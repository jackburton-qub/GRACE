"""
SSR Detection Module — High Performance
----------------------------------------
Detects SSRs (microsatellites) in a genome dictionary:
    genome_dict = {contig_name: sequence}

Performance strategy
--------------------
1. Multiprocessing  — contigs are distributed across all available CPU cores
                      using a process pool. Each worker scans its own batch
                      independently with no shared state.

2. Compiled regex   — one pre-compiled pattern per motif length using a
                      backreference: ([ACGTN]{k})(?:\\1){min_rep-1,}
                      re.finditer() runs in C, far faster than Python loops.

3. Pre-computed     — all possible canonical motifs for each k are computed
   canonical table    upfront (4^k possibilities) and stored in a plain dict.
                      No lru_cache overhead on the hot path.

4. Early skip       — contigs shorter than k * min_rep are skipped immediately.

5. Single uppercase — seq.upper() called once per contig, not per pattern.

6. Chunk-based      — for very large contigs the sequence is scanned in
   scanning           overlapping chunks to keep memory pressure manageable,
                      though for most genomes a full-contig scan is fine.
"""

import re
import os
import multiprocessing as mp
from functools import lru_cache
from itertools import product as iproduct


# ---------------------------------------------------------
# CANONICAL MOTIF
# ---------------------------------------------------------

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def is_compound_motif(motif: str) -> bool:
    """
    Return True if motif is just a shorter unit repeated.
    e.g. CTCT -> True (it's CT*2), ATAT -> True (AT*2),
         CTAG -> False (no shorter period divides it evenly).
    Checks all divisors of len(motif) smaller than len(motif).
    """
    n = len(motif)
    for period in range(1, n):
        if n % period == 0:
            unit = motif[:period]
            if unit * (n // period) == motif:
                return True
    return False


@lru_cache(maxsize=131072)
def canonical_motif(motif: str, level: int = 4) -> str:
    """Compute the canonical form of a motif (cached)."""
    if level == 0:
        return motif
    motif = motif.upper()
    rotations = [motif[i:] + motif[:i] for i in range(len(motif))]
    best = min(rotations)
    if level >= 2:
        rc = reverse_complement(best)
        rc_rots = [rc[i:] + rc[:i] for i in range(len(rc))]
        best = min(best, min(rc_rots))
    return best


def _build_canonical_table(k: int, level: int) -> dict:
    """
    Pre-compute canonical forms for all possible k-mers over ACGT.
    Returns a dict {motif_str: canonical_str}.
    Much faster than calling canonical_motif() on every match at runtime.
    """
    bases = "ACGT"
    table = {}
    for combo in iproduct(bases, repeat=k):
        motif = "".join(combo)
        table[motif] = canonical_motif(motif, level)
    return table


# ---------------------------------------------------------
# PATTERN BUILDER
# ---------------------------------------------------------

def _build_pattern(k: int, min_rep: int) -> re.Pattern:
    """Compile a regex that matches any k-mer repeated >= min_rep times."""
    return re.compile(
        rf"([ACGTN]{{{k}}})(?:\1){{{min_rep - 1},}}",
        re.IGNORECASE,
    )


# ---------------------------------------------------------
# SINGLE-CONTIG WORKER  (must be top-level for pickling)
# ---------------------------------------------------------

def _worker_init():
    """
    Redirect worker stderr to devnull on startup.
    Suppresses BrokenPipeError tracebacks printed by dying workers when the
    pool closes on Windows/Mac — cosmetic noise, not a real error.
    """
    import sys
    import os
    sys.stderr = open(os.devnull, "w")


def _scan_contig(args):
    """
    Scan a single contig for SSRs.
    Args is a tuple for multiprocessing compatibility.
    Patterns and canonical tables are pre-built and passed in to avoid
    rebuilding them for every contig.
    Returns a list of SSR dicts (without ssr_id).
    """
    (
        contig, seq,
        motif_lengths, min_repeats,
        exclude_homopolymers, search_reverse,
        std_level,
        patterns, canonical_tables,   # pre-built, passed in
    ) = args

    seq = seq.upper()
    n   = len(seq)
    results = []

    for k in motif_lengths:
        min_rep = min_repeats.get(k, 3)
        min_len = k * min_rep

        if n < min_len:
            continue

        pattern       = patterns[k]
        canonical_tbl = canonical_tables[k]

        # Forward strand
        for m in pattern.finditer(seq):
            motif_str = m.group(1).upper()

            if exclude_homopolymers and len(set(motif_str)) == 1:
                continue

            if is_compound_motif(motif_str):
                continue

            full_len     = m.end() - m.start()
            repeat_count = full_len // k
            canon        = canonical_tbl.get(motif_str) or canonical_motif(motif_str, std_level)

            results.append({
                "contig":          contig,
                "start":           m.start() + 1,  # 1-based, matches Krait2 convention
                "end":             m.end(),         # 1-based inclusive
                "motif":           motif_str,
                "canonical_motif": canon,
                "repeat_count":    repeat_count,
            })

        # Reverse complement strand (optional)
        if search_reverse:
            rc_seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            for m in pattern.finditer(rc_seq):
                motif_str = m.group(1).upper()

                if exclude_homopolymers and len(set(motif_str)) == 1:
                    continue

                if is_compound_motif(motif_str):
                    continue

                full_len     = m.end() - m.start()
                repeat_count = full_len // k
                canon        = canonical_tbl.get(motif_str) or canonical_motif(motif_str, std_level)

                fwd_start = n - m.end()
                fwd_end   = n - m.start() - 1

                results.append({
                    "contig":          contig,
                    "start":           fwd_start + 1,  # 1-based
                    "end":             fwd_end + 1,     # 1-based inclusive
                    "motif":           motif_str,
                    "canonical_motif": canon,
                    "repeat_count":    repeat_count,
                })

    return results


# ---------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------

def find_ssrs(
    genome_dict: dict,
    motif_lengths: tuple = (2, 3, 4, 5, 6),
    min_repeats: dict = None,
    exclude_homopolymers: bool = True,
    search_reverse: bool = False,
    motif_standardisation_level: int = 4,
    progress_callback=None,
    n_workers: int = None,
    min_contig_len: int = 300,
) -> list:
    """
    Detect SSRs in the genome using multiprocessing and compiled regex.

    Args:
        genome_dict: dict of {contig_name: sequence}
        motif_lengths: tuple of motif sizes to search (e.g. (2, 3, 4))
        min_repeats: dict of {motif_length: min_repeat_count}; defaults to 3 for all
        exclude_homopolymers: skip single-nucleotide repeat motifs (e.g. AAAA)
        search_reverse: also scan reverse complement strand
        motif_standardisation_level: 0-4, controls canonical motif computation
        progress_callback: optional callable(done, total) for UI progress updates
        n_workers: number of worker processes (default: all available CPU cores)
        min_contig_len: skip contigs shorter than this (default 300bp).
                        Contigs this short cannot yield usable SSR markers
                        and skipping them dramatically reduces overhead for
                        fragmented assemblies.

    Returns:
        List of dicts with fields:
            ssr_id, contig, start, end, motif, canonical_motif, repeat_count
    """
    if min_repeats is None:
        min_repeats = {k: 3 for k in motif_lengths}

    if n_workers is None:
        n_workers = max(1, os.cpu_count() or 1)

    contigs       = list(genome_dict.items())
    total_contigs = len(contigs)

    # Filter out contigs that are too short to yield usable SSR markers
    # and build argument tuples for the remaining contigs
    skipped = sum(1 for _, seq in contigs if len(seq) < min_contig_len)
    contigs = [(c, s) for c, s in contigs if len(s) >= min_contig_len]
    total_contigs = len(contigs)

    if skipped:
        import sys
        print(
            f"[ssr_detection] Skipped {skipped:,} contigs shorter than "
            f"{min_contig_len}bp. Scanning {total_contigs:,} contigs.",
            file=sys.stderr,
        )

    if total_contigs == 0:
        return []

    # Pre-build patterns and canonical tables once — not once per contig
    patterns         = {k: _build_pattern(k, min_repeats.get(k, 3)) for k in motif_lengths}
    canonical_tables = {k: _build_canonical_table(k, motif_standardisation_level) for k in motif_lengths}

    task_args = [
        (
            contig, seq,
            motif_lengths, min_repeats,
            exclude_homopolymers, search_reverse,
            motif_standardisation_level,
            patterns, canonical_tables,
        )
        for contig, seq in contigs
    ]

    ssrs = []

    # Use multiprocessing for large genomes, single-process for small ones
    # (process pool overhead isn't worth it for < ~100 contigs)
    # Only fire the progress callback every 1% of contigs to avoid
    # hammering Streamlit with up to 250,000 re-renders on large genomes
    update_interval = max(1, total_contigs // 100)

    if n_workers > 1 and total_contigs >= 10:
        # Tune chunksize: larger chunks reduce IPC overhead for many small contigs
        # Larger chunksize = fewer IPC round trips across process boundary
        # For 250k contigs on 4 cores this gives chunks of ~3000 contigs per batch
        chunksize = max(50, min(5000, total_contigs // (n_workers * 4)))

        # Use spawn context — prevents deadlock on re-run in PyInstaller frozen apps
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_workers, initializer=_worker_init) as pool:
            completed = 0
            for contig_results in pool.imap(
                _scan_contig, task_args, chunksize=chunksize
            ):
                ssrs.extend(contig_results)
                completed += 1
                if progress_callback and (
                    completed % update_interval == 0 or completed == total_contigs
                ):
                    progress_callback(completed, total_contigs)
    else:
        for i, args in enumerate(task_args):
            ssrs.extend(_scan_contig(args))
            if progress_callback and (
                (i + 1) % update_interval == 0 or i + 1 == total_contigs
            ):
                progress_callback(i + 1, total_contigs)

    # Sort by contig then position, then assign sequential IDs
    ssrs.sort(key=lambda x: (x["contig"], x["start"]))
    for i, s in enumerate(ssrs, start=1):
        s["ssr_id"] = i

    return ssrs