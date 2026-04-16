"""
Primer Specificity Module (BLAST-based in-silico PCR)
-----------------------------------------------------
Checks primer pair specificity by BLASTing primers against the genome
and identifying all amplicons that would be produced.

Primer FASTA IDs must use "|" as delimiter:
    >SSR{ssr_id}|L   (left primer)
    >SSR{ssr_id}|R   (right primer)

This is produced automatically by primers_to_blast_fasta() in primer_design.py.
"""

import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import List, Dict, Any, Optional


# ---------------------------------------------------------------------------
# Parameter dataclasses
# ---------------------------------------------------------------------------

@dataclass
class BlastParams:
    task: str = "blastn-short"
    word_size: int = 7
    evalue: float = 0.001
    reward: int = 1
    penalty: int = -3
    gapopen: int = 5
    gapextend: int = 2
    dust: str = "no"
    soft_masking: str = "false"


@dataclass
class SpecificityParams:
    min_product: int = 100
    max_product: int = 300
    max_mismatches: int = 3
    min_align_len: int = 15
    min_identity: float = 80.0
    allow_offtarget_amplicons: bool = False
    max_offtarget_amplicon_size: int = 1000
    pass_mode: str = "amplicons"
    # pass_mode options:
    #   "amplicons" — PASS only if exactly 1 in-range amplicon is formed (default, classic PCR logic)
    #   "hits"      — PASS only if each primer maps to exactly 1 genomic location (stricter)
    max_total_hits: int = 1
    # max_total_hits: when pass_mode="hits", the maximum number of genomic
    # locations each individual primer may match. Default 1 = unique mapping only.


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

PRESET_SPECIFICITY = {
    "Single-locus": SpecificityParams(
        min_product=100,
        max_product=300,
        max_mismatches=0,
        min_align_len=18,
        min_identity=100.0,
        allow_offtarget_amplicons=False,
        max_offtarget_amplicon_size=800,
        pass_mode="hits",
        max_total_hits=1,
    ),
    "Strict": SpecificityParams(
        min_product=100,
        max_product=300,
        max_mismatches=2,
        min_align_len=18,
        min_identity=90.0,
        allow_offtarget_amplicons=False,
        max_offtarget_amplicon_size=800,
    ),
    "Balanced": SpecificityParams(
        min_product=100,
        max_product=300,
        max_mismatches=3,
        min_align_len=16,
        min_identity=85.0,
        allow_offtarget_amplicons=False,
        max_offtarget_amplicon_size=1000,
    ),
    "Lenient": SpecificityParams(
        min_product=80,
        max_product=400,
        max_mismatches=4,
        min_align_len=14,
        min_identity=80.0,
        allow_offtarget_amplicons=True,
        max_offtarget_amplicon_size=1500,
    ),
}

PRESET_BLAST = {
    "Single-locus": BlastParams(
        task="blastn-short",
        word_size=7,
        evalue=0.001,
        reward=1,
        penalty=-3,
        gapopen=5,
        gapextend=2,
        dust="no",
        soft_masking="false",
    ),
    "Strict": BlastParams(
        task="blastn-short",
        word_size=7,
        evalue=0.01,
        reward=1,
        penalty=-3,
        gapopen=5,
        gapextend=2,
        dust="no",
        soft_masking="false",
    ),
    "Balanced": BlastParams(
        task="blastn-short",
        word_size=9,
        evalue=0.1,
        reward=1,
        penalty=-2,
        gapopen=5,
        gapextend=2,
        dust="yes",
        soft_masking="true",
    ),
    "Lenient": BlastParams(
        task="blastn",
        word_size=11,
        evalue=1.0,
        reward=1,
        penalty=-2,
        gapopen=5,
        gapextend=2,
        dust="yes",
        soft_masking="true",
    ),
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _run_cmd(cmd: List[str], cwd: Optional[str] = None) -> str:
    """Run a shell command, raise RuntimeError on failure, return stdout."""
    kwargs = {}
    if sys.platform == "win32":
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW

    result = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        **kwargs,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    return result.stdout


def _write_primers_fasta(primer_results: List[Dict[str, Any]], path: str) -> None:
    """
    Write primers to FASTA using "|" delimiter for safe ID parsing.
    Format: >SSR{ssr_id}|L  and  >SSR{ssr_id}|R
    """
    with open(path, "w") as fh:
        for i, p in enumerate(primer_results):
            ssr_id = p.get("ssr_id", f"pair{i}")
            fh.write(f">SSR{ssr_id}|L\n{p['left_primer']}\n")
            fh.write(f">SSR{ssr_id}|R\n{p['right_primer']}\n")


def _parse_blast_tab(path: str) -> List[Dict[str, Any]]:
    """Parse BLAST tabular output (outfmt 6 + sstrand)."""
    hits = []
    if not os.path.exists(path):
        return hits

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) != 13:
                continue  # skip malformed lines

            (
                qseqid, sseqid, pident, length,
                mismatch, gapopen, qstart, qend,
                sstart, send, evalue, bitscore, sstrand,
            ) = cols

            hits.append({
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": float(pident),
                "length": int(length),
                "mismatch": int(mismatch),
                "gapopen": int(gapopen),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "sstrand": sstrand.strip(),  # strip trailing whitespace/newlines
            })
    return hits


def _pair_hits(
    hits: List[Dict[str, Any]],
    specificity: SpecificityParams,
) -> Dict[str, Dict[str, Any]]:
    """
    Pair left and right primer hits into amplicons and assess specificity.

    Expects query IDs in the format: SSR{ssr_id}|L  and  SSR{ssr_id}|R
    Uses "|" as the delimiter to safely split SSR ID from strand label.

    pass_mode="amplicons" (default):
        PASS = exactly 1 in-range amplicon formed by the primer pair.
        Multiple primers hitting multiple genomic locations is tolerated
        provided only one valid amplicon is produced.

    pass_mode="hits":
        PASS = each individual primer maps to at most max_total_hits
        locations AND exactly 1 in-range amplicon is formed.
        This is stricter: a primer that matches elsewhere in the genome
        (even without forming an amplicon) will cause a FAIL.
    """
    by_query: Dict[str, List[Dict]] = {}
    for h in hits:
        by_query.setdefault(h["qseqid"], []).append(h)

    # Recover unique SSR IDs by splitting on "|"
    ssr_ids = set()
    for q in by_query.keys():
        if "|" in q:
            ssr_ids.add(q.rsplit("|", 1)[0])  # e.g. "SSR42|L" → "SSR42"

    results: Dict[str, Dict[str, Any]] = {}

    def _good(h, sp):
        """Return True if a BLAST hit passes the quality thresholds."""
        return (
            h["length"]   >= sp.min_align_len
            and h["mismatch"] <= sp.max_mismatches
            and h["pident"]   >= sp.min_identity
        )

    def _unique_locations(hit_list):
        """Count distinct (contig, strand, start) locations in a hit list."""
        return {(h["sseqid"], h["sstrand"], h["sstart"]) for h in hit_list}

    for ssr_key in ssr_ids:
        left_hits_all  = by_query.get(f"{ssr_key}|L", [])
        right_hits_all = by_query.get(f"{ssr_key}|R", [])

        left_hits  = [h for h in left_hits_all  if _good(h, specificity)]
        right_hits = [h for h in right_hits_all if _good(h, specificity)]

        left_locations  = _unique_locations(left_hits)
        right_locations = _unique_locations(right_hits)
        n_left_hits     = len(left_locations)
        n_right_hits    = len(right_locations)

        # Find all in-range amplicons.
        # Guard against O(n²) explosion on highly repetitive genomes —
        # if either primer has an unusually large number of hits, cap the
        # search to avoid hanging. Pairs with >200 hits per side are
        # almost certainly non-specific and will FAIL anyway.
        MAX_HITS_PER_SIDE = 200
        left_search  = left_hits[:MAX_HITS_PER_SIDE]
        right_search = right_hits[:MAX_HITS_PER_SIDE]
        amplicons = []
        for lh in left_search:
            for rh in right_search:
                if lh["sseqid"] != rh["sseqid"]:
                    continue
                if lh["sstrand"] == rh["sstrand"]:
                    continue
                coords = [lh["sstart"], lh["send"], rh["sstart"], rh["send"]]
                product = abs(max(coords) - min(coords)) + 1
                if specificity.min_product <= product <= specificity.max_product:
                    amplicons.append({
                        "contig": lh["sseqid"],
                        "product_size": product,
                        "left_hit": lh,
                        "right_hit": rh,
                    })

        n_amplicons = len(amplicons)

        # ── Determine PASS/FAIL ──────────────────────────────
        if specificity.pass_mode == "hits":
            # Strict: each primer must map to <= max_total_hits locations
            # AND exactly one valid amplicon must be formed
            hits_ok = (
                n_left_hits  <= specificity.max_total_hits
                and n_right_hits <= specificity.max_total_hits
            )
            amplicon_ok = n_amplicons == 1
            status = "PASS" if (hits_ok and amplicon_ok) else "FAIL"

        else:
            # Default "amplicons" mode: only the number of formed amplicons matters
            if n_amplicons == 0:
                status = "FAIL"
            elif not specificity.allow_offtarget_amplicons:
                status = "PASS" if n_amplicons == 1 else "FAIL"
            else:
                ok = all(
                    a["product_size"] <= specificity.max_offtarget_amplicon_size
                    for a in amplicons
                )
                status = "PASS" if ok else "FAIL"

        # Strip "SSR" prefix to recover numeric ssr_id for merging
        numeric_id_str = ssr_key[3:] if ssr_key.startswith("SSR") else ssr_key
        try:
            numeric_id = int(numeric_id_str)
        except ValueError:
            numeric_id = numeric_id_str

        results[numeric_id] = {
            "ssr_id":             numeric_id,
            "specificity_status": status,
            "pass_mode":          specificity.pass_mode,
            "n_left_hits":        n_left_hits,
            "n_right_hits":       n_right_hits,
            "n_amplicons":        n_amplicons,
            "amplicons":          amplicons,
            "left_hits":          left_hits,
            "right_hits":         right_hits,
        }

    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def check_specificity_blast(
    genome_fasta: str,
    primer_results: List[Dict[str, Any]],
    blast_params: Optional[BlastParams] = None,
    specificity_params: Optional[SpecificityParams] = None,
    blast_bin_dir: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Run BLAST-based specificity check on primer pairs.

    Args:
        genome_fasta: path to genome FASTA file (used to build BLAST db)
        primer_results: list of dicts with keys: ssr_id, left_primer, right_primer
        blast_params: BlastParams dataclass (defaults to BlastParams())
        specificity_params: SpecificityParams dataclass (defaults to SpecificityParams())
        blast_bin_dir: optional path to BLAST binaries directory

    Returns:
        List of primer result dicts with added keys:
            specificity_status, amplicons, left_hits, right_hits
    """
    if blast_params is None:
        blast_params = BlastParams()
    if specificity_params is None:
        specificity_params = SpecificityParams()

    blastn = "blastn"
    makeblastdb = "makeblastdb"
    if blast_bin_dir:
        blastn = os.path.join(blast_bin_dir, "blastn")
        makeblastdb = os.path.join(blast_bin_dir, "makeblastdb")

    with tempfile.TemporaryDirectory() as tmpdir:
        db_prefix = os.path.join(tmpdir, "genome_db")
        primers_fa = os.path.join(tmpdir, "primers.fasta")
        out_tab = os.path.join(tmpdir, "blast_out.tsv")

        # Build BLAST database
        _run_cmd([
            makeblastdb,
            "-in", genome_fasta,
            "-dbtype", "nucl",
            "-out", db_prefix,
        ])

        # Write primer FASTA
        _write_primers_fasta(primer_results, primers_fa)

        # Run BLAST
        _run_cmd([
            blastn,
            "-task", blast_params.task,
            "-query", primers_fa,
            "-db", db_prefix,
            "-out", out_tab,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand",
            "-word_size", str(blast_params.word_size),
            "-evalue", str(blast_params.evalue),
            "-reward", str(blast_params.reward),
            "-penalty", str(blast_params.penalty),
            "-gapopen", str(blast_params.gapopen),
            "-gapextend", str(blast_params.gapextend),
            "-dust", blast_params.dust,
            "-soft_masking", blast_params.soft_masking,
        ])

        hits_raw = _parse_blast_tab(out_tab)

    hits = hits_raw
    paired = _pair_hits(hits, specificity_params)

    # Merge specificity results back into primer_results by ssr_id
    out: List[Dict[str, Any]] = []
    for p in primer_results:
        ssr_id = p.get("ssr_id")
        spec = paired.get(ssr_id)
        if spec is None:
            merged = dict(p)
            merged["specificity_status"] = "UNKNOWN"
            merged["amplicons"] = []
            merged["left_hits"] = []
            merged["right_hits"] = []
        else:
            merged = dict(p)
            merged.update({
                "specificity_status": spec["specificity_status"],
                "amplicons": spec["amplicons"],
                "left_hits": spec["left_hits"],
                "right_hits": spec["right_hits"],
            })
        out.append(merged)

    return out