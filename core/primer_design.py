"""
Primer Design Module (backend-safe)
-----------------------------------

Specificity is NOT integrated here by design.
Primers can be exported in BLAST-ready FASTA format.

Note on FASTA IDs for BLAST compatibility:
    Uses "|" as delimiter: SSR{id}|L and SSR{id}|R
    This avoids ambiguity when contig names or ssr_ids contain underscores.

Performance note:
    Multiprocessing is intentionally disabled. Spawning child processes from
    inside a PyQt6 QThread on Windows causes deadlocks due to Qt's event loop
    initialisation in the spawned process. primer3 is a C extension and runs
    fast enough in a single QThread for all practical SSR dataset sizes.

GBS mode:
    Genotype-by-sequencing mode uses shorter flanks and smaller product sizes
    to ensure amplicons fit within Illumina read lengths. Primers are tailed
    with standard Illumina adapter sequences for direct use in amplicon
    sequencing workflows.
"""

import primer3
import os
import re


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

STRICT_PRESET = {
    "PRIMER_MIN_SIZE": 20,
    "PRIMER_OPT_SIZE": 22,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_MIN_TM": 59.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 61.0,
    "PRIMER_MIN_GC": 45.0,
    "PRIMER_MAX_GC": 55.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_MAX_SELF_ANY": 4.0,
    "PRIMER_MAX_SELF_END": 1.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 4.0,
    "PRIMER_PAIR_MAX_COMPL_END": 1.0,
    "PRIMER_MAX_HAIRPIN_TH": 12.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 1.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
}

RECOMMENDED_PRESET = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 27,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 30.0,
    "PRIMER_MAX_GC": 70.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_END": 3.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
}

RELAXED_PRESET = {
    "PRIMER_MIN_SIZE": 16,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 30,
    "PRIMER_MIN_TM": 50.0,
    "PRIMER_OPT_TM": 58.0,
    "PRIMER_MAX_TM": 65.0,
    "PRIMER_MIN_GC": 20.0,
    "PRIMER_MAX_GC": 80.0,
    "PRIMER_MAX_POLY_X": 6,
    "PRIMER_MAX_SELF_ANY": 12.0,
    "PRIMER_MAX_SELF_END": 6.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 12.0,
    "PRIMER_PAIR_MAX_COMPL_END": 6.0,
    "PRIMER_MAX_HAIRPIN_TH": 40.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 4.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
}

# GBS preset — optimised for Illumina amplicon sequencing
# Short flanks + small products ensure amplicons fit within 300bp paired-end reads.
# Slightly stricter Tm range for consistent multiplexed PCR.
GBS_PRESET = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_MAX_SELF_ANY": 6.0,
    "PRIMER_MAX_SELF_END": 2.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 6.0,
    "PRIMER_PAIR_MAX_COMPL_END": 2.0,
    "PRIMER_MAX_HAIRPIN_TH": 16.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 1.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
}

PRESETS = {
    "strict":      STRICT_PRESET,
    "recommended": RECOMMENDED_PRESET,
    "relaxed":     RELAXED_PRESET,
    "gbs":         GBS_PRESET,
}

# Default parameters for GBS mode
GBS_DEFAULTS = {
    "flank":            50,     # short flanks to keep product small
    "product_min":      80,
    "product_max":      200,    # must fit in Illumina 300bp paired-end
    "num_pairs":        2,
}

# Standard Illumina adapter tails for M13-tailed primer approach
# Forward primer gets M13F tail, reverse gets M13R tail
# These are added to the 5' end of each primer for adapter ligation
ILLUMINA_TAILS = {
    "forward": "TGTAAAACGACGGCCAGT",   # M13F (-21) universal tail
    "reverse": "CAGGAAACAGCTATGAC",    # M13R universal tail
}

INT_KEYS = {
    "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE", "PRIMER_MAX_POLY_X",
}
FLOAT_KEYS = {
    "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
    "PRIMER_MIN_GC", "PRIMER_MAX_GC",
    "PRIMER_MAX_SELF_ANY", "PRIMER_MAX_SELF_END",
    "PRIMER_PAIR_MAX_COMPL_ANY", "PRIMER_PAIR_MAX_COMPL_END",
    "PRIMER_MAX_HAIRPIN_TH", "PRIMER_PAIR_MAX_DIFF_TM",
    "PRIMER_SALT_MONOVALENT", "PRIMER_DNA_CONC",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalize_types(opts):
    fixed = {}
    for key, value in opts.items():
        if key in INT_KEYS:
            fixed[key] = int(value)
        elif key in FLOAT_KEYS:
            fixed[key] = float(value)
        else:
            fixed[key] = value
    return fixed


def calc_gc(seq):
    seq = seq.upper()
    gc  = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


def calc_tm(seq):
    return round(primer3.calc_tm(seq), 2)


_NN_DG = {
    "AA": -1.0,  "AT": -0.88, "TA": -0.58, "CA": -1.45,
    "GT": -1.44, "CT": -1.28, "GA": -1.30, "CG": -2.17,
    "GC": -2.24, "GG": -1.84, "AC": -1.44, "AG": -1.28,
    "TC": -1.30, "TG": -1.45, "TT": -1.0,  "CC": -1.84,
}


def calc_3prime_dg(seq):
    tetramer = seq[-4:].upper()
    total    = sum(abs(_NN_DG.get(tetramer[i:i+2], -1.0)) for i in range(3))
    return round(total, 2)


# ---------------------------------------------------------------------------
# Low-complexity flank filter
# ---------------------------------------------------------------------------

# Pre-compiled patterns for fast repeat detection
_DI_REPEAT  = re.compile(r'(.{2})\1{3,}', re.IGNORECASE)   # dinucleotide ×4+
_TRI_REPEAT = re.compile(r'(.{3})\1{2,}', re.IGNORECASE)   # trinucleotide ×3+
_MONO_RUN   = re.compile(r'(.)\1{5,}',    re.IGNORECASE)   # mononucleotide ×6+


def _is_low_complexity(seq: str) -> bool:
    """
    Return True if a sequence is dominated by simple repeats and would
    produce unreliable primers. Checks for:
      - Mononucleotide runs of 6+ (e.g. AAAAAA)
      - Dinucleotide repeats of 4+ units (e.g. CTCTCTCT)
      - Trinucleotide repeats of 3+ units (e.g. ATGATGATG)
    """
    if not seq or len(seq) < 8:
        return True
    s = seq.upper()
    if _MONO_RUN.search(s):
        return True
    if _DI_REPEAT.search(s):
        return True
    if _TRI_REPEAT.search(s):
        return True
    return False


def _flank_is_usable(template: str, target_start: int, target_len: int,
                     min_flank: int = 18) -> bool:
    """
    Return True if both flanking sequences are long enough and not
    low-complexity. SSRs that fail this check are skipped before
    calling Primer3, saving significant compute time.
    """
    left_flank  = template[:target_start]
    right_flank = template[target_start + target_len:]

    if len(left_flank) < min_flank or len(right_flank) < min_flank:
        return False
    if _is_low_complexity(left_flank[-30:]):   # check last 30bp of left flank
        return False
    if _is_low_complexity(right_flank[:30]):   # check first 30bp of right flank
        return False
    return True


# ---------------------------------------------------------------------------
# Core design function
# ---------------------------------------------------------------------------

def _design_one(ssr, template, left_bound, target_start, target_len,
                product_size_range, opts, num_pairs,
                gbs_mode=False, add_adapters=False):
    """Design primers for a single SSR. Returns (success_list, failed_list)."""

    primer3_input = {
        "SEQUENCE_ID":               f"{ssr['contig']}|{ssr['ssr_id']}",
        "SEQUENCE_TEMPLATE":         template,
        "SEQUENCE_TARGET":           [target_start, target_len],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_size_range)],
        # Exclude the SSR itself from primer binding — prevents primers
        # landing within the repeat, which causes the low-complexity problem
        "SEQUENCE_EXCLUDED_REGION":  [[target_start, target_len]],
    }

    result    = primer3.design_primers(primer3_input, opts)
    success   = []
    found_any = False

    for rank in range(num_pairs):
        left_key = f"PRIMER_LEFT_{rank}_SEQUENCE"
        if left_key not in result:
            break
        found_any = True
        left_seq  = result[left_key]
        right_seq = result[f"PRIMER_RIGHT_{rank}_SEQUENCE"]

        # Add Illumina adapter tails in GBS mode
        left_seq_with_tail  = left_seq
        right_seq_with_tail = right_seq
        if gbs_mode and add_adapters:
            left_seq_with_tail  = ILLUMINA_TAILS["forward"] + left_seq
            right_seq_with_tail = ILLUMINA_TAILS["reverse"] + right_seq

        rec = {
            "ssr_id":          ssr["ssr_id"],
            "pair_rank":       rank,
            "contig":          ssr["contig"],
            "start":           ssr["start"],
            "end":             ssr["end"],
            "motif":           ssr["motif"],
            "canonical_motif": ssr.get("canonical_motif"),
            "repeat_count":    ssr["repeat_count"],
            "left_primer":     left_seq,
            "right_primer":    right_seq,
            "product_size":    result[f"PRIMER_PAIR_{rank}_PRODUCT_SIZE"],
            "left_gc":         calc_gc(left_seq),
            "right_gc":        calc_gc(right_seq),
            "left_tm":         calc_tm(left_seq),
            "right_tm":        calc_tm(right_seq),
            "left_3end_dg":    calc_3prime_dg(left_seq),
            "right_3end_dg":   calc_3prime_dg(right_seq),
            "gbs_mode":        gbs_mode,
        }

        # Store tailed sequences separately if GBS mode with adapters
        if gbs_mode and add_adapters:
            rec["left_primer_tailed"]  = left_seq_with_tail
            rec["right_primer_tailed"] = right_seq_with_tail

        # Store genomic feature if available
        if "genomic_feature" in ssr:
            rec["genomic_feature"] = ssr["genomic_feature"]

        success.append(rec)

    return (success, [] if found_any else [ssr])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def design_primers_for_all_ssrs(
    genome,
    ssr_list,
    flank=100,
    product_size_range=(100, 300),
    preset="recommended",
    primer_opts=None,
    num_pairs=1,
    progress_callback=None,
    gbs_mode=False,
    add_adapters=False,
):
    """
    Design primers for a list of SSRs using primer3.

    Runs entirely in the calling thread — no multiprocessing. Spawning child
    processes from inside a PyQt6 QThread on Windows causes deadlocks, and
    primer3 (a C extension) is fast enough in a single thread for all
    practical SSR dataset sizes.

    Args:
        genome: dict of {contig_name: sequence}
        ssr_list: list of SSR dicts
        flank: bp of flanking sequence on each side of the SSR
        product_size_range: (min, max) tuple for amplicon size
        preset: one of "strict", "recommended", "relaxed", "gbs"
        primer_opts: dict of primer3 parameter overrides
        num_pairs: number of primer pairs to return per SSR (1–5)
        progress_callback: optional callable(done, total)
        gbs_mode: if True, applies GBS-specific design logic and adds
                  SEQUENCE_EXCLUDED_REGION to keep primers out of the repeat
        add_adapters: if True (GBS mode only), appends Illumina M13 tails
                      to primer sequences for direct use in amplicon sequencing

    Returns:
        dict with keys:
            "success": list of primer result dicts
            "failed":  list of SSR dicts for which no primers were found
            "skipped": list of SSR dicts skipped due to low-complexity flanks
    """
    opts = PRESETS[preset].copy()

    if primer_opts:
        opts.update(primer_opts)

    opts = normalize_types(opts)

    num_pairs = max(1, min(int(num_pairs), 5))
    opts["PRIMER_NUM_RETURN"] = num_pairs

    total     = len(ssr_list)
    success   = []
    failed    = []
    skipped   = []   # low-complexity flanks — not worth BLASTing

    # Minimum primer size determines minimum usable flank
    min_flank = opts.get("PRIMER_MIN_SIZE", 18)

    update_interval = max(1, total // 100)

    for i, ssr in enumerate(ssr_list):
        contig      = ssr["contig"]
        start       = ssr["start"]
        end         = ssr["end"]
        full_seq    = genome[contig]
        left_bound  = max(0, (start - 1) - flank)
        right_bound = min(len(full_seq), end + flank)
        template    = full_seq[left_bound:right_bound]
        target_start = (start - 1) - left_bound
        target_len   = end - (start - 1)

        # Skip SSRs with low-complexity or too-short flanking sequence
        # This prevents Primer3 being called on repeat-dominated flanks
        # which would produce non-specific primers anyway
        if not _flank_is_usable(template, target_start, target_len, min_flank):
            skipped.append(ssr)
            if progress_callback and (
                (i + 1) % update_interval == 0 or i + 1 == total
            ):
                progress_callback(i + 1, total)
            continue

        s, f = _design_one(
            ssr, template, left_bound,
            target_start, target_len,
            product_size_range, opts, num_pairs,
            gbs_mode=gbs_mode,
            add_adapters=add_adapters,
        )
        success.extend(s)
        failed.extend(f)

        if progress_callback and (
            (i + 1) % update_interval == 0 or i + 1 == total
        ):
            progress_callback(i + 1, total)

    return {"success": success, "failed": failed, "skipped": skipped}


def primers_to_blast_fasta(primer_results):
    """
    Convert primer results into a BLAST-ready FASTA string.
    Uses bare primer sequences (no adapter tails) for BLAST.
    """
    lines = []
    for rec in primer_results:
        ssr_id = rec["ssr_id"]
        lines.append(f">SSR{ssr_id}|L")
        lines.append(rec["left_primer"])
        lines.append(f">SSR{ssr_id}|R")
        lines.append(rec["right_primer"])
    return "\n".join(lines) + "\n" if lines else ""


def primers_to_gbs_fasta(primer_results):
    """
    Export GBS primers with Illumina adapter tails if present.
    Falls back to bare sequences if tails were not added.
    For use in ordering primers for amplicon sequencing workflows.
    """
    lines = []
    for rec in primer_results:
        ssr_id = rec["ssr_id"]
        left   = rec.get("left_primer_tailed",  rec["left_primer"])
        right  = rec.get("right_primer_tailed", rec["right_primer"])
        lines.append(f">SSR{ssr_id}|L|tailed")
        lines.append(left)
        lines.append(f">SSR{ssr_id}|R|tailed")
        lines.append(right)
    return "\n".join(lines) + "\n" if lines else ""