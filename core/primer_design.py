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

amplicon mode:
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

# RECOMMENDED_PRESET - Used for Capillary Electrophoresis mode (default)
# Standard parameters for SSR genotyping with fragment analysis
RECOMMENDED_PRESET = {
    # Core primer constraints - MATCHED TO KRAIT2 ACTUAL SETTINGS
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,      # Krait2 screenshot shows 20, not 22 as stated in paper
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_MIN_TM": 52.0,
    "PRIMER_OPT_TM": 58.0,
    "PRIMER_MAX_TM": 60.0,
    "PRIMER_MIN_GC": 30.0,
    "PRIMER_OPT_GC": 40.0,      # Krait2 uses 40%
    "PRIMER_MAX_GC": 60.0,
    
    # Secondary structure
    "PRIMER_MAX_POLY_X": 5,
    "PRIMER_GC_CLAMP": 2,
    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_END": 3.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    
    # Pair constraints
    "PRIMER_PAIR_MAX_DIFF_TM": 5.0,
    
    # Thermodynamic parameters - OLD DEFAULTS to match Krait2
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 0.0,
    "PRIMER_DNTP_CONC": 0.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_SALT_CORRECTIONS": 0,
    "PRIMER_TM_FORMULA": 0,
}

# amplicon_PRESET - Used for Amplicon Sequencing mode
# Optimized for Illumina amplicon sequencing
# Stricter Tm range and shorter products to fit 300bp paired-end reads
amplicon_PRESET = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_GC_CLAMP": 2,
    "PRIMER_MAX_SELF_ANY": 6.0,
    "PRIMER_MAX_SELF_END": 2.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 6.0,
    "PRIMER_PAIR_MAX_COMPL_END": 2.0,
    "PRIMER_MAX_HAIRPIN_TH": 16.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 1.0,
    # Thermodynamic parameters
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 0.0,
    "PRIMER_DNTP_CONC": 0.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_SALT_CORRECTIONS": 0,
    "PRIMER_TM_FORMULA": 0,
}

PRESETS = {
    "recommended": RECOMMENDED_PRESET,
    "amplicon":    amplicon_PRESET,
}

INT_KEYS = {
    "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE", "PRIMER_MAX_POLY_X",
    "PRIMER_GC_CLAMP", "PRIMER_SALT_CORRECTIONS", "PRIMER_TM_FORMULA",
}
FLOAT_KEYS = {
    "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
    "PRIMER_MIN_GC", "PRIMER_MAX_GC",
    "PRIMER_MAX_SELF_ANY", "PRIMER_MAX_SELF_END",
    "PRIMER_PAIR_MAX_COMPL_ANY", "PRIMER_PAIR_MAX_COMPL_END",
    "PRIMER_MAX_HAIRPIN_TH", "PRIMER_PAIR_MAX_DIFF_TM",
    "PRIMER_SALT_MONOVALENT", "PRIMER_SALT_DIVALENT", "PRIMER_DNTP_CONC",
    "PRIMER_DNA_CONC",
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


# Nearest-neighbor thermodynamic parameters (kcal/mol)
# These are NEGATIVE values - more negative = stronger binding
_NN_DG = {
    "AA": -1.0,  "AT": -0.88, "TA": -0.58, "CA": -1.45,
    "GT": -1.44, "CT": -1.28, "GA": -1.30, "CG": -2.17,
    "GC": -2.24, "GG": -1.84, "AC": -1.44, "AG": -1.28,
    "TC": -1.30, "TG": -1.45, "TT": -1.0,  "CC": -1.84,
}


def calc_3prime_dg(seq):
    """
    Calculate 3' end stability using nearest-neighbor thermodynamics.
    
    Evaluates the last 5 nucleotides (4 dinucleotide pairs) of the primer
    and sums their free energy (ΔG) values.
    
    Returns:
        ΔG in kcal/mol (NEGATIVE values)
        More negative = more stable binding
        Less negative = weaker binding
    
    Optimal range: -9 to -6 kcal/mol
        - Too stable (< -9 kcal/mol): Risk of non-specific priming
        - Too weak (> -6 kcal/mol): Poor primer extension by polymerase
    
    Note: Uses simplified nearest-neighbor parameters. For more accurate
    calculations accounting for salt concentration, use primer3.calc_end_stability()
    """
    if len(seq) < 5:
        return 0.0
    
    pentamer = seq[-5:].upper()
    # Sum the ΔG values (keep them NEGATIVE - do NOT use abs())
    total = sum(_NN_DG.get(pentamer[i:i+2], -1.0) for i in range(4))
    return round(total, 2)


def calc_3prime_dg_accurate(seq, mv_conc=50.0, dv_conc=0.0):
    """
    Calculate 3' end stability using Primer3's thermodynamic engine.
    
    This is more accurate than calc_3prime_dg() as it accounts for:
    - Salt concentration effects
    - Temperature-dependent parameters
    - Terminal base penalties
    
    Args:
        seq: Primer sequence
        mv_conc: Monovalent cation concentration in mM (default 50.0)
        dv_conc: Divalent cation concentration in mM (default 0.0)
    
    Returns:
        ΔG in kcal/mol (NEGATIVE values)
    """
    # primer3.calc_end_stability returns a ThermoResult object
    # Extract the dg (delta G) attribute
    result = primer3.calc_end_stability(seq, mv_conc, dv_conc)
    return round(result.dg, 2)


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
                amplicon_mode=False, add_adapters=False):
    """Design primers for a single SSR. Returns (success_list, failed_list)."""

    primer3_input = {
        "SEQUENCE_ID":               f"{ssr['contig']}|{ssr['ssr_id']}",
        "SEQUENCE_TEMPLATE":         template,
        "SEQUENCE_TARGET":           [target_start, target_len],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_size_range)],
        # Keep stricter masking: exclude the SSR itself from primer binding
        # This prevents primers landing within the repeat 
        "SEQUENCE_EXCLUDED_REGION":  [[target_start, target_len]],
    }

    result    = primer3.design_primers(primer3_input, opts)
    success   = []
    found_any = False

    # Get salt concentrations for accurate 3' stability calculation
    mv_conc = opts.get("PRIMER_SALT_MONOVALENT", 50.0)
    dv_conc = opts.get("PRIMER_SALT_DIVALENT", 0.0)

    for rank in range(num_pairs):
        left_key = f"PRIMER_LEFT_{rank}_SEQUENCE"
        if left_key not in result:
            break
        found_any = True
        left_seq  = result[left_key]
        right_seq = result[f"PRIMER_RIGHT_{rank}_SEQUENCE"]

        # Add Illumina adapter tails in amplicon mode
        left_seq_with_tail  = left_seq
        right_seq_with_tail = right_seq
        if amplicon_mode and add_adapters:
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
            "left_tm":         round(result[f"PRIMER_LEFT_{rank}_TM"], 2),
            "right_tm":        round(result[f"PRIMER_RIGHT_{rank}_TM"], 2),
            # Use standard thermodynamic ΔG calculation (negative values, -6 to -9 kcal/mol)
            "left_3end_dg":    calc_3prime_dg(left_seq),
            "right_3end_dg":   calc_3prime_dg(right_seq),
            "amplicon_mode":   amplicon_mode,
        }

        # Store tailed sequences separately if amplicon mode with adapters
        if amplicon_mode and add_adapters:
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
    product_size_range=(100, 350),
    preset="recommended",
    primer_opts=None,
    num_pairs=1,
    progress_callback=None,
    amplicon_mode=False,
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
        preset: one of "strict", "recommended", "relaxed", "amplicon"
        primer_opts: dict of primer3 parameter overrides
        num_pairs: number of primer pairs to return per SSR (1–5)
        progress_callback: optional callable(done, total)
        amplicon_mode: if True, applies amplicon-specific design logic and adds
                  SEQUENCE_EXCLUDED_REGION to keep primers out of the repeat
        add_adapters: if True (amplicon mode only), appends Illumina M13 tails
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
            amplicon_mode=amplicon_mode,
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


def primers_to_amplicon_fasta(primer_results):
    """
    Export amplicon primers with Illumina adapter tails if present.
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


def calculate_primer_quality_score(primer_rec, mode="balanced"):
    """
    Calculate a quality score for a primer pair based on multiple criteria.
    Higher score = better quality.
    
    Args:
        primer_rec: Dictionary containing primer pair data
        mode: Scoring strategy - "balanced", "amplicon", or "capillary"
    
    Returns:
        Dictionary with score, rank category, and score breakdown
    
    Scoring criteria:
        - 3' end stability (optimal range: -8.5 to -6.5 kcal/mol)
        - Tm matching between primers (closer = better)
        - GC content balance (45-55% optimal)
        - Product size (smaller often better for amplicon sequencing)
        - Tm in optimal range (57-60°C)
    """
    score = 0
    breakdown = {}
    
    # Extract values
    left_dg = primer_rec.get("left_3end_dg", 0)
    right_dg = primer_rec.get("right_3end_dg", 0)
    left_tm = primer_rec.get("left_tm", 0)
    right_tm = primer_rec.get("right_tm", 0)
    left_gc = primer_rec.get("left_gc", 0)
    right_gc = primer_rec.get("right_gc", 0)
    product_size = primer_rec.get("product_size", 0)
    
    # 1. 3' End Stability (up to 20 points)
    # Optimal range: -8.5 to -6.5 kcal/mol
    left_dg_score = 0
    right_dg_score = 0
    
    if -8.5 <= left_dg <= -6.5:
        left_dg_score = 10
    elif -9.0 <= left_dg <= -6.0:
        left_dg_score = 7
    elif -10.0 <= left_dg <= -5.5:
        left_dg_score = 4
    else:
        left_dg_score = 0
    
    if -8.5 <= right_dg <= -6.5:
        right_dg_score = 10
    elif -9.0 <= right_dg <= -6.0:
        right_dg_score = 7
    elif -10.0 <= right_dg <= -5.5:
        right_dg_score = 4
    else:
        right_dg_score = 0
    
    score += left_dg_score + right_dg_score
    breakdown["3end_stability"] = left_dg_score + right_dg_score
    
    # 2. Tm Matching (up to 15 points)
    # Primers should have similar Tm
    tm_diff = abs(left_tm - right_tm)
    if tm_diff < 1.0:
        tm_match_score = 15
    elif tm_diff < 2.0:
        tm_match_score = 10
    elif tm_diff < 3.0:
        tm_match_score = 5
    else:
        tm_match_score = 0
    
    score += tm_match_score
    breakdown["tm_matching"] = tm_match_score
    
    # 3. Tm in Optimal Range (up to 15 points)
    # Optimal: 57-60°C
    left_tm_score = 0
    right_tm_score = 0
    
    if 57 <= left_tm <= 60:
        left_tm_score = 7.5
    elif 55 <= left_tm <= 62:
        left_tm_score = 5
    elif 52 <= left_tm <= 65:
        left_tm_score = 2
    
    if 57 <= right_tm <= 60:
        right_tm_score = 7.5
    elif 55 <= right_tm <= 62:
        right_tm_score = 5
    elif 52 <= right_tm <= 65:
        right_tm_score = 2
    
    score += left_tm_score + right_tm_score
    breakdown["tm_range"] = left_tm_score + right_tm_score
    
    # 4. GC Content Balance (up to 20 points)
    # Optimal: 45-55%
    left_gc_score = 0
    right_gc_score = 0
    
    if 45 <= left_gc <= 55:
        left_gc_score = 10
    elif 40 <= left_gc <= 60:
        left_gc_score = 7
    elif 35 <= left_gc <= 65:
        left_gc_score = 4
    
    if 45 <= right_gc <= 55:
        right_gc_score = 10
    elif 40 <= right_gc <= 60:
        right_gc_score = 7
    elif 35 <= right_gc <= 65:
        right_gc_score = 4
    
    score += left_gc_score + right_gc_score
    breakdown["gc_content"] = left_gc_score + right_gc_score
    
    # 5. Product Size (up to 15 points)
    # Mode-specific scoring
    if mode == "amplicon":
        # Amplicon sequencing: smaller is better (80-150bp ideal)
        if 80 <= product_size <= 150:
            size_score = 15
        elif 150 < product_size <= 200:
            size_score = 10
        elif 200 < product_size <= 250:
            size_score = 5
        else:
            size_score = 0
    elif mode == "capillary":
        # Capillary: medium sizes often better (150-300bp)
        if 150 <= product_size <= 300:
            size_score = 15
        elif 100 <= product_size < 150 or 300 < product_size <= 350:
            size_score = 10
        else:
            size_score = 5
    else:  # balanced
        # Balanced: moderate preference for smaller
        if 100 <= product_size <= 200:
            size_score = 15
        elif 200 < product_size <= 300:
            size_score = 10
        else:
            size_score = 5
    
    score += size_score
    breakdown["product_size"] = size_score
    
    # 6. Pair Rank Bonus (up to 15 points)
    # Primer3 ranks pairs by quality - rank 0 is best
    pair_rank = primer_rec.get("pair_rank", 0)
    if pair_rank == 0:
        rank_score = 15
    elif pair_rank == 1:
        rank_score = 10
    elif pair_rank == 2:
        rank_score = 5
    else:
        rank_score = 0
    
    score += rank_score
    breakdown["primer3_rank"] = rank_score
    
    # Total possible: 100 points
    # Determine quality category
    if score >= 80:
        category = "Excellent"
    elif score >= 65:
        category = "Good"
    elif score >= 50:
        category = "Acceptable"
    else:
        category = "Poor"
    
    return {
        "quality_score": round(score, 1),
        "quality_category": category,
        "score_breakdown": breakdown,
    }


def add_quality_scores_to_primers(primer_results, mode="balanced"):
    """
    Add quality scores to a list of primer results.
    
    Args:
        primer_results: List of primer result dictionaries
        mode: Scoring mode - "balanced", "amplicon", or "capillary"
    
    Returns:
        List of primer results with added quality_score and quality_category fields
    """
    scored_primers = []
    for primer in primer_results:
        primer_copy = primer.copy()
        score_data = calculate_primer_quality_score(primer, mode=mode)
        primer_copy["quality_score"] = score_data["quality_score"]
        primer_copy["quality_category"] = score_data["quality_category"]
        # Optionally store breakdown for detailed analysis
        primer_copy["score_breakdown"] = score_data["score_breakdown"]
        scored_primers.append(primer_copy)
    
    return scored_primers