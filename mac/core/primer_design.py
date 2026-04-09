"""
Primer Design Module (backend-safe)
-----------------------------------
Provides strict, recommended, and relaxed primer3 presets.
Allows full override of any primer3 parameter.

Specificity is NOT integrated here by design.
Primers can be exported in BLAST-ready FASTA format.

Note on FASTA IDs for BLAST compatibility:
    Uses "|" as delimiter: SSR{id}|L and SSR{id}|R
    This avoids ambiguity when contig names or ssr_ids contain underscores.
"""

import primer3
import os
import multiprocessing as mp


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


PRESETS = {
    "strict": STRICT_PRESET,
    "recommended": RECOMMENDED_PRESET,
    "relaxed": RELAXED_PRESET,
}


INT_KEYS = {
    "PRIMER_MIN_SIZE",
    "PRIMER_OPT_SIZE",
    "PRIMER_MAX_SIZE",
    "PRIMER_MAX_POLY_X",
}

FLOAT_KEYS = {
    "PRIMER_MIN_TM",
    "PRIMER_OPT_TM",
    "PRIMER_MAX_TM",
    "PRIMER_MIN_GC",
    "PRIMER_MAX_GC",
    "PRIMER_MAX_SELF_ANY",
    "PRIMER_MAX_SELF_END",
    "PRIMER_PAIR_MAX_COMPL_ANY",
    "PRIMER_PAIR_MAX_COMPL_END",
    "PRIMER_MAX_HAIRPIN_TH",
    "PRIMER_PAIR_MAX_DIFF_TM",
    "PRIMER_SALT_MONOVALENT",
    "PRIMER_DNA_CONC",
}


def normalize_types(opts):
    """Cast primer3 option values to correct Python types."""
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
    """Return GC content as a percentage (0–100), rounded to 2 dp."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


def calc_tm(seq):
    """Return Tm in degrees Celsius using primer3, rounded to 2 dp."""
    return round(primer3.calc_tm(seq), 2)


# Nearest-neighbour ΔG for the terminal dinucleotide (SantaLucia 1998,
# 50 mM NaCl, 37°C), kcal/mol. Fast dict lookup replaces calling
# primer3.calcHeterodimer() for every primer (32,000+ calls on large datasets).
_NN_DG = {
    "AA": -1.0,  "AT": -0.88, "TA": -0.58, "CA": -1.45,
    "GT": -1.44, "CT": -1.28, "GA": -1.30, "CG": -2.17,
    "GC": -2.24, "GG": -1.84, "AC": -1.44, "AG": -1.28,
    "TC": -1.30, "TG": -1.45, "TT": -1.0,  "CC": -1.84,
}


def calc_3prime_dg(seq):
    """
    Return the 3' end stability as a positive kcal/mol value by summing the
    nearest-neighbour ΔG for the three dinucleotide steps in the last 4 bases
    (SantaLucia 1998). Higher = more stable 3' end (matches Krait2 convention).
    Fast dict lookup, no primer3 call.
    """
    tetramer = seq[-4:].upper()
    total = sum(
        abs(_NN_DG.get(tetramer[i:i+2], -1.0))
        for i in range(3)
    )
    return round(total, 2)


def _design_one(args):
    """
    Worker function for multiprocessing — designs primers for a single SSR.
    Must be top-level for pickling compatibility.
    """
    (
        ssr, template, left_bound,
        target_start, target_len,
        product_size_range, opts, num_pairs,
    ) = args

    primer3_input = {
        "SEQUENCE_ID":           f"{ssr['contig']}|{ssr['ssr_id']}",
        "SEQUENCE_TEMPLATE":     template,
        "SEQUENCE_TARGET":       [target_start, target_len],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_size_range)],
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
        success.append({
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
        })

    return (success, [] if found_any else [ssr])


def design_primers_for_all_ssrs(
    genome,
    ssr_list,
    flank=100,
    product_size_range=(100, 300),
    preset="recommended",
    primer_opts=None,
    num_pairs=1,
    progress_callback=None,
):
    """
    Design primers for a list of SSRs using primer3.

    Args:
        genome: dict of {contig_name: sequence}
        ssr_list: list of SSR dicts (must include contig, start, end, ssr_id, motif, repeat_count)
        flank: bp of flanking sequence on each side of the SSR
        product_size_range: (min, max) tuple for amplicon size
        preset: one of "strict", "recommended", "relaxed"
        primer_opts: dict of primer3 parameter overrides (applied on top of preset)
        num_pairs: number of primer pairs to return per SSR (1–5); pair 0 is top-ranked
        progress_callback: optional callable(done, total) for UI progress updates

    Returns:
        dict with keys:
            "success": list of primer result dicts (includes pair_rank column)
            "failed": list of SSR dicts for which no primers were found
    """
    opts = PRESETS[preset].copy()

    if primer_opts:
        opts.update(primer_opts)

    opts = normalize_types(opts)

    # Clamp num_pairs and tell Primer3 how many pairs to return
    num_pairs = max(1, min(int(num_pairs), 5))
    opts["PRIMER_NUM_RETURN"] = num_pairs

    total = len(ssr_list)
    n_workers = max(1, os.cpu_count() or 1)
    update_interval = max(1, total // 100)

    # Build task args — extract template sequences here so workers are stateless
    task_args = []
    for ssr in ssr_list:
        contig   = ssr["contig"]
        start    = ssr["start"]
        end      = ssr["end"]
        full_seq = genome[contig]
        # Convert back to 0-based for slicing (coordinates stored as 1-based)
        left_bound  = max(0, (start - 1) - flank)
        right_bound = min(len(full_seq), end + flank)
        template    = full_seq[left_bound:right_bound]
        target_start = (start - 1) - left_bound  # 0-based offset into template
        target_len   = end - (start - 1)
        task_args.append((
            ssr, template, left_bound,
            target_start, target_len,
            product_size_range, opts, num_pairs,
        ))

    success = []
    failed  = []

    def _process_results(results_iter):
        completed = 0
        for result_batch in results_iter:
            s, f = result_batch
            success.extend(s)
            failed.extend(f)
            completed += 1
            if progress_callback and (
                completed % update_interval == 0 or completed == total
            ):
                progress_callback(completed, total)

    if n_workers > 1 and total >= 50:
        chunksize = max(10, min(500, total // (n_workers * 4)))
        # Use spawn context — prevents deadlock on re-run in PyInstaller frozen apps
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_workers) as pool:
            _process_results(pool.imap(_design_one, task_args, chunksize=chunksize))
    else:
        _process_results(_design_one(a) for a in task_args)

    return {"success": success, "failed": failed}


def primers_to_blast_fasta(primer_results):
    """
    Convert primer results into a BLAST-ready FASTA string.

    Uses "|" as delimiter in sequence IDs to avoid conflicts with
    underscores in contig names or SSR IDs.

    Format: >SSR{ssr_id}|L  and  >SSR{ssr_id}|R
    """
    lines = []
    for rec in primer_results:
        ssr_id = rec["ssr_id"]
        left_id = f"SSR{ssr_id}|L"
        right_id = f"SSR{ssr_id}|R"

        lines.append(f">{left_id}")
        lines.append(rec["left_primer"])
        lines.append(f">{right_id}")
        lines.append(rec["right_primer"])

    return "\n".join(lines) + "\n" if lines else ""