"""
amplicon_reference.py
---------------------
Generates expected amplicon sequences and performs GC analysis.
"""

from typing import List, Dict, Any, Optional

_COMP = str.maketrans("ACGT", "TGCA")

def _rc(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]

def extract_amplicon(
    genome: dict,
    primer: Dict[str, Any],
    search_flank: int = 150,
) -> Optional[str]:
    """
    Extract the expected amplicon sequence for a primer pair.
    """
    contig = primer.get("contig", "")
    ssr_start = primer.get("start", 0) - 1  # 0-based
    ssr_end = primer.get("end", 0)

    if contig not in genome:
        return None

    seq = genome[contig].upper()
    s_start = max(0, ssr_start - search_flank)
    s_end = min(len(seq), ssr_end + search_flank)
    region = seq[s_start:s_end]

    fwd = primer.get("left_primer", "").upper()
    rev = primer.get("right_primer", "").upper()
    rev_rc = _rc(rev)

    fwd_pos = region.find(fwd)
    rev_pos = region.rfind(rev_rc)

    if fwd_pos == -1 or rev_pos == -1 or rev_pos <= fwd_pos:
        return None

    return region[fwd_pos : rev_pos + len(rev_rc)]

# Alias for the UI call
extract_amplicon_robust = extract_amplicon

def analyse_amplicon_gc(amplicon_seq: str) -> dict:
    """Return GC percentage and a risk flag."""
    if not amplicon_seq:
        return {"gc_percent": 0.0, "flag": "Unknown"}
    gc = (amplicon_seq.count('G') + amplicon_seq.count('C')) / len(amplicon_seq) * 100
    if gc < 30:
        flag = "Low GC (risk of poor amplification)"
    elif gc > 70:
        flag = "High GC (risk of poor amplification)"
    elif gc < 35 or gc > 65:
        flag = "Marginal GC"
    else:
        flag = "Good"
    return {"gc_percent": round(gc, 1), "flag": flag}

def amplicons_to_fasta(
    pass_primers: List[Dict[str, Any]],
    genome: dict,
    search_flank: int = 150,
):
    """Export expected amplicon sequences as FASTA."""
    lines = []
    n_found = 0
    n_miss = 0

    for p in pass_primers:
        amplicon = extract_amplicon(genome, p, search_flank)
        if not amplicon:
            n_miss += 1
            continue

        ssr_id = p.get("ssr_id", "?")
        contig = p.get("contig", "")
        motif = p.get("canonical_motif", p.get("motif", ""))
        repeats = p.get("repeat_count", "?")
        prod = p.get("product_size", len(amplicon))

        header = f">SSR{ssr_id} contig={contig} motif={motif} repeats={repeats} product={prod}bp"
        lines.append(header)
        lines.append(amplicon)
        n_found += 1

    return "\n".join(lines) + "\n" if lines else "", n_found, n_miss