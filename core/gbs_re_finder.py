"""
gbs_re_finder.py
----------------
Restriction enzyme site scanning and fragment analysis for GBS-RE.
Optimized with per‑contig fragment caching.
"""

import re
from typing import List, Dict, Any, Tuple, Optional

ENZYMES = {
    "PstI": "CTGCAG",
    "MspI": "CCGG",
    "ApeKI": "GCWGC",
    "MseI": "TTAA",
    "SbfI": "CCTGCAGG",
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
    "TaqI": "TCGA",
    "NlaIII": "CATG",
    "None": "",
}

DEGEN = {
    'W': '[AT]', 'S': '[CG]', 'M': '[AC]', 'K': '[GT]',
    'R': '[AG]', 'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}


def expand_degenerate(seq: str) -> str:
    pattern = ''
    for base in seq.upper():
        pattern += DEGEN.get(base, base)
    return pattern


def compile_enzyme_regex(enzyme_seq: str) -> re.Pattern:
    if not enzyme_seq:
        return None
    return re.compile(expand_degenerate(enzyme_seq), re.IGNORECASE)


def find_cut_sites(sequence: str, enzyme_regex: re.Pattern) -> List[int]:
    if enzyme_regex is None:
        return []
    return [m.start() for m in enzyme_regex.finditer(sequence)]


def build_fragments(contig_seq: str, cut_sites: List[int]) -> List[Tuple[int, int]]:
    fragments = []
    sites = [0] + sorted(cut_sites) + [len(contig_seq)]
    for i in range(len(sites) - 1):
        start = sites[i]
        end = sites[i+1]
        if end - start > 0:
            fragments.append((start, end))
    return fragments


def find_ssrs_in_fragments(
    ssrs: List[Dict[str, Any]],
    genome: Dict[str, str],
    enzyme1: str,
    enzyme2: Optional[str] = None,
    min_fragment_size: int = 100,
    max_fragment_size: int = 400,
    min_distance_from_end: int = 20,
    progress_callback: Optional[callable] = None,
) -> Dict[str, Any]:
    """
    Identify SSRs that lie within restriction fragments of desired size.
    Optimized: fragments are pre‑computed once per contig.
    """
    regex1 = compile_enzyme_regex(ENZYMES.get(enzyme1, ""))
    regex2 = compile_enzyme_regex(ENZYMES.get(enzyme2, "")) if enzyme2 and enzyme2 != "None" else None

    contig_fragments = {}
    contig_total_frags = 0
    contig_passing_frags = 0
    passing_fragment_sizes = []

    ssrs_by_contig: Dict[str, List[Dict]] = {}
    for ssr in ssrs:
        contig = ssr.get("contig")
        if contig and contig in genome:
            ssrs_by_contig.setdefault(contig, []).append(ssr)

    for contig, seq in genome.items():
        if contig not in ssrs_by_contig:
            continue

        cuts = find_cut_sites(seq, regex1)
        if regex2:
            cuts.extend(find_cut_sites(seq, regex2))
        cuts = sorted(set(cuts))

        all_frags = build_fragments(seq, cuts)
        contig_total_frags += len(all_frags)

        passing_frags = []
        for start, end in all_frags:
            frag_len = end - start
            if min_fragment_size <= frag_len <= max_fragment_size:
                passing_frags.append((start, end))
                passing_fragment_sizes.append(frag_len)

        contig_passing_frags += len(passing_frags)
        contig_fragments[contig] = passing_frags

    qualified_ssrs = []
    total_ssrs = len(ssrs)

    for idx, ssr in enumerate(ssrs):
        if progress_callback and idx % 100 == 0:
            progress_callback(idx + 1, total_ssrs)

        contig = ssr.get("contig")
        ssr_start = ssr.get("start", 0) - 1
        ssr_end = ssr.get("end", 0)

        fragments = contig_fragments.get(contig, [])
        for frag_start, frag_end in fragments:
            if (ssr_start >= frag_start + min_distance_from_end and
                ssr_end <= frag_end - min_distance_from_end):
                qualified_ssrs.append({
                    **ssr,
                    "fragment_start": frag_start + 1,
                    "fragment_end": frag_end,
                    "fragment_size": frag_end - frag_start,
                })
                break

    if progress_callback:
        progress_callback(total_ssrs, total_ssrs)

    return {
        "qualified_ssrs": qualified_ssrs,
        "total_fragments": contig_total_frags,
        "passing_fragments": contig_passing_frags,
        "passing_fragment_sizes": passing_fragment_sizes,
        "n_qualified": len(qualified_ssrs),
    }