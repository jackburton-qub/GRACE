"""
gbs_enzyme_check.py — Deep diagnostic version
"""
import re
import os
from typing import List, Dict, Any, Optional

LOG_FILE = os.path.join(os.getcwd(), "enzyme_deep_debug.log")

def log_debug(msg: str):
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

# Clear log
with open(LOG_FILE, "w", encoding="utf-8") as f:
    f.write("=== ENZYME DEEP DEBUG LOG ===\n")

DEGEN = {
    'W': '[AT]', 'S': '[CG]', 'M': '[AC]', 'K': '[GT]',
    'R': '[AG]', 'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}

ENZYMES = {
    "PstI":     {"site": "CTGCAG"},
    "MspI":     {"site": "CCGG"},
    "ApeKI":    {"site": "GCWGC"},
    "MseI":     {"site": "TTAA"},
    "SbfI":     {"site": "CCTGCAGG"},
    "EcoRI":    {"site": "GAATTC"},
    "HindIII":  {"site": "AAGCTT"},
    "TaqI":     {"site": "TCGA"},
    "NlaIII":   {"site": "CATG"},
    "PstI-MspI": {"site": "CTGCAG|CCGG"},
}

_COMP = str.maketrans("ACGT", "TGCA")

def _rc(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]

def expand_degenerate(seq: str) -> str:
    pattern = ''
    for base in seq.upper():
        pattern += DEGEN.get(base, base)
    return pattern

def compile_enzyme_regex(enzyme_name: str) -> re.Pattern:
    if enzyme_name not in ENZYMES:
        raise ValueError(f"Unknown enzyme: {enzyme_name}")
    site_str = ENZYMES[enzyme_name]["site"]
    if '|' in site_str:
        patterns = [expand_degenerate(p) for p in site_str.split('|')]
        return re.compile('|'.join(patterns), re.IGNORECASE)
    else:
        return re.compile(expand_degenerate(site_str), re.IGNORECASE)

def extract_amplicon_robust(
    genome: Dict[str, str],
    primer: Dict[str, Any],
    search_flank: int = 500,
) -> Optional[str]:
    contig = primer.get("contig", "")
    ssr_start = primer.get("start", 0) - 1
    ssr_end = primer.get("end", 0)
    product_size = primer.get("product_size", 200)
    ssr_id = primer.get("ssr_id", "?")

    if contig not in genome:
        log_debug(f"FAIL: Contig '{contig}' not in genome for SSR {ssr_id}")
        return None

    seq = genome[contig].upper()
    margin = max(search_flank, product_size + 2000)
    s_start = max(0, ssr_start - margin)
    s_end = min(len(seq), ssr_end + margin)
    region = seq[s_start:s_end]

    fwd = primer.get("left_primer", "").upper()
    rev = primer.get("right_primer", "").upper()
    rev_rc = _rc(rev)

    fwd_pos = region.find(fwd)
    rev_pos = region.rfind(rev_rc)

    if fwd_pos == -1 or rev_pos == -1 or rev_pos <= fwd_pos:
        # Try alternative orientation
        fwd_rc = _rc(fwd)
        fwd_pos2 = region.find(fwd_rc)
        rev_pos2 = region.rfind(rev)
        if fwd_pos2 != -1 and rev_pos2 != -1 and rev_pos2 > fwd_pos2:
            fwd_pos, rev_pos = fwd_pos2, rev_pos2
            rev_rc = rev
        else:
            return None

    amplicon = region[fwd_pos : rev_pos + len(rev_rc)]
    return amplicon

def check_enzyme_sites(primers, genome, enzyme, search_flank=500):
    log_debug(f"=== Enzyme check: {len(primers)} primers, enzyme={enzyme}, flank={search_flank} ===")
    if enzyme not in ENZYMES:
        raise ValueError(f"Unknown enzyme: {enzyme}")
    
    site_regex = compile_enzyme_regex(enzyme)
    log_debug(f"Regex pattern: {site_regex.pattern}")
    
    passed, failed, skipped = [], [], []
    fail_count = 0

    for i, p in enumerate(primers):
        ssr_id = p.get("ssr_id", "?")
        amplicon = extract_amplicon_robust(genome, p, search_flank)
        
        if amplicon is None:
            skipped.append({"ssr_id": ssr_id, "reason": "Extraction failed"})
            continue
        
        if site_regex.search(amplicon):
            passed.append({"ssr_id": ssr_id, "amplicon_length": len(amplicon)})
        else:
            failed.append({"ssr_id": ssr_id, "amplicon_length": len(amplicon), "reason": f"No {enzyme} site"})
            fail_count += 1
            if fail_count <= 10:
                log_debug(f"\n--- FAIL #{fail_count}: SSR {ssr_id} ---")
                log_debug(f"Amplicon length: {len(amplicon)} bp")
                log_debug(f"Amplicon sequence (first 200 bp):\n{amplicon[:200]}")
                # Also check for simple substring matches
                for site in ENZYMES[enzyme]["site"].split('|'):
                    if site.upper() in amplicon.upper():
                        log_debug(f"  WARNING: Simple string match for '{site}' FOUND, but regex failed!")
                    else:
                        log_debug(f"  Simple string match for '{site}' NOT found.")
    
    log_debug(f"=== Complete: {len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped ===\n")
    return {
        "enzyme": enzyme, "n_total": len(primers), "n_passed": len(passed),
        "n_failed": len(failed), "n_skipped": len(skipped),
        "passed_loci": passed, "failed_loci": failed, "skipped_loci": skipped,
        "recommendation": ""
    }