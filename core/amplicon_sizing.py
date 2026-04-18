"""
gbs_allele_sizes.py
-------------------
Calculates expected allele size ranges and filters based on platform limits.
"""

from typing import List, Dict, Any

PLATFORMS = [
    {
        "name": "Illumina MiSeq 150bp SE",
        "read_length": 150,
        "min_amplicon": 70,
        "max_amplicon": 130,
        "paired": False,
    },
    {
        "name": "Illumina MiSeq 300bp PE",
        "read_length": 300,
        "min_amplicon": 80,
        "max_amplicon": 450,
        "paired": True,
    },
    {
        "name": "Illumina MiSeq v3 600bp",
        "read_length": 300,
        "min_amplicon": 100,
        "max_amplicon": 550,
        "paired": True,
    },
    {
        "name": "Oxford Nanopore (MinION)",
        "read_length": 10000,
        "min_amplicon": 100,
        "max_amplicon": 2000,
        "paired": False,
    },
    {
        "name": "PacBio SMRT",
        "read_length": 15000,
        "min_amplicon": 150,
        "max_amplicon": 3000,
        "paired": False,
    },
]


def calculate_allele_sizes(
    primer: Dict[str, Any],
    min_repeats: int = 3,
    max_repeats: int = 60,
) -> Dict[str, Any]:
    motif_len = len(primer.get("motif", "AT"))
    ref_repeats = primer.get("repeat_count", 10)
    product_size = primer.get("product_size", 150)

    repeat_bp = ref_repeats * motif_len
    fixed_bp = max(0, product_size - repeat_bp)

    min_rep = max(1, min_repeats)
    max_rep = max(min_rep, max_repeats)

    min_size = fixed_bp + min_rep * motif_len
    max_size = fixed_bp + max_rep * motif_len
    ref_size = fixed_bp + ref_repeats * motif_len

    compatible = []
    for p in PLATFORMS:
        if p["min_amplicon"] <= min_size and p["max_amplicon"] >= max_size:
            compatible.append(p["name"])

    recommended = compatible[0] if compatible else "Long-read sequencing required"

    return {
        "ssr_id": primer.get("ssr_id"),
        "motif": primer.get("motif", ""),
        "motif_length": motif_len,
        "ref_repeat_count": ref_repeats,
        "ref_product_size": ref_size,
        "fixed_flanking_bp": fixed_bp,
        "min_allele_size": min_size,
        "max_allele_size": max_size,
        "allele_size_range": max_size - min_size,
        "compatible_platforms": compatible,
        "recommended_platform": recommended,
    }


def summarise_pool(
    pass_primers: List[Dict[str, Any]],
    min_repeats: int = 3,
    max_repeats: int = 60,
    platform_name: str = "Illumina MiSeq 300bp PE",
    min_amplicon: int = None,
    max_amplicon: int = None,
) -> Dict[str, Any]:
    """
    Summarise allele size requirements across a full primer pool.
    
    If min_amplicon or max_amplicon are provided, they override platform defaults.
    """
    if not pass_primers:
        return {
            "per_locus": [],
            "pool_max_size": 0,
            "pool_min_size": 0,
            "recommended_platform": "—",
            "platform_breakdown": {},
            "n_loci": 0,
            "n_too_short": 0,
            "n_too_long": 0,
            "filtered_loci": [],
            "platform_min": 0,
            "platform_max": 0,
        }

    per_locus = [calculate_allele_sizes(p, min_repeats, max_repeats) for p in pass_primers]

    # Determine size limits
    platform = next((p for p in PLATFORMS if p["name"] == platform_name), PLATFORMS[1])
    min_allowed = min_amplicon if min_amplicon is not None else platform["min_amplicon"]
    max_allowed = max_amplicon if max_amplicon is not None else platform["max_amplicon"]

    # Flag loci outside size limits
    too_short = []
    too_long = []
    for r in per_locus:
        if r["min_allele_size"] < min_allowed:
            too_short.append(r["ssr_id"])
        if r["max_allele_size"] > max_allowed:
            too_long.append(r["ssr_id"])

    filtered_out = set(too_short) | set(too_long)
    valid_loci = [r for r in per_locus if r["ssr_id"] not in filtered_out]

    pool_max = max((r["max_allele_size"] for r in valid_loci), default=0)
    pool_min = min((r["min_allele_size"] for r in valid_loci), default=0)

    platform_counts = {p["name"]: 0 for p in PLATFORMS}
    for r in valid_loci:
        for pname in r["compatible_platforms"]:
            if pname in platform_counts:
                platform_counts[pname] += 1

    return {
        "per_locus": per_locus,
        "pool_max_size": pool_max,
        "pool_min_size": pool_min,
        "recommended_platform": platform_name if valid_loci else "Long-read sequencing required",
        "platform_breakdown": {k: v for k, v in platform_counts.items() if v > 0},
        "n_loci": len(per_locus),
        "n_too_short": len(too_short),
        "n_too_long": len(too_long),
        "filtered_loci": list(filtered_out),
        "platform_min": min_allowed,
        "platform_max": max_allowed,
    }