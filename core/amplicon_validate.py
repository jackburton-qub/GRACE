"""
amplicon_validate.py
--------------------
Validation orchestrator for amplicon sequencing panels.
No restriction enzyme check. No size conflict filter.
"""
from typing import List, Dict, Any, Optional
from core.multiplex_utils import check_pool_compatibility, optimise_pool
from core.amplicon_sizing import summarise_pool
from core.adapter_tags import TagConfig

def validate_panel(
    primers: List[Dict[str, Any]],
    genome: Dict[str, str],
    platform_name: str,
    tag_config: Optional[TagConfig] = None,
    min_repeats: int = 3,
    max_repeats: int = 60,
    max_3prime_score: int = 3,
    max_tm_spread: float = 5.0,
    min_amplicon: int = None,
    max_amplicon: int = None,
    auto_optimise: bool = False,
) -> Dict[str, Any]:
    if not primers:
        return {"error": "No primers provided", "n_input": 0}

    n_total = len(primers)

    # 1. Size/platform check
    size_result = summarise_pool(
        primers,
        min_repeats=min_repeats,
        max_repeats=max_repeats,
        platform_name=platform_name,
        min_amplicon=min_amplicon,
        max_amplicon=max_amplicon,
    )
    size_failed = set(size_result["filtered_loci"])

    # 2. Multiplex check (dimer + Tm only)
    #    We call the original function but then discard size conflict flags.
    multiplex_result = check_pool_compatibility(
        primers,
        max_3prime_score=max_3prime_score,
        max_tm_spread=max_tm_spread,
        size_conflict_window=2,   # value ignored because we override flags
        tag_config=tag_config,
    )

    # Override flagged_ids to ONLY include dimer risks
    dimer_flagged = set()
    for d in multiplex_result["dimer_risks"]:
        dimer_flagged.add(d["ssr_id_a"])
        dimer_flagged.add(d["ssr_id_b"])
    multiplex_result["flagged_ids"] = dimer_flagged

    # Combine flags
    all_flagged = size_failed | dimer_flagged
    clean_primers = [p for p in primers if p["ssr_id"] not in all_flagged]

    # Optimise if requested
    optimise_result = None
    final_pool = clean_primers
    if auto_optimise and clean_primers:
        optimise_result = optimise_pool(
            clean_primers,
            max_3prime_score=max_3prime_score,
            size_conflict_window=10000,  # effectively disable in optimiser too
            max_tm_spread=max_tm_spread,
            tag_config=tag_config,
        )
        # Also override any size conflict flags from optimiser output
        final_pool = optimise_result["recommended_pool"]

    issues_summary = {
        "dimer_risks": len(multiplex_result["dimer_risks"]),
        "size_filtered_short": size_result["n_too_short"],
        "size_filtered_long": size_result["n_too_long"],
    }

    return {
        "n_input": n_total,
        "n_clean": len(clean_primers),
        "n_final": len(final_pool),
        "issues_summary": issues_summary,
        "multiplex_result": multiplex_result,
        "size_result": size_result,
        "optimise_result": optimise_result,
        "final_pool": final_pool,
        "clean_pool": clean_primers,
        "flagged_ids": list(all_flagged),
        "compatibility_score": multiplex_result["compatibility_score"],
        "platform": platform_name,
    }