"""
capillary_multiplex.py
----------------------
Advanced dye assignment and size binning for capillary multiplex panels.
Fills panels sequentially, maximizing density and using all dyes.
"""

from typing import List, Dict, Any, Tuple, Optional

DYE_SETS = {
    "1 Dye (FAM)": ["FAM"],
    "2 Dyes (FAM, VIC)": ["FAM", "VIC"],
    "3 Dyes (FAM, VIC, NED)": ["FAM", "VIC", "NED"],
    "4 Dyes (FAM, VIC, NED, PET)": ["FAM", "VIC", "NED", "PET"],
}


def _overlaps(bin1: Tuple[int, int], bin2: Tuple[int, int], buffer: int = 0) -> bool:
    lo1, hi1 = bin1
    lo2, hi2 = bin2
    return not (hi1 + buffer < lo2 or hi2 + buffer < lo1)


def assign_dyes_advanced(
    loci: List[Dict[str, Any]],
    dye_set: List[str],
    min_bin_spacing: int = 10,
    max_loci_per_dye: Optional[int] = None,
    allow_overlap_within_dye: bool = False,
    overlap_tolerance: int = 0,
) -> Tuple[List[Dict[str, Any]], List[str], Dict[str, List[Tuple[int, int]]]]:
    """
    Assign loci to dyes within a single panel.
    Returns (assignments, unassigned, dye_bins).
    """
    sorted_loci = sorted(
        loci,
        key=lambda x: (x["max_allele_size"] - x["min_allele_size"]),
        reverse=True
    )
    dye_bins = {dye: [] for dye in dye_set}
    unassigned = []

    for locus in sorted_loci:
        loc_min = locus["min_allele_size"]
        loc_max = locus["max_allele_size"]
        ssr_id = locus["ssr_id"]

        candidates = []
        for dye in dye_set:
            if max_loci_per_dye and len(dye_bins[dye]) >= max_loci_per_dye:
                continue
            conflict = False
            for (b_min, b_max, _) in dye_bins[dye]:
                if allow_overlap_within_dye:
                    overlap = min(loc_max, b_max) - max(loc_min, b_min)
                    if overlap > overlap_tolerance:
                        conflict = True
                        break
                else:
                    if _overlaps((loc_min, loc_max), (b_min, b_max), min_bin_spacing):
                        conflict = True
                        break
            if not conflict:
                candidates.append(dye)

        if candidates:
            best_dye = min(candidates, key=lambda d: len(dye_bins[d]))
            dye_bins[best_dye].append((loc_min, loc_max, ssr_id))
        else:
            unassigned.append(ssr_id)

    assignments = []
    for dye in dye_set:
        for (b_min, b_max, ssr_id) in dye_bins[dye]:
            assignments.append({
                "ssr_id": ssr_id,
                "dye": dye,
                "min_size": b_min,
                "max_size": b_max,
            })

    return assignments, unassigned, {dye: [(bmin, bmax) for (bmin, bmax, _) in bins] for dye, bins in dye_bins.items()}


def assign_to_panels_sequential(
    loci: List[Dict[str, Any]],
    dye_set: List[str],
    max_panels: int = 1,
    min_bin_spacing: int = 5,
    max_loci_per_dye: Optional[int] = None,
    allow_overlap: bool = False,
    overlap_tolerance: int = 0,
) -> Dict[str, Any]:
    """
    Fill Panel 1 completely, then Panel 2, etc.
    """
    remaining_loci = list(loci)  # copy
    panels_assignments = []
    panels_dye_bins = []
    all_unassigned = []

    for panel_idx in range(1, max_panels + 1):
        if not remaining_loci:
            break

        assignments, unassigned, dye_bins = assign_dyes_advanced(
            remaining_loci, dye_set,
            min_bin_spacing=min_bin_spacing,
            max_loci_per_dye=max_loci_per_dye,
            allow_overlap_within_dye=allow_overlap,
            overlap_tolerance=overlap_tolerance,
        )

        if not assignments:
            # No loci could be placed in this panel; stop trying further panels
            all_unassigned.extend([loc["ssr_id"] for loc in remaining_loci])
            break

        # Tag assignments with panel number
        for a in assignments:
            a["panel"] = panel_idx
        panels_assignments.append(assignments)
        panels_dye_bins.append(dye_bins)

        # Remove assigned loci from remaining
        assigned_ids = {a["ssr_id"] for a in assignments}
        remaining_loci = [loc for loc in remaining_loci if loc["ssr_id"] not in assigned_ids]
        all_unassigned.extend(unassigned)

    # Build result
    all_assignments = []
    panels_dict = {}
    for i, assigns in enumerate(panels_assignments, start=1):
        panels_dict[i] = assigns
        all_assignments.extend(assigns)

    return {
        "assignments": all_assignments,
        "panels": panels_dict,
        "n_panels": len(panels_assignments),
        "n_total": len(loci),
        "n_assigned": len(all_assignments),
        "unassigned": all_unassigned,
    }


def suggest_best_panel(
    loci: List[Dict[str, Any]],
    dye_set: List[str],
    min_spacing_options: List[int] = None,
    max_per_dye_options: List[int] = None,
) -> Dict[str, Any]:
    if min_spacing_options is None:
        min_spacing_options = [5, 8, 10, 12, 15]
    if max_per_dye_options is None:
        max_per_dye_options = [10, 15, 20, 25, 0]

    best_result = None
    best_n_assigned = -1

    for spacing in min_spacing_options:
        for max_per_dye in max_per_dye_options:
            assignments, unassigned, _ = assign_dyes_advanced(
                loci=loci,
                dye_set=dye_set,
                min_bin_spacing=spacing,
                max_loci_per_dye=max_per_dye if max_per_dye > 0 else None,
                allow_overlap_within_dye=False,
            )
            if len(assignments) > best_n_assigned:
                best_n_assigned = len(assignments)
                best_result = {
                    "assignments": assignments,
                    "unassigned": unassigned,
                    "n_assigned": len(assignments),
                    "n_total": len(loci),
                    "best_params": {"min_bin_spacing": spacing, "max_loci_per_dye": max_per_dye},
                }

    return best_result


def filter_unique_loci(primers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if not primers:
        return []
    import pandas as pd
    df = pd.DataFrame(primers)
    if "pair_rank" not in df.columns:
        return primers
    df = df.sort_values("pair_rank").drop_duplicates(subset=["ssr_id"], keep="first")
    return df.to_dict(orient="records")