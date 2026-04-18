"""
ld_filter.py
------------
Linkage disequilibrium filter: thin markers by physical distance.
"""

from typing import List, Dict, Any


def thin_markers_by_distance(
    markers: List[Dict[str, Any]],
    min_distance_bp: int,
    contig_key: str = "contig",
    position_key: str = "start",
) -> List[Dict[str, Any]]:
    """
    Return a subset of markers where no two on the same contig are closer than min_distance_bp.

    Args:
        markers: list of marker dicts (must contain contig and position)
        min_distance_bp: minimum allowed distance between markers on same contig
        contig_key: dict key for contig/chromosome name
        position_key: dict key for marker position (1‑based)

    Returns:
        Thinned list of markers.
    """
    if not markers or min_distance_bp <= 0:
        return markers

    # Group by contig
    by_contig: Dict[str, List[Dict]] = {}
    for m in markers:
        contig = m.get(contig_key)
        if contig is None:
            continue
        by_contig.setdefault(contig, []).append(m)

    kept = []
    for contig, contig_markers in by_contig.items():
        # Sort by position
        sorted_markers = sorted(contig_markers, key=lambda x: x.get(position_key, 0))
        last_kept_pos = -min_distance_bp  # ensure first marker is kept
        for m in sorted_markers:
            pos = m.get(position_key, 0)
            if pos - last_kept_pos >= min_distance_bp:
                kept.append(m)
                last_kept_pos = pos

    # Add any markers that lacked contig info (keep them all)
    for m in markers:
        if m.get(contig_key) is None:
            kept.append(m)

    return kept