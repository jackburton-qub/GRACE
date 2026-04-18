"""
gbs_multiplex.py
----------------
Multiplex primer pool compatibility checker for GBS workflows.

Supports tagged primers (adapters) for accurate dimer prediction.
"""

import itertools
from typing import List, Dict, Any, Optional, Set
from collections import defaultdict

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


def _score_3prime(seq_a: str, seq_b: str, n: int = 6) -> int:
    """Score 3' end complementarity."""
    end_a = seq_a[-n:].upper()
    end_b = _rc(seq_b[-n:])
    score = 0
    for a, b in zip(reversed(end_a), reversed(end_b)):
        if a == b:
            score += 1
        else:
            break
    return score


def check_pool_compatibility(
    primers: List[Dict[str, Any]],
    max_3prime_score: int = 3,
    max_tm_spread: float = 5.0,
    size_conflict_window: int = 2,
    tag_config: Optional[Any] = None,   # TagConfig object
) -> Dict[str, Any]:
    """
    Check a list of primer pairs for multiplex pool compatibility.

    If tag_config is provided, uses tagged sequences for dimer checking.
    """
    n = len(primers)
    if n == 0:
        return _empty_result()

    # Apply tags if configured
    if tag_config and tag_config.include_in_dimer_check:
        working_primers = [tag_config.apply_tags(p) for p in primers]
        left_key = "left_primer_tagged"
        right_key = "right_primer_tagged"
    else:
        working_primers = primers
        left_key = "left_primer"
        right_key = "right_primer"

    # Build flat list of all individual primers
    all_primers = []
    for p in working_primers:
        all_primers.append({
            "ssr_id": p["ssr_id"],
            "direction": "F",
            "seq": p[left_key].upper(),
        })
        all_primers.append({
            "ssr_id": p["ssr_id"],
            "direction": "R",
            "seq": p[right_key].upper(),
        })

    # 3' dimer check
    dimer_risks = []
    for a, b in itertools.combinations(all_primers, 2):
        if a["ssr_id"] == b["ssr_id"]:
            continue
        score = _score_3prime(a["seq"], b["seq"])
        if score >= max_3prime_score:
            dimer_risks.append({
                "ssr_id_a": a["ssr_id"],
                "dir_a": a["direction"],
                "ssr_id_b": b["ssr_id"],
                "dir_b": b["direction"],
                "score": score,
            })

    # Size conflict check (using original primers for product size)
    size_conflicts = []
    size_index = {}
    for p in primers:
        sz = p.get("product_size", 0)
        for existing_id, existing_sz in size_index.items():
            if abs(sz - existing_sz) <= size_conflict_window:
                size_conflicts.append({
                    "ssr_id_a": existing_id,
                    "ssr_id_b": p["ssr_id"],
                    "size_a": existing_sz,
                    "size_b": sz,
                })
        size_index[p["ssr_id"]] = sz

    # Tm analysis
    all_tms = [p["left_tm"] for p in primers] + [p["right_tm"] for p in primers]
    tm_min = min(all_tms)
    tm_max = max(all_tms)
    tm_spread = round(tm_max - tm_min, 2)

    # Flagged IDs
    flagged_ids = set()
    for r in dimer_risks:
        flagged_ids.add(r["ssr_id_a"])
        flagged_ids.add(r["ssr_id_b"])
    for c in size_conflicts:
        flagged_ids.add(c["ssr_id_a"])
        flagged_ids.add(c["ssr_id_b"])

    n_clean = n - len(flagged_ids)

    # Compatibility score
    max_pairs = n * (n - 1)
    dimer_pen = len(dimer_risks) / max(1, max_pairs)
    conflict_pen = len(size_conflicts) / max(1, n)
    tm_pen = max(0.0, (tm_spread - max_tm_spread) / 20.0)
    score = max(0.0, 1.0 - dimer_pen * 2 - conflict_pen - tm_pen)

    # Recommendation
    issues = []
    if dimer_risks:
        issues.append(f"{len(dimer_risks)} dimer risk(s)")
    if size_conflicts:
        issues.append(f"{len(size_conflicts)} size conflict(s)")
    if tm_spread > max_tm_spread:
        issues.append(f"Tm spread {tm_spread:.1f}°C")

    if not issues:
        recommendation = f"Pool of {n} pairs looks compatible. Tm spread: {tm_spread:.1f}°C."
    else:
        recommendation = f"{len(flagged_ids)} of {n} pairs flagged: {'; '.join(issues)}."

    return {
        "dimer_risks": dimer_risks,
        "size_conflicts": size_conflicts,
        "tm_spread": tm_spread,
        "tm_min": round(tm_min, 2),
        "tm_max": round(tm_max, 2),
        "flagged_ids": flagged_ids,
        "n_flagged": len(flagged_ids),
        "n_clean": n_clean,
        "n_total": n,
        "compatibility_score": round(score, 3),
        "recommendation": recommendation,
    }


def _empty_result() -> dict:
    return {
        "dimer_risks": [], "size_conflicts": [],
        "tm_spread": 0.0, "tm_min": 0.0, "tm_max": 0.0,
        "flagged_ids": set(), "n_flagged": 0, "n_clean": 0,
        "n_total": 0, "compatibility_score": 1.0,
        "recommendation": "No primers to check.",
    }


def optimise_pool(
    primers: List[Dict[str, Any]],
    max_3prime_score: int = 3,
    size_conflict_window: int = 2,
    max_tm_spread: float = 5.0,
    tag_config: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Find the largest compatible subset of primer pairs using greedy conflict removal.
    """
    if not primers:
        return {
            "recommended_pool": [],
            "removed_ids": set(),
            "removal_reasons": {},
            "n_kept": 0,
            "n_removed": 0,
            "final_tm_spread": 0.0,
            "summary": "No primers to optimise.",
        }

    primer_by_id = {p["ssr_id"]: p for p in primers}
    active_ids = set(primer_by_id.keys())

    # Apply tags for dimer check
    if tag_config and tag_config.include_in_dimer_check:
        working_primers = [tag_config.apply_tags(p) for p in primers]
        left_key = "left_primer_tagged"
        right_key = "right_primer_tagged"
    else:
        working_primers = primers
        left_key = "left_primer"
        right_key = "right_primer"

    # Build conflict adjacency
    conflicts: Dict[int, set] = {sid: set() for sid in active_ids}
    removal_reasons: Dict[int, str] = {}

    all_seqs = []
    for p in working_primers:
        all_seqs.append((p["ssr_id"], "F", p[left_key].upper()))
        all_seqs.append((p["ssr_id"], "R", p[right_key].upper()))

    for i in range(len(all_seqs)):
        for j in range(i + 1, len(all_seqs)):
            id_a, dir_a, seq_a = all_seqs[i]
            id_b, dir_b, seq_b = all_seqs[j]
            if id_a == id_b:
                continue
            score = _score_3prime(seq_a, seq_b)
            if score >= max_3prime_score:
                conflicts[id_a].add(id_b)
                conflicts[id_b].add(id_a)

    # Size conflicts
    size_map = {}
    for p in primers:
        sz = p.get("product_size", 0)
        for existing_id, existing_sz in size_map.items():
            if abs(sz - existing_sz) <= size_conflict_window:
                conflicts[p["ssr_id"]].add(existing_id)
                conflicts[existing_id].add(p["ssr_id"])
        size_map[p["ssr_id"]] = sz

    # Greedy removal
    removed_ids = set()
    while True:
        max_conflicts = 0
        worst_id = None
        for sid in active_ids:
            n = len(conflicts[sid] & active_ids)
            if n > max_conflicts:
                max_conflicts = n
                worst_id = sid
            elif n == max_conflicts and worst_id is not None:
                if primer_by_id[sid].get("product_size", 0) < primer_by_id[worst_id].get("product_size", 0):
                    worst_id = sid
        if max_conflicts == 0:
            break
        active_ids.discard(worst_id)
        removed_ids.add(worst_id)
        removal_reasons[worst_id] = f"Conflicted with {max_conflicts} other pair(s)"

    # Tm spread trim
    remaining = [primer_by_id[sid] for sid in active_ids]
    if len(remaining) > 1:
        import statistics
        all_tms = [p["left_tm"] for p in remaining] + [p["right_tm"] for p in remaining]
        tm_spread = max(all_tms) - min(all_tms)
        if tm_spread > max_tm_spread:
            while len(remaining) > 1:
                tms = [p["left_tm"] for p in remaining] + [p["right_tm"] for p in remaining]
                median_tm = statistics.median(tms)
                tm_spread = max(tms) - min(tms)
                if tm_spread <= max_tm_spread:
                    break
                worst_p = max(remaining, key=lambda p: max(abs(p["left_tm"]-median_tm), abs(p["right_tm"]-median_tm)))
                remaining.remove(worst_p)
                active_ids.discard(worst_p["ssr_id"])
                removed_ids.add(worst_p["ssr_id"])
                removal_reasons[worst_p["ssr_id"]] = f"Tm outlier"

    recommended_pool = [primer_by_id[sid] for sid in active_ids]
    n_kept = len(recommended_pool)
    n_removed = len(removed_ids)

    final_tms = [p["left_tm"] for p in recommended_pool] + [p["right_tm"] for p in recommended_pool] if recommended_pool else [0.0]
    final_tm_spread = round(max(final_tms) - min(final_tms), 2)

    summary = f"Recommended pool: {n_kept:,} pairs ({n_removed:,} removed). Final Tm spread: {final_tm_spread:.1f}°C."

    return {
        "recommended_pool": recommended_pool,
        "removed_ids": removed_ids,
        "removal_reasons": removal_reasons,
        "n_kept": n_kept,
        "n_removed": n_removed,
        "final_tm_spread": final_tm_spread,
        "summary": summary,
    }