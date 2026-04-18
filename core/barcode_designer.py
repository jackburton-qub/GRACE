"""
barcode_designer.py
-------------------
Generate balanced barcode sets for sample multiplexing.
"""

import random
from typing import List, Dict, Any


def hamming_distance(seq1: str, seq2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def validate_barcode_set(barcodes: List[str], min_distance: int = 3) -> Dict[str, Any]:
    """Check a set of barcodes for conflicts."""
    issues = []
    for i, bc1 in enumerate(barcodes):
        # Homopolymer check
        for base in set(bc1):
            if base * 4 in bc1:
                issues.append(f"Barcode {bc1} has homopolymer run")
        # GC check
        gc = (bc1.count('G') + bc1.count('C')) / len(bc1) * 100
        if gc < 30 or gc > 70:
            issues.append(f"Barcode {bc1} GC% = {gc:.1f} (outside 30-70%)")

    conflicts = []
    for i, bc1 in enumerate(barcodes):
        for bc2 in barcodes[i+1:]:
            if hamming_distance(bc1, bc2) < min_distance:
                conflicts.append((bc1, bc2))

    return {
        "valid": len(issues) == 0 and len(conflicts) == 0,
        "issues": issues,
        "conflicts": conflicts,
    }


def generate_barcodes(n: int, length: int = 8, min_distance: int = 3) -> List[str]:
    """
    Generate a set of n barcodes of given length with minimum Hamming distance.
    Simple greedy algorithm.
    """
    bases = ['A', 'C', 'G', 'T']
    barcodes = []
    attempts = 0
    while len(barcodes) < n and attempts < 10000:
        candidate = ''.join(random.choices(bases, k=length))
        # Basic filters
        if any(base * 3 in candidate for base in bases):
            continue
        gc = (candidate.count('G') + candidate.count('C')) / length * 100
        if gc < 40 or gc > 60:
            continue
        # Check distance
        if all(hamming_distance(candidate, bc) >= min_distance for bc in barcodes):
            barcodes.append(candidate)
        attempts += 1
    return barcodes