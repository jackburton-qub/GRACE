"""
gbs_tags.py
-----------
Handle adapter tags and barcodes for GBS primer design.

GBS primers typically include 5' adapter tails for Illumina sequencing.
These tags affect dimer formation and must be included in multiplex checks.
"""

from typing import List, Dict, Any, Optional

# Common Illumina adapter sequences
ILLUMINA_P5 = "AATGATACGGCGACCACCGAGATCTACAC"
ILLUMINA_P7 = "CAAGCAGAAGACGGCATACGAGAT"

# Standard Illumina TruSeq universal adapter (for single-index)
TRUSEQ_UNIVERSAL = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

# Standard indexing primer (reverse)
TRUSEQ_INDEX = "CAAGCAGAAGACGGCATACGAGAT[index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

# Common inline barcode lengths
BARCODE_LENGTHS = [4, 5, 6, 8, 10]


class TagConfig:
    """Container for tag configuration."""
    def __init__(
        self,
        forward_tag: str = "",
        reverse_tag: str = "",
        include_in_dimer_check: bool = True,
        include_in_amplicon_export: bool = False,
    ):
        self.forward_tag = forward_tag.upper() if forward_tag else ""
        self.reverse_tag = reverse_tag.upper() if reverse_tag else ""
        self.include_in_dimer_check = include_in_dimer_check
        self.include_in_amplicon_export = include_in_amplicon_export

    @property
    def has_tags(self) -> bool:
        return bool(self.forward_tag or self.reverse_tag)

    def apply_tags(self, primer: Dict[str, Any]) -> Dict[str, Any]:
        """
        Return a copy of the primer dict with tagged sequences.
        Tags are added to 5' ends: tag + primer sequence.
        """
        tagged = primer.copy()
        if self.forward_tag:
            tagged["left_primer_tagged"] = self.forward_tag + primer.get("left_primer", "")
        else:
            tagged["left_primer_tagged"] = primer.get("left_primer", "")

        if self.reverse_tag:
            tagged["right_primer_tagged"] = self.reverse_tag + primer.get("right_primer", "")
        else:
            tagged["right_primer_tagged"] = primer.get("right_primer", "")

        return tagged

    def get_dimer_check_sequences(self, primer: Dict[str, Any]) -> tuple:
        """Return (left_seq, right_seq) to use for dimer checking."""
        if self.include_in_dimer_check and self.has_tags:
            left = primer.get("left_primer_tagged", primer.get("left_primer", ""))
            right = primer.get("right_primer_tagged", primer.get("right_primer", ""))
        else:
            left = primer.get("left_primer", "")
            right = primer.get("right_primer", "")
        return left, right


# Presets for common configurations
PRESETS = {
    "Illumina TruSeq (single index)": {
        "forward_tag": TRUSEQ_UNIVERSAL,
        "reverse_tag": TRUSEQ_INDEX.replace("[index]", "NNNNNN"),
        "description": "Standard Illumina single-index adapters",
    },
    "Illumina Nextera": {
        "forward_tag": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        "reverse_tag": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
        "description": "Nextera transposase adapters",
    },
    "Custom (no tags)": {
        "forward_tag": "",
        "reverse_tag": "",
        "description": "Primers only, no adapter tails",
    },
}


def validate_barcode_set(
    barcodes: List[str],
    min_hamming: int = 3,
    max_homopolymer: int = 2,
    balanced_gc: bool = True,
) -> Dict[str, Any]:
    """
    Validate a set of sample barcodes for multiplexing.

    Args:
        barcodes: list of barcode sequences
        min_hamming: minimum Hamming distance between any pair
        max_homopolymer: maximum allowed consecutive identical bases
        balanced_gc: if True, warn about GC extremes (<30% or >70%)

    Returns:
        dict with validation results
    """
    n = len(barcodes)
    issues = []
    warnings = []
    conflicts = []

    # Check each barcode individually
    for bc in barcodes:
        # Homopolymer runs
        run = 1
        max_run = 1
        for i in range(1, len(bc)):
            if bc[i] == bc[i-1]:
                run += 1
                max_run = max(max_run, run)
            else:
                run = 1
        if max_run > max_homopolymer:
            issues.append(f"Barcode {bc}: homopolymer run of {max_run} (>{max_homopolymer})")

        # GC content
        if balanced_gc:
            gc = (bc.count('G') + bc.count('C')) / len(bc) * 100
            if gc < 30:
                warnings.append(f"Barcode {bc}: low GC ({gc:.1f}%)")
            elif gc > 70:
                warnings.append(f"Barcode {bc}: high GC ({gc:.1f}%)")

    # Check pairwise Hamming distances
    for i in range(n):
        for j in range(i+1, n):
            bc1, bc2 = barcodes[i], barcodes[j]
            if len(bc1) != len(bc2):
                warnings.append(f"Barcodes have different lengths: {bc1} ({len(bc1)}) vs {bc2} ({len(bc2)})")
            dist = sum(1 for a, b in zip(bc1, bc2) if a != b)
            if dist < min_hamming:
                conflicts.append({
                    "barcode_a": bc1,
                    "barcode_b": bc2,
                    "hamming_distance": dist,
                })

    passed = len(issues) == 0 and len(conflicts) == 0

    return {
        "passed": passed,
        "n_barcodes": n,
        "issues": issues,
        "warnings": warnings,
        "conflicts": conflicts,
        "recommendation": (
            "Barcode set is valid." if passed
            else f"Found {len(issues)} issues and {len(conflicts)} conflicts. Review above."
        ),
    }