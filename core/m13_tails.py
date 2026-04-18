"""
m13_tails.py
------------
M13 universal tail sequences for capillary electrophoresis.
"""

# Standard M13 tails
M13_FORWARD = "TGTAAAACGACGGCCAGT"
M13_REVERSE = "CAGGAAACAGCTATGACC"

# PIG-tail (added to 5' of reverse primer to promote non-templated adenylation)
PIG_TAIL = "GTTTCTT"

def add_m13_tails(primer: dict) -> dict:
    """Return a copy of the primer dict with M13 tails added to 5' ends."""
    tailed = primer.copy()
    tailed["left_primer_tailed"] = M13_FORWARD + primer.get("left_primer", "")
    tailed["right_primer_tailed"] = PIG_TAIL + M13_REVERSE + primer.get("right_primer", "")
    return tailed