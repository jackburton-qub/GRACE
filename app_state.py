"""
app_state.py
------------
Single dataclass holding all application data.
Lives for the entire session — no serialisation, no stale keys.
All panels read from and write to this object directly.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class AppState:
    # ── Genome ────────────────────────────────────────────
    genome:          Optional[dict]  = None   # {contig_name: sequence}
    genome_path:     Optional[str]   = None   # path to FASTA file on disk
    genome_filename: Optional[str]   = None

    # ── SSR Detection ─────────────────────────────────────
    ssrs:               Optional[list] = None  # list of SSR dicts
    selected_ssrs:      Optional[list] = None  # user-selected subset
    ssr_detection_time: Optional[float] = None

    # ── Primer Design ─────────────────────────────────────
    primer_results:          Optional[list] = None  # all designed primers
    filtered_primer_results: Optional[list] = None  # after quality filters
    product_min:             int = 100
    product_max:             int = 250

    # ── Specificity ───────────────────────────────────────
    specificity_results: Optional[list] = None
    blast_raw_rows:      Optional[list] = None
    # blast_bin_dir is no longer stored here — BLAST binaries are bundled
    # with the app and resolved automatically by get_blast_bin_dir() in
    # core/primer_specificity_blast.py.

    # ── Run version stamps ────────────────────────────────
    # Incremented each time a step completes — downstream panels
    # can check these to detect stale data before running
    ssr_version:     int = 0
    primer_version:  int = 0
    blast_version:   int = 0

    def clear_downstream_of_genome(self):
        """Call when a new genome is loaded."""
        self.ssrs               = None
        self.selected_ssrs      = None
        self.ssr_detection_time = None
        self.ssr_version        = 0
        self.clear_downstream_of_ssrs()

    def clear_downstream_of_ssrs(self):
        """Call when new SSRs are detected."""
        self.primer_results          = None
        self.filtered_primer_results = None
        self.primer_version          = 0
        self.clear_downstream_of_primers()

    def clear_downstream_of_primers(self):
        """Call when new primers are designed."""
        self.specificity_results = None
        self.blast_raw_rows      = None
        self.blast_version       = 0

    @property
    def has_genome(self) -> bool:
        return self.genome is not None

    @property
    def has_ssrs(self) -> bool:
        return self.ssrs is not None and len(self.ssrs) > 0

    @property
    def has_primers(self) -> bool:
        return self.primer_results is not None and len(self.primer_results) > 0

    @property
    def has_specificity(self) -> bool:
        return self.specificity_results is not None and len(self.specificity_results) > 0

    @property
    def blast_primers(self) -> list:
        """Returns filtered primers if available, otherwise full set."""
        if self.filtered_primer_results:
            return self.filtered_primer_results
        return self.primer_results or []