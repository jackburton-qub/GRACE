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

    # ── GFF Annotation (optional) ─────────────────────────
    gff_path:        Optional[str]    = None   # path to GFF/GTF file on disk
    gff_filename:    Optional[str]    = None
    gff_features:    Optional[object] = None   # parsed GFFIndex (built lazily)
    chrom_names:     Optional[dict]   = None   # contig accession → chromosome display name

    # ── SSR Detection ─────────────────────────────────────
    ssrs:               Optional[list]  = None  # list of SSR dicts
    selected_ssrs:      Optional[list]  = None  # user-selected subset
    ssr_detection_time: Optional[float] = None

    # ── Primer Design ─────────────────────────────────────
    primer_results:          Optional[list] = None  # all designed primers
    filtered_primer_results: Optional[list] = None  # after quality filters
    product_min:             int  = 100
    product_max:             int  = 250

    # Workflow mode: "capillary", "amplicon", or "gbs_re" (future)
    workflow_mode:           str  = "capillary"

    # ── Specificity ───────────────────────────────────────
    specificity_results: Optional[list] = None
    blast_raw_rows:      Optional[list] = None

    # ── Run version stamps ────────────────────────────────
    ssr_version:     int = 0
    primer_version:  int = 0
    blast_version:   int = 0

    # ---------------------------------------------------------
    # Clear helpers
    # ---------------------------------------------------------

    def clear_downstream_of_genome(self):
        """Call when a new genome is loaded."""
        self.ssrs               = None
        self.selected_ssrs      = None
        self.ssr_detection_time = None
        self.ssr_version        = 0
        self.clear_downstream_of_ssrs()

    def clear_gff(self):
        """Call when GFF is cleared or a new one loaded."""
        self.gff_path     = None
        self.gff_filename = None
        self.gff_features = None
        self.chrom_names  = None

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

    # ---------------------------------------------------------
    # Properties
    # ---------------------------------------------------------

    @property
    def has_genome(self) -> bool:
        return self.genome is not None

    @property
    def has_gff(self) -> bool:
        return self.gff_path is not None

    @property
    def has_chrom_names(self) -> bool:
        return bool(self.chrom_names)

    @property
    def has_ssrs(self) -> bool:
        return self.ssrs is not None and len(self.ssrs) > 0

    @property
    def has_primers(self) -> bool:
        return self.primer_results is not None and len(self.primer_results) > 0

    @property
    def has_specificity(self) -> bool:
        return (self.specificity_results is not None and
                len(self.specificity_results) > 0)

    @property
    def blast_primers(self) -> list:
        """Returns filtered primers if available, otherwise full set."""
        if self.filtered_primer_results:
            return self.filtered_primer_results
        return self.primer_results or []