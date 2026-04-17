"""
gff_parser.py
-------------
Parses GFF3 / GTF annotation files and builds a fast interval index
for classifying genomic positions as:

    exon        — overlaps an exon feature
    CDS         — within a CDS but not an exon
    intron      — within a gene but no exon/CDS overlap
    intergenic  — outside all annotated gene features

Also extracts chromosome name mappings from region features so that
contig accessions (e.g. CP007106.1) can be displayed as chromosome
names (e.g. 2L, chrX) in visualisations.

Performance
-----------
Uses a sorted interval list + binary search (bisect) per contig.
No third-party dependencies required.

Feature priority (highest wins):
    exon > CDS > intron > intergenic
"""

import bisect
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

_Interval = Tuple[int, int, str]
_PRIORITY = {"exon": 0, "CDS": 1, "intron": 2, "intergenic": 3}


class GFFIndex:
    """
    Fast interval index for genomic feature lookup.
    Also stores chromosome name mappings extracted from region features.
    """

    def __init__(self):
        self._intervals: Dict[str, List[_Interval]] = defaultdict(list)
        self._starts:    Dict[str, List[int]]        = defaultdict(list)
        # Maps contig accession → chromosome name (e.g. "CP007106.1" → "2L")
        # Only populated when GFF has explicit chromosome annotations
        self.chrom_names: Dict[str, str] = {}
        self.n_features: int = 0
        self.n_contigs:  int = 0

    def _add(self, contig: str, start: int, end: int, feature: str):
        self._intervals[contig].append((start, end, feature))
        self.n_features += 1

    def build(self):
        for contig in self._intervals:
            self._intervals[contig].sort(key=lambda x: x[0])
            self._starts[contig] = [iv[0] for iv in self._intervals[contig]]
        self.n_contigs = len(self._intervals)

    def classify(self, contig: str, start: int, end: int) -> str:
        """
        Classify a genomic position (1-based inclusive).
        Returns: 'exon', 'CDS', 'intron', or 'intergenic'.
        """
        if contig not in self._intervals:
            return "intergenic"

        q_start = start - 1
        q_end   = end

        intervals = self._intervals[contig]
        starts    = self._starts[contig]

        right = bisect.bisect_right(starts, q_end)
        best  = "intergenic"
        best_priority = _PRIORITY["intergenic"]

        for i in range(right - 1, -1, -1):
            iv_start, iv_end, iv_type = intervals[i]
            if iv_end <= q_start:
                break
            priority = _PRIORITY.get(iv_type, 99)
            if priority < best_priority:
                best          = iv_type
                best_priority = priority
                if best_priority == 0:
                    break

        return best

    def get_display_name(self, contig: str) -> str:
        """Return chromosome name if known, otherwise the contig name."""
        return self.chrom_names.get(contig, contig)


# ---------------------------------------------------------------------------
# Chromosome name extraction helpers
# ---------------------------------------------------------------------------

def _parse_attributes_gff3(attr_str: str) -> dict:
    """Parse GFF3 attribute string into a dict."""
    attrs = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs


def _parse_attributes_gtf(attr_str: str) -> dict:
    """Parse GTF attribute string into a dict."""
    attrs = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        # GTF: key "value" or key value
        m = re.match(r'(\S+)\s+"([^"]*)"', part)
        if m:
            attrs[m.group(1)] = m.group(2)
        else:
            tokens = part.split()
            if len(tokens) >= 2:
                attrs[tokens[0]] = tokens[1]
    return attrs


def _extract_chrom_name_gff3(attrs: dict) -> Optional[str]:
    """
    Extract a human-readable chromosome name from GFF3 region attributes.
    NCBI GFF3 files use patterns like:
      Name=2L  or  chromosome=2L  or  Dbxref=taxon:...;Name=chrX
    """
    # Direct Name attribute (most common in NCBI files)
    name = attrs.get("Name", "")
    if name and _looks_like_chrom_name(name):
        return name

    # chromosome attribute
    chrom = attrs.get("chromosome", "")
    if chrom and _looks_like_chrom_name(chrom):
        return chrom

    # Ensembl uses 'ID=chromosome:2L'
    id_val = attrs.get("ID", "")
    if id_val.startswith("chromosome:"):
        cname = id_val[len("chromosome:"):]
        if _looks_like_chrom_name(cname):
            return cname

    return None


def _looks_like_chrom_name(name: str) -> bool:
    """
    Return True if a name looks like a chromosome name rather than
    an accession number.
    Chromosome names: 2L, 2R, 3L, 3R, X, Y, 4, Chr1, chrX, LG1, MT etc
    Accession numbers: CP007106.1, KI914433.1, CM071426.1 etc
    """
    # Reject if looks like an accession (letters + digits + dot + digit)
    if re.match(r'^[A-Z]{1,4}\d+\.\d+$', name):
        return False
    # Reject very long names (probably IDs not chr names)
    if len(name) > 20:
        return False
    # Accept common chromosome name patterns
    if re.match(r'^(chr|chromosome|chrom|lg|linkage_group|mt|mito|mitochondr)',
                name, re.IGNORECASE):
        return True
    # Accept short alphanumeric names like 2L, 3R, X, Y, 4, I, II, III
    if re.match(r'^(\d{1,3}[LRS]?|[IXVY]{1,6}|[XYZ]\d*)$', name, re.IGNORECASE):
        return True
    return False


# ---------------------------------------------------------------------------
# GFF3 parser
# ---------------------------------------------------------------------------

def _parse_gff3(path: str, index: GFFIndex, progress_callback=None):
    target_features = {"gene", "exon", "CDS"}
    gene_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    n_lines = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid   = parts[0]
            feature = parts[2]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue

            # Extract chromosome names from region features
            if feature == "region":
                attrs   = _parse_attributes_gff3(parts[8])
                genome_type = attrs.get("genome", attrs.get("genome-type", ""))
                if "chromosome" in genome_type.lower() or attrs.get("chromosome"):
                    chrom_name = _extract_chrom_name_gff3(attrs)
                    if chrom_name and chrom_name != seqid:
                        index.chrom_names[seqid] = chrom_name

            elif feature == "gene":
                gene_intervals[seqid].append((start - 1, end))

            elif feature in target_features:
                index._add(seqid, start - 1, end, feature)

            n_lines += 1
            if progress_callback and n_lines % 100_000 == 0:
                progress_callback(n_lines)

    _derive_introns(index, gene_intervals)


# ---------------------------------------------------------------------------
# GTF parser
# ---------------------------------------------------------------------------

def _parse_gtf(path: str, index: GFFIndex, progress_callback=None):
    target_features = {"gene", "exon", "CDS"}
    gene_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    n_lines = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid   = parts[0]
            feature = parts[2]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue

            if feature == "gene":
                gene_intervals[seqid].append((start - 1, end))
            elif feature in target_features:
                index._add(seqid, start - 1, end, feature)

            n_lines += 1
            if progress_callback and n_lines % 100_000 == 0:
                progress_callback(n_lines)

    _derive_introns(index, gene_intervals)


# ---------------------------------------------------------------------------
# Intron derivation
# ---------------------------------------------------------------------------

def _derive_introns(index: GFFIndex,
                    gene_intervals: Dict[str, List[Tuple[int, int]]]):
    """
    Add gene-body intervals at intron priority. Since exon/CDS have higher
    priority in classify(), they always win when a position overlaps both.
    The intron interval only activates in the gaps between exons.
    """
    for contig, genes in gene_intervals.items():
        for g_start, g_end in genes:
            index._add(contig, g_start, g_end, "intron")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_gff_index(path: str, progress_callback=None) -> GFFIndex:
    """
    Parse a GFF3 or GTF file and return a ready-to-query GFFIndex.
    The index also contains chromosome name mappings in index.chrom_names.

    Args:
        path: path to .gff, .gff3, or .gtf file
        progress_callback: optional callable(n_lines_parsed)

    Returns:
        GFFIndex ready for classify() and get_display_name() calls
    """
    index = GFFIndex()
    lower = path.lower()
    if lower.endswith(".gtf"):
        _parse_gtf(path, index, progress_callback)
    else:
        _parse_gff3(path, index, progress_callback)
    index.build()
    return index


def annotate_ssrs(ssrs: list, gff_index: GFFIndex) -> list:
    """
    Add 'genomic_feature' and 'chrom_name' keys to each SSR dict.

    Args:
        ssrs: list of SSR dicts (must have contig, start, end)
        gff_index: built GFFIndex

    Returns:
        Same list with 'genomic_feature' and 'chrom_name' added.
    """
    for ssr in ssrs:
        ssr["genomic_feature"] = gff_index.classify(
            ssr["contig"], ssr["start"], ssr["end"]
        )
        ssr["chrom_name"] = gff_index.get_display_name(ssr["contig"])
    return ssrs