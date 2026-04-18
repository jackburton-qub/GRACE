# GRACE
### Genomic Repeat Analysis and Characterisation Engine

GRACE is a standalone, offline, cross-platform desktop application for SSR (microsatellite) marker development. It integrates SSR detection, Primer3-based primer design, BLAST-based specificity validation, genomic annotation, and lab-ready export in a single graphical interface — no command line required.

---

## Features

- **SSR Detection** — fast regex-based scanning across all motif lengths (di- to hexanucleotide), with full canonical motif standardisation
- **SSR Summary** — karyotype view showing exact SSR positions along each chromosome, motif type breakdown, and genomic feature distribution
- **Primer Design** — Primer3-based design with low-complexity flank filtering, quality filters, and two design modes:
  - **Standard mode** — optimised for capillary electrophoresis / fragment analysis
  - **GBS mode** — optimised for Illumina amplicon sequencing (genotype-by-sequencing)
- **Specificity Checking** — bundled BLAST+ binaries, no external installation required
- **GFF Annotation** — optional GFF3/GTF annotation support for exon/intron/intergenic classification and chromosome name mapping
- **GBS Tools** — multiplex pool compatibility checker, amplicon sequence exporter, and sequencing platform recommender
- **Final Report** — full pipeline summary with charts, filtering funnel, and PDF export

---

## Installation

### Pre-built binaries

Download the latest release for your platform from the [Releases](../../releases) page. No installation required — extract the folder and run `GRACE.exe` (Windows) or open `GRACE.app` (Mac).

BLAST+ binaries are bundled with the application. No separate BLAST installation is needed.

### Running from source

```bash
git clone https://github.com/jackburton-qub/GRACE.git
cd GRACE
pip install -r requirements.txt
python main.py
```

---

## Workflow

```
01  Load Genome        — FASTA/FASTQ + optional GFF3/GTF annotation
02  Detect SSRs        — scan for perfect repeats, annotate with genomic features
03  SSR Summary        — visualise distribution across chromosomes
04  Design Primers     — Standard or GBS mode
05  Check Specificity  — BLAST primers against the genome
06  BLAST Results      — view and filter results
07  Final Report       — pipeline summary, charts, PDF/CSV/FASTA export
08  GBS Tools          — multiplex checker, amplicon FASTA, platform recommender
```

---

## GBS Mode

GRACE supports genotype-by-sequencing (SSR-GBS) workflows for Illumina amplicon sequencing:

- Short product sizes (80–200bp) compatible with 300bp paired-end reads
- Tighter Tm range for consistent multiplexed PCR
- Illumina M13 adapter tails appended to primer sequences
- Low-complexity flank filtering to remove repeat-adjacent primers
- Multiplex pool compatibility checker for cross-primer dimer detection
- Expected amplicon sequence export for allele-calling tools (STRait Razor, HipSTR)
- Sequencing platform recommender based on allele size range

---

## Dependencies

| Package | Purpose |
|---------|---------|
| PyQt6 | Desktop GUI framework |
| primer3-py | Primer design |
| pandas | Data handling |
| numpy | Numerical operations |
| reportlab | PDF report generation |
| BLAST+ | Bundled — specificity checking |

---

## Supported file formats

| Format | Input | Output |
|--------|-------|--------|
| FASTA (.fa, .fasta, .fna) | ✓ Genome | |
| FASTQ (.fastq, .fq) | ✓ Genome | |
| GFF3 (.gff, .gff3) | ✓ Annotation | |
| GTF (.gtf) | ✓ Annotation | |
| CSV | | ✓ SSRs, primers, results |
| FASTA | | ✓ Primers (BLAST-ready, GBS-tailed, amplicons) |
| PDF | | ✓ Final report |
| JSON | ✓ Session restore | ✓ Session save |

---

## Licence

MIT — see [LICENSE](LICENSE) for details.
