# GRACE
### Genomic Repeat Analysis and Characterisation Engine

GRACE is a standalone, offline, cross‑platform desktop application for microsatellite (SSR) marker development. It combines SSR detection, Primer3‑based primer design, BLAST specificity filtering, genomic annotation, and lab‑ready exports in a single graphical interface — no command line needed.

---

## Features

- **SSR Detection** – fast regex‑based scanning for di‑ to hexanucleotide repeats, with canonical motif standardisation and optional adjacent‑motif filtering.
- **SSR Summary** – interactive karyotype view showing SSR positions along chromosomes, motif type distribution, and genomic feature breakdown (exon/intron/intergenic) when a GFF is loaded.
- **Primer Design** – Primer3‑powered design with quality filters. Two parameter presets:
  - **Capillary Electrophoresis** – optimised for fragment analysis (M13 tails, dye multiplexing).
  - **Amplicon Sequencing** – optimised for Illumina amplicon panels (P5/P7 tails, tight Tm, short products).
- **Specificity Checking** – bundled BLAST+ binaries; no external installation required.
- **GFF Annotation** – optional GFF3/GTF support for exon/intron/intergenic classification and chromosome name mapping.
- **Three Complete Workflows**
  - **Amplicon Sequencing** – validate multiplex primer pools for dimer risks, amplicon size, and GC content. Optimise panels and export tailed primers with amplicon reference sequences.
  - **Capillary Electrophoresis** – assign dyes (FAM, VIC, NED, PET) and size bins. Supports multiple panels (injections) and exports ready‑to‑order primers with M13 tails.
  - **GBS‑RE (Restriction‑Enzyme GBS)** – discover SSRs located within user‑defined restriction fragments (e.g. PstI‑MspI). Compare enzyme combinations, generate barcodes, and export marker lists.
- **Linkage Disequilibrium Filter** – thin markers by physical distance across all workflows.
- **Chromosome View** – visualise marker positions on chromosomes or contigs throughout the application.
- **Project Summary** – dynamic, tabbed summary of all completed steps with embedded charts. Export a comprehensive PDF report.

---

## Installation

### Pre‑built Binaries
Download the latest release for your platform from the [Releases](https://github.com/jackburton-qub/GRACE/releases) page. Extract the folder and run `GRACE.exe` (Windows) or open `GRACE.app` (macOS). BLAST+ is bundled – no separate installation required.

### Running from Source
```bash
git clone https://github.com/jackburton-qub/GRACE.git
cd GRACE
pip install -r requirements.txt
python main.py
Workflow
text
01  Load Genome        – FASTA/FASTQ + optional GFF3/GTF annotation
02  Detect SSRs        – scan for perfect repeats, annotate with genomic features
03  SSR Summary        – visualise distribution and statistics
04  Design Primers     – Capillary or Amplicon Sequencing mode
05  Check Specificity  – BLAST primers against the genome
06  BLAST Results      – review and filter by specificity
07  Project Summary    – complete pipeline overview with PDF export
08  Amplicon Tools     – multiplex validation and panel optimisation
09  Capillary Tools    – dye assignment and size binning
10  GBS‑RE Tools       – restriction fragment marker discovery
Dependencies
Package	Purpose
PyQt6	Desktop GUI framework
primer3‑py	Primer design
pandas	Data handling
numpy	Numerical operations
reportlab	PDF report generation
BLAST+	Bundled – specificity checking
Supported File Formats
Format	Input	Output
FASTA (.fa, .fasta, .fna)	✓	
FASTQ (.fastq, .fq)	✓	
GenBank (.gb, .gbff, .gbk)	✓	
GFF3 (.gff, .gff3)	✓	
GTF (.gtf)	✓	
CSV		✓
FASTA (primers, amplicons)		✓
PDF		✓
JSON (session)	✓	✓
Licence
MIT – see LICENSE for details.