# GRACE
### Genomic Repeat Analysis and Characterisation Engine

GRACE is a standalone, offline, cross-platform desktop application for SSR (microsatellite) marker development. It integrates SSR detection, Primer3-based primer design, BLAST-based specificity validation, and lab-ready CSV export in a single graphical interface — no command line required.

---

## Features

- **SSR Detection** — compiled regex-based scanning with multiprocessing, handles highly fragmented assemblies efficiently
- **Primer Design** — full Primer3 integration with user-configurable parameters
- **Quality Filtering** — 3′ end stability (SantaLucia 1998) and GC clamp filters before BLAST
- **Specificity Checking** — BLAST-based in silico PCR against your local genome
- **Offline** — works entirely offline; supports unpublished and proprietary genome assemblies
- **Cross-platform** — Windows, macOS, Linux

---

## Installation

### Pre-built binaries (recommended)

Download the latest release for your platform from the [Releases](../../releases) page. No installation required — extract the folder and run `GRACE.exe` (Windows) or `GRACE` (Mac/Linux).

### Requirements

- **BLAST+** must be installed separately and either placed on your system PATH or the bin directory specified within the application.
  - Download: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

### Running from source

```bash
git clone https://github.com/jackburton-qub/GRACE.git
cd GRACE
pip install -r requirements.txt
python main.py
```

Python 3.12 required.

---

## Building from source

```bash
pip install pyinstaller
python build.py
```

Output will be in `dist/GRACE/`.

---

## Usage

1. **Load Genome** — select a FASTA file (.fa / .fasta / .fna)
2. **Detect SSRs** — configure motif lengths and repeat thresholds, run detection
3. **Design Primers** — set Primer3 parameters, design for all SSRs or a selected subset
4. **Check Specificity** — run BLAST against your genome to validate primer specificity
5. **BLAST Results** — view and download PASS primers as a lab-ready CSV

Session data can be saved and restored via the Home page.

---

## Citation

If you use GRACE in your research, please cite:

> Burton, J.P. & [Supervisor name] (2026) GRACE: an integrated offline desktop application for microsatellite detection, primer design, and BLAST-based specificity validation. *Molecular Ecology Resources*. DOI: [XXXX]

See also [CITATION.cff](CITATION.cff) for machine-readable citation information.

---

## Validation

GRACE was validated on the sugar kelp genome (*Saccharina latissima*; NCBI accession GCA_034768055.1; 251,166 contigs):

- SSR detection concordant with Krait2 at identical 1-based coordinates
- 278 primer pairs passed single-locus BLAST specificity (trinucleotide SSRs)
- 12 primer pairs showed complete sequence identity with Krait2-validated designs

---

## Dependencies

| Package | Purpose |
|---------|---------|
| PyQt6 | Desktop GUI framework |
| primer3-py | Primer design |
| pandas | Data handling |
| numpy | Numerical operations |
| BLAST+ (external) | Specificity checking |

---

## Licence

MIT — see [LICENSE](LICENSE) for details.

---

## Author

Jack Patrick Burton  
School of Biological Sciences, Queen's University Belfast  
[QUB email]
