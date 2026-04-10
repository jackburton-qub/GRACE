# GRACE
### Genomic Repeat Analysis and Characterisation Engine

GRACE is a standalone, offline, cross-platform desktop application for SSR (microsatellite) marker development. It integrates SSR detection, Primer3-based primer design, BLAST-based specificity validation, and lab-ready CSV export in a single graphical interface — no command line required.

---

## Installation

### Pre-built binaries

Download the latest release for your platform from the [Releases](../../releases) page. No installation required — extract the folder and run `GRACE.exe` (Windows) or `GRACE` (Mac/Linux).

### Requirements

- **BLAST+** must be installed separately and either placed on your system PATH or the bin directory specified within the application.
  - Download: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

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

