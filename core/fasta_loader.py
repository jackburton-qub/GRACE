"""
FASTA / FASTQ Loader Module (backend-safe)
------------------------------------------
Provides clean functions for loading FASTA and FASTQ files into a dictionary.

Supports:
- Local file paths (string)          — reads line by line, no full-file buffer
- File-like objects (BytesIO, StringIO)
- Streamlit uploaded files

Performance notes
-----------------
- File paths are read line-by-line to avoid loading the full file into memory
- The IUPAC filter uses str.translate() with a precomputed delete-table
  which runs in C and is ~10x faster than a character-by-character generator
- Sequence chunks are joined once per contig rather than concatenated
"""

import io
import string


# ------------------------------------------------------------------
# Precompute a translate table that strips any character NOT in the
# IUPAC nucleotide alphabet.  str.translate() runs in C.
# ------------------------------------------------------------------
_IUPAC = set("ACGTNRYSWKMBDHV-")
_DELETE_TABLE = str.maketrans("", "", "".join(
    c for c in (string.ascii_letters + string.digits + string.punctuation + " ")
    if c.upper() not in _IUPAC
))


def _clean(line: str) -> str:
    """Strip non-IUPAC characters from a sequence line using C-level translate."""
    return line.upper().translate(_DELETE_TABLE)


def load_fasta(handle) -> dict:
    """
    Load a FASTA file from a file path or file-like object.
    Returns a dict: {sequence_id: sequence}.

    Handles:
    - str path      → opens and reads line-by-line (memory efficient for large files)
    - bytes file-like (e.g. Streamlit upload, BytesIO) → decodes to text, iterates lines
    - str file-like  (e.g. StringIO) → iterates lines directly
    """

    # ------------------------------------------------------------------
    # Normalise the input into an iterable of str lines
    # ------------------------------------------------------------------
    if isinstance(handle, str):
        # File path — open directly and stream line by line
        file_obj = open(handle, "r", encoding="utf-8", errors="replace")
        close_after = True
    else:
        content = handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        file_obj = io.StringIO(content)
        close_after = True

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------
    sequences  = {}
    current_id = None
    chunks     = []          # accumulate sequence chunks per contig

    try:
        for raw_line in file_obj:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous contig
                if current_id is not None:
                    sequences[current_id] = "".join(chunks)
                # Start new contig — first token after ">"
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(_clean(line))

        # Save last contig
        if current_id is not None:
            sequences[current_id] = "".join(chunks)

    finally:
        if close_after:
            file_obj.close()

    return sequences


def load_fastq(handle) -> dict:
    """
    Load a FASTQ file from a file path or file-like object.
    Returns a dict: {sequence_id: sequence}.

    FASTQ format is 4 lines per record:
      Line 1: @sequence_id [description]
      Line 2: sequence
      Line 3: + (separator, ignored)
      Line 4: quality scores (ignored)
    """

    # ------------------------------------------------------------------
    # Normalise the input into an iterable of str lines (same as load_fasta)
    # ------------------------------------------------------------------
    if isinstance(handle, str):
        file_obj = open(handle, "r", encoding="utf-8", errors="replace")
        close_after = True
    else:
        content = handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        file_obj = io.StringIO(content)
        close_after = True

    sequences = {}

    try:
        lines = iter(file_obj)
        for header_line in lines:
            header_line = header_line.strip()
            if not header_line:
                continue
            if not header_line.startswith("@"):
                continue  # skip malformed / unexpected lines

            seq_id   = header_line[1:].split()[0]
            seq_line = next(lines, "").strip()
            next(lines, None)   # "+" separator — discard
            next(lines, None)   # quality scores  — discard

            if seq_id and seq_line:
                sequences[seq_id] = _clean(seq_line)

    finally:
        if close_after:
            file_obj.close()

    return sequences


def load_sequence_file(path_or_handle) -> dict:
    """
    Auto-detect FASTA vs FASTQ by file extension (for paths) or first
    non-empty character (for file-like objects), then dispatch accordingly.
    Returns a dict: {sequence_id: sequence}.
    """
    if isinstance(path_or_handle, str):
        lower = path_or_handle.lower()
        if lower.endswith((".fastq", ".fq")):
            return load_fastq(path_or_handle)
        return load_fasta(path_or_handle)

    # File-like object — peek at the first character
    content = path_or_handle.read()
    if isinstance(content, bytes):
        content = content.decode("utf-8", errors="replace")

    first_char = content.lstrip()[:1]
    buf = io.StringIO(content)
    if first_char == "@":
        return load_fastq(buf)
    return load_fasta(buf)