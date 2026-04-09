"""
specificity_panel.py — Specificity Check
Single scroll area. Runs BLAST. Stores results. Shows pass count only.
"""
import os, sys, copy, shutil

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox, QComboBox,
    QProgressBar, QScrollArea, QLineEdit,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_SIZE_LARGE, FONT_SIZE_SMALL, PANEL_PADDING,
)


def _lbl(text, tip=None):
    w = QLabel(text)
    w.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
    if tip: w.setToolTip(tip)
    return w


class BLASTWorker(QThread):
    finished = pyqtSignal(list)
    error    = pyqtSignal(str)

    def __init__(self, genome_path, primers, blast_params, specificity_params, blast_bin_dir):
        super().__init__()
        self.genome_path        = genome_path
        self.primers            = primers
        self.blast_params       = blast_params
        self.specificity_params = specificity_params
        self.blast_bin_dir      = blast_bin_dir

    def run(self):
        try:
            _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            if _root not in sys.path:
                sys.path.insert(0, _root)
            from core.primer_specificity_blast import check_specificity_blast
            results = check_specificity_blast(
                genome_fasta=self.genome_path,
                primer_results=self.primers,
                blast_params=self.blast_params,
                specificity_params=self.specificity_params,
                blast_bin_dir=self.blast_bin_dir,
            )
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))



def _flatten_blast_hits(spec_results):
    """Flatten nested BLAST hit lists into a flat list of dicts."""
    import pandas as pd
    rows = []
    for rec in spec_results:
        ssr_id = rec.get("ssr_id")
        for side, hits in [("Forward", rec.get("left_hits", [])),
                           ("Reverse", rec.get("right_hits", []))]:
            for h in (hits or []):
                rows.append({
                    "ssr_id": ssr_id, "primer_side": side,
                    "qseqid": h.get("qseqid"), "sseqid": h.get("sseqid"),
                    "pident": h.get("pident"), "length": h.get("length"),
                    "mismatch": h.get("mismatch"), "gapopen": h.get("gapopen"),
                    "qstart": h.get("qstart"), "qend": h.get("qend"),
                    "sstart": h.get("sstart"), "send": h.get("send"),
                    "evalue": h.get("evalue"), "bitscore": h.get("bitscore"),
                    "sstrand": h.get("sstrand"),
                })
    cols = ["ssr_id","primer_side","qseqid","sseqid","pident","length",
            "mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sstrand"]
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)

class SpecificityPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker = None
        self._build_ui()

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        outer.addWidget(scroll)

        content = QWidget()
        scroll.setWidget(content)
        L = QVBoxLayout(content)
        L.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        L.setSpacing(16)

        # Title
        title = QLabel("Specificity Check")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("BLAST primers against the genome to verify unique binding")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # Primer queue
        queue_group = QGroupBox("Primer Queue")
        qg = QVBoxLayout(queue_group)
        self.queue_label = QLabel("No primers available")
        self.queue_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        qg.addWidget(self.queue_label)
        L.addWidget(queue_group)

        # BLAST path
        blast_group = QGroupBox("BLAST+ Location")
        bg = QVBoxLayout(blast_group)
        self.blast_detected_label = QLabel("")
        self.blast_detected_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        bg.addWidget(self.blast_detected_label)
        path_row = QHBoxLayout()
        path_row.addWidget(_lbl("BLAST+ bin directory"))
        self.blast_path_input = QLineEdit()
        self.blast_path_input.setPlaceholderText(r"C:\Program Files\NCBI\blast-2.16.0+\bin  (leave blank if on PATH)")
        self.blast_path_input.setToolTip("Leave blank if blastn is on your system PATH.\nOn Windows, paste the path to your BLAST+ bin folder.")
        self.blast_path_input.textChanged.connect(self._on_path_changed)
        path_row.addWidget(self.blast_path_input)
        bg.addLayout(path_row)
        L.addWidget(blast_group)

        # Settings
        settings_group = QGroupBox("Specificity Settings")
        sg = QGridLayout(settings_group)
        sg.setSpacing(10); sg.setContentsMargins(12, 12, 12, 12)

        sg.addWidget(_lbl("Pass/fail criterion",
            "Hits: each primer must map to exactly 1 genomic location.\n"
            "Amplicons: exactly 1 in-range PCR product must be formed."), 0, 0)
        self.pass_mode = QComboBox()
        self.pass_mode.addItem("Hits — each primer maps to exactly 1 genomic location", "hits")
        self.pass_mode.addItem("Amplicons — exactly 1 in-range PCR product formed", "amplicons")
        sg.addWidget(self.pass_mode, 0, 1, 1, 3)

        sg.addWidget(_lbl("Max locations per primer",
            "Maximum number of genomic locations a primer may map to and still pass (hits mode)."), 1, 0)
        self.max_hits = QSpinBox(); self.max_hits.setRange(1, 10); self.max_hits.setValue(1)
        sg.addWidget(self.max_hits, 1, 1)

        sg.addWidget(_lbl("Max mismatches",
            "Maximum mismatches allowed in a BLAST alignment to count as a hit."), 1, 2)
        self.max_mismatches = QSpinBox(); self.max_mismatches.setRange(0, 10); self.max_mismatches.setValue(0)
        sg.addWidget(self.max_mismatches, 1, 3)

        sg.addWidget(_lbl("Min alignment length (bp)",
            "Minimum alignment length to count as a significant hit."), 2, 0)
        self.min_align_len = QSpinBox(); self.min_align_len.setRange(5, 30); self.min_align_len.setValue(18)
        sg.addWidget(self.min_align_len, 2, 1)

        sg.addWidget(_lbl("Min % identity",
            "Minimum percentage identity for a hit to be counted."), 2, 2)
        self.min_identity = QDoubleSpinBox(); self.min_identity.setRange(50, 100); self.min_identity.setValue(100); self.min_identity.setSingleStep(1)
        sg.addWidget(self.min_identity, 2, 3)

        sg.addWidget(_lbl("BLAST task",
            "blastn-short: optimised for short sequences like primers (recommended).\nblastn: standard nucleotide BLAST."), 3, 0)
        self.blast_task = QComboBox()
        self.blast_task.addItem("blastn-short (recommended for primers)", "blastn-short")
        self.blast_task.addItem("blastn", "blastn")
        sg.addWidget(self.blast_task, 3, 1, 1, 3)

        sg.addWidget(_lbl("Word size",
            "BLAST word size. 7 is recommended for short primer sequences."), 4, 0)
        self.word_size = QSpinBox(); self.word_size.setRange(4, 28); self.word_size.setValue(7)
        sg.addWidget(self.word_size, 4, 1)

        sg.addWidget(_lbl("E-value",
            "Maximum E-value for a BLAST hit to be reported.\n0.001 is standard for primer specificity checking."), 4, 2)
        self.evalue = QDoubleSpinBox(); self.evalue.setRange(0.0001, 1.0); self.evalue.setValue(0.001); self.evalue.setSingleStep(0.001); self.evalue.setDecimals(4)
        sg.addWidget(self.evalue, 4, 3)

        sg.addWidget(_lbl("Dust filter",
            "DUST masks low-complexity regions before BLAST alignment.\n"
            "Recommended: OFF for primer specificity — SSR-flanking primers often contain\n"
            "repetitive sequence which DUST will mask, causing them to return no hits\n"
            "and be classified as UNKNOWN rather than PASS or FAIL."), 5, 0)
        self.dust_filter = QComboBox()
        self.dust_filter.addItem("Off (recommended for primers)", "no")
        self.dust_filter.addItem("On", "yes")
        sg.addWidget(self.dust_filter, 5, 1)

        sg.addWidget(_lbl("Max off-target amplicon (bp)",
            "In amplicons mode, the maximum size of an off-target product to count as a failure."), 6, 0)
        self.max_offtarget = QSpinBox(); self.max_offtarget.setRange(100, 5000); self.max_offtarget.setValue(3000); self.max_offtarget.setSingleStep(50)
        sg.addWidget(self.max_offtarget, 6, 1)

        sg.setColumnStretch(4, 1)
        L.addWidget(settings_group)

        # Run
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Run Specificity Check")
        self.run_btn.setObjectName("primary")
        self.run_btn.clicked.connect(self._run)
        run_row.addWidget(self.run_btn); run_row.addStretch()
        L.addLayout(run_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setFixedHeight(6)
        self.progress_bar.setRange(0, 0)  # indeterminate
        L.addWidget(self.progress_bar)

        self.status_label = QLabel("")
        self.status_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(self.status_label)

        L.addStretch()
        self._refresh()

    # ---------------------------------------------------------
    # PATH
    # ---------------------------------------------------------
    def _on_path_changed(self, text):
        if text.strip():
            self.state.blast_bin_dir = text.strip()

    # ---------------------------------------------------------
    # RUN
    # ---------------------------------------------------------
    def _run(self):
        if not self.state.has_genome:
            self._set_status("Load a genome first", WARNING); return
        if not self.state.has_primers:
            self._set_status("Run primer design first", WARNING); return
        if not self.state.genome_path:
            self._set_status("Genome file path not available — reload genome on Home page", WARNING); return

        blast_bin_dir = self.state.blast_bin_dir or None
        if not blast_bin_dir and not shutil.which("blastn"):
            self._set_status("BLAST+ not found. Install BLAST+ and enter the bin directory above.", ERROR)
            return

        # Snapshot primers — local variable, immune to any state changes during run
        _primers = list(self.state.blast_primers)
        if not _primers:
            self._set_status("No primers to BLAST", WARNING); return

        _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
        if _root not in sys.path: sys.path.insert(0, _root)
        from core.primer_specificity_blast import PRESET_SPECIFICITY, PRESET_BLAST

        blast_params = copy.copy(PRESET_BLAST["Single-locus"])
        blast_params.task      = self.blast_task.currentData()
        blast_params.word_size = self.word_size.value()
        blast_params.evalue    = self.evalue.value()
        blast_params.dust      = self.dust_filter.currentData()

        sp = copy.copy(PRESET_SPECIFICITY["Single-locus"])
        sp.min_product                 = self.state.product_min
        sp.max_product                 = self.state.product_max
        sp.max_mismatches              = self.max_mismatches.value()
        sp.min_align_len               = self.min_align_len.value()
        sp.min_identity                = self.min_identity.value()
        sp.max_offtarget_amplicon_size = self.max_offtarget.value()
        sp.pass_mode                   = self.pass_mode.currentData()
        sp.max_total_hits              = self.max_hits.value()

        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self._set_status(f"Running BLAST on {len(_primers):,} primers...", TEXT_SECONDARY)
        self.mw.set_status("Running BLAST...")

        self._worker = BLASTWorker(
            self.state.genome_path, _primers, blast_params, sp, blast_bin_dir
        )
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, results):
        self.state.specificity_results = results
        self.state.blast_raw_rows      = _flatten_blast_hits(results)
        self.state.blast_version      += 1

        n_pass = sum(1 for r in results if r.get("specificity_status") == "PASS")
        n_fail = len(results) - n_pass

        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self._set_status(
            f"BLAST complete — {n_pass:,} PASS, {n_fail:,} FAIL. Navigate to BLAST Results to view.",
            SUCCESS
        )
        self.mw.set_status(f"Specificity check complete — {n_pass:,} primers passed")
        self.mw.on_step_complete(3)

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        if "not found" in msg.lower() or "no such file" in msg.lower():
            self._set_status("BLAST+ executables not found. Check the bin directory path above.", ERROR)
        else:
            self._set_status(f"BLAST failed: {msg}", ERROR)
        self.mw.set_status("Specificity check failed")

    def _set_status(self, msg, color):
        self.status_label.setText(msg)
        self.status_label.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        # Restore persisted BLAST path
        if self.state.blast_bin_dir:
            self.blast_path_input.setText(self.state.blast_bin_dir)

        # Auto-detect BLAST
        blastn = shutil.which("blastn")
        if blastn:
            self.blast_detected_label.setText(f"BLAST+ detected on PATH — {blastn}")
            self.blast_detected_label.setStyleSheet(f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.blast_detected_label.setText("BLAST+ not detected on PATH — enter bin directory below or install from ncbi.nlm.nih.gov/blast")
            self.blast_detected_label.setStyleSheet(f"color: {WARNING}; font-size: {FONT_SIZE_SMALL}pt;")

        # Queue info
        if self.state.has_primers:
            _filtered = self.state.filtered_primer_results
            _all      = self.state.primer_results or []
            n_blast   = len(_filtered) if _filtered else len(_all)
            n_total   = len(_all)
            if n_blast == n_total:
                self.queue_label.setText(f"{n_blast:,} primer pairs queued for BLAST — no quality filters applied")
            else:
                self.queue_label.setText(
                    f"{n_blast:,} primer pairs queued "
                    f"({n_total - n_blast:,} excluded by quality filters on Primer Design page)"
                )
            self.queue_label.setStyleSheet(f"color: {TEXT_PRIMARY}; font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.queue_label.setText("No primers available — run Primer Design first")
            self.queue_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")

        # Show last run result if available
        if self.state.has_specificity and not self.status_label.text():
            n_pass = sum(1 for r in self.state.specificity_results if r.get("specificity_status") == "PASS")
            self._set_status(
                f"Last run: {n_pass:,} primers passed — navigate to BLAST Results to view, or run again.",
                TEXT_SECONDARY
            )

    def on_show(self):
        self._refresh()