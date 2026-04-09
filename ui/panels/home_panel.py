"""
home_panel.py — Load Genome + Session Management
"""
import os, sys, json

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QFileDialog, QMessageBox,
)
from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)

SAVEABLE_KEYS = [
    "ssrs", "primer_results", "specificity_results",
    "genome_filename", "product_min", "product_max",
    "blast_bin_dir",
    "filtered_primer_results",
]


def _state_to_dict(state):
    """Serialise saveable state to a plain dict."""
    payload = {}
    for k in SAVEABLE_KEYS:
        v = getattr(state, k, None)
        if v is not None:
            if hasattr(v, "to_dict"):          # DataFrame
                payload[k] = v.to_dict(orient="records")
            else:
                payload[k] = v
    return payload


def _dict_to_state(state, payload):
    """Restore state from a plain dict."""
    import pandas as pd
    for k in SAVEABLE_KEYS:
        if k not in payload:
            continue
        v = payload[k]
        setattr(state, k, v)
    # blast_raw_rows lives on state but isn't in SAVEABLE_KEYS above
    # — it gets rebuilt from specificity_results on next view
    # Version stamps need resetting based on what was restored
    if state.ssrs:
        state.ssr_version = 1
    if state.primer_results:
        state.filtered_primer_results = state.primer_results  # restore filtered to full set
        state.primer_version = 1
    if state.specificity_results:
        state.blast_version = 1





# ── Genome load worker ────────────────────────────────────
class GenomeLoadWorker(QThread):
    finished = pyqtSignal(dict, str)
    error    = pyqtSignal(str)

    def __init__(self, path):
        super().__init__()
        self.path = path

    def run(self):
        try:
            _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            if _root not in sys.path:
                sys.path.insert(0, _root)
            from core.fasta_loader import load_sequence_file
            genome = load_sequence_file(self.path)
            self.finished.emit(genome, self.path)
        except Exception as e:
            self.error.emit(str(e))


# ── Panel ─────────────────────────────────────────────────
class HomePanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker = None
        self._selected_path = None
        self._build_ui()

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        from PyQt6.QtWidgets import QScrollArea
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
        title = QLabel("Load Genome")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Select a FASTA or FASTQ file to begin — .fa / .fasta / .fna / .fastq / .fq supported")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # ── Load genome ───────────────────────────────────
        load_group = QGroupBox("Genome File")
        lg = QVBoxLayout(load_group)
        lg.setSpacing(10)

        path_row = QHBoxLayout()
        self.path_label = QLabel("No file selected")
        self.path_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.path_label.setWordWrap(True)
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.setFixedWidth(100)
        self.browse_btn.clicked.connect(self._browse)
        path_row.addWidget(self.path_label, 1)
        path_row.addWidget(self.browse_btn)
        lg.addLayout(path_row)

        self.load_btn = QPushButton("Load Genome")
        self.load_btn.setObjectName("primary")
        self.load_btn.setEnabled(False)
        self.load_btn.clicked.connect(self._load)
        lg.addWidget(self.load_btn)

        self.load_status = QLabel("")
        self.load_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        lg.addWidget(self.load_status)
        L.addWidget(load_group)

        # ── Genome summary ────────────────────────────────
        self.summary_group = QGroupBox("Genome Summary")
        self.summary_group.setVisible(False)
        self.summary_layout = QGridLayout(self.summary_group)
        self.summary_layout.setSpacing(12)
        L.addWidget(self.summary_group)

        # ── Session management ────────────────────────────
        session_group = QGroupBox("Session Management")
        sg = QVBoxLayout(session_group)

        save_restore_row = QHBoxLayout()

        # Save
        save_col = QVBoxLayout()
        save_lbl = QLabel("Save session")
        save_lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        save_lbl.setStyleSheet(f"color: {TEXT_PRIMARY};")
        save_desc = QLabel("Write all results to a JSON file")
        save_desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.save_btn = QPushButton("Save session...")
        self.save_btn.clicked.connect(self._save_session)
        self.save_btn.setToolTip("Saves SSR results, primer results, and BLAST results to a JSON file.\nThe genome file itself is not saved — you will need to reload it.")
        save_col.addWidget(save_lbl)
        save_col.addWidget(save_desc)
        save_col.addWidget(self.save_btn)
        save_restore_row.addLayout(save_col)

        save_restore_row.addSpacing(24)

        # Restore
        restore_col = QVBoxLayout()
        restore_lbl = QLabel("Restore session")
        restore_lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        restore_lbl.setStyleSheet(f"color: {TEXT_PRIMARY};")
        restore_desc = QLabel("Load results from a previously saved JSON file")
        restore_desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.restore_btn = QPushButton("Restore session...")
        self.restore_btn.clicked.connect(self._restore_session)
        self.restore_btn.setToolTip("Load a previously saved session.\nLoad your genome first, then restore to continue where you left off.")
        restore_col.addWidget(restore_lbl)
        restore_col.addWidget(restore_desc)
        restore_col.addWidget(self.restore_btn)
        save_restore_row.addLayout(restore_col)

        save_restore_row.addStretch()
        sg.addLayout(save_restore_row)

        self.session_status = QLabel("")
        self.session_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.session_status.setWordWrap(True)
        sg.addWidget(self.session_status)
        L.addWidget(session_group)

        # ── Reset ─────────────────────────────────────────
        self.reset_btn = QPushButton("Reset — clear all data and start fresh")
        self.reset_btn.setStyleSheet(f"color: {ERROR}; border-color: {ERROR};")
        self.reset_btn.clicked.connect(self._reset)
        self.reset_btn.setVisible(False)
        L.addWidget(self.reset_btn)

        L.addStretch()
        self._refresh()

    # ---------------------------------------------------------
    # GENOME LOAD
    # ---------------------------------------------------------
    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select sequence file", "",
            "Sequence files (*.fa *.fasta *.fna *.fastq *.fq);;FASTA files (*.fa *.fasta *.fna);;FASTQ files (*.fastq *.fq);;All files (*)"
        )
        if path:
            self._selected_path = path
            self.path_label.setText(path)
            self.path_label.setStyleSheet(f"color: {TEXT_PRIMARY}; font-size: {FONT_SIZE_SMALL}pt;")
            self.load_btn.setEnabled(True)

    def _load(self):
        if not self._selected_path:
            return
        self.load_btn.setEnabled(False)
        self.browse_btn.setEnabled(False)
        self.load_status.setText("Loading...")
        self.load_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.mw.set_status("Loading genome...")
        self._worker = GenomeLoadWorker(self._selected_path)
        self._worker.finished.connect(self._on_load_done)
        self._worker.error.connect(self._on_load_error)
        self._worker.start()

    def _on_load_done(self, genome, path):
        self.state.clear_downstream_of_genome()
        self.state.genome          = genome
        self.state.genome_path     = path
        self.state.genome_filename = os.path.basename(path)
        self.load_status.setText(f"Loaded — {len(genome):,} contigs")
        self.load_status.setStyleSheet(f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;")
        self.browse_btn.setEnabled(True)
        self.mw.set_status(f"Genome loaded — {len(genome):,} contigs")
        self.mw.on_step_complete(0)
        self._refresh()

    def _on_load_error(self, msg):
        self.load_status.setText(f"Error: {msg}")
        self.load_status.setStyleSheet(f"color: {ERROR}; font-size: {FONT_SIZE_SMALL}pt;")
        self.load_btn.setEnabled(True)
        self.browse_btn.setEnabled(True)
        self.mw.set_status("Failed to load genome")

    # ---------------------------------------------------------
    # SESSION SAVE / RESTORE
    # ---------------------------------------------------------
    def _save_session(self):
        default = os.path.join(
            os.path.dirname(self.state.genome_path) if self.state.genome_path else os.path.expanduser("~"),
            "grace_session.json"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save session", default, "JSON files (*.json)"
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(_state_to_dict(self.state), f, default=str, indent=2)
            self._set_session_status(f"Session saved to {path}", SUCCESS)
            self.mw.set_status(f"Session saved")
        except Exception as e:
            self._set_session_status(f"Failed to save: {e}", ERROR)

    def _restore_session(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Restore session", os.path.expanduser("~"), "JSON files (*.json)"
        )
        if not path:
            return
        self._load_session_file(path)

    def _load_session_file(self, path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            _dict_to_state(self.state, payload)
            parts = []
            if self.state.ssrs:             parts.append(f"{len(self.state.ssrs):,} SSRs")
            if self.state.primer_results:   parts.append(f"{len(self.state.primer_results):,} primer pairs")
            if self.state.specificity_results: parts.append(f"{len(self.state.specificity_results):,} specificity results")
            summary = ", ".join(parts) if parts else "no result data"
            self._set_session_status(f"Session restored — {summary}", SUCCESS)
            # Update sidebar step indicators
            if self.state.has_genome:    self.mw.on_step_complete(0)
            if self.state.has_ssrs:      self.mw.on_step_complete(1)
            if self.state.has_primers:   self.mw.on_step_complete(2)
            if self.state.has_specificity: self.mw.on_step_complete(3)
            self.mw.refresh_sidebar()
            self.mw.set_status("Session restored")
            self._refresh()
        except Exception as e:
            self._set_session_status(f"Failed to restore: {e}", ERROR)

    def _set_session_status(self, msg, color):
        self.session_status.setText(msg)
        self.session_status.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # RESET
    # ---------------------------------------------------------
    def _reset(self):
        reply = QMessageBox.question(
            self, "Reset", "Clear all data and start fresh?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.state.clear_downstream_of_genome()
            self.state.genome = self.state.genome_path = self.state.genome_filename = None
            self._selected_path = None
            self.path_label.setText("No file selected")
            self.path_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
            self.load_btn.setEnabled(False)
            self.load_status.setText("")
            self.mw.set_status("Reset")
            self.mw.refresh_sidebar()
            self._refresh()

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        has = self.state.has_genome
        self.summary_group.setVisible(has)
        self.reset_btn.setVisible(has)



        if has:
            genome   = self.state.genome
            total_bp = sum(len(s) for s in genome.values())
            while self.summary_layout.count():
                item = self.summary_layout.takeAt(0)
                if item.widget(): item.widget().deleteLater()
            metrics = [
                ("File",              self.state.genome_filename or "—"),
                ("Contigs",           f"{len(genome):,}"),
                ("Total length",      f"{total_bp:,} bp"),
                ("Shortest contig",   f"{min(len(s) for s in genome.values()):,} bp"),
                ("Longest contig",    f"{max(len(s) for s in genome.values()):,} bp"),
            ]
            for row, (label, value) in enumerate(metrics):
                lbl = QLabel(label)
                lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
                val = QLabel(value)
                val.setFont(QFont(FONT_MONO, FONT_SIZE_NORMAL))
                val.setStyleSheet(f"color: {TEXT_PRIMARY};")
                self.summary_layout.addWidget(lbl, row, 0)
                self.summary_layout.addWidget(val, row, 1)

    def on_show(self):
        self._refresh()