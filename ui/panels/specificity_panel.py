"""
specificity_panel.py — Specificity Check
Single scroll area. Runs BLAST. Stores results. Shows pass count only.
BLAST binaries are now bundled with the app — no path input needed.
"""
import os, sys, copy, shutil

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox, QComboBox,
    QProgressBar, QScrollArea,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_SIZE_LARGE, FONT_SIZE_SMALL, PANEL_PADDING,
    BG_MID, BG_LIGHT, BORDER,
)


def _lbl(text, tip=None):
    w = QLabel(text)
    w.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
    if tip: w.setToolTip(tip)
    return w


class BLASTWorker(QThread):
    finished = pyqtSignal(list)
    error    = pyqtSignal(str)

    def __init__(self, genome_path, primers, blast_params, specificity_params):
        super().__init__()
        self.genome_path        = genome_path
        self.primers            = primers
        self.blast_params       = blast_params
        self.specificity_params = specificity_params

    def run(self):
        try:
            from core.primer_specificity_blast import check_specificity_blast
            results = check_specificity_blast(
                genome_fasta=self.genome_path,
                primer_results=self.primers,
                blast_params=self.blast_params,
                specificity_params=self.specificity_params,
                # blast_bin_dir omitted — resolved automatically from bundled blast/
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

        # BLAST status (bundled — read only, no path input)
        blast_group = QGroupBox("BLAST+ Status")
        bg = QVBoxLayout(blast_group)
        self.blast_detected_label = QLabel("")
        self.blast_detected_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        bg.addWidget(self.blast_detected_label)
        L.addWidget(blast_group)

        # Primer queue
        queue_group = QGroupBox("Primer Queue")
        qg = QVBoxLayout(queue_group)
        self.queue_label = QLabel("No primers available")
        self.queue_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        qg.addWidget(self.queue_label)
        L.addWidget(queue_group)

        # Mode selector
        mode_group = QGroupBox("Specificity Mode")
        mg = QVBoxLayout(mode_group)
        mg.setSpacing(8)

        mode_row = QHBoxLayout()
        mode_row.setSpacing(12)

        self._btn_amplicons = QPushButton("Amplicons")
        self._btn_amplicons.setCheckable(True)
        self._btn_amplicons.setChecked(True)
        self._btn_amplicons.setMinimumHeight(60)
        self._btn_amplicons.setStyleSheet(f"""
            QPushButton {{
                background: {BG_MID};
                border: 1px solid {BORDER};
                border-radius: 8px;
                padding: 12px 20px;
                font-weight: bold;
                color: {TEXT_PRIMARY};
                text-align: left;
            }}
            QPushButton:checked {{
                background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
                    stop:0 {ACCENT}33, stop:1 {ACCENT}11);
                border: 2px solid {ACCENT};
                color: {ACCENT};
            }}
            QPushButton:hover {{
                border: 1px solid {ACCENT}88;
                background: {BG_LIGHT};
            }}
        """)

        self._btn_hits = QPushButton("Hits")
        self._btn_hits.setCheckable(True)
        self._btn_hits.setMinimumHeight(60)
        self._btn_hits.setStyleSheet(self._btn_amplicons.styleSheet())

        self._mode_desc = QLabel(
            "Exactly 1 in-range PCR product must be formed\n"
            "Practical validation · Tolerates weak off-target binding"
        )
        self._mode_desc.setWordWrap(True)
        self._mode_desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt; padding: 4px 0;")

        self._btn_amplicons.toggled.connect(lambda checked: self._on_mode_toggled(checked, "amplicons"))
        self._btn_hits.toggled.connect(lambda checked: self._on_mode_toggled(checked, "hits"))

        mode_row.addWidget(self._btn_amplicons)
        mode_row.addWidget(self._btn_hits)
        mode_row.addStretch()
        mg.addLayout(mode_row)
        mg.addWidget(self._mode_desc)

        L.addWidget(mode_group)

        # Settings
        settings_group = QGroupBox("BLAST Settings")
        sg = QGridLayout(settings_group)
        sg.setSpacing(10); sg.setContentsMargins(12, 12, 12, 12)

        sg.addWidget(_lbl("Max locations per primer",
            "Hits mode: Maximum genomic locations each primer may map to and still pass.\n"
            "Amplicons mode: Not used (primers can bind multiple locations if only 1 amplicon forms)."), 0, 0)
        self.max_hits = QSpinBox(); self.max_hits.setRange(1, 10); self.max_hits.setValue(1)
        self.max_hits.setEnabled(False)  # Disabled by default (Amplicons mode is default)
        sg.addWidget(self.max_hits, 0, 1)

        sg.addWidget(_lbl("Max mismatches",
            "Maximum mismatches allowed in a BLAST alignment to count as a hit."), 0, 2)
        self.max_mismatches = QSpinBox(); self.max_mismatches.setRange(0, 10); self.max_mismatches.setValue(2)
        sg.addWidget(self.max_mismatches, 0, 3)

        sg.addWidget(_lbl("Min alignment length (bp)",
            "Minimum alignment length to count as a significant hit."), 1, 0)
        self.min_align_len = QSpinBox(); self.min_align_len.setRange(5, 30); self.min_align_len.setValue(17)
        sg.addWidget(self.min_align_len, 1, 1)

        sg.addWidget(_lbl("Min % identity",
            "Minimum percentage identity for a hit to be counted."), 1, 2)
        self.min_identity = QDoubleSpinBox(); self.min_identity.setRange(50, 100); self.min_identity.setValue(90.0); self.min_identity.setSingleStep(1)
        sg.addWidget(self.min_identity, 1, 3)

        sg.addWidget(_lbl("BLAST task",
            "blastn-short: optimised for short sequences like primers (recommended).\nblastn: standard nucleotide BLAST."), 2, 0)
        self.blast_task = QComboBox()
        self.blast_task.addItem("blastn-short (recommended for primers)", "blastn-short")
        self.blast_task.addItem("blastn", "blastn")
        sg.addWidget(self.blast_task, 2, 1, 1, 3)

        sg.addWidget(_lbl("Word size",
            "BLAST word size. 7 is recommended for short primer sequences."), 3, 0)
        self.word_size = QSpinBox(); self.word_size.setRange(4, 28); self.word_size.setValue(7)
        sg.addWidget(self.word_size, 3, 1)

        sg.addWidget(_lbl("E-value",
            "Maximum E-value for a BLAST hit to be reported.\n0.001 is now standard for both amplicon and hits mode."), 3, 2)
        self.evalue = QDoubleSpinBox(); self.evalue.setDecimals(4); self.evalue.setRange(0.0001, 1.0); self.evalue.setSingleStep(0.001); self.evalue.setValue(0.001)
        sg.addWidget(self.evalue, 3, 3)

        sg.addWidget(_lbl("Dust filter",
            "DUST masks low-complexity regions before BLAST alignment.\n"
            "Recommended: OFF for primer specificity — SSR-flanking primers often contain\n"
            "repetitive sequence which DUST will mask, causing them to return no hits\n"
            "even when they are actually present in the genome."), 4, 0)
        self.dust_filter = QComboBox()
        self.dust_filter.addItem("OFF (recommended for SSR primers)", "no")
        self.dust_filter.addItem("ON", "yes")
        sg.addWidget(self.dust_filter, 4, 1, 1, 3)

        sg.addWidget(_lbl("Max off-target amplicon size (bp)",
            "If a primer pair produces off-target amplicons, reject if any exceed this size."), 5, 0)
        self.max_offtarget = QSpinBox(); self.max_offtarget.setRange(100, 10000); self.max_offtarget.setSingleStep(100); self.max_offtarget.setValue(1000)
        sg.addWidget(self.max_offtarget, 5, 1, 1, 3)

        L.addWidget(settings_group)

        # Run button
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Run BLAST")
        self.run_btn.setFixedHeight(40)
        self.run_btn.clicked.connect(self._run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setFixedHeight(40)
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        run_row.addWidget(self.run_btn)
        run_row.addWidget(self.cancel_btn)
        run_row.addStretch()
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
    # MODE TOGGLE HANDLER
    # ---------------------------------------------------------
    def _on_mode_toggled(self, checked, mode):
        """Handle mode button toggles and update settings."""
        if not checked:
            return
        
        # Update button states and description
        if mode == "amplicons":
            self._btn_hits.setChecked(False)
            self._mode_desc.setText(
                "Exactly 1 in-range PCR product must be formed\n"
                "Practical validation · Tolerates weak off-target binding"
            )
            # Disable max_hits control - not relevant in amplicons mode
            self.max_hits.setEnabled(False)
        else:  # hits
            self._btn_amplicons.setChecked(False)
            self._mode_desc.setText(
                "Each primer maps to exactly 1 genomic location\n"
                "Strict validation · Requires perfect single-locus specificity"
            )
            # Enable max_hits control - relevant in hits mode
            self.max_hits.setEnabled(True)
        
        # Update settings to match preset
        from core.primer_specificity_blast import PRESET_SPECIFICITY, PRESET_BLAST
        
        preset_name = "Amplicons" if mode == "amplicons" else "Hits"
        spec_preset = PRESET_SPECIFICITY[preset_name]
        blast_preset = PRESET_BLAST[preset_name]
        
        # Update UI to match preset values
        self.max_mismatches.setValue(spec_preset.max_mismatches)
        self.min_align_len.setValue(spec_preset.min_align_len)
        self.min_identity.setValue(spec_preset.min_identity)
        self.max_offtarget.setValue(spec_preset.max_offtarget_amplicon_size)
        self.max_hits.setValue(spec_preset.max_total_hits)
        self.evalue.setValue(blast_preset.evalue)
        self.word_size.setValue(blast_preset.word_size)

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

        # Snapshot primers — local variable, immune to any state changes during run
        _primers = list(self.state.blast_primers)
        if not _primers:
            self._set_status("No primers to BLAST", WARNING); return

        from core.primer_specificity_blast import PRESET_SPECIFICITY, PRESET_BLAST

        # Determine which preset to use based on current mode (from button states)
        current_mode = "amplicons" if self._btn_amplicons.isChecked() else "hits"
        preset_name = "Amplicons" if current_mode == "amplicons" else "Hits"

        blast_params = copy.copy(PRESET_BLAST[preset_name])
        blast_params.task      = self.blast_task.currentData()
        blast_params.word_size = self.word_size.value()
        blast_params.evalue    = self.evalue.value()
        blast_params.dust      = self.dust_filter.currentData()

        sp = copy.copy(PRESET_SPECIFICITY[preset_name])
        sp.min_product                 = self.state.product_min
        sp.max_product                 = self.state.product_max
        sp.max_mismatches              = self.max_mismatches.value()
        sp.min_align_len               = self.min_align_len.value()
        sp.min_identity                = self.min_identity.value()
        sp.max_offtarget_amplicon_size = self.max_offtarget.value()
        sp.pass_mode                   = current_mode
        sp.max_total_hits              = self.max_hits.value()

        self.run_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True)
        self._set_status(f"Running BLAST on {len(_primers):,} primers...", TEXT_SECONDARY)
        self.mw.set_status("Running BLAST...")

        self._worker = BLASTWorker(
            self.state.genome_path, _primers, blast_params, sp
        )
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, results):
        print(f"DEBUG BLAST: got {len(results)} results")
        self.state.specificity_results = results
        print(f"DEBUG BLAST: state.specificity_results len = {len(self.state.specificity_results)}")
        self.state.blast_raw_rows = _flatten_blast_hits(results)
        self.state.blast_version += 1
        print(f"DEBUG BLAST: calling mw.refresh_sidebar")
        self.mw.refresh_sidebar()

        n_pass = sum(1 for r in results if r.get("specificity_status") == "PASS")
        n_fail = len(results) - n_pass

        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status(
            f"BLAST complete — {n_pass:,} PASS, {n_fail:,} FAIL. Navigate to BLAST Results to view.",
            SUCCESS
        )
        self.mw.set_status(f"Specificity check complete — {n_pass:,} primers passed")
        self.mw.on_step_complete(4)
        self._cleanup_worker()

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        if "not found" in msg.lower() or "no such file" in msg.lower():
            self._set_status("BLAST+ executables not found — please report this issue.", ERROR)
        else:
            self._set_status(f"BLAST failed: {msg}", ERROR)
        self.mw.set_status("Specificity check failed")
        self._cleanup_worker()

    def _cancel(self):
        if self._worker is not None:
            self._worker.terminate()
            self._worker.wait()
            self._cleanup_worker()
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status("Cancelled", WARNING)
        self.mw.set_status("Specificity check cancelled")

    def _cleanup_worker(self):
        """Disconnect signals and null large refs to prevent idle crashes."""
        if self._worker is not None:
            try:
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            self._worker.primers            = None
            self._worker.blast_params       = None
            self._worker.specificity_params = None
            self._worker = None

    def _set_status(self, msg, color):
        self.status_label.setText(msg)
        self.status_label.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        # Show bundled BLAST status
        from core.primer_specificity_blast import get_blast_bin_dir
        blast_dir = get_blast_bin_dir()
        exe_suffix = ".exe" if sys.platform == "win32" else ""
        blastn_path = os.path.join(blast_dir, f"blastn{exe_suffix}") if blast_dir else ""

        if blast_dir and os.path.isfile(blastn_path):
            self.blast_detected_label.setText(f"✓ Bundled BLAST+ ready — {blastn_path}")
            self.blast_detected_label.setStyleSheet(f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;")
        elif shutil.which("blastn"):
            found = shutil.which("blastn")
            self.blast_detected_label.setText(f"✓ BLAST+ detected on system PATH — {found}")
            self.blast_detected_label.setStyleSheet(f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.blast_detected_label.setText("✗ BLAST+ not found — please report this issue")
            self.blast_detected_label.setStyleSheet(f"color: {ERROR}; font-size: {FONT_SIZE_SMALL}pt;")

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