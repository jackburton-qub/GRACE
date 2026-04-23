"""
primer_panel.py — Primer Design
Single scroll area. Clean layout. Filters appear after design.
Supports Capillary Electrophoresis and Amplicon Sequencing modes.
"""
import os, sys, time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog, QTabWidget,
    QCheckBox, QFrame, QApplication,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)

TABLE_DISPLAY_LIMIT          = 10_000
PRIMER_DESIGN_WARN_THRESHOLD = 50_000


def _lbl(text, tip=None):
    w = QLabel(text)
    w.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
    if tip: w.setToolTip(tip)
    return w


class PrimerWorker(QThread):
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(dict, float)
    error    = pyqtSignal(str)

    def __init__(self, genome, ssr_list, params):
        super().__init__()
        self.genome   = genome
        self.ssr_list = ssr_list
        self.params   = params

    def run(self):
        try:
            from core.primer_design import design_primers_for_all_ssrs
            start = time.time()
            res = design_primers_for_all_ssrs(
                genome=self.genome,
                ssr_list=self.ssr_list,
                flank=self.params["flank"],
                product_size_range=(self.params["product_min"], self.params["product_max"]),
                preset=self.params["preset"],
                primer_opts=self.params["primer_opts"],
                num_pairs=self.params["num_pairs"],
                progress_callback=lambda d, t: self.progress.emit(d, t),
                amplicon_mode=self.params["amplicon_mode"],
            )
            self.finished.emit(res, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


class PrimerPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._worker = None
        self._last_rendered_version = -1
        self._amplicon_mode = False
        self._filter_updating = False
        self._build_ui()

    def _on_mode_toggled(self, checked, mode):
        if not checked:
            return
        if mode == "capillary":
            self._btn_amplicon.setChecked(False)
            self._amplicon_mode = False
            self.state.workflow_mode = "capillary"
            self._mode_desc.setText(
                "Fragment analysis · M13 tails · Dye multiplexing\n"
                "Product sizes 100–350bp · Flank 100bp"
            )
            self.flank.setValue(100)
            self.product_min.setValue(100)
            self.product_max.setValue(250)
            self.tm_min.setValue(52.0)
            self.tm_opt.setValue(58.0)
            self.tm_max.setValue(60.0)
            self.max_tm_diff.setValue(2.0)
            self.max_poly_x.setValue(4)
        else:  # amplicon
            self._btn_capillary.setChecked(False)
            self._amplicon_mode = True
            self.state.workflow_mode = "amplicon"
            self._mode_desc.setText(
                "Illumina amplicon sequencing · P5/P7 tails\n"
                "Product sizes 80–200bp · Flank 50bp · Tighter Tm range"
            )
            self.flank.setValue(50)
            self.product_min.setValue(80)
            self.product_max.setValue(200)
            self.tm_min.setValue(58.0)
            self.tm_opt.setValue(60.0)
            self.tm_max.setValue(62.0)
            self.max_tm_diff.setValue(1.0)
            self.max_poly_x.setValue(3)
        QApplication.processEvents()

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

        title = QLabel("Primer Design")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Design PCR primers flanking each SSR using Primer3")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # Mode selector
        mode_group = QGroupBox("Design Mode")
        mg = QVBoxLayout(mode_group)
        mg.setSpacing(8)

        mode_row = QHBoxLayout()
        mode_row.setSpacing(12)

        self._btn_capillary = QPushButton("Capillary Electrophoresis")
        self._btn_capillary.setCheckable(True)
        self._btn_capillary.setChecked(True)
        self._btn_capillary.setMinimumHeight(60)
        self._btn_capillary.setStyleSheet(f"""
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

        self._btn_amplicon = QPushButton("Amplicon Sequencing")
        self._btn_amplicon.setCheckable(True)
        self._btn_amplicon.setMinimumHeight(60)
        self._btn_amplicon.setStyleSheet(self._btn_capillary.styleSheet())

        self._mode_desc = QLabel(
            "Fragment analysis · M13 tails · Dye multiplexing\n"
            "Product sizes 100–350bp · Flank 100bp"
        )
        self._mode_desc.setWordWrap(True)
        self._mode_desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt; padding: 4px 0;")

        self._btn_capillary.toggled.connect(lambda checked: self._on_mode_toggled(checked, "capillary"))
        self._btn_amplicon.toggled.connect(lambda checked: self._on_mode_toggled(checked, "amplicon"))

        mode_row.addWidget(self._btn_capillary)
        mode_row.addWidget(self._btn_amplicon)
        mode_row.addStretch()
        mg.addLayout(mode_row)
        mg.addWidget(self._mode_desc)

        L.addWidget(mode_group)

        # Settings tabs
        tabs = QTabWidget()

        basic = QWidget()
        bg = QGridLayout(basic)
        bg.setSpacing(10); bg.setContentsMargins(12, 12, 12, 12)

        bg.addWidget(_lbl("Min product size (bp)"), 0, 0)
        self.product_min = QSpinBox(); self.product_min.setRange(50, 1000); self.product_min.setValue(100)
        bg.addWidget(self.product_min, 0, 1)
        bg.addWidget(_lbl("Max product size (bp)"), 0, 2)
        self.product_max = QSpinBox(); self.product_max.setRange(50, 2000); self.product_max.setValue(250)
        bg.addWidget(self.product_max, 0, 3)

        bg.addWidget(_lbl("Flank size (bp)"), 1, 0)
        self.flank = QSpinBox(); self.flank.setRange(20, 1000); self.flank.setValue(100)
        bg.addWidget(self.flank, 1, 1)
        bg.addWidget(_lbl("Pairs per SSR"), 1, 2)
        self.num_pairs = QSpinBox(); self.num_pairs.setRange(1, 5); self.num_pairs.setValue(2)
        bg.addWidget(self.num_pairs, 1, 3)

        bg.addWidget(_lbl("Min primer size (bp)"), 2, 0)
        self.primer_min_size = QSpinBox(); self.primer_min_size.setRange(10, 30); self.primer_min_size.setValue(18)
        bg.addWidget(self.primer_min_size, 2, 1)
        bg.addWidget(_lbl("Optimal primer size (bp)"), 2, 2)
        self.primer_opt_size = QSpinBox(); self.primer_opt_size.setRange(10, 35); self.primer_opt_size.setValue(22)
        bg.addWidget(self.primer_opt_size, 2, 3)
        bg.addWidget(_lbl("Max primer size (bp)"), 3, 0)
        self.primer_max_size = QSpinBox(); self.primer_max_size.setRange(15, 40); self.primer_max_size.setValue(25)
        bg.addWidget(self.primer_max_size, 3, 1)

        bg.addWidget(_lbl("Min Tm (°C)"), 4, 0)
        self.tm_min = QDoubleSpinBox(); self.tm_min.setRange(40, 75); self.tm_min.setValue(52); self.tm_min.setSingleStep(0.5)
        bg.addWidget(self.tm_min, 4, 1)
        bg.addWidget(_lbl("Optimal Tm (°C)"), 4, 2)
        self.tm_opt = QDoubleSpinBox(); self.tm_opt.setRange(40, 75); self.tm_opt.setValue(58); self.tm_opt.setSingleStep(0.5)
        bg.addWidget(self.tm_opt, 4, 3)
        bg.addWidget(_lbl("Max Tm (°C)"), 5, 0)
        self.tm_max = QDoubleSpinBox(); self.tm_max.setRange(40, 75); self.tm_max.setValue(60); self.tm_max.setSingleStep(0.5)
        bg.addWidget(self.tm_max, 5, 1)

        bg.addWidget(_lbl("Min GC (%)"), 6, 0)
        self.gc_min = QDoubleSpinBox(); self.gc_min.setRange(0, 100); self.gc_min.setValue(30)
        bg.addWidget(self.gc_min, 6, 1)
        bg.addWidget(_lbl("Max GC (%)"), 6, 2)
        self.gc_max = QDoubleSpinBox(); self.gc_max.setRange(0, 100); self.gc_max.setValue(60)
        bg.addWidget(self.gc_max, 6, 3)
        bg.setColumnStretch(4, 1)
        tabs.addTab(basic, "Basic Settings")

        adv = QWidget()
        ag = QGridLayout(adv); ag.setSpacing(10); ag.setContentsMargins(12, 12, 12, 12)
        ag.addWidget(_lbl("Max poly-X"), 0, 0)
        self.max_poly_x = QSpinBox(); self.max_poly_x.setRange(1, 10); self.max_poly_x.setValue(4)
        ag.addWidget(self.max_poly_x, 0, 1)
        ag.addWidget(_lbl("Max self-complementarity"), 0, 2)
        self.max_self_any = QDoubleSpinBox(); self.max_self_any.setRange(0, 30); self.max_self_any.setValue(8)
        ag.addWidget(self.max_self_any, 0, 3)
        ag.addWidget(_lbl("Max 3' self-complementarity"), 1, 0)
        self.max_self_end = QDoubleSpinBox(); self.max_self_end.setRange(0, 20); self.max_self_end.setValue(3)
        ag.addWidget(self.max_self_end, 1, 1)
        ag.addWidget(_lbl("Max pair complementarity"), 1, 2)
        self.max_pair_any = QDoubleSpinBox(); self.max_pair_any.setRange(0, 30); self.max_pair_any.setValue(8)
        ag.addWidget(self.max_pair_any, 1, 3)
        ag.addWidget(_lbl("Max pair 3' complementarity"), 2, 0)
        self.max_pair_end = QDoubleSpinBox(); self.max_pair_end.setRange(0, 20); self.max_pair_end.setValue(3)
        ag.addWidget(self.max_pair_end, 2, 1)
        ag.addWidget(_lbl("Max hairpin Tm (°C)"), 2, 2)
        self.max_hairpin = QDoubleSpinBox(); self.max_hairpin.setRange(0, 60); self.max_hairpin.setValue(24)
        ag.addWidget(self.max_hairpin, 2, 3)
        ag.addWidget(_lbl("Max Tm difference (°C)"), 3, 0)
        self.max_tm_diff = QDoubleSpinBox(); self.max_tm_diff.setRange(0, 10); self.max_tm_diff.setValue(2)
        ag.addWidget(self.max_tm_diff, 3, 1)
        ag.setColumnStretch(4, 1)
        tabs.addTab(adv, "Advanced Primer3 Parameters")
        L.addWidget(tabs)

        # Run
        run_group = QGroupBox("Run")
        rg = QVBoxLayout(run_group)
        self.ssr_source_label = QLabel("")
        self.ssr_source_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        rg.addWidget(self.ssr_source_label)
        run_row = QHBoxLayout()
        self.run_all_btn = QPushButton("Design for all SSRs")
        self.run_all_btn.setObjectName("primary")
        self.run_all_btn.clicked.connect(lambda: self._run("all"))
        self.run_sel_btn = QPushButton("Design for selected SSRs")
        self.run_sel_btn.clicked.connect(lambda: self._run("selected"))
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        run_row.addWidget(self.run_all_btn); run_row.addWidget(self.run_sel_btn)
        run_row.addWidget(self.cancel_btn); run_row.addStretch()
        rg.addLayout(run_row)
        self.run_status = QLabel("")
        self.run_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        rg.addWidget(self.run_status)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False); self.progress_bar.setFixedHeight(6)
        rg.addWidget(self.progress_bar)
        L.addWidget(run_group)

        # Results
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        res_layout = QVBoxLayout(self.results_group)

        self.metrics_label = QLabel("")
        self.metrics_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        res_layout.addWidget(self.metrics_label)

        filter_group = QGroupBox("Quality Filters")
        fg = QGridLayout(filter_group); fg.setSpacing(8); fg.setContentsMargins(10, 10, 10, 10)
        fg.addWidget(_lbl("Min 3' stability (kcal/mol)"), 0, 0)
        self.dg_min = QDoubleSpinBox(); self.dg_min.setRange(0, 20); self.dg_min.setValue(0); self.dg_min.setSingleStep(0.5)
        self.dg_min.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.dg_min, 0, 1)
        fg.addWidget(_lbl("Max 3' stability (kcal/mol)"), 0, 2)
        self.dg_max = QDoubleSpinBox(); self.dg_max.setRange(0, 20); self.dg_max.setValue(20); self.dg_max.setSingleStep(0.5)
        self.dg_max.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.dg_max, 0, 3)
        fg.addWidget(_lbl("Min GC clamp"), 1, 0)
        self.clamp_min = QSpinBox(); self.clamp_min.setRange(0, 5); self.clamp_min.setValue(0)
        self.clamp_min.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.clamp_min, 1, 1)
        fg.addWidget(_lbl("Max GC clamp"), 1, 2)
        self.clamp_max = QSpinBox(); self.clamp_max.setRange(0, 5); self.clamp_max.setValue(5)
        self.clamp_max.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.clamp_max, 1, 3)
        self.filter_status = QLabel("")
        self.filter_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        fg.addWidget(self.filter_status, 2, 0, 1, 4)
        fg.setColumnStretch(4, 1)
        res_layout.addWidget(filter_group)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        res_layout.addWidget(self.table)

        dl_row = QHBoxLayout()
        self.download_csv_btn = QPushButton("Download primer table (CSV)")
        self.download_csv_btn.clicked.connect(self._download_csv)
        self.download_fasta_btn = QPushButton("Download BLAST-ready FASTA")
        self.download_fasta_btn.clicked.connect(self._download_fasta)
        dl_row.addWidget(self.download_csv_btn)
        dl_row.addWidget(self.download_fasta_btn)
        dl_row.addStretch()
        res_layout.addLayout(dl_row)

        L.addWidget(self.results_group)
        L.addStretch()
        self._refresh()

    def _run(self, mode):
        if not self.state.has_genome:
            self._set_status("Load a genome first", WARNING); return
        if not self.state.has_ssrs:
            self._set_status("Run SSR detection first", WARNING); return
        if mode == "selected":
            ssr_list = self.state.selected_ssrs
            if not ssr_list:
                self._set_status("No SSRs selected", WARNING); return
        else:
            ssr_list = self.state.ssrs

        if len(ssr_list) > PRIMER_DESIGN_WARN_THRESHOLD:
            from PyQt6.QtWidgets import QMessageBox
            msg = QMessageBox(self)
            msg.setWindowTitle("Large Dataset Warning")
            msg.setText(
                f"You are about to design primers for {len(ssr_list):,} SSRs.\n\n"
                f"This may take a long time.\n"
                f"Consider selecting a subset on the SSR Detection page first."
            )
            msg.addButton("Continue anyway", QMessageBox.ButtonRole.AcceptRole)
            cancel = msg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            msg.exec()
            if msg.clickedButton() == cancel:
                return

        primer_opts = {
            "PRIMER_MIN_SIZE":           self.primer_min_size.value(),
            "PRIMER_OPT_SIZE":           self.primer_opt_size.value(),
            "PRIMER_MAX_SIZE":           self.primer_max_size.value(),
            "PRIMER_MIN_TM":             self.tm_min.value(),
            "PRIMER_OPT_TM":             self.tm_opt.value(),
            "PRIMER_MAX_TM":             self.tm_max.value(),
            "PRIMER_MIN_GC":             self.gc_min.value(),
            "PRIMER_MAX_GC":             self.gc_max.value(),
            "PRIMER_MAX_POLY_X":         self.max_poly_x.value(),
            "PRIMER_MAX_SELF_ANY":       self.max_self_any.value(),
            "PRIMER_MAX_SELF_END":       self.max_self_end.value(),
            "PRIMER_PAIR_MAX_COMPL_ANY": self.max_pair_any.value(),
            "PRIMER_PAIR_MAX_COMPL_END": self.max_pair_end.value(),
            "PRIMER_MAX_HAIRPIN_TH":     self.max_hairpin.value(),
            "PRIMER_PAIR_MAX_DIFF_TM":   self.max_tm_diff.value(),
        }
        params = {
            "flank":        self.flank.value(),
            "product_min":  self.product_min.value(),
            "product_max":  self.product_max.value(),
            "num_pairs":    self.num_pairs.value(),
            "primer_opts":  primer_opts,
            "preset":       "amplicon" if self._amplicon_mode else "recommended",
            "amplicon_mode": self._amplicon_mode,
        }
        self.state.product_min = params["product_min"]
        self.state.product_max = params["product_max"]
        self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        mode_tag = " [Amplicon]" if self._amplicon_mode else " [Capillary]"
        self._set_status(f"Designing primers for {len(ssr_list):,} SSRs{mode_tag}...", TEXT_SECONDARY)
        self.mw.set_status(f"Running Primer3{mode_tag}...")
        self._worker = PrimerWorker(self.state.genome, ssr_list, params)
        self._worker.progress.connect(self._on_progress, Qt.ConnectionType.QueuedConnection)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Designing: {done:,} / {total:,}", TEXT_SECONDARY)
        QApplication.processEvents()

    def _on_done(self, res, elapsed):
        primers   = res["success"]
        n_failed  = len(res["failed"])
        n_skipped = len(res.get("skipped", []))
        self.state.clear_downstream_of_primers()
        self.state.primer_results          = primers
        self.state.filtered_primer_results = primers
        self.state.primer_version         += 1
        self.state.workflow_mode           = "capillary" if not self._amplicon_mode else "amplicon"

        if primers:
            import pandas as pd
            df = pd.DataFrame(primers)
            if "left_3end_dg" in df.columns:
                all_dg = pd.concat([df["left_3end_dg"], df["right_3end_dg"]]).dropna()
                self._filter_updating = True
                for w in [self.dg_min, self.dg_max, self.clamp_min, self.clamp_max]:
                    w.blockSignals(True)
                self.dg_min.setValue(0.0)
                self.dg_max.setValue(min(20.0, round(float(all_dg.max()) + 1, 1)))
                self.clamp_min.setValue(0)
                self.clamp_max.setValue(5)
                for w in [self.dg_min, self.dg_max, self.clamp_min, self.clamp_max]:
                    w.blockSignals(False)
                self._filter_updating = False

        self._apply_filters()
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False); self.progress_bar.setVisible(False)
        msg = f"{len(primers):,} primer pairs designed in {elapsed:.1f}s"
        if n_skipped: msg += f" — {n_skipped:,} skipped"
        if n_failed:  msg += f" — {n_failed:,} failed"
        self._set_status(msg, SUCCESS)
        self.mw.set_status(msg); self.mw.on_step_complete(3)
        self._cleanup_worker(); self._refresh()

    def _on_error(self, msg):
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False); self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", ERROR)
        self.mw.set_status("Primer design failed")
        self._cleanup_worker()

    def _cancel(self):
        if self._worker:
            self._worker.terminate(); self._worker.wait(); self._cleanup_worker()
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False); self.progress_bar.setVisible(False)
        self._set_status("Cancelled", WARNING)
        self.mw.set_status("Primer design cancelled")

    def _cleanup_worker(self):
        if self._worker:
            try:
                self._worker.progress.disconnect()
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            self._worker = None

    def _set_status(self, msg, color):
        self.run_status.setText(msg)
        self.run_status.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    def _on_filter_changed(self):
        if self._filter_updating: return
        self._filter_updating = True
        try:
            self._apply_filters()
            self._populate_table()
        finally:
            self._filter_updating = False

    def _apply_filters(self):
        if not self.state.primer_results:
            return
        import pandas as pd

        try:
            dg_min   = float(self.dg_min.value())
            dg_max   = float(self.dg_max.value())
            clamp_lo = int(self.clamp_min.value())
            clamp_hi = int(self.clamp_max.value())
        except Exception:
            return

        if dg_min > dg_max:
            dg_min, dg_max = dg_max, dg_min
        if clamp_lo > clamp_hi:
            clamp_lo, clamp_hi = clamp_hi, clamp_lo

        df = pd.DataFrame(self.state.primer_results)
        total = len(df)

        # 3' stability filter
        if "left_3end_dg" in df.columns and "right_3end_dg" in df.columns:
            mask = (
                df["left_3end_dg"].between(dg_min, dg_max) &
                df["right_3end_dg"].between(dg_min, dg_max)
            )
            df = df[mask].copy()

        # GC clamp filter (only if primers exist and DataFrame not empty)
        if "left_primer" in df.columns and "right_primer" in df.columns and not df.empty:
            def _gc_clamp(seq):
                try:
                    s = str(seq) if seq is not None else ""
                    return sum(1 for b in s[-5:].upper() if b in "GC")
                except Exception:
                    return 0

            # Use .loc to avoid SettingWithCopyWarning
            df = df.copy()
            df["left_gc_clamp"] = df["left_primer"].apply(_gc_clamp)
            df["right_gc_clamp"] = df["right_primer"].apply(_gc_clamp)

            mask = (
                df["left_gc_clamp"].between(clamp_lo, clamp_hi) &
                df["right_gc_clamp"].between(clamp_lo, clamp_hi)
            )
            df = df[mask].copy()

        # Update filtered results
        if df.empty:
            self.state.filtered_primer_results = self.state.primer_results
            self.filter_status.setText("No primers pass current filters — using full set for BLAST")
            self.filter_status.setStyleSheet(f"color: {WARNING}; font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.state.filtered_primer_results = df.to_dict(orient="records")
            removed = total - len(df)
            if removed > 0:
                self.filter_status.setText(f"{removed:,} pairs removed — {len(df):,} remaining")
                self.filter_status.setStyleSheet(f"color: {WARNING}; font-size: {FONT_SIZE_SMALL}pt;")
            else:
                self.filter_status.setText(f"All {len(df):,} pairs pass filters")
                self.filter_status.setStyleSheet(f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;")

    def _refresh(self):
        if not self.state.has_genome or not self.state.has_ssrs:
            self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        else:
            self.run_all_btn.setEnabled(True)
            n_sel = len(self.state.selected_ssrs) if self.state.selected_ssrs else 0
            self.ssr_source_label.setText(f"{len(self.state.ssrs):,} SSRs available" + (f" — {n_sel:,} selected" if n_sel else ""))
        if not self.state.has_primers:
            self.results_group.setVisible(False); return
        if self.state.primer_version == self._last_rendered_version: return
        self._last_rendered_version = self.state.primer_version
        self.results_group.setVisible(True)
        self._apply_filters(); self._populate_table()

    def _populate_table(self):
        if not self.state.has_primers: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        total = len(self.state.primer_results)
        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df
        metrics = f"{total:,} primer pairs designed | {len(primers):,} pass filters | {df['ssr_id'].nunique():,} SSRs covered"
        if truncated: metrics += f" | showing first {TABLE_DISPLAY_LIMIT:,}"
        self.metrics_label.setText(metrics)

        COLS = {
            "ssr_id": "SSR ID", "pair_rank": "Pair", "contig": "Contig",
            "motif": "Motif", "left_primer": "Forward", "right_primer": "Reverse",
            "product_size": "Size", "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC", "right_gc": "Rev GC",
        }

        # Add Chromosome column if mapping exists
        show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
        if show_chromosome:
            COLS["chromosome"] = "Chromosome"

        display_cols = [c for c in COLS if c in display_df.columns]

        # Insert Chromosome after Contig
        if show_chromosome and "chromosome" in display_cols and "contig" in display_cols:
            display_cols.remove("chromosome")
            contig_idx = display_cols.index("contig")
            display_cols.insert(contig_idx + 1, "chromosome")

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df))
        self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([COLS.get(c, c) for c in display_cols])

        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        for ci, col in enumerate(display_cols):
            if col == "contig":
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(ci, 140)
            elif col == "chromosome":
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(ci, 120)
            elif col in ("left_primer", "right_primer"):
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(ci, 200)
            else:
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.ResizeToContents)

        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                if col == "chromosome":
                    contig = display_df.iat[row_idx, display_df.columns.get_loc("contig")]
                    val = self.state.get_display_name(contig)
                else:
                    val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                text = f"{val:.2f}" if isinstance(val, float) else ("" if val is None else str(val))
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(text))
        self.table.setSortingEnabled(True)

    def _download_csv(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save primer table", "primer_table.csv", "CSV files (*.csv)")
        if not path: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        if hasattr(self.state, 'get_display_name'):
            df["Chromosome"] = df["contig"].apply(lambda c: self.state.get_display_name(c))
        df.rename(columns={
            "ssr_id": "SSR ID", "pair_rank": "Pair", "contig": "Contig",
            "motif": "Motif", "repeat_count": "Repeat count",
            "left_primer": "Forward", "right_primer": "Reverse",
            "product_size": "Size", "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC", "right_gc": "Rev GC",
            "left_3end_dg": "Fwd 3' ΔG", "right_3end_dg": "Rev 3' ΔG",
        }).to_csv(path, index=False, encoding="utf-8-sig")
        self.mw.set_status(f"Saved to {path}")

    def _download_fasta(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save BLAST FASTA", "primers_blast_ready.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        from core.primer_design import primers_to_blast_fasta
        primers = self.state.filtered_primer_results or self.state.primer_results
        with open(path, "w") as f: f.write(primers_to_blast_fasta(primers))
        self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()