"""
primer_panel.py — Primer Design
Single scroll area. Clean layout. Filters appear after design.
Supports Standard and GBS (Genotype-by-Sequencing) modes.
"""
import os, sys, time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog, QTabWidget,
    QCheckBox, QFrame,
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


# ---------------------------------------------------------------------------
# Mode toggle button — big obvious clickable card
# ---------------------------------------------------------------------------

class ModeButton(QFrame):
    """A large clickable card that acts as a mode selector button."""

    def __init__(self, title: str, description: str, tag: str, parent=None):
        super().__init__(parent)
        self.tag      = tag
        self._active  = False
        self._title   = title
        self._desc    = description
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedHeight(80)
        self.setSizePolicy(
            __import__('PyQt6.QtWidgets', fromlist=['QSizePolicy']).QSizePolicy.Policy.Expanding,
            __import__('PyQt6.QtWidgets', fromlist=['QSizePolicy']).QSizePolicy.Policy.Fixed,
        )
        self._build()
        self._apply_style()

    def _build(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 10, 16, 10)
        layout.setSpacing(4)

        self._title_lbl = QLabel(self._title)
        self._title_lbl.setFont(QFont(FONT_UI, FONT_SIZE_NORMAL, QFont.Weight.Bold))

        self._desc_lbl = QLabel(self._desc)
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL - 1))

        layout.addWidget(self._title_lbl)
        layout.addWidget(self._desc_lbl)

    def _apply_style(self):
        if self._active:
            self.setStyleSheet(f"""
                ModeButton {{
                    background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
                        stop:0 {ACCENT}33, stop:1 {ACCENT}11);
                    border: 2px solid {ACCENT};
                    border-radius: 8px;
                }}
            """)
            self._title_lbl.setStyleSheet(f"color: {ACCENT}; background: transparent;")
            self._desc_lbl.setStyleSheet(f"color: {TEXT_PRIMARY}; background: transparent;")
        else:
            self.setStyleSheet(f"""
                ModeButton {{
                    background: {BG_MID};
                    border: 1px solid {BORDER};
                    border-radius: 8px;
                }}
                ModeButton:hover {{
                    border: 1px solid {ACCENT}88;
                    background: {BG_LIGHT};
                }}
            """)
            self._title_lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; background: transparent;")
            self._desc_lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; background: transparent;")

    def set_active(self, active: bool):
        self._active = active
        self._apply_style()

    def mousePressEvent(self, event):
        self.parent().parent()  # just to keep reference
        # Find ModeSelector parent and call select
        p = self.parent()
        while p is not None:
            if hasattr(p, '_on_mode_btn_clicked'):
                p._on_mode_btn_clicked(self.tag)
                break
            p = p.parent()
        super().mousePressEvent(event)


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

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
                gbs_mode=self.params["gbs_mode"],
                add_adapters=self.params["add_adapters"],
            )
            self.finished.emit(res, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class PrimerPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._worker = None
        self._last_rendered_version = -1
        self._gbs_mode = False
        self._build_ui()

    def _on_mode_btn_clicked(self, tag: str):
        """Called by ModeButton when clicked."""
        self._gbs_mode = (tag == "gbs")
        self._btn_standard.set_active(tag == "standard")
        self._btn_gbs.set_active(tag == "gbs")
        self._gbs_options.setVisible(self._gbs_mode)
        if self._gbs_mode:
            self.flank.setValue(50)
            self.product_min.setValue(80)
            self.product_max.setValue(200)
            self.tm_min.setValue(58.0)
            self.tm_opt.setValue(60.0)
            self.tm_max.setValue(62.0)
            self.max_tm_diff.setValue(1.0)
            self.max_poly_x.setValue(3)
        else:
            self.flank.setValue(100)
            self.product_min.setValue(100)
            self.product_max.setValue(250)
            self.tm_min.setValue(52.0)
            self.tm_opt.setValue(58.0)
            self.tm_max.setValue(60.0)
            self.max_tm_diff.setValue(2.0)
            self.max_poly_x.setValue(4)

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
        title = QLabel("Primer Design")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Design PCR primers flanking each SSR using Primer3")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # ── Mode selector — big card buttons ─────────────
        mode_group = QGroupBox("Design Mode")
        mg = QVBoxLayout(mode_group)
        mg.setSpacing(8)

        mode_row = QHBoxLayout()
        mode_row.setSpacing(12)

        self._btn_standard = ModeButton(
            "◉  Standard Mode",
            "Capillary electrophoresis / fragment analysis\n"
            "Product sizes 100–350bp · Flank 100bp",
            "standard",
        )
        self._btn_gbs = ModeButton(
            "◎  GBS Mode",
            "Genotype-by-Sequencing · Illumina amplicon sequencing\n"
            "Product sizes 80–200bp · Flank 50bp · Tighter Tm range",
            "gbs",
        )
        self._btn_standard.set_active(True)
        mode_row.addWidget(self._btn_standard)
        mode_row.addWidget(self._btn_gbs)
        mg.addLayout(mode_row)

        # GBS-specific options
        self._gbs_options = QWidget()
        gbs_row = QHBoxLayout(self._gbs_options)
        gbs_row.setContentsMargins(4, 4, 0, 0)
        self._add_adapters_cb = QCheckBox(
            "Append Illumina M13 adapter tails to primer sequences"
        )
        self._add_adapters_cb.setChecked(True)
        self._add_adapters_cb.setToolTip(
            "Adds M13F/M13R universal tails to the 5' end of each primer.\n"
            "Required for standard Illumina amplicon sequencing workflows.\n"
            "Bare sequences (without tails) are always used for BLAST."
        )
        gbs_row.addWidget(self._add_adapters_cb)
        gbs_row.addStretch()
        self._gbs_options.setVisible(False)
        mg.addWidget(self._gbs_options)
        L.addWidget(mode_group)

        # ── Settings tabs ─────────────────────────────────
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

        # ── Run ───────────────────────────────────────────
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

        # ── Results ───────────────────────────────────────
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
        self.download_gbs_btn = QPushButton("Download GBS FASTA (with adapter tails)")
        self.download_gbs_btn.clicked.connect(self._download_gbs_fasta)
        self.download_gbs_btn.setVisible(False)
        dl_row.addWidget(self.download_csv_btn)
        dl_row.addWidget(self.download_fasta_btn)
        dl_row.addWidget(self.download_gbs_btn)
        dl_row.addStretch()
        res_layout.addLayout(dl_row)

        L.addWidget(self.results_group)
        L.addStretch()
        self._refresh()

    # ---------------------------------------------------------
    # RUN
    # ---------------------------------------------------------
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
            "preset":       "gbs" if self._gbs_mode else "recommended",
            "gbs_mode":     self._gbs_mode,
            "add_adapters": self._add_adapters_cb.isChecked() if self._gbs_mode else False,
        }
        self.state.product_min = params["product_min"]
        self.state.product_max = params["product_max"]
        self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        mode_tag = " [GBS]" if self._gbs_mode else ""
        self._set_status(f"Designing primers for {len(ssr_list):,} SSRs{mode_tag}...", TEXT_SECONDARY)
        self.mw.set_status(f"Running Primer3{mode_tag}...")
        self._worker = PrimerWorker(self.state.genome, ssr_list, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Designing: {done:,} / {total:,}", TEXT_SECONDARY)

    def _on_done(self, res, elapsed):
        primers   = res["success"]
        n_failed  = len(res["failed"])
        n_skipped = len(res.get("skipped", []))
        self.state.clear_downstream_of_primers()
        self.state.primer_results          = primers
        self.state.filtered_primer_results = primers
        self.state.primer_version         += 1
        self.state.gbs_mode                = self._gbs_mode

        if primers:
            import pandas as pd
            df = pd.DataFrame(primers)
            if "left_3end_dg" in df.columns:
                all_dg = pd.concat([df["left_3end_dg"], df["right_3end_dg"]]).dropna()
                for w in [self.dg_min, self.dg_max, self.clamp_min, self.clamp_max]:
                    w.blockSignals(True)
                self.dg_min.setValue(0.0)
                self.dg_max.setValue(min(20.0, round(float(all_dg.max()) + 1, 1)))
                self.clamp_min.setValue(0); self.clamp_max.setValue(5)
                for w in [self.dg_min, self.dg_max, self.clamp_min, self.clamp_max]:
                    w.blockSignals(False)

        self._apply_filters()
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False); self.progress_bar.setVisible(False)
        msg = f"{len(primers):,} primer pairs designed in {elapsed:.1f}s"
        if n_skipped: msg += f" — {n_skipped:,} skipped (low-complexity)"
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
            self._worker.genome = None; self._worker.ssr_list = None
            self._worker = None

    def _set_status(self, msg, color):
        self.run_status.setText(msg)
        self.run_status.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # FILTERS
    # ---------------------------------------------------------
    def _on_filter_changed(self):
        self._apply_filters(); self._populate_table()

    def _apply_filters(self):
        if not self.state.primer_results: return
        import pandas as pd
        df = pd.DataFrame(self.state.primer_results)
        if "left_3end_dg" in df.columns:
            dg_min = self.dg_min.value(); dg_max = self.dg_max.value()
            clamp_lo = self.clamp_min.value(); clamp_hi = self.clamp_max.value()
            mask = df["left_3end_dg"].between(dg_min, dg_max) & df["right_3end_dg"].between(dg_min, dg_max)
            df = df[mask].copy()
            def _gc(seq): return sum(1 for b in seq[-5:].upper() if b in "GC")
            df["left_gc_clamp"]  = df["left_primer"].apply(_gc)
            df["right_gc_clamp"] = df["right_primer"].apply(_gc)
            df = df[df["left_gc_clamp"].between(clamp_lo, clamp_hi) & df["right_gc_clamp"].between(clamp_lo, clamp_hi)].copy()
        total = len(self.state.primer_results)
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

    # ---------------------------------------------------------
    # TABLE
    # ---------------------------------------------------------
    def _refresh(self):
        if not self.state.has_genome or not self.state.has_ssrs:
            self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        else:
            self.run_all_btn.setEnabled(True)
            n_sel = len(self.state.selected_ssrs) if self.state.selected_ssrs else 0
            self.ssr_source_label.setText(
                f"{len(self.state.ssrs):,} SSRs available" + (f" — {n_sel:,} selected" if n_sel else "")
            )
        if not self.state.has_primers:
            self.results_group.setVisible(False); return
        if self.state.primer_version == self._last_rendered_version: return
        self._last_rendered_version = self.state.primer_version
        self.results_group.setVisible(True)
        has_tails = self.state.primer_results and "left_primer_tailed" in self.state.primer_results[0]
        self.download_gbs_btn.setVisible(has_tails)
        self._apply_filters(); self._populate_table()

    def _populate_table(self):
        if not self.state.has_primers: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        total = len(self.state.primer_results)
        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df
        metrics = (
            f"{total:,} primer pairs designed   |   {len(primers):,} pass quality filters   |   "
            f"{df['ssr_id'].nunique():,} SSRs covered"
        )
        if truncated:
            metrics += f"   |   showing first {TABLE_DISPLAY_LIMIT:,} — download CSV for full results"
        self.metrics_label.setText(metrics)
        COLS = {
            "ssr_id": "SSR ID", "pair_rank": "Pair rank", "contig": "Contig",
            "motif": "Motif", "left_primer": "Forward primer", "right_primer": "Reverse primer",
            "product_size": "Product (bp)", "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC%", "right_gc": "Rev GC%",
            "left_3end_dg": "Fwd 3' ΔG", "right_3end_dg": "Rev 3' ΔG",
            "genomic_feature": "Feature",
        }
        display_cols = [c for c in COLS if c in display_df.columns]
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df)); self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([COLS[c] for c in display_cols])
        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        for ci, col in enumerate(display_cols):
            if col == "contig":
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(ci, 160)
            elif col in ("left_primer", "right_primer"):
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(ci, 200)
            else:
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.ResizeToContents)
        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                text = f"{val:.2f}" if isinstance(val, float) else ("" if val is None else str(val))
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(text))
        self.table.setSortingEnabled(True)

    # ---------------------------------------------------------
    # DOWNLOADS
    # ---------------------------------------------------------
    def _download_csv(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save primer table", "primer_table.csv", "CSV files (*.csv)")
        if not path: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        pd.DataFrame(primers).rename(columns={
            "ssr_id": "SSR ID", "pair_rank": "Pair rank", "contig": "Contig",
            "motif": "Motif", "repeat_count": "Repeat count",
            "left_primer": "Forward primer", "right_primer": "Reverse primer",
            "product_size": "Product size (bp)", "left_tm": "Forward Tm (°C)",
            "right_tm": "Reverse Tm (°C)", "left_gc": "Forward GC (%)",
            "right_gc": "Reverse GC (%)", "left_3end_dg": "Forward 3' stability (kcal/mol)",
            "right_3end_dg": "Reverse 3' stability (kcal/mol)",
            "genomic_feature": "Genomic feature",
            "left_primer_tailed": "Forward primer (with adapter tail)",
            "right_primer_tailed": "Reverse primer (with adapter tail)",
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

    def _download_gbs_fasta(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save GBS FASTA", "primers_gbs_tailed.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        from core.primer_design import primers_to_gbs_fasta
        primers = self.state.filtered_primer_results or self.state.primer_results
        with open(path, "w") as f: f.write(primers_to_gbs_fasta(primers))
        self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()