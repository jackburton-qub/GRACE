"""
primer_panel.py — Primer Design
Single scroll area. Clean layout. Filters appear after design.
"""
import os, sys, time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog, QTabWidget,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)

# Cap table display to prevent main-thread freeze on large primer sets.
# Full data is always in state and available for export/BLAST.
TABLE_DISPLAY_LIMIT = 10_000


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
        self.genome = genome; self.ssr_list = ssr_list; self.params = params

    def run(self):
        try:
            _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            if _root not in sys.path: sys.path.insert(0, _root)
            from core.primer_design import design_primers_for_all_ssrs
            start = time.time()
            res = design_primers_for_all_ssrs(
                genome=self.genome, ssr_list=self.ssr_list,
                flank=self.params["flank"],
                product_size_range=(self.params["product_min"], self.params["product_max"]),
                preset="recommended", primer_opts=self.params["primer_opts"],
                num_pairs=self.params["num_pairs"],
                progress_callback=lambda d, t: self.progress.emit(d, t),
            )
            self.finished.emit(res, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


class PrimerPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state; self.mw = main_window; self._worker = None
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
        title = QLabel("Primer Design")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Design PCR primers flanking each SSR using Primer3")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title); L.addWidget(sub)

        # Settings tabs
        tabs = QTabWidget()

        # Basic tab
        basic = QWidget()
        bg = QGridLayout(basic); bg.setSpacing(10); bg.setContentsMargins(12,12,12,12)
        bg.addWidget(_lbl("Min product size (bp)", "Minimum amplicon size.\nTypically 100–150 bp for capillary electrophoresis."), 0, 0)
        self.product_min = QSpinBox(); self.product_min.setRange(50,1000); self.product_min.setValue(100)
        bg.addWidget(self.product_min, 0, 1)
        bg.addWidget(_lbl("Max product size (bp)", "Maximum amplicon size.\nKeep under 350 bp for fragment analysis."), 0, 2)
        self.product_max = QSpinBox(); self.product_max.setRange(50,2000); self.product_max.setValue(250)
        bg.addWidget(self.product_max, 0, 3)
        bg.addWidget(_lbl("Flank size (bp)", "Flanking sequence provided to Primer3 on each side of the SSR.\n100 bp is usually sufficient."), 1, 0)
        self.flank = QSpinBox(); self.flank.setRange(50,1000); self.flank.setValue(100)
        bg.addWidget(self.flank, 1, 1)
        bg.addWidget(_lbl("Pairs per SSR", "Number of primer pairs per SSR.\nPair 0 is Primer3's top-ranked choice."), 1, 2)
        self.num_pairs = QSpinBox(); self.num_pairs.setRange(1,5); self.num_pairs.setValue(2)
        bg.addWidget(self.num_pairs, 1, 3)
        bg.addWidget(_lbl("Min primer size (bp)"), 2, 0)
        self.primer_min_size = QSpinBox(); self.primer_min_size.setRange(10,30); self.primer_min_size.setValue(18)
        bg.addWidget(self.primer_min_size, 2, 1)
        bg.addWidget(_lbl("Optimal primer size (bp)"), 2, 2)
        self.primer_opt_size = QSpinBox(); self.primer_opt_size.setRange(10,35); self.primer_opt_size.setValue(22)
        bg.addWidget(self.primer_opt_size, 2, 3)
        bg.addWidget(_lbl("Max primer size (bp)"), 3, 0)
        self.primer_max_size = QSpinBox(); self.primer_max_size.setRange(15,40); self.primer_max_size.setValue(25)
        bg.addWidget(self.primer_max_size, 3, 1)
        bg.addWidget(_lbl("Min Tm (°C)", "Minimum melting temperature.\n52–60°C is standard for SSR genotyping."), 4, 0)
        self.tm_min = QDoubleSpinBox(); self.tm_min.setRange(40,75); self.tm_min.setValue(52); self.tm_min.setSingleStep(0.5)
        bg.addWidget(self.tm_min, 4, 1)
        bg.addWidget(_lbl("Optimal Tm (°C)"), 4, 2)
        self.tm_opt = QDoubleSpinBox(); self.tm_opt.setRange(40,75); self.tm_opt.setValue(58); self.tm_opt.setSingleStep(0.5)
        bg.addWidget(self.tm_opt, 4, 3)
        bg.addWidget(_lbl("Max Tm (°C)"), 5, 0)
        self.tm_max = QDoubleSpinBox(); self.tm_max.setRange(40,75); self.tm_max.setValue(60); self.tm_max.setSingleStep(0.5)
        bg.addWidget(self.tm_max, 5, 1)
        bg.addWidget(_lbl("Min GC (%)", "Minimum GC content. 30–65% is standard."), 6, 0)
        self.gc_min = QDoubleSpinBox(); self.gc_min.setRange(0,100); self.gc_min.setValue(30)
        bg.addWidget(self.gc_min, 6, 1)
        bg.addWidget(_lbl("Max GC (%)"), 6, 2)
        self.gc_max = QDoubleSpinBox(); self.gc_max.setRange(0,100); self.gc_max.setValue(60)
        bg.addWidget(self.gc_max, 6, 3)
        bg.setColumnStretch(4, 1)
        tabs.addTab(basic, "Basic Settings")

        # Advanced tab
        adv = QWidget()
        ag = QGridLayout(adv); ag.setSpacing(10); ag.setContentsMargins(12,12,12,12)
        ag.addWidget(_lbl("Max poly-X", "Maximum run of a single nucleotide within a primer."), 0, 0)
        self.max_poly_x = QSpinBox(); self.max_poly_x.setRange(1,10); self.max_poly_x.setValue(4)
        ag.addWidget(self.max_poly_x, 0, 1)
        ag.addWidget(_lbl("Max self-complementarity", "Maximum self-complementarity score (any position)."), 0, 2)
        self.max_self_any = QDoubleSpinBox(); self.max_self_any.setRange(0,30); self.max_self_any.setValue(8)
        ag.addWidget(self.max_self_any, 0, 3)
        ag.addWidget(_lbl("Max 3' self-complementarity", "Keep low to prevent primer dimers."), 1, 0)
        self.max_self_end = QDoubleSpinBox(); self.max_self_end.setRange(0,20); self.max_self_end.setValue(3)
        ag.addWidget(self.max_self_end, 1, 1)
        ag.addWidget(_lbl("Max pair complementarity"), 1, 2)
        self.max_pair_any = QDoubleSpinBox(); self.max_pair_any.setRange(0,30); self.max_pair_any.setValue(8)
        ag.addWidget(self.max_pair_any, 1, 3)
        ag.addWidget(_lbl("Max pair 3' complementarity", "Keep low to prevent primer dimer extension."), 2, 0)
        self.max_pair_end = QDoubleSpinBox(); self.max_pair_end.setRange(0,20); self.max_pair_end.setValue(3)
        ag.addWidget(self.max_pair_end, 2, 1)
        ag.addWidget(_lbl("Max hairpin Tm (°C)"), 2, 2)
        self.max_hairpin = QDoubleSpinBox(); self.max_hairpin.setRange(0,60); self.max_hairpin.setValue(24)
        ag.addWidget(self.max_hairpin, 2, 3)
        ag.addWidget(_lbl("Max Tm difference (°C)", "Maximum Tm difference between forward and reverse primers."), 3, 0)
        self.max_tm_diff = QDoubleSpinBox(); self.max_tm_diff.setRange(0,10); self.max_tm_diff.setValue(2)
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
        self.run_all_btn.setToolTip("Design primers for every SSR detected.")
        self.run_sel_btn = QPushButton("Design for selected SSRs")
        self.run_sel_btn.clicked.connect(lambda: self._run("selected"))
        self.run_sel_btn.setToolTip("Design primers only for SSRs selected in the SSR Detection table.")
        run_row.addWidget(self.run_all_btn); run_row.addWidget(self.run_sel_btn); run_row.addStretch()
        rg.addLayout(run_row)
        self.run_status = QLabel("")
        self.run_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        rg.addWidget(self.run_status)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False); self.progress_bar.setFixedHeight(6)
        rg.addWidget(self.progress_bar)
        L.addWidget(run_group)

        # Results — hidden until design runs
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        res_layout = QVBoxLayout(self.results_group)

        self.metrics_label = QLabel("")
        self.metrics_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        res_layout.addWidget(self.metrics_label)

        # Quality filters — inside results, hidden until design
        filter_group = QGroupBox("Quality Filters")
        fg = QGridLayout(filter_group); fg.setSpacing(8); fg.setContentsMargins(10,10,10,10)
        fg.addWidget(_lbl("Min 3' stability (kcal/mol)",
            "Minimum 3' end stability (ΔG of last 4 bases).\nHigher = stronger 3' binding. Set to 0 to disable."), 0, 0)
        self.dg_min = QDoubleSpinBox(); self.dg_min.setRange(0,20); self.dg_min.setValue(0); self.dg_min.setSingleStep(0.5)
        self.dg_min.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.dg_min, 0, 1)
        fg.addWidget(_lbl("Max 3' stability (kcal/mol)", "Maximum 3' end stability."), 0, 2)
        self.dg_max = QDoubleSpinBox(); self.dg_max.setRange(0,20); self.dg_max.setValue(20); self.dg_max.setSingleStep(0.5)
        self.dg_max.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.dg_max, 0, 3)
        fg.addWidget(_lbl("Min GC clamp", "Minimum G/C bases in last 5 positions. 1–2 is optimal."), 1, 0)
        self.clamp_min = QSpinBox(); self.clamp_min.setRange(0,5); self.clamp_min.setValue(0)
        self.clamp_min.valueChanged.connect(self._on_filter_changed)
        fg.addWidget(self.clamp_min, 1, 1)
        fg.addWidget(_lbl("Max GC clamp", "Maximum G/C bases in last 5 positions. More than 3 increases mispriming risk."), 1, 2)
        self.clamp_max = QSpinBox(); self.clamp_max.setRange(0,5); self.clamp_max.setValue(5)
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
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        res_layout.addWidget(self.table)

        dl_row = QHBoxLayout()
        self.download_csv_btn = QPushButton("Download primer table (CSV)")
        self.download_csv_btn.clicked.connect(self._download_csv)
        self.download_fasta_btn = QPushButton("Download BLAST-ready FASTA")
        self.download_fasta_btn.clicked.connect(self._download_fasta)
        dl_row.addWidget(self.download_csv_btn); dl_row.addWidget(self.download_fasta_btn); dl_row.addStretch()
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
                self._set_status("No SSRs selected — use 'Design for all' or select rows on SSR Detection", WARNING); return
        else:
            ssr_list = self.state.ssrs

        primer_opts = {
            "PRIMER_MIN_SIZE": self.primer_min_size.value(),
            "PRIMER_OPT_SIZE": self.primer_opt_size.value(),
            "PRIMER_MAX_SIZE": self.primer_max_size.value(),
            "PRIMER_MIN_TM": self.tm_min.value(), "PRIMER_OPT_TM": self.tm_opt.value(), "PRIMER_MAX_TM": self.tm_max.value(),
            "PRIMER_MIN_GC": self.gc_min.value(), "PRIMER_MAX_GC": self.gc_max.value(),
            "PRIMER_MAX_POLY_X": self.max_poly_x.value(),
            "PRIMER_MAX_SELF_ANY": self.max_self_any.value(), "PRIMER_MAX_SELF_END": self.max_self_end.value(),
            "PRIMER_PAIR_MAX_COMPL_ANY": self.max_pair_any.value(), "PRIMER_PAIR_MAX_COMPL_END": self.max_pair_end.value(),
            "PRIMER_MAX_HAIRPIN_TH": self.max_hairpin.value(), "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff.value(),
        }
        params = {
            "flank": self.flank.value(), "product_min": self.product_min.value(),
            "product_max": self.product_max.value(), "num_pairs": self.num_pairs.value(),
            "primer_opts": primer_opts,
        }
        self.state.product_min = params["product_min"]
        self.state.product_max = params["product_max"]
        self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        self._set_status(f"Designing primers for {len(ssr_list):,} SSRs...", TEXT_SECONDARY)
        self.mw.set_status("Running Primer3...")
        self._worker = PrimerWorker(self.state.genome, ssr_list, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Designing: {done:,} / {total:,}", TEXT_SECONDARY)

    def _on_done(self, res, elapsed):
        primers = res["success"]; n_failed = len(res["failed"])
        self.state.clear_downstream_of_primers()
        self.state.primer_results = primers
        self.state.filtered_primer_results = primers

        # Set filter defaults from data range
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
        self.progress_bar.setVisible(False)
        msg = f"{len(primers):,} primer pairs designed in {elapsed:.1f}s"
        if n_failed: msg += f" — {n_failed} SSRs failed"
        self._set_status(msg, SUCCESS)
        self.mw.set_status(msg); self.mw.on_step_complete(2)
        self._cleanup_worker()
        self._refresh()

    def _on_error(self, msg):
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", ERROR)
        self.mw.set_status("Primer design failed")
        self._cleanup_worker()

    def _cleanup_worker(self):
        """Disconnect signals and release large refs from the worker thread.
        Prevents Qt from holding the genome dict and SSR list alive after the
        thread finishes, which causes idle crashes after ~20 minutes."""
        if self._worker is not None:
            try:
                self._worker.progress.disconnect()
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            self._worker.genome   = None
            self._worker.ssr_list = None
            self._worker = None

    def _set_status(self, msg, color):
        self.run_status.setText(msg)
        self.run_status.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # FILTERS
    # ---------------------------------------------------------
    def _on_filter_changed(self):
        self._apply_filters()
        self._populate_table()

    def _apply_filters(self):
        if not self.state.primer_results:
            return
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
        self.results_group.setVisible(True)
        self._apply_filters()
        self._populate_table()

    def _populate_table(self):
        if not self.state.has_primers: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        total = len(self.state.primer_results)
        self.metrics_label.setText(
            f"{total:,} primer pairs designed   |   {len(primers):,} pass quality filters   |   {df['ssr_id'].nunique():,} SSRs covered"
        )
        COLS = {
            "ssr_id":"SSR ID","pair_rank":"Pair rank","contig":"Contig","motif":"Motif",
            "left_primer":"Forward primer","right_primer":"Reverse primer",
            "product_size":"Product size (bp)",
            "left_tm":"Forward Tm (°C)","right_tm":"Reverse Tm (°C)",
            "left_gc":"Forward GC (%)","right_gc":"Reverse GC (%)",
            "left_3end_dg":"Forward 3' stability","right_3end_dg":"Reverse 3' stability",
        }
        display_cols = [c for c in COLS if c in df.columns]

        # Cap display rows to avoid freezing the main thread on large datasets.
        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df
        if truncated:
            self.metrics_label.setText(
                self.metrics_label.text() +
                f"   |   showing first {TABLE_DISPLAY_LIMIT:,} rows — download CSV for full results"
            )

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df)); self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([COLS[c] for c in display_cols])
        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                text = f"{val:.2f}" if isinstance(val, float) else ("" if val is None else str(val))
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(text))
        self.table.resizeColumnsToContents()
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
            "ssr_id":"SSR ID","pair_rank":"Pair rank","contig":"Contig","motif":"Motif",
            "repeat_count":"Repeat count","left_primer":"Forward primer","right_primer":"Reverse primer",
            "product_size":"Product size (bp)","left_tm":"Forward Tm (°C)","right_tm":"Reverse Tm (°C)",
            "left_gc":"Forward GC (%)","right_gc":"Reverse GC (%)",
            "left_3end_dg":"Forward 3' stability (kcal/mol)","right_3end_dg":"Reverse 3' stability (kcal/mol)",
        }).to_csv(path, index=False, encoding="utf-8-sig")
        self.mw.set_status(f"Saved to {path}")

    def _download_fasta(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save BLAST FASTA", "primers_blast_ready.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
        if _root not in sys.path: sys.path.insert(0, _root)
        from core.primer_design import primers_to_blast_fasta
        primers = self.state.filtered_primer_results or self.state.primer_results
        with open(path, "w") as f:
            f.write(primers_to_blast_fasta(primers))
        self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()