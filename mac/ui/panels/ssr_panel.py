"""
ssr_panel.py — SSR Detection
Single scroll area. Table has fixed height. No layout tricks.
"""
import os, sys, time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QCheckBox, QComboBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)


class SSRWorker(QThread):
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(list, float)
    error    = pyqtSignal(str)

    def __init__(self, genome, params):
        super().__init__()
        self.genome = genome
        self.params = params

    def run(self):
        try:
            _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            if _root not in sys.path:
                sys.path.insert(0, _root)
            from core.ssr_detection import find_ssrs
            start = time.time()
            ssrs = find_ssrs(
                self.genome,
                motif_lengths=self.params["motif_lengths"],
                min_repeats=self.params["min_repeats"],
                exclude_homopolymers=self.params["exclude_homopolymers"],
                search_reverse=self.params["search_reverse"],
                motif_standardisation_level=self.params["std_level"],
                min_contig_len=self.params["min_contig_len"],
                progress_callback=lambda d, t: self.progress.emit(d, t),
            )
            self.finished.emit(ssrs, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


class SSRPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker = None
        self._build_ui()

    def _build_ui(self):
        # One scroll area, everything inside
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        outer.addWidget(scroll)

        content = QWidget()
        scroll.setWidget(content)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        layout.setSpacing(16)

        # Title
        title = QLabel("SSR Detection")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Scan your genome for simple sequence repeats")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        layout.addWidget(title)
        layout.addWidget(sub)

        # Settings
        settings_group = QGroupBox("Detection Settings")
        sg = QGridLayout(settings_group)
        sg.setSpacing(10)

        lbl_motif = QLabel("Motif lengths")
        lbl_motif.setToolTip("Select which repeat unit sizes to search for.\nDi- (2bp) and tri- (3bp) nucleotide repeats are most informative\nfor population genetics and marker development.")
        sg.addWidget(lbl_motif, 0, 0)
        motif_row = QHBoxLayout()
        self.motif_checks = {}
        for k in [2, 3, 4, 5, 6]:
            cb = QCheckBox(f"{k}-mer")
            cb.setChecked(k in [2, 3])
            cb.stateChanged.connect(self._update_repeat_visibility)
            self.motif_checks[k] = cb
            motif_row.addWidget(cb)
        motif_row.addStretch()
        motif_w = QWidget(); motif_w.setLayout(motif_row)
        sg.addWidget(motif_w, 0, 1, 1, 3)

        lbl_contig = QLabel("Min contig length (bp)")
        lbl_contig.setToolTip("Skip contigs shorter than this length.\nContigs under ~300bp cannot yield usable SSR markers.\nIncreasing this speeds up scanning on fragmented assemblies.")
        sg.addWidget(lbl_contig, 1, 0)
        self.min_contig = QSpinBox()
        self.min_contig.setRange(0, 100000); self.min_contig.setValue(300); self.min_contig.setSingleStep(50)
        sg.addWidget(self.min_contig, 1, 1)

        lbl_std = QLabel("Motif standardisation")
        lbl_std.setToolTip("Controls how equivalent motifs are grouped.\nLevel 4 groups all rotations and strand orientations,\ne.g. AC, CA, TG, GT are all reported as AC.")
        sg.addWidget(lbl_std, 1, 2)
        self.std_level = QComboBox()
        for lvl, desc in [(0,"0 — None"),(1,"1 — Smallest rotation"),(2,"2 — + reverse complement"),(3,"3 — Extended"),(4,"4 — Full canonical (recommended)")]:
            self.std_level.addItem(desc, lvl)
        self.std_level.setCurrentIndex(4)
        sg.addWidget(self.std_level, 1, 3)

        self.excl_homo = QCheckBox("Exclude homopolymers")
        self.excl_homo.setChecked(True)
        self.excl_homo.setToolTip("Exclude single-nucleotide runs (e.g. AAAAAAA).\nHomopolymers are uninformative as genetic markers. Recommended: on.")
        self.search_rev = QCheckBox("Search reverse complement strand")
        self.search_rev.setChecked(False)
        self.search_rev.setToolTip("Also scan the reverse complement of each contig.\nMay produce duplicate entries for palindromic SSRs.")
        sg.addWidget(self.excl_homo,  2, 0, 1, 2)
        sg.addWidget(self.search_rev, 2, 2, 1, 2)

        rep_lbl = QLabel("Min repeat counts")
        rep_lbl.setToolTip("The motif must occur this many times consecutively to be reported.\nRecommended: 6 for di-nucleotide, 5 for tri- and longer.")
        sg.addWidget(rep_lbl, 3, 0)
        self.min_repeats = {}
        self.min_repeat_labels = {}
        rep_row = QHBoxLayout(); rep_row.setSpacing(8)
        for k in [2, 3, 4, 5, 6]:
            lbl = QLabel(f"{k}-mer:")
            lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
            sp = QSpinBox(); sp.setRange(2, 20); sp.setValue(6 if k == 2 else 5)
            sp.setMinimumWidth(70); sp.setFixedHeight(28)
            self.min_repeats[k] = sp
            self.min_repeat_labels[k] = lbl
            rep_row.addWidget(lbl); rep_row.addWidget(sp)
        rep_row.addStretch()
        rep_w = QWidget(); rep_w.setLayout(rep_row)
        sg.addWidget(rep_w, 3, 1, 1, 3)
        layout.addWidget(settings_group)

        # Run
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Run SSR Detection")
        self.run_btn.setObjectName("primary")
        self.run_btn.clicked.connect(self._run)
        self.status_label = QLabel("")
        self.status_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        run_row.addWidget(self.run_btn); run_row.addWidget(self.status_label); run_row.addStretch()
        layout.addLayout(run_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False); self.progress_bar.setFixedHeight(6)
        layout.addWidget(self.progress_bar)

        # Results
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        res_layout = QVBoxLayout(self.results_group)

        self.metrics_label = QLabel("")
        self.metrics_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        res_layout.addWidget(self.metrics_label)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        res_layout.addWidget(self.table)

        sel_row = QHBoxLayout()
        self.sel_all_btn   = QPushButton("Select all")
        self.desel_all_btn = QPushButton("Deselect all")
        self.sel_all_btn.clicked.connect(lambda: self.table.selectAll())
        self.desel_all_btn.clicked.connect(lambda: self.table.clearSelection())
        self.sel_count_label = QLabel("")
        self.sel_count_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        self.download_btn = QPushButton("Download SSR results (CSV)")
        self.download_btn.clicked.connect(self._download_csv)
        sel_row.addWidget(self.sel_all_btn); sel_row.addWidget(self.desel_all_btn)
        sel_row.addWidget(self.sel_count_label); sel_row.addStretch()
        sel_row.addWidget(self.download_btn)
        res_layout.addLayout(sel_row)

        layout.addWidget(self.results_group)
        layout.addStretch()

        self._update_repeat_visibility()
        self._refresh()

    def _update_repeat_visibility(self):
        for k in [2, 3, 4, 5, 6]:
            v = self.motif_checks[k].isChecked()
            self.min_repeat_labels[k].setVisible(v)
            self.min_repeats[k].setVisible(v)

    def _run(self):
        if not self.state.has_genome:
            self._set_status("Load a genome first", WARNING); return
        motif_lengths = tuple(k for k, cb in self.motif_checks.items() if cb.isChecked())
        if not motif_lengths:
            self._set_status("Select at least one motif length", WARNING); return
        params = {
            "motif_lengths": motif_lengths,
            "min_repeats": {k: self.min_repeats[k].value() for k in motif_lengths},
            "exclude_homopolymers": self.excl_homo.isChecked(),
            "search_reverse": self.search_rev.isChecked(),
            "std_level": self.std_level.currentData(),
            "min_contig_len": self.min_contig.value(),
        }
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        self._set_status("Scanning...", TEXT_SECONDARY)
        self.mw.set_status("Running SSR detection...")
        self._worker = SSRWorker(self.state.genome, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Scanning contig {done:,} of {total:,}...", TEXT_SECONDARY)

    def _on_done(self, ssrs, elapsed):
        self.state.clear_downstream_of_ssrs()
        self.state.ssrs = ssrs; self.state.ssr_detection_time = elapsed; self.state.selected_ssrs = None
        self.run_btn.setEnabled(True); self.progress_bar.setVisible(False)
        self._set_status(f"Found {len(ssrs):,} SSRs in {elapsed:.1f}s", SUCCESS)
        self.mw.set_status(f"SSR detection complete — {len(ssrs):,} SSRs found")
        self.mw.on_step_complete(1)
        self._refresh()

    def _on_error(self, msg):
        self.run_btn.setEnabled(True); self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", ERROR)
        self.mw.set_status("SSR detection failed")

    def _set_status(self, msg, color):
        self.status_label.setText(msg)
        self.status_label.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    def _refresh(self):
        if not self.state.has_ssrs:
            self.results_group.setVisible(False); return
        self.results_group.setVisible(True)
        import pandas as pd
        ssrs = self.state.ssrs; df = pd.DataFrame(ssrs)
        t = self.state.ssr_detection_time
        self.metrics_label.setText(
            f"{len(ssrs):,} SSRs   |   {df['motif'].nunique():,} unique motifs   |   "
            f"{df['contig'].nunique():,} contigs" + (f"   |   {t:.1f}s" if t else "")
        )
        cols = ["ssr_id","contig","start","end","motif","canonical_motif","repeat_count"]
        col_labels = {"ssr_id":"SSR ID","contig":"Contig","start":"Start (bp)","end":"End (bp)",
                      "motif":"Motif","canonical_motif":"Canonical motif","repeat_count":"Repeat count"}
        display_cols = [c for c in cols if c in df.columns]
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(df)); self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([col_labels.get(c,c) for c in display_cols])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        for row_idx in range(len(df)):
            for col_idx, col in enumerate(display_cols):
                val = df.iat[row_idx, df.columns.get_loc(col)]
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(str(val)))
        self.table.setSortingEnabled(True)

    def _on_selection_changed(self):
        rows = set(i.row() for i in self.table.selectedItems())
        self.sel_count_label.setText(f"{len(rows):,} SSRs selected" if rows else "")
        self.state.selected_ssrs = [self.state.ssrs[r] for r in sorted(rows)] if rows and self.state.has_ssrs else None

    def _download_csv(self):
        if not self.state.has_ssrs: return
        path, _ = QFileDialog.getSaveFileName(self, "Save SSR results", "ssr_detection_results.csv", "CSV files (*.csv)")
        if path:
            import pandas as pd
            pd.DataFrame(self.state.ssrs).rename(columns={
                "ssr_id":"SSR ID","contig":"Contig","start":"Start (bp)","end":"End (bp)",
                "motif":"Motif","canonical_motif":"Canonical motif","repeat_count":"Repeat count"
            }).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()