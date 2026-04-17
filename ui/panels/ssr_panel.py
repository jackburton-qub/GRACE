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

TABLE_DISPLAY_LIMIT = 10_000


# ---------------------------------------------------------------------------
# SSR detection worker
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# GFF annotation worker
# ---------------------------------------------------------------------------

class GFFAnnotationWorker(QThread):
    # annotated ssrs, gff_index, chrom_names dict
    finished = pyqtSignal(list, object, dict)
    error    = pyqtSignal(str)

    def __init__(self, ssrs, gff_path, gff_features):
        super().__init__()
        self.ssrs         = ssrs
        self.gff_path     = gff_path
        self.gff_features = gff_features

    def run(self):
        try:
            from core.gff_parser import build_gff_index, annotate_ssrs
            if self.gff_features is None:
                index = build_gff_index(self.gff_path)
            else:
                index = self.gff_features
            annotated = annotate_ssrs(self.ssrs, index)
            self.finished.emit(annotated, index, dict(index.chrom_names))
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class SSRPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker     = None
        self._gff_worker = None
        self._last_rendered_version = -1
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
        lbl_motif.setToolTip(
            "Select which repeat unit sizes to search for.\n"
            "Di- (2bp) and tri- (3bp) nucleotide repeats are most informative\n"
            "for population genetics and marker development."
        )
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
        lbl_contig.setToolTip(
            "Skip contigs shorter than this length.\n"
            "Contigs under ~300bp cannot yield usable SSR markers."
        )
        sg.addWidget(lbl_contig, 1, 0)
        self.min_contig = QSpinBox()
        self.min_contig.setRange(0, 100000)
        self.min_contig.setValue(300)
        self.min_contig.setSingleStep(50)
        sg.addWidget(self.min_contig, 1, 1)

        lbl_std = QLabel("Motif standardisation")
        lbl_std.setToolTip(
            "Controls how equivalent motifs are grouped.\n"
            "Level 4 groups all rotations and strand orientations."
        )
        sg.addWidget(lbl_std, 1, 2)
        self.std_level = QComboBox()
        for lvl, desc in [
            (0, "0 — None"),
            (1, "1 — Smallest rotation"),
            (2, "2 — + reverse complement"),
            (3, "3 — Extended"),
            (4, "4 — Full canonical (recommended)"),
        ]:
            self.std_level.addItem(desc, lvl)
        self.std_level.setCurrentIndex(4)
        sg.addWidget(self.std_level, 1, 3)

        self.excl_homo = QCheckBox("Exclude homopolymers")
        self.excl_homo.setChecked(True)
        self.excl_homo.setToolTip(
            "Exclude single-nucleotide runs (e.g. AAAAAAA).\n"
            "Homopolymers are uninformative as genetic markers."
        )
        self.search_rev = QCheckBox("Search reverse complement strand")
        self.search_rev.setChecked(False)
        sg.addWidget(self.excl_homo,  2, 0, 1, 2)
        sg.addWidget(self.search_rev, 2, 2, 1, 2)

        rep_lbl = QLabel("Min repeat counts")
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
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        self.status_label = QLabel("")
        self.status_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        run_row.addWidget(self.run_btn)
        run_row.addWidget(self.cancel_btn)
        run_row.addWidget(self.status_label)
        run_row.addStretch()
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
        sel_row.addWidget(self.sel_all_btn)
        sel_row.addWidget(self.desel_all_btn)
        sel_row.addWidget(self.sel_count_label)
        sel_row.addStretch()
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

    # ---------------------------------------------------------
    # RUN
    # ---------------------------------------------------------
    def _run(self):
        if not self.state.has_genome:
            self._set_status("Load a genome first", WARNING); return
        motif_lengths = tuple(k for k, cb in self.motif_checks.items() if cb.isChecked())
        if not motif_lengths:
            self._set_status("Select at least one motif length", WARNING); return
        params = {
            "motif_lengths":        motif_lengths,
            "min_repeats":          {k: self.min_repeats[k].value() for k in motif_lengths},
            "exclude_homopolymers": self.excl_homo.isChecked(),
            "search_reverse":       self.search_rev.isChecked(),
            "std_level":            self.std_level.currentData(),
            "min_contig_len":       self.min_contig.value(),
        }
        self.run_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        self._set_status("Scanning...", TEXT_SECONDARY)
        self.mw.set_status("Running SSR detection...")
        self._worker = SSRWorker(self.state.genome, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _cancel(self):
        if self._worker is not None:
            self._worker.terminate()
            self._worker.wait()
            self._cleanup_worker()
        if self._gff_worker is not None:
            self._gff_worker.terminate()
            self._gff_worker.wait()
            self._gff_worker = None
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status("Cancelled", WARNING)
        self.mw.set_status("SSR detection cancelled")

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Scanning contig {done:,} of {total:,}...", TEXT_SECONDARY)

    def _on_done(self, ssrs, elapsed):
        self.state.clear_downstream_of_ssrs()
        self.state.ssrs             = ssrs
        self.state.ssr_detection_time = elapsed
        self.state.selected_ssrs    = None
        self.state.ssr_version     += 1
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status(f"Found {len(ssrs):,} SSRs in {elapsed:.1f}s", SUCCESS)
        self.mw.set_status(f"SSR detection complete — {len(ssrs):,} SSRs found")
        self.mw.on_step_complete(1)
        self._cleanup_worker()
        self._refresh()
        if self.state.has_gff and ssrs:
            self._run_gff_annotation(ssrs)

    def _run_gff_annotation(self, ssrs):
        self._set_status(
            f"Found {len(ssrs):,} SSRs — annotating with GFF features...",
            TEXT_SECONDARY
        )
        self.mw.set_status("Annotating SSRs with genomic features...")
        self._gff_worker = GFFAnnotationWorker(
            ssrs, self.state.gff_path, self.state.gff_features
        )
        self._gff_worker.finished.connect(self._on_gff_done)
        self._gff_worker.error.connect(self._on_gff_error)
        self._gff_worker.start()

    def _on_gff_done(self, annotated_ssrs, gff_index, chrom_names):
        self.state.ssrs         = annotated_ssrs
        self.state.gff_features = gff_index
        self.state.chrom_names  = chrom_names if chrom_names else {}
        self.state.ssr_version += 1
        n = len(annotated_ssrs)
        chrom_note = f" — {len(chrom_names):,} chromosome names mapped" if chrom_names else ""
        self._set_status(f"Found {n:,} SSRs — genomic features annotated ✓{chrom_note}", SUCCESS)
        self.mw.set_status(f"SSR annotation complete — {n:,} SSRs annotated")
        if self._gff_worker:
            try:
                self._gff_worker.finished.disconnect()
                self._gff_worker.error.disconnect()
            except Exception:
                pass
            self._gff_worker = None
        self._refresh()

    def _on_gff_error(self, msg):
        self._set_status(f"SSRs found but GFF annotation failed: {msg}", WARNING)
        self.mw.set_status("GFF annotation failed — SSRs available without feature labels")
        if self._gff_worker:
            try:
                self._gff_worker.finished.disconnect()
                self._gff_worker.error.disconnect()
            except Exception:
                pass
            self._gff_worker = None

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", ERROR)
        self.mw.set_status("SSR detection failed")
        self._cleanup_worker()

    def _cleanup_worker(self):
        if self._worker is not None:
            try:
                self._worker.progress.disconnect()
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            self._worker.genome = None
            self._worker.params = None
            self._worker = None

    def _set_status(self, msg, color):
        self.status_label.setText(msg)
        self.status_label.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # REFRESH / TABLE
    # ---------------------------------------------------------
    def _refresh(self):
        if not self.state.has_ssrs:
            self.results_group.setVisible(False); return
        if self.state.ssr_version == self._last_rendered_version:
            return
        self._last_rendered_version = self.state.ssr_version
        self.results_group.setVisible(True)

        import pandas as pd
        ssrs = self.state.ssrs
        df   = pd.DataFrame(ssrs)
        t    = self.state.ssr_detection_time
        truncated  = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df

        metrics = (
            f"{len(ssrs):,} SSRs   |   {df['motif'].nunique():,} unique motifs   |   "
            f"{df['contig'].nunique():,} contigs" + (f"   |   {t:.1f}s" if t else "")
        )
        if truncated:
            metrics += f"   |   showing first {TABLE_DISPLAY_LIMIT:,} rows — download CSV for full results"
        self.metrics_label.setText(metrics)

        has_annotation = "genomic_feature" in df.columns
        cols = ["ssr_id", "contig", "start", "end", "motif", "canonical_motif", "repeat_count"]
        if has_annotation:
            cols.append("genomic_feature")

        col_labels = {
            "ssr_id":          "SSR ID",
            "contig":          "Contig",
            "start":           "Start (bp)",
            "end":             "End (bp)",
            "motif":           "Motif",
            "canonical_motif": "Canonical motif",
            "repeat_count":    "Repeat count",
            "genomic_feature": "Genomic feature",
        }
        display_cols = [c for c in cols if c in df.columns]

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df))
        self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([col_labels.get(c, c) for c in display_cols])

        # Column resize modes — contig gets fixed interactive width,
        # others resize to content, last column stretches
        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        for col_idx, col in enumerate(display_cols):
            if col == "contig":
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(col_idx, 160)
            elif col in ("left_primer", "right_primer"):
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(col_idx, 200)
            else:
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.ResizeToContents)

        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                val  = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                item = QTableWidgetItem(str(val))
                if col_idx == 0:
                    item.setData(
                        Qt.ItemDataRole.UserRole,
                        int(display_df.iat[row_idx, display_df.columns.get_loc("ssr_id")])
                    )
                self.table.setItem(row_idx, col_idx, item)

        self.table.setSortingEnabled(True)

    def _on_selection_changed(self):
        rows = set(i.row() for i in self.table.selectedItems())
        if not rows or not self.state.has_ssrs:
            self.sel_count_label.setText("")
            self.state.selected_ssrs = None
            return
        ssr_id_set = set()
        for row in rows:
            item = self.table.item(row, 0)
            if item is not None:
                uid = item.data(Qt.ItemDataRole.UserRole)
                if uid is not None:
                    ssr_id_set.add(uid)
        if ssr_id_set:
            ssr_by_id = {s["ssr_id"]: s for s in self.state.ssrs}
            self.state.selected_ssrs = [
                ssr_by_id[i] for i in sorted(ssr_id_set) if i in ssr_by_id
            ]
        else:
            self.state.selected_ssrs = None
        count = len(self.state.selected_ssrs) if self.state.selected_ssrs else 0
        self.sel_count_label.setText(f"{count:,} SSRs selected" if count else "")

    def _download_csv(self):
        if not self.state.has_ssrs: return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save SSR results", "ssr_detection_results.csv", "CSV files (*.csv)"
        )
        if path:
            import pandas as pd
            pd.DataFrame(self.state.ssrs).rename(columns={
                "ssr_id":          "SSR ID",
                "contig":          "Contig",
                "start":           "Start (bp)",
                "end":             "End (bp)",
                "motif":           "Motif",
                "canonical_motif": "Canonical motif",
                "repeat_count":    "Repeat count",
                "genomic_feature": "Genomic feature",
            }).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()