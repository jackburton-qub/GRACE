"""
results_panel.py — BLAST Results
Pure display. Single scroll area. Reads state only.
"""
import os, sys

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QScrollArea, QFileDialog,
    QCheckBox,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)

# Cap table rows to prevent main-thread freeze on large datasets.
TABLE_DISPLAY_LIMIT = 10_000
RAW_TABLE_LIMIT     = 5_000   # raw BLAST hits can be enormous — cap lower



def _clean_df(df):
    """Remove nested/messy columns before CSV export."""
    import pandas as pd
    exclude = {"amplicons", "left_hits", "right_hits", "pass_mode"}
    cols = [c for c in df.columns if c not in exclude]
    out = df[cols].copy()
    # Rename to human-readable
    rename = {
        "ssr_id": "SSR ID", "pair_rank": "Pair rank",
        "specificity_status": "Status", "contig": "Contig",
        "motif": "Motif", "canonical_motif": "Canonical motif",
        "repeat_count": "Repeat count", "start": "Start (bp)", "end": "End (bp)",
        "left_primer": "Forward primer", "right_primer": "Reverse primer",
        "product_size": "Product size (bp)",
        "left_tm": "Forward Tm (C)", "right_tm": "Reverse Tm (C)",
        "left_gc": "Forward GC (%)", "right_gc": "Reverse GC (%)",
        "left_3end_dg": "Forward 3' stability (kcal/mol)",
        "right_3end_dg": "Reverse 3' stability (kcal/mol)",
        "left_gc_clamp": "Forward GC clamp", "right_gc_clamp": "Reverse GC clamp",
    }
    return out.rename(columns={k: v for k, v in rename.items() if k in out.columns})

class ResultsPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
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
        L = QVBoxLayout(content)
        L.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        L.setSpacing(16)

        # Title
        title = QLabel("BLAST Results")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Specificity check results — primers verified against the genome")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title); L.addWidget(sub)

        # Summary
        self.summary_group = QGroupBox("Summary")
        sg = QGridLayout(self.summary_group)
        sg.setSpacing(12)
        self.metric_labels = {}
        for col, key in enumerate(["PASS", "FAIL", "UNKNOWN", "Pass rate"]):
            val_lbl = QLabel("—")
            val_lbl.setFont(QFont(FONT_MONO, 18, QFont.Weight.Bold))
            val_lbl.setStyleSheet(f"color: {ACCENT};")
            key_lbl = QLabel(key)
            key_lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
            key_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            val_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            sg.addWidget(val_lbl, 0, col)
            sg.addWidget(key_lbl, 1, col)
            self.metric_labels[key] = val_lbl

        # Pass/fail bar
        self.bar_label = QLabel("")
        self.bar_label.setFixedHeight(8)
        self.bar_label.setStyleSheet(f"background: {TEXT_SECONDARY}; border-radius: 3px;")
        sg.addWidget(self.bar_label, 2, 0, 1, 4)
        L.addWidget(self.summary_group)

        # Results table
        results_group = QGroupBox("Primer Results")
        rg = QVBoxLayout(results_group)

        # Filter checkbox
        filter_row = QHBoxLayout()
        self.pass_only_cb = QCheckBox("Show PASS primers only")
        self.pass_only_cb.setChecked(False)
        self.pass_only_cb.stateChanged.connect(self._populate_table)
        filter_row.addWidget(self.pass_only_cb)
        filter_row.addStretch()
        rg.addLayout(filter_row)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        rg.addWidget(self.table)

        dl_row = QHBoxLayout()
        self.download_full_btn = QPushButton("Download full results (CSV)")
        self.download_full_btn.clicked.connect(self._download_full)
        self.download_pass_btn = QPushButton("Download PASS primers (CSV)")
        self.download_pass_btn.clicked.connect(self._download_pass)
        dl_row.addWidget(self.download_full_btn)
        dl_row.addWidget(self.download_pass_btn)
        dl_row.addStretch()
        rg.addLayout(dl_row)
        L.addWidget(results_group)

        # Raw BLAST hits
        raw_group = QGroupBox("Raw BLAST Hits")
        raw_layout = QVBoxLayout(raw_group)
        raw_lbl = QLabel("One row per primer-genome alignment.")
        raw_lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        raw_layout.addWidget(raw_lbl)

        self.raw_table = QTableWidget()
        self.raw_table.setAlternatingRowColors(True)
        self.raw_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.raw_table.setSortingEnabled(True)
        self.raw_table.setFixedHeight(300)
        raw_layout.addWidget(self.raw_table)

        self.download_raw_btn = QPushButton("Download raw BLAST hits (CSV)")
        self.download_raw_btn.clicked.connect(self._download_raw)
        raw_layout.addWidget(self.download_raw_btn)
        L.addWidget(raw_group)

        L.addStretch()

        # No-results placeholder
        self.no_results_label = QLabel("No BLAST results yet — run Specificity Check first.")
        self.no_results_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self.no_results_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        L.addWidget(self.no_results_label)

        self._refresh()

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        has = self.state.has_specificity
        self.summary_group.setVisible(has)
        self.no_results_label.setVisible(not has)

        if not has:
            self._clear_tables()
            return

        # Skip re-render if data hasn't changed since last render
        if self.state.blast_version == self._last_rendered_version:
            return
        self._last_rendered_version = self.state.blast_version

        import pandas as pd
        spec_df   = pd.DataFrame(self.state.specificity_results)
        passed    = spec_df[spec_df["specificity_status"] == "PASS"]
        failed    = spec_df[spec_df["specificity_status"] == "FAIL"]
        unknown   = spec_df[spec_df["specificity_status"] == "UNKNOWN"]
        total     = max(len(spec_df), 1)
        pass_rate = 100 * len(passed) / total

        self.metric_labels["PASS"].setText(str(len(passed)))
        self.metric_labels["PASS"].setStyleSheet(f"color: {SUCCESS}; font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["FAIL"].setText(str(len(failed)))
        self.metric_labels["FAIL"].setStyleSheet(f"color: {ERROR}; font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["UNKNOWN"].setText(str(len(unknown)))
        self.metric_labels["UNKNOWN"].setStyleSheet(f"color: {WARNING}; font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["Pass rate"].setText(f"{pass_rate:.1f}%")
        self.metric_labels["Pass rate"].setStyleSheet(f"color: {ACCENT}; font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")

        # Colour bar
        bar_p = int(pass_rate); bar_f = int(100 * len(failed) / total)
        bar_u = max(0, 100 - bar_p - bar_f)
        self.bar_label.setStyleSheet(
            f"background: qlineargradient(x1:0, y1:0, x2:1, y2:0, "
            f"stop:0 {SUCCESS}, stop:{bar_p/100:.3f} {SUCCESS}, "
            f"stop:{bar_p/100:.3f} {ERROR}, stop:{(bar_p+bar_f)/100:.3f} {ERROR}, "
            f"stop:{(bar_p+bar_f)/100:.3f} {WARNING}, stop:1 {WARNING}); "
            f"border-radius: 3px;"
        )

        self._populate_table()
        self._populate_raw_table()

    def _populate_table(self):
        if not self.state.has_specificity:
            return
        import pandas as pd
        from PyQt6.QtGui import QColor
        spec_df = pd.DataFrame(self.state.specificity_results)
        if self.pass_only_cb.isChecked():
            df = spec_df[spec_df["specificity_status"] == "PASS"].copy()
        else:
            df = spec_df.copy()

        COLS = {
            "ssr_id":             "SSR ID",
            "pair_rank":          "Pair rank",
            "specificity_status": "Status",
            "contig":             "Contig",
            "motif":              "Motif",
            "repeat_count":       "Repeat count",
            "left_primer":        "Forward primer",
            "right_primer":       "Reverse primer",
            "product_size":       "Product size (bp)",
            "left_tm":            "Forward Tm (°C)",
            "right_tm":           "Reverse Tm (°C)",
            "left_gc":            "Forward GC (%)",
            "right_gc":           "Reverse GC (%)",
            "left_3end_dg":       "Forward 3' stability",
            "right_3end_dg":      "Reverse 3' stability",
        }
        exclude = {"amplicons", "left_hits", "right_hits", "pass_mode"}
        display_cols = [c for c in COLS if c in df.columns and c not in exclude]

        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df))
        self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([COLS[c] for c in display_cols])

        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                text = f"{val:.2f}" if isinstance(val, float) else ("" if val is None else str(val))
                item = QTableWidgetItem(text)
                if col == "specificity_status":
                    if text == "PASS":
                        item.setForeground(QColor(SUCCESS))
                    elif text == "FAIL":
                        item.setForeground(QColor(ERROR))
                self.table.setItem(row_idx, col_idx, item)

        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)

        if truncated:
            # Update metrics label to note truncation
            current = self.metric_labels["Pass rate"].toolTip()
            self.table.setToolTip(
                f"Showing first {TABLE_DISPLAY_LIMIT:,} rows — download CSV for full results"
            )

    def _populate_raw_table(self):
        raw = self.state.blast_raw_rows
        if raw is None:
            return
        import pandas as pd
        if not isinstance(raw, pd.DataFrame):
            raw = pd.DataFrame(raw)
        if raw.empty:
            return

        COL_RENAME = {
            "ssr_id": "SSR ID", "primer_side": "Direction",
            "qseqid": "Query ID", "sseqid": "Subject contig",
            "pident": "% Identity", "length": "Alignment length",
            "mismatch": "Mismatches", "gapopen": "Gap opens",
            "qstart": "Q start", "qend": "Q end",
            "sstart": "S start", "send": "S end",
            "evalue": "E-value", "bitscore": "Bit score",
            "sstrand": "Strand",
        }
        display_cols = [c for c in COL_RENAME if c in raw.columns]

        truncated = len(raw) > RAW_TABLE_LIMIT
        display_raw = raw.iloc[:RAW_TABLE_LIMIT] if truncated else raw

        self.raw_table.setSortingEnabled(False)
        self.raw_table.setRowCount(len(display_raw))
        self.raw_table.setColumnCount(len(display_cols))
        self.raw_table.setHorizontalHeaderLabels([COL_RENAME[c] for c in display_cols])

        for row_idx in range(len(display_raw)):
            for col_idx, col in enumerate(display_cols):
                val = display_raw.iat[row_idx, display_raw.columns.get_loc(col)]
                text = f"{val:.4f}" if isinstance(val, float) else ("" if val is None else str(val))
                self.raw_table.setItem(row_idx, col_idx, QTableWidgetItem(text))

        self.raw_table.resizeColumnsToContents()
        self.raw_table.setSortingEnabled(True)

        if truncated:
            self.raw_table.setToolTip(
                f"Showing first {RAW_TABLE_LIMIT:,} rows — download CSV for full results"
            )

    def _clear_tables(self):
        self.table.setRowCount(0)
        self.raw_table.setRowCount(0)

    # ---------------------------------------------------------
    # DOWNLOADS
    # ---------------------------------------------------------
    def _get_display_df(self):
        if not self.state.has_specificity:
            return None
        import pandas as pd
        return pd.DataFrame(self.state.specificity_results)

    def _clean_df(self, df):
        """Remove internal columns and rename for clean CSV output."""
        exclude = {"amplicons", "left_hits", "right_hits", "pass_mode",
                   "left_gc_clamp", "right_gc_clamp"}
        rename = {
            "ssr_id":             "SSR ID",
            "pair_rank":          "Pair rank",
            "specificity_status": "Status",
            "contig":             "Contig",
            "start":              "Start (bp)",
            "end":                "End (bp)",
            "motif":              "Motif",
            "canonical_motif":    "Canonical motif",
            "repeat_count":       "Repeat count",
            "left_primer":        "Forward primer",
            "right_primer":       "Reverse primer",
            "product_size":       "Product size (bp)",
            "left_tm":            "Forward Tm (°C)",
            "right_tm":           "Reverse Tm (°C)",
            "left_gc":            "Forward GC (%)",
            "right_gc":           "Reverse GC (%)",
            "left_3end_dg":       "Forward 3' stability (kcal/mol)",
            "right_3end_dg":      "Reverse 3' stability (kcal/mol)",
        }
        cols = [col for col in df.columns if col not in exclude]
        return df[cols].rename(columns=rename)

    def _download_full(self):
        df = self._get_display_df()
        if df is None: return
        path, _ = QFileDialog.getSaveFileName(self, "Save full results", "specificity_results.csv", "CSV files (*.csv)")
        if path:
            self._clean_df(df).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def _download_pass(self):
        df = self._get_display_df()
        if df is None: return
        path, _ = QFileDialog.getSaveFileName(self, "Save PASS primers", "pass_primers.csv", "CSV files (*.csv)")
        if path:
            self._clean_df(df[df["specificity_status"] == "PASS"]).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def _download_raw(self):
        raw = self.state.blast_raw_rows
        if raw is None: return
        import pandas as pd
        if not isinstance(raw, pd.DataFrame): raw = pd.DataFrame(raw)
        path, _ = QFileDialog.getSaveFileName(self, "Save raw BLAST hits", "blast_raw_hits.csv", "CSV files (*.csv)")
        if path:
            raw.to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()