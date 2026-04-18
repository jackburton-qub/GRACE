"""
gbs_re_panel.py — GBS-RE Marker Discovery
Identify SSRs within restriction fragments for traditional GBS.
Includes LD filter, chromosome view, and multi-enzyme comparison.
"""

import os
import math
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QFileDialog, QSpinBox,
    QComboBox, QCheckBox, QApplication, QTextEdit, QProgressBar,
    QSplitter, QTabWidget, QListWidget, QListWidgetItem, QSizePolicy
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QRectF, QPointF
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)

def _lbl(text, tip=None):
    w = QLabel(text)
    w.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
    if tip:
        w.setToolTip(tip)
    return w


class ChromosomeView(QWidget):
    """Visualize marker positions across contigs."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self._markers = []          # list of (contig, start, end) for all markers
        self._filtered_markers = [] # after LD filter (if applied)
        self._show_filtered = False
        self._contig_lengths = {}
        self.setMinimumHeight(200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, all_markers, filtered_markers=None, contig_lengths=None):
        self._markers = all_markers
        self._filtered_markers = filtered_markers if filtered_markers is not None else []
        self._contig_lengths = contig_lengths or {}
        self.update()

    def set_show_filtered(self, show):
        self._show_filtered = show
        self.update()

    def paintEvent(self, event):
        if not self._markers and not self._filtered_markers:
            return

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        margin_left, margin_right, margin_top, margin_bottom = 100, 20, 20, 30

        # Group markers by contig
        by_contig = {}
        markers_to_show = self._filtered_markers if self._show_filtered else self._markers
        for m in markers_to_show:
            contig = m.get("contig", "")
            if contig not in by_contig:
                by_contig[contig] = []
            by_contig[contig].append(m)

        if not by_contig:
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No markers to display")
            return

        # Get contig lengths if available, else estimate from max position
        contigs = list(by_contig.keys())
        contig_lengths = {}
        for contig in contigs:
            if contig in self._contig_lengths:
                contig_lengths[contig] = self._contig_lengths[contig]
            else:
                max_pos = max(m.get("start", 0) for m in by_contig[contig])
                contig_lengths[contig] = max_pos + 1000

        # Sort contigs by length (largest first) or alphabetically
        contigs = sorted(contigs, key=lambda c: contig_lengths.get(c, 0), reverse=True)
        # Take top 8 contigs to keep view manageable
        contigs = contigs[:8]

        n_contigs = len(contigs)
        row_height = (h - margin_top - margin_bottom) / max(n_contigs, 1)
        bar_height = row_height * 0.6

        painter.setPen(QPen(QColor(TEXT_SECONDARY), 1))
        painter.setFont(QFont(FONT_MONO, 8))

        for i, contig in enumerate(contigs):
            y = margin_top + i * row_height + (row_height - bar_height) / 2
            contig_len = contig_lengths[contig]

            # Draw contig label
            painter.drawText(5, int(y + bar_height/2 + 4), contig[:15])

            # Draw chromosome bar
            bar_x = margin_left
            bar_w = w - margin_left - margin_right
            painter.setBrush(QBrush(QColor(BG_MID)))
            painter.setPen(QPen(QColor(BORDER), 1))
            painter.drawRect(bar_x, int(y), int(bar_w), int(bar_height))

            # Draw markers
            markers = by_contig[contig]
            painter.setBrush(QBrush(QColor(ACCENT)))
            painter.setPen(Qt.PenStyle.NoPen)
            for m in markers:
                pos = m.get("start", 0)
                if contig_len > 0:
                    rel_x = bar_x + (pos / contig_len) * bar_w
                    painter.drawEllipse(QPointF(rel_x, y + bar_height/2), 3, 3)

            # Draw scale tick
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(bar_x, int(y + bar_height + 12), "0")
            painter.drawText(bar_x + bar_w - 30, int(y + bar_height + 12), f"{contig_len/1000:.1f} kb")


class GBSREWorker(QThread):
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, ssrs, genome, params):
        super().__init__()
        self.ssrs = ssrs
        self.genome = genome
        self.params = params

    def run(self):
        try:
            from core.gbs_re_finder import find_ssrs_in_fragments
            result = find_ssrs_in_fragments(
                self.ssrs,
                self.genome,
                enzyme1=self.params["enzyme1"],
                enzyme2=self.params.get("enzyme2"),
                min_fragment_size=self.params["min_frag"],
                max_fragment_size=self.params["max_frag"],
                min_distance_from_end=self.params["min_end_dist"],
                progress_callback=lambda done, total: self.progress.emit(done, total),
            )
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class MultiEnzymeWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)

    def __init__(self, ssrs, genome, enzyme_combos, min_frag, max_frag, min_end_dist):
        super().__init__()
        self.ssrs = ssrs
        self.genome = genome
        self.enzyme_combos = enzyme_combos
        self.min_frag = min_frag
        self.max_frag = max_frag
        self.min_end_dist = min_end_dist

    def run(self):
        try:
            from core.gbs_re_finder import find_ssrs_in_fragments
            results = []
            for combo in self.enzyme_combos:
                self.progress.emit(f"Testing {combo['display']}...")
                enzyme1 = combo["enzyme1"]
                enzyme2 = combo.get("enzyme2")
                res = find_ssrs_in_fragments(
                    self.ssrs,
                    self.genome,
                    enzyme1=enzyme1,
                    enzyme2=enzyme2,
                    min_fragment_size=self.min_frag,
                    max_fragment_size=self.max_frag,
                    min_distance_from_end=self.min_end_dist,
                )
                # Compute median fragment size
                frag_sizes = res.get("passing_fragment_sizes", [])
                median_size = sorted(frag_sizes)[len(frag_sizes)//2] if frag_sizes else 0
                results.append({
                    "display": combo["display"],
                    "enzyme1": enzyme1,
                    "enzyme2": enzyme2,
                    "n_qualified": res["n_qualified"],
                    "passing_fragments": res["passing_fragments"],
                    "median_fragment": median_size,
                })
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))


class BarcodeWorker(QThread):
    finished = pyqtSignal(list)
    error = pyqtSignal(str)

    def __init__(self, n_barcodes, length, min_dist):
        super().__init__()
        self.n_barcodes = n_barcodes
        self.length = length
        self.min_dist = min_dist

    def run(self):
        try:
            from core.barcode_designer import generate_barcodes
            barcodes = generate_barcodes(self.n_barcodes, self.length, self.min_dist)
            self.finished.emit(barcodes)
        except Exception as e:
            self.error.emit(str(e))


class GBSREPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw = main_window
        self._worker = None
        self._barcode_worker = None
        self._multi_worker = None
        self._qualified_ssrs = []
        self._all_markers = []      # for chromosome view
        self._contig_lengths = {}   # from genome
        self._build_ui()
        self._compute_contig_lengths()

    def _compute_contig_lengths(self):
        if self.state.genome:
            self._contig_lengths = {c: len(s) for c, s in self.state.genome.items()}

    def _toggle_ld_filter(self, checked):
        self._ld_distance_spin.setEnabled(checked)
        if hasattr(self, '_chromosome_view'):
            self._chromosome_view.set_show_filtered(checked and self._ld_filter_cb.isChecked())
            self._chromosome_view.update()

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

        title = QLabel("GBS‑RE Marker Discovery")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Identify SSRs within restriction fragments for traditional GBS")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # Enzyme selection
        enzyme_group = QGroupBox("Enzyme Selection")
        eg = QVBoxLayout(enzyme_group)
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Enzyme 1:"))
        self._enzyme1_combo = QComboBox()
        self._enzyme1_combo.addItems(["PstI", "MspI", "ApeKI", "MseI", "SbfI", "EcoRI", "HindIII", "TaqI", "NlaIII"])
        self._enzyme1_combo.setCurrentText("PstI")
        row1.addWidget(self._enzyme1_combo)

        row1.addWidget(QLabel("Enzyme 2 (optional):"))
        self._enzyme2_combo = QComboBox()
        self._enzyme2_combo.addItems(["None", "MspI", "ApeKI", "MseI", "TaqI", "NlaIII"])
        row1.addWidget(self._enzyme2_combo)
        eg.addLayout(row1)

        L.addWidget(enzyme_group)

        # Fragment size filter
        frag_group = QGroupBox("Fragment Size Filter")
        fg = QVBoxLayout(frag_group)
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Min fragment (bp):"))
        self._min_frag_spin = QSpinBox()
        self._min_frag_spin.setRange(50, 1000)
        self._min_frag_spin.setValue(100)
        row2.addWidget(self._min_frag_spin)

        row2.addWidget(QLabel("Max fragment (bp):"))
        self._max_frag_spin = QSpinBox()
        self._max_frag_spin.setRange(100, 2000)
        self._max_frag_spin.setValue(400)
        row2.addWidget(self._max_frag_spin)

        row2.addWidget(QLabel("Min distance from ends (bp):"))
        self._min_end_dist_spin = QSpinBox()
        self._min_end_dist_spin.setRange(0, 100)
        self._min_end_dist_spin.setValue(20)
        row2.addWidget(self._min_end_dist_spin)
        fg.addLayout(row2)

        # LD filter
        ld_row = QHBoxLayout()
        self._ld_filter_cb = QCheckBox("Apply LD filter (thin markers by distance)")
        self._ld_filter_cb.setChecked(False)
        self._ld_filter_cb.toggled.connect(self._toggle_ld_filter)
        ld_row.addWidget(self._ld_filter_cb)

        ld_row.addWidget(QLabel("Min distance (bp):"))
        self._ld_distance_spin = QSpinBox()
        self._ld_distance_spin.setRange(1000, 1000000)
        self._ld_distance_spin.setValue(10000)
        self._ld_distance_spin.setEnabled(False)
        ld_row.addWidget(self._ld_distance_spin)
        ld_row.addStretch()
        fg.addLayout(ld_row)

        L.addWidget(frag_group)

        # Action buttons
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Discover GBS Markers")
        self._run_btn.setObjectName("primary")
        self._run_btn.clicked.connect(self._run_discovery)
        btn_row.addWidget(self._run_btn)

        self._compare_btn = QPushButton("Compare Enzymes")
        self._compare_btn.clicked.connect(self._show_comparison_dialog)
        btn_row.addWidget(self._compare_btn)
        btn_row.addStretch()
        L.addLayout(btn_row)

        self._progress_bar = QProgressBar()
        self._progress_bar.setVisible(False)
        self._progress_bar.setFixedHeight(6)
        L.addWidget(self._progress_bar)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        L.addWidget(self._status_label)

        # Tabs for results and chromosome view
        self._tab_widget = QTabWidget()
        self._tab_widget.setVisible(False)

        # Results tab
        results_tab = QWidget()
        rt_layout = QVBoxLayout(results_tab)
        self._summary_label = QLabel("")
        rt_layout.addWidget(self._summary_label)

        self._results_table = QTableWidget()
        self._results_table.setAlternatingRowColors(True)
        self._results_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._results_table.setSortingEnabled(True)
        self._results_table.setFixedHeight(300)
        rt_layout.addWidget(self._results_table)

        export_layout = QHBoxLayout()
        self._export_csv_btn = QPushButton("Export Markers (CSV)")
        self._export_csv_btn.clicked.connect(self._export_csv)
        export_layout.addWidget(self._export_csv_btn)

        self._export_bed_btn = QPushButton("Export BED File")
        self._export_bed_btn.clicked.connect(self._export_bed)
        export_layout.addWidget(self._export_bed_btn)
        rt_layout.addLayout(export_layout)

        self._tab_widget.addTab(results_tab, "Markers")

        # Chromosome view tab
        chrom_tab = QWidget()
        ct_layout = QVBoxLayout(chrom_tab)
        self._chromosome_view = ChromosomeView()
        ct_layout.addWidget(self._chromosome_view)
        chrom_controls = QHBoxLayout()
        self._show_filtered_cb = QCheckBox("Show LD-filtered markers")
        self._show_filtered_cb.setChecked(False)
        self._show_filtered_cb.toggled.connect(self._chromosome_view.set_show_filtered)
        chrom_controls.addWidget(self._show_filtered_cb)
        chrom_controls.addStretch()
        ct_layout.addLayout(chrom_controls)
        self._tab_widget.addTab(chrom_tab, "Chromosome View")

        L.addWidget(self._tab_widget)

        # Barcode section
        self._barcode_group = QGroupBox("Sample Barcodes")
        self._barcode_group.setVisible(False)
        bg = QVBoxLayout(self._barcode_group)

        barcode_controls = QHBoxLayout()
        barcode_controls.addWidget(QLabel("Number of samples:"))
        self._n_samples_spin = QSpinBox()
        self._n_samples_spin.setRange(1, 384)
        self._n_samples_spin.setValue(96)
        barcode_controls.addWidget(self._n_samples_spin)

        barcode_controls.addWidget(QLabel("Barcode length:"))
        self._barcode_len_spin = QSpinBox()
        self._barcode_len_spin.setRange(6, 12)
        self._barcode_len_spin.setValue(8)
        barcode_controls.addWidget(self._barcode_len_spin)

        barcode_controls.addWidget(QLabel("Min Hamming distance:"))
        self._min_dist_spin = QSpinBox()
        self._min_dist_spin.setRange(2, 6)
        self._min_dist_spin.setValue(3)
        barcode_controls.addWidget(self._min_dist_spin)

        self._gen_barcode_btn = QPushButton("Generate Barcodes")
        self._gen_barcode_btn.clicked.connect(self._generate_barcodes)
        barcode_controls.addWidget(self._gen_barcode_btn)
        bg.addLayout(barcode_controls)

        self._barcode_text = QTextEdit()
        self._barcode_text.setReadOnly(True)
        self._barcode_text.setFixedHeight(150)
        bg.addWidget(self._barcode_text)

        self._export_barcodes_btn = QPushButton("Export Barcodes (CSV)")
        self._export_barcodes_btn.clicked.connect(self._export_barcodes)
        bg.addWidget(self._export_barcodes_btn)

        L.addWidget(self._barcode_group)
        L.addStretch()

    def _show_comparison_dialog(self):
        from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QCheckBox, QProgressBar
        dlg = QDialog(self)
        dlg.setWindowTitle("Select Enzyme Combinations to Compare")
        dlg.resize(400, 300)
        layout = QVBoxLayout(dlg)

        label = QLabel("Choose enzyme combinations:")
        layout.addWidget(label)

        # Common single enzymes
        single_enzymes = ["PstI", "MspI", "ApeKI", "MseI", "SbfI", "EcoRI", "HindIII", "TaqI", "NlaIII"]
        # Common double digests
        double_combos = [("PstI", "MspI"), ("PstI", "ApeKI"), ("MspI", "ApeKI"), ("SbfI", "MspI")]
        checkboxes = {}
        for enz in single_enzymes:
            cb = QCheckBox(enz)
            cb.setChecked(enz == "PstI")
            checkboxes[("single", enz)] = cb
            layout.addWidget(cb)
        for e1, e2 in double_combos:
            display = f"{e1} + {e2}"
            cb = QCheckBox(display)
            cb.setChecked(False)
            checkboxes[("double", e1, e2)] = cb
            layout.addWidget(cb)

        progress = QProgressBar()
        progress.setVisible(False)
        layout.addWidget(progress)

        btn_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btn_box.accepted.connect(dlg.accept)
        btn_box.rejected.connect(dlg.reject)
        layout.addWidget(btn_box)

        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        # Gather selected combos
        combos = []
        for key, cb in checkboxes.items():
            if cb.isChecked():
                if key[0] == "single":
                    combos.append({"display": key[1], "enzyme1": key[1], "enzyme2": None})
                else:
                    combos.append({"display": f"{key[1]}+{key[2]}", "enzyme1": key[1], "enzyme2": key[2]})

        if not combos:
            return

        # Run comparison
        self._multi_worker = MultiEnzymeWorker(
            self.state.ssrs,
            self.state.genome,
            combos,
            self._min_frag_spin.value(),
            self._max_frag_spin.value(),
            self._min_end_dist_spin.value()
        )
        progress.setVisible(True)
        self._multi_worker.progress.connect(lambda msg: progress.setFormat(msg))
        self._multi_worker.finished.connect(lambda res: self._show_comparison_results(res, dlg))
        self._multi_worker.error.connect(lambda err: self._status_label.setText(f"Error: {err}"))
        self._multi_worker.start()

    def _show_comparison_results(self, results, dialog):
        dialog.close()
        # Display in a new dialog with table
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton
        dlg = QDialog(self)
        dlg.setWindowTitle("Enzyme Comparison Results")
        dlg.resize(700, 400)
        layout = QVBoxLayout(dlg)
        table = QTableWidget()
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Enzyme(s)", "Qualified SSRs", "Passing Fragments", "Median Fragment", "Rank"])
        table.setRowCount(len(results))
        # Sort by qualified SSRs descending
        results = sorted(results, key=lambda x: x["n_qualified"], reverse=True)
        for i, r in enumerate(results):
            table.setItem(i, 0, QTableWidgetItem(r["display"]))
            table.setItem(i, 1, QTableWidgetItem(str(r["n_qualified"])))
            table.setItem(i, 2, QTableWidgetItem(str(r["passing_fragments"])))
            table.setItem(i, 3, QTableWidgetItem(f"{r['median_fragment']} bp"))
            table.setItem(i, 4, QTableWidgetItem(str(i+1)))
        table.resizeColumnsToContents()
        layout.addWidget(table)
        btn = QPushButton("Close")
        btn.clicked.connect(dlg.accept)
        layout.addWidget(btn)
        dlg.exec()

    def _run_discovery(self):
        if not self.state.has_ssrs:
            self._status_label.setText("No SSRs detected. Run SSR detection first.")
            return
        if not self.state.has_genome:
            self._status_label.setText("Genome not loaded.")
            return

        self._run_btn.setEnabled(False)
        self._progress_bar.setVisible(True)
        self._progress_bar.setRange(0, len(self.state.ssrs))
        self._progress_bar.setValue(0)
        self._status_label.setText("Scanning restriction fragments...")
        QApplication.processEvents()

        params = {
            "enzyme1": self._enzyme1_combo.currentText(),
            "enzyme2": self._enzyme2_combo.currentText() if self._enzyme2_combo.currentText() != "None" else None,
            "min_frag": self._min_frag_spin.value(),
            "max_frag": self._max_frag_spin.value(),
            "min_end_dist": self._min_end_dist_spin.value(),
        }

        self._worker = GBSREWorker(self.state.ssrs, self.state.genome, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_discovery_done)
        self._worker.error.connect(self._on_discovery_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self._progress_bar.setValue(done)
        self._status_label.setText(f"Scanning SSRs: {done:,} / {total:,}")
        QApplication.processEvents()

    def _on_discovery_done(self, result):
        qualified = result["qualified_ssrs"]
        # Build all markers list for chromosome view
        self._all_markers = []
        for ssr in self.state.ssrs:
            self._all_markers.append({"contig": ssr.get("contig"), "start": ssr.get("start", 0)})

        if self._ld_filter_cb.isChecked():
            from core.ld_filter import thin_markers_by_distance
            qualified = thin_markers_by_distance(qualified, self._ld_distance_spin.value())
            result["qualified_ssrs"] = qualified
            result["n_qualified"] = len(qualified)

        self._qualified_ssrs = qualified
        self._progress_bar.setVisible(False)
        self._run_btn.setEnabled(True)
        self._status_label.setText("Discovery complete.")
        self._tab_widget.setVisible(True)
        self._barcode_group.setVisible(True)

        self._summary_label.setText(
            f"✅ {result['n_qualified']} SSRs found in {result['passing_fragments']} passing fragments "
            f"(out of {len(self.state.ssrs)} total SSRs)"
        )

        self._populate_table(self._qualified_ssrs)
        # Update chromosome view
        filtered_markers = [{"contig": m.get("contig"), "start": m.get("start")} for m in qualified]
        self._chromosome_view.set_data(self._all_markers, filtered_markers, self._contig_lengths)
        self._chromosome_view.set_show_filtered(self._ld_filter_cb.isChecked())

    def _on_discovery_error(self, msg):
        self._progress_bar.setVisible(False)
        self._run_btn.setEnabled(True)
        self._status_label.setText(f"Error: {msg}")

    def _populate_table(self, ssrs):
        if not ssrs:
            self._results_table.setRowCount(0)
            return
        self._results_table.setRowCount(len(ssrs))
        cols = ["SSR ID", "Contig", "Motif", "Fragment Start", "Fragment End", "Fragment Size"]
        self._results_table.setColumnCount(len(cols))
        self._results_table.setHorizontalHeaderLabels(cols)
        for i, s in enumerate(ssrs):
            self._results_table.setItem(i, 0, QTableWidgetItem(str(s.get("ssr_id"))))
            self._results_table.setItem(i, 1, QTableWidgetItem(s.get("contig", "")))
            self._results_table.setItem(i, 2, QTableWidgetItem(s.get("motif", "")))
            self._results_table.setItem(i, 3, QTableWidgetItem(str(s.get("fragment_start", ""))))
            self._results_table.setItem(i, 4, QTableWidgetItem(str(s.get("fragment_end", ""))))
            self._results_table.setItem(i, 5, QTableWidgetItem(str(s.get("fragment_size", ""))))
        self._results_table.resizeColumnsToContents()

    def _generate_barcodes(self):
        n = self._n_samples_spin.value()
        length = self._barcode_len_spin.value()
        min_dist = self._min_dist_spin.value()

        self._gen_barcode_btn.setEnabled(False)
        self._barcode_worker = BarcodeWorker(n, length, min_dist)
        self._barcode_worker.finished.connect(self._on_barcodes_generated)
        self._barcode_worker.error.connect(self._on_barcode_error)
        self._barcode_worker.start()

    def _on_barcodes_generated(self, barcodes):
        self._gen_barcode_btn.setEnabled(True)
        self._barcode_text.setText("\n".join(f"{i+1}\t{bc}" for i, bc in enumerate(barcodes)))

    def _on_barcode_error(self, msg):
        self._gen_barcode_btn.setEnabled(True)
        self._status_label.setText(f"Barcode error: {msg}")

    def _export_csv(self):
        if not self._qualified_ssrs:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Markers", "gbs_markers.csv", "CSV (*.csv)")
        if not path:
            return
        import pandas as pd
        pd.DataFrame(self._qualified_ssrs).to_csv(path, index=False)
        self.mw.set_status(f"Markers exported to {path}")

    def _export_bed(self):
        if not self._qualified_ssrs:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export BED", "gbs_fragments.bed", "BED (*.bed)")
        if not path:
            return
        with open(path, 'w') as f:
            for s in self._qualified_ssrs:
                f.write(f"{s['contig']}\t{s['fragment_start']-1}\t{s['fragment_end']}\tSSR{s['ssr_id']}\n")
        self.mw.set_status(f"BED file exported to {path}")

    def _export_barcodes(self):
        text = self._barcode_text.toPlainText()
        if not text.strip():
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Barcodes", "barcodes.txt", "Text (*.txt)")
        if not path:
            return
        with open(path, 'w') as f:
            f.write(text)
        self.mw.set_status(f"Barcodes exported to {path}")

    def _rebuild(self):
        has_genome = self.state.has_genome
        has_ssrs = self.state.has_ssrs
        self._run_btn.setEnabled(has_genome and has_ssrs)
        if not has_genome:
            self._status_label.setText("Load a genome first.")
        elif not has_ssrs:
            self._status_label.setText("Run SSR detection first.")
        else:
            self._status_label.setText(f"{len(self.state.ssrs):,} SSRs ready for GBS‑RE discovery.")

    def on_show(self):
        self._compute_contig_lengths()
        self._rebuild()