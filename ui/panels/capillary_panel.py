"""
capillary_panel.py — Capillary Fragment Analysis Panel
Design multiplex panels for ABI genetic analyzers.
Includes LD filter and chromosome view.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QFileDialog, QSpinBox,
    QComboBox, QCheckBox, QToolButton, QApplication,
    QTabWidget, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QPointF
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)
from core.capillary_multiplex import DYE_SETS, filter_unique_loci


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
        self._markers = []
        self._filtered_markers = []
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
        markers_to_show = self._filtered_markers if self._show_filtered else self._markers
        if not markers_to_show:
            painter = QPainter(self)
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No markers to display")
            return

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        margin_left, margin_right, margin_top, margin_bottom = 100, 20, 20, 30

        by_contig = {}
        for m in markers_to_show:
            contig = m.get("contig", "")
            if contig not in by_contig:
                by_contig[contig] = []
            by_contig[contig].append(m)

        contigs = list(by_contig.keys())
        contig_lengths = {}
        for contig in contigs:
            if contig in self._contig_lengths:
                contig_lengths[contig] = self._contig_lengths[contig]
            else:
                max_pos = max(m.get("start", 0) for m in by_contig[contig])
                contig_lengths[contig] = max_pos + 1000

        contigs = sorted(contigs, key=lambda c: contig_lengths.get(c, 0), reverse=True)[:8]
        n_contigs = len(contigs)
        row_height = (h - margin_top - margin_bottom) / max(n_contigs, 1)
        bar_height = row_height * 0.6

        painter.setPen(QPen(QColor(TEXT_SECONDARY), 1))
        painter.setFont(QFont(FONT_MONO, 8))

        for i, contig in enumerate(contigs):
            y = margin_top + i * row_height + (row_height - bar_height) / 2
            contig_len = contig_lengths[contig]

            painter.drawText(5, int(y + bar_height/2 + 4), contig[:15])

            bar_x = margin_left
            bar_w = w - margin_left - margin_right
            painter.setBrush(QBrush(QColor(BG_MID)))
            painter.setPen(QPen(QColor(BORDER), 1))
            painter.drawRect(bar_x, int(y), int(bar_w), int(bar_height))

            painter.setBrush(QBrush(QColor(ACCENT)))
            painter.setPen(Qt.PenStyle.NoPen)
            for m in by_contig[contig]:
                pos = m.get("start", 0)
                if contig_len > 0:
                    rel_x = bar_x + (pos / contig_len) * bar_w
                    painter.drawEllipse(QPointF(rel_x, y + bar_height/2), 3, 3)

            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(bar_x, int(y + bar_height + 12), "0")
            painter.drawText(bar_x + bar_w - 30, int(y + bar_height + 12), f"{contig_len/1000:.1f} kb")


class CollapsibleBox(QWidget):
    def __init__(self, title="", parent=None):
        super().__init__(parent)
        self._is_expanded = False
        self._title = title
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.toggle_btn = QToolButton()
        self.toggle_btn.setStyleSheet(f"""
            QToolButton {{
                background: {BG_MID};
                border: 1px solid {BORDER};
                border-radius: 4px;
                padding: 8px 12px;
                font-weight: bold;
                color: {TEXT_PRIMARY};
                text-align: left;
            }}
            QToolButton:hover {{ background: {BG_LIGHT}; }}
        """)
        self.toggle_btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.toggle_btn.setArrowType(Qt.ArrowType.RightArrow)
        self.toggle_btn.setText(self._title)
        self.toggle_btn.clicked.connect(self._on_toggle)

        self.content_area = QWidget()
        self.content_area.setVisible(False)
        content_layout = QVBoxLayout(self.content_area)
        content_layout.setContentsMargins(12, 8, 12, 8)
        self.content_layout = content_layout

        layout.addWidget(self.toggle_btn)
        layout.addWidget(self.content_area)

    def _on_toggle(self):
        self.set_expanded(not self._is_expanded)

    def set_expanded(self, expanded: bool):
        self._is_expanded = expanded
        self.toggle_btn.setArrowType(Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow)
        self.content_area.setVisible(expanded)

    def addWidget(self, widget):
        self.content_layout.addWidget(widget)

    def addLayout(self, layout):
        self.content_layout.addLayout(layout)


class CapillaryWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, primers, params):
        super().__init__()
        self.primers = primers
        self.params = params

    def run(self):
        try:
            from core.capillary_multiplex import assign_dyes_advanced, DYE_SETS
            from core.amplicon_sizing import calculate_allele_sizes

            unique_primers = filter_unique_loci(self.primers)

            loci = []
            for p in unique_primers:
                if self.params.get("use_narrow_range", True):
                    ref_count = p.get("repeat_count", 10)
                    min_rep = max(1, ref_count - 2)
                    max_rep = ref_count + 2
                else:
                    min_rep = self.params["min_repeats"]
                    max_rep = self.params["max_repeats"]

                sizing = calculate_allele_sizes(p, min_repeats=min_rep, max_repeats=max_rep)
                loci.append({
                    "ssr_id": p["ssr_id"],
                    "min_allele_size": sizing["min_allele_size"],
                    "max_allele_size": sizing["max_allele_size"],
                    "primer": p,
                })

            dye_set_name = self.params["dye_set"]
            dye_list = DYE_SETS[dye_set_name]

            result = assign_dyes_advanced(
                loci,
                dye_list,
                min_bin_spacing=self.params["min_buffer"],
                max_loci_per_dye=self.params["max_per_dye"],
                allow_overlap_within_dye=self.params["allow_overlap"],
                overlap_tolerance=self.params["overlap_tolerance"],
            )

            locus_map = {loc["ssr_id"]: loc["primer"] for loc in loci}
            for a in result["assignments"]:
                a.update(locus_map[a["ssr_id"]])

            result["unique_loci_used"] = len(loci)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class SuggestWorker(QThread):
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, primers, params):
        super().__init__()
        self.primers = primers
        self.params = params

    def run(self):
        try:
            from core.capillary_multiplex import suggest_best_panel, DYE_SETS, filter_unique_loci
            from core.amplicon_sizing import calculate_allele_sizes

            unique_primers = filter_unique_loci(self.primers)

            loci = []
            for p in unique_primers:
                if self.params.get("use_narrow_range", True):
                    ref_count = p.get("repeat_count", 10)
                    min_rep = max(1, ref_count - 2)
                    max_rep = ref_count + 2
                else:
                    min_rep = self.params["min_repeats"]
                    max_rep = self.params["max_repeats"]

                sizing = calculate_allele_sizes(p, min_repeats=min_rep, max_repeats=max_rep)
                loci.append({
                    "ssr_id": p["ssr_id"],
                    "min_allele_size": sizing["min_allele_size"],
                    "max_allele_size": sizing["max_allele_size"],
                    "primer": p,
                })

            dye_set_name = self.params["dye_set"]
            dye_list = DYE_SETS[dye_set_name]

            best = suggest_best_panel(
                loci,
                dye_list,
                min_spacing_options=self.params["spacings"],
                max_per_dye_options=self.params["max_per_dye_options"],
            )

            locus_map = {loc["ssr_id"]: loc["primer"] for loc in loci}
            for a in best["assignments"]:
                a.update(locus_map[a["ssr_id"]])

            best["suggested"] = True
            best["unique_loci_used"] = len(loci)
            self.finished.emit(best)
        except Exception as e:
            self.error.emit(str(e))


class CapillaryPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw = main_window
        self._worker = None
        self._last_version = -1
        self._pass_primers = []
        self._result = None
        self._all_markers = []
        self._contig_lengths = {}
        self._build_ui()
        self._compute_contig_lengths()

    def _compute_contig_lengths(self):
        if self.state.genome:
            self._contig_lengths = {c: len(s) for c, s in self.state.genome.items()}

    def _toggle_overlap_tolerance(self, checked):
        self._overlap_tolerance_spin.setEnabled(checked)

    def _toggle_ld_filter(self, checked):
        self._ld_distance_spin.setEnabled(checked)

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

        title = QLabel("Capillary Fragment Analysis")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Design multiplex panels for ABI genetic analyzers")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        self._placeholder = QLabel(
            "Capillary Tools are only available when primers were designed in Capillary Electrophoresis mode.\n"
            "Go to Primer Design, select 'Capillary Electrophoresis', and run primer design first."
        )
        self._placeholder.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._placeholder.setWordWrap(True)
        L.addWidget(self._placeholder)

        self._main_group = QGroupBox("Configuration")
        self._main_group.setVisible(False)
        mg = QVBoxLayout(self._main_group)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Dye Set:"))
        self._dye_set_combo = QComboBox()
        self._dye_set_combo.addItems(list(DYE_SETS.keys()))
        self._dye_set_combo.setCurrentText("4 Dyes (FAM, VIC, NED, PET)")
        row1.addWidget(self._dye_set_combo)
        row1.addStretch()
        mg.addLayout(row1)

        row_range = QHBoxLayout()
        self._use_narrow_cb = QCheckBox("Use narrow allele range (±2 repeats around reference)")
        self._use_narrow_cb.setChecked(True)
        self._use_narrow_cb.setToolTip("Reduces overlap by using conservative allele size estimates")
        row_range.addWidget(self._use_narrow_cb)
        mg.addLayout(row_range)

        row_repeats = QHBoxLayout()
        row_repeats.addWidget(QLabel("Min expected repeats:"))
        self._min_repeats_spin = QSpinBox()
        self._min_repeats_spin.setRange(1, 20)
        self._min_repeats_spin.setValue(3)
        row_repeats.addWidget(self._min_repeats_spin)

        row_repeats.addWidget(QLabel("Max expected repeats:"))
        self._max_repeats_spin = QSpinBox()
        self._max_repeats_spin.setRange(10, 200)
        self._max_repeats_spin.setValue(20)
        row_repeats.addWidget(self._max_repeats_spin)
        row_repeats.addStretch()
        mg.addLayout(row_repeats)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Min bin spacing (bp):"))
        self._min_buffer_spin = QSpinBox()
        self._min_buffer_spin.setRange(0, 50)
        self._min_buffer_spin.setValue(5)
        self._min_buffer_spin.setToolTip("Minimum gap between non‑overlapping bins")
        row2.addWidget(self._min_buffer_spin)

        row2.addWidget(QLabel("Max loci per dye:"))
        self._max_per_dye_spin = QSpinBox()
        self._max_per_dye_spin.setRange(0, 50)
        self._max_per_dye_spin.setValue(20)
        self._max_per_dye_spin.setToolTip("Limit number of loci per dye (0 = unlimited)")
        row2.addWidget(self._max_per_dye_spin)
        row2.addStretch()
        mg.addLayout(row2)

        row3 = QHBoxLayout()
        self._allow_overlap_cb = QCheckBox("Allow overlapping bins within same dye")
        self._allow_overlap_cb.setChecked(False)
        self._allow_overlap_cb.toggled.connect(self._toggle_overlap_tolerance)
        row3.addWidget(self._allow_overlap_cb)

        row3.addWidget(QLabel("Overlap tolerance (bp):"))
        self._overlap_tolerance_spin = QSpinBox()
        self._overlap_tolerance_spin.setRange(0, 50)
        self._overlap_tolerance_spin.setValue(5)
        self._overlap_tolerance_spin.setEnabled(False)
        row3.addWidget(self._overlap_tolerance_spin)
        row3.addStretch()
        mg.addLayout(row3)

        row4 = QHBoxLayout()
        self._add_m13_cb = QCheckBox("Add M13 tails to exported primers")
        self._add_m13_cb.setChecked(True)
        row4.addWidget(self._add_m13_cb)
        mg.addLayout(row4)

        self._ld_filter_cb = QCheckBox("Apply LD filter (thin markers by distance)")
        self._ld_filter_cb.setChecked(False)
        self._ld_filter_cb.toggled.connect(self._toggle_ld_filter)
        mg.addWidget(self._ld_filter_cb)

        ld_row = QHBoxLayout()
        ld_row.addWidget(QLabel("Min distance (bp):"))
        self._ld_distance_spin = QSpinBox()
        self._ld_distance_spin.setRange(1000, 1000000)
        self._ld_distance_spin.setValue(10000)
        self._ld_distance_spin.setEnabled(False)
        ld_row.addWidget(self._ld_distance_spin)
        ld_row.addStretch()
        mg.addLayout(ld_row)

        btn_layout = QHBoxLayout()
        self._run_btn = QPushButton("Design Multiplex Panel")
        self._run_btn.setObjectName("primary")
        self._run_btn.clicked.connect(self._run_assignment)
        btn_layout.addWidget(self._run_btn)

        self._suggest_btn = QPushButton("Suggest Best Panel")
        self._suggest_btn.clicked.connect(self._run_suggestion)
        btn_layout.addWidget(self._suggest_btn)

        btn_layout.addStretch()
        mg.addLayout(btn_layout)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        mg.addWidget(self._status_label)

        L.addWidget(self._main_group)

        self._tab_widget = QTabWidget()
        self._tab_widget.setVisible(False)

        results_tab = QWidget()
        rt_layout = QVBoxLayout(results_tab)
        self._summary_label = QLabel("")
        self._summary_label.setWordWrap(True)
        self._summary_label.setStyleSheet(f"font-weight: bold; color: {TEXT_PRIMARY};")
        rt_layout.addWidget(self._summary_label)

        self._assign_table = QTableWidget()
        self._assign_table.setAlternatingRowColors(True)
        self._assign_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._assign_table.setSortingEnabled(True)
        self._assign_table.setFixedHeight(400)
        rt_layout.addWidget(self._assign_table)

        export_layout = QHBoxLayout()
        self._export_csv_btn = QPushButton("Export Panel (CSV)")
        self._export_csv_btn.clicked.connect(self._export_csv)
        export_layout.addWidget(self._export_csv_btn)
        self._export_genemapper_btn = QPushButton("Export GeneMapper Panel")
        self._export_genemapper_btn.clicked.connect(self._export_genemapper)
        export_layout.addWidget(self._export_genemapper_btn)
        self._export_genalex_btn = QPushButton("Export GenAlEx Format")
        self._export_genalex_btn.clicked.connect(self._export_genalex)
        export_layout.addWidget(self._export_genalex_btn)
        export_layout.addStretch()
        rt_layout.addLayout(export_layout)

        self._tab_widget.addTab(results_tab, "Assignments")

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
        L.addStretch()

    def _run_assignment(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return

        self._run_btn.setEnabled(False)
        self._suggest_btn.setEnabled(False)
        self._status_label.setText("Assigning dyes...")
        QApplication.processEvents()

        max_per_dye = self._max_per_dye_spin.value()
        params = {
            "dye_set": self._dye_set_combo.currentText(),
            "min_buffer": self._min_buffer_spin.value(),
            "max_per_dye": max_per_dye if max_per_dye > 0 else None,
            "allow_overlap": self._allow_overlap_cb.isChecked(),
            "overlap_tolerance": self._overlap_tolerance_spin.value(),
            "min_repeats": self._min_repeats_spin.value(),
            "max_repeats": self._max_repeats_spin.value(),
            "use_narrow_range": self._use_narrow_cb.isChecked(),
        }

        self._worker = CapillaryWorker(self._pass_primers, params)
        self._worker.finished.connect(self._on_assignment_done)
        self._worker.error.connect(self._on_assignment_error)
        self._worker.start()

    def _run_suggestion(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return

        self._run_btn.setEnabled(False)
        self._suggest_btn.setEnabled(False)
        self._status_label.setText("Trying parameter combinations...")
        QApplication.processEvents()

        params = {
            "dye_set": self._dye_set_combo.currentText(),
            "min_repeats": self._min_repeats_spin.value(),
            "max_repeats": self._max_repeats_spin.value(),
            "use_narrow_range": self._use_narrow_cb.isChecked(),
            "spacings": [5, 8, 10, 12, 15],
            "max_per_dye_options": [10, 15, 20, 25, 0],
        }

        self._worker = SuggestWorker(self._pass_primers, params)
        self._worker.finished.connect(self._on_suggestion_done)
        self._worker.error.connect(self._on_assignment_error)
        self._worker.start()

    def _on_assignment_done(self, result):
        self._result = result
        self._run_btn.setEnabled(True)
        self._suggest_btn.setEnabled(True)
        self._status_label.setText("Assignment complete.")
        self._tab_widget.setVisible(True)

        assignments = result["assignments"]
        if self._ld_filter_cb.isChecked():
            from core.ld_filter import thin_markers_by_distance
            assignments = thin_markers_by_distance(assignments, self._ld_distance_spin.value())
            result["assignments"] = assignments
            result["n_assigned"] = len(assignments)

        unique_total = result.get("unique_loci_used", result["n_total"])
        self._summary_label.setText(
            f"✅ {result['n_assigned']} of {unique_total} unique loci assigned. "
            f"{result['n_unassigned']} unassigned."
        )
        self._populate_table(assignments)

        self._all_markers = []
        for p in self._pass_primers:
            self._all_markers.append({"contig": p.get("contig"), "start": p.get("start", 0)})
        filtered_markers = [{"contig": m.get("contig"), "start": m.get("start")} for m in assignments]
        self._chromosome_view.set_data(self._all_markers, filtered_markers, self._contig_lengths)
        self._chromosome_view.set_show_filtered(self._ld_filter_cb.isChecked())

    def _on_suggestion_done(self, result):
        self._result = result
        self._run_btn.setEnabled(True)
        self._suggest_btn.setEnabled(True)
        self._status_label.setText("Suggestion complete.")
        self._tab_widget.setVisible(True)

        assignments = result["assignments"]
        if self._ld_filter_cb.isChecked():
            from core.ld_filter import thin_markers_by_distance
            assignments = thin_markers_by_distance(assignments, self._ld_distance_spin.value())
            result["assignments"] = assignments
            result["n_assigned"] = len(assignments)

        params = result.get("best_params", {})
        unique_total = result.get("unique_loci_used", result["n_total"])
        self._summary_label.setText(
            f"✅ Best panel: {result['n_assigned']} of {unique_total} loci assigned.\n"
            f"Parameters: spacing={params.get('min_bin_spacing')} bp, "
            f"max per dye={params.get('max_loci_per_dye') or 'unlimited'}"
        )
        self._populate_table(assignments)

        self._all_markers = []
        for p in self._pass_primers:
            self._all_markers.append({"contig": p.get("contig"), "start": p.get("start", 0)})
        filtered_markers = [{"contig": m.get("contig"), "start": m.get("start")} for m in assignments]
        self._chromosome_view.set_data(self._all_markers, filtered_markers, self._contig_lengths)
        self._chromosome_view.set_show_filtered(self._ld_filter_cb.isChecked())

    def _on_assignment_error(self, msg):
        self._run_btn.setEnabled(True)
        self._suggest_btn.setEnabled(True)
        self._status_label.setText(f"Error: {msg}")

    def _populate_table(self, assignments):
        if not assignments:
            return
        self._assign_table.setRowCount(len(assignments))
        cols = ["SSR ID", "Dye", "Min Size", "Max Size", "Forward Primer", "Reverse Primer"]
        self._assign_table.setColumnCount(len(cols))
        self._assign_table.setHorizontalHeaderLabels(cols)
        for i, a in enumerate(assignments):
            self._assign_table.setItem(i, 0, QTableWidgetItem(str(a["ssr_id"])))
            self._assign_table.setItem(i, 1, QTableWidgetItem(a["dye"]))
            self._assign_table.setItem(i, 2, QTableWidgetItem(str(a["min_size"])))
            self._assign_table.setItem(i, 3, QTableWidgetItem(str(a["max_size"])))
            fwd = a.get("left_primer", "")
            rev = a.get("right_primer", "")
            if self._add_m13_cb.isChecked():
                from core.m13_tails import M13_FORWARD, PIG_TAIL, M13_REVERSE
                fwd = M13_FORWARD + fwd
                rev = PIG_TAIL + M13_REVERSE + rev
            self._assign_table.setItem(i, 4, QTableWidgetItem(fwd))
            self._assign_table.setItem(i, 5, QTableWidgetItem(rev))
        self._assign_table.resizeColumnsToContents()

    def _export_csv(self):
        if not self._result:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Panel", "capillary_panel.csv", "CSV (*.csv)")
        if not path:
            return
        import pandas as pd
        df = pd.DataFrame(self._result["assignments"])
        df.to_csv(path, index=False)
        self.mw.set_status(f"Panel exported to {path}")

    def _export_genemapper(self):
        self._status_label.setText("GeneMapper export not yet implemented")

    def _export_genalex(self):
        self._status_label.setText("GenAlEx export not yet implemented")

    def _get_pass_primers(self):
        spec = self.state.specificity_results or []
        primers = self.state.primer_results or []
        pass_ids = {r["ssr_id"] for r in spec if r.get("specificity_status") == "PASS"}
        return [p for p in primers if p.get("ssr_id") in pass_ids]

    def _rebuild(self):
        mode = self.state.workflow_mode
        is_capillary = (mode == "capillary")
        has_spec = self.state.has_specificity

        should_show = is_capillary and has_spec
        self._placeholder.setVisible(not should_show)
        self._main_group.setVisible(should_show)
        self._tab_widget.setVisible(False)

        if should_show:
            self._pass_primers = self._get_pass_primers()
            unique_count = len({p["ssr_id"] for p in self._pass_primers})
            self._status_label.setText(f"{unique_count} PASS loci ready.")

    def on_show(self):
        self._compute_contig_lengths()
        self._rebuild()