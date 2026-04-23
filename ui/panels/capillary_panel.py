"""
capillary_panel.py — Capillary Fragment Analysis Panel
Design multiplex panels for ABI genetic analyzers.
Supports multiple panels (injections) with per‑panel tabs and full dye sets.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QFileDialog, QSpinBox,
    QComboBox, QCheckBox, QToolButton, QApplication,
    QTabWidget, QSizePolicy, QGridLayout,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QPointF
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)
from core.capillary_multiplex import DYE_SETS, filter_unique_loci, assign_dyes_advanced, assign_to_panels_sequential


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
            from core.amplicon_sizing import calculate_allele_sizes

            unique_primers = filter_unique_loci(self.primers)

            loci = []
            for p in unique_primers:
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

            max_panels = self.params.get("max_panels", 1)

            if max_panels is None or max_panels <= 1:
                assignments, unassigned, dye_bins = assign_dyes_advanced(
                    loci, dye_list,
                    min_bin_spacing=self.params["min_buffer"],
                    max_loci_per_dye=self.params.get("max_per_dye"),
                    allow_overlap_within_dye=self.params.get("allow_overlap", False),
                    overlap_tolerance=self.params.get("overlap_tolerance", 0),
                )
                for a in assignments:
                    a["panel"] = 1
                result = {
                    "assignments": assignments,
                    "panels": {1: assignments},
                    "n_panels": 1,
                    "n_total": len(loci),
                    "n_assigned": len(assignments),
                    "unassigned": unassigned,
                }
            else:
                result = assign_to_panels_sequential(
                    loci, dye_list,
                    max_panels=max_panels,
                    min_bin_spacing=self.params["min_buffer"],
                    max_loci_per_dye=self.params.get("max_per_dye"),
                    allow_overlap=self.params.get("allow_overlap", False),
                    overlap_tolerance=self.params.get("overlap_tolerance", 0),
                )

            locus_map = {loc["ssr_id"]: loc["primer"] for loc in loci}
            for a in result["assignments"]:
                a.update(locus_map[a["ssr_id"]])

            result["unique_loci_used"] = len(loci)
            self.finished.emit(result)
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
        mg.setSpacing(10)

        # ---- Basic Settings ----
        basic_group = QGroupBox("Basic Settings")
        basic_layout = QGridLayout(basic_group)
        basic_layout.setVerticalSpacing(8)

        basic_layout.addWidget(QLabel("Dye Set:"), 0, 0)
        self._dye_set_combo = QComboBox()
        self._dye_set_combo.addItems(list(DYE_SETS.keys()))
        self._dye_set_combo.setCurrentText("4 Dyes (FAM, VIC, NED, PET)")
        basic_layout.addWidget(self._dye_set_combo, 0, 1)

        basic_layout.addWidget(QLabel("Min bin spacing (bp):"), 1, 0)
        self._min_buffer_spin = QSpinBox()
        self._min_buffer_spin.setRange(0, 50)
        self._min_buffer_spin.setValue(5)
        self._min_buffer_spin.setToolTip("Minimum gap between non‑overlapping bins")
        basic_layout.addWidget(self._min_buffer_spin, 1, 1)

        basic_layout.addWidget(QLabel("Max loci per dye:"), 2, 0)
        self._max_per_dye_spin = QSpinBox()
        self._max_per_dye_spin.setRange(0, 50)
        self._max_per_dye_spin.setValue(20)
        self._max_per_dye_spin.setToolTip("Limit number of loci per dye (0 = unlimited)")
        basic_layout.addWidget(self._max_per_dye_spin, 2, 1)

        mg.addWidget(basic_group)

        # ---- Advanced Settings ----
        adv_group = QGroupBox("Advanced Settings")
        adv_layout = QGridLayout(adv_group)
        adv_layout.setVerticalSpacing(8)

        adv_layout.addWidget(QLabel("Min expected repeats:"), 0, 0)
        self._min_repeats_spin = QSpinBox()
        self._min_repeats_spin.setRange(1, 20)
        self._min_repeats_spin.setValue(3)
        adv_layout.addWidget(self._min_repeats_spin, 0, 1)

        adv_layout.addWidget(QLabel("Max expected repeats:"), 0, 2)
        self._max_repeats_spin = QSpinBox()
        self._max_repeats_spin.setRange(10, 200)
        self._max_repeats_spin.setValue(40)
        adv_layout.addWidget(self._max_repeats_spin, 0, 3)

        self._allow_overlap_cb = QCheckBox("Allow overlapping bins within same dye")
        self._allow_overlap_cb.setChecked(False)
        self._allow_overlap_cb.toggled.connect(self._toggle_overlap_tolerance)
        adv_layout.addWidget(self._allow_overlap_cb, 1, 0, 1, 2)

        adv_layout.addWidget(QLabel("Overlap tolerance (bp):"), 1, 2)
        self._overlap_tolerance_spin = QSpinBox()
        self._overlap_tolerance_spin.setRange(0, 50)
        self._overlap_tolerance_spin.setValue(5)
        self._overlap_tolerance_spin.setEnabled(False)
        adv_layout.addWidget(self._overlap_tolerance_spin, 1, 3)

        mg.addWidget(adv_group)

        # ---- Multiple Panels ----
        panel_group = QGroupBox("Multiple Panels (Injections)")
        panel_layout = QVBoxLayout(panel_group)

        self._split_panels_cb = QCheckBox("Split into multiple panels")
        self._split_panels_cb.setChecked(True)
        panel_layout.addWidget(self._split_panels_cb)

        panel_opts = QHBoxLayout()
        panel_opts.addWidget(QLabel("Number of panels:"))
        self._num_panels_spin = QSpinBox()
        self._num_panels_spin.setRange(1, 10)
        self._num_panels_spin.setValue(2)
        panel_opts.addWidget(self._num_panels_spin)
        panel_opts.addStretch()
        panel_layout.addLayout(panel_opts)

        mg.addWidget(panel_group)

        # ---- LD Filter ----
        ld_group = QGroupBox("Linkage Disequilibrium Filter")
        ld_layout = QVBoxLayout(ld_group)
        self._ld_filter_cb = QCheckBox("Apply LD filter (thin markers by distance)")
        self._ld_filter_cb.setChecked(False)
        self._ld_filter_cb.toggled.connect(self._toggle_ld_filter)
        ld_layout.addWidget(self._ld_filter_cb)

        ld_row = QHBoxLayout()
        ld_row.addWidget(QLabel("Min distance (bp):"))
        self._ld_distance_spin = QSpinBox()
        self._ld_distance_spin.setRange(1000, 1000000)
        self._ld_distance_spin.setValue(10000)
        self._ld_distance_spin.setEnabled(False)
        ld_row.addWidget(self._ld_distance_spin)
        ld_row.addStretch()
        ld_layout.addLayout(ld_row)

        mg.addWidget(ld_group)

        # ---- M13 Tails ----
        self._add_m13_cb = QCheckBox("Add M13 tails to exported primers")
        self._add_m13_cb.setChecked(True)
        mg.addWidget(self._add_m13_cb)

        # ---- Action Buttons ----
        btn_layout = QHBoxLayout()
        self._run_btn = QPushButton("Design Multiplex Panel")
        self._run_btn.setObjectName("primary")
        self._run_btn.clicked.connect(self._run_assignment)
        btn_layout.addWidget(self._run_btn)
        btn_layout.addStretch()
        mg.addLayout(btn_layout)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        mg.addWidget(self._status_label)

        L.addWidget(self._main_group)

        self._main_tab_widget = QTabWidget()
        self._main_tab_widget.setVisible(False)
        L.addWidget(self._main_tab_widget)

        L.addStretch()

    def _run_assignment(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return

        self._run_btn.setEnabled(False)
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
        }

        if self._split_panels_cb.isChecked():
            params["max_panels"] = self._num_panels_spin.value()
        else:
            params["max_panels"] = None

        self._worker = CapillaryWorker(self._pass_primers, params)
        self._worker.finished.connect(self._on_assignment_done)
        self._worker.error.connect(self._on_assignment_error)
        self._worker.start()

    def _apply_ld_filter(self, assignments):
        if not self._ld_filter_cb.isChecked():
            return assignments
        from core.ld_filter import thin_markers_by_distance
        return thin_markers_by_distance(assignments, self._ld_distance_spin.value())

    def _on_assignment_done(self, result):
        self._result = result
        self.state.capillary_result = result
        self._run_btn.setEnabled(True)
        self._status_label.setText("Assignment complete.")
        self._display_results(result)

    def _display_results(self, result):
        self._main_tab_widget.clear()
        self._main_tab_widget.setVisible(True)

        panels = result.get("panels", {1: result["assignments"]})
        n_panels = result.get("n_panels", 1)

        # Summary tab with export button
        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        unique_total = result.get("unique_loci_used", result["n_total"])
        panel_text = f" across {n_panels} panel{'s' if n_panels != 1 else ''}" if n_panels > 1 else ""
        summary_label = QLabel(
            f"✅ {result['n_assigned']} of {unique_total} unique loci assigned{panel_text}. "
            f"{len(result.get('unassigned', []))} unassigned."
        )
        summary_label.setWordWrap(True)
        summary_label.setStyleSheet(f"font-weight: bold; color: {TEXT_PRIMARY};")
        summary_layout.addWidget(summary_label)

        export_all_btn = QPushButton("Export All Panels (CSV)")
        export_all_btn.clicked.connect(lambda: self._export_all_panels(result))
        summary_layout.addWidget(export_all_btn, alignment=Qt.AlignmentFlag.AlignRight)
        summary_layout.addStretch()
        self._main_tab_widget.addTab(summary_tab, "Summary")

        # Per-panel tabs
        for panel_idx in sorted(panels.keys()):
            assignments = self._apply_ld_filter(panels[panel_idx])
            if not assignments:
                continue

            panel_tab = QWidget()
            panel_layout = QVBoxLayout(panel_tab)

            table = QTableWidget()
            table.setAlternatingRowColors(True)
            table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            table.setSortingEnabled(True)

            cols = ["SSR ID", "Dye", "Min Size", "Max Size", "Motif", "Contig"]
            show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
            if show_chromosome:
                cols.append("Chromosome")
            cols.extend(["Forward Primer", "Reverse Primer"])

            table.setColumnCount(len(cols))
            table.setHorizontalHeaderLabels(cols)
            table.setRowCount(len(assignments))

            for i, a in enumerate(assignments):
                col_idx = 0
                table.setItem(i, col_idx, QTableWidgetItem(str(a["ssr_id"])))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(a["dye"]))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a["min_size"])))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a["max_size"])))
                col_idx += 1
                motif = a.get("motif", "")
                table.setItem(i, col_idx, QTableWidgetItem(motif))
                col_idx += 1
                contig = a.get("contig", "")
                table.setItem(i, col_idx, QTableWidgetItem(contig))
                col_idx += 1
                if show_chromosome:
                    chrom = self.state.get_display_name(contig)
                    table.setItem(i, col_idx, QTableWidgetItem(chrom))
                    col_idx += 1
                fwd = a.get("left_primer", "")
                rev = a.get("right_primer", "")
                if self._add_m13_cb.isChecked():
                    from core.m13_tails import M13_FORWARD, PIG_TAIL, M13_REVERSE
                    fwd = M13_FORWARD + fwd
                    rev = PIG_TAIL + M13_REVERSE + rev
                table.setItem(i, col_idx, QTableWidgetItem(fwd))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(rev))

            table.resizeColumnsToContents()
            panel_layout.addWidget(table)
            self._main_tab_widget.addTab(panel_tab, f"Panel {panel_idx}")

        # Chromosome View tab
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
        self._main_tab_widget.addTab(chrom_tab, "Chromosome View")

        self._all_markers = []
        for p in self._pass_primers:
            self._all_markers.append({"contig": p.get("contig"), "start": p.get("start", 0)})
        filtered_markers = [{"contig": m.get("contig"), "start": m.get("start")} for m in self._apply_ld_filter(result["assignments"])]
        self._chromosome_view.set_data(self._all_markers, filtered_markers, self._contig_lengths)
        self._chromosome_view.set_show_filtered(self._ld_filter_cb.isChecked())

    def _export_all_panels(self, result):
        all_assignments = result.get("assignments", [])
        if not all_assignments:
            self._status_label.setText("No assignments to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export All Panels", "capillary_all_panels.csv", "CSV (*.csv)")
        if not path:
            return
        import pandas as pd
        df = pd.DataFrame(all_assignments)
        if hasattr(self.state, 'get_display_name') and "contig" in df.columns:
            df["Chromosome"] = df["contig"].apply(lambda c: self.state.get_display_name(c))
        # Reorder columns for better readability
        col_order = ["panel", "ssr_id", "dye", "min_size", "max_size", "motif", "contig"]
        if "Chromosome" in df.columns:
            col_order.append("Chromosome")
        col_order.extend(["left_primer", "right_primer"])
        df = df[[c for c in col_order if c in df.columns]]
        df.to_csv(path, index=False)
        self.mw.set_status(f"All panels exported to {path}")

    def _on_assignment_error(self, msg):
        self._run_btn.setEnabled(True)
        self._status_label.setText(f"Error: {msg}")

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
        self._main_tab_widget.setVisible(False)

        if should_show:
            self._pass_primers = self._get_pass_primers()
            unique_count = len({p["ssr_id"] for p in self._pass_primers})
            self._status_label.setText(f"{unique_count} PASS loci ready.")

    def on_show(self):
        self._compute_contig_lengths()
        self._rebuild()