"""
amplicon_panel.py — Amplicon Sequencing Panel
Validate multiplex compatibility and optimise your amplicon panel.
Includes LD filter and chromosome view.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QFileDialog, QSpinBox,
    QDoubleSpinBox, QCheckBox, QComboBox, QToolButton,
    QFrame, QApplication, QProgressBar, QGridLayout, QLineEdit,
    QSizePolicy, QTabWidget,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QPointF
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)
from core.amplicon_sizing import PLATFORMS


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


class HistogramWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._sizes = []
        self.setMinimumHeight(180)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, sizes):
        self._sizes = sizes
        self.update()

    def paintEvent(self, event):
        if not self._sizes:
            painter = QPainter(self)
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No data")
            return

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        margin_left, margin_right, margin_top, margin_bottom = 50, 30, 25, 35
        chart_w = w - margin_left - margin_right
        chart_h = h - margin_top - margin_bottom

        try:
            import numpy as np
            counts, bins = np.histogram(self._sizes, bins=15)
        except ImportError:
            min_sz, max_sz = min(self._sizes), max(self._sizes)
            if min_sz == max_sz:
                bins = [min_sz, max_sz+1]
                counts = [len(self._sizes)]
            else:
                n_bins = 15
                bin_width = (max_sz - min_sz) / n_bins
                bins = [min_sz + i * bin_width for i in range(n_bins + 1)]
                counts = [0] * n_bins
                for sz in self._sizes:
                    for i in range(n_bins):
                        if bins[i] <= sz < bins[i+1]:
                            counts[i] += 1
                            break
                    else:
                        if sz == max_sz:
                            counts[-1] += 1

        max_count = max(counts) if counts else 1
        bin_width = chart_w / len(counts)

        painter.setPen(QPen(QColor(TEXT_SECONDARY), 1))
        painter.drawLine(margin_left, h - margin_bottom, w - margin_right, h - margin_bottom)
        painter.drawLine(margin_left, margin_top, margin_left, h - margin_bottom)

        for i, count in enumerate(counts):
            if count == 0:
                continue
            bar_h = (count / max_count) * chart_h
            x = margin_left + i * bin_width
            y = h - margin_bottom - bar_h
            painter.setBrush(QBrush(QColor(ACCENT)))
            painter.setPen(QPen(QColor(ACCENT), 1))
            painter.drawRect(int(x), int(y), int(bin_width - 1), int(bar_h))

        painter.setPen(QColor(TEXT_SECONDARY))
        painter.setFont(QFont(FONT_MONO, 8))
        painter.drawText(margin_left - 10, h - margin_bottom + 15, f"{bins[0]:.0f}")
        painter.drawText(w - margin_right - 25, h - margin_bottom + 15, f"{bins[-1]:.0f}")
        painter.drawText(margin_left - 35, margin_top + 10, f"{max_count}")
        painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL))
        painter.drawText(margin_left + chart_w//2 - 30, h - 5, "Amplicon Size (bp)")
        painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        painter.drawText(margin_left, margin_top - 5, "Size Distribution")


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


class ValidateWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, primers, genome, params):
        super().__init__()
        self.primers = primers
        self.genome = genome
        self.params = params
        self._is_cancelled = False

    def cancel(self):
        self._is_cancelled = True

    def run(self):
        try:
            from core.amplicon_validate import validate_panel
            from core.adapter_tags import TagConfig

            tag_config = TagConfig(
                forward_tag=self.params["forward_tag"],
                reverse_tag=self.params["reverse_tag"],
                include_in_dimer_check=self.params["include_tags_dimer"],
                include_in_amplicon_export=False,
            )

            self.progress.emit("Validating panel...")
            result = validate_panel(
                self.primers,
                self.genome,
                platform_name=self.params["platform"],
                tag_config=tag_config,
                min_repeats=self.params["min_repeats"],
                max_repeats=self.params["max_repeats"],
                max_3prime_score=self.params["max_3prime"],
                max_tm_spread=self.params["max_tm_spread"],
                min_amplicon=self.params["min_amplicon"],
                max_amplicon=self.params["max_amplicon"],
                auto_optimise=False,
            )
            if self._is_cancelled:
                return
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class AmpliconPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw = main_window
        self._worker = None
        self._last_version = -1
        self._pass_primers = []
        self._validation_result = None
        self._all_markers = []
        self._contig_lengths = {}
        self._build_ui()
        self._compute_contig_lengths()

    def _compute_contig_lengths(self):
        if self.state.genome:
            self._contig_lengths = {c: len(s) for c, s in self.state.genome.items()}

    def _on_tag_preset_changed(self, preset: str):
        presets = {
            "Illumina TruSeq": ("AATGATACGGCGACCACCGAGATCTACAC",
                               "CAAGCAGAAGACGGCATACGAGAT"),
            "Illumina Nextera": ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                                "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"),
        }
        if preset == "Custom":
            self._fwd_tag_edit.setReadOnly(False)
            self._rev_tag_edit.setReadOnly(False)
        else:
            self._fwd_tag_edit.setReadOnly(True)
            self._rev_tag_edit.setReadOnly(True)
            fwd, rev = presets[preset]
            self._fwd_tag_edit.setText(fwd)
            self._rev_tag_edit.setText(rev)

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

        title = QLabel("Amplicon Panel Validation")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Validate multiplex compatibility for your amplicon panel")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        self._placeholder = QLabel(
            "Amplicon Tools are only available when primers were designed in Amplicon Sequencing mode.\n"
            "Go to Primer Design, select 'Amplicon Sequencing', and run primer design first."
        )
        self._placeholder.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._placeholder.setWordWrap(True)
        L.addWidget(self._placeholder)

        self._main_group = QGroupBox("Configuration")
        self._main_group.setVisible(False)
        mg = QVBoxLayout(self._main_group)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Sequencing Platform:"))
        self._platform_combo = QComboBox()
        self._platform_combo.addItems([p["name"] for p in PLATFORMS])
        self._platform_combo.setCurrentText("Illumina MiSeq 300bp PE")
        self._platform_combo.currentTextChanged.connect(self._update_platform_limits)
        row1.addWidget(self._platform_combo)
        row1.addStretch()
        mg.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Min amplicon (bp):"))
        self._min_amp_spin = QSpinBox()
        self._min_amp_spin.setRange(30, 1000)
        self._min_amp_spin.setValue(50)
        row2.addWidget(self._min_amp_spin)

        row2.addWidget(QLabel("Max amplicon (bp):"))
        self._max_amp_spin = QSpinBox()
        self._max_amp_spin.setRange(50, 2000)
        self._max_amp_spin.setValue(450)
        row2.addWidget(self._max_amp_spin)

        self._platform_hint = QLabel("(Suggested: 50–450 bp)")
        self._platform_hint.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL-1}pt;")
        row2.addWidget(self._platform_hint)
        row2.addStretch()
        mg.addLayout(row2)

        row3a = QHBoxLayout()
        row3a.addWidget(QLabel("Adapter tag preset:"))
        self._tag_preset_combo = QComboBox()
        self._tag_preset_combo.addItems(["Illumina TruSeq", "Illumina Nextera", "Custom"])
        self._tag_preset_combo.currentTextChanged.connect(self._on_tag_preset_changed)
        row3a.addWidget(self._tag_preset_combo)
        row3a.addStretch()
        mg.addLayout(row3a)

        row3b = QHBoxLayout()
        self._include_tags_dimer_cb = QCheckBox("Include adapter tags in dimer analysis?")
        self._include_tags_dimer_cb.setChecked(False)
        self._include_tags_dimer_cb.setToolTip(
            "If checked, dimer risks are calculated using full primer+tag.\n"
            "Uncheck to only use genomic portions (recommended)."
        )
        row3b.addWidget(self._include_tags_dimer_cb)

        row3b.addWidget(QLabel("Forward tag:"))
        self._fwd_tag_edit = QLineEdit()
        self._fwd_tag_edit.setText("AATGATACGGCGACCACCGAGATCTACAC")
        self._fwd_tag_edit.setMinimumWidth(200)
        self._fwd_tag_edit.setReadOnly(True)
        row3b.addWidget(self._fwd_tag_edit)

        row3b.addWidget(QLabel("Reverse tag:"))
        self._rev_tag_edit = QLineEdit()
        self._rev_tag_edit.setText("CAAGCAGAAGACGGCATACGAGAT")
        self._rev_tag_edit.setMinimumWidth(200)
        self._rev_tag_edit.setReadOnly(True)
        row3b.addWidget(self._rev_tag_edit)
        mg.addLayout(row3b)

        btn_layout = QHBoxLayout()
        self._validate_btn = QPushButton("Validate Panel")
        self._validate_btn.setObjectName("primary")
        self._validate_btn.setMinimumHeight(40)
        self._validate_btn.clicked.connect(self._run_validation)
        btn_layout.addWidget(self._validate_btn)

        self._cancel_btn = QPushButton("Cancel")
        self._cancel_btn.setEnabled(False)
        self._cancel_btn.clicked.connect(self._cancel_worker)
        btn_layout.addWidget(self._cancel_btn)
        btn_layout.addStretch()
        mg.addLayout(btn_layout)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        mg.addWidget(self._status_label)

        self._advanced_box = CollapsibleBox("Advanced Parameters")
        adv_widget = QWidget()
        adv_layout = QGridLayout(adv_widget)
        adv_layout.setVerticalSpacing(8)

        adv_layout.addWidget(_lbl("Max 3' dimer score:"), 0, 0)
        self._max_3prime_spin = QSpinBox()
        self._max_3prime_spin.setRange(2, 10)
        self._max_3prime_spin.setValue(5)
        self._max_3prime_spin.setToolTip("Higher = more tolerant of 3' complementarity")
        adv_layout.addWidget(self._max_3prime_spin, 0, 1)

        adv_layout.addWidget(_lbl("Max Tm spread (°C):"), 0, 2)
        self._max_tm_spread_spin = QDoubleSpinBox()
        self._max_tm_spread_spin.setRange(1.0, 20.0)
        self._max_tm_spread_spin.setValue(5.0)
        self._max_tm_spread_spin.setSingleStep(0.5)
        adv_layout.addWidget(self._max_tm_spread_spin, 0, 3)

        adv_layout.addWidget(_lbl("Min expected repeats:"), 1, 0)
        self._min_repeats_spin = QSpinBox()
        self._min_repeats_spin.setRange(1, 20)
        self._min_repeats_spin.setValue(3)
        adv_layout.addWidget(self._min_repeats_spin, 1, 1)

        adv_layout.addWidget(_lbl("Max expected repeats:"), 1, 2)
        self._max_repeats_spin = QSpinBox()
        self._max_repeats_spin.setRange(10, 200)
        self._max_repeats_spin.setValue(45)
        adv_layout.addWidget(self._max_repeats_spin, 1, 3)

        self._ld_filter_cb = QCheckBox("Apply LD filter (thin markers by distance)")
        self._ld_filter_cb.setChecked(False)
        self._ld_filter_cb.toggled.connect(self._toggle_ld_filter)
        adv_layout.addWidget(self._ld_filter_cb, 2, 0, 1, 2)

        adv_layout.addWidget(_lbl("Min distance (bp):"), 2, 2)
        self._ld_distance_spin = QSpinBox()
        self._ld_distance_spin.setRange(1000, 1000000)
        self._ld_distance_spin.setValue(10000)
        self._ld_distance_spin.setEnabled(False)
        adv_layout.addWidget(self._ld_distance_spin, 2, 3)

        adv_layout.setColumnStretch(4, 1)
        self._advanced_box.addWidget(adv_widget)
        mg.addWidget(self._advanced_box)

        L.addWidget(self._main_group)

        self._histogram_box = CollapsibleBox("Amplicon Size Distribution")
        self._histogram_widget = HistogramWidget()
        self._histogram_box.addWidget(self._histogram_widget)
        L.addWidget(self._histogram_box)

        self._tab_widget = QTabWidget()
        self._tab_widget.setVisible(False)

        results_tab = QWidget()
        rt_layout = QVBoxLayout(results_tab)
        self._summary_label = QLabel("")
        self._summary_label.setWordWrap(True)
        self._summary_label.setStyleSheet(f"font-weight: bold; color: {TEXT_PRIMARY};")
        rt_layout.addWidget(self._summary_label)

        cards_layout = QHBoxLayout()
        self._cards = {}
        for name, color in [("Total Loci", ACCENT), ("Passed", SUCCESS),
                            ("Dimer Risks", WARNING), ("Size Filtered", ERROR)]:
            card = self._make_card(name, "—", color)
            cards_layout.addWidget(card)
            self._cards[name] = card
        rt_layout.addLayout(cards_layout)

        self._issues_box = CollapsibleBox("Detailed Issues")
        issues_widget = QWidget()
        issues_layout = QVBoxLayout(issues_widget)
        self._dimer_table = self._create_issue_table("Dimer Risks")
        issues_layout.addWidget(self._dimer_table)
        self._size_filter_table = self._create_issue_table("Size Filtered Loci")
        issues_layout.addWidget(self._size_filter_table)
        self._gc_warning_table = self._create_issue_table("GC Content Warnings")
        issues_layout.addWidget(self._gc_warning_table)
        self._issues_box.addWidget(issues_widget)
        rt_layout.addWidget(self._issues_box)

        self._pool_table_label = QLabel("Clean Pool (Passed All Checks)")
        self._pool_table_label.setFont(QFont(FONT_UI, FONT_SIZE_NORMAL, QFont.Weight.Bold))
        self._pool_table_label.setStyleSheet(f"color: {ACCENT}; margin-top: 16px;")
        rt_layout.addWidget(self._pool_table_label)

        self._pool_table = QTableWidget()
        self._pool_table.setAlternatingRowColors(True)
        self._pool_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._pool_table.setSortingEnabled(True)
        self._pool_table.setFixedHeight(300)
        rt_layout.addWidget(self._pool_table)

        export_layout = QHBoxLayout()
        self._export_pool_btn = QPushButton("Export Pool (CSV)")
        self._export_pool_btn.clicked.connect(self._export_pool)
        export_layout.addWidget(self._export_pool_btn)
        self._export_fasta_btn = QPushButton("Export Amplicon Reference (FASTA)")
        self._export_fasta_btn.clicked.connect(self._export_amplicon_fasta)
        export_layout.addWidget(self._export_fasta_btn)
        self._export_gc_btn = QPushButton("Export GC Analysis (CSV)")
        self._export_gc_btn.clicked.connect(self._export_gc_analysis)
        export_layout.addWidget(self._export_gc_btn)
        export_layout.addStretch()
        rt_layout.addLayout(export_layout)

        self._tab_widget.addTab(results_tab, "Results")

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

        self._update_platform_limits()

    def _make_card(self, label, value, color):
        w = QWidget()
        w.setStyleSheet(f"background: {BG_MID}; border: 1px solid {BORDER}; border-radius: 6px;")
        l = QVBoxLayout(w)
        l.setContentsMargins(10, 6, 10, 6)
        val_lbl = QLabel(value)
        val_lbl.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
        val_lbl.setStyleSheet(f"color: {color}; background: transparent;")
        lbl = QLabel(label)
        lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL-1}pt; background: transparent;")
        l.addWidget(val_lbl)
        l.addWidget(lbl)
        w.val_lbl = val_lbl
        return w

    def _create_issue_table(self, title):
        container = QWidget()
        layout = QVBoxLayout(container)
        lbl = QLabel(title)
        lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        layout.addWidget(lbl)
        table = QTableWidget()
        table.setAlternatingRowColors(True)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setFixedHeight(120)
        layout.addWidget(table)
        return container

    def _update_card(self, name, value):
        if name in self._cards:
            self._cards[name].val_lbl.setText(str(value))

    def _update_platform_limits(self):
        platform_name = self._platform_combo.currentText()
        platform = next((p for p in PLATFORMS if p["name"] == platform_name), PLATFORMS[1])
        self._min_amp_spin.setValue(platform["min_amplicon"])
        self._max_amp_spin.setValue(platform["max_amplicon"])
        self._platform_hint.setText(f"(Suggested: {platform['min_amplicon']}–{platform['max_amplicon']} bp)")

    def _cancel_worker(self):
        if self._worker:
            self._worker.cancel()
            self._worker.quit()
            self._worker.wait()
            self._worker = None
        self._validate_btn.setEnabled(True)
        self._cancel_btn.setEnabled(False)
        self._status_label.setText("Cancelled")

    def _run_validation(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return
        if not self.state.genome:
            self._status_label.setText("Genome not loaded.")
            return

        self._validate_btn.setEnabled(False)
        self._cancel_btn.setEnabled(True)
        self._status_label.setText("Validating panel...")
        QApplication.processEvents()

        params = {
            "platform": self._platform_combo.currentText(),
            "forward_tag": self._fwd_tag_edit.text().upper(),
            "reverse_tag": self._rev_tag_edit.text().upper(),
            "include_tags_dimer": self._include_tags_dimer_cb.isChecked(),
            "min_amplicon": self._min_amp_spin.value(),
            "max_amplicon": self._max_amp_spin.value(),
            "max_3prime": self._max_3prime_spin.value(),
            "max_tm_spread": self._max_tm_spread_spin.value(),
            "min_repeats": self._min_repeats_spin.value(),
            "max_repeats": self._max_repeats_spin.value(),
            "auto_optimise": False,
        }

        self._worker = ValidateWorker(self._pass_primers, self.state.genome, params)
        self._worker.progress.connect(lambda msg: self._status_label.setText(msg))
        self._worker.finished.connect(self._on_validation_done)
        self._worker.error.connect(self._on_validation_error)
        self._worker.start()

    def _on_validation_done(self, result):
        self._validation_result = result
        self._validate_btn.setEnabled(True)
        self._cancel_btn.setEnabled(False)
        self._status_label.setText("Validation complete.")
        self._tab_widget.setVisible(True)

        total_unique = self._count_unique_loci(self._pass_primers)
        clean_pool = self._filter_unique_loci(result["clean_pool"])

        if self._ld_filter_cb.isChecked():
            from core.ld_filter import thin_markers_by_distance
            clean_pool = thin_markers_by_distance(clean_pool, self._ld_distance_spin.value())

        passed_unique = len(clean_pool)

        size_filtered_ids = set(result["size_result"]["filtered_loci"])
        size_filtered_unique = len({p["ssr_id"] for p in self._pass_primers if p["ssr_id"] in size_filtered_ids})

        dimer_flagged_ids = set()
        for d in result["multiplex_result"]["dimer_risks"]:
            dimer_flagged_ids.add(d["ssr_id_a"])
            dimer_flagged_ids.add(d["ssr_id_b"])
        dimer_unique = len({p["ssr_id"] for p in self._pass_primers if p["ssr_id"] in dimer_flagged_ids})

        self._update_card("Total Loci", str(total_unique))
        self._update_card("Passed", str(passed_unique))
        self._update_card("Dimer Risks", str(dimer_unique))
        self._update_card("Size Filtered", str(size_filtered_unique))

        self._summary_label.setText(f"✅ {passed_unique} of {total_unique} loci passed all checks.")
        self._populate_issue_tables(result)
        self._populate_pool_table(clean_pool)

        sizes = [p.get("product_size", 0) for p in clean_pool if p.get("product_size")]
        self._histogram_widget.set_data(sizes)

        self._all_markers = []
        for p in self._pass_primers:
            self._all_markers.append({"contig": p.get("contig"), "start": p.get("start", 0)})
        filtered_markers = [{"contig": m.get("contig"), "start": m.get("start")} for m in clean_pool]
        self._chromosome_view.set_data(self._all_markers, filtered_markers, self._contig_lengths)
        self._chromosome_view.set_show_filtered(self._ld_filter_cb.isChecked())

    def _on_validation_error(self, msg):
        self._validate_btn.setEnabled(True)
        self._cancel_btn.setEnabled(False)
        self._status_label.setText(f"Error: {msg}")

    def _populate_issue_tables(self, result):
        def _populate_table(table_widget, data, headers, row_mapper):
            table_widget.setRowCount(0)
            if not data:
                return
            table_widget.setRowCount(len(data))
            table_widget.setColumnCount(len(headers))
            table_widget.setHorizontalHeaderLabels(headers)
            for i, item in enumerate(data):
                row_data = row_mapper(item)
                for j, val in enumerate(row_data):
                    table_widget.setItem(i, j, QTableWidgetItem(str(val)))
            table_widget.resizeColumnsToContents()

        dimer_risks = result["multiplex_result"]["dimer_risks"]
        table = self._dimer_table.findChild(QTableWidget)
        _populate_table(
            table, dimer_risks,
            ["SSR A", "Dir A", "SSR B", "Dir B", "Score"],
            lambda r: (r["ssr_id_a"], r["dir_a"], r["ssr_id_b"], r["dir_b"], r["score"])
        )

        size_filtered = result["size_result"]["filtered_loci"]
        table = self._size_filter_table.findChild(QTableWidget)
        _populate_table(
            table, [{"ssr_id": sid} for sid in size_filtered],
            ["SSR ID"],
            lambda x: (x["ssr_id"],)
        )

        gc_warnings = []
        from core.amplicon_reference import extract_amplicon_robust, analyse_amplicon_gc
        for p in self._pass_primers:
            amplicon = extract_amplicon_robust(self.state.genome, p, search_flank=300)
            if amplicon:
                gc_info = analyse_amplicon_gc(amplicon)
                if gc_info["flag"] not in ("Good",):
                    gc_warnings.append({
                        "ssr_id": p["ssr_id"],
                        "gc_percent": gc_info["gc_percent"],
                        "flag": gc_info["flag"]
                    })
        table = self._gc_warning_table.findChild(QTableWidget)
        _populate_table(
            table, gc_warnings,
            ["SSR ID", "GC %", "Flag"],
            lambda x: (x["ssr_id"], x["gc_percent"], x["flag"])
        )

    def _count_unique_loci(self, pool: list) -> int:
        if not pool:
            return 0
        return len({p.get("ssr_id") for p in pool if p.get("ssr_id") is not None})

    def _filter_unique_loci(self, pool: list) -> list:
        if not pool:
            return []
        import pandas as pd
        df = pd.DataFrame(pool)
        if "pair_rank" not in df.columns:
            return pool
        df = df.sort_values("pair_rank").drop_duplicates(subset=["ssr_id"], keep="first")
        return df.to_dict(orient="records")

    def _populate_pool_table(self, pool):
        if not pool:
            self._pool_table.setRowCount(0)
            return
        import pandas as pd
        df = pd.DataFrame(pool)
        cols = {
            "ssr_id": "SSR ID",
            "pair_rank": "Rank",
            "contig": "Contig",
            "motif": "Motif",
            "product_size": "Size (bp)",
            "left_tm": "Fwd Tm",
            "right_tm": "Rev Tm",
            "left_primer": "Forward Primer",
            "right_primer": "Reverse Primer",
        }
        display_cols = [c for c in cols if c in df.columns]
        self._pool_table.setSortingEnabled(False)
        self._pool_table.setRowCount(len(df))
        self._pool_table.setColumnCount(len(display_cols))
        self._pool_table.setHorizontalHeaderLabels([cols[c] for c in display_cols])
        for i, row in df.iterrows():
            for j, col in enumerate(display_cols):
                val = row[col]
                text = f"{val:.2f}" if isinstance(val, float) else str(val)
                self._pool_table.setItem(i, j, QTableWidgetItem(text))
        self._pool_table.resizeColumnsToContents()
        self._pool_table.setSortingEnabled(True)

    def _export_pool(self):
        if not self._validation_result:
            return
        pool = self._validation_result.get("clean_pool", [])
        if not pool:
            self._status_label.setText("No loci to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Pool", "amplicon_pool.csv", "CSV (*.csv)")
        if not path:
            return
        import pandas as pd
        pd.DataFrame(pool).to_csv(path, index=False)
        self.mw.set_status(f"Pool exported to {path}")

    def _export_amplicon_fasta(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return
        if not self.state.genome:
            self._status_label.setText("Genome not loaded.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Amplicon Reference FASTA", "amplicon_reference.fasta", "FASTA files (*.fasta *.fa)")
        if not path:
            return
        try:
            from core.amplicon_reference import amplicons_to_fasta
            fasta_str, n_found, n_miss = amplicons_to_fasta(self._pass_primers, self.state.genome, search_flank=300)
            with open(path, 'w') as f:
                f.write(fasta_str)
            self._status_label.setText(f"Amplicon FASTA saved ({n_found} loci, {n_miss} skipped)")
            self.mw.set_status(f"Amplicon reference exported to {path}")
        except Exception as e:
            self._status_label.setText(f"Export error: {e}")

    def _export_gc_analysis(self):
        if not self._pass_primers:
            self._status_label.setText("No PASS primers available.")
            return
        if not self.state.genome:
            self._status_label.setText("Genome not loaded.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export GC Analysis", "gc_analysis.csv", "CSV (*.csv)")
        if not path:
            return
        try:
            from core.amplicon_reference import extract_amplicon_robust, analyse_amplicon_gc
            import pandas as pd
            rows = []
            for p in self._pass_primers:
                amplicon = extract_amplicon_robust(self.state.genome, p, search_flank=300)
                if amplicon:
                    gc_info = analyse_amplicon_gc(amplicon)
                    rows.append({
                        "ssr_id": p["ssr_id"],
                        "gc_percent": gc_info["gc_percent"],
                        "flag": gc_info["flag"]
                    })
            if rows:
                pd.DataFrame(rows).to_csv(path, index=False)
                self._status_label.setText(f"GC analysis exported ({len(rows)} loci)")
                self.mw.set_status(f"GC analysis exported to {path}")
            else:
                self._status_label.setText("No amplicons could be extracted for GC analysis.")
        except Exception as e:
            self._status_label.setText(f"Export error: {e}")

    def _get_pass_primers(self):
        spec = self.state.specificity_results or []
        primers = self.state.primer_results or []
        pass_ids = {r["ssr_id"] for r in spec if r.get("specificity_status") == "PASS"}
        return [p for p in primers if p.get("ssr_id") in pass_ids]

    def _rebuild(self):
        mode = self.state.workflow_mode
        is_amplicon = (mode == "amplicon")
        has_spec = self.state.has_specificity

        should_show = is_amplicon and has_spec
        self._placeholder.setVisible(not should_show)
        self._main_group.setVisible(should_show)
        self._tab_widget.setVisible(False)

        if should_show:
            self._pass_primers = self._get_pass_primers()
            self._status_label.setText(f"{self._count_unique_loci(self._pass_primers)} PASS loci ready.")

    def on_show(self):
        self._compute_contig_lengths()
        if self.state.blast_version != self._last_version:
            self._last_version = self.state.blast_version
            self._rebuild()