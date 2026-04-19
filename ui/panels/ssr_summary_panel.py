"""
ssr_summary_panel.py — SSR Summary Panel
With fallback chromosome name assignment.
"""

import os
from collections import Counter
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QTabWidget, QSizePolicy,
    QGridLayout,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, SUCCESS, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BORDER,
)

CHART_COLORS = [
    "#4C9BE8", "#56C596", "#F4A261", "#E76F51",
    "#9B89C4", "#48CAE4", "#F6C90E", "#E07BA5",
]


class PieChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data = {}
        self.setMinimumHeight(200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data):
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        size = min(w, h) - 40
        cx, cy = w // 2, h // 2
        total = sum(self._data.values())
        if total == 0:
            return
        start_angle = 0
        legend_x = cx + size // 2 + 20
        legend_y = cy - (len(self._data) * 20) // 2
        for i, (label, count) in enumerate(self._data.items()):
            angle = int(360 * count / total)
            color = QColor(CHART_COLORS[i % len(CHART_COLORS)])
            painter.setBrush(QBrush(color))
            painter.setPen(QPen(Qt.GlobalColor.white, 2))
            painter.drawPie(cx - size//2, cy - size//2, size, size, start_angle * 16, angle * 16)
            painter.setPen(QPen(QColor(TEXT_PRIMARY), 1))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL))
            pct = (count / total) * 100
            painter.drawText(legend_x, legend_y + i * 20, f"{label}: {count:,} ({pct:.1f}%)")
            start_angle += angle


class VerticalBarChart(QWidget):
    """Vertical bar chart with X‑axis labels and Y‑axis counts."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data = []  # list of (label, count)
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data):
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        margin_left, margin_right, margin_top, margin_bottom = 50, 30, 30, 50
        chart_w = w - margin_left - margin_right
        chart_h = h - margin_top - margin_bottom

        max_count = max(c for _, c in self._data) if self._data else 1
        n = len(self._data)
        if n == 0:
            return

        bar_width = max(10, chart_w // n - 4)

        # Draw axes
        painter.setPen(QPen(QColor(TEXT_SECONDARY), 1))
        painter.drawLine(margin_left, h - margin_bottom, w - margin_right, h - margin_bottom)
        painter.drawLine(margin_left, margin_top, margin_left, h - margin_bottom)

        # Y‑axis labels (0, max/2, max)
        painter.setFont(QFont(FONT_MONO, 8))
        painter.drawText(margin_left - 25, h - margin_bottom + 4, "0")
        mid_y = margin_top + chart_h // 2
        painter.drawText(margin_left - 40, mid_y + 4, f"{max_count//2}")
        painter.drawText(margin_left - 40, margin_top + 10, f"{max_count}")

        # Draw bars
        for i, (label, count) in enumerate(self._data):
            x = margin_left + i * (bar_width + 4) + 4
            bar_h = int((count / max_count) * chart_h)
            y = h - margin_bottom - bar_h
            color = QColor(CHART_COLORS[i % len(CHART_COLORS)])
            painter.setBrush(QBrush(color))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawRect(x, y, bar_width, bar_h)

            # X‑axis label
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.setFont(QFont(FONT_UI, 8))
            display_label = str(label)[:12] + ("…" if len(str(label)) > 12 else "")
            painter.drawText(x, h - margin_bottom + 18, display_label)


class SSRSummaryPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw = main_window
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

        title = QLabel("Genomic distribution and statistics of detected SSRs")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        L.addWidget(title)

        self._tab_widget = QTabWidget()
        L.addWidget(self._tab_widget)

        # Statistics Tab
        stats_tab = QWidget()
        self._build_stats_tab(stats_tab)
        self._tab_widget.addTab(stats_tab, "Statistics")

        # Genomic Distribution Tab
        dist_tab = QWidget()
        self._build_distribution_tab(dist_tab)
        self._tab_widget.addTab(dist_tab, "Genomic Distribution")

    def _build_stats_tab(self, tab):
        layout = QVBoxLayout(tab)
        layout.setSpacing(16)

        # Overview cards
        self._overview_grid = QGridLayout()
        layout.addLayout(self._overview_grid)

        # SSR Count by Sequence (vertical bars)
        layout.addWidget(QLabel("SSR Count by Sequence — Top 15 of total sequences"))
        self._contig_bar = VerticalBarChart()
        layout.addWidget(self._contig_bar)

        # Pie charts row (side‑by‑side)
        pie_row = QHBoxLayout()
        pie_row.setSpacing(20)

        # Motif type pie chart (left)
        motif_widget = QWidget()
        motif_layout = QVBoxLayout(motif_widget)
        motif_layout.addWidget(QLabel("Motif Type Distribution"))
        self._motif_pie = PieChart()
        motif_layout.addWidget(self._motif_pie)
        pie_row.addWidget(motif_widget)

        # Genomic feature pie chart (right)
        feature_widget = QWidget()
        feature_layout = QVBoxLayout(feature_widget)
        feature_layout.addWidget(QLabel("Genomic Feature Distribution"))
        self._feature_pie = PieChart()
        feature_layout.addWidget(self._feature_pie)
        pie_row.addWidget(feature_widget)

        layout.addLayout(pie_row)

        # Store reference to feature widget to hide/show later
        self._feature_widget = feature_widget

        # Top motifs vertical bar chart
        layout.addWidget(QLabel("Top 15 Canonical Motifs"))
        self._top_motifs_bar = VerticalBarChart()
        layout.addWidget(self._top_motifs_bar)

    def _build_distribution_tab(self, tab):
        layout = QVBoxLayout(tab)
        layout.setSpacing(16)

        self._warning_banner = QLabel(
            "Genome is fragmented. Showing contigs in arbitrary order. "
            "Load a GFF with chromosome annotations for proper grouping."
        )
        self._warning_banner.setStyleSheet(f"background: {WARNING}22; color: {WARNING}; padding: 8px; border-radius: 4px;")
        self._warning_banner.setWordWrap(True)
        self._warning_banner.setVisible(False)
        layout.addWidget(self._warning_banner)

        layout.addWidget(QLabel("SSR Count by Chromosome / Contig (Top 15)"))
        self._chrom_bar = VerticalBarChart()
        layout.addWidget(self._chrom_bar)

        layout.addWidget(QLabel("All Chromosomes / Contigs"))
        self._dist_table = QTableWidget()
        self._dist_table.setAlternatingRowColors(True)
        self._dist_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._dist_table.setSortingEnabled(True)
        layout.addWidget(self._dist_table)

    def _refresh(self):
        if not self.state.has_ssrs:
            return

        ssrs = self.state.ssrs
        genome = self.state.genome
        has_gff = self.state.gff_features is not None

        # --- FALLBACK: Ensure chrom_name is populated ---
        gff_index = self.state.gff_features
        if gff_index and hasattr(gff_index, 'chrom_names'):
            chrom_names = gff_index.chrom_names
            for ssr in ssrs:
                if "chrom_name" not in ssr:
                    contig = ssr.get("contig", "")
                    ssr["chrom_name"] = chrom_names.get(contig, contig)

        # --- Diagnostic print (remove after testing) ---
        if ssrs:
            sample = ssrs[0]
            print("Sample SSR dict keys:", list(sample.keys()))
            print("chrom_name present?", "chrom_name" in sample)
            print("chrom_name value:", sample.get("chrom_name"))
            print("contig:", sample.get("contig"))

        total_ssrs = len(ssrs)
        genome_size = sum(len(s) for s in genome.values()) if genome else 0
        density = (total_ssrs / (genome_size / 1_000_000)) if genome_size > 0 else 0

        motifs = [ssr.get("canonical_motif", ssr.get("motif", "")) for ssr in ssrs]
        motif_counts = Counter(motifs)
        unique_motifs = len(motif_counts)

        repeat_counts = [ssr.get("repeat_count", 0) for ssr in ssrs]
        mean_repeats = sum(repeat_counts) / total_ssrs if total_ssrs else 0

        # Overview cards
        for i in reversed(range(self._overview_grid.count())):
            w = self._overview_grid.itemAt(i).widget()
            if w:
                w.deleteLater()

        cards = [
            ("Total SSRs", f"{total_ssrs:,}"),
            ("SSR density", f"{density:.2f} / Mb"),
            ("Unique motifs", f"{unique_motifs:,}"),
            ("Mean repeat count", f"{mean_repeats:.1f}"),
        ]
        for col, (label, value) in enumerate(cards):
            card = QWidget()
            card.setStyleSheet(f"background: {BG_MID}; border: 1px solid {BORDER}; border-radius: 6px; padding: 10px;")
            cl = QVBoxLayout(card)
            vl = QLabel(value)
            vl.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
            vl.setStyleSheet(f"color: {ACCENT};")
            cl.addWidget(vl)
            cl.addWidget(QLabel(label))
            self._overview_grid.addWidget(card, 0, col)

        # SSR Count by Sequence
        contig_counts = Counter(ssr.get("contig", "?") for ssr in ssrs)
        top_contigs = contig_counts.most_common(15)
        self._contig_bar.set_data(top_contigs)

        # Motif type distribution
        motif_types = {"Di": 0, "Tri": 0, "Tetra": 0, "Penta": 0, "Hexa": 0}
        for m in motifs:
            length = len(m)
            if length == 2:
                motif_types["Di"] += 1
            elif length == 3:
                motif_types["Tri"] += 1
            elif length == 4:
                motif_types["Tetra"] += 1
            elif length == 5:
                motif_types["Penta"] += 1
            elif length >= 6:
                motif_types["Hexa"] += 1
        self._motif_pie.set_data(motif_types)

        # Remove any motif type that has zero count (not scanned)
        motif_types = {k: v for k, v in motif_types.items() if v > 0}
        self._motif_pie.set_data(motif_types)

        # Genomic feature distribution (only if GFF loaded)
        if has_gff:
            features = [ssr.get("genomic_feature", "intergenic") for ssr in ssrs]
            feature_counts = Counter(features)
            feature_data = {
                "Exon": feature_counts.get("exon", 0),
                "Intron": feature_counts.get("intron", 0),
                "Intergenic": feature_counts.get("intergenic", 0),
            }
            self._feature_pie.set_data(feature_data)
            self._feature_widget.setVisible(True)
        else:
            self._feature_widget.setVisible(False)

        # Top motifs
        top_motifs = motif_counts.most_common(15)
        self._top_motifs_bar.set_data(top_motifs)

        # Genomic Distribution Tab
        counts = {}
        for ssr in ssrs:
            display = ssr.get("chrom_name", ssr.get("contig", "?"))
            counts[display] = counts.get(display, 0) + 1

        has_chrom_names = any(ssr.get("chrom_name") for ssr in ssrs)
        many_contigs = len(counts) > 50
        self._warning_banner.setVisible(not has_chrom_names and many_contigs)

        top_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:15]
        self._chrom_bar.set_data(top_items)

        all_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        self._dist_table.setRowCount(len(all_items))
        self._dist_table.setColumnCount(2)
        self._dist_table.setHorizontalHeaderLabels(["Chromosome / Contig", "SSR Count"])
        for i, (name, cnt) in enumerate(all_items):
            self._dist_table.setItem(i, 0, QTableWidgetItem(str(name)))
            self._dist_table.setItem(i, 1, QTableWidgetItem(str(cnt)))
        self._dist_table.resizeColumnsToContents()

    def on_show(self):
        self._refresh()