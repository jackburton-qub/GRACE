"""
ssr_summary_panel.py — SSR Summary & Explorer
Distribution charts, karyotype view, and statistics after SSR detection.

Chromosome names come from the GFF annotation only — we never guess
which contigs are chromosomes based on size or naming patterns.
If no GFF is loaded, all contigs are shown as contigs.
"""
import re
from collections import Counter, defaultdict

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QScrollArea, QSizePolicy,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import (
    QFont, QPainter, QColor, QPen, QBrush, QPainterPath, QFontMetrics,
)

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)

CHART_COLOURS = [
    "#4C9BE8", "#56C596", "#F4A261", "#E76F51",
    "#9B89C4", "#48CAE4", "#F6C90E", "#E07BA5",
    "#6BCB77", "#FF6B6B", "#4ECDC4", "#FFE66D",
]

FEATURE_COLOURS = {
    "exon":       "#56C596",
    "CDS":        "#4C9BE8",
    "intron":     "#F4A261",
    "intergenic": "#9B89C4",
    "unknown":    "#666688",
}

CHROM_COLOUR  = "#2A2A4A"
CHROM_BORDER  = "#4A4A7A"


# ---------------------------------------------------------------------------
# Karyotype widget
# ---------------------------------------------------------------------------

class KaryotypeWidget(QWidget):
    """
    Horizontal ideogram — one bar per sequence, SSR positions as ticks.
    Label shows chromosome name if known, else contig accession.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._chroms  = []   # (display_name, length, [(pos, feature), ...])
        self._max_len = 1
        self._title   = ""
        self.setMinimumHeight(60)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)

    def set_data(self, chroms: list, title: str = ""):
        self._chroms  = chroms
        self._max_len = max((c[1] for c in chroms), default=1)
        self._title   = title
        self.setMinimumHeight(max(200, len(chroms) * 46 + 60))
        self.update()

    def paintEvent(self, event):
        if not self._chroms:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        w, h    = self.width(), self.height()
        pad_l   = 110
        pad_r   = 100
        pad_t   = 28
        row_h   = 46
        bar_h   = 18
        chart_w = w - pad_l - pad_r
        bp_per_px = self._max_len / chart_w if chart_w > 0 else 1

        # Title
        if self._title:
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL - 1, QFont.Weight.Bold))
            painter.drawText(pad_l, pad_t - 8, self._title)

        for i, (name, length, positions) in enumerate(self._chroms):
            y_c   = pad_t + i * row_h + row_h // 2
            bar_y = y_c - bar_h // 2
            bar_w = max(4, int(length / bp_per_px))

            # Label
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.setFont(QFont(FONT_MONO, max(7, FONT_SIZE_SMALL - 1)))
            fm  = QFontMetrics(painter.font())
            lbl = name[:18]
            painter.drawText(pad_l - fm.horizontalAdvance(lbl) - 6,
                             y_c + fm.ascent() // 2, lbl)

            # Chromosome/contig bar
            path = QPainterPath()
            path.addRoundedRect(pad_l, bar_y, bar_w, bar_h,
                                bar_h // 2, bar_h // 2)
            painter.fillPath(path, QBrush(QColor(CHROM_COLOUR)))
            painter.strokePath(path, QPen(QColor(CHROM_BORDER), 1))

            # SSR ticks — binned to pixel resolution
            if positions:
                n_bins   = max(1, bar_w)
                bin_hits = defaultdict(list)
                for pos, feat in positions:
                    b = min(int(pos / length * n_bins), n_bins - 1)
                    bin_hits[b].append(feat or "unknown")
                for b, feats in bin_hits.items():
                    col = FEATURE_COLOURS.get(feats[0], "#4C9BE8")
                    x   = pad_l + b
                    painter.setPen(QPen(QColor(col), 1))
                    painter.drawLine(x, bar_y + 2, x, bar_y + bar_h - 2)

            # Right-side annotation
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.setFont(QFont(FONT_MONO, max(6, FONT_SIZE_SMALL - 2)))
            painter.drawText(
                pad_l + bar_w + 6, y_c + 4,
                f"{len(positions):,} SSRs   {_fmt_bp(length)}"
            )

        # Feature colour legend
        has_feats = any(
            feat and feat != "unknown"
            for _, _, positions in self._chroms
            for _, feat in positions
        )
        if has_feats:
            lx = pad_l
            ly = pad_t + len(self._chroms) * row_h + 10
            painter.setFont(QFont(FONT_UI, max(7, FONT_SIZE_SMALL - 1)))
            for feat, col in FEATURE_COLOURS.items():
                if feat == "unknown":
                    continue
                painter.setBrush(QBrush(QColor(col)))
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawRect(lx, ly, 10, 10)
                painter.setPen(QColor(TEXT_SECONDARY))
                painter.drawText(lx + 14, ly + 9, feat.capitalize())
                lx += 100

        painter.end()


# ---------------------------------------------------------------------------
# Bar chart widget
# ---------------------------------------------------------------------------

class BarChartWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data  = []
        self._title = ""
        self.setMinimumHeight(200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data: list, title: str = ""):
        self._title = title
        self._data  = []
        for i, item in enumerate(data):
            if len(item) == 3:
                self._data.append(item)
            else:
                label, value = item
                self._data.append((label, value, CHART_COLOURS[i % len(CHART_COLOURS)]))
        self.update()

    def paintEvent(self, event):
        if not self._data:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        pl, pr, pt, pb = 55, 12, 28, 52
        cw = w - pl - pr; ch = h - pt - pb
        max_val = max(v for _, v, _ in self._data) or 1
        n = len(self._data)
        gap = max(2, cw // (n * 5)); bw = max(4, (cw - gap * (n + 1)) // n)

        if self._title:
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
            painter.drawText(pl, pt - 8, self._title)

        painter.setPen(QPen(QColor(BORDER), 1))
        painter.drawLine(pl, pt, pl, pt + ch)
        painter.drawLine(pl, pt + ch, pl + cw, pt + ch)
        painter.setFont(QFont(FONT_MONO, max(6, FONT_SIZE_SMALL - 2)))
        for frac in [0, 0.25, 0.5, 0.75, 1.0]:
            yv = int(max_val * frac); yp = pt + ch - int(ch * frac)
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(0, yp + 4, pl - 4, 12, Qt.AlignmentFlag.AlignRight, _fmt_num(yv))
            if frac > 0:
                painter.setPen(QPen(QColor(BORDER), 1, Qt.PenStyle.DotLine))
                painter.drawLine(pl, yp, pl + cw, yp)
        for i, (label, value, colour) in enumerate(self._data):
            x = pl + gap + i * (bw + gap); bh = int(ch * value / max_val); y = pt + ch - bh
            path = QPainterPath(); path.addRoundedRect(x, y, bw, bh, 3, 3)
            painter.fillPath(path, QBrush(QColor(colour)))
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.save()
            painter.translate(x + bw // 2, pt + ch + 6); painter.rotate(40)
            painter.drawText(0, 0, str(label)[:20]); painter.restore()
        painter.end()


# ---------------------------------------------------------------------------
# Pie chart widget
# ---------------------------------------------------------------------------

class PieChartWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data  = []
        self.setMinimumHeight(200); self.setMinimumWidth(280)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data: list):
        self._data = data; self.update()

    def paintEvent(self, event):
        if not self._data: return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        pad = 16; lw = 160
        pw  = min(w - lw - pad * 3, h - pad * 2)
        px  = pad; py = pad + (h - pw - pad * 2) // 2
        total = sum(v for _, v, _ in self._data) or 1
        angle = 90 * 16
        for _, value, colour in self._data:
            span = int(360 * 16 * value / total)
            painter.setBrush(QBrush(QColor(colour)))
            painter.setPen(QPen(QColor(BG_MID), 2))
            painter.drawPie(px, py, pw, pw, angle, span)
            angle += span
        lx = px + pw + pad; ly = py
        painter.setFont(QFont(FONT_UI, max(7, FONT_SIZE_SMALL - 1)))
        for label, value, colour in self._data:
            pct = value / total * 100
            painter.setBrush(QBrush(QColor(colour))); painter.setPen(Qt.PenStyle.NoPen)
            painter.drawRoundedRect(lx, ly, 10, 10, 2, 2)
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.drawText(lx + 14, ly + 9, f"{label} ({pct:.1f}%)")
            ly += 20
        painter.end()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_num(n):
    if n >= 1_000_000: return f"{n/1_000_000:.1f}M"
    if n >= 1_000:     return f"{n/1_000:.1f}k"
    return str(n)


def _fmt_bp(n):
    if n >= 1_000_000: return f"{n/1_000_000:.1f}Mb"
    if n >= 1_000:     return f"{n/1_000:.0f}kb"
    return f"{n}bp"


def _build_karyotype_data(ssrs, genome, chrom_names):
    """
    Build per-sequence SSR position data.
    display_name = chromosome name if known from GFF, else contig accession.
    Sorted by sequence length descending.
    """
    by_contig = defaultdict(list)
    for s in ssrs:
        feat = s.get("genomic_feature")
        by_contig[s["contig"]].append((s["start"], feat))

    result = []
    for contig, seq in genome.items():
        display = chrom_names.get(contig, contig) if chrom_names else contig
        result.append((display, len(seq), by_contig.get(contig, [])))

    result.sort(key=lambda x: -x[1])
    return result


def _build_contig_bar_data(ssrs, genome, chrom_names, top_n=30):
    """
    Build bar chart data for contig SSR counts.
    Uses display names where available.
    Shows top_n by count plus 'Other'.
    """
    counts = Counter(s["contig"] for s in ssrs)
    items  = sorted(counts.items(), key=lambda x: -x[1])
    top    = items[:top_n]
    other  = sum(n for _, n in items[top_n:])

    result = []
    for contig, n in top:
        display = chrom_names.get(contig, contig) if chrom_names else contig
        # Truncate long accession numbers for display
        if len(display) > 12 and "." in display:
            display = display[:8] + "…"
        result.append((display, n, CHART_COLOURS[len(result) % len(CHART_COLOURS)]))

    if other:
        result.append((f"Other ({len(items)-top_n})", other, "#666688"))

    return result


def _build_motif_distribution(ssrs):
    labels = {2: "Di", 3: "Tri", 4: "Tetra", 5: "Penta", 6: "Hexa"}
    counts = Counter(len(s["motif"]) for s in ssrs)
    return [(labels.get(k, f"{k}bp"), counts[k], CHART_COLOURS[i % len(CHART_COLOURS)])
        for i, k in enumerate(sorted(counts))]


def _build_feature_distribution(ssrs):
    if not ssrs or "genomic_feature" not in ssrs[0]:
        return []
    counts = Counter(s.get("genomic_feature", "unknown") for s in ssrs)
    order  = ["exon", "CDS", "intron", "intergenic", "unknown"]
    return [(f.capitalize(), counts[f], FEATURE_COLOURS.get(f, "#888888"))
            for f in order if f in counts]


def _build_top_motifs(ssrs, n=10):
    counts = Counter(s.get("canonical_motif", s["motif"]) for s in ssrs)
    return [(m, c, CHART_COLOURS[i % len(CHART_COLOURS)])
            for i, (m, c) in enumerate(counts.most_common(n))]


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class SSRSummaryPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._last_version = -1
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
        self._L = QVBoxLayout(content)
        self._L.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        self._L.setSpacing(16)

        # Title
        title = QLabel("SSR Summary")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Genomic distribution and statistics of detected SSRs")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._L.addWidget(title); self._L.addWidget(sub)

        self._placeholder = QLabel("Run SSR Detection first to see the summary.")
        self._placeholder.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._L.addWidget(self._placeholder)

        # Overview stats
        self._stats_group = QGroupBox("Overview")
        self._stats_group.setVisible(False)
        self._stats_layout = QGridLayout(self._stats_group)
        self._stats_layout.setSpacing(12)
        self._L.addWidget(self._stats_group)

        # GFF annotation status note
        self._gff_note = QLabel("")
        self._gff_note.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._gff_note.setWordWrap(True)
        self._gff_note.setVisible(False)
        self._L.addWidget(self._gff_note)

        # Karyotype / positional view
        self._karyo_group = QGroupBox("")
        self._karyo_group.setVisible(False)
        kg = QVBoxLayout(self._karyo_group)
        self._karyo_desc = QLabel("")
        self._karyo_desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._karyo_desc.setWordWrap(True)
        kg.addWidget(self._karyo_desc)
        self._karyo_widget = KaryotypeWidget()
        kg.addWidget(self._karyo_widget)
        self._L.addWidget(self._karyo_group)

        # Top contigs bar chart (when too many sequences to show individually)
        self._contig_group = QGroupBox("SSR Count by Sequence (Top 30)")
        self._contig_group.setVisible(False)
        cg = QVBoxLayout(self._contig_group)
        self._contig_chart = BarChartWidget()
        cg.addWidget(self._contig_chart)
        self._L.addWidget(self._contig_group)

        # Motif type + feature pies
        pie_row = QHBoxLayout()

        self._motif_group = QGroupBox("Motif Type Distribution")
        self._motif_group.setVisible(False)
        mg = QVBoxLayout(self._motif_group)
        self._motif_chart = PieChartWidget()
        mg.addWidget(self._motif_chart)
        pie_row.addWidget(self._motif_group)

        self._feature_group = QGroupBox("Genomic Feature Distribution")
        self._feature_group.setVisible(False)
        fg = QVBoxLayout(self._feature_group)
        self._feature_note = QLabel(
            "Load a GFF annotation on the Home page to enable genomic feature classification."
        )
        self._feature_note.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        self._feature_note.setWordWrap(True)
        self._feature_chart = PieChartWidget()
        self._feature_chart.setVisible(False)
        fg.addWidget(self._feature_note); fg.addWidget(self._feature_chart)
        pie_row.addWidget(self._feature_group)

        pie_w = QWidget(); pie_w.setLayout(pie_row)
        self._L.addWidget(pie_w)

        # Top motifs
        self._topmotif_group = QGroupBox("Top 10 Canonical Motifs")
        self._topmotif_group.setVisible(False)
        tm = QVBoxLayout(self._topmotif_group)
        self._topmotif_chart = BarChartWidget()
        tm.addWidget(self._topmotif_chart)
        self._L.addWidget(self._topmotif_group)

        # Navigate
        nav_row = QHBoxLayout()
        self._design_btn = QPushButton("Design primers for all SSRs →")
        self._design_btn.setObjectName("primary")
        self._design_btn.setVisible(False)
        self._design_btn.clicked.connect(lambda: self.mw.navigate_to(3))
        nav_row.addWidget(self._design_btn); nav_row.addStretch()
        self._L.addLayout(nav_row)
        self._L.addStretch()

    def _rebuild(self):
        if not self.state.has_ssrs:
            self._placeholder.setVisible(True)
            for g in [self._stats_group, self._karyo_group, self._contig_group,
                      self._motif_group, self._feature_group, self._topmotif_group]:
                g.setVisible(False)
            self._gff_note.setVisible(False)
            self._design_btn.setVisible(False)
            return

        self._placeholder.setVisible(False)
        ssrs        = self.state.ssrs
        genome      = self.state.genome or {}
        chrom_names = self.state.chrom_names or {}
        has_gff     = self.state.has_gff
        n_chrom_mapped = len(chrom_names)

        # ── Overview stats ────────────────────────────────
        while self._stats_layout.count():
            item = self._stats_layout.takeAt(0)
            if item.widget(): item.widget().deleteLater()

        total_bp    = sum(len(s) for s in genome.values()) if genome else 0
        density     = len(ssrs) / (total_bp / 1_000_000) if total_bp else 0
        mean_repeat = sum(s["repeat_count"] for s in ssrs) / len(ssrs)
        n_contigs   = len(genome)

        stats = [
            ("Total SSRs",        f"{len(ssrs):,}"),
            ("SSR density",       f"{density:.1f} / Mb"),
            ("Unique motifs",     f"{len(set(s['motif'] for s in ssrs)):,}"),
            ("Mean repeat count", f"{mean_repeat:.1f}"),
            ("Sequences",         f"{n_contigs:,}"),
            ("With SSRs",         f"{len(set(s['contig'] for s in ssrs)):,}"),
            ("Annotation",        "Yes ✓" if "genomic_feature" in ssrs[0] else "None loaded"),
            ("Chr names mapped",  f"{n_chrom_mapped:,}" if n_chrom_mapped else "None (load GFF)"),
        ]
        for i, (label, value) in enumerate(stats):
            lbl = QLabel(label)
            lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
            val = QLabel(value)
            val.setFont(QFont(FONT_MONO, FONT_SIZE_NORMAL, QFont.Weight.Bold))
            val.setStyleSheet(f"color: {TEXT_PRIMARY};")
            row = i // 4; col = (i % 4) * 2
            self._stats_layout.addWidget(lbl, row, col)
            self._stats_layout.addWidget(val, row, col + 1)
        self._stats_group.setVisible(True)

        # ── GFF note ──────────────────────────────────────
        if not has_gff:
            self._gff_note.setText(
                "💡 Tip: Load a GFF annotation file on the Home page to enable chromosome name "
                "labels, genomic feature classification (exon/intron/intergenic), and richer "
                "visualisations throughout GRACE."
            )
            self._gff_note.setVisible(True)
        elif n_chrom_mapped > 0:
            self._gff_note.setText(
                f"✓ GFF loaded — {n_chrom_mapped:,} chromosome name{'s' if n_chrom_mapped!=1 else ''} "
                f"mapped from annotation (e.g. {next(iter(chrom_names.values()))})."
            )
            self._gff_note.setVisible(True)
        else:
            self._gff_note.setText(
                "GFF loaded but no chromosome name mappings found in region features. "
                "Sequences are shown by their accession numbers."
            )
            self._gff_note.setVisible(True)

        # ── Positional view ───────────────────────────────
        n_seqs = len(genome)

        if n_seqs <= 60:
            # Show full karyotype for all sequences
            karyo_data = _build_karyotype_data(ssrs, genome, chrom_names)
            if n_chrom_mapped > 0:
                title = "SSR Positions Along Each Chromosome"
                desc  = (
                    f"Each bar represents a sequence scaled by length. "
                    f"Chromosome names from GFF annotation. "
                    f"Vertical ticks show exact SSR positions."
                    + (" Tick colour = genomic feature." if "genomic_feature" in ssrs[0] else "")
                )
            else:
                title = "SSR Positions Along Each Sequence"
                desc  = (
                    f"Each bar represents a contig/scaffold scaled by length. "
                    f"Vertical ticks show exact SSR positions. "
                    f"Load a GFF to show chromosome names instead of accession numbers."
                )
            self._karyo_group.setTitle(title)
            self._karyo_desc.setText(desc)
            self._karyo_widget.set_data(karyo_data)
            self._karyo_group.setVisible(True)
            self._contig_group.setVisible(False)
        else:
            # Too many sequences for individual bars — show top-N bar chart
            self._karyo_group.setVisible(False)
            contig_data = _build_contig_bar_data(ssrs, genome, chrom_names)
            self._contig_chart.set_data(contig_data)
            self._contig_group.setTitle(
                f"SSR Count by Sequence — Top 30 of {n_seqs:,} total sequences"
            )
            self._contig_group.setVisible(True)

        # ── Motif type pie ────────────────────────────────
        self._motif_chart.set_data(_build_motif_distribution(ssrs))
        self._motif_group.setVisible(True)

        # ── Feature pie ───────────────────────────────────
        feature_data = _build_feature_distribution(ssrs)
        self._feature_group.setVisible(True)
        if feature_data:
            self._feature_note.setVisible(False)
            self._feature_chart.setVisible(True)
            self._feature_chart.set_data(feature_data)
        else:
            self._feature_note.setVisible(True)
            self._feature_chart.setVisible(False)

        # ── Top motifs ────────────────────────────────────
        self._topmotif_chart.set_data(_build_top_motifs(ssrs))
        self._topmotif_group.setVisible(True)

        self._design_btn.setVisible(True)

    def on_show(self):
        if self.state.ssr_version != self._last_version:
            self._last_version = self.state.ssr_version
            self._rebuild()