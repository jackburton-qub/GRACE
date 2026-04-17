"""
report_panel.py — Final Report
Full pipeline summary with charts, statistics, and PDF export.
Shows the complete journey from raw SSRs to validated primer pairs.
"""
import os, sys
from collections import Counter, defaultdict

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QScrollArea, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QFileDialog,
    QFrame, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont, QPainter, QColor, QPen, QBrush, QPainterPath

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)

TABLE_DISPLAY_LIMIT = 10_000

CHART_COLOURS = [
    "#4C9BE8", "#56C596", "#F4A261", "#E76F51",
    "#9B89C4", "#48CAE4", "#F6C90E", "#E07BA5",
]

FEATURE_COLOURS = {
    "exon":       "#56C596",
    "CDS":        "#4C9BE8",
    "intron":     "#F4A261",
    "intergenic": "#9B89C4",
    "unknown":    "#666688",
}


# ---------------------------------------------------------------------------
# Mini chart widgets (self-contained, no external deps)
# ---------------------------------------------------------------------------

class MiniPieChart(QWidget):
    def __init__(self, title="", parent=None):
        super().__init__(parent)
        self._data  = []
        self._title = title
        self.setMinimumHeight(180)
        self.setMinimumWidth(260)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data):
        """data: list of (label, value, colour)"""
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data: return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        pad  = 12
        lw   = 140
        pw   = min(w - lw - pad * 3, h - pad * 2 - 20)
        px   = pad
        py   = pad + 18 + max(0, (h - pw - pad * 2 - 18) // 2)

        if self._title:
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
            painter.drawText(px, pad + 12, self._title)

        total = sum(v for _, v, _ in self._data) or 1
        angle = 90 * 16
        for _, value, colour in self._data:
            span = int(360 * 16 * value / total)
            painter.setBrush(QBrush(QColor(colour)))
            painter.setPen(QPen(QColor(BG_MID), 2))
            painter.drawPie(px, py, pw, pw, angle, span)
            angle += span

        lx = px + pw + pad
        ly = py
        painter.setFont(QFont(FONT_UI, max(7, FONT_SIZE_SMALL - 1)))
        for label, value, colour in self._data:
            pct = value / total * 100
            painter.setBrush(QBrush(QColor(colour)))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawRoundedRect(lx, ly, 9, 9, 2, 2)
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.drawText(lx + 13, ly + 8, f"{label}: {value:,} ({pct:.1f}%)")
            ly += 19
        painter.end()


class MiniBarChart(QWidget):
    def __init__(self, title="", parent=None):
        super().__init__(parent)
        self._data  = []
        self._title = title
        self.setMinimumHeight(180)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_data(self, data):
        """data: list of (label, value, colour)"""
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data: return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        pl, pr, pt, pb = 50, 12, 28, 48
        cw = w - pl - pr
        ch = h - pt - pb
        max_val = max(v for _, v, _ in self._data) or 1
        n   = len(self._data)
        gap = max(2, cw // (n * 5))
        bw  = max(4, (cw - gap * (n + 1)) // n)

        if self._title:
            painter.setPen(QColor(TEXT_PRIMARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
            painter.drawText(pl, pt - 8, self._title)

        painter.setPen(QPen(QColor(BORDER), 1))
        painter.drawLine(pl, pt, pl, pt + ch)
        painter.drawLine(pl, pt + ch, pl + cw, pt + ch)

        painter.setFont(QFont(FONT_MONO, max(6, FONT_SIZE_SMALL - 2)))
        for frac in [0, 0.5, 1.0]:
            yv = int(max_val * frac)
            yp = pt + ch - int(ch * frac)
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.drawText(0, yp + 4, pl - 4, 12, Qt.AlignmentFlag.AlignRight, _fmt_num(yv))
            if frac > 0:
                painter.setPen(QPen(QColor(BORDER), 1, Qt.PenStyle.DotLine))
                painter.drawLine(pl, yp, pl + cw, yp)

        for i, (label, value, colour) in enumerate(self._data):
            x  = pl + gap + i * (bw + gap)
            bh = int(ch * value / max_val)
            y  = pt + ch - bh
            path = QPainterPath()
            path.addRoundedRect(x, y, bw, bh, 3, 3)
            painter.fillPath(path, QBrush(QColor(colour)))
            painter.setPen(QColor(TEXT_SECONDARY))
            painter.setFont(QFont(FONT_MONO, max(6, FONT_SIZE_SMALL - 2)))
            painter.save()
            painter.translate(x + bw // 2, pt + ch + 6)
            painter.rotate(35)
            painter.drawText(0, 0, str(label)[:18])
            painter.restore()
        painter.end()


class FunnelWidget(QWidget):
    """
    Visualises the filtering pipeline as a funnel:
    SSRs detected → Primers designed → Primers filtered → BLAST PASS
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self._stages = []
        self.setMinimumHeight(120)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_stages(self, stages):
        """stages: list of (label, count, colour)"""
        self._stages = stages
        self.update()

    def paintEvent(self, event):
        if not self._stages: return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h  = self.width(), self.height()
        n     = len(self._stages)
        pad   = 16
        bh    = (h - pad * 2) // n - 6
        max_v = max(s[1] for s in self._stages) or 1

        for i, (label, count, colour) in enumerate(self._stages):
            frac  = count / max_v
            bw    = max(60, int((w - pad * 2) * frac))
            x     = pad + (w - pad * 2 - bw) // 2
            y     = pad + i * (bh + 6)

            path = QPainterPath()
            path.addRoundedRect(x, y, bw, bh, 4, 4)
            painter.fillPath(path, QBrush(QColor(colour + "CC")))
            painter.strokePath(path, QPen(QColor(colour), 1))

            painter.setPen(QColor(TEXT_PRIMARY))
            painter.setFont(QFont(FONT_UI, FONT_SIZE_SMALL - 1, QFont.Weight.Bold))
            painter.drawText(x + 8, y + bh // 2 + 5, f"{label}: {count:,}")

            if i < n - 1:
                pct = self._stages[i + 1][1] / count * 100 if count else 0
                painter.setPen(QColor(TEXT_SECONDARY))
                painter.setFont(QFont(FONT_MONO, max(6, FONT_SIZE_SMALL - 2)))
                painter.drawText(w - 60, y + bh + 2, f"↓ {pct:.1f}%")
        painter.end()


# ---------------------------------------------------------------------------
# PDF worker
# ---------------------------------------------------------------------------

class PDFWorker(QThread):
    finished = pyqtSignal(str)
    error    = pyqtSignal(str)

    def __init__(self, pass_primers, state, output_path):
        super().__init__()
        self.pass_primers = pass_primers
        self.state        = state
        self.output_path  = output_path

    def run(self):
        try:
            _write_pdf(self.pass_primers, self.state, self.output_path)
            self.finished.emit(self.output_path)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_num(n):
    if n >= 1_000_000: return f"{n/1_000_000:.1f}M"
    if n >= 1_000:     return f"{n/1_000:.1f}k"
    return str(n)


def _stat_card(label, value, colour=None):
    w = QWidget()
    w.setStyleSheet(f"background: {BG_MID}; border: 1px solid {BORDER}; border-radius: 6px;")
    l = QVBoxLayout(w); l.setContentsMargins(12, 8, 12, 8); l.setSpacing(2)
    lbl = QLabel(label)
    lbl.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL - 1}pt; background: transparent; border: none;")
    val = QLabel(value)
    val.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
    val.setStyleSheet(f"color: {colour or ACCENT}; background: transparent; border: none;")
    l.addWidget(val); l.addWidget(lbl)
    return w


def _section_divider(title):
    w = QWidget()
    l = QHBoxLayout(w); l.setContentsMargins(0, 8, 0, 4); l.setSpacing(8)
    lbl = QLabel(title)
    lbl.setFont(QFont(FONT_UI, FONT_SIZE_NORMAL, QFont.Weight.Bold))
    lbl.setStyleSheet(f"color: {ACCENT};")
    line = QFrame(); line.setFrameShape(QFrame.Shape.HLine)
    line.setStyleSheet(f"color: {BORDER};")
    l.addWidget(lbl); l.addWidget(line, 1)
    return w


# ---------------------------------------------------------------------------
# PDF generation
# ---------------------------------------------------------------------------

def _write_pdf(pass_primers, state, path):
    try:
        from reportlab.lib.pagesizes import A4
        _write_pdf_reportlab(pass_primers, state, path)
    except ImportError:
        _write_txt_fallback(pass_primers, state, path.replace(".pdf", ".txt"))
        raise RuntimeError(
            "reportlab is not installed — saved as plain text instead.\n"
            "Install with: pip install reportlab"
        )


def _write_pdf_reportlab(pass_primers, state, path):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable,
    )
    import datetime

    doc    = SimpleDocTemplate(path, pagesize=A4,
                               leftMargin=2*cm, rightMargin=2*cm,
                               topMargin=2*cm, bottomMargin=2*cm)
    styles = getSampleStyleSheet()
    acc    = colors.HexColor("#4C9BE8")

    t_style  = ParagraphStyle("T",  parent=styles["Title"],   textColor=acc,  fontSize=20, spaceAfter=4)
    h1_style = ParagraphStyle("H1", parent=styles["Heading1"], textColor=acc, fontSize=13, spaceBefore=12, spaceAfter=4)
    h2_style = ParagraphStyle("H2", parent=styles["Heading2"], textColor=colors.HexColor("#888888"), fontSize=10, spaceBefore=8, spaceAfter=2)
    body     = ParagraphStyle("B",  parent=styles["Normal"],  fontSize=9, leading=13)
    small    = ParagraphStyle("S",  parent=styles["Normal"],  fontSize=7, textColor=colors.HexColor("#888888"))

    story = []

    # Header
    story.append(Paragraph("GRACE — Final SSR Marker Report", t_style))
    story.append(Paragraph("Genomic Repeat Analysis and Characterisation Engine", small))
    story.append(HRFlowable(width="100%", color=acc, thickness=1))
    story.append(Spacer(1, 0.3*cm))

    # Run info
    story.append(Paragraph("Run Information", h1_style))
    info = [
        ["Genome file",    state.genome_filename or "—"],
        ["Annotation",     state.gff_filename or "None loaded"],
        ["Design mode",    "GBS (Genotype-by-Sequencing)" if getattr(state, "gbs_mode", False) else "Standard"],
        ["Report date",    datetime.datetime.now().strftime("%Y-%m-%d %H:%M")],
    ]
    story.append(Table(info, colWidths=[4*cm, 13*cm], style=TableStyle([
        ("FONTSIZE", (0,0),(-1,-1), 9),
        ("TEXTCOLOR", (0,0),(0,-1), colors.HexColor("#888888")),
        ("BOTTOMPADDING", (0,0),(-1,-1), 4),
    ])))
    story.append(Spacer(1, 0.3*cm))

    # Pipeline summary
    story.append(Paragraph("Pipeline Summary", h1_style))
    n_ssrs    = len(state.ssrs) if state.ssrs else 0
    n_primers = len(state.primer_results) if state.primer_results else 0
    n_filtered = len(state.filtered_primer_results) if state.filtered_primer_results else n_primers
    n_spec    = len(state.specificity_results) if state.specificity_results else 0
    n_pass    = len(pass_primers)
    pass_rate = n_pass / n_spec * 100 if n_spec else 0

    pipeline = [
        ["Stage", "Count", "Notes"],
        ["SSRs detected",            f"{n_ssrs:,}",     "Raw SSRs from genome scan"],
        ["Primer pairs designed",    f"{n_primers:,}",  "After low-complexity filter"],
        ["After quality filters",    f"{n_filtered:,}", "3' stability and GC clamp filters"],
        ["Submitted to BLAST",       f"{n_spec:,}",     "Specificity checked against genome"],
        ["PASS (specific markers)",  f"{n_pass:,}",     f"{pass_rate:.1f}% pass rate"],
    ]
    story.append(Table(pipeline, colWidths=[5*cm, 3*cm, 9*cm], style=TableStyle([
        ("FONTSIZE",      (0,0),(-1,-1), 9),
        ("FONTNAME",      (0,0),(-1,0), "Helvetica-Bold"),
        ("BACKGROUND",    (0,0),(-1,0), acc),
        ("TEXTCOLOR",     (0,0),(-1,0), colors.white),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [colors.white, colors.HexColor("#F4F8FF")]),
        ("GRID",          (0,0),(-1,-1), 0.25, colors.HexColor("#DDDDDD")),
        ("BOTTOMPADDING", (0,0),(-1,-1), 4),
        ("FONTNAME",      (0,5),(0,5), "Helvetica-Bold"),
        ("TEXTCOLOR",     (1,5),(1,5), colors.HexColor("#56C596")),
    ])))
    story.append(Spacer(1, 0.4*cm))

    # SSR breakdown
    if state.ssrs:
        story.append(Paragraph("SSR Detection Results", h1_style))
        from collections import Counter
        motif_labels = {2:"Dinucleotide",3:"Trinucleotide",4:"Tetranucleotide",
                        5:"Pentanucleotide",6:"Hexanucleotide"}
        motif_counts = Counter(len(s["motif"]) for s in state.ssrs)
        motif_data   = [[motif_labels.get(k,f"{k}bp"), f"{v:,}", f"{v/n_ssrs*100:.1f}%"]
                        for k, v in sorted(motif_counts.items())]
        motif_table  = [["Motif type", "Count", "% of SSRs"]] + motif_data
        story.append(Table(motif_table, colWidths=[5*cm, 3*cm, 4*cm], style=TableStyle([
            ("FONTSIZE",      (0,0),(-1,-1), 9),
            ("FONTNAME",      (0,0),(-1,0), "Helvetica-Bold"),
            ("BACKGROUND",    (0,0),(-1,0), colors.HexColor("#2A2A4A")),
            ("TEXTCOLOR",     (0,0),(-1,0), colors.white),
            ("ROWBACKGROUNDS",(0,1),(-1,-1), [colors.white, colors.HexColor("#F4F8FF")]),
            ("GRID",          (0,0),(-1,-1), 0.25, colors.HexColor("#DDDDDD")),
            ("BOTTOMPADDING", (0,0),(-1,-1), 4),
        ])))

        if state.ssrs and "genomic_feature" in state.ssrs[0]:
            story.append(Spacer(1, 0.2*cm))
            story.append(Paragraph("Genomic Feature Distribution (SSRs)", h2_style))
            feat_counts = Counter(s.get("genomic_feature","unknown") for s in state.ssrs)
            feat_data   = [[f.capitalize(), f"{v:,}", f"{v/n_ssrs*100:.1f}%"]
                           for f, v in sorted(feat_counts.items(), key=lambda x:-x[1])]
            feat_table  = [["Feature", "SSR count", "% of total"]] + feat_data
            story.append(Table(feat_table, colWidths=[4*cm, 3*cm, 4*cm], style=TableStyle([
                ("FONTSIZE",      (0,0),(-1,-1), 9),
                ("FONTNAME",      (0,0),(-1,0), "Helvetica-Bold"),
                ("BACKGROUND",    (0,0),(-1,0), colors.HexColor("#2A2A4A")),
                ("TEXTCOLOR",     (0,0),(-1,0), colors.white),
                ("ROWBACKGROUNDS",(0,1),(-1,-1), [colors.white, colors.HexColor("#F4F8FF")]),
                ("GRID",          (0,0),(-1,-1), 0.25, colors.HexColor("#DDDDDD")),
                ("BOTTOMPADDING", (0,0),(-1,-1), 4),
            ])))
        story.append(Spacer(1, 0.3*cm))

    # PASS primers table
    story.append(Paragraph("PASS Primer Pairs", h1_style))
    story.append(Paragraph(
        f"{n_pass:,} primer pairs passed specificity checking. "
        "Pair rank 0 is Primer3's top-ranked pair per SSR.",
        body
    ))
    story.append(Spacer(1, 0.2*cm))

    has_feat = pass_primers and "genomic_feature" in pass_primers[0]
    headers  = ["SSR ID","Contig","Motif","Pair","Forward (5'→3')","Reverse (5'→3')","Size","Fwd Tm","Rev Tm"]
    col_w    = [1.2*cm, 3*cm, 1.4*cm, 0.8*cm, 4.8*cm, 4.8*cm, 1.2*cm, 1.2*cm, 1.2*cm]
    if has_feat:
        headers.append("Feature"); col_w.append(1.8*cm)

    rows = [headers]
    for p in pass_primers[:4000]:
        row = [
            str(p.get("ssr_id","")),
            str(p.get("contig",""))[-22:],
            str(p.get("canonical_motif", p.get("motif",""))),
            str(p.get("pair_rank",0)),
            str(p.get("left_primer","")),
            str(p.get("right_primer","")),
            str(p.get("product_size","")),
            f"{p.get('left_tm',0):.1f}",
            f"{p.get('right_tm',0):.1f}",
        ]
        if has_feat: row.append(str(p.get("genomic_feature","")))
        rows.append(row)

    if len(pass_primers) > 4000:
        rows.append(["...",f"({len(pass_primers)-4000:,} more — export CSV for full data)","","","","","","",""])

    story.append(Table(rows, colWidths=col_w, repeatRows=1, style=TableStyle([
        ("FONTSIZE",       (0,0),(-1,-1), 7),
        ("FONTNAME",       (0,0),(-1,0), "Helvetica-Bold"),
        ("BACKGROUND",     (0,0),(-1,0), acc),
        ("TEXTCOLOR",      (0,0),(-1,0), colors.white),
        ("ROWBACKGROUNDS", (0,1),(-1,-1), [colors.white, colors.HexColor("#F4F8FF")]),
        ("GRID",           (0,0),(-1,-1), 0.25, colors.HexColor("#DDDDDD")),
        ("BOTTOMPADDING",  (0,0),(-1,-1), 3),
        ("TOPPADDING",     (0,0),(-1,-1), 3),
        ("FONTNAME",       (4,1),(5,-1), "Courier"),
    ])))

    story.append(Spacer(1, 0.5*cm))
    story.append(HRFlowable(width="100%", color=colors.HexColor("#DDDDDD"), thickness=0.5))
    story.append(Paragraph(
        "Generated by GRACE — Genomic Repeat Analysis and Characterisation Engine. "
        "Primers designed using Primer3. Specificity verified by BLAST+ against the reference genome.",
        small
    ))
    doc.build(story)


def _write_txt_fallback(pass_primers, state, path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("GRACE — Final SSR Marker Report\n" + "="*60 + "\n\n")
        f.write(f"Genome: {state.genome_filename or '—'}\n")
        f.write(f"PASS primers: {len(pass_primers):,}\n\n")
        f.write("\t".join(["SSR_ID","Contig","Motif","Pair","Forward","Reverse",
                            "Product","Fwd_Tm","Rev_Tm","Feature"]) + "\n")
        for p in pass_primers:
            f.write("\t".join([
                str(p.get("ssr_id","")), str(p.get("contig","")),
                str(p.get("canonical_motif", p.get("motif",""))),
                str(p.get("pair_rank",0)), str(p.get("left_primer","")),
                str(p.get("right_primer","")), str(p.get("product_size","")),
                f"{p.get('left_tm',0):.1f}", f"{p.get('right_tm',0):.1f}",
                str(p.get("genomic_feature","")),
            ]) + "\n")


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class ReportPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._worker       = None
        self._last_version = -1
        self._pass_primers = []
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
        title = QLabel("Final Report")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        sub = QLabel("Complete pipeline summary — from genome scan to validated SSR markers")
        sub.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title); L.addWidget(sub)

        self._placeholder = QLabel("Complete Specificity Check first to generate the final report.")
        self._placeholder.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(self._placeholder)

        # ── Key metrics cards ─────────────────────────────
        self._cards_widget = QWidget()
        self._cards_widget.setVisible(False)
        cards_row = QHBoxLayout(self._cards_widget)
        cards_row.setSpacing(10)
        self._card_ssrs     = _stat_card("SSRs Detected", "—", ACCENT)
        self._card_primers  = _stat_card("Primers Designed", "—", ACCENT)
        self._card_pass     = _stat_card("PASS Markers", "—", SUCCESS)
        self._card_rate     = _stat_card("Pass Rate", "—", SUCCESS)
        self._card_motif    = _stat_card("Top Motif", "—", ACCENT)
        for c in [self._card_ssrs, self._card_primers, self._card_pass,
                  self._card_rate, self._card_motif]:
            cards_row.addWidget(c)
        L.addWidget(self._cards_widget)

        # ── Pipeline funnel ───────────────────────────────
        self._funnel_group = QGroupBox("Filtering Pipeline")
        self._funnel_group.setVisible(False)
        fg = QVBoxLayout(self._funnel_group)
        self._funnel = FunnelWidget()
        fg.addWidget(self._funnel)
        L.addWidget(self._funnel_group)

        # ── SSR breakdown charts ──────────────────────────
        self._ssr_charts_group = QGroupBox("SSR Detection Breakdown")
        self._ssr_charts_group.setVisible(False)
        scg = QHBoxLayout(self._ssr_charts_group)
        self._motif_pie    = MiniPieChart("Motif Type Distribution")
        self._feature_pie  = MiniPieChart("Genomic Feature Distribution")
        scg.addWidget(self._motif_pie)
        scg.addWidget(self._feature_pie)
        L.addWidget(self._ssr_charts_group)

        # ── PASS primer breakdown charts ──────────────────
        self._pass_charts_group = QGroupBox("PASS Primer Breakdown")
        self._pass_charts_group.setVisible(False)
        pcg = QHBoxLayout(self._pass_charts_group)
        self._pass_motif_pie   = MiniPieChart("Motif Types (PASS)")
        self._pass_feature_pie = MiniPieChart("Genomic Features (PASS)")
        self._product_bar      = MiniBarChart("Product Size Distribution")
        pcg.addWidget(self._pass_motif_pie)
        pcg.addWidget(self._pass_feature_pie)
        pcg.addWidget(self._product_bar)
        L.addWidget(self._pass_charts_group)

        # ── Chromosome coverage ───────────────────────────
        self._chrom_group = QGroupBox("Coverage per Chromosome")
        self._chrom_group.setVisible(False)
        cg = QVBoxLayout(self._chrom_group)
        self._chrom_note = QLabel("")
        self._chrom_note.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        cg.addWidget(self._chrom_note)
        self._chrom_bar = MiniBarChart("PASS markers per chromosome")
        cg.addWidget(self._chrom_bar)
        L.addWidget(self._chrom_group)

        # ── Export ────────────────────────────────────────
        self._export_group = QGroupBox("Export")
        self._export_group.setVisible(False)
        eg = QHBoxLayout(self._export_group)
        self._pdf_btn   = QPushButton("Export PDF Report")
        self._pdf_btn.setObjectName("primary")
        self._pdf_btn.clicked.connect(self._export_pdf)
        self._csv_btn   = QPushButton("Export CSV")
        self._csv_btn.clicked.connect(self._export_csv)
        self._fasta_btn = QPushButton("Export FASTA")
        self._fasta_btn.clicked.connect(self._export_fasta)
        self._gbs_export_btn = QPushButton("Export GBS FASTA (tailed)")
        self._gbs_export_btn.clicked.connect(self._export_gbs_fasta)
        self._gbs_export_btn.setVisible(False)
        eg.addWidget(self._pdf_btn)
        eg.addWidget(self._csv_btn)
        eg.addWidget(self._fasta_btn)
        eg.addWidget(self._gbs_export_btn)
        eg.addStretch()
        self._export_status = QLabel("")
        self._export_status.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        eg.addWidget(self._export_status)
        L.addWidget(self._export_group)

        # ── Table ─────────────────────────────────────────
        self._table_group = QGroupBox("PASS Primer Pairs")
        self._table_group.setVisible(False)
        tg = QVBoxLayout(self._table_group)
        self._metrics_label = QLabel("")
        self._metrics_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        tg.addWidget(self._metrics_label)
        self._table = QTableWidget()
        self._table.setAlternatingRowColors(True)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        self._table.setFixedHeight(500)
        tg.addWidget(self._table)
        L.addWidget(self._table_group)

        L.addStretch()

    # ---------------------------------------------------------
    # REBUILD
    # ---------------------------------------------------------
    def _rebuild(self):
        if not self.state.has_specificity:
            self._placeholder.setVisible(True)
            for w in [self._cards_widget, self._funnel_group, self._ssr_charts_group,
                      self._pass_charts_group, self._chrom_group,
                      self._export_group, self._table_group]:
                w.setVisible(False)
            return

        self._placeholder.setVisible(False)
        spec    = self.state.specificity_results or []
        primers = self.state.primer_results or []
        ssrs    = self.state.ssrs or []

        pass_ids = {r["ssr_id"] for r in spec if r.get("specificity_status") == "PASS"}
        self._pass_primers = [p for p in primers if p.get("ssr_id") in pass_ids]

        n_ssrs     = len(ssrs)
        n_primers  = len(primers)
        n_filtered = len(self.state.filtered_primer_results or primers)
        n_spec     = len(spec)
        n_pass     = len(pass_ids)
        pass_rate  = n_pass / n_spec * 100 if n_spec else 0

        motif_counts = Counter(
            p.get("canonical_motif", p.get("motif", "?")) for p in self._pass_primers
        )
        top_motif = motif_counts.most_common(1)[0][0] if motif_counts else "—"

        # ── Key metric cards ──────────────────────────────
        def _update_card(card, value):
            card.findChild(QLabel).setText(value)

        # Update value labels (first child = value label)
        cards = [
            (self._card_ssrs,    f"{n_ssrs:,}"),
            (self._card_primers, f"{n_primers:,}"),
            (self._card_pass,    f"{n_pass:,}"),
            (self._card_rate,    f"{pass_rate:.1f}%"),
            (self._card_motif,   top_motif),
        ]
        for card, val in cards:
            labels = card.findChildren(QLabel)
            if labels: labels[0].setText(val)
        self._cards_widget.setVisible(True)

        # ── Funnel ────────────────────────────────────────
        self._funnel.set_stages([
            ("SSRs Detected",        n_ssrs,     "#4C9BE8"),
            ("Primers Designed",     n_primers,  "#9B89C4"),
            ("After Quality Filter", n_filtered, "#F4A261"),
            ("BLAST Tested",         n_spec,     "#48CAE4"),
            ("PASS Markers",         n_pass,     "#56C596"),
        ])
        self._funnel_group.setVisible(True)

        # ── SSR motif pie ─────────────────────────────────
        if ssrs:
            mlabels = {2:"Di",3:"Tri",4:"Tetra",5:"Penta",6:"Hexa"}
            mc = Counter(len(s["motif"]) for s in ssrs)
            self._motif_pie.set_data([
                (mlabels.get(k,f"{k}bp"), v, CHART_COLOURS[i % len(CHART_COLOURS)])
                for i,(k,v) in enumerate(sorted(mc.items()))
            ])

            if "genomic_feature" in ssrs[0]:
                fc = Counter(s.get("genomic_feature","unknown") for s in ssrs)
                self._feature_pie.set_data([
                    (f.capitalize(), v, FEATURE_COLOURS.get(f,"#888888"))
                    for f,v in sorted(fc.items(), key=lambda x:-x[1])
                ])
            else:
                self._feature_pie.set_data([("No GFF loaded", 1, BORDER)])

            self._ssr_charts_group.setVisible(True)

        # ── PASS primer charts ────────────────────────────
        if self._pass_primers:
            pm = Counter(p.get("canonical_motif", p.get("motif","?")) for p in self._pass_primers)
            top_n = pm.most_common(8)
            self._pass_motif_pie.set_data([
                (m, c, CHART_COLOURS[i % len(CHART_COLOURS)])
                for i,(m,c) in enumerate(top_n)
            ])

            if "genomic_feature" in self._pass_primers[0]:
                pfc = Counter(p.get("genomic_feature","unknown") for p in self._pass_primers)
                self._pass_feature_pie.set_data([
                    (f.capitalize(), v, FEATURE_COLOURS.get(f,"#888888"))
                    for f,v in sorted(pfc.items(), key=lambda x:-x[1])
                ])
            else:
                self._pass_feature_pie.set_data([("No annotation", 1, BORDER)])

            # Product size bins
            size_bins = defaultdict(int)
            for p in self._pass_primers:
                sz = p.get("product_size", 0)
                if sz < 100:   size_bins["<100bp"]   += 1
                elif sz < 150: size_bins["100–150bp"] += 1
                elif sz < 200: size_bins["150–200bp"] += 1
                elif sz < 250: size_bins["200–250bp"] += 1
                elif sz < 300: size_bins["250–300bp"] += 1
                else:          size_bins[">300bp"]    += 1
            bin_order = ["<100bp","100–150bp","150–200bp","200–250bp","250–300bp",">300bp"]
            self._product_bar.set_data([
                (b, size_bins[b], CHART_COLOURS[i % len(CHART_COLOURS)])
                for i,b in enumerate(bin_order) if size_bins[b]
            ])
            self._pass_charts_group.setVisible(True)

        # ── Chromosome coverage ───────────────────────────
        chrom_counts = Counter(p.get("contig","") for p in self._pass_primers)
        if len(chrom_counts) <= 50:
            top_chroms = chrom_counts.most_common(20)
            self._chrom_bar.set_data([
                (c[-12:], n, CHART_COLOURS[i % len(CHART_COLOURS)])
                for i,(c,n) in enumerate(top_chroms)
            ])
            self._chrom_note.setText(
                f"PASS markers distributed across {len(chrom_counts):,} "
                f"sequence{'s' if len(chrom_counts)!=1 else ''}. "
                f"Showing top {min(20,len(chrom_counts))}."
            )
            self._chrom_group.setVisible(True)

        # ── GBS export button ─────────────────────────────
        has_tails = self._pass_primers and "left_primer_tailed" in self._pass_primers[0]
        self._gbs_export_btn.setVisible(has_tails)

        self._export_group.setVisible(True)
        self._populate_table()
        self._table_group.setVisible(True)

    def _populate_table(self):
        primers   = self._pass_primers
        truncated = len(primers) > TABLE_DISPLAY_LIMIT
        display   = primers[:TABLE_DISPLAY_LIMIT] if truncated else primers

        self._metrics_label.setText(
            f"{len(primers):,} PASS primer pairs" +
            (f"   |   showing first {TABLE_DISPLAY_LIMIT:,} — export for full data" if truncated else "")
        )

        has_feature = primers and "genomic_feature" in primers[0]
        COLS = {
            "ssr_id": "SSR ID", "pair_rank": "Pair", "contig": "Contig",
            "motif": "Motif", "left_primer": "Forward primer",
            "right_primer": "Reverse primer", "product_size": "Product (bp)",
            "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC%", "right_gc": "Rev GC%",
        }
        if has_feature:
            COLS["genomic_feature"] = "Feature"

        display_cols = list(COLS.keys())
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(display))
        self._table.setColumnCount(len(display_cols))
        self._table.setHorizontalHeaderLabels(list(COLS.values()))

        header = self._table.horizontalHeader()
        header.setStretchLastSection(True)
        for ci, col in enumerate(display_cols):
            if col == "contig":
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self._table.setColumnWidth(ci, 160)
            elif col in ("left_primer","right_primer"):
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.Interactive)
                self._table.setColumnWidth(ci, 200)
            else:
                header.setSectionResizeMode(ci, QHeaderView.ResizeMode.ResizeToContents)

        for row_idx, p in enumerate(display):
            for col_idx, col in enumerate(display_cols):
                val  = p.get(col, "")
                text = f"{val:.2f}" if isinstance(val, float) else ("" if val is None else str(val))
                self._table.setItem(row_idx, col_idx, QTableWidgetItem(text))

        self._table.setSortingEnabled(True)

    # ---------------------------------------------------------
    # EXPORTS
    # ---------------------------------------------------------
    def _export_pdf(self):
        if not self._pass_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Export PDF", "GRACE_report.pdf", "PDF files (*.pdf)")
        if not path: return
        self._pdf_btn.setEnabled(False)
        self._set_export_status("Generating PDF...", TEXT_SECONDARY)
        self._worker = PDFWorker(self._pass_primers, self.state, path)
        self._worker.finished.connect(self._on_pdf_done)
        self._worker.error.connect(self._on_pdf_error)
        self._worker.start()

    def _on_pdf_done(self, path):
        self._pdf_btn.setEnabled(True)
        self._set_export_status(f"PDF saved to {path}", SUCCESS)
        self.mw.set_status(f"Report saved to {path}")
        self._worker = None

    def _on_pdf_error(self, msg):
        self._pdf_btn.setEnabled(True)
        self._set_export_status(f"PDF error: {msg}", ERROR)
        self._worker = None

    def _export_csv(self):
        if not self._pass_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", "GRACE_pass_primers.csv", "CSV files (*.csv)")
        if not path: return
        try:
            import pandas as pd
            pd.DataFrame(self._pass_primers).rename(columns={
                "ssr_id":"SSR ID","pair_rank":"Pair rank","contig":"Contig",
                "motif":"Motif","canonical_motif":"Canonical motif",
                "repeat_count":"Repeat count","left_primer":"Forward primer",
                "right_primer":"Reverse primer","product_size":"Product size (bp)",
                "left_tm":"Forward Tm (°C)","right_tm":"Reverse Tm (°C)",
                "left_gc":"Forward GC (%)","right_gc":"Reverse GC (%)",
                "left_3end_dg":"Forward 3' stability (kcal/mol)",
                "right_3end_dg":"Reverse 3' stability (kcal/mol)",
                "genomic_feature":"Genomic feature",
            }).to_csv(path, index=False, encoding="utf-8-sig")
            self._set_export_status(f"CSV saved to {path}", SUCCESS)
            self.mw.set_status(f"CSV saved to {path}")
        except Exception as e:
            self._set_export_status(f"CSV error: {e}", ERROR)

    def _export_fasta(self):
        if not self._pass_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Export FASTA", "GRACE_pass_primers.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        try:
            from core.primer_design import primers_to_blast_fasta
            with open(path, "w") as f:
                f.write(primers_to_blast_fasta(self._pass_primers))
            self._set_export_status(f"FASTA saved to {path}", SUCCESS)
            self.mw.set_status(f"FASTA saved to {path}")
        except Exception as e:
            self._set_export_status(f"FASTA error: {e}", ERROR)

    def _export_gbs_fasta(self):
        if not self._pass_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Export GBS FASTA", "GRACE_pass_gbs_tailed.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        try:
            from core.primer_design import primers_to_gbs_fasta
            with open(path, "w") as f:
                f.write(primers_to_gbs_fasta(self._pass_primers))
            self._set_export_status(f"GBS FASTA saved to {path}", SUCCESS)
            self.mw.set_status(f"GBS FASTA saved to {path}")
        except Exception as e:
            self._set_export_status(f"GBS FASTA error: {e}", ERROR)

    def _set_export_status(self, msg, color):
        self._export_status.setText(msg)
        self._export_status.setStyleSheet(f"color: {color}; font-size: {FONT_SIZE_SMALL}pt;")

    def on_show(self):
        if self.state.blast_version != self._last_version:
            self._last_version = self.state.blast_version
            self._rebuild()