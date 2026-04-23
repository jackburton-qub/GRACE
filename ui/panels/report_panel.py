"""
report_panel.py — Project Summary & PDF Export
Enhanced UI with detailed tables and comprehensive PDF report.
Includes Panel column for Capillary results.
"""

import os
from collections import Counter
from datetime import datetime
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QFileDialog, QMessageBox,
    QTabWidget,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BORDER,
)

from ui.panels.ssr_summary_panel import PieChart, VerticalBarChart

CHART_COLORS = [
    "#4C9BE8", "#56C596", "#F4A261", "#E76F51",
    "#9B89C4", "#48CAE4", "#F6C90E", "#E07BA5",
]


class ReportPanel(QWidget):
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
        self._main_layout = QVBoxLayout(content)
        self._main_layout.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        self._main_layout.setSpacing(16)

        title = QLabel("Project Summary")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {ACCENT};")
        self._main_layout.addWidget(title)

        desc = QLabel("Review your results and generate a comprehensive PDF report.")
        desc.setStyleSheet(f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;")
        desc.setWordWrap(True)
        self._main_layout.addWidget(desc)

        btn_layout = QHBoxLayout()
        self._pdf_btn = QPushButton("📄 Generate PDF Report")
        self._pdf_btn.setObjectName("primary")
        self._pdf_btn.setMinimumHeight(40)
        self._pdf_btn.clicked.connect(self._generate_pdf)
        btn_layout.addWidget(self._pdf_btn)
        btn_layout.addStretch()
        self._main_layout.addLayout(btn_layout)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        self._main_layout.addWidget(self._status_label)

        self._tab_widget = QTabWidget()
        self._main_layout.addWidget(self._tab_widget)

        self._main_layout.addStretch()

    def _refresh(self):
        while self._tab_widget.count() > 0:
            self._tab_widget.removeTab(0)

        self._add_overview_tab()

        all_primers = self.state.primer_results
        if all_primers:
            self._add_primer_tab(all_primers, "All Primers")

        pass_primers = self._get_blast_pass_primers()
        if pass_primers:
            self._add_primer_tab(pass_primers, "PASS Primers")

        if self.state.amplicon_validation_result:
            self._add_amplicon_tab()

        if self.state.capillary_result:
            self._add_capillary_tab()

        if self.state.gbs_re_markers:
            self._add_gbs_re_tab()

        if not any([self.state.has_ssrs, all_primers]):
            self._status_label.setText("No data available. Complete previous steps first.")

    def _get_blast_pass_primers(self):
        spec_results = self.state.specificity_results
        all_primers = self.state.primer_results
        if not spec_results or not all_primers:
            return []
        pass_keys = set()
        for r in spec_results:
            if r.get("specificity_status") == "PASS":
                ssr_id = r.get("ssr_id")
                pair_rank = r.get("pair_rank")
                if ssr_id is not None and pair_rank is not None:
                    pass_keys.add((ssr_id, pair_rank))
        return [p for p in all_primers if (p.get("ssr_id"), p.get("pair_rank")) in pass_keys]

    def _add_overview_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setSpacing(16)

        if self.state.has_ssrs:
            ssrs = self.state.ssrs
            total = len(ssrs)
            motifs = [s.get("canonical_motif", s.get("motif", "")) for s in ssrs]
            motif_counts = Counter(motifs)

            group = QGroupBox("SSR Detection Summary")
            gl = QVBoxLayout(group)

            metrics_layout = QHBoxLayout()
            metrics = [
                ("Total SSRs", f"{total:,}"),
                ("Unique motifs", f"{len(motif_counts):,}"),
                ("Contigs with SSRs", f"{len(set(s.get('contig') for s in ssrs)):,}"),
            ]
            for label, value in metrics:
                card = self._make_metric_card(label, value)
                metrics_layout.addWidget(card)
            gl.addLayout(metrics_layout)

            motif_types = self._get_motif_type_counts(ssrs)
            pie = PieChart()
            pie.set_data(motif_types)
            pie.setMinimumHeight(180)
            gl.addWidget(QLabel("Motif Type Distribution"))
            gl.addWidget(pie)

            top_motifs = motif_counts.most_common(10)
            bar = VerticalBarChart()
            bar.set_data(top_motifs)
            bar.setMinimumHeight(200)
            gl.addWidget(QLabel("Top 10 Canonical Motifs"))
            gl.addWidget(bar)

            layout.addWidget(group)

        all_primers = self.state.primer_results
        pass_primers = self._get_blast_pass_primers()
        if all_primers:
            group = QGroupBox("Primer Design Summary")
            gl = QVBoxLayout(group)
            metrics_layout = QHBoxLayout()
            metrics = [
                ("Total pairs", f"{len(all_primers):,}"),
                ("PASS pairs", f"{len(pass_primers):,}" if pass_primers else "—"),
                ("Unique loci", f"{len(set(p.get('ssr_id') for p in all_primers)):,}"),
            ]
            for label, value in metrics:
                card = self._make_metric_card(label, value)
                metrics_layout.addWidget(card)
            gl.addLayout(metrics_layout)
            layout.addWidget(group)

        layout.addStretch()
        self._tab_widget.addTab(tab, "Overview")

    def _get_motif_type_counts(self, ssrs):
        counts = Counter()
        for ssr in ssrs:
            motif = ssr.get("canonical_motif", ssr.get("motif", ""))
            length = len(motif)
            if length == 2:
                counts["Di"] += 1
            elif length == 3:
                counts["Tri"] += 1
            elif length == 4:
                counts["Tetra"] += 1
            elif length == 5:
                counts["Penta"] += 1
            elif length >= 6:
                counts["Hexa"] += 1
        return {k: v for k, v in counts.items() if v > 0}

    def _add_primer_tab(self, primers, title):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        group = QGroupBox(f"{title} ({len(primers):,} pairs)")
        gl = QVBoxLayout(group)

        table = QTableWidget()
        table.setAlternatingRowColors(True)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSortingEnabled(True)
        cols = ["SSR ID", "Pair Rank", "Contig", "Motif", "Product Size", "Fwd Tm", "Rev Tm", "Fwd Primer", "Rev Primer"]
        show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
        if show_chromosome:
            cols.insert(2, "Chromosome")

        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)

        display_primers = primers[:100]
        table.setRowCount(len(display_primers))
        for i, p in enumerate(display_primers):
            table.setItem(i, 0, QTableWidgetItem(str(p.get("ssr_id", ""))))
            table.setItem(i, 1, QTableWidgetItem(str(p.get("pair_rank", ""))))
            contig = p.get("contig", "")
            if show_chromosome:
                table.setItem(i, 2, QTableWidgetItem(self.state.get_display_name(contig)))
                table.setItem(i, 3, QTableWidgetItem(contig))
                table.setItem(i, 4, QTableWidgetItem(p.get("motif", "")))
                table.setItem(i, 5, QTableWidgetItem(str(p.get("product_size", ""))))
                table.setItem(i, 6, QTableWidgetItem(f"{p.get('left_tm', 0):.1f}"))
                table.setItem(i, 7, QTableWidgetItem(f"{p.get('right_tm', 0):.1f}"))
                table.setItem(i, 8, QTableWidgetItem(p.get("left_primer", "")))
                table.setItem(i, 9, QTableWidgetItem(p.get("right_primer", "")))
            else:
                table.setItem(i, 2, QTableWidgetItem(contig))
                table.setItem(i, 3, QTableWidgetItem(p.get("motif", "")))
                table.setItem(i, 4, QTableWidgetItem(str(p.get("product_size", ""))))
                table.setItem(i, 5, QTableWidgetItem(f"{p.get('left_tm', 0):.1f}"))
                table.setItem(i, 6, QTableWidgetItem(f"{p.get('right_tm', 0):.1f}"))
                table.setItem(i, 7, QTableWidgetItem(p.get("left_primer", "")))
                table.setItem(i, 8, QTableWidgetItem(p.get("right_primer", "")))
        table.resizeColumnsToContents()
        gl.addWidget(table)

        if len(primers) > 100:
            gl.addWidget(QLabel(f"Showing first 100 of {len(primers):,} primer pairs."))

        layout.addWidget(group)
        self._tab_widget.addTab(tab, title)

    def _add_amplicon_tab(self):
        result = self.state.amplicon_validation_result
        if not result:
            return
        tab = QWidget()
        layout = QVBoxLayout(tab)

        group1 = QGroupBox("Validation Summary")
        gl1 = QVBoxLayout(group1)
        clean_pool = result.get("clean_pool", [])
        issues = result.get("issues_summary", {})
        metrics_layout = QHBoxLayout()
        metrics = [
            ("Input loci", str(result.get("n_input", 0))),
            ("Passed checks", str(result.get("n_clean", 0))),
            ("Final pool", str(result.get("n_final", len(clean_pool)))),
            ("Dimer risks", str(issues.get("dimer_risks", 0))),
            ("Size filtered", str(issues.get("size_filtered_short", 0) + issues.get("size_filtered_long", 0))),
            ("Compatibility score", f"{result.get('compatibility_score', 0):.2f}"),
        ]
        for label, value in metrics:
            card = self._make_metric_card(label, value)
            metrics_layout.addWidget(card)
        gl1.addLayout(metrics_layout)
        layout.addWidget(group1)

        if clean_pool:
            group2 = QGroupBox(f"Clean Pool ({len(clean_pool)} loci)")
            gl2 = QVBoxLayout(group2)
            table = QTableWidget()
            cols = ["SSR ID", "Contig", "Motif", "Product Size", "Fwd Tm", "Rev Tm", "Fwd Primer", "Rev Primer"]
            show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
            if show_chromosome:
                cols.insert(1, "Chromosome")
            table.setColumnCount(len(cols))
            table.setHorizontalHeaderLabels(cols)
            table.setRowCount(min(len(clean_pool), 50))
            for i, p in enumerate(clean_pool[:50]):
                table.setItem(i, 0, QTableWidgetItem(str(p.get("ssr_id", ""))))
                contig = p.get("contig", "")
                if show_chromosome:
                    table.setItem(i, 1, QTableWidgetItem(self.state.get_display_name(contig)))
                    table.setItem(i, 2, QTableWidgetItem(contig))
                    table.setItem(i, 3, QTableWidgetItem(p.get("motif", "")))
                    table.setItem(i, 4, QTableWidgetItem(str(p.get("product_size", ""))))
                    table.setItem(i, 5, QTableWidgetItem(f"{p.get('left_tm', 0):.1f}"))
                    table.setItem(i, 6, QTableWidgetItem(f"{p.get('right_tm', 0):.1f}"))
                    table.setItem(i, 7, QTableWidgetItem(p.get("left_primer", "")))
                    table.setItem(i, 8, QTableWidgetItem(p.get("right_primer", "")))
                else:
                    table.setItem(i, 1, QTableWidgetItem(contig))
                    table.setItem(i, 2, QTableWidgetItem(p.get("motif", "")))
                    table.setItem(i, 3, QTableWidgetItem(str(p.get("product_size", ""))))
                    table.setItem(i, 4, QTableWidgetItem(f"{p.get('left_tm', 0):.1f}"))
                    table.setItem(i, 5, QTableWidgetItem(f"{p.get('right_tm', 0):.1f}"))
                    table.setItem(i, 6, QTableWidgetItem(p.get("left_primer", "")))
                    table.setItem(i, 7, QTableWidgetItem(p.get("right_primer", "")))
            table.resizeColumnsToContents()
            gl2.addWidget(table)
            layout.addWidget(group2)

        layout.addStretch()
        self._tab_widget.addTab(tab, "Amplicon")

    def _add_capillary_tab(self):
        result = self.state.capillary_result
        if not result:
            return
        tab = QWidget()
        layout = QVBoxLayout(tab)

        assignments = result.get("assignments", [])
        unassigned = result.get("unassigned", [])
        n_total = result.get("n_total", 0)
        n_assigned = result.get("n_assigned", len(assignments))

        group1 = QGroupBox("Dye Assignment Summary")
        gl1 = QVBoxLayout(group1)
        metrics_layout = QHBoxLayout()
        metrics = [
            ("Total loci", str(n_total)),
            ("Assigned", str(n_assigned)),
            ("Unassigned", str(len(unassigned))),
            ("Dyes used", str(len(set(a.get("dye", "") for a in assignments)))),
        ]
        for label, value in metrics:
            card = self._make_metric_card(label, value)
            metrics_layout.addWidget(card)
        gl1.addLayout(metrics_layout)
        layout.addWidget(group1)

        if assignments:
            group2 = QGroupBox(f"Assignments ({len(assignments)} loci)")
            gl2 = QVBoxLayout(group2)
            table = QTableWidget()
            cols = ["SSR ID", "Dye", "Min Size", "Max Size", "Range", "Contig", "Fwd Primer", "Rev Primer"]
            has_panels = any("panel" in a for a in assignments)
            if has_panels:
                cols.insert(0, "Panel")
            show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
            if show_chromosome:
                cols.insert(cols.index("Contig") + 1, "Chromosome")
            table.setColumnCount(len(cols))
            table.setHorizontalHeaderLabels(cols)
            table.setRowCount(min(len(assignments), 50))
            for i, a in enumerate(assignments[:50]):
                col_idx = 0
                if has_panels:
                    table.setItem(i, col_idx, QTableWidgetItem(str(a.get("panel", ""))))
                    col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a.get("ssr_id", ""))))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(a.get("dye", "")))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a.get("min_size", ""))))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a.get("max_size", ""))))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(str(a.get("max_size", 0) - a.get("min_size", 0))))
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
                table.setItem(i, col_idx, QTableWidgetItem(fwd))
                col_idx += 1
                table.setItem(i, col_idx, QTableWidgetItem(rev))
            table.resizeColumnsToContents()
            gl2.addWidget(table)
            layout.addWidget(group2)

        layout.addStretch()
        self._tab_widget.addTab(tab, "Capillary")

    def _add_gbs_re_tab(self):
        markers = self.state.gbs_re_markers
        passing = self.state.gbs_re_passing_frags
        tab = QWidget()
        layout = QVBoxLayout(tab)

        group1 = QGroupBox("GBS‑RE Marker Discovery Summary")
        gl1 = QVBoxLayout(group1)
        metrics_layout = QHBoxLayout()
        metrics = [
            ("Qualified SSRs", f"{len(markers):,}"),
            ("Passing fragments", f"{passing:,}"),
            ("Unique contigs", f"{len(set(m.get('contig') for m in markers)):,}"),
        ]
        for label, value in metrics:
            card = self._make_metric_card(label, value)
            metrics_layout.addWidget(card)
        gl1.addLayout(metrics_layout)
        layout.addWidget(group1)

        if markers:
            group2 = QGroupBox(f"Markers ({len(markers)} loci)")
            gl2 = QVBoxLayout(group2)
            table = QTableWidget()
            cols = ["SSR ID", "Contig", "Fragment Start", "Fragment End", "Fragment Size", "Motif"]
            show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
            if show_chromosome:
                cols.insert(1, "Chromosome")
            table.setColumnCount(len(cols))
            table.setHorizontalHeaderLabels(cols)
            table.setRowCount(min(len(markers), 50))
            for i, m in enumerate(markers[:50]):
                table.setItem(i, 0, QTableWidgetItem(str(m.get("ssr_id", ""))))
                contig = m.get("contig", "")
                if show_chromosome:
                    table.setItem(i, 1, QTableWidgetItem(self.state.get_display_name(contig)))
                    table.setItem(i, 2, QTableWidgetItem(contig))
                    table.setItem(i, 3, QTableWidgetItem(str(m.get("fragment_start", ""))))
                    table.setItem(i, 4, QTableWidgetItem(str(m.get("fragment_end", ""))))
                    table.setItem(i, 5, QTableWidgetItem(str(m.get("fragment_size", ""))))
                    table.setItem(i, 6, QTableWidgetItem(m.get("motif", "")))
                else:
                    table.setItem(i, 1, QTableWidgetItem(contig))
                    table.setItem(i, 2, QTableWidgetItem(str(m.get("fragment_start", ""))))
                    table.setItem(i, 3, QTableWidgetItem(str(m.get("fragment_end", ""))))
                    table.setItem(i, 4, QTableWidgetItem(str(m.get("fragment_size", ""))))
                    table.setItem(i, 5, QTableWidgetItem(m.get("motif", "")))
            table.resizeColumnsToContents()
            gl2.addWidget(table)
            layout.addWidget(group2)

        layout.addStretch()
        self._tab_widget.addTab(tab, "GBS‑RE")

    def _make_metric_card(self, label, value):
        card = QWidget()
        card.setStyleSheet(f"background: {BG_MID}; border: 1px solid {BORDER}; border-radius: 6px; padding: 10px;")
        cl = QVBoxLayout(card)
        vl = QLabel(value)
        vl.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
        vl.setStyleSheet(f"color: {ACCENT};")
        cl.addWidget(vl)
        cl.addWidget(QLabel(label))
        return card

    def _generate_pdf(self):
        try:
            import reportlab
        except ImportError:
            QMessageBox.warning(self, "Missing Library", "ReportLab is not installed. Run: pip install reportlab")
            return

        if not self.state.has_ssrs:
            QMessageBox.warning(self, "No Data", "No SSR data available.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save PDF Report", "GRACE_report.pdf", "PDF files (*.pdf)")
        if not path:
            return

        if os.path.exists(path):
            try:
                with open(path, 'a') as f:
                    pass
            except PermissionError:
                QMessageBox.warning(self, "File In Use", f"'{os.path.basename(path)}' is open in another program.\nPlease close it and try again.")
                return

        try:
            from reportlab.lib.pagesizes import A4, landscape
            from reportlab.lib import colors
            from reportlab.lib.units import inch, mm
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
            from reportlab.platypus import (
                SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                PageBreak, PageTemplate, Frame, NextPageTemplate
            )
            from reportlab.pdfgen import canvas

            # ------------------------------------------------------------------
            # Custom page template with header and footer
            # ------------------------------------------------------------------
            def header_footer(canvas_obj, doc):
                canvas_obj.saveState()
                w, h = landscape(A4)

                # Header line
                canvas_obj.setStrokeColor(colors.HexColor(ACCENT))
                canvas_obj.setLineWidth(1)
                canvas_obj.line(doc.leftMargin, h - 0.6*inch, w - doc.rightMargin, h - 0.6*inch)

                # Header text
                canvas_obj.setFont("Helvetica-Bold", 14)
                canvas_obj.setFillColor(colors.HexColor(ACCENT))
                canvas_obj.drawString(doc.leftMargin, h - 0.4*inch, "GRACE Project Report")

                # Footer
                canvas_obj.setFont("Helvetica", 8)
                canvas_obj.setFillColor(colors.grey)
                canvas_obj.drawString(doc.leftMargin, 0.4*inch, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
                canvas_obj.drawRightString(w - doc.rightMargin, 0.4*inch, f"Page {doc.page}")

                canvas_obj.restoreState()

            doc = SimpleDocTemplate(path, pagesize=landscape(A4),
                                    leftMargin=0.6*inch, rightMargin=0.6*inch,
                                    topMargin=0.8*inch, bottomMargin=0.8*inch)
            
            # Use the custom template
            frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
            template = PageTemplate(id='header_footer', frames=[frame], onPage=header_footer)
            doc.addPageTemplates([template])

            styles = getSampleStyleSheet()
            story = []

            # Custom styles
            title_style = ParagraphStyle(
                'Title', parent=styles['Title'],
                fontSize=18, textColor=colors.HexColor(ACCENT),
                alignment=TA_LEFT, spaceAfter=12
            )
            heading2_style = ParagraphStyle(
                'Heading2', parent=styles['Heading2'],
                fontSize=14, textColor=colors.HexColor(ACCENT),
                spaceBefore=16, spaceAfter=8, keepWithNext=True
            )
            heading3_style = ParagraphStyle(
                'Heading3', parent=styles['Heading3'],
                fontSize=12, textColor=colors.HexColor(ACCENT),
                spaceBefore=12, spaceAfter=6, keepWithNext=True
            )

            # ------------------------------------------------------------------
            # SSR Summary Section
            # ------------------------------------------------------------------
            story.append(Paragraph("SSR Detection Summary", heading2_style))
            ssrs = self.state.ssrs
            total = len(ssrs)
            motifs = [s.get("canonical_motif", s.get("motif", "")) for s in ssrs]
            motif_counts = Counter(motifs)

            data = [
                ["Total SSRs", f"{total:,}"],
                ["Unique motifs", f"{len(motif_counts):,}"],
                ["Contigs with SSRs", f"{len(set(s.get('contig') for s in ssrs)):,}"],
            ]
            t = Table(data, colWidths=[2.5*inch, 2.5*inch])
            t.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                ('FONTSIZE', (0,0), (-1,0), 11),
                ('BOTTOMPADDING', (0,0), (-1,0), 8),
                ('BACKGROUND', (0,1), (-1,-1), colors.HexColor("#F5F5F5")),
                ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
            ]))
            story.append(t)
            story.append(Spacer(1, 0.3*inch))

            # Motif distribution
            motif_types = self._get_motif_type_counts(ssrs)
            if motif_types:
                story.append(Paragraph("Motif Type Distribution", heading3_style))
                motif_data = [["Motif Type", "Count", "Percentage"]]
                for k, v in motif_types.items():
                    motif_data.append([k, f"{v:,}", f"{v/total*100:.1f}%"])
                t = Table(motif_data, colWidths=[1.5*inch, 1.5*inch, 1.5*inch])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0,0), (-1,0), 10),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
                ]))
                story.append(t)
                story.append(Spacer(1, 0.2*inch))

            # Top motifs
            top_motifs = motif_counts.most_common(10)
            if top_motifs:
                story.append(Paragraph("Top 10 Canonical Motifs", heading3_style))
                top_data = [["Motif", "Count"]]
                for m, c in top_motifs:
                    top_data.append([m, f"{c:,}"])
                t = Table(top_data, colWidths=[2*inch, 2*inch])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0,0), (-1,0), 10),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
                ]))
                story.append(t)

            story.append(NextPageTemplate('header_footer'))
            story.append(PageBreak())

            # ------------------------------------------------------------------
            # Helper for data tables
            # ------------------------------------------------------------------
            def make_data_table(headers, rows, col_widths, header_color=colors.HexColor(ACCENT)):
                data = [headers] + rows
                t = Table(data, colWidths=col_widths, repeatRows=1)
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), header_color),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0,0), (-1,0), 9),
                    ('BOTTOMPADDING', (0,0), (-1,0), 6),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor("#F9F9F9")]),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#DDDDDD")),
                    ('FONTSIZE', (0,1), (-1,-1), 8),
                ]))
                return t

            # ------------------------------------------------------------------
            # PASS Primers
            # ------------------------------------------------------------------
            pass_primers = self._get_blast_pass_primers()
            if pass_primers:
                story.append(Paragraph("PASS Primers (after BLAST)", heading2_style))
                show_chrom = hasattr(self.state, 'get_display_name') and self.state.chrom_names
                headers = ["SSR ID", "Rank", "Contig"]
                if show_chrom:
                    headers.append("Chromosome")
                headers.extend(["Motif", "Size", "Fwd Tm", "Rev Tm", "Fwd Primer", "Rev Primer"])
                rows = []
                for p in pass_primers:
                    contig = p.get("contig", "")
                    row = [str(p.get("ssr_id", "")), str(p.get("pair_rank", "")), contig[:20]]
                    if show_chrom:
                        row.append(self.state.get_display_name(contig))
                    row.extend([p.get("motif", ""), str(p.get("product_size", "")),
                                f"{p.get('left_tm', 0):.1f}", f"{p.get('right_tm', 0):.1f}",
                                p.get("left_primer", ""), p.get("right_primer", "")])
                    rows.append(row)
                col_widths = [0.5*inch, 0.35*inch, 0.9*inch]
                if show_chrom:
                    col_widths.append(0.7*inch)
                col_widths.extend([0.5*inch, 0.35*inch, 0.45*inch, 0.45*inch, 1.6*inch, 1.6*inch])
                story.append(make_data_table(headers, rows, col_widths))
                story.append(NextPageTemplate('header_footer'))
                story.append(PageBreak())

            # ------------------------------------------------------------------
            # Amplicon Section
            # ------------------------------------------------------------------
            if self.state.amplicon_validation_result:
                result = self.state.amplicon_validation_result
                story.append(Paragraph("Amplicon Panel Validation", heading2_style))
                clean_pool = result.get("clean_pool", [])
                issues = result.get("issues_summary", {})
                amplicon_data = [
                    ["Metric", "Value"],
                    ["Input loci", str(result.get("n_input", 0))],
                    ["Passed checks", str(result.get("n_clean", 0))],
                    ["Final pool", str(result.get("n_final", len(clean_pool)))],
                    ["Dimer risks", str(issues.get("dimer_risks", 0))],
                    ["Size filtered", str(issues.get("size_filtered_short", 0) + issues.get("size_filtered_long", 0))],
                    ["Compatibility score", f"{result.get('compatibility_score', 0):.2f}"],
                ]
                t = Table(amplicon_data, colWidths=[2*inch, 2*inch])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('BOTTOMPADDING', (0,0), (-1,0), 8),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
                ]))
                story.append(t)
                story.append(Spacer(1, 0.2*inch))

                if clean_pool:
                    story.append(Paragraph("Clean Pool Primers", heading3_style))
                    show_chrom = hasattr(self.state, 'get_display_name') and self.state.chrom_names
                    headers = ["SSR ID", "Contig"]
                    if show_chrom:
                        headers.append("Chromosome")
                    headers.extend(["Motif", "Size", "Fwd Tm", "Rev Tm", "Fwd Primer", "Rev Primer"])
                    rows = []
                    for p in clean_pool:
                        contig = p.get("contig", "")
                        row = [str(p.get("ssr_id", "")), contig[:20]]
                        if show_chrom:
                            row.append(self.state.get_display_name(contig))
                        row.extend([p.get("motif", ""), str(p.get("product_size", "")),
                                    f"{p.get('left_tm', 0):.1f}", f"{p.get('right_tm', 0):.1f}",
                                    p.get("left_primer", ""), p.get("right_primer", "")])
                        rows.append(row)
                    col_widths = [0.5*inch, 0.9*inch]
                    if show_chrom:
                        col_widths.append(0.7*inch)
                    col_widths.extend([0.5*inch, 0.35*inch, 0.45*inch, 0.45*inch, 1.6*inch, 1.6*inch])
                    story.append(make_data_table(headers, rows, col_widths))
                story.append(NextPageTemplate('header_footer'))
                story.append(PageBreak())

            # ------------------------------------------------------------------
            # Capillary Section
            # ------------------------------------------------------------------
            if self.state.capillary_result:
                result = self.state.capillary_result
                story.append(Paragraph("Capillary Panel Results", heading2_style))
                assignments = result.get("assignments", [])
                cap_data = [
                    ["Metric", "Value"],
                    ["Total loci", str(result.get("n_total", 0))],
                    ["Assigned", str(result.get("n_assigned", len(assignments)))],
                    ["Unassigned", str(len(result.get("unassigned", [])))],
                    ["Dyes used", str(len(set(a.get("dye", "") for a in assignments)))],
                ]
                t = Table(cap_data, colWidths=[2*inch, 2*inch])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('BOTTOMPADDING', (0,0), (-1,0), 8),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
                ]))
                story.append(t)
                story.append(Spacer(1, 0.2*inch))

                if assignments:
                    story.append(Paragraph("Dye Assignments", heading3_style))
                    show_chrom = hasattr(self.state, 'get_display_name') and self.state.chrom_names
                    has_panels = any("panel" in a for a in assignments)
                    headers = []
                    if has_panels:
                        headers.append("Panel")
                    headers.extend(["SSR ID", "Dye", "Min", "Max", "Range", "Motif", "Contig"])
                    if show_chrom:
                        headers.append("Chromosome")
                    headers.extend(["Fwd Primer", "Rev Primer"])
                    rows = []
                    for a in assignments:
                        contig = a.get("contig", "")
                        row = []
                        if has_panels:
                            row.append(str(a.get("panel", "")))
                        row.extend([
                            str(a.get("ssr_id", "")),
                            a.get("dye", ""),
                            str(a.get("min_size", "")),
                            str(a.get("max_size", "")),
                            str(a.get("max_size", 0) - a.get("min_size", 0)),
                            a.get("motif", ""),
                            contig[:20],
                        ])
                        if show_chrom:
                            row.append(self.state.get_display_name(contig))
                        row.extend([a.get("left_primer", ""), a.get("right_primer", "")])
                        rows.append(row)
                    col_widths = []
                    if has_panels:
                        col_widths.append(0.4*inch)
                    col_widths.extend([0.5*inch, 0.5*inch, 0.35*inch, 0.35*inch, 0.35*inch, 0.5*inch, 0.9*inch])
                    if show_chrom:
                        col_widths.append(0.7*inch)
                    col_widths.extend([1.5*inch, 1.5*inch])
                    story.append(make_data_table(headers, rows, col_widths))
                story.append(NextPageTemplate('header_footer'))
                story.append(PageBreak())

            # ------------------------------------------------------------------
            # GBS-RE Section
            # ------------------------------------------------------------------
            if self.state.gbs_re_markers:
                markers = self.state.gbs_re_markers
                story.append(Paragraph("GBS‑RE Marker Discovery", heading2_style))
                gbs_data = [
                    ["Metric", "Value"],
                    ["Qualified SSRs", f"{len(markers):,}"],
                    ["Passing fragments", f"{self.state.gbs_re_passing_frags:,}"],
                    ["Unique contigs", f"{len(set(m.get('contig') for m in markers)):,}"],
                ]
                t = Table(gbs_data, colWidths=[2*inch, 2*inch])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.HexColor(ACCENT)),
                    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
                    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                    ('BOTTOMPADDING', (0,0), (-1,0), 8),
                    ('BACKGROUND', (0,1), (-1,-1), colors.white),
                    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CCCCCC")),
                ]))
                story.append(t)
                story.append(Spacer(1, 0.2*inch))

                if markers:
                    story.append(Paragraph("Qualified Markers", heading3_style))
                    show_chrom = hasattr(self.state, 'get_display_name') and self.state.chrom_names
                    headers = ["SSR ID", "Contig"]
                    if show_chrom:
                        headers.append("Chromosome")
                    headers.extend(["Fragment Start", "Fragment End", "Size", "Motif"])
                    rows = []
                    for m in markers:
                        contig = m.get("contig", "")
                        row = [str(m.get("ssr_id", "")), contig[:20]]
                        if show_chrom:
                            row.append(self.state.get_display_name(contig))
                        row.extend([
                            str(m.get("fragment_start", "")),
                            str(m.get("fragment_end", "")),
                            str(m.get("fragment_size", "")),
                            m.get("motif", ""),
                        ])
                        rows.append(row)
                    col_widths = [0.5*inch, 1.0*inch]
                    if show_chrom:
                        col_widths.append(0.7*inch)
                    col_widths.extend([0.7*inch, 0.7*inch, 0.5*inch, 0.7*inch])
                    story.append(make_data_table(headers, rows, col_widths))

            doc.build(story)
            self._status_label.setText(f"PDF report saved to {path}")
            self.mw.set_status("PDF report generated")

        except Exception as e:
            QMessageBox.critical(self, "PDF Generation Failed", str(e))

    def on_show(self):
        self._refresh()