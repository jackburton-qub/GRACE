"""
main_window.py
--------------
QMainWindow — sidebar navigation + stacked panel area.
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QListWidgetItem, QStackedWidget,
    QLabel, QStatusBar, QSizePolicy, QFrame,
    QButtonGroup, QRadioButton,
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QFont

from ui.style import (
    STYLESHEET, build_stylesheet,
    SIDEBAR_WIDTH, WINDOW_MIN_W, WINDOW_MIN_H,
    ACCENT, SUCCESS, MUTED, TEXT_SECONDARY,
    BG_MID, BG_LIGHT, BG_DARK, BORDER, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_SMALL,
    FONT_PRESETS, FONT_SIZE_PRESET_DEFAULT,
)
from ui.panels.home_panel        import HomePanel
from ui.panels.ssr_panel         import SSRPanel
from ui.panels.ssr_summary_panel import SSRSummaryPanel
from ui.panels.primer_panel      import PrimerPanel
from ui.panels.specificity_panel import SpecificityPanel
from ui.panels.results_panel     import ResultsPanel
from ui.panels.report_panel      import ReportPanel
from ui.panels.amplicon_panel    import AmpliconPanel
from ui.panels.capillary_panel   import CapillaryPanel
from ui.panels.gbs_re_panel      import GBSREPanel


STEPS = [
    (0,  "01  Load Genome"),
    (1,  "02  Detect SSRs"),
    (2,  "03  SSR Summary"),
    (3,  "04  Design Primers"),
    (4,  "05  Check Specificity"),
    (5,  "06  BLAST Results"),
    (6,  "07  Final Report"),
    (7,  "08  Amplicon Tools"),
    (8,  "09  Capillary Tools"),
    (9,  "10  GBS‑RE Tools"),
]


class MainWindow(QMainWindow):
    def __init__(self, state):
        super().__init__()
        self.state = state
        self.setWindowTitle("GRACE — Genomic Repeat Analysis and Characterisation Engine")
        # Increased minimum height to accommodate 10 sidebar items
        self.setMinimumSize(WINDOW_MIN_W, max(WINDOW_MIN_H, 800))
        self.setStyleSheet(STYLESHEET)
        self._build_ui()
        self._connect_signals()
        self.refresh_sidebar()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Sidebar ──────────────────────────────────────
        sidebar_container = QWidget()
        sidebar_container.setFixedWidth(SIDEBAR_WIDTH)
        sidebar_container.setStyleSheet(
            f"background-color: {BG_MID}; border-right: 1px solid {BORDER};"
        )
        sidebar_layout = QVBoxLayout(sidebar_container)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        sidebar_layout.setSpacing(0)

        header = QWidget()
        header.setFixedHeight(64)
        header.setStyleSheet(
            f"background-color: {BG_MID}; border-bottom: 1px solid {BORDER}; padding: 0;"
        )
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)
        header_layout.setSpacing(2)

        title_label = QLabel("GRACE")
        title_label.setFont(QFont(FONT_UI, 13, QFont.Weight.Bold))
        title_label.setStyleSheet(f"color: {ACCENT};")

        sub_label = QLabel("Genomic Repeat Analysis")
        sub_label.setFont(QFont(FONT_MONO, 7))
        sub_label.setStyleSheet(f"color: {TEXT_SECONDARY};")

        header_layout.addWidget(title_label)
        header_layout.addWidget(sub_label)

        self.sidebar = QListWidget()
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.sidebar.setSpacing(2)

        for _, label in STEPS:
            item = QListWidgetItem(label)
            item.setSizeHint(QSize(SIDEBAR_WIDTH, 48))
            item.setFont(QFont(FONT_UI, FONT_SIZE_NORMAL))
            self.sidebar.addItem(item)

        # Font size selector
        font_size_widget = QWidget()
        font_size_widget.setStyleSheet(
            f"background: transparent; border-top: 1px solid {BORDER}; padding: 8px 0 4px 0;"
        )
        font_size_layout = QVBoxLayout(font_size_widget)
        font_size_layout.setContentsMargins(16, 8, 16, 4)
        font_size_layout.setSpacing(6)

        font_size_label = QLabel("TEXT SIZE")
        font_size_label.setStyleSheet(
            f"color: {MUTED}; font-size: 7pt; font-weight: 600; "
            f"letter-spacing: 0.5px; border: none;"
        )
        font_size_layout.addWidget(font_size_label)

        btn_row = QWidget()
        btn_row.setStyleSheet("background: transparent; border: none;")
        btn_row_layout = QHBoxLayout(btn_row)
        btn_row_layout.setContentsMargins(0, 0, 0, 0)
        btn_row_layout.setSpacing(4)

        self._font_btn_group = QButtonGroup(self)
        self._font_preset = FONT_SIZE_PRESET_DEFAULT

        for preset in ("Small", "Medium", "Large"):
            btn = QRadioButton(preset)
            btn.setStyleSheet(f"""
                QRadioButton {{
                    color: {TEXT_SECONDARY}; font-size: 8pt;
                    spacing: 4px; border: none; background: transparent;
                }}
                QRadioButton:checked {{ color: {ACCENT}; }}
                QRadioButton::indicator {{
                    width: 10px; height: 10px;
                    border: 1px solid {MUTED}; border-radius: 5px;
                    background: transparent;
                }}
                QRadioButton::indicator:checked {{
                    background: {ACCENT}; border-color: {ACCENT};
                }}
            """)
            if preset == FONT_SIZE_PRESET_DEFAULT:
                btn.setChecked(True)
            self._font_btn_group.addButton(btn)
            btn_row_layout.addWidget(btn)

        font_size_layout.addWidget(btn_row)

        version_label = QLabel("v1.0.0")
        version_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        version_label.setStyleSheet(f"color: {MUTED}; font-size: 7pt; padding: 8px;")

        sidebar_layout.addWidget(header)
        sidebar_layout.addWidget(self.sidebar)
        sidebar_layout.addStretch()
        sidebar_layout.addWidget(font_size_widget)
        sidebar_layout.addWidget(version_label)

        # ── Main panel stack ─────────────────────────────
        self.stack = QStackedWidget()

        self.home_panel        = HomePanel(self.state, self)
        self.ssr_panel         = SSRPanel(self.state, self)
        self.ssr_summary_panel = SSRSummaryPanel(self.state, self)
        self.primer_panel      = PrimerPanel(self.state, self)
        self.specificity_panel = SpecificityPanel(self.state, self)
        self.results_panel     = ResultsPanel(self.state, self)
        self.report_panel      = ReportPanel(self.state, self)
        self.amplicon_panel    = AmpliconPanel(self.state, self)
        self.capillary_panel   = CapillaryPanel(self.state, self)
        self.gbs_re_panel      = GBSREPanel(self.state, self)

        for panel in [
            self.home_panel, self.ssr_panel, self.ssr_summary_panel,
            self.primer_panel, self.specificity_panel, self.results_panel,
            self.report_panel, self.amplicon_panel, self.capillary_panel,
            self.gbs_re_panel,
        ]:
            self.stack.addWidget(panel)

        root.addWidget(sidebar_container)
        root.addWidget(self.stack, 1)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    def _connect_signals(self):
        self.sidebar.currentRowChanged.connect(self._on_tab_changed)
        self.sidebar.setCurrentRow(0)
        self._font_btn_group.buttonClicked.connect(self._on_font_size_changed)

    def _on_tab_changed(self, index):
        self.stack.setCurrentIndex(index)
        panel = self.stack.currentWidget()
        if hasattr(panel, "on_show"):
            panel.on_show()

    def _on_font_size_changed(self, btn):
        preset = btn.text()
        if preset == self._font_preset:
            return
        self._font_preset = preset
        self.setStyleSheet(build_stylesheet(preset))

    def set_status(self, message: str):
        self.status_bar.showMessage(message)

    def refresh_sidebar(self):
        from PyQt6.QtGui import QColor
        s = self.state
        mode = s.workflow_mode

        is_amplicon = (mode == "amplicon")
        is_capillary = (mode == "capillary")

        statuses = [
            s.has_genome,       # 0 Load Genome
            s.has_ssrs,         # 1 Detect SSRs
            s.has_ssrs,         # 2 SSR Summary
            s.has_primers,      # 3 Design Primers
            s.has_specificity,  # 4 Check Specificity
            s.has_specificity,  # 5 BLAST Results
            s.has_specificity,  # 6 Final Report
            s.has_specificity and is_amplicon,   # 7 Amplicon Tools
            s.has_specificity and is_capillary,  # 8 Capillary Tools
            s.has_genome and s.has_ssrs,         # 9 GBS‑RE Tools
        ]

        for i, (_, label) in enumerate(STEPS):
            item = self.sidebar.item(i)
            if item is None:
                continue
            is_current = self.sidebar.currentRow() == i
            done = statuses[i]

            # Dim inactive tool steps
            if i == 7 and not is_amplicon:
                item.setForeground(QColor(MUTED))
            elif i == 8 and not is_capillary:
                item.setForeground(QColor(MUTED))
            elif done:
                item.setForeground(QColor(SUCCESS))
            elif is_current:
                item.setForeground(QColor(ACCENT))
            else:
                item.setForeground(QColor(TEXT_SECONDARY))

    def navigate_to(self, step: int):
        self.sidebar.setCurrentRow(step)

    def on_step_complete(self, step: int):
        self.refresh_sidebar()
        self.set_status(f"Step {step + 1} complete")