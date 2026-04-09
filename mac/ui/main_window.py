"""
main_window.py
--------------
QMainWindow — sidebar navigation + stacked panel area.
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QListWidgetItem, QStackedWidget,
    QLabel, QStatusBar, QSizePolicy, QFrame,
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QFont

from ui.style import (
    STYLESHEET, SIDEBAR_WIDTH, WINDOW_MIN_W, WINDOW_MIN_H,
    ACCENT, SUCCESS, MUTED, TEXT_SECONDARY,
    BG_MID, BG_LIGHT, BG_DARK, BORDER, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_SMALL,
)
from ui.panels.home_panel        import HomePanel
from ui.panels.ssr_panel         import SSRPanel
from ui.panels.primer_panel      import PrimerPanel
from ui.panels.specificity_panel import SpecificityPanel
from ui.panels.results_panel     import ResultsPanel


STEPS = [
    (0, "01  Load Genome"),
    (1, "02  Detect SSRs"),
    (2, "03  Design Primers"),
    (3, "04  Check Specificity"),
    (4, "05  BLAST Results"),
]


class MainWindow(QMainWindow):
    def __init__(self, state):
        super().__init__()
        self.state = state
        self.setWindowTitle("GRACE — Genomic Repeat Analysis and Characterisation Engine")
        self.setMinimumSize(WINDOW_MIN_W, WINDOW_MIN_H)
        self.setStyleSheet(STYLESHEET)

        self._build_ui()
        self._connect_signals()
        self.refresh_sidebar()

    # ---------------------------------------------------------
    # UI BUILD
    # ---------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Sidebar ──────────────────────────────────────
        sidebar_container = QWidget()
        sidebar_container.setFixedWidth(SIDEBAR_WIDTH)
        sidebar_container.setStyleSheet(f"background-color: {BG_MID}; border-right: 1px solid {BORDER};")
        sidebar_layout = QVBoxLayout(sidebar_container)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        sidebar_layout.setSpacing(0)

        # App name header in sidebar
        header = QWidget()
        header.setFixedHeight(64)
        header.setStyleSheet(f"background-color: {BG_MID}; border-bottom: 1px solid {BORDER}; padding: 0;")
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)
        header_layout.setSpacing(2)

        title_label = QLabel("GRACE")
        title_font = QFont(FONT_UI, 13, QFont.Weight.Bold)
        title_label.setFont(title_font)
        title_label.setStyleSheet(f"color: {ACCENT};")

        sub_label = QLabel("Genomic Repeat Analysis")
        sub_font = QFont(FONT_MONO, 7)
        sub_label.setFont(sub_font)
        sub_label.setStyleSheet(f"color: {TEXT_SECONDARY};")

        header_layout.addWidget(title_label)
        header_layout.addWidget(sub_label)

        # Step list
        self.sidebar = QListWidget()
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.sidebar.setSpacing(2)

        for _, label in STEPS:
            item = QListWidgetItem(label)
            item.setSizeHint(QSize(SIDEBAR_WIDTH, 52))
            font = QFont(FONT_UI, FONT_SIZE_NORMAL)
            item.setFont(font)
            self.sidebar.addItem(item)

        # Version label at bottom of sidebar
        version_label = QLabel("v1.0.0")
        version_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        version_label.setStyleSheet(f"color: {MUTED}; font-size: 7pt; padding: 8px;")

        sidebar_layout.addWidget(header)
        sidebar_layout.addWidget(self.sidebar)
        sidebar_layout.addStretch()
        sidebar_layout.addWidget(version_label)

        # ── Main panel stack ─────────────────────────────
        self.stack = QStackedWidget()

        self.home_panel        = HomePanel(self.state, self)
        self.ssr_panel         = SSRPanel(self.state, self)
        self.primer_panel      = PrimerPanel(self.state, self)
        self.specificity_panel = SpecificityPanel(self.state, self)
        self.results_panel     = ResultsPanel(self.state, self)

        self.stack.addWidget(self.home_panel)
        self.stack.addWidget(self.ssr_panel)
        self.stack.addWidget(self.primer_panel)
        self.stack.addWidget(self.specificity_panel)
        self.stack.addWidget(self.results_panel)

        root.addWidget(sidebar_container)
        root.addWidget(self.stack, 1)

        # ── Status bar ───────────────────────────────────
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    # ---------------------------------------------------------
    # SIGNALS
    # ---------------------------------------------------------
    def _connect_signals(self):
        self.sidebar.currentRowChanged.connect(self._on_tab_changed)
        self.sidebar.setCurrentRow(0)

    def _on_tab_changed(self, index):
        self.stack.setCurrentIndex(index)
        # Notify panel it's being shown — lets panels refresh their display
        panel = self.stack.currentWidget()
        if hasattr(panel, "on_show"):
            panel.on_show()

    # ---------------------------------------------------------
    # PUBLIC API — called by panels to update UI
    # ---------------------------------------------------------
    def set_status(self, message: str):
        self.status_bar.showMessage(message)

    def refresh_sidebar(self):
        """Update sidebar item colours to reflect current data state."""
        from PyQt6.QtGui import QColor
        s = self.state

        statuses = [
            s.has_genome,
            s.has_ssrs,
            s.has_primers,
            s.has_specificity,
            s.has_specificity,
        ]

        for i, (_, label) in enumerate(STEPS):
            item = self.sidebar.item(i)
            is_current = self.sidebar.currentRow() == i
            done = statuses[i]

            if done:
                item.setForeground(QColor(SUCCESS))
            elif is_current:
                item.setForeground(QColor(ACCENT))
            else:
                item.setForeground(QColor(TEXT_SECONDARY))

    def navigate_to(self, step: int):
        """Navigate to a specific step (0-indexed)."""
        self.sidebar.setCurrentRow(step)

    def on_step_complete(self, step: int):
        """Called by panels when a step completes — refreshes sidebar and advances."""
        self.refresh_sidebar()
        self.set_status(f"Step {step + 1} complete")