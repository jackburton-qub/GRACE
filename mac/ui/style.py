"""
style.py
--------
GRACE colour palette and Qt stylesheet.
IBM Plex Sans + IBM Plex Mono — GitHub dark palette with Irish green accent.
"""

# ── Palette ───────────────────────────────────────────────
BG_DARK    = "#0d1117"   # main window background
BG_MID     = "#161b22"   # panel / sidebar background
BG_LIGHT   = "#21262d"   # card / input / button background
BG_HOVER   = "#30363d"   # hover state

ACCENT     = "#2ea043"   # primary green accent
ACCENT_DIM = "#238636"   # darker green for borders / pressed
ACCENT_BRIGHT = "#3fb950" # bright green for highlights

SUCCESS    = "#3fb950"   # PASS
ERROR      = "#f85149"   # FAIL
WARNING    = "#d29922"   # amber
MUTED      = "#484f58"   # very muted

TEXT_PRIMARY   = "#e6edf3"   # main text
TEXT_SECONDARY = "#8b949e"   # labels / help text
TEXT_MONO      = "#79c0ff"   # monospace data values

BORDER        = "#21262d"    # standard border
BORDER_ACTIVE = "#2ea043"    # focused / active border

# ── Fonts ─────────────────────────────────────────────────
FONT_UI   = "IBM Plex Sans"
FONT_MONO = "IBM Plex Mono"
FONT_SIZE_NORMAL = 9
FONT_SIZE_SMALL  = 8
FONT_SIZE_LARGE  = 11
FONT_SIZE_TITLE  = 14

# ── Dimensions ───────────────────────────────────────────
SIDEBAR_WIDTH  = 210
WINDOW_MIN_W   = 1200
WINDOW_MIN_H   = 800
PANEL_PADDING  = 24
CARD_RADIUS    = 6
BUTTON_RADIUS  = 6
BUTTON_HEIGHT  = 32

# ── Stylesheet ────────────────────────────────────────────
STYLESHEET = f"""
/* ── Global ── */
QMainWindow, QWidget {{
    background-color: {BG_DARK};
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", "Segoe UI", system-ui, sans-serif;
    font-size: {FONT_SIZE_NORMAL}pt;
}}

/* ── Sidebar ── */
QListWidget#sidebar {{
    background-color: {BG_MID};
    border: none;
    border-right: 1px solid {BORDER};
    outline: none;
    padding: 4px 0;
}}
QListWidget#sidebar::item {{
    padding: 10px 16px 10px 20px;
    border-left: 3px solid transparent;
    color: {TEXT_SECONDARY};
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_NORMAL}pt;
}}
QListWidget#sidebar::item:selected {{
    background-color: {BG_LIGHT};
    border-left: 3px solid {ACCENT};
    color: {TEXT_PRIMARY};
}}
QListWidget#sidebar::item:hover:!selected {{
    background-color: {BG_HOVER};
    color: {TEXT_PRIMARY};
    border-left: 3px solid {MUTED};
}}

/* ── Buttons ── */
QPushButton {{
    background-color: {BG_LIGHT};
    color: {TEXT_PRIMARY};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 16px;
    min-height: {BUTTON_HEIGHT}px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_NORMAL}pt;
    font-weight: 500;
}}
QPushButton:hover {{
    background-color: {BG_HOVER};
    border-color: {ACCENT};
    color: {TEXT_PRIMARY};
}}
QPushButton:pressed {{
    background-color: {ACCENT_DIM};
    border-color: {ACCENT_DIM};
}}
QPushButton:disabled {{
    color: {MUTED};
    border-color: {BORDER};
    background-color: {BG_MID};
}}
QPushButton#primary {{
    background-color: {ACCENT_DIM};
    border-color: {ACCENT_DIM};
    color: #ffffff;
    font-weight: 600;
}}
QPushButton#primary:hover {{
    background-color: {ACCENT};
    border-color: {ACCENT};
}}
QPushButton#primary:pressed {{
    background-color: #196127;
    border-color: #196127;
}}
QPushButton#success {{
    background-color: #0f2d16;
    border-color: {SUCCESS};
    color: {SUCCESS};
    font-weight: 500;
}}
QPushButton#success:hover {{
    background-color: #1a4d26;
}}

/* ── Inputs ── */
QLineEdit, QTextEdit {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 10px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    selection-background-color: {ACCENT_DIM};
    min-height: 28px;
}}
QLineEdit:focus, QTextEdit:focus {{
    border-color: {ACCENT};
    outline: none;
}}
QSpinBox, QDoubleSpinBox {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 4px 8px;
    padding-right: 20px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    min-height: 28px;
    min-width: 80px;
}}
QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {ACCENT};
}}
QSpinBox::up-button, QDoubleSpinBox::up-button {{
    subcontrol-origin: border;
    subcontrol-position: top right;
    width: 18px;
    border-left: 1px solid {BORDER};
    background: {BG_HOVER};
    border-top-right-radius: {BUTTON_RADIUS}px;
}}
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    subcontrol-origin: border;
    subcontrol-position: bottom right;
    width: 18px;
    border-left: 1px solid {BORDER};
    background: {BG_HOVER};
    border-bottom-right-radius: {BUTTON_RADIUS}px;
}}
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {{
    background: {ACCENT_DIM};
}}
QComboBox {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 10px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    min-height: 28px;
    selection-background-color: {ACCENT_DIM};
}}
QComboBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{
    border: none;
    width: 24px;
}}
QComboBox QAbstractItemView {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: 0px;
    selection-background-color: {ACCENT_DIM};
    color: {TEXT_PRIMARY};
    padding: 2px;
}}

/* ── Labels ── */
QLabel {{
    font-family: "{FONT_UI}", sans-serif;
    color: {TEXT_PRIMARY};
    background: transparent;
}}

/* ── Group box ── */
QGroupBox {{
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    margin-top: 14px;
    padding: 12px 16px 16px 16px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_SMALL}pt;
    font-weight: 600;
    color: {TEXT_SECONDARY};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    padding: 0 6px;
    color: {TEXT_SECONDARY};
    background-color: {BG_DARK};
    font-size: {FONT_SIZE_SMALL}pt;
    font-weight: 600;
    letter-spacing: 0.5px;
    text-transform: uppercase;
}}

/* ── Tables ── */
QTableWidget, QTableView {{
    background-color: {BG_MID};
    alternate-background-color: {BG_DARK};
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    gridline-color: {BORDER};
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_SMALL}pt;
    selection-background-color: #1f3d2a;
    selection-color: {TEXT_PRIMARY};
}}
QTableWidget::item, QTableView::item {{
    padding: 5px 10px;
    border: none;
}}
QTableWidget::item:selected, QTableView::item:selected {{
    background-color: #1f3d2a;
    color: {TEXT_PRIMARY};
}}
QHeaderView {{
    background-color: {BG_MID};
}}
QHeaderView::section {{
    background-color: {BG_LIGHT};
    color: {TEXT_SECONDARY};
    border: none;
    border-right: 1px solid {BORDER};
    border-bottom: 1px solid {BORDER};
    padding: 6px 10px;
    font-family: "{FONT_UI}", sans-serif;
    font-weight: 600;
    font-size: {FONT_SIZE_SMALL}pt;
    letter-spacing: 0.3px;
}}
QHeaderView::section:first {{
    border-top-left-radius: {CARD_RADIUS}px;
}}

/* ── Progress bar ── */
QProgressBar {{
    background-color: {BG_LIGHT};
    border: none;
    border-radius: 3px;
    height: 6px;
    text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{
    background-color: {ACCENT};
    border-radius: 3px;
}}

/* ── Scrollbars ── */
QScrollBar:vertical {{
    background: {BG_MID};
    width: 8px;
    border: none;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {BG_HOVER};
    border-radius: 4px;
    min-height: 24px;
    margin: 2px;
}}
QScrollBar::handle:vertical:hover {{
    background: {MUTED};
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
QScrollBar:horizontal {{
    background: {BG_MID};
    height: 8px;
    border: none;
    margin: 0;
}}
QScrollBar::handle:horizontal {{
    background: {BG_HOVER};
    border-radius: 4px;
    min-width: 24px;
    margin: 2px;
}}
QScrollBar::handle:horizontal:hover {{
    background: {MUTED};
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{ width: 0; }}

/* ── Tab widget ── */
QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    background: {BG_MID};
    top: -1px;
}}
QTabBar::tab {{
    background: {BG_DARK};
    color: {TEXT_SECONDARY};
    padding: 7px 18px;
    border: 1px solid {BORDER};
    border-bottom: none;
    border-top-left-radius: {BUTTON_RADIUS}px;
    border-top-right-radius: {BUTTON_RADIUS}px;
    margin-right: 2px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_SMALL}pt;
    font-weight: 500;
}}
QTabBar::tab:selected {{
    background: {BG_MID};
    color: {TEXT_PRIMARY};
    border-bottom: 2px solid {ACCENT};
}}
QTabBar::tab:hover:!selected {{
    background: {BG_LIGHT};
    color: {TEXT_PRIMARY};
}}

/* ── Checkbox ── */
QCheckBox {{
    color: {TEXT_PRIMARY};
    spacing: 8px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_NORMAL}pt;
}}
QCheckBox::indicator {{
    width: 15px;
    height: 15px;
    border: 1px solid {BORDER};
    border-radius: 3px;
    background: {BG_LIGHT};
}}
QCheckBox::indicator:checked {{
    background-color: {ACCENT_DIM};
    border-color: {ACCENT};
    image: none;
}}
QCheckBox::indicator:hover {{
    border-color: {ACCENT};
}}

/* ── Status bar ── */
QStatusBar {{
    background-color: {BG_MID};
    border-top: 1px solid {BORDER};
    color: {TEXT_SECONDARY};
    font-family: "{FONT_MONO}", monospace;
    font-size: {FONT_SIZE_SMALL}pt;
    padding: 3px 10px;
}}

/* ── Tooltip ── */
QToolTip {{
    background-color: {BG_LIGHT};
    color: {TEXT_PRIMARY};
    border: 1px solid {BORDER};
    padding: 6px 10px;
    border-radius: {BUTTON_RADIUS}px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {FONT_SIZE_SMALL}pt;
}}

/* ── Scroll area ── */
QScrollArea {{
    border: none;
    background: transparent;
}}
QScrollArea > QWidget > QWidget {{
    background: transparent;
}}

/* ── Message box ── */
QMessageBox {{
    background-color: {BG_MID};
    color: {TEXT_PRIMARY};
}}
QMessageBox QPushButton {{
    min-width: 80px;
}}

/* ── File dialog ── */
QFileDialog {{
    background-color: {BG_MID};
    color: {TEXT_PRIMARY};
}}

/* ── Splitter ── */
QSplitter::handle {{
    background: {BORDER};
}}
QSplitter::handle:hover {{
    background: {ACCENT_DIM};
}}
"""