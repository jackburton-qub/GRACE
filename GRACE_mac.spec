# GRACE_mac.spec
# --------------
# PyInstaller build spec for GRACE (PyQt6 version) — macOS build.
#
# Produces dist/GRACE/GRACE.app
#
# Build with:
#   python3.12 build.py
#
# Or directly:
#   pyinstaller GRACE_mac.spec
#
# Before building, ensure BLAST+ binaries are in blast/mac/:
#   blast/mac/blastn
#   blast/mac/makeblastdb

import sys
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules, collect_data_files

# ---------------------------------------------------------------------------
# Collect packages
# ---------------------------------------------------------------------------

qt_datas,       qt_binaries,       qt_hiddenimports       = collect_all("PyQt6")
primer3_datas,  primer3_binaries,  primer3_hiddenimports  = collect_all("primer3")
reportlab_datas,reportlab_binaries,reportlab_hiddenimports= collect_all("reportlab")

pandas_hiddenimports = collect_submodules("pandas")
numpy_hiddenimports  = collect_submodules("numpy")

# ---------------------------------------------------------------------------
# Application source files
# ---------------------------------------------------------------------------

app_datas = [
    ("main.py",      "."),
    ("app_state.py", "."),
    ("ui",           "ui"),
    ("core",         "core"),
    ("assets",       "assets"),
    ("blast/mac",    "blast/mac"),
]

# ---------------------------------------------------------------------------
# Hidden imports
# ---------------------------------------------------------------------------

all_hidden = (
    qt_hiddenimports
    + primer3_hiddenimports
    + reportlab_hiddenimports
    + pandas_hiddenimports
    + numpy_hiddenimports
    + [
        "multiprocessing",
        "multiprocessing.pool",
        "multiprocessing.managers",
        "multiprocessing.reduction",
        "multiprocessing.resource_tracker",
        "multiprocessing.popen_fork",
        "multiprocessing.popen_forkserver",
        "multiprocessing.popen_spawn_posix",
        "multiprocessing.popen_spawn_win32",
        "subprocess",
        "tempfile",
        "socket",
        "threading",
        "dataclasses",
        "re",
        "io",
        "json",
        "copy",
        "shutil",
        "itertools",
        "functools",
        "bisect",
        "collections",
        "reportlab",
        "reportlab.lib",
        "reportlab.lib.pagesizes",
        "reportlab.lib.styles",
        "reportlab.lib.units",
        "reportlab.lib.colors",
        "reportlab.platypus",
        "reportlab.pdfgen",
    ]
)

# ---------------------------------------------------------------------------
# Combined datas and binaries
# ---------------------------------------------------------------------------

all_datas = (
    app_datas
    + qt_datas
    + primer3_datas
    + reportlab_datas
)

all_binaries = (
    qt_binaries
    + primer3_binaries
    + reportlab_binaries
)

# ---------------------------------------------------------------------------
# ANALYSIS
# ---------------------------------------------------------------------------

a = Analysis(
    ["main.py"],
    pathex=["."],
    binaries=all_binaries,
    datas=all_datas,
    hiddenimports=all_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "matplotlib",
        "IPython",
        "jupyter",
        "notebook",
        "scipy",
        "sklearn",
        "tensorflow",
        "torch",
        "streamlit",
        "webview",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GRACE",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    windowed=True,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="GRACE",
)

app = BUNDLE(
    coll,
    name="GRACE.app",
    icon=None,
    bundle_identifier="com.grace.app",
    info_plist={
        "NSPrincipalClass":               "NSApplication",
        "NSHighResolutionCapable":        True,
        "NSRequiresAquaSystemAppearance": False,
        "CFBundleShortVersionString":     "1.0.0",
        "CFBundleName":                   "GRACE",
        "CFBundleDisplayName":            "GRACE",
    },
)
