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
# (chmod +x both files on Mac before building)

import sys
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules, collect_data_files

# ---------------------------------------------------------------------------
# Collect packages PyInstaller cannot fully analyse automatically
# ---------------------------------------------------------------------------

qt_datas, qt_binaries, qt_hiddenimports = collect_all("PyQt6")
primer3_datas, primer3_binaries, primer3_hiddenimports = collect_all("primer3")
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
    ("blast/mac",    "blast/mac"),   # bundled BLAST+ binaries for macOS
]

# ---------------------------------------------------------------------------
# Hidden imports
# ---------------------------------------------------------------------------

all_hidden = (
    qt_hiddenimports
    + primer3_hiddenimports
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
    ]
)

# ---------------------------------------------------------------------------
# Combined datas and binaries
# ---------------------------------------------------------------------------

all_datas    = app_datas + qt_datas + primer3_datas
all_binaries = qt_binaries + primer3_binaries

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
    upx=False,           # UPX disabled — corrupts binaries on Mac too
    console=False,
    windowed=True,
    icon="grace_icon.icns",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,           # UPX disabled
    upx_exclude=[],
    name="GRACE",
)

# Mac .app bundle
app = BUNDLE(
    coll,
    name="GRACE.app",
    icon="grace_icon.icns",
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
