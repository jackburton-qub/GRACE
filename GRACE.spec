# GRACE.spec
# ----------
# PyInstaller build spec for GRACE (PyQt6 version) — Windows build.
#
# Produces a single folder in dist/GRACE/ containing GRACE.exe
#
# Build with:
#   python build.py
#
# Or directly:
#   pyinstaller GRACE.spec
#
# Before building, ensure BLAST+ binaries are in blast/windows/:
#   blast/windows/blastn.exe
#   blast/windows/blastn.exe.manifest
#   blast/windows/makeblastdb.exe
#   blast/windows/makeblastdb.exe.manifest
#   blast/windows/ncbi-vdb-md.dll
#   blast/windows/nghttp2.dll

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
    ("main.py",          "."),
    ("app_state.py",     "."),
    ("ui",               "ui"),
    ("core",             "core"),
    ("assets",           "assets"),
    ("blast/windows",    "blast/windows"),
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
        "NSPrincipalClass":           "NSApplication",
        "NSHighResolutionCapable":    True,
        "CFBundleShortVersionString": "1.0.0",
    },
)
