# -*- mode: python ; coding: utf-8 -*-
#
# PyInstaller recipe of the simulation_explorer executables, which don't need Python installed.
# The executables are not kept in the repository: .github/workflows/build-simulation-explorer.yml
# builds them for Windows, Linux and macOS and attaches them to each LoKI-GM release.
#
# To build them locally (in a clean virtual environment, so that only what is needed is included):
#     pip install -r requirements-build.txt
#     pyinstaller --noconfirm --clean simulation_explorer.spec
# The result is in dist/: a single file on Windows and Linux (simulation_explorer.exe /
# simulation_explorer), and an app bundle on macOS (simulation_explorer.app).

import os
import subprocess
import sys

NAME = 'simulation_explorer'


def linux_dependencies(paths):
    # Full paths of the libraries needed by the given ones, directly or not (ldd resolves the whole tree)
    found = set()
    for path in paths:
        output = subprocess.run(['ldd', path], capture_output=True, text=True).stdout
        for line in output.splitlines():
            target = line.partition('=>')[2].split()
            if target and target[0].startswith('/'):
                found.add(os.path.realpath(target[0]))
    return found


def linux_binaries(binaries):
    # Qt's GTK3 theme plugin would bring the GTK of the build system (GTK, pango, cairo...), which
    # can clash with the GTK modules and themes of the user's system; without it Qt uses its own
    # style. So the plugin is left out, together with the libraries that only it needs.
    # libstdc++ and libgcc_s are left out too: every system has them, and a copy older than the
    # system's breaks the system libraries loaded later (e.g. graphics drivers: "GLIBCXX_... not found")
    gtk3 = [src for dest, src, kind in binaries if dest.endswith('platformthemes/libqgtk3.so')]
    removed = {os.path.realpath(src) for src in gtk3} | linux_dependencies(gtk3)
    removed -= linux_dependencies(src for dest, src, kind in binaries if os.path.realpath(src) not in removed)
    return [(dest, src, kind) for dest, src, kind in binaries
            if os.path.realpath(src) not in removed
            and not os.path.basename(dest).startswith(('libstdc++.so', 'libgcc_s.so'))]

a = Analysis(
    ['simulation_explorer.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    # matplotlib backends: Qt for the window, and those used to save the figure from its toolbar
    # (they are imported only when saving, so they would not be found otherwise)
    hooksconfig={'matplotlib': {'backends': ['QtAgg', 'Agg', 'PDF', 'SVG', 'PS']}},
    runtime_hooks=[],
    excludes=['tkinter'],
    noarchive=False,
    optimize=0,
)
if sys.platform.startswith('linux'):
    a.binaries = linux_binaries(a.binaries)
pyz = PYZ(a.pure)

if sys.platform == 'darwin':
    # macOS: an app bundle, which needs the one-folder mode
    exe = EXE(
        pyz,
        a.scripts,
        [],
        exclude_binaries=True,
        name=NAME,
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,
        console=False,
        disable_windowed_traceback=False,
        argv_emulation=False,
        target_arch=None,
        codesign_identity=None,
        entitlements_file=None,
    )
    coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name=NAME)
    app = BUNDLE(
        coll,
        name=f'{NAME}.app',
        icon=None,
        bundle_identifier=None,
        info_plist={'NSHighResolutionCapable': True},
    )
else:
    # Windows and Linux: a single file. Without console (console=False), so that on Windows no
    # black window opens behind the window; the terminal mode (--cli) then needs Python there
    exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.datas,
        [],
        name=NAME,
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,  # UPX compression gains little with Qt and makes antivirus false positives more likely
        upx_exclude=[],
        runtime_tmpdir=None,
        console=False,
        disable_windowed_traceback=False,
        argv_emulation=False,
        target_arch=None,
        codesign_identity=None,
        entitlements_file=None,
    )
