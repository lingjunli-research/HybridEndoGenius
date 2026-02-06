# -*- mode: python ; coding: utf-8 -*-

import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)
a = Analysis(
    ['GUI.py'],
    pathex=[],
    binaries=[],
    datas=[('Endogenius_Thumbnail.ico', '.')],
    hiddenimports=['tarfile'],
    hookspath=['C:\\Users\\lawashburn\\Documents\\EndoGeniusDistributions\\EndoGenius_v2.0.1_final\\hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='EndoGenius',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['EndoGenius.ico'],
)
