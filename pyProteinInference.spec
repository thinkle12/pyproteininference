# -*- mode: python ; coding: utf-8 -*-
import glob, os
from pathlib import Path
from PyInstaller.utils.hooks import collect_all

import nicegui
import pyopenms

# use nicegui to find the location of its libraries
nicegui_location = f'{Path(nicegui.__file__).parent}'
pyopenms_location = f'{Path(pyopenms.__file__).parent}'
print("Found nicegui location as: " + nicegui_location)
print("Found pyopenms location as: " + pyopenms_location)
VENV_LOCATION = "./venv"
pyopenms_binaries = glob.glob(pyopenms_location + "/*.dylib") + glob.glob(pyopenms_location + "/*.dll")
print("Found pyopenms binaries to include:" + ",".join(pyopenms_binaries))
datas = [(nicegui_location, 'nicegui')]
binaries = [(dylib, '.') for dylib in pyopenms_binaries]
hiddenimports = []
tmp_ret = collect_all('pyopenms')
datas += tmp_ret[0]
binaries += tmp_ret[1]
hiddenimports += tmp_ret[2]


a = Analysis(
    ['scripts/protein_inference_gui.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='pyProteinInference',
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
)
