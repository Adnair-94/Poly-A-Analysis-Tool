#!/usr/bin/env bash
# Build PolyA_Tool as a self-contained executable on macOS or Linux.
# Build on the operating system that will run the executable.

set -euo pipefail

echo
echo "=== Step 1/4: Checking Python ==="
command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 is not on PATH. Install Python 3.10 or newer."
    exit 1
}
python3 --version

echo
echo "=== Step 2/4: Installing dependencies ==="
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt

echo
echo "=== Step 3/4: Running core-science tests ==="
python3 -m pip install pytest
python3 -m pytest -q

echo
echo "=== Step 4/4: Building PolyA_Tool ==="
rm -rf build dist PolyA_Tool.spec

python3 -m PyInstaller \
    --noconfirm \
    --clean \
    --onefile \
    --windowed \
    --name PolyA_Tool \
    --collect-all PySide6 \
    --collect-all matplotlib \
    --collect-submodules openpyxl \
    --collect-submodules PIL \
    --hidden-import polya_core \
    polya_gui.py

echo
echo "=== Build complete ==="
echo "Linux output: dist/PolyA_Tool"
echo "macOS output: dist/PolyA_Tool.app or dist/PolyA_Tool"
