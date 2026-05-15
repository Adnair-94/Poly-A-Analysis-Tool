#!/usr/bin/env bash
# ============================================================
# Build PolyA_Tool as a single self-contained executable.
# Run on macOS or Linux. No Python install needed on the target machine.
# ============================================================

set -e

echo
echo "=== Step 1/3: Checking Python ==="
command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 is not on PATH. Install Python 3.10 or newer."
    exit 1
}
python3 --version

echo
echo "=== Step 2/3: Installing dependencies ==="
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt

echo
echo "=== Step 3/3: Building PolyA_Tool ==="
rm -rf build dist PolyA_Tool.spec

python3 -m PyInstaller \
    --onefile \
    --windowed \
    --name PolyA_Tool \
    --collect-all PySide6 \
    --collect-all matplotlib \
    --collect-submodules openpyxl \
    --hidden-import polya_core \
    polya_gui.py

echo
echo "=== Build complete! ==="
echo
echo "Executable is at:  dist/PolyA_Tool"
echo "(On macOS this is dist/PolyA_Tool.app)"
echo
echo "Fully self-contained. Copy to any compatible machine and run."
echo
