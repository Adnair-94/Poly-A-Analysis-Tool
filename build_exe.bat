@echo off
REM ============================================================
REM Build PolyA_Tool.exe as a single self-contained executable.
REM Run this on Windows. No Python install needed on the target PC.
REM ============================================================

setlocal

echo.
echo === Step 1/3: Checking Python ===
where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python is not installed or not on PATH.
    echo Install Python 3.10 or newer from https://www.python.org/downloads/
    echo Make sure to tick "Add Python to PATH" during installation.
    pause
    exit /b 1
)
python --version

echo.
echo === Step 2/3: Installing dependencies ===
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if errorlevel 1 (
    echo ERROR: Failed to install dependencies.
    pause
    exit /b 1
)

echo.
echo === Step 3/3: Building PolyA_Tool.exe ===
REM Clean previous build artefacts so old files don't get bundled
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist PolyA_Tool.spec del PolyA_Tool.spec

python -m PyInstaller ^
    --onefile ^
    --windowed ^
    --name PolyA_Tool ^
    --collect-all PySide6 ^
    --collect-all matplotlib ^
    --collect-submodules openpyxl ^
    --hidden-import polya_core ^
    polya_gui.py

if errorlevel 1 (
    echo ERROR: PyInstaller build failed.
    pause
    exit /b 1
)

echo.
echo === Build complete! ===
echo.
echo Single executable is at:  dist\PolyA_Tool.exe
echo.
echo This .exe is fully self-contained. Copy it to any Windows PC
echo (Windows 10/11) and double-click to run. No Python needed.
echo.
pause
