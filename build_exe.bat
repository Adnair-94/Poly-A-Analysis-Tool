@echo off
REM ============================================================
REM Build PolyA_Tool.exe as a single self-contained executable.
REM Run on Windows. The target PC does not need Python installed.
REM ============================================================

setlocal

echo.
echo === Step 1/4: Checking Python ===
where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python is not installed or not on PATH.
    echo Install Python 3.10 or newer and select "Add Python to PATH".
    pause
    exit /b 1
)
python --version

echo.
echo === Step 2/4: Installing dependencies ===
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if errorlevel 1 (
    echo ERROR: Dependency installation failed.
    pause
    exit /b 1
)

echo.
echo === Step 3/4: Running core-science tests ===
python -m pip install pytest
python -m pytest -q
if errorlevel 1 (
    echo ERROR: Tests failed. The executable was not built.
    pause
    exit /b 1
)

echo.
echo === Step 4/4: Building PolyA_Tool.exe ===
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist PolyA_Tool.spec del PolyA_Tool.spec

python -m PyInstaller ^
    --noconfirm ^
    --clean ^
    --onefile ^
    --windowed ^
    --name PolyA_Tool ^
    --collect-all PySide6 ^
    --collect-all matplotlib ^
    --collect-submodules openpyxl ^
    --collect-submodules PIL ^
    --hidden-import polya_core ^
    polya_gui.py

if errorlevel 1 (
    echo ERROR: PyInstaller build failed.
    pause
    exit /b 1
)

echo.
echo === Build complete ===
echo Executable: dist\PolyA_Tool.exe
echo.
echo Copy the single EXE to a Windows 10/11 PC and run it.
echo The target PC does not need Python, R, RStudio, or packages.
echo.
pause
