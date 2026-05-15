# Poly(A) Heterogeneity Analysis Tool — Desktop App

Local desktop port of the Shiny app. Same scientific behaviour:
T1 digest → internal species library → BioPharma Finder peak annotation →
class-resolved and report-style length distributions, with mode-mass QC,
AWDI, and Excel/CSV/TIFF exports.

## File overview

| File | Purpose |
|------|---------|
| `polya_core.py`    | All science (digest, mass, matching, QC, plot data). Standalone, no GUI dependency. |
| `polya_gui.py`     | PySide6 + matplotlib desktop interface. |
| `requirements.txt` | Python dependencies. |
| `build_exe.bat`    | Windows one-click build → single `.exe`. |
| `build_exe.sh`     | macOS / Linux one-click build. |

## Run from source (development)

```
pip install -r requirements.txt
python polya_gui.py
```

## Build a single executable (deployment)

### Windows

Double-click `build_exe.bat`, or in a Command Prompt:

```
build_exe.bat
```

Output: `dist\PolyA_Tool.exe` — a single file (~150–250 MB).

Copy it to any Windows 10/11 PC and double-click. **No Python or other
software needs to be installed on the target PC.**

### macOS / Linux

```
chmod +x build_exe.sh
./build_exe.sh
```

Output: `dist/PolyA_Tool` (Linux) or `dist/PolyA_Tool.app` (macOS).

## Expected behaviour on first launch of the .exe

- **5–15 second delay** before the window appears (PyInstaller `--onefile` mode
  unpacks bundled libraries into a temp directory on first run; subsequent
  launches use the cached copy and start faster).
- **Windows SmartScreen** may show "Windows protected your PC" for unsigned
  executables. Click **More info → Run anyway**. This is standard for any
  unsigned binary; the only way to remove it is to buy a code-signing
  certificate.
- The window has the same controls as the Shiny sidebar, organised into
  collapsible groups on the left, with the same four output tabs on the right
  (Summary, Observed mass distribution, Class-resolved tail distribution,
  Total report-style distribution).

## If you want faster startup

`--onefile` is convenient but slow to launch. For an instant-launch version,
edit `build_exe.bat` and change `--onefile` to `--onedir`. You'll get a
folder `dist\PolyA_Tool\` containing `PolyA_Tool.exe` plus its dependencies.
Zip the folder to share. Launch time drops from ~10 s to ~1 s.

## Workflow

1. Click **Upload sample Excel files** and select one or more BioPharma
   Finder outputs (`.xlsx`). The files need columns for Monoisotopic Mass,
   Apex RT, and Relative Abundance — flexible column names are auto-detected.
2. Edit the sequence box, sample ID, range, chemistry, and tolerances as
   needed. Defaults are the SR22 validation sequence with poly(A) range
   110–140 nt and 10 ppm matching.
3. Click **Annotate and Plot**.
4. Inspect tabs; export via the buttons at the bottom of the sidebar.

## Notes

- The `Manual/report mode: nearest-ladder classification` checkbox controls
  whether unmatched peaks are placed on the nearest sequence-derived ladder
  position for report-style plots. This reproduces the manual / RNA Forge
  plotting style. Strict mass-confirmed assignments are still flagged
  separately in the export.
- The `Mode-mass QC` block prevents class/report plots from rendering if
  the most abundant peak doesn't match the generated ladder within tolerance,
  which catches wrong U chemistry, wrong sequence context, or wrong range.
- All numeric outputs (Excel/CSV) carry the exact same column structure as
  the Shiny app's exports.
