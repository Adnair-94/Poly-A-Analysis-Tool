# Poly(A) Heterogeneity Analysis Tool — Desktop App v2.1

Standalone PySide6 desktop application for sequence-derived RNase T1 poly(A)
fragment annotation of BioPharma Finder LC–MS output.

The application builds an internal theoretical species library from a supplied
RNA or DNA sequence, matches observed monoisotopic masses, reports tail-length
and composition distributions, performs mode-mass QC and AWDI calculation, and
exports Excel, CSV and publication-resolution plots.

## Version 2.1 updates

### RNA and DNA sequence input

The app now accepts:

- **RNA sequence, 5′→3′**
- **DNA coding/non-template strand, 5′→3′** — converted by T→U
- **DNA template/antisense strand, 5′→3′** — reverse-complemented, then T→U

The exact RNA sequence used for the RNase T1 digest is shown in a read-only
preview and written to the Excel export.

> For DNA input, supply the **transcribed region ending at the actual
> transcription or linearisation endpoint**. Downstream vector sequence would
> alter the inferred T1 fragment and its automatic 3′ terminus.

### User-selectable termini

The generated tail-fragment library now supports independent terminal choices:

**5′ terminus**

- Auto from T1 digest
- 5′-OH
- 5′-phosphate

**3′ terminus**

- Auto from T1 digest
- 3′-OH
- 3′-phosphate
- 2′,3′-cyclic phosphate

Auto mode retains the previous rules:

- internal RNase T1 fragment: 5′-OH / 3′-phosphate
- transcript-terminal fragment: 5′-OH / 3′-OH

Manual choices override the detected tail-positive fragment. The selected and
automatic termini are both retained in exported tables.

A one-phosphate difference corresponds to **79.9663309 Da**. The annotation
output now flags observed/theoretical errors close to this value and advises the
user to review the selected terminal state.

### C/U placement

C/U variants can be generated:

- 3′ of the A-run — legacy behaviour
- 5′ of the A-run
- at both positions

Moving an identical C/U composition from one side of the A-run to the other does
**not** change intact monoisotopic mass. When both positions are generated they
are therefore isobaric and are reported as positionally ambiguous unless other
evidence resolves them. Terminal chemistry is treated independently from C/U
placement.

### TIFF or JPEG export

The plot export ZIP can contain either:

- TIFF at 600 dpi
- JPEG at 300 dpi with user-selectable quality

## Repository files

| File | Purpose |
|---|---|
| `polya_core.py` | Sequence conversion, T1 digestion, terminal mass model, library generation, matching, QC, AWDI and plot-data functions. |
| `polya_gui.py` | PySide6 desktop GUI and TIFF/JPEG/Excel/CSV exports. |
| `tests/test_polya_core.py` | Regression tests for DNA conversion, termini, C/U placement and the 79.966 Da diagnostic. |
| `requirements.txt` | Runtime and build dependencies. |
| `requirements-dev.txt` | Runtime dependencies plus pytest. |
| `build_exe.bat` | Tested workflow for building one Windows EXE. |
| `build_exe.sh` | macOS/Linux build workflow. |
| `.github/workflows/build-windows.yml` | GitHub Actions build and downloadable Windows artifact. |
| `UPDATE_NOTES.md` | Scientific and implementation rationale for this release. |

## Run from source

```bash
python -m pip install -r requirements.txt
python polya_gui.py
```

## Run tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

## Build a standalone Windows executable

On Windows, double-click:

```text
build_exe.bat
```

Output:

```text
dist\PolyA_Tool.exe
```

The target Windows 10/11 PC does not need Python, R, RStudio or packages.
The executable may trigger Windows SmartScreen until it is code-signed.

## Build with GitHub Actions

Open **Actions → Build Windows executable → Run workflow**. After the job
finishes, download the `PolyA_Tool-Windows` artifact from the workflow run.

## Input data

Sample Excel files must contain columns corresponding to:

- Monoisotopic Mass
- Apex RT
- Relative Abundance

Column-name variants are detected after case and punctuation normalization.

## Analytical cautions

- Intact mass can establish composition only to the extent supported by the
  theoretical library and mass tolerance. It cannot distinguish isobaric
  positional variants by itself.
- A manual terminus override should be documented in the exported settings.
- The mode-mass QC should be reviewed before using class-resolved or total
  distributions.
- The nearest-ladder report mode is useful for plotting but is not equivalent
  to a strict mass-confirmed composition assignment.
