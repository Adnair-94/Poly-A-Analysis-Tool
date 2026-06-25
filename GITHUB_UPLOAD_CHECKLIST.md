# GitHub upload checklist

## Replace existing files

- `polya_core.py`
- `polya_gui.py`
- `README.md`
- `requirements.txt`
- `build_exe.bat`
- `build_exe.sh`
- `.gitignore`

## Add new files

- `requirements-dev.txt`
- `UPDATE_NOTES.md`
- `tests/test_polya_core.py`
- `.github/workflows/build-windows.yml`

## Suggested commit message

```text
Add DNA/RNA input conversion, selectable termini, C/U placement and JPEG export
```

## Validation before merging

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python polya_gui.py
```

Expected test result: `9 passed`.

## Build artifact

After pushing, open the GitHub **Actions** tab and run **Build Windows
executable**. Download `PolyA_Tool-Windows` after the workflow completes.
