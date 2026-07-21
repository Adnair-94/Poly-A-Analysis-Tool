Poly(A) plot export correction

Replace these two files in the repository root:
  1. polya_gui.py
  2. pyi_rth_plot_style.py

Keep the current workflow:
  .github/workflows/build-windows-exe.yml

Corrections:
  - Redraws exports on fresh figures instead of saving the live Qt canvas.
  - Single-sample plots are exactly 7 x 5 inches.
  - TIFF and JPEG/JPG are both 600 dpi.
  - Expected single-sample dimensions: 4200 x 3000 pixels.
  - Removes bbox_inches='tight' so dimensions are reproducible.
  - Axis titles: Arial Bold 12 pt.
  - Tick labels: Arial Bold 9 pt.
  - Legends and statistics: Arial Bold 9 pt.
  - JPEG uses quality control, no chroma subsampling, and optimization.
  - Runtime hook registers Arial but no longer forces every text element to 12 pt.

After replacing both files, commit them and run:
  Actions > Build Windows PolyA Tool EXE > Run workflow
