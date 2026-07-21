"""PyInstaller runtime hook that registers Arial from Windows.

The application itself controls font sizes and export dimensions. This hook
only makes the locally installed Arial family visible to matplotlib; it does
not monkey-patch Figure.savefig or override plot-specific text sizes.

No font files are bundled or redistributed.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
from matplotlib import font_manager


def _register_windows_arial() -> None:
    if os.name != "nt":
        return

    fonts_dir = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    for filename in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
        font_path = fonts_dir / filename
        if font_path.is_file():
            try:
                font_manager.fontManager.addfont(str(font_path))
            except Exception:
                pass


_register_windows_arial()

matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.weight": "bold",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
    }
)
