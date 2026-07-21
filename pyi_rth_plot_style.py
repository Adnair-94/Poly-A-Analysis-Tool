"""PyInstaller runtime hook for consistent publication-quality plot typography.

This hook runs before the desktop application imports matplotlib. It registers
Arial from the local Windows font directory and applies explicit font sizes to
all figure text immediately before export.

No font files are bundled or redistributed.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
from matplotlib import font_manager
from matplotlib.figure import Figure

AXIS_FONT_SIZE = 12
LEGEND_FONT_SIZE = 10
ANNOTATION_FONT_SIZE = 10


def _register_windows_arial() -> None:
    """Register Arial from the user's Windows installation, when available."""
    if os.name != "nt":
        return

    fonts_dir = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    for filename in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
        font_path = fonts_dir / filename
        if font_path.is_file():
            try:
                font_manager.fontManager.addfont(str(font_path))
            except Exception:
                # Matplotlib will use its configured fallback if registration fails.
                pass


def _set_text_style(text, size: int) -> None:
    try:
        text.set_fontfamily("Arial")
        text.set_fontweight("bold")
        text.set_fontsize(size)
        text.set_color("black")
    except Exception:
        pass


def _apply_export_typography(figure: Figure) -> None:
    """Apply explicit typography immediately before saving a figure."""
    for axis in figure.get_axes():
        _set_text_style(axis.xaxis.label, AXIS_FONT_SIZE)
        _set_text_style(axis.yaxis.label, AXIS_FONT_SIZE)
        _set_text_style(axis.title, AXIS_FONT_SIZE)

        for label in list(axis.get_xticklabels()) + list(axis.get_yticklabels()):
            _set_text_style(label, AXIS_FONT_SIZE)

        legend = axis.get_legend()
        if legend is not None:
            for label in legend.get_texts():
                _set_text_style(label, LEGEND_FONT_SIZE)
            _set_text_style(legend.get_title(), LEGEND_FONT_SIZE)

        # Statistical summary boxes and any other axis-level text.
        for text in axis.texts:
            _set_text_style(text, ANNOTATION_FONT_SIZE)

    for legend in getattr(figure, "legends", []):
        for label in legend.get_texts():
            _set_text_style(label, LEGEND_FONT_SIZE)

    for text in figure.texts:
        _set_text_style(text, ANNOTATION_FONT_SIZE)

    try:
        figure.canvas.draw()
    except Exception:
        pass


_register_windows_arial()

matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": AXIS_FONT_SIZE,
        "font.weight": "bold",
        "axes.labelsize": AXIS_FONT_SIZE,
        "axes.labelweight": "bold",
        "axes.titlesize": AXIS_FONT_SIZE,
        "axes.titleweight": "bold",
        "xtick.labelsize": AXIS_FONT_SIZE,
        "ytick.labelsize": AXIS_FONT_SIZE,
        "legend.fontsize": LEGEND_FONT_SIZE,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
    }
)

_original_savefig = Figure.savefig


def _savefig_with_publication_style(self: Figure, *args, **kwargs):
    _apply_export_typography(self)
    return _original_savefig(self, *args, **kwargs)


Figure.savefig = _savefig_with_publication_style
