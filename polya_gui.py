"""Poly(A) Heterogeneity Analysis Tool - PySide6 desktop interface."""

from __future__ import annotations

import os
import sys
import tempfile
import traceback
import zipfile
from datetime import date
from typing import Dict, Optional

import matplotlib

# Publication plot specification. Exported single-sample figures are rendered
# independently of the Qt canvas at exactly 7 x 5 inches and 600 dpi.
FONT_FAMILY = "Arial"
AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 9
TITLE_SIZE = 10
LEGEND_SIZE = 9
ANNOTATION_SIZE = 9
EXPORT_DPI = 600
EXPORT_SINGLE_FIGSIZE = (7.0, 5.0)

matplotlib.rcParams.update(
    {
        "font.family": FONT_FAMILY,
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.weight": "bold",
        "font.size": TICK_LABEL_SIZE,
        "axes.labelsize": AXIS_LABEL_SIZE,
        "axes.labelweight": "bold",
        "axes.titlesize": TITLE_SIZE,
        "axes.titleweight": "bold",
        "xtick.labelsize": TICK_LABEL_SIZE,
        "ytick.labelsize": TICK_LABEL_SIZE,
        "legend.fontsize": LEGEND_SIZE,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
    }
)

import logging

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

import numpy as np
import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

import polya_core as pc

PLOT_FONT = {"family": FONT_FAMILY, "weight": "bold", "color": "black"}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _style_axes(ax, x_label_angle: float = 0):
    """Apply the Cas9-style publication typography to one matplotlib axis."""
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.85)

    ax.tick_params(
        axis="both",
        colors="black",
        width=0.45,
        length=4,
        labelsize=TICK_LABEL_SIZE,
    )

    for label in ax.get_xticklabels():
        label.set_fontfamily(FONT_FAMILY)
        label.set_fontweight("bold")
        label.set_fontsize(TICK_LABEL_SIZE)
        label.set_color("black")
        if x_label_angle:
            label.set_rotation(x_label_angle)
            label.set_ha("right" if x_label_angle else "center")
            label.set_va("center" if x_label_angle >= 80 else "top")

    for label in ax.get_yticklabels():
        label.set_fontfamily(FONT_FAMILY)
        label.set_fontweight("bold")
        label.set_fontsize(TICK_LABEL_SIZE)
        label.set_color("black")

    ax.xaxis.label.set_fontproperties(
        {"family": FONT_FAMILY, "weight": "bold", "size": AXIS_LABEL_SIZE}
    )
    ax.yaxis.label.set_fontproperties(
        {"family": FONT_FAMILY, "weight": "bold", "size": AXIS_LABEL_SIZE}
    )
    ax.grid(False)


def _length_breaks(values, x_limits=None, break_by=None):
    if x_limits is not None and len(x_limits) == 2 and all(np.isfinite(x_limits)):
        lo = int(np.floor(min(x_limits)))
        hi = int(np.ceil(max(x_limits)))
    else:
        array = np.asarray([value for value in values if pd.notna(value) and np.isfinite(value)])
        if len(array) == 0:
            return None
        lo = int(np.floor(array.min()))
        hi = int(np.ceil(array.max()))
    if hi < lo:
        return None
    span = hi - lo
    if break_by is not None and np.isfinite(break_by) and break_by > 0:
        step = int(break_by)
    else:
        step = 2 if span <= 30 else (5 if span <= 70 else 10)
    return list(range(lo, hi + 1, step))


def _length_limits(values):
    array = np.asarray([value for value in values if pd.notna(value) and np.isfinite(value)])
    if len(array) == 0:
        return None
    lo, hi = int(np.floor(array.min())), int(np.ceil(array.max()))
    if hi - lo <= 25:
        return lo, hi
    return int(np.floor((lo - 2) / 5) * 5), int(np.ceil((hi + 3) / 5) * 5)


def _plot_limits_from_data(data, length_col="PlotLength", custom_limits=None):
    if (
        custom_limits is not None
        and len(custom_limits) == 2
        and all(np.isfinite(custom_limits))
        and custom_limits[1] > custom_limits[0]
    ):
        return tuple(custom_limits)
    if data is None or len(data) == 0 or length_col not in data.columns:
        return None
    if "PlotMin" in data.columns and "PlotMax" in data.columns:
        lo = pd.to_numeric(data["PlotMin"], errors="coerce").dropna()
        hi = pd.to_numeric(data["PlotMax"], errors="coerce").dropna()
        if len(lo) and len(hi):
            lo_value, hi_value = float(lo.iloc[0]), float(hi.iloc[0])
            if np.isfinite(lo_value) and np.isfinite(hi_value) and hi_value > lo_value:
                return lo_value, hi_value
    return _length_limits(data[length_col])


def _set_empty(ax, message: str):
    ax.clear()
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        family=FONT_FAMILY,
        fontweight="bold",
        color="black",
        fontsize=TITLE_SIZE,
        wrap=True,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def draw_mass_plot(fig: Figure, plot_data: pd.DataFrame):
    fig.clear()
    if plot_data is None or len(plot_data) == 0:
        _set_empty(fig.add_subplot(111), "No filtered data to plot.")
        return

    class_order = [
        "A_only",
        "C_containing",
        "U_containing",
        "mixed_CU",
        "length_only",
        "unassigned",
    ]
    samples = plot_data["SampleLabel"].unique() if "SampleLabel" in plot_data else [None]
    count = len(samples)
    columns = 1 if count == 1 else 2
    rows = int(np.ceil(count / columns))
    axes = [fig.add_subplot(rows, columns, index + 1) for index in range(count)]

    handles, labels, seen = [], [], set()
    for ax, sample in zip(axes, samples):
        subset = plot_data if sample is None else plot_data[plot_data["SampleLabel"] == sample]
        present = [value for value in class_order if value in set(subset["SpeciesClass"].astype(str))]
        for species_class in present:
            data = subset[subset["SpeciesClass"] == species_class]
            ax.bar(
                data["Mass"],
                data["Fractional.Abundance"],
                width=0.1,
                color=pc.CLASS_COLOURS.get(species_class, "#999999"),
                edgecolor="black",
                linewidth=0.25,
                zorder=2,
            )
            if species_class not in seen:
                seen.add(species_class)
                handles.append(
                    Rectangle(
                        (0, 0),
                        1,
                        1,
                        facecolor=pc.CLASS_COLOURS.get(species_class, "#999999"),
                        edgecolor="black",
                        linewidth=0.5,
                    )
                )
                labels.append(pc.class_to_label(species_class))
        ax.set_xlabel("Monoisotopic mass (Da)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if sample is not None and count > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax)

    if handles:
        axes[0].legend(
            handles,
            labels,
            loc="upper right",
            frameon=True,
            edgecolor="black",
            facecolor="white",
            prop={"family": FONT_FAMILY, "weight": "bold", "size": LEGEND_SIZE},
        )
    fig.tight_layout(pad=0.9)


def draw_class_plot(fig: Figure, class_data: pd.DataFrame, x_limits=None, x_break_by=None, x_label_angle=45):
    fig.clear()
    if class_data is None or len(class_data) == 0:
        _set_empty(fig.add_subplot(111), "No class-resolved data after current filters.")
        return

    class_order = ["A_only", "C_containing", "U_containing", "mixed_CU"]
    present = [value for value in class_order if value in set(class_data["SpeciesClass"].astype(str))]
    samples = class_data["SampleLabel"].unique()
    count = len(samples)
    columns = 1 if count == 1 else 2
    rows = int(np.ceil(count / columns))
    axes = [fig.add_subplot(rows, columns, index + 1) for index in range(count)]

    handles, labels, seen = [], [], set()
    for ax, sample in zip(axes, samples):
        subset = class_data[class_data["SampleLabel"] == sample]
        if len(subset) == 0:
            _set_empty(ax, "No data for sample")
            continue
        limits = _plot_limits_from_data(subset, "PlotLength", custom_limits=x_limits) if count == 1 else x_limits
        number_of_classes = max(len(present), 1)
        bar_width = 0.74 / number_of_classes
        for index, species_class in enumerate(present):
            data = subset[subset["SpeciesClass"] == species_class]
            if len(data) == 0:
                continue
            offset = (index - (number_of_classes - 1) / 2) * bar_width
            ax.bar(
                data["PlotLength"] + offset,
                data["RelAbundance"],
                width=bar_width,
                color=pc.CLASS_COLOURS.get(species_class, "#999999"),
                edgecolor="black",
                linewidth=0.25,
                zorder=2,
            )
            if species_class not in seen:
                seen.add(species_class)
                handles.append(
                    Rectangle(
                        (0, 0),
                        1,
                        1,
                        facecolor=pc.CLASS_COLOURS.get(species_class, "#999999"),
                        edgecolor="black",
                        linewidth=0.5,
                    )
                )
                labels.append(pc.class_to_label(species_class))
        ax.set_xlabel("Length (nt)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if limits is not None:
            ax.set_xlim(limits)
        breaks = _length_breaks(subset["PlotLength"], x_limits=limits, break_by=x_break_by)
        if breaks:
            ax.set_xticks(breaks)
        if count > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax, x_label_angle=x_label_angle)

    if handles:
        axes[0].legend(
            handles,
            labels,
            loc="upper right",
            frameon=True,
            edgecolor="black",
            facecolor="white",
            prop={"family": FONT_FAMILY, "weight": "bold", "size": LEGEND_SIZE},
        )
    fig.tight_layout(pad=0.9)


def draw_total_plot(fig: Figure, total_data: pd.DataFrame, threshold_pct=10, x_limits=None, x_break_by=None, x_label_angle=45):
    fig.clear()
    if total_data is None or len(total_data) == 0:
        _set_empty(fig.add_subplot(111), "No total distribution data after current filters.")
        return

    samples = total_data["SampleLabel"].unique()
    count = len(samples)
    columns = 1 if count == 1 else 2
    rows = int(np.ceil(count / columns))
    axes = [fig.add_subplot(rows, columns, index + 1) for index in range(count)]
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=pc.CLASS_COLOURS["above_10"], edgecolor="black", linewidth=0.5),
        Rectangle((0, 0), 1, 1, facecolor=pc.CLASS_COLOURS["below_10"], edgecolor="black", linewidth=0.5),
    ]
    labels = [f"Individual bin ≥{threshold_pct}%", f"Individual bin <{threshold_pct}%"]

    for ax, sample in zip(axes, samples):
        subset = total_data[total_data["SampleLabel"] == sample].copy()
        if len(subset) == 0:
            _set_empty(ax, "No data")
            continue
        limits = _plot_limits_from_data(subset, "PlotLength", custom_limits=x_limits) if count == 1 else x_limits
        below = subset[subset["ThresholdClass"] == "below_10"]
        above = subset[subset["ThresholdClass"] == "above_10"]
        ax.bar(
            below["PlotLength"],
            below["TotalNorm"],
            width=0.85,
            color=pc.CLASS_COLOURS["below_10"],
            edgecolor="black",
            linewidth=0.25,
            zorder=2,
        )
        ax.bar(
            above["PlotLength"],
            above["TotalNorm"],
            width=0.85,
            color=pc.CLASS_COLOURS["above_10"],
            edgecolor="black",
            linewidth=0.25,
            zorder=2,
        )
        ax.set_xlabel("Length (nt)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=AXIS_LABEL_SIZE)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if limits is not None:
            ax.set_xlim(limits)
        breaks = _length_breaks(subset["PlotLength"], x_limits=limits, break_by=x_break_by)
        if breaks:
            ax.set_xticks(breaks)

        if subset["TotalRaw"].sum() > 0:
            mode_value = float(subset.loc[subset["TotalRaw"].idxmax(), "PlotLength"])
            mean_value = float((subset["PlotLength"] * subset["TotalRaw"]).sum() / subset["TotalRaw"].sum())
            left = subset[subset["PlotLength"] <= mode_value]
            right = subset[subset["PlotLength"] >= mode_value]
            sd_left = (
                float(np.sqrt((left["TotalRaw"] * (left["PlotLength"] - mode_value) ** 2).sum() / left["TotalRaw"].sum()))
                if left["TotalRaw"].sum() > 0
                else np.nan
            )
            sd_right = (
                float(np.sqrt((right["TotalRaw"] * (right["PlotLength"] - mode_value) ** 2).sum() / right["TotalRaw"].sum()))
                if right["TotalRaw"].sum() > 0
                else np.nan
            )
            stats = (
                f"Mode = {round(mode_value)} nt\nMean = {round(mean_value)} nt\n"
                f"SD-L = {round(sd_left)} nt\nSD-R = {round(sd_right)} nt"
            )
            x_position = limits[1] - 0.04 * (limits[1] - limits[0]) if limits else subset["PlotLength"].max()
            ax.text(
                x_position,
                82,
                stats,
                ha="right",
                va="top",
                family=FONT_FAMILY,
                fontweight="bold",
                fontsize=ANNOTATION_SIZE,
                bbox=dict(facecolor="white", edgecolor="black", linewidth=0.35, pad=3),
            )
        if count > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax, x_label_angle=x_label_angle)

    axes[0].legend(
        handles,
        labels,
        loc="upper right",
        frameon=True,
        edgecolor="black",
        facecolor="white",
        prop={"family": FONT_FAMILY, "weight": "bold", "size": LEGEND_SIZE},
    )
    fig.tight_layout(pad=0.9)


def _publication_figsize(sample_count: int) -> tuple[float, float]:
    """Return a fixed publication canvas; expand only for multi-sample facets."""
    sample_count = max(int(sample_count), 1)
    if sample_count == 1:
        return EXPORT_SINGLE_FIGSIZE
    rows = int(np.ceil(sample_count / 2))
    return 10.0, 4.8 * rows


# ---------------------------------------------------------------------------
# Table helper
# ---------------------------------------------------------------------------


def populate_table(widget: QTableWidget, dataframe: pd.DataFrame, max_rows: Optional[int] = None):
    if dataframe is None or len(dataframe) == 0:
        widget.clear()
        widget.setRowCount(1)
        widget.setColumnCount(1)
        widget.setHorizontalHeaderLabels(["Result"])
        widget.setItem(0, 0, QTableWidgetItem("No data to display"))
        return
    data = dataframe.head(max_rows) if max_rows is not None else dataframe
    widget.clear()
    widget.setRowCount(len(data))
    widget.setColumnCount(len(data.columns))
    widget.setHorizontalHeaderLabels([str(column) for column in data.columns])
    for row_index in range(len(data)):
        for column_index, column in enumerate(data.columns):
            value = data.iloc[row_index, column_index]
            if pd.isna(value):
                text = ""
            elif isinstance(value, (float, np.floating)):
                text = f"{value:.6g}"
            else:
                text = str(value)
            item = QTableWidgetItem(text)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            widget.setItem(row_index, column_index, item)
    widget.resizeColumnsToContents()
    widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"Poly(A) Heterogeneity Analysis Tool v{pc.APP_VERSION}")
        self.resize(1550, 980)
        self.uploaded_files: list[str] = []
        self.generated: Optional[pc.GeneratedLibrary] = None
        self.results: Dict[str, dict] = {}
        self.qc_table: Optional[pd.DataFrame] = None
        self._build_ui()
        self._wire_signals()
        self._apply_demo("SR22")

    # ------------------------------------------------------------------ UI

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer = QHBoxLayout(central)
        outer.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Horizontal)
        outer.addWidget(splitter)

        controls_holder = QScrollArea()
        controls_holder.setWidgetResizable(True)
        controls_holder.setMinimumWidth(480)
        controls_inner = QWidget()
        controls_holder.setWidget(controls_inner)
        controls = QVBoxLayout(controls_inner)
        controls.setContentsMargins(10, 10, 10, 10)
        controls.setSpacing(8)

        controls.addWidget(self._build_files_group())
        controls.addWidget(self._build_sequence_group())
        controls.addWidget(self._build_terminal_group())
        controls.addWidget(self._build_range_group())
        controls.addWidget(self._build_cu_group())
        controls.addWidget(self._build_filter_group())
        controls.addWidget(self._build_match_group())
        controls.addWidget(self._build_report_length_group())
        controls.addWidget(self._build_axis_export_group())
        controls.addWidget(self._build_awdi_group())

        self.warnings_label = QLabel()
        self.warnings_label.setWordWrap(True)
        self.warnings_label.setStyleSheet(
            "background:#fff3cd; border:1px solid #ffecb5; color:#664d03; padding:8px; font-weight:bold;"
        )
        self.warnings_label.setVisible(False)
        controls.addWidget(self.warnings_label)

        self.run_btn = QPushButton("Annotate and Plot")
        self.run_btn.setStyleSheet(
            "QPushButton { background:#2563eb; color:white; font-weight:bold; padding:10px; }"
            "QPushButton:hover { background:#1d4ed8; }"
        )
        controls.addWidget(self.run_btn)
        controls.addWidget(self._build_export_group())
        controls.addStretch(1)
        splitter.addWidget(controls_holder)

        self.tabs = QTabWidget()
        splitter.addWidget(self.tabs)
        splitter.setSizes([490, 1060])
        self._build_summary_tab()
        self._build_mass_tab()
        self._build_class_tab()
        self._build_total_tab()

    def _build_files_group(self):
        group = QGroupBox("Sample files")
        layout = QVBoxLayout(group)
        self.files_btn = QPushButton("Upload sample Excel files (.xlsx)")
        layout.addWidget(self.files_btn)
        self.files_label = QLabel("No files selected")
        self.files_label.setWordWrap(True)
        layout.addWidget(self.files_label)
        return group

    def _build_sequence_group(self):
        group = QGroupBox("Sequence input and conversion")
        form = QFormLayout(group)
        self.sample_id_in = QLineEdit("SR22")
        form.addRow("Sample ID:", self.sample_id_in)

        self.sequence_in = QTextEdit()
        self.sequence_in.setMinimumHeight(105)
        form.addRow("Input sequence (5'→3'):", self.sequence_in)

        self.input_molecule_in = QComboBox()
        self.input_molecule_in.addItem("RNA sequence (5'→3')", "rna")
        self.input_molecule_in.addItem("DNA coding/non-template strand (5'→3')", "dna_coding")
        self.input_molecule_in.addItem("DNA template/antisense strand (5'→3')", "dna_template")
        form.addRow("Input molecule/strand:", self.input_molecule_in)

        self.sequence_mode_in = QComboBox()
        self.sequence_mode_in.addItem("Digest transcribed sequence using RNase T1 rules", "digest_after_G")
        self.sequence_mode_in.addItem("Input is already the T1 tail fragment", "tail_fragment")
        form.addRow("Analysis input mode:", self.sequence_mode_in)

        buttons = QHBoxLayout()
        self.dna_to_rna_btn = QPushButton("DNA → RNA preview")
        self.rna_to_dna_btn = QPushButton("RNA → DNA preview")
        buttons.addWidget(self.dna_to_rna_btn)
        buttons.addWidget(self.rna_to_dna_btn)
        wrapper = QWidget()
        wrapper.setLayout(buttons)
        form.addRow("Conversion helpers:", wrapper)

        self.converted_rna_preview = QTextEdit()
        self.converted_rna_preview.setReadOnly(True)
        self.converted_rna_preview.setMinimumHeight(70)
        form.addRow("RNA used for T1 digest:", self.converted_rna_preview)

        self.demo_in = QComboBox()
        self.demo_in.addItem("Keep current", "none")
        self.demo_in.addItem("SR22", "SR22")
        self.demo_in.addItem("CSP split tail", "CSP")
        form.addRow("Load example sequence:", self.demo_in)

        self.chem_in = QComboBox()
        self.chem_in.addItem("Canonical U", "canonical_U")
        self.chem_in.addItem("N1-methylpseudouridine at U positions", "N1mPseudoU")
        form.addRow("Nucleotide chemistry:", self.chem_in)

        self.u_shift_in = QDoubleSpinBox()
        self.u_shift_in.setRange(-100, 100)
        self.u_shift_in.setDecimals(4)
        self.u_shift_in.setSingleStep(0.0001)
        form.addRow("Advanced U-residue offset (Da):", self.u_shift_in)
        return group

    def _build_terminal_group(self):
        group = QGroupBox("Tail-fragment termini")
        form = QFormLayout(group)
        self.five_prime_in = QComboBox()
        self.five_prime_in.addItem("Auto from T1 digest", "auto")
        self.five_prime_in.addItem("5'-OH", "5OH")
        self.five_prime_in.addItem("5'-phosphate", "5p")
        form.addRow("5' terminus:", self.five_prime_in)

        self.three_prime_in = QComboBox()
        self.three_prime_in.addItem("Auto from T1 digest", "auto")
        self.three_prime_in.addItem("3'-OH", "3OH")
        self.three_prime_in.addItem("3'-phosphate", "3p")
        self.three_prime_in.addItem("2',3'-cyclic phosphate", "3cyclicp")
        form.addRow("3' terminus:", self.three_prime_in)

        note = QLabel(
            "Auto retains 5'-OH/3'-phosphate for an internal T1 fragment and 5'-OH/3'-OH for the transcript-terminal fragment. Manual choices override the detected tail fragment only."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color:#444;")
        form.addRow(note)
        return group

    def _build_range_group(self):
        group = QGroupBox("Range and tail detection")
        form = QFormLayout(group)
        self.auto_range_in = QCheckBox("Auto range around detected A-run(s)")
        form.addRow(self.auto_range_in)
        self.expected_a_in = QSpinBox()
        self.expected_a_in.setRange(0, 10000)
        self.expected_a_in.setValue(122)
        form.addRow("Expected poly(A) length:", self.expected_a_in)
        self.range_min_in = QSpinBox()
        self.range_min_in.setRange(0, 10000)
        self.range_min_in.setValue(110)
        self.range_max_in = QSpinBox()
        self.range_max_in.setRange(0, 10000)
        self.range_max_in.setValue(140)
        form.addRow("Range minimum:", self.range_min_in)
        form.addRow("Range maximum:", self.range_max_in)
        self.range_basis_in = QComboBox()
        self.range_basis_in.addItem("Report length used in plots", "report_length")
        self.range_basis_in.addItem("A-count only", "a_count")
        form.addRow("Min/max refers to:", self.range_basis_in)
        self.min_a_run_in = QSpinBox()
        self.min_a_run_in.setRange(1, 1000)
        self.min_a_run_in.setValue(10)
        form.addRow("Minimum A-run for detection:", self.min_a_run_in)
        self.auto_left_in = QSpinBox()
        self.auto_left_in.setRange(0, 1000)
        self.auto_left_in.setValue(10)
        self.auto_right_in = QSpinBox()
        self.auto_right_in.setRange(0, 1000)
        self.auto_right_in.setValue(15)
        form.addRow("Auto range: nt below detected:", self.auto_left_in)
        form.addRow("Auto range: nt above detected:", self.auto_right_in)
        return group

    def _build_cu_group(self):
        group = QGroupBox("C/U variants in generated library")
        form = QFormLayout(group)
        self.max_c_in = QSpinBox()
        self.max_c_in.setRange(0, 10)
        self.max_c_in.setValue(1)
        self.max_u_in = QSpinBox()
        self.max_u_in.setRange(0, 10)
        self.max_u_in.setValue(1)
        form.addRow("Maximum C residues:", self.max_c_in)
        form.addRow("Maximum U residues:", self.max_u_in)

        self.consider_c_in = QCheckBox("Consider C-containing species")
        self.consider_c_in.setChecked(True)
        self.consider_u_in = QCheckBox("Consider U-containing species")
        self.consider_u_in.setChecked(False)
        self.allow_mixed_in = QCheckBox("Allow mixed C/U species")
        form.addRow(self.consider_c_in)
        form.addRow(self.consider_u_in)
        form.addRow(self.allow_mixed_in)

        self.variant_placement_in = QComboBox()
        self.variant_placement_in.addItem("3' of the A-run (legacy)", "three_prime_of_a")
        self.variant_placement_in.addItem("5' of the A-run", "five_prime_of_a")
        self.variant_placement_in.addItem("Generate both positional possibilities", "both")
        form.addRow("C/U placement:", self.variant_placement_in)

        row = QHBoxLayout()
        self.show_a_in = QCheckBox("Poly(A)")
        self.show_a_in.setChecked(True)
        self.show_c_in = QCheckBox("with C")
        self.show_c_in.setChecked(True)
        self.show_u_in = QCheckBox("with U")
        self.show_u_in.setChecked(True)
        self.show_cu_in = QCheckBox("with C/U")
        self.show_cu_in.setChecked(True)
        for checkbox in (self.show_a_in, self.show_c_in, self.show_u_in, self.show_cu_in):
            row.addWidget(checkbox)
        wrapper = QWidget()
        wrapper.setLayout(row)
        form.addRow("Classes in plots:", wrapper)
        return group

    def _build_filter_group(self):
        group = QGroupBox("Filter input peaks")
        form = QFormLayout(group)
        self.rt_min_in = QDoubleSpinBox()
        self.rt_min_in.setRange(0, 1000)
        self.rt_min_in.setDecimals(2)
        self.rt_max_in = QDoubleSpinBox()
        self.rt_max_in.setRange(0, 1000)
        self.rt_max_in.setDecimals(2)
        self.rt_max_in.setValue(100)
        form.addRow("RT minimum (min):", self.rt_min_in)
        form.addRow("RT maximum (min):", self.rt_max_in)
        self.rel_min_in = QDoubleSpinBox()
        self.rel_min_in.setRange(0, 100)
        self.rel_max_in = QDoubleSpinBox()
        self.rel_max_in.setRange(0, 100)
        self.rel_max_in.setValue(100)
        form.addRow("Relative abundance min (%):", self.rel_min_in)
        form.addRow("Relative abundance max (%):", self.rel_max_in)
        self.mass_min_in = QDoubleSpinBox()
        self.mass_min_in.setRange(0, 1e7)
        self.mass_min_in.setDecimals(1)
        self.mass_min_in.setValue(5000)
        self.mass_max_in = QDoubleSpinBox()
        self.mass_max_in.setRange(0, 1e7)
        self.mass_max_in.setDecimals(1)
        self.mass_max_in.setValue(75000)
        form.addRow("Mass minimum (Da):", self.mass_min_in)
        form.addRow("Mass maximum (Da):", self.mass_max_in)
        return group

    def _build_match_group(self):
        group = QGroupBox("Matching and QC")
        form = QFormLayout(group)
        self.tol_unit_in = QComboBox()
        self.tol_unit_in.addItems(["ppm", "Da"])
        form.addRow("Direct mass-match unit:", self.tol_unit_in)
        self.tol_value_in = QDoubleSpinBox()
        self.tol_value_in.setRange(0.001, 10000)
        self.tol_value_in.setDecimals(3)
        self.tol_value_in.setValue(10)
        form.addRow("Direct match tolerance:", self.tol_value_in)
        self.fill_in = QCheckBox("Nearest-library fill outside direct tolerance")
        self.fill_in.setChecked(True)
        form.addRow(self.fill_in)
        self.fill_tol_in = QDoubleSpinBox()
        self.fill_tol_in.setRange(0, 1000)
        self.fill_tol_in.setDecimals(2)
        self.fill_tol_in.setValue(1)
        form.addRow("Nearest-fill tolerance (Da):", self.fill_tol_in)
        self.length_tol_in = QDoubleSpinBox()
        self.length_tol_in.setRange(0, 10000)
        self.length_tol_in.setDecimals(1)
        self.length_tol_in.setValue(165)
        form.addRow("Length-only tolerance (Da):", self.length_tol_in)
        self.legacy_in = QCheckBox("Manual/report mode: nearest-ladder classification")
        self.legacy_in.setChecked(True)
        form.addRow(self.legacy_in)
        self.qc_enable_in = QCheckBox("Mode-mass QC blocks class/report plots")
        self.qc_enable_in.setChecked(True)
        form.addRow(self.qc_enable_in)
        self.qc_unit_in = QComboBox()
        self.qc_unit_in.addItems(["Da", "ppm"])
        form.addRow("Mode QC unit:", self.qc_unit_in)
        self.qc_tol_in = QDoubleSpinBox()
        self.qc_tol_in.setRange(0.001, 1000)
        self.qc_tol_in.setDecimals(3)
        self.qc_tol_in.setValue(1)
        form.addRow("Mode QC tolerance:", self.qc_tol_in)
        self.dup_rule_in = QComboBox()
        self.dup_rule_in.addItem("Use dominant duplicate entry", "dominant")
        self.dup_rule_in.addItem("Sum duplicate entries", "sum")
        form.addRow("Duplicate entries:", self.dup_rule_in)
        return group

    def _build_report_length_group(self):
        group = QGroupBox("Report length convention")
        form = QFormLayout(group)
        self.report_basis_in = QComboBox()
        self.report_basis_in.addItem("A-count only", "a_count")
        self.report_basis_in.addItem("A-count + fixed T1 context", "context")
        self.report_basis_in.addItem("Full generated fragment length", "full_fragment")
        self.report_basis_in.addItem("A-count + manual shift", "manual_shift")
        self.report_basis_in.setCurrentIndex(1)
        form.addRow("Length shown in plots:", self.report_basis_in)
        self.manual_shift_in = QSpinBox()
        self.manual_shift_in.setRange(-500, 500)
        form.addRow("Manual shift (nt):", self.manual_shift_in)
        return group

    def _build_axis_export_group(self):
        group = QGroupBox("Plot display and image export")
        form = QFormLayout(group)
        self.use_custom_x_in = QCheckBox("Manually set x-axis range")
        form.addRow(self.use_custom_x_in)
        self.xmin_in = QSpinBox()
        self.xmin_in.setRange(0, 100000)
        self.xmin_in.setValue(50)
        self.xmax_in = QSpinBox()
        self.xmax_in.setRange(1, 100000)
        self.xmax_in.setValue(100)
        form.addRow("X-axis minimum (nt):", self.xmin_in)
        form.addRow("X-axis maximum (nt):", self.xmax_in)
        self.xbreak_in = QSpinBox()
        self.xbreak_in.setRange(1, 100)
        self.xbreak_in.setValue(1)
        form.addRow("X-axis tick interval (nt):", self.xbreak_in)
        self.xangle_in = QSpinBox()
        self.xangle_in.setRange(0, 90)
        self.xangle_in.setValue(45)
        form.addRow("X-axis label angle (degrees):", self.xangle_in)
        self.threshold_in = QSpinBox()
        self.threshold_in.setRange(0, 100)
        self.threshold_in.setValue(10)
        form.addRow("Individual-bin display threshold (%):", self.threshold_in)

        self.plot_format_in = QComboBox()
        self.plot_format_in.addItem("TIFF (600 dpi, 7 × 5 in)", "tiff")
        self.plot_format_in.addItem("JPEG/JPG (600 dpi, 7 × 5 in)", "jpeg")
        form.addRow("Export image format:", self.plot_format_in)
        self.jpeg_quality_in = QSpinBox()
        self.jpeg_quality_in.setRange(50, 100)
        self.jpeg_quality_in.setValue(95)
        form.addRow("JPEG quality:", self.jpeg_quality_in)
        return group

    def _build_awdi_group(self):
        group = QGroupBox("AWDI")
        form = QFormLayout(group)
        self.awdi_basis_in = QComboBox()
        self.awdi_basis_in.addItem("A-only assigned species", "A_only")
        self.awdi_basis_in.addItem("All assigned species", "all_assigned")
        self.awdi_basis_in.addItem("A-only + C-containing species", "report_default")
        form.addRow("AWDI basis:", self.awdi_basis_in)
        return group

    def _build_export_group(self):
        group = QGroupBox("Export")
        layout = QVBoxLayout(group)
        self.exp_excel_btn = QPushButton("Download summary + annotated data (Excel)")
        self.exp_csv_btn = QPushButton("Download annotated data (CSV)")
        self.exp_lib_btn = QPushButton("Download generated species library (CSV)")
        self.exp_digest_btn = QPushButton("Download T1 digest table (CSV)")
        self.exp_plots_btn = QPushButton("Download plots (TIFF + ZIP)")
        for button in (
            self.exp_excel_btn,
            self.exp_csv_btn,
            self.exp_lib_btn,
            self.exp_digest_btn,
            self.exp_plots_btn,
        ):
            layout.addWidget(button)
        return group

    def _build_summary_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.addWidget(QLabel("<b>Summary</b>"))
        self.summary_table = QTableWidget()
        layout.addWidget(self.summary_table, 2)
        layout.addWidget(QLabel("<b>Mode-mass QC</b>"))
        self.mode_qc_table_w = QTableWidget()
        layout.addWidget(self.mode_qc_table_w, 2)
        self.qc_banner = QLabel()
        self.qc_banner.setWordWrap(True)
        self.qc_banner.setVisible(False)
        layout.addWidget(self.qc_banner)
        layout.addWidget(QLabel("<b>Annotation preview</b>"))
        self.annot_table = QTableWidget()
        layout.addWidget(self.annot_table, 3)
        layout.addWidget(QLabel("<b>Generated species preview</b>"))
        self.lib_table = QTableWidget()
        layout.addWidget(self.lib_table, 2)
        layout.addWidget(QLabel("<b>Detected T1 tail blocks</b>"))
        self.tail_block_table_w = QTableWidget()
        layout.addWidget(self.tail_block_table_w, 1)
        self.tabs.addTab(widget, "Summary")

    def _build_mass_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        self.mass_fig = Figure(figsize=(8, 5), dpi=100)
        self.mass_canvas = FigureCanvas(self.mass_fig)
        layout.addWidget(self.mass_canvas)
        self.tabs.addTab(widget, "Observed mass distribution")

    def _build_class_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        self.class_fig = Figure(figsize=(8, 5), dpi=100)
        self.class_canvas = FigureCanvas(self.class_fig)
        layout.addWidget(self.class_canvas)
        self.tabs.addTab(widget, "Class-resolved tail distribution")

    def _build_total_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        self.total_fig = Figure(figsize=(8, 5), dpi=100)
        self.total_canvas = FigureCanvas(self.total_fig)
        layout.addWidget(self.total_canvas)
        self.tabs.addTab(widget, "Total report-style distribution")

    def _wire_signals(self):
        self.files_btn.clicked.connect(self.on_pick_files)
        self.run_btn.clicked.connect(self.on_run)
        self.demo_in.currentIndexChanged.connect(self._on_demo_changed)
        self.dna_to_rna_btn.clicked.connect(self._refresh_sequence_preview)
        self.rna_to_dna_btn.clicked.connect(self._show_rna_to_dna_preview)
        self.sequence_in.textChanged.connect(self._refresh_sequence_preview)
        self.input_molecule_in.currentIndexChanged.connect(self._refresh_sequence_preview)
        self.plot_format_in.currentIndexChanged.connect(self._update_export_button_text)

        self.exp_excel_btn.clicked.connect(self.on_export_excel)
        self.exp_csv_btn.clicked.connect(self.on_export_annotated_csv)
        self.exp_lib_btn.clicked.connect(self.on_export_lib_csv)
        self.exp_digest_btn.clicked.connect(self.on_export_digest_csv)
        self.exp_plots_btn.clicked.connect(self.on_export_plots)

        for widget in self.findChildren(QCheckBox):
            widget.toggled.connect(self._refresh_warnings)
        for widget in self.findChildren(QComboBox):
            widget.currentIndexChanged.connect(self._refresh_warnings)
        for widget in self.findChildren(QSpinBox):
            widget.valueChanged.connect(self._refresh_warnings)
        for widget in self.findChildren(QDoubleSpinBox):
            widget.valueChanged.connect(self._refresh_warnings)
        for widget in self.findChildren(QLineEdit):
            widget.textChanged.connect(self._refresh_warnings)
        self._update_export_button_text()

    # ------------------------------------------------------------ interaction

    def _on_demo_changed(self, *_):
        demo = self.demo_in.currentData()
        if demo and demo != "none":
            self._apply_demo(demo)
            self.demo_in.blockSignals(True)
            self.demo_in.setCurrentIndex(0)
            self.demo_in.blockSignals(False)

    def _apply_demo(self, name: str):
        if name == "CSP":
            self.sample_id_in.setText("CSP")
            self.sequence_in.setPlainText(pc.CSP_SEQ)
            self.expected_a_in.setValue(70)
            self.range_min_in.setValue(60)
            self.range_max_in.setValue(85)
        else:
            self.sample_id_in.setText("SR22")
            self.sequence_in.setPlainText(pc.SR22_SEQ)
            self.expected_a_in.setValue(122)
            self.range_min_in.setValue(110)
            self.range_max_in.setValue(140)
        self.input_molecule_in.setCurrentIndex(0)
        self.sequence_mode_in.setCurrentIndex(0)
        self._refresh_sequence_preview()
        self._refresh_warnings()

    def _refresh_sequence_preview(self, *_):
        try:
            rna = pc.convert_sequence_to_rna(
                self.sequence_in.toPlainText(), self.input_molecule_in.currentData()
            )
            self.converted_rna_preview.setPlainText(rna)
            self.converted_rna_preview.setStyleSheet("color:#111; background:#f7fff7;")
        except Exception as exc:
            self.converted_rna_preview.setPlainText(str(exc))
            self.converted_rna_preview.setStyleSheet("color:#8b0000; background:#fff0f0;")
        self._refresh_warnings()

    def _show_rna_to_dna_preview(self):
        try:
            rna = pc.convert_sequence_to_rna(
                self.sequence_in.toPlainText(), self.input_molecule_in.currentData()
            )
            QMessageBox.information(self, "RNA → DNA", pc.rna_to_dna(rna))
        except Exception as exc:
            QMessageBox.warning(self, "Conversion failed", str(exc))

    def _update_export_button_text(self, *_):
        format_name = "TIFF" if self.plot_format_in.currentData() == "tiff" else "JPEG"
        self.exp_plots_btn.setText(f"Download plots ({format_name} + ZIP)")
        self.jpeg_quality_in.setEnabled(self.plot_format_in.currentData() == "jpeg")

    def _refresh_warnings(self, *_):
        if not hasattr(self, "warnings_label"):
            return
        warnings = []
        input_type = self.input_molecule_in.currentData()
        if input_type == "dna_template":
            warnings.append(
                "DNA template/antisense input will be reverse-complemented before T→U conversion. Confirm that the supplied strand is written 5'→3'."
            )
        if input_type in ("dna_coding", "dna_template") and self.sequence_mode_in.currentData() == "digest_after_G":
            warnings.append(
                "For DNA input, supply the transcribed region ending at the actual transcription/linearisation endpoint. Downstream vector sequence would change the T1 fragment and automatic 3' terminus."
            )
        if self.five_prime_in.currentData() != "auto" or self.three_prime_in.currentData() != "auto":
            warnings.append(
                "A manual terminal override is active for tail-positive fragments. The selected terminal state changes theoretical mass independently of C/U position."
            )
        if self.variant_placement_in.currentData() == "both":
            warnings.append(
                "Both C/U positions are being generated. They are isobaric when composition and termini are identical, so intact mass alone will report positional ambiguity."
            )
        if self.range_max_in.value() < self.range_min_in.value():
            warnings.append("Range maximum is below range minimum.")
        if self.mass_max_in.value() < self.mass_min_in.value():
            warnings.append("Mass maximum is below mass minimum.")
        if self.rt_max_in.value() < self.rt_min_in.value():
            warnings.append("RT maximum is below RT minimum.")
        self.warnings_label.setText("<br>".join(f"• {message}" for message in warnings))
        self.warnings_label.setVisible(bool(warnings))

    def on_pick_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select BioPharma Finder Excel files",
            "",
            "Excel files (*.xlsx *.xls)",
        )
        if files:
            self.uploaded_files = files
            self.files_label.setText("\n".join(os.path.basename(path) for path in files))

    def _effective_max_c(self):
        return self.max_c_in.value() if self.consider_c_in.isChecked() else 0

    def _effective_max_u(self):
        return self.max_u_in.value() if self.consider_u_in.isChecked() else 0

    def _axis_settings(self):
        limits = None
        if self.use_custom_x_in.isChecked():
            limits = (self.xmin_in.value(), self.xmax_in.value())
        return {
            "limits": limits,
            "break_by": self.xbreak_in.value(),
            "label_angle": self.xangle_in.value(),
        }

    def on_run(self):
        if not self.uploaded_files:
            QMessageBox.warning(self, "No sample files", "Upload at least one sample Excel file.")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self.generated = pc.generate_species_library(
                sample_id=self.sample_id_in.text(),
                seq=self.sequence_in.toPlainText(),
                chemistry=self.chem_in.currentData(),
                u_mass_shift=self.u_shift_in.value(),
                sequence_mode=self.sequence_mode_in.currentData(),
                sequence_input_type=self.input_molecule_in.currentData(),
                use_auto_ranges=self.auto_range_in.isChecked(),
                expected_a_length=self.expected_a_in.value(),
                range_min=self.range_min_in.value(),
                range_max=self.range_max_in.value(),
                range_input_basis=self.range_basis_in.currentData(),
                report_length_basis=self.report_basis_in.currentData(),
                manual_report_shift=self.manual_shift_in.value(),
                min_a_run=self.min_a_run_in.value(),
                auto_range_left=self.auto_left_in.value(),
                auto_range_right=self.auto_right_in.value(),
                max_c=self._effective_max_c(),
                max_u=self._effective_max_u(),
                allow_mixed=self.allow_mixed_in.isChecked(),
                variant_placement=self.variant_placement_in.currentData(),
                five_prime_terminus=self.five_prime_in.currentData(),
                three_prime_terminus=self.three_prime_in.currentData(),
            )

            results = {}
            for path in self.uploaded_files:
                name = os.path.basename(path)
                raw = pc.preprocess_uploaded_excel(path)
                filtered = pc.filter_peak_data(
                    raw,
                    rt_range=(self.rt_min_in.value(), self.rt_max_in.value()),
                    abundance_range=(self.rel_min_in.value(), self.rel_max_in.value()),
                    mass_range=(self.mass_min_in.value(), self.mass_max_in.value()),
                )
                if len(filtered):
                    annotated = pc.match_to_library(
                        filtered,
                        self.generated.SpeciesLibrary,
                        tolerance_value=self.tol_value_in.value(),
                        tolerance_unit=self.tol_unit_in.currentText(),
                        fill_unassigned=self.fill_in.isChecked(),
                        fill_tolerance_da=self.fill_tol_in.value(),
                        length_only_tolerance_da=self.length_tol_in.value(),
                        legacy_nearest_ladder=self.legacy_in.isChecked(),
                    )
                    annotated = pc.add_report_length_fields(
                        annotated,
                        report_length_basis=self.report_basis_in.currentData(),
                        manual_report_shift=self.manual_shift_in.value(),
                    )
                    annotated["File"] = name
                    annotated["SampleLabel"] = pc.clean_sample_label(name)
                else:
                    annotated = pd.DataFrame()
                awdi = pc.calculate_awdi_from_annotated(
                    annotated, basis=self.awdi_basis_in.currentData()
                )
                results[name] = {
                    "raw": raw,
                    "filtered": filtered,
                    "annotated": annotated,
                    "awdi": awdi,
                }

            self.results = results
            self.qc_table = pc.build_mode_mass_qc_table(
                self.results,
                enabled=self.qc_enable_in.isChecked(),
                tolerance_value=self.qc_tol_in.value(),
                tolerance_unit=self.qc_unit_in.currentText(),
            )
            self._update_tables()
            self._update_plots()
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Analysis failed",
                f"{exc}\n\n{traceback.format_exc(limit=6)}",
            )
        finally:
            QApplication.restoreOverrideCursor()

    def _combined_annotated(self):
        frames = [
            result.get("annotated")
            for result in self.results.values()
            if result.get("annotated") is not None and len(result.get("annotated"))
        ]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def _update_tables(self):
        summary = pc.build_summary_table(self.results, self.awdi_basis_in.currentData())
        populate_table(self.summary_table, summary)
        populate_table(self.mode_qc_table_w, self.qc_table)

        annotated = self._combined_annotated()
        preview_columns = [
            "Mass",
            "Apex.RT",
            "Fractional.Abundance",
            "TailBlockID",
            "DisplayLabel",
            "SpeciesClassLabel",
            "A_Count",
            "ReportLength",
            "ReportLengthBasis",
            "C_Count",
            "U_Count",
            "VariantPlacementLabel",
            "FivePrimeTerminus",
            "ThreePrimeTerminus",
            "Termini",
            "TheoreticalMass",
            "MassErrorDa",
            "MassErrorPpm",
            "AnnotationStatus",
            "PhosphateOffsetFlag",
            "PhosphateOffsetMessage",
        ]
        preview_columns = [column for column in preview_columns if column in annotated.columns]
        populate_table(self.annot_table, annotated[preview_columns] if len(annotated) else annotated, max_rows=200)

        library = self.generated.SpeciesLibrary if self.generated is not None else pd.DataFrame()
        library_columns = [
            "TailBlockID",
            "DisplayLabel",
            "SpeciesClassLabel",
            "A_Count",
            "C_Count",
            "U_Count",
            "VariantPlacementLabel",
            "FivePrimeTerminus",
            "ThreePrimeTerminus",
            "Termini",
            "TheoreticalMass",
            "GeneratedSequence",
        ]
        library_columns = [column for column in library_columns if column in library.columns]
        populate_table(self.lib_table, library[library_columns] if len(library) else library, max_rows=200)
        populate_table(
            self.tail_block_table_w,
            self.generated.TailBlocks if self.generated is not None else pd.DataFrame(),
            max_rows=100,
        )

        if self.qc_enable_in.isChecked() and not pc.mode_qc_all_pass(self.qc_table):
            self.qc_banner.setText(pc.mode_qc_failure_message(self.qc_table))
            self.qc_banner.setStyleSheet(
                "background:#f8d7da; border:1px solid #f5c2c7; color:#842029; padding:8px; font-weight:bold;"
            )
            self.qc_banner.setVisible(True)
        else:
            self.qc_banner.setVisible(False)

    def _update_plots(self):
        annotated = self._combined_annotated()
        if len(annotated) == 0:
            for figure, canvas in (
                (self.mass_fig, self.mass_canvas),
                (self.class_fig, self.class_canvas),
                (self.total_fig, self.total_canvas),
            ):
                figure.clear()
                _set_empty(figure.add_subplot(111), "No filtered data to plot.")
                canvas.draw()
            return

        draw_mass_plot(self.mass_fig, annotated)
        self.mass_canvas.draw()

        if self.qc_enable_in.isChecked() and not pc.mode_qc_all_pass(self.qc_table):
            message = pc.mode_qc_failure_message(self.qc_table)
            for figure, canvas in (
                (self.class_fig, self.class_canvas),
                (self.total_fig, self.total_canvas),
            ):
                figure.clear()
                _set_empty(figure.add_subplot(111), message)
                canvas.draw()
            return

        active_classes = []
        if self.show_a_in.isChecked():
            active_classes.append("A_only")
        if self.show_c_in.isChecked():
            active_classes.append("C_containing")
        if self.show_u_in.isChecked():
            active_classes.append("U_containing")
        if self.show_cu_in.isChecked():
            active_classes.append("mixed_CU")

        class_data = pc.build_class_plot_data(
            annotated, duplicate_rule=self.dup_rule_in.currentData()
        )
        if len(class_data) and active_classes:
            class_data = class_data[class_data["SpeciesClass"].isin(active_classes)]
        axis = self._axis_settings()
        draw_class_plot(
            self.class_fig,
            class_data,
            x_limits=axis["limits"],
            x_break_by=axis["break_by"],
            x_label_angle=axis["label_angle"],
        )
        self.class_canvas.draw()

        class_total = pc.build_class_plot_data(
            annotated,
            include_unassigned=True,
            duplicate_rule=self.dup_rule_in.currentData(),
        )
        total_data = pc.build_total_plot_data_from_class(
            class_total, threshold_pct=self.threshold_in.value()
        )
        draw_total_plot(
            self.total_fig,
            total_data,
            threshold_pct=self.threshold_in.value(),
            x_limits=axis["limits"],
            x_break_by=axis["break_by"],
            x_label_angle=axis["label_angle"],
        )
        self.total_canvas.draw()

    # ---------------------------------------------------------------- export

    def _check_results(self):
        if not self.results:
            QMessageBox.warning(self, "No results", "Run 'Annotate and Plot' first.")
            return False
        return True

    def _stem(self):
        return pc.safe_sample_id(self.sample_id_in.text() or "Sample")

    def _settings_dataframe(self):
        if self.generated is None:
            return pd.DataFrame()
        settings = {
            "AppVersion": pc.APP_VERSION,
            "GeneratedDate": str(date.today()),
            "SampleID": self.sample_id_in.text(),
            "OriginalInputSequence": self.generated.OriginalInputSequence,
            "InputSequenceType": self.generated.InputSequenceType,
            "ConvertedRNASequence": self.generated.ConvertedRNASequence,
            "SequenceMode": self.sequence_mode_in.currentData(),
            "Chemistry": self.chem_in.currentData(),
            "ExtraUResidueOffsetDa": self.u_shift_in.value(),
            "FivePrimeTerminusChoice": self.five_prime_in.currentData(),
            "ThreePrimeTerminusChoice": self.three_prime_in.currentData(),
            "VariantPlacement": self.variant_placement_in.currentData(),
            "RangeMin": self.range_min_in.value(),
            "RangeMax": self.range_max_in.value(),
            "RangeBasis": self.range_basis_in.currentData(),
            "MaxC": self._effective_max_c(),
            "MaxU": self._effective_max_u(),
            "AllowMixedCU": self.allow_mixed_in.isChecked(),
            "ToleranceUnit": self.tol_unit_in.currentText(),
            "MassTolerance": self.tol_value_in.value(),
            "FillToleranceDa": self.fill_tol_in.value(),
            "LengthOnlyToleranceDa": self.length_tol_in.value(),
            "ReportLengthBasis": self.report_basis_in.currentData(),
            "ManualShift": self.manual_shift_in.value(),
            "IndividualBinThresholdPct": self.threshold_in.value(),
            "PlotExportFormat": self.plot_format_in.currentData(),
            "PlotExportDPI": EXPORT_DPI,
            "SingleSamplePlotWidthIn": EXPORT_SINGLE_FIGSIZE[0],
            "SingleSamplePlotHeightIn": EXPORT_SINGLE_FIGSIZE[1],
            "AxisLabelFont": f"{FONT_FAMILY} Bold {AXIS_LABEL_SIZE} pt",
            "TickLabelFont": f"{FONT_FAMILY} Bold {TICK_LABEL_SIZE} pt",
            "JPEGQuality": self.jpeg_quality_in.value(),
        }
        return pd.DataFrame({"Setting": list(settings.keys()), "Value": list(settings.values())})

    def on_export_annotated_csv(self):
        if not self._check_results():
            return
        annotated = self._combined_annotated()
        if len(annotated) == 0:
            QMessageBox.information(self, "Nothing to save", "No annotated data.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save annotated data",
            f"{self._stem()}_PolyA_annotated_data_{date.today()}.csv",
            "CSV (*.csv)",
        )
        if path:
            annotated.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_lib_csv(self):
        if self.generated is None:
            QMessageBox.warning(self, "No library", "Run 'Annotate and Plot' first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save species library",
            f"{self._stem()}_PolyA_generated_internal_standards_{date.today()}.csv",
            "CSV (*.csv)",
        )
        if path:
            self.generated.SpeciesLibrary.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_digest_csv(self):
        if self.generated is None:
            QMessageBox.warning(self, "No digest", "Run 'Annotate and Plot' first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save T1 digest",
            f"{self._stem()}_PolyA_T1_digest_table_{date.today()}.csv",
            "CSV (*.csv)",
        )
        if path:
            self.generated.DigestTable.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_excel(self):
        if not self._check_results():
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save summary + data",
            f"{self._stem()}_PolyA_summary_and_data_{date.today()}.xlsx",
            "Excel (*.xlsx)",
        )
        if not path:
            return
        annotated = self._combined_annotated()
        summary = pc.build_summary_table(self.results, self.awdi_basis_in.currentData())
        class_data = (
            pc.build_class_plot_data(annotated, duplicate_rule=self.dup_rule_in.currentData())
            if len(annotated)
            else pd.DataFrame()
        )
        class_total = (
            pc.build_class_plot_data(
                annotated,
                include_unassigned=True,
                duplicate_rule=self.dup_rule_in.currentData(),
            )
            if len(annotated)
            else pd.DataFrame()
        )
        total_data = (
            pc.build_total_plot_data_from_class(class_total, threshold_pct=self.threshold_in.value())
            if len(class_total)
            else pd.DataFrame()
        )
        sequence_conversion = pc.sequence_conversion_table(
            self.sequence_in.toPlainText(), self.input_molecule_in.currentData()
        )
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            summary.to_excel(writer, sheet_name="Summary", index=False)
            if self.qc_table is not None:
                self.qc_table.to_excel(writer, sheet_name="Mode_Mass_QC", index=False)
            annotated.to_excel(writer, sheet_name="Annotated_Data", index=False)
            self.generated.SpeciesLibrary.to_excel(writer, sheet_name="Species_Library", index=False)
            self.generated.DigestTable.to_excel(writer, sheet_name="T1_Digest", index=False)
            self.generated.TailBlocks.to_excel(writer, sheet_name="Tail_Blocks", index=False)
            sequence_conversion.to_excel(writer, sheet_name="Sequence_Conversion", index=False)
            if len(class_data):
                class_data.to_excel(writer, sheet_name="Plot_Data_Class", index=False)
            if len(total_data):
                total_data.to_excel(writer, sheet_name="Plot_Data_Total", index=False)
            self._settings_dataframe().to_excel(writer, sheet_name="Filter_Settings", index=False)
        QMessageBox.information(self, "Saved", path)

    def _build_publication_export_figures(self):
        """Redraw plots on fixed-size canvases for reproducible exports.

        Saving the live Qt figures makes the physical figure size depend on the
        user's window and screen geometry. A fresh Figure gives every
        single-sample export the same 7 x 5 inch, 600 dpi specification used by
        the Cas9 reference plots.
        """
        annotated = self._combined_annotated()
        if len(annotated) == 0:
            return []

        sample_count = (
            int(annotated["SampleLabel"].nunique())
            if "SampleLabel" in annotated.columns
            else 1
        )
        figure_size = _publication_figsize(sample_count)
        axis = self._axis_settings()

        mass_figure = Figure(figsize=figure_size, dpi=100, facecolor="white")
        draw_mass_plot(mass_figure, annotated)

        active_classes = []
        if self.show_a_in.isChecked():
            active_classes.append("A_only")
        if self.show_c_in.isChecked():
            active_classes.append("C_containing")
        if self.show_u_in.isChecked():
            active_classes.append("U_containing")
        if self.show_cu_in.isChecked():
            active_classes.append("mixed_CU")

        class_data = pc.build_class_plot_data(
            annotated, duplicate_rule=self.dup_rule_in.currentData()
        )
        if len(class_data) and active_classes:
            class_data = class_data[class_data["SpeciesClass"].isin(active_classes)]

        class_figure = Figure(figsize=figure_size, dpi=100, facecolor="white")
        draw_class_plot(
            class_figure,
            class_data,
            x_limits=axis["limits"],
            x_break_by=axis["break_by"],
            x_label_angle=axis["label_angle"],
        )

        class_total = pc.build_class_plot_data(
            annotated,
            include_unassigned=True,
            duplicate_rule=self.dup_rule_in.currentData(),
        )
        total_data = pc.build_total_plot_data_from_class(
            class_total, threshold_pct=self.threshold_in.value()
        )
        total_figure = Figure(figsize=figure_size, dpi=100, facecolor="white")
        draw_total_plot(
            total_figure,
            total_data,
            threshold_pct=self.threshold_in.value(),
            x_limits=axis["limits"],
            x_break_by=axis["break_by"],
            x_label_angle=axis["label_angle"],
        )

        return [
            (mass_figure, "Mass_annotation_plot"),
            (class_figure, "Class_resolved_length_plot"),
            (total_figure, "Total_distribution_plot"),
        ]

    def _save_figure(self, figure: Figure, path: str):
        """Save one fixed-size publication figure at 600 dpi."""
        image_format = self.plot_format_in.currentData()
        common = {
            "dpi": EXPORT_DPI,
            "facecolor": "white",
            "edgecolor": "white",
        }
        if image_format == "jpeg":
            figure.savefig(
                path,
                format="jpeg",
                pil_kwargs={
                    "quality": self.jpeg_quality_in.value(),
                    "subsampling": 0,
                    "optimize": True,
                },
                **common,
            )
        else:
            figure.savefig(path, format="tiff", **common)

    def on_export_plots(self):
        if not self._check_results():
            return
        if self.qc_enable_in.isChecked() and not pc.mode_qc_all_pass(self.qc_table):
            QMessageBox.warning(
                self,
                "QC failed",
                "Mode-mass QC failed. Plot export is blocked until the mode peak matches the generated ladder.",
            )
            return

        image_format = self.plot_format_in.currentData()
        label = "TIFF" if image_format == "tiff" else "JPEG"
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save publication plots ZIP",
            f"{self._stem()}_PolyA_plots_{label}_{date.today()}.zip",
            "ZIP (*.zip)",
        )
        if not path:
            return

        extension = "tiff" if image_format == "tiff" else "jpg"
        figure_items = self._build_publication_export_figures()
        if not figure_items:
            QMessageBox.information(self, "Nothing to save", "No plot data are available.")
            return

        with tempfile.TemporaryDirectory() as temporary_directory:
            base = self._stem()
            paths = []
            for figure, plot_name in figure_items:
                filename = f"{base}_{plot_name}.{extension}"
                output_path = os.path.join(temporary_directory, filename)
                self._save_figure(figure, output_path)
                paths.append(output_path)
                figure.clear()

            with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
                for output_path in paths:
                    archive.write(output_path, os.path.basename(output_path))

        QMessageBox.information(
            self,
            "Saved",
            f"{path}\n\nPublication format: 7 × 5 inches per single-sample plot, {EXPORT_DPI} dpi.",
        )


def main():
    app = QApplication(sys.argv)
    app.setFont(QFont("Arial", 10))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
