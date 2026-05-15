"""
Poly(A) Heterogeneity Analysis Tool - desktop GUI.

PySide6 + matplotlib port of the Shiny app. All scientific behaviour lives in polya_core.
"""

from __future__ import annotations

import os
import sys
import io
import zipfile
import tempfile
import traceback
from datetime import date
from typing import Optional, Dict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("QtAgg")
# Font fallback - use Arial on Windows, otherwise fall through to DejaVu Sans
matplotlib.rcParams["font.family"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.weight"] = "bold"
import logging
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QGroupBox, QLabel, QLineEdit, QTextEdit, QComboBox, QCheckBox, QSpinBox,
    QDoubleSpinBox, QPushButton, QTabWidget, QTableWidget, QTableWidgetItem,
    QFileDialog, QMessageBox, QScrollArea, QHeaderView, QFormLayout, QFrame,
    QSizePolicy,
)

import polya_core as pc


# ---------------------------------------------------------------------------
# Plotting (matplotlib) - mirrors the ggplot2 report style
# ---------------------------------------------------------------------------

PLOT_FONT = {"family": "Arial", "weight": "bold", "color": "black"}


def _style_axes(ax, x_label_angle: float = 0):
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.85)
    ax.tick_params(axis="both", colors="black", width=0.45, length=4)
    for lab in ax.get_xticklabels():
        lab.set_fontfamily("Arial")
        lab.set_fontweight("bold")
        lab.set_color("black")
        if x_label_angle:
            lab.set_rotation(x_label_angle)
            lab.set_ha("right" if x_label_angle != 0 else "center")
            lab.set_va("center" if x_label_angle >= 80 else "top")
    for lab in ax.get_yticklabels():
        lab.set_fontfamily("Arial")
        lab.set_fontweight("bold")
        lab.set_color("black")
    ax.xaxis.label.set_fontproperties({"family": "Arial", "weight": "bold", "size": 12})
    ax.yaxis.label.set_fontproperties({"family": "Arial", "weight": "bold", "size": 12})
    ax.grid(False)


def _length_breaks(x, x_limits=None, break_by=None):
    if x_limits is not None and len(x_limits) == 2 and all(np.isfinite(x_limits)):
        lo = int(np.floor(min(x_limits)))
        hi = int(np.ceil(max(x_limits)))
    else:
        xv = np.asarray([v for v in x if pd.notna(v) and np.isfinite(v)])
        if len(xv) == 0:
            return None
        lo = int(np.floor(xv.min()))
        hi = int(np.ceil(xv.max()))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi < lo:
        return None
    span = hi - lo
    if break_by is not None and np.isfinite(break_by) and break_by > 0:
        by = int(break_by)
    else:
        by = 2 if span <= 30 else (5 if span <= 70 else 10)
    return list(range(lo, hi + 1, by))


def _length_limits(x):
    xv = np.asarray([v for v in x if pd.notna(v) and np.isfinite(v)])
    if len(xv) == 0:
        return None
    lo, hi = int(np.floor(xv.min())), int(np.ceil(xv.max()))
    span = hi - lo
    if span <= 25:
        return (lo, hi)
    return (int(np.floor((lo - 2) / 5) * 5), int(np.ceil((hi + 3) / 5) * 5))


def _plot_limits_from_data(data, length_col="PlotLength", custom_limits=None):
    if custom_limits is not None and len(custom_limits) == 2 and all(np.isfinite(custom_limits)) and custom_limits[1] > custom_limits[0]:
        return tuple(custom_limits)
    if data is None or len(data) == 0 or length_col not in data.columns:
        return None
    if "PlotMin" in data.columns and "PlotMax" in data.columns:
        lo = pd.to_numeric(data["PlotMin"], errors="coerce").dropna()
        hi = pd.to_numeric(data["PlotMax"], errors="coerce").dropna()
        if len(lo) and len(hi):
            lo_v, hi_v = float(lo.iloc[0]), float(hi.iloc[0])
            if np.isfinite(lo_v) and np.isfinite(hi_v) and hi_v > lo_v:
                return (lo_v, hi_v)
    return _length_limits(data[length_col])


def _set_empty(ax, message: str):
    ax.clear()
    ax.text(0.5, 0.5, message, ha="center", va="center",
            family="Arial", fontweight="bold", color="black", fontsize=11, wrap=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def draw_mass_plot(fig: Figure, plot_data: pd.DataFrame):
    fig.clear()
    if plot_data is None or len(plot_data) == 0:
        ax = fig.add_subplot(111)
        _set_empty(ax, "No filtered data to plot.")
        return

    class_order = ["A_only", "C_containing", "U_containing", "mixed_CU", "length_only", "unassigned"]
    samples = plot_data["SampleLabel"].unique() if "SampleLabel" in plot_data.columns else [None]

    n = len(samples)
    cols = 1 if n == 1 else 2
    rows = int(np.ceil(n / cols))
    axes = [fig.add_subplot(rows, cols, i + 1) for i in range(n)]

    legend_handles = []
    legend_labels = []
    seen = set()

    for ax, sample in zip(axes, samples):
        sub = plot_data if sample is None else plot_data[plot_data["SampleLabel"] == sample]
        present = [c for c in class_order if c in set(sub["SpeciesClass"].astype(str))]
        for cls in present:
            d = sub[sub["SpeciesClass"] == cls]
            ax.bar(d["Mass"], d["Fractional.Abundance"], width=0.1,
                   color=pc.CLASS_COLOURS.get(cls, "#999999"),
                   edgecolor="black", linewidth=0.25, zorder=2)
            if cls not in seen:
                seen.add(cls)
                lbl = pc.class_to_label(cls)
                legend_handles.append(Rectangle((0, 0), 1, 1,
                                                facecolor=pc.CLASS_COLOURS.get(cls, "#999999"),
                                                edgecolor="black", linewidth=0.5))
                legend_labels.append(lbl)
        ax.set_xlabel("Monoisotopic mass (Da)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if sample is not None and n > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax)

    if legend_handles:
        axes[0].legend(legend_handles, legend_labels, loc="upper right",
                       frameon=True, edgecolor="black", facecolor="white",
                       prop={"family": "Arial", "weight": "bold", "size": 9})
    fig.tight_layout()


def draw_class_plot(fig: Figure, class_data: pd.DataFrame,
                    x_limits=None, x_break_by=None, x_label_angle=45):
    fig.clear()
    if class_data is None or len(class_data) == 0:
        ax = fig.add_subplot(111)
        _set_empty(ax, "No class-resolved data after current filters.")
        return

    class_order = ["A_only", "C_containing", "U_containing", "mixed_CU"]
    present = [c for c in class_order if c in set(class_data["SpeciesClass"].astype(str))]

    samples = class_data["SampleLabel"].unique()
    n = len(samples)
    cols = 1 if n == 1 else 2
    rows = int(np.ceil(n / cols))
    axes = [fig.add_subplot(rows, cols, i + 1) for i in range(n)]

    legend_handles, legend_labels, seen = [], [], set()

    for ax, sample in zip(axes, samples):
        sub = class_data[class_data["SampleLabel"] == sample]
        if len(sub) == 0:
            _set_empty(ax, "No data for sample")
            continue
        lims = _plot_limits_from_data(sub, "PlotLength", custom_limits=x_limits) if n == 1 else x_limits
        # Side-by-side bars per length value
        n_classes = max(len(present), 1)
        bar_w = 0.74 / n_classes
        for i, cls in enumerate(present):
            d = sub[sub["SpeciesClass"] == cls]
            if len(d) == 0:
                continue
            offset = (i - (n_classes - 1) / 2) * bar_w
            ax.bar(d["PlotLength"] + offset, d["RelAbundance"], width=bar_w,
                   color=pc.CLASS_COLOURS.get(cls, "#999999"),
                   edgecolor="black", linewidth=0.25, zorder=2)
            if cls not in seen:
                seen.add(cls)
                legend_handles.append(Rectangle((0, 0), 1, 1,
                                                facecolor=pc.CLASS_COLOURS.get(cls, "#999999"),
                                                edgecolor="black", linewidth=0.5))
                legend_labels.append(pc.class_to_label(cls))
        ax.set_xlabel("Length (nt)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if lims is not None:
            ax.set_xlim(lims)
        breaks = _length_breaks(sub["PlotLength"], x_limits=lims, break_by=x_break_by)
        if breaks:
            ax.set_xticks(breaks)
        if n > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax, x_label_angle=x_label_angle)

    if legend_handles:
        axes[0].legend(legend_handles, legend_labels, loc="upper right",
                       frameon=True, edgecolor="black", facecolor="white",
                       prop={"family": "Arial", "weight": "bold", "size": 9})
    fig.tight_layout()


def draw_total_plot(fig: Figure, total_data: pd.DataFrame, threshold_pct=10,
                    x_limits=None, x_break_by=None, x_label_angle=45):
    fig.clear()
    if total_data is None or len(total_data) == 0:
        ax = fig.add_subplot(111)
        _set_empty(ax, "No total distribution data after current filters.")
        return

    samples = total_data["SampleLabel"].unique()
    n = len(samples)
    cols = 1 if n == 1 else 2
    rows = int(np.ceil(n / cols))
    axes = [fig.add_subplot(rows, cols, i + 1) for i in range(n)]

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=pc.CLASS_COLOURS["above_10"], edgecolor="black", linewidth=0.5),
        Rectangle((0, 0), 1, 1, facecolor=pc.CLASS_COLOURS["below_10"], edgecolor="black", linewidth=0.5),
    ]
    legend_labels = [f"\u2265{threshold_pct}%", f"<{threshold_pct}%"]

    for ax, sample in zip(axes, samples):
        sub = total_data[total_data["SampleLabel"] == sample].copy()
        if len(sub) == 0:
            _set_empty(ax, "No data")
            continue
        lims = _plot_limits_from_data(sub, "PlotLength", custom_limits=x_limits) if n == 1 else x_limits
        above = sub[sub["ThresholdClass"] == "above_10"]
        below = sub[sub["ThresholdClass"] == "below_10"]
        ax.bar(below["PlotLength"], below["TotalNorm"], width=0.85,
               color=pc.CLASS_COLOURS["below_10"], edgecolor="black", linewidth=0.25, zorder=2)
        ax.bar(above["PlotLength"], above["TotalNorm"], width=0.85,
               color=pc.CLASS_COLOURS["above_10"], edgecolor="black", linewidth=0.25, zorder=2)
        ax.set_xlabel("Length (nt)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylabel("Relative abundance (%)", fontdict=PLOT_FONT, fontsize=12)
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        if lims is not None:
            ax.set_xlim(lims)
        breaks = _length_breaks(sub["PlotLength"], x_limits=lims, break_by=x_break_by)
        if breaks:
            ax.set_xticks(breaks)
        # Stats annotation
        if sub["TotalRaw"].sum() > 0:
            mode_val = float(sub.loc[sub["TotalRaw"].idxmax(), "PlotLength"])
            mean_val = float((sub["PlotLength"] * sub["TotalRaw"]).sum() / sub["TotalRaw"].sum())
            left = sub[sub["PlotLength"] <= mode_val]
            right = sub[sub["PlotLength"] >= mode_val]
            sd_l = float(np.sqrt((left["TotalRaw"] * (left["PlotLength"] - mode_val) ** 2).sum() /
                                  left["TotalRaw"].sum())) if left["TotalRaw"].sum() > 0 else np.nan
            sd_r = float(np.sqrt((right["TotalRaw"] * (right["PlotLength"] - mode_val) ** 2).sum() /
                                  right["TotalRaw"].sum())) if right["TotalRaw"].sum() > 0 else np.nan
            stats_text = (f"Mode = {round(mode_val)} nt\nMean = {round(mean_val)} nt\n"
                          f"SD-L = {round(sd_l)} nt\nSD-R = {round(sd_r)} nt")
            x_pos = (lims[1] - 0.04 * (lims[1] - lims[0])) if lims else sub["PlotLength"].max()
            ax.text(x_pos, 82, stats_text, ha="right", va="top",
                    family="Arial", fontweight="bold", fontsize=9,
                    bbox=dict(facecolor="white", edgecolor="black", linewidth=0.35, pad=3))
        if n > 1:
            ax.set_title(str(sample), fontdict=PLOT_FONT, fontsize=11)
        _style_axes(ax, x_label_angle=x_label_angle)

    axes[0].legend(legend_handles, legend_labels, loc="upper right",
                   frameon=True, edgecolor="black", facecolor="white",
                   prop={"family": "Arial", "weight": "bold", "size": 9})
    fig.tight_layout()


# ---------------------------------------------------------------------------
# Helper: pandas -> QTableWidget
# ---------------------------------------------------------------------------

def populate_table(widget: QTableWidget, df: pd.DataFrame, max_rows: Optional[int] = None):
    if df is None or len(df) == 0:
        widget.clear()
        widget.setRowCount(1)
        widget.setColumnCount(1)
        widget.setHorizontalHeaderLabels(["Result"])
        widget.setItem(0, 0, QTableWidgetItem("No data to display"))
        return
    if max_rows is not None:
        df = df.head(max_rows)
    widget.clear()
    widget.setRowCount(len(df))
    widget.setColumnCount(len(df.columns))
    widget.setHorizontalHeaderLabels([str(c) for c in df.columns])
    for i in range(len(df)):
        for j, col in enumerate(df.columns):
            val = df.iloc[i, j]
            if pd.isna(val):
                txt = ""
            elif isinstance(val, float):
                txt = f"{val:.6g}"
            else:
                txt = str(val)
            item = QTableWidgetItem(txt)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            widget.setItem(i, j, item)
    widget.resizeColumnsToContents()


# ---------------------------------------------------------------------------
# Main Window
# ---------------------------------------------------------------------------


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Poly(A) Heterogeneity Analysis Tool")
        self.resize(1500, 950)

        self.uploaded_files: list[str] = []
        self.generated: Optional[pc.GeneratedLibrary] = None
        self.results: Dict[str, dict] = {}
        self.qc_table: Optional[pd.DataFrame] = None

        self._build_ui()
        self._wire_demo_loader()
        self._apply_demo("SR22")  # initialise defaults

    # ---------------------------------------------------------------- UI build

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer = QHBoxLayout(central)
        outer.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Horizontal)
        outer.addWidget(splitter)

        # ---------- left side: scrollable controls ----------
        controls_holder = QScrollArea()
        controls_holder.setWidgetResizable(True)
        controls_holder.setMinimumWidth(440)
        controls_inner = QWidget()
        controls_holder.setWidget(controls_inner)
        v = QVBoxLayout(controls_inner)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(8)

        v.addWidget(self._build_files_group())
        v.addWidget(self._build_sequence_group())
        v.addWidget(self._build_range_group())
        v.addWidget(self._build_cu_group())
        v.addWidget(self._build_filter_group())
        v.addWidget(self._build_match_group())
        v.addWidget(self._build_report_length_group())
        v.addWidget(self._build_axis_group())
        v.addWidget(self._build_awdi_group())

        self.warnings_label = QLabel()
        self.warnings_label.setWordWrap(True)
        self.warnings_label.setStyleSheet(
            "background:#fff3cd; border:1px solid #ffecb5; color:#664d03; "
            "padding:8px; font-weight:bold;"
        )
        self.warnings_label.setVisible(False)
        v.addWidget(self.warnings_label)

        self.run_btn = QPushButton("Annotate and Plot")
        self.run_btn.setStyleSheet(
            "QPushButton { background:#2563eb; color:white; font-weight:bold; padding:10px; }"
            "QPushButton:hover { background:#1d4ed8; }"
        )
        self.run_btn.clicked.connect(self.on_run)
        v.addWidget(self.run_btn)

        v.addWidget(self._build_export_group())
        v.addStretch(1)

        splitter.addWidget(controls_holder)

        # ---------- right side: tabs ----------
        self.tabs = QTabWidget()
        splitter.addWidget(self.tabs)
        splitter.setSizes([460, 1040])

        self._build_summary_tab()
        self._build_mass_tab()
        self._build_class_tab()
        self._build_total_tab()

        # Hook up live update of warnings (PySide6 needs one type per findChildren call)
        for w in self.findChildren(QCheckBox):
            w.toggled.connect(self._refresh_warnings)
        for w in self.findChildren(QComboBox):
            w.currentIndexChanged.connect(self._refresh_warnings)
        for w in self.findChildren(QSpinBox):
            w.valueChanged.connect(self._refresh_warnings)
        for w in self.findChildren(QDoubleSpinBox):
            w.valueChanged.connect(self._refresh_warnings)
        for w in self.findChildren(QLineEdit):
            w.textChanged.connect(self._refresh_warnings)
        self._refresh_warnings()

    # ---------- group builders ----------

    def _build_files_group(self):
        g = QGroupBox("Sample files")
        l = QVBoxLayout(g)
        self.files_btn = QPushButton("Upload sample Excel files (.xlsx)")
        self.files_btn.clicked.connect(self.on_pick_files)
        l.addWidget(self.files_btn)
        self.files_label = QLabel("No files selected")
        self.files_label.setWordWrap(True)
        l.addWidget(self.files_label)
        return g

    def _build_sequence_group(self):
        g = QGroupBox("Internal standard generation")
        f = QFormLayout(g)
        f.setLabelAlignment(Qt.AlignLeft)

        self.sample_id_in = QLineEdit("SR22")
        f.addRow("Sample ID:", self.sample_id_in)

        self.sequence_in = QTextEdit()
        self.sequence_in.setMinimumHeight(110)
        self.sequence_in.setPlainText(pc.SR22_SEQ)
        f.addRow("Sequence (5'→3'):", self.sequence_in)

        self.sequence_mode_in = QComboBox()
        self.sequence_mode_in.addItem("Digest input using RNase T1 rules", "digest_after_G")
        self.sequence_mode_in.addItem("Input is already the T1 tail fragment", "tail_fragment")
        f.addRow("Sequence input type:", self.sequence_mode_in)

        self.demo_in = QComboBox()
        self.demo_in.addItem("Keep current", "none")
        self.demo_in.addItem("SR22", "SR22")
        self.demo_in.addItem("CSP split tail", "CSP")
        f.addRow("Load example sequence:", self.demo_in)

        self.chem_in = QComboBox()
        self.chem_in.addItem("Canonical U", "canonical_U")
        self.chem_in.addItem("N1-methylpseudouridine for U positions", "N1mPseudoU")
        f.addRow("Nucleotide chemistry:", self.chem_in)

        self.u_shift_in = QDoubleSpinBox()
        self.u_shift_in.setRange(-100, 100)
        self.u_shift_in.setDecimals(4)
        self.u_shift_in.setSingleStep(0.0001)
        self.u_shift_in.setValue(0)
        f.addRow("Advanced: extra U-residue offset (Da):", self.u_shift_in)
        return g

    def _build_range_group(self):
        g = QGroupBox("Range and detection")
        f = QFormLayout(g)

        self.auto_range_in = QCheckBox("Auto range around detected A-run(s)")
        f.addRow(self.auto_range_in)

        self.expected_a_in = QSpinBox()
        self.expected_a_in.setRange(0, 10000); self.expected_a_in.setValue(122)
        f.addRow("Expected poly(A) length:", self.expected_a_in)

        self.range_min_in = QSpinBox(); self.range_min_in.setRange(0, 10000); self.range_min_in.setValue(110)
        self.range_max_in = QSpinBox(); self.range_max_in.setRange(0, 10000); self.range_max_in.setValue(140)
        f.addRow("Range minimum:", self.range_min_in)
        f.addRow("Range maximum:", self.range_max_in)

        self.range_basis_in = QComboBox()
        self.range_basis_in.addItem("Report length used in plots", "report_length")
        self.range_basis_in.addItem("A-count only", "a_count")
        f.addRow("Min/max refers to:", self.range_basis_in)

        self.min_a_run_in = QSpinBox(); self.min_a_run_in.setRange(1, 1000); self.min_a_run_in.setValue(10)
        f.addRow("Min A-run for tail detection:", self.min_a_run_in)

        self.auto_left_in = QSpinBox(); self.auto_left_in.setRange(0, 1000); self.auto_left_in.setValue(10)
        self.auto_right_in = QSpinBox(); self.auto_right_in.setRange(0, 1000); self.auto_right_in.setValue(15)
        f.addRow("Auto range: nt below detected:", self.auto_left_in)
        f.addRow("Auto range: nt above detected:", self.auto_right_in)
        return g

    def _build_cu_group(self):
        g = QGroupBox("C/U species in generated library")
        f = QFormLayout(g)

        self.max_c_in = QSpinBox(); self.max_c_in.setRange(0, 10); self.max_c_in.setValue(1)
        self.max_u_in = QSpinBox(); self.max_u_in.setRange(0, 10); self.max_u_in.setValue(1)
        f.addRow("Max C residues:", self.max_c_in)
        f.addRow("Max U residues:", self.max_u_in)

        self.consider_c_in = QCheckBox("Consider C-containing species during annotation"); self.consider_c_in.setChecked(True)
        self.consider_u_in = QCheckBox("Consider U-containing species during annotation"); self.consider_u_in.setChecked(False)
        self.allow_mixed_in = QCheckBox("Allow mixed C/U species when both are enabled"); self.allow_mixed_in.setChecked(False)
        f.addRow(self.consider_c_in)
        f.addRow(self.consider_u_in)
        f.addRow(self.allow_mixed_in)

        plot_classes_row = QHBoxLayout()
        self.show_a_in = QCheckBox("Poly(A)"); self.show_a_in.setChecked(True)
        self.show_c_in = QCheckBox("with C"); self.show_c_in.setChecked(True)
        self.show_u_in = QCheckBox("with U"); self.show_u_in.setChecked(True)
        self.show_cu_in = QCheckBox("with C/U"); self.show_cu_in.setChecked(True)
        for cb in (self.show_a_in, self.show_c_in, self.show_u_in, self.show_cu_in):
            plot_classes_row.addWidget(cb)
        wrap = QWidget(); wrap.setLayout(plot_classes_row)
        f.addRow("Classes in plots:", wrap)
        return g

    def _build_filter_group(self):
        g = QGroupBox("Filter input peaks")
        f = QFormLayout(g)

        self.rt_min_in = QDoubleSpinBox(); self.rt_min_in.setRange(0, 1000); self.rt_min_in.setValue(0); self.rt_min_in.setDecimals(2)
        self.rt_max_in = QDoubleSpinBox(); self.rt_max_in.setRange(0, 1000); self.rt_max_in.setValue(100); self.rt_max_in.setDecimals(2)
        f.addRow("RT min (min):", self.rt_min_in)
        f.addRow("RT max (min):", self.rt_max_in)

        self.rel_min_in = QDoubleSpinBox(); self.rel_min_in.setRange(0, 100); self.rel_min_in.setValue(0)
        self.rel_max_in = QDoubleSpinBox(); self.rel_max_in.setRange(0, 100); self.rel_max_in.setValue(100)
        f.addRow("Rel. abundance min (%):", self.rel_min_in)
        f.addRow("Rel. abundance max (%):", self.rel_max_in)

        self.mass_min_in = QDoubleSpinBox(); self.mass_min_in.setRange(0, 1e7); self.mass_min_in.setValue(5000); self.mass_min_in.setDecimals(1)
        self.mass_max_in = QDoubleSpinBox(); self.mass_max_in.setRange(0, 1e7); self.mass_max_in.setValue(75000); self.mass_max_in.setDecimals(1)
        f.addRow("Mass min (Da):", self.mass_min_in)
        f.addRow("Mass max (Da):", self.mass_max_in)
        return g

    def _build_match_group(self):
        g = QGroupBox("Matching tolerances")
        f = QFormLayout(g)

        self.tol_unit_in = QComboBox(); self.tol_unit_in.addItems(["ppm", "Da"])
        f.addRow("Direct mass-match unit:", self.tol_unit_in)

        self.tol_value_in = QDoubleSpinBox(); self.tol_value_in.setRange(0.001, 10000); self.tol_value_in.setDecimals(3); self.tol_value_in.setValue(10)
        f.addRow("Direct mass-match tolerance:", self.tol_value_in)

        self.fill_in = QCheckBox("Nearest-library fill outside direct tolerance"); self.fill_in.setChecked(True)
        f.addRow(self.fill_in)

        self.fill_tol_in = QDoubleSpinBox(); self.fill_tol_in.setRange(0, 1000); self.fill_tol_in.setDecimals(2); self.fill_tol_in.setValue(1)
        f.addRow("Nearest-fill tolerance (Da):", self.fill_tol_in)

        self.length_tol_in = QDoubleSpinBox(); self.length_tol_in.setRange(0, 10000); self.length_tol_in.setDecimals(1); self.length_tol_in.setValue(165)
        f.addRow("Length-only tol. (Da, min 165):", self.length_tol_in)

        self.legacy_in = QCheckBox("Manual/report mode: nearest-ladder classification"); self.legacy_in.setChecked(True)
        f.addRow(self.legacy_in)

        self.qc_enable_in = QCheckBox("Mode-mass QC blocks plots if mode doesn't match ladder"); self.qc_enable_in.setChecked(True)
        f.addRow(self.qc_enable_in)

        self.qc_unit_in = QComboBox(); self.qc_unit_in.addItems(["Da", "ppm"])
        f.addRow("Mode QC unit:", self.qc_unit_in)

        self.qc_tol_in = QDoubleSpinBox(); self.qc_tol_in.setRange(0.001, 1000); self.qc_tol_in.setDecimals(3); self.qc_tol_in.setValue(1)
        f.addRow("Mode QC tolerance:", self.qc_tol_in)

        self.dup_rule_in = QComboBox()
        self.dup_rule_in.addItem("Use dominant entry", "dominant")
        self.dup_rule_in.addItem("Sum entries", "sum")
        f.addRow("Duplicate entries:", self.dup_rule_in)
        return g

    def _build_report_length_group(self):
        g = QGroupBox("Report length convention")
        f = QFormLayout(g)
        self.report_basis_in = QComboBox()
        self.report_basis_in.addItem("A-count only", "a_count")
        self.report_basis_in.addItem("A-count + fixed T1 context", "context")
        self.report_basis_in.addItem("Full generated fragment length", "full_fragment")
        self.report_basis_in.addItem("A-count + manual shift", "manual_shift")
        self.report_basis_in.setCurrentIndex(1)
        f.addRow("Length shown in plots:", self.report_basis_in)

        self.manual_shift_in = QSpinBox(); self.manual_shift_in.setRange(-500, 500); self.manual_shift_in.setValue(0)
        f.addRow("Manual shift (nt):", self.manual_shift_in)
        return g

    def _build_axis_group(self):
        g = QGroupBox("Length-plot display")
        f = QFormLayout(g)
        self.use_custom_x_in = QCheckBox("Manually set x-axis range")
        f.addRow(self.use_custom_x_in)
        self.xmin_in = QSpinBox(); self.xmin_in.setRange(0, 100000); self.xmin_in.setValue(50)
        self.xmax_in = QSpinBox(); self.xmax_in.setRange(1, 100000); self.xmax_in.setValue(100)
        f.addRow("X-axis minimum (nt):", self.xmin_in)
        f.addRow("X-axis maximum (nt):", self.xmax_in)
        self.xbreak_in = QSpinBox(); self.xbreak_in.setRange(1, 100); self.xbreak_in.setValue(1)
        f.addRow("X-axis tick interval (nt):", self.xbreak_in)
        self.xangle_in = QSpinBox(); self.xangle_in.setRange(0, 90); self.xangle_in.setValue(45)
        f.addRow("X-axis label angle (deg):", self.xangle_in)
        self.threshold_in = QSpinBox(); self.threshold_in.setRange(0, 100); self.threshold_in.setValue(10)
        f.addRow("Report plot threshold (%):", self.threshold_in)
        return g

    def _build_awdi_group(self):
        g = QGroupBox("AWDI")
        f = QFormLayout(g)
        self.awdi_basis_in = QComboBox()
        self.awdi_basis_in.addItem("A-only assigned species", "A_only")
        self.awdi_basis_in.addItem("All assigned species", "all_assigned")
        self.awdi_basis_in.addItem("A-only + C-containing species", "report_default")
        f.addRow("AWDI basis:", self.awdi_basis_in)
        return g

    def _build_export_group(self):
        g = QGroupBox("Export")
        l = QVBoxLayout(g)
        self.exp_excel_btn = QPushButton("Download summary + annotated (Excel)")
        self.exp_csv_btn = QPushButton("Download annotated data (CSV)")
        self.exp_lib_btn = QPushButton("Download species library (CSV)")
        self.exp_digest_btn = QPushButton("Download T1 digest table (CSV)")
        self.exp_plots_btn = QPushButton("Download plots (TIFF + ZIP)")
        for btn, fn in [(self.exp_excel_btn, self.on_export_excel),
                        (self.exp_csv_btn, self.on_export_annotated_csv),
                        (self.exp_lib_btn, self.on_export_lib_csv),
                        (self.exp_digest_btn, self.on_export_digest_csv),
                        (self.exp_plots_btn, self.on_export_plots)]:
            btn.clicked.connect(fn)
            l.addWidget(btn)
        return g

    # ---------- tab builders ----------

    def _build_summary_tab(self):
        w = QWidget(); v = QVBoxLayout(w)
        v.addWidget(QLabel("<b>Summary</b>"))
        self.summary_table = QTableWidget(); v.addWidget(self.summary_table, 2)
        v.addWidget(QLabel("<b>Mode-mass QC</b>"))
        self.mode_qc_table_w = QTableWidget(); v.addWidget(self.mode_qc_table_w, 2)
        self.qc_banner = QLabel(); self.qc_banner.setWordWrap(True); self.qc_banner.setVisible(False)
        v.addWidget(self.qc_banner)
        v.addWidget(QLabel("<b>Annotation preview</b>"))
        self.annot_table = QTableWidget(); v.addWidget(self.annot_table, 3)
        v.addWidget(QLabel("<b>Generated species (preview)</b>"))
        self.lib_table = QTableWidget(); v.addWidget(self.lib_table, 2)
        v.addWidget(QLabel("<b>Detected T1 tail blocks</b>"))
        self.tail_block_table_w = QTableWidget(); v.addWidget(self.tail_block_table_w, 1)
        self.tabs.addTab(w, "Summary")

    def _build_mass_tab(self):
        w = QWidget(); v = QVBoxLayout(w)
        self.mass_fig = Figure(figsize=(8, 5), dpi=100)
        self.mass_canvas = FigureCanvas(self.mass_fig)
        v.addWidget(self.mass_canvas)
        self.tabs.addTab(w, "Observed mass distribution")

    def _build_class_tab(self):
        w = QWidget(); v = QVBoxLayout(w)
        self.class_fig = Figure(figsize=(8, 5), dpi=100)
        self.class_canvas = FigureCanvas(self.class_fig)
        v.addWidget(self.class_canvas)
        self.tabs.addTab(w, "Class-resolved tail distribution")

    def _build_total_tab(self):
        w = QWidget(); v = QVBoxLayout(w)
        self.total_fig = Figure(figsize=(8, 5), dpi=100)
        self.total_canvas = FigureCanvas(self.total_fig)
        v.addWidget(self.total_canvas)
        self.tabs.addTab(w, "Total report-style distribution")

    # ---------------------------------------------------------------- helpers

    def _wire_demo_loader(self):
        self.demo_in.currentIndexChanged.connect(lambda _: self._apply_demo(self.demo_in.currentData()))

    def _apply_demo(self, key):
        if key == "SR22":
            self.sample_id_in.setText("SR22")
            self.sequence_in.setPlainText(pc.SR22_SEQ)
            self.sequence_mode_in.setCurrentIndex(0)
            self.auto_range_in.setChecked(False)
            self.expected_a_in.setValue(122)
            self.range_min_in.setValue(110); self.range_max_in.setValue(140)
        elif key == "CSP":
            self.sample_id_in.setText("CSP")
            self.sequence_in.setPlainText(pc.CSP_SEQ)
            self.sequence_mode_in.setCurrentIndex(0)
            self.auto_range_in.setChecked(True)
            self.expected_a_in.setValue(70)
            self.range_min_in.setValue(60); self.range_max_in.setValue(85)

    def _effective_max_c(self):
        return self.max_c_in.value() if self.consider_c_in.isChecked() else 0

    def _effective_max_u(self):
        return self.max_u_in.value() if self.consider_u_in.isChecked() else 0

    def _effective_allow_mixed(self):
        return self.allow_mixed_in.isChecked() and self.consider_c_in.isChecked() and self.consider_u_in.isChecked()

    def _axis_settings(self):
        use_custom = self.use_custom_x_in.isChecked()
        xmin, xmax = self.xmin_in.value(), self.xmax_in.value()
        valid = use_custom and xmax > xmin
        return {
            "limits": (xmin, xmax) if valid else None,
            "break_by": self.xbreak_in.value(),
            "label_angle": self.xangle_in.value(),
        }

    def _refresh_warnings(self, *args):
        warnings = []
        if self.auto_range_in.isChecked():
            warnings.append("Auto range is ON: explicit min/max are ignored.")
        if self.tol_unit_in.currentText() == "Da" and self.tol_value_in.value() > 2:
            warnings.append("Mass tolerance is in Da and >2 Da. ppm mode is usually safer.")
        if self.chem_in.currentData() == "N1mPseudoU" and abs(self.u_shift_in.value()) > 0:
            warnings.append("N1mΨ chemistry already includes the U shift. Keep extra offset at 0 unless calibrating.")
        if self.report_basis_in.currentData() == "a_count":
            warnings.append("Report plots use A-count only; legacy plots usually want 'A-count + fixed T1 context'.")
        if not self.consider_c_in.isChecked() and self.max_c_in.value() > 0:
            warnings.append("C candidates disabled; Max C is ignored.")
        if not self.consider_u_in.isChecked() and self.max_u_in.value() > 0:
            warnings.append("U candidates disabled; Max U is ignored.")
        if self.allow_mixed_in.isChecked() and (not self.consider_c_in.isChecked() or not self.consider_u_in.isChecked()):
            warnings.append("Mixed C/U requires both C and U toggles enabled.")
        if self.qc_enable_in.isChecked():
            warnings.append("Mode-mass QC is ON: plots blocked unless mode peak matches the generated ladder.")
        if self.use_custom_x_in.isChecked() and self.xmax_in.value() <= self.xmin_in.value():
            warnings.append("Custom x-axis: max must be larger than min. Falling back to automatic limits.")

        if warnings:
            self.warnings_label.setText("\u2022 " + "\n\u2022 ".join(warnings))
            self.warnings_label.setVisible(True)
        else:
            self.warnings_label.setVisible(False)

    # ---------------------------------------------------------------- actions

    def on_pick_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select Excel files",
                                                 filter="Excel Files (*.xlsx)")
        if files:
            self.uploaded_files = files
            self.files_label.setText(f"{len(files)} file(s) selected:\n" + "\n".join(os.path.basename(f) for f in files))

    def on_run(self):
        try:
            seq = self.sequence_in.toPlainText()
            sample_id = self.sample_id_in.text() or "Sample"

            self.generated = pc.generate_species_library(
                sample_id=sample_id,
                seq=seq,
                chemistry=self.chem_in.currentData(),
                u_mass_shift=self.u_shift_in.value(),
                sequence_mode=self.sequence_mode_in.currentData(),
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
                allow_mixed=self._effective_allow_mixed(),
            )

            if not self.uploaded_files:
                QMessageBox.warning(self, "No files", "No Excel files selected. Library generated; please upload sample files to annotate.")
                populate_table(self.lib_table, self.generated.SpeciesLibrary[[
                    "DisplayLabel", "SpeciesClassLabel", "A_Count", "C_Count", "U_Count",
                    "FixedContextLength", "GeneratedSequence", "TheoreticalMass"]], max_rows=200)
                populate_table(self.tail_block_table_w, self.generated.TailBlocks)
                return

            self.results = {}
            lib = self.generated.SpeciesLibrary
            for fp in self.uploaded_files:
                raw = pc.preprocess_uploaded_excel(fp)
                filtered = raw[
                    (raw["Apex.RT"] >= self.rt_min_in.value()) & (raw["Apex.RT"] <= self.rt_max_in.value()) &
                    (raw["Fractional.Abundance"] >= self.rel_min_in.value()) & (raw["Fractional.Abundance"] <= self.rel_max_in.value()) &
                    (raw["Mass"] >= self.mass_min_in.value()) & (raw["Mass"] <= self.mass_max_in.value())
                ].copy()
                if len(filtered) == 0:
                    self.results[os.path.basename(fp)] = {"raw": raw, "filtered": filtered, "annotated": None, "awdi": float("nan")}
                    continue
                ann = pc.match_to_library(
                    filtered=filtered, lib=lib,
                    tolerance_value=self.tol_value_in.value(),
                    tolerance_unit=self.tol_unit_in.currentText(),
                    fill_unassigned=self.fill_in.isChecked(),
                    fill_tolerance_da=self.fill_tol_in.value(),
                    length_only_tolerance_da=self.length_tol_in.value(),
                    legacy_nearest_ladder=self.legacy_in.isChecked(),
                )
                ann = pc.add_report_length_fields(ann,
                                                  report_length_basis=self.report_basis_in.currentData(),
                                                  manual_report_shift=self.manual_shift_in.value())
                ann["File"] = os.path.basename(fp)
                ann["SampleLabel"] = pc.clean_sample_label(fp)
                awdi = pc.calculate_awdi_from_annotated(ann, basis=self.awdi_basis_in.currentData())
                self.results[os.path.basename(fp)] = {"raw": raw, "filtered": filtered, "annotated": ann, "awdi": awdi}

            self._update_displays()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"{e}\n\n{traceback.format_exc()}")

    def _combined_annotated(self):
        frames = [r["annotated"] for r in self.results.values() if r.get("annotated") is not None]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def _update_displays(self):
        # QC
        self.qc_table = pc.build_mode_mass_qc_table(
            self.results, enabled=self.qc_enable_in.isChecked(),
            tolerance_value=self.qc_tol_in.value(),
            tolerance_unit=self.qc_unit_in.currentText(),
        )
        qc_pass = pc.mode_qc_all_pass(self.qc_table) if self.qc_enable_in.isChecked() else True

        # Summary
        summary = pc.build_summary_table(self.results, self.awdi_basis_in.currentData())
        populate_table(self.summary_table, summary)
        populate_table(self.mode_qc_table_w, self.qc_table)

        if self.qc_enable_in.isChecked():
            if qc_pass:
                self.qc_banner.setText("Mode-mass QC passed. Assignment plots enabled.")
                self.qc_banner.setStyleSheet("background:#d1e7dd; border:1px solid #badbcc; color:#0f5132; padding:8px; font-weight:bold;")
            else:
                self.qc_banner.setText("Mode-mass QC failed. Assignment plots blocked.\n" + pc.mode_qc_failure_message(self.qc_table))
                self.qc_banner.setStyleSheet("background:#f8d7da; border:1px solid #f5c2c7; color:#842029; padding:8px; font-weight:bold; white-space:pre-wrap;")
            self.qc_banner.setVisible(True)
        else:
            self.qc_banner.setVisible(False)

        # Tail blocks + library preview
        populate_table(self.tail_block_table_w, self.generated.TailBlocks)
        lib_preview = self.generated.SpeciesLibrary[[
            "SampleID", "TailBlockID", "DisplayLabel", "SpeciesClassLabel", "A_Count",
            "PrefixContextLength", "SuffixContextLength", "FixedContextLength",
            "C_Count", "U_Count", "GeneratedSequence", "Composition", "Termini",
            "TheoreticalMass", "PositionInterpretation"]]
        populate_table(self.lib_table, lib_preview, max_rows=100)

        # Annotation preview
        ann_all = self._combined_annotated()
        if len(ann_all):
            keep = ["File", "Mass", "Apex.RT", "Fractional.Abundance", "TailBlockID",
                    "DisplayLabel", "SpeciesClassLabel", "A_Count", "ReportLength",
                    "ReportLengthBasis", "C_Count", "U_Count", "TheoreticalMass",
                    "MassErrorDa", "MassErrorPpm", "AnnotationStatus", "CandidateCount"]
            keep = [c for c in keep if c in ann_all.columns]
            populate_table(self.annot_table, ann_all[keep], max_rows=60)
        else:
            populate_table(self.annot_table, pd.DataFrame())

        # Plots
        axis = self._axis_settings()
        if not qc_pass and self.qc_enable_in.isChecked():
            msg = "Mode-mass QC failed.\n" + pc.mode_qc_failure_message(self.qc_table)
            for fig, canvas in [(self.mass_fig, self.mass_canvas),
                                (self.class_fig, self.class_canvas),
                                (self.total_fig, self.total_canvas)]:
                fig.clear()
                ax = fig.add_subplot(111); _set_empty(ax, msg)
                canvas.draw()
            return

        if len(ann_all) == 0:
            for canvas in (self.mass_canvas, self.class_canvas, self.total_canvas):
                canvas.figure.clear()
                ax = canvas.figure.add_subplot(111); _set_empty(ax, "No filtered data to plot.")
                canvas.draw()
            return

        # Mass plot
        draw_mass_plot(self.mass_fig, ann_all)
        self.mass_canvas.draw()

        # Class plot
        active_classes = []
        if self.show_a_in.isChecked(): active_classes.append("A_only")
        if self.show_c_in.isChecked(): active_classes.append("C_containing")
        if self.show_u_in.isChecked(): active_classes.append("U_containing")
        if self.show_cu_in.isChecked(): active_classes.append("mixed_CU")
        cd = pc.build_class_plot_data(ann_all, duplicate_rule=self.dup_rule_in.currentData())
        if len(cd) > 0 and active_classes:
            cd = cd[cd["SpeciesClass"].isin(active_classes)]
        draw_class_plot(self.class_fig, cd, x_limits=axis["limits"], x_break_by=axis["break_by"], x_label_angle=axis["label_angle"])
        self.class_canvas.draw()

        # Total plot
        cd_total = pc.build_class_plot_data(ann_all, include_unassigned=True, duplicate_rule=self.dup_rule_in.currentData())
        td = pc.build_total_plot_data_from_class(cd_total, threshold_pct=self.threshold_in.value())
        draw_total_plot(self.total_fig, td, threshold_pct=self.threshold_in.value(),
                        x_limits=axis["limits"], x_break_by=axis["break_by"], x_label_angle=axis["label_angle"])
        self.total_canvas.draw()

    # ---------------------------------------------------------------- exports

    def _check_results(self) -> bool:
        if not self.results:
            QMessageBox.warning(self, "No results", "Run 'Annotate and Plot' first.")
            return False
        return True

    def _stem(self):
        return pc.safe_sample_id(self.sample_id_in.text() or "Sample")

    def on_export_annotated_csv(self):
        if not self._check_results():
            return
        ann = self._combined_annotated()
        if len(ann) == 0:
            QMessageBox.information(self, "Nothing to save", "No annotated data."); return
        path, _ = QFileDialog.getSaveFileName(self, "Save annotated data",
                                              f"{self._stem()}_PolyA_annotated_data_{date.today()}.csv",
                                              filter="CSV (*.csv)")
        if path:
            ann.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_lib_csv(self):
        if self.generated is None:
            QMessageBox.warning(self, "No library", "Run 'Annotate and Plot' first."); return
        path, _ = QFileDialog.getSaveFileName(self, "Save species library",
                                              f"{self._stem()}_PolyA_generated_internal_standards_{date.today()}.csv",
                                              filter="CSV (*.csv)")
        if path:
            self.generated.SpeciesLibrary.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_digest_csv(self):
        if self.generated is None:
            QMessageBox.warning(self, "No digest", "Run 'Annotate and Plot' first."); return
        path, _ = QFileDialog.getSaveFileName(self, "Save T1 digest",
                                              f"{self._stem()}_PolyA_T1_digest_table_{date.today()}.csv",
                                              filter="CSV (*.csv)")
        if path:
            self.generated.DigestTable.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", path)

    def on_export_excel(self):
        if not self._check_results():
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save summary + data",
                                              f"{self._stem()}_PolyA_summary_and_data_{date.today()}.xlsx",
                                              filter="Excel (*.xlsx)")
        if not path:
            return
        ann = self._combined_annotated()
        summary = pc.build_summary_table(self.results, self.awdi_basis_in.currentData())
        cd = pc.build_class_plot_data(ann, duplicate_rule=self.dup_rule_in.currentData()) if len(ann) else pd.DataFrame()
        cd_total = pc.build_class_plot_data(ann, include_unassigned=True, duplicate_rule=self.dup_rule_in.currentData()) if len(ann) else pd.DataFrame()
        td = pc.build_total_plot_data_from_class(cd_total, threshold_pct=self.threshold_in.value()) if len(cd_total) else pd.DataFrame()
        settings = pd.DataFrame({
            "Setting": ["SampleID", "InputSequence", "SequenceMode", "Chemistry",
                        "RangeMin", "RangeMax", "RangeBasis", "MaxC", "MaxU",
                        "ToleranceUnit", "MassTolerance", "FillToleranceDa",
                        "ReportLengthBasis", "ManualShift", "ThresholdPct"],
            "Value": [self.sample_id_in.text(), self.generated.InputSequence,
                      self.sequence_mode_in.currentData(), self.chem_in.currentData(),
                      self.range_min_in.value(), self.range_max_in.value(),
                      self.range_basis_in.currentData(),
                      self._effective_max_c(), self._effective_max_u(),
                      self.tol_unit_in.currentText(), self.tol_value_in.value(),
                      self.fill_tol_in.value(), self.report_basis_in.currentData(),
                      self.manual_shift_in.value(), self.threshold_in.value()],
        })
        with pd.ExcelWriter(path, engine="openpyxl") as xw:
            summary.to_excel(xw, sheet_name="Summary", index=False)
            if self.qc_table is not None:
                self.qc_table.to_excel(xw, sheet_name="Mode_Mass_QC", index=False)
            ann.to_excel(xw, sheet_name="Annotated_Data", index=False)
            self.generated.SpeciesLibrary.to_excel(xw, sheet_name="Species_Library", index=False)
            self.generated.DigestTable.to_excel(xw, sheet_name="T1_Digest", index=False)
            self.generated.TailBlocks.to_excel(xw, sheet_name="Tail_Blocks", index=False)
            if len(cd):
                cd.to_excel(xw, sheet_name="Plot_Data_Class", index=False)
            if len(td):
                td.to_excel(xw, sheet_name="Plot_Data_Total", index=False)
            settings.to_excel(xw, sheet_name="Filter_Settings", index=False)
        QMessageBox.information(self, "Saved", path)

    def on_export_plots(self):
        if not self._check_results():
            return
        if self.qc_enable_in.isChecked() and not pc.mode_qc_all_pass(self.qc_table):
            QMessageBox.warning(self, "QC failed",
                                "Mode-mass QC failed. Plot export is blocked until the mode peak matches the generated ladder.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save plots ZIP",
                                              f"{self._stem()}_PolyA_plots_{date.today()}.zip",
                                              filter="ZIP (*.zip)")
        if not path:
            return
        ann = self._combined_annotated()
        with tempfile.TemporaryDirectory() as tmp:
            base = self._stem()
            paths = []
            self.mass_fig.savefig(os.path.join(tmp, f"{base}_Mass_annotation_plot.tiff"),
                                   dpi=600, format="tiff", bbox_inches="tight")
            paths.append(os.path.join(tmp, f"{base}_Mass_annotation_plot.tiff"))
            self.class_fig.savefig(os.path.join(tmp, f"{base}_Class_resolved_length_plot.tiff"),
                                    dpi=600, format="tiff", bbox_inches="tight")
            paths.append(os.path.join(tmp, f"{base}_Class_resolved_length_plot.tiff"))
            self.total_fig.savefig(os.path.join(tmp, f"{base}_Total_distribution_plot.tiff"),
                                    dpi=600, format="tiff", bbox_inches="tight")
            paths.append(os.path.join(tmp, f"{base}_Total_distribution_plot.tiff"))
            with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
                for p in paths:
                    z.write(p, os.path.basename(p))
        QMessageBox.information(self, "Saved", path)


def main():
    app = QApplication(sys.argv)
    app.setFont(QFont("Arial", 10))
    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
