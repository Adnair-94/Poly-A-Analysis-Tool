"""
Poly(A) Heterogeneity Analysis Tool - core science module.

Python port of the R / Shiny workflow (T1 digest internal-library workflow FINAL v6).
All scientific behaviour mirrors the original R script:
  - RNase T1 cleaves after G.
  - Internal T1 fragments modelled as 5'-OH / 3'-phosphate.
  - Terminal T1 fragment modelled as 5'-OH / 3'-OH.
  - A_Count is the modelled continuous-A length in the detected tail block.
  - ReportLength is recomputed reactively from current settings.
  - "Poly(A) with C/U" is composition-level only; position is not inferred.

Mass conventions (monoisotopic):
  5'-OH / 3'-OH       = sum(residues) + H2O - HPO3
  5'-OH / 3'-phosphate = sum(residues) + H2O
"""

from __future__ import annotations

import os
import re
import unicodedata
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Visual constants (mirrors R CLASS_COLOURS / CLASS_LABELS)
# ---------------------------------------------------------------------------

CLASS_COLOURS = {
    "A_only":       "#008B8B",
    "C_containing": "#A64CE6",
    "U_containing": "#F97316",
    "mixed_CU":     "#64748B",
    "ambiguous":    "#9CA3AF",
    "length_only":  "#D9D9D9",
    "unassigned":   "#D9D9D9",
    "above_10":     "#008B8B",
    "below_10":     "#D9D9D9",
}

CLASS_LABELS = {
    "A_only":       "Poly(A)",
    "C_containing": "Poly(A) with C",
    "U_containing": "Poly(A) with U",
    "mixed_CU":     "Poly(A) with C/U",
    "ambiguous":    "Ambiguous",
    "length_only":  "Nearest length only",
    "unassigned":   "Unassigned",
}

# ---------------------------------------------------------------------------
# Mass constants
# ---------------------------------------------------------------------------

RNA_RESIDUE_MASS = {
    "A": 329.05252,
    "C": 305.04129,
    "G": 345.04744,
    "U": 306.02530,
}
HPO3_MASS = 79.9663309
H2O_MASS = 18.0105647
N1M_PSEUDOURIDINE_SHIFT = 14.0156501  # CH2 shift vs U

# ---------------------------------------------------------------------------
# Default validation sequences
# ---------------------------------------------------------------------------

SR22_SEQ = (
    "GUUUCUUCACAUUCUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAA"
)
CSP_SEQ = (
    "GCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCATATGACTAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGA"
)

# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def _normalise_names(cols):
    """Match R's tolower + replace punctuation/space runs with '.'."""
    out = []
    for c in cols:
        c2 = str(c).lower()
        c2 = re.sub(r"[\W_]+", ".", c2, flags=re.UNICODE)
        out.append(c2)
    return out


def _first_matching_col(nms, patterns):
    for p in patterns:
        rx = re.compile(p)
        for n in nms:
            if rx.search(n):
                return n
    return None


def clean_sample_label(path: str) -> str:
    base = os.path.splitext(os.path.basename(path))[0]
    base = re.sub(r"[-_.]+", " ", base)
    base = re.sub(r"\s+", " ", base).strip()
    return base


def safe_sample_id(x: Optional[str], default: str = "Sample") -> str:
    if x is None or not str(x).strip():
        return default
    s = unicodedata.normalize("NFKD", str(x)).encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s if s else default


def class_to_label(x):
    if isinstance(x, str):
        return CLASS_LABELS.get(x, x)
    return [CLASS_LABELS.get(v, v) for v in x]


def report_length_basis_label(basis: str) -> str:
    return {
        "a_count": "A-count only",
        "context": "A-count + fixed T1 context",
        "full_fragment": "Full generated fragment length",
        "manual_shift": "A-count + manual shift",
    }.get(str(basis), str(basis))


# ---------------------------------------------------------------------------
# BioPharma Finder file parsing
# ---------------------------------------------------------------------------


def preprocess_uploaded_excel(filepath: str) -> pd.DataFrame:
    df = pd.read_excel(filepath, sheet_name=0)
    df.columns = _normalise_names(df.columns)
    mass_col = _first_matching_col(df.columns, [r"^monoisotopic\.mass$", r"monoisotopic.*mass", r"^mass$"])
    rt_col = _first_matching_col(df.columns, [r"^apex\.rt$", r"apex.*rt", r"retention.*time", r"^rt$"])
    rel_col = _first_matching_col(df.columns, [r"^relative\.abundance$", r"relative.*abundance",
                                               r"fractional.*abundance", r"rel.*abund"])
    if mass_col is None or rt_col is None or rel_col is None:
        raise ValueError("Excel file must contain Monoisotopic Mass, Relative Abundance, and Apex RT columns.")
    out = pd.DataFrame({
        "Mass": pd.to_numeric(df[mass_col], errors="coerce"),
        "Apex.RT": pd.to_numeric(df[rt_col], errors="coerce"),
        "Fractional.Abundance": pd.to_numeric(df[rel_col], errors="coerce"),
    })
    return out.dropna().reset_index(drop=True)


# ---------------------------------------------------------------------------
# Sequence and mass logic
# ---------------------------------------------------------------------------


def normalise_rna_sequence(seq: str) -> str:
    if seq is None or not str(seq).strip():
        raise ValueError("Sequence input is empty.")
    s = str(seq).upper()
    s = re.sub(r"[^ACGTU]", "", s)
    s = s.replace("T", "U")
    if not s:
        raise ValueError("Sequence contains no valid A/C/G/U bases after cleanup.")
    return s


def composition_string(seq: str) -> str:
    parts = []
    for b in "ACGU":
        n = seq.count(b)
        if n > 0:
            parts.append(f"{b}{n}")
    return "".join(parts)


def mass_table_for_chemistry(chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> dict:
    m = dict(RNA_RESIDUE_MASS)
    if chemistry == "N1mPseudoU":
        m["U"] = m["U"] + N1M_PSEUDOURIDINE_SHIFT
    if u_mass_shift and np.isfinite(u_mass_shift):
        m["U"] = m["U"] + float(u_mass_shift)
    return m


def calculate_rna_mass(seq: str, terminal_type: str = "3OH",
                       chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> float:
    s = normalise_rna_sequence(seq)
    masses = mass_table_for_chemistry(chemistry, u_mass_shift)
    mass = sum(masses[b] for b in s)
    if terminal_type == "3p":
        return mass + H2O_MASS
    if terminal_type == "3OH":
        return mass + H2O_MASS - HPO3_MASS
    if terminal_type == "3cyclicp":
        return mass
    raise ValueError("Unsupported terminal_type. Use 3OH, 3p, or 3cyclicp.")


def find_a_runs(seq: str) -> pd.DataFrame:
    rows = []
    for m in re.finditer(r"A+", seq):
        rows.append({
            "RunStart": m.start() + 1,
            "RunEnd": m.end(),
            "RunLength": m.end() - m.start(),
        })
    return pd.DataFrame(rows, columns=["RunStart", "RunEnd", "RunLength"])


def digest_rnase_t1(seq: str, chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> pd.DataFrame:
    s = normalise_rna_sequence(seq)
    n = len(s)
    g_positions = [i + 1 for i, ch in enumerate(s) if ch == "G"]
    cut_ends = [p for p in g_positions if p < n]
    frag_ends = list(dict.fromkeys(cut_ends + [n])) or [n]

    rows = []
    start = 1
    for i, end in enumerate(frag_ends, start=1):
        if end < start:
            continue
        frag_seq = s[start - 1:end]
        terminal_type = "3p" if end < n and s[end - 1] == "G" else "3OH"
        mass = calculate_rna_mass(frag_seq, terminal_type, chemistry, u_mass_shift)
        runs = find_a_runs(frag_seq)
        max_run = int(runs["RunLength"].max()) if len(runs) else 0
        rows.append({
            "FragmentNumber": i,
            "FragmentID": f"T1_F{i}",
            "Start": start,
            "End": end,
            "FragmentSequence": frag_seq,
            "FragmentLength": len(frag_seq),
            "TerminalType": terminal_type,
            "Termini": "5'-OH / 3'-phosphate" if terminal_type == "3p" else "5'-OH / 3'-OH",
            "Composition": composition_string(frag_seq),
            "MaxContinuousA": max_run,
            "TheoreticalMass": round(mass, 6),
        })
        start = end + 1
    return pd.DataFrame(rows)


def digest_tail_fragment_input(seq: str, chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> pd.DataFrame:
    s = normalise_rna_sequence(seq)
    mass = calculate_rna_mass(s, "3OH", chemistry, u_mass_shift)
    runs = find_a_runs(s)
    max_run = int(runs["RunLength"].max()) if len(runs) else 0
    return pd.DataFrame([{
        "FragmentNumber": 1,
        "FragmentID": "T1_tail_input",
        "Start": 1,
        "End": len(s),
        "FragmentSequence": s,
        "FragmentLength": len(s),
        "TerminalType": "3OH",
        "Termini": "5'-OH / 3'-OH",
        "Composition": composition_string(s),
        "MaxContinuousA": max_run,
        "TheoreticalMass": round(mass, 6),
    }])


def detect_tail_blocks(digest_df: pd.DataFrame, min_a_run: int = 10) -> pd.DataFrame:
    rows = []
    for _, frag in digest_df.iterrows():
        frag_seq = frag["FragmentSequence"]
        runs = find_a_runs(frag_seq)
        if len(runs) == 0:
            continue
        runs = runs[runs["RunLength"] >= min_a_run]
        if len(runs) == 0:
            continue
        for j, (_, run) in enumerate(runs.iterrows(), start=1):
            rs, re_end = int(run["RunStart"]), int(run["RunEnd"])
            prefix = frag_seq[:rs - 1]
            suffix = frag_seq[re_end:]
            rows.append({
                "TailBlockNumber": len(rows) + 1,
                "TailBlockID": f"{frag['FragmentID']}_Ablock{j}",
                "SourceFragmentID": frag["FragmentID"],
                "SourceFragmentNumber": frag["FragmentNumber"],
                "SourceFragmentSequence": frag_seq,
                "SourceFragmentLength": frag["FragmentLength"],
                "TerminalType": frag["TerminalType"],
                "Termini": frag["Termini"],
                "RunStartInFragment": rs,
                "RunEndInFragment": re_end,
                "RunStartInInput": int(frag["Start"]) + rs - 1,
                "RunEndInInput": int(frag["Start"]) + re_end - 1,
                "DetectedALength": int(run["RunLength"]),
                "PrefixBeforeA": prefix,
                "SuffixAfterA": suffix,
            })
    return pd.DataFrame(rows)


def variant_class(c: int, u: int) -> str:
    if c == 0 and u == 0:
        return "A_only"
    if c > 0 and u == 0:
        return "C_containing"
    if c == 0 and u > 0:
        return "U_containing"
    return "mixed_CU"


def species_display_label(a: int, c: int, u: int) -> str:
    if c == 0 and u == 0:
        return f"A{a}"
    parts = []
    if c > 0:
        parts.append("C" if c == 1 else f"C{c}")
    if u > 0:
        parts.append("U" if u == 1 else f"U{u}")
    return f"A{a} with {'/'.join(parts)}"


def generate_variant_combos(max_c: int = 1, max_u: int = 1, allow_mixed: bool = True) -> pd.DataFrame:
    max_c = max(0, int(max_c))
    max_u = max(0, int(max_u))
    rows = []
    for c in range(max_c + 1):
        for u in range(max_u + 1):
            if not allow_mixed and c > 0 and u > 0:
                continue
            rows.append({"C_Count": c, "U_Count": u})
    df = pd.DataFrame(rows)
    if len(df):
        df["_order"] = df["C_Count"] + df["U_Count"]
        df = df.sort_values(["_order", "C_Count", "U_Count"]).drop(columns="_order").reset_index(drop=True)
    return df


@dataclass
class GeneratedLibrary:
    InputSequence: str
    DigestTable: pd.DataFrame
    TailBlocks: pd.DataFrame
    SpeciesLibrary: pd.DataFrame


def generate_species_library(
    sample_id: str,
    seq: str,
    chemistry: str = "canonical_U",
    u_mass_shift: float = 0.0,
    sequence_mode: str = "digest_after_G",
    use_auto_ranges: bool = False,
    expected_a_length: Optional[float] = None,
    range_min: Optional[float] = None,
    range_max: Optional[float] = None,
    range_input_basis: str = "report_length",
    report_length_basis: str = "context",
    manual_report_shift: float = 0.0,
    min_a_run: int = 10,
    auto_range_left: int = 10,
    auto_range_right: int = 15,
    max_c: int = 1,
    max_u: int = 1,
    allow_mixed: bool = True,
) -> GeneratedLibrary:
    clean_seq = normalise_rna_sequence(seq)
    if sequence_mode == "tail_fragment":
        digest_df = digest_tail_fragment_input(clean_seq, chemistry, u_mass_shift)
    else:
        digest_df = digest_rnase_t1(clean_seq, chemistry, u_mass_shift)

    tail_blocks = detect_tail_blocks(digest_df, min_a_run=min_a_run)
    if len(tail_blocks) == 0:
        raise ValueError("No tail-positive T1 fragment detected. Lower the minimum A-run threshold or check the sequence.")

    sample_id = (sample_id or "").strip() or "Sample"
    combos = generate_variant_combos(max_c, max_u, allow_mixed)

    rows = []
    block_rows = []
    for _, tb in tail_blocks.iterrows():
        detected_n = int(tb["DetectedALength"])
        prefix = tb["PrefixBeforeA"]
        suffix = tb["SuffixAfterA"]
        prefix_len = len(prefix)
        suffix_len = len(suffix)
        fixed_ctx = prefix_len + suffix_len

        rmin = pd.to_numeric(range_min, errors="coerce")
        rmax = pd.to_numeric(range_max, errors="coerce")
        rmin_ok = pd.notna(rmin)
        rmax_ok = pd.notna(rmax)

        if use_auto_ranges or not rmin_ok or not rmax_ok or rmax < rmin:
            lo = max(0, detected_n - auto_range_left)
            hi = detected_n + auto_range_right
            expected_n = detected_n
        else:
            basis = str(range_input_basis) if str(range_input_basis) in ("a_count", "report_length") else "report_length"
            range_shift = 0
            if basis == "report_length":
                rlb = str(report_length_basis) if str(report_length_basis) in ("a_count", "context", "full_fragment", "manual_shift") else "context"
                if rlb in ("context", "full_fragment"):
                    range_shift = fixed_ctx
                elif rlb == "manual_shift":
                    tmp = pd.to_numeric(manual_report_shift, errors="coerce")
                    range_shift = float(tmp) if pd.notna(tmp) else 0
                else:
                    range_shift = 0
            lo = max(0, int(np.floor(rmin - range_shift)))
            hi = max(lo, int(np.ceil(rmax - range_shift)))
            expected_n = int(expected_a_length) if expected_a_length is not None and np.isfinite(expected_a_length) else detected_n

        expected_seq = prefix + ("A" * expected_n) + suffix
        expected_mass = calculate_rna_mass(expected_seq, tb["TerminalType"], chemistry, u_mass_shift)

        block_rows.append({
            "SampleID": sample_id,
            "TailBlockID": tb["TailBlockID"],
            "SourceFragmentID": tb["SourceFragmentID"],
            "SourceFragmentSequence": tb["SourceFragmentSequence"],
            "DetectedALength": detected_n,
            "ExpectedALength": expected_n,
            "ARangeMin": lo,
            "ARangeMax": hi,
            "InputRangeMin": rmin if rmin_ok else np.nan,
            "InputRangeMax": rmax if rmax_ok else np.nan,
            "RangeInputBasis": "A-count" if str(range_input_basis) == "a_count" else "Report length",
            "ExpectedAOnlyMass": round(expected_mass, 6),
            "PrefixBeforeA": prefix,
            "SuffixAfterA": suffix,
            "PrefixContextLength": prefix_len,
            "SuffixContextLength": suffix_len,
            "FixedContextLength": fixed_ctx,
            "TerminalType": tb["TerminalType"],
            "Termini": tb["Termini"],
        })

        for a in range(lo, hi + 1):
            for _, combo in combos.iterrows():
                c = int(combo["C_Count"])
                u = int(combo["U_Count"])
                generated_seq = prefix + ("A" * a) + ("C" * c) + ("U" * u) + suffix
                mass = calculate_rna_mass(generated_seq, tb["TerminalType"], chemistry, u_mass_shift)
                cls = variant_class(c, u)
                rows.append({
                    "SampleID": sample_id,
                    "TailBlockID": tb["TailBlockID"],
                    "SourceFragmentID": tb["SourceFragmentID"],
                    "SpeciesID": f"{sample_id}_{tb['TailBlockID']}_A{a}_C{c}_U{u}",
                    "DisplayLabel": species_display_label(a, c, u),
                    "SpeciesClass": cls,
                    "SpeciesClassLabel": class_to_label(cls),
                    "A_Count": int(a),
                    "AssignedLength": int(a),
                    "C_Count": int(c),
                    "U_Count": int(u),
                    "NonA_Count": int(c + u),
                    "TotalTailModelResidues": int(a + c + u),
                    "PrefixContextLength": prefix_len,
                    "SuffixContextLength": suffix_len,
                    "FixedContextLength": fixed_ctx,
                    "GeneratedSequence": generated_seq,
                    "FragmentLength": len(generated_seq),
                    "Composition": composition_string(generated_seq),
                    "PositionInterpretation": (
                        "A-only" if (c + u) == 0
                        else "Composition-level only; C/U position not inferred"
                    ),
                    "TerminalType": tb["TerminalType"],
                    "Termini": tb["Termini"],
                    "NucleotideChemistry": chemistry,
                    "TheoreticalMass": round(mass, 6),
                    "ExpectedALength": expected_n,
                    "ExpectedAOnlyMass": round(expected_mass, 6),
                    "PlotMin": lo,
                    "PlotMax": hi,
                    "InputRangeMin": rmin if rmin_ok else np.nan,
                    "InputRangeMax": rmax if rmax_ok else np.nan,
                    "RangeInputBasis": "A-count" if str(range_input_basis) == "a_count" else "Report length",
                })

    return GeneratedLibrary(
        InputSequence=clean_seq,
        DigestTable=digest_df,
        TailBlocks=pd.DataFrame(block_rows),
        SpeciesLibrary=pd.DataFrame(rows),
    )


# ---------------------------------------------------------------------------
# Report-length helpers
# ---------------------------------------------------------------------------


def add_report_length_fields(annotated: pd.DataFrame,
                             report_length_basis: str = "context",
                             manual_report_shift: float = 0.0) -> pd.DataFrame:
    if annotated is None or len(annotated) == 0:
        return annotated
    d = annotated.copy()
    basis = str(report_length_basis) if str(report_length_basis) in (
        "a_count", "context", "full_fragment", "manual_shift") else "context"
    ms = pd.to_numeric(manual_report_shift, errors="coerce")
    ms = float(ms) if pd.notna(ms) else 0.0

    a_count = pd.to_numeric(d.get("A_Count"), errors="coerce")
    ctx = pd.to_numeric(d["FixedContextLength"], errors="coerce") if "FixedContextLength" in d else pd.Series(0, index=d.index)
    ctx = ctx.fillna(0)
    non_a = pd.to_numeric(d["NonA_Count"], errors="coerce") if "NonA_Count" in d else pd.Series(0, index=d.index)
    non_a = non_a.fillna(0)
    frag_len = pd.to_numeric(d["FragmentLength"], errors="coerce") if "FragmentLength" in d else pd.Series(np.nan, index=d.index)

    if basis == "a_count":
        report_length = a_count
        applied_shift = pd.Series(0, index=d.index)
    elif basis == "context":
        report_length = a_count + ctx
        applied_shift = ctx
    elif basis == "full_fragment":
        report_length = frag_len
        applied_shift = ctx + non_a
    else:
        report_length = a_count + ms
        applied_shift = pd.Series(ms, index=d.index)

    base_min = pd.to_numeric(d.get("PlotMin"), errors="coerce")
    base_max = pd.to_numeric(d.get("PlotMax"), errors="coerce")

    d["ReportLength"] = report_length.round().astype("Int64")
    d["ReportLengthBasis"] = report_length_basis_label(basis)
    d["ReportLengthShift"] = applied_shift
    d["ReportPlotMin"] = base_min + applied_shift
    d["ReportPlotMax"] = base_max + applied_shift
    return d


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------


_CLASS_PRIORITY = {
    "A_only": 1, "C_containing": 2, "U_containing": 3, "mixed_CU": 4,
    "length_only": 98, "unassigned": 99,
}


def _class_priority(cls):
    return _CLASS_PRIORITY.get(str(cls), 50)


def match_to_library(filtered: pd.DataFrame, lib: pd.DataFrame,
                     tolerance_value: float = 10, tolerance_unit: str = "ppm",
                     fill_unassigned: bool = True, fill_tolerance_da: float = 35,
                     length_only_tolerance_da: float = 165,
                     legacy_nearest_ladder: bool = True) -> pd.DataFrame:
    if len(filtered) == 0 or len(lib) == 0:
        return pd.DataFrame()

    effective_len_tol = max(float(length_only_tolerance_da) if np.isfinite(length_only_tolerance_da) else 165, 165)

    lib_arr_mass = lib["TheoreticalMass"].to_numpy()
    out_rows = []

    for _, obs in filtered.iterrows():
        obs_mass = float(obs["Mass"])
        err_da = obs_mass - lib_arr_mass
        abs_da = np.abs(err_da)
        err_ppm = 1e6 * err_da / lib_arr_mass
        abs_ppm = np.abs(err_ppm)
        in_tol = abs_ppm <= tolerance_value if tolerance_unit == "ppm" else abs_da <= tolerance_value
        cand_idx = np.where(in_tol)[0]

        status = None
        if len(cand_idx) == 0 and fill_unassigned:
            nearest = int(np.argmin(abs_da))
            if np.isfinite(abs_da[nearest]) and abs_da[nearest] <= fill_tolerance_da:
                cand_idx = np.array([nearest])
                status = "nearest_fill"
            elif np.isfinite(abs_da[nearest]) and abs_da[nearest] <= effective_len_tol:
                status = "length_only"
            else:
                status = "unassigned"
        elif len(cand_idx) == 0:
            nearest = int(np.argmin(abs_da))
            if np.isfinite(abs_da[nearest]) and abs_da[nearest] <= effective_len_tol:
                status = "length_only"
            else:
                status = "unassigned"
        else:
            status = "unique" if len(cand_idx) == 1 else "ambiguous"

        if status in ("length_only", "unassigned"):
            nearest = int(np.argmin(abs_da))
            n = lib.iloc[nearest]
            err_d = obs_mass - n["TheoreticalMass"]
            err_p = 1e6 * err_d / n["TheoreticalMass"]

            if legacy_nearest_ladder:
                # Legacy/manual mode: place observed peak on the nearest ladder species
                # and flag the assignment as nearest_ladder.
                row = obs.to_dict()
                row.update({
                    "TailBlockID": n["TailBlockID"],
                    "SourceFragmentID": n["SourceFragmentID"],
                    "SpeciesID": n["SpeciesID"],
                    "DisplayLabel": n["DisplayLabel"],
                    "SpeciesClass": n["SpeciesClass"],
                    "SpeciesClassLabel": n["SpeciesClassLabel"],
                    "A_Count": n["A_Count"],
                    "AssignedLength": n["AssignedLength"],
                    "C_Count": n["C_Count"],
                    "U_Count": n["U_Count"],
                    "NonA_Count": n["NonA_Count"],
                    "TotalTailModelResidues": n["TotalTailModelResidues"],
                    "PrefixContextLength": n["PrefixContextLength"],
                    "SuffixContextLength": n["SuffixContextLength"],
                    "FixedContextLength": n["FixedContextLength"],
                    "GeneratedSequence": n["GeneratedSequence"],
                    "FragmentLength": n["FragmentLength"],
                    "Composition": n["Composition"],
                    "PositionInterpretation": (
                        "Legacy nearest-ladder assignment for report plotting; strict mass tolerance "
                        "not met, so review mass error before using as confirmed chemistry"
                    ),
                    "TerminalType": n["TerminalType"],
                    "Termini": n["Termini"],
                    "NucleotideChemistry": n["NucleotideChemistry"],
                    "TheoreticalMass": n["TheoreticalMass"],
                    "ExpectedALength": n["ExpectedALength"],
                    "ExpectedAOnlyMass": n["ExpectedAOnlyMass"],
                    "MassErrorDa": round(err_d, 6),
                    "MassErrorPpm": round(err_p, 3),
                    "AnnotationStatus": "nearest_ladder",
                    "CandidateCount": 1,
                    "CandidateList": f"{n['DisplayLabel']} [nearest ladder] ({round(err_d, 4)} Da)",
                    "NearestSpecies": n["DisplayLabel"],
                    "NearestTheoreticalMass": n["TheoreticalMass"],
                    "NearestErrorDa": round(err_d, 6),
                    "NearestErrorPpm": round(err_p, 3),
                    "PlotMin": n["PlotMin"],
                    "PlotMax": n["PlotMax"],
                })
                out_rows.append(row)
                continue

            length_class = "length_only" if status == "length_only" else "unassigned"
            length_label = "Nearest length only" if status == "length_only" else "Unassigned"
            row = obs.to_dict()
            row.update({
                "TailBlockID": n["TailBlockID"],
                "SourceFragmentID": n["SourceFragmentID"],
                "SpeciesID": np.nan,
                "DisplayLabel": n["DisplayLabel"],
                "SpeciesClass": length_class,
                "SpeciesClassLabel": length_label,
                "A_Count": n["A_Count"],
                "AssignedLength": n["AssignedLength"],
                "C_Count": np.nan,
                "U_Count": np.nan,
                "NonA_Count": np.nan,
                "TotalTailModelResidues": np.nan,
                "PrefixContextLength": n["PrefixContextLength"],
                "SuffixContextLength": n["SuffixContextLength"],
                "FixedContextLength": n["FixedContextLength"],
                "GeneratedSequence": np.nan,
                "FragmentLength": n["FragmentLength"],
                "Composition": np.nan,
                "PositionInterpretation": (
                    "Nearest length only for report plotting; composition not assigned"
                    if status == "length_only" else
                    "Outside length-only tolerance; composition not assigned"
                ),
                "TerminalType": n["TerminalType"],
                "Termini": n["Termini"],
                "NucleotideChemistry": n["NucleotideChemistry"],
                "TheoreticalMass": np.nan,
                "ExpectedALength": n["ExpectedALength"],
                "ExpectedAOnlyMass": n["ExpectedAOnlyMass"],
                "MassErrorDa": np.nan,
                "MassErrorPpm": np.nan,
                "AnnotationStatus": status,
                "CandidateCount": 0,
                "CandidateList": np.nan,
                "NearestSpecies": n["DisplayLabel"],
                "NearestTheoreticalMass": n["TheoreticalMass"],
                "NearestErrorDa": round(err_d, 6),
                "NearestErrorPpm": round(err_p, 3),
                "PlotMin": n["PlotMin"],
                "PlotMax": n["PlotMax"],
            })
            out_rows.append(row)
            continue

        # Assigned: unique or ambiguous - pick best by abs Da, then class priority, then NonA
        cand = lib.iloc[cand_idx].copy()
        cand["MassErrorDa"] = obs_mass - cand["TheoreticalMass"]
        cand["AbsMassErrorDa"] = cand["MassErrorDa"].abs()
        cand["MassErrorPpm"] = 1e6 * cand["MassErrorDa"] / cand["TheoreticalMass"]
        cand["ClassPriority"] = cand["SpeciesClass"].map(_class_priority)
        cand = cand.sort_values(["AbsMassErrorDa", "ClassPriority", "NonA_Count"])
        best = cand.iloc[0]

        cand_strs = [f"{r['DisplayLabel']} [{r['TailBlockID']}] ({round(r['MassErrorDa'], 4)} Da)"
                     for _, r in cand.head(5).iterrows()]
        candidate_list = "; ".join(cand_strs)

        row = obs.to_dict()
        row.update({
            "TailBlockID": best["TailBlockID"],
            "SourceFragmentID": best["SourceFragmentID"],
            "SpeciesID": best["SpeciesID"],
            "DisplayLabel": best["DisplayLabel"],
            "SpeciesClass": best["SpeciesClass"],
            "SpeciesClassLabel": best["SpeciesClassLabel"],
            "A_Count": best["A_Count"],
            "AssignedLength": best["AssignedLength"],
            "C_Count": best["C_Count"],
            "U_Count": best["U_Count"],
            "NonA_Count": best["NonA_Count"],
            "TotalTailModelResidues": best["TotalTailModelResidues"],
            "PrefixContextLength": best["PrefixContextLength"],
            "SuffixContextLength": best["SuffixContextLength"],
            "FixedContextLength": best["FixedContextLength"],
            "GeneratedSequence": best["GeneratedSequence"],
            "FragmentLength": best["FragmentLength"],
            "Composition": best["Composition"],
            "PositionInterpretation": best["PositionInterpretation"],
            "TerminalType": best["TerminalType"],
            "Termini": best["Termini"],
            "NucleotideChemistry": best["NucleotideChemistry"],
            "TheoreticalMass": best["TheoreticalMass"],
            "ExpectedALength": best["ExpectedALength"],
            "ExpectedAOnlyMass": best["ExpectedAOnlyMass"],
            "MassErrorDa": round(best["MassErrorDa"], 6),
            "MassErrorPpm": round(best["MassErrorPpm"], 3),
            "AnnotationStatus": status,
            "CandidateCount": len(cand),
            "CandidateList": candidate_list,
            "NearestSpecies": best["DisplayLabel"],
            "NearestTheoreticalMass": best["TheoreticalMass"],
            "NearestErrorDa": round(best["MassErrorDa"], 6),
            "NearestErrorPpm": round(best["MassErrorPpm"], 3),
            "PlotMin": best["PlotMin"],
            "PlotMax": best["PlotMax"],
        })
        out_rows.append(row)

    return pd.DataFrame(out_rows)


# ---------------------------------------------------------------------------
# Summary and AWDI
# ---------------------------------------------------------------------------


def calculate_awdi_from_annotated(annotated: pd.DataFrame, basis: str = "A_only") -> float:
    if annotated is None or len(annotated) == 0:
        return float("nan")
    d = annotated[(annotated["AnnotationStatus"] != "unassigned")
                  & annotated["ExpectedAOnlyMass"].apply(lambda x: pd.notna(x) and np.isfinite(x))].copy()
    if basis == "A_only":
        d = d[d["SpeciesClass"] == "A_only"]
    elif basis == "report_default":
        d = d[d["SpeciesClass"].isin(["A_only", "C_containing"])]
    total_w = d["Fractional.Abundance"].sum()
    if len(d) == 0 or total_w <= 0:
        return float("nan")
    w = d["Fractional.Abundance"] / total_w
    return float(((d["Mass"] - d["ExpectedAOnlyMass"]) / d["ExpectedAOnlyMass"]).pow(2).mul(w).sum())


def class_pct(d: pd.DataFrame, cls: str) -> float:
    assigned = d[d["AnnotationStatus"] != "unassigned"]
    total = assigned["Fractional.Abundance"].sum()
    if total <= 0:
        return float("nan")
    return 100 * assigned.loc[assigned["SpeciesClass"] == cls, "Fractional.Abundance"].sum() / total


def unassigned_pct(d: pd.DataFrame) -> float:
    total = d["Fractional.Abundance"].sum()
    if total <= 0:
        return float("nan")
    return 100 * d.loc[d["AnnotationStatus"] == "unassigned", "Fractional.Abundance"].sum() / total


def build_summary_table(results: dict, awdi_basis: str) -> pd.DataFrame:
    rows = []
    for name, res in results.items():
        d = res.get("annotated")
        if d is None or len(d) == 0:
            rows.append({"File": name, "FilteredSpecies": 0, "AssignedSpecies": 0})
            continue
        assigned = d[d["AnnotationStatus"] != "unassigned"]
        most = assigned.iloc[assigned["Fractional.Abundance"].argmax()] if len(assigned) else None
        report_lengths = assigned["ReportLength"] if "ReportLength" in assigned else assigned.get("AssignedLength")
        rows.append({
            "File": name,
            "AWDI": round(res.get("awdi", float("nan")), 6) if res.get("awdi") == res.get("awdi") else None,
            "AWDIBasis": awdi_basis,
            "FilteredSpecies": len(d),
            "AssignedSpecies": len(assigned),
            "UniqueAssigned": int((d["AnnotationStatus"] == "unique").sum()),
            "Ambiguous": int((d["AnnotationStatus"] == "ambiguous").sum()),
            "NearestFilled": int((d["AnnotationStatus"] == "nearest_fill").sum()),
            "NearestLadder": int((d["AnnotationStatus"] == "nearest_ladder").sum()),
            "LengthOnly": int((d["AnnotationStatus"] == "length_only").sum()),
            "Unassigned": int((d["AnnotationStatus"] == "unassigned").sum()),
            "Unassigned % of total filtered abundance": round(unassigned_pct(d), 2),
            "Percent Poly(A)": round(class_pct(d, "A_only"), 2),
            "Percent Poly(A) with C": round(class_pct(d, "C_containing"), 2),
            "Percent Poly(A) with U": round(class_pct(d, "U_containing"), 2),
            "Percent Poly(A) with C/U": round(class_pct(d, "mixed_CU"), 2),
            "MostAbundantSpecies": most["DisplayLabel"] if most is not None else None,
            "MostAbundantClass": most["SpeciesClassLabel"] if most is not None else None,
            "MostAbundantA_Count": most["A_Count"] if most is not None else None,
            "MostAbundantReportLength": (most["ReportLength"] if most is not None and "ReportLength" in most else
                                         (most["AssignedLength"] if most is not None else None)),
            "ReportLengthBasis": most["ReportLengthBasis"] if most is not None and "ReportLengthBasis" in most else None,
            "MedianMass": round(float(d["Mass"].median()), 4) if len(d) else None,
            "MedianA_Count": round(float(assigned["A_Count"].median()), 1) if len(assigned) else None,
            "MedianReportLength": round(float(report_lengths.median()), 1) if len(assigned) and report_lengths is not None else None,
            "A_CountRange": f"{int(assigned['A_Count'].min())}\u2013{int(assigned['A_Count'].max())}" if len(assigned) else None,
            "ReportLengthRange": (f"{int(report_lengths.min())}\u2013{int(report_lengths.max())}"
                                  if len(assigned) and report_lengths is not None else None),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Mode-mass QC
# ---------------------------------------------------------------------------


def build_mode_mass_qc_table(results: dict, enabled: bool = True,
                             tolerance_value: float = 2, tolerance_unit: str = "Da") -> pd.DataFrame:
    unit = tolerance_unit if tolerance_unit in ("Da", "ppm") else "Da"
    tol = pd.to_numeric(tolerance_value, errors="coerce")
    if not pd.notna(tol) or tol < 0:
        tol = 20 if unit == "ppm" else 1

    rows = []
    for name, res in results.items():
        d = res.get("annotated")
        if d is None or len(d) == 0:
            rows.append({
                "File": name, "ModeQCEnabled": bool(enabled),
                "ModeQCPass": not bool(enabled),
                "ModeObservedMass": None, "ModeApexRT": None, "ModeRelativeAbundance": None,
                "ModeNearestSpecies": None, "ModeSpeciesClass": None,
                "ModeReportLength": None, "ModeAnnotationStatus": None,
                "ModeTheoreticalMass": None,
                "ModeMassErrorDa": None, "ModeAbsMassErrorDa": None,
                "ModeMassErrorPpm": None, "ModeAbsMassErrorPpm": None,
                "ModeQCTolerance": tol, "ModeQCToleranceUnit": unit,
                "ModeQCMessage": ("FAIL: no annotated filtered data available." if enabled else "QC disabled."),
            })
            continue

        mode = d.iloc[d["Fractional.Abundance"].argmax()]
        theoretical = pd.to_numeric(mode.get("TheoreticalMass"), errors="coerce")
        if not pd.notna(theoretical):
            theoretical = pd.to_numeric(mode.get("NearestTheoreticalMass"), errors="coerce")

        obs_mass = pd.to_numeric(mode.get("Mass"), errors="coerce")
        if pd.notna(theoretical) and pd.notna(obs_mass):
            err_da = float(obs_mass - theoretical)
            err_ppm = 1e6 * err_da / float(theoretical)
        else:
            err_da = np.nan
            err_ppm = np.nan
        abs_da = abs(err_da) if np.isfinite(err_da) else np.nan
        abs_ppm = abs(err_ppm) if np.isfinite(err_ppm) else np.nan

        pass_ = True
        msg = "QC disabled."
        if enabled:
            if not (np.isfinite(obs_mass) if pd.notna(obs_mass) else False) or not (
                    np.isfinite(theoretical) if pd.notna(theoretical) else False):
                pass_ = False
                msg = "FAIL: mode peak could not be compared with the generated ladder. Check sequence/range/library settings."
            elif unit == "ppm":
                pass_ = np.isfinite(abs_ppm) and abs_ppm <= tol
                msg = (f"PASS: mode peak matches generated ladder within {tol} ppm." if pass_ else
                       f"FAIL: mode peak differs from nearest generated ladder by {round(abs_ppm, 2)} ppm / "
                       f"{round(abs_da, 3)} Da. Check U chemistry, sequence context, termini, or range.")
            else:
                pass_ = np.isfinite(abs_da) and abs_da <= tol
                msg = (f"PASS: mode peak matches generated ladder within {tol} Da." if pass_ else
                       f"FAIL: mode peak differs from nearest generated ladder by {round(abs_da, 3)} Da / "
                       f"{round(abs_ppm, 2)} ppm. Check U chemistry, sequence context, termini, or range.")

        rows.append({
            "File": name, "ModeQCEnabled": bool(enabled), "ModeQCPass": bool(pass_),
            "ModeObservedMass": round(float(obs_mass), 6) if pd.notna(obs_mass) else None,
            "ModeApexRT": float(mode.get("Apex.RT")) if pd.notna(mode.get("Apex.RT")) else None,
            "ModeRelativeAbundance": float(mode.get("Fractional.Abundance")) if pd.notna(mode.get("Fractional.Abundance")) else None,
            "ModeNearestSpecies": str(mode.get("DisplayLabel")) if pd.notna(mode.get("DisplayLabel")) else None,
            "ModeSpeciesClass": str(mode.get("SpeciesClassLabel")) if pd.notna(mode.get("SpeciesClassLabel")) else None,
            "ModeReportLength": (float(mode.get("ReportLength")) if "ReportLength" in mode and pd.notna(mode.get("ReportLength"))
                                 else (float(mode.get("AssignedLength")) if pd.notna(mode.get("AssignedLength")) else None)),
            "ModeAnnotationStatus": str(mode.get("AnnotationStatus")) if pd.notna(mode.get("AnnotationStatus")) else None,
            "ModeTheoreticalMass": round(float(theoretical), 6) if pd.notna(theoretical) else None,
            "ModeMassErrorDa": round(err_da, 6) if np.isfinite(err_da) else None,
            "ModeAbsMassErrorDa": round(abs_da, 6) if np.isfinite(abs_da) else None,
            "ModeMassErrorPpm": round(err_ppm, 3) if np.isfinite(err_ppm) else None,
            "ModeAbsMassErrorPpm": round(abs_ppm, 3) if np.isfinite(abs_ppm) else None,
            "ModeQCTolerance": tol, "ModeQCToleranceUnit": unit,
            "ModeQCMessage": msg,
        })
    return pd.DataFrame(rows)


def mode_qc_all_pass(qc_table: pd.DataFrame) -> bool:
    if qc_table is None or len(qc_table) == 0 or "ModeQCPass" not in qc_table.columns:
        return False
    return bool(qc_table["ModeQCPass"].all())


def mode_qc_failure_message(qc_table: pd.DataFrame) -> str:
    if qc_table is None or len(qc_table) == 0:
        return "Mode-mass QC failed: no QC table available."
    failed = qc_table[qc_table["ModeQCPass"] != True]  # noqa: E712
    if len(failed) == 0:
        return "Mode-mass QC passed."
    return "\n".join(sorted(set(failed["ModeQCMessage"].astype(str))))


# ---------------------------------------------------------------------------
# Plot-data builders
# ---------------------------------------------------------------------------


def build_class_plot_data(annotated: pd.DataFrame, include_unassigned: bool = False,
                          duplicate_rule: str = "dominant") -> pd.DataFrame:
    status_keep = {"unique", "ambiguous", "nearest_fill", "nearest_ladder"}
    length_col = "ReportLength" if "ReportLength" in annotated.columns else "AssignedLength"

    d = annotated if include_unassigned else annotated[
        annotated["AnnotationStatus"].isin(status_keep) & annotated[length_col].notna()]
    if len(d) == 0:
        return pd.DataFrame()

    d = d.copy()
    d["PlotLength"] = pd.to_numeric(d[length_col], errors="coerce")
    d["PlotMinForReport"] = pd.to_numeric(d.get("ReportPlotMin", d.get("PlotMin")), errors="coerce")
    d["PlotMaxForReport"] = pd.to_numeric(d.get("ReportPlotMax", d.get("PlotMax")), errors="coerce")

    d["SpeciesLabel"] = d["SpeciesClass"].map(lambda x: class_to_label(x))
    if "SampleLabel" not in d.columns:
        d["SampleLabel"] = d["File"].map(clean_sample_label)

    d = d[d["PlotLength"].apply(lambda x: pd.notna(x) and np.isfinite(x))]
    d = d[d["Fractional.Abundance"].apply(lambda x: pd.notna(x) and np.isfinite(x))]
    if len(d) == 0:
        return pd.DataFrame()

    agg_fn = (lambda x: x.sum()) if duplicate_rule == "sum" else (lambda x: x.max())
    grouped = d.groupby(["File", "SampleLabel", "TailBlockID", "PlotLength", "SpeciesClass", "SpeciesLabel"],
                        dropna=False, as_index=False).agg(
        RelAbundance=("Fractional.Abundance", agg_fn),
        PlotMin=("PlotMinForReport", "first"),
        PlotMax=("PlotMaxForReport", "first"),
    )
    return grouped


def build_total_plot_data_from_class(class_data: pd.DataFrame, threshold_pct: float = 10) -> pd.DataFrame:
    if class_data is None or len(class_data) == 0:
        return pd.DataFrame()
    total = class_data.groupby(["File", "SampleLabel", "PlotLength"], as_index=False).agg(
        TotalRaw=("RelAbundance", "sum"),
        PlotMin=("PlotMin", "first"),
        PlotMax=("PlotMax", "first"),
    )
    max_data = total.groupby(["File", "SampleLabel"], as_index=False).agg(MaxRaw=("TotalRaw", "max"))
    total = total.merge(max_data, on=["File", "SampleLabel"], how="left")
    total["TotalNorm"] = np.where((total["MaxRaw"] > 0) & total["MaxRaw"].notna(),
                                   100 * total["TotalRaw"] / total["MaxRaw"], 0)
    total["ThresholdClass"] = np.where(total["TotalNorm"] >= threshold_pct, "above_10", "below_10")
    return total
