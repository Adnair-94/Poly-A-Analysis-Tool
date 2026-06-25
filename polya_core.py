"""
Poly(A) Heterogeneity Analysis Tool - core science module.

Version 2.1 adds:
  - Explicit RNA, DNA coding-strand, and DNA template-strand input handling.
  - Independent 5' and 3' terminal chemistry for generated tail fragments.
  - C/U placement at the 5' side, 3' side, or both sides of the A-run.
  - A diagnostic for the characteristic ~79.966 Da phosphate-state offset.

RNase T1 defaults are retained:
  - Internal T1 fragments: 5'-OH / 3'-phosphate.
  - Transcript-terminal T1 fragment: 5'-OH / 3'-OH.

Mass conventions use monoisotopic RNA residue masses. The neutral 5'-OH /
3'-OH oligonucleotide mass is sum(residues) + H2O - HPO3. Each terminal
phosphate adds HPO3. A 2',3'-cyclic phosphate adds HPO3 - H2O relative to
3'-OH.
"""

from __future__ import annotations

import os
import re
import unicodedata
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

APP_VERSION = "2.1.0"

# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------

CLASS_COLOURS = {
    "A_only": "#008B8B",
    "C_containing": "#A64CE6",
    "U_containing": "#F97316",
    "mixed_CU": "#64748B",
    "ambiguous": "#9CA3AF",
    "length_only": "#D9D9D9",
    "unassigned": "#D9D9D9",
    "above_10": "#008B8B",
    "below_10": "#D9D9D9",
}

CLASS_LABELS = {
    "A_only": "Poly(A)",
    "C_containing": "Poly(A) with C",
    "U_containing": "Poly(A) with U",
    "mixed_CU": "Poly(A) with C/U",
    "ambiguous": "Ambiguous",
    "length_only": "Nearest length only",
    "unassigned": "Unassigned",
}

INPUT_SEQUENCE_LABELS = {
    "rna": "RNA sequence (5'→3')",
    "dna_coding": "DNA coding/non-template strand (5'→3')",
    "dna_template": "DNA template/antisense strand (5'→3')",
}

FIVE_PRIME_LABELS = {
    "5OH": "5'-OH",
    "5p": "5'-phosphate",
}

THREE_PRIME_LABELS = {
    "3OH": "3'-OH",
    "3p": "3'-phosphate",
    "3cyclicp": "2',3'-cyclic phosphate",
}

VARIANT_PLACEMENT_LABELS = {
    "three_prime_of_a": "3' of the A-run",
    "five_prime_of_a": "5' of the A-run",
    "both": "Both 5' and 3' possibilities",
    "none": "A-only",
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
N1M_PSEUDOURIDINE_SHIFT = 14.0156501

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
    """Match the former R workflow: lower-case and replace punctuation runs."""
    out = []
    for c in cols:
        c2 = str(c).lower()
        c2 = re.sub(r"[\W_]+", ".", c2, flags=re.UNICODE)
        out.append(c2)
    return out


def _first_matching_col(nms, patterns):
    for pattern in patterns:
        rx = re.compile(pattern)
        for name in nms:
            if rx.search(name):
                return name
    return None


def clean_sample_label(path: str) -> str:
    base = os.path.splitext(os.path.basename(path))[0]
    base = re.sub(r"[-_.]+", " ", base)
    return re.sub(r"\s+", " ", base).strip()


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
# BioPharma Finder parsing and filtering
# ---------------------------------------------------------------------------


def preprocess_uploaded_excel(filepath: str) -> pd.DataFrame:
    df = pd.read_excel(filepath, sheet_name=0)
    df.columns = _normalise_names(df.columns)
    mass_col = _first_matching_col(
        df.columns,
        [r"^monoisotopic\.mass$", r"monoisotopic.*mass", r"^mass$"],
    )
    rt_col = _first_matching_col(
        df.columns,
        [r"^apex\.rt$", r"apex.*rt", r"retention.*time", r"^rt$"],
    )
    rel_col = _first_matching_col(
        df.columns,
        [r"^relative\.abundance$", r"relative.*abundance", r"fractional.*abundance", r"rel.*abund"],
    )
    if mass_col is None or rt_col is None or rel_col is None:
        raise ValueError(
            "Excel file must contain Monoisotopic Mass, Relative Abundance, and Apex RT columns."
        )
    out = pd.DataFrame(
        {
            "Mass": pd.to_numeric(df[mass_col], errors="coerce"),
            "Apex.RT": pd.to_numeric(df[rt_col], errors="coerce"),
            "Fractional.Abundance": pd.to_numeric(df[rel_col], errors="coerce"),
        }
    )
    return out.dropna().reset_index(drop=True)


def filter_peak_data(
    data: pd.DataFrame,
    rt_range: Sequence[float] = (0, 100),
    abundance_range: Sequence[float] = (0, 100),
    mass_range: Sequence[float] = (0, np.inf),
) -> pd.DataFrame:
    if data is None or len(data) == 0:
        return pd.DataFrame(columns=["Mass", "Apex.RT", "Fractional.Abundance"])
    d = data.copy()
    mask = (
        d["Apex.RT"].between(float(rt_range[0]), float(rt_range[1]), inclusive="both")
        & d["Fractional.Abundance"].between(
            float(abundance_range[0]), float(abundance_range[1]), inclusive="both"
        )
        & d["Mass"].between(float(mass_range[0]), float(mass_range[1]), inclusive="both")
    )
    return d.loc[mask].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Sequence conversion
# ---------------------------------------------------------------------------


def _clean_base_string(seq: str) -> str:
    if seq is None or not str(seq).strip():
        raise ValueError("Sequence input is empty.")
    s = str(seq).upper()
    s = re.sub(r"[^ACGTU]", "", s)
    if not s:
        raise ValueError("Sequence contains no valid A/C/G/T/U bases after cleanup.")
    return s


def reverse_complement_dna(seq: str) -> str:
    s = _clean_base_string(seq).replace("U", "T")
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def convert_sequence_to_rna(seq: str, input_type: str = "rna") -> str:
    """
    Convert a user sequence into the RNA sequence that is digested in silico.

    input_type:
      rna          RNA supplied 5'→3'. T is accepted and converted to U.
      dna_coding   Coding/non-template DNA supplied 5'→3'. T→U.
      dna_template Template/antisense DNA supplied 5'→3'. Reverse complement,
                   then T→U.
    """
    kind = str(input_type or "rna")
    aliases = {
        "RNA": "rna",
        "dna": "dna_coding",
        "DNA": "dna_coding",
        "coding": "dna_coding",
        "template": "dna_template",
        "antisense": "dna_template",
    }
    kind = aliases.get(kind, kind)
    if kind not in INPUT_SEQUENCE_LABELS:
        raise ValueError("input_type must be 'rna', 'dna_coding', or 'dna_template'.")

    clean = _clean_base_string(seq)
    if kind == "dna_template":
        return reverse_complement_dna(clean).replace("T", "U")
    return clean.replace("T", "U")


def normalise_rna_sequence(seq: str) -> str:
    """Backward-compatible RNA cleanup used by existing callers."""
    return convert_sequence_to_rna(seq, "rna")


def rna_to_dna(seq: str) -> str:
    return normalise_rna_sequence(seq).replace("U", "T")


def sequence_conversion_table(seq: str, input_type: str) -> pd.DataFrame:
    cleaned = _clean_base_string(seq)
    converted = convert_sequence_to_rna(seq, input_type)
    return pd.DataFrame(
        [
            {
                "InputType": INPUT_SEQUENCE_LABELS.get(input_type, input_type),
                "CleanedInputSequence": cleaned,
                "ConvertedRNASequence": converted,
                "ConvertedLength": len(converted),
            }
        ]
    )


# ---------------------------------------------------------------------------
# Sequence composition and mass logic
# ---------------------------------------------------------------------------


def composition_string(seq: str) -> str:
    s = normalise_rna_sequence(seq)
    parts = []
    for base in "ACGU":
        count = s.count(base)
        if count > 0:
            parts.append(f"{base}{count}")
    return "".join(parts)


def mass_table_for_chemistry(chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> dict:
    masses = dict(RNA_RESIDUE_MASS)
    if chemistry == "N1mPseudoU":
        masses["U"] += N1M_PSEUDOURIDINE_SHIFT
    if u_mass_shift and np.isfinite(u_mass_shift):
        masses["U"] += float(u_mass_shift)
    return masses


def normalise_five_prime_terminus(value: str) -> str:
    aliases = {
        "5'-OH": "5OH",
        "5oh": "5OH",
        "OH": "5OH",
        "5' phosphate": "5p",
        "5'-phosphate": "5p",
        "5phosphate": "5p",
        "phosphate": "5p",
    }
    v = aliases.get(str(value), str(value))
    if v not in FIVE_PRIME_LABELS:
        raise ValueError("5' terminus must be '5OH' or '5p'.")
    return v


def normalise_three_prime_terminus(value: str) -> str:
    aliases = {
        "3'-OH": "3OH",
        "3oh": "3OH",
        "OH": "3OH",
        "3' phosphate": "3p",
        "3'-phosphate": "3p",
        "3phosphate": "3p",
        "phosphate": "3p",
        "2',3'-cyclic phosphate": "3cyclicp",
        "3' cyclic phosphate": "3cyclicp",
        "cyclic": "3cyclicp",
    }
    v = aliases.get(str(value), str(value))
    if v not in THREE_PRIME_LABELS:
        raise ValueError("3' terminus must be '3OH', '3p', or '3cyclicp'.")
    return v


def termini_label(five_prime: str, three_prime: str) -> str:
    f = normalise_five_prime_terminus(five_prime)
    t = normalise_three_prime_terminus(three_prime)
    return f"{FIVE_PRIME_LABELS[f]} / {THREE_PRIME_LABELS[t]}"


def terminal_mass_adjustment(five_prime: str = "5OH", three_prime: str = "3OH") -> float:
    """Mass adjustment relative to sum of RNA residue masses."""
    f = normalise_five_prime_terminus(five_prime)
    t = normalise_three_prime_terminus(three_prime)
    adjustment = H2O_MASS - HPO3_MASS  # 5'-OH / 3'-OH
    if f == "5p":
        adjustment += HPO3_MASS
    if t == "3p":
        adjustment += HPO3_MASS
    elif t == "3cyclicp":
        adjustment += HPO3_MASS - H2O_MASS
    return adjustment


def calculate_rna_mass(
    seq: str,
    terminal_type: str = "3OH",
    chemistry: str = "canonical_U",
    u_mass_shift: float = 0.0,
    five_prime: str = "5OH",
    three_prime: Optional[str] = None,
) -> float:
    """
    Calculate neutral monoisotopic RNA mass.

    terminal_type is retained for backward compatibility and is interpreted as
    the 3' terminus unless three_prime is supplied explicitly.
    """
    s = normalise_rna_sequence(seq)
    masses = mass_table_for_chemistry(chemistry, u_mass_shift)
    three = terminal_type if three_prime is None else three_prime
    return sum(masses[b] for b in s) + terminal_mass_adjustment(five_prime, three)


def find_a_runs(seq: str) -> pd.DataFrame:
    s = normalise_rna_sequence(seq)
    rows = []
    for match in re.finditer(r"A+", s):
        rows.append(
            {
                "RunStart": match.start() + 1,
                "RunEnd": match.end(),
                "RunLength": match.end() - match.start(),
            }
        )
    return pd.DataFrame(rows, columns=["RunStart", "RunEnd", "RunLength"])


# ---------------------------------------------------------------------------
# RNase T1 digestion
# ---------------------------------------------------------------------------


def _digest_row(
    fragment_number: int,
    fragment_id: str,
    start: int,
    end: int,
    fragment_sequence: str,
    five_prime: str,
    three_prime: str,
    chemistry: str,
    u_mass_shift: float,
) -> dict:
    runs = find_a_runs(fragment_sequence)
    max_run = int(runs["RunLength"].max()) if len(runs) else 0
    mass = calculate_rna_mass(
        fragment_sequence,
        chemistry=chemistry,
        u_mass_shift=u_mass_shift,
        five_prime=five_prime,
        three_prime=three_prime,
    )
    return {
        "FragmentNumber": fragment_number,
        "FragmentID": fragment_id,
        "Start": start,
        "End": end,
        "FragmentSequence": fragment_sequence,
        "FragmentLength": len(fragment_sequence),
        "AutoFivePrimeTerminus": five_prime,
        "AutoThreePrimeTerminus": three_prime,
        "FivePrimeTerminus": five_prime,
        "ThreePrimeTerminus": three_prime,
        "TerminalType": three_prime,  # backward compatibility
        "Termini": termini_label(five_prime, three_prime),
        "TerminalOverrideApplied": False,
        "Composition": composition_string(fragment_sequence),
        "MaxContinuousA": max_run,
        "TheoreticalMass": round(mass, 6),
    }


def digest_rnase_t1(seq: str, chemistry: str = "canonical_U", u_mass_shift: float = 0.0) -> pd.DataFrame:
    s = normalise_rna_sequence(seq)
    n = len(s)
    g_positions = [i + 1 for i, ch in enumerate(s) if ch == "G"]
    cut_ends = [position for position in g_positions if position < n]
    fragment_ends = list(dict.fromkeys(cut_ends + [n])) or [n]

    rows = []
    start = 1
    for index, end in enumerate(fragment_ends, start=1):
        if end < start:
            continue
        fragment_sequence = s[start - 1 : end]
        three_prime = "3p" if end < n and s[end - 1] == "G" else "3OH"
        rows.append(
            _digest_row(
                index,
                f"T1_F{index}",
                start,
                end,
                fragment_sequence,
                "5OH",
                three_prime,
                chemistry,
                u_mass_shift,
            )
        )
        start = end + 1
    return pd.DataFrame(rows)


def digest_tail_fragment_input(
    seq: str, chemistry: str = "canonical_U", u_mass_shift: float = 0.0
) -> pd.DataFrame:
    s = normalise_rna_sequence(seq)
    return pd.DataFrame(
        [
            _digest_row(
                1,
                "T1_tail_input",
                1,
                len(s),
                s,
                "5OH",
                "3OH",
                chemistry,
                u_mass_shift,
            )
        ]
    )


def detect_tail_blocks(digest_df: pd.DataFrame, min_a_run: int = 10) -> pd.DataFrame:
    rows = []
    for _, fragment in digest_df.iterrows():
        fragment_sequence = fragment["FragmentSequence"]
        runs = find_a_runs(fragment_sequence)
        if len(runs) == 0:
            continue
        runs = runs[runs["RunLength"] >= int(min_a_run)]
        for block_index, (_, run) in enumerate(runs.iterrows(), start=1):
            run_start = int(run["RunStart"])
            run_end = int(run["RunEnd"])
            prefix = fragment_sequence[: run_start - 1]
            suffix = fragment_sequence[run_end:]
            rows.append(
                {
                    "TailBlockNumber": len(rows) + 1,
                    "TailBlockID": f"{fragment['FragmentID']}_Ablock{block_index}",
                    "SourceFragmentID": fragment["FragmentID"],
                    "SourceFragmentNumber": fragment["FragmentNumber"],
                    "SourceFragmentSequence": fragment_sequence,
                    "SourceFragmentLength": fragment["FragmentLength"],
                    "AutoFivePrimeTerminus": fragment.get("AutoFivePrimeTerminus", "5OH"),
                    "AutoThreePrimeTerminus": fragment.get(
                        "AutoThreePrimeTerminus", fragment.get("TerminalType", "3OH")
                    ),
                    "FivePrimeTerminus": fragment.get("FivePrimeTerminus", "5OH"),
                    "ThreePrimeTerminus": fragment.get(
                        "ThreePrimeTerminus", fragment.get("TerminalType", "3OH")
                    ),
                    "TerminalType": fragment.get("ThreePrimeTerminus", fragment.get("TerminalType", "3OH")),
                    "Termini": fragment["Termini"],
                    "TerminalOverrideApplied": bool(fragment.get("TerminalOverrideApplied", False)),
                    "RunStartInFragment": run_start,
                    "RunEndInFragment": run_end,
                    "RunStartInInput": int(fragment["Start"]) + run_start - 1,
                    "RunEndInInput": int(fragment["Start"]) + run_end - 1,
                    "DetectedALength": int(run["RunLength"]),
                    "PrefixBeforeA": prefix,
                    "SuffixAfterA": suffix,
                }
            )
    return pd.DataFrame(rows)


def _resolve_terminus_choice(auto_value: str, selected_value: str, prime: str) -> str:
    selected = str(selected_value or "auto")
    if selected == "auto":
        return auto_value
    if prime == "5":
        return normalise_five_prime_terminus(selected)
    return normalise_three_prime_terminus(selected)


def apply_tail_terminal_overrides(
    digest_df: pd.DataFrame,
    tail_fragment_ids: Sequence[str],
    five_prime_choice: str = "auto",
    three_prime_choice: str = "auto",
    chemistry: str = "canonical_U",
    u_mass_shift: float = 0.0,
) -> pd.DataFrame:
    d = digest_df.copy()
    ids = set(str(x) for x in tail_fragment_ids)
    for idx, row in d.iterrows():
        if str(row["FragmentID"]) not in ids:
            continue
        auto_five = normalise_five_prime_terminus(row.get("AutoFivePrimeTerminus", "5OH"))
        auto_three = normalise_three_prime_terminus(
            row.get("AutoThreePrimeTerminus", row.get("TerminalType", "3OH"))
        )
        effective_five = _resolve_terminus_choice(auto_five, five_prime_choice, "5")
        effective_three = _resolve_terminus_choice(auto_three, three_prime_choice, "3")
        override_applied = effective_five != auto_five or effective_three != auto_three
        d.at[idx, "FivePrimeTerminus"] = effective_five
        d.at[idx, "ThreePrimeTerminus"] = effective_three
        d.at[idx, "TerminalType"] = effective_three
        d.at[idx, "Termini"] = termini_label(effective_five, effective_three)
        d.at[idx, "TerminalOverrideApplied"] = override_applied
        d.at[idx, "TheoreticalMass"] = round(
            calculate_rna_mass(
                row["FragmentSequence"],
                chemistry=chemistry,
                u_mass_shift=u_mass_shift,
                five_prime=effective_five,
                three_prime=effective_three,
            ),
            6,
        )
    return d


# ---------------------------------------------------------------------------
# Species-library generation
# ---------------------------------------------------------------------------


def variant_class(c: int, u: int) -> str:
    if c == 0 and u == 0:
        return "A_only"
    if c > 0 and u == 0:
        return "C_containing"
    if c == 0 and u > 0:
        return "U_containing"
    return "mixed_CU"


def species_display_label(a: int, c: int, u: int, placement: str = "none") -> str:
    if c == 0 and u == 0:
        return f"A{a}"
    parts = []
    if c > 0:
        parts.append("C" if c == 1 else f"C{c}")
    if u > 0:
        parts.append("U" if u == 1 else f"U{u}")
    side = "5′ of A" if placement == "five_prime_of_a" else "3′ of A"
    return f"A{a} with {'/'.join(parts)} ({side})"


def generate_variant_combos(max_c: int = 1, max_u: int = 1, allow_mixed: bool = True) -> pd.DataFrame:
    max_c = max(0, int(max_c))
    max_u = max(0, int(max_u))
    rows = []
    for c_count in range(max_c + 1):
        for u_count in range(max_u + 1):
            if not allow_mixed and c_count > 0 and u_count > 0:
                continue
            rows.append({"C_Count": c_count, "U_Count": u_count})
    df = pd.DataFrame(rows)
    if len(df):
        df["_order"] = df["C_Count"] + df["U_Count"]
        df = (
            df.sort_values(["_order", "C_Count", "U_Count"])
            .drop(columns="_order")
            .reset_index(drop=True)
        )
    return df


def _variant_placements(c: int, u: int, placement_setting: str) -> list[str]:
    if c + u == 0:
        return ["none"]
    setting = str(placement_setting or "three_prime_of_a")
    if setting == "both":
        return ["five_prime_of_a", "three_prime_of_a"]
    if setting not in ("five_prime_of_a", "three_prime_of_a"):
        raise ValueError("variant_placement must be 'five_prime_of_a', 'three_prime_of_a', or 'both'.")
    return [setting]


def _place_variant(prefix: str, a_count: int, c_count: int, u_count: int, suffix: str, placement: str) -> str:
    a_run = "A" * int(a_count)
    motif = ("C" * int(c_count)) + ("U" * int(u_count))
    if placement == "five_prime_of_a":
        return prefix + motif + a_run + suffix
    return prefix + a_run + motif + suffix


@dataclass
class GeneratedLibrary:
    OriginalInputSequence: str
    InputSequenceType: str
    InputSequence: str
    ConvertedRNASequence: str
    DigestTable: pd.DataFrame
    TailBlocks: pd.DataFrame
    SpeciesLibrary: pd.DataFrame


def generate_species_library(
    sample_id: str,
    seq: str,
    chemistry: str = "canonical_U",
    u_mass_shift: float = 0.0,
    sequence_mode: str = "digest_after_G",
    sequence_input_type: str = "rna",
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
    variant_placement: str = "three_prime_of_a",
    five_prime_terminus: str = "auto",
    three_prime_terminus: str = "auto",
) -> GeneratedLibrary:
    original_seq = str(seq)
    clean_rna = convert_sequence_to_rna(seq, sequence_input_type)

    if sequence_mode == "tail_fragment":
        digest_df = digest_tail_fragment_input(clean_rna, chemistry, u_mass_shift)
    elif sequence_mode == "digest_after_G":
        digest_df = digest_rnase_t1(clean_rna, chemistry, u_mass_shift)
    else:
        raise ValueError("sequence_mode must be 'digest_after_G' or 'tail_fragment'.")

    initial_blocks = detect_tail_blocks(digest_df, min_a_run=min_a_run)
    if len(initial_blocks) == 0:
        raise ValueError(
            "No tail-positive T1 fragment detected. Lower the minimum A-run threshold or check the sequence."
        )

    digest_df = apply_tail_terminal_overrides(
        digest_df,
        initial_blocks["SourceFragmentID"].tolist(),
        five_prime_choice=five_prime_terminus,
        three_prime_choice=three_prime_terminus,
        chemistry=chemistry,
        u_mass_shift=u_mass_shift,
    )
    tail_blocks = detect_tail_blocks(digest_df, min_a_run=min_a_run)

    sample_id = (sample_id or "").strip() or "Sample"
    combos = generate_variant_combos(max_c, max_u, allow_mixed)

    species_rows = []
    block_rows = []
    for _, tail_block in tail_blocks.iterrows():
        detected_n = int(tail_block["DetectedALength"])
        prefix = str(tail_block["PrefixBeforeA"])
        suffix = str(tail_block["SuffixAfterA"])
        prefix_len = len(prefix)
        suffix_len = len(suffix)
        fixed_context = prefix_len + suffix_len

        rmin = pd.to_numeric(range_min, errors="coerce")
        rmax = pd.to_numeric(range_max, errors="coerce")
        rmin_ok = pd.notna(rmin)
        rmax_ok = pd.notna(rmax)

        if use_auto_ranges or not rmin_ok or not rmax_ok or rmax < rmin:
            lo = max(0, detected_n - int(auto_range_left))
            hi = detected_n + int(auto_range_right)
            expected_n = detected_n
        else:
            basis = range_input_basis if range_input_basis in ("a_count", "report_length") else "report_length"
            range_shift = 0.0
            if basis == "report_length":
                report_basis = (
                    report_length_basis
                    if report_length_basis in ("a_count", "context", "full_fragment", "manual_shift")
                    else "context"
                )
                if report_basis in ("context", "full_fragment"):
                    range_shift = fixed_context
                elif report_basis == "manual_shift":
                    manual = pd.to_numeric(manual_report_shift, errors="coerce")
                    range_shift = float(manual) if pd.notna(manual) else 0.0
            lo = max(0, int(np.floor(float(rmin) - range_shift)))
            hi = max(lo, int(np.ceil(float(rmax) - range_shift)))
            expected_n = (
                int(expected_a_length)
                if expected_a_length is not None and np.isfinite(expected_a_length)
                else detected_n
            )

        effective_five = normalise_five_prime_terminus(tail_block["FivePrimeTerminus"])
        effective_three = normalise_three_prime_terminus(tail_block["ThreePrimeTerminus"])
        expected_seq = prefix + ("A" * expected_n) + suffix
        expected_mass = calculate_rna_mass(
            expected_seq,
            chemistry=chemistry,
            u_mass_shift=u_mass_shift,
            five_prime=effective_five,
            three_prime=effective_three,
        )

        block_rows.append(
            {
                "SampleID": sample_id,
                "TailBlockID": tail_block["TailBlockID"],
                "SourceFragmentID": tail_block["SourceFragmentID"],
                "SourceFragmentSequence": tail_block["SourceFragmentSequence"],
                "DetectedALength": detected_n,
                "ExpectedALength": expected_n,
                "ARangeMin": lo,
                "ARangeMax": hi,
                "InputRangeMin": rmin if rmin_ok else np.nan,
                "InputRangeMax": rmax if rmax_ok else np.nan,
                "RangeInputBasis": "A-count" if range_input_basis == "a_count" else "Report length",
                "ExpectedAOnlyMass": round(expected_mass, 6),
                "PrefixBeforeA": prefix,
                "SuffixAfterA": suffix,
                "PrefixContextLength": prefix_len,
                "SuffixContextLength": suffix_len,
                "FixedContextLength": fixed_context,
                "AutoFivePrimeTerminus": tail_block["AutoFivePrimeTerminus"],
                "AutoThreePrimeTerminus": tail_block["AutoThreePrimeTerminus"],
                "FivePrimeTerminus": effective_five,
                "ThreePrimeTerminus": effective_three,
                "TerminalType": effective_three,
                "Termini": termini_label(effective_five, effective_three),
                "TerminalOverrideApplied": bool(tail_block["TerminalOverrideApplied"]),
                "VariantPlacementSetting": variant_placement,
            }
        )

        for a_count in range(lo, hi + 1):
            for _, combo in combos.iterrows():
                c_count = int(combo["C_Count"])
                u_count = int(combo["U_Count"])
                for placement in _variant_placements(c_count, u_count, variant_placement):
                    generated_seq = _place_variant(
                        prefix, a_count, c_count, u_count, suffix, placement
                    )
                    mass = calculate_rna_mass(
                        generated_seq,
                        chemistry=chemistry,
                        u_mass_shift=u_mass_shift,
                        five_prime=effective_five,
                        three_prime=effective_three,
                    )
                    species_class = variant_class(c_count, u_count)
                    placement_code = {
                        "none": "A",
                        "five_prime_of_a": "P5",
                        "three_prime_of_a": "P3",
                    }[placement]
                    species_rows.append(
                        {
                            "SampleID": sample_id,
                            "TailBlockID": tail_block["TailBlockID"],
                            "SourceFragmentID": tail_block["SourceFragmentID"],
                            "SpeciesID": (
                                f"{sample_id}_{tail_block['TailBlockID']}_A{a_count}_"
                                f"C{c_count}_U{u_count}_{placement_code}_"
                                f"{effective_five}_{effective_three}"
                            ),
                            "DisplayLabel": species_display_label(
                                a_count, c_count, u_count, placement
                            ),
                            "SpeciesClass": species_class,
                            "SpeciesClassLabel": class_to_label(species_class),
                            "A_Count": int(a_count),
                            "AssignedLength": int(a_count),
                            "C_Count": int(c_count),
                            "U_Count": int(u_count),
                            "NonA_Count": int(c_count + u_count),
                            "TotalTailModelResidues": int(a_count + c_count + u_count),
                            "PrefixContextLength": prefix_len,
                            "SuffixContextLength": suffix_len,
                            "FixedContextLength": fixed_context,
                            "VariantPlacement": placement,
                            "VariantPlacementLabel": VARIANT_PLACEMENT_LABELS[placement],
                            "GeneratedSequence": generated_seq,
                            "FragmentLength": len(generated_seq),
                            "Composition": composition_string(generated_seq),
                            "PositionInterpretation": (
                                "A-only"
                                if c_count + u_count == 0
                                else (
                                    f"C/U modelled {VARIANT_PLACEMENT_LABELS[placement]}; "
                                    "intact mass alone does not establish positional isomerism"
                                )
                            ),
                            "AutoFivePrimeTerminus": tail_block["AutoFivePrimeTerminus"],
                            "AutoThreePrimeTerminus": tail_block["AutoThreePrimeTerminus"],
                            "FivePrimeTerminus": effective_five,
                            "ThreePrimeTerminus": effective_three,
                            "TerminalType": effective_three,
                            "Termini": termini_label(effective_five, effective_three),
                            "TerminalOverrideApplied": bool(tail_block["TerminalOverrideApplied"]),
                            "TerminalMassAdjustment": round(
                                terminal_mass_adjustment(effective_five, effective_three), 6
                            ),
                            "NucleotideChemistry": chemistry,
                            "TheoreticalMass": round(mass, 6),
                            "ExpectedALength": expected_n,
                            "ExpectedAOnlyMass": round(expected_mass, 6),
                            "PlotMin": lo,
                            "PlotMax": hi,
                            "InputRangeMin": rmin if rmin_ok else np.nan,
                            "InputRangeMax": rmax if rmax_ok else np.nan,
                            "RangeInputBasis": (
                                "A-count" if range_input_basis == "a_count" else "Report length"
                            ),
                        }
                    )

    return GeneratedLibrary(
        OriginalInputSequence=original_seq,
        InputSequenceType=sequence_input_type,
        InputSequence=clean_rna,
        ConvertedRNASequence=clean_rna,
        DigestTable=digest_df,
        TailBlocks=pd.DataFrame(block_rows),
        SpeciesLibrary=pd.DataFrame(species_rows),
    )


# ---------------------------------------------------------------------------
# Report-length helpers
# ---------------------------------------------------------------------------


def add_report_length_fields(
    annotated: pd.DataFrame,
    report_length_basis: str = "context",
    manual_report_shift: float = 0.0,
) -> pd.DataFrame:
    if annotated is None or len(annotated) == 0:
        return annotated
    d = annotated.copy()
    basis = (
        report_length_basis
        if report_length_basis in ("a_count", "context", "full_fragment", "manual_shift")
        else "context"
    )
    manual = pd.to_numeric(manual_report_shift, errors="coerce")
    manual = float(manual) if pd.notna(manual) else 0.0

    a_count = pd.to_numeric(d.get("A_Count"), errors="coerce")
    context = (
        pd.to_numeric(d["FixedContextLength"], errors="coerce")
        if "FixedContextLength" in d
        else pd.Series(0, index=d.index)
    ).fillna(0)
    non_a = (
        pd.to_numeric(d["NonA_Count"], errors="coerce")
        if "NonA_Count" in d
        else pd.Series(0, index=d.index)
    ).fillna(0)
    fragment_length = (
        pd.to_numeric(d["FragmentLength"], errors="coerce")
        if "FragmentLength" in d
        else pd.Series(np.nan, index=d.index)
    )

    if basis == "a_count":
        report_length = a_count
        applied_shift = pd.Series(0, index=d.index)
    elif basis == "context":
        report_length = a_count + context
        applied_shift = context
    elif basis == "full_fragment":
        report_length = fragment_length
        applied_shift = context + non_a
    else:
        report_length = a_count + manual
        applied_shift = pd.Series(manual, index=d.index)

    base_min = pd.to_numeric(d.get("PlotMin"), errors="coerce")
    base_max = pd.to_numeric(d.get("PlotMax"), errors="coerce")
    d["ReportLength"] = report_length.round().astype("Int64")
    d["ReportLengthBasis"] = report_length_basis_label(basis)
    d["ReportLengthShift"] = applied_shift
    d["ReportPlotMin"] = base_min + applied_shift
    d["ReportPlotMax"] = base_max + applied_shift
    return d


# ---------------------------------------------------------------------------
# Matching and terminal-state diagnostics
# ---------------------------------------------------------------------------

_CLASS_PRIORITY = {
    "A_only": 1,
    "C_containing": 2,
    "U_containing": 3,
    "mixed_CU": 4,
    "length_only": 98,
    "unassigned": 99,
}


def _class_priority(species_class):
    return _CLASS_PRIORITY.get(str(species_class), 50)


def phosphate_offset_diagnostic(error_da: float, tolerance_da: float = 2.0) -> tuple[bool, str, str]:
    if error_da is None or not np.isfinite(error_da):
        return False, "", ""
    if abs(abs(float(error_da)) - HPO3_MASS) > float(tolerance_da):
        return False, "", ""
    if error_da < 0:
        direction = "observed_lower"
        message = (
            "Observed mass is approximately one HPO3 unit (79.966 Da) lower than the "
            "selected theoretical species. Review whether one terminal phosphate should be OH."
        )
    else:
        direction = "observed_higher"
        message = (
            "Observed mass is approximately one HPO3 unit (79.966 Da) higher than the "
            "selected theoretical species. Review whether one terminal OH should be phosphate."
        )
    return True, direction, message


def _candidate_label(row: pd.Series, error_da: float) -> str:
    placement = row.get("VariantPlacementLabel", "")
    termini = row.get("Termini", "")
    extra = ", ".join(x for x in (str(placement), str(termini)) if x and x != "nan")
    return f"{row['DisplayLabel']} [{row['TailBlockID']}; {extra}] ({error_da:.4f} Da)"


def _copy_library_row(row: pd.Series) -> dict:
    return {column: row[column] for column in row.index}


def match_to_library(
    filtered: pd.DataFrame,
    lib: pd.DataFrame,
    tolerance_value: float = 10,
    tolerance_unit: str = "ppm",
    fill_unassigned: bool = True,
    fill_tolerance_da: float = 35,
    length_only_tolerance_da: float = 165,
    legacy_nearest_ladder: bool = True,
) -> pd.DataFrame:
    if filtered is None or lib is None or len(filtered) == 0 or len(lib) == 0:
        return pd.DataFrame()

    effective_length_tolerance = max(
        float(length_only_tolerance_da) if np.isfinite(length_only_tolerance_da) else 165,
        165,
    )
    lib_masses = pd.to_numeric(lib["TheoreticalMass"], errors="coerce").to_numpy(dtype=float)
    output_rows = []

    for _, observed in filtered.iterrows():
        observed_mass = float(observed["Mass"])
        error_da_all = observed_mass - lib_masses
        absolute_da = np.abs(error_da_all)
        error_ppm_all = 1e6 * error_da_all / lib_masses
        absolute_ppm = np.abs(error_ppm_all)
        in_tolerance = (
            absolute_ppm <= float(tolerance_value)
            if tolerance_unit == "ppm"
            else absolute_da <= float(tolerance_value)
        )
        candidate_indices = np.where(in_tolerance)[0]

        status = None
        if len(candidate_indices) == 0 and fill_unassigned:
            nearest_index = int(np.nanargmin(absolute_da))
            if np.isfinite(absolute_da[nearest_index]) and absolute_da[nearest_index] <= float(fill_tolerance_da):
                candidate_indices = np.array([nearest_index])
                status = "nearest_fill"
            elif np.isfinite(absolute_da[nearest_index]) and absolute_da[nearest_index] <= effective_length_tolerance:
                status = "length_only"
            else:
                status = "unassigned"
        elif len(candidate_indices) == 0:
            nearest_index = int(np.nanargmin(absolute_da))
            status = (
                "length_only"
                if np.isfinite(absolute_da[nearest_index])
                and absolute_da[nearest_index] <= effective_length_tolerance
                else "unassigned"
            )
        else:
            status = "unique" if len(candidate_indices) == 1 else "ambiguous"

        if status in ("length_only", "unassigned"):
            nearest_index = int(np.nanargmin(absolute_da))
            nearest = lib.iloc[nearest_index]
            nearest_error_da = float(observed_mass - nearest["TheoreticalMass"])
            nearest_error_ppm = 1e6 * nearest_error_da / float(nearest["TheoreticalMass"])
            flag, direction, diagnostic = phosphate_offset_diagnostic(nearest_error_da)

            if legacy_nearest_ladder:
                row = observed.to_dict()
                row.update(_copy_library_row(nearest))
                row.update(
                    {
                        "MassErrorDa": round(nearest_error_da, 6),
                        "MassErrorPpm": round(nearest_error_ppm, 3),
                        "AnnotationStatus": "nearest_ladder",
                        "CandidateCount": 1,
                        "CandidateList": _candidate_label(nearest, nearest_error_da) + " [nearest ladder]",
                        "NearestSpecies": nearest["DisplayLabel"],
                        "NearestTheoreticalMass": nearest["TheoreticalMass"],
                        "NearestErrorDa": round(nearest_error_da, 6),
                        "NearestErrorPpm": round(nearest_error_ppm, 3),
                        "PhosphateOffsetFlag": flag,
                        "PhosphateOffsetDirection": direction,
                        "PhosphateOffsetMessage": diagnostic,
                    }
                )
                output_rows.append(row)
                continue

            row = observed.to_dict()
            row.update(_copy_library_row(nearest))
            row.update(
                {
                    "SpeciesID": np.nan,
                    "SpeciesClass": "length_only" if status == "length_only" else "unassigned",
                    "SpeciesClassLabel": (
                        "Nearest length only" if status == "length_only" else "Unassigned"
                    ),
                    "C_Count": np.nan,
                    "U_Count": np.nan,
                    "NonA_Count": np.nan,
                    "GeneratedSequence": np.nan,
                    "Composition": np.nan,
                    "PositionInterpretation": (
                        "Nearest length only; composition not assigned"
                        if status == "length_only"
                        else "Outside length-only tolerance; composition not assigned"
                    ),
                    "TheoreticalMass": np.nan,
                    "MassErrorDa": np.nan,
                    "MassErrorPpm": np.nan,
                    "AnnotationStatus": status,
                    "CandidateCount": 0,
                    "CandidateList": np.nan,
                    "NearestSpecies": nearest["DisplayLabel"],
                    "NearestTheoreticalMass": nearest["TheoreticalMass"],
                    "NearestErrorDa": round(nearest_error_da, 6),
                    "NearestErrorPpm": round(nearest_error_ppm, 3),
                    "PhosphateOffsetFlag": flag,
                    "PhosphateOffsetDirection": direction,
                    "PhosphateOffsetMessage": diagnostic,
                }
            )
            output_rows.append(row)
            continue

        candidates = lib.iloc[candidate_indices].copy()
        candidates["MassErrorDa"] = observed_mass - candidates["TheoreticalMass"]
        candidates["AbsMassErrorDa"] = candidates["MassErrorDa"].abs()
        candidates["MassErrorPpm"] = 1e6 * candidates["MassErrorDa"] / candidates["TheoreticalMass"]
        candidates["ClassPriority"] = candidates["SpeciesClass"].map(_class_priority)
        candidates = candidates.sort_values(
            ["AbsMassErrorDa", "ClassPriority", "NonA_Count", "VariantPlacement"],
            kind="stable",
        )
        best = candidates.iloc[0]
        candidate_list = "; ".join(
            _candidate_label(candidate, float(candidate["MassErrorDa"]))
            for _, candidate in candidates.head(8).iterrows()
        )
        best_error_da = float(best["MassErrorDa"])
        best_error_ppm = float(best["MassErrorPpm"])
        flag, direction, diagnostic = phosphate_offset_diagnostic(best_error_da)

        row = observed.to_dict()
        row.update(_copy_library_row(best))
        row.update(
            {
                "MassErrorDa": round(best_error_da, 6),
                "MassErrorPpm": round(best_error_ppm, 3),
                "AnnotationStatus": status,
                "CandidateCount": len(candidates),
                "CandidateList": candidate_list,
                "NearestSpecies": best["DisplayLabel"],
                "NearestTheoreticalMass": best["TheoreticalMass"],
                "NearestErrorDa": round(best_error_da, 6),
                "NearestErrorPpm": round(best_error_ppm, 3),
                "PhosphateOffsetFlag": flag,
                "PhosphateOffsetDirection": direction,
                "PhosphateOffsetMessage": diagnostic,
            }
        )
        output_rows.append(row)

    return pd.DataFrame(output_rows)


# ---------------------------------------------------------------------------
# Summary and AWDI
# ---------------------------------------------------------------------------


def calculate_awdi_from_annotated(annotated: pd.DataFrame, basis: str = "A_only") -> float:
    if annotated is None or len(annotated) == 0:
        return float("nan")
    d = annotated[
        (annotated["AnnotationStatus"] != "unassigned")
        & annotated["ExpectedAOnlyMass"].apply(lambda x: pd.notna(x) and np.isfinite(x))
    ].copy()
    if basis == "A_only":
        d = d[d["SpeciesClass"] == "A_only"]
    elif basis == "report_default":
        d = d[d["SpeciesClass"].isin(["A_only", "C_containing"])]
    total_weight = d["Fractional.Abundance"].sum()
    if len(d) == 0 or total_weight <= 0:
        return float("nan")
    weights = d["Fractional.Abundance"] / total_weight
    return float(
        (((d["Mass"] - d["ExpectedAOnlyMass"]) / d["ExpectedAOnlyMass"]) ** 2 * weights).sum()
    )


def class_pct(data: pd.DataFrame, species_class: str) -> float:
    assigned = data[data["AnnotationStatus"] != "unassigned"]
    total = assigned["Fractional.Abundance"].sum()
    if total <= 0:
        return float("nan")
    return 100 * assigned.loc[
        assigned["SpeciesClass"] == species_class, "Fractional.Abundance"
    ].sum() / total


def unassigned_pct(data: pd.DataFrame) -> float:
    total = data["Fractional.Abundance"].sum()
    if total <= 0:
        return float("nan")
    return 100 * data.loc[
        data["AnnotationStatus"] == "unassigned", "Fractional.Abundance"
    ].sum() / total


def build_summary_table(results: dict, awdi_basis: str) -> pd.DataFrame:
    rows = []
    for name, result in results.items():
        data = result.get("annotated")
        if data is None or len(data) == 0:
            rows.append({"File": name, "FilteredSpecies": 0, "AssignedSpecies": 0})
            continue
        assigned = data[data["AnnotationStatus"] != "unassigned"]
        most = assigned.iloc[assigned["Fractional.Abundance"].argmax()] if len(assigned) else None
        report_lengths = (
            assigned["ReportLength"] if "ReportLength" in assigned else assigned.get("AssignedLength")
        )
        rows.append(
            {
                "File": name,
                "AWDI": (
                    round(result.get("awdi", float("nan")), 6)
                    if result.get("awdi") == result.get("awdi")
                    else None
                ),
                "AWDIBasis": awdi_basis,
                "FilteredSpecies": len(data),
                "AssignedSpecies": len(assigned),
                "UniqueAssigned": int((data["AnnotationStatus"] == "unique").sum()),
                "Ambiguous": int((data["AnnotationStatus"] == "ambiguous").sum()),
                "NearestFilled": int((data["AnnotationStatus"] == "nearest_fill").sum()),
                "NearestLadder": int((data["AnnotationStatus"] == "nearest_ladder").sum()),
                "LengthOnly": int((data["AnnotationStatus"] == "length_only").sum()),
                "Unassigned": int((data["AnnotationStatus"] == "unassigned").sum()),
                "PhosphateOffsetFlags": int(
                    pd.Series(data.get("PhosphateOffsetFlag", False)).fillna(False).astype(bool).sum()
                ),
                "Unassigned % of total filtered abundance": round(unassigned_pct(data), 2),
                "Percent Poly(A)": round(class_pct(data, "A_only"), 2),
                "Percent Poly(A) with C": round(class_pct(data, "C_containing"), 2),
                "Percent Poly(A) with U": round(class_pct(data, "U_containing"), 2),
                "Percent Poly(A) with C/U": round(class_pct(data, "mixed_CU"), 2),
                "MostAbundantSpecies": most["DisplayLabel"] if most is not None else None,
                "MostAbundantClass": most["SpeciesClassLabel"] if most is not None else None,
                "MostAbundantA_Count": most["A_Count"] if most is not None else None,
                "MostAbundantReportLength": (
                    most["ReportLength"]
                    if most is not None and "ReportLength" in most
                    else (most["AssignedLength"] if most is not None else None)
                ),
                "ReportLengthBasis": (
                    most["ReportLengthBasis"]
                    if most is not None and "ReportLengthBasis" in most
                    else None
                ),
                "SelectedTermini": most["Termini"] if most is not None else None,
                "VariantPlacement": (
                    most.get("VariantPlacementLabel") if most is not None else None
                ),
                "MedianMass": round(float(data["Mass"].median()), 4),
                "MedianA_Count": round(float(assigned["A_Count"].median()), 1) if len(assigned) else None,
                "MedianReportLength": (
                    round(float(report_lengths.median()), 1)
                    if len(assigned) and report_lengths is not None
                    else None
                ),
                "A_CountRange": (
                    f"{int(assigned['A_Count'].min())}–{int(assigned['A_Count'].max())}"
                    if len(assigned)
                    else None
                ),
                "ReportLengthRange": (
                    f"{int(report_lengths.min())}–{int(report_lengths.max())}"
                    if len(assigned) and report_lengths is not None
                    else None
                ),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Mode-mass QC
# ---------------------------------------------------------------------------


def build_mode_mass_qc_table(
    results: dict,
    enabled: bool = True,
    tolerance_value: float = 2,
    tolerance_unit: str = "Da",
) -> pd.DataFrame:
    unit = tolerance_unit if tolerance_unit in ("Da", "ppm") else "Da"
    tolerance = pd.to_numeric(tolerance_value, errors="coerce")
    if not pd.notna(tolerance) or tolerance < 0:
        tolerance = 20 if unit == "ppm" else 1

    rows = []
    for name, result in results.items():
        data = result.get("annotated")
        if data is None or len(data) == 0:
            rows.append(
                {
                    "File": name,
                    "ModeQCEnabled": bool(enabled),
                    "ModeQCPass": not bool(enabled),
                    "ModeObservedMass": None,
                    "ModeApexRT": None,
                    "ModeRelativeAbundance": None,
                    "ModeNearestSpecies": None,
                    "ModeSpeciesClass": None,
                    "ModeReportLength": None,
                    "ModeAnnotationStatus": None,
                    "ModeTheoreticalMass": None,
                    "ModeMassErrorDa": None,
                    "ModeAbsMassErrorDa": None,
                    "ModeMassErrorPpm": None,
                    "ModeAbsMassErrorPpm": None,
                    "ModeQCTolerance": tolerance,
                    "ModeQCToleranceUnit": unit,
                    "ModePhosphateOffsetFlag": False,
                    "ModeQCMessage": (
                        "FAIL: no annotated filtered data available." if enabled else "QC disabled."
                    ),
                }
            )
            continue

        mode = data.iloc[data["Fractional.Abundance"].argmax()]
        theoretical = pd.to_numeric(mode.get("TheoreticalMass"), errors="coerce")
        if not pd.notna(theoretical):
            theoretical = pd.to_numeric(mode.get("NearestTheoreticalMass"), errors="coerce")
        observed_mass = pd.to_numeric(mode.get("Mass"), errors="coerce")
        if pd.notna(theoretical) and pd.notna(observed_mass):
            error_da = float(observed_mass - theoretical)
            error_ppm = 1e6 * error_da / float(theoretical)
        else:
            error_da = np.nan
            error_ppm = np.nan
        absolute_da = abs(error_da) if np.isfinite(error_da) else np.nan
        absolute_ppm = abs(error_ppm) if np.isfinite(error_ppm) else np.nan
        phosphate_flag, _, phosphate_message = phosphate_offset_diagnostic(error_da)

        passed = True
        message = "QC disabled."
        if enabled:
            if not np.isfinite(observed_mass) or not np.isfinite(theoretical):
                passed = False
                message = (
                    "FAIL: mode peak could not be compared with the generated ladder. "
                    "Check sequence, range, and library settings."
                )
            elif unit == "ppm":
                passed = np.isfinite(absolute_ppm) and absolute_ppm <= tolerance
                message = (
                    f"PASS: mode peak matches generated ladder within {tolerance} ppm."
                    if passed
                    else (
                        f"FAIL: mode peak differs by {absolute_ppm:.2f} ppm / {absolute_da:.3f} Da. "
                        "Check sequence conversion, U chemistry, termini, C/U placement, or range."
                    )
                )
            else:
                passed = np.isfinite(absolute_da) and absolute_da <= tolerance
                message = (
                    f"PASS: mode peak matches generated ladder within {tolerance} Da."
                    if passed
                    else (
                        f"FAIL: mode peak differs by {absolute_da:.3f} Da / {absolute_ppm:.2f} ppm. "
                        "Check sequence conversion, U chemistry, termini, C/U placement, or range."
                    )
                )
            if phosphate_flag and not passed:
                message += " " + phosphate_message

        rows.append(
            {
                "File": name,
                "ModeQCEnabled": bool(enabled),
                "ModeQCPass": bool(passed),
                "ModeObservedMass": round(float(observed_mass), 6) if pd.notna(observed_mass) else None,
                "ModeApexRT": float(mode.get("Apex.RT")) if pd.notna(mode.get("Apex.RT")) else None,
                "ModeRelativeAbundance": (
                    float(mode.get("Fractional.Abundance"))
                    if pd.notna(mode.get("Fractional.Abundance"))
                    else None
                ),
                "ModeNearestSpecies": (
                    str(mode.get("DisplayLabel")) if pd.notna(mode.get("DisplayLabel")) else None
                ),
                "ModeSpeciesClass": (
                    str(mode.get("SpeciesClassLabel"))
                    if pd.notna(mode.get("SpeciesClassLabel"))
                    else None
                ),
                "ModeReportLength": (
                    float(mode.get("ReportLength"))
                    if "ReportLength" in mode and pd.notna(mode.get("ReportLength"))
                    else (
                        float(mode.get("AssignedLength"))
                        if pd.notna(mode.get("AssignedLength"))
                        else None
                    )
                ),
                "ModeAnnotationStatus": (
                    str(mode.get("AnnotationStatus"))
                    if pd.notna(mode.get("AnnotationStatus"))
                    else None
                ),
                "ModeTheoreticalMass": (
                    round(float(theoretical), 6) if pd.notna(theoretical) else None
                ),
                "ModeMassErrorDa": round(error_da, 6) if np.isfinite(error_da) else None,
                "ModeAbsMassErrorDa": round(absolute_da, 6) if np.isfinite(absolute_da) else None,
                "ModeMassErrorPpm": round(error_ppm, 3) if np.isfinite(error_ppm) else None,
                "ModeAbsMassErrorPpm": (
                    round(absolute_ppm, 3) if np.isfinite(absolute_ppm) else None
                ),
                "ModeQCTolerance": tolerance,
                "ModeQCToleranceUnit": unit,
                "ModePhosphateOffsetFlag": phosphate_flag,
                "ModeQCMessage": message,
            }
        )
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


def build_class_plot_data(
    annotated: pd.DataFrame,
    include_unassigned: bool = False,
    duplicate_rule: str = "dominant",
) -> pd.DataFrame:
    if annotated is None or len(annotated) == 0:
        return pd.DataFrame()
    status_keep = {"unique", "ambiguous", "nearest_fill", "nearest_ladder"}
    length_col = "ReportLength" if "ReportLength" in annotated.columns else "AssignedLength"
    data = (
        annotated
        if include_unassigned
        else annotated[
            annotated["AnnotationStatus"].isin(status_keep) & annotated[length_col].notna()
        ]
    )
    if len(data) == 0:
        return pd.DataFrame()

    data = data.copy()
    data["PlotLength"] = pd.to_numeric(data[length_col], errors="coerce")
    data["PlotMinForReport"] = pd.to_numeric(
        data.get("ReportPlotMin", data.get("PlotMin")), errors="coerce"
    )
    data["PlotMaxForReport"] = pd.to_numeric(
        data.get("ReportPlotMax", data.get("PlotMax")), errors="coerce"
    )
    data["SpeciesLabel"] = data["SpeciesClass"].map(class_to_label)
    if "SampleLabel" not in data.columns:
        data["SampleLabel"] = data["File"].map(clean_sample_label)
    data = data[
        data["PlotLength"].apply(lambda x: pd.notna(x) and np.isfinite(x))
        & data["Fractional.Abundance"].apply(lambda x: pd.notna(x) and np.isfinite(x))
    ]
    if len(data) == 0:
        return pd.DataFrame()

    aggregate_function = "sum" if duplicate_rule == "sum" else "max"
    return data.groupby(
        ["File", "SampleLabel", "TailBlockID", "PlotLength", "SpeciesClass", "SpeciesLabel"],
        dropna=False,
        as_index=False,
    ).agg(
        RelAbundance=("Fractional.Abundance", aggregate_function),
        PlotMin=("PlotMinForReport", "first"),
        PlotMax=("PlotMaxForReport", "first"),
    )


def build_total_plot_data_from_class(
    class_data: pd.DataFrame, threshold_pct: float = 10
) -> pd.DataFrame:
    if class_data is None or len(class_data) == 0:
        return pd.DataFrame()
    total = class_data.groupby(["File", "SampleLabel", "PlotLength"], as_index=False).agg(
        TotalRaw=("RelAbundance", "sum"),
        PlotMin=("PlotMin", "first"),
        PlotMax=("PlotMax", "first"),
    )
    maxima = total.groupby(["File", "SampleLabel"], as_index=False).agg(
        MaxRaw=("TotalRaw", "max")
    )
    total = total.merge(maxima, on=["File", "SampleLabel"], how="left")
    total["TotalNorm"] = np.where(
        (total["MaxRaw"] > 0) & total["MaxRaw"].notna(),
        100 * total["TotalRaw"] / total["MaxRaw"],
        0,
    )
    total["ThresholdClass"] = np.where(
        total["TotalNorm"] >= float(threshold_pct), "above_10", "below_10"
    )
    return total
