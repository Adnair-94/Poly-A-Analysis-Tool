import math

import pandas as pd

import polya_core as pc


def test_sequence_conversion_modes():
    assert pc.convert_sequence_to_rna("ATGC", "dna_coding") == "AUGC"
    assert pc.convert_sequence_to_rna("GCAT", "dna_template") == "AUGC"
    assert pc.convert_sequence_to_rna("AUGC", "rna") == "AUGC"
    assert pc.rna_to_dna("AUGC") == "ATGC"


def test_terminal_mass_differences_are_hpo3():
    sequence = "ACGUAAAA"
    oh_oh = pc.calculate_rna_mass(sequence, five_prime="5OH", three_prime="3OH")
    oh_p = pc.calculate_rna_mass(sequence, five_prime="5OH", three_prime="3p")
    p_oh = pc.calculate_rna_mass(sequence, five_prime="5p", three_prime="3OH")
    p_p = pc.calculate_rna_mass(sequence, five_prime="5p", three_prime="3p")

    assert math.isclose(oh_p - oh_oh, pc.HPO3_MASS, abs_tol=1e-9)
    assert math.isclose(p_oh - oh_oh, pc.HPO3_MASS, abs_tol=1e-9)
    assert math.isclose(p_p - oh_p, pc.HPO3_MASS, abs_tol=1e-9)


def test_cyclic_phosphate_mass_convention():
    sequence = "ACGUAAAA"
    oh_oh = pc.calculate_rna_mass(sequence, five_prime="5OH", three_prime="3OH")
    cyclic = pc.calculate_rna_mass(sequence, five_prime="5OH", three_prime="3cyclicp")
    assert math.isclose(cyclic - oh_oh, pc.HPO3_MASS - pc.H2O_MASS, abs_tol=1e-9)


def _library(placement="three_prime_of_a", five="auto", three="auto", max_c=1, max_u=2):
    return pc.generate_species_library(
        sample_id="Test",
        seq="GCUUAAAAAAAAAA",
        sequence_mode="tail_fragment",
        sequence_input_type="rna",
        use_auto_ranges=False,
        expected_a_length=10,
        range_min=10,
        range_max=10,
        range_input_basis="a_count",
        report_length_basis="a_count",
        min_a_run=5,
        max_c=max_c,
        max_u=max_u,
        allow_mixed=True,
        variant_placement=placement,
        five_prime_terminus=five,
        three_prime_terminus=three,
    )


def test_variant_position_changes_sequence_not_intact_mass():
    lib_5 = _library("five_prime_of_a").SpeciesLibrary
    lib_3 = _library("three_prime_of_a").SpeciesLibrary

    row_5 = lib_5[(lib_5.C_Count == 1) & (lib_5.U_Count == 2)].iloc[0]
    row_3 = lib_3[(lib_3.C_Count == 1) & (lib_3.U_Count == 2)].iloc[0]

    assert row_5.GeneratedSequence != row_3.GeneratedSequence
    assert math.isclose(row_5.TheoreticalMass, row_3.TheoreticalMass, abs_tol=1e-9)


def test_manual_terminal_override_changes_generated_ladder():
    lib_oh = _library("three_prime_of_a", three="3OH", max_c=0, max_u=0)
    lib_p = _library("three_prime_of_a", three="3p", max_c=0, max_u=0)
    mass_oh = float(lib_oh.SpeciesLibrary.iloc[0].TheoreticalMass)
    mass_p = float(lib_p.SpeciesLibrary.iloc[0].TheoreticalMass)
    assert math.isclose(mass_p - mass_oh, pc.HPO3_MASS, abs_tol=2e-6)
    assert lib_p.TailBlocks.iloc[0].TerminalOverrideApplied


def test_phosphate_offset_diagnostic():
    flag, direction, message = pc.phosphate_offset_diagnostic(-pc.HPO3_MASS)
    assert flag
    assert direction == "observed_lower"
    assert "phosphate" in message.lower()


def test_exact_matching_and_report_length():
    generated = _library("three_prime_of_a", max_c=0, max_u=0)
    theoretical = float(generated.SpeciesLibrary.iloc[0].TheoreticalMass)
    observed = pd.DataFrame(
        {"Mass": [theoretical], "Apex.RT": [12.3], "Fractional.Abundance": [100.0]}
    )
    annotated = pc.match_to_library(observed, generated.SpeciesLibrary, tolerance_value=5)
    annotated = pc.add_report_length_fields(annotated, "a_count")
    assert annotated.iloc[0].AnnotationStatus == "unique"
    assert int(annotated.iloc[0].ReportLength) == 10
    assert not bool(annotated.iloc[0].PhosphateOffsetFlag)


def test_minus_hpo3_observation_is_flagged():
    generated = _library("three_prime_of_a", three="3p", max_c=0, max_u=0)
    theoretical = float(generated.SpeciesLibrary.iloc[0].TheoreticalMass)
    observed = pd.DataFrame(
        {
            "Mass": [theoretical - pc.HPO3_MASS],
            "Apex.RT": [12.3],
            "Fractional.Abundance": [100.0],
        }
    )
    annotated = pc.match_to_library(
        observed,
        generated.SpeciesLibrary,
        tolerance_value=5,
        fill_tolerance_da=1,
        length_only_tolerance_da=165,
        legacy_nearest_ladder=True,
    )
    assert bool(annotated.iloc[0].PhosphateOffsetFlag)
    assert annotated.iloc[0].PhosphateOffsetDirection == "observed_lower"


def test_both_placements_are_explicitly_isobaric():
    generated = _library("both")
    subset = generated.SpeciesLibrary[
        (generated.SpeciesLibrary.C_Count == 1) & (generated.SpeciesLibrary.U_Count == 2)
    ]
    assert set(subset.VariantPlacement) == {"five_prime_of_a", "three_prime_of_a"}
    assert subset.TheoreticalMass.nunique() == 1
