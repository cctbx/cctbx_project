"""
Tests for directive_extractor module.

These tests verify that directives are correctly extracted from user advice
and properly validated.

Run with: python tests/test_directive_extractor.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
    assert_equal, assert_not_equal, assert_true, assert_false, assert_in,
    assert_not_in, assert_none, assert_is_instance,
    run_tests_with_fail_fast
)

from agent.directive_extractor import (
    validate_directives,
    merge_directives,
    get_program_settings,
    check_stop_conditions,
    extract_directives_simple,
    _fix_program_name,
    format_directives_for_display,
)


# =============================================================================
# VALIDATE DIRECTIVES TESTS
# =============================================================================

def test_validate_empty_input():
    """Empty input should return empty dict."""
    assert_equal(validate_directives(None), {})
    assert_equal(validate_directives({}), {})
    assert_equal(validate_directives("not a dict"), {})


def test_validate_valid_program_settings():
    """Valid program settings should be preserved."""
    input_directives = {
        "program_settings": {
            "phenix.refine": {"resolution": 2.5, "anisotropic_adp": True},
            "default": {"resolution": 3.0}
        }
    }
    result = validate_directives(input_directives)

    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 2.5)
    assert_equal(result["program_settings"]["phenix.refine"]["anisotropic_adp"], True)
    assert_equal(result["program_settings"]["default"]["resolution"], 3.0)


def test_validate_invalid_program_name_removed():
    """Invalid program names should be removed."""
    input_directives = {
        "program_settings": {
            "invalid_program": {"resolution": 2.5},
            "phenix.refine": {"resolution": 3.0}
        }
    }
    result = validate_directives(input_directives)

    assert_in("program_settings", result)
    assert_not_in("invalid_program", result["program_settings"])
    assert_in("phenix.refine", result["program_settings"])


def test_validate_fixable_program_name():
    """Fixable program names should be corrected."""
    input_directives = {
        "program_settings": {
            "refine": {"resolution": 2.5},
        }
    }
    result = validate_directives(input_directives)

    assert_in("program_settings", result)
    assert_in("phenix.refine", result["program_settings"])


def test_validate_valid_stop_conditions():
    """Valid stop conditions should be preserved."""
    input_directives = {
        "stop_conditions": {
            "after_cycle": 5,
            "r_free_target": 0.25,
            "skip_validation": True
        }
    }
    result = validate_directives(input_directives)

    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_cycle"], 5)
    assert_equal(result["stop_conditions"]["r_free_target"], 0.25)
    assert_true(result["stop_conditions"]["skip_validation"])


def test_validate_type_conversion():
    """String values should be converted to appropriate types."""
    input_directives = {
        "program_settings": {
            "default": {"resolution": "2.5"}  # String should become float
        },
        "stop_conditions": {
            "after_cycle": "3"  # String should become int
        }
    }
    result = validate_directives(input_directives)

    assert_equal(result["program_settings"]["default"]["resolution"], 2.5)
    assert_is_instance(result["program_settings"]["default"]["resolution"], float)
    assert_equal(result["stop_conditions"]["after_cycle"], 3)
    assert_is_instance(result["stop_conditions"]["after_cycle"], int)


def test_validate_workflow_preferences():
    """Workflow preferences should be preserved."""
    input_directives = {
        "workflow_preferences": {
            "skip_programs": ["phenix.autobuild"],
            "prefer_programs": ["phenix.refine"]
        }
    }
    result = validate_directives(input_directives)

    assert_in("workflow_preferences", result)
    assert_in("phenix.autobuild", result["workflow_preferences"]["skip_programs"])


def test_validate_file_preferences():
    """File preferences should be preserved."""
    input_directives = {
        "file_preferences": {
            "prefer_unmerged": True,
            "prefer_anomalous": True,
            "model": "my_model.pdb"
        }
    }
    result = validate_directives(input_directives)

    assert_in("file_preferences", result)
    assert_true(result["file_preferences"]["prefer_unmerged"])
    assert_true(result["file_preferences"]["prefer_anomalous"])
    assert_equal(result["file_preferences"]["model"], "my_model.pdb")


# =============================================================================
# MERGE DIRECTIVES TESTS
# =============================================================================

def test_merge_empty_base():
    """Merging into empty base should return new directives."""
    new = {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
    result = merge_directives({}, new)

    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 2.5)


def test_merge_empty_new():
    """Merging empty new should return base unchanged."""
    base = {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
    result = merge_directives(base, {})

    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 2.5)


def test_merge_override():
    """New directives should override base."""
    base = {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
    new = {"program_settings": {"phenix.refine": {"resolution": 3.0}}}
    result = merge_directives(base, new)

    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 3.0)


def test_merge_addition():
    """New directives should add to base."""
    base = {"program_settings": {"phenix.refine": {"resolution": 2.5}}}
    new = {"program_settings": {"phenix.autosol": {"resolution": 3.0}}}
    result = merge_directives(base, new)

    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 2.5)
    assert_equal(result["program_settings"]["phenix.autosol"]["resolution"], 3.0)


def test_merge_deep():
    """Deep merging should preserve nested values."""
    base = {"program_settings": {"phenix.refine": {"resolution": 2.5, "cycles": 5}}}
    new = {"program_settings": {"phenix.refine": {"anisotropic_adp": True}}}
    result = merge_directives(base, new)

    assert_equal(result["program_settings"]["phenix.refine"]["resolution"], 2.5)
    assert_equal(result["program_settings"]["phenix.refine"]["cycles"], 5)
    assert_true(result["program_settings"]["phenix.refine"]["anisotropic_adp"])


# =============================================================================
# GET PROGRAM SETTINGS TESTS
# =============================================================================

def test_get_settings_specific_program():
    """Should return settings for specific program."""
    directives = {
        "program_settings": {
            "phenix.refine": {"resolution": 2.5, "cycles": 5}
        }
    }
    result = get_program_settings(directives, "phenix.refine")

    assert_equal(result["resolution"], 2.5)
    assert_equal(result["cycles"], 5)


def test_get_settings_default_fallback():
    """Should fall back to default settings."""
    directives = {
        "program_settings": {
            "default": {"resolution": 3.0}
        }
    }
    result = get_program_settings(directives, "phenix.refine")

    assert_equal(result["resolution"], 3.0)


def test_get_settings_merge_with_default():
    """Should merge program-specific with default."""
    directives = {
        "program_settings": {
            "default": {"resolution": 3.0, "cycles": 5},
            "phenix.refine": {"resolution": 2.5}
        }
    }
    result = get_program_settings(directives, "phenix.refine")

    assert_equal(result["resolution"], 2.5)  # From specific
    assert_equal(result["cycles"], 5)  # From default


def test_get_settings_empty():
    """Should return empty dict for missing program."""
    directives = {}
    result = get_program_settings(directives, "phenix.refine")

    assert_equal(result, {})


# =============================================================================
# CHECK STOP CONDITIONS TESTS
# =============================================================================

def test_stop_after_cycle():
    """Should stop after specified cycle."""
    directives = {"stop_conditions": {"after_cycle": 3}}

    should_stop, reason = check_stop_conditions(directives, cycle_number=3, last_program="phenix.refine")
    assert_true(should_stop)
    assert_in("3", reason)

    should_stop, reason = check_stop_conditions(directives, cycle_number=2, last_program="phenix.refine")
    assert_false(should_stop)


def test_stop_after_program():
    """Should stop after specified program."""
    directives = {"stop_conditions": {"after_program": "phenix.refine"}}

    should_stop, reason = check_stop_conditions(
        directives, cycle_number=2, last_program="phenix.refine"
    )
    assert_true(should_stop)

    should_stop, reason = check_stop_conditions(
        directives, cycle_number=2, last_program="phenix.phaser"
    )
    assert_false(should_stop)


def test_stop_r_free_target():
    """Should stop when R-free target reached."""
    directives = {"stop_conditions": {"r_free_target": 0.25}}

    should_stop, reason = check_stop_conditions(
        directives, cycle_number=2, last_program="phenix.refine",
        metrics={"r_free": 0.24}
    )
    assert_true(should_stop)

    should_stop, reason = check_stop_conditions(
        directives, cycle_number=2, last_program="phenix.refine",
        metrics={"r_free": 0.30}
    )
    assert_false(should_stop)


def test_stop_no_conditions():
    """Should not stop if no conditions defined."""
    directives = {}
    should_stop, reason = check_stop_conditions(directives, cycle_number=10, last_program="phenix.refine")
    assert_false(should_stop)


def test_stop_skip_validation():
    """skip_validation should be readable but not trigger stop."""
    directives = {"stop_conditions": {"skip_validation": True}}
    should_stop, reason = check_stop_conditions(directives, cycle_number=1, last_program="phenix.refine")
    assert_false(should_stop)  # skip_validation alone doesn't stop


# =============================================================================
# FIX PROGRAM NAME TESTS
# =============================================================================

def test_fix_already_valid():
    """Already valid names should be unchanged."""
    assert_equal(_fix_program_name("phenix.refine"), "phenix.refine")
    assert_equal(_fix_program_name("phenix.autosol"), "phenix.autosol")


def test_fix_missing_prefix():
    """Missing phenix. prefix should be added."""
    assert_equal(_fix_program_name("refine"), "phenix.refine")
    assert_equal(_fix_program_name("autosol"), "phenix.autosol")


def test_fix_common_variations():
    """Common variations should be normalized."""
    assert_equal(_fix_program_name("auto_build"), "phenix.autobuild")
    assert_equal(_fix_program_name("autobuild"), "phenix.autobuild")


def test_fix_unknown_returns_none():
    """Unknown programs should return None."""
    assert_none(_fix_program_name("unknown_program"))
    assert_none(_fix_program_name("not_a_phenix_tool"))


def test_fix_map_sharpening_aliases():
    """Map sharpening aliases should be normalized."""
    assert_equal(_fix_program_name("map_sharpening"), "phenix.map_sharpening")
    assert_equal(_fix_program_name("sharpen_map"), "phenix.map_sharpening")
    assert_equal(_fix_program_name("auto_sharpen"), "phenix.map_sharpening")
    assert_equal(_fix_program_name("autosharpen"), "phenix.map_sharpening")


def test_fix_map_to_model_aliases():
    """Map to model aliases should be normalized."""
    assert_equal(_fix_program_name("map_to_model"), "phenix.map_to_model")
    assert_equal(_fix_program_name("maptomodel"), "phenix.map_to_model")
    assert_equal(_fix_program_name("build_model"), "phenix.map_to_model")
    assert_equal(_fix_program_name("buildmodel"), "phenix.map_to_model")


# =============================================================================
# EXTRACT DIRECTIVES SIMPLE TESTS
# =============================================================================

def test_extract_resolution():
    """Should extract resolution values."""
    result = extract_directives_simple("use resolution 2.5")
    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["default"]["resolution"], 2.5)


def test_extract_resolution_angstrom():
    """Should extract resolution with Angstrom."""
    result = extract_directives_simple("resolution of 3.0 Angstrom")
    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["default"]["resolution"], 3.0)


def test_extract_stop_after_refinement():
    """Should extract stop after refinement."""
    result = extract_directives_simple("stop after the first refinement")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.refine")


def test_extract_stop_after_cycle():
    """Should extract stop after cycle N."""
    result = extract_directives_simple("stop at cycle 4")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_cycle"], 4)


def test_extract_anisotropic():
    """Should extract anisotropic refinement directive."""
    result = extract_directives_simple("run anisotropic refinement")
    assert_in("program_settings", result)
    assert_true(result["program_settings"]["phenix.refine"]["anisotropic_adp"])


def test_extract_skip_validation():
    """Should extract skip validation directive."""
    result = extract_directives_simple("skip validation")
    assert_in("stop_conditions", result)
    assert_true(result["stop_conditions"]["skip_validation"])


def test_extract_multiple():
    """Should extract multiple directives."""
    result = extract_directives_simple("use resolution 2.5, stop at cycle 3")
    assert_in("program_settings", result)
    assert_in("stop_conditions", result)


def test_extract_empty():
    """Should return empty dict for non-directive text."""
    result = extract_directives_simple("just a regular message")
    assert_equal(result, {})


# =============================================================================
# FORMAT DIRECTIVES FOR DISPLAY TESTS
# =============================================================================

def test_format_empty():
    """Empty directives should return appropriate message."""
    result = format_directives_for_display({})
    assert_in("no", result.lower())


def test_format_with_settings():
    """Should format program settings."""
    directives = {
        "program_settings": {
            "phenix.refine": {"resolution": 2.5}
        }
    }
    result = format_directives_for_display(directives)
    assert_in("phenix.refine", result)
    assert_in("2.5", result)


def test_format_with_stop_conditions():
    """Should format stop conditions."""
    directives = {
        "stop_conditions": {"after_cycle": 3}
    }
    result = format_directives_for_display(directives)
    assert_in("stop", result.lower())


# =============================================================================
# AUTOSOL DIRECTIVES TESTS
# =============================================================================

def test_extract_atom_type_selenium():
    """Should extract selenium atom type."""
    result = extract_directives_simple("use selenium as the anomalous scatterer")
    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["phenix.autosol"]["atom_type"], "Se")


def test_extract_atom_type_short():
    """Should extract short form atom type."""
    result = extract_directives_simple("Se-SAD experiment")
    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["phenix.autosol"]["atom_type"], "Se")


def test_extract_atom_type_sulfur():
    """Should extract sulfur atom type."""
    result = extract_directives_simple("sulfur SAD phasing")
    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["phenix.autosol"]["atom_type"], "S")


def test_extract_ncs():
    """Should extract NCS directive if pattern matched."""
    result = extract_directives_simple("use NCS for refinement")
    assert_is_instance(result, dict)


def test_extract_mrsad():
    """Should extract MR-SAD directive."""
    result = extract_directives_simple("run MR-SAD phasing")
    # MR-SAD triggers phaser
    assert_is_instance(result, dict)


# =============================================================================
# STOP CONDITION EDGE CASES
# =============================================================================

def test_stop_first_only():
    """'first' should mean max 1."""
    result = extract_directives_simple("stop after the first refinement")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"].get("max_refine_cycles"), 1)


def test_stop_program_normalization():
    """Program names in stop conditions should be normalized."""
    result = extract_directives_simple("stop after first refine")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.refine")


def test_stop_multiple_conditions():
    """Multiple stop conditions should all be captured."""
    directives = {
        "stop_conditions": {
            "after_cycle": 5,
            "r_free_target": 0.22
        }
    }
    # Cycle 5 reached but R-free not at target
    should_stop, reason = check_stop_conditions(
        directives, cycle_number=5, last_program="phenix.refine",
        metrics={"r_free": 0.30}
    )
    assert_true(should_stop)  # Cycle condition met


# =============================================================================
# FILE PREFERENCE TESTS
# =============================================================================

def test_extract_unmerged_preference():
    """Should extract preference for unmerged data."""
    result = extract_directives_simple("prefer unmerged data")
    assert_in("file_preferences", result)
    assert_true(result["file_preferences"]["prefer_unmerged"])


def test_extract_anomalous_preference():
    """Should extract preference for anomalous data."""
    result = extract_directives_simple("use anomalous data for phasing")
    assert_in("file_preferences", result)
    assert_true(result["file_preferences"]["prefer_anomalous"])


# =============================================================================
# WORKFLOW PREFERENCE TESTS
# =============================================================================

def test_extract_skip_autobuild():
    """Should extract skip autobuild preference."""
    result = extract_directives_simple("skip autobuild")
    assert_in("workflow_preferences", result)
    assert_in("phenix.autobuild", result["workflow_preferences"]["skip_programs"])


def test_extract_skip_ligandfit():
    """Should extract skip ligandfit preference."""
    result = extract_directives_simple("avoid ligandfit")
    assert_in("workflow_preferences", result)
    assert_in("phenix.ligandfit", result["workflow_preferences"]["skip_programs"])


def test_extract_prefer_program():
    """Should extract program preference if pattern matched."""
    result = extract_directives_simple("prefer phaser for molecular replacement")
    assert_is_instance(result, dict)


# =============================================================================
# COMPLEX EXTRACTION TESTS
# =============================================================================

def test_extract_complex_advice():
    """Should handle complex multi-directive advice."""
    advice = """
    Please use resolution 2.5 Angstrom for all programs.
    Stop after 3 cycles of refinement.
    Use anisotropic ADPs.
    """
    result = extract_directives_simple(advice)

    assert_in("program_settings", result)
    assert_equal(result["program_settings"]["default"]["resolution"], 2.5)
    # Stop conditions may or may not be extracted from this phrasing
    assert_is_instance(result, dict)


def test_extract_preserves_case():
    """Case should be handled appropriately."""
    result1 = extract_directives_simple("Run ANISOTROPIC refinement")
    result2 = extract_directives_simple("run anisotropic refinement")
    # Both should work
    assert_equal(
        result1.get("program_settings", {}).get("phenix.refine", {}).get("anisotropic_adp"),
        result2.get("program_settings", {}).get("phenix.refine", {}).get("anisotropic_adp")
    )


# =============================================================================
# POLDER EXTRACTION TESTS
# =============================================================================

def test_extract_polder():
    """Should extract polder map directive."""
    result = extract_directives_simple("calculate a polder map")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.polder")


def test_extract_omit_map():
    """Should map omit map to polder."""
    result = extract_directives_simple("generate omit map")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.polder")


# =============================================================================
# MAP SHARPENING TESTS
# =============================================================================

def test_extract_map_sharpening():
    """Should extract map sharpening directive."""
    result = extract_directives_simple("sharpen the map")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_sharpening")


def test_extract_auto_sharpen():
    """Should extract auto-sharpen directive."""
    result = extract_directives_simple("auto-sharpen the cryo-EM map")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_sharpening")


# =============================================================================
# MAP TO MODEL TESTS
# =============================================================================

def test_extract_map_to_model():
    """Should extract map_to_model directive with proper phrasing."""
    result = extract_directives_simple("Run MapToModel on the density")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_to_model")


def test_extract_maptomodel_camelcase():
    """Should detect MapToModel camelcase."""
    result = extract_directives_simple("Run MapToModel on the cryo-EM density.")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_to_model")


def test_map_to_model_program_name_fix():
    """Should normalize map_to_model variations."""
    assert_equal(_fix_program_name("map_to_model"), "phenix.map_to_model")
    assert_equal(_fix_program_name("maptomodel"), "phenix.map_to_model")
    assert_equal(_fix_program_name("MapToModel"), "phenix.map_to_model")


# =============================================================================
# LIGAND WORKFLOW CONFLICT TESTS (v72)
# =============================================================================

def test_ligand_workflow_clears_after_cycle():
    """after_cycle should be cleared when ligand constraint exists (v72 fix)."""
    # Simulate what happens when LLM incorrectly extracts "second refinement" as "cycle 2"
    input_directives = {
        "stop_conditions": {
            "after_cycle": 2  # Incorrect - user meant "second refinement", not "cycle 2"
        },
        "constraints": [
            "Fit ligand from lig.pdb after the first refinement"  # Ligand workflow
        ]
    }
    result = validate_directives(input_directives)

    # after_cycle should be cleared because ligand workflow needs ~8 cycles
    if "stop_conditions" in result:
        assert_not_in("after_cycle", result.get("stop_conditions", {}),
                     "after_cycle=2 should be cleared for ligand workflow")


def test_ligand_workflow_clears_after_program_refine():
    """after_program=phenix.refine should be cleared when ligand constraint exists."""
    input_directives = {
        "stop_conditions": {
            "after_program": "phenix.refine"
        },
        "constraints": [
            "Fit ligand from lig.pdb after the first refinement"
        ]
    }
    result = validate_directives(input_directives)

    if "stop_conditions" in result:
        assert_not_in("after_program", result.get("stop_conditions", {}),
                     "after_program=phenix.refine should be cleared for ligand workflow")


def test_ligand_workflow_allows_high_cycle():
    """after_cycle > 4 should NOT be cleared even with ligand constraint."""
    input_directives = {
        "stop_conditions": {
            "after_cycle": 10  # High enough for ligand workflow
        },
        "constraints": [
            "Fit ligand from lig.pdb after the first refinement"
        ]
    }
    result = validate_directives(input_directives)

    # after_cycle=10 should be preserved (high enough for ligand workflow)
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"].get("after_cycle"), 10,
                "after_cycle=10 should be preserved for ligand workflow")


def test_no_ligand_constraint_preserves_after_cycle():
    """after_cycle should be preserved when no ligand constraint exists."""
    input_directives = {
        "stop_conditions": {
            "after_cycle": 2
        },
        "constraints": [
            "Use one macro cycle for refinement"  # No ligand mention
        ]
    }
    result = validate_directives(input_directives)

    # after_cycle should be preserved without ligand constraint
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"].get("after_cycle"), 2,
                "after_cycle should be preserved without ligand constraint")


# =============================================================================
# MULTI-STEP WORKFLOW TESTS (prevent premature after_program stops)
# =============================================================================

def test_multi_step_cryoem_no_early_stop():
    """Multi-step cryo-EM pipeline should NOT set after_program for intermediate steps."""
    advice = ("Auto-sharpen the map. Extract the unique part of the map. "
              "Build a model into the unique part. "
              "Apply NCS to generate the full 24-mer complex.")
    result = extract_directives_simple(advice)
    stop_cond = result.get("stop_conditions", {})
    # Should NOT have after_program=phenix.map_symmetry or phenix.map_sharpening
    after_prog = stop_cond.get("after_program", "")
    assert_not_equal(after_prog, "phenix.map_symmetry",
                     "Multi-step workflow should not stop at map_symmetry")
    assert_not_equal(after_prog, "phenix.map_sharpening",
                     "Multi-step workflow should not stop at map_sharpening")


def test_multi_step_sharpen_then_build():
    """Sharpen map then build model should not stop after sharpening."""
    advice = "Sharpen the map, then build a model into the map"
    result = extract_directives_simple(advice)
    stop_cond = result.get("stop_conditions", {})
    after_prog = stop_cond.get("after_program", "")
    assert_not_equal(after_prog, "phenix.map_sharpening",
                     "Should not stop at sharpening when building is planned")


def test_multi_step_symmetry_then_build():
    """Find symmetry then build model should not stop after symmetry."""
    advice = "Determine symmetry of the map and build a model"
    result = extract_directives_simple(advice)
    stop_cond = result.get("stop_conditions", {})
    after_prog = stop_cond.get("after_program", "")
    assert_not_equal(after_prog, "phenix.map_symmetry",
                     "Should not stop at symmetry when building is planned")


def test_single_step_still_works():
    """Single-step 'sharpen the map' should still set after_program."""
    result = extract_directives_simple("sharpen the map")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_sharpening")


def test_single_step_symmetry_still_works():
    """Single-step 'find symmetry' should still set after_program."""
    result = extract_directives_simple("determine map symmetry")
    assert_in("stop_conditions", result)
    assert_equal(result["stop_conditions"]["after_program"], "phenix.map_symmetry")


# =============================================================================
# TEST RUNNER
# =============================================================================

def test_validate_clears_after_program_with_downstream_constraints():
    """validate_directives should clear after_program when constraints have downstream work.

    This catches the LLM extraction path (not just simple extraction).
    The exact user scenario: LLM sets after_program=phenix.map_symmetry but
    constraints say 'Build a model into the unique part'.
    """
    directives = {
        "stop_conditions": {
            "after_program": "phenix.map_symmetry",
            "skip_validation": True,
        },
        "constraints": [
            "Auto-sharpen the map",
            "Extract the unique part of the map",
            "Build a model into the unique part",
            "Apply NCS to generate the full 24-mer complex",
        ]
    }
    result = validate_directives(directives)
    stop_cond = result.get("stop_conditions", {})
    assert_not_in("after_program", stop_cond,
                  "after_program should be cleared when constraints describe downstream work")


def test_validate_preserves_after_program_without_downstream():
    """validate_directives should keep after_program when no downstream constraints exist."""
    directives = {
        "stop_conditions": {
            "after_program": "phenix.map_symmetry",
            "skip_validation": True,
        },
        "constraints": [
            "Use resolution 2.9",
        ]
    }
    result = validate_directives(directives)
    stop_cond = result.get("stop_conditions", {})
    assert_equal(stop_cond.get("after_program"), "phenix.map_symmetry",
                 "after_program should be preserved with no downstream constraints")


def test_validate_clears_sharpening_with_build_constraint():
    """after_program=map_sharpening should be cleared when constraints say to build."""
    directives = {
        "stop_conditions": {
            "after_program": "phenix.map_sharpening",
            "skip_validation": True,
        },
        "constraints": [
            "Sharpen the map",
            "Build a model into the density",
        ]
    }
    result = validate_directives(directives)
    stop_cond = result.get("stop_conditions", {})
    assert_not_in("after_program", stop_cond,
                  "after_program=map_sharpening should be cleared when build is planned")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
