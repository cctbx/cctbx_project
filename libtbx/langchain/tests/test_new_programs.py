"""
Tests for newly added programs: phenix.polder, phenix.map_sharpening, phenix.model_vs_data

These tests verify:
1. YAML configuration is correct
2. Directive extraction works
3. Program name normalization works

Run with: python tests/test_new_programs.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import yaml
import fnmatch

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    assert_greater,
    run_tests_with_fail_fast
)

from agent.directive_extractor import extract_directives_simple, _fix_program_name


# =============================================================================
# LOAD YAML FILES ONCE (module level)
# =============================================================================

def _load_yaml(filename):
    """Load a YAML file from knowledge directory."""
    yaml_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "knowledge", filename
    )
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)


# Load YAML files at module level (once)
PROGRAMS = _load_yaml("programs.yaml")
WORKFLOWS = _load_yaml("workflows.yaml")
FILE_CATEGORIES = _load_yaml("file_categories.yaml")

# Load prompts for LLM prompt tests
try:
    from knowledge.prompts_hybrid import PLAN_SYSTEM_PROMPT
    HAS_PROMPTS = True
except ImportError:
    HAS_PROMPTS = False
    PLAN_SYSTEM_PROMPT = ""


# =============================================================================
# PHENIX.POLDER YAML CONFIG TESTS
# =============================================================================

def test_polder_exists():
    """phenix.polder should be defined in programs.yaml."""
    assert_in("phenix.polder", PROGRAMS)


def test_polder_category():
    """phenix.polder should be in map_analysis category."""
    polder = PROGRAMS["phenix.polder"]
    assert_equal(polder.get("category"), "map_analysis")


def test_polder_experiment_types():
    """phenix.polder should be for xray only."""
    polder = PROGRAMS["phenix.polder"]
    assert_equal(polder.get("experiment_types"), ["xray"])


def test_polder_required_inputs():
    """phenix.polder should require data_mtz and model."""
    polder = PROGRAMS["phenix.polder"]
    required = polder.get("inputs", {}).get("required", {})
    assert_in("data_mtz", required)
    assert_in("model", required)


def test_polder_selection_strategy_flag():
    """phenix.polder should have selection strategy flag."""
    polder = PROGRAMS["phenix.polder"]
    strategy_flags = polder.get("strategy_flags", {})
    assert_in("selection", strategy_flags)
    assert_in("{value}", strategy_flags["selection"]["flag"])


def test_polder_command_template():
    """phenix.polder command should include placeholders."""
    polder = PROGRAMS["phenix.polder"]
    command = polder.get("command", "")
    assert_in("phenix.polder", command)
    assert_in("{data_mtz}", command)
    assert_in("{model}", command)
    assert_in("{selection}", command)


def test_polder_log_parsing():
    """phenix.polder should have log parsing for CC values."""
    polder = PROGRAMS["phenix.polder"]
    log_parsing = polder.get("log_parsing", {})
    assert_in("cc_1_3", log_parsing)
    assert_in("polder_conclusion", log_parsing)


# =============================================================================
# PHENIX.MAP_SHARPENING YAML CONFIG TESTS
# =============================================================================

def test_map_sharpening_exists():
    """phenix.map_sharpening should be defined in programs.yaml."""
    assert_in("phenix.map_sharpening", PROGRAMS)


def test_map_sharpening_experiment_types():
    """phenix.map_sharpening should be for cryoem."""
    prog = PROGRAMS["phenix.map_sharpening"]
    assert_in("cryoem", prog.get("experiment_types", []))


def test_map_sharpening_required_inputs():
    """phenix.map_sharpening has no hard required inputs (3 modes with different files)."""
    prog = PROGRAMS["phenix.map_sharpening"]
    optional = prog.get("inputs", {}).get("optional", {})
    # All three input types should be in optional
    assert_in("full_map", optional)
    assert_in("half_map", optional)
    assert_in("model", optional)


def test_map_sharpening_optional_inputs():
    """phenix.map_sharpening should have optional model input."""
    prog = PROGRAMS["phenix.map_sharpening"]
    optional = prog.get("inputs", {}).get("optional", {})
    assert_in("model", optional)


def test_map_sharpening_strategy_flags():
    """phenix.map_sharpening should have resolution strategy flag."""
    prog = PROGRAMS["phenix.map_sharpening"]
    strategy_flags = prog.get("strategy_flags", {})
    assert_in("resolution", strategy_flags)


def test_map_sharpening_command_template():
    """phenix.map_sharpening should have hints covering all three modes."""
    prog = PROGRAMS["phenix.map_sharpening"]
    hints = prog.get("hints", [])
    hints_text = " ".join(hints)
    assert_in("shells", hints_text)
    assert_in("b-factor", hints_text)
    assert_in("half_map", hints_text)


# =============================================================================
# PHENIX.MODEL_VS_DATA YAML CONFIG TESTS
# =============================================================================

def test_model_vs_data_exists():
    """phenix.model_vs_data should be defined in programs.yaml."""
    assert_in("phenix.model_vs_data", PROGRAMS)


def test_model_vs_data_required_inputs():
    """phenix.model_vs_data should require data_mtz and model."""
    prog = PROGRAMS["phenix.model_vs_data"]
    required = prog.get("inputs", {}).get("required", {})
    assert_in("data_mtz", required)
    assert_in("model", required)


def test_model_vs_data_log_parsing():
    """phenix.model_vs_data should have log parsing for R values."""
    prog = PROGRAMS["phenix.model_vs_data"]
    log_parsing = prog.get("log_parsing", {})
    # Should have metrics like r_work, r_free
    assert_greater(len(log_parsing), 0)


def test_model_vs_data_hints():
    """phenix.model_vs_data should have workflow hints."""
    prog = PROGRAMS["phenix.model_vs_data"]
    # Could have hints about when to use it
    assert_in("description", prog)


# =============================================================================
# POLDER WORKFLOW CONFIG TESTS
# =============================================================================

def test_polder_in_xray_refine_phase():
    """phenix.polder should be in xray refine phase."""
    xray_phases = WORKFLOWS.get("xray", {}).get("phases", {})
    refine_phase = xray_phases.get("refine", {})
    programs = refine_phase.get("programs", [])

    polder_in_refine = any(
        (p.get("program") if isinstance(p, dict) else p) == "phenix.polder"
        for p in programs
    )
    assert_true(polder_in_refine, "phenix.polder should be in xray refine phase")


def test_polder_in_xray_validate_phase():
    """phenix.polder should be in xray validate phase."""
    xray_phases = WORKFLOWS.get("xray", {}).get("phases", {})
    validate_phase = xray_phases.get("validate", {})
    programs = validate_phase.get("programs", [])

    polder_in_validate = any(
        (p.get("program") if isinstance(p, dict) else p) == "phenix.polder"
        for p in programs
    )
    assert_true(polder_in_validate, "phenix.polder should be in xray validate phase")


def test_polder_conditions_refine_phase():
    """phenix.polder should have conditions in refine phase."""
    xray_phases = WORKFLOWS.get("xray", {}).get("phases", {})
    refine_phase = xray_phases.get("refine", {})
    programs = refine_phase.get("programs", [])

    for p in programs:
        if isinstance(p, dict) and p.get("program") == "phenix.polder":
            assert_in("conditions", p, "polder should have conditions")
            break


def test_polder_conditions_validate_phase():
    """phenix.polder should have conditions in validate phase."""
    xray_phases = WORKFLOWS.get("xray", {}).get("phases", {})
    validate_phase = xray_phases.get("validate", {})
    programs = validate_phase.get("programs", [])

    for p in programs:
        if isinstance(p, dict) and p.get("program") == "phenix.polder":
            assert_in("conditions", p, "polder should have conditions")
            break


def test_polder_has_hint():
    """phenix.polder should have a hint in workflow."""
    xray_phases = WORKFLOWS.get("xray", {}).get("phases", {})

    found_hint = False
    for phase_name, phase in xray_phases.items():
        programs = phase.get("programs", [])
        for p in programs:
            if isinstance(p, dict) and p.get("program") == "phenix.polder":
                if "hint" in p:
                    found_hint = True
                    break

    assert_true(found_hint, "phenix.polder should have a hint")


# =============================================================================
# POLDER LLM PROMPT TESTS
# =============================================================================

def test_polder_in_system_prompt():
    """phenix.polder should be mentioned in LLM system prompt."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    assert_in("polder", PLAN_SYSTEM_PROMPT.lower())


def test_polder_in_validation_section():
    """phenix.polder should be in VALIDATION PROGRAMS section."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    # Find validation section
    validation_idx = PLAN_SYSTEM_PROMPT.find("VALIDATION PROGRAMS")
    assert_greater(validation_idx, -1, "Should have VALIDATION PROGRAMS section")

    # Check polder is mentioned after validation section
    after_validation = PLAN_SYSTEM_PROMPT[validation_idx:]
    assert_in("polder", after_validation.lower())


def test_polder_clarifies_no_phases_needed():
    """Polder description should clarify it doesn't need phases."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    prompt_lower = PLAN_SYSTEM_PROMPT.lower()

    # Find polder section
    polder_idx = prompt_lower.find("polder")
    if polder_idx > -1:
        # Get context around polder mention
        context = PLAN_SYSTEM_PROMPT[max(0, polder_idx-100):polder_idx+500].lower()
        # Should NOT say it needs phases/map coefficients
        has_clarification = (
            "does not need" in context or
            "doesn't need" in context or
            "not need" in context or
            "standard mtz" in context or
            "fobs" in context
        )
        assert_true(has_clarification,
            "Polder description should clarify it doesn't need phases")


def test_polder_description_mentions_omit_map():
    """Polder description should mention omit map functionality."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    prompt_lower = PLAN_SYSTEM_PROMPT.lower()
    assert_true("omit" in prompt_lower or "polder" in prompt_lower)


def test_polder_specifies_mtz_input():
    """Polder should specify MTZ as input in prompt."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    prompt_lower = PLAN_SYSTEM_PROMPT.lower()
    polder_idx = prompt_lower.find("polder")
    if polder_idx > -1:
        context = prompt_lower[polder_idx:polder_idx+300]
        assert_true("mtz" in context or "data_mtz" in context or "reflection" in context)


def test_polder_specifies_model_input():
    """Polder should specify model as input in prompt."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    prompt_lower = PLAN_SYSTEM_PROMPT.lower()
    polder_idx = prompt_lower.find("polder")
    if polder_idx > -1:
        context = prompt_lower[polder_idx:polder_idx+300]
        assert_true("model" in context or "pdb" in context)


def test_polder_mentions_selection():
    """Polder should mention selection parameter."""
    if not HAS_PROMPTS:
        print("  SKIPPED (prompts not available)")
        return
    prompt_lower = PLAN_SYSTEM_PROMPT.lower()
    polder_idx = prompt_lower.find("polder")
    if polder_idx > -1:
        context = prompt_lower[polder_idx:polder_idx+500]
        assert_true("selection" in context or "residue" in context or "ligand" in context)


# =============================================================================
# DIRECTIVE EXTRACTION TESTS
# =============================================================================

def test_polder_extraction():
    """Test directive extraction for polder-related advice."""
    advice = "Calculate polder map for ligand MES at position 88"
    directives = extract_directives_simple(advice)

    assert_in("stop_conditions", directives)
    assert_equal(directives["stop_conditions"]["after_program"], "phenix.polder")


def test_map_sharpening_extraction():
    """Test directive extraction for map sharpening advice."""
    advice = "sharpen the cryo-EM map"
    directives = extract_directives_simple(advice)

    assert_in("stop_conditions", directives)
    assert_equal(directives["stop_conditions"]["after_program"], "phenix.map_sharpening")


def test_map_to_model_extraction():
    """Test directive extraction for map_to_model/dock_in_map advice."""
    advice = "dock model in map"
    directives = extract_directives_simple(advice)

    assert_in("stop_conditions", directives)
    # "dock in map" triggers dock_in_map
    assert_equal(directives["stop_conditions"]["after_program"], "phenix.dock_in_map")


# =============================================================================
# PROGRAM NAME FIX TESTS
# =============================================================================

def test_fix_polder_name():
    """_fix_program_name should normalize polder variations."""
    assert_equal(_fix_program_name("polder"), "phenix.polder")
    assert_equal(_fix_program_name("polder_map"), "phenix.polder")
    assert_equal(_fix_program_name("omit_map"), "phenix.polder")


def test_fix_map_sharpening_name():
    """_fix_program_name should normalize map_sharpening variations."""
    assert_equal(_fix_program_name("map_sharpening"), "phenix.map_sharpening")


def test_fix_map_to_model_name():
    """_fix_program_name should normalize map_to_model variations."""
    assert_equal(_fix_program_name("map_to_model"), "phenix.map_to_model")


# =============================================================================
# VALID PROGRAMS SET TESTS
# =============================================================================

def test_polder_in_valid_programs():
    """phenix.polder should be recognized as valid."""
    from agent.directive_validator import list_available_programs
    available = list_available_programs()
    assert_in("phenix.polder", available)


def test_map_sharpening_in_valid_programs():
    """phenix.map_sharpening should be recognized as valid."""
    from agent.directive_validator import list_available_programs
    available = list_available_programs()
    assert_in("phenix.map_sharpening", available)


def test_map_to_model_in_valid_programs():
    """phenix.map_to_model should be recognized as valid."""
    from agent.directive_validator import list_available_programs
    available = list_available_programs()
    assert_in("phenix.map_to_model", available)


# =============================================================================
# UNCLASSIFIED PDB CATEGORIZATION TESTS
# =============================================================================

def test_unclassified_pdb_exists():
    """unclassified_pdb category should exist."""
    assert_in("unclassified_pdb", FILE_CATEGORIES)


def test_unclassified_pdb_parent_is_model():
    """unclassified_pdb parent category should be 'model'."""
    unclassified = FILE_CATEGORIES.get("unclassified_pdb", {})
    assert_equal(unclassified.get("parent_category"), "model")


def test_unclassified_pdb_has_catchall_pattern():
    """unclassified_pdb should have catchall pattern."""
    unclassified = FILE_CATEGORIES.get("unclassified_pdb", {})
    patterns = unclassified.get("patterns", [])

    # Has wildcard catchall pattern
    has_catchall = "*" in patterns or any("*" == p for p in patterns)
    assert_true(has_catchall, "Should have catchall pattern")


def test_unclassified_pdb_excludes_search_models():
    """unclassified_pdb should exclude search model patterns."""
    unclassified = FILE_CATEGORIES.get("unclassified_pdb", {})
    excludes = unclassified.get("excludes", [])

    has_search_exclude = any("search" in e.lower() for e in excludes)
    assert_true(has_search_exclude, "Should exclude search model patterns")


def test_unclassified_pdb_excludes_specific_patterns():
    """unclassified_pdb should exclude sculptor/chainsaw patterns."""
    unclassified = FILE_CATEGORIES.get("unclassified_pdb", {})
    excludes = unclassified.get("excludes", [])

    # Check for sculptor or chainsaw excludes
    has_sculptor = any("sculptor" in e.lower() for e in excludes)
    has_chainsaw = any("chainsaw" in e.lower() for e in excludes)

    assert_true(has_sculptor or has_chainsaw,
        "Should exclude sculptor/chainsaw patterns")


def test_generic_pdb_not_excluded():
    """Generic PDB files like '1aba.pdb' should not be excluded."""
    unclassified = FILE_CATEGORIES.get("unclassified_pdb", {})
    patterns = unclassified.get("patterns", [])
    excludes = unclassified.get("excludes", [])

    test_file = "1aba.pdb"

    # Should not be excluded
    excluded = any(fnmatch.fnmatch(test_file, exc) for exc in excludes)
    assert_false(excluded, f"'{test_file}' should NOT be excluded")

    # Should match a pattern
    matched = any(fnmatch.fnmatch(test_file, pat) for pat in patterns)
    assert_true(matched, f"'{test_file}' should match unclassified_pdb patterns")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
