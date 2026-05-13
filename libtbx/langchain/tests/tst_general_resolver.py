"""
Unit tests for the general after_program resolver (v115.10).

Tests _detect_actions, _resolve_after_program, and end-to-end
behavior through extract_directives_simple and
_apply_workflow_intent_fallback.
"""

from __future__ import absolute_import, division, print_function
import sys
import os

# Setup import paths
sys.path.insert(0, os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))

try:
    from libtbx.langchain.agent.directive_extractor import (
        _detect_actions,
        _resolve_after_program,
        _apply_workflow_intent_fallback,
        extract_directives_simple,
        _ACTION_TABLE,
    )
except ImportError:
    from agent.directive_extractor import (
        _detect_actions,
        _resolve_after_program,
        _apply_workflow_intent_fallback,
        extract_directives_simple,
        _ACTION_TABLE,
    )


# =====================================================================
# Helper
# =====================================================================

def _actions(advice):
    """Shortcut: return action names from advice."""
    return [a for a, _ in _detect_actions(advice.lower())]


def _resolve(advice, initial_ap=None):
    """Shortcut: run resolver, return after_program."""
    d = {}
    if initial_ap:
        d["stop_conditions"] = {"after_program": initial_ap}
    _resolve_after_program(d, advice.lower())
    return d.get("stop_conditions", {}).get("after_program")


# =====================================================================
# TEST 1: ACTION_TABLE structure
# =====================================================================

def test_action_table_structure():
    """Every action has xray, cryoem, and keywords fields."""
    print("Test: action_table_structure")
    for action, info in _ACTION_TABLE.items():
        assert "xray" in info, \
            "%s missing 'xray'" % action
        assert "cryoem" in info, \
            "%s missing 'cryoem'" % action
        assert "keywords" in info, \
            "%s missing 'keywords'" % action
        assert isinstance(info["keywords"], list), \
            "%s keywords not a list" % action
        assert len(info["keywords"]) > 0, \
            "%s has no keywords" % action
        # At least one program must be non-None
        assert info["xray"] or info["cryoem"], \
            "%s has no programs" % action
    print("  PASS")


# =====================================================================
# TEST 2: _detect_actions — single actions
# =====================================================================

def test_detect_single_actions():
    """Each action is detected by its keywords."""
    print("Test: detect_single_actions")
    cases = [
        ("run xtriage", ["analyze"]),
        ("run mtriage", ["analyze"]),
        ("run phaser", ["solve"]),
        ("molecular replacement", ["solve"]),
        ("sad phasing", ["phase"]),
        ("predict and build", ["predict_build"]),
        ("predict a model", ["predict"]),
        ("build a model", ["build"]),
        ("autobuild", ["build"]),
        ("density modification", ["denmod"]),
        ("sharpen the map", ["sharpen"]),
        ("refine the model", ["refine"]),
        ("fit ligand", ["ligandfit"]),
        ("validate the structure", ["validate"]),
        ("polder map", ["polder"]),
        ("dock in map", ["dock"]),
        ("find symmetry", ["symmetry"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 3: _detect_actions — multiple actions
# =====================================================================

def test_detect_multiple_actions():
    """Multiple distinct actions detected in order."""
    print("Test: detect_multiple_actions")
    cases = [
        ("density modify and build a model",
         ["denmod", "build"]),
        ("run phaser and refine",
         ["solve", "refine"]),
        ("fit ligand and refine",
         ["ligandfit", "refine"]),
        ("sharpen, density modify, build",
         ["sharpen", "denmod", "build"]),
        ("refine and validate",
         ["refine", "validate"]),
        ("run phaser then refine",
         ["solve", "refine"]),
        ("density modify and fit the ligand",
         ["denmod", "ligandfit"]),
        ("sharpen and apply symmetry",
         ["sharpen", "symmetry"]),
        ("improve the map and build a model",
         ["denmod", "build"]),
        ("run phaser and pick waters",
         ["solve", "refine"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 4: _detect_actions — predict/build compound rule
# =====================================================================

def test_predict_build_compound():
    """predict + build merge into predict_build."""
    print("Test: predict_build_compound")
    cases = [
        ("predict and build a model", ["predict_build"]),
        ("predict, build, refine", ["predict_build", "refine"]),
        ("predict a model", ["predict"]),
        ("build a model", ["build"]),
        # Direct keyword
        ("predict_and_build", ["predict_build"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 5: _detect_actions — negation detection
# =====================================================================

def test_negation_detection():
    """Negated actions are excluded."""
    print("Test: negation_detection")
    cases = [
        ("don't build a model", []),
        ("do density modification but don't build a model",
         ["denmod"]),
        ("refine but do not validate", ["refine"]),
        ("run phaser, no refinement", ["solve"]),
        ("skip validation, just refine", ["refine"]),
        ("never refine, just validate", ["validate"]),
        # Non-negated control
        ("build a model", ["build"]),
        ("refine and validate", ["refine", "validate"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 6: _detect_actions — word boundary protection
# =====================================================================

def test_word_boundaries():
    """Single-word keywords don't match inside other words."""
    print("Test: word_boundaries")
    cases = [
        # "solve" should NOT match inside "resolve"
        ("resolve_cryo_em and build",
         ["denmod", "build"]),
        # "predict" should NOT match inside "unpredictable"
        ("unpredictable results need refine",
         ["refine"]),
        # "build" should NOT match inside "rebuild" (word boundary)
        # But "rebuild" may match via multi-word "build" keyword
        # because re.search(r'\bbuild\b', 'rebuild') → no match
        # Note: "rebuild the model" → no "build" match
        # Hyphenated keywords use substring matching
        ("create a density-modified map", ["denmod"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 7: _detect_actions — no false positives
# =====================================================================

def test_no_false_positives():
    """Unrelated text produces no actions."""
    print("Test: no_false_positives")
    cases = [
        ("hello world", []),
        ("the crystal has P21 symmetry", []),
        ("upload my files", []),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 8: _resolve_after_program — multi-step clearing
# =====================================================================

def test_resolve_multi_step_clear():
    """Multiple actions without stop → clear after_program."""
    print("Test: resolve_multi_step_clear")
    cases = [
        ("density modify and build a model",
         "phenix.resolve_cryo_em", None),
        ("run phaser and refine",
         "phenix.phaser", None),
        ("fit ligand and refine",
         "phenix.ligandfit", None),
        ("sharpen and density modify",
         "phenix.map_sharpening", None),
        ("refine and validate",
         "phenix.refine", None),
    ]
    for advice, initial, expected in cases:
        got = _resolve(advice, initial)
        assert got == expected, \
            "'%s' (init=%s) → %s, expected %s" % (
                advice, initial, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 9: _resolve_after_program — multi-step + stop → last
# =====================================================================

def test_resolve_multi_step_stop_last():
    """Multiple actions + stop → after_program = last."""
    print("Test: resolve_multi_step_stop_last")
    cases = [
        ("run phaser, refine, and stop",
         "phenix.phaser", "phenix.refine"),
        ("sharpen, density modify, build, and stop",
         "phenix.map_sharpening", "phenix.predict_and_build"),
        ("predict, build, refine, and stop",
         None, "phenix.refine"),
    ]
    for advice, initial, expected in cases:
        got = _resolve(advice, initial)
        assert got == expected, \
            "'%s' (init=%s) → %s, expected %s" % (
                advice, initial, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 10: _resolve_after_program — single + stop → set
# =====================================================================

def test_resolve_single_stop_set():
    """Single action + stop → set after_program (overriding LLM)."""
    print("Test: resolve_single_stop_set")
    cases = [
        # LLM missed it
        ("density modify and stop", None,
         "phenix.autobuild_denmod"),  # xray default
        ("create a density-modified map and stop", None,
         "phenix.resolve_cryo_em"),  # cryo indicator
        ("run xtriage and stop", None,
         "phenix.xtriage"),
        ("fit ligand and stop", None,
         "phenix.ligandfit"),
        ("run mtriage and stop", None,
         "phenix.mtriage"),
        # LLM set WRONG program — resolver overrides
        ("create a density-modified map and stop",
         "phenix.real_space_refine",
         "phenix.resolve_cryo_em"),
        ("run xtriage and stop",
         "phenix.refine",
         "phenix.xtriage"),
    ]
    for advice, initial, expected in cases:
        got = _resolve(advice, initial)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 11: _resolve_after_program — single, no stop → leave as-is
# =====================================================================

def test_resolve_single_no_stop_unchanged():
    """Single action, no stop → leave after_program as-is."""
    print("Test: resolve_single_no_stop_unchanged")
    # Already set → stays
    got = _resolve("solve the structure", "phenix.phaser")
    assert got == "phenix.phaser", \
        "Expected phenix.phaser, got %s" % got
    # Not set → stays None
    got2 = _resolve("solve the structure", None)
    assert got2 is None, \
        "Expected None, got %s" % got2
    print("  PASS")


# =====================================================================
# TEST 12: _resolve_after_program — stop negation
# =====================================================================

def test_resolve_stop_negation():
    """'don't stop' does not trigger stop logic."""
    print("Test: resolve_stop_negation")
    got = _resolve("don't stop after refinement", "phenix.refine")
    assert got == "phenix.refine", \
        "Expected phenix.refine (unchanged), got %s" % got
    print("  PASS")


# =====================================================================
# TEST 13: _resolve_after_program — experiment type inference
# =====================================================================

def test_resolve_experiment_type():
    """Cryo-EM indicators produce cryo-EM programs."""
    print("Test: resolve_experiment_type")
    # cryo-EM indicators
    cryo_cases = [
        ("sharpen and stop", "phenix.map_sharpening"),
        ("create a density-modified map and stop",
         "phenix.resolve_cryo_em"),
        ("cryo-em density modify and stop",
         "phenix.resolve_cryo_em"),
    ]
    for advice, expected in cryo_cases:
        got = _resolve(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    # xray default
    xray_cases = [
        ("density modify and stop",
         "phenix.autobuild_denmod"),
        ("refine and stop",
         "phenix.refine"),
    ]
    for advice, expected in xray_cases:
        got = _resolve(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % (
        len(cryo_cases) + len(xray_cases)))


# =====================================================================
# TEST 14: stop_after_predict
# =====================================================================

def test_stop_after_predict():
    """'predict and stop' sets stop_after_predict=True."""
    print("Test: stop_after_predict")
    d = {}
    _resolve_after_program(d, "predict a model and stop")
    _pab = "phenix.predict_and_build"
    sap = d.get("program_settings", {}).get(
        _pab, {}).get("stop_after_predict", False)
    ap = d.get("stop_conditions", {}).get("after_program")
    assert sap is True, \
        "stop_after_predict should be True, got %s" % sap
    assert ap == _pab, \
        "after_program should be predict_and_build, got %s" % ap
    # predict_build should NOT set stop_after_predict
    d2 = {}
    _resolve_after_program(d2, "predict and build and stop")
    sap2 = d2.get("program_settings", {}).get(
        _pab, {}).get("stop_after_predict", False)
    assert sap2 is False, \
        "predict_build should NOT set stop_after_predict"
    print("  PASS")


# =====================================================================
# TEST 15: start_with_program set for multi-step
# =====================================================================

def test_start_with_program():
    """Multi-step sets start_with_program to first action."""
    print("Test: start_with_program")
    d = {}
    _resolve_after_program(d, "run phaser and refine")
    swp = d.get("stop_conditions", {}).get(
        "start_with_program")
    assert swp == "phenix.phaser", \
        "start_with_program should be phenix.phaser, got %s" % swp
    # With stop: start_with still set
    d2 = {}
    _resolve_after_program(d2, "run phaser, refine, and stop")
    swp2 = d2.get("stop_conditions", {}).get(
        "start_with_program")
    assert swp2 == "phenix.phaser", \
        "start_with_program should be phenix.phaser, got %s" % swp2
    print("  PASS")


# =====================================================================
# TEST 16: Preserved overlays (model_is_placed, MR-SAD, validation)
# =====================================================================

def test_preserved_overlays():
    """model_is_placed, MR-SAD, and validation-only still work."""
    print("Test: preserved_overlays")
    # model_is_placed
    d = {}
    _apply_workflow_intent_fallback(d, "fit atp ligand")
    assert d.get("workflow_preferences", {}).get(
        "model_is_placed") is True, \
        "model_is_placed should be True for 'fit atp'"
    # MR-SAD
    d2 = {}
    _apply_workflow_intent_fallback(d2, "mr-sad phasing")
    assert d2.get("workflow_preferences", {}).get(
        "use_mr_sad") is True, \
        "use_mr_sad should be True for 'mr-sad'"
    # Validation-only
    d3 = {}
    _apply_workflow_intent_fallback(d3, "validate this structure")
    assert d3.get("workflow_preferences", {}).get(
        "wants_validation_only") is True, \
        "wants_validation_only should be True"
    # model_is_placed NOT set without ligand signals
    d4 = {}
    _apply_workflow_intent_fallback(d4, "refine the model")
    assert not d4.get("workflow_preferences", {}).get(
        "model_is_placed", False), \
        "model_is_placed should NOT be set for 'refine'"
    print("  PASS")


# =====================================================================
# TEST 17: End-to-end through extract_directives_simple
# =====================================================================

def test_end_to_end():
    """Full rules path produces correct after_program."""
    print("Test: end_to_end")
    cases = [
        ("density modify the map and build an initial model",
         None, "multi-step cleared"),
        ("create a density-modified map and stop",
         "phenix.resolve_cryo_em", "single + stop"),
        ("run phaser and refine",
         None, "multi-step cleared"),
        ("run phaser, refine, and stop",
         "phenix.refine", "last + stop"),
        ("fit ligand (ATP) and stop",
         "phenix.ligandfit", "ligand + stop"),
        ("run xtriage and stop",
         "phenix.xtriage", "analyze + stop"),
        ("run mtriage and stop",
         "phenix.mtriage", "mtriage + stop"),
    ]
    for advice, expected_ap, desc in cases:
        d = extract_directives_simple(advice)
        got = d.get("stop_conditions", {}).get("after_program")
        assert got == expected_ap, \
            "'%s' (%s) → %s, expected %s" % (
                advice, desc, got, expected_ap)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 18: Formerly handled by continuation_indicators
# =====================================================================

def test_replaces_continuation_indicators():
    """Scenarios formerly caught by continuation_indicators."""
    print("Test: replaces_continuation_indicators")
    cases = [
        # "then" between programs
        ("run phaser then refine", ["solve", "refine"]),
        # "followed by"
        ("phaser followed by refinement", ["solve", "refine"]),
        # "and then"
        ("fit ligand and then refine", ["ligandfit", "refine"]),
        # "afterwards"
        ("run mtriage, afterwards density modify",
         ["analyze", "denmod"]),
        # Edge: "then" + stop → should still set after_program
        ("run xtriage then stop", ["analyze"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    # Edge case: "run xtriage then stop" should SET after_program
    # (old system wrongly treated "then" as continuation)
    got_ap = _resolve("run xtriage then stop")
    assert got_ap == "phenix.xtriage", \
        "xtriage then stop → %s (expected xtriage)" % got_ap
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 19: Formerly handled by downstream_tasks
# =====================================================================

def test_replaces_downstream_tasks():
    """Scenarios formerly caught by downstream_tasks."""
    print("Test: replaces_downstream_tasks")
    cases = [
        ("refine and place the ligand",
         ["refine", "ligandfit"]),
        ("density modify and build a model",
         ["denmod", "build"]),
        ("sharpen then map to model",
         ["sharpen", "build"]),
        # "apply symmetry" → 1 action (correct — 1 program)
        ("find symmetry and apply it",
         ["symmetry"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 20: Formerly handled by multi_program_patterns
# =====================================================================

def test_replaces_multi_program_patterns():
    """Scenarios formerly caught by multi_program_patterns."""
    print("Test: replaces_multi_program_patterns")
    cases = [
        ("sharpen and density modify",
         ["sharpen", "denmod"]),
        ("dock in map and refine",
         ["dock", "refine"]),
        ("fit the ligand and refine",
         ["ligandfit", "refine"]),
        ("sharpen and apply symmetry",
         ["sharpen", "symmetry"]),
    ]
    for advice, expected in cases:
        got = _actions(advice)
        assert got == expected, \
            "'%s' → %s, expected %s" % (advice, got, expected)
    print("  PASS (%d cases)" % len(cases))


# =====================================================================
# TEST 21: Preprocessed text — negated actions in Special Instructions
# =====================================================================

def test_preprocessed_text_negation():
    """Actions in 'Do not proceed to X, Y' are not falsely detected."""
    print("Test: preprocessed_text_negation")
    # Exact preprocessed text from production bug
    preprocessed = (
        "1. Experiment Type: cryo-EM (half-maps provided)\n"
        "2. Primary Goal: Perform density modification of the "
        "cryo-EM half-maps and then stop (no further steps).\n"
        "5. Special Instructions:\n"
        "   * Do not proceed to model building, refinement, "
        "or any downstream processing after density modification."
    ).lower()
    got = _resolve(preprocessed)
    assert got == "phenix.resolve_cryo_em", \
        "Expected resolve_cryo_em, got %s" % got
    # Multi-step preprocessed text — build NOT negated
    preprocessed2 = (
        "2. Primary Goal: Density modify the map and build "
        "an initial model.\n"
        "5. Special Instructions: None"
    ).lower()
    got2 = _resolve(preprocessed2)
    assert got2 is None, \
        "Expected None (multi-step cleared), got %s" % got2
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """Run all tests."""
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    # Simple runner when tst_utils not available
    test_fns = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)
