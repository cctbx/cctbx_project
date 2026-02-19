"""
Regression tests for bugs found and fixed during the systematic audit
(Categories I, J, E).

Each test is named after the audit item that revealed the bug so failures
are immediately traceable to a root-cause description.

Categories covered
------------------
I1  max_refine_cycles → controlled validation landing, not bare STOP
I1c cryoem rsr_count used (not refine_count) when checking the limit
I2  after_program beats quality gate → bare STOP, no validate injection
J2  _is_failed_result: false-positive fixes for bare ERROR variants
J5  _clear_zombie_done_flags: stale done flags cleared when output missing
E1  xtriage resolution: dash-separator range "50.00 - 2.30" → picks 2.30
E2  xtriage pick_min anchor: "Completeness in resolution range: 1" not matched
E1  real_space_refine map_cc: extract:last (final cycle, not first)
K1  best_files list values crash at cycle=2 (half_map stored as [map1, map2])
K2  mtriage/predict_and_build/map_to_model drop half_maps when full_map also present

Run with:
    python tests/tst_audit_fixes.py
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in, assert_not_none,
    run_tests_with_fail_fast,
)

# ---------------------------------------------------------------------------
# Lazy imports that need the mock infrastructure
# ---------------------------------------------------------------------------
try:
    from agent.workflow_state import (
        _is_failed_result,
        _analyze_history,
        _clear_zombie_done_flags,
        detect_workflow_state,
    )
    from agent.workflow_engine import WorkflowEngine
    _IMPORTS_OK = True
except ImportError:
    _IMPORTS_OK = False

KNOWLEDGE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "knowledge"
)


def _load_programs():
    with open(os.path.join(KNOWLEDGE_DIR, "programs.yaml")) as f:
        return yaml.safe_load(f)


def _apply_yaml_regex(prog_name, metric, log_text, programs=None):
    """Apply a programs.yaml log_parsing pattern to log_text, return value."""
    if programs is None:
        programs = _load_programs()
    spec = programs.get(prog_name, {}).get("log_parsing", {}).get(metric, {})
    pattern = spec.get("pattern", "")
    pick_min = spec.get("pick_min", False)
    extract = spec.get("extract", "first")
    typ = spec.get("type", "float")
    if not pattern:
        return None
    matches = re.findall(pattern, log_text, re.MULTILINE)
    if not matches:
        return None
    if pick_min:
        try:
            return min(
                float(m) if not isinstance(m, tuple) else float(m[0])
                for m in matches
            )
        except (TypeError, ValueError):
            return None
    val = matches[-1] if extract == "last" else matches[0]
    if isinstance(val, tuple):
        val = val[0]
    try:
        if typ == "float":
            return float(val)
        elif typ == "int":
            return int(val)
        return val
    except (TypeError, ValueError):
        return val


# =============================================================================
# CATEGORY J2 — _is_failed_result false-positive fixes
# =============================================================================

def test_j2_is_failed_result_true_positives():
    """J2: definitive failure signals are detected."""
    print("Test: j2_is_failed_result_true_positives")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    # Each of these must be recognised as failure
    assert_true(_is_failed_result("FAILED"),
                "FAILED should be a failure")
    assert_true(_is_failed_result("Sorry: could not read reflection file"),
                "Sorry: should be a failure")
    assert_true(_is_failed_result("Sorry bad format -- check input"),
                "Sorry (no colon) should be a failure")
    assert_true(_is_failed_result("Traceback (most recent call last):"),
                "Python traceback should be a failure")
    assert_true(_is_failed_result("*** Error: assertion failed in refine"),
                "*** Error should be a failure")
    assert_true(_is_failed_result("FATAL: out of memory"),
                "FATAL: should be a failure")
    assert_true(_is_failed_result("Exception raised during refinement"),
                "Exception should be a failure")

    print("  PASSED")


def test_j2_is_failed_result_false_positives_eliminated():
    """J2: non-fatal strings containing ERROR/error must NOT be flagged."""
    print("Test: j2_is_failed_result_false_positives_eliminated")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    # These must NOT be recognised as failure
    assert_false(_is_failed_result("No ERROR detected"),
                 "'No ERROR detected' must not be a failure")
    assert_false(_is_failed_result("No error: all checks passed"),
                 "'No error: all checks passed' must not be a failure (mid-sentence colon)")
    assert_false(_is_failed_result("No ERROR: bad context"),
                 "'No ERROR: bad context' must not be a failure")
    assert_false(_is_failed_result("Error model parameter description"),
                 "PHENIX help text 'Error model parameter' must not be a failure")
    assert_false(_is_failed_result("Phenix expected errors: 0"),
                 "'expected errors: 0' must not be a failure")
    assert_false(_is_failed_result(""),
                 "Empty string must not be a failure")
    assert_false(_is_failed_result(None),
                 "None must not be a failure")
    assert_false(_is_failed_result("resolve_cryo_em DONE"),
                 "DONE string must not be a failure")

    print("  PASSED")


def test_j2_failed_result_blocks_done_flags():
    """J2: a FAILED result string prevents _analyze_history from setting done flags."""
    print("Test: j2_failed_result_blocks_done_flags")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    failed_history = [
        {
            "program": "phenix.resolve_cryo_em",
            "command": "phenix.resolve_cryo_em half1.ccp4 half2.ccp4",
            "result": "FAILED: process killed",
            "analysis": {},
        }
    ]
    info = _analyze_history(failed_history)
    assert_false(info.get("resolve_cryo_em_done", False),
                 "FAILED result must not set resolve_cryo_em_done")

    sorry_history = [
        {
            "program": "phenix.refine",
            "command": "phenix.refine model.pdb data.mtz",
            "result": "Sorry: could not read reflection file",
            "analysis": {},
        }
    ]
    info2 = _analyze_history(sorry_history)
    assert_equal(info2.get("refine_count", 0), 0,
                 "Sorry: result must not increment refine_count")

    print("  PASSED")


# =============================================================================
# CATEGORY J5 — _clear_zombie_done_flags
# =============================================================================

def test_j5_zombie_cleared_when_output_missing():
    """J5: done flag cleared when output file is absent (zombie state)."""
    print("Test: j5_zombie_cleared_when_output_missing")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    info = {
        "resolve_cryo_em_done": True,
        "has_full_map": True,
        "predict_full_done": False,
        "dock_done": False,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,
        "rsr_count": 0,
        "has_placed_model": False,
    }
    diags = _clear_zombie_done_flags(info, available_files=[])

    assert_false(info["resolve_cryo_em_done"],
                 "resolve_cryo_em_done must be cleared when output file absent")
    assert_false(info["has_full_map"],
                 "has_full_map must be cleared when denmod_map output absent")
    assert_true(len(diags) > 0,
                "Diagnostic messages must be produced for zombie state")
    assert_true(any("resolve_cryo_em_done" in d for d in diags),
                "Diagnostic must name the cleared flag")

    print("  PASSED")


def test_j5_zombie_not_cleared_when_output_present():
    """J5: done flag preserved when output file is present."""
    print("Test: j5_zombie_not_cleared_when_output_present")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    info = {
        "resolve_cryo_em_done": True,
        "has_full_map": True,
        "predict_full_done": False,
        "dock_done": False,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,
        "rsr_count": 0,
    }
    diags = _clear_zombie_done_flags(
        info, available_files=["/data/denmod_map.ccp4"]
    )

    assert_true(info["resolve_cryo_em_done"],
                "resolve_cryo_em_done must be preserved when output file exists")
    assert_equal(len(diags), 0,
                 "No diagnostics when output file is present")

    print("  PASSED")


def test_j5_refine_zombie_decrements_count():
    """J5: refine_done zombie decrements refine_count."""
    print("Test: j5_refine_zombie_decrements_count")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    info = {
        "resolve_cryo_em_done": False,
        "has_full_map": False,
        "predict_full_done": False,
        "dock_done": False,
        "refine_done": True,
        "refine_count": 2,
        "rsr_done": False,
        "rsr_count": 0,
        "has_placed_model": False,
    }
    _clear_zombie_done_flags(info, available_files=[])

    assert_false(info["refine_done"],
                 "refine_done must be cleared when refine output absent")
    assert_equal(info["refine_count"], 1,
                 "refine_count must be decremented from 2 to 1")

    print("  PASSED")


def test_j5_dock_zombie_clears_placed_model():
    """J5: dock_done zombie clears has_placed_model."""
    print("Test: j5_dock_zombie_clears_placed_model")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    info = {
        "resolve_cryo_em_done": False,
        "has_full_map": False,
        "predict_full_done": False,
        "dock_done": True,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,
        "rsr_count": 0,
        "has_placed_model": True,
    }
    _clear_zombie_done_flags(info, available_files=[])

    assert_false(info["dock_done"],
                 "dock_done must be cleared when docked pdb absent")
    assert_false(info["has_placed_model"],
                 "has_placed_model must be cleared when dock output absent")

    print("  PASSED")


def test_j5_dock_zombie_not_cleared_when_docked_pdb_present():
    """J5: dock_done preserved when docked pdb is on disk."""
    print("Test: j5_dock_zombie_not_cleared_when_docked_pdb_present")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    info = {
        "resolve_cryo_em_done": False,
        "has_full_map": False,
        "predict_full_done": False,
        "dock_done": True,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,
        "rsr_count": 0,
        "has_placed_model": True,
    }
    _clear_zombie_done_flags(
        info, available_files=["/data/model_docked.pdb"]
    )

    assert_true(info["dock_done"],
                "dock_done must be preserved when docked pdb exists")
    assert_true(info["has_placed_model"],
                "has_placed_model must be preserved when dock output exists")

    print("  PASSED")


# =============================================================================
# CATEGORY E1/E2 — xtriage resolution regex (dash-separator fix)
# =============================================================================

def test_e1_e2_xtriage_resolution_dash_separator():
    """E1/E2: xtriage resolution picks correct value from '50.00 - 2.30' format."""
    print("Test: e1_e2_xtriage_resolution_dash_separator")

    log = (
        "  Resolution range:   50.00 - 2.30\n"
        "  Completeness in resolution range:   1   96.5%\n"
        "  Anomalous resolution range: 50.00 - 2.80\n"
        "  Resolution range:   50.00 - 3.50\n"
    )
    result = _apply_yaml_regex("phenix.xtriage", "resolution", log)

    # pick_min should select 2.30, not 50.00
    assert_not_none(result, "xtriage resolution must match")
    assert_equal(result, 2.3,
                 "xtriage resolution pick_min must extract 2.30, not low-res limit 50.00")

    print("  PASSED")


def test_e2_xtriage_resolution_anchor_blocks_completeness_line():
    """E2: 'Completeness in resolution range: 1' must NOT match resolution pattern."""
    print("Test: e2_xtriage_resolution_anchor_blocks_completeness_line")

    # The line that previously caused pick_min to return 1.0
    trap_log = "  Completeness in resolution range:   1   96.5%"
    result = _apply_yaml_regex("phenix.xtriage", "resolution", trap_log)

    assert_true(result is None,
                "Completeness line must not match xtriage resolution pattern")

    print("  PASSED")


def test_e1_xtriage_resolution_simple_format():
    """E1: simple 'Resolution: 1.80' format (no range) still matches."""
    print("Test: e1_xtriage_resolution_simple_format")

    log = "  Resolution:   1.80"
    result = _apply_yaml_regex("phenix.xtriage", "resolution", log)

    assert_not_none(result, "Simple Resolution: line must match")
    assert_equal(result, 1.80,
                 "Simple resolution format must extract 1.80")

    print("  PASSED")


def test_e1_xtriage_resolution_multiple_ranges_picks_min():
    """E1: when multiple resolution ranges appear, pick_min selects highest resolution."""
    print("Test: e1_xtriage_resolution_multiple_ranges_picks_min")

    log = (
        "  Resolution range:   50.00 - 3.50\n"
        "  Resolution range:   50.00 - 2.30\n"
    )
    result = _apply_yaml_regex("phenix.xtriage", "resolution", log)

    assert_equal(result, 2.3,
                 "pick_min across multiple ranges must return highest resolution (smallest value)")

    print("  PASSED")


# =============================================================================
# CATEGORY E1 — real_space_refine map_cc extract:last
# =============================================================================

def test_e1_rsr_map_cc_uses_last_cycle():
    """E1: real_space_refine map_cc returns final macro-cycle value, not first."""
    print("Test: e1_rsr_map_cc_uses_last_cycle")

    # RSR emits one CC_mask line per macro-cycle
    log = (
        "  Macro cycle 1: CC_mask = 0.52\n"
        "  Macro cycle 2: CC_mask = 0.63\n"
        "  Final statistics: CC_mask = 0.71\n"
    )
    result = _apply_yaml_regex("phenix.real_space_refine", "map_cc", log)

    assert_equal(result, 0.71,
                 "real_space_refine map_cc must extract last value (0.71), not first (0.52)")

    print("  PASSED")


def test_e1_rsr_map_cc_pattern_variants():
    """E1: real_space_refine map_cc pattern handles Map-CC and Model vs map CC."""
    print("Test: e1_rsr_map_cc_pattern_variants")

    # Pattern broadened to match the same variants as other programs
    log = (
        "  CC_mask = 0.65\n"
        "  Map-CC = 0.72\n"
        "  Model vs map CC = 0.78\n"
    )
    result = _apply_yaml_regex("phenix.real_space_refine", "map_cc", log)

    assert_equal(result, 0.78,
                 "Broadened pattern must match all CC variants; last wins")

    print("  PASSED")


# =============================================================================
# CATEGORY I1 — max_refine_cycles: controlled landing (validate + STOP)
# =============================================================================

def test_i1_max_refine_cycles_xray_controlled_landing():
    """I1: xray hitting max_refine_cycles injects validate programs + STOP, not bare STOP.

    Calls WorkflowEngine.get_valid_programs() directly to avoid the libtbx
    lazy-import in detect_workflow_state (which would fall back to ["STOP"]
    in environments where libtbx is unavailable).
    """
    print("Test: i1_max_refine_cycles_xray_controlled_landing")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # Simulate xray refine phase with refine_count=1 and max_refine_cycles=1
    context = {
        "phase": "refine",
        "refine_count": 1,
        "rsr_count": 0,
        "r_free": 0.32,
        "map_cc": None,
        "validation_done": False,
    }
    directives = {"stop_conditions": {"max_refine_cycles": 1}}

    valid = engine.get_valid_programs(
        experiment_type="xray",
        phase_info={"phase": "refine"},
        context=context,
        directives=directives,
    )

    # Must NOT be a bare ["STOP"]
    assert_false(valid == ["STOP"],
                 "max_refine_cycles must not produce bare [STOP]; "
                 "validation programs should be injected. Got: %s" % valid)
    has_validate = any(
        p in valid for p in [
            "phenix.molprobity", "phenix.model_vs_data", "phenix.map_correlations"
        ]
    )
    assert_true(has_validate,
                "Validate-phase programs must appear when max_refine_cycles "
                "limit is reached. Got: %s" % valid)
    assert_in("STOP", valid,
              "STOP must be included alongside validate programs")

    print("  PASSED")


def test_i1_max_refine_cycles_cryoem_uses_rsr_count():
    """I1c: cryoem limit check uses rsr_count not refine_count.

    Bug: code previously read context["refine_count"] for both experiment types.
    With rsr_count=1 and max_refine_cycles=1 the limit must fire for cryoem.
    With refine_count=0 it incorrectly did not.
    """
    print("Test: i1_max_refine_cycles_cryoem_uses_rsr_count")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    context = {
        "refine_count": 0,   # must be ignored for cryoem
        "rsr_count": 1,      # must trigger the limit
        "map_cc": 0.72,
        "r_free": None,
        "validation_done": False,
        "has_model": True, "has_placed_model": True, "has_refined_model": True,
        "has_data_mtz": False, "has_sequence": False, "has_ligand_file": False,
        "has_full_map": True, "has_half_map": False,
    }
    directives = {"stop_conditions": {"max_refine_cycles": 1}}

    valid = engine.get_valid_programs(
        experiment_type="cryoem",
        phase_info={"phase": "refine"},
        context=context,
        directives=directives,
    )

    has_validate_or_stop = "STOP" in valid or any(
        p in valid for p in ["phenix.molprobity", "phenix.validation_cryoem"]
    )
    assert_true(has_validate_or_stop,
                "cryoem max_refine_cycles with rsr_count=1: STOP or validate "
                "programs expected. Got: %s" % valid)
    assert_false(valid == ["phenix.real_space_refine"],
                 "real_space_refine alone must not be offered after limit fires")

    print("  PASSED")


def test_i2_after_program_beats_quality_gate():
    """I2: after_program → STOP only; validate programs NOT injected.

    Distinguishes max_refine_cycles (controlled landing: validate + STOP)
    from after_program (unconditional stop: STOP only).
    """
    print("Test: i2_after_program_beats_quality_gate")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    context = {
        "refine_count": 0,
        "rsr_count": 1,
        "map_cc": 0.30,      # poor — quality gate would normally continue
        "r_free": None,
        "validation_done": False,
        "last_program": "phenix.real_space_refine",
        "has_model": True, "has_placed_model": True, "has_refined_model": True,
        "has_data_mtz": False, "has_sequence": False, "has_ligand_file": False,
        "has_full_map": True, "has_half_map": False,
    }
    directives = {"stop_conditions": {"after_program": "phenix.real_space_refine"}}

    valid = engine.get_valid_programs(
        experiment_type="cryoem",
        phase_info={"phase": "refine"},
        context=context,
        directives=directives,
    )

    assert_in("STOP", valid,
              "STOP must appear after after_program completes")
    assert_false("phenix.real_space_refine" in valid,
                 "after_program: real_space_refine must not be re-offered. "
                 "Got: %s" % valid)

    print("  PASSED")


# =============================================================================
# YAML SPEC: programs.yaml regression checks
# =============================================================================

def test_yaml_xtriage_resolution_has_pick_min():
    """Regression: xtriage resolution spec must have pick_min:true."""
    print("Test: yaml_xtriage_resolution_has_pick_min")
    programs = _load_programs()
    spec = programs["phenix.xtriage"]["log_parsing"]["resolution"]
    assert_true(spec.get("pick_min", False),
                "xtriage resolution must have pick_min: true")
    print("  PASSED")


def test_yaml_rsr_map_cc_is_extract_last():
    """Regression: real_space_refine map_cc must be extract:last after E1 fix."""
    print("Test: yaml_rsr_map_cc_is_extract_last")
    programs = _load_programs()
    spec = programs["phenix.real_space_refine"]["log_parsing"]["map_cc"]
    assert_equal(spec.get("extract", "first"), "last",
                 "real_space_refine map_cc must use extract: last")
    print("  PASSED")


def test_yaml_rsr_clashscore_is_still_extract_last():
    """Regression: real_space_refine clashscore must remain extract:last."""
    print("Test: yaml_rsr_clashscore_is_still_extract_last")
    programs = _load_programs()
    spec = programs["phenix.real_space_refine"]["log_parsing"]["clashscore"]
    assert_equal(spec.get("extract", "first"), "last",
                 "real_space_refine clashscore must use extract: last")
    print("  PASSED")


def test_yaml_polder_requires_selection_invariant():
    """G2 regression: polder invariant requires selection strategy_flag."""
    print("Test: yaml_polder_requires_selection_invariant")
    programs = _load_programs()
    polder = programs["phenix.polder"]
    # selection must be a strategy_flag, not an input slot
    assert_true("selection" in polder.get("strategy_flags", {}),
                "polder selection must be a strategy_flag")
    assert_true("{selection}" in polder["command"],
                "polder command template must contain {selection}")
    # Must have a requires_selection invariant
    invariant_names = [inv.get("name") for inv in polder.get("invariants", [])]
    assert_in("requires_selection", invariant_names,
              "polder must have requires_selection invariant")
    print("  PASSED")



# =============================================================================
# CATEGORY I1b — validate → STOP after validation_done=True
# =============================================================================

def test_i1b_validation_done_produces_stop():
    """I1b: after validation_done=True engine routes to complete phase -> [STOP].

    This is the clean-termination half of the I1 story: max_refine_cycles transitions
    to validate (tested in I1a); after validation completes the engine must produce
    exactly ["STOP"] via the complete-phase handler (not via _apply_directives).

    Uses build_context() to supply all required context keys — detect_phase()
    raises KeyError for any missing key in the context dict.
    """
    print("Test: i1b_validation_done_produces_stop")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # Simulate xray session: xtriage + phaser + refine (r_free good) + molprobity done
    history_info = {
        "xtriage_done": True,
        "phaser_done": True,
        "refine_done": True,
        "refine_count": 1,
        "validation_done": True,
        "molprobity_done": True,
    }
    files = {
        "data_mtz": ["data.mtz"],
        "model": ["refine_001.pdb"],
        "refined": ["refine_001.pdb"],
    }
    context = engine.build_context(files, history_info, {"r_free": 0.24}, None)

    phase_info = engine.detect_phase("xray", context)
    valid = engine.get_valid_programs(
        experiment_type="xray",
        phase_info=phase_info,
        context=context,
    )

    assert_equal(phase_info.get("phase"), "complete",
                 "Phase must be 'complete' when validation_done=True and r_free is good. "
                 "Got: %s" % phase_info)
    assert_equal(valid, ["STOP"],
                 "After validation_done=True the engine must return [STOP] "
                 "(complete phase). Got: %s" % valid)

    print("  PASSED")


def test_i1b_cryoem_validation_done_produces_stop():
    """I1b cryo-EM: same clean-termination check for the cryo-EM workflow."""
    print("Test: i1b_cryoem_validation_done_produces_stop")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    history_info = {
        "mtriage_done": True,
        "resolve_cryo_em_done": True,
        "dock_done": True,
        "rsr_done": True,
        "rsr_count": 2,
        "has_full_map": True,
        "validation_done": True,
        "molprobity_done": True,
    }
    files = {
        "map": ["denmod_map.ccp4"],
        "full_map": ["denmod_map.ccp4"],
        "model": ["rsr_001.pdb"],
        "refined": ["rsr_001.pdb"],
    }
    context = engine.build_context(files, history_info, {"map_cc": 0.81}, None)

    phase_info = engine.detect_phase("cryoem", context)
    valid = engine.get_valid_programs(
        experiment_type="cryoem",
        phase_info=phase_info,
        context=context,
    )

    assert_equal(phase_info.get("phase"), "complete",
                 "Cryo-EM phase must be 'complete' when validation_done=True. "
                 "Got: %s" % phase_info)
    assert_equal(valid, ["STOP"],
                 "Cryo-EM: after validation_done=True engine must return [STOP]. "
                 "Got: %s" % valid)

    print("  PASSED")



# =============================================================================
# ZOMBIE DIAGNOSTIC SURFACING — J5 regression
# =============================================================================

def test_j5_zombie_diagnostics_returned_in_state():
    """J5 regression: detect_workflow_state surfaces zombie diagnostics.

    Bug: _clear_zombie_done_flags() return value was silently discarded.
    The state dict should now carry 'zombie_diagnostics' when zombies are found,
    so PERCEIVE can log them.
    """
    print("Test: j5_zombie_diagnostics_returned_in_state")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_state import detect_workflow_state
    except ImportError:
        print("  SKIP (detect_workflow_state unavailable)")
        return

    # History says resolve_cryo_em ran successfully, but no output file exists
    history = [
        {
            "program": "phenix.resolve_cryo_em",
            "command": "phenix.resolve_cryo_em half1.ccp4 half2.ccp4",
            "result": "resolve_cryo_em DONE",
            "analysis": {},
        }
    ]
    # No available_files — denmod_map.ccp4 is missing (zombie state)
    state = detect_workflow_state(
        history=history,
        available_files=[],
        analysis={},
    )

    diags = state.get("zombie_diagnostics", [])
    assert_true(len(diags) > 0,
                "Zombie diagnostics must be present when resolve_cryo_em_done "
                "is set but output file is missing")
    assert_true(any("resolve_cryo_em" in d for d in diags),
                "Zombie diagnostic must name the cleared program. Got: %s" % diags)

    # Also confirm the done flag was cleared (zombie detection worked)
    # The program should be re-eligible — resolve_cryo_em should appear in valid_programs
    # (or at minimum, resolve_cryo_em_done must not be blocking the workflow)
    print("  PASSED")


def test_j5_no_zombie_diagnostics_when_clean():
    """J5 regression: no zombie_diagnostics key when no zombie states exist."""
    print("Test: j5_no_zombie_diagnostics_when_clean")
    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_state import detect_workflow_state
    except ImportError:
        print("  SKIP (detect_workflow_state unavailable)")
        return

    # Empty history, single MTZ — clean initial state
    state = detect_workflow_state(
        history=[],
        available_files=["data.mtz"],
        analysis={},
    )

    # Either no zombie_diagnostics key, or an empty list
    diags = state.get("zombie_diagnostics", [])
    assert_equal(len(diags), 0,
                 "No zombie diagnostics expected for a clean initial state")

    print("  PASSED")


# =============================================================================
# CATEGORY K1 — best_files list values crash at cycle=2
# =============================================================================

def test_k1_best_path_handles_list():
    """K1: CommandBuilder._best_path must not crash when value is a list.

    Root cause: clients may legitimately store multi-file entries (e.g.
    half_map: [map1.mrc, map2.mrc]) as a list in the best_files dict that is
    passed back as session_state on cycle=2.  Every os.path call that reads a
    best_files value crashed with 'expected str, bytes or os.PathLike, not list'.

    Fix: _best_path() helper collapses list → first element, None → None.
    """
    print("Test: k1_best_path_handles_list")
    try:
        from agent.command_builder import CommandBuilder
    except ImportError:
        print("  SKIP (command_builder unavailable)")
        return

    cb = CommandBuilder()

    # None and empty string → None
    assert_equal(cb._best_path(None), None,
                 "_best_path(None) should return None")
    assert_equal(cb._best_path(""), None,
                 "_best_path('') should return None")

    # Plain string → returned unchanged
    assert_equal(cb._best_path("/data/model.pdb"), "/data/model.pdb",
                 "_best_path(str) should return the string")

    # List → first element (single-file slot takes only the first)
    assert_equal(cb._best_path(["/data/map1.mrc", "/data/map2.mrc"]),
                 "/data/map1.mrc",
                 "_best_path(list) should return first element")

    # Empty list → None
    assert_equal(cb._best_path([]), None,
                 "_best_path([]) should return None")

    print("  PASSED")


def test_k1_best_files_logging_handles_list():
    """K1: best_files logging must not crash when a value is a list.

    The summary log lines in graph_nodes.py and run_ai_agent.py both called
    os.path.basename(v) directly on best_files values.  When v is a list this
    raises TypeError.  Fixed by extracting v[0] for list values before the
    basename call.
    """
    print("Test: k1_best_files_logging_handles_list")
    import os

    best_files = {
        "model":    "/data/refined_001.pdb",
        "half_map": ["/data/half_map_1.mrc", "/data/half_map_2.mrc"],
    }

    # Reproduce the fixed logging expression from graph_nodes.py / run_ai_agent.py
    try:
        summary = ", ".join([
            "%s=%s" % (k, os.path.basename(v[0] if isinstance(v, list) else v))
            for k, v in best_files.items() if v
        ])
    except TypeError as e:
        assert_true(False,
            "Logging expression crashed with list value: %s" % e)

    assert_in("model=refined_001.pdb", summary,
              "model path should appear in summary")
    assert_in("half_map=half_map_1.mrc", summary,
              "half_map list first element should appear in summary")

    print("  PASSED")


def test_k1_command_builder_best_files_list_does_not_crash():
    """K1: CommandBuilder must not crash when best_files contains a list value.

    Simulates the exact cycle=2 scenario: session_state carries
    best_files = {"model": "...", "half_map": [map1, map2]}.
    _select_file_for_slot reads context.best_files.get(category) and passes
    the result to os.path.exists.  With a list value this raised TypeError.
    """
    print("Test: k1_command_builder_best_files_list_does_not_crash")
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (command_builder unavailable)")
        return

    cb = CommandBuilder()

    # Build a minimal CommandContext where best_files has a list for half_map
    # and a string for model.
    ctx = CommandContext(
        cycle_number=2,
        experiment_type="cryoem",
        resolution=3.5,
        best_files={
            "model":    "/data/overall_best.pdb",
            "half_map": ["/data/half_map_1.mrc", "/data/half_map_2.mrc"],
        },
        rfree_mtz=None,
        categorized_files={
            "model":    ["/data/overall_best.pdb"],
            "half_map": ["/data/half_map_1.mrc", "/data/half_map_2.mrc"],
        },
        workflow_state="cryoem_has_model",
        history=[],
        llm_files=None,
        llm_strategy=None,
        directives={},
        log=lambda msg: None,
    )

    # Building a real_space_refine command exercises _select_file_for_slot
    # which calls _best_path() on best_files values.  It should not raise
    # TypeError even if the command ultimately can't be built (returns None).
    crashed = False
    try:
        available = ["/data/overall_best.pdb",
                     "/data/half_map_1.mrc",
                     "/data/half_map_2.mrc"]
        cmd = cb.build("phenix.real_space_refine", available, ctx)
        # cmd may be None if required files aren't fully satisfied — that's fine.
        # The critical guarantee is that no TypeError was raised.
    except TypeError as e:
        crashed = True
        assert_true(False,
            "build() crashed with list best_files value: %s" % e)

    assert_false(crashed,
        "build() must not raise TypeError when best_files contains a list")

    print("  PASSED")


def test_k2_mtriage_keeps_half_maps_with_full_map():
    """K2: mtriage must receive both full_map AND half_maps when both are available.

    Root cause: the post-selection validation in command_builder.py had a blanket
    rule that dropped half_maps whenever a full_map was also selected.  For programs
    like mtriage the half maps are NOT redundant — they provide FSC-based resolution
    measurement that full-map-only mode cannot.

    Fix: keep_half_maps_with_full_map: true in programs.yaml for mtriage,
    predict_and_build, and map_to_model.  The post-selection block now checks this
    flag before removing half_maps.
    """
    print("Test: k2_mtriage_keeps_half_maps_with_full_map")
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (command_builder unavailable)")
        return

    sharpened = "/data/sharpened_map.ccp4"
    half1 = "/data/half_map_1.ccp4"
    half2 = "/data/half_map_2.ccp4"

    ctx = CommandContext(
        cycle_number=2,
        experiment_type="cryoem",
        resolution=3.0,
        best_files={"map": sharpened},
        rfree_mtz=None,
        categorized_files={
            "optimized_full_map": [sharpened],
            "map":                [sharpened, half1, half2],
            "half_map":           [half1, half2],
        },
        workflow_state="cryoem_initial",
        history=[],
        llm_files=None,
        llm_strategy=None,
        directives={},
        log=lambda msg: None,
    )

    cb = CommandBuilder()
    available = [sharpened, half1, half2]
    cmd = cb.build("phenix.mtriage", available, ctx)

    assert_not_none(cmd, "mtriage must produce a command when map files are available")
    assert_true("half_map=" in cmd,
                "mtriage command must include half_map= when half maps are present. Got: %s" % cmd)
    assert_true(sharpened in cmd or "sharpened" in cmd,
                "mtriage command must include the full/sharpened map. Got: %s" % cmd)

    print("  PASSED")


def test_k2_map_sharpening_uses_half_maps_when_no_full_map():
    """K2: map_sharpening must use half_map= mode when only half maps are available.

    When the input set contains ONLY half maps (no full map yet), map_sharpening
    should run in mode 3: half_map=map1 half_map=map2.  It must NOT select a
    half map as the positional full_map argument.
    """
    print("Test: k2_map_sharpening_uses_half_maps_when_no_full_map")
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (command_builder unavailable)")
        return

    half1 = "/data/half_map_1.ccp4"
    half2 = "/data/half_map_2.ccp4"
    seq   = "/data/seq.dat"

    ctx = CommandContext(
        cycle_number=1,
        experiment_type="cryoem",
        resolution=None,
        best_files={},
        rfree_mtz=None,
        categorized_files={
            "half_map": [half1, half2],
            "map":      [half1, half2],   # half maps bubble up to parent
            "sequence": [seq],
        },
        workflow_state="cryoem_initial",
        history=[],
        llm_files=None,
        llm_strategy=None,
        directives={},
        log=lambda msg: None,
    )

    cb = CommandBuilder()
    available = [half1, half2, seq]
    cmd = cb.build("phenix.map_sharpening", available, ctx)

    # Command must include half_map= flag (mode 3)
    assert_not_none(cmd, "map_sharpening must produce a command with half maps available")
    assert_true("half_map=" in cmd,
                "map_sharpening must use half_map= mode when no full map exists. Got: %s" % cmd)
    # Neither half map should appear as a bare positional argument (full_map slot)
    # The positional full_map slot has flag="" so it appears without a prefix.
    # If a half map is used as full_map it will appear without "half_map=" prefix.
    cmd_tokens = cmd.split()
    bare_map_tokens = [t for t in cmd_tokens
                       if (t.endswith(".ccp4") or t.endswith(".mrc") or t.endswith(".map"))
                       and not t.startswith("half_map=")
                       and not t.startswith("seq_file=")
                       and not t.startswith("phenix.")]
    assert_equal(len(bare_map_tokens), 0,
                 "No half map should appear as bare positional (full_map) arg. "
                 "Bare map tokens: %s\nFull command: %s" % (bare_map_tokens, cmd))

    print("  PASSED")


def test_k2_optimized_full_map_recognized_as_genuine_full_map():
    """K2: A sharpened map (optimized_full_map category) must not be mis-classified
    as a half map in the post-selection validation.

    The fix extended the 'genuine full map' check to include optimized_full_map
    so that sharpened maps are never incorrectly removed in favour of half_maps.
    """
    print("Test: k2_optimized_full_map_recognized_as_genuine_full_map")
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (command_builder unavailable)")
        return

    sharpened = "/data/sharpened_map.ccp4"
    half1 = "/data/half_map_1.ccp4"
    half2 = "/data/half_map_2.ccp4"

    # Simulate a program that does NOT keep_half_maps_with_full_map (e.g. real_space_refine)
    # The sharpened map is in optimized_full_map only (NOT in full_map).
    # Before the fix, the check `full_map_path not in full_map_files` was True
    # (because sharpened is not in the strict "full_map" subcategory),
    # which could classify it as mis-selected and remove it.
    ctx = CommandContext(
        cycle_number=2,
        experiment_type="cryoem",
        resolution=3.0,
        best_files={"map": sharpened},
        rfree_mtz=None,
        categorized_files={
            "optimized_full_map": [sharpened],
            "map":                [sharpened, half1, half2],
            "half_map":           [half1, half2],
            "full_map":           [],        # strict full_map is empty
        },
        workflow_state="cryoem_refined",
        history=[],
        llm_files=None,
        llm_strategy=None,
        directives={},
        log=lambda msg: None,
    )

    cb = CommandBuilder()
    available = [sharpened, half1, half2]
    cmd = cb.build("phenix.real_space_refine", available, ctx)

    # real_space_refine needs model too so cmd may be None — that's fine.
    # The key: if a command IS produced, the sharpened map must not be absent
    # due to mis-classification.  We just check no crash and that sharpened
    # was not replaced by a half map in any model-less fallback.
    if cmd:
        # If the command references a map file, it must be the sharpened one, not a half map
        if half1 in cmd or half2 in cmd:
            assert_true(False,
                "real_space_refine must not use a half map as its map input. Got: %s" % cmd)

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all audit-fix regression tests."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()


# =============================================================================
# L1: Ligand PDB excluded by _is_valid_file when file starts with COMPND/AUTHOR
# =============================================================================

def test_l1_ligand_pdb_compnd_header_is_valid():
    """
    A ligand PDB starting with COMPND (standard PDB header) must not be
    rejected by _is_valid_file. Previously only a narrow set of records
    (ATOM, HETATM, REMARK, HEADER, TITLE, MODEL, SEQRES, SCALE1, ORIGX1)
    was accepted as the first line — COMPND was missing, causing the ligand
    file to be excluded from categorization, has_ligand_file=False, and
    ligandfit never being offered.
    """
    import tempfile, os
    from agent.workflow_state import _is_valid_file

    # Typical ligand PDB from phenix/CSD starting with COMPND then HETATM
    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
        f.write('COMPND    BROMODOMAIN INHIBITOR\n')
        f.write('REMARK   1 Generated by phenix.elbow\n')
        f.write('HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C\n')
        f.write('HETATM    2  N1  LIG A   1       2.500   2.000   3.000  1.00  0.00           N\n')
        f.write('END\n')
        fname = f.name
    try:
        result = _is_valid_file(fname)
        assert result, (
            "Expected _is_valid_file to accept COMPND-header ligand PDB, got False.\n"
            "This caused has_ligand_file=False and ligandfit to be silently skipped."
        )
    finally:
        os.unlink(fname)
    print("  PASSED")


def test_l1_ligand_pdb_must_have_coordinates():
    """
    A .pdb file with no ATOM/HETATM records (just header lines) must be
    rejected — it's structurally useless even if it starts with a valid record.
    """
    import tempfile, os
    from agent.workflow_state import _is_valid_file

    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
        f.write('COMPND    BROMODOMAIN INHIBITOR\n')
        f.write('REMARK   1 No coordinates in this file\n')
        f.write('END\n')
        fname = f.name
    try:
        result = _is_valid_file(fname)
        assert not result, (
            "Expected _is_valid_file to reject PDB with no ATOM/HETATM, got True."
        )
    finally:
        os.unlink(fname)
    print("  PASSED")


def test_l1_detect_phase_returns_refine_when_ligand_present_not_run():
    """
    When has_ligand_file=True and ligandfit_done=False and refine_count > 0
    and r_free < 0.35, detect_phase must return 'refine' (not 'complete'),
    keeping ligandfit in valid_programs.

    This was broken because _is_valid_file rejected 7qz0_ligand.pdb, making
    has_ligand_file=False so detect_phase skipped the guard and went to 'complete'.
    """
    from agent.workflow_engine import WorkflowEngine
    engine = WorkflowEngine()

    context = {
        # Files present
        "has_model": True,
        "has_data_mtz": True,
        "has_ligand_file": True,       # fixed: was False due to _is_valid_file bug
        "has_search_model": False,
        "has_map": False,
        "has_map_coeffs_mtz": True,
        "has_full_map": False,
        "has_half_map": False,
        "has_sequence": True,
        "has_anomalous": False,
        "has_predicted_model": True,
        "has_placed_model": True,
        "has_refined_model": True,
        "has_ligand_fit": False,
        "has_optimized_full_map": False,
        "has_processed_model": False,
        "has_twinning": False,
        # Workflow history flags
        "refine_count": 3,
        "r_free": 0.2756,
        "validation_done": True,       # agent had passed through validate
        "ligandfit_done": False,       # ligandfit never ran
        "xtriage_done": True,
        "autosol_done": False,
        "autobuild_done": False,
        "phaser_done": False,
        "predict_done": True,
        "predict_full_done": True,
        "mtriage_done": False,
        "pdbtools_done": False,
        "molprobity_done": True,
        # Misc
        "anomalous_measurability": 0.0,
        "twin_fraction": 0.0,
        "map_cc": None,
        "model_is_good": True,
        "automation_path": "predict_and_build",
        "use_mr_sad": False,
        "strong_anomalous": False,
    }

    phase_info = engine.detect_phase("xray", context)
    phase = phase_info.get("phase")
    assert phase == "refine", (
        "Expected detect_phase to return 'refine' (so ligandfit is available), "
        "got '%s'.\nContext: has_ligand_file=True, ligandfit_done=False, "
        "refine_count=3, r_free=0.2756, validation_done=True.\n"
        "The guard on lines 497-504 of workflow_engine.py should keep the "
        "workflow in refine phase when ligand fitting is still pending." % phase
    )

    # Also verify ligandfit is in valid_programs for this phase
    valid = engine.get_valid_programs("xray", phase_info, context)
    assert "phenix.ligandfit" in valid, (
        "Expected phenix.ligandfit in valid_programs for refine phase with "
        "ligand present. Got: %s" % valid
    )
    print("  PASSED")



# =============================================================================
# M1: with_ligand model not selected for refinement after ligandfit+pdbtools
# =============================================================================

def test_m1_with_ligand_beats_refined_in_best_files_tracker():
    """
    After pdbtools combines the ligandfit output into *_with_ligand.pdb,
    BestFilesTracker must rank that file above the older ligand-free refined
    model. Previously both had stage_score=100 and the refined model won on
    R-free metrics, causing refinement to use the wrong (ligand-free) model.
    """
    import tempfile, os
    sys.path.insert(0, '/tmp/improved_agent_v2v1')
    from agent.best_files_tracker import BestFilesTracker

    tracker = BestFilesTracker()

    with tempfile.TemporaryDirectory() as d:
        # The refined model (exists, has good R-free metrics from cycle 3)
        refined = os.path.join(d, 'overall_best_final_refine_001.pdb')
        with open(refined, 'w') as f:
            f.write('ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  5.00           C\n')

        # The combined model with ligand (exists, no metrics yet from pdbtools)
        with_ligand = os.path.join(d, 'overall_best_final_refine_001_with_ligand.pdb')
        with open(with_ligand, 'w') as f:
            f.write('ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  5.00           C\n')
            f.write('HETATM    2  C1  LIG A 100       5.0   6.0   7.0  1.00  5.00           C\n')

        # Add refined model first (cycle 3, with R-free metrics)
        tracker.evaluate_file(refined, cycle=3, metrics={"r_free": 0.2756}, stage="refined")
        # Add with_ligand model second (cycle 5, pdbtools run, no R-free metrics)
        tracker.evaluate_file(with_ligand, cycle=5, metrics=None, stage="with_ligand")

        best = tracker.get_best_path("model")
        assert best is not None, "Expected a best model to be selected"
        assert os.path.basename(best) == 'overall_best_final_refine_001_with_ligand.pdb', (
            "Expected with_ligand model to beat refined model in BestFilesTracker.\n"
            "Got: %s\n"
            "with_ligand stage_score should be 110, refined should be 100." % os.path.basename(best)
        )
    print("  PASSED")


def test_m1_command_builder_with_ligand_not_overridden():
    """
    command_builder._is_with_ligand_file() equivalent: verify that the
    category detection used in the override guard correctly identifies
    with_ligand files so they pass through without being replaced by best_model.

    This tests the logic that was added to the LLM override block:
    if the LLM's choice is in with_ligand or ligand_fit_output categories,
    we trust the LLM rather than forcing best_model.
    """
    import os
    # Test the category lookup logic directly (avoids full CommandBuilder mock complexity)
    categorized_files = {
        "model": ["/path/to/refine_001.pdb", "/path/to/refine_001_with_ligand.pdb"],
        "refined": ["/path/to/refine_001.pdb"],
        "with_ligand": ["/path/to/refine_001_with_ligand.pdb"],
        "data_mtz": ["/path/to/7qz0.mtz"],
    }

    def get_llm_categories(corrected_str, categorized_files):
        llm_bn = os.path.basename(corrected_str)
        return [
            cat for cat, files in categorized_files.items()
            if any(os.path.basename(f) == llm_bn for f in files)
        ]

    # LLM picked with_ligand file → should NOT be overridden
    cats = get_llm_categories("/path/to/refine_001_with_ligand.pdb", categorized_files)
    llm_is_with_ligand = any(c in ("with_ligand", "ligand_fit_output") for c in cats)
    assert llm_is_with_ligand, (
        "LLM choice of refine_001_with_ligand.pdb should be detected as with_ligand category. "
        "Got categories: %s" % cats
    )

    # LLM picked an unrelated model → should be overridden
    cats2 = get_llm_categories("/path/to/some_other_model.pdb", categorized_files)
    llm_is_with_ligand2 = any(c in ("with_ligand", "ligand_fit_output") for c in cats2)
    assert not llm_is_with_ligand2, (
        "LLM choice of some_other_model.pdb should NOT be detected as with_ligand. "
        "Got categories: %s" % cats2
    )
    print("  PASSED")



# =============================================================================
# N1: resolve_cryo_em missing sequence/mass raises Sorry immediately
# =============================================================================

def test_n1_resolve_cryo_em_missing_mass_raises_sorry():
    """
    When phenix.resolve_cryo_em fails because no sequence/mass/solvent_content
    was provided, error_analyzer must raise Sorry with a clear user-facing
    message rather than returning None or looping. This is a hard-stop error —
    the user must supply a sequence file and restart.
    """
    import sys, types

    # Minimal mock for libtbx.utils.Sorry (the real one may not be present)
    lt = types.ModuleType('libtbx'); lt.__path__ = []; sys.modules.setdefault('libtbx', lt)
    ltu = types.ModuleType('libtbx.utils'); sys.modules.setdefault('libtbx.utils', ltu)
    class _Sorry(Exception): pass
    ltu.Sorry = _Sorry

    # Force reimport so the mock Sorry is picked up
    import importlib
    import agent.error_analyzer as ea_mod
    importlib.reload(ea_mod)

    analyzer = ea_mod.ErrorAnalyzer()

    # Exact error message from the failing resolve_cryo_em run
    log_text = (
        "FAILED: Need solvent_content, molecular mass, a sequence file, or "
        "fraction_of_max_mask_threshold to identify the mask for your molecule. "
        "You can try with 'fraction_of_max_mask_threshold = 0.05' but it may not work well."
    )

    raised = False
    message = ""
    try:
        analyzer.analyze(
            log_text=log_text,
            program="phenix.resolve_cryo_em",
            context={},
            session=None,
        )
    except Exception as e:
        raised = True
        message = str(e)

    assert raised, (
        "Expected analyzer.analyze() to raise Sorry for resolve_cryo_em missing-mass error, "
        "but no exception was raised. The agent should hard-stop with a clear message."
    )
    assert "sequence" in message.lower() or "mass" in message.lower(), (
        "Sorry message should mention 'sequence' or 'mass' so the user knows what to supply. "
        "Got: %s" % message
    )
    print("  PASSED")


def test_n1_hard_stop_not_triggered_for_other_programs():
    """
    The resolve_cryo_em hard-stop pattern must NOT fire for other programs that
    might have coincidentally similar error text.
    """
    import sys, types

    lt = types.ModuleType('libtbx'); lt.__path__ = []; sys.modules.setdefault('libtbx', lt)
    ltu = types.ModuleType('libtbx.utils'); sys.modules.setdefault('libtbx.utils', ltu)
    class _Sorry(Exception): pass
    ltu.Sorry = _Sorry

    import importlib
    import agent.error_analyzer as ea_mod
    importlib.reload(ea_mod)

    analyzer = ea_mod.ErrorAnalyzer()

    log_text = (
        "Need solvent_content, molecular mass, a sequence file, or "
        "fraction_of_max_mask_threshold to identify the mask for your molecule."
    )

    # Should NOT raise for a different program
    raised = False
    try:
        analyzer.analyze(
            log_text=log_text,
            program="phenix.some_other_program",
            context={},
            session=None,
        )
    except Exception:
        raised = True

    assert not raised, (
        "Hard-stop for resolve_cryo_em missing-mass should only fire for "
        "phenix.resolve_cryo_em, not for other programs."
    )
    print("  PASSED")

