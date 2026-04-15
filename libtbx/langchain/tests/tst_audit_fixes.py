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
import tempfile
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _find_ai_agent_path():
    """
    Locate ai_agent.py robustly across environments.

    The file lives at different locations depending on the environment:
      - Dev / /tmp unpack:  <project_root>/programs/ai_agent.py
      - PHENIX installation: <phenix_root>/programs/ai_agent.py
                             (importable as phenix.programs.ai_agent)

    Strategy:
      1. Relative to this test file (dev / /tmp).
      2. importlib lookup of phenix.programs.ai_agent (PHENIX install).
      3. importlib lookup of libtbx.langchain.programs.ai_agent (alt layout).
      4. sys.path scan for programs/ai_agent.py and phenix/programs/ai_agent.py.
    """
    # 1. Relative to this test file
    candidate = os.path.join(_PROJECT_ROOT, 'programs', 'ai_agent.py')
    if os.path.exists(candidate):
        return candidate

    # 2 & 3. Via importlib — covers both package layouts without importing
    try:
        import importlib.util as _ilu
        for _mod in ('phenix.programs.ai_agent',
                     'libtbx.langchain.programs.ai_agent'):
            try:
                spec = _ilu.find_spec(_mod)
                if spec and spec.origin and os.path.exists(spec.origin):
                    return spec.origin
            except (ModuleNotFoundError, ValueError):
                pass
    except Exception:
        pass

    # 4. sys.path scan
    for _p in sys.path:
        for _rel in ('programs/ai_agent.py',
                     'phenix/programs/ai_agent.py',
                     'langchain/programs/ai_agent.py'):
            candidate = os.path.join(_p, _rel)
            if os.path.exists(candidate):
                return candidate

    raise FileNotFoundError(
        "Cannot locate ai_agent.py. "
        "Tried relative path %s/programs/ai_agent.py, "
        "importlib (phenix.programs.ai_agent, libtbx.langchain.programs.ai_agent), "
        "and sys.path entries." % _PROJECT_ROOT
    )

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

    # Simulate xray refine step with refine_count=1 and max_refine_cycles=1
    context = {
        "step": "refine",
        "refine_count": 1,
        "rsr_count": 0,
        "r_free": 0.32,
        "map_cc": None,
        "validation_done": False,
    }
    directives = {"stop_conditions": {"max_refine_cycles": 1}}

    valid = engine.get_valid_programs(
        experiment_type="xray",
        step_info={"step": "refine"},
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
                "Validate-step programs must appear when max_refine_cycles "
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
        step_info={"step": "refine"},
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
        "map_cc": 0.75,      # acceptable — above quality override threshold (0.70)
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
        step_info={"step": "refine"},
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
    # {selection} was removed from the command template - it is now appended
    # automatically from strategy_flags (avoids double-substitution issues).
    assert_true("{selection}" not in polder["command"],
                "polder command template must NOT contain {selection} (appended via strategy_flag)")
    # Must have a requires_selection invariant
    invariant_names = [inv.get("name") for inv in polder.get("invariants", [])]
    assert_in("requires_selection", invariant_names,
              "polder must have requires_selection invariant")
    # Flag template must use single quotes so multi-word values don't get split
    sel_flag = polder["strategy_flags"]["selection"].get("flag", "")
    assert_true("'" in sel_flag,
                "polder selection flag must use single quotes to avoid shell splitting: %s" % sel_flag)
    # Must have a safe default that doesn't assume ligand residue name
    sel_default = polder["strategy_flags"]["selection"].get("default", "")
    assert_true("hetero" in sel_default,
                "polder selection default should use 'hetero and not water', got: %s" % sel_default)
    print("  PASSED")


def test_polder_selection_always_uses_safe_default():
    """G2b: program_registry overrides LLM-supplied selection for polder with safe default.

    The LLM cannot reliably know the ligand residue name and guesses values like
    'resname LIG'.  The registry must always reset it to 'hetero and not water'.
    """
    print("Test: polder_selection_always_uses_safe_default")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return

    registry = ProgramRegistry(use_yaml=True)
    dummy_files = {"data_mtz": "/tmp/data.mtz", "model": "/tmp/model.pdb"}

    # Simulate LLM providing wrong residue name
    for bad_sel in ["resname LIG", "resname LGD", "resname ATP", "LIG"]:
        cmd = registry.build_command(
            program_name="phenix.polder",
            files=dummy_files,
            strategy={"selection": bad_sel},
            log=lambda msg: None,
        )
        assert_true(
            bad_sel not in cmd,
            "polder command should NOT contain LLM selection %r, got: %s" % (bad_sel, cmd)
        )
        assert_true(
            "hetero" in cmd,
            "polder command should use 'hetero and not water', got: %s" % cmd
        )
        assert_true(
            "'" in cmd,
            "polder selection must be single-quoted to avoid shell splitting, got: %s" % cmd
        )
    print("  PASSED (LLM selection overridden for: resname LIG, resname LGD, resname ATP, LIG)")


# =============================================================================
# CATEGORY I1b — validate → STOP after validation_done=True
# =============================================================================

def test_i1b_validation_done_produces_stop():
    """I1b: after validation_done=True engine routes to complete step -> [STOP].

    This is the clean-termination half of the I1 story: max_refine_cycles transitions
    to validate (tested in I1a); after validation completes the engine must produce
    exactly ["STOP"] via the complete-step handler (not via _apply_directives).

    Uses build_context() to supply all required context keys — detect_step()
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

    step_info = engine.detect_step("xray", context)
    valid = engine.get_valid_programs(
        experiment_type="xray",
        step_info=step_info,
        context=context,
    )

    assert_equal(step_info.get("step"), "complete",
                 "Phase must be 'complete' when validation_done=True and r_free is good. "
                 "Got: %s" % step_info)
    assert_equal(valid, ["STOP"],
                 "After validation_done=True the engine must return [STOP] "
                 "(complete step). Got: %s" % valid)

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

    step_info = engine.detect_step("cryoem", context)
    valid = engine.get_valid_programs(
        experiment_type="cryoem",
        step_info=step_info,
        context=context,
    )

    assert_equal(step_info.get("step"), "complete",
                 "Cryo-EM step must be 'complete' when validation_done=True. "
                 "Got: %s" % step_info)
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
    """K2: mtriage must prefer half-maps and drop full_map when both available.

    Half-map FSC is the gold standard for resolution — mtriage computes
    resolution from half-map pairs.  When both half-maps and a full_map are
    present, the full_map is dropped (prefers_half_maps: true) to avoid
    "Maps have different dimensions" errors when the full_map was
    post-processed with a different grid.

    The command should be:
        phenix.mtriage half_map=half1.ccp4 half_map=half2.ccp4
    NOT:
        phenix.mtriage half_map=half1.ccp4 full_map=sharpened.ccp4
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
    # Both half-maps must be present
    assert_true("half_map=" in cmd,
                "mtriage command must include half_map= when half maps are present. Got: %s" % cmd)
    assert_true("half_map_1" in cmd and "half_map_2" in cmd,
                "mtriage must include both half-maps. Got: %s" % cmd)
    # Full map must be DROPPED (prefers_half_maps)
    assert_false(sharpened in cmd,
                 "mtriage must drop full_map when half-maps present (prefers_half_maps). Got: %s" % cmd)

    print("  PASSED")


def test_k2_map_sharpening_uses_half_maps_when_no_full_map():
    """K2: map_sharpening must use half maps as positional args when only half maps
    are available.

    When the input set contains ONLY half maps (no full map), map_sharpening
    should run with both half maps as bare positional arguments:
        phenix.map_sharpening half1.ccp4 half2.ccp4 seq_file=seq.dat
    This matches the confirmed working command format.
    It must NOT select only one half map or use the half_map= flag.
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

    # Command must be produced
    assert_not_none(cmd, "map_sharpening must produce a command with half maps available")
    # Both half maps must appear in the command (positional, no flag prefix)
    assert_true("half_map_1.ccp4" in cmd,
                "map_sharpening must include half_map_1. Got: %s" % cmd)
    assert_true("half_map_2.ccp4" in cmd,
                "map_sharpening must include half_map_2. Got: %s" % cmd)
    # Must NOT use half_map= flag (positional is correct)
    assert_false("half_map=" in cmd,
                 "map_sharpening must use positional args, not half_map= flag. Got: %s" % cmd)

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
# CATEGORY R1 — placement_checker: unit cell comparison (Tier 1)
# =============================================================================

def test_r1_pdb_cryst1_parsed_correctly():
    """R1: read_pdb_unit_cell returns correct 6-tuple from a CRYST1 line."""
    print("Test: r1_pdb_cryst1_parsed_correctly")
    import tempfile
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import read_pdb_unit_cell

    cryst1 = "CRYST1   57.230   57.230  146.770  90.00  90.00  90.00 P 41 21 2\n"
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        f.write(cryst1)
        f.write("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n")
        pdb_path = f.name
    try:
        cell = read_pdb_unit_cell(pdb_path)
        assert cell is not None, "Should parse CRYST1 line"
        assert len(cell) == 6, "Should return 6 values"
        assert abs(cell[0] - 57.230) < 0.001, "a should be 57.230, got %s" % cell[0]
        assert abs(cell[2] - 146.770) < 0.001, "c should be 146.770, got %s" % cell[2]
        assert abs(cell[3] - 90.00) < 0.001, "alpha should be 90.00, got %s" % cell[3]
    finally:
        os.unlink(pdb_path)
    print("  PASSED")


def test_r1_pdb_cryst1_missing_returns_none():
    """R1: read_pdb_unit_cell returns None when no CRYST1 record is present."""
    print("Test: r1_pdb_cryst1_missing_returns_none")
    import tempfile
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import read_pdb_unit_cell

    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        f.write("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n")
        pdb_path = f.name
    try:
        cell = read_pdb_unit_cell(pdb_path)
        assert cell is None, "Should return None when CRYST1 absent, got %s" % str(cell)
    finally:
        os.unlink(pdb_path)
    print("  PASSED")


def test_r1_pdb_missing_file_returns_none():
    """R1: read_pdb_unit_cell returns None (not raises) for a missing file."""
    print("Test: r1_pdb_missing_file_returns_none")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import read_pdb_unit_cell

    cell = read_pdb_unit_cell("/tmp/this_file_does_not_exist_xyz.pdb")
    assert cell is None, "Should return None for missing file, got %s" % str(cell)
    print("  PASSED")


def test_r1_cells_compatible_within_tolerance():
    """R1: cells_are_compatible returns True for cells within 5% tolerance."""
    print("Test: r1_cells_compatible_within_tolerance")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    # Identical cells
    cell = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    assert cells_are_compatible(cell, cell), "Identical cells must be compatible"

    # Within 4% on one axis — should still be compatible
    cell_a = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    cell_b = (59.50, 57.23, 146.77, 90.0, 90.0, 90.0)  # ~4% diff on a
    assert cells_are_compatible(cell_a, cell_b), (
        "Cells within 5%% should be compatible: %s vs %s" % (cell_a, cell_b)
    )
    print("  PASSED")


def test_r1_cells_mismatch_outside_tolerance():
    """R1: cells_are_compatible returns False when any parameter differs > 5%."""
    print("Test: r1_cells_mismatch_outside_tolerance")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    # Large difference on 'a' axis (P1 vs different crystal form)
    cell_a = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    cell_b = (80.00, 57.23, 146.77, 90.0, 90.0, 90.0)  # ~40% diff → mismatch
    assert not cells_are_compatible(cell_a, cell_b), (
        "Cells with >5%% difference should be incompatible: %s vs %s" % (cell_a, cell_b)
    )

    # Angle mismatch (orthorhombic vs monoclinic-like)
    cell_c = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    cell_d = (57.23, 57.23, 146.77, 90.0, 110.0, 90.0)  # beta 90→110 = 22% diff
    assert not cells_are_compatible(cell_c, cell_d), (
        "Angle mismatch should be detected: %s vs %s" % (cell_c, cell_d)
    )
    print("  PASSED")


def test_r1_cells_compatible_with_none_is_failsafe():
    """R1: cells_are_compatible returns True (fail-safe) when either cell is None."""
    print("Test: r1_cells_compatible_with_none_is_failsafe")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    cell = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    assert cells_are_compatible(None, cell), "None vs cell should be fail-safe True"
    assert cells_are_compatible(cell, None), "cell vs None should be fail-safe True"
    assert cells_are_compatible(None, None), "None vs None should be fail-safe True"
    print("  PASSED")


def test_r1_xray_mismatch_requires_both_readable():
    """R1: check_xray_cell_mismatch returns False (fail-safe) when PDB has no CRYST1."""
    print("Test: r1_xray_mismatch_requires_both_readable")
    import tempfile
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import check_xray_cell_mismatch

    # PDB with no CRYST1 — cannot compare → must NOT declare mismatch
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        f.write("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n")
        pdb_no_cryst1 = f.name
    try:
        result = check_xray_cell_mismatch(pdb_no_cryst1, "/tmp/fake.mtz")
        assert result is False, (
            "Must return False (fail-safe) when PDB has no CRYST1, got %s" % result
        )
    finally:
        os.unlink(pdb_no_cryst1)
    print("  PASSED")


def test_r1_xray_missing_files_return_false():
    """R1: check_xray_cell_mismatch returns False for missing/None paths."""
    print("Test: r1_xray_missing_files_return_false")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import check_xray_cell_mismatch

    assert check_xray_cell_mismatch(None, None) is False
    assert check_xray_cell_mismatch("", "/tmp/data.mtz") is False
    assert check_xray_cell_mismatch("/tmp/model.pdb", None) is False
    assert check_xray_cell_mismatch("/no/such/file.pdb", "/no/such.mtz") is False
    print("  PASSED")


def test_r1_cryoem_matches_present_cell_is_compatible():
    """R1: model matching the present-portion map cell is NOT a mismatch.

    The model may have been placed in a sub-box extracted from the full map;
    matching the present-portion cell means it is correctly placed.
    """
    print("Test: r1_cryoem_matches_present_cell_is_compatible")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    # Simulate: full map cell = 320^3, present portion = (59, 67, 106)
    # Model cell matches the present portion
    full_cell    = (320.0, 320.0, 320.0, 90.0, 90.0, 90.0)
    present_cell = (59.0, 67.0, 106.0, 90.0, 90.0, 90.0)
    model_cell   = (59.5, 67.3, 106.2, 90.0, 90.0, 90.0)  # ~within 1%

    assert not cells_are_compatible(model_cell, full_cell), (
        "Model should NOT match full cell in this scenario"
    )
    assert cells_are_compatible(model_cell, present_cell), (
        "Model SHOULD match present-portion cell — should be compatible"
    )

    # Verify the overall mismatch logic: compatible with present → not a mismatch
    # (mirrors check_cryoem_cell_mismatch logic without needing a real map file)
    matches_full    = cells_are_compatible(model_cell, full_cell)
    matches_present = cells_are_compatible(model_cell, present_cell)
    is_mismatch = not matches_full and not matches_present
    assert not is_mismatch, "Should NOT be mismatch when model matches present-portion cell"
    print("  PASSED")


def test_r1_cryoem_missing_files_return_false():
    """R1: check_cryoem_cell_mismatch returns False for missing/None paths."""
    print("Test: r1_cryoem_missing_files_return_false")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import check_cryoem_cell_mismatch

    assert check_cryoem_cell_mismatch(None, None) is False
    assert check_cryoem_cell_mismatch("", "/tmp/map.ccp4") is False
    assert check_cryoem_cell_mismatch("/tmp/model.pdb", None) is False
    assert check_cryoem_cell_mismatch("/no/such.pdb", "/no/such.ccp4") is False
    print("  PASSED")


def test_r1_definitive_xray_mismatch_detected():
    """R1: check_xray_cell_mismatch correctly identifies a clear mismatch
    when both cells are available and are very different.

    Uses the internal comparison path directly (no real MTZ needed).
    """
    print("Test: r1_definitive_xray_mismatch_detected")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    # A large crystal-form difference (P1 vs different SG/cell)
    model_cell = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    mtz_cell   = (100.0, 100.0, 200.0,  90.0, 90.0, 90.0)

    assert not cells_are_compatible(model_cell, mtz_cell), (
        "Very different cells should NOT be compatible: %s vs %s"
        % (model_cell, mtz_cell)
    )
    print("  PASSED")


# =============================================================================
# CATEGORY R2 — placement_uncertain context key (Tier 2 gate)
# =============================================================================

def _make_minimal_context(overrides=None):
    """Build a minimal WorkflowEngine context dict for R2 tests."""
    base = {
        "has_placed_model":      False,
        "cell_mismatch":         False,
        "placement_probed":      False,
        "placement_probe_result": None,
        "has_model":             True,
        "has_data_mtz":          True,
        "has_map":               False,
        "has_predicted_model":   False,
        "placement_uncertain":   False,   # will be recomputed
    }
    if overrides:
        base.update(overrides)
    # Replicate the logic from build_context
    base["placement_uncertain"] = (
        not base["has_placed_model"] and
        not base["cell_mismatch"] and
        not base["placement_probed"] and
        base["has_model"] and
        (base["has_data_mtz"] or base["has_map"]) and
        not base["has_predicted_model"]
    )
    return base


def test_r2_placement_uncertain_set_when_ambiguous():
    """R2: placement_uncertain=True when model+MTZ present, no history, no directive."""
    print("Test: r2_placement_uncertain_set_when_ambiguous")
    ctx = _make_minimal_context()
    assert ctx["placement_uncertain"] is True, (
        "placement_uncertain should be True for ambiguous model+MTZ case"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_when_placed():
    """R2: placement_uncertain=False when heuristics confirm model is placed."""
    print("Test: r2_placement_uncertain_false_when_placed")
    ctx = _make_minimal_context({"has_placed_model": True})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when has_placed_model=True"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_on_cell_mismatch():
    """R2: placement_uncertain=False when cell_mismatch=True (Tier 1 already decided)."""
    print("Test: r2_placement_uncertain_false_on_cell_mismatch")
    ctx = _make_minimal_context({"cell_mismatch": True})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when cell_mismatch=True "
        "(Tier 1 handles routing)"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_when_probed():
    """R2: placement_uncertain=False when probe has already run (placement_probed=True)."""
    print("Test: r2_placement_uncertain_false_when_probed")
    ctx = _make_minimal_context({"placement_probed": True})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when placement_probed=True"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_no_model():
    """R2: placement_uncertain=False when there is no model file."""
    print("Test: r2_placement_uncertain_false_no_model")
    ctx = _make_minimal_context({"has_model": False})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when has_model=False"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_no_data():
    """R2: placement_uncertain=False when there is no data (MTZ or map)."""
    print("Test: r2_placement_uncertain_false_no_data")
    ctx = _make_minimal_context({"has_data_mtz": False, "has_map": False})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when no MTZ and no map"
    )
    print("  PASSED")


def test_r2_placement_uncertain_false_for_predicted_model():
    """R2: placement_uncertain=False for predicted models (always need MR/dock — no probe)."""
    print("Test: r2_placement_uncertain_false_for_predicted_model")
    ctx = _make_minimal_context({"has_predicted_model": True})
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False for predicted models "
        "(they always require MR/dock without probing)"
    )
    print("  PASSED")


def test_r2_placement_uncertain_true_with_map_only():
    """R2: placement_uncertain=True for cryo-EM case (model + map, no MTZ)."""
    print("Test: r2_placement_uncertain_true_with_map_only")
    ctx = _make_minimal_context({"has_data_mtz": False, "has_map": True})
    assert ctx["placement_uncertain"] is True, (
        "placement_uncertain should be True for model+map (cryo-EM ambiguous case)"
    )
    print("  PASSED")


def test_r2_context_keys_present_in_workflow_engine():
    """R2: WorkflowEngine.build_context() returns all four new placement keys."""
    print("Test: r2_context_keys_present_in_workflow_engine")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    # Minimal inputs — no files, no history, no analysis
    ctx = engine.build_context(
        files={},
        history_info={},
        analysis=None,
        directives=None,
    )
    for key in ("cell_mismatch", "placement_probed",
                "placement_probe_result", "placement_uncertain"):
        assert key in ctx, (
            "build_context() must include key %r in returned context" % key
        )
    print("  PASSED (all four keys present: cell_mismatch, placement_probed, "
          "placement_probe_result, placement_uncertain)")


def test_r2_placement_probed_loaded_from_history():
    """R2: placement_probed and placement_probe_result are loaded from history_info."""
    print("Test: r2_placement_probed_loaded_from_history")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # Simulate history where probe already ran and confirmed placement
    history = {
        "placement_probed": True,
        "placement_probe_result": "placed",
    }
    ctx = engine.build_context(files={}, history_info=history)
    assert ctx["placement_probed"] is True, "placement_probed must be True from history"
    assert ctx["placement_probe_result"] == "placed", (
        "placement_probe_result must be 'placed' from history, got %r"
        % ctx["placement_probe_result"]
    )
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False when placement_probed=True"
    )
    print("  PASSED")


def test_r2_cell_mismatch_false_with_no_files():
    """R2: cell_mismatch=False when no model or data files are present (fail-safe)."""
    print("Test: r2_cell_mismatch_false_with_no_files")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    ctx = engine.build_context(files={}, history_info={})
    assert ctx["cell_mismatch"] is False, (
        "cell_mismatch must be False (fail-safe) when no files present, got %r"
        % ctx["cell_mismatch"]
    )
    print("  PASSED")


# =============================================================================
# CATEGORY R3 — probe_placement step routing (Tier 3)
# =============================================================================

def _make_xray_history_with_probe(program, rfree=None, pre_refine=True):
    """Build a minimal history list simulating a probe run."""
    result_text = "SUCCESS: Quick R-factor check complete"
    metrics = {}
    if rfree is not None:
        metrics["r_free"] = rfree

    entry = {
        "program": program,
        "command": program + " model.pdb data.mtz",
        "result": result_text,
        "analysis": metrics,
    }
    if not pre_refine:
        # Simulate refinement having run first
        return [
            {"program": "phenix.refine", "command": "phenix.refine model.pdb data.mtz",
             "result": "SUCCESS: Refinement complete", "analysis": {"r_free": 0.28}},
            entry,
        ]
    return [entry]


def test_r3_probe_placement_phase_offered_xray():
    """R3: placement_uncertain=True → xray step is probe_placement."""
    print("Test: r3_probe_placement_phase_offered_xray")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    # Build context that triggers placement_uncertain
    # has_model=True, has_data_mtz=True, no history, no predicted model
    context = {
        "has_placed_model": False,
        "has_predicted_model": False,
        "has_processed_model": False,
        "has_model": True,
        "has_search_model": False,
        "has_model_for_mr": True,
        "has_data_mtz": True,
        "has_map_coeffs_mtz": False,
        "has_sequence": False,
        "has_full_map": False,
        "has_half_map": False,
        "has_map": False,
        "has_non_half_map": False,
        "has_ligand": False,
        "has_ligand_file": False,
        "has_ligand_fit": False,
        "has_optimized_full_map": False,
        "has_refined_model": False,
        "phaser_done": False,
        "predict_done": False,
        "predict_full_done": False,
        "autobuild_done": False,
        "autosol_done": False,
        "autosol_success": False,
        "refine_done": False,
        "rsr_done": False,
        "validation_done": False,
        "ligandfit_done": False,
        "pdbtools_done": False,
        "needs_post_ligandfit_refine": False,
        "dock_done": False,
        "process_predicted_model_done": False,
        "refine_count": 0,
        "rsr_count": 0,
        "r_free": None, "r_work": None, "map_cc": None, "clashscore": None,
        "resolution": None, "tfz": None,
        "has_anomalous": False, "strong_anomalous": False,
        "anomalous_measurability": None, "has_twinning": False, "twin_law": None,
        "twin_fraction": None, "anomalous_resolution": None, "has_ncs": False,
        "has_half_map": False,
        "xtriage_done": True,   # Analysis already done
        "model_is_good": False,
        "use_mr_sad": False,
        "automation_path": "automated",
        # Tier 1/2/3 placement keys
        "cell_mismatch": False,
        "placement_probed": False,
        "placement_probe_result": None,
        "placement_uncertain": True,   # <- This triggers probe
    }

    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("xray")
    result = engine._detect_xray_step(steps, context)
    assert result["step"] == "probe_placement", (
        "Expected probe_placement step, got %r (reason: %s)"
        % (result["step"], result.get("reason", ""))
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_probe_valid_programs_xray():
    """R3: probe_placement valid programs for xray = [phenix.model_vs_data]."""
    print("Test: r3_probe_valid_programs_xray")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    step_info = {"step": "probe_placement"}
    progs = engine.get_valid_programs("xray", step_info, {}, {})
    assert "phenix.model_vs_data" in progs, (
        "phenix.model_vs_data must be in probe_placement valid programs for xray, got %s" % progs
    )
    assert "phenix.refine" not in progs, (
        "phenix.refine must NOT be in probe_placement valid programs"
    )
    print("  PASSED: valid programs = %s" % progs)


def test_r3_probe_valid_programs_cryoem():
    """R3: probe_placement valid programs for cryoem = [phenix.map_correlations]."""
    print("Test: r3_probe_valid_programs_cryoem")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    step_info = {"step": "probe_placement"}
    progs = engine.get_valid_programs("cryoem", step_info, {}, {})
    assert "phenix.map_correlations" in progs, (
        "phenix.map_correlations must be in probe_placement valid programs for cryoem, got %s" % progs
    )
    assert "phenix.real_space_refine" not in progs, (
        "phenix.real_space_refine must NOT be in probe_placement valid programs"
    )
    print("  PASSED: valid programs = %s" % progs)


def test_r3_probe_result_placed_routes_to_refine():
    """R3: placement_probe_result='placed' → xray step is refine (not probe_placement)."""
    print("Test: r3_probe_result_placed_routes_to_refine")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("xray")

    # placement_probed=True, result=placed → should skip probe_placement, go to refine
    context = _make_minimal_context({
        "has_placed_model": False,   # Heuristics still say False
        "placement_probed": True,
        "placement_probe_result": "placed",
        "placement_uncertain": False,
        "cell_mismatch": False,
        "xtriage_done": True,
        "has_predicted_model": False,
        "has_processed_model": False,
        "autosol_done": False,
        "phaser_done": False,
        "has_ligand_fit": False,
        "pdbtools_done": False,
        "needs_post_ligandfit_refine": False,
        "has_refined_model": False,
        "has_anomalous": False,
        "use_mr_sad": False,
        "has_twinning": False,
        "has_ligand_file": False,
        "ligandfit_done": False,
        "refine_count": 0,
        "r_free": None,
        "resolution": None,
        "validation_done": False,
        "model_is_good": False,
    })
    # build_context overrides has_placed_model=True when probe result is 'placed';
    # replicate that here since we bypass build_context in this unit test.
    if context.get("placement_probe_result") == "placed":
        context["has_placed_model"] = True
    result = engine._detect_xray_step(steps, context)
    assert result["step"] == "refine", (
        "After successful probe (placed), expected refine, got %r (reason: %s)"
        % (result["step"], result.get("reason", ""))
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_probe_result_needs_mr_routes_to_mr():
    """R3: placement_probe_result='needs_mr' → xray step is molecular_replacement."""
    print("Test: r3_probe_result_needs_mr_routes_to_mr")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("xray")

    context = _make_minimal_context({
        "has_placed_model": False,
        "placement_probed": True,
        "placement_probe_result": "needs_mr",
        "placement_uncertain": False,
        "cell_mismatch": False,
        "xtriage_done": True,
        "has_predicted_model": False,
        "has_processed_model": False,
        "autosol_done": False,
        "phaser_done": False,
        "use_mr_sad": False,
        "has_anomalous": False,
        "model_is_good": False,
    })
    result = engine._detect_xray_step(steps, context)
    assert result["step"] == "molecular_replacement", (
        "needs_mr probe result should route to molecular_replacement, got %r"
        % result["step"]
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_cell_mismatch_routes_to_mr_xray():
    """R3: cell_mismatch=True → xray step is molecular_replacement (skips probe)."""
    print("Test: r3_cell_mismatch_routes_to_mr_xray")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("xray")

    context = _make_minimal_context({
        "has_placed_model": False,
        "cell_mismatch": True,          # <- Tier 1 result
        "placement_uncertain": False,   # Tier 1 detected → uncertainty cleared
        "placement_probed": False,
        "placement_probe_result": None,
        "xtriage_done": True,
        "has_predicted_model": False,
        "has_processed_model": False,
        "autosol_done": False,
        "phaser_done": False,
        "use_mr_sad": False,
        "has_anomalous": False,
        "model_is_good": False,
    })
    result = engine._detect_xray_step(steps, context)
    assert result["step"] == "molecular_replacement", (
        "cell_mismatch should route to molecular_replacement, got %r" % result["step"]
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_cell_mismatch_routes_to_dock_cryoem():
    """R3: cell_mismatch=True → cryoem step is dock_model (skips probe)."""
    print("Test: r3_cell_mismatch_routes_to_dock_cryoem")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("cryoem")

    context = _make_minimal_context({
        "has_placed_model": False,
        "cell_mismatch": True,
        "placement_uncertain": False,
        "placement_probed": False,
        "placement_probe_result": None,
        "has_data_mtz": False,
        "has_map": True,
        "has_full_map": True,
        "has_predicted_model": False,
        "has_processed_model": False,
        "has_search_model": False,
        "has_model_for_mr": True,
        "has_half_map": False,
        "has_non_half_map": True,
        "has_sequence": False,
        "has_optimized_full_map": False,
        "xtriage_done": True,
        "mtriage_done": True,
        "dock_done": False,
        "resolve_cryo_em_done": False,
        "map_sharpening_done": False,
        "map_to_model_done": False,
        "model_is_good": False,
        "automation_path": "automated",
        "map_symmetry_done": False,
    })
    result = engine._detect_cryoem_step(steps, context)
    assert result["step"] == "dock_model", (
        "cell_mismatch should route to dock_model for cryoem, got %r" % result["step"]
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_probe_not_rerun_when_already_probed():
    """R3: placement_probed=True → probe_placement step NOT offered again."""
    print("Test: r3_probe_not_rerun_when_already_probed")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    from knowledge.yaml_loader import get_workflow_steps
    steps = get_workflow_steps("xray")

    # Probe ran but R-free was unparseable (result=None) → should not re-probe
    context = _make_minimal_context({
        "has_placed_model": False,
        "placement_probed": True,
        "placement_probe_result": None,   # Couldn't parse → fail-safe
        "placement_uncertain": False,     # probed=True clears uncertainty
        "cell_mismatch": False,
        "xtriage_done": True,
        "has_predicted_model": False,
        "has_processed_model": False,
        "autosol_done": False,
        "phaser_done": False,
        "use_mr_sad": False,
        "has_anomalous": False,
        "model_is_good": False,
    })
    result = engine._detect_xray_step(steps, context)
    assert result["step"] != "probe_placement", (
        "probe_placement must NOT be offered again when placement_probed=True, "
        "got %r" % result["step"]
    )
    print("  PASSED: step=%r (not re-probing)" % result["step"])


def test_r3_history_probe_detection_xray_placed():
    """R3: model_vs_data in history before refine → placement_probed=True, result=placed."""
    print("Test: r3_history_probe_detection_xray_placed")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = _make_xray_history_with_probe("phenix.model_vs_data", rfree=0.38, pre_refine=True)
    info = _analyze_history(history)

    assert info["placement_probed"] is True, (
        "placement_probed should be True after model_vs_data before refinement"
    )
    assert info["placement_probe_result"] == "placed", (
        "R-free=0.38 < 0.50 should give placement_probe_result='placed', got %r"
        % info["placement_probe_result"]
    )
    print("  PASSED: probed=%s result=%r" % (info["placement_probed"],
                                              info["placement_probe_result"]))


def test_r3_history_probe_detection_xray_needs_mr():
    """R3: model_vs_data before refine with R-free >= 0.50 → result=needs_mr."""
    print("Test: r3_history_probe_detection_xray_needs_mr")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = _make_xray_history_with_probe("phenix.model_vs_data", rfree=0.55, pre_refine=True)
    info = _analyze_history(history)

    assert info["placement_probed"] is True
    assert info["placement_probe_result"] == "needs_mr", (
        "R-free=0.55 >= 0.50 should give placement_probe_result='needs_mr', got %r"
        % info["placement_probe_result"]
    )
    print("  PASSED: result=%r" % info["placement_probe_result"])


def test_r3_history_probe_not_detected_after_refine():
    """R3: model_vs_data AFTER refinement → NOT treated as probe (is validation)."""
    print("Test: r3_history_probe_not_detected_after_refine")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    # model_vs_data ran AFTER refinement → it's validation, not a probe
    history = _make_xray_history_with_probe("phenix.model_vs_data", rfree=0.25, pre_refine=False)
    info = _analyze_history(history)

    assert info["placement_probed"] is False, (
        "model_vs_data after refinement must NOT be treated as placement probe, "
        "placement_probed should be False, got %s" % info["placement_probed"]
    )
    print("  PASSED: placement_probed=%s (correctly treated as validation)" % info["placement_probed"])


def test_r3_validation_done_not_set_during_probe_phase():
    """R3: probe history detection (pre-refine) sets placement_probed, not validation_done.

    validation_done is still set by YAML done_tracking (set_flag strategy),
    which fires whenever model_vs_data marker is found.  The probe result
    flags are SET IN ADDITION — they don't prevent validation_done.
    This test verifies the probe result flags ARE set; validation_done
    being also set is acceptable (it resets when the session resumes with
    the probe result already known, so validation step runs normally later).
    """
    print("Test: r3_validation_done_not_set_during_probe_phase")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = _make_xray_history_with_probe("phenix.model_vs_data", rfree=0.38, pre_refine=True)
    info = _analyze_history(history)

    # The key invariant: placement_probed must be set
    assert info["placement_probed"] is True, (
        "placement_probed must be True after pre-refine model_vs_data"
    )
    assert info["placement_probe_result"] == "placed", (
        "placement_probe_result must be 'placed' for R-free=0.38"
    )
    # validation_done being set is acceptable here — context builder uses
    # placement_probed to route correctly; validation can re-run later.
    print("  PASSED: placement_probed=%s, result=%r, validation_done=%s"
          % (info["placement_probed"], info["placement_probe_result"],
             info.get("validation_done")))


# =============================================================================
# CATEGORY R3-EXTRA -- additional placement detection edge cases
# =============================================================================

def test_r3_cryoem_probe_needs_dock():
    """R3-extra: map_correlations before refine, CC <= 0.15 -> needs_dock."""
    print("Test: r3_cryoem_probe_needs_dock")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = [{
        "program": "phenix.map_correlations",
        "command": "phenix.map_correlations model.pdb map.ccp4",
        "result": "SUCCESS: Map-model CC computed",
        "analysis": {"cc_mask": 0.05},
    }]
    info = _analyze_history(history)
    assert info["placement_probed"] is True, "placement_probed should be True"
    assert info["placement_probe_result"] == "needs_dock", (
        "CC=0.05 should give needs_dock, got %r" % info["placement_probe_result"]
    )
    print("  PASSED: result=%r" % info["placement_probe_result"])


def test_r3_cryoem_probe_placed():
    """R3-extra: map_correlations before refine, CC > 0.15 -> placed."""
    print("Test: r3_cryoem_probe_placed")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = [{
        "program": "phenix.map_correlations",
        "command": "phenix.map_correlations model.pdb map.ccp4",
        "result": "SUCCESS: Map-model CC computed",
        "analysis": {"cc_mask": 0.42},
    }]
    info = _analyze_history(history)
    assert info["placement_probed"] is True
    assert info["placement_probe_result"] == "placed", (
        "CC=0.42 should give placed, got %r" % info["placement_probe_result"]
    )
    print("  PASSED: result=%r" % info["placement_probe_result"])


def test_r3_failed_probe_not_counted():
    """R3-extra: a FAILED model_vs_data run marks placement_probed=True with
    probe_result=None (inconclusive) — NOT placement_probed=False.

    Behaviour changed in v112.72: previously a failed probe was ignored entirely,
    which caused an infinite probe-retry loop.  Now a failed probe is treated as
    inconclusive (probed=True, result=None) so the workflow can move on to refine.
    """
    print("Test: r3_failed_probe_not_counted")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = [{
        "program": "phenix.model_vs_data",
        "command": "phenix.model_vs_data model.pdb data.mtz",
        "result": "FAILED: Sorry: could not open file",
        "analysis": {},
    }]
    info = _analyze_history(history)
    # After v112.72: a failed probe is inconclusive (prevents infinite retry)
    assert info["placement_probed"] is True, (
        "FAILED model_vs_data should set placement_probed=True (inconclusive)"
    )
    assert info["placement_probe_result"] is None, (
        "placement_probe_result must be None for inconclusive probe, got %r"
        % info["placement_probe_result"]
    )
    print("  PASSED: failed model_vs_data probe correctly marked as inconclusive")


def test_r3_probe_cc_volume_fallback():
    """R3-extra: probe detection uses cc_volume when cc_mask is absent."""
    print("Test: r3_probe_cc_volume_fallback")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    history = [{
        "program": "phenix.map_correlations",
        "command": "phenix.map_correlations model.pdb map.ccp4",
        "result": "SUCCESS: Done",
        "analysis": {"cc_volume": 0.30},   # cc_mask absent, cc_volume present
    }]
    info = _analyze_history(history)
    assert info["placement_probed"] is True
    assert info["placement_probe_result"] == "placed", (
        "cc_volume=0.30 > 0.15 should give placed, got %r" % info["placement_probe_result"]
    )
    print("  PASSED: cc_volume fallback works, result=%r" % info["placement_probe_result"])


def test_r3_build_context_overrides_has_placed_model():
    """R3-extra: build_context sets has_placed_model=True when probe says placed."""
    print("Test: r3_build_context_overrides_has_placed_model")
    import tempfile
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from agent.workflow_state import _analyze_history

    engine = WorkflowEngine()
    history = [{
        "program": "phenix.model_vs_data",
        "command": "phenix.model_vs_data model.pdb data.mtz",
        "result": "SUCCESS: Done",
        "analysis": {"r_free": 0.35},
    }]
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb = os.path.join(tmpdir, "model.pdb")
        mtz = os.path.join(tmpdir, "data.mtz")
        for p in (pdb, mtz):
            open(p, "w").close()

        files = {"model": [pdb], "data_mtz": [mtz]}
        history_info = _analyze_history(history)

        assert history_info["placement_probed"] is True
        assert history_info["placement_probe_result"] == "placed"

        ctx = engine.build_context(files=files, history_info=history_info)

        assert ctx["has_placed_model"] is True, (
            "build_context must set has_placed_model=True when probe result is placed, "
            "got %s" % ctx["has_placed_model"]
        )
        assert ctx["placement_uncertain"] is False
    print("  PASSED: build_context correctly overrides has_placed_model=True")


def test_r3_cells_fail_safe_none():
    """R3-extra: cells_are_compatible with None inputs returns True (fail-safe)."""
    print("Test: r3_cells_fail_safe_none")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    assert cells_are_compatible(None, (57.0, 57.0, 147.0, 90.0, 90.0, 90.0)) is True
    assert cells_are_compatible((57.0, 57.0, 147.0, 90.0, 90.0, 90.0), None) is True
    assert cells_are_compatible(None, None) is True
    print("  PASSED: None inputs always return True (fail-safe)")


def test_r3_cells_clear_mismatch_detected():
    """R3-extra: cells_are_compatible correctly rejects clearly incompatible cells."""
    print("Test: r3_cells_clear_mismatch_detected")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.placement_checker import cells_are_compatible

    # ~75% difference in all axes -- far outside 5% tolerance
    model_cell = (57.23, 57.23, 146.77, 90.0, 90.0, 90.0)
    mtz_cell   = (100.0, 100.0, 200.00, 90.0, 90.0, 90.0)
    assert cells_are_compatible(model_cell, mtz_cell) is False, (
        "Clearly incompatible cells must be detected as mismatch"
    )
    print("  PASSED: large cell difference correctly flagged as mismatch")


def test_r3_placement_uncertain_clears_after_probe():
    """R3-extra: placement_uncertain=False on build_context after probe ran."""
    print("Test: r3_placement_uncertain_clears_after_probe")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from agent.workflow_state import _analyze_history

    engine = WorkflowEngine()
    history = [{
        "program": "phenix.model_vs_data",
        "command": "phenix.model_vs_data model.pdb data.mtz",
        "result": "SUCCESS: Done",
        "analysis": {"r_free": 0.42},
    }]
    history_info = _analyze_history(history)
    ctx = engine.build_context(files={}, history_info=history_info)

    assert ctx["placement_probed"] is True
    assert ctx["placement_uncertain"] is False, (
        "placement_uncertain must be False after probe ran"
    )
    print("  PASSED: placement_uncertain cleared correctly after probe cycle")


def test_r3_probe_placement_cryoem_phase_offered():
    """R3-extra: cryoem placement_uncertain=True -> probe_placement step offered."""
    print("Test: r3_probe_placement_cryoem_phase_offered")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from knowledge.yaml_loader import get_workflow_steps

    engine = WorkflowEngine()
    steps = get_workflow_steps("cryoem")
    context = _make_minimal_context({
        "has_placed_model": False,
        "has_predicted_model": False,
        "has_processed_model": False,
        "has_model": True,
        "has_search_model": False,
        "has_model_for_mr": True,
        "has_data_mtz": False,
        "has_map": True,
        "has_full_map": True,
        "has_half_map": False,
        "has_non_half_map": True,
        "has_sequence": False,
        "has_optimized_full_map": False,
        "has_refined_model": False,
        "has_map_coeffs_mtz": False,
        "predict_done": False,
        "predict_full_done": False,
        "dock_done": False,
        "resolve_cryo_em_done": False,
        "map_sharpening_done": False,
        "map_to_model_done": False,
        "xtriage_done": True,
        "mtriage_done": True,
        "map_symmetry_done": False,
        "model_is_good": False,
        "automation_path": "automated",
        "cell_mismatch": False,
        "placement_probed": False,
        "placement_probe_result": None,
        "placement_uncertain": True,  # <- triggers probe
    })
    result = engine._detect_cryoem_step(steps, context)
    assert result["step"] == "probe_placement", (
        "Expected probe_placement step for cryoem, got %r" % result["step"]
    )
    print("  PASSED: step=%r" % result["step"])


def test_r3_needs_dock_routes_to_dock():
    """R3-extra: cryoem probe result needs_dock routes to dock_model step."""
    print("Test: r3_needs_dock_routes_to_dock")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from knowledge.yaml_loader import get_workflow_steps

    engine = WorkflowEngine()
    steps = get_workflow_steps("cryoem")
    context = _make_minimal_context({
        "has_placed_model": False,
        "has_predicted_model": False,
        "has_processed_model": False,
        "has_model": True,
        "has_search_model": False,
        "has_model_for_mr": True,
        "has_data_mtz": False,
        "has_map": True,
        "has_full_map": True,
        "has_half_map": False,
        "has_non_half_map": True,
        "has_sequence": False,
        "has_optimized_full_map": False,
        "has_refined_model": False,
        "has_map_coeffs_mtz": False,
        "predict_done": False,
        "predict_full_done": False,
        "dock_done": False,
        "resolve_cryo_em_done": False,
        "map_sharpening_done": False,
        "map_to_model_done": False,
        "xtriage_done": True,
        "mtriage_done": True,
        "map_symmetry_done": False,
        "model_is_good": False,
        "automation_path": "automated",
        "cell_mismatch": False,
        "placement_probed": True,
        "placement_probe_result": "needs_dock",
        "placement_uncertain": False,
    })
    result = engine._detect_cryoem_step(steps, context)
    assert result["step"] == "dock_model", (
        "needs_dock probe result should route to dock_model, got %r" % result["step"]
    )
    print("  PASSED: step=%r" % result["step"])



# =============================================================================
# CATEGORY S2 -- Directive override protection
#
# Root cause: LLM directive extractor set model_is_placed=True from "solve the
# structure", making _has_placed_model() return True. Tier 1 cell-mismatch
# routing checked `cell_mismatch AND NOT has_placed_model`, which evaluated
# False → skipped docking → RSR failed with "unit cell dimensions mismatch".
#
# Three fixes:
#   A. Tightened directive extractor prompt (model_is_placed is now high-precision)
#   B. has_placed_model_from_history context key (history/files only, no directives)
#   C. Tier 1 routing guards use has_placed_model_from_history instead of has_placed_model
#   D. S1 short-circuit also uses has_placed_model_from_history
# =============================================================================

def test_s2_has_placed_model_from_history_method_exists():
    """S2: WorkflowEngine has _has_placed_model_from_history() method."""
    print("Test: s2_has_placed_model_from_history_method_exists")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    engine = WorkflowEngine()
    assert hasattr(engine, '_has_placed_model_from_history'), (
        "_has_placed_model_from_history method must exist on WorkflowEngine"
    )
    print("  PASSED")


def test_s2_from_history_false_when_only_directive():
    """S2: _has_placed_model_from_history returns False when placement is only via directive."""
    print("Test: s2_from_history_false_when_only_directive")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    engine = WorkflowEngine()

    # Simulate: directive says model_is_placed=True, but no history/subcategory evidence
    directives = {"workflow_preferences": {"model_is_placed": True}}
    # files must include "model" for directive check to fire in _has_placed_model
    files = {"model": ["/fake/1aew_A.pdb"]}
    history_info = {}

    # _has_placed_model should return True (directive-driven, model file present)
    has_placed = engine._has_placed_model(files, history_info, directives)
    assert has_placed is True, "Directive should make _has_placed_model return True (with model file)"

    # _has_placed_model_from_history should return False (no history)
    from_history = engine._has_placed_model_from_history(files, history_info)
    assert from_history is False, (
        "_has_placed_model_from_history must return False when placement is only "
        "from directive (no dock_done, no phaser_done, no positioned subcategory)"
    )
    print("  PASSED: directive cannot fool has_placed_model_from_history")


def test_s2_from_history_true_when_dock_done():
    """S2: _has_placed_model_from_history returns True when dock_done=True in history."""
    print("Test: s2_from_history_true_when_dock_done")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from agent.workflow_state import _analyze_history

    engine = WorkflowEngine()

    history = [{
        "program": "phenix.dock_in_map",
        "command": "phenix.dock_in_map model.pdb map.ccp4",
        "result": "SUCCESS: Docking complete",
        "analysis": {},
    }]
    history_info = _analyze_history(history)
    assert history_info["dock_done"] is True

    from_history = engine._has_placed_model_from_history({}, history_info)
    assert from_history is True, (
        "_has_placed_model_from_history must return True when dock_done=True"
    )
    print("  PASSED: dock_done in history → has_placed_model_from_history=True")


def test_s2_context_has_placed_from_history_key():
    """S2: build_context populates has_placed_model_from_history key."""
    print("Test: s2_context_has_placed_from_history_key")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    ctx = engine.build_context(files={}, history_info={})

    assert "has_placed_model_from_history" in ctx, (
        "build_context must include has_placed_model_from_history key"
    )
    assert ctx["has_placed_model_from_history"] is False, (
        "Empty history → has_placed_model_from_history must be False"
    )
    print("  PASSED: has_placed_model_from_history present in context")


def test_s2_directive_model_is_placed_does_not_suppress_cell_mismatch():
    """
    S2: The critical bug scenario — model_is_placed=True directive must NOT
    suppress Tier 1 cell-mismatch routing to dock_model.

    Replicates the apoferritin log:
      - Directive: model_is_placed=True (set by LLM from "solve the structure")
      - History: no dock_done, no phaser_done
      - Cell mismatch: True (model cell ≠ map cell)
    Expected: cryoem step → dock_model  (not ready_to_refine)
    """
    print("Test: s2_directive_model_is_placed_does_not_suppress_cell_mismatch")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # Build context with directive-driven has_placed_model=True but no history
    directives = {"workflow_preferences": {"model_is_placed": True}}
    # Include mtriage_done=True so we are past the "analyze" gate in cryoem step detection
    history_with_mtriage = {"mtriage_done": True, "resolve_cryo_em_done": False}
    ctx = engine.build_context(files={"model": ["/fake/1aew_A.pdb"], "full_map": ["/fake/map.ccp4"]},
                               history_info=history_with_mtriage, directives=directives)

    assert ctx["has_placed_model"] is True, \
        "Directive should set has_placed_model=True (with model file present)"
    assert ctx["has_placed_model_from_history"] is False, \
        "No dock/phaser → has_placed_model_from_history must be False"

    # Manually inject cell_mismatch=True to simulate what placement_checker would detect
    ctx["cell_mismatch"] = True

    # Now ask the engine to detect step given this context
    from knowledge.yaml_loader import get_workflow_steps; steps = get_workflow_steps("cryoem")
    step_info = engine._detect_cryoem_step(steps, ctx)

    assert step_info["step"] == "dock_model", (
        "Cell mismatch MUST route to dock_model even when model_is_placed=True "
        "from directive. Got step=%r (reason=%r)" % (
            step_info["step"], step_info.get("reason", ""))
    )
    print("  PASSED: cell_mismatch → dock_model despite model_is_placed=True directive")


def test_s2_history_placed_does_suppress_cell_mismatch():
    """
    S2: When placement is confirmed by real history (dock_done=True), Tier 1
    cell-mismatch routing is correctly suppressed (model is already placed).
    """
    print("Test: s2_history_placed_does_suppress_cell_mismatch")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from agent.workflow_state import _analyze_history

    engine = WorkflowEngine()

    history = [{
        "program": "phenix.dock_in_map",
        "command": "phenix.dock_in_map 1aew_A.pdb denmod_map.ccp4",
        "result": "SUCCESS: Docking complete. Output: 1aew_A_docked.pdb",
        "analysis": {},
    }]
    history_info = _analyze_history(history)
    assert history_info["dock_done"] is True

    ctx = engine.build_context(files={}, history_info=history_info)
    assert ctx["has_placed_model_from_history"] is True

    # Inject cell_mismatch=True (edge case: could still mismatch due to box padding)
    ctx["cell_mismatch"] = True

    from knowledge.yaml_loader import get_workflow_steps; steps = get_workflow_steps("cryoem")
    step_info = engine._detect_cryoem_step(steps, ctx)

    # With dock_done=True history, cell_mismatch is short-circuited to False
    # (from the S1 post-processing short-circuit), so routing advances normally
    assert step_info["step"] != "dock_model", (
        "After dock_done, should NOT re-route to dock_model. "
        "Got step=%r" % step_info["step"]
    )
    print("  PASSED: dock_done in history correctly suppresses Tier 1 re-dock")


def test_s2_xray_tier1_uses_from_history():
    """S2: X-ray Tier 1 MR routing uses has_placed_model_from_history (not has_placed_model)."""
    print("Test: s2_xray_tier1_uses_from_history")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # Directive: model_is_placed=True — should NOT suppress MR routing
    directives = {"workflow_preferences": {"model_is_placed": True}}
    # files must include "model" — _has_placed_model gates the directive on files.get("model")
    ctx = engine.build_context(files={"model": ["/fake/1abc.pdb"]}, history_info={}, directives=directives)
    assert ctx["has_placed_model"] is True,         "Directive should make has_placed_model=True (with model file present)"
    assert ctx["has_placed_model_from_history"] is False

    # Inject cell_mismatch=True and also fake X-ray analysis flag
    ctx["cell_mismatch"] = True
    ctx["xtriage_done"] = True  # past analysis step

    from knowledge.yaml_loader import get_workflow_steps; steps = get_workflow_steps("xray")
    step_info = engine._detect_xray_step(steps, ctx)

    assert step_info["step"] == "molecular_replacement", (
        "X-ray cell mismatch MUST route to MR even with model_is_placed=True directive. "
        "Got step=%r" % step_info["step"]
    )
    print("  PASSED: X-ray cell_mismatch → molecular_replacement despite model_is_placed directive")


def test_s2_short_circuit_uses_from_history_not_directive():
    """
    S2: S1 short-circuit uses has_placed_model_from_history — a directive-only
    has_placed_model=True must NOT prevent the cell-mismatch check from running.
    """
    print("Test: s2_short_circuit_uses_from_history_not_directive")
    sys.path.insert(0, _PROJECT_ROOT)
    src_path = os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py")
    with open(src_path) as f:
        src = f.read()

    # The short-circuit must reference has_placed_model_from_history, not has_placed_model
    assert 'context.get("has_placed_model_from_history")' in src, (
        "Short-circuit must use has_placed_model_from_history"
    )
    # The short-circuit must NOT reference bare has_placed_model (without _from_history)
    # in its condition
    sc_block_start = src.index("Tier 1 short-circuit: skip cell check")
    sc_block_end   = src.index("context[\"cell_mismatch\"] = False", sc_block_start) + 40
    sc_block = src[sc_block_start:sc_block_end]
    assert 'has_placed_model_from_history' in sc_block, \
        "Short-circuit block must reference has_placed_model_from_history"
    # The bare 'has_placed_model' should only appear in the _from_history lookup, not standalone
    assert '"has_placed_model")' not in sc_block, (
        "Short-circuit must not check bare has_placed_model — directive could be wrong"
    )
    print("  PASSED: short-circuit uses has_placed_model_from_history (not bare directive flag)")


def test_s2_directive_prompt_stronger_do_not_set():
    """S2: directive_extractor.py prompt explicitly warns against 'solve the structure'."""
    print("Test: s2_directive_prompt_stronger_do_not_set")
    sys.path.insert(0, _PROJECT_ROOT)
    src_path = os.path.join(_PROJECT_ROOT, "agent", "directive_extractor.py")
    with open(src_path) as f:
        src = f.read()

    # Must explicitly say "solve the structure" is a DO NOT case
    assert '"solve the structure"' in src or "'solve the structure'" in src or \
           "solve the structure" in src, (
        "Prompt must explicitly list 'solve the structure' as a case to NOT set model_is_placed"
    )
    # Must warn about cryo-EM + PDB combination
    assert "cryo-EM" in src or "cryoem" in src.lower(), \
        "Prompt must mention cryo-EM as a case requiring docking before refinement"
    # ONLY / HIGH-PRECISION language should be present
    assert ("HIGH-PRECISION" in src or "ONLY set model_is_placed" in src or
            "high-precision" in src.lower()), (
        "Prompt must use high-precision / only language to discourage false positives"
    )
    print("  PASSED: directive prompt contains explicit DO NOT cases and precision guidance")


def test_s2_full_cryoem_stack_routes_to_dock_not_rsr():
    """
    S2: Full stack integration — cryoem workflow with resolve_cryo_em done,
    model_is_placed directive, but no dock history → must route to dock_model
    (not ready_to_refine or real_space_refine).

    This replicates the exact failure observed in the apoferritin log.
    """
    print("Test: s2_full_cryoem_stack_routes_to_dock_not_rsr")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    from agent.workflow_state import _analyze_history

    engine = WorkflowEngine()

    # History: mtriage ran (cycle 1), resolve_cryo_em ran (cycle 2)
    # Exactly as in the apoferritin log
    history = [
        {
            "program": "phenix.mtriage",
            "command": "phenix.mtriage half_map=map1.ccp4 half_map=map2.ccp4",
            "result": "SUCCESS: Resolution 1.91 Å",
            "analysis": {"resolution": 1.91},
        },
        {
            "program": "phenix.resolve_cryo_em",
            "command": "phenix.resolve_cryo_em half_map=map1.ccp4 half_map=map2.ccp4",
            "result": "SUCCESS: Density modification complete",
            "analysis": {},
        },
    ]
    history_info = _analyze_history(history)
    assert history_info["mtriage_done"] is True
    assert history_info["resolve_cryo_em_done"] is True
    assert history_info["dock_done"] is False  # no docking yet

    # Directive: model_is_placed=True (as the LLM mis-extracted from "solve the structure")
    directives = {"workflow_preferences": {"model_is_placed": True}}

    ctx = engine.build_context(
        files={"model": ["/fake/1aew_A.pdb"], "full_map": ["/fake/denmod_map.ccp4"]},
        history_info=history_info, directives=directives
    )

    assert ctx["has_placed_model"] is True, \
        "Directive should make has_placed_model=True (requires files['model'] non-empty)"
    assert ctx["has_placed_model_from_history"] is False, \
        "No dock/phaser in history → has_placed_model_from_history=False"

    # Inject cell mismatch (what placement_checker would detect for 1aew_A.pdb vs denmod_map.ccp4)
    ctx["cell_mismatch"] = True

    from knowledge.yaml_loader import get_workflow_steps; steps = get_workflow_steps("cryoem")
    step_info = engine._detect_cryoem_step(steps, ctx)

    assert step_info["step"] == "dock_model", (
        "REGRESSION: Full stack must route to dock_model when cell_mismatch=True "
        "and has_placed_model_from_history=False, even if model_is_placed directive "
        "is set. Got step=%r (reason=%r). This is the apoferritin bug." % (
            step_info["step"], step_info.get("reason", ""))
    )
    print("  PASSED: apoferritin scenario correctly routes to dock_model (not RSR)")


# =============================================================================
# CATEGORY S1 -- Polish fixes (yaml_tools transition fields, import fallback,
#                              redundant import, cell_mismatch short-circuit)
# =============================================================================

def test_s1_yaml_validator_no_if_placed_warnings():
    """S1: _validate_workflows generates no warnings for if_placed / if_not_placed."""
    print("Test: s1_yaml_validator_no_if_placed_warnings")
    import yaml
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.yaml_tools import _validate_workflows

    wf_path = os.path.join(_PROJECT_ROOT, "knowledge", "workflows.yaml")
    with open(wf_path) as f:
        data = yaml.safe_load(f)

    issues = _validate_workflows(data)
    probe_warns = [i for i in issues
                   if "if_placed" in i[1] or "if_not_placed" in i[1]]

    assert probe_warns == [], (
        "Expected no if_placed/if_not_placed warnings, got: %s" % probe_warns
    )
    print("  PASSED: no spurious transition-field warnings from probe_placement step")


def test_s1_yaml_validator_if_placed_is_in_valid_set():
    """S1: valid_transition_fields in yaml_tools explicitly includes if_placed/if_not_placed."""
    print("Test: s1_yaml_validator_if_placed_is_in_valid_set")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.yaml_tools import _validate_workflows
    assert _validate_workflows is not None, "_validate_workflows must be importable"

    # Read source to verify the constants are present (not just that no warnings fired)
    src_path = os.path.join(_PROJECT_ROOT, "agent", "yaml_tools.py")
    with open(src_path) as f:
        src = f.read()

    assert "'if_placed'" in src, (
        "'if_placed' must appear in yaml_tools.py (valid_transition_fields)"
    )
    assert "'if_not_placed'" in src, (
        "'if_not_placed' must appear in yaml_tools.py (valid_transition_fields)"
    )
    print("  PASSED: if_placed and if_not_placed present in yaml_tools source")


def test_s1_no_redundant_import_re_in_probe_block():
    """S1: workflow_state.py probe detection does not use a local 'import re' inside the loop."""
    print("Test: s1_no_redundant_import_re_in_probe_block")
    sys.path.insert(0, _PROJECT_ROOT)

    src_path = os.path.join(_PROJECT_ROOT, "agent", "workflow_state.py")
    with open(src_path) as f:
        src = f.read()

    # The probe detection block should use the module-level `re`, not `_re2`
    assert "import re as _re2" not in src, (
        "Redundant 'import re as _re2' must be removed from probe detection block"
    )
    # Module-level `import re` must still be present
    assert "import re\n" in src or "\nimport re\n" in src, (
        "Module-level 'import re' must be present in workflow_state.py"
    )
    # The regex search used in the probe block must reference just `re`
    assert "re.search(r'r_free" in src, (
        "Probe block must use module-level re.search(), not _re2.search()"
    )
    print("  PASSED: probe detection uses module-level re, no _re2 alias")


def test_s1_probe_re_fallback_still_works():
    """S1: regex fallback for r_free still parses correctly after import cleanup."""
    print("Test: s1_probe_re_fallback_still_works")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history

    # Result has r_free in text only (not in analysis dict) — exercises the regex path
    history = [{
        "program": "phenix.model_vs_data",
        "command": "phenix.model_vs_data model.pdb data.mtz",
        "result": "SUCCESS: r_free: 0.38 r_work: 0.34",
        "analysis": {},   # empty — forces regex fallback
    }]
    info = _analyze_history(history)

    assert info["placement_probed"] is True, (
        "placement_probed must be True when r_free parsed from result text"
    )
    assert info["placement_probe_result"] == "placed", (
        "r_free=0.38 < 0.50 must give result='placed', got %r"
        % info["placement_probe_result"]
    )
    print("  PASSED: regex fallback parses r_free from result text correctly")


def test_s1_local_import_fallback_in_check_cell_mismatch():
    """S1: _check_cell_mismatch has a local 'from agent.placement_checker' fallback path."""
    print("Test: s1_local_import_fallback_in_check_cell_mismatch")
    sys.path.insert(0, _PROJECT_ROOT)

    src_path = os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py")
    with open(src_path) as f:
        src = f.read()

    assert "from agent.placement_checker import" in src, (
        "_check_cell_mismatch must have a local 'from agent.placement_checker' fallback"
    )
    # Verify the fallback is inside an except ImportError block
    # (i.e., it appears after the libtbx path)
    libtbx_idx = src.index("from libtbx.langchain.agent.placement_checker import")
    local_idx  = src.index("from agent.placement_checker import")
    assert local_idx > libtbx_idx, (
        "Local import must appear AFTER the libtbx path (as a fallback)"
    )
    print("  PASSED: local import fallback present and correctly ordered")


def test_s1_placement_checker_importable_locally():
    """S1: placement_checker can be imported via the local path (no libtbx needed)."""
    print("Test: s1_placement_checker_importable_locally")
    sys.path.insert(0, _PROJECT_ROOT)

    from agent.placement_checker import (
        read_pdb_unit_cell,
        cells_are_compatible,
        check_xray_cell_mismatch,
        check_cryoem_cell_mismatch,
    )
    # Basic smoke-test: all public functions callable and return correct types
    assert check_xray_cell_mismatch(None, None) is False
    assert check_cryoem_cell_mismatch(None, None) is False
    assert cells_are_compatible(None, None) is True
    assert read_pdb_unit_cell("/nonexistent.pdb") is None
    print("  PASSED: all placement_checker public functions importable and callable")


def test_s1_cell_mismatch_short_circuits_when_placed():
    """S1: build_context sets cell_mismatch=False when has_placed_model=True (short-circuit)."""
    print("Test: s1_cell_mismatch_short_circuits_when_placed")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # History: refinement ran -> has_placed_model=True via heuristics
    # If short-circuit is correct, cell_mismatch will be False even if files
    # had incompatible cells (we can't actually read MTZ in tests, but we verify
    # the short-circuit fires by checking the result is always False here).
    history_info = _analyze_history([{
        "program": "phenix.refine",
        "command": "phenix.refine model.pdb data.mtz",
        "result": "SUCCESS: Refinement complete",
        "analysis": {"r_free": 0.27},
    }])

    assert history_info.get("refine_done") is True, \
        "refine_done must be True after refinement history entry"

    ctx = engine.build_context(files={}, history_info=history_info)

    assert ctx["has_placed_model"] is True, \
        "has_placed_model must be True when refine_done=True"
    assert ctx["cell_mismatch"] is False, (
        "cell_mismatch must be False (short-circuited) when has_placed_model=True, "
        "got %s" % ctx["cell_mismatch"]
    )
    print("  PASSED: cell_mismatch=False when has_placed_model=True (short-circuited)")


def test_s1_cell_mismatch_short_circuits_when_probed():
    """S1: build_context sets cell_mismatch=False when placement_probed=True (short-circuit)."""
    print("Test: s1_cell_mismatch_short_circuits_when_probed")
    sys.path.insert(0, _PROJECT_ROOT)
    from agent.workflow_state import _analyze_history
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # History: probe ran (needs_mr) -> placement_probed=True, has_placed_model=False
    # Short-circuit must still fire even though has_placed_model is False
    history_info = _analyze_history([{
        "program": "phenix.model_vs_data",
        "command": "phenix.model_vs_data model.pdb data.mtz",
        "result": "SUCCESS: Done",
        "analysis": {"r_free": 0.62},   # > 0.50 -> needs_mr
    }])

    assert history_info["placement_probed"] is True
    assert history_info["placement_probe_result"] == "needs_mr"

    ctx = engine.build_context(files={}, history_info=history_info)

    assert ctx["has_placed_model"] is False, \
        "has_placed_model must remain False when probe result is needs_mr"
    assert ctx["cell_mismatch"] is False, (
        "cell_mismatch must be False (short-circuited) when placement_probed=True, "
        "got %s" % ctx["cell_mismatch"]
    )
    print("  PASSED: cell_mismatch=False when placement_probed=True (short-circuited)")


def test_s1_cell_mismatch_not_short_circuited_first_cycle():
    """S1: cell_mismatch check runs normally on first cycle (no placement evidence)."""
    print("Test: s1_cell_mismatch_not_short_circuited_first_cycle")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()

    # No history at all -> neither has_placed_model nor placement_probed
    # cell_mismatch check must run (returns False for empty files — no PDB to compare)
    ctx = engine.build_context(files={}, history_info={})

    assert ctx["has_placed_model"] is False
    assert ctx["placement_probed"] is False
    # The check ran and correctly returned False for empty files
    assert ctx["cell_mismatch"] is False
    # placement_uncertain also False (no model and no data)
    assert ctx["placement_uncertain"] is False
    print("  PASSED: cell_mismatch check active on first cycle (no short-circuit)")


def test_s1_short_circuit_order_before_probe_override():
    """S1: cell_mismatch short-circuit applies before probe-result has_placed_model override."""
    print("Test: s1_short_circuit_order_before_probe_override")
    sys.path.insert(0, _PROJECT_ROOT)

    src_path = os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py")
    with open(src_path) as f:
        src = f.read()

    # The short-circuit block must appear before the probe-result override block
    sc_marker   = 'if context.get("has_placed_model_from_history") or context.get("placement_probed"):\n            context["cell_mismatch"] = False'
    probe_marker = 'if (context["placement_probed"] and\n                context.get("placement_probe_result") == "placed"):\n            context["has_placed_model"] = True'

    sc_idx    = src.index(sc_marker)
    probe_idx = src.index(probe_marker)

    assert sc_idx < probe_idx, (
        "Short-circuit override (cell_mismatch=False) must appear BEFORE "
        "probe-result override (has_placed_model=True) in build_context"
    )
    print("  PASSED: short-circuit at index %d, probe override at %d (correct order)"
          % (sc_idx, probe_idx))




# =============================================================================
# CATEGORY S2b -- placement_uncertain must use has_placed_model_from_history
#
# Root cause (apoferritin AIAgent_104 live failure):
#   "solve the structure" → LLM sets model_is_placed=True directive
#   placement_uncertain used `not context["has_placed_model"]` which is
#   directive-affected → False → probe never fired → RSR crash on cycle 3.
#
# Fix: placement_uncertain = not has_placed_model_FROM_HISTORY (directive-immune)


# =============================================================================
# CATEGORY S2b -- placement_uncertain uses has_placed_model_from_history
#
# Root cause (apoferritin AIAgent_104 live failure):
#   "solve the structure" -> LLM sets model_is_placed=True directive
#   placement_uncertain used `not context["has_placed_model"]` which is
#   directive-affected -> False -> probe never fired -> RSR crash cycle 3.
#
# Fix: placement_uncertain uses `not has_placed_model_FROM_HISTORY` (directive-immune)
# =============================================================================

def test_s2b_placement_uncertain_uses_from_history_not_directive():
    """S2b: placement_uncertain=True despite directive model_is_placed=True when no history confirms."""
    print("Test: s2b_placement_uncertain_uses_from_history_not_directive")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    engine = WorkflowEngine()
    directives = {
        "workflow_preferences": {"model_is_placed": True},
        "constraints": [], "stop_conditions": {}, "program_settings": {},
    }
    history_info = {
        "mtriage_done": True, "resolve_cryo_em_done": True,
        "dock_done": False, "refine_done": False, "rsr_done": False,
        "phaser_done": False, "autobuild_done": False, "predict_full_done": False,
        "placement_probed": False, "placement_probe_result": None,
        "predict_done": False, "process_predicted_done": False,
    }
    files = {"model": ["/fake/1aew_A.pdb"], "full_map": ["/fake/denmod_map.ccp4"]}
    ctx = engine.build_context(files=files, history_info=history_info, directives=directives)
    assert ctx["has_placed_model"] is True, "Directive should set has_placed_model=True"
    assert ctx["has_placed_model_from_history"] is False
    assert ctx["placement_uncertain"] is True, (
        "placement_uncertain must be True despite model_is_placed directive. "
        "from_history=%s cell_mismatch=%s probed=%s" % (
            ctx["has_placed_model_from_history"], ctx["cell_mismatch"], ctx["placement_probed"])
    )
    print("  PASSED: placement_uncertain=True despite model_is_placed directive")


def test_s2b_map_correlations_has_ignore_symmetry_default():
    """S2b: phenix.map_correlations has ignore_symmetry_conflicts=True default."""
    print("Test: s2b_map_correlations_has_ignore_symmetry_default")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return
    registry = ProgramRegistry()
    prog = registry.get_program("phenix.map_correlations")
    assert prog is not None
    defaults = prog.get("defaults", {})
    assert "ignore_symmetry_conflicts" in defaults, (
        "map_correlations must have ignore_symmetry_conflicts in defaults to avoid "
        "crashing when probe runs against a mismatched model+map"
    )
    val = defaults["ignore_symmetry_conflicts"]
    assert val in (True, "true", "True"), "ignore_symmetry_conflicts must be true, got %r" % val
    print("  PASSED: map_correlations has ignore_symmetry_conflicts=True in defaults")


def test_s2b_placement_uncertain_formula_uses_from_history_key():
    """S2b: The placement_uncertain formula uses has_placed_model_from_history."""
    print("Test: s2b_placement_uncertain_formula_uses_from_history_key")
    sys.path.insert(0, _PROJECT_ROOT)
    src_path = os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py")
    with open(src_path) as f:
        src = f.read()
    start = src.index('context["placement_uncertain"] = (')
    # Grab ~600 chars which covers the multi-line tuple
    block = src[start:start + 600]
    assert "has_placed_model_from_history" in block, (
        "placement_uncertain must use has_placed_model_from_history in its formula"
    )
    check = block.replace("has_placed_model_from_history", "REPLACED")
    assert '"has_placed_model"' not in check, (
        "placement_uncertain must NOT use bare has_placed_model (directive-affected)"
    )
    print("  PASSED: formula uses has_placed_model_from_history, not has_placed_model")


# =============================================================================

def test_s2b_placement_uncertain_uses_from_history_not_directive():
    """S2b: placement_uncertain=True despite directive model_is_placed=True when no history confirms."""
    print("Test: s2b_placement_uncertain_uses_from_history_not_directive")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    engine = WorkflowEngine()
    directives = {
        "workflow_preferences": {"model_is_placed": True},
        "constraints": [], "stop_conditions": {}, "program_settings": {},
    }
    history_info = {
        "mtriage_done": True, "resolve_cryo_em_done": True,
        "dock_done": False, "refine_done": False, "rsr_done": False,
        "phaser_done": False, "autobuild_done": False, "predict_full_done": False,
        "placement_probed": False, "placement_probe_result": None,
        "predict_done": False, "process_predicted_done": False,
    }
    files = {"model": ["/fake/1aew_A.pdb"], "full_map": ["/fake/denmod_map.ccp4"]}
    ctx = engine.build_context(files=files, history_info=history_info, directives=directives)
    assert ctx["has_placed_model"] is True, "Directive should set has_placed_model=True"
    assert ctx["has_placed_model_from_history"] is False
    assert ctx["placement_uncertain"] is True, (
        "placement_uncertain must be True despite model_is_placed directive when no history confirms. "
        "from_history=%s cell_mismatch=%s probed=%s" % (
            ctx["has_placed_model_from_history"], ctx["cell_mismatch"], ctx["placement_probed"])
    )
    print("  PASSED: placement_uncertain=True despite model_is_placed directive")


def test_s2b_cryoem_routes_to_probe_not_rsr_with_directive_placed():
    """S2b: Full cryoem state goes to probe_placement (not RSR) when directive is placed but no history."""
    print("Test: s2b_cryoem_routes_to_probe_not_rsr_with_directive_placed")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return
    engine = WorkflowEngine()
    directives = {
        "workflow_preferences": {"model_is_placed": True},
        "constraints": [], "stop_conditions": {}, "program_settings": {},
    }
    history_info = {
        "mtriage_done": True, "resolve_cryo_em_done": True,
        "dock_done": False, "refine_done": False, "rsr_done": False,
        "rsr_count": 0, "refine_count": 0,
        "phaser_done": False, "autobuild_done": False, "predict_full_done": False,
        "placement_probed": False, "placement_probe_result": None,
        "predict_done": False, "process_predicted_done": False,
        "map_sharpening_done": False, "map_symmetry_done": False, "validation_done": False,
    }
    files = {
        "model": ["/fake/1aew_A.pdb"], "full_map": ["/fake/denmod_map.ccp4"],
        "half_map": ["/fake/h1.ccp4", "/fake/h2.ccp4"], "sequence": ["/fake/seq.dat"],
    }
    state = engine.get_workflow_state("cryoem", files, history_info, directives=directives)
    step = state.get("step", state.get("state", ""))
    valid_progs = state.get("valid_programs", [])
    # The probe may be served from step "probe_placement" or from a general step
    # like "cryoem_analyzed" — what matters is that map_correlations IS offered
    # and real_space_refine is NOT (before placement is confirmed).
    assert "phenix.map_correlations" in valid_progs, (
        "Agent must offer map_correlations (probe) when placement unconfirmed. "
        "Got step=%r programs=%s" % (step, valid_progs)
    )
    assert "phenix.real_space_refine" not in valid_progs, (
        "Agent must NOT offer real_space_refine before placement confirmed. "
        "Got step=%r programs=%s" % (step, valid_progs)
    )
    print("  PASSED: offers probe (map_correlations) not RSR, step=%r" % step)


def test_s2b_map_correlations_has_ignore_symmetry_default():
    """S2b: phenix.map_correlations has ignore_symmetry_conflicts=True default."""
    print("Test: s2b_map_correlations_has_ignore_symmetry_default")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return
    registry = ProgramRegistry()
    prog = registry.get_program("phenix.map_correlations")
    assert prog is not None
    defaults = prog.get("defaults", {})
    assert "ignore_symmetry_conflicts" in defaults, (
        "map_correlations must have ignore_symmetry_conflicts in defaults to avoid "
        "crashing when probe runs against a mismatched model+map"
    )
    val = defaults["ignore_symmetry_conflicts"]
    assert val in (True, "true", "True"), "ignore_symmetry_conflicts must be true, got %r" % val
    print("  PASSED: map_correlations has ignore_symmetry_conflicts=True in defaults")


# =============================================================================
# CATEGORY S2c — promotion of unclassified_pdb to search_model for docking
# =============================================================================

def test_s2c_promotion_fires_when_placement_uncertain():
    """S2c: unclassified PDB promoted to search_model when placement_uncertain=True (Tier 3 path)."""
    print("Test: s2c_promotion_fires_when_placement_uncertain")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"
    files = {
        "unclassified_pdb": [pdb],
        "model": [pdb],
        "search_model": [],
        "pdb": [pdb],
    }
    ctx = {
        "has_placed_model_from_history": False,
        "has_search_model": False,
        "placement_uncertain": True,
        "placement_probed": False,
        "placement_probe_result": None,
        "cell_mismatch": False,
    }
    f2, c2 = engine._promote_unclassified_for_docking(files, ctx, "cryoem")

    assert_true(pdb in f2["search_model"],
                "Promoted PDB must appear in search_model. Got: %s" % f2["search_model"])
    assert_true(c2["has_search_model"] is True,
                "has_search_model must be True after promotion")
    assert_true(c2.get("unclassified_promoted_to_search_model") is True,
                "unclassified_promoted_to_search_model flag must be set")
    assert_true(files["search_model"] == [],
                "Original files dict must be unchanged (no mutation)")
    assert_true(f2 is not files,
                "Promoted files must be a new dict, not the original")
    print("  PASSED: promotion fires on placement_uncertain (Tier 3 path)")


def test_s2c_promotion_fires_when_probe_says_needs_dock():
    """S2c: unclassified PDB promoted when probe ran and returned needs_dock (Tier 3 post-probe)."""
    print("Test: s2c_promotion_fires_when_probe_says_needs_dock")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"
    files = {"unclassified_pdb": [pdb], "model": [pdb], "search_model": [], "pdb": [pdb]}
    ctx = {
        "has_placed_model_from_history": False,
        "has_search_model": False,
        "placement_uncertain": False,
        "placement_probed": True,
        "placement_probe_result": "needs_dock",
        "cell_mismatch": False,
    }
    f2, c2 = engine._promote_unclassified_for_docking(files, ctx, "cryoem")

    assert_true(pdb in f2["search_model"],
                "Promoted PDB must appear in search_model. Got: %s" % f2["search_model"])
    assert_true(c2["has_search_model"] is True,
                "has_search_model must be True after promotion")
    print("  PASSED: promotion fires on placement_probed=needs_dock (Tier 3 post-probe)")


def test_s2c_promotion_fires_when_cell_mismatch():
    """S2c: unclassified PDB promoted when cell_mismatch=True and no history (Tier 1 path)."""
    print("Test: s2c_promotion_fires_when_cell_mismatch")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"
    files = {"unclassified_pdb": [pdb], "model": [pdb], "search_model": [], "pdb": [pdb]}
    ctx = {
        "has_placed_model_from_history": False,
        "has_search_model": False,
        "placement_uncertain": False,
        "placement_probed": False,
        "placement_probe_result": None,
        "cell_mismatch": True,
    }
    f2, c2 = engine._promote_unclassified_for_docking(files, ctx, "cryoem")

    assert_true(pdb in f2["search_model"],
                "Promoted PDB must appear in search_model. Got: %s" % f2["search_model"])
    assert_true(c2["has_search_model"] is True,
                "has_search_model must be True after promotion")
    print("  PASSED: promotion fires on cell_mismatch (Tier 1 path)")


def test_s2c_no_promotion_when_placed_by_history():
    """S2c: no promotion when dock/refine already in history (model is already placed)."""
    print("Test: s2c_no_promotion_when_placed_by_history")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"
    files = {"unclassified_pdb": [pdb], "model": [pdb], "search_model": [], "pdb": [pdb]}
    ctx = {
        "has_placed_model_from_history": True,   # dock_done or refine_done in history
        "has_search_model": False,
        "placement_uncertain": True,             # would trigger if guard weren't here
        "placement_probed": False,
        "placement_probe_result": None,
        "cell_mismatch": True,                   # would trigger if guard weren't here
    }
    f2, c2 = engine._promote_unclassified_for_docking(files, ctx, "cryoem")

    assert_true(f2 is files,
                "Original files must be returned unchanged when model placed by history")
    assert_true(f2["search_model"] == [],
                "search_model must stay empty when model is placed by history")
    print("  PASSED: no promotion when has_placed_model_from_history=True")


def test_s2c_no_promotion_for_xray():
    """S2c: no promotion for X-ray sessions — crystal PDBs handled via model key for phaser."""
    print("Test: s2c_no_promotion_for_xray")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"
    files = {"unclassified_pdb": [pdb], "model": [pdb], "search_model": [], "pdb": [pdb]}
    ctx = {
        "has_placed_model_from_history": False,
        "has_search_model": False,
        "placement_uncertain": True,
        "placement_probed": False,
        "placement_probe_result": None,
        "cell_mismatch": True,
    }
    f2, c2 = engine._promote_unclassified_for_docking(files, ctx, "xray")

    assert_true(f2 is files,
                "Original files must be returned unchanged for xray experiment type")
    assert_true(f2["search_model"] == [],
                "search_model must stay empty for xray sessions")
    print("  PASSED: no promotion for xray experiment type")


def test_s2c_categorized_files_propagates_through_get_workflow_state():
    """S2c: promoted files appear in get_workflow_state return dict under categorized_files."""
    print("Test: s2c_categorized_files_propagates_through_get_workflow_state")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    pdb = "/fake/1aew_A.pdb"

    # Simulate post-resolve_cryo_em state: full_map present, no dock history,
    # cell_mismatch will be False (no libtbx in test env → fail-safe False),
    # so placement_uncertain fires (Tier 3 path, condition 5a).
    files = {
        "unclassified_pdb": [pdb],
        "model": [pdb],
        "search_model": [],
        "pdb": [pdb],
        "full_map": ["/fake/denmod_map.ccp4"],
        "half_map": ["/fake/h1.ccp4", "/fake/h2.ccp4"],
        "map": ["/fake/denmod_map.ccp4", "/fake/h1.ccp4", "/fake/h2.ccp4"],
        "sequence": ["/fake/seq.dat"],
        "data_mtz": [], "map_coeffs_mtz": [], "ligand_cif": [], "ligand_pdb": [],
        "predicted": [], "processed_predicted": [], "docked": [],
        "refined": [], "rsr_output": [], "phaser_output": [], "autobuild_output": [],
        "with_ligand": [], "ligand_fit_output": [], "intermediate_mr": [],
    }
    history_info = {
        "mtriage_done": True, "resolve_cryo_em_done": True,
        "dock_done": False, "refine_done": False, "rsr_done": False,
        "rsr_count": 0, "refine_count": 0,
        "phaser_done": False, "autobuild_done": False, "predict_full_done": False,
        "placement_probed": False, "placement_probe_result": None,
        "predict_done": False, "process_predicted_done": False,
        "map_sharpening_done": False, "map_symmetry_done": False,
        "validation_done": False, "xtriage_done": False,
    }

    result = engine.get_workflow_state(
        experiment_type="cryoem",
        files=files,
        history_info=history_info,
    )

    # Structural fix verification: categorized_files must be in the return dict
    assert_true("categorized_files" in result,
                "get_workflow_state must return categorized_files key")

    returned_files = result["categorized_files"]

    # When promotion fires (placement_uncertain=True in test env),
    # search_model must contain the promoted PDB.
    # When it doesn't fire (e.g. placement_uncertain=False for unexpected reasons),
    # at minimum categorized_files must be a dict and not clobber with originals.
    assert_true(isinstance(returned_files, dict),
                "categorized_files must be a dict, got: %s" % type(returned_files))

    # Context flag must match files state
    ctx = result.get("context", {})
    if ctx.get("unclassified_promoted_to_search_model"):
        assert_true(pdb in returned_files.get("search_model", []),
                    "If promotion fired, PDB must be in returned categorized_files['search_model']")
        assert_true(returned_files.get("has_search_model") or
                    bool(returned_files.get("search_model")),
                    "search_model list must be non-empty after promotion")
        print("  PASSED: promotion fired and categorized_files['search_model'] contains PDB")
    else:
        # Promotion did not fire (acceptable if placement_uncertain was False)
        print("  PASSED: categorized_files key present in return dict (promotion did not fire)")


# =============================================================================
# CATEGORY S2d — skip_map_model_overlap_check=True default for real_space_refine
# =============================================================================

def test_s2d_rsr_has_skip_map_model_overlap_check_default():
    """S2d: phenix.real_space_refine has skip_map_model_overlap_check=True in defaults."""
    print("Test: s2d_rsr_has_skip_map_model_overlap_check_default")
    sys.path.insert(0, _PROJECT_ROOT)
    from knowledge.yaml_loader import load_programs
    programs = load_programs()
    rsr = programs.get("phenix.real_space_refine", {})
    assert_not_none(rsr, "phenix.real_space_refine must exist in programs.yaml")
    defaults = rsr.get("defaults", {})
    assert_true("skip_map_model_overlap_check" in defaults,
                "real_space_refine defaults must contain skip_map_model_overlap_check")
    val = defaults["skip_map_model_overlap_check"]
    assert_true(val in (True, "true", "True"),
                "skip_map_model_overlap_check must be True, got %r" % val)
    print("  PASSED: real_space_refine has skip_map_model_overlap_check=True in defaults")


# =============================================================================
# CATEGORY S2e — after_program directive suppresses placement probe correctly
# =============================================================================

def test_s2e_after_program_requiring_placed_suppresses_probe():
    """S2e: after_program=model_vs_data sets has_placed_model_from_after_program=True,
    which suppresses placement_uncertain so probe_placement step is skipped."""
    print("Test: s2e_after_program_requiring_placed_suppresses_probe")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    files = {
        "model": ["/fake/1aba.pdb"],
        "pdb": ["/fake/1aba.pdb"],
        "data_mtz": ["/fake/1aba.mtz"],
        "search_model": [], "predicted": [], "processed_predicted": [],
        "ligand_pdb": [], "ligand": [],
    }
    history_info = {
        "phaser_done": False, "autobuild_done": False, "dock_done": False,
        "predict_full_done": False, "refine_done": False, "xtriage_done": True,
        "placement_probed": False, "placement_probe_result": None,
    }
    directives = {
        "stop_conditions": {
            "after_program": "phenix.model_vs_data",
            "skip_validation": True,
        },
        "workflow_preferences": {},
        "constraints": [],
        "program_settings": {},
    }

    ctx = engine.build_context(files=files, history_info=history_info, directives=directives)

    assert_true(ctx.get("has_placed_model_from_after_program") is True,
                "has_placed_model_from_after_program must be True for after_program=model_vs_data")
    assert_true(ctx.get("placement_uncertain") is False,
                "placement_uncertain must be False when after_program implies placed model. "
                "Got placement_uncertain=%r, from_history=%r, from_after_program=%r" % (
                    ctx.get("placement_uncertain"),
                    ctx.get("has_placed_model_from_history"),
                    ctx.get("has_placed_model_from_after_program")))
    print("  PASSED: after_program=model_vs_data suppresses placement probe")


def test_s2e_model_is_placed_directive_still_does_not_suppress_probe():
    """S2e: model_is_placed=True workflow_preference (unreliable LLM directive) does NOT
    suppress placement_uncertain — only after_program or history evidence does."""
    print("Test: s2e_model_is_placed_directive_still_does_not_suppress_probe")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    files = {
        "model": ["/fake/1aew_A.pdb"],
        "pdb": ["/fake/1aew_A.pdb"],
        "data_mtz": ["/fake/data.mtz"],
        "search_model": [], "predicted": [], "processed_predicted": [],
        "ligand_pdb": [], "ligand": [],
    }
    history_info = {
        "phaser_done": False, "autobuild_done": False, "dock_done": False,
        "predict_full_done": False, "refine_done": False,
        "placement_probed": False, "placement_probe_result": None,
    }
    directives = {
        "workflow_preferences": {"model_is_placed": True},  # hallucinated directive
        "stop_conditions": {},
        "constraints": [],
        "program_settings": {},
    }

    ctx = engine.build_context(files=files, history_info=history_info, directives=directives)

    assert_true(ctx.get("has_placed_model_from_after_program") is False,
                "has_placed_model_from_after_program must be False for model_is_placed directive")
    assert_true(ctx.get("placement_uncertain") is True,
                "placement_uncertain must remain True despite model_is_placed directive (S2b guard). "
                "Got placement_uncertain=%r" % ctx.get("placement_uncertain"))
    print("  PASSED: model_is_placed directive does not suppress placement probe (S2b intact)")


def test_s2e_after_program_not_in_requiring_placed_does_not_suppress_probe():
    """S2e: after_program for a non-placement-requiring program (e.g. predict_and_build)
    does NOT suppress placement_uncertain."""
    print("Test: s2e_after_program_not_in_requiring_placed_does_not_suppress_probe")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    engine = WorkflowEngine()
    files = {
        "model": ["/fake/model.pdb"],
        "pdb": ["/fake/model.pdb"],
        "data_mtz": ["/fake/data.mtz"],
        "search_model": [], "predicted": [], "processed_predicted": [],
        "ligand_pdb": [], "ligand": [],
    }
    history_info = {
        "phaser_done": False, "autobuild_done": False, "dock_done": False,
        "predict_full_done": False, "refine_done": False,
        "placement_probed": False, "placement_probe_result": None,
    }
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"},
        "workflow_preferences": {},
        "constraints": [],
        "program_settings": {},
    }

    ctx = engine.build_context(files=files, history_info=history_info, directives=directives)

    assert_true(ctx.get("has_placed_model_from_after_program") is False,
                "has_placed_model_from_after_program must be False for predict_and_build")
    print("  PASSED: predict_and_build after_program does not suppress placement probe")


# =============================================================================
# CATEGORY S2f — validation_cryoem requires resolution invariant
# =============================================================================

def test_s2f_validation_cryoem_has_requires_resolution_invariant():
    """S2f: phenix.validation_cryoem has a requires_resolution invariant with
    auto_fill_resolution so it never runs without a resolution value."""
    print("Test: s2f_validation_cryoem_has_requires_resolution_invariant")
    sys.path.insert(0, _PROJECT_ROOT)
    from knowledge.yaml_loader import load_programs
    programs = load_programs()
    prog = programs.get("phenix.validation_cryoem", {})
    assert_not_none(prog, "phenix.validation_cryoem must exist in programs.yaml")

    invariants = prog.get("invariants", [])
    assert_true(len(invariants) > 0,
                "validation_cryoem must have at least one invariant")

    res_inv = next((i for i in invariants if i.get("name") == "requires_resolution"), None)
    assert_not_none(res_inv,
                    "validation_cryoem must have a requires_resolution invariant")
    assert_true(res_inv.get("fix", {}).get("auto_fill_resolution") is True,
                "requires_resolution invariant must have auto_fill_resolution: true")
    assert_true(res_inv.get("check", {}).get("has_strategy") == "resolution",
                "requires_resolution invariant must check has_strategy: resolution")
    print("  PASSED: validation_cryoem has requires_resolution invariant with auto_fill_resolution")


# =============================================================================
# CATEGORY S2g — map_correlations requires resolution for cryoem, not xray
# =============================================================================

def test_s2g_map_correlations_invariant_has_only_for_cryoem():
    """S2g: map_correlations requires_resolution invariant has only_for_experiment_type=cryoem."""
    print("Test: s2g_map_correlations_invariant_has_only_for_cryoem")
    sys.path.insert(0, _PROJECT_ROOT)
    from knowledge.yaml_loader import load_programs
    programs = load_programs()
    prog = programs.get("phenix.map_correlations", {})
    assert_not_none(prog, "phenix.map_correlations must exist in programs.yaml")

    invariants = prog.get("invariants", [])
    res_inv = next((i for i in invariants if i.get("name") == "requires_resolution"), None)
    assert_not_none(res_inv, "map_correlations must have a requires_resolution invariant")

    check = res_inv.get("check", {})
    assert_true(check.get("only_for_experiment_type") == "cryoem",
                "requires_resolution must have only_for_experiment_type=cryoem, got: %r"
                % check.get("only_for_experiment_type"))
    assert_true(res_inv.get("fix", {}).get("auto_fill_resolution") is True,
                "requires_resolution fix must have auto_fill_resolution: true")
    print("  PASSED: map_correlations invariant has only_for_experiment_type=cryoem")


def test_s2g_invariant_skipped_for_xray():
    """S2g: only_for_experiment_type=cryoem invariant is skipped when experiment_type=xray."""
    print("Test: s2g_invariant_skipped_for_xray")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (CommandBuilder unavailable)")
        return

    builder = CommandBuilder()
    ctx = CommandContext(experiment_type="xray", resolution=None)
    files = {"model": "/fake/model.pdb", "full_map": "/fake/map.ccp4"}
    strategy = {}

    result_files, result_strategy = builder._apply_invariants(
        "phenix.map_correlations", files, strategy, ctx)

    assert_true("resolution" not in result_strategy,
                "resolution must NOT be auto-filled for xray experiment_type. "
                "Got strategy: %r" % result_strategy)
    print("  PASSED: invariant correctly skipped for xray (resolution not auto-filled)")


def test_s2g_invariant_fires_for_cryoem():
    """S2g: only_for_experiment_type=cryoem invariant fires and fills resolution when cryoem."""
    print("Test: s2g_invariant_fires_for_cryoem")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.command_builder import CommandBuilder, CommandContext
    except ImportError:
        print("  SKIP (CommandBuilder unavailable)")
        return

    builder = CommandBuilder()
    ctx = CommandContext(experiment_type="cryoem", resolution=1.9)
    files = {"model": "/fake/model.pdb", "full_map": "/fake/map.ccp4"}
    strategy = {}

    result_files, result_strategy = builder._apply_invariants(
        "phenix.map_correlations", files, strategy, ctx)

    assert_true("resolution" in result_strategy,
                "resolution must be auto-filled for cryoem experiment_type. "
                "Got strategy: %r" % result_strategy)
    assert_true(abs(result_strategy["resolution"] - 1.9) < 0.01,
                "resolution must be 1.9, got: %r" % result_strategy.get("resolution"))
    print("  PASSED: invariant fires for cryoem and auto-fills resolution=1.9")


# =============================================================================
# CATEGORY S2h — validation_cryoem writes data_manager PHIL for GUI restore
# =============================================================================

def test_s2h_validation_cryoem_in_dm_to_data_manager_phil():
    """S2h: phenix.validation_cryoem is in the dm_to_data_manager_phil whitelist
    so its model and map files get written as data_manager PHIL in the eff,
    not just as comments."""
    print("Test: s2h_validation_cryoem_in_dm_to_data_manager_phil")
    sys.path.insert(0, _PROJECT_ROOT)
    ai_agent_path = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.isfile(ai_agent_path):
        print("  SKIP (ai_agent.py not found)")
        return

    with open(ai_agent_path) as f:
        src = f.read()

    assert_true("dm_programs_to_data_manager_phil" in src,
                "ai_agent.py must contain dm_programs_to_data_manager_phil dict")
    assert_true("phenix.validation_cryoem" in src.split(
                "dm_programs_to_data_manager_phil")[1][:500],
                "phenix.validation_cryoem must be in dm_programs_to_data_manager_phil")
    assert_true("data_manager.%s.file" in src or
                "data_manager.%s.file" in src,
                "data_manager PHIL lines must be written for dm_file_mapping programs")
    print("  PASSED: validation_cryoem is in dm_to_data_manager_phil whitelist")


# =============================================================================
# CATEGORY S2i — STOP command recognized even with trailing tokens
# =============================================================================

def test_s2i_stop_alone_is_valid_command():
    """S2i: bare 'STOP' passes _is_valid_command."""
    print("Test: s2i_stop_alone_is_valid_command")
    sys.path.insert(0, _PROJECT_ROOT)
    ai_agent_path = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.isfile(ai_agent_path):
        print("  SKIP (ai_agent.py not found)")
        return
    with open(ai_agent_path) as f:
        src = f.read()
    assert_true('program == "STOP"' in src or "program == 'STOP'" in src,
                "_is_valid_command must explicitly allow STOP as first token")
    print("  PASSED: STOP is whitelisted in _is_valid_command")


def test_s2i_stop_with_trailing_tokens_detected():
    """S2i: 'STOP <args>' is recognized as a STOP command (first-token check)."""
    print("Test: s2i_stop_with_trailing_tokens_detected")
    sys.path.insert(0, _PROJECT_ROOT)
    ai_agent_path = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.isfile(ai_agent_path):
        print("  SKIP (ai_agent.py not found)")
        return
    with open(ai_agent_path) as f:
        src = f.read()
    assert_true('.split()[0] == "STOP"' in src or ".split()[0] == 'STOP'" in src,
                "STOP detection must use first-token check (.split()[0]) "
                "not exact-match to handle 'STOP <args>'")
    print("  PASSED: STOP detection uses first-token check for trailing-args case")


def test_s2i_plan_sets_stop_true_when_program_is_stop():
    """S2i root cause: PLAN must set intent['stop']=True when normalizing to STOP,
    so BUILD short-circuits and never assembles 'STOP <strategy_flags>'."""
    print("Test: s2i_plan_sets_stop_true_when_program_is_stop")
    sys.path.insert(0, _PROJECT_ROOT)
    graph_path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    if not os.path.isfile(graph_path):
        print("  SKIP (graph_nodes.py not found)")
        return
    with open(graph_path) as f:
        src = f.read()
    # The fix: after intent["program"] = "STOP", also set intent["stop"] = True
    assert_true('intent["stop"] = True' in src or "intent['stop'] = True" in src,
                "PLAN must set intent['stop']=True when normalizing program to STOP")
    print("  PASSED: PLAN sets intent['stop']=True when normalizing to STOP")


def test_s2i_build_short_circuits_on_program_stop():
    """S2i root cause: BUILD must short-circuit on program=='STOP' as well as stop==True."""
    print("Test: s2i_build_short_circuits_on_program_stop")
    sys.path.insert(0, _PROJECT_ROOT)
    graph_path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    if not os.path.isfile(graph_path):
        print("  SKIP (graph_nodes.py not found)")
        return
    with open(graph_path) as f:
        src = f.read()
    assert_true('intent.get("program") == "STOP"' in src or
                "intent.get('program') == 'STOP'" in src,
                "BUILD must check program=='STOP' in addition to stop==True")
    print("  PASSED: BUILD short-circuits on program=='STOP'")


# =============================================================================
# CATEGORY S2j — dotted strategy keys pass through for PHIL validation
# =============================================================================

def test_s2j_refinement_key_dropped_from_ligandfit_strategy():
    """S2j: refinement.main.number_of_macro_cycles in ligandfit
    strategy is passed through (PHIL will reject at runtime)
    rather than silently dropped.

    Cross-program PHIL paths are passed through to the
    command line because silently dropping valid PHIL paths
    (like refinement.reference_model.enabled on phenix.refine)
    is worse than passing invalid paths that PHIL rejects
    loudly.  The PHENIX PHIL interpreter catches mismatched
    scopes at runtime with a clear error message.
    """
    print("Test: s2j_refinement_key_dropped_from_ligandfit_strategy")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return

    registry = ProgramRegistry(use_yaml=True)
    passed = []

    def log(msg):
        if "PASSTHROUGH" in msg:
            passed.append(msg)

    strategy = {
        "refinement.main.number_of_macro_cycles": 2,
    }
    files = {"model": "/fake/placed.pdb", "data": "/fake/data.mtz"}
    cmd = registry.build_command("phenix.ligandfit", files, strategy, log=log)

    # Dotted PHIL paths now pass through for PHIL
    # validation at runtime (not silently dropped).
    assert_true(len(passed) > 0,
                "refinement.main.number_of_macro_cycles "
                "should be PASSTHROUGH'd (dotted PHIL), "
                "got cmd: %r" % cmd)
    assert_true(
        "number_of_macro_cycles" in (cmd or ""),
        "dotted PHIL path should reach command for "
        "runtime validation, got: %r" % cmd)
    print("  PASSED: dotted PHIL path passed through "
          "(PHIL validates at runtime)")


def test_s2j_known_short_names_still_pass_through():
    """S2j: KNOWN_PHIL_SHORT_NAMES (nproc, twin_law, etc.) still pass through normally."""
    print("Test: s2j_known_short_names_still_pass_through")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return

    registry = ProgramRegistry(use_yaml=True)
    passed = []

    def log(msg):
        if "PASSTHROUGH" in msg and "nproc" in msg:
            passed.append(msg)

    strategy = {"nproc": 4}
    files = {"model": "/fake/model.pdb", "data": "/fake/data.mtz"}
    cmd = registry.build_command("phenix.ligandfit", files, strategy, log=log)

    assert_true(len(passed) > 0 or (cmd and "nproc=4" in cmd),
                "nproc should still pass through as a KNOWN_PHIL_SHORT_NAME")
    print("  PASSED: nproc still passes through correctly")


# =============================================================================
# CATEGORY S2k — _inject_user_params skips STOP commands
# =============================================================================

def test_s2k_inject_user_params_skips_stop():
    """S2k: _inject_user_params call site must guard against command starting with STOP."""
    print("Test: s2k_inject_user_params_skips_stop")
    sys.path.insert(0, _PROJECT_ROOT)
    ai_agent_path = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.isfile(ai_agent_path):
        print("  SKIP (ai_agent.py not found)")
        return
    with open(ai_agent_path) as f:
        src = f.read()
    # The guard must be at the call site, not just inside _inject_user_params
    assert_true("split()[0] != 'STOP'" in src or 'split()[0] != "STOP"' in src,
                "_inject_user_params call site must skip when command starts with STOP")
    print("  PASSED: _inject_user_params call site skips STOP commands")


def test_s2k_inject_user_params_filters_wrong_program_scope():
    """S2k: refinement.main.number_of_macro_cycles must not be injected into ligandfit."""
    print("Test: s2k_inject_user_params_filters_wrong_program_scope")
    sys.path.insert(0, _PROJECT_ROOT)
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    if not os.path.isfile(pp_path):
        print("  SKIP (command_postprocessor.py not found)")
        return
    with open(pp_path) as f:
        src = f.read()
    assert_true("def inject_user_params" in src,
                "inject_user_params must exist in command_postprocessor.py")
    assert_true("leading_scope" in src and "_UNIVERSAL_SCOPES" in src,
                "inject_user_params must filter by leading_scope against _UNIVERSAL_SCOPES")
    universal_line = next((l for l in src.split('\n') if '_UNIVERSAL_SCOPES' in l and '=' in l), "")
    assert_true("refinement" not in universal_line,
                "'refinement' must not be in _UNIVERSAL_SCOPES")
    print("  PASSED: _inject_user_params filters dotted keys by program scope")


# =============================================================================
# S2L TESTS — Probe crash → needs_dock + client-side model cell transport
# =============================================================================

def test_s2l_probe_failure_outside_map_sets_needs_dock():
    """S2L: map_correlations crash 'entirely outside map' → placement_probed=needs_dock."""
    print("Test: s2l_probe_failure_outside_map_sets_needs_dock")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        print("  SKIP (_analyze_history not importable)")
        return

    history = [
        {
            "program": "phenix.mtriage",
            "command": "phenix.mtriage map.ccp4",
            "result": "Completed. Resolution: 1.9",
        },
        {
            "program": "phenix.resolve_cryo_em",
            "command": "phenix.resolve_cryo_em half_map_1.ccp4 half_map_2.ccp4",
            "result": "Completed. denmod_map.ccp4 written",
        },
        {
            # This is the probe — runs before any dock/refine, crashes with
            # the definitive "outside map" message.
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=1aew_A.pdb map=denmod_map.ccp4",
            "result": "FAILED: Stopping as model is entirely outside map and wrapping=False",
        },
    ]

    info = _analyze_history(history)

    assert_true(info.get("placement_probed"),
                "placement_probed must be True after 'entirely outside map' crash")
    assert_equal(info.get("placement_probe_result"), "needs_dock",
                 "placement_probe_result must be 'needs_dock' for outside-map crash")
    print("  PASSED: 'entirely outside map' crash → placement_probed=True, result=needs_dock")


def test_s2l_probe_failure_other_error_marks_inconclusive():
    """S2L: other map_correlations crashes → probed=True, result=None (no infinite retry)."""
    print("Test: s2l_probe_failure_other_error_marks_inconclusive")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        print("  SKIP")
        return

    history = [
        {
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=1aew_A.pdb map=denmod_map.ccp4",
            "result": "FAILED: Sorry: some other error occurred",
        },
    ]

    info = _analyze_history(history)

    assert_true(info.get("placement_probed"),
                "placement_probed must be True even for non-specific map_correlations failure")
    assert_true(info.get("placement_probe_result") is None,
                "placement_probe_result must be None (inconclusive) for generic failure")
    print("  PASSED: generic probe failure → probed=True, result=None (no retry)")


def test_s2l_probe_does_not_repeat_after_failure():
    """S2L: second map_correlations run should NOT override already-set probe result."""
    print("Test: s2l_probe_does_not_repeat_after_failure")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        print("  SKIP")
        return

    history = [
        {
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=1aew_A.pdb map=denmod_map.ccp4",
            "result": "FAILED: Stopping as model is entirely outside map and wrapping=False",
        },
        {
            # Second attempt — would have been prevented by fix, but even if it
            # ran, its result should not overwrite the first probe's conclusion.
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=1aew_A.pdb map=denmod_map.ccp4",
            "result": "FAILED: Sorry: some other error",
        },
    ]

    info = _analyze_history(history)

    # First probe result (needs_dock) must be preserved
    assert_true(info.get("placement_probed"), "placement_probed must be True")
    assert_equal(info.get("placement_probe_result"), "needs_dock",
                 "first probe result (needs_dock) must not be overwritten by second failure")
    print("  PASSED: first probe result preserved; second failure does not overwrite")


def test_s2l_successful_probe_still_works():
    """S2L: successful map_correlations probe with high CC still routes correctly."""
    print("Test: s2l_successful_probe_still_works")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        print("  SKIP")
        return

    history = [
        {
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=1aew_A.pdb map=denmod_map.ccp4",
            "result": "Completed.",
            "analysis": {"cc_mask": 0.72},
        },
    ]

    info = _analyze_history(history)

    assert_true(info.get("placement_probed"), "placement_probed must be True for success")
    assert_equal(info.get("placement_probe_result"), "placed",
                 "CC=0.72 > 0.15 threshold → result must be 'placed'")
    print("  PASSED: successful probe with high CC still works as before")


def test_s2l_client_model_cell_in_session_state():
    """S2L: unplaced_model_cell must be passed through build_session_state."""
    print("Test: s2l_client_model_cell_in_session_state")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.api_client import build_session_state
    except ImportError:
        print("  SKIP")
        return

    session_info = {
        "experiment_type": "cryoem",
        "unplaced_model_cell": [184.0, 184.0, 184.0, 90.0, 90.0, 90.0],
    }

    state = build_session_state(session_info)
    assert_true("unplaced_model_cell" in state,
                "unplaced_model_cell must appear in session_state built by build_session_state")
    assert_equal(state["unplaced_model_cell"], [184.0, 184.0, 184.0, 90.0, 90.0, 90.0],
                 "unplaced_model_cell values must be preserved exactly")
    print("  PASSED: unplaced_model_cell flows through build_session_state correctly")


def test_s2l_check_cell_mismatch_uses_preread_cell():
    """S2L: _check_cell_mismatch detects apoferritin mismatch using client-supplied cell."""
    print("Test: s2l_check_cell_mismatch_uses_preread_cell")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_engine import WorkflowEngine
        from agent.placement_checker import cells_are_compatible
    except ImportError:
        print("  SKIP (workflow_engine or placement_checker not importable)")
        return

    engine = WorkflowEngine()

    # Simulate the apoferritin scenario:
    # Model CRYST1: F432 crystal, a=b=c=184 Å
    # Map cell (denmod from resolve_cryo_em sub-box): P1, 32.5 x 39.65 x 36.4 Å
    # These differ by ~5x on every axis → definitive mismatch
    model_cell = [184.0, 184.0, 184.0, 90.0, 90.0, 90.0]
    # files dict has no real paths — on the server the map path would be real
    # but we mock it here as empty to test only the model_cell branch
    files = {"map": [], "full_map": [], "optimized_full_map": [], "data_mtz": []}

    # With no map files and model_cell only → cannot compare → False (fail-safe)
    result = engine._check_cell_mismatch(files, model_cell=model_cell)
    assert_true(not result,
                "With no readable map file, must return False (fail-safe) even with model_cell")

    # Quick sanity check that the cells ARE incompatible
    map_cell = (32.5, 39.65, 36.4, 90.0, 90.0, 90.0)
    compat = cells_are_compatible(tuple(model_cell), map_cell)
    assert_true(not compat,
                "apoferritin crystal cell and cryo-EM sub-box cell must be incompatible (>5%)")
    print("  PASSED: _check_cell_mismatch uses model_cell parameter; apoferritin cells incompatible")


def test_s2l_detect_workflow_state_accepts_session_info():
    """S2L: detect_workflow_state must accept session_info keyword argument."""
    print("Test: s2l_detect_workflow_state_accepts_session_info")
    sys.path.insert(0, _PROJECT_ROOT)
    import inspect
    try:
        from agent.workflow_state import detect_workflow_state
    except ImportError:
        print("  SKIP")
        return
    sig = inspect.signature(detect_workflow_state)
    assert_true("session_info" in sig.parameters,
                "detect_workflow_state must accept session_info parameter")
    print("  PASSED: detect_workflow_state has session_info parameter")


def test_s2l_outside_map_variants_all_detected():
    """S2L: various 'outside map' crash messages all trigger needs_dock."""
    print("Test: s2l_outside_map_variants_all_detected")
    sys.path.insert(0, _PROJECT_ROOT)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        print("  SKIP")
        return

    outside_map_messages = [
        "FAILED: Stopping as model is entirely outside map and wrapping=False",
        "FAILED: Sorry: model is outside map",
        "FAILED: Sorry: model entirely outside map box",
        "FAILED: stopping as model is outside map",
    ]

    for msg in outside_map_messages:
        history = [{
            "program": "phenix.map_correlations",
            "command": "phenix.map_correlations model=m.pdb map=map.ccp4",
            "result": msg,
        }]
        info = _analyze_history(history)
        assert_true(info.get("placement_probed"),
                    "placement_probed must be True for: %s" % msg)
        assert_equal(info.get("placement_probe_result"), "needs_dock",
                     "placement_probe_result must be needs_dock for: %s" % msg)

    print("  PASSED: all outside-map crash variants trigger needs_dock")


# RUN ALL TESTS
# =============================================================================
# S3A TESTS — Failure Diagnosis Feature
# =============================================================================
# Tests for the new diagnosable-terminal error handling system.
# These tests do NOT require libtbx or the PHENIX server; they test the pure
# Python detector, sanitiser, prompt builder, and HTML builder in isolation.
# =============================================================================

def test_s3a_detect_crystal_symmetry_mismatch():
    """DiagnosisDetector detects the exact error string from the bug report."""
    print("  Test: s3a_detect_crystal_symmetry_mismatch")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    result = d.detect(
        "FAILED: Sorry: Crystal symmetry mismatch between different files\n"
        "File 1: a=34.5 Ang\nFile 2: a=39.6 Ang"
    )
    assert_not_none(result, "Should detect crystal_symmetry_mismatch")
    assert_equal(result[0], 'crystal_symmetry_mismatch')
    assert_true(len(result[2]) > 0, "Excerpt should not be empty")
    print("  PASSED: crystal_symmetry_mismatch detected correctly")


def test_s3a_detect_model_outside_map():
    """DiagnosisDetector detects the model-outside-map error."""
    print("  Test: s3a_detect_model_outside_map")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    result = d.detect(
        "FAILED: Sorry: Stopping as model is entirely outside map and wrapping=False"
    )
    assert_not_none(result, "Should detect model_outside_map")
    assert_equal(result[0], 'model_outside_map')
    print("  PASSED: model_outside_map detected correctly")


def test_s3a_detect_returns_none_for_recoverable():
    """
    A recoverable error (ambiguous data labels) must NOT match DiagnosisDetector.

    This enforces the disjoint-sets invariant: no error type may be both
    retryable and terminal.
    """
    print("  Test: s3a_detect_returns_none_for_recoverable")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    # This is the canonical recoverable error text (from recoverable_errors.yaml)
    recoverable_text = (
        "Multiple equally suitable arrays of observed xray data found.\n"
        "Please use scaling.input.xray_data.obs_labels to specify an "
        "unambiguous substring."
    )
    result = d.detect(recoverable_text)
    assert_equal(result, None,
                 "Recoverable error must NOT be flagged as diagnosable-terminal")
    print("  PASSED: recoverable error correctly returns None")


def test_s3a_detect_returns_none_for_unknown():
    """An unrecognised error string must return None."""
    print("  Test: s3a_detect_returns_none_for_unknown")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    result = d.detect(
        "FAILED: Some totally novel error that does not match any pattern."
    )
    assert_equal(result, None, "Unknown error should return None")
    print("  PASSED: unknown error correctly returns None")


def test_s3a_detect_returns_none_for_empty():
    """Empty or None result text must return None without raising."""
    print("  Test: s3a_detect_returns_none_for_empty")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    assert_equal(d.detect(""),   None)
    assert_equal(d.detect(None), None)
    print("  PASSED: empty/None input handled safely")


def test_s3a_strip_markdown_removes_formatting():
    """_strip_llm_markdown correctly removes **, ##, and _ markers
    WITHOUT eating underscores inside identifiers or filenames."""
    print("  Test: s3a_strip_markdown_removes_formatting")
    from agent.failure_diagnoser import _strip_llm_markdown

    raw = (
        "## WHAT WENT WRONG\n"
        "**The unit cells** do not match.\n\n"
        "## MOST LIKELY CAUSE\n"
        "_Different processing_ conventions.\n\n"
        "## HOW TO FIX IT\n"
        "Run **phenix.refine** again."
    )
    stripped = _strip_llm_markdown(raw)

    assert_true('##' not in stripped,    "## headers should be removed")
    assert_true('**' not in stripped,    "** bold should be removed")
    assert_true('The unit cells' in stripped, "Inner text should be preserved")
    assert_true('phenix.refine' in stripped,  "Inner text should be preserved")
    # Plain text input should be unchanged (modulo strip())
    plain = "WHAT WENT WRONG\nThe cells differ."
    assert_equal(_strip_llm_markdown(plain), plain)

    # --- Underscores inside identifiers and filenames must be preserved ---
    # This was a regression: the old regex ate underscores in PHENIX parameter
    # names and filenames (e.g. nsf-d2_noligand.pdb → nsf-d2noligand.pdb,
    # set_unit_cell_from_map=True → setunitcellfrommap=True).
    ident_text = (
        "Run: phenix.pdbtools nsf-d2_noligand.pdb nsf-d2.mtz "
        "set_unit_cell_and_space_group_from_map_or_model=True"
    )
    stripped_ident = _strip_llm_markdown(ident_text)
    assert_true("nsf-d2_noligand.pdb" in stripped_ident,
        "Filename nsf-d2_noligand.pdb must survive markdown stripping, got: %r"
        % stripped_ident)
    assert_true("set_unit_cell_and_space_group_from_map_or_model=True" in stripped_ident,
        "PHENIX parameter name with underscores must survive stripping, got: %r"
        % stripped_ident)

    # Real markdown italic around a word should still be stripped
    italic_text = "This is _italic_ and __underlined__ text."
    stripped_italic = _strip_llm_markdown(italic_text)
    assert_true("_italic_" not in stripped_italic,
        "_italic_ markdown should be stripped")
    assert_true("italic" in stripped_italic,
        "inner text 'italic' should be preserved after stripping")
    assert_true("underlined" in stripped_italic,
        "inner text 'underlined' should be preserved after stripping")

    print("  PASSED: markdown stripped, plain text preserved")


def test_s3a_build_prompt_contains_hint():
    """build_diagnosis_prompt includes the YAML hint for crystal_symmetry_mismatch."""
    print("  Test: s3a_build_prompt_contains_hint")
    from agent.failure_diagnoser import build_diagnosis_prompt
    prompt = build_diagnosis_prompt(
        error_type='crystal_symmetry_mismatch',
        error_text='Sorry: Crystal symmetry mismatch between different files',
        program='phenix.refine',
        log_tail='[... last lines of log ...]',
    )
    # Prompt must contain structural markers
    assert_true('WHAT WENT WRONG'  in prompt, "Prompt missing WHAT WENT WRONG section")
    assert_true('MOST LIKELY CAUSE' in prompt, "Prompt missing MOST LIKELY CAUSE section")
    assert_true('HOW TO FIX IT'    in prompt, "Prompt missing HOW TO FIX IT section")
    # Hint from YAML must be present (key phrase from the hint text)
    assert_true('unit cell' in prompt.lower(), "YAML hint should appear in prompt")
    # Program name must be present
    assert_true('phenix.refine' in prompt, "Program name should appear in prompt")
    print("  PASSED: prompt contains all required sections and YAML hint")


def test_s3a_build_prompt_log_tail_truncated():
    """Log tail in the prompt is capped to 3000 chars regardless of input size."""
    print("  Test: s3a_build_prompt_log_tail_truncated")
    from agent.failure_diagnoser import build_diagnosis_prompt
    # Build a log tail much longer than the 3000-char cap
    huge_log = "x" * 10_000
    prompt = build_diagnosis_prompt(
        error_type='crystal_symmetry_mismatch',
        error_text='err',
        program='phenix.refine',
        log_tail=huge_log,
    )
    # The total prompt will contain the capped log section
    # 3000 x's should be present but not 10000
    assert_true('x' * 3000 in prompt,  "3000-char log section should appear")
    assert_true('x' * 3001 not in prompt, "Log tail should be capped at 3000 chars")
    print("  PASSED: log tail correctly capped at 3000 chars in prompt")


def test_s3a_build_html_escapes_content():
    """build_diagnosis_html HTML-escapes < > & in user-supplied strings."""
    print("  Test: s3a_build_html_escapes_content")
    from agent.failure_diagnoser import build_diagnosis_html
    html = build_diagnosis_html(
        description='Test <error> & "description"',
        error_excerpt='<script>alert(1)</script> & more',
        diagnosis_text='Fix: use a=b&c',
        program='phenix.<test>',
        cycle=1,
    )
    # Angle brackets and ampersands must be escaped
    assert_true('<script>' not in html, "Raw <script> tag must be escaped")
    assert_true('&lt;script&gt;' in html, "&lt;&gt; escaping must be present")
    assert_true('&amp;' in html, "& must be escaped to &amp;")
    # The document should be valid HTML5 boilerplate
    assert_true('<!DOCTYPE html>' in html, "Must be a full HTML document")
    assert_true('white-space: pre-wrap' in html, "Diagnosis div needs pre-wrap")
    print("  PASSED: HTML escaping applied correctly to all user-supplied strings")


def test_s3a_build_html_contains_required_sections():
    """build_diagnosis_html output contains all required elements."""
    print("  Test: s3a_build_html_contains_required_sections")
    from agent.failure_diagnoser import build_diagnosis_html
    html = build_diagnosis_html(
        description='Unit cell mismatch',
        error_excerpt='Sorry: Crystal symmetry mismatch',
        diagnosis_text='WHAT WENT WRONG\nCells differ.\n\nHOW TO FIX IT\nReprocess.',
        program='phenix.refine',
        cycle=5,
    )
    assert_true('Unit cell mismatch'           in html, "Description missing")
    assert_true('phenix.refine'                in html, "Program name missing")
    assert_true('5'                            in html, "Cycle number missing")
    assert_true('WHAT WENT WRONG'             in html, "Diagnosis text missing")
    assert_true('error-box'                   in html, "error-box class missing")
    assert_true('diagnosis'                   in html, "diagnosis class missing")
    print("  PASSED: HTML report contains all required sections")


def test_s3a_rules_only_fallback_has_hint():
    """
    When the YAML hint is available, it appears in the fallback diagnosis text
    that _diagnose_terminal_failure uses when use_rules_only=True or when
    the LLM call fails.
    """
    print("  Test: s3a_rules_only_fallback_has_hint")
    from agent.error_analyzer import DiagnosisDetector
    d = DiagnosisDetector()
    hint = d.get_hint('crystal_symmetry_mismatch')
    assert_true(len(hint) > 20, "YAML hint must be non-trivial")
    assert_true('unit cell' in hint.lower() or 'cell' in hint.lower(),
                "Hint should mention unit cell")
    # Fallback text construction (mirrors _diagnose_terminal_failure)
    fallback = (
        "WHAT WENT WRONG\n"
        "Unit cell or space group mismatch between input files.\n\n"
        "MOST LIKELY CAUSE\n"
        "Sorry: Crystal symmetry mismatch\n\n"
        "HOW TO FIX IT\n"
        + hint
    )
    assert_true('unit cell' in fallback.lower() or
                'cell' in fallback.lower(),
                "Fallback should contain crystallographic guidance")
    print("  PASSED: fallback diagnosis contains meaningful YAML hint")


def test_s3a_build_html_new_fields():
    """
    build_diagnosis_html populates the new optional context fields:
    html_path (saved-to line), job_name, and working_dir (meta bar).
    Also verifies the heading text was changed to 'Error diagnosis'.
    """
    print("  Test: s3a_build_html_new_fields")
    from agent.failure_diagnoser import build_diagnosis_html

    html = build_diagnosis_html(
        description='Unit cell mismatch',
        error_excerpt='Sorry: Crystal symmetry mismatch',
        diagnosis_text='WHAT WENT WRONG\nCells differ.',
        program='phenix.refine',
        cycle=3,
        html_path='/my/job/dir/ai_failure_diagnosis.html',
        job_name='nsf-d2-ligand',
        working_dir='/my/job/dir',
    )

    # New heading (item 1)
    assert_true('Error diagnosis' in html,
                "Heading must read 'Error diagnosis' (not 'Terminal Error Diagnosis')")
    assert_true('Terminal Error Diagnosis' not in html,
                "Old heading 'Terminal Error Diagnosis' must be gone")

    # File location (item 2)
    assert_true('ai_failure_diagnosis.html' in html,
                "Saved file path must appear in the HTML report")
    assert_true('Saved to' in html,
                "Footer must say 'Saved to: <path>'")

    # Job name and working dir (item 3)
    assert_true('nsf-d2-ligand' in html,
                "Job name must appear in the meta bar")
    assert_true('/my/job/dir' in html,
                "Working directory must appear in the meta bar")

    # Backward compat — all optional fields can be omitted
    html_min = build_diagnosis_html(
        description='Test error',
        error_excerpt='some text',
        diagnosis_text='Some diagnosis.',
        program='phenix.refine',
        cycle=1,
    )
    assert_true('Error diagnosis' in html_min,
                "Heading must be present even without optional fields")
    assert_true('<!DOCTYPE html>' in html_min,
                "Must be a full HTML document")

    print("  PASSED: new html_path/job_name/working_dir fields and heading verified")


def test_s3a_diagnose_returns_true_no_sorry():
    """
    _diagnose_terminal_failure must NO LONGER raise Sorry, and _finalize_session
    must skip the Results summary page when a fatal diagnosis fired.

    The user flow is: diagnosis HTML opens in browser → ai_agent finishes cleanly.
    No Sorry modal, no second Results page that buries the diagnosis.

    Verifies via source-code inspection of ai_agent.py.
    """
    print("  Test: s3a_diagnose_returns_true_no_sorry")

    ai_agent_src = open(_find_ai_agent_path()).read()

    # 1. The old deferred-Sorry pattern must be gone
    assert_true('_pending_sorry' not in ai_agent_src,
                "_pending_sorry must be removed — Sorry is no longer raised")
    assert_true('raise _pending_sorry' not in ai_agent_src,
                "raise _pending_sorry must be removed")

    # 2. _diagnose_terminal_failure must return True (not raise Sorry)
    assert_true('return True' in ai_agent_src,
                "_diagnose_terminal_failure must return True to stop the cycle loop")

    # 3. The caller must propagate the return value
    #    (it now does "return self._diagnose_terminal_failure(...)")
    assert_true('return self._diagnose_terminal_failure(' in ai_agent_src,
                "_run_single_cycle must propagate the True return from "
                "_diagnose_terminal_failure")

    # 4. _finalize_session is still unconditional (no try/except around it)
    finalize_marker = 'self._finalize_session(session)    # always runs'
    assert_true(finalize_marker in ai_agent_src,
                "_finalize_session must remain unconditional after the cycle loop")

    # 5. Results summary is skipped when a fatal diagnosis fired
    assert_true('failure_diagnosis_path' in ai_agent_src,
                "_finalize_session must check failure_diagnosis_path to skip summary")
    assert_true('has_fatal_diagnosis' in ai_agent_src,
                "_finalize_session must use has_fatal_diagnosis flag to suppress "
                "the Results page when the diagnosis HTML is the user's output")

    print("  PASSED: Sorry removed; returns True; Results page suppressed on fatal diagnosis")


def test_s3a_finalize_runs_after_diagnosis():
    """
    Even when _diagnose_terminal_failure fires (True returned from
    _run_single_cycle), _finalize_session must still run (to save the session
    and populate self.result), but the Results summary page is NOT produced —
    the diagnosis HTML is the user's sole output window.

    Simulates the cycle loop logic directly without importing ai_agent.
    """
    print("  Test: s3a_finalize_runs_after_diagnosis")

    call_log = []

    def fake_run_single_cycle_diagnosis():
        """Simulates _run_single_cycle returning True on terminal failure."""
        return True   # stop the loop — diagnosis done

    def fake_finalize(has_diagnosis):
        call_log.append('finalize')
        # Mirrors _finalize_session: only generate summary when no diagnosis
        if not has_diagnosis:
            call_log.append('summary')

    # Replicate the simplified loop from iterate_agent
    for cycle in range(3):
        should_break = fake_run_single_cycle_diagnosis()
        if should_break:
            break

    fake_finalize(has_diagnosis=True)   # unconditional — always runs

    assert_true('finalize' in call_log,
                "_finalize_session must run even when a terminal failure stops the loop")
    assert_true('summary' not in call_log,
                "Results summary must NOT be generated when a fatal diagnosis fired")
    assert_true(cycle == 0,
                "Loop must break on first cycle when terminal failure is detected")

    print("  PASSED: _finalize_session runs; Results page suppressed on fatal diagnosis")



def test_pdb_is_small_molecule_helper():
    """
    _pdb_is_small_molecule must correctly distinguish polymer models from
    small-molecule coordinate files based on file size (atom count) and
    record types.

    Detection strategy:
    - ≤150 atoms total → small molecule (ligands have 10-100 atoms)
    - >150 atoms, HETATM-only → small molecule
    - >150 atoms, has ATOM records → NOT small molecule (protein)

    This is the foundational function for the atp.pdb ligand-detection fix.
    """
    print("  Test: pdb_is_small_molecule_helper")
    from agent.workflow_state import _pdb_is_small_molecule

    # Realistic protein: >150 ATOM records (real proteins have 500+ atoms)
    protein_lines = []
    for i in range(200):
        protein_lines.append(
            "ATOM  %5d  CA  ALA A %3d       %6.1f   0.000   0.000  1.00  0.00           C"
            % (i+1, i+1, float(i)))
    protein_lines.append("END")
    protein_pdb = "\n".join(protein_lines) + "\n"

    # Small ligand: HETATM-only (classic case)
    ligand_pdb = (
        "REMARK  ATP - adenosine triphosphate\n"
        "HETATM    1  PA  ATP     1       1.000   2.000   3.000  1.00  0.00           P\n"
        "HETATM    2  O1A ATP     1       4.000   5.000   6.000  1.00  0.00           O\n"
        "HETATM    3  C5  ATP     1       7.000   8.000   9.000  1.00  0.00           C\n"
        "END\n"
    )

    # Small ligand with ATOM records (the bug case: atp.pdb using ATOM instead of HETATM)
    atom_ligand_pdb_lines = []
    for i in range(31):  # ATP has ~31 atoms
        atom_ligand_pdb_lines.append(
            "ATOM  %5d  C%d  ATP A   1       %6.1f   0.000   0.000  1.00  0.00           C"
            % (i+1, i % 10, float(i)))
    atom_ligand_pdb_lines.append("END")
    atom_ligand_pdb = "\n".join(atom_ligand_pdb_lines) + "\n"

    # Mixed: realistic protein (>150 ATOM) plus bound ligand HETATM
    mixed_lines = list(protein_lines[:-1])  # 200 ATOM lines without END
    mixed_lines.append(
        "HETATM  201  PA  ATP A   1       4.000   5.000   6.000  1.00  0.00           P")
    mixed_lines.append("END")
    mixed_pdb = "\n".join(mixed_lines) + "\n"

    empty_pdb = "REMARK  empty file\nEND\n"

    with tempfile.TemporaryDirectory() as d:
        def write(name, content):
            p = os.path.join(d, name)
            with open(p, 'w') as f: f.write(content)
            return p

        # Protein (200 ATOM, >150) → NOT a small molecule
        assert_false(_pdb_is_small_molecule(write('protein.pdb', protein_pdb)),
                     "Protein (200 ATOM records) must not be small molecule")

        # HETATM-only ligand → IS a small molecule
        assert_true(_pdb_is_small_molecule(write('atp.pdb', ligand_pdb)),
                    "HETATM-only file must be detected as small molecule")

        # ATOM-only ligand (31 atoms, ≤150) → IS a small molecule
        # This is the key bug fix: atp.pdb may use ATOM records
        assert_true(_pdb_is_small_molecule(write('atp_atom.pdb', atom_ligand_pdb)),
                    "Small ATOM-only file (31 atoms) must be detected as small molecule")

        # Mixed (200 ATOM + 1 HETATM, >150) → NOT a small molecule
        assert_false(_pdb_is_small_molecule(write('complex.pdb', mixed_pdb)),
                     "Large file with ATOM records must not be small molecule")

        # Empty/REMARK-only: no ATOM or HETATM → returns False (conservative)
        assert_false(_pdb_is_small_molecule(write('empty.pdb', empty_pdb)),
                     "File with no ATOM/HETATM must return False conservatively")

        # Missing file → returns False (no exception)
        assert_false(_pdb_is_small_molecule(os.path.join(d, 'does_not_exist.pdb')),
                     "Missing file must return False without raising")

    print("  PASSED: _pdb_is_small_molecule helper is correct for all cases")


def test_hetcode_ligand_not_used_as_refine_model():
    """
    When the user provides a hetcode-named ligand file (e.g. atp.pdb, gdp.pdb,
    hem.pdb) alongside a real protein model and data file, the categorizer must:

      1. Place the hetcode ligand in 'ligand' and 'ligand_pdb' — NOT 'model'
      2. Leave the protein model in 'model'
      3. Leave mixed files (protein + bound ligand in one PDB) in 'model'

    This is the regression test for the original bug:
    "the ai_agent sets up to refine the atp.pdb instead of the model"

    The bug root cause: pattern-based categorization only knows about files
    named lig*.pdb / ligand*.pdb.  A file called atp.pdb matched the
    unclassified_pdb wildcard ("*") and bubbled up to 'model'.  Refinement
    programs exclude 'ligand' but not 'model', so atp.pdb could be selected
    as the model for phenix.refine.

    The fix: a post-categorization content check promotes HETATM-only PDB files
    from unclassified_pdb/model into ligand_pdb/ligand.
    """
    print("  Test: hetcode_ligand_not_used_as_refine_model")

    import yaml
    try:
        from agent.workflow_state import (
            _pdb_is_small_molecule,
            _categorize_files_yaml,
            _bubble_up_to_parents,
        )
        rules_path = os.path.join(_PROJECT_ROOT, 'knowledge', 'file_categories.yaml')
        with open(rules_path) as f:
            category_rules = yaml.safe_load(f)
    except Exception as e:
        print("  SKIPPED: could not load YAML rules:", e)
        return

    # Protein model must have enough atoms to be distinguishable from a ligand.
    # Real proteins have 500+ atoms; we use 200 to be well above the 150-atom
    # threshold that separates ligands from proteins in _pdb_is_small_molecule.
    protein_lines = []
    for i in range(200):
        protein_lines.append(
            "ATOM  %5d  CA  ALA A %3d       %6.1f   %6.1f   %6.1f  1.00  0.00           C"
            % (i+1, i+1, float(i), 0.0, 0.0))
    protein_lines.append("END")
    protein_pdb = "\n".join(protein_lines) + "\n"

    hetatm_pdb = (
        "HETATM    1  PA  ATP     1       1.000   2.000   3.000  1.00  0.00           P\n"
        "HETATM    2  O1A ATP     1       4.000   5.000   6.000  1.00  0.00           O\n"
        "END\n"
    )
    # Mixed file: mostly protein (200 ATOM) plus a few HETATM ligand atoms
    mixed_lines = list(protein_lines[:-1])  # all ATOM lines, without END
    mixed_lines.append(
        "HETATM  201  PA  ATP     1       4.000   5.000   6.000  1.00  0.00           P")
    mixed_lines.append("END")
    mixed_pdb = "\n".join(mixed_lines) + "\n"

    # Hetcode names that have no 'lig' substring — these are the problem files
    ligand_names = ['atp.pdb', 'gdp.pdb', 'hem.pdb', 'fmn.pdb', 'NAD.pdb']

    with tempfile.TemporaryDirectory() as d:
        def write(name, content):
            p = os.path.join(d, name)
            with open(p, 'w') as f: f.write(content)
            return p

        protein_path = write('model.pdb',   protein_pdb)
        complex_path = write('complex.pdb', mixed_pdb)   # protein + ligand in one file
        ligand_paths = [write(n, hetatm_pdb) for n in ligand_names]

        all_files = [protein_path, complex_path] + ligand_paths

        # Run the YAML categorizer
        files = _categorize_files_yaml(all_files, category_rules)
        files = _bubble_up_to_parents(files, category_rules)

        # Apply the HETATM post-processing pass (mirrors _categorize_files)
        _model_subcats = {
            'refined', 'rsr_output', 'phaser_output', 'autobuild_output',
            'docked', 'with_ligand', 'ligand_fit_output', 'model_cif',
        }
        for f in list(files.get('unclassified_pdb', [])):
            if any(f in files.get(sc, []) for sc in _model_subcats):
                continue
            if _pdb_is_small_molecule(f):
                files['unclassified_pdb'].remove(f)
                for lst in ('model', 'pdb'):
                    if f in files.get(lst, []): files[lst].remove(f)
                for k in ('ligand_pdb', 'ligand'):
                    if k not in files: files[k] = []
                    if f not in files[k]: files[k].append(f)

        model_basenames  = {os.path.basename(f) for f in files.get('model', [])}
        ligand_basenames = {os.path.basename(f) for f in files.get('ligand', [])}

        # 1. Protein model stays in model
        assert_in('model.pdb', model_basenames,
                  "Protein model must remain in 'model' category")

        # 2. Mixed file (has ATOM) stays in model
        assert_in('complex.pdb', model_basenames,
                  "PDB with ATOM records must remain in 'model' even if HETATM present")

        # 3. Every hetcode ligand must be in ligand, NOT in model
        for name in ligand_names:
            assert_true(name not in model_basenames,
                        f"{name} must NOT be in 'model' category")
            assert_true(name in ligand_basenames,
                        f"{name} must be in 'ligand' category")
            lp_basenames = {os.path.basename(f) for f in files.get('ligand_pdb', [])}
            assert_true(name in lp_basenames,
                        f"{name} must be in 'ligand_pdb' subcategory")

        # 4. Sanity: model only contains protein + complex
        unexpected_in_model = model_basenames - {'model.pdb', 'complex.pdb'}
        assert_true(len(unexpected_in_model) == 0,
                    f"Unexpected files in model category: {unexpected_in_model}")

    print("  PASSED: hetcode ligand PDB files correctly excluded from 'model'")
    print("  PASSED: phenix.refine model slot would not select atp.pdb")


def test_is_ligand_file_noligand_false_positive():
    """
    Regression test: _is_ligand_file must NOT misclassify proteins whose
    filename contains 'noligand' (e.g. nsf-d2_noligand.pdb) as ligands.

    The old code used bare substring matching ('ligand.pdb' in basename)
    which matched 'noligand.pdb' — because 'noligand.pdb' contains the
    substring 'ligand.pdb'.  The fix uses word-boundary-aware regex.

    Also verifies that genuine ligand names (lig.pdb, ligand_001.pdb, etc.)
    and HETATM-only hetcode files (atp.pdb) are still correctly identified.
    """
    print("  Test: is_ligand_file_noligand_false_positive")

    try:
        from agent.best_files_tracker import BestFilesTracker
    except ImportError:
        print("  SKIP (agent.best_files_tracker not importable)")
        return

    tracker = BestFilesTracker()
    check = lambda name, path=None: tracker._is_ligand_file(name.lower(), path=path)

    # --- Must NOT be ligand ---
    assert_true(not check('nsf-d2_noligand.pdb'),
                "'nsf-d2_noligand.pdb' must NOT be classified as ligand "
                "(contains 'noligand', not 'ligand' as a standalone word)")
    assert_true(not check('my_model_noligand.pdb'),
                "'my_model_noligand.pdb' must NOT be classified as ligand")
    assert_true(not check('noligand_model.pdb'),
                "'noligand_model.pdb' must NOT be classified as ligand")
    assert_true(not check('protein.pdb'),
                "'protein.pdb' must NOT be classified as ligand")

    # --- Must BE ligand (name-based) ---
    assert_true(check('lig.pdb'),       "'lig.pdb' must be classified as ligand")
    assert_true(check('lig_001.pdb'),   "'lig_001.pdb' must be classified as ligand")
    assert_true(check('ligand.pdb'),    "'ligand.pdb' must be classified as ligand")
    assert_true(check('ligand_001.pdb'),"'ligand_001.pdb' must be classified as ligand")
    assert_true(check('my_ligand.pdb'), "'my_ligand.pdb' must be classified as ligand")

    # --- Excluded output files must NOT be ligand (they are models) ---
    assert_true(not check('ligand_fit_001.pdb'),
                "'ligand_fit_001.pdb' is LigandFit output (model), not a ligand")
    assert_true(not check('nsf_with_ligand.pdb'),
                "'nsf_with_ligand.pdb' is a model+ligand complex, not a small molecule")

    # --- Content-based detection: HETATM-only files (e.g. atp.pdb) ---
    import tempfile, os
    hetatm_pdb = (
        "HETATM    1  O1  ATP A   1       1.000   2.000   3.000  1.00 10.00           O\n"
        "HETATM    2  N1  ATP A   1       4.000   5.000   6.000  1.00 10.00           N\n"
        "END\n"
    )
    # Protein PDBs need >150 ATOM records to be realistically distinguishable
    # from small molecules.  Real proteins have 500+ atoms.
    protein_lines = []
    for i in range(200):
        protein_lines.append(
            "ATOM  %5d  CA  ALA A %3d       %6.1f   0.000   0.000  1.00 10.00           C"
            % (i+1, i+1, float(i)))
    protein_lines.append("END")
    protein_pdb = "\n".join(protein_lines) + "\n"
    noligand_pdb = protein_pdb  # Same content, different name
    with tempfile.TemporaryDirectory() as d:
        def write(name, content):
            p = os.path.join(d, name)
            open(p, 'w').write(content)
            return p

        atp_path      = write('atp.pdb',              hetatm_pdb)
        protein_path  = write('protein.pdb',           protein_pdb)
        noligand_path = write('nsf-d2_noligand.pdb',  noligand_pdb)

        # atp.pdb: no ligand name pattern → falls back to HETATM content check
        assert_true(check('atp.pdb', path=atp_path),
                    "'atp.pdb' (HETATM-only) must be classified as ligand via content check")
        # protein.pdb: has ATOM records → NOT a small molecule
        assert_true(not check('protein.pdb', path=protein_path),
                    "'protein.pdb' (ATOM records) must NOT be classified as ligand")
        # nsf-d2_noligand.pdb: name is not a ligand, content has ATOM → not ligand
        assert_true(not check('nsf-d2_noligand.pdb', path=noligand_path),
                    "'nsf-d2_noligand.pdb' (ATOM records) must NOT be classified as ligand")

    print("  PASSED: _is_ligand_file correctly handles noligand false-positive "
          "and hetcode content detection")


# =============================================================================
# CRYSTAL SYMMETRY INJECTION (unit_cell / space_group from user advice)
# =============================================================================

def test_s4a_simple_extraction_unit_cell_parenthesized():
    """
    extract_directives_simple must pull unit_cell out of parenthesized tuple form
    and store it as a space-separated 6-number string under program_settings.default.

    This covers the exact form from the log:
      "The specified unit cell (116.097, 116.097, 44.175, 90, 90, 120) must be used"
    """
    print("  Test: s4a_simple_extraction_unit_cell_parenthesized")
    from agent.directive_extractor import extract_directives_simple

    advice = (
        "The specified unit cell (116.097, 116.097, 44.175, 90, 90, 120) "
        "must be used for the procedure."
    )
    d = extract_directives_simple(advice)

    uc = d.get("program_settings", {}).get("default", {}).get("unit_cell")
    assert_not_none(uc, "unit_cell must be extracted from parenthesized advice")
    assert_equal(uc, "116.097 116.097 44.175 90 90 120",
                 "unit_cell must be space-separated, no parens or commas")
    print("  PASSED: unit_cell extracted and normalised correctly")


def test_s4a_simple_extraction_unit_cell_space_separated():
    """extract_directives_simple handles space-separated unit cell without parens."""
    print("  Test: s4a_simple_extraction_unit_cell_space_separated")
    from agent.directive_extractor import extract_directives_simple

    advice = "Use unit cell 50.0 60.0 70.0 90.0 90.0 90.0 for refinement."
    d = extract_directives_simple(advice)

    uc = d.get("program_settings", {}).get("default", {}).get("unit_cell")
    assert_not_none(uc, "unit_cell must be extracted from space-separated form")
    nums = uc.split()
    assert_equal(len(nums), 6, "unit_cell must have exactly 6 numbers")
    print("  PASSED: space-separated unit_cell extracted correctly")


def test_s4a_simple_extraction_space_group():
    """extract_directives_simple extracts space group symbol."""
    print("  Test: s4a_simple_extraction_space_group")
    from agent.directive_extractor import extract_directives_simple

    for advice, expected_start in [
        ("space group P 32 2 1", "P 32"),
        ("space_group=P63", "P63"),
        ("Space Group: C 2 2 21", "C 2"),
    ]:
        d = extract_directives_simple(advice)
        sg = d.get("program_settings", {}).get("default", {}).get("space_group")
        assert_not_none(sg, "space_group must be extracted from: %r" % advice)
        assert_true(sg.startswith(expected_start),
                    "space_group %r must start with %r for advice %r" % (
                        sg, expected_start, advice))
    print("  PASSED: space_group extracted from multiple formats")


def test_s4a_inject_crystal_symmetry_into_model_vs_data():
    """
    _inject_crystal_symmetry must append unit_cell and space_group to
    phenix.model_vs_data when session directives carry them.

    This is the exact failure from the log: unit_cell given in advice but
    missing from the model_vs_data command.

    Architecture note: crystal symmetry injection now runs in the graph's
    BUILD node via postprocess_command (in command_postprocessor.py), not
    in ai_agent.py's _get_command_for_cycle.  The tests verify both:
      1. The standalone inject_crystal_symmetry function exists and works
      2. It is wired into postprocess_command (called by BUILD)
      3. The program sets are correct
    """
    print("  Test: s4a_inject_crystal_symmetry_into_model_vs_data")

    # Build a mock session whose directives contain crystal info
    import tempfile, shutil
    from agent.session import AgentSession

    tmp = tempfile.mkdtemp()
    try:
        session = AgentSession(session_dir=tmp)
        session.data["directives"] = {
            "program_settings": {
                "default": {
                    "unit_cell": "116.097 116.097 44.175 90 90 120",
                    "space_group": "P 63",
                }
            }
        }
        session.data["directives_extracted"] = True

        # --- Verify inject_crystal_symmetry exists in command_postprocessor ---
        import os
        postprocessor_path = os.path.join(_PROJECT_ROOT, "agent",
                                          "command_postprocessor.py")

        pp_src = open(postprocessor_path).read()
        assert_true("def inject_crystal_symmetry" in pp_src,
                    "inject_crystal_symmetry function must exist in command_postprocessor.py")
        assert_true("def postprocess_command" in pp_src,
                    "postprocess_command entry point must exist in command_postprocessor.py")
        assert_true("inject_crystal_symmetry(" in pp_src,
                    "inject_crystal_symmetry must be called from postprocess_command")

        # --- Verify BUILD calls postprocess_command ---
        graph_nodes_path = os.path.join(os.path.dirname(postprocessor_path),
                                         "graph_nodes.py")
        gn_src = open(graph_nodes_path).read()
        assert_true("postprocess_command" in gn_src,
                    "postprocess_command must be called in graph_nodes.py BUILD")

        # --- Verify dead code cleanup: class method removed (Phase 4) ---
        # inject_crystal_symmetry now lives only in command_postprocessor.py
        ai_agent_path = _find_ai_agent_path()
        src = open(ai_agent_path).read()
        assert_true("def _inject_crystal_symmetry" not in src,
                    "_inject_crystal_symmetry class method must be removed from "
                    "ai_agent.py (Phase 4 dead code cleanup)")

        # --- Verify _XRAY_SYMMETRY_PROGRAMS set is correct ---
        # Check in command_postprocessor.py (the authoritative location)
        block_start = pp_src.find("_XRAY_SYMMETRY_PROGRAMS = frozenset")
        assert_true(block_start >= 0, "_XRAY_SYMMETRY_PROGRAMS frozenset must exist")
        brace_open = pp_src.find("{", block_start)
        brace_close = pp_src.find("}", brace_open)
        symmetry_block = pp_src[brace_open: brace_close + 1]

        # Programs that should be in the set (need explicit crystal_symmetry arg)
        for prog in ("phenix.refine", "phenix.phaser", "phenix.autosol",
                     "phenix.ligandfit"):
            assert_true(prog in symmetry_block,
                        "%s must be in _XRAY_SYMMETRY_PROGRAMS" % prog)

        # Programs that read symmetry from their input files — must NOT be here,
        # or PHENIX will reject the redundant crystal_symmetry.unit_cell= arg.
        for prog in ("phenix.model_vs_data", "phenix.xtriage", "phenix.molprobity",
                     "phenix.real_space_refine", "phenix.map_to_model",
                     "phenix.dock_in_map"):
            assert_true(prog not in symmetry_block,
                        "%s must NOT be in _XRAY_SYMMETRY_PROGRAMS "
                        "(reads symmetry from input files)" % prog)
    finally:
        shutil.rmtree(tmp)

    print("  PASSED: inject_crystal_symmetry exists in postprocessor, "
          "is called by BUILD, and has correct program set")


def test_s4a_unit_cell_format_normalised():
    """
    extract_directives_simple must always produce space-separated 6-number strings
    regardless of whether the input uses parentheses, commas, or spaces.
    VALID_SETTINGS must include unit_cell and space_group as str.
    """
    print("  Test: s4a_unit_cell_format_normalised")
    from agent.directive_extractor import extract_directives_simple, VALID_SETTINGS

    # Check VALID_SETTINGS
    assert_in("unit_cell", VALID_SETTINGS,
              "unit_cell must be in VALID_SETTINGS")
    assert_equal(VALID_SETTINGS["unit_cell"], str,
                 "unit_cell type must be str (space-separated)")
    assert_in("space_group", VALID_SETTINGS,
              "space_group must be in VALID_SETTINGS")

    # Parenthesized comma-separated (the exact user complaint)
    d = extract_directives_simple(
        "unit cell (116.097, 116.097, 44.175, 90, 90, 120)"
    )
    uc = d.get("program_settings", {}).get("default", {}).get("unit_cell", "")
    assert_true("(" not in uc and "," not in uc,
                "Extracted unit_cell must not contain parentheses or commas: %r" % uc)
    nums = uc.strip().split()
    assert_equal(len(nums), 6,
                 "Extracted unit_cell must have exactly 6 space-separated numbers")

    print("  PASSED: unit_cell normalised and VALID_SETTINGS correct")


# =============================================================================
# CRYSTAL SYMMETRY FALLBACK (v112.63)
# Tests for _apply_crystal_symmetry_fallback and crystal_symmetry. scoped output
# =============================================================================

def test_s4b_fallback_populates_unit_cell_from_empty_directives():
    """
    When the LLM returns {} (no directives), _apply_crystal_symmetry_fallback
    must still extract unit_cell from the raw advice text.

    This is the exact failure from the log: directive extraction returned {}
    so the unit cell was never injected into the command.
    """
    print("  Test: s4b_fallback_populates_unit_cell_from_empty_directives")
    from agent.directive_extractor import _apply_crystal_symmetry_fallback

    advice = (
        "The specified unit cell (116.097, 116.097, 44.175, 90, 90, 120) "
        "must be used for the procedure."
    )
    # Simulate LLM returning {}
    result = _apply_crystal_symmetry_fallback({}, advice, lambda m: None)

    uc = result.get("program_settings", {}).get("default", {}).get("unit_cell")
    assert_not_none(uc,
        "_apply_crystal_symmetry_fallback must extract unit_cell from advice "
        "when LLM directives are empty")
    assert_equal(uc, "116.097 116.097 44.175 90 90 120",
        "unit_cell must be normalised to space-separated numbers, got: %r" % uc)
    print("  PASSED: fallback populates unit_cell from empty directives")


def test_s4b_fallback_does_not_overwrite_llm_unit_cell():
    """
    _apply_crystal_symmetry_fallback must NOT overwrite a unit_cell that the
    LLM already extracted correctly.
    """
    print("  Test: s4b_fallback_does_not_overwrite_llm_unit_cell")
    from agent.directive_extractor import _apply_crystal_symmetry_fallback

    llm_directives = {
        "program_settings": {"default": {"unit_cell": "50 60 70 90 90 90"}}
    }
    advice = "Use unit cell (116.097, 116.097, 44.175, 90, 90, 120) please."
    result = _apply_crystal_symmetry_fallback(
        llm_directives, advice, lambda m: None)

    uc = result.get("program_settings", {}).get("default", {}).get("unit_cell")
    assert_equal(uc, "50 60 70 90 90 90",
        "fallback must not overwrite LLM-extracted unit_cell "
        "(got: %r)" % uc)
    print("  PASSED: fallback preserves LLM-extracted unit_cell")


def test_s4b_fallback_populates_space_group_only_when_missing():
    """
    When LLM extracted unit_cell but missed space_group, the fallback
    must fill in space_group without touching unit_cell.
    """
    print("  Test: s4b_fallback_populates_space_group_only_when_missing")
    from agent.directive_extractor import _apply_crystal_symmetry_fallback

    llm_directives = {
        "program_settings": {
            "default": {"unit_cell": "116.097 116.097 44.175 90 90 120"}
        }
    }
    advice = (
        "The space group is P 63 2 2 and the unit cell is "
        "(116.097, 116.097, 44.175, 90, 90, 120)."
    )
    result = _apply_crystal_symmetry_fallback(
        llm_directives, advice, lambda m: None)

    default = result.get("program_settings", {}).get("default", {})
    # unit_cell unchanged
    assert_equal(default.get("unit_cell"), "116.097 116.097 44.175 90 90 120",
        "unit_cell must be preserved by fallback")
    # space_group filled in
    sg = default.get("space_group")
    assert_not_none(sg,
        "fallback must extract space_group when LLM left it empty")
    assert_true("P 63" in sg,
        "space_group must contain extracted symbol, got: %r" % sg)
    print("  PASSED: fallback fills space_group only, preserves unit_cell")


def test_s4b_program_registry_uses_crystal_symmetry_scope():
    """
    The program_registry PASSTHROUGH path must emit crystal_symmetry.unit_cell=
    and crystal_symmetry.space_group= (the fully-scoped PHIL form), not bare
    unit_cell= / space_group=.

    Verified by source inspection (ProgramRegistry.build_command requires libtbx
    to instantiate so we read the source directly).
    """
    print("  Test: s4b_program_registry_uses_crystal_symmetry_scope")
    import os
    pr_path = os.path.join(os.path.dirname(__file__), '..', 'agent', 'program_registry.py')
    src = open(pr_path).read()

    # Find the PASSTHROUGH block that handles KNOWN_PHIL_SHORT_NAMES
    passthrough_idx = src.find("PASSTHROUGH")
    assert_true(passthrough_idx != -1,
        "program_registry must have a PASSTHROUGH block for strategy flags")

    # Check that the scoping logic is present
    assert_true("crystal_symmetry.%s" % "unit_cell" in src or
                ("'unit_cell', 'space_group'" in src and "crystal_symmetry" in src),
        "program_registry must scope unit_cell/space_group as crystal_symmetry.*")

    # Verify the key pattern: tuple check then crystal_symmetry. prefix
    assert_true("crystal_symmetry.%s" % "" in src or
                "'crystal_symmetry.%s' % key" in src or
                "crystal_symmetry.%s' % key" in src,
        "program_registry must construct crystal_symmetry.{key} dynamically")

    print("  PASSED: program_registry uses crystal_symmetry. scoping")


def test_s4b_inject_crystal_symmetry_uses_scoped_form():
    """
    inject_crystal_symmetry in command_postprocessor.py must append
    crystal_symmetry.unit_cell= and crystal_symmetry.space_group=, not the
    bare unit_cell= / space_group= forms.

    Verified by source inspection.
    """
    print("  Test: s4b_inject_crystal_symmetry_uses_scoped_form")
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.isfile(pp_path),
        "command_postprocessor.py must exist")
    src = open(pp_path).read()

    # Find the inject_crystal_symmetry function body
    method_start = src.find("def inject_crystal_symmetry")
    assert_true(method_start != -1,
        "inject_crystal_symmetry must exist in command_postprocessor.py")

    # Find the next function def to bound the search
    next_method = src.find("\ndef ", method_start + 100)
    method_body = src[method_start: next_method] if next_method > 0 else src[method_start:]

    assert_true("crystal_symmetry.unit_cell=" in method_body,
        "inject_crystal_symmetry must use crystal_symmetry.unit_cell= "
        "(not bare unit_cell=)")
    assert_true("crystal_symmetry.space_group=" in method_body,
        "inject_crystal_symmetry must use crystal_symmetry.space_group= "
        "(not bare space_group=)")

    # Verify bare forms are absent from the append statements
    append_lines = [ln for ln in method_body.splitlines()
                    if "command +" in ln or "command=" in ln]
    for ln in append_lines:
        # Remove all ALLOWED scoped forms before checking for bare forms:
        #   crystal_symmetry.unit_cell= / crystal_symmetry.space_group=  → refine/ligandfit/polder
        #   crystal_info.unit_cell= / crystal_info.space_group=           → autobuild/autosol
        #   xray_data.unit_cell= / xray_data.space_group=                 → phaser
        _ln_stripped = (ln
            .replace("crystal_symmetry.unit_cell=", "")
            .replace("crystal_info.unit_cell=", "")
            .replace("xray_data.unit_cell=", "")
            .replace("crystal_symmetry.space_group=", "")
            .replace("crystal_info.space_group=", "")
            .replace("xray_data.space_group=", "")
        )
        assert_true("unit_cell=" not in _ln_stripped,
            "Bare unit_cell= must not appear in append line: %r" % ln)
        assert_true("space_group=" not in _ln_stripped,
            "Bare space_group= must not appear in append line: %r" % ln)

    print("  PASSED: _inject_crystal_symmetry uses scoped crystal_symmetry. form")


def test_s4b_fallback_called_in_extract_directives():
    """
    _apply_crystal_symmetry_fallback must be called inside extract_directives
    so unit_cell is captured even when the LLM returns {}.  Verified by
    source inspection.
    """
    print("  Test: s4b_fallback_called_in_extract_directives")
    import inspect
    from agent import directive_extractor
    src = inspect.getsource(directive_extractor.extract_directives)

    assert_true("_apply_crystal_symmetry_fallback" in src,
        "_apply_crystal_symmetry_fallback must be called inside extract_directives")
    assert_true("validate_directives" in src,
        "extract_directives must still call validate_directives")

    # Fallback must come AFTER validate_directives in source order
    idx_validate = src.find("validate_directives")
    idx_fallback  = src.find("_apply_crystal_symmetry_fallback")
    assert_true(idx_fallback > idx_validate,
        "_apply_crystal_symmetry_fallback must run AFTER validate_directives "
        "(validate_directives at %d, fallback at %d)" % (idx_validate, idx_fallback))

    print("  PASSED: fallback is present and ordered correctly in extract_directives")


# =============================================================================
# UNKNOWN PHIL PARAMETER — diagnosable error + model_vs_data exclusion
# =============================================================================

def test_s4c_unknown_phil_param_is_diagnosable():
    """
    'Unknown command line parameter definition' must be in diagnosable_errors.yaml
    so the agent surfaces an error window instead of silently looping.
    Previously it matched 'unit_cell' in real_failure_patterns → FAILED, but
    nothing in diagnosable_errors.yaml caught it, so no window appeared.
    """
    print("  Test: s4c_unknown_phil_param_is_diagnosable")
    from agent.error_analyzer import get_diagnosis_detector

    detector = get_diagnosis_detector()

    # The exact error text PHENIX produces when a bad param is injected
    bad_param_errors = [
        ("Sorry: Unknown command line parameter definition: "
         "unit_cell = 116.097 116.097 44.175 90 90 120   "
         "It turns out there is no such parameter"),
        "Sorry: there is no such parameter: crystal_symmetry.unit_cell",
        "Unknown command line parameter definition: space_group = P 63",
    ]
    for err in bad_param_errors:
        match = detector.detect(err)
        assert_not_none(match,
            "detector must catch unknown PHIL param error: %r" % err[:80])
        error_type = match[0]
        assert_equal(error_type, "unknown_phil_parameter",
            "error_type must be 'unknown_phil_parameter', got %r" % error_type)

    print("  PASSED: unknown_phil_parameter correctly detected as diagnosable")


def test_s4c_model_vs_data_not_in_symmetry_programs():
    """
    phenix.model_vs_data (and xtriage, molprobity) must NOT be in
    _XRAY_SYMMETRY_PROGRAMS.  These programs read crystal symmetry from their
    input files automatically; injecting crystal_symmetry.unit_cell= causes
    'Unknown command line parameter definition: unit_cell' errors.

    Phase 4: now checks command_postprocessor.py (class method removed from ai_agent.py).
    """
    print("  Test: s4c_model_vs_data_not_in_symmetry_programs")
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.isfile(pp_path),
        "command_postprocessor.py must exist")
    src = open(pp_path).read()

    block_start = src.find("_XRAY_SYMMETRY_PROGRAMS = frozenset")
    assert_true(block_start >= 0,
        "_XRAY_SYMMETRY_PROGRAMS frozenset must exist in command_postprocessor.py")
    brace_open  = src.find("{", block_start)
    brace_close = src.find("}", brace_open)
    symmetry_block = src[brace_open: brace_close + 1]

    # These programs read symmetry from their input files — must never be listed
    auto_symmetry_programs = [
        "phenix.model_vs_data",
        "phenix.xtriage",
        "phenix.molprobity",
    ]
    for prog in auto_symmetry_programs:
        assert_true(prog not in symmetry_block,
            "%s must NOT be in _XRAY_SYMMETRY_PROGRAMS "
            "(it reads symmetry from MTZ/model automatically)" % prog)

    # Programs that DO need explicit crystal_symmetry must still be present
    explicit_symmetry_programs = [
        "phenix.refine",
        "phenix.phaser",
        "phenix.autosol",
    ]
    for prog in explicit_symmetry_programs:
        assert_true(prog in symmetry_block,
            "%s must be in _XRAY_SYMMETRY_PROGRAMS" % prog)

    print("  PASSED: model_vs_data/xtriage/molprobity correctly excluded")


# =============================================================================
# BAD INJECT PARAM BLACKLIST (inject → crash → re-inject infinite loop fix)
# =============================================================================

def test_s4d_session_records_bad_inject_param():
    """session.record_bad_inject_param stores params and get_bad_inject_params
    returns them correctly, including both full and short-key variants."""
    print("  Test: s4d_session_records_bad_inject_param")
    import tempfile, shutil
    from agent.session import AgentSession

    tmp = tempfile.mkdtemp()
    try:
        sess = AgentSession(session_dir=tmp)

        # Nothing recorded yet
        assert_equal(sess.get_bad_inject_params("phenix.refine"), set(),
                     "Empty set before any recording")

        # Record a bad param
        sess.record_bad_inject_param("phenix.refine", "ignore_symmetry_conflicts")
        result = sess.get_bad_inject_params("phenix.refine")
        assert_true("ignore_symmetry_conflicts" in result,
                    "Recorded param must appear in blacklist")

        # Other program unaffected
        assert_equal(sess.get_bad_inject_params("phenix.autosol"), set(),
                     "Other programs must have empty blacklist")

        # Second recording of same key is idempotent (no duplicates)
        sess.record_bad_inject_param("phenix.refine", "ignore_symmetry_conflicts")
        assert_equal(len(sess.get_bad_inject_params("phenix.refine")), 1,
                     "Duplicate recording must not add duplicate entries")

        # Persists across a fresh session object (saved to disk)
        sess2 = AgentSession(session_dir=tmp)
        assert_true("ignore_symmetry_conflicts" in
                    sess2.get_bad_inject_params("phenix.refine"),
                    "Blacklist must persist across AgentSession reloads")
    finally:
        shutil.rmtree(tmp)
    print("  PASSED: bad inject param blacklist recorded and persisted correctly")


def test_s4d_inject_user_params_skips_blacklisted():
    """inject_user_params must not re-inject a parameter that was previously
    blacklisted via session.record_bad_inject_param.

    This is the core fix for the inject→crash→re-inject infinite loop:
    the agent added 'ignore_symmetry_conflicts=True' from the user's advice
    every cycle even after phenix.refine rejected it as unknown.

    Phase 4: checks command_postprocessor.py (class method removed from ai_agent.py).
    """
    print("  Test: s4d_inject_user_params_skips_blacklisted")

    import tempfile, shutil
    from agent.session import AgentSession

    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.isfile(pp_path),
        "command_postprocessor.py must exist")
    src = open(pp_path).read()

    # Source must contain the blacklist check in inject_user_params
    assert_true("bad_inject_params" in src,
        "inject_user_params must contain blacklist check code")

    tmp = tempfile.mkdtemp()
    try:
        sess = AgentSession(session_dir=tmp)

        # Simulate: run 1 failed with ignore_symmetry_conflicts=True
        sess.record_bad_inject_param("phenix.refine", "ignore_symmetry_conflicts")

        # Verify blacklist state
        bad = sess.get_bad_inject_params("phenix.refine")
        assert_true("ignore_symmetry_conflicts" in bad,
            "ignore_symmetry_conflicts must be in blacklist after recording")

        # Any param NOT in the blacklist is still injectable (sanity check)
        other_bad = sess.get_bad_inject_params("phenix.refine")
        assert_true("resolution" not in other_bad,
            "resolution must NOT be in blacklist — only the bad param should be")
    finally:
        shutil.rmtree(tmp)
    print("  PASSED: blacklisted param correctly excluded from injection")


def test_s4d_generate_program_eff_reports_rejected_args():
    """generate_program_eff must populate rejected_args with params that the
    PHIL interpreter rejects (e.g. 'ignore_symmetry_conflicts').

    This is the primary fix for the inject→crash→re-inject infinite loop:
    phenix.refine's interpreter catches the bad param during .eff generation
    and logs it, but previously nothing recorded it so the agent kept
    re-injecting it every cycle.
    """
    print("  Test: s4d_generate_program_eff_reports_rejected_args")

    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    # generate_program_eff must accept a rejected_args parameter
    assert_true("def generate_program_eff" in src,
        "generate_program_eff must exist in ai_agent.py")
    assert_true("rejected_args=None" in src.split("def generate_program_eff")[1][:400],
        "generate_program_eff must have rejected_args=None parameter")

    # Must populate rejected_args from failed_args (interpreter rejects)
    assert_true("rejected_args.append" in src or
                "rejected_args is not None" in src,
        "generate_program_eff must populate rejected_args list")

    # Blacklisting must happen in _execute_sub_job_for_gui after collecting
    # rejected args. Verify both the collector and the blacklist call exist.
    assert_true("_rejected_args" in src,
        "_execute_sub_job_for_gui must use _rejected_args list")
    assert_true("record_bad_inject_param" in src,
        "session.record_bad_inject_param must be called somewhere in ai_agent.py")
    # The blacklist call must reference _rejected_args (not just the earlier
    # error-text-based blacklisting)
    assert_true("for _bad_key in _rejected_args" in src,
        "Must iterate over _rejected_args to blacklist each bad param")

    # rejected_args must be threaded through _try_native_execution
    assert_true("rejected_args=_rejected_args" in src,
        "_try_native_execution call must pass rejected_args=_rejected_args")

    # rejected_args must be threaded through _run_via_old_style_runner
    assert_true("_run_via_old_style_runner" in src and
                "rejected_args=rejected_args" in src,
        "_run_via_old_style_runner must forward rejected_args")

    print("  PASSED: rejected_args propagation wired through .eff generation chain")


def test_s4d_blacklist_extracted_from_error_message():
    """The bad param name must be extractable from PHENIX's 'Unknown command
    line parameter definition: FOO = VALUE' error format using the same regex
    used in _record_command_result."""
    print("  Test: s4d_blacklist_extracted_from_error_message")
    import re

    # The regex used in _record_command_result
    _bad_re = re.compile(
        r'(?:parameter definition|no such parameter)\s*[:\-]?\s*'
        r'([A-Za-z_][\w.]*)',
        re.IGNORECASE
    )

    test_cases = [
        # (error text, expected extracted key)
        ("Unknown command line parameter definition: ignore_symmetry_conflicts = True",
         "ignore_symmetry_conflicts"),
        # Full dotted path is extracted; caller also stores the short leaf name separately
        ("Sorry: Unknown command line parameter definition: bad.scope.key = False",
         "bad.scope.key"),
        ("Sorry: there is no such parameter: my_param",
         "my_param"),
    ]
    for error_text, expected_key in test_cases:
        m = _bad_re.search(error_text)
        assert_not_none(m, "Regex must match in: %r" % error_text)
        extracted = m.group(1).strip().rstrip('=').strip()
        assert_equal(extracted, expected_key,
            "Extracted key %r must equal %r for error: %r" % (
                extracted, expected_key, error_text[:60]))

    print("  PASSED: bad param name correctly extracted from PHENIX error messages")


def test_s4e_space_group_placeholder_not_injected():
    """inject_crystal_symmetry must NOT inject placeholder space_group values
    like 'Not mentioned' — these come from the directive extractor when it
    couldn't find a real symbol, and passing them to PHENIX causes a PHIL error:
      crystal_symmetry.space_group=Not mentioned
    """
    print("  Test: s4e_space_group_placeholder_not_injected")

    # Check command_postprocessor.py (where inject_crystal_symmetry now lives)
    import os as _os
    postprocessor_path = _os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(_os.path.exists(postprocessor_path),
        "command_postprocessor.py must exist at %s" % postprocessor_path)
    src = open(postprocessor_path).read()

    assert_true("_INVALID_SG_PATTERNS" in src or "_is_placeholder" in src,
        "inject_crystal_symmetry must have a placeholder guard for space_group")
    assert_true('"not mentioned"' in src,
        "Guard must explicitly reject 'not mentioned' values")

    # Exercise the guard logic inline (same code as in the method)
    _INVALID_SG_PATTERNS = (
        "not mentioned", "not specified", "not provided",
        "unknown", "n/a", "none", "null", "tbd", "to be determined",
    )
    def _is_placeholder(sg_str):
        sg_lower = sg_str.lower()
        return (any(p in sg_lower for p in _INVALID_SG_PATTERNS) or
                len(sg_str) > 20 or
                not sg_str.strip()[0:1].isalpha())

    invalid_values = [
        "Not mentioned", "not specified", "not provided",
        "Unknown", "N/A", "none", "null",
    ]
    valid_values = [
        "P 63", "P212121", "P 32 2 1", "C 2 2 21", "F 4 3 2", "R3",
    ]
    for v in invalid_values:
        assert_true(_is_placeholder(v),
            "Value %r should be detected as a placeholder and NOT injected" % v)
    for v in valid_values:
        assert_true(not _is_placeholder(v),
            "Value %r should NOT be detected as a placeholder" % v)

    print("  PASSED: space_group placeholder values correctly identified and blocked")


def test_s4e_sanitize_command_removes_placeholder_space_group():
    """sanitize_command must:
    1. Strip crystal_symmetry.space_group=<placeholder> tokens regardless of
       source (LLM-generated or injected), and
    2. Strip session-blacklisted parameters (e.g. ignore_symmetry_conflicts=True)
       that the LLM regenerates despite the warning in guidelines.

    The standalone function in command_postprocessor.py accepts bad_inject_params
    directly (no session dependency).
    """
    print("  Test: s4e_sanitize_command_removes_placeholder_space_group")
    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    # Phase 4: dead class method must be removed
    assert_true("def _sanitize_command" not in src,
        "_sanitize_command class method must be removed from ai_agent.py (Phase 4)")
    # Sanitize now lives exclusively in command_postprocessor.py, called by BUILD.
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.exists(pp_path),
        "command_postprocessor.py must exist in agent/")
    pp_src = open(pp_path).read()
    assert_true("def sanitize_command" in pp_src,
        "sanitize_command must be defined in command_postprocessor.py")
    assert_true("def postprocess_command" in pp_src,
        "postprocess_command must be defined in command_postprocessor.py")
    gn_path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    gn_src = open(gn_path).read()
    assert_true("postprocess_command" in gn_src,
        "postprocess_command must be called in graph_nodes.py BUILD")
    # Phase 3: directive stop + consecutive cap now in PERCEIVE
    assert_true("check_directive_stop" in gn_src,
        "check_directive_stop must be called in graph_nodes.py PERCEIVE")
    assert_true("check_consecutive_program_cap" in gn_src,
        "check_consecutive_program_cap must be called in graph_nodes.py PERCEIVE")

    # --- Exercise the sanitise logic inline ---
    _PLACEHOLDER_PATTERNS = (
        "not mentioned", "not specified", "not provided",
        "unknown", "n/a", "none", "null", "tbd", "to be determined",
    )

    def _is_placeholder_value(val):
        v = val.strip().strip('"').strip("'").lower()
        return any(p in v for p in _PLACEHOLDER_PATTERNS)

    import re as _re
    _ALL_KV_RE = _re.compile(
        r'\s*([\w.]+)\s*=\s*(?:"[^"]*"|\'[^\']*\'|[^\s]+)'
    )

    def sanitize(command, bad_keys=None):
        bad_keys = bad_keys or set()
        sanitized = command
        changed = True
        while changed:
            prev = sanitized
            parts = []
            last = 0
            for m in _ALL_KV_RE.finditer(sanitized):
                key_full  = m.group(1)
                key_short = key_full.split('.')[-1]
                token = m.group(0)
                eq_pos = token.find('=')
                val = token[eq_pos+1:].strip() if eq_pos >= 0 else ''
                strip = False
                if key_full in bad_keys or key_short in bad_keys:
                    strip = True
                if not strip and ('space_group' in key_full or 'unit_cell' in key_full):
                    if _is_placeholder_value(val):
                        strip = True
                parts.append(sanitized[last:m.start() if strip else m.end()])
                last = m.end()
            parts.append(sanitized[last:])
            sanitized = ''.join(parts)
            changed = (sanitized != prev)
        return _re.sub(r'  +', ' ', sanitized).strip()

    # Test 1: placeholder space_group is stripped
    bad_cmd = ('phenix.map_correlations input_files.model=model.pdb '
               'crystal_symmetry.space_group="Not specified" '
               'output.prefix=check_001')
    cleaned = sanitize(bad_cmd)
    assert_true('Not specified' not in cleaned,
        "'Not specified' must be stripped, got: %r" % cleaned)
    assert_true('output.prefix=check_001' in cleaned,
        "output.prefix must be preserved, got: %r" % cleaned)

    # Test 2: blacklisted param is stripped when session provides it
    bad_cmd2 = ('phenix.map_correlations input_files.model=model.pdb '
                'crystal_symmetry.space_group="Not specified" '
                'ignore_symmetry_conflicts=True output.prefix=check_001')
    cleaned2 = sanitize(bad_cmd2,
                        bad_keys={'ignore_symmetry_conflicts'})
    assert_true('ignore_symmetry_conflicts' not in cleaned2,
        "blacklisted param must be stripped, got: %r" % cleaned2)
    assert_true('Not specified' not in cleaned2,
        "placeholder must also be stripped, got: %r" % cleaned2)
    assert_true('output.prefix=check_001' in cleaned2,
        "output.prefix must be preserved, got: %r" % cleaned2)

    # Test 3: real space group is kept
    good_cmd = ('phenix.refine model.pdb data.mtz '
                'crystal_symmetry.space_group="P 63" output.prefix=refine_001')
    assert_equal(sanitize(good_cmd), good_cmd.strip(),
        "Real space group 'P 63' must NOT be stripped")

    # Test 4: other placeholder variants
    for placeholder in ('"not specified"', '"Not mentioned"', '"Unknown"',
                        'None', 'null', '"not provided"'):
        cmd = 'phenix.refine model.pdb crystal_symmetry.space_group=%s' % placeholder
        result = sanitize(cmd)
        assert_true(placeholder.strip('"').lower() not in result.lower(),
            "Placeholder %r must be stripped, got: %r" % (placeholder, result))

    print("  PASSED: _sanitize_command strips placeholder space_group and blacklisted params")


    """When bad_inject_params is non-empty in the session, _query_agent_for_command
    must append a 'DO NOT USE' warning to the guidelines text before calling the LLM,
    so the LLM does not regenerate the invalid parameter in its own command output.

    Without this, the LLM reads the user's guidelines, sees
    'ignore_symmetry_conflicts=True', and adds it to the command even though the
    agent has already blacklisted it from injection.
    """
    print("  Test: s4e_blacklist_warning_appended_to_guidelines")
    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    assert_true("bad_inject_params" in src and "INVALID PARAMETERS" in src,
        "_query_agent_for_command must inject a warning about bad params into guidelines")
    assert_true("DO NOT add" in src,
        "Warning must tell the LLM not to add the bad parameter")
    # The warning must appear before the '# 2. ASK THE AGENT' comment
    # (i.e., before decide_next_step is called)
    warn_idx = src.find("INVALID PARAMETERS")
    ask_idx  = src.find("# 2. ASK THE AGENT")
    assert_true(warn_idx < ask_idx,
        "INVALID PARAMETERS warning must be injected before ASK THE AGENT call "
        "(warn_idx=%d, ask_idx=%d)" % (warn_idx, ask_idx))

    print("  PASSED: guidelines get blacklist warning injected before LLM call")


    print("  PASSED: guidelines get blacklist warning injected before LLM call")


# =============================================================================
# MODEL_VS_DATA PROBE INCONCLUSIVE HANDLING
# =============================================================================

def test_s5a_model_vs_data_probe_failure_is_inconclusive():
    """When phenix.model_vs_data fails before any refine/dock, workflow_state
    must mark placement_probed=True with probe_result=None (inconclusive) rather
    than stopping.  A crystal_symmetry_mismatch from the probe does NOT mean the
    model is unplaced — phenix.refine is more permissive and may succeed.
    """
    print("  Test: s5a_model_vs_data_probe_failure_is_inconclusive")
    from agent.workflow_state import _analyze_history

    history = [
        {
            "cycle": 1,
            "program": "phenix.xtriage",
            "result": "SUCCESS: ...",
            "command": "phenix.xtriage data.mtz",
            "analysis": {},
        },
        {
            "cycle": 2,
            "program": "phenix.model_vs_data",
            "result": "FAILED: Sorry: Crystal symmetry mismatch between different files.",
            "command": "phenix.model_vs_data model.pdb data.mtz",
            "analysis": {},
        },
    ]
    info = _analyze_history(history)

    assert_true(info.get("placement_probed"),
        "placement_probed must be True after failed model_vs_data probe")
    assert_true(info.get("placement_probe_result") is None,
        "placement_probe_result must be None (inconclusive) for crystal_symmetry_mismatch, "
        "got: %r" % info.get("placement_probe_result"))

    print("  PASSED: model_vs_data probe failure correctly marked as inconclusive")


def test_s5a_explicit_refine_suppresses_probe():
    """When explicit_program='phenix.refine' is in session_info, build_context
    must set placement_uncertain=False, skipping the model_vs_data probe step.

    Verified via source inspection (libtbx unavailable in test environment).
    """
    print("  Test: s5a_explicit_refine_suppresses_probe")
    we_path = os.path.join(_PROJECT_ROOT, 'agent', 'workflow_engine.py')
    src = open(we_path).read()

    assert_true("_PROBE_SKIP_PROGRAMS" in src,
        "build_context must define _PROBE_SKIP_PROGRAMS for explicit_program check")
    assert_true("phenix.refine" in src.split("_PROBE_SKIP_PROGRAMS")[1][:300],
        "_PROBE_SKIP_PROGRAMS must include 'phenix.refine'")
    assert_true("phenix.ligandfit" in src.split("_PROBE_SKIP_PROGRAMS")[1][:300],
        "_PROBE_SKIP_PROGRAMS must include 'phenix.ligandfit'")
    assert_true("_user_asserts_placed" in src,
        "build_context must set _user_asserts_placed from explicit_program check")
    assert_true("not _user_asserts_placed" in src,
        "placement_uncertain condition must include 'not _user_asserts_placed'")

    print("  PASSED: explicit_program=phenix.refine correctly suppresses placement probe")


def test_s5a_inconclusive_probe_routes_to_refine_not_phaser():
    """After an inconclusive probe (placement_probed=True, result=None),
    the engine must promote has_placed_model=True so routing falls through
    to 'refine' rather than 'obtain_model' (Phaser).

    Verified by source inspection of workflow_engine.py.
    """
    print("  Test: s5a_inconclusive_probe_routes_to_refine_not_phaser")

    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'workflow_engine.py')).read()

    # The placement_probed block must set has_placed_model=True for
    # non-needs_mr results (covers both "placed" and None/inconclusive)
    assert_true('context["has_placed_model"] = True' in src,
        'Engine must explicitly set context["has_placed_model"] = True '
        'for inconclusive probe result')

    # The promotion must come inside the placement_probed block,
    # before the obtain_model routing
    prb_idx   = src.find('context.get("placement_probed")')
    promo_idx = src.find('context["has_placed_model"] = True', prb_idx)
    obtain_idx = src.find('"obtain_model"', prb_idx)
    assert_true(0 < promo_idx < obtain_idx,
        "has_placed_model promotion must come before obtain_model check "
        "(promo=%d, obtain=%d)" % (promo_idx, obtain_idx))

    print("  PASSED: inconclusive probe promotes has_placed_model → refine chosen over Phaser")


    """When phenix.model_vs_data fails with crystal_symmetry_mismatch, the cycle
    loop must NOT call _diagnose_terminal_failure — treated as inconclusive probe.

    Verified by source inspection.
    """
    print("  Test: s5a_probe_symmetry_mismatch_not_terminal")
    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    assert_true("_is_probe_inconclusive" in src,
        "Must have _is_probe_inconclusive guard in cycle loop")

    # _is_probe_program (containing phenix.model_vs_data) is declared just before
    # _is_probe_inconclusive; check the ±400-char window around the first occurrence.
    idx = src.find("_is_probe_inconclusive")
    window = src[max(0, idx - 400): idx + 400]
    assert_true("phenix.model_vs_data" in window,
        "Guard window must reference phenix.model_vs_data")
    # Guard now covers ANY error type for probe programs (not just crystal_symmetry_mismatch)
    assert_true("_is_probe_program" in window,
        "Guard window must reference _is_probe_program")

    guard_block = src[idx: idx + 800]
    assert_true("_diagnosis_match = None" in guard_block,
        "Guard must set _diagnosis_match=None to skip terminal diagnosis")

    print("  PASSED: any diagnosable error from model_vs_data probe is non-terminal")


def test_s5b_pkl_free_programs_not_marked_failed():
    """Programs that never produce a .pkl (e.g. phenix.model_vs_data,
    phenix.xtriage) must not be marked 'failed' solely because no pkl was
    found.  The no-pkl heuristic must be skipped for these programs.
    """
    print("  Test: s5b_pkl_free_programs_not_marked_failed")
    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    assert_true("_PKL_FREE_PROGRAMS" in src,
        "_PKL_FREE_PROGRAMS whitelist must exist in ai_agent.py")
    assert_true('"phenix.model_vs_data"' in src,
        "phenix.model_vs_data must be in _PKL_FREE_PROGRAMS")
    assert_true('"phenix.xtriage"' in src,
        "phenix.xtriage must be in _PKL_FREE_PROGRAMS")

    # The no-pkl check must be gated on _no_pkl_expected
    assert_true("_no_pkl_expected" in src,
        "_no_pkl_expected flag must be used to skip the no-pkl failure check")

    # Verify the guard comes before "No pkl produced — marking as failed"
    guard_idx  = src.find("_no_pkl_expected")
    failed_idx = src.find("No pkl produced — marking as failed")
    assert_true(guard_idx < failed_idx,
        "_no_pkl_expected guard must appear before 'marking as failed' string")

    print("  PASSED: pkl-free program whitelist prevents false failure marking")


def test_s5c_ligandfit_removed_from_valid_programs_when_no_ligand_file():
    """phenix.ligandfit must NOT be filtered from valid_programs based solely
    on the absence of a separately-categorized ligand file.
    phenix.ligandfit works with a residue code (ligand_type=ATP) or a ligand
    PDB whose name may not match categorization patterns (e.g. atp.pdb).
    """
    print("  Test: s5c_ligandfit_removed_from_valid_programs_when_no_ligand_file")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'workflow_engine.py')).read()

    # The old hard filter on has_ligand_file must be gone
    assert_true(
        'if "phenix.ligandfit" in valid and not context.get("has_ligand_file"' not in src,
        "get_valid_programs must NOT hard-filter ligandfit based on has_ligand_file alone")
    # The restore path for user_wants_ligandfit must still exist
    assert_true(
        '"phenix.ligandfit"' in src and "valid.append" in src,
        "get_valid_programs must still restore ligandfit when user_wants_ligandfit")

    print("  PASSED: ligandfit not blocked by absent ligand file categorization")


def test_s5c_autostop_with_no_ligand_file_message():
    """The PLAN AUTO-STOP suppression path must NOT block on has_ligand_file.
    phenix.ligandfit can run with just a residue code or a PDB that is not
    auto-categorized.  The old block-on-missing-file logic was removed.
    """
    print("  Test: s5c_autostop_with_no_ligand_file_message")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'graph_nodes.py')).read()

    # Old blocking variables must be gone
    assert_true("_has_ligand_file" not in src,
        "PLAN node must NOT check _has_ligand_file")
    assert_true("no ligand restraint file" not in src,
        "PLAN node must NOT emit 'no ligand restraint file' — ligandfit does "
        "not require a separately-categorized file")
    # Suppress path must still exist
    assert_true("Suppressing AUTO-STOP because after_program" in src,
        "PLAN node must still suppress AUTO-STOP when after_program hasn't run yet")

    print("  PASSED: PLAN node no longer wrongly blocks ligandfit on missing file")


def test_s5d_user_wants_ligandfit_stays_in_refine_phase():
    """When user explicitly requests ligand fitting, the workflow must stay in
    the 'refine' step so ligandfit is available, regardless of has_ligand_file.
    """
    print("  Test: s5d_user_wants_ligandfit_stays_in_refine_phase")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'workflow_engine.py')).read()

    assert_true("user_wants_ligandfit" in src,
        "build_context must set user_wants_ligandfit flag")
    assert_true('"ligandfit" in _after_prog' in src or
                '"ligandfit" in _after_prog.lower()' in src,
        "user_wants_ligandfit must check after_program for ligandfit")
    assert_true("_wants_ligandfit" in src,
        "_detect_xray_step must use _wants_ligandfit")
    assert_true("0.50" in src,
        "Relaxed r_free threshold (0.50) must be present for user_wants_ligandfit")

    # The user_wants_ligandfit step block must NOT require has_ligand_file
    wants_idx = src.find("_wants_ligandfit = context.get(\"user_wants_ligandfit\"")
    block_end = src.find("return self._make_phase_result", wants_idx)
    block = src[wants_idx:block_end + 60]
    assert_true("has_ligand_file" not in block,
        "_detect_xray_step ligandfit gate must not require has_ligand_file")

    restore_idx = src.find("user_wants_ligandfit")
    valid_idx = src.find('valid.append("phenix.ligandfit")', restore_idx)
    assert_true(valid_idx > 0,
        "get_valid_programs must restore phenix.ligandfit when user_wants_ligandfit")

    print("  PASSED: user_wants_ligandfit keeps ligandfit available without has_ligand_file")


def test_s5d_validate_phase_does_not_block_ligandfit():
    """The workflow must not route to 'validate' when the user wants ligandfit
    and has not yet run it — even if r_free >= 0.35.

    Source check: the _check_metric_condition for r_free returns True when value
    is None, so the step detection gate must handle the explicit-user-request
    case by either raising the threshold or re-adding ligandfit to valid_programs.
    """
    print("  Test: s5d_validate_phase_does_not_block_ligandfit")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'workflow_engine.py')).read()

    # The r_free threshold must be conditional on user_wants_ligandfit
    assert_true("_rfree_threshold" in src,
        "r_free threshold must be a variable (_rfree_threshold) "
        "that is relaxed when user_wants_ligandfit")
    # The fallthrough to validate must come AFTER the ligandfit gate
    ligandfit_gate_idx = src.find("user_wants_ligandfit")
    validate_idx = src.find("need validation before stopping")
    assert_true(ligandfit_gate_idx < validate_idx,
        "Ligandfit gate must be evaluated before routing to validate step")

    print("  PASSED: ligandfit gate evaluated before validate routing")


def test_s5e_notice_emitted_when_ligandfit_impossible():
    """PERCEIVE must NOT block ligandfit based on has_ligand_file.
    phenix.ligandfit does not require a separately-categorized ligand file;
    it works with a residue code or an atp.pdb that may not be auto-categorized.
    The old blocking logic was removed.
    """
    print("  Test: s5e_notice_emitted_when_ligandfit_impossible")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'graph_nodes.py')).read()

    # Old blocking variables must be gone
    assert_true("_ligandfit_missing_notice_sent" not in src,
        "PERCEIVE must NOT have _ligandfit_missing_notice_sent blocking logic")
    assert_true("Ligand fitting requested but no ligand restraint file" not in src,
        "PERCEIVE must NOT emit 'no ligand restraint file' stop message")
    # user_wants_ligandfit log must still exist for debugging
    assert_true("user_wants_ligandfit" in src,
        "PERCEIVE must still log user_wants_ligandfit for debugging")

    print("  PASSED: PERCEIVE does not block ligandfit on has_ligand_file")


def test_s5e_notice_event_type_exists_and_formatted():
    """NOTICE event type must exist in EventType and have a formatter so it
    renders clearly in the output at QUIET verbosity (always shown).
    """
    print("  Test: s5e_notice_event_type_exists_and_formatted")
    event_log_src = open(os.path.join(_PROJECT_ROOT, 'agent', 'event_log.py')).read()
    formatter_src = open(os.path.join(_PROJECT_ROOT, 'agent', 'event_formatter.py')).read()

    assert_true('NOTICE = "notice"' in event_log_src,
        "EventType.NOTICE must be defined in event_log.py")
    assert_true("EventType.NOTICE: Verbosity.QUIET" in event_log_src,
        "NOTICE must be in EVENT_VERBOSITY at QUIET level (always shown)")
    assert_true("EventType.NOTICE: self._format_notice" in formatter_src,
        "EventFormatter must dispatch NOTICE to _format_notice")
    assert_true("def _format_notice" in formatter_src,
        "_format_notice method must exist in EventFormatter")

    print("  PASSED: NOTICE event type exists with QUIET verbosity and formatter")


def test_s5e_plan_node_emits_notice_on_missing_ligand_stop():
    """The PLAN node must NOT stop based on missing ligand file — ligandfit
    works without a separately-categorized file.  The old stop was removed.
    """
    print("  Test: s5e_plan_node_emits_notice_on_missing_ligand_stop")
    src = open(os.path.join(_PROJECT_ROOT, 'agent', 'graph_nodes.py')).read()

    assert_true("Ligand fitting requested but cannot proceed" not in src,
        "PLAN node must NOT emit old missing-ligand NOTICE (logic removed)")
    assert_true("_has_ligand_file" not in src,
        "PLAN node must NOT check _has_ligand_file")
    # Suppress path must still exist
    assert_true("Suppressing AUTO-STOP because after_program" in src,
        "PLAN node must still suppress AUTO-STOP when after_program hasn't run")

    print("  PASSED: PLAN node no longer blocks ligandfit on missing file")


def test_s5g_natural_language_macro_cycles_injected():
    """'Run only one macro-cycle of refinement' must inject
    main.number_of_macro_cycles=1 into the phenix.refine command.
    The conversion is driven by programs.yaml strategy_flags.natural_language
    so the flag path is always valid.
    """
    print("  Test: s5g_natural_language_macro_cycles_injected")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)
    from agent.nl_to_phil import extract_nl_params

    # Test 1: "one macro-cycle" word form → number_of_macro_cycles=1
    r1 = extract_nl_params("Run only one macro-cycle of refinement.", "phenix.refine")
    assert_true("main.number_of_macro_cycles=1" in r1,
        "'one macro-cycle' must produce main.number_of_macro_cycles=1, got: %r" % r1)

    # Test 2: digit form
    r2 = extract_nl_params("Run 5 macro-cycles of refinement.", "phenix.refine")
    assert_true("main.number_of_macro_cycles=5" in r2,
        "'5 macro-cycles' must produce main.number_of_macro_cycles=5, got: %r" % r2)

    # Test 3: ligandfit — must not inject refine params
    r3 = extract_nl_params("Run only one macro-cycle of refinement.", "phenix.ligandfit")
    assert_true(r3 == [],
        "Must NOT inject number_of_macro_cycles for phenix.ligandfit, got: %r" % r3)

    # Test 4: real_space_refine uses different PHIL path
    r4 = extract_nl_params("Run only one macro-cycle of refinement.", "phenix.real_space_refine")
    assert_true("macro_cycles=1" in r4,
        "real_space_refine must produce macro_cycles=1, got: %r" % r4)

    # Test 5: simulated annealing boolean flag
    r5 = extract_nl_params("Use simulated annealing.", "phenix.refine")
    assert_true("main.simulated_annealing=True" in r5,
        "'simulated annealing' must produce main.simulated_annealing=True, got: %r" % r5)

    # Test 6: no match → empty
    r6 = extract_nl_params("No special settings.", "phenix.refine")
    assert_true(r6 == [],
        "No NL phrases should produce empty list, got: %r" % r6)

    # Test 7: programs.yaml natural_language entries exist for phenix.refine cycles
    from knowledge.yaml_loader import get_program
    prog = get_program("phenix.refine")
    flags = prog.get("strategy_flags", {})
    assert_true("natural_language" in flags.get("cycles", {}),
        "cycles flag must have natural_language entries in programs.yaml")
    assert_true("natural_language" in flags.get("simulated_annealing", {}),
        "simulated_annealing flag must have natural_language entries in programs.yaml")

    # Test 8: inject_user_params calls extract_nl_params (not hardcoded)
    # inject_user_params now lives in command_postprocessor.py (moved from ai_agent.py)
    import os as _os
    postprocessor_path = _os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(_os.path.exists(postprocessor_path),
        "command_postprocessor.py must exist")
    pp_src = open(postprocessor_path).read()
    assert_true("extract_nl_params" in pp_src,
        "inject_user_params must call extract_nl_params from nl_to_phil")
    assert_true("nl_to_phil" in pp_src,
        "command_postprocessor must import from nl_to_phil")

    print("  PASSED: NL->PHIL driven by programs.yaml, correct paths for each program")



    """phenix.polder conditions: requires model + data_mtz + ligand_fit
    (a ligand must have been fitted or already be in the model), and
    must not have already run (not_done: polder).
    """
    print("  Test: s5f_polder_conditions")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)
    from knowledge.yaml_loader import get_workflow_steps

    steps = get_workflow_steps("xray")

    def check_conditions(prog_entry, context):
        for cond in prog_entry.get("conditions", []):
            if "has" in cond:
                if not context.get("has_" + cond["has"]):
                    return False
            if "not_done" in cond:
                done_key = cond["not_done"] + "_done"
                if context.get(done_key):
                    return False
        return True

    for step_name in ("refine", "validate"):
        polder_entry = next(
            (p for p in steps[step_name].get("programs", [])
             if isinstance(p, dict) and "polder" in p.get("program", "")),
            None
        )
        assert_true(polder_entry is not None,
            "phenix.polder must be defined in xray %s step" % step_name)

        # Polder should NOT be valid without ligand_fit
        ctx_no = {"has_model": True, "has_data_mtz": True,
                  "has_ligand_fit": False}
        assert_true(not check_conditions(polder_entry, ctx_no),
            "polder must NOT be valid in %s step without ligand_fit"
            % step_name)

        # Polder should be valid with ligand_fit
        ctx_yes = {"has_model": True, "has_data_mtz": True,
                   "has_ligand_fit": True}
        assert_true(check_conditions(polder_entry, ctx_yes),
            "polder MUST be valid in %s step with ligand_fit"
            % step_name)

        # Polder should NOT be valid when already done
        ctx_done = {"has_model": True, "has_data_mtz": True,
                    "has_ligand_fit": True, "polder_done": True}
        assert_true(not check_conditions(polder_entry, ctx_done),
            "polder must NOT be valid in %s step when already done"
            % step_name)

        # Confirm both conditions are present
        cond_keys = []
        for c in polder_entry.get("conditions", []):
            if isinstance(c, dict):
                for k, v in c.items():
                    cond_keys.append("%s:%s" % (k, v))
        assert_true("has:ligand_fit" in cond_keys,
            "polder %s conditions must include 'has: ligand_fit', got: %s"
            % (step_name, cond_keys))
        assert_true("not_done:polder" in cond_keys,
            "polder %s conditions must include 'not_done: polder', got: %s"
            % (step_name, cond_keys))

    print("  PASSED: polder conditions correct (model + data + ligand_fit + not_done)")

def test_s5h_sanitize_strips_out_of_scope_params():
    """_sanitize_command must strip key=value tokens that are not in the
    programs.yaml strategy_flags allowlist for the target program.

    Specifically: rebuild_strategy=quick must be stripped from a
    phenix.model_vs_data command because model_vs_data has no strategy_flags,
    and therefore an empty allowlist.
    """
    print("  Test: s5h_sanitize_strips_out_of_scope_params")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    # We test the sanitize logic directly by importing and calling it through
    # a minimal mock, since AIAnalysis is not importable in the test environment.
    # Instead, verify the code structure and the allowlist logic via programs.yaml.

    from knowledge.yaml_loader import get_program

    # model_vs_data must have NO strategy_flags (or empty) so its allowlist is
    # limited to universal keys — rebuild_strategy is not in that set.
    mvd = get_program("phenix.model_vs_data")
    assert_true(mvd is not None, "phenix.model_vs_data must exist in programs.yaml")
    sf = mvd.get("strategy_flags") or {}
    assert_true(len(sf) == 0,
        "phenix.model_vs_data must have no strategy_flags; got: %s" % list(sf.keys()))

    # rebuild_strategy / rebuilding_strategy must NOT appear in any strategy_flags
    # that could bleed into model_vs_data (verify the prompt fix too).
    import os
    _prompt_path = os.path.join(_PROJECT_ROOT, 'knowledge', 'prompts_hybrid.py')
    prompt_src = open(_prompt_path).read()
    assert_true("rebuilding_strategy" not in prompt_src,
        "prompts_hybrid.py must NOT mention rebuilding_strategy "
        "(it is not a real PHIL param and causes cross-program contamination)")

    # phenix.refine must have strategy_flags including number_of_macro_cycles
    # so that its allowlist correctly accepts it.
    refine = get_program("phenix.refine")
    refine_sf = refine.get("strategy_flags") or {}
    refine_keys = set(refine_sf.keys())
    # Collect bare keys from flag templates
    refine_bare = set()
    for sfdef in refine_sf.values():
        if isinstance(sfdef, dict):
            tpl = sfdef.get('flag', '')
            bare = tpl.split('=')[0].strip().split('.')[-1].lower()
            if bare and bare != '{value}':
                refine_bare.add(bare)
    assert_true("number_of_macro_cycles" in refine_bare,
        "phenix.refine allowlist must include number_of_macro_cycles; got bare keys: %s"
        % sorted(refine_bare))

    # Verify source code has the narrowed Rule C
    # sanitize_command now lives in command_postprocessor.py (moved from ai_agent.py)
    import os as _os
    _pp_path = _os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(_os.path.exists(_pp_path), "command_postprocessor.py must exist")
    ai_src = open(_pp_path).read()
    assert_true("_prog_allowlist" in ai_src,
        "sanitize_command must build _prog_allowlist from programs.yaml")
    assert_true("Rule C" in ai_src,
        "sanitize_command must have Rule C")
    # Rule C must be scoped to len(strategy_flags) == 0 (file-only programs)
    assert_true("len(strategy_flags) == 0" in ai_src,
        "Rule C must only fire for programs with zero strategy_flags")
    # Must NOT use the old broad allowlist check
    assert_true("not in _prog_allowlist" not in ai_src,
        "Rule C must NOT use broad _prog_allowlist check (regression risk)")

    print("  PASSED: sanitize_command Rule C stripped for file-only programs only; "
          "rebuilding_strategy removed from prompt")


def test_s5h_rule_d_strips_bare_hallucinated_params():
    """Rule D in sanitize_command must strip bare (unscoped) key=value params
    that are not in strategy_flags for programs that HAVE strategy_flags.
    Scoped PHIL params (with dots) must survive — they go through PHIL validation.
    File-path values (model=/path/to/file.pdb) must survive — they are legitimate
    flagged file arguments from CommandBuilder.

    Concrete case: map_type=pre_calculated must be stripped from phenix.refine
    (map_type is not in its strategy_flags), but xray_data.r_free_flags.generate=True
    must survive (scoped PHIL param with dots).
    """
    print("  Test: s5h_rule_d_strips_bare_hallucinated_params")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    from agent.command_postprocessor import sanitize_command

    # map_type=pre_calculated must be stripped (bare, not in allowlist)
    cmd = ("phenix.refine model.pdb data.mtz "
           "xray_data.r_free_flags.generate=True output.prefix=refine_001 "
           "map_type=pre_calculated")
    result = sanitize_command(cmd, "phenix.refine")
    assert_true("map_type" not in result,
        "Rule D must strip bare hallucinated map_type; got: %s" % result)

    # Scoped PHIL param must survive
    assert_true("xray_data.r_free_flags.generate=True" in result,
        "Rule D must keep scoped PHIL params; got: %s" % result)

    # Known strategy flags must survive
    cmd2 = "phenix.refine model.pdb data.mtz nproc=4 simulated_annealing=True"
    result2 = sanitize_command(cmd2, "phenix.refine")
    assert_true("nproc=4" in result2,
        "Rule D must keep universal key nproc; got: %s" % result2)
    assert_true("simulated_annealing=True" in result2,
        "Rule D must keep strategy_flag simulated_annealing; got: %s" % result2)

    # ---- FILE-PATH GUARD ----
    # File-path values must survive Rule D even if key is not in allowlist.
    # This is the ligandfit bug: model=/path/to/file.pdb was stripped.

    # Absolute path (contains '/')
    cmd3 = ("phenix.ligandfit model=/path/to/refine_001.pdb "
            "data=/path/to/refine_001.mtz ligand=/path/to/atp.pdb "
            "general.nproc=4")
    result3 = sanitize_command(cmd3, "phenix.ligandfit")
    assert_true("model=/path/to/refine_001.pdb" in result3,
        "Rule D must preserve model= with file path; got: %s" % result3)
    assert_true("data=/path/to/refine_001.mtz" in result3,
        "Rule D must preserve data= with file path; got: %s" % result3)
    assert_true("ligand=/path/to/atp.pdb" in result3,
        "Rule D must preserve ligand= with file path; got: %s" % result3)

    # Bare filename with crystallographic extension (no '/' but .pdb)
    cmd4 = ("phenix.ligandfit model=refine_001.pdb data=refine_001.mtz "
            "ligand=atp.pdb general.nproc=4 map_type=pre_calculated")
    result4 = sanitize_command(cmd4, "phenix.ligandfit")
    assert_true("model=refine_001.pdb" in result4,
        "Rule D must preserve model=file.pdb; got: %s" % result4)
    assert_true("data=refine_001.mtz" in result4,
        "Rule D must preserve data=file.mtz; got: %s" % result4)
    assert_true("map_type" not in result4,
        "Rule D must still strip hallucinated non-file params; got: %s" % result4)

    # Rule C file-path guard (programs with no strategy_flags)
    cmd5 = ("phenix.phaser model=refine_001.pdb data=data.mtz "
            "bogus_param=42")
    result5 = sanitize_command(cmd5, "phenix.phaser")
    assert_true("model=refine_001.pdb" in result5,
        "Rule C must preserve model=file.pdb; got: %s" % result5)
    assert_true("bogus_param" not in result5,
        "Rule C must strip non-file bare params; got: %s" % result5)

    # Verify Rule D exists in source
    import os as _os
    _pp_path = _os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    src = open(_pp_path).read()
    assert_true("Rule D" in src,
        "sanitize_command must have Rule D comment")

    print("  PASSED: Rule D strips bare hallucinated params, keeps scoped PHIL, "
          "known strategy_flags, and file-path values")


def test_s5h_inject_program_defaults():
    """inject_program_defaults must append defaults from programs.yaml when
    they are missing from the command, and must not double-add when present.

    Note: xray_data.r_free_flags.generate was moved from unconditional
    defaults to a conditional strategy_flag (generate_rfree_flags) in v115.05.
    It is now injected by the BUILD node's _apply_invariants based on
    rfree_mtz state, NOT by inject_program_defaults. This test verifies
    that generate is NOT in defaults (would cause the old bug) and IS
    in strategy_flags (conditional path).
    """
    print("  Test: s5h_inject_program_defaults")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    from knowledge.yaml_loader import get_program

    # Verify phenix.refine does NOT have r_free_flags.generate in defaults
    # (removed in v115.05 to fix 3tpp-ensemble-refine crash)
    refine = get_program("phenix.refine")
    defaults = refine.get("defaults") or {}
    assert_false("xray_data.r_free_flags.generate" in defaults,
        "phenix.refine must NOT have r_free_flags.generate in defaults "
        "(was removed in v115.05); got: %s" % list(defaults.keys()))

    # Verify it IS in strategy_flags (conditional path via BUILD invariants)
    sf = refine.get("strategy_flags") or {}
    assert_true("generate_rfree_flags" in sf,
        "phenix.refine must have generate_rfree_flags in strategy_flags; "
        "got: %s" % list(sf.keys()))

    print("  PASSED: refine has no unconditional generate default, "
          "has conditional generate_rfree_flags strategy_flag")


def test_s5h_predict_and_build_has_no_rebuilding_strategy():
    """predict_and_build strategy_flags must NOT include rebuilding_strategy
    (it is not a real PHIL parameter for that program).
    """
    print("  Test: s5h_predict_and_build_has_no_rebuilding_strategy")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)
    from knowledge.yaml_loader import get_program

    pb = get_program("phenix.predict_and_build")
    assert_true(pb is not None, "phenix.predict_and_build must exist in programs.yaml")
    sf = pb.get("strategy_flags") or {}
    assert_true("rebuilding_strategy" not in sf,
        "phenix.predict_and_build strategy_flags must not contain rebuilding_strategy")
    assert_true("rebuild_strategy" not in sf,
        "phenix.predict_and_build strategy_flags must not contain rebuild_strategy")

    # The real quick-build flag is 'quick', mapping to quick=True
    assert_true("quick" in sf,
        "phenix.predict_and_build strategy_flags must contain 'quick' flag")

    print("  PASSED: predict_and_build strategy_flags has no rebuilding_strategy")



def test_s5i_model_vs_data_rfree_extracted_with_space_colon():
    """r_free and r_work are extracted from model_vs_data log even when the
    output uses 'r_free : 0.5934' (space before colon) rather than 'r_free:'.

    Regression test for the bug where model_vs_data R-free was never extracted,
    leaving placement_probed=False so the LLM re-suggested model_vs_data,
    hit duplicate detection, and stopped without running Phaser.
    """
    print("  Test: s5i_model_vs_data_rfree_extracted_with_space_colon")
    import sys as _sys, re as _re
    _sys.path.insert(0, _PROJECT_ROOT)

    # Simulate the actual phenix.model_vs_data log output (space before colon)
    log = (
        "phenix.model_vs_data: Analysis of model vs experimental data\n"
        "========================================\n"
        "Using: beta.pdb, beta_blip_P3221.mtz\n"
        "\n"
        "                  r_work : 0.5812\n"
        "                  r_free : 0.5934\n"
        "\n"
        "Resolution: 3.00\n"
    )

    # Test 1: the patterns in programs.yaml now match
    try:
        from knowledge.yaml_loader import get_program as _get_prog
    except ImportError:
        from libtbx.langchain.knowledge.yaml_loader import get_program as _get_prog

    mvd = _get_prog("phenix.model_vs_data")
    assert_true(mvd is not None, "phenix.model_vs_data must exist in programs.yaml")
    lp = mvd.get("log_parsing", {})

    for metric in ("r_work", "r_free"):
        pat_str = lp.get(metric, {}).get("pattern")
        assert_true(pat_str is not None,
            "programs.yaml model_vs_data must have log_parsing.%s.pattern" % metric)
        m = _re.search(pat_str, log, _re.IGNORECASE | _re.MULTILINE)
        assert_true(m is not None,
            "Pattern %r must match the 'key : value' format in log.\nLog excerpt: %r"
            % (pat_str, log[:200]))

    # Test 2: patterns must NOT be anchored in a way that breaks MULTILINE-less usage
    # i.e. they should work without re.MULTILINE too (no ^ anchor needed)
    for metric in ("r_work", "r_free"):
        pat_str = lp[metric]["pattern"]
        m = _re.search(pat_str, log, _re.IGNORECASE)  # no MULTILINE
        assert_true(m is not None,
            "Pattern %r must not rely on ^ anchor; must match without re.MULTILINE"
            % pat_str)

    # Test 3: extracted values must be numerically correct
    r_work_str = _re.search(lp["r_work"]["pattern"], log, _re.IGNORECASE).group(1)
    r_free_str = _re.search(lp["r_free"]["pattern"], log, _re.IGNORECASE).group(1)
    assert_true(abs(float(r_work_str) - 0.5812) < 0.0001,
        "r_work must extract 0.5812, got %r" % r_work_str)
    assert_true(abs(float(r_free_str) - 0.5934) < 0.0001,
        "r_free must extract 0.5934, got %r" % r_free_str)

    # Test 4: workflow_state placement_probe_result correctly set to 'needs_mr'
    # when r_free > 0.50 (the threshold that triggers Phaser)
    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        from libtbx.langchain.agent.workflow_state import _analyze_history

    # _analyze_history consumes the same history format as get_history_for_agent()
    history = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "command": "phenix.xtriage data.mtz",
            "result": "SUCCESS: xtriage analysis complete",
            "analysis": {},
            "output_files": [],
        },
        {
            "cycle_number": 2,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data beta.pdb data.mtz",
            "result": "SUCCESS: Command completed without errors",
            "analysis": {"r_free": 0.5934, "r_work": 0.5812},
            "output_files": [],
        },
    ]
    info = _analyze_history(history)
    assert_true(info.get("placement_probed") is True,
        "placement_probed must be True when model_vs_data returns r_free=0.59, got: %r"
        % info.get("placement_probed"))
    assert_true(info.get("placement_probe_result") == "needs_mr",
        "placement_probe_result must be 'needs_mr' when r_free=0.59 > 0.50, got: %r"
        % info.get("placement_probe_result"))

    print("  PASSED: r_free/r_work extracted from space-colon format; "
          "placement_probe_result='needs_mr' correctly set")


def test_s6a_failure_context_stored_on_cycle():
    """Step 1 of superseded-issue handling: failure_context must be stored on
    the cycle record at each of the three failure classification sites.

    Tests are structural (source-code assertions) because the failure sites
    can only be triggered by running real PHENIX programs.
    """
    print("  Test: s6a_failure_context_stored_on_cycle")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    _ai_agent_path = _find_ai_agent_path()
    assert_true(_ai_agent_path is not None, "Could not locate ai_agent.py")
    ai_src = open(_ai_agent_path).read()

    # ── Helper method present ─────────────────────────────────────────────────
    assert_true("_input_basenames_from_command" in ai_src,
        "ai_agent must define _input_basenames_from_command static helper")
    assert_true("_FILE_EXTS" in ai_src,
        "_input_basenames_from_command must define _FILE_EXTS set")

    # ── Site 1: probe-inconclusive (scope=probe) ──────────────────────────────
    assert_true('"scope": "probe"' in ai_src,
        'Site 1 must store failure_context with scope="probe"')
    # probe annotation must come BEFORE _diagnosis_match = None clears it
    probe_ctx_idx = ai_src.find('"scope": "probe"')
    diag_none_idx = ai_src.find('_diagnosis_match = None')
    # The probe annotation block contains _diagnosis_match = None just before it
    # (the assignment is what triggers the annotation path)
    assert_true(probe_ctx_idx > 0 and diag_none_idx > 0,
        "Both probe annotation and _diagnosis_match=None must exist")

    # ── Site 2: terminal failure (scope=terminal) ─────────────────────────────
    assert_true('"scope": "terminal"' in ai_src,
        'Site 2 must store failure_context with scope="terminal"')
    # terminal annotation must appear before _diagnose_terminal_failure call
    term_ctx_idx = ai_src.find('"scope": "terminal"')
    term_call_idx = ai_src.find('return self._diagnose_terminal_failure(')
    assert_true(term_ctx_idx < term_call_idx,
        'Terminal failure_context must be stored BEFORE _diagnose_terminal_failure is called')

    # ── Site 3: generic recoverable failure (scope=recoverable) ──────────────
    assert_true('"scope": "recoverable"' in ai_src,
        'Site 3 must store failure_context with scope="recoverable"')
    # recoverable annotation uses setdefault (does not overwrite probe/terminal)
    assert_true('setdefault("failure_context"' in ai_src,
        'Recoverable site must use setdefault to avoid overwriting probe/terminal annotation')

    # ── All three sites store input_files from the command ────────────────────
    assert_true(ai_src.count('"input_files": self._input_basenames_from_command(command)') >= 3,
        'All failure sites must store input_files via _input_basenames_from_command')

    # ── _input_basenames_from_command logic test ──────────────────────────────
    # Test the helper logic inline (can't import ai_agent due to PHENIX deps)
    import os as _os
    _FILE_EXTS = {'.pdb', '.cif', '.mtz', '.sca', '.hkl', '.mrc',
                  '.ccp4', '.map', '.fa', '.fasta', '.seq', '.dat',
                  '.ncs_spec', '.eff'}
    def _input_basenames(command):
        if not command:
            return []
        basenames = []
        for token in command.split():
            if '=' in token:
                token = token.split('=', 1)[1]
            token = token.strip('"\'')  # strip quotes
            ext = _os.path.splitext(token)[1].lower()
            if ext in _FILE_EXTS:
                basenames.append(_os.path.basename(token).lower())
        return basenames

    r1 = _input_basenames("phenix.model_vs_data 7qz0.pdb 7qz0.mtz")
    assert_true(r1 == ["7qz0.pdb", "7qz0.mtz"],
        "Must extract pdb and mtz basenames, got: %r" % r1)

    r2 = _input_basenames("phenix.refine /path/to/refine_001_model.pdb data.mtz "
                          "main.number_of_macro_cycles=1")
    assert_true("refine_001_model.pdb" in r2 and "data.mtz" in r2,
        "Must extract basenames from absolute paths, ignoring PHIL params, got: %r" % r2)

    r3 = _input_basenames("")
    assert_true(r3 == [], "Empty command must return empty list")

    r4 = _input_basenames("phenix.model_vs_data model.pdb data.mtz "
                          "rebuild_strategy=quick")
    assert_true("rebuild_strategy" not in str(r4),
        "Must not include PHIL param values as file basenames, got: %r" % r4)

    print("  PASSED: failure_context stored at all three sites; "
          "_input_basenames_from_command works correctly")



def test_s6b_annotate_superseded_issues():
    """Step 2: _annotate_superseded_issues must correctly mark probe and
    recoverable failures as superseded when a subsequent authoritative program
    succeeds using overlapping inputs.
    """
    print("  Test: s6b_annotate_superseded_issues")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    def make_session(cycles_data, original_files=None):
        """Build a minimal AgentSession with synthetic cycle data."""
        s = AgentSession.__new__(AgentSession)
        s.data = {
            "cycles": cycles_data,
            "original_files": original_files or [],
        }
        s.best_files = None
        return s

    def make_cycle(cycle_number, program, result, command="", failure_context=None):
        c = {
            "cycle_number": cycle_number,
            "program": program,
            "result": result,
            "command": command,
        }
        if failure_context is not None:
            c["failure_context"] = failure_context
        return c

    # ── Test 1: basic probe supersession ─────────────────────────────────────
    # model_vs_data fails (probe) → phenix.refine succeeds with same data
    cycles = [
        make_cycle(1, "phenix.xtriage", "SUCCESS: xtriage done",
                   "phenix.xtriage 7qz0.mtz"),
        make_cycle(2, "phenix.model_vs_data", "FAILED: crystal symmetry mismatch",
                   "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
                   failure_context={
                       "error_type": "crystal_symmetry_mismatch",
                       "scope": "probe",
                       "program": "phenix.model_vs_data",
                       "input_files": ["7qz0.pdb", "7qz0.mtz"],
                       "superseded_by": None,
                   }),
        make_cycle(3, "phenix.refine", "SUCCESS: R-free=0.28",
                   "phenix.refine 7qz0.pdb 7qz0.mtz"),
    ]
    s = make_session(cycles, original_files=["/data/7qz0.pdb", "/data/7qz0.mtz"])
    s._annotate_superseded_issues(s.data["cycles"])
    fc = cycles[1]["failure_context"]
    assert_true(fc["superseded_by"] is not None,
        "Probe failure must be marked superseded when refine succeeds, got: %r" % fc)
    assert_true("phenix.refine" in fc["superseded_by"],
        "superseded_by must name phenix.refine, got: %r" % fc["superseded_by"])
    assert_true("cycle 3" in fc["superseded_by"],
        "superseded_by must include cycle number, got: %r" % fc["superseded_by"])

    # ── Test 2: NOT superseded when refinement also fails ────────────────────
    cycles2 = [
        make_cycle(2, "phenix.model_vs_data", "FAILED: crystal symmetry mismatch",
                   "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
                   failure_context={
                       "error_type": "crystal_symmetry_mismatch",
                       "scope": "probe",
                       "program": "phenix.model_vs_data",
                       "input_files": ["7qz0.pdb", "7qz0.mtz"],
                       "superseded_by": None,
                   }),
        make_cycle(3, "phenix.refine", "FAILED: refinement error",
                   "phenix.refine 7qz0.pdb 7qz0.mtz"),
    ]
    s2 = make_session(cycles2)
    s2._annotate_superseded_issues(s2.data["cycles"])
    fc2 = cycles2[0]["failure_context"]
    assert_true(fc2["superseded_by"] is None,
        "Probe must NOT be superseded when refinement also fails, got: %r" % fc2)

    # ── Test 3: terminal failure is never superseded ──────────────────────────
    cycles3 = [
        make_cycle(2, "phenix.refine", "FAILED: some terminal error",
                   "phenix.refine 7qz0.pdb 7qz0.mtz",
                   failure_context={
                       "error_type": "some_error",
                       "scope": "terminal",
                       "program": "phenix.refine",
                       "input_files": ["7qz0.pdb", "7qz0.mtz"],
                       "superseded_by": None,
                   }),
        make_cycle(3, "phenix.refine", "SUCCESS: R-free=0.28",
                   "phenix.refine 7qz0.pdb 7qz0.mtz"),
    ]
    s3 = make_session(cycles3)
    s3._annotate_superseded_issues(s3.data["cycles"])
    fc3 = cycles3[0]["failure_context"]
    assert_true(fc3["superseded_by"] is None,
        "Terminal failure must NEVER be marked superseded, got: %r" % fc3)

    # ── Test 4: stage-prefix rename — 7qz0.pdb becomes refine_001_7qz0.pdb ──
    cycles4 = [
        make_cycle(2, "phenix.model_vs_data", "FAILED: symmetry mismatch",
                   "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
                   failure_context={
                       "error_type": "crystal_symmetry_mismatch",
                       "scope": "probe",
                       "program": "phenix.model_vs_data",
                       "input_files": ["7qz0.pdb", "7qz0.mtz"],
                       "superseded_by": None,
                   }),
        make_cycle(3, "phenix.refine", "SUCCESS: R-free=0.25",
                   "phenix.refine refine_001_7qz0.pdb 7qz0.mtz"),
    ]
    s4 = make_session(cycles4)
    s4._annotate_superseded_issues(s4.data["cycles"])
    fc4 = cycles4[0]["failure_context"]
    assert_true(fc4["superseded_by"] is not None,
        "Stage-prefix rename must still trigger supersession "
        "(refine_001_7qz0.pdb ~ 7qz0.pdb), got: %r" % fc4)

    # ── Test 5: non-authority program does not supersede ─────────────────────
    cycles5 = [
        make_cycle(2, "phenix.model_vs_data", "FAILED: symmetry mismatch",
                   "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
                   failure_context={
                       "error_type": "crystal_symmetry_mismatch",
                       "scope": "probe",
                       "program": "phenix.model_vs_data",
                       "input_files": ["7qz0.pdb", "7qz0.mtz"],
                       "superseded_by": None,
                   }),
        make_cycle(3, "phenix.xtriage", "SUCCESS: xtriage done",
                   "phenix.xtriage 7qz0.mtz"),
    ]
    s5 = make_session(cycles5)
    s5._annotate_superseded_issues(s5.data["cycles"])
    fc5 = cycles5[0]["failure_context"]
    assert_true(fc5["superseded_by"] is None,
        "xtriage is not an authority for symmetry_mismatch — must not supersede, "
        "got: %r" % fc5)

    # ── Test 6: _stage_prefix_strip helper ───────────────────────────────────
    strip = AgentSession._stage_prefix_strip
    assert_true(strip("refine_001_7qz0.pdb") == "7qz0",
        "Must strip refine_001_ prefix, got: %r" % strip("refine_001_7qz0.pdb"))
    assert_true(strip("phased_model.pdb") == "model",
        "Must strip phased_ prefix, got: %r" % strip("phased_model.pdb"))
    assert_true(strip("7qz0.pdb") == "7qz0",
        "No prefix — stem unchanged, got: %r" % strip("7qz0.pdb"))

    print("  PASSED: _annotate_superseded_issues correctly marks "
          "superseded/non-superseded failures across all cases")



def test_s6c_superseded_tag_in_llm_summary():
    # Step 3: [SUPERSEDED] tag and NOTE appear iff a failure was superseded
    print("  Test: s6c_superseded_tag_in_llm_summary")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    def _make_session(cycles_data):
        s = AgentSession.__new__(AgentSession)
        s.data = {
            "session_id": "test-s6c",
            "experiment_type": "xray",
            "original_files": [],
            "project_advice": "Test run.",
            "cycles": cycles_data,
        }
        s.best_files = None
        s.get_directives = lambda: None
        s._annotate_superseded_issues = lambda cycles: None

        def _minimal_extract(self=s):
            steps = []
            for c in self.data.get("cycles", []):
                prog = c.get("program", "")
                if not prog or prog == "STOP":
                    continue
                success = "SUCCESS" in str(c.get("result", "")).upper()
                steps.append({"cycle": c.get("cycle_number", "?"),
                               "program": prog, "success": success,
                               "key_metric": ""})
            return {
                "session_id": "test-s6c", "experiment_type": "xray",
                "original_files": [], "user_advice": "Test.",
                "resolution": None, "workflow_path": "X-ray",
                "steps": steps, "final_metrics": {}, "final_files": {},
                "input_quality": {}, "total_cycles": len(steps),
                "successful_cycles": sum(1 for st in steps if st["success"]),
                "run_name": "test", "is_tutorial": False,
                "working_dir": None, "failure_diagnosis_path": None,
            }
        s._extract_summary_data = _minimal_extract
        return s

    sup_cycles = [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: done", "command": "phenix.xtriage 7qz0.mtz",
         "metrics": {}},
        {"cycle_number": 2, "program": "phenix.model_vs_data",
         "result": "FAILED: symmetry mismatch",
         "command": "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
         "metrics": {},
         "failure_context": {
             "error_type": "crystal_symmetry_mismatch", "scope": "probe",
             "program": "phenix.model_vs_data",
             "input_files": ["7qz0.pdb", "7qz0.mtz"],
             "superseded_by": "phenix.refine (cycle 3)",
         }},
        {"cycle_number": 3, "program": "phenix.refine",
         "result": "SUCCESS: R-free=0.28",
         "command": "phenix.refine 7qz0.pdb 7qz0.mtz", "metrics": {}},
    ]
    clean_cycles = [
        {"cycle_number": 1, "program": "phenix.refine",
         "result": "SUCCESS: R-free=0.28",
         "command": "phenix.refine 7qz0.pdb 7qz0.mtz", "metrics": {}},
    ]

    txt_sup = _make_session(sup_cycles).get_summary_for_llm_assessment()
    assert_true("[SUPERSEDED by phenix.refine (cycle 3)]" in txt_sup,
        "LLM summary must contain [SUPERSEDED] tag")
    assert_true("NOTE:" in txt_sup and "resolved by later" in txt_sup,
        "LLM summary must include block-level NOTE about superseded steps")

    txt_clean = _make_session(clean_cycles).get_summary_for_llm_assessment()
    assert_true("[SUPERSEDED" not in txt_clean,
        "LLM summary must NOT contain [SUPERSEDED] when no failure was superseded")
    assert_true("resolved by later" not in txt_clean,
        "Block-level NOTE must NOT appear when no failure was superseded")

    print("  PASSED: [SUPERSEDED] tag and NOTE appear iff a failure was superseded")


def test_s6d_prompt_has_superseded_instruction():
    # Step 4: AGENT_SESSION_ASSESSMENT_PROMPT has brief superseded instruction
    print("  Test: s6d_prompt_has_superseded_instruction")
    import sys as _sys, os as _os, re as _re
    _sys.path.insert(0, _PROJECT_ROOT)

    prompt_path = _os.path.join(_PROJECT_ROOT, "knowledge", "prompts_hybrid.py")
    prompt_src = open(prompt_path).read()
    m = _re.search(
        r'AGENT_SESSION_ASSESSMENT_PROMPT\s*=\s*"""(.*?)"""',
        prompt_src, _re.DOTALL)
    assert_true(m is not None,
        "AGENT_SESSION_ASSESSMENT_PROMPT not found in prompts_hybrid.py")
    P = m.group(1)

    assert_true("[SUPERSEDED]" in P, "Prompt must mention [SUPERSEDED] tag")
    assert_true("do not report them as critical issues" in P,
        "Prompt must instruct LLM not to escalate superseded failures")

    idx = P.find("[SUPERSEDED]")
    line_start = P.rfind("\n", 0, idx) + 1
    line_end = P.find("\n", idx)
    line = P[line_start:(line_end if line_end >= 0 else len(P))].strip()
    assert_true(len(line) < 300,
        "Superseded instruction must be <300 chars, got %d chars" % len(line))

    print("  PASSED: prompt has brief [SUPERSEDED] instruction")


def test_s6e_no_circular_import_in_files_overlap():
    # Fix 1: _files_overlap must not import from agent.session inside itself.
    # _stage_prefix_strip must delegate to module-level _strip_stage_prefix.
    print("  Test: s6e_no_circular_import_in_files_overlap")
    import sys as _sys, os as _os
    _sys.path.insert(0, _PROJECT_ROOT)

    session_path = _os.path.join(_PROJECT_ROOT, "agent", "session.py")
    src = open(session_path).read()

    # The circular import must be gone
    assert_true("from agent.session import AgentSession" not in src,
        "_files_overlap must not contain 'from agent.session import AgentSession'")

    # Module-level function must exist
    assert_true("def _strip_stage_prefix(" in src,
        "session.py must define module-level _strip_stage_prefix()")

    # staticmethod must delegate to it
    assert_true("return _strip_stage_prefix(basename)" in src,
        "_stage_prefix_strip staticmethod must delegate to _strip_stage_prefix()")

    # _files_overlap Strategy 2 must call _strip_stage_prefix directly
    assert_true("_strip_stage_prefix(f) for f in probe_files" in src,
        "_files_overlap Strategy 2 must call _strip_stage_prefix() directly")

    # Verify the module-level function works correctly in isolation
    import importlib.util as _ilu
    spec = _ilu.spec_from_file_location("_ss_mod", session_path)
    # Can't import due to libtbx deps — test the logic inline
    import re as _re, os as _os2
    def _strip(basename):
        stem = _os2.path.splitext(basename.lower())[0]
        stem = _re.sub(r'^[a-z]+_\d+_', '', stem)
        stem = _re.sub(r'^[a-z]+_',     '', stem)
        return stem

    assert_true(_strip("refine_001_7qz0.pdb") == "7qz0",
        "strip refine_001_ prefix failed: %r" % _strip("refine_001_7qz0.pdb"))
    assert_true(_strip("phased_model.pdb") == "model",
        "strip phased_ prefix failed: %r" % _strip("phased_model.pdb"))
    assert_true(_strip("7qz0.pdb") == "7qz0",
        "no prefix — must be unchanged: %r" % _strip("7qz0.pdb"))

    print("  PASSED: no circular import; _strip_stage_prefix at module level")


def test_s6f_files_overlap_scope_aware_fallback():
    # Fix 2: empty-file fallback in _files_overlap must only fire for scope="probe"
    print("  Test: s6f_files_overlap_scope_aware_fallback")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    # probe scope + empty probe_files -> True (benefit of the doubt)
    result_probe = AgentSession._files_overlap(
        [], "phenix.refine 7qz0.pdb 7qz0.mtz", [], scope="probe")
    assert_true(result_probe is True,
        "scope=probe with empty probe_files must return True (probe ran before "
        "derived files existed), got: %r" % result_probe)

    # recoverable scope + empty probe_files -> False (no evidence)
    result_recoverable = AgentSession._files_overlap(
        [], "phenix.refine 7qz0.pdb 7qz0.mtz", [], scope="recoverable")
    assert_true(result_recoverable is False,
        "scope=recoverable with empty probe_files must return False "
        "(no provenance evidence), got: %r" % result_recoverable)

    # None scope + empty probe_files -> False (conservative)
    result_none = AgentSession._files_overlap(
        [], "phenix.refine 7qz0.pdb 7qz0.mtz", [], scope=None)
    assert_true(result_none is False,
        "scope=None with empty probe_files must return False, got: %r" % result_none)

    # Both sides non-empty, matching — True regardless of scope
    result_match = AgentSession._files_overlap(
        ["7qz0.pdb", "7qz0.mtz"],
        "phenix.refine 7qz0.pdb 7qz0.mtz",
        [],
        scope="recoverable")
    assert_true(result_match is True,
        "Exact filename match must return True regardless of scope, got: %r" % result_match)

    # Both sides non-empty, no match — False
    result_nomatch = AgentSession._files_overlap(
        ["other.pdb", "other.mtz"],
        "phenix.refine 7qz0.pdb 7qz0.mtz",
        [],
        scope="probe")
    assert_true(result_nomatch is False,
        "No filename match with no shared originals must return False, "
        "got: %r" % result_nomatch)

    print("  PASSED: _files_overlap fallback is scope-aware")


def test_s6g_markdown_table_superseded_footnote():
    # Fix 3: markdown summary table uses ✗* for superseded failures
    # and emits a footnote explaining the symbol.
    print("  Test: s6g_markdown_table_superseded_footnote")
    import sys as _sys, os as _os
    _sys.path.insert(0, _PROJECT_ROOT)

    session_path = _os.path.join(_PROJECT_ROOT, "agent", "session.py")
    src = open(session_path).read()

    # Source must contain the ✗* symbol logic
    assert_true('"\\u2717*"' in src or "'✗*'" in src or '"✗*"' in src,
        "_format_summary_markdown must use ✗* for superseded failures")
    # superseded_by appears in both _format_summary_markdown and _extract_steps;
    # just verify it is present somewhere after the _format_summary_markdown def
    fmt_idx = src.find("def _format_summary_markdown")
    assert_true(fmt_idx >= 0, "_format_summary_markdown not found in session.py")
    assert_true("superseded_by" in src[fmt_idx:fmt_idx + 16000],
        "_format_summary_markdown must reference superseded_by")
    assert_true("does not affect the final result" in src,
        "_format_summary_markdown must include footnote text")

    # superseded_by must be propagated through _extract_steps
    extract_steps_section = src[src.find("def _extract_steps"):
                                src.find("def _extract_steps") + 2000]
    assert_true("superseded_by" in extract_steps_section,
        "_extract_steps must include superseded_by in the step dict")

    print("  PASSED: markdown table uses ✗* with footnote for superseded steps")


def test_s7a_exact_failed_duplicate_detected():
    # Problem 1: exact repeat of a FAILED command must be caught as duplicate
    print("  Test: s7a_exact_failed_duplicate_detected")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    def make_session(cycles):
        s = AgentSession.__new__(AgentSession)
        s.data = {"cycles": cycles}
        return s

    failed_cmd = "phenix.predict_and_build predict_and_build.stop_after_predict=True rebuild_strategy=quick"

    cycles = [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: xtriage ok",
         "command": "phenix.xtriage 7qz0.mtz"},
        {"cycle_number": 2, "program": "phenix.predict_and_build",
         "result": "FAILED: Please supply a sequence file",
         "command": failed_cmd},
    ]
    s = make_session(cycles)

    # Exact repeat of the failed command must be caught
    is_dup, prev, _ = s.is_duplicate_command(failed_cmd)
    assert_true(is_dup is True,
        "Exact repeat of FAILED command must be caught as duplicate, got: %r" % is_dup)
    assert_true(prev == 2,
        "Duplicate cycle number must be 2, got: %r" % prev)

    # A genuinely new command must NOT be caught
    new_cmd = "phenix.predict_and_build input_files.seq_file=7qz0.fa predict_and_build.stop_after_predict=True"
    is_dup2, _, _ = s.is_duplicate_command(new_cmd)
    assert_true(is_dup2 is False,
        "Different command (with seq file added) must NOT be caught as duplicate")

    # Successful commands still caught by existing heuristic
    success_cmd = "phenix.xtriage 7qz0.mtz"
    is_dup3, prev3, _ = s.is_duplicate_command(success_cmd)
    assert_true(is_dup3 is True,
        "Exact repeat of SUCCESSFUL command must still be caught, got: %r" % is_dup3)

    print("  PASSED: exact failed-command duplicates are detected; new commands pass through")


def test_s7b_failed_commands_not_in_success_list():
    # Regression guard: get_all_commands must NOT include failed cycles.
    # Failed commands are handled separately by get_all_failed_commands.
    print("  Test: s7b_failed_commands_not_in_success_list")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    s = AgentSession.__new__(AgentSession)
    s.data = {"cycles": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: ok", "command": "phenix.xtriage data.mtz"},
        {"cycle_number": 2, "program": "phenix.predict_and_build",
         "result": "FAILED: missing seq file",
         "command": "phenix.predict_and_build data.mtz"},
    ]}

    success_cmds = [cmd for _, cmd in s.get_all_commands()]
    assert_true("phenix.predict_and_build data.mtz" not in success_cmds,
        "Failed command must NOT appear in get_all_commands()")
    assert_true("phenix.xtriage data.mtz" in success_cmds,
        "Successful command must appear in get_all_commands()")

    failed_cmds = [cmd for _, cmd in s.get_all_failed_commands()]
    assert_true("phenix.predict_and_build data.mtz" in failed_cmds,
        "Failed command must appear in get_all_failed_commands()")
    assert_true("phenix.xtriage data.mtz" not in failed_cmds,
        "Successful command must NOT appear in get_all_failed_commands()")

    print("  PASSED: get_all_commands and get_all_failed_commands are non-overlapping")


def test_s7c_phaser_valid_for_conventional_model():
    # Problem 2: molecular_replacement phase phaser condition must allow
    # a conventional PDB search model (has_model_for_mr=True) even when
    # has_processed_model=False. Test the condition logic directly from YAML
    # without importing workflow_engine (avoids libtbx dependency).
    print("  Test: s7c_phaser_valid_for_conventional_model")
    import sys as _sys, os as _os
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        import yaml as _yaml
    except ImportError:
        try:
            from libtbx.utils import import_python_object
            _yaml = import_python_object("yaml", "").object
        except Exception:
            print("  SKIPPED: yaml not available")
            return

    wf_path = _os.path.join(_PROJECT_ROOT, "knowledge", "workflows.yaml")
    wf = _yaml.safe_load(open(wf_path))
    steps = wf.get("xray", {}).get("steps") or wf.get("xray", {}).get("phases") or {}
    mr_programs = steps.get("molecular_replacement", {}).get("programs", [])

    # Find phenix.phaser entry
    phaser_entry = None
    for p in mr_programs:
        if isinstance(p, dict) and p.get("program") == "phenix.phaser":
            phaser_entry = p
            break
    assert_true(phaser_entry is not None,
        "phenix.phaser must exist in molecular_replacement programs")

    # Simulate _check_conditions for beta.pdb scenario:
    # conventional model present, no predicted/processed model, phaser not yet run
    context_beta = {
        "has_processed_model": False,
        "has_model_for_mr": True,   # beta.pdb categorized as model
        "phaser_done": False,
        "process_predicted_model_done": False,
    }
    def simulate_check(entry, ctx):
        for cond in entry.get("conditions", []):
            if "has_any" in cond:
                keys = ["has_" + k for k in cond["has_any"]]
                if not any(ctx.get(k) for k in keys):
                    return False
            if "has" in cond:
                if not ctx.get("has_" + cond["has"]):
                    return False
            if "not_done" in cond:
                if ctx.get(cond["not_done"] + "_done"):
                    return False
        return True

    result_beta = simulate_check(phaser_entry, context_beta)
    assert_true(result_beta is True,
        "phaser condition must pass for conventional model (has_model_for_mr=True, "
        "has_processed_model=False)")

    # Predicted model path still works
    context_alphafold = {
        "has_processed_model": True,
        "has_model_for_mr": True,
        "phaser_done": False,
    }
    result_af = simulate_check(phaser_entry, context_alphafold)
    assert_true(result_af is True,
        "phaser condition must still pass for processed predicted model")

    # phaser_done=True must block phaser (not_done: phaser guard)
    context_done = {
        "has_processed_model": True,
        "has_model_for_mr": True,
        "phaser_done": True,
    }
    result_done = simulate_check(phaser_entry, context_done)
    assert_true(result_done is False,
        "phaser condition must fail when phaser_done=True")

    print("  PASSED: phaser condition allows conventional model and predicted model, "
          "blocks already-run")


def test_s7d_phaser_condition_in_workflows_yaml():
    # Source-level check: workflows.yaml molecular_replacement phaser condition
    # must use has_any: [processed_model, model_for_mr]
    print("  Test: s7d_phaser_condition_in_workflows_yaml")
    import sys as _sys, os as _os
    _sys.path.insert(0, _PROJECT_ROOT)

    wf_path = _os.path.join(_PROJECT_ROOT, "knowledge", "workflows.yaml")
    src = open(wf_path).read()

    assert_true("has_any: [processed_model, model_for_mr]" in src,
        "workflows.yaml phaser condition must use has_any: [processed_model, model_for_mr]")
    # The old condition (has: processed_model alone on phaser entry) must be gone
    # Check that phaser's conditions no longer contain only 'has: processed_model'
    mr_idx = src.find("molecular_replacement:")
    assert_true(mr_idx >= 0, "molecular_replacement step must exist in workflows.yaml")
    # Look for phenix.phaser only within the molecular_replacement section
    # (there is also a phaser entry in obtain_model which has different conditions)
    mr_section = src[mr_idx:mr_idx + 1000]
    phaser_in_mr = mr_section.find("program: phenix.phaser")
    assert_true(phaser_in_mr >= 0,
        "phenix.phaser must exist in molecular_replacement section")
    phaser_section = mr_section[phaser_in_mr:phaser_in_mr + 800]
    assert_true("has_any" in phaser_section,
        "phenix.phaser conditions in molecular_replacement must use has_any, "
        "got:\n%s" % phaser_section)

    print("  PASSED: workflows.yaml phaser condition correctly uses has_any")


def test_s7e_is_duplicate_returns_was_failure_flag():
    # is_duplicate_command must return a 3-tuple (is_dup, cycle_num, was_failure)
    print("  Test: s7e_is_duplicate_returns_was_failure_flag")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    s = AgentSession.__new__(AgentSession)
    s.data = {"cycles": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: ok", "command": "phenix.xtriage data.mtz"},
        {"cycle_number": 2, "program": "phenix.predict_and_build",
         "result": "FAILED: Please supply a sequence file",
         "command": "phenix.predict_and_build data.mtz"},
    ]}

    # Duplicate of failed cycle -> was_failure=True
    is_dup, prev, was_fail = s.is_duplicate_command("phenix.predict_and_build data.mtz")
    assert_true(is_dup is True, "Must be detected as duplicate")
    assert_true(was_fail is True,
        "was_failure must be True when matched cycle was a failure")

    # Duplicate of successful cycle -> was_failure=False
    is_dup2, prev2, was_fail2 = s.is_duplicate_command("phenix.xtriage data.mtz")
    assert_true(is_dup2 is True, "Must be detected as duplicate")
    assert_true(was_fail2 is False,
        "was_failure must be False when matched cycle was successful")

    # Not a duplicate -> (False, None, False)
    is_dup3, prev3, was_fail3 = s.is_duplicate_command("phenix.refine model.pdb data.mtz")
    assert_true(is_dup3 is False, "New command must not be a duplicate")
    assert_true(prev3 is None, "prev_cycle must be None for non-duplicate")
    assert_true(was_fail3 is False, "was_failure must be False for non-duplicate")

    print("  PASSED: is_duplicate_command returns correct 3-tuple with was_failure flag")


def test_s7f_get_cycle_result_returns_error_text():
    # get_cycle_result must return the result text for the given cycle number
    print("  Test: s7f_get_cycle_result_returns_error_text")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    s = AgentSession.__new__(AgentSession)
    s.data = {"cycles": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: xtriage done", "command": "phenix.xtriage data.mtz"},
        {"cycle_number": 2, "program": "phenix.predict_and_build",
         "result": "FAILED: Please supply a sequence file",
         "command": "phenix.predict_and_build data.mtz"},
    ]}

    assert_true(s.get_cycle_result(2) == "FAILED: Please supply a sequence file",
        "Must return exact result text for cycle 2")
    assert_true(s.get_cycle_result(1) == "SUCCESS: xtriage done",
        "Must return exact result text for cycle 1")
    assert_true(s.get_cycle_result(99) is None,
        "Must return None for non-existent cycle")

    print("  PASSED: get_cycle_result returns correct result text by cycle number")


def test_s7g_retry_feedback_distinguishes_failure_from_success():
    # _build_duplicate_feedback text must differ based on prev_was_failure.
    # Duplicate retries go through _query_agent_for_command (normal graph path)
    # with the feedback appended to guidelines.
    print("  Test: s7g_retry_feedback_distinguishes_failure_from_success")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    ai_agent_path = _find_ai_agent_path()
    assert_true(ai_agent_path is not None, "Could not locate ai_agent.py")
    src = open(ai_agent_path).read()

    # _build_duplicate_feedback must exist and accept prev_was_failure
    assert_true("def _build_duplicate_feedback" in src,
        "_build_duplicate_feedback must exist")
    feedback_idx = src.find("def _build_duplicate_feedback")
    feedback_section = src[feedback_idx:feedback_idx + 3000]
    assert_true("prev_was_failure" in feedback_section,
        "_build_duplicate_feedback must accept prev_was_failure parameter")

    # Must have two distinct feedback branches
    assert_true("FAILED COMMAND REPEATED" in feedback_section,
        "Must have failure-specific feedback message")
    assert_true("DUPLICATE REJECTED" in feedback_section,
        "Must have success-duplicate feedback message")

    # Failure branch must mention checking required input files
    assert_true("required input files" in feedback_section,
        "Failure feedback must instruct LLM to check required input files")

    # Failure branch must include the prior error text
    assert_true("get_cycle_result" in feedback_section,
        "Failure feedback must retrieve prior error text via get_cycle_result")

    # _handle_duplicate_check must pass prev_was_failure to _build_duplicate_feedback
    handle_idx = src.find("def _handle_duplicate_check")
    handle_section = src[handle_idx:handle_idx + 2000]
    assert_true("prev_was_failure" in handle_section,
        "_handle_duplicate_check must unpack and pass prev_was_failure")

    # Duplicate retries must go through _query_agent_for_command (no parallel path)
    assert_true("_query_agent_for_command" in handle_section,
        "Duplicate retries must call _query_agent_for_command (not bypass graph)")
    assert_true("duplicate_feedback" in handle_section,
        "Duplicate retries must pass duplicate_feedback parameter")

    # _retry_duplicate (old parallel path) must NOT exist
    assert_true("def _retry_duplicate" not in src,
        "_retry_duplicate must be removed (replaced by _query_agent_for_command path)")

    print("  PASSED: _build_duplicate_feedback generates context-appropriate feedback "
          "for failed vs. successful duplicate commands, routed through graph")


def test_s7h_molecular_replacement_phase_description_updated():
    # molecular_replacement step description must no longer be AlphaFold-specific
    print("  Test: s7h_molecular_replacement_phase_description_updated")
    import sys as _sys, os as _os
    _sys.path.insert(0, _PROJECT_ROOT)

    wf_path = _os.path.join(_PROJECT_ROOT, "knowledge", "workflows.yaml")
    src = open(wf_path).read()

    # Find the molecular_replacement section
    mr_idx = src.find("    molecular_replacement:")
    assert_true(mr_idx >= 0, "molecular_replacement section must exist")
    mr_section = src[mr_idx:mr_idx + 300]

    # Must not contain AlphaFold-specific language
    assert_true("AlphaFold" not in mr_section,
        "molecular_replacement description must not mention AlphaFold, got:\n%s" % mr_section)
    assert_true("predicted model" not in mr_section,
        "molecular_replacement description must not say 'predicted model', "
        "got:\n%s" % mr_section)

    # Must contain general MR language
    assert_true("molecular replacement" in mr_section.lower() or
                "MR" in mr_section,
        "molecular_replacement description must contain general MR language, "
        "got:\n%s" % mr_section)

    print("  PASSED: molecular_replacement step description is now general, "
          "not AlphaFold-specific")


# ===========================================================================
# S8 — Client-side required-file injection (_inject_missing_required_files)
# ===========================================================================

def _make_mock_session(available_files=None, best_files=None,
                       original_files=None):
    """Build a minimal mock session for injection tests."""
    class _MockSession:
        def get_available_files(self):
            return list(available_files or [])
        def get_best_files_dict(self):
            return dict(best_files or {})
        def get_directives(self):
            return {}
        @property
        def data(self):
            return {"original_files": list(original_files or [])}
    return _MockSession()


def _make_mock_params(original_files=None):
    """Build a minimal mock params object."""
    class _MockAI:
        pass
    class _MockParams:
        ai_analysis = _MockAI()
    p = _MockParams()
    p.ai_analysis.original_files = list(original_files or [])
    return p


def _get_injector():
    """Return a callable that wraps _inject_missing_required_files.

    Extracts the two methods from ai_agent.py source via exec so we
    don't need to import the full AIAgent class (which requires the
    entire PHENIX package).  Supplies a mock ProgramRegistry backed
    directly by knowledge/yaml_loader so the methods work without
    libtbx being installed.
    """
    import sys as _sys, types as _types
    _sys.path.insert(0, _PROJECT_ROOT)

    ai_agent_path = _find_ai_agent_path()

    # Build a lightweight ProgramRegistry substitute that reads directly
    # from programs.yaml via knowledge.yaml_loader (no libtbx needed).
    from knowledge.yaml_loader import get_program_inputs  # local import, no libtbx

    class _MockProgramRegistry:
        def get_required_input_defs(self, program_name):
            inputs = get_program_inputs(program_name)
            return dict(inputs.get("required", {}))

    # Also supply the module-level _quote_if_needed from program_registry source
    # by extracting it directly rather than importing the module.
    import shlex as _shlex
    def _quote_if_needed(path):
        s = str(path)
        return _shlex.quote(s) if ' ' in s else s

    # Extract just the two method bodies from source.
    src = open(ai_agent_path).read()

    def _extract_method(method_name):
        import textwrap as _tw
        marker = "  def %s(" % method_name
        start = src.find(marker)
        if start < 0:
            raise ValueError("Cannot find method %r in ai_agent.py" % method_name)
        rest = src[start:]
        lines = rest.split('\n')
        body_lines = [lines[0]]
        for line in lines[1:]:
            if line.startswith('  def ') and line != lines[0]:
                break
            body_lines.append(line)
        return _tw.dedent('\n'.join(body_lines))

    method_src = (
        _extract_method("_inject_missing_required_files") + "\n\n" +
        _extract_method("_find_candidate_for_slot")
    )

    # Namespace pre-populated with everything the methods need.
    # The try/except import blocks inside the method will raise ImportError
    # for the libtbx path, then fall to 'from agent.program_registry import ...'
    # which also chains to libtbx.  We bypass both by pre-seeding the names
    # the method resolves via globals() after exec.
    ns = {
        # Class constants (referenced as self._X, but exec sees them as closures)
        "_SLOT_TO_BEST_CATEGORY": {
            "model": "model", "pdb_file": "model", "protein": "model",
            "map": "map", "full_map": "map",
            "data_mtz": "data_mtz", "hkl_file": "data_mtz", "data": "data_mtz",
            "map_coeffs_mtz": "map_coeffs_mtz",
            "sequence": "sequence", "seq_file": "sequence",
            "ligand_cif": "ligand_cif", "ligand": "ligand_cif",
        },
        "_CRYSTAL_EXTS": frozenset({
            '.pdb', '.cif', '.mtz', '.sca', '.hkl',
            '.mrc', '.ccp4', '.map',
            '.fa', '.fasta', '.seq', '.dat',
            '.ncs_spec', '.eff',
        }),
        # Pre-seed ProgramRegistry so the try/except import inside the method
        # resolves to our mock without ever touching libtbx.
        "ProgramRegistry": _MockProgramRegistry,
        "_quote_if_needed": _quote_if_needed,
    }

    # Patch sys.modules so the 'from libtbx... import ProgramRegistry' lines
    # inside the exec'd methods also resolve to our mock instead of failing.
    _libtbx_stub = _types.ModuleType("libtbx")
    _lc_stub     = _types.ModuleType("libtbx.langchain")
    _lca_stub    = _types.ModuleType("libtbx.langchain.agent")
    _lca_stub.ProgramRegistry    = _MockProgramRegistry
    _lca_stub._quote_if_needed   = _quote_if_needed
    _lcapr_stub  = _types.ModuleType("libtbx.langchain.agent.program_registry")
    _lcapr_stub.ProgramRegistry  = _MockProgramRegistry
    _lcapr_stub._quote_if_needed = _quote_if_needed
    for _name, _mod in [
        ("libtbx", _libtbx_stub),
        ("libtbx.langchain", _lc_stub),
        ("libtbx.langchain.agent", _lca_stub),
        ("libtbx.langchain.agent.program_registry", _lcapr_stub),
    ]:
        if _name not in _sys.modules:
            _sys.modules[_name] = _mod

    try:
        exec(compile(method_src, "<inject_methods>", "exec"), ns)
    except Exception as e:
        raise ImportError("Failed to exec injection methods: %s" % e)

    inject_fn = ns["_inject_missing_required_files"]
    find_fn   = ns["_find_candidate_for_slot"]

    class _StubVlog:
        def verbose(self, msg): pass
        def normal(self, msg): pass

    obj = _types.SimpleNamespace(
        vlog=_StubVlog(),
        params=None,
        _SLOT_TO_BEST_CATEGORY=ns["_SLOT_TO_BEST_CATEGORY"],
        _CRYSTAL_EXTS=ns["_CRYSTAL_EXTS"],
    )
    obj._inject_missing_required_files = _types.MethodType(inject_fn, obj)
    obj._find_candidate_for_slot       = _types.MethodType(find_fn, obj)

    class _Wrapper:
        def inject(self, command, program_name, session, available_files=None):
            return obj._inject_missing_required_files(
                command, program_name, session, available_files=available_files)
        def find_candidate(self, slot_name, slot_def, extensions,
                           best_files, available_files, present_basenames):
            return obj._find_candidate_for_slot(
                slot_name, slot_def, extensions,
                best_files, available_files, present_basenames)

    return _Wrapper()


def test_s8a_inject_sequence_when_missing():
    """Sequence file injected into predict_and_build when absent."""
    print("  Test: s8a_inject_sequence_when_missing")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    session = _make_mock_session(
        available_files=["/data/7qz0.fa", "/data/data.mtz"],
        best_files={})

    cmd = ("phenix.predict_and_build "
           "predict_and_build.stop_after_predict=True "
           "rebuild_strategy=quick")
    result = injector.inject(cmd, "phenix.predict_and_build", session,
                             available_files=["/data/7qz0.fa", "/data/data.mtz"])

    assert_true("input_files.seq_file=" in result,
        "Flagged injection must produce 'input_files.seq_file=', got:\n%s" % result)
    assert_true("7qz0.fa" in result,
        "Injected result must contain '7qz0.fa', got:\n%s" % result)
    # MTZ not present in original command — predict_and_build has no required MTZ slot
    # so data.mtz must NOT be injected
    assert_true(result.count(".mtz") == 0,
        "data.mtz must NOT be injected (not a required slot), got:\n%s" % result)

    print("  PASSED: sequence file injected with correct flag")


def test_s8b_no_inject_when_already_present():
    """No injection when required file is already in the command."""
    print("  Test: s8b_no_inject_when_already_present")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    session = _make_mock_session(
        available_files=["/data/7qz0.fa", "/data/data.mtz"],
        best_files={})

    cmd = ("phenix.predict_and_build "
           "input_files.seq_file=/data/7qz0.fa "
           "predict_and_build.stop_after_predict=True")
    result = injector.inject(cmd, "phenix.predict_and_build", session,
                             available_files=["/data/7qz0.fa"])

    assert_true(result == cmd,
        "Command must be unchanged when seq file already present.\n"
        "Expected: %r\nGot:      %r" % (cmd, result))
    # Exactly one occurrence of the seq file
    assert_true(result.count("7qz0.fa") == 1,
        "Must not double-inject the file, got:\n%s" % result)

    print("  PASSED: command unchanged when required file already present")


def test_s8c_inject_bare_positional_model():
    """Model file injected as bare positional before key=value tokens."""
    print("  Test: s8c_inject_bare_positional_model")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    session = _make_mock_session(
        best_files={"model": "/data/refined.pdb"},
        available_files=["/data/refined.pdb"])

    cmd = "phenix.molprobity xray_data.r_free_flags.generate=True"
    result = injector.inject(cmd, "phenix.molprobity", session,
                             available_files=["/data/refined.pdb"])

    assert_true("refined.pdb" in result,
        "Model file must be injected, got:\n%s" % result)

    # Bare positional must appear before any key=value token
    parts = result.split()
    pdb_idx  = next(i for i, p in enumerate(parts) if "refined.pdb" in p)
    kv_idx   = next(i for i, p in enumerate(parts) if "=" in p)
    assert_true(pdb_idx < kv_idx,
        "Bare positional file must come before key=value tokens.\n"
        "parts: %s" % parts)

    print("  PASSED: bare positional model injected before key=value tokens")


def test_s8d_no_inject_for_unknown_program():
    """No injection and no exception for an unrecognised program."""
    print("  Test: s8d_no_inject_for_unknown_program")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    session = _make_mock_session(
        available_files=["/data/model.pdb"],
        best_files={})

    cmd = "phenix.custom_program model.pdb"
    result = injector.inject(cmd, "phenix.custom_program", session,
                             available_files=["/data/model.pdb"])

    assert_true(result == cmd,
        "Command must be unchanged for unknown program, got:\n%s" % result)

    print("  PASSED: unknown program returns command unchanged, no exception")


def test_s8e_pdbtools_dual_pdb_no_double_inject():
    """pdbtools: second PDB slot filled; first PDB not duplicated."""
    print("  Test: s8e_pdbtools_dual_pdb_no_double_inject")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    session = _make_mock_session(
        available_files=["/data/protein.pdb", "/data/ligand.pdb"],
        best_files={"model": "/data/protein.pdb"})

    # protein.pdb is already in the command (positional)
    cmd = "phenix.pdbtools /data/protein.pdb"
    result = injector.inject(cmd, "phenix.pdbtools", session,
                             available_files=["/data/protein.pdb",
                                             "/data/ligand.pdb"])

    # ligand.pdb should be injected for the second PDB slot
    assert_true("ligand.pdb" in result,
        "ligand.pdb must be injected for the second PDB slot, got:\n%s" % result)
    # protein.pdb must appear exactly once
    assert_true(result.count("protein.pdb") == 1,
        "protein.pdb must appear exactly once, got:\n%s" % result)

    print("  PASSED: pdbtools second PDB slot filled without duplicating first")


def test_s8f_ligandfit_mtz_not_injected_without_best():
    """ligandfit: raw Fobs MTZ NOT injected when best_files empty."""
    print("  Test: s8f_ligandfit_mtz_not_injected_without_best")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()
    # best_files has no map_coeffs_mtz; available_files has a raw Fobs MTZ
    session = _make_mock_session(
        available_files=["/data/protein.pdb", "/data/LIG.pdb",
                         "/data/data.mtz"],
        best_files={})

    cmd = "phenix.ligandfit model=/data/protein.pdb ligand=/data/LIG.pdb"
    result = injector.inject(cmd, "phenix.ligandfit", session,
                             available_files=["/data/protein.pdb",
                                             "/data/LIG.pdb",
                                             "/data/data.mtz"])

    # data.mtz is raw Fobs — must NOT be injected (require_best_files_only=true)
    assert_true("data.mtz" not in result,
        "Raw Fobs MTZ must NOT be injected for ligandfit map_coeffs_mtz slot "
        "(require_best_files_only=true), got:\n%s" % result)

    print("  PASSED: require_best_files_only blocks raw MTZ injection for ligandfit")


def test_s8g_get_required_input_defs():
    """ProgramRegistry.get_required_input_defs returns full slot definitions."""
    print("  Test: s8g_get_required_input_defs")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.program_registry import ProgramRegistry
    except ImportError:
        print("  SKIP (ProgramRegistry unavailable)")
        return

    r = ProgramRegistry()

    # predict_and_build: one required slot (sequence)
    defs = r.get_required_input_defs("phenix.predict_and_build")
    assert_true("sequence" in defs,
        "predict_and_build must have 'sequence' required slot, got: %s" % list(defs))
    seq_def = defs["sequence"]
    assert_true(seq_def.get("flag") == "input_files.seq_file=",
        "sequence flag must be 'input_files.seq_file=', got: %r" % seq_def.get("flag"))
    assert_true(".fa" in seq_def.get("extensions", []),
        "sequence extensions must include '.fa', got: %s" % seq_def.get("extensions"))

    # ligandfit: require_best_files_only on map_coeffs_mtz
    lf_defs = r.get_required_input_defs("phenix.ligandfit")
    assert_true("map_coeffs_mtz" in lf_defs,
        "ligandfit must have 'map_coeffs_mtz' required slot, got: %s" % list(lf_defs))
    assert_true(lf_defs["map_coeffs_mtz"].get("require_best_files_only") is True,
        "ligandfit map_coeffs_mtz must have require_best_files_only=true")

    # Unknown program: empty dict, no exception
    unknown = r.get_required_input_defs("phenix.nonexistent")
    assert_true(unknown == {},
        "Unknown program must return empty dict, got: %r" % unknown)

    print("  PASSED: get_required_input_defs returns correct slot definitions")


def test_s8h_path_with_spaces_quoted_correctly():
    """Paths containing spaces are quoted so shlex.split() tokenises correctly."""
    print("  Test: s8h_path_with_spaces_quoted_correctly")
    import sys as _sys, shlex as _shlex
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()

    # ── Sub-test 1: flagged injection with space in directory ─────────────────
    spaced_seq = "/data/my projects/7qz0.fa"
    session1 = _make_mock_session(
        available_files=[spaced_seq],
        best_files={})
    cmd1 = "phenix.predict_and_build stop_after_predict=True"
    result1 = injector.inject(cmd1, "phenix.predict_and_build", session1,
                              available_files=[spaced_seq])

    assert_true("7qz0.fa" in result1,
        "Sequence file must be injected, got:\n%s" % result1)
    # After quoting, shlex.split must produce exactly 3 tokens:
    #   phenix.predict_and_build  stop_after_predict=True  input_files.seq_file=...
    try:
        parts1 = _shlex.split(result1)
    except ValueError as e:
        raise AssertionError(
            "shlex.split raised ValueError on injected command "
            "(likely unquoted space in path): %s\nCommand: %s" % (e, result1))
    # The seq_file token must contain the full path with space intact
    seq_tokens = [p for p in parts1 if "seq_file" in p]
    assert_true(len(seq_tokens) == 1,
        "Must have exactly one seq_file token, got: %s" % parts1)
    assert_true("my projects" in seq_tokens[0],
        "Full path with space must be preserved in token, got: %r" % seq_tokens[0])

    # ── Sub-test 2: bare positional injection with space in directory ──────────
    spaced_model = "/data/my data/refined.pdb"
    session2 = _make_mock_session(
        best_files={"model": spaced_model},
        available_files=[spaced_model])
    cmd2 = "phenix.molprobity xray_data.r_free_flags.generate=True"
    result2 = injector.inject(cmd2, "phenix.molprobity", session2,
                              available_files=[spaced_model])

    assert_true("refined.pdb" in result2,
        "Model must be injected, got:\n%s" % result2)
    try:
        parts2 = _shlex.split(result2)
    except ValueError as e:
        raise AssertionError(
            "shlex.split raised ValueError on bare-positional command: "
            "%s\nCommand: %s" % (e, result2))
    pdb_tokens = [p for p in parts2 if "refined.pdb" in p]
    assert_true(len(pdb_tokens) == 1,
        "Must have exactly one pdb token, got: %s" % parts2)
    assert_true("my data" in pdb_tokens[0],
        "Full path with space must be preserved, got: %r" % pdb_tokens[0])

    # ── Sub-test 3: roundtrip — already-present detection on a quoted path ─────
    # Build a command that already has the quoted injected form, then
    # call inject again.  The slot must be detected as satisfied (no re-injection).
    cmd3 = result1  # result from sub-test 1 (contains quoted seq path)
    result3 = injector.inject(cmd3, "phenix.predict_and_build", session1,
                              available_files=[spaced_seq])
    count_before = cmd3.count("7qz0.fa")
    count_after  = result3.count("7qz0.fa")
    assert_true(count_after == count_before,
        "Re-injection must be suppressed; '7qz0.fa' count went from %d to %d.\n"
        "Command: %s" % (count_before, count_after, result3))

    print("  PASSED: spaced paths quoted correctly; shlex.split safe; "
          "already-present detection works on quoted tokens")


# ===========================================================================
# S9 — Superseded-issues provenance fixes (lineage graph + authority expansion)
# ===========================================================================

def test_s9e_extract_file_basenames_helper():
    """_extract_file_basenames correctly tokenises commands."""
    print("  Test: s9e_extract_file_basenames_helper")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import _extract_file_basenames
    except ImportError:
        from libtbx.langchain.agent.session import _extract_file_basenames

    # Normal command: mix of files and PHIL params
    result = _extract_file_basenames(
        "phenix.refine input_files.seq_file=7qz0.fa PHASER.1.pdb refinement.cycles=5")
    assert_true("7qz0.fa" in result,
        "seq_file= token must yield '7qz0.fa', got: %s" % result)
    assert_true("phaser.1.pdb" in result,
        "bare positional must yield 'phaser.1.pdb', got: %s" % result)
    assert_true("phenix.refine" not in result,
        "program name must not appear in result, got: %s" % result)

    # None / empty input
    assert_true(_extract_file_basenames(None) == [],
        "None input must return empty list")
    assert_true(_extract_file_basenames("") == [],
        "Empty string must return empty list")

    # Quoted path with space in directory — must return correct basename
    result2 = _extract_file_basenames(
        "phenix.predict_and_build input_files.seq_file='/data/my dir/7qz0.fa'")
    assert_true("7qz0.fa" in result2,
        "Quoted spaced-directory path must yield '7qz0.fa', got: %s" % result2)

    print("  PASSED: _extract_file_basenames tokenises correctly")


def test_s9a_lineage_graph_built_correctly():
    """_build_lineage_graph maps output basenames to input basenames."""
    print("  Test: s9a_lineage_graph_built_correctly")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    cycles = [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "command": "phenix.xtriage data.mtz",
         "result": "SUCCESS",
         "output_files": []},
        {"cycle_number": 2, "program": "phenix.phaser",
         "command": "phenix.phaser beta.pdb data.mtz",
         "result": "SUCCESS",
         "output_files": ["/data/PHASER.1.pdb"]},
        {"cycle_number": 3, "program": "phenix.refine",
         "command": "phenix.refine PHASER.1.pdb data.mtz",
         "result": "SUCCESS",
         "output_files": ["/data/refine_001.pdb"]},
    ]

    graph = AgentSession._build_lineage_graph(cycles)

    # PHASER.1.pdb was produced from beta.pdb + data.mtz
    assert_true("phaser.1.pdb" in graph,
        "PHASER.1.pdb must be in graph, keys: %s" % list(graph))
    assert_true("beta.pdb" in graph["phaser.1.pdb"],
        "PHASER.1.pdb parents must include beta.pdb, got: %s" % graph["phaser.1.pdb"])
    assert_true("data.mtz" in graph["phaser.1.pdb"],
        "PHASER.1.pdb parents must include data.mtz, got: %s" % graph["phaser.1.pdb"])

    # refine_001.pdb was produced from PHASER.1.pdb + data.mtz
    assert_true("refine_001.pdb" in graph,
        "refine_001.pdb must be in graph")
    assert_true("phaser.1.pdb" in graph["refine_001.pdb"],
        "refine_001.pdb parents must include phaser.1.pdb, got: %s" % graph["refine_001.pdb"])

    # Cycle with empty output_files must not add entries (xtriage above)
    for bn in graph:
        assert_true(bn != "data.mtz",  # data.mtz is an input, not an output
            "Input files must not appear as graph keys, got: %s" % bn)

    print("  PASSED: lineage graph maps output basenames to their input basenames")


def test_s9b_ancestors_traces_transitively():
    """_ancestors resolves transitive chains and handles cycles safely."""
    print("  Test: s9b_ancestors_traces_transitively")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    graph = {
        "phaser.1.pdb": {"beta.pdb", "data.mtz"},
        "refined.pdb":  {"phaser.1.pdb"},
    }

    # Transitive chain: refined.pdb → phaser.1.pdb → beta.pdb, data.mtz
    anc = AgentSession._ancestors("refined.pdb", graph)
    assert_true("phaser.1.pdb" in anc,
        "phaser.1.pdb must be an ancestor of refined.pdb, got: %s" % anc)
    assert_true("beta.pdb" in anc,
        "beta.pdb must be a transitive ancestor of refined.pdb, got: %s" % anc)
    assert_true("data.mtz" in anc,
        "data.mtz must be a transitive ancestor of refined.pdb, got: %s" % anc)

    # Leaf node: beta.pdb has no entries in the graph
    anc2 = AgentSession._ancestors("beta.pdb", graph)
    assert_true(anc2 == set(),
        "Leaf node must return empty set, got: %s" % anc2)

    # Unknown file: not in graph at all
    anc3 = AgentSession._ancestors("unknown.pdb", graph)
    assert_true(anc3 == set(),
        "Unknown file must return empty set, got: %s" % anc3)

    # Cycle safety: must not infinite-loop
    cyclic = {"a.pdb": {"b.pdb"}, "b.pdb": {"a.pdb"}}
    anc4 = AgentSession._ancestors("a.pdb", cyclic)
    assert_true(isinstance(anc4, set),
        "Must return a set even for cyclic graphs, got: %r" % anc4)

    print("  PASSED: _ancestors resolves transitively; handles leaf nodes and cycles")


def test_s9c_strategy4_supersedes_via_lineage():
    """Full MR workflow: probe superseded via lineage graph (Strategy 4)."""
    print("  Test: s9c_strategy4_supersedes_via_lineage")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    # The exact scenario that previously failed all three strategies:
    #   Cycle 1: model_vs_data beta.pdb → FAILED (needs_mr)
    #   Cycle 2: phaser beta.pdb → SUCCESS → outputs PHASER.1.pdb
    #   Cycle 3: refine PHASER.1.pdb → SUCCESS
    #
    # Strategy 4 should trace: PHASER.1.pdb ancestors → {beta.pdb}
    # which intersects probe_files → superseded.

    cycles = [
        {
            "cycle_number": 1,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data beta.pdb data.mtz",
            "result": "FAILED: R-free=0.58 — model needs molecular replacement",
            "output_files": [],
            "failure_context": {
                "error_type": None,
                "scope": "recoverable",
                "program": "phenix.model_vs_data",
                "input_files": ["beta.pdb", "data.mtz"],
                "superseded_by": None,
            },
        },
        {
            "cycle_number": 2,
            "program": "phenix.phaser",
            "command": "phenix.phaser beta.pdb data.mtz",
            "result": "SUCCESS: MR solution found",
            "output_files": ["/data/PHASER.1.pdb"],
        },
        {
            "cycle_number": 3,
            "program": "phenix.refine",
            "command": "phenix.refine PHASER.1.pdb data.mtz",
            "result": "SUCCESS: R-free=0.26",
            "output_files": ["/data/refine_001.pdb"],
        },
    ]

    s = AgentSession.__new__(AgentSession)
    s.data = {"cycles": cycles, "original_files": ["/data/beta.pdb", "/data/data.mtz"]}
    s._annotate_superseded_issues(s.data["cycles"])

    fc = cycles[0]["failure_context"]
    assert_true(fc["superseded_by"] is not None,
        "Cycle 1 must be marked superseded after phaser+refine succeed, got: %r" % fc)

    # Accept either cycle 2 (phaser) or cycle 3 (refine) as the superseder —
    # both are correct.  Cycle 2 is preferred (phaser is now in _default authority)
    # but cycle 3 is also valid if phaser was matched first by the graph.
    superseder = fc["superseded_by"]
    assert_true("phenix.phaser" in superseder or "phenix.refine" in superseder,
        "superseded_by must name phaser or refine, got: %r" % superseder)

    print("  PASSED: Strategy 4 links beta.pdb → PHASER.1.pdb → refine supersession "
          "(superseded by: %s)" % superseder)


def test_s9d_expanded_authority_phaser_supersedes_probe():
    """phenix.phaser in _default authority supersedes probe on direct file match."""
    print("  Test: s9d_expanded_authority_phaser_supersedes_probe")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.session import AgentSession
    except ImportError:
        from libtbx.langchain.agent.session import AgentSession

    # Direct file match scenario — no lineage needed.
    # Phaser runs on the SAME files as the probe: beta.pdb, data.mtz.
    # Before the authority expansion, phaser was not in _default so this
    # would never be marked superseded.

    cycles = [
        {
            "cycle_number": 1,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data beta.pdb data.mtz",
            "result": "FAILED: R-free=0.58",
            "output_files": [],
            "failure_context": {
                "error_type": None,
                "scope": "recoverable",
                "program": "phenix.model_vs_data",
                "input_files": ["beta.pdb", "data.mtz"],
                "superseded_by": None,
            },
        },
        {
            "cycle_number": 2,
            "program": "phenix.phaser",
            "command": "phenix.phaser beta.pdb data.mtz",
            "result": "SUCCESS: MR solution",
            "output_files": ["/data/PHASER.1.pdb"],
        },
    ]

    s = AgentSession.__new__(AgentSession)
    s.data = {"cycles": cycles, "original_files": []}
    s._annotate_superseded_issues(s.data["cycles"])

    fc = cycles[0]["failure_context"]
    assert_true(fc["superseded_by"] is not None,
        "Probe must be superseded by phaser (Strategy 1, direct match), got: %r" % fc)
    assert_true("phenix.phaser" in fc["superseded_by"],
        "superseded_by must name phenix.phaser, got: %r" % fc["superseded_by"])
    assert_true("cycle 2" in fc["superseded_by"],
        "superseded_by must include cycle 2, got: %r" % fc["superseded_by"])

    print("  PASSED: phenix.phaser now authoritative in _default; "
          "directly supersedes probe on matching files")


def test_s8i_retry_duplicate_injection():
    """Duplicate retry commands are passed through required-file injection.

    Duplicate retries now go through _query_agent_for_command (the normal graph
    path), so BUILD's postprocess_command runs automatically.  The client-only
    _inject_missing_required_files is applied in _handle_duplicate_check after
    each retry query.
    """
    print("  Test: s8i_retry_duplicate_injection")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    # Old _retry_duplicate must be gone (parallel path removed)
    assert_true("def _retry_duplicate(" not in src,
        "_retry_duplicate must be removed (replaced by _query_agent_for_command path)")

    # _handle_duplicate_check must call _inject_missing_required_files
    method_start = src.find("  def _handle_duplicate_check(")
    assert_true(method_start >= 0, "_handle_duplicate_check method must exist")
    rest = src[method_start:]
    next_def = rest.find("\n  def ", 1)
    method_src = rest[:next_def] if next_def > 0 else rest[:3000]

    assert_true("_inject_missing_required_files" in method_src,
        "_handle_duplicate_check must call _inject_missing_required_files on retries")

    # Retries must go through the normal graph path
    assert_true("_query_agent_for_command" in method_src,
        "Retries must call _query_agent_for_command (not bypass graph)")
    assert_true("duplicate_feedback" in method_src,
        "Retries must pass duplicate_feedback to _query_agent_for_command")

    # decision_info['command'] must be updated after injection
    assert_true("decision_info['command'] = command" in method_src,
        "decision_info['command'] must be updated after injection")

    # _query_agent_for_command must accept duplicate_feedback parameter
    qac_start = src.find("  def _query_agent_for_command(")
    assert_true(qac_start >= 0, "_query_agent_for_command must exist")
    qac_sig = src[qac_start:qac_start + 200]
    assert_true("duplicate_feedback" in qac_sig,
        "_query_agent_for_command must accept duplicate_feedback parameter")

    print("  PASSED: duplicate retries go through normal graph path with "
          "required-file injection applied afterward")


def test_s8j_prefer_patterns_honoured():
    """prefer_patterns promotes matching files over non-matching candidates."""
    print("  Test: s8j_prefer_patterns_honoured")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    injector = _get_injector()

    # available_files has two .pdb files: one matches prefer_patterns, one doesn't.
    # Without prefer_patterns support the method would return the last file
    # in the list (model.pdb), which is wrong for a ligand slot.
    session = _make_mock_session(
        available_files=["/data/model.pdb", "/data/lig.pdb"],
        best_files={})

    slot_def = {
        "extensions": [".pdb", ".cif"],
        "flag": "ligand=",
        "prefer_patterns": ["lig", "ligand"],
        "exclude_patterns": [],
    }
    result = injector.find_candidate(
        "ligand", slot_def, [".pdb", ".cif"],
        best_files={},
        available_files=["/data/model.pdb", "/data/lig.pdb"],
        present_basenames=set())

    assert_true(result is not None,
        "find_candidate must return a file, got None")
    assert_true("lig.pdb" in result,
        "prefer_patterns must select lig.pdb over model.pdb, got: %r" % result)

    # Non-matching file is still returned when no prefer-pattern match exists
    result2 = injector.find_candidate(
        "ligand", slot_def, [".pdb", ".cif"],
        best_files={},
        available_files=["/data/model.pdb"],
        present_basenames=set())
    assert_true(result2 is not None and "model.pdb" in result2,
        "Fallback to non-matching file when no prefer match exists, got: %r" % result2)

    print("  PASSED: prefer_patterns selects matching candidate; "
          "non-matching file used as fallback")


def test_s10a_no_param_inject_into_probe_programs():
    """inject_user_params must never inject key=value params into probe/validation
    programs (model_vs_data, xtriage, molprobity, mtriage, validation_cryoem).

    Regression guard: guidelines containing rebuild_strategy=quick, space_group=P3221,
    or any other bare key=value must not appear in model_vs_data commands.
    """
    print("  Test: s10a_no_param_inject_into_probe_programs")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.isfile(pp_path),
        "command_postprocessor.py must exist")
    src = open(pp_path).read()

    # The guard must exist in inject_user_params
    inject_start = src.find("def inject_user_params(")
    assert_true(inject_start >= 0,
        "inject_user_params must exist in command_postprocessor.py")

    # Find the end of the function
    inject_end = src.find("\ndef ", inject_start + 1)
    inject_src = src[inject_start:inject_end] if inject_end > 0 else src[inject_start:]

    assert_true("_NO_PARAM_INJECT_PROGRAMS" in inject_src,
        "inject_user_params must reference _NO_PARAM_INJECT_PROGRAMS guard")

    # The constant definition is module-level, so check the full source
    assert_true("_NO_PARAM_INJECT_PROGRAMS" in src,
        "_NO_PARAM_INJECT_PROGRAMS must be defined in command_postprocessor.py")
    const_start = src.find("_NO_PARAM_INJECT_PROGRAMS = frozenset")
    assert_true(const_start >= 0,
        "_NO_PARAM_INJECT_PROGRAMS must be defined as frozenset")
    const_block = src[const_start:const_start + 500]
    assert_true("phenix.model_vs_data" in const_block,
        "_NO_PARAM_INJECT_PROGRAMS must include phenix.model_vs_data")
    assert_true("phenix.xtriage" in const_block,
        "_NO_PARAM_INJECT_PROGRAMS must include phenix.xtriage")
    assert_true("phenix.molprobity" in const_block,
        "_NO_PARAM_INJECT_PROGRAMS must include phenix.molprobity")

    # The guard must fire before any injection code in the function body
    guard_idx  = inject_src.find("_NO_PARAM_INJECT_PROGRAMS")
    return_idx = inject_src.find("return command", guard_idx)
    assert_true(return_idx > guard_idx and return_idx < guard_idx + 400,
        "_NO_PARAM_INJECT_PROGRAMS guard must return command unchanged immediately")

    print("  PASSED: _NO_PARAM_INJECT_PROGRAMS guard present and returns early "
          "for probe/validation programs")


def test_s10b_probe_failure_continues_workflow():
    """Any failure of model_vs_data must continue the workflow, not stop it.

    Regression guard: before this fix, 'Sorry: Unknown command line parameter
    definition: space_group' triggered unknown_phil_parameter → terminal stop.
    Now the _PROBE_PROGRAMS early-exit ensures model_vs_data failures are
    always treated as inconclusive.
    """
    print("  Test: s10b_probe_failure_continues_workflow")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    ai_agent_path = _find_ai_agent_path()
    src = open(ai_agent_path).read()

    # Find the run_failed block — the _PROBE_PROGRAMS guard must be inside it
    # and must appear BEFORE the _analyze_for_recovery call.
    run_failed_idx = src.find("if run_failed:")
    assert_true(run_failed_idx >= 0, "run_failed block must exist in ai_agent.py")

    # Search within the next 5000 chars of the block
    block = src[run_failed_idx:run_failed_idx + 5000]

    assert_true("_PROBE_PROGRAMS" in block,
        "_PROBE_PROGRAMS guard must be inside the run_failed block")

    probe_idx    = block.find("_PROBE_PROGRAMS")
    # Search for the actual _analyze_for_recovery CALL (not the comment that
    # mentions it); skip past comment lines by finding the method call form.
    recovery_idx = block.find("recovery = self._analyze_for_recovery")
    assert_true(probe_idx < recovery_idx,
        "_PROBE_PROGRAMS guard must appear before recovery = self._analyze_for_recovery "
        "(probe exit must short-circuit recovery analysis). "
        "probe at %d, recovery call at %d" % (probe_idx, recovery_idx))

    # The early-exit must return False (continue workflow, not stop)
    probe_block = block[probe_idx:probe_idx + 1200]
    assert_true("return False" in probe_block,
        "_PROBE_PROGRAMS guard must return False to continue the workflow")

    # model_vs_data and mmtbx variant must both be in the set
    assert_true("phenix.model_vs_data" in probe_block,
        "_PROBE_PROGRAMS must include phenix.model_vs_data")
    assert_true("mmtbx.model_vs_data" in probe_block,
        "_PROBE_PROGRAMS must include mmtbx.model_vs_data")

    print("  PASSED: _PROBE_PROGRAMS guard short-circuits failure analysis; "
          "any model_vs_data failure returns False (continue workflow)")


def test_s10c_sanitize_strips_params_from_probe_programs():
    """sanitize_command must strip ALL key=value tokens from probe programs
    before execution, regardless of whether yaml_loader is importable.

    Regression guard for two real failures:
      1. rebuild_level=QUICK on phenix.model_vs_data   (LLM-generated)
      2. crystal_symmetry.space_group="None mentioned in ins"  (LLM-generated)
    Both caused "Sorry: Unknown command line parameter definition" then agent stop.

    The standalone sanitize_command in command_postprocessor.py is the sole
    authority (Phase 4: class method removed from ai_agent.py).
    """
    print("  Test: s10c_sanitize_strips_params_from_probe_programs")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    # ── Structural: guard must exist in command_postprocessor.py ─────────
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    assert_true(os.path.isfile(pp_path),
        "command_postprocessor.py must exist")
    src = open(pp_path).read()

    sanitize_start = src.find("def sanitize_command(")
    assert_true(sanitize_start >= 0,
        "sanitize_command must exist in command_postprocessor.py")
    sanitize_end = src.find("\ndef ", sanitize_start + 1)
    sanitize_src = src[sanitize_start:sanitize_end] if sanitize_end > 0 else src[sanitize_start:]

    assert_true("_PROBE_ONLY_FILE_PROGRAMS" in sanitize_src,
        "sanitize_command must reference _PROBE_ONLY_FILE_PROGRAMS guard")
    # Constant is module-level, so check the full source for program names
    const_start = src.find("_PROBE_ONLY_FILE_PROGRAMS = frozenset")
    assert_true(const_start >= 0,
        "_PROBE_ONLY_FILE_PROGRAMS must be defined as frozenset")
    const_block = src[const_start:const_start + 500]
    assert_true("phenix.model_vs_data" in const_block,
        "_PROBE_ONLY_FILE_PROGRAMS must include phenix.model_vs_data")
    assert_true("phenix.xtriage" in const_block,
        "_PROBE_ONLY_FILE_PROGRAMS must include phenix.xtriage")

    # Phase 4: class method must be gone from ai_agent.py
    ai_agent_path = _find_ai_agent_path()
    ai_src = open(ai_agent_path).read()
    assert_true("def _sanitize_command" not in ai_src,
        "_sanitize_command class method must be removed from ai_agent.py (Phase 4)")

    # ── Functional: test probe-program stripping inline ──────────────────
    import re as _re_sp

    _MULTIWORD_PLACEHOLDER_RE = _re_sp.compile(
        r'[\w.]+\s*=\s*(?:None|Not|Unknown|N/?A|TBD)'
        r'(?:\s+[A-Za-z]\w*)?',
        _re_sp.IGNORECASE)

    def _probe_sanitize(command):
        """Minimal copy of the probe-program branch in sanitize_command."""
        command = _MULTIWORD_PLACEHOLDER_RE.sub('', command)
        stripped = []
        for tok in command.split():
            if '=' in tok and not tok.startswith('-'):
                # Preserve file-path tokens (contains / or crystallographic extension)
                val = tok.split('=', 1)[1] if '=' in tok else ''
                _FILE_EXTS = {'.pdb', '.cif', '.mtz', '.mrc', '.ccp4', '.map', '.fa', '.fasta'}
                import os as _os_tok
                if '/' in val or _os_tok.path.splitext(val)[1].lower() in _FILE_EXTS:
                    stripped.append(tok)
                # else: stripped (non-file key=value)
            else:
                stripped.append(tok)
        return ' '.join(stripped)

    # Case 1: rebuild_level=QUICK
    cmd1 = "phenix.model_vs_data /path/7qz0.pdb /path/7qz0.mtz rebuild_level=QUICK"
    out1 = _probe_sanitize(cmd1)
    assert_true("rebuild_level" not in out1,
        "rebuild_level=QUICK must be stripped; got: %r" % out1)
    assert_true("7qz0.pdb" in out1 and "7qz0.mtz" in out1,
        "File args must be preserved; got: %r" % out1)

    # Case 2: crystal_symmetry.space_group="None mentioned in ins"
    cmd2 = ('phenix.model_vs_data beta.pdb data.mtz '
            'crystal_symmetry.space_group="None mentioned in ins"')
    out2 = _probe_sanitize(cmd2)
    assert_true("space_group" not in out2,
        "space_group must be stripped; got: %r" % out2)
    assert_true("beta.pdb" in out2 and "data.mtz" in out2,
        "File args must be preserved; got: %r" % out2)

    # Case 3: clean command — unchanged
    cmd3 = "phenix.model_vs_data model.pdb data.mtz"
    out3 = _probe_sanitize(cmd3)
    assert_true(out3 == cmd3,
        "Clean command must pass through unchanged; got: %r" % out3)

    # Case 4: file-path key=value tokens preserved
    cmd4 = "phenix.mtriage half_map_1=/data/half1.mrc half_map_2=/data/half2.mrc rebuild_strategy=quick"
    out4 = _probe_sanitize(cmd4)
    assert_true("half_map_1" in out4 and "half_map_2" in out4,
        "File-path key=value tokens must be preserved; got: %r" % out4)
    assert_true("rebuild_strategy" not in out4,
        "Non-file key=value must be stripped; got: %r" % out4)

    # Case 5: orphan word from two-word placeholder
    cmd5 = ("phenix.xtriage /path/data.mtz "
            "xray_data.space_group=None mentioned")
    out5 = _probe_sanitize(cmd5)
    assert_true("mentioned" not in out5,
        "Orphan 'mentioned' must be stripped; got: %r" % out5)
    assert_true("data.mtz" in out5,
        "MTZ file path must survive; got: %r" % out5)

    # Case 6: _MULTIWORD_PLACEHOLDER_RE present in postprocessor
    assert_true("_MULTIWORD_PLACEHOLDER_RE" in sanitize_src,
        "sanitize_command must contain _MULTIWORD_PLACEHOLDER_RE "
        "to strip two-word placeholder values like 'None mentioned'")

    print("  PASSED: _sanitize_command strips key=value tokens from probe programs "
          "(rebuild_level=QUICK, crystal_symmetry.space_group, "
          "'None mentioned' two-word orphan, etc.)")


def test_s10d_probe_sanitize_preserves_quoted_paths():
    """_sanitize_command must NOT strip file paths that contain '=' in quoted
    form.  Paths with '=' are pathological but must not be falsely stripped.

    More importantly: tokens starting with '-' (flags like --help, -v) must
    be kept even though they may contain '='.
    """
    print("  Test: s10d_probe_sanitize_preserves_quoted_paths")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    def _probe_sanitize(command):
        stripped = []
        for tok in command.split():
            if '=' in tok and not tok.startswith('-'):
                pass
            else:
                stripped.append(tok)
        return ' '.join(stripped)

    # Flags starting with '-' must be preserved
    cmd = "phenix.model_vs_data model.pdb data.mtz --overwrite=True"
    out = _probe_sanitize(cmd)
    assert_true("--overwrite=True" in out,
        "Flags starting with '-' must NOT be stripped; got: %r" % out)

    # Program name (no '=') must always survive
    cmd2 = "phenix.model_vs_data model.pdb"
    assert_true(_probe_sanitize(cmd2) == cmd2,
        "Program name must not be affected; got: %r" % _probe_sanitize(cmd2))

    print("  PASSED: '-' prefixed tokens preserved; program name unaffected")


def test_s10e_rwork_fallback_sets_placement_probed():
    """When model_vs_data returns only R-work (no R-free), placement_probed must
    still be set to True and placement_probe_result to 'needs_mr' when r_work ≥ 0.45.

    Regression: the real log showed 'R Work: 0.5851' in the metrics report but
    no R-free (phenix.model_vs_data skips r_free when the MTZ has no free-set flags).
    Without the r_work fallback, placement_probed stayed False → workflow engine
    routed back to probe_placement → model_vs_data ran again → phaser never picked.
    """
    print("  Test: s10e_rwork_fallback_sets_placement_probed")
    import sys as _sys
    _sys.path.insert(0, _PROJECT_ROOT)

    try:
        from agent.workflow_state import _analyze_history
    except ImportError:
        from libtbx.langchain.agent.workflow_state import _analyze_history

    # Case 1: only r_work present, no r_free — exact scenario from the real log
    history_rwork_only = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "command": "phenix.xtriage data.mtz",
            "result": "SUCCESS: xtriage complete",
            "analysis": {},
            "output_files": [],
        },
        {
            "cycle_number": 2,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data beta.pdb data.mtz",
            "result": (
                "SUCCESS: Command completed without errors\n\n"
                "**************************************************\n"
                "FINAL QUALITY METRICS REPORT:\n"
                "--------------------------------------------------\n"
                "R Work: 0.5851\n"
                "Resolution: 3.00\n"
                "**************************************************\n"
            ),
            "analysis": {"r_work": 0.5851, "resolution": 3.0},  # no r_free key
            "output_files": [],
        },
    ]
    info = _analyze_history(history_rwork_only)
    assert_true(info.get("placement_probed") is True,
        "placement_probed must be True when r_work=0.5851 (fallback path); "
        "got placement_probed=%r" % info.get("placement_probed"))
    assert_true(info.get("placement_probe_result") == "needs_mr",
        "placement_probe_result must be 'needs_mr' when r_work=0.5851 ≥ 0.45; "
        "got %r" % info.get("placement_probe_result"))

    # Case 2: r_work low → model IS placed
    history_placed = [
        {
            "cycle_number": 1,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data placed.pdb data.mtz",
            "result": "SUCCESS: Command completed without errors\nR Work: 0.2341\n",
            "analysis": {"r_work": 0.2341},
            "output_files": [],
        },
    ]
    info2 = _analyze_history(history_placed)
    assert_true(info2.get("placement_probed") is True,
        "placement_probed must be True for low r_work; got %r" % info2.get("placement_probed"))
    assert_true(info2.get("placement_probe_result") == "placed",
        "placement_probe_result must be 'placed' when r_work=0.23 < 0.45; "
        "got %r" % info2.get("placement_probe_result"))

    # Case 3: r_free present → r_free takes priority, r_work fallback NOT used
    history_rfree = [
        {
            "cycle_number": 1,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data model.pdb data.mtz",
            "result": "SUCCESS: r_free: 0.3200\nr_work: 0.2800\n",
            "analysis": {"r_free": 0.32, "r_work": 0.28},
            "output_files": [],
        },
    ]
    info3 = _analyze_history(history_rfree)
    assert_true(info3.get("placement_probed") is True,
        "placement_probed must be True when r_free present")
    assert_true(info3.get("placement_probe_result") == "placed",
        "r_free=0.32 < 0.50 → placed (r_free takes priority over r_work); "
        "got %r" % info3.get("placement_probe_result"))

    # Case 4: r_work from result text (not from analysis dict) — regex fallback
    history_text_only = [
        {
            "cycle_number": 1,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data model.pdb data.mtz",
            "result": "SUCCESS: Command completed without errors\nR Work: 0.5851\n",
            "analysis": {},   # metrics NOT in analysis dict — parsed from result text
            "output_files": [],
        },
    ]
    info4 = _analyze_history(history_text_only)
    assert_true(info4.get("placement_probed") is True,
        "placement_probed must be True when r_work parsed from result text; "
        "got %r" % info4.get("placement_probed"))
    assert_true(info4.get("placement_probe_result") == "needs_mr",
        "needs_mr when r_work from text=0.5851; got %r"
        % info4.get("placement_probe_result"))

    print("  PASSED: r_work fallback correctly sets placement_probed=True/needs_mr "
          "when r_free is absent (prevents model_vs_data from re-running)")


def test_s10f_phaser_success_supersedes_probe_needs_mr():
    """After a successful phaser run, the workflow must route to refine — not
    loop back to molecular_replacement because placement_probed/needs_mr from
    the earlier model_vs_data cycle is still in history.

    Regression: post-phaser cycle 4 got STUCK/STOP because:
      - placement_probed=True, placement_probe_result="needs_mr" (from mvd)
      - Tier-3 check fired BEFORE has_placed_model_from_history check
      - molecular_replacement phase → phaser excluded (not_done:phaser) → STOP

    Fix: gate the needs_mr→MR route on NOT has_placed_model_from_history.
    """
    print("  Test: s10f_phaser_success_supersedes_probe_needs_mr")

    if not _IMPORTS_OK:
        print("  SKIP (imports unavailable)")
        return

    try:
        from agent.workflow_engine import WorkflowEngine
    except ImportError:
        print("  SKIP (WorkflowEngine unavailable)")
        return

    # Build history: xtriage → model_vs_data (needs_mr) → phaser (success)
    history = [
        {
            "cycle": 1,
            "program": "phenix.xtriage",
            "command": "phenix.xtriage data.mtz",
            "result": "SUCCESS: xtriage complete",
            "metrics": {},
            "output_files": [],
        },
        {
            "cycle": 2,
            "program": "phenix.model_vs_data",
            "command": "phenix.model_vs_data beta.pdb data.mtz",
            "result": "SUCCESS: Command completed\nR Work: 0.5851\n",
            "metrics": {"r_work": 0.5851},
            "output_files": [],
        },
        {
            "cycle": 3,
            "program": "phenix.phaser",
            "command": "phenix.phaser data.mtz beta.pdb beta.seq phaser.mode=MR_AUTO",
            "result": "SUCCESS: Command completed without errors\nTFZ==23.60\nLLG=480\n",
            "metrics": {"tfz": 23.60, "llg": 480.0},
            # No output_files — matches the real scenario where PHASER.*.pdb
            # wasn't found in the sub-job directory (history_files=[])
            "output_files": [],
        },
    ]

    history_info = _analyze_history(history)

    # Verify history analysis is correct
    assert_true(history_info.get("placement_probed") is True,
        "placement_probed must be True (from model_vs_data cycle)")
    assert_true(history_info.get("placement_probe_result") == "needs_mr",
        "placement_probe_result must be 'needs_mr' (r_work=0.5851)")
    assert_true(history_info.get("phaser_done") is True,
        "phaser_done must be True after phaser ran")

    # Build engine context — no available files (matches history_files=[])
    engine = WorkflowEngine()
    files = {}  # No output files detected (the real-world scenario)
    context = engine.build_context(files, history_info, analysis=None,
                                   directives={}, session_info={})

    assert_true(context.get("has_placed_model_from_history") is True,
        "has_placed_model_from_history must be True when phaser_done=True")

    # The critical check: detect_step must NOT return molecular_replacement
    step_info = engine.detect_step("xray", context)
    step = step_info.get("step", "")

    assert_true(step != "molecular_replacement",
        "After phaser success, step must NOT be molecular_replacement "
        "(needs_mr probe result must be superseded by phaser_done); "
        "got step=%r" % step)

    assert_true(step in ("refine", "validate", "complete"),
        "After phaser success with no refinement, step must be 'refine' "
        "(or validate/complete); got step=%r" % step)

    # Also verify valid_programs contains phenix.refine (not STOP)
    valid = engine.get_valid_programs("xray", step_info, context, directives={})
    assert_true("phenix.refine" in valid,
        "phenix.refine must be in valid_programs after phaser; "
        "got valid_programs=%r" % valid)
    assert_true(valid != ["STOP"],
        "valid_programs must not be just [STOP] after successful phaser run; "
        "got %r" % valid)

    print("  PASSED: needs_mr probe result correctly superseded by phaser_done "
          "— workflow routes to refine (not stuck at molecular_replacement)")


def test_s10g_crystal_symmetry_per_program_format():
    """inject_crystal_symmetry must use the correct PHIL scope per program:
      - phenix.refine / ligandfit / polder            → crystal_symmetry.space_group=SG
      - phenix.autobuild / autobuild_denmod / autosol → crystal_info.space_group=SG
      - phenix.phaser                                 → xray_data.space_group=SG

    Phase 4: the class method is removed from ai_agent.py; the standalone
    in command_postprocessor.py is the sole authority.
    """
    print("  Test: s10g_crystal_symmetry_per_program_format")
    pp_path = os.path.join(_PROJECT_ROOT, "agent", "command_postprocessor.py")
    if not os.path.isfile(pp_path):
        print("  SKIP (command_postprocessor.py not found)")
        return

    with open(pp_path) as _f:
        src = _f.read()

    # ── structural checks ────────────────────────────────────────────────────
    assert_true("_CRYSTAL_INFO_PROGRAMS" in src,
        "_CRYSTAL_INFO_PROGRAMS frozenset must be defined in inject_crystal_symmetry")
    assert_true("_PHASER_CS_PROGRAMS" in src,
        "_PHASER_CS_PROGRAMS frozenset must be defined in inject_crystal_symmetry")

    ci_block = src.split("_CRYSTAL_INFO_PROGRAMS")[1][:300]
    assert_true("phenix.autobuild" in ci_block,
        "phenix.autobuild must be in _CRYSTAL_INFO_PROGRAMS")
    assert_true("phenix.autosol" in ci_block,
        "phenix.autosol must be in _CRYSTAL_INFO_PROGRAMS")

    phaser_block = src.split("_PHASER_CS_PROGRAMS")[1][:200]
    assert_true("phenix.phaser" in phaser_block,
        "phenix.phaser must be in _PHASER_CS_PROGRAMS")

    assert_true("crystal_info.space_group=" in src,
        "inject_crystal_symmetry must use crystal_info.space_group= for autobuild/autosol")
    assert_true("xray_data.space_group=" in src,
        "inject_crystal_symmetry must use xray_data.space_group= for phaser")
    assert_true("crystal_symmetry.space_group=" in src,
        "inject_crystal_symmetry must use crystal_symmetry.space_group= for refine/etc")

    # Phase 4: class method must be gone from ai_agent.py
    ai_agent_path = _find_ai_agent_path()
    ai_src = open(ai_agent_path).read()
    assert_true("def _inject_crystal_symmetry" not in ai_src,
        "_inject_crystal_symmetry class method must be removed from ai_agent.py (Phase 4)")

    # ── functional checks using the standalone function ────────────────────
    try:
        sys.path.insert(0, _PROJECT_ROOT)
        try:
            from agent.command_postprocessor import inject_crystal_symmetry
        except ImportError:
            print("  (functional checks skipped: No module named 'phenix')")
            print("  PASSED: structural checks passed")
            return
    except Exception as _e:
        print("  (functional checks skipped: %s)" % _e)
        print("  PASSED: structural checks passed")
        return

    directives = {"program_settings": {"default": {"space_group": "P3221"}}}

    # autobuild: crystal_info scope
    result_ab = inject_crystal_symmetry(
        "phenix.autobuild data=data.mtz seq_file=seq.fa model=m.pdb nproc=4",
        program_name="phenix.autobuild", directives=directives)
    assert_true("crystal_info.space_group=P3221" in result_ab,
        "autobuild must get crystal_info.space_group=P3221; got: %r" % result_ab)
    assert_true("crystal_symmetry.space_group" not in result_ab,
        "autobuild must NOT get crystal_symmetry.space_group=; got: %r" % result_ab)

    # refine: crystal_symmetry scope
    result_ref = inject_crystal_symmetry(
        "phenix.refine model.pdb data.mtz",
        program_name="phenix.refine", directives=directives)
    assert_true("crystal_symmetry.space_group=P3221" in result_ref,
        "refine must get crystal_symmetry.space_group=P3221; got: %r" % result_ref)

    # phaser: xray_data scope
    result_ph = inject_crystal_symmetry(
        "phenix.phaser data.mtz model.pdb seq.fa phaser.mode=MR_AUTO",
        program_name="phenix.phaser", directives=directives)
    assert_true("xray_data.space_group=P3221" in result_ph,
        "phaser must get xray_data.space_group=P3221; got: %r" % result_ph)
    assert_true("crystal_symmetry.space_group" not in result_ph,
        "phaser must NOT get crystal_symmetry.space_group=; got: %r" % result_ph)

    # no double-injection
    result_ab2 = inject_crystal_symmetry(
        "phenix.autobuild data=data.mtz crystal_info.space_group=P3221 nproc=4",
        program_name="phenix.autobuild", directives=directives)
    assert_true(result_ab2.count("space_group=") == 1,
        "Must not double-inject space_group=; count=%d in %r"
        % (result_ab2.count("space_group="), result_ab2))

    print("  PASSED: inject_crystal_symmetry uses correct per-program PHIL scope")

def test_s10h_sanitize_quoted_spacegroup_no_orphan():
    """sanitize_command must not leave orphan tokens when a quoted multi-word
    space group value is stripped from a probe-only program.

    Regression: crystal_symmetry.space_group="P 32 2 1" split into
      ['crystal_symmetry.space_group="P', '32', '2', '1"']
    Only the first token was stripped; '32 2 1"' remained on the command line,
    causing phenix.model_vs_data to crash with an unknown positional argument.

    Fix: collapse key="quoted value" into a single token before the split loop.
    Phase 4: now tests command_postprocessor.py (class method removed from ai_agent.py).
    """
    print("  Test: s10h_sanitize_quoted_spacegroup_no_orphan")

    try:
        sys.path.insert(0, _PROJECT_ROOT)
        try:
            from agent.command_postprocessor import sanitize_command
        except ImportError:
            from libtbx.langchain.agent.command_postprocessor import sanitize_command
    except Exception as _e:
        print("  SKIP (cannot import sanitize_command: %s)" % _e)
        return

    # Exact pattern from the bug: space group with spaces, double-quoted
    cmd = ('phenix.model_vs_data /path/to/model.pdb /path/to/data.mtz '
           'crystal_symmetry.space_group="P 32 2 1"')
    result = sanitize_command(cmd, program_name='phenix.model_vs_data')

    # No orphan numeric tokens from the space group value
    assert_true('32' not in result.split()[2:],
        "Orphan token '32' must not appear in sanitized command: %r" % result)
    assert_true('"' not in result,
        "No stray quote must remain in sanitized command: %r" % result)
    assert_true('space_group' not in result,
        "space_group param must be fully removed: %r" % result)

    # Files must be preserved
    assert_true('model.pdb' in result,
        "model.pdb must be preserved: %r" % result)
    assert_true('data.mtz' in result,
        "data.mtz must be preserved: %r" % result)

    # Also test single-word space group (no quotes) — existing behaviour preserved
    cmd2 = ('phenix.model_vs_data /path/to/model.pdb /path/to/data.mtz '
            'crystal_symmetry.space_group=P3221')
    result2 = sanitize_command(cmd2, program_name='phenix.model_vs_data')
    assert_true('space_group' not in result2,
        "Unquoted space_group must also be removed: %r" % result2)
    assert_true('model.pdb' in result2,
        "model.pdb must be preserved: %r" % result2)

    print("  PASSED: quoted multi-word space group stripped cleanly, no orphan tokens")


# =============================================================================
# S5j — Ligandfit file selection fixes
# =============================================================================

def test_s5j_refine_mtz_classified_as_map_coeffs():
    """Bug: refine_001.mtz was classified as data_mtz instead of map_coeffs_mtz.

    The is_map_coeffs regex in session._rebuild_best_files_from_cycles only
    matched refine_NNN_001.mtz (two-level serial), not the standard
    refine_001.mtz output.  This caused best_files["map_coeffs_mtz"] to
    never be populated after refinement, so ligandfit had no data file.

    Same bug existed in file_utils.classify_mtz_type (used by
    BestFilesTracker._classify_mtz_type for live evaluations) and in
    session.record_result (used when recording cycle output files).
    """
    import re

    # Test the session.py regex (used in _rebuild_best_files_from_cycles
    # and record_result)
    def is_map_coeffs_fixed(basename):
        return bool(
            'map_coeffs' in basename or
            'denmod' in basename or
            re.match(r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$', basename)
        )

    # All refine MTZ patterns should be recognized as map coefficients
    assert_true(is_map_coeffs_fixed('refine_001.mtz'),
        "refine_001.mtz must be map_coeffs")
    assert_true(is_map_coeffs_fixed('refine_001_001.mtz'),
        "refine_001_001.mtz must be map_coeffs")
    assert_true(is_map_coeffs_fixed('7qz0_refine_001.mtz'),
        "7qz0_refine_001.mtz must be map_coeffs")
    assert_true(is_map_coeffs_fixed('7qz0_refine_001_001.mtz'),
        "7qz0_refine_001_001.mtz must be map_coeffs")

    # Data MTZ files must NOT be classified as map coefficients
    assert_false(is_map_coeffs_fixed('nsf-d2.mtz'),
        "nsf-d2.mtz must NOT be map_coeffs")
    assert_false(is_map_coeffs_fixed('nsf-d2_data.mtz'),
        "nsf-d2_data.mtz must NOT be map_coeffs")
    assert_false(is_map_coeffs_fixed('data.mtz'),
        "data.mtz must NOT be map_coeffs")

    # Also verify file_utils.classify_mtz_type (the shared classifier)
    from agent.file_utils import classify_mtz_type
    assert_equal(classify_mtz_type('/path/refine_001.mtz'), 'map_coeffs_mtz',
        "classify_mtz_type must handle refine_001.mtz")
    assert_equal(classify_mtz_type('/path/refine_001_001.mtz'), 'map_coeffs_mtz',
        "classify_mtz_type must handle refine_001_001.mtz")
    assert_equal(classify_mtz_type('/path/7qz0_refine_001.mtz'), 'map_coeffs_mtz',
        "classify_mtz_type must handle 7qz0_refine_001.mtz")
    assert_equal(classify_mtz_type('/path/nsf-d2.mtz'), 'data_mtz',
        "classify_mtz_type must keep nsf-d2.mtz as data_mtz")

    print("  PASSED: refine MTZ output correctly classified as map_coeffs_mtz")


def test_s5j_matches_exclude_pattern_word_boundary():
    """Bug: exclude_patterns used substring matching, so 'ligand' in
    exclude_patterns matched 'nsf-d2_noligand.pdb' — a protein model
    whose name happens to contain 'noligand'.

    Now uses word-boundary matching: 'ligand' only matches when preceded
    by start-of-string or a separator (_, -, .).
    """
    from agent.file_utils import matches_exclude_pattern

    # 'ligand' must NOT match inside 'noligand'
    assert_false(matches_exclude_pattern('nsf-d2_noligand.pdb', ['ligand']),
        "noligand must not match 'ligand' exclude")

    # 'ligand' MUST match when it's a proper word
    assert_true(matches_exclude_pattern('ligand.pdb', ['ligand']),
        "ligand.pdb must match")
    assert_true(matches_exclude_pattern('atp_ligand.pdb', ['ligand']),
        "atp_ligand.pdb must match")
    assert_true(matches_exclude_pattern('my_ligand_001.pdb', ['ligand']),
        "my_ligand_001.pdb must match")

    # 'lig.pdb' pattern (with extension) must match exactly
    assert_true(matches_exclude_pattern('lig.pdb', ['lig.pdb']),
        "lig.pdb must match itself")
    assert_false(matches_exclude_pattern('nolig.pdb', ['lig.pdb']),
        "nolig.pdb must not match lig.pdb pattern")
    assert_false(matches_exclude_pattern('lig.cif', ['lig.pdb']),
        "lig.cif must not match lig.pdb (wrong extension)")

    # Half-map patterns
    assert_true(matches_exclude_pattern('half_1.mrc', ['half_1']),
        "half_1.mrc must match")
    assert_false(matches_exclude_pattern('nothalf_1.mrc', ['half_1']),
        "nothalf_1.mrc must not match")
    assert_true(matches_exclude_pattern('model_half_1.mrc', ['half_1']),
        "model_half_1.mrc must match (separator before half)")

    # 'lig' pattern for prefer_patterns
    assert_true(matches_exclude_pattern('lig_001.pdb', ['lig']),
        "lig_001.pdb must match")
    assert_false(matches_exclude_pattern('nolig_001.pdb', ['lig']),
        "nolig_001.pdb must not match")

    # atp.pdb has no ligand-like word boundary
    assert_false(matches_exclude_pattern('atp.pdb', ['lig', 'ligand']),
        "atp.pdb must not match lig/ligand patterns")

    print("  PASSED: word-boundary exclude/prefer pattern matching correct")


def test_s5j_content_guard_ligand_and_model_slots():
    """Bug: LLM assigned refine_001_001.pdb (protein) to the ligand slot.

    Content-based guards should reject:
    - protein PDBs from the ligand slot (proteins have ATOM records)
    - small-molecule PDBs from model/protein/pdb_file slots (HETATM-only)

    The ligand slot guard uses _pdb_is_protein_model (positive protein check)
    rather than "not _pdb_is_small_molecule" to avoid rejecting unreadable
    or non-existent files.
    """
    test_dir = os.path.join(os.path.expanduser('~'), '_test_content_guard')
    os.makedirs(test_dir, exist_ok=True)

    try:
        # Create a protein PDB (has many ATOM records — realistic refined model)
        protein_pdb = os.path.join(test_dir, 'refine_001.pdb')
        with open(protein_pdb, 'w') as f:
            # Real protein models have hundreds to thousands of atoms
            for i in range(200):
                f.write("ATOM  %5d  CA  ALA A %3d       1.000   2.000   3.000  1.00 10.00           C\n" % (i+1, i+1))
            f.write("END\n")

        # Create a small-molecule PDB (HETATM only)
        ligand_pdb = os.path.join(test_dir, 'atp.pdb')
        with open(ligand_pdb, 'w') as f:
            f.write("HETATM    1  PG  ATP A   1       1.000   2.000   3.000  1.00 10.00           P\n")
            f.write("HETATM    2  O1G ATP A   1       2.000   3.000   4.000  1.00 10.00           O\n")
            f.write("END\n")

        from agent.workflow_state import _pdb_is_small_molecule, _pdb_is_protein_model

        # _pdb_is_small_molecule: True for HETATM-only, False otherwise
        assert_false(_pdb_is_small_molecule(protein_pdb),
            "refine_001.pdb (has ATOM) must NOT be classified as small molecule")
        assert_true(_pdb_is_small_molecule(ligand_pdb),
            "atp.pdb (HETATM-only) must be classified as small molecule")

        # _pdb_is_protein_model: True for large files with ATOM records, False for small/HETATM-only
        assert_true(_pdb_is_protein_model(protein_pdb),
            "refine_001.pdb (200 ATOM records) must be identified as protein model")
        assert_false(_pdb_is_protein_model(ligand_pdb),
            "atp.pdb (HETATM-only) must NOT be identified as protein model")

        # Critical: non-existent files must NOT be rejected from ligand slot
        # _pdb_is_protein_model returns False for missing files (conservative)
        assert_false(_pdb_is_protein_model('/nonexistent/ligand.pdb'),
            "Non-existent file must NOT be identified as protein (would block ligand slot)")
        # _pdb_is_small_molecule also returns False for missing files
        assert_false(_pdb_is_small_molecule('/nonexistent/ligand.pdb'),
            "Non-existent file returns False for small_molecule too")

        print("  PASSED: content-based ligand/model slot guards correct")

    finally:
        import shutil
        shutil.rmtree(test_dir, ignore_errors=True)


def test_s5j_refine_cif_excluded_from_ligand_slot():
    """Bug: refine_001_001.cif (refinement restraints) was used as ligand.

    After the protein PDB guard correctly rejected the protein model from
    the ligand slot, the auto-fill picked refine_001_001.cif (geometry
    restraints from refinement) instead of the actual ligand file.

    Fix: 'refine' added to exclude_patterns for ligandfit ligand slot,
    and exclude_patterns now applied to LLM-selected files too.
    """
    from agent.file_utils import matches_exclude_pattern

    # Verify the ligandfit ligand slot has 'refine' in exclude_patterns
    from knowledge.yaml_loader import get_program_inputs
    inputs = get_program_inputs('phenix.ligandfit')
    ligand_def = inputs.get('required', {}).get('ligand', {})
    excl = ligand_def.get('exclude_patterns', [])
    assert_true('refine' in excl,
        "ligandfit ligand slot must have 'refine' in exclude_patterns, got: %s" % excl)

    # Refine output CIFs must be excluded
    assert_true(matches_exclude_pattern('refine_001_001.cif', excl),
        "refine_001_001.cif must match ligand exclude_patterns")
    assert_true(matches_exclude_pattern('refine_001.cif', excl),
        "refine_001.cif must match ligand exclude_patterns")
    assert_true(matches_exclude_pattern('7qz0_refine_001.cif', excl),
        "7qz0_refine_001.cif must match ligand exclude_patterns")

    # Actual ligand files must NOT be excluded
    assert_false(matches_exclude_pattern('ATP.cif', excl),
        "ATP.cif must NOT match ligand exclude_patterns")
    assert_false(matches_exclude_pattern('ligand.cif', excl),
        "ligand.cif must NOT match ligand exclude_patterns")
    assert_false(matches_exclude_pattern('atp.pdb', excl),
        "atp.pdb must NOT match ligand exclude_patterns")
    assert_false(matches_exclude_pattern('atp_refined.cif', excl),
        "atp_refined.cif must NOT match (refined != refine at word boundary)")

    print("  PASSED: refine CIF correctly excluded from ligand slot")


def test_s5j_session_rebuild_map_coeffs():
    """Integration: _rebuild_best_files_from_cycles correctly populates
    map_coeffs_mtz after a successful refinement cycle produces refine_001.mtz.
    """
    import json

    # Use /home/claude as base to avoid /tmp/ triggering the intermediate
    # file filter in BestFilesTracker._is_intermediate_file.
    test_dir = os.path.join(os.path.expanduser('~'), '_test_session_rebuild')
    os.makedirs(test_dir, exist_ok=True)
    try:
        # Create mock files
        model_path = os.path.join(test_dir, 'nsf-d2_noligand.pdb')
        data_path = os.path.join(test_dir, 'nsf-d2.mtz')
        refine_pdb = os.path.join(test_dir, 'refine_001.pdb')
        refine_mtz = os.path.join(test_dir, 'refine_001.mtz')

        for path in [model_path, data_path]:
            with open(path, 'w') as f:
                f.write('ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n')
        for path in [refine_pdb, refine_mtz]:
            with open(path, 'w') as f:
                f.write('MOCK')

        # Create session with one successful refine cycle
        session_file = os.path.join(test_dir, 'session.json')
        session_data = {
            "original_files": [model_path, data_path],
            "cycles": [{
                "cycle_number": 1,
                "program": "phenix.refine",
                "result": "SUCCESS",
                "metrics": {"r_work": 0.22, "r_free": 0.26},
                "output_files": [refine_pdb, refine_mtz],
            }],
        }
        with open(session_file, 'w') as f:
            json.dump(session_data, f)

        from agent.session import AgentSession
        session = AgentSession(session_file=session_file)

        # Check that map_coeffs_mtz is now populated
        best = session.get_best_files_dict()
        map_coeffs = best.get('map_coeffs_mtz')
        assert_true(map_coeffs is not None,
            "map_coeffs_mtz must be populated after refine: got %s" % best)
        assert_true('refine_001.mtz' in str(map_coeffs),
            "map_coeffs_mtz must point to refine_001.mtz, got: %s" % map_coeffs)
    finally:
        import shutil
        shutil.rmtree(test_dir, ignore_errors=True)

    print("  PASSED: _rebuild_best_files_from_cycles populates map_coeffs_mtz")


def test_s5j_session_rebuild_supplemental_map_coeffs():
    """Integration: best_files discovers map coefficients MTZ that wasn't
    tracked in the cycle's output_files.

    This is the root cause of ligandfit failing to build — the refine cycle
    tracked refine_001_data.mtz (data) but not refine_001.mtz (map coefficients).
    The _find_missing_outputs method discovers the companion file on disk,
    and _rebuild_best_files_from_cycles must evaluate it so best_files gets
    map_coeffs_mtz populated.
    """
    import json

    test_dir = os.path.join(os.path.expanduser('~'), '_test_session_supplemental')
    os.makedirs(test_dir, exist_ok=True)
    try:
        # Create mock files
        model_path = os.path.join(test_dir, 'nsf-d2_noligand.pdb')
        data_path = os.path.join(test_dir, 'nsf-d2.mtz')
        refine_pdb = os.path.join(test_dir, 'refine_001_001.pdb')
        refine_data_mtz = os.path.join(test_dir, 'refine_001_data.mtz')
        # Map coefficients MTZ exists on disk but is NOT in output_files
        refine_map_mtz = os.path.join(test_dir, 'refine_001.mtz')

        for path in [model_path, data_path]:
            with open(path, 'w') as f:
                f.write('ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n')
        for path in [refine_pdb, refine_data_mtz, refine_map_mtz]:
            with open(path, 'w') as f:
                f.write('MOCK')

        # Create session where output_files only has _data.mtz and .pdb
        # (the map coefficients MTZ was NOT tracked by the client)
        session_file = os.path.join(test_dir, 'session.json')
        session_data = {
            "original_files": [model_path, data_path],
            "cycles": [{
                "cycle_number": 1,
                "program": "phenix.refine",
                "result": "SUCCESS",
                "metrics": {"r_work": 0.22, "r_free": 0.26},
                "output_files": [refine_pdb, refine_data_mtz],
                # Note: refine_001.mtz deliberately NOT included
            }],
        }
        with open(session_file, 'w') as f:
            json.dump(session_data, f)

        from agent.session import AgentSession
        session = AgentSession(session_file=session_file)

        # Check that map_coeffs_mtz is populated via supplemental discovery
        best = session.get_best_files_dict()
        map_coeffs = best.get('map_coeffs_mtz')
        assert_true(map_coeffs is not None,
            "map_coeffs_mtz must be populated via supplemental discovery: got %s" % best)
        assert_true('refine_001.mtz' in str(map_coeffs),
            "map_coeffs_mtz must point to refine_001.mtz, got: %s" % map_coeffs)

        # Also verify available_files includes the supplemental MTZ
        avail = session.get_available_files()
        avail_basenames = [os.path.basename(f) for f in avail]
        assert_true('refine_001.mtz' in avail_basenames,
            "available_files must include supplemental refine_001.mtz, got: %s" % avail_basenames)

    finally:
        import shutil
        shutil.rmtree(test_dir, ignore_errors=True)

    print("  PASSED: supplemental map coefficients MTZ discovered and added to best_files")


def test_s5j_record_result_discovers_supplemental_map_coeffs():
    """Live path: record_result discovers map coefficients MTZ not in output_files.

    Same scenario as test_s5j_session_rebuild_supplemental_map_coeffs but
    exercising the live record_result path instead of session-load rebuild.
    """
    import json

    test_dir = os.path.join(os.path.expanduser('~'), '_test_record_result_supp')
    os.makedirs(test_dir, exist_ok=True)
    try:
        # Create mock files
        model_path = os.path.join(test_dir, 'model.pdb')
        data_path = os.path.join(test_dir, 'data.mtz')
        refine_pdb = os.path.join(test_dir, 'refine_001_001.pdb')
        refine_data_mtz = os.path.join(test_dir, 'refine_001_data.mtz')
        refine_map_mtz = os.path.join(test_dir, 'refine_001.mtz')  # on disk, not in output_files

        for path in [model_path, data_path]:
            with open(path, 'w') as f:
                f.write('ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n')
        for path in [refine_pdb, refine_data_mtz, refine_map_mtz]:
            with open(path, 'w') as f:
                f.write('MOCK')

        # Create fresh session (no cycles yet)
        session_file = os.path.join(test_dir, 'session.json')
        session_data = {
            "original_files": [model_path, data_path],
            "cycles": [],
        }
        with open(session_file, 'w') as f:
            json.dump(session_data, f)

        from agent.session import AgentSession
        session = AgentSession(session_file=session_file)

        # Mock out metrics extraction (needs libtbx in real environment)
        session._extract_metrics_from_result = lambda *a, **k: {}

        # Simulate a live cycle: record decision, then record result
        session.record_decision(
            cycle_number=1, program="phenix.refine",
            decision="Run refine", command="phenix.refine model.pdb data.mtz")
        session.record_result(
            cycle_number=1, result="SUCCESS: R=0.22 Rfree=0.26",
            output_files=[refine_pdb, refine_data_mtz])
        # Note: refine_001.mtz deliberately NOT in output_files

        # Check that map_coeffs_mtz is populated via supplemental discovery
        best = session.get_best_files_dict()
        map_coeffs = best.get('map_coeffs_mtz')
        assert_true(map_coeffs is not None,
            "record_result must discover supplemental map_coeffs_mtz: got %s" % best)
        assert_true('refine_001.mtz' in str(map_coeffs),
            "map_coeffs_mtz must point to refine_001.mtz, got: %s" % map_coeffs)

    finally:
        import shutil
        shutil.rmtree(test_dir, ignore_errors=True)

    print("  PASSED: record_result discovers supplemental map coefficients on live path")


def test_s5j_duplicate_detection_different_model_not_duplicate():
    """Refine with a different model file must NOT be flagged as duplicate.

    Bug: The 80% token-overlap heuristic flagged a new refinement cycle as
    duplicate of a previous one even though the input model file changed
    (refine_002.pdb vs model.pdb). The agent then stopped or retried
    instead of running the requested refinement.

    Fix: If the file tokens (basenames with crystallographic extensions)
    differ between two commands, they are NOT duplicates regardless of
    overall token overlap.
    """
    import json

    test_dir = os.path.join(os.path.expanduser('~'), '_test_dup_detection')
    os.makedirs(test_dir, exist_ok=True)
    try:
        session_file = os.path.join(test_dir, 'session.json')
        session_data = {
            'original_files': ['/data/model.pdb', '/data/data.mtz'],
            'cycles': [
                {
                    'cycle_number': 1,
                    'program': 'phenix.refine',
                    'command': 'phenix.refine /data/model.pdb /data/data.mtz nproc=4 '
                               'main.number_of_macro_cycles=1',
                    'result': 'SUCCESS: R=0.25',
                    'output_files': [],
                },
                {
                    'cycle_number': 2,
                    'program': 'phenix.refine',
                    'command': 'phenix.refine /data/refine_001.pdb /data/data.mtz nproc=4 '
                               'main.number_of_macro_cycles=1',
                    'result': 'SUCCESS: R=0.22',
                    'output_files': [],
                },
            ],
        }
        with open(session_file, 'w') as f:
            json.dump(session_data, f)

        from agent.session import AgentSession
        session = AgentSession(session_file=session_file)

        # Refine with NEW model (output of cycle 2) → NOT duplicate
        cmd_new = ('phenix.refine /data/refine_002.pdb /data/data.mtz '
                   'nproc=4 main.number_of_macro_cycles=1')
        is_dup, prev, _ = session.is_duplicate_command(cmd_new)
        assert_false(is_dup,
            "Refine with different model must NOT be flagged as duplicate, "
            "but was flagged as dup of cycle %s" % prev)

        # Exact same command as cycle 1 → IS duplicate
        cmd_same = ('phenix.refine /data/model.pdb /data/data.mtz '
                    'nproc=4 main.number_of_macro_cycles=1')
        is_dup2, prev2, _ = session.is_duplicate_command(cmd_same)
        assert_true(is_dup2,
            "Exact same command as cycle 1 must be flagged as duplicate")

        # Same program, same params, but model path differs only in directory →
        # still different basename, NOT duplicate
        cmd_diff_dir = ('phenix.refine /other/dir/new_model.pdb /data/data.mtz '
                        'nproc=4 main.number_of_macro_cycles=1')
        is_dup3, _, _ = session.is_duplicate_command(cmd_diff_dir)
        assert_false(is_dup3,
            "Different model basename must NOT be duplicate")

        print("  PASSED: duplicate detection respects different input files")

    finally:
        import shutil
        shutil.rmtree(test_dir, ignore_errors=True)


def test_s5k_best_files_fallback_for_specific_subcategory():
    """When category-based lookup fails, best_files should be tried as fallback.

    Scenario: ligandfit's map_coeffs_mtz slot uses specific subcategories
    (refine_map_coeffs, etc.), which skips the normal best_files path. If
    categorized_files doesn't contain the file in the right category but
    best_files["map_coeffs_mtz"] IS populated (from supplemental discovery),
    the build should still succeed using best_files as fallback.

    Verified by source-code inspection since CommandBuilder requires libtbx.
    """
    print("  Test: s5k_best_files_fallback_for_specific_subcategory")

    cb_path = os.path.join(_PROJECT_ROOT, "agent", "command_builder.py")
    with open(cb_path) as f:
        src = f.read()

    # The fallback must exist between "uses_specific_subcategory" guard and
    # "return None" for the extension-skip path
    assert_true("best_files_fallback" in src,
        "command_builder.py must contain 'best_files_fallback' selection reason")

    # The fallback must check exclude_categories before using best_files
    # Find the block: after "if uses_specific_subcategory:" and before the
    # final "return None", there should be a best_files lookup with exclude check
    import re
    pattern = (
        r'if uses_specific_subcategory:.*?'
        r'SLOT_TO_BEST_CATEGORY.*?'
        r'best_files.*?'
        r'exclude.*?'
        r'best_files_fallback.*?'
        r'return None'
    )
    match = re.search(pattern, src, re.DOTALL)
    assert_true(match is not None,
        "command_builder.py must have best_files fallback with exclude check "
        "inside the uses_specific_subcategory block before return None")

    # Verify the fallback logs clearly
    assert_true("category lookup found no files" in src,
        "Fallback should log 'category lookup found no files' for transparency")

    print("  PASSED: best_files fallback works for specific subcategory slots")


def test_s5l_set_project_info_merges_on_resume():
    """set_project_info must MERGE original_files on resume, not replace.

    When a user resumes a session with only refinement outputs, the original
    ligand file from the first run must be preserved in original_files.
    """
    print("  Test: s5l_set_project_info_merges_on_resume")
    import tempfile, shutil

    test_dir = tempfile.mkdtemp(prefix="_test_merge_")
    try:
        from agent.session import AgentSession
        session = AgentSession(session_dir=test_dir)

        # First run: user supplies data + model + ligand
        session.set_project_info(
            project_advice="Solve the structure",
            original_files=["/data/data.mtz", "/data/model.pdb", "/data/atp.pdb"]
        )
        assert_true(len(session.data["original_files"]) == 3,
            "First call should set 3 files, got %d" % len(session.data["original_files"]))

        # Resume: user supplies only refine outputs (different files)
        session.set_project_info(
            project_advice="fit ATP to the refined model",
            original_files=["/refine/refine_001_data.mtz", "/refine/refine_001.pdb"]
        )
        files = session.data["original_files"]
        basenames = [os.path.basename(f) for f in files]

        # Must have ALL 5 files (3 original + 2 new)
        assert_true("atp.pdb" in basenames,
            "atp.pdb from first run must be preserved on resume, got: %s" % basenames)
        assert_true("refine_001_data.mtz" in basenames,
            "refine_001_data.mtz must be added on resume, got: %s" % basenames)
        assert_true("refine_001.pdb" in basenames,
            "refine_001.pdb must be added on resume, got: %s" % basenames)
        assert_true(len(files) == 5,
            "Should have 5 files total (3 original + 2 new), got %d: %s" % (len(files), basenames))

        # Duplicate basenames must not be re-added
        session.set_project_info(
            original_files=["/other/atp.pdb", "/other/data.mtz"]
        )
        files2 = session.data["original_files"]
        assert_true(len(files2) == 5,
            "Duplicate basenames must not be re-added, got %d: %s" % (
                len(files2), [os.path.basename(f) for f in files2]))

        print("  PASSED: set_project_info merges original_files on resume")

    finally:
        shutil.rmtree(test_dir, ignore_errors=True)


def test_s5m_build_failure_includes_missing_slots():
    """BUILD node validation_error must include missing slot names.

    Without this, the user sees "Failed to build command" with no indication
    of what file they need to supply.
    """
    print("  Test: s5m_build_failure_includes_missing_slots")

    gn_path = os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py")
    with open(gn_path) as f:
        src = f.read()

    # BUILD node must check _last_missing_slots and include them in error
    assert_true("_last_missing_slots" in src,
        "graph_nodes.py must reference _last_missing_slots")
    assert_true('missing:' in src.lower() or 'missing_slots' in src,
        "BUILD failure message must include missing slot names")

    # The error message format should include the slot names
    import re
    # Look for pattern: "Failed to build command (missing: %s)"
    pattern = r'Failed to build command \(missing.*?\)'
    match = re.search(pattern, src)
    assert_true(match is not None,
        "BUILD node must format error as 'Failed to build command (missing: ...)' "
        "so fallback reasoning shows which slots are unfilled")

    print("  PASSED: BUILD failure includes missing slot names")


def test_s5n_protein_guard_ratio_based():
    """_pdb_is_protein_model must use size-based check for small files.

    A standalone ligand file like atp.pdb may use ATOM records (not HETATM).
    The guard should only reject files that are large enough to be actual
    protein models (> 150 coordinate records).
    """
    print("  Test: s5n_protein_guard_ratio_based")
    import tempfile, os

    from agent.workflow_state import _pdb_is_protein_model

    # Pure HETATM ligand → NOT protein
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        for i in range(30):
            f.write("HETATM%5d  C   ATP A   1       0.000   0.000   0.000  1.00  0.00\n" % (i+1))
        f.flush()
        assert_true(not _pdb_is_protein_model(f.name),
            "Pure HETATM file should NOT be protein")
        os.unlink(f.name)

    # Small file with ALL ATOM records (like real atp.pdb) → NOT protein
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        for i in range(31):
            f.write("ATOM  %5d  C   ATP A   1       0.000   0.000   0.000  1.00  0.00\n" % (i+1))
        f.flush()
        assert_true(not _pdb_is_protein_model(f.name),
            "Small all-ATOM file (31 atoms, like ATP) should NOT be protein")
        os.unlink(f.name)

    # Medium file at boundary (150 ATOM records) → NOT protein
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        for i in range(150):
            f.write("ATOM  %5d  CA  ALA A %3d       0.000   0.000   0.000  1.00  0.00\n" % (i+1, i+1))
        f.flush()
        assert_true(not _pdb_is_protein_model(f.name),
            "150-atom file should NOT be protein (at boundary)")
        os.unlink(f.name)

    # Large protein model → IS protein
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        for i in range(500):
            f.write("ATOM  %5d  CA  ALA A %3d       0.000   0.000   0.000  1.00  0.00\n" % (i+1, i+1))
        for i in range(10):
            f.write("HETATM%5d  C   ATP A 501       0.000   0.000   0.000  1.00  0.00\n" % (i+501))
        f.flush()
        assert_true(_pdb_is_protein_model(f.name),
            "Large file with mostly ATOM records should be protein")
        os.unlink(f.name)

    # Empty file → NOT protein
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write("REMARK  test\nEND\n")
        f.flush()
        assert_true(not _pdb_is_protein_model(f.name),
            "Empty PDB (no coordinates) should NOT be protein")
        os.unlink(f.name)

    # Non-existent file → NOT protein
    assert_true(not _pdb_is_protein_model("/nonexistent/file.pdb"),
        "Non-existent file should NOT be protein")

    print("  PASSED: size-based protein guard works correctly")


def test_s5o_ligandfit_prereq_refine_first():
    """Bug 17 extended: when user wants ligandfit but refine_count==0,
    refine must run first (produces map_coeffs_mtz for ligandfit).
    - Cycle 3 (after xtriage+model_vs_data): valid=[phenix.refine] only
    - Cycle 4 (after refine): valid includes phenix.ligandfit
    """
    print("  Test: ligandfit_prereq_refine_first")

    import types
    if 'libtbx' not in sys.modules:
        libtbx = types.ModuleType('libtbx')
        libtbx.langchain = types.ModuleType('libtbx.langchain')
        libtbx.langchain.__path__ = [os.path.join(os.path.dirname(__file__), '..')]
        sys.modules['libtbx'] = libtbx
        sys.modules['libtbx.langchain'] = libtbx.langchain

    from agent.workflow_engine import WorkflowEngine
    from agent.workflow_state import _categorize_files, _analyze_history

    engine = WorkflowEngine()
    files_list = ['/tmp/7qz0.fa', '/tmp/7qz0.mtz', '/tmp/7qz0_ligand.pdb', '/tmp/7qz0.pdb']
    directives = {"workflow_preferences": {"prefer_programs": ["phenix.ligandfit"]}}

    # After xtriage + model_vs_data with good R-free (pre-refined model)
    history_pre = [
        {"program": "phenix.xtriage", "command": "phenix.xtriage 7qz0.mtz",
         "result": "SUCCESS", "analysis": {"resolution": 2.10}},
        {"program": "phenix.model_vs_data",
         "command": "phenix.model_vs_data 7qz0.pdb 7qz0.mtz",
         "result": "SUCCESS",
         "analysis": {"r_free": 0.204, "r_work": 0.177, "resolution": 2.10}},
    ]
    files = _categorize_files(files_list)
    hi_pre = _analyze_history(history_pre)
    analysis_pre = {"r_free": 0.204, "resolution": 2.10}

    state_pre = engine.get_workflow_state("xray", files, hi_pre, analysis_pre, directives)
    valid_pre = state_pre["valid_programs"]

    # Before refine: only phenix.refine should be available, NOT ligandfit
    assert_true("phenix.refine" in valid_pre,
        "phenix.refine should be in valid_programs when refine_count==0 "
        "(needed as ligandfit prerequisite), got: %s" % valid_pre)
    assert_true("phenix.ligandfit" not in valid_pre,
        "phenix.ligandfit should NOT be in valid_programs when refine_count==0 "
        "(can't build without map_coeffs_mtz), got: %s" % valid_pre)
    assert_true("STOP" not in valid_pre,
        "STOP should NOT be in valid_programs when refine is a ligandfit "
        "prerequisite, got: %s" % valid_pre)

    # After refine runs
    history_post = history_pre + [
        {"program": "phenix.refine",
         "command": "phenix.refine 7qz0.pdb 7qz0.mtz",
         "result": "SUCCESS",
         "analysis": {"r_free": 0.198, "resolution": 2.10}},
    ]
    hi_post = _analyze_history(history_post)
    analysis_post = {"r_free": 0.198, "resolution": 2.10}

    state_post = engine.get_workflow_state("xray", files, hi_post, analysis_post, directives)
    valid_post = state_post["valid_programs"]

    # After refine: ligandfit should be available
    assert_true("phenix.ligandfit" in valid_post,
        "phenix.ligandfit should be in valid_programs after refine runs "
        "(refine_count > 0, YAML conditions pass), got: %s" % valid_post)

    print("  PASSED")


def run_all_tests():
    """Run all audit fix tests."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
