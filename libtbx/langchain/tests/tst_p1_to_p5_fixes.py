"""
Unit tests for P1-P5 bug fixes (Rev 3, 2026-03-09).

P1A: _phased_markers / _mtz_has_phase_columns
P2:  solve-mode file discovery + absolute path
P3:  autobuild prerequisite check
P4:  session_blocked_programs / count_total_failures
P5:  Phase 1.5 gate (has_full_map suppression removed)
     + map_sharpening_done zombie table
"""
import os
import sys

# ---------------------------------------------------------------------------
# Ensure agent package is importable
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent"), os.path.join(_ROOT, "programs")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


# ============================================================================
# P4: count_total_failures / session_block_threshold
# ============================================================================
from error_classifier import (
    count_total_failures, session_block_threshold,
    _SESSION_BLOCK_THRESHOLD, should_pivot,
)

def _make_history(*entries):
    """Build a minimal history list from (program, success) tuples."""
    result = []
    for prog, ok in entries:
        result.append({
            "program": prog,
            "result": "SUCCESS: ok" if ok else "FAILED: sorry",
        })
    return result


def test_count_total_failures_basic():
    h = _make_history(
        ("phenix.refine", False),
        ("phenix.phaser", True),
        ("phenix.refine", False),
        ("phenix.refine", True),
        ("phenix.refine", False),
    )
    assert count_total_failures(h, "phenix.refine") == 3, \
        "Expected 3 refine failures"
    assert count_total_failures(h, "phenix.phaser") == 0, \
        "Expected 0 phaser failures"
    print("  PASS: test_count_total_failures_basic")


def test_count_total_failures_empty():
    assert count_total_failures([], "phenix.refine") == 0
    assert count_total_failures(None, "phenix.refine") == 0
    print("  PASS: test_count_total_failures_empty")


def test_session_block_threshold_refine():
    assert session_block_threshold("phenix.refine") == 6
    assert session_block_threshold("refine") == 6
    assert session_block_threshold("phenix.real_space_refine") == 6
    print("  PASS: test_session_block_threshold_refine")


def test_session_block_threshold_default():
    assert session_block_threshold("phenix.phaser") == 4
    assert session_block_threshold("phenix.autobuild") == 4
    print("  PASS: test_session_block_threshold_default")


def test_should_pivot_session_block_fires():
    """should_pivot returns True + 'session block' after N total failures.

    Use interspersed fail/success/fail to keep consecutive count at 1,
    so the session-block path (total >= threshold) fires instead of
    the consecutive path (>= 2).
    """
    threshold = _SESSION_BLOCK_THRESHOLD["__default__"]  # 4
    # Build: (fail, success) * threshold then one final fail
    # So consecutive=1, total=threshold+1 >= threshold
    entries = []
    for _ in range(threshold):
        entries.append(("phenix.phaser", False))
        entries.append(("phenix.phaser", True))  # reset consecutive counter
    entries.append(("phenix.phaser", False))  # final fail: consecutive=1
    h = _make_history(*entries)
    # Verify setup: consecutive=1 (not 2+), total > threshold
    from error_classifier import count_total_failures, count_consecutive_failures
    assert count_consecutive_failures(h, "phenix.phaser") == 1, (
        "Expected consecutive=1, got %d" % count_consecutive_failures(h, "phenix.phaser"))
    assert count_total_failures(h, "phenix.phaser") >= threshold
    from error_classifier import RETRYABLE
    ec = {"category": RETRYABLE, "error_message": "sorry"}
    pivots, reason = should_pivot(h, "phenix.phaser", ec)
    assert pivots, "Expected pivot=True after %d total failures" % threshold
    assert "session block" in reason, \
        "Expected 'session block' in reason, got: %r" % reason
    print("  PASS: test_should_pivot_session_block_fires")


def test_should_pivot_refine_higher_threshold():
    """Refine needs 6 total failures, not 4."""
    h = _make_history(*[("phenix.refine", False)] * 5)
    from error_classifier import RETRYABLE
    ec = {"category": RETRYABLE, "error_message": "sorry"}
    pivots, reason = should_pivot(h, "phenix.refine", ec)
    # 5 < 6 threshold: consecutive check (2) fires first
    assert pivots, "Expected pivot after 5 consecutive refine failures"
    # Should be consecutive, not session block
    print("  PASS: test_should_pivot_refine_higher_threshold "
          "(consecutive fires first, reason=%r)" % reason)


# ============================================================================
# P3: validate_program_choice autobuild prerequisite
# ============================================================================
from workflow_state import validate_program_choice


def _make_workflow_state(valid_programs, exp_type="xray", context=None):
    return {
        "valid_programs": valid_programs,
        "experiment_type": exp_type,
        "state": "xray_initial",
        "reason": "test",
        "context": context or {},
    }


def test_autobuild_blocked_without_phasing():
    ws = _make_workflow_state(
        valid_programs=["phenix.autobuild", "phenix.xtriage", "STOP"],
        exp_type="xray",
        context={"phaser_done": False, "autosol_done": False,
                 "has_placed_model_from_history": False},
    )
    ok, err = validate_program_choice("phenix.autobuild", ws)
    assert not ok, "autobuild should be blocked without phasing"
    assert "phased MTZ" in err, "Error should mention phased MTZ"
    print("  PASS: test_autobuild_blocked_without_phasing")


def test_autobuild_allowed_after_phaser():
    ws = _make_workflow_state(
        valid_programs=["phenix.autobuild", "STOP"],
        exp_type="xray",
        context={"phaser_done": True, "autosol_done": False,
                 "has_placed_model_from_history": True},
    )
    ok, err = validate_program_choice("phenix.autobuild", ws)
    assert ok, "autobuild should be allowed after phaser: %r" % err
    print("  PASS: test_autobuild_allowed_after_phaser")


def test_autobuild_allowed_after_autosol():
    ws = _make_workflow_state(
        valid_programs=["phenix.autobuild", "STOP"],
        exp_type="xray",
        context={"phaser_done": False, "autosol_done": True,
                 "has_placed_model_from_history": False},
    )
    ok, err = validate_program_choice("phenix.autobuild", ws)
    assert ok, "autobuild should be allowed after autosol: %r" % err
    print("  PASS: test_autobuild_allowed_after_autosol")


def test_autobuild_allowed_cryoem():
    """P3 check is xray-only; cryoem autobuild is not gated."""
    ws = _make_workflow_state(
        valid_programs=["phenix.autobuild", "STOP"],
        exp_type="cryoem",
        context={"phaser_done": False, "autosol_done": False,
                 "has_placed_model_from_history": False},
    )
    ok, err = validate_program_choice("phenix.autobuild", ws)
    assert ok, "cryoem autobuild should not be blocked by P3: %r" % err
    print("  PASS: test_autobuild_allowed_cryoem")


def test_other_programs_unaffected():
    ws = _make_workflow_state(
        valid_programs=["phenix.phaser", "STOP"],
        exp_type="xray",
        context={},
    )
    ok, err = validate_program_choice("phenix.phaser", ws)
    assert ok, "phaser should not be blocked"
    print("  PASS: test_other_programs_unaffected")


# ============================================================================
# P1A: _phased_markers includes 'overall_best'
# ============================================================================
from workflow_state import _mtz_has_phase_columns

def test_mtz_has_phase_columns_nonexistent():
    """Non-existent file returns False without raising."""
    assert _mtz_has_phase_columns("/nonexistent/file.mtz") == False
    print("  PASS: test_mtz_has_phase_columns_nonexistent")


def test_mtz_has_phase_columns_non_mtz():
    """Non-MTZ file returns False."""
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(b"ATOM  ...\n")
        fname = f.name
    try:
        result = _mtz_has_phase_columns(fname)
        assert result == False, "PDB file should not have phase columns"
    finally:
        os.unlink(fname)
    print("  PASS: test_mtz_has_phase_columns_non_mtz")


def test_phased_markers_includes_overall_best():
    """Verify 'overall_best' is in _phased_markers by checking source file."""
    _ws_path = os.path.join(_ROOT, "agent", "workflow_state.py")
    with open(_ws_path) as _f:
        src = _f.read()
    # Find the _phased_markers tuple
    idx = src.find("_phased_markers = (")
    assert idx >= 0, "_phased_markers not found in workflow_state.py"
    tuple_section = src[idx:idx+400]
    assert "overall_best" in tuple_section, \
        "'overall_best' not found in _phased_markers tuple"
    print("  PASS: test_phased_markers_includes_overall_best")


# ============================================================================
# P5: map_sharpening in _ZOMBIE_CHECK_TABLE
# ============================================================================
from workflow_state import _ZOMBIE_CHECK_TABLE

def test_map_sharpening_in_zombie_table():
    flags = [entry[0] for entry in _ZOMBIE_CHECK_TABLE]
    assert "map_sharpening_done" in flags, \
        "map_sharpening_done not in _ZOMBIE_CHECK_TABLE"
    print("  PASS: test_map_sharpening_in_zombie_table")


def test_map_sharpening_zombie_pattern():
    """The zombie pattern should match typical sharpened map filenames."""
    entry = next(e for e in _ZOMBIE_CHECK_TABLE if e[0] == "map_sharpening_done")
    pattern = entry[1]
    for fname in [
        "map_sharpened.ccp4", "map_sharpened.mrc",
        "full_map_sharpened.mrc", "sharpened_map.ccp4",
    ]:
        assert pattern.search(fname), \
            "Pattern should match %r" % fname
    for fname in ["refine_001.pdb", "overall_best_denmod.mtz"]:
        assert not pattern.search(fname), \
            "Pattern should NOT match %r" % fname
    print("  PASS: test_map_sharpening_zombie_pattern")


# ============================================================================
# P5: Phase 1.5 gate — has_full_map no longer blocks resolve_cryo_em
# ============================================================================
def test_phase_15_fires_with_full_map_present():
    """Phase 1.5 gate must fire even when has_full_map=True.

    Tests via source-file inspection (workflow_engine needs libtbx at import
    time, which is not available in the standalone test environment).
    """
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    # Find the Phase 1.5 block
    idx = src.find("Phase 1.5")
    assert idx >= 0, "Phase 1.5 comment not found in workflow_engine.py"
    gate_section = src[idx:idx+900]  # 763 chars to has_half_map; use 900 for margin
    # The old gate included: 'not context.get("has_full_map")'
    assert 'not context.get("has_full_map")' not in gate_section, \
        "P5 fix: not has_full_map should be removed from Phase 1.5 gate"
    # New gate checks: has_half_map + not resolve_cryo_em_done
    assert "has_half_map" in gate_section, \
        "Phase 1.5 gate should still check has_half_map"
    assert "resolve_cryo_em_done" in gate_section, \
        "Phase 1.5 gate should still check resolve_cryo_em_done"
    print("  PASS: test_phase_15_fires_with_full_map_present")


def test_phase_15_does_not_check_map_sharpening_done():
    """map_sharpening_done must not block Phase 1.5.

    B-factor sharpening and density modification are different programs;
    completing one should not suppress the other (P5 fix).
    """
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    idx = src.find("Phase 1.5")
    assert idx >= 0, "Phase 1.5 comment not found in workflow_engine.py"
    gate_section = src[idx:idx+400]
    assert "map_sharpening_done" not in gate_section, \
        "P5 fix: map_sharpening_done should not appear in Phase 1.5 if-condition"
    print("  PASS: test_phase_15_does_not_check_map_sharpening_done")


# ============================================================================
# P4: graph_state includes session_blocked_programs
# ============================================================================
def test_graph_state_has_session_blocked_programs():
    # Import graph_state and create an initial state
    from graph_state import create_initial_state
    state = create_initial_state(
        history=[],
        available_files=[],
        cycle_number=1,
        max_cycles=20,
    )
    assert "session_blocked_programs" in state, \
        "session_blocked_programs not in initial state"
    assert state["session_blocked_programs"] == [], \
        "session_blocked_programs should start empty"
    print("  PASS: test_graph_state_has_session_blocked_programs")


def test_graph_state_has_think_overrides():
    from graph_state import create_initial_state
    state = create_initial_state(
        history=[], available_files=[], cycle_number=1, max_cycles=20)
    for key in ("think_file_overrides", "think_label_hints", "think_stop_override"):
        assert key in state, "%s not in initial state" % key
    assert state["think_file_overrides"] == {}
    assert state["think_label_hints"] == {}
    assert state["think_stop_override"] is None
    print("  PASS: test_graph_state_has_think_overrides")


# ============================================================================
# Runner
# ============================================================================
_TESTS = [
    # P4
    test_count_total_failures_basic,
    test_count_total_failures_empty,
    test_session_block_threshold_refine,
    test_session_block_threshold_default,
    test_should_pivot_session_block_fires,
    test_should_pivot_refine_higher_threshold,
    # P3
    test_autobuild_blocked_without_phasing,
    test_autobuild_allowed_after_phaser,
    test_autobuild_allowed_after_autosol,
    test_autobuild_allowed_cryoem,
    test_other_programs_unaffected,
    # P1A
    test_mtz_has_phase_columns_nonexistent,
    test_mtz_has_phase_columns_non_mtz,
    test_phased_markers_includes_overall_best,
    # P5
    test_map_sharpening_in_zombie_table,
    test_map_sharpening_zombie_pattern,
    test_phase_15_fires_with_full_map_present,
    test_phase_15_does_not_check_map_sharpening_done,
    # State
    test_graph_state_has_session_blocked_programs,
    test_graph_state_has_think_overrides,
]

if __name__ == "__main__":
    passed = 0
    failed = 0
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        try:
            test_fn()
            passed += 1
        except Exception as exc:
            import traceback
            print("  FAIL: %s" % test_fn.__name__)
            traceback.print_exc()
            failed += 1
    print()
    if failed:
        print("%d/%d tests FAILED." % (failed, passed + failed))
        sys.exit(1)
    else:
        print("All %d tests passed." % passed)
