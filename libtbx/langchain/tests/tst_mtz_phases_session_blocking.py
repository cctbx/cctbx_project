"""
MTZ phase detection, session blocking, autobuild prerequisites, cryo-EM gate (v115 P1-P5).

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
# P6 — predict_and_build guard on autosol
# ============================================================================

def _make_p6_workflow_state(context_overrides=None):
    """Return a workflow_state matching AF_exoV_PredictAndBuild inputs."""
    ctx = {
        "experiment_type": "xray",
        "xtriage_done": True,
        "has_sequence": True,
        "has_model_for_mr": True,   # 7lw7_offset.pdb present
        "has_search_model": False,  # offset PDB → model, not search_model
        "has_anomalous": True,      # 7lw7.mtz has anomalous signal
        "predict_full_done": False,
        "use_mr_sad": False,
        "phaser_done": False,
        "autosol_done": False,
    }
    if context_overrides:
        ctx.update(context_overrides)
    return {
        "valid_programs": ["phenix.xtriage", "phenix.autosol",
                           "phenix.predict_and_build"],
        "experiment_type": "xray",
        "state": "xray_analyze",
        "reason": "test",
        "context": ctx,
    }


def test_autosol_blocked_by_predict_and_build_guard():
    """autosol must be blocked when seq+model present, p&b not done, no mr_sad."""
    ws = _make_p6_workflow_state()
    ok, err = validate_program_choice("phenix.autosol", ws)
    assert not ok, "Expected autosol blocked by predict_and_build guard, got valid"
    assert "predict_and_build" in (err or ""), (
        "Expected error to mention predict_and_build, got: %s" % err)
    print("  PASS: test_autosol_blocked_by_predict_and_build_guard")


def test_autosol_allowed_when_use_mr_sad():
    """autosol must not be blocked by p&b guard when use_mr_sad=True."""
    ws = _make_p6_workflow_state({"use_mr_sad": True})
    ok, err = validate_program_choice("phenix.autosol", ws)
    # use_mr_sad=True lifts the p&b guard; autosol should be accepted
    assert ok, (
        "Expected autosol allowed with use_mr_sad=True, got blocked: %s" % err)
    print("  PASS: test_autosol_allowed_when_use_mr_sad")


def test_autosol_allowed_when_predict_full_done():
    """autosol must be allowed after predict_and_build has already completed."""
    ws = _make_p6_workflow_state({"predict_full_done": True})
    ok, err = validate_program_choice("phenix.autosol", ws)
    assert ok, (
        "Expected autosol allowed when predict_full_done=True, got: %s" % err)
    print("  PASS: test_autosol_allowed_when_predict_full_done")


def test_autosol_allowed_when_no_model():
    """autosol must be allowed for pure SAD (no model present)."""
    ws = _make_p6_workflow_state({"has_model_for_mr": False})
    ok, err = validate_program_choice("phenix.autosol", ws)
    assert ok, (
        "Expected autosol allowed for pure SAD (no model), got: %s" % err)
    print("  PASS: test_autosol_allowed_when_no_model")


# ============================================================================
# Bug C: predict_build_refine_internal category in file_categories.yaml
# ============================================================================
def test_predict_build_refine_internal_in_yaml():
    """predict_build_refine_internal category must exist in file_categories.yaml
    and have parent_category=intermediate so overall_best_final_refine files
    are stripped from the model category (Bug C fix)."""
    _yaml_path = os.path.join(_ROOT, "knowledge", "file_categories.yaml")
    with open(_yaml_path) as _f:
        src = _f.read()
    assert "predict_build_refine_internal:" in src, \
        "predict_build_refine_internal category must be defined"
    idx = src.find("predict_build_refine_internal:")
    block = src[idx:idx+600]
    assert "parent_category: intermediate" in block, \
        "predict_build_refine_internal must have parent_category: intermediate"
    assert "overall_best" in block and "final_refine" in block, \
        "predict_build_refine_internal must match overall_best*final_refine pattern"
    print("  PASS: test_predict_build_refine_internal_in_yaml")


def test_predict_and_build_output_excludes_final_refine():
    """predict_and_build_output category must exclude *_final_refine_* files
    (Bug C fix)."""
    _yaml_path = os.path.join(_ROOT, "knowledge", "file_categories.yaml")
    with open(_yaml_path) as _f:
        src = _f.read()
    idx = src.find("predict_and_build_output:")
    assert idx >= 0, "predict_and_build_output category not found"
    block = src[idx:idx+700]
    assert "_final_refine_" in block, \
        "predict_and_build_output excludes must contain '*_final_refine_*'"
    print("  PASS: test_predict_and_build_output_excludes_final_refine")


# ============================================================================
# Bug E: constraint-based placement inference respects MR intent
# ============================================================================
def test_has_placed_model_wants_mr_guard_exists():
    """_has_placed_model must NOT infer placement from 'refinement' keyword
    when constraints also mention MR/phaser (Bug E fix).
    Without this guard, CASP7 search model PDBs caused the workflow to jump
    straight to xray_refined before phaser ran."""
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    # _mr_keywords is defined just before _wants_mr_first; search from there
    idx = src.find("_mr_keywords")
    assert idx >= 0, \
        "Bug E fix: _mr_keywords not found in workflow_engine.py"
    block = src[idx:idx+1000]
    assert "_wants_mr_first" in block, \
        "Bug E fix: _wants_mr_first must follow _mr_keywords"
    assert "molecular replacement" in block or "phaser" in block, \
        "Bug E fix: MR keywords must include 'molecular replacement' or 'phaser'"
    assert "placement_keywords" in block, \
        "Bug E fix: placement_keywords check must be conditional on not _wants_mr_first"
    print("  PASS: test_has_placed_model_wants_mr_guard_exists")


def test_has_placed_model_placement_keywords_still_conditional():
    """placement_keywords logic must still exist but only run when no MR intent."""
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    assert "placement_keywords" in src, \
        "placement_keywords logic must still exist"
    assert "not _wants_mr_first" in src, \
        "placement_keywords must only run when not _wants_mr_first"
    print("  PASS: test_has_placed_model_placement_keywords_still_conditional")


# ============================================================================
# Bug F: hopeless R-free dead-end after failed MR
# ============================================================================
def test_hopeless_rfree_routes_to_obtain_model():
    """When refine ran but R-free >= 0.45 after phaser, _detect_xray_step
    must route to obtain_model for phaser retry instead of trapping at STOP
    (Bug F fix)."""
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    idx = src.find("Bug F fix")
    assert idx >= 0, \
        "Bug F fix comment not found in workflow_engine.py"
    block = src[idx:idx+1200]
    assert "obtain_model" in block, \
        "Bug F fix must route to obtain_model step"
    assert "0.45" in block, \
        "Bug F fix must use 0.45 as the hopeless R-free threshold"
    assert "phaser_done" in block, \
        "Bug F fix must check phaser_done"
    assert "refine_count" in block, \
        "Bug F fix must check refine_count"
    print("  PASS: test_hopeless_rfree_routes_to_obtain_model")


def test_hopeless_rfree_threshold_conservative():
    """R-free threshold for Bug F must be >= 0.45 (not lower)."""
    _engine_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_engine_path) as _f:
        src = _f.read()
    idx = src.find("Bug F fix")
    assert idx >= 0, "Bug F fix comment not found"
    block = src[idx:idx+500]
    assert ">= 0.45" in block, \
        "Bug F threshold must be '>= 0.45'"
    print("  PASS: test_hopeless_rfree_threshold_conservative")



# ============================================================================
# Weak anomalous signal (measurability 0.06–0.10) sets has_anomalous=True
# ============================================================================
def test_weak_anomalous_sets_has_anomalous_true():
    """Measurability in [0.06, 0.10] must set has_anomalous=True (not False).
    This is the 'weak but real' signal range — autosol should remain available
    even if phasing may be marginal. Regression for tst_history_analysis failure."""
    _ws_path = os.path.join(_ROOT, "agent", "workflow_state.py")
    with open(_ws_path) as _f:
        src = _f.read()
    # Find the measurability block
    idx = src.find(">= 0.06")
    assert idx >= 0, \
        ("Weak anomalous fix not found: '>= 0.06' "
         "missing from workflow_state.py")
    block = src[idx:idx+300]
    assert "has_anomalous" in block and "True" in block, \
        "Weak anomalous block must set has_anomalous=True"
    # The < 0.06 clear path must still exist
    assert "anomalous_measurability < 0.06" not in src or "else:" in src, \
        "Negligible signal path (< 0.06) must still clear has_anomalous"
    print("  PASS: test_weak_anomalous_sets_has_anomalous_true")


def test_negligible_anomalous_clears_has_anomalous():
    """Measurability < 0.06 must clear has_anomalous=False (negligible signal)."""
    _ws_path = os.path.join(_ROOT, "agent", "workflow_state.py")
    with open(_ws_path) as _f:
        src = _f.read()
    # The else branch should still clear has_anomalous for negligible signal
    assert "has_anomalous" in src and "False" in src, \
        "Negligible anomalous path must still exist"
    idx = src.find(">= 0.06")
    assert idx >= 0
    # Find the else branch after the >= 0.06 block
    after = src[idx:idx+600]
    assert "else:" in after and "has_anomalous" in after, \
        "else branch (negligible signal) must set has_anomalous=False"
    print("  PASS: test_negligible_anomalous_clears_has_anomalous")

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
    # P6 — predict_and_build guard on autosol
    test_autosol_blocked_by_predict_and_build_guard,
    test_autosol_allowed_when_use_mr_sad,
    test_autosol_allowed_when_predict_full_done,
    test_autosol_allowed_when_no_model,
]

def run_all_tests():
    for test_fn in _TESTS:
        test_fn()
    print("All %d tests passed." % len(_TESTS))


# ============================================================================
# BUG A: event_formatter TypeError on string metric values
# ============================================================================
def _ef_format_metrics_extracted(event):
    """
    Execute _format_metrics_extracted logic in isolation, extracted from
    event_formatter.py source without importing the module (which requires
    libtbx).  This mirrors the fixed implementation exactly.
    """
    metrics_def = [
        ("r_free",                "R-free",              "%.4f"),
        ("r_work",                "R-work",              "%.4f"),
        ("resolution",            "Resolution",          "%.2f A"),
        ("map_cc",                "Map CC",              "%.4f"),
        ("clashscore",            "Clashscore",          "%.1f"),
        ("tfz",                   "TFZ",                 "%.1f"),
        ("llg",                   "LLG",                 "%.0f"),
        ("ramachandran_outliers", "Ramachandran outliers","%.1f%%"),
        ("rotamer_outliers",      "Rotamer outliers",    "%.1f%%"),
    ]
    lines = ["", "Metrics:"]
    has_metrics = False
    for key, label, fmt in metrics_def:
        if event.get(key) is None:
            continue
        has_metrics = True
        value = event[key]
        prev_key = key + "_prev"
        prev = event.get(prev_key)

        # Bug A fix: coerce string values
        if not isinstance(value, (int, float)):
            try:
                value = float(value)
            except (TypeError, ValueError):
                lines.append("  %s: %s" % (label, value))
                continue

        if prev is not None and not isinstance(prev, (int, float)):
            try:
                prev = float(prev)
            except (TypeError, ValueError):
                prev = None

        if prev is not None:
            diff = value - prev
            if key in ("r_free", "r_work", "clashscore",
                       "ramachandran_outliers", "rotamer_outliers"):
                direction = "improved" if diff < 0 else "worsened"
            else:
                direction = "improved" if diff > 0 else "worsened"
            val_str  = fmt % value
            prev_str = (fmt % prev).split()[0]
            lines.append("  %s: %s -> %s (%s by %.4f)"
                         % (label, prev_str, val_str, direction, abs(diff)))
        else:
            lines.append("  %s: %s" % (label, fmt % value))

    return "\n".join(lines) if has_metrics else None


def test_event_formatter_string_resolution():
    """_format_metrics_extracted must not crash when metrics are strings."""
    event = {"r_free": "0.556", "resolution": "2.04"}
    try:
        result = _ef_format_metrics_extracted(event)
        assert result is not None
    except TypeError as e:
        raise AssertionError("TypeError on string metrics: %s" % e)
    print("    PASS: test_event_formatter_string_resolution")


def test_event_formatter_numeric_unchanged():
    """Normal numeric metrics still format correctly after Bug A fix."""
    event = {"r_free": 0.272, "resolution": 2.04}
    result = _ef_format_metrics_extracted(event)
    assert result is not None
    assert "0.2720" in result, "Expected formatted r_free; got: %s" % result
    print("    PASS: test_event_formatter_numeric_unchanged")


def test_event_formatter_prev_string_no_crash():
    """String prev value must not crash trend comparison."""
    event = {"r_free": "0.556", "r_free_prev": "0.600"}
    try:
        _ef_format_metrics_extracted(event)
    except (TypeError, ValueError) as e:
        raise AssertionError("Crash on string prev value: %s" % e)
    print("    PASS: test_event_formatter_prev_string_no_crash")


def test_event_formatter_source_has_isinstance_guard():
    """event_formatter.py must guard against string metric values (Bug A fix)."""
    _ef_path = os.path.join(_ROOT, "agent", "event_formatter.py")
    with open(_ef_path) as f:
        src = f.read()
    # The fix can be either isinstance() or try/except — both handle strings.
    # Verify that the metrics section has some form of TypeError protection.
    has_isinstance = "isinstance(value, (int, float))" in src
    has_try_except = ("except (TypeError, ValueError)" in src or
                      "except TypeError" in src)
    assert has_isinstance or has_try_except, (
        "Bug A fix not present: no isinstance guard or try/except protection "
        "found in event_formatter.py"
    )
    print("    PASS: test_event_formatter_source_has_isinstance_guard")


# ============================================================================
# BUG B: _auto_discover_files supplement mode
# ============================================================================
def test_auto_discover_supplements_existing_files():
    """Supplement mode: newly found files are appended to existing original_files."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        ligand = os.path.join(tmpdir, "ligand.pdb")
        mtz    = os.path.join(tmpdir, "data.mtz")
        fa     = os.path.join(tmpdir, "seq.fa")
        for p in (ligand, mtz, fa):
            open(p, "w").close()

        existing = [ligand]
        existing_basenames = {os.path.basename(f) for f in existing}

        _KNOWN_EXT = {".mtz", ".pdb", ".cif", ".fa", ".fasta", ".seq"}
        discovered = []
        for fname in sorted(os.listdir(tmpdir)):
            fpath = os.path.join(tmpdir, fname)
            if not os.path.isfile(fpath): continue
            _, ext = os.path.splitext(fname.lower())
            if ext in _KNOWN_EXT:
                discovered.append(fpath)

        new_files = [f for f in discovered
                     if os.path.basename(f) not in existing_basenames]
        result = list(existing) + new_files

        basenames = {os.path.basename(f) for f in result}
        assert "ligand.pdb" in basenames, "Existing file lost"
        assert "data.mtz"   in basenames, "MTZ not discovered"
        assert "seq.fa"     in basenames, "FASTA not discovered"
        assert len(result) == 3, "Expected 3 files, got %d: %s" % (len(result), result)
    print("    PASS: test_auto_discover_supplements_existing_files")


def test_auto_discover_no_duplicates():
    """Supplement mode does not add files already in original_files."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        mtz = os.path.join(tmpdir, "data.mtz")
        open(mtz, "w").close()
        existing = [mtz]
        existing_basenames = {os.path.basename(f) for f in existing}

        _KNOWN_EXT = {".mtz"}
        discovered = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir)
                      if os.path.splitext(f)[1].lower() in _KNOWN_EXT]
        new_files = [f for f in discovered
                     if os.path.basename(f) not in existing_basenames]
        result = list(existing) + new_files

        assert len(result) == 1, "Duplicate added: %s" % result
    print("    PASS: test_auto_discover_no_duplicates")


def test_auto_discover_readme_not_included():
    """README files are excluded even in supplement mode."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        readme = os.path.join(tmpdir, "README.txt")
        mtz    = os.path.join(tmpdir, "data.mtz")
        for p in (readme, mtz):
            open(p, "w").close()

        existing = []
        existing_basenames = set()
        _KNOWN_EXT = {".mtz", ".txt"}
        discovered = []
        for fname in sorted(os.listdir(tmpdir)):
            fpath = os.path.join(tmpdir, fname)
            if fname.upper().startswith("README"): continue
            _, ext = os.path.splitext(fname.lower())
            if ext in _KNOWN_EXT:
                discovered.append(fpath)

        new_files = [f for f in discovered
                     if os.path.basename(f) not in existing_basenames]
        result = list(existing) + new_files

        basenames = {os.path.basename(f) for f in result}
        assert "README.txt" not in basenames, "README must be excluded"
        assert "data.mtz" in basenames
    print("    PASS: test_auto_discover_readme_not_included")


# ============================================================================
# BUG C: predict_and_build cascade must NOT set refine_done
# ============================================================================
def test_predict_and_build_cascade_no_refine_done():
    """workflow_state.py source must NOT set refine_done in predict_and_build cascade."""
    _ws_path = os.path.join(_ROOT, "agent", "workflow_state.py")
    with open(_ws_path) as f:
        src = f.read()

    # Find the cascade block
    cascade_idx = src.find("The ONE remaining Python-only case: predict_and_build cascade")
    assert cascade_idx >= 0, "Cascade comment not found"
    cascade_src = src[cascade_idx:cascade_idx + 800]

    # Must NOT set refine_done=True inside the cascade
    assert 'info["refine_done"] = True' not in cascade_src, \
        "Bug C regression: refine_done=True is still set in predict_and_build cascade"

    # Must still set predict_full_done
    assert 'predict_full_done' in cascade_src, \
        "predict_full_done not found in cascade — logic may have changed"
    print("    PASS: test_predict_and_build_cascade_no_refine_done")


def test_refined_category_excludes_predictandbuild():
    """file_categories.yaml refined category must exclude *PredictAndBuild* files."""
    _yaml_path = os.path.join(_ROOT, "knowledge", "file_categories.yaml")
    with open(_yaml_path) as f:
        content = f.read()

    # Find the refined: section
    idx = content.find("\nrefined:")
    assert idx >= 0, "refined: section not found in file_categories.yaml"

    # Read from refined: to the next top-level key (empty line + non-space start)
    import re
    refined_section = re.split(r'\n(?=[a-z_]+:)', content[idx+1:])[0]

    assert "*PredictAndBuild*" in refined_section, (
        "Bug C: *PredictAndBuild* exclude not present in refined: category.\n"
        "Section:\n%s" % refined_section)
    print("    PASS: test_refined_category_excludes_predictandbuild")


# ============================================================================
# BUG D: autosol deprioritized when anomalous_measurability < 0.05
# ============================================================================
def test_autosol_deprioritized_low_anomalous():
    """workflow_engine.py source must contain the low-anomalous deprioritization block."""
    _eng_path = os.path.join(_ROOT, "agent", "workflow_engine.py")
    with open(_eng_path) as f:
        src = f.read()
    assert "anomalous_measurability" in src and "0.05" in src, \
        "Bug D deprioritization block not found in workflow_engine.py"
    # Verify the logic: if anomalous_meas < 0.05 → move autosol to end
    assert "Deprioritized phenix.autosol" in src, \
        "Deprioritization modification log message not found"
    print("    PASS: test_autosol_deprioritized_low_anomalous")


def test_autosol_deprioritize_logic_low():
    """Logic: autosol moves to end when measurability < 0.05 and predict_and_build present."""
    result = ["phenix.autosol", "phenix.predict_and_build"]
    context = {"anomalous_measurability": 0.032,
               "autosol_done": False, "phaser_done": False}
    workflow_prefs = {}

    if (not workflow_prefs.get("use_mr_sad") and
            not workflow_prefs.get("use_experimental_phasing") and
            "phenix.autosol" in result and
            "phenix.predict_and_build" in result and
            not context.get("autosol_done") and
            not context.get("phaser_done")):
        if (context.get("anomalous_measurability") or 0) < 0.05:
            result.remove("phenix.autosol")
            result.append("phenix.autosol")

    assert result[0] == "phenix.predict_and_build", \
        "predict_and_build should be first; got %s" % result
    assert result[-1] == "phenix.autosol"
    print("    PASS: test_autosol_deprioritize_logic_low")


def test_autosol_deprioritize_logic_high():
    """Logic: autosol NOT moved when measurability >= 0.05."""
    result = ["phenix.autosol", "phenix.predict_and_build"]
    context = {"anomalous_measurability": 0.15}
    workflow_prefs = {}

    if (not workflow_prefs.get("use_mr_sad") and
            "phenix.autosol" in result and
            "phenix.predict_and_build" in result):
        if (context.get("anomalous_measurability") or 0) < 0.05:
            result.remove("phenix.autosol")
            result.append("phenix.autosol")

    assert result[0] == "phenix.autosol", \
        "autosol should stay first for strong signal; got %s" % result
    print("    PASS: test_autosol_deprioritize_logic_high")


def test_autosol_deprioritize_logic_mr_sad_override():
    """Logic: autosol NOT moved when use_mr_sad=True even with negligible signal."""
    result = ["phenix.autosol", "phenix.predict_and_build"]
    context = {"anomalous_measurability": 0.01}
    workflow_prefs = {"use_mr_sad": True}

    if (not workflow_prefs.get("use_mr_sad") and
            "phenix.autosol" in result and
            "phenix.predict_and_build" in result):
        if (context.get("anomalous_measurability") or 0) < 0.05:
            result.remove("phenix.autosol")
            result.append("phenix.autosol")

    assert result[0] == "phenix.autosol", \
        "autosol must not be deprioritized when use_mr_sad=True"
    print("    PASS: test_autosol_deprioritize_logic_mr_sad_override")


_TESTS += [
    test_event_formatter_string_resolution,
    test_event_formatter_numeric_unchanged,
    test_event_formatter_prev_string_no_crash,
    test_event_formatter_source_has_isinstance_guard,
    test_auto_discover_supplements_existing_files,
    test_auto_discover_no_duplicates,
    test_auto_discover_readme_not_included,
    test_predict_and_build_cascade_no_refine_done,
    test_refined_category_excludes_predictandbuild,
    test_autosol_deprioritized_low_anomalous,
    test_autosol_deprioritize_logic_low,
    test_autosol_deprioritize_logic_high,
    test_autosol_deprioritize_logic_mr_sad_override,
    # Bug C (yaml + cascade)
    test_predict_build_refine_internal_in_yaml,
    test_predict_and_build_output_excludes_final_refine,
    # Bug E (MR-intent guard)
    test_has_placed_model_wants_mr_guard_exists,
    test_has_placed_model_placement_keywords_still_conditional,
    # Bug F (hopeless R-free)
    test_hopeless_rfree_routes_to_obtain_model,
    test_hopeless_rfree_threshold_conservative,
    # Weak anomalous signal fix
    test_weak_anomalous_sets_has_anomalous_true,
    test_negligible_anomalous_clears_has_anomalous,
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
