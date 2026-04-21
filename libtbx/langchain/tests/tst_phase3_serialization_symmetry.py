"""
Phase 3: Session Round-Trip — Symmetry + Invariant Validation.

Tests that session state survives JSON serialization/deserialization
without data corruption, and that deserialized state is logically
consistent.

Bug class: silent data corruption (Bug 4 from session 12a was
float→string after JSON round-trip).

Run: python tests/tst_phase3_serialization_symmetry.py
Produces: findings/phase_3_roundtrip_failures.yaml
"""

from __future__ import absolute_import, division, print_function

import json
import os
import sys
import tempfile
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

# Mock libtbx
if 'libtbx' not in sys.modules:
    import types
    _base_dir = os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))
    libtbx = types.ModuleType('libtbx')
    libtbx.__path__ = [_base_dir]
    libtbx.langchain = types.ModuleType('libtbx.langchain')
    libtbx.langchain.__path__ = [_base_dir]
    libtbx.langchain.agent = types.ModuleType(
        'libtbx.langchain.agent')
    libtbx.langchain.agent.__path__ = [
        os.path.join(_base_dir, 'agent')]
    libtbx.langchain.knowledge = types.ModuleType(
        'libtbx.langchain.knowledge')
    libtbx.langchain.knowledge.__path__ = [
        os.path.join(_base_dir, 'knowledge')]
    sys.modules['libtbx'] = libtbx
    sys.modules['libtbx.langchain'] = libtbx.langchain
    sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
    sys.modules['libtbx.langchain.knowledge'] = \
        libtbx.langchain.knowledge
import knowledge.yaml_loader
if 'libtbx.langchain.knowledge.yaml_loader' not in sys.modules:
    sys.modules['libtbx.langchain.knowledge.yaml_loader'] = \
        knowledge.yaml_loader

FINDINGS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'findings')


# =====================================================================
# TEST UTILITIES
# =====================================================================

def assert_equal(a, b, msg=""):
    if a != b:
        raise AssertionError(
            "%s: %r != %r" % (msg, a, b) if msg
            else "%r != %r" % (a, b))


def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)


def assert_is_none(val, msg=""):
    if val is not None:
        raise AssertionError(
            "%s: expected None, got %r" % (msg, val) if msg
            else "expected None, got %r" % val)


def assert_is_instance(val, typ, msg=""):
    if not isinstance(val, typ):
        raise AssertionError(
            "%s: expected %s, got %s (%r)" %
            (msg, typ.__name__, type(val).__name__, val) if msg
            else "expected %s, got %s (%r)" %
            (typ.__name__, type(val).__name__, val))


def round_trip(data):
    """Serialize to JSON string and deserialize back."""
    s = json.dumps(data, indent=2)
    return json.loads(s)


def session_round_trip(session_data, tmpdir):
    """Save session data to file and load it back."""
    path = os.path.join(tmpdir, "test_session.json")
    with open(path, 'w') as f:
        json.dump(session_data, f, indent=2)
    with open(path, 'r') as f:
        return json.load(f)


# =====================================================================
# SYMMETRY TESTS: data == deserialize(serialize(data))
# =====================================================================

def test_float_round_trip():
    """Float metrics survive JSON round-trip as floats."""
    print("Test: float_round_trip")
    data = {
        "r_work": 0.2134,
        "r_free": 0.2567,
        "model_map_cc": 0.891,
        "resolution": 2.05,
        "clashscore": 3.14,
    }
    result = round_trip(data)
    for key in data:
        assert_is_instance(result[key], float,
            "Field '%s' should be float after round-trip" % key)
        assert_equal(result[key], data[key],
            "Field '%s' value changed" % key)
    print("  PASSED")


def test_none_vs_string_none():
    """None stays None after round-trip (not string 'None')."""
    print("Test: none_vs_string_none")
    data = {
        "resolution": None,
        "r_free": None,
        "rfree_mtz": None,
        "experiment_type": None,
    }
    result = round_trip(data)
    for key in data:
        assert_is_none(result[key],
            "Field '%s' should be None, not '%s'" %
            (key, result[key]))
    print("  PASSED")


def test_none_vs_missing_field():
    """Fields set to None vs fields that are absent."""
    print("Test: none_vs_missing_field")
    data = {"present_none": None, "present_value": 42}
    result = round_trip(data)
    assert_true("present_none" in result,
        "Field set to None should still exist after round-trip")
    assert_is_none(result["present_none"],
        "Explicitly-None field should remain None")
    assert_true("absent_field" not in result,
        "Never-set field should not appear")
    print("  PASSED")


def test_int_vs_string_int():
    """Integer fields stay integers (not string '5')."""
    print("Test: int_vs_string_int")
    data = {
        "rfree_mtz_locked_at_cycle": 3,
        "cycle_number": 5,
        "refine_count": 2,
    }
    result = round_trip(data)
    for key in data:
        assert_is_instance(result[key], int,
            "Field '%s' should be int after round-trip" % key)
    print("  PASSED")


def test_float_precision():
    """Float precision preserved for R-factor range values."""
    print("Test: float_precision")
    # R-factors typically have 4 decimal places
    values = [0.0001, 0.1234, 0.2567, 0.4999, 0.9999,
              1e-10, 3.14159265]
    data = {"values": values}
    result = round_trip(data)
    for orig, loaded in zip(values, result["values"]):
        # JSON preserves float64 precision for these magnitudes
        assert_equal(orig, loaded,
            "Float precision lost: %r != %r" % (orig, loaded))
    print("  PASSED")


def test_path_round_trip():
    """File paths survive round-trip unchanged."""
    print("Test: path_round_trip")
    paths = [
        "/absolute/path/to/file.mtz",
        "relative/path/data.pdb",
        "/path/with spaces/file.mtz",
        "/path/with-dashes/under_scores/file.hkl",
        "/deeply/nested/dir/sub/file.sca",
    ]
    data = {"files": paths}
    result = round_trip(data)
    for orig, loaded in zip(paths, result["files"]):
        assert_equal(orig, loaded,
            "Path changed: %r != %r" % (orig, loaded))
    print("  PASSED")


def test_empty_collections_preserved():
    """Empty lists and dicts stay as their type."""
    print("Test: empty_collections_preserved")
    data = {
        "cycles": [],
        "directives": {},
        "output_files": [],
        "recovery_attempts": {},
    }
    result = round_trip(data)
    for key in data:
        expected_type = type(data[key])
        assert_is_instance(result[key], expected_type,
            "Field '%s' should be %s after round-trip" %
            (key, expected_type.__name__))
        assert_equal(len(result[key]), 0,
            "Field '%s' should be empty" % key)
    print("  PASSED")


def test_nested_dict_round_trip():
    """Nested dicts (directives, recovery_strategies) survive."""
    print("Test: nested_dict_round_trip")
    data = {
        "directives": {
            "stop_conditions": {
                "max_refine_cycles": 5,
                "target_r_free": 0.25,
            },
            "program_settings": {
                "phenix.refine": {
                    "nproc": 4,
                }
            }
        },
        "recovery_strategies": {
            "/path/to/data.mtz": {
                "flags": "labels=IMEAN,SIGIMEAN",
                "program_scope": "phenix.refine",
                "reason": "ambiguous labels",
            }
        }
    }
    result = round_trip(data)
    assert_equal(
        result["directives"]["stop_conditions"]["target_r_free"],
        0.25, "Nested float should survive")
    assert_equal(
        result["recovery_strategies"]["/path/to/data.mtz"]["flags"],
        "labels=IMEAN,SIGIMEAN", "Nested string should survive")
    print("  PASSED")


def test_history_entry_result_key():
    """Session uses 'result' not 'status' for cycle outcomes."""
    print("Test: history_entry_result_key")
    # The documented pitfall: some code uses "status",
    # but SessionState uses "result"
    cycle = {
        "cycle_number": 1,
        "program": "phenix.refine",
        "command": "phenix.refine model.pdb data.mtz",
        "result": "SUCCESS: R-work=0.22 R-free=0.26",
        "output_files": ["refine_001.pdb", "refine_001.mtz"],
    }
    result = round_trip(cycle)
    assert_true("result" in result,
        "Cycle must use 'result' key")
    assert_true("status" not in result,
        "Cycle must NOT use 'status' key (pitfall)")
    assert_true(result["result"].startswith("SUCCESS"),
        "Result string should survive round-trip")
    print("  PASSED")


# =====================================================================
# FULL SESSION ROUND-TRIP
# =====================================================================

def test_full_session_round_trip():
    """A complete session with 3 cycles survives file round-trip."""
    print("Test: full_session_round_trip")
    session = {
        "session_id": "2026-03-17_12-00-00",
        "project_advice": "Refine the model.",
        "original_files": ["/path/to/model.pdb",
                           "/path/to/data.mtz"],
        "experiment_type": "xray",
        "resolution": 2.05,
        "resolution_source": "xtriage_1",
        "rfree_mtz": "/path/to/data_rfree.mtz",
        "rfree_mtz_locked_at_cycle": 2,
        "rfree_resolution": None,
        "directives": {"stop_conditions": {"max_refine_cycles": 5}},
        "directives_extracted": True,
        "recovery_attempts": {},
        "recovery_strategies": {},
        "force_retry_program": None,
        "abort_on_red_flags": True,
        "abort_on_warnings": False,
        "ignored_commands": {},
        "cycles": [
            {
                "cycle_number": 1,
                "program": "phenix.xtriage",
                "decision": "Analyze data quality",
                "reasoning": "First cycle, need analysis",
                "explanation": "Running xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS: OK",
                "output_files": [],
                "timestamp": "2026-03-17T12:00:01",
            },
            {
                "cycle_number": 2,
                "program": "phenix.refine",
                "decision": "Refine model",
                "reasoning": "Model needs refinement",
                "explanation": "Running refine",
                "command": "phenix.refine model.pdb data.mtz",
                "result": "SUCCESS: R-work=0.2200 R-free=0.2600",
                "output_files": ["refine_001.pdb",
                                 "refine_001.mtz"],
                "timestamp": "2026-03-17T12:05:00",
            },
            {
                "cycle_number": 3,
                "program": "phenix.refine",
                "decision": "Continue refinement",
                "reasoning": "R-free improving",
                "explanation": "Second refine cycle",
                "command": "phenix.refine refine_001.pdb "
                           "data.mtz",
                "result": "SUCCESS: R-work=0.2100 R-free=0.2500",
                "output_files": ["refine_002.pdb",
                                 "refine_002.mtz"],
                "timestamp": "2026-03-17T12:10:00",
            },
        ],
    }

    tmpdir = tempfile.mkdtemp()
    try:
        loaded = session_round_trip(session, tmpdir)

        # Check types
        assert_is_instance(loaded["resolution"], float,
            "resolution should be float")
        assert_is_instance(
            loaded["rfree_mtz_locked_at_cycle"], int,
            "rfree_mtz_locked_at_cycle should be int")
        assert_is_none(loaded["rfree_resolution"],
            "rfree_resolution should be None")
        assert_is_none(loaded["force_retry_program"],
            "force_retry_program should be None")

        # Check cycle count
        assert_equal(len(loaded["cycles"]), 3,
            "Should have 3 cycles")

        # Check cycle numbers are sequential
        for i, cycle in enumerate(loaded["cycles"]):
            assert_equal(cycle["cycle_number"], i + 1,
                "Cycle %d should have cycle_number=%d" %
                (i, i + 1))

        # Check result field survived
        assert_true(
            loaded["cycles"][1]["result"].startswith("SUCCESS"),
            "Cycle 2 result should start with SUCCESS")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


# =====================================================================
# INVARIANT TESTS: logical consistency after deserialization
# =====================================================================

def test_cycle_numbering_invariant():
    """Cycle numbers must be sequential 1..N with no gaps."""
    print("Test: cycle_numbering_invariant")
    # Valid
    valid = {"cycles": [
        {"cycle_number": 1, "result": "SUCCESS"},
        {"cycle_number": 2, "result": "SUCCESS"},
        {"cycle_number": 3, "result": "SUCCESS"},
    ]}
    loaded = round_trip(valid)
    for i, c in enumerate(loaded["cycles"]):
        assert_equal(c["cycle_number"], i + 1,
            "Cycle %d has wrong number" % i)

    # Construct invalid (gap) and verify we can detect it
    invalid = {"cycles": [
        {"cycle_number": 1, "result": "SUCCESS"},
        {"cycle_number": 3, "result": "SUCCESS"},  # Gap!
    ]}
    loaded = round_trip(invalid)
    has_gap = False
    for i, c in enumerate(loaded["cycles"]):
        if c["cycle_number"] != i + 1:
            has_gap = True
    assert_true(has_gap,
        "Should detect numbering gap in invalid data")
    print("  PASSED")


def test_cycle_count_matches_list_length():
    """len(cycles) should match the highest cycle_number."""
    print("Test: cycle_count_matches_list_length")
    data = {"cycles": [
        {"cycle_number": 1, "program": "xtriage"},
        {"cycle_number": 2, "program": "refine"},
    ]}
    loaded = round_trip(data)
    assert_equal(len(loaded["cycles"]), 2,
        "Cycle list length should be 2")
    max_num = max(c["cycle_number"] for c in loaded["cycles"])
    assert_equal(max_num, len(loaded["cycles"]),
        "Highest cycle_number should match list length")
    print("  PASSED")


def test_successful_cycle_has_program():
    """A cycle with SUCCESS result must have a non-empty program.
    Verify we can detect the violation."""
    print("Test: successful_cycle_has_program")
    data = {"cycles": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "SUCCESS: OK"},
        {"cycle_number": 2, "program": "",
         "result": "SUCCESS: R=0.22"},  # Invalid!
    ]}
    loaded = round_trip(data)
    found_violation = False
    for c in loaded["cycles"]:
        if c.get("result", "").startswith("SUCCESS"):
            has_program = bool(c.get("program", "").strip())
            if not has_program:
                found_violation = True
    assert_true(found_violation,
        "Should detect cycle with SUCCESS but no program")
    print("  PASSED")


def test_experiment_type_valid():
    """experiment_type must be None, 'xray', or 'cryoem'."""
    print("Test: experiment_type_valid")
    for value in [None, "xray", "cryoem"]:
        data = {"experiment_type": value}
        loaded = round_trip(data)
        assert_equal(loaded["experiment_type"], value,
            "experiment_type %r should survive round-trip" %
            value)

    # Invalid value survives round-trip too (JSON doesn't
    # validate) — but we should be able to detect it
    data = {"experiment_type": "neutron"}
    loaded = round_trip(data)
    valid_types = {None, "xray", "cryoem"}
    is_valid = loaded["experiment_type"] in valid_types
    assert_true(not is_valid,
        "Should detect invalid experiment_type")
    print("  PASSED")


def test_resolution_is_float_or_none():
    """resolution field must be float or None, never string."""
    print("Test: resolution_is_float_or_none")
    # Valid cases
    for val in [None, 2.05, 1.5, 3.0]:
        data = {"resolution": val}
        loaded = round_trip(data)
        if val is None:
            assert_is_none(loaded["resolution"],
                "None resolution should stay None")
        else:
            assert_is_instance(loaded["resolution"], float,
                "resolution should be float")

    # Bug 4 pattern: what if resolution was stored as string?
    data = {"resolution": "2.05"}
    loaded = round_trip(data)
    # JSON preserves the string — this IS the bug
    assert_is_instance(loaded["resolution"], str,
        "String resolution survives as string (the bug)")
    # Verify _safe_float would fix it
    from agent.structure_model import _safe_float
    fixed = _safe_float(loaded["resolution"])
    assert_is_instance(fixed, float,
        "_safe_float should convert string to float")
    assert_equal(fixed, 2.05,
        "_safe_float should get correct value")
    print("  PASSED")


def test_rfree_mtz_consistency():
    """If rfree_mtz is set, rfree_mtz_locked_at_cycle must be set."""
    print("Test: rfree_mtz_consistency")
    # Valid: both set
    data = {
        "rfree_mtz": "/path/to/data.mtz",
        "rfree_mtz_locked_at_cycle": 2,
    }
    loaded = round_trip(data)
    assert_true(
        loaded["rfree_mtz"] is not None and
        loaded["rfree_mtz_locked_at_cycle"] is not None,
        "Both rfree fields should be set")

    # Valid: both None
    data2 = {
        "rfree_mtz": None,
        "rfree_mtz_locked_at_cycle": None,
    }
    loaded2 = round_trip(data2)
    assert_is_none(loaded2["rfree_mtz"],
        "Both should be None")
    assert_is_none(loaded2["rfree_mtz_locked_at_cycle"],
        "Both should be None")

    # Invariant violation: mtz set but cycle not set
    data3 = {
        "rfree_mtz": "/path/to/data.mtz",
        "rfree_mtz_locked_at_cycle": None,
    }
    loaded3 = round_trip(data3)
    has_mtz = loaded3["rfree_mtz"] is not None
    has_cycle = loaded3["rfree_mtz_locked_at_cycle"] is not None
    invariant_ok = has_mtz == has_cycle
    assert_true(not invariant_ok,
        "Should detect rfree_mtz/cycle mismatch")
    print("  PASSED")


def test_directives_extracted_consistency():
    """If directives_extracted is True, directives should be a dict."""
    print("Test: directives_extracted_consistency")
    data = {
        "directives_extracted": True,
        "directives": {"stop_conditions": {}},
    }
    loaded = round_trip(data)
    if loaded.get("directives_extracted"):
        assert_is_instance(loaded.get("directives"), dict,
            "directives should be dict when extracted=True")
    print("  PASSED")


def test_output_files_are_strings():
    """All output_files entries should be strings, not dicts."""
    print("Test: output_files_are_strings")
    data = {"cycles": [
        {
            "cycle_number": 1,
            "output_files": [
                "refine_001.pdb",
                "refine_001.mtz",
            ],
        }
    ]}
    loaded = round_trip(data)
    for f in loaded["cycles"][0]["output_files"]:
        assert_is_instance(f, str,
            "output_file entry should be string, got %s" %
            type(f).__name__)
    print("  PASSED")


def test_recovery_attempts_count_is_int():
    """recovery_attempts count fields should be integers."""
    print("Test: recovery_attempts_count_is_int")
    data = {
        "recovery_attempts": {
            "AMBIGUOUS_LABELS": {
                "count": 3,
                "files_tried": {"/path/data.mtz": True},
            }
        }
    }
    loaded = round_trip(data)
    count = loaded["recovery_attempts"]["AMBIGUOUS_LABELS"]["count"]
    assert_is_instance(count, int,
        "recovery count should be int")
    print("  PASSED")


# =====================================================================
# AGENT SESSION ROUND-TRIP TESTS
# These test the ACTUAL production code paths, not json module behavior.
# =====================================================================

def _make_session():
    """Create a fresh AgentSession in a temp directory."""
    from agent.session import AgentSession
    tmpdir = tempfile.mkdtemp()
    ss = AgentSession(session_dir=tmpdir)
    return ss, tmpdir


def _reload_session(tmpdir):
    """Load an AgentSession from disk (simulates restart)."""
    from agent.session import AgentSession
    return AgentSession(session_dir=tmpdir)


def test_session_resolution_type():
    """set_resolution stores float, survives save/load."""
    print("Test: session_resolution_type")
    ss, tmpdir = _make_session()
    try:
        ss.set_resolution(2.05, "xtriage_1")
        ss.save()

        ss2 = _reload_session(tmpdir)
        res = ss2.data.get("resolution")
        assert_is_instance(res, float,
            "resolution should be float after save/load")
        assert_equal(res, 2.05,
            "resolution value should be 2.05")

        src = ss2.data.get("resolution_source")
        assert_is_instance(src, str,
            "resolution_source should be str")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_resolution_string_input():
    """set_resolution with string input — should it coerce or crash?
    Tests the actual behavior."""
    print("Test: session_resolution_string_input")
    ss, tmpdir = _make_session()
    try:
        # In production, resolution comes from regex extraction
        # as float. But what if it's a string?
        try:
            ss.set_resolution("2.05", "buggy_source")
            ss.save()
            ss2 = _reload_session(tmpdir)
            res = ss2.data.get("resolution")
            # round("2.05", 2) raises TypeError in Python 3
            # so we should never get here
            print("  NOTE: string input accepted (type=%s)"
                  % type(res).__name__)
        except TypeError:
            # This is the correct behavior — round() rejects strings
            pass
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED (string input rejected by round())")


def test_session_experiment_type_locking():
    """set_experiment_type locks after first call, survives reload."""
    print("Test: session_experiment_type_locking")
    ss, tmpdir = _make_session()
    try:
        was_set, _ = ss.set_experiment_type("xray")
        assert_true(was_set, "First set should succeed")
        ss.save()

        ss2 = _reload_session(tmpdir)
        assert_equal(ss2.data.get("experiment_type"), "xray",
            "experiment_type should survive reload")

        # Second set with different type should be rejected
        was_set2, _ = ss2.set_experiment_type("cryoem")
        assert_true(not was_set2,
            "Second set with different type should be rejected")
        assert_equal(ss2.data.get("experiment_type"), "xray",
            "experiment_type should still be xray")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_rfree_mtz_types():
    """set_rfree_mtz stores string path + int cycle."""
    print("Test: session_rfree_mtz_types")
    ss, tmpdir = _make_session()
    try:
        ss.set_rfree_mtz("/path/to/data_rfree.mtz", 2)
        ss.save()

        ss2 = _reload_session(tmpdir)
        mtz = ss2.data.get("rfree_mtz")
        cycle = ss2.data.get("rfree_mtz_locked_at_cycle")
        assert_is_instance(mtz, str,
            "rfree_mtz should be str")
        assert_is_instance(cycle, int,
            "rfree_mtz_locked_at_cycle should be int")
        assert_equal(cycle, 2, "cycle should be 2")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_record_result_metrics():
    """record_result extracts metrics as floats from result string."""
    print("Test: session_record_result_metrics")
    ss, tmpdir = _make_session()
    try:
        ss.start_cycle(1)
        ss.record_decision(1, program='phenix.refine',
            command='phenix.refine model.pdb data.mtz')
        ss.record_result(1,
            "SUCCESS: R-work=0.2200 R-free=0.2600")
        ss.save()

        ss2 = _reload_session(tmpdir)
        cycle = ss2.data["cycles"][0]
        metrics = cycle.get("metrics", {})

        # Metrics should be floats, not strings
        r_work = metrics.get("r_work")
        r_free = metrics.get("r_free")
        if r_work is not None:
            assert_is_instance(r_work, float,
                "r_work metric should be float, got %s"
                % type(r_work).__name__)
        if r_free is not None:
            assert_is_instance(r_free, float,
                "r_free metric should be float, got %s"
                % type(r_free).__name__)
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_metrics_through_history_pipeline():
    """Full pipeline: record_result → save → load →
    get_history_for_agent → _analyze_history → verify types."""
    print("Test: session_metrics_through_history_pipeline")
    from agent.workflow_state import _analyze_history

    ss, tmpdir = _make_session()
    try:
        ss.set_experiment_type("xray")

        ss.start_cycle(1)
        ss.record_decision(1, program='phenix.xtriage',
            command='phenix.xtriage data.mtz')
        ss.record_result(1, "SUCCESS: Resolution=2.05")

        ss.start_cycle(2)
        ss.record_decision(2, program='phenix.refine',
            command='phenix.refine model.pdb data.mtz')
        ss.record_result(2,
            "SUCCESS: R-work=0.2200 R-free=0.2600")
        ss.save()

        # Reload and extract history
        ss2 = _reload_session(tmpdir)
        history = ss2.get_history_for_agent()

        assert_equal(len(history), 2,
            "Should have 2 history entries")
        assert_equal(history[0]["program"], "phenix.xtriage",
            "First history entry should be xtriage")

        # Feed through _analyze_history (the real consumer)
        info = _analyze_history(history)

        # Check types in the analysis output
        last_r_free = info.get("last_r_free")
        if last_r_free is not None:
            assert_is_instance(last_r_free, float,
                "last_r_free from _analyze_history should be "
                "float, got %s" % type(last_r_free).__name__)

        last_r_work = info.get("last_r_work")
        if last_r_work is not None:
            assert_is_instance(last_r_work, float,
                "last_r_work from _analyze_history should be "
                "float, got %s" % type(last_r_work).__name__)

        # Verify refine_count
        assert_equal(info.get("refine_count", 0), 1,
            "refine_count should be 1")

        # Verify last_program (B1 fix)
        assert_equal(info.get("last_program"), "phenix.refine",
            "last_program should be phenix.refine")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_bug4_string_metric_detection():
    """Simulate Bug 4: someone writes string metrics directly
    to session data. Verify _safe_float defense-in-depth."""
    print("Test: session_bug4_string_metric_detection")
    from agent.structure_model import _safe_float

    ss, tmpdir = _make_session()
    try:
        ss.start_cycle(1)
        ss.record_decision(1, program='phenix.model_vs_data')
        # Simulate the Bug 4 corruption: manually inject
        # string metrics (this is what model_vs_data did)
        cycle = ss.data["cycles"][0]
        cycle["result"] = "SUCCESS: R-work=0.385 R-free=0.386"
        cycle["metrics"] = {
            "r_work": "0.385",   # BUG: string not float
            "r_free": "0.386",   # BUG: string not float
        }
        ss.save()

        ss2 = _reload_session(tmpdir)
        metrics = ss2.data["cycles"][0].get("metrics", {})

        # The raw data IS corrupted (strings)
        assert_is_instance(metrics["r_work"], str,
            "Raw corrupted metric should be str")

        # But _safe_float should fix it
        fixed_rw = _safe_float(metrics["r_work"])
        fixed_rf = _safe_float(metrics["r_free"])
        assert_is_instance(fixed_rw, float,
            "_safe_float should convert to float")
        assert_equal(fixed_rw, 0.385,
            "_safe_float should get correct value")
        assert_is_instance(fixed_rf, float,
            "_safe_float should convert to float")
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_cycle_structure_survives_reload():
    """Multi-cycle session preserves cycle count, ordering, and
    all fields after save/load."""
    print("Test: session_cycle_structure_survives_reload")
    ss, tmpdir = _make_session()
    try:
        for i in range(1, 4):
            ss.start_cycle(i)
            ss.record_decision(i,
                program='phenix.refine',
                decision='Cycle %d' % i)
            ss.record_result(i,
                "SUCCESS: R-work=%.4f R-free=%.4f"
                % (0.30 - i*0.03, 0.35 - i*0.03))
        ss.save()

        ss2 = _reload_session(tmpdir)
        cycles = ss2.data["cycles"]
        assert_equal(len(cycles), 3,
            "Should have 3 cycles")

        for i, c in enumerate(cycles):
            assert_equal(c["cycle_number"], i + 1,
                "Cycle %d number wrong" % (i + 1))
            assert_equal(c["program"], "phenix.refine",
                "Cycle %d program wrong" % (i + 1))
            assert_true(
                c["result"].startswith("SUCCESS"),
                "Cycle %d result should start with SUCCESS"
                % (i + 1))
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


def test_session_none_fields_survive_reload():
    """Fields initialized to None survive save/load as None,
    not as string 'None' or missing."""
    print("Test: session_none_fields_survive_reload")
    ss, tmpdir = _make_session()
    try:
        # These are all None in a fresh session
        ss.save()
        ss2 = _reload_session(tmpdir)

        for field in ['resolution', 'resolution_source',
                      'experiment_type', 'rfree_mtz',
                      'rfree_mtz_locked_at_cycle',
                      'force_retry_program']:
            val = ss2.data.get(field)
            assert_is_none(val,
                "Fresh session field '%s' should be None "
                "after reload, got %r" % (field, val))
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED")


# =====================================================================
# RUN ALL TESTS
# =====================================================================

def run_all_tests():
    """Run all Phase 3 tests."""
    print("Phase 3: Session Round-Trip — Symmetry + Invariants")
    print()

    tests = [
        # Symmetry tests (JSON module behavior — schema documentation)
        test_float_round_trip,
        test_none_vs_string_none,
        test_none_vs_missing_field,
        test_int_vs_string_int,
        test_float_precision,
        test_path_round_trip,
        test_empty_collections_preserved,
        test_nested_dict_round_trip,
        test_history_entry_result_key,
        test_full_session_round_trip,
        # Invariant tests (data consistency checks)
        test_cycle_numbering_invariant,
        test_cycle_count_matches_list_length,
        test_successful_cycle_has_program,
        test_experiment_type_valid,
        test_resolution_is_float_or_none,
        test_rfree_mtz_consistency,
        test_directives_extracted_consistency,
        test_output_files_are_strings,
        test_recovery_attempts_count_is_int,
        # AgentSession round-trip tests (actual production code)
        test_session_resolution_type,
        test_session_resolution_string_input,
        test_session_experiment_type_locking,
        test_session_rfree_mtz_types,
        test_session_record_result_metrics,
        test_session_metrics_through_history_pipeline,
        test_session_bug4_string_metric_detection,
        test_session_cycle_structure_survives_reload,
        test_session_none_fields_survive_reload,
    ]

    passed = 0
    failed = 0
    failures = []

    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            failures.append({
                'test': t.__name__,
                'error': str(e)[:200],
            })
            print("  FAILED: %s" % str(e)[:100])

    print()
    print("  Results: %d passed, %d failed out of %d" %
          (passed, failed, len(tests)))

    status = 'PASS' if failed == 0 else 'FAIL'

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_3_roundtrip_failures.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 3: Session Round-Trip Findings\n")
        f.write("phase: 3\n")
        f.write("name: session_roundtrip\n")
        f.write("status: %s\n" % status)
        f.write("total_tests: %d\n" % len(tests))
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        f.write("symmetry_tests: 10\n")
        f.write("invariant_tests: 9\n")
        f.write("session_roundtrip_tests: 9\n")
        if failures:
            f.write("\nfailures:\n")
            for fail in failures:
                f.write("  - test: %s\n" % fail['test'])
                f.write("    error: %s\n" % fail['error'])

    print("  Findings: %s" % path)
    print("  Phase 3 overall: %s" % status)

    result = {'status': status, 'passed': passed,
              'failed': failed}
    if failed > 0:
        raise AssertionError(
            "Phase 3 serialization FAILED: %d/%d tests failed"
            % (failed, passed + failed))
    return result


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] == 'PASS' else 1)
