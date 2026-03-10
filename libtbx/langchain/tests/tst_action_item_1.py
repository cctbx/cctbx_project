"""
Tests for Action Item 1 (P1B):
  validate() checks think_file_overrides existence before BUILD consumes them.

Following the established test pattern in this codebase, the check logic is
replicated inline rather than importing graph_nodes directly (graph_nodes has
unconditional libtbx imports that are not available in the test environment).

Seven cases:
  1. Override path found by exact match -> no error, overrides preserved
  2. Override path found by basename match -> no error, overrides preserved
  3. Override path not in available_files -> error set, overrides cleared,
     attempt incremented
  4. Empty / None / [] overrides -> no-op, no error
  5. Override with None/empty path value -> silently skipped
  6. List of paths all missing -> error
  7. List of paths with any missing -> error (all paths must be present)
"""
import sys, os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


# ---------------------------------------------------------------------------
# Inline replication of the P1B check block from validate()
# (mirrors the code in agent/graph_nodes.py :: validate())
# ---------------------------------------------------------------------------

def _increment_attempt(state):
    """Replicated from graph_nodes._increment_attempt."""
    attempt = {
        "intent":   state.get("intent", {}),
        "command":  state.get("command", ""),
        "error":    state.get("validation_error", ""),
    }
    previous = list(state.get("previous_attempts", []))
    return {
        **state,
        "attempt_number":    state.get("attempt_number", 0) + 1,
        "previous_attempts": previous + [attempt],
    }


def _check_file_overrides(state):
    """
    Replicated P1B check from validate() in graph_nodes.py.

    Returns (state, failed) where failed=True means an invalid override
    was found and state already has validation_error set + think_file_overrides
    cleared + attempt_number incremented.
    """
    available_files   = state.get("available_files", [])
    available_set     = set(available_files)
    available_basenames = {os.path.basename(f): f for f in available_set}

    _think_overrides = state.get("think_file_overrides") or {}
    if not _think_overrides:
        return state, False

    _bad_overrides = {}
    for _cat, _path in _think_overrides.items():
        if not _path:
            continue
        _paths = _path if isinstance(_path, list) else [_path]
        for _p in _paths:
            _bn = os.path.basename(_p)
            if _p not in available_set and _bn not in available_basenames:
                _bad_overrides[_cat] = _p
                break

    if not _bad_overrides:
        return state, False

    _detail = "; ".join(
        "category '%s' -> '%s'" % (c, p)
        for c, p in sorted(_bad_overrides.items())
    )
    error = (
        "THINK suggested file_override(s) not found in available_files: "
        "%s. Verify the path is correct and the file has been produced "
        "by a previous step." % _detail
    )
    state = {**state, "think_file_overrides": {}, "validation_error": error}
    state = _increment_attempt(state)
    return state, True


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _state(overrides, available_files):
    return {
        "available_files":    available_files,
        "think_file_overrides": overrides,
        "previous_attempts":  [],
        "attempt_number":     0,
        "intent":             {},
        "command":            "phenix.refine model.pdb data.mtz",
        "validation_error":   None,
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_exact_path_ok():
    """Override path in available_files by exact match -> no error."""
    available = ["/work/model.pdb", "/work/data.mtz", "/work/rebuilt.pdb"]
    state = _state({"model": "/work/rebuilt.pdb"}, available)
    result, failed = _check_file_overrides(state)

    assert not failed, "Should not fail for exact path match"
    assert result.get("validation_error") is None
    assert result.get("think_file_overrides") == {"model": "/work/rebuilt.pdb"}, \
        "Overrides should be preserved on success"
    assert result.get("attempt_number", 0) == 0
    print("PASS test_exact_path_ok")


def test_basename_match_ok():
    """Override path found by basename even when directory differs -> no error."""
    available = [
        "/work/model.pdb",
        "/work/data.mtz",
        "/some/other/path/rebuilt.pdb",   # different directory
    ]
    state = _state({"model": "/think/suggested/rebuilt.pdb"}, available)
    result, failed = _check_file_overrides(state)

    assert not failed, "Should not fail for basename match"
    assert result.get("validation_error") is None
    assert result.get("think_file_overrides") == {"model": "/think/suggested/rebuilt.pdb"}
    print("PASS test_basename_match_ok")


def test_missing_path_triggers_error():
    """Override path not in available_files -> error, overrides cleared, attempt++."""
    available = ["/work/model.pdb", "/work/data.mtz"]
    state = _state({"model": "/work/non_existent.pdb"}, available)
    result, failed = _check_file_overrides(state)

    assert failed, "Should fail for missing path"

    err = result.get("validation_error", "")
    assert err, "validation_error must be set"
    assert "model" in err,            "Error must name the category: %r" % err
    assert "non_existent.pdb" in err, "Error must name the path: %r" % err
    assert "available_files" in err,  "Error must mention available_files: %r" % err

    assert result.get("think_file_overrides") == {}, \
        "think_file_overrides must be cleared, got: %r" % result.get("think_file_overrides")
    assert result.get("attempt_number", 0) == 1, \
        "attempt_number must increment to 1"

    prev = result.get("previous_attempts", [])
    assert len(prev) == 1, "Expected 1 previous_attempt entry, got %d" % len(prev)
    assert prev[0].get("error") == err, \
        "previous_attempts error must match validation_error"
    print("PASS test_missing_path_triggers_error")


def test_empty_overrides_noop():
    """Empty / None / list overrides -> no-op, no error."""
    available = ["/work/model.pdb", "/work/data.mtz"]
    for overrides in ({}, None, []):
        state = _state(overrides, available)
        result, failed = _check_file_overrides(state)
        assert not failed, "Empty overrides %r should not fail" % (overrides,)
        assert result.get("validation_error") is None
    print("PASS test_empty_overrides_noop")


def test_null_path_value_skipped():
    """An override with a None/empty path value is silently skipped."""
    available = ["/work/model.pdb", "/work/data.mtz"]
    state = _state({"model": None, "sequence": ""}, available)
    result, failed = _check_file_overrides(state)
    assert not failed, "None/empty paths should be skipped"
    assert result.get("validation_error") is None
    print("PASS test_null_path_value_skipped")


def test_list_path_value_all_missing():
    """Override with a list of paths: all missing -> error."""
    available = ["/work/model.pdb", "/work/data.mtz"]
    state = _state({"model": ["/work/missing_a.pdb", "/work/missing_b.pdb"]}, available)
    result, failed = _check_file_overrides(state)
    assert failed, "Should fail when list paths are all missing"
    assert result.get("think_file_overrides") == {}
    print("PASS test_list_path_value_all_missing")


def test_list_path_value_any_missing_triggers_error():
    """Override with a list where first path exists but second does not -> error.

    The check iterates all paths and flags the category if ANY path is missing.
    A list override is only valid when all paths in it are present.
    """
    available = ["/work/model.pdb", "/work/data.mtz", "/work/rebuilt.pdb"]
    state = _state({"model": ["/work/rebuilt.pdb", "/work/missing.pdb"]}, available)
    result, failed = _check_file_overrides(state)
    # Second path is missing -> category is flagged
    assert failed, "Should fail when any path in list is missing"
    assert result.get("think_file_overrides") == {}
    print("PASS test_list_path_value_any_missing_triggers_error")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def run_all_tests():
    for t in [
        test_exact_path_ok,
        test_basename_match_ok,
        test_missing_path_triggers_error,
        test_empty_overrides_noop,
        test_null_path_value_skipped,
        test_list_path_value_all_missing,
        test_list_path_value_any_missing_triggers_error,
    ]:
        t()
    print("All 7 tests passed.")

if __name__ == "__main__":
    tests = [
        test_exact_path_ok,
        test_basename_match_ok,
        test_missing_path_triggers_error,
        test_empty_overrides_noop,
        test_null_path_value_skipped,
        test_list_path_value_all_missing,
        test_list_path_value_any_missing_triggers_error,
    ]
    failures = []
    for t in tests:
        try:
            t()
        except Exception as e:
            import traceback
            print("FAIL %s: %s" % (t.__name__, e))
            traceback.print_exc()
            failures.append(t.__name__)
    print()
    if failures:
        print("FAILED: %s" % ", ".join(failures))
        sys.exit(1)
    else:
        print("All %d tests passed." % len(tests))
