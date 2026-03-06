"""
Backward compatibility tests using frozen client fixtures.

These fixtures represent REAL payloads from specific client versions.
They must NEVER be updated — they are the "old clients" we test against.

Each test loads a fixture, runs it through the server's normalization
and version checking, and verifies no crash and correct behavior.

When running in a full PHENIX environment, the tests also run the
complete perceive→plan→build→validate→output pipeline.

Run with:
    PYTHONPATH=. python tests/tst_old_client_compat.py
"""


import glob
import json
import os
import sys

# Add parent directory to path for imports
sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import assert_equal
from tests.tst_utils import assert_true
from tests.tst_utils import assert_false
from tests.tst_utils import assert_in
from tests.tst_utils import run_tests_with_fail_fast

# =========================================================================
# PHENIX/cctbx Linter "Silencer"
# =========================================================================
(assert_equal, assert_true, assert_false, assert_in, run_tests_with_fail_fast)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_FIXTURES_DIR = os.path.join(_PROJECT_ROOT, "tests", "fixtures")

try:
    from agent.contract import (
        normalize_session_info, check_client_version,
        get_deprecation_warnings, SESSION_INFO_FIELDS,
        CURRENT_PROTOCOL_VERSION,
    )
    CONTRACT_AVAILABLE = True
except ImportError:
    CONTRACT_AVAILABLE = False


def _load_fixtures():
    """Load all fixture JSON files."""
    pattern = os.path.join(_FIXTURES_DIR, "client_v*", "*.json")
    fixtures = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        data["_fixture_path"] = path
        fixtures.append(data)
    return fixtures


# ---------------------------------------------------------------------------
# Test: All fixtures load and have required structure
# ---------------------------------------------------------------------------

def test_fixtures_valid_structure():
    """Every fixture must have the expected top-level keys."""
    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return
    assert_true(len(fixtures) >= 2, "Expected at least 2 fixtures, found %d" % len(fixtures))

    required_keys = {"client_version", "scenario", "decide_next_step_args", "expected"}
    required_args = {"log_content", "history", "files", "session_resolution",
                     "session_info", "abort_on_red_flags", "abort_on_warnings"}

    for fix in fixtures:
        path = os.path.basename(fix["_fixture_path"])
        missing = required_keys - set(fix.keys())
        assert_true(not missing, "%s missing top-level keys: %s" % (path, missing))

        args = fix["decide_next_step_args"]
        missing_args = required_args - set(args.keys())
        assert_true(not missing_args,
                    "%s missing args: %s" % (path, missing_args))

    print("  PASSED: %d fixtures have valid structure" % len(fixtures))


# ---------------------------------------------------------------------------
# Test: normalize_session_info fills missing fields for old clients
# ---------------------------------------------------------------------------

def test_normalize_fills_missing_fields():
    """
    For each fixture, normalize_session_info must add defaults for any
    missing SESSION_INFO_FIELDS without crashing or modifying existing fields.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    all_field_names = {name for name, _, _, _ in SESSION_INFO_FIELDS}
    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return

    for fix in fixtures:
        path = os.path.basename(fix["_fixture_path"])
        si = dict(fix["decide_next_step_args"]["session_info"])  # copy
        original_keys = set(si.keys())

        si = normalize_session_info(si)

        # All registered fields must now exist
        for field_name in all_field_names:
            assert_true(field_name in si,
                        "%s: field '%s' missing after normalize" % (path, field_name))

        # Original fields must be unchanged
        original_si = fix["decide_next_step_args"]["session_info"]
        for key in original_keys:
            if key in all_field_names:
                assert_equal(si[key], original_si[key],
                             "%s: field '%s' was modified by normalize" % (path, key))

    print("  PASSED: normalize fills missing fields for all %d fixtures" % len(fixtures))


# ---------------------------------------------------------------------------
# Test: Unknown/junk fields are preserved
# ---------------------------------------------------------------------------

def test_junk_fields_preserved():
    """
    Fixtures with unknown fields (e.g. junk_field_from_future) must not
    have them stripped or cause errors.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return
    junk_tested = False

    for fix in fixtures:
        si = dict(fix["decide_next_step_args"]["session_info"])
        junk_keys = {k for k in si.keys()
                     if k not in {n for n, _, _, _ in SESSION_INFO_FIELDS}}

        if junk_keys:
            si = normalize_session_info(si)
            for jk in junk_keys:
                assert_true(jk in si,
                            "Junk field '%s' was stripped by normalize" % jk)
            junk_tested = True

    assert_true(junk_tested,
                "No fixture has junk fields — add one for proper testing")
    print("  PASSED: Junk/unknown fields preserved through normalization")


# ---------------------------------------------------------------------------
# Test: Version check accepts all fixtures
# ---------------------------------------------------------------------------

def test_version_check_accepts_all():
    """
    All fixtures should be accepted by check_client_version
    (since MIN_SUPPORTED is 1 and all fixtures are >= v1).
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return

    for fix in fixtures:
        path = os.path.basename(fix["_fixture_path"])
        si = dict(fix["decide_next_step_args"]["session_info"])
        si = normalize_session_info(si)

        result = check_client_version(si)
        assert_true(result is None,
                    "%s: should be accepted but got rejection: %s" % (path, result))

    print("  PASSED: All %d fixtures accepted by version check" % len(fixtures))


# ---------------------------------------------------------------------------
# Test: Deprecation warnings for old clients
# ---------------------------------------------------------------------------

def test_deprecation_warnings_for_old_clients():
    """
    Fixtures with protocol version < CURRENT should get deprecation warnings.
    Fixtures at CURRENT should get none.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return
    warned = 0
    clean = 0

    for fix in fixtures:
        path = os.path.basename(fix["_fixture_path"])
        si = dict(fix["decide_next_step_args"]["session_info"])
        si = normalize_session_info(si)

        warnings = get_deprecation_warnings(si)
        client_v = si.get("client_protocol_version", 1)

        if client_v < CURRENT_PROTOCOL_VERSION:
            assert_true(len(warnings) >= 1,
                        "%s (v%d): expected deprecation warning" % (path, client_v))
            warned += 1
        else:
            assert_equal(len(warnings), 0,
                         "%s (v%d): should have no warnings" % (path, client_v))
            clean += 1

    assert_true(warned > 0, "Need at least one old-client fixture to test warnings")
    print("  PASSED: %d fixtures got warnings, %d clean" % (warned, clean))


# ---------------------------------------------------------------------------
# Test: Full pipeline (only in PHENIX environment)
# ---------------------------------------------------------------------------

def test_full_pipeline_if_available():
    """
    If the full PHENIX graph is importable, run each fixture through
    perceive→plan→build→validate→output and check expectations.

    Skips gracefully outside PHENIX.
    """
    try:
        from agent.graph_nodes import perceive, plan, build, validate, output_node
    except ImportError:
        print("  SKIP (full pipeline not importable — not in PHENIX environment)")
        return

    fixtures = _load_fixtures()
    if not fixtures:
        print("  SKIP (no fixtures in %s)" % _FIXTURES_DIR)
        return
    passed = 0

    for fix in fixtures:
        path = os.path.basename(fix["_fixture_path"])
        args = fix["decide_next_step_args"]
        expected = fix["expected"]

        # Build initial state matching what the REST handler creates
        state = {
            "log_text": args["log_content"],
            "history": args["history"],
            "available_files": args["files"],
            "guidelines": args.get("guidelines", ""),
            "session_resolution": args["session_resolution"],
            "session_info": dict(args["session_info"]),
            "abort_on_red_flags": args["abort_on_red_flags"],
            "abort_on_warnings": args["abort_on_warnings"],
            "cycle_number": len(args["history"]) + 1,
            "directives": args["session_info"].get("directives", {}),
            "bad_inject_params": args["session_info"].get("bad_inject_params", {}),
        }

        try:
            state = perceive(state)
            state = plan(state)
            state = build(state)
            state = validate(state)
            state = output_node(state)
        except Exception as e:
            assert_true(False, "%s: Pipeline crashed: %s" % (path, e))

        if expected.get("should_not_crash"):
            pass  # We got here, so it didn't crash

        if expected.get("stop_not_set"):
            assert_false(state.get("stop"),
                         "%s: Unexpected stop=%s reason=%s" %
                         (path, state.get("stop"), state.get("stop_reason")))

        command = state.get("command", "")
        if expected.get("command_not_empty"):
            assert_true(bool(command), "%s: Empty command" % path)

        if expected.get("program_in"):
            program = state.get("program") or (state.get("intent", {}).get("program"))
            assert_in(program, expected["program_in"],
                      "%s: program=%s not in %s" % (path, program, expected["program_in"]))

        passed += 1

    print("  PASSED: Full pipeline ran successfully for %d fixtures" % passed)


# =========================================================================
# Run
# =========================================================================

def run_all_tests():
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
