"""Tests for v116.10 Phase 2 — protocol version invariants in contract.py.

The bug: contract.py had CURRENT_PROTOCOL_VERSION=3 but registered fields
existed at version 4 (plan_has_pending_stages) and version 5 (asu_copies).
This drift meant `get_deprecation_warnings` would tell v3 clients they
were "current" when in fact they were missing newer fields.

The fix:
1. Bumped CURRENT_PROTOCOL_VERSION from 3 to 5.
2. Added contract.validate_contract() which returns (ok, errors) for
   any violation of the protocol invariants.
3. These tests assert that validate_contract() reports ok at every
   invariant we care about, plus probe the function with synthetic
   drift scenarios.

Run from `tests/` directory or in PHENIX environment.  The contract
module has no libtbx dependencies so it imports cleanly standalone.
"""

from __future__ import absolute_import, division, print_function

import os
import sys


# --- Path setup -----------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_AGENT_DIR = os.path.normpath(os.path.join(_HERE, "..", "agent"))
if _AGENT_DIR not in sys.path:
    sys.path.insert(0, _AGENT_DIR)

try:
    from libtbx.langchain.agent import contract
except ImportError:
    import contract


# =====================================================================
# SECTION A: Live contract.py invariants
# =====================================================================

def test_contract_passes_validation():
    """The shipped contract.py must satisfy its own invariants."""
    print("Test: contract_passes_validation")
    ok, errors = contract.validate_contract()
    assert ok, (
        "Contract validation failed; errors:\n  %s" % "\n  ".join(errors))
    print("  PASS")


def test_current_protocol_version_at_least_max_field():
    """CURRENT_PROTOCOL_VERSION must cover every registered field."""
    print("Test: current_protocol_version_at_least_max_field")
    max_field = max(ver for _, _, ver, _ in contract.SESSION_INFO_FIELDS)
    assert contract.CURRENT_PROTOCOL_VERSION >= max_field, (
        "CURRENT_PROTOCOL_VERSION (%d) must be >= max field version (%d). "
        "Bump the constant in contract.py."
        % (contract.CURRENT_PROTOCOL_VERSION, max_field))
    print("  PASS")


def test_min_supported_at_most_current():
    """MIN_SUPPORTED_PROTOCOL_VERSION must be <= CURRENT."""
    print("Test: min_supported_at_most_current")
    assert (contract.MIN_SUPPORTED_PROTOCOL_VERSION
            <= contract.CURRENT_PROTOCOL_VERSION), (
        "MIN_SUPPORTED (%d) > CURRENT (%d) — no client can satisfy this."
        % (contract.MIN_SUPPORTED_PROTOCOL_VERSION,
           contract.CURRENT_PROTOCOL_VERSION))
    print("  PASS")


def test_min_supported_at_least_one():
    """MIN_SUPPORTED must be >= 1 (versions start at 1)."""
    print("Test: min_supported_at_least_one")
    assert contract.MIN_SUPPORTED_PROTOCOL_VERSION >= 1, (
        "MIN_SUPPORTED (%d) < 1 is undefined."
        % contract.MIN_SUPPORTED_PROTOCOL_VERSION)
    print("  PASS")


def test_current_protocol_version_is_8():
    """v120.2 bumped CURRENT to 8.

    Change-detector: this documents the current protocol version as a
    conscious decision.  v116.10 Phase 2 set it to 5; v120 bumped it to 6
    for the plan_current_unrun_lead_program field (Option 2a reactive-
    deviation hold), then to 7 for the input_mtz_has_rfree field
    (client-extracted R-free presence), then to 8 for the mtz_rfree_map
    field (per-file R-free map, v120.2 parity fix).
    If a future phase bumps to 9 or higher, update this test (or remove it).
    """
    print("Test: current_protocol_version_is_8")
    assert contract.CURRENT_PROTOCOL_VERSION == 8, (
        "v120.2 set CURRENT_PROTOCOL_VERSION to 8; "
        "current value is %d. If this is intentional (e.g. v9 added), "
        "update this test." % contract.CURRENT_PROTOCOL_VERSION)
    print("  PASS")


# =====================================================================
# SECTION B: validate_contract() probing
# =====================================================================

def test_validate_returns_tuple_of_correct_shape():
    """validate_contract returns (bool, list[str])."""
    print("Test: validate_returns_tuple_of_correct_shape")
    result = contract.validate_contract()
    assert isinstance(result, tuple) and len(result) == 2, (
        "validate_contract should return a 2-tuple, got %r" % (result,))
    ok, errors = result
    assert isinstance(ok, bool), "First element should be bool, got %r" % type(ok)
    assert isinstance(errors, list), "Second element should be list, got %r" % type(errors)
    assert all(isinstance(e, str) for e in errors), (
        "errors must be a list of strings; got %r" % errors)
    print("  PASS")


def test_validate_detects_current_below_max_field():
    """Synthetic test: simulate CURRENT < max field via monkey-patching.

    Verifies the function detects drift, not just that the current
    contract happens to be valid.
    """
    print("Test: validate_detects_current_below_max_field")
    original_current = contract.CURRENT_PROTOCOL_VERSION
    try:
        # Force drift: set CURRENT below the max field version
        contract.CURRENT_PROTOCOL_VERSION = 1
        ok, errors = contract.validate_contract()
        assert not ok, (
            "validate_contract should report failure when CURRENT < max field, "
            "but returned ok=%r errors=%r" % (ok, errors))
        # The error message should mention the issue
        joined = " ".join(errors).lower()
        assert "current_protocol_version" in joined or "current" in joined, (
            "Error message should mention CURRENT_PROTOCOL_VERSION; got: %r"
            % errors)
    finally:
        contract.CURRENT_PROTOCOL_VERSION = original_current
    print("  PASS")


def test_validate_detects_min_above_current():
    """Synthetic test: simulate MIN > CURRENT."""
    print("Test: validate_detects_min_above_current")
    original_min = contract.MIN_SUPPORTED_PROTOCOL_VERSION
    original_current = contract.CURRENT_PROTOCOL_VERSION
    try:
        # Force invariant violation: MIN > CURRENT
        contract.MIN_SUPPORTED_PROTOCOL_VERSION = original_current + 1
        ok, errors = contract.validate_contract()
        assert not ok, (
            "validate_contract should report failure when MIN > CURRENT, "
            "but returned ok=%r errors=%r" % (ok, errors))
        joined = " ".join(errors).lower()
        assert "min" in joined and "current" in joined, (
            "Error message should mention MIN and CURRENT; got: %r" % errors)
    finally:
        contract.MIN_SUPPORTED_PROTOCOL_VERSION = original_min
    print("  PASS")


def test_validate_detects_min_below_one():
    """Synthetic test: MIN < 1 is invalid."""
    print("Test: validate_detects_min_below_one")
    original_min = contract.MIN_SUPPORTED_PROTOCOL_VERSION
    try:
        contract.MIN_SUPPORTED_PROTOCOL_VERSION = 0
        ok, errors = contract.validate_contract()
        assert not ok, (
            "validate_contract should reject MIN=0, but returned "
            "ok=%r errors=%r" % (ok, errors))
    finally:
        contract.MIN_SUPPORTED_PROTOCOL_VERSION = original_min
    print("  PASS")


def test_validate_succeeds_with_current_above_max_field():
    """Synthetic test: CURRENT > max field is fine (anticipating v6+).

    Even though no v6 fields exist yet, bumping CURRENT to 6 in
    advance is permitted (e.g. reserving the version).
    """
    print("Test: validate_succeeds_with_current_above_max_field")
    original_current = contract.CURRENT_PROTOCOL_VERSION
    try:
        contract.CURRENT_PROTOCOL_VERSION = 99
        ok, errors = contract.validate_contract()
        assert ok, (
            "validate_contract should accept CURRENT > max field; "
            "got ok=%r errors=%r" % (ok, errors))
    finally:
        contract.CURRENT_PROTOCOL_VERSION = original_current
    print("  PASS")


# =====================================================================
# SECTION C: Field registry invariants
# =====================================================================

def test_field_versions_are_positive():
    """All registered fields have version >= 1."""
    print("Test: field_versions_are_positive")
    bad = [(name, ver) for name, _, ver, _
           in contract.SESSION_INFO_FIELDS if ver < 1]
    assert not bad, (
        "These fields have invalid versions (must be >= 1): %s" % bad)
    print("  PASS")


def test_field_names_are_unique():
    """No two SESSION_INFO_FIELDS share a name."""
    print("Test: field_names_are_unique")
    names = [name for name, _, _, _ in contract.SESSION_INFO_FIELDS]
    duplicates = [n for n in set(names) if names.count(n) > 1]
    assert not duplicates, "Duplicate field names: %s" % duplicates
    print("  PASS")


def test_field_entries_are_4_tuples():
    """Each SESSION_INFO_FIELDS entry is (name, default, version, desc)."""
    print("Test: field_entries_are_4_tuples")
    for entry in contract.SESSION_INFO_FIELDS:
        assert isinstance(entry, tuple) and len(entry) == 4, (
            "Field entry must be a 4-tuple; got %r" % (entry,))
        name, _default, version, desc = entry
        assert isinstance(name, str), "field name must be str, got %r" % (entry,)
        assert isinstance(version, int), "field version must be int, got %r" % (entry,)
        assert isinstance(desc, str), "field description must be str, got %r" % (entry,)
    print("  PASS")


# =====================================================================
# SECTION D: normalize_session_info safety
# =====================================================================

def test_normalize_does_not_strip_unknown_fields():
    """normalize_session_info adds defaults but doesn't remove unknown fields.

    This guarantees forward compatibility: if a newer client sends a v6 field
    we don't know about yet, the server preserves it on the way through (in
    case downstream code needs it).
    """
    print("Test: normalize_does_not_strip_unknown_fields")
    s = {"some_future_field": "futureval",
         "another_future_field": [1, 2, 3]}
    result = contract.normalize_session_info(s)
    assert "some_future_field" in result and result["some_future_field"] == "futureval"
    assert "another_future_field" in result and result["another_future_field"] == [1, 2, 3]
    print("  PASS")


def test_normalize_fills_defaults_for_missing_fields():
    """normalize_session_info fills in all registered defaults."""
    print("Test: normalize_fills_defaults_for_missing_fields")
    s = {}
    result = contract.normalize_session_info(s)
    for name, default, _ver, _desc in contract.SESSION_INFO_FIELDS:
        assert name in result, "Missing default for field %s" % name
        # For mutable types, the default must be a copy not the same object
        if isinstance(default, (dict, list)):
            assert result[name] is not default, (
                "Mutable default for %s was not copied "
                "(this could cause cross-request contamination)" % name)
            assert result[name] == default, (
                "Default value mismatch for %s: got %r, want %r"
                % (name, result[name], default))
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        try:
            from tests.tst_utils import run_tests_with_fail_fast
        except ImportError:
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
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


if __name__ == "__main__":
    _standalone_runner()
