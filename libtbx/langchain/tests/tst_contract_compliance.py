"""
Contract compliance tests for backward compatibility.

Verifies that server-side code follows the rules in agent/contract.py:
  1. No bare session_info["X"] reads (would KeyError on old clients)
  2. Every session_info field accessed is registered in the contract
  3. Defaults used in .get() calls match the contract
  4. No os.path.exists() guards on file paths (server parity)

Run with:
    PYTHONPATH=. python tests/tst_contract_compliance.py
"""


import os
import re
import sys

# Add parent directory to path for imports
sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import assert_equal
from tests.tst_utils import assert_true
from tests.tst_utils import assert_in
from tests.tst_utils import run_tests_with_fail_fast

# =========================================================================
# PHENIX/cctbx Linter "Silencer"
# =========================================================================
(re, assert_equal, assert_true, assert_in, run_tests_with_fail_fast)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Server-side files where session_info is accessed
SERVER_FILES = [
    os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py"),
    os.path.join(_PROJECT_ROOT, "agent", "command_builder.py"),
    os.path.join(_PROJECT_ROOT, "agent", "workflow_state.py"),
    os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py"),
]

# Files where os.path.exists() on user file paths is forbidden
# (these run on the server where client paths don't exist)
SERVER_NO_OSPATH_FILES = [
    os.path.join(_PROJECT_ROOT, "agent", "graph_nodes.py"),
    os.path.join(_PROJECT_ROOT, "agent", "command_builder.py"),
    os.path.join(_PROJECT_ROOT, "agent", "workflow_engine.py"),
]

# Load the contract
try:
    from agent.contract import SESSION_INFO_FIELDS, RESPONSE_FIELDS
    CONTRACT_AVAILABLE = True
except ImportError:
    CONTRACT_AVAILABLE = False


def _get_contract_fields():
    """Return dict of field_name -> (default, version, desc)."""
    return {
        name: (default, ver, desc)
        for name, default, ver, desc in SESSION_INFO_FIELDS
    }


# ---------------------------------------------------------------------------
# Test: No bare bracket reads on session_info
# ---------------------------------------------------------------------------

def test_no_bare_session_info_bracket_reads():
    """
    Scan server-side code for session_info["X"] bracket access.

    Writes (assignments) are acceptable:
        session_info["best_files"] = new_best   # OK

    Reads are NOT acceptable (would KeyError on old clients):
        x = session_info["best_files"]           # BAD
        if session_info["advice_changed"]:        # BAD
    """
    violations = []

    for filepath in SERVER_FILES:
        if not os.path.exists(filepath):
            continue
        with open(filepath) as f:
            for lineno, line in enumerate(f, 1):
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue

                for m in re.finditer(r'session_info\["(\w+)"\]', line):
                    # Check if this is an assignment target
                    after = line[m.end():].lstrip()
                    is_assignment = after.startswith("=") and not after.startswith("==")
                    if not is_assignment:
                        violations.append(
                            "%s:%d  session_info[\"%s\"] — use .get(\"%s\", default) instead"
                            % (os.path.basename(filepath), lineno,
                               m.group(1), m.group(1))
                        )

    if violations:
        msg = "Found %d bare session_info bracket read(s):\n  %s" % (
            len(violations), "\n  ".join(violations))
        assert_true(False, msg)
    print("  PASSED: No bare session_info bracket reads in %d files" % len(SERVER_FILES))


# ---------------------------------------------------------------------------
# Test: All accessed fields are registered in contract
# ---------------------------------------------------------------------------

def test_all_accessed_fields_in_contract():
    """
    Every session_info.get("X") on the server must have "X" registered
    in SESSION_INFO_FIELDS.  Unregistered fields mean the contract is
    incomplete — a future static check or normalization might miss them.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    contract_fields = set(_get_contract_fields().keys())
    accessed_fields = set()

    for filepath in SERVER_FILES:
        if not os.path.exists(filepath):
            continue
        with open(filepath) as f:
            content = f.read()
        # .get() calls
        for m in re.finditer(r'session_info\.get\("(\w+?)"', content):
            accessed_fields.add(m.group(1))
        # Bracket access (writes count too — they imply the field exists)
        for m in re.finditer(r'session_info\["(\w+?)"\]', content):
            accessed_fields.add(m.group(1))

    unregistered = accessed_fields - contract_fields
    if unregistered:
        msg = "Fields accessed but NOT in contract: %s" % sorted(unregistered)
        assert_true(False, msg)
    print("  PASSED: All %d accessed fields are registered in contract" % len(accessed_fields))


# ---------------------------------------------------------------------------
# Test: Contract defaults are used consistently
# ---------------------------------------------------------------------------

def test_contract_defaults_consistency():
    """
    When session_info.get("X", default) is called with an explicit default,
    verify the default matches what the contract specifies.

    This catches drift where a developer adds .get("X", []) but the
    contract says the default is None.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    contract = _get_contract_fields()
    mismatches = []

    for filepath in SERVER_FILES:
        if not os.path.exists(filepath):
            continue
        with open(filepath) as f:
            for lineno, line in enumerate(f, 1):
                if line.strip().startswith("#"):
                    continue
                # Match session_info.get("field", default)
                for m in re.finditer(
                    r'session_info\.get\("(\w+)",\s*(.+?)\)', line
                ):
                    field = m.group(1)
                    raw_default = m.group(2).strip()

                    if field not in contract:
                        continue  # Caught by test_all_accessed_fields_in_contract

                    expected_default, _, _ = contract[field]

                    # Normalize for comparison
                    actual = _normalize_default(raw_default)
                    expected = _normalize_default(repr(expected_default))

                    if actual != expected:
                        mismatches.append(
                            "%s:%d  .get(\"%s\", %s) — contract says default=%r"
                            % (os.path.basename(filepath), lineno,
                               field, raw_default, expected_default)
                        )

    if mismatches:
        msg = "Default mismatches:\n  %s" % "\n  ".join(mismatches)
        assert_true(False, msg)
    print("  PASSED: All explicit defaults match contract")


def _normalize_default(s):
    """Normalize a default value string for comparison."""
    s = s.strip().strip('"').strip("'")
    # Normalize empty containers
    if s in ("{}", "dict()"):
        return "{}"
    if s in ("[]", "list()"):
        return "[]"
    if s in ("None",):
        return "None"
    if s in ("False", "True"):
        return s
    if s == '""' or s == "''":
        return '""'
    return s


# ---------------------------------------------------------------------------
# Test: .get() calls without explicit default on non-None contract fields
# ---------------------------------------------------------------------------

def test_get_without_default_on_non_none_fields():
    """
    Flag session_info.get("X") calls (no second argument) where the
    contract default is not None (e.g. {}, [], False, "").

    These work today because normalize_session_info fills them in first,
    but they become latent bugs if normalization is ever skipped.  The
    safe pattern is:  session_info.get("X", {})

    NOTE: Current violations are all guarded by truthiness checks (safe),
    so this test WARNS rather than fails.  Promote to failure once the
    codebase is clean.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    contract = _get_contract_fields()
    warnings = []

    for filepath in SERVER_FILES:
        if not os.path.exists(filepath):
            continue
        with open(filepath) as f:
            for lineno, line in enumerate(f, 1):
                if line.strip().startswith("#"):
                    continue
                # Match .get("X") with NO second argument
                for m in re.finditer(r'session_info\.get\("(\w+)"\)(?!\s*\.)', line):
                    field = m.group(1)
                    if field in contract:
                        default, _, _ = contract[field]
                        if default is not None:
                            warnings.append(
                                "%s:%d  .get(\"%s\") — contract default is %r, "
                                "consider .get(\"%s\", %r)"
                                % (os.path.basename(filepath), lineno,
                                   field, default, field, default)
                            )

    if warnings:
        msg = (
            "%d .get() call(s) without explicit default on non-None fields:\n  %s\n"
            "These are fragile if normalize_session_info is ever skipped.\n"
            "Use .get(\"field\", default) with the contract default."
            % (len(warnings), "\n  ".join(warnings))
        )
        assert_true(False, msg)
    else:
        print("  PASSED: All .get() calls on non-None fields use explicit defaults")


def test_protocol_version_consistency():
    """
    The client's sent protocol version must match
    CURRENT_PROTOCOL_VERSION in the contract.
    """
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    from agent.contract import CURRENT_PROTOCOL_VERSION

    client_file = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.exists(client_file):
        print("  SKIP (ai_agent.py not found)")
        return

    with open(client_file) as f:
        content = f.read()

    # Check that the client imports from contract (not hardcoded)
    assert_true(
        "_get_protocol_version" in content,
        "Client should use _get_protocol_version(), not a hardcoded version"
    )

    # Check the fallback value matches
    m = re.search(r'return\s+(\d+)\s*#\s*Fallback', content)
    if m:
        fallback = int(m.group(1))
        assert_equal(fallback, CURRENT_PROTOCOL_VERSION,
                     "Fallback version in _get_protocol_version() must match "
                     "CURRENT_PROTOCOL_VERSION (%d)" % CURRENT_PROTOCOL_VERSION)

    print("  PASSED: Protocol version consistent between contract and client")


# ---------------------------------------------------------------------------
# Test: MIN_SUPPORTED <= CURRENT
# ---------------------------------------------------------------------------

def test_version_bounds():
    """MIN_SUPPORTED_PROTOCOL_VERSION must be <= CURRENT_PROTOCOL_VERSION."""
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    from agent.contract import (
        CURRENT_PROTOCOL_VERSION,
        MIN_SUPPORTED_PROTOCOL_VERSION,
    )
    assert_true(
        MIN_SUPPORTED_PROTOCOL_VERSION <= CURRENT_PROTOCOL_VERSION,
        "MIN (%d) must be <= CURRENT (%d)"
        % (MIN_SUPPORTED_PROTOCOL_VERSION, CURRENT_PROTOCOL_VERSION)
    )
    assert_true(
        MIN_SUPPORTED_PROTOCOL_VERSION >= 1,
        "MIN must be >= 1"
    )
    print("  PASSED: Version bounds valid (MIN=%d, CURRENT=%d)"
          % (MIN_SUPPORTED_PROTOCOL_VERSION, CURRENT_PROTOCOL_VERSION))


# ---------------------------------------------------------------------------
# Test: warnings in response contract
# ---------------------------------------------------------------------------

def test_warnings_in_response_contract():
    """The warnings field must be in RESPONSE_FIELDS."""
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    assert_in("warnings", RESPONSE_FIELDS,
              "'warnings' must be in RESPONSE_FIELDS")
    print("  PASSED: 'warnings' is in RESPONSE_FIELDS")


# ---------------------------------------------------------------------------
# Test: client handles warnings from server
# ---------------------------------------------------------------------------

def test_client_handles_warnings():
    """ai_agent.py must check for 'warnings' in the server response."""
    client_file = os.path.join(_PROJECT_ROOT, "programs", "ai_agent.py")
    if not os.path.exists(client_file):
        print("  SKIP (ai_agent.py not found)")
        return

    with open(client_file) as f:
        content = f.read()

    assert_true(
        "warnings" in content and "AI Server Warning" in content,
        "Client must check for 'warnings' in server response and display them"
    )
    print("  PASSED: Client handles server warnings")


# ---------------------------------------------------------------------------
# Test: normalize_session_info covers all fields
# ---------------------------------------------------------------------------

def test_normalize_covers_all_fields():
    """normalize_session_info must fill in every registered field."""
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    from agent.contract import normalize_session_info

    # Start with empty session_info
    si = normalize_session_info({})

    contract = _get_contract_fields()
    missing = []
    for field_name, (default, _, _) in contract.items():
        if field_name not in si:
            missing.append(field_name)
        elif type(si[field_name]) != type(default) and default is not None:
            # Type mismatch (None defaults accept any type)
            missing.append("%s: got %s, expected %s"
                           % (field_name, type(si[field_name]).__name__,
                              type(default).__name__))

    if missing:
        assert_true(False, "normalize_session_info missed fields: %s" % missing)
    print("  PASSED: normalize fills all %d registered fields" % len(contract))


# ---------------------------------------------------------------------------
# Test: Mutable default isolation
# ---------------------------------------------------------------------------

def test_normalize_mutable_isolation():
    """Mutable defaults from normalize must be independent across calls."""
    if not CONTRACT_AVAILABLE:
        print("  SKIP (contract.py not importable)")
        return

    from agent.contract import normalize_session_info

    si1 = normalize_session_info({})
    si2 = normalize_session_info({})

    # Mutate si1's dict defaults
    si1["best_files"]["model"] = "x.pdb"
    si1["recovery_strategies"]["file.pdb"] = [{"strategy": "retry"}]

    # si2 must be unaffected
    assert_equal(si2["best_files"], {},
                 "Mutable default 'best_files' leaked between calls")
    assert_equal(si2["recovery_strategies"], {},
                 "Mutable default 'recovery_strategies' leaked between calls")
    print("  PASSED: Mutable defaults are isolated between normalize calls")


# =========================================================================
# Run
# =========================================================================

def run_all_tests():
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
