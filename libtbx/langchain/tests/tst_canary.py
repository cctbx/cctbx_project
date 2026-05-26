"""K_H3a: Startup canary metadata regression tests.

v119.H3a.  Verifies the deployed server's agent_build matches
what tests/canary_expected.json pins.  Catches:
  - Wrong VERSION deployed (built from wrong branch)
  - Tables in core/llm.py edited without VERSION bump (fingerprint
    diverges from expected)
  - H2 injection broken (no agent_build in response)
  - Operator drift between VERSION file and canary_expected.json

Test organization (matches plan rev 2 §4.2):
  §4.2.1  Config-file tests          (2 tests)
  §4.2.2  Probe-mechanics tests      (3 tests)
  §4.2.3  Version-pinning tests      (2 tests)
  §4.2.4  Fingerprint-pinning tests  (2 tests)
  §4.2.5  End-to-end pinning         (1 test)

Total: 10 tests.

The probe sends an intentionally-malformed request through
phenix_ai.run_ai_agent.run() (missing required fields).  That
fails validate_request, flows through the error path
(_build_error_response), gets agent_build injected by H2's
post-processor, and returns.  No LLM call, no LangGraph traversal.

Failure severity per plan rev 2 §4.3 (graduated, Tom's choice c):
  - Hard fail (test failure):  version mismatch, malformed
    fingerprint, missing agent_build, missing config file
  - Soft warn (PASS + stderr): fingerprint drift when version
    matches (may be an intentional mid-version retirement update)
"""
from __future__ import absolute_import, division, print_function

import json
import os
import sys
import time


# ---------- Import / file-discovery helpers ------------------------

def _try_import_run_ai_agent():
    """Try to import run_ai_agent.run; return (run_fn, error_msg).

    Returns (callable, None) on success; (None, msg) on failure.
    Graceful skip pattern matches K_H1's _try_import for sandbox
    isolation.

    The production module path is `phenix.phenix_ai.run_ai_agent`
    (the phenix-project source tree).  We try the libtbx-style
    prefix as a fallback for forward compatibility but the canonical
    path is the first try.
    """
    try:
        from phenix.phenix_ai.run_ai_agent import run
        return run, None
    except ImportError:
        pass
    try:
        from libtbx.langchain.phenix_ai.run_ai_agent import run
        return run, None
    except ImportError:
        pass
    try:
        from phenix_ai.run_ai_agent import run
        return run, None
    except ImportError as e:
        return None, str(e)


def _try_import_canary_utils():
    """Try to import canary_utils; return (module, error_msg)."""
    try:
        from libtbx.langchain.tests import canary_utils
        return canary_utils, None
    except ImportError:
        pass
    try:
        # When run from tests/ directory (sibling import)
        here = os.path.dirname(os.path.abspath(__file__))
        if here not in sys.path:
            sys.path.insert(0, here)
        import canary_utils
        return canary_utils, None
    except ImportError as e:
        return None, str(e)


def _try_import_build_info():
    """Try to import core._build_info; return (module, error_msg)."""
    try:
        from libtbx.langchain.core import _build_info
        return _build_info, None
    except ImportError:
        pass
    try:
        from core import _build_info
        return _build_info, None
    except ImportError as e:
        return None, str(e)


def _version_file_path():
    """Return absolute path to VERSION at the langchain/ root."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "VERSION"))


def _read_version_file():
    """Return contents of VERSION file, or None if missing."""
    path = _version_file_path()
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            return f.read().strip()
    except Exception:
        return None


# ---------- Probe helpers ------------------------------------------

# The probe request shape (per plan rev 2 §4.1, Gemini Q2).
# CRITICAL detail: an "intentionally malformed" request is harder
# to construct than rev 1 anticipated.  REQUEST_V2_REQUIRED lists
# ["api_version", "files", "cycle_number"], but ALL of these have
# schema defaults that apply_request_defaults fills in before
# validate_request runs.  So omitting required fields does NOT
# trigger a validation failure.
#
# Instead, we send `files` with the wrong TYPE (string, not list).
# apply_request_defaults only fills in MISSING fields -- it does
# not repair existing-but-wrong-typed fields.  So validate_request
# catches the type mismatch and returns errors.  run() then routes
# through _build_error_response, which flows through
# _build_group_args_response and gets agent_build injected per H2.
#
# The "client_command" marker is operator-visible in server logs
# to label the malformed-on-purpose request (per Gemini Q2).
PROBE_REQUEST = {
    "api_version": "2.0",
    "files": "CANARY_PING_invalid_type",  # intentionally str, not list
    "cycle_number": 1,
    "client_command": "CANARY_PING",
}


def _send_probe():
    """Send the canary probe through run(); return parsed agent_build.

    Returns (response_dict, agent_build_dict, elapsed_s, error_msg).
    On failure, response_dict and agent_build_dict may be None;
    error_msg is populated.
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        return None, None, 0.0, "could not import run(): %s" % err

    t0 = time.time()
    try:
        result = run(api_version="2.0",
                     request_json=json.dumps(PROBE_REQUEST))
    except Exception as e:
        elapsed = time.time() - t0
        return None, None, elapsed, "run() raised: %s" % e
    elapsed = time.time() - t0

    # result is a group_args; we want response_json (the wire form)
    try:
        response_json = getattr(result, "response_json", None)
        if not response_json:
            return None, None, elapsed, "result has no response_json"
        response = json.loads(response_json)
    except Exception as e:
        return None, None, elapsed, "could not parse response_json: %s" % e

    agent_build = response.get("agent_build")
    return response, agent_build, elapsed, None


# =====================================================================
# §4.2.1: Config-file tests (2)
# =====================================================================

def test_canary_expected_exists():
    """canary_expected.json is present at the expected location."""
    cu, err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % err)
        return
    path = cu.canary_config_path()
    assert os.path.isfile(path), (
        "canary_expected.json not found at %s" % path)
    print("  PASS: test_canary_expected_exists")


def test_canary_expected_has_required_keys():
    """Config loads with both required keys, both str-typed."""
    cu, err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % err)
        return
    cfg, load_err = cu.load_canary_expected()
    assert load_err is None, (
        "load_canary_expected returned error: %s" % load_err)
    assert isinstance(cfg, dict)
    assert "agent_version" in cfg
    assert "defaults_fingerprint" in cfg
    assert isinstance(cfg["agent_version"], str)
    assert isinstance(cfg["defaults_fingerprint"], str)
    print("  PASS: test_canary_expected_has_required_keys")


# =====================================================================
# §4.2.2: Probe-mechanics tests (3)
# =====================================================================

def test_invalid_request_produces_error_response():
    """The probe request goes through run() and returns an error response."""
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err
    assert response is not None
    # The error path populates response["error"] with the validation
    # failure message
    error_field = response.get("error")
    assert error_field, (
        "Expected response.error to be populated for malformed "
        "request; got %r" % error_field)
    print("  PASS: test_invalid_request_produces_error_response")


def test_error_response_contains_agent_build():
    """The error response has agent_build with all three subfields."""
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    _response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err
    assert isinstance(agent_build, dict), (
        "agent_build missing or not a dict: %r" % agent_build)
    for key in ("version", "defaults_fingerprint", "started_at"):
        assert key in agent_build, (
            "agent_build missing subfield %r; got %r"
            % (key, agent_build))
    print("  PASS: test_error_response_contains_agent_build")


def test_error_response_is_fast():
    """Probe completes quickly; slow probe means an LLM call snuck in.

    Threshold: 3 seconds, generous for sandbox/PHENIX environment
    variability (langchain cold-start imports, heavily-loaded build
    machines, etc.).  An LLM call would take 5-30 seconds, so 3s
    still cleanly detects an accidental LLM leak into the error path.
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    _response, _agent_build, elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err
    assert elapsed < 3.0, (
        "Probe took %.2fs; expected <3.0s.  Possible LLM call "
        "leaked into the error path." % elapsed)
    print("  PASS: test_error_response_is_fast (%.3fs)" % elapsed)


# =====================================================================
# §4.2.3: Version-pinning tests (2)
# =====================================================================

def test_agent_version_matches_canary_expected():
    """Deployed agent_build.version equals canary_expected agent_version.

    HARD FAIL on mismatch (wrong build deployed).
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    cu, cu_err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % cu_err)
        return
    cfg, load_err = cu.load_canary_expected()
    assert load_err is None, "config load failed: %s" % load_err

    _response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err

    expected = cfg["agent_version"]
    deployed = agent_build.get("version")
    assert deployed == expected, (
        "VERSION MISMATCH: canary_expected.json pins agent_version=%r "
        "but deployed server reports version=%r.  Either the wrong "
        "build was deployed, or canary_expected.json is stale."
        % (expected, deployed))
    print("  PASS: test_agent_version_matches_canary_expected (%s)" % deployed)


def test_version_file_matches_canary_expected():
    """VERSION file content matches canary_expected agent_version.

    Catches the case where the operator edited canary_expected.json
    but forgot to update VERSION (or vice versa) -- visible before
    the server is started.
    """
    cu, err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % err)
        return
    cfg, load_err = cu.load_canary_expected()
    assert load_err is None, "config load failed: %s" % load_err

    version = _read_version_file()
    assert version is not None, (
        "VERSION file not found at %s" % _version_file_path())
    expected = cfg["agent_version"]
    assert version == expected, (
        "VERSION/canary drift: VERSION file contains %r but "
        "canary_expected.json pins agent_version=%r.  Update both "
        "files together." % (version, expected))
    print("  PASS: test_version_file_matches_canary_expected (%s)" % version)


# =====================================================================
# §4.2.4: Fingerprint-pinning tests (2)
# =====================================================================

def test_fingerprint_matches_canary_expected():
    """Deployed fingerprint matches canary_expected.  SOFT WARN on drift.

    Fingerprint drift can be intentional (e.g., a model retirement
    landed mid-version without bumping VERSION).  This test prints
    a clear warning to stderr but does NOT fail -- the K-suite
    stays green.  Operators reviewing logs see the warning.
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    cu, cu_err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % cu_err)
        return
    cfg, load_err = cu.load_canary_expected()
    assert load_err is None, "config load failed: %s" % load_err

    _response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err

    expected = cfg["defaults_fingerprint"]
    deployed = agent_build.get("defaults_fingerprint", "")
    if deployed == expected:
        print("  PASS: test_fingerprint_matches_canary_expected")
    else:
        # Soft warning -- print but do not fail
        print("  WARN: fingerprint drift -- canary_expected pins %s "
              "but deployed reports %s.  May be intentional "
              "(mid-version retirement update)."
              % (expected, deployed),
              file=sys.stderr)
        print("  PASS: test_fingerprint_matches_canary_expected "
              "(soft-warn: fingerprint drift)")


def test_fingerprint_is_well_formed():
    """Deployed fingerprint starts with 'sha256:' and has 64 hex chars.

    HARD FAIL on malformed fingerprint (unambiguously broken).
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    _response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err
    fingerprint = agent_build.get("defaults_fingerprint", "")
    assert fingerprint.startswith("sha256:"), (
        "fingerprint must start with 'sha256:', got %r"
        % fingerprint[:32])
    hex_part = fingerprint[len("sha256:"):]
    assert len(hex_part) == 64, (
        "fingerprint hex part must be 64 chars, got %d (%r)"
        % (len(hex_part), hex_part))
    import re
    assert re.match(r"^[0-9a-f]{64}$", hex_part), (
        "fingerprint hex must be all lowercase hex, got %r" % hex_part)
    print("  PASS: test_fingerprint_is_well_formed")


# =====================================================================
# §4.2.5: End-to-end pinning (1)
# =====================================================================

def test_canary_pin_full_consistency():
    """End-to-end: load expected, send probe, assert version + fingerprint shape.

    The "if this passes, the deploy is the right build" assertion.
    """
    run, err = _try_import_run_ai_agent()
    if run is None:
        print("  SKIP: cannot import run_ai_agent.run (%s)" % err)
        return
    cu, cu_err = _try_import_canary_utils()
    if cu is None:
        print("  SKIP: cannot import canary_utils (%s)" % cu_err)
        return
    cfg, load_err = cu.load_canary_expected()
    assert load_err is None, "config load failed: %s" % load_err

    _response, agent_build, _elapsed, probe_err = _send_probe()
    assert probe_err is None, "probe failed: %s" % probe_err

    # Version must match exactly
    assert agent_build["version"] == cfg["agent_version"], (
        "version pin failed: deployed=%r expected=%r"
        % (agent_build["version"], cfg["agent_version"]))
    # Fingerprint must be well-formed (may or may not match,
    # see soft-warn test above)
    fp = agent_build["defaults_fingerprint"]
    assert fp.startswith("sha256:"), (
        "deployed fingerprint malformed: %r" % fp)
    assert len(fp) == len("sha256:") + 64, (
        "deployed fingerprint wrong length: %r" % fp)
    # And the VERSION file must agree (the canary_expected.json /
    # VERSION lockstep enforced in
    # test_version_file_matches_canary_expected)
    version = _read_version_file()
    assert version == cfg["agent_version"], (
        "VERSION/canary lockstep broken: VERSION=%r expected=%r"
        % (version, cfg["agent_version"]))

    print("  PASS: test_canary_pin_full_consistency")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §4.2.1: Config-file (2)
    test_canary_expected_exists()
    test_canary_expected_has_required_keys()
    # §4.2.2: Probe mechanics (3)
    test_invalid_request_produces_error_response()
    test_error_response_contains_agent_build()
    test_error_response_is_fast()
    # §4.2.3: Version pinning (2)
    test_agent_version_matches_canary_expected()
    test_version_file_matches_canary_expected()
    # §4.2.4: Fingerprint pinning (2)
    test_fingerprint_matches_canary_expected()
    test_fingerprint_is_well_formed()
    # §4.2.5: End-to-end (1)
    test_canary_pin_full_consistency()


if __name__ == "__main__":
    run_all_tests()
