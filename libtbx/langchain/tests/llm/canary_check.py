"""v119.H3b: Startup canary orchestrator.

Combines two operational checks into a single operator-facing
tool that runs after deploy:

  1. agent_build pinning verification (no LLM).
     Reads canary_expected.json + live get_agent_build_info() and
     asserts they match.  Version mismatch -> hard fail; fingerprint
     drift -> soft warning (may be intentional mid-version retirement
     update).

  2. LLM provider reachability smoke test (one cheap call).
     Picks a provider (Google preferred, OpenAI fallback) and
     sends a minimal directive extraction probe through the same
     production code path the agent uses.  Asserts only that the
     LLM returned a parseable non-empty dict (per Gemini guardrail
     Q3).  Wrapped in a strict 30-second timeout via
     concurrent.futures.ThreadPoolExecutor to prevent hung calls
     from blocking deployment pipelines (Gemini's additional
     critical recommendation).

The LLM probe is the same as the "canary" Scenario registered in
tst_directive_extraction.build_scenarios() -- you can also run it
through the standard suite runner for diagnostic-level output:

    phenix.python tests/llm/run_llm_tests.py --scenario canary

Usage:
    phenix.python tests/llm/canary_check.py

Exit codes (per plan rev 2 §5.4, Tom's choice c -- graduated):
    0  All checks passed (or PASS with stderr fingerprint warning)
    1  Version mismatch -- wrong build deployed
    2  Configuration error -- missing canary_expected.json, no API key
    3  LLM probe failed -- API key set but call errored or timed out
"""
from __future__ import absolute_import, division, print_function

import concurrent.futures
import os
import sys
import time


# Make canary_utils (sibling of tests/) and framework (in tests/llm/)
# importable from this script's location.
_HERE = os.path.dirname(os.path.abspath(__file__))    # tests/llm/
_TESTS_DIR = os.path.dirname(_HERE)                   # tests/
for p in (_TESTS_DIR, _HERE):
    if p not in sys.path:
        sys.path.insert(0, p)


PROBE_TEXT = "Run phenix.refine with resolution 2.5"
# v119.H5.1: probe text includes one concrete settable parameter
# (resolution=2.5) so the validator preserves it.  An earlier
# version used "Run phenix.refine on the model", which OpenAI
# faithfully extracted as {"program_settings": {"phenix.refine": {}}}
# — semantically correct (no parameters mentioned, so none
# extracted) but stripped to {} by validate_directives because
# the inner dict is empty.  Google tended to over-extract and
# returned non-empty, hiding the issue.  Adding one concrete
# parameter makes the probe robust across providers without
# changing the canary's purpose (provider reachability /
# JSON-formatting sanity, not extraction quality).
# Timeout for the LLM probe.  The probe triggers intent
# classification PLUS the main directive extraction (two LLM
# round-trips), each subject to network variance.  5s was too
# tight in practice; 30s gives headroom for cold-network latency
# while still cleanly bounding hung calls (the actual failure
# mode this timeout protects against -- provider unreachable,
# DNS issue, etc.).
PROBE_TIMEOUT_S = 30.0


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _import_canary_utils():
    """Return (module, error_msg).  Never raises."""
    try:
        from libtbx.langchain.tests import canary_utils
        return canary_utils, None
    except ImportError:
        pass
    try:
        import canary_utils
        return canary_utils, None
    except ImportError as e:
        return None, str(e)


def _import_build_info():
    """Return (module, error_msg).  Never raises."""
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


def _import_framework_call_fn():
    """Return (call_fn, error_msg).  Never raises.

    The framework's call_directive_extractor is the
    production-faithful entry point used by the LLM test suite.
    Using it here ensures the canary probe exercises the same
    code path as the registered canary Scenario.
    """
    try:
        from framework import call_directive_extractor
        return call_directive_extractor, None
    except ImportError as e:
        return None, str(e)


def _pick_provider():
    """Return provider name or None per Tom's choice 3.

    Preference order: Google -> OpenAI.  Anthropic is excluded
    because the H1 RAG branch for Anthropic isn't wired through
    get_llm_and_embeddings; claiming Anthropic works here would
    mislead the operator.
    """
    if os.environ.get("GOOGLE_API_KEY"):
        return "google"
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"
    return None


def _check_agent_build(expected):
    """Verify deployed agent_build matches expected.

    Returns (version_ok, fingerprint_ok, info_dict, error_msg).
    error_msg is None on success; populated if the live build
    info could not be read at all.
    """
    bi, err = _import_build_info()
    if bi is None:
        return False, False, None, "could not import _build_info: %s" % err
    try:
        info = bi.get_agent_build_info()
    except Exception as e:
        return False, False, None, "get_agent_build_info raised: %s" % e
    version_ok = info.get("version") == expected["agent_version"]
    fingerprint_ok = (info.get("defaults_fingerprint")
                      == expected["defaults_fingerprint"])
    return version_ok, fingerprint_ok, info, None


def _run_probe(provider, call_fn):
    """Run the canary LLM probe with a strict 30-second timeout.

    Returns (passed, why, elapsed_s, error_or_none).

    The ThreadPoolExecutor context manager (`with`) ensures the
    executor is cleaned up cleanly on TimeoutError, per Gemini's
    minor implementation guardrail.  Note: on timeout, the worker
    thread that's actually waiting on the LLM HTTP call may still
    be alive in the background; Python doesn't support thread
    interruption.  The LLM SDK clients have their own internal
    HTTP timeouts (typically 60-120s), so the worker will
    eventually drain.  The user-visible result -- the canary
    fails fast with a clear message -- happens at the 30s mark
    regardless.

    Threshold rationale: the probe triggers TWO LLM round-trips
    (intent classification + main directive extraction).  Each
    round-trip has variable network latency (1-4s typical, up to
    10-15s under cold-network conditions).  30s gives sufficient
    headroom for normal latency variance while still bounding
    truly hung calls (provider unreachable, DNS issue) at a
    tolerable wait time.
    """
    def _do_call():
        return call_fn(PROBE_TEXT, provider)

    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
        future = ex.submit(_do_call)
        try:
            raw_output, parsed, error, _captured = future.result(
                timeout=PROBE_TIMEOUT_S)
        except concurrent.futures.TimeoutError:
            elapsed = time.time() - t0
            return (False,
                    "LLM call exceeded %.1fs timeout" % PROBE_TIMEOUT_S,
                    elapsed,
                    "timeout")
    elapsed = time.time() - t0
    if error is not None:
        return (False, "LLM call errored", elapsed, error)
    if not isinstance(parsed, dict):
        return (False,
                "expected dict, got %s" % type(parsed).__name__,
                elapsed, None)
    if not parsed:
        return (False,
                "expected non-empty dict; LLM returned {}",
                elapsed, None)
    return (True,
            "directives dict has %d keys" % len(parsed),
            elapsed, None)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print("PHENIX AI Agent -- Startup Canary")
    print("=" * 50)

    # (1) Load expected pins
    cu, cu_err = _import_canary_utils()
    if cu is None:
        print("CONFIG ERROR: could not import canary_utils: %s"
              % cu_err, file=sys.stderr)
        return 2

    expected, load_err = cu.load_canary_expected()
    if load_err is not None:
        print("CONFIG ERROR: %s" % load_err, file=sys.stderr)
        return 2

    # (2) Check agent_build
    version_ok, fingerprint_ok, info, build_err = _check_agent_build(expected)
    if build_err is not None:
        print("BUILD ERROR: %s" % build_err, file=sys.stderr)
        return 1

    print()
    print("Build metadata:")
    print("  expected version:     %s" % expected["agent_version"])
    print("  deployed version:     %s   %s"
          % (info.get("version", "?"),
             "[OK]" if version_ok else "[FAIL]"))
    print("  expected fingerprint: %s..."
          % expected["defaults_fingerprint"][:24])
    print("  deployed fingerprint: %s...   %s"
          % (info.get("defaults_fingerprint", "")[:24],
             "[OK]" if fingerprint_ok else "[WARN]"))
    print("  started_at:           %s" % info.get("started_at", "?"))

    if not version_ok:
        print()
        print("FAIL: deployed version does not match canary "
              "expected.  Wrong build deployed.", file=sys.stderr)
        return 1
    if not fingerprint_ok:
        print()
        print("WARN: fingerprint drift -- canary_expected pins "
              "%s but deployed reports %s.  May be intentional "
              "(mid-version retirement update)."
              % (expected["defaults_fingerprint"],
                 info.get("defaults_fingerprint", "")),
              file=sys.stderr)
        # Not a hard fail per Tom's choice (c).  Continue to LLM probe.

    # (3) LLM probe
    print()
    print("LLM probe:")
    provider = _pick_provider()
    if provider is None:
        print()
        print("CONFIG ERROR: no API key available.  Set "
              "GOOGLE_API_KEY or OPENAI_API_KEY and re-run.",
              file=sys.stderr)
        return 2

    call_fn, fn_err = _import_framework_call_fn()
    if call_fn is None:
        print()
        print("BUILD ERROR: could not import framework "
              "call_directive_extractor: %s" % fn_err,
              file=sys.stderr)
        return 1

    print("  selected provider:    %s" % provider)
    print("  probe:                %r" % PROBE_TEXT)
    print("  timeout:              %.1fs" % PROBE_TIMEOUT_S)

    passed, why, elapsed, error = _run_probe(provider, call_fn)
    status = "OK" if passed else "FAIL"
    print("  result:               %s (%s)" % (status, why))
    print("  elapsed:              %.2fs" % elapsed)

    if not passed:
        print()
        if error == "timeout":
            print("FAIL: LLM call timed out -- provider may be "
                  "experiencing degradation.", file=sys.stderr)
        else:
            print("FAIL: LLM probe did not return a parseable "
                  "non-empty directives dict.", file=sys.stderr)
            if error and error != "timeout":
                # Print only the first line of error to keep stderr clean
                first_line = error.split("\n")[0] if error else ""
                if first_line:
                    print("  error: %s" % first_line, file=sys.stderr)
        return 3

    print()
    print("Status: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
