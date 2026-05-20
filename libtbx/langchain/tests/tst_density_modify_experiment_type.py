"""Tests for v118.9: experiment-type-conditional program canonicalization.

The v118.9 plan adds a code-side validator that corrects after_program
when the LLM directive extractor picks a program that's canonical for
the wrong experiment type (e.g., autobuild_denmod for cryo-EM data
instead of resolve_cryo_em).

This file holds the K_DENMOD test suite (13 tests) that locks the
validator's behavior:

  Bug case (the reason the validator exists):
    K1   — autobuild_denmod picked for cryo-EM input → corrected
            to resolve_cryo_em
    K10  — full preprocessor-style advice with cryo-EM signals
            (the actual Tom-bug input)

  Mirror case (symmetry insurance):
    K2   — resolve_cryo_em picked for X-ray input → corrected
            to autobuild_denmod

  No-change cases (validator must be a no-op):
    K3   — autobuild_denmod for X-ray (LLM got it right)
    K4   — resolve_cryo_em for cryo-EM (LLM got it right)
    K5   — autobuild_denmod for ambiguous experiment type
    K6   — phenix.refine (not in reprints table)
    K7   — empty stop_conditions
    K8   — stop_conditions without after_program
    K9   — raw_advice=None (single-input path)

  Ordering/placement verification (Gemini Gap A response):
    K11  — Tom's case path: grounding bypass + correction
    K12  — Alternate path: no bypass + grounding drops bad name
            (validator must not "resurrect" the dropped name)

  Diagnostic transparency (Gemini Gap C response):
    K13  — Correction emits [DIRECTIVE_CORRECTION] log marker AND
            stderr line AND _corrected_from sidecar field

All tests use the mock-LLM pattern from
tst_after_program_fill_from_raw.py: monkey-patch _call_llm to return
canned JSON, then call extract_directives and assert on the result.
This makes the suite deterministic and runs without API access.
"""

from __future__ import absolute_import, division, print_function

import io
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOTS = (
    os.path.abspath(os.path.join(_HERE, "..")),
    os.path.abspath(_HERE),
)
for _root in _ROOTS:
    if _root not in sys.path:
        sys.path.insert(0, _root)

for mod in ("libtbx", "libtbx.langchain", "libtbx.langchain.agent"):
    if mod not in sys.modules:
        sys.modules[mod] = types.ModuleType(mod)

if "agent.program_registry" not in sys.modules:
    pr = types.ModuleType("agent.program_registry")
    class _Registry:
        def __init__(self, *a, **k): pass
        def get(self, *a, **k): return None
        def all_programs(self): return []
        def get_all_program_keys(self): return []
    pr.ProgramRegistry = _Registry
    sys.modules["agent.program_registry"] = pr
    sys.modules["libtbx.langchain.agent.program_registry"] = pr

if "agent.intent_classifier" not in sys.modules:
    intent_mod = types.ModuleType("agent.intent_classifier")
    def _classify_intent_stub(text):
        return {"intent": "tutorial", "confidence": "medium"}
    intent_mod.classify_intent = _classify_intent_stub
    sys.modules["agent.intent_classifier"] = intent_mod
    sys.modules["libtbx.langchain.agent.intent_classifier"] = intent_mod

try:
    from libtbx.langchain.agent import directive_extractor as de
except ImportError:
    from agent import directive_extractor as de


# =====================================================================
# Test inputs
# =====================================================================

# Preprocessor-shaped advice matching the format Tom's preprocessor
# produces.  Includes the explicit Experiment Type field and the
# half-map file extensions that the validator's regex picks up.
_PROCESSED_CRYOEM = (
    "1. **Input Files Found**: 7mjs_23883_H_1.ccp4, "
    "7mjs_23883_H_2.ccp4, 7mjs_23883_H.fa\n\n"
    "2. **Experiment Type**: cryo-EM (inferred from half-maps)\n\n"
    "3. **Primary Goal**: Perform density modification.\n\n"
    "4. **Key Parameters**:\nNone\n\n"
    "5. **Program Parameters**:\nNone\n\n"
    "6. **Special Instructions**:\nNone\n\n"
    "7. **Stop Condition**: None"
)

_PROCESSED_XRAY = (
    "1. **Input Files Found**: model.pdb, data.mtz\n\n"
    "2. **Experiment Type**: X-ray crystallography\n\n"
    "3. **Primary Goal**: Perform density modification.\n\n"
    "4. **Key Parameters**:\nNone\n\n"
    "5. **Program Parameters**:\nNone\n\n"
    "6. **Special Instructions**:\nNone\n\n"
    "7. **Stop Condition**: None"
)

# Ambiguous: no experiment-type signal, no specific file extensions
_PROCESSED_AMBIGUOUS = (
    "1. **Input Files Found**: data.txt\n\n"
    "2. **Primary Goal**: Perform some operation.\n\n"
    "3. **Stop Condition**: None"
)


def _with_mock_llm(response_str):
    """Replace _call_llm with one that returns a fixed string."""
    def _fake(prompt, provider, model, log):
        return response_str
    de._call_llm = _fake


def _quiet_log(_msg):
    """Default no-op log; replaced when a test inspects log output."""
    pass


# =====================================================================
# K1 — The bug case: autobuild_denmod for cryo-EM gets corrected
# =====================================================================

def test_K1_autobuild_denmod_for_cryoem_corrected():
    """K1: LLM picks phenix.autobuild_denmod for cryo-EM input.
    Validator corrects to phenix.resolve_cryo_em.

    This is the production bug case (Tom's runs 237/239)."""
    print("Test: K1_autobuild_denmod_for_cryoem_corrected")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod", '
        '"skip_validation": true}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.resolve_cryo_em", (
        "Expected after_program=phenix.resolve_cryo_em (corrected from "
        "autobuild_denmod), got %r" % ap)
    print("  PASS")


# =====================================================================
# K2 — The mirror case: resolve_cryo_em for X-ray gets corrected
# =====================================================================

def test_K2_resolve_cryo_em_for_xray_corrected():
    """K2: LLM picks phenix.resolve_cryo_em for X-ray input.
    Validator corrects to phenix.autobuild_denmod.

    Mirror case of K1.  Not yet observed in production but the
    validator handles both directions symmetrically."""
    print("Test: K2_resolve_cryo_em_for_xray_corrected")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.resolve_cryo_em", '
        '"skip_validation": true}}')
    result = de.extract_directives(
        _PROCESSED_XRAY, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.autobuild_denmod", (
        "Expected after_program=phenix.autobuild_denmod (corrected from "
        "resolve_cryo_em), got %r" % ap)
    print("  PASS")


# =====================================================================
# K3 / K4 — LLM picked correctly: no change
# =====================================================================

def test_K3_autobuild_denmod_for_xray_unchanged():
    """K3: LLM picks phenix.autobuild_denmod for X-ray input.
    Validator is no-op (LLM got it right)."""
    print("Test: K3_autobuild_denmod_for_xray_unchanged")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_XRAY, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.autobuild_denmod", (
        "Expected after_program=phenix.autobuild_denmod (no change), "
        "got %r" % ap)
    assert "_corrected_from" not in sc, (
        "Expected no _corrected_from sidecar (no correction occurred), "
        "got %r" % sc.get("_corrected_from"))
    print("  PASS")


def test_K4_resolve_cryo_em_for_cryoem_unchanged():
    """K4: LLM picks phenix.resolve_cryo_em for cryo-EM input.
    Validator is no-op (LLM got it right)."""
    print("Test: K4_resolve_cryo_em_for_cryoem_unchanged")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.resolve_cryo_em"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.resolve_cryo_em", (
        "Expected after_program=phenix.resolve_cryo_em (no change), "
        "got %r" % ap)
    assert "_corrected_from" not in sc, (
        "Expected no _corrected_from sidecar, got %r"
        % sc.get("_corrected_from"))
    print("  PASS")


# =====================================================================
# K5 — Ambiguous experiment type: no correction
# =====================================================================

def test_K5_ambiguous_experiment_type_no_correction():
    """K5: Ambiguous advice (no cryo-EM and no X-ray signals).
    Validator declines to act."""
    print("Test: K5_ambiguous_experiment_type_no_correction")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_AMBIGUOUS, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    # Validator declines on ambiguous experiment.  Whether the
    # after_program survives downstream is up to other passes; we
    # only assert the validator didn't rewrite it.
    if ap is not None:
        assert ap == "phenix.autobuild_denmod", (
            "Expected after_program unchanged (validator declined on "
            "ambiguous experiment), got %r" % ap)
    assert "_corrected_from" not in sc, (
        "Expected no _corrected_from sidecar on ambiguous case, got %r"
        % sc.get("_corrected_from"))
    print("  PASS")


# =====================================================================
# K6 — Other after_program values: no change
# =====================================================================

def test_K6_other_after_program_untouched():
    """K6: phenix.refine (not in PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE).
    Validator is no-op."""
    print("Test: K6_other_after_program_untouched")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.refine"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="refine and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.refine", (
        "Expected after_program=phenix.refine unchanged, got %r" % ap)
    assert "_corrected_from" not in sc
    print("  PASS")


# =====================================================================
# K7 / K8 — Empty / missing stop_conditions cases
# =====================================================================

def test_K7_empty_stop_conditions_no_crash():
    """K7: LLM returns empty dict.  Validator must not crash."""
    print("Test: K7_empty_stop_conditions_no_crash")
    _with_mock_llm('{}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    # No assertion on result content; just verifying no exception
    # was raised.  Validator no-ops on missing stop_conditions.
    assert isinstance(result, dict)
    print("  PASS")


def test_K8_missing_after_program_no_change():
    """K8: stop_conditions present but no after_program field.
    Validator is no-op."""
    print("Test: K8_missing_after_program_no_change")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    # No after_program to correct — and our validator should not
    # have invented one.  (Other code paths may fill it in from raw
    # advice; that's not the validator's concern.)
    # We assert that the validator itself did not record a correction.
    assert "_corrected_from" not in sc, (
        "Validator should not record correction when no after_program "
        "to correct, got %r" % sc.get("_corrected_from"))
    print("  PASS")


# =====================================================================
# K9 — raw_advice=None (single-input path)
# =====================================================================

def test_K9_raw_advice_none_uses_user_advice():
    """K9: raw_advice is None (single-input path).  Validator
    still corrects using user_advice alone."""
    print("Test: K9_raw_advice_none_uses_user_advice")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice=None,
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.resolve_cryo_em", (
        "Expected after_program=phenix.resolve_cryo_em (corrected via "
        "user_advice signals), got %r" % ap)
    print("  PASS")


# =====================================================================
# K10 — End-to-end with Tom's actual bug case
# =====================================================================

def test_K10_end_to_end_density_modify_and_stop():
    """K10: The actual Tom-bug input — preprocessed advice exactly
    matching what the production preprocessor emits for
    'density modify and stop' on cryo-EM half-maps."""
    print("Test: K10_end_to_end_density_modify_and_stop")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod", '
        '"skip_validation": true}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.resolve_cryo_em", (
        "Expected after_program=phenix.resolve_cryo_em, got %r"
        % sc.get("after_program"))
    assert sc.get("stop_after_requested") is True, (
        "Expected stop_after_requested=True preserved, got %r"
        % sc.get("stop_after_requested"))
    assert sc.get("skip_validation") is True, (
        "Expected skip_validation=True preserved, got %r"
        % sc.get("skip_validation"))
    print("  PASS")


# =====================================================================
# K11 — Grounding-bypass interaction (Tom's path, Gemini Gap A)
# =====================================================================

def test_K11_grounding_bypass_then_correction():
    """K11: stop_after_requested=True (triggers grounding bypass)
    + after_program=autobuild_denmod + cryo-EM signals.

    The grounding check is bypassed because the LLM asserted
    user-stop intent (see _validate_after_program_grounded:1697).
    Our validator then runs and corrects the wrong program.

    This is Tom's actual path.  The test verifies that the
    bypass and the correction work together end-to-end."""
    print("Test: K11_grounding_bypass_then_correction")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.resolve_cryo_em", (
        "Expected after_program corrected to phenix.resolve_cryo_em "
        "(grounding bypass let the wrong name survive for us to "
        "correct), got %r" % sc.get("after_program"))
    # The _corrected_from sidecar should be present with the
    # documented schema
    cf = sc.get("_corrected_from")
    assert cf is not None, (
        "Expected _corrected_from sidecar to be present after "
        "correction, got %r" % cf)
    assert cf.get("from") == "phenix.autobuild_denmod"
    assert cf.get("to") == "phenix.resolve_cryo_em"
    assert cf.get("reason") == "experiment_type_mismatch"
    assert cf.get("experiment_type") == "cryoem"
    print("  PASS")


# =====================================================================
# K12 — No grounding bypass + ungrounded name (Gemini Gap A)
# =====================================================================

def test_K12_no_bypass_grounding_drops_bad_name():
    """K12: stop_after_requested=False (no bypass) +
    after_program=autobuild_denmod + cryo-EM signals BUT the
    advice does NOT contain 'autobuild_denmod' as a literal.

    The grounding check fires (because no stop_after_requested
    bypass).  autobuild_denmod isn't in the advice text, so
    grounding drops it as fabrication.  Our validator then has
    no after_program to correct.

    Verifies that the validator placement (AFTER grounding) is
    intentional: when grounding drops the wrong name, our
    validator should NOT inappropriately resurrect it."""
    print("Test: K12_no_bypass_grounding_drops_bad_name")
    _with_mock_llm(
        '{"stop_conditions": {"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify",  # no "stop" trigger
    )
    sc = result.get("stop_conditions", {})
    # Grounding should have dropped after_program because:
    #  - stop_after_requested is not True (no bypass)
    #  - "autobuild_denmod" doesn't appear in user_advice or raw_advice
    ap = sc.get("after_program")
    # Either the stop_conditions block is gone entirely OR
    # after_program is removed.  Both indicate grounding's drop
    # was respected.
    assert ap is None, (
        "Expected after_program dropped by grounding (the wrong "
        "name should not have been resurrected by the validator), "
        "got %r" % ap)
    # _corrected_from should NOT be set — there was nothing to
    # correct.
    assert "_corrected_from" not in sc, (
        "Validator should not record correction when grounding "
        "already dropped the after_program, got %r"
        % sc.get("_corrected_from"))
    print("  PASS")


# =====================================================================
# K13 — Diagnostic emission (Gemini Gap C)
# =====================================================================

def test_K13_correction_emits_diagnostic_markers():
    """K13: When a correction is made, verify ALL three diagnostic
    channels emit the [DIRECTIVE_CORRECTION] marker:
      - log_func receives the message
      - stderr receives the message
      - directives.stop_conditions._corrected_from sidecar is set

    Gemini Gap C: silent mutations are bad.  Every correction must
    be transparently visible."""
    print("Test: K13_correction_emits_diagnostic_markers")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')

    # Capture log messages
    log_messages = []
    def _capture_log(msg):
        log_messages.append(msg)

    # Capture stderr
    real_stderr = sys.stderr
    captured_stderr = io.StringIO()
    sys.stderr = captured_stderr
    try:
        result = de.extract_directives(
            _PROCESSED_CRYOEM, provider="google",
            log_func=_capture_log,
            raw_advice="density modify and stop",
        )
    finally:
        sys.stderr = real_stderr

    stderr_text = captured_stderr.getvalue()

    # Channel 1: log_func received the message
    correction_logs = [m for m in log_messages
                       if "[DIRECTIVE_CORRECTION]" in m]
    assert len(correction_logs) == 1, (
        "Expected exactly one [DIRECTIVE_CORRECTION] log message, "
        "got %d: %r" % (len(correction_logs), correction_logs))
    log_msg = correction_logs[0]
    assert "phenix.autobuild_denmod" in log_msg, (
        "Expected log to mention source program, got %r" % log_msg)
    assert "phenix.resolve_cryo_em" in log_msg, (
        "Expected log to mention target program, got %r" % log_msg)
    assert "cryoem" in log_msg, (
        "Expected log to mention experiment type, got %r" % log_msg)

    # Channel 2: stderr received the message
    assert "[DIRECTIVE_CORRECTION]" in stderr_text, (
        "Expected stderr to contain [DIRECTIVE_CORRECTION] marker, "
        "got %r" % stderr_text[:500])
    assert "phenix.autobuild_denmod" in stderr_text
    assert "phenix.resolve_cryo_em" in stderr_text

    # Channel 3: sidecar field
    sc = result.get("stop_conditions", {})
    cf = sc.get("_corrected_from")
    assert cf is not None, (
        "Expected _corrected_from sidecar, got None")
    # Individual key assertions (not strict dict equality) so future
    # additions to the sidecar (e.g., a timestamp or validator_version)
    # don't break this test.  We only assert the documented fields
    # are present with correct values.
    assert cf.get("from") == "phenix.autobuild_denmod", (
        "Expected sidecar from=phenix.autobuild_denmod, got %r"
        % cf.get("from"))
    assert cf.get("to") == "phenix.resolve_cryo_em", (
        "Expected sidecar to=phenix.resolve_cryo_em, got %r"
        % cf.get("to"))
    assert cf.get("reason") == "experiment_type_mismatch", (
        "Expected sidecar reason=experiment_type_mismatch, got %r"
        % cf.get("reason"))
    assert cf.get("experiment_type") == "cryoem", (
        "Expected sidecar experiment_type=cryoem, got %r"
        % cf.get("experiment_type"))

    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    tests = [
        test_K1_autobuild_denmod_for_cryoem_corrected,
        test_K2_resolve_cryo_em_for_xray_corrected,
        test_K3_autobuild_denmod_for_xray_unchanged,
        test_K4_resolve_cryo_em_for_cryoem_unchanged,
        test_K5_ambiguous_experiment_type_no_correction,
        test_K6_other_after_program_untouched,
        test_K7_empty_stop_conditions_no_crash,
        test_K8_missing_after_program_no_change,
        test_K9_raw_advice_none_uses_user_advice,
        test_K10_end_to_end_density_modify_and_stop,
        test_K11_grounding_bypass_then_correction,
        test_K12_no_bypass_grounding_drops_bad_name,
        test_K13_correction_emits_diagnostic_markers,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            print("  ERROR: %s: %s" % (type(e).__name__, e))
            failed += 1
    print("")
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        raise AssertionError("%d K_DENMOD tests failed" % failed)


if __name__ == "__main__":
    run_all_tests()
