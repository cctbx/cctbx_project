"""Tests for v118.9 + v119.H18: experiment-type-conditional program
canonicalization.

The v118.9 plan adds a code-side validator that corrects after_program
when the LLM directive extractor picks a program that's canonical for
the wrong experiment type (e.g., autobuild_denmod for cryo-EM data
instead of resolve_cryo_em).

v119.H18 extends the validator with file-extension-based detection
as the PRIMARY signal (text-based detection remains as fallback).
Triggered by the AF_7mjs regression: user wrote "density modify and
stop" with cryo-EM half-map inputs, the preprocessor produced terse
output with no experiment-type marker, and v118.9's text-only
detector returned None ("decline to act"), leaving the wrong
program intact.

This file holds the K_DENMOD test suite (20 tests) that locks the
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

  H18: file-extension-based primary signal:
    K14  — AF_7mjs failure verbatim: terse advice + cryo-EM files
            → corrected via files (source=files)
    K15  — Mirror: terse advice + X-ray files → corrected via files
    K16  — Backward compat: cryo-EM advice + NO files → text
            fallback still works (source=text)
    K17  — Pre-H18 bug-path: terse advice + NO files → declines
            (no signal, no correction)
    K18  — Files-win on conflict: cryo-EM files + advice says
            X-ray → files win, OVERRIDDEN logged
    K19  — Mixed-input drift (Pitfall 1): both .mtz and .ccp4 in
            inventory → defers to text, [..._MIXED] logged
    K20  — plan_generator uses shared helper (single source of
            truth refactor)

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
# H18 input fixtures
# =====================================================================
# These mimic the actual preprocessor output for the AF_7mjs failure
# case: terse advice with NO explicit experiment-type marker.  The
# pre-H18 detector returns None for these inputs (no cryo-EM / no
# X-ray tokens).  H18 fixes the case by inspecting the file inventory
# directly.

_PROCESSED_TERSE = (
    "Primary Goal: Perform density modification and stop.\n\n"
    "Stop Condition: stop after density modification completes"
)

# Cryo-EM file inventory (matches AF_7mjs production run):
_AF7MJS_FILES = [
    "/path/7mjs_23883_H_1.ccp4",
    "/path/7mjs_23883_H_2.ccp4",
    "/path/sub_02_resolve_cryo_em/rcm_0/initial_map.ccp4",
    "/path/sub_02_resolve_cryo_em/rcm_0/denmod_map.ccp4",
    "/path/7mjs_23883_H.fa",
]

# Pure X-ray file inventory (mirror case):
_XRAY_FILES = [
    "/path/data.mtz",
    "/path/model.pdb",
    "/path/sequence.fa",
]

# Mixed inventory (Pitfall 1: dirty-directory accumulation):
_MIXED_FILES = [
    "/path/old_data.mtz",      # from an earlier cycle
    "/path/current_map.ccp4",  # current cryo-EM input
]


# =====================================================================
# K14 — AF_7mjs failure case (file-based primary signal)
# =====================================================================

def test_K14_h18_af7mjs_terse_advice_cryoem_files():
    """K14: The AF_7mjs failure verbatim — terse advice with no
    experiment-type tokens, plus cryo-EM file inventory.

    Pre-H18: text-only detector returns None ("ambiguous"); no
    correction; bug persists (after_program stays autobuild_denmod
    and the LLM downstream picks predict_and_build, overriding
    user intent).

    H18: file-based detector identifies cryo-EM from .ccp4
    extensions; correction fires; after_program → resolve_cryo_em
    with source=files in the sidecar."""
    print("Test: K14_h18_af7mjs_terse_advice_cryoem_files")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod", '
        '"skip_validation": true}}')
    result = de.extract_directives(
        _PROCESSED_TERSE, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
        original_files=_AF7MJS_FILES,
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.resolve_cryo_em", (
        "Expected H18 file-based correction; got %r"
        % sc.get("after_program"))
    cf = sc.get("_corrected_from") or {}
    assert cf.get("source") == "files", (
        "Expected source=files in sidecar; got %r" % cf.get("source"))
    assert cf.get("experiment_type") == "cryoem"
    ev = cf.get("evidence") or {}
    assert ".ccp4" in ev.get("cryoem_exts", []), (
        "Expected .ccp4 in evidence; got %r" % ev)
    print("  PASS")


# =====================================================================
# K15 — Mirror: terse advice + X-ray files
# =====================================================================

def test_K15_h18_mirror_terse_advice_xray_files():
    """K15: Symmetric mirror.  LLM picks resolve_cryo_em for X-ray
    data; file-based detector identifies xray from .mtz; correction
    flips to autobuild_denmod."""
    print("Test: K15_h18_mirror_terse_advice_xray_files")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.resolve_cryo_em", '
        '"skip_validation": true}}')
    result = de.extract_directives(
        _PROCESSED_TERSE, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
        original_files=_XRAY_FILES,
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.autobuild_denmod"
    cf = sc.get("_corrected_from") or {}
    assert cf.get("source") == "files"
    assert cf.get("experiment_type") == "xray"
    print("  PASS")


# =====================================================================
# K16 — Backward compat: no files → text fallback still works
# =====================================================================

def test_K16_h18_backward_compat_text_fallback():
    """K16: When original_files is not provided, the validator must
    fall back to text-only detection — the pre-H18 behavior.

    This protects K1-K13 from regression.  Here we re-prove that
    the text path still corrects when given preprocessor-shaped
    advice that contains the cryo-EM markers.  Equivalent to K1 but
    asserts source=text in the sidecar."""
    print("Test: K16_h18_backward_compat_text_fallback")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
        # no original_files
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.resolve_cryo_em"
    cf = sc.get("_corrected_from") or {}
    assert cf.get("source") == "text", (
        "Expected source=text when no files passed; got %r"
        % cf.get("source"))
    print("  PASS")


# =====================================================================
# K17 — Pre-H18 bug path: terse advice + no files → declines
# =====================================================================

def test_K17_h18_no_signal_declines():
    """K17: Demonstrates that H18's improvement REQUIRES files to be
    threaded through.  Without files AND without text experiment-
    type markers, the validator correctly declines (the pre-H18 bug
    path is preserved when files aren't supplied).

    Why this matters: the LLM-pipeline integration must pass
    original_files for the AF_7mjs fix to activate.  This test
    pins that contract."""
    print("Test: K17_h18_no_signal_declines")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')
    result = de.extract_directives(
        _PROCESSED_TERSE, provider="google",
        log_func=_quiet_log,
        raw_advice="density modify and stop",
        original_files=None,  # explicit None
    )
    sc = result.get("stop_conditions", {})
    # No correction → after_program unchanged
    assert sc.get("after_program") == "phenix.autobuild_denmod", (
        "Expected no correction without files+text signal; got %r"
        % sc.get("after_program"))
    assert "_corrected_from" not in sc
    print("  PASS")


# =====================================================================
# K18 — Files-win on conflict (with OVERRIDDEN telemetry)
# =====================================================================

def test_K18_h18_files_win_on_conflict():
    """K18: When files and text disagree, files win (per Gemini's
    H18 review policy).  Telemetry requirement: log must include
    OVERRIDDEN annotation for audit.

    Justification: files are the hard physical boundary — passing
    .ccp4 to an X-ray-only program crashes regardless of what the
    user wrote.  Text-based detection has known false-positive
    vectors ('mad' in 'modify', 'sad' in arbitrary sentences)."""
    print("Test: K18_h18_files_win_on_conflict")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')

    log_messages = []
    def _capture_log(msg):
        log_messages.append(msg)

    # User text says X-ray (preprocessor marker); files are cryo-EM
    result = de.extract_directives(
        _PROCESSED_XRAY, provider="google",
        log_func=_capture_log,
        raw_advice="this is x-ray data, please run density modification",
        original_files=["/p/map.ccp4"],
    )
    sc = result.get("stop_conditions", {})
    # Files won: correction fired in the cryo-EM direction
    assert sc.get("after_program") == "phenix.resolve_cryo_em", (
        "Files should override text; got %r"
        % sc.get("after_program"))
    cf = sc.get("_corrected_from") or {}
    assert cf.get("source") == "files"
    assert cf.get("text_signal_overridden") == "xray", (
        "Expected text_signal_overridden=xray in sidecar; got %r"
        % cf.get("text_signal_overridden"))

    # Telemetry must mark the override
    correction_logs = [m for m in log_messages
                       if "[DIRECTIVE_CORRECTION]" in m
                       and "_MIXED" not in m]
    assert correction_logs, "Expected DIRECTIVE_CORRECTION log"
    assert "OVERRIDDEN" in correction_logs[0], (
        "Log must include OVERRIDDEN annotation; got %r"
        % correction_logs[0])
    print("  PASS")


# =====================================================================
# K19 — Mixed-input drift defers to text + emits [..._MIXED]
# =====================================================================

def test_K19_h18_mixed_input_defers_to_text():
    """K19: Pitfall 1 mitigation.  Long-running sessions can
    accumulate both .mtz and .ccp4 in original_files (per
    Session.set_project_info MERGE semantics).  File-based
    detection returns None for mixed input; we defer to text.

    A [DIRECTIVE_CORRECTION_MIXED] line MUST be logged so any
    future drift is visible in audit logs immediately."""
    print("Test: K19_h18_mixed_input_defers_to_text")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.autobuild_denmod"}}')

    log_messages = []
    def _capture_log(msg):
        log_messages.append(msg)

    # Text says cryo-em; files are mixed (.mtz + .ccp4)
    result = de.extract_directives(
        _PROCESSED_CRYOEM, provider="google",
        log_func=_capture_log,
        raw_advice="density modify and stop",
        original_files=_MIXED_FILES,
    )
    sc = result.get("stop_conditions", {})
    # Correction still fires via text fallback
    assert sc.get("after_program") == "phenix.resolve_cryo_em"
    cf = sc.get("_corrected_from") or {}
    assert cf.get("source") == "text", (
        "Expected source=text when files are mixed; got %r"
        % cf.get("source"))

    # Mixed-input audit marker
    mixed_logs = [m for m in log_messages
                  if "[DIRECTIVE_CORRECTION_MIXED]" in m]
    assert len(mixed_logs) == 1, (
        "Expected exactly one [DIRECTIVE_CORRECTION_MIXED] line; "
        "got %d: %r" % (len(mixed_logs), mixed_logs))
    assert ".mtz" in mixed_logs[0]
    assert ".ccp4" in mixed_logs[0]
    print("  PASS")


# =====================================================================
# K20 — plan_generator uses the shared helper (single SoT)
# =====================================================================

def test_K20_h18_plan_generator_uses_shared_helper():
    """K20: plan_generator._build_context delegates experiment-type
    inference to the shared infer_experiment_type_from_files helper
    (Item 3 of the H18 plan).

    Ensures plan_generator and directive_extractor cannot drift
    apart on detection behavior — both call the same function.

    Behavioral pin: mixed input yields None (semantically tighter
    than pre-H18 last-write-wins, which silently picked the
    LAST-seen extension).  In practice plan inputs are never mixed,
    but this test pins the new contract."""
    print("Test: K20_h18_plan_generator_uses_shared_helper")
    try:
        from libtbx.langchain.agent.file_utils import (
            infer_experiment_type_from_files)
        from libtbx.langchain.agent.plan_generator import _build_context
    except ImportError:
        from agent.file_utils import infer_experiment_type_from_files
        from agent.plan_generator import _build_context

    # Cryo-EM input: both must report cryoem
    files = _AF7MJS_FILES
    ctx = _build_context(available_files=files)
    helper_type, _ = infer_experiment_type_from_files(files)
    assert ctx["experiment_type"] == "cryoem"
    assert helper_type == "cryoem"

    # X-ray input: both must report xray
    ctx_x = _build_context(available_files=_XRAY_FILES)
    helper_x, _ = infer_experiment_type_from_files(_XRAY_FILES)
    assert ctx_x["experiment_type"] == "xray"
    assert helper_x == "xray"

    # Mixed input: both must report None (new tighter semantics)
    helper_mix, ev = infer_experiment_type_from_files(_MIXED_FILES)
    assert helper_mix is None, (
        "Helper must return None for mixed input; got %r"
        % helper_mix)
    assert ev["is_mixed"] is True
    print("  PASS")


# =====================================================================
# K21 — v117.2 fallback also threads original_files (H18.2)
# =====================================================================
#
# Production failure path discovered AFTER H18 and H18.1 shipped:
# the AF_7mjs run with LLM extraction produced
# stop_conditions={"stop_after_requested": True} but NO after_program.
# This triggered the v117.2 fallback at directive_extractor.py:783
# which called _resolve_after_program WITHOUT original_files.  The
# resolver defaulted _exp="xray" via the text-only heuristic and
# mapped denmod → phenix.autobuild_denmod.  The downstream
# _apply_workflow_intent_fallback ran with original_files but saw
# preprocessed advice with "Stop Condition: None" — so its
# _is_stop_after_requested returned False, and the n==1+no-stop
# branch left the buggy after_program in place.
#
# H18.2 fix: the v117.2 callsite now passes original_files=
# original_files to _resolve_after_program, mirroring the H18 policy
# at the OTHER two _resolve_after_program callsites.

def test_K21_h18_2_v117_2_fallback_threads_files():
    """K21: when the LLM emits stop_after_requested=True but no
    after_program, the v117.2 fallback path must use files-first
    experiment-type detection (mirroring H18's policy at the other
    two _resolve_after_program callsites).  Production failure
    reproduction; was the AF_7mjs cycle-3 predict_and_build bug."""
    print("Test: K21_h18_2_v117_2_fallback_threads_files")
    try:
        from libtbx.langchain.agent.directive_extractor import (
            extract_directives,
        )
    except ImportError:
        from agent.directive_extractor import extract_directives

    # Reproduce Tom's AF_7mjs case exactly:
    # - The LLM (mocked here) returns stop_after_requested=True but
    #   omits after_program (the exact production behavior observed
    #   in the runtime tracer output).
    # - Preprocessed advice has "Stop Condition: None"
    # - Raw advice has "density modify and stop"
    # - File inventory is cryo-EM (.ccp4 only)

    PROCESSED_ADVICE = (
        "1. **Input Files Found**: 7mjs_23883_H_1.ccp4, "
        "7mjs_23883_H_2.ccp4, 7mjs_23883_H.fa\n\n"
        "2. **Experiment Type**: cryo-EM (inferred from "
        "half-maps)\n\n"
        "3. **Primary Goal**: Perform density modification.\n\n"
        "4. **Key Parameters**: None\n\n"
        "5. **Program Parameters**: None\n\n"
        "6. **Special Instructions**: None\n\n"
        "7. **Stop Condition**: None"
    )
    RAW_ADVICE = "User instructions:\ndensity modify and stop"
    ORIGINAL_FILES = _AF7MJS_FILES

    # Reproduce the exact LLM output Tom's production trace showed:
    # stop_after_requested=True with NO after_program.  Use the
    # standard _with_mock_llm pattern that returns a JSON string
    # (matching the production protocol).
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}'
    )
    try:
        result = extract_directives(
            user_advice=PROCESSED_ADVICE,
            raw_advice=RAW_ADVICE,
            original_files=ORIGINAL_FILES,
            provider="google",
        )
    finally:
        # Restore — other tests will install their own mocks
        pass

    sc = (result or {}).get("stop_conditions", {})
    after_prog = sc.get("after_program")
    assert after_prog == "phenix.resolve_cryo_em", (
        "K21 production reproduction FAILED.  v117.2 fallback didn't "
        "thread original_files, so denmod defaulted to xray → "
        "phenix.autobuild_denmod.\n"
        "Got: after_program=%r\n"
        "Expected: phenix.resolve_cryo_em\n"
        "Full stop_conditions: %s" % (after_prog, sc))
    assert sc.get("stop_after_requested") is True
    print("  PASS")


def test_K22_h18_2_v117_2_fallback_xray_files():
    """K22: mirror of K21 for X-ray inputs.  When LLM omits
    after_program and files are X-ray (.mtz), v117.2 fallback must
    produce phenix.autobuild_denmod (the correct X-ray choice)."""
    print("Test: K22_h18_2_v117_2_fallback_xray_files")
    try:
        from libtbx.langchain.agent.directive_extractor import (
            extract_directives,
        )
    except ImportError:
        from agent.directive_extractor import extract_directives

    PROCESSED_ADVICE = (
        "1. **Input Files Found**: data.mtz, model.pdb\n\n"
        "2. **Experiment Type**: X-ray\n\n"
        "3. **Primary Goal**: Perform density modification.\n\n"
        "7. **Stop Condition**: None"
    )
    RAW_ADVICE = "User instructions:\ndensity modify and stop"
    XRAY_FILES = ["data.mtz", "model.pdb"]

    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}'
    )
    result = extract_directives(
        user_advice=PROCESSED_ADVICE,
        raw_advice=RAW_ADVICE,
        original_files=XRAY_FILES,
        provider="google",
    )

    sc = (result or {}).get("stop_conditions", {})
    after_prog = sc.get("after_program")
    assert after_prog == "phenix.autobuild_denmod", (
        "K22 FAILED.  Got: after_program=%r (expected "
        "phenix.autobuild_denmod for X-ray inputs)" % after_prog)
    print("  PASS")


def test_K23_h18_2_v117_2_no_files_falls_back_to_text():
    """K23: backward compat — when original_files is None, the
    v117.2 fallback must fall back to text-only detection (pre-H18
    behavior).  Confirms files=None doesn't break the existing
    text path."""
    print("Test: K23_h18_2_v117_2_no_files_falls_back_to_text")
    try:
        from libtbx.langchain.agent.directive_extractor import (
            extract_directives,
        )
    except ImportError:
        from agent.directive_extractor import extract_directives

    # Raw advice with strong cryo-EM text signal — text-only
    # detection should produce phenix.resolve_cryo_em
    RAW_ADVICE = ("User instructions:\nperform cryo-em density "
                  "modification and stop")
    PROCESSED_ADVICE = (
        "1. **Input Files Found**: None\n\n"
        "2. **Experiment Type**: cryo-EM\n\n"
        "3. **Primary Goal**: Density modify.\n\n"
        "7. **Stop Condition**: None"
    )

    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}'
    )
    # original_files=None — must not crash, must use text path
    result = extract_directives(
        user_advice=PROCESSED_ADVICE,
        raw_advice=RAW_ADVICE,
        original_files=None,
        provider="google",
    )

    sc = (result or {}).get("stop_conditions", {})
    after_prog = sc.get("after_program")
    assert after_prog == "phenix.resolve_cryo_em", (
        "K23 FAILED.  Got: after_program=%r (expected "
        "phenix.resolve_cryo_em from text-only cryo-em signal)"
        % after_prog)
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
        # v119.H18 additions
        test_K14_h18_af7mjs_terse_advice_cryoem_files,
        test_K15_h18_mirror_terse_advice_xray_files,
        test_K16_h18_backward_compat_text_fallback,
        test_K17_h18_no_signal_declines,
        test_K18_h18_files_win_on_conflict,
        test_K19_h18_mixed_input_defers_to_text,
        test_K20_h18_plan_generator_uses_shared_helper,
        test_K21_h18_2_v117_2_fallback_threads_files,
        test_K22_h18_2_v117_2_fallback_xray_files,
        test_K23_h18_2_v117_2_no_files_falls_back_to_text,
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
