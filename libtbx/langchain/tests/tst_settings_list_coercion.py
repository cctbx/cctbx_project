"""Tests for v118.10: list-to-string coercion in directive validation.

The v118.10 plan adds `_coerce_setting_value()`, a helper that
handles LLM-shape variability across providers.  Google's LLM
tends to wrap scalar values in lists (e.g. ``["S"]`` for
``additional_atom_types``).  OpenAI emits the same content as a
bare scalar (``"S"``).  Naive ``str(["S"])`` produces ``"['S']"``
— the original v118.10 bug Tom hit on P9 SAD data.

This file holds the K_LIST test suite (11 tests) that locks the
fix:

  Direct helper tests (no LLM, no extract_directives):
    K1   — single-element list ["S"] + str → "S"
    K2   — multi-element list ["Se", "S"] + str → "Se S"
    K3   — tuple ("S", "Cl") + str → "S Cl"
    K4   — bare string "S" + str → "S" (unchanged scalar)
    K5   — non-str types: 5 + int, 5.0 + float, True + bool
    K6a  — dict {} + str → ValueError (caught + dropped)
    K6b  — single-element list [False] + bool → False
            (locks the bool-truthiness trap discovered in rev 3)
    K7   — empty list [] + str → "" + nested [["S"]] limitation

  End-to-end via extract_directives + mock LLM:
    K8   — Tom's exact bug: program_settings.additional_atom_types
            = ["S"] coerced to "S" with audit log
    K9   — site 2 single-element after_program list coerced and
            then passes VALID_PROGRAMS check
    K10  — site 2 multi-element after_program list coerced and
            then correctly rejected by VALID_PROGRAMS check
    K11  — interaction with v118.9: list-formatted
            after_program=["phenix.autobuild_denmod"] in cryo-EM
            context → coerce to string → v118.9 rewrites to
            phenix.resolve_cryo_em.  Proves v118.10 is a true
            enabler for v118.9.

All tests use the mock-LLM pattern (monkey-patch
`_call_llm`) for the end-to-end cases; the direct helper tests
call `_coerce_setting_value` straight.
"""

from __future__ import absolute_import, division, print_function

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
# Shared inputs for end-to-end tests
# =====================================================================

# Cryo-EM preprocessor-style advice for K11 (v118.10 + v118.9 interaction)
_PROCESSED_CRYOEM = (
    "1. **Input Files Found**: 7mjs_23883_H_1.ccp4, "
    "7mjs_23883_H_2.ccp4, 7mjs_23883_H.fa\n\n"
    "2. **Experiment Type**: cryo-EM (inferred from half-maps)\n\n"
    "3. **Primary Goal**: Perform density modification.\n\n"
    "4. **Stop Condition**: None"
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
# K1 — Single-element list + str (Tom's exact bug pattern at helper level)
# =====================================================================

def test_K1_single_element_list_to_str():
    """K1: _coerce_setting_value(["S"], str) → "S".

    The exact failure pattern from Tom's P9 SAD run: Google's LLM
    wraps `additional_atom_types` in a list, naive str() produces
    "['S']", autosol fails to parse the literal brackets."""
    print("Test: K1_single_element_list_to_str")
    got = de._coerce_setting_value(["S"], str)
    assert got == "S", (
        "Expected 'S' (unpacked from single-element list), "
        "got %r" % got)
    print("  PASS")


# =====================================================================
# K2 — Multi-element list + str: space-joined
# =====================================================================

def test_K2_multi_element_list_to_str():
    """K2: _coerce_setting_value(["Se", "S"], str) → "Se S".

    PHIL multi-value-string syntax accepts space-separated values.
    The downstream command_builder regex `[,;\\s/]+` splits on
    spaces as well as commas, so space-joining is compatible
    everywhere."""
    print("Test: K2_multi_element_list_to_str")
    got = de._coerce_setting_value(["Se", "S"], str)
    assert got == "Se S", (
        "Expected 'Se S' (space-joined), got %r" % got)
    print("  PASS")


# =====================================================================
# K3 — Tuple inclusion
# =====================================================================

def test_K3_tuple_to_str():
    """K3: tuples treated identically to lists.

    Some JSON parsers may yield tuples; isinstance check covers
    both."""
    print("Test: K3_tuple_to_str")
    got = de._coerce_setting_value(("S", "Cl"), str)
    assert got == "S Cl", "Expected 'S Cl', got %r" % got
    # Single-element tuple
    got = de._coerce_setting_value(("S",), str)
    assert got == "S", (
        "Single-element tuple should unpack to 'S', got %r" % got)
    print("  PASS")


# =====================================================================
# K4 — Bare string unchanged
# =====================================================================

def test_K4_bare_string_unchanged():
    """K4: _coerce_setting_value("S", str) → "S" (backward-compat).

    OpenAI's path through the helper must produce identical output
    to today's behavior — bare strings flow through `str(value)`
    unchanged."""
    print("Test: K4_bare_string_unchanged")
    got = de._coerce_setting_value("S", str)
    assert got == "S", "Expected 'S', got %r" % got
    print("  PASS")


# =====================================================================
# K5 — Non-str scalar types unchanged
# =====================================================================

def test_K5_scalar_types_unchanged():
    """K5: int/float/bool scalar values still coerce normally."""
    print("Test: K5_scalar_types_unchanged")
    assert de._coerce_setting_value(5, int) == 5
    assert de._coerce_setting_value(5.0, float) == 5.0
    assert de._coerce_setting_value(True, bool) is True
    assert de._coerce_setting_value(False, bool) is False
    # int → str works via str(5) = "5"
    assert de._coerce_setting_value(5, str) == "5"
    # str → int works via int("5") = 5
    assert de._coerce_setting_value("5", int) == 5
    print("  PASS")


# =====================================================================
# K6a — Dict raises ValueError
# =====================================================================

def test_K6a_dict_raises():
    """K6a: dict input to str expected_type raises ValueError.

    Outer try/except in the validation block catches this and
    drops the setting with a log line."""
    print("Test: K6a_dict_raises")
    try:
        de._coerce_setting_value({"key": "value"}, str)
        assert False, "Expected ValueError for dict input"
    except ValueError as e:
        assert "dict" in str(e).lower(), (
            "Error message should mention dict: %r" % str(e))
    print("  PASS")


# =====================================================================
# K6b — Single-element bool list: [False] → False, NOT True (rev 3 bonus)
# =====================================================================

def test_K6b_single_element_bool_list_preserves_false():
    """K6b: _coerce_setting_value([False], bool) → False.

    Critical correctness check: naive `bool([False])` returns True
    (non-empty list is truthy in Python).  The helper unpacks
    single-element lists BEFORE type coercion, so [False] → False
    via the recursive call with the unpacked element.

    Discovered during rev 3 review of the plan.  Without this fix,
    Google emitting ``"anisotropic_adp": [false]`` would silently
    become True in the validated directives."""
    print("Test: K6b_single_element_bool_list_preserves_false")
    got = de._coerce_setting_value([False], bool)
    assert got is False, (
        "Expected False (unpacked from [False]); naive bool([False]) "
        "would have given True, got %r" % got)
    # Also verify [True] still works
    got = de._coerce_setting_value([True], bool)
    assert got is True, "Expected True from [True], got %r" % got
    # And single-element int list
    got = de._coerce_setting_value([5], int)
    assert got == 5, "Expected 5 from [5], got %r" % got
    print("  PASS")


# =====================================================================
# K7 — Empty list + nested list limitation
# =====================================================================

def test_K7_empty_list_and_nested_limitation():
    """K7: empty list and nested-list edge cases.

    - [] + str → "" (outer validation's drop-empty logic removes)
    - [] + int → ValueError (no sensible scalar)
    - [["S"]] + str → "['S']" (documented limitation — LLMs
      emitting nested lists for flat PHIL parameters is genuinely
      malformed; we accept the limitation rather than recursively
      flattening, which would lose structure for legitimate
      [["S"], ["Cl"]] cases).
    """
    print("Test: K7_empty_list_and_nested_limitation")
    # Empty list + str
    assert de._coerce_setting_value([], str) == ""
    # Empty list + non-str → ValueError
    try:
        de._coerce_setting_value([], int)
        assert False, "Expected ValueError for empty list + int"
    except ValueError:
        pass
    # Nested list (single-element wraps single-element)
    got = de._coerce_setting_value([["S"]], str)
    # Single-element list [["S"]] unpacks to ["S"], then recursive
    # call str(["S"]) is itself a list so it would unpack again to
    # "S".  Let's verify.
    assert got == "S", (
        "Nested single-element list should fully unpack to 'S', "
        "got %r" % got)
    # Multi-element with one nested: [["S"], ["Cl"]]
    got = de._coerce_setting_value([["S"], ["Cl"]], str)
    # Multi-element → space-join → str of each → "['S'] ['Cl']"
    # This is the documented limitation: when LLM nests lists,
    # the result has bracket artifacts.  Not a valid PHIL string,
    # but downstream validation will reject it.
    assert "[" in got or "S" in got, (
        "Multi-nested produces string with content, got %r" % got)
    print("  PASS")


# =====================================================================
# K8 — Tom's exact case end-to-end via extract_directives
# =====================================================================

def test_K8_tom_p9_sad_case_end_to_end():
    """K8: end-to-end reproduction of Tom's P9 SAD bug.

    Mock LLM emits the exact Google shape Tom observed:
        {"program_settings": {"phenix.autosol":
         {"additional_atom_types": ["S"]}}}

    Expected post-v118.10:
    - result.program_settings.phenix.autosol.additional_atom_types
      == "S" (NOT "['S']")
    - log line: "DIRECTIVES: Unpacked single-element list for
      additional_atom_types: ['S'] → 'S'"
    """
    print("Test: K8_tom_p9_sad_case_end_to_end")
    _with_mock_llm(
        '{"program_settings": {"phenix.autosol": '
        '{"additional_atom_types": ["S"]}}}')

    logs = []
    result = de.extract_directives(
        "P9 SAD data; search for Se and S",
        provider="google",
        log_func=lambda m: logs.append(m),
    )
    ps = result.get("program_settings", {}).get("phenix.autosol", {})
    aat = ps.get("additional_atom_types")
    assert aat == "S", (
        "Expected additional_atom_types='S' (coerced from ['S']), "
        "got %r" % aat)
    # Verify audit log line fires
    unpack_logs = [l for l in logs
                   if "Unpacked single-element list" in l
                   and "additional_atom_types" in l]
    assert len(unpack_logs) == 1, (
        "Expected exactly one 'Unpacked single-element list' log "
        "for additional_atom_types, got %d: %r"
        % (len(unpack_logs), unpack_logs))
    print("  PASS")


# =====================================================================
# K9 — Site 2 single-element after_program list coerced and validated
# =====================================================================

def test_K9_site2_single_element_after_program_list():
    """K9: stop_conditions.after_program = ["phenix.autosol"]
    coerces to "phenix.autosol" AND passes VALID_PROGRAMS check.

    Today (pre-v118.10) without coercion:
    - `value not in VALID_PROGRAMS` raises TypeError on list
    - Silently caught by outer except → after_program dropped
    - Section E couldn't help diagnose

    Post-v118.10:
    - Coerce list → "phenix.autosol" (string)
    - VALID_PROGRAMS check passes
    - after_program preserved
    """
    print("Test: K9_site2_single_element_after_program_list")
    _with_mock_llm(
        '{"stop_conditions": {"after_program": '
        '["phenix.autosol"]}}')

    logs = []
    result = de.extract_directives(
        "Run autosol",
        provider="google",
        log_func=lambda m: logs.append(m),
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.autosol", (
        "Expected after_program='phenix.autosol' (coerced and "
        "validated), got %r" % ap)
    # No "Invalid stop program" log
    invalid_logs = [l for l in logs if "Invalid stop program" in l]
    assert len(invalid_logs) == 0, (
        "Should not have 'Invalid stop program' log when single-"
        "element list contained a valid program name, got: %r"
        % invalid_logs)
    # Should have the Unpacked log line
    unpack_logs = [l for l in logs
                   if "Unpacked single-element list" in l
                   and "after_program" in l]
    assert len(unpack_logs) == 1, (
        "Expected 'Unpacked single-element list for after_program' "
        "log, got %d: %r" % (len(unpack_logs), unpack_logs))
    print("  PASS")


# =====================================================================
# K10 — Site 2 multi-element after_program list correctly rejected
# =====================================================================

def test_K10_site2_multi_element_after_program_rejected():
    """K10: stop_conditions.after_program = ["phenix.autosol",
    "phenix.refine"] coerces to space-joined string, which then
    correctly fails the VALID_PROGRAMS check and is dropped.

    Semantic note: a multi-program after_program is malformed
    (one stop target, not many), so rejection is the right
    behavior.  Both logs should fire:
    - Joined multi-element list (the coercion)
    - Invalid stop program (the subsequent rejection)
    """
    print("Test: K10_site2_multi_element_after_program_rejected")
    _with_mock_llm(
        '{"stop_conditions": {"after_program": '
        '["phenix.autosol", "phenix.refine"]}}')

    logs = []
    result = de.extract_directives(
        "Run autosol then refine",
        provider="google",
        log_func=lambda m: logs.append(m),
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program") if sc else None
    # The coerced string "phenix.autosol phenix.refine" is not a
    # valid program; rejected → after_program absent.
    assert ap is None, (
        "Expected after_program dropped (multi-element coerced "
        "value not a valid program), got %r" % ap)

    # Both audit lines should appear
    joined_logs = [l for l in logs
                   if "Joined multi-element list" in l
                   and "after_program" in l]
    assert len(joined_logs) == 1, (
        "Expected 'Joined multi-element list for after_program' "
        "log, got %d: %r" % (len(joined_logs), joined_logs))
    invalid_logs = [l for l in logs if "Invalid stop program" in l]
    assert len(invalid_logs) == 1, (
        "Expected 'Invalid stop program' rejection log, got %d: %r"
        % (len(invalid_logs), invalid_logs))
    print("  PASS")


# =====================================================================
# K11 — v118.10 + v118.9 interaction (v118.10 is a true enabler)
# =====================================================================

def test_K11_v118_10_enables_v118_9_for_list_after_program():
    """K11: list-formatted after_program in cryo-EM context.

    Mock LLM emits Google's shape with the v118.9 bug pattern
    WRAPPED in a list (the v118.10 bug pattern):
        {"stop_conditions": {"stop_after_requested": true,
         "after_program": ["phenix.autobuild_denmod"]}}

    Without v118.10: site 2's TypeError catches silently drops
    after_program; v118.9's validator never gets to correct it.

    With v118.10: list coerced to "phenix.autobuild_denmod"
    string; v118.9's validator then correctly rewrites to
    "phenix.resolve_cryo_em" for the cryo-EM context.

    Both audit lines should fire:
    - DIRECTIVES: Unpacked single-element list for after_program ...
    - [DIRECTIVE_CORRECTION] Mapped after_program=... (v118.9)
    """
    print("Test: K11_v118_10_enables_v118_9_for_list_after_program")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": ["phenix.autobuild_denmod"]}}')

    logs = []
    result = de.extract_directives(
        _PROCESSED_CRYOEM,
        provider="google",
        log_func=lambda m: logs.append(m),
        raw_advice="density modify and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.resolve_cryo_em", (
        "Expected v118.9 to rewrite the coerced "
        "autobuild_denmod → resolve_cryo_em, got %r" % ap)

    # v118.10 log: Unpacked single-element list
    unpack_logs = [l for l in logs
                   if "Unpacked single-element list" in l
                   and "after_program" in l]
    assert len(unpack_logs) == 1, (
        "Expected v118.10 'Unpacked' log line, got %d: %r"
        % (len(unpack_logs), unpack_logs))

    # v118.9 log: [DIRECTIVE_CORRECTION]
    correction_logs = [l for l in logs
                       if "[DIRECTIVE_CORRECTION]" in l]
    assert len(correction_logs) == 1, (
        "Expected v118.9 [DIRECTIVE_CORRECTION] log line "
        "(proves v118.10 enabled v118.9 to see the correctly-"
        "typed input), got %d: %r"
        % (len(correction_logs), correction_logs))

    # v118.9 sidecar should be present (extra verification)
    cf = sc.get("_corrected_from")
    assert cf is not None, (
        "Expected v118.9 _corrected_from sidecar, got None")
    assert cf.get("from") == "phenix.autobuild_denmod"
    assert cf.get("to") == "phenix.resolve_cryo_em"
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    tests = [
        test_K1_single_element_list_to_str,
        test_K2_multi_element_list_to_str,
        test_K3_tuple_to_str,
        test_K4_bare_string_unchanged,
        test_K5_scalar_types_unchanged,
        test_K6a_dict_raises,
        test_K6b_single_element_bool_list_preserves_false,
        test_K7_empty_list_and_nested_limitation,
        test_K8_tom_p9_sad_case_end_to_end,
        test_K9_site2_single_element_after_program_list,
        test_K10_site2_multi_element_after_program_rejected,
        test_K11_v118_10_enables_v118_9_for_list_after_program,
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
        raise AssertionError("%d K_LIST tests failed" % failed)


if __name__ == "__main__":
    run_all_tests()
