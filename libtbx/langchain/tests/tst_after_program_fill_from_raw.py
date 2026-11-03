"""Tests for v117.2 fix: when LLM emits stop_after_requested=True but
omits after_program, fill it in from raw advice via the regex resolver.

The interaction this fix addresses:

  v117 Step 1 introduced raw-advice dual-input extraction.  The LLM
  uses the AUTHORITY paragraph to honor raw intent over processed
  intent, but occasionally emits the stop flag without specifying
  after_program — observed on openai with raw='refine and stop' +
  processed advice containing 'Stop Condition: None'.

  Without after_program, workflow_engine._apply_directives' gated
  wipe code is unreachable (gated on `if after_program:`).  The
  user's stop intent is signaled but cannot be acted on.

  v117.2 closes this gap: when stop_after_requested=True and
  after_program is missing, call _resolve_after_program() on the
  raw advice to fill it in via _ACTION_TABLE.

Boundary cases (K1-K5):
  K1: flag True + after_program missing + raw has stop intent
      -> after_program filled in from raw
  K2: flag True + after_program ALREADY set -> no change (LLM choice
      wins, fix doesn't fire)
  K3: flag True + after_program missing + raw has NO stop intent
      -> no change (fix doesn't fire, can't trust the flag)
  K4: flag False/absent + after_program missing -> no change
  K5: raw_advice None (single-input path) -> falls back to user_advice
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


def _with_mock_llm(response_str):
    """Decorator: replace _call_llm with one that returns a fixed string."""
    def _fake(prompt, provider, model, log):
        return response_str
    de._call_llm = _fake


_PROCESSED_C1_REFINE = (
    "1. Input Files Found: model.pdb, data.mtz\n\n"
    "2. Experiment Type: X-ray crystallography\n\n"
    "3. Primary Goal: Refine the model and then validate using "
    "MolProbity to assess geometry.\n\n"
    "4. Key Parameters:\n"
    "   - Resolution limit: 2.0\n\n"
    "5. Program Parameters: None\n\n"
    "6. Special Instructions: None\n\n"
    "7. Stop Condition: None"
)


def test_K1_flag_true_no_after_program_filled_from_raw():
    """K1: LLM emits flag only; raw='refine and stop'.
    v117.2 fills in after_program=phenix.refine."""
    print("Test: K1_flag_true_no_after_program_filled_from_raw")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}')
    result = de.extract_directives(
        _PROCESSED_C1_REFINE, provider="openai",
        log_func=lambda x: None,
        raw_advice="refine and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    sar = sc.get("stop_after_requested")
    assert ap == "phenix.refine", (
        "Expected after_program=phenix.refine, got %r" % ap)
    assert sar is True, "Expected stop_after_requested=True, got %r" % sar
    print("  PASS")


def test_K2_llm_set_after_program_preserved():
    """K2: LLM emits both signals; fix is no-op (after_program is
    already set, the new code path's guard skips it)."""
    print("Test: K2_llm_set_after_program_preserved")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true, '
        '"after_program": "phenix.refine"}}')
    result = de.extract_directives(
        _PROCESSED_C1_REFINE, provider="openai",
        log_func=lambda x: None,
        raw_advice="refine and stop",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.refine", (
        "Expected after_program=phenix.refine, got %r" % ap)
    print("  PASS")


def test_K3_flag_true_no_stop_intent_in_raw_no_fill():
    """K3: LLM emits flag but raw='just refine the model' has NO stop
    intent.  Guard prevents resolver invocation; after_program stays
    missing (the LLM's flag is possibly hallucinated; we don't amplify)."""
    print("Test: K3_flag_true_no_stop_intent_in_raw_no_fill")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}')
    result = de.extract_directives(
        "Just refine the model.", provider="openai",
        log_func=lambda x: None,
        raw_advice="Just refine the model.",
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    # The fix did not fill in (no stop intent in raw).  The flag's
    # fate after that is the caller's problem — but the new code path
    # didn't make anything worse.
    assert ap is None, (
        "Expected after_program=None (fix did not fire), got %r" % ap)
    print("  PASS")


def test_K4_flag_absent_no_change():
    """K4: LLM does NOT set the flag.  Fix is no-op."""
    print("Test: K4_flag_absent_no_change")
    _with_mock_llm(
        '{"program_settings": {"default": {"resolution": 2.0}}}')
    result = de.extract_directives(
        _PROCESSED_C1_REFINE, provider="openai",
        log_func=lambda x: None,
        raw_advice="refine and stop",
    )
    sc = result.get("stop_conditions", {})
    # No stop_conditions block at all — fix didn't fire.
    assert sc.get("after_program") is None, (
        "Expected no after_program; got %r" % sc.get("after_program"))
    print("  PASS")


def test_K5_no_raw_advice_falls_back_to_user_advice():
    """K5: caller didn't pass raw_advice (single-input path).
    Fix falls back to user_advice.  If user_advice has stop intent,
    the resolver runs on it.

    Example: user_advice = 'refine and stop' (a caller using the
    single-input path with bare imperative)."""
    print("Test: K5_no_raw_advice_falls_back_to_user_advice")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}')
    result = de.extract_directives(
        "refine and stop", provider="openai",
        log_func=lambda x: None,
        # raw_advice intentionally omitted (None)
    )
    sc = result.get("stop_conditions", {})
    ap = sc.get("after_program")
    assert ap == "phenix.refine", (
        "Expected fallback to user_advice; after_program=phenix.refine, "
        "got %r" % ap)
    print("  PASS")


def test_K6_xtriage_raw_fills_in_xtriage_program():
    """K6: raw='run xtriage and stop' produces after_program=phenix.xtriage.
    Companion to K1 confirming the resolver handles other action verbs."""
    print("Test: K6_xtriage_raw_fills_in_xtriage_program")
    _with_mock_llm(
        '{"stop_conditions": {"stop_after_requested": true}}')
    result = de.extract_directives(
        "Some preprocessed prose without explicit stop.",
        provider="openai", log_func=lambda x: None,
        raw_advice="run xtriage and stop",
    )
    sc = result.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.xtriage", (
        "Expected phenix.xtriage, got %r" % sc.get("after_program"))
    print("  PASS")


def run_all_tests():
    tests = [
        test_K1_flag_true_no_after_program_filled_from_raw,
        test_K2_llm_set_after_program_preserved,
        test_K3_flag_true_no_stop_intent_in_raw_no_fill,
        test_K4_flag_absent_no_change,
        test_K5_no_raw_advice_falls_back_to_user_advice,
        test_K6_xtriage_raw_fills_in_xtriage_program,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            print("  FAIL: %s" % e)
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
