"""Tests for the v117.1 fix: grounding guardrail must skip
after_program drop when stop_after_requested=True is set.

The interaction the fix addresses:
  - v116.19a grounding guardrail drops after_program when the program
    name isn't grounded in user_advice (Failure 1: pure fabrication;
    Failure 2: preprocessed + no stop intent + no imperative).
  - v117 Step 1 added stop_after_requested — the LLM's explicit
    user-stop assertion based on raw_advice + AUTHORITY paragraph.
  - When both signals fire (LLM extracts after_program AND
    stop_after_requested=True), the flag IS the grounding; the
    guardrail's heuristic check is redundant and harmful.

This file documents the fix and the boundaries:
  - K1: LLM extracts both signals -> after_program preserved
  - K2: LLM extracts only after_program (no flag) -> still dropped
  - K3: pure fabrication scenario -> still dropped
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

try:
    from libtbx.langchain.agent.directive_extractor import (
        _validate_after_program_grounded)
except ImportError:
    from agent.directive_extractor import _validate_after_program_grounded


def _make_logger():
    msgs = []
    def log(m):
        msgs.append(m)
    return log, msgs


def test_K1_llm_sets_both_signals_keeps_after_program():
    """v117.1 fix: when LLM sets after_program AND stop_after_requested,
    after_program is preserved even if not literally in user_advice."""
    print("Test: K1_llm_sets_both_signals_keeps_after_program")
    log, msgs = _make_logger()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.resolve_cryo_em",
            "stop_after_requested": True,
            "skip_validation": True,
        }
    }
    processed = ("1. Input Files Found: model.pdb, data.mtz\n"
                 "Primary Goal: Perform density modification, "
                 "then build and refine.\n"
                 "Stop Condition: None")
    out = _validate_after_program_grounded(directives, processed, log)
    ap = (out.get("stop_conditions") or {}).get("after_program")
    assert ap == "phenix.resolve_cryo_em", (
        "after_program should be preserved when stop_after_requested=True; "
        "got %r" % ap)
    keep_msgs = [m for m in msgs if "Keeping" in m]
    assert keep_msgs, "expected a 'Keeping' log message"
    print("  PASS")


def test_K2_llm_sets_only_after_program_drops_it():
    """Boundary: LLM sets after_program but NOT stop_after_requested.
    The grounding check fires as before -> drop."""
    print("Test: K2_llm_sets_only_after_program_drops_it")
    log, _ = _make_logger()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.real_space_refine",
            # NO stop_after_requested
        }
    }
    af7mjs_advice = ("1. Input Files Found: 7mjs.fa, 7mjs_1.ccp4\n"
                     "3. Primary Goal: Run PredictAndBuild and refine.\n"
                     "7. Stop Condition: None")
    out = _validate_after_program_grounded(directives, af7mjs_advice, log)
    ap = (out.get("stop_conditions") or {}).get("after_program")
    assert ap is None, (
        "after_program should be dropped without stop_after_requested; "
        "got %r" % ap)
    print("  PASS")


def test_K3_pure_fabrication_no_flag_drops():
    """Boundary: pure fabrication scenario (program name absent from
    advice, no flag) -> still dropped."""
    print("Test: K3_pure_fabrication_no_flag_drops")
    log, _ = _make_logger()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.something_completely_made_up",
        }
    }
    out = _validate_after_program_grounded(
        directives, "Run molecular replacement and refine.", log)
    ap = (out.get("stop_conditions") or {}).get("after_program")
    assert ap is None, "pure fabrication should be dropped; got %r" % ap
    print("  PASS")


def test_K4_stop_after_requested_false_drops():
    """Boundary: stop_after_requested=False (not absent) -> grounding
    check still fires, after_program dropped if not grounded."""
    print("Test: K4_stop_after_requested_false_drops")
    log, _ = _make_logger()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.something_made_up",
            "stop_after_requested": False,
        }
    }
    out = _validate_after_program_grounded(
        directives, "Run refinement.", log)
    ap = (out.get("stop_conditions") or {}).get("after_program")
    assert ap is None, (
        "after_program should be dropped when stop_after_requested=False; "
        "got %r" % ap)
    print("  PASS")


def test_K5_grounded_program_no_flag_kept():
    """Sanity: program name literally in advice + no flag -> kept
    (existing v116.19a behavior, unchanged by fix)."""
    print("Test: K5_grounded_program_no_flag_kept")
    log, _ = _make_logger()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine",
        }
    }
    out = _validate_after_program_grounded(
        directives, "Run phenix.refine on the model.", log)
    ap = (out.get("stop_conditions") or {}).get("after_program")
    assert ap == "phenix.refine", (
        "after_program should be kept (grounded); got %r" % ap)
    print("  PASS")


def run_all_tests():
    tests = [
        test_K1_llm_sets_both_signals_keeps_after_program,
        test_K2_llm_sets_only_after_program_drops_it,
        test_K3_pure_fabrication_no_flag_drops,
        test_K4_stop_after_requested_false_drops,
        test_K5_grounded_program_no_flag_kept,
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
