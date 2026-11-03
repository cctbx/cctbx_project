"""Tests for `_format_directives_for_prompt` after_program rendering
under the v117 stop_refactor architecture (stop_after_requested-gated).

Under this architecture, `after_program` carries two distinct kinds
of intent and the renderer emits DIFFERENT prompt content for each:

  (a) USER-EXPLICIT STOP — set by the directive extractor when the
      raw advice contains "stop after X" / "X and stop" / etc.
      Identified by `stop_conditions.stop_after_requested = True`.
      Renders the Phase 6a 4-line "Stop target:" block.

  (b) PLAN-PROGRESSION HINT — emitted per-stage by plan_to_directives
      to ensure the stage's program runs.  Identified by absence of
      the `stop_after_requested` flag.  Renders NO after_program
      content (the workflow_engine min-run guarantee prioritizes the
      target in valid_programs; telling the LLM "stop after X" would
      be misleading).

The "**Stop Conditions:**" header itself is emitted only when there
is non-empty content under it (the body-then-header pattern in
_format_directives_for_prompt).

Test groups:
  K1-K3 — adjacent stop_conditions fields (after_cycle, max_refine,
          r_free_target, etc.) still render correctly
  N1-N5 — stop_after_requested gating:
    N1: True → Phase 6a "Stop target:" block emitted
    N2: True → schema-correct in-list / not-in-list lines emitted
    N3: False/absent → no after_program content emitted
    N4: True + plan-progression hint mixed → only stop_after_requested
        controls the rendering
    N5: anti-regression — the OLD "CRITICAL: You MUST run X" wording
        must NOT appear in any active emission path
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


_HERE = os.path.dirname(os.path.abspath(__file__))
_KNOWLEDGE_DIR = os.path.normpath(os.path.join(_HERE, "..", "knowledge"))
if _KNOWLEDGE_DIR not in sys.path:
    sys.path.insert(0, _KNOWLEDGE_DIR)


def _install_stubs():
    """Install stub modules for libtbx imports."""

    class _StubProgramRegistry:
        def __init__(self):
            pass

    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")
    pr = _ensure("libtbx.langchain.agent.program_registry")
    pr.ProgramRegistry = _StubProgramRegistry


try:
    from libtbx.langchain.knowledge.prompts_hybrid import (
        _format_directives_for_prompt)
except ImportError:
    _install_stubs()
    try:
        from libtbx.langchain.knowledge.prompts_hybrid import (
            _format_directives_for_prompt)
    except ImportError:
        from prompts_hybrid import _format_directives_for_prompt


# =====================================================================
# K1-K3: ADJACENT FIELDS STILL RENDER
# =====================================================================

def test_K1_after_cycle_still_works():
    """after_cycle is independent of stop_after_requested gating."""
    print("Test: K1_after_cycle_still_works")
    directives = {"stop_conditions": {"after_cycle": 5}}
    output = _format_directives_for_prompt(directives)
    assert "**Stop Conditions:**" in output, (
        "Stop Conditions header missing: %r" % output)
    assert "Stop after cycle 5" in output, (
        "after_cycle line missing: %r" % output)
    print("  PASS")


def test_K2_max_refine_cycles_and_targets_render():
    """max_refine_cycles and target metrics still render."""
    print("Test: K2_max_refine_cycles_and_targets_render")
    directives = {"stop_conditions": {
        "max_refine_cycles": 8,
        "r_free_target": 0.25,
        "map_cc_target": 0.70,
        "skip_validation": True,
    }}
    output = _format_directives_for_prompt(directives)
    assert "**Stop Conditions:**" in output, (
        "Header missing: %r" % output)
    assert "Maximum 8 refinement cycles" in output, (
        "max_refine_cycles missing")
    assert "Target R-free: 0.250" in output, "r_free_target missing"
    assert "Target map CC: 0.70" in output, "map_cc_target missing"
    assert "Validation can be skipped before stopping" in output, (
        "skip_validation missing")
    print("  PASS")


def test_K3_no_stop_conditions_renders_no_header():
    """When directives has no stop_conditions, no header is emitted."""
    print("Test: K3_no_stop_conditions_renders_no_header")
    output = _format_directives_for_prompt({})
    assert "**Stop Conditions:**" not in output, (
        "Header should not appear for empty directives")
    output_none = _format_directives_for_prompt(None)
    assert "**Stop Conditions:**" not in output_none
    print("  PASS")


# =====================================================================
# N1-N5: stop_after_requested GATING
# =====================================================================

def test_N1_stop_after_requested_true_emits_phase_6a_block():
    """When stop_after_requested=True, the Phase 6a 4-line block
    is emitted with the program name substituted."""
    print("Test: N1_stop_after_requested_true_emits_phase_6a_block")
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
        "stop_after_requested": True,
    }}
    output = _format_directives_for_prompt(directives)
    assert "**Stop Conditions:**" in output, "Header missing"
    # Line 1 — simple stop confirmation
    assert "Stop after phenix.refine completes" in output, (
        "Stop-after line missing: %r" % output)
    # Line 2 — Stop target framing
    assert ("**Stop target: phenix.refine.**" in output
            or "Stop target: phenix.refine." in output), (
        "Stop target line missing: %r" % output)
    # Lines 3-4 — in-list / not-in-list guidance
    assert "in VALID PROGRAMS this cycle, choose it" in output, (
        "in-list guidance missing: %r" % output)
    assert "NOT in VALID PROGRAMS" in output, (
        "not-in-list guidance missing: %r" % output)
    assert ("Never pick a program outside VALID PROGRAMS"
            in output), (
        "outside-list prohibition missing: %r" % output)
    print("  PASS")


def test_N2_stop_after_requested_substitutes_program_name():
    """Different programs render different lines correctly."""
    print("Test: N2_stop_after_requested_substitutes_program_name")
    for prog in ["phenix.refine", "phenix.predict_and_build",
                 "phenix.resolve_cryo_em", "phenix.ligandfit"]:
        directives = {"stop_conditions": {
            "after_program": prog,
            "stop_after_requested": True,
        }}
        output = _format_directives_for_prompt(directives)
        assert "Stop after %s completes" % prog in output, (
            "Stop-after line not substituted for %s" % prog)
        assert "Stop target: %s." % prog in output, (
            "Stop target line not substituted for %s" % prog)
    print("  PASS")


def test_N3_plan_progression_hint_emits_no_after_program_content():
    """When stop_after_requested is absent or False, after_program
    is a plan-progression hint and the renderer emits NOTHING for
    it.  The "**Stop Conditions:**" header itself is suppressed if
    there is no other content to render."""
    print("Test: N3_plan_progression_hint_emits_no_after_program_content")
    # Case 1: after_program alone (no stop_after_requested) — header SUPPRESSED
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
    }}
    output = _format_directives_for_prompt(directives)
    assert "Stop after phenix.refine completes" not in output, (
        "after_program content should not render without "
        "stop_after_requested=True")
    assert "Stop target: phenix.refine" not in output, (
        "Stop target line should not render without "
        "stop_after_requested=True")
    # Header should also be suppressed because there's no other body
    assert "**Stop Conditions:**" not in output, (
        "Header should be suppressed when after_program is the only "
        "stop_condition and stop_after_requested is False")

    # Case 2: explicit False
    directives2 = {"stop_conditions": {
        "after_program": "phenix.refine",
        "stop_after_requested": False,
    }}
    output2 = _format_directives_for_prompt(directives2)
    assert "Stop after phenix.refine completes" not in output2
    assert "Stop target: phenix.refine" not in output2
    print("  PASS")


def test_N4_mixed_after_program_and_adjacent_fields_render_correctly():
    """When after_program is a plan hint (no stop_after_requested) but
    other stop_conditions are present, those OTHER fields render under
    the header but after_program produces no content."""
    print("Test: N4_mixed_after_program_and_adjacent_fields_render_correctly")
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
        # NO stop_after_requested
        "after_cycle": 3,
        "r_free_target": 0.22,
    }}
    output = _format_directives_for_prompt(directives)
    # Header appears because there is non-empty content
    assert "**Stop Conditions:**" in output
    assert "Stop after cycle 3" in output
    assert "Target R-free: 0.220" in output
    # But the after_program "Stop target:" lines do NOT appear
    assert "Stop target: phenix.refine" not in output, (
        "after_program plan hint should not emit Stop target line")
    assert "Stop after phenix.refine completes" not in output, (
        "after_program plan hint should not emit stop-after line")
    print("  PASS")


def test_N5_old_critical_wording_does_not_appear_in_output():
    """Anti-regression: the OLD 'CRITICAL: You MUST run X' wording
    must NOT appear in any rendering, even when stop_after_requested
    is True.  The new Phase 6a wording replaces it.

    If a future edit accidentally re-introduces the old wording, this
    test catches it."""
    print("Test: N5_old_critical_wording_does_not_appear_in_output")
    for stop_after in [True, False]:
        directives = {"stop_conditions": {
            "after_program": "phenix.refine",
            "stop_after_requested": stop_after,
        }}
        output = _format_directives_for_prompt(directives)
        # These are the OLD lines that v117 Phase 6a replaced
        assert "CRITICAL: You MUST run" not in output, (
            "Old 'CRITICAL: You MUST run' wording leaked back into "
            "output (stop_after_requested=%s)" % stop_after)
        assert "Do NOT keep running refinement cycles" not in output, (
            "Old 'Do NOT keep running' wording leaked back into "
            "output (stop_after_requested=%s)" % stop_after)
    print("  PASS")


# =====================================================================
# Sanity edge cases
# =====================================================================

def test_empty_directives():
    """{} → no Stop Conditions section."""
    print("Test: empty_directives")
    output = _format_directives_for_prompt({})
    assert "Stop Conditions" not in output
    print("  PASS")


def test_none_directives():
    """None → renderer handles gracefully."""
    print("Test: none_directives")
    output = _format_directives_for_prompt(None)
    assert "Stop Conditions" not in output
    print("  PASS")


def test_no_unresolved_placeholders():
    """No {placeholder} survives substitution."""
    print("Test: no_unresolved_placeholders")
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
        "stop_after_requested": True,
    }}
    output = _format_directives_for_prompt(directives)
    # Naive check: no obvious leftover format placeholders
    for pattern in ["{after_prog}", "{program}", "%(prog)s"]:
        assert pattern not in output, (
            "Unresolved placeholder %r in output: %r" % (pattern, output))
    print("  PASS")


def run_all_tests():
    tests = [
        test_K1_after_cycle_still_works,
        test_K2_max_refine_cycles_and_targets_render,
        test_K3_no_stop_conditions_renders_no_header,
        test_N1_stop_after_requested_true_emits_phase_6a_block,
        test_N2_stop_after_requested_substitutes_program_name,
        test_N3_plan_progression_hint_emits_no_after_program_content,
        test_N4_mixed_after_program_and_adjacent_fields_render_correctly,
        test_N5_old_critical_wording_does_not_appear_in_output,
        test_empty_directives,
        test_none_directives,
        test_no_unresolved_placeholders,
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
