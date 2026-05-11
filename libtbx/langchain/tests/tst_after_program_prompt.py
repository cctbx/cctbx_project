"""Tests for v116.10 Phase 6a — LLM prompt alignment for after_program.

The bug: in `_format_directives_for_prompt`, the after_program block
included these lines:

    "- **CRITICAL: You MUST run X before stopping. If it's in VALID
      PROGRAMS, choose it NOW.**"
    "- Do NOT keep running refinement cycles - run X instead!"

The weak "If it's in VALID PROGRAMS" qualifier was undermined by
the strong "Do NOT...run X instead!" directive.  The LLM followed
the strong directive and picked after_program even when it was not
in VALID PROGRAMS — triggering the after_program_not_available
STOP defense at graph_nodes.py:2208.

The fix: target-not-now framing that explicitly tells the LLM:
- after_program is a STOP TARGET, not a now-directive
- Pick from VALID PROGRAMS only
- When the target is not available yet, pick a prerequisite

These tests verify:
(1) The new wording is present
(2) The old (bad) wording is gone
(3) The function still handles the absent-after_program case
(4) Substitution into after_prog placeholders happens correctly
(5) Functional integrity (empty/None inputs don't crash)

These are STATIC regression tests against re-introduction of the
bad text.  They do NOT verify that the LLM actually behaves
correctly with the new prompt — that requires end-to-end
verification against the tutorial corpus (run_29_openai,
run_29_ollama) and is documented in PHASE_6A_SUMMARY.md.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# --- Path setup + libtbx stubs ---------------------------------------------
# prompts_hybrid.py imports from libtbx.langchain.agent.program_registry.
# We stub it so the module can be imported standalone.

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
# SECTION A: The Phase 6a fix — target-not-now framing
# =====================================================================

def test_new_wording_present_in_after_program_block():
    """The new target-not-now wording must appear in the prompt."""
    print("Test: new_wording_present_in_after_program_block")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)

    expected_phrases = [
        "Stop target: phenix.predict_and_build.",
        "If phenix.predict_and_build is in VALID PROGRAMS",
        "If phenix.predict_and_build is NOT in VALID PROGRAMS",
        "The workflow will reach phenix.predict_and_build",
        "Never pick a program outside VALID PROGRAMS",
    ]
    for phrase in expected_phrases:
        assert phrase in output, (
            "Phase 6a wording missing: '%s'\nFull output:\n%s"
            % (phrase, output))
    print("  PASS")


def test_old_wording_removed():
    """The pre-Phase-6a bad lines must be gone."""
    print("Test: old_wording_removed")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)

    forbidden_phrases = [
        # The strong-MUST line that caused LLM to pick out-of-list:
        "CRITICAL: You MUST run",
        # The doubly-bad imperative:
        "Do NOT keep running refinement cycles",
        # The "run X instead!" exhortation:
        "run %s instead!" % "phenix.predict_and_build",
    ]
    for phrase in forbidden_phrases:
        assert phrase not in output, (
            "Old Phase-6a wording leaked back in: '%s'\nFull output:\n%s"
            % (phrase, output))
    print("  PASS")


# =====================================================================
# SECTION B: Substitution correctness
# =====================================================================

def test_after_prog_substitution_complete():
    """All `%s` placeholders in the new wording resolve to after_prog.

    Rather than count occurrences (brittle if the wording changes),
    we verify each templated phrase rendered with the program name.
    """
    print("Test: after_prog_substitution_complete")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)

    # Each templated phrase must appear with the substituted name.
    # If any %s slot fails to render, one of these phrases breaks.
    templated_phrases = [
        "Stop after phenix.predict_and_build completes",
        "Stop target: phenix.predict_and_build.",
        "If phenix.predict_and_build is in VALID PROGRAMS",
        "If phenix.predict_and_build is NOT in VALID PROGRAMS",
        "reach phenix.predict_and_build when its inputs",
    ]
    for phrase in templated_phrases:
        assert phrase in output, (
            "Templated phrase did not render: '%s'\nFull output:\n%s"
            % (phrase, output))

    # And no unresolved placeholders leaked through.
    assert "%s" not in output, (
        "Unresolved %%s placeholder in output:\n%s" % output)
    assert "%d" not in output, (
        "Unresolved %%d placeholder in output:\n%s" % output)
    print("  PASS")


def test_substitution_with_different_program_names():
    """Verify substitution works for various program names, not just predict_and_build."""
    print("Test: substitution_with_different_program_names")
    for prog in ["phenix.refine", "phenix.molprobity",
                 "phenix.real_space_refine", "phenix.ligandfit"]:
        directives = {"stop_conditions": {"after_program": prog}}
        output = _format_directives_for_prompt(directives)
        # Verify the core templated phrases render for this name
        for phrase in [
            "Stop after %s completes" % prog,
            "Stop target: %s." % prog,
            "If %s is in VALID PROGRAMS" % prog,
            "If %s is NOT in VALID PROGRAMS" % prog,
        ]:
            assert phrase in output, (
                "Phrase '%s' did not render for %s\nOutput:\n%s"
                % (phrase, prog, output))
        # No unresolved placeholders for this program either
        assert "%s" not in output, (
            "Unresolved %%s for %s in output:\n%s" % (prog, output))
    print("  PASS")


# =====================================================================
# SECTION C: Adjacent fields unaffected
# =====================================================================

def test_after_cycle_still_works():
    """`after_cycle` rendering must be unchanged by Phase 6a."""
    print("Test: after_cycle_still_works")
    directives = {"stop_conditions": {"after_cycle": 5}}
    output = _format_directives_for_prompt(directives)
    assert "Stop after cycle 5" in output, (
        "after_cycle rendering broken\nOutput:\n%s" % output)
    # Should NOT have any after_program content
    assert "Stop target" not in output, (
        "after_cycle-only directives should not produce Stop target line\n"
        "Output:\n%s" % output)
    assert "VALID PROGRAMS" not in output, (
        "after_cycle-only directives should not mention VALID PROGRAMS\n"
        "Output:\n%s" % output)
    print("  PASS")


def test_both_after_cycle_and_after_program():
    """Both stop conditions can coexist."""
    print("Test: both_after_cycle_and_after_program")
    directives = {
        "stop_conditions": {
            "after_cycle": 10,
            "after_program": "phenix.refine",
        }
    }
    output = _format_directives_for_prompt(directives)
    assert "Stop after cycle 10" in output, (
        "after_cycle missing when both present\nOutput:\n%s" % output)
    assert "Stop target: phenix.refine." in output, (
        "Stop target line missing when both present\nOutput:\n%s" % output)
    assert "If phenix.refine is in VALID PROGRAMS" in output, (
        "in-list guidance missing when both present\nOutput:\n%s" % output)
    print("  PASS")


def test_other_stop_condition_fields_preserved():
    """r_free_target, map_cc_target, etc. still render."""
    print("Test: other_stop_condition_fields_preserved")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.refine",
            "r_free_target": 0.22,
            "map_cc_target": 0.75,
            "max_refine_cycles": 5,
            "skip_validation": True,
        }
    }
    output = _format_directives_for_prompt(directives)
    expected = [
        "Target R-free: 0.220",
        "Target map CC: 0.75",
        "Maximum 5 refinement cycles",
        "Validation can be skipped before stopping",
        "Stop target: phenix.refine.",  # Phase 6a content also present
    ]
    for phrase in expected:
        assert phrase in output, (
            "Expected phrase missing: '%s'\nOutput:\n%s" % (phrase, output))
    print("  PASS")


# =====================================================================
# SECTION D: Edge cases
# =====================================================================

def test_no_after_program_field():
    """No after_program → no Phase-6a block at all."""
    print("Test: no_after_program_field")
    directives = {"stop_conditions": {"after_cycle": 3}}
    output = _format_directives_for_prompt(directives)
    assert "Stop target" not in output, (
        "Stop target line appeared without after_program\nOutput:\n%s"
        % output)
    assert "VALID PROGRAMS" not in output, (
        "VALID PROGRAMS appeared without after_program\nOutput:\n%s"
        % output)
    print("  PASS")


def test_empty_directives():
    """Empty directives returns empty string."""
    print("Test: empty_directives")
    result = _format_directives_for_prompt({})
    assert result == "", (
        "Empty directives should return '', got %r" % result)
    print("  PASS")


def test_none_directives():
    """None directives returns empty string."""
    print("Test: none_directives")
    result = _format_directives_for_prompt(None)
    assert result == "", (
        "None directives should return '', got %r" % result)
    print("  PASS")


def test_no_stop_conditions():
    """Directives without stop_conditions section work."""
    print("Test: no_stop_conditions")
    directives = {"program_settings": {"phenix.refine": {"strategy": "tls"}}}
    output = _format_directives_for_prompt(directives)
    assert "Stop target" not in output, (
        "Stop target appeared without stop_conditions\nOutput:\n%s" % output)
    assert "phenix.refine: strategy=tls" in output, (
        "program_settings rendering broken\nOutput:\n%s" % output)
    print("  PASS")


# =====================================================================
# SECTION E: Output structural integrity
# =====================================================================

def test_output_is_well_formed():
    """No unresolved placeholders, no Python repr leakage."""
    print("Test: output_is_well_formed")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)
    # No Python format placeholders
    assert "%s" not in output, "Unresolved %%s in output:\n%s" % output
    assert "%d" not in output, "Unresolved %%d in output:\n%s" % output
    # No dict reprs leaked: a dict repr contains both '{' and '}' AND
    # typically a ':' between them.  We check for both braces present
    # together as the signal.
    has_brace_pair = "{" in output and "}" in output
    assert not has_brace_pair, (
        "Possible dict repr leakage (output contains both '{' and '}'):"
        "\n%s" % output)
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """cctbx-style runner."""
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
