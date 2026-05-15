"""Tests for `_format_directives_for_prompt` after_program rendering.

This test file asserts the CURRENT wording in
`knowledge/prompts_hybrid.py` `_format_directives_for_prompt`,
which renders `stop_conditions.after_program` as three lines:

    - Stop after <PROG> completes
    - **CRITICAL: You MUST run <PROG> before stopping. If it's in
      VALID PROGRAMS, choose it NOW.**
    - Do NOT keep running refinement cycles - run <PROG> instead!

The rendering is correct and intentional given the current architecture:

* The "If it's in VALID PROGRAMS" qualifier scopes the directive to
  programs the workflow engine has actually offered.
* v116.17 ensures `valid_programs` is NOT incorrectly wiped to
  `["STOP"]` at the validate step when an after_program directive
  is set with after_program_done=True, so the qualifier is
  enforceable: when validation is pending, the LLM sees validation
  programs in VALID PROGRAMS and does not need to pick the (already
  done) after_program.

Historical note: there was an earlier proposal ("Phase 6a") to
replace this wording with "Stop target: <PROG>." framing.  The
v116.x family did not adopt that change.  These tests document
what IS in the code, not what could be there.

These are STATIC regression tests against accidental wording
changes.  They do NOT verify LLM behavior — that requires end-to-
end runs against the tutorial corpus.
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
# SECTION A: Current after_program rendering
# =====================================================================

def test_current_wording_present_in_after_program_block():
    """The current three-line directive block must appear.

    Each rendered line must be present verbatim with the after_program
    name substituted.  These assertions lock the existing wording so
    accidental edits to `_format_directives_for_prompt` are caught.
    """
    print("Test: current_wording_present_in_after_program_block")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)

    expected_phrases = [
        # Line 1: the basic "stop after X completes"
        "Stop after phenix.predict_and_build completes",
        # Line 2: the CRITICAL "must run before stopping" directive
        # — note the in-list qualifier and the bold-italic markdown
        "CRITICAL: You MUST run phenix.predict_and_build before stopping.",
        "If it's in VALID PROGRAMS, choose it NOW.",
        # Line 3: the "don't keep refining" anti-loop directive
        "Do NOT keep running refinement cycles - run "
        "phenix.predict_and_build instead!",
    ]
    for phrase in expected_phrases:
        assert phrase in output, (
            "Current wording missing: '%s'\nFull output:\n%s"
            % (phrase, output))
    print("  PASS")


def test_no_unresolved_placeholders_in_after_program_block():
    """No '%s' or '%d' should leak into the rendered output.

    All three template lines use `%s` for the program name; if any
    formatting step fails the placeholder leaks through.  Catch that.
    """
    print("Test: no_unresolved_placeholders_in_after_program_block")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)
    assert "%s" not in output, \
        "Unresolved %%s in output:\n%s" % output
    assert "%d" not in output, \
        "Unresolved %%d in output:\n%s" % output
    print("  PASS")


def test_phase6a_wording_not_present():
    """Phase 6a wording was discussed but not adopted.

    If someone is running these tests against a future codebase
    where the Phase 6a rewording HAS been applied, this test will
    fail and they should switch to the Phase 6a test file
    (preserved in the v116.10 phase6a workspace).  Until then,
    these phrases must NOT be present.
    """
    print("Test: phase6a_wording_not_present")
    directives = {
        "stop_conditions": {"after_program": "phenix.predict_and_build"}
    }
    output = _format_directives_for_prompt(directives)

    phase6a_phrases = [
        "Stop target: phenix.predict_and_build.",
        "If phenix.predict_and_build is NOT in VALID PROGRAMS",
        "The workflow will reach phenix.predict_and_build",
        "Never pick a program outside VALID PROGRAMS",
    ]
    for phrase in phase6a_phrases:
        assert phrase not in output, (
            "Phase 6a wording detected — either Phase 6a has been "
            "applied (switch to the Phase 6a test file) or this "
            "test file is out of date.\n"
            "Phrase found: '%s'\nFull output:\n%s" % (phrase, output))
    print("  PASS")


# =====================================================================
# SECTION B: Substitution correctness across program names
# =====================================================================

def test_substitution_with_different_program_names():
    """All three template lines must substitute every program name.

    Verifies that `%s` substitution happens on each of the three
    rendered lines for a variety of canonical program names.
    """
    print("Test: substitution_with_different_program_names")
    for prog in [
        "phenix.refine",
        "phenix.molprobity",
        "phenix.real_space_refine",
        "phenix.ligandfit",
        "phenix.predict_and_build",
    ]:
        directives = {"stop_conditions": {"after_program": prog}}
        output = _format_directives_for_prompt(directives)
        for phrase in [
            "Stop after %s completes" % prog,
            "CRITICAL: You MUST run %s before stopping." % prog,
            "run %s instead!" % prog,
        ]:
            assert phrase in output, (
                "Phrase '%s' did not render for %s\nOutput:\n%s"
                % (phrase, prog, output))
        assert "%s" not in output, (
            "Unresolved %%s for %s in output:\n%s" % (prog, output))
    print("  PASS")


def test_in_list_qualifier_present():
    """The 'If it's in VALID PROGRAMS' qualifier is present.

    This qualifier is what makes the directive scope-correct: the
    LLM should only pick after_program when it is actually
    available.  v116.17 ensures the workflow engine doesn't
    incorrectly remove validation programs from VALID PROGRAMS
    at the validate step, so this qualifier is enforceable.
    """
    print("Test: in_list_qualifier_present")
    directives = {
        "stop_conditions": {"after_program": "phenix.refine"}
    }
    output = _format_directives_for_prompt(directives)
    assert "If it's in VALID PROGRAMS" in output, (
        "The in-list qualifier 'If it's in VALID PROGRAMS' is missing.\n"
        "Output:\n%s" % output)
    print("  PASS")


# =====================================================================
# SECTION C: Adjacent stop_conditions fields unaffected
# =====================================================================

def test_after_cycle_still_works():
    """`after_cycle` rendering renders independently of after_program."""
    print("Test: after_cycle_still_works")
    directives = {"stop_conditions": {"after_cycle": 5}}
    output = _format_directives_for_prompt(directives)
    assert "Stop after cycle 5" in output, (
        "after_cycle rendering broken\nOutput:\n%s" % output)
    # Should NOT have any after_program content
    assert "CRITICAL: You MUST run" not in output, (
        "after_cycle-only directives should not produce the "
        "after_program CRITICAL line\nOutput:\n%s" % output)
    print("  PASS")


def test_both_after_cycle_and_after_program():
    """Both stop conditions coexist; both render."""
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
    assert "Stop after phenix.refine completes" in output, (
        "Stop-after-program line missing when both present\nOutput:\n%s"
        % output)
    assert "CRITICAL: You MUST run phenix.refine before stopping." in output, (
        "after_program CRITICAL line missing when both present\n"
        "Output:\n%s" % output)
    print("  PASS")


def test_other_stop_condition_fields_preserved():
    """r_free_target, map_cc_target, max_refine_cycles, skip_validation
    all render alongside after_program.
    """
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
        "Stop after phenix.refine completes",
        "CRITICAL: You MUST run phenix.refine before stopping.",
        "Target R-free: 0.220",
        "Target map CC: 0.75",
        "Maximum 5 refinement cycles",
        "Validation can be skipped before stopping",
    ]
    for phrase in expected:
        assert phrase in output, (
            "Expected phrase missing: '%s'\nOutput:\n%s"
            % (phrase, output))
    print("  PASS")


# =====================================================================
# SECTION D: Edge cases
# =====================================================================

def test_no_after_program_field():
    """No after_program → no after_program block at all."""
    print("Test: no_after_program_field")
    directives = {"stop_conditions": {"after_cycle": 3}}
    output = _format_directives_for_prompt(directives)
    assert "CRITICAL: You MUST run" not in output, (
        "after_program CRITICAL line appeared without after_program\n"
        "Output:\n%s" % output)
    assert "Do NOT keep running refinement cycles" not in output, (
        "after_program anti-loop line appeared without after_program\n"
        "Output:\n%s" % output)
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
    """Directives without stop_conditions section work without errors."""
    print("Test: no_stop_conditions")
    directives = {"program_settings": {"phenix.refine": {"strategy": "tls"}}}
    output = _format_directives_for_prompt(directives)
    assert "CRITICAL: You MUST run" not in output, (
        "after_program CRITICAL line appeared without stop_conditions\n"
        "Output:\n%s" % output)
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
    # No dict reprs leaked
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
