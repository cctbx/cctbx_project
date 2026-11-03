"""Tests for v116.11 Stop Condition false-positive fix.

Bug: directive_extractor's `_resolve_after_program` (v115.10)
interpreted the preprocessor-inserted "Stop Condition: None" header
as user stop intent (matching `\\bstop\\b`).  Combined with detection
of "refine" / "build" / "predict" in the Primary Goal text, this
caused AF_7mjs's plan to set:

    after_program: phenix.real_space_refine
    skip_validation: True
    start_with_program: phenix.predict_and_build

…which caused the plan generator to skip Stages 1 (data assessment)
and 2 (map improvement) and start at Stage 3.

The expected behavior (verified by the working March 30, 2026 run)
is to produce empty or minimal directives, letting the full cryo-EM
plan generate (5 stages, starting with mtriage).

Root cause: `_strip_preprocessor_stop_condition` at line 44
required the header to start with literal `stop condition` and
didn't handle the markdown bold + numbered list prefix the
preprocessor actually produces ("7. **Stop Condition**: None").

The fix has three logical changes:

  1. Strengthen the regex that recognizes preprocessor formatting,
     applied consistently at all THREE locations where it appears:
       - `_strip_preprocessor_stop_condition` (the "Stop Condition"
         and "Primary Goal" strip regexes)
       - The `_PREPROCESSOR_SIGNATURES` detection inside the same
         function (gates Primary Goal stripping)
       - The pre-existing `_is_preprocessed` check in
         `extract_directives_simple` (gates `tutorial_patterns`)
     The strengthened regex handles markdown bold (`**Header**:`),
     numbered list prefixes (`7. Header:`), and bullet markers
     (`- Header:`).  Root-cause fix.

  2. Defense-in-depth strip inside `_resolve_after_program`
     before the `\\bstop\\b` check.  Under normal flow this is
     redundant with #1, but catches future code paths that
     bypass the upstream strip.

  3. Suppress `start_with_program` writes from
     `_resolve_after_program` when advice is preprocessed.
     Multi-action signals from descriptive preprocessor prose
     ("Run PredictAndBuild ... dock/trim ... rebuild ... refine")
     should not be treated as user "skip prerequisites"
     prescription.

These tests use the libtbx-stub pattern so they run without a
real PHENIX install.  Tests that target specific resolver
behavior call `_resolve_after_program` directly to isolate from
intent_classifier, which is a separate, legitimate extraction
path that can also write `after_program`.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# =============================================================================
# Path setup
# =============================================================================

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))

for d in (_PROJECT_ROOT, os.path.join(_PROJECT_ROOT, "agent")):
    if d not in sys.path:
        sys.path.insert(0, d)


# =============================================================================
# libtbx stubs (mirror the pattern used by other test files)
# =============================================================================

def _install_stubs():
    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")
    _ensure("libtbx.langchain.knowledge")

    # Stub intent_classifier — returns None so the patterns/resolver
    # paths run unencumbered.  Real intent_classifier behavior is
    # exercised in integration; this file targets the resolver bug.
    ic = _ensure("libtbx.langchain.agent.intent_classifier")
    ic.classify_intent = lambda advice: None


try:
    from libtbx.langchain.agent.directive_extractor import (
        extract_directives_simple, _resolve_after_program,
        _strip_preprocessor_stop_condition)
except ImportError:
    _install_stubs()
    try:
        from libtbx.langchain.agent.directive_extractor import (
            extract_directives_simple, _resolve_after_program,
            _strip_preprocessor_stop_condition)
    except ImportError:
        from directive_extractor import (
            extract_directives_simple, _resolve_after_program,
            _strip_preprocessor_stop_condition)


# =============================================================================
# Tests
# =============================================================================

# The exact processed-advice format the preprocessor produces.
# This is what `extract_directives_simple` actually sees in production.
AF_7MJS_PROCESSED_ADVICE = """1. **Input Files Found**: 7mjs_23883_H.fa, 7mjs_23883_H_1.ccp4, 7mjs_23883_H_2.ccp4, 7mjs_23883_H_full.ccp4

2. **Experiment Type**: cryo-EM (model building/refinement with PredictAndBuild using half-maps and AlphaFold predictions)

3. **Primary Goal**: Run PredictAndBuild (Cryo-EM) for the Fab heavy chain using the sequence and two half-maps, with PDB templates disabled; dock/trim the predicted model, rebuild the problematic loop (residues ~100-120), and refine the model. Optionally use supplied predicted models via "carry_on" when running offline.

4. **Key Parameters**:
- Wavelength: None
- Resolution limit: Not specified (auto-detected from half-maps)
- Number of expected sites: None
- Heavy atom type: None
- Space group: None

5. **Program Parameters**:
use_templates=False
carry_on=True

6. **Special Instructions**:
- Note: Maps/models are boxed and offset from original positions (no action required).

7. **Stop Condition**: None
"""


def test_af7mjs_stop_condition_none_does_not_set_after_program():
    """Primary regression: 'Stop Condition: None' header must not
    trigger after_program assignment.

    Pre-fix, extract_directives_simple returned:
        {"stop_conditions": {"after_program": "phenix.real_space_refine",
                             "skip_validation": True}}
    Post-fix, the returned directives must NOT contain after_program.
    """
    print("Test: af7mjs_stop_condition_none_does_not_set_after_program")
    directives = extract_directives_simple(AF_7MJS_PROCESSED_ADVICE)
    sc = directives.get("stop_conditions", {})

    assert "after_program" not in sc, (
        "AF_7mjs regression: after_program should NOT be set when the "
        "advice only contains 'Stop Condition: None' header.  Got: %r"
        % sc)
    assert "skip_validation" not in sc, (
        "AF_7mjs regression: skip_validation leaked.  Got: %r" % sc)
    print("  PASS")


def test_af7mjs_internal_flag_does_not_leak_to_output():
    """`_set_by_pattern` is an internal marker used by
    tutorial_patterns / denmod_patterns to enable intent-driven
    cleanup.  It must not appear in returned directives.  The
    AF_7mjs path doesn't write this marker (tutorial_patterns are
    skipped for preprocessed advice), but this test guards against
    future regressions if pattern paths get triggered on preprocessed
    advice."""
    print("Test: af7mjs_internal_flag_does_not_leak_to_output")
    directives = extract_directives_simple(AF_7MJS_PROCESSED_ADVICE)
    sc = directives.get("stop_conditions", {})
    assert "_set_by_pattern" not in sc, (
        "_set_by_pattern leaked into output: %r" % sc)
    print("  PASS")


def test_stop_condition_header_with_value_triggers_resolver():
    """v116.x stop_refactor: 'Stop Condition: <real value>' is now
    treated as a positive stop signal by _is_stop_after_requested.
    This is intentional — the field documents user intent ("the user's
    README says stop after phaser"), and the new architecture wants
    that intent captured.

    Pre-stop_refactor, the resolver had inline logic that stripped
    "Stop Condition:" lines and only checked bare \\bstop\\b in the
    remaining text — so a structured Stop Condition value was ignored.
    Post-stop_refactor, the centralized helper _is_stop_after_requested
    explicitly recognizes _STOP_CONDITION_VALUE (anything other than
    None/not specified/N/A/null) as a positive signal.

    This test documents the new behavior so future edits don't
    accidentally revert it.
    """
    print("Test: stop_condition_header_with_value_triggers_resolver")
    advice_lower = """primary goal: do molecular replacement and refine the model.

stop condition: after phaser.
"""
    directives = {}
    _resolve_after_program(directives, advice_lower)
    sc = directives.get("stop_conditions", {})
    # The resolver detects [solve, refine] actions plus the stop signal
    # from "Stop Condition: after phaser." (a value, not None).  Under
    # n>1+stop, after_program = last action = refine.
    assert "after_program" in sc, (
        "after_program should be set when 'Stop Condition: <value>' "
        "signals user-stop intent.  Got: %r" % sc)
    # stop_after_requested should also be True (the resolver sets it
    # whenever it sets after_program from a positive stop signal).
    assert sc.get("stop_after_requested") is True, (
        "stop_after_requested should be True when Stop Condition has "
        "a real value.  Got: %r" % sc)
    print("  PASS")


def test_legitimate_and_stop_still_fires():
    """Sanity: my fix must not break the resolver's detection of real
    stop intent.  Calls _resolve_after_program directly to isolate
    from intent_classifier (which can also write after_program through
    its own legitimate path)."""
    print("Test: legitimate_and_stop_still_fires")
    directives = {}
    _resolve_after_program(directives, "refine and stop")
    sc = directives.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.refine", (
        "Legitimate 'refine and stop' broke: %r" % sc)
    assert sc.get("skip_validation") is True
    print("  PASS")


def test_multi_action_with_legitimate_stop_picks_last():
    """Multi-action + real 'and stop' still picks last action.
    Calls _resolve_after_program directly to isolate from
    intent_classifier."""
    print("Test: multi_action_with_legitimate_stop_picks_last")
    directives = {}
    _resolve_after_program(directives, "phaser, refine, and stop")
    sc = directives.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.refine", (
        "Multi-action stop broke: %r" % sc)
    assert sc.get("start_with_program") == "phenix.phaser"
    print("  PASS")


def test_multi_action_no_stop_clears_after_program():
    """Multi-action without stop sets start_with_program but does
    not write after_program.  Calls _resolve_after_program directly."""
    print("Test: multi_action_no_stop_clears_after_program")
    directives = {}
    _resolve_after_program(directives, "phaser and refine")
    sc = directives.get("stop_conditions", {})
    assert "after_program" not in sc, (
        "Multi-action no-stop should not write after_program: %r" % sc)
    # start_with_program should be the first action (not preprocessed)
    assert sc.get("start_with_program") == "phenix.phaser"
    print("  PASS")


def test_negation_dont_stop_still_works():
    """The 'don't stop' negation guard must continue to work after
    our header-stripping change."""
    print("Test: negation_dont_stop_still_works")
    # "refine" is a single action; "don't stop" negates the stop.
    # _resolve_after_program should NOT set after_program from the
    # action resolver (n=1, no stop → leave as-is).
    # NOTE: tutorial_patterns may set after_program=phenix.refine
    # independently; we only test that _resolve_after_program's
    # negation guard works.  Use _resolve_after_program directly.
    directives = {}
    _resolve_after_program(directives, "don't stop after refinement")
    sc = directives.get("stop_conditions", {})
    assert "after_program" not in sc, (
        "Negation guard broken: 'don't stop after refinement' should "
        "not set after_program; got %r" % sc)
    print("  PASS")


def test_stop_condition_in_middle_of_sentence_is_ok():
    """If 'stop' appears mid-sentence (not as a header), it still
    counts as a stop signal."""
    print("Test: stop_condition_in_middle_of_sentence_is_ok")
    # "I want to refine and stop" — real stop intent, real action
    # ("refine" matches \brefine\b with word boundary)
    directives = {}
    _resolve_after_program(directives, "i want to refine and stop")
    sc = directives.get("stop_conditions", {})
    # Single action (refine) + has_stop → set after_program
    assert sc.get("after_program") == "phenix.refine", (
        "Real stop intent not detected: %r" % sc)
    print("  PASS")


def test_af7mjs_no_start_with_program_for_preprocessed():
    """v116.11 critical: preprocessed advice with multi-action prose
    must NOT set start_with_program.

    The planner uses start_with_program as a skip-prerequisites
    directive.  When the preprocessor describes a workflow with
    multiple actions ("PredictAndBuild ... dock/trim ... rebuild
    ... refine"), it's descriptive prose, not user prescription —
    setting start_with would cause the planner to skip mtriage +
    denmod stages.

    This test reproduces the AF_7mjs symptom directly."""
    print("Test: af7mjs_no_start_with_program_for_preprocessed")
    directives = extract_directives_simple(AF_7MJS_PROCESSED_ADVICE)
    sc = directives.get("stop_conditions", {})
    assert "start_with_program" not in sc, (
        "AF_7mjs regression: start_with_program should NOT be "
        "set from preprocessed multi-action prose.  Setting it "
        "causes the planner to skip prerequisite stages.  Got: %r"
        % sc)
    print("  PASS")


def test_real_user_prose_still_gets_start_with_program():
    """Sanity: real user prose like 'run phaser and refine' (no
    preprocessor signature headers) must still produce
    start_with_program.  This is the v115.10 design intent for
    imperative multi-step user requests.  Calls
    _resolve_after_program directly."""
    print("Test: real_user_prose_still_gets_start_with_program")
    directives = {}
    _resolve_after_program(directives, "run phaser and refine")
    sc = directives.get("stop_conditions", {})
    assert sc.get("start_with_program") == "phenix.phaser", (
        "Real user prose 'run phaser and refine' must produce "
        "start_with_program=phenix.phaser; got %r" % sc)
    print("  PASS")


def test_strip_preprocessor_handles_markdown_bold():
    """Upstream fix: _strip_preprocessor_stop_condition must strip
    markdown-wrapped 'Stop Condition' headers.

    The advice preprocessor produces output with markdown bold
    markers around field names (e.g., '**Stop Condition**: None').
    Pre-fix, the strip regex required the line to start with
    'stop condition' literally, missing all variants the preprocessor
    actually produces.
    """
    print("Test: strip_preprocessor_handles_markdown_bold")
    variants = [
        "Stop Condition: None",
        "**Stop Condition**: None",
        "7. **Stop Condition**: None",
        "  stop condition: stop after autosol",
        "- **Stop Condition**: Stop after refinement",
    ]
    for inp in variants:
        out = _strip_preprocessor_stop_condition(inp)
        assert "stop" not in out.lower(), (
            "Failed to strip '%s' — got %r" % (inp, out))
    print("  PASS (%d variants stripped)" % len(variants))


def test_strip_preprocessor_preserves_user_prose():
    """Upstream fix must not strip real user advice that happens to
    contain 'stop condition' as part of natural prose."""
    print("Test: strip_preprocessor_preserves_user_prose")
    preserved = [
        "please stop condition processing",
        "the stop condition is map quality > 0.7",
        "stop refining when CC plateaus",
    ]
    for inp in preserved:
        out = _strip_preprocessor_stop_condition(inp)
        assert out == inp, (
            "Should NOT have stripped '%s' — got %r" % (inp, out))
    print("  PASS (%d prose preserved)" % len(preserved))


def test_strip_preprocessor_handles_primary_goal_markdown():
    """The Primary Goal strip path was also strengthened — verify it
    handles the markdown-wrapped variant in preprocessed contexts."""
    print("Test: strip_preprocessor_handles_primary_goal_markdown")
    # Must include a preprocessor signature so Primary Goal stripping
    # is enabled (the function's heuristic).
    preprocessed = """1. **Input Files Found**: foo.mtz

3. **Primary Goal**: do refinement
"""
    out = _strip_preprocessor_stop_condition(preprocessed)
    # 'primary goal' should be stripped
    assert "primary goal" not in out.lower(), (
        "Primary Goal not stripped: %r" % out)
    print("  PASS")


def test_extract_directives_simple_recognizes_markdown_preprocessor():
    """The pre-existing _is_preprocessed check in
    extract_directives_simple was using a regex that required
    signature headers to start with the literal word (no list marker,
    no markdown).  Pre-fix, AF_7mjs-style markdown-formatted advice
    would NOT be recognized as preprocessed, causing
    tutorial_patterns to incorrectly run on descriptive prose.

    This test verifies the strengthened regex by checking the strip
    function's preprocessor-detection behavior directly:
    `_strip_preprocessor_stop_condition` only strips 'Primary Goal:'
    when it detects preprocessor signatures.  If markdown signatures
    are correctly recognized, Primary Goal gets stripped.
    """
    print("Test: extract_directives_simple_recognizes_markdown_preprocessor")
    # Markdown-formatted signatures (the AF_7mjs format)
    advice = """1. **Input Files Found**: half_maps.ccp4

3. **Primary Goal**: Map to model using cryo-EM half-maps.
"""
    out = _strip_preprocessor_stop_condition(advice)
    # If signature detection recognized this as preprocessed,
    # Primary Goal should be stripped.  Pre-fix, the old regex
    # wouldn't match the markdown signature → not preprocessed →
    # Primary Goal preserved.
    assert "primary goal" not in out.lower(), (
        "Markdown preprocessor signatures must be recognized so that "
        "Primary Goal stripping fires.  Got: %r" % out)
    print("  PASS")


def test_strip_alone_handles_af7mjs_stop_condition_header():
    """Smaller regression check: the upstream strip function alone
    must remove the AF_7mjs-style "**Stop Condition**: None" header.

    This isolates the upstream-strip fix from the rest of the
    extraction pipeline.  If this test fails, the regex in
    _strip_preprocessor_stop_condition has regressed.
    """
    print("Test: strip_alone_handles_af7mjs_stop_condition_header")
    # Exact format the preprocessor produces (from ok.txt, March 30)
    advice = "7. **Stop Condition**: None\n"
    out = _strip_preprocessor_stop_condition(advice)
    assert "stop" not in out.lower(), (
        "AF_7mjs regression: '7. **Stop Condition**: None' must be "
        "stripped by _strip_preprocessor_stop_condition.  Got: %r"
        % out)
    print("  PASS")


# =============================================================================
# Test runner (matches the project's tst_*.py convention)
# =============================================================================

def run_all_tests():
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
