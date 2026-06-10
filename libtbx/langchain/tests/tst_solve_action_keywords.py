"""K_H14_ITEM_1: solve action keyword cleanup (v119.H14).

Covers the phaser false-positive surfaced by run_39_openai batch
analysis of the 1029B-sad__rules_only_stop dataset.

The rules-only directive extractor's _ACTION_TABLE["solve"] entry
treated the goal phrase "solve the structure" as a synonym for
"do molecular replacement with phaser."  When a README contained
BOTH that phrase AND another action ("Stop after refinement"),
the multi-action branch in _apply_workflow_intent_fallback set
start_with_program=phenix.phaser — forcing phaser into the
workflow even on SAD/MAD datasets where the correct method is
autosol.

H14 fix: remove "solve the structure" and "solve structure" from
the keyword list.  Explicit method keywords ("molecular
replacement", "phaser", "mr ") still match — only the ambiguous
goal phrase is removed.

This file pins:
  - The two removed keywords are gone (source scan)
  - _detect_actions no longer matches the goal phrase
  - Explicit method keywords still match
  - The 1029B-sad README case produces correct directives
  - Negative cases: legitimate MR + stop still gets start_with_program

Total: 13 tests.
"""
from __future__ import absolute_import, division, print_function

import os
import sys


# =====================================================================
# Import helpers — sandbox-friendly
# =====================================================================

def _try_import_de():
    """Import the directive_extractor module itself (not just funcs)."""
    try:
        from libtbx.langchain.agent import directive_extractor as de
        return de, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from agent import directive_extractor as de
        return de, None
    except ImportError as e:
        return None, str(e)


def _read_source(rel_path):
    """Read a source file from the ship layout for source-scan tests."""
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    path = os.path.join(parent, rel_path)
    with open(path, 'r') as f:
        return f.read()


# =====================================================================
# §A: Source-scan — the removed keywords must NOT be back
# =====================================================================

def test_source_scan_solve_keywords_lack_goal_phrase():
    """Pin the H14 fix at source level: the goal phrases must not
    appear in the `solve` action's keyword list."""
    src = _read_source("agent/directive_extractor.py")

    # Locate the solve action block.
    # We can't just grep for the strings because they may legitimately
    # appear in comments elsewhere — locate the action table entry
    # specifically.
    solve_start = src.find('"solve": {')
    assert solve_start != -1, "Could not find solve action in _ACTION_TABLE"

    # Find the closing brace of the solve entry
    solve_end = src.find('    },', solve_start)
    assert solve_end != -1, "Could not find end of solve action"
    solve_block = src[solve_start:solve_end]

    # The keywords list inside the block must not contain the goal
    # phrases as items (note the quoting matters — they could appear
    # inside the H14 comment explaining WHY they were removed)
    keywords_start = solve_block.find('"keywords":')
    assert keywords_start != -1, "solve action must have keywords"
    keywords_section = solve_block[keywords_start:]

    # Strip out comment lines (lines starting with whitespace + #)
    # so we only look at actual code
    keywords_no_comments = "\n".join(
        line for line in keywords_section.split("\n")
        if not line.lstrip().startswith("#")
    )

    assert '"solve the structure"' not in keywords_no_comments, (
        '"solve the structure" must not appear as a `solve` action '
        'keyword (H14 fix — see batch analysis of 1029B-sad)')
    assert '"solve structure"' not in keywords_no_comments, (
        '"solve structure" must not appear as a `solve` action '
        'keyword (H14 fix — see batch analysis of 1029B-sad)')
    print("  PASS: test_source_scan_solve_keywords_lack_goal_phrase")


def test_source_scan_solve_keywords_keep_explicit_methods():
    """The cleanup must not over-correct: explicit method keywords
    must still trigger `solve`."""
    src = _read_source("agent/directive_extractor.py")
    solve_start = src.find('"solve": {')
    solve_end = src.find('    },', solve_start)
    solve_block = src[solve_start:solve_end]
    keywords_start = solve_block.find('"keywords":')
    keywords_section = solve_block[keywords_start:]
    keywords_no_comments = "\n".join(
        line for line in keywords_section.split("\n")
        if not line.lstrip().startswith("#")
    )

    # These method names MUST still be present
    assert '"molecular replacement"' in keywords_no_comments
    assert '"phaser"' in keywords_no_comments
    assert '"mr "' in keywords_no_comments
    print("  PASS: test_source_scan_solve_keywords_keep_explicit_methods")


# =====================================================================
# §B: _detect_actions behavior
# =====================================================================

def test_detect_actions_goal_phrase_alone_no_solve():
    """v119.H14: 'solve the structure' alone (no method mentioned)
    does NOT trigger the solve action.  Pre-H14 it did, which on
    SAD/MAD data caused phaser to be injected as start_with_program."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    actions = de._detect_actions("solve the structure by standard procedures")
    action_names = [a[0] for a in actions]
    assert "solve" not in action_names, (
        "Pre-H14 'solve the structure' triggered solve action. "
        "Got actions: %r" % action_names)
    print("  PASS: test_detect_actions_goal_phrase_alone_no_solve")


def test_detect_actions_goal_phrase_with_stop_no_solve():
    """The 1029B-sad README pattern: goal phrase + 'Stop after
    refinement'.  Should detect ONLY refine, not solve."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    advice = ("solve the structure by standard procedures. "
              "anomalous scatterer is selenium. stop after refinement")
    actions = de._detect_actions(advice)
    action_names = [a[0] for a in actions]
    assert "solve" not in action_names, (
        "Pre-H14 this matched solve via 'solve the structure'. "
        "Got: %r" % action_names)
    assert "refine" in action_names, (
        "'refinement' should still match refine action. Got: %r"
        % action_names)
    print("  PASS: test_detect_actions_goal_phrase_with_stop_no_solve")


def test_detect_actions_explicit_mr_still_triggers_solve():
    """Explicit 'molecular replacement' MUST still trigger solve."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    actions = de._detect_actions(
        "use molecular replacement to find the model")
    action_names = [a[0] for a in actions]
    assert "solve" in action_names, (
        "Explicit 'molecular replacement' must trigger solve. "
        "Got: %r" % action_names)
    print("  PASS: test_detect_actions_explicit_mr_still_triggers_solve")


def test_detect_actions_explicit_phaser_still_triggers_solve():
    """Explicit 'phaser' MUST still trigger solve."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    actions = de._detect_actions("run phaser to place the model")
    action_names = [a[0] for a in actions]
    assert "solve" in action_names, (
        "Explicit 'phaser' must trigger solve. Got: %r" % action_names)
    print("  PASS: test_detect_actions_explicit_phaser_still_triggers_solve")


def test_detect_actions_mr_token_still_triggers_solve():
    """Explicit 'mr ' (with trailing space to avoid false-positive
    inside 'remove' etc.) MUST still trigger solve."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    actions = de._detect_actions("do mr  with the search model")
    action_names = [a[0] for a in actions]
    assert "solve" in action_names
    print("  PASS: test_detect_actions_mr_token_still_triggers_solve")


# =====================================================================
# §C: End-to-end directive resolution
# =====================================================================

def test_e2e_1029B_sad_pattern_no_start_with_phaser():
    """End-to-end test of Tom's case: the 1029B-sad README pattern
    must produce after_program=phenix.refine WITHOUT
    start_with_program=phenix.phaser."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    advice_lower = (
        "solve the structure by standard procedures. "
        "anomalous scatterer is selenium. "
        "stop after refinement"
    )
    directives = {}
    de._apply_workflow_intent_fallback(directives, advice_lower)

    sc = directives.get("stop_conditions", {})
    assert sc.get("after_program") == "phenix.refine", (
        "after_program should be phenix.refine (set by n==1+stop "
        "branch).  Got: %r" % sc)
    assert "start_with_program" not in sc, (
        "Pre-H14 this incorrectly set start_with_program=phenix.phaser "
        "from 'solve the structure' triggering n>1 branch.  "
        "Got: %r" % sc)
    print("  PASS: test_e2e_1029B_sad_pattern_no_start_with_phaser")


def test_e2e_explicit_mr_plus_stop_sets_start_with():
    """The cleanup must not over-correct: explicit MR + stop SHOULD
    still set start_with_program=phenix.phaser (legitimate two-action
    user request)."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    advice_lower = "run molecular replacement, then refine and stop"
    directives = {}
    de._apply_workflow_intent_fallback(directives, advice_lower)

    sc = directives.get("stop_conditions", {})
    assert sc.get("start_with_program") == "phenix.phaser", (
        "Explicit MR + stop SHOULD set start_with_program. "
        "Got: %r" % sc)
    assert sc.get("after_program") == "phenix.refine"
    print("  PASS: test_e2e_explicit_mr_plus_stop_sets_start_with")


# =====================================================================
# §D: Additional regression cases — verify H14 doesn't break other
# README phrasings.
# =====================================================================

def test_e2e_other_solve_goal_phrasings():
    """Tom-style README goal phrasings beyond '1029B-sad__rules_only_stop'.
    All of these should NOT set start_with_program=phenix.phaser since
    they don't explicitly request MR.

    NOTE: a README that says "solve the structure BY <method>" should
    set start_with_program=<method>'s program if <method> is a
    method keyword — that's the correct behavior.  This test only
    covers cases without an explicit method keyword.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    goal_phrasings = [
        # The 1029B-sad case
        "solve the structure by standard procedures. stop after refinement",
        # Same idea, different surrounding text
        "the goal is to solve the structure. stop after refinement",
        "solve structure. stop after refinement",
        # Goal phrase nested inside larger advice
        "this tutorial demonstrates how to solve the structure. "
        "stop after refinement",
        # Only "solve" with one refinement (no stop signal in the sense
        # of "stop after X" — just a single refinement cycle)
        # Note: "with one refinement" doesn't trigger _is_stop_after_requested
        "solve the structure with one refinement",
    ]
    for advice in goal_phrasings:
        directives = {}
        de._apply_workflow_intent_fallback(directives, advice.lower())
        sc = directives.get("stop_conditions", {})
        # The critical check: no phaser start_with_program (which was
        # the H14 bug — adding the wrong method to a goal-only phrase)
        start = sc.get("start_with_program", "")
        assert start != "phenix.phaser", (
            "Goal-phrasing %r should NOT set start_with_program=phaser. "
            "Got: %r" % (advice, sc))
    print("  PASS: test_e2e_other_solve_goal_phrasings")


def test_e2e_solve_with_explicit_method_uses_method():
    """When 'solve the structure' is paired with an explicit method
    keyword (autosol, phaser, predict_and_build, etc.), the method
    drives — not the bare goal verb.

    This is the GOOD outcome of H14: removing the goal-phrase
    match means the method keyword now wins.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # solve + experimental phasing → autosol
    actions = de._detect_actions(
        "solve the structure by experimental phasing. stop after refinement")
    names = [a[0] for a in actions]
    assert "phase" in names, (
        "'experimental phasing' should trigger `phase` (autosol). "
        "Got: %r" % names)
    assert "solve" not in names, (
        "Bare 'solve the structure' should not trigger `solve`. "
        "Got: %r" % names)

    # solve + phaser → solve (phaser)
    actions = de._detect_actions(
        "solve the structure with phaser. stop after refinement")
    names = [a[0] for a in actions]
    assert "solve" in names, (
        "Explicit 'phaser' should trigger `solve`. "
        "Got: %r" % names)

    # solve + autosol → phase
    actions = de._detect_actions(
        "solve the structure using autosol. stop after refinement")
    names = [a[0] for a in actions]
    assert "phase" in names, (
        "Explicit 'autosol' should trigger `phase`. "
        "Got: %r" % names)
    print("  PASS: test_e2e_solve_with_explicit_method_uses_method")


def test_e2e_explicit_mr_at_start_of_advice():
    """Variant of test_e2e_explicit_mr_plus_stop_sets_start_with —
    the explicit MR keyword can appear in different positions."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    advice_variants = [
        "phaser then refine and stop",
        "use phaser to solve, then refine and stop",
        "molecular replacement followed by refinement, stop after",
        # mr (with trailing space) should also work
        "mr  then refine and stop",
    ]
    for advice in advice_variants:
        directives = {}
        de._apply_workflow_intent_fallback(directives, advice.lower())
        sc = directives.get("stop_conditions", {})
        assert sc.get("start_with_program") == "phenix.phaser", (
            "Advice with explicit MR keyword %r should set "
            "start_with_program=phenix.phaser. Got: %r"
            % (advice, sc))
    print("  PASS: test_e2e_explicit_mr_at_start_of_advice")


def test_e2e_workflow_engine_solve_intent_unaffected():
    """v119.H14 modifies _ACTION_TABLE only.  The separate
    classify_intent path (which uses agent/intent_classifier.py)
    is not affected.  Source-scan based.
    """
    src = _read_source("agent/directive_extractor.py")

    # Verify the intent='solve' branches reference classify_intent
    assert "classify_intent" in src, (
        "directive_extractor.py should still import classify_intent "
        "from agent.intent_classifier")
    # The two _intent == "solve" branches should still exist
    assert src.count('_intent == "solve"') >= 2, (
        "Both _intent == 'solve' branches (at lines ~822 and ~4604) "
        "should still be present")
    print("  PASS: test_e2e_workflow_engine_solve_intent_unaffected")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: Source scans (2)
    test_source_scan_solve_keywords_lack_goal_phrase()
    test_source_scan_solve_keywords_keep_explicit_methods()
    # §B: _detect_actions behavior (5)
    test_detect_actions_goal_phrase_alone_no_solve()
    test_detect_actions_goal_phrase_with_stop_no_solve()
    test_detect_actions_explicit_mr_still_triggers_solve()
    test_detect_actions_explicit_phaser_still_triggers_solve()
    test_detect_actions_mr_token_still_triggers_solve()
    # §C: End-to-end (2)
    test_e2e_1029B_sad_pattern_no_start_with_phaser()
    test_e2e_explicit_mr_plus_stop_sets_start_with()
    # §D: Additional regression cases (4)
    test_e2e_other_solve_goal_phrasings()
    test_e2e_solve_with_explicit_method_uses_method()
    test_e2e_explicit_mr_at_start_of_advice()
    test_e2e_workflow_engine_solve_intent_unaffected()


if __name__ == "__main__":
    print("K_H14_ITEM_1: solve action keyword cleanup (v119.H14)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H14_ITEM_1 complete.")
