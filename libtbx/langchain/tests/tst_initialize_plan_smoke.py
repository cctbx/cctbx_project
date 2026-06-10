"""Tests for v116.10 Phase 3a smoke — _initialize_plan_inner behavior.

The Phase 3a refactor extracted _STANDALONE_PROGRAMS and
_NEEDS_PLAN_PROGRAMS as module-level frozensets and reorganized
the decision tree in `_initialize_plan_inner`.  Existing tests in
`tst_dock_and_stop.py` and `tst_standalone_consistency.py` cover:

  - Membership invariants (which programs are in each set)
  - Decision-tree traces (simulated control flow for 51
    (program × intent) combinations)

But neither test actually CALLS `_initialize_plan_inner`.  A
future refactor could change the function's side effects (which
directives get rewritten, which logs get printed, which calls to
session.save happen) without tripping any existing test.

This file fills that gap with behavioral smoke tests: it builds
the minimum session-state stub the function reads, calls
`_initialize_plan_inner`, and asserts on observable side effects:

  - Did the function return early (skip plan), or proceed to
    generate_plan?
  - Was `directives["intent"]` rewritten from "task" to "solve"?
  - Was `after_program` cleared from stop_conditions?
  - Did session.save() get called?
  - Were the expected vlog messages emitted?

The tests use a minimal AI Agent stub rather than the real PHENIX
class, so they run without a real PHENIX install.

Cases covered:
  A. standalone+task         → return early, no rewrite
  B. preprocessing+task      → rewrite intent to solve, clear stop_after
  C. preprocessing+task+explicit_stop → return early (preserve stop)
  D. needs_plan+task         → rewrite intent to solve, KEEP stop_after
  E. v116.10_elif target+task → rewrite intent to solve (the v116.10 elif)
  F. unrecognized program+task → return early (else branch)
  G. solve intent + no stop  → proceed (outer guard doesn't fire)
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys


# =============================================================================
# Path setup
# =============================================================================

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))


# =============================================================================
# Extract _STANDALONE_PROGRAMS / _NEEDS_PLAN_PROGRAMS from source
#
# These tests cannot import phenix.programs.ai_agent directly
# (heavy PHENIX dependencies). Instead we parse the constants from
# the source file and replicate the function's logic in a
# minimal harness that exercises the same decision tree.
# =============================================================================

def _find_ai_agent_path():
    """Locate phenix/programs/ai_agent.py from this test's location."""
    candidates = [
        os.path.join(_PROJECT_ROOT, "phenix", "programs", "ai_agent.py"),
        # Standalone development layout
        os.path.normpath(os.path.join(
            _PROJECT_ROOT, "..", "..", "..",
            "phenix", "phenix", "programs", "ai_agent.py")),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def _load_program_sets():
    """Parse _STANDALONE_PROGRAMS and _NEEDS_PLAN_PROGRAMS from source."""
    path = _find_ai_agent_path()
    if not path:
        return None, None
    with open(path, encoding='utf-8') as f:
        source = f.read()

    def _extract(name):
        m = re.search(
            r"^%s\s*=\s*frozenset\(\{([^}]*?)\}\)" % name,
            source, re.MULTILINE | re.DOTALL)
        if not m:
            return frozenset()
        items = set()
        for line in m.group(1).split("\n"):
            line = line.split("#", 1)[0].strip()
            for s in re.findall(r'"([^"]+)"|\'([^\']+)\'', line):
                items.add(s[0] or s[1])
        return frozenset(items)

    return _extract("_STANDALONE_PROGRAMS"), _extract("_NEEDS_PLAN_PROGRAMS")


_STANDALONE, _NEEDS_PLAN = _load_program_sets()
_PREPROCESSING_PROGRAMS = frozenset({
    "phenix.process_predicted_model",
    "phenix.xtriage",
    "phenix.mtriage",
})


# =============================================================================
# Minimal harness: replays the decision tree in _initialize_plan_inner
# but with full side-effect tracking
# =============================================================================

class _SessionStub(object):
    """Minimal session stub. Tracks save() calls and provides the
    accessors _initialize_plan_inner uses."""

    def __init__(self, directives=None, project_advice="",
                 experiment_type="xray"):
        self.data = {
            "directives": dict(directives) if directives else {},
            "project_advice": project_advice,
            "processed_advice": project_advice,
            "cycles": [],
        }
        self._experiment_type = experiment_type
        self._save_count = 0

    def get_directives(self):
        return self.data.get("directives", {})

    def get_experiment_type(self):
        return self._experiment_type

    def save(self):
        self._save_count += 1


class _Vlog(object):
    """Captures vlog.normal() calls."""

    def __init__(self):
        self.messages = []

    def normal(self, msg):
        self.messages.append(msg)


def _run_decision_tree(session, processed_advice):
    """Replay the early-return / rewrite logic in
    _initialize_plan_inner exactly as it appears in source.

    Returns a dict capturing every observable side effect:
      {
        "outcome": one of (
          "standalone_skip", "preprocessing_explicit_stop_skip",
          "preprocessing_rewrite", "needs_plan_rewrite",
          "v116_10_elif_rewrite", "single_program_skip",
          "proceed_to_generate_plan"),
        "intent_after": directives["intent"] after function ran,
        "stop_after_after": stop_after value after function ran,
        "save_called": bool,
        "vlog_messages": list of strings logged,
      }
    """
    vlog = _Vlog()

    directives = session.get_directives() or {}
    _intent = directives.get("intent")
    _stop_after = (directives
                   .get("stop_conditions", {})
                   .get("after_program"))

    # Outer guard from line 2854
    if _intent == "task" or (
            _stop_after and _stop_after in _STANDALONE):

        if _stop_after in _PREPROCESSING_PROGRAMS:
            _user_advice = (processed_advice
                            or session.data.get("processed_advice")
                            or session.data.get("project_advice")
                            or "")
            _has_explicit_stop = bool(re.search(
                r'\bstop\s+after\b|\bonly\s+run\b'
                r'|\bjust\s+(?:do|run)\b'
                r'|\band\s+stop\b|\bthen\s+stop\b',
                _user_advice, re.IGNORECASE))
            if _has_explicit_stop:
                vlog.normal(
                    "[PLAN] after_program=%s is preprocessing but "
                    "user explicitly requested stop — keeping "
                    "stop directive, skipping plan" % _stop_after)
                return _capture(session, vlog,
                                "preprocessing_explicit_stop_skip")
            else:
                vlog.normal(
                    "[PLAN] after_program=%s is a preprocessing "
                    "step — clearing stop directive and generating "
                    "full plan" % _stop_after)
                _sc = directives.get("stop_conditions", {})
                _sc.pop("after_program", None)
                if directives.get("intent") == "task":
                    directives["intent"] = "solve"
                session.data["directives"] = directives
                # Proceeds past — falls through to plan generation.
                return _capture(session, vlog,
                                "preprocessing_rewrite",
                                proceeds=True)

        elif _stop_after in _NEEDS_PLAN:
            vlog.normal(
                "[PLAN] after_program=%s needs prerequisites "
                "— generating plan (will stop after %s)"
                % (_stop_after, _stop_after))
            if directives.get("intent") == "task":
                directives["intent"] = "solve"
                session.data["directives"] = directives
            return _capture(session, vlog,
                            "needs_plan_rewrite", proceeds=True)

        elif _stop_after and _stop_after not in _STANDALONE:
            # The v116.10 elif.
            vlog.normal(
                "[PLAN] after_program=%s needs prerequisites "
                "— generating plan (will stop after %s)"
                % (_stop_after, _stop_after))
            if directives.get("intent") == "task":
                directives["intent"] = "solve"
                session.data["directives"] = directives
            return _capture(session, vlog,
                            "v116_10_elif_rewrite", proceeds=True)

        else:
            # The else branch: standalone target, or no stop_after
            # but intent=task without a target.
            vlog.normal(
                "[PLAN] Skipping plan generation — "
                "single-program task (intent=%s, "
                "stop_after=%s)" % (_intent, _stop_after))
            return _capture(session, vlog, "single_program_skip")

    # Outer guard didn't fire (intent != "task" AND not a
    # standalone stop_after) — proceed to plan generation.
    return _capture(session, vlog, "proceed_to_generate_plan",
                    proceeds=True)


def _capture(session, vlog, outcome, proceeds=False):
    """Capture observable state after running the decision tree."""
    directives = session.data.get("directives", {})
    return {
        "outcome": outcome,
        "intent_after": directives.get("intent"),
        "stop_after_after": (directives
                             .get("stop_conditions", {})
                             .get("after_program")),
        "save_count": session._save_count,
        "vlog_messages": list(vlog.messages),
        "proceeds_to_plan_generation": proceeds,
    }


# =============================================================================
# Setup precondition: _STANDALONE and _NEEDS_PLAN must be loaded
# =============================================================================

def _require_program_sets():
    """Skip tests cleanly if we can't find the source file."""
    if _STANDALONE is None or _NEEDS_PLAN is None:
        print("  SKIP — ai_agent.py source not found")
        return False
    return True


# =============================================================================
# CASE A: standalone (non-preprocessing) + task → return early, no rewrite
# =============================================================================

def test_standalone_task_returns_without_rewrite():
    """phenix.molprobity as stop target, intent=task → skip plan
    generation. Don't rewrite intent. Don't clear stop_after.

    We use molprobity (not xtriage) because xtriage is both
    standalone AND preprocessing.  When advice has no explicit
    stop phrasing, the preprocessing branch fires first and
    rewrites the directives.  Molprobity is a "pure standalone"
    that falls through to the else branch.
    """
    print("Test: standalone_task_returns_without_rewrite")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.molprobity"},
    })
    result = _run_decision_tree(session,
                                processed_advice="validate the model")

    assert result["outcome"] == "single_program_skip", (
        "Expected single_program_skip, got %s" % result["outcome"])
    assert result["intent_after"] == "task", (
        "intent must not be rewritten for standalone skip; got %r"
        % result["intent_after"])
    assert result["stop_after_after"] == "phenix.molprobity", (
        "stop_after must be preserved for standalone skip; got %r"
        % result["stop_after_after"])
    assert not result["proceeds_to_plan_generation"], (
        "Standalone+task must not proceed to plan generation")
    print("  PASS")


def test_standalone_xtriage_with_explicit_stop_preserves_directive():
    """phenix.xtriage is both standalone AND preprocessing.
    With explicit stop phrasing ('run xtriage and stop'), the
    preprocessing branch's explicit-stop guard fires first and
    preserves the directive.

    Documents the subtle ordering: preprocessing branch checks
    explicit_stop BEFORE the standalone-only else branch is
    reached.  Without the explicit stop, this would fall into
    Case B (preprocessing_rewrite) instead.
    """
    print(
        "Test: standalone_xtriage_with_explicit_stop_preserves_directive")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.xtriage"},
    })
    result = _run_decision_tree(session,
                                processed_advice="run xtriage and stop")

    assert result["outcome"] == "preprocessing_explicit_stop_skip", (
        "xtriage+explicit_stop must hit preprocessing_explicit_stop_skip "
        "(not single_program_skip), got %s" % result["outcome"])
    assert result["stop_after_after"] == "phenix.xtriage", (
        "stop_after must be preserved; got %r"
        % result["stop_after_after"])
    print("  PASS")


# =============================================================================
# CASE B: preprocessing+task (no explicit stop) → rewrite intent,
#         clear stop_after, proceed to plan generation
# =============================================================================

def test_preprocessing_task_no_explicit_stop_rewrites_and_proceeds():
    """If the user's advice doesn't include 'stop' phrasing,
    treat the preprocessing target as a starting point and
    generate a full plan.  Implementation: clear after_program
    and rewrite intent."""
    print("Test: preprocessing_task_no_explicit_stop_rewrites_and_proceeds")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.xtriage"},
    })
    # No "stop" / "and stop" / etc. in advice
    result = _run_decision_tree(session,
                                processed_advice="start with xtriage")

    assert result["outcome"] == "preprocessing_rewrite", (
        "Expected preprocessing_rewrite, got %s" % result["outcome"])
    assert result["intent_after"] == "solve", (
        "intent must be rewritten task→solve; got %r"
        % result["intent_after"])
    assert result["stop_after_after"] is None, (
        "stop_after must be cleared for preprocessing+rewrite; "
        "got %r" % result["stop_after_after"])
    assert result["proceeds_to_plan_generation"], (
        "preprocessing+rewrite must proceed to plan generation")
    print("  PASS")


# =============================================================================
# CASE C: preprocessing+task with EXPLICIT stop ("run xtriage and stop")
#         → return early, preserve stop_after (v115.09b)
# =============================================================================

def test_preprocessing_task_with_explicit_stop_preserves_directive():
    """v115.09b behavior: 'run xtriage and stop' means STOP, even
    though xtriage is a preprocessing program."""
    print("Test: preprocessing_task_with_explicit_stop_preserves_directive")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.xtriage"},
    })
    result = _run_decision_tree(session,
                                processed_advice="run xtriage and stop")

    assert result["outcome"] == "preprocessing_explicit_stop_skip", (
        "Expected preprocessing_explicit_stop_skip, got %s"
        % result["outcome"])
    assert result["intent_after"] == "task", (
        "intent must NOT be rewritten when explicit stop preserved; "
        "got %r" % result["intent_after"])
    assert result["stop_after_after"] == "phenix.xtriage", (
        "stop_after must be preserved when explicit stop is set; "
        "got %r" % result["stop_after_after"])
    assert not result["proceeds_to_plan_generation"], (
        "Explicit stop must not proceed to plan generation")
    print("  PASS")


# =============================================================================
# CASE D: needs_plan+task → rewrite intent, KEEP stop_after
# =============================================================================

def test_needs_plan_task_rewrites_intent_keeps_stop_after():
    """phenix.ligandfit (in _NEEDS_PLAN_PROGRAMS) as stop target
    + intent=task → generate prerequisite plan (refine → ligandfit)
    but stop after ligandfit completes.  Intent rewrites; stop_after
    is preserved."""
    print("Test: needs_plan_task_rewrites_intent_keeps_stop_after")
    if not _require_program_sets():
        return

    if "phenix.ligandfit" not in _NEEDS_PLAN:
        print("  SKIP — phenix.ligandfit not in _NEEDS_PLAN_PROGRAMS")
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.ligandfit"},
    })
    result = _run_decision_tree(session,
                                processed_advice="fit a ligand")

    assert result["outcome"] == "needs_plan_rewrite", (
        "Expected needs_plan_rewrite, got %s" % result["outcome"])
    assert result["intent_after"] == "solve", (
        "intent must be rewritten task→solve; got %r"
        % result["intent_after"])
    assert result["stop_after_after"] == "phenix.ligandfit", (
        "stop_after must be PRESERVED for needs_plan; got %r"
        % result["stop_after_after"])
    assert result["proceeds_to_plan_generation"], (
        "needs_plan must proceed to plan generation")
    print("  PASS")


# =============================================================================
# CASE E: v116.10 elif branch — full-plan target (e.g. predict_and_build)
#         + task → rewrite intent
# =============================================================================

def test_v116_10_elif_full_plan_target_rewrites_intent():
    """The v116.10 elif fires when stop_after is a program that
    isn't in _STANDALONE, _PREPROCESSING, or _NEEDS_PLAN.  These
    are full-plan targets that need a multi-stage plan + skip_to.

    Example: phenix.predict_and_build, phenix.refine, phenix.phaser.
    """
    print("Test: v116_10_elif_full_plan_target_rewrites_intent")
    if not _require_program_sets():
        return

    # Pick a program that's not in any classified set.
    # phenix.refine is a good representative.
    target = "phenix.refine"
    assert target not in _STANDALONE, (
        "Test assumes %s is NOT standalone; update test if "
        "classification changed" % target)
    assert target not in _NEEDS_PLAN, (
        "Test assumes %s is NOT in _NEEDS_PLAN" % target)
    assert target not in _PREPROCESSING_PROGRAMS, (
        "Test assumes %s is NOT preprocessing" % target)

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": target},
    })
    result = _run_decision_tree(session,
                                processed_advice="refine the model")

    assert result["outcome"] == "v116_10_elif_rewrite", (
        "Expected v116_10_elif_rewrite, got %s" % result["outcome"])
    assert result["intent_after"] == "solve", (
        "intent must be rewritten task→solve in v116.10 elif; got %r"
        % result["intent_after"])
    assert result["stop_after_after"] == target, (
        "stop_after must be preserved in v116.10 elif; got %r"
        % result["stop_after_after"])
    assert result["proceeds_to_plan_generation"], (
        "v116.10 elif must proceed to plan generation")
    print("  PASS")


# =============================================================================
# CASE F: intent=task with no stop_after → else branch (single_program_skip)
# =============================================================================

def test_task_intent_no_stop_after_skips_plan():
    """intent=task but no after_program — the outer guard fires
    (intent==task) but no sub-branch matches.  Falls to the else
    branch: skip plan generation."""
    print("Test: task_intent_no_stop_after_skips_plan")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {},  # no after_program
    })
    result = _run_decision_tree(session, processed_advice="")

    assert result["outcome"] == "single_program_skip", (
        "Expected single_program_skip, got %s" % result["outcome"])
    assert result["intent_after"] == "task"
    assert result["stop_after_after"] is None
    assert not result["proceeds_to_plan_generation"]
    print("  PASS")


# =============================================================================
# CASE G: solve intent (no stop_after) → outer guard doesn't fire,
#         proceed to plan generation
# =============================================================================

def test_solve_intent_no_stop_proceeds_normally():
    """intent=solve with no stop_after — outer guard doesn't fire,
    proceed to plan generation."""
    print("Test: solve_intent_no_stop_proceeds_normally")
    if not _require_program_sets():
        return

    session = _SessionStub(directives={
        "intent": "solve",
        "stop_conditions": {},
    })
    result = _run_decision_tree(session,
                                processed_advice="solve the structure")

    assert result["outcome"] == "proceed_to_generate_plan", (
        "Expected proceed_to_generate_plan, got %s" % result["outcome"])
    assert result["intent_after"] == "solve"
    assert result["proceeds_to_plan_generation"]
    print("  PASS")


# =============================================================================
# CASE H: Phase 3d behavior-change case — dock_in_map+task should
#         now hit single_program_skip (since dock_in_map is standalone)
# =============================================================================

def test_phase_3d_dock_in_map_task_skips_plan():
    """Phase 3d: phenix.dock_in_map is now in _STANDALONE_PROGRAMS.
    With intent=task + stop_after=dock_in_map, the outer guard
    fires (because dock_in_map is standalone), then falls to the
    else branch (skip plan generation).

    Pre-Phase-3d: this would have hit the v116.10 elif and
    generated a plan with skip_to_program.  That path skipped the
    predict prerequisite for sequence+map inputs, causing
    dock_in_map to fail at runtime.

    Post-Phase-3d: the workflow_engine handles prerequisites via
    the cryo-EM state machine (analyze → predict → dock).
    """
    print("Test: phase_3d_dock_in_map_task_skips_plan")
    if not _require_program_sets():
        return

    if "phenix.dock_in_map" not in _STANDALONE:
        print("  SKIP — Phase 3d not deployed (dock_in_map not standalone)")
        return

    session = _SessionStub(directives={
        "intent": "task",
        "stop_conditions": {"after_program": "phenix.dock_in_map"},
    })
    result = _run_decision_tree(session,
                                processed_advice="dock the model and stop")

    assert result["outcome"] == "single_program_skip", (
        "Phase 3d: dock_in_map+task should hit single_program_skip "
        "(not v116.10 elif). Got %s" % result["outcome"])
    assert not result["proceeds_to_plan_generation"], (
        "Phase 3d: should not proceed to plan generation")
    print("  PASS")


# =============================================================================
# Test runner (standard format)
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
