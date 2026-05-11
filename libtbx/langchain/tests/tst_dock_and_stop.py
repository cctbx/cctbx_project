"""Tests for v116.10 Phase 3d — behavior change for
phenix.map_symmetry and phenix.dock_in_map.

Phase 3a was a pure refactor: it extracted _STANDALONE_PROGRAMS
as a module-level constant without changing membership.

Phase 3d is an EXPLICIT BEHAVIOR CHANGE: it adds map_symmetry and
dock_in_map to _STANDALONE_PROGRAMS.  Pre-Phase-3d these were
"full-plan targets" (caught by the v116.10 elif), which generated
a plan and used skip_to_program to fast-forward.  Post-Phase-3d
they skip plan generation and route through workflow_engine
instead.

This file captures the behavior change in three ways:

1. test_phase_3d_membership: documents that the two programs are
   now in _STANDALONE_PROGRAMS (would fail on pre-Phase-3d code).

2. test_phase_3d_decision_tree_changes: traces the decision tree
   for (program × intent) combinations and asserts the expected
   pre/post differences.

3. test_phase_3d_does_not_affect_other_programs: regression
   guard — Phase 3d's behavior change is LIMITED to map_symmetry
   and dock_in_map.  Every other program's decision-tree
   behavior is unchanged.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))


# --- Reference: pre-Phase-3d _STANDALONE_PROGRAMS (6 programs) ---
_PRE_PHASE_3D_STANDALONE = frozenset({
    "phenix.xtriage",
    "phenix.mtriage",
    "phenix.molprobity",
    "phenix.model_vs_data",
    "phenix.map_correlations",
    "phenix.map_sharpening",
})


# --- Reference: programs Phase 3d reclassifies ---
_PHASE_3D_RECLASSIFIED = frozenset({
    "phenix.map_symmetry",
    "phenix.dock_in_map",
})


# --- Other inline constants from _initialize_plan_inner ---
_PREPROCESSING_PROGRAMS = frozenset({
    "phenix.process_predicted_model",
    "phenix.xtriage",
    "phenix.mtriage",
})


def _find_ai_agent_py():
    """Find phenix/programs/ai_agent.py across common PHENIX layouts.

    Tries (in order):
      1. Standalone development layout: tests/ next to phenix/programs/
      2. PHENIX install layout: cctbx_project/libtbx/langchain/tests/ at
         siblings of phenix/phenix/programs/
      3. libtbx.env_config.find_dist_path('phenix')
      4. $PHENIX environment variable
    """
    tried = []

    # 1. Standalone development layout
    cand = os.path.join(_PROJECT_ROOT, "phenix", "programs", "ai_agent.py")
    tried.append(cand)
    if os.path.exists(cand):
        return cand

    # 2. PHENIX install layout: walk up to <modules>/, over to phenix/phenix/
    cand = os.path.normpath(os.path.join(
        _PROJECT_ROOT, "..", "..", "..",
        "phenix", "phenix", "programs", "ai_agent.py"))
    tried.append(cand)
    if os.path.exists(cand):
        return cand

    # 3. libtbx environment (when running under phenix.python)
    try:
        from libtbx.env_config import unpickle
        env = unpickle()
        phenix_dist = env.find_dist_path('phenix', default=None)
        if phenix_dist:
            for sub in ("phenix/programs/ai_agent.py",
                        "programs/ai_agent.py"):
                cand = os.path.join(phenix_dist, sub)
                tried.append(cand)
                if os.path.exists(cand):
                    return cand
    except Exception as e:
        tried.append("libtbx.env_config: %s" % e)

    # 4. $PHENIX environment variable
    phenix_env = os.environ.get("PHENIX")
    if phenix_env:
        cand = os.path.join(
            phenix_env, "modules", "phenix", "phenix",
            "programs", "ai_agent.py")
        tried.append(cand)
        if os.path.exists(cand):
            return cand

    raise RuntimeError(
        "Cannot locate ai_agent.py. Tried:\n  %s"
        % "\n  ".join(tried))


def _load_standalone():
    """Extract _STANDALONE_PROGRAMS from ai_agent.py source."""
    ai_agent_path = _find_ai_agent_py()
    with open(ai_agent_path) as f:
        source = f.read()
    m = re.search(
        r"^_STANDALONE_PROGRAMS\s*=\s*frozenset\(\{([^}]*?)\}\)",
        source, re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError("Cannot find _STANDALONE_PROGRAMS")
    items = set()
    for line in m.group(1).split("\n"):
        line = line.split("#", 1)[0].strip()
        for s in re.findall(r'"([^"]+)"|\'([^\']+)\'', line):
            items.add(s[0] or s[1])
    return frozenset(items)


def _load_needs_plan():
    """Extract _NEEDS_PLAN_PROGRAMS from ai_agent.py source."""
    ai_agent_path = _find_ai_agent_py()
    with open(ai_agent_path) as f:
        source = f.read()
    m = re.search(
        r"^_NEEDS_PLAN_PROGRAMS\s*=\s*frozenset\(\{([^}]*?)\}\)",
        source, re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError("Cannot find _NEEDS_PLAN_PROGRAMS")
    items = set()
    for line in m.group(1).split("\n"):
        line = line.split("#", 1)[0].strip()
        for s in re.findall(r'"([^"]+)"|\'([^\']+)\'', line):
            items.add(s[0] or s[1])
    return frozenset(items)


def _trace_decision_tree(intent, stop_after, standalone_set,
                         needs_plan_set=None):
    """Simulate the decision tree in _initialize_plan_inner.

    Returns one of: 'preprocessing', 'needs_plan', 'v116_10_elif',
    'else_return', 'outer_skipped'.
    """
    if needs_plan_set is None:
        needs_plan_set = frozenset(
            {"phenix.polder", "phenix.ligandfit"})
    if intent == "task" or (stop_after and
                             stop_after in standalone_set):
        if stop_after in _PREPROCESSING_PROGRAMS:
            return "preprocessing"
        elif stop_after in needs_plan_set:
            return "needs_plan"
        elif stop_after and stop_after not in standalone_set:
            return "v116_10_elif"
        else:
            return "else_return"
    else:
        return "outer_skipped"


# =====================================================================
# Test 1: Membership change
# =====================================================================

def test_phase_3d_reclassified_programs_are_standalone():
    """map_symmetry and dock_in_map are in _STANDALONE_PROGRAMS
    after Phase 3d.

    This test FAILS on pre-Phase-3d code (where they were full-plan
    targets).  Documents the intentional behavior change.
    """
    print("Test: phase_3d_reclassified_programs_are_standalone")
    standalone = _load_standalone()
    missing = _PHASE_3D_RECLASSIFIED - standalone
    assert not missing, (
        "Phase 3d should reclassify these as standalone but "
        "they're missing from _STANDALONE_PROGRAMS: %s\n"
        "(If you reverted Phase 3d, also remove this test.)"
        % sorted(missing))
    print("  PASS")


def test_phase_3d_preserves_original_six_standalone_programs():
    """The 6 pre-Phase-3a/3d standalone programs are still
    standalone.

    Regression guard: Phase 3d should ADD two programs, never
    remove any.
    """
    print("Test: phase_3d_preserves_original_six_standalone_programs")
    standalone = _load_standalone()
    missing = _PRE_PHASE_3D_STANDALONE - standalone
    assert not missing, (
        "Phase 3d should preserve the original 6 standalone "
        "programs but these are missing: %s" % sorted(missing))
    print("  PASS")


# =====================================================================
# Test 2: Decision-tree changes
# =====================================================================

def test_phase_3d_decision_tree_changes_for_reclassified_programs():
    """For each Phase 3d reclassified program and each intent,
    document the expected pre→post decision-tree change.

    Pre-Phase-3d behavior:
      - intent=task   → v116_10_elif  (plan generated)
      - intent=solve  → outer_skipped (plan generated normally)

    Post-Phase-3d behavior:
      - intent=task   → else_return   (plan skipped)
      - intent=solve  → else_return   (plan skipped)
    """
    print(
      "Test: phase_3d_decision_tree_changes_for_reclassified_programs")
    standalone = _load_standalone()

    expected = {
        ("phenix.map_symmetry", "task"):  ("v116_10_elif", "else_return"),
        ("phenix.map_symmetry", "solve"): ("outer_skipped", "else_return"),
        ("phenix.dock_in_map",  "task"):  ("v116_10_elif", "else_return"),
        ("phenix.dock_in_map",  "solve"): ("outer_skipped", "else_return"),
    }

    failures = []
    for (prog, intent), (pre_expected, post_expected) in expected.items():
        pre_actual = _trace_decision_tree(
            intent, prog, _PRE_PHASE_3D_STANDALONE)
        post_actual = _trace_decision_tree(
            intent, prog, standalone)
        if pre_actual != pre_expected:
            failures.append(
                "%s/intent=%s: pre-Phase-3d expected %r, got %r"
                % (prog, intent, pre_expected, pre_actual))
        if post_actual != post_expected:
            failures.append(
                "%s/intent=%s: post-Phase-3d expected %r, got %r"
                % (prog, intent, post_expected, post_actual))

    assert not failures, (
        "Phase 3d decision-tree behavior mismatch:\n  %s"
        % "\n  ".join(failures))
    print("  PASS — 4 (program × intent) traces match expected pre/post")


# =====================================================================
# Test 3: Phase 3d does NOT affect any other program
# =====================================================================

def test_phase_3d_does_not_affect_other_programs():
    """For every program that's NOT in _PHASE_3D_RECLASSIFIED, the
    decision-tree behavior must be the same pre/post Phase 3d.

    Regression guard: Phase 3d's blast radius is limited to
    map_symmetry and dock_in_map.  Every other program is
    unchanged.
    """
    print("Test: phase_3d_does_not_affect_other_programs")
    standalone_post = _load_standalone()

    # Comprehensive list of programs that appear in any
    # classification or commonly in directives
    all_programs = [
        # Original 6 standalone
        "phenix.xtriage", "phenix.mtriage", "phenix.molprobity",
        "phenix.model_vs_data", "phenix.map_correlations",
        "phenix.map_sharpening",
        # _NEEDS_PLAN_PROGRAMS
        "phenix.polder", "phenix.ligandfit",
        # Full-plan targets
        "phenix.predict_and_build", "phenix.phaser",
        "phenix.refine", "phenix.real_space_refine",
        "phenix.autosol", "phenix.autobuild",
        "phenix.autobuild_denmod", "phenix.resolve_cryo_em",
        # Preprocessing
        "phenix.process_predicted_model",
        # None / empty stop_after
        None,
    ]

    failures = []
    for prog in all_programs:
        if prog in _PHASE_3D_RECLASSIFIED:
            continue  # Skip the two programs we expect to change
        for intent in ["task", "solve", None]:
            pre = _trace_decision_tree(
                intent, prog, _PRE_PHASE_3D_STANDALONE)
            post = _trace_decision_tree(
                intent, prog, standalone_post)
            if pre != post:
                failures.append(
                    "%s/intent=%s: pre=%r, post=%r"
                    % (prog, intent, pre, post))

    assert not failures, (
        "Phase 3d should NOT affect these programs but did:\n  %s"
        % "\n  ".join(failures))
    n_programs = len(all_programs) - len(_PHASE_3D_RECLASSIFIED)
    print("  PASS — %d programs × 3 intents = %d traces unchanged"
          % (n_programs, n_programs * 3))


# =====================================================================
# Test 4: Sequence + map case for dock_in_map
# =====================================================================

def test_phase_3d_motivating_case_documented():
    """The motivating case for Phase 3d: 'dock and stop' with
    sequence + map (no model).

    Pre-Phase-3d:
      1. Plan generated: [analyze: mtriage, obtain_model:
         predict_and_build, obtain_model: dock_in_map] (or similar
         structure depending on plan generator)
      2. skip_to_program(dock_in_map) marks predict_and_build SKIPPED
      3. dock_in_map runs without a model → FAILS at runtime

    Post-Phase-3d:
      1. Plan generation SKIPPED (dock_in_map now in standalone)
      2. workflow_engine handles the session:
         - state=analyze, valid_programs=[mtriage] → LLM picks mtriage
         - state=obtain_model, valid_programs has predict_and_build
           (has sequence + map) → LLM picks predict_and_build
         - state=obtain_model, after predict_and_build runs, model
           now exists → dock_in_map valid → LLM picks dock_in_map
         - Stop after dock_in_map

    This test documents the difference; the actual workflow_engine
    behavior is verified by integration tests against the tutorial
    corpus, not here.
    """
    print("Test: phase_3d_motivating_case_documented")
    # This test is documentation-only.  It verifies that the
    # classification supports the motivating case.
    standalone = _load_standalone()
    assert "phenix.dock_in_map" in standalone, (
        "Phase 3d motivating case requires dock_in_map in standalone")
    print("  PASS — classification supports the motivating case "
          "(integration test in tutorial corpus)")


# =====================================================================
# Runner
# =====================================================================

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
