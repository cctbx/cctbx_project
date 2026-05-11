"""Tests for v116.10 Phase 3a — consistency between
ai_agent._STANDALONE_PROGRAMS and directive_extractor._ACTION_TABLE.

Before Phase 3a, the set of "standalone" programs (programs that
plan generation can short-circuit on) was duplicated in two places
inside `_initialize_plan_inner`.  Drift between those two tuples,
or between the duplicated tuples and the resolver's view of which
programs are one-shot targets, was kept in sync by manual review.

Phase 3a extracts _STANDALONE_PROGRAMS to a module-level frozenset
and adds these tests to enforce three invariants:

  1. Every program in _STANDALONE_PROGRAMS is a known PHENIX program
     (catches typos).

  2. Every single-program target in _ACTION_TABLE falls into exactly
     one classification: _STANDALONE_PROGRAMS, _NEEDS_PLAN_PROGRAMS,
     or the documented "full-plan" set (programs that need a multi-
     stage plan template).  If someone adds a new action with a new
     target program, this test fires until the developer decides
     where it goes.

  3. The clearly one-shot actions in _ACTION_TABLE
     ("analyze", "validate", "sharpen") MUST resolve to programs in
     _STANDALONE_PROGRAMS.  Drift-catcher: if someone "fixes"
     _STANDALONE_PROGRAMS by removing one of these, this test fails
     loudly.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# --- Path setup + libtbx stubs ---------------------------------------------
# Both ai_agent.py (the client) and directive_extractor.py (the server)
# import from libtbx.  We stub the minimum needed for the constants to be
# importable.

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
for d in (os.path.join(_PROJECT_ROOT, "agent"),
          os.path.join(_PROJECT_ROOT, "phenix", "programs")):
    if d not in sys.path:
        sys.path.insert(0, d)


def _install_stubs():
    class _Stub:
        def __init__(self):
            pass

    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")
    _ensure("libtbx.langchain.knowledge")

    yl = _ensure("libtbx.langchain.knowledge.yaml_loader")
    yl.get_workflow_steps = lambda et: {}
    yl.get_workflow_targets = lambda et, m: None
    yl.get_metric_threshold = lambda et, m: None
    yl.get_program = lambda p: None
    # get_all_programs returns the set of known programs — needed by
    # the typo-catch test.  We provide a comprehensive list that covers
    # everything we expect _STANDALONE_PROGRAMS and _ACTION_TABLE to
    # reference.
    yl.get_all_programs = lambda: {
        "phenix.xtriage", "phenix.mtriage",
        "phenix.molprobity", "phenix.map_correlations",
        "phenix.model_vs_data", "phenix.map_sharpening",
        "phenix.map_symmetry", "phenix.dock_in_map",
        "phenix.predict_and_build", "phenix.phaser",
        "phenix.autosol", "phenix.autobuild",
        "phenix.autobuild_denmod", "phenix.resolve_cryo_em",
        "phenix.refine", "phenix.real_space_refine",
        "phenix.ligandfit", "phenix.polder",
        "phenix.process_predicted_model",
    }

    pr = _ensure("libtbx.langchain.agent.program_registry")
    pr.ProgramRegistry = _Stub


# --- Module loading helpers ------------------------------------------------

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
    # _PROJECT_ROOT = .../modules/cctbx_project/libtbx/langchain/
    # Want:          .../modules/phenix/phenix/programs/ai_agent.py
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
            # find_dist_path may return modules/phenix or modules/phenix/phenix
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


def _load_action_table():
    """Load _ACTION_TABLE from directive_extractor in a few ways."""
    try:
        from libtbx.langchain.agent.directive_extractor import (
            _ACTION_TABLE)
        return _ACTION_TABLE
    except ImportError:
        pass
    try:
        from agent.directive_extractor import _ACTION_TABLE
        return _ACTION_TABLE
    except ImportError:
        from directive_extractor import _ACTION_TABLE
        return _ACTION_TABLE


def _load_standalone_programs():
    """Load _STANDALONE_PROGRAMS from the ai_agent client module."""
    # ai_agent.py is heavy (9000+ lines), but the constant lives near
    # the top.  Importing the full module needs all its dependencies,
    # which we don't fully stub.  Instead, parse the file and extract
    # the constant definitions via regex.
    import re
    ai_agent_path = _find_ai_agent_py()

    with open(ai_agent_path) as f:
        source = f.read()

    # Extract the _STANDALONE_PROGRAMS and _NEEDS_PLAN_PROGRAMS
    # definitions using a tolerant regex.  Each is a frozenset(...)
    # spanning multiple lines until a matching parenthesis.
    def _extract_frozenset(name):
        # Match `name = frozenset({...})` with possibly newlines inside the {}
        pattern = re.compile(
            r"^" + re.escape(name) + r"\s*=\s*frozenset\(\{([^}]*?)\}\)",
            re.MULTILINE | re.DOTALL)
        m = pattern.search(source)
        if not m:
            raise RuntimeError(
                "Cannot find %s in %s" % (name, ai_agent_path))
        body = m.group(1)
        # Eval the body as a tuple/set literal of strings
        items = set()
        for line in body.split("\n"):
            # Strip comments
            line = line.split("#", 1)[0].strip()
            # Match string literals
            for s in re.findall(r'"([^"]+)"|\'([^\']+)\'', line):
                # Match groups: one of the two is non-empty
                items.add(s[0] or s[1])
        return frozenset(items)

    return _extract_frozenset("_STANDALONE_PROGRAMS"), \
           _extract_frozenset("_NEEDS_PLAN_PROGRAMS")


# Try to install stubs and load the constants
try:
    from libtbx.langchain.agent.directive_extractor import _ACTION_TABLE
    _PHENIX_AVAILABLE = True
except ImportError:
    _install_stubs()
    _PHENIX_AVAILABLE = False


# These programs MUST need a multi-stage plan template (not standalone
# and not in _NEEDS_PLAN_PROGRAMS).  Programs in this set hit the
# v116.10 elif in _initialize_plan_inner — a plan is generated and
# skip_to_program fast-forwards to the target stage.
_FULL_PLAN_TARGETS = frozenset({
    "phenix.refine",
    "phenix.real_space_refine",
    "phenix.phaser",
    "phenix.autosol",
    "phenix.predict_and_build",
    "phenix.autobuild",
    "phenix.autobuild_denmod",
    "phenix.resolve_cryo_em",
    # v116.10 Phase 3d: phenix.map_symmetry and phenix.dock_in_map
    # were previously here (full-plan targets).  Phase 3d
    # reclassified them as standalone in ai_agent.py because the
    # workflow_engine path handles them more cleanly than the
    # plan + skip_to_program path (which SKIPPED prerequisites for
    # sequence-only inputs).  They no longer belong here.
})


# =====================================================================
# SECTION A: _STANDALONE_PROGRAMS is well-formed
# =====================================================================

def test_standalone_programs_are_known():
    """Every entry in _STANDALONE_PROGRAMS is a known PHENIX program.

    Catches typos like 'phenix.xtraige'.
    """
    print("Test: standalone_programs_are_known")
    try:
        from libtbx.langchain.knowledge.yaml_loader import get_all_programs
    except ImportError:
        from knowledge.yaml_loader import get_all_programs

    standalone, _needs_plan = _load_standalone_programs()
    known = set(get_all_programs())

    unknown = standalone - known
    assert not unknown, (
        "_STANDALONE_PROGRAMS contains unknown programs: %s\n"
        "(typo? deprecated? missing from YAML?)"
        % sorted(unknown))
    print("  PASS")


def test_needs_plan_programs_are_known():
    """Every entry in _NEEDS_PLAN_PROGRAMS is a known PHENIX program."""
    print("Test: needs_plan_programs_are_known")
    try:
        from libtbx.langchain.knowledge.yaml_loader import get_all_programs
    except ImportError:
        from knowledge.yaml_loader import get_all_programs

    _standalone, needs_plan = _load_standalone_programs()
    known = set(get_all_programs())

    unknown = needs_plan - known
    assert not unknown, (
        "_NEEDS_PLAN_PROGRAMS contains unknown programs: %s"
        % sorted(unknown))
    print("  PASS")


def test_standalone_and_needs_plan_are_disjoint():
    """A program can't be both standalone and a needs-plan template."""
    print("Test: standalone_and_needs_plan_are_disjoint")
    standalone, needs_plan = _load_standalone_programs()
    overlap = standalone & needs_plan
    assert not overlap, (
        "These programs are in BOTH _STANDALONE_PROGRAMS and "
        "_NEEDS_PLAN_PROGRAMS: %s. A program must be in exactly one."
        % sorted(overlap))
    print("  PASS")


# =====================================================================
# SECTION B: _ACTION_TABLE targets are classified
# =====================================================================

def test_action_table_targets_are_classified():
    """Every single-program target in _ACTION_TABLE must be in one of
    _STANDALONE_PROGRAMS, _NEEDS_PLAN_PROGRAMS, or _FULL_PLAN_TARGETS.

    Surfaces drift: if someone adds a new action with a new target,
    the test fires until they explicitly classify it.
    """
    print("Test: action_table_targets_are_classified")
    action_table = _load_action_table()
    standalone, needs_plan = _load_standalone_programs()

    targets = set()
    for action, info in action_table.items():
        for exp in ("xray", "cryoem"):
            prog = info.get(exp)
            if prog:
                targets.add(prog)

    classified = standalone | needs_plan | _FULL_PLAN_TARGETS
    unclassified = targets - classified
    assert not unclassified, (
        "These _ACTION_TABLE targets are not classified as standalone, "
        "needs-plan, or full-plan: %s.\n"
        "Decide which bucket each belongs to:\n"
        "  - _STANDALONE_PROGRAMS in ai_agent.py "
        "(one-shot, user-uploaded inputs)\n"
        "  - _NEEDS_PLAN_PROGRAMS in ai_agent.py "
        "(needs prereq, has plan template)\n"
        "  - _FULL_PLAN_TARGETS here in this test "
        "(needs full multi-stage plan)"
        % sorted(unclassified))
    print("  PASS")


# =====================================================================
# SECTION C: Drift-catchers — clearly one-shot actions stay standalone
# =====================================================================

def test_one_shot_actions_must_be_standalone():
    """Programs that _ACTION_TABLE treats as clear one-shot tasks
    ('analyze', 'validate', 'sharpen') MUST be in _STANDALONE_PROGRAMS.

    Drift-catcher: if someone removes molprobity from
    _STANDALONE_PROGRAMS while leaving 'validate' in _ACTION_TABLE,
    this test fails loudly.
    """
    print("Test: one_shot_actions_must_be_standalone")
    action_table = _load_action_table()
    standalone, _needs_plan = _load_standalone_programs()

    required = set()
    for action_name in ("analyze", "validate", "sharpen"):
        info = action_table.get(action_name, {})
        for exp in ("xray", "cryoem"):
            prog = info.get(exp)
            if prog:
                required.add(prog)

    missing = required - standalone
    assert not missing, (
        "These _ACTION_TABLE one-shot targets MUST be in "
        "_STANDALONE_PROGRAMS:\n  %s\n"
        "These actions ('analyze', 'validate', 'sharpen') are "
        "semantically standalone — the user gives the inputs, the "
        "program runs once, the workflow stops."
        % sorted(missing))
    print("  PASS")


# =====================================================================
# SECTION D: Phase 3b — Decision-tree behavior for the v116.10 elif
# =====================================================================
#
# The v116.10 elif in _initialize_plan_inner generates a plan when
# stop_after is a program that:
#   - is truthy
#   - is NOT in _preprocessing_programs
#   - is NOT in _NEEDS_PLAN_PROGRAMS
#   - is NOT in _STANDALONE_PROGRAMS
#
# Testing _initialize_plan_inner directly requires heavy scaffolding
# (self.vlog, self.params, session, etc.).  These tests instead exercise
# the SET-MEMBERSHIP LOGIC the elif depends on: every full-plan target
# (predict_and_build, refine, phaser, etc.) falls through the chain
# correctly to the elif rather than being captured by an earlier branch.

# The preprocessing programs match the inline set in _initialize_plan_inner
_PREPROCESSING_PROGRAMS = frozenset({
    "phenix.process_predicted_model",
    "phenix.xtriage",
    "phenix.mtriage",
})


def test_full_plan_targets_hit_the_v116_10_elif():
    """Programs in _FULL_PLAN_TARGETS must NOT be in _preprocessing,
    _STANDALONE, or _NEEDS_PLAN — so they reach the v116.10 elif.

    If predict_and_build, refine, phaser, etc. accidentally end up
    in _STANDALONE_PROGRAMS, the v116.10 elif won't fire, and
    "predict and stop" or "refine and stop" sessions will skip
    plan generation and AUTO-STOP at workflow init.
    """
    print("Test: full_plan_targets_hit_the_v116_10_elif")
    standalone, needs_plan = _load_standalone_programs()

    misclassified = []
    for prog in _FULL_PLAN_TARGETS:
        if prog in _PREPROCESSING_PROGRAMS:
            misclassified.append(
                "%s: in _PREPROCESSING_PROGRAMS (would clear "
                "after_program — wrong for full-plan targets)" % prog)
        if prog in standalone:
            misclassified.append(
                "%s: in _STANDALONE_PROGRAMS (would skip plan "
                "generation — wrong; full-plan targets need a plan)"
                % prog)
        if prog in needs_plan:
            misclassified.append(
                "%s: in _NEEDS_PLAN_PROGRAMS (would use the specific "
                "polder/ligandfit branch — wrong template for "
                "full-plan targets)" % prog)

    assert not misclassified, (
        "These full-plan targets are misclassified:\n  %s\n"
        "They must fall through to the v116.10 elif so a full "
        "multi-stage plan is generated and skip_to_program can "
        "fast-forward to the target."
        % "\n  ".join(misclassified))
    print("  PASS")


def test_standalone_programs_short_circuit_correctly():
    """Programs in _STANDALONE_PROGRAMS must NOT be in
    _preprocessing or _NEEDS_PLAN_PROGRAMS.

    The outer condition `_stop_after in _STANDALONE_PROGRAMS`
    triggers entering the block; the inner chain then checks
    preprocessing first, needs-plan second, and the v116.10 elif
    third.  A standalone program would mis-fire on the preprocessing
    or needs-plan branch if it also appeared there.
    """
    print("Test: standalone_programs_short_circuit_correctly")
    standalone, needs_plan = _load_standalone_programs()

    misclassified = []
    for prog in standalone:
        if prog in _PREPROCESSING_PROGRAMS and prog not in (
                "phenix.xtriage", "phenix.mtriage"):
            # xtriage and mtriage are deliberately in both
            # (standalone for "run xtriage and stop", preprocessing
            # for "use xtriage as analysis step").  Other overlap
            # would be an error.
            misclassified.append(
                "%s: in both _STANDALONE_PROGRAMS and "
                "_PREPROCESSING_PROGRAMS (only xtriage/mtriage "
                "may overlap)" % prog)
        if prog in needs_plan:
            misclassified.append(
                "%s: in both _STANDALONE_PROGRAMS and "
                "_NEEDS_PLAN_PROGRAMS (mutually exclusive)" % prog)

    assert not misclassified, (
        "Standalone-program classification errors:\n  %s"
        % "\n  ".join(misclassified))
    print("  PASS")


# =====================================================================
# SECTION E: Documentation invariants
# =====================================================================

def test_standalone_is_a_frozenset():
    """The constant is a frozenset (immutable + supports set ops)."""
    print("Test: standalone_is_a_frozenset")
    standalone, _needs_plan = _load_standalone_programs()
    # Our extraction returns frozenset; this confirms the source uses
    # frozenset rather than set or tuple.  Reading the source directly
    # to confirm:
    import re
    try:
        ai_agent_path = _find_ai_agent_py()
    except RuntimeError:
        print("  SKIP (ai_agent.py not found)")
        return
    with open(ai_agent_path) as f:
        source = f.read()
    assert re.search(
        r"^_STANDALONE_PROGRAMS\s*=\s*frozenset\(",
        source, re.MULTILINE), (
        "_STANDALONE_PROGRAMS should be a frozenset for safety "
        "(prevents accidental mutation)")
    print("  PASS")


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
