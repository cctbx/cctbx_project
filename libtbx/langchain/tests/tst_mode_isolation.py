"""Only `standard` mode may reach the summariser.

WHAT DEFECT THIS EXISTS FOR

Not a defect that has happened -- a trap that is one edit away.

`run_ai_analysis.run()` summarises whenever `something_to_analyze` is
true, and `determine_input_source`, which sets that flag, **has no
knowledge of `analysis_mode`**: it returns True for any supplied log.
So `run()` will summarise anything that reaches it.

The only thing preventing the four agent modes from summarising is a
dispatch decision made two files earlier: `run_job_locally` branches on
`analysis_mode` and each agent mode returns to its own entry point
before the standard path is reached.

**A future maintainer who adds a mode and forgets that branch gets the
summariser by default** -- including defect 6, which produced 382,721
characters from a 156,951-character log, reproducibly, and 2,026,595
from 312,011.

This was found while checking whether deleting the summariser would
break the agent. It would not. But the safety is structural rather than
stated, and structural safety that nobody has written down is one
refactor from gone.

WHAT IT CHECKS

Source-level, against the Phenix tree copy in `reference_source/`.
It cannot run the dispatcher, so it asserts the two properties the
dispatcher depends on:

  * every agent mode has its own branch in `run_job_locally`, and that
    branch returns;
  * none of the four agent entry points in `run_ai_analysis` reaches
    `summarize_log`, `get_log_info` or `_custom_log_chunker`.

WHERE IT LOOKS, AND WHY THAT MATTERS

`PHENIX_SOURCE_DIR` if set -- the live tree, which is the only thing
worth guarding. Otherwise `reference_source/`, which is a **snapshot
copied 2026-08-16** and will report the same answer forever no matter
what the tree does.

The first version of this file used the snapshot silently. **That is
false assurance: ten green tests that cannot fail.** Worse than no test,
and the same class of error as a skip reported as a pass. The run now
says which source it read, and says loudly when it read the snapshot.

It also walks the **transitive** call closure. The first version checked
only each entry point's own body -- a helper two calls deep could have
summarised and passed. (Checked: none does. But the test should not
depend on that having been checked by hand once.)
"""
from __future__ import absolute_import, division, print_function

import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

# --- module resolution -------------------------------------------------
# In the Phenix tree these modules live at
# libtbx/langchain/analysis/; standalone they sit beside the tests.
# Guideline 3 wants the dual-path form, and without it all four new
# suites failed on import in the installed layout -- checked by
# installing into a simulated tree, not by reading.
_ANALYSIS = os.path.join(os.path.dirname(_HERE), "analysis")
if os.path.isdir(_ANALYSIS):
    sys.path.insert(0, _ANALYSIS)
sys.path.insert(0, os.path.dirname(_HERE))
SNAPSHOT = os.path.join(os.path.dirname(_HERE), "reference_source")
LIVE = os.environ.get("PHENIX_SOURCE_DIR")


def source_root():
    """The live tree if we were told where it is; the snapshot otherwise.

    Returns (root, is_live). Callers must report is_live -- a green run
    against the snapshot proves nothing about the shipped code.
    """
    if LIVE and os.path.isdir(os.path.join(LIVE, "phenix_ai")):
        return LIVE, True
    if LIVE:
        print("  NOTE: PHENIX_SOURCE_DIR=%r has no phenix_ai/; "
              "falling back to the snapshot" % LIVE)
    return SNAPSHOT, False

AGENT_MODES = (
    ("agent_session", "run_agent_session_analysis"),
    ("advice_preprocessing", "run_advice_preprocessing"),
    ("directive_extraction", "run_directive_extraction"),
    ("failure_diagnosis", "run_failure_diagnosis"),
)

SUMMARISER = ("summarize_log", "get_log_info", "_custom_log_chunker")


def _require(env_name, path, what):
    """House pattern: SKIP when the variable is unset, FAIL LOUDLY when
    it is set but wrong.

    Without the second half a typo in PHIL_CORPUS_DIR silently falls
    back and the suite reports PASSED with everything skipped -- which
    the harness renders as a green tick. A skip is not a pass, and a
    skip that looks like a pass is worse than a failure.
    """
    if os.environ.get(env_name) and not os.path.isdir(path):
        raise SystemExit(
            "%s is set to %r but %s is not there.\n"
            "Unset it to skip these tests, or point it at the right "
            "directory." % (env_name, os.environ[env_name], what))

PASS, FAIL, SKIPPED = [], [], []


def ok(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print("  %s: %s%s" % ("PASS" if condition else "FAIL", name,
                          "  -- " + detail if detail and not condition else ""))


def function_body(path, name):
    lines = open(path, errors="replace").read().split("\n")
    start = None
    for i, line in enumerate(lines):
        if line.startswith("def %s(" % name):
            start = i
            break
    if start is None:
        return None
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if lines[j].startswith("def "):
            end = j
            break
    return "\n".join(lines[start:end])


def reaches(path, entry):
    """Every function `entry` can reach within this module.

    Body-only inspection would miss a helper two calls deep. Nothing in
    the current tree does that, but the guard should not rely on that
    having been true on the day it was written.
    """
    import ast
    tree = ast.parse(open(path, errors="replace").read())
    funcs = {n.name: n for n in ast.walk(tree)
             if isinstance(n, ast.FunctionDef)}
    seen = set()

    def walk(name):
        if name in seen or name not in funcs:
            return
        seen.add(name)
        for node in ast.walk(funcs[name]):
            if isinstance(node, ast.Call):
                func = node.func
                called = (func.id if isinstance(func, ast.Name)
                          else getattr(func, "attr", None))
                if called:
                    walk(called)

    walk(entry)
    return seen


def t_every_agent_mode_branches_and_returns(ai_analysis):
    body = function_body(ai_analysis, "run_job_locally")
    if body is None:
        SKIPPED.append("run_job_locally not found")
        print("  SKIP: run_job_locally not found")
        return
    for mode, _entry in AGENT_MODES:
        pattern = (r"if analysis_mode == ['\"]%s['\"]:\s*\n\s*return " % mode)
        ok("run_job_locally branches and returns for %s" % mode,
           re.search(pattern, body) is not None)

    # the standard path must come after all four
    idx = [body.find("'%s'" % m) for m, _ in AGENT_MODES]
    std = body.find("# Standard log analysis mode")
    ok("the standard path is reached only after every agent branch",
       std == -1 or all(0 <= i < std for i in idx),
       "standard at %d, branches at %s" % (std, idx))


def t_no_agent_entry_point_summarises(run_ai_analysis):
    for mode, entry in AGENT_MODES:
        if function_body(run_ai_analysis, entry) is None:
            SKIPPED.append("%s not found" % entry)
            print("  SKIP: %s not found" % entry)
            continue
        closure = reaches(run_ai_analysis, entry)
        hits = sorted(closure & set(SUMMARISER))
        ok("%s does not reach the summariser (transitively)" % mode,
           not hits, "found %s" % hits)


def t_the_search_actually_works(run_ai_analysis):
    """The control, and it had to change when the patch landed.

    It used to assert that `standard` mode DOES reach the summariser --
    true before the patch and the whole point of the patch afterwards.
    **The test encoded the pre-patch state as an invariant**, so
    applying the change turned the control red for the right reason,
    which is the wrong signal.

    A control must prove the SEARCH works, not that a particular defect
    is still present. So: `run()` must reach something it certainly
    calls."""
    if function_body(run_ai_analysis, "run") is None:
        SKIPPED.append("run() not found")
        print("  SKIP: run() not found")
        return
    closure = reaches(run_ai_analysis, "run")
    # Must be a function DEFINED in this module: reaches() walks the
    # module's own call graph, so an imported name such as
    # determine_input_source is invisible to it and made a false
    # control on the first attempt.
    ok("the reachability search finds a function run() certainly calls",
       "analyze_summary_with_rag" in closure,
       "closure had %d entries: %s" % (len(closure), sorted(closure)[:6]))
    ok("the search does not return everything indiscriminately",
       len(closure) < 40, "%d entries" % len(closure))


_require("PHENIX_SOURCE_DIR",
         os.path.join(os.environ.get("PHENIX_SOURCE_DIR", ""), "phenix_ai")
         if os.environ.get("PHENIX_SOURCE_DIR") else SNAPSHOT,
         "phenix_ai/ under the tree it points at")


def t_no_blocking_summary_gate_remains(run_ai_analysis):
    """The single-call path produces no summary. Any reachable gate that
    RETURNS or RAISES on `not log_info.summary` therefore rejects every
    request.

    This exists because exactly that happened: the gate in run() was
    corrected and the one inside analyze_summary_with_rag was not, so
    the first real run died with "Sorry: No summary to analyze". One
    gate was changed, four sites needed it, and nothing checked.

    A gate that ASSIGNS to log_info.summary is fine -- that is how
    summary_file gets filled from the report -- so exit and assignment
    are distinguished rather than lumped together.
    """
    import ast
    src = open(run_ai_analysis, errors="replace").read()
    tree = ast.parse(src)
    funcs = {n.name: n for n in ast.walk(tree)
             if isinstance(n, ast.FunctionDef)}

    def walk(name, seen=None):
        seen = seen or set()
        if name in seen or name not in funcs:
            return seen
        seen.add(name)
        for node in ast.walk(funcs[name]):
            if isinstance(node, ast.Call):
                func = node.func
                nm = (func.id if isinstance(func, ast.Name)
                      else getattr(func, "attr", None))
                if nm:
                    walk(nm, seen)
        return seen

    reachable = walk("run") | {"analyze_summary_with_rag"}
    blocking = []
    for fname in sorted(reachable):
        if fname not in funcs:
            continue
        for node in ast.walk(funcs[fname]):
            if not isinstance(node, ast.If):
                continue
            seg = ast.get_source_segment(src, node) or ""
            head = seg.split("\n")[0].strip()
            if "log_info.summary" not in head:
                continue
            if "analysis_payload" in seg or not head.startswith("if not"):
                continue
            if any(isinstance(n, (ast.Return, ast.Raise))
                   for n in ast.walk(node)):
                blocking.append("%s: %s" % (fname, head))
    ok("no reachable gate blocks on a missing summary", not blocking,
       "%s" % blocking)


def t_no_database_gate_blocks_the_report_path(ai_analysis, run_ai):
    """Retrieval is gone; nothing may still demand its database.

    Observed in production with anthropic:
        "'standard' analysis mode needs a local RAG database for
         provider 'anthropic' and none was found"
    The run was refused for want of a database the code no longer
    reads. Removing retrieval from analyzer.py was one change; FIVE
    places gated on that database, and one of them sat three lines
    after another edit in the same file.
    """
    import ast
    src = open(ai_analysis, errors="replace").read()

    ok("'standard' is treated as an LLM-only mode",
       re.search(r"_LLM_ONLY_MODES\s*=\s*\{[^}]*'standard'", src,
                 re.S) is not None)

    # No raise may be guarded by a missing database on this path.
    tree = ast.parse(src)
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        seg = ast.get_source_segment(src, node) or ""
        head = seg.split("\n")[0]
        if not re.search(r"ai_db_dir|has_database|isdir\(", head):
            continue
        if "install_ai_tools" in seg and "raise" in seg:
            offenders.append(head.strip()[:70])
    ok("no missing-database check raises on the report path",
       not offenders, "%s" % offenders)

    # Behaviour, not source text: no `return` inside run() may be
    # guarded by a database error. Grepping for the strings I happened
    # to write is error class 3, which has already produced three false
    # failures in this package.
    run_tree = ast.parse(open(run_ai, errors="replace").read())
    run_src = open(run_ai, errors="replace").read()
    aborts = []
    for node in ast.walk(run_tree):
        if not (isinstance(node, ast.FunctionDef) and node.name == "run"):
            continue
        for inner in ast.walk(node):
            if not isinstance(inner, ast.If):
                continue
            seg = ast.get_source_segment(run_src, inner) or ""
            head = seg.split("\n")[0]
            if "db_error" not in head:
                continue
            if any(isinstance(x, ast.Return) for x in ast.walk(inner)):
                aborts.append(head.strip()[:60])
    ok("a missing database does not abort run()", not aborts,
       "%s" % aborts)


def run_all_tests():
    root, is_live = source_root()
    if is_live:
        print("  source: LIVE TREE at %s" % root)
    else:
        print("  source: SNAPSHOT at reference_source/ (copied 2026-08-16)")
        print("  *** This guards a frozen copy, NOT the shipped code.")
        print("  *** Set PHENIX_SOURCE_DIR to check the tree that ships.")
    ai = os.path.join(root, "programs", "ai_analysis.py")
    raa = os.path.join(root, "phenix_ai", "run_ai_analysis.py")
    if not (os.path.exists(ai) and os.path.exists(raa)):
        print("  SKIP: reference_source/ not present; nothing checked")
        print("\n0 passed, 0 failed, 1 skipped, 1 total")
        print("*** 1 SKIPPED, and a skip is not a pass.")
        return True

    t_every_agent_mode_branches_and_returns(ai)
    t_no_agent_entry_point_summarises(raa)
    t_the_search_actually_works(raa)
    t_no_blocking_summary_gate_remains(raa)
    t_no_database_gate_blocks_the_report_path(ai, raa)

    total = len(PASS) + len(FAIL) + len(SKIPPED)
    print("\n%d passed, %d failed, %d skipped, %d total"
          % (len(PASS), len(FAIL), len(SKIPPED), total))
    if SKIPPED:
        print("*** %d SKIPPED, and a skip is not a pass." % len(SKIPPED))
    if FAIL:
        print("*** FAILED: %s" % ", ".join(FAIL))
    return not FAIL


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
