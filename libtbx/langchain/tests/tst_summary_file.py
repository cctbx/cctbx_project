"""Tests for the deterministic summary -- option (c).

WHAT THIS EXISTS FOR

After the single-call patch, `summary_file` held a byte-identical copy
of the report: 40 of 40 outputs across two passes. That was my choice --
the cheapest thing that kept `display_summary_html`, the history record
and `load_existing_analysis` working -- and it is semantically wrong. A
summary that is the analysis is not a summary.

The fix fills `summary_file` from the LOG, deterministically, with no
LLM involved.

WHY NOT SPLIT THE REPORT INSTEAD

Measured across the twenty reports: 14 open their first section with
`### 1.` and 6 with `**1.`. That is the same two-format variation that
broke the `**3` scrape and cost this project a fortnight. Building a
second consumer of LLM output formatting would repeat defect 5
knowingly.

WHAT THE FIRST DRAFT OF THIS WORK GOT WRONG, CHECKED BEFORE CODING

* `orientation_strip.py` carried its own `load_registry` and `resolve`,
  duplicating `program_identity`. They agreed on 10 of 10 probes -- by
  luck, not construction, and those two have already diverged once in
  this project at 19/20 against a measured 20/20.
* Its paths did not resolve once installed: it looked for
  `HERE/vendor/` and `HERE/knowledge/programs.yaml`, neither of which
  exists at `langchain/analysis/`. It would have silently fallen back to
  a ten-name registry and failed to import its extractor.
* Its closing line promised `[the full analysis is being prepared]`,
  which is wrong in a file written after the analysis exists.

Each has a test below.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ANALYSIS = os.path.join(os.path.dirname(_HERE), "analysis")
if os.path.isdir(_ANALYSIS):
    sys.path.insert(0, _ANALYSIS)
sys.path.insert(0, os.path.dirname(_HERE))


def _require(env_name, path, what):
    if os.environ.get(env_name) and not os.path.isdir(path):
        raise SystemExit(
            "%s is set to %r but %s is not there.\nUnset it to skip these "
            "tests, or point it at the right directory."
            % (env_name, os.environ[env_name], what))


def _corpus_root():
    env = os.environ.get("PHIL_CORPUS_DIR")
    for cand in ([env, os.path.dirname((env or "").rstrip(os.sep))] +
                 [os.path.join(os.path.dirname(_HERE), "test_corpus")]):
        if cand and os.path.isdir(os.path.join(cand, "corpus_gui")):
            return cand
    return os.path.join(os.path.dirname(_HERE), "test_corpus")


CORPUS = os.path.join(_corpus_root(), "corpus_gui")
_require("PHIL_CORPUS_DIR", CORPUS, "corpus_gui")

PASS, FAIL, SKIPPED = [], [], []


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  %s: %s%s" % ("PASS" if cond else "FAIL", name,
                          "  -- " + detail if detail and not cond else ""))


def logs():
    return sorted(f for f in os.listdir(CORPUS) if f.endswith(".log"))


# ------------------------------------------------------------------ tests

def t_every_log_yields_a_summary(mod):
    """The old summariser produced nothing useful on some logs and 2 MB
    on others. This one must produce something bounded on all twenty."""
    empty, oversize = [], []
    for name in logs():
        text = open(os.path.join(CORPUS, name), errors="replace").read()
        out = mod.summary_for_log(text, name)
        if not out.strip():
            empty.append(name)
        if len(out) > 4000:
            oversize.append((name, len(out)))
    ok("every corpus log yields a non-empty summary", not empty,
       "%s" % empty[:3])
    ok("no summary exceeds 4 KB", not oversize, "%s" % oversize[:3])


def t_summary_is_smaller_than_the_log(mod):
    """It is a summary. The old one returned 382,721 characters from a
    156,951-character log, reproducibly."""
    bad = []
    for name in logs():
        text = open(os.path.join(CORPUS, name), errors="replace").read()
        out = mod.summary_for_log(text, name)
        if len(out) >= len(text):
            bad.append((name, len(out), len(text)))
    ok("every summary is smaller than its log", not bad, "%s" % bad[:2])


def t_no_llm_derived_value_reaches_it(mod):
    """Asserted on the AST with docstrings stripped, as in
    tst_program_identity -- an earlier guard of this kind failed on
    prose that named the very thing it forbids."""
    import ast
    tree = ast.parse(open(mod.__file__.replace(".pyc", ".py"),
                          errors="replace").read())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef)):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body.pop(0)
    code = ast.dump(tree)
    forbidden = ("summarize_log", "get_log_info", "invoke",
                 "processed_log_dict", "analysis_payload")
    hits = [w for w in forbidden if w in code]
    ok("no LLM-derived value reaches the summary", not hits,
       "found %s" % hits)


def t_it_does_not_duplicate_the_resolver(mod):
    """It defined its own load_registry and resolve, duplicating
    program_identity -- the error class this whole change cites as its
    reason for reusing tested code."""
    import ast
    src = open(mod.__file__.replace(".pyc", ".py"), errors="replace").read()
    tree = ast.parse(src)
    defined = {n.name for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef)}
    dupes = defined & {"load_registry", "resolve", "_canonical", "_strict"}
    ok("it does not define its own program resolver", not dupes,
       "defines %s" % sorted(dupes))
    ok("it imports program_identity instead",
       "program_identity" in src)


def t_it_agrees_with_the_resolver(mod):
    """The corpus invariant. Whatever it names must be what
    program_identity names -- they diverged once already, 19/20 against
    a measured 20/20."""
    import program_identity as ident
    reg = os.path.join(os.path.dirname(_HERE), "knowledge",
                       "programs.yaml")
    bad = []
    for name in logs():
        text = open(os.path.join(CORPUS, name), errors="replace").read()
        want = ident.resolve_program(
            log_text=text, log_path=name,
            registry_path=reg if os.path.exists(reg) else None).name
        out = mod.summary_for_log(text, name)
        if want and want not in out:
            bad.append((name, want))
    ok("the summary names the program the resolver resolves", not bad,
       "%s" % bad[:3])


def t_it_promises_nothing_pending(mod):
    """The on-screen form ends '[the full analysis is being prepared]'.
    The saved file is written AFTER the analysis exists."""
    text = open(os.path.join(CORPUS, logs()[0]), errors="replace").read()
    saved = mod.summary_for_log(text, logs()[0])
    ok("the saved summary does not promise a pending analysis",
       "being prepared" not in saved and "follows in" not in saved,
       "%r" % saved[-90:])


def t_it_resolves_its_paths_from_the_package(mod):
    """It looked for HERE/vendor and HERE/knowledge, neither of which
    exists at langchain/analysis/. It would have shipped falling back to
    a ten-name registry with no extractor."""
    # Check BEHAVIOUR, not source text. The first version of these two
    # assertions grepped for strings I had written myself and failed a
    # correct implementation twice -- error class 3, in the test written
    # to catch a path bug.
    ok("the registry resolves into the package, not beside the module",
       os.path.basename(os.path.dirname(os.path.dirname(mod.KNOWLEDGE)))
       != "analysis",
       "KNOWLEDGE=%s" % mod.KNOWLEDGE)
    ok("the resolved registry file exists", os.path.exists(mod.KNOWLEDGE),
       "KNOWLEDGE=%s" % mod.KNOWLEDGE)


def t_the_extractor_is_imported_not_vendored(mod):
    """The tree already ships log_structure_extractor at
    libtbx/langchain/ -- run_all_tests.py registers its suite and names
    that path. Vendoring a second copy would put two modules of one name
    on sys.path, which that file's own comment warns silently shadow
    each other."""
    # Behaviour again: whatever it bound must be a usable extractor,
    # and it must not have come from a copy sitting inside this package.
    ok("the bound extractor exposes scan()", hasattr(mod.X, "scan"))
    here = os.path.dirname(os.path.abspath(mod.__file__))
    ok("no second copy of the extractor ships beside the module",
       not os.path.exists(os.path.join(here,
                                       "log_structure_extractor.py")))


def t_no_python_repr_leaks_into_the_summary(mod):
    """A ControlSkip carries `what` and `reason` and has NO `text`, so a
    getattr fallback chain ended at the object and printed
    "ControlSkip(line=13829, ...)" into a user-facing file. Present on 4
    of 20 corpus logs and on a 649 KB one, and shipped because I never
    read the output."""
    import re
    bad = []
    for name in logs():
        text = open(os.path.join(CORPUS, name), errors="replace").read()
        out = mod.summary_for_log(text, name)
        if re.search(r"\w+\(line=\d+", out) or " object at 0x" in out:
            bad.append(name)
    ok("no Python repr appears in any summary", not bad, "%s" % bad[:3])


def t_nothing_is_claimed_as_written(mod):
    """The slot was labelled "wrote:" and printed completion records --
    "Job complete", "Finished: Fri Jun 19 15:22:26 2026". On 13 of 20
    logs it told the reader a timestamp was a file. The extractor
    exposes no output-file structure, so nothing may claim one."""
    bad = []
    for name in logs():
        text = open(os.path.join(CORPUS, name), errors="replace").read()
        out = mod.summary_for_log(text, name)
        if "wrote:" in out:
            bad.append(name)
    ok("no summary claims a file was written", not bad, "%s" % bad[:3])


def t_large_log_is_handled(mod):
    """The corpus tops out at 305 KB. A real 649 KB AutoBuild log found
    both defects above on its first run, which is why this exists: the
    summary must stay bounded and fast as the log grows."""
    big = os.environ.get("PHIL_BIG_LOG")
    if not big or not os.path.exists(big):
        SKIPPED.append("no large log supplied")
        print("  SKIP: set PHIL_BIG_LOG to a multi-hundred-KB log")
        return
    import time
    text = open(big, errors="replace").read()
    start = time.time()
    out = mod.summary_for_log(text, os.path.basename(big))
    elapsed = time.time() - start
    print("      %d chars in, %d chars out, %.2f s"
          % (len(text), len(out), elapsed))
    ok("a large log still yields a bounded summary",
       0 < len(out) < 4000, "%d chars" % len(out))
    ok("the summary is built in under 2 s", elapsed < 2.0,
       "%.2f s" % elapsed)


def run_all_tests():
    if not os.path.isdir(CORPUS):
        print("  SKIP: corpus not found. Set PHIL_CORPUS_DIR.")
        print("\n0 passed, 0 failed, 1 skipped, 1 total")
        print("*** 1 SKIPPED, and a skip is not a pass.")
        return True
    try:
        import orientation_strip as mod
    except ImportError as error:
        print("  FAIL: orientation_strip not importable -- %s" % error)
        print("\n0 passed, 1 failed, 0 skipped, 1 total")
        return False
    if not hasattr(mod, "summary_for_log"):
        print("  FAIL: summary_for_log() does not exist yet")
        print("\n0 passed, 1 failed, 0 skipped, 1 total")
        return False

    for fn in (t_every_log_yields_a_summary,
               t_summary_is_smaller_than_the_log,
               t_no_llm_derived_value_reaches_it,
               t_it_does_not_duplicate_the_resolver,
               t_it_agrees_with_the_resolver,
               t_it_promises_nothing_pending,
               t_it_resolves_its_paths_from_the_package,
               t_the_extractor_is_imported_not_vendored,
               t_no_python_repr_leaks_into_the_summary,
               t_nothing_is_claimed_as_written,
               t_large_log_is_handled):
        fn(mod)

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
