"""Tests for program identity resolution -- Step 1 of IMPLEMENTATION.md.

WHAT DEFECT THIS EXISTS FOR

Two defects, both measured on real runs:

  defect 3  On the default server path the program name is never derived
            from the file name, so content patterns decide it. Wrong on
            6 of 20 corpus logs. Seen live: an xtriage log identified as
            phenix.phaser via the substring 'molecular replacement'.

  defect 4  ai_analysis.py:1067 writes "Log file for phenix.<derived>"
            into the model's prompt. On 20 of 20 corpus logs the derived
            name is malformed -- phenix.xtriage_1, phenix.refine_3 --
            and names a program that does not exist.

WHAT REPLACES THEM

  caller  ->  banner  ->  file name  ->  refuse

Measured at 20 of 20, zero wrong, zero unknown on both corpus shapes.
The banner and the file name are independent signals and agree on all 9
logs where both fire.

THE GUARD THAT MATTERS

Not "run it twenty times and get one answer" -- that is vacuous once the
mechanism is deterministic, and the variance it was meant to catch came
from an LLM scrape that this design does not consult. The guard is that
NO LLM-DERIVED VALUE CAN REACH THE PROGRAM NAME, plus registry
membership, plus explicit refusal.

PROVENANCE IS NOT DECORATION. A file name is INFERRED: registry
membership proves a name exists, not that it is right, and a user can
rename a file. Callers that need certainty must be able to tell.
"""
from __future__ import absolute_import, division, print_function

import os
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

def _corpus_root():
    """Resolve the corpus root whether PHIL_CORPUS_DIR points at the
    parent (test_corpus/) or at one shape inside it (test_corpus/corpus/).

    The harness sets it to the shape directory. The first version of
    this file assumed the parent, so it looked for
    test_corpus/corpus/corpus and skipped BOTH corpus checks -- the
    corpus-level invariant, the most important test here, silently did
    not run. That is exactly the defect found in tst_api.py and
    tst_verifier.py, reproduced in the suite written to fix a different
    one. A skip is not a pass.
    """
    env = os.environ.get("PHIL_CORPUS_DIR")
    candidates = []
    if env:
        candidates += [env, os.path.dirname(env.rstrip(os.sep))]
    candidates.append(os.path.join(os.path.dirname(_HERE), "test_corpus"))
    for cand in candidates:
        if cand and os.path.isdir(os.path.join(cand, "corpus_gui")):
            return cand
    return candidates[-1]


CORPUS = _corpus_root()
KNOWLEDGE = os.path.join(os.path.dirname(_HERE), "knowledge", "programs.yaml")

TRUTH = {
    "PHASER": "phenix.phaser", "autobuild": "phenix.autobuild",
    "autosol": "phenix.autosol", "find_reference": "phenix.find_reference",
    "ligandfit": "phenix.ligandfit", "molprobity": "phenix.molprobity",
    "predict_and_build": "phenix.predict_and_build",
    "refine": "phenix.refine", "xtriage": "phenix.xtriage",
}


def expected(stem):
    for key in sorted(TRUTH, key=len, reverse=True):
        if stem.lower().startswith(key.lower()):
            return TRUTH[key]
    return None



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

PASS = []
FAIL = []
SKIP = []


def ok(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print("  %s: %s%s" % ("PASS" if condition else "FAIL", name,
                          "  -- " + detail if detail and not condition
                          else ""))


# ------------------------------------------------------------------ tests

def t_registry_is_loaded_not_hardcoded(mod):
    """A hardcoded list drifts. The first draft of the resolver in
    orientation_strip.py held 15 names against programs.yaml's 23 and
    would have refused a `polder` log as 'not a Phenix program'."""
    reg = mod.load_registry(KNOWLEDGE)
    ok("registry loads from programs.yaml", len(reg) >= 23,
       "got %d names" % len(reg))
    ok("registry contains programs the hardcoded list missed",
       all(p in reg for p in ("phenix.polder", "phenix.pdbtools",
                              "phenix.map_to_model")))


def t_caller_wins_and_is_authoritative(mod):
    r = mod.resolve_program(caller="phenix.refine", log_text="",
                            log_path="xtriage_1_p9.log", registry_path=KNOWLEDGE)
    ok("a caller-supplied name wins over the file name",
       r.name == "phenix.refine")
    ok("caller is marked authoritative", r.source == "caller" and
       r.authoritative is True)


def t_banner_beats_filename(mod):
    log = "Starting phenix.molprobity\non Sat by terwill\n"
    r = mod.resolve_program(caller=None, log_text=log,
                            log_path="refine_3_esterase.log",
                            registry_path=KNOWLEDGE)
    ok("the log's own banner beats the file name",
       r.name == "phenix.molprobity")
    ok("banner is marked authoritative", r.source == "banner" and
       r.authoritative is True)


def t_filename_is_inferred_not_authoritative(mod):
    """Registry membership proves a name EXISTS, not that it is RIGHT.
    A user can rename a file. 'Never a guess' was too strong."""
    r = mod.resolve_program(caller=None, log_text="no banner here",
                            log_path="refine_3_esterase.log",
                            registry_path=KNOWLEDGE)
    ok("a run number is stripped from the file name",
       r.name == "phenix.refine")
    ok("file name is marked INFERRED, not authoritative",
       r.source == "filename" and r.authoritative is False)


def t_refusal_is_explicit(mod):
    r = mod.resolve_program(caller=None, log_text="nothing useful",
                            log_path="wibble_flarn.log",
                            registry_path=KNOWLEDGE)
    ok("an unresolvable name refuses rather than guessing", r.name is None)
    ok("refusal carries source='none'", r.source == "none")
    ok("refusal is not authoritative", r.authoritative is False)


def t_never_emits_a_nonexistent_program(mod):
    """Defect 4. The shipped code produced phenix.xtriage_1 on 20 of 20."""
    reg = mod.load_registry(KNOWLEDGE)
    bad = []
    for shape in ("corpus", "corpus_gui"):
        d = os.path.join(CORPUS, shape)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            text = open(os.path.join(d, f), errors="replace").read()
            r = mod.resolve_program(caller=None, log_text=text, log_path=f,
                                    registry_path=KNOWLEDGE)
            if r.name is not None and r.name not in reg:
                bad.append((f, r.name))
    ok("no emitted name is outside the registry", not bad, str(bad[:3]))


def t_corpus_invariant_twenty_of_twenty(mod):
    """The corpus-level invariant for this step."""
    for shape in ("corpus", "corpus_gui"):
        d = os.path.join(CORPUS, shape)
        if not os.path.isdir(d):
            SKIP.append("%s not present" % shape)
            print("  SKIP: %s not present" % shape)
            continue
        right = wrong = refused = 0
        for f in sorted(os.listdir(d)):
            text = open(os.path.join(d, f), errors="replace").read()
            r = mod.resolve_program(caller=None, log_text=text, log_path=f,
                                    registry_path=KNOWLEDGE)
            exp = expected(os.path.splitext(f)[0])
            if r.name == exp:
                right += 1
            elif r.name is None:
                refused += 1
            else:
                wrong += 1
        ok("%s: 20 right, 0 wrong, 0 refused" % shape,
           (right, wrong, refused) == (20, 0, 0),
           "got %d right, %d wrong, %d refused" % (right, wrong, refused))


def t_no_llm_value_can_reach_the_name(mod):
    """The guard that actually matters. Asserted on the source, because
    the failure mode is a future maintainer wiring a summary-derived
    string back in -- which is what the shipped code does today."""
    # Inspect CODE, not prose. The first version grepped raw source and
    # failed on the module's own docstring, which names
    # processed_log_dict while explaining that it does NOT use it. That
    # is error class 3 -- matching against text we wrote ourselves --
    # occurring inside the guard written to prevent a different defect.
    import ast
    path = mod.__file__.replace(".pyc", ".py")
    tree = ast.parse(open(path, errors="replace").read())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef)):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body.pop(0)
    code = ast.dump(tree)

    forbidden = ("processed_log_dict", "phenix_program", "summarize_log",
                 "get_log_info", "invoke")
    hits = [w for w in forbidden if w in code]
    ok("the resolver never consults an LLM-derived value", not hits,
       "found %s in executable code" % hits)


_require('PHIL_CORPUS_DIR', CORPUS, 'the corpus')


def run_all_tests():
    try:
        import program_identity as mod
    except ImportError as error:
        print("  FAIL: program_identity does not exist yet -- %s" % error)
        print("\n0 passed, 1 failed, 0 skipped, 1 total")
        print("*** Implementation has not been written.")
        return False

    for fn in (t_registry_is_loaded_not_hardcoded,
               t_caller_wins_and_is_authoritative,
               t_banner_beats_filename,
               t_filename_is_inferred_not_authoritative,
               t_refusal_is_explicit,
               t_never_emits_a_nonexistent_program,
               t_corpus_invariant_twenty_of_twenty,
               t_no_llm_value_can_reach_the_name):
        fn(mod)

    # House summary format -- VERIFY.py parses this line and reported
    # "NO SUMMARY" until it matched. A suite that does not announce its
    # own results the way the harness reads them is a suite the harness
    # cannot see.
    total = len(PASS) + len(FAIL) + len(SKIP)
    print("\n%d passed, %d failed, %d skipped, %d total"
          % (len(PASS), len(FAIL), len(SKIP), total))
    if SKIP:
        print("*** %d SKIPPED, and a skip is not a pass: %s"
              % (len(SKIP), ", ".join(SKIP)))
    if FAIL:
        print("*** FAILED: %s" % ", ".join(FAIL))
    return not FAIL


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
