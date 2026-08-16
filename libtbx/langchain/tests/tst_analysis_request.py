"""Tests for the single-call analysis request -- Step 2.

WHAT DEFECTS THIS EXISTS FOR

  defect 1  The analysis step reads a ~1 KB summary, not the log.
            6 of 7 logs retain under 10%.
  defect 2  _custom_log_chunker drops everything after `Citations`
            when "Files are in the directory" appears. 7 of 20 logs;
            59% of a 304 KB log.
  defect 4  A fabricated program name is written into the prompt.
            20 of 20 logs named a program that does not exist.
  defect 6  The summariser can emit more than it was given, and does so
            reproducibly: 382,721 chars from a 156,951-char log, one
            line 4,328 times; 2,026,595 chars from 312,011, in a single
            line of 1.88 million characters.

The request built here contains the log, whole and unmodified. There is
no summariser and no chunker to exhibit any of them.

THE INVARIANT, AND WHY IT IS NOT A RATIO

An earlier draft of this test asserted `len(payload) >= 0.99 *
len(log)`. A reviewer pointed out that mangled or partly duplicated
input passes that. The assertion here is that **the log appears in the
payload exactly once, byte for byte** -- checked by count and by hash,
so neither truncation nor duplication nor mutation can pass.
"""
from __future__ import absolute_import, division, print_function

import hashlib
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
    env = os.environ.get("PHIL_CORPUS_DIR")
    for cand in ([env, os.path.dirname((env or "").rstrip(os.sep))] +
                 [os.path.join(os.path.dirname(_HERE), "test_corpus")]):
        if cand and os.path.isdir(os.path.join(cand, "corpus_gui")):
            return cand
    return os.path.join(os.path.dirname(_HERE), "test_corpus")


CORPUS = os.path.join(_corpus_root(), "corpus_gui")
HAVE_CORPUS = os.path.isdir(CORPUS)
KNOWLEDGE = os.path.join(os.path.dirname(_HERE), "knowledge", "programs.yaml")


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


# ------------------------------------------------------------------ tests

def t_log_appears_exactly_once_byte_for_byte(mod, ident):
    """Defects 1 and 2 in one assertion, on every corpus log."""
    bad = []
    for name in sorted(os.listdir(CORPUS)):
        log = open(os.path.join(CORPUS, name), errors="replace").read()
        req = mod.build_analysis_request(
            log_text=log,
            identity=ident.resolve_program(log_text=log, log_path=name,
                                           registry_path=KNOWLEDGE))
        if req.payload.count(log) != 1:
            bad.append((name, req.payload.count(log)))
    ok("the log appears exactly once in every payload", not bad,
       "%s" % bad[:3])


def t_hash_of_the_embedded_log_matches_the_source(mod, ident):
    """Count alone would miss mutation inside a single copy."""
    log = open(os.path.join(CORPUS, "refine_3_esterase.log"),
               errors="replace").read()
    req = mod.build_analysis_request(
        log_text=log, identity=ident.resolve_program(
            log_text=log, log_path="refine_3_esterase.log",
            registry_path=KNOWLEDGE))
    ok("the embedded log hashes to the source",
       req.log_sha256 == hashlib.sha256(log.encode("utf-8",
                                                   "replace")).hexdigest())
    ok("the payload actually contains that block", req.log_sha256 in
       hashlib.sha256(mod.extract_log_block(req.payload).encode(
           "utf-8", "replace")).hexdigest() or
       mod.extract_log_block(req.payload) == log)


def t_the_largest_log_is_not_truncated(mod, ident):
    """Defect 2's worst case: 59% of this log never reached the model."""
    name = "predict_and_build_2_bromodomain.log"
    path = os.path.join(CORPUS, name)
    if not os.path.exists(path):
        SKIPPED.append("largest log absent")
        print("  SKIP: %s not present" % name)
        return
    log = open(path, errors="replace").read()
    req = mod.build_analysis_request(
        log_text=log, identity=ident.resolve_program(
            log_text=log, log_path=name, registry_path=KNOWLEDGE))
    ok("the 304 KB log reaches the payload whole",
       req.payload.count(log) == 1 and len(log) > 300000,
       "log %d chars" % len(log))


def t_no_fabricated_program_name(mod, ident):
    """Defect 4. The shipped header named a nonexistent program on 20
    of 20; an unresolved identity must produce no header at all."""
    import yaml
    registry = set(yaml.safe_load(open(KNOWLEDGE))) | {"phenix.find_reference"}
    bad = []
    for name in sorted(os.listdir(CORPUS)):
        log = open(os.path.join(CORPUS, name), errors="replace").read()
        idn = ident.resolve_program(log_text=log, log_path=name,
                                    registry_path=KNOWLEDGE)
        req = mod.build_analysis_request(log_text=log, identity=idn)
        import re
        for found in re.findall(r"phenix\.[a-z_0-9]+", req.header or ""):
            if found not in registry:
                bad.append((name, found))
    ok("no payload header names a program outside the registry", not bad,
       "%s" % bad[:3])


def t_unresolved_identity_yields_no_header(mod, ident):
    idn = ident.resolve_program(log_text="nothing", log_path="wibble.log",
                                registry_path=KNOWLEDGE)
    req = mod.build_analysis_request(log_text="nothing", identity=idn)
    ok("an unresolved program produces no header line", not req.header)
    ok("and the payload says so rather than inventing one",
       "not stated" in req.payload.lower())


def t_payload_is_never_larger_than_the_log_plus_a_bound(mod, ident):
    """Defect 6's guard, applied to the thing that replaces it. The
    summariser emitted 2,026,595 chars from 312,011. Nothing here can
    grow the input, and this asserts it rather than assuming it."""
    worst = 0
    for name in sorted(os.listdir(CORPUS)):
        log = open(os.path.join(CORPUS, name), errors="replace").read()
        req = mod.build_analysis_request(
            log_text=log, identity=ident.resolve_program(
                log_text=log, log_path=name, registry_path=KNOWLEDGE))
        worst = max(worst, len(req.payload) - len(log))
    # Bound set from the measurement, not from a round number. Observed
    # worst case is 863 chars across the corpus; 8000 was 9x that and
    # would have passed a payload carrying an extra 7 KB of anything.
    ok("payload never exceeds the log by more than the instruction block",
       worst < 1200, "worst overhead %d chars" % worst)


def t_verify_catches_a_corrupted_request(mod, ident):
    """verify() is the guard between a construction bug and a report
    built on the wrong text. It was written and never exercised -- dead
    code in the position of a safety check is worse than no check."""
    log = "some log text\nwith lines\n"
    idn = ident.resolve_program(log_text="Starting phenix.refine\n",
                                log_path="refine_1.log",
                                registry_path=KNOWLEDGE)
    req = mod.build_analysis_request(log_text=log, identity=idn)
    good, problems = mod.verify(req)
    ok("verify() passes a well-formed request", good, "%s" % problems)

    req.payload = req.payload.replace("with lines", "with LINES")
    bad, problems = mod.verify(req)
    ok("verify() catches a mutated log block", not bad and problems)

    req2 = mod.build_analysis_request(log_text=log, identity=idn)
    req2.payload = req2.payload + "\n" + mod.LOG_OPEN + "\n"
    bad2, problems2 = mod.verify(req2)
    ok("verify() catches a duplicated log block", not bad2 and problems2)


def t_no_summariser_or_retriever_is_reachable(mod):
    """Defects 1, 5 and 6 cannot occur if the components are absent.
    Asserted on the AST with docstrings stripped -- the first version of
    a guard like this in tst_program_identity failed on prose that named
    the very thing it forbids."""
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
    forbidden = ("summarize_log", "_custom_log_chunker", "retriever",
                 "load_persistent_db", "processed_log_dict")
    hits = [w for w in forbidden if w in code]
    ok("no summariser, chunker or retriever is reachable", not hits,
       "found %s" % hits)


_require('PHIL_CORPUS_DIR', CORPUS, 'corpus_gui')


def run_all_tests():
    if not HAVE_CORPUS:
        # The corpus is not shipped into the Phenix tree: it is 20 real
        # logs carrying absolute home paths from two machines (see the
        # plan's privacy review).  Point PHIL_CORPUS_DIR at it.
        print("  SKIP: corpus not found. Set PHIL_CORPUS_DIR to the")
        print("        directory holding corpus/ and corpus_gui/.")
        print("\n0 passed, 0 failed, 1 skipped, 1 total")
        print("*** 1 SKIPPED, and a skip is not a pass.")
        return True
    try:
        import analysis_request as mod
        import program_identity as ident
    except ImportError as error:
        print("  FAIL: analysis_request does not exist yet -- %s" % error)
        print("\n0 passed, 1 failed, 0 skipped, 1 total")
        return False

    t_log_appears_exactly_once_byte_for_byte(mod, ident)
    t_hash_of_the_embedded_log_matches_the_source(mod, ident)
    t_the_largest_log_is_not_truncated(mod, ident)
    t_no_fabricated_program_name(mod, ident)
    t_unresolved_identity_yields_no_header(mod, ident)
    t_payload_is_never_larger_than_the_log_plus_a_bound(mod, ident)
    t_verify_catches_a_corrupted_request(mod, ident)
    t_no_summariser_or_retriever_is_reachable(mod)

    total = len(PASS) + len(FAIL) + len(SKIPPED)
    print("\n%d passed, %d failed, %d skipped, %d total"
          % (len(PASS), len(FAIL), len(SKIPPED), total))
    if SKIPPED:
        print("*** %d SKIPPED, and a skip is not a pass" % len(SKIPPED))
    if FAIL:
        print("*** FAILED: %s" % ", ".join(FAIL))
    return not FAIL


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
