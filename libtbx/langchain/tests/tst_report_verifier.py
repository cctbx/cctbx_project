"""Tests for the report verifier -- Step 3.

WHAT THIS EXISTS FOR

Phase 0a swept twenty Arm E reports and found **zero fabricated
figures**. That is a guard against a problem we do not currently have,
and it stays cheap as models change. The value is in the rarity: a
verifier that fires often teaches users which warnings to ignore, so
**the false-positive rate is the thing to measure, not the true-positive
rate.**

THE CASE THAT MADE ROUNDING NON-OPTIONAL

In the first real report examined, two of the figures were derived
rather than quoted:

    "Resolution (1.74 Å)"          from  1.74434
    "Twin fraction estimates (<0.03)"  a bound over four values
                                       whose maximum is 0.029

**Both correct. A naive check flags both.** Two false flags on the first
report examined is exactly the failure mode above, so rounding and
bounds are requirements here, not refinements.

WHAT IT DOES NOT REACH

The wrong-phasing recommendation of defect 1 passes every check in this
file, and so does "the enantiomorphic space groups I 4 and I 41" -- a
wrong domain claim in fluent language, in a real shipped report. **The
badge says "Data-source checks passed", never "Analysis verified".**
Even "the numbers are checked" overstates it: if 0.26 appears labelled
R-work when the log says R-free, the number exists and the claim is
wrong.
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


def _root():
    env = os.environ.get("PHIL_CORPUS_DIR")
    for cand in ([env, os.path.dirname((env or "").rstrip(os.sep))] +
                 [os.path.join(os.path.dirname(_HERE), "test_corpus")]):
        if cand and os.path.isdir(os.path.join(cand, "corpus_gui")):
            return cand
    return os.path.join(os.path.dirname(_HERE), "test_corpus")


CORPUS = os.path.join(_root(), "corpus_gui")
REPORTS = os.path.join(os.path.dirname(_HERE), "study_reports")
CAPTURED = os.path.join(os.path.dirname(_HERE), "captured_runs")
REPORTS = os.environ.get("PHIL_REPORTS_DIR", REPORTS)
CAPTURED = os.environ.get("PHIL_CAPTURED_DIR", CAPTURED)
HAVE_DATA = os.path.isdir(REPORTS) and os.path.isdir(CORPUS)


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


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  %s: %s%s" % ("PASS" if cond else "FAIL", name,
                          "  -- " + detail if detail and not cond else ""))


# ------------------------------------------------------------------ tests

def t_exact_figures_are_supported(mod):
    log = "Resolution range: 28.4872 1.74434\nR-free 0.2601\n"
    flags = mod.check_numbers("R-free is 0.2601.", log)
    ok("a figure quoted verbatim raises no flag", not flags, "%s" % flags)


def t_rounding_is_accepted(mod):
    """1.74 from 1.74434 -- observed in a real report."""
    log = "Resolution range: 28.4872 1.74434\n"
    flags = mod.check_numbers("The resolution is 1.74 A.", log)
    ok("a rounded figure raises no flag", not flags, "%s" % flags)


def t_a_derived_bound_is_accepted(mod):
    """'<0.03' over 0.00, 0.008, 0.029, 0.022 -- also observed."""
    log = ("L-test: 0.00\nBritton: 0.008\nH-test: 0.029\ncum: 0.022\n")
    flags = mod.check_numbers(
        "Twin fraction estimates (<0.03) are insignificant.", log)
    ok("a bound over log values raises no flag", not flags, "%s" % flags)


def t_a_fabricated_figure_is_flagged(mod):
    """The control. Without this the suite cannot tell a working check
    from one that never fires."""
    log = "Resolution range: 28.4872 1.74434\n"
    flags = mod.check_numbers("The R-free was 0.187.", log)
    ok("a figure absent from the log IS flagged", len(flags) == 1,
       "%s" % flags)


def t_years_and_versions_are_not_treated_as_measurements(mod):
    log = "Starting phenix.refine\non Mon Jun 15 19:21:26 2026 by terwill\n"
    flags = mod.check_numbers(
        "The run took place in 2026 and used 3 macro-cycles.", log)
    ok("a year present in the log is not flagged", "2026" not in str(flags))


def t_false_positive_rate_on_real_reports(mod):
    """The measurement that decides whether this can be exposed to
    users. Run over every real report in the package against its own
    log. Every flag here is a candidate false positive: Phase 0a's
    sweep found zero fabricated figures in twenty reports."""
    pairs = []
    for name in sorted(os.listdir(REPORTS)):
        if not name.endswith("__B.txt") and not name.endswith("__E.txt"):
            continue
        stem = name.rsplit("__", 1)[0]
        log = os.path.join(CORPUS, stem + ".log")
        if os.path.exists(log):
            pairs.append((log, os.path.join(REPORTS, name)))
    if not pairs:
        SKIPPED.append("no report/log pairs found")
        print("  SKIP: no report/log pairs found")
        return

    total_flags = flagged_reports = 0
    for log_path, rep_path in pairs:
        flags = mod.check_numbers(
            open(rep_path, errors="replace").read(),
            open(log_path, errors="replace").read())
        total_flags += len(flags)
        flagged_reports += 1 if flags else 0
    rate = flagged_reports / float(len(pairs))
    print("      %d reports, %d flagged, %d flags total, rate %.0f%%"
          % (len(pairs), flagged_reports, total_flags, 100 * rate))
    ok("false-positive rate on real reports is under 25%", rate < 0.25,
       "%.0f%% of reports flagged" % (100 * rate))


def t_the_two_known_derived_figures_do_not_flag(mod):
    """The specific case from the first real report examined."""
    log_path = os.path.join(CORPUS, "xtriage_1_p9.log")
    rep_path = os.path.join(CAPTURED, "xtriage_1_p9.analysis.md")
    if not (os.path.exists(log_path) and os.path.exists(rep_path)):
        SKIPPED.append("xtriage capture absent")
        print("  SKIP: xtriage capture absent")
        return
    flags = mod.check_numbers(open(rep_path, errors="replace").read(),
                              open(log_path, errors="replace").read())
    bad = [f for f in flags if "1.74" in str(f) or "0.03" in str(f)]
    ok("neither known derived figure is flagged", not bad, "%s" % bad)


def t_program_existence_check(mod):
    """check_program_names BLOCKS rather than flags, and was written
    without a test -- the same defect found in analysis_request.verify()
    one turn earlier. One line of it would have caught `phenix.royal_`
    on its first appearance."""
    registry = {"phenix.refine", "phenix.xtriage"}
    ok("a real program name is accepted",
       not mod.check_program_names("Run phenix.refine next.", registry))
    ok("an invented program name is caught",
       mod.check_program_names("Try phenix.royal_ instead.", registry)
       == ["phenix.royal_"])
    ok("several are all reported",
       len(mod.check_program_names("phenix.wibble and phenix.flarn",
                                   registry)) == 2)


def t_summarise_wording(mod):
    """The operator-facing line. The badge must never say more than
    'Data-source checks passed' -- overstating it would repeat, at
    product level, the error found in the reports."""
    ok("a clean report reports passing", mod.summarise([], []) ==
       "Data-source checks passed.")
    ok("a blocked report says BLOCK and names the program",
       mod.summarise([], ["phenix.royal_"]).startswith("BLOCK")
       and "phenix.royal_" in mod.summarise([], ["phenix.royal_"]))
    ok("the summary never claims the analysis was verified",
       "verified" not in mod.summarise([], []).lower()
       and "verified" not in mod.summarise([1], []).lower())
    ok("a block outranks a flag",
       mod.summarise([1, 2], ["phenix.x"]).startswith("BLOCK"))


def t_held_out_false_positive_rate(mod):
    """The 10% measured on the study reports is IN-SAMPLE: the checker
    was tuned by looking at what it flagged there. The shipped-pipeline
    reports in captured_runs/ were never consulted while tuning."""
    pairs = []
    for f in sorted(os.listdir(CAPTURED)):
        if ".analysis" not in f:
            continue
        stem = f.split(".analysis")[0]
        log = os.path.join(CORPUS, stem + ".log")
        if os.path.exists(log):
            pairs.append((log, os.path.join(CAPTURED, f)))
    if not pairs:
        SKIPPED.append("no held-out reports")
        print("  SKIP: no held-out reports")
        return
    flagged = sum(1 for lp, rp in pairs
                  if mod.check_numbers(open(rp, errors="replace").read(),
                                       open(lp, errors="replace").read()))
    rate = flagged / float(len(pairs))
    print("      held-out: %d reports, %d flagged, rate %.0f%%"
          % (len(pairs), flagged, 100 * rate))
    ok("held-out rate is under 30%", rate < 0.30,
       "%.0f%% -- in-sample was 10%%" % (100 * rate))


_require('PHIL_CORPUS_DIR', CORPUS, 'corpus_gui')
_require('PHIL_REPORTS_DIR', REPORTS, 'the reports')
_require('PHIL_CAPTURED_DIR', CAPTURED, 'the captures')


def run_all_tests():
    if not HAVE_DATA:
        print("  SKIP: reports or corpus not found. Set PHIL_CORPUS_DIR")
        print("        and PHIL_REPORTS_DIR (and PHIL_CAPTURED_DIR).")
        print("\n0 passed, 0 failed, 1 skipped, 1 total")
        print("*** 1 SKIPPED, and a skip is not a pass.")
        return True
    try:
        import report_verifier as mod
    except ImportError as error:
        print("  FAIL: report_verifier does not exist yet -- %s" % error)
        print("\n0 passed, 1 failed, 0 skipped, 1 total")
        return False

    for fn in (t_exact_figures_are_supported, t_rounding_is_accepted,
               t_a_derived_bound_is_accepted, t_a_fabricated_figure_is_flagged,
               t_years_and_versions_are_not_treated_as_measurements,
               t_the_two_known_derived_figures_do_not_flag,
               t_program_existence_check, t_summarise_wording,
               t_false_positive_rate_on_real_reports,
               t_held_out_false_positive_rate):
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
