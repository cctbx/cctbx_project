"""
B1 — LLM Decision Test: stop_after_requested Extraction (v117)

Five advice strings, each exercising the v117 stop_refactor pattern:
the LLM (or regex backstop) must produce
`stop_conditions.stop_after_requested = True` plus a sensible
`after_program`.

This is the EXTRACTION layer counterpart to B2 (planning layer).
B2 verifies the LLM picks STOP when the gated wipe puts STOP in
valid_programs; B1 verifies the gate gets set in the first place.

The fifth scenario — 'Refine. Stop.' (period-separated) — exercises
the v117 Step B period-stop patch added to _POSITIVE_STOP_AFTER_PATTERNS.

All five are RELIABILITY tests at threshold 0.8.  The regex backstop
should make these very reliable; below 0.8 indicates either the
extractor LLM is regressing or the regex helper has a bug.

Run:
    phenix.python tests/llm/tst_b1_stop_after_requested.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import (Scenario, make_directive_extraction_run_fn)


# =====================================================================
# Expected-output helpers
# =====================================================================

def _density_modify_family_with_stop(directives):
    """B1.1 specific: 'density modify and stop' is ambiguous about
    which density-modification program is meant without a cryo-EM/X-ray
    context cue.  Accept any density-modification-family choice as long
    as stop_after_requested=True is set.
    """
    sc = (directives or {}).get("stop_conditions") or {}
    sar = sc.get("stop_after_requested")
    ap = sc.get("after_program")
    if sar is not True:
        return (False, "stop_after_requested=%r (expected True). "
                       "Got stop_conditions=%r" % (sar, sc))
    if ap is None:
        return (True, "stop_after_requested=True, after_program omitted "
                      "(downstream resolver handles)")
    density_family = (
        "phenix.resolve_cryo_em", "phenix.resolve",
        "phenix.density_modification", "phenix.parrot",
        "phenix.solvent_modify",
    )
    if ap in density_family:
        return (True, "stop_after_requested=True, after_program=%s" % ap)
    return (False, "stop_after_requested=True but after_program=%r "
                   "not in density-modification family" % ap)


def _stop_after_requested_with_program(expected_program):
    """Builder: assert stop_after_requested=True AND after_program
    matches the expected value.
    """
    def check(directives):
        sc = (directives or {}).get("stop_conditions") or {}
        sar = sc.get("stop_after_requested")
        ap = sc.get("after_program")
        if sar is not True:
            return (False, "stop_after_requested=%r (expected True). "
                           "Got stop_conditions=%r" % (sar, sc))
        if ap != expected_program:
            return (False, "after_program=%r (expected %r). "
                           "stop_after_requested correctly True."
                           % (ap, expected_program))
        return (True, "stop_after_requested=True, after_program=%s" % ap)
    return check


def _stop_after_requested_any_program(directives):
    """For inputs where after_program selection is ambiguous (e.g. the
    program isn't directly named in the advice), only assert that
    stop_after_requested is True.
    """
    sc = (directives or {}).get("stop_conditions") or {}
    if sc.get("stop_after_requested") is True:
        return (True, "stop_after_requested=True (after_program=%s)"
                      % sc.get("after_program"))
    return (False, "stop_after_requested=%r" % sc.get("stop_after_requested"))


# =====================================================================
# Scenarios
# =====================================================================

def build_scenarios():
    """B1 — five stop_after_requested reliability scenarios."""
    return [
        Scenario(
            name="b1_density_modify_and_stop",
            description=(
                "Canonical case: 'density modify and stop' must set "
                "stop_after_requested=True with a density-modification-"
                "family after_program (resolve_cryo_em / resolve / "
                "density_modification / parrot are all acceptable; "
                "the raw advice doesn't specify cryo-EM vs X-ray)."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input="density modify and stop",
            expected_fn=_density_modify_family_with_stop,
        ),
        Scenario(
            name="b1_refine_and_stop",
            description=(
                "'Refine the model and stop' must set "
                "stop_after_requested=True with after_program=phenix.refine."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input="Refine the model and stop.",
            expected_fn=_stop_after_requested_with_program("phenix.refine"),
        ),
        Scenario(
            name="b1_xtriage_and_stop",
            description=(
                "'Run xtriage and stop' must set stop_after_requested=True "
                "with after_program=phenix.xtriage."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input="Run xtriage and stop.",
            expected_fn=_stop_after_requested_with_program("phenix.xtriage"),
        ),
        Scenario(
            name="b1_period_stop_refine",
            description=(
                "Period-separated stop phrasing: 'Refine. Stop.'  This "
                "exercises the v117 Step B period-stop patch added to "
                "_POSITIVE_STOP_AFTER_PATTERNS (the bare \\bstop\\b "
                "match that the new helper has to cover via the "
                "r'\\.\\s+stop\\b' pattern)."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input="Refine. Stop.",
            expected_fn=_stop_after_requested_with_program("phenix.refine"),
        ),
        Scenario(
            name="b1_just_run_mtriage",
            description=(
                "'Just run mtriage' — the 'just run X' phrasing is one "
                "of the 11 positive patterns in _is_stop_after_requested. "
                "Must set stop_after_requested=True."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input="Just run mtriage.",
            expected_fn=_stop_after_requested_any_program,
        ),
    ]


# =====================================================================
# Standalone runner
# =====================================================================

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    from framework import (AVAILABLE_PROVIDERS, build_summary,
                           create_log_dir, run_scenario_against_providers,
                           write_run_summary)
    import time

    if not AVAILABLE_PROVIDERS:
        sys.stderr.write(
            "No LLM provider API keys configured.\n")
        return 2

    providers = list(AVAILABLE_PROVIDERS)
    scenarios = build_scenarios()
    run_one = make_directive_extraction_run_fn()
    log_dir = create_log_dir(phase="phase1-b1")

    print()
    print("LLM Decision Tests — B1 stop_after_requested (v117)")
    print("=" * 65)
    print("Available providers: %s" % ", ".join(providers))
    print("Log directory:       %s" % log_dir)
    print()

    all_verdicts = []
    t0 = time.time()
    for scenario in scenarios:
        verdicts = run_scenario_against_providers(
            scenario, providers, run_one, log_dir=log_dir)
        all_verdicts.extend(verdicts)
    wall_time = time.time() - t0

    summary = build_summary(
        all_verdicts, providers_tested=providers,
        providers_skipped=[], wall_time_s=wall_time,
        phase="phase1-b1")
    write_run_summary(log_dir, summary)

    print()
    print("Summary: %d PASS, %d FAIL, %d ERROR, %d SKIP"
          % (summary["totals"]["passed"], summary["totals"]["failed"],
             summary["totals"]["errored"], summary["totals"]["skipped"]))
    print("Log:                 %s" % log_dir)

    return 0 if summary["totals"]["failed"] == 0 and \
                summary["totals"]["errored"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
