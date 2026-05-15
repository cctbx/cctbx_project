"""
LLM Decision Tests — Directive Extraction (D2, Phase 1)

Tests 6 directive-extraction scenarios against each available LLM
provider, using the framework's early-termination runner.

Run from a parent directory with framework.py importable, e.g.:

    cd $PHENIX/modules/cctbx_project/libtbx/langchain
    phenix.python tests/llm/tst_directive_extraction.py

Or via the top-level CLI runner:

    phenix.python tests/llm/run_llm_tests.py \\
        --suite directive_extraction

Each scenario classifies as either:
  - "capability": ≥1 pass over max_runs = success.  The LLM is shown
                  to be CAPABLE of the right answer.
  - "reliability": ≥(threshold × max_runs) passes = success.  The LLM
                   RELIABLY produces the right answer.

The AF_7mjs Bug A regression is a reliability test at threshold 0.9.
The other five are capability tests.
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Make the framework importable when invoked directly
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import (
    Scenario, load_fixture,
    make_directive_extraction_run_fn)


# =====================================================================
# Expected-output helper functions
# =====================================================================
#
# Each helper returns a (passed: bool, why: str) tuple given the
# parsed directives dict.  "why" is a short explanation for the log.

def _no_fabricated_after_program(directives):
    """AF_7mjs regression check: no fabricated after_program or
    prefer_programs."""
    stop = (directives or {}).get("stop_conditions") or {}
    wf = (directives or {}).get("workflow_preferences") or {}
    if stop.get("after_program"):
        return (False, "fabricated after_program=%s" % stop["after_program"])
    if wf.get("prefer_programs"):
        return (False,
                "fabricated prefer_programs=%s" % wf["prefer_programs"])
    return (True, "no after_program, no prefer_programs")


def _after_program_equals(expected):
    """Builder for "after_program must equal X"."""
    def check(directives):
        actual = ((directives or {}).get("stop_conditions") or {}).get(
            "after_program")
        if actual == expected:
            return (True, "after_program=%s" % actual)
        return (False, "expected after_program=%s, got %r"
                       % (expected, actual))
    return check


def _has_resolution_and_space_group(expected_res, expected_sg):
    """Resolution and space_group are emitted under program_settings.

    The prompt's CRITICAL rule says they go under "default" scope, but
    LLMs sometimes put them under specific program scopes.  Accept
    EITHER as long as both values appear somewhere in program_settings.
    """
    def check(directives):
        ps = (directives or {}).get("program_settings") or {}
        if not isinstance(ps, dict):
            return (False, "program_settings is not a dict: %r" % type(ps))

        found_res = None
        found_sg = None

        # Search every scope (default + per-program scopes)
        for scope_name, scope in ps.items():
            if not isinstance(scope, dict):
                continue
            if "resolution" in scope and found_res is None:
                found_res = scope["resolution"]
            if "space_group" in scope and found_sg is None:
                found_sg = scope["space_group"]

        # Normalize for comparison
        try:
            res_match = (found_res is not None and
                         abs(float(found_res) - expected_res) < 0.01)
        except (TypeError, ValueError):
            res_match = False
        sg_match = (found_sg is not None and
                    _normalize_sg(found_sg) == _normalize_sg(expected_sg))

        if res_match and sg_match:
            return (True, "resolution=%s, space_group=%s"
                           % (found_res, found_sg))
        problems = []
        if not res_match:
            problems.append("expected resolution=%s, got %r"
                            % (expected_res, found_res))
        if not sg_match:
            problems.append("expected space_group=%s, got %r"
                            % (expected_sg, found_sg))
        return (False, "; ".join(problems))
    return check


def _normalize_sg(sg):
    """Collapse whitespace in space group strings ('P 32 2 1' == 'P3221')."""
    if sg is None:
        return ""
    return "".join(str(sg).split()).upper()


def _skip_programs_contains(expected):
    """workflow_preferences.skip_programs must contain expected."""
    def check(directives):
        wf = (directives or {}).get("workflow_preferences") or {}
        skip = wf.get("skip_programs") or []
        if not isinstance(skip, list):
            return (False, "skip_programs is not a list: %r" % type(skip))
        if expected in skip:
            return (True, "skip_programs contains %s (got %r)"
                           % (expected, skip))
        return (False, "expected %s in skip_programs, got %r"
                       % (expected, skip))
    return check


def _wants_validation_only(directives):
    """workflow_preferences.wants_validation_only must be truthy."""
    wf = (directives or {}).get("workflow_preferences") or {}
    val = wf.get("wants_validation_only")
    if val:
        return (True, "wants_validation_only=%r" % val)
    return (False, "wants_validation_only not set; "
                   "workflow_preferences=%r" % wf)


def _use_mr_sad(directives):
    """workflow_preferences.use_mr_sad must be truthy."""
    wf = (directives or {}).get("workflow_preferences") or {}
    val = wf.get("use_mr_sad")
    if val:
        return (True, "use_mr_sad=%r" % val)
    return (False, "use_mr_sad not set; workflow_preferences=%r" % wf)


# =====================================================================
# Scenarios
# =====================================================================

def build_scenarios():
    """Build the list of D2 scenarios.

    Wrapped in a function so fixtures are loaded lazily — avoids
    issues with import order in test discovery.
    """
    return [
        Scenario(
            name="af7mjs_bug_a_regression",
            description=(
                "AF_7mjs preprocessed advice — directives must NOT "
                "contain a fabricated after_program or prefer_programs. "
                "v116.19a regression test (pre-fix rate was ~7%; "
                "post-fix should be ~100%)."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.9,
            max_runs=5,
            input=load_fixture("af7mjs_preprocessed.txt"),
            expected_fn=_no_fabricated_after_program,
        ),
        Scenario(
            name="explicit_stop_after_phaser",
            description=(
                "Advice with explicit 'Stop after running phaser' "
                "should produce after_program=phenix.phaser."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input=load_fixture("explicit_stop_after_phaser.txt"),
            expected_fn=_after_program_equals("phenix.phaser"),
        ),
        Scenario(
            name="resolution_and_space_group",
            description=(
                "Advice with explicit resolution=2.5 and space_group="
                "P 32 2 1 should extract both into program_settings."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input=load_fixture("resolution_and_space_group.txt"),
            expected_fn=_has_resolution_and_space_group(2.5, "P 32 2 1"),
        ),
        Scenario(
            name="skip_programs",
            description=(
                "'Don't run phenix.mtriage' should produce "
                "skip_programs containing phenix.mtriage."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input=load_fixture("skip_programs.txt"),
            expected_fn=_skip_programs_contains("phenix.mtriage"),
        ),
        Scenario(
            name="validation_only",
            description=(
                "'Validate this structure' as primary goal should "
                "produce wants_validation_only=True."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input=load_fixture("validation_only.txt"),
            expected_fn=_wants_validation_only,
        ),
        Scenario(
            name="mr_sad_workflow",
            description=(
                "'Use MR-SAD phasing' should produce use_mr_sad=True."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input=load_fixture("mr_sad_workflow.txt"),
            expected_fn=_use_mr_sad,
        ),
    ]


# =====================================================================
# Standalone runner — useful for debugging this file directly
# =====================================================================

def main(argv=None):
    """Run all D2 scenarios against all available providers."""
    if argv is None:
        argv = sys.argv[1:]

    from framework import (
        AVAILABLE_PROVIDERS, build_summary,
        create_log_dir, run_scenario_against_providers,
        write_run_summary)
    import time

    if not AVAILABLE_PROVIDERS:
        sys.stderr.write(
            "No LLM provider API keys configured.  Set at least one of\n"
            "  GOOGLE_API_KEY, OPENAI_API_KEY, ANTHROPIC_API_KEY\n")
        return 2

    providers = list(AVAILABLE_PROVIDERS)
    scenarios = build_scenarios()
    run_one = make_directive_extraction_run_fn()
    log_dir = create_log_dir(phase="phase1-d2")

    print()
    print("LLM Decision Tests — Directive Extraction (Phase 1, D2)")
    print("=" * 65)
    print("Available providers: %s" % ", ".join(providers))
    skipped = [p for p in ("google", "openai", "anthropic")
               if p not in providers]
    if skipped:
        print("Skipped providers:   %s (no API key)" % ", ".join(skipped))
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
        providers_skipped=skipped, wall_time_s=wall_time,
        phase="phase1-d2")
    write_run_summary(log_dir, summary)

    print()
    print("Summary: %d PASS, %d FAIL, %d ERROR, %d SKIP"
          % (summary["totals"]["passed"], summary["totals"]["failed"],
             summary["totals"]["errored"], summary["totals"]["skipped"]))
    print("Total LLM calls: %d" % summary["totals"]["total_llm_calls"])
    print("Wall time:       %.1fs" % wall_time)
    print("Log:             %s" % log_dir)

    return 0 if summary["totals"]["failed"] == 0 and \
                summary["totals"]["errored"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
