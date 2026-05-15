"""
Top-level CLI runner for LLM decision tests.

Usage:
    phenix.python tests/llm/run_llm_tests.py
    phenix.python tests/llm/run_llm_tests.py --provider google
    phenix.python tests/llm/run_llm_tests.py --suite directive_extraction
    phenix.python tests/llm/run_llm_tests.py --scenario af7mjs_bug_a_regression
    phenix.python tests/llm/run_llm_tests.py --max-runs 1   # smoke test

Exit codes:
    0 = all tests passed (or all SKIP due to missing API keys)
    1 = at least one FAIL or ERROR
    2 = configuration error (no providers, bad argument)
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import sys
import time

# Make the framework importable
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import (
    AVAILABLE_PROVIDERS, build_summary,
    create_log_dir, run_scenario_against_providers,
    write_run_summary,
    make_directive_extraction_run_fn,
    make_planning_run_fn)


# Map suite name -> (scenarios-builder, run_one_fn).
# Phase 3 will add: "advice_preprocessing": { ... }
SUITES = {
    "directive_extraction": {
        "module": "tst_directive_extraction",
        "build_scenarios": None,   # filled in lazily
        "run_one_fn": None,
    },
    "planning": {
        "module": "tst_planning",
        "build_scenarios": None,
        "run_one_fn": None,
    },
}


def _load_suite(name):
    """Lazily import the suite's module and populate SUITES[name]."""
    if name not in SUITES:
        raise ValueError("Unknown suite: %s.  Known: %s"
                         % (name, ", ".join(sorted(SUITES))))
    info = SUITES[name]
    if info["build_scenarios"] is None:
        try:
            mod = __import__(info["module"], fromlist=["build_scenarios"])
        except ImportError as e:
            raise ImportError("Could not import suite '%s' (module '%s'): %s"
                              % (name, info["module"], e))
        info["build_scenarios"] = mod.build_scenarios

    if info["run_one_fn"] is None:
        if name == "directive_extraction":
            info["run_one_fn"] = make_directive_extraction_run_fn()
        elif name == "planning":
            info["run_one_fn"] = make_planning_run_fn()
    return info


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Run LLM-based decision tests for the PHENIX agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Exit codes:
  0  all PASS (or all SKIP for missing API keys)
  1  at least one FAIL or ERROR
  2  configuration error (no providers, bad argument)

Examples:
  Run everything with all configured providers:
    %(prog)s

  Smoke test with one run per scenario:
    %(prog)s --max-runs 1

  Debug a single scenario:
    %(prog)s --scenario af7mjs_bug_a_regression --verbose

  One provider only:
    %(prog)s --provider openai
""")
    parser.add_argument(
        "--suite", default=None,
        help="Comma-separated suite names.  Default: all suites.  "
             "Available suites: directive_extraction (Phase 1), "
             "planning (Phase 2).")
    parser.add_argument(
        "--provider", default=None,
        help="Single provider (google, openai, anthropic).  "
             "Default: all available.")
    parser.add_argument(
        "--scenario", default=None,
        help="Run only this scenario (debugging aid).")
    parser.add_argument(
        "--max-runs", type=int, default=None,
        help="Override max_runs for all scenarios.")
    parser.add_argument(
        "--log-dir", default=None,
        help="Where to write logs.  Default: tests/llm/logs/.")
    parser.add_argument(
        "--no-log", action="store_true",
        help="Skip writing log files (stdout only).")
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print full per-run progress to stdout.")
    args = parser.parse_args(argv)

    # Resolve providers
    if args.provider:
        if args.provider not in ("google", "openai", "anthropic"):
            sys.stderr.write("Unknown provider: %s\n" % args.provider)
            return 2
        if args.provider not in AVAILABLE_PROVIDERS:
            sys.stderr.write(
                "Provider %s requested but %s_API_KEY not set.\n"
                % (args.provider, args.provider.upper()))
            return 2
        providers = [args.provider]
    else:
        providers = list(AVAILABLE_PROVIDERS)

    if not providers:
        sys.stderr.write(
            "No LLM provider API keys configured.  Set at least one of\n"
            "  GOOGLE_API_KEY, OPENAI_API_KEY, ANTHROPIC_API_KEY\n")
        return 2

    # Resolve suites
    if args.suite:
        suite_names = [s.strip() for s in args.suite.split(",") if s.strip()]
    else:
        suite_names = list(SUITES.keys())

    unknown = [s for s in suite_names if s not in SUITES]
    if unknown:
        sys.stderr.write("Unknown suites: %s.  Known: %s\n"
                         % (", ".join(unknown), ", ".join(sorted(SUITES))))
        return 2

    # Set up log directory
    log_dir = None
    if not args.no_log:
        log_dir = create_log_dir(base_dir=args.log_dir, phase="run")

    # Header.  Distinguish two categories:
    #   skipped_no_key  — providers with no API key (truly unavailable)
    #   deselected      — providers available but not chosen via --provider
    skipped_no_key = [p for p in ("google", "openai", "anthropic")
                      if p not in AVAILABLE_PROVIDERS]
    deselected = [p for p in AVAILABLE_PROVIDERS if p not in providers]

    print()
    print("LLM Decision Tests")
    print("=" * 65)
    print("Providers selected:  %s" % ", ".join(providers))
    if deselected:
        print("Providers available but not selected (use --provider to include): "
              "%s" % ", ".join(deselected))
    if skipped_no_key:
        print("Providers unavailable: %s (no API key)"
              % ", ".join(skipped_no_key))
    print("Log directory:       %s"
          % (log_dir if log_dir else "<none, --no-log set>"))
    print()

    # Run suites
    all_verdicts = []
    found_scenario_anywhere = False
    t0 = time.time()
    for suite_name in suite_names:
        try:
            suite = _load_suite(suite_name)
        except (ValueError, ImportError) as e:
            sys.stderr.write("ERROR: %s\n" % e)
            return 2

        scenarios = suite["build_scenarios"]()

        # Filter to single scenario if requested.  If not in this suite,
        # just skip — we may find it in a later suite.  Validate at the
        # end that it was found somewhere.
        if args.scenario:
            scenarios = [s for s in scenarios if s.name == args.scenario]
            if not scenarios:
                # Don't error here; the scenario may exist in another suite.
                continue
            found_scenario_anywhere = True

        # Apply max_runs override
        if args.max_runs is not None:
            for s in scenarios:
                s.max_runs = args.max_runs

        print("Suite: %s (%d scenarios)" % (suite_name, len(scenarios)))
        print()

        run_one = suite["run_one_fn"]
        for scenario in scenarios:
            verdicts = run_scenario_against_providers(
                scenario, providers, run_one,
                log_dir=log_dir, verbose=args.verbose)
            all_verdicts.extend(verdicts)
        print()

    if args.scenario and not found_scenario_anywhere:
        sys.stderr.write("Scenario %r not found in any of the requested "
                         "suites: %s\n"
                         % (args.scenario, ", ".join(suite_names)))
        return 2

    wall_time = time.time() - t0

    # Write summary.  providers_skipped reports only "no API key" providers;
    # deselected providers are intentional user choices and not failures.
    summary = build_summary(
        all_verdicts, providers_tested=providers,
        providers_skipped=skipped_no_key, wall_time_s=wall_time,
        phase="run")
    if log_dir is not None:
        write_run_summary(log_dir, summary)

    # Final output
    totals = summary["totals"]
    print("Summary: %d PASS, %d FAIL, %d ERROR, %d SKIP"
          % (totals["passed"], totals["failed"],
             totals["errored"], totals["skipped"]))

    # Phase 2: surface individual-run fail / err counts across all
    # scenarios.  A reliability test can pass overall while still
    # exhibiting individual failures (e.g. 4/5 at threshold 0.8 = PASS
    # but 1 individual fail); the verdict alone hides that.
    indiv_failed = sum(
        len([o for o in v.outcomes if not o.passed and o.error is None])
        for v in all_verdicts)
    indiv_errored = sum(
        len([o for o in v.outcomes if not o.passed and o.error is not None])
        for v in all_verdicts)
    if indiv_failed > 0 or indiv_errored > 0:
        parts = []
        if indiv_failed > 0:
            parts.append("%d individual run%s failed"
                         % (indiv_failed, "s" if indiv_failed != 1 else ""))
        if indiv_errored > 0:
            parts.append("%d errored"
                         % indiv_errored)
        print("         (%s)" % ", ".join(parts))

    print("Total LLM calls: %d" % totals["total_llm_calls"])
    print("Wall time:       %s" % _format_duration(wall_time))
    if log_dir:
        print("Log:             %s/" % log_dir)
    print()

    # Exit code
    if totals["failed"] > 0 or totals["errored"] > 0:
        return 1
    return 0


def _format_duration(seconds):
    """Format seconds as Xm Ys or just Ys."""
    if seconds < 60:
        return "%.1fs" % seconds
    minutes = int(seconds // 60)
    secs = seconds - minutes * 60
    return "%dm %.0fs" % (minutes, secs)


if __name__ == "__main__":
    sys.exit(main())
