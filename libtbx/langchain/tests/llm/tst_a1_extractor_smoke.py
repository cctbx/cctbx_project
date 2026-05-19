"""
A1 — LLM Decision Test: Extractor Smoke Test (v117)

Three short advice strings, each exercising a different code path
through the extractor.  Each is a CAPABILITY test (one pass over
max_runs = success), threshold not applicable.

This is the smoke test that confirms the v117 Step A brace-fix
hasn't broken the LLM extraction path for standalone callers.
Production already used .replace() to dodge the brace bug (v116.20
defensive code), but tests and debug scripts exercise .format() so
A1 confirms both paths produce usable output.

Run from a parent directory with framework.py importable:

    cd $PHENIX/modules/cctbx_project/libtbx/langchain
    phenix.python tests/llm/tst_a1_extractor_smoke.py

Or via the top-level CLI runner (after registering the suite in
run_llm_tests.py — see suite_registration.py.diff for the patch):

    phenix.python tests/llm/run_llm_tests.py --suite a1_smoke
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

def _produces_usable_dict(directives):
    """The extractor returned a dict with at least one section.

    'Usable' is loose: we only assert the return value is a dict and
    contains at least one expected top-level key.  Specific intent
    assertions are covered by B1 and C1.
    """
    if not isinstance(directives, dict):
        return (False, "expected dict, got %r" % type(directives).__name__)
    expected_keys = (
        "program_settings", "stop_conditions",
        "workflow_preferences", "file_preferences", "intent",
    )
    present = [k for k in expected_keys if k in directives]
    if not present:
        return (False, "no expected top-level keys (got: %r)"
                       % list(directives.keys()))
    return (True, "directives keys: %s" % ", ".join(present))


def _has_density_modify_or_resolve_program(directives):
    """A1.1: 'density modify and stop' — expect a directive that
    indicates the LLM understood density modification as a routing
    intent.  Multiple acceptable signals:
      - after_program is a density-modification-family program
      - program_settings mentions resolve_cryo_em or similar
      - stop_after_requested=True (the LLM made a stop assertion)
    Any one signal is sufficient.
    """
    if not isinstance(directives, dict):
        return (False, "not a dict")

    sc = directives.get("stop_conditions") or {}
    ap = sc.get("after_program")
    sar = sc.get("stop_after_requested")

    density_family = (
        "phenix.resolve_cryo_em", "phenix.resolve",
        "phenix.density_modification", "phenix.parrot",
        "phenix.solvent_modify",
    )
    if ap in density_family:
        return (True, "after_program=%s (density-modification family)" % ap)

    if sar is True:
        return (True, "stop_after_requested=True (LLM stop assertion)")

    flat = repr(directives).lower()
    if "resolve_cryo_em" in flat or "resolve cryo" in flat:
        return (True, "directives mention resolve_cryo_em")
    if "density" in flat and "stop" in flat:
        return (True, "directives mention density and stop")
    return (False,
            "no density/resolve/stop signal in directives: %r" % directives)


def _has_resolution_2_8(directives):
    """A1.3: 'use resolution 2.8 in autosol' — expect resolution
    setting somewhere in program_settings.
    """
    if not isinstance(directives, dict):
        return (False, "not a dict")
    ps = directives.get("program_settings") or {}
    # Look anywhere in program_settings (default scope or per-program)
    for scope, params in ps.items():
        if isinstance(params, dict):
            res = params.get("resolution")
            if res is not None:
                # Accept 2.8 or close; LLMs sometimes round
                try:
                    if abs(float(res) - 2.8) < 0.05:
                        return (True, "resolution=%s under %s" % (res, scope))
                except (TypeError, ValueError):
                    pass
    return (False, "no resolution=2.8 found in program_settings: %r" % ps)


# =====================================================================
# Scenarios
# =====================================================================

def build_scenarios():
    """A1 scenarios — three smoke-test inputs."""
    return [
        Scenario(
            name="a1_smoke_density_modify_and_stop",
            description=(
                "Short imperative: 'density modify and stop'.  "
                "Extractor must return a usable dict with some "
                "indication of density modification routing."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input="density modify and stop",
            expected_fn=_has_density_modify_or_resolve_program,
        ),
        Scenario(
            name="a1_smoke_ligand_fitting",
            description=(
                "Short imperative: 'fit the atp ligand into the map'.  "
                "Extractor must return a usable dict."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input="Fit the ATP ligand into the map.",
            expected_fn=_produces_usable_dict,
        ),
        Scenario(
            name="a1_smoke_resolution_setting",
            description=(
                "Numeric parameter: 'use resolution 2.8 in autosol'.  "
                "Extractor must produce resolution=2.8 under "
                "program_settings."),
            decision_point="directive_extraction",
            test_type="capability",
            max_runs=5,
            input="Use resolution 2.8 in autosol.",
            expected_fn=_has_resolution_2_8,
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
    log_dir = create_log_dir(phase="phase1-a1-smoke")

    print()
    print("LLM Decision Tests — A1 Extractor Smoke (v117)")
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
        phase="phase1-a1-smoke")
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
