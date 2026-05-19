"""
C1 — LLM Decision Test: Raw-Authority Extraction (v117)

The critical reliability test for v117 Step 1.  Three variants of
the original openai bug, where the preprocessor mangled short
imperative advice into multi-action prose with 'Stop Condition: None'.

Per the AUTHORITY paragraph in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW,
the extractor LLM should defer to the raw user instruction for
intent fields when raw and processed disagree.

C1 passing reliably at threshold 0.8 is the primary green-light
that v117 Step 1 delivers value at the LLM layer (not just via
the regex backstop).  It is also the prerequisite for evaluating
Phase 2B (preprocessor deprecation) — if the LLM isn't actually
using the AUTHORITY paragraph, dropping the preprocessor wouldn't
be safe.

This file uses the framework extension `call_directive_extractor_with_raw`
and `make_directive_extraction_with_raw_run_fn` (see
framework_extensions.py).  Each scenario's `input` is a
(user_advice, raw_advice) tuple instead of a string.

Provider should be 'openai' to exercise the path that originally
failed; google/anthropic also worth covering.

Run:
    phenix.python tests/llm/tst_c1_raw_authority.py
    phenix.python tests/llm/tst_c1_raw_authority.py --provider openai
"""

from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import Scenario

# v117 added two functions to framework.py for dual-input
# extraction.  Apply v117_framework_dot_py.diff before running this test.
from framework import (call_directive_extractor_with_raw,
                       make_directive_extraction_with_raw_run_fn)


# =====================================================================
# Expected-output helpers
# =====================================================================

def _raw_authority_check(expected_program):
    """Assert that the extractor honored the RAW instruction's stop
    intent, producing after_program=<expected> and
    stop_after_requested=True — despite the processed advice's
    'Stop Condition: None'.
    """
    def check(directives):
        sc = (directives or {}).get("stop_conditions") or {}
        sar = sc.get("stop_after_requested")
        ap = sc.get("after_program")
        if sar is True and ap == expected_program:
            return (True, "raw authority honored: stop_after_requested=True, "
                          "after_program=%s" % ap)
        # Failure mode A: extractor followed the processed advice
        # ('Stop Condition: None' → no stop)
        if sar is not True:
            return (False, "extractor missed raw stop intent: "
                           "stop_after_requested=%r (the processed "
                           "advice's 'Stop Condition: None' likely won)"
                           % sar)
        # Failure mode B: extractor set flag but wrong program
        return (False, "stop_after_requested=True but after_program=%r "
                       "(expected %r)" % (ap, expected_program))
    return check


def _raw_authority_check_density_family(directives):
    """Variant for c1_density_modify_raw_authority.

    Raw='density modify and stop' is genuinely ambiguous about which
    Phenix density-modification program the user means (cryo-EM
    resolve_cryo_em, X-ray resolve, density_modification, parrot, ...).
    Without a cryo-EM/X-ray context cue in raw, the LLM may pick any
    of these reasonably.  Assert only:
      - stop_after_requested=True
      - after_program is set AND looks like a density-modification-family
        program (or omitted, falling back to the regex backstop)
    """
    sc = (directives or {}).get("stop_conditions") or {}
    sar = sc.get("stop_after_requested")
    ap = sc.get("after_program")
    if sar is not True:
        return (False, "extractor missed raw stop intent: "
                       "stop_after_requested=%r" % sar)
    # stop_after_requested is True. Now check after_program is sensible
    # OR is missing (in which case the regex backstop's after_program
    # selection takes over downstream).
    if ap is None:
        return (True, "stop_after_requested=True, after_program left "
                      "for downstream resolver")
    density_family = (
        "phenix.resolve_cryo_em", "phenix.resolve",
        "phenix.density_modification", "phenix.parrot",
        "phenix.solvent_modify",
    )
    if ap in density_family:
        return (True, "stop_after_requested=True, after_program=%s "
                      "(density-modification family)" % ap)
    return (False, "stop_after_requested=True but after_program=%r "
                   "is not a density-modification-family program" % ap)


# =====================================================================
# Scenarios
# =====================================================================

# Inputs are (user_advice, raw_advice) tuples.
#
# user_advice: what the preprocessor produced (the mangled, multi-action
#   version with 'Stop Condition: None' — the actual production output
#   from openai on the corresponding raw advice).
#
# raw_advice: the original user instruction.

_C1_1_INPUT = (
    # user_advice (preprocessor output, openai-mangled)
    "1. Input Files Found: model.pdb, data.mtz, map.ccp4\n"
    "\n"
    "2. Experiment Type: cryo-EM\n"
    "\n"
    "3. Primary Goal: Perform density modification, then build and "
    "refine the model.\n"
    "\n"
    "4. Key Parameters:\n"
    "   - Resolution limit: 3.0\n"
    "\n"
    "5. Program Parameters: None\n"
    "\n"
    "6. Special Instructions: None\n"
    "\n"
    "7. Stop Condition: None",
    # raw_advice (user's actual instruction)
    "density modify and stop",
)

_C1_2_INPUT = (
    "1. Input Files Found: model.pdb, data.mtz\n"
    "\n"
    "2. Experiment Type: X-ray crystallography\n"
    "\n"
    "3. Primary Goal: Refine the model and then validate using "
    "MolProbity to assess geometry.\n"
    "\n"
    "4. Key Parameters:\n"
    "   - Resolution limit: 2.0\n"
    "\n"
    "5. Program Parameters: None\n"
    "\n"
    "6. Special Instructions: None\n"
    "\n"
    "7. Stop Condition: None",
    "refine and stop",
)

_C1_3_INPUT = (
    "1. Input Files Found: data.mtz\n"
    "\n"
    "2. Experiment Type: X-ray crystallography\n"
    "\n"
    "3. Primary Goal: Analyze the data with xtriage to check for "
    "twinning, then refine the model.\n"
    "\n"
    "4. Key Parameters:\n"
    "   - Resolution limit: 2.5\n"
    "\n"
    "5. Program Parameters: None\n"
    "\n"
    "6. Special Instructions: None\n"
    "\n"
    "7. Stop Condition: None",
    "run xtriage and stop",
)


def build_scenarios():
    """C1 — three raw-authority reliability scenarios."""
    return [
        Scenario(
            name="c1_density_modify_raw_authority",
            description=(
                "Raw='density modify and stop', processed expanded to "
                "multi-action with 'Stop Condition: None'.  The "
                "AUTHORITY paragraph must cause the extractor to "
                "honor raw's stop intent.  After_program assertion is "
                "permissive — any density-modification-family program "
                "(resolve_cryo_em / resolve / density_modification / "
                "parrot) is accepted, since the raw advice doesn't "
                "specify cryo-EM vs X-ray context."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_C1_1_INPUT,
            expected_fn=_raw_authority_check_density_family,
        ),
        Scenario(
            name="c1_refine_raw_authority",
            description=(
                "Raw='refine and stop', processed expanded to "
                "multi-action with 'Stop Condition: None'.  Same "
                "raw-authority pattern as c1_density_modify."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_C1_2_INPUT,
            expected_fn=_raw_authority_check("phenix.refine"),
        ),
        Scenario(
            name="c1_xtriage_raw_authority",
            description=(
                "Raw='run xtriage and stop', processed expanded to "
                "multi-action with 'Stop Condition: None'.  Same "
                "raw-authority pattern."),
            decision_point="directive_extraction",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_C1_3_INPUT,
            expected_fn=_raw_authority_check("phenix.xtriage"),
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
    run_one = make_directive_extraction_with_raw_run_fn()
    log_dir = create_log_dir(phase="phase1-c1")

    print()
    print("LLM Decision Tests — C1 Raw Authority (v117)")
    print("=" * 65)
    print("Available providers: %s" % ", ".join(providers))
    print("Log directory:       %s" % log_dir)
    if "openai" not in providers:
        print("WARN: openai not in providers — C1 was designed for the "
              "openai-specific bug; results from other providers are "
              "informative but not the primary signal.")
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
        phase="phase1-c1")
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
