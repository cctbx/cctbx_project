"""
Framework for LLM-based decision tests of the PHENIX agent.

Provides:
- Scenario / RunOutcome / Verdict dataclasses
- Provider detection (AVAILABLE_PROVIDERS) from API key env vars
- Early-termination runner that runs each scenario up to max_runs times
  but stops as soon as the verdict is certain (capability: first pass;
  reliability: first time the threshold is provably met or unreachable)
- Logging helpers (timestamped log directory, per-scenario detail logs,
  run_summary.json)

This framework is callable directly from individual test modules
(tst_*.py) and from the top-level CLI runner (run_llm_tests.py).
It does NOT depend on the PHENIX agent's deterministic test machinery
(run_tests_with_fail_fast etc.); the deterministic suite and the LLM
suite are kept separate by design.

See PHASE1_PLAN.md §3.1 for the design rationale.
"""

from __future__ import absolute_import, division, print_function

import datetime
import json
import math
import os
import sys
import time
import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional


# ----------------------------------------------------------------------
# Provider detection
# ----------------------------------------------------------------------

# Tuple, not list, because this should not be mutated after import.
_KNOWN_PROVIDERS = ("google", "openai", "anthropic")


def _detect_available_providers():
    """Return list of providers whose API key env var is set.

    Side-effect-free; called once at module import to populate
    AVAILABLE_PROVIDERS.
    """
    out = []
    for provider in _KNOWN_PROVIDERS:
        env_var = "%s_API_KEY" % provider.upper()
        if os.environ.get(env_var):
            out.append(provider)
    return out


AVAILABLE_PROVIDERS = _detect_available_providers()


# ----------------------------------------------------------------------
# Core dataclasses
# ----------------------------------------------------------------------

@dataclass
class Scenario:
    """A single test scenario.

    Attributes
    ----------
    name : str
        Short identifier (used as log filename component).
    description : str
        One-line human-readable summary.
    decision_point : str
        Which agent decision point is being tested
        ("directive_extraction" for Phase 1).
    test_type : str
        "capability" (≥1 pass = success) or "reliability"
        (≥threshold × max_runs passes = success).
    threshold : float
        Used by reliability tests only.  Default 0.8.
    max_runs : int
        Hard cap on number of LLM calls per (scenario, provider).
        Early termination usually stops well before this.
    input : Any
        Scenario-specific input.  For D2 this is the user_advice string
        passed to extract_directives.
    expected_fn : Callable
        Function (parsed_output) -> (passed: bool, why: str).
        Called once per LLM call; the "passed" boolean drives the
        verdict logic.
    """
    name: str
    description: str
    decision_point: str
    test_type: str
    threshold: float = 0.8
    max_runs: int = 5
    input: Any = None
    expected_fn: Callable = None


@dataclass
class RunOutcome:
    """Result of a single LLM call within a scenario."""
    run_index: int                  # 0-based
    elapsed_s: float
    passed: bool                    # from expected_fn; False on error
    why: str                        # explanation from expected_fn or error msg
    raw_output: str                 # full LLM response text (or "" on error)
    parsed: Any                     # parsed dict, or None
    error: Optional[str]            # None on success; populated on exception


@dataclass
class Verdict:
    """Aggregated result across all runs of one (scenario, provider)."""
    scenario_name: str
    provider: str
    result: str                     # "PASS", "FAIL", "ERROR", "SKIP"
    n_runs_executed: int
    n_passed: int
    n_failed: int
    n_errored: int
    early_stop_reason: Optional[str]
    outcomes: List[RunOutcome] = field(default_factory=list)
    avg_elapsed_s: float = 0.0
    test_type: str = ""             # mirror for log output
    threshold: float = 0.0


# ----------------------------------------------------------------------
# Pass / fail thresholds
# ----------------------------------------------------------------------

def _required_passes(scenario):
    """Number of passes needed for a reliability test to succeed."""
    return int(math.ceil(scenario.threshold * scenario.max_runs))


# ----------------------------------------------------------------------
# Early-termination runner
# ----------------------------------------------------------------------

def run_with_early_termination(scenario, provider, run_one_fn,
                               progress_fn=None):
    """Run scenario up to max_runs times, stopping early when the
    verdict is certain.

    Verdict logic:

      capability test:
        PASS at first successful run.
        FAIL if 0 successes after max_runs (or ERROR if every run errored).

      reliability test:
        PASS when n_passed >= ceil(threshold * max_runs).
        FAIL when n_passed + (remaining_runs + n_errored) < required passes.

    ERROR vs FAIL handling:
        A "passed=False" outcome counts as either a semantic FAILURE
        (LLM gave wrong answer; n_failed increments) or an ERROR
        (exception, parse failure, etc.; n_errored increments).

        For early-fail unreachability, errors are treated as
        INDETERMINATE — they don't count toward the failure budget.
        Rationale: a single transient API error shouldn't cause us
        to give up on the entire scenario.  In practice this means
        we keep running until either:
          (a) we accumulate enough semantic failures that even if
              all remaining runs pass we couldn't reach threshold;
          (b) we hit max_runs.

        Equivalently: the "best case" for the remaining budget is
        n_passed + (remaining_runs + n_errored_so_far).  If even that
        can't reach the required count, we stop.  This is unchanged
        from a no-errors-yet world where remaining_runs IS the
        best-case budget.

    Final verdict result:
        - n_passed >= required (capability: ≥1; reliability: ≥ceil)
          → PASS
        - n_runs == n_errored (every run errored, no semantic info)
          → ERROR
        - otherwise → FAIL

    Parameters
    ----------
    scenario : Scenario
    provider : str
    run_one_fn : Callable[[Scenario, str, int], RunOutcome]
        Function that performs one LLM call and returns its outcome.
        The int argument is the run_index.
    progress_fn : Optional[Callable[[RunOutcome], None]]
        Called after each run with the new outcome.  Useful for live
        progress printing.

    Returns
    -------
    Verdict
    """
    outcomes = []
    n_pass = 0
    n_fail = 0       # passed=False, error is None (semantic mismatch)
    n_err = 0        # passed=False, error is not None
    early_stop_reason = None

    if scenario.max_runs < 1:
        # Degenerate: nothing to run.  Return a synthetic FAIL.
        return Verdict(
            scenario_name=scenario.name,
            provider=provider,
            result="FAIL",
            n_runs_executed=0,
            n_passed=0,
            n_failed=0,
            n_errored=0,
            early_stop_reason="max_runs=%d < 1" % scenario.max_runs,
            outcomes=[],
            avg_elapsed_s=0.0,
            test_type=scenario.test_type,
            threshold=scenario.threshold if scenario.test_type == "reliability"
                                           else 0.0,
        )

    required_for_pass = _required_passes(scenario)

    for i in range(scenario.max_runs):
        outcome = run_one_fn(scenario, provider, i)
        outcomes.append(outcome)
        if progress_fn:
            progress_fn(outcome)

        if outcome.passed:
            n_pass += 1
        else:
            if outcome.error is not None:
                n_err += 1
            else:
                n_fail += 1

        # Check early-termination conditions

        if scenario.test_type == "capability":
            # First pass clinches it.
            if n_pass >= 1:
                early_stop_reason = (
                    "capability satisfied at run %d" % (i + 1))
                break
            # No early-fail for capability — keep trying until max_runs.

        elif scenario.test_type == "reliability":
            remaining = scenario.max_runs - (i + 1)
            # Early pass: reached required count
            if n_pass >= required_for_pass:
                early_stop_reason = (
                    "reliability satisfied: %d/%d passes (>= %d needed)"
                    % (n_pass, i + 1, required_for_pass))
                break
            # Early fail: even if all remaining runs pass AND every error
            # so far had really been a pass, we still couldn't reach
            # the threshold.  Errors are treated as indeterminate
            # (not as failures) for this check, so a few transient API
            # errors don't trigger premature give-up.
            best_case_total = n_pass + remaining + n_err
            if best_case_total < required_for_pass:
                early_stop_reason = (
                    "reliability unreachable: %d passes, %d semantic "
                    "failures after %d runs; need %d total but best-case "
                    "is %d"
                    % (n_pass, n_fail, i + 1, required_for_pass,
                       best_case_total))
                break

        else:
            raise ValueError(
                "Unknown test_type %r for scenario %s"
                % (scenario.test_type, scenario.name))

    # Determine final verdict
    n_runs = len(outcomes)
    if scenario.test_type == "capability":
        if n_pass >= 1:
            result = "PASS"
        elif n_runs > 0 and n_err == n_runs:
            result = "ERROR"
        else:
            result = "FAIL"
    else:  # reliability
        if n_pass >= required_for_pass:
            result = "PASS"
        elif n_runs > 0 and n_err == n_runs:
            result = "ERROR"
        else:
            result = "FAIL"

    avg_elapsed = (
        sum(o.elapsed_s for o in outcomes) / n_runs if n_runs > 0 else 0.0)

    return Verdict(
        scenario_name=scenario.name,
        provider=provider,
        result=result,
        n_runs_executed=n_runs,
        n_passed=n_pass,
        n_failed=n_fail,
        n_errored=n_err,
        early_stop_reason=early_stop_reason,
        outcomes=outcomes,
        avg_elapsed_s=avg_elapsed,
        test_type=scenario.test_type,
        threshold=scenario.threshold if scenario.test_type == "reliability"
                                       else 0.0,
    )


# ----------------------------------------------------------------------
# Production-faithful LLM invocation for D2 (directive extraction)
# ----------------------------------------------------------------------

def call_directive_extractor(user_advice, provider):
    """Call extract_directives() in a production-faithful way.

    Returns
    -------
    (raw_output, parsed_directives, error_msg, captured_log)

    raw_output : str
        The raw text the framework reconstructs from the captured log
        (the agent's _call_llm logs the response indirectly via
         log_func messages; we capture and concatenate them).
        Empty string on error.
    parsed_directives : dict or None
        The directives dict returned by extract_directives, or None
        on error.
    error_msg : str or None
        Exception message if extract_directives raised, else None.
    captured_log : list of str
        Every log line emitted by extract_directives during the call.
        Useful for diagnosis when a run fails.

    Note: extract_directives does NOT return the raw LLM response
    directly; it returns the parsed/validated directives dict.  To get
    the raw response text we'd need to monkey-patch _call_llm, which is
    more invasive than needed for Phase 1.  Instead we capture the log
    output (which includes some response context) and report the parsed
    directives as the primary signal.  Tests assert on the parsed
    directives, which is what production code consumes anyway.
    """
    try:
        # Use the production import path; fall back to direct path for
        # standalone test invocation.
        try:
            from libtbx.langchain.agent.directive_extractor import (
                extract_directives)
        except ImportError:
            from agent.directive_extractor import extract_directives
    except ImportError as e:
        return ("", None, "ImportError: %s" % e, [])

    captured = []

    def log_fn(msg):
        captured.append(msg)

    try:
        directives = extract_directives(
            user_advice, provider=provider, log_func=log_fn)
        # Convert the captured log into a "raw_output" proxy for the
        # detail log.  The directives dict is the primary signal.
        raw_output = "\n".join(captured)
        return (raw_output, directives, None, captured)
    except Exception as e:
        tb = traceback.format_exc()
        # Include short traceback in error message for diagnosis
        return ("\n".join(captured), None,
                "%s: %s\n%s" % (type(e).__name__, e, tb[-500:]),
                captured)


def make_directive_extraction_run_fn(call_fn=call_directive_extractor):
    """Build a run_one_fn for directive-extraction scenarios.

    call_fn is injectable for testing the framework itself.  Default
    is the production-faithful caller.
    """

    def run_one(scenario, provider, run_index):
        t0 = time.time()
        raw_output, parsed, error, _captured = call_fn(
            scenario.input, provider)
        elapsed = time.time() - t0

        if error is not None:
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="LLM call errored: %s" % error.split("\n")[0],
                raw_output=raw_output,
                parsed=None,
                error=error,
            )

        if parsed is None:
            # Shouldn't happen given the contract of extract_directives
            # (it returns {} on empty advice, not None) but defend.
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="extract_directives returned None",
                raw_output=raw_output,
                parsed=None,
                error="parsed is None",
            )

        try:
            passed, why = scenario.expected_fn(parsed)
        except Exception as e:
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="expected_fn raised: %s" % e,
                raw_output=raw_output,
                parsed=parsed,
                error="expected_fn raised: %s: %s" % (type(e).__name__, e),
            )

        return RunOutcome(
            run_index=run_index,
            elapsed_s=elapsed,
            passed=bool(passed),
            why=str(why),
            raw_output=raw_output,
            parsed=parsed,
            error=None,
        )

    return run_one


# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------

def create_log_dir(base_dir=None, phase="phase1"):
    """Create a timestamped log directory.  Returns its path.

    Structure:
        <base_dir>/<timestamp>-<phase>/
            run_summary.json
            details/
    """
    if base_dir is None:
        # Default to a "logs" directory next to this framework.py file,
        # rather than resolving relative to the current working
        # directory.  This means logs always land in the same place
        # regardless of where the user invokes the test runner from.
        # Override with LLM_TEST_LOG_DIR or pass base_dir explicitly.
        base_dir = os.environ.get(
            "LLM_TEST_LOG_DIR",
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "logs"))
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    log_dir = os.path.join(base_dir, "%s-%s" % (timestamp, phase))
    os.makedirs(os.path.join(log_dir, "details"), exist_ok=True)
    return log_dir


def write_scenario_detail(log_dir, verdict):
    """Write per-scenario per-provider plaintext detail log."""
    if log_dir is None:
        return
    fname = "%s__%s.log" % (verdict.scenario_name, verdict.provider)
    path = os.path.join(log_dir, "details", fname)
    with open(path, "w") as f:
        f.write("=== Scenario: %s ===\n" % verdict.scenario_name)
        f.write("Provider: %s\n" % verdict.provider)
        f.write("Test type: %s" % verdict.test_type)
        if verdict.test_type == "reliability":
            f.write(", threshold=%.2f" % verdict.threshold)
        f.write("\n")
        f.write("Result: %s\n" % verdict.result)
        f.write("Runs: %d (passed=%d, failed=%d, errored=%d)\n"
                % (verdict.n_runs_executed, verdict.n_passed,
                   verdict.n_failed, verdict.n_errored))
        f.write("Avg elapsed: %.2fs\n" % verdict.avg_elapsed_s)
        if verdict.early_stop_reason:
            f.write("Early stop: %s\n" % verdict.early_stop_reason)
        f.write("\n")

        for outcome in verdict.outcomes:
            f.write("--- Run %d ---\n" % (outcome.run_index + 1))
            f.write("Elapsed: %.2fs\n" % outcome.elapsed_s)
            f.write("Passed:  %s\n" % outcome.passed)
            f.write("Why:     %s\n" % outcome.why)
            if outcome.error:
                # Full error (may include traceback up to 500 chars from
                # call_directive_extractor); preserve it for diagnosis.
                f.write("Error:\n")
                for line in outcome.error.splitlines():
                    f.write("  %s\n" % line)
            if outcome.parsed is not None:
                f.write("Parsed:  %s\n"
                        % json.dumps(outcome.parsed,
                                     indent=2,
                                     default=str)[:2000])
            if outcome.raw_output:
                preview = outcome.raw_output[:500]
                if len(outcome.raw_output) > 500:
                    preview += "...[truncated]"
                f.write("Raw output (first 500 chars):\n%s\n" % preview)
            f.write("\n")


def write_run_summary(log_dir, summary):
    """Write run_summary.json."""
    if log_dir is None:
        return
    path = os.path.join(log_dir, "run_summary.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=2, default=str)


def build_summary(verdicts, providers_tested, providers_skipped,
                  wall_time_s, phase="phase1"):
    """Build the run_summary.json structure from a list of Verdicts."""
    # Group verdicts by scenario name
    by_scenario = {}
    for v in verdicts:
        by_scenario.setdefault(v.scenario_name, []).append(v)

    scenarios_out = []
    for name in sorted(by_scenario):
        vs = by_scenario[name]
        # Pull test_type / threshold from first verdict; they're all same.
        first = vs[0]
        scenarios_out.append({
            "name": name,
            "test_type": first.test_type,
            "threshold": first.threshold,
            "results": [
                {
                    "provider": v.provider,
                    "verdict": v.result,
                    "n_passed": v.n_passed,
                    "n_failed": v.n_failed,
                    "n_errored": v.n_errored,
                    "n_runs_executed": v.n_runs_executed,
                    "avg_elapsed_s": v.avg_elapsed_s,
                    "early_stop_reason": v.early_stop_reason,
                }
                for v in vs
            ],
        })

    totals = {
        "passed": sum(1 for v in verdicts if v.result == "PASS"),
        "failed": sum(1 for v in verdicts if v.result == "FAIL"),
        "errored": sum(1 for v in verdicts if v.result == "ERROR"),
        "skipped": sum(1 for v in verdicts if v.result == "SKIP"),
        "total_llm_calls": sum(v.n_runs_executed for v in verdicts),
        "wall_time_s": wall_time_s,
    }

    return {
        "phase": phase,
        "timestamp": datetime.datetime.now().isoformat(),
        "providers_tested": providers_tested,
        "providers_skipped": providers_skipped,
        "scenarios": scenarios_out,
        "totals": totals,
    }


# ----------------------------------------------------------------------
# Fixture loading
# ----------------------------------------------------------------------

# Resolve fixture directory relative to this file.
_FIXTURE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "fixtures")


def load_fixture(name):
    """Read a fixture file by name (no path components).

    Looks in tests/llm/fixtures/.  Returns the file's text content.
    """
    if "/" in name or "\\" in name:
        raise ValueError("fixture name must not contain path components: %r"
                         % name)
    path = os.path.join(_FIXTURE_DIR, name)
    if not os.path.exists(path):
        raise IOError("fixture not found: %s" % path)
    with open(path, "r") as f:
        return f.read()


def load_state_fixture(path):
    """Load a captured workflow-state fixture.

    Phase 1 doesn't use this (D2 fixtures are plain text), but Phase 2
    will.  Provided here so Phase 2 inherits the helper.

    Returns (state_inputs, expected, metadata) as dicts.
    """
    with open(path, "r") as f:
        data = json.load(f)
    return (
        data.get("state_inputs", {}),
        data.get("expected", {}),
        data.get("metadata", {}),
    )


# ----------------------------------------------------------------------
# Top-level scenario runner
# ----------------------------------------------------------------------

def run_scenario_against_providers(scenario, providers, run_one_fn,
                                   log_dir=None, verbose=False,
                                   progress_print=True):
    """Run one scenario against every provider in the list.

    Returns a list of Verdicts, one per provider.
    """
    verdicts = []
    for provider in providers:
        if provider not in AVAILABLE_PROVIDERS:
            verdicts.append(Verdict(
                scenario_name=scenario.name,
                provider=provider,
                result="SKIP",
                n_runs_executed=0,
                n_passed=0,
                n_failed=0,
                n_errored=0,
                early_stop_reason="%s_API_KEY not set" % provider.upper(),
                outcomes=[],
                avg_elapsed_s=0.0,
                test_type=scenario.test_type,
                threshold=(scenario.threshold
                           if scenario.test_type == "reliability" else 0.0),
            ))
            if progress_print:
                print("  [%s]   %-40s  SKIP  (no API key)"
                      % (provider, scenario.name))
            continue

        def _progress(outcome):
            if verbose:
                status = "PASS" if outcome.passed else "FAIL"
                if outcome.error:
                    status = "ERROR"
                print("       run %d  %s  (%.1fs)  %s"
                      % (outcome.run_index + 1, status, outcome.elapsed_s,
                         outcome.why))

        verdict = run_with_early_termination(
            scenario, provider, run_one_fn, progress_fn=_progress)
        verdicts.append(verdict)

        if log_dir is not None:
            write_scenario_detail(log_dir, verdict)

        if progress_print:
            early = ""
            if verdict.early_stop_reason and verdict.result in ("PASS",
                                                                 "FAIL"):
                early = "  *early stop*"
            tag = "%d/%d" % (verdict.n_passed, verdict.n_runs_executed)
            print("  [%s]   %-40s  %s  %-5s (%.1fs)%s"
                  % (provider, scenario.name, verdict.result, tag,
                     verdict.avg_elapsed_s, early))

    return verdicts


# ----------------------------------------------------------------------
# v117 additions — dual-input directive extraction (raw + processed)
# ----------------------------------------------------------------------
#
# v117 Step 1 added a `raw_advice` parameter to extract_directives().
# When raw_advice differs from user_advice, the extractor uses
# DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW (dual-input + AUTHORITY paragraph).
# Scenarios that exercise this path use a (user_advice, raw_advice)
# tuple as input rather than a string, and need the run_one_fn below.
#
# Used by tst_c1_raw_authority.py.  The original single-input path
# (call_directive_extractor + make_directive_extraction_run_fn above)
# remains unchanged.

def call_directive_extractor_with_raw(input_tuple, provider):
    """Call extract_directives() with the v117 raw_advice parameter.

    Production path (v117+):
        extract_directives(user_advice, provider=..., raw_advice=...)

    When raw_advice differs from user_advice, the dual-input prompt
    (DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW) is selected and the LLM
    sees both inputs plus the AUTHORITY paragraph stating raw is
    authoritative for intent.

    Parameters
    ----------
    input_tuple : tuple of (user_advice, raw_advice)
        user_advice : str
            The "processed" advice (typically what the preprocessor
            produces; in tests, the post-preprocessing text).
        raw_advice : str or None
            The original raw user instruction.  If None or equal to
            user_advice, the extractor falls back to the single-input
            prompt (backward compat).

    Returns
    -------
    (raw_output, parsed_directives, error_msg, captured_log)
        Same shape as `call_directive_extractor`.
    """
    user_advice, raw_advice = input_tuple

    try:
        try:
            from libtbx.langchain.agent.directive_extractor import (
                extract_directives)
        except ImportError:
            from agent.directive_extractor import extract_directives
    except ImportError as e:
        return ("", None, "ImportError: %s" % e, [])

    captured = []

    def log_fn(msg):
        captured.append(msg)

    try:
        directives = extract_directives(
            user_advice,
            provider=provider,
            log_func=log_fn,
            raw_advice=raw_advice,
        )
        raw_output = "\n".join(captured)
        return (raw_output, directives, None, captured)
    except Exception as e:
        tb = traceback.format_exc()
        return ("\n".join(captured), None,
                "%s: %s\n%s" % (type(e).__name__, e, tb[-500:]),
                captured)


def make_directive_extraction_with_raw_run_fn(
        call_fn=call_directive_extractor_with_raw):
    """Build a run_one_fn for dual-input directive-extraction scenarios.

    Mirrors `make_directive_extraction_run_fn`, but scenario.input is
    a (user_advice, raw_advice) tuple instead of a string.
    """
    def run_one(scenario, provider, run_index):
        t0 = time.time()
        raw_output, parsed, error, _captured = call_fn(
            scenario.input, provider)
        elapsed = time.time() - t0

        if error is not None:
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="LLM call errored: %s" % error.split("\n")[0],
                raw_output=raw_output,
                parsed=None,
                error=error,
            )

        if parsed is None:
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="extract_directives returned None",
                raw_output=raw_output,
                parsed=None,
                error="parsed is None",
            )

        try:
            passed, why = scenario.expected_fn(parsed)
        except Exception as e:
            return RunOutcome(
                run_index=run_index,
                elapsed_s=elapsed,
                passed=False,
                why="expected_fn raised: %s" % e,
                raw_output=raw_output,
                parsed=parsed,
                error="expected_fn raised: %s: %s" % (type(e).__name__, e),
            )

        return RunOutcome(
            run_index=run_index,
            elapsed_s=elapsed,
            passed=passed,
            why=why,
            raw_output=raw_output,
            parsed=parsed,
            error=None,
        )

    return run_one


# ----------------------------------------------------------------------
# v117.2 stubs — placeholders for names referenced by other modules
# ----------------------------------------------------------------------
#
# These stubs satisfy imports without changing behavior, letting the
# framework load cleanly even when callers reference functions that
# haven't yet been implemented.  Each stub documents what a real
# implementation should do; replace in place when implementing.
#
def validate_planning_state(state):
    """Stub: validate a planning-suite scenario's state dict.

    Until a real implementation is added, this is a no-op — every state
    is accepted.  A real implementation should raise ValueError if the
    state dict is missing required keys or has malformed values (per
    PHASE2_PLAN_v2.md §4 required-fields checklist).

    Called from build_scenarios() in tst_planning.py to catch
    malformed scenarios early.  Until Phase 2A defines the planning
    state schema, scenarios are accepted as-is.
    """
    return None


def is_stop_intent(text):
    """Stub: detect whether a string expresses user-stop intent.

    Until a real implementation is added, returns False unconditionally.
    This is conservative — no text is classified as having stop intent,
    which means planning scenarios depending on stop-intent semantics
    will see a 'no stop' world by default.

    Replace with a real implementation when Phase 2A planning lands.
    The agent's directive extractor has a related helper
    (_is_stop_after_requested) that detects 11 positive patterns in raw
    user advice; a real implementation here could import that helper
    directly, or duplicate the patterns if the planning suite needs a
    framework-local copy.
    """
    return False


def is_api_availability_error(exception):
    """Stub: classifier for LLM provider API availability errors.

    Returns True if the exception represents a transient provider
    issue (rate limit, timeout, 503, etc.) that should be retried or
    cause the test to SKIP rather than FAIL.  Returns False otherwise.

    Until a real implementation is added, returns False unconditionally
    — exceptions are not classified as API-availability errors, so
    the original error propagates normally (same as before this stub
    existed).  This preserves whatever default error handling the
    caller had.

    Replace with a real implementation that checks for provider-
    specific error types or HTTP status codes.
    """
    return False


# tst_b2_gated_wipe_planning.py and some downstream runners (e.g.
# run_llm_tests.py in trees that anticipate Phase 2A) import these
# names from framework.py.  Until Phase 2A merges them as real
# implementations (per PHASE2_PLAN_v2.md §2 and §6), these stubs
# satisfy the import without behavioral change for Phase 1 tests.
#
# Calling either stub raises NotImplementedError with a clear pointer
# to PHASE2_PLAN_v2.md.  When Phase 2A lands, replace these stubs in
# place with the real functions; no caller changes required.


def call_planning_llm(*args, **kwargs):
    """Phase 2A stub.  Replace with real implementation per
    PHASE2_PLAN_v2.md §2 ('Production-faithful LLM invocation') when
    Phase 2A is merged."""
    raise NotImplementedError(
        "call_planning_llm is a Phase 2A stub; Phase 2A "
        "(planning suite) is not yet merged.  See "
        "PHASE2_PLAN_v2.md for the implementation plan.  "
        "Phase 1 suites (directive_extraction, a1_smoke, "
        "b1_stop_after_requested, c1_raw_authority) do not "
        "require this function and run normally without it.")


def make_planning_run_fn(call_fn=None):
    """Phase 2A stub.  Replace with real implementation per
    PHASE2_PLAN_v2.md §2 + §6 when Phase 2A is merged.

    Returns a run_one callable that produces an ERRORED RunOutcome
    instead of raising, so the planning suite loads cleanly and each
    scenario reports an error rather than crashing the entire runner.
    Other suites in the same invocation (e.g. directive_extraction)
    complete normally.
    """
    _msg = ("make_planning_run_fn is a Phase 2A stub; Phase 2A "
            "(planning suite) is not yet merged.  See "
            "PHASE2_PLAN_v2.md for the implementation plan.  "
            "Other suites (directive_extraction, a1_smoke, "
            "b1_stop_after_requested, c1_raw_authority) do not "
            "require this function and run normally without it.")

    def _stub_run_one(scenario, provider, run_index):
        return RunOutcome(
            run_index=run_index,
            elapsed_s=0.0,
            passed=False,
            why=("Phase 2A planning infrastructure not yet merged "
                 "(stub in framework.py)"),
            raw_output="",
            parsed=None,
            error=_msg,
        )
    return _stub_run_one


# ----------------------------------------------------------------------
# Sanity-check entry point (run framework.py standalone)
# ----------------------------------------------------------------------

if __name__ == "__main__":
    print("LLM test framework module — Phase 1")
    print("Available providers: %s" % (AVAILABLE_PROVIDERS or "<none>"))
    print("Fixture dir:         %s" % _FIXTURE_DIR)
    if not AVAILABLE_PROVIDERS:
        print()
        print("No API keys configured.  Set at least one of:")
        for p in _KNOWN_PROVIDERS:
            print("  %s_API_KEY" % p.upper())
        sys.exit(2)
    sys.exit(0)
