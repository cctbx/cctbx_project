# LLM Decision Tests

Tests that exercise specific LLM-driven decisions in the PHENIX agent.
Separate from the deterministic test suite because they:

- require API keys (Google / OpenAI / Anthropic)
- are slow (~4-8 min for a full Phase 1+2 run with two providers)
- cost real money (~$0.10-0.50 per run depending on providers)
- are non-deterministic (LLM outputs vary; we test pass rates with
  early termination)

## What's tested

**Phase 1: Directive Extraction (D2)**

Six scenarios exercising the directive extractor's LLM call.  Tests
whether structured directives (after_program, skip_programs,
wants_validation_only, etc.) are correctly extracted from user
advice.

| Scenario | Type | Threshold |
|----------|------|-----------|
| `af7mjs_bug_a_regression` | reliability | 0.9 |
| `explicit_stop_after_phaser` | capability | — |
| `resolution_and_space_group` | capability | — |
| `skip_programs` | capability | — |
| `validation_only` | capability | — |
| `mr_sad_workflow` | capability | — |

**Phase 2: Planning (D4)**

Eight scenarios exercising the PLAN node's LLM call.  Tests whether
the LLM picks the right next program given a workflow state.  All
reliability tests at threshold 0.8.

| Scenario | Asserts |
|----------|---------|
| `validate_pending_dont_stop` | Picks a validation program (not STOP) when R-free target met but validation not run |
| `plateau_with_validation_done` | Picks STOP when plateau + validation complete |
| `dont_stop_on_plateau_when_validation_pending` | Picks validation (not STOP) even when plateau, if validation pending |
| `cryoem_stop_target_with_validation` | Picks STOP when cryo-EM CC target met + validation done |
| `directives_force_stop_after_phaser` | Picks STOP when `after_program` directive applies |
| `directives_skip_program_honored` | Avoids program in `skip_programs` directive |
| `cryoem_next_step_after_placement` | Picks a valid non-STOP next step after `predict_and_build` |
| `dont_retry_same_failure` | Picks different program or strategy after a failure |

**Capability** tests pass at first success; the LLM is shown to be
capable of the right answer.  **Reliability** tests pass at
≥(threshold × max_runs) successes; the LLM is shown to be reliable.

When a reliability scenario passes overall but individual runs
failed (e.g. 4/5 PASS at threshold 0.8 = 1 individual fail), the
per-scenario line shows `[N fail]` so the failure is visible.

## Quick start

```sh
# Set at least one API key
export GOOGLE_API_KEY=...
# (and/or)
export OPENAI_API_KEY=...
export ANTHROPIC_API_KEY=...

cd $PHENIX/modules/cctbx_project/libtbx/langchain

# Run all tests
phenix.python tests/llm/run_llm_tests.py
```

## Framework self-tests (no API keys needed)

`tests/llm/tst_framework_self.py` contains 17 deterministic tests
that verify the framework's early-termination logic.  No real LLM
calls, no API keys required.  Run them to confirm the framework
itself works before invoking the API-driven tests:

```sh
phenix.python tests/llm/tst_framework_self.py
```

Output: `Framework self-test: 17 / 17 passed`.  If any fail, the
framework has a bug; fix that before running the API-driven tests.

## Common invocations

```sh
# All providers, all suites (Phase 1 + Phase 2 = 14 scenarios)
phenix.python tests/llm/run_llm_tests.py

# One suite only
phenix.python tests/llm/run_llm_tests.py --suite directive_extraction
phenix.python tests/llm/run_llm_tests.py --suite planning

# Smoke test (1 run per scenario, ~30s)
phenix.python tests/llm/run_llm_tests.py --max-runs 1

# Single provider
phenix.python tests/llm/run_llm_tests.py --provider openai

# Debug one scenario
phenix.python tests/llm/run_llm_tests.py \
    --scenario validate_pending_dont_stop --verbose

# Full output to stdout, no log file
phenix.python tests/llm/run_llm_tests.py --verbose --no-log
```

Exit codes:
- `0` — all PASS (or all SKIP for missing API keys)
- `1` — at least one FAIL or ERROR
- `2` — configuration error (no providers, invalid arg)

## Reading the logs

Each invocation creates a timestamped directory under
`tests/llm/logs/`:

```
tests/llm/logs/2026-05-15-15-44-12-phase1/
  run_summary.json                       # high-level: verdict per scenario
  details/
    af7mjs_bug_a_regression__google.log  # full per-run output
    af7mjs_bug_a_regression__openai.log
    ... etc
```

- `run_summary.json` is machine-readable.  Useful for tracking pass
  rates across weekly runs.
- `details/*.log` files are plaintext.  Each shows every run's input,
  raw LLM output (first 500 chars), parsed directives, and verdict
  reason.  Look here first when diagnosing a FAIL.

## Adding a new scenario

1. **Add a fixture file** in `fixtures/`.  For D2, it's preprocessed
   advice text (~300-1000 bytes).  Match the format of existing
   fixtures: numbered sections with `**Header**: value`.

2. **Write an expectation function** in `tst_directive_extraction.py`.
   Signature is `(parsed_directives) -> (passed: bool, why: str)`.
   See the existing helpers for examples.

3. **Add the Scenario** to `build_scenarios()` in
   `tst_directive_extraction.py`.  Pick `test_type="capability"`
   unless you specifically care about reliability rate.

4. **Run it**:
   ```sh
   phenix.python tests/llm/run_llm_tests.py --scenario your_new_name
   ```

## ERROR vs FAIL

The runner distinguishes:

- **FAIL** — the LLM responded but its output was wrong.  The model
  is making bad decisions for this scenario.
- **ERROR** — the LLM call itself failed (API timeout, rate limit,
  malformed JSON the parser can't recover from).  We never got a
  semantic answer.

If ANY runs of a scenario succeed, the verdict is PASS or FAIL based
on the success count.  Only when ALL runs error out is the verdict
ERROR.  This makes "the API is broken" easy to distinguish from
"the LLM is making wrong choices".

## Early termination

To save time and money, the runner stops each scenario as soon as
the verdict is certain:

- **Capability test**: PASS at first successful run.  FAIL only after
  `max_runs` (default 5) failures.
- **Reliability test**: PASS at first time `n_passed >= ceil(threshold
  × max_runs)`.  FAIL when remaining runs can't reach the threshold.

For `max_runs=5, threshold=0.8`: PASS at 4 successes, FAIL at 2
failures.

For the AF_7mjs reliability test (`max_runs=5, threshold=0.9`): PASS
at 5 successes, FAIL at 1 failure.

In the common case (LLM gets things right on the first try), most
scenarios finish in 1 LLM call instead of 5.

## Fixture capture (Phase 2/3)

`tools/capture_session_fixture.py` extracts and sanitizes real state
from production `agent_session.json` files for use in Phase 2 (planning)
and Phase 3 (preprocessing) tests.

Phase 1 doesn't need it — D2 fixtures are simple advice text.  See the
tool's docstring for usage.

## Cost and time

Realistic full Phase 1 run (most scenarios pass first try):

- Wall time: ~1-2 min
- LLM calls: ~10-15
- Cost: ~$0.01-0.05 depending on providers

Worst case (everything fails):

- Wall time: ~5 min
- LLM calls: 30 (6 scenarios × 5 runs × 1 provider)
- Cost: ~$0.10

Cheap enough to run weekly and before releases.  NOT registered in
`tests/run_all_tests.py` — these are not for every-commit CI.
