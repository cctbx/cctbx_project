# CHANGELOG — v116

## Version 116.12 (Cryo-EM Validation-Aware Auto-Stop)

### Summary

Two coordinated fixes for the AF_7mjs Stage 5 (validation) skip
bug. The agent ran the correct 5-stage plan but auto-stopped after
Stage 4 (refinement, CC=0.78) without running validation
(phenix.molprobity).

Root cause: asymmetry between X-ray and cryo-EM stop logic in
`agent/metrics_analyzer.py`. The X-ray path correctly defers
auto-stop until validation runs; the cryo-EM path didn't. The two
paths drifted out of sync.

### Modified Files (2 production files)

| File | Phases | Changes |
|------|--------|---------|
| `agent/metrics_analyzer.py` | Fix #1 | +25 / -4 lines. Mirror the X-ray validation_done check in the cryo-EM `latest_cc > 0.70` branch. Recognizes `phenix.molprobity` and `phenix.validation_cryoem`. When validation isn't done, sets `should_stop=False` and `suggest_validation=True` instead of stopping immediately. |
| `agent/graph_nodes.py` | Fix #2 | +41 / -0 lines. Two additions in the PLAN node's auto-stop chain: (a) new `elif` between `plan_has_pending_stages` and the final `else` — suppresses AUTO-STOP when `workflow_state.step_info.step == "validate"` AND `validation_done == False`; (b) diagnostic context dump in the AUTO-STOP path logging `plan_has_pending_stages`, `step`, `validation_done`, `after_program`, `experiment_type` to `debug_log`. |

### Bug Details

#### Symptom

```
[STATE]  cryoem_refined
  Model refined, need validation - Goal: Ensure model meets quality standards (CC: 0.783)
[STOP]
  SUCCESS: Map-model CC (0.783) above target (0.70)
No command was generated for program 'unknown'. Ending session.
```

AF_7mjs (cryo-EM, half-maps + sequence, PredictAndBuild) had a
plan with 5 stages. After Stage 4 (`phenix.real_space_refine`)
produced `map_cc=0.7832`, the agent should have advanced to Stage
5 (`phenix.molprobity`). Instead it auto-stopped. The state
machine even reported "Model refined, need validation" but the
threshold-based stop fired anyway.

Session JSON confirmed:
- Plan had 5 stages (correct)
- Stages 1-3 complete (correct)
- Stage 4 active, `cycles_used=2`, `success_criteria: model_map_cc > 0.7` met
- Stage 5 `status="pending"`, `programs=["phenix.molprobity"]`, **never ran**
- Cycle 5 reasoning: `"Automatically stopping: SUCCESS: Map-model CC (0.783) above target (0.70). Map CC: 0.783 - TARGET REACHED"`
- Cycle 5 program: empty (None); F7 fallback printed `'unknown'`

#### Root cause analysis

The X-ray and cryo-EM stop paths in `analyze_metrics_trend` should
be symmetric on the "success → validation" question. They were
not.

**X-ray path (lines 365-385) — correct:**
```python
if latest_r_free < success_threshold:
    validation_done = any(
        m.get("program") in ("phenix.molprobity", "phenix.model_vs_data")
        for m in metrics_history
    )
    if validation_done:
        result["should_stop"] = True
    else:
        result["should_stop"] = False
        result["suggest_validation"] = True
```

**Cryo-EM path (lines 506-511) — bug:**
```python
if latest_cc > 0.70:
    result["should_stop"] = True   # No validation check
    return result
```

The cryo-EM path was probably written first and the validation_done
check was added to the X-ray path later. The two paths drifted.

### Three logical changes

| Change | Location | What |
|--------|----------|------|
| **Fix #1: validation_done check** | `metrics_analyzer.py` `_analyze_cryoem_trend`, `latest_cc > 0.70` branch | Mirror X-ray pattern. When validation not done, `should_stop=False` + `suggest_validation=True` + trend_summary `"Map-model CC: X.XXX - TARGET REACHED, VALIDATE BEFORE STOPPING"`. Accepts `phenix.molprobity` and `phenix.validation_cryoem` as valid cryo-EM validation programs. |
| **Fix #2a: validate-step suppression elif** | `graph_nodes.py` `plan()` auto-stop chain, between `plan_has_pending_stages` and `else` | Defense-in-depth. Suppresses AUTO-STOP when `workflow_state.step_info.step == "validate"` AND `not validation_done`. Catches cases where `plan_has_pending_stages` doesn't fire. |
| **Fix #2b: diagnostic context dump** | `graph_nodes.py` `plan()` else (AUTO-STOP) branch | Logs `plan_has_pending_stages`, `step`, `validation_done`, `after_program`, `experiment_type` to `debug_log` before stopping. Makes future regressions debuggable. |

### Priority order in the elif chain (after Fix #2)

```
1. advice_changed       (user provided new advice on resume)
2. after_program        (user-specified stop target)
3. user_wants_ligandfit (ligandfit prerequisite pending)
4. plan_has_pending_stages (expert plan has pending stages)
5. step=="validate" AND !validation_done  ← NEW (defense-in-depth)
6. else AUTO-STOP (now with diagnostic context dump)
```

Rationale: user/plan intent (1-4) wins over engine-derived intent
(5). The new check is the last line of defense.

### Tests (2 new files, 27 new test cases)

| File | Phase | Cases | Covers |
|------|-------|-------|--------|
| `tst_cryoem_stop_validation.py` (S27) | Fix #1 | 16 (31 assertions) | AF_7mjs regression, both validation programs recognized, X-ray validation programs explicitly don't count, CC just-above threshold, CC exactly at threshold (strict `>`), validation_done with EXCESSIVE / PLATEAU paths, X-ray path unchanged |
| `tst_plan_autostop_validation_suppression.py` (S28) | Fix #2 | 11 (28 assertions) | New elif fires when step=="validate" + !validation_done; doesn't fire when validation_done=True; doesn't fire when step != "validate"; priority order with other suppressors; diagnostic dump appears only on AUTO-STOP; dump records correct values |

Both test files use the libtbx-stub pattern so they run standalone
without a PHENIX install.

Registered in `tests/run_all_tests.py` as suites S27 and S28.

### Verification

```
=== Fix #1 tests ===          (pre-fix:)
31 passed, 0 failed           18 passed, 13 failed

=== Fix #2 tests ===          (pre-fix:)
28 passed, 0 failed           10 passed, 13 failed

=== Regression check (existing v116.11 suites) ===
tst_stop_condition_false_positive: 15 passed, 0 failed
tst_general_resolver:              21 passed, 0 failed
tst_dock_and_stop:                  5 passed, 0 failed
tst_standalone_consistency:         8 passed, 0 failed
tst_initialize_plan_smoke:          9 passed, 0 failed
tst_file_encoding:                  7 passed, 0 failed
tst_program_requirements:          40 passed, 0 failed
```

### Why TWO fixes if Fix #1 is the primary

Fix #1 alone should prevent the bug (don't set `should_stop=True`
until validation is done). Fix #2 adds:

1. **Unsolved-mystery insurance.** We could not determine from the
   AF_7mjs session JSON why `plan_has_pending_stages` didn't
   suppress the AUTO-STOP. The function should have returned True
   against the final plan state. Fix #2's diagnostic logging is
   specifically designed to catch this on future runs.

2. **Defense-in-depth on a class of bugs.** Other code paths
   might set `should_stop=True` outside `_analyze_cryoem_trend`
   (PLATEAU, EXCESSIVE, future additions). Fix #2 protects
   against any future "stop fires while validate is pending"
   regression.

3. **No cost.** Fix #2's new elif only fires when the existing
   chain would have stopped AND `step="validate"` AND
   `validation_done=False`. If Fix #1 works, Fix #2's new elif
   never fires in practice.

### Compatibility

- **X-ray workflows unchanged.** Fix #1 only modifies the cryo-EM
  branch. The X-ray pattern it mirrors was already there.
- **Workflows with validation already done**: `should_stop=True`
  fires as before. No behavior change.
- **No directives or contract changes.** No client protocol updates.
- **AF_7mjs and similar cryo-EM tutorials**: workflow now correctly
  advances to validation stage instead of stopping prematurely.

### Known limitations (not blocking AF_7mjs)

| Limitation | Impact |
|------------|--------|
| **Cosmetic 'unknown' in F7 fallback** | "No command was generated for program 'unknown'" comes from `_run_single_cycle`: `decision_info.get('program', '') or 'unknown'`. When `program=None` (stop intent), `None or 'unknown'` = `'unknown'`. Cleaner: show `'STOP'` or `stop_reason`. Not blocking; session correctly ends. |
| **Unknown why `plan_has_pending_stages` didn't suppress** | Function should return True against AF_7mjs session JSON. Possibilities: stale session.data at check time, transport-layer propagation issue, or deployed version differs from analyzed version. Fix #2's diagnostic logging will surface this on future occurrences. |
| **Map-CC trend_summary format change** | Before: `"Map-model CC: 0.620 → 0.780 (+25.8% last cycle)"`. After (validation done): `"Map-model CC: 0.780 - ABOVE TARGET"`. After (validation not done): `"Map-model CC: 0.780 - TARGET REACHED, VALIDATE BEFORE STOPPING"`. UI log parsers that match the trend format may need to be updated. Mirrors X-ray's pre-existing pattern. |
| **MTZ warning in cryo-EM workflows** | Pre-existing diagnostic: `"WARNING: Refinement completed but best_files[map_coeffs_mtz] is EMPTY"`. Cryo-EM doesn't produce MTZ — the slot doesn't apply. Should be gated on `experiment_type == "xray"`. Not blocking. |

## Version 116.11 (AF_7mjs Stop Condition Fix + Test Cleanup)

### Summary

Two coordinated changes shipped together:

1. **Directive extractor fix** for an AF_7mjs (and likely other
   cryo-EM tutorial) regression: the preprocessor-inserted
   `**Stop Condition**: None` header was being misinterpreted as
   real user stop intent, causing the planner to skip prerequisite
   stages and produce wrong plans.

2. **Test file cleanup** — removed `libtbx.find_unused_imports`
   warnings from 5 test files by switching to real symbol references
   (no `# noqa` shortcuts).

The directive extractor fix was the trigger; the test cleanup
shipped in the same cycle since it was in flight.

### Modified Files (1 production file, 5 test files)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | +98 / -17 lines. Three logical changes (see Bug Details below): (1) strengthened the regex used at all three places that detect preprocessor formatting — handles markdown bold (`**Header**:`), numbered list prefixes (`7. Header:`), and bullet markers (`- Header:`); (2) defense-in-depth Stop Condition strip inside `_resolve_after_program` before the `\bstop\b` check; (3) suppression of `start_with_program` writes from `_resolve_after_program` when advice is preprocessor output, since multi-action signals in descriptive Primary Goal prose are not user prescription. |
| `tests/tst_file_encoding.py` | Switched from `from yaml_loader import _load_yaml_file` to `from knowledge import yaml_loader` + `hasattr(yaml_loader, '_load_yaml_file')` availability probe. Real symbol reference, no `# noqa`. |
| `tests/tst_general_resolver.py` | Removed unused `import re` (was at top, never referenced). |
| `tests/tst_skip_to_program.py` | Removed unused `STAGE_FAILED` from the plan_schema import tuple. |
| `tests/tst_standalone_consistency.py` | Switched from importing `_ACTION_TABLE` directly to importing the `directive_extractor` module and checking `directive_extractor is not None` as an availability probe. The actual symbol is loaded later by `_load_action_table()` which has its own fallback chain. |
| `tests/tst_user_advice_filter.py` | Merged the two import attempts (libtbx and local) into a single try/except chain. Eliminated the dead `_real` alias that existed only to satisfy the linter. |

### Bug Details

#### v116.11 root cause: latent strip-regex weakness

`_strip_preprocessor_stop_condition` at line 44 of
`directive_extractor.py` exists specifically to remove
preprocessor-inserted headers at the entry of both extraction
paths (`extract_directives` and `extract_directives_simple`). Its
regex required headers to start with the literal word:

```python
r'^[ \t]*stop\s+condition\s*:\s*[^\n]*$'
```

But the advice-preprocessor LLM produces output with **markdown
bold and numbered list prefixes** like:

```
7. **Stop Condition**: None
```

The strip never matched this format. Before v115.10, no code
path did a bare `\bstop\b` search on the un-stripped text, so the
bug was latent. v115.10's new `_resolve_after_program` introduced
the first sensitive consumer.

#### Symptom: AF_7mjs three-stage plan

When the agent ran AF_7mjs with the v115.10 directive extractor,
the resolver matched `\bstop\b` against the un-stripped
`**Stop Condition**` header, set `_has_stop=True`, and combined
with the multi-action detection in the Primary Goal text
("Run PredictAndBuild ... rebuild ... refine") produced
directives like:

```yaml
after_program: phenix.real_space_refine
skip_validation: True
start_with_program: phenix.predict_and_build
```

The planner interpreted `start_with_program` as a skip-prerequisites
directive and produced a 3-stage plan starting at Stage 3 instead
of the correct 5-stage cryo-EM plan starting with `mtriage`. The
March 30, 2026 working run produced the correct 5-stage output;
something between then and May 11 (v115.10 deployment) exposed
the latent bug.

#### Three logical changes in directive_extractor.py

| Change | Locations | What |
|--------|-----------|------|
| **Strengthen preprocessor-format regex** | `_strip_preprocessor_stop_condition` (Stop Condition strip, Primary Goal strip, and `_PREPROCESSOR_SIGNATURES` detection) + the pre-existing `_is_preprocessed` check at line 2860 in `extract_directives_simple` | Adds optional numbered list prefix (`[\d]+\.`), bullet markers (`[-*]`), and markdown bold (`\**`) wrappers around the header keyword. Three locations, identical pattern, for consistency. |
| **Defense-in-depth Stop Condition strip** | Top of `_resolve_after_program`, before the `\bstop\b` check | Strips the Stop Condition header again before computing `_has_stop`. Under normal flow this is redundant with the upstream strip, but catches future code paths that bypass `_strip_preprocessor_stop_condition`. |
| **Suppress `start_with_program` for preprocessed advice** | Multi-action branch of `_resolve_after_program` | New `_is_preprocessed` check inside the resolver. When the advice is detected as preprocessor output, the multi-action branch no longer writes `start_with_program`. Reasoning: in preprocessed advice, multi-action mentions come from descriptive Primary Goal prose ("Run PredictAndBuild ... dock/trim ... rebuild ... refine"), not user prescription. Setting `start_with_program` from descriptive prose makes the planner skip prerequisite stages. Real user prose like "run phaser and refine" is unaffected — no preprocessor signatures, so `start_with_program=phenix.phaser` still gets set as before. |

### Tests (1 new file, 15 tests)

| File | Tests | Covers |
|------|-------|--------|
| `tst_stop_condition_false_positive.py` (S26) | 15 | Primary AF_7mjs regression; markdown header variants (5); upstream-strip preserves user prose (3); Primary Goal strip strengthened; `_resolve_after_program` directly for resolver-isolated behavior (negation, legitimate stop, multi-action with/without stop, real user prose for start_with); start_with suppression for preprocessed advice. |

Registered in `tests/run_all_tests.py` as suite S26.

The test file uses the libtbx-stub pattern so it runs without a
real PHENIX install. Tests that target specific resolver behavior
call `_resolve_after_program` directly to isolate from
`intent_classifier`, which is a separate, legitimate extraction
path that can also write `after_program` through its own
mechanism.

### Verification

- **Reproducer**: full AF_7mjs processed advice now produces empty
  or minimal directives (matching the March 30 working behavior).
  No `after_program` write, no `start_with_program` write.
- **Real PHENIX**: 15/15 tests pass via `phenix.python
  tst_stop_condition_false_positive.py` after the test-environment
  sensitivity fix (`test_only_primary_goal_extracted_no_stop` was
  initially too strict — it didn't account for `intent_classifier`
  legitimately writing `after_program` for simple "Refine the
  model." advice; rewrote the test to call `_strip_preprocessor_stop_condition`
  in isolation).
- **No regressions**: all existing test suites pass
  (`tst_general_resolver` 21/21, `tst_dock_and_stop` 5/5,
  `tst_program_requirements` 40/40, `tst_standalone_consistency`
  8/8, `tst_initialize_plan_smoke` 9/9, `tst_file_encoding` 7/7).

### Self-Review Findings

The fix went through five review iterations before settling on the
current shape. Each iteration found a problem with the previous
diagnosis:

1. **First diagnosis (wrong shape):** Strip the header inside
   `_resolve_after_program` only. This is a symptom fix — it
   treats the bug at the consumer, not the cause. The right fix
   point was the upstream `_strip_preprocessor_stop_condition`.

2. **Second diagnosis (incomplete):** Mark resolver writes with
   `_set_by_pattern=True` so the intent-driven clearing block
   could remove them. But the clearing block only handles
   `after_program`, not `start_with_program`. The actual
   problematic write on AF_7mjs was `start_with_program`.

3. **Third diagnosis (right shape):** Suppress `start_with_program`
   writes when advice is preprocessed. Simpler than mark-and-clear;
   addresses the root semantic confusion (descriptive prose ≠ user
   prescription).

4. **Fourth iteration:** Found that the pre-existing
   `_is_preprocessed` regex at line 2860 of `extract_directives_simple`
   was inconsistent with my upstream fix — it still used the old
   regex that didn't handle markdown. Strengthened for consistency.

5. **Fifth iteration:** Test failure in real PHENIX revealed that
   `intent_classifier` writes `after_program` through a separate
   path. Rewrote affected tests to use `_resolve_after_program`
   directly when targeting resolver-specific behavior.

### Known limitations (not blocking AF_7mjs)

| Limitation | Impact |
|------------|--------|
| **Conditional stop in Special Instructions** | Text like "Stop the analysis if R-free > 0.5" in Special Instructions could still false-positive `_has_stop=True`. AF_7mjs's actual Special Instructions don't contain stop words. Pre-existing fragility. |
| **DRY violation across 3 `_PREPROC_SIGS` locations** | Same signature list and regex appears at three locations in `directive_extractor.py`. Cosmetic; future refactor should extract to a module-level helper. |
| **Line 2860 `_is_preprocessed` covered only indirectly** | My test suite verifies the strengthened regex via the strip function (which uses the same pattern). The line 2860 instance is mechanically a copy of tested code, but a direct test would be cleaner. |

### Test Cleanup Details

The unused-imports fix used **real symbol references**, not
`# noqa: F401` suppression. The `libtbx.find_unused_imports`
tool requires that imported names actually be referenced; a
`noqa` comment doesn't satisfy it.

The general pattern was: where a test imports a symbol only to
verify the module is available, switch to importing the module
itself and using `hasattr()` or `is not None` as the availability
probe. This is both more honest about what the test depends on
and less brittle if the symbol gets renamed.

---

## Version 116.10 (Reliability and Classification Cleanup)

### Summary

Six bugs in advice parsing, LLM prompting, program selection, and
plan classification, plus drift-detection machinery for the wire
contract and the client's plan-generation logic. The cleanup was
scoped from a single batch-tutorial finding ("predict and stop"
sessions were AUTO-STOPping) and grew as adjacent bugs surfaced
during investigation. Five files modified, 89 new tests across
7 new test files plus one augmented file. Rollout order in
testing: 1 → 6a → 4b → 2 → 6b → 3a → 3d.

### Modified Files (5 production files)

| File | Phases | Changes |
|------|--------|---------|
| `agent/rules_selector.py` | 1 | `_apply_user_advice` `stop_condition_patterns` extended with `"and stop"`, `"then stop"`, `", stop"` so "predict and stop" stops being collapsed to `["STOP"]`. |
| `agent/workflow_engine.py` | 4b, 6b | Phase 4b: new `_filter_programs_missing_data_inputs` strips `xtriage` when `not has_data_mtz and not has_phased_data_mtz`, strips `mtriage` when no map; matching guards added to `_check_program_prerequisites` for defense in depth. Phase 6b: top of `_detect_xray_step` routes sequence-only sessions directly to `obtain_model` (state becomes `xray_analyzed`). |
| `agent/contract.py` | 2 | `CURRENT_PROTOCOL_VERSION` 3 → 5 (matches highest field version after silent drift). New `validate_contract()` returns `(ok, errors)` enforcing three invariants: CURRENT ≥ max field version, MIN ≤ CURRENT, MIN ≥ 1. |
| `knowledge/prompts_hybrid.py` | 6a | `_format_directives_for_prompt` after-program lines replaced: "CRITICAL: You MUST run X" → "Stop target: X. If in VALID PROGRAMS choose it; if not, choose an appropriate prerequisite. Never pick outside VALID PROGRAMS." |
| `phenix/programs/ai_agent.py` | 2, 3a, 3d | Phase 2: `_get_protocol_version()` fallback `return 3` → `return 5`. Phase 3a: `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS` extracted as module-level frozensets near line 815; three call sites in `_initialize_plan_inner` updated to reference them (pure refactor, 51 program×intent traces verified identical). Phase 3d: `phenix.map_symmetry` and `phenix.dock_in_map` added to `_STANDALONE_PROGRAMS` (behavior change: skips plan generation, routes through workflow_engine state machine instead). |

### Bug Details

| Bug | Phase | Symptom | Root Cause | Fix |
|-----|-------|---------|------------|-----|
| 1 | 1 | "predict and stop" sessions immediately AUTO-STOP at session start | `_apply_user_advice` saw `stop` keyword without matching any of `stop after / stop when / stop once / stop if / stop condition / stop at` → fell through to "treat as immediate stop" → `valid_programs` collapsed to `["STOP"]` | Extended `stop_condition_patterns` with `"and stop"`, `"then stop"`, `", stop"` |
| 2 | 6a | LLM picks `predict_and_build` before `xtriage` has run, hitting prerequisite failure | Prompt told LLM "CRITICAL: You MUST run X before stopping. If it's in VALID PROGRAMS, choose it NOW." — overrode prerequisite logic | Reframed as "Stop target: X. If in VALID PROGRAMS choose it; if not, choose an appropriate prerequisite" |
| 3 | 4b | xtriage offered as valid program even when no `.mtz` is uploaded; session stalls | YAML lists xtriage at analyze step unconditionally; no data-input filter | New `_filter_programs_missing_data_inputs` strips xtriage when no .mtz, mtriage when no map; matching prereq guard in `_check_program_prerequisites` for defense in depth |
| 4 | 2 | New v3/v4/v5 fields added to `SESSION_INFO_FIELDS` while `CURRENT_PROTOCOL_VERSION` stuck at 3 — silent drift | No invariant check between `CURRENT_PROTOCOL_VERSION` and registered field versions | Bumped CURRENT to 5; added `validate_contract()` invariant function; augmented existing `tst_contract_compliance.py` with drift-detection test |
| 5 | 6b | Sequence-only X-ray sessions stuck at `xray_initial` (analyze step); xtriage filtered, no fallback | `_detect_xray_step` returned analyze for `not xtriage_done` regardless of whether xtriage could even run | Routing branch at top of `_detect_xray_step`: when no .mtz and has sequence → `obtain_model` (state becomes `xray_analyzed`); `predict_and_build` available via YAML `has: sequence` condition |
| 6 | 3a | `_initialize_plan_inner` standalone-programs list duplicated inline at two call sites; comment admitted the duplication; classification couldn't be tested | Original implementation kept the two tuples in sync by hand | Extracted `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS` as module-level frozensets; new `tst_standalone_consistency.py` enforces alignment with `directive_extractor._ACTION_TABLE` |
| 7 | 3d | "dock and stop" with sequence + map (no model) generates a plan, `skip_to_program(dock_in_map)` marks predict_and_build SKIPPED, dock_in_map then fails at runtime | `dock_in_map` was a "full-plan target" (caught by the v116.10 elif); skip_to_program incorrectly removed the predict prerequisite | Reclassified `phenix.dock_in_map` and `phenix.map_symmetry` as `_STANDALONE_PROGRAMS`; `_initialize_plan_inner` now skips plan generation, lets `workflow_engine` handle prerequisites via state machine (analyze → predict → dock) |

### Tests (7 new files, 1 augmented, 89 new test cases)

| File | Phase | Cases | Covers |
|------|-------|-------|--------|
| `tst_user_advice_filter.py` | 1 (S13) | 16 | "X and stop" / "X then stop" / "X, stop" preserved; recognized stop-conditions still respected |
| `tst_after_program_prompt.py` | 6a (S14) | 12 | Prompt uses target-not-now framing; aggressive "MUST" language removed |
| `tst_data_input_filter.py` | 4b (S15) | 22 | xtriage filtered when no .mtz; mtriage filtered when no map; prereq guard handles directive re-injection |
| `tst_protocol_version.py` | 2 (S16) | 15 | `validate_contract()` invariants via synthetic probes (drift detection, bounds checking) |
| `tst_contract_compliance.py` | 2 (augmented) | +1 | New `test_contract_validate_passes` adds drift-detection entry point |
| `tst_sequence_only_routing.py` | 6b (S17) | 10 | Sequence-only sessions route to `obtain_model`; doesn't fire when data is present or sequence absent |
| `tst_standalone_consistency.py` | 3a (S18) | 8 | `_STANDALONE_PROGRAMS` aligned with `_ACTION_TABLE`; drift detection for new actions |
| `tst_dock_and_stop.py` | 3d (S19) | 5 | Phase 3d behavior change documented; blast-radius regression guard (48 unrelated traces unchanged) |

Registered in `tests/run_all_tests.py` as suites S13–S19.

Run with: `python3 tests/run_all_tests.py` or individually
`python3 tests/tst_user_advice_filter.py`, etc.

### Verification

| Fix | Verification |
|-----|--------------|
| Bug 1 (user advice "X and stop") | ✅ 16 unit tests pass; 6/16 fail on pre-fix code |
| Bug 2 (after-program prompt) | ✅ 12 unit tests pass; 6/12 fail on pre-fix |
| Bug 3 (data-input filter) | ✅ 22 unit tests pass; 18/22 fail on pre-fix; integration smoke test confirms xtriage absent from valid_programs when no .mtz |
| Bug 4 (protocol hygiene) | ✅ 15 + 1 unit tests pass; 8/15 fail on pre-fix; verified `tst_contract_compliance.py:test_protocol_version_consistency` now passes (was failing because fallback `return 3` ≠ CURRENT `5`) |
| Bug 5 (sequence-only routing) | ✅ 10 unit tests pass; 4/10 fail on pre-fix; end-to-end integration test confirms state becomes `xray_analyzed` with `valid_programs=[phenix.predict_and_build]` |
| Bug 6 (Phase 3a refactor) | ✅ 8 consistency tests pass; behavioral equivalence verified by 51 (program × intent) decision-tree traces — all identical pre/post |
| Bug 7 (Phase 3d behavior change) | ✅ 5 unit tests pass; 3/5 fail on pre-fix; blast-radius check confirms 48 unrelated traces unchanged. ⏳ Integration test against tutorial corpus pending. |

### Self-Review Finding (Phase 3a/3d split)

Phase 3 was initially delivered as a single phase that bundled the
de-duplication refactor with two new programs added to
`_STANDALONE_PROGRAMS` (silently changing behavior). Self-review
caught this:

- Refactor PRs that change behavior should be split. "Cleaner"
  behavior is still a behavior change.
- The fix was to split into Phase 3a (verified by 51 trace-equivalence
  checks) and Phase 3d (explicit behavior change, with `tst_dock_and_stop.py`
  including a 48-trace blast-radius regression guard).

This is recorded in `OVERVIEW.md` under Active Development as a
principle for future cleanups.

### Known Issues (deferred to future work)

| Issue | Description |
|-------|-------------|
| Phase 3d integration verification | Decision-tree traces verified, but no tutorial currently exercises "dock and stop with sequence + map". A tutorial that does would close the verification gap. Rollback is one line if a regression appears. |
| Systematic input-availability audit | Phase 4b only filters xtriage and mtriage. A systematic "every program declares its inputs" audit is unfinished. The next program added with hard input requirements should be added to the filter at the same time. |
| Phase 6a "appropriate prerequisite" language | Shipped Option A (ship as-is) over Options B (add hint) and C (surface `rules_priority`). If integration testing shows the LLM picks the wrong prerequisite, revisit. |
| Redundant prediction-only allowance | The v116.10 allowance in `_check_program_prerequisites` is mostly redundant after Phase 6b but retained as defense in depth. A future cleanup could remove it if edge cases (e.g., `xtriage_done=True` with no data) are verified to reach `predict_and_build` through other paths. |

### Post-Phase-5 Addition: CC Key Extraction Fix (S20)

#### Summary

A successful cryo-EM workflow that produced metrics and a structure
model displayed the misleading banner "SESSION STOPPED -
INCOMPLETE / No structure model available." Root-caused to a key
naming inconsistency in `_generate_structure_report` and patched
with a defensive multi-key lookup.

#### Bug Details

| Symptom | Successful cryo-EM run reports "SESSION STOPPED - INCOMPLETE" |
| Root Cause | `_generate_structure_report` looked up CC under `map_model_cc` (the form programs PRINT in logs), but the cycle metrics dict stores it under `model_map_cc` (the canonical storage key, see schema at `ai_agent.py:9083-9088`). The v115.05 author wrote the lookup against the printed form. For cryo-EM workflows without R-free, this meant `_best_cc=None` → `_metrics_good=False` → stopped-report path fires despite success. |
| Fix | Defensive multi-key lookup at lines 8794-8797 and 8822-8823: checks `model_map_cc` (canonical), `map_model_cc` (legacy/unnormalized), `map_cc` (workflow_engine variant), `cc_mask` (real_space_refine output). Works regardless of which form actually appears in the dict. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `phenix/programs/ai_agent.py` | 8794-8797, 8822-8823 | Two surgical edits adding `model_map_cc` to lookup lists |

#### Tests

`tst_cc_key_extraction.py` (11 tests, suite S20): 2 source-level
invariants, 5 behavioral key-variant tests, 3 edge cases, 1
motivating-case test reproducing the Tom-reported scenario.

#### Self-Review Note

Initially framed as "an open-and-shut typo." That was sloppy.
`map_model_cc` is a legitimate spelling — it appears as a phenix
program name and in regex-captured raw log lines. The bug is real
but it's "wrong form of the name" not "typo": the author used the
printed form rather than the stored form. The fix's defensive
multi-key lookup handles both forms regardless.

### Post-Phase-5 Addition: File Encoding Fix (S21)

#### Summary

A user on Chinese-locale Windows reported a `UnicodeDecodeError:
'gbk' codec can't decode byte 0x94` crash at startup while loading
`programs.yaml`. Root cause: Python's `open()` without explicit
`encoding=` falls back to the system code page (`gbk` for
Chinese-locale Windows, `cp1252` for English, `cp932` for
Japanese), which crashes on any non-ASCII content in a UTF-8 file.

The fix was applied in two passes: initially to the 4 files
directly implicated by the crash trace, then to the **entire
v116.10 codebase** (63 additional files surfaced by a recursive
grep audit).

#### Bug Details

| Symptom | `UnicodeDecodeError: 'gbk' codec can't decode byte 0x94 in position 987` on PHENIX startup, Chinese-locale Windows |
| Root Cause | `open(path)` without `encoding=` uses `locale.getpreferredencoding()`. On non-UTF-8 Windows locales, any UTF-8 file with non-ASCII content (curly quotes, em dashes, accented letters) crashes the read. |
| Fix | Add `encoding='utf-8'` to every text-mode `open()` call. YAML, JSON, PDB, and report files are all UTF-8 by spec, so this is correct behavior, not a workaround. |

#### Files Modified

**First pass (the reported crash site + immediate neighbors):**

| File | Sites |
|------|-------|
| `knowledge/yaml_loader.py` | 1 (the reported crash) |
| `phenix/programs/ai_agent.py` | 3 |
| `agent/directive_validator.py` | 1 |
| 4 test files | 11 |

**Second pass (codebase-wide audit, 305 sites across 63 files):**

| Subdirectory | Offenders | Files |
|--------------|-----------|-------|
| `agent/` | 53 | 22 |
| `knowledge/` | 3 | 3 |
| `utils/` | 1 | 1 |
| `tests/` | 248 | 37 |

Patches applied: `open(path)` → `open(path, encoding='utf-8')`;
`open(path, 'r')` → `open(path, 'r', encoding='utf-8')`;
`open(path, errors='X')` → `open(path, encoding='utf-8',
errors='X')`. Binary opens, `Popen`, `urlopen`, and `.open()`
method calls were correctly excluded.

#### Tests

`tst_file_encoding.py` (7 tests, suite S21):

- 3 per-file source scans (`yaml_loader`, `ai_agent`, `directive_validator`)
- 2 directory-scan tests (`test_all_production_code_uses_utf8`,
  `test_all_test_code_uses_utf8`) — walk the entire production and
  test trees automatically, so new files added in the future are
  covered without test updates
- 2 empirical tests (write/read UTF-8, the user's exact byte pattern)

The scanner uses Python's `tokenize` module rather than regex
heuristics, so docstrings and string-literal occurrences of
`open(` are correctly excluded.

#### Verification

| Check | Result |
|-------|--------|
| Re-audit of patched files (tokenize-based) | 0 offenders |
| Python syntax (all 63 files) | All parse |
| Line-count drift | 0 (every patch is in-place) |
| 1:1 substitution (added = removed) | 306 = 306 |
| `tst_file_encoding.py` on patched | 7/7 pass |
| `tst_file_encoding.py` on unpatched | Correctly flags 302 |

#### User Workaround (No Code Change)

If a user can't deploy the patched files immediately, setting
`PYTHONUTF8=1` in their shell before launching phenix forces
Python 3.7+ to use UTF-8 for all file I/O regardless of system
locale.

#### Self-Review Note

The 305 patches are mechanically uniform but lack per-patch
explanatory comments (the first encoding fix added a 5-line
comment to each file; the automated patcher for the codebase-wide
pass did not). The rationale lives in this CHANGELOG entry rather
than per-patch comments. Future maintainers cleaning up the
codebase should know the `encoding='utf-8'` convention is
load-bearing for non-UTF-8 Windows locales.

### Post-Phase-5 Addition: Ligand Workflow Restart Fix (Phase 6c, S22)

#### Summary

The nsf-d2-ligand tutorial restarted its reasoning at cycle 3:
after a successful refinement (cycle 1) and a successful ligand
fit (cycle 2), cycle 3 reported `STATE: xray_initial / Need to
analyze data quality first` and the LLM produced the reasoning
*"As this is the first refinement step, I will set
'generate_rfree_flags=true'."* The same pattern repeated at
cycles 4 and 5. Three interacting issues were diagnosed; two
fixes resolve all three.

#### Bug Details

| Bug | Symptom | Root Cause | Fix |
|-----|---------|------------|-----|
| **1** | Every cycle reports `STATE: xray_initial` even after refine and ligandfit succeeded | `_detect_xray_step` returns "analyze" any time `xtriage_done=False`, regardless of downstream progress. When the user supplies a pre-placed model + `start_with_program=phenix.refine` (which skips xtriage), the state never advances | Added a `past_analysis` check that mirrors `_detect_cryoem_step`: if any downstream program has completed (`refine_done`, `ligandfit_done`, `phaser_done`, etc.) or a positioned model exists on disk, advance past analyze |
| **2** | After ligandfit succeeds, the "best model" pointer stays on the previous unliganded refined model | `best_files_tracker` scores `ligand_fit_output` (stage 105) below `refined` (stage 100 + R-free contribution ~22 = ~122). LigandFit's output (no R-free metrics) loses the scoring contest, so the LLM is told the model is still the unliganded one | Extended the existing `with_ligand` metric inheritance to also cover `ligand_fit_output`. One-condition change: `stage == "with_ligand"` → `stage in ("with_ligand", "ligand_fit_output")` |
| **3** | LLM mechanically re-applies first-cycle directives (`generate_rfree_flags=true`) at cycles 3, 4, 5 | Observable consequence of Bugs 1 and 2 — when the agent presents a stale worldview (state=xray_initial + unliganded model as "current"), the LLM correctly applies first-cycle directives because the worldview says we're at cycle 1 | No separate fix needed; resolves automatically once the LLM sees correct state (Bug 1) and current model (Bug 2) |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `agent/workflow_engine.py` | 1040-1086 | Replaced 3-line early return with 45-line `past_analysis` check (mirrors `_detect_cryoem_step` at lines 1283-1299) |
| `agent/best_files_tracker.py` | 564-583 | Extended condition `stage == "with_ligand"` to `stage in ("with_ligand", "ligand_fit_output")`; existing inheritance machinery handles the rest |

#### Tests

`tst_ligand_workflow_restart.py` (14 tests, suite S22):

- **Section A** (8 tests): `past_analysis` check — verifies each
  downstream flag (`refine_done`, `ligandfit_done`,
  `phaser_done`, etc.) correctly advances past analyze; verifies
  fresh cycle 1 still routes to analyze; reproduces the exact
  nsf-d2-ligand cycle 3 context and confirms it returns
  `combine_ligand`.
- **Section B** (2 tests): regression — Phase 6b sequence-only
  precedence still fires; fresh data+model upload still goes to
  analyze first.
- **Section C** (4 tests): `ligand_fit_output` metric
  inheritance — ligand_fit_1.pdb becomes best with inherited
  R-free; existing `with_ligand` inheritance still works;
  explicit metrics on ligand_fit_output are respected (no
  inheritance); end-to-end refine→ligandfit→refine tracks correctly.

Uses the same libtbx-stub pattern as `tst_sequence_only_routing.py`
(Phase 6b tests), so workflow_engine routing tests run without a
real PHENIX install.

#### Verification

| Check | Result |
|-------|--------|
| Python syntax (both patched files) | Parses |
| `tst_ligand_workflow_restart.py` (14 new tests) | 14/14 pass |
| `tst_sequence_only_routing.py` (Phase 6b, 10 tests) | 10/10 pass (regression) |
| `tst_data_input_filter.py` (Phase 4b, 22 tests) | 22/22 pass (regression) |
| Empirical bug reproduction | Pre-fix: ligand_fit_1.pdb does NOT become best (score 122 stays on refined). Post-fix: ligand_fit_1.pdb becomes best (score 127 = stage 105 + inherited R-free 22) |
| nsf-d2-ligand cycle 3 routing | Returns `combine_ligand` (correct flow through pdbtools) |

#### Workflow Shape Change

A behavioral change worth flagging: the workflow shape is
different after this fix.

| Before | After |
|--------|-------|
| refine → ligandfit → refine → refine → refine (3 refines without progress until guard fires) | refine → ligandfit → combine_ligand (pdbtools) → refine (canonical flow) |

Both reach acceptable R-free, but the post-fix flow is the
canonical sequence with the combine_ligand step properly
executing. Tutorial expectations that encode the old shape
("expect 3 refine cycles after ligandfit") would need updating.

#### Self-Review Note

Both fixes are non-novel — each mirrors a proven pattern already
in the codebase:

- Fix 1 mirrors `_detect_cryoem_step`'s `past_analysis` exactly.
  The cryo-EM version has a comment explaining why the check
  exists: *"This handles tutorials that skip mtriage and also
  prevents the workflow from getting stuck in 'analyze' if
  mtriage's done flag wasn't set."* The X-ray version simply
  lacked the equivalent.
- Fix 2 extends an existing two-line metric-inheritance pattern
  to one additional stage. A one-character-ish change.

A new principle surfaced from this round: **when fixing
analogous bugs in parallel code paths (X-ray vs cryo-EM,
with_ligand vs ligand_fit_output), check whether the fix
pattern already exists nearby**. Both bugs had close analogs
that simply hadn't been extended to cover the affected case.

The highest residual risk is that this delivery was tested at
source level but not via a live tutorial run. Tom should re-run
the nsf-d2-ligand tutorial to confirm: (a) R-free still reaches
~0.21 or better, (b) the reasoning report no longer shows "first
refinement step" 4 times, (c) no tutorial expectations break.

### Post-Phase-5 Addition: Tier 1 Follow-Up Tests (S23 + S24)

#### Summary

Two new test suites that close verification gaps flagged in the
v116.10 review. Both are pure tests (no production code changes);
they verify existing behavior that previously had only indirect
coverage.

#### Tests Added

**S23 — `tst_initialize_plan_smoke.py` (9 tests)**

`_initialize_plan_inner` had decision-tree traces in
`tst_dock_and_stop.py` and `tst_standalone_consistency.py`, but
neither test actually CALLED the function. A future refactor
could change side effects (directive rewrites, stop_after
clearing, plan-generation gating) without tripping any existing
test.

The new tests assert on the function's observable side effects
for 9 scenarios covering each branch of the decision tree:

| Case | Scenario | Expected outcome |
|------|----------|------------------|
| A | standalone (non-preprocessing) + task | single_program_skip; no rewrite |
| B | preprocessing + task (no explicit stop) | rewrite intent; clear stop_after; proceed |
| C | preprocessing + task + "and stop" advice | preprocessing_explicit_stop_skip; preserve directive |
| D | needs_plan + task | rewrite intent; KEEP stop_after; proceed |
| E | v116.10 elif full-plan target | rewrite intent; preserve stop_after; proceed |
| F | task intent + no stop_after | single_program_skip (else branch) |
| G | solve intent + no stop_after | proceed_to_generate_plan |
| H | Phase 3d dock_in_map + task | single_program_skip (post-Phase-3d) |
| Extra | xtriage + explicit_stop | preprocessing_explicit_stop_skip (xtriage is BOTH preprocessing and standalone — preprocessing branch fires first) |

The "extra" case surfaced from a test development misclassification
that was caught by the actual control-flow trace — exactly the
kind of subtlety this layer of testing is meant to catch.

The harness parses `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS`
from source (same pattern as `tst_dock_and_stop.py`), then replays
the decision tree in a minimal harness. The full function can't
be imported because of PHENIX dependencies; the harness mirrors
the control flow exactly so a refactor that changes side effects
requires the test to be updated in parallel.

**S24 — `tst_phase3d_motivating_tutorial.py` (3 tests)**

Closes the **highest-priority verification gap** from the v116.10
review: Phase 3d's behavior change was verified only by
decision-tree traces; no tutorial exercised the dock-and-stop
path with sequence + map.

Three layers:

1. **Decision-tree** (unconditional): asserts dock_in_map + task
   routes to `single_program_skip` (Phase 3d post-fix).
2. **YAML expectations format** (unconditional):
   `get_expectations_yaml_entry()` returns the snippet to add to
   `tutorial_expectations.yaml`; test verifies it parses and
   contains the expected program ordering.
3. **Session-based** (skip-aware): given a real session.json
   fixture (via `DOCK_AND_STOP_FIXTURE` env var or
   `tests/fixtures/dock_and_stop_session.json`), asserts on
   actual program ordering and the absence of the anti-pattern
   `"no model file found"`.

A `docs/PHASE_3D_TUTORIAL_README.md` walks through wiring up the
synthetic tutorial data (sequence + map) and the session fixture.

#### Files Modified

| File | Lines |
|------|-------|
| `tests/tst_initialize_plan_smoke.py` (new) | 631 |
| `tests/tst_phase3d_motivating_tutorial.py` (new) | 367 |
| `tests/run_all_tests.py` (S23 + S24 registered) | +58 |
| `docs/PHASE_3D_TUTORIAL_README.md` (new) | ~100 |

#### Verification

| Check | Result |
|-------|--------|
| Both new test files: Python syntax | Parses |
| `tst_initialize_plan_smoke.py` | 9/9 pass |
| `tst_phase3d_motivating_tutorial.py` | 3/3 pass (one skip-aware) |
| `tst_dock_and_stop.py` (existing, regression) | 5/5 pass |

#### Tier 1 status

With S23 and S24 shipped, all four Tier 1 items from the
v116.10 follow-up plan are now complete:

| Item | Status |
|------|--------|
| 1.1 Phase 5 docs amendment for CC key fix | Done (covered by prior CHANGELOG amendments) |
| 1.2 Multi-clause stop test cases | Handled by Tom directly |
| 1.3 Phase 3a behavioral smoke test | Done (S23) |
| 1.4 Phase 3d motivating-case tutorial | Done (S24) |

### Post-Phase-5 Addition: GUI Tutorial Override Fix

#### Summary

When running a tutorial whose README mentions phenix programs
not yet available in the AI Agent (e.g., `groel-dock-refine`
mentions a currently-unsupported program), the GUI was blocking
Run with a `Sorry` dialog even after the user provided their own
files and advice that overrode the README.

#### Bug Details

| Symptom | "This tutorial requires programs not yet available" dialog blocks Run, even when user has provided their own files and instructions |
| Root Cause | `AIAgent.py:validate_params()` has two consecutive checks. The first correctly distinguishes user-provided input from README fallback. The second (the tutorial-blocking block) fires unconditionally based on `tut.can_run` regardless of whether the user was using the README at all. |
| Fix | One condition added to the tutorial-blocking guard: only fire when `not has_files and not has_advice`. The banner at the top of the page still shows the informational "tutorial requires X" warning — that's fine and was kept. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `wxGUI2/Programs/AIAgent.py` | 69-86 | Added `and not has_files and not has_advice` to the blocking guard |

#### Scenario Verification

| Scenario | files | advice | tut.can_run | Pre-fix | Post-fix |
|----------|-------|--------|-------------|---------|----------|
| Tutorial, no own input | No | No | False | BLOCKED ✓ | BLOCKED ✓ |
| Tutorial + own files (bug case) | Yes | No | False | BLOCKED (bug) | PROCEEDS ✓ |
| Tutorial + own advice | No | Yes | False | BLOCKED (bug) | PROCEEDS ✓ |
| Tutorial + both | Yes | Yes | False | BLOCKED (bug) | PROCEEDS ✓ |
| Non-tutorial, no input | No | No | n/a | First block raises | First block raises ✓ |
| Tutorial with supported programs | any | any | True | PROCEEDS ✓ | PROCEEDS ✓ |

Three pre-fix bug scenarios fixed; no correct behavior changed.

#### Verification

GUI code (wxPython, requires a display to run). Verified by:
- Python syntax check
- Manual scenario tracing through 7 user paths

User-confirmed working: Tom retested groel-dock-refine after
deploying the fix and the dialog no longer blocks Run.

#### Self-Review Note

This fix was discovered during integration testing of Phase 6c
(ligand workflow restart), not from a planned cycle item. The
banner display (informational) was deliberately left in place so
users still see what the README is asking for; only the
blocking dialog was conditioned on user-override.

### Post-Phase-5 Addition: Declarative Program Requirements (Tier 2.1, S25)

#### Summary

Added an optional `requirements:` block to entries in
`programs.yaml`. When present, the block is evaluated by a new
filter pass inside `_filter_programs_missing_data_inputs`.
Programs without the block are unaffected. Initial scope: one
declaration (`phenix.autobuild`) to plug a known gap where the
program could be picked and crash at runtime with `MTZ lacks
phase/FOM columns`.

This is the result of Tier 2.1 design work, reviewed externally
by Gemini and refined in the `PATH_Y_DESIGN.md` document.

#### Bug Details

| Symptom | LLM picks `phenix.autobuild` in sessions without phase information; command runs and fails with `MTZ lacks phase/FOM columns` after building setup. One wasted cycle plus a confusing runtime error. |
| Root Cause | Autobuild had an explanation path (`explain_unavailable_program`) telling the LLM why it shouldn't pick autobuild, but no actual filter. The LLM could ignore the explanation. The state machine alone didn't gate autobuild because the relevant step has multiple valid programs and the conditions in `workflows.yaml` weren't expressive enough. |
| Fix | Declarative `requirements:` block on `phenix.autobuild` evaluated by new `_check_requirements()` parser. Filter integrates inside `_filter_programs_missing_data_inputs` after the existing xtriage/mtriage checks. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `agent/workflow_engine.py` | +135 | New method `_check_requirements()` (~100 lines including grammar docs); new clause-keyword set `_KNOWN_REQUIREMENT_CLAUSES`; integration block inside `_filter_programs_missing_data_inputs` (~22 lines) |
| `knowledge/programs.yaml` | +4 | `requirements:` block on `phenix.autobuild` entry |
| `tests/run_all_tests.py` | +22 | Register S25 (Program Requirements) |

#### Declarative Schema (v1)

```yaml
phenix.<program>:
  # ... existing fields ...
  requirements:    # NEW — optional
    requires:
      - <clause>
      - <clause>
```

Grammar (closed, boolean-only):

| Clause | Semantics |
|--------|-----------|
| `has: <name>` | `context["has_<name>"]` truthy |
| `has_any: [<name>, ...]` | At least one `context["has_<n>"]` truthy |
| `not_has: <name>` | Not `context["has_<name>"]` |
| `done: <name>` | `context["<name>_done"]` truthy |
| `not_done: <name>` | Not `context["<name>_done"]` |

The grammar is intentionally closed. No `if/then`, no nested
logic, no metric comparisons. Per Pitfall 10 of the design
doc: if a requirement needs control flow, it goes in Python
(Mechanism 3 or 4), not here.

#### autobuild Declaration

```yaml
phenix.autobuild:
  requirements:
    requires:
      - has_any: [data_mtz, phased_data_mtz]
      - has_any: [model, placed_model, phased_data_mtz]
```

`has_phased_data_mtz` appears in both clauses because a phased
MTZ satisfies both "x-ray data" and "phase information"
independently. So a phased MTZ alone is sufficient. Other
qualifying scenarios: raw data + model, raw data + placed
model (after phaser). Filtered scenarios: data alone, model
alone, empty session.

#### Tests (30 new test cases)

| Section | Tests | Description |
|---------|-------|-------------|
| A. Parser unit tests | 17 | Each clause type evaluated in isolation; AND of clauses; malformed input handling; unknown-keyword warning + strict-mode raise |
| B. Filter integration | 6 | autobuild scenarios: phased alone (keep), data+model (keep), data+placed_model (keep), data alone (filter), model alone (filter), empty session (filter) |
| C. Backward compat | 7 | Programs without `requirements:` unaffected; existing xtriage/mtriage filters still work; order of operations correct; mixed program lists handled |

Run with `phenix.python tests/tst_program_requirements.py`.

#### Order of Operations

The new filter integrates into the existing pipeline at step 4:

1. State machine (Mechanism 1) — decides current step
2. YAML step conditions (Mechanism 2) — filters programs within step
3. `_filter_programs_missing_data_inputs` — existing xtriage/mtriage checks
4. **`_check_requirements` — NEW declarative filter (this delivery)**
5. `_check_program_prerequisites` (Mechanism 4) — runs for directive additions

`requirements:` can only remove programs from `valid_programs`,
never add them.

#### Forward-Looking Guidance

Per the design doc, `requirements:` is **the preferred path
for new program filtering rules going forward**. The four
existing hand-coded mechanisms aren't being grown.
Documentation rules:

- Add `requirements:` ONLY when the existing mechanisms don't
  cover the filtering need (avoid duplication)
- The grammar is closed; extension requires documented
  justification per Pitfall 10
- The block is boolean-only; file-content checks (e.g., "MTZ
  has FOM columns") must be exposed as context flags elsewhere
  (per Pitfall 16)
- Before adding a new declaration, complete the context-flag
  audit checklist (per Pitfall 15)

See `PATH_Y_DESIGN.md` in the design package for the full
rationale, 16 pitfalls with mitigations, and three rejected
alternatives.

#### Verification

| Check | Result |
|-------|--------|
| Python syntax (workflow_engine.py, run_all_tests.py) | parse OK |
| New test suite (40 tests) | 40 passed, 0 failed |
| Backward-compat: existing `tst_data_input_filter` regression | xtriage/mtriage filters still pass |
| Context-flag audit for autobuild's 4 referenced flags | All set in `_build_context()` before filter runs |
| `tst_scenario_tracer.py` (S5A) | Passes after Tier 2.1.1 fix below |
| `tst_phase4_history_flags.py` (all_read_flags_are_written) | Passes after Tier 2.1.1 fix below |

#### Post-Deployment Fix (Tier 2.1.1): `any_of` clause + phantom flag

After initial deployment, two test failures surfaced that required
the Tier 2.1 implementation to be revised:

**Failure 1 — S5A in `tst_scenario_tracer.py`:**

The simulator reported `autobuild not available after autosol`. The
v1 rule required `has_any: [model, placed_model, phased_data_mtz]`
in its second clause, but in the simulator (and some real session
paths) `autosol_done=True` is set without `has_phased_data_mtz` being
updated on the same cycle. The existing explanation code at
`explain_unavailable_program` uses `phaser_done OR autosol_done OR
has_placed_model_from_history` as the gate; my v1 rule diverged from
this and was stricter.

**Failure 2 — phantom `has_` flag in `tst_phase4_history_flags.py`:**

The static flag scanner regex matches `context.get("has_X")` to
extract flag names. The v1 implementation had `context.get("has_" +
clause["has"])` patterns, which the regex extracted as just `has_`
(captured up to the first closing quote). This produced a false
"flag read but not written" warning.

**Fixes shipped in Tier 2.1.1:**

1. **Grammar extension: `any_of` clause type** (Pitfall 10 justification
   provided in workflow_engine.py docstring). Where `has_any: [...]`
   takes a list of names and auto-prefixes `has_`, `any_of: [...]`
   takes a list of sub-clauses (each is itself a clause: `has`,
   `done`, `not_has`, `has_any`, or nested `any_of`). This allows
   mixing flag families:

   ```yaml
   - any_of:
       - has: phased_data_mtz
       - done: autosol
       - done: phaser
   ```

2. **Refactor key construction in `_check_requirements`** to use
   intermediate variables (matching the pattern already used in
   `_check_conditions`). Replaces `context.get("has_" + name)` with
   `key = "has_" + name; context.get(key)` so the static scanner
   doesn't trip on a literal `"has_"` inside a `context.get(...)`.

3. **Update `phenix.autobuild` requirements**: replace the
   inline-only `has_any` second clause with an `any_of` covering
   both flag families:

   ```yaml
   requirements:
     requires:
       - has_any: [data_mtz, phased_data_mtz]
       - any_of:
           - has: phased_data_mtz
           - has: model
           - has: placed_model
           - has: placed_model_from_history
           - done: phaser
           - done: autosol
   ```

   This now exactly mirrors the existing explanation pattern.

4. **Test suite expanded from 30 → 40**: added 6 tests for the
   new `any_of` clause type (positive/negative/nested/empty/malformed/
   AND-combined), plus 4 regression tests for the S5A scenario
   (autobuild kept after autosol_done, phaser_done,
   has_placed_model_from_history; still filtered when no history at all).

## Version 115.09b (GUI Fixes + Ligand Workflow + Bug 1)

### Summary

Production fixes from GUI testing and run 25 analysis. Three categories:
explicit_program loop guard, ligand-fitting workflow (6 fixes), and
GUI-mode corrections. 7 files modified.

### Modified Files

| File | Changes |
|------|---------|
| `agent/graph_nodes.py` | Bug 1: explicit_program done-flag guard prevents LLM-driven program loops. Removed dead .sca-only proactive check (replaced by diagnosable error). |
| `agent/directive_extractor.py` | Ligand-fit placement: `model_is_placed=True` for "fit ligand" advice. Ligand-fit `after_program` clearing: removes LLM-set `after_program` for multi-step ligand workflows. |
| `agent/workflow_engine.py` | Post-ligandfit exemption: defers `after_program_done` during combine/refine steps. Combine_ligand guard: forces pdbtools-only in combine step. Debug tracing for combine_ligand routing. |
| `agent/command_builder.py` | `phaser_sad.atom_type` interception: converts to `additional_atom_types`, prevents Se→S override in autosol. |
| `knowledge/programs.yaml` | Polder: `requires_resolution` invariant with `auto_fill_resolution` + `xray_data.high_resolution` strategy flag. Pdbtools: `exclude_patterns: [pose]` for ligand input (prevents selecting pose files over final model). Merged duplicate hints blocks. |
| `knowledge/diagnosable_errors.yaml` | `missing_crystal_symmetry`: "No unit cell info available" + "Cell and/or symmetry not specified" patterns for xtriage with .sca data. |
| `phenix/programs/ai_agent.py` | Fix A: `_has_explicit_stop` regex guard on preprocessing stop override (xtriage/mtriage). GUI auto-discovery skip: in GUI mode with user-selected files, skip supplement to prevent PHENIX artifacts from changing workflow intent. |

### Bug Details

| Bug | Symptom | Root Cause | Fix |
|-----|---------|------------|-----|
| Bug 1 (explicit_program loop) | bgal_denmod 3x resolve_cryo_em, hipip-refine 3x refine, emd_6123 3x map_to_model | `explicit_program` injection at graph_nodes.py line 691 bypasses done-flag check | Check done flag before injecting; skip if program already completed |
| Fix A (preprocessing stop) | "run mtriage and stop" generates full 5-stage plan | `_preprocessing_programs` unconditionally clears `after_program` | `_has_explicit_stop` regex check before clearing |
| Fix B (p9-xtriage) | xtriage fails twice on .sca without cell | Proactive file-type check can't distinguish p9 from gene-5-mad | Reactive: diagnosable_errors.yaml entry for "No unit cell" |
| GUI discovery | nsf-d2-ligand picks up pdb_sequences.fa → predict_and_build | Auto-discovery supplements with PHENIX project artifacts | Skip supplement in GUI mode with user-selected files |
| Ligand placement | "fit ATP" advice → predict_and_build instead of ligandfit | `model_is_placed` not set; model treated as unplaced | Rules overlay sets `model_is_placed=True` for ligand-fit signals |
| Ligand after_program | Ligandfit→pdbtools→refine→polder workflow blocked at each step | LLM sets `after_program` to different program each run (ligandfit/refine/polder); each choice blocks a different step | Clear `after_program` when ligand-fit signals detected; plan template drives workflow |
| Combine_ligand stuck | `combine_ligand` step gets STOP instead of pdbtools | `_apply_directives` wipes valid_programs via `after_program_done` | Post-ligandfit exemption + combine_ligand guard forces pdbtools-only |
| Post-ligandfit refine | LLM picks molprobity instead of refine after pdbtools | `_apply_directives` re-adds ligandfit or polder to front of list | Three-way branch: combine→pdbtools-only; refine→prioritize refine; done→STOP |
| phaser_sad.atom_type | Autosol uses S instead of Se | LLM adds `phaser_sad.atom_type=S` which overrides primary atom_type | Intercept and convert to `additional_atom_types` |
| Polder resolution | Polder crashes without resolution limit | No `auto_fill_resolution` invariant for polder | Added `requires_resolution` invariant + `xray_data.high_resolution` strategy flag |
| Pdbtools pose file | Pdbtools uses `ligand_fit_1_pose_5.pdb` instead of `ligand_fit_1.pdb` | `prefer_patterns: [ligand_fit]` matches both final and pose files | `exclude_patterns: [pose]` on ligand input slot |

### Verification

| Fix | Status |
|-----|--------|
| Bug 1 (bgal_denmod/rules_only) | ✅ mtriage→resolve_cryo_em→predict_and_build→RSR, CC=0.68 |
| Bug 1 (hipip-refine/rules_only) | ✅ refine→xtriage→molprobity→auto-stop, R-free=0.223 |
| Fix A (xtriage and stop) | ✅ Stops after xtriage in GUI |
| Fix B (p9-xtriage) | ✅ 1 cycle, clean diagnosis with helpful message |
| GUI discovery skip | ✅ Deployed, prevents pdb_sequences.fa pickup |
| Ligand workflow | ✅ xtriage→refine→ligandfit→pdbtools→refine→polder→molprobity |
| Polder resolution | ⏳ Deployed to server, awaiting verification |
| Pdbtools pose exclusion | ⏳ Deployed, awaiting verification |


## Version 115.09 (Tutorial Routing Fixes)

### Summary

Four routing bugs that prevented tutorials from completing. Fixes span
cryo-EM workflow progression, data-limitation detection, validation-only
routing, and MR-SAD intent recognition. Intent detection lives in the
directive extractor layer (LLM prompt + rules-based fallback); the
routing engine stays deterministic.

### Modified Files (5 production files)

| File | Fixes | Changes |
|------|-------|---------|
| `agent/workflow_engine.py` | 1, 3, 4 | Fix 1: `map_sharpening_done`, `map_symmetry_done`, `has_optimized_full_map` added to `past_analysis` gate in `_detect_cryoem_step`. Fix 3: `wants_validation_only` read from directives in `build_context`; validation shortcut in `_detect_xray_step` (guarded by `has_model` AND `has_data_mtz`). Fix 4: `force_mr` flag set in `build_context` when `use_mr_sad` + `has_model` + `not has_search_model` + `not phaser_done`; routing override in `_detect_xray_step` before placement probes; MR-SAD guard in `get_valid_programs` updated to check `force_mr`. |
| `agent/graph_nodes.py` | 2 | `.sca/.hkl`-only data detection in `perceive()` with abort_message. Guarded by `not _has_model` AND `not _has_sequence` AND no `unit_cell` in directives. Follows established early-stop pattern (stop + stop_reason + abort_message + intent dict). |
| `agent/directive_extractor.py` | 3, 4 | LLM prompt: `wants_validation_only` added to `workflow_preferences` schema with extraction guidance. MR-SAD prompt strengthened with concrete examples. `_validate_directives`: `wants_validation_only` added to allowed boolean keys. `extract_directives_simple`: rules-based fallback for 6 validation signals and 6 MR-SAD patterns. |
| `agent/plan_generator.py` | 3 | `wants_validation_only` added to `_build_context` initial dict and propagated from `workflow_preferences` to template selection context. |
| `knowledge/plan_templates.yaml` | 3 | `validate_existing` template: `applicable_when` requires xray + `wants_validation_only` + `has_search_model` + `model_is_placed` (priority 60, beats `refine_placed_ligand` at 55). Two stages: data_assessment (xtriage) + validation (model_vs_data, molprobity). |

### Bug Details

| Bug | Tutorial | Symptom | Root Cause | Fix |
|-----|----------|---------|------------|-----|
| 1 | bgal_denmod, apoferritin_denmod, ion_channel_denmod | Cryo-EM stops after 1 cycle (map_sharpening or resolve_cryo_em) | `map_sharpening_done` missing from `past_analysis` gate → `detect_step` returns "analyze" → no programs → STOP | 3 new flags in `past_analysis` including file-presence fallback |
| 2 | p9-xtriage | xtriage fails on unmerged .sca without cell dimensions | Data limitation: unmerged scalepack needs unit cell + space group provided interactively | Early detection in `perceive()` with helpful abort message |
| 3 | pka-validate | Runs predict_and_build instead of molprobity despite "Analysis only" advice | No validation-only routing path; directive extractor doesn't capture validation intent | `wants_validation_only` directive + validation shortcut + plan template |
| 4 | lysozyme-MRSAD | R-free stuck at 0.54; phaser never runs despite "MR-SAD" advice | PDB classified as `model` not `search_model`; LLM ignores `use_mr_sad` prompt; no rules-based MR-SAD fallback | `force_mr` flag + MR-SAD rules fallback + prompt strengthening |

### Review Fixes (2 bugs found during review)

| Issue | Fix |
|-------|-----|
| Fix 2 false positive: gene-5-mad (merged .sca + sequence) falsely aborted | Added `not _has_sequence` guard to `.sca-only` condition |
| Fix 3 false positive: "analysis only" matched xtriage-only and cryo-EM tutorials | Removed from rules-based signals; LLM prompt handles the nuance |

### Deployment Fix (discovered via run 19 verification)

Run 19 (OpenAI, post-deployment) showed all four fixes not firing: directives
contained no `wants_validation_only` or `use_mr_sad` despite code being
deployed and bytecache cleared.

**Root cause**: The rules-based intent patterns were only in
`extract_directives_simple()` — the fallback path. OpenAI runs use the
LLM path in `extract_directives()`, which calls the LLM, gets directives
without these flags, and returns. The rules-based patterns never execute.

**Fix**: Extracted patterns into `_apply_workflow_intent_fallback()` shared
helper. Called as **post-LLM overlay** in `extract_directives()` (after LLM
returns, before final return) AND from `extract_directives_simple()`.
Rules always run last and always win for routing flags. The LLM cannot
override them because it never sets them (0/240 extractions).

**Future**: Replace overlay with centralized `_DIRECTIVE_SCHEMA` and
registry-driven `_merge_tiered` merge (v115.10, see `docs/directive_merge_plan.md`).

### Production Verification Fixes (runs 20-22)

| Issue | Discovery | Root Cause | Fix |
|-------|-----------|------------|-----|
| Fix 3: `_is_valid_file` rejects valid 3dnd.pdb | Run 21 `[GATE]` output: `has_model=False` despite 3225 ATOM records | PDB scan limit of 500 lines; 3dnd.pdb has 546 header lines before first ATOM record | `workflow_state.py`: increased scan limit to 2000 |
| Fix 3: `has_data_mtz=False` for phased MTZ | Run 20 analysis: completed structures with phase columns go to `phased_data_mtz` | Validation shortcut only checked `has_data_mtz` | `workflow_engine.py`: added `has_phased_data_mtz` to context + validation shortcut |
| Fix 1: `map_sharpening_done` regex wrong | Diagnostic: regex matches "sharpened" but output is `auto_sharpen_A.ccp4` | Zombie check table regex too narrow | `workflow_state.py`: broadened regex from "sharpened" to "sharpen" |
| Fix 2: `.sca-only` check blocked by sequence guard | Run 19: p9-xtriage has `seq.dat` like gene-5-mad | Proactive file-type approach can't distinguish the two | **Deferred** — needs reactive approach (diagnosable_errors.yaml) |

### Additional files modified (production verification)

| File | Changes |
|------|---------|
| `agent/workflow_state.py` | `_is_valid_file`: PDB scan limit 500→2000. `_ZOMBIE_CHECK_TABLE`: map_sharpening regex "sharpened"→"sharpen". |
| `agent/workflow_engine.py` | `has_phased_data_mtz` added to `build_context` initial dict + validation shortcut condition. |

### Known Issues (deferred to v115.10)

| Issue | Description |
|-------|-------------|
| Preprocessing stop override | `ai_agent.py` line 2761 unconditionally clears `after_program` for xtriage/mtriage even when user explicitly says "stop". Fix: check `_has_explicit_stop` regex before clearing. |
| CIF model categorization | `model.cif` files go to generic `cif` category, not `model`. Affects model-dependent routing for CIF-format structures. |
| Fix 2 redesign | p9-xtriage `.sca-only` proactive check fails (has sequence like gene-5-mad). Need reactive: xtriage "no unit cell" failure as diagnosable terminal error. |
| Directive extraction refactor | Replace post-LLM overlay with `_DIRECTIVE_SCHEMA` registry + `_merge_tiered`. See `docs/directive_merge_plan.md`. |


## Version 115.08 (Phased File Detection + Systematic Testing Framework)

### Summary

Four critical fixes for phased file detection (F1–F4) that caused rab3a-refine
and nsf-d2-refine to fail across all 5 modes. One additional bug fix (B1:
last_program) found by the testing framework. A 10-phase systematic testing
framework exercising file categorization, routing, command building, error
classification, and LLM resilience across 32 tutorials.

### Modified Files (main code — 3 files)

| File | Change | Details |
|------|--------|---------|
| `agent/workflow_state.py` | F1–F4 | Content-based phased detection replacing filename markers. New: `_PHASE_COLUMN_CACHE` (iotbx column-type check cached by abspath+mtime), `_ascii_phase_heuristic()` (ASCII header check for .hkl files), `_has_phase_columns_cached()` (cached wrapper). `_resolve_phased_promotions()` moved from hardcoded-only to shared post-processing (F2). Category-exclusivity enforcement — files appear in one of `data_mtz` or `phased_data_mtz`, not both (F3). Conditional removal from `data_mtz` only when other non-phased data files remain (F4). `_IGNORED_FORMATS` dict for .cv/.xplor/.param/.eff extensions. |
| `agent/workflow_engine.py` | B1, diagnostics | B1 fix: added `"last_program": history_info.get("last_program")` to `build_context()` dict literal — the `after_program` directive check at line 1799 previously always saw None. `[GATE]` diagnostic logging controlled by `PHENIX_AGENT_DIAG_VALID_PROGRAMS=1` env var — prints context flags, detect_step result, and STOP reasons to aid routing debugging. |
| `agent/graph_nodes.py` | WARNING event | WARNING event emitted for each ignored format file in `perceive()` — informs user that .cv/.xplor files are not supported for auto-fill. |

### New Test Files (12 files)

| File | Phase | Tests | Purpose |
|------|-------|-------|---------|
| `tst_file_categorization.py` | — | 42 | Unit tests for v115.08 phased detection fixes |
| `tst_phase0_static_audit.py` | S0 | 5 | Parse check, bare except scan, import fallback check |
| `tst_phase1_contract_gaps.py` | S1 | 128 | AST-based coverage map of 4 key modules |
| `tst_phase2_path_consistency.py` | S2 | 10 | YAML vs hardcoded categorization path diff |
| `tst_phase3_serialization_symmetry.py` | S3 | 28 | JSON symmetry + invariants + AgentSession round-trip |
| `tst_phase4_history_flags.py` | S4 | 8 | Flag writer/reader consistency, B1 fix verification |
| `tst_phase5_error_classification.py` | S8 | 7 | 3 error classifiers, pattern overlap, self-detection |
| `tst_phase6_category_consumer.py` | S5 | 3 | input_priorities + fallback_categories alignment |
| `tst_phase7_routing_simulation.py` | S6 | 32 | 3-cycle routing simulation with real production code |
| `tst_phase8_command_building.py` | S7 | 15 | CommandBuilder.build() for 7 tutorials + 2 edge cases |
| `tst_phase9_llm_perturbation.py` | S9 | 17 | Filename/program/parameter/truncation/empty perturbation |
| `run_all_tests.py` | — | — | Updated with S0–S9 registration |

### New Documentation

| File | Purpose |
|------|---------|
| `docs/PHASE_REVIEW_REPORT.md` | Detailed report of 54 review fixes across all 10 phases |

### Bug Details

| Bug | Severity | Description | Fix |
|-----|----------|-------------|-----|
| F1 | CRITICAL | Filename-marker detection checked 'phased' but files were named 'phases' | Content-based detection (iotbx + ASCII heuristic) |
| F2 | CRITICAL | `_resolve_phased_promotions()` only ran on hardcoded path, not YAML (production) | Moved to shared post-processing |
| F3 | CRITICAL | YAML put files in both `data_mtz` and `phased_data_mtz` | Exclusivity enforcement |
| F4 | CRITICAL | Removing from `data_mtz` could delete the only data file | Conditional removal (only when alternatives exist) |
| B1 | LOW | `last_program` never transferred from `_analyze_history()` to `build_context()` | 1-line addition to `build_context()` dict |

### Test Framework Review (54 fixes across 9 phases)

The testing framework itself was reviewed across multiple rounds. Key
classes of defects found and fixed in the test scripts:

- **5 tests that could never fail** (Phases 2, 6, 7, 9) — status was always
  PASS/PARTIAL but raise checked for FAIL. Fixed with whitelist patterns.
- **3 tautological assertions** (Phases 3, 9) — `assert True`, `assert X or True`.
- **6 phantom/wrong names** in boundary function sets (Phase 1).
- **5 missing exception isolation** patterns (Phases 2, 6, 7, 8) — one crash
  killed all remaining tutorials.
- **5 wrong `best_files` keys** (Phases 8, 9) — `'data'` should be `'data_mtz'`.
- **9 new AgentSession round-trip tests** (Phase 3) — replaced 18 tautological
  JSON identity tests with real production code exercise.

See `docs/PHASE_REVIEW_REPORT.md` for full details.


## Version 115.07 (Run 15b Bug Fixes — Phase 3)

### Summary

Four bugs identified from analysis of run 15b (OpenAI, 371 runs, 42
tutorials). Fixes address a graph crash from uncoerced JSON metrics, a
cryo-EM half-map pair detection gap, a terminal crash loop, and an LLM
hallucination pattern. Two issues deferred as known limitations.

**Run 15b headline:** 9→18 GOOD tutorials (doubled), failure rate 39%→34%.
Phase 3 targets remaining failures: lysozyme-refine (77%), apoferritin
denmod_dock (57%), 7rpq (never refines), AF_7n8i (5/10 fail).

### Modified files (4)

| File | Bug | Changes |
|------|-----|---------|
| `agent/structure_model.py` | 4 | **Numeric coercion**: added `_coerce_numerics()` function — coerces all float/int fields after `from_dict()` deserialization. Added `_safe_float()` guards at all 6 arithmetic and formatting sites in `_detect_problems()`, `get_summary()`, and `_format_report()`. Progress entries also coerced. Prevents `TypeError: unsupported operand type(s) for -: 'str' and 'float'` when model_vs_data stores metrics as strings in session JSON. |
| `agent/workflow_state.py` | 6 | **Two-tier half-map pair detection**: relaxed `_categorize_files()` heuristic. Tier 1 (new): when exactly 2 `full_map` files form a `_1/_2` pair, promote to `half_map` without requiring a companion full map. Tier 2 (existing): when ≥3 `full_map` files include a pair, promote only if a companion remains. Fixes AF_7n8i where `box_1.ccp4`/`box_2.ccp4` were the only maps. |
| `agent/phil_validator.py` | 3 | **Blocked params**: added `phenix.resolve_cryo_em` to `_BLOCKED_PARAMS` with 4 entries: `mask_atoms`, `output.prefix`, `d_min`, `main.number_of_macro_cycles`. These are LLM-hallucinated params copied from mtriage log output that cause RuntimeError. |
| `knowledge/diagnosable_errors.yaml` | 1 | **Terminal diagnosis**: added `unknown_chemical_element` — detects "Unknown chemical element type" errors (bad PDB columns 77-78) and stops immediately with actionable hint instead of crash-looping. |

Also modified (config):

| File | Bug | Changes |
|------|-----|---------|
| `knowledge/programs.yaml` | 3 | Removed `mask_atoms` from `phenix.resolve_cryo_em` `strategy_flags` whitelist. Only `resolution` passes through. |

### Bug details

| Bug | Tutorial | Symptom | Root cause | Fix |
|-----|----------|---------|-----------|-----|
| 1 | lysozyme-refine | 77% fail, no metric | PDB has AU atom missing element columns → refine crashes → agent retries 5× | Terminal diagnosis: stop with hint |
| 3 | apoferritin_denmod_dock | 57% fail (LLM modes) | LLM copies `mask_atoms=True` from mtriage log → `strategy.mask_atoms_atom_radius="True"` RuntimeError | Whitelist removes it; `_BLOCKED_PARAMS` double-blocks it |
| 4 | 7rpq_AF_reference | 0% fail but stops at R=0.386 | model_vs_data stores `r_work="0.385"` (string) → cycle 3 graph crashes on `r_free - r_work` | `_safe_float()` at 5 sites in metric_evaluator (true root cause: USE_YAML_METRICS=True); belt-and-suspenders in metrics_analyzer (7 sites), structure_model, graph_nodes, event_formatter, kb_tags, workflow_state._analyze_history |
| 5 | lowres_restraints | 75% fail (restraint params stripped) | PHIL validator strips `reference_model.file`, `ncs.*`, `secondary_structure.*`, `ramachandran_restraints` — 7 of 8 params lost | Hierarchical prefix whitelist + path resolution + reference model exclusion |
| 6 | AF_7n8i | 5/10 modes fail at mtriage | `box_1.ccp4`/`box_2.ccp4` not recognized as half-maps (old heuristic required ≥3 full_map files) | Tier 1: exactly 2 matching files → promote |
| 7 | 1aba-polder | polder fails (missing: model) | Protein model with ligand atoms misclassified as `ligand_pdb` → polder can't find `model` slot | File-size fallback: PDB >10KB in ligand_pdb rescued to model |
| 8 | apoferritin_denmod_dock (rules_only) | dock_in_map loops 3x, RSR never offered | `emd-20026_auto_sharpen_A.ccp4` excluded from `full_map` by `*_a.*` pattern → `has_full_map=False` → RSR blocked | Orphan-map promotion: map files not in any subcategory → promoted to full_map |

### Bug 5 details — Reference model restraints for phenix.refine

**Problem:** phenix.refine has a massive PHIL parameter tree, but the
strategy_flags whitelist had only 13 entries. Advanced crystallographic
restraint parameters (reference model, NCS, secondary structure, Ramachandran)
were correctly extracted by the LLM from the README but stripped by PHIL
validation. The tutorial requires all of these to produce a meaningful result.

**Fix A — Hierarchical prefix whitelist** (`programs.yaml` + `phil_validator.py`):
Added `allowed_phil_prefixes` to phenix.refine: `reference_model.`,
`secondary_structure`, `ncs.`, `ramachandran`. The PHIL validator now allows
any strategy key containing one of these as a case-insensitive substring.
This covers both short forms (`ncs.type`) and full PHIL paths
(`refinement.pdb_interpretation.ncs.type`). Added `ramachandran_restraints`
as an individual strategy_flag (standalone param, not under a namespace).

**Fix B — Path resolution** (`program_registry.py`):
Added a pre-pass in `build_command()` that detects strategy values ending
in file extensions (`.pdb`, `.cif`, `.params`, `.eff`, etc.) and resolves
them to absolute paths. Builds a basename→path lookup from the command's
input files and the working directory. Handles the LLM writing
`reference_model.file=4pf4.pdb` (relative) by resolving to the full path.

**Fix C — Reference model exclusion** (3-tier):
- Tier 1 (`command_builder.py`): When `reference_model.file` is in the strategy,
  exclude that file from primary model selection. Already existed (lines 380-394),
  but was blocked because Fix A wasn't in place — now unblocked.
- Tier 2 (`workflow_state.py`): When ≥2 PDB files are in the `model` category
  and one matches `reference|homolog|template|restraint|high.res` (but NOT
  agent output prefixes like `refine_`, `autobuild_`), reclassify to
  `reference_model`.
- Tier 3: LLM extraction of `reference_model.file=` from README feeds Tier 1.

**Fix D — Strategy rewrites** (`graph_nodes.py`):
5 new `_STRATEGY_REWRITES` entries normalize verbose PHIL paths to shortest
valid forms (e.g. `refinement.pdb_interpretation.ncs.type` → `ncs.type`).

### Deferred issues

| Tutorial | Symptom | Why deferred |
|----------|---------|-------------|
| lysozyme-MRSAD | 72% fail (multi-wavelength label confusion) | Needs wavelength selection feature (Phase 4) |

### Other changes this session

| File | Change |
|------|--------|
| `agent/metrics_analyzer.py` | Added `_safe_float()` at all 7 numeric read sites in `derive_metrics_from_history()`, `_analyze_xray_trend()`, `_analyze_cryoem_trend()`, `get_latest_resolution()`, `get_best_r_free()`, `get_latest_r_free()`, `get_latest_map_cc()`. Root cause of Bug 4 crash: JSON round-tripping turns floats to strings. |
| `agent/metric_evaluator.py` | Added `_safe_float()` at 5 arithmetic sites: `_analyze_xray_trend()` r_free extraction, `_analyze_cryoem_trend()` CC extraction, `is_significant_improvement()`, `calculate_improvement_rate()`, `is_plateau()`. TRUE root cause of Bug 4: `USE_YAML_METRICS=True` routes through this file, not `metrics_analyzer.py`. Added import fallback for standalone testing. |
| `agent/event_formatter.py` | Wrapped compact metrics formatting (L920-935) with float coercion + try/except. Prevents `str - float` crash in METRICS_EXTRACTED event rendering. |
| `agent/kb_tags.py` | `_trend_tags()` now coerces all R-free trend values to float via try/except before arithmetic (diffs, total_drop). |
| `agent/workflow_state.py` | `_analyze_history()` L1841-1850: all metric reads coerced via local `_sf()` helper. Orphan-map promotion: map files in parent `map` but not in any subcategory (`full_map`, `half_map`, `optimized_full_map`) promoted to `full_map`. Fixes apoferritin_denmod_dock rules_only. |
| `agent/graph_nodes.py` | Removed stale `strategy.mask_atoms` → `mask_atoms` rewrite for resolve_cryo_em (now a no-op since mask_atoms is in `_BLOCKED_PARAMS`). |
| `knowledge/plan_schema.py` | `record_stage_cycle()` catch-up: when agent runs ahead of plan tracker (e.g. ligandfit during refine stage), advance through intermediate stages. Overshoot guard: verify `new_curr` matches program before counting. |
| `phenix_ai/local_agent.py` | Added `client_version=self._get_client_version()` + `_get_client_version()` method for parity with RemoteAgent. |
| `phenix_ai/remote_agent.py` | Added `request["settings"]["verbosity"]` and `events=parsed.get("events", [])` for parity with LocalAgent. |
| `tests/tst_phase3_bug5.py` | 38 test functions / 104 assertions covering all Phase 3 + Bug 5 + Bug 8 fixes, metric_evaluator coercion, orphan-map promotion, kb_tags string handling. |
| `tests/tst_phil_validation.py` | Fixed stale `test_rewrite_resolve_cryo_em_mask_atoms` (mask_atoms now blocked, not allowed). |
| `tests/tst_audit_fixes.py` | 3 stale tests updated: `test_k2_mtriage` (prefers_half_maps), `test_k2_map_sharpening` (positional), `test_s5h_inject_program_defaults` (generate not in defaults). |
| `tests/tst_backward_compat.py` | 3 new parity tests (26 total): `test_local_remote_agent_settings_parity`, `test_local_remote_agent_return_parity`, `test_local_agent_full_roundtrip`. |
| `docs/ARCHITECTURE.md` | Updated PHIL validation section (prefix whitelist + path resolution), added half-map Tier 1, reference model categorizer, `_coerce_numerics`, plan_schema catch-up documentation. |
| `docs/DEVELOPER_GUIDE.md` | RULE 4 expanded (prefix whitelist guidance), RULE 9 added (numeric type safety in StructureModel). |

### Test results

| Suite | Tests |
|-------|-------|
| `tst_structure_model.py` | 77/77 |
| `tst_plan_schema.py` | 53/53 |
| `tst_backward_compat.py` | 26/26 |
| `tst_command_builder.py` | 22/22 |
| `tst_event_system.py` | 13/13 |
| `tst_phase3_bug5.py` | 104/104 |
| `tst_phil_validation.py` | 15/15 |
| **Total** | **310** |
| **Total** | **191/191** |


## Version 115.06 (Cryo-EM Half-Map Handling — Phase 2)

### Summary

Three interconnected fixes addressing the #1 failure source in cryo-EM
tutorials: programs receiving only 1 half-map instead of 2. Root cause was
that the LLM assigns one half-map to a `multiple:true` slot, marking it as
"filled", so auto-fill skips it. Affected mtriage (109 failures),
resolve_cryo_em (49), and map_sharpening (46) — totaling ~200 wasted cycles.

### Modified files (3)

| File | Changes |
|------|---------|
| `agent/command_builder.py` | **Supplement logic**: after LLM files + auto-fill, a new loop checks every `multiple:true` slot. If the LLM filled it with 1 file but the category has 2+, missing files are backfilled using `_find_file_for_slot()`. Uses `os.path.realpath()` for symlink robustness. |
| `knowledge/programs.yaml` | **mtriage**: changed from `keep_half_maps_with_full_map: true` to `prefers_half_maps: true` — drops full_map when half-maps present (eliminates "Maps have different dimensions" error from mismatched grids). **map_sharpening**: added `prefers_half_maps: true` + changed half_map flag from `"half_map="` to `""` (positional args, matching the working command `phenix.map_sharpening h1.ccp4 h2.ccp4 seq_file=seq.fa`). |
| `knowledge/workflows.yaml` | **map_sharpening**: reverted all 4 entries from `has: full_map` to `has_any: [full_map, half_map]` — map_sharpening works with either input type. |

### Data flow: half-map supplement

```
LLM: {half_map=file2.ccp4}         ← one file
Corrector: → half_map_1.ccp4       ← still one file
selected_files["half_map"] = "half_map_1.ccp4"

Auto-fill: half_map already filled → SKIP

NEW SUPPLEMENT:
  half_map is multiple:true AND in selected_files
  _find_file_for_slot returns [half_map_1.ccp4, half_map_2.ccp4]
  len(2) > len(1) → supplement fires
  selected_files["half_map"] = ["half_map_1.ccp4", "half_map_2.ccp4"]

Dedup: prefers_half_maps → drop full_map, keep both half-maps
Command: phenix.mtriage half_map=h1.ccp4 half_map=h2.ccp4  ← SUCCESS
```

### Expected impact

- mtriage: gets 2 half-maps (no full_map) → no dimension mismatch
- resolve_cryo_em: gets 2 half-maps (supplemented) → "Need 2 half maps" eliminated
- map_sharpening: gets 2 half-maps positionally → "use b-factor" error eliminated
- Affects all 9 cryo-EM tutorials that were failing in run 14


## Version 115.05 (Run 13/14 Bug Fixes — Phase 1)

### Summary

Ten bugs identified from analysis of run 13 (expanded tutorial set: 33
README + 35 SOLVE across 25+ tutorials) and verified in run 14 (253 OpenAI
+ 52 Ollama runs). Fixes address R-free flag mismatches, cryo-EM half-map
dedup, map_sharpening workflow routing, docking for unplaced models, space
group parsing, PHIL parameter passthrough, and data label selection.

### Modified files (8)

| File | Bugs | Changes |
|------|------|---------|
| `knowledge/programs.yaml` | 4/5, 2/3, 7, 10 | Removed `generate=True` from refine template+defaults; added `generate_rfree_flags` strategy flag; added `prefers_half_maps: true` to resolve_cryo_em; added `fallback_categories: [model]` to dock_in_map; added `strict_strategy_flags: true` to process_predicted_model |
| `knowledge/workflows.yaml` | 1, 2/3, 7 | resolve_cryo_em: added `has: half_map` to optimize_map entry; dock_in_map: added `model` to has_any in dock_model+refine steps; added dock_in_map to cryo-EM refine step |
| `knowledge/recoverable_errors.yaml` | 9 | Added `phenix.model_vs_data` to data_label_parameters |
| `knowledge/prompts_hybrid.py` | 4/5 | Updated R-free error recovery prompt (was telling LLM "command already includes generate=True") |
| `agent/command_builder.py` | 2/3, 6 | `prefers_half_maps` dedup branch; `multiple:true` append logic with normalization; cctbx-first space group validation with prefix shortening |
| `agent/program_registry.py` | 10 | `strict_strategy_flags` check before KNOWN_PHIL_SHORT_NAMES passthrough |
| `agent/graph_nodes.py` | 7 | Unplaced model guard: redirects to dock_in_map when CC<0.10 and dock_done=False; quality floor prefers dock_in_map for unplaced models |
| `docs/RUN_13_BUG_ANALYSIS_AND_PLAN.md` | — | Analysis document |

### Bug details

| Bug | Issue | Fix | Verified |
|-----|-------|-----|----------|
| 1 | map_sharpening offered with only half-maps → crash | `has: full_map` condition (later reverted in v115.06 to `has_any`) | ✓ bgal, actin_sharpen eliminated |
| 2/3 | resolve_cryo_em gets 1 half-map (dedup drops second) | `prefers_half_maps` flag + append logic | ✓ apoferritin_dock CC=0.496 |
| 4/5 | R-free generate=True on every refinement | Removed from template; conditional via strategy flag | ✓ 3tpp R=0.167 (was CRASH) |
| 6 | Space group "P4 for refinement wit" | cctbx.sgtbx progressive prefix shortening | ✓ |
| 7 | groel: dock_in_map never offered | fallback_categories + refine step + CC guard | ✓ dock_in_map now runs |
| 8 | beta-blip hard MR case | No fix needed | — |
| 9 | "Multiple equally suitable arrays" | model_vs_data added to data_label_parameters | ✓ |
| 10 | process_predicted_model PHIL error | strict_strategy_flags blocks passthrough | ✓ b-CA-mr 0 failures |

### Ollama vs OpenAI (run 14)

Both providers use the same server-side code. Key finding: Ollama achieves
comparable quality (36% vs 39% fail rate) and dramatically outperforms OpenAI
on AF_POMGNT2 (R=0.287 vs 0.553) because Ollama's LLM correctly runs
process_predicted_model before molecular replacement.


## Version 115.04 (ASU Copy Count Tracking)

### Summary

Tracks the number of copies of the search model in the asymmetric unit (ASU)
and passes this automatically to Phaser as `search_copies` each cycle.
Copies are sourced from user directives (always wins) or xtriage log analysis.

### Modified files (5)

| File | Changes |
|------|---------|
| `programs/ai_agent.py` | `_extract_copies_from_directives()` method added; copies injected in all 3 directive extraction paths (rules-only, LLM, fallback); xtriage n_copies update from `history_record` |
| `agent/api_client.py` | `asu_copies` added to both `build_session_state()` and `build_request_v2()` whitelists so value survives client→server round-trip |
| `agent/graph_nodes.py` | `_fallback_extract_metrics()` extracts `n_copies` from xtriage log; BUILD node reads `session_info["asu_copies"]` and injects `component_copies` into strategy (skipped if LLM already set it); sanity bound 1–30 enforced |
| `phenix_ai/run_ai_agent.py` | `asu_copies` forwarded from `session_state` → `session_info` each cycle; serialized into `metadata` for persistence |
| `tests/tst_copies_tracking.py` | New test: directive path, xtriage path, priority (directive wins), range clamping, integration with `run_ai_agent` metadata |

### Data flow

**Directive path** (user says "4 copies in the ASU"):
```
_extract_copies_from_directives() → session.data["asu_copies"] = 4
  → session_info["asu_copies"] → api_client whitelists
    → run_ai_agent session_state → BUILD node → phaser.search_copies=4
```

**Xtriage path** (xtriage log: "Best guess : 4 copies"):
```
_fallback_extract_metrics() → log_analysis["n_copies"] = 4
  → history_record["asu_copies"] → session.data["asu_copies"]
    (only if not already set by directive)
      → next cycle: session_info → BUILD → phaser.search_copies=4
```

### Tests

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_copies_tracking.py` | New | Directive extraction, xtriage extraction, priority, bounds |


## Version 115.03 (Tutorial Run Bug Fixes: Bugs A–F + Anomalous Signal)

### Summary

Six bugs identified from log analysis of 16 tutorial runs across CASP7,
a2u-globulin-mr, a2u-globulin-rebuild, and advanced-xtriage.  The most
impactful fixes are: the `_has_placed_model` MR-keyword guard (Bug E)
that prevented search models from being misrouted to `xray_refined`; the
hopeless R-free retry path (Bug F) that recovered from failed MR solutions
instead of stopping; and the anomalous signal detection fix that correctly
handles weak signal (xtriage `has_anomalous=True` with `measurability` in
`[0.05, 0.10]`) without overwriting it with a measurability-based clear.
Bugs A, B, C (cascade), and D were already present in the tree.

### Modified files (4)

| File | Bug | Changes |
|------|-----|---------|
| `agent/workflow_state.py` | Anomalous | `has_anomalous` detection: `elif has_anomalous` promoted to independent `if` (was silently skipped when `anomalous_resolution` present but ≥ 6.0); `>= 0.06` weak-signal branch added; `< 0.06` clear no longer fires when `analysis["has_anomalous"]` is explicitly `True` |
| `agent/workflow_engine.py` | Bug E | `_has_placed_model`: `_wants_mr_first` guard scans constraints for MR keywords (`"molecular replacement"`, `"phaser"`, `" mr "`, etc.); if found, skips placement-keyword inference so search models are not misidentified as placed |
| `agent/workflow_engine.py` | Bug F | `_detect_xray_step` Step 2e: when `r_free >= 0.45` after MR + ≥1 refinement cycle, routes to `obtain_model` instead of falling through to `validate`/`complete` → STOP |
| `knowledge/file_categories.yaml` | Bug C (yaml) | New `predict_build_refine_internal` intermediate category matching `*overall_best*final_refine*.pdb` and `*overall_best*refine_[0-9][0-9][0-9]*.pdb`; excludes `*_final_refine_*` from `predict_and_build_output` |

### AgentState TypedDict (graph_state.py)

Eight fields added to `AgentState` to support thinking-mode features.
`create_initial_state` gains one new parameter (`session_blocked_programs`).
Server-side test `tst_thinking_defense.py` expected counts updated:
H02 `37` → `45` (field count), H03 `19` → `20` (param count).
Run `patch_server_tests.py <langchain_root>` to apply.

### Tests

| File | Tests added | Covers |
|------|-------------|--------|
| `tests/tst_p1_to_p5_fixes.py` | +2 | Weak anomalous detection; negligible anomalous clear |
| `patch_server_tests.py` | — | Updates `tst_thinking_defense.py` H02/H03 expected counts |

### Bugs confirmed already in tree (no new code needed)

| Bug | Root cause | Fix location |
|-----|-----------|--------------|
| Bug A | `event_formatter.py` `fmt % value` TypeError on string metrics | `agent/event_formatter.py` try/except guard |
| Bug B | `_auto_discover_files` skipped when `original_files` non-empty | `programs/ai_agent.py` supplement-mode discovery |
| Bug C (cascade) | `predict_and_build` cascade set `refine_done=True` prematurely | `agent/workflow_state.py` cascade comment + guard |
| Bug D | autosol deprioritized when `anomalous_measurability < 0.05` | `agent/workflow_engine.py` deprioritization logic |


## Version 115.02 (Report Analysis Fixes & GUI Update)

### Summary

Fixes from report analysis: program alias resolution, plan generation
guard for single-program tasks, solve-mode file discovery fix,
unsupported program handling (Sorry → warning), polder file routing,
predict_and_build label injection, and GUI alignment.

### Modified files (7)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | `phenix.resolve` → `phenix.resolve_cryo_em` alias in `_fix_program_name()`; LLM prompt negative example for resolve |
| `agent/intent_classifier.py` | Added `resolve_cryo_em` pattern to `phenix.resolve_cryo_em` aliases |
| `agent/sanity_checker.py` | Added `_check_autobuild_phib()` — warns after 2+ autobuild PHIB failures |
| `agent/thinking_agent.py` | Enhancement B: `error_classification` and `failure_count` forwarded to thinking context |
| `programs/ai_agent.py` | Fix 1: `_load_result_file_exclusions()` + filter in `_auto_discover_files()`; plan skip for task intent; solve-mode README preprocessing bypass; expanded `_README_IGNORE_PROGRAMS` with utility programs (keeps `raise Sorry` for unsupported core programs like ensemble_refinement) |
| `knowledge/programs.yaml` | Polder `data_mtz` input_priorities (exclude map_coeffs_mtz); `predict_and_build` obs_labels strategy_flag |
| `wxGUI2/Programs/AIAgent.py` | Expanded ignore set for utility programs (keeps `raise Sorry` + blocking banner for unsupported core programs) |


## Version 115.01 (Intent Classification & Infrastructure Fixes)

### Summary

Intent classification system for dual-run evaluation (tutorial vs solve
mode), plus 6 infrastructure fixes addressing file routing bugs discovered
by the first dual-run evaluation.  Expected additional cycle waste
reduction: ~15-20% from preventing technical crashes.

### New files (3)

| File | Lines | Purpose |
|------|-------|---------|
| `agent/intent_classifier.py` | 433 | Classifies user advice into solve/solve_constrained/task/tutorial |
| `knowledge/tutorial_expectations.yaml` | 218 | Expected behavior for all 21 tutorials (answer key) |
| `knowledge/solve_run_info.yaml` | 94 | Minimal hints for solve-mode runs |

### Modified files (5)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | Intent classification integration; intent-driven stop behavior; `_set_by_pattern` flag for distinguishing pattern-inferred vs user-explicit stops; intent simplified to string before return |
| `agent/advice_preprocessor.py` | Early return bypass for "Solve the structure by standard procedures" READMEs (skips LLM preprocessing) |
| `agent/workflow_state.py` | Fix I1: Multi-array MTZ detection (`_detect_mtz_arrays`, `_read_mtz_array_labels`, `_classify_mtz_arrays`); Fix I3: `phased_data_mtz` post-categorization scan for autosol outputs |
| `agent/phil_validator.py` | Fix I2: `_BLOCKED_PARAMS` dict; blocked-param check before allowed-flags check (strips `input_map_file` from autobuild) |
| `agent/graph_nodes.py` | Fix I1: MTZ label injection with SAD/MR ranking rule; Fix I4: `map_sharpening` resolution injection; Fix I5: `_estimate_anomalous_sites`/`_read_sequence` + sites injection + `atom_type` forwarding from directives + `rebuild_in_place=False` when R-free > 0.45 |

### Knowledge base changes (1)

| File | Changes |
|------|---------|
| `knowledge/programs.yaml` | Fix I1: `obs_labels` for xtriage and autosol; Fix I3: autosol output type `data_mtz` → `phased_data_mtz` + denmod pattern; Fix I8: `ligand_cif` optional input + input_priorities + command template for refine |

### Evaluation infrastructure (2)

| File | Changes |
|------|---------|
| `tests/analyze_tutorial_runs.py` | Dual-run support: `run_type` dimension (solve/readme), per-type output files, updated directory parsing |
| `solve_readmes/` | 21 solve-mode README.txt files |

### Tests

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_intent_classifier.py` | 32 | Intent classification (4 categories) |
| Total (all test files) | 106 | Fixes 2-6, Enhancement A, Intent, Infrastructure |


## Version 115.00 (Failure Handling & Stop Condition Fixes)

### Summary

Seven fixes and one enhancement addressing the top causes of
wasted agent cycles, based on analysis of 21 tutorials × 5 modes.
The core changes eliminate fabricated stop conditions (Fix 2),
add PHIL parameter validation (Fix 4), implement a tiered failure
response with error classification (Fix 3), make duplicate
detection parameter-aware (Fix 5), rescue uncategorized map files
(Fix 6), and inject structured error context for LLM
self-correction (Enhancement A).  Expected cycle waste reduction:
41% → ~23% (20 cycles saved across analyzed sessions).

### New files (2)

| File | Lines | Purpose |
|------|-------|---------|
| `agent/phil_validator.py` | 142 | PHIL strategy validation against programs.yaml strategy_flags |
| `agent/error_classifier.py` | 410 | Error classification (TERMINAL/PHIL_ERROR/LABEL_ERROR/RETRYABLE) + tiered pivot logic |

### Modified files (6)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | Fix 2: Strip fabricated "Stop Condition:" and "Goal:" headers from preprocessor output; remove "Stop Condition:" parser block; strengthen LLM extraction prompt with negative examples; guard tutorial_patterns against preprocessor text |
| `agent/graph_nodes.py` | Fix 3: Error classification in PERCEIVE, pivot logic in PLAN (before rules_only); Fix 4: PHIL validation wrapper + call in both build paths; Fix 5: Parameter-aware duplicate check in validate(); Enhancement A: Error context injection into planning prompt |
| `agent/workflow_state.py` | Fix 6: Map extension safety net in _categorize_files(); fixed _load_category_rules() ImportError handling |
| `knowledge/thinking_prompts.py` | Enhancement A: Error classification context in expert prompt ("=== PREVIOUS FAILURE ===" section) |
| `knowledge/programs.yaml` | Fix 4: Added reference_model_enabled, reference_model_use_starting to phenix.refine; added mask_atoms to phenix.resolve_cryo_em; updated comment |
| `docs/reference/ARCHITECTURE.md` | Added phil_validator.py and error_classifier.py to Key Files table |

### Tests (4 files, 68 tests)

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_stop_condition_fix.py` | 20 | Fix 2: stripping, context-awareness, real user commands preserved |
| `tests/tst_phil_validation.py` | 15 | Fix 4: valid/invalid params, rewrite regression, PLAN cases |
| `tests/tst_error_classifier.py` | 24 | Fix 3: terminal/PHIL/label/retryable errors, pivot logic, failure counting |
| `tests/tst_fix_integration.py` | 9 | Cross-fix interactions: Fix 4→3, Fix 3+5, full AF_exoV chain |

Run with: `python3 tests/tst_stop_condition_fix.py && python3 tests/tst_phil_validation.py && python3 tests/tst_error_classifier.py && python3 tests/tst_fix_integration.py`

### Items completed in v115.01/v115.02

| Item | Status |
|------|--------|
| Fix 1: Clean tutorial set | ✓ Done (v115.02) — `_load_result_file_exclusions` in ai_agent.py |
| Fix 7: CASP7 MTZ recognition | Partially mitigated — auto-discover finds files |
| Enhancement B: Model upgrade | ✓ Done (v115.02) — error context forwarded in thinking_agent.py |
| advice_preprocessor.py | ✓ Done (v115.01) — solve-mode bypass |
