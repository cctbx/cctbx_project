# CHANGELOG — v115

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

### Documentation (2)

| File | Purpose |
|------|---------|
| `docs/project/INFRASTRUCTURE_ISSUES.md` | 8 issues identified from dual-run evaluation |
| `docs/project/FIX_PLAN_INFRASTRUCTURE.md` | Detailed fix plan for I1-I5, I8 |
| `docs/project/DUAL_RUN_EVALUATION_PLAN.md` | Framework design for tutorial vs solve evaluation |

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

### New files (4)

| File | Lines | Purpose |
|------|-------|---------|
| `agent/phil_validator.py` | 142 | PHIL strategy validation against programs.yaml strategy_flags |
| `agent/error_classifier.py` | 410 | Error classification (TERMINAL/PHIL_ERROR/LABEL_ERROR/RETRYABLE) + tiered pivot logic |
| `docs/guides/FIX2_ADVICE_PREPROCESSOR_SPEC.md` | 43 | Spec for defense-in-depth changes to advice_preprocessor.py |
| `docs/guides/FIX7_CASP7_MTZ_SPEC.md` | 48 | Spec for CASP7 MTZ recognition (client-side work needed) |

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

### Remaining items

| Item | Status | Notes |
|------|--------|-------|
| Fix 1: Clean tutorial set | Not code | Remove result files from tutorial directories |
| Fix 7: CASP7 MTZ recognition | Spec only | Client-side file discovery issue |
| Enhancement B: Model upgrade for recovery | Needs thinking_agent.py | Upgrade to thinking prompt for recovery decisions |
| advice_preprocessor.py | Spec only | Remove or constrain "Stop Condition" field at source |
