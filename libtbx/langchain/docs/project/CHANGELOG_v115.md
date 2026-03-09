# CHANGELOG — v115 (Failure Handling & Stop Condition Fixes)

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
