# CHANGELOG — v115

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
