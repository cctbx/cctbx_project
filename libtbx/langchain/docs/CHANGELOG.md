# CHANGELOG — v115

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
