# HANDOFF — PHENIX AI Agent Session

**Date:** 2026-03-16
**Upload this file at the start of the next session along with the code archive.**

---

## What This Project Is

The PHENIX AI Agent is a reactive engine that runs crystallography/cryo-EM
programs (phenix.refine, phenix.autosol, etc.) iteratively. Each cycle:
client sends session_info → server runs a 7-node LangGraph
(PERCEIVE→THINK→PLAN→BUILD→VALIDATE→FALLBACK→OUTPUT) → returns one command.
The client executes it, collects the log, and calls again.

Code lives in `agent/`, `knowledge/`, `programs/`, `phenix_ai/`, `tests/`.
The archive has ~50 Python files, ~100K lines total.

---

## Current State

**Run 17** (30/42 tutorials, OpenAI, our fixes deployed):
- 10 GOOD / 9 OK / 2 POOR / 9 NOMET
- 24% failure rate (down from 34% in Run 15 baseline)
- 4 fixes verified working (Bugs 3, 4-partial, 6, 7)

**Fixes implemented this session (7 bugs, 10 code files):**

| Bug | File(s) | Status |
|-----|---------|--------|
| 1 | `diagnosable_errors.yaml` | Not yet in Run 17 |
| 3 | `programs.yaml`, `phil_validator.py` | ✅ Verified: 57%→0% fail |
| 4 | `metrics_analyzer.py`, `structure_model.py` | ⚠ Crash persists (see below) |
| 5 | `programs.yaml`, `phil_validator.py`, `program_registry.py`, `graph_nodes.py`, `workflow_state.py` | Not yet in Run 17 |
| 6 | `workflow_state.py` | ✅ Verified: 26%→2% fail |
| 7 | `workflow_state.py` | ✅ Verified: polder runs |
| cleanup | `graph_nodes.py` | Done |

---

## CRITICAL: Two Unfixed Bugs for Next Session

### Bug 4 Continued — `metric_evaluator.py` needs `_safe_float`

**Symptom:** 7rpq crashes at cycle 3 with `unsupported operand type(s) for -: 'str' and 'float'` despite `metrics_analyzer.py` fix being deployed.

**Root cause:** `USE_YAML_METRICS=True` (line 103 of `graph_nodes.py`) routes
`analyze_metrics_trend()` through `metric_evaluator.py::analyze_refinement_trend()`
instead of the hardcoded path in `metrics_analyzer.py`. The YAML evaluator
(`metric_evaluator.py`) has **zero** `_safe_float` coercion — it reads raw
`m.get("r_free")` at line 352 and passes it to `calculate_improvement_rate()`
(line 226/252) and `is_plateau()` (line 284) which do subtraction.

The `derive_metrics_from_history()` in `metrics_analyzer.py` correctly coerces
values, but `metric_evaluator.py` re-reads them from the metrics_history dicts
without coercion, and also reads from the `metrics` key (not `analysis`) via
`_analyze_history()` in `workflow_state.py` lines 1839-1842 which does NO coercion.

**Fix needed in `metric_evaluator.py`:**
- Add `_safe_float()` function (same as in metrics_analyzer.py)
- Coerce at line 352: `r_free = _safe_float(m.get("r_free"))`
- Coerce at line 392: `previous = _safe_float(r_free_values[-2])` (belt-and-suspenders)
- Coerce in `calculate_improvement_rate()` (lines 226, 229, 252, 255)
- Coerce in `is_plateau()` (lines 284, 286)
- Also coerce CC values in `_analyze_cryoem_trend()` (same pattern)

**Also add belt-and-suspenders in:**
- `graph_nodes.py` lines 788/791 (event emission)
- `event_formatter.py` line 925 (r_free - r_free_prev)
- `workflow_state.py` lines 1841-1850 (_analyze_history metric reads)
- `kb_tags.py` lines 239, 256, 268 (r_free arithmetic in tag generation)

### Orphan-Map Promotion — `workflow_state.py`

**Symptom:** apoferritin_denmod_dock rules_only mode loops dock_in_map 3x,
never offers real_space_refine. LLM modes work (CC=0.817).

**Root cause:** `emd-20026_auto_sharpen_A.ccp4` is categorized into `map`
(parent) but excluded from `full_map` (subcategory) by the YAML pattern
`*_a.*` which is intended for half-map suffixes but false-positives on the
`_A` in `auto_sharpen_A`. Since `has_full_map=False`, `real_space_refine`
(which has `requires_full_map: true`) is never offered.

LLM modes work because they run mtriage first → mtriage output creates a
resolution entry → workflow advances through a different path.

**Fix:** Add post-processing in `_categorize_files()`: any file in the `map`
parent category that is NOT in any map subcategory (`full_map`, `half_map`,
`optimized_full_map`) should be promoted to `full_map`. This mirrors the
existing orphan-PDB promotion (pdb → model) at lines 461-498.

---

## Test Suite

291 tests across 7 suites, all passing:
- tst_phase3_bug5: 85 (covers all session fixes)
- tst_structure_model: 77
- tst_plan_schema: 53
- tst_backward_compat: 26
- tst_command_builder: 22
- tst_phil_validation: 15
- tst_event_system: 13

---

## Archives

Both at `/mnt/user-data/outputs/`:
- `phenix_ai_agent_updates.tar.gz` (39 files, full agent code)
- `phenix_phase1_fixes.tar.gz` (25 files, changed files only)

These do NOT yet include the metric_evaluator.py fix or orphan-map fix.

---

## Run Results on Disk

- `r15b_openai/` — Run 15 baseline (report JSONs)
- `r16_openai/` — Run 16 pre-deploy (report JSONs)
- `r17_openai/` — Run 17 partial (22 tut, report JSONs)
- `r17_full/` — Run 17 fuller (30 tut, report JSONs)
- `r17_logs/` — Run 17 workflow logs (7rpq + apoferritin_denmod_dock)
- `r16_degraded/` — Run 16 degraded tutorial logs

---

## Still Missing from Run 17

13 tutorials not yet run: lysozyme-refine (Bug 1), lowres_restraints (Bug 5),
hipip-refine (regression guard), groel_dock_refine, gene-5-mad, 1029B-sad,
ion_channel_denmod, lysozyme-MRSAD, model-building-scripting, nsf-d2-ligand,
nsf-d2-refine, groel_map_symmetry, if5a-textal.

---

## Deliverables Produced This Session

- `TestingStatus_Run17.pptx` — 5-slide presentation
- `PROJECT_STATUS.md` — comprehensive project status
- `DEFERRED_CHANGES_PLAN.md` — prioritized issue list
- `RUN16_VERIFICATION.md` — verification run commands
- `BUG5_FIX_PLAN.md` — reference model restraints plan
- `PHASE3_FIX_PLAN.md` — Phase 3 bug analysis
- `DebuggingCycle_Run14_to_Run15.pptx` — before/after slide
