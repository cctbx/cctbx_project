# PHENIX AI Agent - Changelog

## Version 102 (January 2025)

### Bug Fix: Retry on Model Overload (503 UNAVAILABLE)

**Problem**: When Gemini returns a 503 UNAVAILABLE error ("The model is overloaded"), the summarization would fail immediately without retrying.

**Solution**: 
1. Added "503" and "unavailable" to rate limit indicators in `rate_limit_handler.py`
2. Added retry logic to `summarize_log()` with exponential backoff (2s, 4s, 8s) for rate limit and overload errors

### Behavior

When model is overloaded:
```
Summarizing log file (using cheap model)...
Model overloaded/rate limited, waiting 2s before retry...
Summarization retry 2/3...
Model overloaded/rate limited, waiting 4s before retry...
Summarization retry 3/3...
[success or final error]
```

### Files Changed

- `agent/rate_limit_handler.py` - Added "503" and "unavailable" to rate limit indicators
- `phenix_ai/run_ai_analysis.py` - Added retry logic to `summarize_log()` with exponential backoff

---

## Version 101 (January 2025)

### Bug Fix: HTML Summary Table Formatting After Failed Steps

**Problem**: When a step failed, the error message in the "Key Metric" column could contain newlines or pipe characters, breaking the markdown table formatting for all subsequent rows.

Example broken output:
```
| 3 | ligandfit | ✗ | FAILED: Sorry: Sorry LigandFit failed
Please... | | 4 | pdbtools | ✓ | Completed |
```

**Solution**: Sanitize the `key_metric` field before adding to table:
- Replace newlines with spaces
- Replace pipe characters (`|`) with dashes
- Collapse multiple spaces
- Truncate to 60 characters max

### Files Changed

- `agent/session.py` - Sanitize key_metric in `_format_summary_as_markdown()`

---

## Version 100 (January 2025)

### Bug Fix: Ligand File Detection Pattern

**Problem**: `7qz0_ligand.pdb` wasn't detected as a ligand file because patterns only matched `ligand_*.pdb` (ligand at start), not `*_ligand.pdb` (ligand at end).

**Solution**: Added patterns `*_ligand.pdb` and `*_lig.pdb` to `ligand_pdb` category.

### Bug Fix: predict_and_build Map Coefficients for Ligandfit

**Problem**: After `predict_and_build`, ligandfit was using `*_refinement.mtz` which doesn't contain map coefficients. The correct file is `*_map_coeffs.mtz` with labels `FP PHIFP`.

**Solution**:
1. Added new `predict_and_build_mtz` category for `*map_coeffs*.mtz` files
2. Updated ligandfit input priorities to include `predict_and_build_mtz`
3. Added `denmod_mtz` and `predict_and_build_mtz` to `specific_subcategories` in command builder so category-based selection takes priority over `best_files`
4. Added invariant to switch labels to `"FP PHIFP"` when using map_coeffs files

### Files Changed

- `knowledge/file_categories.yaml`:
  - Added `*_ligand.pdb` and `*_lig.pdb` patterns to `ligand_pdb`
  - Added `predict_and_build_mtz` category for map coefficient files
- `knowledge/programs.yaml`:
  - Updated ligandfit mtz priorities: `[denmod_mtz, predict_and_build_mtz, refined_mtz]`
  - Added invariant for predict_and_build labels (`FP PHIFP`)
- `agent/command_builder.py`:
  - Added `denmod_mtz` and `predict_and_build_mtz` to `specific_subcategories`
- `tests/test_file_categorization.py`:
  - Added `test_ligand_file_patterns`
  - Added `test_predict_and_build_mtz_detection`

---

## Version 99 (January 2025)

### Feature: maximum_automation=False Now Works for X-ray

**Previously**: `maximum_automation=False` (stepwise mode) only affected cryo-EM workflows, forcing `stop_after_predict=True` for `predict_and_build`.

**Now**: Stepwise mode also applies to X-ray workflows. When `maximum_automation=False`:
- `predict_and_build` will use `stop_after_predict=True` in states: `xray_initial`, `xray_placed`
- This gives users more control over the workflow with intermediate checkpoints
- User can then run `process_predicted_model` → `phaser` → `refine` separately

### Usage

```bash
# Stepwise mode - more control with intermediate checkpoints
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"
```

### Workflow Comparison

**Automated (maximum_automation=True, default)**:
```
xray_initial → xtriage → predict_and_build(full) → xray_refined
```

**Stepwise (maximum_automation=False)**:
```
xray_initial → xtriage → predict_and_build(stop_after_predict)
                              ↓
              process_predicted_model → phaser → refine → xray_refined
```

### Files Changed

- `agent/graph_nodes.py` - Extended stepwise mode handling to X-ray states
- `agent/docs_tools.py` - Updated workflow documentation diagrams
- `agent/workflow_state.py` - Updated stepwise hint message
- `docs/README.md` - Added automation modes section and quick start example
- `tests/test_integration.py` - Added `test_xray_stepwise_forces_stop_after_predict`
- `tests/test_workflow_state.py` - Added X-ray stepwise tests

---

## Version 98 (January 2025)

### Bug Fix: predict_and_build Counts as Refinement for ligandfit

**Problem**: After running `phenix.predict_and_build`, `phenix.ligandfit` was unavailable because `refine_count=0`:
```
PERCEIVE: phenix.ligandfit unavailable: refine_count=0 does not satisfy condition '> 0'
```

But predict_and_build includes internal refinement cycles and produces the same outputs (map coefficients) that ligandfit needs.

**Solution**: When a full `predict_and_build` run completes (not stopped early), increment `refine_count` so downstream programs like `ligandfit` know there's a refined model with map coefficients.

### Files Changed

- `agent/workflow_state.py` - Increment `refine_count` for successful full predict_and_build runs

---

## Version 97 (January 2025)

### Bug Fix: PredictAndBuild Output Categorized as Model (Not Search Model)

**Problem**: Files like `PredictAndBuild_0_overall_best.pdb` were incorrectly categorized as `search_model` instead of `model`. This happened because:
1. The file contains "predict" in the name
2. Multiple categories had `excludes: ["*predict*"]`
3. The file fell through to the `predicted` category (parent: `search_model`)

This caused the sanity check to fail with:
```
PERCEIVE: RED FLAG [search_model_not_positioned]: Cannot refine: search model found but not yet positioned
```

**Solution**: Added new `predict_and_build_output` category that specifically matches `PredictAndBuild_*_overall_best*.pdb` files and categorizes them as `model` (positioned, ready for refinement).

Also added exclusions to the `predicted` category to ensure these files don't get double-categorized.

### Files Changed

- `knowledge/file_categories.yaml` - Added `predict_and_build_output` category, added exclusions to `predicted`
- `tests/test_file_categorization.py` - Added test for PredictAndBuild output categorization

---

## Version 96 (January 2025)

### Bug Fix: LLM Slot Alias Mapping for MTZ Files

**Problem**: When the LLM requested an MTZ file using the slot name `data` (e.g., `data=PredictAndBuild_0_overall_best_refinement.mtz`), but the program defined the input slot as `mtz`, the LLM's file choice was ignored and the wrong file was auto-selected.

**Example Debug Output (Before)**:
```
BUILD: LLM requested files: {model=..., data=PredictAndBuild_0_overall_best_refinement.mtz}
BUILD: Skipping best_files for mtz (program needs specific subcategory)
BUILD: Auto-filled mtz=PredictAndBuild_0_refinement_cycle_2.extended_r_free.mtz  # WRONG!
```

**Solution**: Added `SLOT_ALIASES` mapping that translates common LLM slot names to canonical program input names:
- `data` → `mtz`
- `pdb` → `model`
- `seq_file` → `sequence`
- etc.

**Example Debug Output (After)**:
```
BUILD: LLM requested files: {model=..., data=PredictAndBuild_0_overall_best_refinement.mtz}
BUILD: Mapped LLM slot 'data' to 'mtz'
BUILD: Using LLM-selected file for mtz
```

### Files Changed

- `agent/command_builder.py` - Added `SLOT_ALIASES` dict and updated LLM file processing to use aliases
- `tests/test_command_builder.py` - Added tests for slot alias mapping

---

## Version 94 (January 2025)

### New Feature: Explain Why Programs Are Unavailable

When a program like `phenix.ligandfit` is not available, the debug output now explains WHY:

```
PERCEIVE: Valid programs: phenix.refine, phenix.polder, STOP
PERCEIVE: phenix.ligandfit unavailable: missing required file: ligand_file
```

This helps diagnose issues when expected programs aren't offered.

### Possible Explanations

- `missing required file: ligand_file` - No ligand file (.pdb/.cif) detected
- `missing required file: sequence` - No sequence file (.fa/.seq) detected
- `already completed: ligandfit` - Program already ran (not_done condition failed)
- `r_free=0.40 does not satisfy condition '< 0.35'` - Metric threshold not met
- `refine_count=0 does not satisfy condition '> 0'` - Needs refinement first
- `run_once program already executed` - Program like xtriage already ran

### Files Changed

- `agent/workflow_engine.py` - Added `explain_unavailable_program()` method, added `unavailable_explanations` to workflow state
- `agent/graph_nodes.py` - Added debug logging for unavailable programs

---

## Version 93 (January 2025)

### New Feature: Density Modification Workflow for X-ray

- **Added `phenix.autobuild_denmod` to X-ray workflow**: Before ligand fitting, the agent can now run density modification using `phenix.autobuild maps_only=True`. This creates improved map coefficients (`overall_best_denmod_map_coeffs.mtz` with FWT/PHFWT labels) for better ligand fitting.

### Workflow Changes

The X-ray refine phase now includes:
1. `phenix.refine` (preferred) - standard refinement
2. `phenix.autobuild` - when R-free stuck above threshold
3. **`phenix.autobuild_denmod`** (NEW) - density modification before ligandfit
4. `phenix.ligandfit` - fit ligand when model is good enough
5. `phenix.polder` - omit map calculation

### Technical Changes

- **Prompt clarification**: Added warning that `predict_and_build` is NOT for density modification
- **New file category**: `denmod_mtz` for density-modified MTZ files
- **Ligandfit label switching**: Automatically uses `FWT PHFWT` labels when denmod MTZ is selected
- **Done flag tracking**: Added `autobuild_denmod_done` flag

### Documentation Updates

- **docs/OVERVIEW.md**: Updated workflow example with autobuild_denmod, updated done flags table
- **docs/guides/ADDING_PROGRAMS.md**: Added autobuild_denmod_done to flags table, added note about refine_count/rsr_count
- **docs/guides/TESTING.md**: Updated test table with counts and added "Key Tests for Recent Fixes" section
- **tests/run_all_tests.py**: Updated docstring with current test list and key tests

### Files Changed

- `knowledge/workflows.yaml` - Added autobuild_denmod to refine phase
- `knowledge/programs.yaml` - Added denmod_labels invariant to ligandfit
- `knowledge/file_categories.yaml` - Added denmod_mtz category
- `knowledge/prompts_hybrid.py` - Added autobuild_denmod description, warned about predict_and_build
- `agent/workflow_state.py` - Added autobuild_denmod_done flag and detection
- `agent/command_builder.py` - Added file_matches invariant handling
- `docs/OVERVIEW.md` - Updated workflow examples and done flags
- `docs/guides/ADDING_PROGRAMS.md` - Updated done flags table
- `docs/guides/TESTING.md` - Updated test descriptions
- `tests/run_all_tests.py` - Updated docstring

---

## Version 72 (January 2025)

### Bug Fixes

- **Fixed cycle count to exclude STOP cycles (v92)**: The session summary was reporting "Cycles: 5 (4 successful)" when 4 programs ran and cycle 5 was just STOP. Now STOP cycles are excluded from the count, so it correctly reports "Cycles: 4 (4 successful)".

### Tests Added

- `test_stop_cycle_excluded_from_count` - Verifies STOP cycles are excluded from total_cycles and successful_cycles
- `test_cryoem_done_flags` - Verifies done flags are set for cryo-EM programs (resolve_cryo_em_done, map_sharpening_done, map_to_model_done, dock_done)

### Documentation Updated

- `docs/OVERVIEW.md` - Added "Program Execution Controls" section documenting `not_done` conditions and `run_once` flags with complete table of done flags
- `docs/guides/ADDING_PROGRAMS.md` - Added "Available done flags" table showing all manually-tracked done flags
- `docs/guides/TESTING.md` - Updated test count table with new test files

### Files Changed

- `agent/session.py` - Exclude STOP/None/unknown from total_cycles and successful_cycles
- `agent/session_tools.py` - Exclude STOP cycles from print_session_summary()
- `tests/test_session_summary.py` - Added test_stop_cycle_excluded_from_count
- `tests/test_workflow_state.py` - Added test_cryoem_done_flags

---

## Version 71 (January 2025)

### Bug Fixes

- **Fixed predict_and_build running without resolution for X-ray (v91)**: When user provides `program_settings` for predict_and_build (e.g., `rebuilding_strategy=Quick`), the program was being added to valid_programs even before xtriage ran. This caused predict_and_build to run without resolution, forcing `stop_after_predict=True` (prediction only). Now `_check_program_prerequisites` requires xtriage_done (X-ray) or mtriage_done (cryo-EM) before adding predict_and_build from program_settings.

- **Fixed resolution requirement for predict_and_build full workflow**: The command builder now correctly requires resolution for BOTH X-ray and cryo-EM when `stop_after_predict=False`. If resolution is not available, it forces `stop_after_predict=True` with a message suggesting to run xtriage/mtriage first.

### Root Cause

The workflow engine's `_apply_directives` was adding `predict_and_build` to valid_programs whenever the user had `program_settings` for it:

```python
# Before: Always allowed predict_and_build
if program == "phenix.predict_and_build":
    return True  # Always allow - worst case it does prediction-only
```

This bypassed the normal workflow phase ordering (xtriage → obtain_model).

### The Fix

```python
# After: Require xtriage/mtriage to be done first
if program == "phenix.predict_and_build":
    if phase_name in ("obtain_model", "molecular_replacement", "dock_model"):
        return True  # Let the phase conditions handle it
    if context.get("xtriage_done") or context.get("mtriage_done"):
        return True
    return False  # Don't add to early phases
```

### Correct Workflow Now

1. xtriage runs first → extracts resolution
2. predict_and_build runs with resolution → full workflow (prediction + MR + building)

### Files Changed

- `agent/workflow_engine.py` - Fixed `_check_program_prerequisites` to require xtriage/mtriage
- `agent/command_builder.py` - Fixed `_apply_invariants` to require resolution for full workflow

---

## Version 70 (January 2025)

### Bug Fixes

- **Fixed programs running repeatedly without stopping (v90)**: Multiple programs were missing `not_done` conditions in their workflow definitions, allowing the LLM to choose them repeatedly even after successful completion. Added protection to all one-time-run programs.

### Programs Now Protected from Re-runs

**X-ray workflow:**
- `phenix.predict_and_build` - `not_done: predict_full`
- `phenix.phaser` - `not_done: phaser` (in both obtain_model and molecular_replacement phases)

**Cryo-EM workflow:**
- `phenix.predict_and_build` - `not_done: predict`
- `phenix.dock_in_map` - `not_done: dock`
- `phenix.map_to_model` - `not_done: map_to_model`
- `phenix.resolve_cryo_em` - `not_done: resolve_cryo_em` (in all 4 phases where it appears)
- `phenix.map_sharpening` - `not_done: map_sharpening` (in all 4 phases where it appears)

**Previously protected (unchanged):**
- `phenix.autobuild`, `phenix.autosol`, `phenix.ligandfit`, `phenix.map_symmetry`, `phenix.process_predicted_model`

### Programs Intentionally Without Protection (run multiple times)

- `phenix.refine` / `phenix.real_space_refine` - Iterative refinement
- `phenix.molprobity` / `phenix.model_vs_data` / `phenix.validation_cryoem` - Validation after each cycle
- `phenix.polder` - May run for different sites
- `phenix.pdbtools` - May add multiple ligands
- `phenix.xtriage` / `phenix.mtriage` - Already have `run_once: true` in programs.yaml

### Files Changed

- `knowledge/workflows.yaml` - Added `not_done` conditions to 11 program entries
- `agent/workflow_state.py` - Added done flags for `map_to_model`, `resolve_cryo_em`, `map_sharpening`

---

## Version 69 (January 2025)

### Bug Fixes

- **Fixed predict_and_build forcing `stop_after_predict=True` for X-ray (v89)**: The command builder was forcing `stop_after_predict=True` whenever resolution wasn't in the context, even for X-ray where predict_and_build can read resolution directly from the MTZ file. This caused predict_and_build to run in prediction-only mode repeatedly (3 retry cycles due to duplicate detection). Now `stop_after_predict=True` is only forced for cryo-EM without resolution.

### Root Cause

The command builder had this logic:
```python
if program == "phenix.predict_and_build":
    if not context.resolution and "stop_after_predict" not in strategy:
        strategy["stop_after_predict"] = True
```

This applied to BOTH X-ray and cryo-EM, but:
- For X-ray: predict_and_build can determine resolution from the MTZ automatically
- For cryo-EM: mtriage should run first to get resolution

The fix restricts this to cryo-EM only:
```python
if context.experiment_type == "cryoem":
    if not context.resolution and "stop_after_predict" not in strategy:
        strategy["stop_after_predict"] = True
```

### Why 3 LLM Calls Were Happening

1. LLM chose predict_and_build
2. Command builder forced `stop_after_predict=True` (same as previous run)
3. Validate detected duplicate command → `validation_error`
4. Graph looped back to Plan (retry 1)
5. Same result → retry 2
6. Same result → retry 3 → fallback

### Files Changed

- `agent/command_builder.py` - Only force `stop_after_predict=True` for cryo-EM

---

## Version 68 (January 2025)

### Bug Fixes

- **Fixed ligandfit using input MTZ instead of refined MTZ (v88)**: When ligandfit requires `refined_mtz` (an MTZ with map coefficients), and no refined MTZ exists (because refinement failed), the code was falling back to extension-based matching and selecting the input MTZ. Now when a program requires a specific subcategory (like `refined_mtz`), extension-based fallback is disabled.

- **Fixed refine_count incrementing for failed refinements (v88)**: Previously, `refine_count` was incremented regardless of whether refinement succeeded or failed. This caused workflow conditions like `refine_count > 0` (used by ligandfit) to pass even when no successful refinement had occurred. Now only successful refinements increment the count.

- **Fixed rsr_count incrementing for failed RSR (v88)**: Same fix applied to `rsr_count` for `phenix.real_space_refine`.

### Root Cause

The session showed:
```
Categorized files: model=1, search_model=15, mtz=1, sequence=10, ...
```
No `refined_mtz` category because refinement failed. But `refine_count` was 2 (from failed attempts), so ligandfit was allowed. Then command_builder fell back to extension matching and selected the input MTZ.

### Files Changed

- `agent/command_builder.py` - Added check to prevent extension-based fallback when specific subcategory required
- `agent/workflow_state.py` - Added success checking for refine_count and rsr_count in `_analyze_history()`
- `tests/test_workflow_state.py` - Added `test_failed_refine_not_counted`

---

## Version 67 (January 2025)

### Bug Fixes

- **Fixed stop_after_predict=True being suggested for X-ray (v87)**: The prompt was incorrectly telling the LLM to use `stop_after_predict=True` for ANY stepwise workflow, but this should only apply to cryo-EM stepwise workflows. For X-ray, `predict_and_build` should run the full workflow (prediction → MR → building). Added experiment_type check so the guidance only appears for cryo-EM.

- **Clarified predict_and_build documentation in prompts**: Added explicit notes explaining that by default `predict_and_build` runs the FULL workflow, and `stop_after_predict=True` should only be used for cryo-EM stepwise.

### Root Cause

When user says "stop after PredictAndBuild", the LLM was confusing:
- "Stop the agent workflow after predict_and_build completes" (correct interpretation)
- "Set stop_after_predict=True" (incorrect - this only runs prediction, skipping MR and building)

The prompt was adding `NOTE: Use predict_and_build with strategy: {"stop_after_predict": true}` for stepwise workflows without checking if it's cryo-EM or X-ray.

### Files Changed

- `knowledge/prompts_hybrid.py` - Added experiment_type check for stop_after_predict guidance; clarified predict_and_build documentation

---

## Version 66 (January 2025)

### Critical Bug Fixes

- **Fixed predicted model incorrectly becoming best model after refinement (v85)**: When `phenix.refine` runs, the directory scan picks up ALL files (including pre-existing ones like `PredictAndBuild_0_predicted_model_processed.pdb`). Previously, `record_result()` blindly applied `stage="refined"` to ALL PDB files in `output_files`, causing the predicted model to get a higher score than the actual refined model. Now `record_result()` only applies the program stage to files whose basename matches expected output patterns (e.g., only files containing "refine" get `stage="refined"`).

- **Fixed PHASER models getting inflated scores from refinement metrics**: Similar to above - `PHASER.1.pdb` in refinement's `output_files` was getting `stage="refined"` and metrics from phenix.refine, inflating its score from 70 to 132. Now handled by the same pattern matching fix.

- **Fixed STOP not available after user's workflow completes (v84)**: When user directives specify `after_program` (e.g., "stop after refinement"), STOP is now added to valid_programs after that program completes. Previously STOP was only available after validation, forcing unwanted extra cycles.

### Consistency Fixes (v86)

- **Synchronized pattern matching across all three locations**:
  - `session.py:_rebuild_best_files_from_cycles()` - Rebuild from saved session
  - `session.py:record_result()` - Real-time recording
  - `session_tools.py:rebuild_best_files()` - Manual rebuild tool

- **Added missing program-to-stage mappings**:
  - `session.py:_infer_stage_from_program()` - Added `predict_and_build`, `ligandfit`, `pdbtools`
  - `session_tools.py:infer_stage()` - Added `process_predicted_model`

- **Added missing filename patterns**:
  - `session_tools.py` - Added `processed_predicted` and `autobuild_output` patterns

### Bug Details

**The predicted model bug (v85)**:
- Root cause: `output_files` from directory scanning includes ALL files in the working directory, not just files created by the program
- `record_result()` was giving ALL PDB files the program's stage (e.g., "refined" for phenix.refine)
- This caused `PredictAndBuild_0_predicted_model_processed.pdb` to get `stage="refined"` with the refinement metrics, giving it score 133 vs PHASER.1.pdb's score 132
- Fix: Pattern matching in `record_result()` now matches `_rebuild_best_files_from_cycles()` - only files with matching basenames get the program stage

### Testing

- Added test `test_predicted_model_not_promoted_by_refine` to verify predicted models don't get wrongly promoted when they appear in refine's output_files
- Added test `test_phaser_model_not_promoted_by_refine_metrics` to verify PHASER models don't get refinement metrics
- Added test `test_stop_added_after_after_program_completes` to verify STOP is available after after_program completes

### Files Changed

- `agent/session.py` - Fixed `record_result()` to pattern-match PDB filenames before applying stage; added missing programs to `_infer_stage_from_program()`
- `agent/session_tools.py` - Added `process_predicted_model` to `infer_stage()`; added `processed_predicted` and `autobuild_output` filename patterns
- `agent/workflow_engine.py` - Added after_program completion check in `_apply_directives()`
- `tests/test_best_files_tracker.py` - Added predicted model and PHASER model promotion tests
- `tests/test_workflow_state.py` - Added STOP after after_program test

---

## Version 65 (January 2025)

### Major Bug Fixes

- **Fixed map files getting wrong stage in best_files (v79)**: All map files from resolve_cryo_em were incorrectly getting `stage=optimized_full_map` instead of being classified by filename patterns. This caused `initial_map.ccp4` to win over `denmod_map.ccp4`. Now map files are classified based on their basename: `initial*` → intermediate_map, `denmod*` → optimized_full_map, `sharp*` → sharpened, `half*` → half_map.

- **Fixed rebuild function missing MTZ/PDB stage patterns (v80)**: The `_rebuild_best_files_from_cycles()` function was missing:
  - Phased MTZ detection (`*phased*`, `*phases*`, `*solve*` → `phased_mtz`)
  - Generic MTZ fallback (→ `mtz` stage)
  - `processed_predicted` pattern for process_predicted_model outputs
  - `autobuild_output` pattern for autobuild outputs

- **Fixed autobuild_done set even on failure (v77)**: Previously `autobuild_done=True` was set just because autobuild appeared in history, regardless of success/failure. Now checks for "FAIL", "SORRY", "ERROR" in result before marking done, allowing the agent to try alternatives when autobuild fails.

- **Fixed LLM file suggestions with wrong extension (v77)**: Added extension validation when LLM suggests files for input slots. Now rejects files with wrong extension (e.g., `.ccp4` for model slot which expects `.pdb/.cif`).

- **Allow programs from directives even when workflow state is "past" that phase (v81)**: When user has `program_settings` for a program (e.g., `phenix.predict_and_build`), that program is now added to `valid_programs` even if the workflow state thinks we're past that phase. This allows users to explicitly request earlier-phase programs.

### New Features

- **Programs from program_settings added to valid_programs**: If user directives include settings for a specific program, that program is automatically added to the list of valid programs (subject to prerequisite checks). This respects user intent when they've configured a program they want to run.

- **_check_program_prerequisites() helper**: New method in WorkflowEngine that centralizes prerequisite checking for programs being added via directives. Checks:
  - Refinement programs need a model to refine
  - Ligandfit needs prior refinement (for map coefficients)
  - predict_and_build is always allowed (worst case: prediction-only)

### Testing

- Added tests for directive-based program addition in test_workflow_state.py:
  - `test_program_settings_adds_program_to_valid`
  - `test_program_settings_prioritizes_program`
  - `test_program_settings_respects_prerequisites`
  - `test_default_program_settings_ignored`

### Files Changed

- `agent/workflow_engine.py` - Added `_check_program_prerequisites()`, enhanced `_apply_directives()`
- `agent/workflow_state.py` - Fixed autobuild_done to check for success
- `agent/session.py` - Fixed map file stage assignment, complete rebuild function
- `agent/command_builder.py` - Added LLM file extension validation
- `tests/test_workflow_state.py` - Added directive program tests

---

## Version 64 (January 2025)

### Major Bug Fixes

- **Fixed best_files rebuild from cycle history (v69)**: Previously, `best_files` was persisted in session.json and could get stale when cycles were removed. Now `Session.load()` always rebuilds `best_files` from the cycle history, ensuring consistency. Removed redundant `best_files.evaluate_file()` call from `_track_output_files()` that was overwriting good evaluations with `metrics=None`.

- **Fixed smart stage assignment in rebuild (v71)**: When rebuilding best_files, stage was blindly applied to all files in a cycle's output_files. Now only applies program-specific stage (e.g., "refined") to files whose basename matches expected patterns (e.g., contains "refine"). This prevents PHASER.1.pdb from incorrectly getting `stage=refined` when it appears in a refinement cycle's output_files.

- **Fixed after_cycle directive for ligand workflows (v72)**: When user says "stop after second refinement" with a ligand workflow, LLM incorrectly extracted `after_cycle: 2`. Extended `_fix_ligand_workflow_conflict()` to also clear `after_cycle <= 4` when ligand constraints are present, since ligand workflows need ~8 cycles minimum.

- **Fixed ligandfit using wrong MTZ file (v73)**: phenix.ligandfit needs an MTZ with map coefficients (2FOFCWT, PH2FOFCWT) from refinement, but was getting the original data MTZ. Added logic to skip generic `best_files["mtz"]` when `input_priorities` specifies a specific subcategory like `refined_mtz`.

- **Fixed cryo-EM dock_in_map using initial_map instead of denmod_map (v74)**: After resolve_cryo_em runs, dock_in_map was selecting `initial_map.ccp4` (intermediate) instead of `denmod_map.ccp4` (density-modified output). Added `optimized_full_map` category with score 100, `intermediate_map` with score 5, and proper pattern matching for `denmod_map`, `density_modified`, etc.

- **Fixed pdbtools output naming**: Added `fixes.output_name` to pdbtools configuration to generate output filenames like `{protein_base}_with_ligand.pdb`, ensuring the combined model is properly categorized for downstream programs.

### New Categories & Scoring

**Map Categories (v74):**
| Stage | Score | Description |
|-------|-------|-------------|
| optimized_full_map | 100 | denmod_map, density_modified, sharpened |
| sharpened | 90 | Sharpened maps |
| full_map | 50 | Regular full reconstructions |
| half_map | 10 | Half-maps for FSC |
| intermediate_map | 5 | initial_map (resolve_cryo_em intermediate) |

### Testing

- Added tests for optimized_full_map scoring and classification
- Added tests for intermediate_map low priority
- Added tests for ligand workflow after_cycle clearing
- Added tests for docked model bubbling to model category

### Files Changed

- `agent/session.py` - Added `_rebuild_best_files_from_cycles()`, always rebuild on load
- `agent/session_tools.py` - Updated `rebuild_best_files()` with smart stage assignment
- `agent/command_builder.py` - Skip best_files for specific subcategories
- `agent/directive_extractor.py` - Extended ligand workflow fix for after_cycle
- `agent/best_files_tracker.py` - Added `intermediate_map` stage and scoring
- `knowledge/file_categories.yaml` - Added `optimized_full_map`, `intermediate_map` categories
- `knowledge/programs.yaml` - Updated resolve_cryo_em outputs, pdbtools output naming
- `programs/ai_agent.py` - Added `placed_model`/`docked` to valuable_output_patterns
- `tests/test_best_files_tracker.py` - New map scoring tests
- `tests/test_directive_extractor.py` - New ligand workflow tests
- `tests/test_file_categorization.py` - New cryo-EM map categorization tests

---

## Version 63 (January 2025)

### Major Bug Fixes

- **Fixed predict_and_build output tracking**: Output files in `CarryOn/` directories were being incorrectly excluded as "intermediate files". Added `valuable_output_patterns` list that overrides intermediate exclusions for important outputs like `*_predicted_model*.pdb`, `*overall_best*.pdb`, `*_processed*.pdb`.

- **Fixed workflow state detection for predicted models**: The `_has_placed_model()` function was incorrectly returning `True` when user directives mentioned `phenix.refine` AND any PDB file existed - even if that PDB was a search_model (unpositioned). Now properly checks that PDB files are not in `search_model`, `predicted`, or `processed_predicted` categories before considering the model "placed".

- **Fixed pdbtools file selection**: Added `input_priorities` to `phenix.pdbtools` in programs.yaml to properly select refined model + ligand instead of predicted model. Added support for `priority_patterns` and `prefer_patterns` in command_builder.py.

- **Fixed ligand workflow directive conflict**: When user wants "refine, fit ligand, refine again", the directive extractor was incorrectly setting `after_program: phenix.refine` which stopped at the first refinement. Added `_fix_ligand_workflow_conflict()` post-processing that clears `after_program` when ligand-related constraints are present.

### New Features

- **Automatic Safety Documentation Generator** (`generate_safety_docs.py`): Script that scans the codebase and generates a comprehensive table of all safety checks, validations, and post-processing corrections. Run with `python generate_safety_docs.py > docs/SAFETY_CHECKS.md`.

- **Simplified Verbosity Levels**: Reduced from 4 levels to 3 levels (quiet/normal/verbose). The `debug` level has been removed; `verbose` now includes all detailed output including file selection, LLM traces, and internal state. For backwards compatibility, `debug` is accepted as input but treated as `verbose`.

- **File existence retry mechanism**: Added retry logic (3 attempts, 0.1s delay) in `resolve_file_path()` to handle race conditions where log mentions a file before it's fully written to disk.

- **Session file fsync**: Added explicit `fsync()` call when saving session files to ensure data is written to disk immediately.

### Documentation

- **New SAFETY_CHECKS.md**: Auto-generated documentation of all 70+ safety checks categorized by type:
  - Sanity Checks (Pre-Execution): 20
  - Directive Validation (Post-LLM): 7
  - File Validation: 4
  - Workflow State Validation: 8
  - Command Building Validation: 3
  - Input Validation: 29
  - Post-Processing Corrections: 4

- **New PROGRAM_CONFIG_ROBUSTNESS.md**: Implementation plan for making program configuration more robust with sensible defaults, validation warnings, and dry_run_file_selection mode.

### Files Changed

- `programs/ai_agent.py` - CarryOn fix, simplified verbosity, fsync
- `agent/workflow_engine.py` - Fixed `_has_placed_model()` to exclude search_models
- `agent/command_builder.py` - Added priority_patterns/prefer_patterns support
- `agent/best_files_tracker.py` - CarryOn fix with valuable_output_patterns
- `agent/directive_extractor.py` - Enhanced LLM prompt, added `_fix_ligand_workflow_conflict()`
- `agent/session.py` - Added fsync to save()
- `phenix_ai/log_parsers.py` - Added file existence retry in resolve_file_path()
- `knowledge/programs.yaml` - Added input_priorities for pdbtools
- `generate_safety_docs.py` - New documentation generator
- `docs/SAFETY_CHECKS.md` - New auto-generated safety documentation
- `docs/implementation/PROGRAM_CONFIG_ROBUSTNESS.md` - New implementation plan

---

## Version 42 (January 2025)

### Testing Infrastructure

- **Converted all tests to cctbx-style**: Migrated 8 test files (300+ tests) from `unittest.TestCase` to plain functions with fail-fast behavior
  - Matches PHENIX/cctbx testing conventions
  - First assertion failure stops with full traceback
  - Simpler syntax without class wrappers

- **New test utilities module** (`tests/test_utils.py`):
  - 20+ assert helper functions (`assert_equal`, `assert_in`, `assert_true`, etc.)
  - `run_tests_with_fail_fast()` for cctbx-style test execution
  - Supports both plain functions and TestCase classes for gradual migration

- **New testing documentation** (`docs/guides/TESTING.md`):
  - Complete guide to writing and running tests
  - Migration guide from unittest to cctbx style
  - Best practices and conventions

### Bug Fixes

- **`sanitize_for_transport()` now handles all types**: Previously converted dicts/lists to string representation. Now recursively sanitizes nested structures while preserving their types.
  - Dicts → sanitized dicts
  - Lists → sanitized lists
  - None/int/float/bool → passed through unchanged
  - Strings → sanitized (control chars, tabs, markers removed)

- **`encode_for_rest()` handles non-string input**: Added type checking to JSON-encode dicts/lists before REST encoding, preventing AttributeError on dict input.

- **`validate_directives()` preserves `file_preferences`**: Added support for boolean preferences (`prefer_anomalous`, `prefer_unmerged`, `prefer_merged`) which were previously dropped.

- **`_fix_program_name()` expanded aliases**: Added mappings for:
  - `sharpen_map`, `auto_sharpen` → `phenix.map_sharpening`
  - `build_model`, `buildmodel` → `phenix.map_to_model`

### New Directive Patterns

- **Atom type extraction**: "use selenium", "Se-SAD", "sulfur SAD" → sets `atom_type` in autosol settings
- **File preferences**: "use anomalous data", "prefer unmerged data" → sets `file_preferences`
- **Workflow preferences**: "skip autobuild", "avoid ligandfit" → adds to `skip_programs` list
- **Stop after refinement**: "stop after refinement" now works (previously required "stop after THE FIRST refinement")

### Files Changed

- `agent/transport.py` - Fixed `sanitize_for_transport()` and `encode_for_rest()`
- `agent/directive_extractor.py` - Added new patterns, fixed `validate_directives()`, expanded `_fix_program_name()`
- `tests/test_utils.py` - New assert helpers and test runner
- `tests/test_*.py` - Converted to cctbx-style (8 files)
- `docs/guides/TESTING.md` - New testing documentation
- `docs/README.md` - Updated testing section

---

## Version 41 (January 2025)

### Enhancements
- **Enhanced session summary**: Improved AI-generated summary to include:
  - **Key Output Files**: Now shows file name, type, and descriptive text explaining what each file contains (e.g., "Refined atomic model (X-ray)", "Structure factors with R-free flags and map coefficients")
  - **Key Metrics**: Enhanced prompt requests specific metric values and names; added extraction for Ramachandran outliers, rotamer outliers, and MolProbity score
  - Output files table now formatted with File/Type/Description columns
  - LLM prompt explicitly requests formatted metrics list with values

- **Multi-step workflow support**: Added `start_with_program` directive for handling sequences like "run polder then refine"
  - When user specifies "calculate polder map and then run refinement", the system extracts `start_with_program: phenix.polder`
  - This tells the workflow "run this program first, then continue with normal workflow"
  - Different from `after_program` which means "run this and stop"
  - Cleaner semantics than a `required_programs` list

- **Fixed directive override behavior**: Safer attempt-based strategy
  - First attempt (attempt_number=0): Honor user's directive value (respect explicit request)
  - Retry attempts (attempt_number>0): Trust LLM's interpretation (it may be correcting syntax)
  - This is safer than always trusting one or the other:
    - User's explicit request gets a fair chance first
    - If it fails, LLM can try to fix potential syntax issues
  - Example: User says `selection=solvent molecule MES 88` (invalid syntax)
    - Attempt 0: Uses user's value → fails
    - Attempt 1: LLM interprets as `selection=resname MES and resseq 88` → succeeds

- **Fixed fallback program selection**: Fallback now respects `start_with_program` directive
  - Previously, if LLM failed 3 times, fallback would pick the first valid program (often xtriage)
  - Now fallback prioritizes `start_with_program` if set by directive

### Bug Fixes
- **phenix.polder workflow integration**: Fixed issue where LLM incorrectly assumed polder needs map coefficients from prior refinement
  - Added polder to PROGRAM REFERENCE in LLM system prompt with clear documentation that it works with standard MTZ data (Fobs + R-free flags)
  - Added polder to both `refine` and `validate` phases in workflows.yaml
  - Explicit clarification: "does NOT need pre-calculated map coefficients or phases"

- **Generic PDB file categorization**: Fixed critical bug where unclassified PDB files (e.g., `1aba.pdb`) were being categorized as `search_model` instead of `model`
  - Changed `unclassified_pdb.parent_category` from `search_model` to `model` in file_categories.yaml
  - Added `*search*`, `*sculptor*`, `*chainsaw*` to excludes list to prevent search model files from being miscategorized
  - This ensures generic PDB files are treated as positioned models ready for refinement
  - Files explicitly named as search models (e.g., `search_model.pdb`, `template.pdb`) correctly go to `search_model` category
  - Previously, providing a simple PDB file would cause the workflow to think Phaser/MR was needed
  - Now the workflow correctly recognizes the model is already placed and allows refinement/validation programs

### New Tests
- Added workflow configuration tests for polder (TestPolderWorkflowConfig)
- Added LLM prompt tests for polder (TestPolderLLMPrompt)
- Added file categorization tests (TestUnclassifiedPDBCategorization)
- Tests verify:
  - Polder is in correct workflow phases
  - Prompt clarifies polder doesn't need phases
  - Generic PDB files categorize as `model` not `search_model`

### Files Changed
- `knowledge/prompts_hybrid.py` - Added phenix.polder to VALIDATION PROGRAMS section
- `knowledge/workflows.yaml` - Added phenix.polder to xray refine and validate phases
- `knowledge/file_categories.yaml` - Changed unclassified_pdb parent_category to 'model', added excludes
- `agent/session.py` - Enhanced `_get_final_output_files()` with descriptions, added `_describe_output_file()`, enhanced `_extract_final_metrics()` with more metrics, updated LLM summary prompt
- `agent/directive_extractor.py` - Added `start_with_program` extraction for multi-step workflows
- `agent/directive_validator.py` - Added attempt-based override behavior, `validate_intent()` now accepts `attempt_number`
- `agent/workflow_engine.py` - Added `start_with_program` handling in `_apply_directives()`
- `agent/graph_nodes.py` - Pass `attempt_number` to `validate_intent()`, fallback respects `start_with_program`
- `tests/test_new_programs.py` - Added TestPolderWorkflowConfig, TestPolderLLMPrompt, TestUnclassifiedPDBCategorization
- `tests/test_workflow_state.py` - Fixed test_dock_in_map_option to use clear search model filename
- `docs/guides/USER_DIRECTIVES.md` - Added `start_with_program` docs, attempt-based override docs, multi-step example

---

## Version 40 (January 2025)

### New Features
- **USER_REQUEST_INVALID event**: When user requests a program that's not available (e.g., "run phenix.xxx"), the agent now displays a prominent warning explaining why the request can't be fulfilled and what will run instead
- Warning is shown at QUIET verbosity level (always visible)
- Distinguishes between "unknown program" and "wrong workflow state"

### Files Changed
- `agent/event_log.py` - Added USER_REQUEST_INVALID event type
- `agent/event_formatter.py` - Added formatter for prominent warning display
- `agent/graph_nodes.py` - Emit event when user request detected as invalid

---

## Version 39 (January 2025)

### Bug Fixes
- **Event transport plumbing**: Fixed events not flowing through in two edge cases:
  1. Single-shot mode via `run_job_on_server` - events now decoded from server response
  2. API result retrieval via `get_results_as_JSON()` - events now serialized in output_files

### Files Changed
- `programs/ai_agent.py` - Added events serialization in `_build_output_files_from_history`
- `programs/ai_agent.py` - Added events decoding in `run_job_on_server`

---

## Version 38 (January 2025)

### Event System Phase 4: Display Integration
- Added `verbosity` parameter to `phenix.ai_agent` command
- Integrated EventFormatter for consistent output formatting
- Added `_display_cycle_events()` method for event rendering
- Legacy fallback when events not available

### Files Changed
- `programs/ai_agent.py` - Verbosity parameter, EventFormatter integration

---

## Version 37 (January 2025)

### Event System Phase 3: Transport Integration
- Events included in v2 API response schema
- LocalAgent and RemoteAgent parse events from responses
- Events stored in history_record for persistence

### Files Changed
- `phenix_ai/run_ai_agent.py` - Include events in response
- `phenix_ai/local_agent.py` - Parse events from response
- `agent/api_schema.py` - Updated response schema

---

## Version 36 (January 2025)

### Event System Phase 2: Decision Point Instrumentation
- All graph nodes now emit structured events
- Full LLM reasoning captured without truncation
- File selection reasons tracked

### Files Changed
- `agent/graph_nodes.py` - Event emission in perceive, plan, build nodes

---

## Version 34 (January 2025)

### Event System Phase 1: Core Infrastructure
- Created EventLog class for structured logging
- Created EventFormatter for human-readable output
- Defined 17 event types with verbosity levels
- LangGraph state compatibility (list of dicts)

### New Files
- `agent/event_log.py` - EventLog class, EventType constants
- `agent/event_formatter.py` - EventFormatter class

---

## Version 33 (January 2025)

### Cleanup and Production Hardening
- Removed deprecated state.md files
- Removed redundant backup files
- Fixed program registration after import changes
- Updated test suites for new structure

---

## Version 32 (January 2025)

### Pattern Centralization
- Moved all regex patterns to `knowledge/patterns.yaml`
- Created PatternManager for centralized access
- Updated log_parsers.py to use PatternManager

### New Files
- `knowledge/patterns.yaml` - Centralized regex patterns
- `agent/pattern_manager.py` - Pattern loading and compilation

---

## Version 31 (January 2025)

### Unified Command Builder
- Single CommandBuilder class for all programs
- Reads program definitions from YAML
- Consistent file selection across all programs
- Strategy flags and defaults from YAML

### Files Changed
- `agent/command_builder.py` - Complete rewrite

---

## Version 30 (January 2025)

### File Categorization Consolidation
- Centralized file categorization in `file_categorization.py`
- Semantic categories: model vs search_model distinction
- Categories defined in `file_categories.yaml`

### New Files
- `knowledge/file_categories.yaml`
- `agent/file_categorization.py` - Centralized categorization

---

## Version 29 (January 2025)

### BestFilesTracker
- New class to track best file of each type across cycles
- Scores based on metrics (R-free, resolution)
- R-free flag locking after first refinement

### New Files
- `agent/best_files_tracker.py`

---

## Version 28 (January 2025)

### YAML Configuration System
- Programs defined in `programs.yaml`
- Workflows defined in `workflows.yaml`
- Metrics defined in `metrics.yaml`
- Transport rules defined in `transport.yaml`

### New Files
- `knowledge/programs.yaml`
- `knowledge/workflows.yaml`
- `knowledge/metrics.yaml`
- `knowledge/transport.yaml`
- `knowledge/yaml_loader.py`

---

## Version 25-27 (December 2024)

### User Directives System
- Natural language directive parsing
- Stop conditions: "stop after X", "stop when metric < Y"
- Workflow preferences: "skip program", "prefer program"
- Four-layer stop condition checking

### New Files
- `agent/directive_extractor.py`
- `agent/directive_validator.py`
- `docs/guides/USER_DIRECTIVES.md`

---

## Earlier Versions

### Initial Development (2024)
- LangGraph pipeline architecture
- LLM integration (Claude, Gemini)
- Rules-only fallback mode
- Local and remote execution modes
- Session tracking and history
- Sanity checking system
