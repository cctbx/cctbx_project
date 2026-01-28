# PHENIX AI Agent - Changelog

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
