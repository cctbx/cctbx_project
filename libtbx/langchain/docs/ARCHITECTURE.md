# PHENIX AI Agent Architecture

## System Overview

The PHENIX AI Agent is an automated crystallographic workflow system
that operates at two levels:

1. **Strategic planner** (v114) — produces a multi-stage plan at session
   start, evaluates progress at stage gates after each cycle, retreats
   when strategies fail, and generates crystallographer-level commentary.
2. **Reactive execution engine** — analyzes logs, decides the next
   program, executes it, tracks results, and repeats until the structure
   is solved.

The strategic planner communicates with the reactive engine through
directives — the same interface a human user would use. This means the
reactive agent's safety checks always apply and the planner can be
disabled without code changes.

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                              CLIENT                                      │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                         ai_agent.py                                 ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │   Session   │  │    Agent    │  │   Command   │  │    Log     │ ││
│  │  │   Tracker   │  │  Interface  │  │  Executor   │  │   Parser   │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  │        │                │                │                │        ││
│  │        ▼                ▼                ▼                ▼        ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │ BestFiles   │  │ Local/Remote│  │  easy_run   │  │ Metrics    │ ││
│  │  │ Tracker     │  │   Agent     │  │  subprocess │  │ Extraction │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  └─────────────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ v2 JSON Request
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                    STRATEGIC PLANNER (v114)                               │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐  ┌──────────────────┐ │
│  │   Plan     │  │   Gate     │  │ Structure  │  │  Explanation     │ │
│  │ Generator  │  │ Evaluator  │  │   Model    │  │  Engine          │ │
│  └──────┬─────┘  └──────┬─────┘  └──────┬─────┘  └────────┬─────────┘ │
│         │               │               │                  │          │
│  ┌──────┴───────┐  ┌────┴──────┐  ┌─────┴─────┐           │          │
│  │ Plan Schema  │  │ Hypothesis│  │ Validation │           │          │
│  │ + Templates  │  │ Evaluator │  │ History    │           │          │
│  └──────────────┘  └───────────┘  └───────────┘           │          │
│                                                            │          │
│              directives + advice + commentary ─────────────┘          │
└─────────────────────────────┬───────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         DECISION ENGINE                                  │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                       run_ai_agent.py                               ││
│  │  ┌─────────────────────────────────────────────────────────────────┐││
│  │  │                      LangGraph                                  │││
│  │  │  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐           │││
│  │  │  │PERCEIVE │─▶│  THINK  │─▶│  PLAN   │─▶│  BUILD  │──┐       │││
│  │  │  └─────────┘  └─────────┘  └────┬────┘  └─────────┘  │       │││
│  │  │                                 ▲                     ▼       │││
│  │  │                                 └─── retry < 3 ──VALIDATE     │││
│  │  │                                                       │       │││
│  │  │                                          ┌──────────┐         │││
│  │  │                                          │ FALLBACK │──▶OUTPUT│││
│  │  │                                          └──────────┘         │││
│  │  └─────────────────────────────────────────────────────────────────┘││
│  └─────────────────────────────────────────────────────────────────────┘│
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                        Knowledge Layer                              ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │   YAML      │  │   Rules     │  │  Templates  │  │  Workflow  │ ││
│  │  │  Programs   │  │  Selector   │  │   Builder   │  │   State    │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  └─────────────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────────┘
```

## Unified Request Flow

Both LocalAgent and RemoteAgent use the same v2 JSON API and **identical transport encoding**:

```
┌─────────────────────────────────────────────────────────────────┐
│                         Client Side                              │
├─────────────────────────────────────────────────────────────────┤
│  LocalAgent / RemoteAgent                                        │
│       │                                                          │
│       ▼                                                          │
│  build_request_v2()                                              │
│       │                                                          │
│       ▼                                                          │
│  prepare_request_for_transport()                                 │
│       │                                                          │
│       ├── 1. sanitize_request() [YAML-driven]                    │
│       │       ├── Remove ZZxxZZ markers                          │
│       │       ├── Truncate long quoted strings                   │
│       │       ├── Replace tabs with spaces                       │
│       │       └── Remove control characters                      │
│       │                                                          │
│       ├── 2. json.dumps()                                        │
│       │                                                          │
│       └── 3. encode_for_rest() [text_as_simple_string]           │
│                    │                                             │
│     ┌──────────────┴──────────────┐                              │
│     │                             │                              │
│  LocalAgent:                 RemoteAgent:                        │
│  process_request_from_      send encoded to                      │
│  transport() locally        REST server                          │
│     │                             │                              │
│     └──────────────┬──────────────┘                              │
└────────────────────│────────────────────────────────────────────┘
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                      run_ai_agent.py                             │
├─────────────────────────────────────────────────────────────────┤
│  process_request_from_transport()                                │
│       │                                                          │
│       ├── 1. decode_from_rest() [simple_string_as_text]          │
│       └── 2. json.loads()                                        │
│                    │                                             │
│                    ▼                                             │
│              LangGraph Execution                                 │
│                    │                                             │
│                    ▼                                             │
│  prepare_response_for_transport()                                │
│       │                                                          │
│       ├── 0. inject_agent_build()  [v119.H2; success+error]      │
│       ├── 1. sanitize_response()                                 │
│       ├── 2. json.dumps()                                        │
│       └── 3. encode_for_rest()                                   │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Decision**: LocalAgent performs the full encode/decode roundtrip locally, even though it could skip encoding. This ensures that transport bugs are caught during local testing rather than only appearing in production server scenarios.

**v119.H2 injection point**: `inject_agent_build()` runs at the
top of `_build_group_args_response` before serialization, on
both the success and error paths.  Every response carries
`agent_build = {version, defaults_fingerprint, started_at}` —
see KDD 23.

## Component Responsibilities

### Client Components

#### ai_agent.py
Main entry point and execution loop:
- Manages the iterative workflow loop (prepare state → call graph → execute → repeat)
- Coordinates GUI callbacks and session management
- Handles command execution via PHENIX subprocess
- Post-execution result routing (`_handle_execution_result`, `_handle_failed_execution`)
- Duplicate command detection and retry (via `_handle_duplicate_check` with graph re-query)
- Client-only file injection (`_inject_missing_required_files`, needs `os.path.exists`)
- Result file exclusion for tutorials (`_load_result_file_exclusions`) (v115)
- Plan generation skip for task intent (v115)
- Solve-mode README preprocessing bypass (v115)
- Unsupported program detection in READMEs with expanded ignore set (v115)

**Not responsible for** (moved to graph in v112.66–112.69):
- Command sanitization, user param injection, crystal symmetry injection (→ BUILD)
- Stop decisions: hard stops (after_cycle, metrics targets) in PERCEIVE; after_program minimum-run guarantee in PLAN (or hard stop when `stop_after_requested=True`, v116.x); consecutive-program cap (→ PERCEIVE)
- Duplicate retries bypass (deleted `_retry_duplicate`, now uses graph)

**GUI mode execution caveat (v112.78):** In GUI mode, `_execute_sub_job_for_gui`
runs each program in a dedicated subdirectory (`sub_NN_program/`) and restores CWD
to the parent agent directory afterward.  Any code that runs after the sub-job
returns MUST NOT use `os.getcwd()` to locate output files — it will point to the
parent directory, not the sub-job output.  `_execute_sub_job_for_gui` returns
`gui_output_dir` as the 4th element of its return tuple; callers must use this
for file scanning (`_record_command_result`, `_track_output_files`) while keeping
`os.getcwd()` only for writing the agent's own log files.

#### command_postprocessor.py
Server-safe command transforms called by the BUILD node (new in v112.66):
- `sanitize_command()` — Rules A–D: strip placeholders, blacklisted params, hallucinated cross-program params, bare unscoped params. Rule B2 (v112.72): validates `space_group=` values and strips non-space-group words. Rule B3 (v120): strips an LLM-authored `unit_cell=`/`space_group=` that conflicts with the reflection file's own crystal symmetry.
- `inject_user_params()` — append user key=value params missing from command (scope-matched for dotted keys, strategy_flags-validated for bare keys)
- `inject_crystal_symmetry()` — append unit_cell/space_group from directives (validates with `_is_valid_space_group()`; v120: skips injecting a directive cell/space group that conflicts with the reflection file's crystal symmetry)
- `inject_program_defaults()` — append defaults from programs.yaml if missing (safety net)
- `postprocess_command()` — single entry point calling all four in order
- `_is_valid_space_group()` — validates that a value is a plausible space group symbol (H-M notation or IT number 1-230), rejecting English words like "determination"

All functions take explicit data arguments (no `self`/class dependencies), making
them callable from both the graph (server-side) and `ai_agent.py` (client-side
replay path).

**Rule D consistency (v112.70):** `inject_user_params` now validates bare (undotted)
keys against the program's `strategy_flags` allowlist before injection, mirroring
Rule D in `sanitize_command`. Without this, `sanitize_command` would strip a
hallucinated param like `d_min=2.5` and then `inject_user_params` would re-add it
from the user advice text.

**Strategy-flag alias awareness (v112.75):** `inject_user_params` builds an alias
map from `strategy_flags` key→flag mappings.  When `wavelength` maps to
`autosol.lambda={value}`, the alias leaf is `lambda`.  The duplicate check now
verifies both the bare key (`wavelength`) AND the alias leaf (`lambda`) against the
command string.  This prevents re-injection of `wavelength=0.9792` when
`autosol.lambda=0.9792` is already present — the prior check only looked for
"wavelength" and missed the aliased form.

**Autosol atom_type dedup (v112.75):** `postprocess_command` validates that
`autosol.atom_type` and `mad_ha_add_list` differ after injection.  When both are
the same element (e.g., both `S`), the duplicate `mad_ha_add_list` is stripped to
prevent the secondary scatterer from being silently lost.

**Heavier-atom-wins rule (v112.76):** After the dedup check, a deterministic
validation compares atomic numbers (Z) of `atom_type` and `mad_ha_add_list`.  If
the primary has lower Z than the secondary, they are swapped — the heavier element
is always the primary scatterer in SAD/MAD.  Uses `_ANOMALOUS_Z` table covering
27 common anomalous scatterers.  Skips swap when either element is unknown (do no
harm).  Handles multi-element `mad_ha_add_list` by swapping with the heaviest
secondary.

**`phaser_sad.atom_type` interception (v115.09b):** The LLM sometimes uses the
dotted PHIL path `phaser_sad.atom_type=S` to tell autosol's internal phaser to
search for a secondary atom. This bypasses the `strategy_flags` system and
overrides the primary `atom_type`, causing autosol to use S instead of Se. The
command builder's `_build_strategy()` now intercepts this key: if the value
differs from the primary `atom_type` and `additional_atom_types` isn't already
set, it's converted to `additional_atom_types` (preserving the LLM's intent).
If it matches the primary or `additional_atom_types` already exists, it's
stripped. Six scenarios tested.

**Rule D design tension (v112.77):** Rule D ("fail closed") strips bare params not
in a program's `strategy_flags` allowlist.  This is safe against hallucination but
blocks legitimate LLM error recovery — e.g., `rebuild_in_place=False` correctly
identified by the LLM as the fix for a sequence mismatch was silently stripped
from autobuild.  The catch-all blacklist (v112.76) handles the reverse case (bad
params that *cause* errors) but Rule D has no feedback loop for good params that
get stripped.  The current mitigation is to expand `strategy_flags` for programs
where recovery params are known.  Autobuild expanded from 3 to 6 flags
(`rebuild_in_place`, `n_cycle_build_max`, `maps_only` added in v112.77).  A
future "warn but keep" mode could rely on PHIL validation + catch-all blacklist
as safety nets.

**Unit-cell hallucination guard (v120):** Two coordinated fixes stop a
hallucinated crystal symmetry from reaching `phenix.refine` (which aborts with
"Working unit cell is not compatible").  The trigger was the directive
extractor's few-shot prompt: it contained a realistic example cell
(`116.097 116.097 44.175 90 90 120`) and space group, which the LLM copied into
directives even when the user advice stated none (observed in beta-blip,
AIAgent_8 — the bad cell sent refine into a phaser→refine→phaser loop).  *Fix A
(`directive_extractor.py`):* the prompt's example values are replaced with
syntactically-invalid placeholders (`A_LENGTH B_LENGTH C_LENGTH ...`,
`SPACE_GROUP_SYMBOL`) plus explicit "extract ONLY if the user states one"
guards, so a leak fails safe (no six numbers → injection skipped).  *Fix B
(`command_postprocessor.py`):* a data-consistency guard reads the reflection
file's own crystal symmetry via iotbx (`any_reflection_file(...)
.as_miller_arrays()[0].crystal_symmetry()`) and treats it as authoritative on
**both** command paths — `inject_crystal_symmetry` skips injecting a *directive*
cell/space group that conflicts, and `sanitize_command` Rule B3 strips an
*LLM-authored* cell/space group already in the command that conflicts.  The
original failure took the LLM-authored path, so the injection-side guard alone
was insufficient; both paths are now covered.  Agreement tolerances:
`_CELL_LENGTH_TOL_FRAC=0.005` (0.5% on a,b,c), `_CELL_ANGLE_TOL_DEG=0.5`.  Helper
`_find_reflection_file` returns `"AMBIGUOUS"` when a command has multiple
reflection files (→ skip the check, do no harm); any iotbx read failure also
falls back to current behavior.  Verified end-to-end on the real beta-blip
`PHASER.1.mtz`: the bogus 116.097 cell is stripped while the matching `P 32 2 1`
space group is correctly kept.

#### Session Tracker (session.py)
Persists workflow state across cycles:
- Experiment type (X-ray/cryo-EM)
- Resolution
- R-free MTZ path
- Cycle history
- **User Directives** (extracted structured instructions)
- **Strategy memory** (v113) — accumulated scientific understanding from
  the thinking agent, persisted via `session.data["strategy_memory"]` and
  round-tripped through `build_session_state()` each cycle
- **Supplemental file discovery** — `_rebuild_best_files_from_cycles` (session load)
  and `record_result` (live path) call `_find_missing_outputs` to discover companion
  output files (e.g., map coefficients MTZ) and evaluate them through the best_files
  tracker. This ensures `best_files["map_coeffs_mtz"]` is populated even when the
  client doesn't track all output files.
- **Duplicate detection** — `is_duplicate_command()` uses a two-tier check: exact
  match against all prior commands, then 80%-token-overlap against successful
  commands. The overlap heuristic now compares file tokens (basenames with
  crystallographic extensions) separately — different input files means a different
  computation, regardless of parameter overlap.

#### Directive System
Extracts and enforces user instructions:
- `directive_extractor.py`: Parses natural language → structured JSON
- `intent_classifier.py`: Classifies advice into solve/solve_constrained/task/tutorial (v115)
- `directive_validator.py`: Validates LLM decisions against directives
- Integrated into graph nodes for consistent enforcement
- See [USER_DIRECTIVES.md](../guides/USER_DIRECTIVES.md) for details

#### BestFilesTracker (best_files_tracker.py)
Tracks highest-quality files:
- Scores files using YAML-based criteria
- Maintains current best by category (model, data_mtz, map_coeffs_mtz, etc.)
- Provides files for server decisions
- `STAGE_TO_PARENT` maps stage names to parent categories (e.g.,
  `refine_map_coeffs` → `map_coeffs_mtz`, `original_data_mtz` → `data_mtz`)

#### file_utils.py
Shared file classification and pattern matching:
- `classify_mtz_type()` — classifies MTZ files as `data_mtz` or `map_coeffs_mtz`
  based on filename patterns
- `matches_exclude_pattern()` — word-boundary-aware pattern matching for
  `exclude_patterns` and `prefer_patterns` in slot definitions. Patterns match at
  start-of-string or after separators (`_`, `-`, `.`), preventing false positives
  like "noligand" matching "ligand"
- `infer_experiment_type_from_files()` — file-extension-based experiment-type
  inference (v119.H18).  Returns `(type, evidence_dict)` where `type` is
  `"xray"` (any `.mtz/.sca/.hkl`, no cryo-EM extensions), `"cryoem"` (any
  `.mrc/.ccp4/.map`, no X-ray extensions), or `None` (mixed/empty).  Used by
  `directive_extractor._apply_experiment_type_program_reprints` and
  `_resolve_after_program` as the primary signal for after_program correction
  (files-win on conflict with text-based detection).  Also used by
  `plan_generator._build_context` (single source of truth — see §41).

#### mtz_inspector.py (v119.H16)
MTZ content inspection and per-program obs_labels selection:
- `inspect_mtz(path)` — reads MTZ column structure via `iotbx.mtz` and returns
  a dict describing each column group (anomalous intensities, merged intensities,
  amplitudes, R-free flags, PHIB/FOM pairs).  Side-effect-free; never raises
  (returns `{"error": ...}` on file-read failures).
- `select_obs_labels_for(program, mtz_info)` — applies a per-program preference
  policy (merged for MR, anomalous for SAD) and returns the program-specific
  obs_labels string to inject, or `None` when no ambiguity exists.
- `has_ambiguous_arrays(mtz_info)` — diagnostic helper, True iff inspection
  found multiple legitimate observation arrays.
- `_PROGRAM_PREFERENCES` — declarative table mapping each program
  (`phenix.xtriage`, `phenix.autosol`, `phenix.phaser`, `phenix.predict_and_build`)
  to its preferred array type and fallback chain.
- Consumed by `command_builder._apply_invariants()` via the
  `auto_fill_obs_labels: true` YAML invariant.  See §39.

#### workflow_state.py
File content analysis and categorization:
- `_pdb_is_small_molecule()` — reads first 8KB, returns True for HETATM-only PDB
  files (ligands, cofactors). Used to reject small molecules from model slots.
- `_pdb_is_protein_model()` — positive protein check, returns True for PDB files
  with ATOM records. Used to reject protein models from ligand slots. Returns
  False for non-existent files (safe for use as rejection filter).
- `_detect_mtz_arrays()` — detects MTZ files with multiple observation arrays
  (e.g., merged Iobs and anomalous I(+)/I(-)). Returns array labels and a
  preferred default (merged for MR, anomalous for SAD). (v115)
- `phased_data_mtz` post-categorization — scans data_mtz files for autosol
  output patterns and promotes to `phased_data_mtz` category. (v115)

#### Agent Interface (LocalAgent/RemoteAgent)
Both agents use identical interface, v2 JSON format, **and transport encoding**:
- `LocalAgent`: Full encode/decode roundtrip, then calls `run_ai_agent.run()`
- `RemoteAgent`: Encodes and sends to REST server
- Same `prepare_request_for_transport()` for encoding
- Same `history_record` response format
- Transport module: `agent/transport.py`
- Configuration: `knowledge/transport.yaml`

**Error propagation (v112.78):** `RemoteAgent._send_request()` is wrapped in a
generic `except Exception` handler.  Since `Sorry` inherits from `Exception`,
fatal server errors (e.g., daily usage limit from `rest/__init__.py`) were
silently caught, logged, and returned as `None`.  Fix: `except Sorry: raise`
before the generic handler lets fatal errors propagate to the GUI.

#### Log Parsers (phenix_ai/log_parsers.py)
Extracts metrics and output files from program log output:

**Key functions:**
- `extract_all_metrics(log_text, program)`: Main entry point for metric extraction
- `detect_program(log_text)`: Identifies which program generated a log
- `extract_output_files(log_text, working_dir)`: Finds output file paths
- `format_metrics_report(metrics)`: Generates "FINAL QUALITY METRICS REPORT"

**Note**: Metric extraction patterns are defined in `programs.yaml` and loaded via `metric_patterns.py`. Hardcoded extractors in log_parsers.py only handle complex cases (tables, multi-line context) that YAML patterns can't express.

**YAML `log_parsing` pattern conventions** (see `programs.yaml`):

| Field | Values | Usage |
|-------|--------|-------|
| `extract` | `first` (default), `last` | `last` required when a program emits one metric line per cycle (RSR, autobuild, phaser) |
| `pick_min` | `true` | After collecting all matches, take the minimum value. Used for xtriage `resolution` (multiple lines; smallest value = highest resolution) |
| `type` | `float`, `int`, `string`, `boolean` | Controls return type |

**Important audited patterns** (fixes applied in v112.14):
- `phenix.xtriage` `resolution`: uses `pick_min: true`. Pattern anchored with `^\s*` and negative lookbehind to prevent "Completeness in resolution range: 1" from matching. Skip group `(?:[0-9.]+\s*[-]\s*)?` handles "50.00 - 2.30" format.
- `phenix.real_space_refine` `map_cc`: uses `extract: last`. RSR emits one CC_mask line per macro-cycle; `last` gives the final value.

**Adding new program support:**
See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide.

**Future direction:** Newer PHENIX programs built on `ProgramTemplate`
expose `results_as_json()` which returns structured metrics without
log parsing. As programs adopt this, PERCEIVE can read JSON results
directly with `log_parsing` as fallback. See "Potential improvements"
in the Future Directions section.

### Centralized Configuration Modules

The agent uses a centralized YAML-driven architecture that reduces program configuration from 7 files to 2-3 files.

#### Metric Patterns (knowledge/metric_patterns.py)
Centralizes metric extraction patterns from programs.yaml:

```
programs.yaml          metric_patterns.py
     │                        │
     ▼                        ▼
log_parsing: ────────▶ extract_metrics_for_program()
  pattern:                    │
  type:              ┌────────┴────────┐
  display_name:      ▼                 ▼
                log_parsers.py    session.py
```

**Key functions:**
- `get_all_metric_patterns()`: Load all patterns from YAML (cached)
- `extract_metrics_for_program(log_text, program)`: Extract metrics using YAML patterns
- `format_metric_value(program, metric, value)`: Format using YAML config

#### Program Registration (knowledge/program_registration.py)
Reads `done_tracking` blocks from programs.yaml for workflow done flags:

```
programs.yaml             program_registration.py
     │                            │
     ▼                            ▼
done_tracking:           get_program_done_flag_map()  ← ALL programs
  flag: "xxx_done"       get_trackable_programs()     ← strategy: "run_once"
  strategy: "run_once"            │
  # or "count"           ┌────────┴────────┐
  # or "set_flag"        ▼                 ▼
                 workflow_state.py   workflow_engine.py
                 (_set_done_flags)   (context building)
```

**Key functions:**
- `get_program_done_flag_map()`: All programs → done flag names (from YAML)
- `get_trackable_programs()`: Programs with `strategy: "run_once"` (filtered from valid list after completion)

**`_is_program_already_done` scope (v112.75):** `_apply_directives` in
`workflow_engine.py` calls `_is_program_already_done()` before re-adding programs
from `program_settings` directives.  This function now checks two conditions:
(1) `run_once` programs whose done flag is set (original), and (2) any non-count
program with a program-specific done flag (e.g., `autosol_done` contains `autosol`).
Without check (2), `_apply_directives` could re-add completed programs like autosol
at the front of `valid_programs`, causing the LLM to re-run them.  Count-strategy
programs (refine, rsr, phaser) are excluded since they intentionally repeat.
Shared flags (e.g., `validation_done`) are excluded since the flag name doesn't
contain the program short name.

#### Summary Display (knowledge/summary_display.py)
Configures session summary display from metrics.yaml:

```
metrics.yaml           summary_display.py
     │                        │
     ▼                        ▼
summary_display: ─────▶ format_quality_table_rows()
  quality_table:              │
  step_metrics:        ┌──────┴──────┐
                       ▼             ▼
                session.py     session.py
             (Final Quality)  (Steps table)
```

**Key functions:**
- `get_quality_table_config()`: Load quality table row configs
- `get_step_metrics_config()`: Load per-program step metric configs
- `format_quality_table_rows()`: Format Final Quality table
- `format_step_metric()`: Format a step's key metric

#### Program Validator (agent/program_validator.py)
Validates that programs are fully configured:
```bash
python agent/program_validator.py phenix.map_symmetry  # Check one program
python agent/program_validator.py --all                 # Check all programs
```

### Decision Engine Components

#### run_ai_agent.py
Single entry point for all requests:
- Parses v2 JSON requests
- Validates request schema
- Creates graph initial state
- Invokes LangGraph
- Returns v2 JSON response

#### LangGraph Pipeline
Decision-making workflow:

```
PERCEIVE ──┬──▶ THINK ─▶ PLAN ─▶ BUILD ─▶ VALIDATE ──┬──▶ OUTPUT ─▶ END
           │              │        │         │    │    │
           │              │        │         │    │    │
        Analyze    Expert Select   Build   Validate │  Format
        inputs   reasoning program command response │  decision
           │    (optional)                    │     │
           │                        retry <3 │     │ retry >=3
           │                            ┌────┘     │
           │                            ▼          ▼
           │                           PLAN     FALLBACK ─▶ OUTPUT ─▶ END
           │
           └──▶ OUTPUT ─▶ END  (if red flag abort)
```

#### Decision Flow Architecture

The decision flow follows a clean, layered architecture where each component has a single responsibility:

```
┌─────────────────────────────────────────────────────────────────┐
│                     SESSION START                                │
│  directive_extractor.py extracts structured directives          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 1: WORKFLOW ENGINE                         │
│  workflow_engine.py determines valid_programs                    │
│                                                                  │
│  Input: phase, context, directives                              │
│  Output: valid_programs list (includes STOP if appropriate)     │
│                                                                  │
│  Responsibilities:                                               │
│  - Get programs from workflow phase definition                   │
│  - Add directive-required programs (after_program target,        │
│    as a MIN-RUN guarantee — see Layer 4 note on v112.78)        │
│  - Add start_with_program to front of list (multi-step flows)   │
│  - Add STOP if skip_validation=true                             │
│  - When after_program target has run:                            │
│      stop_after_requested=True  → WIPE to [STOP] (v116.x)        │
│      stop_after_requested=False → do nothing, plan continues     │
│  - Apply workflow preferences (skip/prefer)                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│             LAYER 1.5: EXPERT REASONING (THINK, v113)            │
│  agent/thinking_agent.py — active at basic/advanced/expert       │
│                                                                  │
│  Input: log_text, history, strategy_memory                      │
│  Output: enriched user_advice with expert guidance              │
│                                                                  │
│  Responsibilities:                                               │
│  - Analyze program logs with domain expertise                   │
│  - Detect crystallographic issues (twinning, stalling, etc.)    │
│  - Inject guidance into user_advice for PLAN                    │
│  - Can recommend STOP (sets command="STOP")                     │
│  - Forward error_classification/failure_count to context (v115) │
│                                                                  │
│  NOT responsible for:                                            │
│  - Program selection (that's PLAN's job)                        │
│  - Parameter settings (presents evidence, not commands)         │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 2: LLM DECISION (PLAN)                     │
│  graph_nodes.py plan() function                                  │
│                                                                  │
│  Input: valid_programs, context, history                        │
│  Output: chosen program (must be in valid_programs)             │
│                                                                  │
│  Responsibilities:                                               │
│  - LLM selects best program from valid_programs                 │
│  - Simple validation: is choice in valid_programs?              │
│  - If invalid choice, use first valid program                   │
│                                                                  │
│  NOT responsible for:                                            │
│  - Directive logic (already handled in Layer 1)                 │
│  - Stop condition evaluation (handled in Layer 4)               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 3: BUILD & EXECUTE                         │
│  graph_nodes.py build(), ai_agent.py execute                    │
│                                                                  │
│  Input: chosen program, files, context                          │
│  Output: command result, output files                           │
│                                                                  │
│  v115 guards (before command assembly):                         │
│  - PHIL validation (strip unrecognized/blocked params)          │
│  - MTZ label injection (obs_labels for multi-array MTZ)         │
│  - AutoSol sites estimation (from sequence)                     │
│  - map_sharpening resolution injection (from mtriage)           │
│  - AutoBuild rebuild_in_place=False when R-free > 0.45          │
│  - Polder selection forwarding from directives                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 4: POST-EXECUTION CHECK                    │
│  ai_agent.py _run_single_cycle() - SINGLE PLACE                 │
│                                                                  │
│  Input: completed program, cycle number, directives             │
│  Output: should_stop (bool), stop_reason                        │
│                                                                  │
│  Responsibilities:                                               │
│  - Check if after_cycle condition met (hard stop)               │
│  - Check if metric targets met — r_free, map_cc (hard stop)    │
│  - This is the ONLY place hard stop conditions are evaluated    │
│                                                                  │
│  NOTE (v112.78): after_program is intentionally NOT a hard      │
│  stop here.  It is a minimum-run guarantee: PLAN suppresses     │
│  auto-stop until the target program has run, but the LLM        │
│  decides when to actually stop.  This prevents premature        │
│  termination on multi-goal requests where the directive          │
│  extractor (or plan-to-directives, per stage) can only name     │
│  one program.                                                    │
│                                                                  │
│  HISTORICAL NOTE (v116.x routing): For a period between          │
│  v112.78 and v116.x, the consumer side in                       │
│  workflow_engine._apply_directives drifted from this design     │
│  and wiped valid_programs to [STOP] once the after_program      │
│  target had run.  v116.x distinguishes user-explicit stops      │
│  ("refine and stop") from plan-progression hints (per-stage     │
│  after_program from plan_to_directives) via a new directive     │
│  field, stop_conditions.stop_after_requested.  Set by the       │
│  directive extractor (agent/directive_extractor.py)             │
│  only when the user/README explicitly requested a stop-after    │
│  condition.  Workflow_engine wipes for True, does nothing for   │
│  False.  Prompts_hybrid.py also gates on the flag and emits     │
│  "Stop target: X" prompt lines only for True.  See section      │
│  "Stop-after directive routing (v116.x)" below.                 │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Principles:**

1. **Single Responsibility**: Each layer handles one aspect of decision-making
2. **Clear Data Flow**: Directives flow from extraction → workflow engine → validated programs
3. **Post-Execution Stop**: Stop conditions are only checked after program execution, not during planning
4. **Simple Validation Gate**: The plan() function validates that LLM choice is in valid_programs list

#### Knowledge Layer

**YAML Programs** (`knowledge/programs.yaml`):
- Program definitions with input/output specifications
- `strategy_flags` — valid PHIL parameters per program (used by phil_validator)
- `input_priorities` — file category preferences per input slot
- `obs_labels` for multi-array MTZ handling (xtriage, autosol, phaser, predict_and_build)
- `ligand_cif` optional input for refine
- Autosol output typed as `phased_data_mtz` (not `data_mtz`)
- Invariants (rules that must be satisfied): auto_fill_resolution for
  programs that need it (real_space_refine, dock_in_map, map_to_model,
  polder, map_correlations, validation_cryoem, predict_and_build)
- `exclude_patterns` on input slots for file filtering (e.g., pdbtools
  ligand slot excludes `pose` to prevent selecting individual LigandFit
  pose files instead of the final combined model)
- Scoring rules

**Rules Selector** (`rules_selector.py`):
- Matches workflow state to applicable programs
- Filters by experiment type
- Applies priority rules

**Command Builder** (`command_builder.py`):
- **Single entry point** for all command generation
- 4-stage pipeline:
  1. `_select_files()` - Priority-based file selection
  2. `_build_strategy()` - Strategy with output_prefix
  3. `_apply_invariants()` - Auto-fill resolution, R-free flags
  4. `_assemble_command()` - Final command string
- Uses `CommandContext` dataclass for all parameters
- Replaces fragmented logic from template_builder and graph_nodes
- v115.05: Phaser copies injection reads `log_analysis["n_copies"]`
  from the current cycle (same-cycle), not just `session_info["asu_copies"]`
  from the previous cycle (1-cycle delay). Eliminates the off-by-one
  where the first post-Phaser autobuild missed the copies count.

**Template Builder** (`template_builder.py`):
- Legacy interface (delegates to CommandBuilder)
- YAML-to-slot mapping
- Invariant checking

**Workflow State** (`workflow_state.py`):
- Determines current workflow position
- Categorizes available files
- v115.05: `_anti_ligand_patterns` excludes filenames matching
  `no_ligand` from ligand file classification. Orphaned PDB files
  in the `pdb` category that don't match any model subcategory
  (refined, phaser_output, etc.) and aren't ligands or intermediates
  are promoted to `model` — fixes user-supplied input PDBs (e.g.
  `1aba.pdb`) not appearing as `has_model=True` in PERCEIVE
- Half-map pair detection: numbered CCP4 file pairs (e.g.
  `map_1.ccp4`, `map_2.ccp4`) are promoted to `half_map` when a
  companion full map exists alongside them. Without a companion,
  they stay in `full_map` (prevents false positives on
  resolve_cryo_em segmented outputs). Explicit `half` in the
  filename always takes priority.
- MR solution detection: PDB files with `mr_solution` in the name
  are categorized as `phaser_output`, which triggers
  `model_is_placed` in the plan context (skips unnecessary MR).
  The pattern requires the `mr_` prefix to avoid matching
  `nmr_solution`.
- Sigma-A MTZ handling: files with `sigmaa` in the name stay in
  `data_mtz` (they contain both Fobs and map coefficients). They
  are NOT moved to `map_coeffs_mtz` because `phenix.refine`'s
  `data_mtz` slot has `exclude_categories: [map_coeffs_mtz]` —
  dual membership would prevent refine from finding its input data.
- Detects experiment type

## LLM Provider Architecture

The agent supports five LLM providers: **google**, **openai**, **ollama**,
**anthropic** (Claude), and **portkey** (the Portkey gateway fronting Azure
OpenAI). The canonical list lives in `core/llm.SUPPORTED_PROVIDERS` (single
source of truth — `graph_nodes.py` and `programs/ai_agent.py` import it rather
than keeping their own copies; see §"Provider Routing — Single Source of Truth"
below). The provider is set via `params.communication.provider` (or the
`LLM_PROVIDER` env var) and flows consistently through every LLM call point —
there is no cross-provider fallback.

### LLM Call Points

The system uses three distinct LLM roles, each with provider-specific model defaults:

| Activity | Function | ollama | google | openai | anthropic | portkey |
|---|---|---|---|---|---|---|
| **Agent decisions** (PLAN node) | `get_planning_llm()` → `get_expensive_llm()` | qwen3:32b (json_mode) | gemini-2.5-pro | gpt-5.5 | claude-opus-4-7 | gpt-5 (Azure) |
| **Log summarization** | `get_cheap_llm()` | qwen2.5:7b | gemini-2.5-flash-lite | gpt-5.4-nano | claude-sonnet-4-6 | gpt-5 (Azure) |
| **RAG analysis** | `get_expensive_llm()` | qwen3:32b | gemini-2.5-pro | gpt-5.5 | claude-opus-4-7 | gpt-5 (Azure) |
| **Directive extraction** | `call_llm_simple()` | direct ollama HTTP | gemini-2.5-flash-lite | (provider default) | claude-sonnet-4-6 | gpt-5 (Azure) |

Model defaults come from the five role tables in `core/llm.py` and can be
overridden per provider via env var (`OLLAMA_LLM_MODEL`, `ANTHROPIC_LLM_MODEL`,
`PHENIX_PORTKEY_MODEL`).

### Provider Flow

```
params.communication.provider = "<provider>"
         │
         ├─► LocalAgent.decide_next_step()
         │       └─► request["settings"]["provider"] = "<provider>"
         │            └─► run_ai_agent.py: create_initial_state(provider=...)
         │                 └─► graph_nodes.py PLAN: get_planning_llm(provider)
         │                      └─► get_expensive_llm(provider)
         │
         ├─► _extract_directives()
         │       ├─► ollama: forces run_on_server=False (local ollama not on server)
         │       └─► run_directive_extraction(provider=...)
         │            └─► call_llm_simple(provider=...)
         │
         └─► run_ai_analysis.run() [summarization + RAG]
                 └─► setup_llms(provider=...)
                      ├─► expensive_llm = get_expensive_llm(provider)
                      ├─► cheap_llm = get_cheap_llm(provider)
                      └─► embeddings = provider-specific embeddings
```

### Key Files

| File | Role |
|------|------|
| `core/llm.py` | `get_llm_and_embeddings()`, `get_expensive_llm()`, `get_cheap_llm()` — model creation |
| `agent/graph_nodes.py` | `get_planning_llm()` — cached LLM for PLAN node decisions |
| `agent/api_client.py` | `call_llm_simple()` — lightweight direct calls (directive extraction) |
| `utils/run_utils.py` | `setup_llms()` — creates all three LLM roles for analysis; also the shared `_safe_float()` metric-coercion helper (v120) |
| `phenix_ai/local_agent.py` | Reads provider from `params.communication.provider`, passes to request |
| `phenix_ai/run_ai_agent.py` | Reads provider from settings, passes to `create_initial_state()` |

### Ollama-Specific Behavior

- **LLM-only modes** (v115.05): All analysis modes that don't need the
  RAG database — `directive_extraction`, `advice_preprocessing`,
  `failure_diagnosis`, `agent_session` — are routed to local execution
  when `run_on_server=False` or when `provider=ollama`. Only `standard`
  mode (log analysis with knowledge-base retrieval) requires the server
  database. This prevents concurrent runs from queueing at the Phenix
  server for preprocessing when a local LLM is available.
- **Planning LLM** uses `json_mode=True` (sets `format="json"` in ChatOllama)
  because ollama models need explicit JSON formatting; Google and OpenAI handle
  structured output without this flag
- **GPU requirement**: `rest_server.requires_gpu = True` is set for ollama to
  ensure server dispatch targets GPU-capable machines

### Adding a New Provider

`SUPPORTED_PROVIDERS` in `core/llm.py` is the single source of truth; every
routing point below either reads it or branches on the provider name. The
complete checklist (learned across the v120 portkey/anthropic work):

1. Add the provider to `SUPPORTED_PROVIDERS` in `core/llm.py`.
2. Add model defaults to the role tables in `core/llm.py`
   (`DECISION/RAG/RAG_EMBEDDING/EXPENSIVE_MODEL_DEFAULTS`) and, if it takes an
   env-var model override, to `_PROVIDER_MODEL_ENV_OVERRIDES`.
3. Add the provider branch to `get_llm_and_embeddings()` and `get_expensive_llm()`
   in `core/llm.py`. (`get_cheap_llm()` needs no change — it falls through to
   the factory for any provider without a `CHEAP_MODEL_DEFAULTS` entry.)
4. Add `_call_<provider>_llm()` + dispatch in `agent/api_client.py`.
5. Add the provider branch in `directive_extractor._call_llm` (handler-select,
   both import paths) and `_call_llm_fallback`.
6. Add a rate-limit handler in `agent/rate_limit_handler.py` (a dedicated one,
   not a borrowed `get_openai_handler`, so metrics stay unpolluted).
7. Add the key/config check in `graph_nodes.validate_provider()` (execution-time
   only — never assert keys at import) and the handler-select at the second
   `graph_nodes` LLM site; add the handler in `thinking_agent._get_rate_handler`.
8. Add the provider to the `provider` PHIL enums in `programs/ai_agent.py` and
   `programs/ai_analysis.py` (the GUI dropdown regenerates from these).
9. If the provider can do embeddings (RAG): add it to `run_utils.validate_api_keys`,
   `run_utils.get_db_dir_for_provider`, `analysis/summarizer.py` chunk sizes,
   `run_query_docs.py`, and the database-build tools
   `command_line/rebuild_ai_database.py` + `update_ai_database.py`. If it has no
   native embeddings (like anthropic), wire chat only and let
   `_delegate_embeddings_for_nonnative()` handle the embeddings side.
10. Add the runtime deps to `phenix_ai/install_ai_tools.csh`.

### Provider Routing — Single Source of Truth (v120 P2)

Before v120, `SUPPORTED_PROVIDERS` was duplicated as independent literals in
`agent/graph_nodes.py` and `programs/ai_agent.py` (the latter as a hardcoded
fallback list). The two could drift — a provider added to one but not the
other would validate in one layer and be rejected in another. v120 centralizes
the list in `core/llm.py`; both consumers import it via the standard two-path
`try/except ImportError` fallback. The `ai_agent.py` literal was removed
entirely: its consumer is an *error-advice* helper, so on import failure it
degrades to the provider-agnostic rules-only suggestion rather than raising a
second exception over the original error. A runtime import-identity test
(`tst_supported_providers_single_source.py`) pins that `graph_nodes` re-exports
the same object and that no literal survives outside `core/llm.py`.

### Portkey + Anthropic Specifics (v120 P3)

- **Portkey** is OpenAI-SDK-compatible. The langchain path
  (`get_llm_and_embeddings`) uses `ChatOpenAI`/`OpenAIEmbeddings` pointed at
  `PORTKEY_BASE_URL` with `PORTKEY_AZURE_API_KEY`; the raw-SDK call path
  (`api_client._call_portkey_llm`, `directive_extractor` fallback) uses the
  Portkey SDK with `provider="azure-openai"`. Both read env at call time and
  raise a clear `ValueError` (never a bare `KeyError`) when a var is missing.
  The env check precedes the SDK import so a misconfigured-but-installed
  deployment gets the friendly message.
- **Asymmetric kwargs** are handled by `sanitize_llm_kwargs(provider, kwargs)`,
  keyed per provider: portkey (Azure gpt-5) renames `max_tokens →
  max_completion_tokens` and drops `temperature`; anthropic always keeps
  `max_tokens` (`ChatAnthropic` requires it) but `temperature` is
  **model-dependent** — the Claude 4.6/4.7+ reasoning family
  (`claude-sonnet-4-6`, `claude-opus-4-7`, future 4.8+) rejects
  `temperature`/`top_p`/`top_k` with HTTP 400, so both `sanitize_llm_kwargs`
  and the `ChatAnthropic` constructor drop them for that family (detected by
  `anthropic_model_rejects_sampling_params()`, a model-name family check).
  Older Claude models still accept `temperature` and keep it. The factory's
  `ChatAnthropic(...)` only passes `temperature=` when the model accepts it.
- **Embeddings.** Anthropic has no native embeddings endpoint. The factory
  always returns a working chat LLM; the embeddings object is constructed
  lazily (no network call), so the end-of-run `agent_session` assessment —
  which constructs but never queries embeddings — is safe even with a
  no-access key. For anthropic, `_delegate_embeddings_for_nonnative()`
  delegates to OpenAI/Google when a key is present (emitting a
  `[PROVIDER_DELEGATION]` marker), else returns `None`; it never raises on
  construction. A real embeddings *query* failure (only via `standard`/RAG)
  surfaces a clear execution-time error echoing the raw upstream message,
  not a `"Unauthorized"` string match.

### FORCE_NO_AI_SERVER (v120 P1)

`programs/ai_agent.py::run()` checks `FORCE_NO_AI_SERVER`: when set to the
literal `1` (whitespace-stripped) it forces `communication.run_on_server=False`,
placed after `set_defaults()` and before session-management dispatch so it
covers all execution paths (v2-API, `iterate_agent`,
`run_job_on_server_or_locally`). The same override is mirrored in
`programs/ai_analysis.py::run()` (placed before its
`run_job_on_server_or_locally` dispatch), so the env var also forces local
execution for `phenix.ai_analysis` across every analysis_mode (standard,
agent_session, advice_preprocessing, directive_extraction, failure_diagnosis).
The shipped PHIL default stays `True` — this is the surgical local-only switch,
not a baseline change. No effect on `predict_and_build`/`predict_model` run
separately. The flag is **absolute**: both dispatchers
(`run_job_on_server_or_locally` in `ai_agent.py` and `ai_analysis.py`) re-check
it and refuse to submit to the server, even via the no-local-database fallback
that `standard` mode would otherwise take. When local execution is impossible
(standard mode, no local RAG database for the provider), the dispatcher raises
`Sorry` with actionable guidance instead of silently contacting the server; the
LLM-only modes (directive_extraction, advice_preprocessing, failure_diagnosis,
agent_session) never need the database and always run local. Because provider
travels in the request `settings` and is consumed only in shared
`agent/`+`core/` code, local and server execution run identical branches; the
only difference is whether information crosses machines. The server simply also
needs the provider's keys in its environment (deployment, not code).

## Data Flow

### Request Building (Client)

```python
# Both agents build requests the same way
session_state = {
    "resolution": 2.5,
    "experiment_type": "xray",
    "rfree_mtz": "/path/to/data.mtz",
    "best_files": {"model": "/path/to/model.pdb"},
    "bad_inject_params": {"phenix.refine": ["ignore_symmetry_conflicts"]},
    "advice_changed": False,
    "unplaced_model_cell": None,
    "strategy_memory": {},  # Thinking agent state (v113)
}

request = build_request_v2(
    files=available_files,
    cycle_number=5,
    log_content=last_log,
    history=cycle_history,
    session_state=session_state,
    user_advice=guidelines,
    ...
)
```

**`bad_inject_params` flow** (added in v112.66; static-seed in v120, see §45):
Parameters that previously caused PHENIX errors (learned) — plus
statically-known-invalid pairs (v120) — are blacklisted and propagated to BUILD
for stripping:

```
ai_agent.py: session.get_all_bad_inject_params()   # v120: merged static ∪ learned
    → session_info["bad_inject_params"]            # (was session.data[...], learned-only)
    → build_session_state() → session_state
    → build_request_v2() → normalized_session_state
    → transport encode/decode
    → run_ai_agent.py → create_initial_state(bad_inject_params=...)
    → graph state["bad_inject_params"]
    → BUILD: postprocess_command(bad_inject_params=set(...))
```

**`session_info` field plumbing contract (parity-critical).**  Any value that
`ai_agent.py` puts into the per-cycle `session_info` dict and that server-side code
(PERCEIVE / BUILD / the gate) needs to read must be carried at EVERY hop of the
shared transport path, or it is silently dropped:

```
ai_agent.py session_info["X"]
  → build_session_state()        [api_client.py]   session_info → session_state
  → build_request_v2()           [api_client.py]   session_state → normalized_session_state
                                                    (the WIRE WHITELIST — an explicit
                                                     per-field copy; unlisted keys are dropped)
  → create_request()             [api_schema.py]   dict(session_state), no filter
  → apply_request_defaults()     [api_schema.py]   additive only (fills defaults,
                                                    never strips; the session_state
                                                    "subfields" list is descriptive,
                                                    NOT a whitelist)
  → transport encode/decode      [transport.py]    whole-request JSON, no per-key filter
  → run_ai_agent.py              session_state → session_info  (explicit per-field map-back)
  → graph_nodes.py PERCEIVE      session_info.get("X")
```

This is parity-critical because **LocalAgent and RemoteAgent share this entire
path**: LocalAgent deliberately runs the same `build_request_v2()` + encode/decode
roundtrip (see "Key Design Decision" above).  The single chokepoint is the
`build_request_v2()` wire whitelist — whatever it omits is dropped IDENTICALLY for
local and server, so a missing field doesn't violate parity but does make the
feature inert.  A new `session_info` field therefore requires, at minimum: an entry
in `build_session_state()`, the `build_request_v2()` whitelist, and the
`run_ai_agent.py` map-back.  If it is read in `graph_nodes.py` via the literal
`session_info.get("X")` (not the `_si` alias), it must ALSO be registered in
`agent/contract.py::SESSION_INFO_FIELDS` (with `CURRENT_PROTOCOL_VERSION` bumped to
≥ the field's version), which `tst_contract_compliance.py` enforces.

**Plan-driven program fields + v120 transport repair.**  Three `session_info`
fields let the plan steer `valid_programs` in PERCEIVE:
`plan_has_pending_stages` (suppress AUTO-STOP while the plan has work),
`plan_next_stage_programs` (offer the next pending stage's programs when the engine
would only offer STOP), and `plan_current_unrun_lead_program` (v120 Option 2a —
offer the current ACTIVE stage's un-run lead program, e.g. `phenix.autobuild`, even
when other non-STOP programs are valid; see §38.5).  Before v120, none of the three
fully round-tripped: `plan_has_pending_stages` was carried by
`build_session_state()`/`run_ai_agent.py` but NOT the `build_request_v2()` wire
whitelist, and the other two were carried nowhere — so all three were dropped on
the wire for BOTH local and server (parity intact, but the plan-injection behavior
was inert).  v120 plumbs all three through `build_session_state()`,
`build_request_v2()`, and `run_ai_agent.py`, verified end-to-end through the real
`create_request`/`apply_request_defaults` and a JSON transport roundtrip.
`plan_current_unrun_lead_program` is registered in `agent/contract.py` as a v6
field (default `""`), which bumped `CURRENT_PROTOCOL_VERSION` 5 → 6.

**Tri-state field + per-file map + `is not None` guards (v120 → v120.2,
`input_mtz_has_rfree` / `mtz_rfree_map`).**  Most `session_info` fields are plumbed
with truthy guards (`if session_info.get("X"):`), which is fine for fields whose
falsy value is also their "absent" meaning.  The R-free fields are the exception.
The scalar `input_mtz_has_rfree` is **tri-state** (`True` = the original input MTZ
has an R-free array, `False` = confirmed none, `None` = undetermined), and `False`
is a *meaningful* value distinct from `None`: a truthy guard would drop a
confirmed-`False` and flip the generate decision.  v120.2 adds
`mtz_rfree_map` — a `{basename: bool}` map covering every local MTZ the client
could inspect (inputs **and** intermediate outputs like `PHASER.1.mtz`) — because
the scalar describes only the original input, while refinement often runs on a
derived output (see the "R-free generate guard" subsection).  An empty map is also
meaningful (distinct from absent).  Every hop therefore uses an explicit
`is not None` guard.  Both fields are client-extracted (the server cannot read the
MTZ) and registered in `agent/contract.py`; `mtz_rfree_map` (v8) bumped
`CURRENT_PROTOCOL_VERSION` 7 → 8 (the scalar was the v7 bump).  The round-trips —
the scalar's `True`/`False`/`None` and the map (intact, including `False` entries
and the empty-map case) — are pinned by `tst_input_mtz_has_rfree_plumbing.py`.

**Error pattern expansion (v112.75):** The learning system originally only triggered
on "unknown command line parameter" and "no such parameter" errors.  PHIL
boolean-type errors ("True or False value expected, scope.path.param="value" found")
are now also caught — the full PHIL path and all components ≥ 6 characters are
blacklisted.  Without this, `inject_user_params` could loop indefinitely injecting
a parameter that causes a type mismatch rather than an "unknown parameter" error
(observed as 9 wasted cycles in run 107).

**Catch-all injection blacklist (v112.76):** A supervisor pattern that handles
"unknown unknowns" — PHIL error formats not covered by pattern-based learning.
`postprocess_command` surfaces the list of params added by inject_* steps via
`return_injected=True` (opt-in, default False for backward compat).  The replay
path in ai_agent.py stores this list on `session.data["last_injected_params"]`.
`_update_inject_fail_streak()` tracks consecutive same-program failures via error
fingerprint (first 120 chars, normalized, digits stripped).  After N=2 failures
with matching fingerprints and non-empty injected list, all injected params are
blacklisted.  Recovery retries (`force_retry_program`) are excluded.

```
_get_command_for_cycle()
  ├─ session.data["last_injected_params"] = []        (init)
  ├─ _query_agent_for_command()                        (server path)
  │    └─ graph BUILD node
  │         └─ postprocess_command(return_injected=True)
  │              └─ state["last_injected_params"] = [...]
  │    └─ history_record["last_injected_params"] → session.data
  └─ postprocess_command(return_injected=True)         (replay path)
       └─ session.data["last_injected_params"] = [...]

_record_command_result()
  └─ _update_inject_fail_streak(program, error, session)
       ├─ SUCCESS → clear streak
       └─ FAILURE → fingerprint match?
            ├─ yes → count++ → count>=2? → blacklist injected → reset
            └─ no  → reset streak to count=1
```

**`best_files` value type contract:** Most entries are plain strings (`category → path`). However, multi-file entries such as `half_map` may be stored as a list:

```python
"best_files": {
    "model":    "/path/to/model.pdb",       # single file → str
    "half_map": ["/path/map1.mrc", "/path/map2.mrc"],  # two files → list
}
```

All server-side code that reads `best_files` values must handle both types. Use `CommandBuilder._best_path(value)` to safely extract a single path from either a string or a list. Do **not** pass `best_files` values directly to `os.path` functions — this caused a cycle=2 crash (`TypeError: expected str, bytes or os.PathLike, not list`) that only surfaced when `session_state` was re-sent from a client that stored `half_map` as a list.

**Programs that need both `full_map` and `half_map` simultaneously:** `phenix.map_to_model` accepts both a full map AND half maps at once. It is marked `keep_half_maps_with_full_map: true` in `programs.yaml`. The post-selection validation in `CommandBuilder._select_files()` checks this flag before removing half maps.

**Programs that prefer half-maps over full maps:** `phenix.resolve_cryo_em`, `phenix.mtriage`, and `phenix.map_sharpening` are marked `prefers_half_maps: true`. When both `full_map` and `half_map` are selected, the dedup drops the `full_map` and keeps the half-maps. This is scientifically correct (half-map FSC is the gold standard for resolution) and avoids "Maps have different dimensions" errors when the full map has different grid dimensions from post-processing.

**Supplement logic for `multiple:true` slots:** LLMs often assign only one file to a `multiple:true` slot (e.g., `half_map=file2.ccp4`). Auto-fill skips the slot because it's already "filled". A supplement loop after auto-fill checks each `multiple:true` slot and backfills missing files from the category. Uses `os.path.realpath()` for symlink robustness.

**Half-map pair detection (v115.07):** Some cryo-EM half-maps use `_1/_2` suffixes instead of containing "half" in the name (e.g., `7n8i_24237_box_1.ccp4`). The categorizer's post-processing detects these pairs using a two-tier heuristic:
- *Tier 1*: Exactly 2 `full_map` files whose basenames differ only by a trailing `_1/_2` → promote both to `half_map`. No companion full map needed (these ARE the only maps).
- *Tier 2*: ≥3 `full_map` files with a `_1/_2` pair AND at least one companion remaining → promote the pair, keep the companion.

**Reference model categorizer (v115.07):** When ≥2 PDB files are categorized as `model`, one may be a reference model for restraints (not a model to refine). A post-processing heuristic checks filenames for keywords (`reference`, `homolog`, `template`, `restraint`, `high_res`) and reclassifies the first match to `reference_model`. Agent output files (`refine_*`, `autobuild_*`, etc.) are excluded from this check. Only one file is reclassified per run. This works with the Tier 1 exclusion in `command_builder.py` (line 380), which excludes `reference_model.file` from primary model selection when the strategy dict names it explicitly.

**Orphan-map promotion (v115.07+):** Map files that end up in the `map` parent category but not in any subcategory (`full_map`, `half_map`, `optimized_full_map`) are promoted to `full_map`. This mirrors the orphan-PDB → model promotion. Root cause: the YAML `full_map` excludes list has `*_a.*` and `*_b.*` (for half-map suffixes) which false-positive on filenames like `emd-20026_auto_sharpen_A.ccp4`. Without promotion, `has_full_map=False` and programs with `requires_full_map: true` (e.g. `real_space_refine`) are never offered.

### Request Processing (Server)

```python
# run_ai_agent.py processes all requests uniformly
def run(request_json):
    request = _parse_request(request_json)
    response = _process_request(request)
    return _build_group_args_response(response)
```

### Response Handling (Client)

```python
# Both agents receive the same response format
result = agent.decide_next_step(...)
history_record = result  # Contains program, command, reasoning, etc.
command = history_record["command"]
```

## File Organization

See the [README.md](../README.md#directory-structure) for the complete directory tree.

Key directories:
- `agent/` — Core agent logic (session, graph nodes, command builder,
  workflow engine, structure model, gate evaluator, plan generator,
  hypothesis evaluator, validation history, etc.)
- `knowledge/` — YAML configuration files and supporting Python
  modules (plan schema, plan templates, explanation prompts,
  thinking prompts, expert KB, etc.)
- `phenix_ai/` — Runtime entry points (local/remote agent, log parsers)
- `programs/` — PHENIX program integration (main entry point `ai_agent.py`)
- `analysis/` — Post-run log analysis and session evaluation
- `core/` — LLM provider abstraction
- `tests/` — 55+ test files with 1400+ tests

## Key Design Decisions

### 1. Unified v2 JSON API

**Rationale**: Consistency between local and remote execution.

**Implementation**:
- Both agents use `build_request_v2()` to create requests
- Single `run_ai_agent.run(request_json)` entry point
- Identical response format for both paths

### 2. Client-Server Separation

**Rationale**: Allows server-side updates without client reinstallation.

**Implementation**:
- Client gathers inputs, executes commands, and handles results
- Server makes all decisions (program selection, command construction, stop logic)
- All command post-processing (sanitize, inject params, inject symmetry, inject
  defaults) runs in the graph's BUILD node (server-side)
- Only `_inject_missing_required_files` remains client-side (needs `os.path.exists`)
- Clean API boundary via JSON protocol

### 3. Configuration Over Code

**Rationale**: Domain experts can modify behavior without code changes.

**Implementation**:
- Program definitions in YAML
- Scoring rules in YAML
- Workflow transitions in YAML

### 4. Session State Tracking

**Rationale**: Ensures consistency across workflow cycles.

**Implementation**:
- R-free MTZ locked after first refinement
- Experiment type locked after detection
- Best files tracked with quality scores
- All state transmitted in `session_state`

### 5. LangGraph Pipeline

**Rationale**: Modular, testable decision-making.

**Implementation**:
- Discrete nodes (PERCEIVE, THINK, PLAN, BUILD, VALIDATE, OUTPUT)
- State passed between nodes
- Easy to add/modify nodes

### 5a. Single Source of Truth (v112.66–112.69 refactor)

**Rationale**: The agent had accumulated overlapping decision layers — "should we
stop?" was answered in 5 places and "fix the command" in 2 places. This caused
real bugs: 20× refine loops, premature stops, crystal symmetry mismatches.

**Implementation**: The graph is the single authority for decisions:
- **"Should we stop?"** → PERCEIVE + PLAN only
- **"Fix the command"** → BUILD only (via `postprocess_command()`)
- **"Execute and recover"** → `ai_agent.py` only

`ai_agent.py` is a pure execution loop: prepare state → call graph → execute
command → handle result. No parallel decision paths remain. Duplicate retries
flow through the graph (via `duplicate_feedback` parameter) rather than bypassing
it.

**Command post-processing pipeline** (runs in BUILD after command construction):

```
postprocess_command():
  1. sanitize_command()         — Rules A–D (strip bad params)
  2. inject_user_params()       — from user advice text
  3. inject_crystal_symmetry()  — from directives
  4. inject_program_defaults()  — from programs.yaml defaults
```

**Sanitization rules:**

| Rule | Fires when | Strips |
|------|-----------|--------|
| Probe | Program in `_PROBE_ONLY_FILE_PROGRAMS` | All key=value except file paths |
| A | Always | Blacklisted params (`bad_inject_params`) |
| B | key contains `space_group` or `unit_cell` | Placeholder values ("Not specified", etc.) |
| C | Program has zero strategy_flags | Any key=value not in `_UNIVERSAL_KEYS` |
| D | Program HAS strategy_flags, key has no dots | Bare params not in program's allowlist |

**Rule D consistency (v112.70):** `inject_user_params` mirrors Rule D's allowlist
check for bare keys extracted from user advice. Without this, a param stripped by
Rule D (e.g. `d_min=2.5` for autobuild) would be re-injected from the advice text.

### 6. Automation Modes

**Rationale**: Different users need different levels of control over the workflow.

**Implementation**:
Two modes controlled by `maximum_automation` parameter:

| Mode | Setting | `predict_and_build` Behavior |
|------|---------|------------------------------|
| **Automated** | `maximum_automation=True` (default) | Runs complete workflow (prediction → MR → building) |
| **Stepwise** | `maximum_automation=False` | Stops after prediction only |

The `automation_path` is set in workflow state and propagated to all decision points:

```python
# workflow_engine.py
context["automation_path"] = "automated" if maximum_automation else "stepwise"

# workflow_engine.py (_check_program_prerequisites)
# Block after full workflow has completed (both automated and stepwise)
if context.get("predict_full_done"):
    return False
# Block in stepwise mode after prediction-only step
if automation_path == "stepwise" and context.get("predict_done"):
    return False
```

**Stepwise workflow example**:
```
xtriage → predict_and_build(stop_after_predict) → process_predicted_model → phaser → refine
```

### 7. Fallback Program Tracking

**Rationale**: When LLM planning fails 3 times, the fallback mechanism must correctly report which program it selected.

**Implementation**:
- Fallback node now sets `state["program"]` when returning a command
- Response builder uses `state["program"]` if available, falls back to `intent["program"]`
- Prevents mismatch where PLAN shows one program but command is different

**Fallback diagnostics (v112.70):**
When the fallback cannot produce a command, it now provides specific stop reasons:

| Stop Reason | Meaning |
|------------|---------|
| `cannot_build_any_program` | No program could be built (missing required inputs) |
| `build_failures_and_duplicates` | Some programs can't build, rest are duplicates |
| `all_commands_duplicate` | All built commands are actual duplicates of prior runs |

The `abort_message` includes per-program diagnostics showing exactly which required
input slots could not be filled (via `CommandBuilder._last_missing_slots`).

### 8. Thinking Agent (v113)

**Rationale**: The planning LLM receives a structured prompt optimized for
program selection. It does not see raw program logs or have the context
to recognize crystallographic phenomena like twinning, anomalous signal
quality, or refinement convergence. A second "expert" LLM call fills
this gap by analyzing program output with domain expertise and providing
strategic guidance that PLAN can incorporate.

**Implementation**:

The THINK node sits between PERCEIVE and PLAN. It is controlled by
`thinking_level` (PHIL choice parameter: `none`/`basic`/`advanced`/`expert`,
default `expert`). At the default `expert` level, the THINK node runs
with full context (equivalent to `advanced` inside the graph).
Backward compatible: `use_thinking_agent=True` maps to `thinking_level=basic`.

```
PERCEIVE → THINK → PLAN → BUILD → VALIDATE → OUTPUT
               │
               ├─ thinking_level=none? → pass-through
               ├─ Not strategic program? → skip
               ├─ thinking_level=basic? → log analysis + LLM
               └─ thinking_level=advanced/expert? → validation + KB + metadata + LLM
                    (expert adds planning layer outside the graph, in ai_agent.py)
```

**When THINK engages** (`should_think()` in `agent/thinking_agent.py`):

| Condition | Rationale |
|-----------|-----------|
| After xtriage, phaser, autosol, autobuild | Key scientific decision points |
| After refine, ligandfit | Convergence and ligand assessment |
| After real_space_refine, map_to_model, predict_and_build | Cryo-EM strategic programs |
| After any failure | May need strategy change |
| R-free stalled 3+ cycles | Convergence problem |

THINK does NOT engage on the first cycle (no program output yet),
after non-strategic programs (e.g., validation-only or utility programs),
or when `thinking_level=none`.

**Four-tier thinking_level** (v113.10, v114):

The `thinking_level` parameter controls how much intelligence the
agent applies per cycle. Each level is additive. A separate parameter,
`use_rules_only=True`, is orthogonal — it replaces the LLM in PLAN
with deterministic selection but does not affect THINK, BUILD, VALIDATE,
or any safety checks.

**Level: `none`**

The THINK node is a complete pass-through (`run_think_node` returns
immediately when `thinking_level` is falsy). The pipeline is effectively
PERCEIVE → PLAN → BUILD → VALIDATE → OUTPUT. PLAN still calls the
LLM for program selection (one LLM call per cycle), unless
`use_rules_only=True` is set, in which case there are zero LLM calls.

The between-cycle loop in `ai_agent.py` skips all expert-mode
operations: no plan generation, no gate evaluation, no hypothesis
testing, no cycle commentary, no structure report.

**Level: `basic`**

The THINK node activates with a shallow context. `_build_thinking_context`
assembles the "basic context" only:

- Log sections extracted from the last program's output (per-program
  keyword extraction via `log_section_extractor.py`)
- Current metrics from `log_analysis` (R-free, map CC, resolution, etc.)
- R-free trend collected across cycle history
- Brief history summary (last 5 cycles with inputs/outputs)
- Recent failures list (last 5 failures with program, error, command)
- Strategy memory (accumulated observations from prior cycles)

Basic mode does NOT run structural validation, does NOT query the expert
knowledge base, does NOT build file metadata, and does NOT create or
update the Structure Model or Validation History.

An LLM call is made via `build_thinking_prompt()` with this context,
and the response is parsed into an assessment (action, confidence,
analysis, guidance, concerns). If guidance is produced, it is prepended
to `user_advice` as `[Expert assessment] ...` so the PLAN node's LLM
sees it. Strategy memory is updated and persisted via `session_info`.

The `should_think()` gate applies at all non-`none` levels: THINK only
engages after strategic programs (categories: analysis, model_building,
refinement, ligand — loaded from `programs.yaml`), after failures, or
when R-free is stalled (3+ cycles without improvement via
`StrategyMemory.metrics_stalled()`). THINK skips on the first cycle
(no program output yet).

The between-cycle loop is unchanged from `none`.

Net effect: two LLM calls per cycle when THINK engages (think + plan),
with log-analysis-quality reasoning injected into the planning prompt.

**Level: `advanced`**

Everything from `basic`, plus four subsystems activate in
`_build_thinking_context` (gated on `thinking_level in ("advanced", "expert")`):

*Phase A — Structural validation:*
`_run_structural_validation()` runs headless `run_validation()` on the
current best model (from `best_files["model"]`). Produces Ramachandran
outlier counts, clashscore, rotamer analysis, bonds/angles RMSD, and
model contents (chains, ligands, waters, ions). The result is formatted
into a compact text report by `format_validation_report()` and included
in the THINK prompt context.

*Phase B — Structure Model + Validation History (v114):*
A `StructureModel` object is created or restored from state. It is
updated every cycle from ground-truth validation results (not LLM
reasoning) via `update_from_validation()`, `update_from_xtriage()`,
and `update_from_phaser()`. Tracks data characteristics (resolution,
space group, twinning, anomalous), model state, R-free trajectory with
annotations, strategy blacklist, and hypotheses. The `ValidationHistory`
stores per-cycle snapshots with `get_metric_series()` for trend analysis
and `get_phase_start_metrics()` for the monotonic progress gate. Both
persist across cycles and session resume via state serialization. The
Structure Model summary replaces the raw validation report in the prompt,
providing richer cross-cycle structural knowledge. Current problems are
extracted via `get_current_problems()`.

*Phase C — File metadata:*
`build_file_metadata()` creates a metadata entry for the validated model
(validation stats, metrics, program, cycle). Stored in state and included
in the prompt context so the LLM knows quality characteristics of
available files.

*Phase D — File inventory:*
`_build_thinking_context` populates a `file_inventory` string from
`workflow_state["categorized_files"]`, listing filenames grouped by
category (models, reflection data, sequences, maps, ligands). This
appears as an `=== AVAILABLE FILES ===` section in the THINK prompt,
giving the LLM visibility into what files are available and how they
are classified — particularly useful when sigma-A MTZ files or MR
solution PDBs need special handling.

*Phase D — Expert knowledge base:*
`_query_knowledge_base()` runs IDF-weighted tag matching against a
56-entry crystallographic knowledge base. Takes the current validation
results, metrics, workflow stage, program, resolution, experiment type,
R-free trend, and xtriage results, and returns matching rules as text.
These rules encode domain expertise like "if clashscore > 40 after
autobuild, consider real-space refinement" or "if twinning detected,
add twin law to refine strategy."

Additionally, a hypothesis prompt is built from the Structure Model
(if hypotheses exist) and included in the THINK prompt, so the LLM
can propose new hypotheses or comment on existing ones. Hypotheses
extracted from the LLM response are added to the Structure Model.

The between-cycle loop in `ai_agent.py` is unchanged from `none`/`basic`:
no plan generation, no gate evaluation, no reports. The Structure Model
exists and accumulates knowledge, but it is only used within the THINK
node to inform the expert assessment. It feeds the GUI's Expert
Assessment display panel.

Net effect: two LLM calls per cycle when THINK engages, but with
dramatically richer context. The PLAN LLM receives enriched advice
that reflects actual structural validation, domain rules, and
accumulated structural knowledge.

**Level: `expert`** (default)

Everything from `advanced`, plus the entire Strategic Planner layer
activates. Inside the graph, the THINK node runs identically to
`advanced` — the value `expert` is mapped to `advanced` by
`create_initial_state()` since all planning operations live outside
the graph (see "How `expert` maps to `advanced`" below).

The additions are all in the between-cycle loop in `ai_agent.py`:

*At session start:*
- `_initialize_plan()` runs. The `PlanGenerator` selects from 12
  pre-defined templates based on experiment type, available files,
  resolution, and anomalous atoms. Template selection is deterministic.
  The LLM only customizes parameters within template bounds.
- The plan is a `StructurePlan` containing `StageDef` phases, each
  with programs, success criteria, gate conditions, fallbacks, and
  skip conditions.

*Before each cycle's graph invocation:*
- `plan_to_directives()` translates the current plan stage into the
  directives format: `prefer_programs`, `after_program`,
  `program_settings`. These are merged into the same directives dict
  that user advice produces, so the reactive engine treats plan
  directives identically to user directives.
- `plan_has_pending_stages` is computed and passed into `session_info`,
  which the PLAN node uses to suppress auto-stop and override LLM
  STOP decisions.
- A Strategy Hash is computed; changes trigger `advice_changed` in
  the reactive agent.

*After each cycle:*
- Gate evaluation: `GateEvaluator.evaluate()` checks the plan's
  success criteria for the current stage. Purely deterministic (no
  LLM). Returns one of: continue, advance, retreat, fallback, skip,
  or stop. Success hysteresis (1.5% buffer) prevents oscillation.
  Five anti-oscillation safeguards protect retreat logic.
- On advance: plan moves to next stage, `check_plan_revision`
  adjusts downstream stages, `generate_stage_summary` produces
  a transition summary.
- On retreat: failed strategy is blacklisted in the Structure Model
  (never retried), plan resets to retreat target phase,
  `advice_changed` is set.
- Hypothesis evaluation: `evaluate_hypotheses()` manages the lifecycle
  (proposed → testing → pending → confirmed/refuted/abandoned).
  Single active hypothesis budget. Verification latency via
  `test_cycles_remaining`.
- Cycle commentary: `generate_cycle_commentary()` produces a
  template-based (no LLM) crystallographer-level comment.
- Model placement gate: detects when model fits data (CC > 0.3 or
  R-free < 0.50) and locks `model_is_placed` — suppresses
  destructive programs and fast-forwards plan past MR/phasing stages.

*At session finalization:*
- `generate_final_report()` / `generate_stopped_report()` — LLM call.
- `structure_report.html` with inline SVG trajectory chart.
- `structure_determination_report.txt`.
- `session_summary.json` with metrics, stage outcomes, hypotheses.

**How `expert` maps to `advanced` in the graph**

The LangGraph pipeline only knows three thinking levels: `none`, `basic`,
and `advanced`. When the user sets `thinking_level=expert`,
`create_initial_state()` maps it to `advanced` for graph execution. This
is deliberate: the graph nodes (PERCEIVE, THINK, PLAN, BUILD, VALIDATE)
behave identically at `advanced` and `expert`. Everything that
distinguishes `expert` from `advanced` lives in the outer loop in
`ai_agent.py`, gated on the original `thinking_level == "expert"` check
from the PHIL parameter.

**`use_rules_only` as an orthogonal axis**

Setting `use_rules_only=True` replaces the LLM call in the PLAN node
with `RulesSelector` (deterministic first-valid-program logic via
`_mock_plan`). This change is confined to PLAN — everything else runs
identically:

- PERCEIVE: file categorization, workflow state, metrics, sanity checks
  — all unchanged.
- THINK: still controlled by `thinking_level`, which defaults to
  `expert` (mapped to `advanced` in the graph). The THINK LLM call
  is independent of `use_rules_only`. To suppress all LLM calls, set
  both `use_rules_only=True` and `thinking_level=none`.
- BUILD: runs identically, but the `intent` from `_mock_plan` has empty
  `files` and `strategy` dicts, so the builder relies entirely on
  auto-selection rather than LLM file/strategy hints.
- VALIDATE: workflow validation, file checks, duplicate detection — all
  unchanged.
- Between-cycle operations: unchanged (none at `thinking_level=none`;
  full planning at `thinking_level=expert`).

**Unintended rules fallback: LLM-unavailable observability**

The `use_rules_only=True` mode above is *intentional* and clearly announced.
The agent also falls back to `_mock_plan`/`RulesSelector` *unintentionally*
when the requested LLM provider cannot be reached: `plan` calls
`_handle_llm_failure`, which (below `MAX_LLM_FAILURES`) returns `_mock_plan`.
That fallback used to be silent, so a run could complete on rules without the
user realizing the LLM never contributed. v120 makes it observable, with a
deliberate split that survives both local and server execution:

- **Signal carrier — the `events` list.** `_handle_llm_failure` emits a
  structured `EventType.NOTICE` event with a machine-readable discriminator
  `notice_kind="llm_unavailable"` (plus `provider` and a human `message`), via
  `_emit` (which copies the list — no in-place mutation). `events` is a plain
  `List[Dict]` in `AgentState` (no reducer), preserved by every downstream node
  (`build`/`validate`/`fallback`/`output` all return `{**state}` or never
  clobber it), and is a member of `contract.RESPONSE_FIELDS`, so it crosses the
  REST boundary identically — `run_ai_agent._build_response_from_state` passes
  `events` into both `create_response` and `create_stop_response`. No change to
  the contract, schema, server dispatcher, event log, or formatter was needed.
- **Per-cycle vs per-run scoping.** Graph state is per-cycle (each cycle is a
  fresh `graph.invoke`), so there is no run-level flag in graph state. Within a
  single invoke, emission is guarded once by `llm_notice_emitted_this_invoke`
  (the `validate→plan` retry edge can call the handler more than once per
  cycle) — a key deliberately scoped to the invoke and **not** named like the
  driver's run-level key, to avoid scope confusion. The NOTICE is therefore
  emitted each cycle the LLM is unavailable (honest per-cycle reporting).
- **Driver aggregation.** `ai_agent.print_history_record` →
  `_detect_llm_unavailable_from_events` scans each cycle's
  `history_record["events"]` for the discriminator and sets the run-level
  `session.data["llm_ever_unavailable"]` (which persists across cycles). The
  hook is guarded `if session is not None` (the display-only/replay path calls
  `print_history_record` without a session) and never raises.
  `_finalize_session` prints an `IMPORTANT: THIS RUN DID NOT USE THE LLM`
  banner, prefixes the Results-tab summary, and exposes `llm_ever_unavailable`
  on the result object for the GUI.
- **Visibility is double-covered.** At normal/verbose verbosity the NOTICE
  renders inside the cycle box via `EventFormatter`; at quiet verbosity (where
  the formatter's quiet path omits NOTICE) the same text still appears because
  `_handle_llm_failure` also logs it to `debug_log`. So the per-cycle notice is
  visible at every verbosity through at least one channel, and driver detection
  (from the structured event) is verbosity-independent.

This is observability only — it does not change program selection. Known gap
(documented for a possible follow-up): a pure display-only replay of a saved
session that never re-runs cycles will not show the end-of-run banner, because
session cycle records do not persist `events`; closing that would require a
`session.py` schema change.

**Log section extraction** (`agent/log_section_extractor.py`):

PHENIX logs can be megabytes. The extractor uses per-program keyword tables
to pull out scientifically informative sections within a 3500-character budget.
Sections are priority-ordered so the most important information is always included
first:

| Program | Priority sections (first = highest) |
|---------|--------------------------------------|
| phenix.xtriage | Twinning → Anomalous signal → Wilson statistics → Space group |
| phenix.phaser | MR scores (TFZ/LLG) → Packing → Search strategy |
| phenix.autosol | Phasing stats (FOM/BAYES-CC) → Density modification |
| phenix.autobuild | Building progress → Map quality |
| phenix.refine | R-factors → Geometry → Twinning → Difference map |

For each keyword match, a window of 5 lines before and 15 lines after is
captured and merged with overlapping windows. Unknown programs fall back to
the last 100 lines.

**LLM response format**: The expert LLM returns structured JSON:

```json
{
  "analysis": "Strong anomalous signal to 3.0 Å ...",
  "confidence": "high|medium|low|hopeless",
  "action": "guide_step|let_run|stop|pivot",
  "guidance": "Run autosol with atom_type=Se ...",
  "data_quality": "good",
  "phasing_strategy": "SAD",
  "concerns": ["weak anomalous"],
  "alternatives": ["MR with AlphaFold as backup"]
}
```

**How guidance reaches PLAN**: The assessment's `guidance` string is prepended
to `user_advice` as `[Expert assessment] ...`. This injects domain reasoning
through the existing advice channel — PLAN and BUILD both see it naturally
without any new interface. The expert presents evidence and reasoning;
PLAN makes the final decision.

**Strategy memory** (`agent/strategy_memory.py`):

Persists across cycles via `session_info`, tracking data quality observations,
phasing strategy, R-free history, concerns, and decisions. Used by
`should_think()` to detect R-free stalling and by the prompt builder to
give the expert context about prior assessments. Capped fields prevent
unbounded growth (10 concerns, 20 R-free values, 30 programs, 10 decisions).

**Graceful degradation**: If the thinking LLM fails (rate limit, timeout,
bad parse), the THINK node logs the error and returns state unchanged. PLAN
runs normally as if THINK were disabled. The workflow never crashes due to
a thinking failure.

**Files**:

| File | Purpose |
|------|---------|
| `agent/thinking_agent.py` | Core module: `should_think()`, `run_think_node()`, context assembly, LLM call, validation/KB/metadata integration |
| `agent/strategy_memory.py` | `StrategyMemory` class with serialization, `metrics_stalled()` |
| `agent/log_section_extractor.py` | `extract_sections()` with per-program keyword tables |
| `agent/validation_inspector.py` | Headless structural validation (v113.10) |
| `agent/format_validation.py` | Compact text formatter for validation results (v113.10) |
| `agent/file_metadata.py` | Per-file metadata tracking and queries (v113.10) |
| `agent/kb_tags.py` | Context tag derivation for KB lookup (v113.10) |
| `agent/phil_validator.py` | PHIL strategy validation against programs.yaml strategy_flags (v115) |
| `agent/error_classifier.py` | Error classification + tiered failure response (v115) |
| `agent/intent_classifier.py` | 4-way intent classification: solve/solve_constrained/task/tutorial (v115) |
| `agent/sanity_checker.py` | Pre-execution sanity checks: experiment type, model existence, PHIB guard (v115) |
| `knowledge/thinking_prompts.py` | Prompt builder: `build_thinking_prompt()`, `parse_assessment()` |
| `knowledge/kb_loader.py` | YAML KB loader with IDF-weighted tag matching (v113.10) |
| `knowledge/expert_knowledge_base_v2.yaml` | 56 expert-reviewed entries (v113.10) |
| `knowledge/tutorial_expectations.yaml` | Per-tutorial expected behavior for dual-run evaluation (v115) |
| `knowledge/solve_run_info.yaml` | Minimal hints for solve-mode runs (v115) |

**Data flow** (see per-level detail above):

```
GUI/CLI → params.ai_analysis.thinking_level
  → LocalAgent/RemoteAgent →
      request["settings"]["thinking_level"]
  → run_ai_agent.py →
      create_initial_state(thinking_level=...)
      ("expert" maps to "advanced" for graph)
  → graph: perceive → think → plan → build →
      validate → output
  → response: strategy_memory, expert_assessment,
      structure_model, validation_history
  → ai_agent.py → session.data (persisted)
  → event_formatter: Expert Assessment block

When thinking_level=expert (additional, in ai_agent.py):
  Session start: _initialize_plan()
  Before graph:  plan_to_directives()
  After graph:   GateEvaluator.evaluate()
                 evaluate_hypotheses()
                 generate_cycle_commentary()
  Transitions:   generate_stage_summary()
  Session end:   generate_final/stopped_report()
  Output:        session_summary.json +
                 structure_report.html +
                 structure_determination_report.txt
```

**Advanced mode prompt example** (v113.10):

```
CYCLE 5 | xray | workflow: refinement
Last program: phenix.refine
Current metrics: R-free=0.248, R-work=0.219

=== VALIDATION (cycle 5: phenix.refine) ===
R-work=0.219 R-free=0.248 (prev: 0.295, start: 0.421)
Bonds=0.0060 Angles=0.82
Rama: 97.2% fav, 0.0% outlier  Clashscore: 3.2
Ligands: ATP (A/301)  Waters: 187
R-free trend: 0.421 -> 0.295 -> 0.248 (3 cycles)

=== FILE METADATA ===
Best model: model_005.pdb (A,B; 1x ATP; 187 waters; R-free=0.248)

=== RELEVANT EXPERT RULES ===
[stop_001] R-free plateau after 3 cycles...

Strategy memory: [cycles 1-3: R-free 0.42 -> 0.30]
--- PROGRAM OUTPUT ---
[extracted log sections...]
```

Context budget: ~2550 chars total (~350 validation, ~50 trend, ~1200 KB
rules, ~150 metadata, ~600 memory, ~200 context), leaving ~950 headroom
under the 3500-char limit.

**Follow-up items** (v113.10):

| ID | What | Priority |
|----|------|----------|
| F4 | Real-structure test suite (10-20 structures) | 1 — Highest |
| F1 | Ligand RSCC computation (stub in place) | 2 |
| F2 | Difference density peak search (stub in place) | 3 |
| F3 | Filename → metadata file selection migration | 4 |
| F5 | KB tag derivation tuning after test suite | 5 |
| F6 | Persist validation reports in session state | Optional |
| F7 | Add validation_result/file_metadata to TypedDict | Optional |

**Risk register** (v113.10):

| Risk | Mitigation |
|------|------------|
| Validation is slow | Runs only in advanced mode |
| Validation crashes | Never-raises pattern throughout |
| KB entries are wrong | 56 focused entries, 45 high-confidence |
| Context window too full | ~2550 chars, 950 headroom |
| Advanced mode untested on real structures | Needs F4 test suite |

### 9. Stop-after Routing (v117 stop_refactor)

**Rationale**: The `after_program` directive carries two distinct
kinds of intent that earlier code treated identically:

1. **User-explicit stop** — the user wrote "X and stop" / "stop after X"
   and means it.  The workflow should terminate after X completes.
2. **Plan-progression hint** — the planner injected `after_program`
   per-stage as a min-run target.  The workflow should advance to the
   next stage after X completes, NOT terminate.

Pre-v117, both were treated as "terminate after X" and v112.78 wiped
`valid_programs` to `["STOP"]` whenever `after_program_done`.  This
caused regressions like AF_7mjs where the planner's stage hint was
interpreted as a user-explicit stop.

v117 introduces `stop_conditions.stop_after_requested` as a boolean
flag that distinguishes the two cases.  The flag is set by:

- The directive extractor's `_is_stop_after_requested` helper —
  detects 11 positive patterns including "X and stop", "stop after X",
  "only run X", "just run X", "Stop Condition: <real value>", and the
  period-stop pattern "Refine. Stop." added in v117 Step B.
- The LLM directly — when the JSON schema sets it (per the
  `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` schema documentation).

When the flag is True AND `after_program_done` is True,
`_apply_directives` wipes `valid_programs` to `["STOP"]`.  When the
flag is False/absent, `_apply_directives` does nothing — the plan
progresses naturally.

The prompt-formatting layer (`_format_directives_for_prompt` in
`knowledge/prompts_hybrid.py`) follows the same gate: when the flag
is True, emit the Phase 6a "Stop target: X" 4-line block; when False,
emit nothing about `after_program`.

**Why important for future work**: anyone modifying the after_program
handling must preserve this two-case distinction.  Re-introducing the
v112.78 "wipe everywhere" pattern would resurrect AF_7mjs.  The unit
tests in `tst_validate_step_after_program_guard.py` (M1-M7) and
`tst_after_program_prompt.py` (N1-N5) lock both halves of the
behavior.

### 10. Raw-Advice-Authoritative Directive Extraction (v117 Step 1)

**Rationale**: The advice preprocessor LLM occasionally mangles short
imperative input — most reproducibly, openai expanding
`"density modify and stop"` into multi-action prose with
`"Stop Condition: None"`.  The downstream directive extractor saw only
the mangled output and missed the user's stop intent.

v117 Step 1 gives the directive extractor both inputs:

- `user_advice` — preprocessed text (file lists, parameter
  normalization, experiment type inference)
- `raw_advice` — original user instruction, when different

A new prompt variant `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` includes
an AUTHORITY paragraph telling the LLM that raw is the source of
truth for intent fields (`after_program`, `start_with_program`,
`stop_after_requested`) and processed fills gaps for files,
parameters, and experiment type.

The dual-input path is activated by `extract_directives()` when
`raw_advice` is supplied and differs from `user_advice`.  When raw and
processed match (or raw is None), the original single-input prompt
runs unchanged — backward compatible.

Plumbing: client sets `session.data["raw_advice"]` →
`programs/ai_agent.py::_extract_directives` reads it →
`programs/ai_analysis.py` PHIL `user_advice_raw` plumbs it across the
client-server boundary → `phenix_ai/run_ai_analysis.py::run_directive_extraction`
passes it through to `extract_directives()`.

**Why important for future work**: any caller that constructs
directive-extraction requests should pass both `user_advice` and
`raw_advice`.  Tests verifying the dual-input path are in
`tst_extract_raw_advice.py`.

### 11. Grounding ∧ stop_after_requested Interaction (v117.1)

**Rationale**: The v116.19a grounding guardrail
(`_validate_after_program_grounded`) drops `after_program` when the
program name doesn't appear literally in the user advice, treating it
as a "fabrication" from descriptive prose.  This was correct in
v116.x but conflicted with v117 Step 1: when the AUTHORITY paragraph
correctly led the LLM from "density modify" (in raw) to
`phenix.resolve_cryo_em` (a learned mapping, not a literal substring),
the guardrail dropped the program.

The C1 LLM test (`tst_c1_raw_authority.py`) surfaced this
interaction: both providers produced `stop_after_requested=True` plus
the LLM-mapped program, and the guardrail uniformly dropped the
program.

v117.1 resolves the conflict: when the LLM sets BOTH `after_program`
AND `stop_after_requested=True` on the same `stop_conditions` block,
the flag IS the grounding signal — the LLM made an explicit user-stop
assertion based on raw advice and the AUTHORITY paragraph, which is
strictly stronger evidence than the heuristic literal-substring check.
The guardrail's check is skipped in that one case.

AF_7mjs is preserved by construction: AF_7mjs's preprocessed advice
has no positive stop signal, so `stop_after_requested` is not set,
and the guardrail fires as before.

**Why important for future work**: when modifying either the
grounding guardrail or the `stop_after_requested` flag handling,
preserve this interaction.  The boundary tests in
`tst_grounding_stop_after_requested.py` (K1-K5) lock the behavior:
K1 (both signals → kept), K2 (flag absent → dropped),
K3 (pure fabrication → dropped), K4 (flag False → dropped),
K5 (grounded program no flag → kept).

### 12. Filling after_program from Raw Advice (v117.2)

**Rationale**: After v117.1 fixed the grounding guardrail, a C1 LLM
test surfaced a different gap: openai's extractor occasionally emits
`stop_after_requested=True` but omits `after_program` for inputs like
`"refine and stop"` (raw) with mangled processed advice
(`"Stop Condition: None"`).

Without `after_program`, `workflow_engine._apply_directives`' gated
wipe path is unreachable — the wipe is gated on `if after_program:`.
The user's stop intent is signaled but cannot be acted on.

v117.2 closes this gap: in `extract_directives()`, after the existing
LLM-extraction post-processing, when `stop_after_requested=True` and
`after_program` is missing, call `_resolve_after_program()` against
the raw advice.  The resolver uses `_ACTION_TABLE` to map action
verbs (`refine`, `xtriage`, `density modify`, ...) to programs.

The fix has three guards:
1. Only fires when the LLM left `after_program` unset — preserves
   LLM choices.
2. Only fires when raw advice contains stop intent per
   `_is_stop_after_requested` — guards against amplifying a possibly-
   hallucinated flag.
3. Uses `raw_advice` if available; falls back to `user_advice` for
   callers using the single-input path.

**Why important for future work**: the accumulated post-LLM processing
in `extract_directives()` now has a defined order (grounding →
symmetry → fill_in_from_raw → ollama-empty → intent merge → intent
override → workflow intent fallback).  The v117.2 fill-in is
positioned after grounding (so it can't be undone) and before intent
overrides (so they don't accidentally clear what the resolver just
confirmed).  Boundary tests in `tst_after_program_fill_from_raw.py`
(K1-K6) lock the interactions.

### 13. Extended Stop-Intent Phrasing Recognition (v117.3)

**Rationale**: After v117.2 fixed the C1 "missing after_program"
gap, a separate failure remained in the existing
`tst_directive_extraction.py::explicit_stop_after_phaser` test:
openai's extractor consistently produced `after_program=phenix.phaser`
but failed to set `stop_after_requested=True` when the user's
phrasing was `"Stop the workflow immediately after phaser completes"`.
None of the v117.2 regex backstops or imperative markers recognized
this phrasing, so v117.1's flag-exemption couldn't fire and the
grounding guardrail's Failure 2 dropped the directive.

v117.3 closes the gap with three small additive changes to data
structures `extract_directives` already consumes:

- **`_IMPERATIVE_STOP_MARKERS`** gains 5 new entries: `"stop the
  workflow"`, `"immediately after"`, `"after it completes"`,
  `"after it finishes"`, `"is the last step"`.  These are
  substring-matched within a 300-char window around the program
  name, used by the grounding guardrail's Failure 2 path.

- **`_POSITIVE_STOP_AFTER_PATTERNS`** gains 2 new contiguous-phrase
  regex patterns: `\bstop\s+the\s+workflow\b` and
  `\bis\s+the\s+last\s+step\b`.  These feed
  `_is_stop_after_requested`, which gates v117.1's grounding-flag
  exemption and v117.2's fill-in-from-raw.

- **`DIRECTIVE_EXTRACTION_PROMPT`** schema docs for
  `stop_after_requested` gain 4 new recognized phrasings and 2
  program-neutral few-shot examples.  The dual-input prompt
  inherits via `.replace()`.

**Why important for future work**: v117.3 demonstrates a pattern
that may recur as new LLM-extraction phrasings surface: bare addition
to existing data structures (the marker tuple, the pattern tuple, the
prompt docstring) without introducing new pipeline steps or new state
fields.  Boundary tests in `tst_extended_stop_phrasings.py` (K1, K2,
K4, K5, K6, K7a/b, K8a/b) lock the behaviors.

**Cross-sentence safety**: both new regex patterns are contiguous
multi-word phrases ("stop the workflow", "is the last step"), making
cross-sentence false matches structurally impossible.  Over-permissive
"after X completes/finishes" patterns were considered and rejected
during plan review after the K5 false-positive verification.

**Why the imperative markers added entries that DON'T also appear in
the regex patterns list**:  the markers list is window-bounded (300
chars near program name), so `"immediately after"` is safer there
than as a global regex.  Adding `"immediately after"` as a global
regex pattern would fire on descriptive prose like "the workflow
continues immediately after prediction."  The asymmetric placement
is deliberate.

### 14. Diagnostic Markers for LLM Pipeline Failures (v118 Section E)

**Rationale**: Pre-v118, when an LLM call in the preprocessor or
directive extractor failed (network error, JSON parse failure,
retired model 404, malformed response), the failure was caught
inside `_call_llm` and the code fell through to the rules-only
backstop.  This was correct behavior for resilience but produced
no diagnostic trail: from outside the function, a successful
rules-only fallback looked identical to a successful LLM call.

v118 Section E adds two stderr markers:

- **`[DIRECTIVE_EXTRACTION_FAILED]`** at
  `agent/directive_extractor.py` — fires when the directive-
  extraction LLM call raises any exception.  Message includes
  the provider, model name, and full exception text.

- **`[ADVICE_PREPROCESSING_FAILED]`** at
  `agent/advice_preprocessor.py` — analogous marker for
  preprocessor LLM failures.

Both markers dual-write to the log function AND stderr (the
stderr write is wrapped in `try/except` so logging can never
break the fallback path).

**Why important for future work**: Section E proved essential
for diagnosing v118 Section 8 (the gemini-2.0-flash retirement).
Without Section E, the server-side 404 would have been buried in
server stderr with no marker the client could grep.  With
Section E, the marker name itself was the search key, and the
root cause was visible from client logs alone.  v119's planned
"server stderr → client diagnostic_messages" feature is a direct
continuation of this principle.

### 15. Diagnostic Layer Split: User Intent vs Effective Runtime (v118 Section C-prime)

**Rationale**: Pre-v118, the displayed "Extracted Directives:"
block in client logs was overloaded: it showed the LLM's
extracted directives but ALSO the post-grounding and
post-fallback overlays applied by `_apply_directives`,
`_validate_after_program_grounded`, and friends.  When provider
asymmetry caused different behavior (e.g., OpenAI vs Google on
the same advice), the client log alone couldn't tell whether
the LLM extracted differently or whether the runtime overlays
differed.

v118 Section C-prime splits the diagnostic output into:

- **`directives_user_intent`** — the LLM's pure extraction
  output, before any runtime overlays
- **`directives_effective_runtime`** — what the workflow
  actually consumed after all overlays applied

The two are displayed in adjacent sections of the cycle log.
When they match, the LLM-vs-rules path is straightforward.
When they differ, the difference is exactly the runtime
overlay's contribution, which the user can compare against
documented overlay rules to diagnose.

**Why important for future work**: this is the pattern any
future "what did the LLM actually emit" question relies on.
Section 9's `_corrected_from` sidecar (see Section 21 below)
follows the same principle of preserving provenance through
runtime transformations.

### 16. PHIL Namespace Healing (v118 Section B)

**Rationale**: LLMs sometimes emit PHIL paths with wrong
namespace prefixes — e.g., `data_manager.r_free_flags.generate`
where the canonical PHIL path is bare
`generate_rfree_flags=True` for phenix.refine, or
`refinement.refine.strategy` where the canonical is
`refine.strategy`.  Pre-v118, these wrong-namespace paths
silently flowed to phenix programs and were rejected by PHIL
parsing, producing "unknown parameter" errors with no
attribution to the LLM as the source.

v118 Section B adds:

- **`PHIL_NAMESPACE_TRANSLATIONS`** static table mapping known-
  bad PHIL paths to their canonical forms
- **`_heal_namespaced_phil_keys()`** helper applied in
  `validate_directives` per-program — for each setting key
  the LLM emitted, check the translation table and rewrite
  if matched
- **`BAD_PHIL_NAMESPACE_PREFIXES`** companion table — drop
  (with log line) any key whose prefix matches a known-bad
  pattern that has no canonical translation

Section B's `tst_phil_namespace_cleaner.py` (14 tests) locks
the translation table.

**Why important for future work**: the translation table is
finite and grows with each observed LLM mistake.  v119's
deferred work item §3.1 considers replacing the static table
with dynamic `master_phil` resolution at extract time, which
would eliminate drift between the table and actual program
schemas.

### 17. BUILD experiment_type Threading + R-free Auto-fill (v118 Section F)

**Rationale**: `agent/command_builder.py:_select_files` had
two pre-v118 brittlenesses:

1. **Hardcoded "xray" assumption in the cycle-1 fast path**:
   when `experiment_type` wasn't yet locked (because xtriage
   hadn't run on cycle 1), the file selector defaulted to
   X-ray file slots.  For cryo-EM tutorials where the planner
   skipped xtriage and jumped directly to phenix.refine, this
   produced wrong file selections.

2. **No R-free generate-flag injection**: when the user
   supplied an MTZ with no R-free flags and asked for
   refinement, phenix.refine would error.  The right command
   line is
   `phenix.refine ... xray_data.r_free_flags.generate=True`,
   but Section F's predecessors had no logic to detect
   "first refinement + no rfree_mtz locked" and add the flag.

v118 Section F adds:

- **`experiment_type` threading** through `_select_files`:
  the parameter is passed explicitly from cycle-1 plan
  context (where it's known from input file extensions even
  if not yet locked in session state).  Replaces the hardcoded
  default with a true value.
- **Auto-fill of `xray_data.r_free_flags.generate=True`** when
  `_select_files` detects refinement-without-rfree-mtz on a
  first-refinement cycle.  BUILD log line:
  `BUILD: First refinement - will generate R-free flags (no rfree_mtz locked)`.
  (v120 hardens this: generation is now also suppressed when the
  selected input data MTZ already carries an R-free array — see
  "R-free generate guard" under the recovery-errors section.)

Section F is the **first v118 section that touches server-side
code**.  After Section F, the v118 ledger formally distinguishes
client-side from server-side changes; server-side changes
require both file deployment AND server process restart.

**Why important for future work**: Section F's "experiment
type isn't known yet" handling is a workaround for the deeper
issue that `session.set_experiment_type()` only locks after the
first program returns.  v119 deferred work item §3.3 considers
inferring experiment type from input file extensions at session
creation time, which would obviate Section F's threading hack.

### 18. Optional Dependency Resilience + Environment Probe (v118 Section G + 6.7)

**Rationale**: The `rag/vector_store.py` module conditionally
imports `chromadb` via `langchain_chroma`.  Pre-v118, this
import was wrapped in `except ImportError`.  Production
deployments started failing with `TypeError` from the
opentelemetry-proto / protobuf chain (a dependency of chromadb
that fails on version mismatches with a TypeError, not an
ImportError).  Result: the LLM-features-disabled fallback path
didn't fire, and tutorial runs failed at module load time.

v118 Section G:

- Broadens the exception class from `except ImportError` to
  `except Exception` at the chromadb import site
- Logs the specific exception text so the operator can see
  which dependency version is conflicting
- Falls through to LLM-features-disabled path on any chromadb
  load failure

v118 Section 6.7 (rev 3) adds **`tests/tst_dependencies.py`** —
an environment-readiness test that runs a **true runtime probe**
(`from langchain_chroma import Chroma`) rather than just
attempting to import top-level packages.  This catches the
exact deep-import failure modes Section G was designed for.

Together, Section G and 6.7 are the first v118 sections verified
on both Mac and Linux production environments.  v118.G also
established the **stub-module isolation pattern** for K-tests
(see `tests/tst_optional_dep_resilience.py`): the test
synthesizes a stub `chromadb` module that raises the protobuf
TypeError on import, exercising Section G's broadened exception
handler without needing a real broken chromadb install.

**Why important for future work**: the "catch Exception not
just ImportError" lesson generalizes to any lazy-import gate
in the codebase.  v119 deferred work item §3.4 calls for an
audit of other lazy-import sites.

### 19. Default Model Names + Retired-Model Detection (v118 Section 8)

**Rationale**: Google retired `gemini-2.0-flash` on
2026-05-20.  Production server's directive extraction began
silently 404'ing because three independent hardcoded
model-default tables in the codebase still pointed at the
retired name:

1. `core/llm.py` (already pointed at current model)
2. `agent/api_client.py` (at TWO sites: lines 655 and 673)
3. `agent/directive_extractor.py` (fallback path)

v118 Section 8 rev 3 bumps all three to
`gemini-2.5-flash-lite`.  Section 8 rev 1 fixed only the
directive_extractor — the api_client occurrences were missed,
and the bug persisted.  Rev 2 needed a patch script; rev 3
shipped the full files.  This experience generated the v119
cadence convention: every string-replacement starts with
`grep -rn 'STRING' $PHENIX/modules/cctbx_project/libtbx/langchain/`.

Section 8 is **client + SERVER**: the server's long-running
process has the directive_extractor and api_client modules
loaded at boot.  Code changes don't take effect until
the server restarts.  This experience also informs v119's
planned operational guards (server `/version` endpoint +
startup canary that issues a 1-token dummy inference call to
catch retired-model 404 at boot).

**Why important for future work**: v119 deferred work item
§3.5 calls for centralizing all model defaults in
`core/llm.DEFAULT_MODELS` with a single
`default_model_for_provider()` accessor.  Then the next
retirement requires one string change rather than three.

### 20. Experiment-Type-Conditional Program Canonicalization (v118 Section 9)

**Rationale**: Tom's runs 237 and 239 ("density modify and
stop" on cryo-EM half-maps) produced
`after_program=phenix.autobuild_denmod` (the X-ray density
modification program) instead of `phenix.resolve_cryo_em`
(the cryo-EM equivalent).  Because `autobuild_denmod` never
runs in a cryo-EM workflow, the stop guard never fired and
the workflow continued through all 5 stages.

The bug was at the LLM-extraction layer: the directive-
extraction prompt at lines 443-445 had two parallel rules
(cryo-EM density mod, X-ray density mod) both matching the
literal phrase "density modify".  The LLM picked the X-ray
default (likely recency bias plus the X-ray rule appearing
second) even when the preprocessed advice clearly stated
`Experiment Type: cryo-EM`.

v118 Section 9 adds three layers:

1. **Prompt restructure** (lines 443-462): replaces the two
   parallel rules with an explicit decision tree —
   "CHECK THE EXPERIMENT TYPE in the preprocessed advice
   FIRST" — and a default-to-cryo-EM rule for ambiguous cases.
2. **Module-level `PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE` table**
   with two symmetric entries:
   `("cryoem", "phenix.autobuild_denmod"): "phenix.resolve_cryo_em"`
   and the X-ray mirror.  Read by
   `_apply_experiment_type_program_reprints()` which runs
   AFTER `_validate_after_program_grounded`.  Placement is
   intentional: grounding bypasses when
   `stop_after_requested=True`, so the wrong program survives
   intact for correction.  Moving the validator before
   grounding would create the opposite bug.
3. **`[DIRECTIVE_CORRECTION]` log marker** + `_corrected_from`
   sidecar field in `stop_conditions` recording
   `{from, to, reason, experiment_type}` for traceability.

K_DENMOD test suite (13 tests) covers the bug case, mirror
case, no-change cases, ambiguous-experiment cases, grounding-
bypass interaction, and diagnostic emission verification.

**Why important for future work**: Section 9 demonstrates the
"table-driven canonicalization" pattern: a declarative
mapping table consumed by a generic validator, with new entries
addable in one line.  v119 may extend the table conservatively
if production logs reveal more analogous LLM mistakes; pairs
that do semantically-different things (xtriage vs mtriage)
should NOT be added.

### 21. List-to-String Coercion in Directive Validation (v118 Section 10)

**Rationale**: Tom's P9 SAD run with Google provider produced
`additional_atom_types=['S']` (the Python repr of a single-
element list) instead of the canonical `S`.  autosol rejected
the bracketed string.  The same input on OpenAI produced `S`
correctly.

Root cause: provider-specific JSON shape mismatch.

```
OpenAI:  "additional_atom_types": "S"        ← string
Google:  "additional_atom_types": ["S"]       ← list of strings
```

`validate_directives` called `expected_type(value)` with
`expected_type=str`:

```python
str("S")     →  "S"        # OpenAI: correct
str(["S"])   →  "['S']"    # Google: Python repr - WRONG
```

v118 Section 10 adds **`_coerce_setting_value(value, expected_type)`**
helper applied at TWO call sites in `validate_directives`:

- **Site 1 (program_settings, ~line 2104)** — Tom's bug.
- **Site 2 (stop_conditions, ~line 2209)** — latent
  silent-drop: the pre-coerce check
  `value not in VALID_PROGRAMS` raises
  `TypeError("unhashable type: 'list'")` on list input,
  swallowed by the outer `except (ValueError, TypeError)`.
  Today Google emitting `after_program=["phenix.autosol"]`
  silently drops the directive.

Helper behavior:

- **Single-element list/tuple of any type**: unpack and
  recurse.  Type-preserving for bool/int/float/str.  Fixes the
  Python truthiness trap where `bool([False]) == True` would
  silently flip a boolean.
- **Multi-element list/tuple + str expected_type**: space-join
  (PHIL multi-value-string syntax).
- **Multi-element + non-str**: raise ValueError → outer
  except → logged and dropped.
- **Dict + str**: raise ValueError → same.
- **All other cases**: `expected_type(value)` as before
  (backward-compat for OpenAI's scalar path).

Emits one of two log markers:

- `DIRECTIVES: Unpacked single-element list for X: Y → Z`
  (the JSON-shape-artifact case)
- `DIRECTIVES: Joined multi-element list for X: Y → Z`
  (the content-ambiguity case)

K_LIST test suite (12 tests) including K11 which verifies the
**v118.10 + v118.9 interaction chain**: a list-formatted
`after_program=["phenix.autobuild_denmod"]` in cryo-EM context
is coerced by v118.10 to a string, which v118.9 then rewrites
to `phenix.resolve_cryo_em`.  Without v118.10 the list would
have been silently dropped before v118.9 could correct it.

**Why important for future work**: Section 10's helper
documents a **future-proofing convention** in its docstring —
if a future PHIL field legitimately needs a list value,
declare it with `expected_type=list` in `VALID_SETTINGS` and
the helper passes lists through unchanged.  This structurally
prevents future maintainers from breaking the helper by
adding a list-bearing field with the wrong type annotation.

A latent analogous bug remains in
`file_preferences` validation (line ~2371 of
`directive_extractor.py`): `bool(value)` on
`prefer_anomalous=[False]` would silently return True.
v118.10 doesn't touch that block because Tom's bug is in
`program_settings`.  v119 deferred work item §3.6 calls for
extending the helper to that block.

### 22. Centralized LLM Model Defaults (v119.H1)

**Rationale**: Before v119, model names were repeated as default
arguments at each call site: `model="gpt-4o-mini"` in three
places of `agent/api_client.py`, `model="gemini-2.5-flash-lite"`
in `agent/directive_extractor.py`, and so on.  Each retirement
update required hunting through the codebase for every
occurrence.  Worse, the **default model for a given
(provider, role)** was not a queryable property of the system —
it was implicit in scattered keyword arguments.

v119.H1 establishes five DEFAULTS tables in `core/llm.py`:

- `DECISION_MODEL_DEFAULTS` — what extractor/planner/intent
  classifier use.
- `RAG_MODEL_DEFAULTS` — what the RAG retrieval LLM uses.
- `RAG_EMBEDDING_DEFAULTS` — embedding model names.
- `EXPENSIVE_MODEL_DEFAULTS` — heavyweight reasoning calls.
- `CHEAP_MODEL_DEFAULTS` — lightweight classification calls.

Plus a `default_model_for_provider(provider, role="DECISION")`
helper that all call sites use.  `agent/api_client.py` and
`agent/directive_extractor.py` were refactored to import and
use this helper.

A `RETIRED_MODELS` frozenset names models known to be retired
(e.g., `"gemini-1.0-pro"`, `"gpt-3.5-turbo-16k"`).  When a
request comes in with an explicit model name that's in
`RETIRED_MODELS`, the call site emits a stderr marker
`[DIRECTIVE_EXTRACTION_MODEL_RETIRED]` before falling through to
the default.  This is defense-in-depth against the v118.8 class
of incident (server returning a 404 because the configured model
was retired by the provider).

**Why important for future work**: with one place to update,
retirement-driven model changes become a single-line edit.  The
H2 fingerprint mechanism (KDD 23) intentionally excludes
`RETIRED_MODELS` from its hash so retirement updates don't
trigger canary fingerprint drift.  Future work that needs to
inspect "what model would be used for X" can call
`default_model_for_provider("google", "DECISION")` rather than
reading source.

K_H1 (`tests/tst_default_models.py`, 22 tests) locks the
invariants: every (provider, role) combination resolves, all
DEFAULTS keys are stable strings, RETIRED_MODELS membership
triggers the stderr marker.

### 23. Server Build Metadata Channel — agent_build (v119.H2)

**Rationale**: Before v119, the client had no programmatic way
to know which server build it was talking to.  A wrong-version
server (built from the wrong branch, or with a stale
`core/llm.py`) would behave subtly differently but be
indistinguishable from a correct deploy in client logs.  This
class of incident is hard to catch because individual responses
look fine — only an aggregate "is this the build I deployed?"
question would surface the issue.

v119.H2 establishes an unconditional server→client metadata
channel.  Every server response carries an `agent_build` object:

```json
{
  "agent_build": {
    "version": "119.H3",
    "defaults_fingerprint": "sha256:77bf7421...",
    "started_at": "2026-05-22T18:42:13Z"
  }
}
```

- **`version`** is the contents of `VERSION` (new file at
  `langchain/` root), accessible via `core/_version.py`.  When
  Tom edits VERSION and redeploys, the next response reflects it.
- **`defaults_fingerprint`** is `sha256:<64-hex>` covering all
  five DEFAULTS tables from KDD 22, EXCLUDING `RETIRED_MODELS`.
  Identical fingerprints across responses mean identical tables
  across processes; drift means someone edited `core/llm.py`
  without bumping VERSION.
- **`started_at`** is strict UTC ISO 8601 captured once at module
  load.  Lets the operator tell "is this a fresh process?" from
  the response alone.

Injection happens at the top of `_build_group_args_response` in
`phenix_ai/run_ai_agent.py`, before serialization.  Both success
and error responses get the channel.  The H3a canary explicitly
exercises the error path to confirm the injection fires there.

`core/_build_info.py` exposes two helpers:

- `get_agent_build_info()` — returns the live three-field dict.
- `inject_agent_build(response)` — sets `response["agent_build"]`
  in place.

`knowledge/api_schema.py` gains an `agent_build` schema entry so
the field is documented and validated.

**Why important for future work**: the channel is the
foundation for the H3 canary (KDD 25) and for any future client
that needs to display server-build info to operators
(deferred follow-up §14.5).  K_H2 (24 tests) covers the
fingerprint computation invariants (deterministic across
processes, sensitive to `core/llm.py` edits, insensitive to
`RETIRED_MODELS` edits), the schema entry, and the injection
firing on both success and error paths.

### 24. skip_programs Promotion (v119.H2.1)

**Rationale**: The directive extraction prompt schema offers
two parallel ways to express "don't run program X":

- `program_settings: {phenix.X: {skip: true}}`
- `workflow_preferences: {skip_programs: ["phenix.X"]}`

The LLM, looking at the prompt, often picks the former — it's
a natural reading.  But downstream code (planner, workflow
engine) reads only the latter.  Result: the user said "skip X,"
the LLM correctly captured "skip X" in `program_settings[X].skip`,
and the program ran anyway.  The bug pre-existed v119 and was
surfaced by the existing
`tst_directive_extraction.py::skip_programs` Scenario.

v119.H2.1 adds `_promote_skip_settings_to_skip_programs()` at
the **top** of `validate_directives`, before per-setting
validation runs.  The helper:

1. Scans `program_settings.{program}.skip` for any truthy form
   per `_SKIP_TRUE_VALUES` (12 forms covering `true`, `True`,
   `"yes"`, `"on"`, `1`, etc.).
2. For each truthy hit, appends the program name to
   `workflow_preferences.skip_programs` (creating that path if
   absent).
3. Uses `.pop()` (not read-and-leave) so the LLM's
   `program_settings[X].skip` doesn't get re-interpreted later.
4. Prunes the parent `program_settings[X]` dict if it becomes
   empty after the pop.
5. Deduplicates `skip_programs` (preserves order; the first
   occurrence wins).

The promotion happens **before** per-setting validation so the
LLM's truthy-string forms don't get rejected as
`expected_type=bool` mismatches before promotion can fire.

**Why important for future work**: H2.1's order-of-operations
matters.  Future additions to `validate_directives` that need
to see post-promotion state should run after H2.1.  Future
prompt-consolidation work that closes the
`program_settings[X].skip` ambiguity at the LLM-emission layer
(deferred follow-up §8.1) would let H2.1 be removed, but until
then H2.1 is the load-bearing fix.

K_H2.1 (11 tests in `tests/tst_skip_promotion.py`) covers all
four Gemini guardrails (truthiness, .pop(), pruning, dedup),
idempotency, defensive bail on malformed input, and an
end-to-end assertion against the LLM extraction Scenario.

### 25. Startup Canary Infrastructure (v119.H3)

**Rationale**: With KDD 22's centralization and KDD 23's metadata
channel in place, v119.H3 closes the loop with two independent
consumers that probe distinct failure surfaces:

- **H3a — metadata canary (K-suite)**: tests that `agent_build`
  in a real server response matches operator-pinned expectations
  in `tests/canary_expected.json`.  Catches: wrong branch
  deployed, `core/llm.py` edited without VERSION bump, H2
  injection broken.  No LLM call, no API keys, runs in CI.
- **H3b — LLM smoke canary (operator-invoked)**: tests that the
  deployed environment can actually make an LLM call through
  the production-faithful framework.  Catches: API keys missing
  / wrong / expired, provider returning 404 on configured
  model, langchain dependency version drift, network timeouts.
  Requires API keys; not in CI.

**H3a probe design.**  The probe sends a request through `run()`
that triggers the **error path** (so `agent_build` gets injected
without invoking the LangGraph pipeline).  An earlier design
attempted to trigger this via missing required fields — but
`_parse_request` calls `apply_request_defaults` BEFORE
`validate_request`, which fills in defaults for all "required"
fields.  The current design uses a **type-mismatched** field
(`"files"` set to a string instead of a list).
`apply_request_defaults` only fills in MISSING fields, not
existing-but-wrong-typed ones, so `validate_request` catches the
type error and routes to `_build_error_response`.  A
`client_command: "CANARY_PING"` marker in the probe makes
malformed-on-purpose requests visible in server logs.

**H3b orchestrator design.**  `tests/llm/canary_check.py`
combines (1) the metadata pinning check (using
`get_agent_build_info()` directly, no LLM) with (2) a one-shot
LLM probe through `framework.call_directive_extractor` (the same
production-faithful entry point the LLM test suite uses).  The
probe is **also** registered as a `canary` Scenario in
`tst_directive_extraction.build_scenarios()`, so it can be run
independently via the standard CLI for debugging:
`phenix.python tests/llm/run_llm_tests.py --scenario canary`.

The LLM probe is wrapped in a `concurrent.futures.ThreadPoolExecutor`
context manager with a 30-second timeout.  Per Gemini's
recommendation, the timeout is enforced externally in the
orchestrator (no production-code modification).  The 30-second
threshold accounts for two LLM round-trips (intent classification
+ main directive extraction), each subject to network variance;
shorter timeouts produced false positives during real-network
testing.

**Graduated severity** (Tom's choice c): version mismatch is a
hard fail (wrong build deployed → block deploy); fingerprint
drift is a soft warning (may be intentional, e.g., a retirement
update landed mid-version); malformed fingerprint is a hard fail
(unambiguously broken).  Exit codes from H3b's CLI follow the
same gradient.

**The `canary_expected.json` lockstep**: K_H3a's
`test_version_file_matches_canary_expected` asserts that the
VERSION file and the JSON's `agent_version` field agree.  This
catches the operator-drift case where one is edited but not the
other — visible immediately rather than after a confusing test
failure on PHENIX.

**Why important for future work**: KDD 25 is the load-bearing
piece for deploy-pipeline confidence.  The deferred §14.1
"auto-startup canary" follow-up would have the server itself run
H3 on first request and log to stderr; the manual H3b workflow
covers the same operational need until that lands.

K_H3a (10 tests in `tests/tst_canary.py`) covers: config file
presence and schema, probe mechanics (error path triggers,
agent_build present, probe is fast), version pinning (deployed
matches JSON, JSON matches VERSION file), fingerprint pinning
(hard fail on malformed, soft warn on drift), end-to-end
consistency.

### 26. Preprocessor Telemetry: [STEP_1F] Marker (v119.H4 / H4.1)

**Rationale**: The LLM-based advice preprocessor is on the
live decision path; before deprecating it (Phase 2B), we need
direct evidence that the regex-based scanner extracts the same
file references the LLM extracts.  KDD 26 introduces a
telemetry marker (`[STEP_1F]`) emitted after every successful
preprocessing call that compares the LLM's file-list against
the scanner's file-list on the SAME advice text, capturing the
symmetric difference for later aggregation.

**Implementation**: `phenix_ai/run_ai_analysis.py`'s
`run_advice_preprocessing` calls a new
`agent/raw_advice_scanner.py` module after the LLM
preprocessing completes, then emits one marker line per call:

```
[STEP_1F] preprocessing_metrics preprocessor_mode=llm_comparison
  scanner_version=119.H4.1 llm_files=[...] regex_files=[...]
  in_llm_only=[...] in_regex_only=[...]
```

`scanner_version` is read from the `VERSION` file at runtime
so aggregation can segment telemetry by scanner generation.
Companion `[STEP_1F_FAILED]` marker fires if the metric block
itself raises (wrapped in defensive try/except so the marker
emission never breaks preprocessing).

**H4.1 — Golden-master corpus pinning**: a 31-document corpus
in `tests/step_1f_corpus.json` with hand-labeled file mentions.
K_H4's golden-master test pins scanner recall against this
corpus at **0.9810**, well above the **0.90** trigger
threshold from the v118 next-steps doc §3.7.  Phase 2B
activation is **unblocked** as a result.

**Why important for future work**: KDD 26 is the empirical
basis for the Phase 2B preprocessor deprecation decision.
The marker emission is also the first H-series feature
deliberately designed for a future relay-channel ship (KDD 27)
— H4 emitted to server stderr only; H5 surfaced it
client-side.

### 27. diagnostic_messages Relay Channel (v119.H5 / H5.1)

**Rationale**: Server-side `[STEP_1F]` markers (KDD 26) emit
to the server's stderr.  When the agent dispatches to
`ai.phenix-online.org`, that stderr isn't relayed to the
client — the operator sees the marker locally (when
`run_on_server=False`) but not remotely.  KDD 27 closes the
gap by adding a `diagnostic_messages` list field on
`working_results` that flows through both dispatch modes
identically.

**The uniformity criterion** (introduced in H5 plan rev 5
after a uniformity-violation in rev 4): all actions should
take place identically whether analysis runs on the server or
locally.  Local execution is just "the network round-trip
happens in memory."  Same code paths, same data flow, same
return contract.  Rev 4 had violated this by gating the
client's re-emit on dispatch mode; rev 5 removed the gate
and the helper's stderr write, making local and remote modes
produce identical operator-visible output.

**Architecture**: three layers, each contributing one
invariant:

1. **Engine**: `phenix_ai/run_ai_analysis.py`'s engine
   functions (`run_advice_preprocessing`,
   `run_directive_extraction`, `run_failure_diagnosis`) each
   have a function-local `diagnostic_messages = []` (no
   mutable defaults), threaded through every return path.
   A `_emit_marker(list, str)` helper appends to the list;
   it never writes to stderr and never raises.
2. **Transport (Total Initialization Policy)**: every
   `working_results` from `programs/ai_analysis.py::get_results_from_all`
   gets `diagnostic_messages=[]` regardless of mode.
   `get_results_as_JSON` serializes the field unconditionally;
   `_process_server_success` decodes it unconditionally on
   the client.  Local-mode `_run_*_locally` helpers thread
   markers from the engine into `working_results`.
3. **Client (centralized re-emit)**: a single
   `_relay_diagnostic_messages_to_stderr` helper in
   `programs/ai_analysis.py::run_job_on_server_or_locally`
   re-emits every list entry to stderr after the dispatcher
   returns.  Adding new markers requires NO client-side
   changes — confirmed by H5.1 which added three markers
   without touching `programs/ai_agent.py`.

**Markers covered**:
- v119.H5: `[STEP_1F]`, `[STEP_1F_FAILED]` (from KDD 26's
  preprocessing metric block)
- v119.H5.1: `[ADVICE_PREPROCESSING_FAILED]` (outer
  exception handler of `run_advice_preprocessing`),
  `[DIRECTIVE_EXTRACTION_FAILED]` (outer exception handler
  of `run_directive_extraction` — distinct from the inner
  one in `agent/directive_extractor.py` which still emits
  to stderr only), `[FAILURE_DIAGNOSIS_FAILED]` (outer
  exception handler of `run_failure_diagnosis`).

The H5.1 markers fire only on UNEXPECTED exceptions.
Configuration mistakes (invalid provider, missing API key)
are handled gracefully upstream by `validate_api_keys` and
`setup_llms` and don't reach these handlers.  This is by
design — these markers are catastrophic-failure telemetry,
not configuration-error reporting.

**Backward compatibility (Gemini rev 2 mandate)**: a server
without H5 returns no `diagnostic_messages` field; the
client's defensive unpack returns `[]`; the centralized
re-emit produces no output; the operator's experience is
indistinguishable from the H4 baseline.  Verified by K_H5 §C
`test_old_server_response_handled_defensively` and
`test_corrupted_payload_decoded_defensively`.

**Test strategy**: K_H5 organized into 6 sections at H5.1
(§A helper unit, §B field-in-return, §C backward compat,
§D uniform client re-emit, §E production encode/decode
round-trip, §F extended markers).  §F exception tests use
deterministic monkey-patching (patching `validate_api_keys`
or `setup_llms` to raise) plus Gemini Q2 dual-assertion
(marker present AND debug_log evidence of exception path)
for robustness against future upstream graceful-handling
refactors.

26 tests total; libtbx-dependent tests SKIP cleanly in
sandbox and PASS under PHENIX.

**Why important for future work**: KDD 27 is the
architectural foundation for any future operator-visibility
work.  Adding a new marker is a one-line change at the
engine site — no client-side, no protocol changes, no
backward-compat concerns.  H5.2 (queued, not started) will
extend the channel to three remaining stderr-only markers in
`agent/directive_extractor.py`.

### 28. Canary Probe Robustness Across Providers (v119.H5.1)

**Rationale** (surfaced during H5.1 verification): the H3b
canary probe text `"Run phenix.refine on the model"` produced
false-fail on OpenAI.  OpenAI faithfully extracted no
parameters (the input mentions phenix.refine but no settable
parameters), returning
`{"program_settings": {"phenix.refine": {}}}`.
`validate_directives` (correctly) stripped the empty inner
dict, yielding `{}`.  The canary asserts non-empty, so it
failed even though OpenAI behaved correctly.

Google tended to over-extract (added inferred fields like
`start_with_program`), masking the issue.

**Fix**: change the probe text to `"Run phenix.refine with
resolution 2.5"` — one concrete settable parameter that any
LLM should extract, ensuring the validated dict is non-empty
across providers.  Same canary purpose (provider
reachability), now robust to provider compliance variance.

**Why important for future work**: KDD 28 is a small-scope
example of a broader pattern — canary probes should include
at least one concrete extractable signal, not just an action
verb plus an entity name.  The natural extension for v120 is
to verify the same property holds for other LLM-probing
tests as Phase 2A planning-suite scaffolding lights up.

### 29. Sidecar-Aware Merge + Type-Gated Boolean Coercion (v119.H5.1.1)

**Rationale**: Two latent bugs in
`agent/directive_extractor.py` surfaced during the §2.1
review of v118.9 leftovers.  Both are subtle interactions
between the multi-source directive-extraction architecture
(LLM extraction + simple-extraction regex + experiment-type
correction) and Python's coercion semantics.

**Bug 1 — Stale `_corrected_from` sidecar**: when LLM
extraction returns a wrong-for-experiment-type
`after_program`, `_apply_experiment_type_program_reprints`
corrects it and attaches a `_corrected_from` sidecar
describing the transition.  When `merge_directives` then
combines the corrected directives with simple-extraction
output (which independently re-extracts the original
`after_program` from raw advice), the override-wins
dict-merge had two pathologies:

1. *Correction reverted*: simple-extraction's value matched
   `_corrected_from.from`, the merge silently un-corrected
   the field; sidecar stayed around but no longer described
   a true correction.
2. *Zombie metadata*: simple-extraction (or any other
   override) set `after_program` to a THIRD value (neither
   the corrected `to` nor the original `from`), leaving
   the sidecar with stale `from`/`to` fields that
   contradicted the current value.

**Fix 1**: detect both cases inside `merge_directives`
itself (semantic belongs in the merge helper, not at
callers).  Decision matrix:

| Scenario | Action |
|---|---|
| Override `after_program` matches `_corrected_from.from` | Strip from override; preserve correction; keep sidecar |
| Override `after_program` is a third value, no own sidecar | Let override win; clear base's stale sidecar |
| Override brings its own `_corrected_from` | Standard dict-merge (override's sidecar wins) |
| Override doesn't touch `after_program` | Standard dict-merge |
| Base has no `_corrected_from` | Standard dict-merge (no protection needed) |

Implementation uses dict-comprehension rebuilds to preserve
the existing no-mutation-of-inputs invariant in
`merge_directives` (same pattern as the program_settings
deep-merge above already uses).

**Bug 2 — Boolean list-coercion**: `bool([False])` returns
`True` in Python because non-empty list is truthy.  Two
sites in `validate_directives` did bare `bool(value)` on
boolean preferences:
- `prefer_anomalous`, `prefer_unmerged`, `prefer_merged`
  in `file_preferences`
- `use_experimental_phasing`, `use_molecular_replacement`,
  `use_mr_sad`, `model_is_placed`,
  `wants_validation_only` in `workflow_preferences`

If the LLM emits `[false]` (list-wrapped — happens
occasionally with provider/temperature variance), the
bare `bool()` would silently flip semantics.

**Fix 2**: explicit type-gated inline unwrap at both sites:

```python
if isinstance(_v, list) and len(_v) == 1 and isinstance(_v[0], bool):
    _v = _v[0]
if isinstance(_v, bool):
    valid_prefs[key] = _v
elif isinstance(_v, int) and not isinstance(_v, bool):
    valid_prefs[key] = bool(_v)
elif isinstance(_v, str):
    valid_prefs[key] = bool(_v)
else:
    _log(...)
```

The `isinstance(_v[0], bool)` guard bounds the unwrap to
genuine `[bool]` — `[1]`, `[True, False]`, `["true"]`, etc.
fall through to the drop+log branch.  Future maintainer
adding a list-typed key to either boolean tuple wouldn't
trigger silent flattening (the value would drop+log, which
surfaces the mistake during testing).

Explicit pattern chosen over the existing v118.10
`_coerce_setting_value` helper.  Gemini's plan review
flagged the helper's broader semantics (string-join for
multi-element lists, dict rejection, etc.) as future-
maintenance risk.  The inline pattern is 4 lines that
obviously match the bug they fix.

**Pre-existing wrongness NOT fixed**: `bool("false")`
returns `True` in Python — a separate string-truthy quirk
preserved for non-list inputs.  Separate bug class; if it
manifests in production, separate ship.

**Test strategy**: K_H5_1_1 organized into 2 sections:
§A merge_directives sidecar protection (7 tests), §B
validate_directives boolean list-wrap defense (16 tests).
Both sections include no-mutation-of-inputs invariant
tests verifying the dict-comprehension rebuilds don't
accidentally mutate caller dicts.

23 tests total; no libtbx dependency (`merge_directives`
and `validate_directives` are pure-Python helpers, so all
tests PASS in sandbox and under PHENIX without SKIPs).

**Why important for future work**: KDD 29 demonstrates two
patterns worth applying elsewhere in the cluster:

1. *Semantic-protection helpers belong in the merge layer,
   not at callers*.  Both production callers of
   `merge_directives` benefit from the sidecar protection
   automatically; neither had to be touched.  Future
   merge-class helpers should follow this pattern.
2. *Explicit type-gates for narrow scopes beat
   general-purpose helpers when the cost of helper-coupling
   exceeds the saved lines.*  The `_coerce_setting_value`
   helper is the right tool for its use case
   (`program_settings` and `stop_conditions` value coercion
   with multiple expected types) but a heavier weapon than
   needed for the narrow "[bool] → bool" unwrap.

### 30. Phase 2A Planning-Suite Framework Patterns (v119.H6)

**Context**: H3b shipped `tests/llm/framework.py` with
`call_directive_extractor` implemented and three sibling
functions as stubs.  The 8 planning scenarios in
`tst_planning.py` ERRORed at 0.0s because the stubs raised
`NotImplementedError`.  H6 replaced the stubs with real
implementations, unblocking ~12 minutes of planning-LLM
reliability testing per Mac run.

**Decision**: Four design patterns that together make the
planning framework production-faithful without coupling it
to the production class hierarchy:

#### 30.1 Production parity, not test convenience

`is_stop_intent(intent)` checks BOTH signals that production
checks at `graph_nodes.py:~2168`:
```python
return intent.get("program") == "STOP" or intent.get("stop") is True
```

The strict `is True` (rather than truthy) matches
production exactly.  Subtle but important: a test framework
that accepted `stop: "yes"` or `stop: 1` would be more
permissive than production and could mask future regressions
where production tightens up the check.  The K_H6 §A tests
pin both signals AND the strictness.

#### 30.2 Lazy imports for sandbox testability

The framework's planning functions need `BedrockChatBuilder`,
the rate-limit handler, and PHENIX configuration utilities.
Importing them at module-load time would make the framework
unimportable in sandbox.  Solution: lazy imports inside the
function body, so the test framework loads fine but the
actual planning call requires PHENIX.

Pattern:
```python
def call_planning_llm(state_inputs, provider="google"):
    # Lazy imports — only fire when called, not at module load
    from langchain.agent.graph_nodes import plan as _plan_function
    from libtbx.langchain.agent.rate_limit_handler import (
        get_google_handler)
    ...
```

Sandbox: framework imports OK; K_H6 §D tests SKIP gracefully
because the production dependencies aren't available.

PHENIX: K_H6 §D tests also SKIP — but for a DIFFERENT
reason.  Under PHENIX, the libtbx import succeeds, so the
sandbox `sys.modules` fakes never get a chance to be loaded
by the test.  The mechanic is identical in both
environments; only the gating differs.  K_H6 §D documents
this and pins what it CAN pin (the sandbox version's
mechanic).

#### 30.3 raw_output preservation invariant (Gemini rev-3 critical fix)

When `parse_intent_json` parses the planning LLM's output:
```python
# WRONG (rev 1, rev 2)
try:
    raw_output = ...
    intent = parse_intent_json(raw_output)
except Exception:
    return {"raw_output": ???, ...}  # raw_output undefined!

# RIGHT (rev 3)
raw_output = ...  # ALWAYS bound, BEFORE the try
try:
    intent = parse_intent_json(raw_output)
except Exception:
    return {"raw_output": raw_output, ...}  # always available
```

Hoist `raw_output = ...` BEFORE the try block.  If
`parse_intent_json` raises (JSON parse failure on malformed
LLM output), the exception handler returns the LLM's
original malformed text for diagnosis instead of crashing
the test framework or returning empty string.

K_H6 §D pins this with a test that mocks
`parse_intent_json` to raise, asserts the result contains
the original raw text.

#### 30.4 K-test SKIP-under-PHENIX pattern for monkey-patched modules

`tests/tst_planning_framework.py` §D monkey-patches
`sys.modules['langchain.agent.graph_nodes']` to control
what `call_planning_llm` invokes.  Under sandbox this
works fine (no libtbx in play).  Under PHENIX, the
libtbx-version of the import wins, so the monkey-patch
doesn't take effect.

Resolution: §D tests check for libtbx availability and
SKIP cleanly under PHENIX.  The mechanic is the same in
both environments; sandbox tests pin the contract.  Live
planning tests under PHENIX exercise the actual production
path.

This is the OPPOSITE pattern from H3a / H4 / H5 where
sandbox SKIPs and PHENIX PASSes.  Both patterns are valid
— the question is "where does the import resolution prevent
the test mechanism from working?"

**Why important for future work**: any future framework
function that needs to monkey-patch a production module
should follow §30.4's pattern (SKIP under PHENIX, full
coverage under sandbox).  Production behavior is exercised
by the live LLM tests directly; the K-test pins the
mechanism.

### 31. Scanner-First Activation Patterns (v119.H7)

**Context**: H4.1 measured `raw_advice_scanner`'s recall at
0.9810 against the LLM extraction's output, exceeding the
0.90 threshold from PHASE2_PLAN for Phase 2B activation.
H7 activated Scope B: scanner becomes the primary source
for `extracted_files`, with the LLM extraction kept as
fallback.

**Decision**: Four design patterns informed by Gemini
critique cycles + a critical implementation-time
correction.

#### 31.1 Q2-strict: consumer contract preservation

When swapping the primary file-extraction source, the
output contract of the new source must match the existing
consumer's expectation EXACTLY.  Inspection of
`programs/ai_agent.py:8128-8174` revealed:

- The consumer's own comment at line 8142 names the field
  "Files mentioned in advice" — that's the semantic
  contract.
- Auto-promotion to `original_files` at lines 8157-8158
  means polluted `extracted_files` directly inflates the
  agent's active input set.

Decision: scanner called with `None` hint (NOT
`file_list`).  Preserves "files mentioned in advice text
only" semantics.  Passing `file_list` would have silently
inflated `extracted_files` with workspace files the user
never mentioned — surfacing as "agent analyzed files I
didn't ask about" in production.

**Pattern**: any source-swap should be preceded by
inspection of every downstream consumer to verify contract
compatibility.  Code inspection beats plan inspection.

#### 31.2 Defense-in-depth fallback with libtbx → relative imports

The new primary path can fail; the fallback path must
remain reliable.  Pattern from H7's
`run_ai_analysis.py:1009-1069`:

```python
try:
    try:
        from libtbx.langchain.agent.raw_advice_scanner import (
            scan_files_in_advice)
    except ImportError:
        from agent.raw_advice_scanner import scan_files_in_advice
    extracted_files = scan_files_in_advice(raw_advice, None)
    ...
except Exception as e:
    debug_log.append(f"raw_advice_scanner failed: {e}; ...")
    if processed_advice and processed_advice != raw_advice:
        try:
            try:
                from libtbx.langchain.agent.advice_preprocessor import (
                    extract_files_from_processed_advice)
            except ImportError:
                from agent.advice_preprocessor import (
                    extract_files_from_processed_advice)
            extracted_files = extract_files_from_processed_advice(...)
        except Exception as e2:
            debug_log.append(f"LLM fallback extraction also failed: {e2}; ...")
```

The Gemini rev-2 critical fix: the FALLBACK import path
must ALSO use libtbx → relative.  Rev 1 used a bare libtbx
import in the fallback, which would silently fail in
sandbox (or any environment where libtbx isn't available)
and be swallowed by the outer except handler — a silent
failure of the safety net.

**Pattern**: any "safety net" import must mirror the
primary path's import-fallback pattern.  Inconsistency
defeats the purpose.

#### 31.3 Telemetry-axis independence after primary-source change

H4 designed `[STEP_1F]` to compare LLM extraction vs scanner
extraction.  Pre-H7 this was implicit: `extracted_files` was
the LLM's output, `regex_files_raw` was the scanner's.
Post-H7, `extracted_files` is the scanner's output.

**Trap**: if the [STEP_1F] block reuses `extracted_files`
as the LLM-comparison axis, it would compare
scanner-to-scanner and the recall numbers would always be
inflated.  Caught during implementation review.

Resolution: the [STEP_1F] block now calls
`extract_files_from_processed_advice` SEPARATELY from the
main extraction path:
```python
_llm_extracted = []
try:
    try:
        from libtbx... import extract_files_from_processed_advice
    except ImportError:
        from agent... import extract_files_from_processed_advice
    _llm_extracted = extract_files_from_processed_advice(processed_advice)
except Exception:
    _llm_extracted = []
llm_normalized = _step_1f_normalize(_llm_extracted)  # NOT extracted_files
```

Cost: one extra pure-regex call per request (microseconds).
Benefit: the telemetry's diagnostic value (scanner vs LLM
drift detection) is preserved.

**Pattern**: when changing the SOURCE of a variable that
also drives telemetry, audit every downstream consumer of
that variable to verify the source-change doesn't
invalidate the telemetry's comparison semantics.

K_H7's `test_telemetry_independence_from_extracted_files`
pins this property at the source level — a future refactor
that accidentally reverts to reusing `extracted_files` is
caught immediately.

#### 31.4 Recall metrics with explicit naming and zero-division guards

Two new inline metrics added to `[STEP_1F]`:
- `scanner_recall_against_llm` (`|intersection|/|llm_files|`)
- `llm_recall_against_scanner` (`|intersection|/|regex_files|`)

Explicit semantic names (not `recall` with implicit
direction).  Future readers know the numerator without
checking docs.

Both metrics zero-division-guarded: empty denominator → 1.0
(mathematically sound for recall — 100% of zero targets
captured — and avoids `ZeroDivisionError` in production on
advice with no file mentions).

**Pattern**: any inline metric must (a) have an explicit
semantic name, (b) be zero-division-guarded with a
documented degenerate value, and (c) appear in [STEP_1F]
without modifying existing fields (additive-only).

**Why important for future work**: KDD 31 demonstrates the
discipline needed for primary-source swaps in production
code paths.  Three of the four sub-patterns (Q2-strict,
import-fallback consistency, telemetry-axis independence)
came from defending against silent failure modes.  The
fourth (explicit metric naming) prevents future telemetry
schema drift.

The same discipline applies to any future Phase 2B Scope C
work that removes the LLM preprocessor entirely.  When
`llm_files=[]` always, the recall metrics degrade
gracefully but tagging via `preprocessor_mode=regex_only`
becomes essential for aggregation tooling.

### 32. Template-Literal Allowlist Pattern (v119.H8)

**Context**: An LLM-emitted strategy for `phenix.autobuild_denmod`
was reaching the post-build `sanitize_command` call with
`maps_only=True` stripped out, despite the parameter being a
legitimate part of the program's command template in
`programs.yaml`.  Without `maps_only=True`, `autobuild_denmod`
tries to rebuild the model rather than just producing density-
modified maps — wrong program semantics.

**Decision**: the per-program parameter allowlist used by
`sanitize_command::_load_prog_allowlist` to enforce Rule D
("strip unknown program-specific parameters") must include
**literal `key=value` tokens baked into the command template**,
not just the entries declared in `strategy_flags`.

#### 32.1 Template-literal tokens are program invariants, not LLM strategies

`programs.yaml` distinguishes between two kinds of command-line
parameters for any given program:

- **Strategy-tunable parameters**, declared under
  `strategy_flags` in the program definition.  These are
  parameters the LLM is permitted to set as part of a strategy
  dict (e.g., `nproc`, `resolution`).  The allowlist already
  derived these correctly.
- **Template-literal parameters**, baked directly into the
  `command:` template string (e.g.,
  `phenix.autobuild_denmod ... maps_only=True`).  These are
  program invariants — fixed properties of how the agent invokes
  the program for a given workflow purpose.  The LLM has no
  business changing them; the strategy system has no entry for
  them.

Pre-H8, the allowlist derived only from `strategy_flags`, so
Rule D treated template-literal parameters as "unknown" and
stripped them.  The autobuild_denmod case was the visible
production symptom; the bug class affects any program with
template-literal parameters not also declared in `strategy_flags`.

**Pattern**: when sanitization rules consult a configuration
artifact to decide what to allow, that artifact must enumerate
**all** legitimate parameters — even ones the user/LLM has no
business setting.  An invariant baked into a template is still
a parameter that will appear on the command line; the
sanitizer must know about it.

#### 32.2 Placeholder vs literal in template extraction

The fix walks `prog_def.get('command', '')` with the regex
`r'(?:^|\s)([a-zA-Z_][a-zA-Z0-9_.]*)='` and adds each matched
parameter name (leaf-only, lowercased) to the allowlist.

The `(?:^|\s)` prefix is load-bearing: it requires start-of-
string or whitespace before the identifier, which is precisely
what distinguishes a **literal parameter token** like
` maps_only=True` (preceded by whitespace) from a
**placeholder substitution target** like `{data_mtz}=...`
(preceded by `{`).  The `{` is neither start-of-string nor
whitespace, so the regex doesn't match — and placeholder names
correctly do NOT pollute the allowlist.

This design treats the command-template's text as a parseable
sublanguage with two distinct token classes:

- `KEY=VALUE` at word-boundary → parameter (add KEY's leaf to
  allowlist).
- `{PLACEHOLDER}` inside the template → substitution slot, NOT
  a parameter.

K_H8's `test_bug8_template_literal_extraction_ignores_placeholders`
pins this regex contract.

**Pattern**: in any DSL with both literal tokens and
substitution syntax, the consumer that needs to distinguish
them must anchor on the literal form's boundary properties
(whitespace, punctuation), not just on the identifier shape.

#### 32.3 Program-invariant override at the registry, not the sanitizer

A defense-in-depth concern: now that `maps_only` is in the
allowlist, what stops an LLM from emitting
`strategy.maps_only=False` to override the template's
`maps_only=True`?  Answer: `ProgramRegistry.build_command`
already drops unknown strategy entries — and `maps_only` is
NOT in `phenix.autobuild_denmod`'s `strategy_flags`, so the
LLM's strategy entry IS unknown and IS dropped (with a
"Unknown strategy 'maps_only'" log line).

So the contract is split:

- **`_load_prog_allowlist` / `sanitize_command`** (post-build):
  "preserve template-literal parameters from being stripped as
  unknown."  H8's fix.
- **`ProgramRegistry.build_command`** (build): "drop LLM
  strategy entries that aren't in `strategy_flags`."  Pre-existing
  behavior.

Together they enforce: template invariants reach the command
line; LLM cannot override them via strategy.

K_H8's `test_bug8_adversarial_strategy_override_dropped` (Claude-
reviewer-suggested) pins this two-step invariant.  The test
emits `strategy={"maps_only": False}`, verifies the registry
logs "Unknown strategy 'maps_only'" and drops it, and confirms
the resulting command still contains the template-literal
`maps_only=True` only.

**Pattern**: when defense-in-depth is split across two
processing stages, each stage's tests should assert both its
own contract and the joint outcome.  Don't rely on the reader
to compose the invariants mentally.

**Why important for future work**: KDD 32 surfaces a class of
issue that's invisible until a specific program's specific
parameter is checked — most programs have parameters that pass
through cleanly because their command template happens to fall
into a path the allowlist loader handles.  When adding any new
program with template literals in its command template, run K_H8's
allowlist-loader test against the new program's parameters as a
smoke check.  Generalization candidate: a property test that
verifies the allowlist round-trips for every program in
`programs.yaml`.

### 33. Per-Program PHIL Scope Rewriting (v119.H9)

**Context**: `phenix.predict_and_build` was being given commands
with `crystal_symmetry.unit_cell=...` but the program expects
`crystal_info.unit_cell=...`.  PHENIX rejected with "Some PHIL
parameters are not recognized".  The wrong scope was being
emitted because `agent/program_registry.py` (around line 778)
unconditionally prepends `crystal_symmetry.` to bare unit_cell
and space_group parameters, and the per-program rewrite table
in `knowledge/parameter_fixes.json` — while containing other
entries for `phenix.predict_and_build` (file-name remappings,
`nproc` scopings) — had no scope-rewrite entries to redirect
the `crystal_symmetry.*` parameters to `crystal_info.*` for
this specific program.

**Decision**: the architectural layering for per-parameter
PHIL scope handling is **two-layer with explicit per-program
overrides**:

1. **Layer 1 — uniform default scoping** (`program_registry.py`):
   parameters are prepended with the most common scope
   (`crystal_symmetry.` for unit_cell / space_group).  This is
   correct for the majority of PHENIX programs.
2. **Layer 2 — per-program rewrite table**
   (`knowledge/parameter_fixes.json`): for programs whose PHIL
   scope differs from the default, an explicit entry maps the
   default-scoped parameter name to the program-specific name.

The fix is purely additive at Layer 2 — four new rewrite entries
plus one annotation under the existing `phenix.predict_and_build`
block in the JSON:

```json
"phenix.predict_and_build": {
  // pre-existing entries (file-name remapping, nproc scoping)
  "data_file": "xray_data_file",
  "input_files.data_file": "input_files.xray_data_file",
  "control.nproc": null,
  "prediction.nproc": null,
  // v119.H9 additions
  "_comment_cs": "predict_and_build uses crystal_info.* scope (v119.H9)",
  "crystal_symmetry.space_group": "crystal_info.space_group",
  "crystal_symmetry.unit_cell":   "crystal_info.unit_cell",
  "xray_data.space_group":        "crystal_info.space_group",
  "xray_data.unit_cell":          "crystal_info.unit_cell"
}
```

The `_comment_cs` key is skipped by `fix_program_parameters`
(which ignores keys starting with `_`), so it serves as inline
documentation for the rewrite table.  The `xray_data.*` fallback
entries (Gemini-suggested) handle the case where an LLM emits a
different incorrect scope; both get rewritten to `crystal_info.*`.

#### 33.1 Why a rewrite table, not a fix in the registry

An alternative design would have been to change Layer 1 to
look up the correct scope per program from `programs.yaml` —
removing the need for a separate rewrite table.  H9 deliberately
did not go this route:

- The default-and-override design is **incremental and
  reversible**.  Each program's exception is one JSON entry;
  reverting is one line.  A Layer 1 rewrite would touch
  registry code and risk regressions for programs not yet
  covered by tests.
- The rewrite table is **self-documenting**.  Reading
  `parameter_fixes.json` reveals exactly which programs deviate
  from the default scope, which is itself a maintenance signal.
- The rewrite table is **runtime-checkable** without server
  restart in dev — testing a fix is editing a JSON file.

**Pattern**: when the majority case of a configuration is
correct and a minority needs deviation, prefer a per-exception
table over generalizing the rule.  The exception table is the
contract surface for the deviation.

#### 33.2 K_H9 module-cache invalidation

`parameter_fixes.json` is loaded by `agent/planner.py` on first
call to `get_parameter_fixes()` and cached at module scope via
`_PARAMETER_FIXES = None` (sentinel) → populated dict on first
load.  K_H9's tests need to exercise different file states
(e.g., test 4 verifies idempotence by ensuring no double-rewrite
of an already-correctly-scoped command).  Solution: a
test-module-level `_reload_parameter_fixes()` helper that resets
the cache sentinel back to `None`, forcing the next call to
re-read the JSON file:

```python
def _reload_parameter_fixes():
    """Force re-read of parameter_fixes.json (defeat the module cache)."""
    try:
        from libtbx.langchain.agent import planner as _p
    except ImportError:
        try:
            from agent import planner as _p
        except ImportError:
            return
    _p._PARAMETER_FIXES = None
```

All four K_H9 tests call this at the top.  The mechanism is
deliberately direct: `_PARAMETER_FIXES` is a public-by-convention
module attribute (no leading underscore on the sentinel access
pattern; the leading underscore on the name is hygienic, not a
strict private marker), and resetting it to `None` triggers
`get_parameter_fixes`'s built-in lazy-load path.  No `importlib`
reload is needed — and would be wrong, since reloading the
planner module would invalidate `fix_program_parameters` references
held by other tests.

This pattern generalizes: any K-test that depends on
JSON-or-YAML-loaded data with module-level lazy caching can
expose a similar single-line invalidation helper.

**Pattern**: when a test exercises configuration data that's
lazily cached at module scope, expose an invalidation hook
rather than relying on test-isolation order.  Cached data +
test order dependence = flaky tests.  Direct cache-sentinel
reset is preferable to `importlib.reload` because the latter
also invalidates function references and re-runs module
initialization.

**Why important for future work**: KDD 33 establishes the
contract surface for adding new program-specific PHIL deviations.
When a new PHENIX program is added with non-default scoping
(check by attempting to run a typical command and observing
PHIL rejections), the fix path is: add an entry to
`parameter_fixes.json` + add a K-test in `tst_autosol_bugs.py`
following the K_H9 pattern.  No registry or planner code change
is needed for the common case.

### 34. `exclude_patterns` Word-Boundary Semantics (v119.H10 / H11)

**Context**: AIAgent_62 cycle 7 emitted a `phenix.refine` command
with `refine_001_001.cif` (a model mmCIF from an earlier refine
step) injected as a positional third argument.  PHENIX
interpreted it as a SECOND model and crashed with "wrong number
of models".  The `phenix.refine::ligand_cif` slot in
`programs.yaml` declares `exclude_patterns: [..., "refine_", ...]`
specifically to reject such files; the schema designer's intent
was clear, but the filter wasn't taking effect.

The fix split across two ships because two independent bugs
were stacked.  Both are needed; H10 alone does not fix the
cycle-7 crash.

#### 34.1 Structural application of `exclude_patterns` across all selection paths (H10)

`CommandBuilder._find_file_for_slot` in `agent/command_builder.py`
has multiple file-selection paths (PRIORITY 1 through 4) plus a
refinement pre-population path in `_select_files`.  The in-code
comment that H10 left at line 1279 documents the pre-H10 state:

> Pre-H10, `exclude_patterns` was only consulted at PRIORITY 4
> (extension fallback).  PRIORITY 2 (best_files), PRIORITY 2.5
> (recovery strategies), PRIORITY 3 (category-based), and
> PRIORITY 3.5 (fallback best_files) all bypassed it.

(`_select_files` also has an LLM-selected files path (line 963)
with its own `exclude_patterns` check; whether that check predates
H10 is not marked in the code.  H10's helper covers the auto-fill
paths inside `_find_file_for_slot`, which is where the bug class
lives.)

The AIAgent_62 cycle-7 build picked the bad file via **PRIORITY
3 (category-based)** — one of the paths that bypassed the
filter.  Since cycle 7's LLM did not explicitly select a
ligand_cif file, auto-fill ran and grabbed `refine_001_001.cif`
out of the `ligand_cif` category without checking the slot's
`exclude_patterns`.

**Decision**: define a closure `_matches_exclude_patterns_h10(f)`
at the entry of `_find_file_for_slot` (line 1292) and apply it
at every auto-fill site:

```python
_exclude_patterns_h10 = input_def.get("exclude_patterns", [])
def _matches_exclude_patterns_h10(f):
    if not _exclude_patterns_h10:
        return False
    return matches_exclude_pattern(
        os.path.basename(f), _exclude_patterns_h10)
```

Applied at 7 sites inside `_find_file_for_slot` (PRIORITY 2
best_files at line 1333, PRIORITY 2.5 recovery_strategies at
line 1359, PRIORITY 3 category variants at lines 1406/1417/1437,
PRIORITY 3.5 `require_best_files_only` at line 1472, PRIORITY
3.5 fallback best_files at line 1513) plus 1 site in
`_select_files` (line 698) for refinement pre-population
(`best_files["model"]` honoring exclude_patterns from the model
slot — uses `matches_exclude_pattern` directly rather than the
closure, since `_select_files` doesn't have the closure in
scope).  Total: 8 application sites, 9 `v119.H10` markers
(8 application sites + 1 marker on the helper-definition comment
block).

**Pattern**: when a slot-level filter must apply consistently
across multiple selection algorithms, the filter is the contract;
each algorithm-specific code path must honor it.  Closures
provide the necessary lexical scoping for the input_def-derived
list without requiring helper functions to be passed through
every signature.

(Subsequent helper-consistency refactor — extract module-level
`_file_passes_exclude_patterns(file_path, input_def)` and replace
all 8 sites — is tracked as a deferred follow-up; the H10 closure
form is correct but DRY can be improved.)

#### 34.2 Word-boundary matching semantics in `matches_exclude_pattern` (function-side, intentional)

`agent/file_utils.py::matches_exclude_pattern` uses **word-boundary
regex matching**, not substring matching:

```python
re.search(
    r'(?:^|[_\-\.])' + re.escape(pat_stem) + r'(?=[_\-\.]|$)',
    stem)
```

The pattern must be preceded by start-of-string or one of
`_`/`-`/`.`, AND followed by one of `_`/`-`/`.`/end-of-stem.
(Extension-suffix patterns like `"lig.pdb"` are special-cased:
the basename must end with `.pdb` AND the stem-portion must
satisfy the word-boundary regex against `lig`.)

**Why word-boundary**: prevents false positives like `"ligand"`
matching `"noligand.pdb"`, where the pattern's intent is the
ligand-output naming convention (`ligand.pdb`, `my_ligand.pdb`),
not arbitrary substring containment.  Substring matching would
have semantic gotchas at scale.

**Why this is a design decision, not just a function**: the
semantics constrain the YAML pattern authoring grammar.  Patterns
must be authored TO this semantics.  Patterns authored AGAINST
this semantics (assuming substring) silently fail to match —
exactly what happened with `"refine_"` (trailing underscore
breaks the regex because the next char in the filename is `0`,
not a boundary char).

#### 34.3 YAML pattern authoring grammar (H11)

H11 corrects three YAML patterns that had been authored against
substring semantics:

- `phenix.refine::ligand_cif`: `"refine_"` → `"refine"`
- `phenix.mtriage::full_map`: `"_half"` → `"half"`
- `phenix.real_space_refine::map`: `"_half"` → `"half"`

H11 also embeds a `DESIGN NOTE (exclude_patterns)` block at the
top of `knowledge/programs.yaml` documenting the grammar for
future YAML authors:

- **Rule 1**: Patterns must NOT have leading or trailing
  underscore.  `"refine"` is correct; `"refine_"` and `"_refine"`
  are broken.
- **Rule 2**: To catch alphanumeric-suffix variants where a
  digit immediately follows the pattern stem with no separator
  (e.g., `protein_half1.ccp4`), include explicit suffix patterns
  alongside the bare pattern: `[half, half1, half2, half_1, half_2]`.
  Bare `"half"` does NOT match `protein_half1.ccp4` because `1`
  is not a boundary char — Gemini's plan review caught this and
  required the numeric defense-in-depth patterns be preserved.

**Pattern**: when a function's semantics constrain a configuration
grammar, the constraint must be authored as inline documentation
at the grammar's authoring site.  Function-level docstrings are
necessary but not sufficient — config authors don't open the
function source.

#### 34.4 Semantic-pin tests for sandbox-vs-real-function divergence (H11)

The H10 → H11 cycle exposed a sandboxing failure mode worth
capturing as policy.

The H10 sandbox stub for `matches_exclude_pattern` did substring
matching (the wrong semantics).  H10's K-suite passed in sandbox
because the stub-vs-YAML pattern interaction happened to give
the right answer for the test files.  Against the REAL function,
the same YAML pattern (`"refine_"`) returned False — H10's
filter fired but rejected nothing.  H10 packaged as "fixes the
cycle-7 bug" — it didn't.

**Decision**: every function with non-obvious semantics (regex
matching, parsing, escaping, normalization) that's stubbed in
sandbox should have a paired semantic-pin test that calls the
REAL function with documented expected behavior.
`test_bug11_matches_exclude_pattern_semantics` in
`tst_autosol_bugs.py` is the template — 10 named cases asserting
specific function behavior, with a graceful sandbox-skip fallback
when the real function isn't importable.

**How the divergence is caught going forward**:
- In production environments (Tom's PHENIX), the real function
  is importable; the test runs against it; semantic pins assert
  the documented behavior.  If the function ever changes, the
  test fails.
- In sandbox environments with a stub, the test still runs
  against whatever `matches_exclude_pattern` resolves to.  If
  the stub has substring semantics, the trailing-underscore
  assertion (`assert_false(real_func("refine_001.cif", ["refine_"]))`)
  fails because the substring stub returns True.  The divergence
  surfaces immediately.

**Pattern**: when sandbox stubs exist, the parity check is implicit
in any test that asserts function behavior — provided the test
asserts behavior the stub-vs-real difference would cross.  The
H11 test does this deliberately, including the trailing-underscore
case that specifically distinguishes substring from word-boundary
matching.

(A formal differential parity test — import both stub and real,
assert equal returns on a fixture — is tracked as deferred Item
2.4 in `v119_H11_PLAN_rev2.md`.  It requires a stub-package
artifact that doesn't currently exist as a real PHENIX tree
member; the semantic-pin test is the actionable substitute for
now.)

**Why important for future work**: KDD 34 captures the discipline
needed when configuration grammars constrain function semantics
and vice versa.  The four sub-patterns (structural application of
filters, function-side word-boundary semantics, YAML authoring
grammar, semantic-pin tests) form a unit: applying filters at all
paths only matters if the patterns match; authoring patterns
correctly only matters if filters apply at all paths; and the
sandbox-vs-real semantics drift problem applies anywhere stubs
exist.  Future similar work — for instance, applying a new
`include_patterns` or per-slot file-type rules — should follow
the same four-step discipline.

### 35. Helper Consolidation Pattern (v119.H12)

**Context**: H10 introduced a closure inside `_find_file_for_slot`
that captured `input_def["exclude_patterns"]` and applied a
filter at 7 call sites within that function.  A parallel
`_select_files` call site (line 698) used `matches_exclude_pattern`
directly because the closure wasn't in scope.  Total: 9
v119.H10 markers, 8 application sites, 1 closure + 1 direct
call.

The result was logically correct but structurally noisy: same
intent at 8 sites with slightly different invocation forms.  H12
consolidates these into a single module-level helper.

#### 35.1 When closure → helper is the right consolidation

Three properties of the call-site population made closure → helper
the clean refactor:

- **Single derived input** (the exclude_patterns list) is the
  same at every site within `_find_file_for_slot`.
- **One parallel site outside** the closure's scope
  (`_select_files`) duplicated the logic with a different
  invocation form.
- **No per-site state** beyond the shared input — the helper
  is a pure function of `(file_path, exclude_patterns)`.

When ALL three hold, hoisting to a module-level helper:

1. Eliminates the scope asymmetry (every call site uses the
   same function call).
2. Lets the input be extracted ONCE at the top of the larger
   function (preserves the closure's micro-perf property).
3. Makes the helper unit-testable as a single function.

If any of the three didn't hold — e.g., if `_select_files`
needed a different filter shape, or if the helper depended on
mutable per-site state — closure → helper would force awkward
parameter threading and the closure would be the better
factoring.

**Pattern**: structural refactors of working code should hold
the test suite invariant.  H12's regression gate was Bug 10
+ Bug 11 tests (6 tests directly exercising the exclude_patterns
mechanism) all passing unchanged.  If those had failed,
behavior had drifted; if they passed, the refactor was safe.

#### 35.2 Helper signature: accept the raw input, not the wrapper dict

Gemini's H12 plan review observation: helpers that decouple
from the caller's schema are more reusable.  H12's
`_file_passes_exclude_patterns(file_path, exclude_patterns)`
takes the raw list, not the `input_def` dict that contains it.
This lets `_select_files`'s call site — which works with the
raw list extracted locally — pass it directly without first
re-wrapping into a dict.

**Pattern**: when a helper consumes one field from a larger
schema object, prefer the field type over the schema type as
the parameter.  Schema-typed parameters couple the helper to
the schema's evolution; field-typed parameters are stable.

#### 35.3 Marker convention: preserve history through refactor

H12's refactor sites are marked `# v119.H10/H12` rather than
just `# v119.H12`.  This preserves the audit trail: the
behavior comes from H10's bug fix; the structural form belongs
to H12.  A future reader doing `git blame` or searching for
v119.H10 markers still finds these sites; a future reader
searching for v119.H12 finds them too.

**Pattern**: when refactoring marked code, preserve the
original marker as a suffix or prefix, don't replace it.  Code
markers are commit-shadow metadata that compounds, not the
single most-recent marker.

#### 35.4 Semantic-pin tests for stub-vulnerable functions

The categorizer suite (`tst_file_categorizer.py`, 8 tests, 75
documented assertions) pins every public function in
`agent/file_utils.py` — `classify_mtz_type`, `get_mtz_stage`,
`get_category_for_extension`, `is_mtz_file`, `is_model_file`,
`is_map_file`, `is_sequence_file`.  Generalizes
`test_bug11_matches_exclude_pattern_semantics` from one
function to all of `file_utils.py`.

Why these functions specifically:

- Pre-H11, the sandbox stubbed `matches_exclude_pattern` with
  substring matching, while the real function uses word-boundary
  matching.  The stub-real semantic drift masked a production
  bug for one ship cycle.
- The other functions in `file_utils.py` have similar
  characteristics — small, non-trivial logic that's easy to
  stub badly and hard to verify just by reading.

**Pattern**: when a function has non-obvious semantics
(regex matching, parsing, classification with rule ordering,
documented quirks), it should have a semantic-pin test
that exercises the REAL function against a documented
fixture.  The test's value is not just regression detection;
it's documentation of "what the function actually does"
that's machine-verifiable.

Notable cases the categorizer suite pins:

- `classify_mtz_type` rule ordering — `refine_001.mtz` matches
  both Rule 1 and Rule 2; the test pins that Rule 1 fires first.
- `is_model_file`'s load-bearing quirk: `.pdb` extension
  overrides the `'ligand'`-in-name filter, so `ligand_fit.pdb`
  IS a model file (intended for content-based guards
  downstream).  A future contributor who "simplifies" the
  function would break this if no test pinned the contract.
- Multi-dot extensions, no-extension files (Gemini-suggested
  defensive cases).

**Pattern**: "non-obvious quirks" are the highest-value pinning
targets.  A test that asserts the obvious behavior is less
useful than one that asserts a documented exception.  When you
find yourself writing "BTW, this case is the exception to the
rule", that case should have a test.

### 36. Ollama Provider Robustness Patterns (v119.H13)

**Context**: Tom's `phenix.ai_agent provider=ollama ...` run on
cci-gpu-01 failed with `[DIRECTIVE_EXTRACTION_FAILED]: 404 page
not found`.  Diagnostic isolated two stacked bugs in the
Ollama-provider call sites of `agent/directive_extractor.py`
and `agent/api_client.py`, plus an opportunity to complete the
v118 §3.5 retired-model classification work.

#### 36.1 URL normalization as caller-side robustness

The OpenAI Python SDK appends `/chat/completions` to the
configured `base_url` to construct the actual request URL.
Ollama's OpenAI-compat endpoint lives at
`/v1/chat/completions`, so `base_url` must end in `/v1`.

Pre-H13, the code's default was `"http://localhost:11434/v1"`
but the env-var path
(`os.environ.get("OLLAMA_BASE_URL", default)`) had a silent
trap: if the user set `OLLAMA_BASE_URL=http://localhost:11434`
(without `/v1` — a natural value matching `OLLAMA_HOST`), the
env-var fully overrode the default, including the `/v1` suffix
the default supplied.  Result: request went to
`http://localhost:11434/chat/completions` → Ollama's literal
HTML `404 page not found`.

**Decision**: instead of forcing users to know the SDK's
behavior, the caller normalizes the URL.  The
`normalize_ollama_openai_base_url(base_url)` helper idempotently
ensures `/v1` suffix:

```python
def normalize_ollama_openai_base_url(base_url):
    if not base_url:
        return base_url  # passthrough None and empty string
    stripped = base_url.rstrip('/')
    if stripped.endswith('/v1'):
        return stripped
    return stripped + '/v1'
```

Properties:

- **Idempotent**: calling on an already-normalized URL returns
  the same URL (modulo trailing-slash stripping).
- **Empty-string and None passthrough**: caller's default-handling
  decides what to do.  Without this property, an empty
  `OLLAMA_BASE_URL` would produce `'/v1'` — a useless URL
  fragment that would fail downstream in a confusing way.
- **Single source of truth**: the same helper applies at both
  Ollama call sites, eliminating the silent-trap.

**Crucial NON-application**: `core/llm.py:get_llm_and_embeddings`
uses `ChatOllama` (langchain-ollama's native API), which hits
Ollama's native `/api/chat` endpoint and does NOT want `/v1`.
Applying the helper there would BREAK the native path.  An
explicit comment at that call site documents the asymmetry
and warns future refactors.

**Pattern**: when an environment variable becomes part of a
larger constructed value, the env-var alone is not enough
context to validate correctness.  The construction logic
should be a single helper that both the default and the
env-var path flow through.  This eliminates the silent-trap
where the default works but the env-var (which the user
believes is equivalent) doesn't.

#### 36.2 Centralizing env-var precedence rules

`core/llm.py:get_llm_and_embeddings` correctly honored
`OLLAMA_LLM_MODEL` ("env-var wins, central default falls back")
since v118.  But two other Ollama call sites
(`directive_extractor.py:1219` and `api_client.py:734`) called
`default_model_for_provider("ollama")` unconditionally —
ignoring the env-var.

H1 centralized the model-defaults TABLE.  H13 completes the
job by centralizing the PRECEDENCE RULE:

```python
_PROVIDER_MODEL_ENV_OVERRIDES = {
    "ollama": "OLLAMA_LLM_MODEL",
    # future providers add their env-var name here
}

def resolve_model_for_provider(provider, role="decision"):
    env_var = _PROVIDER_MODEL_ENV_OVERRIDES.get(
        (provider or "").lower().strip())
    if env_var:
        env_value = os.getenv(env_var)
        if env_value:
            return env_value
    return default_model_for_provider(provider, role=role)
```

Empty env-var (`OLLAMA_LLM_MODEL=""`) is treated as unset and
falls through to the central default — protecting operators
from the foot-gun where setting the env-var to `""` would
break their setup with no explanation.

**Pattern**: centralization happens in stages.  Stage 1: extract
the data table (H1 did this for `DEFAULT_MODELS`).  Stage 2:
extract the access patterns / rules over the table (H13 does
this for env-var precedence).  After Stage 2, every consumer
of the table is also a consumer of the rule.  If Stage 1
shipped but Stage 2 didn't, the table acquires inconsistent
access patterns over time as new consumers are added — exactly
the failure H13 found.

The pre-existing inline pattern in `get_llm_and_embeddings`
is left as-is rather than refactored to use the new helper.
Two reasons: (a) it uses `role="rag"` while the helper
defaults to `role="decision"`; (b) it has a sibling
`OLLAMA_EMBED_MODEL` env-var the helper doesn't know about.
Refactoring would require either threading the role arg AND
generalizing the env-var-name lookup, with no behavior change.
A comment at the call site documents the parity (same rule,
different inline implementation).

#### 36.3 Three-way 404 sub-classification

Pre-H13 the existing `_looks_like_retired_model_error` (from
v119.H1) lumped "model not found" together with retirement
phrases.  That worked for the original §3.5 use case (cloud
provider drops a model from its catalog) but misclassifies the
"local Ollama doesn't have this model pulled" case as
retirement.  The operator hint would say "update DEFAULT_MODELS"
— the wrong action.  The right action for the local case is
`ollama pull X`.

H13's `_classify_provider_error(exc, model_name)` returns one of
four marker tags:

| Tag | Trigger | Action |
|---|---|---|
| `MODEL_RETIRED` | 404 + retirement phrase near "model" word | Update DEFAULT_MODELS |
| `MODEL_UNAVAILABLE` | 404 + `model '...' not found` regex | Pull / rename |
| `AUTH_FAILED` | 401 / "unauthorized" / "invalid api key" | Check API key |
| `FAILED` | Anything else | Investigate |

The classification distinctions are operationally meaningful —
each tag maps to a different operator action.  The hint text
in each tag includes the specific action: "Update
`core/llm.DEFAULT_MODELS`", "Run `ollama pull X`", "Check the
provider API key".

**Pattern**: error categories should map to operator actions,
not just error types.  An error type ("404 not found") could
have multiple operator actions depending on context; the
classification should disambiguate.  When a single error type
maps to multiple actions, the classifier should split it.

#### 36.4 Defense-in-depth via class attribute + content matching

The classifier interrogates BOTH the exception's structured
properties AND its stringification:

```python
status_code = (
    getattr(exc, 'status_code', None)
    or getattr(getattr(exc, 'response', None), 'status_code', None))

err_body = ""
if hasattr(exc, 'body') and isinstance(getattr(exc, 'body'), dict):
    err_body = str(exc.body).lower()
elif hasattr(exc, 'message'):
    err_body = str(exc.message).lower()
err_all = (str(exc) + " " + err_body).lower()
```

Defends against two SDK behaviors:

- Some SDKs surface HTTP details on attributes only (`exc.status_code`,
  `exc.response.status_code`) while `__str__` returns an opaque
  class name or short code.  Pure string-matching would miss these.
- Other SDKs put the body in `__str__` directly but not in
  attributes.  Pure attribute interrogation would miss these.

Combined, the classifier handles both shapes.

**Pattern**: when consuming third-party exceptions whose internal
shape varies, redundantly interrogate both structured and
text-based forms.  The cost is small (a few lines) and the
robustness payoff is large (works across SDK versions and
provider variants).

#### 36.5 Co-occurrence rule against semantic false-positives

Naive matching ("404" + "deprecated") would false-positive on
edge-proxy 404 pages mentioning "this endpoint is deprecated;
see new path".  The endpoint is deprecated, not the model;
classifying as `MODEL_RETIRED` would noise the operator alert
channel.

H13's classifier requires the retirement phrase to appear
within 80 characters of the word "model":

```python
for phrase in _RETIREMENT_PHRASES:
    for m in re.finditer(re.escape(phrase), err_all):
        window = err_all[max(0, m.start() - 80) : m.end() + 80]
        if 'model' in window:
            return MODEL_RETIRED, ...
```

The 80-char window is wide enough to capture realistic
retirement messages ("The model X has been deprecated") and
narrow enough to reject co-mentions where the words happen to
appear in the same exception but not in the same semantic
context.

A load-bearing K-test pin
(`test_classify_edge_proxy_404_not_retired`) asserts the
edge-proxy case classifies as `FAILED`, not `RETIRED`.  If a
future refactor weakens the co-occurrence rule, this test
fires loudly.

**Pattern**: pattern-matching error classifiers benefit from
co-occurrence rules over single-phrase rules.  A single phrase
can appear in many contexts; a phrase-pair requirement adds
context-sensitivity for minimal complexity cost.

#### 36.6 H13 known limitations (carried as future work)

- **Other-provider classifier migration**: H13's classifier is
  used by the ollama exception handler in `_call_llm_fallback`
  and by the primary-path `extract_directives` catch.  The
  google/openai/anthropic per-provider catches in
  `_call_llm_fallback` still use H1's older
  `_looks_like_retired_model_error` helper, which doesn't
  distinguish UNAVAILABLE.  Migrating these three providers
  is a candidate H14 cleanup.
- **AUTH_FAILED has no dedicated emit helper**: when the
  classifier returns `AUTH_FAILED`, no marker is written to
  stderr (only the generic log line).  The hint is available
  but unused.  Pattern-consistency would add
  `_emit_auth_failed_marker`; deferred since the immediate
  bug was not auth-related.

### 37. Three Latent Regressions Surfaced by Batch Analysis (v119.H14)

**Context**: Tom's `run_39_openai` batch (520 runs, May 2026)
showed regressions vs `run_25_openai` (720 runs, March 2026):
phenix_sorry rate up from 16.9% to 28.1%, empty_summary rate up
from 24.2% to 36.0%.  `scan_batch_runs.py` (a new ops-side
triage tool that classifies run logs into Tier-1 crashes /
Tier-2 diagnostic markers / Tier-3 state anomalies / Tier-4 soft
anomalies) localized the regressions to three independent code
paths.  H14 fixes all three.

The investigation method itself is worth recording: a controlled
three-way log comparison of the 1029B-sad dataset across
`run_25-solve`, `run_39-solve`, and `run_39-stop` variants
isolated the trigger of the phaser bug to a single-line README
diff.  This kind of "same input minus one line" comparison is
the right diagnostic technique when a batch shows category-level
regressions — it converts a statistical signal into a
mechanistic one.

#### 37.1 Goal phrases vs method requests in `_ACTION_TABLE`

The rules-only directive extractor's `_ACTION_TABLE` maps
keyword matches in user advice to workflow actions.  Pre-H14
the `solve` entry's keyword list was:

```python
"solve": {
    "xray": "phenix.phaser",
    "keywords": ["molecular replacement", "solve the structure",
                 "solve structure", "phaser", "mr "],
}
```

This conflated two distinct semantic categories:

- **Method requests**: `"molecular replacement"`, `"phaser"`,
  `"mr "` — user explicitly asks for MR.
- **Goal phrases**: `"solve the structure"`, `"solve structure"`
  — user expresses a GOAL.  The right METHOD depends on the
  data type (MR for native + search model, autosol for SAD/MAD,
  predict_and_build for native + predicted model).

Including goal phrases in the action's keyword list meant the
extractor inferred "user wants MR" from any README mentioning
the goal.  When that goal phrase appeared alongside another
action keyword (e.g., "Stop after refinement" matching the
`refine` action), the multi-action branch in
`_apply_workflow_intent_fallback` set `start_with_program =
phenix.phaser` — forcing MR into the workflow.  On SAD/MAD
datasets, autosol then ran first anyway (workflow rules), but
the persistent directive meant phaser was retried later cycle
after cycle.

**Decision**: keyword lists in `_ACTION_TABLE` are restricted
to explicit METHOD requests.  Goal phrases (the user's intent
to solve, refine, validate, etc.) belong elsewhere — either in
the `classify_intent` machinery (which separately determines
the user's overall intent independently of method-specific
actions) or in workflow rules that key on data type.

**Property preserved**: when a method keyword IS explicitly
present (`"run phaser to place the model"`,
`"use molecular replacement"`), the `solve` action still
triggers and `start_with_program = phenix.phaser` is still set.
The cleanup removes only the false-positive on goal-only
phrasings.

This decision is small (5 chars removed from a keyword list)
but it pins a SEMANTIC invariant: **method actions match
method requests, not user goals.**  The H14 comment in the
`_ACTION_TABLE` entry records the invariant so future
keyword-list additions can be audited against it.

#### 37.2 Single-emit invariant for diagnostic_messages relay

Pre-v119.H5 markers like `[STEP_1F]` wrote directly to stderr
from inside the preprocessor.  v119.H5 §2.9 introduced the
`diagnostic_messages` channel: the preprocessor appends marker
strings to a list, the dispatcher relays the list to stderr at
the return path.  v119.H5.1 refactored further by moving the
relay INTO the dispatcher (`programs/ai_analysis.py`'s
`run_job_on_server_or_locally`) so all four return paths
(server local, server remote, client local, client remote)
share the same relay implementation.

The H5 → H5.1 refactor left a duplicate.  The pre-refactor
client-side relay in `programs/ai_agent.py` (lines 8087-8103
in the H13.1 ship) survived into H14.  Both relays fired on
every return, doubling every marker.

**Decision**: invariant is **single emit per marker**.  The
relay belongs in ONE place — the dispatcher — because that's
the single chokepoint that handles all dispatch modes
uniformly.  Client-side relay code is removed and a marker
comment records the H14 deletion so future readers don't
re-introduce it from old patterns.

The relay function's docstring already said "this is the
SINGLE, uniform site where operators see them."  H14 brings
the implementation in line with that claim — pre-H14 there
were two sites, and the docstring was wishful by ~one ship.

#### 37.3 Positive shape check for `space_group` values

The pre-H14 `validate_directives()` validator combined:

1. **Sentinel set** (`_SYMMETRY_SENTINELS`): explicit
   case-insensitive match against known placeholder phrases
   (`"not specified"`, `"not provided"`, etc.).
2. **Negative structural checks**: must start with a letter;
   length <= 25 characters.

These two checks together caught most cases but had a
predictable gap class: **prose phrases that happened to start
with a letter and fit in 25 chars** (e.g., `"Solve the
structure"`).  Such phrases passed both checks and reached
PHENIX as bogus PHIL params.

**Decision**: add a POSITIVE shape check
(`_looks_like_space_group`) in addition to the negative ones.
The pattern accepts standard Hermann-Mauguin space-group
symbols and rejects prose.

**Why both negative and positive checks**: the sentinel set is
the cheaper match for the common case (LLMs emit placeholders
verbatim).  The positive shape check is the SOUNDER guarantee
for the uncommon case (LLMs emit prose).  Together they are
both fast (frozenset lookup + single regex) and tight (only
real space-group symbols pass).

**The final regex** (after two review passes):

```python
_HM_FORM_RE = re.compile(
    r'^[PFICRHAB]'
    r'[0-9mcndabehr\s/_\-:]{0,24}'
    r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$',
    re.IGNORECASE
)
```

The alphabet carries crystallographic meaning at every position:

- `[PFICRHAB]` — the 8 Bravais lattice letters
- `0-9` — symmetry-element orders, screw subscripts, digits
- `m c n d a b` — mirror plane and four glide-plane types
- `e` — 2002-ITA double-glide (required for Aem2, Aea2, Cmce,
  Cmme, Ccce; without it 5 of 230 canonical symbols would be
  rejected)
- `h r` — axis-spec suffix letters for rhombohedral
  hexagonal/rhombohedral settings (R3:H, R3:R)
- `space / _ -` — separators (PERCEIVE-style spacing,
  monoclinic slash, inversion-center dash, mmCIF underscore)
- `:` — alternative cell/origin setting prefix (R3:H, P4/n:2,
  P21/c:b)

**Regex evolution — two review passes**:

The first draft pattern was too strict:

```python
r'^[PFICRHAB]\s*-?\d{0,3}(?:\s*[-/_]?\s*\d{0,3}){0,4}'
r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$'
```

It matched only 71/230 canonical symbols.  Self-review against
the full 230-symbol list (per Tom's "must admit all 240 space
groups correctly" instruction) caught this — the strict
digits-only inner pattern rejected all monoclinic slash forms
(P21/c, P21/n), all orthorhombic mirror/glide groups (Pmma,
Pnma, Pbca), all tetragonal mirror forms (P4mm, P4/mmm), and
all cubic high-symmetry groups (Pm-3m, Im-3m, Fd-3m, Ia-3d).
The fix was to switch from a digits-only inner pattern to a
character-class with the symmetry-indicator letters added.

Gemini then reviewed and identified that the self-review-fixed
regex still rejected alternative cell/origin settings:

- `R3:H` / `R3:R` — rhombohedral axis choice
- `P4/n:1` / `P4/n:2` — origin choice (centrosymmetric groups)
- `P21/c:b` — unique-axis cell choice (monoclinic)

These forms appear in real PDB/mmCIF metadata and cctbx tool
output.  The fix was to add `:` to the alphabet (Gemini's
explicit recommendation) plus `h` and `r` (Gemini's example
`R3:H` requires both — colon alone would not have been enough).

**Acknowledged residual limitation** (Gemini Risk B): the
permissive alphabet still admits short English words that
happen to use only HM-alphabet characters: `Panda`, `Fame`,
`Bad`, `Cab`, `Bed`, `Acme` all match.  These are not
realistic LLM outputs for `space_group`, and the pre-H14
negative checks already let them through.  Closing this gap
fully would require enumerating the 230 canonical symbols
(plus their alternative-setting variants) via
`cctbx.sgtbx.space_group_info` — a substantive cctbx
dependency beyond a directive sanity check's scope.  The
limitation is pinned with `test_hm_form_known_limitation_short_words`
so future contributors don't mistake "Panda accepts" for a
bug; the positive guarantee that words with non-alphabet
letters (Phaser, Pizza, Place) ARE rejected is pinned with
`test_hm_form_rejects_words_with_non_alphabet_chars`.

**Sentinel extension policy**: the additions to
`_SYMMETRY_SENTINELS` (the `"not explicitly *"` family plus
common truncated forms like `"not explicitly mentio"`)
recognize a class of LLM output where length-limited models
emit a partial sentinel phrase.  Truncated-form sentinels are
worth carrying because they're cheap (frozenset entries) and
they catch real production cases (Tom's qwen2.5:72b output).

**PHENIX fallback property**: dropping a malformed `space_group`
value is preferable to passing it through, since PHENIX can
usually auto-detect crystal symmetry from data when the
directive lacks an explicit value.  Pass-through of bogus
values produces obscure downstream errors;
`validate_directives` drops with a clear log line.

#### 37.4 H14 known limitations

- **Goal-phrase signal lost from extractor**: removing
  `"solve the structure"` from the `solve` action means the
  rules-only path no longer learns user GOAL from that
  phrasing.  That signal is still available to the
  `classify_intent` path (LLM-driven, used in non-rules-only
  modes).  For rules-only operation, the workflow engine's
  data-type rules suffice to pick the right method without
  the directive layer needing to encode the goal.
- **`asu_copies` → `phaser.search_copies` injection
  unmodified**: the `graph_nodes.py` injection that fires
  whenever phaser runs (using xtriage's Matthews `Best guess`)
  is correct domain behavior — `search_copies` should match
  the ASU copy count.  The 1029B-sad failure was upstream
  (Item 1 above): phaser ran when it shouldn't have, on a
  search model that was already the full ASU.  When phaser
  legitimately runs with a single-chain search model, the
  injection is right.  Defensive instrumentation for the
  multi-chain-search-model case (a `[BUILD]` warning if the
  model file appears multi-chain when `search_copies` > 1) is
  a candidate H15 item, NOT urgent.

#### 37.5 Dual-path validation closure (v119.H14.1)

H14 added the sentinel + Hermann-Mauguin shape check on
`space_group` inside `validate_directives`.  Tom's
2026-05-26 production verification of H14 with ollama showed
the fix did NOT close the bug in production — `space_group=Not
explicitly mentio` still reached the displayed directives even
though `scanner_version=119.H14` confirmed the H14 code was
installed.

**Root cause**: `directive_extractor.py` has TWO paths that
produce a directives dict:

```
                    extract_directives(advice, provider)
                                 │
                ┌────────────────┼────────────────┐
                ↓                ↓                ↓
       LLM returns valid    LLM returns      use_rules_only=True
       parseable JSON       empty/invalid    (bypass LLM)
                │                │                │
                │                ↓                ↓
                │      extract_directives_simple (regex path)
                │                │                │
                │                │                │
                └───→ validate_directives ←───────┘
                                 │              [pre-H14.1: simple
                                 ↓               path returned
                              return              WITHOUT validation]
```

Pre-H14.1, the `extract_directives_simple` path returned its
result directly to the agent.  This worked for the LLM-success
case (where the agent's `extract_directives` wrapped the call
with `validate_directives`), but when the LLM failed and the
fallback was returned, OR when `use_rules_only=True` was set,
validation was skipped.

The simple extractor has its OWN space_group regex (at line
~4350 in directive_extractor.py):

```python
r'space[_ ]group(?:\s*[=:]\s*|\s+(?:is|of|=|:)?\s*)([A-Za-z][A-Za-z0-9 /_-]{1,20})'
```

This pattern is **permissive by design** — it accepts up to
20 chars of `[A-Za-z0-9 /_-]` after a leading letter, so
"Not explicitly mentioned" → captures "Not explicitly mentio"
(21 chars).  Pre-H14.1, the bogus value passed straight
through to the agent.

**Architectural decision in H14.1**: rather than duplicate
the H14 sentinel/shape check inline at the simple extractor's
regex site, make `validate_directives` the canonical
final-sanity step that BOTH paths converge through.  The
simple extractor now ends with:

```python
directives = validate_directives(directives, log)
return directives
```

This means future validator extensions automatically apply to
both paths — no risk of one path drifting from the other.
The pattern is "extract-then-validate" instead of
"validate-inside-each-extractor".

**Prerequisite fix**: `start_with_program` had to be added to
`VALID_STOP_CONDITIONS`.  This was a pre-existing latent bug
— the key was set by `_resolve_after_program` and consumed by
`workflow_engine.py`, `ai_agent.py`, and `ai_analysis.py`,
but missing from the validator's allow-list.  Pre-H14.1 the
gap was invisible because the LLM path's call order was
`extract → validate → fallback-overlay`, so the resolver
added `start_with_program` AFTER validate ran.  In the simple
extractor, the resolver runs INSIDE the extractor (before
return), so the post-H14.1 final-validate would have stripped
the key without this fix — silently breaking H14 Item 1's
fix for the rules-only path.

**Lesson generalized**: when a value's validator and its
extraction sites are decoupled, the canonical pattern is
"all extractors converge through a single validator at the
extract-function boundary" — never "validator runs in one
extractor but not others."  Adding a new extractor adds a
single line (call the validator before return); adding a new
validator rule changes one site (the validator function);
neither requires touching the other.

### 38. Targeted Resume-Reopen for Plan Stages (v119.H15)

**Rationale**: Tom's bromodomain resume failure (run 144)
exposed a three-part interaction.  Bug 1 corrupted
`final_refinement.status` to COMPLETE.  Resume cleared
`gate_stop` but NOT per-stage statuses.  The gate immediately
re-fired "all stages complete" → LLM saw `STATE=complete` →
chose `phenix.polder` against the user's resume directive.

H15 Item 1 prevents the corruption going forward.  Item 2 is
the resume safety net implemented in
`agent/plan_generator.py::reopen_stages_for_directives`.

**Targeted single-stage reopen** (per Gemini's critique of the
original blast-radius proposal):

1. Walk the directives looking for `program_settings.<program>`
   entries.
2. For each program named, find the LATEST completed stage
   whose `programs` list contains that program.
3. Reset ONLY that stage to PENDING.  Reset `cycles_used → 0`
   and clear `started_at`, `completed_at`, `last_result`.
4. Do NOT cascade-reset stages downstream of the reopened one.
5. Earlier stages with the same program are NOT reopened.
6. Skipped stages stay skipped (the user can't un-skip a
   stage by mentioning its program; that requires explicit
   directives).

The "LATEST stage with matching program" rule keeps the blast
radius O(1) regardless of plan size — exactly one stage moves
to PENDING per matched program, even on a 12-stage cryo-EM
pipeline.

**Why latest-stage, not earliest-stage**: when the user
mentions a program in resume advice, the strongest semantic
interpretation is "I want to re-do this program in its most
recent context."  Reopening an earlier stage would discard
the work that was done after it.

**K_H15_ITEM_2** (7 tests in `tests/tst_resume_reopen_stages.py`):
§A Tom's exact scenario; §B no advice change → no-op;
§C directive for unmatched program → no reopen; §D multiple
programs in directives → still O(1); §E skipped stages stay
skipped; §F empty plan / empty directives → no-op, no
exception; §G stage's strategy already honors directive —
H15 design says reopen anyway (semantic: "user re-asserted
intent, so respect that").

**Why important for future work**: §3.3 deferred work item
("infer experiment type from input file extensions at session
creation time") would complement H15 by ensuring a resumed
session's experiment_type lock survives across the gate
re-fire.  H15 fixes the stage-state side; §3.3 would fix
the experiment-type side.

### 38.5 Reactive-Deviation Hold for Un-Run Must-Run Stages (v120)

**Rationale**: H15 (§38) made the plan-deviation catch-up *honest*
(intermediate stages that don't meet criteria are marked FAILED, not
silently COMPLETE), but it did not stop the catch-up from abandoning a
stage that *never ran at all*.  Tom's beta-blip rebuild (AIAgent_21)
exposed the gap: the plan selected `refine_rebuild_placed`, so
`model_rebuilding` (programs `[phenix.autobuild, phenix.refine]`) was the
active stage after refinement.  At cycle 3 the LLM chose
`phenix.molprobity` — reasoning, sensibly, "validate first, then rebuild."
But molprobity matches a LATER stage (validation), so
`record_stage_cycle`'s catch-up advanced PAST `model_rebuilding` (recording
it COMPLETE with `cycles_used == 0`) to validation.  The plan then read
"all stages complete" → STOP.  Autobuild never ran; the user's explicit
rebuild request was lost to a single reactive validation.

This is the architectural "tail wagging the dog" the code comments flag:
the LLM's reactive program choice drives the plan, and choosing a
later-stage program collapses the stages in between.

**The fix (Option 2a — hold-and-sustain, not force-forward).**  Three
coordinated pieces let the LLM take ONE reactive detour without abandoning
an un-run, must-run stage:

1. *`knowledge/plan_schema.py` — the hold.*  `StageDef` gains a
   `reactive_deviations` counter (serialized; tolerant default 0 for
   pre-existing sessions).  In `record_stage_cycle`, when the program
   matches a LATER stage but the CURRENT stage is *must-run-and-un-run*
   (`cycles_used == 0`, its lead program `programs[0]` has not completed in
   history, its success_criteria are not met) and `reactive_deviations < 1`,
   the catch-up is HELD: the stage stays active, no cycle is counted (lead
   never ran), `reactive_deviations` is incremented, and a HELD deviation
   event (`from == to == curr.id`) is recorded for observability.  A log
   line makes it explicit:
   `PLAN_GUARD: Holding model_rebuilding stage while side-cycle (molprobity)
   executes (reactive_deviations: 1/1)`.  On a SECOND deviation the guard is
   exhausted and the normal H15 catch-up proceeds — so there is no infinite
   hold.  The "lead already ran" test counts a variant
   (`phenix.autobuild_denmod`) as the lead (`phenix.autobuild`); it
   deliberately does NOT reuse `_program_matches_phase` for this, because
   that helper returns True on an empty program list and would falsely
   report the lead as run when history is empty.

2. *`programs/ai_agent.py` — surface the lead program.*  A new
   `_get_plan_current_unrun_lead_program(session)` returns the current
   ACTIVE stage's lead program when it is un-run.  This is necessary because
   `_get_plan_next_stage_programs` only inspects `status == "pending"` and so
   misses the active held stage.  The value is passed into `session_info` as
   `plan_current_unrun_lead_program`.

3. *`agent/graph_nodes.py` — offer it in PERCEIVE.*  The existing
   plan-stage injection only fired when `valid_programs` was STOP-only.  A
   new injection adds `plan_current_unrun_lead_program` to `valid_programs`
   *even when other non-STOP programs are already valid* — otherwise, after
   the hold, the engine keeps offering only molprobity and autobuild is
   never presented.  Logged as
   `PERCEIVE: Injected plan lead program <prog> (current stage active and
   un-run)`.

**End-to-end**: cycle 3 molprobity → HELD (stage stays active) → gate
evaluates `model_rebuilding` (skip_if/​success unmet, not exhausted →
continue) → cycle 4 autobuild offered and run → final refine → validation.
The LLM gets its one "look before you leap"; the rebuild still happens.

**Deployment coupling**: the hold's criteria check (`_criteria_met`) lazily
imports `GateEvaluator`.  If that import fails at runtime it returns True
(criteria "met") and the hold does NOT fire — a safe degradation toward
legacy catch-up, but it means `gate_evaluator.py` must be deployed
alongside `plan_schema.py` for the hold to engage.  `graph_nodes.py` is
server-side (PERCEIVE), `ai_agent.py` client-side, `plan_schema.py` shared —
deploy all three together.

**Transport plumbing (parity)**: the `plan_current_unrun_lead_program` field
this section relies on must round-trip the full client→server path
(`build_session_state` → `build_request_v2` wire whitelist → `run_ai_agent`
map-back), or PERCEIVE never sees it and the injection is inert — identically for
local and server.  v120 plumbs it (and the previously-dropped
`plan_has_pending_stages` / `plan_next_stage_programs`); see the "`session_info`
field plumbing contract" under Request Building above.  So the full deploy set for
this feature is `plan_schema.py`, `gate_evaluator.py`, `ai_agent.py`,
`graph_nodes.py`, `contract.py`, `api_client.py`, and `run_ai_agent.py` — all to
BOTH server and Mac.

**Tests**: `tests/tst_reactive_deviation_hold.py` (7, real StructurePlan +
real GateEvaluator: 1st deviation held; 2nd catches up; no hold when lead
ran / criteria met / variant ran; holds with empty history;
serialization round-trip incl. missing-field default) and
`tests/tst_plan_lead_program_offer.py` (8, lead surfacing + PERCEIVE
injection).

**Status**: this is the most invasive plan/LLM-arbitration change to date.
The narrow trigger (un-run + lead-not-run + criteria-unmet + once) and the
safe `_criteria_met` degradation bound the risk, but the next few rebuild
runs should be watched.  Option 1 (broaden the next-stage-program injection
to steer the LLM toward the plan's program directly) remains a documented
fallback if the hold proves too permissive or too strict.

### 39. MTZ Auto-Inspection and obs_labels Injection (v119.H16 / H16.1)

**Rationale**: 88 TIER-1 failures across two batch scans
(run_25, run_39) matched "Sorry: Multiple equally suitable
arrays of observed xray data found", concentrated in
AF_exoV_MRSAD and lysozyme-MRSAD tutorials.  The MTZ files
contained multiple legitimate observation arrays (merged
Iobs AND anomalous I(+)/I(-)) and PHENIX programs error out
unless the user disambiguates via
`scaling.input.xray_data.obs_labels=` (or equivalent
per-program PHIL path).  Pre-H16, the LLM had to remember
to inject this parameter manually for every multi-array
MTZ; when it forgot, the program crashed and the agent
treated it as an unrecoverable error.

H16 adds a three-layer auto-inspection mechanism:

**Layer 1 — `agent/mtz_inspector.py`** (NEW module):

```
inspect_mtz(path) -> {
    "columns": [...],      # raw column list from iotbx.mtz
    "anomalous_intensity_pairs": [("I(+)", "I(-)", "SIGI(+)", "SIGI(-)")],
    "merged_intensity": "IMEAN",
    "merged_amplitude": "FOBS",
    "rfree_flag": "R-free-flags",
    "phib_fom_pair": None,  # populated if columns present
    "error": None,
}
```

`inspect_mtz()` reads via `iotbx.mtz` (cctbx); side-effect-free;
never raises (returns `{"error": ...}` on file-read failures
so the caller can decide).

**Layer 2 — `_PROGRAM_PREFERENCES` table**:

```python
_PROGRAM_PREFERENCES = {
    "phenix.xtriage":          {"prefer": "merged_intensity",
                                "fallback": "merged_amplitude"},
    "phenix.autosol":          {"prefer": "anomalous_intensity",
                                "fallback": "merged_intensity"},
    "phenix.phaser":           {"prefer": "merged_intensity"},
    "phenix.predict_and_build":{"prefer": "merged_intensity"},
}
```

`select_obs_labels_for(program, mtz_info)` reads the table
and returns the program-specific obs_labels string to inject
(or None when no ambiguity exists).

**Layer 3 — YAML invariant + builder hook**: `programs.yaml`
declares `auto_fill_obs_labels: true` on xtriage, autosol,
phaser, and predict_and_build.
`command_builder._apply_invariants()` reads the flag, calls
`inspect_mtz()` on the relevant `data_mtz` input, calls
`select_obs_labels_for()`, and injects the result into the
command before execution.

**Why a separate inspector module, not a method on
CommandBuilder**: `inspect_mtz()` is a side-effect-free
content-reading function that should be unit-testable
without instantiating the builder.  Layer separation also
allows other consumers (e.g., the recovery system) to
inspect MTZs without going through the builder.

**Why per-program preferences**: phaser does MR on merged
intensities; autosol does SAD on anomalous intensities.
Picking the same default for both would force the user to
override one program or the other.  The table records the
canonical choice per program.

**Three-tier test layering** (8 tests in `tst_obs_labels_auto_fill.py`):
- §A: Policy unit tests (pure functions, no I/O — always run)
- §B: MTZ inspector tests (cctbx-dependent, gracefully skip
  in sandbox)
- §C: Builder-integration simulation (mirrors the actual
  `_apply_invariants()` branch without requiring a full
  CommandBuilder instance)
- §D: YAML config validation

H16.1 follow-up: bumped the scanner version pin in
`[STEP_1F]` telemetry to `119.H16.1` so production logs
clearly indicate the obs_labels auto-fill is active.

**Why important for future work**: the inspector + per-program
preference table pattern generalizes.  Other "discovery before
execution" needs (e.g., space-group detection for ambiguous
crystals) could use the same shape.

### 40. Reactive Recovery with `strip_parameter` Resolution (v119.H17 / H17.1)

**Rationale**: Tom's lysozyme-MRSAD cycle 5.  The LLM
correctly identified MR-SAD and ran xtriage → phaser → autosol
successfully (Bayes CC 74.20).  On cycle 5 it called
`phenix.autobuild` with the user-supplied raw anomalous MTZ
(`lyso2001_scala1.mtz`) as `map_file=`.  autobuild requires
`PHIB` phase columns in `input_map_file`/`map_file`;
the program crashed with "Sorry, PHIB is required for
input_map_file"; the agent had no recovery template; the
workflow halted after four retries.

**Why reactive, not proactive**: user-supplied
`map_coeffs.mtz` IS a legitimate input in many workflows
(any post-phasing build).  Pre-emptively rejecting `map_file=`
on all autobuild calls would break the legitimate case.  H17
detects the specific PHIB-required error AFTER it fires and
strips the offending parameter on retry.

**The `strip_parameter` resolution kind** (introduced in H17):
alongside the existing `add_parameter` (adds a flag) and
`select_value` (disambiguates an enum), H17 introduces
`strip_parameter` for the case where the right fix is to
REMOVE an inappropriate flag entirely.  The YAML schema:

```yaml
missing_phib_input_map_file:
  patterns:
    - regex: "PHIB is required for input_map_file"
    - regex: "PHIB.*required.*input_map_file"
  resolution: strip_parameter
  strip_parameters:
    - map_file
    - input_map_file
    - input_files.map_file
  max_retries: 1
  description: "Strip user-supplied map without phases"
```

The `strip_parameters` list names the PHIL paths to strip
(three aliases for the same input slot in this case).

**Strict conjunction detection**: the YAML pattern requires
BOTH the literal phrase "PHIB is required for input_map_file"
AND the column-spec context.  Without the conjunction,
generic "PHIB" mentions in unrelated error messages would
false-positive.

**Two-phase implementation**:

H17 ships the analyzer side: YAML entry, new `strip_flags`
field on `ErrorRecovery`, generic `_resolve_strip_parameter`
handler.

H17.1 ships the executor side (closes Scenario B surfaced
in H17 pre-deploy review).  Three edits to `programs/ai_agent.py`:

1. **`_handle_recovery`** stashes `strip_flags` into
   `session.data["pending_strip_recoveries"]` keyed by program
   name.  Preserves `set_recovery_strategy` for add_parameter
   recoveries.
2. **`_execute_command`** pops any pending entry for the
   current program at the top of the method and applies a
   robust regex to strip each flag prefix.  Emits `[STRIP]`
   log line.  One-shot via `pop()` — entry is consumed on
   use.
3. **`_print_recovery_notice`** shows "Action: Stripping
   [...]" instead of the awkward "Selecting ''" for strip
   recoveries.

**The robust strip regex** (per Gemini's critique on the
H17 plan):

```python
pattern = (r'(?:^|\s)' + re.escape(flag_prefix)
           + r'\s*=\s*(?:"[^"]*"|\'[^\']*\'|\S+)')
```

Handles all the legitimate forms a flag value can take:

```
flag=value          (basic)
flag = value        (PHIL spacing — spaces around =)
flag="path spaces"  (double-quoted with spaces)
flag='path'         (single-quoted)
flag at end of line (no trailing whitespace)
```

A naive `\S+` pattern would corrupt quoted-with-spaces and
miss PHIL-spacing.

**Retroactive side benefit**: `rfree_flags_mismatch` had been
declared in `recoverable_errors.yaml` for months but never
wired in the analyzer (it would never trigger).  H17's
generic `_resolve_strip_parameter` handler retroactively
resolves that error type too.

**End-to-end validation on lysozyme-MRSAD**:
- Cycle 1: xtriage `ambiguous_data_labels` recovery (pre-H17
  mechanism still works)
- Cycle 2: xtriage SUCCESS
- Cycles 3-4: phaser SUCCESS, autosol SUCCESS (Bayes CC 74.20)
- Cycle 5: autobuild FAILED with PHIB error → `[NOTICE]
  DETECTED RECOVERABLE ERROR` with "Action: Stripping
  [map_file, input_map_file, input_files.map_file]"
- Cycle 6: `[STRIP] phenix.autobuild: removed
  'map_file=...' from retry command (recovery:
  missing_phib_input_map_file)` followed by `Running:` with
  no map_file= in command; autobuild SUCCESS.

**K_H17_strip_executor** (9 tests in
`tests/tst_h17_strip_executor.py`) pins the regex pattern.
If the pattern is changed in `ai_agent.py`, the test file's
reference implementation must be updated in lockstep.

**Why important for future work**: the strip_parameter
resolution kind closes one of three reactive recovery
"shapes" (add, select, strip).  The natural fourth shape
is "replace_parameter" (e.g., when a deprecated flag has
a new spelling); not yet motivated by production but easy
to add.

### 41. File-Based Experiment-Type Detection as Primary Signal (v119.H18)

**Rationale**: Tom's AF_7mjs density-modify-and-stop
regression (distinct from §20's bug, surfaced after H17.1
deployed).  User wrote literally "density modify and stop"
with cryo-EM half-map inputs (`.ccp4` files, no `.mtz`).
The extractor produced
`after_program=phenix.autobuild_denmod` (the X-ray
density-modification program).  §20's correction was
supposed to map this to `phenix.resolve_cryo_em` for cryo-EM
data — but its text-only `_detect_experiment_type_signals`
returned None ("ambiguous, decline to act") because the
terse user text "density modify and stop" has neither
cryo-EM nor X-ray tokens.  The LLM downstream then saw the
wrong `after_program` in VALID PROGRAMS, recognized it
couldn't actually run (cryo-EM has no MTZ), and selected
`phenix.predict_and_build` as "the next logical step" —
overriding the user's explicit "stop" instruction.

**Root cause**: §20's correction relied on text-based
detection alone.  For terse advice, text gives no signal,
and the file inventory (the unambiguous evidence) was never
inspected.  §17 documents that `session.set_experiment_type()`
only locks AFTER the first program returns — so at
directive-extraction time, the locked experiment_type is
None and the validator had to infer one.  §3.3 deferred work
item proposes inferring experiment type from input file
extensions at session creation time; H18 is the smallest
first step in that direction.

**New public helper** (`agent/file_utils.py`):

```python
def infer_experiment_type_from_files(files):
    """Returns (type, evidence_dict).
    
    type:    "xray" | "cryoem" | None
    evidence:
        xray_exts:   sorted unique extensions seen (.mtz, .sca, .hkl)
        cryoem_exts: sorted unique extensions seen (.mrc, .ccp4, .map)
        is_mixed:    True iff both types are present
    """
```

Asymmetric semantics match the existing
`_detect_experiment_type_signals` policy: returns
"xray" only if X-ray extensions are present and no cryo-EM
extensions; returns "cryoem" only the reverse; returns None
for mixed or empty inputs.

**Detection priority (locked policy)** — applied identically
in BOTH `_apply_experiment_type_program_reprints` AND
`_resolve_after_program` (the v115.10 post-LLM overlay):

```python
target_type = None
if original_files:
    file_type, evidence = infer_experiment_type_from_files(original_files)
    target_type = file_type

if target_type is None:
    text_type = _detect_experiment_type_signals(combined_advice)
    target_type = text_type

if target_type is None:
    return directives  # decline to act
```

**Files-win on conflict** (policy per Gemini's H18 review):
file extensions are the hard physical boundary (passing
`.ccp4` to an X-ray-only program crashes at parse time
regardless of what the user wrote); text-based detection has
known false-positive vectors ("mad" inside "modify", "sad"
in arbitrary sentences).

**Two correction sites needed the fix**.  H18 rev 1
implemented files-first only in §20's site
(`_apply_experiment_type_program_reprints`).  During the
merge into the existing `tst_density_modify_experiment_type.py`
suite, K14 (the exact AF_7mjs failing case) revealed a
second silent-revert site: `_resolve_after_program` (called
from `_apply_workflow_intent_fallback`, the v115.09 post-LLM
overlay) had its own text-only experiment-type heuristic
that defaulted `_exp="xray"` for terse advice and
unconditionally overrode `after_program` via the
`_ACTION_TABLE` lookup.  H18 rev 2 fixes both sites
identically using the same helper and policy.

**Telemetry (Pitfall 2 — Silent Override Hazard)**: the
enriched `[DIRECTIVE_CORRECTION]` marker records
`source=files|text`, `evidence=['.ccp4', ...]`, and
`text_signal=xray, OVERRIDDEN` (when files override
contrary text).  Example:

```
[DIRECTIVE_CORRECTION] Mapped after_program=phenix.autobuild_denmod
  to phenix.resolve_cryo_em (source=files, target_type=cryoem,
  evidence=['.ccp4'], text_signal=xray, OVERRIDDEN)
```

The `_corrected_from` sidecar on the directives dict gains
matching fields (`source`, `evidence`, `text_signal_overridden`)
for downstream traceability.

**Telemetry (Pitfall 1 — Dirty Directory Poisoning)**:
long-running sessions that MERGE files across cycles (per
`Session.set_project_info`) could accumulate both `.mtz`
and `.ccp4`, pushing detection into the "mixed → None"
state and disabling the correction.  A new
`[DIRECTIVE_CORRECTION_MIXED]` log marker fires when both
types are detected, making accumulated-files drift visible
in audit logs:

```
[DIRECTIVE_CORRECTION_MIXED] Both X-ray (['.mtz']) and
  cryo-EM (['.ccp4']) files present in inventory; deferring
  to text-based detection (text=cryoem)
```

In practice this is a non-issue for the AF_7mjs bug path
because directive extraction is gated by
`directives_extracted=True` and only fires on the initial
extraction, before any cycle has run.  But the telemetry
is defensive.

**Single source of truth**: `agent/plan_generator._build_context`
was refactored to call the new shared helper.  Previously
it had its own private file-extension mapping at
lines 232-236; H18 routes it through
`infer_experiment_type_from_files()` so plan_generator and
directive_extractor cannot drift apart on detection
behavior.  This satisfies §3.3's deferred work direction
in a constrained, surgical way (no session-locking changes).

**PHIL transport for `original_files`**: H18 introduces a
new PHIL parameter `original_files_for_directives` in
`programs/ai_analysis.py`.  The encoder side in
`programs/ai_agent.py::_extract_directives` joins basenames
with commas (avoids the trailing-whitespace issue in full
paths); the decoder side splits and feeds the list to
`run_directive_extraction()` which threads it into
`extract_directives()` and onward to the validators.
Basename-with-comma is a theoretical edge case but PHENIX
conventions don't use comma in filenames; if encountered,
the inference still produces correct extensions (the
comma-fragment without extension is silently ignored).

**Backward compatibility**: `original_files` is optional
with default `None` at every signature change.  Pre-H18
callsites without files get exact pre-H18 behavior —
text-only detection.  The 13 v118.9 §20 tests
(K1-K13) continue to pass unchanged.

**K-test surface**: 7 new tests (K14-K20) extend
`tst_density_modify_experiment_type.py` to a 20-test suite:
- K14: AF_7mjs failure verbatim — terse advice + cryo-EM
  files → corrected via files (source=files)
- K15: Mirror — terse advice + X-ray files
- K16: Backward compat — no files → text fallback still
  corrects (source=text)
- K17: Pre-H18 bug path — no files + terse text → declines
- K18: Files-win on conflict + OVERRIDDEN telemetry
- K19: Mixed input → defers to text + `[..._MIXED]` marker
- K20: `plan_generator._build_context` uses shared helper

**Why important for future work**: H18 implements §3.3's
first step.  The full §3.3 implementation would lock
`session.experiment_type` at session creation from files,
making the correction in `_apply_experiment_type_program_reprints`
unnecessary in the common case.  H18's helper is the
load-bearing piece for that future change; the full
implementation just adds a `session.set_experiment_type()`
call at session-creation time using the same helper.

**Lesson**: when fixing a categorical bug, grep for every
site that does the same kind of inference.  H18 rev 1
fixed `_apply_experiment_type_program_reprints` but the
same flawed text-only logic was re-applying the wrong
default in `_resolve_after_program` downstream.  Without
that second-site fix, the correction was silently reverted
and the AF_7mjs bug persisted.  K14 caught this at merge
time precisely because it tested the production entry
point, not the validator in isolation.

### 42. Dual `master_params` Schemas and PHIL Round-Trip Defense (v119.H18.1)

**Rationale**: H18 deployed clean with all 20 K-tests passing,
but Tom's production re-run of AF_7mjs still ran
`predict_and_build` on cycle 3 instead of stopping after
`resolve_cryo_em`.  The log showed:

```
DIRECTIVES: Extraction failed - Assignment to non-existing
  attribute "ai_analysis.original_files_for_directives"
  File ".../programs/ai_agent.py", line 8513, in _extract_directives
    directive_params.ai_analysis.original_files_for_directives = (
```

H18 had added the PHIL declaration for the new parameter to
`programs/ai_analysis.py` (line 165) but NOT to
`programs/ai_agent.py`'s OWN `master_params` string at
line 144.  The crash was caught by the surrounding try/except
in `_extract_directives` and silently swallowed; directive
extraction returned `{}`; the agent fell back to the default
cryo-EM plan template (mtriage → resolve_cryo_em →
predict_and_build).

**The dual-master_params architecture**.  This codebase has
TWO independent PHIL schemas declared in two different files:

```
programs/ai_agent.py:144           programs/ai_analysis.py:~120
        ↓                                   ↓
   master_params (str)                master_params (str)
        ↓                                   ↓
   libtbx.phil.parse(...)             libtbx.phil.parse(...)
        ↓                                   ↓
   master_phil object                 master_phil object
        ↓                                   ↓
   self.params (in AI agent)          params (in analysis server)
```

These two schemas are NOT transitive.  Adding a parameter to
`ai_analysis.py`'s master_params does NOT add it to
`ai_agent.py`'s schema.  The schemas overlap because both
declare an `ai_analysis` scope, but each declaration is a
separate object.

**Why H18's K-tests didn't catch this**.  The 20 K-tests in
`tst_density_modify_experiment_type.py` call helpers DIRECTLY
(`infer_experiment_type_from_files`,
`_apply_experiment_type_program_reprints`,
`_resolve_after_program`) with Python dicts as inputs.  None
of them go through PHIL parsing → deep-copy → attribute
assignment.  The bug was entirely in the PHIL layer, which
the tests bypassed.

**The fix**: add the missing PHIL declaration to
`programs/ai_agent.py`'s `master_params` string after
`user_advice_raw`, mirroring the declaration already shipped
in `programs/ai_analysis.py`.  14 lines added; no code logic
changed.

**Preventive K-test pattern**:
`tests/tst_h18_1_phil_roundtrip.py` exercises the full
production code path:

1. **Source-grep verification** (sandbox-safe): extract the
   `master_params = """..."""` region from `ai_agent.py` and
   confirm the new parameter name is present.  Runs without
   libtbx, fast, deterministic.
2. **PHIL parse + extract** (PHENIX-only, skips in sandbox):
   call `libtbx.phil.parse()` on the extracted master_params
   string and confirm `extract()` produces an object with the
   new attribute.
3. **Assignment reproduction** (PHENIX-only): deep-copy the
   extracted params and perform the exact assignment that
   crashed in production.  Must succeed.  Verify the value
   round-trips.
4. **Server-side definition present**: confirm
   `programs/ai_analysis.py` still has the matching PHIL
   declaration (don't fix one side and regress the other).
5. **Assignment site anchor**: confirm the assignment is
   still at the expected location.
6. **Cross-file consistency**: same parameter name in both
   files.

The source-grep variant (test 1) is the load-bearing one —
it runs in any environment, no PHENIX required, and catches
the entire class of "I forgot to declare the param in
master_params" bugs at sandbox time.  Verified by
temporarily reverting the H18.1 PHIL block: test 1 failed
with a clear diagnostic message identifying the deploy gap.
Restored — all 6 pass.

**Why important for future work**: the dual-master_params
pattern is permanent (the analysis server runs separately
from the agent and has different PHIL needs).  Every future
client→server PHIL parameter addition must declare in BOTH
files.  The `tst_*_phil_roundtrip.py` K-test pattern should
become the standard for any such addition.

**Generalize**: the source-grep K-test pattern is a defense
against any "schema declaration missing on one side" class
of bug.  It applies to TypedDict subclasses (`AgentState`),
dataclasses with explicit fields (`ErrorRecovery`), and any
other schema-driven type.  When adding a new field to a
schema, add the smallest possible test that grep-checks
each producer file for the new field name.  Fast, no
dependencies, no false positives.

### 43. Three Callsites of `_resolve_after_program` Require Lockstep Updates (v119.H18.2)

**Rationale**: H18 added an optional `original_files` parameter
to `_resolve_after_program()` for files-first experiment-type
detection.  The H18 audit threaded `original_files` through two
callsites — both billed in their surrounding comments as
"experiment-type detection" sites.  A THIRD callsite existed
but was billed differently and slipped through.

Tom re-ran AF_7mjs after H18.1 deployed.  The runtime tracer
revealed: the LLM emitted `{"stop_conditions":
{"stop_after_requested": true}}` with NO `after_program` field.
H18 site 1 (`_apply_experiment_type_program_reprints`)
correctly took its "no after_prog → early return" branch.  But
between H18 site 1 and H18 site 2
(`_apply_workflow_intent_fallback`), a v117.2 fallback at
`directive_extractor.py:783` fired:

```python
# Pre-H18.2: this is the v117.2 fallback path
# (added in v117.2 — predates H18).
_v172_source = raw_advice if raw_advice else user_advice
if _v172_source and _is_stop_after_requested(_v172_source):
    _resolve_after_program(directives, _v172_source.lower())
    #                                                       ↑
    #                                       MISSING: original_files
```

This callsite is in the post-LLM validation pipeline.  Its
purpose: when the LLM signals user-stop intent
(`stop_after_requested=True`) but omits the program field, parse
the raw advice to fill in `after_program`.  Without
`original_files`, the resolver defaulted `_exp="xray"` via the
text-only fallback heuristic (raw advice "density modify and
stop" has no cryo-EM/X-ray tokens) and mapped `denmod` →
`phenix.autobuild_denmod` — overwriting the correct empty state.

The downstream `_apply_workflow_intent_fallback` (H18 site 2)
DID pass `original_files` to `_resolve_after_program`, so the
files-first detection would have produced the right answer.
But by the time it ran, the preprocessed advice contains "Stop
Condition: None", so `_is_stop_after_requested(advice)` returned
False.  The resolver entered the `n==1, no stop → leave as-is`
branch instead of the `n==1, has_stop → override` branch.  The
buggy `after_program` from the v117.2 path persisted.

**Three-callsite architecture (post-H18.2)**:

```
extract_directives()
├── _validate_after_program_grounded()
│
├── _apply_experiment_type_program_reprints()   ← H18 site 1
│       └─ uses original_files (files-first detection)
│
├── v117.2 fallback (lines 783–798)             ← H18.2 site
│       └─ calls _resolve_after_program(..., original_files=original_files)
│
└── _apply_workflow_intent_fallback()           ← H18 site 2
        └─ calls _resolve_after_program(..., original_files=original_files)
```

All three sites now use files-first detection with text as
fallback.

**Why the v117.2 fallback exists**: the LLM extraction occasionally
emits `stop_after_requested=True` without `after_program` (Tom's
production AF_7mjs run is the canonical example).  Without
v117.2's fill-in step, the stop-check downstream is unreachable
(it's gated on `if after_program:`), and the agent runs the plan
template through to completion.  v117.2 closes that gap by
parsing the raw user advice deterministically.

**Why H18 missed it**: the v117.2 callsite is in a code region
billed as "fill in missing field," not "experiment-type
detection."  The internal function call is the same
(`_resolve_after_program`) with the same files-win contract, but
the LABEL was different.  H18's grep audit was driven by labels,
not by callsite enumeration.

**H18.2's fix**: one line — pass `original_files=original_files`
to the v117.2 callsite.  No logic change beyond that.

**The K-test that catches this class of bug** is the
production-faithful end-to-end pattern:

- Mock `_call_llm` to return the EXACT shape production emits
  (`{"stop_conditions": {"stop_after_requested": true}}`).
- Run the full `extract_directives(...)` pipeline.
- Assert on the FINAL state, not on any intermediate function.

K21 in `tst_density_modify_experiment_type.py` follows this
pattern.  When the H18.2 fix is reverted, K21 fails
immediately, identifying the offending callsite by the final
`after_program` value.

**Why important for future work**: the same auditing failure
mode could recur for any future ship that adds an optional
parameter to a shared function.  The defensive practice is to
grep for ALL callsites of the modified function — not just the
callsites in code regions that match the parameter's purpose.

```bash
# Defensive grep before declaring a function-signature change "done":
grep -n "function_name(" *.py | grep -v "^.*def function_name"
```

Then audit each line: does the new parameter affect behavior at
this site?  If so, pass it.  If not, document why.  Don't leave
callsites silently using the legacy default unless that's the
intent.

**Generalize**: any function that takes a behavior-controlling
optional parameter has an implicit lockstep contract across all
its callsites.  Adding the parameter to the function definition
is step 1.  Step 2 is auditing every callsite.  Step 3 is a
production-faithful K-test that drives through the entry-point
contract (not the function under test directly) so any missed
callsite surfaces immediately.

### 44. PHIL Multi-Line String Authoring: the Continuation-Space Rule (v120)

**Symptom**: `phenix.python tst_h18_1_phil_roundtrip.py` failed
with `RuntimeError: Syntax error: expected "=", found "is"
(input line 433)` — and after fixing that, `found "set"`.  The
parse of `programs/ai_agent.py`'s `master_params` string was
aborting partway through.  Nothing was wrong with the parameter
being parsed; the error pointed at words (`is`, `set`) in the
middle of a `.help` description.

**Root cause (verified against the real `libtbx.phil` parser, not
inferred)**: a multi-line **unquoted** `.help` value continued
with a backslash, but the backslash was preceded by a non-space
character:

```
.help = Provider for AI analysis. Ollama is cheapest,\     ← BAD: "...cheapest,\"  (comma then backslash)
     OpenAI is most thorough, ...
```

PHIL's lexer requires a **space before the continuation
backslash** in an unquoted value.  Without it, the continuation
is not recognized; the next line is re-lexed as a fresh
statement, PHIL reads its first token (`OpenAI`, `Normally`) as a
parameter name and then expects `=`, finding the next word
instead (`is`, `set`).  The offending word is a red herring — the
real fault is the missing space.

**The rule, established empirically with the parser:**

| Form | Line ends with | Parses? |
|------|----------------|---------|
| Unquoted multi-line | `<non-space>\` | **NO** — always fails, regardless of continuation content (even a single plain word breaks) |
| Unquoted multi-line | `<space>\`     | yes |
| Quoted multi-line   | `"…"\`         | yes — quoted continuations follow different rules and do **not** need the leading space |

The minimal fix is therefore a **single space before the
backslash**, not quoting and not rewording:

```
.help = Provider for AI analysis. Ollama is cheapest, \    ← GOOD: space before backslash
     OpenAI is most thorough, ...
```

**Scope of the v120 fix**: only THREE unquoted helps were broken
— `provider` (added in v120; the regression that surfaced this),
plus the pre-existing `url` and `port` helps (both continued
`...server\` then `Normally set automatically`).  The sibling
`url_type` and `token` helps were already correct because they
happened to end `. \` (space before backslash).  The fix added
one space to each of the three; the other ~30 multi-line helps in
the two `master_params` strings were left untouched.  Both
`programs/ai_agent.py` and `programs/ai_analysis.py` carry dual
schemas (see §42), so the same three were fixed in both.

**Why the quoted scope `.caption` blocks were never affected**:
the `ai_analysis` scope caption is a multi-line **quoted** string
(`"…"\` per line).  Quoted continuations parse fine without the
leading space, which is why they were never flagged.

**Authoring guidance (avoid this class entirely):**
- Prefer a **single-line** `.help`/`.caption` string when it fits.
- For genuinely long text, either (a) use unquoted continuation
  with a **space before every backslash**, or (b) use a quoted
  multi-line string (`"part one "` / `"part two"` per line).  Pick
  one style per block; do not mix.
- Do **not** assume "it compiles as Python, so it parses as PHIL."
  The triple-quoted `master_params` is just a Python string;
  `py_compile` says nothing about PHIL validity.

**Verification recipe** (works in a plain sandbox — no full
PHENIX needed): `pip install cctbx-base --break-system-packages`
gives a real `libtbx.phil`.  Then extract the `master_params`
triple-quoted block with a regex and call
`libtbx.phil.parse(block, process_includes=True)`.  The parser
stops at the FIRST error, so fix-and-re-parse iteratively until it
passes, then scan for any remaining unquoted multi-line help that
ends in `<non-space>\` to catch latent failures the parser has not
yet reached.  This is exactly how the v120 fix was confirmed
complete (zero remaining unquoted `<non-space>\` continuations in
either file).

### 45. Static Never-Inject Params: Proactive `bad_inject_params` (v120)

**Problem.** The `bad_inject_params` mechanism (§ "bad_inject_params flow")
was purely *reactive*: a parameter is added to the per-program blacklist by
`record_bad_inject_param()` only AFTER PHENIX rejects it once.  So the first
command for a program always included a known-bad param, failed, and only then
taught the agent to skip it.  Concretely (log `ai_agent_503`): a stray
`.ncs_spec` from an earlier map-symmetry/segmentation step was auto-mapped to
`ncs_file=` for `phenix.resolve_cryo_em`, which does not accept a top-level
`ncs_file` — the cycle failed before the reactive blacklist kicked in.

**Fix — seed the same mechanism with known-invalid pairs.** Rather than add a
parallel system, `AgentSession` gains a static table:

```python
_STATIC_BAD_INJECT_PARAMS = {
    "phenix.resolve_cryo_em": ["ncs_file"],
}
```

- `get_bad_inject_params(program)` returns `learned | static`, so the param is
  blacklisted on the FIRST command (no failed cycle needed).
- `get_all_bad_inject_params()` (NEW) returns the full merged
  `{program: [keys]}` dict (static ∪ learned, all programs) for transmission.

**The server-path subtlety (why two changes, not one).** The failing run was a
**server-mode** run; the command was built server-side before reaching any
client injector.  The server only sees what the client transmits.  The client
was sending `session.data["bad_inject_params"]` (learned-only), so the static
entry never crossed the wire and the first fix (local `get_bad_inject_params`
alone) would NOT have prevented the failure.  The send-site in `ai_agent.py`
was changed to transmit `session.get_all_bad_inject_params()`.  The end-to-end
chain, verified against real source:

```
client: session_info["bad_inject_params"] = get_all_bad_inject_params()   # merged
  → run_ai_agent.py: create_initial_state(bad_inject_params=...)           # dict in state
  → graph_nodes.py BUILD: all_bad.get(program) → set                       # per-program set
  → command_postprocessor.sanitize_command Rule A: STRIPS blacklisted token
```

`sanitize_command` Rule A *strips* an existing `key=value` whose key (full OR
short) is blacklisted — it does not merely decline to inject — so an
LLM-emitted `ncs_file=` is removed, not just left un-added.

**Defense in depth.** The merged set is honored at four layers: LLM decision
guidance (`ai_agent.py` guidelines), server build (`postprocess_command` via
`session_info`), local build (`postprocess_command` via `get_bad_inject_params`),
and the client required-file injector (`_inject_missing_required_files`, which
now fetches the blacklist and skips matching slots).

**Robustness (external-review hardening).** The BUILD-node extraction is
defensive against a JSON round-trip that yields `None`/non-dict:
`all_bad = state.get("bad_inject_params") or {}`; `isinstance` guard;
`set(all_bad.get(program) or [])`.  And because `sanitize_command` matches both
the full dotted key and the short key, `_STATIC_BAD_INJECT_PARAMS` carries a
WARNING: future entries that could collide with a *valid* nested PHIL parameter
of the same program must blacklist the full dotted path, not the short name.
(`ncs_file`/`resolve_cryo_em` is safe — the program has no valid `ncs_file` in
any scope.)

**Why not scope the LLM guidance to the active program** (a reasonable
suggestion): at guidelines-construction time the agent has not yet chosen a
program — that is the LLM's job for the cycle — so there is no active program to
filter to.  The ban list is necessarily global.  With a one-entry static map
plus the small learned set the prompt cost is negligible; if the static map ever
grows large, cap/relevance-filter the *combined* list rather than scope per
active program.

**Extending.** New known-invalid pairs are a one-line addition to
`_STATIC_BAD_INJECT_PARAMS`.  Tests: `tst_resolve_cryo_em_ncs_inject.py`
(11 tests), including a `json.dumps`/`loads` integration gate that drives the
hardened extraction and the real `sanitize_command`.

**Strip, not remap.** resolve_cryo_em has NO `ncs_file` parameter (verified
against `phenix.resolve_cryo_em --show_defaults` — no `ncs_file` token in any
scope); the correct way to pass symmetry is `input_files.symmetry_file=...ncs_spec`.
One could imagine *remapping* a stray `ncs_file=` to `symmetry_file=`, but that
would be wrong here: in the failing run the `.ncs_spec` was an auto-discovered
stray, not part of the plan (the decision targeted "optimized full map from
half-maps" and never mentioned symmetry), and resolve_cryo_em estimates symmetry
internally.  So the intended behavior is to DROP the stray param, which the
blacklist does.  The blacklist matches `ncs_file` only — the legitimate
`input_files.symmetry_file` is a distinct key (full and short) and is preserved
(verified against the real `sanitize_command`).

NOTE: an earlier CHANGELOG entry (v119.H14.2) asserted "resolve_cryo_em DOES
accept `ncs_file=`".  That is WRONG (per `--show_defaults`) — there is no
`ncs_file` parameter; symmetry is `input_files.symmetry_file`.  The H14.2 entry
has been corrected in CHANGELOG.

### 46. `_safe_float` Consolidation + sanity_checker String-Metric Fix (v120)

**Two related issues.**

*Bug.* When model placement is unnecessary, `phenix.ai_agent` skips model
rebuilding and runs `phenix.model_vs_data`.  On that branch metric values reach
`metrics_history` as strings (parsed from program output / JSON round-trip), and
`sanity_checker._check_metric_anomalies` did raw arithmetic on them —
`change = curr_rfree - prev_rfree` (and `f"{x:.3f}"` formatting) — crashing with
`TypeError: unsupported operand type(s) for -: 'str' and 'str'` (reported at
sanity_checker line ~500).  There are TWO vulnerable blocks, not one: the R-free
spike check and the Map CC drop check; both now coerce.

Fix in `sanity_checker._check_metric_anomalies`:
- `prev_rfree/curr_rfree/prev_cc/curr_cc = _safe_float(prev.get(...))` before any
  arithmetic or formatting.
- Guards changed from `if prev and curr:` to `if prev is not None and curr is
  not None:`.  The old truthiness test had a latent bug: a legitimate metric
  value of `0.0` is falsy and would wrongly skip the check.

*Consolidation.* `_safe_float()` had been **duplicated** in five agent modules
(`validation_history`, `metric_evaluator`, `metrics_analyzer`, `structure_model`,
`display_data_model`) — all logically identical (`def _safe_float(val)`: return
`float(val)`, or `None` on `None`/`ValueError`/`TypeError`).  Five copies are a
drift hazard.  v120 collapses them to a single definition in
`utils/run_utils.py` (the shared-helpers home alongside `setup_llms`,
`validate_api_keys`, `normalize_none_string`).  All five former definers plus
`sanity_checker` now `from libtbx.langchain.utils.run_utils import _safe_float`.

**Why `run_utils.py` and not a new agent module.** `run_utils` is the
established home for cross-cutting helpers and — importantly — imports nothing
from the `agent` package, so `agent.* → utils.run_utils` introduces no import
cycle.  (An interim `agent/metric_utils.py` was prototyped and discarded in
favor of the existing utils module.)

**Call-site safety.** Removing a module-level `def` and replacing it with an
`import` keeps every bare `_safe_float(...)` call working (the import binds the
same name in the module namespace).  Verified: per-file call-site counts are
unchanged from the originals (3/7/12/45/16 across the five modules), and the
imported object is identity-equal to `run_utils._safe_float`.

Tests: `tst_safe_float_consolidation.py` (8 tests) — exactly one definition;
no agent module re-defines it; each consumer imports the canonical one; behavior
(string `"0.385"`, `None`, `0.0`, non-numeric → `None`, no raise); and
sanity_checker coerces both metric blocks with the `is not None` guard.  A
negative control confirmed the "no duplicate definitions" test fails if a copy
is re-introduced.

### 47. Resolution-Dependent Stop Targets: Don't Drop session_resolution (v120)

**Symptom.** At sub-2A resolution a run could stop after a SINGLE refinement
with "R-free TARGET REACHED" (observed: 1.57A, R-free 0.252, stopped).

**Mechanism.** The stop evaluator's R-free target is resolution-dependent
(`metrics.yaml` `r_free.by_resolution`): the `[1.5, 2.0]` band gives
`acceptable = 0.25`, but the resolution-INDEPENDENT default is `0.28`.  With
`success_threshold = target - 0.02`, the banded path needs R-free < 0.23 (0.252
does NOT qualify) while the default path needs < 0.26 (0.252 DOES) — so dropping
the resolution flips a "keep going" into a false "target reached".

**Root cause.** `perceive()` (`agent/graph_nodes.py`) chose the stop-evaluation
resolution from `get_latest_resolution(metrics_history)`.  That returns None when
a cycle's structured `analysis` dict lacks a resolution — which happens even
though the data resolution is known and carried in `state["session_resolution"]`
(the per-cycle graph alias of `session.data["resolution"]`, which xtriage/mtriage
populate via `set_resolution()` — see §48 for the full data flow).  `derive_metrics_from_history()` only read resolution
from `analysis[...]` and had no result-text fallback (unlike r_free/r_work/etc.),
so history silently carried `resolution=None`.  `get_target("r_free", None)` then
returned 0.28.

**Fix (defense-in-depth, two layers).**
1. `perceive()` prefers the authoritative session value:
   `resolution = _safe_float(state.get("session_resolution")) or
   get_latest_resolution(metrics_history)`.  This matches the
   session_resolution-first pattern already used ~8x elsewhere in graph_nodes.py
   (the perceive() stop-evaluation read was the lone exception).  `_safe_float` coercion matters: the band
   lookup does `lo <= resolution < hi`, which raises `TypeError` if
   `session_resolution` is ever stored as a string; coercion both prevents that
   crash and lets a non-numeric value fall through to the history value.
2. `derive_metrics_from_history()` gains a resolution text-fallback via a
   dedicated `_extract_resolution()` helper, so history retains resolution even
   when `analysis` omits it.  The helper is deliberately NOT a naive
   `_extract_float(..., "resolution")`: real PHENIX output prints
   `Anomalous Resolution: 2.10` *before* the data `Resolution: 1.74`, so a
   first-match regex grabs the anomalous value.  `_extract_resolution` strips the
   anomalous span, takes the high-resolution limit of a `range: lo - hi` form,
   lets the last match win, requires a decimal (so `Completeness in resolution
   range: 1` is not read as 1.0 A), and bounds `0.5 < res < 20`.  It matches the
   structured extractor on real p9 xtriage/autosol/autobuild text.

**Verification.** Driven against the REAL stack (real `yaml_loader` + real
`metrics.yaml` + real `metric_evaluator`): `get_target("r_free", None) = 0.28`,
`get_target("r_free", 1.57) = 0.25`; the None path's trend reads "TARGET
REACHED", the 1.57 path reads "first refinement"; the text-fallback recovers
1.57.  Tests: `tst_rfree_resolution_stop.py` (10) + committed fixture
`tests/data/session_rfree_resolution_bug.json`; a negative control confirms the
perceive() change is required, and a wiring test (the recovered-from-text case
uses the Anomalous-prefix form) fails if the fallback reverts to a naive pattern.

**Downstream impact of the original bug** (why it mattered beyond one extra
cycle): a skipped second refinement means ordered solvent is never added
(`rules_selector` requires `refine_count >= 2`) and a ligand fitted after the
first refine is merged but never refined against the data, leaving strong
unexplained difference density in the (pre-ligand) output map.

**Side effect to be aware of:** supplying the real resolution also activates the
non-YAML fallback's own resolution-dependent target (`_analyze_xray_trend`:
`dynamic_target = clamp(resolution/10, 0.20, 0.30)`), which previously defaulted
to 0.25 when resolution was None — so that path now stops slightly *less* readily
at high resolution and slightly *more* readily at low resolution (both correct);
the default YAML path is unchanged in mechanism.

### 48. Resolution Data Flow: One Value, Three Names, Two Renames

The data resolution is referred to by **three different key names** as it
travels from the session to the graph, with **two renames** along the way.  This
is the single most confusing part of the resolution handling; the names are NOT
independent values — they are the same number relabeled at layer boundaries.

```
LAYER                         KEY NAME                       set / read by
----------------------------  -----------------------------  -----------------------------
1. Session (truth)            session.data["resolution"]     set_resolution() writes it;
                                                             get_resolution() reads it
        |  ai_agent.py: session_resolution = session.get_resolution()
        v
2. Agent call / wire          session_state["resolution"]    build_session_state() renames
                                                             session_resolution ->
                                                             session_state["resolution"]
                                                             (api_client.py)   [RENAME 1]
        |  run_ai_agent.py: create_initial_state(
        |      session_resolution = session_state.get("resolution"))           [RENAME 2]
        v
3. Graph state                state["session_resolution"]    create_initial_state() writes it
                                                             (graph_state.py); perceive(),
                                                             build()/find_resolution(),
                                                             template_builder, command_builder
                                                             read it
```

**The single source of truth is `session.data["resolution"]`**, populated by
`session.set_resolution(value, source)` (e.g. from xtriage/mtriage).
`state["session_resolution"]` is its **per-cycle graph alias**, rebuilt fresh by
`create_initial_state()` every cycle — it is NOT persisted.  Consequences worth
remembering:

- A persisted `agent_session.json` (a dump of `session.data`) contains
  **`resolution`**, not `session_resolution`.  The absence of a
  `session_resolution` key in saved session files is correct, not a bug.
- The two renames are deliberate: the wire/`session_state` layer uses the plain
  name `resolution` (shared with other request fields), while the graph state
  uses the qualified `session_resolution` to distinguish the
  authoritative-session value from per-cycle metric resolutions in history.
- `build_session_state()` only sets `session_state["resolution"]` when the value
  is not None, and `create_initial_state()` defaults `session_resolution=None`,
  so an unknown resolution stays None all the way through (and the stop evaluator
  falls back to history / resolution-independent targets — see §47).

The names are **load-bearing** across `graph_state.py` (TypedDict), `contract.py`,
`api_client.py`, `template_builder.py`, `command_builder.py`, and several tests,
so they are intentionally NOT unified by renaming; this section documents the
mapping instead.

### 49. Unified Resolution Resolver (v120)

`perceive()` (the stop decision) and `build()` (program-parameter auto-fill) both
need "the data resolution", and they used to compute it with **different chains**:

| node | old chain |
|---|---|
| `perceive()` | `session_resolution` -> `get_latest_resolution(metrics_history)` |
| `build()` (nested `find_resolution()`) | `session_resolution` -> `workflow_state["resolution"]` -> `log_analysis["resolution"]` -> previous-command regex |

Two chains over the same state can disagree — e.g. `perceive()` decides "stop" from
a loose history default while `build()` prepares parameters from a different
resolution. v120 replaces both with one module-level helper in `graph_nodes.py`:

```python
def resolve_session_resolution(state, workflow_state=None, metrics_history=None):
    # priority (first hit wins), returns (resolution, source_label):
    #   1. state["session_resolution"]            (coerced; the section-48 alias)
    #   2. workflow_state["resolution"]
    #   3. state["log_analysis"]["resolution"]
    #   4. get_latest_resolution(metrics_history) (perceive's structured source)
    #   5. resolution=NN from a previous command  (last resort)
```

- Each caller passes the optional sources it has in scope: `perceive()` passes
  `metrics_history`; `build()` passes `workflow_state`.  A missing source is
  skipped (so `build()` keeps exactly its original source set — it does not gain
  `metrics_history` — and `perceive()` keeps its metrics source while *also*
  gaining the richer fallbacks).
- Every value is coerced with `_safe_float` and range-checked `0.5 < res < 20`.
  This both prevents a stringized `session_resolution` from crashing the band
  lookup (`lo <= resolution < hi`) and stops the old `find_resolution()` failure
  mode where a raw string like `"2.1"` could be returned into
  `invariant_context["resolution"]` and thence a command parameter.

**Behavior preservation.** `perceive()` is unchanged vs its old
`session_resolution or get_latest_resolution(...)` logic (verified by a 200-case
scan).  `build()` is unchanged vs the old `find_resolution()` except for the
intended string->float coercion and range-guard (a 2000-case scan shows zero true
value divergences; all differences are same-value type coercions or correctly
dropped out-of-range values).  `perceive()` and `build()` now agree whenever the
resolution comes from a source they share (session / workflow_state / log_analysis
/ history).  Tests in `tst_rfree_resolution_stop.py` cover the priority tiers, the
coercion/range guard, and an explicit perceive-vs-build agreement check.

**A third resolver on the live build path.** There is a *third* place resolution
is resolved: `CommandContext.from_state()` in `command_builder.py`.  When
`USE_NEW_COMMAND_BUILDER=True` (the default), `build()` delegates to
`_build_with_new_builder`, which constructs a `CommandContext` via
`from_state()`; that object's `context.resolution` is then consumed by
`round(context.resolution, 1)` and `"%.1f" % context.resolution` in
`_apply_invariants` (the resolution auto-fill).  `from_state()` uses its own
3-source chain (`session_resolution` -> `state["resolution"]` ->
`workflow_state["resolution"]`) -- different again from perceive()/build(), and
it ran UNCOERCED, so a stringized source could reach `round()` and raise
`TypeError: type str doesn't define __round__`.  (The graph_nodes `build()`
unification above only covers the *legacy* `USE_NEW_COMMAND_BUILDER=False` path;
`from_state()` is the live consumer.)

**Consolidated coercion primitive.** Rather than unify the three *source chains*
(they legitimately see different state: perceive() has metrics_history, build()'s
context has none, etc.) or create a `command_builder -> graph_nodes` import
cycle, the shared piece -- the coercion+range-guard *rule* -- is factored into
one primitive, `_coerce_resolution(val)` in `utils/run_utils.py` (no agent
dependencies, so no cycle):

```python
def _coerce_resolution(val):
    v = _safe_float(val)
    if v is not None and 0.5 < v < 20.0:
        return v
    return None
```

Both `resolve_session_resolution` (every tier) and `from_state()` (every source)
now run candidates through it.  `from_state()`'s coercion is behavior-preserving
except for the two intended improvements: a stringized value is coerced to float
(removing the `round()` crash) and a 0.0/negative/out-of-range value is rejected
(the `or` chain skips to the next source, exactly as a falsy value did before).

## Workflow States


### X-ray Crystallography Workflow

```
xray_initial → xtriage → xray_analyzed
                              │
              ┌───────────────┼───────────────┐
              ↓               ↓               ↓
       predict_and_build    phaser         autosol†
              ↓               │               ↓
    xray_has_prediction       │        xray_has_phases
              ↓               │               ↓
  process_predicted_model     │          autobuild
              ↓               │               │
    xray_model_processed ─────┘               │
              ↓                               │
           phaser ────────────────────────────┘
              │
              ├── [anomalous data] → xray_mr_sad
              │                         ↓
              │                 autosol (partpdb_file=PHASER.pdb)
              │                         ↓
              │                  xray_has_phases
              │                         ↓
              │                    autobuild ──┐
              ↓                               │
       xray_has_model                         │
              ↓                               │
         refine (loop) ←──────────────────────┘
              ↓
       xray_refined
         ↓    ↓    ↓
    molprobity refine STOP
              │
              ↓
       [if ligand] ligandfit → pdbtools → refine → polder → validate
```

**Ligand-fitting sub-workflow (v115.09b):** When a ligand file is
present and the user requests fitting, the plan template
`refine_placed_ligand` drives a 6-stage pipeline: data_assessment →
refinement → ligand_fitting → final_refinement → ligand_validation →
validation. Three mechanisms ensure the pipeline completes:

1. **General `after_program` resolver** (`directive_extractor.py`,
   v115.10): Detects "fit ligand" as a single action. If the user
   also mentions "refine" or other actions, after_program is cleared
   (multi-step → plan drives). If the user says "fit ligand and
   stop", after_program is set to ligandfit.  Replaces the per-
   workflow ligand after_program clearing from v115.09b.
2. **`combine_ligand` guard** (`workflow_engine.py`): Forces
   `valid_programs = ["phenix.pdbtools"]` for the combine step,
   regardless of what `_apply_directives` returns.
3. **Post-ligandfit exemption** (`workflow_engine.py`): Defers
   `after_program_done` during combine and refine steps so the
   workflow can complete even if the LLM's `after_program` wasn't
   cleared (defense-in-depth).  After the v116.x stop-after
   routing (see Layer 4 historical note) this is less critical
   than it once was — without `stop_after_requested=True`, the
   workflow_engine does not act on `after_program_done` at all, so
   the exemption's original concern is moot — but the exemption
   is retained: it still does useful work prioritizing refine
   programs at the right cycle, and the combine_ligand guard at
   point (2) still needs to override `_apply_directives` because
   pdbtools must be the sole choice at that step.

**Pose file exclusion (v115.09b):** LigandFit produces individual
pose files (`ligand_fit_1_pose_5.pdb`) alongside the final combined
model (`ligand_fit_1.pdb`). Pose files are filtered at three levels:
(1) `session.get_available_files()` removes `_pose_`/`_pose.` from
the returned list — the primary fix; (2) `workflow_state.py`
categorization excludes pose files from `ligand_fit_output`; (3)
`programs.yaml` `exclude_patterns: [pose]` on pdbtools ligand input.

† Anomalous signal classification (from `_analyze_history` in `workflow_state.py`):
  - `strong_anomalous=True`: measurability > 0.10 or `anomalous_resolution` < 6.0 Å → autosol prioritized
  - `has_anomalous=True` (weak): measurability ≥ 0.06, or xtriage `has_anomalous=True` explicitly → autosol available
  - negligible: measurability < 0.06 and xtriage did not assert `has_anomalous` → `has_anomalous=False`, autosol unavailable
  - v115.05 additional guard: when measurability < 0.05 and `has_anomalous` is not True (i.e., False or absent), autosol is removed from `valid_programs` entirely (prevents the LLM from choosing it even when anomalous data columns exist in the MTZ). The "is not True" check (rather than "== False") also handles the case where `has_anomalous` is absent from the context entirely, as occurs with AlphaFold prediction workflows.
‡ MR-SAD: phaser places model first, then autosol uses it as partpdb_file

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `xray_initial` | Starting point | xtriage |
| `xray_analyzed` | After data analysis — includes `probe_placement` phase (maps to same external state) | predict_and_build, phaser, autosol, model_vs_data (probe only) |
| `xray_has_prediction` | Have AlphaFold model | process_predicted_model |
| `xray_mr_sad` | After phaser + anomalous data (MR-SAD) | autosol (with partpdb_file) |
| `xray_has_phases` | After experimental phasing | autobuild |
| `xray_has_model` | Have placed model | refine |
| `xray_refined` | After refinement **or** validation (`refine` and `validate` phases share this external state; internal `phase_info["phase"]` distinguishes them). If R-free ≥ 0.45 after MR + ≥1 refinement cycle + autobuild completed, routes back to `obtain_model` (retry with different search model). The `autobuild_done` guard (v115.05) ensures the agent tries rebuilding before concluding the MR is wrong. | refine, molprobity, autobuild, STOP |

### Cryo-EM Workflow

```
cryoem_initial → mtriage → cryoem_analyzed
                                │
                ┌───────────────┼───────────────┐
                ↓               ↓               ↓
      predict_and_build    dock_in_map    [half-maps only]
                │               │               ↓
                └───────┬───────┘         resolve_cryo_em
                        ↓                       ↓
                cryoem_has_model          map_sharpening
                        ↓                       │
              real_space_refine (loop) ←────────┘
                        ↓
                 cryoem_refined
                   ↓    ↓    ↓
              molprobity RSR  STOP
```

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `cryoem_initial` | Starting point | mtriage |
| `cryoem_analyzed` | After map analysis — includes `probe_placement` phase (maps to same external state) | predict_and_build, dock_in_map, map_correlations (probe only) |
| `cryoem_has_model` | Half-map optimisation / model check (legacy name — no model yet; `check_map` and `optimize_map` phases) | resolve_cryo_em, map_sharpening |
| `cryoem_docked` | Model docked, ready for first real-space refinement (`ready_to_refine` phase) | real_space_refine |
| `cryoem_refined` | After refinement **or** validation (`refine` and `validate` phases share this external state; internal `phase_info["phase"]` distinguishes them) | real_space_refine, molprobity, STOP |

## Model Placement Detection

When a user supplies an atomic model together with reflection data or a cryo-EM
map but no session history, the agent must decide whether the model is already
positioned in the unit cell / map before choosing the next program.  A three-tier
framework resolves this automatically, with each tier more expensive but more
definitive than the last.

### Tier 1 — Unit cell comparison (free, instant)

`agent/placement_checker.py` reads the unit cell parameters from each file and
compares them with a 5% fractional tolerance.

| Source | Reader |
|---|---|
| PDB CRYST1 record | `read_pdb_unit_cell()` |
| MTZ file | `read_mtz_unit_cell()` (iotbx.mtz, falls back to mtzdump) |
| CCP4/MRC map | `read_map_unit_cells()` — returns **two** cells: full-map and present-portion |

- If any read fails → **fail-safe: no mismatch declared** (workflow falls through to Tier 2)
- Definitive mismatch (> 5% on any parameter) → model cannot be placed →
  route immediately to **MR** (X-ray) or **docking** (cryo-EM)
- Cryo-EM: model is compatible if it matches *either* the full-map cell or the
  present-portion cell (partial maps are common)

### Tier 2 — Existing heuristics (`_has_placed_model`)

Checks history flags (`refine_done`, `dock_done`), file subcategory (`positioned`),
and user directives.  When any of these give a clear signal the framework is done.
When the result is still ambiguous, Tier 3 runs.

**MR-keyword guard:** Before inferring placement from constraint keywords (e.g.
"refinement"), `_has_placed_model` first scans constraints for MR intent keywords
(`"molecular replacement"`, `"phaser"`, `" mr "`, `"autobuild"`, `"place the model"`,
`"molecular_replacement"`).  If any are found, the placement-keyword inference is
skipped — those refinement references are future goals, not evidence of current
placement.  This prevents search-model PDB files from being misidentified as placed
models when the README describes an MR workflow.

The `build_context()` method in `WorkflowEngine` computes a `placement_uncertain`
flag that is `True` exactly when all of the following hold:

- `has_model` and (`has_data_mtz` or `has_map`)
- `has_placed_model` is False (Tier 2 found no evidence)
- `cell_mismatch` is False (Tier 1 found no mismatch)
- `placement_probed` is False (probe has not run yet)
- `has_predicted_model` is False (predicted models always need processing/docking)

### Tier 3 — Diagnostic probe (one program cycle)

When `placement_uncertain` is True the workflow engine routes to the
`probe_placement` phase, which runs a single quick diagnostic:

| Experiment | Program | Threshold | Placed | Not placed |
|---|---|---|---|---|
| X-ray | `phenix.model_vs_data` | R-free < 0.50 | → refine | → molecular_replacement |
| Cryo-EM | `phenix.map_correlations` | CC > 0.15 | → refine | → dock_model |

After the probe cycle, `_analyze_history` in `workflow_state.py` detects the
result **positionally**: the first occurrence of `model_vs_data` or
`map_correlations` that appears *before* any refinement or docking cycle is
identified as the probe.  This requires no schema change to history entries.

`build_context()` sets `placement_probed = True` and
`placement_probe_result = "placed" | "needs_mr" | "needs_dock" | None` (None if
the result could not be parsed — fail-safe: falls through to obtain_model).
When the result is `"placed"`, `build_context` overrides `has_placed_model = True`
so normal refinement routing takes over.

The probe never repeats: `placement_probed = True` in `history_info` prevents
`placement_uncertain` from being set on subsequent cycles.

### Routing summary

```
model + data, no history
        │
        ▼ Tier 1
   Cell mismatch? ─── YES ──▶ MR / docking (skip probe)
        │ NO
        ▼ Tier 2
  Has heuristic evidence? ─── YES ──▶ normal routing
        │ NO (placement_uncertain = True)
        ▼ Tier 3
   probe_placement phase
    (model_vs_data or map_correlations)
        │
   R-free / CC result
        ├── placed    ──▶ refine
        ├── not placed ──▶ MR / docking
        └── unparseable ──▶ obtain_model (fail-safe)
```

### Key files

| File | Role |
|---|---|
| `agent/placement_checker.py` | Unit cell readers and comparison |
| `agent/workflow_engine.py` | `build_context()` new keys; `_detect_*_phase()` routing |
| `agent/workflow_state.py` | `_analyze_history()` probe detection |
| `knowledge/workflows.yaml` | `probe_placement` phase in xray and cryoem |
| `knowledge/programs.yaml` | `done_tracking` for `map_correlations` (pre-existing gap fixed) |
| `agent/yaml_tools.py` | `if_placed` / `if_not_placed` added to `valid_transition_fields` |

### Implementation notes

**Import paths**: `placement_checker` is importable via both `libtbx.langchain.agent.placement_checker`
(production) and `agent.placement_checker` (tests / local dev). `_check_cell_mismatch`
tries both paths — libtbx first, bare `agent` path as fallback — matching the pattern
used throughout the codebase.

**Short-circuit**: `build_context()` computes `cell_mismatch` by calling `_check_cell_mismatch()`
on every cycle. For cryo-EM this involves running `phenix.show_map_info` as a subprocess.
After the context dict is built, a post-processing block overrides `cell_mismatch = False`
when either `has_placed_model` or `placement_probed` is `True` — both conditions mean
placement is already resolved and the check cannot change the outcome. The cell check
still runs unconditionally on the first cycle when placement is genuinely unknown.

**YAML validator**: `if_placed` and `if_not_placed` are registered in
`valid_transition_fields` in `yaml_tools.py` so that `_validate_workflows()` does not
emit spurious "unknown field" warnings for `probe_placement` phases.

### v115.09 Routing Additions

**Cryo-EM `past_analysis` gate** (`_detect_cryoem_step`): The gate that
advances cryo-EM routing past the "analyze" step now includes
`map_sharpening_done`, `map_symmetry_done`, and `has_optimized_full_map`
(file-presence fallback). Previously, tutorials that started with
`map_sharpening` instead of `mtriage` got stuck in "analyze" with no
programs available.

**Validation-only shortcut** (`_detect_xray_step`): When
`wants_validation_only=True` (from directives) and both `has_model` and
data are present (`has_data_mtz` OR `has_phased_data_mtz` — completed
structures often have phase columns in their MTZ), routing jumps directly
from xtriage to the `validate` step. The `validate` step runs
`model_vs_data` first (crystal symmetry sanity check) then `molprobity`.
The corresponding `validate_existing` plan template (priority 60) ensures
the planner selects a 2-stage plan (data_assessment → validation).

**`_is_valid_file` PDB scan limit** (`workflow_state.py`): The Layer 3
structural validity check scans PDB files for ATOM/HETATM records. The
scan limit is 2000 lines (increased from 500 after 3dnd.pdb — 546 header
lines — was rejected as invalid). Since the full file content is already
read into memory, the `any()` short-circuits on first match with no
performance cost.

**`force_mr` flag** (`build_context` + `_detect_xray_step`): When
`use_mr_sad=True` from directives but the PDB is categorized as `model`
(not `search_model`), `force_mr=True` is set. This overrides placement
probes and routes directly to `molecular_replacement`, where phaser's
`has_any: [processed_model, model_for_mr]` condition is satisfied via
`has_model_for_mr=True`. The MR-SAD guard in `get_valid_programs` also
checks `force_mr`, blocking autosol until phaser completes.

**Directive-driven intent**: Both `wants_validation_only` and
`use_mr_sad` are extracted by the directive extractor, not by string
matching in the engine. The engine reads these as boolean flags
from `workflow_preferences`.

**Extraction architecture**: The rules-based patterns in
`_apply_workflow_intent_fallback()` run as a **post-LLM overlay** —
after the LLM returns directives, the overlay applies deterministic
pattern matching and sets routing flags that the LLM missed. This
is necessary because the LLM has never set these flags in 240+
tested extractions (across 61 tutorials × 4 modes). The overlay
also runs in `extract_directives_simple()` for the rules-only path.
Rules always run last and always win for routing flags.

**`map_sharpening_done` regex** (`workflow_state.py`): The zombie check
table and done-flag detection use regex `sharpen.*\.(ccp4|mrc)$` (not
`sharpened`) to match actual `phenix.map_sharpening` output filenames
like `auto_sharpen_A.ccp4` and `bgal_auto_sharpen.ccp4`.

**Future refactor** (v115.10): Replace overlay with centralized
`_DIRECTIVE_SCHEMA` and registry-driven `_merge_tiered` merge where
each field has a declared authority level (RULES or LLM). See
`docs/directive_merge_plan.md`.

**Preprocessing stop override (v115.09b, fixed):** The
`_preprocessing_programs` set (`xtriage`, `mtriage`) previously cleared
`after_program` unconditionally, even when the user explicitly said
"run mtriage and stop." Fixed with `_has_explicit_stop` regex check
before clearing — raw user advice is checked for explicit stop signals
before the preprocessing override fires.

**General `after_program` resolver (v115.10):** Replaces per-workflow
overlays (ligand-fit clearing, denmod stop/clear) with a single
general mechanism.  The resolver uses `_ACTION_TABLE` (14 actions
with keywords mapped to xray/cryoem programs) to detect which
workflow actions the user mentioned, then applies three rules:

- **Multiple actions + "stop"** → `after_program` = last mentioned
  action's program.  "run phaser, refine, and stop" → `after_program
  = phenix.refine`.
- **Multiple actions, no "stop"** → clear `after_program`.  Plan
  template drives the full workflow.  "density modify and build" →
  no `after_program`.
- **Single action + "stop"** → set `after_program` if LLM missed it.
  "density modify and stop" → `after_program = resolve_cryo_em`.

Key features: word-boundary matching for single-word keywords (prevents
"solve" matching inside "resolve"); action-specific negation detection
("don't build" removes the build action); predict+build compound rule
(merged into single `predict_build` action); experiment-type inference
from advice text for xray/cryoem program mapping.

Also removed from `extract_directives_simple`: `continuation_indicators`
(9 patterns), `downstream_tasks` (16 patterns), `multi_program_patterns`
(8 patterns) — all replaced by `_detect_actions()`.  `tutorial_patterns`
and `denmod_patterns` retained for initial extraction; the general
resolver overrides when wrong.

**Stop-after directive routing (v116.x, complement to v115.10):** The
general resolver at v115.10 fixed one side of `after_program`: making
sure user advice was correctly mapped to the right `after_program`
value (or no value, for multi-step requests without an explicit stop).
Plan-driven progression (`plan_to_directives`,
`knowledge/plan_schema.py`) also emits `after_program` per stage as a
min-run hint, by design (see [Plan section](#) above).

The result was that `after_program` carried two distinct kinds of
intent which the consumer (`workflow_engine._apply_directives`) could
not distinguish:

  (a) **USER/README explicit stop** ("refine and stop") — the user
      wants the workflow to stop after this program completes.

  (b) **PLAN progression hint** (per-stage `after_program` from
      `plan_to_directives`) — the plan wants to make sure THIS
      stage's program runs.  No user intent to stop is implied.

For a period between v112.78 and v116.x, the consumer wiped
`valid_programs` to `[STOP]` whenever the named program had run,
treating BOTH cases as a hard stop.  Symptoms:

- Multi-goal user requests (e.g. "refine the model and fit ATP using
  LigandFit. Stop Condition: None") were truncated to just the first
  goal.  The directive extractor or plan-to-directives emits
  `after_program=phenix.refine` for the refinement stage; once refine
  succeeded, the wipe killed ligandfit from `valid_programs` even
  though it was the user's stated second goal.

- Plan-driven progression through stages was disrupted: each stage's
  `after_program` hint caused the workflow to "complete" at every
  stage boundary rather than advancing.

A v112.78-restoration patch ("append STOP, don't wipe") was tried
first — but it left two problems standing.  First, the LLM still got
a misleading "Stop target: X" prompt directive for plan-injected
`after_program`, biasing it toward stopping mid-plan.  Second, for
user-explicit stops the LLM was trusted to pick STOP based on the
prompt; that worked most of the time but offered no structural
guarantee.

**The v116.x architecture distinguishes the two cases at the source
and uses different machinery for each.**  A new directive field,
`stop_conditions.stop_after_requested` (bool), is set by the directive
extractor when — and only when — the user (or README) explicitly
requested a stop-after condition.  Detection lives in
`_is_stop_after_requested()` (`agent/directive_extractor.py`); it
returns True for phrases like `"stop after X"`, `"X and stop"`,
`"only run X"`, `"just (do|run) X"`, and `"Stop Condition: <real
value>"`.  It returns False for:

- The absence of any stop signal
- `"Stop Condition: None"` / `"Stop Condition: not specified"` /
  heading-only `"Stop Conditions:"` (an explicit no-stop signal that
  the previous bare-`\bstop\b` check matched as a false positive —
  the source of Tom's nsf-d2-ligand `Stop Condition: None` bug)
- Explicit negations like `"don't stop"` / `"do not stop"` / `"never
  stop"`

Plan-injected `after_program` from `plan_to_directives` does NOT set
`stop_after_requested` — `plan_to_directives` only sets
`after_program`, leaving the new flag as its default (absent / False).
This is the structural marker that distinguishes plan progression
from user stop.

Two consumers gate on the flag:

1. **`workflow_engine._apply_directives`** (the `elif
   after_program_done:` branch):

   - `stop_after_requested=True` AND target program has run →
     **wipe `valid_programs` to `[STOP]`**.  The user said stop;
     enforce it structurally.  No ambiguity, no risk of the LLM
     running molprobity instead of stopping.
   - `stop_after_requested=False` AND target program has run →
     **do nothing**.  Let the workflow advance to the next plan
     stage / step naturally.  No spurious STOP prompts.

   The min-run guarantee branch (`after_program` not yet done) is
   unchanged: the target program is added to the front of
   `valid_programs` regardless of `stop_after_requested`, because
   we always want the named program to run before we make any
   stop decisions.

2. **`prompts_hybrid.py::_format_directives_for_prompt`** (the
   `**Stop Conditions:**` section):

   - `stop_after_requested=True` → emit the "Stop after X completes" +
     "Stop target: X" v116.10 Phase 6a wording.
   - `stop_after_requested=False` → emit NO `after_program`-related
     prompt lines.  The LLM is never told "stop after X" when the
     workflow is meant to continue past X.  The workflow_engine's
     min-run prioritization of X in `valid_programs` is sufficient
     guidance.

A side effect: the `**Stop Conditions:**` header is now emitted only
when there is non-empty body content under it.  The previous
unconditional emission would have produced an empty header in the
plan-progression case.

**`skip_validation` is no longer special-cased here.** The directive
extractor still sets `skip_validation=true` automatically on any
user-explicit stop (per its "ALWAYS set" rule), but the consumer
does not gate on it.  The wipe enforces the user's explicit stop
intent whether or not validation programs are in the YAML step.  A
`skip_validation`-based carve-out preserving the wipe just for the
validate step was considered and rejected during the v112.78
restoration: it would have re-introduced the premature-stop bug for
any "X and stop" workflow that reaches validate, because
`skip_validation` is set broadly.

**Why "just trust the LLM" wasn't enough on its own:** Under the
v112.78-restoration patch, the prompt's strong directive (`Stop
target: X`) was the mechanism that biased the LLM toward STOP when
the user requested it.  Empirically the LLM respected this most of
the time, but there were two failure modes:

1. The prompt was emitted ALSO for plan-injected `after_program`,
   confusing the LLM into stopping mid-plan.
2. There was no structural backstop if the LLM mis-read the
   directive.

The v116.x design uses the directive flag for both: prompt suppression
prevents the LLM from getting the wrong message, and structural wipe
guarantees the right outcome when the user actually said stop.

**Tracking program success vs attempt (Bug 3 fix):** Distinct from
the stop-routing change but shipped alongside it.  Previously,
`after_program_done` flagged as true whenever the named program
appeared in history, regardless of whether it succeeded.  A failed
`phenix.refine` would mark `after_program_done=True`, then the
(then-wiping, now-wiping-only-on-user-stop) code would prevent any
retry.

`workflow_state.py::_analyze_history` now populates a
`successful_programs` set alongside the existing `programs_run` set.
A program is added to `successful_programs` only if its history
entry doesn't match `_is_failed_result()`.  `workflow_engine`
propagates this through `build_context` and uses it in
`_apply_directives` to override `after_program_done = False` when
the target program is in `programs_run` but not in
`successful_programs`.  The `CONSEC-FAIL` guard at
`graph_nodes.py` line ~1844 still caps retries at 2 consecutive
failures.

**Parameter handling architecture (v115.10 design note):** User parameter
requests ("set resolution to 3", "run 4 cycles of refinement", "box the
map") use a different architecture from the stop/continue resolver above.
Stop/continue is a *routing* decision with a small discrete solution space
(14 actions) where the LLM was systematically wrong — a general resolver
was needed.  Parameters are *value* decisions where the LLM is generally
reliable at extraction, and the solution space is large (hundreds of PHIL
parameters per program).  The parameter system uses three layers:

1. **LLM extraction** → extracts user intent as concept + value
   (e.g., ``resolution: 3``, ``macro_cycles: 4``).  Rules-based
   patterns (``cycle_patterns``) provide fallback for common cases.
2. **``strategy_flags`` in programs.yaml** → maps concept names to
   program-specific PHIL paths (e.g., ``resolution`` →
   ``xray_data.high_resolution={value}`` for refinement vs
   ``crystal_info.resolution={value}`` for predict_and_build).
3. **``command_builder``** → applies mapped parameters to the command.

When a user parameter request fails, the fix is typically to add a
``strategy_flags`` entry in programs.yaml for the missing concept,
not to build new Python extraction code.  The ``auto_fill_resolution``
invariant is an example: it automatically fills the resolution
parameter for programs that need it, using the value from previous
program outputs.

**.sca-only data detection** (`perceive`): When all data files are
`.sca/.hkl` with no `.mtz`, no model, and no sequence, and no
`unit_cell` is provided in directives, `perceive()` emits a helpful
abort message explaining that unit cell parameters must be provided
via the GUI.

---

## Workflow History and Done Flags

`workflow_state.py` maintains two complementary systems for tracking
completed programs:

### _analyze_history: done flag extraction

`_analyze_history(history)` reads the history list and sets done flags
(`refine_done`, `resolve_cryo_em_done`, …) based on marker strings in
each record's `result` field. Key invariants:

- **Failure blocks flags**: `_is_failed_result(result)` is checked
  first. If the result is a failure, all done-flag extraction for that
  record is skipped. Failure signals (in priority order):
  1. Shell exit code (caught at the subprocess layer before this function)
  2. Specific Phenix terminal phrases: `FAILED`, `SORRY:`, `SORRY `,
     `*** ERROR`, `FATAL:`, `TRACEBACK`, `EXCEPTION`
  - Generic `ERROR` words without those prefixes are **not** failure
    signals — Phenix logs contain "Error model parameter", "Expected
    errors: 0", etc. and these must not suppress done flags.

- **predict_and_build cascade**: `predict_full_done` additionally sets
  `refine_done=True` and increments `refine_count` (Python-only; no
  YAML `history_detection` entry) because a full `predict_and_build`
  run includes refinement internally.

### _clear_zombie_done_flags: crash recovery

After `_analyze_history`, `detect_workflow_state` calls
`_clear_zombie_done_flags(history_info, available_files)`.

A **zombie state** occurs when the agent crashed mid-cycle or the user
deleted output files: history records `done_flag=True`, but the output
file is missing from disk. The phase detector sees `done=True` and
skips the program; the workflow becomes stuck.

The function checks each crashable done flag against its expected
output filename pattern:

| Done flag | Expected output | Also clears | Any-PDB fallback |
|-----------|----------------|-------------|-----------------|
| `resolve_cryo_em_done` | `denmod*.ccp4/mrc` | `has_full_map` | No |
| `predict_full_done` | `*_overall_best.pdb` | `has_placed_model` | Yes |
| `dock_done` | `*docked*.pdb`, `*dock*map*.pdb`, `placed_model*.pdb`, `*_placed*.pdb` | `has_placed_model` | Yes |
| `refine_done` | `*_refine_NNN.pdb` | decrements `refine_count` | Yes |
| `rsr_done` | `*_real_space_refined*.pdb` | decrements `rsr_count` | Yes |

**Any-PDB fallback**: when `True`, the presence of any `.pdb` file on disk suppresses
zombie detection even if the pattern doesn't match. Used for model-producing programs
where output may be renamed. `resolve_cryo_em_done` uses `False` because `denmod_map.ccp4`
is referenced by name downstream.

Flags are cleared **in-memory only** — history is never rewritten.
The function returns a `zombie_diagnostics` list which `detect_workflow_state`
attaches to the state dict under the key `"zombie_diagnostics"`. The PERCEIVE
node logs each entry prefixed `"PERCEIVE: ZOMBIE STATE — "`. This makes
crash/restart re-runs self-explaining: the user can see exactly which done flag
was cleared and why the program appears again.

## Best Files Tracking

The agent maintains a **Best Files Tracker** that identifies and tracks the highest-quality file of each type throughout a session, ensuring programs always receive optimal inputs.

### Categories Tracked

| Category | Description | Example |
|----------|-------------|---------|
| `model` | Best atomic model (PDB/mmCIF) | `refine_002_001.pdb` |
| `map` | Best full cryo-EM map | `denmod_map.ccp4` |
| `mtz` | Best reflection data | `data_with_rfree.mtz` |
| `map_coefficients` | Best map coefficients | `refine_001_001.mtz` |
| `sequence` | Sequence file | `sequence.fa` |
| `ligand_cif` | Ligand restraints | `LIG.cif` |

### Scoring System

Each file is scored based on **processing stage** (0-100 points) and **quality metrics** (0-100 points). All scoring parameters are configurable in `knowledge/metrics.yaml` under `best_files_scoring`.

**Model stage scores:** refined=100, rsr_output=100, autobuild_output=80, docked=60, processed_predicted=50, predicted=40, phaser_output=30, pdb=10

**Model metric scores:** R-free (40 pts, linear_inverse), Map CC (30 pts, linear), Clashscore (30 pts, linear_inverse)

**Map stage scores:** optimized_full_map=100, sharpened=90, density_modified=80, full_map=50, half_map=10. Plus resolution bonus (0-30 pts).

### R-free MTZ Locking (X-ray Only)

Once an MTZ with R-free flags is identified, it is **locked for the entire session**. This ensures consistent cross-validation statistics throughout refinement. The locked MTZ is Priority 0 inside `_find_file_for_slot` — but a `data_mtz` preference or LLM hint can pre-fill the slot *before* that point (see the caveat below), so the lock's priority is enforced for the refine data slot by **R-free lock reconciliation** (F1/F2), not by Priority 0 alone.

### File Selection Priority

When building commands, files are selected in this order:

0. **Locked R-free MTZ** — For MTZ inputs in X-ray refinement (highest priority)
1. **Best Files** — From the BestFilesTracker scoring system
2. **Categorized Files** — From workflow_state file categorization
3. **Extension Matching** — Search available_files by extension

> **Caveat (v120.3): preferences fill the slot before Priority 0.**  A
> `directives.file_preferences.data_mtz` or an explicit LLM hint is applied in
> `_select_files` *before* `_find_file_for_slot` injects the Priority-0 lock, so a
> preference pointing at the original (often flagless) input can win over the
> lock.  R-free **lock reconciliation** (below) is the corrective.

### R-free Lock Reconciliation — F1 / F2 (v120.3, widened v120.4)

Two coordinated steps in `command_builder.py` keep the refine data file and the
generate decision consistent with the locked R-free MTZ:

- **F2 — reconcile the data slot to the lock.**  In `build()`, after
  `_select_files` and before `_apply_invariants`, if `rfree_mtz` is locked and the
  data slot holds a *different* file that is **not confirmed to carry flags**, the
  slot is reset to the locked MTZ (logged `BUILD: reconciled … locked R-free MTZ`).
  "Not confirmed flagged" = `_input_mtz_rfree_state` returns `False` (flagless)
  **or** `None` (unverifiable — e.g. a raw scalepack input, which `inspect_mtz`
  cannot read).  A flagged file (`True`) is never overridden, and an `abspath`
  guard (checked first) skips the slot when it already *is* the lock.  v120.3
  reconciled only `False`; **v120.4** widened the trigger to `is not True` so an
  undetermined input (scalepack `p9.sca`, p9-sad AIAgent_95) reconciles too.

- **F1 — generate decision keys off the selected file.**  `_apply_invariants`
  decides `xray_data.r_free_flags.generate` from `_input_mtz_rfree_state(files,
  context)` on the *actually-selected* file (post-F2), not from the mere existence
  of a lock: strip on `True`, add on `False`/undetermined-when-unlocked, with the
  lock as the undetermined fallback.

Net effect: a flagless or unverifiable file selected over the lock is pulled back
to the flagged lock (F2), and the generate flag follows the file that will actually
be refined (F1).  Fail-safe — you never lose R-free flags; worst case you refine
against the locked dataset instead of an unverifiable one.  Pinned by
`tst_rfree_lock_reconciliation.py`.

### Key Files

| File | Role |
|------|------|
| `agent/best_files_tracker.py` | Core tracker class with scoring and STAGE_TO_PARENT mapping |
| `agent/session.py` | Integration with session persistence; supplemental file discovery |
| `agent/template_builder.py` | Uses best_files for command building |
| `agent/file_utils.py` | `classify_mtz_type()` for MTZ classification; `matches_exclude_pattern()` for word-boundary pattern matching |
| `agent/workflow_state.py` | `_pdb_is_small_molecule()` and `_pdb_is_protein_model()` for content-based PDB analysis; MTZ categorization safety net in `_categorize_files()` |
| `knowledge/metrics.yaml` | Scoring configuration (best_files_scoring section) |

### Companion File Discovery

Some clients only track a subset of program output files. The agent discovers
missing companion files in multiple layers:

**Layer 0: `session._discover_cycle_outputs()` (v112.73)** — The foundation
layer.  All other layers depend on files appearing in `available_files`.  This
method resolves output files for any cycle using three strategies:

| Strategy | When it helps |
|---|---|
| Try stored `output_files` paths as-is | Normal operation |
| Resolve relative paths against `_get_session_dir()` | cwd changed between runs |
| Scan `sub_{NN}_{program}/` by cycle number + program name | output_files completely empty |

`_get_session_dir()` returns `os.path.dirname(session_file)` — always known,
never depends on stored file paths.  The `get_available_files()` Step 3
directory scan is also seeded from the session directory, not just from
already-tracked files.  This means even if Steps 1-2 find nothing, the scan
still runs.

Without this layer, all downstream recovery (categorization safety nets,
best_files tracking, etc.) is irrelevant because the files were never in the
working set.

**Intermediate file filtering (v115.09b):** `get_available_files()` applies a
final filter (Step 5) that removes files with `_pose_` or `_pose.` in the
basename. LigandFit produces individual pose files (`ligand_fit_1_pose_5.pdb`)
alongside the final combined model (`ligand_fit_1.pdb`). Without this filter,
`_discover_cycle_outputs` globs `LigandFit_run_1_/*.pdb` and adds all pose
files to `available_files`. The LLM then picks the pose with the highest CC
instead of the final model. Matching is basename-only to avoid directory-name
false positives. Words containing "pose" (decompose, compose, expose) are not
affected because the patterns require an underscore before "pose".

**Layer 1: `session._find_missing_outputs()`** — Runs in
`get_available_files()` after Layer 0.  Derives companion files from
known output file names (e.g., if `refine_001_data.mtz` is found,
looks for `refine_001.mtz`).  Supplements Layer 0's directory scan
with pattern-based inference.

**Layer 2: Best files evaluation (v112.70)** — Both
`_rebuild_best_files_from_cycles` (session load) and `record_result`
(live cycle completion) call `_find_missing_outputs` and evaluate
supplemental files through the best_files tracker. This ensures
`best_files["map_coeffs_mtz"]` is populated even when the client
only tracked `refine_001_data.mtz` in `output_files`. Without this
layer, programs with `require_best_files_only: true` (like
ligandfit's map_coeffs_mtz slot) would fail to build because the
map coefficients MTZ was never evaluated.

**Layer 3: MTZ categorization safety net (v112.71)** — Runs at the
end of `_categorize_files()` after both YAML and hardcoded
categorization paths. Cross-checks every MTZ file against the
authoritative `classify_mtz_type()` regex and corrects three types
of misclassification:

| Failure Mode | Detection | Correction |
|---|---|---|
| File in `data_mtz`, should be `map_coeffs_mtz` | `classify_mtz_type()` returns `map_coeffs_mtz` but file not in that category | Move to `map_coeffs_mtz` + subcategory, remove from `data_mtz` |
| File in BOTH `data_mtz` and `map_coeffs_mtz` | YAML Step 1 extension match + Step 2 pattern match create dual membership | Remove from `data_mtz` (prevents `exclude_categories` rejection) |
| File in `map_coeffs_mtz`, should be `data_mtz` | `classify_mtz_type()` returns `data_mtz` but file in `map_coeffs_mtz` | Move to `data_mtz`, remove from `map_coeffs_mtz` and subcategories |

All corrections are logged at `WARNING` level via Python's logging module.

**Layer 4: Self-contained MTZ classification (v112.73)** — The Layer 3 safety
net depends on `classify_mtz_type()` from `file_utils.py`.  Three deployment
failures can disable it: (a) `get_mtz_stage` not deployed (joint import kills
`classify_mtz_type` too), (b) `file_utils.py` entirely missing, (c) YAML
patterns incomplete.  `_import_mtz_utils()` eliminates all three by **always
returning working functions**.  It tries importing from `file_utils.py` first,
then falls back to inline implementations that embed the refine-output regex
directly in `workflow_state.py`.  No external dependency can break it.

**Layer 5: Principled exclusion rule (v112.73)** — Defense-in-depth in the
command builder for dual-categorization (file in both `data_mtz` and
`map_coeffs_mtz`).  `_should_exclude()` implements the rule: **exclude only
if the file is in an excluded category AND NOT in any desired category**.
With Layer 4 working, dual-categorization is cleaned up before it reaches the
command builder.  This layer catches partial failures where the safety net adds
a file to `map_coeffs_mtz` but doesn't remove it from `data_mtz`.

### MTZ Categorization Diagnostics

Two logging points help diagnose map_coeffs_mtz failures:

1. **`perceive()` node**: After file categorization, logs all MTZ category
   contents: `PERCEIVE: MTZ categories: data_mtz=[...]; map_coeffs_mtz=[...]`.
   Warns if refinement is in history but `map_coeffs_mtz` is empty.

2. **`_categorize_files()` safety net**: Logs `WARNING` when it corrects a
   misclassification (e.g., `MTZ safety net: moved refine_001_001.mtz from
   data_mtz to refine_map_coeffs/map_coeffs_mtz`).

### Intermediate File Filtering

`graph_nodes._filter_intermediate_files()` removes
temporary/intermediate files before categorization.  Runs in the
perceive node after history injection.

**Filtered patterns:**
- Files in `/TEMP`, `/temp`, `/TEMP0/`, `/scratch/` directories
- Files with `EDITED_` or `TEMP_` prefixes

These are internal working files from programs like ligandfit that should never
be used as inputs to other programs.

**`predict_build_refine_internal` category (v115.03):** `predict_and_build`
performs its own internal refinement passes and writes intermediate PDB files
named `overall_best_final_refine_NNN.pdb` and `overall_best_refine_NNN.pdb`.
These are multi-model PDB files (one MODEL record per refinement macro-cycle)
that cannot be passed directly to `phenix.refine` or used as the session model.
The `predict_build_refine_internal` category (`parent_category: intermediate`)
captures these files via glob patterns `*overall_best*final_refine*.pdb` and
`*overall_best*refine_[0-9][0-9][0-9]*.pdb`, ensuring they are never selected
as a model input.  The correct post-predict_and_build output
(`PredictAndBuild_0_overall_best.pdb`) is captured by the existing
`predict_and_build_output` category.

### Best Files Exclusion Check

When `best_files["model"]` is applied as a model override for refinement
programs, it is first checked against the program's `exclude_categories`. If
the best model is in an excluded category (e.g., `ligand_fit_1.pdb` in
`ligand_fit_output` → parent `ligand`), it is skipped and category-based
selection runs instead. This prevents ligand fragments from being used as the
protein model for refinement.

Applied in `command_builder.py` at two points: (1) pre-population of the model
slot from best_files, and (2) LLM override where best_files would normally
take precedence over the LLM's model choice.

### Content-Based File Selection Guards (v112.70)

Three layers of defense prevent files from being assigned to wrong slots:

**Layer 1: `exclude_patterns` (YAML-driven, word-boundary matching)**
Slot definitions in `programs.yaml` specify `exclude_patterns` to reject files
by name. Uses `matches_exclude_pattern()` with word-boundary semantics: `ligand`
matches `atp_ligand.pdb` but NOT `nsf-d2_noligand.pdb`. Applied to both
auto-fill and LLM-selected files.

**Layer 2: Content-based PDB analysis**
After pattern matching, PDB files are checked by content:
- **Model slots** (`model`, `protein`, `pdb_file`): `_pdb_is_small_molecule()`
  rejects HETATM-only files (ligands like `atp.pdb`)
- **Ligand slot**: `_pdb_is_protein_model()` rejects files with ATOM records
  (protein models like `refine_001.pdb`)

Uses `_pdb_is_protein_model()` for the ligand slot rather than
`not _pdb_is_small_molecule()` because the latter returns False for unreadable
files, which would incorrectly reject valid candidates.

**Layer 3: LLM selection validation (v112.70)**
LLM file hint assignments are now validated against the slot's `exclude_patterns`
before acceptance. Previously, the LLM could bypass exclusion rules by explicitly
assigning a file. Now the same guards apply to all file sources:

```
File selection pipeline (all three apply at each stage):
  1. LLM hint → validate against exclude_patterns + content guards
  2. Auto-fill → apply exclude_patterns + content guards + prefer_patterns
  3. Safety net → apply exclude_patterns + content guards + prefer_patterns
```

Applied in:
- `CommandBuilder._find_file_for_slot` (server-side auto-fill)
- `CommandBuilder` LLM hint validation loop
- `_inject_missing_required_files._find_candidate_for_slot` (client-side safety net)

### Ligand Category Content Validation (v112.74)

The YAML categorizer (`file_categories.yaml`) may misclassify protein PDB files
as `ligand_pdb` when broad filename patterns match names like `1aba.pdb` or
`3gx5.pdb`.  A protein with a few HETATM ligand/cofactor atoms is still a
macromolecular model, not a ligand coordinate file.

Post-processing guard in `_categorize_files()` validates every `ligand_pdb`
entry using `_pdb_is_protein_model()`.  If a PDB file has >150 coordinate
records and majority ATOM records, it is a false positive and gets moved from
`ligand_pdb`/`ligand` to `unclassified_pdb`/`pdb`/`model`.  Same defense
pattern as the half-map validation guard.

Only runs when `files_local=True` (file content readable on disk).

| File | ATOM | HETATM | Total | is_protein_model | Action |
|------|------|--------|-------|------------------|--------|
| 1aba.pdb (protein+ligand) | 729 | 20 | 749 | True | Rescue → model |
| atp.pdb (pure ligand) | 0 | 31 | 31 | False | Keep in ligand |
| hem.pdb (cofactor) | 0 | 43 | 43 | False | Keep in ligand |

**Defense layers for PDB ligand classification:**

| Layer | What it catches | Where |
|-------|----------------|-------|
| YAML patterns | Name-based (e.g. `*lig*`) | `_categorize_files_yaml` |
| Hardcoded name check | Word-boundary `lig`/`ligand` | `_categorize_files_hardcoded` |
| Content: unclassified→ligand | HETATM-only in unclassified_pdb | Post-processing (existing) |
| Content: ligand→model | Protein in ligand_pdb | Post-processing (v112.74) |
| best_files_tracker | Content + name check | `_is_ligand_file` |

### Duplicate Command Detection

`session.is_duplicate_command()` prevents the agent from repeating commands.
Two-tier detection:

**Tier 1 — Exact match:** Compares normalized commands against all prior
commands (both successful and failed). Catches verbatim retries.

**Tier 2 — Overlap heuristic:** For same-program commands, computes
basename-level token overlap. >80% overlap flags as duplicate. However, file
tokens (basenames with crystallographic extensions) are compared separately:
if the input files differ, the commands are NOT duplicates regardless of overall
token overlap. This prevents iterative refinement (same program, new model)
from being blocked.

When a duplicate is detected, the agent retries through the normal graph path
with `duplicate_feedback` appended to guidelines (up to 3 attempts).

| File | Role |
|------|------|
| `agent/session.py` | `is_duplicate_command()`, `get_all_commands()`, `get_all_failed_commands()` |
| `programs/ai_agent.py` | `_handle_duplicate_check()`, `_build_duplicate_feedback()` |

## Metrics and Stop Conditions

### Quality Assessment

Metrics are evaluated against YAML-defined thresholds:

| Quality | R-free (2.0Å) | Map CC | Anomalous Measurability |
|---------|---------------|--------|-------------------------|
| Good | < 0.25 | > 0.80 | > 0.10 |
| Acceptable | < 0.30 | > 0.70 | > 0.05 |
| Poor | ≥ 0.30 | ≤ 0.70 | ≤ 0.05 |

### Stop Conditions

The workflow stops when:

1. **Success**: Quality metrics reach targets AND validation done
2. **Plateau**: < 0.5% improvement for 2+ consecutive cycles
3. **Hopeless bailout**: R-free > 0.50 after 1+ refinement cycle, but only if `autobuild_done` (v115.05 — gives autobuild a chance to rebuild missing density before declaring the structure unsolvable)
4. **Hard limit**: 3+ refinement cycles regardless of R-free
5. **Validation Gate**: Must run molprobity before stopping if R-free is good
6. **Clashscore shortcut**: clashscore < 10 with reasonable R-free, but only if `refine_count >= 1` (v115.05 — prevents declaring at-target before any refinement has run)

### Refinement Loop Enforcement

When `_is_at_target()` returns True (conditions 2-4, 6 above),
`get_valid_programs()` actively **removes** `phenix.refine` and
`phenix.real_space_refine` from valid programs in both `validate` and
`refine` phases, and adds `STOP`. This prevents the LLM from selecting
refinement even when it appears as a phase-preferred program.

**Exception:** `needs_post_ligandfit_refine` always allows refinement.
After ligand fitting changes the model, re-refinement is scientifically
required regardless of the current cycle count or R-free value.

The validation gate prevents stopping without validation: if R-free is below the success threshold or 3+ refinement cycles have completed, STOP is removed from valid_programs until validation runs.

## Client-Server Update Model

Understanding which code runs where is essential for planning fixes and knowing
whether a deployed fix reaches existing users without requiring them to reinstall.

### The Core Principle

**Decisions and knowledge → server. Execution and I/O → client.**

The client's job is: receive user input → serialize it → send to server → receive
a command string → run that command string locally. Everything in between happens
on the server.

### Local/Remote Parity Invariant

**LocalAgent MUST produce identical results to RemoteAgent.** Both agents go
through the exact same request preparation pipeline:

```
_query_agent_for_command() — identical for both modes
  → build_session_state(session_info)
  → build_request_v2(files, history, session_state, ...)
  → request["settings"][...] = ...           ← same settings in both
  → prepare_request_for_transport(request)   ← same encoding
  → [LocalAgent: decode locally | RemoteAgent: send over HTTP]
  → run_ai_agent.run(request_json)           ← identical server-side code
  → group_args(history_record, events, ...)  ← same return format
```

The `LocalAgent` intentionally performs the full encode/decode roundtrip
(not a shortcut). This ensures transport bugs are caught during local testing
and that local mode produces byte-identical requests to what the server
would receive.

**When adding a new field or setting, you MUST update both agents:**
- `phenix_ai/local_agent.py` — `decide_next_step()`
- `phenix_ai/remote_agent.py` — `decide_next_step()`

If the agents diverge, users will get different results depending on whether
`run_on_server` is True or False. This is a silent correctness bug that is
extremely difficult to diagnose.

### Execution split diagram

```
┌──────────────────────────────────────────────────────────────────┐
│                        CLIENT (user install)                      │
│                                                                   │
│  run()  ──→  run_job_on_server_or_locally()                      │
│                     │                                             │
│                     ├── run_job_on_server()  ─── serialize ──────┐│
│                     │   _inject_user_params()  ← RUNS CLIENT-SIDE││
│                     │   _run_single_cycle()    ← RUNS CLIENT-SIDE││
│                     │   phenix program execution                  ││
│                     │   GUI components (wxGUI2/)                  ││
│                     │   .eff / .pkl file generation               ││
│                     │                                             ││
│                     └── run_job_locally()  ─── (has local DB) ───┘│
└────────────────────────────────────────────────────────────────── │
                                           HTTP / REST               │
                                               ▼                     │
┌──────────────────────────────────────────────────────────────────┐│
│                        SERVER (your install)                      ││
│                                                                   ││
│  run_ai_agent.run()                                               ││
│     └── LangGraph pipeline                                        ││
│          ├── PERCEIVE  (graph_nodes.py)                           ││
│          ├── PLAN      (graph_nodes.py + prompts_hybrid.py)       ││
│          ├── BUILD     (command_builder.py + program_registry.py) ││
│          ├── VALIDATE  (graph_nodes.py)                           ││
│          └── OUTPUT    (graph_nodes.py)                           ││
│                                                                   ││
│  Knowledge layer                                                  ││
│     ├── knowledge/programs.yaml                                   ││
│     ├── knowledge/workflows.yaml                                  ││
│     ├── knowledge/prompts_hybrid.py                               ││
│     └── agent/workflow_engine.py, workflow_state.py, session.py  ││
└───────────────────────────────────────────────────────────────────┘│
```

### Analysis mode routing (`ai_analysis.py`)

The `ai_analysis.py` module is shared infrastructure for all LLM
interactions. It has five `analysis_mode` values, but only one of
them uses the RAG database:

| Mode | Used by | Needs RAG DB | What it does |
|------|---------|:------------:|--------------|
| `standard` | `phenix.ai_analysis` (standalone) | **Yes** | Analyzes a single PHENIX log file using retrieval-augmented generation against the Phenix knowledge base (documentation, papers, newsletters). Produces a summary + detailed analysis. This is the original `phenix.ai_analysis` program; **the AI Agent does not use this mode**. |
| `directive_extraction` | AI Agent (session start) | No | Parses user advice into structured directives (prefer_programs, after_program, strategy settings). Pure LLM call. |
| `advice_preprocessing` | AI Agent (session start) | No | Reformats a tutorial README into structured guidance the agent can follow. Pure LLM call. |
| `failure_diagnosis` | AI Agent (on terminal error) | No | Produces a three-section diagnosis (what went wrong / cause / fix) from a program's error log. Pure LLM call. |
| `agent_session` | AI Agent (session end) | No | Generates an end-of-run assessment from the session history. Pure LLM call. |

**Why this matters for routing.** The v115.05 `_LLM_ONLY_MODES` set
in `run_job_on_server_or_locally()` routes the four non-standard
modes to local execution when `run_on_server=False` or when
`provider=ollama`. Only `standard` mode is sent to the server,
because it is the only mode that needs the RAG database. Without
this distinction, all five modes would queue at the Phenix server
for every cycle, even when a local LLM is available.

**The agent's THINK node is separate.** The THINK node in the
LangGraph pipeline does its own log analysis using `thinking_prompts.py`
and the expert knowledge base — it does not go through `ai_analysis.py`
at all. The two systems share the same LLM providers but serve
different purposes: THINK produces per-cycle expert assessments that
feed PLAN, while `ai_analysis.py` standard mode produces standalone
summaries for human consumption.

### Observability across the client/server split (v119.H4 lesson)

The local/remote execution split has a non-obvious consequence
for **observability**: stderr and other auxiliary output streams
behave differently depending on which side of the split a piece of
code is executing on. Instrumentation that "works locally" can
silently disappear when the same code runs remotely.

**Local mode** (`run_on_server=False` or `provider=ollama`):
`_run_advice_preprocessing_locally`, `_run_directive_extraction_locally`,
and similar helpers in `ai_analysis.py` call into
`phenix_ai/run_ai_analysis.py` in the **client process**. Anything
those functions write to `sys.stderr` (e.g., diagnostic markers,
exception tracebacks, fault-isolation fallbacks) lands in the
client's terminal stderr, the same stream as the rest of
`phenix.ai_agent`'s output. With `phenix.ai_agent ... 2>&1 | tee
output.log` (or any shell-level capture), the marker is preserved.

**Remote mode**: the same code runs inside the
`ai.phenix-online.org` server process. Its `sys.stderr` is the
server's stderr, not the client's. Successful results are returned
to the client as a structured dict; the server's stderr is **not**
included in that dict. Exceptions are an exception (pun intended):
when the server-side call raises, the traceback IS relayed back as
a string inside the result, which is how the v118 PHIL parse error
became visible to the client in the first place. But normal
diagnostic output from a successful server-side call has no path
to the client log.

**Implications for new instrumentation:**

1. **Server-side markers in normal control flow are invisible to
   the client** unless explicitly threaded through the result dict.
   `sys.stderr.write("[MY_MARKER] ...")` works under
   `run_on_server=False` but disappears under remote mode.
2. **For telemetry that must survive both modes**, the canonical
   pattern is to (a) write to stderr as before, AND (b) include the
   marker value as a field in the `group_args` result so the client
   receives it and can re-emit it on the client side.
3. **For server-only telemetry**, write to a known file path on
   the server host (e.g., a dedicated log file). The operator who
   maintains the server collects this directly; the client doesn't
   need to see it.
4. **For local-only debugging**, plain `sys.stderr.write` or
   `print` is fine. But remember to `force run_on_server=False`
   explicitly when testing — otherwise the run may silently go
   remote and the markers will not appear locally regardless of
   how the code is patched.

**Concrete example (v119.H4 Step 1F):** The Step 1F metric block
in `run_advice_preprocessing` writes its `[STEP_1F]` marker to
`sys.stderr`. This was diagnosed during H4 deployment when a run
that appeared to be local was actually being dispatched to the
remote server, and the markers (correctly emitted on the server's
stderr) were never seen by the client. The marker shows correctly
when `run_on_server=False` is set in PHIL. For broader-coverage
telemetry collection across both modes, future iterations may need
to thread the metric tuple through the result dict (option 2
above) or write to a dedicated server-side log file (option 3).

This concerns observability only — it does not affect the **Local/
Remote Parity Invariant** above. Both modes execute the same code
on identical input and produce identical output via the documented
return contract. The difference is in *how visible* internal
runtime details are to the operator running the client.

### Always server-side (no user action needed)

**A user with yesterday's PHENIX install must be able to connect to
today's server and get correct results.** All of the following changes
are server-side-only and take effect immediately for all users — no
reinstall required. They are safe because the client never reads these
files; it sends data and receives a command string.

If a change requires the client to send new data or handle a new response
format, that is a **protocol change** — see `DEVELOPER_GUIDE.md §8
(Backward Compatibility & Contract)` for the version-bump procedure.

**LLM decision-making — the entire graph**
Any change to how the agent thinks: PERCEIVE, PLAN, BUILD, VALIDATE, ACT nodes,
workflow routing, phase detection, placement logic, error recovery. Users get
this immediately.

**Prompts and knowledge**
`prompts_hybrid.py`, `programs.yaml`, `workflows.yaml`. Any improvement to how
programs are described to the LLM, new invariants, strategy flag fixes, or stop
conditions — effective immediately.

**Command construction**
`command_builder.py`, `program_registry.py`. Changes to passthrough filtering,
invariant application, file selection — all server-side.

**Session and history logic**
`session.py`, `workflow_engine.py`, `workflow_state.py`. How history is
interpreted, context built, S2c promotion, zombie detection — all server-side.

**API schema and transport encoding**
Changes to the request/response format take effect on the server. The client
passes opaque encoded strings through without inspecting them.

### Always client-side (requires user to update their install)

**Post-command injection — `_inject_user_params()`**
This function runs *after* the server returns a command, on the user's machine.
It can corrupt an otherwise-correct server-built command. This is the most
dangerous client-side code path because server fixes alone cannot protect against it.

**The local execution loop**
`_run_single_cycle()`, `_get_command()`, the `for cycle in range(...)` loop,
`iterate_agent`. These control how commands get run, retried, and logged locally.

**Phenix program execution**
`_try_native_execution()`, `_run_easy_run()`, subprocess handling, output
capture. The server does not run phenix programs; the client does.

**GUI components**
`wxGUI2/Programs/DockInMap.py` and all other `wxGUI2` files. Restoration,
display widgets, and result panels are entirely client-side.

**The .eff and .pkl generation**
`generate_program_eff()`, DataManager PHIL mapping. File writing happens on the
client.

**The top-level branching**
`run()`, `run_job_on_server_or_locally()`, `run_job_on_server()`. The logic that
decides whether to call the server at all. A bug here means users might not
reach the server.

### The gray area: `programs/ai_agent.py`

This single file straddles both worlds. The **top half** (through
`run_job_on_server()`) runs client-side. The **bottom half** (when the server
receives `run_on_server=False`) runs server-side. Changes in this file need user
updates if in the client path, and don't if in the server path.

**Practical rule:** if logic can live in `agent/`, `knowledge/`, or
`phenix_ai/run_ai_agent.py` rather than in the client path of `ai_agent.py`,
prefer that location — it turns a required user update into a free server update.
The `_inject_user_params` bug (S2k) is a textbook example: if post-injection
filtering were done server-side before returning the command string, no user
update would ever be needed for that class of bug.

### Practical decision rule for new fixes

When writing a fix, ask: *does this code run before or after
`run_job_on_server()` is called?*

- **Before** (serialization, call site, result handling) → client-side → user must update
- **After** (on the machine that received `run_on_server=False`) → server-side → free

### Summary table

| Component | Runs on | Update required? |
|---|---|---|
| `knowledge/prompts_hybrid.py` | Server | No |
| `knowledge/programs.yaml` | Server | No |
| `knowledge/workflows.yaml` | Server | No |
| `agent/graph_nodes.py` | Server | No |
| `agent/command_builder.py` | Server | No |
| `agent/program_registry.py` | Server | No |
| `agent/workflow_engine.py` | Server | No |
| `agent/workflow_state.py` | Server | No |
| `agent/session.py` | Server | No |
| `phenix_ai/run_ai_agent.py` | Server | No |
| `programs/ai_agent.py` — `_inject_user_params()` | **Client** | **Yes** |
| `programs/ai_agent.py` — `_run_single_cycle()` | **Client** | **Yes** |
| `programs/ai_agent.py` — `run_job_on_server()` | **Client** | **Yes** |
| `programs/ai_agent.py` — `run()` top-level | **Client** | **Yes** |
| `wxGUI2/Programs/DockInMap.py` | **Client** | **Yes** |
| `phenix_ai/remote_agent.py` | **Client** | **Yes** |

---

## Extension Points

### Adding a New Program

See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide. In summary:

1. Add program definition to `knowledge/programs.yaml` (inputs, outputs, log_parsing)
2. Add to appropriate workflow phase in `knowledge/workflows.yaml`
3. Add file categories to `knowledge/file_categories.yaml` (if new file types)
4. Add hardcoded extractor to `phenix_ai/log_parsers.py` (only if YAML patterns insufficient)

### Adding a New Workflow

1. Define states in `workflow_state.py`
2. Add transitions in workflow YAML
3. Update rules selector if needed

### Adding New Session State

1. Add to `session_state` in request schema
2. Map to `session_info` in `_process_request()`
3. Use in graph nodes as needed

## Performance Considerations

### Graph Caching

```python
_cached_graph = None

def _get_graph():
    global _cached_graph
    if _cached_graph is None:
        _cached_graph = build_agent_graph()
    return _cached_graph
```

### Request Size

- Log content can be large (truncate if needed)
- History grows with cycles (summarize old cycles)
- Files list should only include relevant files

### Response Time

- Graph invocation: ~1-5 seconds
- LLM calls: ~2-10 seconds (if used)
- Rules-only mode: <1 second

## Command Builder Architecture

The CommandBuilder provides a unified, single-entry-point system for generating
PHENIX commands. It replaces the previous fragmented approach spread across
multiple modules.

### Pipeline

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         CommandBuilder.build()                               │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────┐  ┌──────────┐  ┌───────────┐  ┌──────────┐  ┌─────────────┐  │
│  │_select_  │─▶│_build_   │─▶│_apply_    │─▶│_assemble_│─▶│ _inject_    │  │
│  │ files()  │  │strategy()│  │invariants()│  │command() │  │ recovery()  │  │
│  └──────────┘  └──────────┘  └───────────┘  └──────────┘  └─────────────┘  │
│       │              │              │              │              │          │
│       ▼              ▼              ▼              ▼              ▼          │
│  Priority order:  Auto-fill:    Auto-fill:    Template-      Append any     │
│  1. LLM hints    - output_pfx  - resolution   based cmd      recovery-     │
│  2. Locked rfree - from hist   - R-free flags  string         sourced       │
│  3. Best files                 - twin_law                     strategy      │
│  4. Categories                                                entries not   │
│  5. Extensions                                                in command    │
│                                                                              │
│  ↕ _apply_recovery_strategies() runs between _build_strategy and            │
│    _apply_invariants, adding recovery flags to the strategy dict.           │
│    Recovery entries are tagged with strategy_sources[key]="recovery".        │
└──────────────────────────────────────────────────────────────────────────────┘
```

The recovery injection step is necessary because `registry.build_command()` only
emits strategy entries matching `strategy_flags` in programs.yaml.  Recovery
params use fully-qualified PHIL paths (e.g., `scaling.input.xray_data.obs_labels`)
that don't match short strategy_flags names.  The post-assembly injection ensures
they reach the final command.  See "Automatic Error Recovery" section for details.

### CommandContext

All parameters are passed via a `CommandContext` dataclass:

```python
@dataclass
class CommandContext:
    cycle_number: int = 1
    experiment_type: str = ""
    resolution: float = None
    best_files: Dict[str, str] = None
    rfree_mtz: str = None
    categorized_files: Dict[str, List[str]] = None
    workflow_state: str = ""
    history: List[Dict] = None
    llm_files: Dict[str, str] = None      # Optional LLM hints
    llm_strategy: Dict[str, Any] = None   # Optional LLM hints
    log: Callable = None
```

### File Selection Priority

1. **LLM hints** - If LLM specified files and they exist
2. **Locked R-free MTZ** - For MTZ slots in X-ray refinement
3. **Best files** - From BestFilesTracker (quality-scored)
4. **Categories** - From `input_priorities` in YAML
5. **Extension fallback** - Most recent file with matching extension

**Best files exclusion:** Before applying `best_files["model"]` as a model
override, the system checks if the file is in a program's `exclude_categories`.
E.g., `ligand_fit_1.pdb` (in `ligand_fit_output` → `ligand`) is excluded from
refine's model slot, falling through to category-based selection.

### Command Validation

The `validate` graph node checks that all file references in a generated command
exist in `available_files`. Output arguments (`output.file_name=X.pdb`,
`output.prefix=Y`) are stripped before extraction since they reference files
that don't exist yet. Duplicate commands and empty commands are also rejected.

### Integration with Graph Nodes

```python
# In graph_nodes.py build() function:
if USE_NEW_COMMAND_BUILDER:
    return _build_with_new_builder(state)

# _build_with_new_builder creates CommandContext from state
# and delegates to CommandBuilder.build()
```

## Automatic Error Recovery

When a PHENIX program fails with a recognized error pattern, the agent can automatically
recover and retry with corrected parameters. This is handled by `agent/error_analyzer.py`
with patterns defined in `knowledge/recoverable_errors.yaml`.

### How It Works

1. **Detection**: After a program fails, the error output is matched against
   `detection_patterns` in `recoverable_errors.yaml`
2. **Extraction**: Regex patterns extract the parameter name and available choices
   from the error message
3. **Resolution**: A strategy selects the best choice based on program type and
   workflow context (e.g., anomalous vs merged data)
4. **Storage**: The resolution is saved as a file-keyed recovery strategy in the
   session (`session.set_recovery_strategy`), and `force_retry_program` is set
5. **Retry**: On the next cycle, `force_retry` bypasses the LLM's program choice,
   the command builder applies recovery flags, and the command is re-executed

### Recovery Parameter Injection (v112.74)

Recovery strategies add fully-qualified PHIL parameters to the command builder's
strategy dict (e.g., `scaling.input.xray_data.obs_labels=I(+)` for ambiguous
data labels).  These must survive two stages to reach the final command:

**Stage 1: build_command template expansion.**  `registry.build_command()` only
emits strategy entries matching `strategy_flags` keys in programs.yaml.  Recovery
params use fully-qualified PHIL paths that don't match short flag names.  Fix:
after `_assemble_command` returns, any strategy entry sourced from recovery
(`strategy_sources[key] == "recovery"`) whose key is not already present in the
command is appended as `key=value`.

**Stage 2: probe-only sanitizer.**  Programs in `_PROBE_ONLY_FILE_PROGRAMS`
(xtriage, model_vs_data, etc.) have all non-file-path key=value tokens stripped.
Data-label selection parameters (`obs_labels`, `labels`, `data_labels`,
`anomalous_labels`, `r_free_flags_labels`) are whitelisted because they come
exclusively from error recovery for ambiguous MTZ arrays.

### Recovery Loop Guard (v112.74)

Before storing a new recovery strategy, the system checks if one already exists
for the same file (`session.get_recovery_strategy(file)`).  If so, the previous
recovery attempt didn't work — re-triggering would create an infinite loop.
The recovery is skipped and the failure falls through to the terminal diagnosis
path.

### Duplicate Check Bypass for Recovery (v112.74)

Recovery retries re-run the same program, often with commands that look like
duplicates of the failed command.  When `forced_retry` is detected in
`decision_info` (or fallback: reasoning text contains recovery marker), the
duplicate check in `_handle_duplicate_check` is skipped entirely.  The
`forced_retry` flag propagates from `state['intent']` through the BUILD node's
return dict to the top-level graph output.

### Currently Handled Errors

- **`ambiguous_data_labels`**: MTZ file contains multiple data arrays (e.g., both
  merged intensities and anomalous pairs). The agent picks the right array based on
  whether the program needs anomalous or merged data.  Resolution kind:
  `select_value`.
- **`ambiguous_experimental_phases`**: MTZ file contains multiple phase arrays.
  The agent selects standard vs anomalous HL coefficients based on context.
  Resolution kind: `select_value`.
- **`missing_phib_input_map_file`** (v119.H17): autobuild was called with
  `map_file=<raw MTZ>` but the MTZ has no PHIB phase columns.  The agent strips
  the `map_file=` / `input_map_file=` / `input_files.map_file=` flags from the
  retry command.  Resolution kind: `strip_parameter`.
- **`rfree_flags_mismatch`** (declared pre-H17; resolved retroactively by H17's
  generic `_resolve_strip_parameter` handler): refine was called with a
  `rfree_mtz=` that has incompatible R-free flags.  The agent strips the flag
  and lets refine generate new flags.  Resolution kind: `strip_parameter`.
- **`rfree_flags_missing`** (v120): refine was given reflection data with **no**
  R-free flag array and neither a flags-bearing MTZ nor
  `xray_data.r_free_flags.generate=True`, so it aborts with "No array of R-free
  flags found".  Observed in AF-bromodomain-ligand (AIAgent_15): refine ran on
  the raw data MTZ directly, failed, and the agent advanced to polder on the
  unrefined protein+ligand model.  The agent forces a retry of `phenix.refine`;
  on the rebuilt retry the server-side BUILD node / `command_builder` re-selects
  the locked R-free MTZ (`context.rfree_mtz`) for the data slot — preserving
  R-free flag continuity (no test-set leakage) — or falls back to
  `generate=True` when no MTZ is locked.  This is the **inverse twin** of
  `rfree_flags_mismatch` (which strips `generate=True` when flags already
  exist).  Resolution kind: `force_retry`.

#### R-free generate guard (v120.2): per-file, parity-safe generate decision

The three R-free invariants the agent must honor are: (1) the data MTZ being
refined has **no** test set → generate one; (2) it **already has** a test set →
keep it, never generate (overwrite); (3) once a set is generated/locked → never
regenerate.

`command_builder.py`'s phenix.refine invariant originally added
`xray_data.r_free_flags.generate=True` whenever `context.rfree_mtz` was not
locked.  But `rfree_mtz` is only locked from a refinement **output**
(`session.set_rfree_mtz`, called post-cycle), never from the **input**.  So a
*first* refinement against pre-flagged input data had no lock yet and wrongly
received `generate=True`, which phenix.refine treats as "create a **new random**
R-free partition", **overwriting** existing flags and silently invalidating
cross-validation (observed on `beta_blip_001.mtz`).  Conversely, in the MR
workflow the first refinement runs against a phaser **output** (`PHASER.1.mtz`)
that has **no** flags, so it *must* generate — omitting it aborts phenix.refine
with `No array of R-free flags found` (observed in AIAgent_35).  The two cases
need **opposite** decisions and both are "first refinement, nothing locked", so
the decision cannot key on the lock alone — it must know whether **the specific
file being refined** has flags.

**The MTZ is examined — once, on the CLIENT.**  The server cannot read client
file paths (`files_local=False`), so the build cannot inspect the data MTZ itself
without breaking server/local parity.  Instead `ai_agent.py` inspects every local
MTZ it can see each cycle — original inputs **and** intermediate outputs like
`PHASER.1.mtz` — via `mtz_inspector.inspect_mtz` (which opens the file with
iotbx and reports the R-free column), and ships a per-file map
`mtz_rfree_map = {basename: bool}` in `session_info`.  (A scalar
`input_mtz_has_rfree`, computed from the original input only, is also shipped as a
fallback; it predates the map and is insufficient on its own because the refined
file is often a derived output, not the original input.)

**The guard never reads the filesystem.**
`CommandBuilder._input_mtz_rfree_state(files, context)` returns a **tri-state**
(`True` / `False` / `None`) from sources that are **all `session_info`-borne**, so
the answer is byte-identical on server and local:
(0) `context.mtz_rfree_map[basename(selected data_mtz)]` — the per-file answer for
the file actually being refined (authoritative);
(1) `context.input_mtz_has_rfree` — the original-input scalar, used only when the
selected file is absent from the map;
(2) `context.mtz_inspection` — precomputed inspection, positive only.
There is **no** `files_local`-gated local read (a v120 source that made the local
build resolve where the server could not — a parity violation; removed in v120.2).

The generate decision (F1, v120.3) keys on the **selected file's** tri-state, with
the lock only as the undetermined fallback — it no longer strips merely because a
lock exists:

- `state is True` → **strip** (positive evidence the selected file has flags);
- `state is None` (undetermined) **and `rfree_mtz` locked** → **strip** (fallback:
  assume the locked dataset's flags; note F2 reconciliation normally resets the
  slot to the lock first, so the selected file usually *is* the flagged lock);
- `state is False` → **generate** (positive evidence it has none);
- `state is None` (undetermined) **and nothing locked** → **generate**.

The last rule (v120.1) restores the robust pre-fix first-refinement default:
when nothing is locked and we cannot prove flags exist, generate.  We **strip
only on positive evidence**, which is what prevents the overwrite; we never omit
generate on an unlocked refine, which is what keeps the MR workflow running even
if the LLM forgets the flag.  BUILD log lines: `BUILD: Stripped
generate_rfree_flags (…)` on the strip path, and `BUILD: First refinement - will
generate R-free flags (…)` on the generate path.

**Plumbing & protocol.**  Both fields travel the shared transport
(`build_session_state` → `build_request_v2` whitelist → `run_ai_agent`
map-back), each with an explicit `is not None` guard so a confirmed-`False`
scalar and an empty map are not dropped.  The live build path
(`graph_nodes._build_with_new_builder`, `USE_NEW_COMMAND_BUILDER=True`)
constructs `CommandContext` directly and passes both fields into all three of its
`CommandContext(...)` constructors.  `mtz_rfree_map` and `input_mtz_has_rfree` are
registered in `contract.py` (protocol v8).  Because the decision depends on no
filesystem state, LocalAgent and RemoteAgent produce identical commands by
construction.

**Parity caveat (dormant dead path).**  `CommandContext.from_state` still calls
`inspect_mtz` to populate `mtz_inspection` — a server-side file read that
contradicts the parity model.  It is currently dormant: `from_state` is NOT on the
live build path (`_build_with_new_builder` builds the context directly), surviving
only in docstrings and the `__main__` demo, and even if reached it degrades to
`None` (inspect_mtz returns `None` on a missing path).  It is a cosmetic
inconsistency, not a live bug, and should be removed in a standalone cleanup so the
codebase has a single rule: the guard never reads files; the client ships the
answer.

The parallel `generate_rfree_flags` injection in `graph_nodes.py` (~4086) is dead
in the default path (guarded `USE_NEW_COMMAND_BUILDER is False`); `command_builder`
is authoritative.  (A separate `validation/` `ProgramValidator` package once
contained an independent R-free injection path — `RefineValidator.prevalidate` →
`add_rfree_generation_if_needed`, which shelled out to `phenix.mtz.dump` and
carried the same server-blindness latent bug — but it was ORPHANED, with no caller
in the live agent path, and was **removed as dead code in v120.2**; see CHANGELOG.)
Pinned by `tests/tst_rfree_generate_guard.py` (16 cases: the three requirements,
LLM-strip, the per-file map cases, map-overrides-scalar, the AIAgent_35 MR
regression, the undetermined-unlocked-generates default, and a parity assertion
that the guard performs no file reads) and
`tests/tst_input_mtz_has_rfree_plumbing.py` (the scalar True/False/None plus the
`mtz_rfree_map` wire round-trip, proving `False` and the map survive intact).

### Resolution kinds (v119.H17, extended v120)

The agent supports these resolution kinds:

- **`select_value`**: disambiguates an enum the LLM left ambiguous (label
  selection, array choice).  Picks a default per a YAML-defined preference policy
  and injects via the recovery parameter mechanism.
- **`strip_parameter`** (v119.H17): removes an inappropriate flag entirely on
  retry.  Detection runs in `error_analyzer.py`; execution wires through
  `programs/ai_agent.py::_execute_command` which applies a robust regex
  (handles PHIL spacing, quoted-with-spaces, single quotes, end-of-line) to
  strip the flag prefix before running the retry command.  Emits `[STRIP]` log
  line.  One-shot via `pop()` — entry is consumed on use so a subsequent retry
  doesn't re-strip.
- **`force_retry`** (v120): re-runs the program with NO flag additions and NO
  flag strips — the fix is delegated entirely to the command rebuild.  The
  forced retry routes through the server-side graph (PLAN emits a retry intent
  with empty `files={}` / `strategy={}`), so `command_builder` does its own file
  selection on the rebuilt command.  Used by `rfree_flags_missing`: the builder
  substitutes the locked R-free MTZ (or falls back to `generate=True`) without
  the recovery having to inject any parameter.  `error_analyzer.py` provides
  `_resolve_force_retry` (returns an `ErrorRecovery` with empty `flags` AND empty
  `strip_flags`, carrying only `retry_program`) plus a `force_retry` marker in
  `_extract_error_info` so `analyze()` doesn't bail before dispatch.  No
  `ai_agent.py` change was needed: `_handle_recovery` already calls
  `set_force_retry_program` for every recovery.

> Note: an `add_parameter` kind (inject a missing flag such as
> `r_free_flags.generate=True`) is *not* implemented.  The
> `rfree_flags_missing` case that would have used it is instead handled by
> `force_retry`, which lets `command_builder`'s existing R-free logic choose
> between the locked MTZ and `generate=True` — covering both situations with one
> entry and no new injection path.

### Configuration

`knowledge/recoverable_errors.yaml` defines:
- **`errors`**: Detection patterns, extraction regexes, resolution strategies, max retries
- **`label_patterns`**: How to classify data labels as anomalous or merged
- **`program_data_preferences`**: Which programs need anomalous vs merged data
- **`context_keywords`**: Workflow hints from user advice (e.g., "SAD" → anomalous)

### Adding a New Recoverable Error

1. Add an entry under `errors:` in `recoverable_errors.yaml` with detection patterns
   and extraction regexes.  Choose a resolution kind: `select_value`,
   `strip_parameter` (v119.H17), or `force_retry` (v120).
2. For `select_value` resolutions: add a resolution handler in
   `error_analyzer.py` (`resolve_error()` method).
3. For `strip_parameter` resolutions: the generic `_resolve_strip_parameter`
   handler reads the `strip_parameters` YAML list directly — no per-error
   handler needed.
4. For `force_retry` resolutions: the generic `_resolve_force_retry` handler
   needs only `resolution: force_retry` (+ optional `retry_program`) in YAML.
   IMPORTANT: also add a marker branch in `_extract_error_info`
   (`return {"resolution": "<kind>"}`) for any new non-extracting kind, or
   `analyze()` returns None at its `if not error_info` guard before dispatch.
5. The fallback node in the graph will automatically use the new pattern.
6. K-test contract: for `strip_parameter` entries, pin the regex against the
   PHIL-spacing/quote variants — see `tests/tst_h17_strip_executor.py` for the
   reference template.  For `force_retry` entries, assert the full `analyze()`
   path returns an `ErrorRecovery` with the right `retry_program` and empty
   `flags`/`strip_flags` — see `tests/tst_rfree_flags_missing.py`.

### Recovery strategies persist intentionally

After a program succeeds with a selected label (e.g., `obs_labels=I(+)`), the
recovery strategy for that MTZ file is NOT cleared.  The ambiguity is a property
of the MTZ file, not the program.  Downstream programs (phaser, refine) that use
the same MTZ also need the label selection.

### Key Files

- `knowledge/recoverable_errors.yaml` — Error patterns and resolution config
- `agent/error_analyzer.py` — Detection, extraction, and resolution logic
- `agent/graph_nodes.py` — Fallback node triggers error analysis on failure
- `agent/command_builder.py` — Recovery param injection after `_assemble_command`
- `agent/command_postprocessor.py` — Label param whitelist in probe-only sanitizer
- `programs/ai_agent.py` — Loop guard, duplicate check bypass, `force_retry` handling

---

## Diagnosable Terminal Errors

Some program failures are categorically unrecoverable — for example, a crystal symmetry
mismatch between input files, or a model entirely outside the cryo-EM map. Rather than
stopping silently or showing a generic error, the agent detects these, asks the LLM to
diagnose the specific cause, and presents a self-contained HTML report to the user.

### Overview

When a recognized terminal failure occurs the agent:

1. Identifies the error type via `DiagnosisDetector` (pattern matching against
   `knowledge/diagnosable_errors.yaml`)
2. Calls the LLM on the server (`analysis_mode=failure_diagnosis`) with the error
   excerpt, a domain hint from the YAML, and the tail of the failing program's log
3. Writes `ai_failure_diagnosis.html` to the job's working directory
4. Opens the report in the user's browser automatically
5. Returns `True` from `_run_single_cycle` to stop the cycle loop cleanly —
   **no `Sorry` exception is raised**; the ai_agent job is considered successful
   because it correctly identified and diagnosed the sub-job failure
6. `_finalize_session` runs unconditionally (saves the session, populates
   `self.result`), but **skips the Results summary page** — the diagnosis HTML
   is the user's sole output window, so a second page would only bury it

### Error flow diagram

```
Program fails
      │
      ▼
DiagnosisDetector.detect(result_text)
      │
      ├── No match → normal failure handling (agent may retry or move on)
      │
      └── Match found: (error_type, description, excerpt)
                │
                ▼
       _diagnose_terminal_failure()
                │
                ├── 1. Read log tail (client-side, 150 lines / 4 000 chars)
                ├── 2. LLM call: analysis_mode=failure_diagnosis
                │         └── fallback: rules-only text if LLM unavailable
                ├── 3. Build & write ai_failure_diagnosis.html
                │         ├── Heading: "Error diagnosis"
                │         ├── Meta bar: program, cycle, job name, working dir
                │         ├── Error excerpt (red box)
                │         ├── AI diagnosis (three paragraphs: what/cause/fix)
                │         └── Footer: full path to the saved file
                ├── 4. Open HTML in browser (load_url)
                ├── 5. Store path in session.data["failure_diagnosis_path"]
                ├── 6. Print diagnosis to CLI log
                └── 7. return True  ← cycle loop breaks; no Sorry raised

_finalize_session()
      │
      ├── Saves session JSON
      ├── Detects failure_diagnosis_path in session.data
      └── Skips _generate_ai_summary() and display_results()
              (Results page suppressed; diagnosis page is the output)
```

### HTML report content

The report (`ai_failure_diagnosis.html`) is self-contained — no CDN dependencies —
so it renders correctly when opened directly from a local path. It includes:

| Section | Content |
|---------|---------|
| Header | "⚠ PHENIX AI Agent — Error diagnosis" |
| Meta bar | Failed program, cycle number, job name, working directory |
| Error | Human-readable error type label + the raw error excerpt |
| AI Diagnosis | Three plain-text paragraphs: WHAT WENT WRONG / MOST LIKELY CAUSE / HOW TO FIX IT |
| Footer | Full path where the file is saved |

### Configuration

`knowledge/diagnosable_errors.yaml` defines each recognizable error type:

```yaml
errors:
  crystal_symmetry_mismatch:
    display_name: "Unit cell or space group mismatch between input files"
    detection_patterns:
      - "crystal symmetry mismatch"
      - "incompatible crystal symmetry"
    diagnosis_hint: |
      Check that all input files were processed from the same crystal form.
      ...
```

Fields:

| Field | Purpose |
|-------|---------|
| `display_name` | Human-readable label shown on the error page and in the log |
| `detection_patterns` | Strings that, if found in the failing program's output, trigger diagnosis |
| `diagnosis_hint` | Domain knowledge injected into the LLM prompt to guide the diagnosis |

### Adding a new diagnosable error

1. Add an entry in `knowledge/diagnosable_errors.yaml` with `display_name`,
   `detection_patterns`, and a `diagnosis_hint`
2. No code changes required — `DiagnosisDetector` and `_diagnose_terminal_failure`
   handle all new entries automatically
3. Add a detection test in `tests/tst_audit_fixes.py` following the `test_s3a_detect_*`
   pattern

### Key files

| File | Role |
|------|------|
| `knowledge/diagnosable_errors.yaml` | Error definitions and LLM hints |
| `agent/error_analyzer.py` | `DiagnosisDetector` — pattern matching and hint lookup |
| `agent/failure_diagnoser.py` | Prompt builder, Markdown sanitiser, HTML report builder |
| `programs/ai_agent.py` | `_diagnose_terminal_failure()` — orchestrates steps 1–7 |
| `tests/tst_audit_fixes.py` | `test_s3a_*` tests covering detection, HTML output, UX flow |

---

## Session Summary Generation

The agent generates structured summaries of completed sessions with optional LLM assessment.

### Summary Components

1. **Header**: Run name derived from the working directory basename
2. **Session line**: Session ID, cycle count (successful / total), and **working directory path**
3. **Input Section**: Files, user advice, experiment type, resolution
4. **Input Data Quality**: Metrics from xtriage/mtriage (resolution, completeness, twinning)
5. **Workflow Path**: High-level description of strategy taken
6. **Steps Performed**: Table of all cycles with key metrics
7. **Final Quality**: Final metrics with quality assessments
8. **Key Output Files**: Output files from the session with relative paths
9. **Failure Diagnosis** *(only on fatal error)*: Path to `ai_failure_diagnosis.html` — shown when the session ended due to a diagnosable terminal failure; see [Diagnosable Terminal Errors](#diagnosable-terminal-errors)
10. **Assessment** (optional): LLM-generated evaluation

> **Note**: When a fatal error diagnosis was produced, the Results summary page is suppressed entirely — `_finalize_session` detects `session.data["failure_diagnosis_path"]` and skips `_generate_ai_summary()`. The diagnosis HTML is the user's sole output window.

### Quality Assessments

Metrics are automatically assessed based on resolution-dependent thresholds:

| Metric | Good | Acceptable | Needs Improvement |
|--------|------|------------|-------------------|
| R-free (2.5Å) | ≤0.25 | ≤0.30 | >0.30 |
| Map CC | ≥0.80 | ≥0.70 | <0.70 |
| Clashscore | ≤5 | ≤10 | >20 |

### LLM Assessment Prompt

When enabled, the LLM evaluates:
1. **Input Data Quality**: Resolution, completeness, issues
2. **Goal and Strategy**: User's goal and agent's approach
3. **Strategy Assessment**: Was it appropriate? Goal achieved?
4. **Current Status**: Ready for deposition?
5. **Next Steps**: Recommendations

### Usage

```python
# Generate summary without LLM
result = session.generate_agent_session_summary(include_llm_assessment=False)
print(result["markdown"])

# Get concise summary for LLM assessment
llm_input = session.get_summary_for_llm_assessment()
```

### Files

- `agent/session.py`: `generate_agent_session_summary()`, `get_summary_for_llm_assessment()`
- `analysis/agent_session_analyzer.py`: LLM assessment integration
- `phenix_ai/run_ai_analysis.py`: `run_agent_session_analysis()`
- `knowledge/prompts_hybrid.py`: Assessment prompt template

---

## Advice Preprocessing

The agent preprocesses user input to create structured, actionable guidance, with
built-in protection against prompt injection attacks and automatic change detection
for mid-session updates.

### Components

| Component | File | Function |
|-----------|------|----------|
| README discovery | `agent/advice_preprocessor.py` | `find_readme_file()` |
| README reading | `agent/advice_preprocessor.py` | `read_readme_file()` |
| Advice gathering | `agent/advice_preprocessor.py` | `gather_raw_advice()` |
| Input sanitization | `agent/advice_preprocessor.py` | `sanitize_advice()` |
| LLM preprocessing | `agent/advice_preprocessor.py` | `preprocess_advice()` |
| Change detection | `programs/ai_agent.py` | `_preprocess_user_advice()` |
| Server routing | `phenix_ai/run_ai_analysis.py` | `run_advice_preprocessing()` |

### Sources

1. **Direct advice**: `project_advice` PHIL parameter
2. **README files**: Found in `input_directory` if specified

### Process

1. **Gather**: Collect advice from all sources
2. **Hash**: Compute MD5 hash of raw advice for change detection
3. **Compare**: Check if advice changed from previous run
4. **Sanitize**: Remove potential prompt injection patterns
5. **Combine**: Merge with clear labeling
6. **Preprocess**: Use LLM to structure (optional)
7. **Store**: Save raw, processed, and hash in session

### Advice Change Detection

The agent detects when advice has changed between runs and automatically reprocesses:

- **Same hash** → reuse cached processed advice
- **Different hash** → reprocess advice AND re-extract directives
- **No new advice** → reuse cached version

This enables mid-session advice updates - users can provide new instructions
(like "stop" or "focus on the ligand") and have them take effect immediately.

#### Extending a Completed Workflow (Q1)

A special case arises when the workflow has already reached the `complete` phase
and the user resumes with new advice. Two independent mechanisms both need to be
active to allow follow-up work:

**Wall 1 — AUTO-STOP in PLAN (pre-existing):**
When `metrics_trend.should_stop` is True, the PLAN node normally terminates
immediately. The `advice_changed` flag suppresses this for exactly one cycle,
letting the LLM plan before reverting to normal termination.

**Wall 2 — `valid_programs = ['STOP']` in PERCEIVE (Q1 fix):**
The `complete` phase returns `['STOP']` from `get_valid_programs`. Even with
Wall 1 down, the LLM is presented with a menu of only STOP and cannot select
any follow-up program. Fix: when `advice_changed=True` AND `phase=='complete'`,
PERCEIVE steps back to the `validate` phase, which contains:

```
['phenix.polder', 'phenix.molprobity', 'phenix.model_vs_data',
 'phenix.map_correlations', 'STOP']
```

After one successful cycle, `advice_changed` is cleared in the post-execution
check and normal termination logic resumes on the next cycle.

**Complete flow on resume with new advice:**

```
_preprocess_user_advice()
  ↓ new hash ≠ stored hash → session.data["advice_changed"] = True

PERCEIVE (graph_nodes.py)
  ↓ phase == 'complete' AND advice_changed
  ↓ valid_programs ← get_valid_programs(exp, {"phase":"validate"}, ctx)
  ↓ phase_info ← {"phase": "validate", "reason": "advice_changed: stepped back"}

PLAN (graph_nodes.py)
  ↓ metrics_trend.should_stop AND advice_changed → skip AUTO-STOP this cycle
  ↓ LLM: sees new advice + validate-phase program menu → chooses phenix.polder

BUILD / VALIDATE / OUTPUT → polder command executed

Post-execution (iterate_agent.py)
  ↓ session.data["advice_changed"] = False
  ↓ next cycle: normal AUTO-STOP and complete-phase logic resume
```

**Usage:**
```bash
# Completed workflow: xtriage → phaser → refine ×3 → ligandfit → molprobity
phenix.ai_agent \
    log_directory=AIAgent_run1 \
    restart_mode=resume \
    project_advice="also run polder on chain B residue 100"
```

**Key design note:** `phenix.polder` intentionally has no `run_once` strategy
in `programs.yaml`. Different residues and ligands may each require a separate
omit map, so `polder_done=True` does not gate polder out of `valid_programs`.

### Session Management Keywords

Two Phil parameters allow inspecting and modifying an existing session without
running new crystallographic cycles.

#### `display_and_stop`

```bash
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=basic
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=detailed
```

Prints the session history (`basic`: one line per cycle; `detailed`: full
reasoning + command), then exits. Populates `self.result` via
`_finalize_session(skip_summary=True)` so `get_results()` / `get_results_as_JSON()`
work identically to a normal run.

#### `remove_last_n`

```bash
phenix.ai_agent log_directory=AIAgent_run1 remove_last_n=2
```

Removes the last N cycles from the session JSON, clears the stale AI summary
(since the cycle set has changed), rebuilds `active_files.json` and `best_files`
from the remaining history, then saves and exits. Useful for pruning a failed
refinement run before re-running.

#### Auto-set `restart_mode=resume`

Both parameters operate on an existing session directory and therefore always
require resume semantics. `run()` automatically forces `restart_mode='resume'`
before `set_defaults()` when either parameter is set, so the user does not need
to remember the flag.

#### `get_results()` safety

`run()` sets `self.result = None` as its very first statement. `get_results()`
uses `getattr(self, 'result', None)` as a defensive fallback. This prevents
`AttributeError` on any early-exit path (session management, red-flag abort,
bad parameters) that would otherwise bypass the normal result-assignment code.

### Input Sanitization

The `sanitize_advice()` function removes potentially malicious patterns:

- Instruction override attempts ("ignore all previous instructions")
- System prompt manipulation (`<system>`, `[system]`)
- Role manipulation ("you are now a...", "act as a...")
- Hidden text (null bytes, control characters)
- Excessive repetition (potential buffer overflow)

```python
from agent.advice_preprocessor import sanitize_advice, is_suspicious

# Check if input contains suspicious patterns
if is_suspicious(user_input):
    clean_input = sanitize_advice(user_input)
```

### README Discovery

Searches for these files (case-insensitive):
- `README`, `README.txt`, `README.dat`, `README.md`
- `notes.txt`, `NOTES.txt`

Files are truncated at `max_readme_chars` (default: 5000).

### LLM Preprocessing

When `preprocess_advice=True` (default), the LLM extracts:

1. **Input Files Found**: Data files mentioned in the text
2. **Experiment Type**: SAD, MAD, MR, cryo-EM, etc.
3. **Primary Goal**: What the user wants to accomplish
4. **Key Parameters**: Wavelength, resolution, sites, heavy atom type
5. **Special Instructions**: Ligands, quality targets, etc.

### Raw Advice Forwarding to Directive Extractor (Step 1)

The preprocessor LLM occasionally mangles short imperative input.  A
reproducible failure: when the user writes `"density modify and stop"`,
the openai-provider preprocessed output expands the request into a
multi-step Primary Goal and writes `Stop Condition: None`.  The user's
stop intent is silently dropped.

Step 1 mitigates this by passing **both** inputs to the LLM directive
extractor — the raw user text AND the preprocessed advice — with an
explicit authority rule:

- **Raw is the source of truth for intent**: `after_program`,
  `start_with_program`, `stop_after_requested`.
- **Processed fills gaps** the raw is silent on: files mentioned only
  in the README, experiment-type hints, parameter ranges, etc.
- **On disagreement, raw wins.**

The extractor selects between two prompts at call time:

| Caller passes | Selector | Prompt used |
|---|---|---|
| `raw_advice=None` (or empty, or equal to user_advice) | falls to else | `DIRECTIVE_EXTRACTION_PROMPT` (single-input — identical to pre-Step-1) |
| `raw_advice` non-empty AND differs from user_advice | dual-input branch | `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` |

The dual-input prompt is derived at module import time by applying
`str.replace()` to the single-input prompt, substituting a
`_SINGLE_INPUT_BLOCK` named fragment for a `_RAW_INPUT_BLOCK` named
fragment that includes the authority paragraph and both
`{raw_advice}` / `{processed_advice}` placeholders.  An import-time
`assert` guards that `_SINGLE_INPUT_BLOCK` matches exactly once, so a
future prompt edit that breaks the derivation fails fast at import
rather than silently producing the wrong prompt.

**Plumbing.** Raw advice flows end-to-end through three paths:

```
Path A (server, normal case):
  ai_agent._extract_directives
    → directive_params.ai_analysis.user_advice_raw = session.data.get("raw_advice") or None
    → ai_analysis._build_server_args encodes user_advice_raw
    → server PHIL parser populates params.ai_analysis.user_advice_raw
    → ai_analysis._run_directive_extraction_locally decodes
    → run_ai_analysis.run_directive_extraction(raw_user_advice=...)
    → directive_extractor.extract_directives(raw_advice=...)
    → dual-input prompt if raw differs from processed

Path B (client-local with run_on_server=False):
  Same as Path A but skips the transport encoding step.

Path C (session-direct, used by tests and standalone scripts):
  caller → AgentSession.extract_directives()
    → reads self.data["raw_advice"] when caller didn't pass it
    → directive_extractor.extract_directives(raw_advice=...)
```

**Backward compatibility.** Old clients send `user_advice_for_directives`
only.  New server's `user_advice_raw` PHIL field defaults to `None`, the
dispatcher passes `None` through, the selector falls back to the
single-input prompt — identical to pre-Step-1 behavior.  Server-first
deployment: clients always lag, no "new client + old server" scenario.

### File Detection Metric (Step 1F)

`run_advice_preprocessing()` emits one `PP_FILE_METRIC` line per call
when the LLM produced a distinct processed output.  The line compares
the LLM-derived `extracted_files` against a regex baseline run directly
on the raw advice, providing observational data for a future decision
about whether the LLM preprocessor is still needed for file detection
or whether a regex pass over raw advice is sufficient.

**Format** (greppable across runs):

```
PP_FILE_METRIC llm=N regex=N overlap=N llm_only=N regex_only=N \
    llm_set=[...] regex_set=[...]
```

For example:

```
PP_FILE_METRIC llm=2 regex=2 overlap=2 llm_only=0 regex_only=0 \
    llm_set=[7qz0.fa,7qz0.mtz] regex_set=[7qz0.fa,7qz0.mtz]
```

`llm_only` is the set of basenames the LLM extracted but the regex
missed; `regex_only` is the reverse.  Both sets are sorted so
cross-run `diff`/`grep -c` produces stable output.

**Skip marker.** When the LLM did not produce a distinct output
(failed, sanitized away, or no API key), the metric is suppressed and
a clearly-marked skip line is emitted instead:

```
PP_FILE_METRIC: skipped (LLM did not produce a distinct processed output — fallback to raw)
```

Without this marker, an analysis that scanned for `llm=0` lines would
conflate "LLM ran and found zero files" with "LLM did not run".

**Implementation files:**

- `agent/advice_preprocessor.py::extract_files_from_raw_advice` —
  regex pathway.  Conservative — only matches tokens with recognized
  crystallographic extensions, requires non-filename trailing
  character (whitespace, comma, quote, etc.) to avoid `file.mtzfoo`-style
  glued tails.  Returns lowercase basenames, deduplicated in
  first-occurrence order.
- `agent/advice_preprocessor.py::_summarize_file_detection_metric` —
  produces the `PP_FILE_METRIC` line.  Underscore-prefixed to signal
  the output format may change.
- `phenix_ai/run_ai_analysis.py::run_advice_preprocessing` — emits
  the metric (or skip marker) into `debug_log` after the existing
  `extract_files_from_processed_advice` call.

### PHIL Parameters

```phil
input_directory = None
  .type = path
  .help = Directory containing input files and optional README

preprocess_advice = True
  .type = bool
  .help = Use LLM to preprocess and clarify user advice

readme_file_patterns = README README.txt README.dat README.md notes.txt
  .type = strings
  .help = Filenames to look for when extracting advice

max_readme_chars = 5000
  .type = int
  .help = Truncate README files longer than this
```

### Usage

```bash
# With direct advice
phenix.ai_agent data.mtz seq.fa project_advice="Solve by MR, R-free < 0.25"

# With README in directory
phenix.ai_agent input_directory=/data/project/

# Both combined
phenix.ai_agent input_directory=/data/project/ project_advice="Prioritize geometry"

# Update advice mid-session (on restart)
phenix.ai_agent session_file=session.json project_advice="stop"
```

### Testing

```bash
python tests/tst_advice_preprocessing.py
```

Tests cover sanitization, README discovery, advice combination, file extraction,
and change detection logic.

---

## Strategic Planner (v114)

The goal-directed layer transforms the reactive single-step agent into a
strategic planner. It sits above the existing reactive loop and communicates
through the directives system. The reactive agent is unchanged.

Activated by `thinking_level=expert`. The `expert` value maps to
`advanced` inside the LangGraph pipeline (the THINK node sees
`advanced` and runs the full validation/KB/Structure Model pipeline).
The planning layer — plan generation, gate evaluation, hypothesis
testing, reports — is gated in `ai_agent.py` around the cycle loop.

### Execution Model

```
Session start:
  1. generate_plan() → StructurePlan (from templates)
  2. plan_to_directives() → directives for reactive agent
  3. Log plan for user visibility

Each cycle:
  1. Reactive agent runs (PERCEIVE → THINK → PLAN → BUILD → ...)
  2. THINK updates Structure Model from validation results
  3. After cycle: GateEvaluator.evaluate()
     → advance | continue | retreat | fallback | skip | stop
  4. On advance/retreat: generate_stage_summary()
  5. On retreat: blacklist_strategy(), compute_hash() →
     advice_changed=True
  6. plan_to_directives() for next cycle

Session end:
  generate_final_report() or generate_stopped_report()
```

### Structure Model (`agent/structure_model.py`)

Maintains ground-truth structural knowledge derived from validation
results. Never from LLM reasoning.

**Key design decisions:**

- Updated from three sources: `update_from_validation()` (every cycle),
  `update_from_xtriage()` (data assessment), `update_from_phaser()`
  (molecular replacement)
- `get_summary(detail_level=)` produces text at three levels: "brief"
  (2-3 lines for cycle display), "normal" (~500 chars for THINK prompt),
  "detailed" (full report)
- `get_current_problems()` returns severity-ordered problem list for
  hypothesis generation
- Strategy Blacklist records tried-and-failed strategies with
  `metrics_at_retreat` to prevent oscillation
- `Hypothesis` data class tracks lifecycle: proposed → testing →
  pending → confirmed/refuted/abandoned. Single active hypothesis
  budget enforced via `get_active_hypothesis()`
- `revalidation_reason` field enables re-examination of confirmed
  hypotheses when evidence weakens

**Numeric coercion (v115.07):** `from_dict()` calls `_coerce_numerics()`
after `_deep_merge()` to convert string values from JSON deserialization
to proper float/int types. Without this, `r_free - r_work` crashes with
`TypeError: unsupported operand type(s) for -: 'str' and 'float'` when
model_vs_data stores metrics as strings. All arithmetic and formatting
sites also use `_safe_float()` as belt-and-suspenders defense.

### Validation History (`agent/validation_history.py`)

Per-cycle validation snapshots persisted to session JSON.

- `record()` stores cycle number, program name, validation result,
  and extracted log metrics
- `get_metric_series(metric_name)` extracts time series for trend
  analysis (used by gate evaluator's monotonic progress gate)
- `get_phase_start_metrics(phase_start_cycle)` retrieves the snapshot
  from when a phase began (used to decide if retreat is warranted)
- `get_at_cycle()` retrieves any specific snapshot
- Serializes as list of dicts in `session.data["validation_history"]`

### Metrics Analyzer (`agent/metrics_analyzer.py`)

Extracts metrics from history and analyzes trends for stop decisions.
Called by PERCEIVE on every graph invocation.

- `derive_metrics_from_history(history)` — reconstructs `metrics_history`
  from the client-side history list. Extracts R-free, R-work, TFZ, LLG,
  resolution, and map-CC from analysis dicts and result text.
- `analyze_metrics_trend()` — detects plateau, success, and excessive
  refinement conditions. Routes to `_analyze_xray_trend()` (R-free) or
  `_analyze_cryoem_trend()` (map-model CC) based on experiment type.

**Numeric coercion (v115.07; consolidated v120):** All numeric values
extracted from history are coerced via `_safe_float()` at read time. JSON
round-tripping between client and server can turn floats into strings (e.g.
`0.385` → `"0.385"`). Without coercion, `previous - latest_r_free` crashes with
TypeError. This was initially diagnosed as the Bug 4 root cause, but see
MetricEvaluator below for the true production crash site.  As of v120
`_safe_float()` is a single shared definition in `utils/run_utils.py` (see
"`_safe_float` consolidation" below); this module imports it rather than
defining its own copy.

### Metric Evaluator (`agent/metric_evaluator.py`)

YAML-driven replacement for the hardcoded trend analysis in metrics_analyzer.
Active when `USE_YAML_METRICS=True` (the default since v115).

- `analyze_trend()` routes to `_analyze_xray_trend()` or
  `_analyze_cryoem_trend()` based on experiment type
- `is_significant_improvement()` / `calculate_improvement_rate()` — compare
  metric values with YAML-defined thresholds
- `is_plateau()` — detect stalled improvement using sliding window
- `get_target()` — resolution-dependent target lookup from metrics.yaml

**Numeric coercion (v115.07+):** This was the TRUE production crash site
for Bug 4. Since `USE_YAML_METRICS=True`, `analyze_metrics_trend()` in
`metrics_analyzer.py` routes to `analyze_refinement_trend()` which calls
`MetricEvaluator.analyze_trend()`. The evaluator re-reads raw values from
`metrics_history` without coercion — `_safe_float()` was added at all five
arithmetic entry points *within this module*: r_free extraction, CC extraction,
`is_significant_improvement`, `calculate_improvement_rate`, `is_plateau`.
(The `_safe_float()` helper itself is shared as of v120 — see
"`_safe_float` consolidation" below. Not to be confused with the five copies
of the *function* that v120 collapsed into one.)

### Plan Schema (`knowledge/plan_schema.py`)

Two data classes define the plan structure:

```
StageDef:
  id, programs, max_cycles, success_criteria,
  gate_conditions, fallbacks, skip_if, directives,
  cycles_used, status (pending/active/done/skipped)

StructurePlan:
  goal, stages, current_stage_index, strategy_hash,
  created_at_cycle, revised_at_cycle, revision_reason
```

Key operations:
- `advance()` / `retreat_to(stage_id)` / `skip_stage(stage_id)`
- `to_directives()` → reactive agent directive dict
- `compute_hash()` → strategy fingerprint (change triggers
  `advice_changed` in the reactive agent)
- `record_stage_cycle()` — counts cycles per stage. When a program
  matches a LATER stage (agent ran ahead), advances through intermediate
  stages marking them complete (v115.07). Overshoot guard: if the target
  stage is SKIPPED, verifies `new_curr` matches the program before counting.

### Plan Templates (`knowledge/plan_templates.yaml`)

Pre-defined plan skeletons:

| Template | Applicable When |
|----------|----------------|
| `mr_refine` | X-ray, has search model |
| `mr_refine_ligand` | X-ray, has search model + ligand |
| `mr_refine_lowres` | X-ray, has search model, resolution > 3.0Å |
| `mr_refine_highres` | X-ray, has search model, resolution < 1.5Å |
| `mr_refine_twinned` | X-ray, has search model, twinned data |
| `refine_placed` | X-ray, model already placed (skip MR) |
| `refine_placed_polder` | X-ray, placed model + `wants_polder` |
| `refine_placed_ligand` | X-ray, placed model + ligand |
| `refine_rebuild_placed` | X-ray, placed model + `wants_rebuild` (advice-driven) |
| `predict_refine` | X-ray, sequence only (no model, no anomalous) |
| `predict_refine_ligand` | X-ray, sequence + ligand (no model, no anomalous) |
| `mr_sad` | X-ray, has search model + anomalous atoms |
| `sad_phasing` | X-ray, anomalous atoms (no search model) |
| `sad_phasing_ligand` | X-ray, anomalous atoms + ligand |
| `validate_existing` | X-ray, `wants_validation_only` + placed model (v115.09) |
| `data_analysis_only` | X-ray, no model, no sequence |
| `cryoem_refine` | Cryo-EM |
| `cryoem_refine_ligand` | Cryo-EM + ligand |
| `cryoem_analysis_only` | Cryo-EM, no model, no sequence |

Templates encode expert crystallographic knowledge. Selection is
deterministic (rule-based). The LLM only customizes parameters
within template bounds (resolution-appropriate thresholds,
ligand-specific settings).

`mr_refine_lowres` relaxes R-free targets and disables ordered
solvent — prevents the common failure mode where the agent adds
waters at low resolution and overfits.

**Advice-driven rebuild** (`refine_rebuild_placed`, v120): plans are
template-driven — the LLM cannot add stages the template lacks — so a placed
model normally selects `refine_placed`, which has NO rebuild stage.  When the
user's advice asks to rebuild (keywords `rebuild`, `re-build`, `autobuild`,
`model building`, etc. — deliberately not bare "build" or "whatever is
necessary"), `wants_rebuild` is set and `refine_rebuild_placed` is selected
instead: `refine_placed` plus a `model_rebuilding` stage
(`[phenix.autobuild, phenix.refine]`).  Its 4-condition match
(`{xray, has_search_model, model_is_placed, wants_rebuild}`) beats
`refine_placed`'s 3.  Two further fixes make the stage actually run rather than
being skipped or abandoned: a placement-skip exception preserves the stage when
the template was explicitly chosen (see ai_agent placement-skip), a command
builder waiver lets `phenix.autobuild` run without a sequence file when a model
is available (autobuild derives the sequence and rebuilds in place), and the
reactive-deviation hold (§38.5) stops a "validate first" detour from collapsing
the plan past the un-run rebuild.

### Plan Generator (`agent/plan_generator.py`)

```python
generate_plan(data_characteristics, user_advice,
              available_models, structure_model) → StructurePlan

plan_to_directives(plan) → Dict  # for reactive agent

check_plan_revision(plan, session_data) → bool
  # Compares strategy hash, sets advice_changed on change

repair_plan(plan, user_directives) → (plan, warnings)
  # Repairs when user directives break plan prerequisites
  # Logs [Plan Repair] or [Plan Conflict]
```

### Gate Evaluator (`agent/gate_evaluator.py`)

Purely deterministic phase evaluation. No LLM.

```python
GateEvaluator.evaluate(plan, structure_model,
                       validation_history,
                       cycle_number) → GateResult

GateResult:
  action: continue | advance | retreat | fallback | skip | stop
  reason: str
  new_stage_id: str (for advance/retreat)
  blacklist_entry: dict (for retreat)
```

**Success hysteresis**: 1.5% buffer on thresholds. Plan says
`r_free: "<0.35"`, gate uses 0.345 for advancement. Prevents
oscillation on noisy single-cycle results.

**Early rebuild gate** (v115.05): `mr_refine` template includes
`r_free > 0.50 after 1 cycles → try_rebuilding`. When triggered,
the gate evaluator returns `action="advance"` to the
`model_rebuilding` stage instead of the previous weak `"fallback"`
hint. This ensures autobuild runs before the agent gives up on a
high-R-free structure (e.g., incomplete AlphaFold model).

**Retreat logic** (5 safeguards):

```
_evaluate_retreat(plan, target_id, structure_model, ...):
  1. Check Strategy Blacklist → refuse if target blacklisted
  2. Check retreat counter → max 2 per phase
  3. Check monotonic progress gate → only if worse than start
  4. Check retreat cooldown → 2+ cycles since last retreat
  5. Check retreat depth → max 1 phase back (unless explicit)

  On retreat:
    → blacklist current strategy
    → rewind plan to target phase
    → recompute strategy hash (triggers advice_changed)
    → log retreat with explanation
```

### Hypothesis Evaluator (`agent/hypothesis_evaluator.py`)

Integrated into the THINK node (advanced mode) and gate evaluator.

```python
evaluate_hypotheses(structure_model, validation_history,
                    cycle_number) → list[HypothesisResult]
  # Manages countdown, evaluation, confirmation/refutation

revalidate_confirmed(structure_model) → list[HypothesisResult]
  # Re-checks confirmed hypotheses for evidence decay

build_hypothesis_prompt(structure_model) → str
  # Generates prompt section for THINK node
  # Enforces single active hypothesis budget
```

**Verification latency**: `test_cycles_remaining` countdown
prevents premature refutation. Example: place Zn²⁺ (cycle N) →
refine (N+1, countdown) → evaluate (N+2, check coordination +
anomalous + B-factor).

**Re-validation triggers**: B-factor > 80 Å², RSCC < 0.5,
R-free spike > 0.02, difference density reappears at site.

### Explanation Engine (`knowledge/explanation_prompts.py`)

Four generators producing crystallographer-level commentary:

| Function | When | Method |
|----------|------|--------|
| `generate_cycle_commentary()` | Every cycle | Template slots from Structure Model |
| `generate_stage_summary()` | Stage transitions | LLM synthesis of multiple cycles |
| `generate_final_report()` | Completion | LLM with Structure Model constraints |
| `generate_stopped_report()` | Early stop | LLM with blacklist + recommendations |

The template-based cycle commentary needs no LLM call. Phase summaries
and reports use the LLM for narrative but constrain it with Structure
Model data (ground truth) to prevent hallucination.

### Integration Points

| Component | Integrates With | How |
|-----------|----------------|-----|
| Structure Model | `thinking_agent.py` | Updated from validation/xtriage/phaser results in THINK node |
| Structure Model | `graph_state.py` | `structure_model` field in AgentState TypedDict |
| Plan Generator | `ai_agent.py` | Called at session start; plan restored on resume |
| Plan directives | `ai_agent.py` | `plan_to_directives()` each cycle before reactive agent |
| Gate Evaluator | `ai_agent.py` | Called after each cycle; handles advance/retreat/stop |
| Strategy Hash | `ai_agent.py` | `check_plan_revision()` sets `advice_changed` |
| Hypothesis | `thinking_agent.py` | `build_hypothesis_prompt()` in THINK prompt |
| Hypothesis | `gate_evaluator.py` | `evaluate_hypotheses()` after test cycles |
| Explanation | `ai_agent.py` | `generate_stage_summary()` at transitions |
| Explanation | `event_formatter.py` | `structure_model_summary` in expert assessment |

---

## Failure Handling & Infrastructure Fixes (v115)

Version 115 addresses cycle waste from repeated failures, invalid
parameters, and file routing bugs. The fixes operate in three layers:
mechanical guards (no LLM needed), structured self-correction (LLM
with error context), and intent classification (route user advice to
the correct execution mode).

### Error Classification

`agent/error_classifier.py` classifies program failures into four
categories that determine the response:

| Category | Example | Response |
|----------|---------|----------|
| TERMINAL | Python traceback, data/model error | Immediate pivot to different program |
| PHIL_ERROR | Unrecognized PHIL parameter | Strip bad param, retry once, then pivot |
| LABEL_ERROR | MTZ column ambiguity | Select correct labels, retry |
| RETRYABLE | Transient timeout | Retry with same params |

The `should_pivot()` function in `graph_nodes.py` uses the classification
plus failure count to decide when to switch programs.

### PHIL Validation

`agent/phil_validator.py` validates LLM-generated strategy parameters
against `strategy_flags` from `programs.yaml` before command building.
Unrecognized parameters are stripped with logging.

Validation order (first match wins):
1. **Blocked params** (`_BLOCKED_PARAMS`): always stripped, even if otherwise
   allowed. E.g., `mask_atoms` for resolve_cryo_em causes RuntimeError.
2. **Exact match** against `strategy_flags` whitelist.
3. **Prefix match** against `allowed_phil_prefixes` (case-insensitive substring).
   E.g., `ncs.type` passes because `"ncs."` is a prefix for phenix.refine.
   This covers entire PHIL namespaces without listing individual params.
4. **Build-pipeline keys** (ligand, output_prefix, etc.).
5. Everything else → stripped.

The `allowed_phil_prefixes` mechanism (v115.07) allows advanced restraint
parameters (NCS, secondary structure, reference model, Ramachandran) to
pass through without individually whitelisting each sub-parameter.

**Path resolution** (v115.07): `program_registry.py :: build_command()`
detects strategy values ending in file extensions (`.pdb`, `.params`, etc.)
and resolves them to absolute paths via basename matching against known
files and the working directory. This handles LLM-generated relative
paths like `reference_model.file=4pf4.pdb`.

### Sanity Checker

`agent/sanity_checker.py` runs pre-execution checks:
- Experiment type stability (didn't change mid-workflow)
- Model exists before refinement (with search_model detection)
- Data exists for experiment type
- Repeated failures detection (4+ identical failures → stop)
- AutoBuild PHIB guard (warns after 2+ phase-related failures)

### Intent Classification

`agent/intent_classifier.py` classifies user advice into four intents
that control scope and stopping behavior:

| Intent | Example | Scope | Stop behavior |
|--------|---------|-------|--------------|
| solve | "solve the structure" | Full workflow | Convergence |
| solve_constrained | "solve with SAD" | Full workflow, method locked | Convergence |
| task | "run xtriage" | Single program | After program |
| tutorial | README with steps | Follow instructions | After described steps |

Intent is integrated into `directive_extractor.py` for both LLM and
rules paths. When intent is `task`, plan generation in `ai_agent.py`
is skipped (no multi-stage plan for a single xtriage run).

### Multi-Array MTZ Handling

When an MTZ file contains both merged (Iobs) and anomalous (I(+)/I(-))
observation arrays, programs crash without explicit labels.

The fix operates in two stages:
1. `workflow_state.py` `_detect_mtz_arrays()` scans MTZ files at
   categorization time using `iotbx.reflection_file_reader`
2. `graph_nodes.py` BUILD phase injects `obs_labels` into the strategy
   with a ranking rule: anomalous labels for SAD/MAD workflows,
   merged labels otherwise

### AutoSol Sites Estimation

`graph_nodes.py` `_estimate_anomalous_sites()` counts anomalous
scatterer sites from sequence files: Met residues for Se-SAD,
Cys+Met for S-SAD. Injected into autosol's `sites=` parameter when
not explicitly provided.

### Dual-Run Evaluation

The evaluation framework compares two run types per tutorial:
- **Solve mode**: minimal hints ("Solve the structure by standard
  procedures"), tests scientific reasoning
- **README mode**: full tutorial instructions, tests instruction-following

Components: `solve_readmes/` (21 files), `tutorial_expectations.yaml`
(per-tutorial answer key), `analyze_tutorial_runs.py` (dual-run
reporter with per-type output).

### Unsupported Program Detection

`ai_agent.py` scans READMEs for `phenix.xxx` tokens and raises `Sorry`
for programs not in the registry. An expanded ignore set covers utility
programs (`map_comparison`, `superpose_models`, `map_box`, `fmodel`,
etc.) that appear in READMEs as optional mentions but are not required.
Core unsupported programs like `ensemble_refinement` correctly block.
The same logic applies in `wxGUI2/Programs/AIAgent.py` for the GUI.

---

## Event System

The agent uses a structured event system for transparent decision logging.

### Event Flow

```
┌────────────────────────────────────────────────────────────────┐
│                        Graph Nodes                              │
│  perceive() → think() → plan() → build() → validate()          │
│       │          │         │          │                        │
│       │          │         │     ┌────┴─────┐                  │
│       │          │         │     │fallback()│                  │
│       │          │         │     └────┬─────┘                  │
│       │          │         │          │                        │
│       │          │         │     output_node()                  │
│       ▼          ▼         ▼          ▼                        │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │               state["events"] (list of dicts)            │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────┐
│              Response Building (run_ai_agent.py)               │
│  response["events"] = state["events"]                          │
│  response["events_as_simple_string"] = json.dumps(events)      │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────┐
│                    Display (ai_agent.py)                        │
│  formatter = EventFormatter(verbosity=params.verbosity)        │
│  output = formatter.format_cycle(events, cycle_number)         │
│  print(output, file=self.logger)                               │
└────────────────────────────────────────────────────────────────┘
```

### Event Types

| Type | Level | Description |
|------|-------|-------------|
| `cycle_start` | quiet | Cycle beginning |
| `cycle_complete` | quiet | Cycle finished |
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `sanity_check` | normal | Red flag or warning detected |
| `program_selected` | normal | Decision with reasoning |
| `program_modified` | normal | Program changed by rules/validation |
| `stop_decision` | normal | Whether to continue |
| `directive_applied` | normal | User directive enforced |
| `user_request_invalid` | quiet | User request unavailable |
| `expert_assessment` | normal | Thinking agent analysis (v113) |
| `files_selected` | verbose | File selection details |
| `file_scored` | verbose | Individual file scoring detail |
| `command_built` | normal | Final command |
| `thought` | verbose | LLM chain-of-thought/reasoning |
| `error` | quiet | Error occurred |
| `warning` | quiet | Non-fatal warning |
| `debug` | verbose | Internal debug information |

### Verbosity Levels

```phil
verbosity = normal
  .type = choice(quiet, normal, verbose)
```

- **quiet**: Errors, warnings, cycle summaries only
- **normal**: Key decisions and metrics (default)
- **verbose**: Full details including file selection, LLM traces, debug info

Note: `debug` is accepted as an alias for `verbose` (3 levels total).

### Implementation Files

- `agent/event_log.py` - EventType, Verbosity, EventLog class
- `agent/event_formatter.py` - EventFormatter class
- `agent/graph_nodes.py` - Event emission (`_emit()` helper)
- `agent/thinking_agent.py` - Expert assessment event emission (v113)

---

## Systematic Testing Framework (v115.08)

A 10-phase bottom-up testing framework that exercises the system boundaries
where bugs hide — file categorization, routing decisions, command building,
error classification, and LLM output resilience. Unlike the existing unit
tests (which test individual functions in isolation), these phases test
cross-module data flow through the real production pipeline.

### Design Rationale

Four code reviews of v115.05–v115.07 found 5 critical bugs that the existing
2,186 unit tests missed. All five were at module boundaries: data flowed
correctly within each module but was corrupted, lost, or misinterpreted at
the handoff. The systematic framework was designed to test these boundaries
directly.

### Phase Architecture

Phases are ordered from static analysis (fast, no I/O) to full pipeline
simulation (slow, creates temp files). Each phase produces a machine-readable
YAML findings file in `findings/`.

| Suite | Phase | What It Tests | Gate Behavior |
|-------|-------|---------------|---------------|
| S0 | Static Audit | Parse check, bare except, import fallbacks | FAIL on blocking issues |
| S1 | Contract Gaps | AST coverage map of 128 functions in 4 modules | FAIL never (informational PARTIAL) |
| S2 | Path Consistency | YAML vs hardcoded categorization for 10 tutorials | FAIL on new unexpected divergence |
| S3 | Session Round-Trip | JSON symmetry + AgentSession save/load/pipeline | FAIL on any test failure |
| S4 | History Flags | Flag writer/reader consistency across modules | FAIL on any test failure |
| S5 | Category-Consumer | input_priorities + fallback_categories alignment | FAIL on unexpected missing input |
| S6 | Routing Simulation | 3-cycle routing for 32 tutorials (real code) | FAIL on unexpected stuck tutorial |
| S7 | Command Building | CommandBuilder.build() for 15 tutorial×program combos | FAIL on any build failure |
| S8 | Error Classification | 3 classifiers × 30+ patterns, overlap detection | FAIL on any test failure |
| S9 | LLM Perturbation | Hallucinated files/programs/params, truncated JSON | FAIL on any test failure |

Phases S6–S8 are skipped in `--quick` mode (~3s saved).

### Whitelist Pattern

Phases 2, 6, and 7 use a whitelist for expected failures that cannot be
fixed in the test environment (e.g., cryo-EM maps with invalid CCP4
headers, known YAML/hardcoded divergences). The pattern:

1. Expected issues are listed in a set (e.g., `_EXPECTED_STUCK`)
2. Each issue is checked against the whitelist
3. Expected issues → PARTIAL status (does not raise)
4. Any NEW unexpected issue → FAIL status (raises AssertionError)

This ensures regressions are caught while known test-env limitations
don't produce false alarms.

### Key Modules Tested

The framework focuses on the 4 highest-risk modules:

- **`workflow_state.py`** — File categorization, history analysis, phased
  detection. 30 functions, 13 with zero coverage (Phase 1).
- **`workflow_engine.py`** — Step detection, program routing, stop decisions.
  32 functions including `_is_at_target` (highest-risk untested).
- **`command_builder.py`** — File selection, command assembly. 29 functions
  including `_find_file_for_slot` (second-highest-risk).
- **`graph_nodes.py`** — LangGraph node functions, LLM intent parsing. 37
  functions.

### Implementation Files

| File | Purpose |
|------|---------|
| `tests/tst_phase0_static_audit.py` | S0: Static analysis |
| `tests/tst_phase1_contract_gaps.py` | S1: Coverage map |
| `tests/tst_phase2_path_consistency.py` | S2: YAML vs hardcoded |
| `tests/tst_phase3_serialization_symmetry.py` | S3: Session round-trip |
| `tests/tst_phase4_history_flags.py` | S4: Flag consistency |
| `tests/tst_phase5_error_classification.py` | S8: Error classifiers |
| `tests/tst_phase6_category_consumer.py` | S5: Category alignment |
| `tests/tst_phase7_routing_simulation.py` | S6: Routing simulation |
| `tests/tst_phase8_command_building.py` | S7: Command building |
| `tests/tst_phase9_llm_perturbation.py` | S9: LLM resilience |
| `findings/*.yaml` | Machine-readable phase outputs |
| `docs/PHASE_REVIEW_REPORT.md` | Review findings (54 fixes) |

### Test Environment Conventions

All phase scripts follow the cctbx test convention: test functions raise
`AssertionError` on failure. The libtbx mock boilerplate (~25 lines) is
applied at the top of each script before importing any agent module.
`_PHASE_COLUMN_CACHE` must be cleared between `_categorize_files` calls.
`ws._mtz_has_phase_columns` is monkeypatched to `lambda f: False` (iotbx
unavailable). CCP4 map files need a valid 1024-byte header to pass
`_is_valid_file()` — use `create_ccp4_map()` from Phase 8.

---

## Future Directions

### Known gaps and active limitations

**Text-only information channel.** The agent's entire understanding
of program results comes from log text and numerical metrics. Density
maps carry spatial information — ligand shape, disorder, connectivity,
unexplained blobs — that text summaries do not capture. This is the
single largest limitation of the current system: an experienced
crystallographer looking at a difference map can immediately see what's
wrong, but the agent can only read the numbers. Adding even basic
spatial awareness (e.g., difference density peak statistics, local
correlation per residue) would significantly improve decision quality
in the refinement and ligand-fitting stages. Adopting
`results_as_json()` from newer PHENIX programs (see Potential
improvements below) would provide richer structured data than log
parsing, but would not address the spatial gap — that requires density
map analysis.

**Hypothesis testing infrastructure.** The hypothesis system (v114) is
architecturally complete — the evaluator, lifecycle management, single-
budget constraint, verification latency, and revalidation logic all
work. However, it does not fire reliably in practice. The prompt only
invites hypothesis proposals when the Structure Model has unresolved
problems or ≥ 2 positive difference density peaks above 4 sigma. In
most tutorial runs the structure either improves steadily (no problems
to trigger a hypothesis) or fails quickly (not enough validation data).
When the invitation does fire, the LLM must return a correctly
structured JSON with `hypothesis`, `test_program`, `confirm_if`, and
`refute_if` fields — which it does not always do. Making this feature
work in real sessions requires either lowering the trigger threshold
(at the risk of spurious hypotheses), improving the prompt to elicit
reliable JSON, or both.

**Single dataset per session.** The agent assumes one dataset. Multi-
crystal merging, serial crystallography data reduction, and ensemble
strategies are not supported. Supporting these would require changes
to the workflow engine (new phases), the file categorizer (dataset
grouping), and possibly a multi-session coordinator.

**Ligand PDB plan selection (fixed v115.05).** `_build_context()`
in `plan_generator.py` had two bugs that prevented correct plan
selection when a ligand is provided as a PDB file (e.g.,
AF_bromodomain_ligand tutorial with `7qz0_ligand.pdb`):

(A) Every `.pdb` file set `has_search_model=True` (line ~200),
including ligand PDB files. A ligand PDB is not a search model — this
falsely triggered MR templates instead of predict-and-build templates.

(B) The ligand name detection (`_ligand_pdb_hints`: "ligand", "lig_",
etc.) only fired when `len(pdb_files) >= 2` (line ~249). When the
ligand PDB was the only PDB file, `has_ligand_code` stayed False and
the plan had no ligandfit stage.

Fix: ligand hints are now checked *during* the file scan, before
setting `has_search_model`. PDB files matching ligand hints set
`has_ligand_code=True` instead of `has_search_model=True`. The
`len(pdb_files) >= 2` guard was removed. With both fixes, the
correct template `predict_refine_ligand` is selected for tutorials
that provide only a ligand PDB + sequence + data.

**Program coverage.** 23 PHENIX programs are registered. Notable gaps
include `ensemble_refinement`, local map sharpening, `map_box`,
`map_comparison`, and `superpose_models`. Deferred items I6
(unsupported programs) and I7 (tar.gz input handling) from the v115
plan are pending workflow engine expansion. Adding programs starts
with YAML definitions but in practice requires iterating on file
categorization guards, error recovery patterns, content-based checks,
and command postprocessor special cases — see the "Adding a New Tool"
discussion in OVERVIEW.md.

**Missing guard in cryo-EM target check.** The `_is_at_target`
clashscore path (line 1334 of `graph_nodes.py`) does not have a
`rsr_count >= 1` guard like the X-ray path does. In practice this
is a theoretical-only gap: cryo-EM programs don't produce clashscore
until after real-space refinement has run, so the unguarded path is
never reached with current tutorials. Worth a future cleanup but low
risk.

**CIF model categorization (v115.09).** When a user provides a
macromolecular model as a `.cif` file (mmCIF format), the YAML
categorizer places it in the generic `cif` category — not `model` or
`search_model`. This causes `has_model=False` in `build_context`,
breaking all model-dependent routing (validation shortcut, placement
probes, refinement). In practice most users provide PDB-format models;
CIF-format ligand restraints are correctly categorized. Fix requires
adding mmCIF model detection to the YAML category rules (checking for
`_atom_site.` loop or similar structural markers).

**Preprocessing stop override (v115.09).** `ai_agent.py` line 2761
has a `_preprocessing_programs` set (`xtriage`, `mtriage`) that
unconditionally clears the `after_program` stop condition and
overrides `intent: task` → `intent: solve`. This prevents stopping
after xtriage/mtriage even when the user explicitly says "run mtriage
and stop." The fix is to check for explicit stop language
(`_has_explicit_stop` regex) before clearing. Located in
`$PHENIX/modules/phenix/phenix/programs/ai_agent.py`, outside the
langchain directory.

### Design tensions

**Rule D ("fail closed") vs LLM error recovery.** The command sanitizer
strips bare parameters not in a program's `strategy_flags` allowlist.
This prevents hallucinated parameters from reaching PHENIX, but also
strips legitimate recovery parameters that the LLM correctly
identifies (e.g., `rebuild_in_place=False` for autobuild sequence
mismatch). The current mitigation is to expand `strategy_flags` for
programs where recovery params are known. A future "warn but keep"
mode could allow unrecognized parameters through if they pass PHIL
validation, relying on the catch-all blacklist (v112.76) as a safety
net for parameters that cause actual failures.

**Plan template rigidity vs expert reasoning.** The plan templates are
deterministic — selected at session start and locked to a phasing
strategy. When the template is wrong (e.g., SAD template locked in
despite anomalous measurability of 0.03), expert reasoning correctly
diagnoses the problem but has no mechanism to override the template.
The v115.05 anomalous gate addresses the specific case of negligible
anomalous signal (measurability < 0.05), but the
general problem — how should the planner respond when the expert says
the strategy is wrong? — remains unsolved. Options include: allowing
the THINK node to flag plan-incompatible evidence that triggers a
plan revision, or adding "escape hatch" gates that the evaluator
checks before entering high-commitment stages.

**Three independent error classification systems.** Error handling
is split across three classifiers that evolved at different times,
have overlapping patterns, and are not aware of each other.

*System 1: `ai_agent.py::_classify_error()` (original).* The oldest
classifier. ~60 hardcoded substring matches (no regex). Two output
categories: `INPUT_ERROR` (agent's fault — don't count) vs
`REAL_FAILURE` (real problem — count it). Used only for deciding
whether a failure enters the cycle history. Knows nothing about the
YAML files.

*System 2: `agent/error_classifier.py::classify_error()` (v115).*
Lives in the graph; PERCEIVE calls it at the start of the next cycle.
Hardcoded regex patterns, more sophisticated. Five categories:
`TERMINAL`, `PHIL_ERROR`, `AMBIGUOUS_PHIL`, `LABEL_ERROR`, `RETRYABLE`.
Extracts details (bad param names, suggestions). Feeds `should_pivot()`
which excludes the failed program from valid_programs. Provides error
context to THINK and PLAN prompts. Also knows nothing about the YAML.

*System 3: YAML-driven (two files, two consumers).*
`recoverable_errors.yaml` + `ErrorAnalyzer` detects errors the agent
can auto-fix (currently: ambiguous data labels, ambiguous experimental
phases). `diagnosable_errors.yaml` + `DiagnosisDetector` detects
terminal errors needing LLM diagnosis (5 error types: crystal symmetry
mismatch, model outside map, SHELX not installed, unknown PHIL
parameter, polymer special position).

The execution order per cycle is:
1. Program runs → result text
2. `ErrorAnalyzer` (YAML recoverable) → if match, set recovery flags
3. `DiagnosisDetector` (YAML diagnosable) → if match, stop run
4. `_classify_error` (hardcoded, ai_agent.py) → classify for history
5. Next cycle: PERCEIVE → `classify_error` (hardcoded,
   error_classifier.py) → feeds THINK/PLAN + pivot

Where they overlap and conflict:

- *Unknown PHIL parameter*: `diagnosable_errors.yaml` says terminal
  (stop, diagnose). `_classify_error` says INPUT_ERROR (agent's fault,
  fixable). `error_classifier.py` says PHIL_ERROR (strip params, retry).
  Three different behaviors for the same error. In practice the
  execution order saves this — DiagnosisDetector fires first and stops
  the run — but the systems disagree on the correct response.

- *Crystal symmetry mismatch*: Present in `diagnosable_errors.yaml`
  but absent from `error_classifier.py`. If the detector misses it
  (e.g., unusual error phrasing), the graph-level classifier has no
  fallback.

- *SHELX not installed*: Same gap — only in the YAML, not in the
  hardcoded classifier.

- *Ambiguous data labels*: `recoverable_errors.yaml` auto-fixes them;
  `error_classifier.py` classifies them as LABEL_ERROR (retryable with
  different approach); `_classify_error` calls them INPUT_ERROR. The
  YAML recovery fires first, so the other two rarely see it, but they
  would disagree if they did.

The practical risk is low because the execution order is stable and the
YAML files are small. But every new error pattern requires checking all
three systems to ensure they agree, and nothing enforces consistency.

*Recommended consolidation path*: Make `error_classifier.py` the single
classifier. Have it load the YAML files as additional pattern sources
rather than duplicating patterns in code. Map its five categories
(TERMINAL, PHIL_ERROR, LABEL_ERROR, AMBIGUOUS_PHIL, RETRYABLE) to the
two needed by `_classify_error` (INPUT_ERROR, REAL_FAILURE) so the
ai_agent.py method becomes a thin wrapper. The YAML files continue to
define the *data* (patterns, hints, recovery strategies) while the
classifier provides the *logic* (matching, extraction, category
assignment). This eliminates the three-way overlap and ensures that
adding a new error pattern requires editing one place.

### Potential improvements

**Structured results via `results_as_json()`.** Newer PHENIX programs
built on `ProgramTemplate` expose a `results_as_json()` method that
returns metrics, output files, and status as structured JSON — no log
parsing needed. The agent currently extracts all metrics by regex-
matching log text, either via YAML `log_parsing` patterns in
`programs.yaml` or hardcoded extractors in `log_parsers.py`. This is
fragile: log format changes silently break extraction, multi-line
patterns are hard to express in YAML, and programs that don't print
clean key=value lines require custom Python parsers.

Switching to `results_as_json()` where available would:
- Eliminate regex fragility for programs that support it
- Provide richer data (e.g., per-residue validation, per-chain
  statistics) that log text summarizes or omits entirely
- Give the THINK node structured input instead of text it must
  re-parse from the LLM's analysis
- Reduce the per-program integration cost (no `log_parsing` YAML
  section needed, no hardcoded extractors)

The migration path is incremental: programs that support
`results_as_json()` can be switched one at a time, with the existing
log-parsing path as a fallback for older programs. The PERCEIVE node
would check for a JSON results file first, then fall back to log text
extraction. The `programs.yaml` entry for each program could add a
`has_json_results: true` flag to signal which path to use.

Not all programs support this yet — legacy programs and those not
built on `ProgramTemplate` will continue to need log parsing. But
as more programs are updated, the regex-heavy extraction path can
be phased out gradually.

**Density map awareness.** Even without full spatial map interpretation,
extracting summary statistics from difference density maps —
peak heights, peak locations relative to the model, local CC per
residue — would give the THINK node evidence for ligand placement,
disorder, and model errors that is currently invisible. This could
plug into the existing expert KB and hypothesis systems.

**Learning from completed sessions.** The agent currently starts fresh
each session. Collecting outcome data from completed runs — which
program sequences solved which types of structures, which strategies
worked at which resolutions — could inform template selection and
strategy recommendations. This would require a session outcome
database and a retrieval mechanism, but the RAG pipeline already
provides the infrastructure for document-grounded retrieval.

**Interactive checkpoints.** The agent currently runs autonomously
or in stepwise mode (stop after prediction for manual inspection).
A middle ground would be structured checkpoints where the agent
presents its assessment and asks the user to confirm or redirect
before committing to an expensive step (e.g., autosol with SAD
phasing, or autobuild after marginal MR). The directive system
already supports this — `after_program` could be extended to
pause-and-ask rather than just suppress auto-stop.

**Broader program integration.** Many PHENIX tools that appear in
tutorial READMEs are not yet registered: `map_box` for extracting
map regions, `superpose_models` for comparing solutions,
`ensemble_refinement` for modeling disorder, and various map
utilities. Each requires YAML definitions plus the inevitable
edge-case iteration, but the pipeline architecture does not need
to change.

**Multi-model and multi-dataset workflows.** Supporting ensemble
strategies (multiple models from different MR solutions),
multi-crystal merging, or comparative analysis across datasets
would require a session model that tracks multiple parallel
branches rather than a single linear cycle history. This is a
significant architectural change but would enable the agent to
handle the more complex structure determination scenarios that
currently require manual intervention.
