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
│       ├── 1. sanitize_response()                                 │
│       ├── 2. json.dumps()                                        │
│       └── 3. encode_for_rest()                                   │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Decision**: LocalAgent performs the full encode/decode roundtrip locally, even though it could skip encoding. This ensures that transport bugs are caught during local testing rather than only appearing in production server scenarios.

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
- Stop decisions: hard stops (after_cycle, metrics targets) in PERCEIVE; after_program minimum-run guarantee in PLAN; consecutive-program cap (→ PERCEIVE)
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
- `sanitize_command()` — Rules A–D: strip placeholders, blacklisted params, hallucinated cross-program params, bare unscoped params. Rule B2 (v112.72): validates `space_group=` values and strips non-space-group words.
- `inject_user_params()` — append user key=value params missing from command (scope-matched for dotted keys, strategy_flags-validated for bare keys)
- `inject_crystal_symmetry()` — append unit_cell/space_group from directives (validates with `_is_valid_space_group()`)
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
│  - Add directive-required programs (after_program target)        │
│  - Add start_with_program to front of list (multi-step flows)   │
│  - Add STOP if skip_validation=true                             │
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
│  extractor can only name one program.                           │
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
- Invariants (rules that must be satisfied)
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

The agent supports three LLM providers: **google**, **openai**, and **ollama**. The
provider is set via `params.communication.provider` (or the `LLM_PROVIDER` env var)
and flows consistently through every LLM call point — there is no cross-provider
fallback.

### LLM Call Points

The system uses three distinct LLM roles, each with provider-specific model defaults:

| Activity | Function | ollama | google | openai |
|---|---|---|---|---|
| **Agent decisions** (PLAN node) | `get_planning_llm()` → `get_expensive_llm()` | qwen3:32b (json_mode) | gemini-2.5-pro | gpt-5 |
| **Log summarization** | `get_cheap_llm()` | qwen2.5:7b | gemini-2.5-flash-lite | gpt-5-nano |
| **RAG analysis** | `get_expensive_llm()` | qwen3:32b | gemini-2.5-pro | gpt-5 |
| **Directive extraction** | `call_llm_simple()` | direct ollama HTTP | gemini-2.0-flash | (provider default) |

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
| `utils/run_utils.py` | `setup_llms()` — creates all three LLM roles for analysis |
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

1. Add provider case to `get_llm_and_embeddings()` in `core/llm.py`
2. Add validation check in `validate_provider()` in `agent/graph_nodes.py`
3. Add `_call_<provider>_llm()` function in `agent/api_client.py`
4. Add rate limit handler in `agent/rate_limit_handler.py` (optional)

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

**`bad_inject_params` flow** (added in v112.66): Parameters that previously
caused PHENIX errors are blacklisted and propagated to BUILD for stripping:

```
ai_agent.py: session.data["bad_inject_params"]
    → session_info["bad_inject_params"]
    → build_session_state() → session_state
    → build_request_v2() → normalized_session_state
    → transport encode/decode
    → run_ai_agent.py → create_initial_state(bad_inject_params=...)
    → graph state["bad_inject_params"]
    → BUILD: postprocess_command(bad_inject_params=set(...))
```

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
- `validation/` — Command validation framework
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
       [if ligand] ligandfit → pdbtools → refine
```

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

**Known issue — preprocessing stop override** (`ai_agent.py` line 2761):
The `_preprocessing_programs` set (`xtriage`, `mtriage`) causes
`after_program` to be unconditionally cleared, even when the user
explicitly says "run mtriage and stop." The intent override also
changes `task` → `solve`. Fix planned for v115.10: add
`_has_explicit_stop` regex check before clearing.

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

Once an MTZ with R-free flags is identified, it is **locked for the entire session**. This ensures consistent cross-validation statistics throughout refinement. The locked MTZ has absolute priority (Priority 0) over all other file selection mechanisms.

### File Selection Priority

When building commands, files are selected in this order:

0. **Locked R-free MTZ** — For MTZ inputs in X-ray refinement (highest priority)
1. **Best Files** — From the BestFilesTracker scoring system
2. **Categorized Files** — From workflow_state file categorization
3. **Extension Matching** — Search available_files by extension

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
  whether the program needs anomalous or merged data.
- **`ambiguous_experimental_phases`**: MTZ file contains multiple phase arrays.
  The agent selects standard vs anomalous HL coefficients based on context.

### Configuration

`knowledge/recoverable_errors.yaml` defines:
- **`errors`**: Detection patterns, extraction regexes, resolution strategies, max retries
- **`label_patterns`**: How to classify data labels as anomalous or merged
- **`program_data_preferences`**: Which programs need anomalous vs merged data
- **`context_keywords`**: Workflow hints from user advice (e.g., "SAD" → anomalous)

### Adding a New Recoverable Error

1. Add an entry under `errors:` in `recoverable_errors.yaml` with detection patterns
   and extraction regexes
2. Add a resolution handler in `error_analyzer.py` (`resolve_error()` method)
3. The fallback node in the graph will automatically use the new pattern

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

**Numeric coercion (v115.07):** All numeric values extracted from history
are coerced via `_safe_float()` at read time. JSON round-tripping between
client and server can turn floats into strings (e.g. `0.385` → `"0.385"`).
Without coercion, `previous - latest_r_free` crashes with TypeError.
This was initially diagnosed as the Bug 4 root cause, but see
MetricEvaluator below for the true production crash site.

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
`metrics_history` without coercion — `_safe_float()` was added at all 5
arithmetic entry points: r_free extraction, CC extraction,
`is_significant_improvement`, `calculate_improvement_rate`, `is_plateau`.

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

Twelve pre-defined plan skeletons:

| Template | Applicable When |
|----------|----------------|
| `mr_refine` | X-ray, has search model |
| `mr_refine_ligand` | X-ray, has search model + ligand |
| `mr_refine_lowres` | X-ray, has search model, resolution > 3.0Å |
| `mr_refine_highres` | X-ray, has search model, resolution < 1.5Å |
| `mr_refine_twinned` | X-ray, has search model, twinned data |
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

See [TRANSPARENCY_LOGGING.md](../project/TRANSPARENCY_LOGGING.md) for full details.

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
