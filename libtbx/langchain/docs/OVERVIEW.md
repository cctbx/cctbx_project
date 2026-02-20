# PHENIX AI Agent V2 - Technical Overview

## Introduction

The PHENIX AI Agent is an intelligent automation system for macromolecular structure determination. It combines:

- **LLM Decision Making** - Uses Claude/Gemini to interpret context and make intelligent program choices
- **Rules-Based Validation** - Enforces workflow constraints and validates all decisions
- **YAML Configuration** - All domain knowledge externalized to editable config files
- **Structured Event Logging** - Transparent decision tracking at multiple verbosity levels

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              USER INTERFACE                                  │
│  phenix.ai_agent original_files="data.mtz seq.fa" project_advice="..."      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            programs/ai_agent.py                              │
│  • Parse parameters                                                          │
│  • Initialize session with BestFilesTracker                                 │
│  • Choose LocalAgent or RemoteAgent                                          │
│  • Display results with EventFormatter                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                    ┌─────────────────┴─────────────────┐
                    ▼                                   ▼
         ┌──────────────────┐                ┌──────────────────┐
         │   LocalAgent     │                │   RemoteAgent    │
         │ (phenix_ai/)     │                │ (phenix_ai/)     │
         │                  │                │                  │
         │ Direct function  │                │ HTTP POST to     │
         │ call to          │                │ PHENIX server    │
         │ run_ai_agent     │                │                  │
         └────────┬─────────┘                └────────┬─────────┘
                  │                                   │
                  └─────────────────┬─────────────────┘
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         phenix_ai/run_ai_agent.py                            │
│  • Build LangGraph state                                                     │
│  • Execute graph nodes (perceive → plan → build → validate → output)        │
│  • Retry loop with fallback on persistent validation failures               │
│  • Build response with events                                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            LangGraph Pipeline                                │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐                  │
│  │ perceive │ → │   plan   │ → │  build   │ → │ validate │                  │
│  └────┬─────┘   └────┬─────┘   └────┬─────┘   └────┬─────┘                  │
│       │              │   ▲          │              │                         │
│       │              │   └──────────│──────────────┘ (retry < 3)            │
│       │              │              │              │                         │
│       │              │              │         ┌────┴─────┐                   │
│       │              │              │         │ fallback │ (retry >= 3)      │
│       │              │              │         └────┬─────┘                   │
│       │              │              │              │                         │
│       │              ▼              ▼              ▼                         │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                       ┌──────────┐                                  │    │
│  │                       │  output  │ → END                            │    │
│  │                       └──────────┘                                  │    │
│  │                      state["events"] (list)                          │    │
│  │  • STATE_DETECTED, METRICS_EXTRACTED, METRICS_TREND                 │    │
│  │  • PROGRAM_SELECTED, USER_REQUEST_INVALID, DIRECTIVE_APPLIED        │    │
│  │  • FILES_SELECTED, COMMAND_BUILT, STOP_DECISION                    │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Core Components

### 1. Graph Nodes (`agent/graph_nodes.py`)

The agent uses a LangGraph pipeline with six nodes:

| Node | Purpose | Key Actions |
|------|---------|-------------|
| **perceive** | Understand current state | Categorize files, detect workflow state, extract metrics, analyze trends |
| **plan** | Decide next action | Call LLM or rules engine, validate against workflow, check user directives |
| **build** | Generate command | Select files using BestFilesTracker, build command string |
| **validate** | Check output | Sanity checks, red flag detection; retry up to 3× on failure |
| **fallback** | Last resort | If validate fails 3×, use mechanical/rules-based program selection |
| **output** | Format response | Package decision, events, and command into response |

The graph has conditional routing: perceive can skip to output (on red flag abort), validate can loop back to plan (retry) or route to fallback (max retries exceeded).

### 2. Workflow Engine (`agent/workflow_engine.py`)

Interprets `workflows.yaml` to determine:
- **Current phase**: X-ray has 9 phases (analyze, obtain_model, molecular_replacement, experimental_phasing, build_from_phases, refine, combine_ligand, validate, complete); cryo-EM has 9 phases (analyze, obtain_model, dock_model, check_map, optimize_map, ready_to_refine, refine, validate, complete)
- **Valid programs** for current phase
- **Transition conditions**
- **Stop criteria** (target reached, plateau detected)

### 3. Command Builder (`agent/command_builder.py`)

Unified command generation that:
- Reads program definitions from `programs.yaml`
- Selects appropriate input files using BestFilesTracker
- Applies strategy flags and defaults
- Tracks WHY each file was selected (for transparency)

### 4. Event System (`agent/event_log.py`, `agent/event_formatter.py`)

Structured logging for transparency:

| Event Type | Verbosity | Description |
|------------|-----------|-------------|
| `cycle_start` | quiet | New cycle beginning |
| `cycle_complete` | quiet | Cycle finished |
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution parsed |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `sanity_check` | normal | Red flag or warning detected |
| `program_selected` | normal | Decision with full reasoning |
| `program_modified` | normal | Program changed by rules/validation |
| `stop_decision` | normal | Whether to continue |
| `directive_applied` | normal | User directive enforced |
| `user_request_invalid` | quiet | User requested unavailable program |
| `files_selected` | verbose | File selection with reasons |
| `file_scored` | verbose | Individual file scoring detail |
| `command_built` | normal | Final command |
| `thought` | verbose | LLM chain-of-thought/reasoning traces |
| `error` | quiet | Error occurred |
| `warning` | quiet | Non-fatal warning |
| `debug` | verbose | Internal debug information |

### 5. File Categorization (`agent/graph_nodes.py`, `knowledge/file_categories.yaml`)

Semantic file classification using rules from `file_categories.yaml`:

| Category | Description | Examples |
|----------|-------------|----------|
| `data_mtz` | Reflection data (Fobs, R-free) | data.mtz, refine_001_data.mtz |
| `map_coeffs_mtz` | Map coefficients (calculated phases) | refine_001_001.mtz, denmod_map_coeffs.mtz |
| `sequence` | Sequence files | seq.fa, protein.fasta |
| `model` | **Positioned** models | refine_001.pdb, phaser_output.pdb |
| `search_model` | Templates/predictions | template.pdb, alphafold.pdb |
| `full_map` | Cryo-EM maps | map.mrc, emd_1234.ccp4 |
| `half_map` | Half-maps | half1.mrc, map_half_a.mrc |
| `ligand` | Ligand files | ligand.pdb, ATP.cif |

**Semantic distinction**: `model` = already positioned in crystal/map; `search_model` = template not yet placed.

### 6. Best Files Tracker (`agent/best_files_tracker.py`)

Tracks the best file of each type across cycles:
- Scores files based on metrics (R-free, resolution, cycle number)
- **Dual MTZ tracking**:
  - `data_mtz`: Locks after first R-free flags generated (consistency for refinement)
  - `map_coeffs_mtz`: Always prefers most recent (maps improve with refinement)
- Provides `best_files["model"]`, `best_files["data_mtz"]`, `best_files["map_coeffs_mtz"]` to CommandBuilder

### 7. User Directives (`agent/directive_extractor.py`)

Parses natural language guidance from `project_advice`:
- **Stop conditions**: "stop after one refinement", "stop when R-free < 0.25"
- **Workflow preferences**: "skip autobuild", "use SAD phasing"
- **Task focus**: "focus on ligand fitting"

### 8. Safety Checks (`agent/sanity_checker.py`, `agent/directive_extractor.py`)

Multiple layers of validation to prevent errors:

| Layer | Location | Examples |
|-------|----------|----------|
| **Sanity Checks** | Pre-execution | No data for workflow, model not positioned, repeated failures |
| **Directive Validation** | Post-LLM | Invalid program names, conflicting stop conditions |
| **Workflow Validation** | State machine | Invalid phase transitions, wrong experiment type |
| **Post-Processing** | After extraction | Ligand workflow conflict resolution |

Key sanity issues that trigger abort:
- `no_data_for_workflow` - Missing data_mtz (X-ray) or map (cryo-EM)
- `search_model_not_positioned` - Trying to refine before MR/docking
- `no_model_for_refine` - No model available for refinement
- `repeated_failures` - Same error 3+ times

See `docs/SAFETY_CHECKS.md` for the complete auto-generated list.

### 9. RAG Pipeline (`rag/`, `utils/query.py`, `analysis/analyzer.py`)

The agent uses Retrieval-Augmented Generation to ground LLM responses in
PHENIX documentation. The pipeline has three stages:

1. **Retrieval**: A Chroma vector store (built from PHENIX docs via
   `run_build_db.py`) returns the top 20 candidate documents using
   embedding similarity.

2. **Reranking**: FlashRank, a local cross-encoder model
   (`ms-marco-MiniLM-L-12-v2`, ~34MB), reranks the 20 candidates and
   selects the top 8 most relevant. This runs on CPU with no API key —
   the `flashrank` package must be installed where analysis runs locally.

3. **Generation**: The reranked documents are passed as context to the
   LLM (Google or OpenAI) which generates the final response.

The reranking retriever is created in `rag/retriever.py` via
`create_reranking_retriever()` and used by both log analysis
(`analysis/analyzer.py`) and documentation queries (`utils/query.py`).

---

## Configuration Files

### programs.yaml

Defines each PHENIX program the agent can run:

```yaml
phenix.refine:
  description: "Crystallographic refinement"
  category: refinement
  experiment_types: [xray]
  
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
        priority_patterns: [refine]  # Prefer refined models
      data_mtz:
        extensions: [.mtz, .sca, .hkl, .sdf]
        flag: ""
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
      - pattern: "*_data.mtz"
        type: data_mtz
      - pattern: "*_refine_*.mtz"
        type: map_coeffs_mtz
    metrics:
      - r_free
      - r_work
      - bonds_rmsd
      - angles_rmsd
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-free"
    r_work:
      pattern: 'R-work\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-work"
    bonds_rmsd:
      pattern: '[Bb]onds?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Bonds RMSD"
    angles_rmsd:
      pattern: '[Aa]ngles?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Angles RMSD"
```

### workflows.yaml

Defines workflow state machines:

```yaml
xray:
  phases:
    analyze:
      description: "Analyze data quality"
      programs:
        - phenix.xtriage
      transitions:
        on_complete: obtain_model
    
    obtain_model:
      description: "Get initial model"
      programs:
        - program: phenix.predict_and_build
          preferred: true
          conditions:
            - has: sequence
            - not_done: predict_full  # Don't re-run after full workflow completes
        - program: phenix.phaser
          conditions:
            - has: search_model
            - not_done: phaser  # MR should only run once
        - program: phenix.autosol
          conditions:
            - has: sequence
            - has: anomalous        # Requires anomalous signal detected by xtriage
            - not_done: autosol
      transitions:
        on_complete: refine
        if_predict_only: molecular_replacement  # Stepwise mode
    
    # MR-SAD: after phaser places model, use anomalous signal with autosol
    experimental_phasing:
      description: "MR-SAD phasing with placed model"
      programs:
        - program: phenix.autosol
          conditions:
            - not_done: autosol
      transitions:
        on_complete: build_from_phases
    
    refine:
      description: "Improve model"
      programs:
        - program: phenix.refine
        - program: phenix.autobuild
          conditions:
            - r_free: "> autobuild_threshold"
            - has: sequence
            - not_done: autobuild
        - program: phenix.autobuild_denmod
          conditions:
            - has: ligand_file
            - has: sequence
            - not_done: autobuild_denmod
            - refine_count: "> 0"    # Must refine first
          hint: "Run density modification before ligand fitting"
        - program: phenix.ligandfit
          conditions:
            - has: ligand_file
            - not_done: ligandfit
            - r_free: "< 0.35"
            - refine_count: "> 0"    # Must refine first
        - program: phenix.polder
          conditions:
            - has: model
            - has: data_mtz
      repeat:
        max_cycles: 4
        until:
          any:
            - r_free: "< target_r_free"
            - condition: plateau
              cycles: 2
              threshold: 0.005
      transitions:
        on_ligandfit: combine_ligand
        on_target_reached: validate
        on_plateau: validate
        on_max_cycles: validate

  targets:
    r_free:
      default: 0.25
      by_resolution:
        - range: [0, 1.5]
          value: 0.20
        - range: [1.5, 2.5]
          value: 0.25
        - range: [2.5, 3.5]
          value: 0.30
```

#### Program Execution Controls

The workflow uses two mechanisms to prevent programs from running repeatedly:

**1. `not_done` conditions** (in workflows.yaml)

Programs can specify `not_done: <flag>` to prevent re-runs:

| Flag | Program | Meaning |
|------|---------|---------|
| `predict_full` | predict_and_build (X-ray) | Full workflow (prediction+MR+building) completed |
| `predict` | predict_and_build (cryo-EM) | Any prediction completed |
| `process_predicted_model` | process_predicted_model | Model processed for MR |
| `phaser` | phaser | Molecular replacement completed |
| `dock` | dock_in_map | Model docked in map |
| `autobuild` | autobuild | Model building completed |
| `autobuild_denmod` | autobuild_denmod | Density modification completed |
| `autosol` | autosol | Experimental phasing completed |
| `ligandfit` | ligandfit | Ligand fitted |
| `resolve_cryo_em` | resolve_cryo_em | Map optimization completed |
| `map_sharpening` | map_sharpening | Map sharpening completed |
| `map_to_model` | map_to_model | De novo model building completed |
| `map_symmetry` | map_symmetry | Map symmetry analysis completed |

**2. `done_tracking` blocks** (in programs.yaml)

Each program's `done_tracking` block defines its workflow done flag and tracking strategy:

- **`strategy: "set_flag"`** (default) — sets a boolean done flag on success. Most programs use this.
- **`strategy: "run_once"`** — sets the done flag AND filters the program from the valid list after first successful run. Used by `phenix.xtriage`, `phenix.mtriage`, and `phenix.map_symmetry`.
- **`strategy: "count"`** — sets the done flag AND increments a counter (e.g., `refine_count`). Used by `phenix.refine`, `phenix.real_space_refine`, and `phenix.phaser`.

All detection is driven by `history_detection.markers` (substring matching). Programs like `phenix.refine` use `exclude_markers: ["real_space"]` to prevent false matches with `phenix.real_space_refine`. The only program requiring Python-only tracking is `phenix.predict_and_build`, which cascades flags across programs.

**Programs that run multiple times** (intentionally):
- `phenix.refine` / `phenix.real_space_refine` - Iterative refinement
- `phenix.molprobity` - Validation after each cycle
- `phenix.polder` - Can run for different sites

### metrics.yaml

Quality thresholds and display configuration:

```yaml
metrics:
  r_free:
    display_name: "R-free"
    format: "{value:.4f}"
    direction: lower_is_better
    thresholds:
      excellent: 0.20
      good: 0.25
      acceptable: 0.30
    
  map_cc:
    display_name: "Map CC"
    format: "{value:.3f}"
    direction: higher_is_better
    thresholds:
      excellent: 0.80
      good: 0.75
      acceptable: 0.70
```

---

## Data Flow

### Single Cycle

```
1. INPUT
   • available_files: [data.mtz, seq.fa, refine_001.pdb]
   • history: [{cycle: 1, program: "xtriage", ...}]
   • log_text: "R-free = 0.28..."
   • user_advice: "continue refinement"

2. PERCEIVE
   • Categorize files → {mtz: [data.mtz], model: [refine_001.pdb]}
   • Detect state → "refine" phase, valid: [refine, ligandfit, molprobity]
   • Extract metrics → {r_free: 0.28, r_work: 0.24}
   • Analyze trend → "improving, no plateau"
   • Emit: STATE_DETECTED, METRICS_EXTRACTED, METRICS_TREND
   
3. PLAN
   • Check directives → no stop condition met
   • Call LLM or rules → phenix.refine
   • Validate → ✓ in valid_programs
   • Emit: PROGRAM_SELECTED (with full reasoning)
   
4. BUILD
   • Select model: refine_001.pdb (from best_files)
   • Select mtz: data.mtz (rfree_locked)
   • Build: "phenix.refine refine_001.pdb data.mtz output.prefix=refine_002"
   • Emit: FILES_SELECTED, COMMAND_BUILT

5. OUTPUT
   • program: "phenix.refine"
   • command: "phenix.refine ..."
   • events: [{type: "state_detected", ...}, ...]
   • stop: false
```

---

## Execution Modes

### Standard Mode (LLM)

```bash
phenix.ai_agent original_files="data.mtz seq.fa"
```

- LLM makes program selection decisions
- Rules engine validates all choices
- Full reasoning captured in events

### Rules-Only Mode

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz seq.fa"
```

- Deterministic program selection (first valid program)
- No LLM calls
- Faster, reproducible for testing

### Dry-Run Mode

```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic
```

- Simulated execution with predefined outcomes
- For testing workflow logic without running PHENIX programs

---

## Error Handling

### User Request Invalid

When user requests an unavailable program:

```
============================================================
  WARNING: Requested program not available
============================================================
  You requested: phenix.refine
  Reason: Not valid in current workflow state 'xray_initial'
  Running instead: phenix.xtriage
  Available programs: phenix.xtriage
  Suggestion: This program requires different conditions.
============================================================
```

The agent:
1. Detects that user mentioned the program in their advice
2. Explains WHY it's not available
3. Suggests what will run instead
4. Always shown (QUIET verbosity level)

### Sanity Checks

Red flags that indicate problems:
- **Experiment type changed** mid-workflow → abort
- **R-free spike** (increased > 0.05) → warning
- **No model** available for refinement → abort
- **Resolution unknown** before refinement → warning

---

## Session Management

The agent provides two mechanisms for inspecting and modifying an existing session
**without running new crystallographic cycles**.

### Viewing a Session

```bash
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=basic
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=detailed
```

`basic` prints a one-line-per-cycle summary table (program, R-free, result).
`detailed` prints full reasoning and command for every cycle.
Both modes populate `self.result` identically to a normal run so GUI calls
(`get_results()`, `get_results_as_JSON()`) work without special cases.

`restart_mode=resume` is automatically set when either session management
parameter is active — no manual flag required.

### Removing Cycles

```bash
phenix.ai_agent log_directory=AIAgent_run1 remove_last_n=2
```

Removes the last N cycles from the session, clears the stale AI summary,
rebuilds `active_files.json` and `best_files` from remaining history,
and saves. Useful for pruning a failed run before re-running.

### Extending a Completed Workflow with New Advice

When a workflow has fully completed, resuming with new `project_advice`
triggers follow-up programs via the Q1 mechanism:

```bash
phenix.ai_agent \
    log_directory=AIAgent_run1 \
    restart_mode=resume \
    project_advice="also run polder on chain B residue 100"
```

1. New advice hash detected → `advice_changed=True`
2. PERCEIVE steps `complete` phase back to `validate` (adds polder, molprobity, etc.)
3. PLAN suppresses AUTO-STOP for one cycle
4. LLM acts on the new advice; after success, normal termination resumes

See [USER_DIRECTIVES.md](guides/USER_DIRECTIVES.md#extending-a-completed-workflow)
and [ARCHITECTURE.md](reference/ARCHITECTURE.md#advice-change-detection) for details.

---

## Testing

### Test Suites (34 files, 32 in runner)

**Standalone (no PHENIX required, 25 in runner):**
- API Schema, Best Files Tracker, Transport, State Serialization
- Command Builder, File Categorization, File Utils
- Session Summary, Session Directives, Session Tools, Audit Fix Regressions
- Advice Preprocessing, Directive Extractor, Directive Validator, Directives Integration
- Event System, Metric Patterns, Pattern Manager
- Program Registration, Summary Display, New Programs
- Error Analyzer, Decision Flow, Phaser Multimodel
- History Analysis, Docs Tools, YAML Tools

**PHENIX-dependent (7 in runner):**
- Workflow State, YAML Config, Sanity Checker
- Metrics Analyzer, Dry Run, Integration, Directives Integration

**Additional test files (not in run_all_tests.py):**
- tst_template.py (template builder), tst_utils.py (assert helpers)

### Running Tests

```bash
python3 tests/run_all_tests.py        # All tests
python3 tests/run_all_tests.py --quick  # Standalone only
python3 tests/tst_event_system.py    # Single suite
```

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v112.31 | **Session management**: `display_and_stop` / `remove_last_n` populate `self.result`; `get_results()` safe before `run()`; `restart_mode` auto-set; **Q1**: resuming with new advice after workflow completion steps back from `complete` to `validate` phase, enabling follow-up programs (polder etc.) |
| v112 | **Steps table metrics**: cycle metrics as primary source; benign warning metrics extraction; ligand typing fix; case-sensitive pattern fix; autobuild_denmod detection; YAML log_parsing for 8 programs |
| v111 | **Summary output fixes**: predict_and_build R-free extraction; ligandfit output in final file list; fallback cycle status check fix |
| v110 | **Stepwise mode**: automation_path controls predict_and_build behavior; fallback program tracking; autobuild scoring equals refined; best files in summary fix |
| v40 | Fixes 12-21: Ligandfit MTZ exclusion, stop condition on failed runs, summary display, predict_and_build resolution handling; USER_REQUEST_INVALID event |
| v39 | Event system plumbing fixes for single-shot mode |
| v38 | Event system Phase 4: display integration |
| v36-37 | Event system Phases 2-3: instrumentation and transport |
| v34 | Event system Phase 1: EventLog and EventFormatter |
| v30-33 | YAML centralization, BestFilesTracker, CommandBuilder unification |

---

## Automation Modes

The agent supports two automation modes controlled by `maximum_automation`:

### Automated Mode (default)

```bash
phenix.ai_agent maximum_automation=True original_files="data.mtz sequence.fa"
```

- `predict_and_build` runs the complete workflow (prediction → MR → building)
- Fewer checkpoints, faster end-to-end processing
- Best for well-understood datasets

### Stepwise Mode

```bash
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"
```

- `predict_and_build` stops after prediction only (`stop_after_predict=True`)
- User can inspect predicted model before proceeding
- Workflow continues: `process_predicted_model` → `phaser` → `refine`
- Best for troubleshooting or when intermediate inspection is needed

The `automation_path` is set in workflow_state and propagated to all decision-making components to ensure consistent behavior throughout the pipeline.
