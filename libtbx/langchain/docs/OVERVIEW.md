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
│  • Execute graph nodes (perceive → plan → build → validate)                 │
│  • Build response with events                                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            LangGraph Pipeline                                │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐                  │
│  │ perceive │ → │   plan   │ → │  build   │ → │ validate │                  │
│  └────┬─────┘   └────┬─────┘   └────┬─────┘   └────┬─────┘                  │
│       │              │              │              │                         │
│       ▼              ▼              ▼              ▼                         │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                      state["events"] (list)                          │    │
│  │  • STATE_DETECTED, METRICS_EXTRACTED, METRICS_TREND                 │    │
│  │  • PROGRAM_SELECTED, USER_REQUEST_INVALID (if applicable)          │    │
│  │  • FILES_SELECTED, COMMAND_BUILT, STOP_DECISION                    │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Core Components

### 1. Graph Nodes (`agent/graph_nodes.py`)

The agent uses a LangGraph pipeline with four main nodes:

| Node | Purpose | Key Actions |
|------|---------|-------------|
| **perceive** | Understand current state | Categorize files, detect workflow state, extract metrics, analyze trends |
| **plan** | Decide next action | Call LLM or rules engine, validate against workflow, check user directives |
| **build** | Generate command | Select files using BestFilesTracker, build command string |
| **validate** | Check output | Sanity checks, red flag detection |

Each node emits structured events that capture its decision-making process.

### 2. Workflow Engine (`agent/workflow_engine.py`)

Interprets `workflows.yaml` to determine:
- **Current phase**: analyze, obtain_model, refine, validate, complete
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
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution parsed |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `program_selected` | normal | Decision with full reasoning |
| `user_request_invalid` | quiet | User requested unavailable program |
| `files_selected` | verbose | File selection with reasons |
| `command_built` | normal | Final command |
| `stop_decision` | normal | Whether to continue |
| `error` | quiet | Error occurred |

### 5. File Categorization (`agent/file_categorization.py`)

Semantic file classification using rules from `file_categories.yaml`:

| Category | Description | Examples |
|----------|-------------|----------|
| `mtz` | Reflection data | data.mtz, scaled.mtz |
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
- **Locks R-free MTZ** after first refinement for consistency
- Provides `best_files["model"]`, `best_files["mtz"]` to CommandBuilder

### 7. User Directives (`agent/directive_extractor.py`)

Parses natural language guidance from `project_advice`:
- **Stop conditions**: "stop after one refinement", "stop when R-free < 0.25"
- **Workflow preferences**: "skip autobuild", "use SAD phasing"
- **Task focus**: "focus on ligand fitting"

---

## Configuration Files

### programs.yaml

Defines each PHENIX program the agent can run:

```yaml
phenix.refine:
  description: "Crystallographic refinement"
  category: refinement
  experiment_types: [xray]
  run_once: false
  
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
        priority_patterns: [refine]  # Prefer refined models
      mtz:
        extensions: [.mtz]
        flag: ""
    optional:
      sequence:
        extensions: [.fa, .fasta]
        flag: "input.xray_data.labels.seq_file="
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
    metrics:
      - r_free
      - r_work
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-free"
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
        - program: phenix.phaser
          conditions:
            - has: search_model
      transitions:
        on_complete: refine
    
    refine:
      description: "Improve model"
      programs:
        - program: phenix.refine
          preferred: true
        - program: phenix.ligandfit
          conditions:
            - has: ligand_file
            - r_free: "< 0.35"
      repeat:
        max_cycles: 4
        until:
          any:
            - r_free: "< target_r_free"
            - condition: plateau
              cycles: 2
              threshold: 0.005
      transitions:
        on_target_reached: validate
        on_plateau: validate

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

## Testing

### Test Suites (20 total)

**Standalone (no PHENIX required):**
- Event System, Transport, Command Builder
- File Categorization, Best Files Tracker
- Directive Extractor/Validator
- Metric Patterns, Program Registration
- Summary Display

**PHENIX-dependent:**
- Workflow State, Sanity Checker
- Metrics Analyzer, Dry Run, Integration

### Running Tests

```bash
python3 tests/run_all_tests.py        # All tests
python3 tests/run_all_tests.py --quick  # Standalone only
python3 tests/test_event_system.py    # Single suite
```

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v40 | Fixes 12-21: Ligandfit MTZ exclusion, stop condition on failed runs, summary display for failed/prediction runs, predict_and_build resolution handling, None program value handling, test import fixes |
| v40 | USER_REQUEST_INVALID event for invalid program requests |
| v39 | Event system plumbing fixes for single-shot mode |
| v38 | Event system Phase 4: display integration |
| v36-37 | Event system Phases 2-3: instrumentation and transport |
| v34 | Event system Phase 1: EventLog and EventFormatter |
| v30-33 | YAML centralization, BestFilesTracker, CommandBuilder unification |
