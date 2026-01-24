# PHENIX AI Agent Logic Documentation

This document describes the decision logic used by the PHENIX AI Agent.

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [YAML Configuration](#yaml-configuration)
3. [Workflow States](#workflow-states)
4. [Program Selection](#program-selection)
5. [Best Files Tracking](#best-files-tracking)
6. [Metrics and Stop Conditions](#metrics-and-stop-conditions)
7. [Operating Modes](#operating-modes)
8. [Modifying Agent Behavior](#modifying-agent-behavior)

---

## Architecture Overview

The PHENIX AI Agent is a YAML-driven workflow automation system for crystallographic and cryo-EM structure determination. It operates as a cyclic graph:

```
PERCEIVE → PLAN → BUILD → VALIDATE → EXECUTE → (loop or stop)
```

### Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            YAML CONFIGURATION                               │
│  ┌───────────────┐ ┌───────────────┐ ┌─────────────┐ ┌───────────────────┐  │
│  │ programs.yaml │ │workflows.yaml │ │metrics.yaml │ │file_categories.yaml│ │
│  │               │ │               │ │             │ │                   │  │
│  │ - inputs      │ │ - phases      │ │ - thresholds│ │ - extensions      │  │
│  │ - outputs     │ │ - transitions │ │ - resolution│ │ - patterns        │  │
│  │ - commands    │ │ - conditions  │ │   dependent │ │ - subcategories   │  │
│  │ - invariants  │ │ - priority_   │ │ - quality   │ │                   │  │
│  │ - input_      │ │   when        │ │   assessment│ │                   │  │
│  │   priorities  │ │               │ │             │ │                   │  │
│  │ - user_advice_│ │               │ │             │ │                   │  │
│  │   keywords    │ │               │ │             │ │                   │  │
│  └───────────────┘ └───────────────┘ └─────────────┘ └───────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            PYTHON ENGINE                                    │
│  ┌─────────────┐   ┌─────────────┐   ┌─────────────┐   ┌─────────────┐     │
│  │ PERCEIVE    │ → │ PLAN        │ → │ BUILD       │ → │ VALIDATE    │     │
│  │             │   │             │   │             │   │             │     │
│  │ - Parse logs│   │ - Get state │   │ - Select    │   │ - Check     │     │
│  │ - Categorize│   │ - Check     │   │   files     │   │   program   │     │
│  │   files     │   │   conditions│   │ - Apply     │   │ - Enforce   │     │
│  │ - Build     │   │ - Recommend │   │   invariants│   │   workflow  │     │
│  │   metrics   │   │   programs  │   │ - Build cmd │   │ - Fallback  │     │
│  └─────────────┘   └─────────────┘   └─────────────┘   └─────────────┘     │
│                                                               │             │
│                                                               ▼             │
│                                                        ┌─────────────┐      │
│                                                        │ EXECUTE     │      │
│                                                        │ or STOP     │      │
│                                                        └─────────────┘      │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Key Design Principles

1. **Configuration over Code** - Crystallography knowledge lives in YAML files, not Python
2. **LLM Optional** - Can run fully deterministically with rules-based selection
3. **Backward Compatible** - Falls back to hardcoded logic if YAML fails
4. **Testable** - Dry run mode allows full workflow testing without PHENIX
5. **Layered Decision Flow** - Clear separation between workflow engine, LLM planning, and stop conditions

### Decision Flow Layers

The agent uses a clean, layered architecture for decision-making:

| Layer | Component | Responsibility |
|-------|-----------|----------------|
| 1 | Workflow Engine | Determine valid_programs, apply directives |
| 2 | LLM/Rules Plan | Select program from valid_programs |
| 3 | Build & Execute | Build command, run program |
| 4 | Post-Execution | Check stop conditions (ONLY place) |

**Key principle**: Stop conditions are only evaluated AFTER program execution, not during planning.

---

## YAML Configuration

All domain knowledge is defined in four YAML files:

### knowledge/programs.yaml

Defines every program the agent can use:

```yaml
phenix.refine:
  description: "Refine atomic model against X-ray data"
  category: refinement
  experiment_types: [xray]

  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
      mtz:
        extensions: [.mtz]
        flag: ""

  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
    metrics:
      - r_work
      - r_free

  command: "phenix.refine {model} {mtz}"

  # Invariants: auto-fix rules applied before command execution
  invariants:
    - name: resolution_when_building
      check:
        any_of:
          - has_strategy: [resolution]
          - strategy_equals: {stop_after_predict: true}
      fix:
        auto_fill_resolution: true

  # Input priorities: which file categories to prefer
  input_priorities:
    model:
      categories: [refined, phaser_output, with_ligand]
      exclude_categories: [predicted]
    mtz:
      categories: [refined_mtz, mtz]

  # User advice keywords: trigger this program when user mentions these
  user_advice_keywords:
    - "refinement"
    - "refine"
```

### knowledge/workflows.yaml

Defines the state machine for each experiment type:

```yaml
xray:
  phases:
    analyze:
      programs:
        - program: phenix.xtriage
      transitions:
        on_complete: obtain_model

    obtain_model:
      programs:
        - program: phenix.predict_and_build
          preferred: true
          conditions:
            - has: sequence
        - program: phenix.autosol
          conditions:
            - has: sequence
            - has: anomalous_data
          priority_when: strong_anomalous  # Prioritize when anomalous signal strong
        - program: phenix.phaser
          conditions:
            - has: search_model

    refine:
      programs:
        - program: phenix.refine
      transitions:
        on_target_reached: validate
        on_plateau: validate

    validate:
      programs:
        - program: phenix.molprobity
      transitions:
        if_quality_acceptable: complete

    complete:
      stop: true
```

### knowledge/metrics.yaml

Defines metrics, thresholds, and quality assessment:

```yaml
r_free:
  description: "Cross-validation R-factor"
  direction: minimize
  thresholds:
    good: 0.25
    acceptable: 0.30
  resolution_dependent:
    - range: [0, 1.5]
      good: 0.20
    - range: [1.5, 2.5]
      good: 0.25
    - range: [2.5, 4.0]
      good: 0.30

anomalous_measurability:
  description: "Measurability of anomalous signal (0-1 scale)"
  direction: maximize
  thresholds:
    good: 0.10       # Strong signal for SAD/MAD
    acceptable: 0.05  # Moderate signal
```

### knowledge/file_categories.yaml

Defines how files are categorized by name and extension:

```yaml
# Primary categories (by extension)
mtz:
  description: "X-ray reflection data files"
  extensions: [".mtz", ".sca", ".hkl", ".sdf"]

pdb:
  description: "Atomic model coordinate files"
  extensions: [".pdb"]

# Subcategories (by filename pattern)
refined:
  description: "X-ray refined models"
  subcategory_of: pdb
  patterns:
    - "*refine*"
  excludes:
    - "*real_space*"
    - "*rsr*"

rsr_output:
  description: "Real-space refined models (cryo-EM)"
  subcategory_of: pdb
  patterns:
    - "*real_space_refined*"
    - "*rsr_*"
```

---

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
              ↓
       xray_has_model
              ↓
         refine (loop) ←──────────────────────┐
              ↓                               │
       xray_refined                           │
         ↓    ↓    ↓                          │
    molprobity refine STOP                    │
              │                               │
              ↓                               │
       [if ligand] ligandfit → pdbtools ──────┘
```

† Prioritized when strong anomalous signal detected (measurability > 0.10)

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `xray_initial` | Starting point | xtriage |
| `xray_analyzed` | After data analysis | predict_and_build, phaser, autosol |
| `xray_has_prediction` | Have AlphaFold model | process_predicted_model |
| `xray_has_phases` | After experimental phasing | autobuild |
| `xray_has_model` | Have placed model | refine |
| `xray_refined` | After refinement | refine, molprobity, autobuild, STOP |

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
| `cryoem_analyzed` | After map analysis | predict_and_build, dock_in_map |
| `cryoem_has_model` | Model in map | real_space_refine |
| `cryoem_refined` | After refinement | real_space_refine, molprobity, STOP |

---

## Program Selection

### Tiered Decision System

| Tier | Type | Description | Example |
|------|------|-------------|---------|
| 1 | Hard Constraints | Must follow | xtriage before phaser |
| 2 | Strong Defaults | Auto-applied | Generate R-free flags on first refine |
| 3 | Soft Guidance | Suggestions | Consider anisotropic B at high res |

### Rules-Based Selection (RulesSelector)

When operating without LLM, the agent uses `RulesSelector` which:

1. Gets valid programs for current workflow state
2. Checks `priority_when` conditions (e.g., strong_anomalous → autosol)
3. Applies `user_advice_keywords` if user gave advice
4. Ranks by priority defined in workflows.yaml
5. Selects highest-priority program with satisfied conditions
6. Auto-selects files using `input_priorities` from programs.yaml
7. Applies `invariants` to auto-fix strategy options

### LLM-Based Selection

When using LLM, the agent:

1. Constructs prompt with workflow state and recommendations
2. LLM chooses from valid programs
3. BUILD node applies invariants and defaults
4. VALIDATE ensures choice is valid
5. Falls back to rules if LLM fails 3 times

---

## Best Files Tracking

The agent maintains a **Best Files Tracker** that identifies and tracks the highest-quality file of each type throughout a session. This ensures programs always receive optimal inputs without relying on filename patterns or recency alone.

### Configuration

Scoring configuration is defined in `knowledge/metrics.yaml` under the `best_files_scoring` section. This follows the "configuration over code" principle, making scoring adjustable without code changes.

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

Each file is scored based on **processing stage** and **quality metrics**. Higher scores are better. All scoring parameters are configurable in `knowledge/metrics.yaml`.

#### Formula Types

The YAML configuration supports three formula types for metric scoring:

| Formula | Description | Use Case |
|---------|-------------|----------|
| `linear` | Higher value is better | map_cc (correlation) |
| `linear_inverse` | Lower value is better | r_free, clashscore, resolution |
| `boolean` | True/False flag | has_rfree_flags |

#### Model Scoring (0-200 points)

**Stage Score (0-100 points):**

| Stage | Score | Description |
|-------|-------|-------------|
| `refined` | 100 | Output from phenix.refine |
| `rsr_output` | 100 | Output from real_space_refine |
| `autobuild_output` | 80 | Output from phenix.autobuild |
| `docked` | 60 | Output from dock_in_map |
| `processed_predicted` | 50 | Processed AlphaFold model |
| `predicted` | 40 | Raw AlphaFold prediction |
| `phaser_output` | 30 | MR solution from Phaser |
| `pdb` | 10 | Unknown/generic PDB |

**Metrics Score (0-100 points):**

| Metric | Max Points | Formula | Best Value | Worst Value |
|--------|------------|---------|------------|-------------|
| R-free | 40 | linear_inverse | 0.20 | 0.40 |
| Map CC | 30 | linear | 1.0 | 0.0 |
| Clashscore | 30 | linear_inverse | 0 | 20 |

**Example:** A refined model (100) with R-free=0.22 (36), map_cc=0.75 (22.5), clashscore=8 (18) = **176.5 points**

#### Map Scoring (0-150 points)

**Stage Score (0-100 points):**

| Stage | Score | Description |
|-------|-------|-------------|
| `optimized_full_map` | 100 | Output from resolve_cryo_em |
| `sharpened` | 90 | Output from map_sharpening |
| `density_modified` | 80 | Density-modified map |
| `full_map` | 50 | User-provided full map |
| `half_map` | 10 | Half-map (not suitable as primary) |

**Resolution Bonus (0-30 points):** linear_inverse formula with best=1.0Å, worst=4.0Å

#### MTZ Scoring (Special Rules)

MTZ files have a unique rule: **the first MTZ with R-free flags locks forever**. This ensures consistent R-free statistics throughout refinement.

| Stage | Score |
|-------|-------|
| `refined_mtz` | 70 |
| `original` | 50 |
| Has R-free flags | +30 |

Once an MTZ with R-free flags is identified, no other MTZ can replace it, regardless of score.

### R-free MTZ Locking (X-ray Only)

In X-ray crystallography, R-free flags are used for cross-validation. Once generated, the **same flags must be used for all subsequent refinements**. Using different R-free flags invalidates the R-free statistic.

The agent handles this with explicit R-free MTZ tracking at the session level:

```
Cycle 5: phenix.refine generates R-free flags
         → Output: refine_001_data.mtz
         → Session locks: rfree_mtz = refine_001_data.mtz

Cycle 6+: All refinements MUST use refine_001_data.mtz
          → Locked MTZ has HIGHEST priority (Priority 0)
          → Overrides best_files, categorized_files, etc.
```

This is tracked in `session.json`:
```json
{
  "rfree_mtz": "/path/to/refine_001_data.mtz",
  "rfree_mtz_locked_at_cycle": 5
}
```

### File Selection Priority

When building commands, the agent selects files in this priority order:

0. **Locked R-free MTZ** - For MTZ inputs in X-ray refinement (highest priority)
1. **Best Files** - From the BestFilesTracker
2. **Categorized Files** - From workflow_state file categorization
3. **Extension Matching** - Search available_files by extension

```
BUILD node file selection for phenix.refine:
┌─────────────────────────────────────────────────────────┐
│ Need: mtz input for phenix.refine                       │
├─────────────────────────────────────────────────────────┤
│ 0. Check session rfree_mtz (LOCKED)                     │
│    → /path/to/refine_001_data.mtz ✓ SELECTED           │
├─────────────────────────────────────────────────────────┤
│ 1-3. (Skipped - locked MTZ takes precedence)            │
└─────────────────────────────────────────────────────────┘
```

### Automatic Updates

The tracker updates automatically in `Session.record_result()`:

1. **After successful cycles**: All output files are evaluated and scored
2. **After validation**: Validation metrics (clashscore, etc.) update the best model's score
3. **Stage inference**: The program name determines the output stage

```python
# Automatic flow in record_result():
phenix.refine succeeds
  → Output: refine_001_001.pdb
  → Stage inferred: "refined" (100 pts)
  → Metrics extracted: r_free=0.25 (30 pts)
  → Total score: 130
  → Becomes best_model (if score > previous)

phenix.molprobity succeeds
  → Metrics extracted: clashscore=8 (18 pts)
  → Updates best_model score: 130 → 148
```

### History Tracking

Every change to "best" is recorded with:
- Old and new file paths
- Old and new scores
- Cycle number when change occurred
- Reason for the change

```python
# Example history entry
BestFileChange(
    category="model",
    old_path="/path/to/docked.pdb",
    new_path="/path/to/refined.pdb",
    old_score=60.0,
    new_score=130.0,
    cycle=3,
    reason="Higher score (130.0 > 60.0)"
)
```

### Session Persistence

Best files state is saved in `session.json` and survives restarts:

```json
{
  "best_files": {
    "best": {
      "model": {
        "path": "/path/to/refine_002_001.pdb",
        "category": "model",
        "stage": "refined",
        "score": 165.5,
        "metrics": {"r_free": 0.22, "clashscore": 8},
        "cycle": 5,
        "reason": "Better metrics"
      },
      "map": {
        "path": "/path/to/denmod_map.ccp4",
        "stage": "optimized_full_map",
        "score": 125.0,
        "cycle": 2
      }
    },
    "mtz_with_rfree_locked": true
  }
}
```

### Debug Logging

The PERCEIVE node logs best files status:

```
PERCEIVE: Best files: model=refine_002_001.pdb, map=denmod_map.ccp4, mtz=data.mtz
```

The BUILD node logs when best files are used:

```
BUILD: Using best_model for slot 'model': refine_002_001.pdb
BUILD: Using best_map for slot 'map': denmod_map.ccp4
```

### Key Implementation Files

| File | Role |
|------|------|
| `agent/best_files_tracker.py` | Core tracker class with scoring |
| `agent/session.py` | Integration with session persistence |
| `agent/template_builder.py` | Uses best_files for command building |
| `agent/graph_nodes.py` | Passes best_files through workflow |

---

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
3. **Excessive**: Too many refinement cycles (default: 5)
4. **Validation Gate**: Must run molprobity before stopping if R-free is good

### Validation Gate

A critical safety mechanism prevents stopping without validation:

- If R-free below success threshold → validation required
- If 3+ refinement cycles completed → validation required
- STOP is removed from valid_programs until validation runs

---

## Operating Modes

### Standard Mode (LLM + YAML)

```bash
phenix.ai_agent original_files="data.mtz sequence.fa"
```

- Uses LLM for program selection
- YAML defines valid programs and constraints
- Best for complex or unusual cases

### Rules-Only Mode (No LLM)

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz sequence.fa"
```

- Fully deterministic, no API calls
- Uses RulesSelector from workflows.yaml
- Fast, reproducible, offline-capable

### Dry Run Mode (Testing)

```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic use_rules_only=True
```

- Uses pre-recorded logs and outputs
- No actual PHENIX programs run
- Tests full workflow without real data

---

## Modifying Agent Behavior

### To Add a New Program

Adding a new program requires only **2-3 files**. See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide.

**Minimal example:**

Edit `knowledge/programs.yaml`:

```yaml
phenix.my_program:
  description: "What it does"
  category: analysis
  experiment_types: [xray]
  run_once: true  # Auto-creates my_program_done flag

  inputs:
    required:
      data:
        extensions: [.mtz]
        flag: ""

  command: "phenix.my_program {data}"

  # Metric extraction patterns (used by log_parsers.py and session.py)
  log_parsing:
    my_metric:
      pattern: 'My Metric[:\s]+([0-9.]+)'
      type: float
```

Edit `knowledge/workflows.yaml`:

```yaml
xray:
  phases:
    analyze:
      programs:
        - program: phenix.my_program
          conditions:
            - not_done: my_program
```

The system automatically handles metric extraction, tracking flags, and summary display.

### To Change Workflow Logic

Edit `knowledge/workflows.yaml`:

```yaml
xray:
  phases:
    my_new_phase:
      programs:
        - program: phenix.my_program
          conditions:
            - has: sequence
          priority_when: strong_anomalous
      transitions:
        on_complete: refine
```

### To Adjust Thresholds

Edit `knowledge/metrics.yaml`:

```yaml
r_free:
  thresholds:
    good: 0.22  # More stringent
    acceptable: 0.28
```

### To Add File Category

Edit `knowledge/file_categories.yaml`:

```yaml
my_output:
  description: "Output from my_program"
  subcategory_of: pdb
  patterns:
    - "*my_program*"
```

### To Validate Changes

```bash
python agent/yaml_tools.py validate
python agent/yaml_tools.py terms --detail full
python agent/yaml_tools.py summary
```

---

## Key Files Reference

| File | Purpose |
|------|---------|
| `knowledge/programs.yaml` | Program definitions |
| `knowledge/workflows.yaml` | Workflow state machines |
| `knowledge/metrics.yaml` | Metric thresholds |
| `knowledge/file_categories.yaml` | File categorization rules |
| `agent/rules_selector.py` | Deterministic program selection |
| `agent/workflow_engine.py` | YAML workflow interpreter, `_apply_directives()` |
| `agent/workflow_state.py` | File categorization, state detection |
| `agent/metric_evaluator.py` | YAML-based quality assessment |
| `agent/template_builder.py` | Command construction with invariants |
| `agent/program_registry.py` | Program info access |
| `agent/dry_run_manager.py` | Test scenario management |
| `agent/directive_extractor.py` | Extract directives from user advice |
| `agent/directive_validator.py` | Apply program settings to intent |
| `agent/graph_nodes.py` | LangGraph nodes including `plan()` |
| `programs/ai_agent.py` | Main agent, post-execution stop check |
| `phenix_ai/log_parsers.py` | Program output parsing |
| `agent/yaml_tools.py` | YAML validation and inspection |
| `agent/session_tools.py` | Session management CLI |
| `agent/docs_tools.py` | Documentation generation |

---

## Testing

### Run All Tests

```bash
python run_all_tests.py
```

### Test Suites

| Suite | Description |
|-------|-------------|
| `test_yaml_config.py` | Configuration validation, invariants, input_priorities, user_advice_keywords |
| `test_workflow_state.py` | State detection, file categorization (including rsr_output) |
| `test_integration.py` | Full workflow paths |
| `test_metrics_analyzer.py` | Quality assessment logic |
| `test_dry_run.py` | Scenario-based testing |
| `test_decision_flow.py` | Decision flow architecture, directive application |
| `test_directive_extractor.py` | Directive extraction from user advice |
| `test_directive_validator.py` | Program settings application |

### Add Test Scenario

Create `tests/scenarios/my_scenario/`:

```
my_scenario/
├── scenario.yaml      # Configuration
├── inputs/            # Initial files
└── steps/
    ├── phenix.program_1/
    │   ├── log.txt
    │   └── outputs/
    └── ...
```

---

## YAML Tools Quick Reference

```bash
# List available YAML files
python agent/yaml_tools.py list

# Validate all configuration
python agent/yaml_tools.py validate

# Display formatted contents
python agent/yaml_tools.py display knowledge/programs.yaml
python agent/yaml_tools.py display knowledge/file_categories.yaml

# Show all defined terms
python agent/yaml_tools.py terms                  # Simple list
python agent/yaml_tools.py terms --detail normal  # With descriptions
python agent/yaml_tools.py terms --detail full    # Full cross-reference

# Show configuration summary
python agent/yaml_tools.py summary

# Compare configurations
python agent/yaml_tools.py compare old_config/ new_config/
```
