# PHENIX AI Agent Logic Documentation

This document describes the decision logic used by the PHENIX AI Agent.

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [YAML Configuration](#yaml-configuration)
3. [Workflow States](#workflow-states)
4. [Program Selection](#program-selection)
5. [Metrics and Stop Conditions](#metrics-and-stop-conditions)
6. [Operating Modes](#operating-modes)
7. [Modifying Agent Behavior](#modifying-agent-behavior)

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

Edit `knowledge/programs.yaml`:

```yaml
phenix.my_program:
  description: "What it does"
  category: analysis
  experiment_types: [xray]

  inputs:
    required:
      data:
        extensions: [.mtz]
        flag: ""

  outputs:
    metrics: [my_metric]

  command: "phenix.my_program {data}"

  # Optional features:
  invariants:
    - name: requires_resolution
      check:
        has_strategy: [resolution]
      fix:
        auto_fill_resolution: true

  input_priorities:
    model:
      categories: [refined]
      exclude_categories: [predicted]

  user_advice_keywords:
    - "my special analysis"
```

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
python yaml_tools.py validate
python yaml_tools.py terms --detail full
python yaml_tools.py summary
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
| `agent/workflow_engine.py` | YAML workflow interpreter |
| `agent/workflow_state.py` | File categorization, state detection |
| `agent/metric_evaluator.py` | YAML-based quality assessment |
| `agent/template_builder.py` | Command construction with invariants |
| `agent/program_registry.py` | Program info access |
| `agent/dry_run_manager.py` | Test scenario management |
| `log_parsers.py` | Program output parsing |
| `yaml_tools.py` | YAML validation and inspection |

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
python yaml_tools.py list

# Validate all configuration
python yaml_tools.py validate

# Display formatted contents
python yaml_tools.py display knowledge/programs.yaml
python yaml_tools.py display knowledge/file_categories.yaml

# Show all defined terms
python yaml_tools.py terms                  # Simple list
python yaml_tools.py terms --detail normal  # With descriptions
python yaml_tools.py terms --detail full    # Full cross-reference

# Show configuration summary
python yaml_tools.py summary

# Compare configurations
python yaml_tools.py compare old_config/ new_config/
```
