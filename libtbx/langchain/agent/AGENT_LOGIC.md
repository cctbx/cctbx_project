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

### Key Design Principles

1. **Configuration over Code** - Crystallography knowledge lives in YAML files, not Python
2. **LLM Optional** - Can run fully deterministically with rules-based selection
3. **Backward Compatible** - Falls back to hardcoded logic if YAML fails
4. **Testable** - Dry run mode allows full workflow testing without PHENIX

---

## YAML Configuration

All domain knowledge is defined in three YAML files:

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
    optional:
      cif:
        extensions: [.cif]
        flag: ""
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
    metrics:
      - r_work
      - r_free
  
  command: "phenix.refine {model} {mtz}"
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
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
```

---

## Workflow States

### X-ray Crystallography Workflow

```
xray_initial → xtriage → xray_analyzed
                              ↓
              [if sequence] predict_and_build → xray_has_prediction
              [if model]    phaser ─────────────────┐
                                                    ↓
                              xray_has_model → refine ←──┐
                                                    ↓    │
                                              xray_refined
                                              ↓    ↓    ↓
                                        molprobity refine STOP
```

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `xray_initial` | Starting point | xtriage |
| `xray_analyzed` | After data analysis | predict_and_build, phaser, autosol |
| `xray_has_prediction` | Have AlphaFold model | process_predicted_model |
| `xray_has_model` | Have placed model | refine |
| `xray_refined` | After refinement | refine, molprobity, autobuild, STOP |

### Cryo-EM Workflow

```
cryoem_initial → mtriage → cryoem_analyzed
                                ↓
                    predict_and_build / dock_in_map
                                ↓
                        cryoem_has_model → real_space_refine ←──┐
                                                    ↓           │
                                              cryoem_refined ───┘
                                              ↓         ↓
                                        molprobity    STOP
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
2. Ranks by priority defined in workflows.yaml
3. Selects highest-priority program with satisfied conditions
4. Auto-selects appropriate files based on extensions
5. Applies strategy options from programs.yaml

### LLM-Based Selection

When using LLM, the agent:

1. Constructs prompt with workflow state and recommendations
2. LLM chooses from valid programs
3. BUILD node applies Tier 2 defaults
4. VALIDATE ensures choice is valid
5. Falls back to rules if LLM fails 3 times

---

## Metrics and Stop Conditions

### Quality Assessment

Metrics are evaluated against YAML-defined thresholds:

| Quality | R-free (2.0Å) | Map CC |
|---------|---------------|--------|
| Good | < 0.25 | > 0.80 |
| Acceptable | < 0.30 | > 0.70 |
| Poor | ≥ 0.30 | ≤ 0.70 |

### Stop Conditions

The workflow stops when:

1. **Success**: Quality metrics reach targets AND validation done
2. **Plateau**: < 0.5% improvement for 2+ consecutive cycles
3. **Excessive**: Too many refinement cycles (default: 5)
4. **Validation Gate**: Must run molprobity before stopping if R-free is good

### Validation Gate

A critical safety mechanism prevents stopping without validation:

- If R-free is below success threshold → validation required
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

Available scenarios:
- `xray_basic` - Standard X-ray workflow
- `cryoem_basic` - Standard cryo-EM workflow

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
  command: "phenix.my_program {data}"
```

### To Change Workflow Logic

Edit `knowledge/workflows.yaml`:

```yaml
xray:
  phases:
    my_new_phase:
      programs:
        - program: phenix.my_program
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

### To Validate Changes

```bash
python yaml_tools.py validate
python yaml_tools.py summary
```

---

## Key Files Reference

| File | Purpose |
|------|---------|
| `knowledge/programs.yaml` | Program definitions |
| `knowledge/workflows.yaml` | Workflow state machines |
| `knowledge/metrics.yaml` | Metric thresholds |
| `agent/rules_selector.py` | Deterministic program selection |
| `agent/workflow_engine.py` | YAML workflow interpreter |
| `agent/metric_evaluator.py` | YAML-based quality assessment |
| `agent/dry_run_manager.py` | Test scenario management |
| `yaml_tools.py` | YAML validation utility |

---

## Testing

### Run All Tests

```bash
python run_all_tests.py
```

### Test Suites

- **Metrics Analyzer** - Quality assessment logic
- **Workflow State** - State detection
- **Integration** - Full workflow paths
- **YAML Config** - Configuration validation
- **Dry Run** - Scenario-based testing

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
