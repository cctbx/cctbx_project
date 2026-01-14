# PHENIX AI Agent Decision-Making Architecture

## Overview

The PHENIX AI Agent automates crystallographic and cryo-EM structure determination workflows. It uses a YAML-driven architecture where domain knowledge (programs, workflows, metrics) is defined in configuration files, while Python code provides a generic execution engine.

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                          YAML CONFIGURATION                         │
│  ┌─────────────┐  ┌─────────────────┐  ┌─────────────────────────┐  │
│  │ programs.yaml│  │ workflows.yaml │  │     metrics.yaml       │  │
│  │             │  │                 │  │                         │  │
│  │ - inputs    │  │ - phases        │  │ - thresholds            │  │
│  │ - outputs   │  │ - transitions   │  │ - resolution-dependent  │  │
│  │ - commands  │  │ - conditions    │  │ - quality assessment    │  │
│  └─────────────┘  └─────────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                         EXECUTION ENGINE                            │
│                                                                     │
│   ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌─────────────┐  │
│   │ PERCEIVE │───►│   PLAN   │───►│  BUILD   │───►│  VALIDATE   │  │
│   └──────────┘    └──────────┘    └──────────┘    └──────┬──────┘  │
│        ▲                                                 │         │
│        │              ┌──────────┐                       │         │
│        │         ┌───►│ FALLBACK │◄──────────────────────┤         │
│        │         │    └──────────┘    (on 3 failures)    │         │
│        │         │                                       │         │
│        │         ▼                                       ▼         │
│   ┌──────────┐◄──────────────────────────────────┌───────────┐    │
│   │  OUTPUT  │                                   │  EXECUTE  │    │
│   └──────────┘◄──────────────────────────────────└───────────┘    │
│        │                                                           │
│        ▼                                                           │
│    [END or LOOP]                                                   │
└─────────────────────────────────────────────────────────────────────┘
```

## Node Responsibilities

| Node | Responsibility |
|------|----------------|
| **PERCEIVE** | Analyze history, categorize files, detect workflow state, extract metrics |
| **PLAN** | Select next program (via LLM or RulesSelector) |
| **BUILD** | Construct command with files and parameters |
| **VALIDATE** | Check program validity, file existence, duplicate prevention |
| **EXECUTE** | Run PHENIX command (or simulate in dry run mode) |
| **FALLBACK** | Mechanical backup when primary selection fails |
| **OUTPUT** | Format results, determine if workflow continues |

## Decision-Making Layers

### Layer 1: Workflow State Detection

The `WorkflowEngine` reads `workflows.yaml` to determine:

- **Current phase** based on history and available files
- **Valid programs** for this phase
- **Transition conditions** to next phase

```python
from agent.workflow_engine import WorkflowEngine

engine = WorkflowEngine()
state = engine.get_workflow_state(history, files, experiment_type)
# Returns: {state_name, valid_programs, phase, reason}
```

### Layer 2: Program Selection

Two modes available:

**LLM Mode** (default):
- Constructs prompt with state, recommendations, constraints
- LLM chooses from valid programs
- Falls back to rules on failure

**Rules Mode** (`use_rules_only=True`):
- `RulesSelector` uses workflows.yaml priorities
- Fully deterministic, no API calls
- Selects highest-priority valid program

```python
from agent.rules_selector import RulesSelector

selector = RulesSelector()
intent = selector.select_next_action(workflow_state, files, metrics_trend)
# Returns: {program, files, strategy, reasoning}
```

### Layer 3: Metrics Evaluation

The `MetricEvaluator` reads `metrics.yaml` to assess:

- **Quality level** (good/acceptable/poor)
- **Improvement** between cycles
- **Plateau detection** for stopping
- **Resolution-dependent thresholds**

```python
from agent.metric_evaluator import MetricEvaluator

evaluator = MetricEvaluator()
quality = evaluator.assess_quality({"r_free": 0.25}, resolution=2.0)
trend = evaluator.analyze_trend(history, "xray", resolution=2.0)
```

### Layer 4: Validation Gate

Critical safety mechanism requiring validation before stopping:

**Conditions requiring validation:**
- R-free below success threshold
- 3+ refinement cycles completed
- Model quality is "good"

**Behavior:**
- STOP removed from valid_programs
- molprobity added at highest priority
- After validation: STOP restored

### Layer 5: Command Building

The `TemplateBuilder` reads `programs.yaml` to:

1. Match files to program inputs by extension
2. Apply required flags from YAML
3. Handle optional parameters
4. Build final command string

## Operating Modes

### Standard Mode

```bash
phenix.ai_agent original_files="data.mtz seq.fa"
```

- LLM selects programs
- YAML defines constraints
- Full automation

### Rules-Only Mode

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz seq.fa"
```

- No LLM/API calls
- Deterministic selection
- Offline-capable

### Dry Run Mode

```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic use_rules_only=True
```

- Simulates program execution
- Uses pre-recorded logs/outputs
- Full workflow testing

## Configuration Files

### programs.yaml Structure

```yaml
phenix.program_name:
  description: "What it does"
  category: analysis|building|refinement|validation
  experiment_types: [xray, cryoem]
  
  inputs:
    required:
      slot_name:
        extensions: [.ext1, .ext2]
        flag: "parameter="
    optional:
      slot_name:
        extensions: [.ext]
        flag: "optional_param="
  
  outputs:
    files:
      - pattern: "*.pdb"
        type: model
    metrics: [metric1, metric2]
  
  command: "phenix.program {slot1} {slot2}"
  
  log_parsing:
    metric_name:
      pattern: 'regex pattern'
      type: float|int|str
```

### workflows.yaml Structure

```yaml
experiment_type:
  description: "Workflow description"
  
  phases:
    phase_name:
      description: "Phase description"
      programs:
        - program: phenix.program
          preferred: true
          conditions:
            - has: file_type
      transitions:
        on_complete: next_phase
        on_target_reached: validate
        on_plateau: validate
    
    complete:
      stop: true
      stop_reasons: ["Success message"]
```

### metrics.yaml Structure

```yaml
metric_name:
  description: "What it measures"
  direction: minimize|maximize
  thresholds:
    good: 0.25
    acceptable: 0.30
  resolution_dependent:
    - range: [0, 1.5]
      good: 0.20
    - range: [1.5, 2.5]
      good: 0.25
```

## Stop Conditions

| Condition | Trigger | Action |
|-----------|---------|--------|
| Success | Metrics at target + validated | STOP |
| Plateau | <0.5% improvement for 2 cycles | STOP |
| Excessive | 5+ refinement cycles | STOP |
| Validation Required | Good metrics, no validation | Force molprobity |
| All Duplicates | No non-duplicate commands | STOP |

## Test Scenarios

Located in `tests/scenarios/`:

| Scenario | Workflow | Programs |
|----------|----------|----------|
| `xray_basic` | xtriage → predict_and_build → refine(×3) → molprobity | 6 steps |
| `cryoem_basic` | mtriage → predict_and_build → real_space_refine(×3) → molprobity | 6 steps |

### Scenario Structure

```
scenario_name/
├── scenario.yaml          # Configuration
├── inputs/                # Initial files
│   ├── data.mtz
│   └── sequence.fa
└── steps/
    ├── phenix.program_1/
    │   ├── log.txt        # Simulated output
    │   └── outputs/       # Output files
    └── ...
```

## Key Classes

| Class | File | Purpose |
|-------|------|---------|
| `WorkflowEngine` | `agent/workflow_engine.py` | YAML workflow interpreter |
| `RulesSelector` | `agent/rules_selector.py` | Deterministic program selection |
| `MetricEvaluator` | `agent/metric_evaluator.py` | YAML-based metrics |
| `ProgramRegistry` | `agent/program_registry.py` | Program definitions |
| `TemplateBuilder` | `agent/template_builder.py` | Command construction |
| `DryRunManager` | `agent/dry_run_manager.py` | Test scenario management |

## Utilities

### yaml_tools.py

```bash
# Validate YAML files
python yaml_tools.py validate

# Display formatted contents
python yaml_tools.py display knowledge/programs.yaml

# Show summary of all configurations
python yaml_tools.py summary
```

### run_all_tests.py

```bash
# Run all test suites
python run_all_tests.py

# Individual suites
python tests/test_dry_run.py
python tests/test_integration.py
```

## Future Enhancements

1. **More test scenarios** - SAD phasing, ligand fitting, twinning
2. **Better plateau detection** - Use multiple metrics
3. **Adaptive thresholds** - Learn from successful workflows
4. **Parallel execution** - Run independent programs simultaneously
5. **Web interface** - Visual workflow monitoring
