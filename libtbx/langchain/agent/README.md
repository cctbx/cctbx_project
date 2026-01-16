# PHENIX AI Agent V2 - YAML-Driven Crystallographic Workflow Automation

## Overview

The PHENIX AI Agent automates crystallographic and cryo-EM structure determination workflows. It uses a **YAML-driven architecture** where domain knowledge is defined in configuration files, while Python code provides a generic execution engine.

### Key Features

1. **YAML Configuration** - All programs, workflows, metrics, and file categories defined in YAML
2. **Metrics Tracking** - Automatic R-free/CC tracking, plateau detection, auto-stop
3. **Workflow State Machine** - Enforced valid program sequences for X-ray and cryo-EM
4. **Rules-Only Mode** - Fully deterministic operation without LLM
5. **Dry Run Testing** - Test workflows without running actual PHENIX programs

## Configuration Files

The agent uses four YAML configuration files in `knowledge/`:

| File | Purpose |
|------|---------|
| `programs.yaml` | Program definitions (inputs, outputs, commands, invariants) |
| `workflows.yaml` | Workflow state machines (phases, transitions, conditions) |
| `metrics.yaml` | Quality metrics and thresholds |
| `file_categories.yaml` | File categorization rules (extensions, patterns) |

### Key Configuration Features

- **invariants** - Auto-fix rules (e.g., auto-fill resolution from context)
- **input_priorities** - File category preferences per program
- **user_advice_keywords** - Keywords to match user requests to programs
- **priority_when** - Conditional program priority (e.g., strong anomalous → autosol)

## Directory Structure

```
improved_agent_v2/
├── agent/
│   ├── graph.py               # LangGraph definition
│   ├── graph_nodes.py         # Graph node implementations
│   ├── graph_state.py         # State schema
│   ├── workflow_engine.py     # YAML workflow interpreter
│   ├── workflow_state.py      # File categorization, state detection
│   ├── rules_selector.py      # Deterministic program selection
│   ├── program_registry.py    # Program info access
│   ├── template_builder.py    # Command construction with invariants
│   ├── metric_evaluator.py    # YAML-based quality assessment
│   ├── metrics_analyzer.py    # Metrics trend analysis
│   ├── dry_run_manager.py     # Test scenario management
│   └── session.py             # Session management
├── knowledge/
│   ├── programs.yaml          # Program definitions
│   ├── workflows.yaml         # Workflow state machines
│   ├── metrics.yaml           # Metric thresholds
│   ├── file_categories.yaml   # File categorization rules
│   ├── yaml_loader.py         # YAML access functions
│   └── prompts_hybrid.py      # LLM prompt templates
├── tests/
│   ├── test_yaml_config.py    # YAML configuration tests
│   ├── test_workflow_state.py # Workflow state tests
│   ├── test_integration.py    # Integration tests
│   ├── test_metrics_analyzer.py
│   ├── test_dry_run.py
│   └── scenarios/             # Dry run test scenarios
├── yaml_tools.py              # YAML validation & inspection utility
├── log_parsers.py             # Program output parsing
├── run_all_tests.py           # Master test runner
├── README.md                  # This file
└── AGENT_LOGIC.md             # Detailed decision-making logic
```

## YAML Tools

The `yaml_tools.py` utility helps manage YAML configuration:

```bash
# List all YAML files
python yaml_tools.py list

# Validate all YAML files
python yaml_tools.py validate

# Display formatted contents
python yaml_tools.py display knowledge/programs.yaml
python yaml_tools.py display knowledge/file_categories.yaml

# Show all defined terms
python yaml_tools.py terms                  # Simple list
python yaml_tools.py terms --detail normal  # With descriptions
python yaml_tools.py terms --detail full    # Full cross-reference

# Show summary of all configuration
python yaml_tools.py summary

# Compare two configurations
python yaml_tools.py compare old/ new/
```

## Running Tests

```bash
# Run all tests
python run_all_tests.py

# Run individual test suites
python tests/test_yaml_config.py
python tests/test_workflow_state.py
python tests/test_integration.py
python tests/test_dry_run.py
```

## X-ray Workflow

```
xray_initial → xtriage → xray_analyzed
                              │
              ┌───────────────┼───────────────┐
              ↓               ↓               ↓
       predict_and_build    phaser         autosol
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

### Program Selection by State

| State | Valid Programs |
|-------|----------------|
| `xray_initial` | xtriage |
| `xray_analyzed` | predict_and_build*, phaser, autosol† |
| `xray_has_prediction` | process_predicted_model |
| `xray_model_processed` | phaser |
| `xray_has_phases` | autobuild |
| `xray_has_model` | refine |
| `xray_refined` | refine, molprobity, autobuild‡, ligandfit, STOP |

\* Default preferred when sequence available
† Prioritized when strong anomalous signal detected (measurability > 0.10)
‡ Available when R-free above autobuild threshold

## Cryo-EM Workflow

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

## Stop Conditions

The workflow stops when:

1. **Success**: Quality metrics reach targets AND validation completed
2. **Plateau**: < 0.5% improvement for 2+ consecutive cycles
3. **Excessive**: Too many refinement cycles (default: 5)
4. **Validation Gate**: Must run molprobity before stopping if metrics are good

### X-ray Thresholds

| Resolution | Good R-free | Acceptable R-free |
|------------|-------------|-------------------|
| < 1.5 Å | 0.20 | 0.25 |
| 1.5-2.5 Å | 0.25 | 0.30 |
| 2.5-3.5 Å | 0.30 | 0.35 |
| > 3.5 Å | 0.35 | 0.40 |

### Cryo-EM Thresholds

| Quality | Map-Model CC |
|---------|--------------|
| Good | > 0.80 |
| Acceptable | > 0.70 |

## Anomalous Data Detection

The agent automatically detects anomalous signal from xtriage output:

- **Metric**: `anomalous_measurability` (0-1 scale)
- **Strong signal**: measurability > 0.10 → prioritizes `phenix.autosol`
- **Moderate signal**: measurability 0.05-0.10 → SAD may be possible
- **Weak signal**: measurability < 0.05 → SAD unlikely

## Operating Modes

### Standard Mode (LLM + YAML)
```bash
phenix.ai_agent original_files="data.mtz sequence.fa"
```

### Rules-Only Mode (No LLM)
```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz sequence.fa"
```

### Dry Run Mode (Testing)
```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic use_rules_only=True
```

## Modifying Agent Behavior

### Add a New Program

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

  # Optional: auto-fix rules
  invariants:
    - name: requires_something
      check:
        has_strategy: [some_param]
      fix:
        auto_fill_resolution: true

  # Optional: file selection preferences
  input_priorities:
    model:
      categories: [refined, phaser_output]
      exclude_categories: [predicted]

  # Optional: user keywords that suggest this program
  user_advice_keywords:
    - "my analysis"
    - "special processing"
```

### Change Workflow Logic

Edit `knowledge/workflows.yaml`:

```yaml
xray:
  phases:
    my_new_phase:
      programs:
        - program: phenix.my_program
          conditions:
            - has: sequence
          priority_when: strong_anomalous  # Optional
      transitions:
        on_complete: refine
```

### Adjust Thresholds

Edit `knowledge/metrics.yaml`:

```yaml
r_free:
  thresholds:
    good: 0.22  # More stringent
    acceptable: 0.28
```

### Add File Category

Edit `knowledge/file_categories.yaml`:

```yaml
my_output:
  description: "Output from my program"
  subcategory_of: pdb
  patterns:
    - "*my_program*"
    - "*my_output*"
```

### Validate Changes

```bash
python yaml_tools.py validate
python yaml_tools.py terms --detail full
```

## Documentation

| Document | Description |
|----------|-------------|
| `README.md` | This file - overview and quick start |
| `AGENT_LOGIC.md` | Detailed decision-making logic and architecture |
