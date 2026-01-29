# PHENIX AI Agent Documentation

## Overview

The PHENIX AI Agent automates macromolecular structure determination for X-ray crystallography and cryo-EM. It combines LLM-based decision making with rules-based validation to intelligently guide workflows from raw data to refined structures.

## Quick Start

```bash
# X-ray structure determination
phenix.ai_agent original_files="data.mtz sequence.fa"

# Cryo-EM structure determination  
phenix.ai_agent original_files="map.mrc sequence.fa"

# With user guidance
phenix.ai_agent original_files="data.mtz model.pdb ligand.pdb" \
    project_advice="Solve the structure and fit the ligand"

# Stop after specific step
phenix.ai_agent original_files="data.mtz seq.fa" \
    project_advice="Stop after one refinement job"

# Stepwise mode - more control with intermediate checkpoints
# Forces predict_and_build to stop after prediction (both X-ray and cryo-EM)
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"

# Rules-only mode (deterministic, no LLM)
phenix.ai_agent use_rules_only=True original_files="data.mtz sequence.fa"

# Control output verbosity
phenix.ai_agent verbosity=verbose original_files="data.mtz sequence.fa"
```

## Documentation Index

### For Users

| Document | Description |
|----------|-------------|
| [guides/USER_DIRECTIVES.md](guides/USER_DIRECTIVES.md) | How to guide the agent with natural language |

### For Developers

| Document | Description |
|----------|-------------|
| [OVERVIEW.md](OVERVIEW.md) | Technical overview, architecture, data flow |
| [reference/ARCHITECTURE.md](reference/ARCHITECTURE.md) | Deep-dive into components and transport |
| [reference/API_REFERENCE.md](reference/API_REFERENCE.md) | V2 JSON API specification |
| [guides/ADDING_PROGRAMS.md](guides/ADDING_PROGRAMS.md) | How to add new PHENIX programs |
| [guides/TESTING.md](guides/TESTING.md) | Testing guide and conventions |
| [SAFETY_CHECKS.md](SAFETY_CHECKS.md) | Auto-generated list of all safety validations |

### Implementation Details

| Document | Description |
|----------|-------------|
| [project/TRANSPARENCY_LOGGING.md](project/TRANSPARENCY_LOGGING.md) | Event system design and implementation |
| [project/THOUGHT_EXPERIMENT.md](project/THOUGHT_EXPERIMENT.md) | Example workflow traces |
| [project/CHANGELOG.md](project/CHANGELOG.md) | Version history |
| [implementation/PROGRAM_CONFIG_ROBUSTNESS.md](implementation/PROGRAM_CONFIG_ROBUSTNESS.md) | Plan for robust program configuration |

---

## Key Features

### 1. YAML-Driven Configuration

All domain knowledge is externalized to editable YAML files:

| File | Purpose |
|------|---------|
| `knowledge/programs.yaml` | Program definitions: inputs, outputs, log parsing |
| `knowledge/workflows.yaml` | Workflow phases, transitions, stop conditions |
| `knowledge/metrics.yaml` | Quality thresholds, display formats |
| `knowledge/file_categories.yaml` | File type categorization rules |
| `knowledge/patterns.yaml` | Regex patterns for metric extraction |

### 2. Intelligent Workflow Management

The agent automatically:
- **Detects experiment type** from input files (MTZ → X-ray, MRC → cryo-EM)
- **Tracks workflow state** and enforces valid program sequences
- **Monitors quality metrics** (R-free, map CC) and detects convergence
- **Stops automatically** when targets reached or improvement plateaus

### 3. Automation Modes

Control workflow granularity with `maximum_automation`:

| Mode | Setting | Behavior |
|------|---------|----------|
| **Automated** | `maximum_automation=True` (default) | Full automation - `predict_and_build` runs the complete workflow |
| **Stepwise** | `maximum_automation=False` | More control - `predict_and_build` stops after prediction |

**Automated workflow (X-ray):**
```
xtriage → predict_and_build(full) → xray_refined → refine cycles → STOP
```

**Stepwise workflow (X-ray):**
```
xtriage → predict_and_build(stop_after_predict) → process_predicted_model → phaser → refine → ...
```

Stepwise mode gives you checkpoints to inspect intermediate results before proceeding.

### 4. User Directives

Natural language guidance for the workflow:

```bash
# Stop conditions
project_advice="Stop after one refinement job"
project_advice="Stop when R-free is below 0.25"

# Workflow preferences
project_advice="Skip autobuild and go straight to refinement"
project_advice="Use SAD phasing instead of molecular replacement"

# Task focus
project_advice="Focus on ligand fitting"
```

### 5. Transparent Decision Making

The agent explains every decision at configurable verbosity:

| Level | Shows |
|-------|-------|
| `quiet` | Errors only, one-line cycle summaries |
| `normal` | Key decisions, metrics, commands (default) |
| `verbose` | Full details: file selection, LLM traces, debug info |

Note: `debug` is an alias for `verbose` (3 levels total).

### 6. Safety Checks

The agent includes 70+ safety checks across multiple layers:

| Category | Count | Examples |
|----------|-------|----------|
| Sanity Checks | 20 | No data for workflow, model not positioned |
| Directive Validation | 7 | Invalid program names, conflicting stops |
| Workflow Validation | 8 | Wrong phase, invalid transitions |
| Post-Processing | 4 | Ligand workflow conflict resolution |

See [SAFETY_CHECKS.md](SAFETY_CHECKS.md) for the complete auto-generated list.

### 7. User Feedback

When the agent can't fulfill a request, it explains why:

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

---

## Typical Workflows

### X-ray Crystallography (8-12 cycles)

```
1. phenix.xtriage           → Analyze data quality, detect twinning
2. phenix.predict_and_build → AlphaFold prediction + model building
3. phenix.refine            → Iterative refinement (3-6 cycles)
4. phenix.ligandfit         → Fit ligand (if provided, when R-free < 0.35)
5. phenix.pdbtools          → Combine protein + ligand
6. phenix.refine            → Final refinement with ligand
7. phenix.molprobity        → Geometry validation
8. STOP                     → Target achieved or plateau detected
```

### Cryo-EM (6-10 cycles)

```
1. phenix.mtriage           → Analyze map, determine resolution
2. phenix.resolve_cryo_em   → Density modification (if half-maps provided)
3. phenix.predict_and_build → AlphaFold prediction
4. phenix.process_predicted_model → Prepare model for docking
5. phenix.dock_in_map       → Dock model into density-modified map
6. phenix.real_space_refine → Iterative refinement (3-6 cycles)
7. phenix.molprobity        → Geometry validation
8. STOP                     → Target achieved or plateau detected
```

---

## Execution Modes

| Mode | Command | Use Case |
|------|---------|----------|
| **Standard** | `phenix.ai_agent ...` | Production use |
| **Rules-only** | `use_rules_only=True` | Deterministic, no LLM |
| **Dry-run** | `dry_run=True dry_run_scenario=xray_basic` | Testing workflows |
| **Server** | via REST API | Distributed systems |

---

## Testing

The test suite uses **cctbx-style testing** with plain functions and fail-fast behavior.

```bash
# Run all tests (300+ tests across 12+ files)
cd improved_agent_v2v1
python tests/run_all_tests.py

# Quick mode (standalone tests only, no PHENIX required)
python tests/run_all_tests.py --quick

# Individual test file
python tests/test_directive_extractor.py

# Tests matching a pattern
python tests/run_all_tests.py --pattern "directive"
```

Key features:
- **Fail-fast**: First assertion failure stops with full traceback
- **Plain functions**: No `unittest.TestCase` classes needed
- **Assert helpers**: `assert_equal()`, `assert_in()`, etc. from `test_utils.py`

See [guides/TESTING.md](guides/TESTING.md) for full documentation.

---

## Command-Line Tools

```bash
# Validate YAML configuration
python3 agent/yaml_tools.py validate

# Display program definitions
python3 agent/yaml_tools.py display programs

# Manage sessions
python3 agent/session_tools.py --show
python3 agent/session_tools.py --remove-last 2

# Generate reference documentation
python3 agent/docs_tools.py > docs/reference/AGENT_LOGIC.md
```

---

## Directory Structure

```
improved_agent_v2/
├── agent/                    # Core agent logic
│   ├── graph_nodes.py        # LangGraph node implementations
│   ├── workflow_engine.py    # YAML workflow interpreter
│   ├── workflow_state.py     # State detection
│   ├── command_builder.py    # Command generation
│   ├── file_categorization.py# Semantic file classification
│   ├── best_files_tracker.py # Track best files per type
│   ├── event_log.py          # Structured event logging
│   ├── event_formatter.py    # Output formatting
│   ├── sanity_checker.py     # Red flag detection
│   ├── rules_selector.py     # Rules-only program selection
│   ├── directive_extractor.py# Parse user directives
│   └── ...
├── knowledge/                # Configuration & domain knowledge
│   ├── programs.yaml         # Program definitions
│   ├── workflows.yaml        # Workflow state machines
│   ├── metrics.yaml          # Quality thresholds
│   ├── file_categories.yaml  # File categorization
│   ├── patterns.yaml         # Regex patterns
│   └── transport.yaml        # Sanitization rules
├── phenix_ai/                # Runtime entry points
│   ├── local_agent.py        # Local execution
│   ├── remote_agent.py       # Server execution
│   └── run_ai_agent.py       # Decision engine
├── programs/                 # PHENIX program integration
│   └── ai_agent.py           # Main entry point
├── tests/                    # Test suites (~700 tests total)
│   ├── run_all_tests.py      # Test runner
│   └── scenarios/            # Dry-run test scenarios
└── docs/                     # This documentation
```

---

## See Also

- [PHENIX Documentation](https://phenix-online.org/documentation/)
- [AlphaFold](https://alphafold.ebi.ac.uk/) - Structure prediction
