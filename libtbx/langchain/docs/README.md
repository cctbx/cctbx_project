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
| [reference/VALIDATION.md](reference/VALIDATION.md) | Validation system documentation |
| [guides/ADDING_PROGRAMS.md](guides/ADDING_PROGRAMS.md) | How to add new PHENIX programs |
| [guides/TESTING.md](guides/TESTING.md) | Testing guide and conventions |
| [SAFETY_CHECKS.md](SAFETY_CHECKS.md) | Auto-generated list of all safety validations |

### Implementation Details

| Document | Description |
|----------|-------------|
| [project/CHANGELOG.md](project/CHANGELOG.md) | Version history and release notes |
| [project/THOUGHT_EXPERIMENT.md](project/THOUGHT_EXPERIMENT.md) | Example workflow traces (v110, updated through v112) |
| [project/TRANSPARENCY_LOGGING.md](project/TRANSPARENCY_LOGGING.md) | Event system design and implementation |


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
| `knowledge/recoverable_errors.yaml` | Error patterns and recovery strategies |
| `knowledge/transport.yaml` | Sanitization rules for API transport |

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
| `quiet` | Errors, warnings, one-line cycle summaries |
| `normal` | Key decisions, metrics, commands (default) |
| `verbose` | Full details: file selection, LLM traces, debug events |

### 6. Safety Checks

The agent includes 75+ safety checks across multiple layers:

| Category | Count | Examples |
|----------|-------|----------|
| Sanity Checks | 20 | No data for workflow, model not positioned |
| Input Validation | 29 | Command syntax, API schema, pattern validation |
| Workflow Validation | 8 | Wrong phase, invalid transitions |
| Directive Validation | 7 | Invalid program names, conflicting stops |
| File Validation | 4 | Intermediate file filtering, YAML validation |
| Command Building | 3 | Invariant checks, file existence |
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
2. phenix.map_symmetry      → Detect point-group symmetry (if applicable)
3. phenix.resolve_cryo_em   → Density modification (if half-maps provided)
4. phenix.predict_and_build → AlphaFold prediction
5. phenix.process_predicted_model → Prepare model for docking
6. phenix.dock_in_map       → Dock model into density-modified map
7. phenix.real_space_refine → Iterative refinement (3-6 cycles)
8. phenix.ligandfit         → Fit ligand (if provided)
9. phenix.molprobity        → Geometry validation
10. STOP                    → Target achieved or plateau detected
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
# Run all tests (745+ tests across 34 suites)
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

### Configuration Management

```bash
# List all YAML configuration files
python3 agent/yaml_tools.py list

# Validate all YAML files for syntax and structural errors
python3 agent/yaml_tools.py validate

# Validate a specific file
python3 agent/yaml_tools.py validate programs.yaml

# Display formatted contents of a YAML file
python3 agent/yaml_tools.py display programs

# Compare two YAML files or directories
python3 agent/yaml_tools.py compare programs.yaml programs_backup.yaml

# Show overview of all configuration
python3 agent/yaml_tools.py summary

# Show all defined terms in the configuration system
python3 agent/yaml_tools.py terms
python3 agent/yaml_tools.py terms --detail full  # With cross-references
```

### Session Management

```bash
# Show current session summary
python3 agent/session_tools.py --show

# Show detailed session info (files, metrics, reasoning for each cycle)
python3 agent/session_tools.py --show --detailed

# Remove last N cycles from session
python3 agent/session_tools.py --remove-last 2

# Reset entire session
python3 agent/session_tools.py --reset

# Dry-run (show what would be done without saving)
python3 agent/session_tools.py --remove-last 3 --dry-run

# Use a specific session directory
python3 agent/session_tools.py --dir /path/to/session --show
```

### Documentation Generation

```bash
# Generate safety checks documentation
python3 agent/generate_safety_docs.py > docs/SAFETY_CHECKS.md
```

### RAG Documentation Database

```bash
# Build vector database from PHENIX documentation
python3 run_build_db.py

# Inspect database contents
python3 run_inspect_db.py

# Query the documentation
python3 run_query_docs.py "How do I set up SAD phasing in phenix.autosol?"
```

### Program Validation

```bash
# Validate a specific program's configuration completeness
python3 agent/program_validator.py phenix.polder

# Validate all programs
python3 agent/program_validator.py --all

# List all configured programs
python3 agent/program_validator.py --list
```

### Pattern Management

```bash
# Validate all metric extraction patterns and run tests
python3 agent/pattern_manager.py
```

### Directive Validation

```bash
# Run directive validator self-test (checks program availability detection)
python3 agent/directive_validator.py
```

### Testing

```bash
# Run all standalone tests (no PHENIX required)
python3 tests/run_all_tests.py --quick

# Run all tests (including PHENIX-dependent)
python3 tests/run_all_tests.py

# Verbose output
python3 tests/run_all_tests.py --verbose

# Individual test file
python3 tests/test_file_utils.py

# Tests matching a pattern
python3 tests/run_all_tests.py --pattern "directive"
```

---

## Directory Structure

```
improved_agent_v2/
├── agent/                      # Core agent logic
│   ├── graph.py                # LangGraph state machine definition
│   ├── graph_state.py          # Agent state type definitions
│   ├── graph_nodes.py          # LangGraph node implementations
│   ├── planner.py              # Agent planning and next-move generation
│   ├── workflow_engine.py      # YAML workflow interpreter
│   ├── workflow_state.py       # State detection from files and history
│   ├── command_builder.py      # Unified command generation
│   ├── template_builder.py     # YAML-driven command templates
│   ├── file_utils.py           # Shared file classification (MTZ, etc.)
│   ├── best_files_tracker.py   # Track best files per type
│   ├── error_analyzer.py       # Automatic error recovery
│   ├── advice_preprocessor.py  # README discovery, advice processing
│   ├── directive_extractor.py  # Parse user directives from advice
│   ├── directive_validator.py  # Pre-validate user requests
│   ├── sanity_checker.py       # Red flag detection
│   ├── rules_selector.py       # Rules-only program selection
│   ├── config_loader.py        # Load decision_config.json thresholds
│   ├── memory.py               # Persistent learned syntax tips
│   ├── metrics_analyzer.py     # Metric trends and convergence
│   ├── metric_evaluator.py     # Metric quality evaluation
│   ├── session.py              # Persistent session tracking (AgentSession)
│   ├── api_client.py           # V2 API request/response building
│   ├── phenix_utils.py         # REST encoding, standalone PHENIX utilities
│   ├── transport.py            # Sanitization and encoding
│   ├── rate_limit_handler.py   # LLM rate limiting with backoff
│   ├── program_registry.py     # YAML program registry
│   ├── pattern_manager.py      # Regex pattern management
│   ├── event_log.py            # Structured event logging
│   ├── event_formatter.py      # Output formatting (verbosity levels)
│   ├── dry_run_manager.py      # Testing: dry-run workflow simulation
│   ├── utils.py                # General utility functions
│   ├── command_templates.json  # Program command templates and file slots
│   ├── decision_config.json    # Tiered decision rules and thresholds
│   ├── parameter_fixes.json    # Wrong→correct parameter name mappings
│   ├── yaml_tools.py           # CLI: YAML validation and inspection
│   ├── session_tools.py        # CLI: Session management
│   ├── docs_tools.py           # CLI: Documentation generation
│   ├── generate_safety_docs.py # CLI: Safety checks documentation
│   ├── generate_logic_doc.py   # CLI: Decision logic documentation
│   └── program_validator.py    # CLI: Program config validation
├── knowledge/                  # Configuration & domain knowledge
│   ├── programs.yaml           # Program definitions
│   ├── workflows.yaml          # Workflow state machines
│   ├── metrics.yaml            # Quality thresholds
│   ├── file_categories.yaml    # File categorization
│   ├── patterns.yaml           # Regex patterns
│   ├── recoverable_errors.yaml # Error recovery patterns
│   ├── transport.yaml          # Sanitization rules
│   ├── api_schema.py           # V2 API schema definitions
│   ├── yaml_loader.py          # Configuration loading
│   ├── metric_patterns.py      # YAML-driven metric extraction
│   ├── phenix_programs.py      # Program discovery and introspection
│   ├── program_registration.py # Program detection from logs
│   ├── prompts.py              # LLM prompts for planning and commands
│   ├── prompts_hybrid.py       # Hybrid planning prompts (rules + LLM)
│   └── summary_display.py      # Quality table formatting
├── phenix_ai/                  # Runtime entry points
│   ├── local_agent.py          # Local execution (same process)
│   ├── remote_agent.py         # Server execution (REST API)
│   ├── run_ai_agent.py         # Decision engine (graph execution)
│   ├── run_ai_analysis.py      # Log analysis (standalone)
│   ├── log_parsers.py          # Log metric extraction
│   └── utilities.py            # Shared runtime utilities
├── programs/                   # PHENIX program integration
│   ├── ai_agent.py             # Main PHENIX entry point
│   └── ai_analysis.py          # Log analysis entry point
├── analysis/                   # Log analysis and post-run analysis
│   ├── analyzer.py             # RAG-based log analysis
│   ├── log_info.py             # High-level log info extraction
│   ├── state_extractor.py      # Project state extraction from logs
│   ├── summarizer.py           # Map-reduce log summarization
│   └── agent_session_analyzer.py # Session performance analysis
├── core/                       # LLM integration
│   ├── llm.py                  # Provider abstraction (Google, OpenAI, etc.)
│   └── types.py                # Core data types (AgentPlan, etc.)
├── commands/                   # Command building framework
│   └── base.py                 # Abstract base class for command builders
├── strategies/                 # Planning strategy framework
│   └── base.py                 # Abstract base class for planning strategies
├── validation/                 # Command validation framework
│   ├── base.py                 # Abstract base class for validators
│   ├── core_validator.py       # Syntax validation and LLM-based fixing
│   ├── mtz_utils.py            # MTZ R-free flag checking
│   ├── phenix_refine.py        # phenix.refine-specific validation
│   └── registry.py             # Validator auto-discovery and lookup
├── rag/                        # RAG (Retrieval-Augmented Generation)
│   ├── document_loader.py      # Document loading and chunking
│   ├── retriever.py            # Vector DB retrieval with reranking
│   └── vector_store.py         # Chroma vector store creation
├── utils/                      # General utilities
│   ├── query.py                # Documentation query interface
│   ├── run_utils.py            # Log parsing and HTML output
│   └── text_processing.py      # Text block extraction helpers
├── phenix_knowledge.py         # Phenix program allow-list and syntax hints
├── phenix_learned_info/        # Persistent learned knowledge
│   └── phenix_learned_memory.json # Learned syntax tips per program
├── run_build_db.py             # CLI: Build RAG documentation database
├── run_inspect_db.py           # CLI: Inspect RAG database contents
├── run_query_docs.py           # CLI: Query PHENIX documentation
├── tst_langchain_tools.py      # Unit tests for external modules
├── tests/                      # Test suites (~745+ tests)
│   ├── run_all_tests.py        # Test runner (32 registered suites)
│   ├── test_utils.py           # Assert helpers (cctbx-style)
│   └── scenarios/              # Dry-run test scenarios
└── docs/                       # Documentation
    ├── README.md               # This file
    ├── OVERVIEW.md             # Technical overview
    ├── SAFETY_CHECKS.md        # Auto-generated safety checks
    ├── guides/                 # How-to guides
    │   ├── USER_DIRECTIVES.md  # Natural language guidance
    │   ├── ADDING_PROGRAMS.md  # Adding new PHENIX programs
    │   └── TESTING.md          # Testing guide
    ├── reference/              # API and logic reference
    │   ├── ARCHITECTURE.md     # Component deep-dive
    │   ├── API_REFERENCE.md    # V2 JSON API spec
    │   └── VALIDATION.md       # Validation system
    ├── project/                # Design history and changelog
    │   ├── CHANGELOG.md        # Version history
    │   ├── THOUGHT_EXPERIMENT.md # Example workflow traces
    │   └── TRANSPARENCY_LOGGING.md # Event system design
```

---

## See Also

- [PHENIX Documentation](https://phenix-online.org/documentation/)
- [AlphaFold](https://alphafold.ebi.ac.uk/) - Structure prediction
