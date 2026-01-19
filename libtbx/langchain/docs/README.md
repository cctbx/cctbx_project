# PHENIX AI Agent - Documentation

## Documents

| Document | Description |
|----------|-------------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | System architecture, components, data flow, CommandBuilder pipeline |
| [API_DOCUMENTATION.md](API_DOCUMENTATION.md) | V2 API request/response schemas, transport encoding |
| [VALIDATION.md](VALIDATION.md) | Sanity checking, red flags, validation gates |

## Quick Reference

### Key Components

| Component | File | Purpose |
|-----------|------|---------|
| CommandBuilder | `agent/command_builder.py` | Unified command generation pipeline |
| Transport | `agent/transport.py` | Request/response encoding for client-server |
| WorkflowState | `agent/workflow_state.py` | File categorization, state detection |
| SanityChecker | `agent/sanity_checker.py` | Validation and red flag detection |
| ProgramRegistry | `agent/program_registry.py` | YAML program definitions access |

### Configuration Files

| File | Purpose |
|------|---------|
| `knowledge/programs.yaml` | Program definitions, inputs, outputs, invariants |
| `knowledge/workflows.yaml` | Workflow state machines, transitions |
| `knowledge/metrics.yaml` | Quality metric thresholds |
| `knowledge/file_categories.yaml` | File categorization rules |
| `knowledge/transport.yaml` | Transport encoding settings |

### Testing

```bash
# Run all tests (12 test suites)
python tests/run_all_tests.py

# Run with verbose output
python tests/run_all_tests.py --verbose

# Skip slow integration tests
python tests/run_all_tests.py --quick
```

## See Also

- Main [AGENT_SUMMARY.md](../AGENT_SUMMARY.md) for overview and workflows
- [YAML_SCORING_PLAN.md](YAML_SCORING_PLAN.md) for planned YAML-based scoring (future)
