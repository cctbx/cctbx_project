# Validation and Sanity Checking

## Overview

The PHENIX AI Agent includes a sanity-checking system that detects when the agent
is in an impossible or nonsensical state. When critical issues are detected, the
agent aborts gracefully with a clear explanation and suggestions for resolution.

## Design Principles

1. **Abort only when truly necessary** - Red flags indicate real problems, not edge cases
2. **User control** - Configurable abort behavior
3. **Clear, actionable messages** - Explain what went wrong and how to fix it
4. **Resumable sessions** - After fixing issues, users can continue the session

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `abort_on_red_flags` | `True` | Abort on critical issues (experiment type change, impossible states) |
| `abort_on_warnings` | `False` | Also abort on warnings (metric anomalies, missing resolution) |

## Red Flag Categories

### Critical Issues (abort_on_red_flags)

| Code | Check | Trigger |
|------|-------|---------|
| `experiment_type_changed` | Experiment type stability | Type differs from initial (xray↔cryoem) |
| `no_model_for_refine` | Model required for refinement | Refine state reached with no PDB files |
| `no_data_for_workflow` | Data required for workflow | No MTZ (xray) or map (cryoem) |
| `repeated_failures` | Same error repeated | Same program fails 3+ times with identical error |

### Warnings (abort_on_warnings)

| Code | Check | Trigger |
|------|-------|---------|
| `resolution_unknown` | Resolution should be established | Entering refine phase without resolution |
| `r_free_spike` | R-free anomaly | R-free increases by >0.15 in one cycle |
| `map_cc_drop` | Map CC anomaly | Map CC drops by >0.3 in one cycle |

## Validation Gates

Before the workflow can stop with success, validation must be completed:

1. **Quality Check Required**: Must run `phenix.molprobity` before stopping
2. **Metrics Must Be Acceptable**: R-free or map CC must meet thresholds
3. **No Red Flags**: No critical issues detected

### Validation Gate Logic

```
if workflow_should_stop:
    if not validation_completed:
        return "Run molprobity first"
    if red_flags:
        return "Cannot stop: " + red_flag_message
    return "STOP"
```

## Implementation

### SanityChecker Class

Located in `agent/sanity_checker.py`:

```python
from agent.sanity_checker import SanityChecker

checker = SanityChecker()

# Run all checks
result = checker.check_state(
    workflow_state="xray_refined",
    experiment_type="xray",
    available_files=[...],
    history=[...],
    session_info={...}
)

if result.red_flags:
    print("Critical issues:", result.red_flags)
if result.warnings:
    print("Warnings:", result.warnings)
```

### YAML Configuration

Sanity checks are defined in `knowledge/workflows.yaml`:

```yaml
sanity_checks:
  critical:
    - name: experiment_type_stable
      message: "Experiment type changed from {initial} to {current}"
      suggestion: "Start a new session for different experiment types"

    - name: model_required_for_refine
      message: "No model file found but refinement is required"
      suggestion: "Check that model building succeeded"

  warnings:
    - name: resolution_established
      message: "Resolution not established before refinement"
      suggestion: "Run xtriage or mtriage first"
```

## Response Format

When validation issues are detected, they appear in the response:

```json
{
  "api_version": "2.0",
  "stop": true,
  "stop_reason": "red_flag",
  "metadata": {
    "red_flags": [
      {
        "code": "experiment_type_changed",
        "message": "Experiment type changed from xray to cryoem",
        "suggestion": "Start a new session for different experiment types"
      }
    ],
    "warnings": []
  }
}
```

## Testing

Run sanity checker tests:

```bash
python tests/test_sanity_checker.py
```

Test scenarios include:
- Experiment type change detection
- Missing model detection
- Repeated failure detection
- Metric spike detection
- Validation gate enforcement
