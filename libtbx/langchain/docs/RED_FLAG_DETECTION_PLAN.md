# Red Flag Detection System - Implementation Plan

## Overview

A sanity-checking system that detects when the PHENIX AI Agent is in an impossible or nonsensical state. When critical issues are detected, the agent aborts gracefully with a clear explanation and suggestions for resolution.

## Design Principles

1. **Abort only when truly necessary** - Red flags indicate real problems, not edge cases
2. **User control** - `abort_on_red_flags=True` (default), `abort_on_warnings=False` (default)
3. **Clear, actionable messages** - Explain what went wrong and how to fix it
4. **Resumable sessions** - After fixing issues, users can continue the session
5. **Clean integration** - Fits naturally into existing YAML-first architecture

---

## User-Facing Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `abort_on_red_flags` | `True` | Abort on critical issues (experiment type change, impossible states) |
| `abort_on_warnings` | `False` | Also abort on warnings (metric anomalies, missing resolution) |

---

## Red Flag Categories

### Critical (abort_on_red_flags)

| Code | Check | Trigger |
|------|-------|---------|
| `experiment_type_changed` | Experiment type stability | Type differs from initial (xray↔cryoem) |
| `no_model_for_refine` | Model required for refinement | Refine state reached with no PDB files |
| `no_data_for_workflow` | Data required for workflow | No MTZ (xray) or map (cryoem) |
| `repeated_failures` | Same error repeated | Same program fails 3+ times with identical error |

### Warning (abort_on_warnings)

| Code | Check | Trigger |
|------|-------|---------|
| `resolution_unknown` | Resolution should be established | Entering refine phase without resolution |
| `multi_sequence_stepwise` | Multi-chain + stepwise mode | Multiple FASTA chains with automated mode disabled |
| `r_free_spike` | R-free anomaly | R-free increases by >0.15 in one cycle |
| `map_cc_drop` | Map CC anomaly | Map CC drops by >0.3 in one cycle |

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     workflows.yaml                          │
│                                                             │
│  sanity_checks:                                             │
│    critical:                                                │
│      - name: experiment_type_stable                         │
│        message: "..."                                       │
│        suggestion: "..."                                    │
│    warning:                                                 │
│      - name: resolution_established                         │
│        ...                                                  │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│              agent/sanity_checker.py (NEW)                  │
│                                                             │
│  class SanityChecker:                                       │
│      def check(context, session) -> SanityResult            │
│          - Loads checks from workflows.yaml                 │
│          - Evaluates each check against context             │
│          - Returns issues with severity/message/suggestion  │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                 agent/graph_nodes.py                        │
│                                                             │
│  perceive():                                                │
│      ...existing logic...                                   │
│      sanity = checker.check(context, session_info)         │
│      if sanity.should_abort:                                │
│          return abort_state(sanity.issues)                  │
│      elif sanity.issues:                                    │
│          log warnings                                       │
│      ...continue...                                         │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                   agent/session.py                          │
│                                                             │
│  Session data now includes:                                 │
│      - experiment_type: str (locked after cycle 1)          │
│      - abort_on_red_flags: bool                             │
│      - abort_on_warnings: bool                              │
└─────────────────────────────────────────────────────────────┘
```

---

## File Changes Summary

| File | Change Type | Description |
|------|-------------|-------------|
| `knowledge/workflows.yaml` | Add section | New `sanity_checks` section with check definitions |
| `agent/sanity_checker.py` | **New file** | SanityChecker class and SanityResult dataclass |
| `agent/session.py` | Extend | Add `experiment_type`, flags to session data |
| `agent/graph_nodes.py` | Modify | Call sanity checker in PERCEIVE, handle abort |
| `agent/workflow_state.py` | Minor | Pass session info through for sanity checking |
| `ai_agent.py` | Extend | Accept and pass through new parameters |
| `knowledge/yaml_loader.py` | Extend | Add `load_sanity_checks()` function |

---

## Implementation Phases

### Phase 1: Infrastructure
**Goal**: Framework in place, no behavior change

#### Step 1.1: Extend session.py
Add to session data structure:
```python
self.data = {
    ...existing fields...
    "experiment_type": None,  # "xray" or "cryoem", locked after cycle 1
    "abort_on_red_flags": True,
    "abort_on_warnings": False,
}
```

Add methods:
```python
def set_experiment_type(self, exp_type):
    """Set experiment type (only if not already set)."""
    if self.data.get("experiment_type") is None:
        self.data["experiment_type"] = exp_type
        self.save()
        return True
    return False

def get_experiment_type(self):
    """Get locked experiment type."""
    return self.data.get("experiment_type")
```

#### Step 1.2: Create sanity_checker.py
```python
"""
Sanity Checker for PHENIX AI Agent.

Detects impossible or nonsensical states and provides
clear error messages with suggestions for resolution.
"""

from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class SanityIssue:
    """A single sanity check failure."""
    severity: str          # "critical" or "warning"
    code: str              # e.g., "experiment_type_changed"
    message: str           # Human-readable description
    suggestion: str        # How to fix it
    details: dict = field(default_factory=dict)  # Additional context

@dataclass
class SanityResult:
    """Result of sanity checking."""
    ok: bool                          # True if no issues
    issues: List[SanityIssue]         # All detected issues
    should_abort: bool                # True if should stop execution
    abort_message: Optional[str]      # Formatted message for user


class SanityChecker:
    """Checks for impossible or nonsensical agent states."""

    def __init__(self):
        self._load_checks()

    def _load_checks(self):
        """Load check definitions from workflows.yaml."""
        # Implementation in Step 1.4
        pass

    def check(self, context, session_info, abort_on_red_flags=True,
              abort_on_warnings=False) -> SanityResult:
        """
        Run all sanity checks against current context.

        Args:
            context: Dict with workflow context (files, history, etc.)
            session_info: Dict with session state (experiment_type, etc.)
            abort_on_red_flags: Whether to abort on critical issues
            abort_on_warnings: Whether to abort on warnings too

        Returns:
            SanityResult with any issues found
        """
        # Implementation in Phase 2
        return SanityResult(ok=True, issues=[], should_abort=False,
                           abort_message=None)
```

#### Step 1.3: Add sanity_checks to workflows.yaml
```yaml
# =============================================================================
# SANITY CHECKS
# =============================================================================
# These checks detect impossible or nonsensical states.
# Critical checks cause abort (if abort_on_red_flags=True)
# Warning checks are logged (and abort if abort_on_warnings=True)

sanity_checks:
  critical:
    # Checks that indicate fundamental problems
    - name: experiment_type_stable
      description: "Experiment type should not change during workflow"
      message: "Experiment type changed from {initial} to {current}"
      suggestion: "This usually indicates spurious files were detected. Check for intermediate files in subdirectories like run_mr/. Remove them and restart, or use session_utils to remove recent cycles."

    - name: model_exists_for_refine
      description: "Cannot refine without a model"
      message: "Refinement requested but no model file available"
      suggestion: "Ensure molecular replacement or model building completed successfully before refinement."

    - name: data_exists_for_workflow
      description: "Cannot proceed without experimental data"
      message: "No {data_type} found for {experiment_type} workflow"
      suggestion: "Provide the required data files: MTZ for X-ray, or MRC/CCP4 map for cryo-EM."

    - name: repeated_failures
      description: "Same program failing repeatedly"
      message: "Program {program} has failed {count} times with: {error}"
      suggestion: "Check the error message and input files. The program may have a fundamental problem with the provided data."

  warning:
    # Checks that indicate potential problems
    - name: resolution_established
      description: "Resolution should be known before refinement"
      message: "Entering refinement phase without established resolution"
      suggestion: "Run xtriage (X-ray) or mtriage (cryo-EM) first to determine data resolution."

    - name: multi_sequence_stepwise
      description: "Multiple sequences may need special handling"
      message: "Multiple sequences detected ({count} chains) with stepwise mode"
      suggestion: "Consider using automated mode (maximum_automation=True) or provide pre-processed multi-chain model."

    - name: r_free_spike
      description: "R-free should not increase dramatically"
      message: "R-free jumped from {old:.3f} to {new:.3f} ({change:+.3f})"
      suggestion: "This may indicate the model was corrupted or wrong data was used. Check the refinement input files."

    - name: map_cc_drop
      description: "Map correlation should not drop dramatically"
      message: "Map CC dropped from {old:.2f} to {new:.2f}"
      suggestion: "The model may have been displaced from the map. Check docking and refinement results."
```

#### Step 1.4: Extend yaml_loader.py
```python
def load_sanity_checks(force_reload=False):
    """
    Load sanity check definitions from workflows.yaml.

    Returns:
        dict: {"critical": [...], "warning": [...]}
    """
    workflows = load_workflows(force_reload)
    return workflows.get("sanity_checks", {"critical": [], "warning": []})
```

#### Step 1.5: Integrate into PERCEIVE node
In `graph_nodes.py`, after workflow_state detection:
```python
# Run sanity checks
from libtbx.langchain.agent.sanity_checker import SanityChecker

checker = SanityChecker()
sanity_context = {
    "experiment_type": workflow_state.get("experiment_type"),
    "initial_experiment_type": state.get("session_info", {}).get("experiment_type"),
    "has_model": bool(workflow_state.get("categorized_files", {}).get("pdb")),
    "state": workflow_state.get("state"),
    "history": history,
    "metrics_history": metrics_history,
    # ... other relevant context
}

sanity_result = checker.check(
    sanity_context,
    session_info=state.get("session_info", {}),
    abort_on_red_flags=state.get("abort_on_red_flags", True),
    abort_on_warnings=state.get("abort_on_warnings", False)
)

# Log any warnings
for issue in sanity_result.issues:
    if issue.severity == "warning":
        state = _log(state, "SANITY WARNING: %s" % issue.message)

# Handle abort
if sanity_result.should_abort:
    state = _log(state, "SANITY ABORT: %s" % sanity_result.abort_message)
    return {
        **state,
        "command": "STOP",
        "stop": True,
        "stop_reason": "red_flag",
        "red_flag_issues": [asdict(i) for i in sanity_result.issues],
        "abort_message": sanity_result.abort_message,
    }
```

#### Step 1.6: Pass parameters through ai_agent.py
Add to params and pass through to graph state:
```python
abort_on_red_flags = getattr(self.params.ai_analysis, 'abort_on_red_flags', True)
abort_on_warnings = getattr(self.params.ai_analysis, 'abort_on_warnings', False)
```

**Phase 1 Deliverable**: Framework in place, all existing tests pass, sanity checker returns ok=True always.

---

### Phase 2: Critical Checks
**Goal**: Implement checks that catch obvious problems

#### Step 2.1: Experiment type stability
```python
def _check_experiment_type_stable(self, context, session_info):
    """Check that experiment type hasn't changed."""
    initial = session_info.get("experiment_type")
    current = context.get("experiment_type")

    if initial and current and initial != current:
        return SanityIssue(
            severity="critical",
            code="experiment_type_changed",
            message=f"Experiment type changed from {initial} to {current}",
            suggestion="This usually indicates spurious files were detected...",
            details={"initial": initial, "current": current}
        )
    return None
```

#### Step 2.2: Model exists for refinement
```python
def _check_model_for_refine(self, context):
    """Check that model exists when entering refine state."""
    state = context.get("state", "")
    has_model = context.get("has_model", False)

    if "refine" in state.lower() and not has_model:
        return SanityIssue(
            severity="critical",
            code="no_model_for_refine",
            message="Refinement requested but no model file available",
            suggestion="Ensure molecular replacement or model building completed...",
        )
    return None
```

#### Step 2.3: Data exists for workflow
```python
def _check_data_exists(self, context):
    """Check that required data exists for experiment type."""
    exp_type = context.get("experiment_type")

    if exp_type == "xray" and not context.get("has_mtz"):
        return SanityIssue(
            severity="critical",
            code="no_data_for_workflow",
            message="No MTZ found for X-ray workflow",
            suggestion="Provide MTZ file with diffraction data",
        )

    if exp_type == "cryoem" and not context.get("has_map"):
        return SanityIssue(
            severity="critical",
            code="no_data_for_workflow",
            message="No map found for cryo-EM workflow",
            suggestion="Provide MRC or CCP4 map file",
        )
    return None
```

#### Step 2.4: Repeated failures
```python
def _check_repeated_failures(self, context):
    """Check for repeated identical failures."""
    history = context.get("history", [])

    # Count recent failures by program
    failure_counts = {}
    for h in history[-10:]:  # Last 10 cycles
        if h.get("result", "").startswith("FAIL"):
            prog = h.get("program", "unknown")
            error = h.get("result", "")[:100]  # First 100 chars
            key = (prog, error)
            failure_counts[key] = failure_counts.get(key, 0) + 1

    for (prog, error), count in failure_counts.items():
        if count >= 3:
            return SanityIssue(
                severity="critical",
                code="repeated_failures",
                message=f"Program {prog} has failed {count} times",
                suggestion="Check the error message and input files...",
                details={"program": prog, "count": count, "error": error}
            )
    return None
```

**Phase 2 Deliverable**: Core safety net catches obvious problems.

---

### Phase 3: Warning Checks
**Goal**: Add checks that warn but don't abort by default

#### Step 3.1: Resolution established
```python
def _check_resolution_established(self, context):
    """Warn if entering refine without resolution."""
    state = context.get("state", "")
    resolution = context.get("resolution")

    if "refine" in state.lower() and resolution is None:
        return SanityIssue(
            severity="warning",
            code="resolution_unknown",
            message="Entering refinement without established resolution",
            suggestion="Run xtriage or mtriage first...",
        )
    return None
```

#### Step 3.2: Multi-sequence handling
```python
def _check_multi_sequence(self, context, session_info):
    """Warn about multiple sequences in stepwise mode."""
    sequence_count = context.get("sequence_chain_count", 1)
    automation = session_info.get("maximum_automation", True)

    if sequence_count > 1 and not automation:
        return SanityIssue(
            severity="warning",
            code="multi_sequence_stepwise",
            message=f"Multiple sequences ({sequence_count} chains) with stepwise mode",
            suggestion="Consider automated mode or pre-processed model",
        )
    return None
```

#### Step 3.3: Metric anomalies
```python
def _check_metric_anomalies(self, context):
    """Warn about dramatic metric changes."""
    metrics = context.get("metrics_history", [])
    if len(metrics) < 2:
        return None

    prev = metrics[-2]
    curr = metrics[-1]
    issues = []

    # R-free spike
    if prev.get("r_free") and curr.get("r_free"):
        change = curr["r_free"] - prev["r_free"]
        if change > 0.15:
            issues.append(SanityIssue(
                severity="warning",
                code="r_free_spike",
                message=f"R-free jumped from {prev['r_free']:.3f} to {curr['r_free']:.3f}",
                suggestion="Model may be corrupted...",
            ))

    # Map CC drop
    if prev.get("map_cc") and curr.get("map_cc"):
        change = prev["map_cc"] - curr["map_cc"]
        if change > 0.3:
            issues.append(SanityIssue(
                severity="warning",
                code="map_cc_drop",
                message=f"Map CC dropped from {prev['map_cc']:.2f} to {curr['map_cc']:.2f}",
                suggestion="Model may be displaced...",
            ))

    return issues if issues else None
```

**Phase 3 Deliverable**: More informative feedback about potential issues.

---

### Phase 4: Message Formatting and Analysis
**Goal**: Make abort messages genuinely helpful

#### Step 4.1: Format abort message
```python
def _format_abort_message(self, issues: List[SanityIssue]) -> str:
    """Format a clear, helpful abort message."""
    lines = [
        "=" * 60,
        "AGENT STOPPED: Sanity check failed",
        "=" * 60,
        "",
    ]

    for i, issue in enumerate(issues, 1):
        lines.append(f"Issue {i}: {issue.code}")
        lines.append(f"  {issue.message}")
        lines.append(f"  Suggestion: {issue.suggestion}")
        if issue.details:
            lines.append(f"  Details: {issue.details}")
        lines.append("")

    lines.append("To resume after fixing the issue:")
    lines.append("  1. Address the problem described above")
    lines.append("  2. Run the agent again with the same session directory")
    lines.append("  3. Or use session_utils --remove-last to undo recent cycles")
    lines.append("")

    return "\n".join(lines)
```

#### Step 4.2: Root cause hints for experiment type change
```python
def _analyze_experiment_type_change(self, context) -> str:
    """Try to identify why experiment type changed."""
    # Look for suspicious files
    files = context.get("all_files", [])
    suspicious = []

    for f in files:
        if "/run_mr/" in f:
            suspicious.append(f"Intermediate file: {f}")
        if f.endswith(".mtz") and "run_mr" in f:
            suspicious.append(f"Spurious MTZ: {f}")

    if suspicious:
        return "Possible cause:\n" + "\n".join(f"  - {s}" for s in suspicious)
    return ""
```

**Phase 4 Deliverable**: Users can understand and fix problems.

---

### Phase 5: Testing
**Goal**: Ensure checks work correctly

#### Test Cases

| Test | Input | Expected |
|------|-------|----------|
| `test_experiment_type_change` | Mock xray→cryoem switch | Critical issue, abort |
| `test_experiment_type_stable` | Consistent xray workflow | No issues |
| `test_no_model_for_refine` | Refine state, no PDB | Critical issue |
| `test_valid_refine` | Refine state with PDB | No issues |
| `test_repeated_failures` | Same error 3x | Critical issue |
| `test_r_free_spike` | R-free +0.20 | Warning |
| `test_normal_r_free_change` | R-free +0.05 | No issues |
| `test_abort_on_warnings_false` | Warning issue | No abort |
| `test_abort_on_warnings_true` | Warning issue | Abort |

**Phase 5 Deliverable**: Robust, tested system.

---

## Session Resumability

Sessions remain resumable after abort:

1. **Session state preserved** - All cycles up to abort are saved
2. **User fixes issue** - Remove spurious files, add missing data, etc.
3. **Restart agent** - With same session directory
4. **Agent continues** - From where it left off (or user removes cycles with session_utils)

The abort is a clean stop, not a crash. The session JSON is valid and complete.

---

## Summary

This system adds a safety net that catches fundamental problems early, provides clear explanations, and allows users to fix issues and continue. It integrates cleanly with the existing YAML-first architecture and doesn't add significant complexity to the core workflow.
