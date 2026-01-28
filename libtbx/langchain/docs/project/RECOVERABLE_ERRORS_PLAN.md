# Recoverable Errors: Implementation Plan

## Overview

This document describes the implementation of automatic error recovery for the PHENIX AI Agent, starting with the "ambiguous data labels" error that occurs when an MTZ file contains multiple suitable data arrays.

## The Problem

When an MTZ file contains multiple data arrays (e.g., both merged intensities and anomalous pairs), PHENIX programs fail with a structured error message:

```
Multiple equally suitable arrays of observed xray data found.

Possible choices:
  /path/to/data.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
  /path/to/data.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-),merged

Please use scaling.input.xray_data.obs_labels
to specify an unambiguous substring of the target label.
```

Currently, this causes the workflow to fail. However, the error is:
1. **Structured** - Contains parseable information
2. **Prescriptive** - Tells us exactly what keyword to use
3. **Solvable** - We can intelligently select the right array based on context

## Design Principles

1. **Extract keyword from error message** - Don't hardcode parameter names; parse them from the error
2. **Use next cycle for retry** - Force retry of the same program with recovery flags
3. **Give up after 3 attempts** - Prevent infinite loops on unrecoverable situations
4. **Key recovery by filename** - Store recovery strategies per-file to avoid cross-contamination
5. **Use existing strategy system** - Inject recovery as strategy parameters, not command string hacking
6. **YAML-driven patterns** - Externalize error signatures for maintainability
7. **Opt-out available** - `auto_recovery=False` for debugging and power users

---

## Architecture

### Component Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        ai_agent.py                               │
│                     _run_single_cycle()                          │
│                            │                                     │
│                    [Program FAILED]                              │
│                            │                                     │
│                            ▼                                     │
│              ┌─────────────────────────┐                        │
│              │     ErrorAnalyzer       │                        │
│              │  analyze(log_text,      │                        │
│              │          program,       │                        │
│              │          context)       │                        │
│              └───────────┬─────────────┘                        │
│                          │                                       │
│           ┌──────────────┼──────────────┐                       │
│           ▼              ▼              ▼                       │
│    [Not Recoverable] [Recoverable]  [Max Retries]               │
│           │              │              │                       │
│           ▼              ▼              ▼                       │
│      Continue        Save to        Report &                    │
│      (normal fail)   Session        Continue                    │
│                          │                                       │
│                          │  • recovery_strategies[file] = flags │
│                          │  • force_retry_program = program     │
│                          │                                       │
└──────────────────────────┼───────────────────────────────────────┘
                           │
                           ▼
                    [Next Cycle]
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                        Planner                                   │
│                                                                  │
│   1. Check session for force_retry_program                      │
│   2. If set, select that program (skip normal planning)         │
│   3. Clear force_retry_program after use                        │
└─────────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                     CommandBuilder                               │
│                                                                  │
│   1. Check recovery_strategies for files being used             │
│   2. If file has recovery flags, merge into strategy            │
│   3. Build command with recovery parameter                      │
└─────────────────────────────────────────────────────────────────┘
```

### New Components

#### 1. ErrorAnalyzer (`agent/error_analyzer.py`)

Core class that:
- Matches error patterns in log text
- Extracts structured information (choices, keywords, affected files)
- Applies context-aware resolution logic
- Returns recovery strategy or None (if `auto_recovery=False`)

#### 2. Recoverable Errors Config (`knowledge/recoverable_errors.yaml`)

YAML configuration defining:
- Error patterns to match
- Extraction regex patterns
- Resolution strategies per error type
- Context-aware selection rules

#### 3. ErrorRecovery Data Class

```python
@dataclass
class ErrorRecovery:
    error_type: str           # e.g., "ambiguous_data_labels"
    affected_file: str        # e.g., "/path/to/data.mtz"
    flags: dict               # e.g., {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"}
    reason: str               # Human-readable explanation
    retry_program: str        # Program to force-retry
    selected_choice: str      # What was selected
    all_choices: list         # All available options
```

---

## Detailed Design

### 1. Error Pattern Matching

**Pattern for ambiguous data labels:**

```yaml
# knowledge/recoverable_errors.yaml
errors:
  ambiguous_data_labels:
    description: "MTZ file contains multiple suitable data arrays"
    
    # Patterns to identify this error (any match triggers)
    detection_patterns:
      - "Multiple equally suitable arrays"
      - "multiple.*suitable.*arrays"
      - "Please use.*to specify an unambiguous substring"
    
    # Regex to extract the keyword name from error message
    # Captures: "scaling.input.xray_data.obs_labels"
    keyword_extraction: 'Please use\s+(\S+)\s*\n?\s*to specify'
    
    # Regex to extract available choices
    # Applied line-by-line to find choice lines
    # Captures: filepath, labels
    choice_extraction: '^\s*(\S+\.mtz):(.+)$'
    
    # How to resolve this error
    resolution: context_based_label_selection
```

### 2. Context-Aware Label Selection

The selection logic considers:

**A. Program Type**

| Program | Data Preference | Reason |
|---------|-----------------|--------|
| phenix.autosol | Anomalous | SAD/MAD phasing needs Bijvoet pairs |
| phenix.hyss | Anomalous | Heavy atom substructure search |
| phenix.phaser (EP mode) | Anomalous | Experimental phasing |
| phenix.phaser (MR mode) | Merged | Molecular replacement |
| phenix.refine | Merged | Standard refinement |
| phenix.xtriage | Either | Analysis (prefer merged) |
| phenix.scale_and_merge | Anomalous if SAD | Depends on workflow |

**B. Workflow Context**

- If `project_advice` contains "SAD", "MAD", "anomalous" → prefer anomalous
- If `project_advice` contains "MR", "molecular replacement" → prefer merged
- If experiment involves `autosol` in history → anomalous workflow

**C. Label Pattern Heuristics**

```yaml
label_patterns:
  anomalous_indicators:
    - '\(\+\)'      # I(+), F(+)
    - '\(-\)'       # I(-), F(-)
    - 'anom'        # I_anom, F_anom
    - 'DANO'        # Anomalous difference
    - 'f_prime'     # f', f''
    - 'bijvoet'
    
  merged_indicators:
    - 'IMEAN'       # Mean intensity
    - 'FMEAN'       # Mean amplitude
    - 'F_obs'       # Observed F
    - 'I_obs'       # Observed I
    - '^F,'         # Simple F
    - '^I,'         # Simple I
    - 'merged'      # Explicitly merged
```

**Selection Algorithm:**

```python
def select_label(choices, program, context):
    """
    Select the best data label based on context.
    
    Returns: (selected_label, reason)
    """
    needs_anomalous = (
        program in ANOMALOUS_PROGRAMS or
        context.get("workflow_type") == "SAD" or
        "anomalous" in context.get("project_advice", "").lower()
    )
    
    anomalous_choices = [c for c in choices if has_anomalous_pattern(c)]
    merged_choices = [c for c in choices if has_merged_pattern(c)]
    
    if needs_anomalous:
        if anomalous_choices:
            return anomalous_choices[0], "Selected anomalous data for phasing workflow"
        else:
            # Fallback: warn but use first available
            return choices[0], "WARNING: No anomalous data found, using first available"
    else:
        if merged_choices:
            return merged_choices[0], "Selected merged data for refinement/MR"
        else:
            return choices[0], "Using first available data array"
```

### 3. Retry Tracking

To prevent infinite loops, track recovery attempts in session:

```python
# In session.data
{
    "recovery_attempts": {
        "ambiguous_data_labels": {
            "count": 1,
            "last_attempt": "2024-01-25T10:30:00",
            "files_tried": {
                "/path/to/data.mtz": ["I_CuKa(+)"]  # Track selections per file
            }
        }
    },
    # File-keyed recovery strategies (persist across cycles)
    "recovery_strategies": {
        "/path/to/data.mtz": {
            "flags": {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"},
            "program_scope": ["phenix.autosol", "phenix.phaser"],  # Optional restriction
            "reason": "Selected anomalous data for MRSAD workflow"
        }
    },
    # Force the planner to retry a specific program
    "force_retry_program": "phenix.autosol"
}
```

**Retry Logic:**

```python
MAX_RECOVERY_ATTEMPTS = 3

def should_attempt_recovery(session, error_type, affected_file):
    attempts = session.data.get("recovery_attempts", {}).get(error_type, {})
    count = attempts.get("count", 0)
    
    if count >= MAX_RECOVERY_ATTEMPTS:
        return False, f"Max recovery attempts ({MAX_RECOVERY_ATTEMPTS}) reached for {error_type}"
    
    # Check if we've already tried all options for this file
    files_tried = attempts.get("files_tried", {})
    tried_for_file = files_tried.get(affected_file, [])
    # Could add logic here to try different selections
    
    return True, None
```

### 4. Integration with Existing Code

#### A. ai_agent.py - `_run_single_cycle()`

```python
def _run_single_cycle(self, cycle_number, session, ...):
    ...
    result_text = self._execute_command(command, ...)
    
    # Check if the run failed
    run_failed = (
        result_text.startswith("FAILED:") or
        result_text.startswith("IGNORED:") or
        "FAILED" in result_text.upper()
    )
    
    if run_failed:
        # NEW: Check for recoverable errors (if auto_recovery enabled)
        auto_recovery = self.params.auto_recovery if hasattr(self.params, 'auto_recovery') else True
        
        if auto_recovery:
            from libtbx.langchain.agent.error_analyzer import ErrorAnalyzer
            
            analyzer = ErrorAnalyzer()
            recovery = analyzer.analyze(
                log_text=result_text,
                program=actual_program,
                context={
                    "project_advice": session.data.get("project_advice", ""),
                    "history": session.data.get("cycles", []),
                    "experiment_type": session.data.get("experiment_type"),
                },
                session=session  # For retry tracking
            )
            
            if recovery:
                # Store file-keyed recovery strategy
                strategies = session.data.setdefault("recovery_strategies", {})
                strategies[recovery.affected_file] = {
                    "flags": recovery.flags,
                    "program_scope": [actual_program],  # At minimum, this program
                    "reason": recovery.reason
                }
                
                # Force retry of the same program
                session.data["force_retry_program"] = actual_program
                session.save()
                
                # Loud notification block
                print("\n" + "=" * 50, file=self.logger)
                print("[NOTICE] DETECTED RECOVERABLE ERROR", file=self.logger)
                print(f"Issue: {recovery.error_type} in {recovery.affected_file}", file=self.logger)
                print(f"Action: Auto-selecting '{recovery.selected_choice}'", file=self.logger)
                print(f"Context: {recovery.reason}", file=self.logger)
                print("=" * 50 + "\n", file=self.logger)
                
                # Continue to next cycle (will retry same program)
                return False
        
        # Not recoverable - normal failure handling
        print("\n[NOTICE] Program failed - stop condition NOT triggered", file=self.logger)
        ...
```

#### B. Planner - Force Retry Check

In the planning phase, check if we should force a specific program:

```python
def plan(state):
    """Plan the next action."""
    session = state.get("session")
    
    # Check for forced retry (from error recovery)
    if session:
        force_program = session.data.pop("force_retry_program", None)
        if force_program:
            # Skip normal planning - retry the failed program
            state = _log(state, f"PLAN: Forced retry of {force_program} (error recovery)")
            return {
                **state,
                "intent": {
                    "program": force_program,
                    "reasoning": "Retrying with recovery flags after recoverable error",
                    "strategy": {}  # Strategy will be merged from recovery_strategies
                }
            }
    
    # Normal planning logic...
```

#### C. CommandBuilder - File-Keyed Strategy Lookup

```python
def build_command(self, program, context, strategy=None, session=None, files_used=None):
    """Build command with file-specific recovery strategies."""
    
    # Check for file-specific recovery strategies
    if session and files_used:
        recovery_strategies = session.data.get("recovery_strategies", {})
        
        for file_path in files_used:
            file_recovery = recovery_strategies.get(file_path)
            if file_recovery:
                # Check program scope (if specified)
                scope = file_recovery.get("program_scope", [])
                if not scope or program in scope:
                    # Merge file-specific flags into strategy
                    flags = file_recovery.get("flags", {})
                    if flags:
                        if strategy:
                            strategy = {**strategy, **flags}
                        else:
                            strategy = dict(flags)
                        
                        self._log(context, f"BUILD: Applying recovery flags for {file_path}: {flags}")
            session.save()
    
    # Continue with normal command building...
```

#### D. Session Data Schema Update

```python
# In session.py, document the new fields
"""
Session data fields for error recovery:

recovery_strategies: dict
    File-keyed recovery strategies that persist across cycles.
    Example: {
        "/path/to/data.mtz": {
            "flags": {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"},
            "program_scope": ["phenix.autosol"],  # Optional
            "reason": "Selected anomalous data for MRSAD"
        }
    }
    
force_retry_program: str
    If set, the planner will force this program on the next cycle,
    bypassing normal planning logic. Used after recoverable errors.
    Cleared after use.
    
recovery_attempts: dict
    Tracks retry counts per error type to prevent infinite loops.
    Example: {
        "ambiguous_data_labels": {
            "count": 1,
            "files_tried": {"/path/to/data.mtz": ["I_CuKa(+)"]}
        }
    }
"""
```

---

## Event System Integration

### Verbosity Strategy: "Loud on Action, Silent on Success"

Recovery events should be transparent but not cluttering.

### 1. Live Stream (During Execution)

When a recoverable error is detected, print a clear notification block:

```
==================================================
[NOTICE] DETECTED RECOVERABLE ERROR
Issue: Ambiguous data arrays in data.mtz
Action: Auto-selecting anomalous array 'I_CuKa(+)'
Context: Program requires anomalous data for phasing
==================================================
```

When the retry runs, the command string itself shows the fix:
```
phenix.autosol data.mtz seq.fa scaling.input.xray_data.obs_labels="I_CuKa(+)"
```

No additional logging needed - the flag is visible proof.

### 2. Event Log

Add a new `RECOVERY` event type (or use existing `NOTICE`):

```python
# In event_log.py
class EventType:
    ...
    RECOVERY = "recovery"  # Auto-recovery from error

# Event payload
{
    "type": "recovery",
    "error_type": "ambiguous_data_labels",
    "affected_file": "/path/to/data.mtz",
    "selected": "I_CuKa(+)",
    "reason": "Selected anomalous data for MRSAD workflow",
    "retry_program": "phenix.autosol"
}
```

### 3. Session Summary

**On Success:** Append a small badge to the successful cycle:

```
============================================================
WORKFLOW SUMMARY
============================================================

| Cycle | Program          | Result                              |
|-------|------------------|-------------------------------------|
| 1     | phenix.xtriage   | Success                             |
| 2     | phenix.autosol   | FAILED (ambiguous data)             |
| 3     | phenix.autosol   | Success (recovered from data error) |
| 4     | phenix.refine    | Success                             |
```

**On Failure (max retries reached):** Show full details:

```
============================================================
ERROR: Recovery Failed
============================================================
Error Type: ambiguous_data_labels
File: /path/to/data.mtz
Attempts: 3
Tried: I_CuKa(+), IMEAN_CuKa, F_CuKa
Suggestion: Manually specify data labels using:
  scaling.input.xray_data.obs_labels="YOUR_CHOICE"
============================================================
```

---

## Opt-Out Flag

### Parameter Definition

```python
# In programs/ai_agent.py parameter definitions
auto_recovery = True
    .type = bool
    .help = "Automatically attempt to recover from certain errors"
    .short_caption = "Auto-recovery"
```

### Behavior When Disabled

```python
if auto_recovery:
    recovery = analyzer.analyze(...)
    if recovery:
        # Apply recovery...
else:
    # Just log the suggestion but don't act
    suggestion = analyzer.get_suggestion(log_text, program)
    if suggestion:
        print(f"\n[INFO] Recoverable error detected but auto_recovery=False", file=self.logger)
        print(f"[INFO] Suggestion: {suggestion}", file=self.logger)
    # Proceed with normal failure handling
```

---

## File Changes Summary

### New Files

| File | Purpose |
|------|---------|
| `agent/error_analyzer.py` | Core error analysis and recovery logic |
| `knowledge/recoverable_errors.yaml` | Error patterns and resolution config |

### Modified Files

| File | Changes |
|------|---------|
| `programs/ai_agent.py` | Add error analysis after failure detection; add `auto_recovery` parameter |
| `agent/graph_nodes.py` | Add force_retry_program check in plan() node |
| `agent/command_builder.py` | Merge file-keyed recovery strategies into command |
| `agent/session.py` | Document new session data fields |
| `agent/event_log.py` | Add RECOVERY event type (optional) |

### Test Files

| File | Purpose |
|------|---------|
| `tests/test_error_analyzer.py` | Unit tests for error detection and recovery |

---

## Implementation Phases

### Phase 1: Core Infrastructure (This PR)

1. Create `error_analyzer.py` with `ErrorAnalyzer` class
2. Create `recoverable_errors.yaml` with ambiguous_data_labels pattern
3. Implement label selection logic (anomalous vs merged)
4. Add retry tracking in session (file-keyed)
5. Integrate into `ai_agent.py` with auto_recovery flag
6. Add force_retry_program handling in Planner
7. Update `CommandBuilder` to apply file-keyed recovery strategies
8. Write unit tests

### Phase 2: Additional Recoverable Errors (Future)

Candidates for future implementation:

| Error | Detection | Recovery |
|-------|-----------|----------|
| Missing resolution | "Please set the resolution" | Run mtriage/xtriage first |
| Space group mismatch | "Space group.*does not match" | May need user input |
| Missing sequence | "No sequence provided" | Check for .fa/.fasta files |
| Memory error | "Out of memory" | Reduce nproc, add memory flags |

### Phase 3: Enhanced Features (Future)

- User notification of recovery actions
- Recovery history in session summary
- Confidence levels for selections
- "Ask user" mode for low-confidence recoveries

---

## Example Walkthrough

### Scenario: MRSAD Workflow with Ambiguous Data

**Cycle 1:** User starts workflow
```
phenix.ai_agent original_files="lyso.mtz lyso.fa" project_advice="Solve using MRSAD"
```

**Cycle 2:** Agent runs autosol, which fails
```
phenix.autosol data.mtz seq.fa
→ FAILED: Multiple equally suitable arrays...
   Choices: IMEAN_CuKa, I_CuKa(+)...
   Please use scaling.input.xray_data.obs_labels
```

**Error Analysis:**
1. ErrorAnalyzer detects "Multiple equally suitable arrays"
2. Extracts keyword: `scaling.input.xray_data.obs_labels`
3. Extracts choices: `["IMEAN_CuKa,...", "I_CuKa(+),..."]`
4. Context: program=autosol, advice mentions "MRSAD"
5. Selection: `I_CuKa(+)` (anomalous data for phasing)
6. Saves to session: `recovery_strategy = {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"}`

**Cycle 3:** Agent retries autosol with fix
```
phenix.autosol data.mtz seq.fa scaling.input.xray_data.obs_labels="I_CuKa(+)"
→ SUCCESS
```

**Session Summary:**
```
Cycle 2: phenix.autosol - FAILED (ambiguous data)
         [Recovery: Selected anomalous data I_CuKa(+) for MRSAD workflow]
Cycle 3: phenix.autosol - SUCCESS
```

---

## Testing Strategy

### Unit Tests

```python
def test_detect_ambiguous_data_error():
    """Test that we correctly identify the error pattern."""
    log = """
    Multiple equally suitable arrays of observed xray data found.
    
    Possible choices:
      /path/data.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
      /path/data.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)
    
    Please use scaling.input.xray_data.obs_labels
    to specify an unambiguous substring of the target label.
    """
    
    analyzer = ErrorAnalyzer()
    result = analyzer._detect_error_type(log)
    
    assert result == "ambiguous_data_labels"

def test_extract_keyword_from_error():
    """Test keyword extraction from error message."""
    # ... test keyword_extraction regex
    
def test_extract_choices_from_error():
    """Test choice extraction from error message."""
    # ... test choice_extraction regex

def test_select_anomalous_for_autosol():
    """Test that autosol gets anomalous data."""
    choices = ["IMEAN_CuKa,SIGIMEAN_CuKa", "I_CuKa(+),SIGI_CuKa(+)"]
    
    selected, reason = select_label(choices, program="phenix.autosol", context={})
    
    assert "I_CuKa(+)" in selected
    assert "anomalous" in reason.lower()

def test_select_merged_for_refine():
    """Test that refine gets merged data."""
    choices = ["IMEAN_CuKa,SIGIMEAN_CuKa", "I_CuKa(+),SIGI_CuKa(+)"]
    
    selected, reason = select_label(choices, program="phenix.refine", context={})
    
    assert "IMEAN" in selected

def test_max_retry_limit():
    """Test that we give up after 3 attempts."""
    session = MockSession()
    session.data["recovery_attempts"] = {
        "ambiguous_data_labels": {"count": 3}
    }
    
    analyzer = ErrorAnalyzer()
    result = analyzer.analyze(log, program, context, session)
    
    assert result is None  # Should not attempt recovery

def test_recovery_strategy_applied():
    """Test that CommandBuilder applies recovery strategy."""
    # ... integration test
```

### Integration Tests

- Full cycle test: fail → recover → succeed
- Multiple recovery attempts
- Different programs get different selections
- Recovery strategy cleared after use

---

## Design Decisions (Resolved)

These questions were raised during review and have been resolved:

### 1. Should we log recovery attempts to the event system?

**Decision: YES, with severity grading**

- **Live Stream:** Log as `NOTICE` or `RECOVERY` event with clear block format
- **Session Summary:** Don't clutter - just append "(recovered from data error)" badge to successful step
- **On Failure:** Show full details of what was tried

### 2. What if the same error occurs with a different MTZ file?

**Decision: Key recovery by filename**

Recovery strategies are stored per-file:
```python
session.data["recovery_strategies"] = {
    "/path/to/data.mtz": {
        "flags": {"obs_labels": "I(+)"},
        "program_scope": ["phenix.autosol"]
    }
}
```

CommandBuilder checks if the file being used has recovery flags before applying them.

### 3. Should recovery be opt-out?

**Decision: YES**

Add `auto_recovery=True` parameter (default enabled). When `False`:
- ErrorAnalyzer logs suggestion but returns None
- Normal failure handling proceeds
- Useful for debugging the error handler itself

### 4. How verbose should recovery logging be?

**Decision: "Loud on Action, Silent on Success"**

- **On Detection:** Print clear notification block (always visible)
- **On Retry:** The command string shows the injected flag - no extra logging needed
- **On Success:** Just the "(recovered)" badge in summary

### 5. How to ensure retry runs the same program?

**Decision: Force retry via session flag**

Set `session.data["force_retry_program"] = "phenix.autosol"`. The Planner checks this first and bypasses normal planning if set.

---

## Timeline Estimate

| Task | Effort |
|------|--------|
| Create error_analyzer.py | 2-3 hours |
| Create recoverable_errors.yaml | 1 hour |
| Integrate into ai_agent.py | 1-2 hours |
| Add force_retry to Planner | 1 hour |
| Update CommandBuilder (file-keyed) | 1-2 hours |
| Write tests | 2-3 hours |
| Documentation | 1 hour |
| **Total** | **9-13 hours** |

---

## References

- Fix 13: Stop condition on failed runs (provides the failure detection hook)
- Fix 15: predict_and_build resolution handling (similar "graceful degradation" pattern)
- CommandBuilder strategy system (existing mechanism we'll leverage)
