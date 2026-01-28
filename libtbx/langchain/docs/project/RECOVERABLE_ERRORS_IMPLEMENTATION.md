# Recoverable Errors: Detailed Implementation Plan

## Executive Summary

This document provides a step-by-step implementation guide for the error recovery system described in `RECOVERABLE_ERRORS_PLAN.md`. The implementation is broken into 8 concrete steps that can be completed incrementally.

---

## Implementation Steps

### Step 1: Create `knowledge/recoverable_errors.yaml`

**Purpose:** Define error patterns, extraction rules, and resolution strategies in a maintainable YAML format.

**File:** `knowledge/recoverable_errors.yaml`

**Key Sections:**
1. `errors:` - Dictionary of error type definitions
2. `label_patterns:` - Patterns for identifying anomalous vs merged data
3. `program_data_preferences:` - Which programs need which data type

**Contents:**
```yaml
# PHENIX AI Agent - Recoverable Error Definitions
# =============================================
# This file defines error patterns that can be automatically recovered from.

errors:
  ambiguous_data_labels:
    description: "MTZ file contains multiple suitable data arrays"
    max_retries: 3
    
    # Patterns to detect this error (any match triggers)
    detection_patterns:
      - "Multiple equally suitable arrays"
      - "Please use.*to specify an unambiguous substring"
    
    # Regex to extract the parameter name from error message
    # Example: "Please use scaling.input.xray_data.obs_labels"
    keyword_extraction: 'Please use\s+(\S+)\s+to specify'
    
    # Regex to extract available choices (file:labels format)
    # Applied line-by-line
    choice_extraction: '^\s*(\S+\.mtz):(.+)$'
    
    resolution: context_based_label_selection

label_patterns:
  anomalous_indicators:
    - '\(\+\)'      # I(+), F(+)
    - '\(-\)'       # I(-), F(-)
    - 'anom'        # I_anom, F_anom
    - 'DANO'        # Anomalous difference
    - 'f_prime'     # f', f''
    - 'bijvoet'
    - 'PANO'        # Phased anomalous
    
  merged_indicators:
    - 'IMEAN'       # Mean intensity
    - 'FMEAN'       # Mean amplitude
    - 'F_obs'       # Observed F
    - 'I_obs'       # Observed I
    - '^F,'         # Simple F at start
    - '^I,'         # Simple I at start
    - 'merged'      # Explicitly merged

program_data_preferences:
  # Programs that need anomalous data (SAD/MAD phasing)
  anomalous:
    - phenix.autosol
    - phenix.hyss
    - phenix.phaser_ep   # Experimental phasing mode
    - phenix.scale_and_merge  # When doing SAD
    
  # Programs that prefer merged data
  merged:
    - phenix.refine
    - phenix.phaser      # MR mode (default)
    - phenix.molrep
    
  # Programs that can use either (prefer merged)
  either:
    - phenix.xtriage
    - phenix.maps
```

---

### Step 2: Create `agent/error_analyzer.py`

**Purpose:** Core class for detecting recoverable errors and determining recovery strategies.

**Key Classes:**
1. `ErrorRecovery` - Dataclass holding recovery information
2. `ErrorAnalyzer` - Main class for error detection and resolution

**File Structure:**
```python
"""
Error Analyzer - Automatic Recovery from Recoverable Errors.

This module detects structured errors in PHENIX log output and determines
appropriate recovery strategies.

Currently supported errors:
- ambiguous_data_labels: Multiple data arrays in MTZ file
"""

from __future__ import absolute_import, division, print_function

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple

# YAML loading
from libtbx.langchain.knowledge.yaml_loader import get_yaml_path, load_yaml_file


@dataclass
class ErrorRecovery:
    """Information needed to recover from an error."""
    error_type: str           # e.g., "ambiguous_data_labels"
    affected_file: str        # e.g., "/path/to/data.mtz"
    flags: Dict[str, str]     # e.g., {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"}
    reason: str               # Human-readable explanation
    retry_program: str        # Program to force-retry
    selected_choice: str      # What was selected
    all_choices: List[str]    # All available options


class ErrorAnalyzer:
    """
    Analyzes program errors and determines recovery strategies.
    
    Usage:
        analyzer = ErrorAnalyzer()
        recovery = analyzer.analyze(
            log_text="...",
            program="phenix.autosol",
            context={"project_advice": "..."},
            session=session
        )
        
        if recovery:
            # Apply recovery.flags and retry recovery.retry_program
    """
    
    def __init__(self):
        """Initialize with configuration from YAML."""
        self._config = self._load_config()
        self._label_patterns = self._config.get("label_patterns", {})
        self._program_prefs = self._config.get("program_data_preferences", {})
    
    def _load_config(self) -> dict:
        """Load recoverable errors configuration."""
        yaml_path = get_yaml_path("recoverable_errors.yaml")
        if yaml_path and os.path.exists(yaml_path):
            return load_yaml_file(yaml_path) or {}
        return {}
    
    def analyze(self, log_text: str, program: str, 
                context: Dict[str, Any], session) -> Optional[ErrorRecovery]:
        """
        Analyze log text for recoverable errors.
        
        Args:
            log_text: Full log/error text from program
            program: Program name that failed
            context: Dict with project_advice, history, experiment_type
            session: Session object for tracking retries
            
        Returns:
            ErrorRecovery if a recovery is possible, None otherwise
        """
        # 1. Detect error type
        error_type = self._detect_error_type(log_text)
        if not error_type:
            return None
        
        # 2. Check retry limits
        can_retry, limit_reason = self._check_retry_limits(
            session, error_type
        )
        if not can_retry:
            self._log_max_retries(error_type, limit_reason)
            return None
        
        # 3. Extract structured information
        error_info = self._extract_error_info(log_text, error_type)
        if not error_info:
            return None
        
        # 4. Determine recovery strategy
        recovery = self._determine_recovery(
            error_type, error_info, program, context
        )
        
        # 5. Update retry tracking in session
        if recovery:
            self._update_retry_tracking(session, error_type, 
                                        recovery.affected_file,
                                        recovery.selected_choice)
        
        return recovery
    
    def get_suggestion(self, log_text: str, program: str) -> Optional[str]:
        """
        Get a human-readable suggestion without attempting recovery.
        
        Used when auto_recovery=False.
        """
        error_type = self._detect_error_type(log_text)
        if not error_type:
            return None
        
        error_info = self._extract_error_info(log_text, error_type)
        if not error_info:
            return None
        
        # Generate suggestion text
        if error_type == "ambiguous_data_labels":
            keyword = error_info.get("keyword", "obs_labels")
            choices = error_info.get("choices", [])
            choices_str = ", ".join(choices[:3])
            return (f"Ambiguous data labels. Add {keyword}=\"YOUR_CHOICE\" "
                    f"where YOUR_CHOICE is one of: {choices_str}")
        
        return None
    
    # =========================================================================
    # DETECTION
    # =========================================================================
    
    def _detect_error_type(self, log_text: str) -> Optional[str]:
        """Detect which recoverable error type (if any) is present."""
        errors = self._config.get("errors", {})
        
        for error_type, error_def in errors.items():
            patterns = error_def.get("detection_patterns", [])
            for pattern in patterns:
                if re.search(pattern, log_text, re.IGNORECASE):
                    return error_type
        
        return None
    
    # =========================================================================
    # EXTRACTION
    # =========================================================================
    
    def _extract_error_info(self, log_text: str, 
                            error_type: str) -> Optional[Dict[str, Any]]:
        """Extract structured information from error message."""
        error_def = self._config.get("errors", {}).get(error_type, {})
        
        if error_type == "ambiguous_data_labels":
            return self._extract_ambiguous_labels_info(log_text, error_def)
        
        return None
    
    def _extract_ambiguous_labels_info(self, log_text: str, 
                                       error_def: dict) -> Optional[Dict[str, Any]]:
        """Extract keyword and choices from ambiguous data labels error."""
        result = {
            "keyword": None,
            "choices": [],
            "affected_file": None,
            "choice_details": []  # [(file, labels), ...]
        }
        
        # Extract keyword name
        keyword_pattern = error_def.get("keyword_extraction", "")
        if keyword_pattern:
            match = re.search(keyword_pattern, log_text)
            if match:
                result["keyword"] = match.group(1)
        
        # Extract choices (file:labels pairs)
        choice_pattern = error_def.get("choice_extraction", "")
        if choice_pattern:
            for line in log_text.split('\n'):
                match = re.match(choice_pattern, line.strip())
                if match:
                    file_path = match.group(1)
                    labels = match.group(2)
                    result["choices"].append(labels)
                    result["choice_details"].append((file_path, labels))
                    if not result["affected_file"]:
                        result["affected_file"] = file_path
        
        if not result["choices"]:
            return None
        
        return result
    
    # =========================================================================
    # RESOLUTION
    # =========================================================================
    
    def _determine_recovery(self, error_type: str, error_info: Dict[str, Any],
                           program: str, context: Dict[str, Any]) -> Optional[ErrorRecovery]:
        """Determine the recovery strategy."""
        if error_type == "ambiguous_data_labels":
            return self._resolve_ambiguous_labels(error_info, program, context)
        
        return None
    
    def _resolve_ambiguous_labels(self, error_info: Dict[str, Any],
                                  program: str, 
                                  context: Dict[str, Any]) -> Optional[ErrorRecovery]:
        """Resolve ambiguous data labels by selecting appropriate array."""
        choices = error_info.get("choices", [])
        keyword = error_info.get("keyword")
        affected_file = error_info.get("affected_file")
        
        if not choices or not keyword:
            return None
        
        # Determine if we need anomalous data
        needs_anomalous = self._needs_anomalous_data(program, context)
        
        # Classify choices
        anomalous_choices = [c for c in choices if self._is_anomalous_label(c)]
        merged_choices = [c for c in choices if self._is_merged_label(c)]
        
        # Select appropriate choice
        if needs_anomalous:
            if anomalous_choices:
                selected = anomalous_choices[0]
                reason = f"Selected anomalous data for {program} (phasing workflow)"
            else:
                # No anomalous available, use first
                selected = choices[0]
                reason = f"WARNING: No anomalous data found, using first available"
        else:
            if merged_choices:
                selected = merged_choices[0]
                reason = f"Selected merged data for {program}"
            else:
                selected = choices[0]
                reason = f"Using first available data array"
        
        # Extract the main label (first column) for the flag value
        label_value = self._extract_main_label(selected)
        
        return ErrorRecovery(
            error_type="ambiguous_data_labels",
            affected_file=affected_file,
            flags={keyword: label_value},
            reason=reason,
            retry_program=program,
            selected_choice=selected,
            all_choices=choices
        )
    
    def _needs_anomalous_data(self, program: str, context: Dict[str, Any]) -> bool:
        """Determine if the program/context needs anomalous data."""
        # Check program preference
        anomalous_programs = self._program_prefs.get("anomalous", [])
        if program in anomalous_programs or program.replace("phenix.", "") in [
            p.replace("phenix.", "") for p in anomalous_programs
        ]:
            return True
        
        # Check project advice
        advice = context.get("project_advice", "").lower()
        if any(kw in advice for kw in ["sad", "mad", "anomalous", "mrsad"]):
            return True
        
        # Check experiment type
        exp_type = context.get("experiment_type")
        if exp_type in ["sad", "mad"]:
            return True
        
        # Check history for phasing programs
        history = context.get("history", [])
        phasing_indicators = ["autosol", "hyss", "phaser_ep"]
        for entry in history:
            if isinstance(entry, dict):
                prog = entry.get("program", "") + entry.get("command", "")
                if any(ind in prog.lower() for ind in phasing_indicators):
                    return True
        
        return False
    
    def _is_anomalous_label(self, label: str) -> bool:
        """Check if label string indicates anomalous data."""
        patterns = self._label_patterns.get("anomalous_indicators", [])
        for pattern in patterns:
            if re.search(pattern, label, re.IGNORECASE):
                return True
        return False
    
    def _is_merged_label(self, label: str) -> bool:
        """Check if label string indicates merged data."""
        patterns = self._label_patterns.get("merged_indicators", [])
        for pattern in patterns:
            if re.search(pattern, label, re.IGNORECASE):
                return True
        return False
    
    def _extract_main_label(self, label_string: str) -> str:
        """Extract the main label from a comma-separated label string.
        
        E.g., "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-)" -> "I_CuKa(+)"
        """
        parts = label_string.split(",")
        return parts[0].strip() if parts else label_string
    
    # =========================================================================
    # RETRY TRACKING
    # =========================================================================
    
    def _check_retry_limits(self, session, error_type: str) -> Tuple[bool, Optional[str]]:
        """Check if we've exceeded retry limits for this error type."""
        error_def = self._config.get("errors", {}).get(error_type, {})
        max_retries = error_def.get("max_retries", 3)
        
        attempts = session.data.get("recovery_attempts", {}).get(error_type, {})
        count = attempts.get("count", 0)
        
        if count >= max_retries:
            return False, f"Max recovery attempts ({max_retries}) reached"
        
        return True, None
    
    def _update_retry_tracking(self, session, error_type: str, 
                               affected_file: str, selected_choice: str):
        """Update session with new retry attempt."""
        recovery_attempts = session.data.setdefault("recovery_attempts", {})
        type_attempts = recovery_attempts.setdefault(error_type, {
            "count": 0,
            "files_tried": {}
        })
        
        type_attempts["count"] = type_attempts.get("count", 0) + 1
        files_tried = type_attempts.setdefault("files_tried", {})
        file_choices = files_tried.setdefault(affected_file, [])
        file_choices.append(selected_choice)
    
    def _log_max_retries(self, error_type: str, reason: str):
        """Log when max retries reached."""
        print(f"\n[WARNING] {error_type}: {reason}")
```

---

### Step 3: Update Session Data Schema

**File:** `agent/session.py`

**Changes:**
1. Document new fields in `_init_new_session()`
2. Add helper method `get_recovery_strategy(file_path)`
3. Add helper method `clear_force_retry()`

**New Fields in `self.data`:**
```python
# Error recovery tracking
"recovery_attempts": {},    # {error_type: {count: N, files_tried: {...}}}
"recovery_strategies": {},  # {file_path: {flags: {...}, program_scope: [...], reason: "..."}}
"force_retry_program": None # Program name to force on next cycle
```

**New Methods:**
```python
def get_recovery_strategy(self, file_path: str) -> Optional[Dict]:
    """Get recovery strategy for a specific file."""
    return self.data.get("recovery_strategies", {}).get(file_path)

def clear_force_retry(self) -> Optional[str]:
    """Clear and return the force_retry_program flag."""
    return self.data.pop("force_retry_program", None)

def set_recovery_strategy(self, file_path: str, flags: Dict, 
                          program: str, reason: str):
    """Set recovery strategy for a file."""
    strategies = self.data.setdefault("recovery_strategies", {})
    strategies[file_path] = {
        "flags": flags,
        "program_scope": [program],
        "reason": reason
    }
```

---

### Step 4: Update `ai_agent.py` - Error Analysis Integration

**File:** `programs/ai_agent.py`

**Changes:**

1. Add `auto_recovery` parameter to phil scope:
```python
auto_recovery = True
  .type = bool
  .help = "Automatically attempt to recover from certain structured errors"
  .short_caption = "Auto-recovery"
```

2. Modify `_run_single_cycle()` to analyze failures:
```python
if run_failed:
    # NEW: Check for recoverable errors
    auto_recovery = getattr(self.params, 'auto_recovery', True)
    
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
            session=session
        )
        
        if recovery:
            # Store file-keyed recovery strategy
            session.set_recovery_strategy(
                recovery.affected_file,
                recovery.flags,
                actual_program,
                recovery.reason
            )
            
            # Force retry of the same program
            session.data["force_retry_program"] = actual_program
            session.save()
            
            # Loud notification block
            self._print_recovery_notice(recovery)
            
            # Continue to next cycle (will retry same program)
            return False
    else:
        # Log suggestion but don't act
        suggestion = analyzer.get_suggestion(result_text, actual_program)
        if suggestion:
            print(f"\n[INFO] Recoverable error detected (auto_recovery=False)", 
                  file=self.logger)
            print(f"[INFO] Suggestion: {suggestion}", file=self.logger)
    
    # Not recoverable - normal failure handling
    print("\n[NOTICE] Program failed - stop condition NOT triggered", 
          file=self.logger)
    ...
```

3. Add helper method for notification:
```python
def _print_recovery_notice(self, recovery):
    """Print formatted recovery notification."""
    print("\n" + "=" * 60, file=self.logger)
    print("[NOTICE] DETECTED RECOVERABLE ERROR", file=self.logger)
    print(f"Issue: {recovery.error_type}", file=self.logger)
    print(f"File: {os.path.basename(recovery.affected_file)}", file=self.logger)
    print(f"Action: Selecting '{recovery.selected_choice}'", file=self.logger)
    print(f"Reason: {recovery.reason}", file=self.logger)
    print("=" * 60 + "\n", file=self.logger)
```

---

### Step 5: Update `graph_nodes.py` - Force Retry in Planner

**File:** `agent/graph_nodes.py`

**Changes:** Modify `plan()` function to check for force_retry_program early.

**Insert at the beginning of `plan()` function (after docstring, before auto-stop check):**

```python
def plan(state):
    """..."""
    
    # 0. Check for forced retry (from error recovery)
    session_info = state.get("session_info", {})
    session = session_info.get("session_object")  # Need to pass session object
    
    if session:
        force_program = session.clear_force_retry()
        if force_program:
            state = _log(state, f"PLAN: Forced retry of {force_program} (error recovery)")
            
            # Emit event for transparency
            state = _emit(state, EventType.PLAN_DECISION,
                          program=force_program,
                          reason="error_recovery",
                          forced=True)
            
            return {
                **state,
                "intent": {
                    "program": force_program,
                    "reasoning": "Retrying with recovery flags after recoverable error",
                    "strategy": {},  # Strategy merged from recovery_strategies in BUILD
                    "files": {},
                    "forced_retry": True
                }
            }
    
    # 1. Check for auto-stop from metrics trend
    ...
```

**Alternative approach** (if session object isn't easily available in state):

Store force_retry in state itself during _query_agent_for_command():
```python
# In session_info passed to state
session_info = {
    ...
    "force_retry_program": session.data.get("force_retry_program"),
}
```

Then in plan():
```python
force_program = session_info.get("force_retry_program")
if force_program:
    # Clear it from session (via callback or flag)
    state = {**state, "clear_force_retry": True}
    ...
```

---

### Step 6: Update `command_builder.py` - Apply Recovery Strategies

**File:** `agent/command_builder.py`

**Changes:**

1. Extend `CommandContext` with recovery info:
```python
@dataclass
class CommandContext:
    ...
    # Recovery strategies (file-keyed)
    recovery_strategies: Dict[str, Dict] = field(default_factory=dict)
    
    @classmethod
    def from_state(cls, state: dict) -> 'CommandContext':
        ...
        # Get recovery strategies from session_info
        recovery_strategies = session_info.get("recovery_strategies", {})
        
        return cls(
            ...
            recovery_strategies=recovery_strategies,
        )
```

2. Modify `_build_strategy()` to merge recovery flags:
```python
def _build_strategy(self, program: str, context: CommandContext) -> Dict[str, Any]:
    """Build strategy flags from context."""
    strategy = {}
    
    # Start with LLM suggestions (if any)
    if context.llm_strategy:
        strategy.update(context.llm_strategy)
    
    # NEW: Merge recovery strategies for files being used
    if context.recovery_strategies:
        # We need to know which files are being used
        # This is called AFTER file selection, so we need to pass files
        # Actually, better to do this in build() after _select_files()
        pass
    
    return strategy
```

3. Better: Modify `build()` to inject recovery flags:
```python
def build(self, program: str, available_files: List[str],
          context: CommandContext) -> Optional[str]:
    ...
    
    # 1. Select files
    files = self._select_files(program, available_files, context)
    if files is None:
        return None
    
    # 2. Build strategy
    strategy = self._build_strategy(program, context)
    
    # 2.5 NEW: Apply file-specific recovery strategies
    strategy = self._apply_recovery_strategies(
        program, files, strategy, context
    )
    
    # 3. Apply invariants
    ...

def _apply_recovery_strategies(self, program: str, files: Dict[str, Any],
                               strategy: Dict[str, Any], 
                               context: CommandContext) -> Dict[str, Any]:
    """Merge recovery strategies for files being used."""
    if not context.recovery_strategies:
        return strategy
    
    # Get all file paths from selected files
    file_paths = set()
    for slot, value in files.items():
        if isinstance(value, str) and os.path.isfile(value):
            file_paths.add(value)
            file_paths.add(os.path.abspath(value))
    
    # Check each recovery strategy
    for file_path, recovery in context.recovery_strategies.items():
        # Check if this file is being used (match by basename or full path)
        file_match = (
            file_path in file_paths or
            os.path.basename(file_path) in [os.path.basename(f) for f in file_paths]
        )
        
        if file_match:
            # Check program scope
            scope = recovery.get("program_scope", [])
            if not scope or program in scope:
                # Merge flags
                flags = recovery.get("flags", {})
                if flags:
                    strategy = {**strategy, **flags}
                    self._log(context, 
                              f"BUILD: Applied recovery flags for {os.path.basename(file_path)}: {flags}")
    
    return strategy
```

---

### Step 7: Pass Recovery Context Through Pipeline

**Files to update:**

1. **`ai_agent.py`** - Include recovery_strategies in session_info:
```python
session_info = {
    ...
    "recovery_strategies": session.data.get("recovery_strategies", {}),
    "force_retry_program": session.data.get("force_retry_program"),
}
```

2. **`graph_nodes.py`** - Handle clearing force_retry after use:
```python
# In plan() after handling force_retry:
if force_program:
    # Tell session to clear the flag
    session.data.pop("force_retry_program", None)
    session.save()
```

3. **`session.py`** - Ensure recovery data persists:
```python
# In save():
# recovery_attempts and recovery_strategies are already in self.data
# so they get serialized automatically
```

---

### Step 8: Create Unit Tests

**File:** `tests/test_error_analyzer.py`

**Test Cases:**

```python
"""
Tests for error_analyzer.py - Recoverable Error Handling.
"""

import unittest
import os
import tempfile
import json

from libtbx.langchain.agent.error_analyzer import ErrorAnalyzer, ErrorRecovery


class TestErrorDetection(unittest.TestCase):
    """Test error type detection."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_detect_ambiguous_data_error(self):
        """Test detection of ambiguous data labels error."""
        log = """
        Multiple equally suitable arrays of observed xray data found.
        
        Possible choices:
          /path/data.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
          /path/data.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)
        
        Please use scaling.input.xray_data.obs_labels
        to specify an unambiguous substring of the target label.
        """
        
        error_type = self.analyzer._detect_error_type(log)
        self.assertEqual(error_type, "ambiguous_data_labels")
    
    def test_no_error_detected(self):
        """Test that normal logs don't trigger false positives."""
        log = "phenix.autosol completed successfully. R-factor = 0.20"
        
        error_type = self.analyzer._detect_error_type(log)
        self.assertIsNone(error_type)


class TestKeywordExtraction(unittest.TestCase):
    """Test extraction of parameter keywords from errors."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_extract_keyword(self):
        """Test keyword extraction from error message."""
        log = """
        Please use scaling.input.xray_data.obs_labels
        to specify an unambiguous substring.
        """
        
        info = self.analyzer._extract_ambiguous_labels_info(
            log, 
            self.analyzer._config.get("errors", {}).get("ambiguous_data_labels", {})
        )
        
        self.assertIsNotNone(info)
        self.assertEqual(info["keyword"], "scaling.input.xray_data.obs_labels")


class TestChoiceExtraction(unittest.TestCase):
    """Test extraction of available choices."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_extract_choices(self):
        """Test extraction of label choices."""
        log = """
        Possible choices:
          /data/test.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
          /data/test.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)
        
        Please use obs_labels to specify.
        """
        
        info = self.analyzer._extract_ambiguous_labels_info(
            log,
            self.analyzer._config.get("errors", {}).get("ambiguous_data_labels", {})
        )
        
        self.assertIsNotNone(info)
        self.assertEqual(len(info["choices"]), 2)
        self.assertIn("IMEAN_CuKa,SIGIMEAN_CuKa", info["choices"])


class TestLabelClassification(unittest.TestCase):
    """Test anomalous vs merged label classification."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_anomalous_label_detection(self):
        """Test detection of anomalous data labels."""
        # These should be detected as anomalous
        anomalous_labels = [
            "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)",
            "F(+),SIGF(+),F(-),SIGF(-)",
            "DANO,SIGDANO",
            "I_anom,SIGI_anom",
        ]
        
        for label in anomalous_labels:
            self.assertTrue(
                self.analyzer._is_anomalous_label(label),
                f"Expected {label} to be anomalous"
            )
    
    def test_merged_label_detection(self):
        """Test detection of merged data labels."""
        merged_labels = [
            "IMEAN_CuKa,SIGIMEAN_CuKa",
            "FMEAN,SIGFMEAN",
            "F_obs,SIGF_obs",
        ]
        
        for label in merged_labels:
            self.assertTrue(
                self.analyzer._is_merged_label(label),
                f"Expected {label} to be merged"
            )


class TestProgramDataPreference(unittest.TestCase):
    """Test program-based data preference."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_autosol_needs_anomalous(self):
        """Test that autosol prefers anomalous data."""
        needs = self.analyzer._needs_anomalous_data(
            "phenix.autosol", {}
        )
        self.assertTrue(needs)
    
    def test_refine_prefers_merged(self):
        """Test that refine prefers merged data."""
        needs = self.analyzer._needs_anomalous_data(
            "phenix.refine", {}
        )
        self.assertFalse(needs)
    
    def test_context_override_sad(self):
        """Test that SAD in project advice triggers anomalous preference."""
        needs = self.analyzer._needs_anomalous_data(
            "phenix.phaser",  # Normally prefers merged
            {"project_advice": "Solve structure using SAD phasing"}
        )
        self.assertTrue(needs)


class TestLabelSelection(unittest.TestCase):
    """Test label selection logic."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_select_anomalous_for_autosol(self):
        """Test that autosol gets anomalous data."""
        choices = [
            "IMEAN_CuKa,SIGIMEAN_CuKa", 
            "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)"
        ]
        
        error_info = {
            "keyword": "obs_labels",
            "choices": choices,
            "affected_file": "/path/data.mtz"
        }
        
        recovery = self.analyzer._resolve_ambiguous_labels(
            error_info, "phenix.autosol", {}
        )
        
        self.assertIsNotNone(recovery)
        self.assertIn("I_CuKa(+)", recovery.selected_choice)
        self.assertIn("anomalous", recovery.reason.lower())
    
    def test_select_merged_for_refine(self):
        """Test that refine gets merged data."""
        choices = [
            "IMEAN_CuKa,SIGIMEAN_CuKa",
            "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)"
        ]
        
        error_info = {
            "keyword": "obs_labels",
            "choices": choices,
            "affected_file": "/path/data.mtz"
        }
        
        recovery = self.analyzer._resolve_ambiguous_labels(
            error_info, "phenix.refine", {}
        )
        
        self.assertIsNotNone(recovery)
        self.assertIn("IMEAN", recovery.selected_choice)


class TestRetryLimits(unittest.TestCase):
    """Test retry limit enforcement."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_max_retry_limit(self):
        """Test that we give up after max attempts."""
        # Create mock session
        class MockSession:
            def __init__(self):
                self.data = {
                    "recovery_attempts": {
                        "ambiguous_data_labels": {"count": 3}
                    }
                }
        
        session = MockSession()
        can_retry, reason = self.analyzer._check_retry_limits(
            session, "ambiguous_data_labels"
        )
        
        self.assertFalse(can_retry)
        self.assertIn("Max", reason)


class TestFullRecoveryFlow(unittest.TestCase):
    """Integration test for full recovery flow."""
    
    def setUp(self):
        self.analyzer = ErrorAnalyzer()
    
    def test_full_analyze_flow(self):
        """Test complete analyze() call."""
        log = """
        Multiple equally suitable arrays of observed xray data found.
        
        Possible choices:
          /data/lyso.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
          /data/lyso.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)
        
        Please use scaling.input.xray_data.obs_labels
        to specify an unambiguous substring of the target label.
        """
        
        class MockSession:
            def __init__(self):
                self.data = {}
        
        session = MockSession()
        
        recovery = self.analyzer.analyze(
            log_text=log,
            program="phenix.autosol",
            context={"project_advice": "MRSAD phasing"},
            session=session
        )
        
        self.assertIsNotNone(recovery)
        self.assertEqual(recovery.error_type, "ambiguous_data_labels")
        self.assertEqual(recovery.retry_program, "phenix.autosol")
        self.assertIn("I_CuKa(+)", recovery.selected_choice)
        self.assertIn("scaling.input.xray_data.obs_labels", recovery.flags)


if __name__ == '__main__':
    unittest.main()
```

---

## Implementation Order

Execute the steps in this order for incremental progress:

| Order | Step | Estimated Time | Dependencies |
|-------|------|----------------|--------------|
| 1 | Create recoverable_errors.yaml | 30 min | None |
| 2 | Create error_analyzer.py | 2-3 hours | Step 1 |
| 3 | Create test_error_analyzer.py | 1-2 hours | Step 2 |
| 4 | Update session.py | 30 min | None |
| 5 | Update ai_agent.py | 1-2 hours | Steps 2, 4 |
| 6 | Update graph_nodes.py | 1 hour | Steps 4, 5 |
| 7 | Update command_builder.py | 1-2 hours | Step 4 |
| 8 | Integration testing | 1-2 hours | All above |

**Total: 8-12 hours**

---

## Verification Checklist

After implementation, verify:

- [ ] ErrorAnalyzer correctly detects ambiguous_data_labels errors
- [ ] Keyword extraction works with different PHENIX error formats
- [ ] Choice extraction handles various MTZ label formats
- [ ] autosol/hyss get anomalous data, refine gets merged
- [ ] SAD/MAD project advice triggers anomalous preference
- [ ] Retry limits prevent infinite loops
- [ ] force_retry_program makes planner retry the same program
- [ ] Recovery flags appear in the rebuilt command
- [ ] auto_recovery=False shows suggestion but doesn't act
- [ ] Session data persists recovery state across cycles
- [ ] All unit tests pass
- [ ] Integration test: fail → recover → succeed works end-to-end
