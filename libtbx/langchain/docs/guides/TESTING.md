# Testing Guide

This document describes the testing infrastructure for the PHENIX AI Agent.

## Overview

The test suite uses a **cctbx-style testing approach** with plain functions and fail-fast behavior, meaning the first assertion failure stops execution with a full traceback. This matches the testing conventions used throughout the PHENIX/cctbx ecosystem.

## Quick Start

```bash
# Run all tests
cd improved_agent_v2v1
python tests/run_all_tests.py

# Run quick mode (standalone tests only, no PHENIX required)
python tests/run_all_tests.py --quick

# Run a single test file
python tests/tst_directive_extractor.py

# Run tests matching a pattern
python tests/run_all_tests.py --pattern "directive"
```

## Test Architecture

### Fail-Fast Behavior

Unlike Python's standard `unittest` which collects all failures, our tests crash immediately on the first failure:

```python
def test_something():
    assert_equal(result, expected)  # If this fails, execution stops here
    assert_true(condition)          # This line never runs if above fails
```

This provides immediate feedback with full tracebacks, making debugging faster.

### Assert Helper Functions

All tests use helper functions from `tests/tst_utils.py` instead of `unittest.TestCase` assertions:

| Function | Purpose |
|----------|---------|
| `assert_equal(a, b)` | Check `a == b` |
| `assert_not_equal(a, b)` | Check `a != b` |
| `assert_true(x)` | Check `x` is truthy |
| `assert_false(x)` | Check `x` is falsy |
| `assert_none(x)` | Check `x is None` |
| `assert_not_none(x)` | Check `x is not None` |
| `assert_in(item, container)` | Check `item in container` |
| `assert_not_in(item, container)` | Check `item not in container` |
| `assert_is(a, b)` | Check `a is b` |
| `assert_is_not(a, b)` | Check `a is not b` |
| `assert_is_instance(obj, cls)` | Check `isinstance(obj, cls)` |
| `assert_greater(a, b)` | Check `a > b` |
| `assert_greater_equal(a, b)` | Check `a >= b` |
| `assert_less(a, b)` | Check `a < b` |
| `assert_less_equal(a, b)` | Check `a <= b` |
| `assert_almost_equal(a, b, places)` | Check floats equal to N decimal places |
| `assert_raises(exc, func, *args)` | Check function raises exception |
| `assert_startswith(s, prefix)` | Check string starts with prefix |
| `assert_endswith(s, suffix)` | Check string ends with suffix |
| `assert_contains(text, substring)` | Check substring in text |

### Test File Structure

Each test file follows this structure:

```python
"""
Tests for module_name.

Run with: python tests/tst_module_name.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    run_tests_with_fail_fast
)

from agent.module_name import function_to_test


# =============================================================================
# SECTION NAME
# =============================================================================

def test_basic_functionality():
    """Description of what this tests."""
    result = function_to_test("input")
    assert_equal(result, "expected")


def test_edge_case():
    """Test edge case handling."""
    result = function_to_test(None)
    assert_none(result)


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
```

### Test Discovery

The `run_tests_with_fail_fast()` function automatically discovers tests:

1. Finds all functions starting with `test_` in the calling module
2. Also finds `unittest.TestCase` classes (for backward compatibility)
3. Runs them in sorted order
4. Stops at first failure with full traceback

## Test Files

| File | Tests | Description |
|------|-------|-------------|
| `tst_directive_extractor.py` | 73 | Directive extraction and validation |
| `tst_transport.py` | 66 | REST transport sanitization |
| `tst_workflow_state.py` | 66 | Workflow state detection, done flags, refine counting, MR-SAD |
| `tst_best_files_tracker.py` | 50 | File tracking, scoring, model stage detection |
| `tst_directive_validator.py` | 48 | Intent validation and modification |
| `tst_new_programs.py` | 44 | YAML config for new programs (polder, map_sharpening, autobuild_denmod) |
| `tst_error_analyzer.py` | 43 | Error recovery system |
| `tst_langchain_tools.py` | 34 | Legacy module tests (core, analysis, RAG, validation, prompts, memory) |
| `tst_sanity_checker.py` | 29 | Sanity check logic |
| `tst_yaml_config.py` | 29 | YAML configuration validation |
| `tst_file_categorization.py` | 26 | File type detection (including denmod_mtz) |
| `tst_api_schema.py` | 26 | API request/response validation |
| `tst_decision_flow.py` | 21 | Directive flow architecture |
| `tst_metrics_analyzer.py` | 19 | Metric extraction and trends |
| `tst_pattern_manager.py` | 17 | Pattern management |
| `tst_directives_integration.py` | 16 | End-to-end directive tests |
| `tst_command_builder.py` | 22 | Unified command generation, MR-SAD partpdb_file |
| `tst_event_system.py` | 13 | Event logging system |
| `tst_program_registration.py` | 13 | Program registry tests |
| `tst_integration.py` | 13 | End-to-end workflow tests |
| `tst_file_utils.py` | 12 | Shared file classification utilities |
| `tst_session_directives.py` | 12 | Session-level directive handling |
| `tst_state_serialization.py` | 12 | State packaging/unpackaging |
| `tst_summary_display.py` | 12 | Summary formatting |
| `tst_advice_preprocessing.py` | 11 | README discovery, advice processing |
| `tst_metric_patterns.py` | 11 | Log parsing patterns |
| `tst_history_analysis.py` | 10 | History analysis, anomalous workflow support |
| `tst_session_summary.py` | 10 | Session summary generation, STOP cycle exclusion |
| `tst_yaml_tools.py` | 9 | YAML validation and inspection |
| `tst_docs_tools.py` | 8 | Documentation generation |
| `tst_dry_run.py` | 8 | Dry run manager functionality |
| `tst_session_tools.py` | 7 | Session management utilities |
| `tst_template.py` | 5 | Template builder |
| `tst_phaser_multimodel.py` | 3 | Phaser multi-model handling |
| `tst_utils.py` | 2 | Assert helpers |

Total: **800 tests across 35 files**

### Key Tests for Recent Fixes

| Test | File | What It Verifies |
|------|------|------------------|
| `tst_automation_path_in_workflow_state` | tst_workflow_state.py | automation_path correctly set in state |
| `tst_stepwise_mode_blocks_predict_and_build_after_prediction` | tst_workflow_state.py | predict_and_build blocked in stepwise mode after prediction |
| `tst_autobuild_beats_earlier_refine_with_better_metrics` | tst_best_files_tracker.py | Autobuild with better R-free becomes best model |
| `tst_yaml_stage_scores` | tst_best_files_tracker.py | autobuild_output has same score as refined (100) |
| `tst_stop_cycle_excluded_from_count` | tst_session_summary.py | STOP cycles not counted in total_cycles |
| `tst_cryoem_done_flags` | tst_workflow_state.py | cryo-EM done flags (resolve_cryo_em, map_sharpening, etc.) |
| `tst_failed_refine_not_counted` | tst_workflow_state.py | Failed refinements don't increment refine_count |
| `tst_mr_sad_after_phaser_with_anomalous` | tst_workflow_state.py | MR-SAD state after phaser + anomalous data |
| `tst_mr_sad_not_triggered_without_anomalous` | tst_workflow_state.py | MR-SAD not triggered without anomalous signal |
| `tst_mr_sad_not_triggered_when_autosol_done` | tst_workflow_state.py | MR-SAD skipped if autosol already ran |
| `tst_mr_sad_directive_prioritizes_phaser` | tst_workflow_state.py | use_mr_sad prioritizes phaser, removes autosol from obtain_model |
| `tst_normal_sad_still_works` | tst_workflow_state.py | Normal SAD workflow unaffected by MR-SAD changes |
| `tst_autosol_has_partial_model_config` | tst_workflow_state.py | autosol YAML has partpdb_file optional input |
| `tst_mr_sad_directive_overrides_no_anomalous` | tst_workflow_state.py | use_mr_sad triggers MR-SAD even without has_anomalous |
| `tst_experimental_phasing_yaml_structure` | tst_workflow_state.py | workflows.yaml experimental_phasing phase structure |
| `tst_predict_and_build_blocked_after_full_completion` | tst_workflow_state.py | predict_and_build not re-offered after full run, even with directives |
| `tst_autosol_partpdb_file_in_command` | tst_command_builder.py | autosol MR-SAD command includes partpdb_file=PHASER.1.pdb |
| `tst_model_scoring_prefers_refined` | tst_best_files_tracker.py | Refined models scored higher than MR output |
| `tst_predicted_model_exclusion` | tst_best_files_tracker.py | Predicted models excluded from model category |

## Writing New Tests

### 1. Create Test File

```bash
touch tests/tst_my_feature.py
```

### 2. Add Boilerplate

```python
"""
Tests for my_feature module.

Run with: python tests/tst_my_feature.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_in,
    run_tests_with_fail_fast
)

from agent.my_feature import my_function


# =============================================================================
# MY FEATURE TESTS
# =============================================================================

def test_my_function_basic():
    """Test basic functionality."""
    result = my_function("input")
    assert_equal(result, "expected")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
```

### 3. Run Tests

```bash
python tests/tst_my_feature.py
```

## Handling Optional Dependencies

Some tests require PHENIX/libtbx. Use runtime checks to skip gracefully:

```python
# Try to import - may fail without libtbx
try:
    from agent.workflow_engine import WorkflowEngine
    HAS_WORKFLOW_ENGINE = True
except ImportError:
    HAS_WORKFLOW_ENGINE = False
    WorkflowEngine = None


def test_workflow_engine_feature():
    """Test that requires workflow engine."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return
    
    engine = WorkflowEngine()
    result = engine.some_method()
    assert_equal(result, expected)
```

## Test Categories

### Unit Tests (No External Dependencies)

These test individual functions in isolation:
- `tst_directive_extractor.py`
- `tst_transport.py`
- `tst_best_files_tracker.py`

### Integration Tests (May Need PHENIX)

These test component interactions:
- `tst_directives_integration.py`
- `tst_decision_flow.py`
- `tst_session_directives.py`

### Scenario Tests (Full Workflows)

Located in `tests/scenarios/`, these test complete workflows:
- `dry_run_xray_basic.py`
- `dry_run_cryoem_basic.py`

## Continuous Integration

The test suite is designed for CI environments:

```bash
# Quick CI run (no PHENIX required)
python tests/run_all_tests.py --quick

# Full CI run (with PHENIX)
python tests/run_all_tests.py
```

Exit codes:
- `0`: All tests passed
- `1`: One or more tests failed

## Debugging Test Failures

When a test fails, you get a full traceback:

```
  test_something ... FAIL
Traceback (most recent call last):
  File "tests/tst_module.py", line 42, in <module>
    run_all_tests()
  File "tests/tst_utils.py", line 357, in run_tests_with_fail_fast
    func()
  File "tests/tst_module.py", line 25, in test_something
    assert_equal(result, expected)
  File "tests/tst_utils.py", line 47, in assert_equal
    raise AssertionError(msg)
AssertionError: expected 'foo', got 'bar'
```

To debug interactively:

```python
# Add this before the failing assertion
import pdb; pdb.set_trace()
```

## Best Practices

1. **One assertion focus per test** - Test one specific behavior
2. **Descriptive names** - `test_extract_resolution_from_angstrom_notation`
3. **Clear docstrings** - Explain what the test verifies
4. **Minimal setup** - Keep tests fast and independent
5. **Section headers** - Group related tests with `# ===` comments
6. **Edge cases** - Test None, empty, boundary conditions

## Migration from unittest

The codebase previously used `unittest.TestCase`. If you encounter old-style tests:

**Before (unittest style):**
```python
class TestSomething(unittest.TestCase):
    def setUp(self):
        self.obj = create_object()
    
    def test_feature(self):
        self.assertEqual(self.obj.method(), "expected")
```

**After (cctbx style):**
```python
def test_feature():
    """Test the feature."""
    obj = create_object()
    assert_equal(obj.method(), "expected")
```

The test runner supports both styles during migration, but new tests should use plain functions.
