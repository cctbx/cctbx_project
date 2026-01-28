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
python tests/test_directive_extractor.py

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

All tests use helper functions from `tests/test_utils.py` instead of `unittest.TestCase` assertions:

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

Run with: python tests/test_module_name.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
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
| `test_directive_extractor.py` | 61 | Directive extraction and validation |
| `test_directive_validator.py` | 48 | Intent validation and modification |
| `test_transport.py` | 66 | REST transport sanitization |
| `test_new_programs.py` | 44 | YAML config for new programs |
| `test_error_analyzer.py` | 34 | Error recovery system |
| `test_decision_flow.py` | 20 | Directive flow architecture |
| `test_directives_integration.py` | 16 | End-to-end directive tests |
| `test_session_directives.py` | 12 | Session-level directive handling |
| `test_best_files_tracker.py` | 45 | File tracking and selection |
| `test_file_categorization.py` | ~30 | File type detection |
| `test_event_system.py` | ~40 | Event logging system |
| `test_yaml_tools.py` | ~20 | YAML validation |
| `test_metric_patterns.py` | ~50 | Log parsing patterns |
| `test_workflow_state.py` | ~60 | Workflow state detection and done flags |
| `test_session_summary.py` | ~20 | Session summary generation |

Total: **350+ tests**

## Writing New Tests

### 1. Create Test File

```bash
touch tests/test_my_feature.py
```

### 2. Add Boilerplate

```python
"""
Tests for my_feature module.

Run with: python tests/test_my_feature.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
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
python tests/test_my_feature.py
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
- `test_directive_extractor.py`
- `test_transport.py`
- `test_best_files_tracker.py`

### Integration Tests (May Need PHENIX)

These test component interactions:
- `test_directives_integration.py`
- `test_decision_flow.py`
- `test_session_directives.py`

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
  File "tests/test_module.py", line 42, in <module>
    run_all_tests()
  File "tests/test_utils.py", line 357, in run_tests_with_fail_fast
    func()
  File "tests/test_module.py", line 25, in test_something
    assert_equal(result, expected)
  File "tests/test_utils.py", line 47, in assert_equal
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
