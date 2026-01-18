#!/usr/bin/env python
"""
Run all tests for the PHENIX AI Agent.

Usage:
    python tests/run_all_tests.py
    python tests/run_all_tests.py --verbose
    python tests/run_all_tests.py --quick  # Skip slow integration tests
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import time
import argparse
import traceback

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_test_module(module_name, run_func, verbose=False):
    """
    Run a test module and return success status.

    Args:
        module_name: Name of the test module
        run_func: Function to call to run tests
        verbose: If True, show full output

    Returns:
        tuple: (success, elapsed_time)
    """
    print(f"\n{'='*60}")
    print(f"Running: {module_name}")
    print('='*60)

    start_time = time.time()

    try:
        run_func()
        elapsed = time.time() - start_time
        print(f"\n✅ {module_name} PASSED ({elapsed:.2f}s)")
        return True, elapsed
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n❌ {module_name} FAILED ({elapsed:.2f}s)")
        if verbose:
            traceback.print_exc()
        else:
            print(f"   Error: {e}")
        return False, elapsed


def main():
    parser = argparse.ArgumentParser(description="Run all AI Agent tests")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show verbose output including tracebacks")
    parser.add_argument("--quick", "-q", action="store_true",
                        help="Skip slow integration tests")
    args = parser.parse_args()

    print("="*60)
    print("PHENIX AI AGENT - TEST SUITE")
    print("="*60)

    results = []
    total_start = time.time()

    # --- API Schema Tests ---
    try:
        from tests.test_api_schema import run_all_tests as run_api_schema_tests
        success, elapsed = run_test_module(
            "test_api_schema", run_api_schema_tests, args.verbose)
        results.append(("API Schema", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_api_schema: {e}")
        results.append(("API Schema", False, 0))

    # --- Best Files Tracker Tests ---
    try:
        from tests.test_best_files_tracker import run_all_tests as run_best_files_tests
        success, elapsed = run_test_module(
            "test_best_files_tracker", run_best_files_tests, args.verbose)
        results.append(("Best Files Tracker", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_best_files_tracker: {e}")
        results.append(("Best Files Tracker", False, 0))

    # --- Workflow State Tests ---
    try:
        from tests.test_workflow_state import run_all_tests as run_workflow_tests
        success, elapsed = run_test_module(
            "test_workflow_state", run_workflow_tests, args.verbose)
        results.append(("Workflow State", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_workflow_state: {e}")
        results.append(("Workflow State", False, 0))

    # --- YAML Config Tests ---
    try:
        from tests.test_yaml_config import run_all_tests as run_yaml_tests
        success, elapsed = run_test_module(
            "test_yaml_config", run_yaml_tests, args.verbose)
        results.append(("YAML Config", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_yaml_config: {e}")
        results.append(("YAML Config", False, 0))

    # --- Sanity Checker Tests ---
    try:
        from tests.test_sanity_checker import run_all_tests as run_sanity_tests
        success, elapsed = run_test_module(
            "test_sanity_checker", run_sanity_tests, args.verbose)
        results.append(("Sanity Checker", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_sanity_checker: {e}")
        results.append(("Sanity Checker", False, 0))

    # --- Metrics Analyzer Tests ---
    try:
        from tests.test_metrics_analyzer import run_all_tests as run_metrics_tests
        success, elapsed = run_test_module(
            "test_metrics_analyzer", run_metrics_tests, args.verbose)
        results.append(("Metrics Analyzer", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_metrics_analyzer: {e}")
        results.append(("Metrics Analyzer", False, 0))

    # --- Dry Run Tests ---
    try:
        from tests.test_dry_run import run_all_tests as run_dry_run_tests
        success, elapsed = run_test_module(
            "test_dry_run", run_dry_run_tests, args.verbose)
        results.append(("Dry Run", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_dry_run: {e}")
        results.append(("Dry Run", False, 0))

    # --- Integration Tests (slow, skip with --quick) ---
    if not args.quick:
        try:
            from tests.test_integration import run_all_tests as run_integration_tests
            success, elapsed = run_test_module(
                "test_integration", run_integration_tests, args.verbose)
            results.append(("Integration", success, elapsed))
        except ImportError as e:
            print(f"⚠️  Could not import test_integration: {e}")
            results.append(("Integration", False, 0))
    else:
        print("\n⏭️  Skipping integration tests (--quick mode)")

    # --- Summary ---
    total_elapsed = time.time() - total_start

    print("\n")
    print("="*60)
    print("TEST SUMMARY")
    print("="*60)

    passed = 0
    failed = 0

    for name, success, elapsed in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"  {name:25s} {status} ({elapsed:.2f}s)")
        if success:
            passed += 1
        else:
            failed += 1

    print("-"*60)
    print(f"  Total: {passed} passed, {failed} failed ({total_elapsed:.2f}s)")
    print("="*60)

    if failed > 0:
        print("\n❌ SOME TESTS FAILED")
        sys.exit(1)
    else:
        print("\n✅ ALL TESTS PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
