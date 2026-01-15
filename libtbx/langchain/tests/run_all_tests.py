#!/usr/bin/env python
"""
Master test runner for PHENIX AI Agent v2.

Runs all test suites:
1. Metrics analyzer tests
2. Workflow state tests
3. Integration tests
4. YAML configuration tests
5. Dry run manager tests
6. Sanity checker tests

Usage:
    python run_all_tests.py

Or run individual test files:
    python tests/test_metrics_analyzer.py
    python tests/test_workflow_state.py
    python tests/test_integration.py
    python tests/test_sanity_checker.py
"""

from __future__ import absolute_import, division, print_function
import sys
import traceback


def run_test_suite(name, test_module):
    """Run a test suite and report results."""
    print("\n" + "=" * 70)
    print("RUNNING: %s" % name)
    print("=" * 70 + "\n")

    try:
        test_module.run_all_tests()
        return True
    except Exception as e:
        print("\n!!! TEST SUITE FAILED !!!")
        print("Error: %s" % str(e))
        traceback.print_exc()
        return False


def main():
    """Run all test suites."""
    print("=" * 70)
    print("PHENIX AI AGENT V2 - TEST SUITE")
    print("=" * 70)

    results = {}

    # Test 1: Metrics Analyzer
    from libtbx.langchain.tests import test_metrics_analyzer
    results["Metrics Analyzer"] = run_test_suite("Metrics Analyzer Tests", test_metrics_analyzer)

    # Test 2: Workflow State
    from libtbx.langchain.tests import test_workflow_state
    results["Workflow State"] = run_test_suite("Workflow State Tests", test_workflow_state)

    # Test 3: Integration
    from libtbx.langchain.tests import test_integration
    results["Integration"] = run_test_suite("Integration Tests", test_integration)

    # Test 4: YAML Configuration
    from libtbx.langchain.tests import test_yaml_config
    results["YAML Config"] = run_test_suite("YAML Configuration Tests", test_yaml_config)

    # Test 5: Dry Run Manager
    from libtbx.langchain.tests import test_dry_run
    results["Dry Run"] = run_test_suite("Dry Run Manager Tests", test_dry_run)

    # Test 6: Sanity Checker
    from libtbx.langchain.tests import test_sanity_checker
    results["Sanity Checker"] = run_test_suite("Sanity Checker Tests", test_sanity_checker)

    # Summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    all_passed = True
    for suite, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print("  %s: %s" % (suite, status))
        if not passed:
            all_passed = False

    print("=" * 70)

    if all_passed:
        print("\n*** ALL TESTS PASSED ***\n")
        return 0
    else:
        print("\n*** SOME TESTS FAILED ***\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
