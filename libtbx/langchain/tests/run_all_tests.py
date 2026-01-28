#!/usr/bin/env python
"""
Run all tests for the PHENIX AI Agent.

Tests use cctbx-style fail-fast behavior: the first assertion failure
raises an uncaught exception with full traceback, rather than collecting
all failures. This makes it easy to identify and debug issues.

Test Suites (standalone - no PHENIX required):
  1. API Schema - Request/response validation
  2. Best Files Tracker - File tracking and scoring
  3. Transport - Encoding/decoding round-trips
  4. State Serialization - State packaging/unpackaging
  5. Command Builder - Unified command generation
  6. File Categorization - File classification logic
  7. Session Summary - Agent session summary generation
  8. Advice Preprocessing - README discovery, advice processing, change detection
  9. Directive Extractor - LLM-based directive extraction from user advice
  10. Directive Validator - Pre-validation of user requests against capabilities
  11. Session Directives - Session directive storage and retrieval
  12. YAML Tools - YAML configuration validation and inspection
  13. Session Tools - Session management utilities
  14. Docs Tools - Documentation generation
  15. Error Analyzer - Error detection and recovery strategies
  16. Decision Flow - Decision flow logic testing
  17. Phaser Multimodel - Phaser multi-model handling
  18. Event System - Event logging and tracking
  19. Metric Patterns - Metric extraction patterns
  20. Pattern Manager - Pattern management
  21. Program Registration - Program registry tests
  22. Summary Display - Summary formatting

Test Suites (require PHENIX environment):
  23. Workflow State - State detection and transitions
  24. YAML Config - YAML configuration validation
  25. Sanity Checker - Sanity check logic
  26. Metrics Analyzer - Metric extraction and trends
  27. Dry Run - Dry run manager functionality
  28. Integration - End-to-end workflow tests
  29. Directives Integration - End-to-end directive system tests

Usage:
    python tests/run_all_tests.py
    python tests/run_all_tests.py --verbose
    python tests/run_all_tests.py --quick  # Skip slow integration tests

Note: Tests requiring PHENIX will be skipped with a warning if libtbx is not available.
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

    # --- Transport Tests ---
    try:
        from tests.test_transport import run_all_tests as run_transport_tests
        success, elapsed = run_test_module(
            "test_transport", run_transport_tests, args.verbose)
        results.append(("Transport", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_transport: {e}")
        results.append(("Transport", False, 0))

    # --- State Serialization Tests ---
    try:
        from tests.test_state_serialization import run_all_tests as run_serialization_tests
        success, elapsed = run_test_module(
            "test_state_serialization", run_serialization_tests, args.verbose)
        results.append(("State Serialization", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_state_serialization: {e}")
        results.append(("State Serialization", False, 0))

    # --- Command Builder Tests ---
    try:
        from tests.test_command_builder import run_all_tests as run_command_builder_tests
        success, elapsed = run_test_module(
            "test_command_builder", run_command_builder_tests, args.verbose)
        results.append(("Command Builder", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_command_builder: {e}")
        results.append(("Command Builder", False, 0))

    # --- File Categorization Tests ---
    try:
        from tests.test_file_categorization import run_all_tests as run_file_categorization_tests
        success, elapsed = run_test_module(
            "test_file_categorization", run_file_categorization_tests, args.verbose)
        results.append(("File Categorization", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_file_categorization: {e}")
        results.append(("File Categorization", False, 0))

    # --- Session Summary Tests ---
    try:
        from tests.test_session_summary import run_all_tests as run_session_summary_tests
        success, elapsed = run_test_module(
            "test_session_summary", run_session_summary_tests, args.verbose)
        results.append(("Session Summary", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_session_summary: {e}")
        results.append(("Session Summary", False, 0))

    # --- Advice Preprocessing Tests ---
    try:
        from tests.test_advice_preprocessing import run_all_tests as run_advice_preprocessing_tests
        success, elapsed = run_test_module(
            "test_advice_preprocessing", run_advice_preprocessing_tests, args.verbose)
        results.append(("Advice Preprocessing", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_advice_preprocessing: {e}")
        results.append(("Advice Preprocessing", False, 0))

    # --- Directive Extractor Tests ---
    try:
        from tests.test_directive_extractor import run_all_tests as run_directive_extractor_tests
        success, elapsed = run_test_module(
            "test_directive_extractor", run_directive_extractor_tests, args.verbose)
        results.append(("Directive Extractor", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_directive_extractor: {e}")
        results.append(("Directive Extractor", False, 0))

    # --- Directive Validator Tests ---
    try:
        from tests.test_directive_validator import run_all_tests as run_directive_validator_tests
        success, elapsed = run_test_module(
            "test_directive_validator", run_directive_validator_tests, args.verbose)
        results.append(("Directive Validator", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_directive_validator: {e}")
        results.append(("Directive Validator", False, 0))

    # --- Session Directives Tests ---
    try:
        from tests.test_session_directives import run_all_tests as run_session_directives_tests
        success, elapsed = run_test_module(
            "test_session_directives", run_session_directives_tests, args.verbose)
        results.append(("Session Directives", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_session_directives: {e}")
        results.append(("Session Directives", False, 0))

    # --- YAML Tools Tests ---
    try:
        from tests.test_yaml_tools import run_all_tests as run_yaml_tools_tests
        success, elapsed = run_test_module(
            "test_yaml_tools", run_yaml_tools_tests, args.verbose)
        results.append(("YAML Tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_yaml_tools: {e}")
        results.append(("YAML Tools", False, 0))

    # --- Session Tools Tests ---
    try:
        from tests.test_session_tools import run_all_tests as run_session_tools_tests
        success, elapsed = run_test_module(
            "test_session_tools", run_session_tools_tests, args.verbose)
        results.append(("Session Tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_session_tools: {e}")
        results.append(("Session Tools", False, 0))

    # --- Docs Tools Tests ---
    try:
        from tests.test_docs_tools import run_all_tests as run_docs_tools_tests
        success, elapsed = run_test_module(
            "test_docs_tools", run_docs_tools_tests, args.verbose)
        results.append(("Docs Tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_docs_tools: {e}")
        results.append(("Docs Tools", False, 0))

    # --- Metric Patterns Tests ---
    try:
        from tests.test_metric_patterns import run_all_tests as run_metric_patterns_tests
        success, elapsed = run_test_module(
            "test_metric_patterns", run_metric_patterns_tests, args.verbose)
        results.append(("Metric Patterns", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_metric_patterns: {e}")
        results.append(("Metric Patterns", False, 0))

    # --- Pattern Manager Tests ---
    try:
        from tests.test_pattern_manager import run_all_tests as run_pattern_manager_tests
        success, elapsed = run_test_module(
            "test_pattern_manager", run_pattern_manager_tests, args.verbose)
        results.append(("Pattern Manager", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_pattern_manager: {e}")
        results.append(("Pattern Manager", False, 0))

    # --- Event System Tests ---
    try:
        from tests.test_event_system import run_all_tests as run_event_system_tests
        success, elapsed = run_test_module(
            "test_event_system", run_event_system_tests, args.verbose)
        results.append(("Event System", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_event_system: {e}")
        results.append(("Event System", False, 0))

    # --- Program Registration Tests ---
    try:
        from tests.test_program_registration import run_all_tests as run_program_registration_tests
        success, elapsed = run_test_module(
            "test_program_registration", run_program_registration_tests, args.verbose)
        results.append(("Program Registration", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_program_registration: {e}")
        results.append(("Program Registration", False, 0))

    # --- Summary Display Tests ---
    try:
        from tests.test_summary_display import run_all_tests as run_summary_display_tests
        success, elapsed = run_test_module(
            "test_summary_display", run_summary_display_tests, args.verbose)
        results.append(("Summary Display", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_summary_display: {e}")
        results.append(("Summary Display", False, 0))

    # --- Directives Integration Tests (can be slow) ---
    if not args.quick:
        try:
            from tests.test_directives_integration import run_all_tests as run_directives_integration_tests
            success, elapsed = run_test_module(
                "test_directives_integration", run_directives_integration_tests, args.verbose)
            results.append(("Directives Integration", success, elapsed))
        except ImportError as e:
            print(f"⚠️  Could not import test_directives_integration: {e}")
            results.append(("Directives Integration", False, 0))
    else:
        print("\n⏭️  Skipping directives integration tests (--quick mode)")

    # --- Error Analyzer Tests ---
    try:
        from tests.test_error_analyzer import run_all_tests as run_error_analyzer_tests
        success, elapsed = run_test_module(
            "test_error_analyzer", run_error_analyzer_tests, args.verbose)
        results.append(("Error Analyzer", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_error_analyzer: {e}")
        results.append(("Error Analyzer", False, 0))

    # --- Decision Flow Tests ---
    try:
        from tests.test_decision_flow import run_all_tests as run_decision_flow_tests
        success, elapsed = run_test_module(
            "test_decision_flow", run_decision_flow_tests, args.verbose)
        results.append(("Decision Flow", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_decision_flow: {e}")
        results.append(("Decision Flow", False, 0))

    # --- Phaser Multimodel Tests ---
    try:
        from tests.test_phaser_multimodel import run_all_tests as run_phaser_multimodel_tests
        success, elapsed = run_test_module(
            "test_phaser_multimodel", run_phaser_multimodel_tests, args.verbose)
        results.append(("Phaser Multimodel", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_phaser_multimodel: {e}")
        results.append(("Phaser Multimodel", False, 0))

    # --- New Programs Tests (polder, map_sharpening, etc.) ---
    try:
        from tests.test_new_programs import run_all_tests as run_new_programs_tests
        success, elapsed = run_test_module(
            "test_new_programs", run_new_programs_tests, args.verbose)
        results.append(("New Programs", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import test_new_programs: {e}")
        results.append(("New Programs", False, 0))

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
