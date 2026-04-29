#!/usr/bin/env python
"""
Run all tests for the PHENIX AI Agent.

Tests use cctbx-style fail-fast behavior: the first assertion failure
raises an uncaught exception with full traceback, rather than collecting
all failures. This makes it easy to identify and debug issues.

Test Suites (standalone - no PHENIX required):
  1. API Schema - Request/response validation
  2. Best Files Tracker - File tracking, scoring, and model stage detection
  3. Transport - Encoding/decoding round-trips
  4. State Serialization - State packaging/unpackaging
  5. Command Builder - Unified command generation
  6. File Categorization - File classification logic
  7. Session Summary - Agent session summary generation (cycle counting)
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
  23. New Programs - YAML config for new programs (polder, map_sharpening, autobuild_denmod)
  24. History Analysis - History analysis, anomalous workflow support (v110)
  25. File Utils - Shared file classification utilities (v110)

Test Suites (require PHENIX environment):
  26. Workflow State - State detection, transitions, done flags, stepwise mode
  27. YAML Config - YAML configuration validation
  28. Sanity Checker - Sanity check logic
  29. Metrics Analyzer - Metric extraction and trends
  30. Dry Run - Dry run manager functionality
  31. Integration - End-to-end workflow tests
  32. Directives Integration - End-to-end directive system tests

  33. Langchain Tools - Legacy module tests (core, analysis, RAG, validation, prompts, memory)
  34. Hardcoded Cleanup - Conformance guards for YAML-driven architecture (v112.10)
  35. v112.13        - Companion files, intermediate filtering, file categorisation
  36. Audit Fix Regressions - Categories I/J/E/G/H: max_refine_cycles landing, zombie
      state detection, _is_failed_result false-positives, xtriage resolution regex,
      real_space_refine map_cc extract strategy (v112 systematic audit)
  37. Autosol Bugs - Wavelength duplication, atom_type swap, autosol re-run,
      heavier-atom-wins, catch-all injection blacklist, rebuild_in_place
      allowlist (v112.75–v112.77)
  38. Thinking Defense - Pre-implementation baseline for Thinking Agent (v113)
  39. Contract Compliance - Contract compliance validation
  40. Explanation Prompts - LLM explanation prompt generation
  41. Gate Evaluator - Gate evaluation logic (advance/skip/retreat/stop)
  42. Hypothesis Evaluator - Hypothesis lifecycle, ligand/metal evaluation
  43. MTZ Crosscheck - MTZ file classification and crosscheck rescue
  44. Old Client Compat - Backward compatibility with older clients
  45. Perceive Stop Checks - Perceive-step stop condition checks
  46. Phase→Stage Rename Audit - Runtime introspection audit
  47. Plan Generator - Plan generation, revision, and repair logic
  48. Plan Schema - Plan schema validation, roundtrips, workflow tests
  49. Shared Code Imports - Shared code import verification
  50. Structure Model - Structure model state, metrics, problem detection
  51. Template - Template system tests
  52. Validation History - Validation history roundtrips, scenarios
  53. Display Data Model - DisplayDataModel properties, HTML report generation
  54. Scenario Tracer - Plan selection, gate evaluation, cycle counting, multi-cycle progression
  55. v115 P1-P5   - v115 bug fixes: session blocking, autobuild prerequisites, MTZ phase
      detection, cryo-EM phase 1.5 gate, solve-mode file discovery; plus Bug A–F fixes:
      event_formatter string metrics, file discovery supplement mode, predict_and_build
      internal refine PDB exclusion, autosol anomalous deprioritization, CASP7 search
      model misclassification, hopeless R-free MR retry, weak anomalous signal detection
      (v115)
  56. v115 N-Bugs    - N1 event_formatter TypeError, N3 file_categories refined exclude
  57. v115 Fix 2     - fabricated stop-condition lines stripped before directive extraction
  58. v115 F1+F3     - abort_message for all stops; max-cycles explicit message
  59. v115 F2+F4+F7  - API key errors; blocked-program naming; empty command stop
  60. v115 F5        - workflow_complete stop reason on rules-only success
  61. v115 F6        - diagnosed_failure stop reason for terminal diagnosis
  62. v115 Follow-up - abort_message for no_workflow_state and after_program_not_available
  63. P1B            - STOP_REASON_CODES, think_stop_override, think_file_overrides
  64. P1B Prompt     - stop_reason table injection, guard, graceful degradation
  65. P1B Action 1   - validate() think_file_overrides existence check

Systematic Testing Framework (v115.08):
  S0. Phase 0: Static Audit - automated gate (parse, bare except, import fallbacks)
  S1. Phase 1: Contract Gaps - coverage map of 128 functions across 4 key modules
  S2. Phase 2: Path Consistency - YAML vs hardcoded categorization (10 tutorials)
  S3. Phase 3: Session Round-Trip - JSON symmetry + invariants + AgentSession (28 tests)
  S4. Phase 4: History Flags - flag writer/reader consistency, 12 dead flags documented
  S5. Phase 6: Category-Consumer - input_priorities + fallback_categories (14 checks)
  S6. Phase 7: Routing Simulation - 32 tutorials x 3 cycles, 9 expected stuck whitelisted
  S7. Phase 8: Command Building - 7 tutorials + 2 edge cases, 15 command builds
  S8. Phase 5: Error Classification - 3 classifiers, 30+ patterns, consistency checks
  S9. Phase 9: LLM Perturbation - filename/program/parameter/truncation/empty (17 tests)
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
        from tests.tst_api_schema import run_all_tests as run_api_schema_tests
        success, elapsed = run_test_module(
            "tst_api_schema", run_api_schema_tests, args.verbose)
        results.append(("API Schema", "tst_api_schema", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_api_schema: {e}")
        results.append(("API Schema", "tst_api_schema", False, 0))

    # --- Best Files Tracker Tests ---
    try:
        from tests.tst_best_files_tracker import run_all_tests as run_best_files_tests
        success, elapsed = run_test_module(
            "tst_best_files_tracker", run_best_files_tests, args.verbose)
        results.append(("Best Files Tracker", "tst_best_files_tracker", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_best_files_tracker: {e}")
        results.append(("Best Files Tracker", "tst_best_files_tracker", False, 0))

    # --- Workflow State Tests ---
    try:
        from tests.tst_workflow_state import run_all_tests as run_workflow_tests
        success, elapsed = run_test_module(
            "tst_workflow_state", run_workflow_tests, args.verbose)
        results.append(("Workflow State", "tst_workflow_state", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_workflow_state: {e}")
        results.append(("Workflow State", "tst_workflow_state", False, 0))

    # --- YAML Config Tests ---
    try:
        from tests.tst_yaml_config import run_all_tests as run_yaml_tests
        success, elapsed = run_test_module(
            "tst_yaml_config", run_yaml_tests, args.verbose)
        results.append(("YAML Config", "tst_yaml_config", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_yaml_config: {e}")
        results.append(("YAML Config", "tst_yaml_config", False, 0))

    # --- Sanity Checker Tests ---
    try:
        from tests.tst_sanity_checker import run_all_tests as run_sanity_tests
        success, elapsed = run_test_module(
            "tst_sanity_checker", run_sanity_tests, args.verbose)
        results.append(("Sanity Checker", "tst_sanity_checker", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_sanity_checker: {e}")
        results.append(("Sanity Checker", "tst_sanity_checker", False, 0))

    # --- Metrics Analyzer Tests ---
    try:
        from tests.tst_metrics_analyzer import run_all_tests as run_metrics_tests
        success, elapsed = run_test_module(
            "tst_metrics_analyzer", run_metrics_tests, args.verbose)
        results.append(("Metrics Analyzer", "tst_metrics_analyzer", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_metrics_analyzer: {e}")
        results.append(("Metrics Analyzer", "tst_metrics_analyzer", False, 0))

    # --- Dry Run Tests ---
    try:
        from tests.tst_dry_run import run_all_tests as run_dry_run_tests
        success, elapsed = run_test_module(
            "tst_dry_run", run_dry_run_tests, args.verbose)
        results.append(("Dry Run", "tst_dry_run", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_dry_run: {e}")
        results.append(("Dry Run", "tst_dry_run", False, 0))

    # --- Integration Tests (slow, skip with --quick) ---
    if not args.quick:
        try:
            from tests.tst_integration import run_all_tests as run_integration_tests
            success, elapsed = run_test_module(
                "tst_integration", run_integration_tests, args.verbose)
            results.append(("Integration", "tst_integration", success, elapsed))
        except ImportError as e:
            print(f"⚠️  Could not import tst_integration: {e}")
            results.append(("Integration", "tst_integration", False, 0))
    else:
        print("\n⏭️  Skipping integration tests (--quick mode)")

    # --- Transport Tests ---
    try:
        from tests.tst_transport import run_all_tests as run_transport_tests
        success, elapsed = run_test_module(
            "tst_transport", run_transport_tests, args.verbose)
        results.append(("Transport", "tst_transport", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_transport: {e}")
        results.append(("Transport", "tst_transport", False, 0))

    # --- State Serialization Tests ---
    try:
        from tests.tst_state_serialization import run_all_tests as run_serialization_tests
        success, elapsed = run_test_module(
            "tst_state_serialization", run_serialization_tests, args.verbose)
        results.append(("State Serialization", "tst_state_serialization", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_state_serialization: {e}")
        results.append(("State Serialization", "tst_state_serialization", False, 0))

    # --- Command Builder Tests ---
    try:
        from tests.tst_command_builder import run_all_tests as run_command_builder_tests
        success, elapsed = run_test_module(
            "tst_command_builder", run_command_builder_tests, args.verbose)
        results.append(("Command Builder", "tst_command_builder", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_command_builder: {e}")
        results.append(("Command Builder", "tst_command_builder", False, 0))

    # --- File Categorization Tests ---
    try:
        from tests.tst_file_categorization import run_all_tests as run_file_categorization_tests
        success, elapsed = run_test_module(
            "tst_file_categorization", run_file_categorization_tests, args.verbose)
        results.append(("File Categorization", "tst_file_categorization", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_file_categorization: {e}")
        results.append(("File Categorization", "tst_file_categorization", False, 0))

    # --- Session Summary Tests ---
    try:
        from tests.tst_session_summary import run_all_tests as run_session_summary_tests
        success, elapsed = run_test_module(
            "tst_session_summary", run_session_summary_tests, args.verbose)
        results.append(("Session Summary", "tst_session_summary", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_session_summary: {e}")
        results.append(("Session Summary", "tst_session_summary", False, 0))

    # --- Advice Preprocessing Tests ---
    try:
        from tests.tst_advice_preprocessing import run_all_tests as run_advice_preprocessing_tests
        success, elapsed = run_test_module(
            "tst_advice_preprocessing", run_advice_preprocessing_tests, args.verbose)
        results.append(("Advice Preprocessing", "tst_advice_preprocessing", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_advice_preprocessing: {e}")
        results.append(("Advice Preprocessing", "tst_advice_preprocessing", False, 0))

    # --- Directive Extractor Tests ---
    try:
        from tests.tst_directive_extractor import run_all_tests as run_directive_extractor_tests
        success, elapsed = run_test_module(
            "tst_directive_extractor", run_directive_extractor_tests, args.verbose)
        results.append(("Directive Extractor", "tst_directive_extractor", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_directive_extractor: {e}")
        results.append(("Directive Extractor", "tst_directive_extractor", False, 0))

    # --- Directive Validator Tests ---
    try:
        from tests.tst_directive_validator import run_all_tests as run_directive_validator_tests
        success, elapsed = run_test_module(
            "tst_directive_validator", run_directive_validator_tests, args.verbose)
        results.append(("Directive Validator", "tst_directive_validator", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_directive_validator: {e}")
        results.append(("Directive Validator", "tst_directive_validator", False, 0))

    # --- Session Directives Tests ---
    try:
        from tests.tst_session_directives import run_all_tests as run_session_directives_tests
        success, elapsed = run_test_module(
            "tst_session_directives", run_session_directives_tests, args.verbose)
        results.append(("Session Directives", "tst_session_directives", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_session_directives: {e}")
        results.append(("Session Directives", "tst_session_directives", False, 0))

    # --- YAML Tools Tests ---
    try:
        from tests.tst_yaml_tools import run_all_tests as run_yaml_tools_tests
        success, elapsed = run_test_module(
            "tst_yaml_tools", run_yaml_tools_tests, args.verbose)
        results.append(("YAML Tools", "tst_yaml_tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_yaml_tools: {e}")
        results.append(("YAML Tools", "tst_yaml_tools", False, 0))

    # --- Session Tools Tests ---
    try:
        from tests.tst_session_tools import run_all_tests as run_session_tools_tests
        success, elapsed = run_test_module(
            "tst_session_tools", run_session_tools_tests, args.verbose)
        results.append(("Session Tools", "tst_session_tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_session_tools: {e}")
        results.append(("Session Tools", "tst_session_tools", False, 0))

    # --- Docs Tools Tests ---
    try:
        from tests.tst_docs_tools import run_all_tests as run_docs_tools_tests
        success, elapsed = run_test_module(
            "tst_docs_tools", run_docs_tools_tests, args.verbose)
        results.append(("Docs Tools", "tst_docs_tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_docs_tools: {e}")
        results.append(("Docs Tools", "tst_docs_tools", False, 0))

    # --- Metric Patterns Tests ---
    try:
        from tests.tst_metric_patterns import run_all_tests as run_metric_patterns_tests
        success, elapsed = run_test_module(
            "tst_metric_patterns", run_metric_patterns_tests, args.verbose)
        results.append(("Metric Patterns", "tst_metric_patterns", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_metric_patterns: {e}")
        results.append(("Metric Patterns", "tst_metric_patterns", False, 0))

    # --- Pattern Manager Tests ---
    try:
        from tests.tst_pattern_manager import run_all_tests as run_pattern_manager_tests
        success, elapsed = run_test_module(
            "tst_pattern_manager", run_pattern_manager_tests, args.verbose)
        results.append(("Pattern Manager", "tst_pattern_manager", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_pattern_manager: {e}")
        results.append(("Pattern Manager", "tst_pattern_manager", False, 0))

    # --- Event System Tests ---
    try:
        from tests.tst_event_system import run_all_tests as run_event_system_tests
        success, elapsed = run_test_module(
            "tst_event_system", run_event_system_tests, args.verbose)
        results.append(("Event System", "tst_event_system", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_event_system: {e}")
        results.append(("Event System", "tst_event_system", False, 0))

    # --- Program Registration Tests ---
    try:
        from tests.tst_program_registration import run_all_tests as run_program_registration_tests
        success, elapsed = run_test_module(
            "tst_program_registration", run_program_registration_tests, args.verbose)
        results.append(("Program Registration", "tst_program_registration", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_program_registration: {e}")
        results.append(("Program Registration", "tst_program_registration", False, 0))

    # --- Summary Display Tests ---
    try:
        from tests.tst_summary_display import run_all_tests as run_summary_display_tests
        success, elapsed = run_test_module(
            "tst_summary_display", run_summary_display_tests, args.verbose)
        results.append(("Summary Display", "tst_summary_display", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_summary_display: {e}")
        results.append(("Summary Display", "tst_summary_display", False, 0))

    # --- Directives Integration Tests (can be slow) ---
    if not args.quick:
        try:
            from tests.tst_directives_integration import run_all_tests as run_directives_integration_tests
            success, elapsed = run_test_module(
                "tst_directives_integration", run_directives_integration_tests, args.verbose)
            results.append(("Directives Integration", "tst_directives_integration", success, elapsed))
        except ImportError as e:
            print(f"⚠️  Could not import tst_directives_integration: {e}")
            results.append(("Directives Integration", "tst_directives_integration", False, 0))
    else:
        print("\n⏭️  Skipping directives integration tests (--quick mode)")

    # --- Error Analyzer Tests ---
    try:
        from tests.tst_error_analyzer import run_all_tests as run_error_analyzer_tests
        success, elapsed = run_test_module(
            "tst_error_analyzer", run_error_analyzer_tests, args.verbose)
        results.append(("Error Analyzer", "tst_error_analyzer", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_error_analyzer: {e}")
        results.append(("Error Analyzer", "tst_error_analyzer", False, 0))

    # --- Decision Flow Tests ---
    try:
        from tests.tst_decision_flow import run_all_tests as run_decision_flow_tests
        success, elapsed = run_test_module(
            "tst_decision_flow", run_decision_flow_tests, args.verbose)
        results.append(("Decision Flow", "tst_decision_flow", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_decision_flow: {e}")
        results.append(("Decision Flow", "tst_decision_flow", False, 0))

    # --- Phaser Multimodel Tests ---
    try:
        from tests.tst_phaser_multimodel import run_all_tests as run_phaser_multimodel_tests
        success, elapsed = run_test_module(
            "tst_phaser_multimodel", run_phaser_multimodel_tests, args.verbose)
        results.append(("Phaser Multimodel", "tst_phaser_multimodel", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_phaser_multimodel: {e}")
        results.append(("Phaser Multimodel", "tst_phaser_multimodel", False, 0))

    # --- New Programs Tests (polder, map_sharpening, etc.) ---
    try:
        from tests.tst_new_programs import run_all_tests as run_new_programs_tests
        success, elapsed = run_test_module(
            "tst_new_programs", run_new_programs_tests, args.verbose)
        results.append(("New Programs", "tst_new_programs", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_new_programs: {e}")
        results.append(("New Programs", "tst_new_programs", False, 0))

    # --- History Analysis Tests (v110 - anomalous workflow support) ---
    try:
        from tests.tst_history_analysis import run_all_tests as run_history_analysis_tests
        success, elapsed = run_test_module(
            "tst_history_analysis", run_history_analysis_tests, args.verbose)
        results.append(("History Analysis", "tst_history_analysis", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_history_analysis: {e}")
        results.append(("History Analysis", "tst_history_analysis", False, 0))

    # --- File Utils Tests (v110 - shared file classification) ---
    try:
        from tests.tst_file_utils import run_all_tests as run_file_utils_tests
        success, elapsed = run_test_module(
            "tst_file_utils", run_file_utils_tests, args.verbose)
        results.append(("File Utils", "tst_file_utils", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_file_utils: {e}")
        results.append(("File Utils", "tst_file_utils", False, 0))

    # --- Langchain Tools Tests (require PHENIX + langchain_core) ---
    try:
        from tests.tst_langchain_tools import run_all_tests as run_langchain_tools_tests
        success, elapsed = run_test_module(
            "tst_langchain_tools", run_langchain_tools_tests, args.verbose)
        results.append(("Langchain Tools", "tst_langchain_tools", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_langchain_tools: {e}")
        results.append(("Langchain Tools", "tst_langchain_tools", False, 0))

    # --- Hardcoded Cleanup Conformance Tests ---
    try:
        from tests.tst_hardcoded_cleanup import run_all_tests as run_hardcoded_cleanup_tests
        success, elapsed = run_test_module(
            "tst_hardcoded_cleanup", run_hardcoded_cleanup_tests, args.verbose)
        results.append(("Hardcoded Cleanup", "tst_hardcoded_cleanup", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_hardcoded_cleanup: {e}")
        results.append(("Hardcoded Cleanup", "tst_hardcoded_cleanup", False, 0))

    # --- v112.13 Fix Tests ---
    try:
        from tests.tst_v112_13_fixes import run_all_tests as run_v112_13_tests
        success, elapsed = run_test_module(
            "tst_v112_13_fixes", run_v112_13_tests, args.verbose)
        results.append(("v112.13 Fixes", "tst_v112_13_fixes", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_v112_13_fixes: {e}")
        results.append(("v112.13 Fixes", "tst_v112_13_fixes", False, 0))

    # --- Audit Fix Regression Tests (Categories I, J, E, G, H) ---
    try:
        from tests.tst_audit_fixes import run_all_tests as run_audit_fixes_tests
        success, elapsed = run_test_module(
            "tst_audit_fixes", run_audit_fixes_tests, args.verbose)
        results.append(("Audit Fixes", "tst_audit_fixes", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_audit_fixes: {e}")
        results.append(("Audit Fixes", "tst_audit_fixes", False, 0))

    # --- Autosol Bugs Tests (v112.75-v112.77) ---
    try:
        from tests.tst_autosol_bugs import run_all_tests as run_autosol_bugs_tests
        success, elapsed = run_test_module(
            "tst_autosol_bugs", run_autosol_bugs_tests, args.verbose)
        results.append(("Autosol Bugs", "tst_autosol_bugs", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_autosol_bugs: {e}")
        results.append(("Autosol Bugs", "tst_autosol_bugs", False, 0))

    # --- Thinking Defense Tests (v113) ---
    try:
        from tests.tst_thinking_defense import run_all_tests as run_thinking_defense_tests
        success, elapsed = run_test_module(
            "tst_thinking_defense", run_thinking_defense_tests, args.verbose)
        results.append(("Thinking Defense", "tst_thinking_defense", success, elapsed))
    except ImportError as e:
        print(f"Could not import tst_thinking_defense: {e}")
        results.append(("Thinking Defense", "tst_thinking_defense", False, 0))

    # --- Strategy Memory Tests (v113) ---
    try:
        from tests.tst_strategy_memory import run_all_tests as run_strategy_memory_tests
        success, elapsed = run_test_module(
            "tst_strategy_memory", run_strategy_memory_tests, args.verbose)
        results.append(("Strategy Memory", "tst_strategy_memory", success, elapsed))
    except ImportError as e:
        print(f"Could not import tst_strategy_memory: {e}")
        results.append(("Strategy Memory", "tst_strategy_memory", False, 0))

    # --- Log Extractor Tests (v113) ---
    try:
        from tests.tst_log_extractor import run_all_tests as run_log_extractor_tests
        success, elapsed = run_test_module(
            "tst_log_extractor", run_log_extractor_tests, args.verbose)
        results.append(("Log Extractor", "tst_log_extractor", success, elapsed))
    except ImportError as e:
        print(f"Could not import tst_log_extractor: {e}")
        results.append(("Log Extractor", "tst_log_extractor", False, 0))

    # --- Thinking Agent Tests (v113) ---
    try:
        from tests.tst_thinking_agent import run_all_tests as run_thinking_agent_tests
        success, elapsed = run_test_module(
            "tst_thinking_agent", run_thinking_agent_tests, args.verbose)
        results.append(("Thinking Agent", "tst_thinking_agent", success, elapsed))
    except ImportError as e:
        print(f"Could not import tst_thinking_agent: {e}")
        results.append(("Thinking Agent", "tst_thinking_agent", False, 0))

    # --- Contract Compliance Tests ---
    try:
        from tests.tst_contract_compliance import run_all_tests as run_contract_compliance_tests
        success, elapsed = run_test_module(
            "tst_contract_compliance", run_contract_compliance_tests, args.verbose)
        results.append(("Contract Compliance", "tst_contract_compliance", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_contract_compliance: {e}")
        results.append(("Contract Compliance", "tst_contract_compliance", False, 0))

    # --- Explanation Prompts Tests ---
    try:
        from tests.tst_explanation_prompts import run_tests as run_explanation_prompts_tests
        success, elapsed = run_test_module(
            "tst_explanation_prompts", run_explanation_prompts_tests, args.verbose)
        results.append(("Explanation Prompts", "tst_explanation_prompts", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_explanation_prompts: {e}")
        results.append(("Explanation Prompts", "tst_explanation_prompts", False, 0))

    # --- Gate Evaluator Tests ---
    try:
        from tests.tst_gate_evaluator import run_tests as run_gate_evaluator_tests
        success, elapsed = run_test_module(
            "tst_gate_evaluator", run_gate_evaluator_tests, args.verbose)
        results.append(("Gate Evaluator", "tst_gate_evaluator", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_gate_evaluator: {e}")
        results.append(("Gate Evaluator", "tst_gate_evaluator", False, 0))

    # --- Hypothesis Evaluator Tests ---
    try:
        from tests.tst_hypothesis_evaluator import run_tests as run_hypothesis_evaluator_tests
        success, elapsed = run_test_module(
            "tst_hypothesis_evaluator", run_hypothesis_evaluator_tests, args.verbose)
        results.append(("Hypothesis Evaluator", "tst_hypothesis_evaluator", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_hypothesis_evaluator: {e}")
        results.append(("Hypothesis Evaluator", "tst_hypothesis_evaluator", False, 0))

    # --- MTZ Crosscheck Tests ---
    try:
        from tests.tst_mtz_crosscheck import TESTS as mtz_crosscheck_tests_list
        from tests.tst_utils import run_tests_with_fail_fast
        success, elapsed = run_test_module(
            "tst_mtz_crosscheck",
            lambda: run_tests_with_fail_fast(mtz_crosscheck_tests_list),
            args.verbose)
        results.append(("MTZ Crosscheck", "tst_mtz_crosscheck", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_mtz_crosscheck: {e}")
        results.append(("MTZ Crosscheck", "tst_mtz_crosscheck", False, 0))

    # --- Old Client Compat Tests ---
    try:
        from tests.tst_old_client_compat import run_all_tests as run_old_client_compat_tests
        success, elapsed = run_test_module(
            "tst_old_client_compat", run_old_client_compat_tests, args.verbose)
        results.append(("Old Client Compat", "tst_old_client_compat", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_old_client_compat: {e}")
        results.append(("Old Client Compat", "tst_old_client_compat", False, 0))

    # --- Perceive Stop Checks Tests ---
    try:
        from tests.tst_perceive_stop_checks import run as run_perceive_stop_checks_tests
        success, elapsed = run_test_module(
            "tst_perceive_stop_checks", run_perceive_stop_checks_tests, args.verbose)
        results.append(("Perceive Stop Checks", "tst_perceive_stop_checks", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_perceive_stop_checks: {e}")
        results.append(("Perceive Stop Checks", "tst_perceive_stop_checks", False, 0))

    # --- Phase→Stage Rename Audit ---
    try:
        from tests.audit_phase_rename import run as run_phase_rename_audit
        success, elapsed = run_test_module(
            "audit_phase_rename", run_phase_rename_audit, args.verbose)
        results.append(("Phase Rename Audit", "audit_phase_rename", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import audit_phase_rename: {e}")
        results.append(("Phase Rename Audit", "audit_phase_rename", False, 0))

    # --- Plan Generator Tests ---
    try:
        from tests.tst_plan_generator import run_tests as run_plan_generator_tests
        success, elapsed = run_test_module(
            "tst_plan_generator", run_plan_generator_tests, args.verbose)
        results.append(("Plan Generator", "tst_plan_generator", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_plan_generator: {e}")
        results.append(("Plan Generator", "tst_plan_generator", False, 0))

    # --- Plan Schema Tests ---
    try:
        from tests.tst_plan_schema import run_tests as run_plan_schema_tests
        success, elapsed = run_test_module(
            "tst_plan_schema", run_plan_schema_tests, args.verbose)
        results.append(("Plan Schema", "tst_plan_schema", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_plan_schema: {e}")
        results.append(("Plan Schema", "tst_plan_schema", False, 0))

    # --- Shared Code Imports Tests ---
    try:
        from tests.tst_shared_code_imports import run_all_tests as run_shared_code_imports_tests
        success, elapsed = run_test_module(
            "tst_shared_code_imports", run_shared_code_imports_tests, args.verbose)
        results.append(("Shared Code Imports", "tst_shared_code_imports", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_shared_code_imports: {e}")
        results.append(("Shared Code Imports", "tst_shared_code_imports", False, 0))

    # --- Structure Model Tests ---
    try:
        from tests.tst_structure_model import run_tests as run_structure_model_tests
        success, elapsed = run_test_module(
            "tst_structure_model", run_structure_model_tests, args.verbose)
        results.append(("Structure Model", "tst_structure_model", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_structure_model: {e}")
        results.append(("Structure Model", "tst_structure_model", False, 0))

    # --- Template Tests ---
    try:
        from tests.tst_template import run_all_tests as run_template_tests
        success, elapsed = run_test_module(
            "tst_template", run_template_tests, args.verbose)
        results.append(("Template", "tst_template", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_template: {e}")
        results.append(("Template", "tst_template", False, 0))

    # --- Validation History Tests ---
    try:
        from tests.tst_validation_history import run_tests as run_validation_history_tests
        success, elapsed = run_test_module(
            "tst_validation_history", run_validation_history_tests, args.verbose)
        results.append(("Validation History", "tst_validation_history", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_validation_history: {e}")
        results.append(("Validation History", "tst_validation_history", False, 0))

    # --- Display Data Model Tests (v114) ---
    try:
        from tests.tst_display_data_model import run_tests as run_display_data_model_tests
        success, elapsed = run_test_module(
            "tst_display_data_model", run_display_data_model_tests, args.verbose)
        results.append(("Display Data Model", "tst_display_data_model", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_display_data_model: {e}")
        results.append(("Display Data Model", "tst_display_data_model", False, 0))

    # --- Scenario Tracer Tests (v114) ---
    try:
        from tests.tst_scenario_tracer import run_all as run_scenario_tracer_tests
        success, elapsed = run_test_module(
            "tst_scenario_tracer", run_scenario_tracer_tests, args.verbose)
        results.append(("Scenario Tracer", "tst_scenario_tracer", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_scenario_tracer: {e}")
        results.append(("Scenario Tracer", "tst_scenario_tracer", False, 0))

    # --- P1–P5 Fix Tests (v115) ---
    try:
        from tests.tst_mtz_phases_session_blocking import run_all_tests as run_p1_to_p5_tests
        success, elapsed = run_test_module(
            "tst_mtz_phases_session_blocking", run_p1_to_p5_tests, args.verbose)
        results.append(("MTZ Phases Session Blocking", "tst_mtz_phases_session_blocking", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_mtz_phases_session_blocking: {e}")
        results.append(("MTZ Phases Session Blocking", "tst_mtz_phases_session_blocking", False, 0))

    try:
        from tests.tst_event_formatter_file_categories import run_all_tests as run_n_bugs_tests
        success, elapsed = run_test_module(
            "tst_event_formatter_file_categories", run_n_bugs_tests, args.verbose)
        results.append(("Event Formatter File Categories", "tst_event_formatter_file_categories", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_event_formatter_file_categories: {e}")
        results.append(("Event Formatter File Categories", "tst_event_formatter_file_categories", False, 0))

    # --- v115 Fix 2: Fabricated stop-condition stripping ---
    try:
        from tests.tst_fabricated_stop_condition import run_all_tests as run_stop_condition_tests
        success, elapsed = run_test_module(
            "tst_fabricated_stop_condition", run_stop_condition_tests, args.verbose)
        results.append(("Fabricated Stop Condition", "tst_fabricated_stop_condition", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_fabricated_stop_condition: {e}")
        results.append(("Fabricated Stop Condition", "tst_fabricated_stop_condition", False, 0))


    # --- v115 F1+F3: abort_message for all stops; max-cycles message ---
    try:
        from tests.tst_abort_message_all_stops import run_all_tests as run_f1_f3_tests
        success, elapsed = run_test_module(
            "tst_abort_message_all_stops", run_f1_f3_tests, args.verbose)
        results.append(("Abort Message All Stops", "tst_abort_message_all_stops", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_abort_message_all_stops: {e}")
        results.append(("Abort Message All Stops", "tst_abort_message_all_stops", False, 0))

    # --- v115 F2+F4+F7: API key errors; blocked-program naming; empty command ---
    try:
        from tests.tst_stop_feedback_api_plan_command import run_all_tests as run_f2_f4_f7_tests
        success, elapsed = run_test_module(
            "tst_stop_feedback_api_plan_command", run_f2_f4_f7_tests, args.verbose)
        results.append(("Stop Feedback API Plan Command", "tst_stop_feedback_api_plan_command", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_stop_feedback_api_plan_command: {e}")
        results.append(("Stop Feedback API Plan Command", "tst_stop_feedback_api_plan_command", False, 0))

    # --- v115 F5: workflow_complete stop on rules-only success ---
    try:
        from tests.tst_workflow_complete_stop import run_all_tests as run_f5_tests
        success, elapsed = run_test_module(
            "tst_workflow_complete_stop", run_f5_tests, args.verbose)
        results.append(("Workflow Complete Stop", "tst_workflow_complete_stop", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_workflow_complete_stop: {e}")
        results.append(("Workflow Complete Stop", "tst_workflow_complete_stop", False, 0))

    # --- v115 F6: diagnosed_failure stop for terminal diagnosis ---
    try:
        from tests.tst_diagnosed_failure_stop import run_all_tests as run_f6_tests
        success, elapsed = run_test_module(
            "tst_diagnosed_failure_stop", run_f6_tests, args.verbose)
        results.append(("Diagnosed Failure Stop", "tst_diagnosed_failure_stop", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_diagnosed_failure_stop: {e}")
        results.append(("Diagnosed Failure Stop", "tst_diagnosed_failure_stop", False, 0))

    # --- v115 Follow-up: abort_message for no_workflow_state and after_program_not_available ---
    try:
        from tests.tst_silent_stop_abort_messages import run_all_tests as run_followup_abort_tests
        success, elapsed = run_test_module(
            "tst_silent_stop_abort_messages", run_followup_abort_tests, args.verbose)
        results.append(("Silent Stop Abort Messages", "tst_silent_stop_abort_messages", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_silent_stop_abort_messages: {e}")
        results.append(("Silent Stop Abort Messages", "tst_silent_stop_abort_messages", False, 0))

    # --- P1B: STOP_REASON_CODES, think_stop_override, think_file_overrides ---
    try:
        from tests.tst_stop_reason_codes import run_all_tests as run_p1b_tests
        success, elapsed = run_test_module(
            "tst_stop_reason_codes", run_p1b_tests, args.verbose)
        results.append(("Stop Reason Codes", "tst_stop_reason_codes", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_stop_reason_codes: {e}")
        results.append(("Stop Reason Codes", "tst_stop_reason_codes", False, 0))

    # --- P1B Prompt: stop_reason table injection and guard ---
    try:
        from tests.tst_stop_reason_prompt import run_all_tests as run_p1b_prompt_tests
        success, elapsed = run_test_module(
            "tst_stop_reason_prompt", run_p1b_prompt_tests, args.verbose)
        results.append(("Stop Reason Prompt", "tst_stop_reason_prompt", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_stop_reason_prompt: {e}")
        results.append(("Stop Reason Prompt", "tst_stop_reason_prompt", False, 0))

    # --- P1B Action Item 1: validate() file override existence check ---
    try:
        from tests.tst_validate_file_overrides import run_all_tests as run_action_item_1_tests
        success, elapsed = run_test_module(
            "tst_validate_file_overrides", run_action_item_1_tests, args.verbose)
        results.append(("Validate File Overrides", "tst_validate_file_overrides", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_validate_file_overrides: {e}")
        results.append(("Validate File Overrides", "tst_validate_file_overrides", False, 0))


    # --- Copies Tracking: ASU copies from xtriage and directives ---
    try:
        from tests.tst_copies_tracking import run as run_copies_tracking_tests
        success, elapsed = run_test_module(
            "tst_copies_tracking", run_copies_tracking_tests, args.verbose)
        results.append(("Copies Tracking", "tst_copies_tracking", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_copies_tracking: {e}")
        results.append(("Copies Tracking", "tst_copies_tracking", False, 0))

    # --- Half-map pairs, sigma-A MTZ, file inventory, anomalous gate ---
    try:
        from tests.tst_halfmap_sigmaa_inventory import run as run_halfmap_sigmaa_tests
        success, elapsed = run_test_module(
            "tst_halfmap_sigmaa_inventory", run_halfmap_sigmaa_tests, args.verbose)
        results.append(("Half-map/SigmaA/Inventory", "tst_halfmap_sigmaa_inventory", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_halfmap_sigmaa_inventory: {e}")
        results.append(("Half-map/SigmaA/Inventory", "tst_halfmap_sigmaa_inventory", False, 0))

    try:
        from tests.tst_fix_verification import run as run_fix_verification_tests
        success, elapsed = run_test_module(
            "tst_fix_verification", run_fix_verification_tests, args.verbose)
        results.append(("Fix Verification", "tst_fix_verification", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_fix_verification: {e}")
        results.append(("Fix Verification", "tst_fix_verification", False, 0))

    # --- Error Classifier Tests (v115) ---
    try:
        from tests.tst_error_classifier import run_all_tests as run_error_classifier_tests
        success, elapsed = run_test_module(
            "tst_error_classifier", run_error_classifier_tests, args.verbose)
        results.append(("Error Classifier", "tst_error_classifier", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_error_classifier: {e}")
        results.append(("Error Classifier", "tst_error_classifier", False, 0))

    # --- PHIL Validation Tests (v115) ---
    try:
        from tests.tst_phil_validation import run_all_tests as run_phil_validation_tests
        success, elapsed = run_test_module(
            "tst_phil_validation", run_phil_validation_tests, args.verbose)
        results.append(("PHIL Validation", "tst_phil_validation", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phil_validation: {e}")
        results.append(("PHIL Validation", "tst_phil_validation", False, 0))

    # --- Phase 3 + Bug 5 Tests (v115.07) ---
    try:
        from tests.tst_phase3_bug5 import run_tests as run_phase3_bug5_tests
        success, elapsed = run_test_module(
            "tst_phase3_bug5", run_phase3_bug5_tests, args.verbose)
        results.append(("Phase 3 + Bug 5", "tst_phase3_bug5", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase3_bug5: {e}")
        results.append(("Phase 3 + Bug 5", "tst_phase3_bug5", False, 0))

    # --- Stop Condition Fix Tests (v115) ---
    try:
        from tests.tst_stop_condition_fix import run_all_tests as run_stop_condition_fix_tests
        success, elapsed = run_test_module(
            "tst_stop_condition_fix", run_stop_condition_fix_tests, args.verbose)
        results.append(("Stop Condition Fix", "tst_stop_condition_fix", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_stop_condition_fix: {e}")
        results.append(("Stop Condition Fix", "tst_stop_condition_fix", False, 0))

    # --- Intent Classifier Tests (v115.01) ---
    try:
        from tests.tst_intent_classifier import run_all_tests as run_intent_classifier_tests
        success, elapsed = run_test_module(
            "tst_intent_classifier", run_intent_classifier_tests, args.verbose)
        results.append(("Intent Classifier", "tst_intent_classifier", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_intent_classifier: {e}")
        results.append(("Intent Classifier", "tst_intent_classifier", False, 0))

    # --- N-Bugs Tests (v115) ---
    try:
        from tests.tst_n_bugs import run_all_tests as run_n_bugs_extra_tests
        success, elapsed = run_test_module(
            "tst_n_bugs", run_n_bugs_extra_tests, args.verbose)
        results.append(("N-Bugs", "tst_n_bugs", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_n_bugs: {e}")
        results.append(("N-Bugs", "tst_n_bugs", False, 0))

    # --- Action Item 1 Tests ---
    try:
        from tests.tst_action_item_1 import run_all_tests as run_action_item_1_tests
        success, elapsed = run_test_module(
            "tst_action_item_1", run_action_item_1_tests, args.verbose)
        results.append(("Action Item 1", "tst_action_item_1", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_action_item_1: {e}")
        results.append(("Action Item 1", "tst_action_item_1", False, 0))

    # --- Backward Compatibility Tests (v115.05/v115.06) ---
    try:
        from tests.tst_backward_compat import run_all_tests as run_backward_compat_tests
        success, elapsed = run_test_module(
            "tst_backward_compat", run_backward_compat_tests, args.verbose)
        results.append(("Backward Compat", "tst_backward_compat", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_backward_compat: {e}")
        results.append(("Backward Compat", "tst_backward_compat", False, 0))

    # ===================================================================
    # SYSTEMATIC TESTING FRAMEWORK (v115.08)
    # ===================================================================
    # Bottom-up testing of system boundaries where bugs hide.
    # These tests exercise real production code paths end-to-end
    # and were designed after 4 code reviews found 3 critical bugs
    # that 2,186 existing unit tests missed.
    # ===================================================================

    # --- Phase 0: Static Audit (automated gate) ---
    try:
        from tests.tst_phase0_static_audit import run_all_tests as run_phase0_tests
        success, elapsed = run_test_module(
            "tst_phase0_static_audit", run_phase0_tests, args.verbose)
        results.append(("Phase 0: Static Audit", "tst_phase0_static_audit", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase0_static_audit: {e}")
        results.append(("Phase 0: Static Audit", "tst_phase0_static_audit", False, 0))

    # --- Phase 1: Contract Gaps (coverage map) ---
    try:
        from tests.tst_phase1_contract_gaps import run_all_tests as run_phase1_tests
        success, elapsed = run_test_module(
            "tst_phase1_contract_gaps", run_phase1_tests, args.verbose)
        results.append(("Phase 1: Contract Gaps", "tst_phase1_contract_gaps", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase1_contract_gaps: {e}")
        results.append(("Phase 1: Contract Gaps", "tst_phase1_contract_gaps", False, 0))

    # --- Phase 2: Path Consistency (YAML vs hardcoded) ---
    try:
        from tests.tst_phase2_path_consistency import run_all_tests as run_phase2_tests
        success, elapsed = run_test_module(
            "tst_phase2_path_consistency", run_phase2_tests, args.verbose)
        results.append(("Phase 2: Path Consistency", "tst_phase2_path_consistency", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase2_path_consistency: {e}")
        results.append(("Phase 2: Path Consistency", "tst_phase2_path_consistency", False, 0))

    # --- Phase 3: Session Round-Trip (symmetry + invariants) ---
    try:
        from tests.tst_phase3_serialization_symmetry import run_all_tests as run_phase3_tests
        success, elapsed = run_test_module(
            "tst_phase3_serialization_symmetry", run_phase3_tests, args.verbose)
        results.append(("Phase 3: Session Round-Trip", "tst_phase3_serialization_symmetry", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase3_serialization_symmetry: {e}")
        results.append(("Phase 3: Session Round-Trip", "tst_phase3_serialization_symmetry", False, 0))

    # --- Phase 4: History Flag Consistency ---
    try:
        from tests.tst_phase4_history_flags import run_all_tests as run_phase4_tests
        success, elapsed = run_test_module(
            "tst_phase4_history_flags", run_phase4_tests, args.verbose)
        results.append(("Phase 4: History Flags", "tst_phase4_history_flags", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase4_history_flags: {e}")
        results.append(("Phase 4: History Flags", "tst_phase4_history_flags", False, 0))

    # --- Phase 5: Error Classification Consistency ---
    try:
        from tests.tst_phase5_error_classification import run_all_tests as run_phase5_tests
        success, elapsed = run_test_module(
            "tst_phase5_error_classification", run_phase5_tests, args.verbose)
        results.append(("Phase 5: Error Classification", "tst_phase5_error_classification", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase5_error_classification: {e}")
        results.append(("Phase 5: Error Classification", "tst_phase5_error_classification", False, 0))

    # --- Phase 6: Category-Consumer Alignment ---
    try:
        from tests.tst_phase6_category_consumer import run_all_tests as run_phase6_tests
        success, elapsed = run_test_module(
            "tst_phase6_category_consumer", run_phase6_tests, args.verbose)
        results.append(("Phase 6: Category-Consumer", "tst_phase6_category_consumer", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase6_category_consumer: {e}")
        results.append(("Phase 6: Category-Consumer", "tst_phase6_category_consumer", False, 0))

    # --- Phase 7: Routing Simulation (32 tutorials x 3 cycles) ---
    if not args.quick:
        try:
            from tests.tst_phase7_routing_simulation import run_all_tests as run_phase7_tests
            success, elapsed = run_test_module(
                "tst_phase7_routing_simulation", run_phase7_tests, args.verbose)
            results.append(("Phase 7: Routing Simulation", "tst_phase7_routing_simulation", success, elapsed))
        except ImportError as e:
            print(f"\u26a0\ufe0f  Could not import tst_phase7_routing_simulation: {e}")
            results.append(("Phase 7: Routing Simulation", "tst_phase7_routing_simulation", False, 0))
    else:
        print("\n  Skipping Phase 7 (routing simulation) in --quick mode")

    # --- Phase 8: Command Building Simulation ---
    if not args.quick:
        try:
            from tests.tst_phase8_command_building import run_all_tests as run_phase8_tests
            success, elapsed = run_test_module(
                "tst_phase8_command_building", run_phase8_tests, args.verbose)
            results.append(("Phase 8: Command Building", "tst_phase8_command_building", success, elapsed))
        except ImportError as e:
            print(f"\u26a0\ufe0f  Could not import tst_phase8_command_building: {e}")
            results.append(("Phase 8: Command Building", "tst_phase8_command_building", False, 0))
    else:
        print("  Skipping Phase 8 (command building) in --quick mode")

    # --- Phase 9: LLM Perturbation Tests ---
    if not args.quick:
        try:
            from tests.tst_phase9_llm_perturbation import run_all_tests as run_phase9_tests
            success, elapsed = run_test_module(
                "tst_phase9_llm_perturbation", run_phase9_tests, args.verbose)
            results.append(("Phase 9: LLM Perturbation", "tst_phase9_llm_perturbation", success, elapsed))
        except ImportError as e:
            print(f"\u26a0\ufe0f  Could not import tst_phase9_llm_perturbation: {e}")
            results.append(("Phase 9: LLM Perturbation", "tst_phase9_llm_perturbation", False, 0))
    else:
        print("  Skipping Phase 9 (LLM perturbation) in --quick mode")

    # --- Summary ---
    total_elapsed = time.time() - total_start

    print("\n")
    print("="*60)
    print("TEST SUMMARY")
    print("="*60)

    passed = 0
    failed = 0
    failed_modules = []

    for name, module, success, elapsed in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"  {name:30s} {status} ({elapsed:.2f}s)")
        if success:
            passed += 1
        else:
            failed += 1
            failed_modules.append((name, module))

    print("-"*60)
    print(f"  Total: {passed} passed, {failed} failed ({total_elapsed:.2f}s)")
    print("="*60)

    if failed > 0:
        print("\n❌ SOME TESTS FAILED")
        print("\nTO RERUN FAILED TESTS:")
        for name, module in failed_modules:
            print(f"  # {name}")
            print(f"  phenix.python $PHENIX/modules/cctbx_project/libtbx/langchain/tests/{module}.py")
        sys.exit(1)
    else:
        print("\n✅ ALL TESTS PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
