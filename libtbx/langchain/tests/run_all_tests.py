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
  S10. General Resolver - ACTION_TABLE, _detect_actions, _resolve_after_program (20 tests, 91 assertions)
  S11. Skip to Program - StructurePlan.skip_to_program() + format_plan_header display (30 tests, v116.10)
  S12. Directive Validator Typo Hints - difflib-based suggestion algorithm and
       'Did you mean?' hints in after_program issues (13 tests, v116.10 Phase 4a)
  S13. User Advice Filter - 'X and stop' / 'X then stop' / 'X, stop' phrasings no longer
       truncate valid_programs to [STOP] (16 tests, v116.10 Phase 1)
  S14. After-Program Prompt Alignment - target-not-now framing replaces
       'CRITICAL: You MUST run X / Do NOT...run X instead!' (12 tests, v116.10 Phase 6a)
  S15. Data-Input Filter - xtriage removed when no .mtz, mtriage removed
       when no map (16 tests, v116.10 Phase 4b)
  S16. Protocol Version Invariants - CURRENT_PROTOCOL_VERSION drift detection
       and validate_contract() function (15 tests, v116.10 Phase 2)
  S17. Sequence-Only Routing - _detect_xray_step routes to obtain_model when
       no .mtz data is present and has_sequence=True (10 tests, v116.10 Phase 6b)
  S18. Standalone Programs Consistency - _STANDALONE_PROGRAMS / _NEEDS_PLAN_PROGRAMS
       de-duplication + drift detection against _ACTION_TABLE
       (8 tests, v116.10 Phase 3a/3b)
  S19. Dock-and-Stop Behavior Change - phenix.map_symmetry and phenix.dock_in_map
       reclassified as standalone (5 tests, v116.10 Phase 3d)
  S20. CC Key Extraction - _generate_structure_report success-detection logic
       recognizes the canonical "model_map_cc" key, fixing misleading
       "SESSION STOPPED - INCOMPLETE" banner on successful cryo-EM runs
       (11 tests, v116.10)
  S21. File Encoding - text-mode open() calls specify encoding='utf-8',
       preventing UnicodeDecodeError crashes on non-UTF-8 system locales
       (Chinese/Japanese Windows in particular) (5 tests, v116.10)
  S22. Ligand Workflow Restart - _detect_xray_step advances past "analyze"
       when downstream programs have run; best_files_tracker treats
       ligand_fit_output as new best by inheriting metrics from prior
       refined model. Fixes nsf-d2-ligand tutorial restarting at cycle 3
       (14 tests, v116.10 Phase 6c)
  S23. Initialize Plan Smoke - behavioral tests for _initialize_plan_inner
       decision tree. Verifies side effects (intent rewrite, stop_after
       clearing, plan-generation gating) for each branch of the function.
       Complements tst_dock_and_stop.py's decision-tree traces with actual
       behavioral assertions (9 tests, v116.10 Tier 1 follow-up)
  S24. Phase 3d Motivating Tutorial - integration-level verification of the
       Phase 3d dock-and-stop fix. Skip-aware against an optional session
       fixture; runs decision-tree + YAML format checks unconditionally.
       Closes the highest-priority verification gap from the v116.10 review
       (3 tests, v116.10 Tier 1 follow-up)
  S25. Program Requirements - declarative `requirements:` schema for
       programs.yaml. Verifies the _check_requirements parser, the filter
       integration (autobuild scenarios), and backward-compat (programs
       without the block are unaffected). Plugs the autobuild gap where
       the program could be picked and crash at runtime
       (40 tests, v116.10 Tier 2.1 + 2.1.1)

  S26. Stop Condition False Positive - v116.11 fix for AF_7mjs regression.
       The preprocessor-inserted "**Stop Condition**: None" header was
       interpreted as a real stop signal by _resolve_after_program,
       causing the planner to skip prerequisite stages.  Tests the
       strengthened _strip_preprocessor_stop_condition regex, the
       defense-in-depth strip in the resolver, and the start_with_program
       suppression for preprocessed advice
       (15 tests, v116.11)

  S27. Cryo-EM Stop Validation - v116.12 Fix #1 for AF_7mjs Stage 5 skip.
       The cryo-EM auto-stop condition (CC > 0.70) lacked the
       validation_done check that the X-ray path has.  This caused
       the agent to stop after refinement without running validation
       (phenix.molprobity).  Tests verify the validation_done check
       in _analyze_cryoem_trend, accept both phenix.molprobity and
       phenix.validation_cryoem as cryo-EM validation, and confirm
       the X-ray path is unaffected
       (16 tests, v116.12 Fix #1)

  S28. PLAN AUTO-STOP Suppression - v116.12 Fix #2 defense-in-depth.
       New elif suppresses AUTO-STOP when workflow_engine reports
       step="validate" with validation_done=False.  Diagnostic context
       dump in the AUTO-STOP path records plan_has_pending_stages,
       step, validation_done, after_program, experiment_type to aid
       future debugging.  Defense-in-depth for cases where
       plan_has_pending_stages doesn't fire
       (11 tests, v116.12 Fix #2)

  S29. Cryo-EM Metric Evaluator Validation - v116.13 actual fix.
       v116.12 patched metrics_analyzer.py but graph_nodes.py sets
       USE_YAML_METRICS=True (default), which delegates trend
       analysis to MetricEvaluator.analyze_trend() in
       metric_evaluator.py — the actual production path.  v116.12's
       fix was in dead code.  v116.13 adds the validation_done check
       to MetricEvaluator's cryo-EM SUCCESS branch, mirroring the
       X-ray pattern (which already had it).  Tests verify the fix
       at the class level, X-ray path unchanged, and validation
       program recognition (phenix.molprobity, phenix.validation_cryoem)
       (15 tests, v116.13)
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

    # --- FORCE_NO_AI_SERVER Tests (v120 Task 1) ---
    try:
        from tests.tst_force_no_ai_server import run_all_tests as run_force_no_ai_server_tests
        success, elapsed = run_test_module(
            "tst_force_no_ai_server", run_force_no_ai_server_tests, args.verbose)
        results.append(("FORCE_NO_AI_SERVER", "tst_force_no_ai_server", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_force_no_ai_server: {e}")
        results.append(("FORCE_NO_AI_SERVER", "tst_force_no_ai_server", False, 0))

    # --- SUPPORTED_PROVIDERS Single Source Tests (v120 Phase 2) ---
    try:
        from tests.tst_supported_providers_single_source import run_all_tests as run_supported_providers_tests
        success, elapsed = run_test_module(
            "tst_supported_providers_single_source", run_supported_providers_tests, args.verbose)
        results.append(("SUPPORTED_PROVIDERS Single Source", "tst_supported_providers_single_source", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_supported_providers_single_source: {e}")
        results.append(("SUPPORTED_PROVIDERS Single Source", "tst_supported_providers_single_source", False, 0))

    # --- Portkey + Anthropic Provider Wiring Tests (v120 Phase 3) ---
    try:
        from tests.tst_v120_providers import run_all_tests as run_v120_providers_tests
        success, elapsed = run_test_module(
            "tst_v120_providers", run_v120_providers_tests, args.verbose)
        results.append(("v120 Providers (Portkey/Anthropic)", "tst_v120_providers", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_v120_providers: {e}")
        results.append(("v120 Providers (Portkey/Anthropic)", "tst_v120_providers", False, 0))

    try:
        from tests.tst_llm_unavailable_notice import run_all_tests as run_llm_unavailable_tests
        success, elapsed = run_test_module(
            "tst_llm_unavailable_notice", run_llm_unavailable_tests, args.verbose)
        results.append(("LLM-unavailable notice", "tst_llm_unavailable_notice", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_llm_unavailable_notice: {e}")
        results.append(("LLM-unavailable notice", "tst_llm_unavailable_notice", False, 0))

    # --- resolve_cryo_em ncs_file never-inject ---
    try:
        from tests.tst_resolve_cryo_em_ncs_inject import run_all_tests as run_resolve_ncs_tests
        success, elapsed = run_test_module(
            "tst_resolve_cryo_em_ncs_inject", run_resolve_ncs_tests, args.verbose)
        results.append(("resolve_cryo_em ncs_file", "tst_resolve_cryo_em_ncs_inject", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_resolve_cryo_em_ncs_inject: {e}")
        results.append(("resolve_cryo_em ncs_file", "tst_resolve_cryo_em_ncs_inject", False, 0))

    # --- portkey embedding model name + silent-fallback guard ---
    try:
        from tests.tst_portkey_embedding_model import run_all_tests as run_portkey_embed_tests
        success, elapsed = run_test_module(
            "tst_portkey_embedding_model", run_portkey_embed_tests, args.verbose)
        results.append(("portkey embedding model", "tst_portkey_embedding_model", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_portkey_embedding_model: {e}")
        results.append(("portkey embedding model", "tst_portkey_embedding_model", False, 0))

    # --- provider lock (provider-specific install) ---
    try:
        from tests.tst_provider_lock import run_all_tests as run_provider_lock_tests
        success, elapsed = run_test_module(
            "tst_provider_lock", run_provider_lock_tests, args.verbose)
        results.append(("provider lock", "tst_provider_lock", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_provider_lock: {e}")
        results.append(("provider lock", "tst_provider_lock", False, 0))

    # --- _safe_float consolidation + sanity_checker metric coercion ---
    try:
        from tests.tst_safe_float_consolidation import run_all_tests as run_safe_float_tests
        success, elapsed = run_test_module(
            "tst_safe_float_consolidation", run_safe_float_tests, args.verbose)
        results.append(("safe_float consolidation", "tst_safe_float_consolidation", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_safe_float_consolidation: {e}")
        results.append(("safe_float consolidation", "tst_safe_float_consolidation", False, 0))

    # --- dropped-resolution premature-stop (rfree target) ---
    try:
        from tests.tst_rfree_resolution_stop import run_all_tests as run_rfree_res_tests
        success, elapsed = run_test_module(
            "tst_rfree_resolution_stop", run_rfree_res_tests, args.verbose)
        results.append(("rfree resolution stop", "tst_rfree_resolution_stop", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_rfree_resolution_stop: {e}")
        results.append(("rfree resolution stop", "tst_rfree_resolution_stop", False, 0))

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

    # --- Raw-Advice Authoritative Extraction Tests (Step 1) ---
    try:
        from tests.tst_extract_raw_advice import run_all_tests as run_extract_raw_advice_tests
        success, elapsed = run_test_module(
            "tst_extract_raw_advice", run_extract_raw_advice_tests, args.verbose)
        results.append(("Raw-Advice Authoritative Extraction", "tst_extract_raw_advice", success, elapsed))
    except ImportError as e:
        print(f"⚠️  Could not import tst_extract_raw_advice: {e}")
        results.append(("Raw-Advice Authoritative Extraction", "tst_extract_raw_advice", False, 0))

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

    # --- General Resolver Tests (v115.10) ---
    try:
        from tests.tst_general_resolver import run_all_tests as run_general_resolver_tests
        success, elapsed = run_test_module(
            "tst_general_resolver", run_general_resolver_tests, args.verbose)
        results.append(("General Resolver", "tst_general_resolver", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_general_resolver: {e}")
        results.append(("General Resolver", "tst_general_resolver", False, 0))

    # --- Skip-to-Program Tests (v116.10) ---
    try:
        from tests.tst_skip_to_program import run_all_tests as run_skip_to_program_tests
        success, elapsed = run_test_module(
            "tst_skip_to_program", run_skip_to_program_tests, args.verbose)
        results.append(("Skip to Program", "tst_skip_to_program", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_skip_to_program: {e}")
        results.append(("Skip to Program", "tst_skip_to_program", False, 0))

    # --- Directive Validator Typo Hints (v116.10 Phase 4a) ---
    try:
        from tests.tst_directive_validator_typo_hints import (
            run_all_tests as run_directive_validator_typo_hints_tests)
        success, elapsed = run_test_module(
            "tst_directive_validator_typo_hints",
            run_directive_validator_typo_hints_tests, args.verbose)
        results.append(("Directive Validator Typo Hints",
                       "tst_directive_validator_typo_hints",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_directive_validator_typo_hints: {e}")
        results.append(("Directive Validator Typo Hints",
                       "tst_directive_validator_typo_hints", False, 0))

    # --- User Advice Filter (v116.10 Phase 1) ---
    try:
        from tests.tst_user_advice_filter import (
            run_all_tests as run_user_advice_filter_tests)
        success, elapsed = run_test_module(
            "tst_user_advice_filter",
            run_user_advice_filter_tests, args.verbose)
        results.append(("User Advice Filter",
                       "tst_user_advice_filter",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_user_advice_filter: {e}")
        results.append(("User Advice Filter",
                       "tst_user_advice_filter", False, 0))

    # --- After-Program Prompt Alignment (v116.10 Phase 6a) ---
    try:
        from tests.tst_after_program_prompt import (
            run_all_tests as run_after_program_prompt_tests)
        success, elapsed = run_test_module(
            "tst_after_program_prompt",
            run_after_program_prompt_tests, args.verbose)
        results.append(("After-Program Prompt Alignment",
                       "tst_after_program_prompt",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_after_program_prompt: {e}")
        results.append(("After-Program Prompt Alignment",
                       "tst_after_program_prompt", False, 0))

    # --- Data-Input Filter (v116.10 Phase 4b) ---
    try:
        from tests.tst_data_input_filter import (
            run_all_tests as run_data_input_filter_tests)
        success, elapsed = run_test_module(
            "tst_data_input_filter",
            run_data_input_filter_tests, args.verbose)
        results.append(("Data-Input Filter",
                       "tst_data_input_filter",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_data_input_filter: {e}")
        results.append(("Data-Input Filter",
                       "tst_data_input_filter", False, 0))

    # --- Protocol Version Invariants (v116.10 Phase 2) ---
    try:
        from tests.tst_protocol_version import (
            run_all_tests as run_protocol_version_tests)
        success, elapsed = run_test_module(
            "tst_protocol_version",
            run_protocol_version_tests, args.verbose)
        results.append(("Protocol Version Invariants",
                       "tst_protocol_version",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_protocol_version: {e}")
        results.append(("Protocol Version Invariants",
                       "tst_protocol_version", False, 0))

    # --- Sequence-Only Routing (v116.10 Phase 6b) ---
    try:
        from tests.tst_sequence_only_routing import (
            run_all_tests as run_sequence_only_routing_tests)
        success, elapsed = run_test_module(
            "tst_sequence_only_routing",
            run_sequence_only_routing_tests, args.verbose)
        results.append(("Sequence-Only Routing",
                       "tst_sequence_only_routing",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_sequence_only_routing: {e}")
        results.append(("Sequence-Only Routing",
                       "tst_sequence_only_routing", False, 0))

    # --- Standalone Programs Consistency (v116.10 Phase 3a/3b) ---
    try:
        from tests.tst_standalone_consistency import (
            run_all_tests as run_standalone_consistency_tests)
        success, elapsed = run_test_module(
            "tst_standalone_consistency",
            run_standalone_consistency_tests, args.verbose)
        results.append(("Standalone Programs Consistency",
                       "tst_standalone_consistency",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_standalone_consistency: {e}")
        results.append(("Standalone Programs Consistency",
                       "tst_standalone_consistency", False, 0))

    # --- Dock-and-Stop Behavior Change (v116.10 Phase 3d) ---
    try:
        from tests.tst_dock_and_stop import (
            run_all_tests as run_dock_and_stop_tests)
        success, elapsed = run_test_module(
            "tst_dock_and_stop",
            run_dock_and_stop_tests, args.verbose)
        results.append(("Dock-and-Stop Behavior Change",
                       "tst_dock_and_stop",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_dock_and_stop: {e}")
        results.append(("Dock-and-Stop Behavior Change",
                       "tst_dock_and_stop", False, 0))

    # --- CC Key Extraction (v116.10) ---
    try:
        from tests.tst_cc_key_extraction import (
            run_all_tests as run_cc_key_extraction_tests)
        success, elapsed = run_test_module(
            "tst_cc_key_extraction",
            run_cc_key_extraction_tests, args.verbose)
        results.append(("CC Key Extraction",
                       "tst_cc_key_extraction",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_cc_key_extraction: {e}")
        results.append(("CC Key Extraction",
                       "tst_cc_key_extraction", False, 0))

    # --- File Encoding (v116.10) ---
    try:
        from tests.tst_file_encoding import (
            run_all_tests as run_file_encoding_tests)
        success, elapsed = run_test_module(
            "tst_file_encoding",
            run_file_encoding_tests, args.verbose)
        results.append(("File Encoding",
                       "tst_file_encoding",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_file_encoding: {e}")
        results.append(("File Encoding",
                       "tst_file_encoding", False, 0))

    # --- Ligand Workflow Restart (v116.10 Phase 6c) ---
    try:
        from tests.tst_ligand_workflow_restart import (
            run_all_tests as run_ligand_workflow_restart_tests)
        success, elapsed = run_test_module(
            "tst_ligand_workflow_restart",
            run_ligand_workflow_restart_tests, args.verbose)
        results.append(("Ligand Workflow Restart",
                       "tst_ligand_workflow_restart",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_ligand_workflow_restart: {e}")
        results.append(("Ligand Workflow Restart",
                       "tst_ligand_workflow_restart", False, 0))

    # --- Initialize Plan Smoke (v116.10 Tier 1 follow-up) ---
    try:
        from tests.tst_initialize_plan_smoke import (
            run_all_tests as run_initialize_plan_smoke_tests)
        success, elapsed = run_test_module(
            "tst_initialize_plan_smoke",
            run_initialize_plan_smoke_tests, args.verbose)
        results.append(("Initialize Plan Smoke",
                       "tst_initialize_plan_smoke",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_initialize_plan_smoke: {e}")
        results.append(("Initialize Plan Smoke",
                       "tst_initialize_plan_smoke", False, 0))

    # --- Phase 3d Motivating Tutorial (v116.10 Tier 1 follow-up) ---
    try:
        from tests.tst_phase3d_motivating_tutorial import (
            run_all_tests as run_phase3d_motivating_tutorial_tests)
        success, elapsed = run_test_module(
            "tst_phase3d_motivating_tutorial",
            run_phase3d_motivating_tutorial_tests, args.verbose)
        results.append(("Phase 3d Motivating Tutorial",
                       "tst_phase3d_motivating_tutorial",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_phase3d_motivating_tutorial: {e}")
        results.append(("Phase 3d Motivating Tutorial",
                       "tst_phase3d_motivating_tutorial", False, 0))

    # --- Program Requirements (v116.10 Tier 2.1 — declarative schema) ---
    try:
        from tests.tst_program_requirements import (
            run_all_tests as run_program_requirements_tests)
        success, elapsed = run_test_module(
            "tst_program_requirements",
            run_program_requirements_tests, args.verbose)
        results.append(("Program Requirements",
                       "tst_program_requirements",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_program_requirements: {e}")
        results.append(("Program Requirements",
                       "tst_program_requirements", False, 0))

    # --- Stop Condition False Positive (v116.11) ---
    try:
        from tests.tst_stop_condition_false_positive import (
            run_all_tests as run_stop_condition_tests)
        success, elapsed = run_test_module(
            "tst_stop_condition_false_positive",
            run_stop_condition_tests, args.verbose)
        results.append(("Stop Condition False Positive",
                       "tst_stop_condition_false_positive",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_stop_condition_false_positive: {e}")
        results.append(("Stop Condition False Positive",
                       "tst_stop_condition_false_positive", False, 0))

    # --- Cryo-EM Stop Validation (v116.12 Fix #1) ---
    try:
        from tests.tst_cryoem_stop_validation import (
            run_all_tests as run_cryoem_stop_tests)
        success, elapsed = run_test_module(
            "tst_cryoem_stop_validation",
            run_cryoem_stop_tests, args.verbose)
        results.append(("Cryo-EM Stop Validation",
                       "tst_cryoem_stop_validation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_cryoem_stop_validation: {e}")
        results.append(("Cryo-EM Stop Validation",
                       "tst_cryoem_stop_validation", False, 0))

    # --- PLAN AUTO-STOP Suppression (v116.12 Fix #2) ---
    try:
        from tests.tst_plan_autostop_validation_suppression import (
            run_all_tests as run_plan_autostop_tests)
        success, elapsed = run_test_module(
            "tst_plan_autostop_validation_suppression",
            run_plan_autostop_tests, args.verbose)
        results.append(("PLAN AUTO-STOP Suppression",
                       "tst_plan_autostop_validation_suppression",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_plan_autostop_validation_suppression: {e}")
        results.append(("PLAN AUTO-STOP Suppression",
                       "tst_plan_autostop_validation_suppression", False, 0))

    # --- Cryo-EM Metric Evaluator Validation (v116.13) ---
    try:
        from tests.tst_cryoem_metric_evaluator_validation import (
            run_all_tests as run_cryoem_metric_eval_tests)
        success, elapsed = run_test_module(
            "tst_cryoem_metric_evaluator_validation",
            run_cryoem_metric_eval_tests, args.verbose)
        results.append(("Cryo-EM Metric Evaluator Validation",
                       "tst_cryoem_metric_evaluator_validation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_cryoem_metric_evaluator_validation: {e}")
        results.append(("Cryo-EM Metric Evaluator Validation",
                       "tst_cryoem_metric_evaluator_validation", False, 0))

    # --- Directive Extractor Sentinel Space Group (v116.15) ---
    try:
        from tests.tst_directive_extractor_sentinel_space_group import (
            run_all_tests as run_sentinel_sg_tests)
        success, elapsed = run_test_module(
            "tst_directive_extractor_sentinel_space_group",
            run_sentinel_sg_tests, args.verbose)
        results.append(("Directive Extractor Sentinel Space Group",
                       "tst_directive_extractor_sentinel_space_group",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_directive_extractor_sentinel_space_group: {e}")
        results.append(("Directive Extractor Sentinel Space Group",
                       "tst_directive_extractor_sentinel_space_group",
                       False, 0))

    # --- Directive Extractor Behavior Regression (v116.15) ---
    try:
        from tests.tst_directive_extractor_behavior_regression import (
            run_all_tests as run_dx_regression_tests)
        success, elapsed = run_test_module(
            "tst_directive_extractor_behavior_regression",
            run_dx_regression_tests, args.verbose)
        results.append(("Directive Extractor Behavior Regression",
                       "tst_directive_extractor_behavior_regression",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_directive_extractor_behavior_regression: {e}")
        results.append(("Directive Extractor Behavior Regression",
                       "tst_directive_extractor_behavior_regression",
                       False, 0))

    # --- Validate-Step After-Program Guard (v116.17) ---
    try:
        from tests.tst_validate_step_after_program_guard import (
            run_all_tests as run_validate_guard_tests)
        success, elapsed = run_test_module(
            "tst_validate_step_after_program_guard",
            run_validate_guard_tests, args.verbose)
        results.append(("Validate-Step After-Program Guard",
                       "tst_validate_step_after_program_guard",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_validate_step_after_program_guard: {e}")
        results.append(("Validate-Step After-Program Guard",
                       "tst_validate_step_after_program_guard",
                       False, 0))

    # --- Directive Extractor Grounding Guardrail (v116.19) ---
    try:
        from tests.tst_directive_extractor_grounding import (
            run_all_tests as run_grounding_tests)
        success, elapsed = run_test_module(
            "tst_directive_extractor_grounding",
            run_grounding_tests, args.verbose)
        results.append(("Directive Extractor Grounding Guardrail",
                       "tst_directive_extractor_grounding",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_directive_extractor_grounding: {e}")
        results.append(("Directive Extractor Grounding Guardrail",
                       "tst_directive_extractor_grounding",
                       False, 0))

    # --- Grounding + stop_after_requested Interaction (v117.1) ---
    try:
        from tests.tst_grounding_stop_after_requested import (
            run_all_tests as run_grounding_flag_tests)
        success, elapsed = run_test_module(
            "tst_grounding_stop_after_requested",
            run_grounding_flag_tests, args.verbose)
        results.append(("Grounding + stop_after_requested Interaction",
                       "tst_grounding_stop_after_requested",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_grounding_stop_after_requested: {e}")
        results.append(("Grounding + stop_after_requested Interaction",
                       "tst_grounding_stop_after_requested",
                       False, 0))

    # --- After-Program Fill-In from Raw Advice (v117.2) ---
    try:
        from tests.tst_after_program_fill_from_raw import (
            run_all_tests as run_after_program_fill_tests)
        success, elapsed = run_test_module(
            "tst_after_program_fill_from_raw",
            run_after_program_fill_tests, args.verbose)
        results.append(("After-Program Fill-In from Raw Advice",
                       "tst_after_program_fill_from_raw",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_after_program_fill_from_raw: {e}")
        results.append(("After-Program Fill-In from Raw Advice",
                       "tst_after_program_fill_from_raw",
                       False, 0))

    # --- Extended Stop-Intent Phrasings (v117.3) ---
    try:
        from tests.tst_extended_stop_phrasings import (
            run_all_tests as run_extended_stop_tests)
        success, elapsed = run_test_module(
            "tst_extended_stop_phrasings",
            run_extended_stop_tests, args.verbose)
        results.append(("Extended Stop-Intent Phrasings",
                       "tst_extended_stop_phrasings",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_extended_stop_phrasings: {e}")
        results.append(("Extended Stop-Intent Phrasings",
                       "tst_extended_stop_phrasings",
                       False, 0))

    # --- Preprocessor File Override (v118.A2) ---
    try:
        from tests.tst_preprocessor_file_override import (
            run_all_tests as run_preprocessor_file_override_tests)
        success, elapsed = run_test_module(
            "tst_preprocessor_file_override",
            run_preprocessor_file_override_tests, args.verbose)
        results.append(("Preprocessor File Override",
                       "tst_preprocessor_file_override",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_preprocessor_file_override: {e}")
        results.append(("Preprocessor File Override",
                       "tst_preprocessor_file_override",
                       False, 0))

    # --- Directive Layer Diagnostics (v118.C-prime) ---
    try:
        from tests.tst_directive_layer_diagnostics import (
            run_all_tests as run_directive_layer_diagnostics_tests)
        success, elapsed = run_test_module(
            "tst_directive_layer_diagnostics",
            run_directive_layer_diagnostics_tests, args.verbose)
        results.append(("Directive Layer Diagnostics",
                       "tst_directive_layer_diagnostics",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_directive_layer_diagnostics: {e}")
        results.append(("Directive Layer Diagnostics",
                       "tst_directive_layer_diagnostics",
                       False, 0))

    # --- PHIL Namespace Cleaner (v118.B) ---
    try:
        from tests.tst_phil_namespace_cleaner import (
            run_all_tests as run_phil_namespace_cleaner_tests)
        success, elapsed = run_test_module(
            "tst_phil_namespace_cleaner",
            run_phil_namespace_cleaner_tests, args.verbose)
        results.append(("PHIL Namespace Cleaner",
                       "tst_phil_namespace_cleaner",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_phil_namespace_cleaner: {e}")
        results.append(("PHIL Namespace Cleaner",
                       "tst_phil_namespace_cleaner",
                       False, 0))

    # --- Extraction-Failure Visibility (v118.E) ---
    try:
        from tests.tst_extraction_failure_visibility import (
            run_all_tests as run_extraction_failure_visibility_tests)
        success, elapsed = run_test_module(
            "tst_extraction_failure_visibility",
            run_extraction_failure_visibility_tests, args.verbose)
        results.append(("Extraction-Failure Visibility",
                       "tst_extraction_failure_visibility",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_extraction_failure_visibility: {e}")
        results.append(("Extraction-Failure Visibility",
                       "tst_extraction_failure_visibility",
                       False, 0))

    # --- BUILD experiment_type + R-free Auto-fill (v118.F) ---
    try:
        from tests.tst_build_experiment_type_and_rfree import (
            run_all_tests as run_build_experiment_type_and_rfree_tests)
        success, elapsed = run_test_module(
            "tst_build_experiment_type_and_rfree",
            run_build_experiment_type_and_rfree_tests, args.verbose)
        results.append(("BUILD experiment_type + R-free Auto-fill",
                       "tst_build_experiment_type_and_rfree",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_build_experiment_type_and_rfree: {e}")
        results.append(("BUILD experiment_type + R-free Auto-fill",
                       "tst_build_experiment_type_and_rfree",
                       False, 0))

    # --- Optional Dependency Resilience (v118.G) ---
    try:
        from tests.tst_optional_dep_resilience import (
            run_all_tests as run_optional_dep_resilience_tests)
        success, elapsed = run_test_module(
            "tst_optional_dep_resilience",
            run_optional_dep_resilience_tests, args.verbose)
        results.append(("Optional Dependency Resilience",
                       "tst_optional_dep_resilience",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_optional_dep_resilience: {e}")
        results.append(("Optional Dependency Resilience",
                       "tst_optional_dep_resilience",
                       False, 0))

    # --- Density Modification Experiment-Type Reprints (v118.9) ---
    try:
        from tests.tst_density_modify_experiment_type import (
            run_all_tests as run_denmod_reprints_tests)
        success, elapsed = run_test_module(
            "tst_density_modify_experiment_type",
            run_denmod_reprints_tests, args.verbose)
        results.append(("Density Modify Experiment Type Reprints",
                       "tst_density_modify_experiment_type",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_density_modify_experiment_type: {e}")
        results.append(("Density Modify Experiment Type Reprints",
                       "tst_density_modify_experiment_type",
                       False, 0))

    # --- H18.1 PHIL round-trip for original_files_for_directives ---
    # Catches the deploy gap where master_params in ai_agent.py is
    # missing a PHIL definition that an assignment site assumes
    # exists.  See tst_h18_1_phil_roundtrip.py module docstring.
    try:
        from tests.tst_h18_1_phil_roundtrip import (
            run_all_tests as run_h18_1_phil_tests)
        success, elapsed = run_test_module(
            "tst_h18_1_phil_roundtrip",
            run_h18_1_phil_tests, args.verbose)
        results.append(("H18.1 PHIL Round-Trip",
                       "tst_h18_1_phil_roundtrip",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_h18_1_phil_roundtrip: {e}")
        results.append(("H18.1 PHIL Round-Trip",
                       "tst_h18_1_phil_roundtrip",
                       False, 0))

    # --- Settings List Coercion (v118.10) ---
    try:
        from tests.tst_settings_list_coercion import (
            run_all_tests as run_list_coercion_tests)
        success, elapsed = run_test_module(
            "tst_settings_list_coercion",
            run_list_coercion_tests, args.verbose)
        results.append(("Settings List Coercion",
                       "tst_settings_list_coercion",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import "
              f"tst_settings_list_coercion: {e}")
        results.append(("Settings List Coercion",
                       "tst_settings_list_coercion",
                       False, 0))

    # --- Environment Dependency Check (v118.6.7) ---
    # Note: this checks the PYTHON ENVIRONMENT, not code behavior.  It
    # will FAIL on any environment that doesn't have the full ai_agent
    # dependency set installed.  Mark it as optional in the registry so
    # contributors with partial envs aren't blocked.
    try:
        from tests.tst_dependencies import (
            run_all_tests as run_dependencies_tests)
        success, elapsed = run_test_module(
            "tst_dependencies",
            run_dependencies_tests, args.verbose)
        results.append(("Environment Dependency Check",
                       "tst_dependencies",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_dependencies: {e}")
        results.append(("Environment Dependency Check",
                       "tst_dependencies",
                       False, 0))

    # --- Default Model Centralization (v119.H1) ---
    # Regression suite for the centralized model-default tables in
    # core/llm.py.  Asserts that every call-site reads via
    # default_model_for_provider() and that no orphan model strings
    # appear in agent/ runtime code (AST scan).  Also baseline-
    # fingerprints the H1 table values so accidental drift surfaces
    # immediately.
    try:
        from tests.tst_default_models import (
            run_all_tests as run_default_models_tests)
        success, elapsed = run_test_module(
            "tst_default_models",
            run_default_models_tests, args.verbose)
        results.append(("Default Model Centralization",
                       "tst_default_models",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_default_models: {e}")
        results.append(("Default Model Centralization",
                       "tst_default_models",
                       False, 0))

    # --- Agent Build Info (v119.H2) ---
    # Verifies the v119.H2 server build metadata channel: VERSION
    # file, defaults fingerprint, agent_build response field, and
    # the single injection point in _build_group_args_response.
    # Includes backwards-compat assertions that existing fields are
    # unchanged.
    try:
        from tests.tst_agent_build_info import (
            run_all_tests as run_agent_build_info_tests)
        success, elapsed = run_test_module(
            "tst_agent_build_info",
            run_agent_build_info_tests, args.verbose)
        results.append(("Agent Build Info",
                       "tst_agent_build_info",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_agent_build_info: {e}")
        results.append(("Agent Build Info",
                       "tst_agent_build_info",
                       False, 0))

    # --- skip_programs Promotion (v119.H2.1) ---
    # Verifies the v119.H2.1 promotion of LLM-emitted
    # program_settings[X].skip=true to
    # workflow_preferences.skip_programs[X], wired into the top
    # of validate_directives in agent/directive_extractor.py.
    # Includes side-effect correctness, idempotency, defensive
    # bail behavior, and an end-to-end assertion that fixes the
    # skip_programs scenario in tests/llm/tst_directive_extraction.py.
    try:
        from tests.tst_skip_promotion import (
            run_all_tests as run_skip_promotion_tests)
        success, elapsed = run_test_module(
            "tst_skip_promotion",
            run_skip_promotion_tests, args.verbose)
        results.append(("skip_programs Promotion",
                       "tst_skip_promotion",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_skip_promotion: {e}")
        results.append(("skip_programs Promotion",
                       "tst_skip_promotion",
                       False, 0))

    # --- Startup Canary (v119.H3a) ---
    # Verifies the deployed server's agent_build matches what
    # canary_expected.json pins.  Catches wrong-build deploys
    # (built from wrong branch), tables in core/llm.py edited
    # without a VERSION bump (fingerprint drift), and operator
    # drift between VERSION and canary_expected.json.
    try:
        from tests.tst_canary import (
            run_all_tests as run_canary_tests)
        success, elapsed = run_test_module(
            "tst_canary",
            run_canary_tests, args.verbose)
        results.append(("Startup Canary",
                       "tst_canary",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_canary: {e}")
        results.append(("Startup Canary",
                       "tst_canary",
                       False, 0))

    # --- Preprocessor Metrics Logger / Step 1F (v119.H4) ---
    # Tests the regex-based file detector module used by the
    # Step 1F metric block in run_advice_preprocessing.
    # Telemetry-gathering for Phase 2B trigger evaluation:
    # records (llm_files, regex_files) tuples in production
    # logs to measure regex_recall = |LLM \u2229 regex|/|LLM|.
    # Scanner is stdlib-only — all 16 tests run in sandbox.
    try:
        from tests.tst_raw_advice_scanner import (
            run_all_tests as run_h4_tests)
        success, elapsed = run_test_module(
            "tst_raw_advice_scanner",
            run_h4_tests, args.verbose)
        results.append(("Preprocessor Metrics",
                       "tst_raw_advice_scanner",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_raw_advice_scanner: {e}")
        results.append(("Preprocessor Metrics",
                       "tst_raw_advice_scanner",
                       False, 0))

    # --- Diagnostic Messages / Relay Channel (v119.H5 + H5.1) ---
    # Verifies the _emit_marker helper and the diagnostic_messages
    # field in run_advice_preprocessing's result.  Closes the
    # "server-side stderr markers invisible to client" gap that
    # was diagnosed during H4 deployment.
    #
    # 26 tests organized into 6 sections:
    #   §A helper unit tests (sandbox + PHENIX)
    #   §B field-in-return tests (PHENIX only)
    #   §C backward compatibility, including corrupted payloads
    #      (sandbox + PHENIX)
    #   §D uniform client re-emit (sandbox + PHENIX, paraphrased ref impl)
    #   §E production encode/decode roundtrip (PHENIX only — exercises
    #      get_results_from_all (Total Init seed) and
    #      Program.get_results_as_JSON directly against the same JSON
    #      wire format _process_server_success decodes)
    #   §F extended markers (v119.H5.1, PHENIX only — verifies
    #      [ADVICE_PREPROCESSING_FAILED], [DIRECTIVE_EXTRACTION_FAILED],
    #      and [FAILURE_DIAGNOSIS_FAILED] markers via deterministic
    #      monkey-patching; also pins Total Init for the new
    #      run_directive_extraction and run_failure_diagnosis paths)
    try:
        from tests.tst_diagnostic_messages import (
            run_all_tests as run_h5_tests)
        success, elapsed = run_test_module(
            "tst_diagnostic_messages",
            run_h5_tests, args.verbose)
        results.append(("Diagnostic Messages",
                       "tst_diagnostic_messages",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_diagnostic_messages: {e}")
        results.append(("Diagnostic Messages",
                       "tst_diagnostic_messages",
                       False, 0))

    # 56. Directive Validation (K_H5_1_1, v119.H5.1.1) -- merge_directives
    #     _corrected_from sidecar protection (Item 3) and validate_directives
    #     boolean list-wrap defense (Item 4).  §2.1 v118.9 leftovers cleanup.
    try:
        from tests.tst_directive_validation import (
            run_all_tests as run_h5_1_1_tests)
        success, elapsed = run_test_module(
            "tst_directive_validation",
            run_h5_1_1_tests, args.verbose)
        results.append(("Directive Validation",
                       "tst_directive_validation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_directive_validation: {e}")
        results.append(("Directive Validation",
                       "tst_directive_validation",
                       False, 0))

    # 57. Planning Framework (K_H6, v119.H6) -- Phase 2A planning suite
    #     framework infrastructure: is_stop_intent, validate_planning_state,
    #     make_planning_run_fn factory, and call_planning_llm raw_output
    #     preservation on parse_intent_json failures.
    try:
        from tests.tst_planning_framework import (
            run_all_tests as run_h6_tests)
        success, elapsed = run_test_module(
            "tst_planning_framework",
            run_h6_tests, args.verbose)
        results.append(("Planning Framework",
                       "tst_planning_framework",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_planning_framework: {e}")
        results.append(("Planning Framework",
                       "tst_planning_framework",
                       False, 0))

    # 58. Phase 2B Activation (K_H7, v119.H7) -- scanner-first file
    #     extraction.  Tests the building blocks composed by
    #     run_advice_preprocessing post-H7: scanner contract
    #     (type/sorting/dedup invariants), scanner-first-with-fallback
    #     decision matrix, recall metrics with zero-division guards,
    #     and [STEP_1F] telemetry format regression.
    try:
        from tests.tst_phase2b_activation import (
            run_all_tests as run_h7_tests)
        success, elapsed = run_test_module(
            "tst_phase2b_activation",
            run_h7_tests, args.verbose)
        results.append(("Phase 2B Activation",
                       "tst_phase2b_activation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_phase2b_activation: {e}")
        results.append(("Phase 2B Activation",
                       "tst_phase2b_activation",
                       False, 0))

    # 59. Solve action keyword cleanup (K_H14_ITEM_1, v119.H14) --
    #     Phaser false-positive surfaced by run_39_openai batch
    #     analysis.  _ACTION_TABLE["solve"] entry treated the goal
    #     phrase "solve the structure" as a synonym for "do MR with
    #     phaser".  When combined with another action + stop, this
    #     forced phaser into the workflow on SAD/MAD datasets where
    #     autosol is the correct method.
    try:
        from tests.tst_solve_action_keywords import (
            run_all_tests as run_h14_item1_tests)
        success, elapsed = run_test_module(
            "tst_solve_action_keywords",
            run_h14_item1_tests, args.verbose)
        results.append(("Solve Action Keywords (H14)",
                       "tst_solve_action_keywords",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_solve_action_keywords: {e}")
        results.append(("Solve Action Keywords (H14)",
                       "tst_solve_action_keywords",
                       False, 0))

    # 60. STEP_1F single-emit (K_H14_ITEM_2, v119.H14) --
    #     Duplicate diagnostic_messages relay surfaced by
    #     run_39_openai batch analysis: 60.6% of runs showed adjacent
    #     duplicate [STEP_1F] preprocessing_metrics lines.  Two relay
    #     sites both wrote to stderr; the client-side duplicate
    #     (programs/ai_agent.py:8087-8103) was removed; the dispatcher
    #     (programs/ai_analysis.py) is the single uniform site.
    try:
        from tests.tst_step1f_single_emit import (
            run_all_tests as run_h14_item2_tests)
        success, elapsed = run_test_module(
            "tst_step1f_single_emit",
            run_h14_item2_tests, args.verbose)
        results.append(("STEP_1F Single-Emit (H14)",
                       "tst_step1f_single_emit",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_step1f_single_emit: {e}")
        results.append(("STEP_1F Single-Emit (H14)",
                       "tst_step1f_single_emit",
                       False, 0))

    # 61. Space-group validation extensions (K_H14_ITEM_3, v119.H14) --
    #     Extends _SYMMETRY_SENTINELS with the "Not explicitly
    #     mentioned" family (including truncated LLM-output forms
    #     like "Not explicitly mentio").  Adds a positive
    #     Hermann-Mauguin shape check (_looks_like_space_group) that
    #     catches prose phrases ("Solve the structure") that the
    #     pre-H14 negative checks let through.
    try:
        from tests.tst_space_group_validation import (
            run_all_tests as run_h14_item3_tests)
        success, elapsed = run_test_module(
            "tst_space_group_validation",
            run_h14_item3_tests, args.verbose)
        results.append(("Space-Group Validation (H14)",
                       "tst_space_group_validation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_space_group_validation: {e}")
        results.append(("Space-Group Validation (H14)",
                       "tst_space_group_validation",
                       False, 0))

    # --- Simple-Extractor Validation Closure Tests (v119.H14.1) ---
    # v119.H14.1: extract_directives_simple now ends with
    # validate_directives, closing the ollama-fallback bypass that
    # Tom's 2026-05-26 xtriage run surfaced (space_group="Not
    # explicitly mentio" reached the directives despite H14 Item 3
    # being installed).  Also fixes a latent VALID_STOP_CONDITIONS
    # gap (start_with_program was set by _resolve_after_program and
    # consumed downstream but wasn't in VALID_STOP_CONDITIONS).
    try:
        from tests.tst_simple_extractor_validation import (
            run_all_tests as run_h14_1_tests)
        success, elapsed = run_test_module(
            "tst_simple_extractor_validation",
            run_h14_1_tests, args.verbose)
        results.append(("Simple Extractor Validation (H14.1)",
                       "tst_simple_extractor_validation",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_simple_extractor_validation: {e}")
        results.append(("Simple Extractor Validation (H14.1)",
                       "tst_simple_extractor_validation",
                       False, 0))

    # --- predict_and_build no-NCS-input config fix (v119.H14.2) ---
    # v119.H14.2: programs.yaml fix — predict_and_build does not
    # accept an NCS-file PHIL parameter.  Pre-H14.2 the yaml
    # declared `flag: "map_model.ncs_file="` for predict_and_build,
    # which the agent dutifully emitted, and the command failed with
    # "Some PHIL parameters are not recognized."  H14.2 removes the
    # entry (Tom's 2026-05-26 1029B-sad ollama-no-PDB run surfaced).
    try:
        from tests.tst_predict_and_build_no_ncs import (
            run_all_tests as run_h14_2_tests)
        success, elapsed = run_test_module(
            "tst_predict_and_build_no_ncs",
            run_h14_2_tests, args.verbose)
        results.append(("predict_and_build no-NCS (H14.2)",
                       "tst_predict_and_build_no_ncs",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_predict_and_build_no_ncs: {e}")
        results.append(("predict_and_build no-NCS (H14.2)",
                       "tst_predict_and_build_no_ncs",
                       False, 0))

    # --- Plan catch-up failure semantics (v119.H15 Item 1) ---
    # H15 Item 1: StructurePlan.advance() now distinguishes
    # criteria-met-COMPLETE vs criteria-unmet-FAILED when called
    # from record_stage_cycle's catch-up path.  Catches Tom's
    # bromodomain bug where final_refinement was silently marked
    # complete despite R-free=0.272 > target=0.25.
    try:
        from tests.tst_plan_catchup_failure_semantics import (
            run_all_tests as run_h15_1_tests)
        success, elapsed = run_test_module(
            "tst_plan_catchup_failure_semantics",
            run_h15_1_tests, args.verbose)
        results.append(("Plan catch-up failure semantics (H15 Item 1)",
                       "tst_plan_catchup_failure_semantics",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_plan_catchup_failure_semantics: {e}")
        results.append(("Plan catch-up failure semantics (H15 Item 1)",
                       "tst_plan_catchup_failure_semantics",
                       False, 0))

    # --- Resume reopen targeted stages (v119.H15 Item 2) ---
    # H15 Item 2: reopen_stages_for_directives() — when the user
    # resumes with new advice, walk directives.program_settings and
    # reopen ONLY the LATEST completed stage matching each program.
    # Single-stage targeted reopen (per Gemini critique of cascade).
    try:
        from tests.tst_resume_reopen_stages import (
            run_all_tests as run_h15_2_tests)
        success, elapsed = run_test_module(
            "tst_resume_reopen_stages",
            run_h15_2_tests, args.verbose)
        results.append(("Resume reopen stages (H15 Item 2)",
                       "tst_resume_reopen_stages",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_resume_reopen_stages: {e}")
        results.append(("Resume reopen stages (H15 Item 2)",
                       "tst_resume_reopen_stages",
                       False, 0))

    # --- Reasoning/command divergence detection (v119.H15 Item 3) ---
    # H15 Item 3: cross_check_reasoning_vs_command() — detect-only
    # telemetry that fires when the LLM's reasoning paragraph
    # references a categorical file ("MTZ from last refinement")
    # but the command uses a different file.  Blocks ONLY when the
    # divergent file doesn't exist on disk (FileNotFoundError).
    try:
        from tests.tst_reasoning_command_divergence import (
            run_all_tests as run_h15_3_tests)
        success, elapsed = run_test_module(
            "tst_reasoning_command_divergence",
            run_h15_3_tests, args.verbose)
        results.append(("Reasoning/command divergence (H15 Item 3)",
                       "tst_reasoning_command_divergence",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_reasoning_command_divergence: {e}")
        results.append(("Reasoning/command divergence (H15 Item 3)",
                       "tst_reasoning_command_divergence",
                       False, 0))

    # --- obs_labels auto-fill for multi-array MTZ (v119.H16) ---
    # H16: cctbx-based MTZ inspection + per-program label-selection
    # policy + auto_fill_obs_labels invariant.  Targets 88 TIER-1
    # "Multiple equally suitable arrays" failures in AF_exoV_MRSAD
    # and lysozyme-MRSAD tutorials.  Three programs covered:
    # refine, phaser, autosol.
    try:
        from tests.tst_obs_labels_auto_fill import (
            run_all_tests as run_h16_tests)
        success, elapsed = run_test_module(
            "tst_obs_labels_auto_fill",
            run_h16_tests, args.verbose)
        results.append(("obs_labels auto-fill (H16)",
                       "tst_obs_labels_auto_fill",
                       success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_obs_labels_auto_fill: {e}")
        results.append(("obs_labels auto-fill (H16)",
                       "tst_obs_labels_auto_fill",
                       False, 0))

    # --- R-free-flags-missing auto-recovery (force_retry) ---
    try:
        from tests.tst_rfree_flags_missing import run_all_tests as run_rfree_missing_tests
        success, elapsed = run_test_module(
            "tst_rfree_flags_missing", run_rfree_missing_tests, args.verbose)
        results.append(("rfree flags missing recovery", "tst_rfree_flags_missing", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_rfree_flags_missing: {e}")
        results.append(("rfree flags missing recovery", "tst_rfree_flags_missing", False, 0))

    # --- classify_mtz_type _001 suffix fix (bare *_001.mtz is data, not map) ---
    try:
        from tests.tst_classify_mtz_type_001 import run_all_tests as run_classify_001_tests
        success, elapsed = run_test_module(
            "tst_classify_mtz_type_001", run_classify_001_tests, args.verbose)
        results.append(("classify_mtz_type _001 fix", "tst_classify_mtz_type_001", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_classify_mtz_type_001: {e}")
        results.append(("classify_mtz_type _001 fix", "tst_classify_mtz_type_001", False, 0))

    # --- phenix.refine Final-vs-Start R-free extraction fix ---
    try:
        from tests.tst_refine_final_rfree import run_all_tests as run_refine_final_rfree_tests
        success, elapsed = run_test_module(
            "tst_refine_final_rfree", run_refine_final_rfree_tests, args.verbose)
        results.append(("refine final rfree", "tst_refine_final_rfree", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_refine_final_rfree: {e}")
        results.append(("refine final rfree", "tst_refine_final_rfree", False, 0))

    # --- spurious stop_after_requested strip (open-ended advice != stop) ---
    try:
        from tests.tst_spurious_stop_after_requested import run_all_tests as run_spurious_stop_tests
        success, elapsed = run_test_module(
            "tst_spurious_stop_after_requested", run_spurious_stop_tests, args.verbose)
        results.append(("spurious stop strip", "tst_spurious_stop_after_requested", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_spurious_stop_after_requested: {e}")
        results.append(("spurious stop strip", "tst_spurious_stop_after_requested", False, 0))

    # --- advice-driven rebuild template selection ---
    try:
        from tests.tst_rebuild_template_selection import run_all_tests as run_rebuild_template_tests
        success, elapsed = run_test_module(
            "tst_rebuild_template_selection", run_rebuild_template_tests, args.verbose)
        results.append(("rebuild template selection", "tst_rebuild_template_selection", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_rebuild_template_selection: {e}")
        results.append(("rebuild template selection", "tst_rebuild_template_selection", False, 0))

    # --- placement-skip preserves explicit rebuild ---
    try:
        from tests.tst_placement_skip_rebuild import run_all_tests as run_placement_skip_tests
        success, elapsed = run_test_module(
            "tst_placement_skip_rebuild", run_placement_skip_tests, args.verbose)
        results.append(("placement-skip rebuild", "tst_placement_skip_rebuild", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_placement_skip_rebuild: {e}")
        results.append(("placement-skip rebuild", "tst_placement_skip_rebuild", False, 0))

    # --- autobuild sequence waiver ---
    try:
        from tests.tst_autobuild_sequence_waiver import run_all_tests as run_autobuild_seq_tests
        success, elapsed = run_test_module(
            "tst_autobuild_sequence_waiver", run_autobuild_seq_tests, args.verbose)
        results.append(("autobuild sequence waiver", "tst_autobuild_sequence_waiver", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_autobuild_sequence_waiver: {e}")
        results.append(("autobuild sequence waiver", "tst_autobuild_sequence_waiver", False, 0))

    # --- Option 2a: reactive-deviation hold ---
    try:
        from tests.tst_reactive_deviation_hold import run_all_tests as run_reactive_hold_tests
        success, elapsed = run_test_module(
            "tst_reactive_deviation_hold", run_reactive_hold_tests, args.verbose)
        results.append(("reactive deviation hold", "tst_reactive_deviation_hold", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_reactive_deviation_hold: {e}")
        results.append(("reactive deviation hold", "tst_reactive_deviation_hold", False, 0))

    # --- Option 2a: plan lead-program offer ---
    try:
        from tests.tst_plan_lead_program_offer import run_all_tests as run_lead_offer_tests
        success, elapsed = run_test_module(
            "tst_plan_lead_program_offer", run_lead_offer_tests, args.verbose)
        results.append(("plan lead program offer", "tst_plan_lead_program_offer", success, elapsed))
    except ImportError as e:
        print(f"\u26a0\ufe0f  Could not import tst_plan_lead_program_offer: {e}")
        results.append(("plan lead program offer", "tst_plan_lead_program_offer", False, 0))

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
