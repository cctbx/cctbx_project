"""
Unit tests for history analysis and anomalous workflow support.

These tests verify:
1. Analysis/metrics key handling in transport (the key mismatch fix)
2. Anomalous metrics extraction from session results
3. History info extraction from workflow_state context
4. Autosol availability when anomalous signal is detected

Run with: python tests/tst_history_analysis.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import tempfile

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    run_tests_with_fail_fast
)


# =============================================================================
# TRANSPORT: ANALYSIS/METRICS KEY HANDLING
# =============================================================================

def test_build_request_v2_handles_analysis_key():
    """build_request_v2 extracts metrics from 'analysis' key (session format)."""
    print("Test: build_request_v2_handles_analysis_key")

    try:
        from agent.api_client import build_request_v2

        # History with 'analysis' key (as returned by session.get_history_for_agent)
        history = [
            {
                "cycle_number": 1,
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS",
                "analysis": {  # <-- 'analysis' key from session
                    "resolution": 2.5,
                    "has_anomalous": True,
                    "anomalous_measurability": 0.15,
                },
                "output_files": [],
            }
        ]

        request = build_request_v2(
            files=["data.mtz"],
            cycle_number=2,
            history=history
        )

        # Check that metrics were extracted from 'analysis' key
        normalized_history = request["history"]
        assert_equal(len(normalized_history), 1)
        metrics = normalized_history[0]["metrics"]

        assert_equal(metrics.get("resolution"), 2.5)
        assert_equal(metrics.get("has_anomalous"), True)
        assert_equal(metrics.get("anomalous_measurability"), 0.15)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


def test_build_request_v2_handles_metrics_key():
    """build_request_v2 also handles 'metrics' key (fallback)."""
    print("Test: build_request_v2_handles_metrics_key")

    try:
        from agent.api_client import build_request_v2

        # History with 'metrics' key directly
        history = [
            {
                "cycle_number": 1,
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS",
                "metrics": {  # <-- 'metrics' key directly
                    "resolution": 2.5,
                    "has_anomalous": True,
                },
                "output_files": [],
            }
        ]

        request = build_request_v2(
            files=["data.mtz"],
            cycle_number=2,
            history=history
        )

        normalized_history = request["history"]
        metrics = normalized_history[0]["metrics"]

        assert_equal(metrics.get("resolution"), 2.5)
        assert_equal(metrics.get("has_anomalous"), True)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


# =============================================================================
# WORKFLOW_STATE: ANALYZE HISTORY
# =============================================================================

def test_analyze_history_extracts_anomalous_from_analysis():
    """_analyze_history extracts anomalous info from 'analysis' key."""
    print("Test: analyze_history_extracts_anomalous_from_analysis")

    try:
        from agent.workflow_state import _analyze_history

        history = [
            {
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS",
                "analysis": {
                    "has_anomalous": True,
                    "anomalous_measurability": 0.15,
                    "anomalous_resolution": 3.5,
                }
            }
        ]

        info = _analyze_history(history)

        assert_true(info.get("has_anomalous"), "has_anomalous should be True")
        assert_true(info.get("strong_anomalous"), "strong_anomalous should be True (measurability > 0.10)")
        assert_equal(info.get("anomalous_measurability"), 0.15)
        assert_equal(info.get("anomalous_resolution"), 3.5)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


def test_analyze_history_extracts_anomalous_from_metrics():
    """_analyze_history extracts anomalous info from 'metrics' key (after transport)."""
    print("Test: analyze_history_extracts_anomalous_from_metrics")

    try:
        from agent.workflow_state import _analyze_history

        # After transport, the key is 'metrics' not 'analysis'
        history = [
            {
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS",
                "metrics": {  # <-- 'metrics' key after transport
                    "has_anomalous": True,
                    "anomalous_measurability": 0.20,
                }
            }
        ]

        info = _analyze_history(history)

        assert_true(info.get("has_anomalous"), "has_anomalous should be True")
        assert_true(info.get("strong_anomalous"), "strong_anomalous should be True")
        assert_equal(info.get("anomalous_measurability"), 0.20)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


def test_analyze_history_extracts_ncs():
    """_analyze_history extracts NCS info from analysis."""
    print("Test: analyze_history_extracts_ncs")

    try:
        from agent.workflow_state import _analyze_history

        history = [
            {
                "program": "phenix.map_symmetry",
                "command": "phenix.map_symmetry map.mrc",
                "result": "SUCCESS",
                "analysis": {
                    "ncs_found": True,
                }
            }
        ]

        info = _analyze_history(history)
        assert_true(info.get("has_ncs"), "has_ncs should be True")

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


def test_analyze_history_weak_anomalous_not_strong():
    """Weak anomalous signal (measurability < 0.10) doesn't set strong_anomalous."""
    print("Test: analyze_history_weak_anomalous_not_strong")

    try:
        from agent.workflow_state import _analyze_history

        history = [
            {
                "program": "phenix.xtriage",
                "analysis": {
                    "has_anomalous": True,
                    "anomalous_measurability": 0.05,  # Weak signal
                }
            }
        ]

        info = _analyze_history(history)

        assert_true(info.get("has_anomalous"), "has_anomalous should be True")
        assert_false(info.get("strong_anomalous"), "strong_anomalous should be False (weak signal)")

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


# =============================================================================
# SESSION: METRICS EXTRACTION
# =============================================================================

def test_session_extracts_anomalous_from_result():
    """Session._extract_metrics_from_result extracts anomalous metrics."""
    print("Test: session_extracts_anomalous_from_result")

    # Mock the YAML loading to avoid libtbx dependency
    import sys
    import types

    # Create mock for libtbx
    if 'libtbx' not in sys.modules:
        libtbx = types.ModuleType('libtbx')
        libtbx.langchain = types.ModuleType('libtbx.langchain')
        libtbx.langchain.knowledge = types.ModuleType('libtbx.langchain.knowledge')
        libtbx.langchain.knowledge.metric_patterns = types.ModuleType('libtbx.langchain.knowledge.metric_patterns')
        libtbx.langchain.knowledge.metric_patterns.extract_metrics_for_program = lambda x, y: {}
        sys.modules['libtbx'] = libtbx
        sys.modules['libtbx.langchain'] = libtbx.langchain
        sys.modules['libtbx.langchain.knowledge'] = libtbx.langchain.knowledge
        sys.modules['libtbx.langchain.knowledge.metric_patterns'] = libtbx.langchain.knowledge.metric_patterns

    from agent.session import AgentSession

    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Simulate xtriage result with anomalous metrics
    result_text = """SUCCESS: Command completed without errors

**************************************************
FINAL QUALITY METRICS REPORT:
--------------------------------------------------
Anomalous Measurability: 0.1500
Anomalous Resolution: 3.50
Has Anomalous: True
Resolution: 2.50
**************************************************
"""

    metrics = session._extract_metrics_from_result(result_text, "phenix.xtriage")

    # Check anomalous metrics were extracted
    assert_equal(metrics.get("anomalous_measurability"), 0.15,
                 f"Expected anomalous_measurability=0.15, got {metrics.get('anomalous_measurability')}")
    assert_equal(metrics.get("anomalous_resolution"), 3.50,
                 f"Expected anomalous_resolution=3.50, got {metrics.get('anomalous_resolution')}")
    assert_true(metrics.get("has_anomalous"),
                f"Expected has_anomalous=True, got {metrics.get('has_anomalous')}")

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


# =============================================================================
# GRAPH_NODES: HISTORY INFO EXTRACTION
# =============================================================================

def test_extract_history_info_uses_workflow_context():
    """_extract_history_info pulls from workflow_state context, not history."""
    print("Test: extract_history_info_uses_workflow_context")

    try:
        from agent.graph_nodes import _extract_history_info

        # State with pre-computed workflow_state context
        state = {
            "cycle_number": 3,
            "workflow_state": {
                "context": {
                    "refine_count": 2,
                    "rsr_count": 1,
                    "twin_law": "-k,-h,-l",
                    "has_ncs": True,
                }
            },
            # Old history format - should be ignored now
            "history": [
                {"program": "phenix.refine"},
                {"program": "phenix.refine"},
            ]
        }

        info = _extract_history_info(state)

        # Should get values from context, not recompute from history
        assert_equal(info["refine_count"], 2)
        assert_equal(info["rsr_count"], 1)
        assert_equal(info["twin_law"], "-k,-h,-l")
        assert_true(info["has_ncs"])
        assert_equal(info["cycle_number"], 3)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


def test_extract_history_info_handles_missing_context():
    """_extract_history_info returns defaults when context is missing."""
    print("Test: extract_history_info_handles_missing_context")

    try:
        from agent.graph_nodes import _extract_history_info

        # State without workflow_state
        state = {
            "cycle_number": 1,
        }

        info = _extract_history_info(state)

        assert_equal(info["refine_count"], 0)
        assert_equal(info["rsr_count"], 0)
        assert_equal(info["twin_law"], None)
        assert_false(info["has_ncs"])
        assert_equal(info["cycle_number"], 1)

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
        else:
            raise


# =============================================================================
# INTEGRATION: ANOMALOUS WORKFLOW
# =============================================================================

def test_anomalous_history_enables_autosol():
    """History with anomalous signal makes autosol available in obtain_model phase."""
    print("Test: anomalous_history_enables_autosol")

    try:
        from agent.workflow_state import _analyze_history
        from agent.workflow_engine import WorkflowEngine

        # Simulate xtriage detected anomalous signal
        history = [
            {
                "program": "phenix.xtriage",
                "command": "phenix.xtriage data.mtz",
                "result": "SUCCESS",
                "analysis": {
                    "resolution": 2.5,
                    "has_anomalous": True,
                    "anomalous_measurability": 0.15,
                }
            }
        ]

        # Analyze history
        history_info = _analyze_history(history)

        # Check anomalous was detected
        assert_true(history_info.get("has_anomalous"), "has_anomalous should be True")
        assert_true(history_info.get("strong_anomalous"), "strong_anomalous should be True")

        # Build context with files
        files = {
            "data_mtz": ["data.mtz"],
            "sequence": ["seq.fa"],
        }

        engine = WorkflowEngine()
        context = engine.build_context(files, history_info, analysis={})

        # Check context has anomalous flag
        assert_true(context.get("has_anomalous"), "Context should have has_anomalous=True")

        # Get workflow state
        state = engine.get_workflow_state("xray", files, history_info, analysis={})

        # Check autosol is in valid programs
        valid_programs = state.get("valid_programs", [])
        assert_in("phenix.autosol", valid_programs,
                  f"autosol should be valid. Valid programs: {valid_programs}")

        print("  PASSED")

    except (ImportError, ModuleNotFoundError) as e:
        if "libtbx" in str(e):
            print(f"  SKIPPED (requires libtbx)")
        else:
            raise


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
