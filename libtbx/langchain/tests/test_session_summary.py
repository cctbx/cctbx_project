"""
Tests for Agent Session Summary generation.

Run with: python tests/test_session_summary.py
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import tempfile
import types

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock langchain_core only if not available
try:
    from langchain_core.prompts import PromptTemplate as _PromptTemplate  # noqa: F401
    del _PromptTemplate  # Not used directly, just checking availability
except ImportError:
    langchain_core = types.ModuleType('langchain_core')
    langchain_core.prompts = types.ModuleType('langchain_core.prompts')
    class MockPromptTemplate:
        def __init__(self, *args, **kwargs):
            pass
    langchain_core.prompts.PromptTemplate = MockPromptTemplate
    sys.modules['langchain_core'] = langchain_core
    sys.modules['langchain_core.prompts'] = langchain_core.prompts

# Now import session
from agent.session import AgentSession


def create_test_session():
    """Create a test session with sample data."""
    # Create temp directory
    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Set up session data
    session.data["project_advice"] = "solve the structure using data to 3 A"
    session.data["original_files"] = ["/data/7qz0.mtz", "/data/7qz0.fa"]
    session.data["experiment_type"] = "xray"
    session.data["resolution"] = 2.1

    # Add cycles
    session.data["cycles"] = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "decision": "Run phenix.xtriage to analyze data quality",
            "command": "phenix.xtriage /data/7qz0.mtz",
            "result": "SUCCESS: Command completed\n\nResolution: 2.10\nCompleteness: 98.0%\nNo Twinning Suspected: True",
            "output_files": [],
        },
        {
            "cycle_number": 2,
            "program": "phenix.predict_and_build",
            "decision": "Generate AlphaFold model from sequence",
            "command": "phenix.predict_and_build ...",
            "result": "SUCCESS: Command completed\n\nPlddt: 95.25",
            "output_files": ["predicted_model.pdb"],
        },
        {
            "cycle_number": 3,
            "program": "phenix.phaser",
            "decision": "Run molecular replacement",
            "command": "phenix.phaser ...",
            "result": "SUCCESS: Command completed\n\nTFZ: 45.2\nSpace Group: P 32 2 1",
            "output_files": ["PHASER.1.pdb", "PHASER.1.mtz"],
        },
        {
            "cycle_number": 4,
            "program": "phenix.refine",
            "decision": "Initial refinement",
            "command": "phenix.refine ...",
            "result": "SUCCESS: Command completed\n\nR Free: 0.3151\nR Work: 0.2368\nClashscore: 9.28",
            "output_files": ["refine_001_001.pdb", "refine_001_001.mtz"],
        },
        {
            "cycle_number": 5,
            "program": "phenix.refine",
            "decision": "Continue refinement",
            "command": "phenix.refine ...",
            "result": "SUCCESS: Command completed\n\nR Free: 0.2812\nR Work: 0.2507\nClashscore: 6.18\nBonds Rmsd: 0.010\nAngles Rmsd: 1.19",
            "output_files": ["refine_002_001.pdb", "refine_002_001.mtz"],
        },
    ]

    return session, temp_dir


def test_extract_summary_data():
    """Test that summary data extraction works."""
    print("Test: extract_summary_data")

    session, temp_dir = create_test_session()

    data = session._extract_summary_data()

    # Check basic fields
    assert data["experiment_type"] == "xray", f"Expected xray, got {data['experiment_type']}"
    assert data["user_advice"] == "solve the structure using data to 3 A"
    assert len(data["original_files"]) == 2
    assert data["total_cycles"] == 5
    assert data["successful_cycles"] == 5

    # Check workflow path
    assert "AlphaFold" in data["workflow_path"] or "molecular replacement" in data["workflow_path"].lower()

    # Check steps
    assert len(data["steps"]) == 5
    assert data["steps"][0]["program"] == "phenix.xtriage"
    assert data["steps"][0]["success"] == True

    # Check final metrics
    assert "r_free" in data["final_metrics"]
    assert abs(data["final_metrics"]["r_free"] - 0.2812) < 0.001

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_input_quality_extraction():
    """Test input quality metrics extraction."""
    print("Test: input_quality_extraction")

    session, temp_dir = create_test_session()

    data = session._extract_summary_data()

    # Check input quality from xtriage
    iq = data.get("input_quality", {})
    assert "resolution" in iq, "Resolution not extracted"
    assert abs(iq["resolution"] - 2.10) < 0.01
    assert "completeness" in iq, "Completeness not extracted"
    assert abs(iq["completeness"] - 98.0) < 0.1
    assert "twinning" in iq, "Twinning not extracted"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_workflow_path_detection():
    """Test workflow path detection for different scenarios."""
    print("Test: workflow_path_detection")

    session, temp_dir = create_test_session()

    # Test X-ray with AlphaFold
    cycles = session.data["cycles"]
    path = session._determine_workflow_path(cycles, "xray")
    assert "AlphaFold" in path or "molecular replacement" in path.lower()

    # Test cryo-EM
    cryoem_cycles = [
        {"program": "phenix.mtriage"},
        {"program": "phenix.dock_in_map"},
        {"program": "phenix.real_space_refine"},
    ]
    path = session._determine_workflow_path(cryoem_cycles, "cryoem")
    assert "Cryo-EM" in path

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_markdown_generation():
    """Test Markdown summary generation."""
    print("Test: markdown_generation")

    session, temp_dir = create_test_session()

    result = session.generate_agent_session_summary(include_llm_assessment=False)

    markdown = result["markdown"]

    # Check structure - header uses ## now for smaller text
    assert "## Phenix AI" in markdown  # Either "Phenix AI Run:" or "Phenix AI Tutorial:"
    assert "## Input" in markdown
    assert "## Workflow Path" in markdown
    assert "## Steps Performed" in markdown
    assert "## Final Quality" in markdown

    # Check content
    assert "7qz0.mtz" in markdown
    assert "solve the structure using data to 3 A" in markdown
    assert "xtriage" in markdown
    assert "refine" in markdown
    assert "0.2812" in markdown  # R-free

    # Check table format
    assert "| Cycle | Program | Result | Key Metric |" in markdown

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_llm_summary_input():
    """Test the concise summary for LLM assessment."""
    print("Test: llm_summary_input")

    session, temp_dir = create_test_session()

    llm_input = session.get_summary_for_llm_assessment()

    # Check key sections
    assert "EXPERIMENT TYPE:" in llm_input
    assert "INPUT FILES:" in llm_input
    assert "USER GOAL:" in llm_input
    assert "WORKFLOW PATH:" in llm_input
    assert "STEPS TAKEN:" in llm_input
    assert "FINAL METRICS:" in llm_input

    # Check content is concise
    lines = llm_input.split("\n")
    assert len(lines) < 50, f"LLM input too long: {len(lines)} lines"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_quality_assessments():
    """Test quality metric assessments."""
    print("Test: quality_assessments")

    session, temp_dir = create_test_session()

    # Test R-free assessment
    assert session._assess_r_free(0.20, 2.0) == "Good"
    assert session._assess_r_free(0.28, 2.0) == "Acceptable"
    assert session._assess_r_free(0.35, 2.0) == "Needs improvement"

    # Resolution-dependent
    assert session._assess_r_free(0.22, 1.0) == "Acceptable"  # Stricter at high res
    assert session._assess_r_free(0.32, 4.0) == "Good"  # Looser at low res

    # Test map CC assessment
    assert session._assess_map_cc(0.85) == "Good"
    assert session._assess_map_cc(0.75) == "Acceptable"
    assert session._assess_map_cc(0.60) == "Needs improvement"

    # Test clashscore assessment
    assert session._assess_clashscore(3) == "Excellent"
    assert session._assess_clashscore(8) == "Good"
    assert session._assess_clashscore(15) == "Acceptable"
    assert session._assess_clashscore(25) == "Needs improvement"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_empty_session():
    """Test handling of empty session."""
    print("Test: empty_session")

    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Should not crash on empty session
    result = session.generate_agent_session_summary(include_llm_assessment=False)

    assert result["markdown"] is not None
    assert "## Phenix AI" in result["markdown"]  # Header uses ## now
    assert result["data"]["total_cycles"] == 0

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_stop_cycle_excluded_from_count():
    """Test that STOP cycles are excluded from cycle counts.

    STOP is a termination signal, not a real program run.
    The summary should report only actual program cycles.
    """
    print("Test: stop_cycle_excluded_from_count")

    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Set up minimal session data
    session.data["project_advice"] = "test"
    session.data["original_files"] = ["/data/test.mtz"]
    session.data["experiment_type"] = "xray"

    # Add 3 real cycles + 1 STOP cycle
    session.data["cycles"] = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "command": "phenix.xtriage test.mtz",
            "result": "SUCCESS: completed",
        },
        {
            "cycle_number": 2,
            "program": "phenix.refine",
            "command": "phenix.refine ...",
            "result": "SUCCESS: R-free=0.30",
        },
        {
            "cycle_number": 3,
            "program": "phenix.refine",
            "command": "phenix.refine ...",
            "result": "FAILED: error occurred",
        },
        {
            "cycle_number": 4,
            "program": "STOP",
            "command": "STOP",
            "result": "Workflow complete",
        },
    ]

    data = session._extract_summary_data()

    # Should count only real programs (3), not STOP
    assert data["total_cycles"] == 3, \
        f"Expected 3 cycles (excluding STOP), got {data['total_cycles']}"

    # Should count only successful real programs (2)
    assert data["successful_cycles"] == 2, \
        f"Expected 2 successful cycles, got {data['successful_cycles']}"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_extract_metrics_includes_anomalous():
    """Test that _extract_metrics_from_result extracts anomalous metrics."""
    print("Test: extract_metrics_includes_anomalous")

    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Simulate xtriage result text with anomalous metrics
    result_text = """SUCCESS: Command completed without errors

**************************************************
FINAL QUALITY METRICS REPORT:
--------------------------------------------------
Anomalous Measurability: 0.1500
Anomalous Resolution: 3.50
Completeness: 98.50
Has Anomalous: True
Resolution: 2.50
**************************************************
"""

    try:
        metrics = session._extract_metrics_from_result(result_text, "phenix.xtriage")
    except Exception as e:
        # May fail without libtbx, but we can still test the regex patterns
        if "libtbx" in str(e):
            print("  SKIPPED (requires libtbx)")
            import shutil
            shutil.rmtree(temp_dir)
            return
        raise

    # Check that anomalous metrics were extracted
    assert metrics.get("anomalous_measurability") == 0.15, \
        f"Expected anomalous_measurability=0.15, got {metrics.get('anomalous_measurability')}"

    assert metrics.get("anomalous_resolution") == 3.50, \
        f"Expected anomalous_resolution=3.50, got {metrics.get('anomalous_resolution')}"

    assert metrics.get("has_anomalous") == True, \
        f"Expected has_anomalous=True, got {metrics.get('has_anomalous')}"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def test_history_preserves_anomalous_across_restart():
    """Test that anomalous metrics are available in history after restart."""
    print("Test: history_preserves_anomalous_across_restart")

    temp_dir = tempfile.mkdtemp()
    session = AgentSession(session_dir=temp_dir)

    # Simulate a session with xtriage that found anomalous signal
    session.data["cycles"] = [
        {
            "cycle_number": 1,
            "program": "phenix.xtriage",
            "command": "phenix.xtriage data.mtz",
            "result": "SUCCESS: Command completed without errors",
            "metrics": {
                "resolution": 2.5,
                "anomalous_measurability": 0.15,
                "anomalous_resolution": 3.5,
                "has_anomalous": True
            }
        }
    ]
    session.save()

    # Reload the session (simulating restart)
    session2 = AgentSession(session_dir=temp_dir)

    # Get history as the agent would see it
    history = session2.get_history_for_agent()

    # Check that history has the anomalous metrics
    assert len(history) == 1, f"Expected 1 history entry, got {len(history)}"
    analysis = history[0].get("analysis", {})

    assert analysis.get("has_anomalous") == True, \
        f"Expected has_anomalous=True in history, got {analysis.get('has_anomalous')}"

    assert analysis.get("anomalous_measurability") == 0.15, \
        f"Expected anomalous_measurability=0.15 in history, got {analysis.get('anomalous_measurability')}"

    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

    print("  PASSED")


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.test_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
