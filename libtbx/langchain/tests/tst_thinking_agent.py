"""Tests for agent/thinking_agent.py (v113 Thinking Agent).

Uses mocked LLM to test should_think(), run_think_node(),
and parse_assessment() integration.
"""

from __future__ import absolute_import, division, print_function
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.thinking_agent import (
  should_think, run_think_node, _build_thinking_context, _update_memory)
from knowledge.thinking_prompts import parse_assessment


# =========================================================================
# Helper: minimal state factory
# =========================================================================

def _make_state(**overrides):
  """Build a minimal valid state for testing."""
  state = {
    "thinking_level": "basic",
    "history": [],
    "log_text": "",
    "log_analysis": {},
    "session_info": {"experiment_type": "xray"},
    "workflow_state": {"state": "initial"},
    "cycle_number": 1,
    "user_advice": "Focus on model quality",
    "strategy_memory": {},
    "provider": "google",
    "debug_log": [],
    "stop": False,
    "stop_reason": None,
    "expert_assessment": {},
  }
  state.update(overrides)
  return state


# =========================================================================
# should_think tests
# =========================================================================

def test_should_think_first_cycle():
  print("Test: should_think_first_cycle")
  # First cycle has no program output — nothing to analyze
  assert should_think({"history": []}) == False
  print("  PASSED")


def test_should_think_after_xtriage():
  print("Test: should_think_after_xtriage")
  state = {"history": [{"program": "phenix.xtriage"}]}
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_after_phaser():
  print("Test: should_think_after_phaser")
  state = {"history": [{"program": "phenix.phaser"}]}
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_after_autosol():
  print("Test: should_think_after_autosol")
  state = {"history": [{"program": "phenix.autosol"}]}
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_after_autobuild():
  print("Test: should_think_after_autobuild")
  state = {"history": [{"program": "phenix.autobuild"}]}
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_routine_refine():
  print("Test: should_think_routine_refine")
  # molprobity is a validation program (non-strategic) — should not trigger thinking
  state = {"history": [{"program": "phenix.molprobity", "status": "OK"}]}
  assert should_think(state) == False
  print("  PASSED")


def test_should_think_after_failure():
  print("Test: should_think_after_failure")
  state = {"history": [{"program": "phenix.refine", "status": "FAILED"}]}
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_stalled_rfree():
  print("Test: should_think_stalled_rfree")
  state = {
    "history": [{"program": "phenix.refine", "status": "OK"}],
    "strategy_memory": {"r_free_history": [0.30, 0.30, 0.301]},
  }
  assert should_think(state) == True
  print("  PASSED")


def test_should_think_improving_rfree():
  print("Test: should_think_improving_rfree")
  # Use a non-strategic program so only the R-free trend is checked.
  # Improving R-free should not trigger thinking.
  state = {
    "history": [{"program": "phenix.molprobity", "status": "OK"}],
    "strategy_memory": {"r_free_history": [0.30, 0.28, 0.26]},
  }
  assert should_think(state) == False
  print("  PASSED")


# =========================================================================
# run_think_node tests
# =========================================================================

def test_disabled_passthrough():
  print("Test: disabled_passthrough")
  state = _make_state(thinking_level=None)
  result = run_think_node(state)
  assert result is state
  print("  PASSED")


def test_enabled_no_llm_graceful():
  print("Test: enabled_no_llm_graceful")
  state = _make_state(
    history=[{"program": "phenix.xtriage", "status": "OK"}])
  result = run_think_node(state)
  # Should not crash, should not stop
  assert not result.get("stop")
  # Should have debug log entry
  debug = result.get("debug_log", [])
  assert any("THINK:" in str(d) for d in debug)
  print("  PASSED")


def test_routine_step_skipped():
  print("Test: routine_step_skipped")
  # molprobity is a validation program (non-strategic category)
  # and should not trigger the thinking agent
  state = _make_state(
    history=[{"program": "phenix.molprobity", "status": "OK"}])
  result = run_think_node(state)
  debug = result.get("debug_log", [])
  assert any("Skipping" in str(d) for d in debug)
  print("  PASSED")


def test_first_cycle_skipped():
  print("Test: first_cycle_skipped")
  state = _make_state(history=[])
  result = run_think_node(state)
  debug = result.get("debug_log", [])
  assert any("Skipping" in str(d) for d in debug)
  print("  PASSED")


# =========================================================================
# _build_thinking_context tests
# =========================================================================

def test_build_context_first_cycle():
  print("Test: build_context_first_cycle")
  state = _make_state()
  ctx = _build_thinking_context(state)
  assert ctx["cycle_number"] == 1
  assert ctx["experiment_type"] == "xray"
  assert ctx["program_name"] == ""
  assert "(first cycle)" in ctx["history_summary"]
  print("  PASSED")


def test_build_context_with_history():
  print("Test: build_context_with_history")
  state = _make_state(
    cycle_number=3,
    history=[
      {"program": "phenix.xtriage", "status": "OK", "cycle_number": 1},
      {"program": "phenix.phaser", "status": "OK", "cycle_number": 2},
    ],
    log_text="TFZ=12.5 LLG=150 SOLU SET",
  )
  ctx = _build_thinking_context(state)
  assert ctx["program_name"] == "phenix.phaser"
  assert "xtriage" in ctx["history_summary"]
  assert len(ctx["log_sections"]) > 0
  print("  PASSED")


# =========================================================================
# _update_memory tests
# =========================================================================

def test_update_memory_basic():
  print("Test: update_memory_basic")
  updated = _update_memory({}, {
    "data_quality": "twinned",
    "phasing_strategy": "SAD",
    "guidance": "Try autosol with Se atoms",
  }, cycle_number=3)
  assert updated["data_quality"] == "twinned"
  assert updated["phasing_strategy"] == "SAD"
  assert len(updated["decisions"]) == 1
  assert updated["decisions"][0][0] == 3  # actual cycle number
  print("  PASSED")


def test_update_memory_preserves_existing():
  print("Test: update_memory_preserves_existing")
  existing = {
    "data_quality": "good",
    "r_free_history": [0.30, 0.28],
    "programs_run": ["phenix.refine"],
    "concerns": [],
    "decisions": [],
    "phasing_strategy": "MR",
  }
  updated = _update_memory(existing, {"concerns": ["ice rings"]})
  assert updated["data_quality"] == "good"
  assert updated["r_free_history"] == [0.30, 0.28]
  assert "ice rings" in updated["concerns"]
  print("  PASSED")


# =========================================================================
# parse_assessment integration
# =========================================================================

def test_parse_valid_json():
  print("Test: parse_valid_json")
  raw = json.dumps({
    "analysis": "Good MR solution",
    "confidence": "high",
    "action": "guide_step",
    "guidance": "Proceed to refinement",
  })
  a = parse_assessment(raw)
  assert a["analysis"] == "Good MR solution"
  assert a["action"] == "guide_step"
  assert a["confidence"] == "high"
  print("  PASSED")


def test_parse_with_fences():
  print("Test: parse_with_fences")
  raw = "```json\n" + json.dumps({
    "analysis": "Test",
    "action": "let_run",
  }) + "\n```"
  a = parse_assessment(raw)
  assert a["analysis"] == "Test"
  print("  PASSED")


def test_parse_invalid_action_corrected():
  print("Test: parse_invalid_action_corrected")
  raw = json.dumps({"action": "destroy", "analysis": "test"})
  a = parse_assessment(raw)
  assert a["action"] == "let_run"
  print("  PASSED")


def test_parse_garbage():
  print("Test: parse_garbage")
  a = parse_assessment("this is not json")
  assert a["action"] == "let_run"
  assert "this is not json" in a["analysis"]
  print("  PASSED")


def test_parse_none():
  print("Test: parse_none")
  a = parse_assessment(None)
  assert a["action"] == "let_run"
  assert a["analysis"] == ""
  print("  PASSED")


# =========================================================================

def run_tests_with_fail_fast():
  tests = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
  passed = 0
  failed = 0
  for t in tests:
    name = t.__name__
    try:
      t()
      passed += 1
      print("  PASS: %s" % name)
    except Exception as e:
      failed += 1
      print("  FAIL: %s — %s" % (name, e))
      break
  print("\n%d passed, %d failed" % (passed, failed))
  return failed == 0


def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  success = run_tests_with_fail_fast()
  sys.exit(0 if success else 1)
