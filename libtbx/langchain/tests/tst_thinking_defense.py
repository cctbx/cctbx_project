#!/usr/bin/env python
"""
Defensive tests for Thinking Agent integration.

These tests lock down the EXISTING behavior of every component
that Phase A will touch. Run before AND after each Phase A step
to verify nothing breaks.

Two categories:
  STANDALONE (A-F, I): Run anywhere, no PHENIX needed.
    Tests graph_state, routing contracts, state serialization,
    user_advice enrichment, session_info extensibility.

  PHENIX-REQUIRED (G, H): Only run in PHENIX environment.
    Tests perceive/plan pipeline, node pass-through behavior.

Run with: python tests/tst_thinking_defense.py
"""

from __future__ import absolute_import, division, print_function
import sys
import os
import copy
import json
import inspect

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
  assert_equal, assert_true, assert_false, assert_in,
  assert_none, assert_is_instance,
  run_tests_with_fail_fast
)

# graph_state imports cleanly (no libtbx dependency)
from agent.graph_state import AgentState, create_initial_state

# graph.py and graph_nodes.py need langgraph/libtbx — guard imports
HAS_GRAPH_NODES = False
try:
  from agent.graph import (
    route_after_perceive, route_after_validate)
  import agent.graph_nodes
  assert agent.graph_nodes is not None
  HAS_GRAPH_NODES = True
except ImportError:
  # Define routing function contracts locally for standalone testing.
  # These mirror the EXACT logic in graph.py — if graph.py changes,
  # these tests catch the mismatch via source verification (test I03).
  def route_after_perceive(state):
    if state.get("stop"):
      return "output"
    return "think"

  def route_after_validate(state):
    if state.get("stop"):
      return "output"
    if state.get("validation_error"):
      attempt = state.get("attempt_number", 0)
      if attempt < 3:
        return "plan"
      else:
        return "fallback"
    return "output"


# =========================================================================
# A. Graph State: fields, defaults, create_initial_state
# =========================================================================

def test_A01_create_initial_state_returns_dict():
  """create_initial_state returns a plain dict."""
  print("Test: A01_create_initial_state_returns_dict")
  state = create_initial_state(available_files=["data.mtz"])
  assert_is_instance(state, dict)
  print("  PASSED")


def test_A02_all_required_fields_present():
  """Every expected field is present in created state."""
  print("Test: A02_all_required_fields_present")
  state = create_initial_state(available_files=["data.mtz"])
  required = [
    "history", "available_files", "cycle_number", "max_cycles",
    "user_advice", "provider", "session_resolution", "session_info",
    "directives", "bad_inject_params",
    "maximum_automation", "use_rules_only",
    "abort_on_red_flags", "abort_on_warnings",
    "thinking_level",
    "expert_assessment", "strategy_memory",
    "log_text", "log_analysis",
    "intent", "command", "validation_error",
    "attempt_number", "previous_attempts",
    "metrics_history", "metrics_trend", "workflow_state",
    "stop", "stop_reason", "fallback_used",
    "red_flag_issues", "abort_message",
    "debug_log",
  ]
  for field in required:
    assert_in(field, state,
      "Missing field '%s' in state" % field)
  print("  PASSED")


def test_A03_default_values_correct():
  """Default values match expected baseline."""
  print("Test: A03_default_values_correct")
  state = create_initial_state(available_files=["x.mtz"])
  assert_equal(state["history"], [])
  assert_equal(state["available_files"], ["x.mtz"])
  assert_equal(state["cycle_number"], 1)
  assert_equal(state["max_cycles"], 20)
  assert_equal(state["user_advice"], "")
  assert_equal(state["provider"], "google")
  assert_none(state["session_resolution"])
  assert_equal(state["session_info"], {})
  assert_equal(state["directives"], {})
  assert_equal(state["bad_inject_params"], {})
  assert_true(state["maximum_automation"])
  assert_false(state["use_rules_only"])
  assert_none(state["thinking_level"])
  assert_equal(state["expert_assessment"], {})
  assert_equal(state["strategy_memory"], {})
  assert_true(state["abort_on_red_flags"])
  assert_false(state["abort_on_warnings"])
  assert_equal(state["log_text"], "")
  assert_equal(state["log_analysis"], {})
  assert_equal(state["intent"], {})
  assert_equal(state["command"], "")
  assert_none(state["validation_error"])
  assert_equal(state["attempt_number"], 0)
  assert_equal(state["previous_attempts"], [])
  assert_false(state["stop"])
  assert_none(state["stop_reason"])
  assert_false(state["fallback_used"])
  assert_equal(state["debug_log"], [])
  print("  PASSED")


def test_A04_all_kwargs_accepted():
  """create_initial_state accepts all documented kwargs."""
  print("Test: A04_all_kwargs_accepted")
  state = create_initial_state(
    available_files=["data.mtz"],
    log_text="some log",
    history=[{"program": "test"}],
    user_advice="advice",
    max_cycles=10,
    cycle_number=3,
    provider="openai",
    maximum_automation=False,
    session_resolution=2.5,
    use_rules_only=True,
    session_info={"experiment_type": "xray"},
    abort_on_red_flags=False,
    abort_on_warnings=True,
    directives={"run_once": ["phenix.xtriage"]},
    bad_inject_params={"phenix.refine": ["strategy"]})
  assert_equal(state["log_text"], "some log")
  assert_equal(state["max_cycles"], 10)
  assert_equal(state["cycle_number"], 3)
  assert_equal(state["provider"], "openai")
  assert_false(state["maximum_automation"])
  assert_equal(state["session_resolution"], 2.5)
  assert_true(state["use_rules_only"])
  assert_false(state["abort_on_red_flags"])
  assert_true(state["abort_on_warnings"])
  print("  PASSED")


def test_A05_create_initial_state_signature_stable():
  """create_initial_state parameter list has not changed."""
  print("Test: A05_create_initial_state_signature_stable")
  sig = inspect.signature(create_initial_state)
  param_names = list(sig.parameters.keys())
  expected = [
    "available_files", "log_text", "history", "user_advice",
    "max_cycles", "cycle_number", "provider",
    "maximum_automation", "session_resolution", "use_rules_only",
    "session_info", "abort_on_red_flags", "abort_on_warnings",
    "directives", "bad_inject_params",
    "thinking_level", "strategy_memory",
  ]
  for name in expected:
    assert_in(name, param_names,
      "Expected parameter '%s' missing from create_initial_state" % name)
  print("  PASSED")


# =========================================================================
# B. Routing Function Contracts
# =========================================================================

def test_B01_route_perceive_normal():
  """route_after_perceive returns 'think' for normal state."""
  print("Test: B01_route_perceive_normal")
  assert_equal(route_after_perceive({"stop": False}), "think")
  print("  PASSED")


def test_B02_route_perceive_stop():
  """route_after_perceive returns 'output' when stop=True."""
  print("Test: B02_route_perceive_stop")
  assert_equal(route_after_perceive({"stop": True}), "output")
  print("  PASSED")


def test_B03_route_perceive_missing_stop():
  """route_after_perceive defaults to 'think' if stop absent."""
  print("Test: B03_route_perceive_missing_stop")
  assert_equal(route_after_perceive({}), "think")
  print("  PASSED")


def test_B04_route_perceive_only_valid_outputs():
  """route_after_perceive only returns 'think' or 'output'."""
  print("Test: B04_route_perceive_only_valid_outputs")
  for stop_val in [True, False, None, 0, 1, ""]:
    result = route_after_perceive({"stop": stop_val})
    assert_in(result, ["think", "output"],
      "Unexpected route %r for stop=%r" % (result, stop_val))
  print("  PASSED")


def test_B05_route_validate_pass():
  """route_after_validate returns 'output' when valid."""
  print("Test: B05_route_validate_pass")
  state = {"stop": False, "validation_error": None}
  assert_equal(route_after_validate(state), "output")
  print("  PASSED")


def test_B06_route_validate_retry():
  """route_after_validate returns 'plan' for retry < 3."""
  print("Test: B06_route_validate_retry")
  state = {"stop": False, "validation_error": "bad",
           "attempt_number": 1}
  assert_equal(route_after_validate(state), "plan")
  print("  PASSED")


def test_B07_route_validate_fallback():
  """route_after_validate returns 'fallback' at attempt >= 3."""
  print("Test: B07_route_validate_fallback")
  state = {"stop": False, "validation_error": "bad",
           "attempt_number": 3}
  assert_equal(route_after_validate(state), "fallback")
  print("  PASSED")


def test_B08_route_validate_stop_overrides():
  """route_after_validate returns 'output' if stop=True."""
  print("Test: B08_route_validate_stop_overrides")
  state = {"stop": True, "validation_error": "bad",
           "attempt_number": 0}
  assert_equal(route_after_validate(state), "output")
  print("  PASSED")


def test_B09_route_validate_only_valid_outputs():
  """route_after_validate only returns valid routes."""
  print("Test: B09_route_validate_only_valid_outputs")
  valid = ["plan", "fallback", "output"]
  cases = [
    {"stop": True},
    {"stop": False, "validation_error": None},
    {"stop": False, "validation_error": "e", "attempt_number": 0},
    {"stop": False, "validation_error": "e", "attempt_number": 2},
    {"stop": False, "validation_error": "e", "attempt_number": 3},
    {"stop": False, "validation_error": "e", "attempt_number": 99},
  ]
  for state in cases:
    result = route_after_validate(state)
    assert_in(result, valid,
      "Unexpected route %r for state %r" % (result, state))
  print("  PASSED")


# =========================================================================
# C. State Serialization and Extensibility
# =========================================================================

def test_C01_state_json_roundtrip():
  """State round-trips through JSON without loss."""
  print("Test: C01_state_json_roundtrip")
  state = create_initial_state(
    available_files=["data.mtz", "model.pdb"],
    user_advice="try MR first",
    session_info={"experiment_type": "xray"})
  text = json.dumps(state)
  restored = json.loads(text)
  assert_equal(restored["user_advice"], "try MR first")
  assert_equal(restored["available_files"],
               ["data.mtz", "model.pdb"])
  assert_equal(restored["session_info"]["experiment_type"],
               "xray")
  print("  PASSED")


def test_C02_state_deepcopy_safe():
  """State can be deep-copied without mutation."""
  print("Test: C02_state_deepcopy_safe")
  state = create_initial_state(
    available_files=["data.mtz"],
    history=[{"program": "phenix.xtriage", "result": "OK"}],
    session_info={"best_files": {"data_mtz": "data.mtz"}})
  state2 = copy.deepcopy(state)
  state2["user_advice"] = "modified"
  state2["history"].append({"program": "new"})
  assert_equal(state["user_advice"], "")
  assert_equal(len(state["history"]), 1)
  print("  PASSED")


def test_C03_session_info_passthrough():
  """session_info dict passes through unchanged."""
  print("Test: C03_session_info_passthrough")
  info = {
    "experiment_type": "xray",
    "best_files": {"data_mtz": "data.mtz"},
    "rfree_mtz": "data_rfree.mtz",
  }
  state = create_initial_state(
    available_files=["data.mtz"],
    session_info=info)
  assert_equal(state["session_info"]["experiment_type"], "xray")
  assert_equal(state["session_info"]["rfree_mtz"], "data_rfree.mtz")
  print("  PASSED")


def test_C04_extra_session_info_keys_survive():
  """session_info preserves unknown keys (future-proofing)."""
  print("Test: C04_extra_session_info_keys_survive")
  info = {
    "experiment_type": "xray",
    "strategy_memory": {"data_quality": "good"},
    "thinking_level": "advanced",
    "totally_new_key": [1, 2, 3],
  }
  state = create_initial_state(
    available_files=["data.mtz"],
    session_info=info)
  assert_in("strategy_memory", state["session_info"])
  assert_in("thinking_level", state["session_info"])
  assert_in("totally_new_key", state["session_info"])
  assert_equal(
    state["session_info"]["strategy_memory"]["data_quality"],
    "good")
  print("  PASSED")


def test_C05_nested_session_info_json_roundtrip():
  """Nested session_info survives JSON round-trip."""
  print("Test: C05_nested_session_info_json_roundtrip")
  info = {
    "experiment_type": "xray",
    "strategy_memory": {
      "phasing_strategy": "MR",
      "r_free_history": [0.45, 0.42, 0.38],
      "concerns": ["borderline twinning"],
    },
  }
  state = create_initial_state(
    available_files=["data.mtz"],
    session_info=info)
  text = json.dumps(state)
  restored = json.loads(text)
  sm = restored["session_info"]["strategy_memory"]
  assert_equal(sm["phasing_strategy"], "MR")
  assert_equal(sm["r_free_history"], [0.45, 0.42, 0.38])
  print("  PASSED")


def test_C06_extra_state_keys_preserved():
  """Adding extra keys to state dict doesn't break it."""
  print("Test: C06_extra_state_keys_preserved")
  state = create_initial_state(available_files=["data.mtz"])
  state["thinking_level"] = "advanced"
  state["expert_assessment"] = {"analysis": "test"}
  state["strategy_memory"] = {"quality": "good"}
  # These survive json round-trip
  text = json.dumps(state)
  restored = json.loads(text)
  assert_equal(restored["thinking_level"], "advanced")
  assert_equal(
    restored["expert_assessment"]["analysis"], "test")
  assert_equal(
    restored["strategy_memory"]["quality"], "good")
  print("  PASSED")


# =========================================================================
# D. user_advice Enrichment Contracts
# =========================================================================

def test_D01_prepend_preserves_original():
  """Prepending expert text keeps original advice."""
  print("Test: D01_prepend_preserves_original")
  state = create_initial_state(
    available_files=["data.mtz"],
    user_advice="user says try MR")
  expert = "[Expert assessment] Data looks good for MR."
  state["user_advice"] = expert + "\n\n" + state["user_advice"]
  assert_in("user says try MR", state["user_advice"])
  assert_in("[Expert assessment]", state["user_advice"])
  print("  PASSED")


def test_D02_empty_enrichment_harmless():
  """Empty prepend doesn't break anything."""
  print("Test: D02_empty_enrichment_harmless")
  state = create_initial_state(
    available_files=["data.mtz"],
    user_advice="original")
  state["user_advice"] = "" + state["user_advice"]
  assert_equal(state["user_advice"], "original")
  print("  PASSED")


def test_D03_long_enrichment_not_truncated():
  """Long expert text is preserved in state."""
  print("Test: D03_long_enrichment_not_truncated")
  state = create_initial_state(available_files=["data.mtz"])
  long_text = "[Expert] " + "x" * 2000
  state["user_advice"] = long_text
  assert_equal(len(state["user_advice"]), 2009)
  print("  PASSED")


def test_D04_multiline_enrichment():
  """Multi-line expert text works."""
  print("Test: D04_multiline_enrichment")
  state = create_initial_state(
    available_files=["data.mtz"],
    user_advice="user advice")
  expert = (
    "[Expert assessment]\n"
    "Resolution 2.2A, completeness 99%.\n"
    "Weak anomalous signal.\n"
    "Try MR first.")
  state["user_advice"] = expert + "\n\n" + state["user_advice"]
  assert_in("Resolution 2.2A", state["user_advice"])
  assert_in("user advice", state["user_advice"])
  # JSON round-trip preserves newlines
  text = json.dumps(state)
  restored = json.loads(text)
  assert_in("Resolution 2.2A", restored["user_advice"])
  print("  PASSED")


# =========================================================================
# E. History / Backpack Extensibility
# =========================================================================

def test_E01_history_with_extra_fields():
  """History records with extra fields are not stripped."""
  print("Test: E01_history_with_extra_fields")
  history = [
    {"program": "phenix.xtriage", "result": "SUCCESS",
     "analysis": {"resolution": 2.5},
     "expert_assessment": {"analysis": "good data"},
     "strategy_memory": {"quality": "good"}},
  ]
  state = create_initial_state(
    available_files=["data.mtz"],
    history=history)
  rec = state["history"][0]
  assert_in("expert_assessment", rec)
  assert_equal(rec["expert_assessment"]["analysis"],
               "good data")
  print("  PASSED")


def test_E02_history_json_roundtrip_with_extras():
  """History with expert fields survives JSON."""
  print("Test: E02_history_json_roundtrip_with_extras")
  history = [
    {"program": "phenix.refine", "result": "SUCCESS",
     "analysis": {"r_free": 0.30},
     "expert_assessment": {
       "action": "continue",
       "analysis": "R-free improving steadily",
       "guidance": "Continue refinement",
     }},
  ]
  state = create_initial_state(
    available_files=["data.mtz"],
    history=history)
  text = json.dumps(state)
  restored = json.loads(text)
  rec = restored["history"][0]
  assert_equal(rec["expert_assessment"]["action"], "continue")
  print("  PASSED")


def test_E03_history_is_copied():
  """State copies history list (not by reference). Fixed in v113."""
  print("Test: E03_history_is_copied")
  history = [{"program": "x", "result": "OK"}]
  state = create_initial_state(
    available_files=["data.mtz"],
    history=history)
  # After v113 fix: state gets a copy, not the original
  assert_false(state["history"] is history,
    "State should copy history, not share reference")
  # But content should be equal
  assert_equal(state["history"], history)
  # Mutating state history should not affect original
  state["history"].append({"program": "y"})
  assert_equal(len(history), 1,
    "Original history must not be mutated")
  print("  PASSED")


# =========================================================================
# F. No-op Pass-Through Contract
# =========================================================================

def test_F01_noop_preserves_all_fields():
  """A no-op function preserves every state field."""
  print("Test: F01_noop_preserves_all_fields")
  def noop(state):
    return state

  state = create_initial_state(
    available_files=["data.mtz", "model.pdb"],
    user_advice="advice text",
    session_info={"experiment_type": "xray"})
  before = json.dumps(state, sort_keys=True)
  result = noop(state)
  after = json.dumps(result, sort_keys=True)
  assert_equal(before, after)
  print("  PASSED")


def test_F02_conditional_noop_pattern():
  """The think stub pattern works as a no-op."""
  print("Test: F02_conditional_noop_pattern")
  def think_stub(state):
    if not state.get("thinking_level"):
      return state
    # Would do thinking here
    return state

  state = create_initial_state(available_files=["data.mtz"])
  before = json.dumps(state, sort_keys=True)
  result = think_stub(state)
  after = json.dumps(result, sort_keys=True)
  assert_equal(before, after)
  print("  PASSED")


def test_F03_enriching_noop_pattern():
  """The enrichment pattern modifies only expected fields."""
  print("Test: F03_enriching_noop_pattern")
  def enrich(state):
    state = dict(state)  # shallow copy
    advice = state.get("user_advice", "")
    state["user_advice"] = "[Expert] test\n\n" + advice
    state["expert_assessment"] = {"analysis": "test"}
    return state

  state = create_initial_state(
    available_files=["data.mtz"],
    user_advice="original")
  result = enrich(state)
  # Modified fields
  assert_in("[Expert] test", result["user_advice"])
  assert_in("original", result["user_advice"])
  assert_equal(result["expert_assessment"]["analysis"], "test")
  # Unmodified fields
  assert_equal(result["available_files"], ["data.mtz"])
  assert_false(result["stop"])
  print("  PASSED")


def test_F04_noop_chain_equivalent():
  """Chaining no-ops doesn't change state."""
  print("Test: F04_noop_chain_equivalent")
  def noop(state):
    if not state.get("thinking_level"):
      return state
    return state

  state = create_initial_state(
    available_files=["data.mtz"],
    user_advice="test")
  before = json.dumps(state, sort_keys=True)
  result = noop(noop(noop(state)))
  after = json.dumps(result, sort_keys=True)
  assert_equal(before, after)
  print("  PASSED")


# =========================================================================
# G. graph.py Source Verification
#    (Catch if routing logic changes out from under us)
# =========================================================================

def test_G01_graph_py_has_route_after_perceive():
  """graph.py defines route_after_perceive."""
  print("Test: G01_graph_py_has_route_after_perceive")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  assert_in("def route_after_perceive", source)
  print("  PASSED")


def test_G02_graph_py_perceive_routes_to_think():
  """graph.py: perceive routes to think (added in Phase A3)."""
  print("Test: G02_graph_py_perceive_routes_to_think")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  assert_in('return "think"', source,
    "route_after_perceive should return 'think'")
  print("  PASSED")


def test_G03_graph_py_has_all_nodes():
  """graph.py registers all expected nodes."""
  print("Test: G03_graph_py_has_all_nodes")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  for node in ["perceive", "think", "plan", "build", "validate",
                "fallback", "output"]:
    assert_in('add_node("%s"' % node, source,
      "Missing node '%s' in graph.py" % node)
  print("  PASSED")


def test_G04_graph_py_has_think_node():
  """graph.py has a think node (added in Phase A3)."""
  print("Test: G04_graph_py_has_think_node")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  assert_in('add_node("think"', source,
    "think node should exist after Phase A3")
  assert_in('add_edge("think", "plan")', source,
    "think→plan edge should exist after Phase A3")
  print("  PASSED")


def test_G05_graph_py_edge_topology():
  """graph.py has expected edge structure."""
  print("Test: G05_graph_py_edge_topology")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  # Key edges that must exist
  assert_in('add_edge("think", "plan")', source)
  assert_in('add_edge("plan", "build")', source)
  assert_in('add_edge("build", "validate")', source)
  assert_in('add_edge("fallback", "output")', source)
  print("  PASSED")


# =========================================================================
# H. graph_state.py Source Verification
# =========================================================================

def test_H01_graph_state_has_thinking_fields():
  """graph_state.py has thinking fields (added in Phase A1)."""
  print("Test: H01_graph_state_has_thinking_fields")
  state_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_state.py")
  with open(state_path) as f:
    source = f.read()
  assert_in("thinking_level", source)
  assert_in("expert_assessment", source)
  assert_in("strategy_memory", source)
  print("  PASSED")


def test_H02_graph_state_field_count():
  """AgentState has exactly the expected number of fields."""
  print("Test: H02_graph_state_field_count")
  # Count fields in the TypedDict
  annotations = getattr(AgentState, '__annotations__', {})
  expected = 45
  actual = len(annotations)
  assert_equal(actual, expected,
    "AgentState has %d fields, expected %d. "
    "If you added thinking fields, update this test." %
    (actual, expected))
  print("  PASSED")


def test_H03_create_initial_state_param_count():
  """create_initial_state has expected parameter count."""
  print("Test: H03_create_initial_state_param_count")
  sig = inspect.signature(create_initial_state)
  # Current: 20 parameters (15 original + thinking_level,
  #   strategy_memory, structure_model, validation_history,
  #   session_blocked_programs)
  expected = 20
  actual = len(sig.parameters)
  assert_equal(actual, expected,
    "create_initial_state has %d params, expected %d. "
    "If you added thinking params, update this test." %
    (actual, expected))
  print("  PASSED")


# =========================================================================
# I. graph_nodes.py Source Verification
# =========================================================================

def test_I01_graph_nodes_has_think_function():
  """graph_nodes.py has think() stub (added in Phase A2)."""
  print("Test: I01_graph_nodes_has_think_function")
  nodes_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_nodes.py")
  with open(nodes_path, errors='ignore') as f:
    source = f.read()
  assert_in("def think(state)", source,
    "think() should exist after Phase A2")
  # Verify it checks thinking_level
  assert_in("thinking_level", source)
  print("  PASSED")


def test_I02_graph_nodes_plan_has_rules_only():
  """plan() function has rules-only early return."""
  print("Test: I02_graph_nodes_plan_has_rules_only")
  nodes_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_nodes.py")
  with open(nodes_path) as f:
    source = f.read()
  assert_in("use_rules_only", source)
  assert_in("_mock_plan", source)
  print("  PASSED")


def test_I03_graph_nodes_plan_has_directive_validation():
  """plan() function has directive validation section."""
  print("Test: I03_graph_nodes_plan_has_directive_validation")
  nodes_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_nodes.py")
  with open(nodes_path) as f:
    source = f.read()
  assert_in("VALIDATE INTENT AGAINST USER DIRECTIVES", source)
  assert_in("validate_intent", source)
  print("  PASSED")


def test_I04_graph_nodes_plan_has_forced_program():
  """plan() has forced program enforcement."""
  print("Test: I04_graph_nodes_plan_has_forced_program")
  nodes_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_nodes.py")
  with open(nodes_path) as f:
    source = f.read()
  assert_in("ENFORCE FORCED PROGRAM", source)
  assert_in("forced_program", source)
  print("  PASSED")


def test_I05_graph_nodes_exports_expected_functions():
  """graph.py imports expected functions from graph_nodes."""
  print("Test: I05_graph_nodes_exports_expected_functions")
  graph_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph.py")
  with open(graph_path) as f:
    source = f.read()
  for fn in ["perceive", "think", "plan", "build", "validate",
             "fallback", "output_node"]:
    assert_in(fn, source,
      "graph.py should import '%s' from graph_nodes" % fn)
  print("  PASSED")


# =========================================================================
# J. ai_agent.py Source Verification
# =========================================================================

def test_J01_ai_agent_has_agent_cycle_callback():
  """ai_agent.py sends agent_cycle callbacks."""
  print("Test: J01_ai_agent_has_agent_cycle_callback")
  agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(agent_path, errors='ignore') as f:
    source = f.read()
  assert_in("agent_cycle", source)
  assert_in('phase="decision"', source)
  assert_in('phase="result"', source)
  print("  PASSED")


def test_J02_ai_agent_has_thinking_plumbing():
  """ai_agent.py has thinking agent plumbing (added in Phase A5)."""
  print("Test: J02_ai_agent_has_thinking_plumbing")
  agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(agent_path, errors='ignore') as f:
    source = f.read()
  assert_in("strategy_memory", source,
    "strategy_memory should be in ai_agent.py after Phase A5")
  assert_in("expert_assessment", source,
    "expert_assessment should be in ai_agent.py after Phase A5")
  print("  PASSED")


# =========================================================================
# K. GUI Source Verification
# =========================================================================

def test_K01_gui_has_thinking_checkbox():
  """AIAgent.py GUI has thinking_level control (Phase A6)."""
  print("Test: K01_gui_has_thinking_checkbox")
  gui_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "wxGUI2", "Programs", "AIAgent.py")
  if not os.path.isfile(gui_path):
    print("  SKIP (AIAgent.py not found)")
    return
  with open(gui_path, errors='ignore') as f:
    source = f.read()
  assert_in("thinking_level", source,
    "thinking_level should be in GUI settings")
  print("  PASSED")


def test_K02_gui_has_agent_cycle_handler():
  """AIAgent.py has _on_agent_cycle callback handler."""
  print("Test: K02_gui_has_agent_cycle_handler")
  gui_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "wxGUI2", "Programs", "AIAgent.py")
  if not os.path.isfile(gui_path):
    print("  SKIP (AIAgent.py not found)")
    return
  with open(gui_path, errors='ignore') as f:
    source = f.read()
  assert_in("_on_agent_cycle", source)
  assert_in("agent_cycle", source)
  print("  PASSED")


def test_K03_gui_has_expert_display():
  """AIAgent.py displays expert assessments (added in Phase A6)."""
  print("Test: K03_gui_has_expert_display")
  gui_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "wxGUI2", "Programs", "AIAgent.py")
  if not os.path.isfile(gui_path):
    print("  SKIP (AIAgent.py not found)")
    return
  with open(gui_path, errors='ignore') as f:
    source = f.read()
  assert_in("expert_assessment", source,
    "expert_assessment should be in GUI after Phase A6")
  assert_in("Expert Review", source,
    "Expert Review display should be in GUI")
  print("  PASSED")


# =========================================================================
# L. Phase B File Absence Tests
#    (These should PASS now and FAIL once Phase B creates the files)
# =========================================================================

def test_L01_thinking_agent_py_exists():
  """agent/thinking_agent.py exists (created in Phase B4)."""
  print("Test: L01_thinking_agent_py_exists")
  path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "thinking_agent.py")
  assert_true(os.path.exists(path),
    "thinking_agent.py should exist after Phase B4")
  from agent.thinking_agent import run_think_node
  # Disabled state passes through
  result = run_think_node({"thinking_level": None})
  assert_false(result.get("stop"))
  print("  PASSED")


def test_L02_strategy_memory_py_exists():
  """agent/strategy_memory.py exists (created in Phase B1)."""
  print("Test: L02_strategy_memory_py_exists")
  path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "strategy_memory.py")
  assert_true(os.path.exists(path),
    "strategy_memory.py should exist after Phase B1")
  from agent.strategy_memory import StrategyMemory
  m = StrategyMemory()
  d = m.to_dict()
  m2 = StrategyMemory.from_dict(d)
  assert_equal(m2.to_dict(), d)
  print("  PASSED")


def test_L03_log_section_extractor_py_exists():
  """agent/log_section_extractor.py exists (created in Phase B2)."""
  print("Test: L03_log_section_extractor_py_exists")
  path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "log_section_extractor.py")
  assert_true(os.path.exists(path),
    "log_section_extractor.py should exist after Phase B2")
  from agent.log_section_extractor import extract_sections
  result = extract_sections("test line", "phenix.unknown")
  assert_true(len(result) > 0, "Should return fallback for unknown program")
  print("  PASSED")


def test_L04_thinking_prompts_py_exists():
  """knowledge/thinking_prompts.py exists (created in Phase B3)."""
  print("Test: L04_thinking_prompts_py_exists")
  path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "knowledge", "thinking_prompts.py")
  assert_true(os.path.exists(path),
    "thinking_prompts.py should exist after Phase B3")
  from knowledge.thinking_prompts import build_thinking_prompt
  sys_msg, user_msg = build_thinking_prompt({"cycle_number": 1})
  assert_true(len(sys_msg) > 100, "System prompt should be substantial")
  assert_true(len(user_msg) > 10, "User prompt should exist")
  print("  PASSED")


def test_L05_history_includes_file_provenance():
  """Think node history must show input/output files per cycle.

  Regression: Expert couldn't detect wrong-file-for-refinement
  because history only showed 'program (SUCCESS)' with no file info.
  """
  print("Test: L05_history_includes_file_provenance")
  from agent.thinking_agent import _summarize_history

  history = [
    {
      "cycle_number": 5,
      "program": "phenix.pdbtools",
      "command": "phenix.pdbtools refine_001.pdb ligand_fit_1.pdb",
      "result": "SUCCESS: OK",
      "output_files": ["/p/refine_001_modified.pdb"],
    },
    {
      "cycle_number": 6,
      "program": "phenix.refine",
      "command": "phenix.refine refine_001.pdb data.mtz",
      "result": "SUCCESS: OK",
      "output_files": ["/p/refine_002.pdb"],
    },
  ]

  summary = _summarize_history(history)

  # Must show input files
  assert_true("refine_001.pdb" in summary,
    "History should show input file basenames")
  assert_true("ligand_fit_1.pdb" in summary,
    "History should show ligand input file")

  # Must show output files
  assert_true("refine_001_modified.pdb" in summary,
    "History should show pdbtools output file")
  assert_true("refine_002.pdb" in summary,
    "History should show refine output file")

  # Must label Input/Output
  assert_true("Input:" in summary,
    "History should label input files")
  assert_true("Output:" in summary,
    "History should label output files")

  print("  PASSED")


def test_L06_system_prompt_mentions_file_provenance():
  """System prompt must instruct expert to check file provenance."""
  print("Test: L06_system_prompt_mentions_file_provenance")
  from knowledge.thinking_prompts import SYSTEM_PROMPT

  assert_true("provenance" in SYSTEM_PROMPT.lower()
              or "file" in SYSTEM_PROMPT.lower(),
    "System prompt should mention file checking")
  assert_true("RECENT HISTORY" in SYSTEM_PROMPT
              or "input" in SYSTEM_PROMPT.lower(),
    "System prompt should reference history section")

  print("  PASSED")


# =========================================================================
# Run
# =========================================================================

def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
