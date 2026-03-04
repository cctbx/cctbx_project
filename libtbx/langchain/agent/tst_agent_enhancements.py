"""
Integration test for the agent enhancement (Phases A-D).

Run from the PHENIX source tree:
  libtbx.python tst_agent_enhancements.py

Tests both with and without PHENIX imports.
All assertions pass in standalone mode (no PHENIX).
Real structural validation tests require PHENIX.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import traceback

# Ensure the project root is on sys.path so
# imports like 'from agent.X' work when running
# as 'python tst_agent_enhancements.py' from any dir.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)
if _PROJECT_ROOT not in sys.path:
  sys.path.insert(0, _PROJECT_ROOT)


def run_tests():
  """Run all integration tests."""
  passed = 0
  failed = 0
  skipped = 0

  def test(name, fn):
    nonlocal passed, failed, skipped
    try:
      result = fn()
      if result == "SKIP":
        print("  SKIP: %s" % name)
        skipped += 1
      else:
        print("  PASS: %s" % name)
        passed += 1
    except Exception as e:
      print("  FAIL: %s — %s" % (name, e))
      traceback.print_exc()
      failed += 1

  print("=" * 60)
  print("Agent Enhancement Integration Tests")
  print("=" * 60)
  print()

  # --- Phase A: Validation ---
  print("Phase A: Validation")
  test("import validation_inspector",
    test_import_validation_inspector)
  test("import format_validation",
    test_import_format_validation)
  test("run_validation graceful fallback",
    test_validation_fallback)
  test("format_validation_report full",
    test_format_report_full)
  test("format_validation_report None input",
    test_format_report_none)
  test("format_validation_trend",
    test_format_trend)
  test("real structural validation",
    test_real_validation)
  print()

  # --- Phase C: File Metadata ---
  print("Phase C: File Metadata")
  test("build_file_metadata",
    test_build_metadata)
  test("find_best_model",
    test_find_best)
  test("get_model_summary",
    test_model_summary)
  test("format_file_metadata_block",
    test_metadata_block)
  print()

  # --- Phase D: Knowledge Base ---
  print("Phase D: Knowledge Base")
  test("KB loader",
    test_kb_loader)
  test("derive_tags refinement",
    test_tags_refinement)
  test("derive_tags stuck",
    test_tags_stuck)
  test("derive_tags cryoem",
    test_tags_cryoem)
  test("KB query refinement",
    test_kb_query_refinement)
  test("KB query stuck (cross-category)",
    test_kb_query_stuck)
  print()

  # --- Integration ---
  print("Integration")
  test("should_think includes refine",
    test_should_think_refine)
  test("metrics access (flat dict)",
    test_metrics_flat)
  test("R-free trend from analysis key",
    test_rfree_trend_analysis_key)
  test("full prompt assembly",
    test_full_prompt)
  test("run_think_node no crash",
    test_run_think_node)
  print()

  # --- Thinking Level Parameter ---
  print("Thinking Level Parameter")
  test("create_initial_state thinking_level",
    test_create_initial_state_thinking_level)
  test("thinking_level basic and advanced",
    test_backward_compat_thinking_agent)
  test("think node routing",
    test_think_node_routing)
  test("basic vs advanced context",
    test_basic_vs_advanced_context)
  print()

  print("=" * 60)
  print(
    "Results: %d passed, %d failed, %d skipped"
    % (passed, failed, skipped)
  )
  print("=" * 60)
  return failed == 0


# =========================================================
# Phase A tests
# =========================================================

def test_import_validation_inspector():
  from agent.validation_inspector import run_validation
  assert callable(run_validation)

def test_import_format_validation():
  from agent.format_validation import (
    format_validation_report,
    format_validation_trend,
  )
  assert callable(format_validation_report)
  assert callable(format_validation_trend)

def test_validation_fallback():
  from agent.validation_inspector import run_validation
  assert run_validation(None) is None
  assert run_validation("") is None
  assert run_validation("/nonexistent.pdb") is None

def test_format_report_full():
  from agent.format_validation import (
    format_validation_report,
  )
  vr = {
    "geometry": {
      "rama_favored": 0.972,
      "rama_outliers": 0.0,
      "rotamer_outliers": 0.013,
      "rotamer_outlier_list": [
        "A/Leu47", "A/Val102",
      ],
      "clashscore": 3.2,
    },
    "model_contents": {
      "chains": ["A", "B"],
      "ligands": [{
        "name": "ATP", "chain": "A",
        "resid": 301, "n_atoms": 31,
      }],
      "waters": 187,
      "ions": [{"name": "MG"}],
    },
  }
  report = format_validation_report(
    vr,
    log_metrics={
      "r_free": 0.248, "r_work": 0.219,
      "bonds_rmsd": 0.006, "angles_rmsd": 0.82,
    },
    cycle_number=5,
    program_name="phenix.refine",
    prev_r_free=0.295,
    start_r_free=0.421,
  )
  assert "R-free=0.248" in report
  assert "Bonds=0.0060" in report
  assert "Rama: 97.2%" in report
  assert "ATP (A/301)" in report
  assert "Waters: 187" in report
  assert len(report) < 500

def test_format_report_none():
  from agent.format_validation import (
    format_validation_report,
  )
  report = format_validation_report(
    None,
    log_metrics={"r_free": 0.30},
    cycle_number=3,
  )
  assert "not available" in report
  assert "R-free=0.300" in report

def test_format_trend():
  from agent.format_validation import (
    format_validation_trend,
  )
  t = format_validation_trend([0.42, 0.30, 0.25])
  assert "0.420" in t
  assert "0.250" in t
  assert "3 cycles" in t
  assert format_validation_trend([0.3]) == ""
  assert format_validation_trend([]) == ""
  assert format_validation_trend(None) == ""

def test_real_validation():
  """Test with real PHENIX validation (requires mmtbx)."""
  try:
    import iotbx.data_manager
    assert iotbx.data_manager is not None
  except ImportError:
    return "SKIP"

  from agent.validation_inspector import run_validation

  # Use a tutorial model if available
  try:
    import libtbx.load_env
    test_dir = libtbx.env.under_dist(
      "phenix_regression", "wizards"
    )
    model = os.path.join(
      test_dir, "p9_refine_001.pdb"
    )
    if not os.path.isfile(model):
      return "SKIP"
  except Exception:
    return "SKIP"

  result = run_validation(model)
  assert result is not None, \
    "Real validation should succeed"
  assert "geometry" in result
  assert "model_contents" in result
  geom = result["geometry"]
  assert geom["rama_favored"] > 0
  assert geom["clashscore"] >= 0
  mc = result["model_contents"]
  assert len(mc["chains"]) > 0


# =========================================================
# Phase C tests
# =========================================================

def test_build_metadata():
  from agent.file_metadata import build_file_metadata
  meta = build_file_metadata(
    "/work/model_005.pdb",
    validation_result={
      "model_contents": {
        "chains": ["A"],
        "ligands": [{"name": "ATP", "chain": "A",
                     "resid": 301, "n_atoms": 31}],
        "waters": 100, "ions": [],
      },
      "geometry": {"clashscore": 3.0},
    },
    log_metrics={"r_free": 0.25},
    program_name="phenix.refine",
    cycle_number=5,
  )
  assert meta["quality"]["r_free"] == 0.25
  assert meta["contents"]["chains"] == ["A"]
  assert meta["provenance"]["cycle"] == 5

def test_find_best():
  from agent.file_metadata import (
    build_file_metadata, find_best_model,
  )
  m1 = build_file_metadata(
    "/w/a.pdb", log_metrics={"r_free": 0.40},
    cycle_number=1,
  )
  m2 = build_file_metadata(
    "/w/b.pdb", log_metrics={"r_free": 0.25},
    cycle_number=2,
  )
  fm = {m1["path"]: m1, m2["path"]: m2}
  best = find_best_model(fm)
  assert "b.pdb" in best

def test_model_summary():
  from agent.file_metadata import (
    build_file_metadata, get_model_summary,
  )
  m = build_file_metadata(
    "/w/m.pdb",
    validation_result={
      "model_contents": {
        "chains": ["A", "B"],
        "ligands": [{"name": "ATP", "chain": "A",
                     "resid": 1, "n_atoms": 10}],
        "waters": 50, "ions": [],
      },
      "geometry": {},
    },
    log_metrics={"r_free": 0.25},
  )
  s = get_model_summary({m["path"]: m}, m["path"])
  assert "A,B" in s
  assert "ATP" in s
  assert "R-free=0.250" in s

def test_metadata_block():
  from agent.file_metadata import (
    build_file_metadata,
    format_file_metadata_block,
  )
  m = build_file_metadata(
    "/w/m.pdb",
    log_metrics={"r_free": 0.25},
    cycle_number=1,
  )
  block = format_file_metadata_block(
    {m["path"]: m}
  )
  assert "FILE METADATA" in block
  assert "Best model" in block
  assert format_file_metadata_block({}) == ""


# =========================================================
# Phase D tests
# =========================================================

def test_kb_loader():
  from knowledge.kb_loader import ExpertKnowledgeBase
  # Try multiple locations for v2 and v1 YAML
  _kb_dir = os.path.join(
    os.path.dirname(__file__), "..", "knowledge",
  )
  candidates = [
    os.path.join(_kb_dir,
      "expert_knowledge_base_v2.yaml"),
    os.path.join(_kb_dir,
      "expert_knowledge_base.yaml"),
    os.path.join(os.path.dirname(__file__),
      "expert_knowledge_base_v2.yaml"),
    os.path.join(os.path.dirname(__file__),
      "expert_knowledge_base.yaml"),
  ]
  kb_path = None
  for c in candidates:
    if os.path.isfile(c):
      kb_path = c
      break
  if kb_path is None:
    return "SKIP"
  kb = ExpertKnowledgeBase(kb_path)
  # v2 has ~56 entries; v1 had ~290.
  # Either version should load successfully.
  assert kb.stats["total"] >= 10, (
    "KB too small: %d entries" % kb.stats["total"]
  )
  # v2 must have the core categories
  if "v2" in kb_path:
    cats = kb.stats["by_category"]
    assert "stopping" in cats, (
      "Missing 'stopping' category in v2 KB"
    )
    assert "phenix_pitfalls" in cats, (
      "Missing 'phenix_pitfalls' in v2 KB"
    )

def test_tags_refinement():
  from agent.kb_tags import derive_tags
  cat, tags = derive_tags(
    resolution=2.5,
    log_metrics={"r_free": 0.42, "r_work": 0.35},
    workflow_stage="refinement",
  )
  assert cat == "refinement"
  assert "r_free_high" in tags
  assert "medium_resolution" in tags

def test_tags_stuck():
  from agent.kb_tags import derive_tags
  _, tags = derive_tags(
    r_free_trend=[0.38, 0.38, 0.38],
    log_metrics={"r_free": 0.38},
    workflow_stage="refinement",
  )
  assert "r_free_stuck" in tags
  assert "plateau" in tags

def test_tags_cryoem():
  from agent.kb_tags import derive_tags
  cat, tags = derive_tags(
    experiment_type="cryoem",
    workflow_stage="cryoem_refinement",
  )
  assert cat == "cryoem"
  assert "cryoem" in tags

def test_kb_query_refinement():
  from agent.thinking_agent import (
    _query_knowledge_base,
  )
  rules = _query_knowledge_base(
    log_metrics={"r_free": 0.30},
    workflow_stage="refinement",
    resolution=2.5,
  )
  assert "EXPERT RULES" in rules

def test_kb_query_stuck():
  from agent.thinking_agent import (
    _query_knowledge_base,
  )
  rules = _query_knowledge_base(
    log_metrics={"r_free": 0.38},
    workflow_stage="refinement",
    r_free_trend=[0.38, 0.38, 0.38],
  )
  assert "EXPERT RULES" in rules
  # Should get failure/stopping entries
  assert "fail_" in rules or "stop_" in rules


# =========================================================
# Integration tests
# =========================================================

def test_should_think_refine():
  from agent.thinking_agent import should_think
  state = {
    "history": [
      {"program": "phenix.refine", "status": "OK"},
    ],
    "strategy_memory": {},
  }
  assert should_think(state) is True
  # First cycle still skips
  assert should_think({
    "history": [], "strategy_memory": {}
  }) is False
  # Program from log_analysis takes priority over
  # history[-1] (which may be the previous cycle)
  state2 = {
    "history": [
      {"program": "phenix.xtriage", "status": "OK"},
    ],
    "log_analysis": {"program": "phenix.refine"},
    "strategy_memory": {},
  }
  assert should_think(state2) is True

def test_metrics_flat():
  """log_analysis is a flat dict, not nested."""
  from knowledge.thinking_prompts import (
    _format_metrics,
  )
  flat = {
    "program": "phenix.refine",
    "r_free": 0.248, "r_work": 0.219,
    "bonds_rmsd": 0.006, "resolution": 2.5,
    "macro_cycles": 3,
  }
  fmt = _format_metrics(flat)
  assert "r_free=0.248" in fmt
  assert "bonds_rmsd=0.006" in fmt
  assert "program=" not in fmt
  assert "macro_cycles=" not in fmt

def test_rfree_trend_analysis_key():
  """History entries use 'analysis' key."""
  from agent.thinking_agent import (
    _collect_r_free_trend,
    _get_prev_r_free,
    _get_start_r_free,
  )
  history = [
    {"analysis": {"program": "xtriage"}},
    {"analysis": {"r_free": 0.42}},
    {"analysis": {"r_free": 0.30}},
    {"analysis": {"r_free": 0.25}},
  ]
  trend = _collect_r_free_trend(history)
  assert trend == [0.42, 0.30, 0.25]
  # history[-1] is the previous cycle (current cycle
  # is NOT in history — it's still being processed)
  assert _get_prev_r_free(history) == 0.25
  # Start R-free skips xtriage (no R-free) and finds
  # the first cycle that actually had one
  assert _get_start_r_free(history) == 0.42
  # Also works with empty/missing
  assert _collect_r_free_trend([]) == []
  assert _get_start_r_free([]) is None
  assert _get_prev_r_free([]) is None
  # Scans backwards past non-refine cycles
  gap_history = [
    {"analysis": {"r_free": 0.35}},
    {"analysis": {"r_free": 0.30}},
    {"analysis": {"program": "ligandfit"}},
  ]
  assert _get_prev_r_free(gap_history) == 0.30

def test_full_prompt():
  from agent.thinking_agent import (
    _build_thinking_context,
  )
  from knowledge.thinking_prompts import (
    build_thinking_prompt,
  )
  state = {
    "history": [
      {"cycle_number": 1,
       "program": "phenix.xtriage",
       "analysis": {"program": "xtriage"}},
      {"cycle_number": 2,
       "program": "phenix.phaser",
       "analysis": {"TFZ": 12.3}},
      {"cycle_number": 3,
       "program": "phenix.refine",
       "analysis": {"r_free": 0.35}},
      {"cycle_number": 4,
       "program": "phenix.refine",
       "analysis": {"r_free": 0.28}},
    ],
    "log_text": "Refinement done.",
    "session_info": {"experiment_type": "xray"},
    "workflow_state": {"state": "refinement"},
    "log_analysis": {
      "program": "phenix.refine",
      "r_free": 0.248, "r_work": 0.219,
      "resolution": 2.5, "bonds_rmsd": 0.006,
    },
    "cycle_number": 5,
    "user_advice": "",
  }
  # Advanced mode: includes KB rules
  ctx = _build_thinking_context(state, "advanced")
  _, msg = build_thinking_prompt(ctx)
  assert "CYCLE 5" in msg
  assert "r_free=0.248" in msg
  assert "bonds_rmsd=0.006" in msg
  # Trend includes 2 from history + current appended
  assert "0.350 -> 0.280 -> 0.248" in msg
  assert "EXPERT RULES" in msg
  assert "program=" not in msg
  # Program comes from log_analysis, not history
  assert ctx["program_name"] == "phenix.refine"
  assert ctx["r_free_trend"] == [0.35, 0.28, 0.248]

  # Basic mode: no KB rules, no validation
  ctx_basic = _build_thinking_context(state, "basic")
  assert ctx_basic["kb_rules_text"] == ""
  assert ctx_basic["validation_result"] is None
  _, msg_basic = build_thinking_prompt(ctx_basic)
  assert "CYCLE 5" in msg_basic
  assert "r_free=0.248" in msg_basic
  # Basic still has R-free trend
  assert "0.350 -> 0.280 -> 0.248" in msg_basic

def test_run_think_node():
  from agent.thinking_agent import run_think_node
  # Advanced level
  state = {
    "thinking_level": "advanced",
    "history": [
      {"program": "phenix.refine",
       "status": "OK", "cycle_number": 1,
       "analysis": {"r_free": 0.35}},
    ],
    "log_text": "",
    "session_info": {"experiment_type": "xray"},
    "workflow_state": {"state": "refinement"},
    "log_analysis": {"r_free": 0.30},
    "cycle_number": 2,
    "user_advice": "",
    "strategy_memory": {},
    "provider": "google",
    "debug_log": [],
    "events": [],
  }
  result = run_think_node(state)
  assert "debug_log" in result

  # None level → pass-through (unchanged)
  state_none = {**state, "thinking_level": None}
  result_none = run_think_node(state_none)
  assert result_none is state_none

  # Basic level
  state_basic = {**state, "thinking_level": "basic"}
  result_basic = run_think_node(state_basic)
  assert "debug_log" in result_basic


# =========================================================
# Thinking level parameter tests
# =========================================================

def test_create_initial_state_thinking_level():
  """create_initial_state handles thinking_level."""
  from agent.graph_state import create_initial_state
  # Default: no thinking
  s = create_initial_state(
    available_files=["a.pdb"]
  )
  assert s["thinking_level"] is None

  # Explicit advanced
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="advanced"
  )
  assert s["thinking_level"] == "advanced"

  # Explicit basic
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="basic"
  )
  assert s["thinking_level"] == "basic"

  # Invalid string → None
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="bogus"
  )
  assert s["thinking_level"] is None

  # Case insensitive
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="ADVANCED"
  )
  assert s["thinking_level"] == "advanced"

def test_backward_compat_thinking_agent():
  """thinking_level=basic explicitly set works correctly."""
  from agent.graph_state import create_initial_state
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="basic"
  )
  assert s["thinking_level"] == "basic"

  # Advanced takes priority over basic
  s = create_initial_state(
    available_files=["a.pdb"],
    thinking_level="advanced"
  )
  assert s["thinking_level"] == "advanced"

def test_think_node_routing():
  """think() node respects thinking_level."""
  # We can't import the real graph_nodes.think()
  # without LangGraph, but we test the core logic
  # via run_think_node which it calls.
  from agent.thinking_agent import run_think_node

  # None → pass-through
  s = {"thinking_level": None, "debug_log": []}
  assert run_think_node(s) is s

  # Empty string → pass-through
  s2 = {"thinking_level": "", "debug_log": []}
  assert run_think_node(s2) is s2

def test_basic_vs_advanced_context():
  """Basic mode omits validation and KB; advanced
  includes them."""
  from agent.thinking_agent import (
    _build_thinking_context,
  )
  state = {
    "history": [
      {"program": "phenix.refine",
       "cycle_number": 1,
       "analysis": {"r_free": 0.30}},
    ],
    "log_text": "Final R-free = 0.25",
    "session_info": {"experiment_type": "xray"},
    "workflow_state": {"state": "refinement"},
    "log_analysis": {
      "program": "phenix.refine",
      "r_free": 0.25, "resolution": 2.0,
    },
    "cycle_number": 2,
    "user_advice": "",
  }

  basic = _build_thinking_context(state, "basic")
  assert basic["kb_rules_text"] == ""
  assert basic["validation_result"] is None
  assert basic["validation_report"] == ""
  assert basic["file_metadata"] == {}
  # But basic still has metrics and trend
  assert basic["metrics"]["r_free"] == 0.25
  assert len(basic["r_free_trend"]) >= 1

  adv = _build_thinking_context(state, "advanced")
  # Advanced attempts KB query (returns text even
  # without validation result, based on tags)
  assert isinstance(adv["kb_rules_text"], str)
  # Advanced has same base metrics
  assert adv["metrics"]["r_free"] == 0.25


if __name__ == "__main__":
  success = run_tests()
  sys.exit(0 if success else 1)
