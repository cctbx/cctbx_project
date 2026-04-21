"""Tests for agent.display_data_model."""

from __future__ import (
  absolute_import, division, print_function,
)

import os
import sys

# Allow running from project root
_this_dir = os.path.dirname(os.path.abspath(__file__))
_project_dir = os.path.dirname(_this_dir)
if _project_dir not in sys.path:
  sys.path.insert(0, _project_dir)

from tests.tst_utils import (
  assert_equal, assert_true, assert_false,
  assert_none,
)
from agent.display_data_model import (
  DisplayDataModel,
  _extract_key_metric,
  _fmt_rfree, _fmt_pct, _fmt_float,
  _format_metric_value,
)


# ── Formatting helpers ─────────────────────────

def test_fmt_rfree():
  assert_equal(_fmt_rfree(0.231), "0.231")
  assert_equal(_fmt_rfree(0.3), "0.300")
  assert_equal(_fmt_rfree(None), None)
  assert_equal(_fmt_rfree("0.25"), "0.250")
  assert_equal(_fmt_rfree("bad"), None)
  print("  PASS: test_fmt_rfree")


def test_fmt_pct():
  assert_equal(_fmt_pct(0.978), "97.8%")
  assert_equal(_fmt_pct(97.8), "97.8%")
  assert_equal(_fmt_pct(0.01), "1.0%")
  assert_equal(_fmt_pct(1.0), "100.0%")
  assert_equal(_fmt_pct(None), None)
  print("  PASS: test_fmt_pct")


def test_fmt_float():
  assert_equal(_fmt_float(3.14159, 2), "3.14")
  assert_equal(_fmt_float(3.14159, 1), "3.1")
  assert_equal(_fmt_float(None), None)
  print("  PASS: test_fmt_float")


def test_format_metric_value():
  # R-free → 3 decimals
  assert_equal(_format_metric_value(
    "R-free", 0.231), "0.231")
  # FOM → 3 decimals
  assert_equal(_format_metric_value(
    "FOM", 0.42), "0.420")
  # Map CC → 3 decimals
  assert_equal(_format_metric_value(
    "Map CC", 0.82), "0.820")
  # TFZ → 1 decimal
  assert_equal(_format_metric_value(
    "TFZ", 14.2), "14.2")
  # Res → 2 decimals
  assert_equal(_format_metric_value(
    "Res", 2.5), "2.50")
  # Clash → 1 decimal
  assert_equal(_format_metric_value(
    "Clash", 4.2), "4.2")
  print("  PASS: test_format_metric_value")


# ── Key metric extraction ──────────────────────

def test_extract_key_metric_refine():
  name, val = _extract_key_metric(
    "phenix.refine", {"r_free": 0.25, "r_work": 0.21})
  assert_equal(name, "R-free")
  assert_equal(val, 0.25)
  print("  PASS: test_extract_key_metric_refine")


def test_extract_key_metric_rsr():
  """real_space_refine must NOT match phenix.refine."""
  name, val = _extract_key_metric(
    "phenix.real_space_refine",
    {"model_map_cc": 0.75})
  assert_equal(name, "Map CC")
  assert_equal(val, 0.75)
  print("  PASS: test_extract_key_metric_rsr")


def test_extract_key_metric_phaser():
  name, val = _extract_key_metric(
    "phenix.phaser", {"tfz": 14.2, "llg": 342.0})
  assert_equal(name, "TFZ")
  assert_equal(val, 14.2)
  print("  PASS: test_extract_key_metric_phaser")


def test_extract_key_metric_autosol():
  name, val = _extract_key_metric(
    "phenix.autosol", {"fom": 0.42, "bayes_cc": 54.2})
  assert_equal(name, "FOM")
  assert_equal(val, 0.42)
  print("  PASS: test_extract_key_metric_autosol")


def test_extract_key_metric_ligandfit():
  name, val = _extract_key_metric(
    "phenix.ligandfit", {"cc": 0.85})
  assert_equal(name, "Ligand CC")
  assert_equal(val, 0.85)
  print("  PASS: test_extract_key_metric_ligandfit")


def test_extract_key_metric_no_metrics():
  name, val = _extract_key_metric(
    "phenix.refine", {})
  assert_equal(name, "")
  assert_none(val)
  print("  PASS: test_extract_key_metric_no_metrics")


def test_extract_key_metric_fallback_rfree():
  """Unknown program should fall back to r_free."""
  name, val = _extract_key_metric(
    "phenix.unknown_program", {"r_free": 0.30})
  assert_equal(name, "R-free")
  assert_equal(val, 0.30)
  print("  PASS: test_extract_key_metric_fallback_rfree")


def test_extract_key_metric_fallback_cc():
  """Unknown program should fall back to model_map_cc
  when no r_free."""
  name, val = _extract_key_metric(
    "phenix.unknown_program",
    {"model_map_cc": 0.80})
  assert_equal(name, "Map CC")
  assert_equal(val, 0.80)
  print("  PASS: test_extract_key_metric_fallback_cc")


# ── from_session ───────────────────────────────

def test_from_session_none():
  ddm = DisplayDataModel.from_session(None)
  assert_equal(ddm.outcome_status, "incomplete")
  assert_equal(ddm.timeline, [])
  print("  PASS: test_from_session_none")


def test_from_session_empty():
  ddm = DisplayDataModel.from_session({})
  assert_equal(ddm.outcome_status, "incomplete")
  assert_equal(ddm.final_metrics, {})
  assert_equal(ddm.stage_outcomes, [])
  print("  PASS: test_from_session_empty")


# ── Outcome ────────────────────────────────────

def _make_mr_session():
  return {
    "experiment_type": "xray",
    "structure_model": {
      "_version": 1,
      "data_characteristics": {
        "resolution": 2.0,
        "space_group": "P212121",
        "mr_tfz": 14.2,
        "mr_llg": 342.0,
      },
      "model_state": {
        "r_free": 0.231,
        "r_work": 0.198,
        "chains": [{"chain_id": "A"}],
        "ligands": [],
        "waters": 88,
        "geometry": {
          "clashscore": 3.2,
          "rama_favored": 0.978,
        },
      },
      "progress": [],
      "strategy_blacklist": [],
    },
    "plan": {
      "goal": "MR + refinement",
      "phases": [
        {"id": "data_assessment",
         "status": "complete",
         "description": "Analyze data",
         "programs": ["phenix.xtriage"],
         "cycles_used": 1},
        {"id": "mr", "status": "complete",
         "description": "Find MR solution",
         "programs": ["phenix.phaser"],
         "cycles_used": 1},
        {"id": "refine", "status": "complete",
         "description": "Refine model",
         "programs": ["phenix.refine"],
         "cycles_used": 3},
      ],
      "_version": 1,
    },
    "best_files": {
      "model": "/path/to/refine_003.pdb",
      "map_coeffs_mtz": "/path/to/refine_003.mtz",
    },
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.xtriage",
       "result": "SUCCESS: OK",
       "metrics": {"resolution": 2.0}},
      {"cycle_number": 2,
       "program": "phenix.phaser",
       "result": "SUCCESS: OK",
       "metrics": {"tfz": 14.2, "llg": 342.0}},
      {"cycle_number": 3,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.350}},
      {"cycle_number": 4,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.285}},
      {"cycle_number": 5,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.231}},
    ],
  }


def test_outcome_determined():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  assert_equal(ddm.outcome_status, "determined")
  assert_true("0.231" in ddm.outcome_message)
  print("  PASS: test_outcome_determined")


def test_outcome_stopped():
  ddm = DisplayDataModel.from_session({
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.refine",
       "result": "FAILED: error"},
    ],
    "structure_model": {
      "_version": 1,
      "model_state": {"r_free": 0.48},
      "data_characteristics": {},
    },
  })
  assert_equal(ddm.outcome_status, "stopped")
  assert_true("0.480" in ddm.outcome_message)
  print("  PASS: test_outcome_stopped")


def test_outcome_cryoem():
  ddm = DisplayDataModel.from_session({
    "experiment_type": "cryoem",
    "structure_model": {
      "_version": 1,
      "model_state": {"model_map_cc": 0.82},
      "data_characteristics": {
        "experiment_type": "cryoem"},
    },
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.82}},
    ],
  })
  assert_equal(ddm.outcome_status, "determined")
  assert_true("0.820" in ddm.outcome_message)
  assert_true("Map CC" in ddm.outcome_message)
  print("  PASS: test_outcome_cryoem")


# ── Final metrics ──────────────────────────────

def test_final_metrics_full():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  m = ddm.final_metrics
  assert_equal(m["R-free"], "0.231")
  assert_equal(m["R-work"], "0.198")
  assert_equal(m["Clashscore"], "3.2")
  assert_equal(m["Ramachandran favored"], "97.8%")
  assert_equal(m["Chains"], "1")
  assert_equal(m["Waters"], "88")
  assert_equal(m["Resolution"], "2.00 Å")
  assert_equal(m["Space group"], "P212121")
  print("  PASS: test_final_metrics_full")


def test_final_metrics_fallback_no_sm():
  """Without structure_model, fall back to cycles."""
  ddm = DisplayDataModel.from_session({
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.25, "r_work": 0.21}},
    ],
  })
  m = ddm.final_metrics
  assert_equal(m["R-free"], "0.250")
  assert_equal(m["R-work"], "0.210")
  print("  PASS: test_final_metrics_fallback_no_sm")


# ── Timeline ───────────────────────────────────

def test_timeline_basic():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  tl = ddm.timeline
  assert_equal(len(tl), 5)
  assert_equal(tl[0].program, "phenix.xtriage")
  assert_equal(tl[0].status, "OK")
  assert_equal(tl[0].key_metric_name, "Res")
  assert_equal(tl[1].key_metric_name, "TFZ")
  assert_equal(tl[2].key_metric_name, "R-free")
  print("  PASS: test_timeline_basic")


def test_timeline_rfree_deltas():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  tl = ddm.timeline
  # Cycle 3: first R-free, no delta
  assert_none(tl[2].key_metric_delta)
  # Cycle 4: delta from 0.350 → 0.285
  assert_true(
    abs(tl[3].key_metric_delta - (-0.065)) < 0.001)
  # Cycle 5: delta from 0.285 → 0.231
  assert_true(
    abs(tl[4].key_metric_delta - (-0.054)) < 0.001)
  print("  PASS: test_timeline_rfree_deltas")


def test_timeline_mapcc_deltas():
  ddm = DisplayDataModel.from_session({
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.65}},
      {"cycle_number": 2,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.75}},
    ],
  })
  tl = ddm.timeline
  assert_none(tl[0].key_metric_delta)
  assert_true(
    abs(tl[1].key_metric_delta - 0.10) < 0.001)
  print("  PASS: test_timeline_mapcc_deltas")


# ── R-free trajectory ──────────────────────────

def test_rfree_trajectory():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  traj = ddm.rfree_trajectory
  assert_equal(len(traj), 3)  # cycles 3,4,5
  assert_equal(traj[0].cycle, 3)
  assert_true(abs(traj[0].value - 0.350) < 0.001)
  assert_equal(traj[2].cycle, 5)
  assert_true(abs(traj[2].value - 0.231) < 0.001)
  print("  PASS: test_rfree_trajectory")


def test_cc_trajectory():
  ddm = DisplayDataModel.from_session({
    "experiment_type": "cryoem",
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.65}},
      {"cycle_number": 2,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.82}},
    ],
  })
  traj = ddm.cc_trajectory
  assert_equal(len(traj), 2)
  assert_true(abs(traj[0].value - 0.65) < 0.001)
  print("  PASS: test_cc_trajectory")


def test_primary_trajectory_xray():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  traj = ddm.primary_trajectory
  assert_true(len(traj) > 0)
  # Should be R-free values
  assert_true(traj[0].value < 1.0)
  print("  PASS: test_primary_trajectory_xray")


def test_primary_trajectory_cryoem():
  ddm = DisplayDataModel.from_session({
    "experiment_type": "cryoem",
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.75}},
    ],
  })
  traj = ddm.primary_trajectory
  assert_equal(len(traj), 1)
  assert_true(abs(traj[0].value - 0.75) < 0.001)
  print("  PASS: test_primary_trajectory_cryoem")


def test_retreat_markers():
  ddm = DisplayDataModel.from_session({
    "structure_model": {
      "_version": 1,
      "model_state": {},
      "data_characteristics": {},
      "progress": [],
      "strategy_blacklist": [
        {"retreat_cycle": 4, "strategy_id": "s1",
         "reason": "stuck"},
      ],
    },
    "cycles": [
      {"cycle_number": 3,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.40}},
      {"cycle_number": 4,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.39}},
      {"cycle_number": 5,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.32}},
    ],
  })
  traj = ddm.rfree_trajectory
  assert_false(traj[0].is_retreat_point)
  assert_true(traj[1].is_retreat_point)
  assert_false(traj[2].is_retreat_point)
  print("  PASS: test_retreat_markers")


# ── Phase outcomes ─────────────────────────────

def test_stage_outcomes():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  po = ddm.stage_outcomes
  assert_equal(len(po), 3)
  assert_equal(po[0].stage_id, "data_assessment")
  assert_equal(po[0].status, "complete")
  assert_equal(po[0].description, "Analyze data")
  assert_equal(po[2].cycles_used, 3)
  print("  PASS: test_stage_outcomes")


# ── Data summary ───────────────────────────────

def test_data_summary_mr():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  ds = ddm.data_summary
  assert_equal(ds["Resolution"], "2.00 Å")
  assert_equal(ds["Space group"], "P212121")
  assert_true("TFZ=14.2" in ds["MR solution"])
  assert_true("LLG=342" in ds["MR solution"])
  print("  PASS: test_data_summary_mr")


def test_data_summary_sad():
  ddm = DisplayDataModel.from_session({
    "structure_model": {
      "_version": 1,
      "model_state": {},
      "data_characteristics": {
        "resolution": 2.5,
        "phasing_fom": 0.42,
        "phasing_bayes_cc": 54.2,
        "sites_found": 8,
      },
    },
    "cycles": [],
  })
  ds = ddm.data_summary
  assert_true("FOM=0.420" in ds["Phasing"])
  assert_true("8 sites" in ds["Phasing"])
  print("  PASS: test_data_summary_sad")


# ── Output files ───────────────────────────────

def test_output_files_string():
  ddm = DisplayDataModel.from_session({
    "best_files": {
      "model": "/path/to/model.pdb",
      "map_coeffs_mtz": "/path/to/map.mtz",
    },
    "cycles": [],
  })
  of = ddm.output_files
  assert_equal(of["model_path"], "/path/to/model.pdb")
  assert_equal(of["map_path"], "/path/to/map.mtz")
  print("  PASS: test_output_files_string")


def test_output_files_list():
  """BestFilesTracker may store lists."""
  ddm = DisplayDataModel.from_session({
    "best_files": {
      "model": ["/path/a.pdb", "/path/b.pdb"],
    },
    "cycles": [],
  })
  of = ddm.output_files
  # Should take last (most recent)
  assert_equal(of["model_path"], "/path/b.pdb")
  print("  PASS: test_output_files_list")


# ── Format helpers ─────────────────────────────

def test_format_cycle_compact():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  tl = ddm.timeline
  s = ddm.format_cycle_compact(tl[0])
  assert_true("xtriage" in s)
  assert_true("Res" in s)
  s2 = ddm.format_cycle_compact(tl[4])
  assert_true("→" in s2)  # delta arrow
  assert_true("0.231" in s2)
  print("  PASS: test_format_cycle_compact")


def test_format_outcome_block():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  rows = ddm.format_outcome_block()
  labels = [r[0] for r in rows]
  assert_true("Status" in labels)
  assert_true("R-free" in labels)
  assert_true("Model" in labels)
  status_val = [v for l, v in rows
                if l == "Status"][0]
  assert_true("Determined" in status_val)
  print("  PASS: test_format_outcome_block")


# ── Run stats ──────────────────────────────────

def test_run_stats():
  ddm = DisplayDataModel.from_session(
    _make_mr_session())
  stats = ddm.run_stats
  assert_equal(stats["n_cycles"], 5)
  assert_equal(stats["n_programs"], 3)
  print("  PASS: test_run_stats")


def test_stop_reason():
  ddm = DisplayDataModel.from_session({
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.44}},
      {"cycle_number": 2,
       "program": "STOP",
       "result": "STOP: stuck",
       "reasoning": "R-free stuck. Try different model."},
    ],
    "structure_model": {
      "_version": 1,
      "model_state": {"r_free": 0.44},
      "data_characteristics": {},
    },
  })
  assert_true("stuck" in ddm.stop_reason)
  # No stop → empty
  ddm2 = DisplayDataModel.from_session(
    _make_mr_session())
  assert_equal(ddm2.stop_reason, "")
  print("  PASS: test_stop_reason")


# ── HTML report ────────────────────────────────

def test_html_report_determined():
  try:
    from knowledge.html_report_template import (
      generate_html_report,
    )
  except ImportError:
    print("  SKIP: test_html_report_determined "
          "(html_report_template not available)")
    return
  html = generate_html_report(_make_mr_session())
  assert_true(len(html) > 100)
  assert_true("Structure Determined" in html)
  assert_true("0.231" in html)
  assert_true("P212121" in html)
  assert_true("svg" in html.lower())
  assert_true("polyline" in html.lower())
  print("  PASS: test_html_report_determined")


def test_html_report_stopped():
  try:
    from knowledge.html_report_template import (
      generate_html_report,
    )
  except ImportError:
    print("  SKIP: test_html_report_stopped")
    return
  html = generate_html_report({
    "structure_model": {
      "_version": 1,
      "model_state": {"r_free": 0.45},
      "data_characteristics": {},
    },
    "cycles": [
      {"cycle_number": 1,
       "program": "STOP",
       "result": "STOP: stuck",
       "reasoning": "R-free stuck above 0.40."},
    ],
  })
  assert_true("Agent Stopped" in html)
  assert_true("stuck" in html)
  print("  PASS: test_html_report_stopped")


def test_html_report_retreat_markers():
  try:
    from knowledge.html_report_template import (
      generate_html_report,
    )
  except ImportError:
    print("  SKIP: test_html_report_retreat_markers")
    return
  html = generate_html_report({
    "structure_model": {
      "_version": 1,
      "model_state": {"r_free": 0.28},
      "data_characteristics": {},
      "progress": [],
      "strategy_blacklist": [
        {"retreat_cycle": 3,
         "strategy_id": "s1",
         "reason": "stuck"},
      ],
    },
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.40}},
      {"cycle_number": 2,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.38}},
      {"cycle_number": 3,
       "program": "phenix.autobuild",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.33}},
      {"cycle_number": 4,
       "program": "phenix.refine",
       "result": "SUCCESS: OK",
       "metrics": {"r_free": 0.28}},
    ],
  })
  assert_true("retreat" in html)
  assert_true("stroke-dasharray" in html)
  print("  PASS: test_html_report_retreat_markers")


def test_html_report_cryoem():
  try:
    from knowledge.html_report_template import (
      generate_html_report,
    )
  except ImportError:
    print("  SKIP: test_html_report_cryoem")
    return
  html = generate_html_report({
    "experiment_type": "cryoem",
    "structure_model": {
      "_version": 1,
      "model_state": {"model_map_cc": 0.82},
      "data_characteristics": {
        "experiment_type": "cryoem",
        "resolution": 3.2,
      },
    },
    "cycles": [
      {"cycle_number": 1,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.65}},
      {"cycle_number": 2,
       "program": "phenix.real_space_refine",
       "result": "SUCCESS: OK",
       "metrics": {"model_map_cc": 0.82}},
    ],
  })
  assert_true("Map CC" in html)
  assert_true("polyline" in html)
  print("  PASS: test_html_report_cryoem")


def test_html_report_empty():
  try:
    from knowledge.html_report_template import (
      generate_html_report,
    )
  except ImportError:
    print("  SKIP: test_html_report_empty")
    return
  html = generate_html_report({})
  assert_true("Incomplete" in html)
  print("  PASS: test_html_report_empty")


# ── Runner ─────────────────────────────────────

def run_tests():
  print("Testing display_data_model...")
  test_fmt_rfree()
  test_fmt_pct()
  test_fmt_float()
  test_format_metric_value()
  test_extract_key_metric_refine()
  test_extract_key_metric_rsr()
  test_extract_key_metric_phaser()
  test_extract_key_metric_autosol()
  test_extract_key_metric_ligandfit()
  test_extract_key_metric_no_metrics()
  test_extract_key_metric_fallback_rfree()
  test_extract_key_metric_fallback_cc()
  test_from_session_none()
  test_from_session_empty()
  test_outcome_determined()
  test_outcome_stopped()
  test_outcome_cryoem()
  test_final_metrics_full()
  test_final_metrics_fallback_no_sm()
  test_timeline_basic()
  test_timeline_rfree_deltas()
  test_timeline_mapcc_deltas()
  test_rfree_trajectory()
  test_cc_trajectory()
  test_primary_trajectory_xray()
  test_primary_trajectory_cryoem()
  test_retreat_markers()
  test_stage_outcomes()
  test_data_summary_mr()
  test_data_summary_sad()
  test_output_files_string()
  test_output_files_list()
  test_format_cycle_compact()
  test_format_outcome_block()
  test_run_stats()
  test_stop_reason()
  test_html_report_determined()
  test_html_report_stopped()
  test_html_report_retreat_markers()
  test_html_report_cryoem()
  test_html_report_empty()
  print("\nAll display_data_model tests passed "
        "(41 tests)")


if __name__ == "__main__":
  run_tests()
