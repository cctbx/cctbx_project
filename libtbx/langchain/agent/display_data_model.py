"""
Display Data Model for PHENIX AI Agent (v114).

Single data provider for all display views:
  - Progress tab (live monitoring)
  - Results tab (post-hoc review)
  - HTML report (shareable document)

All three views read from this model. Metric
formatting, status indicators, file paths, and
stage descriptions are computed once, here.

No imports from phenix.*, wx.*, or langchain.*.
Pure data transformation — safe for agent/ directory.
"""

from __future__ import (
  absolute_import, division, print_function,
)

import logging
import os
import time

logger = logging.getLogger(__name__)


# ── Named tuples / data classes ─────────────────

class CycleEntry(object):
  """One row in the timeline table."""

  __slots__ = (
    "cycle", "program", "status",
    "key_metric_name", "key_metric_value",
    "key_metric_delta", "reasoning",
    "command", "expert_assessment",
    "is_retreat_point",
  )

  def __init__(self, cycle, program, status,
               key_metric_name="",
               key_metric_value=None,
               key_metric_delta=None,
               reasoning="", command="",
               expert_assessment=None,
               is_retreat_point=False):
    self.cycle = cycle
    self.program = program
    self.status = status
    self.key_metric_name = key_metric_name
    self.key_metric_value = key_metric_value
    self.key_metric_delta = key_metric_delta
    self.reasoning = reasoning
    self.command = command
    self.expert_assessment = (
      expert_assessment or {})
    self.is_retreat_point = is_retreat_point

  def to_dict(self):
    return {
      "cycle": self.cycle,
      "program": self.program,
      "status": self.status,
      "key_metric_name": self.key_metric_name,
      "key_metric_value": self.key_metric_value,
      "key_metric_delta": self.key_metric_delta,
      "reasoning": self.reasoning,
      "command": self.command,
      "is_retreat_point": self.is_retreat_point,
    }


class PhaseEntry(object):
  """One row in the stage outcomes list."""

  __slots__ = (
    "stage_id", "description", "status",
    "programs", "cycles_used", "key_metric",
  )

  def __init__(self, stage_id, description, status,
               programs=None, cycles_used=0,
               key_metric=""):
    self.stage_id = stage_id
    self.description = description
    self.status = status
    self.programs = programs or []
    self.cycles_used = cycles_used
    self.key_metric = key_metric


class RfreePoint(object):
  """One point in the R-free trajectory."""

  __slots__ = (
    "cycle", "value", "program",
    "is_retreat_point",
  )

  def __init__(self, cycle, value, program="",
               is_retreat_point=False):
    self.cycle = cycle
    self.value = value
    self.program = program
    self.is_retreat_point = is_retreat_point


# ── Metric formatting helpers ───────────────────

def _fmt_rfree(val):
  """Format R-free/R-work as '0.231'."""
  if val is None:
    return None
  try:
    return "%.3f" % float(val)
  except (ValueError, TypeError):
    return None


def _fmt_pct(val):
  """Format a fraction as percentage '97.4%'."""
  if val is None:
    return None
  try:
    v = float(val)
    # If already in percentage form (>1.0)
    if v > 1.0:
      return "%.1f%%" % v
    return "%.1f%%" % (v * 100)
  except (ValueError, TypeError):
    return None


def _fmt_float(val, decimals=1):
  """Format a float to N decimal places."""
  if val is None:
    return None
  try:
    return ("%%.%df" % decimals) % float(val)
  except (ValueError, TypeError):
    return None


def _safe_float(val):
  """Convert to float or return None."""
  if val is None:
    return None
  try:
    return float(val)
  except (ValueError, TypeError):
    return None


def _format_metric_value(name, value):
  """Format a metric value with appropriate precision
  based on the metric name.

  R-free/R-work: 3 decimals (0.231)
  Map CC, FOM, Ligand CC: 3 decimals (0.820)
  TFZ, LLG, Clash, BAYES-CC: 1 decimal (14.2)
  Resolution: 2 decimals (2.50)
  """
  _3_DECIMAL = {
    "R-free", "R-work", "Map CC", "FOM",
    "Ligand CC",
  }
  _2_DECIMAL = {"Res",}
  if name in _3_DECIMAL:
    return _fmt_float(value, 3)
  elif name in _2_DECIMAL:
    return _fmt_float(value, 2)
  else:
    return _fmt_float(value, 1)


# ── Key metric selection per program ────────────

_PROGRAM_KEY_METRICS = {
  "phenix.xtriage": [
    ("resolution", "Res"),
  ],
  "phenix.mtriage": [
    ("resolution", "Res"),
  ],
  "phenix.phaser": [
    ("tfz", "TFZ"),
    ("llg", "LLG"),
  ],
  "phenix.autosol": [
    ("fom", "FOM"),
    ("bayes_cc", "BAYES-CC"),
  ],
  "phenix.refine": [
    ("r_free", "R-free"),
  ],
  "phenix.real_space_refine": [
    ("model_map_cc", "Map CC"),
  ],
  "phenix.autobuild": [
    ("r_free", "R-free"),
  ],
  "phenix.ligandfit": [
    ("cc", "Ligand CC"),
  ],
  "phenix.molprobity": [
    ("clashscore", "Clash"),
  ],
  "phenix.model_vs_data": [
    ("r_free", "R-free"),
  ],
  "phenix.map_to_model": [
    ("model_map_cc", "Map CC"),
  ],
  "phenix.predict_and_build": [
    ("r_free", "R-free"),
    ("model_map_cc", "Map CC"),
  ],
  "phenix.dock_in_map": [
    ("model_map_cc", "Map CC"),
  ],
  "phenix.validation_cryoem": [
    ("model_map_cc", "Map CC"),
  ],
  "phenix.map_correlations": [
    ("model_map_cc", "Map CC"),
  ],
}

# Metrics that should show delta between cycles
_DELTA_METRICS = {"R-free", "Map CC"}


def _extract_key_metric(program, metrics):
  """Extract the most informative metric for a
  program from its result metrics dict.

  Returns (name, value) or ("", None).
  """
  if not metrics or not isinstance(metrics, dict):
    return ("", None)
  prog = (program or "").lower().replace(
    "phenix.", "")
  # Exact match first, then prefix match
  matched = None
  for key in _PROGRAM_KEY_METRICS:
    key_short = key.replace("phenix.", "")
    if key_short == prog:
      matched = _PROGRAM_KEY_METRICS[key]
      break
  if matched is None:
    # Prefix match (e.g. "refine" matches
    # "phenix.refine" but NOT "real_space_refine")
    for key in _PROGRAM_KEY_METRICS:
      key_short = key.replace("phenix.", "")
      if prog == key_short or (
        prog.startswith(key_short + "_")
        or prog.startswith(key_short + ".")
      ):
        matched = _PROGRAM_KEY_METRICS[key]
        break
  if matched is not None:
    for metric_key, display_name in matched:
      val = _safe_float(metrics.get(metric_key))
      if val is not None:
        return (display_name, val)
  # Fallback: r_free if present
  rf = _safe_float(metrics.get("r_free"))
  if rf is not None:
    return ("R-free", rf)
  # Fallback: model_map_cc if present
  cc = _safe_float(metrics.get("model_map_cc"))
  if cc is not None:
    return ("Map CC", cc)
  return ("", None)


# ── Main class ──────────────────────────────────

class DisplayDataModel(object):
  """Single provider for all display views.

  Takes raw session data and produces a unified
  display object with pre-formatted fields.

  All three views (Progress, Results, HTML) read
  from this model. Metric formatting, status
  indicators, file paths, and stage descriptions
  are computed once, here.

  Never raises — all properties return safe
  defaults on missing data.
  """

  def __init__(self):
    # Raw data (set by from_session)
    self._sm_data = {}
    self._plan_data = {}
    self._cycles = []
    self._best_files = {}
    self._original_files = []
    self._directives = {}
    self._start_time = None
    self._experiment_type = ""

  @classmethod
  def from_session(cls, session_data):
    """Construct from a session data dict.

    Args:
      session_data: dict from session.to_dict()
        or session.data. Must have at least
        'cycles' key.

    Returns:
      DisplayDataModel instance.

    Never raises.
    """
    ddm = cls()
    if not isinstance(session_data, dict):
      return ddm
    try:
      ddm._sm_data = session_data.get(
        "structure_model", {}) or {}
      ddm._plan_data = session_data.get(
        "plan", {}) or {}
      ddm._cycles = session_data.get(
        "cycles", []) or []
      ddm._best_files = session_data.get(
        "best_files", {}) or {}
      ddm._original_files = session_data.get(
        "original_files", []) or []
      ddm._directives = session_data.get(
        "directives", {}) or {}
      ddm._start_time = session_data.get(
        "start_time")
      ddm._experiment_type = session_data.get(
        "experiment_type", "") or ""
      # Fallback: structure model may have it
      if not ddm._experiment_type:
        dc = (ddm._sm_data.get(
          "data_characteristics") or {})
        ddm._experiment_type = dc.get(
          "experiment_type", "") or ""
    except Exception:
      logger.debug(
        "DisplayDataModel.from_session failed",
        exc_info=True,
      )
    return ddm

  # ── Outcome ─────────────────────────────────

  @property
  def outcome_status(self):
    """'determined', 'stopped', 'incomplete'.

    'determined' = agent completed all stages or
      reached target metrics.
    'stopped' = agent stopped before completion
      (gate_stop, red_flag, user request, etc.)
    'incomplete' = no cycles ran or no metrics.
    """
    if not self._cycles:
      return "incomplete"
    # Check plan completion
    stages = self._plan_data.get("stages") or self._plan_data.get("phases") or []
    if stages:
      all_done = all(
        p.get("status") in ("complete", "skipped")
        for p in stages
        if isinstance(p, dict)
      )
      if all_done:
        return "determined"
    # Check if last cycle was a stop
    last = self._cycles[-1] if self._cycles else {}
    result = str(last.get("result", "")).upper()
    program = str(last.get("program", "")).upper()
    if program == "STOP" or "STOP" in result:
      # Stopped — but check if we have good metrics
      # (agent may have stopped because it's done)
      rf = self._get_model_metric("r_free")
      cc = self._get_model_metric("model_map_cc")
      if cc is None:
        cc = self._get_model_metric("map_cc")
      if rf is not None and rf < 0.30:
        return "determined"
      if cc is not None and cc > 0.7:
        return "determined"
      return "stopped"
    if "FAILED" in result:
      return "stopped"
    # Has good metrics?
    rf = self._get_model_metric("r_free")
    cc = self._get_model_metric("model_map_cc")
    if cc is None:
      cc = self._get_model_metric("map_cc")
    if rf is not None and rf < 0.30:
      return "determined"
    if cc is not None and cc > 0.7:
      return "determined"
    # Has any metrics at all? Then it ran but
    # didn't reach targets — "stopped" not
    # "incomplete".
    if rf is not None or cc is not None:
      return "stopped"
    return "incomplete"

  @property
  def outcome_message(self):
    """Human-readable one-line summary."""
    status = self.outcome_status
    rf = self._get_model_metric("r_free")
    cc = self._get_model_metric("model_map_cc")
    if cc is None:
      cc = self._get_model_metric("map_cc")

    if status == "determined":
      if rf is not None:
        return (
          "Structure determined (R-free: %s)"
          % _fmt_rfree(rf)
        )
      if cc is not None:
        return (
          "Structure determined (Map CC: %s)"
          % _fmt_float(cc, 3)
        )
      return "Structure determination complete"
    elif status == "stopped":
      if rf is not None:
        return (
          "Agent stopped — R-free: %s"
          % _fmt_rfree(rf)
        )
      if cc is not None:
        return (
          "Agent stopped — Map CC: %s"
          % _fmt_float(cc, 3)
        )
      return "Agent stopped"
    else:
      return "Incomplete — no results yet"

  @property
  def stop_reason(self):
    """Reasoning from the last cycle when the agent
    stopped. Useful for displaying in the Outcome
    section to explain WHY it stopped.

    Returns empty string if not stopped or no
    reasoning available.
    """
    if self.outcome_status != "stopped":
      return ""
    # Check last cycle for reasoning
    for c in reversed(self._cycles):
      if not isinstance(c, dict):
        continue
      result = str(c.get("result", "")).upper()
      reasoning = c.get("reasoning", "")
      if "STOP" in result or "FAILED" in result:
        if reasoning:
          return reasoning.strip()
      # Also check if last cycle has a stop-related
      # program (STOP or failed)
      prog = (c.get("program") or "").upper()
      if prog == "STOP" and reasoning:
        return reasoning.strip()
    return ""

  @property
  def final_metrics(self):
    """Dict of final metric name → formatted value.

    Keys present only when the metric has a value.
    All values are strings, pre-formatted.

    Sources (in priority order):
    1. StructureModel model_state (v114+)
    2. Last successful cycle's metrics (fallback)
    """
    ms = self._get_model_state()
    geom = ms.get("geometry", {})
    dc = self._get_data_chars()
    result = {}

    # R-factors — try SM first, then scan all cycles
    rf = _safe_float(ms.get("r_free"))
    rw = _safe_float(ms.get("r_work"))
    cc = _safe_float(ms.get("model_map_cc"))
    if rf is None:
      rf = self._best_metric_from_cycles("r_free")
    if rw is None:
      rw = self._best_metric_from_cycles("r_work")
    if cc is None:
      cc = self._best_metric_from_cycles(
        "model_map_cc")
      if cc is None:
        cc = self._best_metric_from_cycles("map_cc")

    if rf is not None:
      result["R-free"] = _fmt_rfree(rf)
    if rw is not None:
      result["R-work"] = _fmt_rfree(rw)
    if cc is not None:
      result["Map CC"] = _fmt_float(cc, 3)

    # Geometry — try SM first, then scan cycles
    # (molprobity cycle metrics as fallback)
    cs_val = _safe_float(geom.get("clashscore"))
    if cs_val is None:
      cs_val = self._best_metric_from_cycles(
        "clashscore")
    cs = _fmt_float(cs_val)
    if cs:
      result["Clashscore"] = cs
    rf_pct_val = _safe_float(
      geom.get("rama_favored"))
    if rf_pct_val is None:
      # Try ramachandran_favored from cycle metrics
      rf_pct_val = self._best_metric_from_cycles(
        "ramachandran_favored")
      # Also try ramachandran_outliers (lower = better)
      if rf_pct_val is None:
        ro_val = self._best_metric_from_cycles(
          "ramachandran_outliers")
        if ro_val is not None:
          result["Ramachandran outliers"] = (
            "%.1f%%" % ro_val)
    rf_pct = _fmt_pct(rf_pct_val)
    if rf_pct:
      result["Ramachandran favored"] = rf_pct
    # Rotamer outliers
    ro_pct = _safe_float(
      geom.get("rotamer_outliers"))
    if ro_pct is None:
      ro_pct = self._best_metric_from_cycles(
        "rotamer_outliers")
    if ro_pct is not None:
      result["Rotamer outliers"] = "%.1f%%" % ro_pct
    # MolProbity score
    mp_score = self._best_metric_from_cycles(
      "molprobity_score")
    if mp_score is not None:
      result["MolProbity score"] = _fmt_float(
        mp_score, 2)

    # Model contents
    chains = ms.get("chains", [])
    if chains:
      result["Chains"] = str(len(chains))
    waters = ms.get("waters", 0)
    if waters:
      result["Waters"] = str(waters)
    ligands = ms.get("ligands", [])
    if ligands:
      lig_names = [
        l.get("name", "?")
        for l in ligands
        if isinstance(l, dict)
      ]
      result["Ligands"] = ", ".join(lig_names)

    # Data characteristics
    res = dc.get("resolution")
    if res is None:
      res = self._best_metric_from_cycles(
        "resolution")
    if res is not None:
      result["Resolution"] = "%s Å" % _fmt_float(
        res, 2)
    sg = dc.get("space_group")
    if sg:
      result["Space group"] = str(sg)

    # Ligand CC (from ligandfit cycles)
    lig_cc = self._best_metric_from_cycles(
      "ligand_cc")
    if lig_cc is not None:
      result["Ligand CC"] = _fmt_float(lig_cc, 3)

    return result

  @property
  def output_files(self):
    """Dict: model_path, map_path, report_path.

    Values are absolute paths (str) or None.
    Handles list values from BestFilesTracker
    by taking the last (most recent) entry.
    """
    bf = self._best_files

    def _resolve(val):
      if isinstance(val, list):
        return val[-1] if val else None
      return val if val else None

    return {
      "model_path": _resolve(bf.get("model")),
      "map_path": _resolve(
        bf.get("map_coeffs_mtz")
        or bf.get("refine_map_coeffs")
      ),
    }

  @property
  def run_stats(self):
    """Dict: n_cycles, n_programs, duration_str,
    run_date_str."""
    n_cycles = len(self._cycles)
    programs = set()
    for c in self._cycles:
      if isinstance(c, dict):
        p = c.get("program", "")
        if p:
          programs.add(p)

    duration = ""
    if self._start_time:
      try:
        elapsed = time.time() - float(
          self._start_time)
        if elapsed > 3600:
          duration = "%.1f hours" % (
            elapsed / 3600)
        elif elapsed > 60:
          duration = "%d minutes" % (
            elapsed / 60)
        else:
          duration = "%d seconds" % elapsed
      except (ValueError, TypeError):
        pass

    return {
      "n_cycles": n_cycles,
      "n_programs": len(programs),
      "duration": duration,
    }

  # ── Timeline ────────────────────────────────

  @property
  def timeline(self):
    """List of CycleEntry for the timeline table.

    One entry per cycle, with the most informative
    metric extracted from the cycle's result.
    """
    entries = []
    prev_metrics = {}  # metric_name → prev value
    for i, c in enumerate(self._cycles):
      if not isinstance(c, dict):
        continue
      cycle_num = c.get("cycle_number", i + 1)
      program = c.get("program", "?")
      result = str(c.get("result", ""))

      # Status
      if "SUCCESS" in result.upper():
        status = "OK"
      elif "FAILED" in result.upper():
        status = "FAILED"
      elif "STOP" in result.upper():
        status = "STOP"
      else:
        status = "?"

      # Key metric from cycle metrics
      metrics = c.get("metrics", {})
      if not isinstance(metrics, dict):
        metrics = {}
      name, val = _extract_key_metric(
        program, metrics)

      # Compute delta for trackable metrics
      delta = None
      if name in _DELTA_METRICS and val is not None:
        prev = prev_metrics.get(name)
        if prev is not None:
          delta = val - prev
        prev_metrics[name] = val

      entries.append(CycleEntry(
        cycle=cycle_num,
        program=program,
        status=status,
        key_metric_name=name,
        key_metric_value=val,
        key_metric_delta=delta,
        reasoning=c.get("reasoning", ""),
        command=c.get("command", ""),
        expert_assessment=c.get(
          "expert_assessment", {}),
      ))
    return entries

  # ── R-free trajectory ───────────────────────

  @property
  def rfree_trajectory(self):
    """List of RfreePoint for X-ray trajectory chart.

    Includes retreat markers from plan progress.
    """
    return self._build_trajectory("r_free")

  @property
  def cc_trajectory(self):
    """List of RfreePoint for cryo-EM trajectory
    chart (model-map CC)."""
    return self._build_trajectory("model_map_cc")

  @property
  def primary_trajectory(self):
    """The trajectory most relevant for this
    experiment type. R-free for X-ray, Map CC
    for cryo-EM."""
    if self._experiment_type == "cryoem":
      traj = self.cc_trajectory
      if traj:
        return traj
    return self.rfree_trajectory

  def _build_trajectory(self, metric_key):
    """Build a trajectory for a specific metric."""
    points = []
    for i, c in enumerate(self._cycles):
      if not isinstance(c, dict):
        continue
      cycle_num = c.get("cycle_number", i + 1)
      program = c.get("program", "")
      metrics = c.get("metrics", {})
      if not isinstance(metrics, dict):
        continue
      val = _safe_float(metrics.get(metric_key))
      if val is not None:
        points.append(RfreePoint(
          cycle=cycle_num,
          value=val,
          program=program,
        ))

    # Mark retreat points from strategy blacklist
    retreat_cycles = set()
    bl = self._sm_data.get(
      "strategy_blacklist", [])
    for entry in bl:
      if isinstance(entry, dict):
        rc = entry.get("retreat_cycle")
        if rc is not None:
          retreat_cycles.add(int(rc))
    for pt in points:
      if pt.cycle in retreat_cycles:
        pt.is_retreat_point = True

    return points

  # ── Stage outcomes ──────────────────────────

  @property
  def stage_outcomes(self):
    """List of PhaseEntry from the plan."""
    stages = self._plan_data.get("stages") or self._plan_data.get("phases") or []
    entries = []
    for p in stages:
      if not isinstance(p, dict):
        continue
      entries.append(PhaseEntry(
        stage_id=p.get("id", "?"),
        description=p.get(
          "description",
          p.get("id", "?").replace(
            "_", " ").title()),
        status=p.get("status", "pending"),
        programs=p.get("programs", []),
        cycles_used=p.get("cycles_used", 0),
      ))
    return entries

  # ── Data characteristics ────────────────────

  @property
  def data_summary(self):
    """Dict of data characteristics for display.

    All values are pre-formatted strings.
    """
    dc = self._get_data_chars()
    result = {}

    res = dc.get("resolution")
    if res is not None:
      result["Resolution"] = "%s Å" % _fmt_float(
        res, 2)
    sg = dc.get("space_group")
    if sg:
      result["Space group"] = str(sg)
    exp = dc.get("experiment_type")
    if exp:
      result["Experiment"] = str(exp)

    # MR metrics
    tfz = _safe_float(dc.get("mr_tfz"))
    if tfz is not None:
      mr_str = "TFZ=%.1f" % tfz
      llg = _safe_float(dc.get("mr_llg"))
      if llg is not None:
        mr_str += ", LLG=%.0f" % llg
      result["MR solution"] = mr_str

    # Phasing metrics
    fom = _safe_float(dc.get("phasing_fom"))
    if fom is not None:
      ph_str = "FOM=%.3f" % fom
      bcc = _safe_float(
        dc.get("phasing_bayes_cc"))
      if bcc is not None:
        ph_str += ", BAYES-CC=%.1f" % bcc
      sites = dc.get("sites_found")
      if sites is not None:
        ph_str += ", %d sites" % int(sites)
      result["Phasing"] = ph_str

    return result

  # ── Formatting helpers for views ────────────

  def format_cycle_compact(self, entry):
    """Format a CycleEntry as a one-line string
    for the Progress tab.

    Example: 'Cycle 3: phenix.refine — R-free:
    0.350 → 0.312'
    """
    parts = ["Cycle %d: %s" % (
      entry.cycle, entry.program)]
    if entry.status and entry.status != "?":
      parts[0] += " — %s" % entry.status
    if (entry.key_metric_name
        and entry.key_metric_value is not None):
      # Format value based on metric type
      val_str = _format_metric_value(
        entry.key_metric_name,
        entry.key_metric_value,
      )
      metric_str = "%s: %s" % (
        entry.key_metric_name, val_str)
      if (entry.key_metric_delta is not None
          and entry.key_metric_name
          in _DELTA_METRICS):
        prev = (
          entry.key_metric_value
          - entry.key_metric_delta
        )
        metric_str = "%s: %s → %s" % (
          entry.key_metric_name,
          _format_metric_value(
            entry.key_metric_name, prev),
          val_str,
        )
      parts.append(metric_str)
    return "  ".join(parts)

  def format_outcome_block(self):
    """Format the Outcome section as multi-line
    text for the Results tab.

    Returns list of (label, value) tuples.
    """
    rows = []
    status = self.outcome_status
    if status == "determined":
      rows.append(("Status", "✓ Structure Determined"))
    elif status == "stopped":
      rows.append(("Status", "⚠ Agent Stopped"))
    else:
      rows.append(("Status", "● Incomplete"))

    for k, v in self.final_metrics.items():
      rows.append((k, v))

    files = self.output_files
    if files.get("model_path"):
      rows.append((
        "Model",
        os.path.basename(files["model_path"]),
      ))
    if files.get("map_path"):
      rows.append((
        "Map",
        os.path.basename(files["map_path"]),
      ))

    stats = self.run_stats
    stats_str = "%d cycles" % stats["n_cycles"]
    if stats["duration"]:
      stats_str += ", %s" % stats["duration"]
    rows.append(("Run", stats_str))

    return rows

  # ── Private helpers ─────────────────────────

  def _get_model_state(self):
    if not isinstance(self._sm_data, dict):
      return {}
    return self._sm_data.get(
      "model_state", {}) or {}

  def _get_data_chars(self):
    if not isinstance(self._sm_data, dict):
      return {}
    return self._sm_data.get(
      "data_characteristics", {}) or {}

  def _get_model_metric(self, key):
    ms = self._get_model_state()
    val = _safe_float(ms.get(key))
    if val is None:
      # Fallback: scan all successful cycles
      val = self._best_metric_from_cycles(key)
    return val

  def _last_cycle_metrics(self):
    """Get metrics from the last successful cycle.

    Scans cycles in reverse for the first one with
    a SUCCESS result and metrics dict.
    """
    for c in reversed(self._cycles):
      if not isinstance(c, dict):
        continue
      result = str(c.get("result", "")).upper()
      if "SUCCESS" in result:
        m = c.get("metrics", {})
        if isinstance(m, dict) and m:
          return m
    return {}

  def _best_metric_from_cycles(self, key):
    """Find the best value of a metric across all
    successful cycles.

    For metrics in LOWER_IS_BETTER (r_free, r_work),
    returns the minimum. For others, returns the
    last non-None value.

    Args:
      key: str metric name.

    Returns:
      float or None.
    """
    _lower = {"r_free", "r_work", "clashscore",
              "ramachandran_outliers",
              "rotamer_outliers",
              "molprobity_score"}
    best = None
    for c in self._cycles:
      if not isinstance(c, dict):
        continue
      result = str(c.get("result", "")).upper()
      if "SUCCESS" not in result:
        continue
      m = c.get("metrics", {})
      if not isinstance(m, dict):
        continue
      v = _safe_float(m.get(key))
      if v is not None:
        if best is None:
          best = v
        elif key in _lower and v < best:
          best = v
        else:
          best = v  # Last seen for non-lower
    return best
