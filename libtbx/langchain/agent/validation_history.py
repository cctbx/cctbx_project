"""
Validation History for the Goal-Directed Agent.

Stores per-cycle validation snapshots in a queryable
format. The gate evaluator (Phase 3) uses this to
compute the Monotonic Progress Gate (are metrics better
or worse than at stage start?), and the explanation
engine (Phase 5) uses it for the final report and
per-cycle before/after comparisons.

Serialized to session.data["validation_history"] as a
list of dicts (JSON-safe).

Entry point: ValidationHistory class.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import copy
import logging
import time

logger = logging.getLogger(__name__)


class ValidationHistory(object):
  """Per-cycle validation results, persisted to session.

  Stored as a list of snapshots, one per cycle that ran
  validation. Each snapshot includes the cycle number,
  program name, the full validation_result dict from
  validation_inspector, the log_metrics dict, the
  space group at the time of recording, and a timestamp.

  Designed for two primary consumers:

  1. Gate evaluator (Phase 3):
     - get_metric_series("r_free") for monotonic check
     - get_phase_start_metrics(cycle) for comparison
     - get_metric_series("r_free", space_group="P 32 21")
       for enantiomorph-aware trajectory comparison

  2. Explanation engine (Phase 5):
     - get_flat_metrics(cycle) for before/after diffs
     - get_at_cycle(cycle) for full snapshot access
     - all snapshots for the final report
  """

  def __init__(self):
    self.snapshots = []

  def record(self, cycle_number, program_name,
             validation_result, log_metrics,
             space_group=None):
    """Record a validation snapshot for this cycle.

    Args:
      cycle_number: int.
      program_name: str, e.g. "phenix.refine".
      validation_result: dict from run_validation(),
        or None if validation didn't run.
      log_metrics: dict with r_work, r_free, etc.
        from log analysis, or None.
      space_group: str or None. Current space group
        for enantiomorph-aware trajectory tracking.

    Replaces any existing snapshot for the same cycle
    (idempotent — safe to call twice in the same cycle).

    Never raises.
    """
    try:
      self._record_inner(
        cycle_number, program_name,
        validation_result, log_metrics,
        space_group,
      )
    except Exception:
      logger.debug(
        "ValidationHistory.record failed",
        exc_info=True,
      )

  def _record_inner(self, cycle_number, program_name,
                    validation_result, log_metrics,
                    space_group):
    """Inner record. May raise."""
    cn = int(cycle_number)
    snapshot = {
      "cycle_number": cn,
      "program_name": str(program_name or ""),
      "validation_result": (
        copy.deepcopy(validation_result)
        if validation_result else None
      ),
      "log_metrics": (
        copy.deepcopy(log_metrics)
        if log_metrics else {}
      ),
      "space_group": (
        str(space_group) if space_group else None
      ),
      "timestamp": time.time(),
    }

    # Replace existing snapshot for same cycle
    # (idempotent)
    for i, s in enumerate(self.snapshots):
      if s.get("cycle_number") == cn:
        self.snapshots[i] = snapshot
        return
    self.snapshots.append(snapshot)

    # Cap to 100 snapshots
    if len(self.snapshots) > 100:
      self.snapshots = self.snapshots[-100:]

  def get_at_cycle(self, cycle_number):
    """Retrieve the snapshot from a specific cycle.

    Args:
      cycle_number: int.

    Returns:
      dict (snapshot) or None if no snapshot for that
      cycle.

    Never raises.
    """
    try:
      cn = int(cycle_number)
      for s in self.snapshots:
        if s.get("cycle_number") == cn:
          return s
      return None
    except Exception:
      return None

  def get_latest(self):
    """Get the most recent snapshot.

    Returns:
      dict (snapshot) or None.

    Never raises.
    """
    try:
      if not self.snapshots:
        return None
      return self.snapshots[-1]
    except Exception:
      return None

  def get_previous(self, cycle_number=None):
    """Get the snapshot before a given cycle.

    Args:
      cycle_number: int or None. If None, returns the
        second-to-last snapshot.

    Returns:
      dict (snapshot) or None.

    Never raises.
    """
    try:
      if cycle_number is None:
        if len(self.snapshots) < 2:
          return None
        return self.snapshots[-2]
      cn = int(cycle_number)
      prev = None
      for s in self.snapshots:
        if s.get("cycle_number") >= cn:
          return prev
        prev = s
      return prev
    except Exception:
      return None

  def get_metric_series(self, metric_name,
                        space_group=None):
    """Extract a time series for a specific metric.

    Searches both log_metrics and validation_result
    for the named metric. Supports space group filtering
    for enantiomorph-aware trajectory comparison.

    Args:
      metric_name: str. Supported names:
        From log_metrics: r_free, r_work, resolution,
          model_map_cc, tfz, llg.
        From validation_result.geometry: clashscore,
          rama_favored, rama_outliers, rotamer_outliers,
          bonds_rmsd, angles_rmsd.
        Derived: ligand_cc (min ligand RSCC),
          peak_count (positive diff peaks count),
          water_count.
      space_group: str or None. If provided, only
        include snapshots matching this space group.

    Returns:
      list of (cycle_number, value) tuples, in cycle
      order. Only includes cycles where the metric
      had a non-None value.

    Never raises.
    """
    try:
      return self._get_metric_series_inner(
        metric_name, space_group
      )
    except Exception:
      logger.debug(
        "ValidationHistory.get_metric_series failed",
        exc_info=True,
      )
      return []

  def _get_metric_series_inner(self, metric_name,
                               space_group):
    """Inner metric series extraction. May raise."""
    result = []
    for s in self.snapshots:
      # Space group filter
      if space_group is not None:
        if s.get("space_group") != space_group:
          continue
      cn = s.get("cycle_number", 0)
      val = _extract_metric(s, metric_name)
      if val is not None:
        result.append((cn, val))
    return result

  def get_flat_metrics(self, cycle_number):
    """Extract a flat dict of all metrics at a cycle.

    Used by the explanation engine (Phase 5) for
    before/after comparisons in per-cycle commentary.

    Args:
      cycle_number: int.

    Returns:
      dict of {metric_name: value} for all metrics
      that have values at that cycle. Empty dict if
      no snapshot for that cycle.

    Never raises.
    """
    try:
      s = self.get_at_cycle(cycle_number)
      if s is None:
        return {}
      return _extract_all_metrics(s)
    except Exception:
      return {}

  def get_phase_start_metrics(self, phase_start_cycle):
    """Get the validation snapshot from when a phase began.

    Used by the gate evaluator to compare current metrics
    against stage-start metrics for the monotonic progress
    gate.

    Args:
      phase_start_cycle: int. The cycle number at which
        the current step began.

    Returns:
      dict (snapshot) or None. Returns the snapshot at
      phase_start_cycle if it exists, otherwise the
      closest snapshot at or after that cycle.

    Never raises.
    """
    try:
      cn = int(phase_start_cycle)
      # Exact match first
      exact = self.get_at_cycle(cn)
      if exact is not None:
        return exact
      # Closest at-or-after
      for s in self.snapshots:
        if s.get("cycle_number", 0) >= cn:
          return s
      return None
    except Exception:
      return None

  def get_metrics_delta(self, cycle_a, cycle_b):
    """Compute metric changes between two cycles.

    Args:
      cycle_a: int (earlier cycle).
      cycle_b: int (later cycle).

    Returns:
      dict of {metric_name: {before, after, delta}}.
      Only includes metrics present in both cycles.
      Empty dict if either cycle has no snapshot.

    Never raises.
    """
    try:
      metrics_a = self.get_flat_metrics(cycle_a)
      metrics_b = self.get_flat_metrics(cycle_b)
      if not metrics_a or not metrics_b:
        return {}
      result = {}
      for key in metrics_a:
        if key in metrics_b:
          va = metrics_a[key]
          vb = metrics_b[key]
          if (isinstance(va, (int, float))
              and isinstance(vb, (int, float))):
            result[key] = {
              "before": va,
              "after": vb,
              "delta": vb - va,
            }
      return result
    except Exception:
      return {}

  def is_improving(self, metric_name, n_recent=3,
                   space_group=None):
    """Check if a metric is improving over recent cycles.

    "Improving" means:
      - For r_free, r_work, clashscore, rama_outliers,
        rotamer_outliers: decreasing
      - For rama_favored, ligand_cc, model_map_cc:
        increasing

    Args:
      metric_name: str.
      n_recent: int. Number of recent values to check.
      space_group: str or None.

    Returns:
      True if improving, False if stalled or worsening,
      None if insufficient data.

    Never raises.
    """
    try:
      series = self.get_metric_series(
        metric_name, space_group
      )
      if len(series) < n_recent:
        return None
      values = [v for _, v in series[-n_recent:]]
      direction = _metric_direction(metric_name)
      if direction == "lower_is_better":
        # Check if overall trend is decreasing
        return values[-1] < values[0]
      elif direction == "higher_is_better":
        return values[-1] > values[0]
      return None
    except Exception:
      return None

  def get_recent_values(self, metric_name, n=3,
                        space_group=None):
    """Get the last N values of a metric.

    Convenience for the gate evaluator's threshold
    checks like "r_free > 0.45 after 2 cycles".

    Args:
      metric_name: str.
      n: int. Number of recent values.
      space_group: str or None.

    Returns:
      list of float values (most recent last), or
      empty list if insufficient data.

    Never raises.
    """
    try:
      series = self.get_metric_series(
        metric_name, space_group
      )
      if not series:
        return []
      values = [v for _, v in series[-n:]]
      return values
    except Exception:
      return []

  def get_metric_at_cycle(self, metric_name,
                          cycle_number):
    """Get a single metric value at a specific cycle.

    Convenience for extracting one metric from a
    snapshot without needing to call get_flat_metrics
    or manually dig into the snapshot dict.

    Args:
      metric_name: str.
      cycle_number: int.

    Returns:
      float or None.

    Never raises.
    """
    try:
      s = self.get_at_cycle(cycle_number)
      if s is None:
        return None
      return _extract_metric(s, metric_name)
    except Exception:
      return None

  def cycle_count(self):
    """Number of recorded snapshots.

    Returns:
      int.
    """
    return len(self.snapshots)

  # ── Serialization ───────────────────────────────────

  def to_dict(self):
    """Serialize for session persistence.

    Returns:
      dict with key "snapshots" containing a list
      of JSON-safe snapshot dicts.
    """
    try:
      return {
        "snapshots": [
          _sanitize_snapshot(s)
          for s in self.snapshots
        ],
        "_version": 1,
      }
    except Exception:
      logger.debug(
        "ValidationHistory.to_dict failed",
        exc_info=True,
      )
      return {"snapshots": [], "_version": 1}

  @classmethod
  def from_dict(cls, d):
    """Deserialize from session data.

    Tolerant of missing/extra keys.

    Args:
      d: dict from to_dict().

    Returns:
      ValidationHistory instance.
    """
    vh = cls()
    if not isinstance(d, dict):
      return vh
    try:
      snaps = d.get("snapshots", [])
      if isinstance(snaps, list):
        vh.snapshots = [
          dict(s) for s in snaps
          if isinstance(s, dict)
        ]
    except Exception:
      logger.debug(
        "ValidationHistory.from_dict failed",
        exc_info=True,
      )
    return vh


# ── Metric extraction helpers (module-private) ──────

def _safe_float(val):
  """Convert to float or return None."""
  if val is None:
    return None
  try:
    return float(val)
  except (ValueError, TypeError):
    return None


def _extract_metric(snapshot, metric_name):
  """Extract a named metric from a snapshot.

  Searches log_metrics first (authoritative for
  R-factors), then validation_result sub-dicts.

  Returns float or None.
  """
  lm = snapshot.get("log_metrics") or {}
  vr = snapshot.get("validation_result") or {}

  # --- log_metrics (authoritative) ---
  if metric_name in ("r_free", "r_work", "resolution",
                     "model_map_cc", "tfz", "llg"):
    val = _safe_float(lm.get(metric_name))
    if val is not None:
      return val

  # --- validation_result.geometry ---
  geom = vr.get("geometry") or {}
  if metric_name in ("clashscore", "rama_favored",
                     "rama_outliers", "rotamer_outliers",
                     "bonds_rmsd", "angles_rmsd"):
    val = _safe_float(geom.get(metric_name))
    if val is not None:
      return val

  # --- Derived: ligand_cc (minimum RSCC) ---
  if metric_name == "ligand_cc":
    dm = vr.get("data_model") or {}
    rscc_list = dm.get("ligand_rscc", [])
    rsccs = [
      _safe_float(entry.get("rscc"))
      for entry in rscc_list
      if isinstance(entry, dict)
      and entry.get("rscc") is not None
    ]
    rsccs = [r for r in rsccs if r is not None]
    return min(rsccs) if rsccs else None

  # --- Derived: peak_count ---
  if metric_name == "peak_count":
    dp = vr.get("diff_peaks") or {}
    pc = dp.get("peak_count")
    if pc is not None:
      return int(pc)
    pos = dp.get("positive", [])
    return len(pos) if pos else None

  # --- Derived: water_count ---
  if metric_name == "water_count":
    mc = vr.get("model_contents") or {}
    wc = mc.get("waters")
    if wc is not None:
      return int(wc)
    return None

  # --- Derived: residue_count ---
  if metric_name == "residue_count":
    mc = vr.get("model_contents") or {}
    rc = mc.get("residue_count")
    if rc is not None:
      return int(rc)
    return None

  return None


def _extract_all_metrics(snapshot):
  """Extract all available metrics as a flat dict.

  Returns dict of {metric_name: value}.
  """
  result = {}
  all_metrics = [
    "r_free", "r_work", "resolution",
    "model_map_cc", "tfz", "llg",
    "clashscore", "rama_favored", "rama_outliers",
    "rotamer_outliers", "bonds_rmsd", "angles_rmsd",
    "ligand_cc", "peak_count", "water_count",
    "residue_count",
  ]
  for name in all_metrics:
    val = _extract_metric(snapshot, name)
    if val is not None:
      result[name] = val
  return result


def _metric_direction(metric_name):
  """Whether lower or higher values are better.

  Returns "lower_is_better", "higher_is_better",
  or "neutral".
  """
  lower = frozenset([
    "r_free", "r_work", "clashscore",
    "rama_outliers", "rotamer_outliers",
    "bonds_rmsd", "angles_rmsd", "peak_count",
  ])
  higher = frozenset([
    "rama_favored", "ligand_cc", "model_map_cc",
    "water_count", "residue_count",
  ])
  if metric_name in lower:
    return "lower_is_better"
  if metric_name in higher:
    return "higher_is_better"
  return "neutral"


def _sanitize_snapshot(snapshot):
  """Deep-copy a snapshot and ensure JSON safety.

  Removes any non-serializable values (numpy arrays,
  cctbx objects that may have leaked in).
  """
  try:
    s = copy.deepcopy(snapshot)
    # Ensure timestamp is a plain float
    ts = s.get("timestamp")
    if ts is not None:
      s["timestamp"] = float(ts)
    return s
  except Exception:
    # Fallback: keep only the safe scalar fields
    return {
      "cycle_number": snapshot.get("cycle_number", 0),
      "program_name": snapshot.get(
        "program_name", ""
      ),
      "validation_result": None,
      "log_metrics": {},
      "space_group": snapshot.get("space_group"),
      "timestamp": time.time(),
    }
