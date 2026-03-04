"""
Compact text formatter for validation results (v113).

Produces a text block suitable for injection into the
thinking agent's context window alongside the existing
log analysis and strategy memory.

Entry points:
  format_validation_report(...) -> str
  format_validation_trend(r_free_history) -> str

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function


def format_validation_report(
  validation_result,
  log_metrics=None,
  cycle_number=None,
  program_name=None,
  prev_r_free=None,
  start_r_free=None,
  max_chars=500,
):
  """Format validation results as compact text.

  Args:
    validation_result: dict from run_validation(),
      or None if validation failed.
    log_metrics: dict from existing log_analysis,
      contains r_work, r_free, bonds, angles.
    cycle_number: int or None.
    program_name: str or None.
    prev_r_free: float or None (from previous cycle).
    start_r_free: float or None (from first cycle).
    max_chars: int, soft limit on report length.

  Returns:
    str. Always returns a string, possibly just a
    header with "(not available)".
  """
  lines = []

  # Header
  header = "=== VALIDATION"
  if cycle_number is not None:
    header += " (cycle %d" % cycle_number
    if program_name:
      header += ": %s" % program_name
    header += ")"
  header += " ==="
  lines.append(header)

  # R-factors from log_metrics (authoritative source)
  if log_metrics:
    r_line = _format_r_factors(
      log_metrics, prev_r_free, start_r_free
    )
    if r_line:
      lines.append(r_line)

  if validation_result is None:
    lines.append("(detailed validation not available)")
    return "\n".join(lines)

  # Geometry
  geom = validation_result.get("geometry")
  if geom:
    lines.extend(_format_geometry(geom, log_metrics))

  # Model contents
  contents = validation_result.get("model_contents")
  if contents:
    lines.extend(_format_model_contents(contents))

  # Ligand RSCC (when available)
  dm = validation_result.get("data_model", {})
  rscc_list = dm.get("ligand_rscc", [])
  if rscc_list:
    lines.extend(_format_ligand_rscc(rscc_list))

  # Difference peaks (when available)
  peaks = validation_result.get("diff_peaks")
  if peaks:
    lines.extend(_format_diff_peaks(peaks))

  report = "\n".join(lines)

  # Truncate if too long
  if len(report) > max_chars:
    report = report[:max_chars - 20] + "\n[truncated]"

  return report


def format_validation_trend(r_free_history):
  """One-line R-free trend from cycle history.

  Args:
    r_free_history: list of float (R-free per cycle).

  Returns:
    str, e.g. "R-free trend: 0.421 -> 0.295 -> 0.248
    (3 cycles)" or "" if insufficient data.
  """
  if not r_free_history or len(r_free_history) < 2:
    return ""
  formatted = " -> ".join(
    "%.3f" % rf for rf in r_free_history
  )
  return "R-free trend: %s (%d cycles)" % (
    formatted, len(r_free_history)
  )


# =========================================================
# Internal formatters
# =========================================================

def _format_r_factors(metrics, prev_rf, start_rf):
  """Format R-factor line with trend."""
  rw = metrics.get("r_work")
  rf = metrics.get("r_free")
  if rw is None and rf is None:
    return None
  parts = []
  if rw is not None:
    parts.append("R-work=%.3f" % rw)
  if rf is not None:
    parts.append("R-free=%.3f" % rf)
  trend = ""
  if rf is not None and prev_rf is not None:
    trend += " (prev: %.3f" % prev_rf
    if start_rf is not None:
      trend += ", start: %.3f" % start_rf
    trend += ")"
  return " ".join(parts) + trend


def _format_geometry(geom, log_metrics=None):
  """Format geometry metrics.

  Prefers bonds/angles from log_metrics (refinement
  log) over computed values, since those are the
  authoritative source.
  """
  lines = []

  # Bonds/angles: prefer log_metrics, fall back to
  # computed values from the restraints manager.
  # NOTE: log_parsers uses "bonds_rmsd"/"angles_rmsd"
  # Use explicit None checks — 0.0 is a valid value.
  b = None
  a = None
  if log_metrics:
    b = log_metrics.get("bonds_rmsd")
    if b is None:
      b = log_metrics.get("bonds")
    a = log_metrics.get("angles_rmsd")
    if a is None:
      a = log_metrics.get("angles")
  if b is None:
    b = geom.get("bonds_rmsd")
  if a is None:
    a = geom.get("angles_rmsd")
  if b is not None and a is not None:
    lines.append("Bonds=%.4f Angles=%.2f" % (b, a))

  # Ramachandran
  rf = geom.get("rama_favored")
  ro = geom.get("rama_outliers")
  if rf is not None:
    line = "Rama: %.1f%% fav" % (rf * 100)
    if ro is not None:
      line += ", %.1f%% outlier" % (ro * 100)
    ol = geom.get("rama_outlier_list", [])
    if ol:
      shown = ol[:5]
      line += " (%s)" % ", ".join(shown)
      if len(ol) > 5:
        line += " +%d more" % (len(ol) - 5)
    lines.append(line)

  # Rotamer outliers
  rot = geom.get("rotamer_outliers")
  if rot is not None:
    line = "Rotamer outliers: %.1f%%" % (rot * 100)
    ol = geom.get("rotamer_outlier_list", [])
    if ol:
      shown = ol[:5]
      line += " (%s)" % ", ".join(shown)
      if len(ol) > 5:
        line += " +%d more" % (len(ol) - 5)
    lines.append(line)

  # Clashscore
  cs = geom.get("clashscore")
  if cs is not None:
    lines.append("Clashscore: %.1f" % cs)

  return lines


def _format_model_contents(contents):
  """Format model contents inventory."""
  lines = []

  # Ligands
  ligs = contents.get("ligands", [])
  if ligs:
    parts = []
    for lig in ligs[:5]:
      parts.append("%s (%s/%s)" % (
        lig["name"], lig["chain"], lig["resid"],
      ))
    line = "Ligands: " + ", ".join(parts)
    if len(ligs) > 5:
      line += " +%d more" % (len(ligs) - 5)
    lines.append(line)
  else:
    lines.append("Ligands: none")

  # Waters
  wc = contents.get("waters", 0)
  if wc > 0:
    lines.append("Waters: %d" % wc)

  # Ions
  ions = contents.get("ions", [])
  if ions:
    names = [i["name"] for i in ions[:5]]
    lines.append("Ions: " + ", ".join(names))

  return lines


def _format_ligand_rscc(rscc_list):
  """Format per-ligand RSCC."""
  lines = []
  for lig in rscc_list[:5]:
    lines.append(
      "  %s (%s/%s): RSCC=%.2f, B=%.0f" % (
        lig["name"], lig["chain"],
        lig["resid"], lig["rscc"],
        lig.get("b_mean", 0),
      )
    )
  return lines


def _format_diff_peaks(peaks):
  """Format difference density peaks."""
  lines = []
  pos = peaks.get("positive", [])
  neg = peaks.get("negative", [])
  if pos:
    parts = []
    for p in pos[:5]:
      parts.append(
        "%.1f-sigma near %s" % (
          p["height"],
          p.get("near_residue", "?"),
        )
      )
    lines.append(
      "Diff peaks (+): " + "; ".join(parts)
    )
  if neg:
    parts = []
    for p in neg[:3]:
      parts.append(
        "%.1f-sigma on %s" % (
          abs(p["height"]),
          p.get("near_residue", "?"),
        )
      )
    lines.append(
      "Diff peaks (-): " + "; ".join(parts)
    )
  return lines
