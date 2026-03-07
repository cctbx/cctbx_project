"""
HTML Report Generator for PHENIX AI Agent (v114).

Template-based HTML report using DisplayDataModel.
No LLM call — all data comes from StructureModel,
StructurePlan, and cycle history.

Includes inline SVG R-free/CC trajectory chart
with retreat markers.

Public API:
  generate_html_report(session_data) -> str
"""

from __future__ import (
  absolute_import, division, print_function,
)

import logging

logger = logging.getLogger(__name__)

# Import DisplayDataModel
try:
  from libtbx.langchain.agent.display_data_model \
    import DisplayDataModel
except ImportError:
  try:
    from agent.display_data_model \
      import DisplayDataModel
  except ImportError:
    DisplayDataModel = None


def generate_html_report(session_data):
  """Generate a complete HTML report from session.

  Args:
    session_data: dict from session.to_dict().

  Returns:
    str: Complete HTML document, or empty string
    if DisplayDataModel is unavailable.

  Never raises.
  """
  if DisplayDataModel is None:
    return ""
  try:
    return _generate_html_report_inner(
      session_data)
  except Exception:
    logger.debug(
      "generate_html_report failed",
      exc_info=True,
    )
    return ""


def _generate_html_report_inner(session_data):
  """Inner HTML report. May raise."""
  ddm = DisplayDataModel.from_session(
    session_data or {})

  parts = []
  parts.append(_HTML_HEAD)

  # Title
  status = ddm.outcome_status
  if status == "determined":
    title_icon = "&#x2705;"
    title_text = "Structure Determined"
  elif status == "stopped":
    title_icon = "&#x26A0;"
    title_text = "Agent Stopped"
  else:
    title_icon = "&#x25CF;"
    title_text = "Run Incomplete"

  parts.append(
    '<div class="header">'
    '<h1>%s %s</h1>'
    '<p class="subtitle">'
    'PHENIX AI Agent &mdash; '
    'Structure Determination Report</p>'
    '</div>' % (title_icon, _esc(title_text))
  )

  # Outcome message
  msg = ddm.outcome_message
  if msg:
    parts.append(
      '<p class="outcome-msg">%s</p>'
      % _esc(msg))

  # Stop reason
  if status == "stopped":
    reason = ddm.stop_reason
    if reason:
      parts.append(
        '<p class="stop-reason">%s</p>'
        % _esc(reason[:500]))

  # === DATA section ===
  data = ddm.data_summary
  if data:
    parts.append('<h2>Data</h2>')
    parts.append('<table class="metrics">')
    for label, value in data.items():
      parts.append(
        '<tr><td class="label">%s</td>'
        '<td>%s</td></tr>'
        % (_esc(label), _esc(value)))
    parts.append('</table>')

  # === FINAL MODEL section ===
  metrics = ddm.final_metrics
  if metrics:
    parts.append('<h2>Final Model</h2>')
    parts.append('<table class="metrics">')
    for label, value in metrics.items():
      parts.append(
        '<tr><td class="label">%s</td>'
        '<td>%s</td></tr>'
        % (_esc(label), _esc(value)))
    parts.append('</table>')

  # === TRAJECTORY CHART ===
  traj = ddm.primary_trajectory
  if len(traj) >= 2:
    is_cc = (
      ddm._experiment_type == "cryoem"
      and ddm.cc_trajectory)
    metric_name = "Map CC" if is_cc else "R-free"
    svg = _build_trajectory_svg(
      traj, metric_name)
    if svg:
      parts.append(
        '<h2>%s Trajectory</h2>' % metric_name)
      parts.append(
        '<div class="chart">%s</div>' % svg)

  # === WORKFLOW section ===
  stages = ddm.stage_outcomes
  timeline = ddm.timeline
  if stages or timeline:
    parts.append('<h2>Workflow</h2>')

  # Run stats
  stats = ddm.run_stats
  stats_parts = []
  if stats["n_cycles"]:
    stats_parts.append(
      "%d cycles" % stats["n_cycles"])
  if stats["n_programs"]:
    stats_parts.append(
      "%d programs" % stats["n_programs"])
  if stats["duration"]:
    stats_parts.append(stats["duration"])
  if stats_parts:
    parts.append(
      '<p class="run-stats">%s</p>'
      % _esc(", ".join(stats_parts)))

  # Stage outcomes
  if stages:
    parts.append(
      '<table class="stages">')
    for p in stages:
      if p.status == "complete":
        icon = "&#x2713;"
        cls = "stage-ok"
      elif p.status == "active":
        icon = "&#x25CF;"
        cls = "stage-active"
      elif p.status == "skipped":
        icon = "&#x2014;"
        cls = "stage-skip"
      elif p.status == "failed":
        icon = "&#x2717;"
        cls = "stage-fail"
      else:
        icon = "&#x25CB;"
        cls = "stage-pending"
      cycle_str = ""
      if p.cycles_used > 0:
        cycle_str = " (%d cycle%s)" % (
          p.cycles_used,
          "s" if p.cycles_used != 1 else "")
      parts.append(
        '<tr class="%s">'
        '<td class="stage-icon">%s</td>'
        '<td>%s%s</td></tr>'
        % (cls, icon,
           _esc(p.description),
           _esc(cycle_str)))
    parts.append('</table>')

  # Timeline table
  if timeline:
    parts.append(
      '<table class="timeline">'
      '<tr><th>Cycle</th><th>Program</th>'
      '<th>Status</th><th>Key Metric</th></tr>')
    for entry in timeline:
      if entry.status == "OK":
        s_icon = "&#x2705;"
      elif entry.status == "FAILED":
        s_icon = "&#x274C;"
      elif entry.status == "STOP":
        s_icon = "&#x23F9;"
      else:
        s_icon = "&#x2022;"
      metric_str = ""
      if (entry.key_metric_name
          and entry.key_metric_value is not None):
        metric_str = "%s: %s" % (
          entry.key_metric_name,
          _fmt_val(entry.key_metric_name,
                   entry.key_metric_value))
      parts.append(
        '<tr><td>%d</td><td>%s</td>'
        '<td>%s</td><td>%s</td></tr>'
        % (entry.cycle,
           _esc(entry.program),
           s_icon,
           _esc(metric_str)))
    parts.append('</table>')

  # === OUTPUT FILES ===
  files = ddm.output_files
  model_path = files.get("model_path")
  map_path = files.get("map_path")
  if model_path or map_path:
    parts.append('<h2>Output Files</h2>')
    parts.append('<table class="metrics">')
    if model_path:
      parts.append(
        '<tr><td class="label">Model</td>'
        '<td class="mono">%s</td></tr>'
        % _esc(model_path))
    if map_path:
      parts.append(
        '<tr><td class="label">Map</td>'
        '<td class="mono">%s</td></tr>'
        % _esc(map_path))
    parts.append('</table>')

  parts.append(_HTML_FOOT)
  return "\n".join(parts)


# ── SVG trajectory chart ──────────────────────

def _build_trajectory_svg(points, metric_name):
  """Build an inline SVG trajectory chart.

  Args:
    points: list of RfreePoint.
    metric_name: "R-free" or "Map CC".

  Returns:
    str: SVG markup, or empty string.
  """
  if len(points) < 2:
    return ""

  values = [p.value for p in points]
  cycles = [p.cycle for p in points]

  # Chart dimensions
  w, h = 500, 200
  pad_l, pad_r = 60, 20
  pad_t, pad_b = 20, 40
  plot_w = w - pad_l - pad_r
  plot_h = h - pad_t - pad_b

  # Value range (with 10% margin)
  v_min = min(values)
  v_max = max(values)
  v_range = v_max - v_min
  if v_range < 0.01:
    v_range = 0.1
    v_min -= 0.05
    v_max += 0.05
  else:
    margin = v_range * 0.1
    v_min -= margin
    v_max += margin
    v_range = v_max - v_min

  c_min = min(cycles)
  c_max = max(cycles)
  c_range = max(c_max - c_min, 1)

  def x_pos(cycle):
    return pad_l + (
      (cycle - c_min) / c_range * plot_w)

  def y_pos(val):
    return pad_t + plot_h - (
      (val - v_min) / v_range * plot_h)

  svg = []
  svg.append(
    '<svg xmlns="http://www.w3.org/2000/svg" '
    'viewBox="0 0 %d %d" '
    'style="max-width:%dpx;width:100%%">'
    % (w, h, w))

  # Background
  svg.append(
    '<rect x="%d" y="%d" width="%d" height="%d" '
    'fill="#fafafa" stroke="#ddd"/>'
    % (pad_l, pad_t, plot_w, plot_h))

  # Y-axis labels (3-4 ticks)
  n_ticks = 4
  for i in range(n_ticks + 1):
    val = v_min + (v_range * i / n_ticks)
    yp = y_pos(val)
    # Grid line
    svg.append(
      '<line x1="%d" y1="%.1f" '
      'x2="%d" y2="%.1f" '
      'stroke="#eee" stroke-width="1"/>'
      % (pad_l, yp, pad_l + plot_w, yp))
    # Label
    svg.append(
      '<text x="%d" y="%.1f" '
      'text-anchor="end" '
      'font-size="11" fill="#666">'
      '%.3f</text>'
      % (pad_l - 5, yp + 4, val))

  # X-axis labels
  for c in cycles:
    xp = x_pos(c)
    svg.append(
      '<text x="%.1f" y="%d" '
      'text-anchor="middle" '
      'font-size="11" fill="#666">%d</text>'
      % (xp, h - pad_b + 20, c))

  # Axis labels
  svg.append(
    '<text x="%d" y="%d" '
    'text-anchor="middle" '
    'font-size="12" fill="#333">Cycle</text>'
    % (pad_l + plot_w // 2, h - 2))
  svg.append(
    '<text x="12" y="%d" '
    'text-anchor="middle" '
    'font-size="12" fill="#333" '
    'transform="rotate(-90, 12, %d)">'
    '%s</text>'
    % (pad_t + plot_h // 2,
       pad_t + plot_h // 2,
       _esc(metric_name)))

  # Retreat markers (vertical dashed lines)
  for pt in points:
    if pt.is_retreat_point:
      xp = x_pos(pt.cycle)
      svg.append(
        '<line x1="%.1f" y1="%d" '
        'x2="%.1f" y2="%d" '
        'stroke="#e74c3c" '
        'stroke-width="1.5" '
        'stroke-dasharray="4,3"/>'
        % (xp, pad_t, xp, pad_t + plot_h))
      svg.append(
        '<text x="%.1f" y="%d" '
        'text-anchor="middle" '
        'font-size="9" fill="#e74c3c">'
        '&#x26A0; retreat</text>'
        % (xp, pad_t - 4))

  # Data line segments (break at retreats)
  segments = []
  current_segment = []
  for pt in points:
    if pt.is_retreat_point and current_segment:
      segments.append(current_segment)
      current_segment = []
    current_segment.append(pt)
  if current_segment:
    segments.append(current_segment)

  colors = ["#2980b9", "#27ae60", "#8e44ad",
            "#e67e22"]
  for seg_idx, segment in enumerate(segments):
    if len(segment) < 2:
      continue
    color = colors[seg_idx % len(colors)]
    path_pts = " ".join(
      "%.1f,%.1f" % (
        x_pos(pt.cycle), y_pos(pt.value))
      for pt in segment
    )
    svg.append(
      '<polyline points="%s" '
      'fill="none" stroke="%s" '
      'stroke-width="2"/>'
      % (path_pts, color))

  # Data points (circles)
  for pt in points:
    xp = x_pos(pt.cycle)
    yp = y_pos(pt.value)
    color = "#e74c3c" if pt.is_retreat_point \
      else "#2980b9"
    svg.append(
      '<circle cx="%.1f" cy="%.1f" r="4" '
      'fill="%s" stroke="white" '
      'stroke-width="1.5"/>'
      % (xp, yp, color))

  svg.append('</svg>')
  return "\n".join(svg)


# ── Helpers ───────────────────────────────────

def _esc(text):
  """HTML-escape a string."""
  if not text:
    return ""
  return (str(text)
          .replace("&", "&amp;")
          .replace("<", "&lt;")
          .replace(">", "&gt;")
          .replace('"', "&quot;"))


def _fmt_val(name, value):
  """Format a metric value for the timeline table."""
  _3d = {"R-free", "R-work", "Map CC", "FOM",
         "Ligand CC"}
  _2d = {"Res"}
  try:
    v = float(value)
    if name in _3d:
      return "%.3f" % v
    elif name in _2d:
      return "%.2f" % v
    else:
      return "%.1f" % v
  except (ValueError, TypeError):
    return str(value)


# ── HTML template fragments ───────────────────

_HTML_HEAD = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PHENIX AI Agent Report</title>
<style>
  body {
    font-family: -apple-system, BlinkMacSystemFont,
      "Segoe UI", Roboto, Helvetica, sans-serif;
    max-width: 700px;
    margin: 40px auto;
    padding: 0 20px;
    color: #333;
    line-height: 1.5;
  }
  .header { margin-bottom: 24px; }
  h1 {
    font-size: 24px;
    margin: 0 0 4px 0;
    color: #222;
  }
  .subtitle {
    color: #888;
    font-size: 14px;
    margin: 0;
  }
  .outcome-msg {
    font-size: 16px;
    color: #555;
    margin: 8px 0 16px 0;
  }
  .stop-reason {
    background: #fff3cd;
    border-left: 4px solid #ffc107;
    padding: 8px 12px;
    margin: 8px 0 16px 0;
    font-size: 14px;
    color: #856404;
  }
  h2 {
    font-size: 18px;
    color: #444;
    border-bottom: 1px solid #eee;
    padding-bottom: 4px;
    margin: 28px 0 12px 0;
  }
  table.metrics {
    border-collapse: collapse;
    margin: 8px 0;
  }
  table.metrics td {
    padding: 3px 16px 3px 0;
    font-size: 14px;
  }
  table.metrics td.label {
    font-weight: 600;
    color: #555;
    white-space: nowrap;
  }
  td.mono {
    font-family: "SF Mono", Monaco, Consolas,
      monospace;
    font-size: 13px;
    word-break: break-all;
  }
  table.stages {
    border-collapse: collapse;
    margin: 8px 0 16px 0;
  }
  table.stages td {
    padding: 4px 8px;
    font-size: 14px;
  }
  td.stage-icon { width: 24px; text-align: center; }
  .stage-ok td { color: #333; }
  .stage-active td { color: #2980b9; font-weight: 600; }
  .stage-fail td { color: #c0392b; }
  .stage-skip td { color: #999; }
  .stage-pending td { color: #bbb; }
  table.timeline {
    border-collapse: collapse;
    margin: 12px 0;
    width: 100%;
  }
  table.timeline th {
    text-align: left;
    font-size: 12px;
    color: #888;
    border-bottom: 1px solid #ddd;
    padding: 4px 12px 4px 0;
  }
  table.timeline td {
    padding: 4px 12px 4px 0;
    font-size: 14px;
    border-bottom: 1px solid #f0f0f0;
  }
  .chart {
    margin: 12px 0;
    text-align: center;
  }
  .run-stats {
    font-size: 14px;
    color: #666;
    margin: 4px 0 12px 0;
  }
  @media print {
    body { margin: 20px; }
  }
</style>
</head>
<body>
"""

_HTML_FOOT = """\
<div style="margin-top:40px;padding-top:12px;
  border-top:1px solid #eee;color:#aaa;
  font-size:12px;">
  Generated by PHENIX AI Agent
</div>
</body>
</html>
"""
