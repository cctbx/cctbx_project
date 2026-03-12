"""
Thinking Agent Prompt Builder (v113).

Builds (system_msg, user_msg) for the expert crystallographer
reasoning LLM call. Imports nothing from agent — pure data
formatting so it can be iterated on without touching any
workflow code.
"""

from __future__ import absolute_import, division, print_function
import json


# =========================================================================
# System prompt — defines the expert persona and output format
# =========================================================================

SYSTEM_PROMPT = """\
You are an expert crystallographer reviewing an automated \
structure determination workflow. You have decades of experience \
with X-ray crystallography, cryo-EM, and the PHENIX software suite.

Your role is to analyze the latest program output and provide \
strategic guidance. You are NOT running programs yourself — you \
are advising the automated agent on what to do next.

IMPORTANT: Present evidence and reasoning. Do NOT issue commands \
or parameter settings directly. The planning agent will make the \
final decision based on your analysis.

IMPORTANT: The "Available programs" line shows which programs \
the workflow engine will allow on the next cycle. Your guidance \
MUST recommend one of these programs. Do not suggest a program \
that is not in the available list — it will confuse the user \
when the agent runs a different program than you recommended.

If a "Plan goal" or "Current plan prefers" line is shown, \
align your guidance with the plan's strategy. If a "User stop \
condition" is shown, respect it in your recommendation.

FILE PROVENANCE: The RECENT HISTORY section shows input and \
output files for each cycle. Always verify that each program \
received the correct input file — a common failure mode is \
refinement running on a stale or wrong model. Flag any \
suspicious file provenance in your concerns list.

Respond with a JSON object (no markdown fences):

{
  "analysis": "2-3 sentence scientific assessment of the current state",
  "confidence": "high|medium|low|hopeless",
  "action": "guide_step|let_run|stop|pivot",
  "guidance": "1-2 sentence advice for the planning agent",
  "data_quality": "brief quality note or empty string",
  "phasing_strategy": "current strategy or empty string",
  "concerns": ["list of concerns, if any"],
  "alternatives": ["list of alternative approaches, if any"]
}

Actions:
- guide_step: Provide specific guidance for the next step
- let_run: Routine step, let the agent proceed normally
- stop: Recommend stopping (convergence reached or hopeless)
- pivot: Recommend changing strategy

Confidence levels:
- high: Clear path forward, good data quality
- medium: Some uncertainty, proceed with caution
- low: Significant problems, may need strategy change
- hopeless: Fundamental issues, consider stopping
"""


# =========================================================================
# User message builder
# =========================================================================

def build_thinking_prompt(context, strategy_memory_dict=None):
  """Build the (system_msg, user_msg) tuple for the thinking LLM.

  Args:
    context: Dict with keys:
      - log_sections: str (from log_section_extractor)
      - program_name: str (last program run)
      - cycle_number: int
      - experiment_type: str ("xray" or "cryoem")
      - workflow_state: str (e.g. "initial", "refinement")
      - metrics: dict (latest metrics)
      - user_advice: str (original user advice)
      - history_summary: str (brief history)
      - validation_report: str (Phase A, optional)
      - r_free_trend: list of float (Phase A, optional)
      - structure_model_summary: str (v114, optional)
      - current_problems: list of dict (v114, optional)
      - hypothesis_prompt: str (Phase 4, optional)
    strategy_memory_dict: Dict from StrategyMemory.to_dict()
      or None.

  Returns:
    (system_msg, user_msg) tuple of strings.
  """
  parts = []

  # Header
  cycle = context.get("cycle_number", "?")
  exp = context.get("experiment_type", "unknown")
  wf = context.get("workflow_state", "unknown")
  parts.append(
    "CYCLE %s | %s | workflow: %s" % (cycle, exp, wf)
  )

  # Valid programs (so guidance aligns with what the
  # agent can actually do next)
  valid_progs = context.get("valid_programs", [])
  non_stop = [p for p in valid_progs if p != "STOP"]
  if non_stop:
    parts.append(
      "Available programs: %s"
      % ", ".join(non_stop)
    )

  # Plan phase context (so guidance aligns with
  # the agent's current strategy)
  plan_phase = context.get("plan_phase", "")
  if plan_phase:
    parts.append(plan_phase)
  plan_goal = context.get("plan_goal", "")
  if plan_goal:
    parts.append("Plan goal: %s" % plan_goal)
  asu_copies = context.get("asu_copies")
  if asu_copies:
    parts.append(
      "ASU copy count: %d (used by Phaser as search_copies)" % asu_copies)
  stop_after = context.get("stop_after", "")
  if stop_after:
    parts.append(
      "User stop condition: stop after %s"
      % stop_after
    )

  # File inventory (v115.05) — show what files are
  # available by category so the expert doesn't
  # hallucinate missing files.
  file_inventory = context.get("file_inventory", "")
  if file_inventory:
    parts.append(
      "\n=== AVAILABLE FILES ===\n%s"
      % file_inventory
    )

  # Program and metrics
  program = context.get("program_name", "")
  if program:
    parts.append("\nLast program: %s" % program)
  metrics = context.get("metrics", {})
  if metrics:
    metrics_str = _format_metrics(metrics)
    if metrics_str:
      parts.append("Current metrics: %s" % metrics_str)

  # Structural validation report (Phase A)
  validation_report = context.get(
    "validation_report", ""
  )
  if validation_report:
    parts.append("\n" + validation_report)

  # R-free trend (Phase A)
  r_free_trend = context.get("r_free_trend", [])
  if r_free_trend and len(r_free_trend) >= 2:
    try:
      from libtbx.langchain.agent.format_validation \
        import format_validation_trend
    except ImportError:
      try:
        from agent.format_validation import (
          format_validation_trend,
        )
      except ImportError:
        format_validation_trend = None
    if format_validation_trend is not None:
      trend_str = format_validation_trend(r_free_trend)
      if trend_str:
        parts.append(trend_str)

  # Structure Model summary (v114)
  sm_summary = context.get(
    "structure_model_summary", ""
  )
  if sm_summary:
    parts.append(
      "\n=== STRUCTURE MODEL ===\n%s" % sm_summary
    )

  # Current problems from Structure Model (v114)
  current_problems = context.get(
    "current_problems", []
  )
  if current_problems:
    prob_lines = []
    for p in current_problems[:5]:
      if isinstance(p, dict):
        prob_lines.append(
          "- %s" % p.get("problem", "")
        )
    if prob_lines:
      parts.append(
        "Problems:\n%s" % "\n".join(prob_lines)
      )

  # Error classification from previous cycle (Fix 3, v115)
  # When the previous cycle failed, provide structured error
  # info so the expert can give recovery guidance.
  error_class = context.get("error_classification")
  if error_class and isinstance(error_class, dict):
    _ecat = error_class.get("category", "")
    if _ecat and _ecat != "NO_ERROR":
      _emsg = error_class.get("error_message", "")
      _esug = error_class.get("suggestion", "")
      _efc = context.get("failure_count", 0)
      _elines = [
        "=== PREVIOUS FAILURE ===",
        "Error: %s — %s" % (_ecat, _emsg),
      ]
      if _esug:
        _elines.append("Suggested fix: %s" % _esug)
      bad = error_class.get("bad_params", [])
      if bad:
        _elines.append(
          "Bad params: %s" % ", ".join(bad))
      if _efc >= 2:
        _elines.append(
          "CONSECUTIVE FAILURES: %d — recommend "
          "pivoting to a different program" % _efc)
      parts.append("\n".join(_elines))

  # File metadata summary (Phase C)
  file_metadata = context.get("file_metadata", {})
  if file_metadata:
    try:
      try:
        from libtbx.langchain.agent.file_metadata \
          import format_file_metadata_block
      except ImportError:
        from agent.file_metadata \
          import format_file_metadata_block
      fm_block = format_file_metadata_block(
        file_metadata
      )
      if fm_block:
        parts.append("\n" + fm_block)
    except ImportError:
      pass

  # Expert knowledge base rules (Phase D)
  kb_rules_text = context.get("kb_rules_text", "")
  if kb_rules_text:
    parts.append("\n" + kb_rules_text)

  # Hypothesis prompt (Phase 4)
  hypothesis_prompt = context.get(
    "hypothesis_prompt", ""
  )
  if hypothesis_prompt:
    parts.append(hypothesis_prompt)

  # Strategy memory summary
  if strategy_memory_dict \
     and isinstance(strategy_memory_dict, dict):
    mem_summary = _format_memory(strategy_memory_dict)
    if mem_summary:
      parts.append(
        "\nStrategy memory: %s" % mem_summary
      )

  # History summary
  history_summary = context.get("history_summary", "")
  if history_summary:
    parts.append(
      "\nHistory: %s" % history_summary[:300]
    )

  # Recent failures — explicit so expert knows what
  # already failed and can recommend alternatives
  recent_failures = context.get(
    "recent_failures", [])
  if recent_failures:
    parts.append(
      "\n--- RECENT FAILURES (try a different "
      "approach) ---")
    for f in recent_failures[-3:]:
      parts.append(
        "Cycle %s: %s FAILED: %s"
        % (f.get("cycle", "?"),
           f.get("program", "?"),
           f.get("error", "?")[:150]))
      cmd = f.get("command", "")
      if cmd:
        parts.append(
          "  Command was: %s" % cmd[:150])
    parts.append(
      "Do NOT recommend the same command. "
      "Suggest different parameters, a "
      "different strategy, or a different "
      "program.")

  # User advice (brief excerpt)
  user_advice = context.get("user_advice", "")
  if user_advice:
    advice_excerpt = user_advice[:200]
    if len(user_advice) > 200:
      advice_excerpt += "..."
    parts.append("\nUser advice: %s" % advice_excerpt)

  # Log sections (the bulk of the content)
  log_sections = context.get("log_sections", "")
  if log_sections:
    parts.append("\n--- PROGRAM OUTPUT ---")
    parts.append(log_sections)
  else:
    parts.append("\n(No program output available)")

  # Assemble
  parts.append(
    "\nProvide your expert assessment as JSON."
  )

  user_msg = "\n".join(parts)

  # P1B: inject stop_reason table into system prompt
  # when errors.yaml is available.
  _table = build_stop_reason_table()
  if _table:
    _system = (
      SYSTEM_PROMPT
      + "\n\nWhen recommending \"stop\", use the "
      "reference list below to classify the "
      "reason.  Put the matching code in the "
      "\"stop_reason_code\" field of your JSON "
      "response.\n\n"
      + _table
    )
  else:
    _system = SYSTEM_PROMPT
  return (_system, user_msg)


def _format_metrics(metrics):
  """Format metrics dict into a brief string.

  The metrics dict may come directly from log_analysis
  (flat dict with extra keys like 'program'). We filter
  to only numeric/meaningful metrics.
  """
  if not metrics:
    return ""

  # Keys to skip — not metrics, just metadata
  _SKIP = frozenset([
    "program", "error_type", "error_message",
    "warning", "output_files", "status",
    "macro_cycles", "detected_program",
  ])

  items = []
  # Prioritize key metrics
  _PRIORITY = [
    "r_free", "r_work", "map_cc", "TFZ", "LLG",
    "FOM", "completeness", "resolution",
    "bonds_rmsd", "angles_rmsd", "clashscore",
  ]
  for key in _PRIORITY:
    val = metrics.get(key)
    if val is not None:
      items.append("%s=%s" % (key, val))
  # Add any remaining numeric values
  priority_set = set(_PRIORITY)
  for key, val in sorted(metrics.items()):
    if key in priority_set or key in _SKIP:
      continue
    if val is not None and len(items) < 10:
      # Only include numeric values
      if isinstance(val, (int, float)):
        items.append("%s=%s" % (key, val))
  return ", ".join(items)


def _format_memory(memory_dict):
  """Format strategy memory dict into a brief string."""
  parts = []
  if memory_dict.get("data_quality"):
    parts.append("quality=%s" % memory_dict["data_quality"])
  if memory_dict.get("phasing_strategy"):
    parts.append("strategy=%s" % memory_dict["phasing_strategy"])
  rfree = memory_dict.get("r_free_history", [])
  if rfree:
    parts.append("r_free_trend=[%s]" % ", ".join(
      "%.3f" % v for v in rfree[-5:]))
  concerns = memory_dict.get("concerns", [])
  if concerns:
    parts.append("concerns: %s" % "; ".join(concerns[:3]))
  decisions = memory_dict.get("decisions", [])
  if decisions:
    last = decisions[-1]
    parts.append("last_decision: cycle %s — %s" % (last[0], last[1]))
  return "; ".join(parts) if parts else ""


def parse_assessment(llm_output):
  """Parse the LLM's JSON response into a structured assessment.

  Tolerant of markdown fences, trailing text, and missing keys.

  Args:
    llm_output: Raw string from the LLM.

  Returns:
    Dict with assessment fields. Always has at least
    'analysis', 'confidence', 'action', 'guidance'.
  """
  defaults = {
    "analysis": "",
    "confidence": "medium",
    "action": "let_run",
    "guidance": "",
    "data_quality": "",
    "phasing_strategy": "",
    "concerns": [],
    "alternatives": [],
    "stop_reason_code": None,
    "file_overrides": {},
  }

  if not llm_output or not isinstance(llm_output, str):
    return defaults

  # Strip markdown fences
  text = llm_output.strip()
  if text.startswith("```"):
    lines = text.split("\n")
    # Remove first and last fence lines
    if lines[0].startswith("```"):
      lines = lines[1:]
    if lines and lines[-1].strip() == "```":
      lines = lines[:-1]
    text = "\n".join(lines)

  # Try to find JSON object
  start = text.find("{")
  end = text.rfind("}")
  if start >= 0 and end > start:
    try:
      parsed = json.loads(text[start:end + 1])
      if isinstance(parsed, dict):
        result = dict(defaults)
        result.update(parsed)
        # Validate action
        if result["action"] not in (
            "guide_step", "let_run", "stop", "pivot"):
          result["action"] = "let_run"
        # Validate confidence
        if result["confidence"] not in (
            "high", "medium", "low", "hopeless"):
          result["confidence"] = "medium"
        # P1B: validate stop_reason_code
        _src = result.get("stop_reason_code")
        if _src is not None:
          try:
            from agent.graph_state import (
              STOP_REASON_CODES)
          except ImportError:
            try:
              from graph_state import (
                STOP_REASON_CODES)
            except ImportError:
              STOP_REASON_CODES = frozenset()
          if _src not in STOP_REASON_CODES:
            result["stop_reason_code"] = None
        else:
          result["stop_reason_code"] = None
        # P1B: normalize file_overrides
        _fo = result.get("file_overrides")
        if not isinstance(_fo, dict):
          result["file_overrides"] = {}
        return result
    except (json.JSONDecodeError, ValueError):
      pass

  # Couldn't parse JSON — use the raw text as analysis
  defaults["analysis"] = text[:500]
  return defaults


# =========================================================================
# P1B: Stop-reason table for system prompt
# =========================================================================

def _load_errors_yaml():
  """Load stop reason codes from knowledge/errors.yaml.

  Returns list of dicts with 'code', 'summary',
  'trigger_evidence', 'negative_constraint' keys.
  Returns [] if the file is missing or unparseable.

  Never raises.
  """
  import os
  try:
    import yaml
  except ImportError:
    return []
  _yaml_dir = os.path.join(
    os.path.dirname(__file__))
  _yaml_path = os.path.join(
    _yaml_dir, "errors.yaml")
  if not os.path.isfile(_yaml_path):
    return []
  try:
    with open(_yaml_path) as f:
      data = yaml.safe_load(f)
    if isinstance(data, dict):
      return data.get("stop_reason_codes", [])
    return []
  except Exception:
    return []


def build_stop_reason_table():
  """Build a text table of stop reason codes for
  injection into the system prompt.

  Returns a multi-line string, or "" if errors.yaml
  is absent or empty.

  Never raises.
  """
  entries = _load_errors_yaml()
  if not entries:
    return ""
  lines = [
    "STOP REASON CODES — use one of these codes in "
    "the stop_reason_code field when action is "
    "\"stop\":",
    "",
  ]
  for entry in entries:
    code = entry.get("code", "")
    summary = entry.get("summary", "")
    trigger = entry.get(
      "trigger_evidence", "").strip()
    negative = entry.get(
      "negative_constraint", "").strip()
    lines.append("  %s — %s" % (code, summary))
    if trigger:
      lines.append(
        "    USE WHEN: %s" % trigger[:200])
    if negative:
      lines.append(
        "    DO NOT USE WHEN: %s"
        % negative[:200])
    lines.append("")
  return "\n".join(lines)
