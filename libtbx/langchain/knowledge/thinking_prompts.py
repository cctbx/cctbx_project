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

CRITICAL: Check file provenance across cycles. The RECENT HISTORY \
section shows which files were used as input and produced as output \
for each step. Verify that each step uses the correct output from \
the previous step. Common errors include:
- Refinement using an old model instead of the latest combined \
  model (e.g. after pdbtools merges protein + ligand)
- Programs receiving data files instead of map coefficient files
- Output files from one step being ignored in favor of older files

IMPORTANT: Present evidence and reasoning. Do NOT issue commands \
or parameter settings directly. The planning agent will make the \
final decision based on your analysis.

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
  parts.append("CYCLE %s | %s | workflow: %s" % (cycle, exp, wf))

  # Program and metrics
  program = context.get("program_name", "")
  if program:
    parts.append("\nLast program: %s" % program)
  metrics = context.get("metrics", {})
  if metrics:
    metrics_str = _format_metrics(metrics)
    if metrics_str:
      parts.append("Current metrics: %s" % metrics_str)

  # Strategy memory summary
  if strategy_memory_dict and isinstance(strategy_memory_dict, dict):
    mem_summary = _format_memory(strategy_memory_dict)
    if mem_summary:
      parts.append("\nStrategy memory: %s" % mem_summary)

  # History summary (includes file provenance per cycle)
  history_summary = context.get("history_summary", "")
  if history_summary:
    parts.append("\n--- RECENT HISTORY (files used/produced) ---")
    parts.append(history_summary[:1200])

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
  parts.append("\nProvide your expert assessment as JSON.")

  user_msg = "\n".join(parts)
  return (SYSTEM_PROMPT, user_msg)


def _format_metrics(metrics):
  """Format metrics dict into a brief string."""
  if not metrics:
    return ""
  items = []
  # Prioritize key metrics
  for key in ["r_free", "r_work", "map_cc", "TFZ", "LLG",
               "FOM", "completeness", "resolution"]:
    val = metrics.get(key)
    if val is not None:
      items.append("%s=%s" % (key, val))
  # Add any remaining
  for key, val in sorted(metrics.items()):
    if key not in ["r_free", "r_work", "map_cc", "TFZ", "LLG",
                    "FOM", "completeness", "resolution"]:
      if val is not None and len(items) < 8:
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
        return result
    except (json.JSONDecodeError, ValueError):
      pass

  # Couldn't parse JSON — use the raw text as analysis
  defaults["analysis"] = text[:500]
  return defaults


# =========================================================================
# Stop-review prompt -- expert reviews a STOP decision
# =========================================================================

STOP_REVIEW_SYSTEM = """\
You are an expert crystallographer reviewing an automated \
structure determination workflow's decision to STOP.

The planning agent has decided to stop the workflow. \
Your job is to evaluate whether stopping is the right call, \
or whether there is a viable path forward that was missed.

Consider:
- Is the stated reason valid and insurmountable?
- Are there workarounds or alternative approaches?
- Could different parameters, programs, or strategies help?
- Is the data quality sufficient to continue?

Respond with a JSON object (no markdown fences):

{
  "agree_with_stop": true or false,
  "analysis": "2-3 sentence assessment of the stop decision",
  "alternatives": ["list of alternative approaches if any"],
  "guidance": "1-2 sentence recommendation for the user"
}
"""


def build_stop_review_prompt(
    stop_reasoning, history_summary, user_advice,
    available_files_summary, metrics_summary=""):
  """Build (system_msg, user_msg) for reviewing a STOP decision.

  Args:
    stop_reasoning: The planning agent's reason for stopping.
    history_summary: Brief summary of cycle history.
    user_advice: Original user instructions.
    available_files_summary: Brief list of available files.
    metrics_summary: Current metrics if any.

  Returns:
    (system_msg, user_msg) tuple.
  """
  parts = []
  parts.append("The planning agent has decided to STOP.")
  parts.append("")
  parts.append("STOP REASONING:")
  parts.append(str(stop_reasoning)[:600])

  if history_summary:
    parts.append("")
    parts.append("WORKFLOW HISTORY:")
    parts.append(str(history_summary)[:400])

  if metrics_summary:
    parts.append("")
    parts.append("CURRENT METRICS: %s" % str(metrics_summary)[:200])

  if user_advice:
    parts.append("")
    parts.append("USER INSTRUCTIONS:")
    parts.append(str(user_advice)[:300])

  if available_files_summary:
    parts.append("")
    parts.append("AVAILABLE FILES:")
    parts.append(str(available_files_summary)[:400])

  parts.append("")
  parts.append(
    "Should the workflow stop? Evaluate and respond as JSON.")

  return (STOP_REVIEW_SYSTEM, "\n".join(parts))


def parse_stop_review(llm_output):
  """Parse the stop-review LLM response.

  Returns dict with 'agree_with_stop', 'analysis',
  'alternatives', 'guidance'.
  """
  defaults = {
    "agree_with_stop": True,
    "analysis": "",
    "alternatives": [],
    "guidance": "",
  }

  if not llm_output or not isinstance(llm_output, str):
    return defaults

  # Strip markdown fences
  text = llm_output.strip()
  if text.startswith("```"):
    lines = text.split("\n")
    if lines[0].startswith("```"):
      lines = lines[1:]
    if lines and lines[-1].strip() == "```":
      lines = lines[:-1]
    text = "\n".join(lines)

  start = text.find("{")
  end = text.rfind("}")
  if start >= 0 and end > start:
    try:
      parsed = json.loads(text[start:end + 1])
      if isinstance(parsed, dict):
        result = dict(defaults)
        result.update(parsed)
        # Ensure boolean
        result["agree_with_stop"] = bool(
          result.get("agree_with_stop", True))
        return result
    except (json.JSONDecodeError, ValueError):
      pass

  defaults["analysis"] = text[:500]
  return defaults


# =========================================================================
# Pre-flight readiness check -- before cycle 1
# =========================================================================

PREFLIGHT_SYSTEM = """\
You are an expert crystallographer reviewing the inputs for an \
automated structure determination workflow BEFORE it begins.

Your job is to check whether all the necessary files and \
information are present to accomplish the user's stated goal. \
Catch problems early so the workflow does not waste cycles \
running programs that will eventually fail due to missing inputs.

Common requirements by goal:
- Molecular replacement: needs reflection data (.mtz), search model (.pdb)
- SAD/MAD phasing: needs reflection data (.mtz), sequence (.fa/.seq), \
possibly heavy atom sites
- Ligand fitting: needs reflection data (.mtz), model (.pdb), AND a \
ligand definition file (.cif or ligand .pdb with HETATM records)
- Refinement: needs reflection data (.mtz) and model (.pdb)
- Cryo-EM model building: needs map (.mrc/.ccp4), possibly model (.pdb)

Respond with a JSON object (no markdown fences):

{
  "ready": true or false,
  "analysis": "Brief assessment of readiness",
  "missing_files": ["list of missing file types, if any"],
  "warnings": ["list of potential issues"],
  "recommendation": "What the user should do if not ready"
}
"""


def build_preflight_prompt(
    user_advice, available_files, experiment_type=""):
  """Build (system_msg, user_msg) for the pre-flight readiness check.

  Args:
    user_advice: User instructions/goal.
    available_files: List of file paths.
    experiment_type: Known experiment type if any.

  Returns:
    (system_msg, user_msg) tuple.
  """
  import os

  parts = []
  parts.append("PRE-FLIGHT READINESS CHECK")
  parts.append("")

  if user_advice:
    parts.append("USER GOAL:")
    parts.append(str(user_advice)[:500])
  else:
    parts.append("USER GOAL: (not specified)")

  if experiment_type:
    parts.append("")
    parts.append("EXPERIMENT TYPE: %s" % experiment_type)

  parts.append("")
  parts.append("AVAILABLE FILES:")
  if available_files:
    for f in available_files[:30]:
      basename = os.path.basename(f)
      # Include extension info
      ext = os.path.splitext(basename)[1].lower()
      parts.append("  - %s  [%s]" % (basename, ext or "no extension"))
    if len(available_files) > 30:
      parts.append("  ... and %d more files" % (
        len(available_files) - 30))
  else:
    parts.append("  (no files provided)")

  parts.append("")
  parts.append(
    "Does the workflow have everything needed to accomplish "
    "the user's goal? Check for missing files and respond "
    "as JSON.")

  return (PREFLIGHT_SYSTEM, "\n".join(parts))


def parse_preflight(llm_output):
  """Parse the pre-flight check LLM response.

  Returns dict with 'ready', 'analysis', 'missing_files',
  'warnings', 'recommendation'.
  """
  defaults = {
    "ready": True,
    "analysis": "",
    "missing_files": [],
    "warnings": [],
    "recommendation": "",
  }

  if not llm_output or not isinstance(llm_output, str):
    return defaults

  # Strip markdown fences
  text = llm_output.strip()
  if text.startswith("```"):
    lines = text.split("\n")
    if lines[0].startswith("```"):
      lines = lines[1:]
    if lines and lines[-1].strip() == "```":
      lines = lines[:-1]
    text = "\n".join(lines)

  start = text.find("{")
  end = text.rfind("}")
  if start >= 0 and end > start:
    try:
      parsed = json.loads(text[start:end + 1])
      if isinstance(parsed, dict):
        result = dict(defaults)
        result.update(parsed)
        result["ready"] = bool(result.get("ready", True))
        return result
    except (json.JSONDecodeError, ValueError):
      pass

  defaults["analysis"] = text[:500]
  return defaults
