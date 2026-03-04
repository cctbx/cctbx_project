"""
Thinking Agent Module (v113).

Expert crystallographer reasoning node for the PHENIX AI Agent.
Provides strategic guidance to the plan node by analyzing program
logs with domain expertise.

Entry point: run_think_node(state) -- called from graph_nodes.think().
"""

from __future__ import absolute_import, division, print_function
import traceback

# --- Imports from B1, B2, B3 ---
# Use try/except for libtbx paths vs standalone testing
try:
  from libtbx.langchain.agent.strategy_memory import StrategyMemory
  from libtbx.langchain.agent.log_section_extractor import extract_sections
  from libtbx.langchain.knowledge.thinking_prompts import (
    build_thinking_prompt, parse_assessment,
    build_stop_review_prompt, parse_stop_review,
    build_preflight_prompt, parse_preflight)
except ImportError:
  from agent.strategy_memory import StrategyMemory
  from agent.log_section_extractor import extract_sections
  from knowledge.thinking_prompts import (
    build_thinking_prompt, parse_assessment,
    build_stop_review_prompt, parse_stop_review,
    build_preflight_prompt, parse_preflight)


# =========================================================================
# Strategic program list -- programs whose output warrants expert review
# =========================================================================

_STRATEGIC_PROGRAMS = frozenset([
  "phenix.mtriage", "mtriage",
  "phenix.xtriage", "xtriage",
  "phenix.phaser", "phaser",
  "phenix.autosol", "autosol",
  "phenix.autobuild", "autobuild",
  "phenix.model_vs_data", "model_vs_data",
  "phenix.real_space_refine", "real_space_refine",
  "phenix.map_to_model", "map_to_model",
  "phenix.predict_and_build", "predict_and_build",
])


def should_think(state):
  """Decide whether the thinking LLM should engage this cycle.

  The thinking LLM is expensive. Only invoke it at:
  - After strategic programs (xtriage, phaser, autosol, autobuild, etc.)
  - After failures
  - When R-free is stalled (3+ cycles without improvement)

  NOT invoked on:
  - First cycle (no program output to analyze yet)
  - Routine refinement steps

  Args:
    state: Graph state dict.

  Returns:
    True if the thinking LLM should be called.
  """
  history = state.get("history", [])

  # First cycle: no program output to analyze yet.
  # The first useful think is AFTER xtriage/first program runs.
  if not history:
    print("  [THINK] should_think: False (no history yet)")
    return False

  last = history[-1]
  program = last.get("program", "")

  # Strategic programs
  if any(p in program for p in _STRATEGIC_PROGRAMS):
    print("  [THINK] should_think: True (strategic program: %s)" % program)
    return True

  # Last cycle failed
  # NOTE: History entries from build_request_v2 use "result" key, not "status".
  # Check both for robustness.
  result_text = last.get("result", last.get("status", ""))
  if result_text and "FAIL" in str(result_text).upper():
    print("  [THINK] should_think: True (last cycle failed: %s)"
          % str(result_text)[:80])
    return True

  # R-free stalled
  memory_dict = state.get("strategy_memory", {})
  if memory_dict:
    memory = StrategyMemory.from_dict(memory_dict)
    if memory.metrics_stalled():
      print("  [THINK] should_think: True (metrics stalled)")
      return True

  print("  [THINK] should_think: False (routine step, program=%s)"
        % program)
  return False


def run_think_node(state):
  """Main entry point -- called from graph_nodes.think().

  Builds context, calls the thinking LLM, parses the response,
  and enriches the state for the plan node.

  Args:
    state: Graph state dict.

  Returns:
    Updated state dict.
  """
  # Safety: if not enabled, return unchanged
  if not state.get("use_thinking_agent"):
    return state

  print("  [THINK] Thinking agent is enabled for cycle %s"
        % state.get("cycle_number", "?"))

  # Check if we should think this cycle
  if not should_think(state):
    state = _log(state, "THINK: Skipping (routine step)")
    return state

  print("  [THINK] Expert reasoning engaged -- calling LLM")
  state = _log(state, "THINK: Expert reasoning engaged")

  try:
    # Build context
    context = _build_thinking_context(state)
    print("  [THINK] Context built: program=%s, experiment=%s, "
          "log_sections=%d chars"
          % (context.get("program_name", "?"),
             context.get("experiment_type", "?"),
             len(context.get("log_sections", ""))))

    # Build prompt
    memory_dict = state.get("strategy_memory", {})
    system_msg, user_msg = build_thinking_prompt(context, memory_dict)
    print("  [THINK] Prompt built: system=%d chars, user=%d chars"
          % (len(system_msg), len(user_msg)))

    # Call LLM
    raw_response = _call_thinking_llm(state, system_msg, user_msg)

    if raw_response is None:
      print("  [THINK] WARNING: LLM call returned None -- "
            "continuing without expert reasoning")
      state = _log(state, "THINK: LLM call failed, continuing without")
      return state

    print("  [THINK] LLM response received: %d chars" % len(raw_response))

    # Parse response
    assessment = parse_assessment(raw_response)
    print("  [THINK] Assessment: action=%s confidence=%s"
          % (assessment["action"], assessment["confidence"]))
    analysis = assessment.get("analysis", "")
    if analysis:
      # Print first 200 chars of analysis to log
      display = analysis.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Analysis: %s" % display)
    guidance = assessment.get("guidance", "")
    if guidance:
      display = guidance.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Guidance: %s" % display)

    state = _log(state,
      "THINK: Assessment -- action=%s confidence=%s"
      % (assessment["action"], assessment["confidence"]))

    # Emit structured event for GUI Events panel
    state = _emit(state,
      action=assessment["action"],
      confidence=assessment["confidence"],
      analysis=assessment.get("analysis", "")[:200])

    # Expert recommends stopping?
    if assessment["action"] == "stop":
      print("  [THINK] Expert recommends STOP")
      state = _log(state, "THINK: Expert recommends STOP")
      cycle = state.get("cycle_number", 1)
      abort_msg = "Expert assessment: %s" % assessment["analysis"][:300]
      return {
        **state,
        "command": "STOP",
        "stop": True,
        "stop_reason": "expert: %s" % assessment["analysis"][:200],
        "abort_message": abort_msg,
        "expert_assessment": assessment,
        "strategy_memory": _update_memory(
          memory_dict, assessment, cycle),
      }

    # Inject guidance into user_advice for PLAN and BUILD
    if guidance:
      enriched_advice = (
        "[Expert assessment] %s\n\n%s"
        % (guidance, state.get("user_advice", "")))
      print("  [THINK] Injected expert guidance into user_advice")
    else:
      enriched_advice = state.get("user_advice", "")

    cycle = state.get("cycle_number", 1)
    return {
      **state,
      "user_advice": enriched_advice,
      "expert_assessment": assessment,
      "strategy_memory": _update_memory(
        memory_dict, assessment, cycle),
    }

  except Exception as e:
    # Never crash the workflow -- degrade gracefully
    print("  [THINK] ERROR (continuing without): %s" % str(e))
    print("  [THINK] %s" % traceback.format_exc()[:500])
    state = _log(state,
      "THINK: Error (continuing without): %s" % str(e))
    state = _log(state, "THINK: %s" % traceback.format_exc()[:500])
    return state


# =========================================================================
# Internal helpers
# =========================================================================

def _build_thinking_context(state):
  """Assemble context dict for the thinking prompt.

  Gathers log sections, metrics, workflow state, and history
  into the format expected by build_thinking_prompt().
  """
  history = state.get("history", [])
  last = history[-1] if history else {}
  program = last.get("program", "")

  # Extract informative log sections
  log_text = state.get("log_text", "")
  log_sections = ""
  if log_text and program:
    log_sections = extract_sections(log_text, program)
  elif log_text:
    # No program name -- use generic tail extraction
    log_sections = extract_sections(log_text, "unknown")

  # Get experiment type and workflow state
  session_info = state.get("session_info", {})
  workflow_state_dict = state.get("workflow_state", {})

  # Extract recent metrics from log_analysis
  log_analysis = state.get("log_analysis", {})
  metrics = log_analysis.get("metrics", {})

  # Build brief history summary
  history_summary = _summarize_history(history)

  return {
    "log_sections": log_sections,
    "program_name": program,
    "cycle_number": state.get("cycle_number", 1),
    "experiment_type": session_info.get("experiment_type", "unknown"),
    "workflow_state": workflow_state_dict.get("state", "unknown"),
    "metrics": metrics,
    "user_advice": state.get("user_advice", ""),
    "history_summary": history_summary,
  }


def _summarize_history(history):
  """Create a detailed text summary of cycle history.

  Includes input files (from command), output files, and result
  so the expert can trace file provenance across cycles — e.g.
  detect when refinement uses the wrong model after pdbtools.
  """
  if not history:
    return "(first cycle)"

  import os
  lines = []
  for h in history[-6:]:  # Last 6 cycles for better context
    cycle_num = h.get("cycle_number", h.get("cycle", "?"))
    prog = h.get("program", "?")
    result = h.get("result", h.get("status", "?"))

    # Determine OK/FAILED
    result_str = str(result) if result else "?"
    if "FAILED" in result_str.upper():
      # Include error reason (compact)
      status = result_str[:80]
    elif "SUCCESS" in result_str.upper():
      status = "OK"
    else:
      status = result_str[:40]

    # Extract input file basenames from command
    command = h.get("command", "")
    input_files = _extract_file_basenames(command)

    # Extract output file basenames
    output_files = h.get("output_files", [])
    output_basenames = []
    if output_files:
      for f in output_files[:4]:
        bn = os.path.basename(f) if f else ""
        if bn and bn.endswith(('.pdb', '.mtz', '.cif', '.ccp4', '.mrc')):
          output_basenames.append(bn)

    # Build compact multi-line entry
    lines.append("  C%s: %s -> %s" % (cycle_num, prog, status))
    if input_files:
      lines.append("      Input: %s" % ", ".join(input_files[:5]))
    if output_basenames:
      lines.append("      Output: %s" % ", ".join(
        output_basenames[:4]))

  return "\n".join(lines)


def _extract_file_basenames(command):
  """Extract file basenames from a phenix command string.

  Parses basenames of .pdb, .mtz, .cif, .ccp4, .mrc, .fa, .seq
  files from the command, ignoring parameter assignments.
  """
  import os
  if not command:
    return []
  basenames = []
  file_extensions = (
    '.pdb', '.mtz', '.cif', '.ccp4', '.mrc', '.fa', '.seq',
    '.map')
  for token in command.split():
    # Skip parameter assignments (key=value where value isn't a file)
    if '=' in token:
      # But handle model=/path/to/file.pdb style
      _val = token.split('=', 1)[-1]
      if any(_val.lower().endswith(ext) for ext in file_extensions):
        basenames.append(os.path.basename(_val))
      continue
    if any(token.lower().endswith(ext) for ext in file_extensions):
      basenames.append(os.path.basename(token))
  return basenames


def _update_memory(memory_dict, assessment, cycle_number=None):
  """Update strategy memory from assessment and return as dict."""
  memory = StrategyMemory.from_dict(memory_dict)
  memory.update(assessment)

  # Record the decision
  if assessment.get("guidance"):
    cycle = cycle_number if cycle_number is not None else (
      len(memory.decisions) + 1)
    memory.record_decision(cycle, assessment["guidance"][:100])

  return memory.to_dict()


def _call_thinking_llm(state, system_msg, user_msg):
  """Call the LLM for expert reasoning.

  Uses the same LLM infrastructure as the plan node.

  Args:
    state: Graph state dict (for provider).
    system_msg: System prompt string.
    user_msg: User prompt string.

  Returns:
    Raw response string, or None on failure.
  """
  provider = state.get("provider", "google")
  print("  [THINK] Calling LLM (provider=%s)" % provider)

  try:
    # Get LLM (reuses cached instance from plan node)
    try:
      from libtbx.langchain.agent.graph_nodes import get_planning_llm
    except ImportError:
      from agent.graph_nodes import get_planning_llm

    llm, error = get_planning_llm(provider)
    if llm is None:
      print("  [THINK] WARNING: get_planning_llm returned None: %s" %
            (error or "unknown error"))
      return None

    try:
      from langchain_core.messages import SystemMessage, HumanMessage
    except ImportError:
      print("  [THINK] WARNING: langchain_core not available")
      # Fallback for standalone testing
      return None

    messages = [
      SystemMessage(content=system_msg),
      HumanMessage(content=user_msg),
    ]

    # Try with rate limit handler
    handler = _get_rate_handler(provider)
    if handler:
      def make_call():
        return llm.invoke(messages)
      response = handler.call_with_retry(make_call, lambda msg: None)
    else:
      response = llm.invoke(messages)

    print("  [THINK] LLM call succeeded")
    return response.content

  except ImportError as e:
    # Missing packages (standalone testing) -- degrade silently
    print("  [THINK] WARNING: Import error in LLM call: %s" % str(e))
    return None
  except Exception as e:
    print("  [THINK] WARNING: LLM call failed: %s" % str(e))
    return None


def _get_rate_handler(provider):
  """Get rate limit handler for the provider, or None."""
  try:
    try:
      from libtbx.langchain.agent.rate_limit_handler import (
        get_google_handler, get_openai_handler, get_anthropic_handler)
    except ImportError:
      from agent.rate_limit_handler import (
        get_google_handler, get_openai_handler, get_anthropic_handler)

    if provider == "google":
      return get_google_handler()
    elif provider == "openai":
      return get_openai_handler()
    elif provider == "anthropic":
      return get_anthropic_handler()
  except ImportError:
    pass
  return None


# =========================================================================
# Stop-decision review (called from ai_agent.py, not from the graph)
# =========================================================================

def review_stop_decision(
    stop_reasoning,
    history,
    user_advice="",
    available_files=None,
    metrics=None,
    provider="google",
):
  """Ask the expert to review a STOP decision.

  Called from ai_agent._run_single_cycle when the planning agent
  decides to stop and use_thinking_agent is enabled.

  Args:
    stop_reasoning: The plan node's reason for stopping (str).
    history: List of cycle history dicts.
    user_advice: User instructions/guidelines.
    available_files: List of file paths.
    metrics: Dict of current metrics (r_free, etc.).
    provider: LLM provider name.

  Returns:
    Dict with 'agree_with_stop', 'analysis', 'alternatives',
    'guidance'. Returns None if the LLM call fails.
  """
  print("  [THINK] Reviewing STOP decision...")

  try:
    # Build history summary
    history_summary = _summarize_history(history or [])

    # Build files summary (basenames only)
    files_summary = ""
    if available_files:
      import os
      basenames = [os.path.basename(f) for f in available_files[:20]]
      files_summary = ", ".join(basenames)
      if len(available_files) > 20:
        files_summary += " (and %d more)" % (
          len(available_files) - 20)

    # Build metrics summary
    metrics_summary = ""
    if metrics and isinstance(metrics, dict):
      items = []
      for k in ["r_free", "r_work", "map_cc", "resolution"]:
        v = metrics.get(k)
        if v is not None:
          items.append("%s=%s" % (k, v))
      metrics_summary = ", ".join(items)

    # Build prompt
    system_msg, user_msg = build_stop_review_prompt(
      stop_reasoning=stop_reasoning,
      history_summary=history_summary,
      user_advice=user_advice,
      available_files_summary=files_summary,
      metrics_summary=metrics_summary,
    )
    print("  [THINK] Stop-review prompt: system=%d chars, "
          "user=%d chars" % (len(system_msg), len(user_msg)))

    # Call LLM (reuse existing infrastructure)
    raw_response = _call_thinking_llm(
      {"provider": provider}, system_msg, user_msg)

    if raw_response is None:
      print("  [THINK] WARNING: Stop-review LLM call "
            "returned None")
      return None

    print("  [THINK] Stop-review response: %d chars"
          % len(raw_response))

    # Parse
    review = parse_stop_review(raw_response)
    agrees = review.get("agree_with_stop", True)
    analysis = review.get("analysis", "")
    alternatives = review.get("alternatives", [])
    guidance = review.get("guidance", "")

    print("  [THINK] Stop review: agree=%s" % agrees)
    if analysis:
      display = analysis.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Stop review analysis: %s" % display)
    if alternatives:
      print("  [THINK] Stop review alternatives: %s"
            % "; ".join(str(a)[:60] for a in alternatives[:3]))
    if guidance:
      display = guidance.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Stop review guidance: %s" % display)

    return review

  except Exception as e:
    print("  [THINK] Stop review ERROR: %s" % str(e))
    print("  [THINK] %s" % traceback.format_exc()[:300])
    return None


# =========================================================================
# Pre-flight readiness check (called from ai_agent.py before cycle 1)
# =========================================================================

def preflight_check(
    user_advice="",
    available_files=None,
    experiment_type="",
    provider="google",
):
  """Check whether the workflow has everything it needs.

  Called from ai_agent.iterate_agent before the cycle loop starts,
  when use_thinking_agent is enabled. This catches missing files
  (like a ligand .cif) before the workflow wastes cycles.

  Args:
    user_advice: User instructions/goal.
    available_files: List of file paths.
    experiment_type: Known experiment type if any.
    provider: LLM provider name.

  Returns:
    Dict with 'ready', 'analysis', 'missing_files', 'warnings',
    'recommendation'. Returns None if LLM call fails.
  """
  print("  [THINK] Running pre-flight readiness check...")

  try:
    # Build prompt
    system_msg, user_msg = build_preflight_prompt(
      user_advice=user_advice,
      available_files=available_files or [],
      experiment_type=experiment_type,
    )
    print("  [THINK] Preflight prompt: system=%d chars, "
          "user=%d chars" % (len(system_msg), len(user_msg)))

    # Call LLM
    raw_response = _call_thinking_llm(
      {"provider": provider}, system_msg, user_msg)

    if raw_response is None:
      print("  [THINK] WARNING: Preflight LLM call "
            "returned None")
      return None

    print("  [THINK] Preflight response: %d chars"
          % len(raw_response))

    # Parse
    result = parse_preflight(raw_response)
    ready = result.get("ready", True)
    analysis = result.get("analysis", "")
    missing = result.get("missing_files", [])
    warnings = result.get("warnings", [])
    recommendation = result.get("recommendation", "")

    if ready:
      print("  [THINK] Preflight: READY")
    else:
      print("  [THINK] Preflight: NOT READY")
    if analysis:
      display = analysis.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Preflight analysis: %s" % display)
    if missing:
      print("  [THINK] Preflight missing: %s"
            % ", ".join(str(m)[:60] for m in missing[:5]))
    if warnings:
      print("  [THINK] Preflight warnings: %s"
            % "; ".join(str(w)[:60] for w in warnings[:3]))
    if recommendation:
      display = recommendation.strip()
      if len(display) > 200:
        display = display[:197] + "..."
      print("  [THINK] Preflight recommendation: %s" % display)

    return result

  except Exception as e:
    print("  [THINK] Preflight ERROR: %s" % str(e))
    print("  [THINK] %s" % traceback.format_exc()[:300])
    return None


def _log(state, msg):
  """Add debug message to state, matching graph_nodes pattern.

  Also prints to console so the message appears in the log file.
  """
  print("  [THINK] %s" % msg)
  debug_log = list(state.get("debug_log", []))
  debug_log.append(msg)
  return {**state, "debug_log": debug_log}


def _emit(state, **kwargs):
  """Emit a structured event, matching graph_nodes pattern."""
  events = list(state.get("events", []))
  events.append({"type": "expert_assessment", **kwargs})
  return {**state, "events": events}
