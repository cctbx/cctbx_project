"""
Thinking Agent Module (v113).

Expert crystallographer reasoning node for the PHENIX AI Agent.
Provides strategic guidance to the plan node by analyzing program
logs with domain expertise.

Entry point: run_think_node(state) — called from graph_nodes.think().
"""

from __future__ import absolute_import, division, print_function
import traceback

# --- Imports from B1, B2, B3 ---
# Use try/except for libtbx paths vs standalone testing
try:
  from libtbx.langchain.agent.strategy_memory import StrategyMemory
  from libtbx.langchain.agent.log_section_extractor import extract_sections
  from libtbx.langchain.knowledge.thinking_prompts import (
    build_thinking_prompt, parse_assessment)
except ImportError:
  from agent.strategy_memory import StrategyMemory
  from agent.log_section_extractor import extract_sections
  from knowledge.thinking_prompts import (
    build_thinking_prompt, parse_assessment)


# =========================================================================
# Strategic program list — programs whose output warrants expert review
# =========================================================================

_STRATEGIC_PROGRAMS = frozenset([
  "phenix.xtriage", "xtriage",
  "phenix.phaser", "phaser",
  "phenix.autosol", "autosol",
  "phenix.autobuild", "autobuild",
  "phenix.refine", "refine",
  "phenix.real_space_refine", "real_space_refine",
  "phenix.map_to_model", "map_to_model",
  "phenix.predict_and_build", "predict_and_build",
  "phenix.ligandfit", "ligandfit",
])


def should_think(state):
  """Decide whether the thinking LLM should engage.

  When thinking_level is set, the user has accepted the
  extra LLM cost. With validation, KB rules, and R-free
  trend now available, the thinking agent provides
  substantial value during refinement — not just after
  strategic milestones.

  Triggers:
  - After any strategic program (including refine)
  - After failures
  - When R-free is stalled (3+ cycles w/o improvement)

  NOT invoked on:
  - First cycle (no program output yet)

  Args:
    state: Graph state dict.

  Returns:
    True if the thinking LLM should be called.
  """
  history = state.get("history", [])

  # First cycle: no program output to analyze yet.
  # The first useful think is AFTER xtriage/first program runs.
  if not history:
    return False

  last = history[-1]

  # Current program: prefer log_analysis (set by
  # perceive for the current cycle's output).
  log_analysis = state.get("log_analysis") or {}
  program = (
    log_analysis.get("program")
    or log_analysis.get("detected_program")
    or last.get("program", "")
  )

  # Strategic programs
  if any(p in program for p in _STRATEGIC_PROGRAMS):
    return True

  # Last cycle failed
  status = last.get("status", "")
  if status and "FAIL" in str(status).upper():
    return True

  # R-free stalled
  memory_dict = state.get("strategy_memory") or {}
  if memory_dict:
    memory = StrategyMemory.from_dict(memory_dict)
    if memory.metrics_stalled():
      return True

  return False


def run_think_node(state):
  """Main entry point — called from graph_nodes.think().

  Builds context, calls the thinking LLM, parses the response,
  and enriches the state for the plan node.

  Respects thinking_level:
    "basic": Log analysis + strategy memory only.
    "advanced": Full pipeline (validation + KB + metadata).

  Args:
    state: Graph state dict.

  Returns:
    Updated state dict.
  """
  # Safety: if not enabled, return unchanged
  thinking_level = state.get("thinking_level")
  if not thinking_level:
    return state

  # Check if we should think this cycle
  if not should_think(state):
    state = _log(state, "THINK: Skipping (routine step)")
    return state

  state = _log(
    state,
    "THINK: Expert reasoning engaged (level=%s)"
    % thinking_level
  )

  try:
    # Build context — level controls depth
    context = _build_thinking_context(
      state, thinking_level
    )

    # Store validation report in state for debugging
    # and potential GUI display (advanced only)
    vr = context.get("validation_report", "")
    if vr:
      state = {**state, "validation_report": vr}
      # Log the most informative line (usually the
      # R-factor line, which is line [1] after header)
      vr_lines = vr.split("\n")
      summary = (
        vr_lines[1][:80] if len(vr_lines) > 1
        else vr_lines[0][:80]
      )
      state = _log(
        state, "THINK: Validation: %s" % summary
      )

    # Build file metadata (Phase C) — advanced only
    if thinking_level == "advanced":
      state = _update_file_metadata(state, context)
      # Sync context so the prompt reflects the
      # current cycle's metadata, not the previous.
      context["file_metadata"] = (
        state.get("file_metadata") or {}
      )

    # Build prompt
    memory_dict = state.get("strategy_memory") or {}
    system_msg, user_msg = build_thinking_prompt(
      context, memory_dict
    )

    # Call LLM
    raw_response = _call_thinking_llm(state, system_msg, user_msg)

    if raw_response is None:
      state = _log(state, "THINK: LLM call failed, continuing without")
      return state

    # Parse response
    assessment = parse_assessment(raw_response)
    state = _log(state,
      "THINK: Assessment — action=%s confidence=%s"
      % (assessment["action"], assessment["confidence"]))

    # Build a validation summary for the event by
    # stripping the "=== VALIDATION ===" header line
    # (redundant inside the Expert Assessment block).
    validation_summary = ""
    vr_text = context.get("validation_report", "")
    if vr_text:
      vr_lines = vr_text.split("\n")
      kept = [
        ln for ln in vr_lines
        if ln.strip() and not ln.startswith("===")
      ]
      validation_summary = "\n".join(kept)

    # R-free trend line (one line, e.g. "0.42 -> 0.30 -> 0.25")
    rfree_trend_text = ""
    rfree_trend = context.get("r_free_trend", [])
    if rfree_trend and len(rfree_trend) >= 2:
      rfree_trend_text = " -> ".join(
        "%.3f" % v for v in rfree_trend
      )

    # Emit structured event for GUI Events panel
    state = _emit(state,
      action=assessment["action"],
      confidence=assessment["confidence"],
      thinking_level=thinking_level,
      analysis=assessment.get("analysis", ""),
      guidance=assessment.get("guidance", ""),
      concerns=assessment.get("concerns", []),
      validation_summary=validation_summary,
      rfree_trend=rfree_trend_text,
      kb_rules_matched=bool(
        context.get("kb_rules_text", "")
      ),
    )

    # Expert recommends stopping?
    if assessment["action"] == "stop":
      state = _log(state, "THINK: Expert recommends STOP")
      cycle = state.get("cycle_number", 1)
      abort_msg = (
        "Expert assessment: %s"
        % assessment["analysis"][:300]
      )
      # Emit stop_decision so the event formatter
      # shows it (plan() will short-circuit and
      # never get to emit its own stop_decision).
      stop_reason = (
        "expert: %s"
        % assessment["analysis"][:200]
      )
      events = list(state.get("events", []))
      events.append({
        "type": "stop_decision",
        "stop": True,
        "reason": stop_reason,
      })
      state = {**state, "events": events}
      return {
        **state,
        "command": "STOP",
        "stop": True,
        "stop_reason": stop_reason,
        "abort_message": abort_msg,
        "expert_assessment": assessment,
        "strategy_memory": _update_memory(
          memory_dict, assessment, cycle),
      }

    # Inject guidance into user_advice for PLAN and BUILD
    guidance = assessment.get("guidance", "")
    if guidance:
      enriched_advice = (
        "[Expert assessment] %s\n\n%s"
        % (guidance, state.get("user_advice", "")))
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
    # Never crash the workflow — degrade gracefully
    state = _log(state,
      "THINK: Error (continuing without): %s" % str(e))
    state = _log(state, "THINK: %s" % traceback.format_exc()[:500])
    return state


# =========================================================================
# Internal helpers
# =========================================================================

def _build_thinking_context(state, thinking_level="advanced"):
  """Assemble context dict for the thinking prompt.

  Gathers log sections, metrics, workflow state, history,
  and (for advanced mode) structural validation results
  into the format expected by build_thinking_prompt().

  Args:
    state: Graph state dict.
    thinking_level: "basic" or "advanced".
      basic: log sections, metrics, history, r_free_trend
      advanced: all of the above plus validation, KB,
        and file metadata.
  """
  history = state.get("history", [])
  last = history[-1] if history else {}

  # Current program: prefer log_analysis (set by
  # perceive for the current cycle) over history[-1]
  # (which may be the PREVIOUS cycle's program if
  # history hasn't been updated yet).
  log_analysis = state.get("log_analysis") or {}
  program = (
    log_analysis.get("program")
    or log_analysis.get("detected_program")
    or last.get("program", "")
  )

  # Extract informative log sections
  log_text = state.get("log_text", "")
  log_sections = ""
  if log_text and program:
    log_sections = extract_sections(log_text, program)
  elif log_text:
    # No program name — use generic tail extraction
    log_sections = extract_sections(log_text, "unknown")

  # Get experiment type and workflow state
  session_info = state.get("session_info") or {}
  workflow_state_dict = state.get("workflow_state") or {}

  # Extract recent metrics from log_analysis
  metrics = log_analysis if log_analysis else {}

  # Build brief history summary
  history_summary = _summarize_history(history)

  # --- R-free trend (both levels) ---
  r_free_trend = _collect_r_free_trend(history)
  current_rf = metrics.get("r_free")
  if current_rf is not None:
    try:
      r_free_trend.append(float(current_rf))
    except (ValueError, TypeError):
      pass

  # --- Basic context (shared by both levels) ---
  context = {
    "log_sections": log_sections,
    "program_name": program,
    "cycle_number": state.get("cycle_number", 1),
    "experiment_type": session_info.get(
      "experiment_type", "unknown"
    ),
    "workflow_state": workflow_state_dict.get(
      "state", "unknown"
    ),
    "metrics": metrics,
    "user_advice": state.get("user_advice", ""),
    "history_summary": history_summary,
    "r_free_trend": r_free_trend,
    # Advanced fields default to empty
    "validation_report": "",
    "validation_result": None,
    "validated_model_path": None,
    "file_metadata": {},
    "kb_rules_text": "",
  }

  # --- Advanced mode additions ---
  if thinking_level == "advanced":
    # Structural validation (Phase A)
    validation_report, validation_result, model_path = (
      _run_structural_validation(
        state, session_info, metrics, program
      )
    )
    context["validation_report"] = validation_report
    context["validation_result"] = validation_result
    context["validated_model_path"] = model_path

    # File metadata from state (Phase C)
    context["file_metadata"] = (
      state.get("file_metadata") or {}
    )

    # Expert knowledge base query (Phase D)
    xtriage_results = _extract_xtriage(history)
    context["kb_rules_text"] = _query_knowledge_base(
      validation_result=validation_result,
      log_metrics=metrics,
      workflow_stage=workflow_state_dict.get(
        "state", "unknown"
      ),
      program_name=program,
      resolution=metrics.get("resolution"),
      experiment_type=session_info.get(
        "experiment_type", "xray"
      ),
      r_free_trend=r_free_trend,
      xtriage_results=xtriage_results,
    )

  return context


def _summarize_history(history):
  """Create a brief text summary of cycle history."""
  if not history:
    return "(first cycle)"
  lines = []
  for h in history[-5:]:  # Last 5 cycles
    prog = h.get("program", "?")
    status = h.get("status", "?")
    lines.append("C%s: %s (%s)" % (
      h.get("cycle_number", "?"), prog, status))
  return "; ".join(lines)


def _run_structural_validation(state, session_info,
                               metrics, program):
  """Run headless validation and return results.

  Extracts the best model path from session_info,
  calls run_validation(), and formats the result as
  a compact text block for the thinking agent prompt.

  Returns:
    tuple of (report_str, validation_result, model_path)
    where validation_result is the raw dict (or None)
    and model_path is the model that was validated.
    Returns ("", None, None) if unavailable.
  """
  try:
    try:
      from libtbx.langchain.agent.validation_inspector \
        import run_validation
      from libtbx.langchain.agent.format_validation \
        import format_validation_report
    except ImportError:
      from agent.validation_inspector \
        import run_validation
      from agent.format_validation \
        import format_validation_report
  except ImportError:
    # Modules not available — skip validation
    return ("", None, None)

  # Get model path from best_files
  best_files = session_info.get("best_files") or {}
  model_path = best_files.get("model")
  if not model_path:
    return ("", None, None)

  # Get data and map coefficients paths
  data_path = best_files.get("data")
  map_coeffs = best_files.get("map_coefficients")
  exp_type = session_info.get(
    "experiment_type", "xray"
  )

  # Run validation (never raises)
  validation_result = run_validation(
    model_path=model_path,
    data_path=data_path,
    map_coeffs_path=map_coeffs,
    experiment_type=exp_type,
  )

  # Get R-free trend for the report
  cycle_number = state.get("cycle_number", 1)
  history = state.get("history", [])
  prev_rf = _get_prev_r_free(history)
  start_rf = _get_start_r_free(history)

  # Format
  report = format_validation_report(
    validation_result=validation_result,
    log_metrics=metrics,
    cycle_number=cycle_number,
    program_name=program,
    prev_r_free=prev_rf,
    start_r_free=start_rf,
  )

  return (report, validation_result, model_path)


def _collect_r_free_trend(history):
  """Collect R-free values from cycle history.

  History entries from get_history_for_agent() use
  key "analysis" (flat dict), not "log_analysis".

  Returns:
    list of float, one per cycle that had R-free.
  """
  trend = []
  for h in history:
    # Try "analysis" first (from session history),
    # then "log_analysis" (from in-graph state)
    la = h.get("analysis") or h.get(
      "log_analysis", {}
    )
    if not la:
      continue
    # Handle flat dict (r_free at top level)
    # or nested dict (r_free under "metrics")
    rf = la.get("r_free")
    if rf is None and isinstance(la.get("metrics"), dict):
      rf = la["metrics"].get("r_free")
    if rf is not None:
      try:
        trend.append(float(rf))
      except (ValueError, TypeError):
        pass
  return trend


def _get_prev_r_free(history):
  """R-free from the most recent cycle that had one.

  Note: history contains completed cycles only.
  The current cycle (whose log we're analyzing) is
  NOT yet in history. Scans backwards to find the
  most recent R-free, skipping cycles like ligandfit
  or polder that don't produce R-free.
  """
  if not history:
    return None
  for entry in reversed(history):
    la = entry.get("analysis") or entry.get(
      "log_analysis", {}
    )
    if not la:
      continue
    rf = la.get("r_free")
    if rf is None and isinstance(
      la.get("metrics"), dict
    ):
      rf = la["metrics"].get("r_free")
    if rf is not None:
      try:
        return float(rf)
      except (ValueError, TypeError):
        pass
  return None


def _get_start_r_free(history):
  """R-free from the first cycle that had R-free.

  Scans forward through history to find the first
  entry with an R-free value. This is typically the
  first refinement cycle, not history[0] (which is
  usually xtriage).
  """
  if not history:
    return None
  for entry in history:
    la = entry.get("analysis") or entry.get(
      "log_analysis", {}
    )
    if not la:
      continue
    rf = la.get("r_free")
    if rf is None and isinstance(
      la.get("metrics"), dict
    ):
      rf = la["metrics"].get("r_free")
    if rf is not None:
      try:
        return float(rf)
      except (ValueError, TypeError):
        pass
  return None


# =========================================================
# Xtriage extraction
# =========================================================

def _extract_xtriage(history):
  """Extract xtriage results from history.

  Xtriage is typically the first cycle. Looks for
  history entries where program is "phenix.xtriage"
  and returns any xtriage-specific fields from the
  analysis dict.

  Returns:
    dict with xtriage fields, or None.
  """
  if not history:
    return None
  for entry in history:
    prog = entry.get("program", "")
    if "xtriage" not in prog:
      continue
    la = entry.get("analysis") or entry.get(
      "log_analysis", {}
    )
    if not la:
      continue
    # Build xtriage result from available fields
    result = {}
    for key in ("twin_fraction", "twin_law",
                "tncs_detected", "anisotropy_detected",
                "ice_rings", "space_group",
                "resolution"):
      if key in la:
        result[key] = la[key]
    if result:
      return result
  return None


# =========================================================
# Expert Knowledge Base (Phase D)
# =========================================================

_kb_instance = None


def _get_kb():
  """Lazy-load the expert knowledge base (singleton).

  Returns ExpertKnowledgeBase instance or None if
  not available.
  """
  global _kb_instance
  if _kb_instance is not None:
    return _kb_instance

  import os
  try:
    try:
      from libtbx.langchain.knowledge.kb_loader \
        import ExpertKnowledgeBase
    except ImportError:
      from knowledge.kb_loader \
        import ExpertKnowledgeBase
  except ImportError:
    return None

  # Find the YAML file relative to this module
  this_dir = os.path.dirname(
    os.path.abspath(__file__)
  )
  # Try v2 first, then fall back to v1
  _CANDIDATE_PATHS = [
    os.path.join(
      this_dir, "..", "knowledge",
      "expert_knowledge_base_v2.yaml",
    ),
    os.path.join(
      this_dir, "..", "knowledge",
      "expert_knowledge_base.yaml",
    ),
    os.path.join(
      this_dir, "data",
      "expert_knowledge_base_v2.yaml",
    ),
    os.path.join(
      this_dir, "data",
      "expert_knowledge_base.yaml",
    ),
  ]
  kb_path = None
  for candidate in _CANDIDATE_PATHS:
    if os.path.isfile(candidate):
      kb_path = candidate
      break
  if kb_path is None:
    return None

  try:
    _kb_instance = ExpertKnowledgeBase(kb_path)
  except Exception:
    return None
  return _kb_instance


def _query_knowledge_base(
  validation_result=None,
  log_metrics=None,
  workflow_stage=None,
  program_name=None,
  resolution=None,
  experiment_type=None,
  r_free_trend=None,
  xtriage_results=None,
):
  """Query the KB and return formatted rules text.

  Returns:
    str, e.g. "=== RELEVANT EXPERT RULES ===\n..."
    or "" if KB not available or no matches.
  """
  kb = _get_kb()
  if kb is None:
    return ""

  try:
    try:
      from libtbx.langchain.agent.kb_tags \
        import derive_tags
    except ImportError:
      from agent.kb_tags import derive_tags
  except ImportError:
    return ""

  category, tags = derive_tags(
    validation_result=validation_result,
    log_metrics=log_metrics,
    workflow_stage=workflow_stage,
    program_name=program_name,
    resolution=resolution,
    experiment_type=experiment_type,
    r_free_trend=r_free_trend,
    xtriage_results=xtriage_results,
  )

  # When diagnostic tags are present, prioritize
  # stopping entries.
  # NOTE: r_free_high (>= 0.40) is deliberately
  # excluded — it is normal in early refinement
  # (cycle 1-3 post-MR). It only becomes diagnostic
  # when paired with plateau/stuck, which have their
  # own entries in this set.
  _DIAGNOSTIC_TAGS = frozenset([
    "r_free_stuck", "plateau", "r_free_very_high",
    "geometry_terrible", "geometry_poor",
    "overfitting",
  ])
  if _DIAGNOSTIC_TAGS & set(tags):
    entries = kb.query(
      category="stopping",
      tags=tags,
      max_results=2,
    )
    # Always add at least 1 entry from all categories
    # so we mix in relevant thresholds/pitfalls, not
    # just stopping advice.
    seen = set(e.get("id") for e in entries)
    extra = kb.query(
      category=None,
      tags=tags,
      max_results=5,
    )
    for e in extra:
      if e.get("id") not in seen:
        entries.append(e)
        seen.add(e.get("id"))
      if len(entries) >= 3:
        break
  else:
    # Try category-filtered first
    entries = kb.query(
      category=category,
      tags=tags,
      max_results=3,
    )
    # Fall back to tag-only if category filter
    # returned nothing (common with v2 KB whose
    # categories differ from workflow stages).
    if not entries:
      entries = kb.query(
        category=None,
        tags=tags,
        max_results=3,
      )

  if not entries:
    return ""

  rules_text = kb.format_for_prompt(
    entries, max_chars=3000
  )
  if not rules_text:
    return ""

  return (
    "=== RELEVANT EXPERT RULES ===\n" + rules_text
  )


def _update_file_metadata(state, context):
  """Build and store file metadata from validation.

  Uses the raw validation_result and metrics from
  _build_thinking_context to create a metadata entry
  for the validated model, and stores it in state.

  Args:
    state: Graph state dict.
    context: Dict from _build_thinking_context.

  Returns:
    Updated state dict.
  """
  model_path = context.get("validated_model_path")
  vr = context.get("validation_result")
  if not model_path:
    return state

  try:
    try:
      from libtbx.langchain.agent.file_metadata \
        import build_file_metadata
    except ImportError:
      from agent.file_metadata \
        import build_file_metadata
  except ImportError:
    return state

  meta = build_file_metadata(
    file_path=model_path,
    validation_result=vr,
    log_metrics=context.get("metrics"),
    program_name=context.get("program_name"),
    cycle_number=context.get("cycle_number"),
  )

  # Merge into existing file_metadata
  fm = dict(state.get("file_metadata") or {})
  fm[model_path] = meta
  return {**state, "file_metadata": fm}


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

  try:
    # Get LLM (reuses cached instance from plan node)
    try:
      from libtbx.langchain.agent.graph_nodes import get_planning_llm
    except ImportError:
      from agent.graph_nodes import get_planning_llm

    llm, error = get_planning_llm(provider)
    if llm is None:
      return None

    try:
      from langchain_core.messages import SystemMessage, HumanMessage
    except ImportError:
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

    return response.content

  except ImportError:
    # Missing packages (standalone testing) — degrade silently
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


def _log(state, msg):
  """Add debug message to state, matching graph_nodes pattern."""
  debug_log = list(state.get("debug_log", []))
  debug_log.append(msg)
  return {**state, "debug_log": debug_log}


def _emit(state, **kwargs):
  """Emit a structured event, matching graph_nodes pattern."""
  events = list(state.get("events", []))
  events.append({"type": "expert_assessment", **kwargs})
  return {**state, "events": events}
