"""
Thinking Agent Module (v113).

Expert crystallographer reasoning node for the PHENIX AI Agent.
Provides strategic guidance to the plan node by analyzing program
logs with domain expertise.

Entry point: run_think_node(state) — called from graph_nodes.think().
"""

from __future__ import absolute_import, division, print_function
import os
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

# --- Goal-directed agent imports (v114) ---
try:
  from libtbx.langchain.agent.structure_model \
    import StructureModel, Hypothesis
  from libtbx.langchain.agent.validation_history \
    import ValidationHistory
except ImportError:
  try:
    from agent.structure_model import (
      StructureModel, Hypothesis,
    )
    from agent.validation_history \
      import ValidationHistory
  except ImportError:
    # Graceful fallback — goal-directed features
    # disabled but agent still works.
    StructureModel = None
    ValidationHistory = None
    Hypothesis = None


# =========================================================
# Strategic program list — programs whose output
# warrants expert review. Loaded from programs.yaml
# (categories: analysis, model_building, refinement,
# ligand). Falls back to a minimal set if YAML
# unavailable.
# =========================================================

_strategic_programs_cache = None


def _get_strategic_programs():
  """Load strategic programs from YAML config.

  Returns a frozenset of program names (with and
  without phenix. prefix) whose output warrants
  expert-level review.

  Cached after first call.
  """
  global _strategic_programs_cache
  if _strategic_programs_cache is not None:
    return _strategic_programs_cache

  strategic_categories = frozenset([
    "analysis", "model_building", "refinement",
    "ligand",
  ])

  try:
    try:
      from libtbx.langchain.knowledge.yaml_loader \
        import get_all_programs
    except ImportError:
      from knowledge.yaml_loader \
        import get_all_programs
    all_progs = get_all_programs()
    names = set()
    for prog in all_progs:
      cat = prog.get("category", "")
      if cat in strategic_categories:
        name = prog.get("name", "")
        if name:
          names.add(name)
          # Also add short name without phenix.
          if name.startswith("phenix."):
            names.add(name[7:])
    if names:
      _strategic_programs_cache = frozenset(names)
      return _strategic_programs_cache
  except Exception:
    pass

  # Fallback: load from programs.yaml directly
  try:
    import os
    import yaml
    yaml_path = os.path.join(
      os.path.dirname(os.path.abspath(__file__)),
      "..", "knowledge", "programs.yaml",
    )
    if os.path.isfile(yaml_path):
      with open(yaml_path) as f:
        data = yaml.safe_load(f)
      if isinstance(data, dict):
        names = set()
        for name, info in data.items():
          if not isinstance(info, dict):
            continue
          cat = info.get("category", "")
          if cat in strategic_categories:
            names.add(name)
            if name.startswith("phenix."):
              names.add(name[7:])
        if names:
          _strategic_programs_cache = frozenset(
            names
          )
          return _strategic_programs_cache
  except Exception:
    pass

  # Minimal fallback — should never reach here
  # in a properly configured installation
  _strategic_programs_cache = frozenset()
  return _strategic_programs_cache


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
  strategic = _get_strategic_programs()
  if strategic and any(
    p in program for p in strategic
  ):
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

    # --- Emit diagnostic for advanced mode so the
    # user sees what the agent's "eyes" found (or
    # why they didn't run), BEFORE the LLM call.
    # This ensures visibility even if the LLM fails.
    if thinking_level in ("advanced", "expert"):
      skip_reason = context.get(
        "validation_skip_reason", ""
      )
      if model_path and not skip_reason:
        diag_parts.append(
          "Validated: %s"
          % os.path.basename(model_path)
        )
      elif skip_reason:
        diag_parts.append(
          "Validation skipped: %s" % skip_reason
        )
      # KB status
      kb_text = context.get("kb_rules_text", "")
      if kb_text:
        # Count rules (each starts with "[")
        n_rules = kb_text.count("[")
        diag_parts.append(
          "%d KB rule(s) matched" % n_rules
        )
      else:
        diag_parts.append("no KB rules matched")
      # File metadata status
      fm = context.get("file_metadata") or {}
      if fm:
        diag_parts.append(
          "%d file(s) tracked" % len(fm)
        )
      # Emit as a visible debug_log entry —
      # the event formatter shows these at NORMAL
      # inside the Expert Assessment block, but
      # _log entries are VERBOSE only. Use _emit
      # with a dedicated type that the formatter
      # can render, or append to the later
      # expert_assessment event. For now, ensure
      # at minimum a _log so verbose users see it.
      diag_msg = (
        "THINK: Advanced context: %s"
        % "; ".join(diag_parts)
      )
      state = _log(state, diag_msg)

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

    # Build file metadata (Phase C) — advanced/expert only
    if thinking_level in ("advanced", "expert"):
      state = _update_file_metadata(state, context)
      # Sync context so the prompt reflects the
      # current cycle's metadata, not the previous.
      context["file_metadata"] = (
        state.get("file_metadata") or {}
      )

    # Persist Structure Model + Validation History
    # to state BEFORE the LLM call so structural
    # knowledge survives even if the LLM fails.
    # (v114: goal-directed agent)
    # Fix: was "advanced" only — expert mode also builds SM
    # so it must persist here too, otherwise session.data
    # never gets structure_model and the Results panel
    # shows "No structure model available / INCOMPLETE".
    if thinking_level in ("advanced", "expert"):
      sm = context.get("structure_model")
      if sm is not None:
        state = {
          **state,
          "structure_model": sm.to_dict(),
        }
      vh = context.get("validation_history")
      if vh is not None:
        state = {
          **state,
          "validation_history": vh.to_dict(),
        }

    # Build prompt
    memory_dict = state.get("strategy_memory") or {}
    system_msg, user_msg = build_thinking_prompt(
      context, memory_dict
    )

    # Call LLM
    raw_response = _call_thinking_llm(state, system_msg, user_msg)

    if raw_response is None:
      state = _log(state,
        "THINK: LLM call failed, continuing without")
      # Still emit the validation / structure data so
      # the user sees what advanced mode found, even
      # though the LLM didn't respond.
      sm_summary = context.get(
        "structure_model_summary", ""
      )
      if thinking_level in ("advanced", "expert") and (
        vr or sm_summary
      ):
        vr_lines_clean = []
        if vr:
          vr_lines_clean = [
            ln for ln in vr.split("\n")
            if ln.strip()
            and not ln.startswith("===")
          ]
        state = _emit(state,
          action="let_run",
          confidence="unknown",
          thinking_level=thinking_level,
          analysis="LLM call failed; validation"
                   " data shown for reference.",
          guidance="",
          concerns=[],
          validation_summary="\n".join(
            vr_lines_clean),
          rfree_trend="",
          kb_rules_matched=bool(
            context.get("kb_rules_text", "")),
          structure_model_summary=context.get(
            "structure_model_summary", ""),
          current_problems=context.get(
            "current_problems", []),
        )
      return state

    # Parse response
    assessment = parse_assessment(raw_response)
    state = _log(state,
      "THINK: Assessment — action=%s confidence=%s"
      % (assessment["action"], assessment["confidence"]))

    # --- Extract hypothesis from LLM response ---
    # If the LLM proposed a hypothesis in its JSON,
    # create a Hypothesis object and add it to the
    # StructureModel. The hypothesis_evaluator will
    # manage its lifecycle from here.
    _extract_hypothesis_from_assessment(
      assessment, context, state
    )

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
      validation_skip_reason=context.get(
        "validation_skip_reason", ""),
      validated_model=os.path.basename(
        context.get("validated_model_path", "") or ""
      ),
      rfree_trend=rfree_trend_text,
      kb_rules_matched=bool(
        context.get("kb_rules_text", "")
      ),
      structure_model_summary=context.get(
        "structure_model_summary", ""),
      current_problems=context.get(
        "current_problems", []),
    )

    # Enrich the assessment dict with advanced-mode data
    # so the GUI progress panel (AIAgent.py) can display
    # structural validation alongside the LLM analysis.
    # The GUI reads expert_assessment from the callback
    # data, which comes from state["expert_assessment"].
    if thinking_level in ("advanced", "expert"):
      assessment["validation_summary"] = (
        validation_summary)
      assessment["validation_skip_reason"] = (
        context.get("validation_skip_reason", ""))
      assessment["validated_model"] = os.path.basename(
        context.get("validated_model_path", "") or "")
      assessment["rfree_trend"] = rfree_trend_text
      assessment["kb_rules_matched"] = bool(
        context.get("kb_rules_text", ""))
      assessment["structure_model_summary"] = (
        context.get("structure_model_summary", ""))
      assessment["current_problems"] = (
        context.get("current_problems", []))

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
      #
      # P1B: use classified stop_reason_code when
      # available; fall back to freeform analysis.
      _src = assessment.get("stop_reason_code")
      if _src:
        stop_reason = _src
      else:
        stop_reason = (
          "expert: %s"
          % assessment["analysis"][:200]
        )
      # P1B: build think_stop_override for output_node
      _think_stop_override = None
      if _src:
        _think_stop_override = {
          "code": _src,
          "analysis": assessment.get(
            "analysis", "")[:300],
        }
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
        "think_stop_override": _think_stop_override,
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

    # P1B: forward file_overrides from expert assessment
    # to BUILD node via think_file_overrides state key.
    _file_overrides = assessment.get(
      "file_overrides", {})
    if not isinstance(_file_overrides, dict):
      _file_overrides = {}

    cycle = state.get("cycle_number", 1)
    return {
      **state,
      "user_advice": enriched_advice,
      "expert_assessment": assessment,
      "think_file_overrides": _file_overrides,
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
    thinking_level: "basic", "advanced", or "expert".
      basic: log sections, metrics, history, r_free_trend
      advanced/expert: all of the above plus validation,
        KB, file metadata, and StructureModel.
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
    "valid_programs": workflow_state_dict.get(
      "valid_programs", []),
    "plan_phase": "",
    "plan_goal": "",
    "stop_after": "",
    "metrics": metrics,
    "user_advice": state.get("user_advice", ""),
    "history_summary": history_summary,
    "r_free_trend": r_free_trend,
    "recent_failures": _collect_recent_failures(
      history),
    # Advanced fields default to empty
    "validation_report": "",
    "validation_result": None,
    "validated_model_path": None,
    "file_metadata": {},
    "kb_rules_text": "",
    "file_inventory": "",
  }

  # File inventory: give the expert a category-grouped
  # view of available files so it can see what inputs
  # exist without guessing from filenames alone.
  _cat_files = workflow_state_dict.get(
    "categorized_files", {})
  if _cat_files:
    # Group into meaningful labels, skip empty categories
    _inv_map = [
      ("Models", [
        "model", "refined", "phaser_output",
        "rsr_output", "predicted",
        "autobuild_output"]),
      ("Sequences", ["sequence"]),
      ("Reflection data", [
        "data_mtz", "phased_data_mtz"]),
      ("Map coefficients", [
        "map_coeffs_mtz", "refine_map_coeffs",
        "denmod_map_coeffs"]),
      ("Maps", [
        "full_map", "half_map"]),
      ("Ligands", [
        "ligand_pdb", "ligand_cif"]),
    ]
    _inv_lines = []
    for label, cats in _inv_map:
      fnames = []
      for cat in cats:
        for fpath in _cat_files.get(cat, []):
          bn = os.path.basename(fpath)
          tag = cat.replace("_", " ")
          entry = "%s (%s)" % (bn, tag)
          if entry not in fnames:
            fnames.append(entry)
      if fnames:
        _inv_lines.append(
          "  %s: %s" % (label, ", ".join(fnames)))
    if _inv_lines:
      context["file_inventory"] = (
        "\n".join(_inv_lines))

  # Plan phase and directives (so guidance aligns
  # with what the agent is trying to accomplish)
  _directives = state.get("directives", {})
  _wf_prefs = _directives.get(
    "workflow_preferences", {})
  _stop_conds = _directives.get(
    "stop_conditions", {})
  # Plan phase from prefer_programs (set by
  # plan_to_directives when expert mode active)
  _prefer = _wf_prefs.get("prefer_programs", [])
  if _prefer:
    context["plan_phase"] = (
      "Current plan prefers: %s"
      % ", ".join(_prefer[:3])
    )
  # Stop condition
  _after = _stop_conds.get("after_program", "")
  if _after:
    context["stop_after"] = _after
  # Plan goal from session_info (set by ai_agent.py)
  _plan_data = (session_info or {}).get("plan")
  if not _plan_data:
    _plan_data = state.get("session_info", {}).get(
      "plan")
  # Not worth full deserialization — just check for
  # a goal string
  if isinstance(_plan_data, dict):
    context["plan_goal"] = (
      _plan_data.get("goal", ""))

  # ASU copy count (copies feature): expose so THINK prompt can reason
  # about whether Phaser found the expected number of copies.
  _asu_copies = (session_info or {}).get("asu_copies")
  if _asu_copies:
    context["asu_copies"] = _asu_copies
  if thinking_level in ("advanced", "expert"):
    # Structural validation (Phase A)
    validation_report, validation_result, model_path, \
      validation_skip = (
        _run_structural_validation(
          state, session_info, metrics, program
        )
      )
    context["validation_report"] = validation_report
    context["validation_result"] = validation_result
    context["validated_model_path"] = model_path
    context["validation_skip_reason"] = validation_skip

    # --- Structure Model + Validation History (v114) ---
    # Restore from state or create fresh. Updated every
    # cycle so structural knowledge accumulates.
    structure_model, validation_hist = (
      _update_structure_model(
        state, validation_result, metrics,
        program, session_info,
      )
    )
    context["structure_model"] = structure_model
    context["validation_history"] = validation_hist

    # Replace ad-hoc validation_report with Structure
    # Model summary for the THINK prompt (richer, and
    # includes cross-cycle knowledge).
    sm_summary = ""
    if structure_model is not None:
      sm_summary = structure_model.get_summary(
        detail_level="normal"
      )
    if sm_summary:
      context["structure_model_summary"] = sm_summary
    else:
      context["structure_model_summary"] = ""

    # Current problems for the THINK prompt
    current_problems = []
    if structure_model is not None:
      current_problems = (
        structure_model.get_current_problems()
      )
    context["current_problems"] = current_problems

    # Hypothesis prompt for the THINK prompt (Phase 4)
    hypothesis_prompt = ""
    if structure_model is not None:
      try:
        try:
          from libtbx.langchain.agent \
            .hypothesis_evaluator import (
              build_hypothesis_prompt,
            )
        except ImportError:
          from agent.hypothesis_evaluator import (
            build_hypothesis_prompt,
          )
        hypothesis_prompt = build_hypothesis_prompt(
          structure_model
        )
      except ImportError:
        pass
      except Exception:
        pass
    context["hypothesis_prompt"] = hypothesis_prompt

    # Error classification from PERCEIVE (Enhancement B, v115)
    # These feed into the "=== PREVIOUS FAILURE ===" section
    # of thinking_prompts.py expert prompt.
    context["error_classification"] = (
      state.get("error_classification") or ""
    )
    context["failure_count"] = (
      state.get("failure_count") or 0
    )

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
  """Create a brief text summary of cycle history.

  Includes input files (extracted from the command string)
  and output files so the expert can detect wrong-file-for-
  refinement errors.
  """
  if not history:
    return "(first cycle)"
  lines = []
  for h in history[-5:]:  # Last 5 cycles
    prog = h.get("program", "?")
    status = h.get("result", h.get("status", "?"))
    cycle = h.get("cycle_number", "?")
    parts = ["C%s: %s (%s)" % (cycle, prog, status)]

    # Input files — basenames from the command string
    cmd = h.get("command", "")
    if cmd:
      import shlex
      try:
        tokens = shlex.split(cmd)
      except ValueError:
        tokens = cmd.split()
      # Skip the program token itself
      input_names = []
      for tok in tokens[1:]:
        # Only include tokens that look like filenames
        if ("." in os.path.basename(tok)
            and not tok.startswith("-")):
          input_names.append(os.path.basename(tok))
      if input_names:
        parts.append(
            "  Input: %s" % ", ".join(input_names))

    # Output files — basenames from output_files list
    out_files = h.get("output_files") or []
    if out_files:
      out_names = [os.path.basename(f) for f in out_files]
      parts.append(
          "  Output: %s" % ", ".join(out_names))

    lines.append("\n".join(parts))
  return "\n".join(lines)


def _collect_recent_failures(history):
  """Collect recent program failures for thinking.

  Returns a list of dicts with program, error, and
  cycle number for the last N failures. This gives
  the thinking agent explicit awareness of what has
  already been tried and failed, so it can recommend
  a different approach.
  """
  failures = []
  if not history:
    return failures
  for h in history[-5:]:
    if not isinstance(h, dict):
      continue
    result = str(h.get("result", ""))
    result_upper = result.upper()
    is_fail = (
      result_upper.startswith("FAIL")
      or ("ERROR" in result_upper
          and "WITHOUT ERROR" not in result_upper))
    if is_fail:
      failures.append({
        "cycle": h.get("cycle_number", "?"),
        "program": h.get("program", "?"),
        "error": result[:200],
        "command": h.get("command", "")[:200],
      })
  return failures


def _run_structural_validation(state, session_info,
                               metrics, program):
  """Run headless validation and return results.

  Extracts the best model path from session_info,
  calls run_validation(), and formats the result as
  a compact text block for the thinking agent prompt.

  Returns:
    tuple of (report_str, validation_result, model_path,
              skip_reason)
    where validation_result is the raw dict (or None),
    model_path is the model that was validated, and
    skip_reason is a short string explaining why
    validation did not run (empty if it succeeded).
    Returns ("", None, None, reason) if unavailable.
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
    return ("", None, None,
            "validation modules not available")

  # Get model path from best_files
  best_files = session_info.get("best_files") or {}
  model_path = best_files.get("model")
  if not model_path:
    return ("", None, None,
            "no model in best_files (too early in run?)")

  if not os.path.isfile(model_path):
    return ("", None, None,
            "model file not found: %s"
            % os.path.basename(model_path))

  # Get data and map coefficients paths.
  # BestFilesTracker categories:
  #   "data_mtz" — reflection data with R-free flags
  #   "map_coeffs_mtz" — map coefficients from refine
  #   "sequence" — sequence file
  data_path = best_files.get("data_mtz")
  map_coeffs = best_files.get("map_coeffs_mtz")
  exp_type = session_info.get(
    "experiment_type", "xray"
  )

  # Get sequence path from best_files or available
  # files (Step 1.3: chain completeness needs it
  # for accurate expected-residue counts).
  seq_path = best_files.get("sequence")
  if not seq_path:
    avail_files = state.get(
      "available_files",
    ) or []
    _SEQ_EXT = (".seq", ".fasta", ".fa", ".pir")
    for f in avail_files:
      if any(str(f).lower().endswith(ext)
             for ext in _SEQ_EXT):
        if os.path.isfile(str(f)):
          seq_path = str(f)
          break

  # Run validation (never raises)
  validation_result = run_validation(
    model_path=model_path,
    data_path=data_path,
    map_coeffs_path=map_coeffs,
    sequence_path=seq_path,
    experiment_type=exp_type,
  )

  if validation_result is None:
    return ("", None, model_path,
            "validation returned no results for %s"
            % os.path.basename(model_path))

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

  return (report, validation_result, model_path, "")


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


def _update_structure_model(state, validation_result,
                           log_metrics, program_name,
                           session_info):
  """Restore or create StructureModel + ValidationHistory.

  Called from _build_thinking_context (advanced mode).
  Restores from state on resume, creates fresh on first
  cycle. Updates from validation results and log metrics.

  Args:
    state: Graph state dict.
    validation_result: dict from run_validation() or None.
    log_metrics: dict with r_work, r_free, etc.
    program_name: str, e.g. "phenix.refine".
    session_info: dict with experiment_type, best_files.

  Returns:
    tuple of (StructureModel, ValidationHistory).
    Returns (None, None) if classes not available.

  Never raises.
  """
  if StructureModel is None or ValidationHistory is None:
    return (None, None)

  try:
    return _update_structure_model_inner(
      state, validation_result, log_metrics,
      program_name, session_info,
    )
  except Exception:
    import logging
    logging.getLogger(__name__).debug(
      "_update_structure_model failed",
      exc_info=True,
    )
    return (None, None)


def _update_structure_model_inner(
  state, validation_result, log_metrics,
  program_name, session_info,
):
  """Inner structure model update. May raise."""
  cycle_number = state.get("cycle_number", 1)
  history = state.get("history", [])

  # --- Restore or create StructureModel ---
  # Check state first (populated if run_ai_agent.py
  # passes structure_model= to create_initial_state).
  # Fall back to session_info (always populated by
  # ai_agent.py, even before run_ai_agent.py is
  # updated to extract it).
  sm_dict = state.get("structure_model")
  if not sm_dict or not isinstance(sm_dict, dict):
    sm_dict = (session_info or {}).get(
      "structure_model"
    )
  if sm_dict and isinstance(sm_dict, dict) \
      and sm_dict.get("_version"):
    sm = StructureModel.from_dict(sm_dict)
  else:
    sm = StructureModel()

  # --- Restore or create ValidationHistory ---
  vh_dict = state.get("validation_history")
  if not vh_dict or not isinstance(vh_dict, dict):
    vh_dict = (session_info or {}).get(
      "validation_history"
    )
  if vh_dict and isinstance(vh_dict, dict) \
      and vh_dict.get("_version"):
    vh = ValidationHistory.from_dict(vh_dict)
  else:
    vh = ValidationHistory()

  # --- Update from xtriage (if first time) ---
  # Xtriage data comes from history (it ran in an
  # earlier cycle). Only update if data_characteristics
  # is still empty (idempotent — don't overwrite on
  # every cycle).
  if sm.data_characteristics.get("resolution") is None:
    xtriage_results = _extract_xtriage(history)
    if xtriage_results:
      sm.update_from_xtriage(xtriage_results)

  # Set experiment type if not already set
  if sm.data_characteristics.get(
    "experiment_type"
  ) is None:
    exp_type = (session_info or {}).get(
      "experiment_type"
    )
    if exp_type:
      sm.data_characteristics[
        "experiment_type"
      ] = exp_type

  # --- Update from phaser (if in history) ---
  # Same pattern: only on first encounter.
  if sm.data_characteristics.get("mr_tfz") is None:
    for entry in history:
      prog = entry.get("program", "")
      if "phaser" not in prog:
        continue
      la = entry.get("analysis") or entry.get(
        "log_analysis", {}
      )
      if not la:
        continue
      phaser_data = {}
      for key in ("tfz", "llg", "space_group",
                   "unit_cell", "n_copies"):
        if key in la:
          phaser_data[key] = la[key]
      if phaser_data:
        sm.update_from_phaser(phaser_data)
        break

  # --- Update from autosol (if in history) ---
  if sm.data_characteristics.get(
    "phasing_fom"
  ) is None:
    for entry in history:
      prog = (entry.get("program") or "")
      if "autosol" not in prog:
        continue
      la = entry.get("analysis") or entry.get(
        "log_analysis", {}
      )
      if not la:
        continue
      autosol_data = {}
      for key in ("fom", "bayes_cc",
                   "sites_found"):
        if key in la:
          autosol_data[key] = la[key]
      if autosol_data:
        sm.update_from_autosol(autosol_data)
        break

  # --- Update from current cycle's validation ---
  sm.update_from_validation(
    validation_result, log_metrics,
    cycle_number, program_name,
  )

  # --- Record validation history snapshot ---
  space_group = sm.data_characteristics.get(
    "space_group"
  )
  vh.record(
    cycle_number, program_name,
    validation_result, log_metrics,
    space_group=space_group,
  )

  return (sm, vh)


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


def _extract_hypothesis_from_assessment(
  assessment, context, state
):
  """Extract a hypothesis from the LLM response.

  If the assessment dict contains a "hypothesis" key
  (from the LLM including it in its JSON response),
  create a Hypothesis object and add it to the
  StructureModel.

  Args:
    assessment: dict from parse_assessment.
    context: dict from _build_thinking_context.
    state: graph state dict (for cycle_number).

  Modifies the StructureModel in context in-place.
  Never raises.
  """
  hyp_text = assessment.get("hypothesis")
  if not hyp_text:
    return

  sm = context.get("structure_model")
  if sm is None:
    return

  try:
    if Hypothesis is None:
      return
    cycle = state.get("cycle_number", 0)

    h = Hypothesis(
      id="hyp_cycle%d" % cycle,
      statement=str(hyp_text)[:200],
      test_program=str(
        assessment.get("test_program", "")
      ),
      test_parameters=(
        assessment.get("test_parameters")
        if isinstance(
          assessment.get("test_parameters"), dict
        ) else {}
      ),
      confirm_if=str(
        assessment.get("confirm_if", "")
      ),
      refute_if=str(
        assessment.get("refute_if", "")
      ),
      status="proposed",
      proposed_at_cycle=cycle,
      test_cycles_remaining=int(
        assessment.get("test_cycles", 1)
      ),
    )

    added = sm.add_hypothesis(h)
    if added:
      import logging as _log_mod
      _log_mod.getLogger(__name__).debug(
        "Hypothesis extracted from LLM: %s",
        h.statement[:80],
      )
  except Exception:
    pass  # Non-critical — hypothesis is optional


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
