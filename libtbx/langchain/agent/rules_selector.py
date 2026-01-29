"""
Rules-Based Program Selector for PHENIX AI Agent.

This module provides deterministic program selection using only
workflow rules from YAML, without requiring an LLM.

The selector can:
1. Select the best program based on current state
2. Choose optimal files for each program
3. Apply appropriate strategy options
4. Provide clear reasoning for each decision

Usage:
    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()
    intent = selector.select_next_action(
        workflow_state=state,
        files=categorized_files,
        metrics_trend=trend,
        history_info=history_info
    )

This enables fully automated workflows without LLM involvement.
"""

from __future__ import absolute_import, division, print_function

import os

# Import dependencies
from libtbx.langchain.agent.workflow_engine import WorkflowEngine
from libtbx.langchain.agent.metric_evaluator import MetricEvaluator
from libtbx.langchain.agent.program_registry import ProgramRegistry


class RulesSelector:
    """
    Deterministic program selector using YAML-defined rules.

    Selects the next program based on:
    1. Current workflow phase
    2. Available files
    3. Metrics trend
    4. Quality targets
    """

    def __init__(self):
        """Initialize the selector with YAML-driven components."""
        self.workflow_engine = WorkflowEngine()
        self.metric_evaluator = MetricEvaluator()
        self.program_registry = ProgramRegistry()

    # =========================================================================
    # MAIN SELECTION
    # =========================================================================

    def select_next_action(self, workflow_state, files, metrics_trend=None,
                          history_info=None, user_advice=None, log_analysis=None):
        """
        Select the next action using rules only (no LLM).

        Args:
            workflow_state: Output from detect_workflow_state()
            files: Categorized files dict
            metrics_trend: Output from analyze_metrics_trend()
            history_info: Analyzed history dict
            user_advice: Optional user guidance string
            log_analysis: Optional log analysis dict with error info

        Returns:
            dict: Intent compatible with graph_nodes.py:
                {
                    program: str or None,
                    files: dict,
                    strategy: dict,
                    reasoning: str,
                    stop: bool,
                    stop_reason: str or None,
                    confidence: str,
                    selection_method: "rules",
                }
        """
        history_info = history_info or {}
        metrics_trend = metrics_trend or {}
        log_analysis = log_analysis or {}

        valid_programs = workflow_state.get("valid_programs", [])
        experiment_type = workflow_state.get("experiment_type", "xray")
        state_name = workflow_state.get("state", "unknown")
        resolution = workflow_state.get("resolution")

        # Extract strategy hints from user advice (e.g., resolution)
        user_strategy_hints = {}
        if user_advice:
            user_strategy_hints = self._extract_user_strategy_hints(user_advice)
            # User-specified resolution overrides detected resolution
            if "resolution" in user_strategy_hints:
                resolution = user_strategy_hints["resolution"]

        # 1. Check for stop conditions
        if self._should_stop(valid_programs, metrics_trend, workflow_state):
            return self._make_stop_intent(metrics_trend, workflow_state)

        # 2. Filter programs if user gave advice
        if user_advice:
            valid_programs = self._apply_user_advice(valid_programs, user_advice)

        # 3. Select best program from valid options
        program, reason = self._select_program(
            valid_programs,
            experiment_type,
            state_name,
            files,
            history_info,
            metrics_trend,
            workflow_state
        )

        if program is None:
            return self._make_stop_intent(
                {"reason": "No valid program available"},
                workflow_state
            )

        # 4. Select files for the program
        selected_files = self._select_files(program, files, workflow_state)

        # 5. Determine strategy options
        strategy = self._select_strategy(
            program,
            experiment_type,
            workflow_state,
            history_info,
            metrics_trend,
            resolution,
            log_analysis
        )

        # Apply any additional user strategy hints
        for key, value in user_strategy_hints.items():
            if key not in strategy or strategy[key] is None:
                strategy[key] = value

        # 6. Validate invariants and apply fixes (via TemplateBuilder)
        # This is optional here since build node also calls it, but ensures
        # the reasoning reflects any fixes applied
        from libtbx.langchain.agent.template_builder import TemplateBuilder
        builder = TemplateBuilder()
        invariant_context = {
            "resolution": resolution,
            "workflow_state": workflow_state,
        }
        selected_files, strategy, warnings = builder.validate_and_fix(
            program, selected_files, strategy,
            context=invariant_context
        )
        # warnings are logged if needed, but reasoning is built with fixed values

        # 7. Build reasoning
        full_reason = self._build_reasoning(
            program, reason, selected_files, strategy,
            workflow_state, metrics_trend
        )

        # Add any invariant warnings to reasoning
        if warnings:
            full_reason += " [Note: %s]" % "; ".join(warnings)

        return {
            "program": program,
            "files": selected_files,
            "strategy": strategy,
            "reasoning": full_reason,
            "stop": False,
            "stop_reason": None,
            "confidence": "high",
            "selection_method": "rules",
        }

    # =========================================================================
    # STOP CONDITIONS
    # =========================================================================

    def _should_stop(self, valid_programs, metrics_trend, workflow_state):
        """Check if we should stop."""
        # Explicit stop from metrics
        if metrics_trend.get("should_stop"):
            return True

        # Only STOP is valid
        if valid_programs == ["STOP"]:
            return True

        # Complete state
        if workflow_state.get("state") == "complete":
            return True

        return False

    def _make_stop_intent(self, metrics_trend, workflow_state):
        """Create a stop intent."""
        reason = metrics_trend.get("reason") or workflow_state.get("reason") or "Workflow complete"
        trend_summary = metrics_trend.get("trend_summary", "")

        return {
            "program": None,
            "files": {},
            "strategy": {},
            "reasoning": "Stopping: %s. %s" % (reason, trend_summary),
            "stop": True,
            "stop_reason": reason,
            "confidence": "high",
            "selection_method": "rules",
        }

    # =========================================================================
    # PROGRAM SELECTION
    # =========================================================================

    def _select_program(self, valid_programs, experiment_type, state_name,
                       files, history_info, metrics_trend, workflow_state):
        """
        Select the best program from valid options.

        Returns:
            tuple: (program_name, reason)
        """
        if not valid_programs:
            return None, "No valid programs"

        # Filter out STOP for now
        candidates = [p for p in valid_programs if p != "STOP"]

        if not candidates:
            return None, "Only STOP available"

        # Single candidate - easy choice
        if len(candidates) == 1:
            return candidates[0], "Only valid program for this state"

        # Multiple candidates - use priority rules
        return self._prioritize_programs(
            candidates, experiment_type, state_name,
            files, history_info, metrics_trend, workflow_state
        )

    def _prioritize_programs(self, candidates, experiment_type, state_name,
                            files, history_info, metrics_trend, workflow_state):
        """
        Prioritize among multiple valid programs.

        Priority rules:
        0. Forced program from after_program directive (highest priority)
        1. Programs with priority_when conditions met (from workflow YAML)
        2. Validation programs when metrics are good
        3. Refinement programs when model needs improvement
        4. Building programs when model is missing
        5. Analysis programs at start
        """
        # Check for forced program from directives (after_program)
        # The workflow_engine puts this at the front of valid_programs
        # We check if it's in candidates and select it immediately
        forced_program = workflow_state.get("forced_program")
        if forced_program and forced_program in candidates:
            return forced_program, "User directive: run %s before stopping" % forced_program

        # First check for programs with priority_when conditions met
        # These are defined in workflows.yaml and checked by workflow_engine
        program_priorities = workflow_state.get("program_priorities", [])
        for prog in program_priorities:
            if prog in candidates:
                return prog, "Priority condition met (from workflow)"

        # Define priority order by context
        if state_name in ("analyze", "xray_initial", "cryoem_initial"):
            # At start: prefer analysis
            priority = [
                "phenix.xtriage", "phenix.mtriage",
                "phenix.predict_and_build", "phenix.phaser",
            ]

        elif state_name in ("validate",):
            # Validation phase: prefer validation tools
            priority = [
                "phenix.molprobity", "phenix.model_vs_data",
                "phenix.refine", "phenix.real_space_refine",
            ]

        elif state_name in ("refine",):
            # Refinement: check if we should validate first
            r_free = metrics_trend.get("r_free_trend", [None])[-1] if metrics_trend.get("r_free_trend") else None
            map_cc = metrics_trend.get("map_cc_trend", [None])[-1] if metrics_trend.get("map_cc_trend") else None

            # If quality is good, prefer validation
            if r_free and r_free < 0.28:
                if "phenix.molprobity" in candidates:
                    return "phenix.molprobity", "Quality good (R-free=%.3f), validating" % r_free
            if map_cc and map_cc > 0.75:
                if "phenix.molprobity" in candidates:
                    return "phenix.molprobity", "Quality good (CC=%.3f), validating" % map_cc

            # Otherwise continue refinement
            priority = [
                "phenix.refine", "phenix.real_space_refine",
                "phenix.molprobity", "phenix.model_vs_data",
            ]

        elif state_name in ("obtain_model", "molecular_replacement"):
            # Model building priority
            # Note: strong_anomalous check is now handled by program_priorities above
            if files.get("sequence") and not files.get("pdb"):
                priority = [
                    "phenix.predict_and_build",
                    "phenix.phaser", "phenix.autobuild",
                ]
            else:
                priority = [
                    "phenix.phaser", "phenix.dock_in_map",
                    "phenix.predict_and_build", "phenix.autobuild",
                ]

        elif "ligand" in state_name.lower():
            # Ligand work
            priority = [
                "phenix.ligandfit", "phenix.elbow",
                "phenix.pdbtools", "phenix.refine",
            ]

        else:
            # Default priority
            priority = [
                "phenix.refine", "phenix.real_space_refine",
                "phenix.molprobity", "phenix.autobuild",
                "phenix.predict_and_build", "phenix.phaser",
            ]

        # Find first candidate in priority list
        for prog in priority:
            if prog in candidates:
                return prog, "Highest priority for %s state" % state_name

        # Fall back to first candidate
        return candidates[0], "First valid program"

    def _apply_user_advice(self, valid_programs, user_advice):
        """Filter programs based on user advice.

        Filters valid_programs based on user preferences. Does not add programs
        that aren't already valid - the workflow state machine controls validity.

        IMPORTANT: This function should only filter when the user gives a SPECIFIC
        IMMEDIATE instruction (e.g., "run refine now", "use phaser"). For multi-step
        workflows (e.g., "refine, fit ligand, then refine again"), the directives
        system handles the sequencing - we should NOT filter here.

        Keywords are now defined in programs.yaml under user_advice_keywords.
        """
        advice_lower = user_advice.lower()

        # Check for specific program mentions (e.g., "use phaser", "run refine")
        # BUT only if it looks like an immediate instruction, not a workflow description
        # Workflow descriptions contain sequencing words or multiple program mentions
        sequencing_words = [
            "then", "after", "first", "next", "finally", "before", "followed by",
            "with", "and then", "including", "include", "plus", "also",
            "workflow", "sequence", "steps", "primary goal", "goal:",
        ]
        is_multi_step = any(word in advice_lower for word in sequencing_words)

        # Also check if multiple programs are mentioned - that's a multi-step workflow
        programs_mentioned = 0
        for prog in valid_programs:
            if prog == "STOP":
                continue
            prog_name = prog.replace("phenix.", "").replace("_", " ").lower()
            if prog_name in advice_lower or prog.lower() in advice_lower:
                programs_mentioned += 1
        if programs_mentioned > 1:
            is_multi_step = True

        if not is_multi_step:
            # Single-step instruction - check for specific program mention
            for prog in valid_programs:
                prog_name = prog.replace("phenix.", "").replace("_", " ")
                # Skip STOP - we handle it specially later to avoid matching "stop" in "stop after..."
                if prog == "STOP":
                    continue
                if prog_name in advice_lower or prog.lower() in advice_lower:
                    return [prog]  # User specifically wants this

            # Check each valid program's keywords from YAML
            for prog in valid_programs:
                keywords = self.program_registry.get_user_advice_keywords(prog)
                if keywords:
                    for kw in keywords:
                        if kw.lower() in advice_lower:
                            return [prog]

        # NOTE: Removed aggressive generic keyword fallbacks (e.g., "refine" -> only refine programs)
        # These caused problems with multi-step workflows like "refine, fit ligand, refine"
        # The directives system and LLM handle complex advice better than simple keyword matching.

        # Only treat "stop" as an immediate stop request if it's at the start
        # or clearly indicates immediate stop (not "stop after X")
        # Words like "stop after", "stop when", "stop condition" are conditions, not requests
        if "stop" in advice_lower:
            # Check if it's a stop condition vs an immediate stop request
            stop_condition_patterns = [
                "stop after",
                "stop when",
                "stop once",
                "stop if",
                "stop condition",
                "stop at",
            ]
            is_stop_condition = any(pattern in advice_lower for pattern in stop_condition_patterns)

            if not is_stop_condition:
                # Looks like an immediate stop request (e.g., "please stop", "stop now")
                return ["STOP"] if "STOP" in valid_programs else valid_programs

        return valid_programs

    def _extract_user_strategy_hints(self, user_advice):
        """
        Extract strategy hints from user advice.

        Parses user advice for specific values like resolution.

        Args:
            user_advice: User-provided guidance string

        Returns:
            dict: Strategy hints to apply (e.g., {"resolution": 3.5})
        """
        import re
        hints = {}

        if not user_advice:
            return hints

        advice_lower = user_advice.lower()

        # Extract resolution
        # Patterns: "resolution=3.5", "resolution 3.5", "3.5 Angstrom", "3.5A resolution"
        resolution_patterns = [
            r'resolution\s*[=:]\s*([0-9.]+)',
            r'resolution\s+([0-9.]+)',
            r'([0-9.]+)\s*(?:angstrom|Å|A)\s*(?:resolution)?',
            r'(?:use|at|to)\s+([0-9.]+)\s*(?:angstrom|Å|A)?',
        ]

        for pattern in resolution_patterns:
            match = re.search(pattern, user_advice, re.IGNORECASE)
            if match:
                try:
                    res = float(match.group(1))
                    # Sanity check: resolution should be between 0.5 and 20
                    if 0.5 <= res <= 20:
                        hints["resolution"] = res
                        break
                except ValueError:
                    pass

        return hints

    # =========================================================================
    # FILE SELECTION
    # =========================================================================

    def _select_files(self, program, files, workflow_state=None):
        """
        Select appropriate files for a program.

        Selects both required and optional inputs when matching files are available.
        For predict_and_build in stepwise mode, excludes data files (MTZ/map).

        Args:
            program: Program name
            files: Categorized files dict
            workflow_state: Optional workflow state dict

        Returns:
            dict: {input_name: filepath, ...}
        """
        selected = {}
        workflow_state = workflow_state or {}
        automation_path = workflow_state.get("automation_path", "automated")
        experiment_type = workflow_state.get("experiment_type", "xray")

        # For predict_and_build in stepwise mode, skip data files
        # This triggers stop_after_predict=True behavior
        skip_data_files = (
            program == "phenix.predict_and_build" and
            automation_path == "stepwise"
        )
        data_inputs = {"data_mtz", "map_coeffs_mtz", "full_map", "half_map", "map", "data"}

        # Get required inputs from registry
        required = self.program_registry.get_required_inputs(program)

        for input_name in required:
            if skip_data_files and input_name in data_inputs:
                continue
            file_list = self._get_files_for_input(input_name, files, program, experiment_type)
            if file_list:
                # Select best file (usually most recent or highest cycle)
                selected[input_name] = self._select_best_file(file_list, input_name)

        # Also include optional inputs when we have matching files
        optional = self.program_registry.get_optional_inputs(program)

        for input_name in optional:
            if skip_data_files and input_name in data_inputs:
                continue
            if input_name not in selected:  # Don't override required
                file_list = self._get_files_for_input(input_name, files, program, experiment_type)
                if file_list:
                    selected[input_name] = self._select_best_file(file_list, input_name)

        return selected

    def _get_files_for_input(self, input_name, files, program=None, experiment_type=None):
        """Map input name to file category.

        Args:
            input_name: Name of the input (e.g., 'model', 'mtz')
            files: Categorized files dict
            program: Program name (for context-specific selection)
            experiment_type: 'xray' or 'cryoem'
        """
        # Check if program has YAML-defined input priorities
        priorities = self.program_registry.get_input_priorities(program, input_name)
        priority_categories = priorities.get("categories", [])
        exclude_categories = priorities.get("exclude_categories", [])

        # If we have YAML-defined priorities, use them
        if priority_categories:
            result_files = []
            for category in priority_categories:
                cat_files = files.get(category, [])
                for f in cat_files:
                    if f not in result_files:
                        # Check exclusions
                        excluded = False
                        for excl_cat in exclude_categories:
                            if f in files.get(excl_cat, []):
                                excluded = True
                                break
                        if not excluded:
                            result_files.append(f)

            # Also add any files from pdb category that aren't in excluded categories
            if input_name == "model" and "pdb" in priority_categories:
                for f in files.get("pdb", []):
                    if f not in result_files:
                        basename = os.path.basename(f).lower()
                        excluded = False
                        for excl_cat in exclude_categories:
                            # Check by category and by pattern
                            if f in files.get(excl_cat, []):
                                excluded = True
                                break
                            if excl_cat == "predicted" and 'predicted' in basename and 'processed' not in basename:
                                excluded = True
                                break
                        if not excluded:
                            result_files.append(f)

            if result_files:
                return result_files

        # Fallback: simple default mappings for programs without input_priorities
        # Most programs should have input_priorities defined in programs.yaml

        # Default model priority (used if no YAML priorities defined)
        if input_name == "model":
            # Generic fallback: prefer refined/placed models over raw predictions
            model_files = []
            for cat in ["with_ligand", "refined", "phaser_output", "rsr_output",
                        "docked", "processed_predicted", "autobuild_output"]:
                model_files.extend(files.get(cat, []))
            # Add generic PDBs, excluding raw predictions
            for f in files.get("pdb", []):
                basename = os.path.basename(f).lower()
                if 'predicted' not in basename or 'processed' in basename:
                    if f not in model_files:
                        model_files.append(f)
            return model_files

        # Default data_mtz priority (for refinement)
        if input_name in ("data_mtz", "data"):
            mtz_files = list(files.get("original_data_mtz", []))
            for f in files.get("data_mtz", []):
                if f not in mtz_files:
                    mtz_files.append(f)
            return mtz_files

        # Map coefficients MTZ priority (for ligand fitting)
        if input_name == "map_coeffs_mtz":
            coeffs_files = list(files.get("refine_map_coeffs", []))
            for f in files.get("denmod_map_coeffs", []):
                if f not in coeffs_files:
                    coeffs_files.append(f)
            for f in files.get("map_coeffs_mtz", []):
                if f not in coeffs_files:
                    coeffs_files.append(f)
            return coeffs_files

        # Simple mappings for other input types
        mapping = {
            "sequence": files.get("sequence", []),
            "map": files.get("full_map", []) + files.get("map", []),
            "full_map": files.get("full_map", []),
            "half_map": files.get("half_map", []),
            "half_map_1": files.get("half_map", [])[:1] if files.get("half_map") else [],
            "half_map_2": files.get("half_map", [])[1:2] if len(files.get("half_map", [])) > 1 else [],
            "ligand_cif": files.get("ligand_cif", []),
            "ligand_pdb": files.get("ligand_pdb", []),
        }

        return mapping.get(input_name, [])

    def _select_best_file(self, file_list, input_type):
        """Select the best file from a list."""
        if not file_list:
            return None

        if len(file_list) == 1:
            return file_list[0]

        # For models, prefer by category then by cycle number
        if input_type == "model":
            # The list is already ordered by priority from _get_files_for_input
            # Within each category, prefer higher cycle numbers

            import re

            def get_model_score(f):
                """Score model files: higher = better."""
                basename = os.path.basename(f).lower()

                # Category scores (higher = better)
                if 'with_ligand' in basename:
                    category_score = 7000
                elif '_refine_' in basename or 'refined' in basename:
                    category_score = 6000
                elif 'phaser' in basename or '.1.pdb' in basename:
                    category_score = 5000
                elif 'real_space_refined' in basename or '_rsr_' in basename:
                    category_score = 4000
                elif 'dock' in basename:
                    category_score = 3000
                elif 'processed' in basename:
                    category_score = 2000
                elif 'autobuild' in basename:
                    category_score = 1500
                else:
                    category_score = 0

                # Extract cycle number for tie-breaking within category
                match = re.search(r'_(\d+)\.pdb', f)
                cycle = int(match.group(1)) if match else 0

                return category_score + cycle

            return max(file_list, key=get_model_score)

        # Default: return first
        return file_list[0]

    # =========================================================================
    # STRATEGY SELECTION
    # =========================================================================

    def _select_strategy(self, program, experiment_type, workflow_state,
                        history_info, metrics_trend, resolution, log_analysis=None):
        """
        Select strategy options for a program.

        Note: This provides basic strategy selection. Complex error recovery
        (like R-free flag generation) is handled by the LLM when available.

        Returns:
            dict: Strategy options to pass to program
        """
        strategy = {}

        # Resolution (round to 1 decimal place)
        if resolution:
            strategy["resolution"] = round(float(resolution), 1)

        # Get cycle number for output prefix
        cycle_number = history_info.get("cycle_number", 1)

        # Program-specific strategy
        if program == "phenix.predict_and_build":
            # Check if stepwise mode
            if workflow_state.get("automation_path") == "stepwise":
                strategy["stop_after_predict"] = True

        elif program == "phenix.refine":
            # Use descriptive prefix: refine_001, refine_002, etc.
            strategy["output_prefix"] = "refine_%03d" % cycle_number

            # Check for twinning
            if history_info.get("twin_law"):
                strategy["twin_law"] = history_info["twin_law"]

            # Check if ncs is present
            if history_info.get("has_ncs"):
                strategy["ncs_constraints"] = True

            # Add ordered solvent in later cycles
            refine_count = history_info.get("refine_count", 0)
            if refine_count >= 2:
                strategy["ordered_solvent"] = True

        elif program == "phenix.real_space_refine":
            # Use descriptive prefix: rsr_001, rsr_002, etc.
            strategy["output_prefix"] = "rsr_%03d" % cycle_number

            # Run longer in later cycles
            rsr_count = history_info.get("rsr_count", 0)
            if rsr_count >= 2:
                strategy["macro_cycles"] = 5

        elif program == "phenix.autobuild":
            # Quick build mode if we have good phases
            tfz = metrics_trend.get("tfz") if metrics_trend else None
            if tfz and tfz > 12:
                strategy["quick"] = True

        return strategy

    # =========================================================================
    # REASONING
    # =========================================================================

    def _build_reasoning(self, program, select_reason, files, strategy,
                        workflow_state, metrics_trend):
        """Build human-readable reasoning string."""
        parts = []

        # State context
        state = workflow_state.get("state", "unknown")
        exp_type = workflow_state.get("experiment_type", "unknown")
        parts.append("State: %s (%s)" % (state, exp_type))

        # Selection reason
        parts.append("Selected %s: %s" % (program, select_reason))

        # Files
        if files:
            file_list = ["%s=%s" % (k, os.path.basename(v)) for k, v in files.items()]
            parts.append("Files: %s" % ", ".join(file_list))

        # Strategy
        if strategy:
            strat_list = ["%s=%s" % (k, v) for k, v in strategy.items()]
            parts.append("Strategy: %s" % ", ".join(strat_list))

        # Metrics context
        if metrics_trend:
            summary = metrics_trend.get("trend_summary", "")
            if summary:
                parts.append("Metrics: %s" % summary)

        return ". ".join(parts)


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

_selector = None

def get_selector():
    """Get global RulesSelector instance."""
    global _selector
    if _selector is None:
        _selector = RulesSelector()
    return _selector


def select_action_by_rules(workflow_state, files, metrics_trend=None,
                           history_info=None, user_advice=None, log_analysis=None):
    """Convenience function for rules-based selection."""
    return get_selector().select_next_action(
        workflow_state, files, metrics_trend, history_info,
        user_advice=user_advice, log_analysis=log_analysis
    )


# =============================================================================
# TESTING
# =============================================================================

if __name__ == "__main__":
    print("Testing RulesSelector...")
    print()

    selector = RulesSelector()

    # Test X-ray initial state
    workflow_state = {
        "state": "analyze",
        "experiment_type": "xray",
        "valid_programs": ["phenix.xtriage"],
        "resolution": 2.0,
    }
    files = {"data_mtz": ["data.mtz"], "sequence": ["seq.fa"]}

    intent = selector.select_next_action(workflow_state, files)
    print("X-ray initial:")
    print("  Program:", intent["program"])
    print("  Reasoning:", intent["reasoning"][:80] + "...")
    print()

    # Test X-ray refinement state
    workflow_state = {
        "state": "refine",
        "experiment_type": "xray",
        "valid_programs": ["phenix.refine", "phenix.molprobity", "STOP"],
        "resolution": 2.0,
    }
    files = {"data_mtz": ["data.mtz"], "pdb": ["model.pdb"], "refined": ["model_refine_001.pdb"]}
    metrics_trend = {
        "r_free_trend": [0.30, 0.28, 0.26],
        "improvement_rate": 7.0,
        "should_stop": False,
    }
    history_info = {"refine_count": 3}

    intent = selector.select_next_action(workflow_state, files, metrics_trend, history_info)
    print("X-ray refining:")
    print("  Program:", intent["program"])
    print("  Strategy:", intent["strategy"])
    print("  Reasoning:", intent["reasoning"][:80] + "...")
    print()

    # Test cryo-EM state
    workflow_state = {
        "state": "refine",
        "experiment_type": "cryoem",
        "valid_programs": ["phenix.real_space_refine", "phenix.molprobity"],
        "resolution": 3.0,
    }
    files = {"full_map": ["map.mrc"], "pdb": ["model.pdb"]}

    intent = selector.select_next_action(workflow_state, files)
    print("Cryo-EM refining:")
    print("  Program:", intent["program"])
    print("  Reasoning:", intent["reasoning"][:80] + "...")
    print()

    # Test stop condition
    workflow_state = {
        "state": "complete",
        "valid_programs": ["STOP"],
    }
    metrics_trend = {"should_stop": True, "reason": "Target reached"}

    intent = selector.select_next_action(workflow_state, {}, metrics_trend)
    print("Stop condition:")
    print("  Stop:", intent["stop"])
    print("  Reason:", intent["stop_reason"])
