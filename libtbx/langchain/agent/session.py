"""
Agent session management - persistent tracking of agent runs.

This module handles:
- Tracking cycle-by-cycle progress
- Persisting session state to JSON
- Checking for duplicate commands
- Generating session summaries

Usage:
    from libtbx.langchain.agent import AgentSession

    session = AgentSession(session_dir='./agent_session')
    session.start_cycle(cycle_number=1)
    session.record_decision(program='phenix.xtriage', decision='...', ...)
    session.record_result(result='Success', output_files=[...])
    session.save()
"""
from __future__ import absolute_import, division, print_function

import os
import re
import json
from datetime import datetime
from langchain_core.prompts import PromptTemplate

from libtbx.langchain.agent.best_files_tracker import BestFilesTracker


class AgentSession:
    """
    Manages persistent state for an agent session across multiple cycles.

    The session is stored as a JSON file and tracks:
    - Project info (advice, original files)
    - Each cycle's decision, command, and result
    - Output files from each cycle
    - A final LLM-generated summary
    """

    def __init__(self, session_dir=None, session_file=None):
        """
        Initialize or load an agent session.

        Args:
            session_dir: Directory to store session file
            session_file: Explicit path to session file (overrides session_dir)
        """
        if session_file:
            self.session_file = session_file
        elif session_dir:
            if not os.path.exists(session_dir):
                os.makedirs(session_dir)
            self.session_file = os.path.join(session_dir, "agent_session.json")
        else:
            self.session_file = "agent_session.json"

        # Load existing session or create new one
        if os.path.exists(self.session_file):
            self.load()
        else:
            self._init_new_session()

    def _init_new_session(self):
        """Initialize a new session."""
        self.data = {
            "session_id": datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
            "project_advice": "",
            "original_files": [],
            "cycles": [],
            "summary": "",
            "resolution": None,  # Single source of truth for resolution
            "resolution_source": None,  # Which program/cycle set it (e.g., "mtriage_1")
            "ignored_commands": {},  # Commands that failed due to input errors
            # Red flag detection settings
            "experiment_type": None,  # "xray" or "cryoem", locked after first detection
            "abort_on_red_flags": True,  # Abort on critical sanity check failures
            "abort_on_warnings": False,  # Also abort on warning-level issues
            # R-free MTZ tracking (X-ray only)
            # Once R-free flags are generated, we must use the same MTZ for all refinements
            "rfree_mtz": None,  # Path to the locked R-free MTZ
            "rfree_mtz_locked_at_cycle": None,  # Cycle when it was locked
        }
        # Initialize best files tracker (not stored in self.data, serialized separately)
        self.best_files = BestFilesTracker()

    def load(self):
        """Load session from file."""
        try:
            with open(self.session_file, 'r') as f:
                self.data = json.load(f)
            # Restore best files tracker
            best_files_data = self.data.get("best_files")
            self.best_files = BestFilesTracker.from_dict(best_files_data)
        except Exception as e:
            print(f"Warning: Could not load session file: {e}")
            self._init_new_session()

    def save(self):
        """Save session to file."""
        try:
            # Include best files tracker in saved data
            self.data["best_files"] = self.best_files.to_dict()
            self.data["best_files_history"] = [h.to_dict() for h in self.best_files.get_history()]
            with open(self.session_file, 'w') as f:
                json.dump(self.data, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save session file: {e}")

    def set_project_info(self, project_advice=None, original_files=None):
        """Set project-level information."""
        if project_advice:
            self.data["project_advice"] = project_advice
        if original_files:
            if isinstance(original_files, str):
                self.data["original_files"] = original_files.split()
            else:
                self.data["original_files"] = list(original_files)
        self.save()

    def get_resolution(self):
        """
        Get the session resolution (single source of truth).

        Returns:
            float or None: Resolution in Angstroms
        """
        return self.data.get("resolution")

    def set_resolution(self, resolution, source, force=False):
        """
        Set the session resolution.

        Args:
            resolution: Resolution value in Angstroms
            source: String describing source (e.g., "mtriage_1", "resolve_cryo_em_3")
            force: If True, always update. If False, only update if better (lower) or not set.

        Returns:
            tuple: (was_updated, message)
        """
        current = self.data.get("resolution")

        if current is None:
            # First time setting resolution
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution set to {resolution:.2f}Å from {source}"

        if force:
            old = current
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution updated from {old:.2f}Å to {resolution:.2f}Å (forced by {source})"

        # Only update if better (lower resolution value = higher resolution data)
        if resolution < current:
            old = current
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution improved from {old:.2f}Å to {resolution:.2f}Å by {source}"
        elif resolution > current + 0.5:
            # Significantly worse - warn but don't update
            return False, f"WARNING: {source} reports {resolution:.2f}Å but session resolution is {current:.2f}Å (not updating)"
        else:
            # Similar value, no update needed
            return False, None

    # =========================================================================
    # EXPERIMENT TYPE (locked after first detection)
    # =========================================================================

    def get_experiment_type(self):
        """
        Get the locked experiment type.

        Returns:
            str or None: "xray" or "cryoem", or None if not yet set
        """
        return self.data.get("experiment_type")

    def set_experiment_type(self, exp_type):
        """
        Set the experiment type (only if not already set).

        The experiment type is locked after first detection to prevent
        spurious files from changing the workflow mid-session.

        Args:
            exp_type: "xray" or "cryoem"

        Returns:
            tuple: (was_set, message)
        """
        current = self.data.get("experiment_type")

        if current is None:
            self.data["experiment_type"] = exp_type
            self.save()
            return True, f"Experiment type set to {exp_type}"

        if current == exp_type:
            return False, None  # Same type, no change needed

        # Type mismatch - this is a red flag situation
        return False, f"Experiment type already locked as {current} (attempted to set {exp_type})"

    # =========================================================================
    # R-FREE MTZ TRACKING (X-ray only)
    # =========================================================================
    # In X-ray crystallography, R-free flags are used for cross-validation.
    # Once generated, the SAME flags must be used for all subsequent refinements.
    # Using different R-free flags invalidates the R-free statistic.
    #
    # This tracking ensures:
    # 1. First refinement that generates R-free flags locks the MTZ
    # 2. All subsequent refinements MUST use this locked MTZ
    # 3. The locked MTZ is visible in session.json for debugging

    def get_rfree_mtz(self):
        """
        Get the locked R-free MTZ path.

        Returns:
            str or None: Path to the MTZ with R-free flags, or None if not yet set
        """
        return self.data.get("rfree_mtz")

    def set_rfree_mtz(self, path, cycle):
        """
        Lock the R-free MTZ. Can only be set once per session.

        This should be called when phenix.refine first generates or uses
        R-free flags. Once locked, all subsequent refinements must use
        this MTZ to maintain valid R-free statistics.

        Args:
            path: Full path to the MTZ file with R-free flags
            cycle: Cycle number when the lock occurred

        Returns:
            tuple: (was_locked, message)
        """
        current = self.data.get("rfree_mtz")

        if current is not None:
            # Already locked - this is expected, not an error
            if current == path:
                return False, None  # Same path, no change
            return False, f"R-free MTZ already locked to {current} (cycle {self.data.get('rfree_mtz_locked_at_cycle')})"

        # Lock it
        self.data["rfree_mtz"] = path
        self.data["rfree_mtz_locked_at_cycle"] = cycle
        self.save()
        return True, f"R-free MTZ locked to {os.path.basename(path)} at cycle {cycle}"

    def is_rfree_mtz_locked(self):
        """Check if an R-free MTZ has been locked."""
        return self.data.get("rfree_mtz") is not None

    # =========================================================================
    # SANITY CHECK SETTINGS
    # =========================================================================

    def get_abort_settings(self):
        """
        Get the abort settings for sanity checking.

        Returns:
            dict: {"abort_on_red_flags": bool, "abort_on_warnings": bool}
        """
        return {
            "abort_on_red_flags": self.data.get("abort_on_red_flags", True),
            "abort_on_warnings": self.data.get("abort_on_warnings", False),
        }

    def set_abort_settings(self, abort_on_red_flags=None, abort_on_warnings=None):
        """
        Set the abort settings for sanity checking.

        Args:
            abort_on_red_flags: If True, abort on critical issues
            abort_on_warnings: If True, also abort on warning-level issues
        """
        if abort_on_red_flags is not None:
            self.data["abort_on_red_flags"] = abort_on_red_flags
        if abort_on_warnings is not None:
            self.data["abort_on_warnings"] = abort_on_warnings
        self.save()

    # =========================================================================
    # IGNORED COMMANDS (Input errors that should not be repeated)
    # =========================================================================

    def record_ignored_command(self, command, reason):
        """
        Record a command that failed due to input error (wrong parameters, etc).

        These commands should not be repeated - they indicate agent mistakes
        that need to be avoided on retry.

        Args:
            command: The full command string that failed
            reason: Brief description of why it failed
        """
        if "ignored_commands" not in self.data:
            self.data["ignored_commands"] = {}
        self.data["ignored_commands"][command] = {
            "reason": reason,
            "timestamp": datetime.now().isoformat()
        }
        self.save()

    def get_ignored_commands(self):
        """
        Get dict of commands that failed due to input errors.

        Returns:
            dict: {command: {"reason": str, "timestamp": str}, ...}
        """
        return self.data.get("ignored_commands", {})

    def is_command_ignored(self, command):
        """
        Check if a command was previously ignored.

        Args:
            command: Command string to check

        Returns:
            tuple: (is_ignored, reason) - reason is None if not ignored
        """
        ignored = self.data.get("ignored_commands", {})
        if command in ignored:
            return True, ignored[command].get("reason", "Unknown")
        return False, None

    def clear_ignored_commands(self):
        """Clear all ignored commands (use when starting fresh)."""
        self.data["ignored_commands"] = {}
        self.save()

    # =========================================================================
    # BEST FILES TRACKING
    # =========================================================================

    def get_best_model(self):
        """
        Get the path to the best model file.

        Returns:
            str or None: Path to best model
        """
        return self.best_files.get_best_path("model")

    def get_best_map(self):
        """
        Get the path to the best map file.

        Returns:
            str or None: Path to best map
        """
        return self.best_files.get_best_path("map")

    def get_best_mtz(self):
        """
        Get the path to the best MTZ file.

        Returns:
            str or None: Path to best MTZ (with R-free flags if available)
        """
        return self.best_files.get_best_path("mtz")

    def get_best_files_dict(self):
        """
        Get all best files as a dict of category -> path.

        Returns:
            dict: {category: path, ...}
        """
        return self.best_files.get_best_dict()

    def get_best_files_summary(self):
        """
        Get a human-readable summary of best files.

        Returns:
            str: Multi-line summary
        """
        return self.best_files.get_summary()

    def evaluate_file_for_best(self, path, cycle, metrics=None, stage=None):
        """
        Manually evaluate a file for the best files tracker.

        Useful for evaluating user-provided input files at session start.

        Args:
            path: File path
            cycle: Cycle number (use 0 for input files)
            metrics: Optional metrics dict
            stage: Optional stage override

        Returns:
            bool: True if file became new best
        """
        return self.best_files.evaluate_file(path, cycle, metrics, stage)

    def update_best_model_metrics(self, metrics, cycle=None):
        """
        Update metrics for the best model and recalculate its score.

        This should be called when validation provides new metrics for
        the current best model (e.g., after molprobity validation).

        Args:
            metrics: Dict with model metrics (r_free, map_cc, clashscore, etc.)
            cycle: Optional cycle number for logging

        Returns:
            tuple: (was_updated, old_score, new_score)
        """
        result = self.best_files.update_best_model_metrics(metrics, cycle)
        if result[0]:  # was_updated
            self.save()
        return result

    def write_history_file(self, output_path=None):
        """
        Write session history as a JSON file for the agent to read.

        Writes the full cycle data so the graph can properly analyze
        what programs have been run.

        Args:
            output_path: Path to write to. If None, creates temp file.

        Returns:
            str: Path to the written file
        """
        import tempfile

        if output_path is None:
            fd, output_path = tempfile.mkstemp(suffix='.json', prefix='session_history_')
            os.close(fd)

        # Write the full session data including all cycles
        # The graph's history loader will extract the cycles
        history_data = {
            "session_id": self.data.get("session_id", ""),
            "experiment_type": self.data.get("experiment_type"),
            "cycles": self.data.get("cycles", []),
        }

        with open(output_path, 'w') as f:
            json.dump(history_data, f, indent=2)

        return output_path

    def start_cycle(self, cycle_number):
        """
        Start a new cycle. Creates the cycle entry if it doesn't exist.

        Args:
            cycle_number: The cycle number (1-indexed)
        """
        # Ensure we have enough cycle entries
        while len(self.data["cycles"]) < cycle_number:
            self.data["cycles"].append({
                "cycle_number": len(self.data["cycles"]) + 1,
                "program": "",
                "decision": "",
                "reasoning": "",
                "explanation": "",
                "command": "",
                "result": "",
                "output_files": [],  # NEW: Track output files
                "timestamp": ""
            })
        self.save()

    def record_decision(self, cycle_number, program=None, decision=None,
                        reasoning=None, explanation=None, command=None):
        """
        Record the decision made for a cycle.

        Args:
            cycle_number: Which cycle to update
            program: The selected program
            decision: Short decision text
            reasoning: Detailed reasoning
            explanation: Full explanation including plan
            command: The generated command
        """
        self.start_cycle(cycle_number)
        cycle = self.data["cycles"][cycle_number - 1]

        if program is not None:
            cycle["program"] = program
        if decision is not None:
            cycle["decision"] = decision
        if reasoning is not None:
            cycle["reasoning"] = reasoning
        if explanation is not None:
            cycle["explanation"] = explanation
        if command is not None:
            cycle["command"] = command

        cycle["timestamp"] = datetime.now().isoformat()
        self.save()

    def record_result(self, cycle_number, result, output_files=None):
        """
        Record the result of running a cycle's command.

        Args:
            cycle_number: Which cycle to update
            result: Result text (success/failure description)
            output_files: List of output file paths produced by the command
        """
        self.start_cycle(cycle_number)
        cycle = self.data["cycles"][cycle_number - 1]
        cycle["result"] = result

        # Store output files
        if output_files:
            # Ensure output_files is a list
            if isinstance(output_files, str):
                output_files = [output_files]
            cycle["output_files"] = list(output_files)

        # Extract and store metrics from result for use by workflow_state
        program = cycle.get("program", "")
        metrics = self._extract_metrics_from_result(result, program)
        if metrics:
            cycle["metrics"] = metrics

            # Update session resolution from key programs
            if metrics.get("resolution"):
                resolution = metrics["resolution"]
                program_lower = program.lower()
                source = f"{program}_{cycle_number}"

                # Programs that establish resolution (first time)
                if "xtriage" in program_lower or "mtriage" in program_lower:
                    updated, msg = self.set_resolution(resolution, source)
                    if msg:
                        cycle["resolution_note"] = msg

                # Programs that might improve resolution (cryo-EM map optimization)
                elif "resolve_cryo_em" in program_lower or "map_sharpening" in program_lower:
                    # These can improve resolution - update if better
                    updated, msg = self.set_resolution(resolution, source)
                    if msg:
                        cycle["resolution_note"] = msg

                # Other programs - just note if significantly different
                else:
                    current = self.get_resolution()
                    if current and abs(resolution - current) > 0.5:
                        cycle["resolution_note"] = f"Program reports {resolution:.2f}Å (session: {current:.2f}Å)"

        # Update best files tracker with output files from successful cycles
        if output_files and "FAILED" not in result.upper():
            model_stage = self._infer_stage_from_program(program)

            # Check if this is an X-ray refinement program (produces MTZ with R-free flags)
            is_xray_refinement = (
                program and
                "refine" in program.lower() and
                "real_space" not in program.lower()  # real_space_refine is cryo-EM
            )

            for f in output_files:
                # Build file-specific metrics
                file_metrics = dict(metrics) if metrics else {}

                # Determine stage based on file type
                # MTZ files from refinement get "refined_mtz" stage
                # PDB/CIF files get the model stage (e.g., "refined")
                if f.lower().endswith('.mtz'):
                    basename = os.path.basename(f).lower()

                    # Check if this MTZ has R-free flags (from refinement or autosol)
                    # Patterns that indicate R-free flags are present:
                    # - *refine* (from phenix.refine)
                    # - *refinement_data* (from autosol)
                    has_rfree_patterns = (
                        'refine' in basename or
                        'refinement_data' in basename
                    )

                    # Patterns that indicate this is NOT suitable for refinement
                    # (phased MTZ files have multiple arrays)
                    is_phased_mtz = (
                        'phased' in basename or
                        'phases' in basename or
                        'solve' in basename
                    )

                    if has_rfree_patterns and not is_phased_mtz:
                        file_stage = "refined_mtz"
                        file_metrics["has_rfree_flags"] = True

                        # Lock the R-free MTZ (first one wins)
                        was_locked, lock_msg = self.set_rfree_mtz(f, cycle_number)
                        if was_locked:
                            cycle["rfree_mtz_locked"] = f
                            if lock_msg:
                                cycle["rfree_mtz_note"] = lock_msg
                    elif is_phased_mtz:
                        # Phased MTZ - not for direct refinement
                        file_stage = "phased_mtz"
                    else:
                        # Generic MTZ
                        file_stage = "mtz"
                else:
                    file_stage = model_stage  # Use inferred stage for models

                # Evaluate each output file - tracker will classify and score
                updated = self.best_files.evaluate_file(
                    path=f,
                    cycle=cycle_number,
                    metrics=file_metrics,
                    stage=file_stage
                )
                if updated:
                    category = self.best_files._classify_category(f)
                    if category:
                        cycle.setdefault("best_file_updates", []).append(category)

        # If this was a validation program, update best model metrics
        # (validation metrics apply to the current best model, not a new file)
        if metrics and program:
            program_lower = program.lower()
            validation_programs = ["molprobity", "validation", "clashscore"]
            is_validation = any(vp in program_lower for vp in validation_programs)

            if is_validation:
                # Extract model-relevant metrics
                model_metrics = {}
                for key in ["clashscore", "rama_favored", "rama_outliers",
                           "rotamer_outliers", "cbeta_outliers"]:
                    if key in metrics:
                        model_metrics[key] = metrics[key]

                if model_metrics:
                    was_updated, old_score, new_score = self.best_files.update_best_model_metrics(
                        model_metrics, cycle_number
                    )
                    if was_updated:
                        cycle["best_model_score_update"] = {
                            "old": old_score,
                            "new": new_score,
                            "reason": "validation metrics"
                        }

        self.save()

    def _infer_stage_from_program(self, program):
        """
        Infer output file stage from program name.

        Args:
            program: Program name (e.g., "phenix.refine", "phenix.real_space_refine")

        Returns:
            str: Stage name for best files tracker
        """
        if not program:
            return None

        program_lower = program.lower()

        # Model stages
        if "refine" in program_lower and "real_space" not in program_lower:
            return "refined"
        if "real_space_refine" in program_lower:
            return "rsr_output"
        if "autobuild" in program_lower:
            return "autobuild_output"
        if "dock_in_map" in program_lower:
            return "docked"
        if "process_predicted_model" in program_lower:
            return "processed_predicted"
        if "phaser" in program_lower:
            return "phaser_output"

        # Map stages
        if "resolve_cryo_em" in program_lower:
            return "optimized_full_map"
        if "map_sharpening" in program_lower or "local_aniso" in program_lower:
            return "sharpened"
        if "density_modification" in program_lower or "dm" in program_lower:
            return "density_modified"

        # MTZ stages
        if "xtriage" in program_lower:
            return "original"  # xtriage doesn't modify data

        return None

    def get_cycle(self, cycle_number):
        """Get data for a specific cycle."""
        if cycle_number <= len(self.data["cycles"]):
            return self.data["cycles"][cycle_number - 1]
        return None

    def get_all_commands(self):
        """
        Get all commands that have been run successfully (for duplicate detection).
        Only includes cycles where the command actually executed successfully.

        Returns:
            list of tuples: [(cycle_number, normalized_command), ...]
        """
        commands = []
        for cycle in self.data["cycles"]:
            # Skip CRASH cycles
            program = cycle.get("program", "")
            if program == "CRASH" or program == "":
                continue

            # Only include successful runs
            result = cycle.get("result", "")
            if not result.startswith("SUCCESS"):
                continue

            cmd = cycle.get("command", "")
            if cmd:
                # Normalize command (remove extra whitespace)
                norm_cmd = " ".join(cmd.strip().split())
                commands.append((cycle.get("cycle_number", 0), norm_cmd))

        return commands

    def cleanup_failed_cycles(self):
        """
        Remove failed/incomplete cycles from session.
        Called at session start to clean up previous crashed runs.

        Returns:
            int: Number of cycles removed
        """
        original_count = len(self.data["cycles"])

        # Keep only cycles that have both a command and a successful result
        self.data["cycles"] = [
            c for c in self.data["cycles"]
            if c.get("program") not in ("CRASH", "")
            and c.get("result", "").startswith("SUCCESS")
        ]

        # Renumber remaining cycles
        for i, cycle in enumerate(self.data["cycles"]):
            cycle["cycle_number"] = i + 1

        removed = original_count - len(self.data["cycles"])
        if removed > 0:
            print(f"Cleaned up {removed} failed/incomplete cycles from previous session")
            self.save()

        return removed

    def is_duplicate_command(self, command):
        """
        Check if a command (or very similar) has already been run successfully.

        Args:
            command: The command to check

        Returns:
            tuple: (is_duplicate: bool, previous_cycle: int or None)
        """
        if not command or command == "No command generated.":
            return False, None

        norm_new = " ".join(command.strip().split())

        for cycle_num, norm_cmd in self.get_all_commands():
            # Skip if this cycle had no real command
            if not norm_cmd:
                continue

            # Exact match
            if norm_cmd == norm_new:
                return True, cycle_num

            # Check if same program with same core parameters
            # (ignore path differences)
            new_parts = set(os.path.basename(p) for p in norm_new.split())
            old_parts = set(os.path.basename(p) for p in norm_cmd.split())

            # If >80% overlap in tokens, consider it a duplicate
            if len(new_parts) > 0 and len(old_parts) > 0:
                overlap = len(new_parts & old_parts) / max(len(new_parts), len(old_parts))
                if overlap > 0.8:
                    return True, cycle_num

        return False, None

    def get_history_for_agent(self):
        """
        Get cycle history in format expected by the agent.

        IMPORTANT: This now includes output_files so the agent knows
        what files were produced by each cycle, and metrics (as 'analysis')
        so workflow_state can access resolution, r_free, etc.

        Returns:
            list: List of dicts compatible with agent history format
        """
        history = []
        for cycle in self.data["cycles"]:
            if cycle.get("command") and cycle.get("result"):
                history_entry = {
                    "cycle_number": cycle.get("cycle_number", 0),
                    "job_id": str(cycle.get("cycle_number", "")),
                    "program": cycle.get("program", ""),
                    "command": cycle.get("command", ""),
                    "result": cycle.get("result", ""),
                    "summary": f"Decision: {cycle.get('decision', '')}\n"
                               f"Command: {cycle.get('command', '')}\n"
                               f"Result: {cycle.get('result', '')}",
                    # FIXED: Return metrics dict (not reasoning string) so workflow_state
                    # can access resolution, r_free, map_cc, etc.
                    "analysis": cycle.get("metrics", {}),
                    # NEW: Include output files so agent knows what was produced
                    "output_files": cycle.get("output_files", []),
                    "next_move": {
                        "command": cycle.get("command", ""),
                        "program": cycle.get("program", "")
                    }
                }
                history.append(history_entry)
        return history

    def get_all_output_files(self):
        """
        Get all output files from all successful cycles.

        Returns:
            list: List of file paths
        """
        all_files = []
        for cycle in self.data["cycles"]:
            if cycle.get("result", "").startswith("SUCCESS"):
                output_files = cycle.get("output_files", [])
                for f in output_files:
                    if f and f not in all_files:
                        all_files.append(f)
        return all_files

    def get_num_cycles(self):
        """Get the number of cycles recorded."""
        return len(self.data["cycles"])

    def format_cycle_summary(self, cycle_number):
        """
        Format a single cycle's information for display.

        Args:
            cycle_number: Which cycle to format

        Returns:
            str: Formatted cycle summary
        """
        cycle = self.get_cycle(cycle_number)
        if not cycle:
            return f"Cycle {cycle_number}: No data"

        lines = [
            f"",
            f"{'='*60}",
            f"CYCLE {cycle_number}",
            f"{'='*60}",
            f"Program: {cycle.get('program', 'N/A')}",
            f"Decision: {cycle.get('decision', 'N/A')}",
        ]

        # Truncate long reasoning
        reasoning = cycle.get('reasoning', 'N/A')
        if len(reasoning) > 300:
            reasoning = reasoning[:300] + "..."
        lines.append(f"Reasoning: {reasoning}")

        lines.extend([
            f"Command: {cycle.get('command', 'N/A')}",
            f"Result: {cycle.get('result', 'N/A')}",
        ])

        # Show output files if any
        output_files = cycle.get('output_files', [])
        if output_files:
            lines.append(f"Output files: {len(output_files)} files")
            for f in output_files[:5]:  # Show first 5
                lines.append(f"  - {os.path.basename(f)}")
            if len(output_files) > 5:
                lines.append(f"  ... and {len(output_files) - 5} more")

        return "\n".join(lines)

    def format_all_cycles(self):
        """
        Format all cycles for display.

        Returns:
            str: Formatted summary of all cycles
        """
        if not self.data["cycles"]:
            return "No cycles recorded."

        lines = [
            "",
            f"{'#'*60}",
            f"AGENT SESSION SUMMARY",
            f"{'#'*60}",
            f"Session ID: {self.data.get('session_id', 'N/A')}",
            f"Project Advice: {self.data.get('project_advice', 'None')}",
            f"Total Cycles: {len(self.data['cycles'])}",
        ]

        for i in range(1, len(self.data["cycles"]) + 1):
            lines.append(self.format_cycle_summary(i))

        if self.data.get("summary"):
            lines.extend([
                "",
                f"{'='*60}",
                "AI SUMMARY OF SESSION",
                f"{'='*60}",
                self.data["summary"]
            ])

        lines.append(f"{'#'*60}")
        return "\n".join(lines)

    def generate_log_for_summary(self):
        """
        Generate a log-like text suitable for LLM summarization.
        Consolidates metrics to show only final/best values.

        Returns:
            str: Log text for summarization
        """
        # Determine experiment type from files and history
        original_files = self.data.get('original_files', [])
        experiment_type = self._detect_experiment_type_for_summary(original_files)

        lines = [
            f"WORKING DIRECTORY:{os.getcwd()}",
            "COMMAND THAT WAS RUN: phenix.run_agent",
            "PHENIX AGENT SESSION LOG",
            f"Experiment Type: {experiment_type}",
            f"Project Advice: {self.data.get('project_advice', 'None')}",
            f"Original Files: {', '.join(original_files)}",
            f"Total Cycles: {len(self.data['cycles'])}",
        ]

        # Add resolution if set
        resolution = self.data.get("resolution")
        if resolution:
            lines.append(f"Resolution: {resolution:.2f} Å")

        lines.append("")

        # Collect key metrics by program (keep only final values)
        metrics_by_program = {}
        final_output_files = []
        programs_run = []

        for cycle in self.data["cycles"]:
            program = cycle.get('program', 'unknown')
            if program and program not in ['STOP', 'unknown', None]:
                if program not in programs_run:
                    programs_run.append(program)

            # Extract metrics from result text
            result = cycle.get('result', '')
            if isinstance(result, str):
                # Parse key metrics from result
                metrics = self._extract_metrics_from_result(result, program)
                if metrics:
                    metrics_by_program[program] = metrics  # Overwrite with latest

            # Collect output files (final cycle's files are most relevant)
            output_files = cycle.get('output_files', [])
            if output_files:
                final_output_files = output_files  # Keep updating to get latest

        # Programs run
        lines.extend([
            "="*40,
            "PROGRAMS RUN",
            "="*40,
            ", ".join(programs_run),
            "",
        ])

        # Consolidated metrics table (final values only)
        lines.extend([
            "="*40,
            "FINAL METRICS (most recent values)",
            "="*40,
        ])

        for program, metrics in metrics_by_program.items():
            if metrics:
                lines.append(f"\n{program}:")
                for key, value in metrics.items():
                    lines.append(f"  {key}: {value}")

        # Brief cycle summary (just program + outcome)
        lines.extend([
            "",
            "="*40,
            "CYCLE HISTORY (brief)",
            "="*40,
        ])

        for cycle in self.data["cycles"]:
            cycle_num = cycle.get('cycle_number', '?')
            program = cycle.get('program', 'N/A')
            result = cycle.get('result', 'N/A')
            # Truncate result to first line / 80 chars
            if isinstance(result, str):
                result_brief = result.split('\n')[0][:80]
            else:
                result_brief = str(result)[:80]
            lines.append(f"Cycle {cycle_num}: {program} - {result_brief}")

        # Final output files
        if final_output_files:
            lines.extend([
                "",
                "="*40,
                "KEY OUTPUT FILES",
                "="*40,
            ])
            for f in final_output_files[:10]:
                lines.append(f"  {os.path.basename(f)}")

        return "\n".join(lines)

    def _detect_experiment_type_for_summary(self, files):
        """
        Detect experiment type from files and history for summary generation.

        Args:
            files: List of file paths

        Returns:
            str: "Cryo-EM" or "X-ray Crystallography"
        """
        # Check if mtriage or real_space_refine was run (definitive cryo-EM indicators)
        for cycle in self.data.get("cycles", []):
            program = cycle.get("program", "").lower()
            if "mtriage" in program or "real_space_refine" in program:
                return "Cryo-EM"

        # Check file extensions
        has_mtz = False
        has_map = False

        for f in files:
            f_lower = f.lower()
            basename = os.path.basename(f_lower)

            # X-ray data files
            if f_lower.endswith(('.mtz', '.sca', '.hkl')):
                has_mtz = True

            # Map files (cryo-EM)
            if f_lower.endswith(('.ccp4', '.mrc', '.map')):
                has_map = True

        if has_map and not has_mtz:
            return "Cryo-EM"
        elif has_mtz:
            return "X-ray Crystallography"
        else:
            return "Unknown"

    def _extract_metrics_from_result(self, result_text, program):
        """Extract key metrics from result text."""
        import re
        metrics = {}

        if not result_text:
            return metrics

        # Common metric patterns
        # IMPORTANT: Use lowercase keys to match workflow_state expectations
        patterns = {
            'r_free': r'R.?free[:\s=]+([0-9.]+)',
            'r_work': r'R.?work[:\s=]+([0-9.]+)',
            'resolution': r'[Rr]esolution[:\s=]+([0-9.]+)',
            'clashscore': r'[Cc]lashscore[:\s=]+([0-9.]+)',
            'bonds_rmsd': r'bonds.?rmsd[:\s=]+([0-9.]+)',
            'angles_rmsd': r'angles.?rmsd[:\s=]+([0-9.]+)',
            'tfz': r'TFZ[:\s=]+([0-9.]+)',
            'llg': r'LLG[:\s=]+([0-9.]+)',
            'map_cc': r'(?:map.?model.?)?CC[:\s=]+([0-9.]+)',
            'completeness': r'[Cc]ompleteness[:\s=]+([0-9.]+)',
        }

        for name, pattern in patterns.items():
            match = re.search(pattern, result_text, re.IGNORECASE)
            if match:
                try:
                    value = float(match.group(1))
                    # Skip obviously wrong values
                    if name == 'resolution' and value > 10:
                        continue  # Likely wrong (e.g., 47.31 is not resolution)
                    metrics[name] = value
                except ValueError:
                    pass

        return metrics

    def generate_summary(self, llm):
        """
        Use LLM to generate a summary of the session (synchronous version).

        Args:
            llm: Language model to use

        Returns:
            str: Generated summary
        """
        log_text = self.generate_log_for_summary()


        # Define the template using a standard string
        template = """You are an expert crystallographer reviewing an automated Phenix structure determination session.

Create a concise final report from the following session log.

IMPORTANT FORMATTING RULES:
- Do NOT repeat metrics multiple times - show only the FINAL values
- Keep the report brief and focused on outcomes
- Use a clean format without excessive duplication

Include these sections:
1. **Summary**: What was accomplished in 2-3 sentences
2. **Final Model Quality**: Key metrics (R-free, R-work, resolution, clashscore) - ONE set of values only
3. **Key Output Files**: The most important final files (refined model, MTZ)
4. **Status**: Whether the workflow completed successfully

SESSION LOG:
<session_log>
{log_text}
</session_log>

FINAL REPORT:"""

        # Create the safe prompt object
        prompt_template = PromptTemplate.from_template(template)

        # Generate the string safely
        # .format handles escaping if log_text contains { or }
        prompt = prompt_template.format(log_text=log_text)

        try:
            response = llm.invoke(prompt)
            summary = response.content.strip()
            self.data["summary"] = summary
            self.save()
            return summary
        except Exception as e:
            error_msg = f"Could not generate summary: {e}"
            print(error_msg)
            return error_msg

    def to_dict(self):
        """Return the session data as a dictionary."""
        return self.data.copy()

    def get_available_files(self):
        """
        Compute the list of available files from session data.

        This is the SINGLE SOURCE OF TRUTH for what files the agent can use.
        It combines:
        - original_files: The user's initial input files
        - output_files from each cycle: Files produced by completed cycles

        When cycles are removed (e.g., via session_utils --remove-last),
        this automatically returns the correct reduced file list.

        Returns:
            list: Absolute paths to all available files
        """
        import os
        files = []
        seen = set()  # Track by basename to avoid duplicates

        # 1. Start with original input files
        for f in self.data.get("original_files", []):
            if f:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                basename = os.path.basename(abs_path)
                if basename not in seen:
                    files.append(abs_path)
                    seen.add(basename)

        # 2. Add output files from each cycle (in order)
        for cycle in self.data.get("cycles", []):
            for f in cycle.get("output_files", []):
                if f:
                    abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                    basename = os.path.basename(abs_path)
                    # Only add if file exists and not already tracked
                    if basename not in seen and os.path.exists(abs_path):
                        files.append(abs_path)
                        seen.add(basename)

        return files

    def reset(self):
        """Reset the session (start fresh)."""
        self._init_new_session()
        self.save()

    # =========================================================================
    # AGENT SESSION SUMMARY
    # =========================================================================

    def generate_agent_session_summary(self, include_llm_assessment=True):
        """
        Generate a structured summary of the agent session.

        This produces a clean Markdown summary with:
        - Input files and user advice
        - Workflow path taken
        - Steps performed with key metrics
        - Final quality assessment
        - Output files

        Args:
            include_llm_assessment: If True, include placeholder for LLM assessment

        Returns:
            dict with:
                - markdown: The formatted Markdown summary
                - data: Structured data for LLM assessment
        """
        # Extract all structured data
        summary_data = self._extract_summary_data()

        # Generate Markdown
        markdown = self._format_summary_markdown(summary_data, include_llm_assessment)

        return {
            "markdown": markdown,
            "data": summary_data,
        }

    def _extract_summary_data(self):
        """
        Extract structured summary data from the session.

        Returns:
            dict with all summary components
        """
        cycles = self.data.get("cycles", [])
        original_files = self.data.get("original_files", [])

        # Detect experiment type
        experiment_type = self.data.get("experiment_type")
        if not experiment_type:
            experiment_type = self._detect_experiment_type_for_summary(original_files)

        # Determine workflow path
        workflow_path = self._determine_workflow_path(cycles, experiment_type)

        # Extract steps with metrics
        steps = self._extract_steps(cycles)

        # Get final metrics (from last successful refinement or analysis)
        final_metrics = self._extract_final_metrics(cycles, experiment_type)

        # Get final output files
        final_files = self._get_final_output_files(cycles)

        # Get input data quality metrics (from xtriage/mtriage)
        input_quality = self._extract_input_quality(cycles, experiment_type)

        return {
            "session_id": self.data.get("session_id", "unknown"),
            "experiment_type": experiment_type,
            "original_files": [os.path.basename(f) for f in original_files],
            "user_advice": self.data.get("project_advice", "None provided"),
            "resolution": self.data.get("resolution"),
            "workflow_path": workflow_path,
            "steps": steps,
            "final_metrics": final_metrics,
            "final_files": final_files,
            "input_quality": input_quality,
            "total_cycles": len(cycles),
            "successful_cycles": sum(1 for c in cycles if "SUCCESS" in str(c.get("result", ""))),
        }

    def _determine_workflow_path(self, cycles, experiment_type):
        """Determine the workflow path taken."""
        programs = [c.get("program", "") for c in cycles if c.get("program")]

        if experiment_type == "cryoem" or "Cryo-EM" in str(experiment_type):
            if "phenix.predict_and_build" in programs:
                return "Cryo-EM model building with AlphaFold prediction"
            elif "phenix.dock_in_map" in programs:
                return "Cryo-EM docking of existing model"
            else:
                return "Cryo-EM refinement"
        else:
            # X-ray
            if "phenix.autosol" in programs:
                return "X-ray experimental phasing (SAD/MAD)"
            elif "phenix.predict_and_build" in programs and "phenix.phaser" in programs:
                return "X-ray molecular replacement with AlphaFold prediction"
            elif "phenix.phaser" in programs:
                return "X-ray molecular replacement"
            elif "phenix.autobuild" in programs:
                return "X-ray automated model building"
            else:
                return "X-ray structure determination"

    def _extract_steps(self, cycles):
        """Extract step-by-step summary."""
        steps = []
        for cycle in cycles:
            program = cycle.get("program", "unknown")
            if program in ["STOP", None, "unknown"]:
                continue

            result = cycle.get("result", "")
            success = "SUCCESS" in str(result)

            # Extract key metric for this step
            key_metric = self._get_key_metric_for_step(cycle, program)

            # Get brief decision (first sentence)
            decision = cycle.get("decision", "")
            if decision:
                decision_brief = decision.split(".")[0] + "."
                if len(decision_brief) > 100:
                    decision_brief = decision_brief[:97] + "..."
            else:
                decision_brief = ""

            steps.append({
                "cycle": cycle.get("cycle_number", "?"),
                "program": program,
                "success": success,
                "key_metric": key_metric,
                "decision_brief": decision_brief,
            })

        return steps

    def _get_key_metric_for_step(self, cycle, program):
        """Get the most important metric for a given step."""
        result = str(cycle.get("result", ""))

        # Program-specific key metrics
        if "xtriage" in program.lower():
            # Look for high-resolution limit (d_min) with specific patterns
            # Try "High resolution limit" first (most specific)
            match = re.search(r'High.?resolution.?limit[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                return f"Resolution: {match.group(1)} Å"
            # Try generic "Resolution" but filter out d_max values (> 15 Å)
            match = re.search(r'Resolution[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                res = float(match.group(1))
                if res < 15:  # d_min is typically < 10 Å, d_max is typically > 20 Å
                    return f"Resolution: {match.group(1)} Å"
            match = re.search(r'Completeness[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                return f"Completeness: {match.group(1)}%"

        elif "mtriage" in program.lower():
            match = re.search(r'd_fsc[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                return f"Resolution (FSC): {match.group(1)} Å"

        elif "predict_and_build" in program.lower():
            match = re.search(r'[Pp]lddt[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                return f"pLDDT: {match.group(1)}"

        elif "phaser" in program.lower():
            match = re.search(r'TFZ[:\s=]+([0-9.]+)', result, re.IGNORECASE)
            if match:
                return f"TFZ: {match.group(1)}"
            # Fallback to space group if TFZ not found
            # Match space group like "P 32 2 1" - letters, numbers, and spaces but stop at newline
            match = re.search(r'Space.?Group[:\s=]+([A-Z][A-Z0-9 ]*)', result, re.IGNORECASE)
            if match:
                sg = match.group(1).strip()
                # Clean up - only keep the space group designation
                sg = sg.split('\n')[0].strip()
                return f"SG: {sg}"

        elif "refine" in program.lower():
            match = re.search(r'R.?[Ff]ree[:\s=]+([0-9.]+)', result)
            if match:
                return f"R-free: {match.group(1)}"

        elif "real_space_refine" in program.lower():
            match = re.search(r'[Mm]ap.?[Cc][Cc][:\s=]+([0-9.]+)', result)
            if match:
                return f"Map CC: {match.group(1)}"

        elif "autobuild" in program.lower():
            match = re.search(r'R.?[Ff]ree[:\s=]+([0-9.]+)', result)
            if match:
                return f"R-free: {match.group(1)}"

        elif "molprobity" in program.lower():
            # Try clashscore first
            match = re.search(r'[Cc]lashscore[:\s=]+([0-9.]+)', result)
            if match:
                clashscore = match.group(1)
                # Also try to get Ramachandran outliers
                rama_match = re.search(r'[Rr]amachandran.?[Oo]utliers?[:\s=]+([0-9.]+)', result)
                if rama_match:
                    return f"Clash: {clashscore}, Rama: {rama_match.group(1)}%"
                return f"Clashscore: {clashscore}"
            # Fallback to Ramachandran
            match = re.search(r'[Rr]amachandran.?[Oo]utliers?[:\s=]+([0-9.]+)', result)
            if match:
                return f"Rama outliers: {match.group(1)}%"

        return ""

    def _extract_final_metrics(self, cycles, experiment_type):
        """Extract final quality metrics from the last relevant cycles."""
        metrics = {}

        # Work backwards through cycles to get most recent values
        for cycle in reversed(cycles):
            result = str(cycle.get("result", ""))
            program = cycle.get("program", "").lower()

            # R-factors (X-ray)
            if "r_free" not in metrics:
                match = re.search(r'R.?[Ff]ree[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["r_free"] = float(match.group(1))

            if "r_work" not in metrics:
                match = re.search(r'R.?[Ww]ork[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["r_work"] = float(match.group(1))

            # Map CC (Cryo-EM)
            if "map_cc" not in metrics:
                match = re.search(r'[Mm]ap.?[Cc][Cc][:\s=]+([0-9.]+)', result)
                if match:
                    metrics["map_cc"] = float(match.group(1))

            # Geometry
            if "clashscore" not in metrics:
                match = re.search(r'[Cc]lashscore[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["clashscore"] = float(match.group(1))

            if "bonds_rmsd" not in metrics:
                match = re.search(r'[Bb]onds.?[Rr]msd[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["bonds_rmsd"] = float(match.group(1))

            if "angles_rmsd" not in metrics:
                match = re.search(r'[Aa]ngles.?[Rr]msd[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["angles_rmsd"] = float(match.group(1))

        # Add quality assessments
        if "r_free" in metrics:
            metrics["r_free_assessment"] = self._assess_r_free(
                metrics["r_free"], self.data.get("resolution"))

        if "map_cc" in metrics:
            metrics["map_cc_assessment"] = self._assess_map_cc(metrics["map_cc"])

        if "clashscore" in metrics:
            metrics["clashscore_assessment"] = self._assess_clashscore(metrics["clashscore"])

        return metrics

    def _assess_r_free(self, r_free, resolution):
        """Assess R-free quality based on resolution."""
        if resolution is None:
            resolution = 2.5  # Default assumption

        # Resolution-dependent thresholds
        if resolution < 1.5:
            good, acceptable = 0.20, 0.25
        elif resolution < 2.5:
            good, acceptable = 0.25, 0.30
        elif resolution < 3.5:
            good, acceptable = 0.30, 0.35
        else:
            good, acceptable = 0.35, 0.40

        if r_free <= good:
            return "Good"
        elif r_free <= acceptable:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _assess_map_cc(self, map_cc):
        """Assess map-model correlation."""
        if map_cc >= 0.80:
            return "Good"
        elif map_cc >= 0.70:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _assess_clashscore(self, clashscore):
        """Assess geometry clashscore."""
        if clashscore <= 5:
            return "Excellent"
        elif clashscore <= 10:
            return "Good"
        elif clashscore <= 20:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _get_final_output_files(self, cycles):
        """Get the most relevant final output files."""
        final_files = []

        # Work backwards to find the last cycle with output files
        for cycle in reversed(cycles):
            output_files = cycle.get("output_files", [])
            if output_files:
                # Prioritize certain file types
                for f in output_files:
                    basename = os.path.basename(f)
                    # Prioritize refined models and data
                    if any(x in basename.lower() for x in
                           ["refine", "overall_best", "real_space_refined"]):
                        if basename.endswith(".pdb"):
                            final_files.insert(0, {"name": basename, "type": "model"})
                        elif basename.endswith(".mtz"):
                            final_files.append({"name": basename, "type": "data"})
                    elif basename.endswith(".pdb") and len(final_files) < 5:
                        final_files.append({"name": basename, "type": "model"})
                    elif basename.endswith(".mtz") and len(final_files) < 5:
                        final_files.append({"name": basename, "type": "data"})

                if final_files:
                    break

        return final_files[:5]  # Limit to 5 files

    def _extract_input_quality(self, cycles, experiment_type):
        """Extract input data quality metrics from xtriage/mtriage."""
        quality = {}

        for cycle in cycles:
            program = cycle.get("program", "").lower()
            result = str(cycle.get("result", ""))

            if "xtriage" in program:
                # Resolution - look for d_min (high-resolution limit)
                # Try "High resolution limit" first (most specific)
                match = re.search(r'High.?resolution.?limit[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                if match:
                    quality["resolution"] = float(match.group(1))
                else:
                    # Fallback to generic "Resolution" but filter out d_max values
                    match = re.search(r'Resolution[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match:
                        res = float(match.group(1))
                        if res < 15:  # d_min is typically < 10 Å
                            quality["resolution"] = res

                # Completeness
                match = re.search(r'Completeness[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                if match:
                    quality["completeness"] = float(match.group(1))

                # Twinning
                if "No Twinning" in result or "twin" not in result.lower():
                    quality["twinning"] = "None detected"
                else:
                    match = re.search(r'Twin.?Fraction[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match and float(match.group(1)) > 0.1:
                        quality["twinning"] = f"Possible (fraction: {match.group(1)})"
                    else:
                        quality["twinning"] = "None detected"

                break  # Only need first xtriage

            elif "mtriage" in program:
                # FSC resolution
                match = re.search(r'd_fsc[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                if match:
                    quality["resolution"] = float(match.group(1))

                # d99
                match = re.search(r'd99[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                if match:
                    quality["d99"] = float(match.group(1))

                break  # Only need first mtriage

        return quality

    def _format_summary_markdown(self, data, include_llm_assessment):
        """Format the summary data as Markdown."""
        lines = []

        # Header - extract date from session_id (format: YYYY-MM-DD_HH-MM-SS)
        session_id = data['session_id']
        # Try to extract just the date part
        if '_' in session_id:
            date_part = session_id.split('_')[0]  # e.g., "2026-01-19"
        else:
            date_part = session_id

        lines.append(f"# Phenix AI Agent Run {date_part}")
        lines.append("")
        lines.append(f"**Session ID:** {session_id}")
        lines.append(f"**Total Cycles:** {data['total_cycles']} ({data['successful_cycles']} successful)")
        if include_llm_assessment:
            lines.append("")
            lines.append("_[STATUS_PLACEHOLDER]_")
        lines.append("")

        # Input Section
        lines.append("## Input")
        lines.append("")
        lines.append(f"- **Files:** {', '.join(data['original_files'])}")
        lines.append(f"- **User Advice:** {data['user_advice']}")
        lines.append(f"- **Experiment Type:** {data['experiment_type']}")
        if data.get("resolution"):
            lines.append(f"- **Resolution:** {data['resolution']:.2f} Å")
        lines.append("")

        # Input Data Quality
        if data.get("input_quality"):
            lines.append("## Input Data Quality")
            lines.append("")
            iq = data["input_quality"]
            if "resolution" in iq:
                lines.append(f"- **Resolution:** {iq['resolution']:.2f} Å")
            if "completeness" in iq:
                lines.append(f"- **Completeness:** {iq['completeness']:.1f}%")
            if "twinning" in iq:
                lines.append(f"- **Twinning:** {iq['twinning']}")
            if "d99" in iq:
                lines.append(f"- **d99:** {iq['d99']:.2f} Å")
            lines.append("")

        # Workflow Path
        lines.append("## Workflow Path")
        lines.append("")
        lines.append(data["workflow_path"])
        lines.append("")

        # Steps Performed
        lines.append("## Steps Performed")
        lines.append("")
        lines.append("| Cycle | Program | Result | Key Metric |")
        lines.append("|-------|---------|--------|------------|")
        for step in data["steps"]:
            result_symbol = "✓" if step["success"] else "✗"
            program_short = step["program"].replace("phenix.", "")
            lines.append(f"| {step['cycle']} | {program_short} | {result_symbol} | {step['key_metric']} |")
        lines.append("")

        # Final Quality
        if data.get("final_metrics"):
            lines.append("## Final Quality")
            lines.append("")
            lines.append("| Metric | Value | Assessment |")
            lines.append("|--------|-------|------------|")
            fm = data["final_metrics"]
            if "r_free" in fm:
                lines.append(f"| R-free | {fm['r_free']:.4f} | {fm.get('r_free_assessment', '')} |")
            if "r_work" in fm:
                lines.append(f"| R-work | {fm['r_work']:.4f} | |")
            if "map_cc" in fm:
                lines.append(f"| Map CC | {fm['map_cc']:.3f} | {fm.get('map_cc_assessment', '')} |")
            if "clashscore" in fm:
                lines.append(f"| Clashscore | {fm['clashscore']:.2f} | {fm.get('clashscore_assessment', '')} |")
            if "bonds_rmsd" in fm:
                lines.append(f"| Bonds RMSD | {fm['bonds_rmsd']:.4f} | |")
            if "angles_rmsd" in fm:
                lines.append(f"| Angles RMSD | {fm['angles_rmsd']:.3f} | |")
            lines.append("")

        # Output Files
        if data.get("final_files"):
            lines.append("## Output Files")
            lines.append("")
            for f in data["final_files"]:
                lines.append(f"- {f['name']} ({f['type']})")
            lines.append("")

        # LLM Assessment placeholder
        if include_llm_assessment:
            lines.append("## Assessment")
            lines.append("")
            lines.append("_[LLM assessment will be inserted here]_")
            lines.append("")

        return "\n".join(lines)

    def get_summary_for_llm_assessment(self):
        """
        Get a concise summary suitable for LLM assessment.

        Returns a text block with key information for the LLM to evaluate:
        - Input data quality
        - Goal/strategy
        - Steps taken
        - Current status
        """
        data = self._extract_summary_data()

        lines = []
        lines.append("=== AGENT SESSION FOR ASSESSMENT ===")
        lines.append("")

        # Input
        lines.append(f"EXPERIMENT TYPE: {data['experiment_type']}")
        lines.append(f"INPUT FILES: {', '.join(data['original_files'])}")
        lines.append(f"USER GOAL: {data['user_advice']}")
        lines.append("")

        # Input quality
        if data.get("input_quality"):
            lines.append("INPUT DATA QUALITY:")
            for k, v in data["input_quality"].items():
                lines.append(f"  {k}: {v}")
            lines.append("")

        # Workflow
        lines.append(f"WORKFLOW PATH: {data['workflow_path']}")
        lines.append("")

        # Steps (brief)
        lines.append("STEPS TAKEN:")
        for step in data["steps"]:
            status = "OK" if step["success"] else "FAILED"
            lines.append(f"  {step['cycle']}. {step['program']} - {status} {step['key_metric']}")
        lines.append("")

        # Final metrics
        if data.get("final_metrics"):
            lines.append("FINAL METRICS:")
            fm = data["final_metrics"]
            if "r_free" in fm:
                lines.append(f"  R-free: {fm['r_free']:.4f} ({fm.get('r_free_assessment', '')})")
            if "r_work" in fm:
                lines.append(f"  R-work: {fm['r_work']:.4f}")
            if "map_cc" in fm:
                lines.append(f"  Map CC: {fm['map_cc']:.3f} ({fm.get('map_cc_assessment', '')})")
            if "clashscore" in fm:
                lines.append(f"  Clashscore: {fm['clashscore']:.2f} ({fm.get('clashscore_assessment', '')})")
            lines.append("")

        lines.append(f"TOTAL CYCLES: {data['total_cycles']} ({data['successful_cycles']} successful)")

        return "\n".join(lines)
