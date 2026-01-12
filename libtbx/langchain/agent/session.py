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
import json
from datetime import datetime
from langchain_core.prompts import PromptTemplate


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
        }

    def load(self):
        """Load session from file."""
        try:
            with open(self.session_file, 'r') as f:
                self.data = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load session file: {e}")
            self._init_new_session()

    def save(self):
        """Save session to file."""
        try:
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

    def write_history_file(self, output_path=None):
        """
        Write session history as a job_*.json file for the agent to read.

        Args:
            output_path: Path to write to. If None, creates temp file.

        Returns:
            str: Path to the written file
        """
        import tempfile

        if output_path is None:
            fd, output_path = tempfile.mkstemp(suffix='.json', prefix='session_history_')
            os.close(fd)

        # Build a history record that looks like job_*.json
        # Include all successful cycles
        successful_cycles = [
            c for c in self.data["cycles"]
            if c.get("result", "").startswith("SUCCESS")
        ]

        if not successful_cycles:
            # No successful cycles, write minimal record
            history = {
                "job_id": "session_0",
                "timestamp": self.data.get("session_id", ""),
                "program": "session",
                "summary": f"Session started. Project advice: {self.data.get('project_advice', 'None')}",
                "analysis": "",
                "error": None,
                "next_move": None
            }
        else:
            # Build summary from all successful cycles
            cycle_summaries = []
            for c in successful_cycles:
                cycle_summaries.append(
                    f"Cycle {c['cycle_number']}: {c['program']} - {c['result']}\n"
                    f"Command: {c['command']}"
                )

            last_cycle = successful_cycles[-1]

            history = {
                "job_id": f"session_{len(successful_cycles)}",
                "timestamp": last_cycle.get("timestamp", ""),
                "program": last_cycle.get("program", ""),
                "summary": f"Session history ({len(successful_cycles)} successful cycles):\n" +
                          "\n\n".join(cycle_summaries),
                "analysis": self.data.get("summary", ""),
                "error": None,
                "next_move": {
                    "command": last_cycle.get("command", ""),
                    "program": last_cycle.get("program", "")
                }
            }

        with open(output_path, 'w') as f:
            json.dump(history, f, indent=2)

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
        
        self.save()

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
