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
    session.record_result(result='Success')
    session.save()
"""
from __future__ import absolute_import, division, print_function

import os
import json
from datetime import datetime


class AgentSession:
    """
    Manages persistent state for an agent session across multiple cycles.

    The session is stored as a JSON file and tracks:
    - Project info (advice, original files)
    - Each cycle's decision, command, and result
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
            "summary": ""
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

    def record_result(self, cycle_number, result):
        """
        Record the result of running a cycle's command.

        Args:
            cycle_number: Which cycle to update
            result: Result text (success/failure description)
        """
        self.start_cycle(cycle_number)
        self.data["cycles"][cycle_number - 1]["result"] = result
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

            # Only include cycles that ran successfully
            result = cycle.get("result", "")
            if not result.startswith("SUCCESS"):
                continue

            cmd = cycle.get("command", "").strip()
            if cmd and cmd != "No command generated.":
                # Normalize: collapse whitespace
                norm_cmd = " ".join(cmd.split())
                commands.append((cycle.get("cycle_number", 0), norm_cmd))
        return commands

    def cleanup_crash_cycles(self):
        """
        Remove CRASH cycles and cycles that never ran successfully from the session.
        Call this at the start of a new run to clean up from previous failures.
        """
        original_count = len(self.data["cycles"])
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
        Get cycle history in format expected by generate_next_move.

        Returns:
            list: List of dicts compatible with run_history format
        """
        history = []
        for cycle in self.data["cycles"]:
            if cycle.get("command") and cycle.get("result"):
                history.append({
                    "job_id": str(cycle.get("cycle_number", "")),
                    "program": cycle.get("program", ""),
                    "summary": f"Decision: {cycle.get('decision', '')}\n"
                               f"Command: {cycle.get('command', '')}\n"
                               f"Result: {cycle.get('result', '')}",
                    "analysis": cycle.get("reasoning", ""),
                    "next_move": {
                        "command": cycle.get("command", ""),
                        "program": cycle.get("program", "")
                    }
                })
        return history

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

        Returns:
            str: Log text for summarization
        """
        lines = [
            f"WORKING DIRECTORY:{os.getcwd()}",
            "COMMAND THAT WAS RUN: phenix.run_agent",
            "PHENIX AGENT SESSION LOG",
            f"Project Advice: {self.data.get('project_advice', 'None')}",
            f"Original Files: {', '.join(self.data.get('original_files', []))}",
            f"Total Cycles: {len(self.data['cycles'])}",
            "",
            "="*40,
            "CYCLE-BY-CYCLE LOG",
            "="*40,
        ]

        for cycle in self.data["cycles"]:
            lines.extend([
                "",
                f"--- Cycle {cycle.get('cycle_number', '?')} ---",
                f"Program: {cycle.get('program', 'N/A')}",
                f"Decision: {cycle.get('decision', 'N/A')}",
                f"Reasoning: {cycle.get('reasoning', 'N/A')}",
                f"Command: {cycle.get('command', 'N/A')}",
                f"Result: {cycle.get('result', 'N/A')}",
            ])

        return "\n".join(lines)

    def generate_summary(self, llm):
        """
        Use LLM to generate a summary of the session (synchronous version).

        Args:
            llm: Language model to use

        Returns:
            str: Generated summary
        """
        log_text = self.generate_log_for_summary()

        prompt = f"""You are an expert crystallographer reviewing an automated Phenix structure determination session.

Summarize the following agent session log. Include:
1. Overall goal and what was accomplished
2. Key decisions and their outcomes (successes and failures)
3. Current state of the structure determination
4. Recommended next steps (if work remains)

Be concise but informative. Use bullet points for clarity.

SESSION LOG:
{log_text}

SUMMARY:"""

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

    def reset(self):
        """Reset the session (start fresh)."""
        self._init_new_session()
        self.save()
