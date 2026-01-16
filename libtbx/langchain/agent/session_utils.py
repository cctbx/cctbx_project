#!/usr/bin/env python
"""
Utility script to manage agent session history.

Usage:
    phenix.python session_utils.py --show                    # Show current session (brief)
    phenix.python session_utils.py --detailed                # Show detailed session info
    phenix.python session_utils.py --remove-last 2          # Remove last 2 cycles
    phenix.python session_utils.py --reset                   # Reset entire session
    phenix.python session_utils.py --dir ai_agent_directory  # Specify session directory

The --detailed (-v) flag shows for each cycle:
    - Full command
    - Reasoning/decision
    - Output files (with existence check)
    - Extracted metrics (R-free, R-work, resolution, etc.)
    - Cumulative file list
"""
from __future__ import absolute_import, division, print_function

import argparse
import os
import sys
import json


def load_session(session_dir):
    """Load session from directory."""
    session_file = os.path.join(session_dir, "agent_session.json")
    if not os.path.exists(session_file):
        print("No session file found at: %s" % session_file)
        return None, session_file

    with open(session_file, 'r') as f:
        data = json.load(f)
    return data, session_file


def save_session(data, session_file):
    """Save session to file."""
    with open(session_file, 'w') as f:
        json.dump(data, f, indent=2)
    print("Session saved to: %s" % session_file)


def show_session(data, detailed=False):
    """Display session summary."""
    print("\n" + "="*60)
    print("AGENT SESSION SUMMARY")
    print("="*60)
    print("Session ID: %s" % data.get("session_id", "N/A"))
    print("Project Advice: %s" % (data.get("project_advice", "None") or "None"))
    print("Total Cycles: %d" % len(data.get("cycles", [])))
    print("-"*60)

    # Show original files
    original_files = data.get("original_files", [])
    print("\nORIGINAL INPUT FILES (%d):" % len(original_files))
    for f in original_files:
        print("  - %s" % f)

    if not detailed:
        # Brief view
        print("\n" + "-"*60)
        for cycle in data.get("cycles", []):
            cycle_num = cycle.get("cycle_number", "?")
            program = cycle.get("program", "N/A")
            command = cycle.get("command", "N/A")
            result = cycle.get("result", "N/A")

            # Truncate long values
            if len(command) > 60:
                command = command[:60] + "..."
            if len(result) > 40:
                result = result[:40] + "..."

            status = "✓" if result.startswith("SUCCESS") else "✗"
            print("%s Cycle %s: %s" % (status, cycle_num, program))
            print("    Command: %s" % command)
            print("    Result: %s" % result)
    else:
        # Detailed view
        print("\n")
        for cycle in data.get("cycles", []):
            show_cycle_detailed(cycle)

    # Show cumulative file list
    print("\n" + "="*60)
    print("CUMULATIVE AVAILABLE FILES")
    print("="*60)
    all_files = get_available_files_from_data(data)
    for i, f in enumerate(all_files, 1):
        print("  %d. %s" % (i, f))
    print("Total: %d files" % len(all_files))

    # Show summary if present
    if data.get("summary"):
        print("\n" + "="*60)
        print("SESSION SUMMARY (AI-generated)")
        print("="*60)
        print(data["summary"])

    print("="*60 + "\n")


def show_cycle_detailed(cycle):
    """Display detailed information for a single cycle."""
    cycle_num = cycle.get("cycle_number", "?")
    program = cycle.get("program", "N/A")
    result = cycle.get("result", "N/A")

    status = "✓" if result.startswith("SUCCESS") else "✗"

    print("=" * 60)
    print("%s CYCLE %s: %s" % (status, cycle_num, program))
    print("=" * 60)

    # Timestamp
    if cycle.get("timestamp"):
        print("Timestamp: %s" % cycle["timestamp"])

    # Decision/Reasoning
    if cycle.get("decision"):
        print("\nDECISION: %s" % cycle["decision"])

    if cycle.get("reasoning"):
        print("\nREASONING:")
        # Wrap long reasoning text
        reasoning = cycle["reasoning"]
        if len(reasoning) > 500:
            reasoning = reasoning[:500] + "..."
        print("  %s" % reasoning)

    # Full command
    print("\nCOMMAND:")
    command = cycle.get("command", "N/A")
    print("  %s" % command)

    # Result
    print("\nRESULT:")
    result_text = cycle.get("result", "N/A")
    if len(result_text) > 300:
        result_text = result_text[:300] + "..."
    print("  %s" % result_text)

    # Output files from this cycle
    output_files = cycle.get("output_files", [])
    if output_files:
        print("\nOUTPUT FILES (%d):" % len(output_files))
        for f in output_files:
            # Check if file exists
            exists = os.path.exists(f) if f else False
            status_mark = "✓" if exists else "✗"
            print("  %s %s" % (status_mark, f))
    else:
        print("\nOUTPUT FILES: None recorded")

    # Extract and show any metrics from the result
    metrics = extract_metrics_from_result(cycle.get("result", ""))
    if metrics:
        print("\nMETRICS EXTRACTED:")
        for key, value in metrics.items():
            print("  %s: %s" % (key, value))

    print("")


def extract_metrics_from_result(result_text):
    """Extract key metrics from result text."""
    import re
    metrics = {}

    if not result_text:
        return metrics

    # Common metric patterns
    patterns = {
        'R-free': r'R.?free[:\s=]+([0-9.]+)',
        'R-work': r'R.?work[:\s=]+([0-9.]+)',
        'Resolution': r'resolution[:\s=]+([0-9.]+)\s*[AÅ]?',
        'Clashscore': r'clashscore[:\s=]+([0-9.]+)',
        'Ramachandran favored': r'ramachandran.{0,20}favored[:\s=]+([0-9.]+)',
        'Ramachandran outliers': r'ramachandran.{0,20}outliers?[:\s=]+([0-9.]+)',
        'Rotamer outliers': r'rotamer.{0,20}outliers?[:\s=]+([0-9.]+)',
        'Bonds RMSD': r'bonds?.{0,10}rmsd[:\s=]+([0-9.]+)',
        'Angles RMSD': r'angles?.{0,10}rmsd[:\s=]+([0-9.]+)',
        'TFZ': r'TFZ[:\s=]+([0-9.]+)',
        'LLG': r'LLG[:\s=]+([0-9.]+)',
        'Map CC': r'(?:map.?model.?)?CC[:\s=]+([0-9.]+)',
        'Completeness': r'completeness[:\s=]+([0-9.]+)',
        'Anomalous signal': r'anomalous.{0,20}signal[:\s=]+([0-9.]+)',
    }

    for name, pattern in patterns.items():
        match = re.search(pattern, result_text, re.IGNORECASE)
        if match:
            try:
                value = float(match.group(1))
                metrics[name] = value
            except ValueError:
                pass

    return metrics


def get_available_files_from_data(data):
    """
    Compute available files from session data.
    Mirrors AgentSession.get_available_files() logic.
    """
    files = []
    seen = set()

    # 1. Original input files
    for f in data.get("original_files", []):
        if f:
            abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
            basename = os.path.basename(abs_path)
            if basename not in seen:
                files.append(abs_path)
                seen.add(basename)

    # 2. Output files from each cycle
    for cycle in data.get("cycles", []):
        for f in cycle.get("output_files", []):
            if f:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                basename = os.path.basename(abs_path)
                if basename not in seen:
                    files.append(abs_path)
                    seen.add(basename)

    return files


def remove_last_cycles(data, n, session_dir=None):
    """Remove last N cycles from session and update active_files.json."""
    original_count = len(data.get("cycles", []))
    if n <= 0:
        print("Nothing to remove (n=%d)" % n)
        return data, 0

    if n >= original_count:
        print("Removing all %d cycles" % original_count)
        data["cycles"] = []
    else:
        data["cycles"] = data["cycles"][:-n]

    # Renumber remaining cycles
    for i, cycle in enumerate(data["cycles"]):
        cycle["cycle_number"] = i + 1

    removed = original_count - len(data["cycles"])
    print("Removed last %d cycle(s). %d cycles remaining." % (removed, len(data["cycles"])))

    # Clear the AI-generated summary since it may reference removed cycles
    if data.get("summary"):
        data["summary"] = ""
        print("Cleared stale session summary (referenced removed cycles)")

    # Rebuild active_files.json to match the new session state
    # This ensures any code reading active_files.json stays in sync
    if session_dir:
        rebuild_active_files(data, session_dir)

    return data, removed


def rebuild_active_files(data, session_dir):
    """
    Rebuild active_files.json from session data.

    This ensures the active_files.json stays in sync with the session
    after cycles are removed or the session is modified.
    """
    import os

    files = []
    seen = set()

    # 1. Original input files
    for f in data.get("original_files", []):
        if f:
            abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
            basename = os.path.basename(abs_path)
            if basename not in seen:
                files.append(abs_path)
                seen.add(basename)

    # 2. Output files from remaining cycles
    for cycle in data.get("cycles", []):
        for f in cycle.get("output_files", []):
            if f:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                basename = os.path.basename(abs_path)
                if basename not in seen and os.path.exists(abs_path):
                    files.append(abs_path)
                    seen.add(basename)

    # Write to active_files.json
    active_files_path = os.path.join(session_dir, "active_files.json")
    try:
        with open(active_files_path, 'w') as f:
            json.dump(files, f, indent=2)
        print("Rebuilt active_files.json with %d files" % len(files))
    except Exception as e:
        print("Warning: Could not write active_files.json: %s" % e)


def reset_session(data, session_dir=None):
    """Reset session to empty state."""
    num_cycles = len(data.get("cycles", []))
    data["cycles"] = []
    data["summary"] = ""

    # Clear resolution (single source of truth from xtriage/mtriage)
    if "resolution" in data:
        del data["resolution"]
    if "resolution_source" in data:
        del data["resolution_source"]

    print("Session reset. Removed %d cycles." % num_cycles)

    # Rebuild active_files.json (will only contain original_files now)
    if session_dir:
        rebuild_active_files(data, session_dir)

    return data


def main():
    parser = argparse.ArgumentParser(description="Manage agent session history")
    parser.add_argument("--dir", "-d", default="ai_agent_directory",
                        help="Session directory (default: ai_agent_directory)")
    parser.add_argument("--show", "-s", action="store_true",
                        help="Show current session summary")
    parser.add_argument("--detailed", "-v", action="store_true",
                        help="Show detailed session info (files, metrics, reasoning for each cycle)")
    parser.add_argument("--remove-last", "-r", type=int, metavar="N",
                        help="Remove last N cycles from session")
    parser.add_argument("--reset", action="store_true",
                        help="Reset entire session (remove all cycles)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without saving")

    args = parser.parse_args()

    # Load session
    data, session_file = load_session(args.dir)
    if data is None:
        return 1

    modified = False

    # Show session (detailed or summary)
    if args.show or args.detailed or (not args.remove_last and not args.reset):
        show_session(data, detailed=args.detailed)

    # Remove last N cycles
    if args.remove_last:
        data, removed = remove_last_cycles(data, args.remove_last, session_dir=args.dir)
        if removed > 0:
            modified = True
            show_session(data, detailed=args.detailed)

    # Reset session
    if args.reset:
        data = reset_session(data, session_dir=args.dir)
        modified = True

    # Save if modified
    if modified and not args.dry_run:
        save_session(data, session_file)
    elif modified and args.dry_run:
        print("(Dry run - changes not saved)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

