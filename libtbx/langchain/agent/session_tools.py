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
    # Exclude STOP cycles from count - STOP is a decision, not a real program run
    real_cycles = [c for c in data.get("cycles", []) if c.get("program") not in ["STOP", None, "unknown"]]
    print("Total Cycles: %d" % len(real_cycles))
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

    # Show best files tracking
    best_files_data = data.get("best_files", {})
    if best_files_data:
        print("\n" + "="*60)
        print("BEST FILES TRACKING")
        print("="*60)
        best_entries = best_files_data.get("best", {})
        if best_entries:
            for category, entry in best_entries.items():
                if entry:
                    path = entry.get("path", "N/A")
                    stage = entry.get("stage", "")
                    cycle = entry.get("cycle", "?")
                    reason = entry.get("reason", "")
                    print("  %s: %s" % (category, os.path.basename(path)))
                    if detailed:
                        print("      Path: %s" % path)
                        print("      Stage: %s, Cycle: %s" % (stage, cycle))
                        if reason:
                            print("      Reason: %s" % reason)
        else:
            print("  (no best files tracked)")

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
        'Resolution': r'(?<!nomalous )(?<!nomalous  )resolution[:\s=]+([0-9.]+)\s*[AÅ]?',
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
    """Remove last N cycles from session and update active_files.json and best_files."""
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

    # Rebuild best_files tracking from remaining cycles
    rebuild_best_files(data)

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


def rebuild_best_files(data):
    """
    Rebuild best_files tracking from remaining cycles.

    After removing cycles, the best_files may reference files from
    removed cycles. This function rebuilds the tracking by replaying
    the remaining cycles through the BestFilesTracker.

    NOTE: This is a simplified version for session_tools. The full
    implementation is in Session._rebuild_best_files_from_cycles().
    """
    try:
        from libtbx.langchain.agent.best_files_tracker import BestFilesTracker
    except ImportError:
        try:
            from agent.best_files_tracker import BestFilesTracker
        except ImportError:
            print("Warning: Could not import BestFilesTracker, skipping best_files rebuild")
            return

    # Create a fresh tracker
    tracker = BestFilesTracker()

    # Helper to infer stage from program
    def infer_stage(program):
        program_lower = program.lower() if program else ""
        if "refine" in program_lower and "real_space" not in program_lower:
            return "refined"
        elif "real_space_refine" in program_lower:
            return "rsr_output"
        elif "phaser" in program_lower:
            return "phaser_output"
        elif "process_predicted_model" in program_lower:
            return "processed_predicted"
        elif "predict_and_build" in program_lower:
            return "predicted"
        elif "autobuild" in program_lower:
            return "autobuild_output"
        elif "dock_in_map" in program_lower:
            return "docked"
        elif "ligandfit" in program_lower:
            return "ligand_fit_output"
        elif "pdbtools" in program_lower:
            return "with_ligand"
        else:
            return None

    # First, evaluate original input files (cycle 0)
    print(f"  Evaluating {len(data.get('original_files', []))} original files (cycle 0):")
    for f in data.get("original_files", []):
        if f and os.path.exists(f):
            updated = tracker.evaluate_file(
                path=f,
                cycle=0,
                metrics=None,
                stage=None
            )
            if updated:
                basename = os.path.basename(f)
                category = tracker._classify_category(f)
                print(f"    -> {basename}: became best {category}")

    # Replay each remaining cycle's output files through the tracker
    for cycle in data.get("cycles", []):
        cycle_num = cycle.get("cycle_number", 1)
        program = cycle.get("program", "")
        result = cycle.get("result", "")
        metrics = cycle.get("metrics", {})
        output_files = cycle.get("output_files", [])

        # Skip failed cycles
        if "FAILED" in result.upper():
            print(f"  Cycle {cycle_num} ({program}): SKIPPED (failed)")
            continue

        # Determine stage from program name
        stage = infer_stage(program)

        print(f"  Cycle {cycle_num} ({program}): stage={stage}, files={len(output_files)}, metrics={list(metrics.keys()) if metrics else []}")

        # Debug: show all output files
        for f in output_files:
            exists = os.path.exists(f) if f else False
            print(f"    File: {os.path.basename(f) if f else 'None'} (exists={exists})")

        for f in output_files:
            if f and os.path.exists(f):
                # Build file-specific metrics
                file_metrics = dict(metrics) if metrics else {}

                # Determine file-specific stage based on file AND program
                # Only apply program-specific stage to files that match expected output patterns
                file_stage = None  # Default: let tracker infer from filename
                basename = os.path.basename(f).lower()

                if f.lower().endswith('.mtz'):
                    if 'refine' in basename or 'refinement_data' in basename:
                        file_stage = "refined_mtz"
                        file_metrics["has_rfree_flags"] = True
                elif f.lower().endswith('.pdb') or f.lower().endswith('.cif'):
                    # Only apply program stage if file basename matches expected output
                    if stage == "refined" and 'refine' in basename:
                        file_stage = "refined"
                    elif stage == "phaser_output" and ('phaser' in basename or basename.startswith('mr_')):
                        file_stage = "phaser_output"
                    elif stage == "predicted" and ('predict' in basename or 'alphafold' in basename):
                        file_stage = "predicted"
                    elif stage == "processed_predicted" and ('processed' in basename or 'trimmed' in basename):
                        file_stage = "processed_predicted"
                    elif stage == "docked" and ('dock' in basename or 'placed' in basename):
                        file_stage = "docked"
                    elif stage == "autobuild_output" and ('autobuild' in basename or 'overall_best' in basename):
                        file_stage = "autobuild_output"
                    elif stage == "ligand_fit_output" and ('ligand_fit' in basename or 'lig_fit' in basename):
                        file_stage = "ligand_fit_output"
                    elif stage == "with_ligand" and 'with_ligand' in basename:
                        file_stage = "with_ligand"
                    elif stage == "rsr_output" and ('real_space' in basename or 'rsr_' in basename):
                        file_stage = "rsr_output"
                    # else: file_stage stays None, tracker will infer from filename

                updated = tracker.evaluate_file(f, cycle=cycle_num, metrics=file_metrics, stage=file_stage)
                if updated:
                    print(f"    -> Updated best: {os.path.basename(f)} (stage={file_stage})")

    # Update the session data with the rebuilt tracker
    data["best_files"] = tracker.to_dict()
    data["best_files_history"] = [h.to_dict() for h in tracker.get_history()]

    # Log what we rebuilt
    best_dict = tracker.get_best_dict()
    if best_dict:
        best_summary = ", ".join(["%s=%s" % (k, os.path.basename(v)) for k, v in best_dict.items() if v])
        print("Rebuilt best_files: %s" % best_summary)
    else:
        print("Rebuilt best_files: (empty)")


def _program_to_stage(program):
    """Map program name to a stage identifier for best_files tracking."""
    program_lower = program.lower() if program else ""

    if "refine" in program_lower:
        return "refined"
    elif "phaser" in program_lower:
        return "phaser_output"
    elif "predict_and_build" in program_lower:
        return "predicted"
    elif "autobuild" in program_lower:
        return "autobuild"
    elif "autosol" in program_lower:
        return "autosol"
    elif "xtriage" in program_lower or "mtriage" in program_lower:
        return "analysis"
    elif "dock_in_map" in program_lower:
        return "docked"
    elif "ligandfit" in program_lower:
        return "ligand_fitted"
    else:
        return "unknown"


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

    # Clear best_files tracking
    if "best_files" in data:
        del data["best_files"]
    if "best_files_history" in data:
        del data["best_files_history"]

    # Clear recovery state
    if "recovery_strategies" in data:
        del data["recovery_strategies"]
    if "recovery_attempts" in data:
        del data["recovery_attempts"]
    if "force_retry_program" in data:
        del data["force_retry_program"]

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

