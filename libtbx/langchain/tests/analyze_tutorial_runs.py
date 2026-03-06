#!/usr/bin/env python3
"""
PHENIX AI Agent — Tutorial Run Analyzer
========================================

Discovers completed agent tutorial runs, extracts structured data from
each session JSON, compares across modes, and produces publication-ready
tables plus raw data for further analysis.

Usage
-----
  python analyze_tutorial_runs.py /path/to/results_root
  python analyze_tutorial_runs.py /path/to/results_root --output report
  python analyze_tutorial_runs.py /path/to/results_root --format csv md json

Directory naming convention
---------------------------
  {tutorial_name}__{mode}     (double underscore preferred)
  {tutorial_name}_{mode}      (single underscore also works)

  where mode is one of:
    rules_only           — deterministic rules selector, no LLM
    llm                  — LLM planning, no thinking agent
    llm_think            — LLM planning + thinking (default)
    llm_think_advanced   — LLM planning + advanced thinking
    llm_think_expert     — LLM planning + expert thinking + plan/gates

  Examples:
    1J4R-ligand__rules_only/
    1J4R-ligand__llm/
    1J4R-ligand__llm_think/
    1J4R-ligand__llm_think_advanced/
    1J4R-ligand__llm_think_expert/

The script searches each directory (recursively) for agent_session.json.
"""

from __future__ import annotations

import argparse
import csv
import glob
import io
import json
import os
import re
import sys
from collections import OrderedDict, defaultdict
from datetime import datetime

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

KNOWN_MODES = [
    "rules_only",
    "llm",
    "llm_think",
    "llm_think_advanced",
    "llm_think_expert",
]

# Canonical mapping: what each mode means internally
MODE_CANONICAL = {
    "rules_only": "rules_only",
    "llm": "llm",
    "llm_think": "llm_think",
    "llm_think_advanced": "llm_think_advanced",
    "llm_think_expert": "llm_think_expert",
}

MODE_LABELS = {
    "rules_only": "Rules",
    "llm": "LLM",
    "llm_think": "LLM+Think",
    "llm_think_advanced": "LLM+Think+Adv",
    "llm_think_expert": "LLM+Think+Expert",
}

# Short labels for tight table columns
MODE_SHORT = {
    "rules_only": "Rules",
    "llm": "LLM",
    "llm_think": "Think",
    "llm_think_advanced": "Adv",
    "llm_think_expert": "Expert",
}

# Metrics where lower is better
LOWER_IS_BETTER = {
    "r_free", "r_work", "clashscore", "bonds_rmsd",
    "angles_rmsd", "ramachandran_outliers",
    "rotamer_outliers", "molprobity_score",
}

# Metrics where higher is better
HIGHER_IS_BETTER = {
    "map_cc", "model_map_cc", "cc_mask", "cc_volume",
    "cc_peaks", "cc_box", "ramachandran_favored",
    "fom", "ligand_cc",
}

METRIC_LABELS = {
    "r_free": "R-free",
    "r_work": "R-work",
    "map_cc": "Map CC",
    "model_map_cc": "Model-Map CC",
    "cc_mask": "CC mask",
    "cc_volume": "CC volume",
    "clashscore": "Clashscore",
    "bonds_rmsd": "Bonds RMSD",
    "angles_rmsd": "Angles RMSD",
    "ramachandran_outliers": "Rama outliers (%)",
    "ramachandran_favored": "Rama favored (%)",
    "rotamer_outliers": "Rotamer outliers (%)",
    "molprobity_score": "MolProbity",
    "fom": "FOM",
    "ligand_cc": "Ligand CC",
    "tfz": "TFZ",
    "llg": "LLG",
}

# ---------------------------------------------------------------------------
#  Discovery
# ---------------------------------------------------------------------------

def find_session_json(directory):
    """Find agent_session.json inside a run directory."""
    p = os.path.join(directory, "agent_session.json")
    if os.path.isfile(p):
        return p
    p = os.path.join(
        directory, "ai_agent_directory",
        "agent_session.json")
    if os.path.isfile(p):
        return p
    matches = glob.glob(
        os.path.join(directory, "**",
                     "agent_session.json"),
        recursive=True)
    if matches:
        return matches[0]
    return None


def parse_directory_name(dirname):
    """Parse a run directory name into (tutorial, mode).

    Supports two naming conventions:
      {tutorial}__{mode}     (double underscore)
      {tutorial}_{mode}      (single underscore)

    Double underscore is tried first (unambiguous).
    Returns (dirname, "unknown") if no mode found.
    """
    base = os.path.basename(dirname.rstrip("/"))

    # Try double-underscore split first
    if "__" in base:
        parts = base.split("__", 1)
        if len(parts) == 2:
            tutorial, mode_str = parts
            if mode_str in KNOWN_MODES:
                return tutorial, mode_str

    # Fall back to single-underscore suffix matching
    # (longest mode first to avoid prefix collisions)
    for mode in sorted(KNOWN_MODES, key=len,
                       reverse=True):
        suffix = "_" + mode
        if base.endswith(suffix):
            tutorial = base[: -len(suffix)]
            return tutorial, mode

    return base, "unknown"


def discover_runs(root_dir):
    """Discover all tutorial run directories."""
    runs = []
    if not os.path.isdir(root_dir):
        print("ERROR: %s is not a directory"
              % root_dir, file=sys.stderr)
        return runs

    for entry in sorted(os.listdir(root_dir)):
        full = os.path.join(root_dir, entry)
        if not os.path.isdir(full):
            continue
        session_path = find_session_json(full)
        if session_path is None:
            continue
        tutorial, mode = parse_directory_name(entry)
        runs.append({
            "tutorial": tutorial,
            "mode": mode,
            "directory": full,
            "dir_name": entry,
            "session_path": session_path,
        })
    return runs


# ---------------------------------------------------------------------------
#  Extraction
# ---------------------------------------------------------------------------

def load_session(path):
    """Load and return the raw session dict."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_run_data(session):
    """Extract all key data from a session dict."""
    cycles = session.get("cycles", [])

    # ── Basic info ────────────────────────────────
    experiment_type = session.get(
        "experiment_type", "unknown")
    resolution = session.get("resolution")
    advice = session.get("project_advice", "")

    # ── Program sequence ──────────────────────────
    program_sequence = []
    for c in cycles:
        prog = c.get("program", "")
        if prog and prog not in [
            None, "unknown", "STOP"
        ]:
            program_sequence.append(prog)

    # ── Per-cycle metrics ─────────────────────────
    cycle_metrics = []
    for c in cycles:
        prog = c.get("program", "")
        if not prog or prog in [
            None, "unknown", "STOP"
        ]:
            continue
        m = c.get("metrics", {})
        result_str = str(c.get("result", ""))
        success = "SUCCESS" in result_str.upper()

        entry = {
            "cycle": c.get("cycle_number", 0),
            "program": prog,
            "success": success,
        }
        # Copy known metrics
        for key in [
            "r_free", "r_work",
            "map_cc", "model_map_cc",
            "cc_mask", "cc_volume",
            "resolution",
            "bonds_rmsd", "angles_rmsd",
            "clashscore",
            "ramachandran_outliers",
            "ramachandran_favored",
            "rotamer_outliers",
            "molprobity_score",
            "tfz", "llg", "fom", "bayes_cc",
            "ligand_cc",
        ]:
            if key in m:
                entry[key] = m[key]

        # Try to pull r_free from result text
        if "r_free" not in entry:
            match = re.search(
                r'R.?[Ff]ree[:\s=]+([0-9.]+)',
                result_str)
            if match:
                try:
                    entry["r_free"] = float(
                        match.group(1))
                except ValueError:
                    pass

        # Expert assessment
        ea = c.get("expert_assessment", {})
        if ea and isinstance(ea, dict):
            if ea.get("analysis"):
                entry["expert_action"] = ea.get(
                    "action", "")
                entry["expert_analysis"] = ea.get(
                    "analysis", "")
                entry["expert_confidence"] = ea.get(
                    "confidence", "")
                entry["expert_guidance"] = ea.get(
                    "guidance", "")

        # Decision reasoning
        decision = c.get("decision", {})
        if isinstance(decision, dict):
            reasoning = decision.get("reasoning", "")
            if reasoning:
                entry["reasoning"] = reasoning

        cycle_metrics.append(entry)

    # ── Final metrics ─────────────────────────────
    final = _extract_final_metrics(cycles, session)

    # ── Counts ────────────────────────────────────
    real_cycles = [
        c for c in cycles
        if c.get("program") not in [
            None, "", "unknown", "STOP"]]
    successful = [
        c for c in real_cycles
        if "SUCCESS" in str(
            c.get("result", "")).upper()]
    failed = [
        c for c in real_cycles
        if "FAILED" in str(
            c.get("result", "")).upper()]

    # ── R-free trajectory ─────────────────────────
    rfree_trajectory = []
    for cm in cycle_metrics:
        if "r_free" in cm:
            rfree_trajectory.append(
                (cm["cycle"], cm["r_free"]))

    # ── Map CC trajectory (cryo-EM) ───────────────
    cc_trajectory = []
    for cm in cycle_metrics:
        for key in ("model_map_cc", "map_cc",
                     "cc_mask"):
            if key in cm:
                cc_trajectory.append(
                    (cm["cycle"], cm[key]))
                break

    # ── Programs used (unique, ordered) ───────────
    programs_used = list(
        OrderedDict.fromkeys(program_sequence))

    # ── Expert assessment summary ─────────────────
    expert_assessments = []
    for cm in cycle_metrics:
        if "expert_analysis" in cm:
            expert_assessments.append({
                "cycle": cm["cycle"],
                "program": cm["program"],
                "action": cm.get(
                    "expert_action", ""),
                "analysis": cm.get(
                    "expert_analysis", ""),
                "confidence": cm.get(
                    "expert_confidence", ""),
                "guidance": cm.get(
                    "expert_guidance", ""),
            })

    # ── Stop info ─────────────────────────────────
    stop_cycle = None
    stop_reasoning = ""
    for c in cycles:
        if c.get("program") == "STOP":
            stop_cycle = c.get("cycle_number")
            stop_reasoning = c.get("reasoning", "")

    # ── Plan info (v114 expert mode) ──────────────
    plan_info = _extract_plan_info(session)

    # ── Structure model (v114) ────────────────────
    sm_info = _extract_structure_model(session)

    # ── Directives ────────────────────────────────
    directives = session.get("directives", {})

    return {
        "experiment_type": experiment_type,
        "resolution": resolution,
        "advice": advice,
        "program_sequence": program_sequence,
        "programs_used": programs_used,
        "cycle_metrics": cycle_metrics,
        "final_metrics": final,
        "total_cycles": len(real_cycles),
        "successful_cycles": len(successful),
        "failed_cycles": len(failed),
        "rfree_trajectory": rfree_trajectory,
        "cc_trajectory": cc_trajectory,
        "expert_assessments": expert_assessments,
        "stop_cycle": stop_cycle,
        "stop_reasoning": stop_reasoning,
        "plan_info": plan_info,
        "structure_model": sm_info,
        "directives": directives,
    }


def _extract_final_metrics(cycles, session=None):
    """Extract final quality metrics.

    Checks StructureModel first (v114), then falls
    back to scanning cycles in reverse.
    """
    metrics = {}

    # Try structure_model first (v114)
    if session:
        sm = session.get("structure_model")
        if sm and isinstance(sm, dict):
            ms = sm.get("model_state", {})
            dc = sm.get("data_characteristics", {})
            for key, src in [
                ("r_free", ms),
                ("r_work", ms),
                ("model_map_cc", ms),
                ("clashscore",
                 ms.get("geometry", {})),
                ("ramachandran_favored",
                 ms.get("geometry", {})),
                ("resolution", dc),
            ]:
                if isinstance(src, dict) and (
                    key in src
                ):
                    try:
                        metrics[key] = float(
                            src[key])
                    except (ValueError, TypeError):
                        pass

    # Fill from cycle history (reverse scan)
    for cycle in reversed(cycles):
        result = str(cycle.get("result", ""))
        cm = cycle.get("metrics", {})

        for key, patterns in [
            ("r_free",
             [r'R.?[Ff]ree[:\s=]+([0-9.]+)']),
            ("r_work",
             [r'R.?[Ww]ork[:\s=]+([0-9.]+)']),
            ("map_cc",
             [r'[Mm]ap.?[Cc][Cc][:\s=]+([0-9.]+)']),
            ("model_map_cc",
             [r'[Mm]odel.?[Mm]ap.?[Cc][Cc]'
              r'[:\s=]+([0-9.]+)']),
            ("cc_mask",
             [r'CC[_ ]?mask\s*[=:]\s*([0-9.]+)']),
            ("clashscore",
             [r'[Cc]lashscore[:\s=]+([0-9.]+)']),
            ("bonds_rmsd",
             [r'[Bb]onds.?[Rr]msd[:\s=]+'
              r'([0-9.]+)']),
            ("angles_rmsd",
             [r'[Aa]ngles.?[Rr]msd[:\s=]+'
              r'([0-9.]+)']),
            ("ramachandran_outliers",
             [r'[Rr]ama(?:chandran)?.?[Oo]utliers?'
              r'[:\s=]+([0-9.]+)']),
            ("ramachandran_favored",
             [r'[Rr]ama(?:chandran)?.?[Ff]avore?d?'
              r'[:\s=]+([0-9.]+)']),
            ("rotamer_outliers",
             [r'[Rr]otamer.?[Oo]utliers?'
              r'[:\s=]+([0-9.]+)']),
            ("molprobity_score",
             [r'[Mm]ol[Pp]robity.?[Ss]core'
              r'[:\s=]+([0-9.]+)']),
            ("fom",
             [r'[Ff][Oo][Mm][:\s=]+([0-9.]+)']),
            ("ligand_cc",
             [r'[Ll]igand.?[Cc][Cc][:\s=]+'
              r'([0-9.]+)']),
        ]:
            if key in metrics:
                continue
            if key in cm:
                try:
                    metrics[key] = float(cm[key])
                except (ValueError, TypeError):
                    pass
                continue
            for pat in patterns:
                match = re.search(pat, result)
                if match:
                    try:
                        metrics[key] = float(
                            match.group(1))
                    except ValueError:
                        pass
                    break
    return metrics


def _extract_plan_info(session):
    """Extract plan phase information (v114)."""
    plan = session.get("plan")
    if not plan or not isinstance(plan, dict):
        return None

    phases = plan.get("phases", [])
    if not phases:
        return None

    phase_list = []
    for p in phases:
        if not isinstance(p, dict):
            continue
        phase_list.append({
            "id": p.get("id", ""),
            "status": p.get("status", "pending"),
            "description": p.get("description", ""),
            "programs": p.get("programs", []),
            "max_cycles": p.get("max_cycles", 0),
            "cycles_used": p.get("cycles_used", 0),
        })

    return {
        "goal": plan.get("goal", ""),
        "template_id": plan.get("template_id", ""),
        "phases": phase_list,
        "n_complete": sum(
            1 for p in phase_list
            if p["status"] == "complete"),
        "n_total": len(phase_list),
    }


def _extract_structure_model(session):
    """Extract structure model summary (v114)."""
    sm = session.get("structure_model")
    if not sm or not isinstance(sm, dict):
        return None
    ms = sm.get("model_state", {})
    dc = sm.get("data_characteristics", {})
    geom = ms.get("geometry", {})
    chains = ms.get("chains", [])
    ligands = ms.get("ligands", [])

    return {
        "r_free": ms.get("r_free"),
        "r_work": ms.get("r_work"),
        "model_map_cc": ms.get("model_map_cc"),
        "resolution": dc.get("resolution"),
        "space_group": dc.get("space_group"),
        "clashscore": geom.get("clashscore"),
        "rama_favored": geom.get(
            "rama_favored",
            geom.get("ramachandran_favored")),
        "n_chains": len(chains),
        "n_ligands": len(ligands),
        "n_waters": ms.get("waters", 0),
        "mr_tfz": dc.get("mr_tfz"),
        "mr_llg": dc.get("mr_llg"),
        "phasing_fom": dc.get("phasing_fom"),
    }


# ---------------------------------------------------------------------------
#  Grouping
# ---------------------------------------------------------------------------

def group_by_tutorial(runs_data):
    """Group extracted run data by tutorial name."""
    groups = defaultdict(dict)
    for info, data in runs_data:
        groups[info["tutorial"]][info["mode"]] = (
            info, data)
    return OrderedDict(sorted(groups.items()))


def _modes_present(groups):
    """Return sorted list of modes actually present."""
    present = set()
    for modes in groups.values():
        present.update(modes.keys())
    return [m for m in KNOWN_MODES if m in present] + (
        ["unknown"] if "unknown" in present else [])


# ---------------------------------------------------------------------------
#  Output: CSV
# ---------------------------------------------------------------------------

def generate_summary_csv(groups):
    """Generate summary CSV."""
    buf = io.StringIO()
    w = csv.writer(buf)
    modes = _modes_present(groups)

    w.writerow([
        "Tutorial", "Mode", "Mode_Label",
        "Experiment", "Resolution",
        "Total_Cycles", "Successful", "Failed",
        "R_free", "R_work",
        "Map_CC", "Model_Map_CC", "CC_mask",
        "Clashscore", "Bonds_RMSD", "Angles_RMSD",
        "Rama_Outliers", "Rama_Favored",
        "MolProbity", "FOM", "Ligand_CC",
        "Programs_Used", "Expert_Assessments",
        "Plan_Template", "Plan_Phases_Complete",
    ])

    for tutorial, mode_dict in groups.items():
        for mode in modes + ["unknown"]:
            if mode not in mode_dict:
                continue
            info, data = mode_dict[mode]
            fm = data["final_metrics"]
            pi = data.get("plan_info") or {}
            w.writerow([
                tutorial,
                mode,
                MODE_LABELS.get(mode, mode),
                data["experiment_type"],
                _fmt(data["resolution"], ".2f"),
                data["total_cycles"],
                data["successful_cycles"],
                data["failed_cycles"],
                _fmt(fm.get("r_free"), ".4f"),
                _fmt(fm.get("r_work"), ".4f"),
                _fmt(fm.get("map_cc"), ".3f"),
                _fmt(fm.get("model_map_cc"), ".3f"),
                _fmt(fm.get("cc_mask"), ".3f"),
                _fmt(fm.get("clashscore"), ".1f"),
                _fmt(fm.get("bonds_rmsd"), ".3f"),
                _fmt(fm.get("angles_rmsd"), ".2f"),
                _fmt(fm.get(
                    "ramachandran_outliers"), ".2f"),
                _fmt(fm.get(
                    "ramachandran_favored"), ".1f"),
                _fmt(fm.get(
                    "molprobity_score"), ".2f"),
                _fmt(fm.get("fom"), ".3f"),
                _fmt(fm.get("ligand_cc"), ".3f"),
                " > ".join(data["programs_used"]),
                len(data["expert_assessments"]),
                pi.get("template_id", ""),
                "%d/%d" % (
                    pi.get("n_complete", 0),
                    pi.get("n_total", 0))
                if pi else "",
            ])

    return buf.getvalue()


def generate_trajectories_csv(groups):
    """Generate per-cycle metrics CSV for plotting."""
    buf = io.StringIO()
    w = csv.writer(buf)
    modes = _modes_present(groups)

    w.writerow([
        "Tutorial", "Mode", "Mode_Label",
        "Cycle", "Program", "Success",
        "R_free", "R_work",
        "Map_CC", "Model_Map_CC", "CC_mask",
        "Expert_Action", "Expert_Confidence",
    ])

    for tutorial, mode_dict in groups.items():
        for mode in modes + ["unknown"]:
            if mode not in mode_dict:
                continue
            _, data = mode_dict[mode]
            for cm in data["cycle_metrics"]:
                w.writerow([
                    tutorial,
                    mode,
                    MODE_LABELS.get(mode, mode),
                    cm["cycle"],
                    cm["program"],
                    cm["success"],
                    _fmt(cm.get("r_free"), ".4f"),
                    _fmt(cm.get("r_work"), ".4f"),
                    _fmt(cm.get("map_cc"), ".3f"),
                    _fmt(cm.get(
                        "model_map_cc"), ".3f"),
                    _fmt(cm.get("cc_mask"), ".3f"),
                    cm.get("expert_action", ""),
                    cm.get("expert_confidence", ""),
                ])

    return buf.getvalue()


# ---------------------------------------------------------------------------
#  Output: Markdown
# ---------------------------------------------------------------------------

def generate_markdown_report(groups):
    """Generate a publication-oriented Markdown report."""
    lines = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    modes = _modes_present(groups)

    lines.append(
        "# PHENIX AI Agent — Tutorial Evaluation Report")
    lines.append("")
    lines.append("Generated: %s" % now)
    lines.append("")
    lines.append("Tutorials: %d" % len(groups))
    mode_strs = [MODE_LABELS.get(m, m) for m in modes]
    lines.append("Modes compared: %s"
                 % ", ".join(mode_strs))
    lines.append("")

    # ── 1. Summary table ──────────────────────────
    lines.append("## 1. Summary Comparison")
    lines.append("")
    _add_summary_table(lines, groups, modes)
    lines.append("")

    # ── 2. Per-tutorial detail ────────────────────
    lines.append("## 2. Per-Tutorial Detail")
    lines.append("")
    for tutorial, mode_dict in groups.items():
        lines.append("### %s" % tutorial)
        lines.append("")

        any_data = next(iter(mode_dict.values()))[1]
        lines.append(
            "**Experiment:** %s  "
            "**Resolution:** %s A"
            % (any_data["experiment_type"],
               _fmt(any_data["resolution"], ".2f")))
        if any_data["advice"]:
            adv = any_data["advice"]
            if len(adv) > 120:
                adv = adv[:117] + "..."
            lines.append("**Advice:** %s" % adv)
        lines.append("")

        # Workflow comparison
        _add_workflow_comparison(
            lines, tutorial, mode_dict, modes)
        lines.append("")

        # Metrics comparison
        _add_metrics_comparison(
            lines, mode_dict, modes)
        lines.append("")

        # Plan info (expert mode)
        for mode in ("llm_think_expert",
                     "llm_think_advanced"):
            if mode in mode_dict:
                _, mdata = mode_dict[mode]
                pi = mdata.get("plan_info")
                if pi and pi.get("phases"):
                    lines.append(
                        "**Plan (%s):** %s"
                        % (MODE_LABELS[mode],
                           pi.get("goal", "")))
                    lines.append("")
                    for p in pi["phases"]:
                        icon = {
                            "complete": "+",
                            "active": ">",
                            "skipped": "-",
                            "failed": "x",
                        }.get(p["status"], "o")
                        lines.append(
                            "- [%s] %s (%d/%d cycles)"
                            % (icon,
                               p.get("description",
                                     p["id"]),
                               p["cycles_used"],
                               p["max_cycles"]))
                    lines.append("")
                    break  # Show plan from first available

        # Expert assessments
        for mode in ("llm_think_expert",
                     "llm_think_advanced",
                     "llm_think", "llm"):
            if mode in mode_dict:
                _, mdata = mode_dict[mode]
                ea_list = mdata["expert_assessments"]
                if ea_list:
                    lines.append(
                        "**Expert Observations "
                        "(%s):**"
                        % MODE_LABELS[mode])
                    lines.append("")
                    for ea in ea_list:
                        action = ea.get(
                            "action", "continue")
                        conf = ea.get(
                            "confidence", "")
                        analysis = ea.get(
                            "analysis", "")
                        if len(analysis) > 150:
                            analysis = (
                                analysis[:147]
                                + "...")
                        lines.append(
                            "- Cycle %d (%s): "
                            "[%s/%s] %s"
                            % (ea["cycle"],
                               ea["program"],
                               action, conf,
                               analysis))
                    lines.append("")
                    break  # Only show for first mode

    # ── 3. Aggregate statistics ───────────────────
    lines.append("## 3. Aggregate Statistics")
    lines.append("")
    _add_aggregate_stats(lines, groups, modes)
    lines.append("")

    return "\n".join(lines)


def _add_summary_table(lines, groups, modes):
    """Add the main summary comparison table."""
    header = ["Tutorial", "Exp", "Res"]
    for mode in modes:
        label = MODE_SHORT.get(mode, mode)
        header.extend([
            "%s Cyc" % label,
            "%s Metric" % label])
    lines.append("| " + " | ".join(header) + " |")
    lines.append(
        "| " + " | ".join(
            ["---"] * len(header)) + " |")

    for tutorial, mode_dict in groups.items():
        any_data = next(iter(mode_dict.values()))[1]
        exp = any_data["experiment_type"]
        is_em = exp in ("cryoem", "Cryo-EM")
        row = [
            tutorial,
            "EM" if is_em else "XR",
            _fmt(any_data["resolution"], ".1f"),
        ]
        for mode in modes:
            if mode in mode_dict:
                _, data = mode_dict[mode]
                fm = data["final_metrics"]
                row.append(str(data["total_cycles"]))
                if is_em:
                    val = (fm.get("model_map_cc")
                           or fm.get("cc_mask")
                           or fm.get("map_cc"))
                    s = _fmt(val, ".3f")
                    if s != "\u2014":
                        s += " CC"
                    row.append(s)
                else:
                    row.append(
                        _fmt(fm.get("r_free"),
                             ".4f"))
            else:
                row.extend(["\u2014", "\u2014"])
        lines.append(
            "| " + " | ".join(row) + " |")


def _add_workflow_comparison(lines, tutorial,
                             mode_dict, modes):
    """Add cycle-by-cycle workflow comparison."""
    max_cyc = 0
    for _, data in mode_dict.values():
        if data["cycle_metrics"]:
            max_cyc = max(
                max_cyc,
                max(cm["cycle"]
                    for cm in data["cycle_metrics"]))

    if max_cyc == 0:
        lines.append("*No cycles recorded.*")
        return

    mode_order = [m for m in modes if m in mode_dict]
    header = ["Cycle"]
    for mode in mode_order:
        label = MODE_SHORT.get(mode, mode)
        header.extend([
            "%s Prog" % label,
            "%s Metric" % label])
    lines.append("| " + " | ".join(header) + " |")
    lines.append(
        "| " + " | ".join(
            ["---"] * len(header)) + " |")

    mode_cycles = {}
    for mode in mode_order:
        _, data = mode_dict[mode]
        mode_cycles[mode] = {
            cm["cycle"]: cm
            for cm in data["cycle_metrics"]}

    for cyc in range(1, max_cyc + 1):
        row = [str(cyc)]
        for mode in mode_order:
            cm = mode_cycles[mode].get(cyc)
            if cm:
                prog = cm["program"].replace(
                    "phenix.", "")
                metric = ""
                for key, fmt in [
                    ("r_free", ".4f"),
                    ("model_map_cc", ".3f"),
                    ("map_cc", ".3f"),
                    ("cc_mask", ".3f"),
                    ("tfz", ".1f"),
                    ("fom", ".2f"),
                    ("ligand_cc", ".2f"),
                ]:
                    if key in cm:
                        metric = _fmt(cm[key], fmt)
                        if key == "tfz":
                            metric = "TFZ=" + metric
                        elif key == "fom":
                            metric = "FOM=" + metric
                        break
                if not cm["success"]:
                    prog = "~~%s~~" % prog
                    metric = "FAIL"
                row.extend([prog, metric])
            else:
                row.extend(["\u2014", "\u2014"])
        lines.append(
            "| " + " | ".join(row) + " |")


def _add_metrics_comparison(lines, mode_dict, modes):
    """Add final metrics comparison."""
    all_keys = set()
    for _, data in mode_dict.values():
        all_keys.update(data["final_metrics"].keys())

    if not all_keys:
        return

    key_order = [k for k in [
        "r_free", "r_work",
        "map_cc", "model_map_cc", "cc_mask",
        "cc_volume",
        "clashscore", "bonds_rmsd", "angles_rmsd",
        "ramachandran_outliers",
        "ramachandran_favored",
        "rotamer_outliers", "molprobity_score",
        "fom", "ligand_cc",
    ] if k in all_keys]

    mode_order = [m for m in modes if m in mode_dict]
    header = ["Metric"] + [
        MODE_SHORT.get(m, m) for m in mode_order]
    lines.append("| " + " | ".join(header) + " |")
    lines.append(
        "| " + " | ".join(
            ["---"] * len(header)) + " |")

    for key in key_order:
        row = [METRIC_LABELS.get(key, key)]
        values = []
        for mode in mode_order:
            _, data = mode_dict[mode]
            v = data["final_metrics"].get(key)
            values.append(v)

        best_idx = _find_best(key, values)

        for i, v in enumerate(values):
            if key.startswith("r_"):
                s = _fmt(v, ".4f")
            elif "cc" in key or key == "fom":
                s = _fmt(v, ".3f")
            else:
                s = _fmt(v, ".2f")
            if (i == best_idx
                    and len(values) > 1
                    and s != "\u2014"):
                s = "**%s**" % s
            row.append(s)

        lines.append(
            "| " + " | ".join(row) + " |")


def _add_aggregate_stats(lines, groups, modes):
    """Add aggregate statistics across tutorials."""
    # Count best-metric wins
    wins = defaultdict(int)
    totals = 0

    for tutorial, mode_dict in groups.items():
        any_data = next(iter(mode_dict.values()))[1]
        is_em = any_data["experiment_type"] in (
            "cryoem", "Cryo-EM")
        if is_em:
            for key in ("model_map_cc", "cc_mask",
                        "map_cc"):
                if any(
                    mode_dict[m][1][
                        "final_metrics"].get(key)
                    for m in mode_dict
                    if m in modes
                ):
                    break
            else:
                key = "model_map_cc"
        else:
            key = "r_free"

        values = {}
        for mode in modes:
            if mode in mode_dict:
                v = mode_dict[mode][1][
                    "final_metrics"].get(key)
                if v is not None:
                    values[mode] = v

        if not values:
            continue
        totals += 1

        best_val = None
        best_mode = None
        for mode, v in values.items():
            if best_val is None:
                best_val = v
                best_mode = mode
            elif (key in LOWER_IS_BETTER
                  and v < best_val):
                best_val = v
                best_mode = mode
            elif (key in HIGHER_IS_BETTER
                  and v > best_val):
                best_val = v
                best_mode = mode
        if best_mode:
            wins[best_mode] += 1

    if totals == 0:
        lines.append(
            "*No comparable metrics found.*")
        return

    lines.append(
        "**Best final metric by mode** "
        "(across %d tutorials):" % totals)
    lines.append("")
    for mode in modes:
        label = MODE_LABELS.get(mode, mode)
        count = wins.get(mode, 0)
        lines.append(
            "- %s: %d/%d" % (label, count, totals))

    # Average cycle count
    lines.append("")
    lines.append("**Average cycle count by mode:**")
    lines.append("")
    for mode in modes:
        counts = []
        for mode_dict in groups.values():
            if mode in mode_dict:
                counts.append(
                    mode_dict[mode][1][
                        "total_cycles"])
        if counts:
            avg = sum(counts) / len(counts)
            label = MODE_LABELS.get(mode, mode)
            lines.append(
                "- %s: %.1f cycles (n=%d)"
                % (label, avg, len(counts)))

    # Expert assessment frequency
    lines.append("")
    lines.append(
        "**Expert assessments per run (average):**")
    lines.append("")
    for mode in modes:
        ea_counts = []
        for mode_dict in groups.values():
            if mode in mode_dict:
                ea_counts.append(
                    len(mode_dict[mode][1][
                        "expert_assessments"]))
        if ea_counts:
            avg = sum(ea_counts) / len(ea_counts)
            label = MODE_LABELS.get(mode, mode)
            lines.append(
                "- %s: %.1f (n=%d)"
                % (label, avg, len(ea_counts)))


# ---------------------------------------------------------------------------
#  Output: JSON
# ---------------------------------------------------------------------------

def generate_json_report(groups):
    """Generate complete JSON export."""
    output = {
        "generated": datetime.now().isoformat(),
        "tutorials": {},
    }

    for tutorial, mode_dict in groups.items():
        tut_data = {}
        for mode, (info, data) in mode_dict.items():
            tut_data[mode] = {
                "directory": info["directory"],
                "session_path": info["session_path"],
                **data,
            }
        output["tutorials"][tutorial] = tut_data

    return json.dumps(output, indent=2, default=str)


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _fmt(value, fmt_spec):
    """Format a numeric value, returning dash for None."""
    if value is None:
        return "\u2014"
    try:
        return format(float(value), fmt_spec)
    except (ValueError, TypeError):
        return str(value)


def _find_best(metric_key, values):
    """Return index of the best value, or None."""
    valid = [
        (i, v) for i, v in enumerate(values)
        if v is not None]
    if len(valid) < 2:
        return None
    if metric_key in LOWER_IS_BETTER:
        return min(valid, key=lambda x: x[1])[0]
    elif metric_key in HIGHER_IS_BETTER:
        return max(valid, key=lambda x: x[1])[0]
    return None


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze PHENIX AI Agent tutorial runs"),
        formatter_class=(
            argparse.RawDescriptionHelpFormatter),
        epilog=__doc__,
    )
    parser.add_argument(
        "results_dir",
        help="Root directory containing run dirs")
    parser.add_argument(
        "--output", "-o", default="report",
        help="Output file prefix (default: report)")
    parser.add_argument(
        "--format", "-f", nargs="+",
        default=["csv", "md", "json"],
        choices=["csv", "md", "json"],
        help="Output formats (default: csv md json)")
    args = parser.parse_args()

    # ── Discover runs ─────────────────────────────
    print("Scanning: %s" % args.results_dir)
    runs = discover_runs(args.results_dir)
    if not runs:
        print("No tutorial runs found.",
              file=sys.stderr)
        sys.exit(1)

    print("Found %d run(s):" % len(runs))
    for r in runs:
        label = MODE_LABELS.get(r["mode"], r["mode"])
        print("  %s  ->  tutorial=%s  mode=%s (%s)"
              % (r["dir_name"], r["tutorial"],
                 r["mode"], label))

    # ── Extract data ──────────────────────────────
    runs_data = []
    for r in runs:
        try:
            session = load_session(r["session_path"])
            data = extract_run_data(session)
            runs_data.append((r, data))
            n = data["total_cycles"]
            fm = data["final_metrics"]
            rf = fm.get("r_free")
            cc = (fm.get("model_map_cc")
                  or fm.get("cc_mask")
                  or fm.get("map_cc"))
            metric_str = ""
            if rf is not None:
                metric_str = "R-free=%.4f" % rf
            elif cc is not None:
                metric_str = "CC=%.3f" % cc
            pi = data.get("plan_info")
            plan_str = ""
            if pi:
                plan_str = " plan=%d/%d" % (
                    pi["n_complete"], pi["n_total"])
            print("  + %s: %d cycles, %s%s"
                  % (r["dir_name"], n, metric_str,
                     plan_str))
        except Exception as e:
            print("  x %s: ERROR - %s"
                  % (r["dir_name"], e),
                  file=sys.stderr)

    if not runs_data:
        print("No valid sessions found.",
              file=sys.stderr)
        sys.exit(1)

    # ── Group by tutorial ─────────────────────────
    groups = group_by_tutorial(runs_data)
    print("\nGrouped into %d tutorial(s):"
          % len(groups))
    for tut, mode_dict in groups.items():
        mode_list = ", ".join(
            MODE_LABELS.get(m, m)
            for m in sorted(mode_dict.keys()))
        print("  %s: %s" % (tut, mode_list))

    # ── Generate outputs ──────────────────────────
    output_dir = os.path.dirname(args.output) or "."
    prefix = os.path.basename(args.output)

    if "csv" in args.format:
        path = os.path.join(
            output_dir, "%s_summary.csv" % prefix)
        with open(path, "w", encoding="utf-8") as f:
            f.write(generate_summary_csv(groups))
        print("\nWrote: %s" % path)

        path = os.path.join(
            output_dir,
            "%s_trajectories.csv" % prefix)
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                generate_trajectories_csv(groups))
        print("Wrote: %s" % path)

    if "md" in args.format:
        path = os.path.join(
            output_dir, "%s.md" % prefix)
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                generate_markdown_report(groups))
        print("Wrote: %s" % path)

    if "json" in args.format:
        path = os.path.join(
            output_dir, "%s.json" % prefix)
        with open(path, "w", encoding="utf-8") as f:
            f.write(generate_json_report(groups))
        print("Wrote: %s" % path)

    print("\nDone.")


if __name__ == "__main__":
    main()
