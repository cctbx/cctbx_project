#!/usr/bin/env python
"""
Generate human-readable documentation of the PHENIX AI Agent logic.

This script extracts and formats all decision logic from:
- decision_config.json (tiered decision architecture)
- workflow_state.py (state machine and transitions)
- command_templates.json (program configurations)
- metrics_analyzer.py (stop conditions)
- prompts_hybrid.py (LLM guidance)

Output: A clear markdown document describing all agent behavior.

Usage:
    python generate_logic_doc.py > AGENT_LOGIC.md
    python generate_logic_doc.py --format text  # Plain text output
"""

import json
import os
import re
import sys
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_decision_config():
    """Load decision_config.json."""
    path = os.path.join(SCRIPT_DIR, "decision_config.json")
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def load_command_templates():
    """Load program definitions from programs.yaml."""
    # First try the knowledge directory relative to SCRIPT_DIR
    yaml_path = os.path.join(os.path.dirname(SCRIPT_DIR), "knowledge", "programs.yaml")
    if not os.path.exists(yaml_path):
        # Try sibling directory
        yaml_path = os.path.join(SCRIPT_DIR, "..", "knowledge", "programs.yaml")

    try:
        import yaml
        with open(yaml_path) as f:
            programs = yaml.safe_load(f)

        # Convert YAML format to the old templates format for compatibility
        templates = {}
        for prog_name, prog_def in programs.items():
            if isinstance(prog_def, dict):
                templates[prog_name] = {
                    "description": prog_def.get("description", ""),
                    "inputs": prog_def.get("inputs", {}),
                    "outputs": prog_def.get("outputs", {}),
                    "category": prog_def.get("category", "other"),
                    "experiment_type": prog_def.get("experiment_types", prog_def.get("experiment_type", [])),
                }
        return templates
    except FileNotFoundError:
        print(f"Warning: Could not find programs.yaml at {yaml_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Warning: Error loading programs.yaml: {e}", file=sys.stderr)
        return {}


def extract_workflow_states():
    """Extract workflow states and transitions from workflow_state.py."""
    path = os.path.join(SCRIPT_DIR, "workflow_state.py")
    with open(path) as f:
        content = f.read()

    states = []

    # Find all state returns in _detect_xray_state
    xray_match = re.search(
        r'def _detect_xray_state\(.*?\):(.*?)(?=\ndef |\Z)',
        content,
        re.DOTALL
    )
    if xray_match:
        xray_code = xray_match.group(1)
        # Extract state definitions
        state_patterns = re.findall(
            r'"state":\s*"([^"]+)".*?"valid_programs":\s*\[([^\]]+)\].*?"reason":\s*"([^"]+)"',
            xray_code,
            re.DOTALL
        )
        for state, programs, reason in state_patterns:
            progs = [p.strip().strip('"\'') for p in programs.split(',')]
            states.append({
                "name": state,
                "type": "xray",
                "valid_programs": progs,
                "reason": reason.split(',')[0]  # First part of reason
            })

    # Find all state returns in _detect_cryoem_state
    cryoem_match = re.search(
        r'def _detect_cryoem_state\(.*?\):(.*?)(?=\ndef |\Z)',
        content,
        re.DOTALL
    )
    if cryoem_match:
        cryoem_code = cryoem_match.group(1)
        state_patterns = re.findall(
            r'"state":\s*"([^"]+)".*?"valid_programs":\s*\[([^\]]+)\].*?"reason":\s*"([^"]+)"',
            cryoem_code,
            re.DOTALL
        )
        for state, programs, reason in state_patterns:
            progs = [p.strip().strip('"\'') for p in programs.split(',')]
            states.append({
                "name": state,
                "type": "cryoem",
                "valid_programs": progs,
                "reason": reason.split(',')[0]
            })

    return states


def extract_stop_conditions():
    """Extract stop conditions from metrics_analyzer.py."""
    path = os.path.join(SCRIPT_DIR, "metrics_analyzer.py")
    with open(path) as f:
        content = f.read()

    conditions = []

    # Find success threshold logic
    if "dynamic_target" in content:
        conditions.append({
            "name": "X-ray Success",
            "description": "R-free below dynamic target (resolution/10, bounded 0.20-0.30)",
            "formula": "R-free < max(0.20, min(0.30, resolution/10)) - 0.02"
        })

    # Find plateau detection
    if "PLATEAU" in content:
        conditions.append({
            "name": "Plateau Detection",
            "description": "Less than 0.5% improvement for 2+ consecutive cycles",
            "formula": "improvement < 0.5% for last 2 cycles"
        })

    # Find excessive refinement
    if "EXCESSIVE" in content:
        conditions.append({
            "name": "Excessive Refinement",
            "description": "Too many consecutive refinement cycles",
            "formula": "consecutive_refines >= 5"
        })

    # Find cryo-EM success
    if "map_cc" in content and "0.75" in content:
        conditions.append({
            "name": "Cryo-EM Success",
            "description": "Map-model correlation above threshold",
            "formula": "map_cc > 0.75"
        })

    return conditions


def extract_file_categorization():
    """Extract file categorization rules from workflow_state.py."""
    path = os.path.join(SCRIPT_DIR, "workflow_state.py")
    with open(path) as f:
        content = f.read()

    rules = []

    # Find _categorize_files function
    match = re.search(
        r'def _categorize_files\(.*?\):(.*?)(?=\ndef |\Z)',
        content,
        re.DOTALL
    )
    if match:
        func_code = match.group(1)

        # Extract patterns
        patterns = [
            ("MTZ files", ".mtz", "X-ray data"),
            ("Sequence files", ".fa, .fasta, .seq", "Protein sequence"),
            ("Map files", ".mrc, .ccp4, .map", "Cryo-EM maps"),
            ("PDB files", ".pdb", "Model coordinates"),
        ]

        for name, ext, desc in patterns:
            rules.append({"name": name, "extensions": ext, "description": desc})

        # PDB subcategories
        if "'phaser'" in func_code:
            rules.append({"name": "Phaser output", "pattern": "'phaser' in filename", "description": "MR solution"})
        if "'refine'" in func_code and "'real_space'" in func_code:
            rules.append({"name": "Refined model", "pattern": "'refine' in filename (not 'real_space')", "description": "X-ray refined"})
        if "'predict'" in func_code:
            rules.append({"name": "Predicted model", "pattern": "'predict' or 'alphafold' in filename", "description": "AlphaFold prediction"})
        if "'autobuild'" in func_code or "'buccaneer'" in func_code:
            rules.append({"name": "Autobuild output", "pattern": "'autobuild', 'buccaneer', 'build' in filename", "description": "Model building output"})
        if "'ligand'" in func_code:
            rules.append({"name": "Ligand CIF", "pattern": ".cif without 'refine'", "description": "Ligand restraints"})

    return rules


def extract_key_thresholds():
    """Extract important threshold values from decision_config.json or code."""
    thresholds = []

    # Try to get from decision_config.json first
    config = load_decision_config()
    if config:
        xray = config.get("thresholds", {}).get("xray", {})
        cryoem = config.get("thresholds", {}).get("cryoem", {})

        # Resolution-dependent thresholds
        res_dep = xray.get("resolution_dependent", {})
        for bin_name, bin_config in res_dep.items():
            if bin_name == "default":
                continue
            range_info = bin_config.get("range", {})
            range_str = ""
            if "max" in range_info and "min" not in range_info:
                range_str = "< %.1fÅ" % range_info["max"]
            elif "min" in range_info and "max" in range_info:
                range_str = "%.1f-%.1fÅ" % (range_info["min"], range_info["max"])
            elif "min" in range_info:
                range_str = "> %.1fÅ" % range_info["min"]

            thresholds.append({
                "name": "%s (%s)" % (bin_name.replace("_", " ").title(), range_str),
                "value": "autobuild=%.2f, good=%.2f, ligandfit=%.2f" % (
                    bin_config.get("autobuild_rfree", 0),
                    bin_config.get("good_model_rfree", 0),
                    bin_config.get("ligandfit_rfree", 0)
                ),
                "description": "Resolution-dependent R-free thresholds"
            })

        # Other thresholds
        if "plateau_improvement_threshold" in xray:
            thresholds.append({
                "name": "Plateau detection",
                "value": "< %.1f%% for %d cycles" % (
                    xray["plateau_improvement_threshold"] * 100,
                    xray.get("plateau_cycles_required", 3)
                ),
                "description": "Trigger plateau warning"
            })

        if "excessive_refine_cycles" in xray:
            thresholds.append({
                "name": "Excessive refinement",
                "value": ">= %d cycles" % xray["excessive_refine_cycles"],
                "description": "Too many consecutive refinement cycles"
            })

        if "success_cc" in cryoem:
            thresholds.append({
                "name": "Cryo-EM success",
                "value": "CC > %.2f" % cryoem["success_cc"],
                "description": "Map-model correlation target"
            })

        return thresholds

    # Fallback to code extraction
    ws_path = os.path.join(SCRIPT_DIR, "workflow_state.py")
    with open(ws_path) as f:
        ws_content = f.read()

    if "0.35" in ws_content:
        thresholds.append({
            "name": "Autobuild threshold",
            "value": "R-free > 0.35 (medium res)",
            "description": "Offer autobuild when model quality is poor"
        })

    return thresholds


def generate_tiered_docs(config):
    """Generate documentation for the tiered decision architecture."""
    lines = []

    if not config:
        return lines

    tiers = config.get("tiers", {})

    # === TIER 1: HARD CONSTRAINTS ===
    lines.append("## Tiered Decision Architecture")
    lines.append("")
    lines.append("The agent uses a three-tier decision system:")
    lines.append("")
    lines.append("| Tier | Type | Description | Override? |")
    lines.append("|------|------|-------------|-----------|")
    lines.append("| 1 | Hard Constraints | Must be followed | No |")
    lines.append("| 2 | Strong Defaults | Applied automatically | Yes (with warning) |")
    lines.append("| 3 | Soft Guidance | Suggestions only | Yes (no warning) |")
    lines.append("")

    # Tier 1
    tier1 = tiers.get("tier1_hard_constraints", {})
    lines.append("### Tier 1: Hard Constraints")
    lines.append("")
    lines.append(tier1.get("description", "Cannot be overridden"))
    lines.append("")

    rules = tier1.get("rules", {})
    if rules:
        lines.append("| Rule | Description | Applies To |")
        lines.append("|------|-------------|------------|")
        for rule_name, rule_info in rules.items():
            desc = rule_info.get("description", rule_name)
            applies = rule_info.get("applies_to", "all")
            lines.append("| %s | %s | %s |" % (rule_name, desc, applies))
        lines.append("")

    # Tier 2
    tier2 = tiers.get("tier2_strong_defaults", {})
    lines.append("### Tier 2: Strong Defaults")
    lines.append("")
    lines.append(tier2.get("description", "Applied automatically; LLM can override"))
    lines.append("")

    lines.append("| Default | Condition | Value | Override Warning |")
    lines.append("|---------|-----------|-------|------------------|")
    for rule_name, rule_info in tier2.items():
        if not isinstance(rule_info, dict) or "applies_to" not in rule_info:
            continue
        condition = rule_info.get("condition", str(rule_info.get("conditions", "")))
        if len(condition) > 40:
            condition = condition[:40] + "..."
        value = rule_info.get("default_value", "")
        warning = rule_info.get("override_warning", "")
        if len(warning) > 30:
            warning = warning[:30] + "..."
        lines.append("| %s | %s | %s | %s |" % (rule_name, condition, value, warning))
    lines.append("")

    # Tier 3
    tier3 = tiers.get("tier3_soft_guidance", {})
    lines.append("### Tier 3: Soft Guidance")
    lines.append("")
    lines.append(tier3.get("description", "Suggestions only"))
    lines.append("")

    suggestions = tier3.get("suggestions", [])
    if suggestions:
        for item in suggestions:
            condition = item.get("condition", "")
            suggestion = item.get("suggestion", "")
            lines.append("- **When** `%s`: %s" % (condition, suggestion))
        lines.append("")

    # === PROGRAM RANKINGS ===
    rankings = config.get("program_rankings", {})
    if rankings:
        lines.append("### Program Rankings by State")
        lines.append("")
        lines.append("The graph recommends programs in priority order:")
        lines.append("")

        for state_name, state_info in rankings.items():
            desc = state_info.get("description", "")
            lines.append("**%s** - %s" % (state_name, desc))
            lines.append("")

            state_rankings = state_info.get("rankings", [])
            if state_rankings:
                lines.append("| Priority | Program | Condition |")
                lines.append("|----------|---------|-----------|")
                for item in state_rankings:
                    prog = item.get("program", "")
                    cond = item.get("condition", "always")
                    if len(cond) > 50:
                        cond = cond[:50] + "..."
                    prio = item.get("priority", "?")
                    lines.append("| %s | %s | %s |" % (prio, prog, cond))
                lines.append("")

    return lines


def generate_markdown(templates, states, stop_conditions, file_rules, thresholds, config=None):
    """Generate markdown documentation."""
    lines = []

    lines.append("# PHENIX AI Agent Logic Documentation")
    lines.append("")
    lines.append("This document describes all decision logic used by the AI agent.")
    lines.append("Generated automatically from source code and configuration files.")
    lines.append("")

    # Table of Contents
    lines.append("## Table of Contents")
    lines.append("1. [Tiered Decision Architecture](#tiered-decision-architecture)")
    lines.append("2. [Workflow States](#workflow-states)")
    lines.append("3. [Program Templates](#program-templates)")
    lines.append("4. [Stop Conditions](#stop-conditions)")
    lines.append("5. [File Categorization](#file-categorization)")
    lines.append("6. [Key Thresholds](#key-thresholds)")
    lines.append("")

    # Tiered Decision Architecture (from config)
    lines.append("---")
    tiered_docs = generate_tiered_docs(config)
    lines.extend(tiered_docs)

    # Workflow States
    lines.append("---")
    lines.append("## Workflow States")
    lines.append("")

    # X-ray states
    lines.append("### X-ray Crystallography Workflow")
    lines.append("")
    lines.append("```")
    lines.append("xray_initial → xtriage → xray_analyzed")
    lines.append("                              ↓")
    lines.append("              [if sequence] predict_and_build → xray_has_prediction")
    lines.append("              [if model]    phaser ─────────────────┐")
    lines.append("                                                    ↓")
    lines.append("                     xray_has_prediction → process_predicted_model")
    lines.append("                                                    ↓")
    lines.append("                                    xray_model_processed → phaser")
    lines.append("                                                              ↓")
    lines.append("                                                    xray_has_model → refine")
    lines.append("                                                              ↓")
    lines.append("                                                       xray_refined")
    lines.append("                                                     ↓    ↓    ↓")
    lines.append("                                              refine  autobuild  STOP")
    lines.append("```")
    lines.append("")

    xray_states = [s for s in states if s["type"] == "xray"]
    if xray_states:
        lines.append("| State | Valid Programs | Description |")
        lines.append("|-------|----------------|-------------|")
        for s in xray_states:
            progs = ", ".join(s["valid_programs"][:3])
            if len(s["valid_programs"]) > 3:
                progs += "..."
            lines.append(f"| {s['name']} | {progs} | {s['reason'][:50]} |")
        lines.append("")

    # Cryo-EM states
    lines.append("### Cryo-EM Workflow")
    lines.append("")
    lines.append("**Automated path** (maximum_automation=True):")
    lines.append("```")
    lines.append("cryoem_initial → mtriage → cryoem_analyzed → predict_and_build(full)")
    lines.append("                                                      ↓")
    lines.append("                                            cryoem_has_model → real_space_refine")
    lines.append("                                                      ↓")
    lines.append("                                               cryoem_refined → STOP")
    lines.append("```")
    lines.append("")
    lines.append("**Stepwise path** (maximum_automation=False):")
    lines.append("```")
    lines.append("cryoem_initial → mtriage → cryoem_analyzed → predict(stop_after_predict)")
    lines.append("                                                      ↓")
    lines.append("                                          process_predicted_model → dock_in_map")
    lines.append("                                                      ↓")
    lines.append("                                            cryoem_has_model → real_space_refine")
    lines.append("```")
    lines.append("")

    cryoem_states = [s for s in states if s["type"] == "cryoem"]
    if cryoem_states:
        lines.append("| State | Valid Programs | Description |")
        lines.append("|-------|----------------|-------------|")
        for s in cryoem_states:
            progs = ", ".join(s["valid_programs"][:3])
            if len(s["valid_programs"]) > 3:
                progs += "..."
            lines.append(f"| {s['name']} | {progs} | {s['reason'][:50]} |")
        lines.append("")

    # Program Templates
    lines.append("---")
    lines.append("## Program Templates")
    lines.append("")

    for prog_name, config in templates.items():
        lines.append(f"### {prog_name}")
        lines.append("")
        lines.append(f"**Description:** {config.get('description', 'N/A')}")
        lines.append("")

        # File slots
        if "file_slots" in config:
            lines.append("**Required Files:**")
            for slot, info in config["file_slots"].items():
                req = "required" if info.get("required") else "optional"
                exts = ", ".join(info.get("extensions", []))
                lines.append(f"- `{slot}`: {exts} ({req})")
            lines.append("")

        # Defaults
        if config.get("defaults"):
            lines.append("**Default Flags:**")
            for key, val in config["defaults"].items():
                lines.append(f"- `{key}={val}`")
            lines.append("")

        # Strategy flags
        if config.get("strategy_flags"):
            lines.append("**Strategy Options:**")
            for flag, options in config["strategy_flags"].items():
                if isinstance(options, dict):
                    if "format" in options:
                        lines.append(f"- `{flag}`: {options['format']}")
                    elif "true" in options:
                        lines.append(f"- `{flag}`: true → `{options['true']}`")
            lines.append("")

        # Hints
        if config.get("hints"):
            lines.append("**Usage Hints:**")
            for hint in config["hints"]:
                lines.append(f"- {hint}")
            lines.append("")

    # Stop Conditions
    lines.append("---")
    lines.append("## Stop Conditions")
    lines.append("")
    lines.append("The agent will recommend stopping when any of these conditions are met:")
    lines.append("")

    for cond in stop_conditions:
        lines.append(f"### {cond['name']}")
        lines.append(f"- **Description:** {cond['description']}")
        lines.append(f"- **Formula:** `{cond['formula']}`")
        lines.append("")

    # File Categorization
    lines.append("---")
    lines.append("## File Categorization")
    lines.append("")
    lines.append("Files are categorized by extension and naming patterns:")
    lines.append("")
    lines.append("| Category | Extensions/Pattern | Description |")
    lines.append("|----------|-------------------|-------------|")
    for rule in file_rules:
        ext = rule.get("extensions", rule.get("pattern", ""))
        lines.append(f"| {rule['name']} | {ext} | {rule['description']} |")
    lines.append("")

    # Key Thresholds
    lines.append("---")
    lines.append("## Key Thresholds")
    lines.append("")
    lines.append("Important decision thresholds used by the agent:")
    lines.append("")
    lines.append("| Threshold | Value | Description |")
    lines.append("|-----------|-------|-------------|")
    for t in thresholds:
        lines.append(f"| {t['name']} | {t['value']} | {t['description']} |")
    lines.append("")

    # Footer
    lines.append("---")
    lines.append("")
    lines.append("## Modifying Agent Behavior")
    lines.append("")
    lines.append("To change agent behavior, edit these files:")
    lines.append("")
    lines.append("| What to Change | File | Location |")
    lines.append("|----------------|------|----------|")
    lines.append("| **Thresholds & Defaults** | `agent/decision_config.json` | tiers, thresholds, program_rankings |")
    lines.append("| Program commands/flags | `agent/command_templates.json` | file_slots, defaults, strategy_flags |")
    lines.append("| Workflow state transitions | `agent/workflow_state.py` | _detect_xray_state(), _detect_cryoem_state() |")
    lines.append("| Stop conditions | `agent/metrics_analyzer.py` | analyze_metrics_trend() |")
    lines.append("| LLM prompts/guidance | `knowledge/prompts_hybrid.py` | SYSTEM_PROMPT, get_planning_prompt() |")
    lines.append("| File categorization | `agent/workflow_state.py` | _categorize_files() |")
    lines.append("")
    lines.append("### Configuration Priority")
    lines.append("")
    lines.append("1. `decision_config.json` - Centralized thresholds and tier definitions (preferred)")
    lines.append("2. `workflow_state.py` - State machine logic (uses config_loader)")
    lines.append("3. `metrics_analyzer.py` - Stop condition detection")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate AI agent logic documentation")
    parser.add_argument("--format", choices=["markdown", "text"], default="markdown",
                       help="Output format (default: markdown)")
    parser.add_argument("--output", "-o", type=str, default=None,
                       help="Output file (default: stdout)")
    args = parser.parse_args()

    # Extract all information
    config = load_decision_config()
    templates = load_command_templates()
    states = extract_workflow_states()
    stop_conditions = extract_stop_conditions()
    file_rules = extract_file_categorization()
    thresholds = extract_key_thresholds()

    # Generate documentation
    doc = generate_markdown(templates, states, stop_conditions, file_rules, thresholds, config)

    # Output
    if args.output:
        with open(args.output, "w") as f:
            f.write(doc)
        print(f"Documentation written to {args.output}", file=sys.stderr)
    else:
        print(doc)


if __name__ == "__main__":
    main()
