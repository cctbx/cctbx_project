"""
Workflow State Utilities for PHENIX AI Agent.

This module provides:
- File categorization by type and purpose
- History analysis (what programs have been run)
- Experiment type detection (X-ray vs Cryo-EM)
- Workflow state detection (delegating to WorkflowEngine)
- Prompt formatting for LLM

The actual workflow logic is defined in:
- knowledge/workflows.yaml (state machine)
- agent/workflow_engine.py (YAML interpreter)
"""

from __future__ import absolute_import, division, print_function
import os
import re


# =============================================================================
# FILE CATEGORIZATION
# =============================================================================

def _categorize_files(available_files):
    """
    Categorize files by type and purpose.

    Returns dict with keys:
        mtz, pdb, sequence, map, ligand_cif, ligand_pdb,
        phaser_output, refined, rsr_output, with_ligand, ligand_fit, predicted,
        full_map, half_map (for cryo-EM half-maps)
    """
    files = {
        "mtz": [],  # Also includes other X-ray data formats (.sca, .hkl, etc.)
        "refined_mtz": [],  # MTZ files output from refinement (have R-free flags)
        "pdb": [],
        "sequence": [],
        "map": [],  # All map files (for backward compatibility)
        "full_map": [],  # Full cryo-EM maps
        "half_map": [],  # Half maps (usually come in pairs)
        "ligand_cif": [],
        "ligand_pdb": [],
        "phaser_output": [],
        "refined": [],  # X-ray refinement output
        "rsr_output": [],  # real_space_refine output (cryo-EM)
        "with_ligand": [],
        "ligand_fit": [],
        "predicted": [],
        "processed_predicted": [],
        "autobuild_output": [],
        "docked": [],  # dock_in_map output
    }

    # X-ray data file extensions (PHENIX can read these)
    xray_data_extensions = ('.mtz', '.sca', '.hkl', '.sdf')

    def is_half_map(basename):
        """Detect if a map file is a half-map based on naming conventions."""
        name_lower = basename.lower()
        name_no_ext = os.path.splitext(name_lower)[0]

        if 'half' in name_lower:
            return True
        if re.search(r'[_\-]?[12]$', name_no_ext):
            return True
        if re.search(r'[_\-][ab]$', name_no_ext):
            return True

        return False

    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)

        # Primary type categorization
        if f_lower.endswith(xray_data_extensions):
            files["mtz"].append(f)
            # Also categorize refined MTZ files (output from phenix.refine)
            # These have R-free flags and should be preferred for subsequent refinement
            if 'refine' in basename or '_001' in basename or '_002' in basename:
                files["refined_mtz"].append(f)
        elif f_lower.endswith('.pdb'):
            files["pdb"].append(f)

            # Subcategorize PDBs by their origin/purpose
            if 'phaser' in basename or basename.startswith('phaser'):
                files["phaser_output"].append(f)

            if 'refine' in basename and 'real_space' not in basename and 'rsr' not in basename:
                files["refined"].append(f)

            # RSR output detection - real_space_refine outputs contain 'real_space_refined'
            # e.g., model_real_space_refined_000.pdb
            if 'real_space_refined' in basename or 'rsr_' in basename or '_rsr' in basename:
                files["rsr_output"].append(f)

            if 'with_ligand' in basename:
                files["with_ligand"].append(f)
            if 'ligand_fit' in basename or 'ligandfit' in basename:
                files["ligand_fit"].append(f)

            if 'predict' in basename or 'alphafold' in basename or 'colabfold' in basename:
                files["predicted"].append(f)

            if 'processed' in basename:
                files["processed_predicted"].append(f)

            is_autobuild = (
                ('autobuild' in basename or 'auto_build' in basename) or
                ('AutoBuild' in os.path.basename(f)) or
                ('overall_best' in basename) or
                ('build' in basename or 'built' in basename) or
                'buccaneer' in basename or
                'arp_warp' in basename or
                ('shelxe' in basename and 'trace' in basename)
            )
            if is_autobuild and 'predict' not in basename:
                files["autobuild_output"].append(f)

            if 'dock' in basename and 'map' in basename:
                files["docked"].append(f)

            if (basename.startswith('lig') and len(basename) < 20) or 'ligand' in basename:
                if not any(x in basename for x in ['ligand_fit', 'ligandfit', 'with_ligand']):
                    files["ligand_pdb"].append(f)

        elif f_lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            files["sequence"].append(f)
        elif f_lower.endswith(('.mrc', '.ccp4', '.map')):
            files["map"].append(f)
            if is_half_map(basename):
                files["half_map"].append(f)
            else:
                files["full_map"].append(f)
        elif f_lower.endswith('.cif'):
            if 'refine' in basename:
                files["pdb"].append(f)
                files["refined"].append(f)
            else:
                files["ligand_cif"].append(f)

    return files


def _detect_experiment_type(files, history_info=None):
    """
    Determine if this is X-ray crystallography or Cryo-EM.

    Logic:
    - If mtriage has been run → cryo-EM (definitive)
    - Has map (full or half) but no MTZ → cryo-EM
    - Has MTZ → X-ray (even if map also present)
    - Neither → unknown (default to X-ray)
    """
    if history_info and history_info.get("mtriage_done"):
        return "cryoem"

    if history_info and history_info.get("rsr_done"):
        return "cryoem"

    has_mtz = bool(files["mtz"])
    has_map = bool(files["map"]) or bool(files["full_map"]) or bool(files["half_map"])

    if has_map and not has_mtz:
        return "cryoem"
    else:
        return "xray"


# =============================================================================
# HISTORY ANALYSIS
# =============================================================================

def _analyze_history(history):
    """
    Extract information about what has been done from history.

    Returns dict with program completion flags and counts.
    """
    info = {
        "programs_run": set(),
        "xtriage_done": False,
        "mtriage_done": False,
        "phaser_done": False,
        "phaser_success": False,
        "predict_done": False,
        "predict_full_done": False,
        "process_predicted_done": False,
        "autobuild_done": False,
        "autosol_done": False,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,
        "rsr_count": 0,
        "ligandfit_done": False,
        "pdbtools_done": False,
        "dock_done": False,
        "validation_done": False,
        "last_program": None,
        "last_r_free": None,
        "last_map_cc": None,
        "last_tfz": None,
        "resolution": None,
        "anomalous_resolution": None,
        "has_anomalous": False,
        "has_twinning": False,
        "twin_law": None,
        "twin_fraction": None,
    }

    if not history:
        return info

    for entry in history:
        prog = ""
        cmd = ""

        if isinstance(entry, str):
            prog = entry.lower()
            cmd = entry.lower()
        elif isinstance(entry, dict):
            prog = entry.get("program", "").lower()
            cmd = entry.get("command", "").lower()

            # Extract metrics
            analysis = entry.get("analysis", {})
            if isinstance(analysis, dict):
                if analysis.get("r_free"):
                    info["last_r_free"] = analysis["r_free"]
                if analysis.get("map_cc"):
                    info["last_map_cc"] = analysis["map_cc"]
                if analysis.get("tfz"):
                    info["last_tfz"] = analysis["tfz"]
                if analysis.get("resolution"):
                    info["resolution"] = analysis["resolution"]
                if analysis.get("anomalous_resolution"):
                    info["anomalous_resolution"] = analysis["anomalous_resolution"]
                    info["has_anomalous"] = True
                elif analysis.get("has_anomalous"):
                    info["has_anomalous"] = analysis["has_anomalous"]
                # Store anomalous measurability for decision making
                if analysis.get("anomalous_measurability"):
                    info["anomalous_measurability"] = analysis["anomalous_measurability"]
                    # Strong anomalous signal if measurability > 0.10
                    if analysis["anomalous_measurability"] > 0.10:
                        info["has_anomalous"] = True
                        info["strong_anomalous"] = True
                # Twinning threshold (0.20) from workflows.yaml shared section
                if analysis.get("twin_law") and analysis.get("twin_fraction"):
                    twin_frac = analysis["twin_fraction"]
                    if twin_frac > 0.20 and not analysis.get("no_twinning_suspected"):
                        info["has_twinning"] = True
                        info["twin_law"] = analysis["twin_law"]
                        info["twin_fraction"] = twin_frac
        else:
            continue

        combined = prog + " " + cmd

        info["programs_run"].add(prog)
        info["last_program"] = prog

        # Program completion flags
        if "xtriage" in combined:
            info["xtriage_done"] = True
        if "mtriage" in combined:
            info["mtriage_done"] = True
        if "validation" in combined or "molprobity" in combined or "holton_geometry" in combined:
            info["validation_done"] = True
        if "phaser" in combined and "real_space" not in combined:
            info["phaser_done"] = True
            if info["last_tfz"] and info["last_tfz"] > 8:
                info["phaser_success"] = True
        if "predict_and_build" in combined:
            info["predict_done"] = True
            if "stop_after_predict=true" not in combined and "stop_after_predict=True" not in combined:
                info["predict_full_done"] = True
        if "process_predicted_model" in combined:
            info["process_predicted_done"] = True
        if "autobuild" in combined:
            info["autobuild_done"] = True
        if "autosol" in combined:
            info["autosol_done"] = True
        if "refine" in combined and "real_space" not in combined:
            info["refine_done"] = True
            info["refine_count"] += 1
        if "real_space_refine" in combined:
            info["rsr_done"] = True
            info["rsr_count"] += 1
        if "ligandfit" in combined:
            info["ligandfit_done"] = True
        if "pdbtools" in combined:
            info["pdbtools_done"] = True
        if "dock_in_map" in combined:
            info["dock_done"] = True

    return info


# =============================================================================
# WORKFLOW STATE DETECTION
# =============================================================================

# Import workflow engine for YAML-driven detection
from libtbx.langchain.agent.workflow_engine import WorkflowEngine


def detect_workflow_state(history, available_files, analysis=None, maximum_automation=True,
                         use_yaml_engine=True):
    """
    Determine current workflow state based on history and files.

    This function delegates to the YAML-driven WorkflowEngine for state detection.

    Args:
        history: List of cycle records from client
        available_files: List of available file paths
        analysis: Current log analysis dict (optional)
        maximum_automation: If True, use fully automated cryo-EM path
        use_yaml_engine: If True, use YAML-driven WorkflowEngine (default: True)

    Returns:
        dict: {
            state: str,              # State name
            experiment_type: str,    # "xray" or "cryoem"
            valid_programs: list,    # Programs allowed in this state
            reason: str,             # Human-readable explanation
            conditions: dict,        # Conditional program availability
            automation_path: str,    # "stepwise" or "automated" (cryo-EM only)
            categorized_files: dict, # Pre-categorized files (full_map, half_map, etc.)
        }
    """
    # Categorize files
    files = _categorize_files(available_files)

    # Analyze history
    history_info = _analyze_history(history)

    # Determine experiment type
    experiment_type = _detect_experiment_type(files, history_info)

    # Use YAML-driven workflow engine
    if use_yaml_engine:
        try:
            engine = WorkflowEngine()
            state = engine.get_workflow_state(experiment_type, files, history_info, analysis)
            state["categorized_files"] = files
            # Set automation_path for both experiment types
            state["automation_path"] = "automated" if maximum_automation else "stepwise"
            return state
        except Exception as e:
            import sys
            print("Warning: YAML workflow engine failed: %s" % e, file=sys.stderr)

    # Fallback: return minimal state (should not happen if YAML is properly configured)
    return {
        "state": "unknown",
        "experiment_type": experiment_type,
        "valid_programs": ["STOP"],
        "reason": "Workflow engine unavailable",
        "conditions": {},
        "automation_path": "automated" if maximum_automation else "stepwise",
        "categorized_files": files,
    }


# =============================================================================
# VALIDATION
# =============================================================================

def validate_program_choice(chosen_program, workflow_state):
    """
    Validate that a program choice is allowed in the current state.

    Args:
        chosen_program: Program the LLM chose
        workflow_state: Dict from detect_workflow_state()

    Returns:
        tuple: (is_valid: bool, error_message: str or None)
    """
    if chosen_program is None:
        return True, None

    if chosen_program == "STOP":
        return True, None

    valid = workflow_state["valid_programs"]

    if chosen_program in valid:
        return True, None

    all_known_programs = [
        "phenix.xtriage", "phenix.mtriage", "phenix.predict_and_build",
        "phenix.phaser", "phenix.refine", "phenix.real_space_refine",
        "phenix.ligandfit", "phenix.pdbtools", "phenix.dock_in_map",
        "phenix.autobuild", "phenix.autosol", "phenix.process_predicted_model",
        "phenix.molprobity", "phenix.resolve_cryo_em", "phenix.map_sharpening",
    ]

    if chosen_program in all_known_programs:
        error = (
            "Program '%s' is not valid in state '%s'. "
            "Valid programs: %s. Reason: %s"
        ) % (
            chosen_program,
            workflow_state["state"],
            ", ".join(valid),
            workflow_state["reason"]
        )
    else:
        error = "Unknown program '%s'. Valid programs: %s" % (chosen_program, ", ".join(valid))

    return False, error


# =============================================================================
# PROMPT FORMATTING
# =============================================================================

def format_workflow_for_prompt(workflow_state):
    """
    Format workflow state for inclusion in LLM prompt.

    Args:
        workflow_state: Output from detect_workflow_state()

    Returns:
        str: Formatted text for prompt
    """
    lines = []

    lines.append("### WORKFLOW STATE: %s" % workflow_state["state"])
    lines.append("Experiment type: %s" % workflow_state["experiment_type"])

    if workflow_state.get("automation_path"):
        lines.append("Automation path: %s" % workflow_state["automation_path"])

    lines.append("")
    lines.append(workflow_state["reason"])
    lines.append("")
    lines.append("**VALID PROGRAMS FOR THIS STATE:**")
    lines.append(", ".join(workflow_state["valid_programs"]))
    lines.append("")
    lines.append("⚠️ You MUST choose a program from the list above, or set \"stop\": true.")
    lines.append("Choosing an invalid program will cause a validation error.")

    if workflow_state.get("conditions"):
        lines.append("")
        lines.append("Conditional availability:")
        for prog, condition in workflow_state["conditions"].items():
            lines.append("  - %s: requires %s" % (prog, condition))

    if (workflow_state["experiment_type"] == "cryoem" and
        workflow_state.get("automation_path") == "stepwise"):

        if workflow_state["state"] == "cryoem_analyzed":
            lines.append("")
            lines.append("NOTE (Stepwise mode): Use predict_and_build with strategy: {\"stop_after_predict\": true}")

    return "\n".join(lines)
