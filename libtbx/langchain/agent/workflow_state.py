"""
Workflow State Machine for PHENIX AI Agent.

Defines valid states and transitions for X-ray and Cryo-EM workflows.
Provides strict enforcement of program choices based on current state.

Key Features:
- X-ray workflow: xtriage → predict/phaser → refine → [ligandfit → pdbtools] → refine
- Cryo-EM workflow: mtriage → predict/dock → real_space_refine → [ligandfit]
- Cryo-EM has two paths controlled by maximum_automation flag:
  - Stepwise (False): More control, intermediate checkpoints
  - Automated (True): Fewer steps, fully automated
- Tiered decision architecture:
  - Tier 1: Hard constraints (cannot override)
  - Tier 2: Strong defaults (LLM can override with warning)
  - Tier 3: Soft guidance (suggestions only)
"""

from __future__ import absolute_import, division, print_function
import os

# Import config_loader (handle both package and direct import)
try:
    from . import config_loader
except ImportError:
    import config_loader


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
    # .mtz - CCP4 MTZ format (most common)
    # .sca - Scalepack format (HKL2000/HKL3000)
    # .hkl - Various formats (SHELX, XDS, etc.)
    # .sdf - structure factor file
    xray_data_extensions = ('.mtz', '.sca', '.hkl', '.sdf')
    
    # Half-map naming patterns (case-insensitive)
    # Common patterns: half_map_1, half_map_2, half1, half2, 
    #                  map_1, map_2, map1, map2, _a, _b, _half1, _half2
    def is_half_map(basename):
        """Detect if a map file is a half-map based on naming conventions."""
        name_lower = basename.lower()
        # Remove extension for pattern matching
        name_no_ext = os.path.splitext(name_lower)[0]
        
        # Explicit half-map patterns
        if 'half' in name_lower:
            return True
        
        # Numbered patterns at end: _1, _2, -1, -2, 1, 2
        # But be careful not to match things like "map_final"
        import re
        # Match patterns like: map_1, map_2, map-1, map-2, map1, map2
        if re.search(r'[_\-]?[12]$', name_no_ext):
            return True
        # Match patterns like: _a, _b, -a, -b (at end of name)
        if re.search(r'[_\-][ab]$', name_no_ext):
            return True
        
        return False
    
    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)
        
        # Primary type categorization
        if f_lower.endswith(xray_data_extensions):
            files["mtz"].append(f)  # All X-ray data goes in "mtz" category
        elif f_lower.endswith('.pdb'):
            files["pdb"].append(f)
            
            # Subcategorize PDBs by their origin/purpose
            if 'phaser' in basename or basename.startswith('phaser'):
                files["phaser_output"].append(f)
            
            # refined: X-ray refinement output (exclude real_space_refine)
            if 'refine' in basename and 'real_space' not in basename and 'rsr' not in basename:
                files["refined"].append(f)
            
            # rsr_output: real_space_refine output (cryo-EM)
            if 'real_space_refined' in basename or '_rsr_' in basename:
                files["rsr_output"].append(f)
            
            if 'with_ligand' in basename:
                files["with_ligand"].append(f)
            if 'ligand_fit' in basename or 'ligandfit' in basename:
                files["ligand_fit"].append(f)
            
            # predicted: AlphaFold/ColabFold predictions (check 'predict' not just 'predicted')
            if 'predict' in basename or 'alphafold' in basename or 'colabfold' in basename:
                files["predicted"].append(f)
            
            # processed_predicted: after process_predicted_model
            if 'processed' in basename:
                files["processed_predicted"].append(f)
            
            # autobuild: model building programs (NOT predict_and_build outputs)
            # Include: autobuild, buccaneer, arp_warp, shelxe builds
            # Also include phenix.autobuild outputs like overall_best.pdb
            # phenix.autobuild creates AutoBuild_run_X_ directories
            is_autobuild = (
                ('autobuild' in basename or 'auto_build' in basename) or
                ('AutoBuild' in os.path.basename(f)) or  # Case-sensitive directory name
                ('overall_best' in basename) or  # phenix.autobuild output
                ('build' in basename or 'built' in basename) or
                'buccaneer' in basename or
                'arp_warp' in basename or
                ('shelxe' in basename and 'trace' in basename)
            )
            # Exclude predict_and_build outputs
            if is_autobuild and 'predict' not in basename:
                files["autobuild_output"].append(f)
            
            # dock_in_map output
            if 'dock' in basename and 'map' in basename:
                files["docked"].append(f)
            
            # Small PDB files with 'lig' in name are likely ligand coordinates
            if (basename.startswith('lig') and len(basename) < 20) or 'ligand' in basename:
                if not any(x in basename for x in ['ligand_fit', 'ligandfit', 'with_ligand']):
                    files["ligand_pdb"].append(f)
                    
        elif f_lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            files["sequence"].append(f)
        elif f_lower.endswith(('.mrc', '.ccp4', '.map')):
            # Categorize maps as full_map or half_map
            files["map"].append(f)  # All maps go here for backward compatibility
            if is_half_map(basename):
                files["half_map"].append(f)
            else:
                files["full_map"].append(f)
        elif f_lower.endswith('.cif'):
            # Distinguish between model CIF (from refinement) and ligand CIF
            # Model CIFs typically have 'refine' in the name (output from phenix.refine)
            # Ligand CIFs are restraint files, typically small with ligand-related names
            if 'refine' in basename:
                # This is a model CIF (mmCIF format from refinement), treat like a PDB
                files["pdb"].append(f)
                files["refined"].append(f)
            else:
                # Assume it's a ligand restraint CIF
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
    # If mtriage was run, it's definitely cryo-EM
    if history_info and history_info.get("mtriage_done"):
        return "cryoem"
    
    # If real_space_refine was run, it's cryo-EM
    if history_info and history_info.get("rsr_done"):
        return "cryoem"
    
    has_mtz = bool(files["mtz"])
    has_map = bool(files["map"]) or bool(files["full_map"]) or bool(files["half_map"])
    
    if has_map and not has_mtz:
        return "cryoem"
    else:
        return "xray"  # Default to X-ray


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
        "phaser_success": False,  # TFZ > 8
        "predict_done": False,
        "predict_full_done": False,  # predict_and_build without stop_after_predict
        "process_predicted_done": False,
        "autobuild_done": False,
        "autosol_done": False,
        "refine_done": False,
        "refine_count": 0,
        "rsr_done": False,  # real_space_refine
        "rsr_count": 0,
        "ligandfit_done": False,
        "pdbtools_done": False,
        "dock_done": False,
        "validation_done": False,  # Whether validation has been run
        "last_program": None,
        "last_r_free": None,
        "last_map_cc": None,
        "last_tfz": None,
        "resolution": None,  # From xtriage or mtriage
        "anomalous_resolution": None,  # Resolution of useful anomalous signal
        "has_anomalous": False,  # Whether anomalous data exists
        "has_twinning": False,  # Whether twinning is detected
        "twin_law": None,  # Twin law operator if twinning detected
        "twin_fraction": None,  # Britton alpha
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
                # Anomalous data from xtriage
                if analysis.get("anomalous_resolution"):
                    info["anomalous_resolution"] = analysis["anomalous_resolution"]
                    info["has_anomalous"] = True
                elif analysis.get("has_anomalous"):
                    info["has_anomalous"] = analysis["has_anomalous"]
                # Twinning from xtriage (only if twin_fraction > 0.20 and no "No twinning" message)
                if analysis.get("twin_law") and analysis.get("twin_fraction"):
                    twin_frac = analysis["twin_fraction"]
                    if twin_frac > 0.20 and not analysis.get("no_twinning_suspected"):
                        info["has_twinning"] = True
                        info["twin_law"] = analysis["twin_law"]
                        info["twin_fraction"] = twin_frac
        else:
            continue
        
        # Normalize: check both prog and cmd for program names
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
            # Check for successful MR (TFZ > 8)
            if info["last_tfz"] and info["last_tfz"] > 8:
                info["phaser_success"] = True
        if "predict_and_build" in combined:
            info["predict_done"] = True
            # Check if it was a full run or just prediction
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
# X-RAY STATE DETECTION
# =============================================================================

def _detect_xray_state(files, history_info, analysis):
    """
    Detect current state in X-ray crystallography workflow.
    
    States:
    - xray_initial: Need xtriage
    - xray_analyzed: Need model (predict or phaser)
    - xray_has_model: Need refinement
    - xray_refined: Can do ligand, more refine, or stop
    - xray_has_ligand: Need pdbtools to combine
    - xray_combined: Final refinement
    - complete: Done
    """
    # Check for various conditions
    has_ligand = bool(files["ligand_cif"] or files["ligand_pdb"])
    has_model = bool(files["pdb"])
    has_sequence = bool(files["sequence"])
    has_search_model = bool([f for f in files["pdb"] if not any(
        x in os.path.basename(f).lower() 
        for x in ['refine', 'ligand', 'phaser', 'autobuild', 'predicted']
    )])
    
    # Check for predicted models (from AlphaFold/predict_and_build)
    has_predicted_model = bool(files["predicted"]) or history_info["predict_done"]
    has_processed_predicted = bool(files["processed_predicted"]) or history_info["process_predicted_done"]
    
    # For X-ray, a model is "placed" only after phaser or autobuild
    # predict_and_build alone does NOT place the model - it just generates an AlphaFold prediction
    has_placed_model = bool(
        files["phaser_output"] or 
        files["refined"] or 
        files["autobuild_output"] or
        history_info["phaser_done"] or 
        history_info["autobuild_done"]
    )
    has_refined = bool(files["refined"]) or history_info["refine_count"] > 0
    has_ligand_fit = bool(files["ligand_fit"]) or history_info["ligandfit_done"]
    has_combined = bool(files["with_ligand"]) or history_info["pdbtools_done"]
    
    # Get current R-free
    r_free = None
    if analysis and analysis.get("r_free"):
        r_free = analysis["r_free"]
    elif history_info["last_r_free"]:
        r_free = history_info["last_r_free"]
    
    # === STATE DETECTION (order matters!) ===
    
    # 1. Initial: Need xtriage
    if not history_info["xtriage_done"]:
        return {
            "state": "xray_initial",
            "experiment_type": "xray",
            "valid_programs": ["phenix.xtriage"],
            "reason": "Need to analyze data quality with xtriage first",
            "conditions": {},
        }
    
    # 2. Have predicted model but need to process and run MR
    if has_predicted_model and not has_placed_model:
        if not has_processed_predicted:
            # Need to process the predicted model first
            return {
                "state": "xray_has_prediction",
                "experiment_type": "xray",
                "valid_programs": ["phenix.process_predicted_model"],
                "reason": "Have AlphaFold prediction, need to process it for molecular replacement",
                "conditions": {},
            }
        else:
            # Have processed model, need phaser
            return {
                "state": "xray_model_processed",
                "experiment_type": "xray",
                "valid_programs": ["phenix.phaser"],
                "reason": "Predicted model processed, need molecular replacement with phaser",
                "conditions": {},
            }
    
    # 3. Analyzed but no model yet: Need to get a model
    if not has_placed_model:
        valid = []
        conditions = {}
        reason_parts = []
        
        # Option 1: AlphaFold prediction (if we have sequence)
        if has_sequence:
            valid.append("phenix.predict_and_build")
        
        # Option 2: Molecular replacement (if we have a search model)
        if has_search_model or has_model:
            valid.append("phenix.phaser")
        
        # Option 3: Experimental phasing with autosol
        # Requires: sequence file
        # Note: We always offer autosol when sequence is available because:
        # - The user may have SAD/MAD data that xtriage didn't explicitly flag
        # - xtriage may not have been run yet or may not have detected anomalous signal
        # - The LLM and user can decide if experimental phasing is appropriate
        anomalous_res = history_info.get("anomalous_resolution")
        has_useful_anomalous = (
            history_info.get("has_anomalous") and 
            anomalous_res is not None and 
            anomalous_res <= 4.0
        )
        
        if has_sequence and not history_info.get("autosol_done"):
            valid.append("phenix.autosol")
            if has_useful_anomalous:
                if anomalous_res <= 2.5:
                    conditions["phenix.autosol"] = "strong_anomalous_signal_%.1fA" % anomalous_res
                    reason_parts.append("strong anomalous signal to %.1fÅ - autosol recommended" % anomalous_res)
                else:
                    conditions["phenix.autosol"] = "anomalous_signal_%.1fA" % anomalous_res
                    reason_parts.append("anomalous signal to %.1fÅ - autosol possible" % anomalous_res)
            else:
                conditions["phenix.autosol"] = "requires_anomalous_data"
                reason_parts.append("autosol available if data has anomalous signal (SAD/MAD)")
        
        if not valid:
            # Can't proceed - no sequence or model
            return {
                "state": "xray_analyzed",
                "experiment_type": "xray",
                "valid_programs": ["STOP"],
                "reason": "STUCK: No sequence file for AlphaFold and no search model for MR. Upload a sequence (.fa) or search model (.pdb)",
                "conditions": {},
            }
        
        reason = "Data analyzed, need to obtain placed model"
        if reason_parts:
            reason += " (" + ", ".join(reason_parts) + ")"
        
        return {
            "state": "xray_analyzed",
            "experiment_type": "xray",
            "valid_programs": valid,
            "reason": reason,
            "conditions": conditions,
        }
    
    # 3b. After autosol: have phases, need autobuild
    # Autosol produces phases but may not produce a complete model
    if history_info.get("autosol_done") and not history_info.get("autobuild_done") and not has_refined:
        if has_sequence:
            return {
                "state": "xray_has_phases",
                "experiment_type": "xray",
                "valid_programs": ["phenix.autobuild"],
                "reason": "Experimental phasing complete, need to build model with autobuild",
                "conditions": {},
            }
    
    # 4. Has placed model but not refined yet
    if not has_refined:
        return {
            "state": "xray_has_model",
            "experiment_type": "xray",
            "valid_programs": ["phenix.refine"],
            "reason": "Have placed model, need initial refinement",
            "conditions": {},
        }
    
    # 4. Ligand workflow: fitted but not combined
    if has_ligand_fit and not has_combined:
        return {
            "state": "xray_has_ligand",
            "experiment_type": "xray",
            "valid_programs": ["phenix.pdbtools"],
            "reason": "Ligand fitted, need to combine with protein model using pdbtools",
            "conditions": {},
        }
    
    # 5. Ligand workflow: combined, need final refine
    if has_combined:
        valid = ["phenix.refine"]
        reason = "Model + ligand combined"
        
        # === VALIDATION GATE FOR STOP (same logic as xray_refined) ===
        # Get resolution-dependent thresholds
        resolution = history_info.get("resolution")
        if resolution:
            thresholds = config_loader.get_resolution_thresholds(resolution)
        else:
            thresholds = config_loader.get_resolution_thresholds(2.0)
        good_model_threshold = thresholds["good_model_rfree"]
        success_threshold = good_model_threshold - 0.02
        
        # Check if validation should be suggested (same criteria as xray_refined)
        suggest_validation = False
        if not history_info.get("validation_done"):
            if r_free and r_free < good_model_threshold:
                suggest_validation = True
            elif history_info["refine_count"] >= 3:
                suggest_validation = True
        
        # Validation required if suggest_validation OR r_free below success threshold
        validation_required_before_stop = False
        if not history_info.get("validation_done"):
            if suggest_validation:
                validation_required_before_stop = True
            elif r_free is not None and r_free < success_threshold:
                validation_required_before_stop = True
        
        if validation_required_before_stop:
            # Add molprobity as required
            valid.insert(0, "phenix.molprobity")
            reason += " - VALIDATION REQUIRED BEFORE STOPPING"
        elif r_free and r_free < 0.30:
            reason += ", R-free (%.3f) looks good" % r_free
        
        # Add STOP - prioritize it if validation is done and model is good
        if not validation_required_before_stop:
            if history_info.get("validation_done") and r_free and r_free <= good_model_threshold:
                # Validation done and model is good - STOP should be first choice
                valid.insert(0, "STOP")
                reason += " - validation complete, consider stopping"
            else:
                valid.append("STOP")
                if r_free and r_free < 0.30:
                    reason += " - consider stopping"
        
        return {
            "state": "xray_combined",
            "experiment_type": "xray",
            "valid_programs": valid,
            "reason": reason,
            "conditions": {},
        }
    
    # 6. Refined: Can do ligand fitting, more refinement, autobuild, or stop
    # === USE CONFIG FOR THRESHOLDS ===
    
    # Get resolution for resolution-dependent logic
    resolution = None
    if analysis and analysis.get("resolution"):
        resolution = analysis["resolution"]
    elif history_info.get("resolution"):
        resolution = history_info["resolution"]
    
    # Get thresholds from config (resolution-dependent)
    if resolution:
        thresholds = config_loader.get_resolution_thresholds(resolution)
    else:
        thresholds = config_loader.get_resolution_thresholds(2.0)  # Default to medium
    
    autobuild_threshold = thresholds["autobuild_rfree"]
    good_model_threshold = thresholds["good_model_rfree"]
    ligandfit_threshold = thresholds["ligandfit_rfree"]
    
    # Determine if model quality is good enough for ligand fitting
    model_is_good = r_free is not None and r_free < ligandfit_threshold
    
    # === FIX 3: AUTOBUILD REQUIRES PHASES ===
    # Autobuild needs phased MTZ (from phaser or refine output)
    has_phases = history_info["phaser_done"] or history_info["refine_count"] > 0
    
    # Check if autobuild has been run - either from history OR from presence of autobuild output files
    autobuild_already_done = (
        history_info.get("autobuild_done", False) or 
        bool(files.get("autobuild_output"))
    )
    
    # === VALIDATION SUGGESTION ===
    suggest_validation = False
    validation_reason = None
    
    if not history_info.get("validation_done"):
        if r_free and r_free < good_model_threshold:
            suggest_validation = True
            validation_reason = "r_free_target_reached"
        elif history_info["refine_count"] >= 3:
            suggest_validation = True
            validation_reason = "multiple_refine_cycles"
    
    # === BUILD CONTEXT FOR CONFIG EVALUATION ===
    # Helper: autobuild should be available if R-free is unknown OR above threshold
    r_free_unknown_or_high = (r_free is None) or (r_free > autobuild_threshold)
    r_free_unknown_or_above_target = (r_free is None) or (r_free > good_model_threshold)
    
    context = {
        "r_free": r_free,
        "resolution": resolution,
        "refine_count": history_info["refine_count"],
        "has_ligand": has_ligand,
        "has_ligand_fit": has_ligand_fit,
        "model_is_good": model_is_good,
        "has_sequence": has_sequence,
        "has_phases": has_phases,
        "autobuild_done": autobuild_already_done,
        "autobuild_threshold": autobuild_threshold,
        "good_model_threshold": good_model_threshold,
        "suggest_validation": suggest_validation,
        "has_twinning": history_info.get("has_twinning", False),
        "twin_law": history_info.get("twin_law"),
        "twin_fraction": history_info.get("twin_fraction"),
        "validation_done": history_info.get("validation_done", False),
        # Helper values for conditions
        "r_free_unknown_or_high": r_free_unknown_or_high,
        "r_free_unknown_or_above_target": r_free_unknown_or_above_target,
    }
    
    # === GET RANKED PROGRAMS FROM CONFIG ===
    rankings = config_loader.get_program_rankings("xray_refined", context)
    
    # Build valid programs list from rankings (only those with conditions met)
    valid = []
    ranked_programs = []
    for item in rankings:
        if item["condition_met"]:
            valid.append(item["program"])
            ranked_programs.append({
                "program": item["program"],
                "rank": item["rank"],
                "reason": item["reason"],
            })
    
    # Ensure we always have at least refine
    if "phenix.refine" not in valid:
        valid.append("phenix.refine")
    
    # === VALIDATION GATE FOR STOP ===
    # Require validation before allowing STOP if:
    # 1. R-free is below success threshold (model is good), OR
    # 2. suggest_validation is True (3+ refine cycles or R-free below target)
    # In either case, validation hasn't been done yet
    validation_required_before_stop = False
    
    # Check if validation is suggested but not done
    if suggest_validation and not history_info.get("validation_done"):
        validation_required_before_stop = True
        validation_gate_reason = validation_reason or "validation_suggested"
    
    # Also check success threshold
    if r_free is not None and not validation_required_before_stop:
        success_threshold = good_model_threshold - 0.02
        if r_free < success_threshold and not history_info.get("validation_done"):
            validation_required_before_stop = True
            validation_gate_reason = "r_free_below_success_threshold"
    
    # Apply validation gate
    if validation_required_before_stop:
        # Ensure molprobity is in valid programs
        if "phenix.molprobity" not in valid:
            valid.insert(0, "phenix.molprobity")  # High priority
            ranked_programs.insert(0, {
                "program": "phenix.molprobity",
                "rank": 1,
                "reason": "Validation required before stopping (%s)" % validation_gate_reason,
            })
        # Remove STOP from valid programs
        if "STOP" in valid:
            valid.remove("STOP")
    
    # Add STOP only if validation not required
    if "STOP" not in valid and not validation_required_before_stop:
        valid.append("STOP")
    
    # Build reason string
    reason = "Model refined"
    if r_free:
        reason += " (R-free: %.3f)" % r_free
    if resolution:
        reason += " at %.1fÅ resolution" % resolution
    if validation_required_before_stop:
        reason += " - VALIDATION REQUIRED BEFORE STOPPING"
    
    # === GET TIER 2 DEFAULTS (RECOMMENDED STRATEGY) ===
    recommended_strategy = config_loader.get_tier2_defaults("phenix.refine", context)
    
    # === GET TIER 3 SUGGESTIONS ===
    soft_suggestions = config_loader.get_tier3_suggestions(context)
    
    # === BUILD HARD CONSTRAINTS LIST ===
    hard_constraints = []
    if history_info["refine_count"] > 0:
        hard_constraints.append("Must use existing R-free flags from previous refinement")
    
    # Add twinning notice if detected
    if history_info.get("has_twinning") and history_info.get("twin_law"):
        reason += ", TWINNING DETECTED (%.0f%% twin fraction)" % (history_info["twin_fraction"] * 100)
    
    return {
        "state": "xray_refined",
        "experiment_type": "xray",
        "valid_programs": valid,
        "reason": reason,
        "conditions": {},  # Legacy field
        
        # === NEW RECOMMENDATION STRUCTURE ===
        "ranked_programs": ranked_programs,
        "recommended_strategy": recommended_strategy,
        "hard_constraints": hard_constraints,
        "soft_suggestions": soft_suggestions,
        
        # === CONTEXT FOR LLM ===
        "resolution": resolution,
        "thresholds": {
            "autobuild_rfree": autobuild_threshold,
            "good_model_rfree": good_model_threshold,
            "ligandfit_rfree": ligandfit_threshold,
        },
        
        # === LEGACY FIELDS (for backward compatibility) ===
        "refine_hints": list(recommended_strategy.keys()),  # Simple list for old code
        "has_twinning": history_info.get("has_twinning", False),
        "twin_law": history_info.get("twin_law"),
        "high_res_suggestions": [s for s in soft_suggestions if "anisotropic" in s.lower() or "hydrogen" in s.lower()],
        "suggest_validation": suggest_validation,
        "validation_reason": validation_reason,
    }


# =============================================================================
# CRYO-EM STATE DETECTION
# =============================================================================

def _detect_cryoem_state(files, history_info, analysis, maximum_automation):
    """
    Detect current state in Cryo-EM workflow.
    
    Two paths controlled by maximum_automation:
    
    Stepwise (False):
      mtriage → predict(stop_after_predict) → process_predicted_model → 
      phaser → autobuild → real_space_refine
      
    Automated (True):
      mtriage → predict(full) → real_space_refine
      
    Note: real_space_refine requires a full map, not half-maps.
    """
    has_model = bool(files["pdb"])
    has_sequence = bool(files["sequence"])
    has_ligand = bool(files["ligand_cif"] or files["ligand_pdb"])
    has_predicted = bool(files["predicted"]) or history_info["predict_done"]
    has_processed = bool(files["processed_predicted"]) or history_info["process_predicted_done"]
    
    # Check if we have a full map (required for real_space_refine)
    # Half-maps alone are not sufficient for real_space_refine
    has_full_map = bool(files["full_map"])
    has_half_maps = len(files["half_map"]) >= 2
    
    # For automated path, model is placed after full predict_and_build or dock
    # For stepwise path, model is placed only after autobuild completes
    if maximum_automation:
        has_placed = bool(
            files["autobuild_output"] or
            files["docked"] or  # dock_in_map output
            history_info["autobuild_done"] or
            history_info["dock_done"] or
            history_info["predict_full_done"]  # Full predict_and_build
        )
    else:
        # Stepwise: only after autobuild (phaser alone doesn't count)
        has_placed = bool(
            files["autobuild_output"] or
            files["docked"] or  # dock_in_map output
            history_info["autobuild_done"] or
            history_info["dock_done"]
        )
    
    # Check if model has been refined (either from history or from file presence)
    has_refined = history_info["rsr_count"] > 0 or bool(files["rsr_output"])
    has_ligand_fit = bool(files["ligand_fit"]) or history_info["ligandfit_done"]
    
    # Get current map-model CC
    map_cc = None
    if analysis and analysis.get("map_cc"):
        map_cc = analysis["map_cc"]
    elif history_info["last_map_cc"]:
        map_cc = history_info["last_map_cc"]
    
    # Get resolution from mtriage analysis or history
    resolution = None
    if analysis and analysis.get("resolution"):
        resolution = analysis["resolution"]
    elif history_info.get("resolution"):
        resolution = history_info["resolution"]
    
    automation_path = "automated" if maximum_automation else "stepwise"
    
    # === STATE DETECTION ===
    
    # 1. Initial: Need mtriage
    if not history_info["mtriage_done"]:
        return {
            "state": "cryoem_initial",
            "experiment_type": "cryoem",
            "valid_programs": ["phenix.mtriage"],
            "reason": "Need to analyze map quality with mtriage first",
            "conditions": {},
            "automation_path": automation_path,
            "resolution": resolution,
        }
    
    # 2. Analyzed: Need to get model
    if not has_placed:
        valid = []
        
        # Path diverges here based on automation setting
        if maximum_automation:
            # Automated: just predict_and_build (full) or dock if model exists
            if has_sequence:
                valid.append("phenix.predict_and_build")
            if has_model:
                valid.append("phenix.dock_in_map")
        else:
            # Stepwise path
            if not has_predicted:
                # Need to get AlphaFold model first
                if has_sequence:
                    valid.append("phenix.predict_and_build")  # with stop_after_predict=True
                if has_model:
                    valid.append("phenix.dock_in_map")
            elif not has_processed:
                # Have predicted model, need to process it
                valid.append("phenix.process_predicted_model")
            elif not history_info["phaser_done"]:
                # Have processed model, need phaser
                valid.append("phenix.phaser")
            elif not history_info["autobuild_done"]:
                # Have phases, need autobuild
                valid.append("phenix.autobuild")
        
        # If only half-maps available, also offer map optimization programs
        # These create a full sharpened map that can be used for real_space_refine later
        if has_half_maps and not has_full_map:
            valid.append("phenix.resolve_cryo_em")
            valid.append("phenix.map_sharpening")
        
        if not valid:
            return {
                "state": "cryoem_analyzed",
                "experiment_type": "cryoem",
                "valid_programs": ["STOP"],
                "reason": "STUCK: No sequence for AlphaFold and no model to dock. Upload a sequence (.fa) or model (.pdb)",
                "conditions": {},
                "automation_path": automation_path,
                "resolution": resolution,
            }
        
        reason = "Map analyzed"
        if not maximum_automation:
            if not has_predicted:
                reason += ", need to obtain AlphaFold model"
            elif not has_processed:
                reason += ", need to process predicted model for MR"
            elif not history_info["phaser_done"]:
                reason += ", need molecular replacement with phaser"
            elif not history_info["autobuild_done"]:
                reason += ", need to build model with autobuild"
        else:
            reason += ", need to obtain/dock model"
        
        if has_half_maps and not has_full_map:
            reason += " (map optimization available: resolve_cryo_em or map_sharpening)"
        
        # Special state names for stepwise path
        if not maximum_automation:
            if has_predicted and not has_processed:
                state_name = "cryoem_has_prediction"
            elif has_processed and not history_info["phaser_done"]:
                state_name = "cryoem_model_processed"
            elif history_info["phaser_done"] and not history_info["autobuild_done"]:
                state_name = "cryoem_has_phases"
            else:
                state_name = "cryoem_analyzed"
        else:
            state_name = "cryoem_analyzed"
        
        return {
            "state": state_name,
            "experiment_type": "cryoem",
            "valid_programs": valid,
            "reason": reason,
            "conditions": {},
            "automation_path": automation_path,
            "resolution": resolution,
        }
    
    # 3. Has model, need real-space refinement
    if not has_refined:
        if has_full_map:
            return {
                "state": "cryoem_has_model",
                "experiment_type": "cryoem",
                "valid_programs": ["phenix.real_space_refine"],
                "reason": "Have model in map, need real-space refinement",
                "conditions": {},
                "automation_path": automation_path,
                "resolution": resolution,
            }
        elif has_half_maps:
            # Only half-maps available - need to create a full map first
            return {
                "state": "cryoem_needs_full_map",
                "experiment_type": "cryoem",
                "valid_programs": ["phenix.resolve_cryo_em", "phenix.map_sharpening"],
                "reason": "Have model but only half-maps available. Need to create optimized full map before real_space_refine. Try resolve_cryo_em first.",
                "conditions": {},
                "automation_path": automation_path,
                "resolution": resolution,
            }
        else:
            return {
                "state": "cryoem_no_map",
                "experiment_type": "cryoem",
                "valid_programs": ["STOP"],
                "reason": "Have model but no maps available.",
                "conditions": {},
                "automation_path": automation_path,
                "resolution": resolution,
            }
    
    # 4. Refined: Can do ligand, more refine, or stop
    valid = []
    if has_full_map:
        valid.append("phenix.real_space_refine")
    conditions = {}
    reason = "Model refined against map"
    if map_cc:
        reason += " (CC: %.3f)" % map_cc
    
    if has_ligand and not has_ligand_fit:
        valid.insert(0, "phenix.ligandfit")
        conditions["phenix.ligandfit"] = "has_ligand_file"
        reason += ", ligand available for fitting"
    
    # === VALIDATION GATE FOR STOP (cryo-EM uses map_cc) ===
    validation_required_before_stop = False
    if map_cc and map_cc > 0.75 and not history_info.get("validation_done"):
        validation_required_before_stop = True
        # Add validation as required
        if "phenix.molprobity" not in valid:
            valid.insert(0, "phenix.molprobity")
        reason += " - VALIDATION REQUIRED BEFORE STOPPING"
    elif map_cc and map_cc > 0.75:
        reason += " - consider stopping (good CC)"
    
    # Only add STOP if validation not required
    if not validation_required_before_stop:
        valid.append("STOP")
    
    return {
        "state": "cryoem_refined",
        "experiment_type": "cryoem",
        "valid_programs": valid,
        "reason": reason,
        "conditions": conditions,
        "automation_path": automation_path,
        "resolution": resolution,
    }


# =============================================================================
# MAIN DETECTION FUNCTION
# =============================================================================

def detect_workflow_state(history, available_files, analysis=None, maximum_automation=True):
    """
    Determine current workflow state based on history and files.
    
    Args:
        history: List of cycle records from client
        available_files: List of available file paths
        analysis: Current log analysis dict (optional)
        maximum_automation: If True, use fully automated cryo-EM path
        
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
    
    # Analyze history first (needed for experiment type detection)
    history_info = _analyze_history(history)
    
    # Determine experiment type (considers both files and history)
    experiment_type = _detect_experiment_type(files, history_info)
    
    # Detect state based on experiment type
    if experiment_type == "cryoem":
        state = _detect_cryoem_state(files, history_info, analysis, maximum_automation)
    else:
        state = _detect_xray_state(files, history_info, analysis)
    
    # Add categorized files to state for use by command builder
    state["categorized_files"] = files
    
    return state


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
    
    # Build helpful error message
    all_known_programs = [
        "phenix.xtriage", "phenix.mtriage", "phenix.predict_and_build",
        "phenix.phaser", "phenix.refine", "phenix.real_space_refine",
        "phenix.ligandfit", "phenix.pdbtools", "phenix.dock_in_map",
        "phenix.autobuild", "phenix.autosol", "phenix.process_predicted_model",
    ]
    
    if chosen_program in all_known_programs:
        error = (
            "Program '%s' is not valid in state '%s'. "
            "Valid programs: %s. "
            "Reason: %s"
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
    
    # Add notes about conditions
    if workflow_state.get("conditions"):
        lines.append("")
        lines.append("Conditional availability:")
        for prog, condition in workflow_state["conditions"].items():
            lines.append("  - %s: requires %s" % (prog, condition))
    
    # Special instructions for stepwise cryo-EM
    if (workflow_state["experiment_type"] == "cryoem" and 
        workflow_state.get("automation_path") == "stepwise"):
        
        if workflow_state["state"] == "cryoem_analyzed":
            lines.append("")
            lines.append("NOTE (Stepwise mode): Use predict_and_build with strategy: {\"stop_after_predict\": true}")
    
    return "\n".join(lines)
