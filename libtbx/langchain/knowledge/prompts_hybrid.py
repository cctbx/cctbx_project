"""
Hybrid Planning Prompts for PHENIX AI Agent.

This module constructs the prompts for the LLM planner, including:
- Workflow state constraints
- Metrics trend information
- Error recovery guidance
- File inventory
- Tiered recommendations (ranked programs, defaults, suggestions)
"""

from __future__ import absolute_import, division, print_function
import json
import os


def format_recommendations_for_prompt(workflow_state):
    """
    Format the tiered recommendations from workflow_state for the LLM prompt.
    
    Args:
        workflow_state: Output from detect_workflow_state() with new recommendation fields
        
    Returns:
        str: Formatted recommendations section for prompt
    """
    if not workflow_state:
        return ""
    
    sections = []
    
    # === RANKED PROGRAMS ===
    ranked_programs = workflow_state.get("ranked_programs", [])
    if ranked_programs:
        sections.append("**Ranked Programs (in order of recommendation):**")
        for item in ranked_programs:
            if item.get("rank"):
                sections.append("  %d. %s - %s" % (item["rank"], item["program"], item["reason"]))
    
    # === RECOMMENDED STRATEGY ===
    recommended_strategy = workflow_state.get("recommended_strategy", {})
    if recommended_strategy:
        sections.append("")
        sections.append("**Recommended Strategy (defaults will be applied if you don't specify):**")
        for key, info in recommended_strategy.items():
            value = info.get("value")
            reason = info.get("reason", "")
            sections.append("  - %s: %s" % (key, value))
            if reason:
                sections.append("    Reason: %s" % reason)
        sections.append("")
        sections.append("You may override these defaults. If you do, your choice will be accepted but logged.")
    
    # === HARD CONSTRAINTS ===
    hard_constraints = workflow_state.get("hard_constraints", [])
    if hard_constraints:
        sections.append("")
        sections.append("**Hard Constraints (cannot override):**")
        for constraint in hard_constraints:
            sections.append("  - %s" % constraint)
    
    # === SOFT SUGGESTIONS ===
    soft_suggestions = workflow_state.get("soft_suggestions", [])
    if soft_suggestions:
        sections.append("")
        sections.append("**Suggestions:**")
        for suggestion in soft_suggestions:
            sections.append("  - %s" % suggestion)
    
    # === THRESHOLDS INFO ===
    thresholds = workflow_state.get("thresholds", {})
    if thresholds:
        sections.append("")
        sections.append("**Current Thresholds (resolution-dependent):**")
        for key, value in thresholds.items():
            # Format nicely
            nice_key = key.replace("_", " ").replace("rfree", "R-free")
            sections.append("  - %s: %.2f" % (nice_key, value))
    
    if sections:
        return "\n### GRAPH RECOMMENDATIONS\n\n" + "\n".join(sections) + "\n"
    
    return ""

SYSTEM_PROMPT = """You are an expert crystallographer and cryo-EM analyst.
Your goal is to determine the next processing step to improve the model or map.

### WORKFLOW RULES

**The workflow state machine tells you EXACTLY what programs are valid.**
You will receive a "WORKFLOW STATE" section that specifies:
1. The current state (e.g., "xray_refined", "cryoem_analyzed")
2. Valid programs for that state
3. Why those programs are valid

**YOU MUST OBEY THE WORKFLOW STATE.** If a program is not in the valid list,
choosing it will cause a validation error and waste a retry.

### GRAPH RECOMMENDATIONS

You will also receive "GRAPH RECOMMENDATIONS" which include:

1. **Ranked Programs**: Programs sorted by recommendation priority. You may choose
   any valid program, but the ranking reflects crystallographic best practices.

2. **Recommended Strategy**: Default values that will be automatically applied.
   You can override these by specifying different values in your strategy.
   If you override defaults, your choice will be accepted but logged as an override.

3. **Hard Constraints**: Rules that cannot be overridden (e.g., R-free flag handling).

4. **Suggestions**: Soft guidance based on current context (e.g., consider validation).

### INPUT DATA
You have access to:
1. **Workflow State**: Current state and valid programs
2. **Graph Recommendations**: Ranked programs and strategy defaults
3. **Metrics Trend**: R-free or CC progression over recent cycles
4. **History**: Results of previous cycles
5. **Analysis**: Metrics from the most recent log
6. **Available Files**: A strict list of files you can use

### OUTPUT FORMAT
You must output a SINGLE JSON object matching this schema:
{
    "program": "phenix.refine",
    "reasoning": "Detailed explanation of why this step is appropriate.",
    "files": {
        "model": "PHASER.1.pdb",
        "data": "data.mtz"
    },
    "strategy": {
        "use_simulated_annealing": true,
        "nproc": 4
    },
    "stop": false,
    "stop_reason": null
}

### PROGRAM REFERENCE

**phenix.xtriage** - Data quality analysis (X-ray)
  Files: {data: .mtz/.sca/.hkl}
  Use: ALWAYS run first on X-ray data
  
**phenix.mtriage** - Map quality analysis (Cryo-EM)
  Files: {full_map: .mrc/.ccp4, half_map: [.mrc/.ccp4, .mrc/.ccp4]}
  Use: ALWAYS run first on cryo-EM maps
  NOTE on map types:
    - full_map: A single full reconstruction (e.g., full_map.ccp4)
    - half_map: A pair of half-maps (e.g., map_1.ccp4, map_2.ccp4 or map_a.ccp4, map_b.ccp4)
    - If BOTH are available, use both: full_map=X half_map=Y half_map=Z
    - If ONLY half-maps available, use: half_map=Y half_map=Z
    - If ONLY full map available, use: full_map=X

**phenix.resolve_cryo_em** - Optimize cryo-EM map clarity (RECOMMENDED)
  Files: {half_map: [.mrc/.ccp4, .mrc/.ccp4], sequence: .fa/.seq (optional)}
  Use: Creates optimized/sharpened full map from half-maps
  Output: *_sharpened.ccp4 (full map suitable for real_space_refine)
  Try this FIRST for map improvement before real_space_refine

**phenix.map_sharpening** - Sharpen cryo-EM map
  Files: {half_map: [.mrc/.ccp4, .mrc/.ccp4], sequence: .fa/.seq (optional)}
  Use: Alternative to resolve_cryo_em for map sharpening
  Output: *_sharpened.ccp4 (full map suitable for real_space_refine)
  
**phenix.predict_and_build** - AlphaFold prediction + building
  Files: {sequence: .fa/.seq/.dat, data: .mtz/.sca (X-ray), full_map: .mrc/.ccp4 (cryo-EM), half_map: [.mrc, .mrc] (cryo-EM)}
  Strategy: {resolution: N, stop_after_predict: true/false}
  IMPORTANT: Set resolution if building (get from xtriage/mtriage)
  For cryo-EM: Use full_map= for single map, or half_map= twice for half-maps, or both
  
**phenix.process_predicted_model** - Prepare AlphaFold model for MR
  Files: {model: .pdb}
  Use: After predict_and_build (X-ray or stepwise cryo-EM) to prepare model for phaser
  
**phenix.phaser** - Molecular replacement
  Files: {data: .mtz/.sca/.hkl, model: .pdb}
  REQUIRES: A model file

**phenix.autosol** - Experimental phasing (SAD/MAD)
  Files: {data: .mtz/.sca/.hkl, sequence: .fa/.seq/.dat}
  Strategy: {atom_type: "Se"/"S"/"Zn"/etc}
  Use: When xtriage reports useful anomalous signal (to ~4Å or better)
  REQUIRES: Data with anomalous signal (I+/I- or F+/F-) AND sequence file
  IMPORTANT: User must specify atom_type (anomalous scatterer) if not obvious
  After autosol: run phenix.autobuild to complete the model

**phenix.refine** - Crystallographic refinement
  Files: {model: .pdb, data: .mtz/.sca/.hkl}
  Strategy: {generate_rfree_flags: true/false, add_waters: true/false, use_simulated_annealing: true/false, nproc: N}
  
  R-FREE FLAGS:
    - FIRST refinement after phaser: set generate_rfree_flags=true
    - ALL subsequent refinements: do NOT set generate_rfree_flags (use existing flags in MTZ)
  
  WATER BUILDING (ordered_solvent):
    - Do NOT add waters on first 1-2 refinement cycles
    - Add waters (add_waters=true) only when: R-free < 0.35 AND resolution <= 3.0Å AND 2+ refine cycles done
  
  SIMULATED ANNEALING:
    - Do NOT use for initial/normal refinement
    - Only use simulated_annealing=true if R-free > 0.40 AND model is stuck

**phenix.real_space_refine** - Cryo-EM refinement
  Files: {model: .pdb, map: .mrc/.ccp4 (FULL MAP REQUIRED)}
  Strategy: {resolution: N} (REQUIRED - get from mtriage)
  NOTE: Cannot use half-maps directly - requires a full/combined map
  
**phenix.ligandfit** - Fit ligand into density
  Files: {model: .pdb, data: .mtz/.sca OR map: .mrc, ligand: .cif/.pdb}
  REQUIRES: Ligand coordinates file
  For cryo-EM: use map_in=file.mrc

**phenix.pdbtools** - Combine PDB files
  Files: {model: .pdb (multiple)}
  Use: Combine protein + ligand after ligandfit

**phenix.dock_in_map** - Dock model into cryo-EM map
  Files: {model: .pdb, map: .mrc}

**phenix.autobuild** - Build model into electron density
  Files: {data: .mtz (PHASED), sequence: .fa/.seq/.dat, model: .pdb (RECOMMENDED)}
  IMPORTANT: Always include model= with the latest refined PDB for best results
  REQUIRES: Phased data with map coefficients (from refine or phaser output)
  If resolution < 2.0Å, set resolution=2.0 to speed up

### VALIDATION PROGRAMS

**phenix.molprobity** - Geometry validation (X-ray and cryo-EM)
  Files: {model: .pdb}
  Good values: Ramachandran outliers <0.5%, rotamer <1%, clashscore <4, MolProbity <1.5
  
**phenix.model_vs_data** - Model vs X-ray data validation
  Files: {model: .pdb, data: .mtz/.sca/.hkl}
  
**phenix.validation_cryoem** - Model vs cryo-EM map validation  
  Files: {model: .pdb, map: .mrc}

**phenix.holton_geometry_validation** - Detailed geometry scoring
  Files: {model: .pdb}
  Overall geometry energy should be close to expected (~67)

### STOP CONDITIONS

Set "stop": true when:
- **X-ray**: R-free < 0.25 (or dynamic target based on resolution)
- **Cryo-EM**: Map-model CC > 0.70
- **Plateau**: Metrics trend shows <0.3% improvement for 3+ cycles
- **Stuck**: Same error repeated 3+ times with no progress
- **No valid programs**: Workflow state shows only STOP is valid

### SPECIAL CASES

**TWINNING**: If xtriage reports twinning (Britton alpha > 20%), include twin_law in refine
  Example: strategy: {"twin_law": "-h,-k,l"}

**HIGH RESOLUTION** (< 1.5Å): Consider anisotropic_adp=true in refine

**VALIDATION**: Run validation programs when:
  - R-free target reached (to confirm model quality)
  - After 3+ refine cycles (to check for problems)
  - If R-free is good but model looks wrong

### IMPORTANT RULES

1. **Trust the workflow state** - it knows what's valid
2. **Watch the metrics trend** - stop if plateau detected
3. **Don't repeat failing commands** - change strategy or try different program
4. **Always set resolution for predict_and_build** if building
5. **Files must exist** - only use files from the inventory
"""


def get_planning_prompt(history, analysis, available_files, previous_attempts=None, 
                        user_advice="", metrics_trend=None, workflow_state=None):
    """
    Constructs the system and user messages for the LLM planner.
    
    Args:
        history: List of previous cycle records
        analysis: Current log analysis dict
        available_files: List of available files
        previous_attempts: List of failed attempts in current cycle
        user_advice: User-provided instructions
        metrics_trend: Output from analyze_metrics_trend()
        workflow_state: Output from detect_workflow_state()
        
    Returns:
        tuple: (system_message, user_message)
    """
    # Start with the provided available files
    files_list = list(available_files)
    
    # Accumulate output files from history
    if history:
        for h in history:
            if isinstance(h, dict):
                output_files = h.get('output_files', [])
                for f in output_files:
                    if f and f not in files_list:
                        files_list.append(f)
    
    # Also add output files from current analysis
    if analysis and isinstance(analysis, dict):
        current_outputs = analysis.get('output_files', [])
        for f in current_outputs:
            if f and f not in files_list:
                files_list.append(f)
    
    # Analyze all available files
    # X-ray data files can be .mtz, .sca (scalepack), .hkl, .sdf
    xray_data_extensions = ('.mtz', '.sca', '.hkl', '.sdf')
    
    pdb_files = [f for f in files_list if f.endswith('.pdb')]
    data_files = [f for f in files_list if f.endswith(xray_data_extensions)]
    seq_files = [f for f in files_list if f.endswith(('.fa', '.fasta', '.seq', '.dat'))]
    map_files = [f for f in files_list if f.endswith(('.mrc', '.ccp4', '.map'))]
    
    # Categorize map files into full maps and half maps
    # Half-map naming patterns: half_map_1, half1, map_1, map_2, _a, _b, etc.
    def is_half_map(filepath):
        import re
        basename = os.path.basename(filepath).lower()
        # Explicit half-map patterns
        if 'half' in basename:
            return True
        # Match patterns like: map_1, map_2, map-1, map-2, map1, map2, _a, _b
        if re.search(r'[_-][12ab]\.', basename):
            return True
        if re.search(r'map[_-]?[12]', basename):
            return True
        return False
    
    full_map_files = [f for f in map_files if not is_half_map(f)]
    half_map_files = [f for f in map_files if is_half_map(f)]
    
    # CIF files need careful categorization:
    # - Model CIFs: from refinement (contain 'refine' in name) - these are NOT ligands
    # - Ligand CIFs: restraint files for ligands (typically small, have 'lig' in name, or are user-provided)
    cif_files = [f for f in files_list if f.endswith('.cif')]
    model_cif_files = []
    ligand_cif_files = []
    for f in cif_files:
        basename = os.path.basename(f).lower()
        # Model CIFs from refinement - NOT ligands
        if 'refine' in basename or 'autobuild' in basename or 'autosol' in basename or 'overall_best' in basename:
            model_cif_files.append(f)
        else:
            # Assume other CIFs are ligand restraint files
            ligand_cif_files.append(f)
    
    # Build file summary
    file_summary = []
    if pdb_files:
        file_summary.append("MODELS (.pdb): %s" % ", ".join(pdb_files))
    else:
        file_summary.append("MODELS (.pdb): NONE AVAILABLE")
    
    # Include model CIFs with models
    if model_cif_files:
        file_summary.append("MODEL CIFs (from refinement): %s" % ", ".join(model_cif_files))
        
    if data_files:
        file_summary.append("DATA (.mtz/.sca/.hkl): %s" % ", ".join(data_files))
    else:
        file_summary.append("DATA (.mtz/.sca/.hkl): NONE")
        
    if seq_files:
        file_summary.append("SEQUENCE (.fa/.seq/.dat): %s" % ", ".join(seq_files))
    else:
        file_summary.append("SEQUENCE: NONE")
    
    # Show map files with proper categorization for cryo-EM
    if full_map_files or half_map_files:
        if full_map_files:
            file_summary.append("FULL MAPS (.mrc/.ccp4): %s" % ", ".join(full_map_files))
        if half_map_files:
            file_summary.append("HALF MAPS (.mrc/.ccp4): %s" % ", ".join(half_map_files))
            if len(half_map_files) == 2:
                file_summary.append("  -> Use with mtriage: half_map=%s half_map=%s" % (half_map_files[0], half_map_files[1]))
            elif len(half_map_files) > 2:
                file_summary.append("  -> WARNING: More than 2 half maps detected, select the matching pair")
        
    if ligand_cif_files:
        file_summary.append("LIGAND RESTRAINTS (.cif): %s" % ", ".join(ligand_cif_files))

    # === WORKFLOW STATE SECTION ===
    workflow_section = ""
    if workflow_state:
        valid_list = ", ".join(workflow_state.get("valid_programs", []))
        workflow_section = """
### WORKFLOW STATE: %s
Experiment type: %s
%s

**VALID PROGRAMS:** %s

You MUST choose from the valid programs above, or set "stop": true.
""" % (
            workflow_state.get("state", "unknown"),
            workflow_state.get("experiment_type", "unknown"),
            workflow_state.get("reason", ""),
            valid_list
        )
        
        # Add automation path note for cryo-EM
        if workflow_state.get("automation_path"):
            workflow_section += "\nAutomation path: %s\n" % workflow_state["automation_path"]
            
            if workflow_state["automation_path"] == "stepwise":
                if "predict_and_build" in valid_list:
                    workflow_section += "NOTE: Use predict_and_build with strategy: {\"stop_after_predict\": true}\n"
        
        # Add recommendations section if available
        recommendations = format_recommendations_for_prompt(workflow_state)
        if recommendations:
            workflow_section += recommendations

    # === METRICS TREND SECTION ===
    metrics_section = ""
    if metrics_trend:
        metrics_section = "\n### REFINEMENT PROGRESS\n"
        metrics_section += metrics_trend.get("trend_summary", "No trend data")
        
        if metrics_trend.get("consecutive_refines", 0) > 0:
            metrics_section += "\nConsecutive X-ray refinement cycles: %d" % metrics_trend["consecutive_refines"]
        if metrics_trend.get("consecutive_rsr", 0) > 0:
            metrics_section += "\nConsecutive real-space refinement cycles: %d" % metrics_trend["consecutive_rsr"]
        
        if metrics_trend.get("recommendation") == "consider_stopping":
            metrics_section += "\n\n*** DIMINISHING RETURNS - Consider stopping or trying different strategy ***"
        
        if metrics_trend.get("should_stop"):
            metrics_section += "\n\n*** STOP RECOMMENDED: %s ***\nSet \"stop\": true in your response." % (
                metrics_trend.get("reason", "Unknown reason")
            )

    # === HISTORY SECTION ===
    history_str = ""
    last_error = None
    last_failed_command = None
    last_failed_program = None
    
    if history:
        for h in history[-5:]:  # Last 5 cycles
            if isinstance(h, dict):
                cycle_num = h.get('cycle_number', '?')
                program = h.get('program', 'unknown')
                result = h.get('result', h.get('summary', 'unknown'))
                error = h.get('error', '')
                command = h.get('command', '')
            elif isinstance(h, str):
                cycle_num = '?'
                program = 'unknown'
                result = h
                error = ''
                command = ''
                if 'FAILED' in h or 'Sorry' in h or 'Error' in h:
                    error = h
            else:
                continue
            
            # Detect failure in result string
            if str(result).startswith("FAILED") or "Sorry" in str(result):
                error = result
            
            if len(str(result)) > 100:
                result = str(result)[:100] + "..."
            
            line = "- Cycle %s: %s" % (cycle_num, program)
            if error:
                line += " [ERROR: %s]" % str(error)[:80]
                last_error = str(error)
                last_failed_command = command
                last_failed_program = program
            else:
                line += " -> %s" % result
            history_str += line + "\n"
    else:
        history_str = "No previous history. This is the FIRST cycle."

    # === RETRY CONTEXT ===
    retry_msg = ""
    if previous_attempts:
        last = previous_attempts[-1]
        retry_msg = """
!!! VALIDATION FAILED - RETRY REQUIRED !!!
Error: %s
Command attempted: %s

YOU MUST FIX THIS. Try a DIFFERENT approach:
- Check if required files exist
- Choose a valid program from the workflow state
- Use different files
""" % (
            last.get('error', 'unknown'),
            last.get('command', 'none')
        )

    # === RUNTIME ERROR CONTEXT ===
    runtime_error_msg = ""
    if last_error and not previous_attempts:
        is_phil_error = any(x in last_error.lower() for x in [
            "phil parameter", "unknown parameter", "unrecognized", "invalid keyword",
            "not a valid", "syntax error", "unexpected"
        ])
        
        if is_phil_error:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED - PHIL SYNTAX ERROR !!!
Program: %s
Error: "%s"

THIS IS A SYNTAX ERROR - DO NOT SWITCH PROGRAMS!
Retry with corrected parameter names.
""" % (last_failed_program or "unknown", last_error[:200])
        else:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED !!!
Program: %s
Error: "%s"

You must adapt:
1. Add/change strategy parameters
2. Try a different valid program
3. Use different files
""" % (last_failed_program or "unknown", last_error[:200])

    # === RESOLUTION HINT ===
    resolution_hint = ""
    if analysis and analysis.get("resolution"):
        resolution_hint = "\nNOTE: Resolution %.2f found. Use this for predict_and_build: \"strategy\": {\"resolution\": %.1f}" % (
            analysis["resolution"], analysis["resolution"]
        )

    # === BUILD USER MESSAGE ===
    user_msg = """
%s
%s
### CURRENT STATUS
User Advice: %s

Log Analysis: %s

### FILE INVENTORY
%s

### HISTORY
%s
%s
%s
%s
Based on the workflow state and available files, what is the next step?
Output JSON only.
""" % (
        workflow_section,
        metrics_section,
        user_advice or "None provided",
        json.dumps(analysis, indent=2) if analysis else "No analysis yet (first run)",
        "\n".join(file_summary),
        history_str,
        runtime_error_msg,
        retry_msg,
        resolution_hint
    )

    return SYSTEM_PROMPT, user_msg
