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

# Import YAML-based program registry
from libtbx.langchain.agent.program_registry import ProgramRegistry


def generate_program_reference_from_yaml():
    """
    Generate PROGRAM REFERENCE section from YAML definitions.

    This can be used to replace the hardcoded program reference
    in SYSTEM_PROMPT when ready for full migration.

    Returns:
        str: Formatted program reference section
    """
    registry = ProgramRegistry(use_yaml=True)

    sections = ["### PROGRAM REFERENCE (Generated from YAML)\n"]

    # Group by category
    categories = {
        "analysis": "DATA ANALYSIS",
        "model_building": "MODEL BUILDING",
        "refinement": "REFINEMENT",
        "ligand": "LIGAND HANDLING",
        "map_optimization": "MAP OPTIMIZATION",
        "validation": "VALIDATION",
        "utility": "UTILITY",
    }

    programs_by_category = {}
    for prog_name in registry.get_all_programs():
        prog = registry.get_program(prog_name)
        category = prog.get("category", "other")
        if category not in programs_by_category:
            programs_by_category[category] = []
        programs_by_category[category].append((prog_name, prog))

    for cat_key, cat_name in categories.items():
        if cat_key not in programs_by_category:
            continue

        sections.append("\n**%s**\n" % cat_name)

        for prog_name, prog in programs_by_category[cat_key]:
            desc = prog.get("description", "")
            sections.append("**%s** - %s" % (prog_name, desc))

            # Format inputs
            inputs = prog.get("inputs", {})
            required = inputs.get("required", {})
            optional = inputs.get("optional", {})

            if required:
                file_parts = []
                for slot, slot_def in required.items():
                    exts = slot_def.get("extensions", [])
                    ext_str = "/".join(exts) if exts else "file"
                    file_parts.append("%s: %s" % (slot, ext_str))
                sections.append("  Files: {%s}" % ", ".join(file_parts))

            # Add hints
            hints = prog.get("hints", [])
            for hint in hints[:2]:  # First 2 hints
                sections.append("  - %s" % hint)

            sections.append("")

    return "\n".join(sections)


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

### USER ADVICE (HIGHEST PRIORITY)

**If the user provides advice, you MUST follow it.**
User advice may include:
- Specific programs to use or avoid
- Strategy parameters (e.g., resolution, number of cycles)
- Workflow preferences (e.g., "use experimental phasing", "try molecular replacement first")
- Troubleshooting hints (e.g., "the model has wrong hand", "try different space group")

**CRITICAL: User advice overrides default program preferences.**

Examples of user advice that MUST be followed:
- "use experimental phasing" or "use SAD phasing" → Choose phenix.autosol (NOT predict_and_build)
- "use molecular replacement" → Choose phenix.phaser
- "skip prediction" → Do NOT use predict_and_build

When user advice conflicts with workflow recommendations:
1. User advice takes precedence over soft recommendations and "preferred" markers
2. User advice cannot override hard constraints (e.g., file availability)
3. If user requests an invalid program, explain why and suggest alternatives

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

**phenix.predict_and_build** - AlphaFold prediction + MR + building
  Files: {sequence: .fa/.seq/.dat, data: .mtz/.sca (X-ray), full_map: .mrc/.ccp4 (cryo-EM), half_map: [.mrc, .mrc] (cryo-EM)}
  Strategy: {resolution: N, stop_after_predict: true/false, rebuilding_strategy: "Quick"/"Standard"}
  IMPORTANT: Set resolution if building (get from xtriage/mtriage)
  NOTE: By default (stop_after_predict=False), this runs the FULL workflow: prediction → molecular replacement → model building
  NOTE: Only set stop_after_predict=True for cryo-EM stepwise workflow where you want just the predicted model
  WARNING: This is NOT a density modification tool! Do NOT use for "density modification" - use phenix.autobuild_denmod instead
  For cryo-EM maps:
    - If you have a FULL MAP: use full_map=filename.ccp4
    - If you ONLY have HALF MAPS: use half_map=file1.ccp4 half_map=file2.ccp4 (NO full_map!)
    - NEVER use a half-map as full_map - half-maps have names like _1.ccp4, _2.ccp4, half1, half2

**phenix.autobuild_denmod** - X-ray density modification (map improvement)
  Files: {data: .mtz (refined MTZ with phases), sequence: .fa/.seq, model: .pdb (optional)}
  Output: overall_best_denmod_map_coeffs.mtz with FWT/PHFWT map coefficients
  Use: Before ligand fitting to improve map quality
  NOTE: This runs autobuild with maps_only=True (no model building, just map improvement)
  NOTE: For ligandfit using this output, set file_info.input_labels="FWT PHFWT"

**phenix.process_predicted_model** - Prepare AlphaFold model for MR
  Files: {model: .pdb}
  Use: After predict_and_build (X-ray or stepwise cryo-EM) to prepare model for phaser

**phenix.phaser** - Molecular replacement
  Files: {data: .mtz/.sca/.hkl, model: .pdb, sequence: .seq/.fa (for composition)}
  REQUIRES: A model file

  Include ALL sequence files for correct solvent content calculation.
  Example: phenix.phaser data.mtz model.pdb seq1.seq seq2.seq phaser.mode=MR_AUTO

**phenix.autosol** - Experimental phasing (SAD/MAD)
  Files: {data: .mtz/.sca/.hkl, sequence: .fa/.seq/.dat}
  Strategy fields: {atom_type: "Se", additional_atom_types: "S", wavelength: 0.9792, sites: 5, resolution: 2.5}
  Use: When xtriage reports useful anomalous signal (to ~4Å or better)
  REQUIRES: Data with anomalous signal (I+/I- or F+/F-) AND sequence file
  CRITICAL - Include these in strategy if user provides them:
    - atom_type: Primary anomalous scatterer (Se, S, Zn, etc) - ONE type only
    - additional_atom_types: Extra atom types to search (e.g., "S" if also looking for sulfur)
    - wavelength: X-ray wavelength in Angstroms (e.g., 0.9792)
    - sites: Expected number of anomalous sites (e.g., 5)
    - resolution: High resolution limit in Angstroms (e.g., 2.5)
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

**phenix.polder** - Calculate Polder omit maps to evaluate ligand/residue placement
  Files: {model: .pdb, data: .mtz}
  Strategy: {selection: "chain A and resseq 88"} (REQUIRED - specify atoms to evaluate)
  IMPORTANT: Polder works with STANDARD MTZ data (Fobs + R-free flags). It does NOT need
  pre-calculated map coefficients or phases - it calculates the omit maps internally.
  Use: Verify ligand placement by checking if density exists for the omitted atoms.
  Output: CC values and conclusion ("likely to show", "inconclusive", "unlikely to show")

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


def _format_directives_for_prompt(directives):
    """
    Format directives for inclusion in the LLM prompt.

    Shows the extracted directives so the LLM knows what was understood
    from the user's advice. This helps the LLM make consistent decisions.

    Args:
        directives: Directives dict

    Returns:
        str: Formatted section for prompt, or empty string if no directives
    """
    if not directives:
        return ""

    lines = []
    lines.append("### EXTRACTED DIRECTIVES (from user advice)")
    lines.append("The following directives were extracted and will be enforced automatically.")
    lines.append("You should be aware of them when planning, but don't need to repeat them in your strategy.")
    lines.append("")

    # Program settings
    prog_settings = directives.get("program_settings", {})
    if prog_settings:
        lines.append("**Program Settings:**")
        for prog, settings in prog_settings.items():
            if settings:
                settings_str = ", ".join("%s=%s" % (k, v) for k, v in settings.items())
                lines.append("- %s: %s" % (prog, settings_str))
        lines.append("")

    # Stop conditions
    stop_cond = directives.get("stop_conditions", {})
    if stop_cond:
        lines.append("**Stop Conditions:**")
        if "after_cycle" in stop_cond:
            lines.append("- Stop after cycle %d" % stop_cond["after_cycle"])
        if "after_program" in stop_cond:
            after_prog = stop_cond["after_program"]
            lines.append("- Stop after %s completes" % after_prog)
            # Add explicit guidance to run the program - make it very clear
            lines.append("- **CRITICAL: You MUST run %s before stopping. If it's in VALID PROGRAMS, choose it NOW.**" % after_prog)
            lines.append("- Do NOT keep running refinement cycles - run %s instead!" % after_prog)
        if "max_refine_cycles" in stop_cond:
            lines.append("- Maximum %d refinement cycles" % stop_cond["max_refine_cycles"])
        if "r_free_target" in stop_cond:
            lines.append("- Target R-free: %.3f" % stop_cond["r_free_target"])
        if "map_cc_target" in stop_cond:
            lines.append("- Target map CC: %.2f" % stop_cond["map_cc_target"])
        if stop_cond.get("skip_validation"):
            lines.append("- Validation can be skipped before stopping")
        lines.append("")

    # Workflow preferences
    workflow_prefs = directives.get("workflow_preferences", {})
    if workflow_prefs:
        lines.append("**Workflow Preferences:**")
        if workflow_prefs.get("skip_programs"):
            lines.append("- Skip: %s" % ", ".join(workflow_prefs["skip_programs"]))
        if workflow_prefs.get("prefer_programs"):
            lines.append("- Prefer: %s" % ", ".join(workflow_prefs["prefer_programs"]))
        if workflow_prefs.get("use_experimental_phasing"):
            lines.append("- Use experimental phasing (SAD/MAD)")
        if workflow_prefs.get("use_molecular_replacement"):
            lines.append("- Use molecular replacement")
        lines.append("")

    # Constraints
    constraints = directives.get("constraints", [])
    if constraints:
        lines.append("**Additional Constraints:**")
        for c in constraints:
            lines.append("- %s" % c)
        lines.append("")

    return "\n".join(lines)


def get_planning_prompt(history, analysis, available_files, previous_attempts=None,
                        user_advice="", metrics_trend=None, workflow_state=None,
                        directives=None, best_files=None):
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
        directives: Structured directives extracted from user advice
        best_files: Dict of {category: path} for best files to use

    Returns:
        tuple: (system_message, user_message)
    """
    # Helper to escape % characters in user-provided strings to avoid format errors
    def escape_percent(s):
        if s is None:
            return ""
        return str(s).replace("%", "%%")

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
        else:
            file_summary.append("FULL MAPS (.mrc/.ccp4): NONE - do NOT use half-maps as full_map!")
        if half_map_files:
            file_summary.append("HALF MAPS (.mrc/.ccp4): %s" % ", ".join(half_map_files))
            if len(half_map_files) == 2:
                file_summary.append("  -> For predict_and_build/mtriage with half-maps only: half_map=%s half_map=%s" % (half_map_files[0], half_map_files[1]))
            elif len(half_map_files) > 2:
                file_summary.append("  -> WARNING: More than 2 half maps detected, select the matching pair")

    if ligand_cif_files:
        file_summary.append("LIGAND RESTRAINTS (.cif): %s" % ", ".join(ligand_cif_files))

    # Add RECOMMENDED FILES section if best_files are available
    # This tells the LLM which files to use for iterative workflows
    if best_files:
        recommended = []
        if best_files.get("model"):
            recommended.append("**USE THIS MODEL:** %s" % os.path.basename(best_files["model"]))
        if best_files.get("mtz"):
            recommended.append("**USE THIS DATA:** %s" % os.path.basename(best_files["mtz"]))
        if best_files.get("map"):
            recommended.append("**USE THIS MAP:** %s" % os.path.basename(best_files["map"]))

        if recommended:
            file_summary.append("")
            file_summary.append(">>> RECOMMENDED FILES FOR NEXT STEP <<<")
            file_summary.extend(recommended)
            file_summary.append("(These are the latest/best files from previous cycles - always use these for refinement!)")

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

            # Only show stop_after_predict guidance for cryo-EM stepwise workflow
            # For X-ray, predict_and_build should run the full workflow (prediction + MR + building)
            experiment_type = workflow_state.get("experiment_type", "")
            if (workflow_state["automation_path"] == "stepwise" and
                experiment_type == "cryoem" and
                "predict_and_build" in valid_list):
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

            line = "- Cycle %s: %s" % (escape_percent(cycle_num), escape_percent(program))
            if error:
                line += " [ERROR: %s]" % escape_percent(str(error)[:80])
                last_error = str(error)
                last_failed_command = command
                last_failed_program = program
            else:
                line += " -> %s" % escape_percent(result)
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
            escape_percent(last.get('error', 'unknown')),
            escape_percent(last.get('command', 'none'))
        )

    # === RUNTIME ERROR CONTEXT ===
    runtime_error_msg = ""
    if last_error and not previous_attempts:
        last_error_lower = last_error.lower()

        # Check for specific error types
        is_phil_error = any(x in last_error_lower for x in [
            "phil parameter", "unknown parameter", "unrecognized", "not recognized",
            "invalid keyword", "not a valid", "syntax error", "unexpected",
            "arguments are not recognized", "ambiguous parameter"
        ])

        is_rfree_error = any(x in last_error_lower for x in [
            "r-free", "rfree", "r_free", "free_flags", "test set",
            "no test reflections", "fraction of test"
        ])

        is_resolution_error = any(x in last_error_lower for x in [
            "resolution", "d_min", "high_resolution"
        ])

        is_file_error = any(x in last_error_lower for x in [
            "file not found", "no such file", "cannot open", "does not exist",
            "missing file", "input file"
        ])

        # Escape % in error messages to prevent format string issues
        safe_error = escape_percent(last_error[:200])
        safe_program = escape_percent(last_failed_program or "unknown")

        if is_phil_error:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED - PHIL SYNTAX ERROR !!!
Program: %s
Error: "%s"

THIS IS A SYNTAX ERROR - DO NOT SWITCH PROGRAMS!
Retry with corrected parameter names.
""" % (safe_program, safe_error)

        elif is_rfree_error:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED - R-FREE FLAG ISSUE !!!
Program: %s
Error: "%s"

THIS IS AN R-FREE FLAG ERROR - DO NOT SWITCH PROGRAMS!
The refinement command already includes xray_data.r_free_flags.generate=True
which should auto-generate R-free flags. Retry refinement - it should work now.
If using a different MTZ file, ensure it has reflection data.
""" % (safe_program, safe_error)

        elif is_resolution_error:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED - RESOLUTION ISSUE !!!
Program: %s
Error: "%s"

THIS IS A RESOLUTION ERROR - DO NOT SWITCH PROGRAMS!
Add or fix the resolution parameter in strategy.
""" % (safe_program, safe_error)

        elif is_file_error:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED - FILE NOT FOUND !!!
Program: %s
Error: "%s"

A required file was not found. Check the FILE INVENTORY and use only files that exist.
""" % (safe_program, safe_error)

        else:
            runtime_error_msg = """
!!! PREVIOUS CYCLE FAILED !!!
Program: %s
Error: "%s"

You must adapt:
1. Add/change strategy parameters
2. Try a different valid program
3. Use different files
""" % (safe_program, safe_error)

    # === RESOLUTION HINT ===
    resolution_hint = ""
    if analysis and analysis.get("resolution"):
        resolution_hint = "\nNOTE: Resolution %.2f found. Use this for predict_and_build: \"strategy\": {\"resolution\": %.1f}" % (
            analysis["resolution"], analysis["resolution"]
        )

    # === BUILD USER MESSAGE ===
    # User advice at the top for prominence
    user_advice_section = ""
    if user_advice and user_advice.strip() and user_advice.strip().lower() != "none provided":
        user_advice_section = """
### USER ADVICE (FOLLOW THIS)
%s

**IMPORTANT**: Extract any specific parameters from the user advice above (wavelength, atom type,
resolution, number of sites, etc.) and include them in your "strategy" field. For example:
- If user mentions wavelength 0.9792 → add to strategy: "wavelength": 0.9792
- If user mentions Se atoms → add to strategy: "atom_type": "Se"
- If user mentions additional S atoms → add to strategy: "additional_atom_types": "S"
- If user mentions 5 sites → add to strategy: "sites": 5
- If user mentions resolution 2.5 Å → add to strategy: "resolution": 2.5

""" % escape_percent(user_advice)

    # Directives section - show extracted structured directives
    directives_section = ""
    if directives:
        directives_section = _format_directives_for_prompt(directives)

    user_msg = """
%s%s%s
### CURRENT STATUS
%s
Log Analysis: %s

### FILE INVENTORY
%s

### HISTORY
%s
%s
%s
%s
Based on the workflow state, user advice, and available files, what is the next step?
Output JSON only.
""" % (
        user_advice_section,
        directives_section,
        workflow_section,
        metrics_section,
        json.dumps(analysis, indent=2) if analysis else "No analysis yet (first run)",
        "\n".join(file_summary),
        history_str,
        runtime_error_msg,
        retry_msg,
        resolution_hint
    )

    return SYSTEM_PROMPT, user_msg


# =============================================================================
# AGENT SESSION ASSESSMENT PROMPT
# =============================================================================

AGENT_SESSION_ASSESSMENT_PROMPT = """You are a senior crystallographer reviewing the results of an automated structure determination workflow.

Analyze the session summary below and provide a brief assessment covering:

1. **Input Data Quality**: What is the quality of the input data? (resolution, completeness, any issues)

2. **Goal and Strategy**: What was the user's goal and what strategy did the agent use to achieve it?

3. **Strategy Assessment**: Was the strategy appropriate for the data and goal? Was the goal achieved?

4. **Current Status**: What is the current state of the structure/analysis? Is it ready for further analysis or deposition?

5. **Next Steps**: What are appropriate next steps? (e.g., more refinement, ligand fitting, validation, deposition)

**IMPORTANT**: Check the "Stop Condition" in the session summary. If the session was a FOCUSED TASK or TUTORIAL
(e.g., "stop after xtriage", "stop after density modification", "this is a focused task"), then:
- The workflow was INTENTIONALLY limited - this is SUCCESS, not failure
- Do NOT suggest the workflow is "stalled", "incomplete", or "failed"
- Assess whether the specific focused task was completed successfully
- For "Next Steps", suggest what the user might do OUTSIDE this automated session

Keep your assessment concise (3-5 sentences per section). Focus on practical insights.

=== SESSION SUMMARY ===
{session_summary}
=== END SESSION SUMMARY ===

Provide your assessment in Markdown format with the headers above."""


def get_agent_session_assessment_prompt():
    """
    Get the prompt template for assessing an agent session.

    Returns:
        str: The assessment prompt template
    """
    return AGENT_SESSION_ASSESSMENT_PROMPT
