"""
Directive Extractor for PHENIX AI Agent.

This module extracts structured directives from natural language user advice
using an LLM. The directives are used to validate and augment LLM planning
decisions throughout the workflow.

The extraction happens once at session start, and the resulting directives
persist across all cycles.

Usage:
    from libtbx.langchain.agent.directive_extractor import extract_directives

    directives = extract_directives(
        user_advice="use resolution 3 in autosol but 2.5 in refine",
        provider="google",
        model="gemini-2.5-flash-lite"   # or None for the agent default
    )

    # Result:
    # {
    #   "program_settings": {
    #     "phenix.autosol": {"resolution": 3.0},
    #     "default": {"resolution": 2.5}
    #   }
    # }
"""

from __future__ import absolute_import, division, print_function

import json
import os
import re
import warnings


# v119.H1: central per-provider default models for the fallback LLM
# path.  Resolves model defaults via core/llm.py rather than
# hardcoding strings.  Fallback import path supports both PHENIX and
# standalone-tests invocation per AI Agent guideline section 3.
try:
    from libtbx.langchain.core.llm import (
        default_model_for_provider,
        # v119.H13
        normalize_ollama_openai_base_url,
        resolve_model_for_provider,
        _classify_provider_error,
    )
except ImportError:
    from core.llm import (
        default_model_for_provider,
        # v119.H13
        normalize_ollama_openai_base_url,
        resolve_model_for_provider,
        _classify_provider_error,
    )


# =============================================================================
# Stop-condition program pattern loader
# =============================================================================

_STOP_DIRECTIVE_PATTERNS_CACHE = None


def _strip_preprocessor_stop_condition(user_advice):
    """Strip LLM-generated structured headers from preprocessed advice.

    The advice preprocessor LLM creates structured output with headers
    like 'Stop Condition:', 'Primary Goal:', etc.  These are LLM
    summaries of the README, not user instructions, and can cause
    false-positive directive extraction.  For example:

    - 'Stop Condition: Stop after molecular replacement' causes a
      fabricated after_program directive (none of the original READMEs
      contain the word 'stop')
    - 'Primary Goal: Density modification of cryo-EM map' triggers
      the tutorial_patterns section to set after_program for denmod

    Real user commands ('stop after refinement', 'run xtriage') are
    natural prose and handled by later extraction patterns.  They
    never appear in the 'Header: value' format that only the
    preprocessor produces.

    'Stop Condition:' is always stripped (never real user input).
    'Goal:' / 'Primary Goal:' are only stripped when the text
    contains other preprocessor signature headers, to avoid
    stripping raw user input that happens to start with 'Goal:'.

    Returns:
        str: advice with preprocessor header lines removed
    """
    if not user_advice:
        return user_advice

    # 'Stop Condition:' is ALWAYS safe to strip — real users
    # write 'stop after X', not 'Stop Condition: stop after X'.
    #
    # v116.11: Updated regex to handle:
    #   - Markdown bold wrapping: '**Stop Condition**: ...'
    #   - Numbered list prefix:   '7. **Stop Condition**: ...'
    #   - Bullet list prefix:     '- Stop Condition: ...'
    # All variations are produced by the advice preprocessor LLM on
    # different runs; the original regex (anchored at line-start with
    # optional whitespace only) missed all but the plainest form.
    cleaned = re.sub(
        # Optional leading whitespace, list markers (1., 7., -, *),
        # markdown bold (**), then 'stop condition', then optional
        # close-bold and arbitrary punctuation before the colon.
        r'^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*'
        r'stop\s+condition'
        r'\**\s*:\s*[^\n]*$',
        '', user_advice,
        flags=re.IGNORECASE | re.MULTILINE,
    )

    # 'Goal:' / 'Primary Goal:' — only strip when we detect
    # preprocessor output (indicated by other signature headers).
    # This avoids stripping raw user input like 'Goal: R-free
    # below 0.25' when preprocess_advice=False.
    #
    # v116.11: Updated regex to handle markdown bold + numbered
    # list prefix, same as the Stop Condition strip above.
    _PREPROCESSOR_SIGNATURES = (
        r'input\s+files?\s+found',
        r'experiment\s+type',
        r'key\s+parameters?',
        r'special\s+instructions?',
    )
    is_preprocessed = any(
        re.search(r'(?i)^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*%s\**\s*:'
                  % sig, cleaned,
                  re.MULTILINE)
        for sig in _PREPROCESSOR_SIGNATURES
    )
    if is_preprocessed:
        cleaned = re.sub(
            r'^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*'
            r'(?:primary\s+)?goal'
            r'\**\s*:\s*[^\n]*$',
            '', cleaned,
            flags=re.IGNORECASE | re.MULTILINE,
        )

    # Collapse any resulting blank-line runs
    cleaned = re.sub(r'\n{3,}', '\n\n', cleaned)
    return cleaned

def _get_stop_directive_patterns():
    """Load program name mappings for stop-condition parsing from YAML.

    Returns list of (pattern, program_name) tuples, sorted by pattern
    length descending so more specific patterns match first (e.g.,
    "map_to_model" before "dock_in_map", "real_space_refine" before "refine").

    Falls back to hardcoded patterns with DeprecationWarning if YAML fails.
    """
    global _STOP_DIRECTIVE_PATTERNS_CACHE
    if _STOP_DIRECTIVE_PATTERNS_CACHE is not None:
        return _STOP_DIRECTIVE_PATTERNS_CACHE

    entries = []
    try:
        try:
            from libtbx.langchain.knowledge.yaml_loader import load_programs
        except ImportError:
            from knowledge.yaml_loader import load_programs
        programs = load_programs()
        for name, defn in programs.items():
            if not isinstance(defn, dict):
                continue
            pats = defn.get("stop_directive_patterns")
            if not pats or not isinstance(pats, list):
                continue
            for pat in pats:
                entries.append((pat, name))
    except Exception:
        pass

    if not entries:
        warnings.warn(
            "Could not load stop_directive_patterns from programs.yaml; "
            "using hardcoded fallback. Add stop_directive_patterns to "
            "programs.yaml to silence this warning.",
            DeprecationWarning,
            stacklevel=2,
        )
        entries = [
            (r'mtriage', 'phenix.mtriage'),
            (r'xtriage', 'phenix.xtriage'),
            (r'phaser', 'phenix.phaser'),
            (r'molecular\s+replacement', 'phenix.phaser'),
            (r'ligand\s*fit', 'phenix.ligandfit'),
            (r'fit.*ligand', 'phenix.ligandfit'),
            (r'refine', 'phenix.refine'),
            (r'autobuild', 'phenix.autobuild'),
            (r'map\s*to\s*model', 'phenix.map_to_model'),
            (r'maptomodel', 'phenix.map_to_model'),
            (r'build.*model.*(?:into|in)\s*(?:the\s*)?map', 'phenix.map_to_model'),
            (r'(?:automated?\s+)?model\s+building.*(?:into|in)\s*(?:the\s*)?map',
             'phenix.map_to_model'),
            (r'dock.*(?:in|into)\s*(?:the\s*)?map', 'phenix.dock_in_map'),
            (r'fit.*model.*(?:to|into)\s*(?:the\s*)?map', 'phenix.dock_in_map'),
            (r'map\s*symmetry', 'phenix.map_symmetry'),
            (r'(?:determin|find|detect).*symmetry', 'phenix.map_symmetry'),
            (r'symmetry.*(?:map|cryo)', 'phenix.map_symmetry'),
        ]

    # Sort by pattern length descending — longer (more specific) patterns
    # are checked first, handling cases like map_to_model vs dock_in_map.
    entries.sort(key=lambda x: len(x[0]), reverse=True)
    _STOP_DIRECTIVE_PATTERNS_CACHE = entries
    return entries


# =============================================================================
# EXTRACTION PROMPT
# =============================================================================

DIRECTIVE_EXTRACTION_PROMPT = """You are analyzing user instructions for a PHENIX crystallography automation system.

Extract structured directives from this user advice. Be precise and extract ONLY what is explicitly stated or clearly implied.

=== USER ADVICE ===
{user_advice}
=== END USER ADVICE ===

Output a JSON object with these sections. Include ONLY sections that have relevant content from the user advice:

1. "program_settings": Program-specific parameters
   - Use exact program names: "phenix.refine", "phenix.autosol", "phenix.autobuild", "phenix.phaser", "phenix.molprobity", "phenix.predict_and_build", "phenix.process_predicted_model", "phenix.real_space_refine", "phenix.dock_in_map", "phenix.map_to_model", "phenix.map_sharpening", "phenix.polder", "phenix.xtriage", "phenix.mtriage", "phenix.map_symmetry", "phenix.ligandfit", "phenix.resolve_cryo_em"
   - IMPORTANT: Do NOT use "phenix.resolve" — the correct name is "phenix.resolve_cryo_em"
   - Use "default" for settings that apply to all programs unless overridden
   - Common parameters and their types:
     * resolution: float (e.g., 2.5)
     * cycles: int (number of refinement macro-cycles)
     * anisotropic_adp: bool (anisotropic B-factors)
     * add_waters: bool (ordered solvent)
     * simulated_annealing: bool
     * atom_type: string (e.g., "Se", "S", "Zn")
     * wavelength: float (e.g., 0.9792)
     * sites: int (number of anomalous sites)
     * twin_law: string (e.g., "-h,-k,l")
     * riding_hydrogens: bool
     * unit_cell: string — space-separated "a b c alpha beta gamma" (e.g., "116.097 116.097 44.175 90 90 120")
     * space_group: string (e.g., "P 32 2 1", "P 1", "C 2 2 21")
     * copies: int (number of copies of the search model in the ASU, e.g. 4 — place under "default")
   - IMPORTANT: unit_cell, space_group, and copies apply to ALL programs — always place them under "default", not under a specific program
   - unit_cell FORMAT: always convert the user's value to a space-separated string of 6 numbers in order a b c alpha beta gamma
     * "(116.097, 116.097, 44.175, 90, 90, 120)" → "116.097 116.097 44.175 90 90 120"
     * "116 116 44 90 90 120" → "116 116 44 90 90 120"
     * Never use parentheses or commas in the extracted string value
   - For phenix.map_sharpening specifically:
     * resolution: float (required for model-based sharpening)
     * sharpening_method: string (e.g., "b-factor", "model_sharpening")
   - For phenix.polder specifically:
     * selection: string (atom selection, e.g., "chain A and resseq 88", "resname LIG")

2. "stop_conditions": When to stop the workflow (ONLY if the user gave an EXPLICIT stop command — see rules below)
   - "after_program": string - Stop after this program completes (e.g., "phenix.xtriage", "phenix.refine", "phenix.phaser", "phenix.ligandfit", "phenix.map_to_model", "phenix.dock_in_map", "phenix.map_sharpening", "phenix.polder")
   - "after_cycle": int - Stop after this cycle number
   - "max_refine_cycles": int - Maximum number of refinement cycles to run
   - "skip_validation": bool - If true, allow stopping without running molprobity
   - "r_free_target": float - Stop when R-free reaches this value
   - "map_cc_target": float - Stop when map correlation reaches this value
   - "stop_after_requested": bool, optional
     Set to true when the user's RAW INSTRUCTION explicitly says to stop
     after a specific program. Recognized phrasings:
       - "stop after X" / "X and stop" / "X, then stop"
       - "stop when ..." / "stop once ..." / "stop if ..." / "stop at ..."
       - "only run X" / "just run X" / "just do X"
       - "Stop Condition: <concrete value>" (a real condition, not "None")
       - "stop the workflow [immediately] after X" / "stop X after it completes"
       - "after X completes/finishes, stop" / "after X is done, stop"
       - "X is the last step" (X is the user-intended terminal program)
       - Any user phrasing that names a target program followed by a
         clear signal that no further programs should run after that one.
     Set false or omit when:
       - No such phrasing in the raw instruction
       - Raw says "don't stop" / "never stop" / "continue past X"
       - Stop Condition reads "None" / "not specified" / "N/A"
     The raw instruction is the authority. If raw says "X and stop" but
     preprocessed says "Stop Condition: None", set this true.
     Examples:
       - "stop the workflow immediately after refinement completes" →
         after_program=phenix.refine, stop_after_requested=true,
         skip_validation=true
       - "<program> is the last step" → after_program=<program>,
         stop_after_requested=true, skip_validation=true

3. "file_preferences": Specific files to use or avoid
   - "model": string - Preferred model file name
   - "sequence": string - Preferred sequence file name
   - "data_mtz": string - Preferred reflection data file (Fobs for refinement)
   - "map_coeffs_mtz": string - Preferred map coefficients file (for ligand fitting)
   - "exclude": list of strings - Files to avoid using

4. "workflow_preferences": High-level workflow choices
   - "skip_programs": list of strings - Programs to never run
   - "prefer_programs": list of strings - Programs to prefer when multiple options exist
   - "use_experimental_phasing": bool - Prefer SAD/MAD over molecular replacement
   - "use_molecular_replacement": bool - Prefer MR over experimental phasing
   - "use_mr_sad": bool - MR-SAD workflow: run phaser first to place model, then autosol with placed model
   - "model_is_placed": bool - The user's model is already positioned in the unit cell (skip MR)
   - "wants_validation_only": bool - User's PRIMARY goal is validation/analysis of an existing model (MolProbity), NOT refinement or MR. Only set when the user explicitly wants validation as the main task, NOT when validation is a final step after refinement.

**CRITICAL: use_experimental_phasing and use_mr_sad**
- ONLY set these if the user EXPLICITLY requests SAD, MAD, experimental phasing, anomalous phasing, or MR-SAD.
- Do NOT infer phasing method from the data files, wavelength values, or atom types.
- The system automatically detects anomalous signal via xtriage and adjusts the workflow accordingly.
- If the user just provides data + sequence/model without mentioning phasing method, leave these unset.
- Examples of when to set: "use SAD phasing", "run MR-SAD", "experimental phasing with Fe", "anomalous phasing"
- Examples of when NOT to set: user provides wavelength/atom_type but doesn't mention SAD/phasing method
- EXAMPLES that MUST set use_mr_sad=true:
  - "MR-SAD using intrinsic sulfur (S-SAD)"
  - "molecular replacement using X.pdb, then run MR-SAD phasing"
  - "molecular replacement followed by SAD/anomalous phasing"
  - "perform molecular replacement, then AutoSol with placed model"

**CRITICAL: wants_validation_only**
- Set wants_validation_only=true ONLY when the user's PRIMARY GOAL is validation or analysis of an existing model.
- Examples of when to set: "model validation", "comprehensive validation", "analysis only", "run MolProbity on this structure", "structure validation and correction"
- Examples of when NOT to set: "refine and then validate", "solve the structure", any tutorial that mentions validation as a FINAL STEP after refinement

5. "constraints": list of strings - Other instructions that don't fit above categories
   - Keep these as clear, actionable statements
   - Examples: "Do not add waters until R-free < 0.30", "Use TLS after initial refinement"

IMPORTANT GUIDELINES:
- Extract actual values, not vague descriptions
- "first refinement" or "one refinement" as the ONLY goal means after_program="phenix.refine" with max_refine_cycles=1
- "stop after N cycles" means after_cycle=N
- "anisotropic" or "anisotropic refinement" means anisotropic_adp=true
- If user says "stop" or "don't run validation", set skip_validation=true
- Program-specific settings override default settings
- If no directives can be extracted, return empty object: {{}}

**CRITICAL: Do NOT confuse wavelength with resolution**
- Wavelength (in Angstroms) is the X-ray beam wavelength, typically 0.5-2.5 Å. Extract as wavelength, NEVER as resolution.
- Resolution is the diffraction limit (d_min), typically 1.5-4.0 Å. Only extract as resolution if explicitly stated as resolution/d_min.
- If the text says "Resolution limit: Not mentioned" or similar, do NOT extract a resolution value from anywhere else in the text.

**CRITICAL: unit_cell and space_group — always use "default" scope**
- When the user specifies a unit cell or space group, put it under "default" so it applies to all programs:
  ```json
  {{"program_settings": {{"default": {{"unit_cell": "116.097 116.097 44.175 90 90 120"}}}}}}
  ```
- Convert any parenthesized comma-separated tuple to a plain space-separated string: strip `(`, `)`, and replace `,` with spaces.
- Example triggers: "use unit cell (a, b, c, α, β, γ)", "the specified unit cell is ...", "space group P 32 2 1"

**CRITICAL: max_refine_cycles vs after_program vs after_cycle**
- max_refine_cycles=N: Limits the NUMBER of refinement jobs to N. The workflow continues normally until refinement, then stops after N refinement jobs.
- after_program="X": FORCES program X to be run IMMEDIATELY, bypassing normal workflow. Only use when user wants to skip directly to a specific program.
- after_cycle=N: Stops after N AGENT CYCLES (each cycle = one program execution). ONLY use when user says "stop after N cycles" with an explicit number.
- "maximum of one refinement" or "at most one refinement" → ONLY set max_refine_cycles=1, do NOT set after_program
- "solve the structure with one refinement" → max_refine_cycles=1 (workflow proceeds normally: xtriage → model → refine)
- "stop after refinement" or "stop after one refinement" → max_refine_cycles=1, skip_validation=true
- "just run refinement" or "only refinement" → after_program="phenix.refine" (skip to refinement immediately)
- Do NOT use after_cycle for "stop after refinement" — that would stop after the first agent cycle (e.g., xtriage), not after refinement.

**CRITICAL: skip_validation RULE**
If the user specifies ANY explicit stop condition (like "stop after X" or "Stop Condition: ..."),
ALWAYS set skip_validation=true. This tells the system that the user knows what they want
and doesn't need automatic validation before stopping. Examples:
- "Stop after running mtriage" → skip_validation=true
- "Stop after density modification" → skip_validation=true
- "Stop Condition: Stop after generating the improved map" → skip_validation=true

**CRITICAL: "don't run X" vs after_program**
If the user says "don't run X" or "skip X", NEVER set after_program to that program.
The after_program must be a program the user WANTS to run, not one they want to skip.
Instead, set after_program to the LAST program the user DOES want to run.
Example: "Run real_space_refine, don't run mtriage" → after_program="phenix.real_space_refine", skip_programs=["phenix.mtriage"]
Example: "Do one refinement, skip validation" → max_refine_cycles=1, skip_validation=true
Example: "Stop after refinement, do not run mtriage" → max_refine_cycles=1, skip_programs=["phenix.mtriage"]
WRONG: "Stop after refinement, do not run mtriage" → after_program="phenix.mtriage" (NEVER do this!)

**CRITICAL: CRYO-EM vs X-RAY PROGRAM SELECTION**:
Different programs are used for cryo-EM and X-ray experiments. Choose correctly based on the experiment type:

Cryo-EM ONLY programs:
- phenix.mtriage: Analyze cryo-EM map quality (NOT phenix.xtriage)
- phenix.map_to_model: Build/rebuild atomic model into cryo-EM map (NOT phenix.autobuild)
- phenix.dock_in_map: Dock existing model into cryo-EM map
- phenix.resolve_cryo_em: Cryo-EM density modification
- phenix.map_sharpening: Sharpen cryo-EM map
- phenix.map_symmetry: Find symmetry in cryo-EM map
- phenix.real_space_refine: Refine model against cryo-EM map (NOT phenix.refine)
- phenix.validation_cryoem: Validate cryo-EM model

X-ray ONLY programs:
- phenix.xtriage: Analyze X-ray data quality (NOT phenix.mtriage)
- phenix.autobuild: Build/rebuild model for X-ray (NOT for cryo-EM)
- phenix.phaser: Molecular replacement
- phenix.autosol: Experimental phasing (SAD/MAD)
- phenix.refine: Reciprocal-space refinement (NOT for cryo-EM)
- phenix.polder: Polder omit maps (X-ray only)

If the advice mentions cryo-EM, maps (.mrc/.ccp4), or real-space methods:
- "model building" or "rebuild model" → phenix.map_to_model (NOT phenix.autobuild)
- "refinement" → phenix.real_space_refine (NOT phenix.refine)
- "data analysis" → phenix.mtriage (NOT phenix.xtriage)
- "density modification" → phenix.resolve_cryo_em (NOT phenix.autobuild_denmod)

**TUTORIAL/PROCEDURE DETECTION**:
If the user gives an EXPLICIT COMMAND to run a single specific program
(e.g., "run xtriage", "just run phaser", "only analyze data quality"),
add appropriate stop_conditions:
- "run xtriage", "check for twinning", "analyze data" → after_program="phenix.xtriage", skip_validation=true
- "run phaser", "test MR", "try molecular replacement" → after_program="phenix.phaser", skip_validation=true
- "run one refinement", "quick refinement test" → after_program="phenix.refine", max_refine_cycles=1, skip_validation=true
- "check data quality", "analyze reflection data" → after_program="phenix.xtriage", skip_validation=true
- "run mtriage", "analyze map quality" → after_program="phenix.mtriage", skip_validation=true
- "map symmetry", "determine symmetry", "find symmetry" → after_program="phenix.map_symmetry", skip_validation=true
- "map sharpening", "sharpen the map", "sharpen map", "automatic sharpening" → after_program="phenix.map_sharpening", skip_validation=true

**CRITICAL: Do NOT fabricate stop conditions from goal descriptions.**
Only populate stop_conditions if the user has EXPLICITLY stated when to stop
using imperative language like "stop after X", "only run X", "just do X".
Tutorial descriptions of goals are NOT stop conditions.  Do not infer a stop
point from the tutorial's stated purpose.  If no explicit stop instruction is
present, do NOT include a stop_conditions section at all.

Examples of goal descriptions that are NOT stop conditions (do NOT extract):
- "The goal is to solve the structure by molecular replacement" → NO stop condition
- "This tutorial demonstrates density modification" → NO stop condition
- "Purpose: rebuild the model using AutoBuild" → NO stop condition
- "We will run MR-SAD to determine the structure" → NO stop condition

Examples of real stop commands (DO extract):
- "Stop after running phaser" → after_program="phenix.phaser"
- "Only run xtriage, nothing else" → after_program="phenix.xtriage"
- "Just do one refinement" → max_refine_cycles=1

If you see "Stop Condition: None" or no stop instruction, do NOT add stop_conditions.

**CRITICAL: model_is_placed — detecting when the model is already positioned**
ONLY set model_is_placed=true when the user EXPLICITLY states their model is already positioned or placed in the unit cell / map. This is a HIGH-PRECISION flag: when in doubt, do NOT set it. The workflow will figure out placement automatically.

Set model_is_placed=true ONLY when the user uses language like:
- "the model is already placed", "the model is already positioned", "skip molecular replacement"
- "skip docking", "skip MR", "the structure is solved", "I already ran phaser/dock_in_map"
- "run refinement on this placed model", "refine this placed model"
- User explicitly confirms no MR/docking step is needed because they placed it themselves

Do NOT set model_is_placed=true in any of these cases:
- User says "solve the structure" (ambiguous — may need MR or docking)
- User says "refine this model" or "run refinement" (the model may still need docking first)
- User says "fit a ligand" or "validate model" (use after_program stop instead)
- User provides a PDB alongside half-maps without explicitly saying the model is placed
- User provides a homology model / starting model for structure determination
- Goal is ambiguous or generic ("run the workflow", "analyze these files")
- User provides a PDB + cryo-EM map combination (model almost certainly needs docking)

IMPORTANT: For cryo-EM workflows where a PDB is provided alongside maps, do NOT set
model_is_placed=true unless the user explicitly says the model has already been docked.
An unplaced PDB + cryo-EM map always requires phenix.dock_in_map before refinement.

**CRITICAL: MR-SAD workflow**
- For MR-SAD experiments, set use_mr_sad=true in workflow_preferences. Do NOT set after_program="phenix.autosol".
- MR-SAD requires phaser to run FIRST to place the model, THEN autosol runs with the placed model.
- Setting after_program="phenix.autosol" would skip phaser, which breaks MR-SAD.
- Only set after_program="phenix.autosol" if the user explicitly says to skip molecular replacement and run autosol standalone.
  (Note: phenix.map_sharpening is the dedicated map sharpening tool - use this for sharpening requests)
- For "density modification", "denmod", or "density modify":
  CHECK THE EXPERIMENT TYPE in the preprocessed advice FIRST.
  * If experiment type is cryo-EM, OR input mentions half-maps, .ccp4, .mrc,
    mtriage, or resolve_cryo_em:
    → after_program="phenix.resolve_cryo_em", skip_validation=true
    (phenix.resolve_cryo_em is the cryo-EM density modification program;
     do NOT confuse with phenix.map_sharpening which is for sharpening)
  * If experiment type is X-ray, OR input mentions .mtz, .hkl, xtriage,
    phaser, autosol, SAD, MAD, or "improve phases":
    → after_program="phenix.autobuild_denmod", skip_validation=true
    (X-ray density modification uses phenix.autobuild with maps_only=True)
  * If experiment type is not explicitly stated and inputs are ambiguous:
    default to phenix.resolve_cryo_em.  Cryo-EM is the more common case
    for half-map inputs; X-ray density modification requires explicit
    X-ray-only signals like .mtz files.
- "MR-SAD", "MR SAD", "MRSAD", "molecular replacement SAD" → Set use_mr_sad=true in workflow_preferences AND use_experimental_phasing=true. Do NOT set after_program="phenix.autosol" because phaser must run first to place the model. The workflow is: phaser → autosol with the phaser output as partpdb_file.
  If user says "stop after autosol" or similar, set after_program="phenix.autosol" AND skip_validation=true.
- "dock in map", "fit model to map" → after_program="phenix.dock_in_map", skip_validation=true
- "map to model", "MapToModel", "build model into map", "automated model building", "rebuild model", "model rebuilding" → after_program="phenix.map_to_model", skip_validation=true
- "polder", "polder map", "omit map", "evaluate ligand placement" → after_program="phenix.polder", skip_validation=true
  (Note: phenix.polder calculates polder omit maps to evaluate ligand/residue placement in density)
- "fit ligand", "ligandfit", "place ligand" → after_program="phenix.ligandfit", skip_validation=true
- If the stop condition mentions generating a specific output file, set skip_validation=true.
- "model validation", "comprehensive validation", "analysis only", "validate this structure", "run MolProbity", "structure validation and correction" → Set wants_validation_only=true in workflow_preferences. Do NOT set after_program. Do NOT set this when validation is mentioned as a final step after refinement.

**CRITICAL: WORKFLOW CONTINUATION INDICATORS**:
Do NOT set after_program stop conditions if the user indicates they want additional steps AFTER a program:
- Words like "then", "later", "afterwards", "and then", "continue", "next", "followed by" indicate MULTI-STEP workflows
- "fit the ligand later" → Do NOT stop early, workflow should continue to ligand fitting
- "run predict_and_build, then refine" → Do NOT stop after predict_and_build
- "build model and add ligand" → Do NOT stop after model building
- "refine then validate" → Do NOT stop after refinement
- If user mentions a downstream task (ligand fitting, validation, additional refinement), do NOT set early stop
- Put downstream tasks in "constraints" instead so the agent knows to do them

**CRITICAL: predict_and_build is ALWAYS INTERMEDIATE**:
When a user says "solve the structure using predict_and_build" or "use PredictAndBuild workflow",
they want the FULL pipeline: prediction → molecular replacement → building → refinement.
Do NOT set after_program="phenix.predict_and_build" for structure solution requests.
Only set after_program="phenix.predict_and_build" if the user EXPLICITLY says "just run the prediction"
or "only predict the model, nothing else".

**CRITICAL: LIGAND FITTING WORKFLOWS**:
When user mentions ligand fitting as part of the workflow:
- "refine, fit ligand, then refine again" → Do NOT set after_program, let workflow complete naturally
- "stop after second refinement" (when ligand fitting is mentioned) → The second refine is AFTER ligandfit, so do NOT set after_program="phenix.refine"
- "one refinement, ligandfit, final refinement" → Do NOT set after_program, this is a complete workflow
- If ligand fitting is mentioned AND refinement after it, the workflow should be: refine → ligandfit → refine → stop
- Put "Fit ligand after first refinement" or similar in constraints, NOT in stop_conditions

Examples of what NOT to do:
- User says "solve structure using PredictAndBuild" → Do NOT set after_program="phenix.predict_and_build" (this is a full workflow, not a single step!)
- User says "run predict_and_build and fit ligand later" → Do NOT set after_program="phenix.predict_and_build"
- User says "try MR then refine" → Do NOT set after_program="phenix.phaser"
- User says "refine, fit ligand, refine again" → Do NOT set after_program="phenix.refine" (this would stop before ligandfit!)
- User says "stop after the second refinement" (with ligand workflow) → Do NOT set after_program="phenix.refine"

Output ONLY valid JSON. No explanation, no markdown code blocks, just the JSON object."""


# v117 Step 1: dual-input prompt for raw-advice-authoritative extraction.
#
# Derived from DIRECTIVE_EXTRACTION_PROMPT by swapping the USER ADVICE
# block for a dual-input block (raw + processed) plus an AUTHORITY
# paragraph telling the LLM that raw is authoritative for intent
# (after_program, start_with_program, stop_after_requested).
#
# Used by extract_directives() when raw_advice is provided and differs
# from the processed advice — i.e., when the preprocessor's transformation
# may have lost or distorted the user's intent. Single-input behavior
# (backward compat) is preserved when raw_advice is None or matches
# user_advice.
_DUAL_INPUT_BLOCK = """=== USER'S RAW INSTRUCTION (authoritative for intent) ===
{raw_advice}
=== END RAW INSTRUCTION ===

=== PREPROCESSED ADVICE (for file mentions, parameters, experiment type) ===
{processed_advice}
=== END PREPROCESSED ADVICE ===

AUTHORITY:
The RAW INSTRUCTION is the source of truth for the user's intent.
This includes which program(s) the user wants to run (after_program,
start_with_program) and whether the user wants the workflow to stop
after a specific program (stop_after_requested).

When the raw and preprocessed advice disagree about intent — for
example, raw says "X and stop" but preprocessed says
"Stop Condition: None", or raw mentions one action but preprocessed
expands to multiple — defer to the raw instruction. The preprocessor
sometimes loses stop signals or hallucinates additional steps when
expanding short commands.

Use the PREPROCESSED ADVICE for fields the raw instruction is silent
on: specific filenames, numeric parameters, experiment type, and
crystallographic details."""

_SINGLE_INPUT_BLOCK = """=== USER ADVICE ===
{user_advice}
=== END USER ADVICE ==="""

DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW = (
    DIRECTIVE_EXTRACTION_PROMPT.replace(
        _SINGLE_INPUT_BLOCK,
        _DUAL_INPUT_BLOCK,
    ))


# =============================================================================
# EXTRACTION FUNCTION
# =============================================================================

def extract_directives(user_advice, provider="google", model=None, log_func=None,
                       use_rules_only=False, raw_advice=None):
    """
    Extract structured directives from user advice using LLM.

    Args:
        user_advice: Natural language user instructions (typically the
            preprocessed advice when called from the server).
        provider: LLM provider ("google", "openai", "anthropic")
        model: Specific model to use (optional, uses default for provider)
        log_func: Optional logging function
        use_rules_only: If True, skip LLM and use simple pattern extraction
        raw_advice: Optional raw user instruction (pre-preprocessing).
            When provided and different from user_advice, the LLM sees
            both: raw is authoritative for intent (after_program,
            start_with_program, stop_after_requested); processed is the
            source for file mentions, parameters, and experiment type.
            When None or equal to user_advice, falls back to the
            single-input prompt (backward compat).

    Returns:
        dict: Extracted directives, or empty dict if extraction fails

    Example:
        >>> directives = extract_directives(
        ...     "use resolution 3 in autosol but 2.5 elsewhere",
        ...     provider="google"
        ... )
        >>> directives["program_settings"]["phenix.autosol"]["resolution"]
        3.0
    """
    def log(msg):
        if log_func:
            log_func(msg)

    # Handle empty or trivial advice
    if not user_advice or not user_advice.strip():
        log("DIRECTIVES: No user advice provided")
        return {}

    # Strip LLM-fabricated "Stop Condition:" lines from preprocessed
    # advice before any extraction path sees them (Fix 2, v115).
    user_advice = _strip_preprocessor_stop_condition(user_advice)

    # Classify intent (v115) — same as simple path.
    intent_result = None
    try:
        try:
            from libtbx.langchain.agent.intent_classifier \
                import classify_intent
        except ImportError:
            from agent.intent_classifier \
                import classify_intent
        intent_result = classify_intent(user_advice)
        if intent_result:
            log("DIRECTIVES: Intent classified as '%s' "
                "(%s)" % (intent_result.get("intent"),
                           intent_result.get("reason", "")[:60]))
    except Exception:
        pass  # intent_classifier not available

    advice_lower = user_advice.lower().strip()
    if advice_lower in ("none", "none provided", "n/a", ""):
        log("DIRECTIVES: User advice is empty/none")
        return {}

    # Check if advice is too short to contain meaningful directives
    if len(advice_lower) < 10:
        log("DIRECTIVES: User advice too short for directive extraction")
        return {}

    # Skip LLM if use_rules_only is set
    if use_rules_only:
        log("DIRECTIVES: Skipping LLM (use_rules_only=True), using simple patterns")
        return extract_directives_simple(user_advice)

    # Build the prompt
    # v116.20: use .replace() instead of .format() because the prompt template
    # contains unescaped JSON examples like {"program_settings": {...}} whose
    # braces would be misinterpreted as format placeholders, causing
    # KeyError when extract_directives is called standalone (e.g. from tests
    # or debug scripts).  Production avoids this because directive extraction
    # routes through the ai_analysis server, but standalone callers trip on
    # it.  .replace() does the same single substitution without invoking
    # str.format()'s placeholder mini-language; behavior-preserving for the
    # production path.
    #
    # v117 Step 1: dual-input form when raw_advice differs from
    # user_advice — gives the extractor LLM the raw user instruction
    # alongside the preprocessed advice, with an AUTHORITY paragraph
    # telling it that raw is authoritative for intent fields.  Falls
    # back to single-input when raw_advice is None or matches
    # user_advice (backward compat).
    if (raw_advice and raw_advice.strip()
            and raw_advice != user_advice):
        prompt = (DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
                  .replace("{raw_advice}", raw_advice)
                  .replace("{processed_advice}", user_advice))
    else:
        prompt = DIRECTIVE_EXTRACTION_PROMPT.replace(
            "{user_advice}", user_advice)

    # Call LLM
    try:
        response = _call_llm(prompt, provider, model, log)
        if not response:
            log("DIRECTIVES: No response from LLM")
            # For ollama, fall back to simple extraction since local models
            # may struggle with complex JSON extraction tasks
            if provider == "ollama":
                log("DIRECTIVES: Falling back to simple pattern extraction for ollama")
                return extract_directives_simple(user_advice)
            return {}

        # Log response length for debugging
        log("DIRECTIVES: Got response from %s (%d chars)" % (provider, len(response)))

        # Parse JSON response
        directives = _parse_json_response(response, log)

        # Log what was parsed
        if not directives:
            log("DIRECTIVES: Parsed to empty dict - LLM may not have found actionable directives")
            # Show first part of response for debugging
            log("DIRECTIVES: Response preview: %s" % response[:300].replace('\n', ' '))
            # For ollama, try simple extraction as fallback since smaller models
            # may not follow JSON formatting instructions well
            if provider == "ollama":
                log("DIRECTIVES: Trying simple pattern extraction as ollama fallback")
                simple_directives = extract_directives_simple(user_advice)
                if simple_directives:
                    log("DIRECTIVES: Simple extraction found: %s" % list(simple_directives.keys()))
                    return simple_directives
        else:
            log("DIRECTIVES: Parsed sections: %s" % list(directives.keys()))

        # Validate and clean directives
        directives = validate_directives(directives, log)

        # v116.19a (Option A): narrow grounding guardrail — drop
        # LLM-fabricated after_program / prefer_programs in two cases:
        #   (1) the program name doesn't appear anywhere in the advice
        #       (pure fabrication, e.g. AF_7mjs's
        #       after_program=phenix.real_space_refine where neither
        #       "real_space_refine" nor any spelling variant appears in
        #       the README); or
        #   (2) the advice came through the preprocessor AND has no
        #       explicit stop intent AND no imperative marker near the
        #       program name (AF_7mjs's misextraction of
        #       phenix.predict_and_build would also fall here, since
        #       "PredictAndBuild" is in the advice but no "stop after"
        #       or "only" or "just" appears near it and the user said
        #       Stop Condition: None).
        # Bare user input like "run xtriage on my data" passes the
        # guardrail (program name is present, advice is not
        # preprocessed) so the directive extractor's documented
        # bare-imperative mappings continue to work.
        directives = _validate_after_program_grounded(
            directives, user_advice, log)

        # v118.9: correct LLM-emitted after_program when it's
        # canonical for the wrong experiment type (e.g., picked
        # autobuild_denmod for cryo-EM data).  See module-level
        # PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE for the mapping table
        # and `_apply_experiment_type_program_reprints` for placement
        # rationale (runs AFTER grounding by design).
        directives = _apply_experiment_type_program_reprints(
            directives, user_advice, raw_advice, log)

        # Regex fallback: ensure unit_cell and space_group are always captured.
        # The LLM sometimes returns an empty dict (or only stop_conditions) even
        # when the advice contains an explicit unit cell.  This pass is
        # deterministic and only fills in fields the LLM left empty.
        directives = _apply_crystal_symmetry_fallback(directives, user_advice, log)

        # v117.2: when the LLM set stop_after_requested=True but did NOT
        # set after_program, fall back to the regex resolver on the raw
        # advice to fill in after_program.  This closes a gap surfaced
        # by C1 c1_refine_raw_authority on openai: the LLM occasionally
        # signals user-stop intent without specifying the target program,
        # which leaves workflow_engine's gated wipe code unreachable
        # (the wipe path is gated on `if after_program:`).  The resolver
        # uses _ACTION_TABLE to map verbs in the raw advice (refine,
        # xtriage, density modify, ...) to programs.
        #
        # Guards:
        #  - Only fires when the LLM left after_program unset.
        #  - Only fires when raw advice contains stop intent (per
        #    _is_stop_after_requested), preventing the resolver from
        #    amplifying a possibly-hallucinated flag.
        #  - Uses raw_advice if available; falls back to user_advice
        #    for callers that don't pass raw_advice (single-input path).
        _sc_v172 = directives.get("stop_conditions") or {}
        if (_sc_v172.get("stop_after_requested") is True
                and not _sc_v172.get("after_program")):
            _v172_source = raw_advice if raw_advice else user_advice
            if _v172_source and _is_stop_after_requested(_v172_source):
                _v172_before = dict(_sc_v172)
                _resolve_after_program(directives, _v172_source.lower())
                _sc_v172_after = directives.get("stop_conditions") or {}
                if _sc_v172_after.get("after_program"):
                    log("DIRECTIVES: Filled in after_program=%s from raw "
                        "advice (LLM set stop_after_requested but omitted "
                        "after_program)" % _sc_v172_after["after_program"])

        # If validation stripped everything for ollama, try simple extraction
        if not directives and provider == "ollama":
            log("DIRECTIVES: Validation emptied directives, trying simple extraction")
            simple_directives = extract_directives_simple(user_advice)
            if simple_directives:
                return simple_directives

        # Store intent classification in directives (v115)
        # Low-confidence means the default fallback fired — not a
        # genuine detection, so don't pollute directives.
        if intent_result and intent_result.get("confidence") != "low":
            directives["intent"] = intent_result

        # Intent-driven stop override (v115) — same logic
        # as in extract_directives_simple.
        _intent = (
            directives.get("intent") or {}).get("intent")

        # Check if user EXPLICITLY said "stop after X"
        # (preserve explicit stops regardless of intent).
        # v116.x: use the centralized helper that handles
        # "Stop Condition: None" correctly.
        _has_explicit_stop = _is_stop_after_requested(user_advice)

        if _intent == "task":
            _task_prog = (
                directives.get("intent", {})
                .get("task_program"))
            if _task_prog:
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                directives["stop_conditions"][
                    "after_program"] = _task_prog
                directives["stop_conditions"][
                    "skip_validation"] = True
                # v116.x: intent="task" is a user-explicit single-
                # program request.  Flag the stop_after_requested so
                # the workflow engine treats it correctly.
                directives["stop_conditions"][
                    "stop_after_requested"] = True
        elif _intent == "solve":
            sc = directives.get("stop_conditions", {})
            if ("after_program" in sc
                    and not _has_explicit_stop):
                log("DIRECTIVES: Intent=solve, removing "
                    "after_program=%s"
                    % sc["after_program"])
                del sc["after_program"]
                sc.pop("skip_validation", None)
                sc.pop("stop_after_requested", None)
                if not sc:
                    directives.pop(
                        "stop_conditions", None)
        elif _intent == "solve_constrained":
            sc = directives.get("stop_conditions", {})
            if ("after_program" in sc
                    and not _has_explicit_stop):
                log("DIRECTIVES: Intent=solve_constrained"
                    ", removing after_program=%s"
                    % sc["after_program"])
                del sc["after_program"]
                sc.pop("skip_validation", None)
                sc.pop("stop_after_requested", None)
                if not sc:
                    directives.pop(
                        "stop_conditions", None)
            _method = (
                directives.get("intent", {})
                .get("method_constraint"))
            if _method:
                if "workflow_preferences" not in directives:
                    directives["workflow_preferences"] = {}
                directives["workflow_preferences"][
                    "method_constraint"] = _method
        # tutorial: keep whatever the LLM set

        # Simplify intent to just the string for
        # downstream consumers
        if isinstance(directives.get("intent"), dict):
            directives["intent"] = directives[
                "intent"].get("intent")

        log("DIRECTIVES: Extracted %d directive sections"
            % len(directives))

        # v115.09: Post-LLM overlay — apply deterministic workflow
        # intent patterns.  The LLM often ignores wants_validation_only
        # and use_mr_sad even when the advice clearly requests them.
        # This ensures the critical routing flags are always set.
        _apply_workflow_intent_fallback(
            directives, user_advice.lower())

        return directives

    except Exception as e:
        log("DIRECTIVES: Extraction failed - %s" % str(e))
        # v118.E (Q3 per Gemini): developer-visibility log line for
        # the silent-fallback case.  Without this, a developer
        # debugging a strange run might assume the LLM extracted
        # the parameters when it was actually the regex backstop
        # saving the cycle.
        log("CRITICAL_FALLBACK: LLM extraction failed; falling "
            "back to rule-based directive extractor.")
        # v118.E (Q1 per Gemini): operator-visibility stderr marker.
        # log_func may be filtered or unwired in some deployments;
        # stderr is always visible.  Wrapped in try/except so the
        # diagnostic can never break the existing fallback path
        # (same safety pattern as C-prime's [DIRECTIVE_LAYERS]
        # diagnostic emission).
        #
        # v119.H13: route through _classify_provider_error so the
        # operator sees one of:
        #   [DIRECTIVE_EXTRACTION_MODEL_RETIRED]      (action: update DEFAULT_MODELS)
        #   [DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE]  (action: pull/rename — Tom's case)
        #   [DIRECTIVE_EXTRACTION_AUTH_FAILED]        (action: check API key)
        # instead of an opaque [DIRECTIVE_EXTRACTION_FAILED] for the
        # primary api_client path.  The fallback _call_llm_fallback
        # path also classifies (at the per-provider exception block);
        # H13 makes the primary path match.
        try:
            # Resolve the model name we attempted, matching the same
            # precedence the api_client path uses, so the marker's
            # hint mentions the actual model name (not a placeholder).
            try:
                attempted_model = (
                    model
                    or resolve_model_for_provider(provider))
            except Exception:
                attempted_model = model or "<unknown>"
            tag, hint = _classify_provider_error(e, attempted_model)
            import sys
            if tag == 'DIRECTIVE_EXTRACTION_MODEL_RETIRED':
                _emit_retired_model_marker(
                    provider, attempted_model, log)
            elif tag == 'DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE':
                _emit_unavailable_model_marker(
                    provider, attempted_model, hint, log)
            # Always also emit the legacy generic FAILED marker so
            # existing log-grep pipelines (which Tom's machine relies
            # on) still see the v118.E format.  The richer marker
            # appears IN ADDITION to FAILED, not in place of it.
            sys.stderr.write(
                "[DIRECTIVE_EXTRACTION_FAILED] LLM extraction "
                "call did not complete (provider=%s): %s. "
                "Falling back to rules-only resolver.\n"
                % (provider, str(e)))
            sys.stderr.flush()
        except Exception:
            pass
        return {}


def _call_llm(prompt, provider, model, log):
    """
    Call the LLM to extract directives.

    Uses the same LLM infrastructure as the main agent with rate limit handling.
    """
    try:
        # Try to use the existing API client infrastructure
        from libtbx.langchain.agent.api_client import call_llm_simple

        # Import rate limit handler - try multiple paths
        handler = None
        try:
            from libtbx.langchain.agent.rate_limit_handler import (
                get_google_handler, get_openai_handler, get_anthropic_handler
            )

            # Select handler based on provider
            if provider == "google":
                handler = get_google_handler()
            elif provider == "openai":
                handler = get_openai_handler()
            elif provider == "anthropic":
                handler = get_anthropic_handler()
        except ImportError:
            try:
                from agent.rate_limit_handler import (
                    get_google_handler, get_openai_handler, get_anthropic_handler
                )
                if provider == "google":
                    handler = get_google_handler()
                elif provider == "openai":
                    handler = get_openai_handler()
                elif provider == "anthropic":
                    handler = get_anthropic_handler()
            except ImportError:
                pass  # No retry logic available

        def log_wrapper(msg):
            log("DIRECTIVES: %s" % msg)

        def make_call():
            return call_llm_simple(
                prompt=prompt,
                provider=provider,
                model=model,
                temperature=0.1,  # Low temperature for consistent extraction
                max_tokens=2000
            )

        if handler:
            return handler.call_with_retry(make_call, log_wrapper)
        else:
            return make_call()

    except ImportError:
        # Fallback: try direct provider calls
        log("DIRECTIVES: Using fallback LLM call")
        return _call_llm_fallback(prompt, provider, model, log)


def _looks_like_retired_model_error(err_str, model_name):
    """Best-effort detection of 'model no longer exists' API responses.

    v119.H1: distinguishes retired-model 404s from generic API
    failures so operators see exactly which file to edit
    (core/llm.py) when a provider drops a model.

    Matches on either:
      - a strong retirement phrase ("no longer available",
        "has been deprecated", "is deprecated", "model not found"),
        which is unambiguous regardless of HTTP status, OR
      - a 404/NOT_FOUND signature combined with the model name or
        a weaker "is not found" phrase.

    False-negatives are fine -- the generic [DIRECTIVE_EXTRACTION_FAILED]
    still fires.  False-positives would mis-route a transient error
    as a retired model -- noisy but recoverable.

    Args:
      err_str: the exception message (str) from a failed API call
      model_name: the model name passed to the failing call

    Returns:
      bool
    """
    if not err_str:
        return False
    lower = err_str.lower()

    # Strong retirement phrases match regardless of HTTP status.
    # Real APIs report retirement variously: 404 with body, 400
    # with message, 200 OK with error in JSON.
    strong_phrases = (
        "no longer available",
        "has been deprecated",
        "is deprecated",
        "model not found",
    )
    for phrase in strong_phrases:
        if phrase in lower:
            return True

    # 404/NOT_FOUND class -- requires extra evidence to avoid false
    # positives on URL typos and similar generic 404s.
    has_404 = (
        ("404" in lower)
        or ("not_found" in lower)
        or ("not found" in lower))
    if has_404:
        if model_name and model_name.lower() in lower:
            return True
        if "is not found" in lower:
            return True
    return False


def _emit_retired_model_marker(provider, model_name, log):
    """Emit a [DIRECTIVE_EXTRACTION_MODEL_RETIRED] stderr marker.

    v119.H1.  Parallel to v118.E's [DIRECTIVE_EXTRACTION_FAILED] /
    [ADVICE_PREPROCESSING_FAILED] markers.  Tells the operator
    exactly which file to edit (core/llm.py) and which constants.

    Args:
      provider: the provider string ("google", "openai", ...)
      model_name: the retired model that was attempted
      log: the existing log function passed through _call_llm_fallback

    Never raises -- the stderr write is wrapped in try/except so a
    broken stderr can never break the fallback path.
    """
    try:
        import sys
        sys.stderr.write(
          "[DIRECTIVE_EXTRACTION_MODEL_RETIRED] "
          "provider=%s model=%r appears to be retired.  "
          "Update the relevant table in "
          "libtbx/langchain/core/llm.py "
          "(DECISION_MODEL_DEFAULTS / RAG_MODEL_DEFAULTS / "
          "RAG_EMBEDDING_DEFAULTS / EXPENSIVE_MODEL_DEFAULTS / "
          "CHEAP_MODEL_DEFAULTS) to a current model and add "
          "%r to RETIRED_MODELS.\n" % (
            provider, model_name, model_name))
        sys.stderr.flush()
    except Exception:
        # Stderr write must never break the fallback path.
        pass
    log("DIRECTIVES: model %s appears retired - see "
        "[DIRECTIVE_EXTRACTION_MODEL_RETIRED] in stderr"
        % model_name)


def _emit_unavailable_model_marker(provider, model_name, hint, log):
    """Emit a [DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE] stderr marker.

    v119.H13.  Parallel to _emit_retired_model_marker.  Distinct
    operational meaning: the model name is unknown to THIS server
    (e.g., not pulled on Ollama, or a typo in OLLAMA_LLM_MODEL).
    The action is to pull/rename, not to update DEFAULT_MODELS.

    Args:
      provider: the provider string
      model_name: the model name that was attempted
      hint: a short operator-action string from
        _classify_provider_error
      log: the existing log function

    Never raises -- the stderr write is wrapped in try/except.
    """
    try:
        import sys
        sys.stderr.write(
          "[DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE] "
          "provider=%s model=%r is unknown to this server.  "
          "%s\n" % (provider, model_name, hint or ""))
        sys.stderr.flush()
    except Exception:
        pass
    log("DIRECTIVES: model %s unavailable on this server - see "
        "[DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE] in stderr"
        % model_name)


def _call_llm_fallback(prompt, provider, model, log):
    """
    Fallback LLM call when main infrastructure not available.
    Includes retry logic for rate limit errors with exponential backoff and decay.

    Property (v119.H1): on any exception from the API call, this
    function emits at most one [DIRECTIVE_EXTRACTION_MODEL_RETIRED]
    marker (when the failure matches a retired-model signature) and
    returns None -- no model-family variant cycling.  Preserve this
    single-call property in future refactors.
    """
    # Import rate limit handler - try multiple paths
    get_google_handler = None
    get_openai_handler = None
    get_anthropic_handler = None

    try:
        from libtbx.langchain.agent.rate_limit_handler import (
            get_google_handler, get_openai_handler, get_anthropic_handler
        )
    except ImportError:
        try:
            # Try relative import for standalone testing
            from agent.rate_limit_handler import (
                get_google_handler, get_openai_handler, get_anthropic_handler
            )
        except ImportError:
            pass  # No retry logic available

    def log_wrapper(msg):
        log("DIRECTIVES: %s" % msg)

    if provider == "google":
        # v119.H1: resolve default once; both SDK paths below share it.
        if model is None:
            model_name = default_model_for_provider("google")
        else:
            model_name = model
        try:
            # Try new google.genai package first (recommended)
            try:
                from google import genai
                from google.genai import types

                # Get API key from environment
                api_key = os.environ.get("GOOGLE_API_KEY")
                if not api_key:
                    log("DIRECTIVES: GOOGLE_API_KEY not set in environment")
                    return None

                client = genai.Client(api_key=api_key)

                def make_call():
                    response = client.models.generate_content(
                        model=model_name,
                        contents=prompt,
                        config=types.GenerateContentConfig(
                            temperature=0.1,
                            max_output_tokens=2000
                        )
                    )
                    return response.text

                if get_google_handler:
                    handler = get_google_handler()
                    return handler.call_with_retry(make_call, log_wrapper)
                else:
                    return make_call()

            except ImportError:
                # Fall back to deprecated google.generativeai if new package not available
                import google.generativeai as genai_old

                # Configure API key from environment
                api_key = os.environ.get("GOOGLE_API_KEY")
                if not api_key:
                    log("DIRECTIVES: GOOGLE_API_KEY not set in environment")
                    return None
                genai_old.configure(api_key=api_key)

                gen_model = genai_old.GenerativeModel(model_name)

                def make_call():
                    response = gen_model.generate_content(
                        prompt,
                        generation_config=genai_old.types.GenerationConfig(
                            temperature=0.1,
                            max_output_tokens=2000
                        )
                    )
                    return response.text

                if get_google_handler:
                    handler = get_google_handler()
                    return handler.call_with_retry(make_call, log_wrapper)
                else:
                    return make_call()

        except Exception as e:
            err_str = str(e)
            # v119.H1: distinguish retired-model 404 from generic
            # API failures so operators see exactly which file to
            # edit.  Property: at most one marker per fallback call,
            # then return None -- no model-family variant cycling.
            if _looks_like_retired_model_error(err_str, model_name):
                _emit_retired_model_marker("google", model_name, log)
            log("DIRECTIVES: Google API call failed - %s" % err_str)
            return None

    elif provider == "openai":
        if model is None:
            model_name = default_model_for_provider("openai")
        else:
            model_name = model
        try:
            import openai

            client = openai.OpenAI()

            def make_call():
                response = client.chat.completions.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}],
                    temperature=0.1,
                    max_tokens=2000
                )
                return response.choices[0].message.content

            if get_openai_handler:
                handler = get_openai_handler()
                return handler.call_with_retry(make_call, log_wrapper)
            else:
                return make_call()

        except Exception as e:
            err_str = str(e)
            if _looks_like_retired_model_error(err_str, model_name):
                _emit_retired_model_marker("openai", model_name, log)
            log("DIRECTIVES: OpenAI API call failed - %s" % err_str)
            return None

    elif provider == "anthropic":
        if model is None:
            model_name = default_model_for_provider("anthropic")
        else:
            model_name = model
        try:
            import anthropic

            client = anthropic.Anthropic()

            def make_call():
                response = client.messages.create(
                    model=model_name,
                    max_tokens=2000,
                    messages=[{"role": "user", "content": prompt}]
                )
                return response.content[0].text

            if get_anthropic_handler:
                handler = get_anthropic_handler()
                return handler.call_with_retry(make_call, log_wrapper)
            else:
                return make_call()

        except Exception as e:
            err_str = str(e)
            if _looks_like_retired_model_error(err_str, model_name):
                _emit_retired_model_marker("anthropic", model_name, log)
            log("DIRECTIVES: Anthropic API call failed - %s" % err_str)
            return None

    elif provider == "ollama":
        # v119.H13 Item B: honor OLLAMA_LLM_MODEL env-var override.
        # Previously called default_model_for_provider("ollama")
        # unconditionally, which ignored the env-var contract that
        # core/llm.py:get_llm_and_embeddings has always honored.
        if model is None:
            model_name = resolve_model_for_provider("ollama")
        else:
            model_name = model
        try:
            # Ollama uses OpenAI-compatible API
            import openai

            # v119.H13 Item A: normalize the base URL.  Setting
            # OLLAMA_BASE_URL=http://localhost:11434 (no /v1 suffix)
            # is reasonable but previously broke OpenAI-SDK path
            # construction → "404 page not found".  The helper
            # idempotently appends /v1.  Default fallback is also
            # the bare host now; the helper supplies /v1 uniformly.
            base_url = normalize_ollama_openai_base_url(
                os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434"))

            client = openai.OpenAI(
                base_url=base_url,
                api_key="ollama"  # Ollama doesn't need a real API key
            )

            def make_call():
                response = client.chat.completions.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}],
                    temperature=0.1,
                    max_tokens=2000
                )
                return response.choices[0].message.content

            # No rate limit handler for ollama (local model)
            return make_call()

        except Exception as e:
            # v119.H13 Item C: route through unified classifier.
            # Three sub-classes: MODEL_RETIRED (update DEFAULT_MODELS),
            # MODEL_UNAVAILABLE (pull/rename), AUTH_FAILED (key issue),
            # FAILED (generic).  Backwards-compatible with v119.H1's
            # _emit_retired_model_marker path.
            err_str = str(e)
            tag, hint = _classify_provider_error(e, model_name)
            if tag == 'DIRECTIVE_EXTRACTION_MODEL_RETIRED':
                _emit_retired_model_marker("ollama", model_name, log)
            elif tag == 'DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE':
                _emit_unavailable_model_marker(
                    "ollama", model_name, hint, log)
            log("DIRECTIVES: Ollama API call failed - %s" % err_str)
            return None

    else:
        log("DIRECTIVES: Unknown provider %s" % provider)
        return None


def _parse_json_response(response, log):
    """
    Parse JSON from LLM response, handling common issues.
    """
    if not response:
        return {}

    # Clean up response
    text = response.strip()

    # Remove markdown code blocks if present
    if text.startswith("```"):
        # Find the end of the opening fence
        first_newline = text.find("\n")
        if first_newline > 0:
            text = text[first_newline + 1:]
        # Remove closing fence
        if text.endswith("```"):
            text = text[:-3]
        text = text.strip()

    # Try to find JSON object in response
    # Look for outermost braces
    start = text.find("{")
    end = text.rfind("}")

    if start >= 0 and end > start:
        text = text[start:end + 1]

    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        log("DIRECTIVES: JSON parse error - %s" % str(e))
        log("DIRECTIVES: Raw response: %s" % text[:200])
        return {}


# =============================================================================
# VALIDATION
# =============================================================================

# Valid program names
VALID_PROGRAMS = {
    "phenix.refine",
    "phenix.autosol",
    "phenix.autobuild",
    "phenix.phaser",
    "phenix.molprobity",
    "phenix.predict_and_build",
    "phenix.process_predicted_model",
    "phenix.real_space_refine",
    "phenix.dock_in_map",
    "phenix.map_to_model",
    "phenix.map_sharpening",
    "phenix.resolve_cryo_em",
    "phenix.polder",
    "phenix.ligandfit",
    "phenix.model_vs_data",
    "phenix.xtriage",
    "phenix.mtriage",
    "phenix.map_symmetry",
}

# Valid setting keys by type
VALID_SETTINGS = {
    "resolution": float,
    "cycles": int,
    "anisotropic_adp": bool,
    "add_waters": bool,
    "simulated_annealing": bool,
    "atom_type": str,
    "additional_atom_types": str,
    "wavelength": float,
    "sites": int,
    "twin_law": str,
    "riding_hydrogens": bool,
    "stop_after_predict": bool,
    "ncs": bool,
    "tls": bool,
    "sharpening_method": str,
    "selection": str,
    "unit_cell": str,    # space-separated "a b c alpha beta gamma"
    "space_group": str,  # e.g. "P 32 2 1"
    "copies": int,       # ASU copy count (e.g. "4 copies of the search model")
}


# v118.B2a: Translation map for known semantic near-misses
# (information-conservation per v118 master plan).
#
# When the preprocessor LLM emits a recognizable INTENT under a wrong
# PHIL namespace, we HEAL it to the bare strategy form rather than
# dropping it.  The bare form is what the planner/BUILD pipeline
# already knows how to translate to the correct namespace per program
# (see knowledge/prompts_hybrid.py:316, Strategy schema).
#
# The most common near-miss observed in production (AIAgent_211,
# testit fresh-run 2026-05-18) is r_free_flags.generate emitted under
# data_manager. or xray_data. namespace prefixes.
#
# Add new entries here as new near-miss patterns are observed in
# production logs.  Conservative by design: only translate when the
# semantic intent CERTAINLY matches the target bare name.
PHIL_NAMESPACE_TRANSLATIONS = {
    # All forms of "generate R-free flags" → bare strategy hint
    "data_manager.r_free_flags.generate":  "generate_rfree_flags",
    "xray_data.r_free_flags.generate":     "generate_rfree_flags",
    "r_free_flags.generate":               "generate_rfree_flags",
}


# v118.B2b: Bad PHIL namespace prefixes.
#
# Any program_settings key starting with one of these prefixes is a
# PHIL namespace identifier that does not belong in program-agnostic
# program_settings.  The cleaner drops them with a log message.
#
# IMPORTANT — program-specific prefixes are NOT included here:
#   autosol.   ligandfit.   phaser.   xtriage.   polder.   etc.
# Those programs are routinely orchestrated and their parameter
# namespaces (e.g. "autosol.atom_type=Se") must be allowed through.
BAD_PHIL_NAMESPACE_PREFIXES = (
    "data_manager.",      # iotbx data_manager scope
    "refinement.",        # internal phenix.refine scope
    "miller_array.",      # data_manager subfield
    "fmodel.",            # data_manager subfield
    "scaling.input.",     # xtriage internal
)


def _heal_namespaced_phil_keys(settings, prog_name, log_func):
    """Heal or drop namespaced PHIL keys in a program_settings dict.

    For each key in `settings`:

      Rule 1 (HEAL): if the key is in PHIL_NAMESPACE_TRANSLATIONS,
        replace it with the bare strategy form.  Preserves the
        user's semantic intent across a near-miss namespace error.

      Rule 2 (DROP): if the key starts with a known-bad PHIL prefix,
        drop it with a log message.  These prefixes identify scopes
        that don't belong in program-agnostic program_settings.

      Rule 3 (KEEP): everything else passes through unchanged.
        Preserves both known params and forward compatibility for
        future bare-name params we don't know about yet.  Also lets
        program-specific dotted keys (e.g. "autosol.atom_type")
        pass through for orchestration.

    Conflict policy: if a translated key would overwrite an existing
    bare key (i.e. both "data_manager.r_free_flags.generate" and
    "generate_rfree_flags" present), the EXISTING bare form is
    authoritative; the translation is logged but does not overwrite.

    Returns a NEW dict with healed/dropped keys applied.  Logs
    every translation and every drop with program name, key,
    value, and reason.
    """
    if not isinstance(settings, dict):
        return settings

    cleaned = {}
    for key, value in settings.items():
        # Rule 1: known semantic translation (heal)
        if key in PHIL_NAMESPACE_TRANSLATIONS:
            new_key = PHIL_NAMESPACE_TRANSLATIONS[key]
            if new_key in settings:
                # Bare form already present; preserve it, don't overwrite.
                log_func(
                    "DIRECTIVES: Discarding namespaced PHIL key in %s: "
                    "%s=%s (bare form %s already present)"
                    % (prog_name, key, value, new_key))
            elif new_key in cleaned:
                # Already translated something else to this bare key.
                log_func(
                    "DIRECTIVES: Skipping duplicate namespaced PHIL "
                    "key in %s: %s=%s (bare form %s already set "
                    "this pass)" % (prog_name, key, value, new_key))
            else:
                log_func(
                    "DIRECTIVES: Healed namespaced PHIL key in %s: "
                    "%s=%s -> %s=%s (bare form is correct strategy "
                    "hint)" % (prog_name, key, value, new_key, value))
                cleaned[new_key] = value
            continue

        # Rule 2: known-bad PHIL namespace prefix (drop)
        bad_prefix = None
        for prefix in BAD_PHIL_NAMESPACE_PREFIXES:
            if key.startswith(prefix):
                bad_prefix = prefix
                break
        if bad_prefix is not None:
            log_func(
                "DIRECTIVES: Dropping namespaced PHIL key from %s: "
                "%s=%s (prefix %r - emit bare name instead)"
                % (prog_name, key, value, bad_prefix))
            continue

        # Rule 3: keep everything else (existing behavior).
        cleaned[key] = value

    return cleaned

# Valid stop condition keys
VALID_STOP_CONDITIONS = {
    "after_program": str,
    "after_cycle": int,
    "max_refine_cycles": int,
    "skip_validation": bool,
    "r_free_target": float,
    "map_cc_target": float,
    # v116.x: gate for stop analysis in
    # workflow_engine._apply_directives.  True iff the directive
    # extractor detected an explicit user-stop intent.  Plan-injected
    # after_program (per-stage hint) does NOT set this flag.
    "stop_after_requested": bool,
    # v119.H14.1: start_with_program — set by _resolve_after_program
    # (called from _apply_workflow_intent_fallback) when the user's
    # advice contains multiple actions plus a stop intent.  Downstream
    # consumers in workflow_engine.py (line ~2288), ai_agent.py
    # (~2816, ~3004), and ai_analysis.py (~160) read this key.
    # Pre-H14.1 it was a latent bug: present everywhere in the codebase
    # but not in VALID_STOP_CONDITIONS, so validate_directives would
    # log "Unknown stop condition start_with_program" and drop it
    # whenever validate_directives was called on a dict that already
    # had it set.  This didn't break the LLM path (the resolver runs
    # AFTER validate, so the key is added post-validation), but it
    # broke when H14.1 added validate_directives at the end of
    # extract_directives_simple — Tom's run_39 rules-only-stop case
    # for 1029B-sad would have its start_with_program=phenix.phaser
    # silently dropped by validate.
    "start_with_program": str,
}


# v118.9: experiment-type-conditional program canonicalization.
#
# Some Phenix programs are canonical for a specific experiment type
# but have a counterpart for the other type that does roughly the
# equivalent thing.  The LLM directive extractor occasionally picks
# the wrong one — e.g. "density modify" on cryo-EM data gets mapped
# to `phenix.autobuild_denmod` (X-ray) instead of
# `phenix.resolve_cryo_em` (cryo-EM).
#
# This table declares known LLM-mistake mappings.  Keys are
# `(target_experiment_type, wrong_program)`; values are the right
# program for that experiment type.
#
# Detection logic lives in `_apply_experiment_type_program_reprints`.
# Adding new entries requires no code changes beyond adding a row
# here and a corresponding K-test in
# `tst_density_modify_experiment_type.py`.
#
# Be conservative when extending: include only pairs that are
# truly canonical-equivalents in different experiment types
# (same operation, different program).  Pairs that do
# semantically-different things (e.g. xtriage vs mtriage — data
# triage vs map triage) should NOT be in this table.
PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE = {
    # Cryo-EM data + LLM picked X-ray density-modification program.
    # Observed in production: runs 237/239 with raw advice
    # "density modify and stop" on half-map cryo-EM data.
    ("cryoem", "phenix.autobuild_denmod"): "phenix.resolve_cryo_em",

    # Mirror case: X-ray data + LLM picked cryo-EM density-modification
    # program.  Not yet observed in production but the validator handles
    # both directions symmetrically.
    ("xray", "phenix.resolve_cryo_em"): "phenix.autobuild_denmod",
}


def _normalize_unit_cell(value):
    """
    Normalise a unit-cell value to a plain space-separated string.

    The LLM may return the value in several formats:
        "(116.097, 116.097, 44.175, 90, 90, 120)"   ← parenthesised tuple
        "116.097, 116.097, 44.175, 90, 90, 120"     ← comma-separated
        "116.097 116.097 44.175 90 90 120"           ← already correct

    All are normalised to "116.097 116.097 44.175 90 90 120".

    Returns the normalised string, or None if the value doesn't parse as
    six numbers.
    """
    if not isinstance(value, str):
        return None
    # Strip outer brackets/parens
    cleaned = value.strip().strip("()[]")
    # Replace commas (with optional surrounding whitespace) with a single space
    cleaned = re.sub(r'\s*,\s*', ' ', cleaned).strip()
    # Verify we have exactly 6 numeric tokens
    tokens = cleaned.split()
    if len(tokens) != 6:
        return None
    try:
        floats = [float(t) for t in tokens]
    except ValueError:
        return None
    # Reconstruct with minimal precision (strip trailing zeros)
    return ' '.join(('%g' % f) for f in floats)


# v116.15: Sentinel values that the advice preprocessor LLM emits when
# the user did NOT specify a value.  The preprocessor uses a standardized
# template ("Space group: None", "Unit cell: Not specified", etc.) and
# these placeholders must never be stored as real directives — they
# pollute every session with meaningless values and confuse the planner.
#
# Used by BOTH `_apply_crystal_symmetry_fallback` (regex path) AND
# `validate_directives` (LLM-extracted path) so the rules cannot drift.
_SYMMETRY_SENTINELS = frozenset([
    "none", "not mentioned", "not specified", "not given", "not provided",
    "not applicable", "unknown", "n/a", "na", "null", "tbd",
    "to be determined", "see above", "auto", "automatic", "default",
    "identification",  # historical: seen in some preprocessor outputs
    "false", "true",   # the preprocessor sometimes outputs literal booleans
    # v119.H14: extended sentinel set for "Not explicitly ..." patterns.
    # The run_39_openai batch (and Tom's xtriage tutorial verification
    # for H13.1) showed qwen2.5:72b and other LLMs emitting these as
    # space_group values when the README didn't specify one.  Includes
    # truncated forms (LLM output sometimes hits length caps mid-phrase,
    # e.g., "Not explicitly mentio").
    "not explicitly mentioned", "not explicitly mentio",
    "not explicitly stated", "not explicitly state",
    "not explicitly given", "not explicitly giv",
    "not explicitly specified", "not explicitly specif",
    "not explicitly defined", "not explicitly defin",
    "not explicitly listed", "not explicitly list",
    "not explicitly noted", "not explicitly not",
])


def _is_symmetry_sentinel(value):
    """Return True when `value` is one of the standardized placeholder
    strings that mean "the user did not specify a value" rather than a
    real crystallographic value.

    Compares case-insensitively against the trimmed value.  Returns
    False for None and non-strings (caller is expected to have
    a string already).
    """
    if not isinstance(value, str):
        return False
    return value.strip().lower() in _SYMMETRY_SENTINELS


# v119.H14: Hermann-Mauguin space-group pattern.
#
# Accepts (case-insensitively): any short token starting with a lattice
# letter [PFICRHAB] followed by characters from the Hermann-Mauguin
# alphabet (digits, m, c, n, d, a, b, e, h, r, spaces, slashes, dashes,
# underscores, colons), with optional trailing "(No. N)" annotation.
#
# All 230 official space groups (International Tables Vol A) match.
#
# Alternative cell/origin settings (via colon-prefixed suffix) also
# match: R3:H vs R3:R (rhombohedral hexagonal-axes vs rhombohedral-axes
# settings — note the 'h' and 'r' in the alphabet), P4/n:1 vs P4/n:2
# (origin choices 1 and 2 for centrosymmetric groups), P21/c:b
# (unique-axis cell choice for monoclinic).  Per Gemini H14 review:
# these appear in PDB/mmCIF metadata and cctbx tool output; rejecting
# them as malformed would drop legitimate user input.
#
# The 'e' letter is required for the 2002-ITA renamings of the
# centered orthorhombic groups Aem2, Aea2, Cmce, Cmme, Ccce.
# The 'h' and 'r' letters are required for rhombohedral axis
# specifiers (R3:H, R3:R) — they don't appear in any of the 230
# standard symbols themselves.
#
# Examples rejected: "Not explicitly mentioned", "Solve the
# structure", "garbage value 42", "" — multi-word prose phrases
# that the pre-H14 negative checks let through.
#
# Acknowledged limitation: the regex is permissive within the
# space-group alphabet — short English words that happen to start
# with a lattice letter and use only alphabet characters would
# technically match: Panda, Fame, Bad, Cab, Can, Bed, Acme, etc.
# These are not realistic LLM outputs for space_group (LLMs emit
# either a valid HM symbol or an obviously-prose sentinel), and
# the pre-H14 negative checks already let them through.  The H14
# regex closes the LARGE class of prose leaks (multi-word, contains
# characters outside the alphabet) without trying to be a full
# space-group validator — that would require knowing the 230
# official symbols, which is beyond the scope of a directive
# sanity check.  See test_hm_form_known_limitation_short_words.
_HM_FORM_RE = re.compile(
    r'^[PFICRHAB]'
    r'[0-9mcndabehr\s/_\-:]{0,24}'
    r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$',
    re.IGNORECASE
)


def _looks_like_space_group(value):
    """Return True if `value` resembles a Hermann-Mauguin space-group
    symbol.  Used by validate_directives() to catch LLM-emitted prose
    (e.g., "Solve the structure") that survives the sentinel check.

    All 230 official Hermann-Mauguin space-group symbols return True.
    Multi-word English phrases and tokens containing characters outside
    the HM alphabet return False.

    Pre-H14 the validator relied on negative checks (sentinel set +
    starts-with-letter + length<=25), which let multi-word phrases
    through if they started with a letter.  This positive shape check
    closes that gap while preserving all real space groups.
    """
    if not isinstance(value, str):
        return False
    return bool(_HM_FORM_RE.match(value.strip()))


# =============================================================================
# v116.19a — Narrow Grounding Guardrail (Option A) for LLM-extracted
#             after_program / prefer_programs
# =============================================================================
#
# The directive-extraction LLM occasionally fabricates a stop_conditions.
# after_program or workflow_preferences.prefer_programs entry from goal
# text alone — e.g. extracting after_program=phenix.real_space_refine
# from "rebuild the loop and refine the model" in a multi-stage tutorial
# whose Stop Condition: field is literally "None".
#
# Observed rate is low (~7% across mixed evidence) but produces
# spurious directives that pollute every subsequent prompt and bias
# the planner.  v116.17 ensures the misextraction is non-fatal at the
# validate-step chokepoint; v116.19a prevents the bad directive from
# being emitted in the first place by validating LLM output against
# the source text.
#
# OPTION A — Narrow guardrail.  Drops a directive only when EITHER:
#   (1) The program name is COMPLETELY absent from user_advice in any
#       canonical spelling variant.  This catches pure fabrications
#       like AF_7mjs's after_program=phenix.real_space_refine (the
#       string "real_space_refine" never appears in the advice).
#       Rule (1) fires regardless of whether the advice came from a
#       preprocessor.
#
#   (2) BOTH of:
#       (a) The advice came through the preprocessor (we detect this
#           by signature headers like "Input Files Found",
#           "Experiment Type", etc.) AND lacks any explicit "stop after"
#           imperative anywhere, indicating the user did not request a
#           stop; AND
#       (b) No imperative marker appears within 300 chars of the
#           program-name occurrence in user_advice.
#
# Rule (2) catches the AF_7mjs class where the program name IS in the
# preprocessed advice ("PredictAndBuild" is part of the Primary Goal)
# but no stop was requested.  It also catches LLM hallucinations where
# a different program name was confabulated.
#
# Both rules deliberately DO NOTHING for raw (non-preprocessed) user
# advice that doesn't fail rule (1).  This preserves the directive
# extractor prompt's documented mappings for bare imperatives like
# "run xtriage" → after_program=phenix.xtriage.  Those phrases never
# arrive in preprocessed form, so rule (2)(a) fails and rule (1) is
# satisfied (program name IS in advice).

# Imperative phrases that indicate the user is scoping or stopping
# the workflow at a specific program.  Case-insensitive substring
# match on a window around the program-name occurrence.
_IMPERATIVE_STOP_MARKERS = (
    "stop after",
    "stop_after",
    "only run",
    "just run",
    "just do",
    "after running",
    "nothing else",
    "skip the rest",
    "skip everything",
    "and then stop",
    "and stop",
    "stop when",
    "stop_condition: stop",  # advice preprocessor sometimes uses this form
    "stop condition: stop",
    # v117.3: phrasings surfaced by the explicit_stop_after_phaser
    # openai failure (0/5 on v117.2).  All are window-bounded substring
    # matches near the program-name occurrence in preprocessed advice.
    # AF_7mjs verified clean against all five (none of these phrases
    # appear in AF_7mjs's stripped advice, so K2 grounding test
    # remains green).  See v117.3 plan §4 for the false-positive
    # analysis.
    "stop the workflow",
    "immediately after",
    "after it completes",
    "after it finishes",
    "is the last step",
)


# Signature headers produced by the advice preprocessor LLM.  When any
# of these appears as a top-of-line header in the advice text, we
# conclude the advice came through the preprocessor (and thus may
# contain a fabricated Stop Condition: None line that was stripped
# before extraction).  Matches the _PREPROC_SIGS list in
# _resolve_after_program (which sees the post-strip advice).
_PREPROCESSOR_SIGNATURE_PATTERNS = (
    r'input\s+files?\s+found',
    r'experiment\s+type',
    r'key\s+parameters?',
    r'special\s+instructions?',
)


def _advice_came_from_preprocessor(user_advice):
    """Return True if user_advice has the structural signature of
    advice-preprocessor output.

    Used by the v116.19a guardrail as gate (2)(a): the narrow rule
    only fires for preprocessed advice, leaving raw user input
    (e.g. ``"run xtriage on my data"``) untouched.
    """
    if not isinstance(user_advice, str) or not user_advice:
        return False
    advice_lower = user_advice.lower()
    for sig in _PREPROCESSOR_SIGNATURE_PATTERNS:
        if re.search(
                r'(?i)^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*%s\**\s*:'
                % sig, advice_lower, re.MULTILINE):
            return True
    return False


def _advice_has_explicit_stop_intent(user_advice):
    """Return True if user_advice anywhere contains an imperative stop
    intent.

    This is the global check that distinguishes "user wanted a stop"
    (in which case we trust the LLM's after_program even on preprocessed
    advice) from "user wanted no stop" (in which case any after_program
    is likely a fabrication from goal verbs).

    The check looks for ``\\bstop after\\b`` and a small set of
    related global imperatives.  More forgiving than the per-program
    window check because we don't need to tie the imperative to a
    specific program name — we only need to know if the user ever
    said "stop".
    """
    if not isinstance(user_advice, str) or not user_advice:
        return False
    advice_lower = user_advice.lower()
    # Use word boundaries to avoid substring matches like
    # "after running into trouble".
    patterns = (
        r'\bstop\s+after\b',
        r'\bstop_after\b',
        r'\bonly\s+run\b',
        r'\bjust\s+run\b',
        r'\bjust\s+do\b',
        r'\bnothing\s+else\b',
        r'\band\s+stop\b',
        r'\bstop\s+when\b',
    )
    for pattern in patterns:
        if re.search(pattern, advice_lower):
            return True
    return False


def _program_name_variants(program):
    """Return all spelling variants of `program` to search for in advice.

    For ``"phenix.real_space_refine"`` returns 5 variants:
        "phenix.real_space_refine"   — canonical
        "real_space_refine"           — bare suffix
        "real space refine"           — underscores → spaces
        "real-space-refine"           — underscores → hyphens
        "RealSpaceRefine"             — CamelCase

    For ``"phenix.predict_and_build"`` returns:
        "phenix.predict_and_build", "predict_and_build",
        "predict and build", "predict-and-build", "PredictAndBuild"

    Returns an empty tuple for None / empty / non-string input.
    """
    if not isinstance(program, str) or not program:
        return ()
    suffix = program
    if suffix.lower().startswith("phenix."):
        suffix = suffix[len("phenix."):]
    variants = [program]
    if suffix != program:
        variants.append(suffix)
    # underscored → spaces
    spaced = suffix.replace("_", " ")
    if spaced != suffix and spaced not in variants:
        variants.append(spaced)
    # underscored → hyphens
    hyphened = suffix.replace("_", "-")
    if hyphened != suffix and hyphened not in variants:
        variants.append(hyphened)
    # underscored → camelcased (e.g. predict_and_build → PredictAndBuild)
    if "_" in suffix:
        camel = "".join(part.capitalize() for part in suffix.split("_"))
        if camel and camel not in variants:
            variants.append(camel)
    return tuple(variants)


def _find_variant_in_text(variants, text, log=None):
    """Find the first (variant, position) pair where `variant` appears
    in `text` on a word boundary.  Case-insensitive.  Returns
    (None, -1) if no variant matches.

    Word-boundary matching avoids matching "refine" inside the
    longer token "real_space_refine".
    """
    text_lower = text.lower()
    for variant in variants:
        # Build a word-boundary pattern.  re.escape handles dots and
        # other regex metacharacters in the variant text.  Word
        # boundaries are anchored at non-word transitions.
        pattern = r'(?<![A-Za-z0-9_])' + re.escape(variant.lower()) + \
                  r'(?![A-Za-z0-9_])'
        match = re.search(pattern, text_lower)
        if match:
            return (variant, match.start())
    return (None, -1)


def _imperative_marker_nearby(text, position, window=300):
    """Return True iff any phrase in `_IMPERATIVE_STOP_MARKERS`
    appears in `text` within `window` characters before or after
    `position`.  Case-insensitive substring match.

    Substring matching (not word-boundary) means ``"just run"`` would
    match ``"just running"``.  This is intentional — false positives
    here cause the guardrail to KEEP a directive that arguably should
    be dropped (a benign failure), while false negatives would DROP a
    legitimate user-requested stop (which we want to avoid).
    """
    if position < 0:
        return False
    lo = max(0, position - window)
    hi = position + window
    window_text = text[lo:hi].lower()
    for marker in _IMPERATIVE_STOP_MARKERS:
        if marker in window_text:
            return True
    return False


def _is_program_grounded(program, user_advice, log=None):
    """v116.19a Option A grounding check.

    Returns True (= keep the directive) UNLESS one of the failure modes
    fires:

      Failure 1 — Pure fabrication: the program name (any variant) is
                  not present in user_advice.  Drop the directive.

      Failure 2 — Preprocessed + no-stop-intent + no imperative near
                  program: the advice came from the preprocessor,
                  the user did NOT request a stop, and no imperative
                  appears near the program name.  Drop the directive.

    Otherwise return True.

    This is the Option A narrow variant of the original v116.19 helper.
    Compared with the original (broad) version:

      * Both versions drop on Failure 1 (pure fabrication).
      * The original ALSO dropped on "program name in advice but no
        imperative nearby", regardless of preprocessor.  Option A
        does NOT drop in that case for raw (non-preprocessed) advice,
        so the directive-extractor's documented bare-imperative
        mappings ("run xtriage" → after_program=phenix.xtriage)
        continue to work.
    """
    if not program or not user_advice:
        # Defensive: if either input is empty, can't make a judgment.
        # Return True (keep) — let other layers handle.
        return True
    variants = _program_name_variants(program)
    if not variants:
        return True
    variant, position = _find_variant_in_text(variants, user_advice, log)

    # Failure 1: pure fabrication — program name absent from advice.
    if variant is None:
        if log:
            log("DIRECTIVES: %s not grounded — program name absent "
                "from user advice (pure fabrication)" % program)
        return False

    # Failure 2: preprocessed + no global stop intent + no imperative
    # marker near the program.
    if (_advice_came_from_preprocessor(user_advice) and
            not _advice_has_explicit_stop_intent(user_advice) and
            not _imperative_marker_nearby(user_advice, position)):
        if log:
            log("DIRECTIVES: %s not grounded — preprocessed advice "
                "with no stop intent and no imperative marker within "
                "300 chars of %r" % (program, variant))
        return False

    return True


def _validate_after_program_grounded(directives, user_advice, log):
    """Drop stop_conditions.after_program and workflow_preferences.
    prefer_programs entries that the LLM fabricated from goal text.

    See `_is_program_grounded` for the grounding criteria (v116.19a
    Option A: drops on pure fabrication OR on preprocessed advice
    with no stop intent and no nearby imperative).

    Returns the (possibly modified) directives dict.  Logs each drop.
    """
    if not isinstance(directives, dict) or not user_advice:
        return directives

    # after_program
    stop_cond = directives.get("stop_conditions")
    if isinstance(stop_cond, dict) and stop_cond.get("after_program"):
        after_prog = stop_cond["after_program"]
        # v117 Step 1 interaction: if the LLM also set
        # stop_after_requested=True on the same stop_conditions block,
        # that flag *is* the grounding signal — the LLM made an
        # explicit user-stop assertion based on the raw advice and the
        # AUTHORITY paragraph in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW.
        # In that case, skip the literal-name grounding check that
        # would otherwise drop after_program as "fabrication".  The
        # AF_7mjs case is preserved: AF_7mjs's preprocessed advice has
        # no positive stop signal, so stop_after_requested would NOT
        # be set, and the grounding check fires as before.
        if stop_cond.get("stop_after_requested") is True:
            if log:
                log("DIRECTIVES: Keeping after_program=%s "
                    "(stop_after_requested=True confirms user stop "
                    "intent; grounding check skipped)" % after_prog)
        elif not _is_program_grounded(after_prog, user_advice, log):
            log("DIRECTIVES: Dropping fabricated "
                "stop_conditions.after_program=%s" % after_prog)
            del stop_cond["after_program"]
            # Clean up empty stop_conditions
            if not stop_cond:
                del directives["stop_conditions"]

    # prefer_programs
    wf_prefs = directives.get("workflow_preferences")
    if isinstance(wf_prefs, dict) and wf_prefs.get("prefer_programs"):
        prefer = wf_prefs["prefer_programs"]
        if isinstance(prefer, list):
            kept = []
            for prog in prefer:
                if _is_program_grounded(prog, user_advice, log):
                    kept.append(prog)
                else:
                    log("DIRECTIVES: Dropping fabricated "
                        "workflow_preferences.prefer_programs "
                        "entry %s" % prog)
            if kept:
                wf_prefs["prefer_programs"] = kept
            else:
                del wf_prefs["prefer_programs"]
                if not wf_prefs:
                    del directives["workflow_preferences"]

    return directives


# v118.9: helpers for experiment-type-conditional program canonicalization.
# See PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE module-level table for the
# mapping data.  See v118_9_PLAN_rev2.md for design rationale and
# the analysis of why this validator runs AFTER
# _validate_after_program_grounded (Gemini's Gap A critique reviewed
# in detail).

def _detect_experiment_type_signals(combined_advice):
    """Return one of 'cryoem', 'xray', or None based on signals in
    the given advice text.

    Uses the same regex as the rules-only resolver in
    extract_directives_simple (around line ~3700).  Factored out so
    the LLM-path validator and the rules-only fallback share one
    canonical detector.

    The cryo-EM signals are: cryo-em, half-map (any spacing/hyphen),
    .ccp4 / .mrc file extensions, mtriage, resolve_cryo_em, full-map.

    The X-ray signals are: x-ray, .mtz / .hkl file extensions,
    xtriage, phaser, autosol, sad, mad, molecular replacement,
    autobuild_denmod.

    Returns None when signals are ambiguous (both fire) or absent
    (neither fires).  The caller MUST handle None as "decline to act".
    """
    if not combined_advice:
        return None
    advice_lower = combined_advice.lower()
    is_cryoem = bool(re.search(
        r'cryo-?em|half.?map|\.mrc\b|\.ccp4\b|mtriage|'
        r'resolve_cryo_em|full.?map',
        advice_lower))
    is_xray = bool(re.search(
        r'x-?ray|\.mtz\b|\.hkl\b|xtriage|phaser|autosol|'
        r'sad|mad|molecular.?replacement|autobuild_denmod',
        advice_lower))
    if is_cryoem and not is_xray:
        return "cryoem"
    if is_xray and not is_cryoem:
        return "xray"
    return None  # ambiguous (both or neither)


def _apply_experiment_type_program_reprints(
        directives, user_advice, raw_advice, log):
    """Correct after_program when the LLM picked a program that's
    canonical for the wrong experiment type.

    Reads PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE (module-level table).
    For each (target_type, wrong_program) -> right_program mapping,
    if:
      - directives' after_program matches wrong_program AND
      - experiment-type detection identifies target_type unambiguously

    ...then replaces with right_program.  Logs every correction via
    the [DIRECTIVE_CORRECTION] marker to both log_func and stderr
    (mirrors Section E's diagnostic safety pattern), and records a
    `_corrected_from` sidecar in the stop_conditions dict.

    Placement note: this validator runs AFTER
    _validate_after_program_grounded.  Why:

      - In the observed-bug case (Tom's runs 237/239), the LLM
        emits stop_after_requested=True alongside the wrong
        after_program.  The grounding check bypasses entirely in
        that case (see _validate_after_program_grounded:1697-1712),
        leaving the wrong program intact for us to correct.

      - In the no-user-intent case (stop_after_requested=False), the
        grounding check fires.  Since the LLM-fabricated wrong
        program name typically doesn't appear in the user advice
        either (user said "density modify", not "autobuild_denmod"),
        grounding drops the after_program.  Our validator then has
        nothing to correct — which is the safer outcome (no
        after_program is better than a wrong one).

      - Reversing the order would create a different bug: our
        validator would correct autobuild_denmod -> resolve_cryo_em,
        then grounding would notice that "resolve_cryo_em" also
        isn't in user_advice (user said "density modify") and
        drop it as fabrication.

    No-op when:
      - directives has no stop_conditions or no after_program
      - after_program is not a key in the reprints table
      - experiment-type signals are ambiguous (None)

    Args:
        directives: post-LLM directives dict (may be modified in
            place)
        user_advice: the preprocessed advice (may contain
            Experiment Type, file extensions)
        raw_advice: the raw user advice (may be None for
            single-input path).  When present, combined with
            user_advice for experiment-type detection to maximize
            signal coverage.
        log: log_func

    Returns:
        the (possibly modified) directives dict
    """
    if not isinstance(directives, dict):
        return directives
    stop_cond = directives.get("stop_conditions") or {}
    after_prog = stop_cond.get("after_program")
    if not after_prog:
        return directives

    # Combine raw + preprocessed to maximize signal coverage.
    # raw_advice may contain file mentions ("half-maps", ".ccp4")
    # that user_advice has already paraphrased; preprocessed
    # advice contains the explicit "Experiment Type: cryo-EM"
    # field.  Both are useful.
    combined = " ".join(filter(None, [raw_advice or "",
                                        user_advice or ""]))
    if not combined.strip():
        return directives

    target_type = _detect_experiment_type_signals(combined)
    if target_type is None:
        return directives  # ambiguous: decline to act

    key = (target_type, after_prog)
    right_prog = PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE.get(key)
    if right_prog is None:
        return directives  # not a known incorrect choice for this type

    # Perform correction
    stop_cond["after_program"] = right_prog
    stop_cond["_corrected_from"] = {
        "from": after_prog,
        "to": right_prog,
        "reason": "experiment_type_mismatch",
        "experiment_type": target_type,
    }
    diag_msg = (
        "[DIRECTIVE_CORRECTION] Mapped after_program=%s to %s "
        "based on Experiment Type: %s"
        % (after_prog, right_prog, target_type))
    if log:
        log(diag_msg)
    # Also write to stderr.  Mirrors Section E's safety pattern:
    # try/except so the diagnostic can never break the correction.
    try:
        import sys
        sys.stderr.write(diag_msg + "\n")
        sys.stderr.flush()
    except Exception:
        pass
    return directives


def _apply_crystal_symmetry_fallback(directives, user_advice, log):
    """
    Regex-based fallback for unit_cell and space_group extraction.

    Called unconditionally after LLM extraction so that an explicit unit cell
    or space group in the user's advice is always captured, even when the LLM
    returns an otherwise empty or incomplete dict.  Only fills in fields that
    the LLM left empty — never overwrites a value the LLM already extracted.

    This is critical because the directive-extraction LLM sometimes returns {}
    for complex advice that mixes stop-conditions with crystal-symmetry data.

    Args:
        directives: Already-validated directive dict (may be {})
        user_advice: Raw advice text to scan
        log: Logging callable

    Returns:
        Updated directives dict (new dict; caller's dict is not mutated)
    """
    if not user_advice:
        return directives

    existing_default = (
        (directives.get("program_settings") or {}).get("default") or {}
    )
    needs_uc = "unit_cell" not in existing_default
    needs_sg = "space_group" not in existing_default
    if not needs_uc and not needs_sg:
        return directives   # both already set — nothing to do

    directives = dict(directives)   # shallow copy so we don't mutate caller's dict

    def _ensure_default():
        """Return (possibly freshly created) default settings dict."""
        if "program_settings" not in directives:
            directives["program_settings"] = {}
        if "default" not in directives["program_settings"]:
            directives["program_settings"]["default"] = {}
        return directives["program_settings"]["default"]

    # ------------------------------------------------------------------
    # unit_cell — match: "unit cell (116.097, 116.097, 44.175, 90, 90, 120)"
    #             or  "unit_cell = 116 116 44 90 90 120"
    #             or  "unit cell: 116.097 116.097 44.175 90 90 120"
    # ------------------------------------------------------------------
    if needs_uc:
        _NUM = r'[-+]?\d+(?:\.\d+)?'
        _SEP = r'[\s,]+'
        _uc_pat = (
            r'unit[_ ]cell'
            r'(?:\s*[=:(\s]\s*|\s+(?:is|of|=|:)?\s*[(]?)'
            r'(' + _NUM + _SEP + _NUM + _SEP + _NUM + _SEP
              + _NUM + _SEP + _NUM + _SEP + _NUM + r')'
        )
        uc_match = re.search(_uc_pat, user_advice, re.IGNORECASE)
        if uc_match:
            nums = re.findall(r'[-+]?\d+(?:\.\d+)?', uc_match.group(1))
            if len(nums) == 6:
                normalized = _normalize_unit_cell(" ".join(nums))
                if normalized:
                    _ensure_default()["unit_cell"] = normalized
                    log("DIRECTIVES: Regex fallback extracted "
                        "unit_cell=%s" % normalized)

    # ------------------------------------------------------------------
    # space_group — match: "space group P 32 2 1"  "space_group=P63"
    # ------------------------------------------------------------------
    if needs_sg:
        _sg_pat = (
            r'space[_ ]group'
            r'(?:\s*[=:]\s*|\s+(?:is|of|=|:)?\s*)'
            r'([A-Za-z][A-Za-z0-9 /_-]{1,20})'
        )
        sg_match = re.search(_sg_pat, user_advice, re.IGNORECASE)
        if sg_match:
            sg_str = sg_match.group(1).strip().rstrip('.')
            # v116.15: Reject sentinel values that mean "not specified".
            # The advice preprocessor LLM emits a standardized template
            # that always includes lines like "Space group: None".  The
            # raw regex above happily matches these and stored the literal
            # sentinel string as a directive, polluting every session
            # with a meaningless space_group="None" entry.  Real space-
            # group symbols (P 32 2 1, P63, C2, etc.) never collide with
            # these natural-language placeholders.
            if sg_str and not _is_symmetry_sentinel(sg_str):
                _ensure_default()["space_group"] = sg_str
                log("DIRECTIVES: Regex fallback extracted "
                    "space_group=%s" % sg_str)
            elif sg_str:
                log("DIRECTIVES: Regex fallback ignored space_group=%r "
                    "(sentinel value, not a real space group)" % sg_str)

    return directives


def _coerce_setting_value(value, expected_type):
    """Convert a directive setting value to its expected type,
    handling LLM-shape variability across providers.

    Google's directive-extraction LLM tends to wrap scalar values
    in lists (e.g. ``["S"]`` for ``additional_atom_types``).
    OpenAI emits the same content as a bare scalar (``"S"``).
    Naive ``str(["S"])`` produces ``"['S']"`` — the original
    v118.10 bug.  This helper normalizes both shapes.

    Behavior:

    * value is list/tuple, length 0:
      - expected_type is str → return ``""``
      - non-str expected_type → raise ValueError
    * value is list/tuple, length 1:
      Unpack the single element and recurse with expected_type.
      Handles ``[False]`` → ``False``, ``[5]`` → ``5``,
      ``["S"]`` → ``"S"`` all correctly (no Python truthiness
      trap from ``bool([False]) == True``).
    * value is list/tuple, length 2+:
      - expected_type is str → join elements with `" "` (PHIL
        multi-value string syntax), strip outer whitespace
      - non-str expected_type → raise ValueError (no sensible
        scalar coercion)
    * value is dict + expected_type is str:
      raise ValueError (no sensible coercion).
    * all other cases:
      apply ``expected_type(value)`` as before (backward-compat).

    Future-proofing convention: this helper coerces lists to
    scalar (str via space-join, others via single-element unpack)
    ONLY for non-list expected_types.  If a future setting
    legitimately needs to be a list, declare it with
    ``expected_type=list`` in VALID_SETTINGS — this function will
    then pass-through the list unchanged.

    Known limitation: deeply nested lists (e.g. ``[["S"]]``)
    unpack to ``["S"]``, then space-join produces ``"['S']"`` —
    back to the original bug shape.  This is accepted because
    LLMs emitting nested lists for flat PHIL parameters is
    genuinely malformed output.

    Quote-stripping NOT performed: LLMs occasionally emit
    elements with internal quotes like ``["'xtriage'"]``.
    Stripping outer quotes would damage legitimate PHIL content
    like ``["O5'"]`` (sugar atom with prime mark) or selection
    strings with internal quotes.  Whitespace stripping is done
    only on the joined multi-element result.

    Args:
        value: the raw LLM-emitted value
        expected_type: one of str, int, float, bool from
            VALID_SETTINGS or VALID_STOP_CONDITIONS

    Returns:
        the value coerced to expected_type

    Raises:
        ValueError or TypeError if coercion isn't sensible
    """
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            if expected_type is str:
                return ""
            raise ValueError(
                "empty list cannot coerce to %s"
                % expected_type.__name__)
        if len(value) == 1:
            # Single-element: unpack and recurse.  Type-preserving
            # for bool/int/float/str.
            return _coerce_setting_value(value[0], expected_type)
        # Multi-element: only str can absorb (space-join is the
        # canonical PHIL multi-value-string syntax)
        if expected_type is str:
            return " ".join(str(x) for x in value).strip()
        raise ValueError(
            "multi-element list cannot coerce to %s: %r"
            % (expected_type.__name__, value))
    if expected_type is str and isinstance(value, dict):
        raise ValueError(
            "dict cannot coerce to str: %r" % (value,))
    return expected_type(value)


# v119.H2.1: values that count as "skip this program" at the
# program_settings level.  LLMs fluctuate between native JSON bool
# and string variants depending on provider/temperature/token-
# boundary effects; this list catches the common permutations so
# Google's "skip": true and OpenAI's "skip": "true" both promote
# deterministically.  Per Gemini guardrail A.
#
# Note on equality semantics: Python treats 1.0 == 1 == True, so
# any of these compare-equal forms also match this membership
# check.  That's fine -- float 1.0 from an LLM still means "truthy".
_SKIP_TRUE_VALUES = (True, "true", "True", "TRUE",
                     "yes", "Yes", "YES",
                     1, "1",
                     "skip", "Skip", "SKIP")


def _promote_skip_settings_to_skip_programs(directives, log=None):
    """Promote program_settings[X].skip=truthy to skip_programs.

    v119.H2.1.  Fixes a pre-existing gap where the LLM correctly
    encodes "don't run X" as program_settings[X].skip=true but the
    downstream consumer reads skip_programs from
    workflow_preferences.

    In-place mutation.  Modifies `directives`.

    For each program_settings entry whose settings dict contains a
    'skip' key with a value in _SKIP_TRUE_VALUES:
      1. Append the program name to workflow_preferences.skip_programs
         (creating the list if needed; de-duplicates).
      2. Pop the 'skip' key from the program's settings dict (so the
         downstream per-setting validation does NOT log "Unknown
         setting skip=True (keeping)").
      3. If the program's settings dict is now empty, remove the
         program entry entirely (matches what the regex fallback
         in extract_directives_simple produces for the same input,
         keeps output tidy).

    IMPORTANT -- call order in validate_directives:
      This helper MUST run BEFORE the per-setting validation loop
      iterates settings.items().  The wiring at the top of
      validate_directives ensures this.  Without that ordering,
      the "Unknown setting skip=True (keeping)" log line would
      fire even when the value gets promoted, polluting logs and
      breaking the no-leftover-log K_H2.1 assertion.

    Defensive bail policy:
      If workflow_preferences or skip_programs exist on the input
      but have unexpected types (workflow_preferences not a dict,
      skip_programs not a list), this helper leaves ALL
      program_settings entries unchanged (including their skip
      keys).  The downstream per-setting validation will then log
      the legacy "Unknown setting skip=True (keeping)" line, which
      is the right signal that the LLM emitted malformed
      workflow_preferences.  Partial processing would be worse
      than no processing in this corner case.

    Liberal on input acceptance, strict on output shape: any value
    in _SKIP_TRUE_VALUES counts as truthy; the promoted result is
    always a list of strings under workflow_preferences.skip_programs.

    Args:
        directives: dict.  Modified in place.
        log: optional logging callable (signature: log(msg_str))

    Returns:
        None.  The mutation is the side effect.
    """
    def _log(msg):
        if log:
            log(msg)

    if not isinstance(directives, dict):
        return
    prog_settings = directives.get("program_settings")
    if not isinstance(prog_settings, dict):
        return

    # Collect promotions first.  Don't mutate prog_settings while
    # iterating its keys.
    to_promote = []  # list of program names whose skip-flag we will pop
    for prog, settings in prog_settings.items():
        if not isinstance(settings, dict):
            continue
        # Filter empty / non-string program names.  If the LLM
        # somehow emits program_settings[""] = {"skip": true},
        # promoting "" to skip_programs would put junk through
        # the downstream pipeline.
        if not prog or not isinstance(prog, str):
            continue
        if "skip" not in settings:
            continue
        if settings["skip"] in _SKIP_TRUE_VALUES:
            to_promote.append(prog)

    if not to_promote:
        return

    # Ensure workflow_preferences.skip_programs exists.  Defensive
    # bail: if either is present but malformed, return without
    # popping any skip keys.  Leaves the LLM's malformed output
    # visible to downstream validation.
    if "workflow_preferences" not in directives:
        directives["workflow_preferences"] = {}
    wf = directives["workflow_preferences"]
    if not isinstance(wf, dict):
        _log("DIRECTIVES: workflow_preferences is not a dict; "
             "skip-promotion skipped (input was malformed)")
        return
    if "skip_programs" not in wf:
        wf["skip_programs"] = []
    elif not isinstance(wf["skip_programs"], list):
        _log("DIRECTIVES: skip_programs is not a list; "
             "skip-promotion skipped (input was malformed)")
        return

    for prog in to_promote:
        settings = prog_settings[prog]
        # Pop (not delete-after-read): cleaner idiomatic mutation,
        # guarantees the key is gone before any downstream code
        # sees it.
        settings.pop("skip", None)
        # Add to skip_programs (de-duplicated).
        if prog not in wf["skip_programs"]:
            wf["skip_programs"].append(prog)
            _log("DIRECTIVES: Promoted %s from program_settings.skip "
                 "to workflow_preferences.skip_programs" % prog)
        else:
            _log("DIRECTIVES: Dropped redundant program_settings.skip "
                 "for %s (already in skip_programs)" % prog)
        # If the program's settings dict is now empty, remove the
        # program entry.
        if not settings:
            del prog_settings[prog]


def validate_directives(directives, log=None):
    """
    Validate and clean extracted directives.

    Removes invalid entries, converts types, and ensures consistency.

    Args:
        directives: Raw extracted directives dict
        log: Optional logging function

    Returns:
        dict: Validated directives
    """
    def _log(msg):
        if log:
            log(msg)

    if not directives or not isinstance(directives, dict):
        return {}

    # v119.H2.1: Promote LLM-emitted program_settings[X].skip=true
    # to workflow_preferences.skip_programs[X] BEFORE per-setting
    # validation runs.  Must run before the per-setting validation
    # loop below so the popped 'skip' key never reaches the
    # "Unknown setting skip=True (keeping)" log line later in this
    # function.  See _promote_skip_settings_to_skip_programs for
    # full rationale.
    _promote_skip_settings_to_skip_programs(directives, log=_log)

    validated = {}

    # Validate program_settings
    if "program_settings" in directives:
        prog_settings = directives["program_settings"]
        if isinstance(prog_settings, dict):
            valid_prog_settings = {}

            for prog, settings in prog_settings.items():
                # Check program name (allow "default")
                if prog != "default" and prog not in VALID_PROGRAMS:
                    # Try to fix common variations
                    fixed_prog = _fix_program_name(prog)
                    if fixed_prog:
                        prog = fixed_prog
                    else:
                        _log("DIRECTIVES: Ignoring unknown program %s" % prog)
                        continue

                if not isinstance(settings, dict):
                    continue

                # v118.B2: Heal/drop namespaced PHIL keys BEFORE
                # type-validation.  Translates known near-misses
                # (e.g. data_manager.r_free_flags.generate ->
                # generate_rfree_flags) and drops keys in bad-PHIL
                # namespace prefixes.  Program-specific dotted keys
                # like "autosol.atom_type" pass through.
                settings = _heal_namespaced_phil_keys(
                    settings, prog, _log)

                valid_settings = {}
                for key, value in settings.items():
                    if key in VALID_SETTINGS:
                        try:
                            expected_type = VALID_SETTINGS[key]
                            # v118.10: coerce list/tuple shapes
                            # (Google emits lists where OpenAI emits
                            # scalars) BEFORE the type cast.  Bare
                            # `expected_type(value)` would call
                            # ``str(["S"])`` which produces
                            # ``"['S']"`` — the original bug.
                            coerced = _coerce_setting_value(
                                value, expected_type)
                            if isinstance(value, (list, tuple)):
                                # Distinct log lines per content:
                                # single-element = JSON-shape artifact;
                                # multi-element = content-ambiguity
                                # signal.  Both worth visibility.
                                if len(value) <= 1:
                                    _log(
                                        "DIRECTIVES: Unpacked "
                                        "single-element list for "
                                        "%s: %r → %r"
                                        % (key, value, coerced))
                                else:
                                    _log(
                                        "DIRECTIVES: Joined "
                                        "multi-element list for "
                                        "%s: %r → %r"
                                        % (key, value, coerced))
                            valid_settings[key] = coerced
                        except (ValueError, TypeError):
                            _log("DIRECTIVES: Invalid value for %s: %s" % (key, value))
                    else:
                        # Keep unknown settings but log them
                        valid_settings[key] = value
                        _log("DIRECTIVES: Unknown setting %s=%s (keeping)" % (key, value))

                # Normalize unit_cell to a plain space-separated string regardless
                # of how the LLM formatted it (parentheses, commas, etc.)
                if "unit_cell" in valid_settings:
                    raw_uc = valid_settings["unit_cell"]
                    normalized = _normalize_unit_cell(str(raw_uc))
                    if normalized:
                        if normalized != raw_uc:
                            _log("DIRECTIVES: Normalised unit_cell %r → %r" % (raw_uc, normalized))
                        valid_settings["unit_cell"] = normalized
                    else:
                        _log("DIRECTIVES: Dropping malformed unit_cell %r (need 6 numbers)" % raw_uc)
                        del valid_settings["unit_cell"]

                # Validate space_group: reject placeholder/non-crystallographic values
                # v116.15: Use the shared _is_symmetry_sentinel helper so the
                # sentinel set stays in sync with _apply_crystal_symmetry_fallback.
                # Structural checks (must start with a letter, length) remain here
                # because they are LLM-output specific (the regex path already
                # guarantees the shape by construction).
                #
                # v119.H14: Added positive Hermann-Mauguin check (_HM_FORM_RE).
                # The pre-H14 negative checks (sentinel + first-letter +
                # length<=25) let multi-word phrases like "Solve the
                # structure" or "From data file" through, since they start
                # with a letter and are <=25 chars.  Such values then reached
                # PHENIX as bogus PHIL params and produced obscure errors.
                # The HM check positively confirms the value looks like a
                # space-group symbol (lattice letter + digits/spaces) before
                # accepting it.  Forms accepted: "P 1 21 1", "P21", "P212121",
                # "P-1", "F222", "P 1 21 1 (No. 4)".  Truly exotic forms
                # (e.g., "R 3 :H" with axis hint) will be dropped — those
                # are rare; PHENIX can usually detect from data.
                if "space_group" in valid_settings:
                    raw_sg = str(valid_settings["space_group"]).strip()
                    _sg_invalid = (
                        not raw_sg
                        or _is_symmetry_sentinel(raw_sg)
                        # Must start with a letter (space group symbol)
                        or not raw_sg[0].isalpha()
                        # Must be short enough to be a real symbol
                        or len(raw_sg) > 25
                        # v119.H14: must match Hermann-Mauguin shape
                        or not _looks_like_space_group(raw_sg)
                    )
                    if _sg_invalid:
                        _log("DIRECTIVES: Dropping invalid space_group value: %r" % raw_sg)
                        del valid_settings["space_group"]

                if valid_settings:
                    valid_prog_settings[prog] = valid_settings

            if valid_prog_settings:
                validated["program_settings"] = valid_prog_settings

        # --- Sanity check: resolution vs wavelength confusion ---
        # X-ray wavelengths (0.5-2.5 Å) can be confused with resolution.
        # Collect all wavelength values across all programs first.
        if "program_settings" in validated:
            all_wavelengths = set()
            for prog, settings in validated["program_settings"].items():
                wl = settings.get("wavelength")
                if wl is not None:
                    all_wavelengths.add(wl)

            for prog, settings in list(validated["program_settings"].items()):
                res = settings.get("resolution")
                if res is not None:
                    remove = False
                    reason = ""
                    # Check if resolution matches any wavelength
                    for wl in all_wavelengths:
                        if abs(res - wl) < 0.01:
                            remove = True
                            reason = ("matches wavelength=%.4f" % wl)
                            break
                    # Resolution < 1.2 Å is almost certainly a wavelength
                    # (very few structures resolve below 1.2 Å)
                    if not remove and res < 1.2:
                        remove = True
                        reason = ("%.4f Å is in X-ray wavelength range, "
                                  "not a plausible resolution" % res)
                    if remove:
                        _log("DIRECTIVES: Removing resolution=%.4f from %s "
                             "(%s)" % (res, prog, reason))
                        del settings["resolution"]

            # Clean up empty settings dicts
            validated["program_settings"] = {
                prog: s for prog, s in validated["program_settings"].items() if s
            }
            if not validated["program_settings"]:
                del validated["program_settings"]

    # Validate stop_conditions
    if "stop_conditions" in directives:
        stop_cond = directives["stop_conditions"]
        if isinstance(stop_cond, dict):
            valid_stop = {}

            for key, value in stop_cond.items():
                if key in VALID_STOP_CONDITIONS:
                    try:
                        expected_type = VALID_STOP_CONDITIONS[key]
                        # v118.10: coerce list/tuple shapes BEFORE
                        # the `value not in VALID_PROGRAMS` check
                        # below — that check raises TypeError on a
                        # list input ("unhashable type: 'list'"),
                        # which the outer except would swallow,
                        # silently dropping the directive.
                        coerced = _coerce_setting_value(
                            value, expected_type)
                        if isinstance(value, (list, tuple)):
                            if len(value) <= 1:
                                _log(
                                    "DIRECTIVES: Unpacked "
                                    "single-element list for "
                                    "%s: %r → %r"
                                    % (key, value, coerced))
                            else:
                                _log(
                                    "DIRECTIVES: Joined "
                                    "multi-element list for "
                                    "%s: %r → %r"
                                    % (key, value, coerced))
                        if key == "after_program":
                            # Validate program name (now safe — the
                            # coerced value is a string).
                            if coerced not in VALID_PROGRAMS:
                                fixed = _fix_program_name(coerced)
                                if fixed:
                                    coerced = fixed
                                else:
                                    _log("DIRECTIVES: Invalid stop program %s" % coerced)
                                    continue
                        valid_stop[key] = coerced
                    except (ValueError, TypeError):
                        _log("DIRECTIVES: Invalid stop condition %s: %s" % (key, value))
                else:
                    _log("DIRECTIVES: Unknown stop condition %s" % key)

            if valid_stop:
                validated["stop_conditions"] = valid_stop

    # Validate file_preferences
    if "file_preferences" in directives:
        file_prefs = directives["file_preferences"]
        if isinstance(file_prefs, dict):
            valid_prefs = {}

            for key, value in file_prefs.items():
                if key in ("model", "sequence", "data_mtz", "map_coeffs_mtz", "map"):
                    # File path preferences
                    if isinstance(value, str) and value:
                        valid_prefs[key] = value
                elif key == "mtz":
                    # Backward compat: "mtz" -> "data_mtz"
                    if isinstance(value, str) and value:
                        valid_prefs["data_mtz"] = value
                elif key == "exclude":
                    # Exclusion list
                    if isinstance(value, list):
                        valid_prefs[key] = [str(v) for v in value if v]
                elif key in ("prefer_anomalous", "prefer_unmerged", "prefer_merged"):
                    # Boolean data type preferences.
                    # v119.H5.1.1: defend against LLM emitting
                    # list-wrapped booleans (e.g. [False] which
                    # bare bool() would incorrectly treat as
                    # truthy because non-empty list is truthy).
                    # Type-gated by isinstance(v[0], bool) so
                    # adding a list-typed key to this tuple in
                    # the future wouldn't trigger silent
                    # flattening — only genuine [bool] unwraps.
                    _v = value
                    if (isinstance(_v, list) and len(_v) == 1
                            and isinstance(_v[0], bool)):
                        _v = _v[0]
                    if isinstance(_v, bool):
                        valid_prefs[key] = _v
                    elif isinstance(_v, int) and not isinstance(_v, bool):
                        # Bare int 0/1 preserved (existing
                        # behaviour).  bool is a subclass of int
                        # in Python, but the case above already
                        # caught it.
                        valid_prefs[key] = bool(_v)
                    elif isinstance(_v, str):
                        # Bare string preserved (existing
                        # behaviour, including the legacy quirk
                        # that bool("false") == True; intentional
                        # — separate bug class, separate ship).
                        valid_prefs[key] = bool(_v)
                    else:
                        _log("DIRECTIVES: Invalid value for "
                             "%s: %r" % (key, value))

            if valid_prefs:
                validated["file_preferences"] = valid_prefs

    # Validate workflow_preferences
    if "workflow_preferences" in directives:
        workflow_prefs = directives["workflow_preferences"]
        if isinstance(workflow_prefs, dict):
            valid_wf = {}

            for key, value in workflow_prefs.items():
                if key in ("skip_programs", "prefer_programs"):
                    if isinstance(value, list):
                        # Validate program names in list
                        valid_list = []
                        for prog in value:
                            if prog in VALID_PROGRAMS:
                                valid_list.append(prog)
                            else:
                                fixed = _fix_program_name(prog)
                                if fixed:
                                    valid_list.append(fixed)
                        if valid_list:
                            valid_wf[key] = valid_list
                elif key in ("use_experimental_phasing", "use_molecular_replacement",
                             "use_mr_sad", "model_is_placed",
                             "wants_validation_only"):
                    # v119.H5.1.1: same list-wrap defense as the
                    # file_preferences boolean block above — see
                    # there for full rationale.
                    _v = value
                    if (isinstance(_v, list) and len(_v) == 1
                            and isinstance(_v[0], bool)):
                        _v = _v[0]
                    if isinstance(_v, bool):
                        valid_wf[key] = _v
                    elif isinstance(_v, int) and not isinstance(_v, bool):
                        valid_wf[key] = bool(_v)
                    elif isinstance(_v, str):
                        valid_wf[key] = bool(_v)
                    else:
                        _log("DIRECTIVES: Invalid value for "
                             "%s: %r" % (key, value))

            if valid_wf:
                validated["workflow_preferences"] = valid_wf

    # Keep constraints as-is (list of strings)
    if "constraints" in directives:
        constraints = directives["constraints"]
        if isinstance(constraints, list):
            valid_constraints = [str(c) for c in constraints if c]
            if valid_constraints:
                validated["constraints"] = valid_constraints

    # POST-PROCESSING: Check for conflicting directives
    # If constraints mention ligand fitting AND stop_conditions has after_program=phenix.refine,
    # this is likely wrong - the user wants the full ligand workflow, not to stop at first refine
    validated = _fix_ligand_workflow_conflict(validated, _log)

    # If constraints describe steps beyond the after_program target, clear it
    # e.g., after_program=phenix.map_symmetry but constraints say "build a model" → don't stop
    validated = _fix_multi_step_workflow_conflict(validated, _log)

    # Fix after_cycle=1 when max_refine_cycles is set
    # LLM often confuses "one refinement cycle" with "one agent cycle"
    validated = _fix_after_cycle_refinement_conflict(validated, _log)

    # Fix contradiction: after_program set to a program that's in skip_programs
    # e.g., "don't run mtriage" but stop_conditions.after_program=phenix.mtriage
    validated = _fix_skip_after_program_conflict(validated, _log)

    return validated


def _fix_ligand_workflow_conflict(directives, log):
    """
    Fix conflicting directives when ligand fitting workflow is detected.

    If user wants ligand fitting after refinement, we should NOT stop at phenix.refine
    or at a low cycle number, because that would stop before the ligand fitting happens.

    A typical ligand workflow is:
      1. xtriage
      2. predict_and_build
      3. process_predicted_model (maybe)
      4. phaser
      5. refine (first)
      6. ligandfit
      7. pdbtools (combine)
      8. refine (second) <- stop here

    So "stop after second refinement" needs at least 8 cycles, not 2.
    """
    constraints = directives.get("constraints", [])
    stop_conditions = directives.get("stop_conditions", {})

    if not stop_conditions:
        return directives

    after_program = stop_conditions.get("after_program", "")
    after_cycle = stop_conditions.get("after_cycle")

    # Check if constraints mention ligand fitting
    ligand_keywords = ["ligand", "ligandfit", "fit ligand", "place ligand", "lig.pdb"]
    has_ligand_constraint = any(
        any(kw in c.lower() for kw in ligand_keywords)
        for c in constraints
    ) if constraints else False

    # If we have ligand constraints AND after_program is phenix.refine,
    # this is likely a conflict - user wants ligand workflow to complete
    if has_ligand_constraint and after_program == "phenix.refine":
        log("DIRECTIVES: Clearing after_program=phenix.refine due to ligand workflow in constraints")
        log("DIRECTIVES: (User wants full ligand fitting workflow, not to stop at first refinement)")

        # v116.x: clear stop_after_requested and skip_validation
        # since the user-stop intent is being overridden by the
        # ligand-workflow conflict.
        del directives["stop_conditions"]["after_program"]
        directives["stop_conditions"].pop("stop_after_requested", None)
        directives["stop_conditions"].pop("skip_validation", None)

        # If stop_conditions is now empty, remove it
        if not directives["stop_conditions"]:
            del directives["stop_conditions"]

    # If we have ligand constraints AND after_cycle is too low (<=4),
    # the LLM probably misinterpreted "second refinement" as "cycle 2"
    # A ligand workflow needs at least 6-8 cycles to complete
    if has_ligand_constraint and after_cycle is not None and after_cycle <= 4:
        log(f"DIRECTIVES: Clearing after_cycle={after_cycle} due to ligand workflow in constraints")
        log("DIRECTIVES: (Ligand workflow needs ~8 cycles; 'second refinement' != 'cycle 2')")

        # Remove the conflicting after_cycle
        del directives["stop_conditions"]["after_cycle"]

        # If stop_conditions is now empty, remove it
        if not directives["stop_conditions"]:
            del directives["stop_conditions"]

    return directives


def _fix_multi_step_workflow_conflict(directives, log):
    """
    Fix conflicting directives when constraints describe a multi-step workflow
    but after_program would stop at an intermediate step.

    For example:
    - after_program=phenix.map_symmetry but constraints say "build a model" → clear after_program
    - after_program=phenix.map_sharpening but constraints say "apply NCS" → clear after_program

    This catches cases where the LLM extraction sets after_program for an intermediate
    step even though the user's constraints describe a longer pipeline.
    """
    constraints = directives.get("constraints", [])
    stop_conditions = directives.get("stop_conditions", {})

    if not stop_conditions or not constraints:
        return directives

    after_program = stop_conditions.get("after_program", "")
    if not after_program:
        return directives

    constraints_text = " ".join(str(c).lower() for c in constraints if c)

    # Define which programs are "intermediate" and what downstream keywords indicate
    # the workflow should continue past them
    intermediate_programs = {
        "phenix.predict_and_build": [
            # predict_and_build is almost always intermediate — the user wants
            # MR + building + refinement after prediction.  Only stop here if the
            # user explicitly asks "just predict" or "only prediction".
            r'solv',            # "solve the structure"
            r'refine',
            r'build',
            r'molecular\s+replacement',
            r'phaser',
            r'MR\b',
            r'model\s+building',
            r'ligand',
            r'structure\s+(?:determination|solution)',
        ],
        "phenix.map_symmetry": [
            r'build\s+(?:a\s+)?model',
            r'model\s+building',
            r'map[\s_]to[\s_]model',
            r'dock.*(?:in|into)',
            r'refine',
            r'sharpen',
            r'density\s+modif',
            r'resolve_cryo_em',
            r'apply\s+(?:ncs|symmetry)',
            r'generate.*(?:complex|assembly|full|complete)',
            r'extract\s+(?:the\s+)?(?:unique|asymmetric)',
        ],
        "phenix.map_sharpening": [
            r'build\s+(?:a\s+)?model',
            r'model\s+building',
            r'map[\s_]to[\s_]model',
            r'dock.*(?:in|into)',
            r'refine',
            r'(?:find|determine|detect)\s+symmetry',
            r'map[\s_]symmetry',
            r'apply\s+(?:ncs|symmetry)',
            r'generate.*(?:complex|assembly|full|complete)',
        ],
        "phenix.mtriage": [
            r'build\s+(?:a\s+)?model',
            r'refine',
            r'dock.*(?:in|into)',
            r'sharpen',
            r'density\s+modif',
            r'resolve_cryo_em',
            r'map[\s_]to[\s_]model',
        ],
        "phenix.xtriage": [
            r'phaser',
            r'molecular\s+replacement',
            r'refine',
            r'build',
            r'ligand',
        ],
        "phenix.phaser": [
            r'refine',
            r'ligand',
            r'build',
            r'validate',
        ],
    }

    downstream_patterns = intermediate_programs.get(after_program, [])
    if not downstream_patterns:
        return directives

    has_downstream = any(
        re.search(pattern, constraints_text, re.IGNORECASE)
        for pattern in downstream_patterns
    )

    if has_downstream:
        log("DIRECTIVES: Clearing after_program=%s because constraints describe downstream work" %
            after_program)
        log("DIRECTIVES: Constraints: %s" % constraints_text[:200])

        del directives["stop_conditions"]["after_program"]
        # v116.x: clear stop_after_requested too — the user-stop
        # intent is being overridden because the constraints
        # contradict it (they describe work AFTER the named program).
        directives["stop_conditions"].pop("stop_after_requested", None)

        # If stop_conditions is now empty (or only has skip_validation), clean up
        remaining = {k: v for k, v in directives["stop_conditions"].items()
                     if k != "skip_validation"}
        if not remaining:
            # Remove skip_validation too since there's no stop to validate
            if "skip_validation" in directives["stop_conditions"]:
                del directives["stop_conditions"]["skip_validation"]
        if not directives["stop_conditions"]:
            del directives["stop_conditions"]

    return directives


def _fix_after_cycle_refinement_conflict(directives, log):
    """
    Fix conflict where LLM sets after_cycle=1 alongside max_refine_cycles.

    When user says "stop after refinement" or "one refinement cycle", the LLM
    sometimes produces both max_refine_cycles=1 AND after_cycle=1. The after_cycle=1
    is wrong — it would stop after the first agent cycle (e.g., xtriage), not after
    the refinement program completes.

    Rules:
    - after_cycle=1 + max_refine_cycles=N → remove after_cycle (keep max_refine_cycles)
    - after_cycle=1 alone (no max_refine_cycles, no after_program) → suspicious,
      only keep if no workflow programs would run before refinement
    """
    stop_conditions = directives.get("stop_conditions", {})
    if not stop_conditions:
        return directives

    after_cycle = stop_conditions.get("after_cycle")
    max_refine = stop_conditions.get("max_refine_cycles")
    after_program = stop_conditions.get("after_program")

    # Case 1: after_cycle + max_refine_cycles — the after_cycle is redundant/wrong
    if after_cycle is not None and max_refine is not None:
        log("DIRECTIVES: Removing after_cycle=%d (conflicts with max_refine_cycles=%d)" %
            (after_cycle, max_refine))
        log("DIRECTIVES: (LLM confused 'one refinement' with 'one agent cycle')")
        del directives["stop_conditions"]["after_cycle"]

    # Case 2: after_cycle=1 alone (no after_program, no max_refine)
    # This almost always means the LLM misinterpreted "stop after refinement"
    # as "stop after 1 cycle". Convert to max_refine_cycles=1.
    elif after_cycle == 1 and after_program is None and max_refine is None:
        log("DIRECTIVES: Converting after_cycle=1 → max_refine_cycles=1")
        log("DIRECTIVES: (after_cycle=1 would stop after xtriage, not after refinement)")
        del directives["stop_conditions"]["after_cycle"]
        directives["stop_conditions"]["max_refine_cycles"] = 1

    return directives


def _fix_skip_after_program_conflict(directives, log):
    """
    Fix contradiction where after_program is set to a program the user wants skipped.

    This catches two cases:
    1. after_program is in skip_programs (explicit)
    2. constraints say "do not run X" or "skip X" and after_program references X

    When user says "don't run mtriage" but the LLM sets after_program=phenix.mtriage,
    this is clearly wrong. The LLM saw "mtriage" mentioned and associated it with the
    stop condition, ignoring the negation.

    Fix: Remove after_program if it conflicts, and try to infer the correct stop
    based on what the user actually wants.
    """
    stop_conditions = directives.get("stop_conditions", {})
    workflow_prefs = directives.get("workflow_preferences", {})
    constraints = directives.get("constraints", [])

    after_program = stop_conditions.get("after_program", "")
    skip_programs = workflow_prefs.get("skip_programs", [])

    if not after_program:
        return directives

    after_normalized = after_program.replace("phenix.", "")

    # Check 1: after_program in skip_programs
    is_skipped = False
    for skip_prog in skip_programs:
        skip_normalized = skip_prog.replace("phenix.", "")
        if after_normalized == skip_normalized:
            is_skipped = True
            break

    # Check 2: constraints mention "do not run X", "skip X", "don't run X"
    if not is_skipped and constraints:
        import re
        skip_patterns = [
            r"(?:do\s+not|don't|skip|no|never|avoid)\s+(?:run(?:ning)?|use|execute)\s+(?:phenix\.)?(\w+)",
            r"(?:phenix\.)?(\w+)\s+(?:should\s+)?(?:not\s+be|is\s+not)\s+(?:run|used|needed)",
        ]
        for constraint in constraints:
            constraint_lower = constraint.lower()
            for pattern in skip_patterns:
                m = re.search(pattern, constraint_lower)
                if m:
                    skipped_prog = m.group(1).replace("phenix.", "")
                    if skipped_prog == after_normalized:
                        is_skipped = True
                        log("DIRECTIVES: Constraint '%s' contradicts after_program=%s" %
                            (constraint, after_program))
                        break
            if is_skipped:
                break

    if is_skipped:
        log("DIRECTIVES: CONFLICT - after_program=%s contradicts skip/constraint!" % after_program)
        log("DIRECTIVES: Removing after_program (user explicitly said don't run this program)")

        del directives["stop_conditions"]["after_program"]
        # v116.x: clear stop_after_requested too — the user-stop
        # intent is being overridden because the user also said
        # don't run this program.
        directives["stop_conditions"].pop("stop_after_requested", None)

        # Try to infer correct stop: if user wanted refinement, set max_refine_cycles
        # Look for refinement-related terms in constraints, program_settings, or
        # the removed after_program's context
        if "max_refine_cycles" not in stop_conditions:
            inferred = False
            prog_settings = directives.get("program_settings", {})
            rsr_settings = prog_settings.get("phenix.real_space_refine", {})
            refine_settings = prog_settings.get("phenix.refine", {})

            # Check program_settings for cycle count
            if rsr_settings.get("cycles") or rsr_settings.get("macro_cycles"):
                log("DIRECTIVES: Inferring max_refine_cycles=1 from RSR settings")
                directives["stop_conditions"]["max_refine_cycles"] = 1
                inferred = True
            elif refine_settings.get("cycles") or refine_settings.get("macro_cycles"):
                log("DIRECTIVES: Inferring max_refine_cycles=1 from refine settings")
                directives["stop_conditions"]["max_refine_cycles"] = 1
                inferred = True

            # Check constraints for refinement-related stop language
            if not inferred and constraints:
                import re
                refine_stop_patterns = [
                    r'stop\s+after\s+(?:the\s+)?(?:single\s+)?(?:macro[- ]?cycle|refinement|refine)',
                    r'(?:one|single|1)\s+(?:macro[- ]?cycle|refinement)',
                    r'(?:only|just)\s+(?:one|a\s+single)\s+(?:refinement|refine)',
                ]
                for constraint in constraints:
                    for pattern in refine_stop_patterns:
                        if re.search(pattern, constraint.lower()):
                            log("DIRECTIVES: Inferring max_refine_cycles=1 from constraint: '%s'" % constraint)
                            directives["stop_conditions"]["max_refine_cycles"] = 1
                            inferred = True
                            break
                    if inferred:
                        break

            # If the removed program was a validation/analysis tool (not refinement),
            # the user likely wanted to stop after refinement completes
            if not inferred:
                removed_normalized = after_normalized
                analysis_programs = {"mtriage", "xtriage", "molprobity", "validation_cryoem"}
                if removed_normalized in analysis_programs:
                    log("DIRECTIVES: Removed after_program was analysis tool (%s)" % after_program)
                    log("DIRECTIVES: Inferring max_refine_cycles=1 (stop after refinement)")
                    directives["stop_conditions"]["max_refine_cycles"] = 1

        # If stop_conditions is now empty (except maybe skip_validation), clean up
        remaining = {k: v for k, v in directives["stop_conditions"].items()
                     if k != "skip_validation"}
        if not remaining:
            # Keep skip_validation if it was set
            if directives["stop_conditions"].get("skip_validation"):
                pass  # Keep the dict with just skip_validation
            else:
                del directives["stop_conditions"]

    return directives


def _fix_program_name(name):
    """
    Try to fix common variations in program names.

    Args:
        name: Potentially incorrect program name

    Returns:
        str: Corrected program name, or None if not fixable
    """
    if not name:
        return None

    name_lower = name.lower().strip()

    # Remove "phenix." prefix if checking without it
    if name_lower.startswith("phenix."):
        name_lower = name_lower[7:]

    # Common mappings - include various case/format variations
    mappings = {
        # Refinement
        "refine": "phenix.refine",
        "refinement": "phenix.refine",

        # Phasing
        "autosol": "phenix.autosol",
        "phaser": "phenix.phaser",

        # Model building - X-ray
        "autobuild": "phenix.autobuild",
        "auto_build": "phenix.autobuild",

        # Model building - Cryo-EM
        "map_to_model": "phenix.map_to_model",
        "maptomodel": "phenix.map_to_model",
        "map-to-model": "phenix.map_to_model",
        "build_model": "phenix.map_to_model",  # Common alternative phrasing
        "buildmodel": "phenix.map_to_model",

        # Docking
        "dock_in_map": "phenix.dock_in_map",
        "dockinmap": "phenix.dock_in_map",
        "dock-in-map": "phenix.dock_in_map",

        # AlphaFold/prediction
        "predict_and_build": "phenix.predict_and_build",
        "predictandbuild": "phenix.predict_and_build",
        "predict-and-build": "phenix.predict_and_build",
        "process_predicted_model": "phenix.process_predicted_model",
        "processpredictedmodel": "phenix.process_predicted_model",

        # Real space refinement
        "real_space_refine": "phenix.real_space_refine",
        "realspacerefine": "phenix.real_space_refine",
        "real-space-refine": "phenix.real_space_refine",
        "rsr": "phenix.real_space_refine",

        # Density modification
        "resolve_cryo_em": "phenix.resolve_cryo_em",
        "resolvecryoem": "phenix.resolve_cryo_em",
        "resolve-cryo-em": "phenix.resolve_cryo_em",
        # LLM often shortens to "resolve" — map to cryo_em
        # (legacy phenix.resolve is not in the agent's
        # program list; all density modification goes
        # through resolve_cryo_em)
        "resolve": "phenix.resolve_cryo_em",
        "autobuild_denmod": "phenix.autobuild_denmod",

        # Map tools
        "map_sharpening": "phenix.map_sharpening",
        "mapsharpening": "phenix.map_sharpening",
        "sharpen_map": "phenix.map_sharpening",  # Alternative phrasing
        "sharpenmap": "phenix.map_sharpening",
        "auto_sharpen": "phenix.map_sharpening",  # Alternative name
        "autosharpen": "phenix.map_sharpening",
        "map_symmetry": "phenix.map_symmetry",
        "mapsymmetry": "phenix.map_symmetry",

        # Validation
        "molprobity": "phenix.molprobity",
        "model_vs_data": "phenix.model_vs_data",
        "modelvsdata": "phenix.model_vs_data",
        "validation_cryoem": "phenix.validation_cryoem",
        "validationcryoem": "phenix.validation_cryoem",
        "holton_geometry_validation": "phenix.holton_geometry_validation",

        # Map analysis
        "polder": "phenix.polder",
        "polder_map": "phenix.polder",
        "omit_map": "phenix.polder",

        # Ligand fitting
        "ligandfit": "phenix.ligandfit",
        "ligand_fit": "phenix.ligandfit",

        # Analysis
        "xtriage": "phenix.xtriage",
        "mtriage": "phenix.mtriage",

        # PDB tools
        "pdbtools": "phenix.pdbtools",
        "pdb_tools": "phenix.pdbtools",
    }

    return mappings.get(name_lower)


# =============================================================================
# MERGING
# =============================================================================

def merge_directives(base, override):
    """
    Merge two directive dicts, with override taking precedence.

    Args:
        base: Base directives dict
        override: Override directives dict (takes precedence)

    Returns:
        dict: Merged directives
    """
    if not base:
        return override or {}
    if not override:
        return base or {}

    merged = {}

    # Merge program_settings (deep merge)
    base_prog = base.get("program_settings", {})
    over_prog = override.get("program_settings", {})
    if base_prog or over_prog:
        merged_prog = dict(base_prog)
        for prog, settings in over_prog.items():
            if prog in merged_prog:
                merged_prog[prog] = {**merged_prog[prog], **settings}
            else:
                merged_prog[prog] = settings
        merged["program_settings"] = merged_prog

    # Merge stop_conditions (override wins)
    base_stop = base.get("stop_conditions", {})
    over_stop = override.get("stop_conditions", {})
    if base_stop or over_stop:
        # v119.H5.1.1: protect after_program corrections.
        # The _corrected_from sidecar describes the CURRENT
        # value of after_program. If override changes
        # after_program:
        #   - to the same value the correction reverted away
        #     from (over_stop[after_program] ==
        #     _corrected_from.from): prevent the un-correction
        #     by stripping from override; sidecar stays valid.
        #   - to any other value, AND override doesn't bring
        #     its own sidecar: let override win for
        #     after_program but clear base's now-stale sidecar
        #     (Gemini blanket-wipe).  Otherwise the dict-merge
        #     would leave a "zombie" sidecar whose `to` field
        #     no longer describes the result.
        # The dict comprehension rebuilds preserve the
        # merge_directives no-mutation-of-inputs invariant
        # (same pattern as the program_settings deep-merge
        # above uses dict()).
        base_cf = base_stop.get("_corrected_from")
        if base_cf and "after_program" in over_stop:
            if over_stop["after_program"] == base_cf.get("from"):
                over_stop = {k: v for k, v in over_stop.items()
                             if k != "after_program"}
            elif "_corrected_from" not in over_stop:
                base_stop = {k: v for k, v in base_stop.items()
                             if k != "_corrected_from"}
        merged["stop_conditions"] = {**base_stop, **over_stop}

    # Merge file_preferences (override wins)
    base_files = base.get("file_preferences", {})
    over_files = override.get("file_preferences", {})
    if base_files or over_files:
        merged["file_preferences"] = {**base_files, **over_files}

    # Merge workflow_preferences (override wins for bools, extend for lists)
    base_wf = base.get("workflow_preferences", {})
    over_wf = override.get("workflow_preferences", {})
    if base_wf or over_wf:
        merged_wf = dict(base_wf)
        for key, value in over_wf.items():
            if key in ("skip_programs", "prefer_programs"):
                # Extend lists
                existing = merged_wf.get(key, [])
                merged_wf[key] = list(set(existing + value))
            else:
                merged_wf[key] = value
        merged["workflow_preferences"] = merged_wf

    # Merge constraints (concatenate)
    base_const = base.get("constraints", [])
    over_const = override.get("constraints", [])
    if base_const or over_const:
        merged["constraints"] = base_const + over_const

    return merged


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_program_settings(directives, program_name):
    """
    Get settings for a specific program, with default fallback.

    Args:
        directives: Directives dict
        program_name: Name of program (e.g., "phenix.refine")

    Returns:
        dict: Settings for the program (merged with defaults)
    """
    if not directives:
        return {}

    prog_settings = directives.get("program_settings", {})

    # Get default settings
    default = dict(prog_settings.get("default", {}))

    # Get program-specific settings (override defaults)
    specific = prog_settings.get(program_name, {})

    return {**default, **specific}


def check_stop_conditions(directives, cycle_number, last_program, metrics=None):
    """
    Check if any stop conditions are met.

    Args:
        directives: Directives dict
        cycle_number: Current cycle number
        last_program: Last program that was run
        metrics: Optional dict with r_free, map_cc, etc.

    Returns:
        tuple: (should_stop: bool, reason: str or None)
    """
    if not directives:
        return False, None

    stop_cond = directives.get("stop_conditions", {})
    if not stop_cond:
        return False, None

    # Check after_cycle
    if "after_cycle" in stop_cond:
        if cycle_number >= stop_cond["after_cycle"]:
            return True, "Reached cycle %d (directive: after_cycle=%d)" % (
                cycle_number, stop_cond["after_cycle"]
            )

    # after_program — intentionally NOT a hard stop (v112.78, Bug 7)
    # ─────────────────────────────────────────────────────────────────
    # Previously, after_program caused an immediate stop when the target
    # program completed.  This broke multi-goal requests: the directive
    # extractor can only name one program, so goals beyond that program
    # were silently dropped.
    #
    # Now after_program is a *minimum-run guarantee* — the PLAN node uses
    # it to suppress auto-stop until the target program has run, but the
    # LLM decides when to actually stop.
    #
    # Hard stops that remain: after_cycle and metrics targets.

    # Check max_refine_cycles
    if "max_refine_cycles" in stop_cond and last_program == "phenix.refine":
        # Would need history to count refine cycles
        # For now, this is handled in the validator
        pass

    # Check metrics targets
    if metrics:
        if "r_free_target" in stop_cond and "r_free" in metrics:
            if metrics["r_free"] <= stop_cond["r_free_target"]:
                return True, "R-free %.3f reached target %.3f" % (
                    metrics["r_free"], stop_cond["r_free_target"]
                )

        if "map_cc_target" in stop_cond and "map_cc" in metrics:
            if metrics["map_cc"] >= stop_cond["map_cc_target"]:
                return True, "Map CC %.3f reached target %.3f" % (
                    metrics["map_cc"], stop_cond["map_cc_target"]
                )

    return False, None


def format_directives_for_display(directives):
    """
    Format directives for human-readable display.

    Args:
        directives: Directives dict

    Returns:
        str: Formatted string
    """
    if not directives:
        return "No directives extracted."

    lines = ["Extracted Directives:"]

    if "program_settings" in directives:
        _ps = directives["program_settings"]
        if isinstance(_ps, dict):
            lines.append("\n  Program Settings:")
            for prog, settings in _ps.items():
                if isinstance(settings, dict):
                    settings_str = ", ".join(
                        "%s=%s" % (k, v)
                        for k, v in settings.items())
                    lines.append(
                        "    %s: %s" % (prog, settings_str))
                else:
                    lines.append(
                        "    %s: %s" % (prog, settings))

    if "stop_conditions" in directives:
        lines.append("\n  Stop Conditions:")
        for key, value in directives["stop_conditions"].items():
            lines.append("    %s: %s" % (key, value))

    if "file_preferences" in directives:
        lines.append("\n  File Preferences:")
        for key, value in directives["file_preferences"].items():
            lines.append("    %s: %s" % (key, value))

    if "workflow_preferences" in directives:
        lines.append("\n  Workflow Preferences:")
        for key, value in directives["workflow_preferences"].items():
            lines.append("    %s: %s" % (key, value))

    if "constraints" in directives:
        lines.append("\n  Constraints:")
        for constraint in directives["constraints"]:
            lines.append("    - %s" % constraint)

    return "\n".join(lines)


# =============================================================================
# SIMPLE EXTRACTION (NO LLM FALLBACK)
# =============================================================================

# =====================================================================
# General after_program resolver (v115.10)
#
# Maps user language to workflow actions, then resolves after_program
# based on how many actions the user mentioned and whether they said
# "stop."  Replaces per-workflow overlays (ligand, denmod) with one
# unified mechanism.
# =====================================================================

# Action → (xray_program, cryoem_program, keywords)
# Keywords are matched longest-first.  Single-word keywords use word-
# boundary matching; multi-word use substring matching.
_ACTION_TABLE = {
    "analyze": {
        "xray": "phenix.xtriage",
        "cryoem": "phenix.mtriage",
        "keywords": [
            "analyze data", "data quality", "data analysis",
            "analyze map", "xtriage", "mtriage", "map quality",
            "map triage",
        ],
    },
    "solve": {
        "xray": "phenix.phaser",
        "cryoem": None,
        "keywords": [
            # v119.H14: removed "solve the structure" and
            # "solve structure" — these are goal phrases, not method
            # requests.  When a README says "Solve the structure..."
            # combined with another action (e.g., "...and stop after
            # refinement"), the multi-action branch was setting
            # start_with_program=phenix.phaser, which forced phaser
            # into the workflow even on SAD/MAD datasets where the
            # correct method is autosol.  The remaining keywords
            # still match explicit method requests ("molecular
            # replacement", "phaser", "mr "), which is the intended
            # semantic of the `solve` action.  See run_39_openai
            # batch analysis (1029B-sad and similar datasets).
            "molecular replacement", "phaser", "mr ",
        ],
    },
    "phase": {
        "xray": "phenix.autosol",
        "cryoem": None,
        "keywords": [
            "experimental phasing", "sad phasing",
            "mad phasing", "autosol", "anomalous phasing",
        ],
    },
    "predict_build": {
        "xray": "phenix.predict_and_build",
        "cryoem": "phenix.predict_and_build",
        "keywords": [
            "predict and build", "predict, build",
            "predict then build", "predict_and_build",
        ],
    },
    "predict": {
        "xray": "phenix.predict_and_build",
        "cryoem": "phenix.predict_and_build",
        "keywords": [
            "predict a model", "predict model", "predict",
        ],
    },
    "build": {
        "xray": "phenix.autobuild",
        "cryoem": "phenix.predict_and_build",
        "keywords": [
            "build a model", "build model", "build the model",
            "model building", "autobuild", "map to model",
            "map_to_model", "build",
        ],
    },
    "denmod": {
        "xray": "phenix.autobuild_denmod",
        "cryoem": "phenix.resolve_cryo_em",
        "keywords": [
            "density modif", "density-modif", "density mod",
            "denmod", "resolve_cryo_em", "resolve cryo",
            "improve map", "improve the map", "improve phases",
        ],
    },
    "sharpen": {
        "xray": None,
        "cryoem": "phenix.map_sharpening",
        "keywords": [
            "map sharpening", "auto-sharpen", "autosharpen",
            "sharpen",
        ],
    },
    "refine": {
        "xray": "phenix.refine",
        "cryoem": "phenix.real_space_refine",
        "keywords": [
            "real_space_refine", "refinement", "refine",
            "pick water", "add water",
        ],
    },
    "ligandfit": {
        "xray": "phenix.ligandfit",
        "cryoem": None,
        "keywords": [
            "ligand fit", "fit ligand", "ligandfit",
            "place ligand", "fit the ligand",
            "place the ligand", "fit atp", "fit nad",
            "fit fad", "fit heme", "add ligand",
            "add the ligand", "dock ligand",
        ],
    },
    "validate": {
        "xray": "phenix.molprobity",
        "cryoem": "phenix.molprobity",
        "keywords": [
            "validation", "validate", "molprobity",
        ],
    },
    "polder": {
        "xray": "phenix.polder",
        "cryoem": None,
        "keywords": [
            "polder map", "omit map", "polder",
        ],
    },
    "dock": {
        "xray": None,
        "cryoem": "phenix.dock_in_map",
        "keywords": [
            "dock in map", "dock into map", "dock_in_map",
        ],
    },
    "symmetry": {
        "xray": None,
        "cryoem": "phenix.map_symmetry",
        "keywords": [
            "map symmetry", "find symmetry",
            "determine symmetry", "apply symmetry",
            "apply ncs", "ncs", "oligomer",
        ],
    },
}

# Pre-compiled negation pattern — checks text immediately
# before a keyword match for negation words.
_NEGATION_BEFORE_RE = re.compile(
    r"(?:don'?t|do\s+not|never|skip|no)\s*$",
    re.IGNORECASE)


def _detect_actions(advice_lower):
    """Detect workflow action mentions in advice, ordered by
    position of first appearance.

    Uses _ACTION_TABLE keywords.  Longest keywords are checked
    first to avoid partial matches ("predict and build" before
    "predict" and "build").  Single-word keywords use word-boundary
    matching to avoid false positives ("solve" inside "resolve").
    Negated actions ("don't build", "skip refinement") are excluded.

    Compound rule: if both "predict" and "build" are detected
    separately, they merge into a single "predict_build" action.

    Returns:
        list of (action_name, position) sorted by position.
    """
    # Build (keyword, action) pairs, longest first
    _all_kw = []
    for action, info in _ACTION_TABLE.items():
        for kw in info["keywords"]:
            _all_kw.append((kw, action))
    _all_kw.sort(key=lambda x: -len(x[0]))

    # Find first non-negated occurrence of each action
    found = {}  # action → position
    for kw, action in _all_kw:
        if action in found:
            continue
        if " " not in kw and "-" not in kw:
            # Single-word: word-boundary matching
            m = re.search(
                r'\b' + re.escape(kw) + r'\b', advice_lower)
            if m:
                before = advice_lower[:m.start()]
                if _NEGATION_BEFORE_RE.search(before):
                    continue  # Negated
                found[action] = m.start()
        else:
            # Multi-word: substring matching
            pos = advice_lower.find(kw)
            if pos >= 0:
                before = advice_lower[:pos]
                if _NEGATION_BEFORE_RE.search(before):
                    continue
                found[action] = pos

    # Compound rule: predict + build → predict_build
    if "predict" in found and "build" in found:
        pos = found["predict"]
        del found["predict"]
        del found["build"]
        if "predict_build" not in found:
            found["predict_build"] = pos
    # If predict_build matched directly, remove individuals
    if "predict_build" in found:
        found.pop("predict", None)
        found.pop("build", None)

    return sorted(found.items(), key=lambda x: x[1])


# ── Stop-after detection ──────────────────────────────────────────
# v116.x: distinguish USER-explicit "stop after X" from plan-injected
# `after_program` per-stage hints.  Used as the gate for stop-analysis
# in workflow_engine._apply_directives.  See ARCHITECTURE.md
# "Stop-after directive routing" for the design rationale.

_POSITIVE_STOP_AFTER_PATTERNS = (
    re.compile(r'\bstop\s+after\b', re.IGNORECASE),
    re.compile(r'\band\s+stop\b', re.IGNORECASE),
    re.compile(r'\bthen\s+stop\b', re.IGNORECASE),
    re.compile(r',\s*stop\b', re.IGNORECASE),
    re.compile(r'\.\s+stop\b', re.IGNORECASE),
    re.compile(r'\bstop\s+when\b', re.IGNORECASE),
    re.compile(r'\bstop\s+once\b', re.IGNORECASE),
    re.compile(r'\bstop\s+if\b', re.IGNORECASE),
    re.compile(r'\bstop\s+at\b', re.IGNORECASE),
    re.compile(r'\bonly\s+run\b', re.IGNORECASE),
    re.compile(r'\bjust\s+(?:do|run)\b', re.IGNORECASE),
    # v117.3: contiguous-phrase patterns for phrasings that
    # _IMPERATIVE_STOP_MARKERS catches near the program name but
    # that ALSO deserve global recognition (so _is_stop_after_requested
    # returns True and the v117.1 grounding-flag exemption fires
    # even when the LLM doesn't set stop_after_requested itself).
    # Both patterns are contiguous multi-word phrases that cannot
    # span sentence boundaries by construction.  Over-permissive
    # "after X completes/finishes" patterns were considered and
    # rejected during v117.3 plan review — see plan §3 Change 2.
    re.compile(r'\bstop\s+the\s+workflow\b', re.IGNORECASE),
    re.compile(r'\bis\s+the\s+last\s+step\b', re.IGNORECASE),
)

_NEGATIVE_STOP_PATTERN = re.compile(
    r"(?:don'?t|do\s+not|never)\s+stop", re.IGNORECASE)

# "Stop Condition: None" / "Stop Condition: not specified" / heading-
# only — explicit NO-stop-condition signal in README/preprocessed
# advice.  Must be stripped before pattern matching so the bare word
# "stop" inside this phrase doesn't trigger a false positive.
_STOP_CONDITION_NONE = re.compile(
    r'\*?\*?stop\s+conditions?\*?\*?\s*[:=]\s*'
    r'(?:none|not\s+specified|n/?a|null|\s*$)',
    re.IGNORECASE)

# "Stop Condition: <real value>" — explicit YES-stop-condition signal.
_STOP_CONDITION_VALUE = re.compile(
    r'\*?\*?stop\s+conditions?\*?\*?\s*[:=]\s*'
    r'(?!none\b|not\s+specified\b|n/?a\b|null\b|\s*$)\S',
    re.IGNORECASE)


def _is_stop_after_requested(advice):
    """Return True iff the user explicitly requested a stop-after condition.

    True iff the raw user advice (or README content) contains an
    explicit stop-after phrasing.  False for the absence of any stop
    signal, for explicit "Stop Condition: None", and for explicit
    negations like "do not stop".

    The bare word "stop" is not sufficient; it must match one of the
    recognised positive patterns.  This is distinct from plan-injected
    `after_program` hints — those do NOT invoke this function and
    therefore stop_after_requested stays False under plan progression.

    Args:
        advice: raw user advice string (may include README content,
                preprocessed sections, etc.).

    Returns:
        True if the advice contains an explicit stop-after directive.
    """
    if not advice:
        return False

    # Explicit negation overrides everything else.
    if _NEGATIVE_STOP_PATTERN.search(advice):
        return False

    # Strip "Stop Condition: None" patterns so they don't trigger
    # bare-"stop" false positives downstream.
    text = _STOP_CONDITION_NONE.sub(' ', advice)

    # Positive: "Stop Condition: <real value>".
    if _STOP_CONDITION_VALUE.search(text):
        return True

    # Positive: recognised stop-after phrasings.
    for pat in _POSITIVE_STOP_AFTER_PATTERNS:
        if pat.search(text):
            return True

    return False


def _resolve_after_program(directives, advice_lower):
    """General after_program resolver.

    Corrects after_program based on how many workflow actions the
    user mentioned and whether they said "stop."

    Rules:
      Multiple actions + stop  → after_program = last mentioned
      Multiple actions + no stop → clear after_program (plan drives)
      Single action + stop     → set after_program if not set
      Otherwise                → leave as-is

    Also sets start_with_program to the first action when multiple
    actions are detected (helps the workflow engine prioritize).

    Called from _apply_workflow_intent_fallback (both LLM and
    rules-only paths).
    """
    # If this is preprocessed text (has "Primary Goal:" header),
    # extract just that line for action detection.  The full
    # preprocessed text contains sections like
    # "Special Instructions: Do not proceed to refinement"
    # that cause false positive action detection — the negation
    # "Do not" is too far from "refinement" for the immediate-
    # predecessor check to catch.
    _action_text = advice_lower
    _goal_match = re.search(
        r'primary\s+goal\s*:\s*(.+?)(?:\n|$)',
        advice_lower)
    if _goal_match:
        _action_text = _goal_match.group(1).strip()

    actions = _detect_actions(_action_text)
    if not actions:
        return  # No actions detected — leave directives as-is

    # Detect stop intent.  Use the dedicated helper that handles
    # "Stop Condition: None" / "not specified" correctly — the
    # bare \bstop\b check used previously had a false positive on
    # advice containing "Stop Condition: None" (e.g. Tom's
    # nsf-d2-ligand case, where the user explicitly said "no stop"
    # but the regex matched the word "Stop" in the heading).
    _has_stop = _is_stop_after_requested(advice_lower)

    # Infer experiment type from advice text (same heuristic
    # as extract_directives_simple).
    _is_cryoem = bool(re.search(
        r'cryo-?em|half.?map|\.mrc|\.map|full.?map|mtriage'
        r'|sharpen|map.sharpening|auto.sharpen'
        r'|resolve.cryo|dock.in.map|map.to.model'
        r'|density.modif\w*\s+map',
        advice_lower))
    _is_xray = bool(re.search(
        r'x-?ray|\.mtz|\.hkl|xtriage|phaser|autosol|'
        r'sad|mad|molecular.?replacement',
        advice_lower))
    if _is_cryoem and not _is_xray:
        _exp = "cryoem"
    else:
        _exp = "xray"  # Default to xray

    # Detect preprocessed advice — used below to suppress
    # start_with_program writes when the multi-action signal
    # comes from descriptive preprocessor prose rather than
    # imperative user advice.
    #
    # v116.11: AF_7mjs regression — the preprocessor produces
    # Primary Goal text like "Run PredictAndBuild ... rebuild
    # the loop and refine the model" describing a workflow.
    # _detect_actions finds multiple actions (predict_build,
    # refine), and pre-fix the resolver set start_with_program=
    # phenix.predict_and_build.  The planner then treats that
    # as a skip-prerequisites directive and skips Stage 1
    # (mtriage) + Stage 2 (denmod), giving a 3-stage plan
    # instead of the correct 5-stage plan.  When the advice
    # is preprocessed, descriptive multi-action prose is the
    # norm, not user prescription; don't infer start_with.
    _PREPROC_SIGS = (
        r'input\s+files?\s+found',
        r'experiment\s+type',
        r'key\s+parameters?',
        r'special\s+instructions?',
    )
    _is_preprocessed = any(
        re.search(r'(?i)^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*%s\**\s*:'
                  % sig, advice_lower, re.MULTILINE)
        for sig in _PREPROC_SIGS
    )

    _sc = directives.get("stop_conditions", {})
    n = len(actions)

    if n > 1:
        if _has_stop:
            # "do A, B and stop" → after_program = B (last)
            last_action = actions[-1][0]
            last_prog = _ACTION_TABLE[last_action].get(
                _exp) or _ACTION_TABLE[last_action].get(
                "xray")
            if last_prog:
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                directives["stop_conditions"][
                    "after_program"] = last_prog
                directives["stop_conditions"][
                    "skip_validation"] = True
                # v116.x: flag user-explicit stop so the consumer
                # (workflow_engine._apply_directives) knows to apply
                # the after_program_done stop analysis.  Plan-injected
                # after_program does NOT set this flag, which is what
                # lets multi-stage plans progress without being cut
                # short at each stage boundary.
                directives["stop_conditions"][
                    "stop_after_requested"] = True
        else:
            # "do A and B" → clear, plan drives
            _sc.pop("after_program", None)
            _sc.pop("skip_validation", None)
            _sc.pop("stop_after_requested", None)
        # Set start_with_program to first action.
        #
        # v116.11: Skip this for preprocessed advice — the
        # multi-action signal there comes from descriptive
        # prose (Primary Goal summary), not user prescription.
        # Setting start_with would cause the planner to skip
        # prerequisite stages (AF_7mjs regression).  Real user
        # prose like "run phaser and refine" still gets
        # start_with=phaser (preprocessed=False).
        if not _is_preprocessed:
            first_action = actions[0][0]
            first_prog = _ACTION_TABLE[first_action].get(
                _exp) or _ACTION_TABLE[first_action].get(
                "xray")
            if first_prog:
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                directives["stop_conditions"][
                    "start_with_program"] = first_prog
    elif n == 1 and _has_stop:
        # "do A and stop" → set after_program to the detected
        # action's program.  Always override — the LLM often
        # sets the wrong program (e.g., real_space_refine instead
        # of resolve_cryo_em for "density modify and stop").
        action_name = actions[0][0]
        prog = _ACTION_TABLE[action_name].get(
            _exp) or _ACTION_TABLE[action_name].get(
            "xray")
        if prog:
            if "stop_conditions" not in directives:
                directives["stop_conditions"] = {}
            directives["stop_conditions"][
                "after_program"] = prog
            directives["stop_conditions"][
                "skip_validation"] = True
            # v116.x: see comment above on stop_after_requested.
            directives["stop_conditions"][
                "stop_after_requested"] = True
    # n == 1, no stop → leave as-is
    # n == 0 → leave as-is

    # Special: "predict" action (not predict_build) implies
    # stop_after_predict=True — the user wants prediction only,
    # not the full predict-and-build workflow.
    # Nested under program name so command_builder merges it
    # into strategy (line 1645: prog_settings[program]).
    # NOTE: must re-read stop_conditions after resolver modified it.
    _final_sc = directives.get("stop_conditions", {})
    if (n >= 1 and actions[-1][0] == "predict"
            and _final_sc.get("after_program")):
        if "program_settings" not in directives:
            directives["program_settings"] = {}
        _pab = "phenix.predict_and_build"
        if _pab not in directives["program_settings"]:
            directives["program_settings"][_pab] = {}
        directives["program_settings"][_pab][
            "stop_after_predict"] = True


def _apply_workflow_intent_fallback(directives, advice_lower):
    """Apply deterministic workflow intent patterns to directives.

    v115.09: Shared by both LLM and rules-based paths.  When the LLM
    extracts directives, it may miss ``wants_validation_only`` or
    ``use_mr_sad`` even when the advice clearly requests them.  This
    function overlays the deterministic pattern matches on top of
    whatever the LLM returned, ensuring critical routing flags are
    always set when the advice text contains unambiguous signals.

    Called from:
      - ``extract_directives`` (post-LLM overlay)
      - ``extract_directives_simple`` (rules-only path)

    Args:
        directives: dict — modified in place.
        advice_lower: str — lowercased user advice text.
    """
    # Validation-only intent.
    # NOTE: "analysis only" deliberately omitted — it matches
    # data-quality analysis (xtriage-only) tutorials, not just
    # model validation.  The LLM prompt handles that nuance.
    _validation_signals = [
        "model validation",
        "structure validation", "comprehensive validation",
        "run molprobity", "validate this structure",
        "validation and correction",
    ]
    if any(sig in advice_lower for sig in _validation_signals):
        if "workflow_preferences" not in directives:
            directives["workflow_preferences"] = {}
        directives["workflow_preferences"][
            "wants_validation_only"] = True

    # MR-SAD intent.
    _mr_sad_patterns = [
        "mr-sad", "mr sad", "mrsad",
        "molecular replacement sad",
        "molecular replacement followed by sad",
        "mr followed by sad",
    ]
    if any(pat in advice_lower for pat in _mr_sad_patterns):
        if "workflow_preferences" not in directives:
            directives["workflow_preferences"] = {}
        directives["workflow_preferences"]["use_mr_sad"] = True
        directives["workflow_preferences"][
            "use_experimental_phasing"] = True

    # Ligand-fitting implies model is placed.
    # You can't fit a ligand into an unplaced model.  When the user
    # says "fit ATP" or "fit ligand", the model must already be in
    # the correct position — no MR or docking needed.
    # Guard: only set model_is_placed when the advice does NOT also
    # mention molecular replacement, phaser, solving, or docking
    # (those indicate the model still needs placement first).
    _ligand_fit_signals = [
        "fit ligand", "fit the ligand", "ligandfit",
        "fit atp", "fit nad", "fit fad", "fit heme",
        "place ligand", "place the ligand",
        "add ligand", "add the ligand",
        "dock ligand", "dock the ligand",
    ]
    _mr_signals = [
        "molecular replacement", "phaser", "solve",
        "mr ", "autosol", "predict",
        "dock in map", "dock_in_map", "dock into",
    ]
    if (any(sig in advice_lower for sig in _ligand_fit_signals) and
            not any(sig in advice_lower for sig in _mr_signals)):
        if "workflow_preferences" not in directives:
            directives["workflow_preferences"] = {}
        directives["workflow_preferences"][
            "model_is_placed"] = True

    # ── General after_program resolver (v115.10) ──────────
    # Replaces per-workflow overlays (ligand after_program
    # clearing, denmod stop/clear) with a single mechanism
    # that counts distinct action mentions and checks for
    # "stop" to determine the correct after_program.
    _resolve_after_program(directives, advice_lower)


def extract_directives_simple(user_advice, log=None):
    """
    Extract directives using simple pattern matching (no LLM).

    This is a fallback when LLM is not available or fails.
    Only handles common, unambiguous patterns.

    Args:
        user_advice: User advice string
        log: optional log function (e.g., ``print`` or a custom
            logger).  Forwarded to ``validate_directives`` at the
            end so sanity-check drops are visible to operators.
            v119.H14.1: defaults to None (silent) for backward
            compatibility with the four external callers in
            ai_agent.py and run_ai_analysis.py that didn't pass
            a logger pre-H14.1.

    Returns:
        dict: Extracted directives (limited), passed through
        ``validate_directives`` for sanity-check normalization.
    """
    if not user_advice:
        return {}

    # Strip LLM-fabricated "Stop Condition:" lines (Fix 2, v115).
    user_advice = _strip_preprocessor_stop_condition(user_advice)

    # Classify intent (v115).
    # The intent controls stopping behavior downstream.
    intent_result = None
    try:
        try:
            from libtbx.langchain.agent.intent_classifier \
                import classify_intent
        except ImportError:
            from agent.intent_classifier \
                import classify_intent
        intent_result = classify_intent(user_advice)
    except Exception:
        pass  # intent_classifier not available

    directives = {}

    # Store intent in directives for downstream use.
    # Low-confidence means the default fallback fired (no pattern
    # matched) — not a genuine detection, so don't pollute directives.
    if intent_result and intent_result.get("confidence") != "low":
        directives["intent"] = intent_result

    advice_lower = user_advice.lower()

    # Extract resolution
    res_patterns = [
        r'resolution\s*[=:]\s*([0-9.]+)',
        r'resolution\s+(?:of\s+)?([0-9.]+)',  # "resolution 3.0" or "resolution of 3.0"
        r'resolution\s+limit\s*[=:]?\s*([0-9.]+)',
        r'([0-9.]+)\s*(?:angstrom|Å|A)\s*resolution',
    ]
    # Note: removed overly broad patterns like "to X Å" and "X Å" standalone
    # which match wavelengths and other Angstrom values

    # Also extract wavelength so we can check for confusion
    wavelength = None
    wl_match = re.search(r'wavelength\s*[=:]\s*([0-9.]+)', user_advice, re.IGNORECASE)
    if wl_match:
        try:
            wavelength = float(wl_match.group(1))
        except ValueError:
            pass

    for pattern in res_patterns:
        match = re.search(pattern, user_advice, re.IGNORECASE)
        if match:
            try:
                res = float(match.group(1))
                if 1.0 <= res <= 20:
                    # Skip if it matches the wavelength (likely confused)
                    if wavelength and abs(res - wavelength) < 0.01:
                        continue
                    if "program_settings" not in directives:
                        directives["program_settings"] = {}
                    directives["program_settings"]["default"] = {"resolution": res}
                    break
            except ValueError:
                pass

    # Check for anisotropic
    if "anisotropic" in advice_lower:
        if "program_settings" not in directives:
            directives["program_settings"] = {}
        if "phenix.refine" not in directives["program_settings"]:
            directives["program_settings"]["phenix.refine"] = {}
        directives["program_settings"]["phenix.refine"]["anisotropic_adp"] = True

    # ==========================================================================
    # UNIT CELL EXTRACTION
    # ==========================================================================
    # Match forms like:
    #   "unit cell (116.097, 116.097, 44.175, 90, 90, 120)"
    #   "unit_cell = 116 116 44 90 90 120"
    #   "the specified unit cell (a, b, c, α, β, γ) must be used"
    #   "unit cell: 116.097 116.097 44.175 90 90 120"
    _NUM = r'[-+]?\d+(?:\.\d+)?'
    _SEP = r'[\s,]+'
    _uc_pattern = (
        r'unit[_ ]cell'
        r'(?:\s*[=:(\s]\s*|\s+(?:is|of|=|:)?\s*[(]?)'    # separator
        r'(' + _NUM + _SEP + _NUM + _SEP + _NUM + _SEP
          + _NUM + _SEP + _NUM + _SEP + _NUM + r')'       # 6 numbers
    )
    uc_match = re.search(_uc_pattern, user_advice, re.IGNORECASE)
    if uc_match:
        # Normalise to 6 space-separated numbers (strip parens/commas)
        raw = uc_match.group(1)
        nums = re.findall(r'[-+]?\d+(?:\.\d+)?', raw)
        if len(nums) == 6:
            uc_str = " ".join(nums)
            if "program_settings" not in directives:
                directives["program_settings"] = {}
            if "default" not in directives["program_settings"]:
                directives["program_settings"]["default"] = {}
            directives["program_settings"]["default"]["unit_cell"] = uc_str

    # ==========================================================================
    # SPACE GROUP EXTRACTION
    # ==========================================================================
    # Match: "space group P 32 2 1", "space_group=P63", "space group: C 2 2 21"
    _sg_pattern = (
        r'space[_ ]group'
        r'(?:\s*[=:]\s*|\s+(?:is|of|=|:)?\s*)'
        r'([A-Za-z][A-Za-z0-9 /_-]{1,20})'  # Herman–Mauguin symbol
    )
    sg_match = re.search(_sg_pattern, user_advice, re.IGNORECASE)
    if sg_match:
        sg_str = sg_match.group(1).strip().rstrip('.')
        if sg_str:
            if "program_settings" not in directives:
                directives["program_settings"] = {}
            if "default" not in directives["program_settings"]:
                directives["program_settings"]["default"] = {}
            directives["program_settings"]["default"]["space_group"] = sg_str

    # ==========================================================================
    # ASU COPY COUNT EXTRACTION
    # ==========================================================================
    # Match forms like:
    #   "4 copies of the search model"        "there are 4 copies in the ASU"
    #   "search for 4 copies"                 "4 copies in the asymmetric unit"
    #   "copies=4"  "copies: 4"               "find 4 copies"
    #   "place 4 copies"                      "4 molecules in the ASU"
    # Guard: only accept integers 1-30 (same sanity bound as xtriage path).
    _copies_patterns = [
        # Explicit number before "copies" or "molecules" with context keywords
        r'(\d+)\s+cop(?:y|ies)\s+(?:of|in)',          # "4 copies of/in"
        r'(\d+)\s+molecules?\s+in\s+(?:the\s+)?asu',  # "4 molecules in the ASU"
        r'(?:search|find|place|look)\s+(?:for\s+)?(\d+)\s+cop',  # "search for 4 copies"
        r'cop(?:y|ies)\s*[=:]\s*(\d+)',               # "copies=4" / "copies: 4"
        r'(?:there\s+are|has?)\s+(\d+)\s+cop',        # "there are 4 copies"
        r'asu\s+(?:contains?|has?)\s+(\d+)\s+cop',    # "ASU contains 4 copies"
    ]
    for _cp_pat in _copies_patterns:
        _cp_m = re.search(_cp_pat, user_advice, re.IGNORECASE)
        if _cp_m:
            try:
                _cp_v = int(_cp_m.group(1))
                if 1 <= _cp_v <= 30:
                    if "program_settings" not in directives:
                        directives["program_settings"] = {}
                    if "default" not in directives["program_settings"]:
                        directives["program_settings"]["default"] = {}
                    directives["program_settings"]["default"]["copies"] = _cp_v
                    break
            except (ValueError, IndexError):
                pass

    # Check for macro_cycles / number of refinement cycles
    cycle_patterns = [
        r'(\d+)\s*(?:macro[- ]?)?cycle',  # "1 macro cycle", "3 cycles"
        r'(?:macro[- ]?)?cycles?\s*[=:]\s*(\d+)',  # "macro_cycles=1", "cycles: 3"
        r'number[_ ]of[_ ]macro[_ ]cycles\s*[=:]\s*(\d+)',  # "number_of_macro_cycles=1"
    ]
    for pattern in cycle_patterns:
        match = re.search(pattern, advice_lower)
        if match:
            try:
                n_cycles = int(match.group(1))
                if 1 <= n_cycles <= 50:
                    if "program_settings" not in directives:
                        directives["program_settings"] = {}
                    # Apply to both refine and real_space_refine
                    for prog in ("phenix.refine", "phenix.real_space_refine"):
                        if prog not in directives["program_settings"]:
                            directives["program_settings"][prog] = {}
                        directives["program_settings"][prog]["macro_cycles"] = n_cycles
                    break
            except ValueError:
                pass

    # ==========================================================================
    # ATOM TYPE EXTRACTION (for autosol)
    # ==========================================================================
    atom_type_patterns = [
        # "use selenium" or "selenium anomalous scatterer"
        (r'\bselenium\b|se[- ]?sad|se[- ]?met|\bse\s+as\b', 'Se'),
        # "use sulfur" or "sulfur SAD"
        (r'\bsulfur\b|s[- ]?sad|\bs\s+as\b', 'S'),
        # "use bromine" or "Br SAD"
        (r'\bbromine\b|br[- ]?sad|\bbr\s+as\b', 'Br'),
        # "use iodine"
        (r'\biodine\b|i[- ]?sad|\bi\s+as\b', 'I'),
        # "use zinc"
        (r'\bzinc\b|zn[- ]?sad|\bzn\s+as\b', 'Zn'),
        # "use iron"
        (r'\biron\b|fe[- ]?sad|\bfe\s+as\b', 'Fe'),
    ]

    for pattern, atom_type in atom_type_patterns:
        if re.search(pattern, advice_lower, re.IGNORECASE):
            if "program_settings" not in directives:
                directives["program_settings"] = {}
            if "phenix.autosol" not in directives["program_settings"]:
                directives["program_settings"]["phenix.autosol"] = {}
            directives["program_settings"]["phenix.autosol"]["atom_type"] = atom_type
            break

    # ==========================================================================
    # FILE PREFERENCES EXTRACTION
    # ==========================================================================
    # Check for anomalous data preference
    if re.search(r'(?:use|prefer|keep)\s+(?:the\s+)?anomalous|anomalous\s+(?:data|signal)', advice_lower):
        if "file_preferences" not in directives:
            directives["file_preferences"] = {}
        directives["file_preferences"]["prefer_anomalous"] = True

    # Check for unmerged data preference
    if re.search(r'(?:use|prefer|keep)\s+(?:the\s+)?unmerged|unmerged\s+data', advice_lower):
        if "file_preferences" not in directives:
            directives["file_preferences"] = {}
        directives["file_preferences"]["prefer_unmerged"] = True

    # ==========================================================================
    # WORKFLOW PREFERENCES EXTRACTION
    # ==========================================================================
    # Check for prefer program patterns (include, use, run)
    prefer_program_patterns = [
        (r'(?:include|use|run|do|perform|add)\s+(?:the\s+)?ligand\s*(?:fit|fitting)', 'phenix.ligandfit'),
        (r'(?:fit|place|add)\s+(?:the\s+)?ligand', 'phenix.ligandfit'),
        (r'with\s+ligand\s*(?:fit|fitting)', 'phenix.ligandfit'),
        (r'ligand\s+(?:should\s+be\s+)?(?:fit|fitted|placed)', 'phenix.ligandfit'),
        (r'(?:include|use|run|do|perform)\s+(?:the\s+)?polder', 'phenix.polder'),
        (r'(?:calculate|compute|generate)\s+polder\s+map', 'phenix.polder'),
    ]

    for pattern, program in prefer_program_patterns:
        if re.search(pattern, advice_lower, re.IGNORECASE):
            if "workflow_preferences" not in directives:
                directives["workflow_preferences"] = {"prefer_programs": []}
            if "prefer_programs" not in directives["workflow_preferences"]:
                directives["workflow_preferences"]["prefer_programs"] = []
            if program not in directives["workflow_preferences"]["prefer_programs"]:
                directives["workflow_preferences"]["prefer_programs"].append(program)

    # ==========================================================================
    # POLDER SELECTION EXTRACTION
    # ==========================================================================
    # Extract atom selection for polder from patterns like:
    #   "polder map for chain A resseq 88"
    #   "polder for chain A and resseq 88"
    #   "polder selection chain A and resseq 88"
    if "polder" in advice_lower:
        _polder_sel = _extract_polder_selection(user_advice)
        if _polder_sel:
            if "program_settings" not in directives:
                directives["program_settings"] = {}
            if "phenix.polder" not in directives["program_settings"]:
                directives["program_settings"]["phenix.polder"] = {}
            directives["program_settings"]["phenix.polder"]["selection"] = _polder_sel

    # Check for skip program patterns
    skip_program_patterns = [
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?autobuild', 'phenix.autobuild'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?ligand\s*fit', 'phenix.ligandfit'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?molprobity', 'phenix.molprobity'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?phaser', 'phenix.phaser'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?mtriage', 'phenix.mtriage'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?xtriage', 'phenix.xtriage'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?(?:real.?space.?)?refine', 'phenix.real_space_refine'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?validation', 'phenix.molprobity'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?map.?sharpening', 'phenix.map_sharpening'),
        (r'(?:skip|no|avoid|(?:don\'?t|do\s+not)\s+(?:run|use))\s+(?:the\s+)?resolve', 'phenix.resolve_cryo_em'),
    ]

    for pattern, program in skip_program_patterns:
        if re.search(pattern, advice_lower, re.IGNORECASE):
            if "workflow_preferences" not in directives:
                directives["workflow_preferences"] = {"skip_programs": []}
            if "skip_programs" not in directives["workflow_preferences"]:
                directives["workflow_preferences"]["skip_programs"] = []
            if program not in directives["workflow_preferences"]["skip_programs"]:
                directives["workflow_preferences"]["skip_programs"].append(program)

    # Check for stop conditions
    if "stop" in advice_lower or "skip" in advice_lower or "don't" in advice_lower or "no valid" in advice_lower:
        if "stop_conditions" not in directives:
            directives["stop_conditions"] = {}

        # "stop after cycle N"
        cycle_match = re.search(r'stop\s+(?:after|at)\s+cycle\s+(\d+)', advice_lower)
        if cycle_match:
            directives["stop_conditions"]["after_cycle"] = int(cycle_match.group(1))

        # "stop after first/one refinement" - also accept "stop after refinement"
        if re.search(r'stop\s+after\s+(?:the\s+)?(?:first|one|1|a\s+single)?\s*refine', advice_lower):
            directives["stop_conditions"]["after_program"] = "phenix.refine"
            # v116.x: this regex matches "stop after" — an explicit user stop.
            directives["stop_conditions"]["stop_after_requested"] = True
            # Only set max_refine_cycles=1 if explicitly "first" or "one"
            if re.search(r'(?:first|one|1|single)\s*refine', advice_lower):
                directives["stop_conditions"]["max_refine_cycles"] = 1

        # "skip validation" or "don't validate" or "no validation"
        if re.search(r"(?:skip|no|don'?t)\s+validat", advice_lower):
            directives["stop_conditions"]["skip_validation"] = True

        # Clean up empty stop_conditions
        if not directives["stop_conditions"]:
            del directives["stop_conditions"]

    # Detect experiment type from context to help with density modification
    is_cryoem = bool(re.search(r'cryo-?em|half.?map|\.mrc|\.map|full.?map|mtriage', advice_lower))
    is_xray = bool(re.search(r'x-?ray|\.mtz|\.hkl|xtriage|phaser|autosol|sad|mad|molecular.?replacement', advice_lower))

    # NOTE (Fix 2, v115): Removed "Stop Condition:" parser that used to
    # live here.  That block parsed LLM-fabricated stop conditions from
    # the advice preprocessor (e.g., "Stop Condition: Stop after
    # molecular replacement").  None of the original READMEs contain
    # the word "stop"; every such line was invented by the preprocessor.
    # Real user stop commands are handled by the explicit patterns above
    # (lines ~2025-2039 and ~2200-2260).  The preprocessor's "Stop
    # Condition:" lines are now stripped by
    # _strip_preprocessor_stop_condition() at entry.

    # v115.10: continuation_indicators, downstream_tasks, and
    # multi_program_patterns removed — replaced by the general
    # after_program resolver in _apply_workflow_intent_fallback
    # which counts distinct action mentions via _ACTION_TABLE.

    # Detect tutorial/procedure patterns that imply stop after
    # specific program.  These provide initial after_program
    # extraction; the general resolver in
    # _apply_workflow_intent_fallback overrides when wrong
    # (e.g., multi-step detected).
    #
    # GUARD (Fix 2, v115): Skip this section when the advice
    # is preprocessor output.  Preprocessor text contains
    # descriptive phrases like "Density modification of
    # cryo-EM map" that match tutorial_patterns but are NOT
    # user commands.  Preprocessor output is identified by
    # signature headers (Input Files Found, Experiment Type,
    # Key Parameters, Special Instructions).
    #
    # v116.11: Updated regex to handle markdown bold + numbered
    # list prefix, matching the strengthening in
    # _strip_preprocessor_stop_condition.  Pre-fix, this regex
    # required headers to start with the literal signature word
    # (no list marker, no markdown), so AF_7mjs-style preprocessor
    # output with "1. **Input Files Found**:" was NOT detected as
    # preprocessed, causing tutorial_patterns to run on
    # descriptive prose.
    _PREPROC_SIGS = (
        r'input\s+files?\s+found',
        r'experiment\s+type',
        r'key\s+parameters?',
        r'special\s+instructions?',
    )
    _is_preprocessed = any(
        re.search(r'(?i)^[ \t]*(?:[\d]+\.\s*|[-*]\s*)?\**\s*%s\**\s*:'
                  % sig, user_advice, re.MULTILINE)
        for sig in _PREPROC_SIGS
    )

    tutorial_patterns = [
        # Xtriage/twinning analysis
        (r'(?:run|check|analyze).*(?:xtriage|twinning)', 'phenix.xtriage'),
        (r'(?:check|look)\s+for\s+twin', 'phenix.xtriage'),
        (r'analyze.*(?:data\s+quality|reflection)', 'phenix.xtriage'),
        # Phaser/MR testing
        (r'(?:run|try|test).*(?:phaser|molecular\s+replacement|MR)', 'phenix.phaser'),
        (r'(?:test|check).*MR\s+solution', 'phenix.phaser'),
        # Mtriage analysis
        (r'(?:run|check).*mtriage', 'phenix.mtriage'),
        (r'analyze.*map\s+quality', 'phenix.mtriage'),
        # Map symmetry analysis
        (r'(?:run|check|find|determine|detect).*(?:map\s*symmetry|symmetry)', 'phenix.map_symmetry'),
        (r'(?:symmetry|ncs).*(?:map|cryo)', 'phenix.map_symmetry'),
        # Map to model (automated model building into cryo-EM maps)
        # Check BEFORE dock_in_map patterns
        (r'map\s*to\s*model', 'phenix.map_to_model'),
        (r'maptomodel', 'phenix.map_to_model'),
        (r'(?:use|run|try).*map\s*to\s*model', 'phenix.map_to_model'),
        (r'(?:automated?\s+)?model\s+building.*(?:into|in)\s*(?:the\s*)?(?:cryo[- ]?em\s+)?(?:density\s+)?map', 'phenix.map_to_model'),
        (r'build.*model.*(?:into|in)\s*(?:the\s*)?(?:density\s+)?map', 'phenix.map_to_model'),
        # Map sharpening
        (r'map\s*sharpening', 'phenix.map_sharpening'),
        (r'sharpen.*(?:the\s+)?map', 'phenix.map_sharpening'),
        (r'(?:run|use|try).*map\s*sharpening', 'phenix.map_sharpening'),
        # Polder omit maps
        (r'polder', 'phenix.polder'),
        (r'polder\s*map', 'phenix.polder'),
        (r'omit\s*map', 'phenix.polder'),
        (r'(?:evaluate|check|verify).*ligand.*(?:placement|density)', 'phenix.polder'),
        (r'(?:evaluate|check|verify).*(?:placement|density).*ligand', 'phenix.polder'),
        # Dock in map (rigid body fitting)
        (r'dock.*(?:in|into)\s+map', 'phenix.dock_in_map'),
        (r'(?:rigid\s+body\s+)?fit.*model.*(?:to|into)\s+map', 'phenix.dock_in_map'),
        # Ligand fitting
        (r'(?:run|try).*ligand\s*fit', 'phenix.ligandfit'),
        (r'(?:fit|place|add).*ligand', 'phenix.ligandfit'),
        (r'ligand\s*fitting', 'phenix.ligandfit'),
    ]

    # Density modification patterns - need to determine experiment type
    denmod_patterns = [
        r'density\s+modif',
        r'denmod',
        r'(?:improve).*(?:map|phases)',  # Removed 'sharpen' - now goes to map_sharpening
        r'one\s+cycle.*(?:density|denmod|map\s+improv)',
    ]

    _allow_patterns = (
        not _is_preprocessed
        or (intent_result
            and intent_result.get("intent")
            == "tutorial"))

    if _allow_patterns:
        for pattern, program in tutorial_patterns:
            if re.search(pattern, advice_lower,
                         re.IGNORECASE):
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                if "after_program" not in directives.get(
                    "stop_conditions", {}):
                    directives["stop_conditions"][
                        "after_program"] = program
                    directives["stop_conditions"][
                        "skip_validation"] = True
                    # v116.x: tutorial pattern match = explicit
                    # stop intent.  Tentative — intent classifier
                    # processing below will pop it along with
                    # after_program if the intent classifier
                    # disagrees.
                    directives["stop_conditions"][
                        "stop_after_requested"] = True
                    # Internal flag: this was inferred
                    # from patterns, not user command
                    directives["stop_conditions"][
                        "_set_by_pattern"] = True
                break

    # Handle density modification separately
    if _allow_patterns and "after_program" not in directives.get("stop_conditions", {}):
        for pattern in denmod_patterns:
            if re.search(pattern, advice_lower, re.IGNORECASE):
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                # Determine program based on experiment type
                if is_cryoem and not is_xray:
                    directives["stop_conditions"]["after_program"] = "phenix.resolve_cryo_em"
                elif is_xray and not is_cryoem:
                    directives["stop_conditions"]["after_program"] = "phenix.autobuild_denmod"
                else:
                    # Ambiguous or unknown - default to cryo-EM as it's more commonly requested
                    # The LLM extraction will do better with context
                    directives["stop_conditions"]["after_program"] = "phenix.resolve_cryo_em"
                directives["stop_conditions"]["skip_validation"] = True
                # v116.x: denmod pattern match = explicit "density
                # modify and stop" intent.  Tentative — intent
                # classifier may override.
                directives["stop_conditions"]["stop_after_requested"] = True
                directives["stop_conditions"]["_set_by_pattern"] = True
                break

    # ── Intent-driven stop behavior (v115) ──────────────────────
    # Override or set stop conditions based on classified intent.
    # Runs AFTER all pattern extraction so it can override.
    _intent = (
        directives.get("intent") or {}).get("intent")

    if _intent == "task":
        # Task: single program, stop after.
        _task_prog = (
            directives.get("intent", {})
            .get("task_program"))
        if _task_prog:
            if "stop_conditions" not in directives:
                directives["stop_conditions"] = {}
            directives["stop_conditions"][
                "after_program"] = _task_prog
            directives["stop_conditions"][
                "skip_validation"] = True
            # v116.x: intent="task" = user explicitly wants a
            # single program then stop.
            directives["stop_conditions"][
                "stop_after_requested"] = True

    elif _intent == "solve":
        # Solve: no artificial stops.  Remove any
        # after_program set by tutorial patterns.
        sc = directives.get("stop_conditions", {})
        if "after_program" in sc and sc.get(
                "_set_by_pattern"):
            del sc["after_program"]
            sc.pop("skip_validation", None)
            sc.pop("stop_after_requested", None)
            sc.pop("_set_by_pattern", None)
            if not sc:
                directives.pop("stop_conditions", None)

    elif _intent == "solve_constrained":
        # Solve-constrained: no artificial stops, but
        # set the method preference.
        sc = directives.get("stop_conditions", {})
        if "after_program" in sc and sc.get(
                "_set_by_pattern"):
            del sc["after_program"]
            sc.pop("skip_validation", None)
            sc.pop("stop_after_requested", None)
            sc.pop("_set_by_pattern", None)
            if not sc:
                directives.pop("stop_conditions", None)
        _method = (
            directives.get("intent", {})
            .get("method_constraint"))
        if _method:
            if "workflow_preferences" not in directives:
                directives["workflow_preferences"] = {}
            directives["workflow_preferences"][
                "method_constraint"] = _method

    # intent=tutorial: keep whatever patterns set
    # (after_program from tutorial_patterns is correct)

    # v115.09 Fix 3+4: Detect validation-only and MR-SAD intent
    # Calls shared helper (also used as post-LLM overlay).
    _apply_workflow_intent_fallback(directives, advice_lower)

    # Extract unit_cell and space_group if mentioned
    _extract_crystal_symmetry_simple(user_advice, directives)

    # Strip internal flags before returning
    sc = directives.get("stop_conditions", {})
    sc.pop("_set_by_pattern", None)
    if "stop_conditions" in directives and not sc:
        del directives["stop_conditions"]

    # Simplify intent: keep only the classification string,
    # not the full internal dict with reason/confidence/etc.
    _intent_dict = directives.get("intent")
    if isinstance(_intent_dict, dict):
        directives["intent"] = _intent_dict.get("intent")

    # v119.H14.1: apply validate_directives as a final sanity pass.
    #
    # Pre-H14.1, extract_directives_simple returned its result without
    # going through validate_directives.  The H14 Item 3 fix
    # (sentinel + Hermann-Mauguin shape check on space_group) lived in
    # validate_directives, so it was BYPASSED in the ollama-fallback
    # path: when the ollama LLM failed to return parseable JSON,
    # extract_directives fell back to extract_directives_simple, whose
    # OWN space-group regex at line ~4350 captured the literal
    # "Not explicitly mentio" (truncated by the regex's {1,20}
    # quantifier from "Not explicitly mentioned" in the preprocessed
    # advice), then returned that dict directly to the agent.  Tom's
    # ollama xtriage tutorial run on 2026-05-26 (gate C) confirmed
    # the bug — the bogus space_group value reached the displayed
    # directives even with H14 installed.
    #
    # The fix is to ensure validate_directives runs on ALL paths,
    # not just the LLM-extracted-JSON path.  This makes
    # validate_directives the canonical "final-sanity" step for every
    # consumer of extract_directives* — future validator extensions
    # automatically apply to both paths.
    #
    # Prerequisite: VALID_STOP_CONDITIONS must include
    # start_with_program (added in H14.1 above) — extract_directives_simple
    # legitimately sets that key via _resolve_after_program when the
    # advice has multiple actions + stop intent, and pre-H14.1
    # validate_directives would have stripped it.
    directives = validate_directives(directives, log)

    return directives


# =============================================================================
# CRYSTAL SYMMETRY EXTRACTION HELPER (shared by simple and LLM paths)
# =============================================================================

def _extract_polder_selection(user_advice):
    """Extract polder atom selection from user advice.

    Matches patterns like:
        "polder map for chain A resseq 88"
        "polder for chain A and resseq 88"
        "polder selection chain A and resseq 88"
        "calculate polder map for resname LIG"

    Returns:
        str or None — the selection string, or None.
    """
    if not user_advice:
        return None

    # Pattern: "for <selection>" after polder mention
    m = re.search(
        r'polder\s+(?:map\s+)?for\s+(.+?)(?:\.|$)',
        user_advice, re.IGNORECASE)
    if m:
        sel = m.group(1).strip()
        if sel and len(sel) > 3:
            return sel

    # Pattern: "selection <selection>" near polder
    m = re.search(
        r'(?:polder|omit)\s+(?:map\s+)?'
        r'selection[=:\s]+["\']?(.+?)["\']?'
        r'(?:\.|$)',
        user_advice, re.IGNORECASE)
    if m:
        sel = m.group(1).strip()
        if sel and len(sel) > 3:
            return sel

    # Pattern: "selection='chain A and resseq 88'"
    m = re.search(
        r'selection\s*=\s*["\']([^"\']+)["\']',
        user_advice, re.IGNORECASE)
    if m:
        return m.group(1).strip()

    return None


def _extract_crystal_symmetry_simple(user_advice, directives):
    """
    Extract unit_cell and space_group from user advice using regex patterns.

    Mutates *directives* in place (adding/merging program_settings.default).

    Handles common formats::

        "unit cell (116.097, 116.097, 44.175, 90, 90, 120)"
        "unit_cell = 116.097 116.097 44.175 90 90 120"
        "use unit cell 116 116 44 90 90 120"
        "space group P 32 2 1"  /  "space_group=P212121"
    """
    NUMBER = r'[0-9]+(?:\.[0-9]*)?'
    SEP    = r'[\s,]+'

    # --- unit cell ---
    UC_PATTERNS = [
        # parenthesised: (a, b, c, al, be, ga)
        r'unit[_ ]?cell\s*[=:]?\s*\(\s*'
        r'(?P<a>{n}){s}(?P<b>{n}){s}(?P<c>{n}){s}'
        r'(?P<al>{n}){s}(?P<be>{n}){s}(?P<ga>{n})\s*\)',
        # plain key=value: unit_cell = a b c al be ga
        r'unit[_ ]?cell\s*[=:]\s*'
        r'(?P<a>{n})\s+(?P<b>{n})\s+(?P<c>{n})\s+'
        r'(?P<al>{n})\s+(?P<be>{n})\s+(?P<ga>{n})',
        # prose: "use unit cell a b c al be ga"
        r'unit\s+cell\s+'
        r'(?P<a>{n})\s+(?P<b>{n})\s+(?P<c>{n})\s+'
        r'(?P<al>{n})\s+(?P<be>{n})\s+(?P<ga>{n})',
    ]
    for raw_pat in UC_PATTERNS:
        pat = raw_pat.format(n=NUMBER, s=SEP)
        m = re.search(pat, user_advice, re.IGNORECASE)
        if m:
            raw = '%s %s %s %s %s %s' % (
                m.group('a'), m.group('b'), m.group('c'),
                m.group('al'), m.group('be'), m.group('ga'))
            normalized = _normalize_unit_cell(raw)
            if normalized:
                if "program_settings" not in directives:
                    directives["program_settings"] = {}
                if "default" not in directives["program_settings"]:
                    directives["program_settings"]["default"] = {}
                directives["program_settings"]["default"]["unit_cell"] = normalized
            break

    # --- space group ---
    m = re.search(r'space[_ ]?group\s*[=:]?\s*([A-Za-z][A-Za-z0-9 _\-]{1,20})',
                  user_advice, re.IGNORECASE)
    if m:
        sg = m.group(1).strip().rstrip('.,;')
        if sg:
            if "program_settings" not in directives:
                directives["program_settings"] = {}
            if "default" not in directives["program_settings"]:
                directives["program_settings"]["default"] = {}
            directives["program_settings"]["default"]["space_group"] = sg
