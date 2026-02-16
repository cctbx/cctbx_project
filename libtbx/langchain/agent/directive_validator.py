"""
Directive Validator for PHENIX AI Agent.

This module provides two types of validation:

1. PRE-VALIDATION (at startup):
   Validates user directives against available capabilities BEFORE starting
   the workflow. If a user requests something we cannot deliver (e.g., a
   program that doesn't exist, or a parameter we don't support), we stop
   immediately with a clear error message.

2. RUNTIME VALIDATION (during execution):
   Applies user directives to modify the LLM's intent during each cycle.
   This includes applying program settings and checking stop conditions.

All program and parameter information is loaded dynamically from the
programs.yaml file, so this validator stays in sync automatically.

Pre-validation usage:
    from agent.directive_validator import validate_directives

    result = validate_directives(
        user_advice="Calculate polder maps using phenix.polder",
        directives=extracted_directives,  # From directive_extractor
    )

    if not result.valid:
        print(f"Cannot proceed: {result.message}")
        for issue in result.issues:
            print(f"  - {issue}")

Runtime validation usage:
    from agent.directive_validator import validate_intent

    result = validate_intent(
        intent={"program": "phenix.refine", "strategy": {}},
        directives={"program_settings": {"default": {"resolution": 2.5}}},
        cycle_number=1,
    )
    # result["validated_intent"] has the modified intent
    # result["modifications"] lists changes made
"""

from __future__ import absolute_import, division, print_function

import copy
import os
import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set, Tuple

# Silence unused import warnings (these are used in type hints)
assert Dict is not None
assert Optional is not None
assert Set is not None


# =============================================================================
# VALIDATION RESULT
# =============================================================================

@dataclass
class ValidationResult:
    """Result of directive validation."""
    valid: bool
    message: str = ""
    issues: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # Detailed breakdown
    unavailable_programs: List[str] = field(default_factory=list)
    unsupported_parameters: List[Tuple[str, str]] = field(default_factory=list)  # (program, param)


# =============================================================================
# DYNAMIC PROGRAM LOADING FROM YAML
# =============================================================================

# Cache for loaded programs
_cached_programs = None
_cached_all_phenix_programs = None


def _load_available_programs():
    """
    Load the set of available programs from programs.yaml via ProgramRegistry.

    Returns:
        Dict mapping program name to its definition including strategy_flags
    """
    global _cached_programs

    if _cached_programs is not None:
        return _cached_programs

    try:
        # Import the program registry which loads from YAML
        try:
            from libtbx.langchain.agent.program_registry import ProgramRegistry
        except ImportError:
            from agent.program_registry import ProgramRegistry

        registry = ProgramRegistry()
        programs = {}

        # Get all program names from the registry
        for prog_name in registry.list_programs():
            prog_def = registry.get_program(prog_name)
            if prog_def:
                programs[prog_name] = {
                    "strategy_flags": list(registry.get_strategy_flags(prog_name).keys()),
                    "description": prog_def.get("description", ""),
                    "category": prog_def.get("category", ""),
                    "experiment_types": prog_def.get("experiment_types", []),
                }

        _cached_programs = programs
        return programs

    except Exception as e:
        # Fallback: try to load YAML directly
        try:
            return _load_programs_from_yaml_directly()
        except Exception:
            # Last resort: return hardcoded minimal set
            print(f"Warning: Could not load programs from registry: {e}")
            return _get_fallback_programs()


def _load_programs_from_yaml_directly():
    """
    Fallback: Load programs directly from YAML file without using ProgramRegistry.
    """
    import yaml

    # Find the YAML file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    yaml_paths = [
        os.path.join(current_dir, "..", "knowledge", "programs.yaml"),
        os.path.join(current_dir, "knowledge", "programs.yaml"),
    ]

    yaml_path = None
    for path in yaml_paths:
        if os.path.exists(path):
            yaml_path = path
            break

    if not yaml_path:
        raise FileNotFoundError("Could not find programs.yaml")

    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)

    programs = {}
    for prog_name, prog_def in data.items():
        if prog_name.startswith("phenix.") and isinstance(prog_def, dict):
            strategy_flags = list(prog_def.get("strategy_flags", {}).keys())
            programs[prog_name] = {
                "strategy_flags": strategy_flags,
                "description": prog_def.get("description", ""),
                "category": prog_def.get("category", ""),
                "experiment_types": prog_def.get("experiment_types", []),
            }

    return programs


def _get_fallback_programs():
    """Minimal fallback list if all else fails."""
    return {
        "phenix.xtriage": {"strategy_flags": [], "description": "Data analysis"},
        "phenix.mtriage": {"strategy_flags": ["resolution"], "description": "Map analysis"},
        "phenix.phaser": {"strategy_flags": [], "description": "Molecular replacement"},
        "phenix.refine": {"strategy_flags": ["resolution", "cycles", "anisotropic_adp", "add_waters"],
                        "description": "Refinement"},
        "phenix.real_space_refine": {"strategy_flags": ["resolution", "cycles"], "description": "Real-space refinement"},
        "phenix.autobuild": {"strategy_flags": [], "description": "Automated model building"},
        "phenix.autosol": {"strategy_flags": ["resolution", "atom_type", "wavelength", "sites"],
                         "description": "Experimental phasing"},
        "phenix.ligandfit": {"strategy_flags": [], "description": "Ligand fitting"},
        "phenix.dock_in_map": {"strategy_flags": ["resolution"], "description": "Dock model in map"},
        "phenix.predict_and_build": {"strategy_flags": [], "description": "AlphaFold prediction"},
        "phenix.molprobity": {"strategy_flags": [], "description": "Validation"},
        "phenix.model_vs_data": {"strategy_flags": [], "description": "Model validation"},
    }


def _get_all_phenix_programs():
    """
    Get a list of ALL known PHENIX programs (not just those in the agent).

    This is used to give better error messages - we can tell users
    "this program exists but isn't in the agent" vs "this program doesn't exist".

    Returns:
        Set of all known PHENIX program names
    """
    global _cached_all_phenix_programs

    if _cached_all_phenix_programs is not None:
        return _cached_all_phenix_programs

    # These are programs that exist in PHENIX but are NOT in the agent workflow
    # We maintain this list to give better error messages
    known_but_unavailable = {
        # Map calculation tools
        "phenix.fem",         # Feature-enhanced maps
        "phenix.maps",        # General map calculation
        "phenix.fmodel",      # Structure factor calculation

        # Model preparation tools
        "phenix.elbow",       # Ligand restraint generation
        "phenix.ready_set",   # Model preparation
        "phenix.reduce",      # Add hydrogens
        "phenix.pdbtools",    # PDB manipulation (some functions)

        # Analysis tools not in agent
        "phenix.find_alt_orig_sym_mate",  # Alternative origins
        "phenix.superpose_pdbs",          # Structure alignment
        "phenix.get_cc_mtz_mtz",          # Map correlation
        "phenix.get_cc_mtz_pdb",          # Model-map correlation

        # Other tools
        "phenix.sculptor",    # Model preparation for MR
        "phenix.ensembler",   # Ensemble preparation
        "phenix.mr_rosetta",  # MR with Rosetta
        "phenix.morph_model", # Model morphing
        "phenix.den_refine",  # DEN refinement
        "phenix.rosetta_refine",  # Rosetta refinement
    }

    # Combine with programs that ARE available
    available = set(_load_available_programs().keys())
    _cached_all_phenix_programs = available | known_but_unavailable

    return _cached_all_phenix_programs


def _normalize_program_name(name: str) -> str:
    """
    Normalize a program name to canonical form.

    Args:
        name: Program name (with or without phenix. prefix)

    Returns:
        Normalized name with phenix. prefix
    """
    name = name.lower().strip()

    # Handle common aliases - include various case/format variations
    aliases = {
        # Refinement
        "rsr": "phenix.real_space_refine",
        "real_space_refine": "phenix.real_space_refine",
        "realspacerefine": "phenix.real_space_refine",
        "real-space-refine": "phenix.real_space_refine",
        "refine": "phenix.refine",
        "refinement": "phenix.refine",

        # Phasing
        "phaser": "phenix.phaser",
        "autosol": "phenix.autosol",

        # Model building - X-ray
        "autobuild": "phenix.autobuild",
        "auto_build": "phenix.autobuild",

        # Model building - Cryo-EM
        "map_to_model": "phenix.map_to_model",
        "maptomodel": "phenix.map_to_model",
        "map-to-model": "phenix.map_to_model",

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

        # Density modification
        "resolve_cryo_em": "phenix.resolve_cryo_em",
        "resolvecryoem": "phenix.resolve_cryo_em",
        "resolve-cryo-em": "phenix.resolve_cryo_em",
        "autobuild_denmod": "phenix.autobuild_denmod",

        # Map tools
        "map_sharpening": "phenix.map_sharpening",
        "mapsharpening": "phenix.map_sharpening",
        "map_symmetry": "phenix.map_symmetry",
        "mapsymmetry": "phenix.map_symmetry",

        # Analysis
        "xtriage": "phenix.xtriage",
        "mtriage": "phenix.mtriage",

        # Validation
        "molprobity": "phenix.molprobity",
        "model_vs_data": "phenix.model_vs_data",
        "modelvsdata": "phenix.model_vs_data",
        "validation_cryoem": "phenix.validation_cryoem",
        "validationcryoem": "phenix.validation_cryoem",
        "holton_geometry_validation": "phenix.holton_geometry_validation",

        # Ligand fitting
        "ligandfit": "phenix.ligandfit",
        "ligand_fit": "phenix.ligandfit",

        # PDB tools
        "pdbtools": "phenix.pdbtools",
        "pdb_tools": "phenix.pdbtools",

        # Known unavailable programs (still normalize for better error messages)
        "polder": "phenix.polder",
        "fem": "phenix.fem",
        "elbow": "phenix.elbow",
        "maps": "phenix.maps",
        "fmodel": "phenix.fmodel",
        "ready_set": "phenix.ready_set",
        "reduce": "phenix.reduce",
    }

    if name in aliases:
        return aliases[name]

    if not name.startswith("phenix."):
        return f"phenix.{name}"

    return name


# =============================================================================
# PARAMETER MAPPING
# =============================================================================

def _get_parameter_aliases():
    """
    Get mapping of user-friendly parameter names to canonical names.

    These are loaded dynamically based on what's in strategy_flags.
    """
    return {
        # Resolution
        "resolution": "resolution",
        "res": "resolution",
        "high_resolution": "resolution",
        "d_min": "resolution",

        # Cycles - map various user terms to our 'cycles' parameter
        "macro_cycles": "cycles",
        "macro-cycles": "cycles",
        "macrocycles": "cycles",
        "cycles": "cycles",
        "refinement_cycles": "cycles",
        "ncycles": "cycles",
        "number_of_cycles": "cycles",

        # B-factors
        "anisotropic": "anisotropic_adp",
        "anisotropic_b": "anisotropic_adp",
        "anisotropic_adp": "anisotropic_adp",
        "aniso": "anisotropic_adp",

        # Waters
        "waters": "add_waters",
        "add_waters": "add_waters",
        "ordered_solvent": "add_waters",
        "solvent": "add_waters",

        # Simulated annealing
        "simulated_annealing": "simulated_annealing",
        "sa": "simulated_annealing",
        "anneal": "simulated_annealing",

        # Hydrogens
        "riding_hydrogens": "riding_hydrogens",
        "hydrogens": "riding_hydrogens",
        "h_atoms": "riding_hydrogens",

        # Twin law
        "twin_law": "twin_law",
        "twinning": "twin_law",
        "twin": "twin_law",

        # R-free
        "rfree": "generate_rfree_flags",
        "r_free": "generate_rfree_flags",
        "generate_rfree": "generate_rfree_flags",

        # Anomalous scattering
        "atom_type": "atom_type",
        "anomalous_scatterer": "atom_type",
        "heavy_atom": "atom_type",
        "wavelength": "wavelength",
        "sites": "sites",
        "num_sites": "sites",
        "number_of_sites": "sites",
    }


# =============================================================================
# TEXT EXTRACTION
# =============================================================================

def _extract_program_references(text: str) -> List[str]:
    """
    Extract program references from user advice text.

    Looks for:
    - Explicit program names: phenix.polder, phenix.refine
    - Implicit references: "run polder", "use autobuild"
    - Goal statements: "calculate polder maps"

    Returns:
        List of normalized program names
    """
    if not text:
        return []

    programs = []
    text_lower = text.lower()

    # Pattern 1: Explicit phenix.X references
    explicit_pattern = r'phenix\.(\w+)'
    for match in re.finditer(explicit_pattern, text_lower):
        prog_name = f"phenix.{match.group(1)}"
        normalized = _normalize_program_name(prog_name)
        if normalized not in programs:
            programs.append(normalized)

    # Pattern 2: "run/use/with X" patterns for known program names
    all_programs = _get_all_phenix_programs()
    # Extract just the short names (without phenix. prefix)
    short_names = {p.replace("phenix.", "") for p in all_programs}

    run_patterns = [
        r'(?:run|use|using|with|via|through)\s+(\w+)',
        r'(\w+)\s+(?:map|maps|analysis|calculation)',
        r'(?:calculate|compute|generate)\s+(\w+)\s+(?:map|maps)',
    ]

    for pattern in run_patterns:
        for match in re.finditer(pattern, text_lower):
            word = match.group(1)
            if word in short_names:
                normalized = _normalize_program_name(word)
                if normalized not in programs:
                    programs.append(normalized)

    return programs


def _extract_parameter_references(text: str) -> List[Tuple[str, str]]:
    """
    Extract parameter references from user advice text.

    Returns:
        List of (parameter_name, context) tuples
    """
    if not text:
        return []

    params = []
    text_lower = text.lower()

    # Pattern: "N macro cycles" or "macro cycles = N"
    macro_cycle_patterns = [
        r'(\d+)\s*(?:macro[- ]?cycles?)',
        r'(?:macro[- ]?cycles?)\s*[=:]\s*(\d+)',
        r'(?:use|with|set)\s+(\d+)\s+(?:macro[- ]?cycles?)',
        r'one\s+macro[- ]?cycle',
        r'single\s+macro[- ]?cycle',
    ]

    for pattern in macro_cycle_patterns:
        if re.search(pattern, text_lower):
            params.append(("macro_cycles", "refinement"))
            break

    # Pattern: resolution references
    if re.search(r'resolution\s*(?:of|=|:)?\s*[\d.]+', text_lower) or \
       re.search(r'[\d.]+\s*(?:angstrom|Å|A)\s+resolution', text_lower):
        params.append(("resolution", "general"))

    # Pattern: anisotropic
    if re.search(r'anisotropic', text_lower):
        params.append(("anisotropic", "refinement"))

    # Pattern: waters/solvent
    if re.search(r'(?:add|ordered)\s+(?:water|solvent)', text_lower):
        params.append(("waters", "refinement"))

    # Pattern: simulated annealing
    if re.search(r'simulated\s+anneal', text_lower):
        params.append(("simulated_annealing", "refinement"))

    return params


# =============================================================================
# MAIN VALIDATION FUNCTION
# =============================================================================

def validate_directives(
    user_advice: str,
    directives: Optional[Dict] = None,
    available_programs: Optional[Dict] = None,
) -> ValidationResult:
    """
    Validate user directives against available capabilities.

    This should be called BEFORE starting the workflow to catch
    impossible requests early.

    Args:
        user_advice: Raw user advice text
        directives: Extracted directives dict (from directive_extractor)
        available_programs: Optional dict of available programs (for testing)

    Returns:
        ValidationResult with valid=False if there are blocking issues
    """
    issues = []
    warnings = []
    unavailable_programs = []
    unsupported_parameters = []

    # Load available programs from YAML
    if available_programs is None:
        available_programs = _load_available_programs()

    available_set = set(available_programs.keys())
    all_phenix = _get_all_phenix_programs()

    # =================================================================
    # 1. CHECK PROGRAM REFERENCES IN USER ADVICE (WARNINGS ONLY)
    # =================================================================
    # Programs mentioned in text are informational - don't block workflow.
    # Only programs in directives structure should block.

    referenced_programs = _extract_program_references(user_advice)

    for prog in referenced_programs:
        if prog in available_set:
            # Program is available - good!
            continue
        elif prog in all_phenix:
            # Program exists but isn't in agent workflow
            # This is just a mention in text - warn but don't block
            unavailable_programs.append(prog)
            warnings.append(
                f"Note: '{prog}' mentioned in advice exists in PHENIX but is not available in the AI agent. "
                f"You can run it manually if needed: {prog} <input_files>"
            )
        else:
            # Program not recognized - might be a typo, just warn
            unavailable_programs.append(prog)
            suggestions = _suggest_similar_programs(prog, available_set)
            if suggestions:
                warnings.append(
                    f"Note: '{prog}' mentioned in advice is not recognized. "
                    f"Did you mean: {', '.join(suggestions)}?"
                )
            # Don't warn about completely unknown programs - might be non-phenix tools

    # =================================================================
    # 2. CHECK DIRECTIVES FOR PROGRAM REFERENCES (BLOCKING)
    # =================================================================

    if directives:
        # Check stop_conditions.after_program
        stop_conditions = directives.get("stop_conditions", {})
        after_prog = stop_conditions.get("after_program")

        if after_prog:
            normalized = _normalize_program_name(after_prog)
            if normalized not in available_set:
                if normalized not in unavailable_programs:
                    unavailable_programs.append(normalized)
                    if normalized in all_phenix:
                        issues.append(
                            f"Requested stop program '{after_prog}' is not available in the agent workflow."
                        )
                    else:
                        issues.append(
                            f"Requested stop program '{after_prog}' is not recognized."
                        )

        # Check workflow_preferences
        workflow_prefs = directives.get("workflow_preferences", {})
        for pref_list in ["skip_programs", "prefer_programs"]:
            for prog in workflow_prefs.get(pref_list, []):
                normalized = _normalize_program_name(prog)
                if normalized not in available_set and normalized not in unavailable_programs:
                    warnings.append(
                        f"Program '{prog}' in {pref_list} is not available in the agent."
                    )

        # Check program_settings for unknown programs
        prog_settings = directives.get("program_settings", {})
        for prog in prog_settings.keys():
            if prog != "default":
                normalized = _normalize_program_name(prog)
                if normalized not in available_set:
                    warnings.append(
                        f"Settings specified for '{prog}' which is not available in the agent."
                    )

    # =================================================================
    # 3. CHECK PARAMETER SUPPORT
    # =================================================================

    param_aliases = _get_parameter_aliases()
    param_refs = _extract_parameter_references(user_advice)

    for param, context in param_refs:
        canonical_param = param_aliases.get(param, param)

        # macro_cycles IS directly controllable via program_settings directives
        # (applied by command_builder._build_strategy)
        if param in ["macro_cycles", "macro-cycles", "macrocycles"]:
            continue  # Supported - no warning needed

        # Check if parameter is supported by relevant programs
        if context == "refinement" and canonical_param != "resolution":
            supported = False
            for prog in ["phenix.refine", "phenix.real_space_refine"]:
                if prog in available_programs:
                    prog_params = available_programs[prog].get("strategy_flags", [])
                    if canonical_param in prog_params:
                        supported = True
                        break

            if not supported:
                unsupported_parameters.append(("phenix.refine", param))
                warnings.append(
                    f"Parameter '{param}' may not be directly supported for refinement."
                )

    # =================================================================
    # BUILD RESULT
    # =================================================================

    if issues:
        message = "Cannot proceed with the requested workflow:\n"
        for i, issue in enumerate(issues, 1):
            message += f"  {i}. {issue}\n"

        return ValidationResult(
            valid=False,
            message=message,
            issues=issues,
            warnings=warnings,
            unavailable_programs=unavailable_programs,
            unsupported_parameters=unsupported_parameters,
        )

    if warnings:
        message = "Workflow can proceed, but note:\n"
        for w in warnings:
            message += f"  - {w}\n"

        return ValidationResult(
            valid=True,
            message=message,
            issues=[],
            warnings=warnings,
            unavailable_programs=[],
            unsupported_parameters=unsupported_parameters,
        )

    return ValidationResult(
        valid=True,
        message="All requested capabilities are available.",
    )


def _suggest_similar_programs(prog: str, available: Set[str]) -> List[str]:
    """
    Suggest similar program names for typos.
    """
    prog_lower = prog.lower().replace("phenix.", "")
    suggestions = []

    for avail in available:
        avail_short = avail.replace("phenix.", "")
        # Simple similarity: shared prefix or substring
        if prog_lower in avail_short or avail_short in prog_lower:
            suggestions.append(avail)
        elif len(set(prog_lower) & set(avail_short)) > len(prog_lower) * 0.5:
            suggestions.append(avail)

    return suggestions[:3]


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def check_program_available(program_name: str) -> Tuple[bool, str]:
    """
    Quick check if a specific program is available.
    """
    available = _load_available_programs()
    normalized = _normalize_program_name(program_name)

    if normalized in available:
        return True, f"'{normalized}' is available."
    elif normalized in _get_all_phenix_programs():
        return False, f"'{normalized}' exists but is not available in the AI agent workflow."
    else:
        return False, f"'{normalized}' is not recognized."


def get_supported_parameters(program_name: str) -> List[str]:
    """
    Get list of supported parameters for a program.
    """
    available = _load_available_programs()
    normalized = _normalize_program_name(program_name)

    if normalized in available:
        return available[normalized].get("strategy_flags", [])
    return []


def list_available_programs() -> List[str]:
    """
    List all programs available in the agent workflow.
    """
    return sorted(_load_available_programs().keys())


def list_unavailable_programs() -> List[str]:
    """
    List known PHENIX programs that are NOT in the agent workflow.
    """
    available = set(_load_available_programs().keys())
    all_known = _get_all_phenix_programs()
    return sorted(all_known - available)


# =============================================================================
# RUNTIME VALIDATION FUNCTIONS
# =============================================================================
# These functions are used during agent execution to apply user directives
# to the LLM's intent, modifying parameters as needed.

def validate_intent(intent, directives, cycle_number=1, history=None, attempt_number=0, log_func=None):
    """
    Validate and modify an LLM intent against user directives.

    This applies program_settings from directives to the intent.

    Behavior depends on attempt_number:
    - attempt_number=0 (first try): Override LLM with directive values (honor user's explicit request)
    - attempt_number>0 (retry): Trust LLM's interpretation (it may be correcting syntax errors)

    Args:
        intent: Dict with 'program', 'strategy', 'files', 'reasoning'
        directives: Dict with 'program_settings', 'stop_conditions', etc.
        cycle_number: Current cycle number
        history: List of past cycle results
        attempt_number: Current attempt within this cycle (0=first, 1+=retry)
        log_func: Optional logging function

    Returns:
        Dict with:
        - validated_intent: Modified intent
        - modifications: List of changes made
        - warnings: List of conflicts detected
    """
    if log_func is None:
        log_func = lambda x: None

    if history is None:
        history = []

    # Make a copy to avoid mutating original
    validated = copy.deepcopy(intent)
    modifications = []
    warnings = []

    # Apply program settings from directives
    if directives:
        mods, warns = _apply_program_settings(
            validated, directives, validated.get("program", ""),
            attempt_number=attempt_number, log_func=log_func
        )
        modifications.extend(mods)
        warnings.extend(warns)

    return {
        "validated_intent": validated,
        "modifications": modifications,
        "warnings": warnings,
    }


def _apply_program_settings(intent, directives, program, attempt_number=0, log_func=None):
    """
    Apply program_settings from directives to intent.

    Behavior depends on attempt_number:
    - attempt_number=0: Override LLM with directive values (honor user's explicit request)
    - attempt_number>0: Trust LLM's interpretation (it may be correcting syntax errors)

    This is safer than always trusting LLM or always trusting directives:
    - First attempt honors the user's explicit request
    - Retries allow LLM to fix potential syntax issues

    Args:
        intent: Intent dict (modified in place)
        directives: Directives dict
        program: Program name
        attempt_number: Current attempt (0=first, 1+=retry)
        log_func: Logging function

    Returns:
        (modifications, warnings) tuple
    """
    if log_func is None:
        log_func = lambda x: None

    modifications = []
    warnings = []

    program_settings = directives.get("program_settings", {})
    if not program_settings:
        return modifications, warnings

    # Get settings for this program (specific overrides default)
    default_settings = program_settings.get("default", {})
    specific_settings = program_settings.get(program, {})

    # Merge: default first, then specific overrides
    combined_settings = {}
    combined_settings.update(default_settings)
    combined_settings.update(specific_settings)

    if not combined_settings:
        return modifications, warnings

    # Ensure strategy dict exists
    if "strategy" not in intent or intent["strategy"] is None:
        intent["strategy"] = {}

    # Apply settings
    # First normalize: for RSR, "cycles" and "macro_cycles" are the same thing
    # Deduplicate to avoid logging the same override twice
    if program == "phenix.real_space_refine":
        if "cycles" in combined_settings and "macro_cycles" in combined_settings:
            # Both present — keep only macro_cycles (the canonical name for RSR)
            del combined_settings["cycles"]
        elif "cycles" in combined_settings:
            combined_settings["macro_cycles"] = combined_settings.pop("cycles")

    applied_normalized = set()  # Track normalized keys to avoid duplicates
    for key, value in combined_settings.items():
        existing = intent["strategy"].get(key)

        if existing is None:
            # Add missing setting from directive
            intent["strategy"][key] = value
            modifications.append(f"Added {key}={value} from directives")
            log_func(f"  Added {key}={value}")
        elif existing != value:
            # LLM chose different value than directive
            if attempt_number == 0:
                # First attempt: honor user's explicit directive
                warnings.append(
                    f"LLM chose {key}={existing}, directive specifies {value} - using directive (first attempt)"
                )
                intent["strategy"][key] = value
                log_func(f"  Override {key}: {existing} -> {value} (honoring directive)")
            else:
                # Retry: trust LLM's interpretation (it may be correcting syntax)
                warnings.append(
                    f"LLM interpreted '{value}' as {key}={existing} - using LLM (retry attempt {attempt_number})"
                )
                log_func(f"  Trusting LLM: {key}={existing} (retry)")
                # Keep LLM's value (don't override)

    return modifications, warnings


def augment_intent_with_directives(intent, directives):
    """
    Add missing settings from directives without overriding existing values.

    This is a lighter version of validate_intent that only fills in gaps.

    Args:
        intent: Intent dict
        directives: Directives dict

    Returns:
        Modified intent (copy)
    """
    result = copy.deepcopy(intent)

    program = result.get("program", "")
    if program == "STOP":
        return result

    program_settings = directives.get("program_settings", {}) if directives else {}
    if not program_settings:
        return result

    # Get settings for this program
    default_settings = program_settings.get("default", {})
    specific_settings = program_settings.get(program, {})

    # Merge: default first, then specific overrides
    combined_settings = {}
    combined_settings.update(default_settings)
    combined_settings.update(specific_settings)

    if not combined_settings:
        return result

    # Ensure strategy dict exists
    if "strategy" not in result or result["strategy"] is None:
        result["strategy"] = {}

    # Only add missing settings
    for key, value in combined_settings.items():
        if key not in result["strategy"]:
            result["strategy"][key] = value

    return result


def get_stop_reason_from_directives(directives, cycle_number, last_program, history=None):
    """
    Check if stop conditions are met and return reason.

    Args:
        directives: Directives dict
        cycle_number: Current cycle number
        last_program: Last program that ran
        history: List of past cycle results

    Returns:
        Stop reason string if should stop, None otherwise
    """
    if not directives:
        return None

    should_stop, reason = _check_stop_conditions(
        directives, cycle_number, last_program, history or [], lambda x: None
    )

    return reason if should_stop else None


def _check_stop_conditions(directives, cycle_number=None, last_program=None, history=None, log_func=None, log=None):
    """
    Check all stop conditions from directives.

    Args:
        directives: Directives dict
        cycle_number: Current cycle number
        last_program: Last program that ran
        history: List of past cycle results
        log_func: Logging function (or use 'log' for compatibility)
        log: Alias for log_func (for backward compatibility)

    Returns:
        (should_stop, reason) tuple
    """
    # Support both log_func and log parameter names
    if log_func is None:
        log_func = log if log is not None else lambda x: None

    if history is None:
        history = []

    stop_conditions = directives.get("stop_conditions", {})
    if not stop_conditions:
        return False, None

    # Check after_cycle
    after_cycle = stop_conditions.get("after_cycle")
    if after_cycle is not None and cycle_number is not None and cycle_number >= after_cycle:
        return True, f"Reached cycle {cycle_number} (stop after cycle {after_cycle})"

    # Check after_program
    after_program = stop_conditions.get("after_program")
    if after_program and last_program == after_program:
        return True, f"Completed {after_program} (stop after this program)"

    # Check R-free target
    r_free_target = stop_conditions.get("r_free_target")
    if r_free_target is not None and history:
        for entry in reversed(history):
            metrics = entry.get("metrics", {})
            r_free = metrics.get("r_free")
            if r_free is not None and r_free <= r_free_target:
                return True, f"R-free {r_free:.3f} reached target {r_free_target}"

    # Check map CC target
    map_cc_target = stop_conditions.get("map_cc_target")
    if map_cc_target is not None and history:
        for entry in reversed(history):
            metrics = entry.get("metrics", {})
            map_cc = metrics.get("map_cc") or metrics.get("cc_mask")
            if map_cc is not None and map_cc >= map_cc_target:
                return True, f"Map CC {map_cc:.3f} reached target {map_cc_target}"

    # Check skip_validation
    if stop_conditions.get("skip_validation"):
        return True, "skip_validation directive is active"

    # Check max program cycles
    should_stop, reason = _check_max_program_cycles(
        directives, last_program, history, log_func
    )
    if should_stop:
        return True, reason

    return False, None


def _check_max_program_cycles(directives, program, history, log_func):
    """
    Check if max cycles for a specific program have been reached.

    Args:
        directives: Directives dict
        program: Program name to check
        history: List of past cycle results
        log_func: Logging function

    Returns:
        (should_stop, reason) tuple
    """
    stop_conditions = directives.get("stop_conditions", {})

    # Check max_refine_cycles
    max_refine = stop_conditions.get("max_refine_cycles")
    if max_refine is not None and program in ("phenix.refine", "phenix.real_space_refine"):
        refine_count = sum(
            1 for h in history
            if h.get("program") in ("phenix.refine", "phenix.real_space_refine")
        )
        if refine_count >= max_refine:
            return True, f"Reached max refinement cycles ({max_refine})"

    # Check max_build_cycles
    max_build = stop_conditions.get("max_build_cycles")
    if max_build is not None and program in ("phenix.autobuild", "phenix.map_to_model"):
        build_count = sum(
            1 for h in history
            if h.get("program") in ("phenix.autobuild", "phenix.map_to_model")
        )
        if build_count >= max_build:
            return True, f"Reached max build cycles ({max_build})"

    return False, None


def format_validation_result(result):
    """
    Format a validation result for display.

    Args:
        result: Dict from validate_intent

    Returns:
        Formatted string
    """
    lines = []

    modifications = result.get("modifications", [])
    warnings = result.get("warnings", [])
    should_stop = result.get("should_stop", False)
    stop_reason = result.get("stop_reason")

    if modifications:
        lines.append("Modifications:")
        for mod in modifications:
            lines.append(f"  - {mod}")

    if warnings:
        lines.append("Warnings:")
        for warn in warnings:
            lines.append(f"  - {warn}")

    if should_stop:
        lines.append(f"Stop: {stop_reason or 'Yes'}")

    if not lines:
        return "No changes from directives"

    return "\n".join(lines)


# =============================================================================
# MAIN (for testing)
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Directive Validator Test")
    print("=" * 60)

    # Test 1: Unavailable program
    print("\n1. Testing unavailable program (phenix.fem):")
    result = validate_directives(
        "Calculate feature-enhanced maps using phenix.fem"
    )
    print(f"   Valid: {result.valid}")
    print(f"   Message: {result.message}")

    # Test 2: Available program (polder is now available!)
    print("\n2. Testing available program (phenix.polder):")
    result = validate_directives(
        "Calculate polder maps using phenix.polder for the ligand"
    )
    print(f"   Valid: {result.valid}")
    print(f"   Message: {result.message}")

    # Test 3: Unsupported parameter
    print("\n3. Testing parameter request (macro_cycles):")
    result = validate_directives(
        "Use one macro cycle in refinement"
    )
    print(f"   Valid: {result.valid}")
    print(f"   Warnings: {result.warnings}")

    # Test 4: Multiple issues
    print("\n4. Testing multiple issues:")
    result = validate_directives(
        "Run phenix.fem and then phenix.maps for general map calculation"
    )
    print(f"   Valid: {result.valid}")
    print(f"   Issues: {result.issues}")

    print("\n" + "=" * 60)
    print("Programs available in agent (from YAML):")
    for prog in list_available_programs():
        params = get_supported_parameters(prog)
        param_str = ", ".join(params[:3]) if params else "(none)"
        if len(params) > 3:
            param_str += "..."
        print(f"  {prog}: [{param_str}]")

    print("\n" + "-" * 60)
    print("Known PHENIX programs NOT in agent:")
    for prog in list_unavailable_programs():
        print(f"  {prog}")

    print("\nTest complete!")

