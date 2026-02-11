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
        model="gemini-2.0-flash"
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
   - For phenix.map_sharpening specifically:
     * resolution: float (required for model-based sharpening)
     * sharpening_method: string (e.g., "b-factor", "model_sharpening")
   - For phenix.polder specifically:
     * selection: string (atom selection, e.g., "chain A and resseq 88", "resname LIG")

2. "stop_conditions": When to stop the workflow
   - "after_program": string - Stop after this program completes (e.g., "phenix.xtriage", "phenix.refine", "phenix.phaser", "phenix.ligandfit", "phenix.map_to_model", "phenix.dock_in_map", "phenix.map_sharpening", "phenix.polder")
   - "after_cycle": int - Stop after this cycle number
   - "max_refine_cycles": int - Maximum number of refinement cycles to run
   - "skip_validation": bool - If true, allow stopping without running molprobity
   - "r_free_target": float - Stop when R-free reaches this value
   - "map_cc_target": float - Stop when map correlation reaches this value

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

**CRITICAL: max_refine_cycles vs after_program**
- max_refine_cycles=N: Limits the NUMBER of refinement jobs to N. The workflow continues normally until refinement, then stops after N refinement jobs.
- after_program="X": FORCES program X to be run IMMEDIATELY, bypassing normal workflow. Only use when user wants to skip directly to a specific program.
- "maximum of one refinement" or "at most one refinement" → ONLY set max_refine_cycles=1, do NOT set after_program
- "solve the structure with one refinement" → max_refine_cycles=1 (workflow proceeds normally: xtriage → model → refine)
- "just run refinement" or "only refinement" → after_program="phenix.refine" (skip to refinement immediately)

**CRITICAL: skip_validation RULE**
If the user specifies ANY explicit stop condition (like "stop after X" or "Stop Condition: ..."),
ALWAYS set skip_validation=true. This tells the system that the user knows what they want
and doesn't need automatic validation before stopping. Examples:
- "Stop after running mtriage" → skip_validation=true
- "Stop after density modification" → skip_validation=true
- "Stop Condition: Stop after generating the improved map" → skip_validation=true

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
If the advice describes a SPECIFIC LIMITED TASK rather than full structure determination, automatically add appropriate stop_conditions:
- "run xtriage", "check for twinning", "analyze data" → after_program="phenix.xtriage", skip_validation=true
- "run phaser", "test MR", "try molecular replacement" → after_program="phenix.phaser", skip_validation=true
- "run one refinement", "quick refinement test" → after_program="phenix.refine", max_refine_cycles=1, skip_validation=true
- "check data quality", "analyze reflection data" → after_program="phenix.xtriage", skip_validation=true
- "run mtriage", "analyze map quality" → after_program="phenix.mtriage", skip_validation=true
- "map symmetry", "determine symmetry", "find symmetry" → after_program="phenix.map_symmetry", skip_validation=true
- "map sharpening", "sharpen the map", "sharpen map", "automatic sharpening" → after_program="phenix.map_sharpening", skip_validation=true
  (Note: phenix.map_sharpening is the dedicated map sharpening tool - use this for sharpening requests)
- For cryo-EM density modification (NOT sharpening): "density modification", "denmod", "resolve_cryo_em" → after_program="phenix.resolve_cryo_em", skip_validation=true
  (Note: phenix.resolve_cryo_em is for density modification, NOT for map sharpening - use phenix.map_sharpening for sharpening)
- For X-ray: "density modification", "improve phases", "denmod" → after_program="phenix.autobuild_denmod", skip_validation=true
  (Note: X-ray density modification uses phenix.autobuild with maps_only=True)
- "MR-SAD", "MR SAD", "MRSAD", "molecular replacement SAD" → Set use_mr_sad=true in workflow_preferences AND use_experimental_phasing=true. Do NOT set after_program="phenix.autosol" because phaser must run first to place the model. The workflow is: phaser → autosol with the phaser output as partpdb_file.
  If user says "stop after autosol" or similar, set after_program="phenix.autosol" AND skip_validation=true.
- "dock in map", "fit model to map" → after_program="phenix.dock_in_map", skip_validation=true
- "map to model", "MapToModel", "build model into map", "automated model building", "rebuild model", "model rebuilding" → after_program="phenix.map_to_model", skip_validation=true
- "polder", "polder map", "omit map", "evaluate ligand placement" → after_program="phenix.polder", skip_validation=true
  (Note: phenix.polder calculates polder omit maps to evaluate ligand/residue placement in density)
- "fit ligand", "ligandfit", "place ligand" → after_program="phenix.ligandfit", skip_validation=true
- Any procedure that ends with a specific analysis step should stop after that step.
- If the stop condition mentions generating a specific output file, set skip_validation=true.

**CRITICAL: WORKFLOW CONTINUATION INDICATORS**:
Do NOT set after_program stop conditions if the user indicates they want additional steps AFTER a program:
- Words like "then", "later", "afterwards", "and then", "continue", "next", "followed by" indicate MULTI-STEP workflows
- "fit the ligand later" → Do NOT stop early, workflow should continue to ligand fitting
- "run predict_and_build, then refine" → Do NOT stop after predict_and_build
- "build model and add ligand" → Do NOT stop after model building
- "refine then validate" → Do NOT stop after refinement
- If user mentions a downstream task (ligand fitting, validation, additional refinement), do NOT set early stop
- Put downstream tasks in "constraints" instead so the agent knows to do them

**CRITICAL: LIGAND FITTING WORKFLOWS**:
When user mentions ligand fitting as part of the workflow:
- "refine, fit ligand, then refine again" → Do NOT set after_program, let workflow complete naturally
- "stop after second refinement" (when ligand fitting is mentioned) → The second refine is AFTER ligandfit, so do NOT set after_program="phenix.refine"
- "one refinement, ligandfit, final refinement" → Do NOT set after_program, this is a complete workflow
- If ligand fitting is mentioned AND refinement after it, the workflow should be: refine → ligandfit → refine → stop
- Put "Fit ligand after first refinement" or similar in constraints, NOT in stop_conditions

Examples of what NOT to do:
- User says "run predict_and_build and fit ligand later" → Do NOT set after_program="phenix.predict_and_build"
- User says "try MR then refine" → Do NOT set after_program="phenix.phaser"
- User says "refine, fit ligand, refine again" → Do NOT set after_program="phenix.refine" (this would stop before ligandfit!)
- User says "stop after the second refinement" (with ligand workflow) → Do NOT set after_program="phenix.refine"

Output ONLY valid JSON. No explanation, no markdown code blocks, just the JSON object."""


# =============================================================================
# EXTRACTION FUNCTION
# =============================================================================

def extract_directives(user_advice, provider="google", model=None, log_func=None,
                       use_rules_only=False):
    """
    Extract structured directives from user advice using LLM.

    Args:
        user_advice: Natural language user instructions
        provider: LLM provider ("google", "openai", "anthropic")
        model: Specific model to use (optional, uses default for provider)
        log_func: Optional logging function
        use_rules_only: If True, skip LLM and use simple pattern extraction

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
    prompt = DIRECTIVE_EXTRACTION_PROMPT.format(user_advice=user_advice)

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

        # If validation stripped everything for ollama, try simple extraction
        if not directives and provider == "ollama":
            log("DIRECTIVES: Validation emptied directives, trying simple extraction")
            simple_directives = extract_directives_simple(user_advice)
            if simple_directives:
                return simple_directives

        log("DIRECTIVES: Extracted %d directive sections" % len(directives))
        return directives

    except Exception as e:
        log("DIRECTIVES: Extraction failed - %s" % str(e))
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


def _call_llm_fallback(prompt, provider, model, log):
    """
    Fallback LLM call when main infrastructure not available.
    Includes retry logic for rate limit errors with exponential backoff and decay.
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
        try:
            # Try new google.genai package first (recommended)
            try:
                from google import genai
                from google.genai import types

                model_name = model or "gemini-2.0-flash"

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

                model_name = model or "gemini-2.0-flash"
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
            log("DIRECTIVES: Google API call failed - %s" % str(e))
            return None

    elif provider == "openai":
        try:
            import openai

            model_name = model or "gpt-4o-mini"
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
            log("DIRECTIVES: OpenAI API call failed - %s" % str(e))
            return None

    elif provider == "anthropic":
        try:
            import anthropic

            model_name = model or "claude-sonnet-4-20250514"
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
            log("DIRECTIVES: Anthropic API call failed - %s" % str(e))
            return None

    elif provider == "ollama":
        try:
            # Ollama uses OpenAI-compatible API
            import openai

            model_name = model or "llama3.2"
            base_url = os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434/v1")

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
            log("DIRECTIVES: Ollama API call failed - %s" % str(e))
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
}

# Valid stop condition keys
VALID_STOP_CONDITIONS = {
    "after_program": str,
    "after_cycle": int,
    "max_refine_cycles": int,
    "skip_validation": bool,
    "r_free_target": float,
    "map_cc_target": float,
}


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

                valid_settings = {}
                for key, value in settings.items():
                    if key in VALID_SETTINGS:
                        try:
                            expected_type = VALID_SETTINGS[key]
                            valid_settings[key] = expected_type(value)
                        except (ValueError, TypeError):
                            _log("DIRECTIVES: Invalid value for %s: %s" % (key, value))
                    else:
                        # Keep unknown settings but log them
                        valid_settings[key] = value
                        _log("DIRECTIVES: Unknown setting %s=%s (keeping)" % (key, value))

                if valid_settings:
                    valid_prog_settings[prog] = valid_settings

            if valid_prog_settings:
                validated["program_settings"] = valid_prog_settings

    # Validate stop_conditions
    if "stop_conditions" in directives:
        stop_cond = directives["stop_conditions"]
        if isinstance(stop_cond, dict):
            valid_stop = {}

            for key, value in stop_cond.items():
                if key in VALID_STOP_CONDITIONS:
                    try:
                        expected_type = VALID_STOP_CONDITIONS[key]
                        if key == "after_program":
                            # Validate program name
                            if value not in VALID_PROGRAMS:
                                fixed = _fix_program_name(value)
                                if fixed:
                                    value = fixed
                                else:
                                    _log("DIRECTIVES: Invalid stop program %s" % value)
                                    continue
                        valid_stop[key] = expected_type(value)
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
                    # Boolean data type preferences
                    valid_prefs[key] = bool(value)

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
                elif key in ("use_experimental_phasing", "use_molecular_replacement"):
                    valid_wf[key] = bool(value)

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

        # Remove the conflicting after_program
        del directives["stop_conditions"]["after_program"]

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

    # Check after_program - normalize names for comparison
    if "after_program" in stop_cond:
        target_program = stop_cond["after_program"]
        # Normalize: remove "phenix." prefix for comparison
        target_normalized = target_program.replace("phenix.", "")
        last_normalized = last_program.replace("phenix.", "") if last_program else ""

        if last_normalized == target_normalized or last_program == target_program:
            return True, "Completed %s (directive: after_program)" % target_program

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
        lines.append("\n  Program Settings:")
        for prog, settings in directives["program_settings"].items():
            settings_str = ", ".join("%s=%s" % (k, v) for k, v in settings.items())
            lines.append("    %s: %s" % (prog, settings_str))

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

def extract_directives_simple(user_advice):
    """
    Extract directives using simple pattern matching (no LLM).

    This is a fallback when LLM is not available or fails.
    Only handles common, unambiguous patterns.

    Args:
        user_advice: User advice string

    Returns:
        dict: Extracted directives (limited)
    """
    if not user_advice:
        return {}

    directives = {}
    advice_lower = user_advice.lower()

    # Extract resolution
    res_patterns = [
        r'resolution\s*[=:]\s*([0-9.]+)',
        r'resolution\s+(?:of\s+)?([0-9.]+)',  # "resolution 3.0" or "resolution of 3.0"
        r'resolution\s+limit\s*[=:]?\s*([0-9.]+)',
        r'([0-9.]+)\s*(?:angstrom|Å|A)\s*resolution',
        r'to\s+([0-9.]+)\s*(?:angstrom|Å|A)',
        r'([0-9.]+)\s*(?:angstrom|Å|A)\b',  # Standalone "3.0 Å"
    ]

    for pattern in res_patterns:
        match = re.search(pattern, user_advice, re.IGNORECASE)
        if match:
            try:
                res = float(match.group(1))
                if 0.5 <= res <= 20:
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

    # Check for skip program patterns
    skip_program_patterns = [
        (r'(?:skip|no|avoid|don\'?t\s+(?:run|use))\s+(?:the\s+)?autobuild', 'phenix.autobuild'),
        (r'(?:skip|no|avoid|don\'?t\s+(?:run|use))\s+(?:the\s+)?ligand\s*fit', 'phenix.ligandfit'),
        (r'(?:skip|no|avoid|don\'?t\s+(?:run|use))\s+(?:the\s+)?molprobity', 'phenix.molprobity'),
        (r'(?:skip|no|avoid|don\'?t\s+(?:run|use))\s+(?:the\s+)?phaser', 'phenix.phaser'),
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

    # Check for explicit "Stop Condition:" section from preprocessed advice
    # This is output from the advice preprocessor and should trigger skip_validation
    if re.search(r'stop\s+condition\s*:', advice_lower, re.IGNORECASE):
        if "stop_conditions" not in directives:
            directives["stop_conditions"] = {}
        # User explicitly specified a stop condition - honor it
        directives["stop_conditions"]["skip_validation"] = True

        # Try to extract program from stop condition text
        stop_match = re.search(r'stop\s+condition\s*:\s*stop\s+after\s+(?:running\s+)?(.+?)(?:\.|$)', advice_lower, re.IGNORECASE)
        if stop_match:
            stop_text = stop_match.group(1).strip()
            # Map common phrases to programs
            # NOTE: Order matters - more specific patterns should come first
            program_mappings = [
                (r'mtriage', 'phenix.mtriage'),
                (r'xtriage', 'phenix.xtriage'),
                (r'phaser', 'phenix.phaser'),
                (r'molecular\s+replacement', 'phenix.phaser'),
                (r'ligand\s*fit', 'phenix.ligandfit'),
                (r'fit.*ligand', 'phenix.ligandfit'),
                (r'refine', 'phenix.refine'),
                (r'autobuild', 'phenix.autobuild'),
                # map_to_model patterns (check BEFORE dock_in_map since "map" appears in both)
                (r'map\s*to\s*model', 'phenix.map_to_model'),
                (r'maptomodel', 'phenix.map_to_model'),
                (r'build.*model.*(?:into|in)\s*(?:the\s*)?map', 'phenix.map_to_model'),
                (r'(?:automated?\s+)?model\s+building.*(?:into|in)\s*(?:the\s*)?map', 'phenix.map_to_model'),
                # dock_in_map patterns (more specific than just "dock")
                (r'dock.*(?:in|into)\s*(?:the\s*)?map', 'phenix.dock_in_map'),
                (r'fit.*model.*(?:to|into)\s*(?:the\s*)?map', 'phenix.dock_in_map'),
                (r'map\s*symmetry', 'phenix.map_symmetry'),
                (r'(?:determin|find|detect).*symmetry', 'phenix.map_symmetry'),
                (r'symmetry.*(?:map|cryo)', 'phenix.map_symmetry'),
            ]
            for pattern, program in program_mappings:
                if re.search(pattern, stop_text, re.IGNORECASE):
                    if "after_program" not in directives["stop_conditions"]:
                        directives["stop_conditions"]["after_program"] = program
                    break

            # Handle density modification based on experiment type
            if "after_program" not in directives.get("stop_conditions", {}):
                if re.search(r'density\s+modif|denmod', stop_text, re.IGNORECASE):
                    if is_cryoem and not is_xray:
                        directives["stop_conditions"]["after_program"] = "phenix.resolve_cryo_em"
                    elif is_xray and not is_cryoem:
                        directives["stop_conditions"]["after_program"] = "phenix.autobuild_denmod"
                    else:
                        # Default to cryo-EM if ambiguous
                        directives["stop_conditions"]["after_program"] = "phenix.resolve_cryo_em"

    # Check for workflow continuation indicators that mean we should NOT set early stops
    # Words like "then", "later", "afterwards" indicate multi-step workflows
    continuation_indicators = [
        r'\bthen\b',
        r'\blater\b',
        r'\bafterwards?\b',
        r'\band\s+then\b',
        r'\bcontinue\b',
        r'\bnext\b',
        r'\bfollowed\s+by\b',
        r'\bafter\s+(?:that|this|which)\b',
        r'\bsubsequent',
    ]

    # Downstream task indicators - if user mentions these, don't stop early
    downstream_tasks = [
        r'ligand\s+(?:fit|fitting|place|placement)',
        r'fit\s+(?:the\s+)?ligand',
        r'add\s+(?:the\s+)?ligand',
        r'place\s+(?:the\s+)?ligand',
        r'validate|validation|molprobity',
        r'refine.*(?:further|more|additional)',
        r'additional\s+refine',
        r'water\s+(?:pick|add|place)',
        # Cryo-EM downstream tasks (indicate steps AFTER an earlier program)
        r'build\s+(?:a\s+)?model',
        r'model\s+building',
        r'map[\s_]to[\s_]model',
        r'apply\s+(?:ncs|symmetry)',
        r'generate\s+(?:the\s+)?(?:full|complete)',
        r'(?:full|complete)\s+(?:complex|assembly|oligomer|multimer)',
        r'(?:real[- ]space\s+)?refine.*(?:model|structure)',
        r'extract\s+(?:the\s+)?(?:unique|asymmetric)',
    ]

    has_continuation = any(re.search(p, advice_lower) for p in continuation_indicators)
    has_downstream_task = any(re.search(p, advice_lower) for p in downstream_tasks)

    # Also detect multi-step workflows by counting distinct program mentions
    # If multiple programs are referenced, this is a pipeline, not a single-step task
    multi_program_patterns = [
        (r'(?:auto[- ]?sharpen|map[\s_]sharpening|sharpen\s+(?:the\s+)?map)', 'sharpening'),
        (r'(?:map[\s_]symmetry|find\s+symmetry|determine\s+symmetry|ncs\s+(?:from|in))', 'symmetry'),
        (r'(?:map[\s_]to[\s_]model|build.*model.*(?:into|in)\s+(?:the\s+)?map)', 'build'),
        (r'(?:dock.*(?:in|into)\s+map|fit\s+model\s+to\s+map)', 'dock'),
        (r'(?:refine|refinement)', 'refine'),
        (r'(?:ligand\s*fit|fit.*ligand|place.*ligand)', 'ligand'),
        (r'(?:density\s+modif|resolve_cryo_em)', 'denmod'),
        (r'(?:apply\s+(?:ncs|symmetry)|generate.*(?:complex|assembly))', 'apply_ncs'),
    ]
    distinct_steps = set()
    for pattern, step_name in multi_program_patterns:
        if re.search(pattern, advice_lower, re.IGNORECASE):
            distinct_steps.add(step_name)
    has_multi_step = len(distinct_steps) >= 2

    # If there are continuation indicators, downstream tasks, or multiple distinct
    # program steps, this is a multi-step workflow - don't set after_program
    if has_continuation or has_downstream_task or has_multi_step:
        # This is a multi-step workflow - extract the first program to start with
        # Common patterns: "run X then Y", "calculate X and then refine"
        first_program_patterns = [
            (r'(?:calculate|run|compute)\s+(?:a\s+)?polder', 'phenix.polder'),
            (r'polder\s+(?:map|omit)', 'phenix.polder'),
            (r'(?:run|do)\s+xtriage', 'phenix.xtriage'),
            (r'(?:run|try)\s+phaser', 'phenix.phaser'),
            (r'(?:run|do)\s+mtriage', 'phenix.mtriage'),
            (r'dock.*(?:in|into)\s+(?:the\s+)?map', 'phenix.dock_in_map'),
            (r'map\s*to\s*model', 'phenix.map_to_model'),
        ]

        for pattern, program in first_program_patterns:
            if re.search(pattern, advice_lower, re.IGNORECASE):
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                # Set as start_with_program - workflow should run this first
                directives["stop_conditions"]["start_with_program"] = program
                break  # Only extract the first matching program
    else:
        # Detect tutorial/procedure patterns that imply stop after specific program
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

        for pattern, program in tutorial_patterns:
            if re.search(pattern, advice_lower, re.IGNORECASE):
                # This looks like a focused task - add stop condition
                if "stop_conditions" not in directives:
                    directives["stop_conditions"] = {}
                # Only set if not already specified
                if "after_program" not in directives.get("stop_conditions", {}):
                    directives["stop_conditions"]["after_program"] = program
                    directives["stop_conditions"]["skip_validation"] = True
                break

        # Handle density modification separately since it depends on experiment type
        if "after_program" not in directives.get("stop_conditions", {}):
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
                    break

    return directives
