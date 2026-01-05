"""
Core command validation functions.

This module handles:
- Syntax validation using --dry-run
- LLM-based command fixing
- General validation utilities

Usage:
    from libtbx.langchain.validation import validate_phenix_command, fix_command_syntax

    is_valid, error = validate_phenix_command(command)
    if not is_valid:
        fixed_command = fix_command_syntax(command, error, llm)
"""
from __future__ import absolute_import, division, print_function
from langchain_core.prompts import PromptTemplate

# Programs known to NOT support --dry-run (tested empirically)
PROGRAMS_WITHOUT_DRY_RUN = {
    'phenix.automr',
    'phenix.average_map_coeffs',
    'phenix.average_maps',
    'phenix.cc_star',
    'phenix.cif_as_pdb',
    'phenix.composite_omit_map',
    'phenix.data_viewer',
    'phenix.density_modification',
    'phenix.doc',
    'phenix.ensembler',
    'phenix.fab_elbow_angle',
    'phenix.feature_enhanced_map',
    'phenix.fem',
    'phenix.fft',
    'phenix.find_peaks_holes',
    'phenix.get_cc_ano',
    'phenix.get_cc_iso',
    'phenix.get_latest_version',
    'phenix.get_patterson_skew',
    'phenix.grow_density',
    'phenix.hyss',
    'phenix.iterative_ss_refine',
    'phenix.kinemage',
    'phenix.ligand_pipeline',
    'phenix.map_comparison',
    'phenix.merging_statistics',
    'phenix.model_idealization',
    'phenix.molprobity',
    'phenix.mrage',
    'phenix.mtz2map',
    'phenix.mtz_as_cif',
    'phenix.multistart_sa',
    'phenix.pdb.b_factor_stats',
    'phenix.pdb_interpretation',
    'phenix.phaser',
    'phenix.phaser_mp',
    'phenix.prime_and_switch_map',
    'phenix.probe',
    'phenix.reduce',
    'phenix.reindex',
    'phenix.rosetta.run_phenix_interface',
    'phenix.rosetta_refine',
    'phenix.sad_data_from_pdb',
    'phenix.sceds',
    'phenix.sculptor',
    'phenix.secondary_structure_restraints',
    'phenix.segment_and_split_map',
    'phenix.simple_ncs_from_pdb',
    'phenix.sisa',
    'phenix.sort_hetatms',
    'phenix.start_coot',
    'phenix.superpose_ligands',
    'phenix.validate_H',
    'phenix.xtriage',
}

def validate_is_phenix_command(command):
    """
    Ensure the command starts with a valid Phenix program name.

    Args:
        command: The command string to validate

    Returns:
        tuple: (is_valid: bool, error_message: str or None)
    """
    if not command:
        return False, "Empty command"

    command = command.strip()

    # Get the first word (program name)
    first_word = command.split()[0] if command.split() else ""

    # Must start with phenix.
    if not first_word.startswith('phenix.'):
        return False, f"Command does not start with 'phenix.': got '{first_word[:50]}...'" if len(first_word) > 50 else f"Command does not start with 'phenix.': got '{first_word}'"

    # Check for obvious garbage (non-ASCII characters in program name)
    try:
        first_word.encode('ascii')
    except UnicodeEncodeError:
        return False, f"Command contains non-ASCII characters: '{first_word[:30]}...'"

    return True, None

def validate_phenix_command(command: str) -> tuple:
    """
    Validates a Phenix command by running with --dry-run flag.
    This validates all parameters without actually running the job.

    Args:
        command: The Phenix command to validate

    Returns:
        tuple: (is_valid: bool, error_message: str)

    Example:
        is_valid, error = validate_phenix_command("phenix.refine model.pdb data.mtz")
        if not is_valid:
            print(f"Validation failed: {error}")
    """
    from libtbx import easy_run

    parts = command.split()
    if not parts:
        return False, "Empty command"

    program = parts[0]

    # Skip validation for programs that don't support --dry-run
    if program in PROGRAMS_WITHOUT_DRY_RUN:
        return True, ""

    # Extract parameters (anything with '=')
    params_to_test = []
    for part in parts[1:]:
        if '=' in part:
            params_to_test.append(part)

    if not params_to_test:
        return True, ""

    # Build test command with --dry-run
    # This will validate parameters without running the actual job
    test_command = f"{program} --dry-run {' '.join(params_to_test)}"

    try:
        result = easy_run.fully_buffered(test_command)

        # Check both stderr and stdout for errors
        stderr_text = "\n".join(result.stderr_lines)
        stdout_text = "\n".join(result.stdout_lines)
        all_output = stderr_text + stdout_text

        # If the error is about --dry-run itself, the program doesn't support it
        # In this case, skip validation (be permissive)
        dry_run_error_indicators = [
            'dry-run',
            'dry_run',
            '--dry-run',
        ]
        all_output_lower = all_output.lower()
        for indicator in dry_run_error_indicators:
            if indicator in all_output_lower:
                # Add to known list for future (in memory only)
                PROGRAMS_WITHOUT_DRY_RUN.add(program)
                return True, ""  # Skip validation, assume command is OK


        print(f"DEBUG VALIDATE: all_output_lower contains: {all_output_lower[:200]}")
        print(f"DEBUG VALIDATE: checking for 'an atomic model is required': {'an atomic model is required' in all_output_lower}")

        # Errors expected in dry-run validation (we don't pass input files)
        # These are NOT real syntax errors - they just mean the program needs files
        benign_dry_run_errors = [
            "an atomic model is required",
            "reflection file",
            "no data file",
            "no model file",
            "no pdb file",
            "no input file",
            "requires a model",
            "requires an mtz",
            "requires reflection",
        ]

        all_output_lower = all_output.lower()
        for benign in benign_dry_run_errors:
            if benign in all_output_lower:
                return True, ""  # Skip validation - texpected without files

        # Common error patterns
        error_patterns = [
            "Sorry:",
            "Error:",
            "ERROR:",
            "Ambiguous parameter",
            "not recognized",
            "Unknown parameter",
            "Invalid parameter",
            "Unrecognized PHIL"
        ]

        for pattern in error_patterns:
            if pattern in stderr_text or pattern in stdout_text:
                # Extract error context
                error_lines = []
                all_lines = result.stderr_lines + result.stdout_lines
                found_error = False
                for line in all_lines:
                    if pattern in line:
                        found_error = True
                    if found_error:
                        error_lines.append(line)
                        if len(error_lines) >= 15:  # Get more context
                            break

                error_msg = "\n".join(error_lines)
                return False, error_msg

        return True, ""

    except Exception as e:
        # If validation itself fails, be permissive
        return True, ""

def fix_ambiguous_parameters(command: str, error_message: str) -> str:
    """
    Mechanically fix common ambiguous parameter errors.
    """
    import re

    # Extract the ambiguous parameter from error
    match = re.search(r'Ambiguous parameter definition:\s*(\w+)\s*=\s*(\S+)', error_message)
    if not match:
        return command

    param_name = match.group(1)
    param_value = match.group(2)

    # Extract best matches from error message - try multiple patterns
    best_matches = []

    # Pattern 1: Lines that are just indented parameter names
    best_matches = re.findall(r'^\s{2,}(\S+\.[\w\.]+)\s*$', error_message, re.MULTILINE)

    # Pattern 2: After "Best matches:"
    if not best_matches:
        section_match = re.search(r'Best matches?:(.*?)(?:\n\n|\Z)', error_message, re.DOTALL)
        if section_match:
            best_matches = re.findall(r'(\S+\.[\w\.]+)', section_match.group(1))

    if not best_matches:
        return command

    # Use the first best match
    full_param = best_matches[0].strip()

    # Replace the ambiguous parameter with the full path
    pattern = rf'\b{re.escape(param_name)}\s*=\s*{re.escape(param_value)}'
    replacement = f'{full_param}={param_value}'
    new_command = re.sub(pattern, replacement, command)

    return new_command

def fix_command_syntax(
    command: str,
    error_message: str,
    rag_chain,
    attempt_number: int = 1
) -> str:
    """
    Uses an LLM to fix a command that failed validation.

    Args:
        command: The failing command
        error_message: The error from validation
        rag_chain: RAG chain for documentation lookup
        attempt_number: Which fix attempt this is

    Returns:
        Fixed command string
    """
    # NEW: Try mechanical fix for ambiguous parameters FIRST
    if "ambiguous parameter" in error_message.lower():
        fixed = fix_ambiguous_parameters(command, error_message)
        if fixed != command:
            return fixed

    # NEW: Try mechanical fix for unrecognized parameters
    if "not recognized" in error_message.lower():
        import re
        # Try to extract and remove the bad parameter
        patterns = [
            r"parameter.*?:\s*(\w+)\s*=",
            r"definition:\s*(\w+)\s*=",
        ]
        for pattern in patterns:
            match = re.search(pattern, error_message, re.IGNORECASE)
            if match:
                bad_param = match.group(1)
                stripped = re.sub(rf'\s*\b{bad_param}\s*=\s*\S+', '', command)
                stripped = ' '.join(stripped.split())
                if stripped != command:
                    return stripped


    program = command.split()[0] if command else ""

    # Get documentation for context
    docs_context = ""
    if rag_chain and program:
        try:
            program_name = program.replace("phenix.", "")
            docs_result = rag_chain.invoke(f"{program_name} command line syntax parameters")
            docs_context = docs_result.get("answer", "")[:2000]
        except Exception as e:
            pass

    template = """You are a Phenix command-line syntax expert.

A command failed validation. Fix it.
I will enclose the command in <failed_command> tags,
and the error message in <error_traceback> tags

**Failed Command:**
<failed_command>{command}</failed_command>

**Error Message:**
<error_traceback>{error_message}</error_traceback>

**Attempt Number:** {attempt_number}

**CRITICAL FIXING RULES:**

1. **"Ambiguous parameter" errors:**
   - The error shows "Best matches:" with specific parameter names
   - You MUST use one of those EXACT matches
   - Example: If error says "resolution = 3.5" is ambiguous with matches:
     - crystal_info.resolution
     - prediction.resolution
   - You should change `resolution=3.5` to `crystal_info.resolution=3.5`

2. **"not recognized" or "Unknown parameter" errors:**
   - The parameter does NOT exist for this program.
   - You MUST REMOVE the invalid parameter entirely.
   - Do NOT try to rename it or guess alternative names.

3. **Input file errors:**
   - Use positional arguments: `phenix.program model.pdb data.mtz`
   - Remove any `file_name=`, `model_file_name=` prefixes.

4. **NEVER use placeholders:**
   - Never use `/path/to/`, `/dev/null`, `model.pdb`, `data.mtz`
   - Use the actual filenames from the command

**Output:**
Return ONLY the fixed command. No explanation, no markdown.
"""

    # Create the template object
    prompt_template = PromptTemplate.from_template(template)

    # Generate the final string safely
    fix_prompt = prompt_template.format(
        command=command,
        error_message=error_message,
        attempt_number=attempt_number
    )


    try:
        from libtbx.langchain.model_setup import get_model
        llm = get_model(expensive=False)
        response = llm.invoke(fix_prompt)

        fixed = response.content if hasattr(response, 'content') else str(response)
        fixed = fixed.strip().strip('`').strip()

        # Remove any markdown formatting
        if fixed.startswith('bash') or fixed.startswith('shell'):
            fixed = fixed.split('\n', 1)[-1].strip()

        # Ensure it starts with the same program
        if program and not fixed.startswith(program):
            # Try to extract command from response
            for line in fixed.split('\n'):
                if line.strip().startswith(program):
                    fixed = line.strip()
                    break
            else:
                # If we can't find a valid command, return original
                return command

        return fixed

    except Exception as e:
        print(f"Error in fix_command_syntax: {e}")
        return command
