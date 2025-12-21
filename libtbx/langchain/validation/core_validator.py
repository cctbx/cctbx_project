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

    fix_prompt = f"""You are a Phenix command-line syntax expert.

A command failed validation. Fix it.

**Failed Command:**
{command}

**Error Message:**
{error_message}

**Attempt Number:** {attempt_number}

**Documentation Context:**
{docs_context}

**CRITICAL FIXING RULES:**

1. **"not recognized" or "Unknown parameter" errors:**
   - The parameter does NOT exist for this program.
   - You MUST REMOVE the invalid parameter entirely.
   - Do NOT try to rename it or guess alternative names.
   - Example: If `nproc=4` is not recognized, just remove it. Don't try `n_processors=4`.

2. **"Ambiguous parameter" errors:**
   - The parameter name is incomplete or matches multiple options.
   - If the error suggests specific matches, use the FIRST match that makes sense.
   - If no good match, REMOVE the parameter.

3. **"atomic model is required" or "reflection file" errors:**
   - These are benign dry-run errors - the command syntax is likely correct.
   - Return the command unchanged.

4. **Input file syntax errors:**
   - Use positional arguments: `phenix.program model.pdb data.mtz`
   - Remove any `file_name=`, `model_file_name=`, `reflection_file_name=` prefixes.

5. **General principle:** When in doubt, REMOVE the problematic parameter rather than guessing.

**Output:**
Return ONLY the fixed command. No explanation, no markdown.
"""

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
