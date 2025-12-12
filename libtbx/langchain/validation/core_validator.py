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


def fix_command_syntax(command: str, error_message: str, llm) -> str:
    """
    Ask the LLM to fix a command based on the specific error message.

    Args:
        command: The failed command
        error_message: The error message from validation
        llm: Language model to use for fixing

    Returns:
        str: The fixed command (or original if fixing fails)

    Example:
        fixed = fix_command_syntax(
            "phenix.refine model.pdb data.mtz twin_law=-h,-k,l",
            "Ambiguous parameter: twin_law",
            llm
        )
        # Returns: "phenix.refine model.pdb data.mtz refinement.main.twin_law=-h,-k,l"
    """
    fix_prompt = f"""You are a Phenix command-line expert. A command has a syntax error.

**Failed Command:**
{command}

**Error Message:**
{error_message}

**Your Task:**
Fix ONLY the syntax error. Do not change the logic or files.
Return ONLY the corrected command string, nothing else. No markdown, no explanation.

**Common Fixes:**
- "Ambiguous parameter" → Use full PHIL path (e.g., `xray_data.r_free_flags.generate=True` not `r_free_flags.generate=True`)
- "Ambiguous parameter" for twin_law → Use `refinement.main.twin_law=...`
- "not recognized" → Check spelling and PHIL hierarchy
- If the error mentions specific valid options, use one of those

**Critical:** Return ONLY the command, one line, no formatting.
"""

    try:
        response = llm.invoke(fix_prompt)
        fixed_command = response.content.strip()

        # Clean up any markdown formatting
        if "```" in fixed_command:
            lines = fixed_command.split('\n')
            for line in lines:
                if line and not line.startswith('```') and not line.startswith('#'):
                    fixed_command = line.strip()
                    break

        return fixed_command

    except Exception as e:
        return command  # Return original if fix fails

