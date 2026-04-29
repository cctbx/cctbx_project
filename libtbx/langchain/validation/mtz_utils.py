"""
MTZ file utilities for validation.

This module handles:
- Checking for R-free flags in MTZ files
- Adding R-free generation parameter when needed

Usage:
    from libtbx.langchain.validation.mtz_utils import mtz_has_rfree_flags

    if not mtz_has_rfree_flags('data.mtz'):
        print("MTZ file lacks R-free flags")
"""
from __future__ import absolute_import, division, print_function


def mtz_has_rfree_flags(mtz_file: str) -> bool:
    """
    Check if an MTZ file contains R-free flags.

    Args:
        mtz_file: Path to the MTZ file

    Returns:
        bool: True if R-free flags are present, False otherwise

    Example:
        if mtz_has_rfree_flags('data.mtz'):
            print("R-free flags found")
        else:
            print("Need to generate R-free flags")
    """
    from libtbx import easy_run

    try:
        # Use phenix.mtz.dump to check contents
        result = easy_run.fully_buffered(f"phenix.mtz.dump {mtz_file}")
        output = "\n".join(result.stdout_lines)

        # Look for R-free-flags column
        # Common names: R-free-flags, FreeR_flag, FREE, RFREE
        rfree_indicators = [
            'R-free-flags',
            'FreeR_flag',
            'FREE',
            'RFREE',
            'R_FREE_FLAG',
            'test_flag'
        ]

        output_lower = output.lower()
        for indicator in rfree_indicators:
            if indicator.lower() in output_lower:
                return True

        return False

    except Exception as e:
        # If we can't check, assume no flags (safer)
        return False


def add_rfree_generation_if_needed(command: str, program: str) -> str:
    """
    For phenix.refine, check if the MTZ file has R-free flags.
    If not, automatically add xray_data.r_free_flags.generate=True

    Args:
        command: The command to check/modify
        program: The program name

    Returns:
        str: The command, possibly with R-free generation added

    Example:
        command = add_rfree_generation_if_needed(
            "phenix.refine model.pdb data.mtz",
            "phenix.refine"
        )
        # If data.mtz lacks R-free flags, returns:
        # "phenix.refine model.pdb data.mtz xray_data.r_free_flags.generate=True"
    """
    if 'phenix.refine' not in program:
        return command

    # Check if command already has r_free_flags parameter
    if 'r_free_flags.generate' in command:
        return command  # Already handled

    # Extract MTZ file from command
    parts = command.split()
    mtz_file = None
    for part in parts:
        if part.endswith('.mtz') and not part.startswith('-'):
            mtz_file = part
            break

    if not mtz_file:
        return command  # No MTZ file found

    # Check if MTZ has R-free flags
    if not mtz_has_rfree_flags(mtz_file):
        # Add the parameter
        command += " xray_data.r_free_flags.generate=True"

    return command
