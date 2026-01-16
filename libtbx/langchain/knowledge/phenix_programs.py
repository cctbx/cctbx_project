"""
Phenix program discovery and information.

This module handles:
- Discovering available Phenix programs
- Getting program keywords via introspection

Usage:
    from libtbx.langchain.knowledge.phenix_programs import get_phenix_program_list, get_keywords_as_phil_string

    programs = get_phenix_program_list()
    keywords = get_keywords_as_phil_string('phenix.refine')
"""
from __future__ import absolute_import, division, print_function

import os


def get_phenix_program_list() -> list:
    """
    Returns a sorted list of all available Phenix programs.
    Scans the bin directory associated with the current environment.

    Returns:
        list: Sorted list of program names (e.g., ['phenix.xtriage', 'phenix.refine', ...])

    Example:
        programs = get_phenix_program_list()
        print(f"Found {len(programs)} Phenix programs")
    """
    # Attempt to find the bin directory
    import libtbx.load_env
    bin_dir = abs(libtbx.env.bin_path)
    if not os.path.isdir(bin_dir):
        # Fallback for some environments
        return []

    programs = []
    for filename in os.listdir(bin_dir):
        # Filter for likely Phenix/CCTBX programs
        if filename.find("development") > -1:
            continue
        if filename.startswith("phenix."):
            programs.append(filename)

    return sorted(programs)


def get_keywords_as_phil_string(program: str) -> str:
    """
    Get valid keywords for a Phenix program via introspection.

    Tries multiple strategies:
    1. --show_defaults (PHIL parameters)
    2. --help (help text)
    3. Run with no arguments (usage message)

    Args:
        program: Name of the Phenix program (e.g., 'phenix.refine')

    Returns:
        str: String containing valid keywords/parameters, or error message

    Example:
        keywords = get_keywords_as_phil_string('phenix.refine')
        print(keywords[:500])  # Print first 500 chars
    """
    from libtbx import easy_run

    # Import phenix_knowledge for usage hints
    try:
        from libtbx.langchain import phenix_knowledge as pk
    except ImportError:
        pk = None

    # Use the external Knowledge Base for hints
    hint_text = ""
    if pk and hasattr(pk, 'USAGE_HINTS') and program in pk.USAGE_HINTS:
        hint_text = f"\n{'-'*40}\n{pk.USAGE_HINTS[program]}\n{'-'*40}\n\n"

    print(f"Introspecting {program} to get valid keywords...")

    # Strategy 1: PHIL Defaults
    cmd_defaults = f"{program} --show_defaults 3"
    result = easy_run.fully_buffered(cmd_defaults)

    if result.return_code == 0 and len(result.stdout_lines) > 10 and \
       "unrecognized argument" not in "\n".join(result.stderr_lines):
        return hint_text + "\n".join(result.stdout_lines)

    # Strategy 2: Help Text
    print(f"  --show_defaults failed. Trying --help...")
    cmd_help = f"{program} --help"
    result_help = easy_run.fully_buffered(cmd_help)
    help_text = "\n".join(result_help.stdout_lines + result_help.stderr_lines)
    if len(help_text) > 100:
        return hint_text + help_text

    # Strategy 3: Naked Call
    print(f"  --help failed. Trying run with no arguments...")
    result_none = easy_run.fully_buffered(f"{program}")
    none_text = "\n".join(result_none.stdout_lines + result_none.stderr_lines)
    if len(none_text) > 100:
        return hint_text + none_text

    return f"ERROR: Could not extract keywords for {program}. It might not be in the path."
