"""
Agent planning utilities.

Retained functions:
  - fix_program_parameters(): Fix common parameter mistakes using JSON config
  - extract_output_files(): Parse output files from summary text

Historical note: This module previously contained generate_next_move(),
construct_command_mechanically(), and other LLM-based planning functions.
Those were superseded by the YAML-driven CommandBuilder + rules_selector
pipeline and removed in v112.10 (Fix 7).
"""
from __future__ import absolute_import, division, print_function

import re
import os
import json


# =============================================================================
# Parameter Fixes
# =============================================================================
_PARAMETER_FIXES = None


def load_parameter_fixes():
    """Load program-specific parameter fixes from JSON.

    The file lives in knowledge/ (PHENIX domain config), not agent/.
    """
    fixes = {}

    # Look in knowledge/ directory (sibling of agent/)
    agent_dir = os.path.dirname(__file__)
    knowledge_dir = os.path.join(os.path.dirname(agent_dir), 'knowledge')
    path = os.path.join(knowledge_dir, 'parameter_fixes.json')

    # Fallback: try libtbx path for installed PHENIX
    if not os.path.exists(path):
        try:
            from libtbx.langchain.knowledge import __file__ as kfile
            path = os.path.join(os.path.dirname(kfile), 'parameter_fixes.json')
        except Exception:
            pass

    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                fixes = json.load(f)
            print(f"Loaded parameter fixes from {path}")
        except Exception as e:
            print(f"Warning: Could not load parameter fixes from {path}: {e}")
    else:
        print(f"Warning: parameter_fixes.json not found")

    return fixes


def get_parameter_fixes():
    """Get parameter fixes, loading from file if needed."""
    global _PARAMETER_FIXES
    if _PARAMETER_FIXES is None:
        _PARAMETER_FIXES = load_parameter_fixes()
    return _PARAMETER_FIXES


def fix_program_parameters(command, program):
    """
    Fix common parameter mistakes for a specific program.

    Uses a JSON config that maps wrong parameter names to correct ones.

    Args:
        command: The generated command string
        program: The program name (e.g., 'phenix.ligandfit')

    Returns:
        Fixed command string
    """
    fixes = get_parameter_fixes()

    if program not in fixes:
        return command

    program_fixes = fixes[program]

    for wrong_param, right_param in program_fixes.items():
        # Skip comment keys
        if wrong_param.startswith('_'):
            continue

        # Skip if right_param is already in the command (avoid double-fixing)
        if right_param and right_param + '=' in command:
            continue

        # Pattern to match parameter=value (handles quoted values too)
        # Use word boundary and negative lookbehind to avoid matching if already prefixed
        # e.g., don't match "input_labels" if "file_info.input_labels" exists
        pattern = rf'(?<![.\w]){re.escape(wrong_param)}=("(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'|\S+)'

        if right_param == "":
            # Remove the parameter name but keep the value
            command = re.sub(pattern, r'\1', command)
        elif right_param is None:
            # Remove entirely
            command = re.sub(pattern, '', command)
        else:
            # Replace with correct parameter name
            command = re.sub(pattern, rf'{right_param}=\1', command)

    # Clean up any double spaces
    command = ' '.join(command.split())

    return command


# =============================================================================
# Output File Extraction
# =============================================================================


def extract_output_files(summary_text):
    """Parses the 'Key Output Files' section from a summary."""
    if not summary_text:
        return []
    files = []
    try:
        if "**Key Output Files:**" in summary_text:
            # Grab text after header
            post = summary_text.split("**Key Output Files:**")[1]
            # Stop at next header or end
            content = post.split("**")[0].strip()
            # Extract filenames (pattern includes paths with /)
            files = re.findall(r'[\w\-\.\/]+\.\w{2,4}', content)
    except Exception as e:
        pass
    return files
