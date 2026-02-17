#!/usr/bin/env python
"""
YAML Tools for PHENIX AI Agent Configuration.

This script provides tools for working with the YAML configuration files:
1. Validate YAML syntax and structure
2. Display formatted contents
3. Compare two sets of YAML files to show differences

Usage:
    # Validate and display all YAML files
    python yaml_tools.py validate

    # Validate a specific file
    python yaml_tools.py validate programs.yaml

    # Compare two directories
    python yaml_tools.py compare /path/to/old /path/to/new

    # Compare specific files
    python yaml_tools.py compare old/programs.yaml new/programs.yaml

    # Display formatted contents
    python yaml_tools.py display programs.yaml

    # Show summary of all files
    python yaml_tools.py summary
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import argparse

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Try to import yaml
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Install with: pip install PyYAML")
    sys.exit(1)


# =============================================================================
# COLORS AND FORMATTING
# =============================================================================

class Colors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    RESET = '\033[0m'

    @classmethod
    def disable(cls):
        """Disable colors (for non-terminal output)."""
        cls.HEADER = ''
        cls.BLUE = ''
        cls.CYAN = ''
        cls.GREEN = ''
        cls.YELLOW = ''
        cls.RED = ''
        cls.BOLD = ''
        cls.DIM = ''
        cls.RESET = ''


def colored(text, color):
    """Apply color to text."""
    return color + str(text) + Colors.RESET


def print_header(text, char='='):
    """Print a header line."""
    width = min(80, len(text) + 4)
    print()
    print(colored(char * width, Colors.BOLD))
    print(colored("  " + text, Colors.BOLD))
    print(colored(char * width, Colors.BOLD))


def print_subheader(text):
    """Print a subheader."""
    print()
    print(colored(">>> " + text, Colors.CYAN + Colors.BOLD))
    print(colored("-" * (len(text) + 4), Colors.CYAN))


def print_success(text):
    """Print success message."""
    print(colored("✓ " + text, Colors.GREEN))


def print_error(text):
    """Print error message."""
    print(colored("✗ " + text, Colors.RED))


def print_warning(text):
    """Print warning message."""
    print(colored("⚠ " + text, Colors.YELLOW))


def print_info(text):
    """Print info message."""
    print(colored("  " + text, Colors.DIM))


# =============================================================================
# YAML LOADING AND VALIDATION
# =============================================================================

def find_yaml_files(directory):
    """Find all YAML files in a directory.

    Searches for YAML files in:
    1. <directory>/knowledge/ (if exists)
    2. <directory>/../knowledge/ (sibling directory, if exists)
    3. <directory> itself
    """
    yaml_files = []

    # Try subdirectory first
    knowledge_dir = os.path.join(directory, "knowledge")
    if os.path.isdir(knowledge_dir):
        directory = knowledge_dir
    else:
        # Try sibling directory (for when script is in agent/ but YAML is in knowledge/)
        parent_dir = os.path.dirname(directory)
        sibling_knowledge = os.path.join(parent_dir, "knowledge")
        if os.path.isdir(sibling_knowledge):
            directory = sibling_knowledge

    if os.path.isdir(directory):
        for f in os.listdir(directory):
            if f.endswith(('.yaml', '.yml')):
                yaml_files.append(os.path.join(directory, f))

    return sorted(yaml_files)


def list_yaml_files(paths=None):
    """
    List all available YAML files with descriptions.

    Args:
        paths: Optional list of directories to search. If empty, uses default location.
    """
    # Determine directory to search
    if paths:
        search_dirs = paths
    else:
        search_dirs = [script_dir]

    print_header("Available YAML Configuration Files")
    print()

    # File descriptions
    file_descriptions = {
        'programs.yaml': 'Program definitions (inputs, outputs, commands, invariants)',
        'workflows.yaml': 'Workflow state machines (phases, transitions, conditions)',
        'metrics.yaml': 'Quality metrics and thresholds',
        'file_categories.yaml': 'File categorization rules (extensions, patterns)',
    }

    total_files = 0

    for search_dir in search_dirs:
        yaml_files = find_yaml_files(search_dir)

        if not yaml_files:
            print(colored("  No YAML files found in: %s" % search_dir, Colors.YELLOW))
            continue

        # Show the actual directory where files were found
        if yaml_files:
            display_dir = os.path.dirname(yaml_files[0])
        else:
            display_dir = search_dir

        print(colored("  Directory: %s" % display_dir, Colors.DIM))
        print()

        for filepath in yaml_files:
            filename = os.path.basename(filepath)
            filesize = os.path.getsize(filepath)

            # Get description
            desc = file_descriptions.get(filename, "")

            # Try to get entry count from file
            data, error = load_yaml_file(filepath)
            entry_count = ""
            if not error and isinstance(data, dict):
                count = len(data)
                entry_count = "(%d entries)" % count

            # Format size
            if filesize < 1024:
                size_str = "%d B" % filesize
            else:
                size_str = "%.1f KB" % (filesize / 1024.0)

            print("  %s" % colored(filename, Colors.CYAN + Colors.BOLD))
            print("    Path: %s" % filepath)
            print("    Size: %s %s" % (size_str, entry_count))
            if desc:
                print("    %s" % colored(desc, Colors.DIM))
            print()

            total_files += 1

    print(colored("  Total: %d YAML file(s)" % total_files, Colors.GREEN))
    print()
    print(colored("  Use 'yaml_tools.py display <file>' to view contents", Colors.DIM))
    print(colored("  Use 'yaml_tools.py validate' to check for errors", Colors.DIM))
    print(colored("  Use 'yaml_tools.py terms' to see all defined terms", Colors.DIM))


def load_yaml_file(filepath):
    """
    Load a YAML file and return (data, error).

    Returns:
        tuple: (data, None) on success, (None, error_message) on failure
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()

        # Try to parse
        data = yaml.safe_load(content)
        return data, None

    except yaml.YAMLError as e:
        # Extract line/column info
        error_msg = str(e)
        if hasattr(e, 'problem_mark'):
            mark = e.problem_mark
            error_msg = "Line %d, Column %d: %s" % (
                mark.line + 1, mark.column + 1, e.problem or "syntax error"
            )
        return None, error_msg

    except IOError as e:
        return None, "Cannot read file: %s" % e


def validate_yaml_structure(data, filepath):
    """
    Validate the structure of a YAML file based on its type.

    Returns:
        list: List of (level, message) tuples where level is 'error', 'warning', or 'info'
    """
    issues = []
    filename = os.path.basename(filepath).lower()

    # Detect file type from name (allow variations like programs_test.yaml)
    if 'program' in filename:
        issues.extend(_validate_programs(data))
    elif 'workflow' in filename:
        issues.extend(_validate_workflows(data))
    elif 'metric' in filename:
        issues.extend(_validate_metrics(data))
    elif 'file_categor' in filename or 'categories' in filename:
        issues.extend(_validate_file_categories(data))
    else:
        # Try to auto-detect based on content
        if any(k.startswith('phenix.') for k in data.keys()):
            issues.extend(_validate_programs(data))
        elif any('phases' in v for v in data.values() if isinstance(v, dict)):
            issues.extend(_validate_workflows(data))
        elif any('direction' in v or 'thresholds' in v for v in data.values() if isinstance(v, dict)):
            issues.extend(_validate_metrics(data))
        elif any('extensions' in v or 'patterns' in v for v in data.values() if isinstance(v, dict)):
            issues.extend(_validate_file_categories(data))

    return issues


def _validate_file_categories(data):
    """Validate file_categories.yaml structure."""
    issues = []

    if not isinstance(data, dict):
        issues.append(('error', "Root must be a dictionary"))
        return issues

    # Valid fields for a category
    valid_category_fields = {
        'description', 'extensions', 'patterns', 'excludes',
        'subcategory_of', 'also_in', 'notes', 'max_basename_length'
    }

    # Track category names for cross-reference validation
    category_names = set(data.keys())

    for cat_name, cat_def in data.items():
        if not isinstance(cat_def, dict):
            issues.append(('error', "Category '%s' must be a dictionary" % cat_name))
            continue

        # Check for unknown fields
        for field in cat_def.keys():
            if field not in valid_category_fields:
                issues.append(('warning', "Category '%s' has unknown field '%s'" % (cat_name, field)))

        # Check required fields
        if 'description' not in cat_def:
            issues.append(('warning', "Category '%s' missing 'description'" % cat_name))

        # Validate extensions are a list
        extensions = cat_def.get('extensions', [])
        if extensions and not isinstance(extensions, list):
            issues.append(('error', "Category '%s' extensions must be a list" % cat_name))

        # Validate patterns are a list
        patterns = cat_def.get('patterns', [])
        if patterns and not isinstance(patterns, list):
            issues.append(('error', "Category '%s' patterns must be a list" % cat_name))

        # Validate excludes are a list
        excludes = cat_def.get('excludes', [])
        if excludes and not isinstance(excludes, list):
            issues.append(('error', "Category '%s' excludes must be a list" % cat_name))

        # Validate subcategory_of reference exists
        parent = cat_def.get('subcategory_of')
        if parent and parent not in category_names:
            issues.append(('warning', "Category '%s' references unknown parent '%s'" % (cat_name, parent)))

        # Validate also_in references exist
        also_in = cat_def.get('also_in', [])
        if also_in:
            if not isinstance(also_in, list):
                issues.append(('error', "Category '%s' also_in must be a list" % cat_name))
            else:
                for ref in also_in:
                    if ref not in category_names:
                        issues.append(('warning', "Category '%s' also_in references unknown '%s'" % (cat_name, ref)))

    return issues

    return issues


def _validate_programs(data):
    """Validate programs.yaml structure."""
    issues = []

    if not isinstance(data, dict):
        issues.append(('error', "Root must be a dictionary"))
        return issues

    # Programs can be under 'programs' key or at root level
    programs = data.get('programs', {})
    if not programs:
        # Check if programs are at root level (keys starting with 'phenix.')
        programs = {k: v for k, v in data.items()
                   if isinstance(v, dict) and (k.startswith('phenix.') or '.' in k)}

    if not programs:
        issues.append(('warning', "No programs defined"))
        return issues

    # Valid top-level fields for a program
    valid_program_fields = {
        'description', 'category', 'experiment_type', 'experiment_types',
        'inputs', 'outputs', 'command', 'log_parsing', 'strategy_options',
        'requires', 'produces', 'hint', 'notes', 'defaults', 'strategy_flags',
        'hints',
        # New fields from steps 1-3
        'invariants',           # Step 1: Auto-fix rules
        'input_priorities',     # Step 2: File category priorities
        'user_advice_keywords', # Step 3: Keywords for user advice matching
        # Additional program-level fields
        'done_tracking',        # Done flag configuration (flag name, run_once, etc.)
        'gui_app_id',           # wxGUI2 app_id for project History (fallback)
        'gui_app_id_cryoem',    # Cryo-EM variant app_id (if different from gui_app_id)
        'stop_directive_patterns',  # Regex patterns for stop-condition parsing
        'requires_full_map',    # Program needs full map, not half-maps
    }

    # Valid fields for input definitions
    valid_input_fields = {
        'extensions', 'flag', 'required', 'optional', 'type', 'description',
        'default', 'pattern', 'multiple', 'exclude_patterns', 'priority_patterns',
        'prefer_patterns', 'from_session'
    }

    # Valid fields for output definitions
    valid_output_fields = {
        'pattern', 'type', 'description', 'multiple'
    }

    # Valid fields for log_parsing
    valid_log_parsing_fields = {
        'pattern', 'type', 'group', 'default', 'fallback_pattern', 'extract',
        # Display/formatting fields used by summary_display.py
        'display_name', 'summary_format',
        # No-match handling fields used by metric_patterns.py
        'no_match_pattern', 'no_match_value'
    }

    # Valid fields for invariants (Step 1)
    valid_invariant_fields = {
        'name', 'description', 'check', 'fix', 'message'
    }

    # Valid fields for invariant checks
    valid_check_fields = {
        'has_file', 'has_strategy', 'strategy_equals', 'any_of', 'all_of',
        'has_file_category', 'not_has_file_category_only'  # File category checks
    }

    # Valid fields for invariant fixes
    valid_fix_fields = {
        'set_strategy', 'set_file', 'remove_file', 'auto_fill_resolution',
        'auto_fill_output_prefix'  # Auto-generate output prefix
    }

    # Valid fields for input_priorities (Step 2)
    valid_input_priority_fields = {
        'categories', 'exclude_categories', 'reason'
    }

    for prog_name, prog_def in programs.items():
        if not isinstance(prog_def, dict):
            issues.append(('error', "Program '%s' must be a dictionary" % prog_name))
            continue

        # Check for unknown top-level fields
        for field in prog_def.keys():
            if field not in valid_program_fields:
                issues.append(('warning', "Program '%s' has unknown field '%s' (typo?)" % (prog_name, field)))

        # Check required fields
        if 'description' not in prog_def:
            issues.append(('warning', "Program '%s' missing 'description'" % prog_name))

        # Validate inputs structure
        inputs = prog_def.get('inputs', {})
        if inputs and isinstance(inputs, dict):
            for input_group, input_defs in inputs.items():
                if input_group not in ('required', 'optional'):
                    # It's an input name directly, check its fields
                    if isinstance(input_defs, dict):
                        for field in input_defs.keys():
                            if field not in valid_input_fields:
                                issues.append(('warning', "Input '%s.%s' has unknown field '%s' (typo?)" % (
                                    prog_name, input_group, field)))
                elif isinstance(input_defs, dict):
                    # It's a required/optional group
                    for input_name, input_def in input_defs.items():
                        if isinstance(input_def, dict):
                            for field in input_def.keys():
                                if field not in valid_input_fields:
                                    issues.append(('warning', "Input '%s.%s.%s' has unknown field '%s' (typo?)" % (
                                        prog_name, input_group, input_name, field)))

        # Validate outputs structure
        outputs = prog_def.get('outputs', {})
        if outputs and isinstance(outputs, dict):
            for out_type, out_defs in outputs.items():
                if out_type not in ('files', 'metrics'):
                    issues.append(('info', "Program '%s' has unusual output type '%s'" % (prog_name, out_type)))
                if isinstance(out_defs, list):
                    for out_def in out_defs:
                        if isinstance(out_def, dict):
                            for field in out_def.keys():
                                if field not in valid_output_fields:
                                    issues.append(('warning', "Output in '%s.%s' has unknown field '%s' (typo?)" % (
                                        prog_name, out_type, field)))

        # Validate log_parsing structure
        log_parsing = prog_def.get('log_parsing', {})
        if log_parsing and isinstance(log_parsing, dict):
            for metric_name, metric_def in log_parsing.items():
                if isinstance(metric_def, dict):
                    for field in metric_def.keys():
                        if field not in valid_log_parsing_fields:
                            issues.append(('warning', "Log parsing '%s.%s' has unknown field '%s' (typo?)" % (
                                prog_name, metric_name, field)))

        # Validate invariants structure (Step 1)
        invariants = prog_def.get('invariants', [])
        if invariants:
            if not isinstance(invariants, list):
                issues.append(('error', "Program '%s' invariants must be a list" % prog_name))
            else:
                for i, inv in enumerate(invariants):
                    if not isinstance(inv, dict):
                        issues.append(('error', "Program '%s' invariant %d must be a dict" % (prog_name, i)))
                        continue

                    # Check invariant fields
                    for field in inv.keys():
                        if field not in valid_invariant_fields:
                            issues.append(('warning', "Invariant '%s[%d]' has unknown field '%s'" % (
                                prog_name, i, field)))

                    # Check required invariant fields
                    if 'name' not in inv:
                        issues.append(('warning', "Invariant '%s[%d]' missing 'name'" % (prog_name, i)))
                    if 'check' not in inv:
                        issues.append(('warning', "Invariant '%s[%d]' missing 'check'" % (prog_name, i)))

                    # Validate check structure
                    check = inv.get('check', {})
                    if isinstance(check, dict):
                        _validate_check_recursive(check, valid_check_fields, prog_name, i, issues)

                    # Validate fix structure
                    fix = inv.get('fix', {})
                    if isinstance(fix, dict):
                        for field in fix.keys():
                            if field not in valid_fix_fields:
                                issues.append(('warning', "Invariant fix '%s[%d]' has unknown field '%s'" % (
                                    prog_name, i, field)))

        # Validate input_priorities structure (Step 2)
        input_priorities = prog_def.get('input_priorities', {})
        if input_priorities:
            if not isinstance(input_priorities, dict):
                issues.append(('error', "Program '%s' input_priorities must be a dict" % prog_name))
            else:
                for input_name, priority_def in input_priorities.items():
                    if not isinstance(priority_def, dict):
                        issues.append(('error', "Input priority '%s.%s' must be a dict" % (
                            prog_name, input_name)))
                        continue

                    for field in priority_def.keys():
                        if field not in valid_input_priority_fields:
                            issues.append(('warning', "Input priority '%s.%s' has unknown field '%s'" % (
                                prog_name, input_name, field)))

                    # Validate categories is a list
                    categories = priority_def.get('categories', [])
                    if categories and not isinstance(categories, list):
                        issues.append(('error', "Input priority '%s.%s' categories must be a list" % (
                            prog_name, input_name)))

                    exclude = priority_def.get('exclude_categories', [])
                    if exclude and not isinstance(exclude, list):
                        issues.append(('error', "Input priority '%s.%s' exclude_categories must be a list" % (
                            prog_name, input_name)))

        # Validate user_advice_keywords (Step 3)
        keywords = prog_def.get('user_advice_keywords', [])
        if keywords:
            if not isinstance(keywords, list):
                issues.append(('error', "Program '%s' user_advice_keywords must be a list" % prog_name))
            else:
                for kw in keywords:
                    if not isinstance(kw, str):
                        issues.append(('warning', "Keyword in '%s' should be a string, got %s" % (
                            prog_name, type(kw).__name__)))

        # Check command template
        command = prog_def.get('command', {})
        if command and isinstance(command, dict):
            if 'base' not in command:
                issues.append(('warning', "Program '%s' command has no base" % prog_name))

    return issues


def _validate_check_recursive(check, valid_fields, prog_name, inv_idx, issues):
    """Recursively validate check structure (handles any_of, all_of)."""
    for field in check.keys():
        if field not in valid_fields:
            issues.append(('warning', "Check in '%s[%d]' has unknown field '%s'" % (
                prog_name, inv_idx, field)))

        # Recurse into any_of / all_of
        if field in ('any_of', 'all_of'):
            sub_checks = check[field]
            if isinstance(sub_checks, list):
                for sub_check in sub_checks:
                    if isinstance(sub_check, dict):
                        _validate_check_recursive(sub_check, valid_fields, prog_name, inv_idx, issues)


def _validate_workflows(data):
    """Validate workflows.yaml structure."""
    issues = []

    if not isinstance(data, dict):
        issues.append(('error', "Root must be a dictionary"))
        return issues

    # Workflows can be under 'workflows' key or at root level
    workflows = data.get('workflows', {})
    if not workflows:
        # Check for workflow-like keys at root level
        workflows = {k: v for k, v in data.items()
                    if isinstance(v, dict) and 'phases' in v}

    if not workflows:
        issues.append(('warning', "No workflows defined"))
        return issues

    # Valid top-level fields for a workflow
    valid_workflow_fields = {
        'description', 'phases', 'targets', 'thresholds', 'automation_paths'
    }

    # Valid fields for a phase
    valid_phase_fields = {
        'description', 'goal', 'programs', 'transitions', 'extracts',
        'repeat', 'stop', 'stop_reasons', 'conditions', 'requires', 'optional'
    }

    # Valid fields for a program entry in a phase
    valid_program_entry_fields = {
        'program', 'preferred', 'conditions', 'strategy', 'hint',
        'required', 'optional', 'requires_resolution', 'strategy_depends_on',
        'priority_when',  # Step 3: Condition for priority selection
    }

    # Valid priority_when conditions
    valid_priority_when_conditions = {
        'strong_anomalous',  # Strong anomalous signal detected
        # Add more conditions here as needed
    }

    # Valid fields for transitions
    valid_transition_fields = {
        'on_complete', 'on_target_reached', 'on_plateau', 'on_max_cycles',
        'on_ligandfit', 'if_predict_only', 'if_quality_acceptable',
        'if_quality_poor', 'if_has_full_map', 'if_needs_optimization', 'else'
    }

    for wf_name, wf_def in workflows.items():
        if not isinstance(wf_def, dict):
            issues.append(('error', "Workflow '%s' must be a dictionary" % wf_name))
            continue

        # Check for unknown workflow fields
        for field in wf_def.keys():
            if field not in valid_workflow_fields:
                issues.append(('warning', "Workflow '%s' has unknown field '%s' (typo?)" % (wf_name, field)))

        # Check for phases
        phases = wf_def.get('phases', {})
        if not phases:
            issues.append(('warning', "Workflow '%s' has no phases" % wf_name))
            continue

        # Check each phase
        for phase_name, phase_def in phases.items():
            if not isinstance(phase_def, dict):
                issues.append(('error', "Phase '%s.%s' must be a dictionary" % (wf_name, phase_name)))
                continue

            # Check for unknown phase fields
            for field in phase_def.keys():
                if field not in valid_phase_fields:
                    issues.append(('warning', "Phase '%s.%s' has unknown field '%s' (typo?)" % (
                        wf_name, phase_name, field)))

            # Check programs list
            programs = phase_def.get('programs', [])
            if isinstance(programs, list):
                for prog_entry in programs:
                    if isinstance(prog_entry, dict):
                        for field in prog_entry.keys():
                            if field not in valid_program_entry_fields:
                                issues.append(('warning', "Program entry in '%s.%s' has unknown field '%s' (typo?)" % (
                                    wf_name, phase_name, field)))

                        # Validate priority_when value
                        priority_when = prog_entry.get('priority_when')
                        if priority_when:
                            if priority_when not in valid_priority_when_conditions:
                                issues.append(('warning', "Program '%s' in '%s.%s' has unknown priority_when '%s'" % (
                                    prog_entry.get('program', '?'), wf_name, phase_name, priority_when)))

            # Check transitions
            transitions = phase_def.get('transitions', {})
            if isinstance(transitions, dict):
                for field in transitions.keys():
                    if field not in valid_transition_fields:
                        issues.append(('warning', "Transition in '%s.%s' has unknown field '%s' (typo?)" % (
                            wf_name, phase_name, field)))

            if 'programs' not in phase_def and not phase_def.get('stop'):
                issues.append(('warning', "Phase '%s.%s' has no programs" % (wf_name, phase_name)))

    return issues


def _validate_metrics(data):
    """Validate metrics.yaml structure."""
    issues = []

    if not isinstance(data, dict):
        issues.append(('error', "Root must be a dictionary"))
        return issues

    # Metrics can be under 'metrics' key or at root level
    metrics = data.get('metrics', {})
    if not metrics:
        # Check for metric-like keys (those with thresholds or direction)
        metrics = {k: v for k, v in data.items()
                  if isinstance(v, dict) and ('thresholds' in v or 'direction' in v)}

    if not metrics:
        issues.append(('warning', "No metrics defined"))
        return issues

    # Valid top-level fields for a metric
    valid_metric_fields = {
        'description', 'direction', 'thresholds', 'resolution_dependent',
        'unit', 'format', 'display_name', 'category', 'notes',
        'type', 'by_resolution', 'extraction', 'interpretation',  # Additional fields used in metrics.yaml
    }

    # Valid fields for thresholds
    valid_threshold_fields = {
        'good', 'acceptable', 'poor', 'target', 'warning', 'critical'
    }

    # Valid fields for resolution_dependent entries
    valid_resolution_fields = {
        'min', 'max', 'range', 'good', 'acceptable', 'poor', 'target'
    }

    # Valid directions
    valid_directions = {'minimize', 'maximize'}

    for metric_name, metric_def in metrics.items():
        if not isinstance(metric_def, dict):
            issues.append(('error', "Metric '%s' must be a dictionary" % metric_name))
            continue

        # Check for unknown metric fields
        for field in metric_def.keys():
            if field not in valid_metric_fields:
                issues.append(('warning', "Metric '%s' has unknown field '%s' (typo?)" % (metric_name, field)))

        # Check direction
        direction = metric_def.get('direction')
        if direction and direction not in valid_directions:
            issues.append(('warning', "Metric '%s' has invalid direction '%s' (should be 'minimize' or 'maximize')" % (
                metric_name, direction)))

        # Check thresholds structure
        thresholds = metric_def.get('thresholds', {})
        if isinstance(thresholds, dict):
            for field in thresholds.keys():
                if field not in valid_threshold_fields:
                    issues.append(('warning', "Threshold in '%s' has unknown field '%s' (typo?)" % (
                        metric_name, field)))

        # Check resolution_dependent structure
        res_dep = metric_def.get('resolution_dependent', [])
        if isinstance(res_dep, list):
            for entry in res_dep:
                if isinstance(entry, dict):
                    for field in entry.keys():
                        if field not in valid_resolution_fields:
                            issues.append(('warning', "Resolution entry in '%s' has unknown field '%s' (typo?)" % (
                                metric_name, field)))

        # Check for missing thresholds
        if 'thresholds' not in metric_def and 'resolution_dependent' not in metric_def:
            issues.append(('info', "Metric '%s' has no thresholds defined" % metric_name))

    return issues


# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================

def display_programs(data):
    """Display programs.yaml contents in a formatted way."""
    # Programs can be under 'programs' key or at root level
    programs = data.get('programs', {})
    if not programs:
        # Check if programs are at root level
        programs = {k: v for k, v in data.items()
                   if isinstance(v, dict) and (k.startswith('phenix.') or '.' in k)}

    print_subheader("Programs (%d total)" % len(programs))
    print()
    print(colored("  This file defines all PHENIX programs the agent can use.", Colors.DIM))
    print(colored("  Each program specifies its required/optional inputs, outputs,", Colors.DIM))
    print(colored("  and how to build commands. Edit this file to add new programs", Colors.DIM))
    print(colored("  or modify existing program definitions.", Colors.DIM))
    print()
    print(colored("  Format: program_name", Colors.CYAN))
    print(colored("          Description of what the program does", Colors.DIM))
    print(colored("          Required: input files that must be provided", Colors.GREEN))
    print(colored("          Optional: input files that may be provided", Colors.DIM))
    print(colored("          Outputs: files/metrics produced by the program", Colors.DIM))

    # Group by category
    by_category = {}
    for name, defn in programs.items():
        category = defn.get('category', 'other')
        if category not in by_category:
            by_category[category] = []
        by_category[category].append((name, defn))

    for category in sorted(by_category.keys()):
        print()
        print(colored("  [%s]" % category.upper(), Colors.YELLOW + Colors.BOLD))

        for name, defn in sorted(by_category[category]):
            desc = defn.get('description', 'No description')[:60]
            exp_types = defn.get('experiment_types', defn.get('experiment_type', ['any']))
            if isinstance(exp_types, str):
                exp_types = [exp_types]

            print("    %s" % colored(name, Colors.CYAN))
            print("      %s" % colored(desc, Colors.DIM))

            # Show inputs - handle nested structure
            inputs = defn.get('inputs', {})
            if inputs:
                required_inputs = []
                optional_inputs = []

                # Check for nested required/optional structure
                if 'required' in inputs and isinstance(inputs['required'], dict):
                    required_inputs = list(inputs['required'].keys())
                if 'optional' in inputs and isinstance(inputs['optional'], dict):
                    optional_inputs = list(inputs['optional'].keys())

                # Also check for flat structure (input_name: {required: true})
                if not required_inputs and not optional_inputs:
                    for input_name, input_def in inputs.items():
                        if isinstance(input_def, dict):
                            if input_def.get('required'):
                                required_inputs.append(input_name)
                            else:
                                optional_inputs.append(input_name)
                        elif input_name not in ('required', 'optional'):
                            optional_inputs.append(input_name)

                if required_inputs:
                    print("      Required: %s" % colored(", ".join(required_inputs), Colors.GREEN))
                if optional_inputs:
                    print("      Optional: %s" % colored(", ".join(optional_inputs), Colors.DIM))

            # Show outputs
            outputs = defn.get('outputs', {})
            if outputs:
                if isinstance(outputs, dict):
                    # Outputs can be {metrics: [...], files: [...]}
                    output_parts = []
                    for out_type, out_list in outputs.items():
                        if out_list:
                            if isinstance(out_list, list):
                                # Extract just the names/patterns
                                items = []
                                for item in out_list[:4]:
                                    if isinstance(item, dict):
                                        # Get pattern or type from dict
                                        items.append(item.get('pattern', item.get('type', '?')))
                                    else:
                                        items.append(str(item))
                                output_parts.append("%s: %s" % (out_type, ", ".join(items)))
                            else:
                                output_parts.append(str(out_type))
                    if output_parts:
                        print("      Outputs: %s" % "; ".join(output_parts))
                elif isinstance(outputs, list):
                    print("      Outputs: %s" % ", ".join(str(o) for o in outputs))

            # Show invariants (Step 1)
            invariants = defn.get('invariants', [])
            if invariants:
                inv_names = [inv.get('name', '?') for inv in invariants if isinstance(inv, dict)]
                print("      Invariants: %s" % colored(", ".join(inv_names), Colors.YELLOW))

            # Show input_priorities (Step 2)
            input_priorities = defn.get('input_priorities', {})
            if input_priorities:
                priority_parts = []
                for input_name, prio_def in input_priorities.items():
                    if isinstance(prio_def, dict):
                        cats = prio_def.get('categories', [])
                        excludes = prio_def.get('exclude_categories', [])
                        if cats:
                            part = "%s: %s" % (input_name, "→".join(cats[:3]))
                            if excludes:
                                part += " (excludes: %s)" % ", ".join(excludes[:3])
                            priority_parts.append(part)
                if priority_parts:
                    print("      Input priorities: %s" % colored("; ".join(priority_parts), Colors.BLUE))

            # Show user_advice_keywords (Step 3)
            keywords = defn.get('user_advice_keywords', [])
            if keywords:
                kw_display = ", ".join(keywords[:4])
                if len(keywords) > 4:
                    kw_display += " (+%d more)" % (len(keywords) - 4)
                print("      User advice keywords: %s" % colored(kw_display, Colors.CYAN))


def display_workflows(data):
    """Display workflows.yaml contents in a formatted way."""
    # Workflows can be under 'workflows' key or at root level
    workflows = data.get('workflows', {})
    if not workflows:
        workflows = {k: v for k, v in data.items()
                    if isinstance(v, dict) and 'phases' in v}

    print_subheader("Workflows (%d total)" % len(workflows))
    print()
    print(colored("  This file defines the workflow state machine for each experiment type.", Colors.DIM))
    print(colored("  Workflows control which programs are valid at each phase and how", Colors.DIM))
    print(colored("  the agent progresses through structure determination.", Colors.DIM))
    print()
    print(colored("  Phases:", Colors.DIM))
    print(colored("    - Each workflow has sequential phases (analyze → build → refine → validate)", Colors.DIM))
    print(colored("    - Each phase lists which programs are valid", Colors.DIM))
    print(colored("    - Transitions define when to move to the next phase", Colors.DIM))
    print()
    print(colored("  Format: phase_name: program1, program2, ... or STOP", Colors.DIM))

    for wf_name, wf_def in workflows.items():
        print()
        print(colored("  [%s]" % wf_name.upper(), Colors.YELLOW + Colors.BOLD))

        phases = wf_def.get('phases', {})
        for phase_name, phase_def in phases.items():
            desc = phase_def.get('description', '')
            programs = phase_def.get('programs', [])

            if phase_def.get('stop'):
                print("    %s: %s" % (
                    colored(phase_name, Colors.GREEN),
                    colored("STOP", Colors.RED + Colors.BOLD)
                ))
            else:
                prog_parts = []
                priority_progs = []
                preferred_prog = None

                for p in programs:
                    if isinstance(p, str):
                        prog_parts.append(p)
                    elif isinstance(p, dict):
                        prog_name = p.get('program', '?')

                        # Check for preferred marker
                        if p.get('preferred'):
                            preferred_prog = prog_name
                            prog_parts.append(prog_name + "*")
                        # Check for priority_when
                        elif p.get('priority_when'):
                            priority_progs.append((prog_name, p['priority_when']))
                            prog_parts.append(prog_name + "†")
                        else:
                            prog_parts.append(prog_name)

                print("    %s: %s" % (
                    colored(phase_name, Colors.CYAN),
                    ", ".join(prog_parts) or "(no programs)"
                ))

                # Show legend if needed
                if preferred_prog or priority_progs:
                    notes = []
                    if preferred_prog:
                        notes.append("* = default preferred")
                    if priority_progs:
                        for prog, cond in priority_progs:
                            notes.append("† %s when %s" % (prog, cond))
                    print("      %s" % colored("; ".join(notes), Colors.DIM))

            if desc:
                print("      %s" % colored(desc, Colors.DIM))


def display_metrics(data):
    """Display metrics.yaml contents in a formatted way."""
    # Metrics can be under 'metrics' key or at root level
    metrics = data.get('metrics', {})
    if not metrics:
        # Skip special keys
        skip_keys = {'improvement', 'quality_rules'}
        metrics = {k: v for k, v in data.items()
                  if isinstance(v, dict) and k not in skip_keys and
                  ('thresholds' in v or 'direction' in v or 'resolution_dependent' in v)}

    print_subheader("Metrics (%d total)" % len(metrics))
    print()
    print(colored("  This file defines quality metrics and their thresholds.", Colors.DIM))
    print(colored("  The agent uses these to assess model quality and decide when to stop.", Colors.DIM))
    print()
    print(colored("  Direction:", Colors.DIM))
    print(colored("    - minimize: lower is better (e.g., R-free, clashscore)", Colors.GREEN + Colors.DIM))
    print(colored("    - maximize: higher is better (e.g., map_cc, completeness)", Colors.BLUE + Colors.DIM))
    print()
    print(colored("  Thresholds:", Colors.DIM))
    print(colored("    - good: Target quality level", Colors.DIM))
    print(colored("    - acceptable: Minimum acceptable quality", Colors.DIM))
    print(colored("    - Some metrics have resolution-dependent thresholds", Colors.DIM))

    for name, defn in sorted(metrics.items()):
        direction = defn.get('direction', 'unknown')
        desc = defn.get('description', '')

        dir_color = Colors.GREEN if direction == 'minimize' else Colors.BLUE
        print("  %s  %s" % (
            colored(name, Colors.CYAN),
            colored("(%s)" % direction, dir_color)
        ))

        # Show thresholds
        thresholds = defn.get('thresholds', {})
        if thresholds:
            good = thresholds.get('good')
            acceptable = thresholds.get('acceptable')
            if good is not None or acceptable is not None:
                parts = []
                if good is not None:
                    parts.append("good: %s" % good)
                if acceptable is not None:
                    parts.append("acceptable: %s" % acceptable)
                print("    Thresholds: %s" % ", ".join(parts))

        # Show resolution-dependent thresholds
        res_thresholds = defn.get('resolution_dependent', [])
        if res_thresholds:
            print("    Resolution-dependent thresholds:")
            for rt in res_thresholds:
                range_str = "%s-%sÅ" % (rt.get('min', 0), rt.get('max', '∞'))
                print("      %s: good=%s, acceptable=%s" % (
                    range_str, rt.get('good'), rt.get('acceptable')
                ))


def display_file(filepath):
    """Display a YAML file's contents in formatted way."""
    # If file doesn't exist, try to find it in the knowledge directory
    if not os.path.exists(filepath):
        # Try in knowledge subdirectory
        knowledge_path = os.path.join(script_dir, "knowledge", filepath)
        if os.path.exists(knowledge_path):
            filepath = knowledge_path
        else:
            # Try to find by basename in knowledge directory
            yaml_files = find_yaml_files(script_dir)
            for yf in yaml_files:
                if os.path.basename(yf) == filepath or os.path.basename(yf) == os.path.basename(filepath):
                    filepath = yf
                    break

    data, error = load_yaml_file(filepath)

    if error:
        print_error("Cannot load %s: %s" % (filepath, error))
        return

    filename = os.path.basename(filepath)
    print_header(filename)

    if filename == 'programs.yaml':
        display_programs(data)
    elif filename == 'workflows.yaml':
        display_workflows(data)
    elif filename == 'metrics.yaml':
        display_metrics(data)
    elif filename == 'file_categories.yaml':
        display_file_categories(data)
    else:
        # Generic display
        print()
        print(yaml.dump(data, default_flow_style=False, indent=2))


def display_file_categories(data):
    """Display file_categories.yaml contents in a formatted way."""
    print_subheader("File Categories (%d total)" % len(data))
    print()
    print(colored("  This file defines how files are categorized based on name/extension.", Colors.DIM))
    print(colored("  Categories are used in input_priorities and workflow state detection.", Colors.DIM))
    print()
    print(colored("  Primary categories: Match by file extension", Colors.DIM))
    print(colored("  Subcategories: Match by filename patterns within a parent category", Colors.DIM))
    print()

    # Group into primary (have extensions) and subcategories
    primary = []
    by_parent = {}

    for cat_name, cat_def in data.items():
        parent = cat_def.get("subcategory_of")
        if parent:
            if parent not in by_parent:
                by_parent[parent] = []
            by_parent[parent].append((cat_name, cat_def))
        elif cat_def.get("extensions"):
            primary.append((cat_name, cat_def))

    for cat_name, cat_def in sorted(primary):
        desc = cat_def.get("description", "")
        exts = cat_def.get("extensions", [])

        print("  %s" % colored(cat_name, Colors.CYAN + Colors.BOLD))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        if exts:
            print("    Extensions: %s" % ", ".join(exts))

        # Show subcategories
        subcats = by_parent.get(cat_name, [])
        if subcats:
            print("    Subcategories:")
            for sub_name, sub_def in sorted(subcats):
                sub_desc = sub_def.get("description", "")[:50]
                patterns = sub_def.get("patterns", [])
                excludes = sub_def.get("excludes", [])

                pattern_str = ""
                if patterns:
                    pattern_str = " patterns: %s" % ", ".join(patterns[:3])
                    if len(patterns) > 3:
                        pattern_str += "..."
                if excludes:
                    pattern_str += " (excludes: %s)" % ", ".join(excludes[:2])

                print("      %s - %s%s" % (
                    colored(sub_name, Colors.CYAN),
                    colored(sub_desc, Colors.DIM),
                    colored(pattern_str, Colors.DIM) if pattern_str else ""
                ))
        print()


# =============================================================================
# COMPARISON FUNCTIONS
# =============================================================================

def compare_yaml_files(file1, file2):
    """
    Compare two YAML files and show differences.

    Returns:
        dict: {added: [], removed: [], changed: []}
    """
    data1, err1 = load_yaml_file(file1)
    data2, err2 = load_yaml_file(file2)

    if err1:
        return {'error': "Cannot load %s: %s" % (file1, err1)}
    if err2:
        return {'error': "Cannot load %s: %s" % (file2, err2)}

    return _compare_dicts(data1, data2, "")


def _compare_dicts(old, new, path):
    """Recursively compare two dictionaries."""
    changes = {'added': [], 'removed': [], 'changed': []}

    if not isinstance(old, dict) or not isinstance(new, dict):
        if old != new:
            changes['changed'].append((path, old, new))
        return changes

    all_keys = set(old.keys()) | set(new.keys())

    for key in all_keys:
        key_path = "%s.%s" % (path, key) if path else key

        if key not in old:
            changes['added'].append((key_path, new[key]))
        elif key not in new:
            changes['removed'].append((key_path, old[key]))
        elif isinstance(old[key], dict) and isinstance(new[key], dict):
            sub_changes = _compare_dicts(old[key], new[key], key_path)
            changes['added'].extend(sub_changes['added'])
            changes['removed'].extend(sub_changes['removed'])
            changes['changed'].extend(sub_changes['changed'])
        elif old[key] != new[key]:
            changes['changed'].append((key_path, old[key], new[key]))

    return changes


def display_comparison(file1, file2):
    """Display comparison between two files."""
    print_header("Comparing: %s vs %s" % (
        os.path.basename(file1), os.path.basename(file2)
    ))

    result = compare_yaml_files(file1, file2)

    if 'error' in result:
        print_error(result['error'])
        return False

    has_changes = False

    # Show additions
    if result['added']:
        has_changes = True
        print()
        print(colored("ADDED (%d):" % len(result['added']), Colors.GREEN + Colors.BOLD))
        for path, value in result['added']:
            print("  + %s" % colored(path, Colors.GREEN))
            if isinstance(value, dict):
                for k, v in list(value.items())[:3]:
                    print("      %s: %s" % (k, _truncate(str(v), 50)))
                if len(value) > 3:
                    print("      ... and %d more" % (len(value) - 3))
            else:
                print("      = %s" % _truncate(str(value), 60))

    # Show removals
    if result['removed']:
        has_changes = True
        print()
        print(colored("REMOVED (%d):" % len(result['removed']), Colors.RED + Colors.BOLD))
        for path, value in result['removed']:
            print("  - %s" % colored(path, Colors.RED))

    # Show changes
    if result['changed']:
        has_changes = True
        print()
        print(colored("CHANGED (%d):" % len(result['changed']), Colors.YELLOW + Colors.BOLD))
        for path, old_val, new_val in result['changed']:
            print("  ~ %s" % colored(path, Colors.YELLOW))
            print("      old: %s" % _truncate(str(old_val), 50))
            print("      new: %s" % _truncate(str(new_val), 50))

    if not has_changes:
        print_success("Files are identical")

    return has_changes


def compare_directories(dir1, dir2):
    """Compare all YAML files between two directories."""
    files1 = {os.path.basename(f): f for f in find_yaml_files(dir1)}
    files2 = {os.path.basename(f): f for f in find_yaml_files(dir2)}

    all_files = set(files1.keys()) | set(files2.keys())

    print_header("Comparing YAML configurations")
    print("  Old: %s" % dir1)
    print("  New: %s" % dir2)

    total_changes = 0

    for filename in sorted(all_files):
        if filename not in files1:
            print()
            print(colored("NEW FILE: %s" % filename, Colors.GREEN + Colors.BOLD))
            total_changes += 1
        elif filename not in files2:
            print()
            print(colored("DELETED FILE: %s" % filename, Colors.RED + Colors.BOLD))
            total_changes += 1
        else:
            if display_comparison(files1[filename], files2[filename]):
                total_changes += 1

    print()
    if total_changes > 0:
        print_warning("Total files with changes: %d" % total_changes)
    else:
        print_success("No differences found")


def _truncate(text, max_len):
    """Truncate text to max length."""
    if len(text) <= max_len:
        return text
    return text[:max_len - 3] + "..."


# =============================================================================
# VALIDATION COMMAND
# =============================================================================

def validate_files(paths):
    """Validate YAML files and report issues."""
    if not paths:
        # Default to knowledge directory
        paths = find_yaml_files(script_dir)

    print_header("YAML Validation Report")

    all_valid = True

    for filepath in paths:
        if not os.path.exists(filepath):
            print_error("File not found: %s" % filepath)
            all_valid = False
            continue

        filename = os.path.basename(filepath)
        print()
        print(colored("Validating: %s" % filename, Colors.CYAN + Colors.BOLD))

        # Load and check syntax
        data, error = load_yaml_file(filepath)

        if error:
            print_error("YAML Syntax Error:")
            print("  %s" % error)
            all_valid = False
            continue

        print_success("Syntax OK")

        # Validate structure
        issues = validate_yaml_structure(data, filepath)

        errors = [msg for level, msg in issues if level == 'error']
        warnings = [msg for level, msg in issues if level == 'warning']
        infos = [msg for level, msg in issues if level == 'info']

        if errors:
            all_valid = False
            print_error("Errors (%d):" % len(errors))
            for msg in errors:
                print("  - %s" % msg)

        if warnings:
            print_warning("Warnings (%d):" % len(warnings))
            for msg in warnings:
                print("  - %s" % msg)

        if infos:
            print_info("Info (%d):" % len(infos))
            for msg in infos:
                print("  - %s" % msg)

        if not errors and not warnings:
            print_success("Structure OK")

    print()
    if all_valid:
        print_success("All files validated successfully!")
    else:
        print_error("Validation found issues")

    return all_valid


# =============================================================================
# SUMMARY COMMAND
# =============================================================================

def show_summary():
    """Show summary of all YAML configuration."""
    yaml_files = find_yaml_files(script_dir)

    print_header("PHENIX AI Agent - YAML Configuration Summary")
    print()
    print(colored("  The AI Agent uses four YAML configuration files:", Colors.DIM))
    print()
    print(colored("  1. programs.yaml       - Defines available programs (inputs, outputs, commands)", Colors.CYAN))
    print(colored("  2. workflows.yaml      - Defines workflow phases and transitions", Colors.CYAN))
    print(colored("  3. metrics.yaml        - Defines quality metrics and thresholds", Colors.CYAN))
    print(colored("  4. file_categories.yaml - Defines file categorization rules", Colors.CYAN))
    print()
    print(colored("  Edit these files to customize agent behavior without changing code.", Colors.DIM))
    print(colored("  Use 'yaml_tools.py validate' to check for errors after editing.", Colors.DIM))
    print(colored("  Use 'yaml_tools.py terms' to see all defined terms and their usage.", Colors.DIM))
    print()
    print(colored("  Key Configuration Features:", Colors.YELLOW + Colors.BOLD))
    print(colored("    • invariants        - Auto-fix rules (e.g., auto-fill resolution)", Colors.DIM))
    print(colored("    • input_priorities  - File category preferences per program", Colors.DIM))
    print(colored("    • user_advice_keywords - Keywords to match user requests", Colors.DIM))
    print(colored("    • priority_when     - Conditional program priority in workflows", Colors.DIM))
    print(colored("    • file_categories   - How files are categorized by name/extension", Colors.DIM))

    for filepath in yaml_files:
        data, error = load_yaml_file(filepath)

        if error:
            print_error("%s: %s" % (os.path.basename(filepath), error))
            continue

        display_file(filepath)

    # Show feature statistics
    print()
    print_header("Feature Statistics")
    _show_feature_stats(yaml_files)

    print()
    print_header("Validation")
    validate_files(yaml_files)


def _show_feature_stats(yaml_files):
    """Show statistics about the new YAML features."""
    stats = {
        'programs_with_invariants': [],
        'programs_with_input_priorities': [],
        'programs_with_keywords': [],
        'phases_with_priority_when': [],
    }

    for filepath in yaml_files:
        data, error = load_yaml_file(filepath)
        if error:
            continue

        filename = os.path.basename(filepath).lower()

        if 'program' in filename:
            # Programs at root level
            programs = {k: v for k, v in data.items()
                       if isinstance(v, dict) and (k.startswith('phenix.') or '.' in k)}

            for prog_name, prog_def in programs.items():
                if prog_def.get('invariants'):
                    stats['programs_with_invariants'].append(prog_name)
                if prog_def.get('input_priorities'):
                    stats['programs_with_input_priorities'].append(prog_name)
                if prog_def.get('user_advice_keywords'):
                    stats['programs_with_keywords'].append(prog_name)

        elif 'workflow' in filename:
            workflows = {k: v for k, v in data.items()
                        if isinstance(v, dict) and 'phases' in v}

            for wf_name, wf_def in workflows.items():
                phases = wf_def.get('phases', {})
                for phase_name, phase_def in phases.items():
                    programs = phase_def.get('programs', [])
                    for p in programs:
                        if isinstance(p, dict) and p.get('priority_when'):
                            stats['phases_with_priority_when'].append(
                                "%s.%s.%s" % (wf_name, phase_name, p.get('program', '?'))
                            )

    print()
    print(colored("  Programs with invariants: %d" % len(stats['programs_with_invariants']), Colors.CYAN))
    if stats['programs_with_invariants']:
        print(colored("    " + ", ".join(stats['programs_with_invariants']), Colors.DIM))

    print(colored("  Programs with input_priorities: %d" % len(stats['programs_with_input_priorities']), Colors.CYAN))
    if stats['programs_with_input_priorities']:
        print(colored("    " + ", ".join(stats['programs_with_input_priorities']), Colors.DIM))

    print(colored("  Programs with user_advice_keywords: %d" % len(stats['programs_with_keywords']), Colors.CYAN))
    if stats['programs_with_keywords']:
        print(colored("    " + ", ".join(stats['programs_with_keywords']), Colors.DIM))

    print(colored("  Programs with priority_when: %d" % len(stats['phases_with_priority_when']), Colors.CYAN))
    if stats['phases_with_priority_when']:
        print(colored("    " + ", ".join(stats['phases_with_priority_when']), Colors.DIM))


# =============================================================================
# MAIN
# =============================================================================

# =============================================================================
# TERMS COMMAND - Show all defined terms and their usage
# =============================================================================

def collect_all_terms(yaml_files):
    """
    Collect all terms defined in the YAML system.

    Returns dict with:
        file_categories: {name: {definition, used_by, description}}
        strategy_keys: {name: {programs, types, description}}
        metrics: {name: {definition, description}}
        workflow_phases: {workflow: {phase: definition}}
        conditions: {name: {used_in, description}}
        check_types: {name: {description, example}}
        fix_types: {name: {description, example}}
    """
    terms = {
        "file_categories": {},
        "strategy_keys": {},
        "metrics": {},
        "workflow_phases": {},
        "conditions": {},
        "priority_conditions": {},
        "check_types": {},
        "fix_types": {},
        "user_advice_keywords": {},
    }

    programs_data = None
    workflows_data = None
    metrics_data = None
    file_categories_data = None

    # Load all YAML files
    for filepath in yaml_files:
        data, error = load_yaml_file(filepath)
        if error:
            continue

        filename = os.path.basename(filepath).lower()

        if 'program' in filename:
            programs_data = data
        elif 'workflow' in filename:
            workflows_data = data
        elif 'metric' in filename:
            metrics_data = data
        elif 'file_categories' in filename or 'categories' in filename:
            file_categories_data = data

    # Extract file categories
    if file_categories_data:
        for cat_name, cat_def in file_categories_data.items():
            terms["file_categories"][cat_name] = {
                "definition": cat_def,
                "description": cat_def.get("description", ""),
                "extensions": cat_def.get("extensions", []),
                "patterns": cat_def.get("patterns", []),
                "subcategory_of": cat_def.get("subcategory_of"),
                "used_by": [],  # Will be filled from programs
            }

    # Extract from programs
    if programs_data:
        programs = {k: v for k, v in programs_data.items()
                   if isinstance(v, dict) and (k.startswith('phenix.') or '.' in k)}

        for prog_name, prog_def in programs.items():
            # Extract strategy keys
            for key, key_def in prog_def.get("strategy_flags", {}).items():
                if key not in terms["strategy_keys"]:
                    terms["strategy_keys"][key] = {
                        "programs": [],
                        "type": key_def.get("type", "unknown"),
                        "description": key_def.get("hint", ""),
                    }
                terms["strategy_keys"][key]["programs"].append(prog_name)

            # Extract input_priorities usage of file categories
            for input_name, prio_def in prog_def.get("input_priorities", {}).items():
                for cat in prio_def.get("categories", []):
                    if cat in terms["file_categories"]:
                        usage = "%s (%s: preferred)" % (prog_name, input_name)
                        if usage not in terms["file_categories"][cat]["used_by"]:
                            terms["file_categories"][cat]["used_by"].append(usage)
                for cat in prio_def.get("exclude_categories", []):
                    if cat in terms["file_categories"]:
                        usage = "%s (%s: excluded)" % (prog_name, input_name)
                        if usage not in terms["file_categories"][cat]["used_by"]:
                            terms["file_categories"][cat]["used_by"].append(usage)

            # Extract invariant check/fix types
            for inv in prog_def.get("invariants", []):
                check = inv.get("check", {})
                _extract_check_types(check, terms["check_types"], prog_name)

                fix = inv.get("fix", {})
                for fix_type in fix.keys():
                    if fix_type not in terms["fix_types"]:
                        terms["fix_types"][fix_type] = {
                            "programs": [],
                            "description": _get_fix_description(fix_type),
                        }
                    if prog_name not in terms["fix_types"][fix_type]["programs"]:
                        terms["fix_types"][fix_type]["programs"].append(prog_name)

            # Extract user_advice_keywords
            for kw in prog_def.get("user_advice_keywords", []):
                if kw not in terms["user_advice_keywords"]:
                    terms["user_advice_keywords"][kw] = {"programs": []}
                terms["user_advice_keywords"][kw]["programs"].append(prog_name)

    # Extract from workflows
    if workflows_data:
        workflows = {k: v for k, v in workflows_data.items()
                    if isinstance(v, dict) and 'phases' in v}

        for wf_name, wf_def in workflows.items():
            terms["workflow_phases"][wf_name] = {}

            for phase_name, phase_def in wf_def.get("phases", {}).items():
                terms["workflow_phases"][wf_name][phase_name] = {
                    "description": phase_def.get("description", ""),
                    "programs": [],
                }

                for prog_entry in phase_def.get("programs", []):
                    if isinstance(prog_entry, dict):
                        prog_name = prog_entry.get("program", "?")
                        terms["workflow_phases"][wf_name][phase_name]["programs"].append(prog_name)

                        # Extract conditions
                        for cond in prog_entry.get("conditions", []):
                            if isinstance(cond, str):
                                cond_name = cond
                            elif isinstance(cond, dict):
                                cond_name = str(cond)
                            else:
                                continue

                            if cond_name not in terms["conditions"]:
                                terms["conditions"][cond_name] = {"used_in": []}
                            terms["conditions"][cond_name]["used_in"].append(
                                "%s.%s.%s" % (wf_name, phase_name, prog_name)
                            )

                        # Extract priority_when
                        priority_when = prog_entry.get("priority_when")
                        if priority_when:
                            if priority_when not in terms["priority_conditions"]:
                                terms["priority_conditions"][priority_when] = {"used_by": []}
                            terms["priority_conditions"][priority_when]["used_by"].append(
                                "%s.%s.%s" % (wf_name, phase_name, prog_name)
                            )
                    elif isinstance(prog_entry, str):
                        terms["workflow_phases"][wf_name][phase_name]["programs"].append(prog_entry)

    # Extract metrics
    if metrics_data:
        for metric_name, metric_def in metrics_data.items():
            if isinstance(metric_def, dict) and ('thresholds' in metric_def or 'direction' in metric_def):
                terms["metrics"][metric_name] = {
                    "description": metric_def.get("description", ""),
                    "direction": metric_def.get("direction", "unknown"),
                    "thresholds": metric_def.get("thresholds", {}),
                }

    return terms


def _extract_check_types(check, check_types, prog_name):
    """Recursively extract check types from invariant check structure."""
    for check_type, value in check.items():
        if check_type not in check_types:
            check_types[check_type] = {
                "programs": [],
                "description": _get_check_description(check_type),
            }
        if prog_name not in check_types[check_type]["programs"]:
            check_types[check_type]["programs"].append(prog_name)

        # Recurse into any_of / all_of
        if check_type in ("any_of", "all_of") and isinstance(value, list):
            for sub_check in value:
                if isinstance(sub_check, dict):
                    _extract_check_types(sub_check, check_types, prog_name)


def _get_check_description(check_type):
    """Get description for a check type."""
    descriptions = {
        "has_file": "Check if any of the named input slots have files",
        "has_strategy": "Check if any of the strategy keys are set",
        "strategy_equals": "Check if a strategy key equals a specific value",
        "any_of": "True if any sub-check passes",
        "all_of": "True if all sub-checks pass",
    }
    return descriptions.get(check_type, "")


def _get_fix_description(fix_type):
    """Get description for a fix type."""
    descriptions = {
        "set_strategy": "Set a strategy key to a specific value",
        "auto_fill_resolution": "Auto-fill resolution from context",
        "set_file": "Set a file input",
        "remove_file": "Remove a file input",
    }
    return descriptions.get(fix_type, "")


def show_terms(yaml_files, detail_level="simple"):
    """
    Show all terms defined in the YAML system.

    Args:
        yaml_files: List of YAML file paths
        detail_level: "simple", "normal", or "full"
    """
    terms = collect_all_terms(yaml_files)

    print_header("PHENIX AI Agent - Term Reference")

    if detail_level == "simple":
        _show_terms_simple(terms)
    elif detail_level == "normal":
        _show_terms_normal(terms)
    else:  # full
        _show_terms_full(terms)


def _show_terms_simple(terms):
    """Show simple list of terms."""
    print()
    print(colored("  FILE CATEGORIES", Colors.YELLOW + Colors.BOLD))
    cats = sorted(terms["file_categories"].keys())
    print("    " + ", ".join(cats))

    print()
    print(colored("  STRATEGY KEYS", Colors.YELLOW + Colors.BOLD))
    keys = sorted(terms["strategy_keys"].keys())
    print("    " + ", ".join(keys))

    print()
    print(colored("  METRICS", Colors.YELLOW + Colors.BOLD))
    metrics = sorted(terms["metrics"].keys())
    print("    " + ", ".join(metrics))

    print()
    print(colored("  WORKFLOW PHASES", Colors.YELLOW + Colors.BOLD))
    for wf_name, phases in terms["workflow_phases"].items():
        phase_names = sorted(phases.keys())
        print("    %s: %s" % (colored(wf_name, Colors.CYAN), ", ".join(phase_names)))

    print()
    print(colored("  CONDITIONS", Colors.YELLOW + Colors.BOLD))
    conds = sorted(terms["conditions"].keys())
    print("    " + ", ".join(conds[:10]))
    if len(conds) > 10:
        print("    ... and %d more" % (len(conds) - 10))

    print()
    print(colored("  PRIORITY CONDITIONS", Colors.YELLOW + Colors.BOLD))
    prio_conds = sorted(terms["priority_conditions"].keys())
    print("    " + ", ".join(prio_conds) if prio_conds else "    (none)")

    print()
    print(colored("  INVARIANT CHECK TYPES", Colors.YELLOW + Colors.BOLD))
    checks = sorted(terms["check_types"].keys())
    print("    " + ", ".join(checks))

    print()
    print(colored("  INVARIANT FIX TYPES", Colors.YELLOW + Colors.BOLD))
    fixes = sorted(terms["fix_types"].keys())
    print("    " + ", ".join(fixes))

    print()
    print(colored("  USER ADVICE KEYWORDS", Colors.YELLOW + Colors.BOLD))
    kws = sorted(terms["user_advice_keywords"].keys())
    print("    " + ", ".join(kws[:8]))
    if len(kws) > 8:
        print("    ... and %d more" % (len(kws) - 8))


def _show_terms_normal(terms):
    """Show terms with descriptions."""
    # File Categories
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  FILE CATEGORIES", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  Categories used in input_priorities and file selection", Colors.DIM))
    print()

    # Group by parent category
    primary = []
    by_parent = {}
    for cat_name, cat_def in terms["file_categories"].items():
        parent = cat_def.get("subcategory_of")
        if parent:
            if parent not in by_parent:
                by_parent[parent] = []
            by_parent[parent].append((cat_name, cat_def))
        else:
            primary.append((cat_name, cat_def))

    for cat_name, cat_def in sorted(primary):
        desc = cat_def.get("description", "")
        exts = cat_def.get("extensions", [])
        print("  %s" % colored(cat_name, Colors.CYAN + Colors.BOLD))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        if exts:
            print("    Extensions: %s" % ", ".join(exts))

        # Show subcategories
        for sub_name, sub_def in sorted(by_parent.get(cat_name, [])):
            sub_desc = sub_def.get("description", "")
            print("      %s - %s" % (colored(sub_name, Colors.CYAN), sub_desc))
        print()

    # Strategy Keys
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  STRATEGY KEYS", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  Parameters that can be passed to programs", Colors.DIM))
    print()

    for key_name, key_def in sorted(terms["strategy_keys"].items()):
        key_type = key_def.get("type", "unknown")
        desc = key_def.get("description", "")
        progs = key_def.get("programs", [])

        print("  %s (%s)" % (colored(key_name, Colors.CYAN + Colors.BOLD), key_type))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        print("    Programs: %s" % ", ".join(progs[:4]))
        if len(progs) > 4:
            print("              ... and %d more" % (len(progs) - 4))
        print()

    # Metrics
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  METRICS", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  Quality metrics tracked by the agent", Colors.DIM))
    print()

    for metric_name, metric_def in sorted(terms["metrics"].items()):
        direction = metric_def.get("direction", "unknown")
        desc = metric_def.get("description", "")
        thresholds = metric_def.get("thresholds", {})

        dir_color = Colors.GREEN if direction == "minimize" else Colors.BLUE
        print("  %s (%s)" % (
            colored(metric_name, Colors.CYAN + Colors.BOLD),
            colored(direction, dir_color)
        ))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        if thresholds:
            thresh_str = ", ".join("%s: %s" % (k, v) for k, v in sorted(thresholds.items()))
            print("    Thresholds: %s" % thresh_str)
        print()

    # Invariant Types
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  INVARIANT CHECK TYPES", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for check_type, check_def in sorted(terms["check_types"].items()):
        desc = check_def.get("description", "")
        progs = check_def.get("programs", [])
        print("  %s" % colored(check_type, Colors.CYAN + Colors.BOLD))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        print("    Used by: %s" % ", ".join(progs))
        print()

    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  INVARIANT FIX TYPES", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for fix_type, fix_def in sorted(terms["fix_types"].items()):
        desc = fix_def.get("description", "")
        progs = fix_def.get("programs", [])
        print("  %s" % colored(fix_type, Colors.CYAN + Colors.BOLD))
        if desc:
            print("    %s" % colored(desc, Colors.DIM))
        print("    Used by: %s" % ", ".join(progs))
        print()


def _show_terms_full(terms):
    """Show full details including cross-references."""
    # Start with normal output
    _show_terms_normal(terms)

    # Add cross-reference sections
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  FILE CATEGORY USAGE (Cross-Reference)", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for cat_name, cat_def in sorted(terms["file_categories"].items()):
        used_by = cat_def.get("used_by", [])
        if used_by:
            print("  %s" % colored(cat_name, Colors.CYAN + Colors.BOLD))
            for usage in used_by:
                print("    - %s" % usage)
            print()

    # Workflow phases with programs
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  WORKFLOW PHASES (Full Detail)", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for wf_name, phases in terms["workflow_phases"].items():
        print(colored("  [%s]" % wf_name.upper(), Colors.YELLOW + Colors.BOLD))

        for phase_name, phase_def in phases.items():
            desc = phase_def.get("description", "")
            progs = phase_def.get("programs", [])

            print("    %s" % colored(phase_name, Colors.CYAN + Colors.BOLD))
            if desc:
                print("      %s" % colored(desc, Colors.DIM))
            if progs:
                print("      Programs: %s" % ", ".join(progs))
        print()

    # Conditions
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  CONDITIONS (Full Detail)", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for cond_name, cond_def in sorted(terms["conditions"].items()):
        used_in = cond_def.get("used_in", [])
        print("  %s" % colored(cond_name, Colors.CYAN + Colors.BOLD))
        for usage in used_in:
            print("    - %s" % usage)
        print()

    # Priority conditions
    if terms["priority_conditions"]:
        print()
        print(colored("=" * 60, Colors.BOLD))
        print(colored("  PRIORITY CONDITIONS", Colors.YELLOW + Colors.BOLD))
        print(colored("=" * 60, Colors.BOLD))
        print()

        for cond_name, cond_def in sorted(terms["priority_conditions"].items()):
            used_by = cond_def.get("used_by", [])
            print("  %s" % colored(cond_name, Colors.CYAN + Colors.BOLD))
            for usage in used_by:
                print("    - %s" % usage)
            print()

    # User advice keywords
    print()
    print(colored("=" * 60, Colors.BOLD))
    print(colored("  USER ADVICE KEYWORDS (Full Detail)", Colors.YELLOW + Colors.BOLD))
    print(colored("=" * 60, Colors.BOLD))
    print()

    for kw, kw_def in sorted(terms["user_advice_keywords"].items()):
        progs = kw_def.get("programs", [])
        print("  %s → %s" % (colored(kw, Colors.CYAN), ", ".join(progs)))


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="YAML Tools for PHENIX AI Agent Configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  list      List all available YAML configuration files
  validate  Check YAML files for syntax and structural errors
  display   Show formatted contents of a YAML file
  compare   Compare two YAML files or directories
  summary   Show overview of all configuration files
  terms     Show all defined terms in the configuration system

Examples:
  %(prog)s list                        # List all available YAML files
  %(prog)s validate                    # Validate all YAML files
  %(prog)s validate programs.yaml      # Validate specific file
  %(prog)s display programs.yaml       # Show formatted contents
  %(prog)s display file_categories.yaml
  %(prog)s compare old/ new/           # Compare two directories
  %(prog)s summary                     # Show summary of all files
  %(prog)s terms                       # Show all defined terms (simple)
  %(prog)s terms --detail normal       # Show terms with descriptions
  %(prog)s terms --detail full         # Show full cross-reference

Detail levels for 'terms' command:
  simple  Just list term names (default)
  normal  Include descriptions and basic info
  full    Full cross-reference showing where each term is used
        """
    )

    parser.add_argument('command', choices=['list', 'validate', 'display', 'compare', 'summary', 'terms'],
                       help='Command to run (see Commands above)')
    parser.add_argument('paths', nargs='*', help='File or directory paths')
    parser.add_argument('--no-color', action='store_true', help='Disable colored output')
    parser.add_argument('--detail', choices=['simple', 'normal', 'full'], default='simple',
                       help='Detail level for terms command: simple (names only), '
                            'normal (with descriptions), full (cross-reference)')

    args = parser.parse_args()

    # Disable colors if requested or not a terminal
    if args.no_color or not sys.stdout.isatty():
        Colors.disable()

    if args.command == 'list':
        list_yaml_files(args.paths)

    elif args.command == 'validate':
        success = validate_files(args.paths)
        sys.exit(0 if success else 1)

    elif args.command == 'display':
        if not args.paths:
            print_error("Please specify a file to display")
            sys.exit(1)
        for path in args.paths:
            display_file(path)

    elif args.command == 'compare':
        if len(args.paths) != 2:
            print_error("Compare requires exactly 2 paths")
            sys.exit(1)

        path1, path2 = args.paths

        if os.path.isdir(path1) and os.path.isdir(path2):
            compare_directories(path1, path2)
        elif os.path.isfile(path1) and os.path.isfile(path2):
            display_comparison(path1, path2)
        else:
            print_error("Both paths must be files or both must be directories")
            sys.exit(1)

    elif args.command == 'summary':
        show_summary()

    elif args.command == 'terms':
        yaml_files = find_yaml_files(script_dir)
        show_terms(yaml_files, args.detail)


if __name__ == "__main__":
    main()
