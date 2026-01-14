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
    """Find all YAML files in a directory."""
    yaml_files = []
    knowledge_dir = os.path.join(directory, "knowledge")

    if os.path.isdir(knowledge_dir):
        directory = knowledge_dir

    if os.path.isdir(directory):
        for f in os.listdir(directory):
            if f.endswith(('.yaml', '.yml')):
                yaml_files.append(os.path.join(directory, f))

    return sorted(yaml_files)


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
    else:
        # Try to auto-detect based on content
        if any(k.startswith('phenix.') for k in data.keys()):
            issues.extend(_validate_programs(data))
        elif any('phases' in v for v in data.values() if isinstance(v, dict)):
            issues.extend(_validate_workflows(data))
        elif any('direction' in v or 'thresholds' in v for v in data.values() if isinstance(v, dict)):
            issues.extend(_validate_metrics(data))

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
        'hints'
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
        'pattern', 'type', 'group', 'default', 'fallback_pattern', 'extract'
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

        # Check command template
        command = prog_def.get('command', {})
        if command and isinstance(command, dict):
            if 'base' not in command:
                issues.append(('warning', "Program '%s' command has no base" % prog_name))

    return issues


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
        'required', 'optional', 'requires_resolution', 'strategy_depends_on'
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
        'unit', 'format', 'display_name', 'category'
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
                prog_names = []
                for p in programs:
                    if isinstance(p, str):
                        prog_names.append(p)
                    elif isinstance(p, dict):
                        prog_names.append(p.get('program', '?'))

                print("    %s: %s" % (
                    colored(phase_name, Colors.CYAN),
                    ", ".join(prog_names) or "(no programs)"
                ))

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
    else:
        # Generic display
        print()
        print(yaml.dump(data, default_flow_style=False, indent=2))


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
    print(colored("  The AI Agent uses three YAML configuration files:", Colors.DIM))
    print()
    print(colored("  1. programs.yaml  - Defines available programs (inputs, outputs, commands)", Colors.CYAN))
    print(colored("  2. workflows.yaml - Defines workflow phases and transitions", Colors.CYAN))
    print(colored("  3. metrics.yaml   - Defines quality metrics and thresholds", Colors.CYAN))
    print()
    print(colored("  Edit these files to customize agent behavior without changing code.", Colors.DIM))
    print(colored("  Use 'yaml_tools.py validate' to check for errors after editing.", Colors.DIM))

    for filepath in yaml_files:
        data, error = load_yaml_file(filepath)

        if error:
            print_error("%s: %s" % (os.path.basename(filepath), error))
            continue

        display_file(filepath)

    print()
    print_header("Validation")
    validate_files(yaml_files)


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="YAML Tools for PHENIX AI Agent Configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s validate                    # Validate all YAML files
  %(prog)s validate programs.yaml      # Validate specific file
  %(prog)s display programs.yaml       # Show formatted contents
  %(prog)s compare old/ new/           # Compare two directories
  %(prog)s summary                     # Show summary of all files
        """
    )

    parser.add_argument('command', choices=['validate', 'display', 'compare', 'summary'],
                       help='Command to run')
    parser.add_argument('paths', nargs='*', help='File or directory paths')
    parser.add_argument('--no-color', action='store_true', help='Disable colored output')

    args = parser.parse_args()

    # Disable colors if requested or not a terminal
    if args.no_color or not sys.stdout.isatty():
        Colors.disable()

    if args.command == 'validate':
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


if __name__ == "__main__":
    main()
