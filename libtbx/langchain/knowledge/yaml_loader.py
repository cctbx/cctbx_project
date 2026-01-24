"""
YAML Configuration Loader for PHENIX AI Agent.

This module loads the YAML configuration files that define:
- programs.yaml: All programs the agent can use
- workflows.yaml: State machines for each experiment type
- metrics.yaml: Quality metrics and thresholds

Usage:
    from libtbx.langchain.knowledge.yaml_loader import (
        load_programs,
        load_workflows,
        load_metrics,
        get_program,
        get_workflow,
        get_metric_threshold,
    )

    programs = load_programs()
    program_info = get_program("phenix.refine")
"""

from __future__ import absolute_import, division, print_function

import os
import re

# Import yaml - required dependency (pip install PyYAML)
try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required. Install with: pip install PyYAML")

# Cache loaded configurations
_PROGRAMS = None
_WORKFLOWS = None
_METRICS = None
_FILE_CATEGORIES = None


# =============================================================================
# FILE LOADING
# =============================================================================

def _get_knowledge_dir():
    """Get the path to the knowledge directory."""
    # Try relative to this file
    this_dir = os.path.dirname(os.path.abspath(__file__))
    knowledge_dir = os.path.join(this_dir, "knowledge")
    if os.path.isdir(knowledge_dir):
        return knowledge_dir

    # Try parent directory
    parent_dir = os.path.dirname(this_dir)
    knowledge_dir = os.path.join(parent_dir, "knowledge")
    if os.path.isdir(knowledge_dir):
        return knowledge_dir

    # Try same directory (if this file is in knowledge/)
    if os.path.exists(os.path.join(this_dir, "programs.yaml")):
        return this_dir

    raise FileNotFoundError("Cannot find knowledge directory with YAML files")


def _load_yaml_file(filename):
    """Load a YAML file from the knowledge directory."""
    knowledge_dir = _get_knowledge_dir()
    filepath = os.path.join(knowledge_dir, filename)

    if not os.path.exists(filepath):
        print(f"Warning: {filename} not found at {filepath}")
        return {}

    with open(filepath, 'r') as f:
        return yaml.safe_load(f) or {}


def load_programs(force_reload=False):
    """
    Load program definitions from programs.yaml.

    Returns:
        dict: Program name -> program definition
    """
    global _PROGRAMS
    if _PROGRAMS is None or force_reload:
        _PROGRAMS = _load_yaml_file("programs.yaml")
    return _PROGRAMS


def load_workflows(force_reload=False):
    """
    Load workflow definitions from workflows.yaml.

    Returns:
        dict: Workflow name -> workflow definition
    """
    global _WORKFLOWS
    if _WORKFLOWS is None or force_reload:
        _WORKFLOWS = _load_yaml_file("workflows.yaml")
    return _WORKFLOWS


def load_sanity_checks(force_reload=False):
    """
    Load sanity check definitions from workflows.yaml.

    Returns:
        dict: {"critical": [...], "warning": [...]}
    """
    workflows = load_workflows(force_reload)
    return workflows.get("sanity_checks", {"critical": [], "warning": []})


def load_metrics(force_reload=False):
    """
    Load metric definitions from metrics.yaml.

    Returns:
        dict: Metric name -> metric definition
    """
    global _METRICS
    if _METRICS is None or force_reload:
        _METRICS = _load_yaml_file("metrics.yaml")
    return _METRICS


def load_file_categories(force_reload=False):
    """
    Load file category definitions from file_categories.yaml.

    Returns:
        dict: Category name -> category definition
    """
    global _FILE_CATEGORIES
    if _FILE_CATEGORIES is None or force_reload:
        _FILE_CATEGORIES = _load_yaml_file("file_categories.yaml")
    return _FILE_CATEGORIES


def reload_all_configs():
    """
    Force reload all YAML configuration files.

    Call this at the start of a new agent session or when configs may have changed.
    This ensures the server picks up any changes to:
    - programs.yaml
    - workflows.yaml
    - metrics.yaml
    - file_categories.yaml
    """
    load_programs(force_reload=True)
    load_workflows(force_reload=True)
    load_metrics(force_reload=True)
    load_file_categories(force_reload=True)


def get_file_category(category_name):
    """
    Get definition for a specific file category.

    Args:
        category_name: Name of the category (e.g., "refined", "rsr_output")

    Returns:
        dict: Category definition or None if not found
    """
    categories = load_file_categories()
    return categories.get(category_name)


def get_all_file_categories():
    """
    Get all file category names.

    Returns:
        list: List of category names
    """
    return list(load_file_categories().keys())


# =============================================================================
# PROGRAM ACCESS
# =============================================================================

def get_program(program_name):
    """
    Get definition for a specific program.

    Args:
        program_name: Name of the program (e.g., "phenix.refine")

    Returns:
        dict: Program definition or None if not found
    """
    programs = load_programs()
    return programs.get(program_name)


def get_all_programs():
    """
    Get list of all known program names.

    Returns:
        list: Program names
    """
    return list(load_programs().keys())


def get_programs_for_experiment(experiment_type):
    """
    Get programs valid for a specific experiment type.

    Args:
        experiment_type: "xray" or "cryoem"

    Returns:
        list: Program names valid for this experiment type
    """
    programs = load_programs()
    result = []
    for name, defn in programs.items():
        exp_types = defn.get("experiment_types", [])
        if experiment_type in exp_types:
            result.append(name)
    return result


def get_programs_by_category(category):
    """
    Get programs in a specific category.

    Args:
        category: Category name (e.g., "refinement", "validation")

    Returns:
        list: Program names in this category
    """
    programs = load_programs()
    result = []
    for name, defn in programs.items():
        if defn.get("category") == category:
            result.append(name)
    return result


def get_program_command_template(program_name):
    """
    Get command template for a program.

    Args:
        program_name: Name of the program

    Returns:
        str: Command template or None
    """
    program = get_program(program_name)
    if program:
        return program.get("command")
    return None


def get_program_inputs(program_name):
    """
    Get input requirements for a program.

    Args:
        program_name: Name of the program

    Returns:
        dict: {"required": {...}, "optional": {...}}
    """
    program = get_program(program_name)
    if program:
        return program.get("inputs", {})
    return {}


def get_program_outputs(program_name):
    """
    Get expected outputs for a program.

    Args:
        program_name: Name of the program

    Returns:
        dict: {"files": [...], "metrics": [...]}
    """
    program = get_program(program_name)
    if program:
        return program.get("outputs", {})
    return {}


def get_all_output_file_patterns():
    """
    Get all output file patterns from all programs.

    This is used by utilities.scan_directory_for_output_files() to find
    output files produced by any PHENIX program.

    Returns:
        list: Unique glob patterns for output files (e.g., ["*_refine_*.pdb", "PHASER.*.pdb"])
    """
    programs = load_programs()
    patterns = set()

    for prog_name, prog_def in programs.items():
        outputs = prog_def.get("outputs", {})
        files = outputs.get("files", [])

        for file_def in files:
            if isinstance(file_def, dict):
                pattern = file_def.get("pattern")
                if pattern:
                    patterns.add(pattern)
            elif isinstance(file_def, str):
                patterns.add(file_def)

    return list(patterns)


def get_program_log_patterns(program_name):
    """
    Get log parsing patterns for a program.

    Args:
        program_name: Name of the program

    Returns:
        dict: metric_name -> pattern definition
    """
    program = get_program(program_name)
    if program:
        return program.get("log_parsing", {})
    return {}


def get_program_description(program_name):
    """
    Get human-readable description of a program.

    Args:
        program_name: Name of the program

    Returns:
        str: Description or empty string
    """
    program = get_program(program_name)
    if program:
        return program.get("description", "")
    return ""


def get_program_hints(program_name):
    """
    Get usage hints for a program.

    Args:
        program_name: Name of the program

    Returns:
        list: List of hint strings
    """
    program = get_program(program_name)
    if program:
        return program.get("hints", [])
    return []


# =============================================================================
# WORKFLOW ACCESS
# =============================================================================

def get_workflow(workflow_name):
    """
    Get definition for a specific workflow.

    Args:
        workflow_name: "xray" or "cryoem"

    Returns:
        dict: Workflow definition or None
    """
    workflows = load_workflows()
    return workflows.get(workflow_name)


def get_workflow_phases(workflow_name):
    """
    Get phases for a workflow.

    Args:
        workflow_name: "xray" or "cryoem"

    Returns:
        dict: Phase name -> phase definition
    """
    workflow = get_workflow(workflow_name)
    if workflow:
        return workflow.get("phases", {})
    return {}


def get_workflow_targets(workflow_name):
    """
    Get quality targets for a workflow.

    Args:
        workflow_name: "xray" or "cryoem"

    Returns:
        dict: Target name -> target definition
    """
    workflow = get_workflow(workflow_name)
    if workflow:
        return workflow.get("targets", {})
    return {}


def get_phase_programs(workflow_name, phase_name):
    """
    Get valid programs for a specific phase.

    Args:
        workflow_name: "xray" or "cryoem"
        phase_name: Phase name (e.g., "refine")

    Returns:
        list: Program definitions for this phase
    """
    phases = get_workflow_phases(workflow_name)
    phase = phases.get(phase_name, {})
    return phase.get("programs", [])


# =============================================================================
# METRICS ACCESS
# =============================================================================

def get_metric(metric_name):
    """
    Get definition for a specific metric.

    Args:
        metric_name: Name of the metric (e.g., "r_free")

    Returns:
        dict: Metric definition or None
    """
    metrics = load_metrics()
    return metrics.get(metric_name)


def get_metric_threshold(metric_name, level="acceptable", resolution=None):
    """
    Get threshold value for a metric.

    Args:
        metric_name: Name of the metric
        level: "good", "acceptable", or "poor"
        resolution: Optional resolution for resolution-dependent thresholds

    Returns:
        float: Threshold value or None
    """
    metric = get_metric(metric_name)
    if not metric:
        return None

    # Check for resolution-dependent thresholds
    if resolution and "by_resolution" in metric:
        for entry in metric["by_resolution"]:
            range_min, range_max = entry.get("range", [0, 999])
            if range_min <= resolution < range_max:
                return entry.get(level)

    # Fall back to default thresholds
    thresholds = metric.get("thresholds", {})
    return thresholds.get(level)


def get_metric_direction(metric_name):
    """
    Get whether higher or lower values are better.

    Args:
        metric_name: Name of the metric

    Returns:
        str: "minimize" or "maximize"
    """
    metric = get_metric(metric_name)
    if metric:
        return metric.get("direction", "minimize")
    return "minimize"


def get_metric_patterns(metric_name):
    """
    Get extraction patterns for a metric.

    Args:
        metric_name: Name of the metric

    Returns:
        list: Regex patterns for extraction
    """
    metric = get_metric(metric_name)
    if metric:
        extraction = metric.get("extraction", {})
        return extraction.get("patterns", [])
    return []


def is_metric_good(metric_name, value, resolution=None):
    """
    Check if a metric value is considered "good".

    Args:
        metric_name: Name of the metric
        value: The value to check
        resolution: Optional resolution for context

    Returns:
        bool: True if value is good
    """
    threshold = get_metric_threshold(metric_name, "good", resolution)
    if threshold is None or value is None:
        return False

    direction = get_metric_direction(metric_name)
    if direction == "minimize":
        return value <= threshold
    else:
        return value >= threshold


def is_metric_acceptable(metric_name, value, resolution=None):
    """
    Check if a metric value is considered "acceptable".

    Args:
        metric_name: Name of the metric
        value: The value to check
        resolution: Optional resolution for context

    Returns:
        bool: True if value is acceptable (includes good)
    """
    threshold = get_metric_threshold(metric_name, "acceptable", resolution)
    if threshold is None or value is None:
        return False

    direction = get_metric_direction(metric_name)
    if direction == "minimize":
        return value <= threshold
    else:
        return value >= threshold


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def extract_metric_from_log(log_text, metric_name, program_name=None):
    """
    Extract a metric value from log text.

    Uses patterns from:
    1. Program-specific patterns (if program_name provided)
    2. Global metric patterns

    Args:
        log_text: Log file content
        metric_name: Name of metric to extract
        program_name: Optional program name for program-specific patterns

    Returns:
        float or str or None: Extracted value
    """
    if not log_text:
        return None

    patterns = []
    extract_method = "first"  # or "last"
    value_type = "float"

    # Get program-specific patterns first
    if program_name:
        prog_patterns = get_program_log_patterns(program_name)
        if metric_name in prog_patterns:
            pattern_def = prog_patterns[metric_name]
            if isinstance(pattern_def, dict):
                patterns.append(pattern_def.get("pattern", ""))
                extract_method = pattern_def.get("extract", "first")
                value_type = pattern_def.get("type", "float")
            else:
                patterns.append(pattern_def)

    # Add global metric patterns
    patterns.extend(get_metric_patterns(metric_name))

    # Try each pattern
    all_matches = []
    for pattern in patterns:
        if not pattern:
            continue
        matches = re.findall(pattern, log_text, re.IGNORECASE)
        all_matches.extend(matches)

    if not all_matches:
        return None

    # Get the appropriate match
    if extract_method == "last":
        raw_value = all_matches[-1]
    else:
        raw_value = all_matches[0]

    # Convert to appropriate type
    if value_type == "float":
        try:
            return float(raw_value)
        except (ValueError, TypeError):
            return None
    elif value_type == "boolean":
        return True  # Pattern matched
    else:
        return str(raw_value).strip()


def get_target_r_free(resolution):
    """
    Get target R-free for a given resolution.

    Args:
        resolution: Resolution in Angstroms

    Returns:
        float: Target R-free value
    """
    return get_metric_threshold("r_free", "acceptable", resolution) or 0.25


def format_metric_value(metric_name, value):
    """
    Format a metric value for display.

    Args:
        metric_name: Name of the metric
        value: Value to format

    Returns:
        str: Formatted string
    """
    metric = get_metric(metric_name)
    if not metric or value is None:
        return str(value)

    fmt = metric.get("format", "{}")
    try:
        return fmt.format(value)
    except Exception:
        return str(value)


# =============================================================================
# VALIDATION
# =============================================================================

def validate_yaml_files():
    """
    Validate that all YAML files are present and well-formed.

    Returns:
        tuple: (is_valid, list of errors)
    """
    errors = []

    # Check programs.yaml
    try:
        programs = load_programs(force_reload=True)
        if not programs:
            errors.append("programs.yaml is empty or missing")
        else:
            # Check required fields
            for name, defn in programs.items():
                if not defn.get("description"):
                    errors.append(f"{name}: missing description")
                if not defn.get("command"):
                    errors.append(f"{name}: missing command template")
                if not defn.get("inputs"):
                    errors.append(f"{name}: missing inputs")
    except Exception as e:
        errors.append(f"Error loading programs.yaml: {e}")

    # Check workflows.yaml
    try:
        workflows = load_workflows(force_reload=True)
        if not workflows:
            errors.append("workflows.yaml is empty or missing")
        else:
            for name, defn in workflows.items():
                if name in ["xray", "cryoem"]:
                    if not defn.get("phases"):
                        errors.append(f"{name} workflow: missing phases")
    except Exception as e:
        errors.append(f"Error loading workflows.yaml: {e}")

    # Check metrics.yaml
    try:
        metrics = load_metrics(force_reload=True)
        if not metrics:
            errors.append("metrics.yaml is empty or missing")
    except Exception as e:
        errors.append(f"Error loading metrics.yaml: {e}")

    return (len(errors) == 0, errors)


# =============================================================================
# MAIN (for testing)
# =============================================================================

if __name__ == "__main__":
    print("Testing YAML loader...")
    print()

    # Validate files
    is_valid, errors = validate_yaml_files()
    if is_valid:
        print("✓ All YAML files valid")
    else:
        print("✗ Validation errors:")
        for e in errors:
            print(f"  - {e}")
    print()

    # Test program loading
    programs = load_programs()
    print(f"Loaded {len(programs)} programs:")
    for name in sorted(programs.keys())[:5]:
        print(f"  - {name}: {get_program_description(name)[:50]}...")
    print()

    # Test workflow loading
    workflows = load_workflows()
    print(f"Loaded {len(workflows)} workflows:")
    for name in workflows:
        if name in ["xray", "cryoem"]:
            phases = get_workflow_phases(name)
            print(f"  - {name}: {len(phases)} phases")
    print()

    # Test metrics loading
    metrics = load_metrics()
    print(f"Loaded {len(metrics)} metrics:")
    for name in ["r_free", "map_cc", "clashscore"]:
        if name in metrics:
            good = get_metric_threshold(name, "good")
            print(f"  - {name}: good threshold = {good}")
    print()

    # Test resolution-dependent thresholds
    print("R-free targets by resolution:")
    for res in [1.5, 2.0, 2.5, 3.0, 3.5]:
        target = get_target_r_free(res)
        print(f"  {res}Å: {target}")
