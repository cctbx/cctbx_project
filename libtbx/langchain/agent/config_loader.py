"""
Configuration loader for PHENIX AI Agent decision-making.

Loads decision_config.json and provides helper functions for accessing
thresholds, evaluating conditions, and getting recommendations.
"""

from __future__ import absolute_import, division, print_function
import json
import os

# Path to config file
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "decision_config.json")

# Cached config
_config = None


def load_config(path=None):
    """Load configuration from JSON file."""
    global _config
    if _config is None or path is not None:
        config_path = path or CONFIG_PATH
        with open(config_path) as f:
            _config = json.load(f)
    return _config


def get_config():
    """Get the loaded configuration (loads if needed)."""
    if _config is None:
        load_config()
    return _config


# =============================================================================
# THRESHOLD ACCESS
# =============================================================================

def get_threshold(experiment_type, key, resolution=None):
    """
    Get a threshold value, optionally resolution-dependent.
    
    Args:
        experiment_type: "xray" or "cryoem"
        key: Threshold key (e.g., "autobuild_rfree", "success_cc")
        resolution: Optional resolution for resolution-dependent thresholds
        
    Returns:
        Threshold value
    """
    config = get_config()
    thresholds = config["thresholds"].get(experiment_type, {})
    
    # Check for resolution-dependent threshold
    if resolution is not None and "resolution_dependent" in thresholds:
        res_dep = thresholds["resolution_dependent"]
        
        # Find matching resolution bin
        for bin_name, bin_config in res_dep.items():
            if bin_name == "default":
                continue
            
            bin_range = bin_config.get("range", {})
            min_res = bin_range.get("min", 0)
            max_res = bin_range.get("max", float("inf"))
            
            if min_res <= resolution < max_res:
                if key in bin_config:
                    return bin_config[key]
        
        # Fall back to default
        if "default" in res_dep and key in res_dep["default"]:
            return res_dep["default"][key]
    
    # Non-resolution-dependent threshold
    if key in thresholds:
        return thresholds[key]
    
    return None


def get_resolution_thresholds(resolution):
    """
    Get all resolution-dependent thresholds for X-ray.
    
    Args:
        resolution: Data resolution in Angstroms
        
    Returns:
        dict with autobuild_rfree, good_model_rfree, ligandfit_rfree
    """
    config = get_config()
    res_dep = config["thresholds"]["xray"].get("resolution_dependent", {})
    
    # Find matching resolution bin
    for bin_name, bin_config in res_dep.items():
        if bin_name == "default":
            continue
        
        bin_range = bin_config.get("range", {})
        min_res = bin_range.get("min", 0)
        max_res = bin_range.get("max", float("inf"))
        
        if min_res <= resolution < max_res:
            return {
                "autobuild_rfree": bin_config.get("autobuild_rfree"),
                "good_model_rfree": bin_config.get("good_model_rfree"),
                "ligandfit_rfree": bin_config.get("ligandfit_rfree"),
                "resolution_bin": bin_name,
            }
    
    # Fall back to default
    default = res_dep.get("default", {})
    return {
        "autobuild_rfree": default.get("autobuild_rfree", 0.35),
        "good_model_rfree": default.get("good_model_rfree", 0.25),
        "ligandfit_rfree": default.get("ligandfit_rfree", 0.35),
        "resolution_bin": "default",
    }


# =============================================================================
# CONDITION EVALUATION
# =============================================================================

def evaluate_condition(condition, context):
    """
    Evaluate a condition against a context dictionary.
    
    Args:
        condition: Can be:
            - A string like "refine_count == 0" or "resolution < 2.0"
            - A string with "and"/"or" like "has_sequence and has_phases"
            - A dict like {"operator": ">=", "value": 2}
            - A dict of conditions (all must be true)
        context: Dict with current values (r_free, resolution, refine_count, etc.)
        
    Returns:
        bool: Whether condition is satisfied
    """
    if condition is None:
        return True
    
    if isinstance(condition, bool):
        return condition
    
    if isinstance(condition, str):
        # Handle "and" conditions
        if " and " in condition:
            parts = condition.split(" and ")
            return all(evaluate_condition(p.strip(), context) for p in parts)
        
        # Handle "or" conditions
        if " or " in condition:
            parts = condition.split(" or ")
            return any(evaluate_condition(p.strip(), context) for p in parts)
        
        # Handle "not" prefix
        if condition.startswith("not "):
            return not evaluate_condition(condition[4:].strip(), context)
        
        # Simple string condition: "key operator value"
        try:
            # Handle comparison operators
            for op in ["==", "!=", ">=", "<=", ">", "<"]:
                if op in condition:
                    parts = condition.split(op)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        value_str = parts[1].strip()
                        
                        # Get actual value from context
                        actual = context.get(key)
                        if actual is None:
                            return False
                        
                        # Parse expected value
                        try:
                            if value_str.lower() == "true":
                                expected = True
                            elif value_str.lower() == "false":
                                expected = False
                            elif value_str.lower() == "null" or value_str.lower() == "none":
                                expected = None
                            elif "." in value_str:
                                expected = float(value_str)
                            else:
                                # Try int, then fall back to looking up in context
                                try:
                                    expected = int(value_str)
                                except ValueError:
                                    # It's a variable reference
                                    expected = context.get(value_str)
                                    if expected is None:
                                        return False
                        except ValueError:
                            expected = value_str
                        
                        # Evaluate
                        if op == "==":
                            return actual == expected
                        elif op == "!=":
                            return actual != expected
                        elif op == ">=":
                            return actual >= expected
                        elif op == "<=":
                            return actual <= expected
                        elif op == ">":
                            return actual > expected
                        elif op == "<":
                            return actual < expected
            
            # Boolean key check (e.g., "has_sequence" or "always")
            if condition.strip() == "always":
                return True
            return bool(context.get(condition.strip()))
            
        except Exception:
            return False
    
    elif isinstance(condition, dict):
        # Check if it's an operator dict or a dict of conditions
        if "operator" in condition:
            # Single condition: {"operator": ">=", "value": 2}
            # This form requires the key to be specified separately
            return True  # Can't evaluate without key
        else:
            # Dict of conditions: all must be true
            for key, cond_spec in condition.items():
                actual = context.get(key)
                
                if isinstance(cond_spec, dict) and "operator" in cond_spec:
                    op = cond_spec["operator"]
                    expected = cond_spec["value"]
                    
                    if actual is None and op != "==" and expected is not None:
                        return False
                    
                    if op == "==" and not (actual == expected):
                        return False
                    elif op == "!=" and not (actual != expected):
                        return False
                    elif op == ">=" and not (actual >= expected):
                        return False
                    elif op == "<=" and not (actual <= expected):
                        return False
                    elif op == ">" and not (actual > expected):
                        return False
                    elif op == "<" and not (actual < expected):
                        return False
                else:
                    # Simple equality check
                    if actual != cond_spec:
                        return False
            
            return True
    
    return False


def evaluate_conditions(conditions, context):
    """
    Evaluate multiple conditions (all must be true).
    
    Args:
        conditions: Dict of {key: condition_spec}
        context: Current state context
        
    Returns:
        bool: Whether all conditions are satisfied
    """
    if not conditions:
        return True
    
    return evaluate_condition(conditions, context)


# =============================================================================
# TIER 2 DEFAULTS
# =============================================================================

def get_tier2_defaults(program, context):
    """
    Get Tier 2 strong defaults that apply to a program given current context.
    
    Args:
        program: Program name (e.g., "phenix.refine")
        context: Dict with current state (r_free, resolution, refine_count, etc.)
        
    Returns:
        dict: {
            "key": {
                "value": default_value,
                "reason": "Explanation string"
            },
            ...
        }
    """
    config = get_config()
    tier2 = config["tiers"]["tier2_strong_defaults"]
    
    defaults = {}
    
    for rule_name, rule in tier2.items():
        # Skip non-rule entries like "description"
        if not isinstance(rule, dict) or "applies_to" not in rule:
            continue
        
        # Check if rule applies to this program
        applies_to = rule.get("applies_to")
        if applies_to and applies_to != program:
            continue
        
        # Check conditions
        condition = rule.get("condition")
        conditions = rule.get("conditions")
        
        condition_met = False
        if condition:
            condition_met = evaluate_condition(condition, context)
        elif conditions:
            condition_met = evaluate_conditions(conditions, context)
        else:
            condition_met = True
        
        if not condition_met:
            continue
        
        # Get the key and value
        key = rule.get("key", rule_name)
        
        if "default_value_from" in rule:
            # Value comes from context
            value = context.get(rule["default_value_from"])
        else:
            value = rule.get("default_value")
        
        if value is None:
            continue
        
        # Format the reason - handle None values
        reason_template = rule.get("reason_template", rule_name)
        try:
            safe_context = {k: (v if v is not None else "N/A") for k, v in context.items()}
            reason = reason_template.format(**safe_context)
        except (KeyError, ValueError, TypeError):
            reason = reason_template
        
        defaults[key] = {
            "value": value,
            "reason": reason,
            "override_warning": rule.get("override_warning", "Overriding default for %s" % key),
        }
    
    return defaults


# =============================================================================
# TIER 3 SUGGESTIONS
# =============================================================================

def get_tier3_suggestions(context):
    """
    Get Tier 3 soft suggestions that apply given current context.
    
    Args:
        context: Dict with current state
        
    Returns:
        list of suggestion strings
    """
    config = get_config()
    tier3 = config["tiers"]["tier3_soft_guidance"]
    
    suggestions = []
    
    for item in tier3.get("suggestions", []):
        condition = item.get("condition")
        
        if evaluate_condition(condition, context):
            suggestion_template = item.get("suggestion", "")
            try:
                safe_context = {k: (v if v is not None else "N/A") for k, v in context.items()}
                suggestion = suggestion_template.format(**safe_context)
            except (KeyError, ValueError, TypeError):
                suggestion = suggestion_template
            suggestions.append(suggestion)
    
    return suggestions


# =============================================================================
# PROGRAM RANKINGS
# =============================================================================

def get_program_rankings(state_name, context):
    """
    Get ranked program recommendations for a workflow state.
    
    Args:
        state_name: Workflow state name (e.g., "xray_refined")
        context: Dict with current state for condition evaluation
        
    Returns:
        list of {
            "program": str,
            "rank": int,
            "reason": str,
            "condition_met": bool
        }
    """
    config = get_config()
    rankings_config = config.get("program_rankings", {})
    
    state_config = rankings_config.get(state_name, {})
    rankings = state_config.get("rankings", [])
    
    result = []
    rank = 1
    
    for item in rankings:
        program = item.get("program")
        condition = item.get("condition")
        reason_template = item.get("reason", program)
        
        # Check condition
        if condition:
            condition_met = evaluate_condition(condition, context)
        else:
            condition_met = True
        
        # Format reason - handle None values
        try:
            # Replace None values with placeholder strings for formatting
            safe_context = {k: (v if v is not None else "N/A") for k, v in context.items()}
            reason = reason_template.format(**safe_context)
        except (KeyError, ValueError, TypeError):
            reason = reason_template
        
        if condition_met:
            result.append({
                "program": program,
                "rank": rank,
                "reason": reason,
                "condition_met": True,
            })
            rank += 1
        else:
            # Include but mark as not met (for transparency)
            result.append({
                "program": program,
                "rank": None,
                "reason": reason,
                "condition_met": False,
            })
    
    return result


def get_recommended_program(state_name, context):
    """
    Get the top recommended program for a workflow state.
    
    Args:
        state_name: Workflow state name
        context: Current state context
        
    Returns:
        tuple: (program_name, reason) or (None, None)
    """
    rankings = get_program_rankings(state_name, context)
    
    for item in rankings:
        if item["condition_met"]:
            return item["program"], item["reason"]
    
    return None, None


# =============================================================================
# VALIDATION CRITERIA
# =============================================================================

def get_validation_criteria():
    """Get validation criteria for model quality assessment."""
    config = get_config()
    return config.get("validation_criteria", {})


def assess_validation_quality(metrics):
    """
    Assess model quality based on validation metrics.
    
    Args:
        metrics: Dict with validation metrics (ramachandran_outliers, clashscore, etc.)
        
    Returns:
        dict: {
            "overall": "good" | "acceptable" | "poor",
            "issues": ["list of specific issues"],
            "details": {metric: {"value": v, "threshold": t, "status": s}}
        }
    """
    config = get_config()
    criteria = config.get("validation_criteria", {}).get("molprobity", {})
    good = criteria.get("good", {})
    poor_mult = criteria.get("poor_multiplier", 2.0)
    
    issues = []
    details = {}
    has_poor = False
    has_acceptable = False
    
    metric_mapping = {
        "ramachandran_outliers": "ramachandran_outliers_max",
        "rotamer_outliers": "rotamer_outliers_max",
        "clashscore": "clashscore_max",
        "rms_bonds": "rms_bonds_max",
        "rms_angles": "rms_angles_max",
        "molprobity_score": "molprobity_score_max",
    }
    
    for metric, threshold_key in metric_mapping.items():
        if metric not in metrics:
            continue
        
        value = metrics[metric]
        threshold = good.get(threshold_key)
        
        if threshold is None:
            continue
        
        poor_threshold = threshold * poor_mult
        
        if value > poor_threshold:
            status = "poor"
            has_poor = True
            issues.append("%s %.2f > %.2f (poor)" % (metric, value, poor_threshold))
        elif value > threshold:
            status = "acceptable"
            has_acceptable = True
            issues.append("%s %.2f > %.2f (acceptable)" % (metric, value, threshold))
        else:
            status = "good"
        
        details[metric] = {
            "value": value,
            "good_threshold": threshold,
            "poor_threshold": poor_threshold,
            "status": status,
        }
    
    if has_poor:
        overall = "poor"
    elif has_acceptable:
        overall = "acceptable"
    else:
        overall = "good"
    
    return {
        "overall": overall,
        "issues": issues,
        "details": details,
    }


# =============================================================================
# HARD CONSTRAINTS (Tier 1)
# =============================================================================

def get_hard_constraints(program=None):
    """
    Get Tier 1 hard constraints.
    
    Args:
        program: Optional program to filter constraints
        
    Returns:
        list of constraint descriptions
    """
    config = get_config()
    tier1 = config["tiers"]["tier1_hard_constraints"]
    
    constraints = []
    for rule_name, rule in tier1.get("rules", {}).items():
        applies_to = rule.get("applies_to", "all")
        
        if program is None or applies_to == "all" or applies_to == program:
            constraints.append(rule.get("description", rule_name))
    
    return constraints
