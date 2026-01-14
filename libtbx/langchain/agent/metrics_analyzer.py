"""
Metrics Analyzer for PHENIX AI Agent.

Extracts metrics from history and analyzes trends to detect:
- Refinement plateaus (R-free stopped improving)
- Success conditions (R-free below target)
- Excessive refinement (too many consecutive cycles)

Design: Metrics are derived from history on each call, not persisted.
This keeps the API clean and ensures metrics always match history.

Uses YAML-driven evaluation via MetricEvaluator.
"""

from __future__ import absolute_import, division, print_function
import re

# Import YAML-driven evaluator
from libtbx.langchain.agent.metric_evaluator import analyze_refinement_trend


# =============================================================================
# METRICS EXTRACTION
# =============================================================================

def derive_metrics_from_history(history):
    """
    Extract metrics from each cycle in history.

    This reconstructs metrics_history from the history list that's
    passed from the client on every call. No new API fields needed.

    Args:
        history: List of cycle records from client. Can be:
            - List of dicts with 'program', 'analysis', 'result' keys
            - List of strings (legacy format)
            - Mixed

    Returns:
        list: [{cycle, program, r_free, r_work, tfz, llg, resolution, map_cc}, ...]
    """
    if not history:
        return []

    metrics_history = []

    for i, entry in enumerate(history):
        cycle_num = i + 1

        if isinstance(entry, str):
            # Legacy string format - try to parse metrics from text
            metrics = _parse_metrics_from_string(entry)
            metrics["cycle"] = cycle_num
            metrics_history.append(metrics)

        elif isinstance(entry, dict):
            metrics = {
                "cycle": entry.get("cycle_number", cycle_num),
                "program": entry.get("program", "unknown"),
                "r_free": None,
                "r_work": None,
                "tfz": None,
                "llg": None,
                "resolution": None,
                "map_cc": None,
            }

            # Extract from analysis dict if present
            analysis = entry.get("analysis", {})
            if isinstance(analysis, dict):
                metrics["r_free"] = analysis.get("r_free")
                metrics["r_work"] = analysis.get("r_work")
                metrics["tfz"] = analysis.get("tfz")
                metrics["llg"] = analysis.get("llg")
                metrics["resolution"] = analysis.get("resolution")
                metrics["map_cc"] = analysis.get("map_cc")

            # Also check result/summary for metrics (fallback)
            result_text = str(entry.get("result", "") or entry.get("summary", ""))

            # Only fill in missing values from result text
            if not metrics["r_free"]:
                metrics["r_free"] = _extract_float(result_text, r"R-free[:\s=]+([0-9.]+)")
            # Try autobuild table format: "SOLUTION  CYCLE     R        RFREE"
            # followed by data row like "1         1      0.21        0.25"
            if not metrics["r_free"]:
                metrics["r_free"] = _extract_autobuild_rfree(result_text)
            if not metrics["r_work"]:
                metrics["r_work"] = _extract_float(result_text, r"R-work[:\s=]+([0-9.]+)")
            if not metrics["tfz"]:
                metrics["tfz"] = _extract_float(result_text, r"TFZ[=:\s]+([0-9.]+)")
            if not metrics["llg"]:
                metrics["llg"] = _extract_float(result_text, r"LLG[=:\s]+([0-9.]+)")
            if not metrics["map_cc"]:
                metrics["map_cc"] = _extract_float(result_text, r"(?:map.model.CC|CC_mask|CC)[:\s=]+([0-9.]+)")

            # Detect program from command if not in entry
            if metrics["program"] == "unknown":
                cmd = entry.get("command", "")
                metrics["program"] = _detect_program_from_command(cmd)

            metrics_history.append(metrics)

    return metrics_history


def _extract_float(text, pattern):
    """Extract float from text using regex pattern."""
    if not text:
        return None
    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except (ValueError, IndexError):
            pass
    return None


def _extract_autobuild_rfree(text):
    """
    Extract R-free from autobuild table format.

    Autobuild outputs a table like:
        SOLUTION  CYCLE     R        RFREE     BUILT   PLACED
        1         1      0.21        0.25      0       0

    We want to extract the RFREE value (0.25 in this example).
    """
    if not text:
        return None

    # Look for the autobuild summary table header
    # Then find the data row and extract the 4th column (RFREE)
    lines = text.split('\n')
    header_idx = None
    rfree_col = None

    for i, line in enumerate(lines):
        # Find header line with RFREE column
        if 'RFREE' in line.upper() and 'SOLUTION' in line.upper():
            header_idx = i
            # Find which column RFREE is in
            parts = line.upper().split()
            for j, part in enumerate(parts):
                if part == 'RFREE':
                    rfree_col = j
                    break
            break

    if header_idx is not None and rfree_col is not None:
        # Look at the next few lines for data rows
        for line in lines[header_idx + 1:header_idx + 5]:
            parts = line.split()
            if len(parts) > rfree_col:
                try:
                    # First column should be a number (SOLUTION)
                    int(parts[0])
                    # Extract RFREE value
                    rfree = float(parts[rfree_col])
                    if 0 < rfree < 1:  # Sanity check
                        return rfree
                except (ValueError, IndexError):
                    continue

    return None


def _detect_program_from_command(cmd):
    """Detect program name from command string."""
    if not cmd:
        return "unknown"
    cmd_lower = cmd.lower()

    programs = [
        "phenix.real_space_refine",  # Check before refine
        "phenix.refine",
        "phenix.phaser",
        "phenix.xtriage",
        "phenix.mtriage",
        "phenix.predict_and_build",
        "phenix.ligandfit",
        "phenix.pdbtools",
        "phenix.autobuild",
        "phenix.autosol",
        "phenix.dock_in_map",
        "phenix.process_predicted_model",
    ]

    for prog in programs:
        if prog in cmd_lower:
            return prog

    return "unknown"


def _parse_metrics_from_string(text):
    """Parse metrics from a string history entry (legacy format)."""
    r_free = _extract_float(text, r"R-free[:\s=]+([0-9.]+)")
    if not r_free:
        r_free = _extract_autobuild_rfree(text)
    return {
        "program": _detect_program_from_command(text),
        "r_free": r_free,
        "r_work": _extract_float(text, r"R-work[:\s=]+([0-9.]+)"),
        "tfz": _extract_float(text, r"TFZ[=:\s]+([0-9.]+)"),
        "llg": _extract_float(text, r"LLG[=:\s]+([0-9.]+)"),
        "resolution": _extract_float(text, r"resolution[:\s=]+([0-9.]+)"),
        "map_cc": _extract_float(text, r"(?:map.model.CC|CC_mask|CC)[:\s=]+([0-9.]+)"),
    }


# =============================================================================
# TREND ANALYSIS
# =============================================================================

def analyze_metrics_trend(metrics_history, resolution=None, experiment_type="xray",
                         use_yaml_evaluator=False):
    """
    Analyze metrics history to detect plateaus and recommend stopping.

    Args:
        metrics_history: List of metrics dicts from derive_metrics_from_history()
        resolution: Data resolution in Angstroms (for dynamic stop threshold)
        experiment_type: "xray" or "cryoem"
        use_yaml_evaluator: If True, use YAML-driven MetricEvaluator

    Returns:
        dict: {
            should_stop: bool,
            reason: str or None,
            trend_summary: str,
            r_free_trend: list of last N R-free values,
            map_cc_trend: list of last N CC values (for cryo-EM),
            improvement_rate: float (% improvement last cycle),
            consecutive_refines: int,
            recommendation: "continue" | "stop" | "consider_stopping",
        }
    """
    # Use YAML-driven evaluator
    if use_yaml_evaluator:
        try:
            return analyze_refinement_trend(metrics_history, experiment_type, resolution)
        except Exception as e:
            import sys
            print("Warning: YAML evaluator failed, using hardcoded: %s" % e, file=sys.stderr)

    result = {
        "should_stop": False,
        "reason": None,
        "trend_summary": "No refinement history",
        "r_free_trend": [],
        "map_cc_trend": [],
        "improvement_rate": 0.0,
        "consecutive_refines": 0,
        "consecutive_rsr": 0,  # real_space_refine count
        "recommendation": "continue",
    }

    if not metrics_history:
        return result

    # Detect experiment type from metrics if not provided
    if experiment_type == "xray":
        # Check if we have map_cc but no r_free - might be cryo-EM
        has_r_free = any(m.get("r_free") for m in metrics_history)
        has_map_cc = any(m.get("map_cc") for m in metrics_history)
        if has_map_cc and not has_r_free:
            experiment_type = "cryoem"

    if experiment_type == "cryoem":
        return _analyze_cryoem_trend(metrics_history, result)
    else:
        return _analyze_xray_trend(metrics_history, resolution, result)


def _analyze_xray_trend(metrics_history, resolution, result):
    """Analyze X-ray refinement trend (R-free based)."""

    # Count consecutive refines at end of history
    # Also count 'unknown' as potentially refine if it has r_free
    consecutive = 0
    for m in reversed(metrics_history):
        prog = m.get("program", "").lower()
        has_r_free = m.get("r_free") is not None
        if "refine" in prog and "real_space" not in prog:
            consecutive += 1
        elif prog == "unknown" and has_r_free:
            # Unknown with R-free is likely a refine cycle
            consecutive += 1
        else:
            break
    result["consecutive_refines"] = consecutive

    # Extract R-free values from refine cycles (or unknown with r_free in refine context)
    # First, identify if we're in a refine-dominated workflow
    refine_count = sum(1 for m in metrics_history if "refine" in m.get("program", "").lower())
    in_refine_context = refine_count > 0

    r_free_values = []
    for m in metrics_history:
        prog = m.get("program", "").lower()
        r_free = m.get("r_free")
        if r_free is None:
            continue

        # Include if explicitly refine
        if "refine" in prog and "real_space" not in prog:
            r_free_values.append(r_free)
        # Or if unknown but in a refine context (has r_free and we've seen refines before)
        elif prog == "unknown" and in_refine_context:
            r_free_values.append(r_free)

    result["r_free_trend"] = r_free_values[-5:]

    if len(r_free_values) < 1:
        result["trend_summary"] = "No refinement cycles with R-free yet"
        return result

    latest_r_free = r_free_values[-1]

    # Calculate dynamic stop threshold based on resolution
    # Rule: R-free ~ resolution/10, but bounded [0.20, 0.30]
    if resolution and resolution > 0:
        dynamic_target = min(0.30, max(0.20, resolution / 10.0))
    else:
        dynamic_target = 0.25  # Default

    # === SUCCESS CHECK (do this first, even with only 1 value) ===
    # 1. SUCCESS: R-free well below target (at least 0.02 below)
    # For borderline cases (within 0.02 of target), let LLM/workflow decide
    # This allows autobuild to be tried when R-free is close but not great
    success_threshold = dynamic_target - 0.02  # e.g., 0.28 for 3Å data

    if latest_r_free < success_threshold:
        # Check if validation has been done - don't auto-stop without validation
        validation_done = any(
            m.get("program") in ("phenix.molprobity", "phenix.model_vs_data")
            for m in metrics_history
        )

        if validation_done:
            result["should_stop"] = True
            result["reason"] = "SUCCESS: R-free (%.3f) well below target (%.2f)" % (latest_r_free, dynamic_target)
            result["recommendation"] = "stop"
            result["trend_summary"] = "R-free: %.3f - BELOW TARGET" % latest_r_free
            return result
        else:
            # Success reached but validation not done - recommend validation, not stop
            result["should_stop"] = False
            result["reason"] = "R-free (%.3f) below target - recommend validation before stopping" % latest_r_free
            result["recommendation"] = "validate"
            result["trend_summary"] = "R-free: %.3f - TARGET REACHED, VALIDATE BEFORE STOPPING" % latest_r_free
            result["suggest_validation"] = True
            return result

    # Borderline case - close to target but could be better
    is_borderline = latest_r_free < dynamic_target and latest_r_free >= success_threshold

    # Build trend summary string
    if len(r_free_values) >= 2:
        trend_parts = ["%0.3f" % r for r in r_free_values[-5:]]
        previous = r_free_values[-2]
        improvement = (previous - latest_r_free) / previous * 100 if previous > 0 else 0
        result["improvement_rate"] = improvement

        # Add trend direction
        if improvement > 1.0:
            direction = "improving"
        elif improvement > 0:
            direction = "slightly improving"
        elif improvement > -1.0:
            direction = "stable"
        else:
            direction = "WORSENING"

        result["trend_summary"] = "R-free: %s (%s, %+.1f%% last cycle)" % (
            " → ".join(trend_parts), direction, improvement
        )

        # Add autobuild suggestion for borderline cases
        if is_borderline:
            result["trend_summary"] += " - near target, could try autobuild"
    else:
        result["trend_summary"] = "R-free: %0.3f (first refinement)" % latest_r_free
        if is_borderline:
            result["trend_summary"] += " - near target, could try autobuild"
        return result

    # === STOP CONDITIONS ===

    # 1. SUCCESS already checked above

    # 2. PLATEAU: Less than 0.3% improvement for 3+ cycles (more patience before giving up)
    if len(r_free_values) >= 4:
        recent_improvements = []
        # Calculate improvements between consecutive cycles (last 3)
        for i in range(len(r_free_values) - 1, 0, -1):
            if i <= len(r_free_values) - 4:  # Only look at last 3 improvements
                break
            if r_free_values[i-1] > 0:
                imp = (r_free_values[i-1] - r_free_values[i]) / r_free_values[i-1] * 100
                recent_improvements.append(imp)

        # All recent improvements below threshold
        if len(recent_improvements) >= 3 and all(imp < 0.3 for imp in recent_improvements):
            result["should_stop"] = True
            result["reason"] = "PLATEAU: R-free improvement stalled (<0.3%% for %d cycles)" % len(recent_improvements)
            result["recommendation"] = "stop"
            return result

    # 3. EXCESSIVE REFINEMENT: 8+ consecutive refines (was 5, now more patience)
    if consecutive >= 8:
        result["should_stop"] = True
        result["reason"] = "EXCESSIVE: %d consecutive refinement cycles" % consecutive
        result["recommendation"] = "stop"
        return result

    # 4. WORSENING: R-free got significantly worse
    if result["improvement_rate"] < -2.0:
        result["recommendation"] = "consider_stopping"
        result["trend_summary"] += " ⚠️"

    # 5. DIMINISHING RETURNS: Warn but don't stop
    elif 0 < result["improvement_rate"] < 1.0 and consecutive >= 3:
        result["recommendation"] = "consider_stopping"

    return result


def _analyze_cryoem_trend(metrics_history, result):
    """Analyze cryo-EM refinement trend (map-model CC based)."""

    # Count consecutive real_space_refine at end of history
    consecutive = 0
    for m in reversed(metrics_history):
        prog = m.get("program", "").lower()
        if "real_space" in prog:
            consecutive += 1
        else:
            break
    result["consecutive_rsr"] = consecutive

    # Extract CC values from real_space_refine cycles
    rsr_metrics = [
        m for m in metrics_history
        if "real_space" in m.get("program", "").lower()
    ]
    cc_values = [m["map_cc"] for m in rsr_metrics if m.get("map_cc") is not None]
    result["map_cc_trend"] = cc_values[-5:]

    if len(cc_values) < 1:
        result["trend_summary"] = "No real-space refinement cycles with CC yet"
        return result

    latest_cc = cc_values[-1]

    # Build trend summary string
    if len(cc_values) >= 2:
        trend_parts = ["%0.3f" % cc for cc in cc_values[-5:]]
        previous = cc_values[-2]
        improvement = (latest_cc - previous) / previous * 100 if previous > 0 else 0
        result["improvement_rate"] = improvement

        result["trend_summary"] = "Map-model CC: %s (%+.1f%% last cycle)" % (
            " → ".join(trend_parts), improvement
        )
    else:
        result["trend_summary"] = "Map-model CC: %0.3f (first refinement)" % latest_cc
        return result

    # === STOP CONDITIONS ===

    # 1. SUCCESS: CC above target (was 0.75, now 0.70)
    if latest_cc > 0.70:
        result["should_stop"] = True
        result["reason"] = "SUCCESS: Map-model CC (%.3f) above 0.70 target" % latest_cc
        result["recommendation"] = "stop"
        return result

    # 2. PLATEAU: Less than 0.3% improvement for 3+ cycles (more patience)
    if len(cc_values) >= 4:
        recent_improvements = []
        for i in range(len(cc_values) - 1, max(0, len(cc_values) - 4), -1):
            if cc_values[i-1] > 0:
                imp = (cc_values[i] - cc_values[i-1]) / cc_values[i-1] * 100
                recent_improvements.append(imp)

        if len(recent_improvements) >= 3 and all(imp < 0.3 for imp in recent_improvements):
            result["should_stop"] = True
            result["reason"] = "PLATEAU: CC improvement stalled (<0.3%% for %d cycles)" % len(recent_improvements)
            result["recommendation"] = "stop"
            return result

    # 3. EXCESSIVE: 8+ consecutive real_space_refine (was 5, now more patience)
    if consecutive >= 8:
        result["should_stop"] = True
        result["reason"] = "EXCESSIVE: %d consecutive refinement cycles" % consecutive
        result["recommendation"] = "stop"
        return result

    # 4. WORSENING
    if result["improvement_rate"] < -2.0:
        result["recommendation"] = "consider_stopping"
        result["trend_summary"] += " ⚠️"

    return result


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_latest_resolution(metrics_history):
    """
    Get the most recent resolution value from metrics history.

    Args:
        metrics_history: List of metrics dicts

    Returns:
        float or None: Resolution in Angstroms
    """
    for m in reversed(metrics_history):
        if m.get("resolution"):
            return m["resolution"]
    return None


def get_best_r_free(metrics_history):
    """
    Get the best (lowest) R-free achieved.

    Args:
        metrics_history: List of metrics dicts

    Returns:
        float or None: Best R-free value
    """
    r_free_values = [m["r_free"] for m in metrics_history if m.get("r_free") is not None]
    if r_free_values:
        return min(r_free_values)
    return None


def get_latest_r_free(metrics_history):
    """
    Get the most recent R-free value.

    Args:
        metrics_history: List of metrics dicts

    Returns:
        float or None: Latest R-free value
    """
    for m in reversed(metrics_history):
        if m.get("r_free"):
            return m["r_free"]
    return None


def get_latest_map_cc(metrics_history):
    """
    Get the most recent map-model CC value.

    Args:
        metrics_history: List of metrics dicts

    Returns:
        float or None: Latest CC value
    """
    for m in reversed(metrics_history):
        if m.get("map_cc"):
            return m["map_cc"]
    return None


def format_trend_for_prompt(trend_analysis):
    """
    Format trend analysis for inclusion in LLM prompt.

    Args:
        trend_analysis: Output from analyze_metrics_trend()

    Returns:
        str: Formatted text for prompt
    """
    lines = []

    lines.append("### REFINEMENT PROGRESS")
    lines.append(trend_analysis["trend_summary"])

    if trend_analysis["consecutive_refines"] > 0:
        lines.append("Consecutive X-ray refinement cycles: %d" % trend_analysis["consecutive_refines"])
    if trend_analysis["consecutive_rsr"] > 0:
        lines.append("Consecutive real-space refinement cycles: %d" % trend_analysis["consecutive_rsr"])

    if trend_analysis["recommendation"] == "consider_stopping":
        lines.append("")
        lines.append("⚠️ DIMINISHING RETURNS - Consider stopping or trying different strategy")

    if trend_analysis["should_stop"]:
        lines.append("")
        lines.append("🛑 STOP RECOMMENDED: %s" % trend_analysis["reason"])
        lines.append("Set \"stop\": true in your response.")

    return "\n".join(lines)
