"""
Metric Evaluator for PHENIX AI Agent.

This module provides YAML-driven metrics evaluation:
- Quality assessment (is model good/acceptable/poor?)
- Improvement detection (is metric improving?)
- Target calculations (resolution-dependent thresholds)
- Plateau detection (has refinement stalled?)

Reads definitions from metrics.yaml instead of hardcoding thresholds.

Usage:
    from libtbx.langchain.agent.metric_evaluator import MetricEvaluator

    evaluator = MetricEvaluator()

    # Check if R-free is acceptable
    is_ok = evaluator.is_acceptable("r_free", 0.25, resolution=2.0)

    # Get quality assessment
    quality = evaluator.assess_quality({"r_free": 0.25, "clashscore": 3})

    # Analyze trend
    trend = evaluator.analyze_trend(metrics_history, "xray")
"""

from __future__ import absolute_import, division, print_function

# Import YAML loader
from libtbx.langchain.knowledge.yaml_loader import (
    get_metric_threshold,
    get_metric_direction,
    is_metric_good,
    is_metric_acceptable,
    load_metrics,
)


class MetricEvaluator:
    """
    YAML-driven metric evaluation.

    Provides quality assessment, improvement detection, and
    trend analysis using thresholds from metrics.yaml.
    """

    def __init__(self):
        """Initialize the evaluator."""
        # Load improvement thresholds from YAML
        metrics = load_metrics()
        improvement_config = metrics.get("improvement", {})
        self.significant_change = improvement_config.get("significant_change", {})
        self.plateau_config = improvement_config.get("plateau", {})

    # =========================================================================
    # THRESHOLD ACCESS
    # =========================================================================

    def get_threshold(self, metric_name, level="acceptable", resolution=None):
        """
        Get threshold for a metric.

        Args:
            metric_name: Name of metric (e.g., "r_free")
            level: "good", "acceptable", or "poor"
            resolution: Optional resolution for resolution-dependent thresholds

        Returns:
            float: Threshold value or None
        """
        return get_metric_threshold(metric_name, level, resolution)

    def get_target(self, metric_name, resolution=None):
        """
        Get target value for a metric (acceptable threshold).

        Args:
            metric_name: Name of metric
            resolution: Optional resolution

        Returns:
            float: Target value
        """
        return self.get_threshold(metric_name, "acceptable", resolution)

    def get_direction(self, metric_name):
        """
        Get whether higher or lower is better.

        Args:
            metric_name: Name of metric

        Returns:
            str: "minimize" or "maximize"
        """
        return get_metric_direction(metric_name)

    # =========================================================================
    # QUALITY ASSESSMENT
    # =========================================================================

    def is_good(self, metric_name, value, resolution=None):
        """
        Check if metric value is "good".

        Args:
            metric_name: Name of metric
            value: The value to check
            resolution: Optional resolution

        Returns:
            bool: True if good
        """
        return is_metric_good(metric_name, value, resolution)

    def is_acceptable(self, metric_name, value, resolution=None):
        """
        Check if metric value is "acceptable".

        Args:
            metric_name: Name of metric
            value: The value to check
            resolution: Optional resolution

        Returns:
            bool: True if acceptable (includes good)
        """
        return is_metric_acceptable(metric_name, value, resolution)

    def get_quality_level(self, metric_name, value, resolution=None):
        """
        Get quality level for a metric value.

        Args:
            metric_name: Name of metric
            value: The value to assess
            resolution: Optional resolution

        Returns:
            str: "good", "acceptable", "poor", or "unknown"
        """
        if value is None:
            return "unknown"

        if self.is_good(metric_name, value, resolution):
            return "good"
        elif self.is_acceptable(metric_name, value, resolution):
            return "acceptable"
        else:
            return "poor"

    def assess_quality(self, metrics, resolution=None):
        """
        Assess overall model quality from multiple metrics.

        Args:
            metrics: Dict of metric_name -> value
            resolution: Optional resolution

        Returns:
            dict: {
                overall: "good" | "acceptable" | "poor",
                details: {metric_name: quality_level, ...},
                issues: [list of problems],
            }
        """
        details = {}
        issues = []

        # Assess each metric
        for name, value in metrics.items():
            if value is None:
                continue

            level = self.get_quality_level(name, value, resolution)
            details[name] = level

            if level == "poor":
                threshold = self.get_threshold(name, "acceptable", resolution)
                direction = self.get_direction(name)
                if direction == "minimize":
                    issues.append("%s (%.3f) above acceptable threshold (%.3f)" % (name, value, threshold))
                else:
                    issues.append("%s (%.3f) below acceptable threshold (%.3f)" % (name, value, threshold))

        # Determine overall quality
        if not details:
            overall = "unknown"
        elif all(level == "good" for level in details.values()):
            overall = "good"
        elif all(level in ("good", "acceptable") for level in details.values()):
            overall = "acceptable"
        else:
            overall = "poor"

        return {
            "overall": overall,
            "details": details,
            "issues": issues,
        }

    # =========================================================================
    # IMPROVEMENT DETECTION
    # =========================================================================

    def is_significant_improvement(self, metric_name, old_value, new_value):
        """
        Check if change represents significant improvement.

        Args:
            metric_name: Name of metric
            old_value: Previous value
            new_value: Current value

        Returns:
            bool: True if significant improvement
        """
        if old_value is None or new_value is None:
            return False

        direction = self.get_direction(metric_name)
        threshold = self.significant_change.get(metric_name, 0.005)

        if direction == "minimize":
            # Lower is better - improvement means new < old
            change = old_value - new_value
        else:
            # Higher is better - improvement means new > old
            change = new_value - old_value

        return change >= threshold

    def calculate_improvement_rate(self, metric_name, old_value, new_value):
        """
        Calculate percentage improvement.

        Args:
            metric_name: Name of metric
            old_value: Previous value
            new_value: Current value

        Returns:
            float: Percentage improvement (positive = better)
        """
        if old_value is None or new_value is None or old_value == 0:
            return 0.0

        direction = self.get_direction(metric_name)

        if direction == "minimize":
            # Lower is better
            return (old_value - new_value) / old_value * 100
        else:
            # Higher is better
            return (new_value - old_value) / old_value * 100

    def is_plateau(self, values, metric_name=None):
        """
        Check if recent values indicate a plateau.

        Args:
            values: List of recent values
            metric_name: Optional metric name for threshold lookup

        Returns:
            bool: True if plateau detected
        """
        cycles_needed = self.plateau_config.get("cycles", 3)
        threshold_pct = self.plateau_config.get("threshold", 0.003) * 100  # Convert to percentage

        if len(values) < cycles_needed + 1:
            return False

        # Calculate improvements for last N cycles
        recent = values[-(cycles_needed + 1):]
        improvements = []

        direction = self.get_direction(metric_name) if metric_name else "minimize"

        for i in range(1, len(recent)):
            if recent[i-1] == 0:
                continue
            if direction == "minimize":
                imp = (recent[i-1] - recent[i]) / recent[i-1] * 100
            else:
                imp = (recent[i] - recent[i-1]) / recent[i-1] * 100
            improvements.append(imp)

        if len(improvements) < cycles_needed:
            return False

        # Plateau if all recent improvements are below threshold
        return all(imp < threshold_pct for imp in improvements[-cycles_needed:])

    # =========================================================================
    # TREND ANALYSIS
    # =========================================================================

    def analyze_trend(self, metrics_history, experiment_type="xray", resolution=None):
        """
        Analyze metrics trend to detect plateaus and recommend stopping.

        This is a YAML-driven replacement for analyze_metrics_trend().

        Args:
            metrics_history: List of metrics dicts
            experiment_type: "xray" or "cryoem"
            resolution: Optional resolution

        Returns:
            dict: Trend analysis result (compatible with metrics_analyzer output)
        """
        result = {
            "should_stop": False,
            "reason": None,
            "trend_summary": "No refinement history",
            "r_free_trend": [],
            "map_cc_trend": [],
            "improvement_rate": 0.0,
            "consecutive_refines": 0,
            "consecutive_rsr": 0,
            "recommendation": "continue",
        }

        if not metrics_history:
            return result

        if experiment_type == "cryoem":
            return self._analyze_cryoem_trend(metrics_history, result)
        else:
            return self._analyze_xray_trend(metrics_history, resolution, result)

    def _analyze_xray_trend(self, metrics_history, resolution, result):
        """Analyze X-ray refinement trend."""

        # Count consecutive refines
        consecutive = 0
        for m in reversed(metrics_history):
            prog = m.get("program", "").lower()
            has_r_free = m.get("r_free") is not None
            if "refine" in prog and "real_space" not in prog:
                consecutive += 1
            elif prog == "unknown" and has_r_free:
                consecutive += 1
            else:
                break
        result["consecutive_refines"] = consecutive

        # Extract R-free values
        r_free_values = []
        for m in metrics_history:
            r_free = m.get("r_free")
            if r_free is not None:
                prog = m.get("program", "").lower()
                if "refine" in prog or prog == "unknown":
                    r_free_values.append(r_free)

        result["r_free_trend"] = r_free_values[-5:]

        if not r_free_values:
            result["trend_summary"] = "No refinement cycles with R-free yet"
            return result

        latest_r_free = r_free_values[-1]

        # Get target from YAML
        target = self.get_target("r_free", resolution) or 0.25
        success_threshold = target - 0.02

        # Check for validation
        validation_done = any(
            m.get("program") in ("phenix.molprobity", "phenix.model_vs_data")
            for m in metrics_history
        )

        # SUCCESS check
        if latest_r_free < success_threshold:
            if validation_done:
                result["should_stop"] = True
                result["reason"] = "SUCCESS: R-free (%.3f) below target (%.2f)" % (latest_r_free, target)
                result["recommendation"] = "stop"
            else:
                result["should_stop"] = False
                result["reason"] = "R-free below target - recommend validation"
                result["recommendation"] = "validate"
                result["suggest_validation"] = True
            result["trend_summary"] = "R-free: %.3f - TARGET REACHED" % latest_r_free
            return result

        # Build trend summary
        if len(r_free_values) >= 2:
            previous = r_free_values[-2]
            improvement = self.calculate_improvement_rate("r_free", previous, latest_r_free)
            result["improvement_rate"] = improvement

            trend_parts = ["%.3f" % r for r in r_free_values[-5:]]
            if improvement > 1.0:
                direction = "improving"
            elif improvement > 0:
                direction = "slightly improving"
            elif improvement > -1.0:
                direction = "stable"
            else:
                direction = "WORSENING"

            result["trend_summary"] = "R-free: %s (%s, %+.1f%%)" % (
                " → ".join(trend_parts), direction, improvement
            )
        else:
            result["trend_summary"] = "R-free: %.3f (first refinement)" % latest_r_free
            return result

        # PLATEAU check using YAML thresholds
        if self.is_plateau(r_free_values, "r_free"):
            result["should_stop"] = True
            result["reason"] = "PLATEAU: R-free improvement stalled"
            result["recommendation"] = "stop"
            return result

        # EXCESSIVE refinement check
        max_cycles = 8
        if consecutive >= max_cycles:
            result["should_stop"] = True
            result["reason"] = "EXCESSIVE: %d consecutive refinement cycles" % consecutive
            result["recommendation"] = "stop"
            return result

        # WORSENING check
        if result["improvement_rate"] < -2.0:
            result["recommendation"] = "consider_stopping"
            result["trend_summary"] += " ⚠️"

        return result

    def _analyze_cryoem_trend(self, metrics_history, result):
        """Analyze cryo-EM refinement trend."""

        # Count consecutive RSR
        consecutive = 0
        for m in reversed(metrics_history):
            prog = m.get("program", "").lower()
            if "real_space" in prog:
                consecutive += 1
            else:
                break
        result["consecutive_rsr"] = consecutive

        # Extract CC values
        cc_values = []
        for m in metrics_history:
            cc = m.get("map_cc")
            if cc is not None:
                prog = m.get("program", "").lower()
                if "real_space" in prog:
                    cc_values.append(cc)

        result["map_cc_trend"] = cc_values[-5:]

        if not cc_values:
            result["trend_summary"] = "No real-space refinement with CC yet"
            return result

        latest_cc = cc_values[-1]

        # Get target from YAML
        target = self.get_target("map_cc") or 0.70

        # SUCCESS check
        if latest_cc > target:
            result["should_stop"] = True
            result["reason"] = "SUCCESS: Map-model CC (%.3f) above target (%.2f)" % (latest_cc, target)
            result["recommendation"] = "stop"
            result["trend_summary"] = "Map CC: %.3f - TARGET REACHED" % latest_cc
            return result

        # Build trend summary
        if len(cc_values) >= 2:
            previous = cc_values[-2]
            improvement = self.calculate_improvement_rate("map_cc", previous, latest_cc)
            result["improvement_rate"] = improvement

            trend_parts = ["%.3f" % cc for cc in cc_values[-5:]]
            result["trend_summary"] = "Map CC: %s (%+.1f%%)" % (
                " → ".join(trend_parts), improvement
            )
        else:
            result["trend_summary"] = "Map CC: %.3f (first refinement)" % latest_cc
            return result

        # PLATEAU check
        if self.is_plateau(cc_values, "map_cc"):
            result["should_stop"] = True
            result["reason"] = "PLATEAU: CC improvement stalled"
            result["recommendation"] = "stop"
            return result

        # EXCESSIVE refinement
        if consecutive >= 8:
            result["should_stop"] = True
            result["reason"] = "EXCESSIVE: %d consecutive RSR cycles" % consecutive
            result["recommendation"] = "stop"
            return result

        return result


# =============================================================================
# SINGLETON AND CONVENIENCE FUNCTIONS
# =============================================================================

_evaluator = None

def get_evaluator():
    """Get global MetricEvaluator instance."""
    global _evaluator
    if _evaluator is None:
        _evaluator = MetricEvaluator()
    return _evaluator


def assess_model_quality(metrics, resolution=None):
    """Convenience function to assess model quality."""
    return get_evaluator().assess_quality(metrics, resolution)


def analyze_refinement_trend(metrics_history, experiment_type="xray", resolution=None):
    """Convenience function to analyze refinement trend."""
    return get_evaluator().analyze_trend(metrics_history, experiment_type, resolution)


# =============================================================================
# TESTING
# =============================================================================

if __name__ == "__main__":
    print("Testing MetricEvaluator...")
    print()

    evaluator = MetricEvaluator()

    # Test thresholds
    print("R-free thresholds:")
    for res in [1.5, 2.0, 2.5, 3.0]:
        target = evaluator.get_target("r_free", res)
        print("  %.1fÅ: %.3f" % (res, target or 0))
    print()

    # Test quality assessment
    print("Quality assessment:")
    metrics = {"r_free": 0.24, "clashscore": 3.5}
    quality = evaluator.assess_quality(metrics, resolution=2.0)
    print("  Metrics:", metrics)
    print("  Overall:", quality["overall"])
    print("  Details:", quality["details"])
    print()

    # Test improvement detection
    print("Improvement detection:")
    print("  R-free 0.28 -> 0.27: significant?",
          evaluator.is_significant_improvement("r_free", 0.28, 0.27))
    print("  R-free 0.28 -> 0.279: significant?",
          evaluator.is_significant_improvement("r_free", 0.28, 0.279))
    print()

    # Test plateau detection
    print("Plateau detection:")
    values = [0.30, 0.29, 0.289, 0.288, 0.287]
    print("  Values:", values)
    print("  Is plateau?", evaluator.is_plateau(values, "r_free"))
    print()

    # Test trend analysis
    print("Trend analysis:")
    history = [
        {"program": "phenix.refine", "r_free": 0.35},
        {"program": "phenix.refine", "r_free": 0.30},
        {"program": "phenix.refine", "r_free": 0.27},
    ]
    trend = evaluator.analyze_trend(history, "xray", resolution=2.0)
    print("  Summary:", trend["trend_summary"])
    print("  Recommendation:", trend["recommendation"])
