"""
Sanity Checker for PHENIX AI Agent.

Detects impossible or nonsensical states and provides clear error
messages with suggestions for resolution.

The checker evaluates two categories of issues:
- Critical: Fundamental problems that should abort (e.g., experiment type changed)
- Warning: Potential problems worth noting (e.g., resolution unknown)

Usage:
    from libtbx.langchain.agent.sanity_checker import SanityChecker

    checker = SanityChecker()
    result = checker.check(context, session_info)

    if result.should_abort:
        print(result.abort_message)
"""

from __future__ import absolute_import, division, print_function

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class SanityIssue:
    """A single sanity check failure."""
    severity: str          # "critical" or "warning"
    code: str              # e.g., "experiment_type_changed"
    message: str           # Human-readable description
    suggestion: str        # How to fix it
    details: Dict[str, Any] = field(default_factory=dict)  # Additional context

    def to_dict(self):
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


@dataclass
class SanityResult:
    """Result of sanity checking."""
    ok: bool                          # True if no critical issues
    issues: List[SanityIssue]         # All detected issues
    should_abort: bool                # True if should stop execution
    abort_message: Optional[str]      # Formatted message for user

    def to_dict(self):
        """Convert to dictionary for JSON serialization."""
        return {
            "ok": self.ok,
            "issues": [i.to_dict() for i in self.issues],
            "should_abort": self.should_abort,
            "abort_message": self.abort_message,
        }


# =============================================================================
# SANITY CHECKER
# =============================================================================

class SanityChecker:
    """
    Checks for impossible or nonsensical agent states.

    Loads check definitions from workflows.yaml and evaluates them
    against the current context.
    """

    def __init__(self):
        """Initialize the checker and load check definitions."""
        self._checks = {"critical": [], "warning": []}
        self._load_checks()

    def _load_checks(self):
        """Load check definitions from workflows.yaml."""
        try:
            from libtbx.langchain.knowledge.yaml_loader import load_sanity_checks
            self._checks = load_sanity_checks()
        except ImportError:
            try:
                import sys
                import os
                script_dir = os.path.dirname(os.path.abspath(__file__))
                parent_dir = os.path.dirname(script_dir)
                if parent_dir not in sys.path:
                    sys.path.insert(0, parent_dir)
                from libtbx.langchain.knowledge.yaml_loader import load_sanity_checks
                self._checks = load_sanity_checks()
            except (ImportError, Exception):
                # Use defaults if YAML not available
                self._checks = {"critical": [], "warning": []}

    def check(self, context: Dict[str, Any], session_info: Dict[str, Any],
              abort_on_red_flags: bool = True,
              abort_on_warnings: bool = False) -> SanityResult:
        """
        Run all sanity checks against current context.

        Args:
            context: Dict with workflow context including:
                - experiment_type: Current detected type ("xray" or "cryoem")
                - has_model: Whether model files exist
                - has_mtz: Whether MTZ files exist
                - has_map: Whether map files exist
                - state: Current workflow state name
                - history: List of cycle history dicts
                - metrics_history: List of metrics from each cycle
                - categorized_files: Dict of files by category

            session_info: Dict with session state including:
                - experiment_type: Locked experiment type (or None)

            abort_on_red_flags: Whether to abort on critical issues
            abort_on_warnings: Whether to abort on warnings too

        Returns:
            SanityResult with any issues found
        """
        issues = []

        # Run critical checks
        issues.extend(self._run_critical_checks(context, session_info))

        # Run warning checks
        issues.extend(self._run_warning_checks(context, session_info))

        # Determine if we should abort
        critical_issues = [i for i in issues if i.severity == "critical"]
        warning_issues = [i for i in issues if i.severity == "warning"]

        should_abort = False
        if critical_issues and abort_on_red_flags:
            should_abort = True
        if warning_issues and abort_on_warnings:
            should_abort = True

        # Format abort message if needed
        abort_message = None
        if should_abort:
            abort_message = self._format_abort_message(issues)

        return SanityResult(
            ok=len(critical_issues) == 0,
            issues=issues,
            should_abort=should_abort,
            abort_message=abort_message,
        )

    # =========================================================================
    # CRITICAL CHECKS
    # =========================================================================

    def _run_critical_checks(self, context: Dict, session_info: Dict) -> List[SanityIssue]:
        """Run all critical-level checks."""
        issues = []

        # Check: Experiment type stability
        issue = self._check_experiment_type_stable(context, session_info)
        if issue:
            issues.append(issue)

        # Check: Model exists for refinement
        issue = self._check_model_for_refine(context)
        if issue:
            issues.append(issue)

        # Check: Data exists for workflow
        issue = self._check_data_exists(context)
        if issue:
            issues.append(issue)

        # Check: Repeated failures
        issue = self._check_repeated_failures(context)
        if issue:
            issues.append(issue)

        return issues

    def _check_experiment_type_stable(self, context: Dict, session_info: Dict) -> Optional[SanityIssue]:
        """Check that experiment type hasn't changed unexpectedly."""
        initial = session_info.get("experiment_type")
        current = context.get("experiment_type")

        # Only check if both are set
        if initial and current and initial != current:
            return SanityIssue(
                severity="critical",
                code="experiment_type_changed",
                message=f"Experiment type changed from {initial} to {current}",
                suggestion="This usually indicates spurious files were detected. "
                          "Check for intermediate files in subdirectories like run_mr/. "
                          "Remove them and restart, or use session_utils --remove-last to undo recent cycles.",
                details={"initial": initial, "current": current}
            )
        return None

    def _check_model_for_refine(self, context: Dict) -> Optional[SanityIssue]:
        """
        Check that a POSITIONED model exists when entering refine state.

        With semantic categories:
        - 'model' = positioned models (phaser_output, refined, docked, etc.)
        - 'search_model' = templates NOT yet positioned (predicted, pdb_template)

        If user has only search_model but requests refinement, give a specific
        error explaining they need to run Phaser/docking first.
        """
        state = context.get("state", "")

        # Only check if we're in a refinement state
        if "refine" not in state.lower():
            return None

        # Check semantic categories
        has_model = context.get("has_model", False)
        has_search_model = context.get("has_search_model", False)

        # Also check categorized_files for more detail
        categorized = context.get("categorized_files", {})
        model_files = categorized.get("model", [])
        search_model_files = categorized.get("search_model", [])

        # Update has_model/has_search_model from categorized_files if available
        if categorized:
            has_model = len(model_files) > 0
            has_search_model = len(search_model_files) > 0

        if not has_model:
            if has_search_model:
                # User has templates but no positioned model - common mistake!
                exp_type = context.get("experiment_type", "unknown")

                if exp_type == "xray":
                    suggestion = (
                        "You have a search model (predicted structure or template) but no positioned model. "
                        "For X-ray crystallography, you must first run molecular replacement (Phaser) "
                        "to position the model in the unit cell before refinement. "
                        "The search model will be used as input for Phaser."
                    )
                elif exp_type == "cryoem":
                    suggestion = (
                        "You have a search model (predicted structure or template) but no positioned model. "
                        "For cryo-EM, you must first run docking (dock_in_map) to position the model "
                        "in the map before refinement. "
                        "The search model will be used as input for docking."
                    )
                else:
                    suggestion = (
                        "You have a search model (predicted structure or template) but no positioned model. "
                        "You must first run molecular replacement (X-ray) or docking (cryo-EM) "
                        "to position the model before refinement."
                    )

                return SanityIssue(
                    severity="critical",
                    code="search_model_not_positioned",
                    message="Cannot refine: search model found but not yet positioned",
                    suggestion=suggestion,
                    details={
                        "search_model_files": [f.split('/')[-1] for f in search_model_files[:3]],
                        "experiment_type": exp_type,
                        "has_search_model": True,
                        "has_model": False,
                    }
                )
            else:
                # No model at all
                return SanityIssue(
                    severity="critical",
                    code="no_model_for_refine",
                    message="Refinement requested but no model file available",
                    suggestion="Ensure molecular replacement or model building completed successfully. "
                              "Check that output PDB files from previous steps exist. "
                              "For X-ray, you may need to run predict_and_build or Phaser first. "
                              "For cryo-EM, you may need to run predict_and_build or dock_in_map first.",
                )

        return None

    def _check_data_exists(self, context: Dict) -> Optional[SanityIssue]:
        """Check that required experimental data exists for the workflow."""
        exp_type = context.get("experiment_type")
        state = context.get("state", "")

        # Only check after initial state
        if state in ("initial", "unknown", ""):
            return None

        if exp_type == "xray" and not context.get("has_mtz"):
            return SanityIssue(
                severity="critical",
                code="no_data_for_workflow",
                message="No reflection data file found for X-ray workflow",
                suggestion="Provide a reflection data file (MTZ, SCA, or HKL format).",
                details={"experiment_type": exp_type, "missing": "reflection_data"}
            )

        if exp_type == "cryoem" and not context.get("has_map"):
            return SanityIssue(
                severity="critical",
                code="no_data_for_workflow",
                message="No map file found for cryo-EM workflow",
                suggestion="Provide an MRC or CCP4 map file.",
                details={"experiment_type": exp_type, "missing": "map"}
            )

        return None

    def _check_repeated_failures(self, context: Dict) -> Optional[SanityIssue]:
        """Check for repeated identical failures."""
        history = context.get("history", [])

        # Count recent failures by program and error
        failure_counts = {}
        for h in history[-10:]:  # Last 10 cycles
            if not isinstance(h, dict):
                continue
            result = h.get("result", "")
            result_upper = result.upper()

            # Check for actual failures - must start with FAIL or contain ERROR
            # but exclude "without errors" which indicates success
            is_failure = False
            if result_upper.startswith("FAIL"):
                is_failure = True
            elif "ERROR" in result_upper and "WITHOUT ERROR" not in result_upper:
                is_failure = True

            if is_failure:
                prog = h.get("program", "unknown")
                # Extract first line of error for grouping
                error_key = result.split('\n')[0][:100] if result else "unknown"
                key = (prog, error_key)
                failure_counts[key] = failure_counts.get(key, 0) + 1

        # Check for 3+ repeated failures
        for (prog, error), count in failure_counts.items():
            if count >= 3:
                return SanityIssue(
                    severity="critical",
                    code="repeated_failures",
                    message=f"Program {prog} has failed {count} times with similar error",
                    suggestion="Check the error message and input files. "
                              "The program may have a fundamental problem with the provided data. "
                              "Consider using session_utils to examine the session history.",
                    details={"program": prog, "count": count, "error_preview": error[:50]}
                )

        return None

    # =========================================================================
    # WARNING CHECKS
    # =========================================================================

    def _run_warning_checks(self, context: Dict, session_info: Dict) -> List[SanityIssue]:
        """Run all warning-level checks."""
        issues = []

        # Check: Resolution established before refinement
        issue = self._check_resolution_established(context)
        if issue:
            issues.append(issue)

        # Check: Metric anomalies
        anomalies = self._check_metric_anomalies(context)
        if anomalies:
            issues.extend(anomalies)

        return issues

    def _check_resolution_established(self, context: Dict) -> Optional[SanityIssue]:
        """Warn if entering refinement without established resolution."""
        state = context.get("state", "")
        resolution = context.get("resolution")

        if "refine" in state.lower() and resolution is None:
            return SanityIssue(
                severity="warning",
                code="resolution_unknown",
                message="Entering refinement phase without established resolution",
                suggestion="Run xtriage (X-ray) or mtriage (cryo-EM) first to determine data resolution. "
                          "Some programs may fail or produce suboptimal results without resolution.",
            )
        return None

    def _check_metric_anomalies(self, context: Dict) -> List[SanityIssue]:
        """Check for dramatic metric changes that indicate problems."""
        metrics_history = context.get("metrics_history", [])
        issues = []

        if len(metrics_history) < 2:
            return issues

        prev = metrics_history[-2]
        curr = metrics_history[-1]

        # R-free spike (X-ray)
        prev_rfree = prev.get("r_free")
        curr_rfree = curr.get("r_free")
        if prev_rfree and curr_rfree:
            change = curr_rfree - prev_rfree
            if change > 0.15:
                issues.append(SanityIssue(
                    severity="warning",
                    code="r_free_spike",
                    message=f"R-free jumped from {prev_rfree:.3f} to {curr_rfree:.3f} ({change:+.3f})",
                    suggestion="This large increase may indicate the model was corrupted, "
                              "wrong data was used, or refinement diverged. "
                              "Check the input files and consider reverting to a previous cycle.",
                    details={"previous": prev_rfree, "current": curr_rfree, "change": change}
                ))

        # Map CC drop (cryo-EM)
        prev_cc = prev.get("map_cc")
        curr_cc = curr.get("map_cc")
        if prev_cc and curr_cc:
            drop = prev_cc - curr_cc
            if drop > 0.3:
                issues.append(SanityIssue(
                    severity="warning",
                    code="map_cc_drop",
                    message=f"Map CC dropped from {prev_cc:.2f} to {curr_cc:.2f} ({-drop:+.2f})",
                    suggestion="The model may have been displaced from the map. "
                              "Check docking and refinement results. "
                              "Consider reverting to a previous cycle.",
                    details={"previous": prev_cc, "current": curr_cc, "drop": drop}
                ))

        return issues

    # =========================================================================
    # MESSAGE FORMATTING
    # =========================================================================

    def _format_abort_message(self, issues: List[SanityIssue]) -> str:
        """Format a clear, helpful abort message."""
        lines = [
            "",
            "=" * 70,
            "AGENT STOPPED: Sanity check failed",
            "=" * 70,
            "",
        ]

        critical = [i for i in issues if i.severity == "critical"]
        warnings = [i for i in issues if i.severity == "warning"]

        if critical:
            lines.append("CRITICAL ISSUES:")
            lines.append("-" * 40)
            for i, issue in enumerate(critical, 1):
                lines.append(f"  {i}. [{issue.code}]")
                lines.append(f"     {issue.message}")
                lines.append(f"     → {issue.suggestion}")
                lines.append("")

        if warnings:
            lines.append("WARNINGS:")
            lines.append("-" * 40)
            for i, issue in enumerate(warnings, 1):
                lines.append(f"  {i}. [{issue.code}]")
                lines.append(f"     {issue.message}")
                lines.append(f"     → {issue.suggestion}")
                lines.append("")

        lines.extend([
            "=" * 70,
            "TO RESUME AFTER FIXING:",
            "  1. Address the issue(s) described above",
            "  2. Run the agent again with the same session directory",
            "  3. Or use: session_utils --remove-last to undo recent cycles",
            "=" * 70,
            "",
        ])

        return "\n".join(lines)
