"""
Best Files Tracker for PHENIX AI Agent.

Maintains knowledge of the highest-quality file of each type throughout
a session. This enables the agent to always use optimal inputs for each
program, rather than relying on filename patterns or recency alone.

Usage:
    tracker = BestFilesTracker()

    # After each cycle, evaluate new output files
    tracker.evaluate_file("/path/to/refined.pdb", cycle=3,
                         metrics={"r_free": 0.22, "clashscore": 8.5},
                         stage="refined")

    # Get best file for use
    best_model = tracker.get_best("model")
    if best_model:
        print(f"Best model: {best_model.path} (score: {best_model.score})")

    # Persist to session
    data = tracker.to_dict()
    # ... save to session.json ...

    # Restore from session
    tracker = BestFilesTracker.from_dict(data)
"""

from __future__ import absolute_import, division, print_function

import os
import re
from datetime import datetime


# =============================================================================
# DATA CLASSES
# =============================================================================

class BestFileEntry:
    """
    Represents the current best file for a category.

    Attributes:
        path: Full path to the file
        category: Category name (e.g., "model", "map")
        stage: Processing stage (e.g., "refined", "docked")
        score: Computed quality score
        metrics: Metrics at time of selection
        cycle: Cycle number when selected as best
        reason: Human-readable explanation
        timestamp: When this became best
    """

    def __init__(self, path, category, stage=None, score=0.0, metrics=None,
                 cycle=0, reason="", timestamp=None):
        self.path = path
        self.category = category
        self.stage = stage or "unknown"
        self.score = score
        self.metrics = metrics or {}
        self.cycle = cycle
        self.reason = reason
        self.timestamp = timestamp or datetime.now().isoformat()

    def to_dict(self):
        """Convert to dictionary for serialization."""
        return {
            "path": self.path,
            "category": self.category,
            "stage": self.stage,
            "score": self.score,
            "metrics": self.metrics,
            "cycle": self.cycle,
            "reason": self.reason,
            "timestamp": self.timestamp,
        }

    @classmethod
    def from_dict(cls, data):
        """Create from dictionary."""
        if not data:
            return None
        return cls(
            path=data.get("path", ""),
            category=data.get("category", ""),
            stage=data.get("stage"),
            score=data.get("score", 0.0),
            metrics=data.get("metrics"),
            cycle=data.get("cycle", 0),
            reason=data.get("reason", ""),
            timestamp=data.get("timestamp"),
        )

    def __repr__(self):
        return (f"BestFileEntry({self.category}: {os.path.basename(self.path)}, "
                f"score={self.score:.1f}, stage={self.stage})")


class BestFileChange:
    """
    Records a change in the best file for a category.

    Used to track the history of how "best" evolved over the session.
    """

    def __init__(self, category, old_path, new_path, old_score, new_score,
                 cycle, reason, timestamp=None):
        self.category = category
        self.old_path = old_path
        self.new_path = new_path
        self.old_score = old_score
        self.new_score = new_score
        self.cycle = cycle
        self.reason = reason
        self.timestamp = timestamp or datetime.now().isoformat()

    def to_dict(self):
        """Convert to dictionary for serialization."""
        return {
            "category": self.category,
            "old_path": self.old_path,
            "new_path": self.new_path,
            "old_score": self.old_score,
            "new_score": self.new_score,
            "cycle": self.cycle,
            "reason": self.reason,
            "timestamp": self.timestamp,
        }

    @classmethod
    def from_dict(cls, data):
        """Create from dictionary."""
        if not data:
            return None
        return cls(
            category=data.get("category", ""),
            old_path=data.get("old_path", ""),
            new_path=data.get("new_path", ""),
            old_score=data.get("old_score", 0.0),
            new_score=data.get("new_score", 0.0),
            cycle=data.get("cycle", 0),
            reason=data.get("reason", ""),
            timestamp=data.get("timestamp"),
        )

    def __repr__(self):
        return (f"BestFileChange({self.category}: "
                f"{os.path.basename(self.old_path)} -> {os.path.basename(self.new_path)}, "
                f"cycle {self.cycle})")


# =============================================================================
# MAIN TRACKER CLASS
# =============================================================================

class BestFilesTracker:
    """
    Tracks the best file of each category throughout a session.

    Categories tracked:
        - model: Best atomic model (PDB/mmCIF)
        - map: Best full cryo-EM map
        - mtz: Best reflection data (with R-free flags)
        - map_coefficients: Best map coefficients
        - sequence: Sequence file
        - ligand_cif: Ligand restraints

    The tracker uses a scoring system based on:
        - Processing stage (refined > docked > predicted, etc.)
        - Quality metrics (R-free, map_cc, clashscore)
        - Special rules (MTZ: earliest with R-free flags wins forever)

    Scoring configuration is loaded from knowledge/metrics.yaml.
    """

    # Categories we track
    CATEGORIES = [
        "model",
        "map",
        "mtz",
        "map_coefficients",
        "sequence",
        "ligand_cif",
    ]

    def __init__(self):
        """Initialize empty tracker."""
        self.best = {}  # category -> BestFileEntry
        self.history = []  # List of BestFileChange
        self._mtz_with_rfree_locked = False  # Special flag for MTZ handling
        self._scoring_config = None  # Loaded from YAML
        self._load_scoring_config()

    # =========================================================================
    # YAML CONFIGURATION LOADING
    # =========================================================================

    def _load_scoring_config(self):
        """Load scoring configuration from YAML."""
        try:
            # Try to load from YAML
            import yaml
            yaml_path = self._find_yaml_path()
            if yaml_path and os.path.exists(yaml_path):
                with open(yaml_path, 'r') as f:
                    data = yaml.safe_load(f)
                self._scoring_config = data.get("best_files_scoring", {})
                if self._scoring_config:
                    return  # Successfully loaded
        except Exception as e:
            # Log but don't fail
            pass

        # Fall back to defaults
        self._scoring_config = self._get_default_scoring()

    def _find_yaml_path(self):
        """Find the metrics.yaml file."""
        # Try multiple possible locations
        possible_paths = [
            os.path.join(os.path.dirname(__file__), "..", "knowledge", "metrics.yaml"),
            os.path.join(os.path.dirname(__file__), "knowledge", "metrics.yaml"),
            "knowledge/metrics.yaml",
        ]

        # Also try libtbx path
        try:
            import libtbx.load_env
            phenix_path = libtbx.env.find_in_repositories("phenix")
            if phenix_path:
                possible_paths.insert(0, os.path.join(
                    phenix_path, "langchain", "knowledge", "metrics.yaml"))
        except Exception:
            pass

        for path in possible_paths:
            abs_path = os.path.abspath(path)
            if os.path.exists(abs_path):
                return abs_path
        return None

    def _get_default_scoring(self):
        """Return hardcoded default scoring as fallback."""
        return {
            "model": {
                "stage_scores": {
                    "refined": 100,
                    "rsr_output": 100,
                    "autobuild_output": 80,
                    "docked": 60,
                    "processed_predicted": 50,
                    "predicted": 40,
                    "phaser_output": 30,
                    "pdb": 10,
                    "_default": 10,
                },
                "metric_scores": {
                    "r_free": {
                        "max_points": 40,
                        "formula": "linear_inverse",
                        "best_value": 0.20,
                        "worst_value": 0.40,
                    },
                    "map_cc": {
                        "max_points": 30,
                        "formula": "linear",
                        "best_value": 1.0,
                        "worst_value": 0.0,
                    },
                    "clashscore": {
                        "max_points": 30,
                        "formula": "linear_inverse",
                        "best_value": 0,
                        "worst_value": 20,
                    },
                },
            },
            "map": {
                "stage_scores": {
                    "optimized_full_map": 100,
                    "sharpened": 90,
                    "density_modified": 80,
                    "full_map": 50,
                    "map": 40,
                    "half_map": 10,
                    "_default": 40,
                },
                "metric_scores": {
                    "resolution": {
                        "max_points": 30,
                        "formula": "linear_inverse",
                        "best_value": 1.0,
                        "worst_value": 4.0,
                    },
                },
            },
            "mtz": {
                "stage_scores": {
                    "refined_mtz": 70,
                    "original": 50,
                    "mtz": 40,
                    "_default": 40,
                },
                "metric_scores": {
                    "has_rfree_flags": {
                        "max_points": 30,
                        "formula": "boolean",
                    },
                },
                "special_rules": {
                    "lock_on_rfree": True,
                },
            },
            "map_coefficients": {
                "stage_scores": {
                    "refined": 80,
                    "_default": 50,
                },
            },
            "sequence": {
                "stage_scores": {
                    "_default": 50,
                },
            },
            "ligand_cif": {
                "stage_scores": {
                    "user_provided": 60,
                    "_default": 50,
                },
            },
        }

    def _get_stage_score(self, category, stage):
        """Get stage score for a category from config."""
        cat_config = self._scoring_config.get(category, {})
        stage_scores = cat_config.get("stage_scores", {})
        if stage in stage_scores:
            return stage_scores[stage]
        return stage_scores.get("_default", 10)

    def _get_metric_score(self, category, metrics):
        """Calculate metric score for a category from config."""
        if not metrics:
            return 0.0

        cat_config = self._scoring_config.get(category, {})
        metric_configs = cat_config.get("metric_scores", {})

        total = 0.0
        for metric_name, config in metric_configs.items():
            value = metrics.get(metric_name)
            if value is not None:
                total += self._apply_formula(value, config)

        return total

    def _apply_formula(self, value, config):
        """
        Apply a scoring formula to a value.

        Supported formulas:
            - linear: Higher value is better
            - linear_inverse: Lower value is better
            - boolean: True gives max_points, False gives 0

        Args:
            value: The metric value
            config: Dict with formula, max_points, best_value, worst_value

        Returns:
            float: Score contribution (0 to max_points)
        """
        formula = config.get("formula", "linear")
        max_points = config.get("max_points", 0)

        if formula == "boolean":
            return max_points if value else 0

        elif formula == "linear":
            # Higher value is better
            best = config.get("best_value", 1.0)
            worst = config.get("worst_value", 0.0)
            if best == worst:
                return 0
            score = max_points * (value - worst) / (best - worst)
            return max(0, min(max_points, score))

        elif formula == "linear_inverse":
            # Lower value is better
            best = config.get("best_value", 0.0)
            worst = config.get("worst_value", 1.0)
            if worst == best:
                return 0
            score = max_points * (worst - value) / (worst - best)
            return max(0, min(max_points, score))

        return 0

    def _get_special_rules(self, category):
        """Get special rules for a category."""
        cat_config = self._scoring_config.get(category, {})
        return cat_config.get("special_rules", {})

    # =========================================================================
    # PUBLIC API
    # =========================================================================

    def evaluate_file(self, path, cycle, metrics=None, stage=None, category=None):
        """
        Evaluate a file and update best if it's better than current.

        Args:
            path: Full path to the file
            cycle: Current cycle number
            metrics: Quality metrics dict (r_free, map_cc, clashscore, etc.)
            stage: Processing stage (e.g., "refined", "docked")
                   If None, will be inferred from filename
            category: File category (e.g., "model", "map")
                      If None, will be inferred from extension

        Returns:
            bool: True if this file became the new best
        """
        if not path or not os.path.basename(path):
            return False

        # Classify file if category/stage not provided
        if category is None:
            category = self._classify_category(path)
        if category is None:
            return False  # Unknown file type

        if stage is None:
            stage = self._classify_stage(path, category)

        # Skip intermediate/temporary files
        if self._is_intermediate_file(path):
            return False

        # Calculate score
        score = self._calculate_score(path, category, stage, metrics)

        # Special handling for MTZ: earliest with R-free flags wins forever
        if category == "mtz":
            return self._evaluate_mtz(path, cycle, metrics, stage, score)

        # For other categories: higher score wins, with recency as tiebreaker
        return self._evaluate_standard(path, category, cycle, metrics, stage, score)

    def get_best(self, category):
        """
        Get the current best file for a category.

        Args:
            category: Category name (e.g., "model", "map")

        Returns:
            BestFileEntry or None
        """
        return self.best.get(category)

    def get_best_path(self, category):
        """
        Get the path to the best file for a category.

        Args:
            category: Category name

        Returns:
            str or None: Path to best file
        """
        entry = self.best.get(category)
        return entry.path if entry else None

    def get_best_dict(self):
        """
        Get all best files as a simple dict of category -> path.

        Returns:
            dict: {category: path, ...}
        """
        return {cat: entry.path for cat, entry in self.best.items()}

    def get_best_entries(self):
        """
        Get all best file entries.

        Returns:
            dict: {category: BestFileEntry, ...}
        """
        return dict(self.best)

    def get_history(self, category=None):
        """
        Get the history of best file changes.

        Args:
            category: Optional - filter to specific category

        Returns:
            list: List of BestFileChange objects
        """
        if category:
            return [h for h in self.history if h.category == category]
        return list(self.history)

    def get_summary(self):
        """
        Get a human-readable summary of current best files.

        Returns:
            str: Multi-line summary
        """
        lines = ["Best Files:"]
        for category in self.CATEGORIES:
            entry = self.best.get(category)
            if entry:
                lines.append(f"  {category}: {os.path.basename(entry.path)} "
                            f"(score={entry.score:.1f}, stage={entry.stage}, "
                            f"cycle={entry.cycle})")
            else:
                lines.append(f"  {category}: (none)")
        return "\n".join(lines)

    def update_metrics(self, category, metrics, cycle=None):
        """
        Update metrics for an existing best file and recalculate score.

        This is useful when metrics arrive later (e.g., validation after refinement).
        The file remains "best" but its score is updated, which may affect
        future comparisons.

        Args:
            category: Category name (e.g., "model")
            metrics: New metrics dict to merge with existing
            cycle: Optional cycle number for logging

        Returns:
            tuple: (was_updated, old_score, new_score)
        """
        entry = self.best.get(category)
        if not entry:
            return False, 0, 0

        old_score = entry.score

        # Merge new metrics with existing
        merged_metrics = dict(entry.metrics)
        merged_metrics.update(metrics)

        # Recalculate score with merged metrics
        new_score = self._calculate_score(entry.path, category, entry.stage, merged_metrics)

        if new_score != old_score:
            # Update the entry
            entry.metrics = merged_metrics
            entry.score = new_score
            if cycle:
                entry.reason = f"Metrics updated in cycle {cycle}: score {old_score:.1f} -> {new_score:.1f}"

            # Record in history
            change = BestFileChange(
                category=category,
                old_path=entry.path,
                new_path=entry.path,  # Same file, just updated score
                old_score=old_score,
                new_score=new_score,
                cycle=cycle or entry.cycle,
                reason=f"Metrics updated: {list(metrics.keys())}",
            )
            self.history.append(change)

            return True, old_score, new_score

        return False, old_score, new_score

    def update_best_model_metrics(self, metrics, cycle=None):
        """
        Convenience method to update metrics for the best model.

        Args:
            metrics: Dict with model metrics (r_free, map_cc, clashscore, etc.)
            cycle: Optional cycle number

        Returns:
            tuple: (was_updated, old_score, new_score)
        """
        return self.update_metrics("model", metrics, cycle)

    # =========================================================================
    # SERIALIZATION
    # =========================================================================

    def to_dict(self):
        """
        Convert tracker state to dictionary for persistence.

        Returns:
            dict: Serializable state
        """
        return {
            "best": {cat: entry.to_dict() for cat, entry in self.best.items()},
            "history": [h.to_dict() for h in self.history],
            "mtz_with_rfree_locked": self._mtz_with_rfree_locked,
        }

    @classmethod
    def from_dict(cls, data):
        """
        Create tracker from dictionary.

        Args:
            data: Dictionary from to_dict()

        Returns:
            BestFilesTracker instance
        """
        tracker = cls()

        if not data:
            return tracker

        # Restore best entries
        best_data = data.get("best", {})
        for category, entry_data in best_data.items():
            entry = BestFileEntry.from_dict(entry_data)
            if entry:
                tracker.best[category] = entry

        # Restore history
        history_data = data.get("history", [])
        for change_data in history_data:
            change = BestFileChange.from_dict(change_data)
            if change:
                tracker.history.append(change)

        # Restore MTZ lock flag
        tracker._mtz_with_rfree_locked = data.get("mtz_with_rfree_locked", False)

        return tracker

    # =========================================================================
    # SCORING
    # =========================================================================

    def _calculate_score(self, path, category, stage, metrics):
        """
        Calculate quality score for a file using YAML configuration.

        Args:
            path: File path
            category: File category
            stage: Processing stage
            metrics: Quality metrics dict

        Returns:
            float: Score (higher is better)
        """
        # Get stage score from config
        stage_score = self._get_stage_score(category, stage)

        # Get metric score from config
        metric_score = self._get_metric_score(category, metrics)

        return stage_score + metric_score

    # =========================================================================
    # EVALUATION LOGIC
    # =========================================================================

    def _evaluate_standard(self, path, category, cycle, metrics, stage, score):
        """
        Standard evaluation: higher score wins, recency breaks ties.

        Returns:
            bool: True if this file became the new best
        """
        current = self.best.get(category)

        # Determine if this is better
        is_better = False
        reason = ""

        if current is None:
            is_better = True
            reason = f"First {category} file"
        elif score > current.score:
            is_better = True
            reason = f"Higher score ({score:.1f} > {current.score:.1f})"
        elif score == current.score and cycle > current.cycle:
            # Tie-breaker: prefer more recent
            is_better = True
            reason = f"Same score, more recent (cycle {cycle} > {current.cycle})"

        if is_better:
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=current)
            return True

        return False

    def _evaluate_mtz(self, path, cycle, metrics, stage, score):
        """
        Special MTZ evaluation: earliest with R-free flags wins forever.

        Once we have an MTZ with R-free flags, we lock to it and never change.
        This ensures consistent R-free statistics throughout refinement.

        Returns:
            bool: True if this file became the new best
        """
        category = "mtz"
        current = self.best.get(category)
        has_rfree = metrics and metrics.get("has_rfree_flags")

        # If we've locked to an MTZ with R-free, never change
        if self._mtz_with_rfree_locked:
            return False

        # If this MTZ has R-free flags, lock to it
        if has_rfree:
            reason = "First MTZ with R-free flags (locked)"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=current)
            self._mtz_with_rfree_locked = True
            return True

        # No R-free flags - use standard evaluation if we don't have anything yet
        if current is None:
            reason = "First MTZ file (no R-free flags yet)"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=None)
            return True

        return False

    def _update_best(self, path, category, stage, score, metrics, cycle, reason,
                    old_entry=None):
        """Update the best file for a category and record history."""
        # Create new entry
        new_entry = BestFileEntry(
            path=path,
            category=category,
            stage=stage,
            score=score,
            metrics=metrics or {},
            cycle=cycle,
            reason=reason,
        )

        # Record history if this is a change
        if old_entry:
            change = BestFileChange(
                category=category,
                old_path=old_entry.path,
                new_path=path,
                old_score=old_entry.score,
                new_score=score,
                cycle=cycle,
                reason=reason,
            )
            self.history.append(change)

        # Update best
        self.best[category] = new_entry

    # =========================================================================
    # FILE CLASSIFICATION
    # =========================================================================

    def _classify_category(self, path):
        """
        Determine file category from extension.

        Args:
            path: File path

        Returns:
            str or None: Category name
        """
        if not path:
            return None

        lower = path.lower()

        if lower.endswith('.pdb') or lower.endswith('.cif') and 'refine' in lower:
            return "model"
        elif lower.endswith('.mtz'):
            return "mtz"
        elif lower.endswith(('.mrc', '.ccp4', '.map')):
            # Distinguish maps from map coefficients (which are MTZ)
            return "map"
        elif lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            return "sequence"
        elif lower.endswith('.cif') and ('lig' in lower or len(os.path.basename(lower)) < 10):
            return "ligand_cif"

        return None

    def _classify_stage(self, path, category):
        """
        Determine processing stage from filename.

        Args:
            path: File path
            category: File category

        Returns:
            str: Stage name
        """
        basename = os.path.basename(path).lower()

        if category == "model":
            # Check for refined models
            if 'refine' in basename and 'real_space' not in basename:
                return "refined"
            if 'real_space_refined' in basename or 'rsr_' in basename:
                return "rsr_output"
            if 'overall_best' in basename or 'autobuild' in basename:
                return "autobuild_output"
            if 'placed' in basename or 'dock' in basename:
                return "docked"
            if 'processed' in basename:
                return "processed_predicted"
            if 'predict' in basename or 'alphafold' in basename:
                return "predicted"
            if 'phaser' in basename:
                return "phaser_output"
            return "pdb"

        elif category == "map":
            if 'denmod' in basename or 'density_mod' in basename:
                return "optimized_full_map"
            if 'sharp' in basename:
                return "sharpened"
            # Check for half-maps
            if re.search(r'[_-][12ab]\.', basename) or 'half' in basename:
                return "half_map"
            return "full_map"

        elif category == "mtz":
            if 'refine' in basename:
                return "refined_mtz"
            return "original"

        return "unknown"

    def _is_intermediate_file(self, path):
        """
        Check if a file is an intermediate that shouldn't be tracked.

        Args:
            path: File path

        Returns:
            bool: True if intermediate/temporary
        """
        # Patterns for intermediate files
        intermediate_patterns = [
            '/run_mr/',           # dock_in_map intermediate directory
            'run_mr.',            # dock_in_map intermediate files
            '_mr.',               # MR intermediate files
            '/AutoBuild_run_',    # autobuild intermediate directory
            'mask.ccp4',          # mtriage mask output
            '/temp/',             # Temporary directories
            '/tmp/',
            '.tmp.',
        ]

        path_check = path.lower()
        for pattern in intermediate_patterns:
            if pattern.lower() in path_check:
                return True

        return False


# =============================================================================
# MODULE-LEVEL FUNCTIONS
# =============================================================================

def create_tracker():
    """Create a new BestFilesTracker instance."""
    return BestFilesTracker()


# =============================================================================
# TESTING SUPPORT
# =============================================================================

if __name__ == "__main__":
    # Simple self-test
    print("BestFilesTracker self-test")
    print("=" * 50)

    tracker = BestFilesTracker()

    # Simulate a workflow
    print("\n1. First model (predicted):")
    tracker.evaluate_file("/path/to/predicted_model.pdb", cycle=1,
                         stage="predicted")
    print(f"   Best model: {tracker.get_best('model')}")

    print("\n2. Docked model:")
    tracker.evaluate_file("/path/to/placed_model.pdb", cycle=2,
                         stage="docked")
    print(f"   Best model: {tracker.get_best('model')}")

    print("\n3. Refined model with metrics:")
    tracker.evaluate_file("/path/to/refine_001_001.pdb", cycle=3,
                         stage="refined",
                         metrics={"r_free": 0.25, "clashscore": 10})
    print(f"   Best model: {tracker.get_best('model')}")

    print("\n4. Another refined model with better metrics:")
    tracker.evaluate_file("/path/to/refine_002_001.pdb", cycle=4,
                         stage="refined",
                         metrics={"r_free": 0.22, "clashscore": 8})
    print(f"   Best model: {tracker.get_best('model')}")

    print("\n5. Map files:")
    tracker.evaluate_file("/path/to/half_map_1.ccp4", cycle=1, stage="half_map")
    print(f"   Best map after half-map: {tracker.get_best('map')}")
    tracker.evaluate_file("/path/to/denmod_map.ccp4", cycle=2,
                         stage="optimized_full_map",
                         metrics={"resolution": 2.5})
    print(f"   Best map after denmod: {tracker.get_best('map')}")

    print("\n6. MTZ with R-free flags (should lock):")
    tracker.evaluate_file("/path/to/data.mtz", cycle=1,
                         stage="original",
                         metrics={"has_rfree_flags": True})
    print(f"   Best MTZ: {tracker.get_best('mtz')}")

    print("\n7. Try to change MTZ (should fail - locked):")
    result = tracker.evaluate_file("/path/to/refine_001_data.mtz", cycle=3,
                                   stage="refined_mtz",
                                   metrics={"has_rfree_flags": True})
    print(f"   Changed: {result}")
    print(f"   Best MTZ: {tracker.get_best('mtz')}")

    print("\n" + "=" * 50)
    print(tracker.get_summary())

    print("\nHistory:")
    for change in tracker.get_history():
        print(f"  {change}")

    print("\nSerialization test:")
    data = tracker.to_dict()
    tracker2 = BestFilesTracker.from_dict(data)
    print(f"  Restored {len(tracker2.best)} best entries")
    print(f"  Restored {len(tracker2.history)} history entries")

    print("\nSelf-test complete!")
