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
from datetime import datetime

# Centralized pattern utilities - handle both PHENIX and standalone imports
try:
    from libtbx.langchain.agent.pattern_manager import is_half_map
except ImportError:
    from agent.pattern_manager import is_half_map


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
    Tracks the best file of each SEMANTIC category throughout a session.

    Categories tracked:
        - model: Best POSITIONED atomic model (ready for refinement)
        - search_model: Best template for MR/docking (NOT yet positioned)
        - map: Best full cryo-EM map
        - data_mtz: Best reflection data with Fobs/R-free (for refinement)
        - map_coeffs_mtz: Best map coefficients (for ligand fitting/visualization)
        - sequence: Sequence file
        - ligand: Best ligand file (coordinates or restraints)

    IMPORTANT DISTINCTIONS:
        - 'model' = Positioned in experimental reference frame, ready for refinement
        - 'search_model' = Template that needs to be placed via Phaser or dock_in_map
        - 'data_mtz' = Measured structure factors (Fobs, R-free) - for refinement
        - 'map_coeffs_mtz' = Calculated map coefficients (phases) - for ligand fitting

    The tracker uses a scoring system based on:
        - Processing stage (refined > docked > predicted, etc.)
        - Quality metrics (R-free, map_cc, clashscore)
        - Special rules (data_mtz: earliest with R-free flags wins forever)
        - Special rules (map_coeffs_mtz: most recent wins - maps improve)

    Scoring configuration is loaded from knowledge/metrics.yaml.
    """

    # Semantic categories we track
    CATEGORIES = [
        "model",           # Positioned models for refinement
        "search_model",    # Templates for MR/docking (NOT positioned)
        "map",
        "data_mtz",        # Reflection data with Fobs/R-free (for refinement)
        "map_coeffs_mtz",  # Calculated map coefficients (for ligand fitting)
        "sequence",
        "ligand",          # Small molecule coordinates/restraints
    ]

    # Mapping from subcategory/stage to parent semantic category
    STAGE_TO_PARENT = {
        # model subcategories (positioned, ready for refinement)
        "refined": "model",
        "rsr_output": "model",
        "phaser_output": "model",
        "autobuild_output": "model",
        "docked": "model",
        "with_ligand": "model",
        "ligand_fit_output": "model",
        "model_cif": "model",
        # search_model subcategories (templates, NOT positioned)
        "predicted": "search_model",
        "processed_predicted": "search_model",
        "pdb_template": "search_model",
        # ligand subcategories
        "ligand_pdb": "ligand",
        "ligand_cif": "ligand",
        # intermediate - NOT tracked (returns None)
        "intermediate_mr": None,
        "autobuild_temp": None,
        "carryover_temp": None,
    }

    def __init__(self):
        """Initialize empty tracker."""
        self.best = {}  # category -> BestFileEntry
        self.history = []  # List of BestFileChange
        self._data_mtz_with_rfree_locked = False  # Special flag for data_mtz handling
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
            # MODEL: Positioned coordinates ready for refinement
            "model": {
                "stage_scores": {
                    "refined": 100,
                    "rsr_output": 100,
                    "with_ligand": 100,      # Tag, same score as refined
                    "ligand_fit_output": 90,
                    "autobuild_output": 100,  # AutoBuild does internal refinement
                    "phaser_output": 70,      # MR output - positioned in unit cell
                    "docked": 60,             # Docked into map
                    "model_cif": 100,         # mmCIF model from refinement
                    "_default": 50,           # Unknown model type
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
            # SEARCH_MODEL: Templates for MR/docking - NOT yet positioned
            "search_model": {
                "stage_scores": {
                    "processed_predicted": 70,  # Best for MR - trimmed
                    "pdb_template": 60,         # PDB homolog - may need sculptor
                    "predicted": 50,            # Raw prediction - needs processing
                    "_default": 40,
                },
                "metric_scores": {
                    "plddt_mean": {
                        "max_points": 30,
                        "formula": "linear",
                        "best_value": 90,
                        "worst_value": 50,
                    },
                },
            },
            # LIGAND: Small molecule coordinates/restraints
            "ligand": {
                "stage_scores": {
                    "ligand_cif": 60,   # Restraints - preferred
                    "ligand_pdb": 50,   # Coordinates only
                    "_default": 40,
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
                    "intermediate_map": 5,  # resolve_cryo_em initial_map - not for downstream use
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
            # DATA_MTZ: Measured structure factors (Fobs, R-free) for refinement
            "data_mtz": {
                "stage_scores": {
                    "original_data_mtz": 70,    # Input data preserved with R-free
                    "phased_data_mtz": 60,      # Has phases but may lack R-free
                    "data_mtz": 50,             # Generic data MTZ
                    "_default": 40,
                },
                "metric_scores": {
                    "has_rfree_flags": {
                        "max_points": 30,
                        "formula": "boolean",
                    },
                },
                "special_rules": {
                    "lock_on_rfree": True,      # Earliest with R-free wins forever
                },
            },
            # MAP_COEFFS_MTZ: Calculated map coefficients for ligand fitting
            "map_coeffs_mtz": {
                "stage_scores": {
                    "refine_map_coeffs": 80,        # Best quality maps from refine
                    "denmod_map_coeffs": 70,        # Density-modified - excellent for ligand
                    "predict_build_map_coeffs": 60, # From predict_and_build
                    "map_coeffs_mtz": 50,           # Generic
                    "_default": 40,
                },
                "special_rules": {
                    "prefer_recent": True,      # Most recent wins - maps improve
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
            # Backward compatibility alias
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
            stage: Processing stage (e.g., "refined", "docked", "predicted")
                   If None, will be inferred from filename
            category: File category (e.g., "model", "search_model", "map")
                      If None, will be inferred from extension/filename
                      or from stage if stage is provided

        Returns:
            bool: True if this file became the new best
        """
        if not path or not os.path.basename(path):
            return False

        # Skip intermediate/temporary files first
        if self._is_intermediate_file(path):
            return False

        # If stage is provided but category is not, try to infer category from stage
        if category is None and stage is not None:
            category = self.STAGE_TO_PARENT.get(stage)

        # Classify file if category still not known
        if category is None:
            category = self._classify_category(path)
        if category is None:
            return False  # Unknown file type

        # Classify stage if not provided
        if stage is None:
            stage = self._classify_stage(path, category)

        # Calculate score
        score = self._calculate_score(path, category, stage, metrics)

        # Special handling for data_mtz: earliest with R-free flags wins forever
        if category == "data_mtz":
            return self._evaluate_data_mtz(path, cycle, metrics, stage, score)

        # Special handling for map_coeffs_mtz: most recent wins (maps improve)
        if category == "map_coeffs_mtz":
            return self._evaluate_map_coeffs_mtz(path, cycle, metrics, stage, score)

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
            "data_mtz_with_rfree_locked": self._data_mtz_with_rfree_locked,
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

        # Restore data_mtz lock flag (with backward compat for old mtz_with_rfree_locked)
        tracker._data_mtz_with_rfree_locked = data.get(
            "data_mtz_with_rfree_locked",
            data.get("mtz_with_rfree_locked", False)  # Backward compat
        )

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
        # SAFETY: If stage wasn't determined but filename clearly indicates refinement,
        # override stage. This handles cases where program context was lost.
        if category == "model" and stage in (None, "model", "_default"):
            basename = os.path.basename(path).lower()
            if 'refine' in basename and 'real_space' not in basename:
                stage = "refined"
            elif 'rsr_' in basename or '_rsr' in basename or 'real_space_refined' in basename:
                stage = "rsr_output"
            elif 'phaser' in basename:
                stage = "phaser_output"
            elif 'autobuild' in basename or 'overall_best' in basename:
                stage = "autobuild_output"

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

    def _evaluate_data_mtz(self, path, cycle, metrics, stage, score):
        """
        Special data_mtz evaluation: earliest with R-free flags wins forever.

        Once we have a data MTZ with R-free flags, we lock to it and never change.
        This ensures consistent R-free statistics throughout refinement.

        IMPORTANT: Resolution-limited R-free flags are NOT locked, as they
        only cover a subset of the resolution range and would cause problems
        for programs like polder that need full-resolution R-free flags.

        Returns:
            bool: True if this file became the new best
        """
        category = "data_mtz"
        current = self.best.get(category)
        has_rfree = metrics and metrics.get("has_rfree_flags")
        is_resolution_limited = metrics and metrics.get("rfree_resolution_limited")

        # If we've locked to a data_mtz with R-free, never change
        if self._data_mtz_with_rfree_locked:
            return False

        # If this MTZ has R-free flags AND is not resolution-limited, lock to it
        if has_rfree and not is_resolution_limited:
            reason = "First data MTZ with R-free flags (locked)"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=current)
            self._data_mtz_with_rfree_locked = True
            return True

        # If MTZ has R-free flags but is resolution-limited, record but don't lock
        if has_rfree and is_resolution_limited:
            if current is None:
                reason = "Data MTZ with resolution-limited R-free flags (NOT locked - needs full resolution)"
                self._update_best(path, category, stage, score, metrics, cycle, reason,
                                old_entry=None)
                return True
            # If we already have a non-locked MTZ, don't replace with a limited one
            return False

        # No R-free flags - use standard evaluation if we don't have anything yet
        if current is None:
            reason = "First data MTZ file (no R-free flags yet)"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=None)
            return True

        return False

    def _evaluate_map_coeffs_mtz(self, path, cycle, metrics, stage, score):
        """
        Special map_coeffs_mtz evaluation: most recent wins.

        Map coefficients improve as refinement progresses, so we always
        prefer the most recent one (higher cycle number).

        Returns:
            bool: True if this file became the new best
        """
        category = "map_coeffs_mtz"
        current = self.best.get(category)

        # If no current best, this becomes best
        if current is None:
            reason = f"First map coefficients MTZ (cycle {cycle})"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=None)
            return True

        # Prefer more recent (higher cycle) - maps improve with refinement
        if cycle > current.cycle:
            reason = f"More recent map coefficients (cycle {cycle} > {current.cycle})"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=current)
            return True

        # Same cycle: prefer higher score (better stage)
        if cycle == current.cycle and score > current.score:
            reason = f"Better map coefficients at cycle {cycle} (score {score:.1f} > {current.score:.1f})"
            self._update_best(path, category, stage, score, metrics, cycle, reason,
                            old_entry=current)
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
        Determine file's SEMANTIC parent category.

        This returns the high-level category (model, search_model, ligand, etc.)
        based on what the file CAN BE USED FOR, not just its extension.

        Args:
            path: File path

        Returns:
            str or None: Semantic category name (model, search_model, map, etc.)
        """
        if not path:
            return None

        lower = path.lower()
        basename = os.path.basename(lower)

        # Skip intermediate files entirely
        if self._is_intermediate_file(path):
            return None

        # PDB/CIF files need semantic classification
        if lower.endswith('.pdb') or lower.endswith('.cif'):
            return self._classify_pdb_cif_category(path)

        # Other file types
        if lower.endswith('.mtz'):
            # Classify MTZ as data_mtz or map_coeffs_mtz based on filename
            return self._classify_mtz_type(path)
        elif lower.endswith(('.mrc', '.ccp4', '.map')):
            return "map"
        elif lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            return "sequence"

        return None

    def _classify_mtz_type(self, path):
        """
        Classify MTZ as data_mtz or map_coeffs_mtz based on filename patterns.

        Delegates to shared file_utils.classify_mtz_type() for consistency.

        Args:
            path: File path

        Returns:
            str: "data_mtz" or "map_coeffs_mtz"
        """
        try:
            from libtbx.langchain.agent.file_utils import classify_mtz_type
        except ImportError:
            from agent.file_utils import classify_mtz_type
        return classify_mtz_type(path)

    def _classify_pdb_cif_category(self, path):
        """
        Classify a PDB/CIF file into its semantic parent category.

        This is the key method that determines whether a coordinate file
        is a 'model' (positioned, ready for refinement) or a 'search_model'
        (template that needs to be placed first).

        Args:
            path: File path

        Returns:
            str: "model", "search_model", or "ligand"
        """
        basename = os.path.basename(path).lower()
        lower = path.lower()

        # First check for ligands (small molecules)
        if self._is_ligand_file(basename):
            return "ligand"

        # Check for SEARCH_MODEL indicators FIRST (templates, NOT positioned)
        # These take priority because we don't want to accidentally send
        # an unpositioned model to refinement
        search_model_indicators = [
            'predict',          # AlphaFold/ESMFold prediction
            'alphafold',
            'colabfold',
            'esmfold',
            'af-',              # AlphaFold ID prefix
            'template',         # PDB template
            'homolog',          # Homologous structure
            'sculptor',         # Sculptor output (still a template)
            'chainsaw',         # Chainsaw output (still a template)
        ]

        # Also check for "processed" + "model" pattern (common for processed predictions)
        if 'processed' in basename and ('model' in basename or 'predict' in basename):
            # But NOT if it also has phaser/refine
            if 'phaser' not in basename and 'refine' not in basename:
                return "search_model"

        # Exclusions for search_model (these indicate a POSITIONED model)
        positioned_indicators = [
            'phaser',           # PHASER output is positioned
            'refine',           # Refined is positioned
            '_rsr',             # RSR output is positioned
            'rsr_',             # RSR output is positioned
            'placed',           # Docked/placed is positioned
            'dock',             # dock_in_map output is positioned
            'autobuild',        # AutoBuild output is positioned
            'overall_best',     # AutoBuild best is positioned
            'built',            # Built model is positioned
        ]

        # Check search_model indicators
        for indicator in search_model_indicators:
            if indicator in basename:
                # Check if any positioned indicator overrides
                if any(pos in basename for pos in positioned_indicators):
                    return "model"
                return "search_model"

        # Check for MODEL indicators (positioned, ready for refinement)
        model_indicators = [
            'refine',           # Output from phenix.refine
            'rsr_',             # Real-space refined
            'real_space_refined',
            'phaser',           # PHASER output (positioned via MR)
            'placed',           # Docked/placed model
            'dock',             # dock_in_map output
            'autobuild',        # AutoBuild output
            'auto_build',
            'overall_best',     # AutoBuild best model
            'built',            # Built model
            'buccaneer',        # Buccaneer output
            'shelxe',           # SHELXE output
            'with_ligand',      # Model with ligand added
            '_liganded',        # Model with ligand added
            'ligand_fit',       # LigandFit output
            'ligandfit',        # LigandFit output
        ]

        for indicator in model_indicators:
            if indicator in basename:
                return "model"

        # CIF files: check if it's a model or ligand restraints
        if path.lower().endswith('.cif'):
            return self._classify_cif_content(path, basename)

        # Default: assume it's a model (conservative - better to try refinement)
        return "model"

    def _is_ligand_file(self, basename):
        """Check if file is a ligand (small molecule)."""
        ligand_patterns = ['lig.pdb', 'lig.cif', 'ligand.pdb', 'ligand.cif',
                          'lig_', 'ligand_', 'restraint', 'constraint']
        not_ligand = ['ligand_fit', 'ligandfit', 'with_ligand', '_liganded']

        is_ligand = any(p in basename for p in ligand_patterns)
        is_excluded = any(p in basename for p in not_ligand)
        is_small_name = len(basename) < 20

        return is_ligand and not is_excluded and is_small_name

    def _classify_cif_content(self, path, basename):
        """
        Classify CIF file based on naming or content.

        Args:
            path: Full file path
            basename: Lowercase basename

        Returns:
            str: "model", "ligand", or default
        """
        # Naming conventions first (fast)
        if any(p in basename for p in ['restraint', 'constraint', 'lig_']):
            return "ligand"
        if any(p in basename for p in ['refine', 'model', 'coord']):
            return "model"

        # Content detection (slower but accurate)
        try:
            with open(path, 'r') as f:
                content = f.read(4096)  # Read first 4KB
                if '_atom_site.' in content or '_atom_site_' in content:
                    return "model"
                if '_chem_comp.' in content or 'data_comp_' in content:
                    return "ligand"
        except Exception:
            pass

        # Default based on file size (restraints are usually small)
        try:
            if os.path.getsize(path) < 50000:  # < 50KB
                return "ligand"
        except Exception:
            pass

        return "model"  # Default to model

    def _classify_stage(self, path, category):
        """
        Determine processing stage/subcategory from filename.

        For PDB files, this returns the specific subcategory within the
        semantic parent (e.g., "refined" within "model", "predicted" within "search_model").

        Args:
            path: File path
            category: Semantic parent category (model, search_model, ligand, etc.)

        Returns:
            str: Stage/subcategory name
        """
        basename = os.path.basename(path).lower()

        if category == "model":
            # Model subcategories (positioned, ready for refinement)
            # Check more specific patterns first!
            if 'with_ligand' in basename or '_liganded' in basename:
                return "with_ligand"
            if 'ligand_fit' in basename or 'ligandfit' in basename:
                return "ligand_fit_output"
            if 'real_space_refined' in basename or 'rsr_' in basename or '_rsr' in basename:
                return "rsr_output"
            if 'refine' in basename and 'real_space' not in basename:
                return "refined"
            if 'overall_best' in basename or 'autobuild' in basename:
                return "autobuild_output"
            if 'placed' in basename or 'dock' in basename:
                return "docked"
            if 'phaser' in basename:
                return "phaser_output"
            if basename.endswith('.cif') and 'refine' in basename:
                return "model_cif"
            # Default for unknown model type
            return "model"

        elif category == "search_model":
            # Search model subcategories (templates, NOT positioned)
            if 'processed' in basename:
                return "processed_predicted"
            if 'template' in basename or 'homolog' in basename:
                return "pdb_template"
            if 'sculptor' in basename or 'chainsaw' in basename:
                return "pdb_template"
            # Default for predictions
            return "predicted"

        elif category == "ligand":
            # Ligand subcategories
            if basename.endswith('.cif'):
                return "ligand_cif"
            return "ligand_pdb"

        elif category == "map":
            # Skip intermediate maps that shouldn't be used as primary outputs
            if 'initial_map' in basename or 'initial' in basename:
                return "intermediate_map"  # Low priority, won't be selected as best
            if 'denmod' in basename or 'density_mod' in basename:
                return "optimized_full_map"
            if 'sharp' in basename:
                return "sharpened"
            # Check for half-maps using centralized pattern
            if is_half_map(basename):
                return "half_map"
            return "full_map"

        elif category == "data_mtz":
            # Data MTZ subcategories (measured Fobs, R-free)
            if '_data.mtz' in basename:
                return "original_data_mtz"
            if '_refinement.mtz' in basename or 'refinement_data' in basename:
                return "original_data_mtz"
            if 'phased' in basename or 'phases' in basename or 'phaser' in basename:
                return "phased_data_mtz"
            return "data_mtz"

        elif category == "map_coeffs_mtz":
            # Map coefficients MTZ subcategories (calculated phases)
            if 'denmod' in basename or 'density_mod' in basename:
                return "denmod_map_coeffs"
            if 'map_coeffs' in basename and 'predict' in basename:
                return "predict_build_map_coeffs"
            if 'overall_best_map_coeffs' in basename:
                return "predict_build_map_coeffs"
            # Default for refine output
            return "refine_map_coeffs"

        return "unknown"

    def _is_intermediate_file(self, path):
        """
        Check if a file is an intermediate that shouldn't be tracked.

        Args:
            path: File path

        Returns:
            bool: True if intermediate/temporary
        """
        basename = os.path.basename(path)

        # Patterns that indicate VALUABLE output files (never skip these)
        valuable_patterns = [
            '_predicted_model',  # predict_and_build main output
            'overall_best',      # predict_and_build best model
            '_processed',        # processed predicted model
            'with_ligand',       # Model combined with ligand
        ]

        # Check if this is a valuable output first
        if any(pat in basename for pat in valuable_patterns):
            return False  # Not intermediate - track it!

        # Patterns for intermediate files
        intermediate_patterns = [
            '/run_mr/',           # dock_in_map intermediate directory
            'run_mr.',            # dock_in_map intermediate files
            '_mr.',               # MR intermediate files
            '/AutoBuild_run_',    # autobuild intermediate directory
            'mask.ccp4',          # mtriage mask output
            '/temp/',             # Temporary directories
            '/tmp/',
            '/TEMP/',             # LigandFit temp directory
            '/TEMP0/',            # LigandFit temp subdirectory
            '.tmp.',
            '/CarryOn/',          # predict_and_build intermediate directory
            '_CarryOn/',          # predict_and_build intermediate directory (alt)
            'reference',          # Reference/template files
            'EDITED',             # Edited intermediate files
            'superposed_predicted_models',  # Alignment intermediates
            'superposed_predicted_untrimmed',  # Intermediate predictions
            '_ELBOW.',            # Elbow geometry files (not fitted ligands)
            'ELBOW.',             # Elbow geometry files
        ]

        path_check = path.lower()
        basename_check = basename.lower()

        for pattern in intermediate_patterns:
            pattern_lower = pattern.lower()
            if pattern_lower in path_check or pattern_lower in basename_check:
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

    print("\n6. Data MTZ with R-free flags (should lock):")
    tracker.evaluate_file("/path/to/data.mtz", cycle=1,
                         category="data_mtz",
                         stage="original_data_mtz",
                         metrics={"has_rfree_flags": True})
    print(f"   Best data_mtz: {tracker.get_best('data_mtz')}")

    print("\n7. Another data MTZ (should NOT change - locked):")
    result = tracker.evaluate_file("/path/to/refine_001_data.mtz", cycle=3,
                                   category="data_mtz",
                                   stage="original_data_mtz",
                                   metrics={"has_rfree_flags": True})
    print(f"   Changed: {result}")
    print(f"   Best data_mtz: {tracker.get_best('data_mtz')}")

    print("\n8. Map coefficients MTZ (should track):")
    tracker.evaluate_file("/path/to/refine_001_001.mtz", cycle=2,
                         category="map_coeffs_mtz",
                         stage="refine_map_coeffs")
    print(f"   Best map_coeffs_mtz: {tracker.get_best('map_coeffs_mtz')}")

    print("\n9. Newer map coefficients MTZ (should update - prefer recent):")
    tracker.evaluate_file("/path/to/refine_002_001.mtz", cycle=3,
                         category="map_coeffs_mtz",
                         stage="refine_map_coeffs")
    print(f"   Best map_coeffs_mtz: {tracker.get_best('map_coeffs_mtz')}")

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
