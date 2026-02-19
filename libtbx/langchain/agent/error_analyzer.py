"""
Error Analyzer - Automatic Recovery from Recoverable Errors.

This module detects structured errors in PHENIX log output and determines
appropriate recovery strategies. It enables the agent to automatically
recover from certain well-defined error conditions without user intervention.

Currently supported errors:
- ambiguous_data_labels: Multiple data arrays in MTZ file

Usage:
    from libtbx.langchain.agent.error_analyzer import ErrorAnalyzer

    analyzer = ErrorAnalyzer()
    recovery = analyzer.analyze(
        log_text="...",
        program="phenix.autosol",
        context={"project_advice": "MRSAD phasing"},
        session=session
    )

    if recovery:
        # recovery.flags contains the fix
        # recovery.retry_program is the program to re-run
        session.set_recovery_strategy(
            recovery.affected_file,
            recovery.flags,
            recovery.retry_program,
            recovery.reason
        )
        session.data["force_retry_program"] = recovery.retry_program

Design Principles:
1. Extract parameter names from error messages (don't hardcode)
2. Use context (program type, project advice) to make smart selections
3. Track retries to prevent infinite loops
4. Key recovery strategies by filename to avoid cross-contamination
"""

from __future__ import absolute_import, division, print_function

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple

try:
    from libtbx.utils import Sorry
except ImportError:
    class Sorry(Exception):
        pass

# Silence unused import warnings (these are used in type hints)
assert Optional is not None
assert Any is not None
assert Tuple is not None

# YAML loading - use same pattern as other knowledge files
try:
    import yaml
except ImportError:
    yaml = None


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ErrorRecovery:
    """
    Information needed to recover from an error.

    Attributes:
        error_type: Type of error (e.g., "ambiguous_data_labels")
        affected_file: Path to the file causing the error
        flags: Dict of parameter flags to add to command
        reason: Human-readable explanation of the recovery
        retry_program: Program to force-retry
        selected_choice: The option that was selected
        all_choices: All available options
        selected_label: The main label extracted (e.g., "FTOXD3")
        selected_label_pair: The full label pair (e.g., "FTOXD3,SIGFTOXD3")
    """
    error_type: str
    affected_file: str
    flags: Dict[str, str]
    reason: str
    retry_program: str
    selected_choice: str
    all_choices: List[str] = field(default_factory=list)
    selected_label: str = ""
    selected_label_pair: str = ""


# =============================================================================
# ERROR ANALYZER
# =============================================================================

class ErrorAnalyzer:
    """
    Analyzes program errors and determines recovery strategies.

    This class:
    1. Detects recoverable errors in log output
    2. Extracts structured information (choices, keywords)
    3. Applies context-aware resolution logic
    4. Tracks retry attempts to prevent infinite loops

    Configuration is loaded from knowledge/recoverable_errors.yaml.
    """

    def __init__(self):
        """Initialize with configuration from YAML."""
        self._config = self._load_config()
        self._label_patterns = self._config.get("label_patterns", {})
        self._program_prefs = self._config.get("program_data_preferences", {})
        self._context_keywords = self._config.get("context_keywords", {})

    def _load_config(self) -> dict:
        """Load recoverable errors configuration from YAML."""
        if yaml is None:
            print("Warning: PyYAML not available, error recovery disabled")
            return {}

        # Find the knowledge directory
        this_dir = os.path.dirname(os.path.abspath(__file__))

        # Try parent/knowledge (agent/ -> knowledge/)
        knowledge_dir = os.path.join(os.path.dirname(this_dir), "knowledge")
        yaml_path = os.path.join(knowledge_dir, "recoverable_errors.yaml")

        if not os.path.exists(yaml_path):
            # Try sibling directory
            yaml_path = os.path.join(this_dir, "..", "knowledge", "recoverable_errors.yaml")
            yaml_path = os.path.normpath(yaml_path)

        if not os.path.exists(yaml_path):
            print(f"Warning: recoverable_errors.yaml not found")
            return {}

        try:
            with open(yaml_path, 'r') as f:
                return yaml.safe_load(f) or {}
        except Exception as e:
            print(f"Warning: Could not load recoverable_errors.yaml: {e}")
            return {}

    # =========================================================================
    # MAIN API
    # =========================================================================

    def analyze(self, log_text: str, program: str,
                context: Dict[str, Any], session) -> Optional[ErrorRecovery]:
        """
        Analyze log text for recoverable errors.

        This is the main entry point for error analysis. It:
        1. Detects the error type (if any)
        2. Checks retry limits
        3. Extracts structured information
        4. Determines recovery strategy
        5. Updates retry tracking

        Args:
            log_text: Full log/error text from the failed program
            program: Program name that failed (e.g., "phenix.autosol")
            context: Dict containing:
                - project_advice: User's project description
                - history: List of previous cycle records
                - experiment_type: "xray", "cryoem", "sad", etc.
            session: Session object for tracking retries

        Returns:
            ErrorRecovery if a recovery is possible, None otherwise
        """
        if not log_text:
            return None

        # 0. Check for hard-stop errors first — these raise Sorry immediately
        #    and cannot be auto-recovered.
        self._check_hard_stop_errors(log_text, program)

        # 1. Detect error type
        error_type = self._detect_error_type(log_text)
        if not error_type:
            return None

        # 2. Check retry limits
        can_retry, limit_reason = self._check_retry_limits(session, error_type)
        if not can_retry:
            self._log_max_retries(error_type, limit_reason)
            return None

        # 3. Extract structured information
        error_info = self._extract_error_info(log_text, error_type)
        if not error_info:
            return None

        # 4. Determine recovery strategy
        recovery = self._determine_recovery(error_type, error_info, program, context)

        # 5. Update retry tracking in session
        if recovery and session:
            self._update_retry_tracking(
                session,
                error_type,
                recovery.affected_file,
                recovery.selected_choice
            )

        return recovery

    def get_suggestion(self, log_text: str, program: str) -> Optional[str]:
        """
        Get a human-readable suggestion without attempting recovery.

        Used when auto_recovery=False to inform the user what could be done.

        Args:
            log_text: Log/error text
            program: Program that failed

        Returns:
            Human-readable suggestion string, or None
        """
        error_type = self._detect_error_type(log_text)
        if not error_type:
            return None

        error_info = self._extract_error_info(log_text, error_type)
        if not error_info:
            return None

        # Generate suggestion text based on error type
        if error_type == "ambiguous_data_labels":
            keyword = error_info.get("keyword", "obs_labels")
            choices = error_info.get("choices", [])
            if choices:
                # Show first few choices
                choices_preview = choices[:3]
                choices_str = ", ".join(
                    self._extract_main_label(c) for c in choices_preview
                )
                if len(choices) > 3:
                    choices_str += ", ..."
                return (
                    f"Ambiguous data labels detected. "
                    f"Add {keyword}=\"YOUR_CHOICE\" to the command, "
                    f"where YOUR_CHOICE is one of: {choices_str}"
                )

        return None

    # =========================================================================
    # DETECTION
    # =========================================================================

    def _detect_error_type(self, log_text: str) -> Optional[str]:
        """
        Detect which recoverable error type (if any) is present.

        Searches log text for patterns defined in the YAML config.
        Returns the first matching error type.
        """
        errors = self._config.get("errors", {})

        for error_type, error_def in errors.items():
            patterns = error_def.get("detection_patterns", [])
            for pattern in patterns:
                try:
                    # Use DOTALL so .* matches newlines (for multi-line patterns)
                    if re.search(pattern, log_text, re.IGNORECASE | re.DOTALL):
                        return error_type
                except re.error:
                    # Invalid regex pattern in config
                    continue

        return None

    # =========================================================================
    # EXTRACTION
    # =========================================================================

    def _extract_error_info(self, log_text: str,
                            error_type: str) -> Optional[Dict[str, Any]]:
        """
        Extract structured information from error message.

        Dispatches to error-type-specific extraction methods.
        """
        error_def = self._config.get("errors", {}).get(error_type, {})

        if error_type == "ambiguous_data_labels":
            return self._extract_ambiguous_labels_info(log_text, error_def, error_type)
        elif error_type == "ambiguous_experimental_phases":
            return self._extract_ambiguous_labels_info(log_text, error_def, error_type)

        # Future error types would be handled here
        return None

    def _extract_ambiguous_labels_info(self, log_text: str,
                                       error_def: dict,
                                       error_type: str = "ambiguous_data_labels") -> Optional[Dict[str, Any]]:
        """
        Extract keyword and choices from ambiguous data labels error.

        Parses error messages like:
            Multiple equally suitable arrays...
            Possible choices:
              /path/data.mtz:IMEAN,SIGIMEAN
              /path/data.mtz:I(+),SIGI(+),I(-),SIGI(-)
            Please use scaling.input.xray_data.obs_labels
            to specify an unambiguous substring.
        """
        result = {
            "keyword": None,
            "choices": [],
            "affected_file": None,
            "choice_details": []  # [(file_path, labels), ...]
        }

        # Extract keyword name (the parameter to use)
        keyword_pattern = error_def.get("keyword_extraction", "")
        if keyword_pattern:
            try:
                match = re.search(keyword_pattern, log_text, re.IGNORECASE | re.MULTILINE)
                if match:
                    result["keyword"] = match.group(1)
            except (re.error, IndexError):
                pass

        # If no keyword found, try common fallbacks
        if not result["keyword"]:
            # Check for common keywords in the text
            common_keywords = [
                "scaling.input.xray_data.obs_labels",
                "miller_array.labels.name",
                "obs_labels",
                "labels.name",
                "labels",
            ]
            for kw in common_keywords:
                if kw in log_text:
                    result["keyword"] = kw
                    break

        # Extract choices (file:labels pairs)
        choice_pattern = error_def.get("choice_extraction", "")
        if choice_pattern:
            for line in log_text.split('\n'):
                line = line.strip()
                try:
                    match = re.match(choice_pattern, line)
                    if match:
                        file_path = match.group(1)
                        labels = match.group(2).strip()
                        result["choices"].append(labels)
                        result["choice_details"].append((file_path, labels))
                        if not result["affected_file"]:
                            result["affected_file"] = file_path
                except (re.error, IndexError):
                    continue

        # Validation: we need at least choices to proceed
        if not result["choices"]:
            return None

        # Default keyword if still not found - depends on error type
        if not result["keyword"]:
            if error_type == "ambiguous_experimental_phases":
                result["keyword"] = "miller_array.labels.name"
            else:
                result["keyword"] = "obs_labels"

        return result

    # =========================================================================
    # RESOLUTION
    # =========================================================================

    def _check_hard_stop_errors(self, log_text: str, program: str) -> None:
        """
        Check for hard-stop errors that cannot be auto-recovered.

        These are errors that require user action (e.g. supplying a missing
        sequence file). We raise Sorry immediately with a clear actionable
        message rather than silently failing or looping.

        Args:
            log_text: Log/error text from the failed program
            program: Program name that failed

        Raises:
            Sorry: If a hard-stop error is detected
        """
        hard_stops = self._config.get("hard_stop_errors", {})

        for error_name, error_def in hard_stops.items():
            # Check program filter (if specified, only match listed programs)
            allowed_programs = error_def.get("programs", [])
            if allowed_programs and program not in allowed_programs:
                continue

            # Check detection patterns
            patterns = error_def.get("detection_patterns", [])
            for pattern in patterns:
                if re.search(pattern, log_text, re.IGNORECASE | re.DOTALL):
                    message = error_def.get("message", "").strip()
                    raise Sorry(
                        "\n\n[AI Agent] Hard stop: %s\n\n%s" % (
                            error_def.get("description", error_name),
                            message
                        )
                    )

    def _determine_recovery(self, error_type: str, error_info: Dict[str, Any],
                           program: str, context: Dict[str, Any]) -> Optional[ErrorRecovery]:
        """
        Determine the recovery strategy based on error type and context.
        """
        if error_type == "ambiguous_data_labels":
            return self._resolve_ambiguous_labels(error_info, program, context)
        elif error_type == "ambiguous_experimental_phases":
            return self._resolve_ambiguous_phases(error_info, program, context)

        # Future error types would be handled here
        return None

    def _resolve_ambiguous_labels(self, error_info: Dict[str, Any],
                                  program: str,
                                  context: Dict[str, Any]) -> Optional[ErrorRecovery]:
        """
        Resolve ambiguous data labels by selecting appropriate array.

        Selection logic:
        1. Determine if program/context needs anomalous data
        2. Classify each choice as anomalous or merged
        3. Select the appropriate choice
        4. Build the recovery with correct parameter
        """
        choices = error_info.get("choices", [])
        keyword = error_info.get("keyword")
        affected_file = error_info.get("affected_file", "unknown")

        if not choices or not keyword:
            return None

        # Determine if we need anomalous data
        needs_anomalous = self._needs_anomalous_data(program, context)

        # Classify choices
        anomalous_choices = [c for c in choices if self._is_anomalous_label(c)]
        merged_choices = [c for c in choices if self._is_merged_label(c)]

        # Select appropriate choice with clear reasoning
        if needs_anomalous:
            if anomalous_choices:
                selected = anomalous_choices[0]
                reason = f"Selected anomalous data for {program} (phasing workflow)"
            else:
                # No anomalous available - warn but proceed
                selected = choices[0]
                reason = f"WARNING: No anomalous data found for {program}, using first available"
        else:
            if merged_choices:
                selected = merged_choices[0]
                reason = f"Selected merged data for {program}"
            elif anomalous_choices:
                # Only anomalous available - that's okay for most programs
                selected = anomalous_choices[0]
                reason = f"Using anomalous data (only option) for {program}"
            else:
                selected = choices[0]
                reason = f"Using first available data array for {program}"

        # Extract the main label for the flag value
        # E.g., "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-)" -> "I_CuKa(+)"
        label_value = self._extract_main_label(selected)

        return ErrorRecovery(
            error_type="ambiguous_data_labels",
            affected_file=affected_file,
            flags={keyword: label_value},
            reason=reason,
            retry_program=program,
            selected_choice=selected,
            all_choices=choices,
            selected_label=label_value,
            selected_label_pair=selected
        )

    def _resolve_ambiguous_phases(self, error_info: Dict[str, Any],
                                  program: str,
                                  context: Dict[str, Any]) -> Optional[ErrorRecovery]:
        """
        Resolve ambiguous experimental phase labels by selecting appropriate HL coefficients.

        For phenix.refine and other refinement programs, we typically want the
        standard HL coefficients (HLAM,HLBM,HLCM,HLDM) not the anomalous ones
        (HLanomA,HLanomB,HLanomC,HLanomD).

        Selection logic:
        1. Check if this is a phasing program that needs anomalous HL
        2. Otherwise prefer standard (non-anomalous) HL coefficients
        """
        choices = error_info.get("choices", [])
        keyword = error_info.get("keyword")
        affected_file = error_info.get("affected_file", "unknown")

        if not choices or not keyword:
            return None

        # Classify HL coefficient choices
        standard_hl = []
        anomalous_hl = []

        for choice in choices:
            if self._is_anomalous_hl(choice):
                anomalous_hl.append(choice)
            else:
                standard_hl.append(choice)

        # Determine if we need anomalous phases
        needs_anomalous = self._needs_anomalous_data(program, context)

        # Select appropriate choice
        if needs_anomalous:
            if anomalous_hl:
                selected = anomalous_hl[0]
                reason = f"Selected anomalous HL coefficients for {program} (phasing workflow)"
            elif standard_hl:
                selected = standard_hl[0]
                reason = f"Using standard HL coefficients for {program} (no anomalous available)"
            else:
                selected = choices[0]
                reason = f"Using first available HL coefficients for {program}"
        else:
            # For refinement, prefer standard (non-anomalous) HL coefficients
            if standard_hl:
                selected = standard_hl[0]
                reason = f"Selected standard HL coefficients for {program}"
            elif anomalous_hl:
                selected = anomalous_hl[0]
                reason = f"Using anomalous HL coefficients for {program} (only option available)"
            else:
                selected = choices[0]
                reason = f"Using first available HL coefficients for {program}"

        # Extract the main label for the flag value
        label_value = self._extract_main_label(selected)

        return ErrorRecovery(
            error_type="ambiguous_experimental_phases",
            affected_file=affected_file,
            flags={keyword: label_value},
            reason=reason,
            retry_program=program,
            selected_choice=selected,
            all_choices=choices,
            selected_label=label_value,
            selected_label_pair=selected
        )

    def _is_anomalous_hl(self, label: str) -> bool:
        """
        Check if HL coefficient label indicates anomalous phases.

        Anomalous HL typically have 'anom' in the name:
        - HLanomA, HLanomB, HLanomC, HLanomD
        """
        label_lower = label.lower()
        return 'anom' in label_lower

    # =========================================================================
    # LABEL CLASSIFICATION
    # =========================================================================

    def _needs_anomalous_data(self, program: str, context: Dict[str, Any]) -> bool:
        """
        Determine if the program/context needs anomalous data.

        Checks in order:
        1. Program type (autosol, hyss, etc. need anomalous)
        2. Project advice keywords (SAD, MAD, anomalous)
        3. Experiment type
        4. History (previous phasing programs)
        """
        # Normalize program name
        prog_normalized = program.lower().replace("phenix.", "")

        # 1. Check program preference from config
        anomalous_programs = self._program_prefs.get("anomalous", [])
        for ap in anomalous_programs:
            ap_normalized = ap.lower().replace("phenix.", "")
            if prog_normalized == ap_normalized:
                return True

        # 2. Check project advice for keywords
        advice = context.get("project_advice", "").lower()
        anomalous_keywords = self._context_keywords.get("anomalous_workflow", [])
        for kw in anomalous_keywords:
            if kw.lower() in advice:
                return True

        # 3. Check experiment type
        exp_type = context.get("experiment_type", "").lower()
        if exp_type in ["sad", "mad", "mrsad"]:
            return True

        # 4. Check history for phasing programs
        history = context.get("history", [])
        phasing_indicators = ["autosol", "hyss", "phaser_ep", "solve"]
        for entry in history:
            if isinstance(entry, dict):
                prog = (entry.get("program", "") + " " + entry.get("command", "")).lower()
                for indicator in phasing_indicators:
                    if indicator in prog:
                        return True

        return False

    def _is_anomalous_label(self, label: str) -> bool:
        """
        Check if label string indicates anomalous data.

        Matches patterns like I(+), F(-), DANO, *anom, etc.
        """
        patterns = self._label_patterns.get("anomalous_indicators", [])
        label_lower = label.lower()

        for pattern in patterns:
            try:
                if re.search(pattern, label, re.IGNORECASE):
                    return True
            except re.error:
                # Simple string match as fallback
                if pattern.lower() in label_lower:
                    return True

        return False

    def _is_merged_label(self, label: str) -> bool:
        """
        Check if label string indicates merged data.

        Matches patterns like IMEAN, FMEAN, F_obs, etc.
        """
        patterns = self._label_patterns.get("merged_indicators", [])
        label_lower = label.lower()

        for pattern in patterns:
            try:
                if re.search(pattern, label, re.IGNORECASE):
                    return True
            except re.error:
                # Simple string match as fallback
                if pattern.lower() in label_lower:
                    return True

        return False

    def _extract_main_label(self, label_string: str) -> str:
        """
        Extract the main label from a comma-separated label string.

        PHENIX typically wants just the first column name.
        E.g., "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)" -> "I_CuKa(+)"
        """
        if not label_string:
            return label_string

        parts = label_string.split(",")
        return parts[0].strip() if parts else label_string

    # =========================================================================
    # RETRY TRACKING
    # =========================================================================

    def _check_retry_limits(self, session, error_type: str) -> Tuple[bool, Optional[str]]:
        """
        Check if we've exceeded retry limits for this error type.

        Returns (can_retry, reason_if_not).
        """
        if session is None:
            # No session = no tracking = allow retry
            return True, None

        error_def = self._config.get("errors", {}).get(error_type, {})
        max_retries = error_def.get("max_retries", 3)

        attempts = session.data.get("recovery_attempts", {}).get(error_type, {})
        count = attempts.get("count", 0)

        if count >= max_retries:
            files_tried = attempts.get("files_tried", {})
            return False, f"Max recovery attempts ({max_retries}) reached. Tried: {files_tried}"

        return True, None

    def _update_retry_tracking(self, session, error_type: str,
                               affected_file: str, selected_choice: str):
        """
        Update session with new retry attempt.

        Tracks:
        - Total count per error type
        - Which files we've tried
        - Which choices we've selected per file
        """
        if session is None:
            return

        recovery_attempts = session.data.setdefault("recovery_attempts", {})
        type_attempts = recovery_attempts.setdefault(error_type, {
            "count": 0,
            "files_tried": {}
        })

        # Increment count
        type_attempts["count"] = type_attempts.get("count", 0) + 1

        # Track which choices we've tried for this file
        files_tried = type_attempts.setdefault("files_tried", {})
        file_choices = files_tried.setdefault(affected_file, [])
        if selected_choice not in file_choices:
            file_choices.append(selected_choice)

    def _log_max_retries(self, error_type: str, reason: str):
        """Log when max retries reached."""
        print(f"\n{'='*60}")
        print(f"[WARNING] RECOVERY LIMIT REACHED")
        print(f"Error type: {error_type}")
        print(f"Reason: {reason}")
        print(f"The agent will not attempt further automatic recovery.")
        print(f"{'='*60}\n")


# =============================================================================
# MODULE-LEVEL CONVENIENCE FUNCTIONS
# =============================================================================

_analyzer_instance = None

def get_analyzer() -> ErrorAnalyzer:
    """Get singleton ErrorAnalyzer instance."""
    global _analyzer_instance
    if _analyzer_instance is None:
        _analyzer_instance = ErrorAnalyzer()
    return _analyzer_instance


def analyze_error(log_text: str, program: str,
                  context: Dict[str, Any], session) -> Optional[ErrorRecovery]:
    """
    Convenience function to analyze an error.

    Equivalent to ErrorAnalyzer().analyze(...).
    """
    return get_analyzer().analyze(log_text, program, context, session)
