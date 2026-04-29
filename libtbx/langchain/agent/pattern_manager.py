"""
Pattern Manager for PHENIX AI Agent.

Provides centralized regex pattern management with:
- Primitive interpolation ({FLOAT} -> actual regex)
- Named pattern lookup
- Compiled regex caching
- Match/search/findall helpers

Usage:
    from libtbx.langchain.agent.pattern_manager import patterns

    # Match a filename
    match = patterns.match("system.refine_output", "refine_001_002.pdb")
    if match:
        cycle = match.group(1)

    # Extract a metric from log text
    r_free = patterns.extract_float("metrics.r_free", log_text)

    # Get raw regex string for use in other contexts
    regex_str = patterns.get_regex_string("metrics.r_free")

Pattern Definition (patterns.yaml):
    Patterns are defined in knowledge/patterns.yaml with:
    - regex: Pattern string with optional {PRIMITIVE} placeholders
    - description: Human-readable explanation
    - examples: Test cases with expected matches
    - groups: Named description of capture groups

Primitive Interpolation:
    Patterns can reference primitives using {NAME} syntax.
    Example: "R-free{OPT_WS}={OPT_WS}({FLOAT})"
    The PatternManager will replace {FLOAT} and {OPT_WS} with actual regex.
"""

from __future__ import absolute_import, division, print_function

import os
import re

# YAML loading
try:
    import yaml
except ImportError:
    yaml = None


class PatternManager:
    """
    Centralized regex pattern management.

    Loads patterns from patterns.yaml, performs primitive interpolation,
    and provides helpers for matching/searching.
    """

    def __init__(self, patterns_file=None):
        """
        Initialize the pattern manager.

        Args:
            patterns_file: Path to patterns.yaml (default: knowledge/patterns.yaml)
        """
        self._primitives = {}
        self._patterns = {}
        self._compiled_cache = {}
        self._loaded = False

        if patterns_file is None:
            # Try multiple locations for flexibility
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            patterns_file = os.path.join(base_dir, "knowledge", "patterns.yaml")

        self._patterns_file = patterns_file

    def _ensure_loaded(self):
        """Lazy load patterns on first use."""
        if not self._loaded:
            self._load_patterns()
            self._loaded = True

    def _load_patterns(self):
        """Load and parse patterns.yaml."""
        if not yaml:
            # Graceful degradation without PyYAML
            import sys
            print("Warning: PyYAML not available, PatternManager disabled",
                  file=sys.stderr)
            return

        if not os.path.exists(self._patterns_file):
            # Graceful degradation if file doesn't exist
            import sys
            print("Warning: patterns.yaml not found at %s" % self._patterns_file,
                  file=sys.stderr)
            return

        with open(self._patterns_file, 'r') as f:
            data = yaml.safe_load(f)

        if not data:
            return

        # Load primitives first
        self._primitives = {}
        for name, defn in data.get("primitives", {}).items():
            if isinstance(defn, dict):
                self._primitives[name] = defn.get("regex", "")
            else:
                self._primitives[name] = str(defn)

        # Load pattern sections
        for section in ["system", "metrics", "filenames", "directives"]:
            section_data = data.get(section, {})
            if not section_data:
                continue
            for name, defn in section_data.items():
                if isinstance(defn, dict):
                    full_name = "%s.%s" % (section, name)
                    self._patterns[full_name] = defn

    def _interpolate(self, pattern):
        """
        Replace {PRIMITIVE} placeholders with actual regex.

        Example: "R-free{OPT_WS}={OPT_WS}({FLOAT})"
              -> "R-free\\s*=\\s*([-+]?(?:\\d+\\.\\d*|\\d*\\.\\d+|\\d+))"

        Args:
            pattern: Pattern string with {PRIMITIVE} placeholders

        Returns:
            Interpolated regex string
        """
        self._ensure_loaded()

        result = pattern
        # Sort by length descending to handle overlapping names
        # (e.g., {FLOAT_POSITIVE} before {FLOAT})
        for name in sorted(self._primitives.keys(), key=len, reverse=True):
            placeholder = "{%s}" % name
            if placeholder in result:
                result = result.replace(placeholder, self._primitives[name])

        return result

    def get_compiled(self, name, flags=0):
        """
        Get compiled regex for a named pattern.

        Args:
            name: Pattern name (e.g., "metrics.r_free")
            flags: Optional re flags (re.IGNORECASE, etc.)

        Returns:
            Compiled regex pattern

        Raises:
            KeyError: If pattern name not found
        """
        self._ensure_loaded()

        # Check cache
        cache_key = (name, flags)
        if cache_key in self._compiled_cache:
            return self._compiled_cache[cache_key]

        if name not in self._patterns:
            raise KeyError("Unknown pattern: %s" % name)

        defn = self._patterns[name]

        # Use pre-compiled if available, otherwise interpolate
        if "compiled" in defn:
            regex_str = defn["compiled"]
        else:
            regex_str = self._interpolate(defn.get("regex", ""))

        compiled = re.compile(regex_str, flags)
        self._compiled_cache[cache_key] = compiled
        return compiled

    def get_regex_string(self, name):
        """
        Get the raw (interpolated) regex string for a pattern.

        Useful when you need to combine patterns or use them in other contexts.

        Args:
            name: Pattern name

        Returns:
            Interpolated regex string
        """
        self._ensure_loaded()

        if name not in self._patterns:
            raise KeyError("Unknown pattern: %s" % name)

        defn = self._patterns[name]
        if "compiled" in defn:
            return defn["compiled"]
        return self._interpolate(defn.get("regex", ""))

    def has_pattern(self, name):
        """Check if a pattern exists."""
        self._ensure_loaded()
        return name in self._patterns

    def match(self, name, text, flags=0):
        """
        Match pattern at start of text.

        Args:
            name: Pattern name
            text: Text to match against
            flags: Optional re flags

        Returns:
            Match object or None
        """
        pattern = self.get_compiled(name, flags)
        return pattern.match(str(text))

    def search(self, name, text, flags=0):
        """
        Search for pattern anywhere in text.

        Args:
            name: Pattern name
            text: Text to search
            flags: Optional re flags

        Returns:
            Match object or None
        """
        pattern = self.get_compiled(name, flags)
        return pattern.search(str(text))

    def findall(self, name, text, flags=0):
        """
        Find all matches of pattern in text.

        Args:
            name: Pattern name
            text: Text to search
            flags: Optional re flags

        Returns:
            List of matches (strings or tuples)
        """
        pattern = self.get_compiled(name, flags)
        return pattern.findall(str(text))

    def extract(self, name, text, group=1, default=None, flags=0):
        """
        Extract a value from text using a pattern.

        Args:
            name: Pattern name
            text: Text to search
            group: Which group to extract (default: 1)
            default: Value if no match
            flags: Optional re flags

        Returns:
            Extracted string or default
        """
        match = self.search(name, text, flags)
        if match:
            try:
                return match.group(group)
            except (IndexError, AttributeError):
                return default
        return default

    def extract_float(self, name, text, group=1, default=None, flags=0):
        """
        Extract a float value from text.

        Args:
            name: Pattern name
            text: Text to search
            group: Which group to extract
            default: Value if no match or conversion fails
            flags: Optional re flags

        Returns:
            Float value or default
        """
        value = self.extract(name, text, group, None, flags)
        if value is not None:
            try:
                return float(value)
            except (ValueError, TypeError):
                pass
        return default

    def extract_last_float(self, name, text, group=1, default=None, flags=0):
        """
        Extract the LAST float value matching a pattern in text.

        This is important for metrics like R-free which may appear multiple
        times in a refinement log (once per macro cycle). We want the final
        value, not the first.

        Args:
            name: Pattern name
            text: Text to search
            group: Which group to extract
            default: Value if no match or conversion fails
            flags: Optional re flags

        Returns:
            Float value from last match, or default
        """
        pattern = self.get_compiled(name, flags)
        matches = list(pattern.finditer(str(text)))
        if matches:
            last_match = matches[-1]
            try:
                return float(last_match.group(group))
            except (ValueError, TypeError, IndexError):
                pass
        return default

    def extract_int(self, name, text, group=1, default=None, flags=0):
        """
        Extract an integer value from text.

        Args:
            name: Pattern name
            text: Text to search
            group: Which group to extract
            default: Value if no match or conversion fails
            flags: Optional re flags

        Returns:
            Integer value or default
        """
        value = self.extract(name, text, group, None, flags)
        if value is not None:
            try:
                return int(value)
            except (ValueError, TypeError):
                pass
        return default

    def list_patterns(self, section=None):
        """
        List available pattern names.

        Args:
            section: Optional section filter (e.g., "metrics")

        Returns:
            List of pattern names
        """
        self._ensure_loaded()

        if section:
            prefix = "%s." % section
            return [n for n in self._patterns.keys() if n.startswith(prefix)]
        return list(self._patterns.keys())

    def list_sections(self):
        """List available pattern sections."""
        self._ensure_loaded()
        sections = set()
        for name in self._patterns.keys():
            if "." in name:
                sections.add(name.split(".")[0])
        return sorted(sections)

    def get_description(self, name):
        """Get the description for a pattern."""
        self._ensure_loaded()
        if name in self._patterns:
            return self._patterns[name].get("description", "")
        return ""

    def get_examples(self, name):
        """Get example test cases for a pattern."""
        self._ensure_loaded()
        if name in self._patterns:
            return self._patterns[name].get("examples", [])
        return []

    def validate_pattern(self, name):
        """
        Validate a pattern against its test cases.

        Args:
            name: Pattern name

        Returns:
            tuple: (passed, failed) lists of test results
        """
        self._ensure_loaded()

        if name not in self._patterns:
            raise KeyError("Unknown pattern: %s" % name)

        defn = self._patterns[name]
        pattern = self.get_compiled(name)

        passed = []
        failed = []

        # Test matches (should match)
        for test in defn.get("test_matches", []):
            if pattern.search(str(test)):
                passed.append(("match", test))
            else:
                failed.append(("match", test, "should match but didn't"))

        # Test non-matches (should not match)
        for test in defn.get("test_non_matches", []):
            if pattern.search(str(test)):
                failed.append(("non_match", test, "should not match but did"))
            else:
                passed.append(("non_match", test))

        # Test examples with expected groups
        for example in defn.get("examples", []):
            if isinstance(example, dict):
                text = example.get("input", "")
                match = pattern.search(str(text))
                if not match:
                    failed.append(("example", text, "no match"))
                else:
                    for key, expected in example.items():
                        if key.startswith("group_"):
                            try:
                                group_num = int(key.split("_")[1])
                                actual = match.group(group_num)
                                if actual != expected:
                                    failed.append(("example", text,
                                        "group %d: expected '%s', got '%s'" %
                                        (group_num, expected, actual)))
                                else:
                                    passed.append(("example", text, "group %d" % group_num))
                            except (IndexError, ValueError) as e:
                                failed.append(("example", text, str(e)))
            elif isinstance(example, str):
                # Simple example - just check it matches
                if pattern.search(example):
                    passed.append(("example", example))
                else:
                    failed.append(("example", example, "no match"))

        return passed, failed

    def validate_all(self):
        """
        Validate all patterns against their test cases.

        Returns:
            dict: {pattern_name: (passed, failed)} for each pattern
        """
        self._ensure_loaded()

        results = {}
        for name in self._patterns.keys():
            results[name] = self.validate_pattern(name)
        return results

    def get_primitive(self, name):
        """Get a primitive regex by name."""
        self._ensure_loaded()
        return self._primitives.get(name, "")

    def list_primitives(self):
        """List available primitive names."""
        self._ensure_loaded()
        return list(self._primitives.keys())


# =============================================================================
# UTILITY FUNCTIONS - Common operations used across the codebase
# =============================================================================

def extract_cycle_number(filename, default=0):
    """
    Extract cycle number from a filename using unified logic.

    This is the SINGLE SOURCE OF TRUTH for cycle extraction.
    All code that needs to extract cycle numbers should use this function
    to ensure consistent behavior across Planner, CommandBuilder, and GraphNodes.

    Priority order:
    1. Program prefix (refine_001, rsr_5) - most reliable
    2. Run directory (run_001) - if present in path
    3. Trailing number (model_001.pdb) - fallback

    Args:
        filename: Filename or path to extract cycle from
        default: Value to return if no cycle found (default: 0)

    Returns:
        int: Extracted cycle number, or default if not found

    Examples:
        >>> extract_cycle_number("refine_003_001.pdb")
        3
        >>> extract_cycle_number("rsr_5_output.pdb")
        5
        >>> extract_cycle_number("model_001.pdb")
        1
        >>> extract_cycle_number("run_2/output.pdb")
        2
    """
    import os
    pm = get_patterns()

    # Get basename for most checks
    basename = os.path.basename(str(filename)).lower()

    # Priority 1: Check for program prefix (most reliable)
    cycle = pm.extract_int("system.cycle_from_prefix", basename, group=1)
    if cycle is not None:
        return cycle

    # Priority 2: Check for run directory in path
    full_path = str(filename).lower()
    cycle = pm.extract_int("system.cycle_from_run_dir", full_path, group=1)
    if cycle is not None:
        return cycle

    # Priority 3: Fallback to trailing number
    cycle = pm.extract_int("system.cycle_from_suffix", basename, group=1)
    if cycle is not None:
        return cycle

    return default


def extract_all_numbers(filename):
    """
    Extract all numbers from a filename as a tuple for secondary sorting.

    This is used when files have the same primary cycle number but differ
    in secondary numbers (e.g., rsr_001_..._000.pdb vs rsr_001_..._002.pdb).

    Args:
        filename: Filename or path

    Returns:
        tuple: All numbers found, as integers, in order

    Examples:
        >>> extract_all_numbers("rsr_001_refined_002.pdb")
        (1, 2)
        >>> extract_all_numbers("refine_003_001.pdb")
        (3, 1)
        >>> extract_all_numbers("model.pdb")
        (0,)
    """
    import os
    basename = os.path.basename(str(filename))
    numbers = re.findall(r'(\d+)', basename)
    if numbers:
        return tuple(int(n) for n in numbers)
    return (0,)


def is_half_map(filename):
    """
    Check if a filename indicates a half-map.

    This is the SINGLE SOURCE OF TRUTH for half-map detection.

    Args:
        filename: Filename to check

    Returns:
        bool: True if filename indicates a half-map
    """
    pm = get_patterns()
    basename = os.path.basename(str(filename)).lower()

    # Use the centralized pattern
    if pm.has_pattern("system.half_map"):
        return pm.search("system.half_map", basename) is not None

    # Fallback if pattern not defined (shouldn't happen)
    return bool(re.search(r'[_-][12ab]\.', basename) or 'half' in basename)


# =============================================================================
# Global Instance (Lazy Loaded)
# =============================================================================

_patterns = None


def get_patterns():
    """Get the global PatternManager instance."""
    global _patterns
    if _patterns is None:
        _patterns = PatternManager()
    return _patterns


class _PatternsProxy:
    """
    Proxy that lazily loads the PatternManager.

    This allows `from libtbx.langchain.agent.pattern_manager import patterns` to work
    without immediately loading the YAML file.
    """

    def __getattr__(self, name):
        return getattr(get_patterns(), name)

    def __repr__(self):
        return "<PatternManager proxy>"


# Convenience instance - can do `from libtbx.langchain.agent.pattern_manager import patterns`
patterns = _PatternsProxy()


# =============================================================================
# Validation Script (when run directly)
# =============================================================================

def main():
    """Validate all patterns when run as script."""
    pm = get_patterns()

    print("=" * 60)
    print("PATTERN VALIDATION")
    print("=" * 60)

    # Validate primitives first
    print("\n--- Primitives ---")
    for name in pm.list_primitives():
        regex = pm.get_primitive(name)
        try:
            re.compile(regex)
            print("  ✓ %s" % name)
        except re.error as e:
            print("  ❌ %s: %s" % (name, e))

    # Validate patterns by section
    all_passed = True
    total_tests = 0
    total_failed = 0

    for section in pm.list_sections():
        print("\n--- %s ---" % section.capitalize())
        for name in pm.list_patterns(section):
            passed, failed = pm.validate_pattern(name)
            total_tests += len(passed) + len(failed)
            total_failed += len(failed)

            if failed:
                print("  ❌ %s:" % name)
                for f in failed:
                    print("     %s" % str(f))
                all_passed = False
            else:
                test_count = len(passed)
                if test_count > 0:
                    print("  ✓ %s (%d tests)" % (name, test_count))
                else:
                    print("  ○ %s (no tests)" % name)

    print("\n" + "=" * 60)
    print("SUMMARY: %d tests, %d failed" % (total_tests, total_failed))
    print("=" * 60)

    return 0 if all_passed else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
