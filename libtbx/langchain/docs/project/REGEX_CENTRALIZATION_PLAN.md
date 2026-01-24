# Implementation Plan: Regex Pattern Centralization

## Executive Summary

Currently, regex patterns are scattered throughout the codebase - in Python files, YAML configs, and inline code. This makes them:
- **Fragile**: A typo in one place doesn't affect others
- **Hard to understand**: Raw regex like `[-+]?(?:\d*\.\d+|\d+)` is cryptic
- **Inconsistent**: Different definitions of "float" in different places
- **Hard to test**: No central place to verify patterns work

This plan introduces **Semantic Pattern Interpolation** - a system where:
1. Primitives (`FLOAT`, `INT`, `WHITESPACE`) are defined once in `patterns.yaml`
2. Complex patterns compose primitives: `R-free{OPT_WS}[=:]{OPT_WS}({FLOAT})`
3. A `PatternManager` handles loading, interpolation, and matching
4. All code uses named patterns instead of inline regex

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    knowledge/patterns.yaml                       │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐  │
│  │ primitives  │  │   system    │  │        metrics          │  │
│  │  FLOAT      │  │ cycle_det   │  │  r_factor, resolution   │  │
│  │  INT        │  │ filename    │  │  map_cc, clashscore     │  │
│  │  OPT_WS     │  │ traceback   │  │                         │  │
│  └─────────────┘  └─────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  agent/pattern_manager.py                        │
│                                                                  │
│  PatternManager:                                                 │
│    - _load_patterns()      # Load YAML                          │
│    - _interpolate(pattern) # Replace {FLOAT} with actual regex  │
│    - get_regex(name)       # Return compiled regex              │
│    - match(name, text)     # Run match, return groups           │
│    - search(name, text)    # Run search, return groups          │
│    - findall(name, text)   # Find all matches                   │
│                                                                  │
│  Global: patterns = PatternManager()                             │
└─────────────────────────────────────────────────────────────────┘
                              │
                    ┌─────────┴─────────┐
                    ▼                   ▼
        ┌───────────────────┐  ┌───────────────────┐
        │  programs.yaml    │  │  Python code      │
        │                   │  │                   │
        │  r_free:          │  │  from agent...    │
        │    pattern: ...   │  │    import patterns│
        │    uses: {FLOAT}  │  │                   │
        └───────────────────┘  │  patterns.match() │
                               └───────────────────┘
```

---

## Phase 1: Audit Existing Patterns

### Step 1.1: Find All Regex in Python Code

```bash
grep -rn "re\.\(compile\|search\|match\|findall\)" agent/*.py knowledge/*.py
```

**Expected locations:**
- `agent/session.py` - Filename parsing, log extraction
- `agent/best_files_tracker.py` - Cycle number extraction
- `agent/workflow_state.py` - File pattern matching
- `knowledge/metric_patterns.py` - Metric extraction from logs
- `knowledge/yaml_loader.py` - Pattern compilation

### Step 1.2: Find All Patterns in YAML

```bash
grep -rn "pattern:" knowledge/*.yaml
```

**Expected locations:**
- `programs.yaml` - Log parsing patterns for each program
- `file_categories.yaml` - Filename patterns (mostly globs, some regex)

### Step 1.3: Categorize Patterns

| Category | Examples | Count (est.) |
|----------|----------|--------------|
| **Primitives** | FLOAT, INT, WHITESPACE | 5-10 |
| **Filename** | refine_XXX_YYY.pdb, PHASER.N.pdb | 10-15 |
| **Metrics** | R-free, resolution, map CC | 20-30 |
| **Log Markers** | Error:, Warning:, Traceback | 5-10 |
| **Program Output** | TFZ score, LLG, clash count | 15-20 |

---

## Phase 2: Create patterns.yaml

### Step 2.1: Define Primitives

```yaml
# knowledge/patterns.yaml

# =============================================================================
# PRIMITIVES - Atomic building blocks
# =============================================================================
primitives:
  # Numbers
  FLOAT:
    regex: "[-+]?(?:\\d*\\.\\d+|\\d+\\.?)"
    description: "Floating point number (1.0, -0.5, 10, 10.)"
    test_matches: ["1.0", "-0.5", "10", "0.123", "10."]
    test_non_matches: ["abc", "1.2.3"]
  
  FLOAT_POSITIVE:
    regex: "\\d*\\.?\\d+"
    description: "Positive float only"
    test_matches: ["1.0", "0.5", "10"]
    test_non_matches: ["-1.0"]
  
  INT:
    regex: "\\d+"
    description: "Positive integer"
    test_matches: ["1", "123", "0"]
    test_non_matches: ["-1", "1.5"]
  
  INT_PADDED:
    regex: "\\d{3}"
    description: "Zero-padded 3-digit integer (001, 012, 123)"
    test_matches: ["001", "012", "123"]
    test_non_matches: ["1", "12", "1234"]
  
  # Whitespace
  WS:
    regex: "\\s+"
    description: "One or more whitespace"
  
  OPT_WS:
    regex: "\\s*"
    description: "Optional whitespace (zero or more)"
  
  # Common tokens
  EQUALS:
    regex: "[=:]"
    description: "Equals sign or colon"
  
  SEPARATOR:
    regex: "[=:\\s]+"
    description: "Separator (equals, colon, or whitespace)"
```

### Step 2.2: Define System Patterns

```yaml
# =============================================================================
# SYSTEM PATTERNS - Used by Python engine for file/directory handling
# =============================================================================
system:
  # Cycle detection from directory names
  run_directory:
    regex: "run_{INT}"
    compiled: "run_(\\d+)"
    description: "Extract cycle number from run directory"
    examples: ["run_001", "run_12", "run_1"]
    groups:
      1: cycle_number
  
  # Refinement output files
  refine_output:
    regex: "refine_{INT_PADDED}_{INT_PADDED}\\.pdb"
    compiled: "refine_(\\d{3})_(\\d{3})\\.pdb"
    description: "Standard phenix.refine output naming"
    examples: ["refine_001_001.pdb", "refine_002_003.pdb"]
    groups:
      1: macro_cycle
      2: micro_cycle
  
  # Real-space refine output
  rsr_output:
    regex: "(?:real_space_refined|rsr).*{INT}"
    description: "Real-space refinement output"
    examples: ["model_real_space_refined_001.pdb", "rsr_cycle_5.pdb"]
  
  # Phaser output
  phaser_output:
    regex: "PHASER\\.{INT}\\.pdb"
    compiled: "PHASER\\.(\\d+)\\.pdb"
    description: "Phaser molecular replacement output"
    examples: ["PHASER.1.pdb", "PHASER.2.pdb"]
    groups:
      1: solution_number
  
  # Error detection
  python_traceback:
    regex: "Traceback \\(most recent call last\\):"
    description: "Python exception traceback marker"
  
  phenix_error:
    regex: "^(?:Error|ERROR|Sorry):{OPT_WS}"
    description: "PHENIX error message prefix"
```

### Step 2.3: Define Metric Patterns

```yaml
# =============================================================================
# METRIC PATTERNS - For extracting values from program logs
# =============================================================================
metrics:
  # R-factors (X-ray)
  r_free:
    regex: "R-?free{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "R-free value from refinement"
    examples: 
      - input: "R-free = 0.2534"
        group_1: "0.2534"
      - input: "Rfree: 0.25"
        group_1: "0.25"
  
  r_work:
    regex: "R-?work{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "R-work value from refinement"
  
  r_factor_pair:
    regex: "R-?work{OPT_WS}{EQUALS}{OPT_WS}({FLOAT}){OPT_WS}R-?free{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "R-work and R-free on same line"
    groups:
      1: r_work
      2: r_free
  
  # Resolution
  resolution:
    regex: "(?:resolution|Resolution){OPT_WS}{EQUALS}{OPT_WS}({FLOAT}){OPT_WS}(?:Å|A|Angstrom)?"
    description: "Resolution in Angstroms"
  
  resolution_range:
    regex: "({FLOAT}){OPT_WS}-{OPT_WS}({FLOAT}){OPT_WS}(?:Å|A)"
    description: "Resolution range (high - low)"
    groups:
      1: high_resolution
      2: low_resolution
  
  # Cryo-EM metrics
  map_cc:
    regex: "(?:map_model_cc|Map-model CC|CC_mask){OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "Map-model correlation coefficient"
  
  fsc_resolution:
    regex: "FSC{OPT_WS}(?:0\\.143|=0\\.5)?{OPT_WS}{EQUALS}?{OPT_WS}({FLOAT}){OPT_WS}(?:Å|A)"
    description: "FSC-based resolution estimate"
  
  # Geometry metrics
  clashscore:
    regex: "(?:Clashscore|clashscore){OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "MolProbity clashscore"
  
  ramachandran_outliers:
    regex: "Ramachandran{OPT_WS}outliers{OPT_WS}{EQUALS}{OPT_WS}({FLOAT}){OPT_WS}%?"
    description: "Percentage of Ramachandran outliers"
  
  # Molecular replacement
  tfz_score:
    regex: "TFZ{OPT_WS}(?:score)?{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "Phaser translation function Z-score"
  
  llg:
    regex: "LLG{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})"
    description: "Log-likelihood gain from Phaser"
```

---

## Phase 3: Create PatternManager

### Step 3.1: Create agent/pattern_manager.py

```python
"""
Pattern Manager for PHENIX AI Agent.

Provides centralized regex pattern management with:
- Primitive interpolation ({FLOAT} -> actual regex)
- Named pattern lookup
- Compiled regex caching
- Match/search/findall helpers

Usage:
    from agent.pattern_manager import patterns
    
    # Match a filename
    match = patterns.match("system.refine_output", "refine_001_002.pdb")
    if match:
        cycle = match.group(1)
    
    # Extract a metric from log text
    r_free = patterns.extract("metrics.r_free", log_text)
"""

from __future__ import absolute_import, division, print_function

import os
import re
from functools import lru_cache

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
        """Initialize the pattern manager."""
        self._primitives = {}
        self._patterns = {}
        self._compiled = {}
        
        if patterns_file is None:
            patterns_file = os.path.join(
                os.path.dirname(os.path.dirname(__file__)),
                "knowledge", "patterns.yaml"
            )
        
        self._load_patterns(patterns_file)
    
    def _load_patterns(self, filepath):
        """Load and parse patterns.yaml."""
        if not yaml:
            raise ImportError("PyYAML required for PatternManager")
        
        if not os.path.exists(filepath):
            # Graceful degradation if file doesn't exist yet
            return
        
        with open(filepath, 'r') as f:
            data = yaml.safe_load(f)
        
        # Load primitives first
        self._primitives = {}
        for name, defn in data.get("primitives", {}).items():
            self._primitives[name] = defn.get("regex", "")
        
        # Load pattern sections
        for section in ["system", "metrics", "filenames"]:
            for name, defn in data.get(section, {}).items():
                full_name = f"{section}.{name}"
                self._patterns[full_name] = defn
    
    def _interpolate(self, pattern):
        """
        Replace {PRIMITIVE} placeholders with actual regex.
        
        Example: "R-free{OPT_WS}={OPT_WS}({FLOAT})"
              -> "R-free\\s*=\\s*([-+]?(?:\\d*\\.\\d+|\\d+))"
        """
        result = pattern
        for name, regex in self._primitives.items():
            placeholder = "{" + name + "}"
            result = result.replace(placeholder, regex)
        return result
    
    @lru_cache(maxsize=256)
    def get_compiled(self, name, flags=0):
        """
        Get compiled regex for a named pattern.
        
        Args:
            name: Pattern name (e.g., "metrics.r_free")
            flags: Optional re flags (re.IGNORECASE, etc.)
        
        Returns:
            Compiled regex pattern
        """
        if name not in self._patterns:
            raise KeyError(f"Unknown pattern: {name}")
        
        defn = self._patterns[name]
        
        # Use pre-compiled if available, otherwise interpolate
        if "compiled" in defn:
            regex_str = defn["compiled"]
        else:
            regex_str = self._interpolate(defn.get("regex", ""))
        
        return re.compile(regex_str, flags)
    
    def get_regex_string(self, name):
        """Get the raw (interpolated) regex string for a pattern."""
        if name not in self._patterns:
            raise KeyError(f"Unknown pattern: {name}")
        
        defn = self._patterns[name]
        if "compiled" in defn:
            return defn["compiled"]
        return self._interpolate(defn.get("regex", ""))
    
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
        return pattern.match(text)
    
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
        return pattern.search(text)
    
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
        return pattern.findall(text)
    
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
            except IndexError:
                return default
        return default
    
    def extract_float(self, name, text, group=1, default=None, flags=0):
        """Extract a float value from text."""
        value = self.extract(name, text, group, None, flags)
        if value is not None:
            try:
                return float(value)
            except ValueError:
                pass
        return default
    
    def extract_int(self, name, text, group=1, default=None, flags=0):
        """Extract an integer value from text."""
        value = self.extract(name, text, group, None, flags)
        if value is not None:
            try:
                return int(value)
            except ValueError:
                pass
        return default
    
    def list_patterns(self, section=None):
        """List available pattern names."""
        if section:
            prefix = f"{section}."
            return [n for n in self._patterns.keys() if n.startswith(prefix)]
        return list(self._patterns.keys())
    
    def get_description(self, name):
        """Get the description for a pattern."""
        if name in self._patterns:
            return self._patterns[name].get("description", "")
        return ""
    
    def validate_pattern(self, name):
        """
        Validate a pattern against its test cases.
        
        Returns:
            (passed, failed) tuple of lists
        """
        if name not in self._patterns:
            raise KeyError(f"Unknown pattern: {name}")
        
        defn = self._patterns[name]
        pattern = self.get_compiled(name)
        
        passed = []
        failed = []
        
        # Test matches
        for test in defn.get("test_matches", []):
            if pattern.search(test):
                passed.append(("match", test))
            else:
                failed.append(("match", test, "should match but didn't"))
        
        # Test non-matches
        for test in defn.get("test_non_matches", []):
            if pattern.search(test):
                failed.append(("non_match", test, "should not match but did"))
            else:
                passed.append(("non_match", test))
        
        # Test examples with expected groups
        for example in defn.get("examples", []):
            if isinstance(example, dict):
                text = example.get("input", "")
                match = pattern.search(text)
                if not match:
                    failed.append(("example", text, "no match"))
                else:
                    for key, expected in example.items():
                        if key.startswith("group_"):
                            group_num = int(key.split("_")[1])
                            actual = match.group(group_num)
                            if actual != expected:
                                failed.append(("example", text, 
                                    f"group {group_num}: expected {expected}, got {actual}"))
                            else:
                                passed.append(("example", text, f"group {group_num}"))
        
        return passed, failed


# Global instance - lazy loaded
_patterns = None

def get_patterns():
    """Get the global PatternManager instance."""
    global _patterns
    if _patterns is None:
        _patterns = PatternManager()
    return _patterns

# Convenience - can do `from agent.pattern_manager import patterns`
class _PatternsProxy:
    """Proxy that lazily loads the PatternManager."""
    def __getattr__(self, name):
        return getattr(get_patterns(), name)

patterns = _PatternsProxy()
```

### Step 3.2: Create Tests

```python
# tests/test_pattern_manager.py

def test_primitive_interpolation():
    """Test that primitives are correctly interpolated."""
    pm = PatternManager()
    
    # Test FLOAT interpolation
    result = pm._interpolate("{FLOAT}")
    assert "\\d" in result
    
    # Test multiple primitives
    result = pm._interpolate("value{OPT_WS}={OPT_WS}{FLOAT}")
    assert "\\s*" in result


def test_pattern_matching():
    """Test pattern matching works."""
    pm = PatternManager()
    
    # Test system patterns
    match = pm.match("system.refine_output", "refine_001_002.pdb")
    assert match is not None
    assert match.group(1) == "001"
    assert match.group(2) == "002"


def test_metric_extraction():
    """Test metric extraction from log text."""
    pm = PatternManager()
    
    log = "Final R-free = 0.2534"
    value = pm.extract_float("metrics.r_free", log)
    assert abs(value - 0.2534) < 0.0001


def test_pattern_validation():
    """Test pattern self-validation."""
    pm = PatternManager()
    
    passed, failed = pm.validate_pattern("primitives.FLOAT")
    assert len(failed) == 0, f"FLOAT validation failed: {failed}"
```

---

## Phase 4: Migrate Existing Code

### Step 4.1: Update metric_patterns.py

**Current:**
```python
def _get_refine_patterns():
    return {
        "r_free": {
            "patterns": [
                r"R-free\s*[=:]\s*([0-9.]+)",
                r"Final R-free\s*=\s*([0-9.]+)",
            ],
            ...
        }
    }
```

**New:**
```python
def _get_refine_patterns():
    from agent.pattern_manager import patterns
    
    return {
        "r_free": {
            "patterns": [
                patterns.get_regex_string("metrics.r_free"),
                patterns.get_regex_string("metrics.r_free_final"),
            ],
            ...
        }
    }
```

### Step 4.2: Update best_files_tracker.py

**Current:**
```python
def _extract_cycle_number(self, path):
    match = re.search(r"_(\d+)\.pdb$", path)
    if match:
        return int(match.group(1))
```

**New:**
```python
def _extract_cycle_number(self, path):
    from agent.pattern_manager import patterns
    
    # Try refine output pattern first
    cycle = patterns.extract_int("system.refine_output", path, group=1)
    if cycle is not None:
        return cycle
    
    # Fallback to generic cycle pattern
    return patterns.extract_int("system.generic_cycle", path, default=0)
```

### Step 4.3: Update programs.yaml

**Current:**
```yaml
phenix.refine:
  log_patterns:
    r_free:
      patterns:
        - 'R-free\s*[=:]\s*([0-9.]+)'
        - 'Final R-free\s*=\s*([0-9.]+)'
```

**New:**
```yaml
phenix.refine:
  log_patterns:
    r_free:
      pattern_refs:
        - "metrics.r_free"
        - "metrics.r_free_final"
      # Or keep inline but use interpolation:
      patterns:
        - 'R-free{OPT_WS}{EQUALS}{OPT_WS}({FLOAT})'
```

---

## Phase 5: Testing and Validation

### Step 5.1: Create Validation Script

```python
# scripts/validate_patterns.py

"""Validate all patterns in patterns.yaml against their test cases."""

from agent.pattern_manager import PatternManager

def main():
    pm = PatternManager()
    
    all_passed = True
    for name in pm.list_patterns():
        passed, failed = pm.validate_pattern(name)
        if failed:
            print(f"❌ {name}:")
            for f in failed:
                print(f"   {f}")
            all_passed = False
        else:
            print(f"✓ {name} ({len(passed)} tests)")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    exit(main())
```

### Step 5.2: Add to Test Suite

Add `test_pattern_manager.py` to the test runner.

---

## Migration Checklist

### Files to Create
- [ ] `knowledge/patterns.yaml` - Pattern definitions
- [ ] `agent/pattern_manager.py` - PatternManager class
- [ ] `tests/test_pattern_manager.py` - Unit tests

### Files to Modify
- [ ] `knowledge/metric_patterns.py` - Use PatternManager for metric patterns
- [ ] `agent/best_files_tracker.py` - Use PatternManager for filename parsing
- [ ] `agent/session.py` - Use PatternManager for log parsing
- [ ] `agent/workflow_state.py` - Use PatternManager where regex used
- [ ] `knowledge/yaml_loader.py` - Add interpolation support when loading programs.yaml

### Patterns to Migrate (Priority Order)

**High Priority (used everywhere):**
1. FLOAT, INT primitives
2. R-free, R-work extraction
3. Resolution extraction
4. Filename cycle extraction

**Medium Priority:**
5. Map CC extraction
6. Clashscore extraction
7. TFZ score extraction
8. Error/warning detection

**Low Priority (can migrate later):**
9. File category patterns (mostly globs)
10. Program-specific specialized patterns

---

## Success Criteria

1. **All existing tests pass** after migration
2. **Pattern validation script** passes with 0 failures
3. **No inline regex** in core Python files (agent/*.py)
4. **Consistent metrics** - same pattern extracts same value everywhere
5. **Self-documenting** - new developers can read patterns.yaml and understand

---

## Timeline Estimate

| Phase | Effort | Dependencies |
|-------|--------|--------------|
| Phase 1: Audit | 1 hour | None |
| Phase 2: patterns.yaml | 2 hours | Phase 1 |
| Phase 3: PatternManager | 2 hours | Phase 2 |
| Phase 4: Migration | 3-4 hours | Phase 3 |
| Phase 5: Testing | 1 hour | Phase 4 |

**Total: ~8-10 hours**

---

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Regex behavior changes | Tests fail | Run full test suite after each migration |
| Performance regression | Slower parsing | Cache compiled patterns (already in design) |
| Missing edge cases | Metrics not extracted | Keep original patterns as fallback initially |
| YAML parsing issues | Crash on load | Graceful degradation if patterns.yaml missing |
