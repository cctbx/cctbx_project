# Implementation Plan: Robust Program Configuration System

## Overview

This document describes the implementation of a system to make program configuration more robust and less error-prone. The system has three main components:

1. **Sensible Defaults** - Automatic inference of `input_priorities` from input names
2. **Validation with Warnings** - Load-time checks that warn (or crash in testing) on suspicious configurations
3. **Dry Run File Selection** - A mode to preview file selection without running programs

---

## Phase 1: Sensible Defaults for Input Priorities

### 1.1 Create New Module: `agent/input_defaults.py`

```python
"""
Default input priorities for common input slot names.

When a program doesn't define explicit input_priorities for an input,
these defaults are used based on the input name and extensions.

This prevents subtle bugs where the wrong file is selected because
the program author forgot to specify input_priorities.
"""

# Mapping of input name patterns to default priorities
# Keys are checked in order - first match wins
# Patterns are matched case-insensitively against input_name

INPUT_NAME_DEFAULTS = {
    # Model/protein inputs - prefer refined models, exclude ligands and templates
    "model": {
        "categories": ["model"],
        "prefer_subcategories": ["refined", "with_ligand", "phaser_output", "docked"],
        "exclude_categories": ["ligand", "search_model", "intermediate", "ligand_fit_output"],
    },
    "protein": {
        "categories": ["model"],
        "prefer_subcategories": ["refined", "with_ligand"],
        "exclude_categories": ["ligand", "search_model", "intermediate", "ligand_fit_output"],
    },
    "pdb_file": {
        "categories": ["model"],
        "prefer_subcategories": ["refined", "with_ligand", "phaser_output"],
        "exclude_categories": ["ligand", "search_model", "intermediate"],
    },
    
    # MTZ/reflection data inputs
    "mtz": {
        "categories": ["mtz"],
        "prefer_subcategories": ["refined_mtz"],
        "exclude_categories": ["phased_mtz"],
    },
    "hkl_file": {
        "categories": ["mtz"],
        "prefer_subcategories": ["refined_mtz"],
        "exclude_categories": ["phased_mtz"],
    },
    "data": {
        "categories": ["mtz"],
        "prefer_subcategories": ["refined_mtz"],
        "exclude_categories": [],
    },
    "xray_data": {
        "categories": ["mtz"],
        "prefer_subcategories": ["refined_mtz"],
        "exclude_categories": ["phased_mtz"],
    },
    
    # Map inputs (cryo-EM)
    "map": {
        "categories": ["map", "full_map"],
        "prefer_subcategories": ["sharpened", "full_map"],
        "exclude_categories": ["half_map"],
    },
    "full_map": {
        "categories": ["full_map", "map"],
        "prefer_subcategories": [],
        "exclude_categories": ["half_map"],
    },
    "half_map": {
        "categories": ["half_map"],
        "prefer_subcategories": [],
        "exclude_categories": ["full_map", "map"],
    },
    
    # Ligand inputs
    "ligand": {
        "categories": ["ligand_fit_output", "ligand"],
        "prefer_subcategories": ["ligand_fit_output"],
        "exclude_categories": ["model", "search_model", "refined"],
    },
    "ligand_file": {
        "categories": ["ligand_fit_output", "ligand"],
        "prefer_subcategories": ["ligand_fit_output"],
        "exclude_categories": ["model", "search_model"],
    },
    
    # Search model / template inputs (for MR, docking)
    "search_model": {
        "categories": ["search_model"],
        "prefer_subcategories": ["processed_predicted", "predicted"],
        "exclude_categories": ["model", "ligand", "refined"],
    },
    "template": {
        "categories": ["search_model"],
        "prefer_subcategories": ["processed_predicted", "predicted"],
        "exclude_categories": ["model", "ligand"],
    },
    
    # Sequence inputs
    "sequence": {
        "categories": ["sequence"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
    "seq_file": {
        "categories": ["sequence"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
}

# Extension-based fallback defaults (when name doesn't match)
# Used when input name is non-standard but extensions give hints
EXTENSION_DEFAULTS = {
    ".pdb": {
        "categories": ["model"],
        "prefer_subcategories": ["refined"],
        "exclude_categories": ["ligand", "intermediate"],
    },
    ".cif": {
        "categories": ["model"],
        "prefer_subcategories": ["refined"],
        "exclude_categories": ["ligand", "intermediate"],
    },
    ".mtz": {
        "categories": ["mtz"],
        "prefer_subcategories": ["refined_mtz"],
        "exclude_categories": [],
    },
    ".mrc": {
        "categories": ["map", "full_map"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
    ".ccp4": {
        "categories": ["map", "full_map"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
    ".fa": {
        "categories": ["sequence"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
    ".fasta": {
        "categories": ["sequence"],
        "prefer_subcategories": [],
        "exclude_categories": [],
    },
}


def get_default_input_priorities(input_name, extensions=None):
    """
    Get default input_priorities for an input based on its name and extensions.
    
    Args:
        input_name: Name of the input slot (e.g., "model", "protein", "mtz")
        extensions: List of file extensions this input accepts (e.g., [".pdb", ".cif"])
    
    Returns:
        dict with keys: categories, prefer_subcategories, exclude_categories
        Returns None if no reasonable default can be inferred.
    """
    input_lower = input_name.lower()
    
    # Try exact match first
    if input_lower in INPUT_NAME_DEFAULTS:
        return INPUT_NAME_DEFAULTS[input_lower].copy()
    
    # Try partial match (e.g., "protein_file" matches "protein")
    for pattern, defaults in INPUT_NAME_DEFAULTS.items():
        if pattern in input_lower or input_lower in pattern:
            return defaults.copy()
    
    # Fall back to extension-based inference
    if extensions:
        for ext in extensions:
            ext_lower = ext.lower()
            if ext_lower in EXTENSION_DEFAULTS:
                return EXTENSION_DEFAULTS[ext_lower].copy()
    
    # No default available
    return None


def can_infer_defaults(input_name, extensions=None):
    """
    Check if defaults can be inferred for this input.
    
    Used by the validator to decide whether to warn about missing input_priorities.
    """
    return get_default_input_priorities(input_name, extensions) is not None
```

### 1.2 Integrate Defaults into Command Builder

**File**: `agent/command_builder.py`

**Changes to `_find_file_for_slot()` method**:

```python
# At the start of PRIORITY 3 section, add:

# Get input_priorities from registry, or infer defaults
priorities = self._registry.get_input_priorities(program, input_name)

if not priorities or not priorities.get("categories"):
    # No explicit priorities - try to infer defaults
    from libtbx.langchain.agent.input_defaults import get_default_input_priorities
    extensions = input_def.get("extensions", [])
    defaults = get_default_input_priorities(input_name, extensions)
    
    if defaults:
        self._log(context, "BUILD: Using default priorities for '%s': categories=%s" % (
            input_name, defaults.get("categories", [])))
        priorities = defaults
```

### 1.3 Unit Tests for Defaults

**File**: `tests/test_input_defaults.py`

```python
"""Tests for input_defaults module."""

import unittest
from libtbx.langchain.agent.input_defaults import (
    get_default_input_priorities,
    can_infer_defaults,
)


class TestInputDefaults(unittest.TestCase):
    
    def test_model_input(self):
        """Model input should prefer refined, exclude ligands."""
        defaults = get_default_input_priorities("model")
        self.assertIn("model", defaults["categories"])
        self.assertIn("refined", defaults["prefer_subcategories"])
        self.assertIn("ligand", defaults["exclude_categories"])
    
    def test_protein_input(self):
        """Protein input should work like model."""
        defaults = get_default_input_priorities("protein")
        self.assertIn("model", defaults["categories"])
        self.assertIn("ligand", defaults["exclude_categories"])
    
    def test_mtz_input(self):
        """MTZ input should prefer refined_mtz."""
        defaults = get_default_input_priorities("mtz")
        self.assertIn("mtz", defaults["categories"])
        self.assertIn("refined_mtz", defaults["prefer_subcategories"])
    
    def test_ligand_input(self):
        """Ligand input should prefer ligand_fit_output, exclude model."""
        defaults = get_default_input_priorities("ligand")
        self.assertIn("ligand_fit_output", defaults["categories"])
        self.assertIn("model", defaults["exclude_categories"])
    
    def test_extension_fallback(self):
        """Unknown name with .pdb extension should get model defaults."""
        defaults = get_default_input_priorities("coords", extensions=[".pdb"])
        self.assertIn("model", defaults["categories"])
    
    def test_no_match_returns_none(self):
        """Completely unknown input should return None."""
        defaults = get_default_input_priorities("foobar", extensions=[".xyz"])
        self.assertIsNone(defaults)
    
    def test_partial_match(self):
        """Partial name match should work (e.g., 'model_file' matches 'model')."""
        defaults = get_default_input_priorities("model_file")
        self.assertIn("model", defaults["categories"])
    
    def test_can_infer(self):
        """can_infer_defaults should return True for known patterns."""
        self.assertTrue(can_infer_defaults("model"))
        self.assertTrue(can_infer_defaults("mtz"))
        self.assertTrue(can_infer_defaults("unknown", [".pdb"]))
        self.assertFalse(can_infer_defaults("foobar", [".xyz"]))


if __name__ == "__main__":
    unittest.main()
```

---

## Phase 2: Validation with Warnings

### 2.1 Create Validator Module: `agent/program_validator.py`

```python
"""
Validator for program definitions in programs.yaml.

Checks for common configuration mistakes and missing fields.
Can be run at load time (with warnings) or in strict mode (crashes on issues).
"""

from __future__ import absolute_import, division, print_function

import os
import sys


class ValidationIssue:
    """A single validation issue."""
    
    LEVELS = ["info", "warning", "error"]
    
    def __init__(self, level, program, message, field=None):
        self.level = level
        self.program = program
        self.message = message
        self.field = field
    
    def __str__(self):
        prefix = f"[{self.level.upper()}] {self.program}"
        if self.field:
            prefix += f".{self.field}"
        return f"{prefix}: {self.message}"


class ProgramValidator:
    """
    Validates program definitions for completeness and correctness.
    
    Usage:
        validator = ProgramValidator(strict=False)
        issues = validator.validate_all(programs_dict)
        validator.report(issues)
    
    In strict mode, warnings are treated as errors and will raise an exception.
    """
    
    def __init__(self, strict=False):
        """
        Initialize validator.
        
        Args:
            strict: If True, warnings are treated as errors (for testing)
        """
        self.strict = strict
        self._issues = []
    
    def validate_all(self, programs):
        """
        Validate all program definitions.
        
        Args:
            programs: Dict of program_name -> program_definition
        
        Returns:
            List of ValidationIssue objects
        """
        self._issues = []
        
        for program_name, program_def in programs.items():
            self._validate_program(program_name, program_def)
        
        return self._issues
    
    def _validate_program(self, name, definition):
        """Validate a single program definition."""
        
        # Check required fields
        if not definition.get("description"):
            self._add_issue("warning", name, "Missing description", "description")
        
        if not definition.get("command"):
            self._add_issue("error", name, "Missing command template", "command")
        
        if not definition.get("experiment_types"):
            self._add_issue("warning", name, 
                "Missing experiment_types - program won't appear in any workflow",
                "experiment_types")
        
        # Check inputs
        inputs = definition.get("inputs", {})
        required_inputs = inputs.get("required", {})
        optional_inputs = inputs.get("optional", {})
        all_inputs = {**required_inputs, **optional_inputs}
        
        if not required_inputs and not optional_inputs:
            self._add_issue("info", name, "No inputs defined - program takes no files")
        
        # Check each input for proper configuration
        input_priorities = definition.get("input_priorities", {})
        
        for input_name, input_def in all_inputs.items():
            self._validate_input(name, input_name, input_def, input_priorities)
        
        # Check outputs
        outputs = definition.get("outputs", {})
        if not outputs.get("files"):
            self._add_issue("info", name,
                "No output patterns defined - file tracking will rely on directory scanning",
                "outputs")
        
        # Check command template has all required slots
        command = definition.get("command", "")
        for input_name in required_inputs:
            slot = "{%s}" % input_name
            if slot not in command:
                self._add_issue("warning", name,
                    f"Required input '{input_name}' not in command template",
                    "command")
    
    def _validate_input(self, program, input_name, input_def, input_priorities):
        """Validate a single input definition."""
        from libtbx.langchain.agent.input_defaults import can_infer_defaults
        
        extensions = input_def.get("extensions", [])
        
        # Check if input has explicit priorities
        has_explicit_priorities = input_name in input_priorities
        can_use_defaults = can_infer_defaults(input_name, extensions)
        
        if not has_explicit_priorities and not can_use_defaults:
            # This is problematic - no priorities and can't infer
            self._add_issue("warning", program,
                f"Input '{input_name}' has no input_priorities and name/extensions "
                f"don't match known patterns - file selection may be unreliable",
                f"inputs.{input_name}")
        
        # Check for potentially ambiguous inputs
        if not has_explicit_priorities and extensions:
            pdb_like = any(ext in [".pdb", ".cif"] for ext in extensions)
            if pdb_like and "ligand" not in input_name.lower():
                # PDB input without explicit priorities - could pick ligand by mistake
                if not input_def.get("exclude_patterns"):
                    self._add_issue("info", program,
                        f"Input '{input_name}' accepts PDB files but has no "
                        f"exclude_patterns - consider excluding ligand patterns",
                        f"inputs.{input_name}")
    
    def _add_issue(self, level, program, message, field=None):
        """Add a validation issue."""
        issue = ValidationIssue(level, program, message, field)
        self._issues.append(issue)
    
    def report(self, issues=None, out=None):
        """
        Print validation report.
        
        Args:
            issues: List of issues (uses self._issues if None)
            out: Output stream (default: sys.stderr)
        
        Raises:
            ValueError in strict mode if there are warnings or errors
        """
        if issues is None:
            issues = self._issues
        
        if out is None:
            out = sys.stderr
        
        if not issues:
            return
        
        # Group by level
        errors = [i for i in issues if i.level == "error"]
        warnings = [i for i in issues if i.level == "warning"]
        infos = [i for i in issues if i.level == "info"]
        
        # Print issues
        for issue in errors + warnings + infos:
            print(str(issue), file=out)
        
        # Summary
        print(f"\nValidation: {len(errors)} errors, {len(warnings)} warnings, "
              f"{len(infos)} info", file=out)
        
        # In strict mode, fail on warnings or errors
        if self.strict and (errors or warnings):
            raise ValueError(
                f"Program validation failed in strict mode: "
                f"{len(errors)} errors, {len(warnings)} warnings"
            )
    
    def has_errors(self):
        """Check if there are any error-level issues."""
        return any(i.level == "error" for i in self._issues)
    
    def has_warnings(self):
        """Check if there are any warning-level issues."""
        return any(i.level == "warning" for i in self._issues)


def validate_programs(programs, strict=False, out=None):
    """
    Convenience function to validate programs and report issues.
    
    Args:
        programs: Dict of program definitions
        strict: If True, raise on warnings/errors
        out: Output stream for report
    
    Returns:
        List of ValidationIssue objects
    """
    validator = ProgramValidator(strict=strict)
    issues = validator.validate_all(programs)
    if issues:
        validator.report(issues, out=out)
    return issues
```

### 2.2 Integrate Validation into YAML Loader

**File**: `knowledge/yaml_loader.py`

Add validation call after loading programs:

```python
def load_programs(validate=True, strict=False):
    """
    Load program definitions from YAML.
    
    Args:
        validate: If True, run validation checks
        strict: If True (and validate=True), fail on warnings
    
    Returns:
        Dict of program definitions
    """
    programs = _load_yaml_file("programs.yaml")
    
    if validate:
        from libtbx.langchain.agent.program_validator import validate_programs
        validate_programs(programs, strict=strict)
    
    return programs
```

### 2.3 Environment Variable for Strict Mode

Add support for `PHENIX_AI_STRICT_VALIDATION=1` environment variable:

```python
# In yaml_loader.py
import os

def _should_use_strict_validation():
    """Check if strict validation is enabled via environment."""
    return os.environ.get("PHENIX_AI_STRICT_VALIDATION", "").lower() in ("1", "true", "yes")
```

### 2.4 Test Configuration for Strict Mode

**File**: `tests/conftest.py` or test runner setup

```python
import os

# Enable strict validation for all tests
os.environ["PHENIX_AI_STRICT_VALIDATION"] = "1"
```

### 2.5 Unit Tests for Validator

**File**: `tests/test_program_validator.py`

```python
"""Tests for program_validator module."""

import unittest
from libtbx.langchain.agent.program_validator import (
    ProgramValidator,
    ValidationIssue,
    validate_programs,
)


class TestProgramValidator(unittest.TestCase):
    
    def test_valid_program(self):
        """Well-configured program should have no issues."""
        programs = {
            "phenix.test": {
                "description": "Test program",
                "command": "phenix.test {model}",
                "experiment_types": ["xray"],
                "inputs": {
                    "required": {
                        "model": {"extensions": [".pdb"]}
                    }
                },
                "input_priorities": {
                    "model": {"categories": ["model"]}
                },
                "outputs": {
                    "files": [{"pattern": "*_out.pdb"}]
                }
            }
        }
        
        validator = ProgramValidator(strict=False)
        issues = validator.validate_all(programs)
        
        errors = [i for i in issues if i.level == "error"]
        warnings = [i for i in issues if i.level == "warning"]
        self.assertEqual(len(errors), 0)
        self.assertEqual(len(warnings), 0)
    
    def test_missing_command(self):
        """Missing command should be an error."""
        programs = {
            "phenix.test": {
                "description": "Test program",
                "experiment_types": ["xray"],
            }
        }
        
        validator = ProgramValidator(strict=False)
        issues = validator.validate_all(programs)
        
        errors = [i for i in issues if i.level == "error"]
        self.assertEqual(len(errors), 1)
        self.assertIn("command", errors[0].message.lower())
    
    def test_unknown_input_without_priorities(self):
        """Unknown input name without priorities should warn."""
        programs = {
            "phenix.test": {
                "description": "Test program",
                "command": "phenix.test {foobar}",
                "experiment_types": ["xray"],
                "inputs": {
                    "required": {
                        "foobar": {"extensions": [".xyz"]}
                    }
                }
            }
        }
        
        validator = ProgramValidator(strict=False)
        issues = validator.validate_all(programs)
        
        warnings = [i for i in issues if i.level == "warning"]
        self.assertTrue(any("foobar" in w.message for w in warnings))
    
    def test_known_input_uses_defaults(self):
        """Known input name without priorities should NOT warn (uses defaults)."""
        programs = {
            "phenix.test": {
                "description": "Test program",
                "command": "phenix.test {model}",
                "experiment_types": ["xray"],
                "inputs": {
                    "required": {
                        "model": {"extensions": [".pdb"]}
                    }
                }
                # Note: no input_priorities - should use defaults
            }
        }
        
        validator = ProgramValidator(strict=False)
        issues = validator.validate_all(programs)
        
        # Should not warn about 'model' because defaults can be inferred
        warnings = [i for i in issues if i.level == "warning"]
        model_warnings = [w for w in warnings if "model" in w.message]
        self.assertEqual(len(model_warnings), 0)
    
    def test_strict_mode_raises(self):
        """Strict mode should raise on warnings."""
        programs = {
            "phenix.test": {
                "description": "Test program",
                # Missing command - will be an error
            }
        }
        
        validator = ProgramValidator(strict=True)
        issues = validator.validate_all(programs)
        
        with self.assertRaises(ValueError):
            validator.report(issues)


if __name__ == "__main__":
    unittest.main()
```

---

## Phase 3: Dry Run File Selection Mode

### 3.1 Add PHIL Parameter

**File**: `programs/ai_agent.py`

Add to `master_params`:

```python
  dry_run_file_selection = False
    .type = bool
    .short_caption = Preview file selection without running programs
    .help = When True, shows which files would be selected for each program \
            without actually executing any commands. Useful for debugging \
            file selection issues. The agent will run through its decision \
            process but stop before execution.
```

### 3.2 Implement Dry Run File Selection Mode

**File**: `programs/ai_agent.py`

In `_run_single_cycle()`, add check before command execution:

```python
def _run_single_cycle(self, cycle, session, session_start_time):
    """Run a single cycle. Returns True if should break out of loop."""
    self.vlog.cycle_header(cycle)
    
    session.start_cycle(cycle)
    self.params.ai_analysis.log_as_simple_string = None
    
    command, decision_info = self._get_command_for_cycle(cycle, session)
    
    # ... existing validation code ...
    
    # DRY RUN FILE SELECTION MODE
    if getattr(self.params.ai_analysis, 'dry_run_file_selection', False):
        self._display_file_selection_preview(command, decision_info, cycle, session)
        session.record_result(cycle, "DRY_RUN_FILE_SELECTION: Preview only, no execution")
        return False  # Continue to next cycle for more previews
    
    # ... rest of existing code ...
```

Add new method:

```python
def _display_file_selection_preview(self, command, decision_info, cycle, session):
    """Display detailed file selection information without executing."""
    self.vlog.normal("\n" + "=" * 60)
    self.vlog.normal("DRY RUN FILE SELECTION - Cycle %d" % cycle)
    self.vlog.normal("=" * 60)
    
    program = decision_info.get('program', 'unknown')
    self.vlog.normal("\nProgram: %s" % program)
    self.vlog.normal("Reasoning: %s" % decision_info.get('reasoning', 'N/A'))
    
    # Show the command that would be run
    self.vlog.normal("\nCommand that would be run:")
    self.vlog.normal("  %s" % command)
    
    # Parse command to show individual file selections
    self.vlog.normal("\nFile selections:")
    # Extract file paths from command
    import shlex
    try:
        parts = shlex.split(command)
        for part in parts[1:]:  # Skip program name
            if os.path.exists(part):
                basename = os.path.basename(part)
                # Try to determine category
                category = self._get_file_category(part, session)
                self.vlog.normal("  - %s [%s]" % (basename, category or "uncategorized"))
            elif "=" in part:
                # Strategy flag
                self.vlog.normal("  - %s [flag]" % part)
    except Exception as e:
        self.vlog.normal("  (Could not parse command: %s)" % e)
    
    # Show what files were available
    active_files = session.get_available_files()
    self.vlog.normal("\nAvailable files (%d total):" % len(active_files))
    for f in active_files[:10]:
        category = self._get_file_category(f, session)
        self.vlog.normal("  - %s [%s]" % (os.path.basename(f), category or "uncategorized"))
    if len(active_files) > 10:
        self.vlog.normal("  ... and %d more" % (len(active_files) - 10))
    
    # Show best files
    best_files = session.get_best_files_dict()
    if best_files:
        self.vlog.normal("\nBest files:")
        for category, path in best_files.items():
            if path:
                self.vlog.normal("  - %s: %s" % (category, os.path.basename(path)))
    
    self.vlog.normal("\n" + "=" * 60)
    self.vlog.normal("(No command executed - dry_run_file_selection=True)")
    self.vlog.normal("=" * 60 + "\n")

def _get_file_category(self, filepath, session):
    """Get the category of a file from session's categorized files."""
    basename = os.path.basename(filepath)
    # This would need access to categorized_files from the last agent call
    # For now, return None - we can enhance this later
    return None
```

### 3.3 Limit Cycles in Dry Run Mode

To prevent infinite loops in dry run mode, add a max cycles check:

```python
# In iterate_agent(), when dry_run_file_selection is True
if getattr(self.params.ai_analysis, 'dry_run_file_selection', False):
    max_preview_cycles = 3  # Only preview first 3 cycles
    self.vlog.normal("DRY RUN FILE SELECTION MODE - previewing up to %d cycles" % max_preview_cycles)
```

---

## Phase 4: Documentation

### 4.1 Update `knowledge/programs.yaml` Header

Add comprehensive documentation at the top of the file:

```yaml
# =============================================================================
# PHENIX AI AGENT - PROGRAM DEFINITIONS
# =============================================================================
#
# This file defines how the AI Agent interacts with PHENIX programs.
#
# ADDING A NEW PROGRAM - CHECKLIST
# --------------------------------
#
# Required fields:
#   [ ] description    - Human-readable description of what the program does
#   [ ] command        - Command template with {slot} placeholders for inputs
#   [ ] experiment_types - List: [xray], [cryoem], or [xray, cryoem]
#   [ ] inputs.required - Dict of required input slots with extensions
#
# Recommended fields:
#   [ ] input_priorities - Category-based file selection for each input
#                          (If omitted, defaults are inferred from input name)
#   [ ] outputs.files   - Patterns for output files (helps with tracking)
#   [ ] category        - Program category (refinement, building, etc.)
#
# Optional fields:
#   [ ] inputs.optional - Dict of optional input slots
#   [ ] strategy_flags  - Command-line flags the agent can set
#   [ ] fixes           - Auto-fill rules (resolution, output_prefix, etc.)
#   [ ] invariants      - Validation rules
#   [ ] log_parsing     - Regex patterns for extracting metrics from logs
#
# INPUT PRIORITIES
# ----------------
# For each input, you can specify:
#   categories:           - File categories to search (e.g., [model, search_model])
#   prefer_subcategories: - Preferred subcategories in order (e.g., [refined, phaser_output])
#   exclude_categories:   - Categories to never use (e.g., [ligand, intermediate])
#   fallback_categories:  - Try these if primary categories are empty
#
# If input_priorities is not specified, defaults are inferred from the input name:
#   - "model", "protein", "pdb" → prefers refined models, excludes ligands
#   - "mtz", "data"             → prefers refined_mtz
#   - "ligand"                  → prefers ligand_fit_output
#   - "search_model", "template"→ prefers processed_predicted
#
# COMMON MISTAKES TO AVOID
# ------------------------
# 1. Forgetting to exclude ligand files from protein/model inputs
#    → Solution: Add exclude_categories: [ligand, ligand_fit_output]
#
# 2. Forgetting to exclude search_model from refinement inputs
#    → Solution: Add exclude_categories: [search_model]
#
# 3. Not specifying prefer_subcategories
#    → Result: May select oldest file instead of best/most recent
#
# 4. PDB inputs without exclude_patterns for ligand patterns
#    → Solution: Add exclude_patterns: [ligand_fit, lig.pdb]
#
# VALIDATION
# ----------
# Run validation with: phenix.ai_agent --validate-programs
# Enable strict mode (fails on warnings) with: PHENIX_AI_STRICT_VALIDATION=1
#
# =============================================================================
```

### 4.2 Create `docs/ADDING_PROGRAMS.md`

```markdown
# Adding Programs to PHENIX AI Agent

This guide explains how to add new PHENIX programs to the AI Agent's repertoire.

## Quick Start

Add your program to `knowledge/programs.yaml`:

```yaml
phenix.my_program:
  description: "Brief description of what this program does"
  category: utility  # refinement, building, analysis, utility, etc.
  experiment_types: [xray, cryoem]  # or just [xray] or [cryoem]
  
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""  # No flag prefix (positional argument)
      data:
        extensions: [.mtz]
        flag: ""
  
  command: "phenix.my_program {model} {data}"
  
  outputs:
    files:
      - pattern: "*_output.pdb"
        type: model
```

## Detailed Configuration

### Input Priorities (Recommended)

For precise control over file selection, add `input_priorities`:

```yaml
  input_priorities:
    model:
      categories: [model]
      prefer_subcategories: [refined, with_ligand]
      exclude_categories: [ligand, search_model]
    data:
      categories: [mtz]
      prefer_subcategories: [refined_mtz]
```

### Strategy Flags

Allow the agent to customize program behavior:

```yaml
  strategy_flags:
    resolution:
      flag: "high_resolution={value}"
      type: float
    cycles:
      flag: "main.number_of_macro_cycles={value}"
      type: int
      default: 3
```

### Output Patterns

Help the agent track output files:

```yaml
  outputs:
    files:
      - pattern: "*_refined.pdb"
        type: model
      - pattern: "*_refined.mtz"
        type: mtz
    metrics:
      - r_free
      - r_work
```

## Testing Your Configuration

1. **Validate syntax**:
   ```bash
   phenix.ai_agent --validate-programs
   ```

2. **Preview file selection**:
   ```bash
   phenix.ai_agent dry_run_file_selection=True ...
   ```

3. **Test with verbose output**:
   ```bash
   phenix.ai_agent verbosity=verbose ...
   ```

## Common Patterns

### Refinement Program
```yaml
phenix.my_refine:
  category: refinement
  experiment_types: [xray]
  inputs:
    required:
      model:
        extensions: [.pdb]
      mtz:
        extensions: [.mtz]
  input_priorities:
    model:
      categories: [model]
      prefer_subcategories: [refined]
      exclude_categories: [ligand, search_model]
    mtz:
      categories: [mtz]
      prefer_subcategories: [refined_mtz]
```

### Utility Program (Combining Files)
```yaml
phenix.combine:
  category: utility
  inputs:
    required:
      protein:
        extensions: [.pdb]
        exclude_patterns: [ligand, lig.]
      ligand:
        extensions: [.pdb]
        prefer_patterns: [ligand_fit]
  input_priorities:
    protein:
      categories: [model]
      exclude_categories: [ligand, ligand_fit_output]
    ligand:
      categories: [ligand_fit_output, ligand]
      exclude_categories: [model]
```
```

---

## Implementation Order

### Step 1: Create `agent/input_defaults.py`
- Implement `INPUT_NAME_DEFAULTS` mapping
- Implement `EXTENSION_DEFAULTS` mapping  
- Implement `get_default_input_priorities()` function
- Implement `can_infer_defaults()` function
- Add unit tests

### Step 2: Integrate defaults into `command_builder.py`
- Modify `_find_file_for_slot()` to use defaults when `input_priorities` is missing
- Add logging when defaults are used
- Test with existing programs

### Step 3: Create `agent/program_validator.py`
- Implement `ValidationIssue` class
- Implement `ProgramValidator` class
- Add validation checks for common issues
- Add unit tests

### Step 4: Integrate validation into `yaml_loader.py`
- Add `validate` and `strict` parameters to `load_programs()`
- Add environment variable support for strict mode
- Configure tests to use strict mode

### Step 5: Add dry_run_file_selection mode
- Add PHIL parameter
- Implement `_display_file_selection_preview()` method
- Add cycle limit for dry run mode
- Test with sample workflows

### Step 6: Documentation
- Update `programs.yaml` header with comprehensive guide
- Create `docs/ADDING_PROGRAMS.md`
- Add examples for common program types

---

## Verification

After implementation, verify:

1. **Existing programs still work**: Run existing test scenarios
2. **Defaults are applied correctly**: Check logs for "Using default priorities" messages
3. **Validation catches issues**: Test with intentionally broken program definitions
4. **Strict mode works in tests**: Ensure tests fail on validation warnings
5. **Dry run shows useful info**: Test dry_run_file_selection with sample files
6. **Documentation is accurate**: Follow the guide to add a test program

---

## Future Enhancements

1. **Auto-generate input_priorities**: Tool to suggest priorities based on program analysis
2. **Visual file selection debugger**: Web UI showing file flow
3. **Program templates**: Pre-defined templates for common program types
4. **Validation in CI**: Automated validation on programs.yaml changes
