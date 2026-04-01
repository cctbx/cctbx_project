# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The Computational Crystallography Toolbox (cctbx) is a large-scale crystallographic computation library developed as the open source component of the Phenix project. The codebase contains 1M+ lines of code organized as a village of independent modules (toolboxes), each providing specific functionality. The architecture combines ISO C++ classes with Python bindings for computational efficiency and flexibility.

## Key Modules (Toolboxes)

The project is organized into specialized modules, each ending in "tbx" (toolbox):

- **cctbx**: Core crystallographic computations (symmetry, unit cells, Miller indices, structure factors)
- **scitbx**: Scientific toolbox (linear algebra, optimization, array operations, mathematical functions)
- **iotbx**: Input/output toolbox (PDB, CIF, MTZ, reflection files)
- **mmtbx**: Macromolecular toolbox (refinement, validation, model building)
- **smtbx**: Small molecule toolbox (includes olex2.refine engine)
- **libtbx**: Library toolbox (build system, utilities, testing infrastructure)
- **rstbx**: Reciprocal space toolbox
- **gltbx**: OpenGL toolbox (graphics)
- **fable**: Fortran to C++ converter
- **boost_adaptbx**: Boost Python bindings
- Additional specialized modules: simtbx, spotfinder, xfel, wxtbx, etc.

## Architecture

- **Hybrid C++/Python**: Core algorithms in C++ for performance, Python bindings for usability
- **Modular design**: Each module is independent with minimal cross-dependencies
- **SCons build system**: Uses SConscript files for compilation
- **libtbx infrastructure**: Provides module configuration (libtbx_refresh.py) and build utilities

## Common Development Commands

source build/setpaths.sh to make correct python and all the commands available
libtbx.scons -j 11 to recompile

### Running Tests

libtbx.run_tests_parallel needs to be executed in empty folder.

```bash
# Run all tests for a specific module (fastest, recommended)
libtbx.run_tests_parallel module=<module_name>

# Examples:
libtbx.run_tests_parallel module=cctbx
libtbx.run_tests_parallel module=mmtbx

# Run tests in a specific directory
libtbx.run_tests_parallel directory=<path>

# Run tests with multiple processes
libtbx.run_tests_parallel module=cctbx nproc=4

# Run a single test file
libtbx.python <module>/regression/tst_example.py
```

### Code Quality Tools

```bash
# Check for common code issues (tabs, trailing whitespace, missing imports, bare excepts)
libtbx.find_clutter

# Fix tabs and trailing whitespace
libtbx.clean_clutter

# Find unused imports
libtbx.find_unused_imports_crude

# Add future division imports
libtbx.add_from_future_import_division
```

### Building

The project uses SCons for building C++ extensions. After modifying C++ code or adding new files:

```bash
# Rebuild (specific commands depend on your build configuration)
libtbx.scons
```

For development installations, the build is typically managed through the bootstrap.py script during initial setup.

## Testing Guidelines

### Test Organization

- Tests located in `<module>/regression/` directories
- Test files named `tst_<name>.py`
- Each module has `run_tests.py` that lists all tests to execute
- **CRITICAL**: New tests MUST be added to `<module>/run_tests.py` or they won't run

### Test Structure

```python
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

def exercise():
  """Brief description of what this test validates."""
  # Test code here
  result = some_function()
  assert approx_equal(result, expected_value)

if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
```

### Test Requirements

- Keep runtime under 30 seconds (60 seconds absolute max)
- Tests should be focused on one functionality
- Use `libtbx.test_utils.approx_equal()` for numerical comparisons
- Prefer generated data over stored files
- Output actual values in failed assertions: `assert a>b, "%f > %f failed" % (a,b)`
- All new functionality requires unit tests
- All bug fixes require regression tests

## Coding Standards

### Python Style

- Use 2-space indentation (not PEP8's 4 spaces)
- Always include `from __future__ import absolute_import, division, print_function`
- Use `//` for integer division, `/` for floating-point division
- Use `print(text, file=log)` not bare `print(text)`
- Never use `import *` (except in `__init__.py` for C++ extensions)
- Never use `isinstance()` - rely on duck typing
- Nested imports preferred over global imports to reduce overhead
- No tabs or trailing whitespace
- No bare `except:` statements - use `except Exception:` at minimum

### Code Organization

- Avoid code duplication - search for existing functionality first
- Use self-explanatory names over short cryptic names
- Functions should generally fit on one page
- Add docstrings and comments for clarity
- Update documentation when code changes
- Use central constants (e.g., `math.pi`, `scitbx.constants.two_pi_sq`)

### Before Committing

1. Run `libtbx.find_clutter` and fix all issues
2. Run `libtbx.clean_clutter` to fix whitespace
3. Run relevant tests with `libtbx.run_tests_parallel`
4. Ensure all tests pass - **fixing broken tests is highest priority**

## Development Workflow

### Before Writing Code

1. **Ask questions** to fully understand the goal - clarify ambiguous tasks
2. **Read existing code** before modifying - understand function purpose, callers, and tests
3. **Check for existing functionality** - search the codebase to avoid duplication
4. **Plan the implementation** as markdown covering:
   - Problem to be solved
   - Overall approach
   - Step-by-step implementation
   - Unit tests for new/modified code
   - Affected existing tests
5. **Before adding dependencies**: Discuss with team - new dependencies are strongly discouraged

### While Writing Code

1. **Break into clearly defined steps** - present each step for review
2. **Run tests after each step** - fix failures before proceeding
3. **Self-review before presenting**:
   - Handle edge cases (None, empty lists, missing keys, invalid types)
   - Cross-platform compatibility (macOS, Linux, Windows)
   - No scoping issues (variable defined where referenced)
   - Match surrounding code style
   - Would you want to fix anything if asked for improvements?
4. **Final audit**:
   - Parse-check all changed files
   - Run test suites that exercise changed code
   - Check line-length violations
   - Verify imports resolve

### After Writing Code

1. **Update documentation** for user-visible or architectural changes
2. **Run code quality tools**: `libtbx.find_clutter`, `libtbx.clean_clutter`
3. **Run tests**: `libtbx.run_tests_parallel module=<affected_module>`
4. **Ensure all tests pass** - fixing broken tests is highest priority

## Common Patterns

### Using Models and Data in Tests

```python
# Use existing models
from mmtbx.regression import model_1yjp
# or
pdb_file = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/pdb/1yjp_h.pdb",
  test=os.path.isfile)

# Generate random structures
from cctbx.development import random_structure

# Embed small test structures as strings
pdb_str = """
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 10.00           C
"""
pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
```

### Checking for Optional Modules

```python
import libtbx.load_env
if not libtbx.env.has_module("phenix"):
  print("phenix tree missing, skipping test")
  return
```

### Common Imports

```python
# Crystallographic core
from cctbx import crystal, miller, xray, sgtbx

# PHENIX utilities
from libtbx.utils import Sorry, null_out
from libtbx import group_args, adopt_init_args

# File I/O
from iotbx import pdb, mtz, cif
from iotbx.data_manager import DataManager

# Testing
from libtbx.test_utils import approx_equal, show_diff, Exception_expected
```

### Crystallographic Conventions

- **Resolution**: Always in Ångströms, stored as `d_min` (not `resolution`)
- **R-factors**: Always fractions (0.25), stored as `r_free` / `r_work` (not `rfree`)
- **Angles**: Degrees
- **B-factors**: Ų
- **Space groups**: Use `sgtbx.space_group_info` objects, not strings
  - Compare with `sg1.group() == sg2.group()`, not string comparison
- **Arrays**: Use `scitbx.array_family.flex`, not NumPy
  - `flex.sum(array)`, not `array.sum()`
  - Boolean selection: `array.select(flags)`, not `array[flags]`

## File Locations

- **Tests**: `<module>/regression/tst_*.py`
- **Test registry**: `<module>/run_tests.py`
- **Module config**: `<module>/libtbx_refresh.py`
- **Build scripts**: `<module>/SConscript`
- **C++ code**: Typically in `<module>/boost_python/` or module root
- **Python code**: Throughout module directories

## LLM Programming Guidelines

**For AI Agent development** (libtbx/langchain subsystem), see comprehensive LLM-specific guidelines:
- `libtbx/langchain/docs/prompts/CCTBX_LLM_PROGRAMMING_GUIDELINES.md` - General cctbx/PHENIX patterns
- `libtbx/langchain/docs/prompts/AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md` - Agent-specific patterns

These documents provide detailed guidance on development process, error handling, common pitfalls, and agent-specific patterns.

### AI Agent Subsystem Specifics

The AI Agent code (libtbx/langchain) has unique requirements:

- **Parameter verification**: Check `knowledge/programs.yaml` strategy_flags before using any command parameter - the sanitizer strips unlisted parameters
- **Import fallbacks**: All imports in `agent/` and `knowledge/` need try/except for both installed PHENIX and standalone contexts:
  ```python
  try:
    from libtbx.langchain.agent.structure_model import StructureModel
  except ImportError:
    from agent.structure_model import StructureModel
  ```
- **Logging patterns**:
  - `agent/` and `knowledge/`: Use `logging.getLogger(__name__)`
  - `programs/ai_agent.py`: Use `self.vlog` (quiet/normal/verbose)
  - Test files: Use `print()`
- **Session state persistence**: Every state field must appear in:
  1. Written to `session.data["key"]`
  2. Read from `session.data["key"]` on resume
  3. Included in `create_initial_state()`
  4. Registered in `agent/contract.py` with default value
- **Prompt templates**: Store in `agent/prompts/` as plain text files, not hardcoded in Python
- **Client-server awareness**: Changes to `agent/` and `knowledge/` are server-side (immediate), changes to `programs/ai_agent.py` are client-side (require user update)

## Critical Patterns & Pitfalls

### None-Safety with dict.get()

**Common bug**: `.get("key", "")` returns `""` when key is missing, but `None` when key exists with value `None`. This crashes on `.lower()`, `.upper()`, `.strip()`:

```python
# WRONG — crashes when value is None
name = params.get("model_name", "").lower()

# RIGHT — handles both missing and None
name = (params.get("model_name") or "").lower()
```

### Platform-Independent Paths

Never use hardcoded `/` or `\` in file paths:

```python
# WRONG
path = data_dir + "/" + "model.pdb"

# RIGHT
path = os.path.join(data_dir, "model.pdb")
```

### R-free: Fraction vs Percentage

R-free is always a fraction (0.25), never a percentage (25):

```python
# WRONG
if r_free < 25:  # This is a bug!

# RIGHT
if r_free < 0.25:
```

### Floating Point Comparisons

Never use `==` for floating point values:

```python
# WRONG
assert result.r_free == 0.25

# RIGHT
from libtbx.test_utils import approx_equal
assert approx_equal(result.r_free, 0.25, eps=1e-4)
```

### PHIL Parameters

Don't invent PHENIX parameters - verify against actual PHIL definitions:

```python
# For AI Agent: check knowledge/programs.yaml strategy_flags
# For PHENIX programs: verify against program PHIL definitions
```

### Error Handling

Use `libtbx.utils.Sorry` for user-facing errors, standard exceptions for programming errors:

```python
from libtbx.utils import Sorry

# User error — clear message
if not os.path.isfile(pdb_file):
  raise Sorry("Model file not found: %s" % pdb_file)

# Programming error
assert space_group is not None, \
  "space_group must be set before refinement"
```

Never use bare `except:` - always catch `Exception` or more specific types. Document exception swallowing:

```python
# Non-critical cleanup — safe to skip on failure
try:
  os.remove(temp_file)
except Exception:
  pass
```

### Subprocess Management

Use `libtbx.easy_run` instead of `subprocess` or `os.system()`:

```python
from libtbx import easy_run

result = easy_run.fully_buffered(
  command="phenix.refine model.pdb data.mtz"
)
if result.return_code != 0:
  raise Sorry("Refinement failed: %s" % result.stderr_lines[-1])
```

## Pre-Commit Checklist

Before committing code changes:

- [ ] All changed files parse cleanly (`python -c "import py_compile; py_compile.compile('<file>')"`)
- [ ] All relevant test suites pass (`libtbx.run_tests_parallel module=<name>`)
- [ ] New tests added to `<module>/run_tests.py`
- [ ] Run `libtbx.find_clutter` and fix all issues
- [ ] Run `libtbx.clean_clutter` to fix whitespace
- [ ] No bare `except:` blocks - use `except Exception:` with explanatory comment
- [ ] No hardcoded path separators - use `os.path.join()`
- [ ] None-safety: `(x.get("key") or "")` not `x.get("key", "")`
- [ ] R-free comparisons use fractions (0.25), not percentages (25)
- [ ] Floating point comparisons use `approx_equal()`, not `==`
- [ ] New `from __future__ import absolute_import, division, print_function`
- [ ] Use `//` for integer division, `/` for floating-point division
- [ ] Documentation updated for user-visible or architectural changes

**For AI Agent code additionally:**
- [ ] New strategy flags added to `knowledge/programs.yaml` before use
- [ ] Import fallbacks present for all `libtbx.langchain` imports
- [ ] New session state fields in `session.data`, `create_initial_state()`, and `contract.py`
- [ ] New YAML files validated with `agent/yaml_tools.py validate`
- [ ] Prompt templates stored in `agent/prompts/`, not hardcoded

## Important Notes

- This is an **installed build** located at `/Users/oleg/Documents/phenix/installs/cif_svn/modules/cctbx_project`
- The project is part of the larger Phenix crystallography package
- Main development branch: `master`
- Community support: cctbxbb@phenix-online.org mailing list
- For development builds, use `bootstrap.py` with `--use-conda` flag for better OS compatibility
- **Fixing broken tests is the highest priority** - broken tests block other developers from committing
