# Guidelines for cctbx/PHENIX Development

You are working on the cctbx/PHENIX codebase — a large
scientific software ecosystem for macromolecular
crystallography and cryo-EM. Follow these guidelines to
produce code that integrates cleanly.

---

## 1. Process

### Before writing code

- Ask enough questions to fully understand the goal. If the
  task is ambiguous, ask for clarification rather than guessing.
- Read existing code before modifying it. Understand what the
  function does, who calls it, and what tests cover it.
- Create a plan as a markdown file. The plan should address:
  - The problem to be solved
  - The overall approach
  - A detailed, step-by-step implementation
  - Unit tests for all new and modified code
  - Which existing tests might be affected
- Show the plan to the programmer and ask for suggestions.
  Recommend that the programmer show the plan to a second LLM
  for additional review. Build the code only after the plan is
  agreed upon.

### While writing code

- Break the implementation into clearly defined steps. At the
  conclusion of each step, present the current state to the
  programmer for review.
- Run all relevant tests after each step. Fix failures before
  proceeding.
- **Self-review before presenting each step.** Before telling
  the programmer a step is done, stop and ask yourself:
  - Did I handle the edge cases? (None values, empty lists,
    missing keys, invalid input types)
  - Does this work on all platforms? (macOS, Linux, Windows)
  - Are there scoping issues? (Variable defined in one
    method but referenced in another)
  - Did I match the style of surrounding code? (Indent
    level, naming conventions, error handling pattern)
  - Would I want to fix anything if the programmer asked
    "any other fixes or additions to this change?"
  This self-review catches the majority of follow-up bugs.
  It is much cheaper to fix issues before presenting than
  to require an additional round-trip.
- After finishing, do a full audit:
  - Parse-check every changed file (`ast.parse()` or
    `python -c "import py_compile; py_compile.compile()"`)
  - Run the test suites that exercise the changed code
  - Check for line-length violations in changed regions
  - Verify imports resolve correctly

### After writing code

- Update documentation for any user-visible or architectural
  changes.
- Update test counts and file inventories where applicable.

---

## 2. Style and Formatting

### Indentation and line width

- 2-space indentation throughout (the standard cctbx
  convention). Some legacy modules may differ — match the
  host file's convention when editing.
- 80-character maximum line width where feasible. Long
  strings, URLs, and regex patterns may exceed this.

### Naming

- Descriptive function and class names:
  `find_peaks_holes()`, not `fph()`.
- Private helpers prefixed with underscore: `_safe_float()`.
- Constants in UPPER_SNAKE_CASE:
  `MAX_RESIDUES_FOR_QUICK_SCAN = 500`.
- Boolean variables read as assertions: `has_model`,
  `is_twinned`, `wants_prediction`.
- Use cctbx naming conventions for crystallographic
  concepts: `d_min` (not `resolution`), `r_free` /
  `r_work` (not `rfree`), `space_group_info` (not `sg`).

### Cleanup

- Remove all trailing whitespace.
- Remove all unused imports. If the libtbx linter flags an
  import that IS used (e.g., re-exported or used only in a
  string context), mark it so the linter stays quiet:

  ```python
  from libtbx.utils import Sorry  # noqa (re-raised)
  assert Sorry is not None  # used in except clause below
  ```

  Or use a tuple reference:

  ```python
  (Sorry,)  # Suppress linter false-positive
  ```

### Docstrings

- Every public function and class gets a docstring.
- Include `Args:`, `Returns:`, and `Never raises.` (if
  applicable) sections.
- For methods that intentionally swallow exceptions, document
  this as a "Never raises" contract and explain why:

  ```python
  def update_from_model(self, ...):
    """Update state from model inspection.

    Never raises — partial updates are acceptable because
    downstream code handles missing fields gracefully.
    """
  ```

### Type hints

Type hints are optional in cctbx — usage varies across
modules.

- Match the host file's convention. If the file you're
  editing uses type hints, add them to new functions. If
  it doesn't, don't introduce them.
- When using hints, import from `typing` for broad Python
  compatibility: `Optional`, `List`, `Dict`, `Tuple`.
- Type hints on function signatures are useful. Type hints
  on every local variable are noise — skip those.

### Code presentation

When presenting code to the programmer:

- Never use placeholder comments like
  `# ... rest of code stays the same ...` inside a code
  block. The programmer will paste your output directly.
- For small changes: show the full modified function.
- For large files: show the exact old text and new text
  so the change is unambiguous.
- If you must show a partial file, clearly mark the
  boundaries with line numbers.

---

## 3. cctbx/PHENIX-Specific Patterns

### Don't invent parameters

LLMs frequently hallucinate PHENIX command-line arguments
or PHIL parameters that don't exist.

- If you are generating or modifying a PHENIX command
  string, verify parameters against the program's PHIL
  definition or its documentation. Do not guess parameter
  names from context.
- If you are calling a cctbx/PHENIX Python API, check the
  actual function signature in the source. Don't assume
  keyword arguments from the function name.
- When in doubt, ask: "Is `twin_law` a valid parameter
  for `phenix.refine`?" rather than assuming.

### PHIL parameters

PHENIX programs use the PHIL (Python Hierarchical
Interchange Language) system for configuration. When
writing or modifying PHIL definitions:

```python
my_param = 10
  .type = int
  .short_caption = My parameter
  .help = Description of what this does. \
          Continuation lines use backslash.
```

- `.type` values: `int`, `float`, `str`, `bool`, `path`,
  `choice`, `strings`, `ints`, `floats`
- For `choice` parameters, mark the default with `*`:
  `mode = fast *thorough exhaustive`
- PHIL help strings use `\` for continuation — do not use
  triple-quoted strings or markdown.
- Test that PHIL parsing works:
  `master_params.fetch(source).extract()`

### Import conventions

cctbx modules are often imported with specific patterns.
Common conventions:

```python
# Miller arrays, crystal symmetry, etc.
from cctbx import crystal, miller, xray

# PHENIX-specific utilities
from libtbx.utils import Sorry, null_out
from libtbx import group_args, adopt_init_args

# File I/O
from iotbx import pdb, mtz, cif
from iotbx.data_manager import DataManager
```

When writing code that may run in multiple contexts
(installed PHENIX vs. development tree), use try/except
for optional imports:

```python
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  ProgramTemplate = None  # Not available outside PHENIX
```

### easy_run and subprocesses

PHENIX uses `libtbx.easy_run` for subprocess management.
Prefer it over raw `subprocess` calls:

```python
from libtbx import easy_run

result = easy_run.fully_buffered(
  command="phenix.refine model.pdb data.mtz"
)
if result.return_code != 0:
  raise Sorry("Refinement failed: %s"
    % result.stderr_lines[-1])
```

Never use `os.system()`. When `easy_run` is not
available, use `subprocess.run()` with explicit
`check=True` or manual return-code checking.

### None-safety on dict.get()

**This is a common bug pattern in Python generally and
in cctbx code that processes JSON or PHIL parameters.**

`.get("key", "")` returns the default `""` when the key
is MISSING, but returns `None` when the key EXISTS with
value `None`. This crashes on `.lower()`, `.upper()`,
`.strip()`:

```python
# WRONG — crashes when value is None
name = params.get("model_name", "").lower()

# RIGHT — handles both missing and None
name = (params.get("model_name") or "").lower()
```

### Raising errors

Use `libtbx.utils.Sorry` for user-facing errors (bad
input, missing files, invalid parameters). Use standard
Python exceptions for programming errors:

```python
from libtbx.utils import Sorry

# User error — clear message, no traceback
if not os.path.isfile(pdb_file):
  raise Sorry(
    "Model file not found: %s" % pdb_file)

# Programming error — should never happen
assert space_group is not None, \
  "space_group must be set before refinement"
```

### Test structure

cctbx uses a fail-fast testing style with plain functions:

```python
from libtbx.test_utils import (
  approx_equal, show_diff, Exception_expected,
)

def exercise_basic():
  result = compute(input_data)
  assert result.r_free < 0.3
  assert approx_equal(result.bonds_rmsd, 0.01,
                       eps=0.005)

def exercise_edge_cases():
  # Test with empty input
  result = compute(None)
  assert result is None

def run():
  exercise_basic()
  exercise_edge_cases()
  print("OK")

if __name__ == "__main__":
  run()
```

Conventions:
- Test file names: `tst_module_name.py`
- Test functions: `exercise_*` or `test_*`
- Entry point: `run()` function
- Success indicator: `print("OK")` at the end
- Use `libtbx.test_utils.approx_equal` for floating
  point comparisons, not `==`
- Use `libtbx.test_utils.show_diff` for string
  comparisons that might fail (produces readable diffs)

### Serialization

When writing classes that persist to JSON or pickle:

- `to_dict()` / `from_dict()` must survive a round-trip
- `from_dict()` must be tolerant of missing keys (the
  data may be from an older version):

  ```python
  @classmethod
  def from_dict(cls, d):
    obj = cls()
    # Use .get() with defaults, not d["key"]
    obj.resolution = d.get("resolution")
    obj.space_group = d.get("space_group", "P1")
    return obj
  ```

- Test the round-trip explicitly:

  ```python
  def exercise_serialization():
    original = MyClass(resolution=2.0)
    d = original.to_dict()
    restored = MyClass.from_dict(d)
    assert approx_equal(
      restored.resolution, 2.0)
  ```

---

## 4. Common Pitfalls

### Pitfall: `.get("key", "")` with None values

See "None-safety on dict.get()" above. Use `or ""`
pattern for any `.lower()`, `.upper()`, `.strip()`.

### Pitfall: hardcoded path separators

PHENIX runs on macOS, Linux, and Windows. Never use
hardcoded `/` or `\` in file paths:

```python
# WRONG
path = data_dir + "/" + "model.pdb"

# RIGHT
path = os.path.join(data_dir, "model.pdb")
```

Also use `os.path.normpath()` when comparing paths and
`os.path.relpath()` with a try/except (it raises
`ValueError` across Windows drive letters).

### Pitfall: R-free as percentage vs fraction

R-free is always stored and compared as a fraction
(0.25), never a percentage (25). Resolution is always in
Ångströms. Angles are in degrees. B-factors are in Å².
If you see a threshold like `r_free < 25`, that's a bug
— it should be `r_free < 0.25`.

Variable naming: use `d_min` for resolution limits (the
crystallographic convention), `r_free` / `r_work` for
R-factors, `map_cc` for map correlation coefficients.

### Pitfall: floating point comparison

Never use `==` for floating point values from
crystallographic calculations:

```python
# WRONG
assert result.r_free == 0.25

# RIGHT
from libtbx.test_utils import approx_equal
assert approx_equal(result.r_free, 0.25, eps=1e-4)
```

### Pitfall: large file I/O

PDB files can be megabytes (large complexes, NMR
ensembles). MTZ files can be hundreds of megabytes
(high-resolution, anomalous data). When reading files
for inspection (not full parsing), limit the read:

```python
def _is_small_molecule(path, max_bytes=32768):
  with open(path, "rb") as f:
    header = f.read(max_bytes)
  # ... inspect header only ...
```

For full parsing, use the appropriate iotbx reader
(`iotbx.pdb.input`, `iotbx.mtz.object`) which handles
memory efficiently.

### Pitfall: forgetting `flex` array semantics

cctbx uses `scitbx.array_family.flex` arrays, not NumPy.
Common gotchas:

- `flex.double()` creates an empty array, not a scalar
- Indexing returns Python objects, not flex scalars
- `flex.sum()` is a function, not a method:
  `flex.sum(array)`, not `array.sum()`
- Boolean selection: `array.select(flags)`, not
  `array[flags]`

### Pitfall: space group conventions

cctbx uses `sgtbx.space_group_info` objects, not strings.
Converting between representations:

```python
from cctbx import sgtbx

# String to object
sgi = sgtbx.space_group_info("P 21 21 21")

# Object to string (various formats)
hall = sgi.type().hall_symbol()
hm = str(sgi)  # Hermann-Mauguin
number = sgi.type().number()
```

Never compare space groups as strings — different
notations represent the same group. Use:

```python
sg1.group() == sg2.group()
```

---

## 5. Error Handling

### "Never raises" functions

Functions that interact with external data (file I/O,
log parsing, model inspection) should follow the "never
raises" pattern when they are non-critical:

```python
def get_resolution_from_log(log_text):
  """Extract resolution from program log.

  Never raises — returns None if parsing fails.
  """
  try:
    return _parse_resolution(log_text)
  except Exception:
    return None
```

The caller handles `None` gracefully. This keeps
pipelines running when one data source produces
unexpected output.

### No bare excepts

Never use bare `except:` — it catches `SystemExit` and
`KeyboardInterrupt`, which prevents clean shutdown and
Ctrl-C. Always catch a specific type:

```python
# WRONG — catches SystemExit, KeyboardInterrupt
try:
  do_something()
except:
  pass

# RIGHT — lets SystemExit/KeyboardInterrupt propagate
try:
  do_something()
except Exception:
  pass
```

### Exception swallowing requires a comment

Even with `except Exception`, never write a bare `pass`
without documenting why the exception is safe to discard.
Every swallowed exception must have either a comment or a
log statement:

```python
# Non-critical cleanup — safe to skip on failure
try:
  os.remove(temp_file)
except Exception:
  pass
```

Or with logging (preferred when debugging may be needed):

```python
try:
  result = compute_map_coefficients(...)
except Exception as e:
  logger.debug("Map coefficient computation "
    "failed: %s" % e)
  result = None
```

---

## 6. Checklist Before Presenting Code

Before telling the programmer a step is done:

- [ ] All changed files parse cleanly
- [ ] All relevant test suites pass
- [ ] No lines > 80 chars in changed regions (where
      feasible)
- [ ] No trailing whitespace
- [ ] No unused imports
- [ ] No bare `except:` — always `except Exception:` or
      more specific
- [ ] Every `except Exception: pass` has a comment or
      log explaining why the exception is safe to discard
- [ ] `(x.get("key") or "")` pattern for any `.lower()`,
      `.upper()`, `.strip()` on dict values that could be
      None
- [ ] `os.path.join()` for all file paths, never string
      concatenation with `/`
- [ ] `approx_equal()` for floating point comparisons in
      tests, never `==`
- [ ] Serialization round-trip test for any new
      `to_dict()` / `from_dict()` pair
- [ ] Documentation updated for user-visible or
      architectural changes
