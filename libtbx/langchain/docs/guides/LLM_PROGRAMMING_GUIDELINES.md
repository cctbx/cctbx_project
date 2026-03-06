# Guidelines for PHENIX AI Agent Development

You are working on the PHENIX AI Agent codebase — an
LLM-driven automation system for macromolecular structure
determination. Follow these guidelines to produce code
that integrates cleanly.

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
    missing keys, resumed sessions, end-of-plan)
  - Does this work in all contexts? (PHENIX environment,
    standalone, GUI mode off, no plan exists)
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
  - Parse-check every changed file (`ast.parse()`)
  - Run the test suites that exercise the changed code
  - Check for 80-character line violations in changed regions
  - Verify imports work in both PHENIX and standalone contexts

### After writing code

- Update documentation (OVERVIEW.md, ARCHITECTURE.md,
  CHANGELOG.md) for any user-visible or architectural changes.
- Update test counts and file inventories in documentation.
- CHANGELOG entries should follow the established format:

  ```
  ## Version NNN.NN (Short Title)

  ### Summary
  One paragraph describing the change.

  ### New files (N)
  | File | Lines | Purpose |
  |------|-------|---------|
  | `path/to/file.py` | 500 | Description |

  ### Modified files (N)
  | File | Changes |
  |------|---------|
  | `path/to/file.py` | What changed |

  ### Tests
  N new tests. Run with: `python3 tests/run_all_tests.py`
  ```

  See `docs/project/CHANGELOG.md` for examples.

---

## 2. Style and Formatting

### Indentation and line width

- 2-space indentation throughout (matching the existing
  codebase). Some older test files use 4-space — match the
  host file's convention when editing.
- 80-character maximum line width where feasible. Long strings,
  URLs, and regex patterns may exceed this.

### Naming

- Descriptive function and class names:
  `_find_diff_peaks()`, not `_fdp()`.
- Private helpers prefixed with underscore: `_safe_float()`.
- Constants in UPPER_SNAKE_CASE: `PHASE_COMPLETE = "complete"`.
- Boolean variables read as assertions: `has_model`,
  `is_twinned`, `wants_prediction`.

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
  def update_from_validation(self, ...):
    """Update model state from validation results.

    Never raises — partial updates are acceptable because
    the model accumulates over many cycles.
    """
  ```

### Type hints

Type hints are optional in this codebase — some modules use
them (`command_builder.py`, `directive_validator.py`,
`error_analyzer.py`, `graph_state.py`) and some don't.

- Match the host file's convention. If the file you're
  editing uses type hints, add them to new functions. If
  it doesn't, don't introduce them.
- When using hints, import from `typing` for Python 3.7+
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
  (the `str_replace` pattern) so the change is
  unambiguous.
- If you must show a partial file, clearly mark the
  boundaries with line numbers.

---

## 3. PHENIX-Specific Patterns

### Don't invent parameters

LLMs frequently hallucinate PHENIX command-line arguments,
PHIL parameters, or Python API signatures that don't exist.

- If you are generating or modifying a PHENIX command string,
  verify parameters against `programs.yaml`. If a parameter
  is not in the YAML `strategy_flags` for that program, do
  not include it unless the programmer explicitly requested
  it.
- If you are calling a PHENIX Python API, check the actual
  function signature in the source code. Don't guess
  keyword arguments from the function name.
- When in doubt, ask: "Is `rebuild_in_place` a valid
  parameter for `phenix.autobuild`?" rather than assuming.

### Logging conventions

The codebase uses three logging patterns in different
contexts:

- `agent/` and `knowledge/` modules: use
  `logging.getLogger(__name__)` — e.g.,
  `logger.debug("update failed", exc_info=True)`
- `programs/ai_agent.py`: use `self.vlog` at the
  appropriate verbosity — `self.vlog.quiet()` for
  errors, `.normal()` for decisions, `.verbose()` for
  debug detail
- Test files: use `print()` —
  `print("  PASS: test_name")`

Don't mix these patterns — no `print()` in `agent/`
modules, no `logger` in test files.

### File reading limits

When reading PDB, CIF, or MTZ files for content inspection
(not full parsing), limit the read to avoid memory spikes
on large structures:

```python
def _pdb_is_small_molecule(path, max_bytes=32768):
  with open(path, "rb") as f:
    header = f.read(max_bytes)
  # ... inspect header only ...
```

The existing guards in `agent/workflow_state.py` use
`max_bytes=32768` (32 KB). Follow this pattern for any new
content-based checks.

### Import fallbacks

The codebase runs in two contexts: inside PHENIX (where
modules live under `libtbx.langchain`) and standalone (where
they're imported directly). All cross-module imports must have
a fallback:

```python
try:
  from libtbx.langchain.agent.structure_model \
    import StructureModel
except ImportError:
  from agent.structure_model import StructureModel
```

Never use bare `from libtbx.langchain.X import Y` without a
fallback `except ImportError` block.

### None-safety on dict.get()

**This is the single most common bug source in this codebase.**

`.get("key", "")` returns the default `""` when the key is
MISSING, but returns `None` when the key EXISTS with value
`None`. This crashes on `.lower()`, `.upper()`, `.strip()`:

```python
# WRONG — crashes when session_info["explicit_program"] is None
prog = session_info.get("explicit_program", "").lower()

# RIGHT — handles both missing and None
prog = (session_info.get("explicit_program") or "").lower()
```

This applies to any dict from client-supplied data
(`session_info`, `history` entries, `directives`, etc.) where
JSON `null` values are common.

### Test structure

Tests use cctbx-style fail-fast behavior, not `unittest`:

```python
from tests.tst_utils import (
  assert_equal, assert_true, assert_false,
  assert_in, run_tests_with_fail_fast,
)

def test_something():
  result = compute(input)
  assert_equal(result, expected, "helpful message")

def test_another():
  assert_true(condition, "what went wrong")

def run_all_tests():
  run_tests_with_fail_fast()

if __name__ == "__main__":
  run_all_tests()
```

Every test file must:
- Add `sys.path.insert(0, project_root)` before imports
- Have a `run_all_tests()` (or `run_tests()`) function
- Be registered in `tests/run_all_tests.py`
- Use assert helpers from `tests/tst_utils.py`, not bare
  `assert`

### File-exists guards in tests

Tests that open files outside `agent/` and `knowledge/`
(e.g., `programs/ai_agent.py`, `wxGUI2/Programs/AIAgent.py`,
`phenix_ai/remote_agent.py`, `rest/__init__.py`) MUST guard
with `os.path.isfile()` because these directories may not be
present in all deployment contexts:

```python
def test_something_in_ai_agent():
  path = os.path.join(PROJECT_ROOT,
    "programs", "ai_agent.py")
  if not os.path.isfile(path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(path) as f:
    source = f.read()
  assert_in("expected_string", source)
```

Files in `agent/` and `knowledge/` are always present — no
guard needed for those.

### Serialization round-trips

Every class with `to_dict()` / `from_dict()` must survive:

```python
original = MyClass(...)
d = original.to_dict()
restored = MyClass.from_dict(d)
d2 = restored.to_dict()
assert json.dumps(d, sort_keys=True) == \
       json.dumps(d2, sort_keys=True)
```

`from_dict()` should be tolerant of missing keys (the session
may be from an older version of the code). Use `.get()` with
sensible defaults, not direct key access.

### Backward compatibility

The agent has a client-server architecture. Server code may be
newer than client code. Rules:

- Never read `session_info["key"]` — always
  `session_info.get("key", default)`.
- Every new `session_info` field must be registered in
  `agent/contract.py` with a default value.
- Shared code in `agent/` must not import server-only
  dependencies (langchain, openai, anthropic, etc.).

See `docs/guides/BACKWARD_COMPATIBILITY.md` for the full
contract system.

### GUI callbacks

When adding new data to GUI callbacks:

- New fields in existing callbacks (e.g., `agent_cycle`) must
  use `getattr(data, 'new_field', default)` on the GUI side —
  old callbacks won't have the field.
- New callback types (e.g., `agent_gate_transition`) must be
  added to the dispatcher in `AIAgent.py` and ignored silently
  by older GUIs.
- Always wrap callback sends in `try/except Exception: pass`
  with a comment explaining why (callbacks are non-critical
  and must never crash the agent).

---

## 4. Common Pitfalls

### Pitfall: `.get("key", "")` with None values

See "None-safety on dict.get()" above. Use `or ""` pattern.

### Pitfall: search window too small for source scanning

Tests that scan source code (e.g., checking that a function
contains a specific string) must account for function length:

```python
# WRONG — function may be longer than 5000 chars
idx = src.find("def my_function")
assert "expected" in src[idx:idx + 5000]

# RIGHT — search to the next function definition
idx = src.find("def my_function")
next_def = src.find("\ndef ", idx + 1)
assert "expected" in src[idx:next_def]
```

### Pitfall: test failures not propagating

`run_tests_with_fail_fast()` raises `AssertionError` on
failure. If your `run_all_tests()` wraps it in a try/except,
failures will be silently swallowed and the test runner will
report PASSED:

```python
# WRONG — swallows test failures (and uses bare except)
def run_all_tests():
  try:
    run_tests_with_fail_fast()
  except Exception:
    pass  # Test failures silently discarded!

# RIGHT — let failures propagate
def run_all_tests():
  run_tests_with_fail_fast()
```

### Pitfall: stale plan data across thinking levels

Session data persists across runs. If a session was created
with `thinking_level=expert` (which generates a plan), then
resumed with `thinking_level=advanced`, the plan data still
exists in `session.data["plan"]`. Code that reads plan data
must check the current thinking level, not just whether the
data exists.

### Pitfall: `result` vs `status` in history entries

History entries use `"result"` (e.g., `"SUCCESS: OK"`), not
`"status"`. Reading `h.get("status")` returns `None`. Always
check both with a fallback:

```python
status = h.get("result", h.get("status", "?"))
```

### Pitfall: forgetting to persist state on resume

Every piece of state that the agent needs across cycles must
be:

1. Written to `session.data["key"]` after the cycle
2. Read from `session.data["key"]` on resume
3. Included in `create_initial_state()` if passed through
   the graph

If you add a new state field, grep for all three patterns and
verify the field appears in each.

### Pitfall: hardcoded path separators

PHENIX runs on macOS, Linux, and Windows. Never use
hardcoded `/` or `\` in file paths:

```python
# WRONG
path = log_dir + "/" + "session.json"

# RIGHT
path = os.path.join(log_dir, "session.json")
```

Also use `os.path.normpath()` when comparing paths and
`os.path.relpath()` with a try/except (it raises
`ValueError` across Windows drive letters).

### Pitfall: R-free as percentage vs fraction

R-free is always stored and compared as a fraction (0.25),
never a percentage (25). Resolution is always in Ångströms.
Angles are in degrees. If you see a threshold like
`r_free < 25`, that's a bug — it should be `r_free < 0.25`.

Variable naming: use `d_min` for resolution limits (the
crystallographic convention), `r_free` / `r_work` for
R-factors, `map_cc` for map correlation coefficients.

---

## 5. YAML Configuration

The agent's domain knowledge is externalized to YAML files in
`knowledge/`. When adding or modifying YAML:

- Run `python agent/yaml_tools.py validate` after changes
- Ensure the YAML validator handles the file's top-level type
  (dict or list — `expert_knowledge_base_v2.yaml` is a list)
- Skip macOS resource fork files (`._*.yaml`) in any file
  discovery code

---

## 6. Error Handling

### "Never raises" functions

Functions that interact with external data (validation results,
log parsing, file I/O) should follow the "never raises"
pattern when they are non-critical:

```python
def update_from_validation(self, ...):
  """..."""
  try:
    self._update_inner(...)
  except Exception:
    logger.debug("update failed", exc_info=True)
```

This keeps the agent running even when one data source
produces unexpected output. The inner function can raise
freely; the outer wrapper logs and continues.

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
# GUI callback — non-critical, must not crash agent
try:
  call_back(message="agent_plan", ...)
except Exception:
  pass
```

Or with logging (preferred when debugging may be needed):

```python
try:
  self._update_inner(...)
except Exception:
  logger.debug("update failed", exc_info=True)
```

---

## 7. Checklist for Pull Requests

Before submitting code:

- [ ] All changed files parse cleanly (`ast.parse()`)
- [ ] All relevant test suites pass
- [ ] New test functions registered in `run_all_tests.py`
- [ ] No lines > 80 chars in changed regions (where feasible)
- [ ] No trailing whitespace
- [ ] No unused imports
- [ ] No bare `except:` — always `except Exception:` or
      more specific
- [ ] Every `except Exception: pass` has a comment or log
      explaining why the exception is safe to discard
- [ ] Import fallbacks for all `libtbx.langchain` imports
- [ ] `(x.get("key") or "")` pattern for any `.lower()`,
      `.upper()`, `.strip()` on dict values from client data
- [ ] File-exists guards for any test that opens files outside
      `agent/` or `knowledge/`
- [ ] Serialization round-trip test for any new `to_dict()` /
      `from_dict()` pair
- [ ] Documentation updated (OVERVIEW.md, ARCHITECTURE.md,
      CHANGELOG.md) for user-visible or architectural changes
