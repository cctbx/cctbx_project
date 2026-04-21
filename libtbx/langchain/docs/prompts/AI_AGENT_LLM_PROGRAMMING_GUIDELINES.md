# AI Agent — LLM Programming Guidelines

Supplement to `CCTBX_LLM_PROGRAMMING_GUIDELINES.md`.
This document covers patterns and pitfalls specific to
the AI Agent codebase (`agent/`, `knowledge/`,
`programs/ai_agent.py`). Read the CCTBX guide first —
everything there applies here too.

---

## 1. Parameter Verification Against programs.yaml

LLMs frequently hallucinate PHENIX command-line
arguments. The agent has a specific defense: the
command sanitizer strips any parameter not in the
program's `strategy_flags` allowlist in
`knowledge/programs.yaml`.

When generating or modifying agent code that builds
commands:

- Check `programs.yaml` for the program's
  `strategy_flags` section. If a parameter is not
  listed, the sanitizer will strip it — your code
  will silently have no effect.
- If you need a new strategy flag, add it to the
  YAML first, then use it in code.
- The YAML also defines `inputs` (required/optional
  file slots), `invariants` (auto-fill rules), and
  `defaults` (always appended). Understand all four
  sections before modifying command building.

Example: to add `rebuild_in_place=False` for
autobuild, you must first verify that
`rebuild_in_place` is in autobuild's `strategy_flags`.
If it isn't, the sanitizer strips it and the parameter
never reaches PHENIX.

---

## 2. Logging Conventions

The codebase uses three logging patterns in different
contexts. Don't mix them.

**`agent/` and `knowledge/` modules:**
Use `logging.getLogger(__name__)`:

```python
logger = logging.getLogger(__name__)
logger.debug("update failed", exc_info=True)
```

**`programs/ai_agent.py`:**
Use `self.vlog` at the appropriate verbosity level:

```python
self.vlog.quiet("Fatal: %s" % error)   # always
self.vlog.normal("Decision: %s" % prog) # default
self.vlog.verbose("File list: %s" % f)  # detail
```

**Test files:**
Use `print()`:

```python
print("  PASS: test_name")
print("  SKIP (ai_agent.py not found)")
```

No `print()` in `agent/` modules.  No `logger` in
test files.

---

## 3. Import Fallbacks

The agent code runs in two contexts: inside PHENIX
(where modules live under `libtbx.langchain`) and
standalone (where they're imported directly). All
cross-module imports in `agent/` and `knowledge/`
must have a fallback:

```python
try:
  from libtbx.langchain.agent.structure_model \
    import StructureModel
except ImportError:
  from agent.structure_model import StructureModel
```

Never use bare `from libtbx.langchain.X import Y`
without a fallback `except ImportError` block.

This is different from the general cctbx pattern
(where optional imports set the name to `None`).
In the agent, both paths must resolve to the same
class — the code that follows uses it unconditionally.

---

## 4. Analysis Mode Routing

`ai_analysis.py` has five `analysis_mode` values, but
the AI Agent only uses four of them. Understanding
which modes need the server is critical:

| Mode | Used by | Needs server |
|------|---------|:------------:|
| `standard` | `phenix.ai_analysis` (standalone) | **Yes** (RAG DB) |
| `directive_extraction` | Agent (session start) | No |
| `advice_preprocessing` | Agent (session start) | No |
| `failure_diagnosis` | Agent (on terminal error) | No |
| `agent_session` | Agent (session end) | No |

The agent **never uses `standard` mode**. That mode is
the standalone `phenix.ai_analysis` program which does
retrieval-augmented generation against the Phenix
knowledge base.

The four agent modes are pure LLM calls routed locally
when `run_on_server=False` or `provider=ollama`. If you
add a new analysis mode, decide whether it needs the
RAG database and add it to `_LLM_ONLY_MODES` in
`run_job_on_server_or_locally()` if it doesn't.

The agent's THINK node does its own log analysis via
`thinking_prompts.py` and the expert knowledge base —
it does not go through `ai_analysis.py` at all.

---

## 5. Three Error Classification Systems

Error handling is split across three independent
systems. When adding a new error pattern, you must
check all three to ensure they agree:

1. **`agent/error_classifier.py`** — called by
   PERCEIVE at the start of the next cycle. Five
   categories: TERMINAL, PHIL_ERROR, AMBIGUOUS_PHIL,
   LABEL_ERROR, RETRYABLE. Feeds `should_pivot()`.

2. **`agent/error_analyzer.py`** — YAML-driven.
   `recoverable_errors.yaml` for auto-fixable errors
   (ambiguous labels). `diagnosable_errors.yaml` for
   terminal errors needing LLM diagnosis.

3. **`programs/ai_agent.py::_classify_error()`** —
   oldest classifier. Maps to INPUT_ERROR (agent's
   fault, don't count) or REAL_FAILURE (count it).
   Used only for cycle history bookkeeping.

See ARCHITECTURE.md "Design tensions" for the full
overlap analysis and recommended consolidation path.

---

## 6. Client-Server Code Path Awareness

Changes to files in `agent/` and `knowledge/` are
server-side — they take effect immediately for all
users. Changes to `programs/ai_agent.py` are
client-side — users must update their install.

When modifying `programs/ai_agent.py`, identify
whether your change is in the **top half** (client
code: `run()`, `_run_single_cycle()`,
`_inject_user_params()`) or the **bottom half**
(server code called when `run_on_server=False`).
Top-half changes require user updates; bottom-half
changes are effectively server-side.

See ARCHITECTURE.md "Client-Server Update Model"
for the full execution split diagram.

---

## 7. Session State Persistence

Every piece of state the agent needs across cycles
must appear in three places:

1. Written to `session.data["key"]` after the cycle
2. Read from `session.data["key"]` on resume
3. Included in `create_initial_state()` if passed
   through the graph

When you add a new state field, grep for all three
patterns and verify the field appears in each. Also
register new `session_info` fields in
`agent/contract.py` with a default value.

---

## 8. CHANGELOG Format

CHANGELOG entries for the agent follow a specific
format:

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
N new tests.
Run with: `python3 tests/run_all_tests.py`
```

See `docs/project/CHANGELOG.md` for examples.

---

## 9. Test Patterns Specific to the Agent

### File-exists guards

Tests that open files outside `agent/` and
`knowledge/` (e.g., `programs/ai_agent.py`,
`wxGUI2/Programs/AIAgent.py`) MUST guard with
`os.path.isfile()` because these directories may
not be present in all deployment contexts:

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

Files in `agent/` and `knowledge/` are always
present — no guard needed.

### Source-scanning search windows

Tests that scan source code for specific patterns
must account for function length:

```python
# WRONG — function may be longer than 5000 chars
idx = src.find("def my_function")
assert "expected" in src[idx:idx + 5000]

# RIGHT — search to the next function definition
idx = src.find("def my_function")
next_def = src.find("\ndef ", idx + 1)
assert "expected" in src[idx:next_def]
```

### Don't swallow test failures

`run_tests_with_fail_fast()` raises `AssertionError`
on failure. If your `run_all_tests()` wraps it in a
try/except, failures are silently swallowed:

```python
# WRONG — test failures discarded
def run_all_tests():
  try:
    run_tests_with_fail_fast()
  except Exception:
    pass

# RIGHT — let failures propagate
def run_all_tests():
  run_tests_with_fail_fast()
```

---

## 10. Checklist (Agent-Specific Additions)

In addition to the CCTBX checklist:

- [ ] New strategy flags added to `programs.yaml`
      before use in command-building code
- [ ] New error patterns checked against all three
      classification systems
- [ ] New `analysis_mode` values added to
      `_LLM_ONLY_MODES` if they don't need the RAG DB
- [ ] Client vs server code path identified for any
      changes to `programs/ai_agent.py`
- [ ] Session state fields appear in `session.data`,
      `create_initial_state()`, and `contract.py`
- [ ] File-exists guards for tests opening files
      outside `agent/` or `knowledge/`
- [ ] Source-scanning tests use function-boundary
      windows, not fixed character counts
- [ ] CHANGELOG entry follows the standard table
      format
- [ ] Documentation updated (OVERVIEW.md,
      ARCHITECTURE.md, DEVELOPER_GUIDE.md) for
      user-visible or architectural changes
