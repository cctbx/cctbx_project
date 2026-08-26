# AI Agent — LLM Programming Guidelines

Supplement to `CCTBX_LLM_PROGRAMMING_GUIDELINES.md`.
This document covers patterns and pitfalls specific to
the AI Agent codebase (`agent/`, `knowledge/`,
`programs/ai_agent.py`). Read the CCTBX guide first —
everything there applies here too.

---

## What changed in this revision

**Section numbering 1–10 is unchanged**, so existing references still
resolve. New material is added as **§0** (before §1) and **§6a**.
**Nine of the ten existing sections changed**; only §8 is untouched.
Every original paragraph, code block and checklist item is carried
over — the changes are additions and one regrouping of the checklist.

| where | change |
|---|---|
| **§0** | **NEW — the Local/Remote Parity Invariant.** It was documented only in `ARCHITECTURE.md`, ~6300 lines in. It is the single easiest rule to break silently, so it now leads. |
| **§6a** | **NEW — crossing the client-server boundary.** The two independent plumbing paths and their two distinct failure modes. Sourced from `DEVELOPER_GUIDE.md` and `ARCHITECTURE.md`. |
| §1 | note that `programs.yaml` is safety-critical configuration, not just a lookup table |
| §2 | **surveyed**: the `getLogger` rule is local to `agent/`+`knowledge/`; five sibling packages use bare `print()`. The split is now documented and the new-package question left open rather than guessed |
| §3 | the package/module name-collision trap |
| §4 | `FORCE_NO_AI_SERVER`; what `standard` mode does when it cannot run locally; parity cross-reference |
| §5 | note that a checker of *output* is not a fourth *run*-error classifier |
| §6 | parity cross-reference; prefer a server-side fix where one exists |
| §7 | pointer to §6a — session state that must reach the server has a second contract, and the touchpoint count that follows |
| §9 | two test patterns: parity, and the PHIL round-trip |
| §10 | checklist extended |

**Why this was needed, measured rather than asserted:** of the nine
identifiers that govern client-server correctness — `local_agent.py`,
`remote_agent.py`, `build_request_v2`, `apply_request_defaults`,
`SESSION_INFO_FIELDS`, `CURRENT_PROTOCOL_VERSION`,
`tst_contract_compliance`, `tst_h18_1_phil_roundtrip`,
`FORCE_NO_AI_SERVER` — **none appeared anywhere in the previous version
of this document.** All nine are in `ARCHITECTURE.md` or
`DEVELOPER_GUIDE.md`. A reader who took this supplement as their guide
to the client/server split had no way to reach them.

---

## 0. The Local/Remote Parity Invariant

**This rule outranks everything else in this document. Read it first.**

> **LocalAgent MUST produce identical results to RemoteAgent.**
> Both modes execute the same code on identical input and produce
> identical output via the documented return contract. **The only
> difference is whether information crosses machines.**

Canonical statement and the execution-split diagram:
`ARCHITECTURE.md` → *Local/Remote Parity Invariant*.

### What "identical" covers — and what it does not

**It covers the model.** The provider travels in the request
`settings` and is consumed only in shared `agent/`+`core/` code, so
local and server execution run identical branches with the **same
provider**. Running locally is *not* a licence to use a different
model, and a difference in output is *not* excused by "it ran locally."
The server simply also needs the provider's keys in its environment —
deployment, not code.

**It covers the request, byte for byte.** `LocalAgent` deliberately
performs the full encode/decode roundtrip rather than short-cutting it,
specifically so that transport bugs surface during local testing and
local mode produces byte-identical requests to what the server would
receive.

**Do not "normalise away" a parity difference.** If two paths produce
different bytes, that is the bug the roundtrip exists to expose.
Filtering timestamps or paths out of a comparison to make it pass
defeats the design.

**It does not cover observability.** Which stderr markers the operator
happens to see differs between modes; that is visibility, not parity.

**It does not cover LLM sampling variation.** Two server runs vary the
same way two local runs do. Sampling variance is orthogonal to parity
and must never be offered as the explanation for a parity failure.

### When you add a field or a setting

**You MUST update both agents:**

- `phenix_ai/local_agent.py` — `decide_next_step()`
- `phenix_ai/remote_agent.py` — `decide_next_step()`

If the agents diverge, users get different results depending on whether
`run_on_server` is True or False. **This is a silent correctness bug
and is extremely difficult to diagnose** — nothing fails, the answer is
just different for half your users.

Then see **§6a**, because updating both agents is necessary and *not
sufficient*.

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

**`programs.yaml` is safety-critical configuration, not a lookup
table.** It is simultaneously the sanitizer's allowlist, the input/slot
table, and the log-parsing pattern set. Two consequences:

- **Never regenerate it by round-tripping through a YAML dumper.** That
  silently strips every comment from the file — the `.help` text and
  the rationale for individual entries — while leaving the parsed
  content intact, so nothing fails and the documentation is gone. Apply
  changes as targeted edits.
- **State which sections a change touches.** A change to one program's
  `log_parsing` cannot affect command construction; a change to
  `strategy_flags` can. Saying which is which turns a review into a
  check.

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

**The rule above is real and it is local to `agent/` + `knowledge/`.
The rest of `libtbx/langchain` does not follow it.** Surveyed
2026-08-15 across `analysis/`, `core/`, `rag/`, `strategies/` and
`utils/`:

    logging.getLogger            0 occurrences
    bare print()                 ~120 occurrences
    print(..., file=out)         2  (utils/run_utils.py)
    strategies/                  neither -- no output calls at all

So there are **two conventions in this tree, split by package**, and
this document previously described only one. Related local idioms
worth knowing: `analysis/summarizer.py` defines a `debug_print()`
helper gated on a debug flag rather than using log levels, and
`utils/run_utils.py` uses the stream form `print(..., file=out)`.

**`print(x, file=stream)` is not forbidden anywhere.** It is an
established pattern — `print(output, file=self.logger)` with
`agent/event_formatter.py` in `DEVELOPER_GUIDE.md`, and
`utils/run_utils.py` above.

**Which convention a NEW package should follow is an open question,
not settled here.** An earlier draft of this section asserted
`logging.getLogger` for new server-side packages, reasoning that a
bare `print()` lands on the server's stdout where the client never
sees it. **That reasoning is unverified**: `analysis/` is server-side
and prints freely, so either those prints are being lost or the
reasoning is wrong, and the survey cannot tell which. The options:

- (a) follow `agent/` + `knowledge/` — `logging.getLogger`
- (b) follow `analysis/` and the rest — bare `print()`, matching five
  sibling packages
- (c) decide by role — `getLogger` for anything on the server-side
  request path, `print()` elsewhere
- (d) the stream form `print(..., file=out)`, as `utils/run_utils.py`
  does

Until this is decided, **match the package you are sitting in**, and
say in the CHANGELOG which convention a new package adopted and why.

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

### Three traps in the fallback form

**Flat intra-package imports are not a fallback.** `from routing
import Router` resolves only when the package's own directory is on
`sys.path`. The fallback form above needs the package's **parent** on
the path instead. Code that works standalone via flat imports will fail
outright once installed, and no test that runs standalone will catch it.

**Test the installed path, not just the standalone one.** The negative
control is to remove the package directory from `PYTHONPATH` and
confirm the suites still import. Without it, green tests only prove
that the path you were already using still works.

**Never name a module the same as its package.** With the parent on
`sys.path`, `pkg/pkg.py` is shadowed by `pkg/` itself and
`from pkg import Thing` raises `ImportError: cannot import name
'Thing' from 'pkg'`. The fallback then resolves to the package on one
path and the module on the other — the same name meaning two different
things depending on how you got there. Also ensure the package has an
`__init__.py`; without one it is a namespace package and is not
importable as a package at all.

---

## 4. Analysis Mode Routing

`ai_analysis.py` has five `analysis_mode` values, but
the AI Agent only uses four of them. Understanding
which modes need the server is critical:

| Mode | Used by | Needs server |
|------|---------|:------------:|
| `standard` | `phenix.ai_analysis` (standalone) | No (v121) |
| `directive_extraction` | Agent (session start) | No |
| `advice_preprocessing` | Agent (session start) | No |
| `failure_diagnosis` | Agent (on terminal error) | No |
| `agent_session` | Agent (session end) | No |

The agent **never uses `standard` mode**. That mode is
the standalone `phenix.ai_analysis` program. **As of
v121 it no longer does retrieval**: it makes one
long-context call on the whole log.

All five modes are now pure LLM calls routed locally
when `run_on_server=False` or `provider=ollama`. If you
add a new analysis mode, decide whether it needs the
RAG database and add it to `_LLM_ONLY_MODES` in
`run_job_on_server_or_locally()` if it doesn't — note
that the set currently contains every declared mode, so
the database guard there is unreachable until one does.

The agent's THINK node does its own log analysis via
`thinking_prompts.py` and the expert knowledge base —
it does not go through `ai_analysis.py` at all.

**Local routing does not relax §0.** No mode needs the database as of
v121, and all run local when asked to. They still run *identical
branches* — the only difference is whether information crosses machines.

**When a mode cannot run locally, it fails loudly.** Preserve that
principle: a silent fallback to the server is the failure this
guideline exists to prevent.

**The concrete instance is gone as of v121.** `standard` mode used to
raise `Sorry` when the provider had no local RAG database. It no longer
retrieves, so it needs no database, and that `Sorry` was refusing runs
for want of a file the code never read — observed with `anthropic`. Any
new mode that does need a database should fail the same loud way.

**`FORCE_NO_AI_SERVER=1` is the local-only switch, and it is
absolute.** Set to the literal `1` (whitespace-stripped), it forces
`run_on_server=False` for both `phenix.ai_agent` and
`phenix.ai_analysis`, across **every** analysis mode. Both dispatchers
re-check it and refuse to submit to the server even via the
no-local-database fallback `standard` mode used to take. The
shipped PHIL default stays `True`; this is a surgical switch, not a
baseline change. **Use it to drive the local side of a parity test**
rather than building a separate harness.

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

**These three classify a *run*.** A checker that classifies *output* —
whether a produced report or command is supported by its inputs — is a
different kind of thing and does not belong in this list. If you add
one, say so in its docstring, so a later consolidation pass does not
sweep it in.

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

**The governing principle:** *decisions and knowledge → server;
execution and I/O → client.* A user with yesterday's PHENIX install
must be able to connect to today's server and get correct results.

**Client-side changes are the expensive ones.** Before adding one, ask
whether the same outcome can be reached server-side. A change that
makes a client-side fix unnecessary — for example, having the producer
of an artifact record what a consumer would otherwise have to infer —
reaches every user immediately, while the client-side version reaches
only those who reinstall.

---

## 6a. Crossing the client-server boundary

**There are two independent plumbing paths, with two different failure
modes. Both are silent.**

### Path 1 — PHIL parameters: the dual-schema mirror

`programs/ai_agent.py` (client) and `programs/ai_analysis.py` (server)
parse **independent** `master_params` schemas. **Adding a parameter to
one does not propagate to the other.**

**Failure mode.** The client does `copy.deepcopy(self.params)` and
assigns the new attribute. If it is not declared in the *client's*
`master_params`:

```
AttributeError: Assignment to non-existing attribute
  "ai_analysis.your_new_field"
```

**In production this is caught by the surrounding try/except and
silently swallowed** — directive extraction returns `{}`, the agent
ignores the parameter, and behaviour drifts in ways that may take a
production run to notice. This is the recorded **H18 → H18.1**
production failure.

**Checklist for a new PHIL parameter flowing client → server:**

1. Declare it in `programs/ai_agent.py:master_params` **and**
   `programs/ai_analysis.py:master_params`.
2. Add a `.help` comment in each noting the cross-file mirror
   requirement.
3. Add a `tst_*_phil_roundtrip.py` K-test exercising
   parse → extract → deepcopy → assign (see §9).
4. Run sandbox tests; the source-grep variant should pass.
5. Deploy and run the live-PHIL variant.
6. **Both files in the same commit.**

### Path 2 — `session_info`: the wire whitelist

Any value `ai_agent.py` puts into the per-cycle `session_info` dict
that server-side code needs to read must be carried at **every** hop:

```
ai_agent.py session_info["X"]
  → build_session_state()      [api_client.py]  ← THE WIRE WHITELIST
  → build_request_v2()         [api_client.py]  (passes session_state through)
  → create_request()           [api_schema.py]
  → apply_request_defaults()   [api_schema.py]  (additive only)
  → transport encode/decode    [transport.py]
  → run_ai_agent.py            (explicit per-field map-back)
  → graph_nodes.py PERCEIVE    session_info.get("X")
```

**The whitelist is `build_session_state()`** — an explicit per-field
copy; unlisted keys never enter `session_state` at all.

**`build_request_v2()` used to hold a SECOND whitelist**, twenty lines
later in the same file, re-enumerating everything the first had just
built. Eleven fields were lost there across two audits, including two
fixes that had been made and marked verified —
`client_protocol_version` and `log_program`, the latter being the whole
purpose of protocol v9, which had therefore never once worked. It now
passes `session_state` through, applying only coercions and derived
flags. **Do not reintroduce an enumeration there.**

`build_request_v2()` remains the **parity chokepoint**: both
`LocalAgent` and `RemoteAgent` call it exactly once and add nothing to
`session_state` afterwards. Anything added in either agent instead
would break `server == local`, and the failure would be invisible in
whichever mode was not being exercised.

**Failure mode, and note how it differs from Path 1.** An omission here
drops the field **identically for local and server**, so **parity stays
intact and the feature is simply inert.** Nothing diverges, nothing
errors, and the code you wrote never runs. In v120 three plan-driven
fields — `plan_has_pending_stages`, `plan_next_stage_programs`,
`plan_current_unrun_lead_program` — round-tripped nowhere for exactly
this reason.

**A new `session_info` field requires, at minimum:** an entry in
`build_session_state()`, an entry in the `build_request_v2()`
whitelist, and the `run_ai_agent.py` map-back. If it is read in
`graph_nodes.py` via the literal `session_info.get("X")` (not the `_si`
alias), it must **also** be registered in
`agent/contract.py::SESSION_INFO_FIELDS`, with
`CURRENT_PROTOCOL_VERSION` bumped to ≥ the field's version.
`tst_contract_compliance.py` enforces this.

**If a change requires the client to send new data or handle a new
response format, that is a protocol change** — see
`DEVELOPER_GUIDE.md §8 (Backward Compatibility & Contract)` for the
version-bump procedure.

### This is now enforced, because documenting it was not enough

The failure above is recorded here for v120 — three plan-driven fields
that "round-tripped nowhere". It then happened four more times, in an
audit of a single week's work:

| field | what was missing |
|---|---|
| `log_program` | registered in the contract, read by the server, carried by neither allow-list |
| `state["warnings"]` | written by `perceive()`, printed by the client, never passed to `create_response` |
| `client_protocol_version` | set by the client, carried by nothing — so the server read the default `1` for **every** client, and `get_deprecation_warnings` fired for everyone on every cycle |
| `mtz_rfree_map` | the v120.2 **parity fix** itself: forwarded into `session_state`, never rebuilt into `session_info`, so `command_builder` saw `None` and the guard it exists to make identical behaved differently server-side |

Plus a rename that no test could see: the per-cycle metrics dict was
translated **to** the declared name `analysis` by `session.py`, then
renamed **back** to `metrics` one hop later, leaving 8 of 9 consumers
reading a key that never arrived.

**`tests/tst_transport_surfaces.py` now checks this mechanically** — for
every field in every declared surface: a producer exists, each
allow-list hop carries it, a consumer exists, and no key is renamed in
flight. It reads the sources as text, so it runs anywhere, and it
**fails rather than skips** when a source is unreadable.

Run it whenever you add or move a field across the boundary. Six
occurrences of one failure mode, all in code whose author had this
section available, is the argument for a check rather than a paragraph.

### On `except` blocks that produce a value

Related, and the reason several of the above stayed invisible:

> **An `except` that swallows a failure while still producing a value
> the agent acts on must record that it did.**

`graph_nodes.py` has ~18 such sites, most of them fine. The one to
learn from converts `float(_m["r_free"])` failing into `None`, which
downstream reads as *"no R-free"* rather than *"R-free unparseable"* —
and those lead to different decisions. If the value feeds a stop check,
a plateau test, or a command, log the swallow.

### The rule of thumb

**Parity (§0) and reachability (§6a) are different properties, and
passing one tells you nothing about the other.** Both agents updated
and the field not whitelisted → parity holds, feature dead. Field
whitelisted and only one agent updated → feature live, results diverge
by mode. Check both.

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

**These three places cover persistence across cycles on one side of the
wire. If the field must also reach server-side code, §6a Path 2
applies as well** — and the whitelist omission is invisible to all
three checks above.

**Count the cost of a new state field before adding one.** Persisting
across cycles costs the three places above; if it must also reach
server-side code, §6a Path 2 adds the `build_request_v2()` whitelist
and the `run_ai_agent.py` map-back, and `agent/contract.py` may add a
sixth. Four to six touchpoints, each silent when missed. That is the
number to weigh — not a reason to avoid persistence, which is what
session resume is built on.

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

See `docs/CHANGELOG.md` for examples.

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

**A skip is not a pass.** A guard that fires because a path is wrong
looks exactly like one that fires because a directory is legitimately
absent. Suites should report what was skipped and why, and a skip count
that grows unexplained should be investigated rather than accepted.

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

The same trap appears in regex form. A body pattern like
`(?:\s{4,}.*\n)+` after a handler runs to the end of the function,
because every following line is also indented — so a check that looks
for a problem *inside* the handler ends up inspecting the whole
function and passing. **Delimit by indentation, not by a greedy
quantifier**, and prove it by planting the violation you expect to
catch.

### PHIL round-trip tests

For any new PHIL parameter crossing the boundary (§6a Path 1), a
source-grep variant that runs in sandbox plus a live-PHIL variant that
skips gracefully:

```python
def test_phil_roundtrip_for_new_param():
    with open("programs/ai_agent.py") as f:
        content = f.read()
    start = content.find('master_params = """')
    end = content.find('master_phil = libtbx.phil.parse')
    assert "your_new_field" in content[start:end]

    try:
        import libtbx.phil, copy
    except ImportError:
        return
    ...
    directive_params.ai_analysis.your_new_field = "test"
    assert directive_params.ai_analysis.your_new_field == "test"
```

See `tests/tst_h18_1_phil_roundtrip.py`.

### Parity tests

Any change touching the request pipeline should be exercised **both
ways** and compared. Drive the local side with `FORCE_NO_AI_SERVER=1`
(§4) rather than a bespoke harness.

- Assert the built request is **byte-identical** between the two paths.
- Assert the return contract matches.
- **Do not** normalise differences away to make the comparison pass.
- Include a negative control: break the local path deliberately and
  confirm the test fails. A parity test that has never failed has not
  been shown to work.

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

Note the distinction: a handler that catches and **reports** a failure
is fine, and so is one where the exception is the expected outcome. The
rule is about a *runner* discarding a suite's result.

---

## 10. Checklist (Agent-Specific Additions)

In addition to the CCTBX checklist:

**Parity and the boundary**

- [ ] **Both agents updated** — `local_agent.py::decide_next_step()`
      and `remote_agent.py::decide_next_step()` (§0)
- [ ] **Exercised both ways** and compared, with
      `FORCE_NO_AI_SERVER=1` driving the local side (§0, §4)
- [ ] No parity difference explained away by "it ran locally", by
      sampling variance, or by normalising the comparison (§0)
- [ ] New PHIL parameter mirrored in **both** `master_params`, with a
      round-trip test, in one commit (§6a Path 1)
- [ ] New `session_info` field present at **every** hop, including the
      `build_request_v2()` whitelist and the `run_ai_agent.py`
      map-back (§6a Path 2)
- [ ] `agent/contract.py::SESSION_INFO_FIELDS` registered and
      `CURRENT_PROTOCOL_VERSION` bumped if read via a literal
      `session_info.get()` (§6a)
- [ ] Protocol change? `DEVELOPER_GUIDE.md §8` version-bump procedure
      followed (§6a)

**Everything else**

- [ ] New strategy flags added to `programs.yaml`
      before use in command-building code (§1)
- [ ] `programs.yaml` edited in place, never regenerated through a YAML
      dumper; which sections a change touches is stated (§1)
- [ ] New error patterns checked against all three
      classification systems (§5)
- [ ] New `analysis_mode` values added to
      `_LLM_ONLY_MODES` if they don't need the RAG DB, and they fail
      loudly rather than falling back to the server (§4)
- [ ] Client vs server code path identified for any
      changes to `programs/ai_agent.py` — and the server-side
      alternative considered first (§6)
- [ ] Session state fields appear in `session.data`,
      `create_initial_state()`, and `contract.py` (§7)
- [ ] Cross-module imports use the dual-path fallback; no flat
      intra-package imports; no module named after its package (§3)
- [ ] Import fallbacks tested with the package directory **off**
      `PYTHONPATH` (§3)
- [ ] Output convention matches the package being edited: `getLogger`
      in `agent/`+`knowledge/`, `print()` elsewhere. A new package
      states its choice in the CHANGELOG (§2)
- [ ] File-exists guards for tests opening files
      outside `agent/` or `knowledge/`; skips reported, not absorbed (§9)
- [ ] Source-scanning tests use function-boundary
      windows, not fixed character counts or greedy quantifiers (§9)
- [ ] CHANGELOG entry follows the standard table
      format (§8)
- [ ] Documentation updated (OVERVIEW.md,
      ARCHITECTURE.md, DEVELOPER_GUIDE.md) for
      user-visible or architectural changes (§10)

---

## Appendix — where the authoritative text lives

This supplement summarises; it does not replace. When the two
disagree, the source document wins.

| topic | authoritative source |
|---|---|
| Local/Remote Parity Invariant, execution split | `ARCHITECTURE.md` — *Local/Remote Parity Invariant*, *Client-Server Update Model* |
| `session_info` plumbing contract, wire whitelist | `ARCHITECTURE.md` — *`session_info` field plumbing contract (parity-critical)* |
| `FORCE_NO_AI_SERVER` | `ARCHITECTURE.md` — *FORCE_NO_AI_SERVER (v120 P1)* |
| dual `master_params`, H18 → H18.1 | `DEVELOPER_GUIDE.md` — *Mirroring PHIL on client and server* |
| protocol version bumps | `DEVELOPER_GUIDE.md §8` — *Backward Compatibility & Contract* |
| error-classifier overlap and consolidation | `ARCHITECTURE.md` — *Design tensions* |
| general cctbx patterns | `CCTBX_LLM_PROGRAMMING_GUIDELINES.md` |
