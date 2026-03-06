# Policy for LLM-Assisted Development in PHENIX/cctbx

This document describes the workflow that PHENIX developers
agree to follow when using LLMs (Claude, Gemini, ChatGPT,
etc.) to write or modify code in the PHENIX/cctbx ecosystem.

The goal is not to restrict LLM use — it is to ensure that
LLM-generated code meets the same standards as hand-written
code: correct, tested, portable, and maintainable.

---

## 1. Session Setup

Before asking the LLM to write any code, provide it with
context:

**Required context (always attach):**

- `CCTBX_LLM_PROGRAMMING_GUIDELINES.md` — the companion
  document to this policy. It contains the coding standards,
  patterns, pitfalls, and checklist that the LLM should
  follow. Attaching it at the start of every session is the
  single highest-leverage step you can take.
- The existing source file(s) to be modified, in full. LLMs
  produce better code when they can see the surrounding
  style, imports, and conventions.

**Recommended context (attach when relevant):**

- `OVERVIEW.md` or `ARCHITECTURE.md` from the relevant
  subsystem — grounds the LLM in the system design so it
  doesn't reinvent existing infrastructure.
- The test file(s) that cover the code being modified —
  helps the LLM understand what's already tested and write
  compatible new tests.
- Error messages or log output that motivated the change —
  gives the LLM concrete evidence of the problem.

**Never rely on the LLM's training data for PHENIX
specifics.** LLMs hallucinate PHIL parameters, API
signatures, and command-line flags that don't exist. Always
provide the actual source code or YAML definitions as
context rather than asking the LLM to recall them.

---

## 2. Plan-Verify-Execute Workflow

All LLM-assisted development follows a three-phase cycle.
Do not skip phases, even for "simple" changes — the bugs
that escape review are always the ones that looked simple.

### Phase 1: Plan

Ask the LLM to produce a markdown plan before writing any
code. The plan must include:

- **Problem statement**: What is broken or missing, with
  specific evidence (error messages, failing test output,
  user report).
- **Approach**: How the fix or feature will work, at a
  level of detail sufficient for another developer to
  review. For scientific code, this includes the
  mathematical formulas or algorithms being implemented.
- **Dependencies and edge cases**: What existing code is
  affected, what inputs could be unusual (None values,
  empty files, Windows paths, resumed sessions), and how
  each is handled.
- **Implementation steps**: A numbered sequence of small,
  testable increments.
- **Test plan**: Which new tests will be written, which
  existing tests might be affected, and how to verify
  correctness.

**Review the plan before proceeding.** Read it yourself.
Optionally, show it to a second LLM and ask: *"What could
go wrong with this approach?"* This cross-check catches
design-level issues before any code is written.

### Phase 2: Implement and verify

The LLM implements the plan one step at a time. After each
step:

1. **The LLM self-reviews** (this is built into the
   guidelines document — the LLM is instructed to check
   for edge cases, scoping issues, and style before
   presenting each step).
2. **You ask**: *"Are there any other fixes or additions
   you would like to make to this change before we move
   on?"* This prompt is the single most effective quality
   gate. It triggers the LLM to re-examine its own work
   and consistently catches bugs that the initial
   implementation misses.
3. **The LLM runs tests** and confirms they pass before
   presenting the step.
4. **You review** the code. Check that it matches the plan,
   follows the coding standards, and makes sense
   scientifically.

Do not batch multiple steps and review at the end. The
cost of finding a bug in step 2 while you're reviewing
step 5 is much higher than catching it immediately.

### Phase 3: Final audit

After all steps are complete, ask the LLM to run a full
audit:

- Parse-check every changed file
- Run all relevant test suites
- Check for line-length violations in changed regions
- Verify the checklist in the guidelines document (Section
  7 of the AI Agent guidelines, Section 6 of the general
  cctbx guidelines)

Then ask: *"Check every changed file for the patterns in
the Common Pitfalls section."* This catches the classes of
bugs that LLMs produce most frequently (None-safety,
path separators, floating point comparison, parameter
hallucination).

---

## 3. What the Programmer Is Responsible For

The LLM is a tool. The programmer is responsible for the
code that gets committed. Specifically:

### Scientific correctness

LLM-generated code is a hypothesis, not a proof. The
programmer must verify:

- **Algorithms match the literature.** If the code
  implements a crystallographic calculation (structure
  factor computation, map coefficient generation, peak
  searching), check the math against the relevant
  reference, not just the LLM's assertion that it's
  correct.
- **Units are consistent.** Distances in Ångströms, angles
  in degrees, R-factors as fractions (0.25, not 25),
  B-factors in Å². An LLM that has seen both conventions
  in its training data may use either one.
- **Coordinate-space operations are cache-aware.** Spatial
  data structures (KD-trees, neighbor lists) must be
  invalidated when the unit cell or symmetry changes.
  LLMs build correct trees but frequently omit the
  invalidation logic.
- **Numerical stability.** Division by zero, log of zero,
  very small denominators in correlation coefficients.
  LLMs rarely add these guards unprompted.

### Platform correctness

PHENIX runs on macOS, Linux, and Windows. The programmer
must verify:

- No hardcoded path separators (`/` or `\`) — only
  `os.path.join()`.
- Explicit `encoding='utf-8'` on file operations where
  the content may contain non-ASCII characters.
- No reliance on Unix-specific behavior (signal handling,
  symlinks, case-sensitive filesystems).

### Test adequacy

The programmer must verify that the LLM's tests actually
test what they claim to:

- **Coverage of edge cases.** Does the test exercise None
  inputs, empty lists, missing keys, and boundary values?
  LLMs tend to write tests for the happy path.
- **Meaningful assertions.** `assert result is not None`
  doesn't verify correctness — it verifies existence.
  Tests should check specific values, ideally against
  known crystallographic results.
- **Floating point comparisons.** Tests must use
  `approx_equal()` (from `libtbx.test_utils`), never
  `==`, for any value derived from a calculation.

### Commit quality

Before committing LLM-generated code:

- Run the full test suite, not just the tests the LLM
  wrote. Changes in shared modules can break downstream
  code that the LLM never saw.
- Review the diff, not just the final file. LLMs
  occasionally make unrelated "cleanup" changes that
  alter behavior.
- Write a commit message that describes the change, not
  the process. "Fix None crash in space group lookup"
  is a good message. "Changes suggested by Claude" is
  not.

---

## 4. Coding Standards (Summary)

The full coding standards are in the
`CCTBX_LLM_PROGRAMMING_GUIDELINES.md` document that gets
fed to the LLM. The programmer should be familiar with
these and enforce them during review:

**Style**: 2-space indentation, 80-character line width,
descriptive names, no trailing whitespace, no unused
imports.

**None-safety**: `(d.get("key") or "").lower()`, never
`d.get("key", "").lower()` — the latter crashes when the
key exists with value `None`.

**Error handling**: No bare `except:` (catches
`SystemExit`). Every `except Exception: pass` gets a
comment explaining why. Non-critical functions use the
"Never raises" pattern with logging.

**Paths**: `os.path.join()` always. Never string
concatenation with `/`.

**Tests**: Fail-fast style, `approx_equal` for floats,
explicit edge-case coverage.

**Serialization**: `to_dict()` / `from_dict()` must
round-trip. `from_dict()` must tolerate missing keys.

---

## 5. When NOT to Use LLMs

LLMs are effective for:

- Implementing well-defined algorithms with clear
  specifications
- Writing test suites for existing code
- Refactoring for consistency (renaming, reformatting,
  extracting functions)
- Debugging with specific error messages as input
- Writing documentation for existing code
- Boilerplate (PHIL definitions, serialization, CLI
  wrappers)

LLMs are unreliable for:

- **Novel scientific algorithms** where correctness
  depends on domain knowledge the LLM may not have.
  Use LLMs to implement the algorithm you've designed,
  not to design it.
- **Performance-critical inner loops** where the LLM
  can't profile or benchmark. Write the inner loop
  yourself; use the LLM for the surrounding
  infrastructure.
- **Security-sensitive code** (authentication, file
  permissions, subprocess sandboxing). LLMs produce
  code that works but may have subtle vulnerabilities.
- **Guessing at APIs they haven't seen.** If the LLM
  doesn't have the source code in context, it will
  invent plausible-looking function signatures that
  don't exist.

In these cases, write the critical code yourself and use
the LLM for the non-critical surrounding code (tests,
documentation, boilerplate).

---

## 6. Prompt Reference Card

Keep these prompts handy during LLM sessions:

| When | Prompt |
|------|--------|
| Start | Attach guidelines + source files |
| Before coding | *"Create a plan as a markdown file"* |
| Plan review | *"What could go wrong with this?"* |
| Each step | *"Any other fixes before we move on?"* |
| Each step | *"Run all relevant tests"* |
| All done | *"Run a full audit of all changes"* |
| Final | *"Check changed files for Common Pitfalls"* |

---

## 7. Checklist for Every Pull Request

Before committing LLM-assisted code, the programmer
confirms:

- [ ] A plan was created and reviewed before coding began
- [ ] Each step was reviewed individually, not batched
- [ ] "Any other fixes?" was asked after each step
- [ ] All changed files parse cleanly
- [ ] The full test suite passes (not just new tests)
- [ ] The diff contains only intentional changes
- [ ] No hardcoded path separators
- [ ] No bare `except:` statements
- [ ] No `d.get("key", "").lower()` on potentially-None
      values
- [ ] `approx_equal()` used for all floating-point
      comparisons in tests
- [ ] Scientific units are correct (Å, fractions, degrees)
- [ ] Documentation updated for user-visible changes
- [ ] Commit message describes the change, not the tool
