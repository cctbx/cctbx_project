# cctbx guidance for coding agents

`CONTRIBUTING.md` holds the project's contribution rules; this file adds
what coding agents most often get wrong. Where the two differ, this file
governs.

## Working rules

- **Understand before changing.** Read the code you are about to touch and
  the tests that cover it, and follow the patterns of the host file. Check
  whether the functionality already exists elsewhere in cctbx before
  writing it again, and search open and closed PRs and issues for prior
  work on the same problem; if a duplicate exists, stop and say so.
- **Solve a real problem.** Every change needs a concrete motivation: a
  failing test, a traceback, a user report, a measured slowdown, a feature
  someone asked for. A review-tool finding is a lead, not a motivation:
  reproduce it first. If asked to "fix some issues" without a specific
  problem, ask what broke before writing code.
- **One problem per change, in the smallest contained diff.** Solve it
  where the cause is, with no new abstractions, options, files, or tests
  the problem does not require. Do not bundle unrelated changes into one
  commit or pull request, and do not propose unrelated refactoring;
  improve neighboring code only where it affects the work at hand.
- **Use test-driven development.** For a bug, write the regression test
  that reproduces it and watch it fail before touching the code; for a
  feature, write the test for the behavior you are adding first. Then
  write the minimal code that makes it pass, then run the changed
  module's suite and the suites of modules that use the changed code
  (`mmtbx` after a `cctbx` or `iotbx` change; CI runs no `mmtbx` tests
  and deletes `chem_data`, so that suite only ever runs locally). Keep
  the regression test, limited to assertions about the new behavior;
  delete any other test written only to check a hypothesis about the
  cause.
- **Verify before reporting done.** Run the code, execute the tests, or
  open the generated file, and report what actually happened. If
  something fails or you skipped verification, say so plainly instead of
  presenting the work as done. Never invent a test result, a measurement,
  or a claim about behavior you did not observe. Run `git fetch` and
  compare with `origin/master` before reporting what has landed upstream.
  When reviewing, verify each finding against the code before reporting
  it, and quote error text and failing test output in full.
- **Ask before installing anything**, system-wide or into the project.
  Say what you want to install and why, then wait for an okay.
- **Commit only on an explicit request to commit, and never push.** Show
  the complete diff first; approving its content ("looks good") is not
  that request, and a request to commit is never a request to push. Get
  the same explicit approval before opening a pull request. Do not run
  `git push` in any form (plain, `--force`, `--set-upstream`, or inside a
  script), and do not let `gh` push for you: the person directing you
  pushes the branch.
- **Measure before optimizing.** Before building an optimization, profile
  a realistic run of a program that uses the code; if the code takes a
  few percent of it, leave it.

## Commits and pull requests

- Write commit messages and tidy the branch as `CONTRIBUTING.md`,
  "Commits and pull requests", describes: a scoped summary line, a short
  body (preferably a bullet list) that says why, one logical change per
  commit. Rewrite only commits that have not left your machine; once the
  branch is pushed, add commits unless the person directing you asks for
  a rebase.
- Do not add a `Co-Authored-By:` trailer or any other co-authorship or
  tool-attribution line, and no private tracker keys (Jira `UR-...`);
  GitHub issue and PR numbers are fine.
- Open the pull request (for example with `gh pr create`) once the
  branch is pushed, and fill in `.github/PULL_REQUEST_TEMPLATE.md` with
  real, specific answers. The "Coding agent" section always applies to
  you: name your model and harness there, not in commit messages. Tick
  only the checklist boxes you verified, and leave the human-review box
  to the person directing you.
- Flag interface changes and other breaking changes to the person
  directing you; announcing them is their call.

## Build and environment

This repo is **not** built with `pip install` or `setup.py` (the
`setup.py` is a stub for conda-forge packaging). The development build
uses `libtbx`'s own workflow:

- `libtbx.configure <modules...>`: run inside the build directory to set
  up the build and write the dispatchers; `bootstrap.py` does this for a
  fresh tree.
- `libtbx.refresh`: regenerate the dispatchers after adding, renaming,
  or removing a `<module>/command_line/<name>.py` entry point, or editing
  a module's `libtbx_refresh.py`.
- `make` in `<root>/build` (which runs `libtbx.scons`): rebuild the C++
  extensions after editing a `.cpp`, `.h`, or SConscript. For a new
  extension, copy the SConscript and wrapper layout of a neighboring
  package in the same module.

Never run module code with a plain `python`. Run the dispatchers by full
path into the build that belongs to this checkout, `<root>/build/bin/`,
where `<root>` is two levels above the `cctbx_project` top level, for
example `<root>/build/bin/libtbx.python`.
Each build imports the modules of its own tree, so another tree's build
tests the wrong code. **Never `source build/setpaths.sh`**: it mutates
the shell's PATH and environment. Bare dispatcher names in the examples
below stand for the full path.

Do not edit `.gitattributes` unless asked. Its thousands of per-file
`-text` entries are an SVN-migration artifact, but they pin line endings:
removing them renormalizes CRLF test data and changes Windows checkouts.

## Module dependencies

Each top-level directory is a module whose `libtbx_config` declares what
it requires. Strictly respect the dependency direction: `libtbx` must not
import `scitbx`; `scitbx` must not import `cctbx`; `cctbx` must not
import `iotbx` or `mmtbx`; `iotbx` must not import `mmtbx`. Breaking this
is a review-block. The rule covers library code; regression tests may
import downstream modules. Existing violations predate the rule: do not
copy them, and do not fix them outside a change about them.

## Coding conventions (cctbx-specific)

Key points that deviate from generic Python style:

- **2-space indentation** throughout Python, not 4. Match the host file if
  it differs.
- **80-character lines** where feasible. Long strings, URLs, and regular
  expressions may run over.
- **Python 3 only**, tracking the current CI matrix (the `python_version`
  list in `.github/workflows/syntax.yml`) rather than a fixed floor. Do
  not add `from __future__` imports; leave existing ones alone.
- **No `import *`** except inside `__init__.py` re-exporting from a C++
  extension.
- **Prefer duck typing and polymorphism to `isinstance()` chains.** A
  check against an abstract base class, or to narrow an argument that
  may legitimately be a string or a list, is fine; a ladder of concrete
  class checks that grows with every new type is not.
- **Type hints are optional.** Match the host file: add them to new
  functions only where the file already uses them, and keep them on
  signatures rather than local variables.
- **numpydoc docstrings** for every new public function and class
  (format and example in `CONTRIBUTING.md`, "Use the cctbx coding style").
- **Imports inside functions** when feasible: top-of-module imports add
  measurable startup overhead, and a package the build may lack can only
  be guarded inline. Small, cheap modules used throughout a file (`os`,
  `sys`, `math`) go at the top.
- **cctbx stands alone.** Do not import modules outside this repository
  (`phenix`, `elbow`, `phaser`, and the like) in library code or tests,
  guarded or not; code that needs an application belongs in that
  application. Existing guarded imports predate this rule: do not copy
  them. The only outside input is data from `chem_data` and
  `phenix_regression`, which may be absent: locate a file with
  `find_in_repositories(..., test=os.path.isfile)` and skip the test, or
  raise `Sorry` in library code, when it returns `None`.
- **Put general helpers in the general module** (`scitbx`, `cctbx`,
  `iotbx`, `libtbx`), not inlined into specialized code.
- **Plain over clever.** Explicit loops and sentinels rather than
  `for/else`, `dict.setdefault` as a lookup, the walrus operator, or
  nested comprehensions. No functions that only forward their arguments
  and no redundant guards; fix the root cause.
- **No new external dependencies** without discussion among the
  developers.
- **Print to `log`, not stdout** in library code: `print(msg, file=log)`,
  and prefer a `show(..., log=log)` method over inline prints. Delete
  debugging prints before committing, or gate them behind a `verbose`
  or `debug` flag that defaults to off.
- **Use `libtbx.utils.Sorry`** for user-facing errors (bad input, missing
  file). Use `assert` / standard exceptions for programmer errors.
- **Use `libtbx.easy_run`**, never `os.system` or raw `subprocess`.
  Neither `call` nor `fully_buffered` fails on its own: assert the exit
  code or call `.raise_if_errors()` (`CONTRIBUTING.md`, "Python do's").
- **Portable paths.** cctbx runs on macOS, Linux, and Windows: build
  paths with `os.path.join`, never a hard-coded `/` or `\`.

### Crystallographic conventions

- A resolution cutoff is `d_min` or `d_max` (Å); in map code `resolution`
  means the nominal map resolution and keeps that name. R-factors are
  fractions (`r_free < 0.25`), never percentages. B-factors are Å².
- Space groups are `sgtbx.space_group_info` objects. **Never compare space
  groups as strings**: two different Hall/HM symbols can represent the
  same group. Compare with `sg1.group() == sg2.group()`.

### `flex` array pitfalls

cctbx uses `scitbx.array_family.flex` arrays, **not NumPy**, pervasively.
The most common mistakes:

- `flex.double()` creates an empty array, not a scalar.
- `flex.sum(arr)` is a free function, not a method.
- Boolean selection: `arr.select(flags)`, not `arr[flags]`.
- Slicing returns a copy: `a[1:3][0] = 99` leaves `a` unchanged. Edit in
  place with `a.set_selected(selection, values)`.
- `flex` has its own random number generator: seed `random.seed(n)` and
  `flex.set_random_seed(n)` both for a reproducible test.

### Programs and PHIL parameters

Build new programs on `libtbx.program_template.ProgramTemplate` and
define their settings as PHIL scopes, as `CONTRIBUTING.md`, "Writing a
command-line program", describes. Two checks agents tend to skip:

- When modifying a program's parameter set, locate its existing
  `master_phil`, `master_phil_str`, or `master_params_str`. Agents
  hallucinate PHIL parameter names easily, so verify against the source,
  not against assumed naming.
- `fetch` silently drops a user parameter that the master does not
  define. Check input with
  `master_phil.fetch(source=user_phil, track_unused_definitions=True)`,
  which returns `(working_phil, unused)`, and assert `unused` is empty.

## Tests

Follow `CONTRIBUTING.md`, "Guidance for Developing Tests", for where a
test lives, how `run_tests.py` registers it, which data to use, and the
template. Points agents get wrong:

- Run tests from an empty scratch directory outside the checkout:
  `run_tests_parallel` refuses to start in a non-empty one, and many
  tests write fixed-name files into the current directory.

  ```sh
  # All tests for one module, in parallel, each in its own subdirectory:
  libtbx.run_tests_parallel module=mmtbx nproc=8 run_in_unique_dirs=True

  # One test file directly (fastest iteration):
  libtbx.python <root>/modules/cctbx_project/mmtbx/regression/tst_foo.py
  ```

- Nothing auto-discovers; a test that is not listed in `run_tests.py`
  never runs.
- If your change makes a test fail, fix it, even if you did not write the
  test. Never comment a test out to get around a failure, never move it to
  the expected-failure or expected-unstable list (entries there cannot
  fail the suite), and never loosen a tolerance without understanding why
  the values moved.
- Before committing, run `libtbx.find_clutter <changed files>`; CI runs
  it over the whole tree and fails on any report. Review flagged imports
  by hand (many are re-exports or side-effecting extension loads) and
  mark deliberate ones as `CONTRIBUTING.md`, "Run find clutter", shows.
