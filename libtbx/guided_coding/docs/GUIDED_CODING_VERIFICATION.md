# GuidedCoding r10 rev16 candidate verification

**Current status (2026-10-01):** r10 rev16 development handoff from the
verified supplied rev14, skipping failed rev15. The Rev16 section separates
the retained setup behavior, bounded local checks and procedural simulations
from the actual-client setup gate, which remains pending. Previous release
observations establish no rev16 result. Live reconciliation, outside review
and Developer acceptance/integration/activation/publication are not complete.

**Historical status (2026-09-29):** a proposed kit, committed locally in the
Developer's `cctbx_project` and trialled on the intended Mac (Rev12
section below); not published and not approved for Worker use.

SUPERSEDED 2026-09-29 (rev12; see the Rev12 section below): "This describes a proposed kit. It has not been installed in the Developer's
Claude Code Desktop app, integrated into `cctbx_project`, or approved for
Worker use."

## Release and source checks

From the canonical `libtbx/guided_coding` root, first compare the downloaded
archive with the independently supplied SHA-256 in the review cover.
Replace the path below with that canonical root and run this whole Bash
block. Each dependent command runs only if the preceding command succeeded:

```bash
cd /absolute/path/to/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev16 candidate' payload/RELEASE
```

Require every listed hash to pass, `VERIFIED complete source`, and the
`r10 rev16` release label. After a successful directory change, the listed-hash check establishes
the checker bytes before running it. The complete-source check compares the entire directory
with the manifest, including unlisted files, and refuses links, hardlinks,
missing files and changed bytes. `-I` prevents a module in the tools
folder or working directory from shadowing Python's standard library;
`-B` avoids writing bytecode into the source tree. Any failure returns nonzero and skips the rest of this block. Do not
continue into another block or dependent preparation after failure. Repeat
this guarded verification before running central tools, integration, or publication. Run the tools with
`python3 -I -B`. These checks establish internal consistency with the
supplied manifest; they do not authenticate the supplier. A separately
stated archive hash and an independent review are still needed.

## CLI version gate and rev11 registration

Only after the whole source-check block succeeds, and before creating the personal skill link or
using `/gc setup` or `/gc <task>`, run:

```bash
python3 -I -B payload/tools/screen_check.py check-claude-version
```

The command looks up `claude` on this shell's PATH, calls `--version`, and
requires a readable version at least 2.1.281. It refuses missing,
unrecognized, failing, and older executables without installing updates.
The cutoff accounts for the rev8 Mac trial's 2.1.268 `AGENTS.md` loading
failure and [Claude Code's documented limitations before 2.1.281](https://code.claude.com/docs/en/memory#when-agentsmd-support-is-unavailable).
Tests use disposable
fake executables: 2.1.268 and 2.1.280 fail; 2.1.281, the Developer's
reported 2.1.284, and 2.2.0 pass. They also cover missing, malformed and
nonzero results. These are author-side tests on Linux, not a live Mac
installation or proof that Claude Code Desktop loads `AGENTS.md`.
The gate trusts the stdout of whichever `claude` executable is found on
PATH. A misleading executable can report a passing version. Inspect the
reported executable path in the live setup; this gate does not authenticate
it or the Desktop app.
The rev9 user-guide command had a defect: when `skills/guided_coding` already
existed as a directory or directory symlink, `ln -s` created a nested link
inside it and returned success. With an existing link to the source, that
nested link invalidated the complete-source inventory. This was reproduced
in a disposable setup. Rev10 replaces the shell `ln` with
`screen_check.py register-skill .`, which rechecks source and CLI, creates
the exact link exclusively, and refuses an occupied path without writing a
nested link. Its unit cases cover fresh registration, a second run, a link
to another package, a real directory, a dangling link, a symlinked skills
parent, and an older CLI. The rev10 author-side Linux suite ran **41 tests:
40 passed, one Mac case-insensitive test skipped**. Its reviewer found a
further case: an absolute configuration path with a missing prefix followed
by `..` let recursive directory creation change the effective path after
preflight. It could create an empty directory in the source before refusing
an occupied destination, or write through a hidden symlinked parent. Rev11
rejects any `..` component before directory creation. Its new unit case
snapshots source files and directories and the other directory across fresh,
occupied, configuration-link and skills-link endpoints. The rev11 author-side
Linux suite ran **42 tests: 41 passed, one Mac case-insensitive test skipped**.
No real Claude Code CLI or Desktop run was made for this rev11 candidate
(true when written; SUPERSEDED 2026-09-29 (rev12; see the Rev12 section below)).
Read-only `/gc help` and `/gc status`, and the personal-link removal
`/gc uninstall`, remain usable if the version gate fails. The one-time
conversational setup is instructed to run the
same command before it creates a link. The registration helper enforces its
own checks even when called directly; there is still no per-project installer.

## Rev7 review findings addressed in rev8

| Reviewer probe | Rev7 result | Rev8 check and result |
| --- | --- | --- |
| Add unlisted `payload/tools/hashlib.py` | Listed-hash command passed and the module executed on `screen_check.py --help` | Complete inventory refuses the extra file before guided work. A disposable CLI test confirms neither tool's direct `--help` loads it. |
| Put companions in `packet/companions/` | Bundle helper accepted them in both domains | Bundle helper refuses before writing; packet verification still passes. It also refuses reserved companion filenames already included in a packet, even if another companion directory is outside. |
| Send source archive without release proof | Reviewer could not establish gate provenance or packet identity | Rev8 review handoff includes source archive, stated archive SHA-256, complete diff from rev7, frozen packet, separate companions, and author-side criterion/runner record. This is evidence to review, not an approval. |

The rev8 suite reported 36 tests, OK, one skipped on Linux with Python 3.12.
The skipped case needs a case-insensitive Mac volume. The new tests cover
clean and changed source, extra code, extra ordinary file and symlink, direct
CLI loading, and nested companions. An extraction check compares the source
archive's files to the candidate and reruns the source verifier. The review
packet records the exact commands, outputs and exit statuses. Author-side
checks are not an independent audit.

## Existing behavior and limits

The personal `~/.claude/skills/guided_coding` symlink is intended to expose
one central, manual `/gc` skill to local projects. A target's current
instructions must adopt the exact Developer–Guide Contract identity and
identify its own project method before GuidedCoding governs that target.
`/gc help`, `status`, `setup` and `uninstall` are control requests; setup and
adoption are conversational instructions, not an executable installer.
An ordinary task in a fresh conversation remains ordinary after older
always-on instructions are retired. These are design claims checked against
the files; they are not live Desktop observations.

The screen checker verifies a frozen evidence packet's listed bytes and
exact inventory. For that packet's identity, calculate SHA-256 of its
`MANIFEST.sha256`; `screen_check.py verify` does not print that hash. The
bundle helper requires its output in the packet's own parent directory and
companions outside the packet. It does not authenticate authors or test
truth. Source verification and evidence verification are distinct.

Before Worker use on the intended Mac, independently test the `/` menu, `/gc`,
`/guided_coding`, an ordinary fresh conversation, target selection, adopted
and unadopted projects, a case-insensitive path probe, and the relevant
project instruction loading (especially an `AGENTS.md` project). The
central source and its complete release packet need independent review and
explicit adoption. No live Mac, server, or Claude Code run is claimed here
(SUPERSEDED 2026-09-29 (rev12; see the Rev12 section below): the Mac CLI and Desktop runs are reported
there; no server run is claimed).

The earlier r09 documentation review and r10 rev3–rev10 review findings
informed this candidate. They do not approve subsequent revisions. Rev13 had 24 regular source
files including its manifest. Rev14 adds the setup workflow, general defaults
and method template, while retaining a single manual skill and no per-project
installer.

## Rev12: `/skills` guidance from the rev11 Mac trial

Rev12 changes documentation and the release label only; no tool, test or
contract byte changed. On 2026-09-29 the rev11 source was trialled on the
intended Mac (macOS, case-insensitive APFS, Claude Code 2.1.284) with an
isolated test configuration and, briefly, the normal configuration for
Desktop. The 42 tests passed with no skip. Typed `/gc` and `/guided_coding`
both ran the skill from the same link in the CLI and in Desktop; ordinary
prompts stayed ordinary; target selection and `AGENTS.md`/`CLAUDE.md`
loading matched Claude Code's documented defaults. Two observations
contradicted the rev11 guide: Desktop's `/skills` listed **no** personal
skill, even after `/reload-skills`, and Desktop's `/` menu listed `gc`
while every transcript (CLI and Desktop) recorded the command as
`/guided_coding`. The guide now tells users to confirm with the `/` menu or
`/gc help`, not `/skills`. Not observed: the interactive CLI `/` menu,
Desktop `/gc <task>` in an unadopted project, Windows. In
`cctbx_project`, `libtbx/tst_guided_coding.py` runs this package's tests
from the shared test suite.

## Rev13: installer-precompiled bytecode

An installer-built PHENIX runs `libtbx.py_compile_all -i` over its modules,
which calls `compileall.compile_dir` and writes
`__pycache__/<name>.<tag>.pyc` beside every `.py` file. Rev12's
complete-source check refused those files, so `/gc` would stop on such an
installation. Rev13's `verify-source` accepts exactly that pattern, and only
beside a **listed** module, and prints a NOTE with the count. Every other
unlisted file is still refused, including bytecode for an unlisted module, a
sourceless `.pyc`, other files in `__pycache__`, bytecode for a listed
name in another directory, and unchecked-hash bytecode (which Python runs
without consulting the source; compileall never writes it by default).
**Limit:** accepted bytecode is checked by name and header type, not by
content. The tools never load it, but an ordinary import does: running the
package tests directly on an installation imports through that cache. The
shared `libtbx/tst_guided_coding.py` therefore runs the package tests on a
temporary copy of the listed files only. The tools never load that bytecode:
`screen_check.py` runs as a script and imports only the standard library,
and `review_bundle.py` now executes `screen_check.py`'s source text instead
of using the caching import loader. Tests: a precompiled synthetic source
passes; seven unlisted-code variants are refused; a crafted bytecode file with
a valid header is loaded by an ordinary cached import (positive control) but
not by either tool. The same three tests fail against the rev12 tools. On
the Mac the real `libtbx.py_compile_all -i` was run on copies: rev12 refused,
rev13 verified with the NOTE.

## Rev14: recoverable, relevant setup

The new instructions live in `payload/SETUP.md` with general proposed
`SETUP_DEFAULTS.md` and an optional `templates/PROJECT_METHOD.md` structure.
SKILL.md routes explicit setup and incomplete/stale task setup to them;
GUIDE.md and WORKER.md route later gaps. Registration prints a `/gc setup`
next step without creating project settings. The user guide covers first
registration, existing projects, optional domain profiles and migration.

Required behavior is evaluated with fictional projects: local-only first
use; changed environment with deferred remote setup; csh/PATH recovery with
hooks retained; conflicting legacy settings/custom authority; a narrowly
authorized saved edit; and an optional domain profile from another owner.
Check that relevant values and their sources survive, current authority and
unrelated settings remain, only necessary questions are asked, and no
unrequested jobs/connections occur. These are conversational checks, not a
claim that prose instructions mechanically enforce each condition.

The registration cue has an executable regression assertion alongside the
existing no-overwrite/no-adoption registration checks. Run the suite from a
clean copy of the manifest-listed source, only after the whole source-check
block above succeeds. This is a separate optional test step, not a command
to append after a failed verification block:

```bash
python3 -I -B -m unittest discover -s tests -p 'tst_*.py' -v
```

The accompanying development evidence records commands, outputs, source
identity, scenario prompts and observations, authorship, and remaining
limits. A fresh Claude Code CLI/Desktop trial is still required to establish
actual skill invocation, instruction loading, and the interactive flow on
the intended machines. No live remote host, real account setup, native
Windows environment, source installation or publication is claimed here.
An outside release reading and the Developer's integration/rollout decisions
remain separate from author-side validation.

## Rev16: isolated stabilization and validation from rev14

Rev16 retains the Rev14 workflow described above. No claim tool, new
installation-file replacement/restoration engine, parallel-job support,
second PHENIX installation, broader assumption-dependent implementation
permission or contract exception is added. WORKER.md's existing
save/apply/test/restore guidance remains; this build does not exercise it
on a live installation or certify arbitrary restoration behavior.

Run the same standalone unittest command above on the final, manifest-listed
source. The original repository wrapper is not supplied in this handoff.
Keep temporary tests under a symlink-free resolved TMPDIR; record platform
skips separately from failures. The source manifest establishes complete
bytes, not correct setup behavior or instruction loading.

The accompanying builder-b evidence supplies requirements-first criteria,
fixture inputs, exact helpers/commands, source identities, file states,
saved procedural responses, raw relevant logs and actual role attribution.
A fake version executable tests registration only. A procedural simulation
using candidate instructions is not an actual /gc invocation. Neither a
PATH version, registration symlink, matching text nor a simulation establishes
Claude Code's actual instruction loading. Observe that in a fresh disposable
intended-client configuration/project with this exact candidate, using the
packet's client-validation handoff, before proposing activation.

The supplied Stage A method and optional PHENIX profile snapshots remain
unchanged context, outside generic source. They are not fresh observations
or permission to replace live settings. Stage B remains deferred and must
not be saved or activated for rev16. The user-reported 2026-10-01 disk cleanup
is context, not a measurement or an edit to those historical snapshots.

Outside review must read the final source and bound evidence; any old
reading belongs only to its old candidate. Prepare materials for the
Developer to deliver, then keep acceptance, integration, activation and
publication separate from successful internal checks.
