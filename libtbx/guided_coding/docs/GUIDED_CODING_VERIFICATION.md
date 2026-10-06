# GuidedCoding verification and limits

This file separates three kinds of evidence: what the supplied source
implements, what earlier execution records report, and what remains
unestablished. A successful checksum, test suite or client-loading check
is not approval to integrate or publish.

## Status and identities

The published implementation baseline is r10 rev16
`enumcheck-20261002T194333Z`, published on 2026-10-02 as the limited opt-in
pilot. Its recorded identities are:

| Item | Identity |
| --- | --- |
| Published commit | `c36887c7f489018af4f91773246ecde32d1b4e24` |
| Source archive SHA-256 | `fce3627e46982f3b5fad1baaad7ada538984e9087a84b2147548b4564e844968` |
| Source manifest SHA-256 | `b86d0490df59f854c214e6771ad8bb44f03069a38612be198edff4ebb11238e0` |
| Contract version / SHA-256 | `2026-09-17` / `ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b` |
| Publication review bundle SHA-256 | `0e4c112d3726026bd397a4e50ac2a3cff79b146883f8ceab269e3a069bc2c6ca` |
| Publication packet identity | `9834e7a5bfdf6b8dd8401972265715f527432ba277b63c4ee5c50f130508f44b` |

The `docs-20261003` revision (published 2026-10-03 as commit
`b0747a4a55f29db3abe04358480d5867e94cb792`, source manifest SHA-256
`5edbec14fe45d3538da533fd995a5d59ddaa47df7629daf27321ab2a1134f1be`)
changed documentation only and kept the literal label `r10 rev16 candidate`.

This revision, r10 rev17 candidate `followups-20261004`, changes
`screen_check.py` (the reading's `Scope:` line, `OUTGOING.txt` bindings and
the `NOT RUN` waiver form), adds `publication_precheck.py` and
`records_history.py` with their tests, and revises the skill entry,
Worker, Guide, Setup, defaults, method template, screens, reviewer brief,
transport, Helper and Roles texts. The contract, the adoption rule and the
SETUP §5 recipe are unchanged. Its commit, archive, manifest, bundle and
packet identities are recorded in its task and publication records
(`phenix/.claude/records/2026-10-04-gc-followups-A/`), not in this file,
which cannot carry its own hash. Earlier reviews and test runs apply to
their recorded source; they do not approve this diff. The literal
release-family check now expects `r10 rev17 candidate`; it does not
establish an exact revision or approval.

A later bounded change (2026-10-06, record
`phenix/.claude/records/2026-10-06-gc-app-version-check/`) made the version
check app-aware (`claude_version`, `check_claude_version`,
`check_app_engine_version` and `APP_ENTRYPOINT` in `screen_check.py`),
rewrote `ClaudeVersionChecks`, and revised the matching sentences in the
skill entry, user guide, architecture table and this file. The label,
contract, adoption rule and setup recipe are unchanged.

The October 2–3 publication, client and installation facts below come from
the supplied controller and outside-review records. The documentation pass
checked the uploaded baseline against its complete manifest and read the
implementation and tests; it did not rerun those sessions, inspect live
remote refs or act as an Outside Reviewer. Raw historical evidence is kept
in the associated task/review records, not duplicated in the generic package.

## Release and source checks

First compare an archive with the independently supplied SHA-256 in its
review/publication handoff. From a trusted extraction, replace the path
below with the canonical source root and run this complete Bash block:

```bash
cd /absolute/path/to/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev17 candidate' payload/RELEASE
```

Every dependency uses `&&`. A failed directory change or listed checksum
must stop before executing the checker; a failed inventory check must stop
the later commands. Do not append unguarded preparation or registration.
After a successful block, separately compare the manifest's SHA-256 and
full release label with the expected handoff identity. Repeat the guarded
source check before central tools and consequential work as the skill
requires.

`verify-source` rejects missing/changed listed bytes, extra files, symlinks
and multiply linked files. The accepted installer-bytecode pattern is an
exception, not authentication of bytecode contents. `-I` isolates imports
from the current directory and `PYTHONPATH`; `-B` suppresses writing caches.
The tool invocation uses source rather than package caches. These controls
establish consistency, not supplier identity or model compliance.

## Trace a claim to the source

For example: **registration refuses an existing destination instead of
creating a nested link.**

1. Open [`register_skill` in `screen_check.py`](../payload/tools/screen_check.py).
   It computes the single `skills/guided_coding` destination, checks
   `exists()` or `is_symlink()`, and calls `fail` if occupied. The creation
   is a direct `os.symlink`, not a recursive `ln -s` command.
2. Read `SkillRegistrationChecks` in
   [`tst_screen_check_v20.py`](../tests/tst_screen_check_v20.py). The fresh/
   second-run case and the other-link/directory/dangling-link case assert
   refusal and absence of a nested link. Other cases cover version and
   parent-path refusals.
3. To establish that those tests ran for a particular revision, inspect
   its captured command, source identity, stdout/stderr and exit status.
   Reading the test source alone establishes its intended assertion,
   not an observed pass.

Other important mappings:

| Claim | Implementation / test source | What it does not establish |
| --- | --- | --- |
| Complete-source consistency | `verify_source`, `regular_files`; `SourceInventoryChecks` | Origin, semantic correctness or authenticated cache contents |
| Claude Code minimum 2.1.281: the PATH `claude` in a Terminal session; the app engine named by `CLAUDE_CODE_EXECPATH` when `CLAUDE_CODE_ENTRYPOINT=claude-desktop`, NOT CHECKED when that version cannot be read | `MIN_CLAUDE_VERSION`, `APP_ENTRYPOINT`, `claude_version`, `check_claude_version`, `check_app_engine_version`; `ClaudeVersionChecks` | That the app marker stays stable, actual instruction loading, or trustworthy executables |
| Screen form and result identity | `blocks`, `developer_view`, `code_identity`, `check`; `ScreenChecks` | Correct criterion, true logs, complete tests or Developer authorization |
| Matching outside reading | `reading_matches`; conditional/waiver/pending controls | Independence or authenticity of the reviewer, satisfaction of a semantic condition |
| Two-domain review transport | `make_bundle` in `review_bundle.py`; `ReviewBundleChecks` | Privacy, evidence completeness or permission to publish |
| Reading scope (exactly one `Scope:` line, counted before its value is checked), outgoing bindings (no repeated key within a block), waiver form | `reading_matches`, `parse_outgoing`, `suite_waiver` in `screen_check.py`; the scope, outgoing and waiver cases in `ScreenChecks` | Git state, the reviewer's independence, or that the quoted words authorize anything |
| Pre-push destinations, settings, bindings, dry run (every Git query's exit checked before its output; a failed configuration or identity query refuses, while Git's no-match answer to the configuration query is accepted), command shape (whitelist of the documented form) | `publication_precheck.py`; `tst_publication_precheck_v20.py` on temporary repositories and with a fake Git on PATH | A decision to push, the remote's later state, settings outside the repository's Git config scopes, any Git command not submitted to the validator, or false acceptances other than the demonstrated ones |
| Read-only task history | `records_history.py`; `tst_records_history_v20.py` | Truth of a record's contents, or any grant, reservation or scheduling |
| Checked setup ordering | `SKILL.md`, `SETUP.md` §5 and its embedded recipe | Mechanical enforcement by Claude Code or reliable obedience by every Guide |
| Test isolation in the repository | [`libtbx/tst_guided_coding.py`](../../tst_guided_coding.py) | Complete source verification of arbitrary unlisted files in the original installation |

There are **118 test methods** in the four shipped package test files
(63 for the checker, 13 for the review bundle, 27 for the pre-push
checks, 15 for the history listing). Subtests exercise additional cases;
a count of 118 is not 118 independent behavioral claims. The checker,
pre-push and history tests of the rev17 release were written by two separate
test-writer sessions from the interface specification, without reading
the implementation; the baseline control (`observations/baseline_control.txt`
in the release's evidence packet: the final checker test file run against
the previous published checker) shows which of them fail there; the
pre-push and history tools have no earlier version to compare with. The
version-check tests of the 2026-10-06 change (`ClaudeVersionChecks`, the
app registration case and `VersionGateScopeChecks`) were written by a Worker
subagent from the Developer's case list in the same session as the
implementation, which it was allowed to read for exact strings; the control
run of that test file against the previous checker is in that change's record
(`evidence-work/base_control.txt`). The setup
shell recipe and native LLM scenarios are not covered by unit methods;
their controls are in the task and review records.

## Recorded checks of the published implementation

| Evidence | Recorded result and reach |
| --- | --- |
| Enumeration correction and follow-up reading | The two `find` producers' failures, including partial listings, stop verification. A/B/B2 did not print VERIFY_OK; the tested driver held apply; protected files were unchanged. Seven earlier recipe controls retained their results. R-F1 was closed. |
| Package and repository checks | Controller reported 45 package tests passing on the Mac. The publication reviewer ran the repository harness on Linux: unit success with one case-insensitive-filesystem skip and a libtbx precompile skip. Skips are not passes of those branches. |
| C1, disposable CLI configuration | Exact enumcheck candidate, CLI 2.1.284: help and status, source guard and command-to-expansion comparison. Reach is CLI loading/help/status, not setup/recovery. |
| 4b, live CLI | Help loaded the expected source and ran the guard. The successful attempt followed an expired-auth failure and user sign-in. |
| D1, Desktop Code tab | Help loaded the personal central link and ran the source guard; transcript reported runtime 2.1.286. Reach is **loading/help only**, not Desktop status or full setup/recovery. |
| Publication r2 | Outside reading: PROCEED for the limited pilot; P-F1/P-F2/P-F3/P-N1 closed, R-F1 retained closed. A clerical description of three test assertions was corrected externally. Genuine publication screen passed; Developer PUBLISH and verified remote pushes were recorded. |
| Mac and shared Linux installation synchronization | Controller reported clean checkouts at the publication commit on the Mac and on one shared Linux installation, source equality between them and repository harness success on both. Synchronized GC files on that installation are not native Claude Code client validation: no Claude client was installed or run there. This was not a PHENIX full-suite run. |
| Tag removal, October 3 | Controller reported the pilot tag absent remotely, its local copy removed after preservation, and the published commit still in upstream history. Repository tags are reserved for releases; use the commit/manifest, not the removed tag. |
| CI for the two published commits (read 2026-10-04) | GitHub Actions `quick` and `clutter_and_syntax` completed with conclusion success for c36887c7 (2026-10-02) and b0747a4a (2026-10-03); the mirror workflow succeeded after each. Read through the public API without authentication; a read, not a rerun, and not a PHENIX suite. |
| Native checks of this candidate, Mac CLI 2.1.284, `claude -p … --permission-mode auto`, candidate loaded through a project-level skill link in a scratch project (2026-10-04) | Four sessions, each printing the guard result (VERIFIED complete source, `r10 rev17 candidate`) for the worktree path: `status` recovered adoption, method and records from files alone; `history` listed the one record read-only and created no file; a setup save with a truncated approved record stopped at verify (FAIL: approved record does not match) and held apply; a setup save with no prior record or ABSENT marker stopped at verify and held apply; the protected files' hashes were unchanged in every session. The CLI title instruction (`/rename GuidedCoding: <task>`) was printed in three of the four sessions and omitted in one. Reach: one client, one mode, scratch targets; the sessions used the Developer's live login and wrote their own client bookkeeping outside the targets. |
| Desktop session rename, Guide's own session (2026-10-04) | In the Guide's own Desktop session (which had loaded the published rev16 skill, not this candidate), the client's rename tool replaced an app-generated title without a prompt: the tool's availability, not this candidate's behaviour. |
| Desktop session titles, pinned snapshot (2026-10-04) | Two Desktop sessions in a scratch project whose project-level skill directory was a byte-identical, read-only copy of an earlier revision of this package (candidate 3 of this release). The skill entry `SKILL.md`, whose session-title instructions these observations exercised, is byte-identical between that revision and the published one; the tools, tests and some texts changed in later revisions, so the observations are of the title behaviour only, not of the package as a whole. Each transcript shows the skill expansion and a passing guard on that snapshot. (a) An app-generated title was replaced by `GuidedCoding: help` through the client's rename tool without a prompt. (b) A title the Developer had set by hand (`KEEP MY TITLE`) was left unchanged: the skill read the session's title, judged it user-set and did not call the rename tool, so the app's approval dialog was not exercised. (d) The Developer reported both titles unchanged after closing the windows. Four earlier Desktop sessions are preserved in the task record and not counted: two met a failed guard because the Guide was editing the worktree, one renamed itself on the previous candidate, one met a failed guard and left a hand-set title untouched. The app's "isn't a command here" notice for a project-level skill appeared each time; it is not evidence either way. Reach: one Desktop client, one snapshot, two sessions. |
| PHENIX test discovery (A7, 2026-10-04) | `phenix.find_program search_type=tests search_text=<function> tests.search_tests_by=function_called` traced `run_autobuild` to its calling tests; the default mode matched test names; a function newer than the static index (dated 2026-05-06) produced no entry and no message; `git_affected_tests=True` saw only uncommitted modifications in the three module directories. Project guidance, not a package feature. |
| App-session version check (2026-10-06, record `phenix/.claude/records/2026-10-06-gc-app-version-check/`) | In the Guide's own Claude app (Code tab) session, with the separately installed CLI removed from PATH for each command only: the candidate checker printed `VERIFIED Claude app's Claude Code engine 2.1.288 …` for the engine named by `CLAUDE_CODE_EXECPATH`; with that variable unset it printed one `NOT CHECKED` line and exited 0; with `CLAUDE_CODE_ENTRYPOINT` unset and no `claude` on PATH it failed as a Terminal session; with the marker unset and the CLI on PATH it verified the PATH `claude` (2.1.284). One real Terminal observation: the PATH CLI 2.1.284, run once non-interactively (`claude -p`, scrubbed environment, a one-off session-start hook passed with `--settings`), reported `CLAUDE_CODE_ENTRYPOINT=sdk-cli`, not `claude-desktop`, so it takes the Terminal path; an interactive terminal session was not observed. Reach: one Mac, one app session, one non-interactive CLI run; nothing uninstalled, no settings changed. |

The full PHENIX server suite and a baseline/candidate roster were **NOT RUN /
NOT PRODUCED for the GC-only pilot publication**, under the Developer's
specific waiver. The publication record states this publicly. Local tests
do not replace that suite, and the waiver does not apply automatically to
a later change.

## Remaining qualifications

- **Native damaged-record branch:** two non-interactive CLI sessions of this
  candidate held apply on a damaged and on a missing recovery record (table
  above). That is one client and one mode on scratch targets; the
  driver/Guide still determines what gets invoked next in any other session.
- **Exact proposal and recovery history:** earlier cases omitted complete
  proposals or approved recovery text. Assisted and later native cases
  showed improvement; they do not erase those failures. A relocated setup
  record was disclosed but not wholly literal to its displayed proposal.
- **Fresh-session context:** disposable auto-memory was written outside the
  target during earlier trials. Later sessions using that configuration
  cannot be called memory-free or dependent only on target records without
  checking those inputs. No broad fresh-recovery guarantee follows.
- **A1 provenance:** an untouched pre-redaction debug log was not located
  by the reported bounded search. Other transcripts and partial excerpts
  remain; the reviewer found no material behavioral claim unsupported solely
  by that loss. The complete original and “only alteration” claim remain
  unverified. This is a provenance limit, not proof of another behavioral
  failure.
- **Client/platform scope:** a passing Terminal version check says nothing
  about the app engine, and the app branch rests on the
  `CLAUDE_CODE_ENTRYPOINT=claude-desktop` marker and `CLAUDE_CODE_EXECPATH`
  as observed on 2026-10-06, whose future stability is not established.
  Earlier TUI observations had continuation/permission-mode
  qualifications. Native Windows, real csh and all-client setup/recovery
  behavior are not established by the final-candidate observations above.
  The repository wrapper explicitly limits the tools to macOS and Linux;
  its Windows skip is not evidence of Windows support.
- **Isolation:** observed file equality and absence of recorded tool calls
  have bounded reach. They do not establish whole-machine immutability,
  all syscalls or unchanged live authentication state.
- **General effectiveness:** this package does not establish a failure
  rate, superiority to unstructured LLM use, or reliability for every future
  Guide. Separate test authorship is recorded, not proof of independence.

## How to check a later documentation or source revision

Keep the previous evidence unchanged and bind the new source separately.
For a docs-only revision, check the diff against the implementation,
command syntax and preserved authority boundaries. Record unchanged
executable/contract bytes explicitly. If instructions or claims change,
assess those changes; a unit suite cannot establish their behavioral effect.

After the source guard succeeds, package tests can be run from a trusted
clean, manifest-listed copy with a resolved, symlink-free temporary path:

```bash
python3 -I -B -m unittest discover -s tests -p 'tst_*.py' -v
```

From an appropriate cctbx environment, the repository wrapper is:

```bash
libtbx.python /absolute/path/to/cctbx_project/libtbx/tst_guided_coding.py
```

These are separate authorized check steps, not commands to continue after
a failed source guard. The wrapper copies listed files into scratch, uses
a resolved temporary directory, runs the package tests, and checks that
listed installed bytes did not change. When libtbx precompilation is
available it exercises that branch on another copy. It skips Windows and
Python 2; it can also skip the real precompile branch. Report each skip.
The wrapper reports existing unlisted files rather than using them for the
test copy, so retain the complete-source check as a separate requirement.

For native validation, identify the exact source, client/runtime, selected
skill path, loaded project instructions and persistent context. Inspect
recorded command-to-expansion linkage, actual operations and resulting
files. A version string, registration link, model-reported path or procedural
simulation alone is insufficient. Matching decoded transcript text is not
proof of file-byte loading. Select observations for the changed behavior
and intended support scope; do not label an earlier candidate's session as
a fresh run of the new one.

The Outside Reviewer evaluates the new source and bound evidence; the
Developer's acceptance, integration, activation and publication decisions
remain distinct. A genuine final screen must cite a qualifying reading of
its own frozen packet. Use the existing contract/project rules for any
waiver; do not infer one from this document.

## Relevant history

Rev10–11 corrected exclusive registration and configuration-path handling.
Rev12 corrected user guidance from the Mac trial. Rev13 added the narrow
installer-bytecode exception and source-only tool loading. Rev14 introduced
recoverable, relevant setup. Failed rev15 parallel-job work was excluded
from the rev16 lineage.

Rev16 stabilized rev14, then corrected fail-fast checksum instructions,
approved-text preservation, explicit comparison failure propagation and
checked file enumeration. A failed `check && echo OK` inside a script was
not a sufficient stopping rule; unchecked producer exits were also not
sufficient. The final published implementation is enumcheck, not an earlier
F1/checked-save/failprop candidate. Earlier reviews belong to those earlier
identities; frozen packets and errata retain that history.

Rev17 (`followups-20261004`) followed three corrections recorded in the
PHENIX project during October 2026: a development tag pushed to a
repository that reserves tags for releases (removed by its maintainer on
2026-10-03), a publication that proceeded on an inferred suite waiver
(corrected by the Developer on 2026-10-04), and installation updates that
were planned after the push rather than with it. It adds no parallel-job
coordination; that remains a separate planned release.

Current client reference material: [skills](https://code.claude.com/docs/en/skills)
and [instruction loading](https://code.claude.com/docs/en/memory).
Vendor documentation describes the client; it is not a GC native test run.
