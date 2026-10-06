# GuidedCoding architecture

GuidedCoding has three distinct parts: centrally stored instructions and
tools, project-owned settings and evidence, and the Claude Code client
that loads the skill and performs operations. A check in one part does
not establish the behavior of the others.

## Source, client and target

The central source is normally `cctbx_project/libtbx/guided_coding/`.
A personal symlink at `~/.claude/skills/guided_coding` points to that
root. With `CLAUDE_CONFIG_DIR`, the same relative registration lives under
that configuration directory. `SKILL.md` declares `name: gc` and
`disable-model-invocation: true`. The recommended user spelling is
`/guided_coding`; `/gc` remains its short name.

The skill selects Claude Code's primary working project unless the user
explicitly names a different target. Reading a central file does not
change the target to the procedure repository. Each task records which
central source identity governs it.

| Location | Responsibility |
| --- | --- |
| `SKILL.md` | Entry, source/target resolution, control requests, adoption and readiness instructions |
| `SOURCE_MANIFEST.sha256` and `payload/RELEASE` | Listed source identities and human-readable revision label |
| `docs/` | User-facing description, implementation map and validation limits |
| `payload/DEVELOPER_GUIDE_CONTRACT.md`, `ROLES.md` | Reserved Developer decisions and role boundaries |
| `payload/GUIDE.md`, `WORKER.md` | Governing task workflow and required evidence |
| `payload/SETUP.md`, `SETUP_DEFAULTS.md` | Relevant setup, recovery and proposed general defaults |
| `payload/screens/`, `templates/` | Decision screens, report/handoff formats and optional method structure |
| `payload/tools/` | Executable source, evidence, screen and bundle checks |
| Target authority and project method | Explicit adoption plus environment, work locations, commands and restrictions |
| Target `.claude/records/` | Durable task/setup proposals, prior files, decisions and evidence |
| Client configuration, transcripts and memory | Client-managed state; not confined by GC to the target records |

The project may already have a method in another location. Its authority
names the canonical path; the template does not impose a second copy.
Optional domain cards supply candidate data with provenance. They are not
global instructions, project adoption or transferable permission.

## What is executable and what is instructed

| Capability | Implementation | Boundary |
| --- | --- | --- |
| Check complete source | `screen_check.py verify-source SOURCE` | Checks inventory and bytes against the supplied manifest, with a narrow bytecode exception; not supplier authentication |
| Check Claude Code minimum | `check-claude-version` | Parses `--version` of the PATH `claude` in a Terminal session, or of the app engine named by `CLAUDE_CODE_EXECPATH` when `CLAUDE_CODE_ENTRYPOINT=claude-desktop` (the marker observed on 2026-10-06; NOT CHECKED, not a failure, when that version cannot be read); not instruction loading, and not the marker's future stability |
| Register personal skill | `register-skill SOURCE` | Exclusive link creation after source/CLI checks; does not adopt or set up a project |
| Help, status, setup, uninstall and task entry | Branches in `SKILL.md` | Instructions followed by the model; no executable dispatcher for these requests |
| Preserve/verify/apply setup | Instructions and Bash example in `SETUP.md` §5 | Guide must construct correct records and obey ordering; no client-enforced transaction |
| Freeze and verify evidence | `freeze DIR`, `verify DIR` | Binds file contents and inventory; does not establish truthful logs or complete evidence |
| Check/present a decision | `check KIND SCREEN`, `present KIND SCREEN` | Checks required form and specified evidence bindings; not a review or a Developer decision |
| Package outside review | `review_bundle.py PACKET COMPANIONS OUTPUT` | Builds two integrity domains with restricted layout; not privacy review or approval authentication |
| Bind a publication proposal | `check publication` reads `OUTGOING.txt`, the reading's `Scope:` line and the `NOT RUN` waiver form | Form and identity binding only; not Git state and not authorization |
| Pre-push Git checks | `publication_precheck.py check REPO OUTGOING [--repository NAME] [--fetch] [--dry-run]`; `vet -- <command>` | Comparison of the live repository with the record; a dry run counts only when Git exits successfully and reports exactly one fast-forward of the approved branch; `vet` accepts only the documented command shape and refuses every other word (Git abbreviates options and reads settings from the environment, so a list of forbidden words cannot be complete). The validator checks the command text submitted to it; it does not intercept Git operations run any other way. `--fetch` and `--dry-run` contact the remote, nothing else does; not a push, not a decision |
| Task history | `records_history.py RECORDS_DIR [--limit N]` | Read-only listing from `JOB_SUMMARY.txt` or `RECORD.md`; grants nothing and runs nothing found in a record |

`GC_PAYLOAD_ROOT` is a shell convention used to address central tools and
payload files. Neither Python tool reads it as a configuration variable.
The explicit command argument determines which source or evidence root
is checked.

## Authority and setup

The contract version is `2026-09-17`, SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`.
Before Guide/Worker instructions govern a target, its current authority
must explicitly adopt that identity and identify its own method. A skill
invocation, link or old release label is insufficient. Control requests
are handled before the task adoption gate; setup can propose adoption.

The Developer retains meaning, value, risk, waiver, acceptance and
consequential-action decisions. A Guide scopes work; Workers perform
bounded tasks; an Outside Reviewer evaluates the supplied evidence at a
gate; the optional Helper explains it for the Developer. A subagent with
its own context is not thereby independent of the Guide who briefs it.
Records distinguish who chose a criterion, wrote a test and ran it.

Setup recovers known settings before asking for missing values. It checks
readiness only for the requested operation. Stale verification does not
erase a known setting. Default setup does not connect to a server or run
a build/test suite. Optional remote gaps need not block independent
local work authorized by the project.

The save sequence requires complete proposals and lossless project recovery
records before protected method/authority changes. Verification compares
approved records with separately regenerated expected text and checks
prior state. The Bash example checks each `find` producer's exit status
before consuming either NUL-separated list. An incomplete listing with
an error cannot count as successful enumeration in that example.

The Guide issues apply only after verify succeeds. Apply guards the
requested files against the saved prior state before its first target
write, then writes and reads back individual files. It does not create
an atomic multi-file update. Failures after writes can leave PARTIAL state;
recovery records must remain available. The recipe has separate `verify`
and `apply` invocations, not a client-enforced verified-state token.

Instruction selection belongs to the client. Setup must preserve the
project's actually loaded authority, including relevant `AGENTS.md`
settings, rather than blindly adding a competing instruction file.
Client auto-memory and transcripts may persist outside the target, so
fresh-session recovery claims require observation of that context.

## Source verification and bytecode

Use the complete guarded block in
[Verification](GUIDED_CODING_VERIFICATION.md#release-and-source-checks).
The trusted checksum command must succeed before the package's checker
executes. The checker then verifies exact listed bytes and complete
inventory, rejecting symlinks, multiply linked files and unlisted files
except the recognized precompile pattern.

Permitted bytecode has a listed `.py` counterpart in the same directory,
a recognized `__pycache__` name and an accepted header flag. Its contents
are not authenticated. `screen_check.py` runs as source; `review_bundle.py`
explicitly loads the sibling checker source text. Documented calls use
`python3 -I -B`. Ordinary test imports can behave differently, which is
why the repository wrapper copies only manifest-listed files before testing.

A source manifest detects changes relative to itself. Its separately
recorded hash binds that manifest to a review or task. Neither mechanism
proves who supplied it. Recheck the source before consequential work;
if it changes, resolve the changed procedure before dependent work continues.

## Evidence and decision screens

Source and task evidence use different manifests:

| Manifest | Contents and identity |
| --- | --- |
| `SOURCE_MANIFEST.sha256` | All listed package files except itself, using `./relative/path`; source identity is this manifest's SHA-256 |
| `MANIFEST.sha256` | All frozen evidence files except itself, using relative names without the source manifest's `./` convention; packet identity is this manifest's SHA-256 |
| Bundle `SHA256SUMS` | Inner packet archive and external companions; the outer archive's hash is supplied separately |

`verify` prints `VERIFIED evidence`, not the packet hash. Compute the
manifest hash separately when recording identity. `freeze` refuses to
replace an existing manifest; a changed packet needs a fresh freeze.

The screen checker requires at most 28 lines and 220 words, specified
headings, no unresolved template placeholders, and the matching final
ACTION. Plan/result screens also need why it matters, a recommendation
with a reason, and the next action.

For a result, the checker checks code-identity fields, a passing recorded
outcome, a matching tested-tree reference, and a report containing the
complete `CHANGE.diff` plus exact-changes/new-tests headings. It does not
execute that check or prove that the reported test authorship is true.

For a full result or publication it binds the reading to the exact packet
and reading-file hash and requires a proceeding verdict. A conditional
verdict needs a matching disposition. WAIVED and PENDING conditions remain
visible; quoted authorization is a record to evaluate, not authentication
of its author. A light result explicitly defers the outside reading to
publication.

The reading carries exactly one `Scope:` line. RESULT accepts
`integration` or `integration and publication`; PUBLICATION accepts
`publication` or `integration and publication`, so one reading can serve
both screens of one frozen packet. A publication packet must contain
`OUTGOING.txt`, naming for each repository the remote, its URL, the base,
the commit, its tree and the `<commit>:refs/heads/<branch>` refspec, and
every outgoing commit must appear in the screen's BATCH section. A suite
file whose first line starts with `SERVER_SUITE: NOT RUN` must carry a
`Waiver (Developer ...):` line followed by `> ` quoted lines. These are
form checks: they bind what is claimed to the packet. Git state is
compared separately by `publication_precheck.py` (which contacts the
remote only with `--fetch` or `--dry-run`), and
authorization remains the Developer's quoted decision.

For publication, the tool additionally requires nonempty
`SERVER_SUITE.txt` and `ROSTER_COMPARISON.txt`. It does not parse their
results or apply the result-specific code-identity checks. If a required
suite is explicitly waived, those files must truthfully record NOT RUN /
NOT PRODUCED and the waiver, with its justification reviewed separately.
File presence alone cannot turn that into a passing suite.

`present` first checks the screen, then formats the Developer view and
hides internal identity strings. The frozen evidence retains the identities.
No screen command integrates code or pushes a repository.

## Review transport

The [transport instructions](../payload/REVIEW_TRANSPORT.md) and
[reviewer brief](../payload/OUTSIDE_REVIEWER_BRIEF.md) govern the handoff.
The frozen packet includes the evidence and proposed report. Separate
companions contain `APPROVALS.md`, `CHECKER.md`, and optional `ERRATA.md`.
The bundle helper rejects nested or duplicated reserved companions and
requires `README_FIRST.md` and `PROOF_SUMMARY.md` in the packet.

Output must be in the packet's own parent directory, not merely somewhere
outside it. The helper uses directory identity and open handles to avoid
specified destination redirections. It accepts an alias to that parent;
it does not promise safety under arbitrary hostile concurrent movement.
The packet and output should not be modified during packaging. Deliver
only a successfully completed and verified bundle.

The cover identifies the bundle and packet separately. A reviewer reply
and later dispositions stay outside the frozen packet. A new packet needs
a reading for its new identity; an earlier HOLD or a synthetic control
cannot supply a genuine PROCEED gate. One reading scoped to both gates
may serve both screens of one packet; a re-staged proposal (new commit,
base or packet) is a new identity and needs a reassessment of the
affected scope bound to it, which may be a short addendum.

## What this version does not add

There is no per-project procedure installer, automatic profile creation,
account hook, `/wrap` command, `gc_claim.py`, parallel-job coordinator or
new installation-restoration engine. The Worker's existing candidate
save/apply/test/restore instructions remain. Worktree separation alone
does not establish separate imported installations or shared-resource safety.
The history listing writes nothing and the pre-push checks change
nothing in the repository (they contact the remote only with `--fetch`
or `--dry-run`): neither
reserves a resource, schedules work, pushes, or updates an installation.

The personal link exposes the current central checkout, not a pinned
release. Updating that checkout therefore needs coordination with active
work and the repository's integration rules. Keep project methods,
existing records and permissions separate from procedure updates. The
published pilot and later documentation changes have distinct source
identities and publication decisions.
