# Guided Coding registration and command reference

Use this reference for exact registration commands, control requests and source checks. For installation and a first task, start with the [User Guide](GUIDED_CODING_USER_GUIDE.md). The [verification record](GUIDED_CODING_VERIFICATION.md) explains what the checks establish and where their evidence stops.

## Register once, then choose a project

Keep the procedure in one reviewed, stable directory, normally
`cctbx_project/libtbx/guided_coding/`. Register that source once per local
Claude Code configuration. You do not install another procedure copy in
every repository.

In ordinary Claude Code, you can request registration with this message,
replacing the example with your actual absolute path:

```text
Please set up GuidedCoding from /absolute/path/to/cctbx_project/libtbx/guided_coding for this machine. Read docs/GUIDED_CODING_COMMAND_REFERENCE.md there. Verify the source and release, check the Claude Code minimum version, inspect any existing personal skill, and register the central link only if the checks pass and the destination is unoccupied. Show what changed, then offer /guided_coding setup for the project I choose, recovering saved settings first. Do not connect to servers, change permission settings, or start a coding task.
```

This message asks Claude to perform the setup steps. It needs the source
path to find the kit from another project. If the skill folder was created
after the conversation started, open a new conversation.
Confirm the command with `/guided_coding help`; the `/skills` listing alone
was not a reliable indicator in the recorded Mac trial.

The personal link is `~/.claude/skills/guided_coding`. The directory name
provides `/guided_coding`; the skill's `name: gc` supplies the short name.
A same-named skill can affect selection, so inspect the resolved source if
the response identifies an unexpected revision. See the current
[Claude Code skill documentation](https://code.claude.com/docs/en/skills#how-a-skill-gets-its-command-name).
A cloud session can obtain its own copy from the source repository. The personal
link created here is for this Claude Code configuration; see the User Guide
for cloud and project-specific use.

## Verify a source archive before using it

Compare the archive's SHA-256 with the independently supplied review or
publication record using a trusted checksum tool. The archive and its own
manifest cannot establish their origin by themselves. Inspect its member
names before extraction; the expected source archive contains
`libtbx/guided_coding/`, without an enclosing release-named folder.

Use an empty staging directory. Replace the paths and run this **whole
Bash block**, not separate commands:

```bash
mkdir /path/to/empty-gc-review &&
tar -xzf /path/to/reviewed-guided-coding.tgz -C /path/to/empty-gc-review &&
cd /path/to/empty-gc-review/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev18 candidate' payload/RELEASE
```

Any failed command skips the rest of the block. In particular, a failed
listed checksum must stop **before executing the package's checker**.
Stop dependent work on failure; a successful command appended later would
not repair the failed check. Require all listed hashes, `VERIFIED complete
source`, and the expected release identity from your handoff. The final
`grep` is only a family-label check, not an exact revision check.

The complete-source check rejects missing, changed and extra source files,
links and multiply linked files. Its narrow exception permits certain
installer-generated `__pycache__` files beside listed Python modules. It
reports the ignored count and checks their names and header type, not their
contents. The two tools do not load that package bytecode when invoked as
documented. Ordinary Python imports, including direct unit-test imports,
can load caches; run those tests on a trusted clean copy.

The commands shown here require Bash and `shasum`, plus a suitable Python 3.
They are not csh syntax. Use the appropriate shell explicitly rather than
pasting Bash environment assignments into csh. Project build/test commands
may have a different, recorded shell requirement.

## One-time registration

After the reviewed source is in its authorized central location, inspect
any existing personal skill path. Then run the complete guarded block with
your canonical source path:

```bash
cd /path/to/cctbx_project/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev18 candidate' payload/RELEASE &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py check-claude-version &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py register-skill .
```

`register-skill` rechecks the source and the Claude Code version itself, then creates exactly
one symlink. It refuses an occupied destination, including an existing
correct link, a real directory or a dangling link. A second invocation is
a refusal, not an upgrade. It also refuses a symlinked configuration or
skills directory. Leave unrelated files and links alone; changing an
existing installation needs its own inspected, authorized edit.

If an absolute `CLAUDE_CONFIG_DIR` is set, the destination is
`CLAUDE_CONFIG_DIR/skills/guided_coding`. The helper rejects `..` components
before creating directories. Use a resolved, ordinary source path: the
complete-source checker refuses symlinks in that path, including a `/var`
alias when `/private/var` is the actual path on a Mac.

The version check requires Claude Code **2.1.281 or newer**. What it checks
depends on the session. In a Terminal session it runs `claude --version`
from this shell's PATH and fails when no `claude` is found. In a Claude app
session, recognized by the environment variable
`CLAUDE_CODE_ENTRYPOINT=claude-desktop` (the marker observed in app sessions
on 2026-10-06; any other or missing value is treated as Terminal), it runs
the app's own Claude Code engine, named by `CLAUDE_CODE_EXECPATH`; a
separately installed `claude` command is not consulted. If the app engine's
version cannot be read, the check prints one `NOT CHECKED` line giving the
reason and continues; a version that is read and is below the minimum still
fails. The check does not update the client or authenticate either
executable. Inspect the reported path; if an update is needed, use your
approved client-update method (in the Claude app, update the app) and retry.
Help, status and uninstall remain available when the version check fails.

Registration makes the command available. It does not set up your project,
connect to a server, reserve an installation or run tests. A server with
only the kit files needs Claude Code before it can run the native command.

## Set up the target project

Open a new conversation in the project you want to change and enter:

```text
/guided_coding setup
```

The primary working directory selects the target unless you explicitly
name another project. Loading the procedure from `cctbx_project` does not
make that repository the target.

The [setup workflow](../payload/SETUP.md) first reads the current method
and its referenced settings, records and migration backups. It proposes
recovered values as defaults with their sources, then asks only for gaps
that matter to the intended work. It does not make you remember a command
that is already saved.

| Area | What setup establishes |
| --- | --- |
| Project and instructions | Which project to change, related repositories, instructions Claude loaded, the existing method file, and whose use of the contract those instructions cover |
| Locations | Source, worktree/scratch locations, durable records and outputs; remote locations only when relevant |
| Environment | Shell, PATH or activation, actual executable/import target, build and test commands with working directories |
| Optional remote work | Host, account, access method, installation, logs, load limits and sharing/reservation rule |
| Decisions | Existing authorization, restrictions, verification still needed, and actions requiring your choice |

A value can be recorded, verified, stale, conflicting, unknown, deferred or
not applicable. A known path needing re-verification should not be erased
and replaced with “unknown.” A tool missing from Bash's PATH might already
be configured under csh. Setup inspects the intended environment; it does
not bypass hooks or rewrite global startup files to make a lookup succeed.

The [general defaults](../payload/SETUP_DEFAULTS.md) and optional
[method template](../payload/templates/PROJECT_METHOD.md) are proposals.
They contain no universal server, PHENIX command or account. Optional
remote work may be deferred while independent authorized local work proceeds.
No build, test suite or remote connection is part of default setup.
Additional actions use existing applicable authorization or a specific
Developer decision.

### Defaults cards

You may supply a separate project card, such as a PHENIX defaults card.
Its account names, paths and commands are candidate settings to reconcile
with your current method. Another developer's grants do not transfer.
A personal card and a shareable project card are different artifacts; do
not publish personal settings as general defaults.

For example, a saved shorthand such as `t96` needs its actual definition,
working directory, shell and applicable load rule. Its name alone is not
a test command. Each developer has their own installation path. A saved
hardware count, test-process cap and build-process count are distinct facts.
The generic GC package supplies none of these PHENIX-specific values.
Relevant adopted values and their provenance go into the project's method,
so a later session does not need the original card attachment.

### Contract adoption and loaded instructions

Before guided work, the project’s current instructions must state that it
uses Developer–Guide Contract version **2026-09-17**, SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`,
and identify its project method. The contract is unchanged from the
published pilot. An existing statement naming this exact version and
SHA-256 can suffice. A skill link, release name or command invocation does not record
your agreement to use the contract.

If adoption is missing, setup proposes a small declaration for your
approval. For example, using the project's actual method path:

```markdown
GuidedCoding applies here only after explicit /guided_coding or /gc invocation.
Adopt Developer–Guide Contract version 2026-09-17,
SHA-256 ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b.
Project method: follow .claude/PROJECT_METHOD.md for this project's
build, test, coordination and publication rules.
```

Review that example before using it. Keep the project’s existing
instruction files; do not create another one just to match the example.
In particular, adding `CLAUDE.local.md` can change whether `AGENTS.md` loads
under Claude Code's default selection. Inspect the active setting rather
than assuming both files load; see the
[client instruction-loading documentation](https://code.claude.com/docs/en/memory).

### Checked saving and recovery

The governing details and shell example are in [SETUP §5](../payload/SETUP.md#5-save-a-usable-method-and-a-recoverable-change).
The intended order is:

1. Show complete proposed files, authority edits and setup-record text.
   Identify any allowed mechanical placeholders before authorization.
2. Write lossless approved texts and prior affected files, or explicit
   absence markers, into project recovery records.
3. Verify those records against the complete approved proposal and prior
   state before protected method/authority edits. A failed comparison or
   failed file enumeration must fail the verification step itself.
4. Only after successful verification, apply the authorized edits with
   prior-state guards and read back each saved file.

The records are writes too: the promise concerns verification **before
protected edits**, not before any write whatsoever. A failure during apply
or readback can occur after some files changed. Stop, report **PARTIAL**,
preserve the records and propose the next authorized recovery action.
The recipe does not make multiple file writes atomic or restore them
automatically. The Guide must follow the ordering; the client does not
mechanically enforce it.

Keep the entire setup record literal except for disclosed placeholders.
A meaningful relocation or other unlisted change to the proposal needs to
be accounted for; it should not be described as an exact save merely because
the method and authority files match.

### Existing projects and migrations

Routine startup checks only the settings needed for the requested activity;
it should not repeat a completed interview. A new host, account, checkout
or shell may need selective re-verification.

Before relying on opt-in behavior, inspect always-loaded instructions that
say to read an old `.claude/WORKER.md` for every task. Migrate only the
relevant startup rule under authorization. Preserve or explicitly account
for every setting, restriction, grant and pending item before retiring an
old profile. Do not delete historical records as part of registration.
An unfinished old change must reach its required boundary before a migration
that would alter its governing procedure.

## Commands and decisions

| Command | Intended effect |
| --- | --- |
| `/guided_coding` | Introduce the procedure and wait for a task |
| `/guided_coding help` | Explain commands and opt-in behavior |
| `/guided_coding status` | Read source/link/adoption/method status and relevant setup gaps |
| `/guided_coding history` | List this project's task records read-only: date, job, state, outcome, record and ticket |
| `/guided_coding setup` | Recover settings and prepare authorized method/adoption changes |
| `/guided_coding <task>` | Start bounded guided work after source, adoption and readiness checks |
| `/guided_coding uninstall` | Remove only the personal symlink resolving to this verified source |
| `/guided_coding auto <goal or job list>` | Prepare, test locally and park several jobs one after another, unattended |
| `/guided_coding auto status` | Morning report: one row per job, approval summaries and how to approve |
| `/guided_coding auto resume` | Continue a queue after an interruption, while its grant lasts |
| `/guided_coding stop` | Stop the active auto queue safely |

The same arguments work with `/gc`. These requests load instructions for
Claude; they are not all Python subcommands. Invoke the skill at the start
of each guided conversation, not before every message.

A task presents **APPROVE / REVISE / STOP** for its plan and
**INTEGRATE / REVISE / DISCARD** for the tested result. Publication has its
own **PUBLISH / HOLD** decision. The approval report includes the complete
tested diff and exact new test code, or an explicit “none added” with the
check used. An Outside Reviewer recommends; the Developer decides.
A tool-permission dialog does not replace those decisions.

For a batch that will be published, the frozen packet also records the
exact outgoing commits (`OUTGOING.txt`) and the suite result or the
Developer's quoted waiver for that batch, so that one outside reading can
answer the integration and the publication question for that packet. If
the proposal changes afterwards, it gets a new identity and a
reassessment; the earlier verdict is not inherited. A publication plan
names each installation to update separately; publication alone updates
none. The accepted push is one explicit `<commit>:refs/heads/<branch>`
with tag following disabled, after `publication_precheck.py` has compared
the live repository with the record; no development or recovery tag is
created.

Use the project's adopted testing and publication rules. Preparing a test,
running it, integration and a push are different actions. The pilot's
one-time server-suite waiver is not permission to skip a later required
suite. Auto mode runs one job at a time; there is no parallel-job
coordination. Follow the current project restrictions and sharing rules.

### Auto mode

Type `/gc auto` yourself (pick it from the command menu) and then paste
the job text after it. A pasted message that merely begins with the
command does not start the skill. Claude turns the request into a job
list, shows one checked start notice and then works without questions.
The start notice gives the absolute path of the queue helper in use, the
stop command and the run boundary.

**What the invocation grants.** It is your advance grant for the listed
jobs only, filed as `GRANT.md` in the queue's record before any job
starts. It covers isolated worktrees and branches from the recorded
baseline, edits within each job's scope, local checks in the test-only
installation, reviewer and checker subagent readings, frozen evidence,
tickets and parking. It waives the per-job plan approval and the
outside reading before parking; both remain available when you review.
It does not cover integration, publication, installation updates, new
server access, accepting an unchosen risk or choices outside a job's
text. It ends when every job is parked or blocked, at a stop, or at the
run boundary (by default 08:00 local time the next morning).

**Each job.** One worker at a time prepares the job in its own worktree,
records the test criterion before any test and runs every build and test
through the queue helper. Each run has a time limit; a run that exceeds
it is ended with its identified processes and recorded as timed out,
which is a failure. Unless the job's grant allows outside contact, the
helper blocks outbound connections and name lookups for its runs. Tests
that would contact outside services are screened out first. The frozen
packet must contain `SCREENING.txt`, one line per test considered:
`test: NAME | outside: WHAT | included` or `excluded`, where WHAT is
exactly `none` or a description of the outside access found. The helper
refuses **Ready for approval** when that file is missing or malformed,
when no test is included, or when an included test has outside access
without the job's network grant. It checks the form, not whether the
screening is complete. A reviewer subagent reads the candidate before
the freeze and a checker subagent reads the frozen packet. The job is
then parked as **Ready for approval** or **Blocked** with one reason and,
where needed, the decision you must make. Master and the working
installations are not changed. The ticket states that the Guide chose
the criterion.

**Morning report.** `/gc auto status` shows one row per job and, under
APPROVAL SUMMARIES, a seven-line block for each ready job: its title,
then the bug, the fix, the test, the criterion, the limits, and a line
saying that approval accepts this exact packet and merges nothing. When
no job is ready, that section says so. The checker checks the form of
the report, not the truth of the summaries. Anything done outside a
job's agreed discretion is asked as its own question before that job can
be approved. You can approve several jobs in one message. Each approval
is recorded with your words and bound to the packet's frozen manifest;
the helper refuses it if the packet has changed.

**Stop.** Send `/gc stop`, or run the helper's stop command shown in the
start notice from any Terminal; it works in any shell, csh included. It
reports **Stopped** only when the registered worker has ended and no
build or test started through the helper, or any descendant it can
identify, is still alive; otherwise it reports **Stopping** and what
remains. Known gaps: a system program that a test starts in its own
session, from a parent that exits within about a second, is not seen;
and a request that a test hands to another system service is not
blocked by the network block.

**Usage limits and `/gc auto resume`.** Recovery after a usage limit is
manual in this release. A limit stops the queue where it is; after the
reset, send `/gc auto resume`. Whether the Claude app retries an
interrupted turn by itself, and what happens to a running worker then,
is untested. The start notice says so. `/gc auto resume` reconciles the
queue, rechecks the test installation and continues the remaining jobs
only while the grant has not ended. It never registers a second worker.
In a new conversation it keeps the queue stopping until you confirm that
the earlier session has ended.

**Not in this release.** There is no `auto integrate`. Integration and
publication of approved jobs are ordinary guided work on your later
instruction, with your PUBLISH decision for the combined packet.

## Ordinary conversations and session titles

Start a fresh conversation without invoking the skill for ordinary work.
Persistent project instructions and auto-memory can still load there;
`/clear` is not a cleanup of those files. Client bookkeeping may also be
written outside the target. See [Claude Code memory](https://code.claude.com/docs/en/memory).

A session's automatic title is not evidence of which skill loaded. Once
the source check passes, the skill names the session
`GuidedCoding: <task>` through the client's rename tool where one exists
(Claude Code Desktop); a title you set yourself is left alone, and
declining the rename is final for that session. On the CLI, which has no
such tool, the skill prints the `/rename GuidedCoding: <task>` instruction
once; in Desktop you can also click the session title. A failed rename
never stops the work. These are client features, not GC commands; see
[sessions](https://code.claude.com/docs/en/sessions#name-your-sessions) and
[Desktop](https://code.claude.com/docs/en/desktop). The
[verification record](GUIDED_CODING_VERIFICATION.md) lists the recorded
title observations and their limits.

## Update, publish or uninstall

Because the personal link points at a checkout, changing the source there
changes what future invocations read. Reconcile active work, preserve a
recoverable baseline and use the project's normal review/integration rules.
A worker that detects changed central source must stop dependent work;
registration does not pin an immutable copy for it.

For publication to `cctbx_project`, tags are reserved for repository releases.
Do not recreate the removed GC pilot tag or create another development tag.
Keep pilot status in these docs and refer to the publication commit and
manifest. Other destinations have their own repository conventions.
Nothing in registration or a successful checker authorizes a push.

To uninstall, use the control request above or ask ordinary Claude Code
to inspect the personal link from this guide's absolute path. Remove only
a confirmed symlink to this exact verified package. Do not remove a real
directory, an unknown link, the central source, project methods or records.
Use a fresh conversation afterward. There is no per-project installer,
automatic profile deletion, automatic lock guard or `/wrap` command.
