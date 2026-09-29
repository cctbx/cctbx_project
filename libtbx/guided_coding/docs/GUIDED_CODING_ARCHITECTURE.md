# GuidedCoding r10 rev11 architecture

The **one central source** is `cctbx_project/libtbx/guided_coding/`. A
personal skill symlink at `~/.claude/skills/guided_coding` points at this
directory. A review archive contains these files directly under
`libtbx/guided_coding/`, without an enclosing candidate directory. Its
`SKILL.md` declares `name: gc` and
`disable-model-invocation: true`, so a person invokes `/gc` or the
directory alias `/guided_coding`; Claude does not choose to start it.
Registration is one personal symlink, optionally created by an ordinary
Claude Code setup request given the central source's absolute path. The
`register-skill` helper verifies the source and CLI, then creates that exact
link exclusively. An existing link (even to this source), directory, dangling
link or symlinked skills directory is refused without a nested link.
An absolute `CLAUDE_CONFIG_DIR` containing `..` is rejected before any
directory creation. This keeps preflight checks and creation on the same
effective path, even if an earlier component is absent.
`/gc help`, `/gc status`, `/gc setup`, and `/gc uninstall` are handled
before the guided-task adoption gate. No libtbx dispatcher or project
installer is required. Setup proposes project adoption and a method;
uninstall removes only a symlink resolving to this verified central
source. Neither command deletes project records.

| Path | Purpose | Read or written when |
| --- | --- | --- |
| `docs/` | README, user guide, architecture and verification | Read from the central source when needed |
| `SKILL.md` | Resolve the target repository and load the central procedure | Only when invoked |
| `SOURCE_MANIFEST.sha256`, complete-source check and `payload/RELEASE` | Identify listed source bytes, reject unlisted source files, and display release label | Check all at start; recheck before central tools and consequential work |
| `screen_check.py check-claude-version` | Require a readable Claude Code CLI at version 2.1.281 or newer | Run after source verification, before personal link registration, project setup, or guided work; help, status, and uninstall remain available |
| `screen_check.py register-skill .` | Recheck source and CLI, then create only the unoccupied personal skill link | One-time registration; refuses existing destinations without changing them |
| `payload/{GUIDE,WORKER,ROLES,DEVELOPER_GUIDE_CONTRACT}.md` | Scope, steps and authority | Read from central source |
| `payload/screens/`, `payload/templates/` | Decision formats and complete approval report | Read from central source |
| `payload/tools/` | Screen, evidence and review-bundle checks | Run from central source |
| Target's `.claude/records/` | Per-repository plans, decisions and proof | Written for that repository only |
| Target's current authority and project method | Explicit adoption of exact contract version and SHA-256; project-specific working rules | Checked before Guide and Worker govern a project; preserve actually loaded `AGENTS.md` instructions |
| Target's `CLAUDE.local.md` and settings | Local machine facts and permissions | Read if present; no automatic GC instruction |

The project method may contain build and test commands, a `t96` command,
and specific servers such as `anaconda.lbl.gov` or `cci-gpu-00.lbl.gov`
when they are relevant and verified. It is a small project-specific
working document, not a global installer configuration or a grant to
run those commands. A non-PHENIX target substitutes its own environment,
build, CI and coordination method; it does not require libtbx. One
personal link exposes `/gc` in local projects, while each project's
contract adoption remains an explicit choice.
In a project that relies on `AGENTS.md`, adding a `CLAUDE.local.md` may
change which instructions Claude Code loads by default. Setup checks the
active instruction sources before proposing an adoption location.
The minimum-version check invokes the `claude` executable found on the
current shell's PATH and refuses missing, malformed, failing or older
versions before setup changes. It does not update Claude Code or verify
the Desktop app's instruction loading. A fresh Desktop and `AGENTS.md`
check is still required for use there.
The version check parses the executable's reported stdout; it cannot prove
the executable's identity. A misleading program on PATH could report an
acceptable version. Review the path shown by the check as part of live setup.

The invocation selects the repository of Claude Code's primary working
directory unless the user names a different repository. Reading a skill
from `cctbx_project` does not change the target to `cctbx_project`. Each
central tool call uses the absolute `payload/tools` path; its output and
the record stay with the target. From the central root run `shasum -a 256
-c SOURCE_MANIFEST.sha256`, then `python3 -I -B
payload/tools/screen_check.py verify-source .`. The first checks the listed
checker bytes; the second checks the full inventory without importing code
from the tools directory. Read `payload/RELEASE` for the label. These
checks do not authenticate who supplied the package; compare the
archive to a trusted review message before accepting it. The target's
change record names the source manifest's SHA-256. If central files
change during a task, the worker stops dependent work rather than
silently switching procedures.

**Activation is a separate check.** The central contract version
`2026-09-17` has SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`.
Before guided work the skill reads the target's applicable current
authority and verifies that it explicitly adopts that exact identity
and defines or identifies its own project method. An `/gc` invocation,
personal symlink, bundled contract, or old record does not establish
adoption. If missing or ambiguous, the skill pauses and proposes a
small project-authority declaration for the Developer to approve; it
does not write one or treat the Worker as controlling automatically.
This keeps procedure code central while requiring a project decision.

The existing r09 contract and five roles remain: the Developer owns
meaning, value, risk, waiver, acceptance, and consequential action; a
Guide scopes and directs; a Worker handles bounded work; an Outside
Reviewer evaluates a frozen packet; a Helper is the Developer's own
optional explainer. No registration step transfers authority. The
Worker checks that tests import the intended checkout, shows a concise
plan, proves the candidate with controls, obtains the required review,
then presents the exact changes for local integration. Publication has
its own suite, review, and decision.

`payload/tools/screen_check.py` requires a full approval report with
**`## Exact tested changes`** and **`## New tests`** headings and the
complete `CHANGE.diff`. Its evidence `verify` command says whether the
frozen packet matches `MANIFEST.sha256`; it does **not** print that
manifest's hash. To record the evidence packet's identity, run
`shasum -a 256 MANIFEST.sha256` inside the packet (or the equivalent
SHA-256 command). The source and evidence manifests are distinct.
The review-bundle tool requires companions outside the packet and permits
output only in the **frozen packet's own parent directory**. It refuses
a nested companions directory before archive creation. It compares the
opened output directory with that
parent by filesystem identity, keeps the handle open through creation,
and confirms the packet still occupies its original directory before
writing. A separately movable directory is refused before any write,
even if it was outside the packet at first. A symlink alias to the
packet's parent can be used; swapping that alias does not redirect the
opened handle. Moving the packet parent under its child is prohibited by
the filesystem while the packet stays there. Concurrent hostile moves
of the packet itself are outside this tool's guarantee; use an attended
workspace without concurrent changes. Live Mac case-variant and
Desktop tests remain necessary before adoption.

Claude Code automatically reads project and personal `CLAUDE.md` files.
An older profile that instructs it to read `.claude/WORKER.md` for every
task overrides the opt-in intent until that line is migrated. Likewise,
skill content persists after invocation inside the conversation; use a
fresh conversation for ordinary work. The skill has no preapproved tools,
and ordinary Claude Code permission rules continue to apply.

This candidate has no r09 `install`/`verify` command, automatic profile
creation step, `/wrap` command, or automatic shared-tree lock guard.
Finish or discard an old open change in its Worker session before a
reviewed profile migration. Use the normal `cctbx_project` controls for
central source changes and coordinate with any working installation.

This candidate retires the per-repository installer and its test, plus
three obsolete profile/migration examples and its acceptance-hash list.
Two checker tests move alongside the central tools. Historical procedure
copies and evidence in previously configured repositories are **not**
deleted by this package; they are migrated separately after inspection.

`/gc` is a Claude Code skill in the slash-command picker, not a new
permission mode next to Auto. Its manual-only content is loaded for each
guided conversation and remains in that conversation; it does not have
to be retyped for each message. Auto and GuidedCoding can be selected
independently. Always-loaded project instructions could turn a procedure
on for every task, but that would change this version's explicit opt-in
contract.
