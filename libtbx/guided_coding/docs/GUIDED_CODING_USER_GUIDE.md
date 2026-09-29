# GuidedCoding user guide — opt in with `/gc` (r10 rev11 candidate)

## The short version

Keep one reviewed copy in `cctbx_project/libtbx/guided_coding`. Register
it **once per local Claude Code configuration**, then use it in any local
repository. You do not install the procedure in every repository.

After the candidate has been reviewed and integrated into your actual
`cctbx_project`, open Claude Code there and send this message, replacing
the path with that checkout's real absolute path:

```text
Please set up GuidedCoding from /absolute/path/to/cctbx_project/libtbx/guided_coding for this Mac. Read docs/GUIDED_CODING_USER_GUIDE.md there. Verify the source and release, run the Claude Code minimum-version check, inspect any existing personal skill, and register the one central skill link only if those checks pass and it is safe. Show what you changed. Do not connect to servers, change permission settings, or start a coding task.
```

This is a normal Claude Code request, **not** a built-in installer. The
absolute path is necessary the first time: Claude Code in another project
cannot guess where your `cctbx_project` checkout lives. Review the path
and any proposed change to an existing skill link. If the `skills/`
directory was created after Claude Code started, use `/reload-skills` or
open a new conversation; confirm the skill appears in `/skills`.

In a project you want to use with GuidedCoding, start a new conversation
in that project's folder and enter `/gc <task>`. If its contract adoption
or method is missing, the skill proposes the missing setup for your
decision before guided work. You can use `/gc setup` first if you want to
prepare the project without starting a task. An ordinary new
conversation without `/gc` remains ordinary.

## Review a candidate archive

The archive contains paths beginning with `libtbx/guided_coding/`. It has
**no enclosing candidate-named directory**. To inspect it safely, create
an empty staging folder of your choice and extract there:

```bash
mkdir -p /path/to/empty-gc-review
tar -xzf /path/to/GuidedCoding_2_0_opt_in_r10_candidate_rev11.tgz -C /path/to/empty-gc-review
cd /path/to/empty-gc-review/libtbx/guided_coding
shasum -a 256 -c SOURCE_MANIFEST.sha256
python3 -I -B payload/tools/screen_check.py verify-source .
cat payload/RELEASE
```

Every listed file should say `OK`, and the complete-source check should
print `VERIFIED complete source`. The release line must identify **r10
rev11**, rather than r08, r09, or another copy. The manifest checks the
bytes listed **by that manifest**; the second check refuses an extra
file or link in the effective source. Neither can establish who supplied the
manifest or authenticate a different checkout. Compare the archive's
SHA-256 to the independently supplied review message before trusting it.
No `libtbx.install_guided_coding` dispatcher, per-repository installer,
or automatic profile creation ships in this candidate. The one-time
personal skill registration is separate from each project's method and
contract adoption.

## One-time registration on your Mac

After release review and an authorized integration place the complete,
reviewed directory at `cctbx_project/libtbx/guided_coding`. Then register
the **one central copy** on your machine (replace the example path with
your real checkout). Inspect any existing skill path first. Registration
refuses an occupied path, including a link already pointing to this source;
it never overwrites or creates a nested link.

```bash
cd /path/to/cctbx_project/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev11 candidate' payload/RELEASE &&
python3 -I -B payload/tools/screen_check.py check-claude-version &&
python3 -I -B payload/tools/screen_check.py register-skill .
```

The `&&` sequence does not create the link if any check fails. Confirm its
four checks: listed files say `OK`, the complete inventory passes, the
release reads `r10 rev11`, and the CLI
check says `VERIFIED Claude Code CLI` for version 2.1.281 or newer. Your
2.1.284 meets that minimum. A missing, unreadable or older CLI stops
registration; run `claude update` yourself and restart before retrying.
The registration helper independently rechecks the complete source and
CLI version before creating the link at the exact destination. A second
run stops with `already exists`, without changing the source or link.
The check uses the `claude` executable on this shell's PATH; it does not
prove which instruction files the Desktop app loads. Verify Desktop and
`AGENTS.md` in fresh sessions separately. The paths in the manifest are
relative to this directory.
The symlink is **personal configuration**, not an install in every repo.
Neither making the link nor running `/gc` checks a shared tree lock.
For a central source update, coordinate the change in `cctbx_project`
under that project's normal repository and installation controls.

If you set an **absolute** `CLAUDE_CONFIG_DIR` with no `..` path components,
the registration helper puts
the `skills/guided_coding` symlink there instead of under `~/.claude`. It
rejects `..` before creating any directory, and refuses a symlinked
configuration or skills directory. The central checkout must exist
and be readable from each local Claude Code session. The link registers the
command once for all local projects; it does not copy the procedure into
those projects.

If the link already exists, inspect its type and destination. Leave a
real directory or a link to another package alone; do not overwrite or
rename it as an automatic upgrade. Show the user exactly what would
change before migrating an older installation.

### Retire old automatic startup instructions

Before relying on the opt-in choice, inspect each repository's existing
`CLAUDE.md`, `CLAUDE.local.md`, `.claude/rules/`, and any user-level
instruction for a line that says to read `.claude/WORKER.md` or to treat
all coding tasks as GuidedCoding. Replace only that startup instruction
with a neutral note such as “Use GuidedCoding only when I type `/gc`.”
Preserve unrelated build paths, machine facts, existing approvals, and
records. Verify the central source and record any reviewed profile edit
separately. Historical `.claude/` records can remain. Remove an older
procedure copy only after checking its records and review references.

### Adopt the contract in each target project

The **procedure stays in one central directory**. A project still needs
its own explicit adoption of the Developer–Guide Contract before `/gc`
can govern work there. Check the project's current authority first (its
applicable `CLAUDE.md`, `CLAUDE.local.md`, an actually loaded `AGENTS.md`,
or another authority the project has designated). It must name contract
version `2026-09-17` and SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`
and separately define or point to the project's own working method.
An older adoption of **this exact contract identity** can satisfy the
gate; an older release label alone does not. The skill link and typing
`/gc` do not adopt a contract for the project.

For a project without adoption, `/gc <task>` stops before guided work and
offers a short declaration for your review. If you authorize adoption,
place it in that project's current authority, preserving existing facts
and grants. For example, if `CLAUDE.local.md` is the chosen authority and
`PROJECT_METHOD.md` is an already reviewed project method, add a small
section like this (use the project's *actual* method path):

```markdown
GuidedCoding is used in this repository only after explicit /gc invocation.
Adopt Developer–Guide Contract version 2026-09-17,
SHA-256 ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b.
Project method: follow PROJECT_METHOD.md for this repository's build,
test, coordination and publication rules.
```

If no project method exists, agree on its content and record it separately
before guided work. A local adoption affects only that checkout; a shared
team adoption belongs in its designated, reviewed project authority.
This is a short project decision, not a copy or installation of the kit.
If the project relies on `AGENTS.md`, do not blindly create
`CLAUDE.local.md`: Claude Code's default instruction loading can then
skip `AGENTS.md`. Inspect the active instruction setting and propose a
location that preserves the existing instructions before saving anything.

### What `/gc setup` asks about

The command reads existing project instructions first. It separates
facts about this checkout from the shared GuidedCoding contract, proposes
the smallest missing method, and shows the exact project-authority text
for your decision. It asks only for information it cannot verify. The
following are **options**, not required fields in every project:

| Project fact | When to record it | Example or question |
| --- | --- | --- |
| Target directory and authority | Always | Which checkout is the target, and which instruction file is actually loaded? Is adoption personal or shared? |
| Local environment, build and tests | If applicable | Exact setup/build/test commands, working directory and expected runtime; do not guess command syntax. |
| `t96` | Only if this project's method uses it | What exact command runs, where, and what does it validate? |
| `anaconda.lbl.gov` | Only if this project uses that server | Login alias, installation and scratch paths, lock rule, and whether a task may use it. |
| `cci-gpu-00.lbl.gov` | Only if this project uses that server | Same questions independently; do not assume anaconda's paths or approvals. |
| Other machines or non-PHENIX tools | Only if relevant | The actual CI, container, build system or test command for this project. |
| Coordination and publication | Always where applicable | Existing lock/branch/review and publication process; preserve any stronger project rule. |

Host names and commands belong to the **project's method**, not the
central skill. Recording a server does not grant permission to connect,
change its installation or publish. Do not put passwords, tokens or SSH
keys in the method. Setup does not run `t96`, build, tests or remote
commands. If a command or path is unknown, mark it unknown and ask before
the first task that needs it. Do not invent a default for non-PHENIX
projects.

For a non-PHENIX tree, start Claude Code in that tree and use `/gc setup`
the same way. Its method may name `pytest`, `make`, a web build, CI or no
remote server at all. The procedure source still lives at the central
`cctbx_project` path and the records for a guided change live in the
**target** tree. The target need not import `libtbx` or have a `.claude`
procedure copy. A repository without an agreed project method can use
ordinary Claude Code until its method is settled.

## A guided task

Open Claude Code in the repository you want to change. Type, for example:

```text
/gc Fix the broken copy_extra regression test.
```

Or type `/guided_coding Fix the broken copy_extra regression test.` The
directory name supplies the longer alias; the skill's `name: gc` supplies
the short command. If either alias does not appear, check `/skills`, the
symlink target, and your Claude Code version. The personal skill is local
to your machine and is unavailable to cloud sessions that cannot read it.

The first task invocation checks the central source and the target's
contract adoption and project method before reading its Guide and Worker
as controlling. It then records the release and source hash in the current
project. The first guided change in a repository receives a short welcome.
`/gc` alone introduces the method and waits for your task; it does not
start a code change. You can ask “How does GuidedCoding work?” at any time
in that guided session.

Other control requests:

| Command | Effect |
| --- | --- |
| `/gc help` | Explain the commands and opt-in behavior. |
| `/gc status` | Read-only check of the central release, personal link and this project's adoption. |
| `/gc setup` | Prepare and, after your adoption decision, save this project's method and declaration. |
| `/gc uninstall` | Remove only the personal `guided_coding` symlink if it points at this exact central copy. The central source and project records stay. |

If you are uninstalling without a working `/gc`, ask ordinary Claude Code
to read this guide from its absolute path and inspect the personal link.
Remove only a confirmed symlink to this package. A real directory, an
unknown link and the central source itself must not be removed by this
instruction. Start a new conversation after uninstalling, since an
already invoked skill remains in the current conversation.

### Where `/gc` appears

`/gc` is a **custom skill in the Slash commands menu**. On Desktop, type
`/` in the prompt or choose **+ → Slash commands**; `/skills` lists the
available skills. It is not an additional mode in the selector containing
**Auto**. Auto controls tool permissions, while `/gc` selects an optional
working procedure; both can be used together. The personal link exposes
the skill in each local project, but you invoke it at the start of **each
new guided conversation**, not on every message and not as a separate
installation per project. Its instructions remain in that conversation.

Making GuidedCoding automatic through always-loaded `CLAUDE.md` would
change this explicit-choice design. This candidate does not do that.
The project adoption declaration says GuidedCoding applies only on
explicit `/gc` invocation.

## An ordinary task

Start a **new conversation** in the repository and describe the work
normally, without `/gc`. The opt-in skill is manual-only, so a fresh
ordinary conversation does not load it. For a later ordinary task in the
same chat, use `/clear` first or start another conversation: Claude Code
keeps an invoked skill's instructions in the current conversation.

## Decisions and permissions

The guided task presents a plan for your approval, a tested result and
full approval report for **INTEGRATE / REVISE / DISCARD**, and a separate
**PUBLISH / HOLD** decision if publication is proposed. The Outside
Reviewer advises; you choose. An “allow command” dialog concerns only
the command shown there. `/gc` does not grant permission to use a remote
server, edit another repository, integrate, contact a reviewer, or push.

The approval report includes headings **`## Exact tested changes`**
and **`## New tests`**, the complete tested diff, and the exact code of
each added test (or an explicit statement that none was added).
`screen_check.py verify` reports whether the frozen evidence still
matches its manifest; to record the evidence packet's identity, separately
run `shasum -a 256 MANIFEST.sha256` inside that packet's directory.
When building an Outside Reviewer bundle, write it **beside the frozen
packet** (in the packet's own parent directory). The bundler refuses an
output directory elsewhere, even if it is outside the packet.

If an older project has an unfinished GuidedCoding change, finish or
discard it in that Worker session before migrating its startup profile.
There is no `/wrap` command in this candidate.
