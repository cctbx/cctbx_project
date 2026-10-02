# GuidedCoding user guide — opt in with `/gc` (r10 rev16 candidate)

## The short version

Keep one reviewed copy in `cctbx_project/libtbx/guided_coding`. Register
it **once per local Claude Code configuration**, then use it in any local
repository. You do not install the procedure in every repository.

After the candidate has been reviewed and integrated into your actual
`cctbx_project`, open Claude Code there and send this message, replacing
the path with that checkout's real absolute path:

```text
Please set up GuidedCoding from /absolute/path/to/cctbx_project/libtbx/guided_coding for this machine. Read docs/GUIDED_CODING_USER_GUIDE.md there. Verify the source and release, run the Claude Code minimum-version check, inspect any existing personal skill, and register the one central skill link only if those checks pass and it is safe. Show what you changed, then offer /gc setup for the project I choose, recovering any saved settings first. Do not connect to servers, change permission settings, or start a coding task.
```

This is a normal Claude Code request, **not** a built-in installer. The
absolute path is necessary the first time: Claude Code in another project
cannot guess where your `cctbx_project` checkout lives. Review the path
and any proposed change to an existing skill link. If the `skills/`
directory was created after Claude Code started, use `/reload-skills` or
open a new conversation. Then confirm the command itself: type `/` (in the
Mac trial the Desktop menu listed `gc`) or run `/gc help`. Do not rely on
`/skills`: in the Mac trial (Claude Code 2.1.284) Desktop's `/skills` did
not list this personal skill even though `/gc` worked.

In a project you want to use with GuidedCoding, start a new conversation
in that project's folder and enter `/gc <task>`. If its contract adoption
or method is missing, or the requested work needs incomplete setup, the
skill walks through relevant gaps using saved values as proposed defaults
before dependent work. You can use `/gc setup` first if you want to
prepare the project without starting a task. An ordinary new
conversation without `/gc` remains ordinary.

## Review a candidate archive

The archive contains paths beginning with `libtbx/guided_coding/`. It has
**no enclosing candidate-named directory**. To inspect it safely, create
a new staging folder of your choice and extract there. First compare the
archive SHA-256 with the independently supplied review message using trusted
tools; do not extract a rejected archive. Replace the paths and run the whole
Bash block below. An existing staging folder is refused:

```bash
mkdir /path/to/empty-gc-review &&
tar -xzf /path/to/GuidedCoding_2_0_opt_in_r10_candidate_rev16_f1-correction-20261001T172033Z.tgz -C /path/to/empty-gc-review &&
cd /path/to/empty-gc-review/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev16 candidate' payload/RELEASE
```

A failed creation, extraction, directory change or check returns nonzero
and skips all later commands in the block. Stop dependent work if it fails;
do not append preparation or registration as unguarded lines.
Every listed file should say `OK`, and the complete-source check should
print `VERIFIED complete source`. The release line must identify **r10
rev16**. The supplied rev14 base and failed rev15 are historical copies. The manifest checks the
bytes listed **by that manifest**; the second check refuses an extra
file or link in the effective source. On an installation whose Python
files were precompiled (for example by an installer's
`libtbx.py_compile_all`), it also prints a NOTE that it ignored the
`__pycache__` bytecode beside listed modules; the tools never load that
bytecode, and any other unlisted file, including other bytecode, is still
refused. Accepted bytecode is checked by name and header type, not by
content: running the package's tests yourself with Python imports through
that cache, so run them only on a trusted copy. Neither can establish who supplied the
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
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev16 candidate' payload/RELEASE &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py check-claude-version &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py register-skill .
```

The `&&` sequence does not create the link if any check fails. Confirm its
four checks: listed files say `OK`, the complete inventory passes, the
release reads `r10 rev16`, and the CLI
check says `VERIFIED Claude Code CLI` for version 2.1.281 or newer. The historically reported
2.1.284 meets that minimum; check the current executable rather than assuming it. A missing, unreadable or older CLI stops
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
change before migrating an older installation. After successful registration,
offer `/gc setup` for the project the user chooses; the registration helper
prints this next step but performs no project adoption or project setup.
If no target has been named, ask for the project directory rather than
assuming the central procedure checkout is the target.

### Retire old automatic startup instructions

Before relying on the opt-in choice, inspect each repository's existing
`CLAUDE.md`, `CLAUDE.local.md`, `.claude/rules/`, and any user-level
instruction for a line that says to read `.claude/WORKER.md` or to treat
all coding tasks as GuidedCoding. Replace only that startup instruction
with a neutral note such as “Use GuidedCoding only when I type `/gc`.”
Preserve unrelated build paths, machine facts, existing approvals, and
records. Follow SETUP.md before retiring old profiles: map each setting,
restriction and pending item into the retained or current method and verify
its saved value and source. Stale verification never justifies erasing a
known value; keep it marked as needing verification. Verify the central source and record any reviewed profile edit
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

### What `/gc setup` does

Setup follows the central [setup workflow](../payload/SETUP.md). The Guide
first finds the current method and its referenced settings, handoffs and
migration backups. It uses recovered values as proposed defaults and tells
you where they came from. It asks only for missing choices or facts it
cannot inspect, in small steps appropriate to this project.

The package's [general defaults](../payload/SETUP_DEFAULTS.md) contain
project-neutral proposals, not machine-specific settings or permission.
You may attach a separate profile supplied by your project, for example a
PHENIX setup profile, and ask the Guide to use it as candidate defaults.
The Guide records its source, adapts account/path choices to your own setup,
and establishes your own permissions and resource limits. Such a profile
is not loaded for unrelated projects and does not supersede your current
project method. Relevant values and provenance are saved in your method so
later sessions do not depend on retaining the attachment.

The same workflow is reached when a guided task needs incomplete setup,
even in an already adopted project. A routine startup reads the saved
method and checks only what the requested operation needs; it does not
repeat the interview. A different machine, account, checkout or execution
shell may require checking the affected entries again. `/gc status` shows
the method location and unresolved setup by activity without changing it.

| Area | What the Guide establishes |
| --- | --- |
| Project and instructions | Target checkout, related repositories and scope, actual loaded authority, canonical method path, personal versus shared adoption |
| Working directories | Source, permitted worktree/scratch locations, durable records/results, and remote working directories only when needed |
| Local environment and commands | Shell, executable/PATH or activation setup, exact build/test commands, working directory, import target and known runtime |
| Optional CI or remote work | Only relevant services/hosts, access method, installation/workspace, commands, log retrieval, approved limits and sharing/lock rules |
| Actions and restrictions | Existing decisions about verification, tests, installation changes, integration and publication; what still needs a decision |

For a non-PHENIX project there are no default PHENIX commands or servers.
If it needs only local tests, remote work is not applicable or is deferred
by choice. A method using a name such as `t96` keeps the name and its actual
definition; the Guide inspects that definition before asking you to recall
it. Every host has its own settings: a recorded hardware count is not an
approved concurrency cap, and one server's paths do not configure another.

Setup distinguishes recorded values, verified facts, unknowns, conflicts,
deferred work, and things that are not applicable. A saved path stays in the
method when it needs verification. It does not become "unknown". A tool
absent from bash's PATH may already be installed and configured in csh;
the Guide reads relevant definitions and verifies the intended execution
environment. It does not bypass hooks or change global startup files to
work around a lookup failure.

Project-specific settings stay in the existing project method, or in
`.claude/PROJECT_METHOD.md` when that is the chosen new location. Existing
shared defaults may be referenced with explicit project overrides. Setup
does not create a new global profile, copy procedure prompts into the
project, or change which instruction files load accidentally. The optional
[method template](../payload/templates/PROJECT_METHOD.md) is guidance, not
a form every project must fill out.

The Guide shows the complete proposed files and any adoption/pointer change
before saving unless that exact bounded edit is already authorized. It then
writes the approved texts and the prior bytes into a setup record and reads
them back before editing anything, saves only the authorized changes, checks
readback and the method pointer, and gives you their exact locations. It
records sources, decisions and unresolved items so the next Guide can
continue. During migration, every old setting/restriction is retained,
mapped into the current method, or explicitly retired by your decision;
old profiles and records are not deleted by setup.

Default setup is local read-only discovery and an authorized configuration
save. It runs no build or suite, opens no server connection, and starts no
coding task. If you request further verification or setup actions, the
Guide reuses applicable authorization or asks for the particular action
still needing it. No credentials are stored in the method. Preparing tests,
running them and pushing are distinct: a no-push request does not prevent
configuration discovery, and a historical publication-only test rule still
needs your decision before a test-only run outside that rule.

Setup ends by saying what is configured, what remains deferred or needs
verification, and the smallest next step. Missing optional remote setup
does not prevent independent authorized local work. A project may continue
ordinary Claude Code work without adopting GuidedCoding.

## A guided task

Open Claude Code in the repository you want to change. Type, for example:

```text
/gc Fix the broken copy_extra regression test.
```

Or type `/guided_coding Fix the broken copy_extra regression test.` The
directory name supplies the longer alias; the skill's `name: gc` supplies
the short command. Both run the same skill; a session transcript records
either as `/guided_coding`. If `/gc help` does not run the skill, check
the symlink target and your Claude Code version (`/skills` may not list
it). The personal skill is local
to your machine and is unavailable to cloud sessions that cannot read it.

The first task invocation checks the central source and the target's
contract adoption and project method before reading its Guide and Worker
as controlling. It then checks setup for the requested operation, recovering
saved values and walking through any relevant gaps using SETUP.md. It then records the release and source hash in the current
project. The first guided change in a repository receives a short welcome.
`/gc` alone introduces the method and waits for your task; it does not
start a code change. You can ask “How does GuidedCoding work?” at any time
in that guided session.

Other control requests:

| Command | Effect |
| --- | --- |
| `/gc help` | Explain the commands and opt-in behavior. |
| `/gc status` | Read-only check of the central release, personal link, project adoption, method location and unresolved setup for relevant activities. |
| `/gc setup` | Recover defaults, walk through relevant environment/directory setup, and save authorized method/adoption changes. |
| `/gc uninstall` | Remove only the personal `guided_coding` symlink if it points at this exact central copy. The central source and project records stay. |

If you are uninstalling without a working `/gc`, ask ordinary Claude Code
to read this guide from its absolute path and inspect the personal link.
Remove only a confirmed symlink to this package. A real directory, an
unknown link and the central source itself must not be removed by this
instruction. Start a new conversation after uninstalling, since an
already invoked skill remains in the current conversation.

### Where `/gc` appears

`/gc` is a **custom skill in the Slash commands menu**. On Desktop, type
`/` in the prompt or choose **+ → Slash commands**; in the Mac trial the
menu listed `gc` (not `guided_coding`), while `/skills` listed no personal
skill. Typing `/gc` or `/guided_coding` worked either way. It is not an additional mode in the selector containing
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
