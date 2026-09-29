---
name: gc
description: Run GuidedCoding for the user's explicitly chosen coding task.
disable-model-invocation: true
---

# GuidedCoding, by explicit request

The user's task is: $ARGUMENTS

This is the one manual entry point to the central GuidedCoding source.
The skill directory is `${CLAUDE_SKILL_DIR}`. Resolve that directory's
canonical absolute path, including any personal-skill symlink; that is
the procedure root. Its `payload/` directory is the procedure source.
The session's project directory is `${CLAUDE_PROJECT_DIR}`. **The target
repository is the repository containing that session directory when this
command was invoked**, unless the user explicitly names another target. Never treat the
procedure repository as the target merely because it holds this skill.

Before guided work, from the canonical procedure root run
`shasum -a 256 -c SOURCE_MANIFEST.sha256`, then run
`python3 -I -B payload/tools/screen_check.py verify-source .`.
The first command checks the listed checker file before it runs; the
second checks the **complete source inventory** and rejects unlisted files,
links and changed bytes. The isolated Python invocation prevents an
unlisted module in `payload/tools/` from loading during verification.
If either fails, stop. Then read `payload/RELEASE` and
`payload/DEVELOPER_GUIDE_CONTRACT.md`.
This candidate requires the release label to say `r10 rev11`; if it does not,
stop and report the different version before starting a guided change.

## Setup and control requests

Treat exactly `help`, `status`, `setup`, `setup <target directory>`, and
`uninstall` as control requests rather than coding tasks. Read
`docs/GUIDED_CODING_USER_GUIDE.md` for their exact scope. `help` gives a
short command list; `status` reads the central identity, personal link and
target adoption without changing anything. Neither needs project adoption.
Before `setup` or an ordinary guided task, run
`python3 -I -B payload/tools/screen_check.py check-claude-version` from the
verified central root. If it fails, stop before writing project adoption or
starting guided work; do not update Claude Code on the user's behalf. The
check requires CLI 2.1.281 or newer on this session's PATH. It does not
establish the Desktop app's own instruction-loading behavior, which needs a
separate live check. `help` and `status` are read-only; `uninstall` may
remove only the verified personal link. All three remain available when
the CLI check fails.

For `setup`, inspect the target's applicable project instructions and
existing method first, including `AGENTS.md` and Claude Code's instruction
loading setting. Do not silently create `CLAUDE.local.md` in a project
that relies on `AGENTS.md`; under the default setting this can stop
Claude Code loading `AGENTS.md`. Propose a short project-specific method
covering the actual checkout, local build and test commands, coordination and
publication rules. Include optional server names, remote installation and
lock paths, and commands such as `t96` only when verified and relevant to
that target. Do not infer their syntax or assume a PHENIX environment in an
unrelated project. Ask for only the facts inspection cannot establish.
Show the exact proposed project method and contract declaration for the
Developer's project-adoption decision. On authorization, update only the
chosen current project authority and method, preserving unrelated
instructions; verify what was saved. `setup` does not run a build, contact
a server, grant tool permissions or begin the coding task. If adoption is
declined, leave the project unchanged and offer ordinary assistance.

For `uninstall`, inspect the personal `skills/guided_coding` entry under
`CLAUDE_CONFIG_DIR` if set, otherwise `~/.claude`. Remove **only** that
entry if it is a symlink whose canonical target is this verified procedure
root; refuse a real directory, a different link target, or an uncertain
path. Never remove the source, project authorities, methods or records.
Explain that already loaded instructions persist in this conversation,
and use a fresh conversation for ordinary work. No project adoption is
needed to uninstall the personal link.

After handling a control request, stop; do not fall through into the
guided-task activation or start work on an unrelated task.
For an ordinary `/gc <task>`, continue with the activation check below.
The contract's project activation rule is a separate gate: inspect the
**target repository's current authority** (including its applicable
`CLAUDE.md`, `CLAUDE.local.md`, an actually loaded `AGENTS.md` if used, and
any project-designated authority). It
must explicitly adopt the contract version **2026-09-17**, SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`,
and define or identify that project's method separately. Compare the
declared identity with the verified central contract bytes. A personal
skill link, `/gc` invocation, old records, a bundled contract, or this
task's approval alone does not establish project adoption. If the
declaration, project method or current authority is absent or uncertain,
pause **before** claiming that Guide and Worker govern the target or
making a guided change. Show the Developer a small proposed adoption
declaration naming the exact contract identity and that project's
specific method; obtain the Developer's project adoption decision.
Only on explicit authorization, record it in the current project
authority chosen by the Developer, preserving existing instructions,
then verify the saved declaration and method. If adoption is declined,
offer ordinary assistance without claiming GuidedCoding controls it.

Once adoption is verified, read `payload/ROLES.md`, `payload/GUIDE.md`
and `payload/WORKER.md`. The Guide establishes scope and roles; the Worker
follows the bounded change procedure.
The procedure is read-only for the task. Recheck both the listed hashes
and the complete inventory before integration, publication, or running a
central tool. Use `python3 -I -B` for later central tool invocations so
local modules cannot shadow the standard library and the tools do not
create unlisted bytecode in the source directory. If verification fails
or its source changes while the task runs, stop dependent work and report
the mismatch.
Record the procedure release and the SHA-256 of `SOURCE_MANIFEST.sha256`
in the target's local change record so a later reviewer knows which exact
instructions and tools were used.

Set `GC_PAYLOAD_ROOT` to the **canonical absolute `payload/` path inside
each shell command** that calls a procedure tool; shell state does not carry
between tool calls. Read screens, templates, briefs and Python tools from
that central directory. Keep each project's state and records in that
project's own `.claude/records/`. Never copy procedure instructions into
the target repository. Do not overwrite its existing records, settings,
profile, permissions, or other project instructions. Check its existing
approval and any revocation before consequential actions.

If there is no task after `/gc`, give the first-use introduction and ask for
the task without starting a change; do not demand project adoption for an
introduction or explanation. If asked how GuidedCoding works, answer
without opening a change. A normal prompt in a fresh conversation does not
invoke this manual-only skill. The instructions remain in *this* conversation
once invoked; a later ordinary task needs a fresh conversation or `/clear`.
Neither this command nor its symlink authorizes installation, integration,
server work, outside delivery or publication.
