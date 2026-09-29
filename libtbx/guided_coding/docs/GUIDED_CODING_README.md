# GuidedCoding — central opt-in r10 rev13 candidate

**Status (2026-09-29):** r10 rev13 lets the complete-source check accept the
bytecode an installer's precompile step writes beside listed modules (the
tools never load it) and keeps refusing any other unlisted code. r10 rev12
changed only the `/skills` guidance, to
match the rev11 Mac trial; rev11 itself was reviewed (no blocking finding)
and trialled on the intended Mac. Rev12 is committed locally in
`cctbx_project`, not published, and needs its own review before
publication.

SUPERSEDED 2026-09-29 (rev12; see the Rev12 section below): "Status: proposed successor to rejected r10 rev10. This revision refuses
configuration paths containing `..` before any directory creation. It has not
been reviewed as a release,
integrated, or used on a real repository task." ("This revision" is rev11.) The two readings
of the **earlier r09 documentation** are addressed in the revised guide;
they are not reviews of r10 rev9. The r10 rev5 review found that a
movable destination can still enter the packet after a preflight check;
the rev6 destination rule addresses that case. Neither earlier verdict
approves rev11. Rev7 added conversational setup and `/gc` control requests.
The rev7 review found that listed-hash checking accepted unlisted runtime
code and that the bundle helper accepted companions inside the packet.
Rev8 checks the entire source inventory and refuses nested companions.
The rev8 Mac trial passed its automatic and CLI skill checks but found that
Claude Code 2.1.268 did not load `AGENTS.md`; after updating, the Developer
reported CLI 2.1.284. Rev9 refuses CLI versions below 2.1.281 before
registration, setup, or guided work. An independent review of rev9 then found
that the documented `ln -s` command could create a nested link inside an
existing destination. Rev10 uses an exclusive registration helper and tests
occupied destinations. The rev10 review then found that a missing prefix
followed by `..` could bypass registration's preflight and write inside the
source or through a linked parent. Rev11 rejects such paths and tests all
four reported cases. Its release packet and Desktop behavior still require
independent review and live checks. (SUPERSEDED 2026-09-29 (rev12; see the Rev12 section below): rev11 was
reviewed and its Desktop behavior trialled; rev12 still needs review.)

Keep **one procedure copy** at `cctbx_project/libtbx/guided_coding/`.
Register it once on your machine with a symlink at
`~/.claude/skills/guided_coding`. In any repository, start a guided task
with `/gc <task>`; the same personal skill should also answer to
`/guided_coding`. An ordinary Claude Code request can register that one
link from a specified reviewed source path. `/gc setup` prepares each
target's method and adoption for the Developer's decision; `/gc status`,
`/gc help`, and `/gc uninstall` provide small control requests. There is
no new libtbx dispatcher or per-project procedure install. Build commands,
`t96`, `anaconda.lbl.gov`, and `cci-gpu-00.lbl.gov` are optional project
facts, not universal global installer options.
A normal prompt in a fresh conversation stays ordinary
**once old profiles that turn GuidedCoding on for every task are retired**.
Claude Code does not automatically invoke this manual-only skill.
Before guided work in a project, its current authority must explicitly
adopt the exact contract version and identity and name its own project
method. The skill pauses and asks for that project decision if absent;
the small declaration is not a copy of the procedure. An ordinary task
does not require adoption.

The procedure still shows a plan, tested result with complete diff and exact
new test code, an integration decision, and a separate publication decision.
Its evidence and records belong to the **repository being changed**. The
skill and all procedure screens and tools live centrally. Any older
per-repository copies can be retired after their always-on startup
instructions are updated without losing local machine facts or grants.
No automatic `CLAUDE.local.md`, per-repository install, automatic lock guard,
or `/wrap` command is part of this version.

The four user documents live together in `docs/`. Read [GUIDED_CODING_USER_GUIDE.md](GUIDED_CODING_USER_GUIDE.md) for setup and normal use,
[GUIDED_CODING_ARCHITECTURE.md](GUIDED_CODING_ARCHITECTURE.md) for paths and authority, and
[GUIDED_CODING_VERIFICATION.md](GUIDED_CODING_VERIFICATION.md) for what this candidate has and has not
been tested against. The skill itself is [SKILL.md](../SKILL.md). A review
archive extracts `libtbx/guided_coding/` directly, so create an empty
staging directory first; the exact commands are in the user guide.
