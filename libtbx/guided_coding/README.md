# GuidedCoding

A procedure for making a change to this codebase with Claude Code that leaves evidence behind, feels simple to use, and asks the developer only for decisions that are genuinely theirs. It is Markdown, JSON and one small Python command - no service, no daemon, no new dependency.

## What it does

A Claude Code session started in a repository where GuidedCoding is installed runs a change as a **Worker**: it investigates, plans, builds in isolated copies with fresh-context subagents, verifies in the live tree with save-and-restore, gets one fresh review, freezes an evidence packet, and stops before integration. The **Developer** decides meaning, scope, acceptance and publication. A separate chat conversation, the **Guide**, plans the project and maintains the procedure; an independent outside reader, the **Helper**, reads proof packets before anything is published. `payload/HELPER.md` defines the roles.

## Installing into a repository

1. Make sure the package is current: `git pull` in cctbx_project, then `libtbx.refresh` once (that creates the `libtbx.guided_coding` command).
2. Open a Claude Code session in the target repository and say: **install guided_coding from libtbx**. The session follows `INSTALL.md`: an interview whose answers become your personal profile, `CLAUDE.local.md`, at the repository root - shown to you before it is written. Then it runs `libtbx.guided_coding install <repository>`.
3. `libtbx.guided_coding verify <repository>` checks every installed file against its manifest.

The command refuses rather than guesses: a symbolic or hard link where it would write, a locally modified procedure file, or a change still open in that repository all stop it before it places anything. Your records, your settings and your profile are yours; `remove` deletes only the package's own files.

A committed project `CLAUDE.md` in the repository is never touched; Claude Code loads it alongside your `CLAUDE.local.md`.

## Starting a change

Open a Claude Code session in the repository and paste, as the first message:

    Read CLAUDE.local.md and follow the workflow skill for everything below.

followed by the problem in plain words. The Worker takes it from there and stops before integration. `payload/skills/workflow/SKILL.md` is the procedure; `payload/DECISIONS.md` says who answers what.

## Where things are

    payload/                  the files installed under <repository>/.claude/
    INSTALL.md                the interview a fresh session runs
    PROFILE_TEMPLATE.md       what CLAUDE.local.md is generated from
    PAYLOAD_MANIFEST.sha256   sha256 of every payload file
    ACCEPTED_PREVIOUS.sha256  every file hash ever released - what install may replace
    ../command_line/guided_coding.py   the command
    ../tst_guided_coding.py            its test (python3, git; no libtbx needed)

## Changing the procedure

Never by hand in an installed repository: the install command detects a modified file and stops. Findings go to the Guide; a revision is a commit to this package, made through a Worker change like any other, and reaches each repository by `git pull`, `libtbx.refresh`, and `libtbx.guided_coding install`. `payload/DEVELOPING_THE_PROCEDURE.md` governs.
