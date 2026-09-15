# CLAUDE.local.md profile template

Filled by the INSTALL.md interview; every <CHOICE:...> is the developer's answer, shown and approved before writing. The generic procedure never changes per site - only this profile does. This file is written as `CLAUDE.local.md` at the repository root: Claude Code's personal, per-repository instruction file, kept out of version control by `.git/info/exclude`. A committed `CLAUDE.md` in the same repository belongs to the project and is loaded alongside it; never write one.

This checkout is `<CHOICE:repository-name>` inside <CHOICE:installation-description>.

## Standing facts
- Repository: `<CHOICE:repository-path>`.
- Environment: <CHOICE:how-the-environment-is-established; for example "source <path>/build/setpaths.sh"; state whether $PHENIX and the interpreter are present in a fresh tool shell, which they are NOT on the developer's Mac>.
- Interpreter: `<CHOICE:interpreter>`. It imports <CHOICE:package-name> from <CHOICE:import-source-as-shown-by-the-check>. <CHOICE:verification-arrangement-sentence>
- Scope: <CHOICE:change-scope; default pure-Python only>.
- Platform caveats: <CHOICE:platform-caveats-or-unknown>; only before/after transitions on the same test count.
- Server suite: <CHOICE:suite-command-and-venue-or-none>; quiet rule: <CHOICE:quiet-rule-or-n/a>.
- Server installation path: <CHOICE:host-installation-path-or-n/a>; canonical workspace and lock path: <CHOICE:workspace-path-or-n/a>; shadow arrangement: <CHOICE:commissioned-and-proven-on-DATE / not commissioned / declined>; installation the shadow links into (SHADOW.md): <CHOICE:read-only execution environment / disposable test installation at PATH / working installation by explicit election, exposure stated>.
- Every line in this profile is revisable in either direction at any time (INSTALL.md section 3b); declining a capability is a recorded posture, not a failure.
- Coding standards document: <CHOICE:standards-doc-or-none>.
- Procedure source of truth: <CHOICE:release-source>.
- Publishing convention: <CHOICE:publish-convention>.
- Remote verification leg: <CHOICE:remote-leg-or-none>.
- Helper policy: <CHOICE:helper-policy; default mandatory with recorded developer waiver>.
- Guide and Helper arrangement on this installation: <CHOICE:who-the-guide-conversation-is; who-the-helper-reader-is; or none>.
- Change records live in `.claude/records/`, one file per change, from `.claude/records/TEMPLATE.md`; machine-local, never pushed. Also there: `HANDOFF.json` (the state file, from `HANDOFF_TEMPLATE.json`) and `known_failure_variation.md`.
- `.claude/` and the root `CLAUDE.local.md` are excluded from git's view by `.git/info/exclude` on this machine; neither may ever appear in `git status --porcelain` - if one does, stop.
- Whole-tree search rule: any change to a value that can reach output starts with a search of EVERY repository under <CHOICE:modules-root> (tracked files only, data-only modules skipped), results reported per repository; this search is a standing grant.
- Data files routed through git-lfs: <CHOICE:lfs-arrangement-or-none; if the binary is absent from the tool shell, say so and name the pointer-hash check as the instrument>.

## If /fix is unavailable on this surface
Begin a change with: "Read CLAUDE.local.md and follow the workflow skill for everything below", followed by the bug report.

## The consequential rule (short form)
Material uncertainty, material consequences, or expanded scope means discussion and plan-mode approval; simple, low-impact work runs immediately.
