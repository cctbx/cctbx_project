# Installing GuidedCoding on a new machine or repository

For a Claude Code session with NO prior setup. The developer says: "install guided_coding from libtbx" (or points you at this file). You are now the installer. Work in the four-line question shape, one question at a time, plain language; every choice below is the DEVELOPER'S - never assume, never skip.

## 1. Find the package and the target

- Locate the package: $PHENIX/modules/cctbx_project/libtbx/guided_coding/ (this file lives in it, beside `payload/`). If $PHENIX is not set, ask the developer for the installation root. The command `libtbx.guided_coding` exists once `libtbx.refresh` has run after the package arrived by `git pull`; check with `libtbx.guided_coding status <repository>` - if the command is missing, ask the developer to run `libtbx.refresh` first.
- Ask: which repository is the procedure being installed into? (Full path. One install per repository; repeat the whole sequence for a second repository.)
- Check the target is a git repository and report its state (branch, clean or not). A dirty tree does not block installation - the procedure files are invisible to git - but say so.

## 2. Establish the machine facts (the step-1 check)

Run and SHOW the developer, never assume:
- Which interpreter runs this project's code (for PHENIX: libtbx.python), and where it imports the project from: `<interpreter> -c "import <package>; print(<package>.__file__)"`. If it imports from the target repository, verification will run in the live tree with save and restore; record that. If from somewhere else, record where - the verification arrangement depends on it.
- Claude Code version, git version, platform.

## 3. The interview (one question at a time, each recorded)

1. Is there a shared server the team publishes to, and what is the publish convention? (For PHENIX: git_ci - pull --rebase then push.)
2. Is there a full regression suite? If yes: the exact command, where it runs (this machine or a named server), the per-run log convention, and - if a server - the login, the quiet-machine rule, the waiting limit, and two more facts the shadow arrangement needs: the path of the project installation ON that server, and where per-change workspaces should live there (one canonical path, used for the lock too). If no: the checklist's suite lines will resolve NOT APPLICABLE, and nothing about servers is ever raised again.
3. Is there a remote verification leg (a Linux/other tree)? Host, tree path, or none.
4. Which coding standards document governs, by path or URL - or none.
5. Are there platform caveats (tests that fail on this OS but pass elsewhere)? If unknown, record unknown - the transitions-only rule protects either way.
6. Scope: which changes may this procedure make? (Default: pure-Python only; anything needing a rebuild is out of scope until the developer widens it.)
7. Helper policy: mandatory with recorded waiver (the default), or stricter.
8. Where do procedure releases come from on this machine? (Default: this package directory; releases arrive by `git pull` of cctbx_project followed by one `libtbx.refresh`, and reach each repository by `libtbx.guided_coding install <repository>`.)

## 3a. Say this before the server questions, and mean it

"You can decline any of this now and add it later, or remove it later; either change is one short session." Every profile choice is revisable in both directions, at a natural boundary, without penalty. A developer who declines the remote leg gets a complete working procedure whose checklist suite lines resolve NOT APPLICABLE by design - not a degraded one.

## 3b. Changing the profile later (either direction)

Rerun only the affected interview questions, regenerate the profile, SHOW the diff, get approval, write the new `CLAUDE.local.md`, then run `libtbx.guided_coding install <repository>` followed by `verify` - the command records the profile's new hash in the manifest; the manifest is never edited by hand. A change in flight completes under the profile it started with; the new posture takes effect at the next change. When a capability is turned off, its commissioning artifacts (the shadow record, the proofs) are KEPT, so turning it back on is re-verification against the current build rather than commissioning from nothing. Turning something off is a recorded posture change, never a failure.

## 4. Generate, install, verify

- Fill PROFILE_TEMPLATE.md (beside this file) with the answers to produce this repository's personal profile, `CLAUDE.local.md`; SHOW it to the developer and get approval before writing anything. Write it at the repository root. Never write a file named `CLAUDE.md`: that name is the project's shared instruction file, which may already exist and is not yours.
- Run `libtbx.guided_coding install <repository>`. It places the payload under `.claude/`, refuses if anything is unsafe (a symbolic or hard link where it would write, a locally modified file, an open change), writes the git exclude lines for `.claude/` and `/CLAUDE.local.md` (machine-local, never committed), records the profile's hash in `.claude/MANIFEST.sha256`, and verifies every file it wrote. Show the developer its capture.
- If the repository had an older GuidedCoding install with a hidden `CLAUDE.md` profile: a packaged profile is renamed to `CLAUDE.local.md` by the command itself; a profile the interview generated or the developer revised is not recognisable as the package's, and the developer runs `libtbx.guided_coding adopt-profile <repository>` to migrate it on their own say-so.
- Run `libtbx.guided_coding verify <repository>` and show the result.
- End with: "Your files are here:" and the full path of every installed file; then the roles pointer (.claude/HELPER.md); then: "To start your first change, open a session in this repository and paste: Read CLAUDE.local.md and follow the workflow skill for everything below - followed by the problem in plain words."

## 5. What you never do during installation

Never modify the developer's source code, never touch git history, never contact any server, and never carry another machine's facts into this profile - every fact here was established or answered above, on this machine, by this developer.
