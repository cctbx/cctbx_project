# How the procedure itself is developed

**r01 - 2026-09-11.** This governs work on GuidedCoding, not work on the codebase. The procedure for fixing a bug is `.claude/skills/workflow/SKILL.md`; this is how that file and its companions come to say what they say. Written because every rule here was real practice for two days before anyone wrote it down, and practice that is not written down is lost when a conversation ends.

## 1. Roles and authority

**The Developer** owns the procedure. Every change to it is his decision; he is also its only user with real work at stake, so his experience of using it is primary evidence, not anecdote.

**The Guide** (an ordinary chat conversation, no hands on any machine) drafts changes, builds and tests packages, dispositions reviews, analyses runs, and writes the evidence documents. It proposes; it never decides. It cannot touch a repository, so every change it makes reaches the Developer as a file he downloads, verifies by hash, and installs himself.

**The Helper** (an independent outside reader - a different model where possible, with no hands and no part in building) reads what the Guide built before it installs and what a Worker proved before it publishes, and answers one question each time: does the evidence support the claim? It finds; it does not build or decide.

**A Worker** (a Claude Code session running a real change) is a *source of findings*, never an author of the procedure. It records what it hit and hands it over. It may not edit a procedure file; the packaged permissions deny the built-in editing tools on those files, interpreter-run code is not stopped by them, and the manifest check at session start and closeout is what detects a change.

**Outside reviewers** - other models, other people - are commissioned by the Developer to read a package and report. They have no authority either; their findings are dispositioned like any other.

## 2. The standing constraint

Prompts, Markdown, and configuration - plus exactly one small program, the `libtbx.guided_coding` command that installs, verifies and removes them, with its own test. No orchestration program beyond that, no harness, no daemon. When a problem can only be solved by writing software that would drive the procedure, the honest answer is that it is out of scope, not a smaller program pretending otherwise. Tools built *by* the procedure for its own use (the roster parser) are changes to the codebase, run through the ordinary workflow, not procedure text.

## 3. One revision at a time

A revision is the unit. It gathers every finding that has accumulated since the last one, is built and tested by the Guide, read by the Helper, committed to `libtbx/guided_coding/` in cctbx_project through an ordinary Worker change, and installed into each repository by the developer with `libtbx.guided_coding install` between changes - a Worker that finds a newer release at stage 0 asks whether to upgrade now or proceed on the installed one, and never upgrades itself. Between revisions the procedure on disk does not move.

Two rules make this safe:

- **Never during a change.** A running Worker has the old text in its context while its subagents would read the new text from disk; half the change would follow each. A change completes on the revision it started with.
- **Never by hand.** The install command regenerates `MANIFEST.sha256` and refuses to place files over a locally modified one, so a hand edit trips a conflict-stop at the next install. A Worker that suggests hand-editing a procedure file is wrong, and this has happened.

## 4. Where findings come from, and what happens to each

Four sources, in rough order of value:

1. **The Developer using it.** A question he could not answer, a prompt he could not judge, an output he could not follow. These are the most valuable findings in the project and produce the rules that matter most.
2. **A Worker's own report.** Procedure findings recorded during a change, dispositioned afterward.
3. **Outside review of a package.** Read cold, against the files, by someone with no stake.
4. **The Guide's own checks** when building or testing a revision.

Every finding gets one of: **accepted** (and appears as operative text), **accepted with change** (say what changed and why), **rejected** (say why, in writing), or **deferred** (to a named ladder entry, not to memory). A finding that is neither in the package nor in `LADDER.md` is not tracked, and saying so out loud is part of the job.

## 5. How a revision is built and proven

1. Edit `payload/` in a working copy of `libtbx/guided_coding/`. There is ONE payload for every repository; the per-site profile is the developer's `CLAUDE.local.md`, never shipped.
2. Regenerate `PAYLOAD_MANIFEST.sha256` with bytecode writing off, and set `payload/RELEASE`.
3. Append the previous release's payload manifest to `ACCEPTED_PREVIOUS.sha256` - **every** prior release, never a hand-picked tail. A hand-maintained tail once produced a false CONFLICT that cost two rounds and an unfair suspicion of a Worker. The command decides what it may replace, retire or migrate ONLY from this file and the payload manifest; the installed manifest in a repository is editable and never proves ownership.
4. **Test in a container before delivering**: run `libtbx/tst_guided_coding.py` (it installs into scratch repositories, upgrades from the previous release, and exercises every refusal), then the Helper reads the kit and runs the command's own test and its probes.
5. Commit through a Worker change in cctbx_project, decisive check = that test, published through the full gate.
6. Each repository upgrades with `libtbx.guided_coding install <repository>` after `git pull` and one `libtbx.refresh`; the command refuses while a change is open there.

The Guide checks what it ships as strictly as it asks a Worker to check a fix. Lessons paid for on 2026-09-14, each now a rule: split every long checksum across two short lines in a paste block, so a wrapped line cannot break it; never give `mktemp` a suffix after the X's on macOS - feed a helper program to Python on stdin rather than through a temporary file; and test every delivered script against the awkward cases before it goes out - a leftover file with the same name, running from the Downloads folder, and both the has-data and no-data branches - not only the happy path. A delivered script that fails on the developer's machine for a reason the container could have shown is a Guide defect, recorded as one. It has failed this twice - a false conflict list, and its own container's bytecode shipped inside two packages - and both are in `RATIONALE.md` rather than forgotten.

## 6. What the Guide may never do

- Edit a repository, run anything on the Developer's machines, or contact a server.
- Apply a procedure change by instructing the Developer to hand-edit a file.
- Ship a package it has not installed and exercised in a container.
- Present a claim it has not checked - including a claim about its own package.
- Decide anything reserved to the Developer: scope, meaning, publication, or whether a trade-off is worth it.
- Quietly improve the Developer's stated reasoning. Corrections sit beside his words, never inside them.

## 7. Adding rules, and not adding them

The procedure grows only from evidence. A rule needs an incident, a review finding, or a measured failure behind it; "it would be safer if" is not enough. Every new rule states what it prevents.

Against growth, three tests applied at every revision:

- **The promise test.** Each retained paragraph serves a stated guarantee, the current trial commitment, or plain legibility. Anything else moves to `RATIONALE.md` or goes.
- **One authoritative home.** A rule lives in exactly one file; others point at it. Duplications are where contradictions breed - one such pair made the roster rule block on a case that could not happen while ignoring the case that did.
- **Rule and reason are separated.** Operative files say what to do; `RATIONALE.md` says why. A reader must be able to tell which sentences bind.

When a rule is removed or changed deliberately - not lost in editing - say so explicitly in `RATIONALE.md`. Condensation and policy change must never be confused.

## 8. How it is measured

Word count is not the measure. The measures are:

- **Questions the Developer could not answer** during a change. The target is zero.
- **Prompts that reached him** that were not decisions. Also zero; a routine prompt is a defect in how a command was constructed.
- **Time waiting for the Developer** versus wall-clock and machine-working time, recorded per change (attention minutes cannot be observed and are not recorded).
- **Defects caught** that a plain session missed, and defects the procedure's own review caught in its output.
- Eventually: **whether a second developer can install it and run a change** without help. That is the only test of the seventh goal, and it has not been run.

## 9. When to stop

Stop adding when the last three changes produced no findings of the Developer-couldn't-answer kind. Stop a revision when its findings are dispositioned, not when the text feels finished. And stop the project's expansion where it stands: the extensions reserved and not built are reserved deliberately, and each would need its own evidence to earn a place.
