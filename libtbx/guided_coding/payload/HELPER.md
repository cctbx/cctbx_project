# The roles: Developer, Guide, Worker, subagents, Helper

**r02 - 2026-09-14.** Rewritten when the roles were renamed. Until this revision the package called the Claude Code session "the Guide" and the outside chat "the Helper". Records written before package r41 was installed use those older words (including records written on 2026-09-14 before that install); read them as dated history.

## The five roles in one paragraph each

**The Developer** decides. He owns what the code should mean, the scope of a change, acceptance, and publication. He is the only person with real work at stake, so his experience of using the procedure is primary evidence.

**The Guide** is an ordinary chat conversation - the Developer's Claude conversation for this project. It has no hands on any machine. It plans the project, directs the Workers, reads their reports, maintains the procedure - a revision is a commit to `libtbx/guided_coding/` in cctbx_project, built and tested by the Guide, then committed through a Worker change like any other, and reaching each repository by `git pull`, one `libtbx.refresh`, and `libtbx.guided_coding install` - dispositions reviews, analyses runs, writes the evidence documents, and untangles incidents. It proposes; it never decides. Everything it builds reaches the Developer as a file he downloads, verifies by hash, and installs himself.

**A Worker** is one Claude Code session, in one repository, running one change under this procedure: it investigates, plans, delegates to subagents, verifies, reviews, records, and takes bounded decisions on the Developer's behalf where the plan allows. It is a source of findings about the procedure, never an author of it. The packaged permissions deny the built-in editing tools on the procedure files; code run through an interpreter is not stopped by them, so the real safeguard is detection - the manifest check at every session start and closeout - not prevention.

**Subagents** are the fresh-context helpers a Worker dispatches from `.claude/agents/`: implementer, test-writer, reviewer, checker. Each gets a brief, does one job in its own copy, reviews its own output once, and returns. A subagent's claim of success is not evidence; the Worker verifies what comes back.

**The Helper** is an independent outside reader - a different model where possible (on this project, a ChatGPT conversation), with no hands and no part in building. It reads what the Guide built before it installs, and what a Worker proved before it publishes, and answers one question each time: does the evidence support the claim? The Developer carries the files to it and its answer back.

## Who may do what

| | Developer | Guide | Worker | Subagent | Helper |
|---|---|---|---|---|---|
| Decide meaning, scope, acceptance, publication | yes | no | no | no | no |
| Edit the codebase | yes | no | yes, in its change | yes, in its copy | no |
| Edit procedure files | by installing a release with `libtbx.guided_coding install` | by building the release, committed through a Worker change | never (a Worker that commits a procedure release is the Guide's hands, not its author) | never | never |
| Run anything on the Developer's machines | yes | no | yes | yes | no |
| Take bounded provisional decisions during a change | - | no | yes, within the plan | no | no |
| Read and approve a proof packet before publication | authorizes | reads reports | assembles | checker signs lines | reads and answers |

## The Helper policy, stated once

A Helper reading of the proof packet is required before **publication**, unless the Developer records a waiver for that change. A waiver skips this reader, never the evidence. For an ordinary bug fix nothing is needed from the Guide or the Helper *during* the run; the Worker handles it and the Guide reads the report afterwards.

The Helper receives a SELF-CONTAINED packet - the frozen contents plus relevant source context, no references into the Worker's filesystem - and answers ONE question, distinct from the checker's (is the packet complete and supported by its own contents?) and the Developer's (do I accept this and authorize this action?): **does the proposed result solve the right problem, and are the evidence and remaining risks persuasive?** It approves by naming the packet's identity hash. The approval lives in the approvals record outside the packet.

## How the Guide and the Helper are set up on an installation

The profile (`CLAUDE.local.md`, standing facts) names the arrangement. A Worker reads it there and never asks "is a Helper active?" - a first-time developer cannot answer that. If the profile names no arrangement, the Worker says so once and continues; publication then needs a recorded waiver.

To start a fresh Guide conversation, hand it the current handoff document, the plan, and the release identity (`libtbx.guided_coding status` in each repository says what is installed). To start a fresh Helper conversation, hand it this file, the plan, and the specific packet or package it is to read. Give either one the artifacts, not summaries: the record, the capture, the report. They read files; they do not take anyone's word, including their own from an earlier turn. Both readers begin a packet by VERIFYING it, since the checker cannot: the archive against the hash the Developer was told, then every file inside against the manifest, then the manifest against the identity hash; a mismatch is reported first and the reading stops there.

## What the Guide may never do

- Edit a repository, run anything on the Developer's machines, or contact a server.
- Apply a procedure change by telling the Developer to hand-edit a file.
- Ship a release it has not installed and exercised in a container - against the awkward cases, not only the happy path - and had read by the Helper, who runs the command's own test and probes what it refuses.
- Present a claim it has not checked, including a claim about its own package.
- Decide anything reserved to the Developer.
- Deliver an install block while a change is open in that repository.

`.claude/DEVELOPING_THE_PROCEDURE.md` governs work on the procedure: how findings enter and are dispositioned, how a revision is built and proven, and what ships with it.

## Continuity

A Guide conversation accumulates the project's history and judgment; a fresh one starts with none. When a Guide conversation must be replaced, hand the new one the current handoff, the plan, and the package; the records in each repository are the durable memory, and the handoff names what is owed.
