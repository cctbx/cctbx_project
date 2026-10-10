# Guided Coding overview

When you ask Claude to change your code, you need more than a description of the fix. You need to understand what changed, see how it was checked and decide whether the result is ready to keep.

Guided Coding organizes that work. You and Claude agree on the task and a way to check it. Claude prepares the change, tests and a report, and helps you obtain a review from a separate conversation. You decide whether to apply the change to your working project and whether to push the commits to a remote repository.

You can use the same procedure in different Git projects. Each project keeps its own build commands, tests, settings and restrictions. The shared kit supplies the procedure and checking tools.

For installation and a first task, start with the [User Guide](GUIDED_CODING_USER_GUIDE.md).

## What a task looks like

Guided Coding uses the same role names as Guided Workflow. The Guide plans and directs the work. A Worker carries out one agreed change. The Reviewer checks it in a separate chat and is called the Outside Reviewer here. An optional Helper explains things to you. In the procedure files you are called the Developer, and the decisions are yours.

Start with a small, useful change in a project you know. Claude investigates the problem and proposes a plan that says what it will do, how it will check the result and what choices it needs from you.

Once you agree on the plan, Claude works through the change and collects the evidence. Its report includes the complete changes and the exact new tests, along with what passed, what was not checked and what remains uncertain.

For a change that needs outside review before it is applied, Claude prepares the review files for you. Use a separate chat, preferably with a different AI assistant. The reviewer examines that version of the work and sends back its findings. Claude addresses those findings before presenting the result for your decision.

A small change can use a lighter procedure only when it cannot change behavior, requirements, interfaces or an important claim, and you accept the worst possible effect stated in the plan. Its outside review can wait until publication. The number of changed lines does not determine which procedure applies.

If a report is hard to follow, give it to a [Guided Workflow Helper](https://www.thomasterwilliger.org/guided_workflow/index.html#panel-explain).

Applying the reviewed change and publishing it are separate decisions. Before publication, the procedure calls for the project's required full test suite, a comparison with a valid baseline run and outside review covering the commits to be sent. A baseline is an earlier result used for comparison. If the full test suite is not run, Claude must obtain your explicit decision to proceed without it; a general request to publish does not count. The outside review cannot be skipped in this version.

## Shared instructions, project settings

The shared kit contains instructions for Claude and small Python programs that check source files, evidence and decision records. It normally lives in `cctbx_project/libtbx/guided_coding/`, or in the `guided_coding/` folder inside the downloaded `GuidedCoding/` folder.

Making the command available creates a personal link to that kit. Setting up a project is a different step: it records how to work on that project and your decision to use the Guided Coding contract. The contract defines the responsibilities and the decisions that remain yours.

Setup reads saved project settings before asking for missing information. It shows you the proposed settings and instruction changes. Before saving, Claude must keep and check copies of the old files and the new text you approved. It must check that the project files have not changed since those copies were made. After saving, it reads the files back to check the result.

The package includes general defaults and a starting template for the project method—the document holding its commands and working environment. An optional setup card can suggest settings for your project. The kit does not supply your accounts, permissions or a personal server profile.

## Starting and returning to work

After making the command available, open Claude Code in the project you choose. First send:

```text
/guided_coding setup
```

In a fresh conversation, send the command and the task together:

```text
/guided_coding Fix the bug described in this issue.
```

Include the issue or bug report in that same message. The [User Guide](GUIDED_CODING_USER_GUIDE.md) explains the full sequence, including the decisions and outside-review handoff.

The shorter command `/gc` is also available. Using the command does not change Claude Code's permission settings. The app's Auto setting can reduce routine permission prompts, but it does not give permission to apply a change or publish it.

Help and status explain the available commands and current setup. History lists saved task summaries without starting or resuming a task. A task records the exact shared source it used, so later readers can identify its instructions and tools.

Start a guided conversation with the command. Without it, Claude does ordinary work unless saved project instructions activate Guided Coding. A new conversation does not erase those instructions or Claude’s saved memory.

## What the checks establish

The checking programs can confirm that files match a recorded version, that an evidence packet has not changed and that a decision screen contains the required information. Publication checks can compare the recorded commits and destination with Git and examine the proposed push command.

Those checks do not establish that a test checks the right thing, that a log is truthful or that a reviewer is independent. Records therefore say who chose the success criteria, wrote the tests and ran them. Another conversation can examine the evidence, but its review is an additional check rather than proof that the change is correct.

Many steps remain instructions Claude must follow. In particular, the setup save is not one indivisible operation. If a write or readback fails, the project can be partly updated. Claude must stop the affected work, describe the state and preserve the recovery copies.

The kit does not automatically reserve shared installations or coordinate parallel jobs. Auto mode runs one job at a time and takes a lock on its own test-only installation; it does not coordinate with other sessions that ignore that lock. Separate Git worktrees give tasks separate working files; they do not by themselves ensure that tests use separate installations or that shared resources are safe to use at the same time. Updating an installation also requires its own authorized plan; permission to publish is not permission to update it.

## Current status and limits

This version is a limited trial (a pilot) of the procedure, used only when you explicitly ask for it. Its release file carries the label `r10 rev18 candidate`. That is a version label used by the startup check, not proof of approval or publication.

This revision adds auto mode: Claude works through several jobs one at a time while you are away and saves each for your approval. It merges and publishes nothing. Recovery after a usage limit is manual, and the app's behaviour at a real usage limit has not been tested.

Recorded observations cover particular Mac Terminal and Claude app sessions. They include help loading and four non-interactive Terminal sessions that checked status, listed task history and included two setup trials that stopped on damaged or missing recovery records. The app observations have narrower scope, including help loading, session titles and the app-aware version check. They do not establish all setup and recovery behavior across clients or platforms.

The project's full server test suite was not run for the earlier Guided Coding-only pilot publication. The user explicitly waived that run for that publication. Local package tests are not a replacement, and the waiver is not a standing exemption.

Guided Coding has not established a failure rate, better fixes than unstructured AI use or reliable compliance by every future Claude session. The [verification document](GUIDED_CODING_VERIFICATION.md) explains the observations and their limits.

## Find the right document

| Document | What it helps you do |
| --- | --- |
| [User Guide](GUIDED_CODING_USER_GUIDE.md) | Install the kit, set up a project and work through a task. |
| [Registration and command reference](GUIDED_CODING_COMMAND_REFERENCE.md) | Find exact manual registration steps and command details. |
| [Architecture](GUIDED_CODING_ARCHITECTURE.md) | Understand where the parts live, who decides and what the tools check. |
| [Verification and limits](GUIDED_CODING_VERIFICATION.md) | Examine the checks, recorded trials, remaining uncertainty and source identities. |
| [Skill entry](../SKILL.md) | Read the instructions loaded when you invoke the command. |
| [Contract](../payload/DEVELOPER_GUIDE_CONTRACT.md) and [roles](../payload/ROLES.md) | Read the responsibilities and decisions reserved for you. |
| [Setup procedure](../payload/SETUP.md) and [general defaults](../payload/SETUP_DEFAULTS.md) | Examine settings recovery and the required save sequence. |
| [Guide](../payload/GUIDE.md) and [Worker](../payload/WORKER.md) | Read the instructions for planning and carrying out the work. |
| [Reviewer brief](../payload/OUTSIDE_REVIEWER_BRIEF.md) and [review transport](../payload/REVIEW_TRANSPORT.md) | Prepare and examine an outside-review bundle. |

## Version notes

This is the rev18 candidate: the rev17 revision with the October 6 app-aware version check, plus auto mode. Its changes and history are in [Verification and limits](GUIDED_CODING_VERIFICATION.md#historical-record).
