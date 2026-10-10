# Guided Coding architecture

Guided Coding keeps the coding procedure in a shared kit. The details of each project are kept in that project. Claude Code brings the two together when you ask it to do a guided task.

That arrangement lets you use the same procedure in several Git projects, each with its own build commands, tests and working environment. A task also leaves a record of the proposed change, the checks and your decisions, so you can examine the result later.

This document explains how those pieces fit together and what the checking tools actually establish. For installation and everyday use, start with the [User Guide](GUIDED_CODING_USER_GUIDE.md).

## The three parts

| Part | What it holds | What it does |
| --- | --- | --- |
| Shared kit | The procedure, checking tools and documentation | Supplies the instructions and tools for guided work. |
| Your project | Project settings, your existing instructions and task records | Says how to work on this project and keeps the evidence from each task. |
| Claude Code | The conversation, loaded instructions and access to tools | Reads the procedure and project settings, then carries out the work. |

These parts have different responsibilities. Checking the kit's files establishes which instructions and tools are present. It does not establish that Claude loaded them correctly or followed every instruction. Similarly, a saved project setting tells Claude which test to run; the task record must show whether that test actually ran.

### The shared kit

The source normally lives in `cctbx_project/libtbx/guided_coding/`. In the website download, it is the `guided_coding/` folder inside `GuidedCoding/`. The kit is shared: a project does not need its own copy of the procedure.

Making the command available creates a symbolic link, which is a filesystem entry pointing to the kit's folder. The usual location is `~/.claude/skills/guided_coding`. If you use a separate Claude Code configuration, the link goes under that configuration instead.

The recommended command is `/guided_coding`. The skill also declares the short name `gc`, so `/gc` is available as a shorter spelling. It is configured for explicit use: Claude should not choose this skill on its own during an ordinary conversation.

### Your project

The project keeps its own instructions and a project method. The method is a document that records practical details such as working directories, build and test commands, server settings and restrictions.

If a project already has a suitable method, Guided Coding uses it at its existing location. The project's instructions identify that location. The supplied template is a starting point when one is needed, rather than a reason to create another copy.

Project instructions also record your decision to use the Guided Coding contract. The contract sets out the responsibilities and decisions reserved for you. Making the command available and agreeing to use the procedure in a project are separate steps.

A setup card can suggest settings for a particular kind of project. Its suggestions retain their source and need to be checked for relevance. A card does not give permission to work on a project or turn its defaults into instructions for every project.

### Claude Code

The coding task normally concerns the Git repository containing Claude Code's working project when you invoke the command. You can explicitly choose a different target.

Reading the shared kit does not change the target to the repository that holds the kit. Each task records both the project being changed and the version of the shared procedure being used.

Claude Code controls which project instructions it loads. Guided Coding therefore has to preserve the instructions actually used by that client, including a relevant `AGENTS.md` if the project uses one. Setup should not create a competing instruction file just because that filename appears in an example.

Claude Code may keep transcripts and memory outside the project's task records. Starting a new conversation does not necessarily remove that context. If a report says a new conversation worked only from saved project files, check what else Claude loaded.

## From a request to a decision

The Guide helps define the task and plan. A Worker carries out an agreed piece of work. Claude Code can serve as the Guide and also perform work, directly or through another conversation or a subagent. A subagent is another AI conversation given a specific assignment.

You keep the decisions about what the change means, whether it is useful, which risks to accept, whether to waive a requirement, and whether to accept the result or take a consequential action such as publication. In the procedure files, you are called the **Developer**.

The normal sequence is:

1. Confirm the shared source and the project settings needed for this task.
2. Agree on a plan and a way to check the result.
3. Make the proposed change and collect the test results.
4. Obtain outside review when the procedure requires it.
5. Present the result for your decision about applying it to the working project.
6. Prepare a separate publication decision if you want to push the commits to a remote repository.

An Outside Reviewer examines the supplied work and evidence. The reviewer does not approve a push on your behalf. An optional Helper can help you understand the material, but its opinion does not replace the required review.

The check used to support accepting the change, and its evidence, must come from a separate conversation or a subagent that did not plan the work. Give that subagent a specific, limited assignment. Rerunning a script can show that its result is repeatable, but does not by itself make the check independent.

The records say who chose the test's success criteria, who wrote the test and who ran it. If none was independent of the planning, the record must say so and the Outside Reviewer must examine whether the test checks the right thing, rather than simply rerunning it.

Whenever the Guide produces evidence for its own plan, both the packet and the review message must say so. That makes outside review especially important: a passing test may still check the wrong thing.

A subagent can provide another reading or a bounded check, but having its own context does not by itself make it independent of the Guide who briefed it. The outside review should come from a separate conversation, using a different model where possible.

In auto mode, your one advance grant replaces the plan approval and the outside review before each job is saved. Claude chooses each job's test criterion and the ticket says so. Nothing is applied to your working project or published until you approve the job and later decide to merge and publish it.

## Programs check some things; instructions govern others

Guided Coding combines instructions for Claude with small checking programs. They create the command link, check the source files and Claude Code version, freeze evidence and check decision screens. Other programs prepare review bundles, compare publication records with Git and list past tasks. Their responsibilities are explained below.

Help, status, setup, uninstall and task requests tell Claude to follow instructions in `SKILL.md`. A program does not take over those requests. Saving setup files also depends on Claude following the required steps.

A program can check that files match a recorded version or that a screen has the required information. It cannot establish that logs are truthful, that the evidence is complete or that an approval was written by the person it names. Claude still has to carry out the procedure correctly, and you still make the decisions.

## Project setup and recovery

Setup first looks for settings already recorded in the project method and referenced cards, records or backups. It asks you for the relevant gaps rather than making you supply the same information again.

A known setting and a recently verified setting are different. An old verification date does not erase a known command or directory. Setup retains the value, its source and what still needs checking.

Setup checks only what your requested task needs. Missing server settings need not block independent local work that the project permits. Default setup discovers local information and prepares an authorized settings save; it does not connect to a server or build or test the project.

Before changing the project method or its instruction file, the Guide must preserve enough information to reconstruct both the approved new text and the prior files. The save has three stages:

1. **Record.** Save the complete approved text and an exact copy of each affected prior file. Record explicitly when a file did not exist.
2. **Verify.** Read those records back and confirm that the new and prior contents can be reconstructed exactly.
3. **Apply.** Check that the current files still match the saved prior state, make only the authorized edits, and read the results back.

The records normally live under the project's `.claude/records/`, or at its designated records location. Copies held only in a conversation are not a substitute.

The example in `SETUP.md` checks the saved records in one call, then changes the project files in another. Before changing any project file, it checks that all affected files still match the saved copies. It writes files individually, so a failure can leave a partly updated project. The record must say so and the recovery copies must remain available.

The example checks each `find` command before using its file list. A partial list produced with an error cannot count as a complete inventory.

## Checking the shared source

A **manifest** lists files and their checksums. A SHA-256 checksum is a fingerprint calculated from a file's contents. Guided Coding uses the manifest's own checksum to identify the exact source for a task or review. These checks detect changes; they do not establish who supplied the files.

First, a trusted checksum command checks the listed files. Only after that succeeds may the kit's own checker run. The checker then verifies the complete inventory, rejecting missing or changed files, symbolic links, files with more than one hard link, and unexpected files, with the bytecode exception below.

Use the complete guarded command block in [Verification](GUIDED_CODING_VERIFICATION.md#release-and-source-checks). Keep its commands together so a failure stops the later commands.

The procedure is read-only during a task. Its source is checked again before consequential work and before running a central tool. If it has changed, dependent work stops until the difference is resolved.

### Python bytecode

Python can leave compiled cache files in `__pycache__` directories. The source checker permits only a recognized cache filename beside a listed Python source file, with a permitted header flag. It does not verify the contents of that cache file.

The documented tool calls use `python3 -I -B`. The `-I` option limits Python's use of local import paths and environment settings; `-B` prevents it from writing new bytecode files. The checker runs from source. The bundle and publication tools explicitly load the sibling checker from its source text.

Ordinary test imports can behave differently, so run the shipped tests on a trusted clean copy of the kit.

### The client version check

This version requires Claude Code 2.1.281 or newer. In Terminal, the check runs the `claude` command found on that session's search path. Failure to read that version stops setup or a guided task.

In a Claude app session recognized by the app marker, the check uses the app's own Claude Code engine. A separately installed Terminal command does not decide whether the app is suitable.

If the app engine's version cannot be read, the check prints one `NOT CHECKED` message and continues. If it reads a version below the minimum, it fails. Neither result establishes which instructions the app loaded. The app marker was observed on October 6, 2026; its behavior in future clients has not been established.

Help, status, history and uninstall remain available when the version check fails.

## Evidence and outside review

A **frozen packet** contains the proposed change, test results, report and other evidence, recorded in a manifest. Freezing lets later checks detect changed files. It does not make the evidence complete or its claims true. The tool refuses to replace an existing manifest; changed evidence needs a new freeze and identity.

The review bundle holds the frozen packet and, separately, the approvals record, checker notes and any corrections. Both groups have their own checksums. This preserves the proposal as reviewed; the reviewer's reply and later responses to findings stay outside the packet too.

The [review transport instructions](../payload/REVIEW_TRANSPORT.md) and [reviewer brief](../payload/OUTSIDE_REVIEWER_BRIEF.md) give the exact layout and reading order. The bundle helper checks the layout and file integrity. It does not review the files for private information or authenticate their authors.

Leave the packet and output files alone during packaging, and deliver only a completed, verified bundle. The exact layout and output location are described in the technical reference below.

### Which version the review covers

A review names the packet it examined and says whether it covers applying the change, publishing it, or both. One review can cover both decisions for the same frozen packet.

If the commit, starting version, packet or review file changes, the review must be reconsidered for that version. A short follow-up may be enough, but it must identify the new version and the decisions it covers. An earlier refusal or a made-up reply used for testing cannot count as a review that says to proceed.

For a change on the full path, outside review precedes the decision to apply it. A change can take the light path only if it cannot change behavior, requirements, interfaces or an important claim, and the plan states the worst effect if it is wrong. You must accept that limit in the plan. The light path defers outside review to publication. Size alone does not decide which path applies.

## Decision screens and publication

A decision screen is a short saved account of the proposed plan, completed work or publication. It keeps the immediate choice readable while pointing to the detailed report and evidence.

The checker verifies required headings, the final action line and the specified records. For a result, it requires a passing recorded check, a matching reference to the tested code, and a report containing the complete change diff. A **diff** shows the exact changes between versions. The report must also have the required sections for exact changes and new tests.

These checks do not run the reported test or establish who wrote it. The records and outside review still need to support those claims.

For a full result or publication, the checker requires a review that identifies the same packet, covers the appropriate decision and says to proceed. If the reviewer makes that conditional, a separate record must address the exact condition. A waived or still-pending condition must remain visible in the decision screen. Recording a person's words does not authenticate them.

Before publication, the procedure calls for the project's full test suite—the complete required set of tests—and a comparison with a valid earlier baseline run. The packet records those results, the outgoing commits and their destination. The screen checker requires the relevant files to be present and nonempty. It does not interpret the test results or perform the result screen's code-identity checks at this stage.

If a required suite was not run, the records must say so and contain your quoted decision to proceed without it. The corresponding comparison cannot be presented as a passing result. The justification for proceeding still needs examination; merely including the files does not make an unrun suite pass.

### Checking Git before a push

The publication tool compares the outgoing record with the live repository. It checks the effective destination, relevant Git settings, the recorded commit and tree, the parent commit and the remote-tracking branch. The tree identifies the files in a Git commit; the remote-tracking branch is Git's local record of the remote branch.

The optional `--fetch` step contacts the remote and refreshes local Git information. The optional `--dry-run` asks Git what a push would do without pushing. It passes only when Git exits successfully and reports exactly one permitted branch update: a fast-forward, which adds the proposed commits without replacing the branch's existing history.

Without those options, the checks do not contact the remote. The tool never pushes the commits itself.

The tool also checks a proposed push command against a narrow permitted form. It refuses other forms because Git's abbreviated options and environment settings make a list of forbidden words insufficient. This checks only the submitted command text; it neither intercepts other Git commands nor gives permission to publish.

### Showing the decision

The `present` command first checks the saved screen, then formats it for you. It hides internal identity strings from that view while leaving them in the saved evidence.

No screen command applies code or pushes a repository. Those operations remain separate steps requiring the applicable decision.

## Updates and working environments

The personal link points to the current contents of the shared kit's folder. It does not hold a fixed copy of an older release. Updating that folder can therefore affect every project using it, and must be coordinated with active tasks and the repository's integration rules.

Keep project methods, permissions and existing records separate from kit updates. Record the source identity used for each task; different procedure versions and publication decisions must not be treated as interchangeable.

This version does not install a separate procedure in each project, create profiles automatically, add an account hook or coordinate parallel jobs. Saving, applying, testing and restoring a candidate still follow the Worker's existing instructions.

A separate Git worktree gives a task separate working files. It does not by itself establish that tests import from a separate installation or that shared resources are safe to use concurrently. The history and publication-check tools do not reserve resources, schedule work or update installations.

## Technical reference

The details below are useful when examining a record, changing the implementation or checking compatibility. They do not need to be memorized to use Guided Coding.

### Where to find the files

| File or location | Purpose |
| --- | --- |
| `SKILL.md` | Entry instructions, source and target selection, control requests, and project-readiness checks. |
| `SOURCE_MANIFEST.sha256` and `payload/RELEASE` | Exact listed source contents and the readable revision label. |
| `docs/` | Reader documentation and implementation details. |
| `payload/DEVELOPER_GUIDE_CONTRACT.md` and `ROLES.md` | Your reserved decisions and the responsibilities of each role. |
| `payload/GUIDE.md` and `WORKER.md` | Task procedure and required evidence. |
| `payload/SETUP.md` and `SETUP_DEFAULTS.md` | Setup and recovery instructions, with proposed general defaults. |
| `payload/screens/` and `templates/` | Decision-screen formats, report and handoff formats, and a starting project method. |
| `payload/tools/` | Executable checking and packaging tools. |
| `auto_queue.py` and `payload/AUTO.md` | The auto-mode queue helper and its procedure. |
| The project's instruction file and method | Its decision to use the contract, working environment, commands and restrictions. |
| The project's `.claude/records/` | Setup proposals, prior file copies, task evidence and decisions. |
| Client configuration, transcripts and memory | State managed by Claude Code, which may be kept outside the project. |

### Registration and contract identity

The personal link is `skills/guided_coding` under `CLAUDE_CONFIG_DIR` when that variable is set, otherwise under `~/.claude`. The registration tool checks the source and client, then creates the link exclusively. It refuses an occupied destination. It does not set up or adopt a project.

The adopted contract for this version is dated `2026-09-17`. Its SHA-256 is:

```
ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b
```

The project's current instructions must adopt that exact contract and identify the project method before Guide and Worker instructions govern a coding task. A command invocation, an existing link or an old release label is insufficient. Control requests are handled before that adoption check; setup can propose adoption.

The app-version branch is selected by `CLAUDE_CODE_ENTRYPOINT=claude-desktop`. It reads the engine path from `CLAUDE_CODE_EXECPATH`. Other or missing entrypoint values select the Terminal branch.

`GC_PAYLOAD_ROOT` is a shell variable used in commands to locate the shared `payload/` directory. Set it to the canonical absolute path in each shell call that uses it. The Python tools do not read it as configuration; explicit command arguments identify the source or evidence to check.

### Tool commands

The tools are in the shared kit's `payload/tools/` directory. Invoke the Python tools with `python3 -I -B` as described above. In this table, capitalized words such as `SOURCE` and `DIR` stand for paths you supply; `KIND` is the type of screen.

| Command | What it does |
| --- | --- |
| `screen_check.py verify-source SOURCE` | Checks the complete source inventory and listed file contents. |
| `screen_check.py check-claude-version` | Checks the relevant Claude Code client version, with the app exception described above. |
| `screen_check.py register-skill SOURCE` | Checks the source and client, then creates the personal command link. |
| `screen_check.py freeze DIR` | Creates a manifest for the evidence files without replacing an existing one. |
| `screen_check.py verify DIR` | Checks the frozen evidence against its manifest. |
| `screen_check.py check KIND SCREEN` | Checks a saved decision screen, with evidence and review arguments when required. |
| `screen_check.py present KIND SCREEN` | Checks the saved screen and formats it for you. |
| `review_bundle.py PACKET COMPANIONS OUTPUT` | Packages the frozen packet and separate companions for outside review. |
| `publication_precheck.py check REPO OUTGOING.txt [--repository NAME] [--fetch] [--dry-run]` | Compares the publication record with Git, optionally fetching or checking a proposed push without sending commits. |
| `publication_precheck.py vet -- <proposed push command>` | Checks that the submitted push command has the permitted form. |
| `records_history.py RECORDS_DIR [--limit N]` | Lists saved task summaries without activating a task. |

Square brackets mark optional arguments. Result and publication screens also need `--evidence DIR`, and `--reading READING` and `--disposition NOTE` where required. The screen type is `plan`, `result`, `publication` or `stop`.

### Source, packet and bundle identities

| Record | What it lists | How its identity is recorded |
| --- | --- | --- |
| `SOURCE_MANIFEST.sha256` | Package files except the manifest itself, with names starting `./`. | The manifest's SHA-256 identifies the source. |
| `MANIFEST.sha256` | Frozen evidence files except the manifest itself, with relative names without that `./` prefix. | The manifest's SHA-256 identifies the packet. |
| Bundle `SHA256SUMS` | The inner packet archive and external companions. | The outer archive's SHA-256 is supplied separately. |

The evidence `verify` command prints `VERIFIED evidence`, not the packet hash. Calculate the manifest's checksum separately when recording the packet identity.

Permitted bytecode names follow the recognized `__pycache__` pattern beside a listed `.py` file. The accepted header flags are `0` or `3`; the cache contents are not authenticated.

### Screens and review records

A saved screen is limited to 28 lines and 220 words. It must use the required headings, contain no unresolved template placeholders, and end with exactly the required `ACTION` line. Plan and result screens also need a statement of why the work matters, a recommendation with a reason, and the next action.

A review has exactly one `Scope:` line:

| Screen | Accepted review scope |
| --- | --- |
| RESULT | `integration` or `integration and publication` |
| PUBLICATION | `publication` or `integration and publication` |

If a review sets conditions, a separate note must say how each was handled. This is the disposition record. Its status is `SATISFIED`, `WAIVED` or `PENDING`; the latter two must be visible in the decision screen.

The frozen review packet includes `README_FIRST.md`, `PROOF_SUMMARY.md` and the proposed report. The separate companions include `APPROVALS.md` and `CHECKER.md`, with `ERRATA.md` when needed. The bundle helper refuses nested companions directories and duplicated reserved companion filenames in the packet.

The archive must be written in the packet's own parent directory; another path to that same directory is allowed. The helper checks directory identity and keeps directories open to limit output redirection. It does not promise safety if another process moves directories while it runs.

### Publication records and commands

`OUTGOING.txt` records, for each repository, the remote name, URL, base commit, outgoing commit, its tree and the destination in the form `<commit>:refs/heads/<branch>`. This form is called a **refspec**: it names the exact commit to send and the branch to receive it. Every outgoing commit must also appear in the screen's `BATCH` section.

Publication requires nonempty `SERVER_SUITE.txt`, which records the full test-suite result or why it was not run, and `ROSTER_COMPARISON.txt`, which records the comparison with the baseline. If the suite file's first line starts with `SERVER_SUITE: NOT RUN`, it must include a `Waiver (Developer ...):` line followed by nonempty quoted lines beginning `> `. These are form checks, not authentication or interpretation of your decision.

The Git checking tool accepts:

```
publication_precheck.py check REPO OUTGOING.txt [--repository NAME] [--fetch] [--dry-run]
publication_precheck.py vet -- <proposed push command>
```

The exact permitted push form is documented in [the tool's source](../payload/tools/publication_precheck.py). It requires `--no-follow-tags`, a remote given by name rather than a URL, and a refspec naming one exact commit and branch. This command-form check does not establish that the named remote is configured. It permits `--dry-run` and `--porcelain` in the documented positions.

The history tool accepts:

```
records_history.py RECORDS_DIR [--limit N]
```

It reads `JOB_SUMMARY.txt` or `RECORD.md`, reports unreadable rows, and does not execute instructions found in a record.
