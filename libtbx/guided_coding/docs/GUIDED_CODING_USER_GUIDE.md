# Guided Coding User Guide

Guided Coding helps you make a software change with Claude and understand the result. You say what you want done. Claude investigates the problem, proposes a plan, makes the change and shows you how it was tested. You decide whether to keep it.

Start with a small task in a project you know. This guide tells you how to install Guided Coding and set up your project so Claude knows which commands to use to build and test your program. It then takes you through your first task.

## Before you begin

These instructions use the [Guided Coding download](https://www.thomasterwilliger.org/guided_coding/). Use Claude Code on your computer and a project managed with Git. You also need Python 3.10 or newer, Bash to run the shell commands, and `shasum` to check the download. In Terminal, `git --version` and `python3 --version` show which versions of Git and Python you have.

On a Mac, open the Claude app, choose the **Code** tab and select **Local**. The app includes Claude Code. If you need the app, follow Anthropic’s [desktop setup instructions](https://code.claude.com/docs/en/desktop-quickstart). Update the app if needed and sign in with an account that includes Claude Code.

You can also use Claude Code in Terminal on macOS or Linux. Follow Anthropic’s [Terminal installation](https://code.claude.com/docs/en/setup) and [sign-in and API key](https://code.claude.com/docs/en/authentication) instructions. If you use an API key, keep using the Terminal window where you set it.

Guided Coding checks for Claude Code 2.1.281 or newer. In the app, it checks the copy of Claude Code built into the app. If it cannot read that version, it prints **NOT CHECKED** and continues. You can continue, but the required version has not been confirmed. In Terminal, it checks the `claude` program that runs when you type that command.

## Using Guided Coding in the cloud or with cctbx_project or PHENIX

A cloud session cannot use the link on your computer. It can obtain its own copy from the [source repository](https://github.com/cctbx/cctbx_project/tree/master/libtbx/guided_coding); this guide does not cover cloud setup.

[cctbx_project](https://github.com/cctbx/cctbx_project) includes the Guided Coding source. If you use your own copy, follow the one-time setup below with the path to its `libtbx/guided_coding` folder instead of downloading the kit.

[PHENIX](https://www.phenix-online.org) includes Guided Coding. Run `phenix.developer` to open the instructions and default settings for using it with PHENIX.

## 1. Download and check the kit

Download [guided_coding.zip](https://www.thomasterwilliger.org/guided_coding/guided_coding.zip) and its [checksum file](https://www.thomasterwilliger.org/guided_coding/guided_coding.zip.sha256) into Downloads. The checksum file contains a fingerprint of the ZIP. The command below checks whether your downloaded ZIP has that fingerprint.

These instructions use the exact filenames shown above. If your browser changes a filename or unzips the ZIP automatically, save an unopened copy named `guided_coding.zip` before continuing.

In Terminal, run this complete command:

```bash
cd ~/Downloads && shasum -a 256 -c guided_coding.zip.sha256
```

The result should end with `guided_coding.zip: OK`. If it does not, stop and find out why before using the download. This checks the file against the checksum; it does not tell you who supplied the files.

Unzip `guided_coding.zip` in Downloads. It creates a folder called `GuidedCoding`. Keep that folder there. If you prefer another permanent location, change the path in the message in the next section. If a `GuidedCoding` folder is already there, follow “Updates and removal” below.

Keep the files inside `GuidedCoding/guided_coding` as they came in the download. Guided Coding checks them for changes, so editing them or adding files can make the check fail. Keep your project settings and task records in your project.

## 2. One-time setup making the /guided_coding command available

Open a Claude Code conversation on your computer. In the app, choose the **Code** tab, select **Local** and choose a folder. Your home folder is a suitable choice for this one-time step. In Terminal, run `claude`. Copy and send the whole message below. It gives Claude the exact setup instructions:

```text
Please set up GuidedCoding from ~/Downloads/GuidedCoding/guided_coding for this Claude Code configuration. Read docs/GUIDED_CODING_COMMAND_REFERENCE.md there. Verify the source and release, check the Claude Code version for this session, and inspect the existing personal skill. If a working shared Guided Coding registration already exists, keep it and show its resolved source. Otherwise register the central link only if the destination is unoccupied. Show what changed. Do not connect to servers, change permission settings, or start a coding task.
```

Claude checks the kit and creates a link to it. This step is called **registration**. The link lets Claude Code find the `/guided_coding` command. If a working link already exists, Claude should keep it and show you where it points. If something else is in the way, Claude should inspect it before proposing a change.

Register once for each Claude Code configuration you use. For the usual setup, that means once on your computer. Registration does not start work on a project. When it is finished, open a new conversation and enter:

```text
/guided_coding help
```

If the command is missing, see “If something goes wrong” below.

## 3. Set up your project

Open a new Claude Code conversation in the folder for the Git project you want to change. Select that folder in the app. In Terminal, use `cd` to enter the folder before running `claude`.

The download includes an optional [general setup card](https://www.thomasterwilliger.org/guided_coding/GENERAL_DEVELOPER_CARD.md). You can attach it in the app, or tell Claude in Terminal to read `~/Downloads/GuidedCoding/GENERAL_DEVELOPER_CARD.md`.

Send:

```text
/guided_coding setup
```

Claude reads your project’s documentation and saved settings. It looks for the commands to build and test your program, the folders where those commands should run, and places to keep temporary files and task records. It should ask you only for information it cannot find. Server settings are needed only if your project uses a server.

Read the proposed settings before approving them. Check that the build and test commands are right and that Claude will run them in the right folders. Claude may also propose a few lines in your project’s instructions saying that Guided Coding’s rules apply there. Read those too. It can keep an existing correct setup.

After you approve, let Claude finish saving and checking the files. It should save the old settings and the text you approved before changing the project’s instructions. If saving fails, some files may already have changed. Ask what was saved and what needs to be put back before continuing.

Setup records how to work on your project. By itself, it does not start a coding task, connect to a server or push changes. Each project has its own settings and uses the same Guided Coding kit. If you use a second copy of the same project, run setup there too.

## 4. Start your first task

Open a new Claude Code conversation in the project you just set up. Put `/guided_coding`, your task and any bug report in one message. Here is an example:

```text
/guided_coding My program crashes when it reads an empty file. Change it so that it prints "The file is empty" and stops. When the file contains data, it should work as before. Add tests for both cases.

[Add details about the program and the crash.]
```

Use your own task and details before sending the message. Keep the command and bug report together. If you send only a bug report, Claude may start ordinary coding work without loading Guided Coding.

If you want to keep the changes on your computer, include `Do not publish` in the task. Start each guided conversation with `/guided_coding`; you do not need it before every reply. `/gc` is a shorter spelling of the same command.

## 5. Read the plan and choose

The plan should explain the problem, what Claude will change, how it will check the result and what is still uncertain. Read it before letting Claude make the change.

Choose **APPROVE** when the plan describes the work you want and a useful way to check it. Choose **REVISE** to ask for a different plan. Choose **STOP** to end the attempt.

Approving the plan lets Claude do that work. Pushing changes to your repository needs a separate decision later. Claude Code may also ask permission to use a tool. Those prompts are separate from your plan and result decisions. In the app, Claude may suggest Auto mode to reduce routine tool prompts when it is available. You still decide whether to accept the plan, keep the change and publish it.

## 6. Ask for outside review when requested

Claude prepares two files for the reviewer: a message and a package containing the work to review. Use a separate chat, preferably with a different AI assistant. The assistant needs to accept file uploads, open a `.tgz` package and review code. Attach both files, then write:

```text
Please review this package.
```

Copy the reviewer’s complete response back into the Guided Coding conversation. Claude should explain what was fixed in response and what still needs your decision. If the proposed change or review package changes after review, Claude should rerun the relevant tests and ask the reviewer to check the changed version. A short follow-up review may be enough.

Outside review may wait until publication only if the change cannot affect behavior, requirements, interfaces or an important claim, and you accept the worst possible effect described in the plan. Correcting formatting without changing meaning can qualify. The reviewer recommends what to do; you decide whether to accept the work.

## 7. Decide whether to keep the result

Read what changed, which tests ran and what remains uncertain. The report should show all the file changes and the full code for any new tests. If no tests were added, it should explain how the result was checked.

Try to explain the change in your own words. Check one important statement in the report. For example, if it says a test passed, look at the test and its recorded result. If you are unsure what the report means, ask for an explanation. You can also give it to a [Guided Workflow Helper](https://www.thomasterwilliger.org/guided_workflow/index.html#panel-explain) in a separate chat.

Choose **INTEGRATE** to put the tested change into your project. Choose **REVISE** to ask for more work. Choose **DISCARD** to abandon the proposed change. The change you accept must be the one that was tested.

Before publication, the procedure calls for a full test suite and an outside review of exactly the changes to be sent. You can choose to skip the test suite by saying so explicitly for that set of changes. In this version, the outside review cannot be skipped. Claude should record your choice and which checks were not done. A **PUBLISH** decision alone does not skip those checks.

The choices are **PUBLISH** and **HOLD**. Check which commits will be sent, to which branch and to which repository. A commit is Git’s saved record of a set of changes. Updating copies of your program on other computers or servers needs a separate update plan.

## Returning to a project

Open a new Claude Code conversation in your project and send `/guided_coding` with the next task. Claude should read the saved settings instead of asking you for them again. If you change computers, shells or project folders, or a saved command fails, it may need to check the relevant settings.

These commands can help:

| Command | What it does |
| --- | --- |
| `/guided_coding help` | Explains the commands. |
| `/guided_coding status` | Shows which kit is loaded and whether it and your project are set up. |
| `/guided_coding history` | Lists the project’s recorded tasks without starting work. |
| `/guided_coding setup` | Reads saved settings and proposes any setup changes you need. |
| `/guided_coding uninstall` | Removes the link that makes this kit available to Claude Code. |
| `/guided_coding auto` and a job list | Works through several jobs one at a time while you are away and saves each for your approval. |
| `/guided_coding auto status` | Shows the report on those jobs and how to approve them. |

For ordinary Claude Code work, open a new conversation without `/guided_coding`. The project’s instructions and Claude’s saved memory may still be loaded.

## Stopping work

You can ask Claude to change direction or stop at any time. To ask it to undo this task, say:

```text
Stop this task and restore the project to its state before this task started. Preserve any work that was already present.
```

Read Claude’s report to see what it restored and what remains unfinished. Do not assume that sending the request was enough to undo the work.

## Running several jobs unattended (auto mode)

Auto mode works through a list of jobs one at a time while you are away, for example overnight. Type `/gc auto` yourself, choosing it from the command menu, and then paste the list of jobs after it. A pasted message that only begins with the command does not start it.

Claude shows one start notice and then works without asking questions. Each job is prepared in its own branch and tested in a separate test-only installation. It is then saved as **Ready for approval**, or as **Blocked** with one reason. Nothing is merged into your main branch, your working installation is not changed and nothing is published. Tests that would contact outside services are left out unless the job allows it.

Starting auto mode lets Claude choose each job's test criterion and skip the plan approval and the outside review before the job is saved. Each ticket says that Claude chose the criterion. You can still ask for an outside review before you approve.

In the morning, send `/gc auto status`. The report has one row per job. For each ready job it gives a short summary: the bug, the fix, the test, the criterion, the limits and what approving means. Approving a job accepts that exact tested change; it merges nothing. You can approve several jobs in one message. Merging and publishing approved jobs is ordinary guided work that needs your later instruction and your PUBLISH decision.

To stop, send `/gc stop`, or run the stop command shown in the start notice in any Terminal window. It works in any shell. Claude reports **Stopped** only when nothing the queue started is still running; otherwise it says what remains.

If you reach a usage limit, the queue stops and does not continue by itself. After the limit resets, send `/gc auto resume`. Keep the computer awake and the lid open: a closed lid still sleeps.

The [registration and command reference](GUIDED_CODING_COMMAND_REFERENCE.md) gives the details and limits.

## Updates and removal

Finish active Guided Coding tasks before updating the kit. Download and check the new ZIP, then unzip it into a new, permanent folder. Keep the old folder in place.

1. Ask Claude to check the new kit before changing the working link.

2. In a conversation using the old kit, send `/guided_coding uninstall` to remove its registration link.

3. Open a new conversation. Send the whole one-time setup message from step 2, replacing `~/Downloads/GuidedCoding/guided_coding` with the actual path to the new `guided_coding` folder.

4. Open another new conversation and send `/guided_coding help` to check that the new kit is available.

Keep the old folder until the new registration works.

To remove the link, use `/guided_coding uninstall`. Claude should confirm that the link points to this kit before removing it. It keeps the kit’s files, your project settings and your task records. Open a new conversation afterward.

## If something goes wrong

**The command is missing, or Claude finds the wrong kit.** Open a new conversation after registration. Ask Claude which folder the registration link points to and whether another project skill uses the same name. It should check those before changing anything.

**Claude Code is too old, or the check says NOT CHECKED.** Update the app or Terminal program you are using, then restart it. The Claude Code version inside the app is different from the app version shown in About Claude. NOT CHECKED means Guided Coding could not read the Claude Code version.

**Claude cannot find a build or test command.** Ask it to read the project’s saved settings and documentation. It should also check the shell in which the command is supposed to run. A command can work in one shell and be unavailable in another. Supply only the missing information.

**Setup fails while saving.** Stop work that uses the new settings. Keep the backup records. Ask Claude which files changed and which ones still need to be restored.

## Further details

The [documentation index](https://www.thomasterwilliger.org/guided_coding/documentation.html) links to explanations of the process, the checking tools and the checks that have been run. The [registration and command reference](https://www.thomasterwilliger.org/guided_coding/user-guide.html#one-time-registration) gives the manual registration instructions, including how to use a different configuration folder.

Guided Coding is still being developed. Tests and outside review can miss mistakes. Claude must still follow the procedure. We have not established that Guided Coding produces better fixes than other approaches.
