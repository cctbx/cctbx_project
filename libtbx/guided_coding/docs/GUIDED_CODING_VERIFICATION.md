# Guided Coding verification and limits

Guided Coding gives you checks and a record you can examine. To assess a result, you need to know what each check establishes and where its evidence stops.

There are three different questions: What does the code check? What happened when it was run? What remains uncertain? Reading a test tells you what it is intended to check. A captured run tells you what happened for a particular source and environment. Neither alone proves that a later task will behave the same way.

This document starts with those practical limits, then gives the checking steps, a guide to the source and the detailed historical record. For everyday use, read the [User Guide](GUIDED_CODING_USER_GUIDE.md).

## What is checked

The package checks its source files against a manifest—a list of files and checksums. A checksum is a fingerprint calculated from file contents. These checks detect a changed or incomplete kit; they do not establish who supplied it.

It can also freeze and verify task evidence, check decision-screen forms, prepare review bundles and compare publication records with Git. A frozen packet is evidence whose contents and inventory have been recorded for later comparison. The [architecture document](GUIDED_CODING_ARCHITECTURE.md) explains how those parts fit together.

The checks have specific limits:

| Check | What it establishes | What still needs examination |
| --- | --- | --- |
| Source verification | Listed files and the complete inventory match the supplied manifest, with a narrow exception for recognized Python cache files. | The supplier, correctness of the procedure and contents of permitted caches. |
| Claude Code version | The relevant client meets the minimum version, or the app branch reports that the version could not be checked. | Which instructions were loaded and whether future clients use the same app marker. |
| Result screen | Required information, recorded code identities and a passing recorded check are present and match. | Whether the test checks the right thing, actually ran and supports the claimed result. |
| Outside-review record | The reply identifies the same packet, covers the required decision and says to proceed. | The reviewer's independence, authenticity of the reply and meaning of any condition. |
| Review bundle | Required files were packaged with checksums in the permitted layout. | Privacy, completeness of the evidence and permission to send or publish it. |
| Publication checks | The outgoing record passes the repository checks performed; the submitted push command has the permitted form. | Your decision to publish, a later change at the remote and commands run outside the validator. |
| Task history | Saved summaries were listed without activating a task. | Truth of those summaries and any permission or resource reservation. |
| Auto-mode queue | Jobs change state only in the permitted order; one worker is registered at a time; runs have time limits and, by default, no outside network; a ready packet carries a well-formed `SCREENING.txt` and still matches its frozen manifest at approval; the start notice and morning report have the required form. | Whether the chosen tests and screening are adequate, whether the approval summaries are true, and what the app does at a usage limit. |

Setup's record, verify and apply sequence is an instruction Claude must follow. It is not a save operation enforced by the client. Verification must succeed before the apply call; a failure after a write can still leave the project partly updated.

A successful check is not permission to apply a change, publish it or update an installation. Those decisions remain separate.

## What the recorded trials support

The records describe tests on disposable configurations and scratch projects, plus some sessions using the developer's normal Mac client. They include successful source checks, command loading and specific refusals. Two Mac Terminal setup trials stopped before applying changes when a recovery record was damaged or missing.

Those observations are useful, but narrow. They do not establish all setup and recovery behavior across clients, operating systems or future sessions. The source-level tests and client observations answer different questions.

The full PHENIX server test suite was **not run** for the earlier Guided Coding-only pilot publication, and its baseline/candidate test comparison was **not produced**. The developer explicitly waived that run for that publication. Local package tests do not replace it, and the waiver does not automatically apply to another change.

We have not established that Guided Coding produces better fixes than ordinary AI coding, how often it fails, or whether every future Guide will follow the instructions reliably.

## Important remaining uncertainty

**Setup and recovery.** Earlier trials omitted complete proposals or approved recovery text. Later trials improved on those failures, but do not erase them. In one trial, the setup record was saved to a different place than proposed and did not match the displayed proposal word for word. The session disclosed this. The two damaged-record trials cover one Mac client, one mode and scratch targets; Claude still determines what happens next in other sessions.

**New conversations.** Earlier clients wrote memory outside the target project. Later sessions using that configuration cannot be described as working only from project records without checking what else was loaded. A new conversation does not necessarily begin with no prior context.

**Missing original log.** An untouched debug log from an earlier trial could not be found in the recorded locations. Transcripts and excerpts remained. The reviewer found no important behavior claim that depended solely on that missing log. We still cannot confirm the complete original or that removing private information was its only alteration. This limits how fully the record can be checked; it does not show another failure in Guided Coding’s behavior.

**Clients and platforms.** A Terminal version check says nothing about the app engine. The app branch uses markers observed on October 6, 2026; their future stability is not established. Earlier interactive Terminal observations had continuation and permission-mode qualifications. Native Windows use, a real csh session and complete setup/recovery behavior across clients have not been established. Attended use on Linux (a real `/gc` session) has not been observed; the kit's tests run there, and an attended Linux trial is a listed follow-up. A Windows skip in a test wrapper is not evidence of Windows support.

**Auto mode.** Auto-mode runs need macOS; elsewhere `run` refuses. Recovery after a usage limit is manual. Whether the app retries an interrupted turn, and what happens to a running worker then, is untested. The stop command cannot see a system program that a test starts in its own session from a parent that exits within about a second. The network block does not cover a request that a test hands to another system service. The screening record and the approval summaries are checked for form only.

**Changes outside the project.** Observations about unchanged files and recorded tool calls apply only to what was examined. They do not establish that every file, operation or login state on the computer stayed the same.

**Independence.** Separate test authorship is recorded, but is not by itself proof of independence. A repeatable test can still check the wrong claim. Records must say who chose the success criteria, wrote the test and ran it.

## Release and source checks

Before extracting an archive, compare it with the SHA-256 checksum supplied independently in its review or publication handoff. Use a trusted extraction for the next steps.

Replace the path below with the canonical source root—the absolute path after resolving symbolic links—and paste this complete block as one command. It works in bash, zsh, csh and tcsh; each line ends with `\`, which continues it on the next line:

```bash
cd /absolute/path/to/libtbx/guided_coding && \
shasum -a 256 -c SOURCE_MANIFEST.sha256 && \
env GC_PAYLOAD_ROOT="`pwd -P`/payload" python3 -I -B payload/tools/screen_check.py verify-source . && \
cat payload/RELEASE && \
grep -Eq '^guided_coding 2\.0 pilot r10 rev19 ' payload/RELEASE
```

Where `shasum` is missing (some Linux systems), replace `shasum -a 256 -c SOURCE_MANIFEST.sha256` with `sha256sum -c SOURCE_MANIFEST.sha256`; it reads the same manifest file.

Each `&&` makes the next command depend on the previous command succeeding. A failed directory change or checksum must stop the block before the kit's checker runs. A failed inventory check must stop the remaining commands. Do not append unguarded preparation or registration commands.

After the block succeeds, separately compare the manifest's own SHA-256 and the full release label with the expected handoff. The final `grep` checks that the release line begins with the `r10 rev19` label (a later release that names rev19 only in its history is refused); it does not identify an exact revision or establish approval.

Repeat the complete source check before running central tools and before consequential work, as required by `SKILL.md`. If the source changes, stop dependent work until the difference has been resolved.

The checker rejects missing or changed listed files, unexpected files, symbolic links and files with more than one hard link. It permits only a recognized installer-generated Python cache pattern. That exception checks the filename and header type, not the cached contents.

The calls use `python3 -I -B`. The `-I` option isolates imports from the current directory and `PYTHONPATH`; `-B` prevents writing caches. The documented tools run from source. These controls check consistency, not supplier identity or Claude's compliance.

## Following a claim back to the code

Consider this claim: **registration refuses an occupied destination rather than placing another link inside it.**

First read [`register_skill` in the checker](../payload/tools/screen_check.py). It constructs one `skills/guided_coding` destination, checks whether it exists or is a symbolic link, and refuses if it is occupied. It creates the link directly with `os.symlink`.

Next read `SkillRegistrationChecks` in [the checker tests](../tests/tst_screen_check_v20.py). The tests cover a second registration attempt and destinations that are directories, other links or dangling links. They check both refusal and the absence of a nested link. Other cases examine version and parent-path refusals.

Finally, to establish that those tests ran for a particular revision, examine the captured command, source identity, output and exit status. The test's source shows its intended check; it is not a record of a passing run.

### Running the package tests

Once the source guard succeeds, the five shipped test files can be run on a trusted clean copy containing only the manifest-listed source. Use a resolved temporary directory without symbolic links. On macOS the default temporary directory is under `/var`, which is a symbolic link, so the command below uses `/tmp` resolved to its real path (`/private/tmp` on macOS):

```bash
env TMPDIR="`cd /tmp && pwd -P`" python3 -I -B -m unittest discover -s tests -p 'tst_*.py' -v
```

It works in bash, zsh, csh and tcsh. Run it from the top folder of that clean copy (the folder that contains `tests`). On macOS the run takes about eight minutes. It is a separate authorized checking step, not a command to continue after a failed source guard. Report every skipped case.

The five files contain 318 test methods: 126 for the checker, 129 for the auto-mode queue helper, 13 for the review bundle, 35 for publication checks and 15 for task history. Some methods contain additional cases. The number of methods is not a count of independent claims.

The default run contacts no outside host: tests of the auto-mode network block use the documentation address 192.0.2.1, where no host exists. Three auto-mode tests would contact a real outside host or name if allowed, so they run only when the environment variable `GC_TEST_NETWORK` is set to `1`, and are otherwise reported as skipped with the reason. When enabled they contact `https://www.rcsb.org` (a DNS lookup and HTTPS), `1.1.1.1` (port 443) and the DNS resolver for `example.com`. A default run on a Mac therefore reports three skipped auto-mode cases. One control makes a single 5-second connection attempt to 192.0.2.1 outside the block, to show that the refusal inside comes from the block; if the surrounding environment already refuses it, the three block tests are skipped as well, with that reason. Other platforms can skip further cases, each with its reason.

Auto-mode runs need macOS: the queue helper's network block uses `sandbox-exec` and its search for a run's processes uses the `ps -E` environment marker. On any other system `auto_queue.py run` refuses with exit 15 before anything starts, with or without a network grant; the helper has no option, environment variable or setting that switches that off. The tests are meant to run on Linux as well. There, and on a Mac with `GC_TEST_SIMULATE_NON_MACOS=1` (a variable only the tests read), they run the helper through the test-only launcher `tests/auto_queue_standin.py`, which accepts the platform and replaces `sandbox-exec` with a stand-in that runs the command unrestricted, inside its own process only; every command the tests run through it is local. The launcher runs commands without the sandbox and without the platform refusal on any system, so it acts only on queues the kit's tests created: before it starts or records anything it refuses with exit 16 unless every queue path the helper would open lies inside one directory directly under the temporary directory that holds the marker file `GC_KIT_TEST_QUEUE_ROOT`, which the tests write there. It resolves those paths as the helper does: the home from `--home` or `GC_AUTO_HOME` (the default home is refused), every `--queue`, the queue named by the home's `ACTIVE` file and, for `init`, the records root and the lock directory. An empty or blank `--home`, `--queue`, `--records`, `--lock-dir` or `GC_AUTO_HOME` is refused too, so that the paths the launcher checks are the ones the helper opens: the helper treats a blank value (for example `--home " "`) as a relative path, falls back to another home for an empty `--home` or `GC_AUTO_HOME` (and to the `ACTIVE` queue for an empty `--queue`), resolves an empty `--records` to the current directory, and refuses an empty `--lock-dir` itself. A real queue is therefore refused. The limits: a marker file written by hand into a directory under the temporary directory is accepted; and the check resolves the queue directory, not what is inside it, so a symlink planted by hand inside a marked test queue (for example its `runs/` pointing to a real queue's `runs/`) lets a run record land in the real queue. The check guards against mistaken use, not deliberate misuse. It checks only those queue paths: other paths, such as `--log`, `--cwd`, `--file`, `--packet`, `--script`, `--grant` and `--jobs`, are not checked, so a run started through it can write its log anywhere. The following is the Linux coverage, shown by a simulated non-macOS run on a Mac and by a real Linux run of the kit suite (anaconda, October 12, 2026: 269 tests, 261 passed, 8 skipped): the tests cover queue states and transitions, init and its refusals, runs, run records and time limits, the probe, approvals and readiness, frozen packets and screening records, the worker registry, the lock, stop, reconcile, cleanup and leftover processes found by group, session or parentage, the installation check, the morning report and the auto-mode screens. Off macOS (real, or simulated with `GC_TEST_SIMULATE_NON_MACOS=1`) these tests are skipped with their reason (revision 13 corrects the earlier count of four): `test_outside_network_blocked_by_default` and `test_probe_network_blocked_by_default` (the real sandbox-exec refusal of outside connections); `test_survivor_setsid_found_by_run_id_marker` and `test_probe_surviving_killed_wrapper_is_found` (the `ps -E` marker); `test_outside_name_lookup_blocked` (the real block of an outside name lookup; off macOS a failed lookup would only show that the host was offline); `test_network_allowed_job` and `test_outside_name_lookup_allowed_job` (the allowed-job companions of the block: off macOS the stand-in never blocks, so their success says nothing about the production control); and `test_probe_timeout_with_detached_child_holding_output` (a duplicate there: the P10 variant runs the same body with the real ps). The last three network tests are also opt-in (`GC_TEST_NETWORK=1`). Off a real Mac (not under the simulation) these are also skipped: `test_module_under_sandbox_exec_fails_once` (needs `/usr/bin/sandbox-exec`), the two Mach-O normalizer tests (need `cc` producing Mach-O), and the macOS subtests of the `ps -E` snapshot tests. The tests of the refusal itself run the real helper, not the launcher, with the non-macOS platform set inside that helper's own process. A default Linux run therefore reports 12 skips: those eight, `test_module_under_sandbox_exec_fails_once`, the two Mach-O normalizer tests (`test_normalizer_build_time_masks_date_and_derived_fields`, `test_normalizer_masks_code_signature_only_at_the_end`) and `test_case_variant_of_packet_refused_on_case_insensitive_volume` (skipped where the filesystem is case-sensitive); the anaconda run of October 13, 2026 (checkpoint 3 source, 314 tests) reported exactly these 12. Two further skips depend on the host on any platform: `test_install_check_classifies_timestamp_only_archive_change` is skipped with "needs ar and ranlib on PATH" where those tools are missing, or with "this ar writes no member dates" where `ar` writes no dates even when asked (`rcU`); and the three bytecode tests of the checker are skipped where a bytecode cache prefix is set. A default Mac run reports 3 skips, the opt-in tests. Default runs contact no outside service; the opt-in tests say what they would contact.

Revision 13 adds three things to what the auto-mode helper checks. First, the process table fails closed: when `ps` cannot be read, the helper refuses with exit 17 and never prints Stopped, Cleaned or a released lock, never passes readiness and shows unknown in the report; tests replace the helper's `ps` with a missing, refused, empty or malformed stand-in in the helper's own process. Second, the kit's own suite checks `ps` first and, inside `sandbox-exec`, reports one failure naming the cause and the supported route (an attended run outside the queue's outer sandbox) instead of many unrelated failures. Third, `install check` can classify a difference as "timestamp-only differences: attended validation needed" from a `NORMALIZED` value (`payload/tools/fingerprint_normalize.py`), tested with real `ar` and `ranlib` archives and, on macOS, a `cc`-built bundle; the installation stays unavailable either way and only an attended validation records a new value. Some of the new tests need macOS; the complete off-macOS skip list is given above (it supersedes an earlier statement here of six macOS skips).

Revision 13 also changes the tests and publication tools. Timing: the suite's wall time and its 20 slowest tests were measured before and after (records of the job). Tests that depended on how fast helper subprocesses start now synchronize instead of sleeping: the stop-during-probe test holds its probe (which ignores SIGTERM) until the stop command has returned Stopping while the starting wrapper is alive (a release file, at most 30 s), the cleanup-during-probe and queue-lock tests do the same, the reused-pid test waits for a later ps start second instead of 2.1 s, and the probe time-out test measures its bound from the run record (`probe_started` to `ended`, at least 5 s and at most 7 s for a 3 s limit) rather than from the wall clock around a helper start-up, and waits at most 5 s for the killed probe group to be reaped. These were each run 20 times in a row and 10 times under induced CPU load on the Mac. What remains a limit: the real timing tests (time limits, TERM/KILL escalation, whole-second run windows) still use wall-clock bounds with stated margins, so a heavily loaded host (the Linux timing sensitivity observed on anaconda) can still fail them; a Linux 20-run repeat is the Guide's later run. Portability: the probe time-out leftover check uses separate `-o pgid= -o command=` flags (procps reads `-o pgid=,command=` as one column). Logs: `libtbx/tst_guided_coding.py` prints the SHA-256 of `SOURCE_MANIFEST.sha256`, the RELEASE line and the Git HEAD when there is one, marked `(dirty: ...)` when the package has changes not in that commit. Publication: `publication_precheck.py outgoing` generates the one-commit OUTGOING.txt block on standard output (notes about push settings go to standard error) and checks it with the same parser; the final PUBLICATION check now applies check-draft's OUTGOING.txt line rule (the seven keys only), the `SERVER_SUITE: <outcome>` first-line form and the companion refusal; a misplaced verified file is accepted from a correction companion only when its bytes equal a file of the same name listed in the frozen manifest and the screen's LIMITS names `ERRATA.md`. The fingerprint normalizer masks a Mach-O code signature only at the end of `__LINKEDIT` and of the file; any other declared range makes it compare the file byte for byte. These rules guard against mistakes, not deliberate edits: the declared ranges come from the file's own load commands, so a file edited on purpose can be given the same normalized value. A normalized match only reports that attended validation is needed; an installation is available only when its exact recorded fingerprint matches, and a new fingerprint is recorded only by an attended validation. Reading fields written as Markdown headings or table rows are compared like other decorated lines. The auto-queue helper's stop now signals every run confirmed by `ps -p` before refusing for an unreadable table; stop, cleanup, lock release and transfer read the table even with no run record; and on macOS an empty or unparsable `ps -E` output refuses.

The setup shell recipe and native Claude scenarios are not covered by those unit methods. Their controls and observations are in the task and review records.

### Checking a later revision

Keep earlier evidence unchanged and identify the new source separately. For a documentation-only change, compare the rewritten claims and command syntax with the implementation. Record which executable and contract files are unchanged. If an instruction or claim changes, assess that change; passing unit tests cannot establish its effect on Claude's behavior.

For a native-client check, identify the exact source, client and runtime, selected skill path, loaded project instructions and persistent context. Examine the actual command expansion, operations and resulting files. A version string, registration link, model-reported path or simulated procedure alone is insufficient. Matching decoded transcript text also does not prove which file bytes loaded.

Choose observations relevant to the changed behavior and intended support. Do not describe an earlier candidate's session as a fresh run of the new one.

An Outside Reviewer examines the new source and its evidence. Accepting the work, applying it, starting to use it and publishing it remain separate decisions. A final decision screen must refer to an actual review that covers its exact frozen packet. Any waiver must follow the existing contract and project rules; this document grants none.

## Historical record

The following account preserves the scope of the earlier records. It is not a new run of those sessions or a fresh inspection of live repositories.

The original documentation pass checked the uploaded baseline against its complete manifest and read the implementation and tests. Publication, client and installation facts came from the records of the session that ran and recorded the work and from outside reviews. Raw historical evidence remains in the associated task and review records, rather than being copied into the general package.

### Published pilot and setup corrections: October 2–3, 2026

The published implementation baseline was rev16, identified as `enumcheck-20261002T194333Z`. A documentation-only update followed on October 3. Its label still said `r10 rev16 candidate`; the label alone did not establish publication.

The setup correction checked the exit status of both `find` commands before using their file lists. Failures, including partial lists, stopped verification. Three test cases did not print `VERIFY_OK`; the script running the checks stopped before applying changes, and the protected files were unchanged. Seven earlier test cases retained their results, and the review finding about file listing was closed.

The session that ran and recorded the work reported 45 package tests passing on the Mac. The publication reviewer ran the repository harness on Linux and reported successful unit tests, with one case-insensitive-filesystem skip and a skipped libtbx precompile branch. A skip is not a pass for that branch.

A disposable Mac Terminal configuration using Claude Code 2.1.284 loaded help and status, ran the source guard and compared the command with its expansion. This established those loading/help/status observations, not setup or recovery.

A live Terminal session also loaded help from the expected source and ran the guard. That successful attempt followed an authentication-expiry failure and the user's sign-in.

A Claude app Code-tab session loaded help through the personal link and ran the source guard. Its transcript reported runtime 2.1.286. Its scope was loading and help, not app status or full setup/recovery.

The final publication review said to proceed with the limited pilot. It closed its four findings; the earlier file-listing finding remained closed. A description of three test assertions was corrected outside the frozen packet. The genuine publication screen passed; the record contains the developer's PUBLISH decision and verified remote pushes.

The session that ran and recorded the work reported clean Mac and shared Linux checkouts at the publication commit, matching source files and successful repository-harness checks on both. No Claude client was installed or run on that Linux installation. Matching kit files there did not establish native client support, and those checks were not a full PHENIX suite.

On October 3, the session that ran and recorded the work reported the development tag absent from the remote and its local copy removed after preservation. The published commit remained in history. That removed tag is not an installation target.

### Follow-up revision: October 4, 2026

Rev17, `followups-20261004`, changes the checker to examine review scope, outgoing-commit records and the required form of a waived test suite. It adds publication and history tools and their tests, and updates the associated procedure texts. The contract, project-adoption rule and setup-save recipe are unchanged.

The changes followed three problems in the development project: a development tag sent to a repository that reserves tags for releases, a publication based on an inferred suite waiver, and installation updates planned after the push. The tag was removed on October 3; the developer corrected the waiver on October 4. This revision adds no parallel-job coordinator.

A public API read on October 4 reported successful GitHub Actions `quick` and `clutter_and_syntax` jobs for both published commits, and a successful mirror workflow after each. It was a read of recorded results, not a rerun or a PHENIX suite.

#### Four Mac Terminal sessions

Four non-interactive Claude Code 2.1.284 sessions used Auto mode and a project-level skill link in a scratch project. Each printed the source guard's success for the worktree path and the rev17 candidate label.

In those sessions, status recovered the adoption declaration, method and records from files, and history listed one record without creating a file. A setup save with truncated approved text stopped at verification; another with no prior record or absence marker also stopped there. Apply was held, and the protected files' hashes were unchanged in every session.

The Terminal title instruction was printed in three sessions and omitted in one. These observations cover one client, one mode and scratch targets. They used the developer's live login and wrote client bookkeeping outside the targets, so they do not establish recovery without any other persistent context.

#### Claude app session titles

The Guide is the session that plans and directs the work. In its own app session, the rename tool replaced an app-generated title without a prompt. That session had loaded the published rev16 skill, rather than the candidate. It showed that the client offered the tool, not how the candidate behaved.

Two other app sessions used an exact copy of an earlier candidate, kept read-only in a scratch project. Its `SKILL.md` title instructions matched those in the recorded published revision, but later changes affected the tools, tests and some texts. These observations cover title behavior on that copy.

Both transcripts showed the skill expansion and a passing source guard. One app-generated title became `GuidedCoding: help` through the rename tool without a prompt. A title set by the developer, `KEEP MY TITLE`, was left unchanged. That second case did not exercise an approval dialog. The developer reported both titles unchanged after closing the windows.

Four earlier app sessions remain in the record but were not counted: two met a failed guard while the worktree was being edited, one renamed itself on the preceding candidate, and one met a failed guard and left a hand-set title untouched. A project-skill notice saying the name was not a command appeared each time; it establishes neither success nor failure.

#### Project-specific test discovery

A PHENIX test lookup found tests calling `run_autobuild`. Its default mode matched test names. A function newer than its May 6, 2026 static index produced no entry or explanatory message, and its changed-file mode considered only uncommitted modifications in three module directories. This was a project-specific tool, not a feature shipped in the general kit. The recorded command was `phenix.find_program search_type=tests search_text=<function> tests.search_tests_by=function_called`; the changed-file option was `git_affected_tests=True`.

### App-aware version check: October 6, 2026

This separately recorded change makes the version check use the app's own engine in a recognized app session. It updates the version tests and corresponding documentation. The release label, contract, adoption rule and setup recipe are unchanged.

In the Guide's own Mac app Code-tab session, the separately installed Terminal command was removed from the search path for each checking command only. The checker read the app engine as version 2.1.288. With the engine-path variable unset, it printed one `NOT CHECKED` line and exited successfully.

With the app marker unset and no Terminal `claude` available, it failed as a Terminal session. With the marker unset and the Terminal command available, it verified that command as version 2.1.284.

One real non-interactive Terminal run used a scrubbed environment and a one-off session-start hook supplied with `--settings`. It reported `CLAUDE_CODE_ENTRYPOINT=sdk-cli`, so it selected the Terminal branch. An interactive Terminal session was not observed.

These checks cover one Mac, one app session and one non-interactive Terminal run. Nothing was uninstalled and no settings were changed.

### Auto mode: October 8–9, 2026

The auto-mode release (release label `r10 rev18 candidate`) adds auto mode: `auto_queue.py`, `payload/AUTO.md`, the `auto_start` and `auto_report` screens and their `screen_check.py` rules, with tests, and routes the auto commands in the skill entry. Standard mode, the contract, roles and setup are unchanged. Revision 9 added the screening record required before a job is ready, the approval summaries in the morning report, manual recovery after a usage limit, a stop command for any shell, and stricter checks of job identifiers. The release-family check now expects `r10 rev18 candidate`. The checks, the hands-off runs and their limits are in the task record `phenix/.claude/records/2026-10-08-gc-auto-release-plan/`. The auto-mode tests were written from the written specification by separate test-writer subagents, which did not read the implementation; the Guide wrote the implementation and ran the checks; reviewer subagents of the same model family re-read each revision. That is separate authorship, not proof of independence.

Revision 12 (October 11, 2026) makes `run` refuse on any system other than macOS, removes the default lock directory so that `init` needs `--lock-dir`, and adds the test-only launcher so that the kit's tests can run on Linux. It followed the shared Linux suite's failures of `libtbx/tst_guided_coding.py` after the auto-mode release; real Linux verification on the server is a separate step.

### Revision 13: October 12–13, 2026

Revision 13 builds on the published auto-mode release. Its release label is `r10 rev19` (the line in `payload/RELEASE`). The work and its checks are in the task record `phenix/.claude/records/2026-10-12-gc-next-revision/`.

New commands and files:

- `screen_check.py check-draft publication DIR` checks an unfrozen publication packet and names the checks it skipped. `payload/PUBLICATION_PACKET.md` lists every file of a publication packet and its format.
- The **ACCEPT** screen (`accept`) records acceptance of an exact frozen candidate while its outside review is pending. Its DECISION must say that master and origin stay unchanged.
- `screen_check.py check` and `present` take `--corrections DIR`: a verified file frozen in the wrong place is accepted from a correction companion only when its bytes equal the file's entry in the frozen manifest, that file has the same name, one frozen file serves one name, and the screen's LIMITS names `ERRATA.md`; `present` names the file read from the companion.
- `publication_precheck.py outgoing` writes the one-commit `OUTGOING.txt` block from the live repository and checks it with the same parser.
- `auto_queue.py transfer` hands jobs over to attended work; they become **Transferred**, a final state.
- `payload/tools/fingerprint_normalize.py` and `install check --expect-normalized` recognize timestamp-only differences after a rebuild. The installation stays unavailable until an attended validation.
- Procedure texts: one definition of a clerical correction (the `ERRATA.md` companion, in `WORKER.md`); when earlier test results carry over after a restage; a plain-words summary and a test table at the top of each report; an estimate of the next update; an optional per-job staging grant and pre-batch checks (`PRECHECK.txt`) in auto mode; the attended route for running the kit's own tests.
- Documentation: separate steps for installing Claude Code, registering once and setting up a project; setup now covers the auto-mode settings; Terminal commands work in bash, zsh, csh and tcsh; `sha256sum` as an alternative to `shasum`, chosen automatically by the `/gc` source guard in `SKILL.md`; session coordination.

Behaviour changes. Some of these refuse what an earlier revision accepted:

- When the process table cannot be read, the auto-mode helper refuses with exit 17 and reports nothing as stopped, cleaned or released (`run`, `stop`, `cleanup`, `lock release`, Ready, Approved, `report` and `transfer`, also when no run record exists).
- The final PUBLICATION check now applies the draft's rules: `OUTGOING.txt` holds only the seven key lines, `SERVER_SUITE.txt` starts with `SERVER_SUITE: <outcome>`, and no companion file may be inside the packet. A packet frozen earlier that passed can now be refused; correct it by freezing a corrected packet. Published packets are not checked again. The undocumented optional `companions:` key in `OUTGOING.txt` is no longer accepted.
- In an outside reading, identical repeated `Packet identity:`, `Verdict:` and `Scope:` lines count as one. Repeats with a different value are refused, including indented, bulleted, quoted or bold ones and ones written as a Markdown heading (`## Verdict: ...`) or as a table row whose first cell is the field name, and so is a `Packet identity:` value that is not a hash.
- `worker start` is refused for a job in a final state.
- The release check (in `SKILL.md` and the documented source blocks) requires the release line to begin `guided_coding 2.0 pilot r10 rev19 `; a line that names `r10 rev19` only later, for example as a predecessor, is refused.
- Inside `sandbox-exec`, the kit's auto-queue tests report one failure, `PROCESS TABLE UNREADABLE`, naming the supported route.

## Technical reference

### Source and test map

| Claim or check | Implementation and test source |
| --- | --- |
| Complete source consistency | `verify_source`, `regular_files`; `SourceInventoryChecks`. |
| Client version and app exception | `MIN_CLAUDE_VERSION`, `APP_ENTRYPOINT`, `claude_version`, `check_claude_version`, `check_app_engine_version`; `ClaudeVersionChecks`. |
| Source and screen commands without a client on the search path | `VersionGateScopeChecks`. |
| Screen form and result identity | `blocks`, `developer_view`, `code_identity`, `check`; `ScreenChecks`. |
| Review identity, conditions and scope | `reading_matches` and the conditional, waived, pending and scope cases in `ScreenChecks`. |
| Outgoing records and suite-waiver form | `parse_outgoing`, `suite_waiver` and their `ScreenChecks` cases. |
| Accept screen and publication draft check | `accept_screen`, `check_draft`; the accept cases in `ScreenChecks` and `DraftPublicationChecks`. |
| Review-bundle layout | `make_bundle` in `review_bundle.py`; `ReviewBundleChecks`. |
| Publication repository and command checks | `publication_precheck.py`; `tst_publication_precheck_v20.py`. |
| Read-only history listing | `records_history.py`; `tst_records_history_v20.py`. |
| Auto-mode queue, grant, worker registration, runs, stop and cleanup | `auto_queue.py`; `tst_auto_queue_v20.py` and its test-only launcher `auto_queue_standin.py`. |
| Auto-mode start notice and morning report | The `auto_start` and `auto_report` rules in `screen_check.py`; `tst_screen_check_v20.py`. |
| Required setup-save ordering | `SKILL.md` and `SETUP.md` section 5, including its shell example. |

The checker requires exactly one `Scope:` value before checking it; identical repeated `Packet identity:`, `Verdict:` and `Scope:` lines count as one and conflicting repeats are refused. An outgoing repository block may not repeat a key. When the suite is marked `NOT RUN`, its quoted waiver is checked for form, not interpreted as genuine permission.

Publication queries check Git's exit status before using its output. Failed configuration or identity queries refuse, even if their output looks right. Git's documented no-match answer to the configuration query is accepted. The `vet` command accepts only the documented push form; it does not intercept other Git commands. These tests cover the demonstrated cases, rather than every possible false acceptance or every setting outside the checked Git configuration.

### Test authorship and controls

The rev17 checker, publication and history tests were written by two separate test-writer sessions from the interface specification, without reading the implementation. A control ran the final checker tests against the earlier published checker to show which tests failed there. Publication and history had no earlier tool versions for that comparison.

The October 6 version tests were written from the developer's case list by a separate AI worker started within the same conversation as the implementation. It was allowed to read the implementation for exact strings. Its earlier-version control is stored separately.

| Record | Location |
| --- | --- |
| Rev17 task and publication evidence | `phenix/.claude/records/2026-10-04-gc-followups-A/` |
| Earlier-checker control | `observations/baseline_control.txt` in that release's evidence packet |
| App-version change evidence | `phenix/.claude/records/2026-10-06-gc-app-version-check/` |
| Earlier version-checker control | `evidence-work/base_control.txt` in that change's record |

These paths identify historical records. They are not additional files included in the general download.

### Repository test wrapper

The source repository has a wrapper, [`libtbx/tst_guided_coding.py`](../../tst_guided_coding.py), which is outside this kit. From a suitable cctbx environment it can be run separately:

```bash
libtbx.python /absolute/path/to/cctbx_project/libtbx/tst_guided_coding.py
```

It copies manifest-listed files into scratch, resolves the temporary directory, runs the package tests and checks that listed installed files did not change. When libtbx precompilation is available, it checks that branch on another copy.

The wrapper skips Windows and Python 2, and may skip real precompilation. Report each skip. It reports unlisted files in the original installation without copying them into the test source; that is why the complete-source check remains a separate requirement.

Download users can run the five shipped test files as described above. They do not need a repository wrapper that their download does not contain.

### Historical source identities

The baseline identities below identify the earlier publication, not this revision or these rewritten documents.

| Item | Recorded identity |
| --- | --- |
| Pilot implementation (first publication) | r10 rev16, `enumcheck-20261002T194333Z`, October 2, 2026 |
| Published commit | `c36887c7f489018af4f91773246ecde32d1b4e24` |
| Source archive SHA-256 | `fce3627e46982f3b5fad1baaad7ada538984e9087a84b2147548b4564e844968` |
| Source manifest SHA-256 | `b86d0490df59f854c214e6771ad8bb44f03069a38612be198edff4ebb11238e0` |
| Contract version | `2026-09-17` |
| Contract SHA-256 | `ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b` |
| Publication review bundle SHA-256 | `0e4c112d3726026bd397a4e50ac2a3cff79b146883f8ceab269e3a069bc2c6ca` |
| Publication packet identity | `9834e7a5bfdf6b8dd8401972265715f527432ba277b63c4ee5c50f130508f44b` |
| Documentation-only revision | `docs-20261003`, October 3, 2026 |
| Documentation commit | `b0747a4a55f29db3abe04358480d5867e94cb792` |
| Documentation source manifest SHA-256 | `5edbec14fe45d3538da533fd995a5d59ddaa47df7629daf27321ab2a1134f1be` |

The commit, archive, manifest, bundle and packet identities of rev17 and of later revisions belong in their own task and publication records. They cannot be derived from the family label, and this document cannot include its own checksum.

Earlier reviews and tests apply to the source they name. They do not approve a later difference. The October 6 app-version change likewise has its own record, despite keeping the release label.

### Earlier development

Revisions 10–11 corrected exclusive registration and configuration paths. Revision 12 improved the instructions after a Mac trial. Revision 13 added the limited installer-cache exception and source-only loading. Revision 14 introduced recoverable setup. Failed parallel-job work from revision 15 was excluded from the revision 16 lineage.

Revision 16 corrected checksum failures that did not stop later commands, preservation of approved text, propagation of comparison failures and checked file enumeration. A failed `check && echo OK` inside a script was not enough to stop it, and unchecked file-listing commands were not sufficient. The published implementation is the final enumeration correction, not an earlier candidate. Frozen packets and errata preserve those distinctions.

### Client reference material

Anthropic's [skills documentation](https://code.claude.com/docs/en/skills) and [instruction-loading documentation](https://code.claude.com/docs/en/memory) describe the client. They are reference material, not a native Guided Coding test run.
