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

**Clients and platforms.** A Terminal version check says nothing about the app engine. The app branch uses markers observed on October 6, 2026; their future stability is not established. Earlier interactive Terminal observations had continuation and permission-mode qualifications. Native Windows use, a real csh session and complete setup/recovery behavior across clients have not been established. A Windows skip in a test wrapper is not evidence of Windows support.

**Changes outside the project.** Observations about unchanged files and recorded tool calls apply only to what was examined. They do not establish that every file, operation or login state on the computer stayed the same.

**Independence.** Separate test authorship is recorded, but is not by itself proof of independence. A repeatable test can still check the wrong claim. Records must say who chose the success criteria, wrote the test and ran it.

## Release and source checks

Before extracting an archive, compare it with the SHA-256 checksum supplied independently in its review or publication handoff. Use a trusted extraction for the next steps.

Replace the path below with the canonical source root—the absolute path after resolving symbolic links—and run this complete Bash block:

```bash
cd /absolute/path/to/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev17 candidate' payload/RELEASE
```

Each `&&` makes the next command depend on the previous command succeeding. A failed directory change or checksum must stop the block before the kit's checker runs. A failed inventory check must stop the remaining commands. Do not append unguarded preparation or registration commands.

After the block succeeds, separately compare the manifest's own SHA-256 and the full release label with the expected handoff. The final `grep` checks a version-family label; it does not identify an exact revision or establish approval.

Repeat the complete source check before running central tools and before consequential work, as required by `SKILL.md`. If the source changes, stop dependent work until the difference has been resolved.

The checker rejects missing or changed listed files, unexpected files, symbolic links and files with more than one hard link. It permits only a recognized installer-generated Python cache pattern. That exception checks the filename and header type, not the cached contents.

The calls use `python3 -I -B`. The `-I` option isolates imports from the current directory and `PYTHONPATH`; `-B` prevents writing caches. The documented tools run from source. These controls check consistency, not supplier identity or Claude's compliance.

## Following a claim back to the code

Consider this claim: **registration refuses an occupied destination rather than placing another link inside it.**

First read [`register_skill` in the checker](../payload/tools/screen_check.py). It constructs one `skills/guided_coding` destination, checks whether it exists or is a symbolic link, and refuses if it is occupied. It creates the link directly with `os.symlink`.

Next read `SkillRegistrationChecks` in [the checker tests](../tests/tst_screen_check_v20.py). The tests cover a second registration attempt and destinations that are directories, other links or dangling links. They check both refusal and the absence of a nested link. Other cases examine version and parent-path refusals.

Finally, to establish that those tests ran for a particular revision, examine the captured command, source identity, output and exit status. The test's source shows its intended check; it is not a record of a passing run.

### Running the package tests

Once the source guard succeeds, the four shipped test files can be run on a trusted clean copy containing only the manifest-listed source. Use a resolved temporary directory without symbolic links:

```bash
python3 -I -B -m unittest discover -s tests -p 'tst_*.py' -v
```

Run this from the top folder of that clean copy. It is a separate authorized checking step, not a command to continue after a failed source guard. Report every skipped case.

The four files contain 118 test methods: 63 for the checker, 13 for the review bundle, 27 for publication checks and 15 for task history. Some methods contain additional cases. The number of methods is not a count of independent claims.

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
| Review-bundle layout | `make_bundle` in `review_bundle.py`; `ReviewBundleChecks`. |
| Publication repository and command checks | `publication_precheck.py`; `tst_publication_precheck_v20.py`. |
| Read-only history listing | `records_history.py`; `tst_records_history_v20.py`. |
| Required setup-save ordering | `SKILL.md` and `SETUP.md` section 5, including its shell example. |

The checker requires exactly one `Scope:` line before checking its value. An outgoing repository block may not repeat a key. When the suite is marked `NOT RUN`, its quoted waiver is checked for form, not interpreted as genuine permission.

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

Download users can run the four shipped test files as described above. They do not need a repository wrapper that their download does not contain.

### Historical source identities

The baseline identities below identify the earlier publication, not the current candidate or these rewritten documents.

| Item | Recorded identity |
| --- | --- |
| Published implementation | r10 rev16, `enumcheck-20261002T194333Z`, October 2, 2026 |
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

The rev17 candidate's commit, archive, manifest, bundle and packet identities belong in its own task and publication records. They cannot be derived from the family label, and this document cannot include its own checksum.

Earlier reviews and tests apply to the source they name. They do not approve a later difference. The October 6 app-version change likewise has its own record, despite keeping the release label.

### Earlier development

Revisions 10–11 corrected exclusive registration and configuration paths. Revision 12 improved the instructions after a Mac trial. Revision 13 added the limited installer-cache exception and source-only loading. Revision 14 introduced recoverable setup. Failed parallel-job work from revision 15 was excluded from the revision 16 lineage.

Revision 16 corrected checksum failures that did not stop later commands, preservation of approved text, propagation of comparison failures and checked file enumeration. A failed `check && echo OK` inside a script was not enough to stop it, and unchecked file-listing commands were not sufficient. The published implementation is the final enumeration correction, not an earlier candidate. Frozen packets and errata preserve those distinctions.

### Client reference material

Anthropic's [skills documentation](https://code.claude.com/docs/en/skills) and [instruction-loading documentation](https://code.claude.com/docs/en/memory) describe the client. They are reference material, not a native Guided Coding test run.
