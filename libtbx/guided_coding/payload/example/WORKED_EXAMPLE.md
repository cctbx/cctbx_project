# Worked example — the walkthrough script for build step 3

A representative small fix, run end to end WITHOUT a real integration. Its purpose is to exercise every stage and every file of the procedure once, on this machine, before the real trial. The example change is deliberately trivial; the procedure is what is being tested.

## The example task

The developer says: "In a small pure-Python utility module of your choosing under libtbx, one function mishandles an edge case (we will pick a real, harmless candidate together at walkthrough time — for the dry run the 'bug' can be a behavior we simply agree to call wrong, such as a formatting function's handling of an empty input). Fix it."

Choosing the target module happens at the walkthrough itself, from the live tree, so the example never encodes an invented fact about the code.

## What must be observed, stage by stage

1. **Describe → plan mode.** The Worker classifies the task as consequential (it changes committed code) and enters plan mode. Observe: no file is edited before approval.
2. **Investigate.** Explore or direct reads locate the function and its callers. Observe: reads only.
3. **Plan.** The Worker drafts the plan, self-reviews it, and presents the revision. The developer MUST edit one detail (any detail, even a wording change) before approving - the point is to exercise the editing mechanism, not to improve the plan. Observe: the plan names diagnosis, change, decisive check, file scope, repair allowance, judgment criterion.
4. **Pin the test.** The test-writer subagent runs in its worktree from the baseline, self-reviews on request, commits `tst_*.py`. The Worker then runs it in the live tree per `.claude/rules/testing.md` and shows it failing for the stated reason. Observe: the worktree was created and the Worker, not the writer, ran the demonstration; the import check printed a path inside the repository.
5. **Construct.** The implementer subagent builds in its own worktree from the same baseline, self-reviews, commits. The Worker checks `git diff --name-only` against the approved scope. Observe: scope check happens before anything else.
6. **Verify.** Full live-tree procedure: state recorded, stash if dirty, baseline fail, candidate pass, restore, porcelain match. Observe: every command lands in the change record; the tree ends exactly as it started.
7. **Review.** The reviewer subagent reports findings with file:line citations; the Worker dispositions each in writing. Observe: the reviewer edited nothing.
8. **Decide.** The Worker presents diff, comparison, findings, limitations, recommendation, and asks one question with options. The developer chooses **discard** — this is a dry run. Observe: the decision options were concrete and the detail was available on request.
9. **Wrap.** `/wrap` reports the record complete, the tree state clean (restored, nothing integrated), this change's own worktrees removed after the preservation tags, and listed.

## Exit criterion

Every observation above held, and the change record for the dry run reads as a complete, truthful account that a fresh session could resume from. Instruction gaps found along the way are fixed in the procedure files (by the developer's word — the deny rules stop the built-in editing tools from touching them (code run by an interpreter is not stopped by these rules - the manifest check is what detects that after the fact)) before build step 4, the outside review.
