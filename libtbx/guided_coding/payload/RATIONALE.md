# Why the rules are what they are

Not needed to follow the procedure. Read this when a rule looks arbitrary, or before proposing to change one. Every entry names the evidence that produced it.

**Only before/after transitions count.** Some registered tests fail on macOS while all pass on Linux. An absolute pass or fail means nothing here; a both-sides failure whose MODE changed does belong to the change.

**The executed test is always the approved bytes.** An existing test file's baseline and candidate commits need not contain the test-writer's additions. Verified by hash before each decisive run.

**A test part added during a repair gets its own discrimination run** against the pre-repair candidate - the only state it was written to catch. Otherwise it is machinery that has only ever been seen to pass.

**Locks are taken with `mkdir`.** Check-then-create is not atomic; two sessions can both pass it.

**"Still running" is a two-pass, four-way question.** Executed in the 2026-09-11 drill: a launcher exited while its job ran on, so a single-PID check reported "not running" on a live job. Pass 1 enumerates the process group by PGID alone - a plain `sleep 30` child carries no token and would be missed by a token filter. Pass 2 hunts escaped descendants by token and working directory. Read pass 1 FIRST: a present terminal marker never overrides a non-empty group, observed twice where a shell blocked waiting on an escapee stayed in the group after the job ended. And the honest limit, executed three times: a descendant that changes group, discards the token and leaves the directory is invisible to both passes. Three clean checks make a retry likely safe, never proven safe. Closing that needs containment, not better detection.

**Exit status is never a verdict.** The suite exits 0 while reporting failures - observed twice. Rosters come from log content.

**`[FAIL]` lines are printed twice per failing test** - the end-of-run replay loop in `parallel.py`, not the retry annotation. Count only the streaming section.

**There is no "timed out" outcome.** `parallel.py` has no timeout mechanism; a hung test prints no result line and surfaces as MISSING. A free-text search for "timeout" matches test FILENAMES and reports healthy tests as timed out. The roster rule once blocked on the impossible case and omitted the real one.

**A suite cannot pass by disappearing.** Newly missing, newly skipped, zero collected, or did-not-run all block like a new failure.

**A fact is not a claim.** Facts originate from the command you caused, the exit the harness observed, and records written where candidate code cannot write. Captured test output is a claim: displayed, never counted. This is provenance, not a security control - a candidate writing arbitrary bytes into the harness's stream is beyond what any reader can distinguish. (From the GuidedCoding project's P3.)

**VERIFIED versus RECORDED.** A checker tests predicates; it cannot test judgments. Writing VERIFIED against a judgment overstates what it knows.

**The shadow arrangement selects what executes; it does not prevent writes.** Symlinks are transparent. The host installation's safety is a separate, profile-recorded posture. Recursive permission changes were tried and removed: they cannot be reversed faithfully and an interruption leaves the installation altered.

**Subagents use absolute paths.** A relative path from the wrong directory wrote into the developer's real file for three minutes on 2026-09-11. Detected because the subagent disclosed it and the restore was verified by content.

**A subagent's deliverable must survive its worktree.** A test-writer's edited file was lost when its worktree was auto-cleaned; only a patch left in the scratchpad saved the work.

**The pre-harvest guard exists** because the suite harvests every `tst*` file anywhere under the module - including scratch copies under `.claude/worktrees/`, which appeared as registered tests (`test_all_parallel.py:35`, no hidden-directory exclusion).

**Procedure files are never hand-edited.** The install command regenerates the manifest and refuses to place a file over a locally modified one, so a hand edit trips a conflict-stop at the next install.

**The next change's baseline is the current integrated state.** Otherwise change two is tested against a state that no longer exists.

**Approvals must be understandable.** A chained command with variables cannot be judged, so approving it is theatre. Prompts that are procedure-mandated get grants; prompts from tools you reached for yourself get removal instead.

**Questions the developer cannot answer are defects in the question.** Most of the procedure's 2026-09-11 revisions came from the developer asking something he could not answer - which repository, whether tests had run, whether an outside reader existed. Each became a rule that keeps the question from reaching anyone again.

**The full server suite is publication-only - an intentional policy change (2026-09-11).** Earlier text dispatched a server baseline suite early in a change. It now runs only for publication. Reason: local integration does not warrant a 144-core suite, and every early-dispatch run in practice cost machine time without changing a decision. The case-specific registered tests still run locally, both sides, before every integrate decision. Recorded here so the change is not mistaken for something lost in condensation.
