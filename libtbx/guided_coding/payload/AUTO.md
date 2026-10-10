# Auto mode — one controller, one worker, one queue

For `/gc auto <goal or job list>`. Specializes GUIDE.md and WORKER.md only where stated here;
everything else in them applies. The contract and ROLES.md govern. Standard `/gc <task>` is
unchanged.

Re-read this file and the queue record (`auto_queue.py status`) at the start of every job, after
a resume, and after the conversation has been summarized.

## What `/gc auto` authorizes

The invocation is the Developer's advance grant (contract §9), filed as `GRANT.md` in the queue
directory before any job starts:

- **Covers**, for the listed jobs only: isolated worktrees and branches from the recorded
  baseline; edits within each job's scope; local checks in the test-only installation; reviewer
  and checker subagent readings; freezing evidence; tickets; parking as Ready for approval.
- **Waives** (a recorded, user-authorized deviation, contract §13): the per-job PLAN approval and
  the full-path outside reading before parking. Both stay available at morning approval.
- **Does not cover**: required tests, new server access, accepting an unchosen risk, behaviour
  choices outside a job's text, integration into master, publication, installation updates,
  removing anything that existed before the queue.
- **Ends** when every job is parked or blocked, at a stop, or at the run boundary (default 08:00
  local the next morning). It is not renewed by a timer. Integration and publication never rely
  on it.

Re-read `GRANT.md`, check `auto_queue.py check-stop` and look for a revocation note before each
job and before each run.

## Start (target: about one minute for a prepared project)

1. The SKILL.md source and version checks.
2. Read the method: test installation path, its fingerprint script, lock directory, queue home,
   disk floor. Missing → every job that needs installation tests is Blocked, with one reason.
3. Turn the request into `JOBS.json`: id, title, scope, `requires` (required local check labels),
   `depends_on`, and `"network": true` only for a job whose grant explicitly authorizes outside
   contact. A job whose scope or checks cannot be stated from the request is Blocked with the
   needed decision; do not guess. `init` creates the queue folder itself; do not create it first.
4. Draft `screens/AUTO_START.md` into the queue's record, check it with
   `screen_check.py present auto_start FILE`, show the output once. The STOP section gives the
   absolute path of the `auto_queue.py` actually in use (the verified procedure root).
5. `auto_queue.py init …`, write `GRANT.md` first, then
   `auto_queue.py install check --script <method's fingerprint script> --expect <recorded value>`
   (unavailable → every installation-test job is Blocked with that reason). Never pass
   `--new-baseline` unattended: recording a new starting state for the test installation is an
   attended step after a deliberate refresh, and the helper accepts it only when the check matches.
6. Keep the computer awake: where the session has the app's keep-awake request, call it with
   `until: session_idle` (it holds across follow-up turns until the session is idle about five
   minutes). The app's own "Keep computer awake while Claude works" setting, when on, covers
   working turns. A closed lid still sleeps; say so in the start notice, and that a usage limit
   stops the queue and resume is manual (`/gc auto resume`). Start job 1.
   Ask no question. A reply of "stop" at any time is honoured.

## Each job, in this order

worker prepares candidate and unfrozen evidence → reviewer reads → any correction gets its checks
rerun → freeze → checker reads the frozen packet → park. A correction after the freeze needs a new
evidence directory and a new freeze; the earlier one is kept and marked superseded.

**Controller (the Guide):**

1. `reconcile`; grant not ended, no stop, lock free, no worker registered; `set JOB Preparing`.
2. Write the front door (GUIDE §4) with these additions: queue and job ids; record the criterion
   with `auto_queue.py criterion` before any test; run `check-stop` before every edit, commit and
   new subprocess and stop at once with a checkpoint if it fails; run every build or test through
   `auto_queue.py run` with `--probe` naming the tested source/build paths; do not freeze; ask the
   Developer nothing — an uncovered choice makes the job Blocked with the decision written down.
3. `install check` again (the installation must still match its fingerprint), `worker start`,
   start one Worker subagent in the background, wait for its notification,
   `worker end`. If the worker cannot be confirmed finished or cancelled, keep the job and queue
   Stopping and start no other worker.
4. Reviewer subagent on the unfrozen candidate and evidence. A defect goes back to a worker (same
   registration rule) and its checks are rerun.
5. `screen_check.py freeze`; checker subagent on the frozen packet. A needed change returns to 4
   with a new evidence directory.
6. Check delivery (GUIDE §6, not a gate): the frozen packet names the recorded candidate tree;
   every required run is current, exit 0, for this candidate and criterion; master, the live index
   and the working installations are unchanged.
7. `set JOB "Ready for approval" --packet DIR` (the helper verifies the whole frozen packet and
   records its manifest identity; it refuses unless 6 holds) or
   `set JOB Blocked --note REASON [--decision WHAT THE DEVELOPER MUST DECIDE]` (plain words, no run
ids; omit --decision when nothing needs deciding). Then the next eligible job; never wait for an
approval.

**Criterion.** The Guide derives it from the requested behaviour and the required project checks,
records it before testing, and never weakens it to obtain a pass. Changing it makes earlier runs
stale and is shown in the ticket with old and new text. The ticket states: "criterion chosen by the
Guide under the auto grant; not seen by you before the work".

**Worker** (WORKER.md, with these differences): worktree and branch from the queue baseline, or
the exact saved predecessor for a recorded dependency; `criterion`; smallest fix; commit under the
project's commit rule; `candidate`; `lock take`; check the candidate out in the test-only
installation (refresh or rebuild when compiled code changed); required checks through `run`;
return the installation to its baseline and verify its fingerprint, or mark it unavailable;
`lock release`; evidence with `CODE_IDENTITY.txt` (a `candidate_id:` line with the id that
`candidate` printed and a `tested_tree:` line for every candidate repository), `RUN_LOG.txt`, `CHANGE.diff`,
`SCREENING.txt` and `APPROVAL_REPORT.md` (criterion first); ticket draft; `JOB_SUMMARY.txt` lines. Return without
freezing. After two repair cycles with no new passing check, stop: the job is Blocked with the
finding.

**Time limits and outside contact.** Every check runs through `run` with a time limit (the job's
own where the criterion states one, else the queue default); a run that exceeds it is ended by the
helper with all its identified processes and recorded as timed out, which is a failure, never a
pass — no question to the Developer is needed for that. Before choosing tests, screen each one for
outside access (web requests, structure or sequence search, PDB or model download, any server
address) and leave such tests out unless the grant authorizes outside contact for that job;
record the result in the packet's `SCREENING.txt`, one line per test considered:
`test: NAME | outside: WHAT | included` (or `excluded`), where WHAT is exactly `none` or what was
found (anything other than exactly `none` counts as outside access). The helper refuses
Ready without it, with no included test, or with an included test that has outside access when
the job has no network grant; it does not judge whether the screening is complete. The
helper also blocks outbound connections and name lookups for every run of a job without
`"network": true`, so an unscreened test that connects directly fails rather than makes contact.
Limit: a request a test hands to another system service (for example macOS URLSession through its
own daemon) is not blocked; screening is the first defence.

**Leftover processes and choices outside discretion.** If a run leaves processes behind (the
helper reports them; `install check` refuses while they use the installation), end them with
`auto_queue.py cleanup JOB` — an ordinary step, run automatically after the job's runs and before
`lock release`, which the helper refuses while any queue-owned process is alive. Never end a
process by hand. A process using the installation that the helper cannot attribute to the queue
keeps the installation unavailable; report it, do not end it. Anything a worker or the controller does outside its
named discretion — ending a process by hand included — is recorded with
`auto_queue.py provisional JOB --note TEXT` and becomes its own question for the Developer; the
helper refuses Approved while any of the job's provisional choices is undecided.

A failed job blocks only that job. An unusable test installation blocks further installation
tests until it is repaired or recreated; the working installation is never borrowed.

## Stop

`/gc stop` here, or from any Terminal (any shell, csh included)
`python3 -I -B <procedure root>/auto_queue.py stop`. Report
**Stopped** only when the helper prints Stopped; otherwise report **Stopping** and what remains.
The helper can stop only what it tracks: the registered worker (which stops at its next
`check-stop`) and processes started through `run`, including their process group. Anything the
worker started outside `run` is not tracked; that is why `run` is mandatory.
Reservations are released only after tracked processes have ended and the installation's
recovery is verified. Saved candidates and pre-existing work are kept.

## Morning: `/gc auto status`

Run `auto_queue.py report`, draft `screens/AUTO_REPORT.md` with that table copied exactly, check
it with `screen_check.py present auto_report FILE --queue DIR`, and show the output. Its
APPROVAL SUMMARIES section (its form is checked; its truth is the Guide's) gives, for each ready job in table order, a seven-line block
with no blank lines between blocks (`Job A: title`, then `Bug:`, `Fix:`, `Test:`, `Criterion:`,
`Limits:`, `Approval:`) — a short approval summary in plain
words — the bug, the fix, the test (failed before, passes after), the criterion the Guide chose,
the material limits and any behaviour change, and that approval accepts this exact packet and
merges nothing; give each ticket path with it.
Show `auto_queue.py workers` (worker periods, overlap: none) and say plainly if a job ran out of
time because the machine slept. Ask every pending provisional choice as its own question
("Do you accept P1 of job A: …?") and record the answer with
`auto_queue.py decide JOB N --adopt|--reject --quote "<his words>"`; never infer it from an
approval of the job. The Developer may approve several jobs in one message. Record each approval verbatim with
`set JOB Approved --quote "<his words>" --manifest <identity recorded at Ready>`; the helper
verifies the packet again and refuses if it changed, so the approval binds to that exact packet
and candidate. It changes no master. Revise starts a new candidate needing approval; discard keeps the record.
An outside review is optional when parking or approving an auto candidate; before integration and
publication the standard procedure applies, including its required outside review and its own
result packet in the standard form (an auto packet naming several repositories is not reused).

## Usage limits and `/gc auto resume`

Recovery after a usage limit is manual in this release (Developer, 2026-10-09): a limit stops the
queue where it is, and the Developer resumes it with `/gc auto resume` after the reset. Nothing
inside a session runs while the account is at a limit; whether the app retries a turn by itself
(Anthropic's documentation, read 2026-10-08, describes an "Auto-continue when limits reset"
checkbox on the Code tab's session-limit card) is the app's behaviour, not this procedure's, and
is untested here, as is what happens to a running worker subagent then. Never tell the Developer
the queue will continue automatically. The start notice says so (the checker requires it to
mention the usage limit, manual resume and `/gc auto resume`), and a Waiting job's row in the
morning report gives manual resume as its next step.

Whenever the controller finds itself continuing after an interruption (a retried turn, a summarized
conversation, `/gc auto resume`): re-read this file, `GRANT.md` and `auto_queue.py status`; run
`reconcile` and `install check`. If a worker is registered, end it in the records only after
confirming that the subagent itself has finished or was cancelled (its own completion
notification or task status in this session). In a new conversation that cannot be confirmed
from records alone (whether an ended session also ends its subagents is not established): keep
the queue Stopping and report, and continue only after the Developer confirms the old session is
gone. Never register a second worker. Then redo the interrupted job from its last checkpoint; `run` itself
refuses (exit 13) while an earlier run of the job has live processes. `wait JOB --until ...`
records a pause the controller knows about (for example a reset time it can see after resuming);
a controller stopped by a limit cannot record one. `/gc auto resume` in a new conversation does
the same and continues the remaining jobs only while the grant has not ended; otherwise it
reports.

## Not in this release

There is no `auto integrate`. Integration and publication of approved candidates are
ordinary guided work: only on the Developer's later instruction naming the jobs, and the push only
on his exact PUBLISH for the combined packet, under the standard WORKER.md rules. The overnight
grant is never used for them.
