# Auto mode — one controller, one worker, one queue

For `/gc auto <goal or job list>`. Specializes GUIDE.md and WORKER.md only where stated here;
everything else in them applies. The contract and ROLES.md govern. Standard `/gc <task>` is
unchanged.

Re-read this file and the queue record (`auto_queue.py status`) at the start of every job, after
a resume, and after the conversation has been summarized.

**Auto-mode runs need macOS.** `auto_queue.py run` blocks outside network with
`/usr/bin/sandbox-exec` and finds a run's processes by their environment marker, which it reads
on macOS only (with `ps -E` until revision 16; since then through psutil, like the whole process
table). The helper and its tools run with Python 3.10 or newer (written `libtbx.python`; see
SKILL.md); the helper also needs psutil, and without it stops with exit 18 ("requires psutil,
which <python> does not have: make a venv and install it there ..."; not `pip install --user`,
which the kit's `-I` ignores). On any
other system `run` refuses with exit 15 ("auto-mode runs need macOS (sandbox-exec and
`ps -E`)") before it reads the queue, writes a run record or starts the probe or CMD, with or
without a network grant. The helper has no option, environment variable or setting that
permits runs off macOS; the other commands work on any system. A queue whose checks must run
elsewhere is not an auto queue. The kit's tests include a test-only launcher,
`tests/auto_queue_standin.py`, that runs the helper without the sandbox and without the platform
refusal so that the portable tests can run on Linux (see "Tests off macOS" below). It is for the
kit's own tests only and refuses real queues (exit 16); a marker written by hand, or a symlink
planted by hand inside a marked test queue, can defeat that, so it guards against mistakes, not
deliberate misuse.

## What `/gc auto` authorizes

The invocation is the Developer's advance grant (contract §9), filed as `GRANT.md` in the queue
directory before any job starts:

- **Covers**, for the listed jobs only: isolated worktrees and branches from the recorded
  baseline; edits within each job's scope; local checks in the test-only installation; reviewer
  and checker subagent readings; freezing evidence; tickets; parking as Ready for approval.
- **Waives** (a recorded, user-authorized deviation, contract §13): the per-job PLAN approval and
  the full-path outside reading before parking. Both stay available at morning approval.
- **Does not cover**: required tests, server access other than a job's staging grant (below),
  accepting an unchosen risk, behaviour choices outside a job's text, integration into master,
  publication, installation updates, removing anything that existed before the queue. This
  replaces the earlier blanket exclusion of "new server access".
- **Ends** when every job is parked or blocked, at a stop, or at the run boundary (default 08:00
  local the next morning). It is not renewed by a timer. Integration and publication never rely
  on it.

Re-read `GRANT.md`, check `auto_queue.py check-stop` and look for a revocation note before each
job and before each run.

**Staging grant (optional, per job).** `GRANT.md` may give a job a staging grant naming, all
three exactly: the host, the isolated workspace root on it, and the staging command form (for
example `rsync -a --delete <worktree>/ host:/scratch/gc-stage/<job>/`). It lets that routine
staging run without a question each time. It covers only that staging, to that workspace, in
that form: not suite runs or any other command on the host, not installation updates, not
integration and not publication; each of those still needs its own authority. A job without a
staging grant has no server access. Reducing the app's permission prompts for that staging is a
separate step: a narrowly scoped settings rule for exactly that command, one rule per job, with
the worktree and workspace named literally (for example
`Bash(rsync -a --delete /Users/me/gc_worktrees/gc-q5-a/ host:/scratch/gc-stage/gc-q5-a/)`), never a
broad `Bash(ssh:*)` or `Bash(ssh host:*)` rule. A proposed rule has no wildcard that could add
options (such as `-e` or `--rsync-path`, which run commands) or a `..` path (which with `--delete`
reaches outside the workspace) (checkpoint 3; replaces an example with wildcards). A settings rule is persistent configuration: the Guide
proposes its exact text and the Developer approves it; the Guide does not write it on its own.

## Start (target: about one minute for a prepared project)

1. The SKILL.md source and version checks.
2. Read the method (the template's "Auto mode" section; `/gc setup` proposes it): test
   installation path, its fingerprint script, lock directory, queue home, disk floor.
   Missing → every job that needs installation tests is Blocked, with one reason.
3. Turn the request into `JOBS.json`: id, title, scope, `requires` (required local check labels),
   `depends_on`, and `"network": true` only for a job whose grant explicitly authorizes outside
   contact. A job whose scope or checks cannot be stated from the request is Blocked with the
   needed decision; do not guess. `init` creates the queue folder itself; do not create it first.
4. Draft `screens/AUTO_START.md` into the queue's record, check it with
   `screen_check.py present auto_start FILE`, paste the output in full into the reply once (tool
   output is not visible to the Developer). The STOP section gives the
   absolute path of the `auto_queue.py` actually in use (the verified procedure root).
   Write its stop line exactly as `screen_check.py stop-command` prints it, run like the other
   tools (revision 18: without an argument it prints the absolute interpreter path of the Python
   running it, `sys.executable`, and refuses when that Python has no psutil; revision 17: a
   path containing a space goes in single quotes, which every supported shell reads as one word).
5. `auto_queue.py init --lock-dir <method's lock directory> …` (there is no default lock: without
   `--lock-dir` init refuses with exit 2 and creates nothing), write `GRANT.md` first, then
   `auto_queue.py install check --script <method's fingerprint script> --expect <recorded value>`
   (unavailable → every installation-test job is Blocked with that reason). Never pass
   `--new-baseline` unattended: recording a new starting state for the test installation is an
   attended step after a deliberate refresh, and the helper accepts it only when the check matches.
   Pass `--expect-normalized <recorded normalized value>` too when the method records one (see
   "Test installation after a rebuild").
   Then the pre-batch checks (after init, so the jobs and the queue exist; checkpoint 3 moved
   them here from step 1). Small preparation, not a new trial: for each job confirm the baseline
   commit, that the reported failure reproduces at that baseline (or that the job is not a fix),
   the platform its checks and expected outputs need, that its checks can run through `run`
   (local, within the time limit, no outside contact unless granted), and the outside-service
   screening of the tests it will use. Write one line per job in `PRECHECK.txt` in the queue's
   record: `job: ID | baseline: COMMIT | reproduces: yes/no/not a fix | platform: WHAT |
   run: yes/no | outside: none/WHAT | start: ok/blocked REASON`. A job failing a check is set
   Blocked (`set JOB Blocked`) with that reason (and the decision needed, if any); it is not
   started.
6. Keep the computer awake: where the session has the app's keep-awake request, call it with
   `until: session_idle` (it holds across follow-up turns until the session is idle about five
   minutes). The app's own "Keep computer awake while Claude works" setting, when on, covers
   working turns. A closed lid still sleeps; say so in the start notice, and that a usage limit
   stops the queue and resume is manual (`/gc auto resume`). Start job 1.
   Ask no question. A reply of "stop" at any time is honoured.

## Each job, in this order

worker prepares candidate and unfrozen evidence → reviewer reads → any correction of substance gets
its checks rerun → freeze → checker reads the frozen packet → park. A correction of substance after
the freeze needs a new evidence directory and a new freeze; the earlier one is kept and marked
superseded. A clerical correction follows WORKER.md's one definition (ERRATA.md companion, no rerun;
a misplaced verified file through `--corrections`; another checker-refused slip refrozen with the
same bytes and a short addendum reading); this replaces the earlier rule that every correction
refreezes.

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
5. `screen_check.py freeze`; checker subagent on the frozen packet. A needed change of substance
   returns to 4 with a new evidence directory; a clerical one goes in `ERRATA.md` and the job goes on
   (replaces "every needed change returns to 4"; post-reading errata are named in the report's
   LIMITS and offered to the reviewer with the next reading, WORKER.md).
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
`SCREENING.txt` and `APPROVAL_REPORT.md` (the template's Summary first, with the criterion in it); ticket draft; `JOB_SUMMARY.txt` lines. Return without
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
keeps the installation unavailable; report it, do not end it. A zombie (a process that has
exited but that its parent has not collected, status zombie (`ps` state `Z`), on Linux a
single-threaded one;
common in containers whose process 1 never collects them) counts as gone: it keeps no run,
group, `stop`, `cleanup` or `lock release` waiting. Anything a worker or the controller does outside its
named discretion — ending a process by hand included — is recorded with
`auto_queue.py provisional JOB --note TEXT` and becomes its own question for the Developer; the
helper refuses Approved while any of the job's provisional choices is undecided.

A failed job blocks only that job. An unusable test installation blocks further installation
tests until it is repaired or recreated; the working installation is never borrowed.

**Unreadable process table.** The helper finds queue processes through psutil (revision 16;
before, `ps`). When psutil cannot list the processes, cannot read a field a decision needs
(status, parent, group, start time; AccessDenied or ZombieProcess), finds a process gone that a
re-check still finds, or leaves out the helper's own process (on macOS also when its own
environment cannot be read), every command that needs the table
refuses with exit 17 ("cannot read the process table ...; whether queue processes are alive is
unknown") and claims nothing: `run` starts nothing, `cleanup` never prints Cleaned, `lock
release` keeps the lock, Ready and Approved are refused, `stop` records the stop but prints
Stopping (never Stopped), and `report` prints its table with **Process table unreadable** and
every job that has runs as unknown, then exits 17. A run already going records the gap
(`process_table_unreadable`), and such a run never counts for readiness. Report the condition;
do not work around it. The project's fingerprint script does the same: no fingerprint when the
table cannot be read (gc_test's `fingerprint.sh` exits 4).

**Test installation after a rebuild.** A rebuild that restores the recorded commits can still
change files through known timestamp fields only (macOS: static archive member dates, including
the date ranlib writes into the index; a library that embeds its build time as text). The
project's fingerprint script may print, besides `FINGERPRINT`, a `NORMALIZED` value computed with
`payload/tools/fingerprint_normalize.py`, which ignores exactly those fields. Verified
automatically: the exact fingerprint (availability needs it, always) and, when exact differs,
whether the normalized value equals the recorded normalized value; then `install check` says
"timestamp-only differences: attended validation needed" and the installation stays unavailable
(exit 12), exactly as for any other difference. Never accepted automatically: the helper never
makes an installation available from the normalized value, and nothing replaces a recorded
fingerprint to clear a block. Attended (the Guide, holding the lock): no process from the
installation, every repository at its recorded commit with no tracked changes, the cause of each
changed file explained, compiled extensions import and smoke checks pass, then a checked save of
the new exact (and normalized) value. The busy check and the lock are unchanged.

## Transfer to attended work

When the Developer decides that a job (or the rest of the queue) continues as ordinary attended
work, record it with `auto_queue.py transfer JOB --to <successor record directory> --note TEXT`
(or `transfer --all ...` for every job that is not final). The successor directory must exist
(create the attended task record first). The job becomes **Transferred**, a final state: it leaves
Preparing, Blocked or wherever it was; any decision it had is no longer offered (the report row
says "transferred to <directory>"; the decision, packet and approval fields move into the job's
transfer record); its history, runs and evidence are kept; a queue whose remaining jobs are all
final counts as finished. The helper refuses while a worker for the job is registered (exit 9) or
a run of it may be alive (exit 13, or 17 when the process table cannot be read). A transferred job
is never reopened in the queue: the successor owns the work, and any later candidate is made and
approved there under the standard procedure. Transfer is allowed under a stop.

## Stop

`/gc stop` here, or from any Terminal (any shell, csh included)
`<interpreter> -I -B <procedure root>/auto_queue.py stop`, with both paths absolute as the
start screen prints them (the interpreter: the absolute interpreter path printed on the start
screen, the Python that runs the kit, which has psutil; never a bare `python3`, which may be
another Python without psutil). The checker accepts any absolute path to an existing executable
file named `python`, `python3`, `python3.N` or `libtbx.python` (revision 18). A path containing
a space is written in single or double quotes (`'/Users/a b/venv/bin/python3' -I -B ... stop`);
`screen_check.py stop-command` prints the line in that form, and the checker refuses the path unquoted. A path containing a
single quote or `!` cannot be written (both are refused). Inside quotes every character must be
printable (revision 18): accented letters and other scripts, as in a home folder's name, are
accepted; a tab, another control character (Ctrl-U, CR, ESC, DEL), a format character such as the
bidi overrides U+202E and U+2066 or the zero-width space U+200B, a line or paragraph separator or a
non-breaking space is refused, because a terminal or line editor can act on it, or hide it, when the
line is pasted. Report
**Stopped** only when the helper prints Stopped; otherwise report **Stopping** and what remains.
The helper can stop only what it tracks: the registered worker (which stops at its next
`check-stop`) and processes started through `run`, including their process group. Anything the
worker started outside `run` is not tracked; that is why `run` is mandatory.
Reservations are released only after tracked processes have ended and the installation's
recovery is verified. Saved candidates and pre-existing work are kept.

## Morning: `/gc auto status`

Run `auto_queue.py report`, draft `screens/AUTO_REPORT.md` with that table copied exactly, check
it with `screen_check.py present auto_report FILE --queue DIR`, and paste the output in full into the
reply, blocked jobs included; an earlier report shown again is labelled superseded. Its
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

## Running the kit's own tests

The kit's own suite (`tests/tst_*_v20.py`, and `libtbx/tst_guided_coding.py`) runs with Python
3.10 or newer that has psutil: `libtbx.python -B -m unittest discover -s tests -p 'tst_*.py' -v` from the
kit root, TMPDIR set to a resolved path. It needs a readable process table, which it reads
through psutil (revision 16): that works inside `sandbox-exec` too, where the setuid `/bin/ps`
was refused, so the suite now runs there. Inside `sandbox-exec` a second profile cannot be
applied (`sandbox_apply: Operation not permitted`), so the tests that run the helper's `run`
(whose network block is a `sandbox-exec` profile) are skipped there with that reason; the
process-table tests run. When psutil cannot
read the table, the auto-queue test module reports ONE failure (`PROCESS TABLE UNREADABLE: ...`)
naming the cause and this route instead of many unrelated failures; it is not skipped silently.
Without psutil the auto-queue tests (and the screen tests that build a queue) are skipped with
the reason, and one test checks that the helper stops with exit 18. Without git on PATH the tests
that need a repository skip with their reason (`git not on PATH: ...`); `publication_precheck.py
check` and `outgoing` then refuse (exit 2, "requires git on PATH"), and `auto_queue.py candidate`
refuses (exit 2). The supported route remains an attended run from an ordinary shell. The
helper's own restrictions are unchanged. A queue job that must run the kit suite needs that exception
stated in advance in its grant; never disable the network protection, silently or otherwise,
to make it run.

## Tests off macOS

The kit's tests (`tests/tst_*_v20.py`, and `libtbx/tst_guided_coding.py` in the shared suite) are
meant to run on Linux. There, and on a Mac when `GC_TEST_SIMULATE_NON_MACOS=1` simulates it
inside the test processes, the tests run the helper through the test-only launcher `tests/auto_queue_standin.py`:
it imports `auto_queue.py` and, in its own process only, accepts the platform and replaces
`sandbox-exec` with a stand-in that runs the command unrestricted. Every command the tests run
through it is local. The installed helper has no hook for it; on a Mac without the simulation the
tests use the real helper and the real `sandbox-exec`. Because the launcher runs commands
unrestricted on any system, it acts only on queues the tests created: before it starts or
records anything it refuses with exit 16 unless every queue path the helper would open lies inside
one directory directly under the temporary directory that holds the marker file
`GC_KIT_TEST_QUEUE_ROOT`, which the tests write there. It finds those paths as the helper does:
the home given by `--home` or `GC_AUTO_HOME` (the default home is refused), every `--queue`, the
queue named by the home's `ACTIVE` file and, for `init`, the records root and the lock directory.
An empty or blank `--home`, `--queue`, `--records`, `--lock-dir` or `GC_AUTO_HOME` is refused too,
so that the paths the launcher checks are the ones the helper opens: the helper treats a blank
value (for example `--home " "`) as a relative path, falls back to another home for an empty
`--home` or `GC_AUTO_HOME` (and to the `ACTIVE` queue for an empty `--queue`), resolves an empty
`--records` to the current directory, and refuses an empty `--lock-dir` itself. A real queue, for example under `.claude/records/` or
`~/GuidedCoding/auto`, is refused. The limits: a marker file written by hand into a directory under
the temporary directory is accepted, and a queue there can then be run; and the check resolves the
queue directory, not what is inside it, so a symlink planted by hand inside a marked test queue
(for example its `runs/` pointing to a real queue's `runs/`) lets a run record land in the real
queue. The check guards against mistaken use, not deliberate misuse. It checks only those queue
paths: other paths, such as `--log`, `--cwd`, `--file`, `--packet`, `--script`, `--grant` and
`--jobs`, are not checked, so a run started through it can write its log anywhere.

The Linux coverage below is shown by a simulated non-macOS run on a Mac and by a real Linux run
of the kit suite (anaconda, October 12, 2026: 269 tests, 261 passed, 8 skipped).

- **Runs on Linux:** queue states and transitions, init and its refusals (including
  the missing `--lock-dir`), runs and their records, time limits, the probe, stale runs, approvals and
  readiness, frozen packets and `SCREENING.txt`, the worker registry, the lock, stop, check-stop,
  reconcile, cleanup and leftover processes found by group, session or parentage, the installation
  check, the morning report and the `auto_start` and `auto_report` screens; and the refusal of
  `run` off macOS, with and without a network grant, by the real helper.
- **Needs macOS (skipped elsewhere, each with its reason):** the complete list is in
  `docs/GUIDED_CODING_VERIFICATION.md`, "Running the package tests": eight tests skipped off
  macOS (real or simulated) and one more off a real Mac. A tenth
  (`test_case_variant_of_packet_refused_on_case_insensitive_volume`) is skipped
  on a case-sensitive filesystem, not for lack of macOS; a default Linux run
  reports all 10, plus the Mach-O real-tool conformance test (and the archive one where `ar` or
  `ranlib` is missing). The archive classification test and the two Mach-O normalizer tests never
  skip (revision 15): they write their archives and bundles themselves and use no `ar`,
  `ranlib` or `cc`; only the two conformance tests use those tools, where present.
  Off macOS the helper reads no run-id markers, as before (Linux `ps` had no `-E`); on macOS an
  unreadable environment of the helper's own process makes the process table unreadable (exit
  17), and another process whose environment cannot be read shows no marker, as with `ps -E`.
- **Opt-in everywhere:** the three tests that would contact a real outside host or name run only
  with `GC_TEST_NETWORK=1`.

## Not in this release

There is no `auto integrate`. Integration and publication of approved candidates are
ordinary guided work: only on the Developer's later instruction naming the jobs, and the push only
on his exact PUBLISH for the combined packet, under the standard WORKER.md rules. The overnight
grant is never used for them.
