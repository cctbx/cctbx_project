# Change record: CHANGEID — one-line title

Dates: opened YYYY-MM-DD, closed YYYY-MM-DD
(The status lives in exactly one place: the generated Status section below.)
(For a PUBLICATION BATCH: the Plan section lists the constituent changes - id, record, rebased commit, evidence archive or packet with hash - and the roster prediction; the Identities section names baseline = shared master and candidate = rebased local master; construction and pinned-test sections say "not applicable - see the constituent records"; the Review section says the same ONLY if every constituent is a full-path change with its own fresh review - otherwise it holds the batch's light-constituent review of EVERY light-path constituent, findings verbatim and dispositions by identifier; the Verification section holds the suite comparison.)

## Plan as approved
(The plan text the developer approved, verbatim: diagnosis, smallest coherent change, decisive check, files in scope, repair allowance, how the result is judged. Note any developer edits to the draft.)

## Identities
- NAMED BASELINE commit (from the approved plan):
- Subagent actual base commits (per subagent, from the conformance check) and transport commits if any, with reasons:
- Pinned test: path, commit, and whether it existed or was written for this change
- BASELINE RUN commit (named baseline + pinned test) - distinct from the named baseline; label both:
- Candidate commit:
- Construction worktrees: implementer branch/path, test-writer branch/path

## Verification runs
(One entry per run, in order. Command, working directory, interpreter path, import-check output, exit status, output excerpt. Include the baseline failure and why it demonstrates the intended problem, the candidate pass, and each regression check. A pass is recorded as a test result with that test's limits.)

## Review
(Reviewer findings verbatim, each with its written disposition: repaired at commit X / held for stage Y with reason / not accepted with reason. Note whether an outside review was also obtained and its outcome.)

## Repairs and allowance
(Each repair round: what failed, its classification — implementation, test construction, setup, environment, meaning — what was done, by whom. State the agreed allowance and how much was used.)

## Pre-integration recheck
(Confirmation that candidate and target matched what was evaluated; what was refreshed if they had moved.)

## Decision and integration
(The developer's decision and the options presented. Integration commit. Post-integration tree state.)

## Limitations and open items
(What this change does not establish; anything held open, and where it goes next.)

## Status (generated, exactly one place)
- Status: OPEN   <- exactly this line, first in this section, `- Status: OPEN` while the change is open and `- Status: CLOSED` once /wrap has run; the installer's open-change check reads this line, so it is never omitted or reworded
(change id; parent change, if this continues a sequence; outcome class; the identities; what has and has not been run; the facts this was generated from. Everything elsewhere in this record is dated history.)

## Timeline
(session start and end; a timestamped line at every stage boundary, every question put to the developer, and every answer; and the three closeout numbers: wall-clock, machine-working, and time waiting for the developer)

## Procedure release
(the release identity this change runs on; chosen at session start, never upgraded mid-change)

## Server suite rosters
(baseline and final, compared as rosters; or NOT APPLICABLE per profile)

## Incident log
(empty, or each incident's five steps and authorized recovery)

## Standing permission grants
(any rules added to .claude/settings.local.json during this change, quoted, with the approval that authorized them)

## Material approvals and unexpected prompts
(scope changes, waivers, recovery authority, publication decisions - always. Routine permission prompts only when UNEXPECTED, or when the developer asked for full telemetry for a trial.)

## Approvals record
(a SEPARATE document OUTSIDE the frozen packet, referencing the packet's identity hash: checker sheet, Helper approval or developer waiver, developer authorization - appended after freezing, never inside)

> Two other files live beside the change records: `HANDOFF.json`, the state file (from `HANDOFF_TEMPLATE.json`), kept live through the change - rewritten at every stage boundary and state change - and finally at closeout; and `known_failure_variation.md` (from `KNOWN_FAILURE_VARIATION_TEMPLATE.md`), the tests whose failures vary from run to run.
