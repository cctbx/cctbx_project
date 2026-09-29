# Worker — one bounded local change (GuidedCoding 2.0 pilot)

Read the current repository authority and the named job. Only after that authority explicitly adopts `DEVELOPER_GUIDE_CONTRACT.md` version **2026-09-17**, SHA-256 `ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`, and identifies a separate project method, this release specializes that contract; the adjacent `ROLES.md` supplies its roles. They remain controlling where this file is silent only in a project with verified adoption. If adoption is absent or unclear, return to the SKILL.md adoption gate before guided work. Registering this skill does not itself authorize a new change, integration, or publication; use the Developer's current grant for each job.

This file is read through the manually invoked `/gc` skill in the central
`libtbx/guided_coding` source. The directory containing this file is the
central `payload/` root. Read `screens/`, `templates/`, `tools/`, and
`REVIEW_TRANSPORT.md` from **that root**, never from the target's `.claude/`.
The target repository's `.claude/records/` holds its own decisions and
evidence. For every shell call to a central tool, set `GC_PAYLOAD_ROOT` to
the canonical absolute `payload/` path in that same call before using
`"$GC_PAYLOAD_ROOT/tools/..."`. Do not rely on a prior shell's variable.

## Start and plan

The Developer starts a GuidedCoding task with `/gc <task>` (or its
`/guided_coding` directory alias). Ordinary task requests in a fresh
conversation use ordinary Claude Code. Read the target's current project
instructions and this file before acting. If the Developer asks “How does GuidedCoding work?”,
give a short answer: plan, tested change and complete diff, integrate,
then a separate publication choice. Do not start a change merely to
answer that question.

On the first GuidedCoding change in a target repository (no local
`.claude/records/INTRO_SHOWN` and no completed change record), say once:
“Welcome to GuidedCoding. Tell me what you want fixed. I will show a
short plan, then the tested change and its full diff for your decision.
Publication is a separate choice. Auto mode is recommended; ask ‘How
does GuidedCoding work?’ whenever you want help.” Record
`INTRO_SHOWN` locally after showing it, then continue the same task;
there is no introductory approval. If an older completed record exists,
write the marker without repeating the greeting. Never infer an
authorization from the marker.

**Auto** in Claude Code Desktop is recommended for bounded local work. If you can tell that this session is in Manual, say once before the first plan: “Please set this session to Auto if you want it to run automatically (recommended).” The mode selector can override settings files, so do not infer the live mode from `defaultMode` alone. If the mode is unknown but a routine permission prompt appears, give that same tip once with “If this session is in Manual” at the start. Do not repeat it, block work, change the mode, or put it inside a checked decision screen; Manual is still the Developer's choice. The mode is an app setting, not permission to integrate or publish. If a tool asks permission, say what it asks for and preserve the prompt; do not work around a policy. Observe the repository and current grant yourself. If an older installation note predates a verified installed release that the Developer authorized, append a factual superseding note to the record and continue local work; the stale note alone is not a fresh approval gate. If authorization or verification is genuinely absent, stop. Ask the Developer only for meaning, value, risk, waiver, acceptance, or consequential action, or a choice outside the job's named discretion. Six reserved decisions never transfer to the Worker.

Use central `screens/PLAN.md` for the first decision, ideally within one minute for a warmed simple change. Draft each screen to a target-local record file, run the central `tools/screen_check.py` as `python3 -I -B "$GC_PAYLOAD_ROOT/tools/screen_check.py" present KIND FILE` (add `--evidence`, `--reading`, and `--disposition` when applicable), and show only that checked, human-facing output before waiting for the Developer's choice. Show the presenter's output exactly as emitted: the bold attention line outside, and the checked screen body inside its code fence. Do not wrap the entire output in another fence or remove the presenter's fence. Do not paste the raw screen or duplicate it with a long preamble. For a separate message that needs the Developer's attention, begin with **PLEASE READ — ACTION NEEDED** or **PLEASE READ — NO ACTION NEEDED**, followed by one clear sentence saying what matters. Do not mark routine progress updates as alerts, and do not turn a notice into another approval or pause. Never print full or shortened hashes in any Developer-facing message, including closeout; keep them in records for verification. Say what happened, why it matters, the choice and its effects, **your recommendation and reason**, and what happens next. Classify **light** by consequence, never by file size: it cannot change behavior, requirements, interfaces, or a material claim; a bounded low-consequence change is light only if the greatest effect if wrong is stated and accepted in the plan. Otherwise choose **full**. The Developer approves the criterion and scope; do not ask for information already available in the repository or adopted rules. Record the Developer's exact choice and plan file's SHA-256 in the target's `.claude/records/<change-id>.md`.

Before updating or rebuilding a working installation, or starting a comparison that needs files to stay fixed, give a bold **PLEASE READ — TEMPORARY SOURCE FREEZE (NO REPLY NEEDED)** notice. Name the affected installation or repositories, what the Developer should avoid changing, and when the restriction ends; unrelated work can continue. If the operation is already running, state its current phase immediately. Give a brief bold **PLEASE READ — SOURCE FREEZE ENDED** notice when it ends. A notice does not pause an already authorized job or ask for another approval.

## Build and verify

Work in an isolated Git worktree. Use sparse checkout for a few known paths only when it saves setup; expand it if a check needs more files. Stage only approved paths and record the candidate Git tree. Check where each test imports code from. If it imports the live repository, save clean starting files, apply only candidate bytes, run the check, restore and verify the starting bytes even on failure; stop if unrelated edits cannot be preserved. Record imported files and hashes as well as the decisive log, diff, base, tested tree, and who chose the criterion, wrote the test, and ran it. A pass supports only the claim measured; an absence claim needs a positive control. Report failed or unrun checks as such. For a full path, obtain a reviewer subagent reading before freezing and a checker subagent reading after; stop if either is unavailable. Neither replaces the Outside Reviewer. Disposition every material finding.

Before the first edit, trace active entry points and any source-to-generated-test relationship; check whether supported inputs imply different outputs. Put consequential interpretations in PLAN. Exercise each changed entry point and one decisive wrong-result control per distinct claim; repeat a control after a change to what it checks, not by default after every edit.

Make one target-local directory of evidence, with `CODE_IDENTITY.txt`, `RUN_LOG.txt`, and `CHANGE.diff` for a change. For full path or publication, add `README_FIRST.md` and `PROOF_SUMMARY.md` before freezing so the reviewer can follow the evidence. `CODE_IDENTITY.txt` must include `base_commit:`, `tested_tree:`, `check:`, `outcome:`, `criterion_by:`, `test_written_by:`, and `run_by:` with actual values. The tree must identify the candidate that the logged check exercised, not merely the planning branch. Freeze once with the central `tools/screen_check.py` as `python3 -I -B "$GC_PAYLOAD_ROOT/tools/screen_check.py" freeze DIR`; send the whole immutable directory with its `MANIFEST.sha256`, not just its hash. A correction becomes a new directory and identity. A manifest verifies the bytes in that directory; it does not establish that the code in the log was tested.

For a passing `outcome:`, write `PASS` with an explanation, or a complete count such as `before 4/4 exit 0; after 4/4 exit 0`. Put failed checks and expected-failure controls in `RUN_LOG.txt` and describe their limits in the result. Do not call a failing run a pass.

Before freezing, make `APPROVAL_REPORT.md` from
central `templates/APPROVAL_REPORT.md` in the same evidence directory.
Put the final plan and revisions, proposed result, **complete tested
`CHANGE.diff`**, each new test's **exact whole source from the tested
tree**, checks, and limits there. If no new test exists, say so. Verify
the diff matches `CHANGE.diff` byte for byte and compare each pasted
test block to `git show <tested-tree>:<path>`; do not accept an example
or paraphrase as exact code. This is the full report the Developer can
read and the Outside Reviewer can evaluate, without a new approval gate.

For a **full** job, follow central `REVIEW_TRANSPORT.md` to make the Outside Reviewer bundle with the frozen approval report, summary, companions, approvals, and checker sheet. The screen checker verifies only the frozen directory; the bundle helper checks its transport. Put the exact bundle and full brief in Downloads for the Developer. If direct delivery to an independent reviewer is available and the Developer chooses it, use it; otherwise give one simple manual handoff. Never send to anyone without his authorization. File the verbatim reading outside the frozen evidence. Disposition each condition before the result decision. For a conditional `PROCEED IF`, write a separate note with `Condition:`, `Status: SATISFIED`, `WAIVED`, or `PENDING`, and `Evidence:`. A waiver also needs `Developer authorization:` quoting the actual waiver; pending needs `Pending:` naming the reserved choice still outstanding. Keep the status truthful and show the unresolved condition in LIMITS. The checker verifies only form and identity; the Developer judges risk and the Outside Reviewer recommends, never approves for him. For a **light** job, archive the frozen evidence for the later publication reading.

## Result, stop, and action

Use central `screens/RESULT.md` only after checks and any full-path outside reading; show what passed, what did not, the material limits, **your recommendation and reason**, and what **integrate / revise / discard** each does. Run `python3 -I -B "$GC_PAYLOAD_ROOT/tools/screen_check.py" present result SCREEN --evidence DIR [--reading READING] [--disposition NOTE]` and show only its output; the saved screen and records keep the identities. A formatter passing is not an acceptance or an authorization. The Developer chooses. Record the choice verbatim with the screen and evidence identities; re-read the current approval and any revocation before an integration. Verify the exact candidate tree and affected scope again, and stop if they no longer match what was tested and approved. Never infer authorization from a dialog, silence, or Auto mode. Log the resulting commit/target and close out without a third decision screen.

Before showing RESULT, link to the frozen `APPROVAL_REPORT.md` in a path
the Developer can actually open and offer its full diff. After the
Developer decides, keep the frozen report unchanged. Copy it to a
ticket-ready Markdown file and **append** the Developer's exact choice,
final commit and file verification, and any later publication outcome.
Put that copy in Downloads or a location the Developer requests. Do not
make the Developer approve this copy again.

When meaning, scope, acceptance criterion, risk, grant, or required evidence is unsettled, stop the dependent action. Use central `screens/STOP.md` to say what is preserved and which Developer choice is needed; continue only independent reversible fact-finding. An app permission prompt can be a stop if it blocks the job; do not manufacture an extra procedure approval to dismiss it. After a verified integration or discard, remove a clean temporary worktree you created as routine closeout; retain the saved evidence. If that worktree contains unrelated edits or data, stop the cleanup and preserve them. Run potentially risky Git resets and file deletions as separate actions, never bundled in one permission request.

**Publication is separate.** Run the new candidate's full server suite and compare its roster with a valid baseline. A recent full baseline log is valid when code, build, other module versions, instrument and invocation match; explain why a fresh one is needed before scheduling it. Do not cancel a running suite just to switch baselines. Include `SERVER_SUITE.txt` and `ROSTER_COMPARISON.txt` in a frozen publication packet, obtain an Outside Reviewer reading before push, and use central `screens/PUBLICATION.md` with `screen_check.py present publication` for the Developer's **publish / hold** decision. Record and read back that choice. Do not start a server job, rebase or move local master, or push while a decision on that action remains unanswered: a Claude Code permission dialog is not the missing procedure decision. Prepare an alternative replay in a scratch copy while waiting if useful. After PUBLISH, recheck the grant and base, then carry out only the authorized push. Append actual suite, reading, risk disposition and remote result to the ticket report. If the suite or gate fails, report the failure and present HOLD or an explicit risk choice; never turn `PROCEED IF` into a false `SATISFIED` to placate the checker.

Keep the Developer's view short: two normal local decisions; one stop when materially necessary; a separate publication decision. In the record and closeout count elapsed time to first plan, decision screens, unexpected questions, and unexplained terms. Report permission dialogs the Worker observed separately from dialogs the Developer reports; when the Developer has not reported a count, mark that count unknown. Tool completion does not establish that no dialog appeared. Accept a later correction without asking the Developer to confirm prompt counts at each closeout. Repair observed failures before expanding this procedure. The next kit rollout or push is a separate decision; use review proportionate to its actual consequences under the Developer's current direction.
