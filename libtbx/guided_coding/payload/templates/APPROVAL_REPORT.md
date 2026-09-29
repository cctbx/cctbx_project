# Approval report — replace this heading with the task and date

**Status: proposed; no Developer acceptance or publication is implied.**
Prepare this file from the **tested tree**, before freezing evidence.
Include it as `APPROVAL_REPORT.md` in the evidence directory, and give
the Developer a usable path to the same file before INTEGRATE. It is the
substantive report the Outside Reviewer can evaluate and the Developer
can approve. The Outside Reviewer never approves on his behalf.

## Problem and intended result

Describe the observed bug, why it matters, and the agreed meaning.

## Decisions so far

Copy every PLAN and STOP screen from this run, including revisions that
changed scope or the criterion; quote the Developer's decisions. Do not include hashes in
the short human explanation; keep identities in the frozen evidence.

## What was built and checked

Copy the proposed RESULT content, including checks, failed controls,
unchanged existing tests, and material limits. Name the tested Git tree
in the audit record. Do not state that a review has passed before it has.

## Exact tested changes

List **every changed file**, then put the **complete `CHANGE.diff`** below
in a fenced `diff` block. Do not replace it with an excerpt, summary,
commit message, or a diff from a different tree. Compare the displayed
diff against the frozen `CHANGE.diff` before asking for approval.

## New tests — exact source in the tested tree

List each new test by path and test name. For each, copy the **whole test
function or added test file** from the tested tree inside a code fence,
including setup, assertions, and decorators. Compare this block with
`git show <tested-tree>:<path>`; the source must match. If no tests were
added, write `None added` and explain the check used instead. A diff
may also show the code, but it does not replace this readable list.

## Review, limits, and recommendation

After the outside reading, attach its verbatim file outside the frozen
packet and disposition each material finding. The proposed report above
remains unchanged. On the final decision screen, name the reviewer
finding and remaining risk; recommend INTEGRATE or DISCARD with a reason.

## Ticket copy after the Developer's decision

Do not change the approved report above. Make a copy and **append** the
Developer's exact decision, the final local commit and file check, any
PUBLISH/HOLD decision and remote outcome, plus unresolved limits. Put
the ticket copy in Downloads (or a location the Developer chooses).
Append the final RESULT and PUBLICATION screens if they were not present
in the proposed report, without altering the text approved above.
The Developer may attach it to the ticket; this copy is a record, not a
new approval request. On REVISE, rebuild and re-review the report for
the new tested tree instead of editing a frozen one.
