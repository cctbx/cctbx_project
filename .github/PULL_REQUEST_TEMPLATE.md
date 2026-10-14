<!--
Keep this short. A few lines is the norm in cctbx: say what changed and
why, then how it was tested. Describe the problem you solved, not just
what you changed. One problem per PR; split unrelated changes into
separate PRs. Delete any section that does not apply.
-->

## Problem
<!-- What broke, what was missing, or what was slow: a traceback, a user
     report, a failing test, a measurement. "Could theoretically fail" or
     "a review tool flagged this" is not a problem statement. Link the
     GitHub issue if there is one; "Fixes #123" closes it on merge. -->

## Changes
<!-- A sentence or a few bullets. Give the reasoning only where the diff
     does not show it. -->

## Tests
<!-- Tests added or changed, and which test suites you ran on which
     platform. For a bug fix, name the regression test and confirm it
     failed before the fix. If something was not run, say so. -->

## Coding agent
<!-- If an AI coding agent wrote any of this PR, say which one and where
     it ran; a change is weighed by what produced it. Write "none" if the
     PR was written by hand. -->

| Field | Value |
|-------|-------|
| Model and version | |
| Harness and version | |

## Checklist
- [ ] Searched open and closed PRs and issues for prior work (link them,
      or say "none found")
- [ ] New or changed tests are registered in the module's `run_tests.py`
- [ ] Docstrings and PHIL `.help` strings cover the new or changed behavior
- [ ] No new external dependency, or it was agreed with the developers
- [ ] Interface or breaking changes are called out under Changes
- [ ] A human has reviewed the complete diff
