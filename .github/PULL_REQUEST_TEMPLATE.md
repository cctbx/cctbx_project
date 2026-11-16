<!--
Keep this short; a few lines is the norm in cctbx. One problem per PR;
split unrelated changes into separate PRs. Write "Not applicable" under
a section that does not apply.
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
<!-- Fill this in if an AI coding agent wrote any of this PR; a change
     is weighed by what produced it. If several agents wrote or fixed
     code, list each, in the same order in the model and harness rows.
     If no agent wrote any of it, write "Not applicable" and delete the
     rest of this section.

     Do not open a PR on an agent's first draft: run the
     adversarial review rounds, with the severity scale, that
     CONTRIBUTING.md, "Commits and pull requests", describes. In the
     reviewer row, write "person" for a human reviewer; for several
     reviewers, list each and mark the one that ran the final round
     "(final)". -->

| Field | Value |
|-------|-------|
| Model and version | |
| Harness | |
| Rounds of adversarial review | |
| Reviewer model, version, and harness | |

<details>
<summary>Final review round</summary>

<!-- Tell the reviewer its output will be public, so it leaves out
     personal instructions, memory, and local paths. Paste its complete
     output between the ~~~~~~~~ lines, which GitHub shows as plain
     text. Edit it only to replace private details, such as tracker
     keys, local paths, credentials, email addresses, code from
     non-public modules, or anything from personal instructions or
     memory, with [redacted], before posting: GitHub keeps earlier
     versions of an edited description. If the output
     makes this description too long for GitHub, have the reviewer
     write a concise summary that keeps every finding and its level,
     and paste that. If the summary is still too long, post it in order
     across PR comments, each in its own ~~~~~~~~ block, and say so
     here. If the output has a line of 8 or more tildes, lengthen both
     fence lines. -->
~~~~~~~~
~~~~~~~~

</details>

<!-- For each final-round finding: its level if the reviewer gave none,
     then fixed (how) or left (why). -->

## Checklist
- [ ] Searched open and closed PRs and issues for prior work (link them,
      or say "none found")
- [ ] New or changed tests are registered in the module's `run_tests.py`
- [ ] Docstrings and PHIL `.help` strings cover the new or changed behavior
- [ ] No new external dependency, or it was agreed with the developers
- [ ] No imports of modules outside this repository; `chem_data` and
      `phenix_regression` are used only for data, with a skip when absent
- [ ] Interface or breaking changes are called out under Changes
- [ ] Agent-written code: the writing agents had `AGENTS.md` and
      `CONTRIBUTING.md` in context, and the reviewing agents had
      `AGENTS.md`
- [ ] A human has reviewed the complete diff
