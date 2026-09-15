# PROOF_SUMMARY template

Every proof packet's FIRST file, named `PROOF_SUMMARY.md`, written by the Worker. It is what a reviewer reads before anything else and it must be readable by someone who has never heard of this change. Plain words, no procedure vocabulary, no packet forensics. Six sections, in this order, nothing before them.

The checklist, the record and the evidence are what this summary is CHECKED AGAINST - they are not the front door. A reviewer who has to reconstruct the story from a completed checklist has been handed the wrong document first.

---

## 1. The problem

What was wrong, who hit it, and what it cost. If there was a user report, name it. If the damage was silent, say so, and give the measured size of it if there is a number. Two short paragraphs at most.

## 2. What was done, in order

The process actually followed: investigated, planned, what the Developer approved, who built what, how it was verified, who reviewed it, what was repaired and why. A numbered list of six to ten lines. This is where a reviewer learns that the test was written before the fix was seen, and by a separate subagent.

## 3. The fix

The change itself, in plain words, then the diff. Name every file touched. State explicitly what is UNCHANGED for runs that worked before.

## 4. The tests

What test was written or extended, where it lives, whether it is registered and runs automatically, and what each part of it checks. If anything was demonstrated by a probe or a throwaway harness rather than by the registered test, say which claims rest on it.

## 5. The results

Every check that was run, with its outcome:
- the decisive test on the old code and on the new code, with the actual values;
- the registered regression tests that reach the changed code, both sides, and what differed;
- the full suite comparison if one was run, by roster, with every category;
- the tree restored and verified, and by what means.

Numbers, not adjectives.

## 6. The verdict

For a change to data (expected outputs, tables, generated files): say what did NOT change (scientific or numeric values) separately from what did (layout, names, text), and where the developer reclassified rows the instrument flagged, keep the instrument's raw count in the sentence.

One of:
- **Complete.** The problem is fixed, the evidence is in this packet, and here is what remains open (list it, or say "nothing").
- **Complete with limits.** As above, naming precisely what was NOT established - platforms not tested, paths never exercised, questions still unanswered.
- **Not complete.** What is missing and what would settle it.

Then, always: the outcome class (published / saved locally / discarded / not reproduced / no candidate), the name of the packet's manifest file (the identity hash is the hash OF that manifest and cannot be written inside the packet it identifies; it lives in the approvals record outside), and the one-line statement of where the change now stands - on this machine only, or on the shared server as of which commit.

---

A summary that cannot be written because the evidence does not support it is a finding, not a writing problem. Say so instead of writing around it.
