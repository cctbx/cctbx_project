# Review and reassessment

> **RECORD — not edited after the fact.** This is the review as written. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

Written after the improvement round, before deciding what comes next.

**State:** 141 tests, 51 negative controls all discriminating, full
validation 10 of 10 over 253 logs, house suite green on Tom's machine.
12 output channels. 657 distinct logs held; the extractor has been run over
all of them.

---

## 1. What the review of this round found

**One defect, and it is the class this project keeps meeting.** A `Sorry:`
block ran to the 40-line cap and stopped there **by luck** — the real end
was a box-drawing border one line later, in the agent's own framed output.
A cap that silently decides a boundary is silent truncation.

Fixed twice over: the boundary rule now also ends at a box border, a
`Skipping`, or a `[NOTICE]`; and **cap-truncation is flagged and printed**
rather than hidden. Across all 657 logs, 31 `Sorry` blocks, median 13 lines,
max 40, and **none now hits the cap**. Control C51 holds it.

**Clean on probing:** path elision does not touch `R-free = 0.23/0.25`,
`ratio a/b/c/d/e`, or short real paths; preamble regions are 36 lines
(cap 90); the ledger still balances; no line is in two states.

## 2. A number moved, and the documents have not caught up

Adding `terminal_failures` to the diagnostic set changes the headline
metric that several documents quote:

| | before | now |
|---|---|---|
| `corpus2/work` ok | 30% | **33%** |
| `corpus2/work` err | 9% | **17%** |
| the 291 fresh logs | — | **37%** |

The err figure nearly doubling is the point of the whole round: failed runs
now yield diagnostic structure because the program's own account of
stopping is finally read.

**`EXTRACTOR_ARCHITECTURE.md`, `MEASUREMENTS.md` and
`EXTRACTOR_REQUIREMENTS_v2.md` still quote the old numbers**, along with
"11 channels" and older test counts. That is document drift of exactly the
kind the project warned about, and it is now the largest open defect.

## 3. The holdout is spent, and I spent a little more of it

Re-measuring reach required reading `corpus2/holdout`. It came back 6%,
unchanged — but **that number is now in-sample**, because it was taken
after the change it measures.

The one-shot holdout has done its job and cannot do it again. **The 291
logs Tom sent are the replacement**, and they are better: never used to
build anything, 99% structural, 99% basic, **37% diagnostic**, 0
misidentifications. Going forward the holdout should be retired in the
documents, with its original numbers preserved as the record of what was
measured when it was fresh.

---

## 3b. A regression I introduced, found by Tom running it

Tom's verdict on the AutoSol report: *"pretty thin output but there are
lots of useful numbers in the log"*, with `Map-model CC: 0.47`, `R: 0.42
Rfree: 0.48` and the refinement block as examples.

**Both of his examples were captured all along. The report hid them.** The
previous round collapsed 342 labeled values into one `BACKGROUND` line to
fix "hard to read" — and the values are where the numbers are. Fixing
readability by hiding the content is not fixing readability.

Three changes:

- **Values are listed again**, grouped by label (342 → 261 distinct on that
  log), showing the last value and the line, **most recent first**.
  Chronology is not a judgement about importance, but a run's last numbers
  are the ones a reader came for, and ascending order buried them under 30
  lines of setup. The AutoSol report now opens on `CC: 0.553`,
  `Model-map CC (around atoms): 0.69`, `Built: 52`, `R: 0.42`.
- **`--label <name>`** pulls one label's whole series:
  `--label "Map-model CC"` gives 0.52 → 0.44 → 0.47 across the run, with
  line numbers. Exact match preferred — `--label R` first matched every
  label containing the letter r, 196 of them.
- **Compound pairs split.** `R:   0.42  Rfree:   0.48` was one value with
  the second number buried; 9% of key:value lines across 554 logs carry two
  or more pairs. Split **only when the pairs account for the whole line**,
  because a naive split turns `Space group: C 1 2 1 (No. 5)` into `C`.

**What is still not captured, and honestly:** `Model is in /path`,
`This model has selenomethionine as MSE`, `Adding waters in refinement` —
prose with no `label: value` shape. Reading them needs either a phrase list
or a file-mention channel; 324 of that log's 3,877 lines mention a file, so
the channel would be mostly noise. Left unread and visible in the unparsed
ledger rather than guessed at.

**145 tests, 53 controls, validation 10 of 10** (148,590 items verified).

## 3c. `--grep`, and a correction to something I told Tom

**I said the unread prose "sits in the unparsed ledger". It does not.**
`unparsed` holds only lines the frozen SCREEN proposed — **737 of that
AutoSol log's 3,877** — and `Model is in /path` matches no screen rule, so
it is in neither the channels nor the ledger. A `--grep` built over
`unparsed`, which is what I had described, would have missed exactly the
lines it exists to find. It searches the raw lines.

The same correction applies to what I told Tom about the consuming
project: `structure.unparsed` does **not** carry those lines either. A
program wanting them has to scan the text.

```
$ ... --grep "Model is in" AutoSol_run_1_1.log
LINES MATCHING 'Model is in' (5)
    2719  [not read] Model is in .../initial_build_cycle_1.pdb
    3705  [not read] Model is in .../full_model_1.pdb
    3778  [not read] Model is in .../full_model_combine.pdb
    3790  [not read] Model is in .../full_model_combine_ref_001.pdb
    3811  [not read] Final model is in  .../phase_and_build_2.pdb
```

Every hit is tagged with what read it, so `[not read]` and
`[labeled_values]` are distinguishable — the difference between "the tool
saw this and filed it" and "the tool never saw it". Searching
`Map-model CC` shows both, including a summary line at 2730 carrying the
same numbers in a form nothing reads.

**Two defects found reviewing it:**

- **Options were consumed by value, not position**, so
  `--grep Sorry Sorry.log` dropped the file and searched nothing.
- **`--summary` silently ignored `--grep` and `--label`.** It now refuses
  with exit 2; silently ignoring a flag is how a user comes to distrust
  the output.

And a guard: `--grep` without the log text says so rather than reporting
zero hits, because a silent false negative is the failure this project
keeps meeting.

**150 tests, 54 controls, validation 10 of 10.**

## 4. Reassessment — what the next step should be

Ranked by what it changes for the consuming project, not by effort.

### (a) Reconcile the documents. Small, and it is now the biggest defect.

Four documents quote superseded reach numbers, channel counts and test
counts. The project's own rule is one current-state document with numbers
attached to their corpus; it is currently violated in four places by my own
changes. **Do this first** — it is an hour, and every later claim rests on
it.

### (b) Take the diagnostic-reach finding to the parent project. Unchanged in importance.

30–37% on known programs, **6% on unseen ones**. An evidence layer built on
stage tables and cycles helps wizard runs and has little to say about the
~60-program long tail. This has been the standing recommendation since P9
and it is still the most consequential thing in the project. It belongs in
Phase 0a's plan *before* its arms are run.

The terminal-failure channel changes this a little for the better — failed
runs went 9% → 17% — but not the shape of the conclusion.

### (c) The raw run-directory shape is still under-measured.

12 of 677 logs. Enough to justify the `PHENIX <program>` signal, not enough
to measure the shape. **Not worth chasing further** unless those logs are
easy to collect: the signal works, and the alternative is inventing
coverage I do not have.

### (d) What I would NOT do next, and why

- **No new parser.** The five gaps Tom found are closed and nothing in the
  657 logs is now both common and unread. The largest unclaimed classes are
  xtriage ASCII tables and PDB `CRYST`/`SCALE` records — real shapes, and
  each would be a vocabulary or a format reader, which is the treadmill.
- **No re-opening the holdout.** It is spent; treat the 291 as the fresh
  set and say so.
- **No screen extension.** Declined once on measurement (+157% candidates
  for ~6 labels); nothing since has changed that arithmetic.

### (e) The question I would put back to Tom

The tool now reads what a run did and what it said about stopping. **The
open question is no longer about the extractor — it is whether the
consuming project can use 33% diagnostic reach.** If Phase 0a needs a
finding on every run, this tool will not supply one on two thirds of them,
and no amount of parsing will: most Phenix programs do not tabulate
anything. That is worth deciding before more effort goes into either side.
