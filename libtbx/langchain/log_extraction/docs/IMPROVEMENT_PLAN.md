# Improvement plan — from the AutoSol run

> **RECORD — not edited after the fact.** This is the improvement plan and what it found. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

Everything below is measured on the log Tom ran
(`AutoSol_run_1_1.log`, 3,877 lines) and against the full corpus. Ordered
by what it costs a reader, not by what is easy.

---

## 0a. Re-measured on 291 fresh logs — the evidence changed

Tom sent 305 logs, **291 of them new to me**. That is a free second test
set, and I ran everything against it before touching the plan.

**The extractor holds up on unseen data:** 0 crashes, **0
misidentifications**, 99% structural reach, 27% diagnostic reach, 19.5 MB
in 3.0 s. Identification names 24% — lower than the 58% on `corpus2`,
because most of these are agent-run logs whose first program line is not a
banner.

**But the failure evidence is completely different from what the project
has been assuming**, and that changes item 0 below:

| marker in PROGRAM output | new 291 (raw captures) | corpus2/work 253 (GUI) |
|---|---|---|
| `Traceback` | **25 (9%)** | 2 (1%) |
| `Sorry:` | **25 (9%)** | 4 (2%) |
| `Please try again.` | **19 (7%)** | 0 |
| `STOPWIZARD` | 2 | 1 |
| any terminal marker | **27 (9%)** | 5 (2%) |

The project has repeated since the corpus was built that *"27 of 30 failed
runs contain no error keyword at all in what the program itself wrote"*.
**That is true of the GUI-shape corpus and false of these raw captures**,
because in the GUI corpus the failure went to stderr and was stripped into
the agent wrapper. Here it is in the log, in full:

```
Traceback (most recent call last):
  ...
    raise Sorry("\n".join(data.err))
Sorry: If previously used R-free flags are available run this command again
with the name of the file containing the original flags as an
additional input. ...
xray_data.r_free_flags.generate=True
...
Please try again.
```

A traceback, the program's own diagnosis, **and the remedy** — all
quotable, all with line numbers. This is the single most valuable thing in
the new material and the tool currently says nothing about it.

**So item 0 is no longer a thin rule built on 3 STOPWIZARD logs.** It is a
terminal-failure block seen in 27 of 291 unseen logs, and one of those 27
(`predict_and_build_2.log`) currently reports a `COMPLETION` record and no
mention of the traceback at all.

---

## ALL FIVE ITEMS ARE DONE

| | item | result |
|---|---|---|
| 0 | terminal-failure channel | **40 of 657 logs** now report the program's own account of stopping; **12 of them also carry a completion record**, and both are shown |
| 1 | `PHENIX <program>` signal | AutoSol identified from line 3; identification 24% → 44% across everything held |
| 2 | preamble as a region | `Texas Engineering Experiment Station` gone; 59 → 57 sections on the AutoSol log |
| 3 | report ordering, path elision | how-it-ended first; 57 sections + 342 labeled values collapsed to one `BACKGROUND` line, `--all` restores them |
| 4 | unclaimed shapes | the six masked shapes shown with counts and a line number each |

**139 tests, 50 controls all moving, full validation 10 of 10** (148,358
items verified, 151 named / 0 wrong).

**A process defect worth recording.** The script that was supposed to add
the nine new tests silently did nothing — its `replace()` target did not
match and it asserted nothing, so it reported success. I then reported "all
131 tests passed" as though the new tests had run. **Two negative controls
caught it**: C48 and C50 disabled features and changed no test result.
Every edit script now asserts that its replacement applied.

**A judgement call, written down rather than assumed.** The non-goal 2.8
guard banned any attribute containing "failure" and fired on the new
`terminal_failures` channel. Loosening a non-goal guard is how a non-goal
dies, so the distinction is now explicit: 2.8 forbids a **verdict** about
the run, not a **channel** of things the program said. The guard bans
verdict-shaped names (`success`, `failed`, `outcome`, `status`, `verdict`,
…) and additionally asserts every observation channel is a list — a verdict
would be a bool or a string.

---

## 0. The finding Tom did not raise, and it is worse than the ones he did

**The run was ABORTED and the report does not say so.**

```
line 3877:  STOPWIZARD IN AutoSol_run_1_/STOPWIZARD
```

That is the last line of the file. What the report prints instead:

```
COMPLETION (1)
    1994  Finished with search...    applies to: search...
```

— a *child* completion at line 1994, with **1,883 more lines after it**. A
reader could reasonably conclude the run finished. It did not: it was
killed.

This is not a new discovery. The project record has carried it since the
corpus was built: *"`autosol_1_gene-5-mad` sits in the good set and ends in
`STOPWIZARD`"*, filed under "the err label is a lower bound". It was
recorded as a corpus-labelling caveat and never turned into a channel. Tom
then handed me that exact program, and the tool stayed silent about the one
thing that mattered.

**Fix: a `terminal_failure` channel** (renamed from `abnormal_end` on the
new evidence). It captures, with line numbers:

- a `Traceback` block, if present;
- the `Sorry:` diagnosis and the multi-line remedy that follows it;
- explicit abort markers (`STOPWIZARD`, `EndOfResolve`);
- and nothing else — an unrecognised ending still falls through to the
  existing "no completion record" silence rather than to a guess.

It reports what the program said. It does **not** conclude the run failed
(non-goal 2.8), and where a log carries both a completion record and a
traceback, it shows both and lets the reader see the conflict.

**Priority 1.** It is the difference between a report that is unhelpful and
one that is misleading.

---

## 1. DONE — it now says `phenix.autosol`, from line 3

```
program:  phenix.autosol   (self_identification, line 3)
```

Tom sent the `grep PHENIX phenix/programs/*.py` output, which settled this
from the **source** rather than from guessing at logs: ~60 programs print

```python
print("\n" + 60*"*" + "\n" + 10*" " + "PHENIX AutoBuild" + ...)
```

and `ai_agent`/`ai_analysis` substitute the name at runtime. The name is
whatever precedes the weekday, which is what makes multi-word names work
(`PHENIX Plan SAD experiment`), along with CamelCase (`AutoSol`) and dotted
(`phenix.ai_agent`) forms.

**Measured on all 657 distinct logs I hold:** the header is present in 36
and the name agrees with the filename in 35.

**The one disagreement is the signal being right and the filename being
wrong.** `ligand_identification_68.log` is a GUI job that runs LigandFit,
and says so three times — `Running LigandFit process 4`, `PHENIX
ligandfit`, `Starting LigandFit with the command:`. That is the **first
demonstrated defect in the filename-derived truth labels**, and it is
exactly why §4.9 forbids the filename as evidence. Recorded as a known
truth-label exception rather than used to weaken the gate.

Identification across everything I hold rises **24% → 44%**; on
`corpus2/work` **147 → 151 named, still 0 wrong**; full validation still
passes 10 of 10. A guard test pins the weekday anchor, because without it
`PHENIX components are copyrighted by:` — in the preamble of every Phenix
log — would become a program name.

## 1b. What remains on the shape itself

```
      3  '               PHENIX autosol  Sun Jul 19 10:20:17 2026'
```

The extractor looks for `Starting <program>` and this shape writes
`PHENIX <program>` instead. Measured across every corpus, in the first 12
lines: **10 correct, 0 wrong, and it is present in only 10 of 386 logs.**

That last number is the real finding. **This is a third log shape the
corpus barely contains** — the raw wizard run-directory log, as opposed to
the GUI shape (253 logs) and the agent-wrapped shape (20). Tom's actual
use is the shape we have 2.6% coverage of.

The signal is in. What is *not* fixed is the corpus gap: the raw wizard
run-directory shape is **12 of 677 logs**. That is enough to justify the
signal and not enough to measure the shape, so it is marked
under-measured in the architecture.

**Note what was NOT used.** The grep lists ~60 program names. Shipping that
list would be a vocabulary — the thing this tool exists to avoid. It was
used only to establish the SHAPE and to check the regex offline; no program
name appears in the module.

---

## 2. "Texas Engineering Experiment Station" — the Phenix preamble is boilerplate

Every Phenix log opens with ~50 lines of developer list, copyright,
third-party components and citation. It is identical in every log and it is
not about the run. It currently yields:

```
      16  http://www.phenix-online.org/
      26  Texas Engineering Experiment Station
```

The second is a copyright continuation line that happens to sit above a
rule, so the form-B section rule fires. Both are correct by the rules and
useless to a reader.

**Fix: recognise the preamble as a region**, the same machinery as the
agent wrapper — bounded, marked, and excluded from `sections` and
`labeled_values` — with its own `source` value so nothing is silently
dropped and the extent is auditable.

**Priority 3.** Small, and it is the first thing on the screen.

---

## 3. The report is hard to read because it is ordered by channel, not by interest

For this log the report prints, in order: 12 phases, 6 decisions, 3 skips,
18 exclusions, 1 cycle, 1 completion, **59 sections**, **342 labeled
values**. The three things a reader wants — *the run was killed*, *it
skipped model-building*, *cycle 2 produced nothing* — are 3rd, 4th and 6th,
and the two noisiest channels are the largest.

Also measured on this log: **long paths dominate.** `Working directory`
appears 4 times, phases are truncated because they carry 106-character
paths, and the most common labeled value is `NOTE` (10).

**Fix:**

1. **A findings header first**: identification, whether it ended normally,
   and a one-line count of what follows.
2. **Order by diagnostic value**: abnormal end → decisions → skips →
   exclusions → cycles → stages/trajectory → phases → completion. Sections
   and labeled values last, collapsed to a count unless `--all`.
3. **Elide long paths** to their basename plus the line number. A
   106-character path is not a finding; the line it is on is.
4. `--all` restores everything, because the truncation must never be the
   only version available.

**Priority 4.** No parser changes, only presentation — and presentation is
what Tom actually hit.

---

## 4. "NOT UNDERSTOOD" is opaque because it is counts without examples

```
NOT UNDERSTOOD  unclaimed=175  generic_only=499  rule_excluded=61
```

Three numbers and no way to judge them. Measured, the 175 unclaimed lines
are dominated by:

```
  30  ----------------------------------------------------------
  29  | #                          | #                      | #
  12  SCALE#      #  #  #        #
  11  | # - #  |  #     |  #     |  #     |  #     |  #     |  #
   4  unit_cell = # # # # # #
   4  CRYST#   #   #   #  # #  # C # # #
```

which is: **xtriage's ASCII tables, PDB CRYST/SCALE records, and bare
rules.** That is a completely different impression from "175 things we
could not read", and it is reassuring rather than alarming.

**Fix: show the shapes, not just the totals** — the commonest unclaimed
patterns with digits masked and a count, and one example line number each,
capped at ~8 rows. The full inventory stays available via
`run_measurements.py`.

**Priority 5.** Cheap, and it turns a number nobody can act on into a list
anybody can scan.

---

## 5. What I am deliberately not proposing

- **No vocabulary.** Not for xtriage's tables, not for CRYST/SCALE. They
  become recognised *shapes* (a pipe-delimited table row, a PDB record) or
  they stay unclaimed and visible.
- **No inference of failure.** `abnormal_end` reports an abort marker the
  program itself wrote, with its line. It does not conclude the run failed
  — non-goal 2.8 stands, and 27 of 30 failed runs still say nothing at all.
- **No re-tuning against this one log.** Every rule above is measured
  across the corpus first, and the raw-run-directory shape gets its own
  logs before its numbers go in any document.

---

## 6. Order of work, and what each is worth

| | change | evidence | risk |
|---|---|---|---|
| 1 | `terminal_failure` channel | **27 of 291 unseen logs**; Sorry+remedy quotable | markers are explicit and measured |
| 2 | `PHENIX <program>` signal | 10 present, 0 wrong | low; positional, not a word list |
| 2b | **get raw run-directory logs** | 2.6% of the corpus | none — it is measurement |
| 3 | preamble as a region | ~50 lines in every Phenix log | low; bounded and marked |
| 4 | report ordering and path elision | this log: 59 sections, 342 values | none; presentation only |
| 5 | unclaimed shapes with examples | 175 lines, 6 shapes | none; presentation only |

**Item 2b is only half answered.** The 305 logs contain **2** raw
wizard run-directory logs (`AutoSol_run_1_1.log`, `LigandFit_run_1_1.log`)
— the shape with the `PHENIX <program>` preamble. Across everything I now
hold, that shape is **12 of 677 logs**. The other 289 new ones are the
agent-wrapped and GUI shapes I already had.

That is enough to justify the `PHENIX <program>` signal (12 present, 0
wrong) and not enough to measure the shape. If raw run-directory logs are
easy to collect, another 20 would let me treat it as a corpus rather than
as two examples. If they are not easy, say so — I will implement the
signal, mark the shape under-measured in the architecture, and stop asking.
