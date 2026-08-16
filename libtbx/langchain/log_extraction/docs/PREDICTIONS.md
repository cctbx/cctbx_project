# PREDICTIONS

> **RECORD — not edited after the fact.** This is the predictions as registered, and how they scored. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

Pre-registered before the phase is implemented, per plan §4. A number
measured only after the fact cannot surprise me; that is the whole point.

**Rule:** a materially missed prediction forces a documented model change
or a new defect entry. Explanation alone does not close it.

**Status of earlier phases:** P0, P1 and P2 were run WITHOUT
pre-registration. That is recorded as a procedural gap in
`REVIEW_P0_P2.md` §4.1, not retrofitted here — writing predictions after
seeing the answers would be worse than admitting the omission.

---

## P3 — phases, control skips, item exclusions

Written after measuring the raw line shapes (that much is legitimate
"measure before asserting") but before writing any parser, and before
looking at what the parsers would produce.

### What is already measured, so not a prediction

On `corpus2/work n=253`: 1,109 lines match `^\s*Running\s+\S` in 65 logs;
381 lines match `^\s*Skipping\b` in 39 logs, splitting into 99
`as`/`because`, 49 dash-form, and **233 matching neither**.

### P3.1 — the two-way skip distinction holds on the complaint case

`find_reference_14` (wave-1) yields **exactly 18 exclusions in 2 reason
groups of 15 and 3, and 0 control skips**. No band: this is the case the
consuming project exists to answer, and an off-by-one is a defect.

### P3.2 — `Starting` produces no phases

**0** phases derive from a `Starting` line, on both corpora. Hard.
An early draft of the requirements had this backwards.

### P3.3 — the 233 unmatched skips all surface

Every one of the 233 `Skipping` lines matching neither form appears in
`unparsed`, none silently dropped. Predicted count **233 ± 0**; the
requirement (§6) is that a `Skipping` line matching neither form is
reported.

> **MISSED. Measured 234.** Band was zero, so this needs a resolution,
> not an explanation.
>
> The single extra line is
> `groel_dock_refine__dock_in_map_4.log:258`:
> `Skipping as CC_mask is < 1/2 min_cc` — a control skip **with no subject**.
> The prediction was made with a loose proxy (`\b(as|because)\b` anywhere
> in the line); the parser requires `<what> as <reason>` with both sides
> non-empty, so this one falls through to unparsed.
>
> **Resolution: the parser is right and the prediction was wrong.** The
> item name is optional for exclusions because a rejected candidate is
> still a rejection; but a control skip with no subject does not say which
> phase did not run, and inventing one — or emitting `what=None` — would
> be exactly the "invents a phase that never existed" failure §4.5 warns
> about. Reported as unparsed is the honest answer for a single line.
>
> **Recorded as a model change:** the two-way skip distinction has an
> asymmetry that was not written down before. An exclusion may be
> anonymous; a control skip may not. Now stated in the module and pinned
> by a test.
>
> The deeper lesson is the one this project keeps re-learning: the loose
> proxy I predicted with and the strict rule I implemented were two
> representations of the same idea, and nothing forced them to agree.

### P3.4 — diagnostic reach finally becomes non-zero

Currently 0% by construction. P3 adds control skips and exclusions to the
diagnostic set (phases are *not* diagnostic — announcing that something
began is not saying what happened).

Predicted: **30 ± 10 of the 253 working logs** gain diagnostic reach.
Reasoning: 39 logs carry a `Skipping` line at all, and some of those
carry only unmatched forms.

Predicted for the error subset: **2 ± 2 of 23**. Failed runs mostly stop
before they skip anything.

### P3.5 — the unclaimed pile falls

Predicted drop in `unparsed/unclaimed`: **1,000 ± 400** from the current
57,906. Reasoning: 1,109 phase lines plus 148 matched skips, minus those
already claimed as section titles (`Running refine_ca ...` occurs 209
times and is also a form-B heading, so the overlap could be large).

The width of this band is honest: I do not know the overlap.

### P3.6 — wave-1 conformance rises to 15

Currently 13 passed / 6 failed / 1 skipped against suite v2. P3 should
turn `the complaint case is exclusions not skips` and `the cryo-EM branch
is named` green. Predicted **15 ± 1 passed**. The remaining failures need
stage tables and cycles (P4) and identification (P7).

### P3.7 — performance stays inside budget

Three more line patterns. Predicted **≥ 5.5 MB/s** corpus-wide, from 6.3.
If it lands below 5.0 the budget is breached and the honest move is to
raise it with a stated reason, not to quietly re-measure.

### P3.8 — a phase may also be a section

`Running refine_ca to refine and make full model for ...` is both a phase
announcement and a form-B section title. Requirements §3 says nothing is
deduplicated across kinds, so both must be emitted. Predicted: **at least
one log where a single line carries both a Phase and a Section**.

This is the one design decision P3 makes that is not forced: SPECIFIC
parsers (phases, skips, exclusions) ignore other parsers' claims, while
the GENERIC channel (`labeled_values`) respects them. Rationale — a
heading that is also a phase is two true facts; a heading that is also
captured as a nameless "fact" is one fact counted twice.

---

## P4 — stage tables and cycles

Written after measuring the raw shapes and before writing any parser.

### Already measured, so not predictions

`corpus2/work n=253`: 72 stage-table headers in 47 logs, 1,179 legend
lines, 103 `Cycle N` lines in 43 logs — of which 75 carry the
`R/Rfree=a/b Built=n Placed=n` form and 28 do not (`Cycle 3 of morphing
chain trace`, `Cycle 2 try 1 of density modification`).

`wave1 agent n=20`: 7 headers in 6 logs, 6 cycle lines in 5 logs.

### P4.1 — the refine trajectory, exactly

`refine_5_beta_blip` yields **28 stages**, first `r_free` **0.5371**,
last **0.4935**, and **exactly three** stages worsening `r_free` by
≥ 0.001: `1_xyzrec`, `1_realsrl2`, `1_adp`. No band; the conformance
suite asserts these.

The trap I have already found by measuring: the table has **29**
consecutive rows and the 29th is named `end` with `r_free` 0.4942. `end`
is a summary line, not a stage. A parser that takes every row gets 29
stages and a final value of 0.4942 — which is precisely the number the
consuming project noted as *hiding* the finding ("the summary 0.4880 →
0.4942 hides it"). **Predicted: `end` is excluded from `stages`.**

> **HIT, all of P4.1-P4.5 exactly.** 28 stages, 0.5371 → 0.4935, and the
> three worsening stages `1_xyzrec`, `1_realsrl2`, `1_adp`.

### P4.2 — the sentinel

`autobuild_4_bromodomain` contains `Cycle 2 R/Rfree=999.90/999.90
Built=0 Placed=0`. Predicted: that cycle is flagged `sentinel=True`,
`Built` is 0, and **no r_work/r_free metric of 999.90 is emitted as a
number a consumer could quote**. 999.90 means "no usable result"; passing
it through as an R-factor is worse than omitting it.

### P4.3 — `R/Rfree=a/b` is one key and two values

Predicted: a cycle line of that form yields **two** metrics, `r_work` and
`r_free`, not one string.

### P4.4 — a bare cycle counter is not a record

The 28 non-metric cycle lines (`Cycle 3 of morphing chain trace`) carry
no metrics. Predicted: they produce **no Cycle item** and are reported as
unparsed. A counter is not a record.

### P4.5 — stage rows are anchored to their header

A global row pattern matches **40** lines in `refine_5_beta_blip` where
only 29 belong to the table. Predicted: rows are taken only as a
contiguous run following a header, so the count is 28 (after P4.1) and
not 40.

### P4.6 — diagnostic reach rises

Currently 33 of 230 ok logs (14%) and 2 of 23 err (9%). Stage tables
appear in 47 logs and cycles in 43, with overlap unknown.

Predicted ok: **75 ± 20 of 230**. Predicted err: **4 ± 3 of 23**.

> **HIT.** Measured 69 of 230 ok (30%) and 2 of 23 err (9%).
> 1,008 stages and 72 cycles in the ok set; 4 stages and 3 cycles in err.

### P4.7 — conformance reaches 18

Currently 15/4/1. P4 should turn `form counts` (needs `stage_table` = 6),
`refine stage trajectory` and `the failed build is visible` green.
Predicted **18 ± 1 passed**, with only identification outstanding.

> **MISSED. Measured 16 passed, 3 failed, 1 skipped.** Two causes, both
> mine, and neither a defect in the parser:
>
> 1. **`the failed build is visible` also asserts a DECISION**
>    (`rebuild_in_place`), which is P5. I predicted from the test's NAME
>    rather than reading its body. The sentinel cycle, the `Built=0` and
>    the model-building control skip all pass; the test fails on a line I
>    had not read.
>
> 2. **`form counts` is a naming mismatch, not a count.** It asserts
>    `forms` contains `stage_table == 6`; the contract names the field
>    `stages`, so `forms` reports `stages: 6`. **The count is right.**
>    `stage_table` is the prototype's label, frozen into the suite before
>    the contract named the field.
>
>    **CORRECTION, after the rename was applied:** the naming was not the
>    only cause. `t_form_counts` also asserts `decisions >= 5`, which is
>    P5. The suite still reports 16 passed after suite v3. I diagnosed a
>    failing test from its first failing assertion and stopped reading —
>    the same mistake as predicting from the test's title, one line
>    further down.
>
> **Resolution:** no model change — the parser is correct on both. One
> item for sign-off (suite v3: `stage_table` → `stages`, count unchanged
> at 6), and a procedural note that a prediction about a test must be
> made from the test body, not its title. That is the same
> representation-boundary failure as P3.3: the name and the assertion are
> two accounts of the same thing and nothing forced them to agree.

### P4.8 — performance

> **HIT.** 6.1 MB/s, above the 6.0 it started at — the header anchor means
> the row pattern runs only inside tables.

Predicted **≥ 5.0 MB/s**, from 6.0. This is the phase most likely to
breach the budget: it adds a header scan plus a row scan. If it lands
below 5.0 I will say so and propose a stated budget change rather than
re-measure until it passes.

---

## P5 — decisions

Written after measuring the raw shapes, before writing any parser.

### Already measured, so not predictions

`corpus2/work n=253`: **1,283 `Setting` lines in 133 logs**, of which only
**208 carry the documented `X=Y as/because Z` form**. The other 1,075
split into shapes that state a value but *no reason*:

```
Setting output.overwrite to True                      129
Setting macro_cycles to  3                             69
Setting default value of  True  for  thorough_denmod   41
Setting chain type to  PROTEIN                         41
Setting up prediction keywords...                      25   (not a decision)
```

`wave1 agent n=20`: 197 `Setting` lines, 27 in the full form, spread over
**5 logs** (predict_and_build 8, autosol_1 6, autosol_2 6, autobuild_4 4,
autobuild_3 3).

> **P5.1, P5.2, P5.3, P5.5, P5.6 all HIT.** 208 decisions and 1,075
> reasonless settings refused, both exact. 0 from `Setting up`. Setting,
> value and reason separate. Diagnostic reach 69/230 (band 65-95) —
> unchanged from P4, because decisions land in logs that already had a
> stage table or a skip. Performance 5.9 MB/s (predicted ≥ 5.5).

### P5.1 — only the reasoned form becomes a Decision

Requirements §4.4 defines a decision as the program changing its own
configuration **and stating the reason** — "the program announcing a
branch". A `Setting` line with no reason is a value without a branch.

Predicted: **208 decisions** from `corpus2/work`, **27** from wave-1, and
the remaining 1,075 reported as unparsed with a named refusal, **not**
fitted with a second pattern.

This follows the P3 precedent exactly: a control skip with no subject was
refused rather than invented, and the 234 unmatched skips were reported
rather than pattern-fitted. The 1,075 reason-less settings are recorded
as a candidate v2 channel with their measurement, for a decision made
deliberately rather than on first sighting.

### P5.2 — `Setting up X` is not a decision

`Setting up prediction keywords...` is a gerund, not an assignment.
Predicted: **0** decisions derive from a `Setting up` line.

### P5.3 — setting, value and reason are captured separately

Predicted: `Setting rebuild_in_place=False as maps_only=True` yields
`setting="rebuild_in_place"`, `value="False"`, `reason="maps_only=True"` —
three fields, not one string.

### P5.4 — conformance reaches 19

Currently 16/3/1. P5 should turn `form counts` green (it needs
`decisions >= 5`; measured 5 wave-1 logs carry a full-form decision —
**exactly on the line, with no margin**) and `the failed build is
visible` green (it needs `rebuild_in_place`, which is present).

Predicted **19 ± 0 passed**, leaving only identification. The zero band
is deliberate: both remaining assertions are checkable in advance and I
have now read both test bodies, which is what P4.7 cost me.

> **MISSED. Measured 18 passed, 1 failed, 1 skipped.**
>
> **Not a model error — an arithmetic error.** 16 passed with 3 failed;
> turning two green gives 18 passed and 1 failed, not 19. The suite has
> 20 tests and one is a SKIP, which is not a pass. I predicted the number
> I wanted rather than adding up the one I had.
>
> Both substantive halves were right: `form counts` and `the failed build
> is visible` both turned green, and identification is the only
> outstanding failure. **Resolution: no model change; the prediction was
> arithmetic and it was wrong.** Recorded because a zero-band miss must be
> resolved rather than explained away, and because "I read both test
> bodies this time" was true and still not sufficient.

### P5.5 — diagnostic reach

Decisions are a diagnostic kind. 133 logs carry a `Setting` line but only
the reasoned form counts.

Predicted ok: **80 ± 15 of 230** (from 69). Predicted err: **2 ± 2**.

### P5.6 — performance

Predicted **≥ 5.5 MB/s**, from 6.1. One more line pattern, cheap.

---

## P6 — completion records

Written after measuring the raw shapes, before writing any parser.

### Already measured, so not predictions

`corpus2/work/ok n=230`, after clustering adjacent markers into events:
70 logs with **no** completion event, 119 with exactly one, 41 with more
than one. `corpus2/work/err n=23`: 21 with none, 2 with four.
`wave1 agent n=20`: 7 none, 7 one, 6 more.

Marker vocabulary and counts across `corpus2/work`: `Finished` 691,
`Job complete` 188, `usr+sys time:` 188, `wall clock time:` 188,
`EXIT STATUS:` 22.

The `Finished` forms split: **548 timestamped** (`Finished: Wed Aug 12
…`), **113 named child** (`Finished with search...`, `Finished with
predict_chain`), 25 bare, 5 other.

> **P6.1, P6.2, P6.4, P6.5, P6.6, P6.7 all HIT.** One record from a
> three-line ending. **113 named-child records, exact.** No API exposes a
> success/failure verdict. Phaser's `EXIT STATUS` stays `unknown`.
> Diagnostic reach unchanged at 69/230 and 2/23. Performance 5.3 MB/s,
> exactly the predicted floor.

### P6.1 — a multi-line ending is one record, not three

The Program Template ends with `Job complete` / `usr+sys time:` /
`wall clock time:` on three consecutive lines. Predicted: that block
yields **one** CompletionRecord, not three. Adjacent markers within 6
lines cluster.

### P6.2 — `Finished with X` names its child

Predicted: **113 records with `applies_to = "X"`** — a named child, not
`top_level` and not `unknown`. This is the only positive evidence of
attribution the corpus offers, and it was found by reading the shapes
rather than assumed.

### P6.3 — `top_level` requires positive evidence

A record is `top_level` only if it is the **last** completion event in
the file and nothing but blank lines follows it. Otherwise `unknown`.

Predicted: **131 ± 10 of 230** ok logs carry a `top_level` record, and
**1 ± 1 of 23** err logs. Every other record is `unknown` or a named
child.

> **MISSED on first measurement (105), then a defect was found and it
> came in at 131 — exactly.**
>
> The shortfall was real: in **26 logs** a child's `Finished with
> process_predicted_model` sat within six lines of the program's own
> `Job complete` / `usr+sys` / `wall clock` ending. Every marker kind
> differed, so the clustering rule merged them — and **the run's own
> ending disappeared into a record attributed to the child**.
>
> `Finished with X` is self-contained and must not absorb what follows.
> Fixed, with a test and control C39. The prediction was right and the
> parser was wrong, which is the outcome the exercise is for.
>
> This is the second clustering defect in one phase: the first (two
> separate `Job complete` lines merging) was caught by the phase's own
> fixture. Proximity clustering is a rule with more edge cases than it
> looks like it has.

The asymmetry is the point: 21 of 23 failed runs have no completion
record at all, and 1 has one that is not top-level.

### P6.4 — an empty list never means failure

Predicted: **0** places in the code where an empty `completion_records`
is turned into a claim about success or failure. This is non-goal 2.8 and
the reason the field is a list of observations rather than a verdict.
Enforced by a test that asserts the API exposes no such interpretation.

### P6.5 — phaser's per-module EXIT STATUS stays unknown

One wave-1 log carries `EXIT STATUS: SUCCESS` many times over. Predicted:
each becomes its own record with `applies_to = unknown`, and **none** is
promoted to `top_level` except at most the last.

### P6.6 — completion records are not diagnostic reach

A record that something finished is basic information, not a statement of
what happened. Predicted: `basic` reach may rise; `diagnostic` reach
**unchanged at 69/230 and 2/23**.

### P6.7 — performance

Predicted **≥ 5.3 MB/s**, from 5.8. One more line pattern over every
line. The margin is thin; if it breaches 5.0 I will say so and propose a
stated budget change rather than re-measure until it passes.

---

## P7 — program identification

Written after measuring the candidate signals, before writing the parser.
This phase has a **hard gate (no log misidentified)** and a **target
(≥10 of 20 wave-1 logs named)** that pull against each other on purpose.

### Signals measured, and two rejected

**Signal A — the line-1 `Starting <program>` banner**, after
`strip_agent_prefix`, with two degenerate tokens blacklisted
(`Starting phenix` with no name, `Starting libtbx.start_process`, the
launcher):

| corpus | named | abstained | **wrong** |
|---|---|---|---|
| corpus2/work n=253 | 146 (58%) | 107 | **0** |
| wave1 agent n=20 | 7 | 13 | **0** |
| wave1 gui n=20 | 7 | 13 | **0** |

**Signal B — REJECTED. The structural hypothesis fails.** Requirements
§4.9 proposed, untested, that "the invoking program writes its own
parameter block first". Tested on the 107 banner-abstainers in
corpus2/work: **6 correct, 50 wrong, 51 no root found.** The roots are
not program names — `scaling` is xtriage's root (27 wrong on its own),
`data_manager` is the shared scope. Salvaging it would need a
root→program map from the PHIL masters, which is non-goal §2.3.

Recorded as a measured dead end, which is what the plan asked for: "worth
trying before any vocabulary scheme."

**Signal C — REJECTED ON PRINCIPLE, not on its numbers.** The
program-emitted `Found phil, …/<stem>.eff` record. On corpus2/work
abstainers: 13 correct, 1 wrong, 93 absent. It looks usable — but the
GUI names the `.eff` after the job, which is also what names the log
file. Using it is **the filename by proxy**, and §4.9's flat prohibition
exists because `ai_analysis` derived the program from the filename and
then printed a confident report about a program that never ran. On wave-1
it is 0 correct and 3 wrong, which is the same weakness showing.

Declining a signal worth 13 logs is the point of having the rule.

### P7.1 — the hard gate

Predicted **0 misidentified** on all three corpora. No band. A wrong
program name is worse than none.

### P7.2 — coverage

Predicted **146 ± 0 named** on corpus2/work (the banner is deterministic
and already measured), **7 of 20** on wave-1 agent-shape, 7 on GUI-shape.

### P7.3 — the wave-1 target will NOT be met, and I am not going to fit to it

The conformance suite requires **≥10 of 20** named. The banner reaches
**7**. Requirements §4.9 records that "discriminative section titles"
reach 9 on wave-1 — **derived from the same 20 logs they were tested on**,
which §4.9 itself flags as guaranteed by construction.

I predict 7, I predict the conformance test fails, and I am proposing a
restated target rather than a vocabulary fitted to the 20 logs it is
scored on. See the sign-off item in the handoff.

### P7.4 — signals are reported

Predicted: every result carries `signals_fired`, so a run that identifies
from one signal where it used to use two fails loudly rather than
narrowing quietly. Abstention is `candidates == []`, never a
low-confidence guess.

### P7.5 — identification still cannot gate extraction

Predicted: the P0 invariant holds unchanged — `scan(text)` and
`scan(text, program_name=…)` produce identical layer-A items on all 253
logs, and the negative control that forces abstention changes nothing.

### P7.6 — performance

Predicted **≥ 5.2 MB/s**, from 5.3. The banner reads six lines.

---

## P8 — measurements, attached only

Written after measuring, before writing the parser.

### Already measured, so not predictions

Plan §2.9 restricts this channel to numbers arriving **attached to a
structure already parsed**: stage rows, cycle lines, and the agent metrics
block marked `source="agent"`. Generic `Key: value` stays in
`labeled_values`, uninterpreted.

A flat stage/cycle channel would emit **9,433 measurements** on
`corpus2/work n=253`.

**The scoring set is `agent_metrics_labels.json`: 75 field-instances
across 14 wave-1 logs, minus the 2 known-defective `Residues Built`
entries (digit strings of an unrelated R-free).**

### P8.1 — the score must be taken on GUI-shape, and it is 24/75

Measured both ways, because this is where the project already made this
exact error once (A2: "67 of 75 appear verbatim", measured against text
that contained the block):

| | attached | labeled_values | neither |
|---|---|---|---|
| agent-shape | 12 | **61** | 2 |
| GUI-shape | 12 | 12 | **51** |

The agent-shape 73/75 is **circular** — 61 of those come from the agent's
own metrics block, i.e. from reading the answer sheet. **The honest number
is 24 of 75, on GUI-shape**, and that is what P8 reports.

Predicted: **24 ± 0** on GUI-shape, **73 ± 0** on agent-shape, and the
handoff quotes the GUI number.

> **MISSED. Measured 31 of 75 on GUI-shape**, with a zero band.
>
> The prediction was taken with a weaker matcher than the one the test
> uses. My probe compared `str(label)` against a set of stage/cycle values
> only; the test normalises numbers on BOTH sides (`2.50` and `2.5` are
> the same number, requirements trap 3) and also splits multi-value
> labeled values into tokens.
>
> **Third instance of the same failure** — P3.3 (a loose `as` proxy),
> P6.3 (a proxy that ignored the named-child rule), now this. Predicting
> with a quick probe and implementing with a careful matcher gives two
> accounts of one number, and nothing forces them to agree.
>
> Worse: the first version of the test itself used **string equality** and
> scored 15 — committing the exact error requirements trap 3 warns about,
> inside the test written to measure it. Corrected to numeric comparison.
>
> **Resolution: no model change; the prediction and the first test were
> both wrong, the parser was right.** 31 of 75 is the number of record.

### P8.2 — the agent block is reproduced and marked

Predicted: **75 of 75** label instances recovered from the agent block on
agent-shape, every one `source="agent"`, and **none** of them also
emitted as a `labeled_value`.

### P8.3 — the gap is the deliverable, not the coverage

The 51 uncovered values are **present in the log** and sit in shapes the
frozen screen does not take:

```
    All-atom Clashscore : 3.56      indent 4, space before the colon
  Clashscore            =   3.56    an `=` assignment
    overall d_min       = 1.744     an `=` assignment
```

Predicted screen-extension arithmetic, for a sign-off item, **not to be
applied in P8**:

| extension | new candidates | 75-label reach (GUI) |
|---|---|---|
| frozen | 41,023 | 24 |
| + any indent, space before colon | +8,055 | — |
| + `Key = value` | **+56,396** | — |
| both | +64,451 | **30** |

Predicted conclusion: **the extension is a bad trade and I will recommend
against it.** It more than doubles the candidate pile (+157%) to move the
label score from 24 to 30. Those 56,396 `=` lines are mostly PHIL
parameter echoes, which the consuming project already reads from its own
routing layer.

### P8.4 — no double emission

Predicted: **0** lines emitted as both a `measurement` and a
`labeled_value`, corpus-wide. Plan §2.9 requires it and nothing tests it
yet.

### P8.5 — reach and performance

Predicted: `diagnostic` reach **unchanged** (a measurement is basic
information). `basic` reach unchanged or +≤2. Performance **≥ 5.2 MB/s**
from 5.6 — the channel re-emits values already parsed, so the cost is
allocation, not scanning.

---

## P9 — the holdout, opened once

**WRITTEN BEFORE THE HOLDOUT IS TOUCHED.** The md5 of this file is
recorded in the build log immediately below, before any holdout log is
read. `corpus2/holdout` has never been scanned by this extractor.

The holdout is **87 logs across 44 whole programs the extractor has never
seen, plus 6 error logs**. Its originals still carry the agent wrapper, so
they must be stripped the same way `corpus2/work` was before scoring;
otherwise the comparison is between different shapes.

**Working-set baselines**, measured today at P8:

| | structural | basic | diagnostic | identified |
|---|---|---|---|---|
| ok n=230 | 97% | 96% | 30% | 54% |
| err n=23 | 96% | 39% | 9% | 96% |

### P9.1 — the hard gate

**0 logs misidentified.** No band, and this is the one that matters. The
44 held programs are, by construction, ones the banner has never seen —
if a degenerate banner form exists that I have not blacklisted, this is
where it appears, and it will produce a wrong name rather than an
abstention. That is the specific risk W4 has been recording since the
plan was written.

### P9.2 — identification coverage

Predicted **50 ± 15%** named on the holdout program set, against 54% on
the working set. The signal is positional, not a vocabulary, so it should
transfer; the band is wide because 44 unseen programs may write banners in
forms I have not met.

### P9.3 — basic-information reach

Predicted **≥ 81%** on the holdout program set — the working set is 96%
and the plan's failure condition is "more than 15 points below".

### P9.4 — diagnostic reach

Predicted **15 ± 15%** on the holdout program set, against 30%. Wide and
low deliberately: diagnostic structure comes from stage tables, cycles,
decisions and skips, and those live in the wizard programs
(refine, autobuild, autosol, predict_and_build), which are heavily
represented in the *working* set. A holdout drawn from programs with no
error logs is biased toward the long tail, where I expect little.

The plan's failure condition is >15 points below 30%, i.e. below 15%. I
am predicting the failure threshold itself, which is honest about how
little I know.

### P9.5 — the error holdout

6 logs. Predicted: **0–2 with a top-level completion record**, matching
the 1-of-23 on the working set. No meaningful reach prediction on n=6.

### P9.6 — the flat-log audit, pass 2

Pass 1 found all 7 zero-candidate logs in the working set were two-column
FSC data tables, not logs at all. Predicted for the holdout: **0–4
zero-candidate logs, and if any exist they are the same class** (data
files, not narrative logs with unreachable facts). A different class
appearing would be the more interesting result and would reopen §3.1.

---

## P9 — SCORED. The holdout was opened once, on the file hashed
## `ea3124b77d9a05a77b874e3a117f78ca`.

Stripped the same way `corpus2/work` was built: 93 logs (87 program
holdout, 6 error holdout), 0 empty.

| | predicted | measured | |
|---|---|---|---|
| P9.1 misidentified | **0**, no band | **0** | **PASS — the hard gate holds** |
| P9.2 identification | 50 ± 15% | **32%** | **MISSED, below band** |
| P9.3 basic reach | ≥ 81% | **93%** | PASS |
| P9.4 diagnostic reach | 15 ± 15% | **6%** | prediction HIT, **plan failure condition BREACHED** |
| P9.5 error holdout top-level | 0–2 of 6 | 1 | HIT |
| P9.6 flat-log audit | 0–4, same class | 2, same class | HIT |

### P9.1 — the hard gate held, and that is the result that matters

**0 of 93 misidentified**, on 44 programs the extractor had never seen.
W4 has been recording since the plan was written that an unblacklisted
degenerate banner form would appear here and produce a wrong name. One
did appear (P9.2) and it produced an **abstention**, not a wrong name.

### P9.2 — the miss, and the defect behind it

**18 of the 87 held-out logs write `**Starting phenix.molprobity **`** —
the banner wrapped in asterisks, a form absent from the entire working
corpus. The extractor abstained on all 18.

That is the whole purpose of holding programs back: a form I could not
have invented, and would not have found by staring at logs I had already
read.

**The 32% is the number of record.** The fix is one character class in
one regex and takes holdout coverage to 53% — but that number is
**in-sample**, measured after tuning to the set, and it is reported only
so the size of the miss is legible. A holdout is opened once.

Abstention breakdown of the 59 unnamed: 18 decorated banners (now fixed),
8 `Starting libtbx.start_process`, 5 `Starting phenix`, 28 with no banner
in the first six lines at all.

### P9.4 — the plan's failure condition is breached, and I am reporting it

Diagnostic reach on the holdout is **6%** against 30% on the working set.
Plan §7 fixed the failure threshold at "more than 15 points below", i.e.
below 15%. **6% fails it.**

I predicted 15 ± 15% and said in the same breath that I was predicting
the failure threshold itself, "which is honest about how little I know".
Being right about a failure does not make it not a failure.

**What it means, stated plainly:** diagnostic structure — stage tables,
cycles, decisions, control skips, exclusions — lives in the wizard
programs. The holdout was drawn from programs with **no error logs**,
which biases it hard toward the long tail: small single-purpose tools
that announce nothing. The extractor is not broken on unseen programs
(structural 98%, basic 93%, identification transferring at 0 wrong); it
is that **unseen programs mostly have no diagnostic structure to find.**

That is a fact about Phenix's long tail, not a fix to make. It bears
directly on the consuming project: an evidence layer built on stage
tables and cycles will help on wizard runs and will have little to say
about the other ~60 programs.

### P9.7 — what would count as failure

Stated before the draw, per plan §7:

- any misidentification → **hard failure**, and the banner rule is wrong
- basic reach below 81% → failure
- diagnostic reach below 15% → failure
- a missed prediction forces a documented model change or a defect entry
