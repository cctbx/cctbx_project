# HANDOFF — Phenix log structure extractor

**Single current-state document.** If it disagrees with anything else it
wins — except numbers, where `MEASUREMENTS.md` wins, because that file is
regenerated from the corpus.

**As of:** close of the build session. Extractor sha256 `aa2487b830737ddb`,
validated against corpus id `corpus2 + wave1 9d8689b3861aeeb2c7b32585aad350c2`.

**Read this first, then `EXTRACTOR_ARCHITECTURE.md`.** Numbers here are
reconciled against `MEASUREMENTS.md`, which is regenerated from the corpus;
if the two ever disagree, MEASUREMENTS.md is right.

---

## 1. Where things stand

Built, installed in `libtbx/langchain/log_extraction/`, registered in the
house test runner, and run by Tom on macOS under `phenix.python`.

| | |
|---|---|
| suite (tiers 1+2, ships in cctbx) | **151 tests, 0 failed**, ~1 s |
| negative controls | **54, all discriminating** |
| wave-1 conformance (suite v4) | **20 passed, 0 failed, 0 skipped** |
| full validation (tier 3, external corpus) | **10 of 10**, 148,590 items |
| house `run_all_tests.py` | 168 suites, 0 failed |

**Channels: 12.** sections · phases · stages · cycles · decisions ·
control_skips · exclusions · measurements · labeled_values ·
completion_records · **terminal_failures** · unparsed.

**Reach** (`MEASUREMENTS.md` is authoritative and regenerated):

| corpus | structural | basic | diagnostic | identified |
|---|---|---|---|---|
| corpus2/work ok n=230 | 97% | 96% | **33%** | 56% |
| corpus2/work err n=23 | 96% | 39% | **17%** | 96% |
| fresh raw captures n=288 | 100% | 100% | **38%** | 24% |
| holdout n=93 (spent) | 98% | 93% | **6%** | 32% |

## 2. What it does

Text in, structure out, every item carrying the line it came from. It reads
what the program says about itself and does not interpret, rank, diagnose or
call a model. Three layers kept apart by the call order: generic observation,
then identification (which may abstain), and no program-specific
interpretation at all.

Run it on a log:

    phenix.python .../log_structure_extractor.py my_run.log
    ... --label "Map-model CC" my_run.log     one label's series
    ... --grep "Model is in" my_run.log       any line, read or not
    ... --all my_run.log                      everything
    ... --summary *.log                       one line each

## 2b. The three commands

    phenix.python langchain/tests/run_all_tests.py          # every run: 151 tests, ~1 s
    phenix.python .../log_extraction/validate.py            # before shipping a change
    phenix.python .../log_extraction/run_controls.py        # occasionally: 54 controls, ~8 min

`validate.py` needs the 21 MB corpus, which is **not** in cctbx — unpack
`extraction_test_data.tar.gz` into `~/Downloads/` or pass `--corpus`. It
**fails** rather than skips if the corpus is absent, and writes a
validation record beside the corpus naming the extractor revision it
passed.

## 3. What is outstanding, none of it in this package

**(a) DONE — written up as `NOTE_FOR_PHASE0A.md`, in the plan's own terms.**
Measured per program family across 253 logs and 63 programs: the plan's 9
families yield event findings on **46%**, everything else on **13%**, and
**48 of 63 programs (62% of the corpus) yield none at all**. Three of the
four commonest programs by log count — xtriage, phaser, real_space_refine —
are 0% on the event axis and **85–100% on measurement**. The note's ask is
not more extraction: it is that the plan state the event/measurement split
as the expected shape, and take the Phase 1 scoping decision against a
corpus containing the long tail. Tom's call from here.

**(b) The raw run-directory shape is under-measured** — 12 of 677 logs.
Enough to justify the `PHENIX <program>` identification signal, not enough
to measure the shape. **Not worth chasing** unless such logs are easy to
collect; the signal works and the alternative is inventing coverage.

**(c) Never executed anywhere but Tom's Mac and this container.** One
Phenix build, one user, two machines, 657 logs. That is the assumption most
likely to be wrong and the cheapest to test.

## 4. What NOT to do next, and why

- **No new parser.** Nothing in 657 logs is both common and unread. The
  largest unclaimed classes are xtriage ASCII tables and PDB `CRYST`/`SCALE`
  records; each would be a vocabulary or a format reader, which is the
  treadmill this tool exists to replace.
- **No phrase list for prose** (`Model is in /path`, `Adding waters in
  refinement`). Measured: gerund-initial lines are 6.1% of all lines across
  166 distinct verbs; file-mentioning lines are 3.3%. Both too broad to be
  findings. `--grep` answers the question without a vocabulary.
- **No re-opening the holdout.** It is spent; the 288 fresh raw captures are
  the replacement set.
- **No screen extension.** Declined on measurement: +157% candidates to move
  the label score by ~6.

## 5. Documents, and which kind each is

**Current-state** — must be true now: `MEASUREMENTS.md` (authoritative for
every number, regenerated), `EXTRACTOR_ARCHITECTURE.md` (start here),
`EXTRACTOR_REQUIREMENTS_v2.md`, `DESIGN.md`, `INSTALL.md`,
`UNPARSED_INVENTORY.md`, this file.

**Records** — true as of their date, never edited: `PREDICTIONS.md`,
`REVIEW_P0_P2.md`, `FLAT_LOG_AUDIT*.md`, `PROCEDURE_REVIEW_REQUEST.md`,
`TEST_TIER_PLAN_v2.md`, `IMPROVEMENT_PLAN.md`, `REVIEW_AND_NEXT.md`. A
pre-registered prediction rewritten after the fact is not a prediction.

## 6. Standing limitations

- **~2.5% of any corpus is not a log** — two-column FSC data tables with a
  `.log` extension. Every reach metric carries that floor.
- **It cannot tell you a run failed.** It reports completion records and the
  program's own account of stopping; 27 of 30 GUI-shape failures said
  nothing at all, though 9% of raw captures do.
- **Identification abstains on ~41%.** Zero wrong across 657 logs; the one
  filename disagreement is the filename being wrong.
- **One Phenix build, one user, two machines.**

---

## 7. The full defect record

Everything below is the running record of the build — every defect with how
it was found. It is a **record**: section numbers below restart from the
cycle they were written in, and nothing here is edited to match later
figures.

---

## 5. Defects found and fixed this cycle

**D1 — the mid-file wrapper rule could never fire.** The first draft
ended a footer block at the next wrapper key and judged mid-file-ness by
how much text *followed* the block. With no following key the block
always ran to EOF, so the "remaining" count was always 0 and every block
was called `agent`. *Caught by:* its own fixture test.

**D2 — the second draft then misjudged short files.** Keying on start
position alone called a 5-line fixture's genuine footer ambiguous. The
rule now uses two signals: **position** (a real footer sits in the back
half — `corpus2` error originals put `ERROR TEXT:` at line 8041 of 8225)
and **repetition** (a real wrapper writes each key at most once, while a
log that *quotes* wrapper blocks repeats them — `p9-sad__ai_agent_1`
carries 8 `FINAL QUALITY METRICS` blocks). Either signal failing means
`unknown`. *Caught by:* the same fixture, after D1's fix.

**D3 — a corpus invariant was wrong, not the corpus.** "No wrapper
regions remain in the working set" failed on two logs. Those two are the
agent's *own* logs, which quote their children's wrapper blocks mid-file
— content, not residue. The invariant that actually catches leakage is
about the *ends* of the file: no agent header, and no `agent`-sourced
footer. *Revision recorded in the test docstring rather than made
quietly.*

**D4 — a feature was implemented and never tested.** Control C4 forced
every item to `source=program` and **the whole suite still passed**.
`find_regions` was covered; the propagation of a region's source onto
the items inside it was not. Two tests added (one fixture, one corpus).
This is the clearest evidence this cycle that the controls earn their
place — no amount of reading the tests would have found it.

**D5 — `bool` is an `int`.** `Phase("x", True)` silently became line 1.
Line validation now rejects bools.

**D6 — the module docstring's own example was Python-3 only**
(`open(path, errors="replace")`), in a module required to run under
`libtbx.python`. Replaced with the bytes-then-decode form.

**D7 — `Section.start` vs `Section.line`.** Requirements §3 names a
section's first line `start`; every other item calls it `line`. Both are
now provided, so the frozen conformance suite and the uniform item API
agree without either moving.

## 5b. Defects found in the P0 REVIEW pass (after the suite was green)

The suite was 32/32 green before this pass. These five are what a
deliberate review found on top of it, which is the argument for keeping
review as a separate step rather than stopping at a green suite.

**D8 — the agent header swallowed the program's own banner. The worst
defect of the cycle.** The agent does not add a line at the top; it
PREFIXES the program's first line — `LOG TEXT: Starting AutoBuild` is
nine characters of agent in front of program content. The first version
marked that whole line `source="agent"`. MEASURED: 15 of the 20
agent-shape logs carry the prefix and **in 9 of them the payload is a
`Starting <program>` banner** — the single best identification signal
there is, buried inside a region labelled "written by the agent". P7
would have had to either ignore the marking or re-derive the banner, and
one of those is a silent misread. The line is now its own region with
`source="program"`, and `strip_agent_prefix()` is exported for parsers.
Controls C10 and C11 added to hold it. Note the wave-1 GUI corpus still
carries this prefix on line 1 of 15 of its 20 logs, so tolerating it is
not optional.

**D9 — a test assertion that could not fail.** `test_screen_recognises_
its_four_shapes` contained `assert A == "verb" or B == "numeric_row"`
where B was separately asserted on the next line. Replaced with exact
single expectations.

**D10 — two invented fixtures asserted the opposite of what they
claimed.** A hand-written pipe-table row was given three numeric cells,
not four, so it screened as `None` while the test claimed it screened as
`numeric_row`. A second attempt, this time a real corpus line, had the
same problem. Both would have passed a reading and failed the data. The
fixture is now a line from `corpus2/work` verified against the rule.

**D11 — the numeric-row rule deviates from the plan's wording, and the
CODE is right.** The plan says "a row of ≥4 whitespace-separated numeric
columns", which reads as a majority rule. MEASURED on `corpus2/work
n=253`: ≥4 numeric tokens gives 46,643 candidates; adding a majority
requirement gives 35,280; majority ignoring separator cells gives
39,181. The 11,363 lines a majority rule drops include genuine table
rows — `| 2.85 - 2.75 | 92.0 | 82.6 | 72.1 | 57.7 |` (pipes outnumber
numbers) and phaser solution records `SOLU TRIAL ENSEMBLE ... EULER
11.606 60.800 213.257 RF -14.7 RFZ 2.50`. Over-selecting prose that
carries four numbers is the cheaper error for a screen whose job is to
over-select. **The plan's wording is what gets corrected, not the code.**

**D12 — control C8's patch target had gone stale.** It silently reported
"patch target not found" after D5 changed the line it patched. The
control runner already treats a skipped control as a failure, which is
why it surfaced; worth noting that the controls need the same
maintenance as the tests.

## 5c. The flat-log audit (plan §3.1, pass 1)

Run before any parser exists, sample fixed before reading: all 7 logs in
`corpus2/work` with zero screened candidates. See `FLAT_LOG_AUDIT.md`.

**All seven are the same thing and none is a log** — two-column FSC data
tables (`1/resolution  CC` then floats to EOF) written to a `.log`
extension. They match no screened shape for a correct reason: the
numeric-row rule needs ≥4 numeric tokens and these rows have 2.

So the review's hypothesised failure mode — a headerless *narrative* log
carrying facts that match nothing — **does not occur in this corpus**.
What the audit found instead is a corpus-composition fact: 2.8% of
`corpus2/work` is data files, and they must be reported as such rather
than counted as extraction failures, or every reach metric carries a
floor it can never close.

## 6. Defects found in the P1 cycle

**D17 — a pure rule was read as a title made of dashes.** The inline
pattern `[-=*#_]{3,}(\S.*?)[-=*#_]{3,}` matches `------------`, taking
the middle dashes for a title. The rule test now runs *first*. *Caught
by:* two fixture tests failing at once.

**D18 — `*` banners were treated as underlines.** MEASURED: of 8,644
form-B candidates, the rule beneath is `-` 5,967 times, `=` 1,818, `*`
824, `#` 35. The `*` rules wrap a block above and below — they are
delimiters — and the titles they produced were the junk
(`hydrogens.refine=riding`, `Unset PHENIX_OVERWRITE_ALL ...`). Form B
now accepts only `-` and `=`. *Rejected first:* a blank-line-above rule,
which would have refused 4,603 candidates including `OUTPUT FILES`
(240), `TRANSLATIONAL NCS` (192) and `WILSON DISTRIBUTION` (126). It
costs more than it saves — measured, not assumed.

**D19 — a blank line above a rule was emitted as a REFUSED title.** It
is not a title we refused; it is not a title. This inflated
`rule_excluded` by **5,764 records** (9,583 → 3,819). *Caught by:* the
paired-shape invariant, which failed on exactly one of the 20 pairs with
a single extra empty item. No fixture would have found it.

**D20 — refusals were ordered wrongly.** A 61-character table header
reported `too_long`, sending a reader after a truncation problem that is
not there. The more specific diagnosis now wins.

**D21 — `run_measurements.py` had drifted into lying.** Its heading says
"FROZEN SCREEN — the denominator", but once P1 started claiming lines it
was reporting post-claim residue under that heading. The screen total is
now counted independently of the unparsed list. *This is the failure the
whole project is built to avoid, in our own instrumentation.*

**D22 — my own meta-test caught my own new test.** The D9 guard (no
assertion may be a bare `A or B` disjunction) fired on a test written
this cycle. It works on the author who wrote it.

## 6b. Open items carried into P2

1. **`corpus2` does not ship the agent-shape originals**, so §5.2's
   paired-shape invariant (190 pairs) cannot run from the shipped
   package alone — P0 tested regions against `wave1 n=20` instead. Fix
   before P6, when `source` is exposed: add `corpus2/agent_shape/` with
   the `ERROR TEXT` block removed, so the pairs exist without unsealing
   the label.
2. ~~The flat-log audit~~ **DONE, see §5c.** Pass 2 runs on the holdout
   at acceptance, same procedure.
3. ~~`claimed_outside_screen`~~ **DONE.** 12,289 lines in
   `corpus2/work`, dominated by form-B titles ending in a colon, which no
   screen rule proposes. Watched, not tuned.
4. ~~Plan §3's numeric-row wording~~ **DONE** in plan v2.1a.
5. ~~P1 will break `t_form_counts`~~ **It did, plus two more.** See §7 —
   and the number is 20, not the 19 I predicted from a looser pattern.
6. **PERFORMANCE.** 27 → 15.7 → 7.5 → 6.0 MB/s against a 5 MB/s budget.
   About 1 MB/s of headroom left. P3 must profile before writing, and if
   it does not fit, the honest move is to raise the budget with a stated
   reason rather than let it slide.
7. **`LABEL_DISTINCT_LIMIT = 50` is now frozen** with its measurement in
   the module. N=12 and N=20 were measured and rejected: they cost real
   per-cycle series (`CC` with 26 distinct values, `Residual` with 22)
   for a modest item-count gain, while the runaway cases sit far above
   (`Target left` 1,999).
8. **The diagnostic reach is 0% and stays 0% until P3/P4.** Any coverage
   claim made before then is a claim about reach, not about value.

## 6p. The final review — two passes

**Pass 1, straight through the package.** Five findings, all fixed:

**D47 — the tests were never registered.** Requirements §8 asks for
`tests/tst_<n>.py` registered in `run_all_tests.py`. The tests existed
from P0 and the runner never did, so anyone following the house
convention would have run the suite and found nothing. Written, plus
`tests/__init__.py`.

**D48 — the shipped CLI crashed on a pipe.** `log_extractor.py --help |
head` raised `BrokenPipeError` and printed a traceback. Guarded.

**D49 — `__pycache__` was shipping in the tarball.** Excluded.

**D50 — three named deliverables did not exist.** Plan §8 lists
`DESIGN.md`, `MEASUREMENTS.md` and `EXTRACTOR_REQUIREMENTS.md` v2.
`REVIEW_P0_P2` §4.2 flagged this at P2 and it stayed open through six
more phases. `DESIGN.md` and `MEASUREMENTS.md` are now written;
**requirements v2 remains the one outstanding deliverable** (see §11).

**Verified clean:** stdlib-only imports (`bisect`, `re`, `sys`); no
f-strings, dataclasses, walrus, type hints, keyword-only arguments,
`yield from`, `raise…from` or `nonlocal`; imports and runs from an
unrelated working directory.

**AN UNTESTED CLAIM, stated plainly:** the module has **never actually
been executed under `libtbx.python`**. Portability is asserted by a
grep-based test, which is evidence and not proof. First thing to check on
a machine that has it.

## 6q. Pass 2 — reviewed as the consumer, not the builder

A different question: not "is the code right" but **"does the output
support the claim the consuming project needs to make?"**

**The flagship complaint case, `find_reference_14`, GUI shape, end to
end:** program identified as `phenix.find_reference` from its banner;
diagnostic reach true; **18 exclusions in exactly 2 reason groups**, the
larger being 15 × `not x-ray and not computational`; every one carrying a
line number, and line 56 reads `Skipping 9GSD.A - not x-ray and not
computational`. The user's own structures — 9GSD, 9GSE, 9GSF, **9GSG** —
are all there by name with their line numbers. A human can check every
part of that claim in the file in under a minute. **That is the whole
purpose of the tool, and it works.**

**The pairing test — which requirements §9 calls "the strongest
acceptance test available".** Never run until now, because it needs
good/failed pairs of the same program and those only existed once
`corpus2` was built. 23 pairs available; on the first 10:

**0 of 10 have identical finding sets.** A failed autobuild loses its
cycles, decisions, control skips and completion record while keeping 9 of
13 sections. A failed resolve_cryo_em drops from 13 phases to 0. A failed
refine loses all 31 stages.

The requirement was: *"If a successful run and a failed run of the same
program produce the same findings, the extractor is not reading the thing
that matters."* It does not, and it is.

**The honest limit, from the same test:** `xtriage_1_err` differs from its
good mate in **1 of 9 channels** — sections 38 → 0, everything else zero
in both. It is a 7-line log. For the shortest failures the only signal the
extractor has is total absence of structure, which is a real finding but a
thin one.

## 6o. The P9 cycle — inventory, audit, holdout

**D45 — a banner form the working corpus does not contain.** 18 of the 87
held-out logs write `**Starting phenix.molprobity **`, the banner wrapped
in asterisks. The extractor abstained on every one. **This is exactly what
holding programs back is for** — a form I could not have invented and
would not have found by re-reading logs I had already read.

Fixed (one character class in one regex), which takes holdout coverage
32% → 53%. **The 32% is the number of record**; the 53% is in-sample,
measured after tuning to the set, and is quoted only so the size of the
miss is legible. A holdout is opened once.

**D46 — the flat-log audit's sample definition changed meaning between
its two passes.** Pass 1 selected `if not scan(text).unparsed`, which at
P0 was equivalent to "zero screened candidates" because no parser existed.
By P9 parsers claim most candidates, so the same expression selects "logs
where nothing remains unexplained" — a different, larger set. Run
naively it picked up four ordinary short failed runs. **The comparison the
audit exists to make would have been worthless.** Corrected to the literal
criterion, which reproduces pass 1's set exactly.

Result, with the corrected criterion: **2 of 93 holdout logs have zero
screened candidates, both the same FSC data-table class as pass 1.** Two
passes, 346 logs, 9 such files, all data tables. The failure mode the
review worried about — a headerless narrative log whose facts match
nothing — **does not occur in either corpus.**

**The unparsed inventory** (`UNPARSED_INVENTORY.md`) is the §7
deliverable: 54,441 unclaimed · 35,089 generic-only · 9,347 refused by 11
named rules, with the 20 commonest unclaimed shapes and a per-program
breakdown. **Reported, not tuned to zero.** The largest unclaimed classes
are numeric table rows (45,094) and bare rules (9,265).

## 6n. Defects found in the P8 REVIEW pass

**D43 — the corpus sampling was quietly unsound, and the guard for it did
not exist.** Three linked findings:

1. **The control runner reported against a red baseline.** At stride 12 it
   printed `baseline 125 passed, 2 failed` and then reported 47 control
   numbers anyway. A control result is a *difference* from the baseline, so
   every one of them was meaningless. The runner now **refuses to run**
   when the baseline is not clean, prints the failing tests, and exits 1.
   Verified by deliberately breaking a test and confirming the refusal.

2. **The stride was being applied to the 20-log wave-1 corpora**, not just
   the 253-log working set. At stride 40 that left one log and four tests
   failed for want of data rather than want of correctness. The stride now
   applies only to the working corpus.

3. **Six tests were existence claims** — "some record names a child", "all
   three refusal reasons occur", "the R-factor lines are captured" — which
   a sampled corpus can legitimately miss. They now always read the whole
   corpus. So is the identification coverage gate: **a gate you can pass by
   sampling is not a gate.**

The suite is now identical at strides 1, 12 and 100. It diverges at 253
(one log), which is the honest limit and is stated rather than papered
over.

**D44 — a dead parameter.** `find_agent_measurements(lines, regions)`
never used `regions`; the caller passed `structure.regions` and the
function ignored it, which reads as though provenance were being consulted
when it was not. Removed.

**New invariant:** every attached measurement must cite its parent stage
or cycle line — 9,433 of them carry no literal evidence of their own (the
name is derived), so the span verifier skips them, and without this they
were entirely unchecked. All 9,433 pass.

**Clean on probing:** the 76 mid-file quoted agent metrics carry 22
plausible names (Resolution, Completeness, Space Group…), all
`source=unknown`.

## 6m. The P8 cycle — measurements

**SIGN-OFF: NO — the screen is not extended. Tom's answer, recorded.**
31 of 75 stands as what an attached-only channel reaches. The arithmetic
below is kept so the decision stays auditable.

**The recommendation, which was accepted:** The 44
labels we do not recover are present in the log, in shapes the frozen
screen refuses:

```
    All-atom Clashscore : 3.56      indent 4, space before the colon
  Clashscore            =   3.56    an `=` assignment
```

| extension | new candidates | 75-label reach |
|---|---|---|
| frozen | 41,023 | 31 |
| + any indent, space before colon | +8,055 | — |
| + `Key = value` | **+56,396** | — |
| both | +64,451 | ~37 |

**It more than doubles the candidate pile (+157%) to move the score by
about six.** Those 56,396 `=` lines are mostly PHIL parameter echoes,
which the consuming project already reads from its own routing layer. Say
the word if you disagree; my recommendation is to leave the screen alone
and record 31 of 75 as what an attached-only channel reaches.

**D39 — the channel emitted nothing at all, silently.** The agent block
opens with the marker, then a `-----` rule, then the metrics, and closes
with `*****`. Treating any rule as the close ended the block on its
opening rule, so `find_agent_measurements` returned zero records on every
log. Caught immediately by a smoke test; it would have been invisible in
the aggregate.

**D40 — the span verifier was made to assert something untrue.**
`Measurement.evidence()` returned the metric name unconditionally, and an
attached measurement's name is DERIVED (`r_free` from `R/Rfree=a/b`) and
appears nowhere on the line. Only the agent block writes a literal
`Name: value`. Attached measurements now claim no literal evidence — the
stage or cycle they came from is already span-verified.

**D41 — a hardcoded `source="agent"` turned uncertainty back into a
claim.** 76 measurements in `corpus2/work` come from the agent's own logs
quoting their children mid-file, where the region source is `unknown`
precisely because we cannot tell the agent wrapping this log from this log
quoting an agent. Hardcoding undid the one rule `source=unknown` exists
for.

**D42 — the negative-control baseline was RED and I nearly reported the
numbers anyway.** At stride 12 the control runner's baseline was 125/2,
because two tests assert that a specific line EXISTS somewhere in the
corpus — and a sampled corpus can legitimately miss it. Every control
number in that run was measured against a red baseline. Existence tests
now ignore the stride. **The runner should refuse to report at all when
the baseline is not clean**; that is an open item for P9.

## 6l. Defects found in the P7 REVIEW pass

**D37 — a latent wrong-name hole, found before it bit.** The banner
pattern captured `(\S.*?)` — everything to end of line. No current
identification was affected, because inside the six-line window every
banner in both corpora is a bare name. But **the corpus proves the unsafe
shape exists**: `Starting AutoBuild with the command:` sits at line ~55
in four logs, and ligandfit writes `Starting CC of ligand as input to
map: 0.714` — a metric wearing a verb, and the very line §4.6 warns
about. Whole-line capture would have produced
`phenix.cc of ligand as input to map: 0.714`; the obvious repair, taking
the first token, would have produced **`phenix.cc`** — a wrong name, which
fails the hard gate.

Requiring a bare single token abstains on both. Verified to change
nothing measurable: 146/253 and 7/20, identical before and after. Control
C42 holds it.

This is the first defect found *before* it produced a wrong output, and
it was found by asking what the pattern could match rather than what it
did match.

**D38 — layer B was not computed last, though its comment said so.**
`identify_program` sat between the cycle parser and the completion parser,
with a comment claiming it ran after every layer-A parser. Nothing used
it, so behaviour was correct and the documentation was not. Moved to the
end of `scan()`, and a test now asserts the ordering in the source —
because the order *is* the enforcement: a parser cannot consult a program
name that does not exist yet.

**Clean on probing:** no blacklisted banner shadows a real one; no
identified name contains whitespace or exceeds 40 characters (26 distinct
names); four abstaining logs have a usable banner just beyond the
six-line window and are deliberately left abstaining rather than widening
it into wizard logs where a child's banner would be reachable.

## 6k. The P7 cycle — identification

**One signal, and two rejected after measurement.**

| corpus | named | abstained | **wrong** |
|---|---|---|---|
| corpus2/work n=253 | **146 (58%)** | 107 | **0** |
| wave1 agent n=20 | 7 | 13 | **0** |
| wave1 gui n=20 | 7 | 13 | **0** |

The signal is the positional `Starting <program>` banner, read through
the agent prefix, with two degenerate forms blacklisted: `Starting
phenix` with no name (15 logs) and `Starting libtbx.start_process` (11,
the launcher). Both abstain **while reporting that the signal fired** —
the difference between a fingerprint that fails loudly and one that rots
quietly.

**REJECTED 1 — the structural hypothesis in §4.9 fails.** "The invoking
program writes its own parameter block first", tested on the 107
banner-abstainers: **6 correct, 50 wrong, 51 with no root at all.** Roots
are not program names — `scaling` is xtriage's (27 wrong by itself) and
`data_manager` is the shared scope. Salvaging it needs a root→program map
from the PHIL masters, which is non-goal §2.3. The plan asked for this to
be tried before any vocabulary scheme; it was, and it died.

**REJECTED 2 — on principle, not on its numbers.** The program-emitted
`Found phil, …/<stem>.eff` record scores 13 correct / 1 wrong / 93 absent
on the corpus2 abstainers. But the GUI names the `.eff` after the job,
which is also what names the log file, so using it is **the filename by
proxy** — and §4.9 forbids the filename because `ai_analysis` derived the
program that way and then printed a confident report about a program that
never ran. On wave-1 it is 0 correct and 3 wrong: the same weakness,
showing. Declining a signal worth 13 logs is what the rule is for.

**Predictions: all six hit**, including the one that predicted its own
conformance failure (P7.3). Predicting a failure and then not quietly
fixing the test to avoid it is the part that mattered.

## 6j. The P6 REVIEW pass

**D36 — a stated fact was wrong and had been propagated for four
documents.** The planning note recorded that phaser writes
`EXIT STATUS: SUCCESS` **per module**, and that one log carried 59 of
them. That claim justified making `completion_records` a list rather than
a scalar, and it was repeated into the plan, the module docstring, the
P6 predictions and the last handoff.

**Re-measured across both corpora: `EXIT STATUS` appears EXACTLY ONCE in
each of the 24 logs that carry it.** The per-module marker is
`Finished: <timestamp>` — 59 in `a2u-globulin-mr__phaser_3`, 34 in
`AF_POMGNT2__phaser_6`, 29 in three more.

**The conclusion was right; the evidence cited for it was not.** A list is
still correct — 59 completion records in one log is exactly the case a
scalar field would have mangled. But the number came from a clustering
probe that counted several marker kinds together and attributed the total
to the one I happened to name. Corrected in the module, the plan and here,
and pinned by a test.

This is the failure the whole procedure is built to catch and it survived
four documents, because nothing re-derives a claim once it has been
written down. Worth raising in the procedure review: **a measurement cited
as justification should be re-run when the code that depends on it is
written**, not just when it is first taken.

**A design consequence to record, not a defect:** completion records claim
their lines, so `wall clock time: 2.69 seconds` is no longer also a
labeled value. Labeled values fell 23,825 → 23,555 and `generic_only`
35,890 → 35,153. The information is in the record's own text; the
specific-parser-wins rule is working as designed. Stated because a
consumer querying `labeled_values` for timing will not find it.

**Clean on probing:** `wall clock time:` never appears mid-sentence (the
unanchored search is safe); all 949 completion records are
`source=program` except 2 in the agent's own logs, and no `top_level`
record is attributed to a non-program source; the 5 non-standard
`Finished` lines are all `Finished (Job PredictAndBuild_0 cycle 2 at …)`,
correctly `unknown`; 74 distinct named children.

**Known limitation, recorded:** `Finished with search...` yields
`applies_to="search..."`, trailing dots included. Cleaning it would be
interpretation; the raw form is what the program wrote.

## 6i. The P6 cycle — completion records

**Six of seven predictions hit**, including 113 named-child records to
the unit. The seventh missed, a defect was found, and the re-measurement
came in **exactly on the predicted number**.

**COMPLETION RECORDS, `corpus2/work n=253`:**

| | top_level | named child | unknown |
|---|---|---|---|
| ok (n=230) | 131 | 92 | 605 |
| err (n=23) | 1 | 21 | 4 |

`top_level` requires positive evidence: the record is the last event in
the file and nothing but blank lines follows. `Finished with X` names its
child — the only positive attribution the corpus offers, found by reading
the shapes rather than assumed. Everything else is `unknown`, because a
child completing is not the run completing and in a wizard log we usually
cannot tell which we are looking at.

**The asymmetry worth carrying to the consuming project:** 21 of 23
failed runs carry no top-level completion record and 21 of their 26
records name a child. Stated as an observation, never as a rule — an
empty list is not a verdict, and non-goal 2.8 is enforced by a test that
asserts the API exposes no `success`/`failure` surface at all.

**D34 — two separate endings merged into one.** Proximity clustering
alone merged two `Job complete` lines two apart. A repeated marker kind
now starts a new event. *Caught by the phase's own fixture.*

**D35 — a child's ending swallowed the run's own.** In **26 logs**
`Finished with process_predicted_model` sat within six lines of the
program's `Job complete` / `usr+sys` / `wall clock` block; every marker
kind differed, so they clustered and the run's own ending was attributed
to the child. `Finished with X` is now self-contained. *Caught by the
missed prediction* — nothing else was looking, and the wrong answer was a
plausible attribution rather than an error.

Two clustering defects in one phase is the lesson: proximity clustering
has more edge cases than it looks like it has, and both were invisible to
a green suite.

## 6h. Defects found in the P5 REVIEW pass

**D32 — a refusal reason that misleads, on 490 of 1,075 lines.** Every
refused `Setting` line carried `setting_states_no_reason`. Measured, that
is the wrong diagnosis for nearly half of them:

```
572  Setting output.overwrite to True        an assignment, genuinely no reason
310  Setting default value of True for X     a narration -- assigns nothing
180  Setting up prediction keywords...       a gerund
 13  Setting macro_cycles=3                  an assignment, no reason
```

A reader auditing the ledger would go hunting for a missing reason on 490
lines that are not assignments at all. Now three named refusals
(`setting_states_no_reason`, `setting_line_assigns_nothing`,
`setting_up_is_not_an_assignment`), with a control holding the split.

**This is the third instance of the same failure** — D20 (a table header
reported as `too_long`), D31 (a stage summary reported as a section-title
refusal), now D32. The pattern: a refusal is easy to get *right* and easy
to get *named wrong*, and nothing about a wrong name shows up as a
failure. Worth watching in P6 rather than discovering again.

**D33 — the portability meta-test had a false positive.** `(?:=` in a
non-capturing group contains the characters `:=`, so the walrus check
fired on correct Python 2-compatible code. Same class as D30: a meta-test
with a false positive is worse than none, because the tempting fix is to
weaken the code it flags.

**Clean on probing:** no line carries more than one unparsed record; all
235 decisions are `source=program`; no decision line is also a section
title; no decision value or reason is mis-split.

## 6g. The P5 cycle

**Predictions: five of six hit exactly**, including 208 decisions and
1,075 reasonless settings refused — both to the unit.

**The miss was arithmetic, not modelling.** I predicted conformance would
reach 19 ± 0; it reached 18. From 16 passed and 3 failed, turning two
green gives 18, not 19 — the suite has 20 tests and one is a *skip*,
which is not a pass. I predicted the number I wanted instead of adding up
the one I had. Both substantive halves were right. Recorded because a
zero-band miss must be resolved rather than explained, and because "I
read both test bodies this time" was true and still not enough.

**THE DESIGN DECISION OF THE CYCLE, and it is a refusal.** 1,283
`Setting` lines exist in 133 logs; only **208 carry a reason**. The other
1,075 state a value and no branch:

```
Setting output.overwrite to True                      129
Setting macro_cycles to  3                             69
Setting default value of  True  for  thorough_denmod   41
Setting up prediction keywords...                      25  (a gerund)
```

Requirements §4.4 defines a decision as the program changing its own
configuration **and stating the reason** — announcing a branch. A value
with no reason is not a branch, so those 1,075 are reported with a named
refusal rather than fitted with a second pattern. Same precedent as the
234 unmatched skips and the subjectless control skip. **They are a
candidate v2 channel**, recorded with their measurement so the decision
is taken deliberately rather than on first sighting — the consuming
project wants `program_chosen` facts, and this is where they would come
from.

**Review probes, all clean:** no decision value or setting is mis-split
(4 distinct values, 14 distinct reasons, 0 odd); no decision line is also
emitted as a labeled value; no reason still contains a stray ` as `,
which would have meant a greedy split.

**The full refusal ledger at P5** (`corpus2/work n=253`) — every line the
extractor declined, with the rule that declined it:

```
label_exceeds_distinct_value_limit  4196
rule_is_a_banner_not_an_underline   1342
section_title_too_long              1427
section_title_is_a_table_header      805
setting_states_no_reason             585
setting_line_assigns_nothing         310
skip_matches_neither_form            234
setting_up_is_not_an_assignment      180
section_title_is_a_numeric_row       170
stage_row_is_the_table_summary        70
cycle_line_carries_no_metrics         28
```

9,347 refusals against 54,441 still unclaimed and 35,890 generic-only.
Nothing is dropped silently.

## 6f. Defects found in the P4 REVIEW pass

**D31 — the stage-table summary row was reported only by accident.** The
`end` row is deliberately not a stage (it carries r_free 0.4942, the
number that hides the finding). But it was *claimed* and emitted nowhere,
which is the silent-drop this project's whole ledger exists to prevent.
It surfaced anyway — because it happens to sit above a rule line, so the
SECTION parser refused it as `section_title_is_a_numeric_row`. A
coincidence, and a misleading diagnosis for anyone reading the ledger.
**77 such rows across the two corpora were relying on that accident.**
Now refused by name (`stage_row_is_the_table_summary`), with a corpus
test asserting every one of them carries that reason.

*How it was found:* a review probe asking "is any line claimed and
emitted nowhere?" — which initially reported 0 and looked clean. Checking
whether the probe was vacuous is what exposed it.

**Clean on probing** (recorded so the next reviewer need not repeat
them): no stage row has a value/column count mismatch; no stage
description leaked from an unrelated block; 35 distinct stage names, all
plausible; no stage or cycle is attributed to a non-program source; every
cycle metric key is one of `r_work`, `r_free`, `Built`, `Placed`,
`Resolution`.

**Corpus observation worth carrying forward:** 25 of 75 cycles across 17
logs are sentinels — a third of all cycle records say *no usable result*.
They cluster in autosol and autobuild. That is a finding for the
consuming project, not a defect here, but it says the sentinel rule earns
its place: without it, 25 nonsense R-factors of 999.90 would be quotable.

## 6e. The P3 cycle — and the first pre-registered predictions

**Predictions were written before any P3 code** (`PREDICTIONS.md`),
closing the gap `REVIEW_P0_P2.md` §4.1 recorded. Seven of eight hit
inside their bands:

| | predicted | measured |
|---|---|---|
| P3.1 complaint case | 18 exclusions, groups 15+3, 0 control skips | exact |
| P3.2 `Starting` yields no phase | 0 | 0 |
| P3.3 unmatched skips | **233 ± 0** | **234 — MISSED** |
| P3.4 diagnostic reach (ok) | 30 ± 10 | 33 |
| P3.4 diagnostic reach (err) | 2 ± 2 | 2 |
| P3.5 unclaimed drop | 1,000 ± 400 | 1,084 |
| P3.6 conformance | 15 ± 1 | 15 |
| P3.7 performance | ≥ 5.5 MB/s | 6.0 |
| P3.8 a line both phase and section | ≥ 1 log | 21 logs |

**The miss earned its keep.** The band was zero, so the exit rule forced
a resolution rather than an explanation. The single extra line is
`Skipping as CC_mask is < 1/2 min_cc` — a control skip **with no
subject**. I predicted with a loose proxy (`as` anywhere in the line);
the parser requires `<what> as <reason>` with both sides present.

Resolution: **the parser is right and the prediction was wrong**, and the
reason is an asymmetry nobody had written down. An exclusion may be
anonymous — a rejected candidate with no name is still a rejection. A
control skip may not: that line does not say *which* phase did not run,
and emitting `what=None` would invent a phase that never existed, the
exact failure the two-way distinction exists to prevent. Now stated in
the module and pinned by a test.

The deeper lesson is the one this project keeps re-learning: the loose
proxy I predicted with and the strict rule I implemented were two
representations of the same idea, and nothing forced them to agree.

**D30 — the D9 meta-test had a false positive.** It flagged
`assert x.reason == "not protein or too few residues."` because the word
`or` sits inside a *string literal*. A meta-test with a false positive is
worse than none, because the tempting fix is to weaken the assertion it
flags. String content is now stripped before the check.

## 6d. Defects found in the P2 cycle

**D26 — my own P1 optimisation was wrong, and only the equivalence test
said so.** The cursor-based section lookup is correct *only if every call
arrives in ascending line order*. P2 added a second pass over the file, so
the cursor was already at the end when the ledger loop restarted at line
1, and every lookup returned `None` for the whole file. Replaced with a
stateless bisect. **Nothing else would have caught this**: the wrong
answer was a plausible `None`, not a crash, and the suite was green.
Control C21 now holds it.

**D27 — the test runner hid its own failures. The worst of the cycle,
because it is a defect in the instrument.** `run_all_tests` caught only
`AssertionError`. A control that made a test raise `KeyError` aborted the
run at test 20 of 70, printed **no counts at all**, and the control
runner reported "suite did not run" — so two controls (C1, C17) looked
like untested features when in fact they were working. In a normal run,
any non-assertion bug would have silently hidden every later test. Every
exception is now a failure and the run continues. *Found by:* a negative
control behaving oddly, not by any test.

**D28 — a control was patching the wrong site.** C16's pattern
`if line_no in claimed:` occurs twice with identical indentation, and
`replace(..., 1)` took the first. It was disabling the claim check inside
the labeled-value parser while claiming to disable it in the ledger. Both
sites now have distinct patterns and separate controls (C16, C22).

**D29 — a test passed with the feature disabled because its fixture could
not exercise it.** `test_P2_a_section_title_is_not_also_a_labeled_value`
used `Processing files:` as the heading — which ends in a colon and is
therefore not a `Key: value` candidate at all, so the claim check it was
testing never came into play. The fixture is now a heading that IS a
candidate (`Space group: P 21 21 21`). This is the same failure class as
D10: a fixture that cannot fail.

## 6c. Defects found in the P1 REVIEW pass (after the suite was green)

**D23 — a section swallowed the agent's footer. A provenance defect.**
A section runs to the next heading or to EOF, so the last program-written
section of a wrapped log extended into the agent block. MEASURED: 16
cases across both corpora, including molprobity's `Summary` section
"containing" the agent's FINAL QUALITY METRICS REPORT. A consumer
grouping items by section would then attribute the agent's numbers to the
program — the exact failure requirements §4.8 exists to prevent. Sections
are now clipped at any provenance boundary. *Caught by:* a review probe
asking whether sections respect regions; no test was looking.

**D24 — the section lookup was quadratic.** `section_of` scanned the
whole section list per line: 389 sections × 11,450 items on the largest
log, and one log carries 419 sections. Replaced with a cursor walk, valid
because sections are produced in line order and do not overlap.

**D25 — `_is_numeric_row` was 38% of scan time.** It ran a regex per
token on every line the first three screen rules did not claim (42,426
calls on the largest log). Now short-circuits on a cheap first-character
test and exits as soon as the outcome is decided.

**Both optimisations are covered by a corpus-level equivalence test**
against naive reference implementations — 416,987 lines and 101,682
section lookups compared. Arguing that an optimisation is safe is not the
same as showing it.

## 7. THREE FROZEN NUMBERS — SIGNED OFF AND APPLIED (the W1 rule, first
exercise)

Plan §10 W1: a frozen conformance number moves only with the measurement
that justifies it, a suite version bump and Tom's sign-off. All three
were proposed with measurements, approved, and applied. The suite is now
**v2**, with this change log in its own docstring so the next reader
finds it before the numbers.

**W1-1 — `t_form_counts` asserts `sections == 15`; correct is 20.**
The prototype recognised only titles embedded in a rule. The Program
Template's own shape is a title on its own line over a rule, which
appears in 225 of 253 `corpus2/work` logs. With both forms, all 20
wave-1 logs carry sections. *Applied:* 15 → 20.

**W1-2 — `t_a_bare_rule_is_not_a_section` uses an ambiguous fixture.**
Its input is `"text\n" + "="*70 + "\nmore text\n"` and it asserts no
section. But that IS form B, and in `corpus2/work` a text line above an
`=` rule is a real heading **1,818 times** (`Starting job`, `OUTPUT
FILES`, `REFINEMENT`). The test means "a rule with nothing above it";
the fixture accidentally supplies something above it. *Applied:* dropped
the leading `"text\n"`. The assertion then tests what it names, and my
own suite has both that test and its complement.

**W1-3 — `t_both_corpus_shapes_yield_the_same_structure` compares raw
counts.** It asserts `len(a.sections) == len(b.sections)`. Requirements
§4.8 asks something narrower: that the shapes agree *once agent-added
content is excluded*. MEASURED across all 14 wrapped pairs: the delta is
exactly −1 every time, and in every case the extra section is
`source="agent"` — the agent's own metrics banner. Excluding
agent-sourced items, the shapes agree exactly. *Applied:* filters
`source != "agent"` before comparing, which is also plan §5.2's stronger
semantic form.

Applied as suite v2. Prototype regression against v2: 15/4/1 (was 16/3/1
against v1) — the tightened suite separates the prototype from the new
module more sharply, as it should.

## 8. Standing limitations

- I have seen all 253 logs in `corpus2/work`. Every number above is
  in-sample. Only `corpus2/holdout` (87 logs across 44 whole programs,
  plus 6 error logs) is out-of-sample, and it opens once.
- 23 error logs in the working set, 2 of which yield diagnostic
  structure. Every claim about failed runs rests on that.
- The `err` label is a lower bound: `autosol_1_gene-5-mad` sits in the
  good set and ends in `STOPWIZARD`.
- One Phenix build, one user, two machines.


---

## 11. Deliverables — complete

`EXTRACTOR_REQUIREMENTS_v2.md` is written. It carries the four contract
additions, the corrected §4.7 (31 of 75 on GUI-shape, and why the
agent-shape 73 is circular), the corrected §4.9 (the structural
hypothesis tested and failed; the filename-by-proxy signal declined), the
three new non-goals, the screen at v2 with the decision not to extend it,
the restated identification target, and a change log with v1's numbers
left attached to wave 1.

**The two contract deviations are now DECLARED rather than silent**, each
with a test that fails if the code changes without the spec:

- **`Measurement.unit` is permanently `None` in v1.** Retained rather
  than removed, because the consuming project built against the stub.
  Populating it means parsing units out of values, which is
  interpretation (non-goal 2.10).
- **`Unparsed`'s v1 `kind` field is `screen_rule`.** `kind` names the
  *channel* on every other item and is what `items()` sorts by, so
  reusing it for "which screen rule matched" was a collision.

## 11b. House test integration (done)

The tests now run and fail **exactly like the samples**, and are registered
in `tests/run_all_tests.py`.

**Protocol matched:** `PROJECT_ROOT` on `sys.path`; `tests.tst_utils`
imported with an inline fallback; `run_all_tests()` calls
`run_tests_with_fail_fast()`, so **the first assertion failure raises an
uncaught exception with a full traceback** and the module exits 1. Verified
by breaking a test deliberately.

**A NAME COLLISION, found while registering.** `tst_log_extractor` is
**already registered** in `run_all_tests.py` — "Log Extractor Tests
(v113)", line 768. Two test modules of one name shadow each other on
`sys.path` with no warning, and whichever is imported first wins. Renamed
throughout:

    log_extractor.py            -> log_structure_extractor.py
    tests/tst_log_extractor.py  -> tests/tst_log_structure_extractor.py

registered as **"Log Structure Extractor"**. Worth checking what the v113
module tests — if it covers the same ground, one of the two should go.

**One deliberate deviation, and why.** The house protocol is fail-fast,
which throws away the magnitude of a failure. A negative control result is
the *difference* a disabled feature makes, and "C1 disabled → 24 passed, 10
failed" says more than "failed". So the module keeps a second entry point,
`tst_log_structure_extractor.py --count`, used **only** by
`run_controls.py`. `run_all_tests()` is unchanged house behaviour.

## 11c. Review of the converted tests

**D51 — a green tick with a third of the value skipped. The most important
finding, and it is the failure this project keeps making.** Without
`PHENIX_LOG_CORPUS` the module printed **"All 130 tests passed"** and the
house runner would show a green tick — while **34 corpus-level invariants
did not run**. Those are most of the value: fixture-only tests encode the
same assumptions as the code, which is how four defects survived 50
passing tests in the consuming project.

It is the third instance of the pattern: `tst_conformance` and the
prototype's own suite read different env vars, so the documented command
silently skipped 4 of 19; the ligandfit suite counted skips as passes and
reported 4/14 where the truth was 3/14.

Fixed: every skip registers itself and the module ends with a loud banner
naming the count and the two variables to set. Silent when everything ran.
Verified for no-corpus, partial-corpus and full-corpus runs.

**D52 — 232 lines of duplicated output.** Each test printed
`Test: <description>` and `  PASSED`; the house runner already prints
`  Running <name>...` and `  PASS: <name>`. Stripped. The 28 results that
carry a *measurement* were kept and reindented — those are the ones worth
reading (13,036 sections, 12,879 claims outside the screen, the full
refusal ledger, "1 log where splitlines() would have been wrong").

**D53 — four imported helpers were never used.** `assert_equal`,
`assert_true`, `assert_false` and `assert_in` were imported from
`tst_utils` with a fallback, and every test uses bare `assert`. Import
narrowed to `run_tests_with_fail_fast`, which is what is actually called.

**Checked and clean:** no absolute paths anywhere in the module, tests or
runners; `PROJECT_ROOT` resolves correctly from `langchain/tests/`;
fail-fast aborts at the first failure with a traceback and exits 1;
`--count` still emits the line `run_controls.py` parses; 130 tests and 47
controls green after the changes.

**Two things to know rather than fix:** the suite takes **~22 s** with the
corpus set, making it one of the slower modules in `run_all_tests.py`;
and `run_tests_with_fail_fast` sorts names in ASCII order, so the
uppercase-prefixed regression tests (`test_D11…`, `test_P4…`) run before
the lowercase ones.

**Samples could not be run here:** all three need `agent/` and
`tests/tst_utils.py`, neither in the bundle. `tst_error_analyzer.py` has no
fallback for `tst_utils` and dies on import; `tst_error_classifier.py` has
one for `tst_utils` but not for `agent.error_classifier`. The runner
handles both as `ImportError` and records them as failed, so nothing is
silently skipped.

**Worth your eye:** `agent/error_classifier.py`'s `classify_error(log_text,
program)` does regex classification over log text — the same treadmill this
extractor exists to replace. If it were rebuilt on `scan()` output it would
cite line numbers instead of matching patterns, which is the whole argument
for the tool.

## 12. What is left, and neither part is mine to do

**1. Put it in the tree and run it under `libtbx.python`.** The module
does not yet live under `libtbx/langchain/` and **has never been executed
under `libtbx.python`**. Portability is asserted by a grep-based test —
no f-strings, dataclasses, walrus, type hints, keyword-only arguments,
stdlib-only imports — which is evidence and not proof. Then run it
against a corpus you generate. Everything measured here comes from 346
logs, one Phenix build, one user, two machines: the assumption most
likely to be wrong and the cheapest to test.

**2. Take the holdout result to the consuming project's plan.**
Diagnostic reach is 30% on the working corpus and **6% on 44 unseen
programs**. An evidence layer built on stage tables and cycles helps on
wizard runs and has little to say about the long tail. If Phase 0a's gate
is to be measured across program families, that number belongs in the
plan before the arms are run, not after.



`EXTRACTOR_REQUIREMENTS.md` v2 is not written. It must carry, at minimum:

- the four contract additions (`source`, `identification`,
  `completion_records`, `labeled_values`) and `items()`
- §4.7 corrected: the label score is **31 of 75 on GUI-shape**, and the
  agent-shape figure is circular
- §4.9 corrected: the structural hypothesis was tested and **failed**
  (6 correct, 50 wrong); the filename-by-proxy signal was declined
- §2 extended with the three new non-goals
- the screen at version 2, and the decision not to extend it further
- a change log, with v1's numbers left attached to wave 1

Everything it needs is measured and recorded across `DESIGN.md`,
`MEASUREMENTS.md`, `PREDICTIONS.md` and this file. It is transcription,
not new work — which is exactly the kind of task that stays undone.
