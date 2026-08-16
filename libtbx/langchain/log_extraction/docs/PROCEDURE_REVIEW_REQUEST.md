# Review request: is our development procedure working?

> **RECORD — not edited after the fact.** This is the procedure review request as sent. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

**What this is.** Not a code review. A request for criticism of the
*procedure* we are using to build a small deterministic tool, prompted by
an observation: this cycle produced an unusually high number of caught
mistakes. We want an outside reading of whether that means the procedure
is working well or badly, and what to change.

**What we need from you.** Concrete suggestions. Section 7 lists the
questions we cannot answer from inside, but do not feel bound by them.

---

## 1. Background

### 1.1 The thing being built

A deterministic reader for PHENIX crystallography log files. Text in,
structure out: sections, phases, stages, cycles, decisions, skips,
exclusions, labeled values, completion records — each carrying the line
number it came from. It does not interpret, rank, diagnose, or call a
model.

It exists because a consuming project needs *checkable* evidence about
what a program run actually did. That project's hard gate is that fewer
than 5% of its claims may be incorrect, because a user cannot verify a
claim about a runtime value they never saw. So the extractor's job is not
to be clever. A missing structure costs a finding; a wrong one costs the
project.

The tool is one Python module, standard library only, roughly 600 lines
so far, with a test suite and a set of negative controls.

### 1.2 The people and the cadence

One human (a crystallographer, author of much of the software this reads)
and an LLM assistant doing the construction. Two other LLMs have been
used as external reviewers of plans.

The cadence is **think → act → review → present**, one cycle per phase,
with a handoff document written each cycle so work can resume after a
context loss.

Standing rules, established before this project from earlier work
together:

- Tests written and shown failing before the implementation, in the same
  cycle.
- Real pass/fail counts. "Tests pass" without numbers is not a result.
- At least one corpus-level invariant per feature — because fixture-only
  tests encode the same assumptions as the code.
- At least one **negative control** per feature: disable it, show the
  suite fails, report both numbers. Without the second number the first
  only shows that the tests agree with themselves.
- **Measure before asserting.** Several plausible ideas have died on
  contact with a five-line measurement.
- Revisions recorded in docstrings, not made quietly.
- Superseded documents deleted, not archived, because a stale duplicate
  once gave a wrong answer about which step came next.

### 1.3 The evidence base

- **Wave 1:** 20 logs, 9 programs, frozen, with every number in the
  requirements stated against its md5. All from a wizard-heavy subset.
- **Wave 2:** 347 original logs from ~107 programs, arriving mid-project.
  Split into a 253-log working set, a sealed holdout of 87 logs across 44
  whole programs, and a sealed answer key.

Wave 2 invalidated a headline claim from wave 1 ("every log yields
structure; none is flat" — true for 20/20 in wave 1, true for 69% in wave
2), and forced a plan rewrite before any code was written.

### 1.4 Where we are

Phase P0 is complete: types, line discipline, wrapper-region detection, a
frozen "candidate screen", a CLI, derived views — and deliberately **no
content parsers yet**. 45 tests pass, 0 fail. 12 negative controls, all of
which move the numbers. Nine further phases are planned.

---

## 2. The defect record

This is the material for your judgment. Every defect found since the
project began, with **which step caught it**.

### 2.1 Found before any code (planning and corpus work)

| # | Defect | Caught by |
|---|---|---|
| A1 | A ground-truth label file recorded two garbage values: `Residues Built: 2336488640192294` is the digit string of an unrelated R-free from elsewhere in the same log | re-measuring a stated number |
| A2 | The claim "67 of 75 labels appear verbatim in program output" was measured against text that *contained* the labels. Honest figure ~64/75 | re-measuring |
| A3 | The corpus-stripping script leaked 486 lines of wrapper into 25 of 29 test logs | comparing against originals that arrived later |
| A4 | **Our own prior conclusion was wrong**: we had reported that failed runs carry the program's own remedy lines. They came from the leaked wrapper. True figure: 3 of 30 failed runs carry any error keyword | new data arriving |
| A5 | The "GUI-shape" corpus still carried an agent prefix on line 1 of 15 of its 20 logs | diffing two corpora |
| A6 | The acceptance suite would **fail a correct extractor** — it asserts a section count that a correct fix changes from 15 to 19 | measurement while planning |
| A7 | Two test suites read different environment variable names, so running with the documented one silently skipped 4 of 19 tests and printed a clean result | reading the code |
| A8 | Plan v1 sequenced the build on an assumption ("logs are structured") that the new corpus contradicts | external LLM review |
| A9 | A proposed scalar `termination` field was incoherent, and wrong in substance: 17% of logs carry more than one completion record, one carries 59 | external review, then measurement |
| A10 | An external reviewer's two recommendations, adopted together, would have silently defeated each other (a generic capture channel would drive the anti-gaming metric to zero by construction) | reviewing the review |
| A11 | A reviewer's proposed filter would have deleted the scarcest signal in the corpus (`Sorry:` lines in program output) | measuring before adopting |

### 2.2 Found during P0 construction

| # | Defect | Caught by |
|---|---|---|
| D1 | A rule's branch could never fire — the condition was computed after a loop that guaranteed it was always zero | its own fixture test |
| D2 | The fix then misjudged short inputs | same fixture, after D1 |
| D3 | A corpus invariant was wrong, not the corpus | corpus test |
| D4 | **A feature was implemented and never tested.** Forcing it off changed nothing and the whole suite still passed | negative control |

### 2.3 Found in the first review pass, after the suite was green (32/32)

| # | Defect | Caught by |
|---|---|---|
| D5 | `bool` is a subclass of `int`, so `True` was silently accepted as line 1 | manual probe |
| D6 | The module's own docstring example used Python-3-only syntax, in a file required to run under Python 2 | grep |
| D7 | The frozen contract names a field `start`; the implementation called it `line` | contract cross-check |
| D8 | **The worst of the cycle.** The wrapper does not add a line at the top — it *prefixes* the program's first line. Marking that whole line as wrapper-written buried the program's own identification banner in 9 of 20 logs | targeted measurement |
| D9 | A test assertion of the form `assert A or B`, where B was separately asserted on the next line — it could not fail | reading the tests |
| D10 | Two invented test fixtures asserted the opposite of what they claimed (a hand-written table row given three numeric cells where the rule requires four). Both would pass a careful reading and fail the data | running the fixture |
| D11 | The implementation deviates from the plan's wording — and **the code is right**. Measurement showed the plan's reading drops 11,363 genuine table rows | measurement |
| D12 | A negative control's patch target had gone stale and it silently reported "target not found" | the control runner |

### 2.4 Found in the second review pass, after the suite was green again (36/36)

| # | Defect | Caught by |
|---|---|---|
| D13 | `str.splitlines()` breaks on bare CR, VT, FF, NEL and U+2028 as well as newline. One corpus log has 100 bare CRs: we reported 539 lines where `grep -n` reports 439, so **every cited line number after the first CR pointed a human at the wrong place** — in a tool whose entire value is that a claim can be checked | probe |
| D14 | The wrapper-prefix rule required a header above it, so the commonest case in the shipped corpus (15 of 20 logs) went unmarked | probe |
| D15 | The span verifier split lines differently from the extractor — a verifier checking against the wrong ruler | the new D13 test |
| D16 | A region test counted "not a header" as a footer, so when a third region kind was added the reported number jumped from 14 to 29 and the test passed while printing a wrong number | reading output |

### 2.5 Tally

**28 defects.** Rough classification:

- **In the code:** 13
- **In the tests or controls:** 6 (D9, D10, D12, D15, D16, and A7)
- **In our own written claims or plans:** 9 (A1, A2, A4, A6, A8, A9, A10, A11, D11)

By catching mechanism:

| Mechanism | Caught |
|---|---|
| Targeted measurement (during planning or review) | 9 |
| Manual review reading code, tests or output | 7 |
| Corpus or fixture tests during construction | 4 |
| External review | 3 |
| Negative controls | 2 |
| New data arriving | 2 |
| Reviewing a review | 1 |

---

## 3. What we think this shows

Offered so you can disagree with it.

**3.1 A green suite is a weak stopping signal.** Twelve defects were
found *after* the suite was green, eight of them in deliberate review
passes. If we had stopped at "45 tests pass" we would have shipped a tool
that cites wrong line numbers.

**3.2 Our tests are as error-prone as our code, and nothing tests the
tests.** Six defects were in test or control code. Two of those (D9, D10)
are invisible to a green suite *by definition*: a test that cannot fail
reports success. The negative controls are the only structural defence,
and D12 showed the controls themselves rot.

**3.3 Defects cluster at representation boundaries.** Every one of D7,
D11, D13, D15 and A5 is a place where two representations of the same
thing had to agree and nothing forced them to: contract prose vs field
name, plan prose vs regex, text vs lines, verifier's line map vs
extractor's line map, corpus A vs corpus B. That is a pattern we did not
design against.

**3.4 The cheapest mechanism found the most.** Nine defects fell to a
measurement that took under five minutes. This is already a standing rule
and it is still the most productive thing we do.

**3.5 Negative controls found what nothing else could.** Only two
defects, but D4 — a feature implemented and never tested — was invisible
to a green suite and to reading. No amount of care would have found it.

**3.6 One defect we did not catch, and could not have.** A4 was a wrong
conclusion we had already reported. It was corrected only because more
data arrived. Nothing in the procedure would have caught it otherwise.
That is the failure mode we are least protected against and most worried
about.

---

## 4. What the procedure currently does not do

- No pre-registered review checklist. Each review is shaped by what was
  just built, so it looks where the builder already suspects.
- No test-of-tests tier beyond the negative controls.
- No boundary contract: nothing requires that two representations of the
  same fact be asserted equal.
- No distinction between self-review and adversarial review.
- No measure of how long a defect survived before detection.
- No rule that a claim written in a handoff document must carry the
  measurement that supports it — though in practice most do.

---

## 5. Constraints on any suggestion

- The builder is an LLM without memory across sessions; handoff documents
  are the continuity mechanism. Anything requiring sustained attention
  across many cycles must be written down or it will not survive.
- The human reviewer is expert in the domain (crystallography, this
  software) and has limited time. His attention should go where domain
  judgment is required, not where mechanism can substitute.
- The corpus is large but seen. Only the sealed holdout is out-of-sample,
  and it can be opened once.
- Cost matters. A procedure that doubles the cycle time needs to catch
  more than it costs.

---

## 6. The counter-argument we want tested

Our defect count went up sharply this cycle. Two readings:

**(a) Construction is getting worse.** More code, more haste, more
mistakes.

**(b) Detection is getting better.** The negative controls, the corpus
invariants and the deliberate review pass are new or newly systematic.
The defects were always there; we are now seeing them.

We believe (b) but we cannot distinguish them from inside, because we
have no baseline defect rate from before these mechanisms existed. If you
can suggest a way to tell (b) from (a) — even a crude one — that would be
the single most useful thing you could give us.

---

## 7. Questions

1. **Is 28 defects in this much work alarming or expected?** What would
   you have expected for a 600-line module with this much specification
   churn?
2. **How do we test our tests?** Is the negative-control mechanism
   sufficient, or is there a better structural answer to D9 and D10 (a
   test that cannot fail; a fixture that asserts the opposite of its
   claim)?
3. **Should review be pre-registered?** A checklist fixed before the
   cycle starts, so the review cannot be shaped by what was built. What
   should be on it?
4. **Is think/act/review/present the right granularity?** Would two
   review passes of different character — one self-review, one
   adversarial with a different reader — be better than one thorough one?
   We effectively did that this cycle by accident.
5. **How do we protect against A4** — a wrong conclusion, already
   reported, that no test will contradict because it is a claim about the
   world rather than about the code?
6. **Should every representation boundary get a mandatory agreement
   test?** It would have caught 5 of 28. It also adds tests nobody reads.
7. **What is the right stopping rule for a review pass?** We currently
   stop when we run out of ideas, which is not a rule.
8. **Are we over-documenting?** Each cycle produces a handoff, and the
   defect record above is itself a document. There is a real risk that
   the documentation becomes the work.

---

## 8. If you want to look at the artefacts

- `EXTRACTOR_BUILD_PLAN_v2.1a.md` — the plan, with its own weakness
  register (§10) and a disposition table for both external reviews (§9)
- `p0/log_extractor.py` — the module; defects are documented in place, in
  the docstring of the thing they affected
- `p0/tests/tst_log_extractor.py` — 45 tests, including one regression
  test per defect in §2.3 and §2.4
- `p0/run_controls.py` — 12 negative controls
- `p0/HANDOFF.md` — the current-state document
- `p0/FLAT_LOG_AUDIT.md` — an audit run deliberately before parsers
  existed, so it could not be shaped by what was built
