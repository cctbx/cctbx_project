# Phenix log extractor — requirements v2

**Supersedes v1.** Delete v1; do not archive it in place.

Every number in v1 was measured against the frozen twenty
(`wave1`, md5 `9d8689b3861aeeb2c7b32585aad350c2`). Those numbers are
**left attached to wave 1** wherever they appear below. Where a v1 claim
was tested and found wrong, the correction is marked **CORRECTED** with
the measurement that replaced it.

---

## 0. Change log, v1 → v2

| | change | why |
|---|---|---|
| §2 | three non-goals added (2.8–2.10) | evidence found during the build |
| §3 | four fields added: `source`, `identification`, `completion_records`, `labeled_values`; `items()`; `section_id` on every item | the table could not express what §4.8 and §4.9 already required |
| §3 | `Unparsed(kind, …)` → `Unparsed(text, line, screen_rule, status, excluded_by)` | `kind` means the channel name on every other item; see §3.2 |
| §3 | `Measurement.unit` documented as permanently `None` in v1 | populating it is interpretation; see §3.3 |
| §4.1 | a second section form recognised | the Program Template's own shape; wave-1 sections 15 → 20 |
| §4.7 | **CORRECTED** — "67 of 75 verbatim" does not reproduce; the score is **31 of 75 on GUI-shape** | v1's figure was measured against text containing the block being scored |
| §4.9 | **CORRECTED** — the untested structural hypothesis was tested and **failed** | 6 correct, 50 wrong, 51 absent |
| §4.9 | coverage target moved off wave 1 to ≥40% on `corpus2/work` | the wave-1 target could not be met without fitting to its own test set |
| §6 | `unparsed` gains three statuses | one number is gameable from both ends |
| §7 | acceptance suite at v4 with a change log of its own | three frozen numbers moved, each with a measurement and sign-off |
| §9 | wave 2 landed: 253 working logs, 93 held out and opened once | |

---

## 1. What this is

Unchanged from v1. A deterministic reader for Phenix log files: text in,
structure out, every item carrying the line it came from. It reads what
the program says about itself. It does not interpret, rank, diagnose, or
call a model.

**Added in v2 — the organising principle**, learned when wave 2 arrived:

> **Headerlessness is a normal state, not a degraded one.** The program
> may be unknown, the section absent, the source unknown and the
> completion unassignable — and a locally supported fact can still be
> extracted and cited. A 40-line headerless log yielding three cited
> facts and no program name is a **successful** extraction.

Three layers, kept apart by the call order rather than by convention:
**A** generic observation, **B** identification, **C** program-specific
interpretation (out of scope). `identify_program` is the last thing
`scan()` does, so no layer-A parser *can* consult a program name.

`scan()` never requires identification, and failure to identify cannot
suppress any otherwise extractable structure. Enforced by a corpus
invariant over 253 logs and by a negative control.

---

## 2. Non-goals

v1's seven stand: it does not interpret or judge; does not rank or score
severity; does not read PHIL masters or defaults dumps; does not call an
LLM; needs no network and no third-party packages; does not guess when it
does not know; does not handle non-Phenix logs.

**2.8 It does not classify a run as failed.** Of 30 failed runs in
`corpus2`, **27 contain no error keyword at all** in what the program
itself wrote. A failure classifier built here would be wrong most of the
time while looking authoritative. `completion_records` reports records
observed; an empty list is not a verdict. Enforced by a test asserting the
API exposes no `success`/`failure` surface.

**2.9 It does not absorb `programs.yaml`'s metric vocabulary in v1.**
`measurements` carries only numbers arriving attached to a structure
already parsed. Enforced by a test that the module contains no metric
names.

**2.10 It does not interpret a labeled value.** No unit normalisation, no
type inference, no mapping of `R-free` / `R Free` / `r_free` onto one
name. That mapping is the treadmill; it belongs to the consumer.

---

## 3. The output contract

```python
scan(text, program_name=None) -> LogStructure
```

| Field | Contents |
|---|---|
| `sections` | `Section(title, line/start, end)` |
| `phases` | `Phase(name, line)` |
| `stages` | `Stage(name, line, metrics, description)` |
| `cycles` | `Cycle(number, line, metrics, sentinel)` |
| `decisions` | `Decision(setting, value, reason, line)` |
| `control_skips` | `ControlSkip(what, reason, line)` |
| `exclusions` | `Exclusion(item, reason, line)` |
| `measurements` | `Measurement(name, value, unit, line, context)` — §3.3 |
| `labeled_values` | **v2** `LabeledValue(label, value, line, repeat_count, end)` |
| `completion_records` | **v2** `CompletionRecord(text, line, applies_to)` |
| `unparsed` | `Unparsed(text, line, screen_rule, status, excluded_by)` — §3.2 |
| `identification` | **v2** `Identification(candidates, signals_fired)` |
| `forms`, `is_flat` | as v1 |
| `regions` | **v2**, diagnostic only, not part of the frozen contract |
| `claimed_outside_screen` | **v2**, the count that says the screen is too narrow |

Every item additionally carries **`source`** (`program` \| `agent` \|
`unknown`) and **`section_id`**, where `None` is normal.

Methods: `metric_moves()`, `exclusion_groups()` as v1, plus **`items()`**
(every item of every kind in line order) and **`reach()`** (three
metrics, always together — §7.2).

### 3.1 Non-negotiable properties

v1's five stand and are each a corpus-level invariant, not an assertion:
every item carries a line number; order is preserved; nothing is
deduplicated across kinds; `scan()` is pure; encoding-safe.

**Added in v2:** line numbers follow **newlines only**, as `grep -n` does.
`str.splitlines()` also breaks on CR, VT, FF, NEL and U+2028; one corpus
log carries 100 bare CRs and `splitlines()` gave it 539 lines where
`grep -n` gives 439.

### 3.2 `Unparsed`, renamed and extended

v1 specified `Unparsed(kind, line, text)`. **`kind` means the channel name
on every other item** (`"sections"`, `"cycles"`), so reusing it for "which
screen rule matched" is a collision. v2 names that field `screen_rule` and
adds `status` and `excluded_by` (§6).

### 3.3 `Measurement.unit` is permanently `None` in v1

The field is in the contract and is `None` on all 9,509 instances.
Populating it means parsing units out of values, which is interpretation
(§2.10). **The field is retained, not removed**, because the consuming
project built against the stub; it is documented as always-`None` rather
than quietly dropped.

---

## 4. Structures to recognise, with the evidence

### 4.1 Sections — TWO forms

v1 recognised only a title inside a rule (`wave1`: 15/20, 73 titles).
**v2 adds the Program Template's own shape**, a title on its own line over
a rule — 225 of 253 logs in `corpus2/work`, 2,109 distinct titles. With
both forms, `wave1` sections are **20/20**.

A bare rule with no title is a separator. Three refusals, each named:
title over 60 characters, title containing `|` (a table header), title
that is a numeric row. Underlines are `-` and `=` only; `*` rules are
banner delimiters, and accepting them added 824 junk titles.

A section stops at a **provenance boundary** — it must not run from
program-written text into the agent's footer.

### 4.2 Stage table

`wave1` 6/20; `corpus2/work` 72 headers in 47 logs. Rows are **anchored to
their header**, not matched globally: the row shape matches 40 lines in
`refine_5_beta_blip` where 29 belong to the table.

**The trailing `end` row is not a stage.** It carries the summary R-free
(0.4942) while the last real stage carries 0.4935 — and 0.4942 is exactly
the number that hides the finding. It is reported as a named refusal, not
dropped.

### 4.3 Cycle lines

`R/Rfree=a/b` is one key and two values. **999.90 is a sentinel** meaning
no usable result: the cycle is flagged and the nonsense R-factors are not
emitted. 25 of 75 cycles in `corpus2/work` are sentinels. A cycle line
with no metrics is a counter, not a record.

### 4.4 Decisions

`Setting X=Y as/because Z` — the program changing its own configuration
**and stating the reason**. `corpus2/work`: 1,283 `Setting` lines, only
**208 carry a reason**. The other 1,075 are reported under three named
refusals and are **not** fitted with a second pattern. They are the
natural home of a future `program_chosen` channel.

### 4.5 Skips — two different things

Unchanged, and confirmed: control skip (`as`/`because`) vs item exclusion
(dash form). The item name is optional. **New in v2 — the asymmetry:** an
exclusion may be anonymous; a **control skip may not**, because a skip
with no subject does not say which phase did not run. 234 `Skipping` lines
in `corpus2/work` match neither form and are reported.

### 4.6 Phases

`Running X` is a phase; `Starting X` is not. Confirmed at scale: 1,109
phase lines in 65 logs, 0 derived from a `Starting` line.

### 4.7 Measurements — **CORRECTED**

v1: "67 of 75 (89%) appear verbatim in program output. Reproducing the 67
is the target."

**That does not reproduce.** It was measured against text that contained
the metrics block being scored. Re-measured:

| | recovered |
|---|---|
| agent-shape | 73 of 75 — **circular, do not quote** |
| **GUI-shape** | **31 of 75** |

61 of the agent-shape 73 come from the agent's own block. **The target is
31 of 75 on GUI-shape**, and the labels remain a test, not patterns to
fit. Two labels are defective (`Residues Built` values that are the digit
string of an unrelated R-free) and are excluded.

The 44 uncovered values are in the log, in shapes the screen refuses
(`All-atom Clashscore : 3.56`, `Clashscore = 3.56`). **Extending the
screen was measured and declined:** +64,451 candidates on 41,023 to move
the score by about six.

### 4.8 What the agent adds

Unchanged in intent. **`source` now has three values, not two** —
`unknown` is real, because a metrics block quoted mid-file by the agent's
own log cannot be attributed. The agent prefixes the program's first line
(`LOG TEXT: Starting AutoBuild`) rather than adding a line; that line is
**program-sourced**, and in 9 of 20 wave-1 logs its payload is the best
identification signal available.

### 4.9 Which program wrote this log? — **CORRECTED**

**The untested hypothesis was tested and failed.** v1: "the invoking
program writes its own parameter block first". On the 107 banner
abstainers in `corpus2/work`: **6 correct, 50 wrong, 51 with no root at
all.** Roots are not program names — `scaling` is xtriage's (27 wrong
alone), `data_manager` is the shared scope. Salvaging it needs a
root→program map from the PHIL masters, which is non-goal §2.3.

**What works: one positional signal.** `Starting <program>` on the first
lines, read through the agent prefix, as a bare token, with two degenerate
forms blacklisted (`Starting phenix` with no name, 15 logs; `Starting
libtbx.start_process`, 11). Decoration is tolerated —
`**Starting phenix.molprobity **` is a form the holdout revealed.

| corpus | named | abstained | wrong |
|---|---|---|---|
| `corpus2/work` n=253 | 58% | 107 | **0** |
| `wave1` n=20 | 7 | 13 | **0** |
| **holdout n=87, unseen programs** | **32%** (see §9) | 59 | **0** |

**A signal was declined on principle.** `Found phil, …/<stem>.eff` scores
13 correct / 1 wrong / 93 absent on the abstainers — but the GUI names the
`.eff` after the job, which also names the log file, so it is the filename
by proxy. v1's flat prohibition on filenames stands.

**Gates, restated.** HARD: no log misidentified, abstention always
allowed — unchanged, and met on 44 unseen programs. TARGET: v1 asked for
≥10 of 20 on wave 1; that cannot be met without fitting a vocabulary to
the same 20 logs it is scored on, which v1 itself flags as guaranteed by
construction. **v2 target: ≥40% named on `corpus2/work`** (met at 58%).
Wave 1 keeps a regression floor of 7.

**The real fix remains outside this tool:** `Starting phenix` with no name
appears in 15 logs. One identifying line from the Program Template would
end the fingerprinting problem permanently.

---

## 5. Known traps

v1's list stands and every one was met in practice. **Added in v2:**
`str.splitlines()` is not newline splitting (§3.1); a refusal that is
right with a *wrong name* never shows up as a failure (three separate
defects); an optimisation must be proved output-identical over the whole
corpus before it is kept.

---

## 6. Reporting what is not understood

`unparsed` is a requirement, not a nicety. **v2 gives it three statuses,
always reported together**, because one number is gameable from both ends:

- `unclaimed` — screened, no parser claimed it
- `generic_only` — recorded as a labeled value, understood as nothing in
  particular
- `rule_excluded` — refused by a **named** rule, reported with the rule

The candidate screen is **frozen** and deliberately over-selective, so
that `unparsed` cannot be reduced by seeing less. It moves only with a
recorded measurement and sign-off; it is at **version 2** (single-character
labels admitted, which recovered 115 previously invisible `R: 0.26
Rfree: 0.29` lines across 33 logs).

`corpus2/work`: 54,441 unclaimed · 35,089 generic-only · 9,347 refused
across 11 named rules. Full breakdown in `UNPARSED_INVENTORY.md`.

---

## 7. Acceptance

### 7.1 The suite

`tst_conformance.py`, black-box, at **version 4** with its own change log.
Three frozen numbers moved, each with a measurement and sign-off:
sections 15 → 20; the bare-rule fixture corrected; the both-shapes test
excludes agent-sourced items; the identification tests read the contract
field names; the coverage target moved to `corpus2`.

**Result: 20 passed, 0 failed, 0 skipped**, both corpus shapes set.

### 7.2 Three reach metrics, always together

Structural reach without diagnostic reach is a flattering number.

| | structural | basic | diagnostic |
|---|---|---|---|
| `corpus2/work` ok n=230 | 97% | 96% | 33% |
| `corpus2/work` err n=23 | 96% | 39% | 17% |
| holdout n=87 | 98% | 93% | **6%** |
| fresh raw captures n=288 | 100% | 100% | 38% |

**~2.5% of both corpora are data files, not logs** (two-column FSC
tables). Every reach metric carries that floor and cannot close it.

### 7.3 Method required of the build

v1's three stand — tests before implementation with real counts, a
corpus-level invariant per feature, a negative control per feature — and
all three held. **Added in v2:** predictions pre-registered before each
phase with tolerance bands, and a missed prediction forces a documented
model change or a defect entry. Three of the most valuable findings came
from missed predictions.

---

## 8. House conventions

Unchanged. **Two are not yet met and are recorded as outstanding:** the
module does not yet live under `libtbx/langchain/`, and **it has never
been executed under `libtbx.python`** — portability is asserted by a
grep-based test, which is evidence and not proof.

---

## 9. Corpus

| | logs | note |
|---|---|---|
| `wave1` agent / GUI | 20 / 20 | frozen, hashed |
| `corpus2/work` | 253 | 230 ok, 23 err |
| `corpus2/holdout` | 93 | 44 unseen programs + 6 err |

**The holdout was opened once**, against predictions hashed
`ea3124b77d9a…` beforehand. Result: **0 misidentified** on 44 unseen
programs; identification 32%; basic reach 93%; **diagnostic reach 6%,
which breaches the plan's stated 15% floor**.

That failure is reported, not softened. Its meaning: diagnostic structure
lives in the wizard programs, and unseen long-tail programs mostly have
none to find. **This is a fact about Phenix's long tail, and it bears
directly on the consuming project** — an evidence layer built on stage
tables and cycles will help on wizard runs and have little to say about
the other ~60 programs.

**Standing limitation:** one Phenix build, one user, two machines.
