# FLAT_LOG_AUDIT — pass 1, before any parser exists

> **RECORD — not edited after the fact.** This is the audit as run before any parser existed. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

**Why now.** Plan §3.1: a log can carry facts and match no screened shape at
all, and broadening the screen until such lines match would defeat the
anti-gaming design. The audit is the alternative test. It only means
something *before* parsers exist, because afterwards I would be reading these
logs already knowing what I built — so P0 is the last honest moment.

**Sample, fixed before anything was read:** every log in `corpus2/work n=253`
with **zero** screened candidates. That is 7 logs, small enough to read in
full rather than sample.

**Pass 2** runs on `corpus2/holdout` at acceptance. The comparison is the
result; a single audit on logs I have already seen is my judgment about my
own corpus.

---

## The 7 logs

| File | Lines |
|---|---|
| `AF_7mjs_H_PredictAndBuild__fsc_half_maps.masked.mtriage.log` | 225 |
| `AF_7mjs_H_PredictAndBuild__map_model_fsc.log` | 296 |
| `AF_7n8i_L_PredictAndBuild__fsc_half_maps.masked.mtriage.log` | 104 |
| `actin_denmod__map_model_fsc.log` | 651 |
| `apoferritin_denmod__map_model_fsc.log` | 826 |
| `model-building-scripting__map_model_fsc.log` | 3526 |
| `rotavirus-autosharpen__map_model_fsc.log` | 68 |

## Finding

**All seven are the same thing, and none of them is a log.** Each is a
two-column numeric data table, most with a `1/resolution    CC` header and
then nothing but pairs of floats to EOF:

```
    1/resolution    CC
    0.046604188     0.679858408
    0.076736282     0.619271119
```

They are FSC curves written to a `.log` extension. There is no narrative, no
banner, no heading, no `Key: value` — nothing for a *log* reader to read,
because they are not records of what a run did.

**They match no screened shape for a correct reason:** the `numeric_row` rule
requires ≥4 numeric tokens and these rows have 2. Widening it to 2 would
match an enormous amount of ordinary prose and text-plot output across the
other 246 logs.

## What this does and does not tell us

**It does not confirm the failure mode the review was worried about.** The
hypothetical was a headerless narrative log — `Reading model foo.pdb / There
are 134 residues` — that carries facts and matches nothing. **No such log
appears in `corpus2/work`.** The zero-candidate class here is a different
thing entirely.

**It does surface a corpus-composition fact worth carrying forward:** 7 of
253 files (2.8%) are data tables rather than logs. They should be *reported
as such*, not counted as extraction failures, or every reach metric carries a
2.8% floor it can never close and the number quietly misleads.

**It leaves a real question open, and the honest answer is "not yet".** A
two-column FSC curve *is* extractable content — a named series with a header.
Emitting it would mean deciding that a resolution/CC table is a measurement
series, which is program-specific interpretation and out of scope for v1
(plan §2.9, §2.10). Recorded here rather than fixed, so that if the holdout
pass finds the same class, the decision is made deliberately with two
measurements behind it instead of on the first sighting.

## Actions

1. **P9 reach metrics report the data-file count separately.** A file with no
   narrative content is not a log the extractor failed to read.
2. **The screen does not move.** Two-column tables stay out of it.
3. **Pass 2 on the holdout**, same procedure, sample fixed before reading.
4. **Open for Tom:** should `phenix.map_model_fsc` output be in a *log*
   corpus at all? If the GUI points `ai_analysis` at these files, that is a
   product question worth knowing about independently of this tool — the
   honest answer for such a file is "this is a data table, not a run log",
   and something has to say it.
