# FLAT_LOG_AUDIT — pass 2, on the holdout

> **RECORD — not edited after the fact.** This is the audit as run on the holdout. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

Pass 1 (`FLAT_LOG_AUDIT.md`) ran before any parser existed, on the
working corpus. This is the comparison it was written for.

**Sample, fixed by the same rule as pass 1:** every log with **zero
screened candidates**.

---

## A defect in the audit itself, found before the result

Pass 1 selected its sample with `if not scan(text).unparsed`. At P0 that
was equivalent to "zero screened candidates", because no parser existed
and every candidate was therefore unclaimed. **By P9 it is not
equivalent**: parsers now claim most candidates, so the same expression
selects "logs where nothing remains unexplained" — a different and much
larger set.

Run naively at P9 it selected 4 error-holdout logs that are ordinary
short failed runs, not flat logs at all:

```
Starting phenix.predict_and_build
on Wed Aug 12 15:16:16 2026 by terwill
==============================================================================
Processing files:
------------------------------------------------------------------------------
  No files found
```

Those have candidates in abundance; they are simply all claimed. **The
sample definition had silently changed meaning between the two passes**,
which would have made the comparison worthless. Corrected to the literal
criterion — zero lines matching the frozen screen — which reproduces pass
1's set exactly (7 of 253) and is what pass 2 uses below.

## The result

| | logs | zero screened candidates |
|---|---|---|
| working set (pass 1) | 253 | 7 (2.8%) |
| **holdout (pass 2)** | **93** | **2 (2.2%)** |

The two holdout logs:

```
AF_7mjs_H_PredictAndBuild__fsc_half_maps.unmasked.mtriage.log
AF_7n8i_L_PredictAndBuild__fsc_half_maps.unmasked.mtriage.log
```

**Same class as pass 1, exactly.** Two-column FSC data tables written to
a `.log` extension — `1/resolution  CC` and then floats to EOF. Not
narrative logs. In pass 1 the masked variants of these same two datasets
appeared; here the unmasked ones do.

## What this settles, and what it does not

**Settles:** the failure mode the review was worried about — a headerless
*narrative* log carrying facts that match no screened shape — **does not
occur in either the working set or the holdout**. Two passes, 346 logs,
9 zero-candidate files, all of them data tables. The frozen screen is not
missing a class of readable log; it is correctly declining files that are
not logs.

**Does not settle:** whether `phenix.map_model_fsc` and `mtriage` FSC
output should be in a log corpus at all. It is the same product question
pass 1 raised. The honest answer for such a file is "this is a data
table, not a run log", and something still has to say it.

**Recorded as a limitation:** 2.2% of the holdout and 2.8% of the working
set are data files. Every reach metric in this project carries that floor
and cannot close it. Stated wherever reach is quoted.
