# Note for the Phase 0a plan — what the extractor's reach implies

**For:** the log-analysis redesign (`phenix.ai_analysis`), Phase 0a.
**From:** the log structure extractor, now built and validated.
**Why now:** Phase 0a's gate is scored **per program family**, and a family
that yields nothing is excluded from Phase 1. The size of that excluded set
is measurable **before** the arms are run, and it is larger than the frozen
20-log corpus can show.

---

## 1. The one number

The plan measured, on the frozen 20: **EVENT findings 9/20, T1 MEASUREMENT
findings 17/20, union 20/20.**

Measured on 253 production logs across **63 programs**:

| | logs | event findings |
|---|---|---|
| the plan's **9 families** | 143 | **46%** |
| **everything else** | 110 | **13%**, across 43 programs |
| all | 253 | 33% |

**48 of the 63 programs — 159 logs, 62% of the corpus — yield no event
finding at all.** Not because the extractor fails on them: because they do
not tabulate anything. Most are 20–300 line tools that announce a result
and stop.

Per family, in the plan's own terms:

| program | logs | event | measurement | union |
|---|---|---|---|---|
| autosol | 9 | **100%** | 100% | 100% |
| map_to_model | 7 | **100%** | 100% | 100% |
| autobuild | 20 | **90%** | 100% | 100% |
| refine | 33 | **87%** | 93% | 93% |
| resolve_cryo_em | 9 | 77% | 77% | 77% |
| predict_and_build | 6 | 33% | 50% | 50% |
| **xtriage** | 28 | **0%** | 96% | 96% |
| **phaser** | 19 | **0%** | 100% | 100% |
| **real_space_refine** | 14 | **0%** | 85% | 85% |
| **model_vs_data** | 14 | **0%** | 100% | 100% |
| **molprobity, mtriage, polder, process_predicted_model, …** | | **0%** | 57–100% | |
| map_model_fsc | 5 | 0% | **0%** | **0%** |

## 2. What this changes about the plan

**(a) The per-family rule bites harder than the 20-log corpus suggests.**
"A 0% family does not fail the project but EXCLUDES that program from Phase
1 scope" was written when the corpus had 9 families. At production scale it
excludes **three of the four commonest programs by log count** — xtriage
(28 logs), phaser (19), real_space_refine (14) — on the *event* axis.

**(b) It does not follow that those families are lost**, and this is the
important half. Their **measurement** reach is 85–100%. xtriage yields no
event finding and 96% measurement; phaser 0% and 100%. The plan already
knows findings come in two families; **at scale the split is not incidental,
it is the shape of the whole corpus.** For most programs the product is
measurements, not events.

**(c) One family yields nothing on either axis.** `map_model_fsc`: 0% and
0%, 5 logs. Those are two-column FSC data tables with a `.log` extension —
**not logs at all**. ~2.5% of any Phenix corpus is this. Worth excluding by
recognition rather than letting it score as a failure.

## 3. What I would put in the plan, concretely

1. **State the event/measurement split as the expected shape, not a
   fallback.** The gate says "≥1 correct materially-useful claim Arm B
   misses". For xtriage, phaser, real_space_refine and most of the long
   tail, that claim will be measurement-derived. If the scoring rubric
   implicitly expects event-shaped findings, three high-volume families
   score 0 for the wrong reason.

2. **Score the gate on a corpus that contains the long tail.** The frozen 20
   has 9 families, all main-line. The pooled ≥70% measured there does not
   predict the pooled figure at production scale, where 43 of 63 programs
   are single-purpose tools. I am not proposing to unfreeze the corpus —
   predictions measured against it stay valid — but the **Phase 1 scoping
   decision** needs the wider number.

3. **Decide the `map_model_fsc` class explicitly.** Recognise-and-decline is
   honest; scoring it as a failed family is not.

4. **Expect the failed-run picture to differ by input shape.** In GUI-shape
   logs the failure went to stderr and was stripped: 5 of 253 carry any
   terminal marker. In raw captures **27 of 288 (9%)** carry a traceback,
   the program's `Sorry:` diagnosis and its remedy. The extractor now reads
   these. **Which shape Arm B is handed changes what there is to find** — the
   same lesson as the agent-header correction, one level down.

## 4. What the extractor now supplies, per finding family

**Event:** stage trajectories with the metric that moved and where
(`R-free 0.4880 → 0.4915 at 1_xyzrec`), sentinel cycles flagged as *no
usable result*, skipped phases with the program's stated reason, exclusion
groups with counts and line numbers, decisions with their reasons, and the
program's own account of stopping — traceback, `Sorry:` and the remedy.

**Measurement:** every `label: value` the run reported, uninterpreted, with
line numbers; `--label X` gives one label's series across a run.

Everything cites a line so a human can check it. Nothing is interpreted,
ranked or diagnosed — that stays the LLM's job, which is what the
findings-first redesign asked for.

## 5. What I am not claiming

- These are **availability** numbers — whether a finding of that shape
  exists — not whether it is **material**. Materiality is a property of
  (finding, run state) and is the scoring question, not this one.
- 33% is not a ceiling on the product. It is a ceiling on *event* findings
  from what programs write today. The union with measurements is 93%.
- One Phenix build, one user, two machines. 253 logs.
