# Review of the work so far — P0 through P2

> **RECORD — not edited after the fact.** This is the state of the build at P2. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

Not a summary. A deliberate look across the three cycles for things that
are wrong, missing, or true-but-flattering. Written after re-running
every probe below, not from memory.

**Current state:** 72 tests pass, 0 fail. 22 negative controls, all
moving the numbers. Wave-1 conformance 13/6/1 against suite v2. Six of
eleven output channels are still empty by design.

---

## 1. What is solid

Checked, not assumed:

- **The contract surface is complete.** All eleven item channels, all
  four contract additions (`source`, `identification`,
  `completion_records`, `labeled_values`), plus `regions`,
  `claimed_outside_screen`, and every method — `metric_moves`,
  `exclusion_groups`, `items`, `reach`, `unparsed_counts` — exist and
  respond. Nothing in plan §1 is unimplemented.
- **Purity, encoding-safety, span verification, hint-independence** all
  hold at corpus scale (253 logs, 138,543 items).
- **Nothing is deduplicated across kinds**, as requirements §3 demands: a
  labeled value is also a `generic_only` unparsed record, deliberately.
- **The refactor discipline held.** Every optimisation this cycle was
  proved output-identical over the whole corpus before being kept — the
  quadratic-claim fix below changed 0 of 138,543 items.

## 2. What the review found wrong

### 2.1 The screen silently drops R-factors — the most valuable find

`_SCREEN_KV_RE` requires `[A-Za-z]` followed by `{1,48}` more characters,
so **a single-character label cannot match**. MEASURED on
`corpus2/work n=253`: 124 such lines, 115 of them `R:`, and they look
like this:

```
R: 0.26  Rfree: 0.29
R:   0.27  Rfree:   0.31
```

The R-factor and R-free pair is the single most diagnostically important
measurement in crystallography, and the extractor cannot currently see
it in these lines.

There is a second, quieter problem underneath. The screen's own comment
says "label 1-48 chars". `[A-Za-z]` plus `{1,48}` is **2 to 49**. The
implementation has never matched its documented intent — the same
prose-versus-code mismatch as D11, in a different file.

Three readings measured:

| reading | candidates | vs frozen |
|---|---|---|
| frozen `{1,48}` → 2–49 chars | 40,899 | — |
| `{0,47}` → 1–48 chars (matches the prose) | 40,686 | **−213** |
| `{0,48}` → 1–49 chars (adds short labels, keeps long) | 41,023 | **+124** |

Note the trap: the reading that matches the prose **loses 213
candidates**, because it also tightens the upper bound. Only the third
adds without subtracting.

**The screen is frozen, so this needs sign-off** — the module says
changing it changes every unparsed number ever measured. *Proposal:*
adopt `{0,48}`, and correct the comment to say 1–49. Known limitation to
record either way: `R: 0.26  Rfree: 0.29` is two facts on one line, and
v1 will capture it as `label="R", value="0.26  Rfree: 0.29"`. Splitting
compound lines is interpretation and belongs later, if at all.

### 2.2 A quadratic loop in the labeled-value parser

The claim loop sat *inside* the per-distinct-value loop, re-walking every
occurrence once per distinct value — O(distinct × occurrences) to produce
a result identical to one pass. A label with 50 distinct values over
1,000 lines did 50,000 dictionary writes to claim 1,000 lines. Fixed;
verified output-identical across all 253 logs; 6.0 → 6.3 MB/s.

### 2.3 Two universal invariants existed only on paper

Plan §5.2 lists **order preservation** and the **performance budget** as
universal invariants. Neither had a test. The budget lived in a
measurement script — a budget nothing enforces is a hope, and this
project has spent three cycles watching that number fall (27 → 15.7 →
7.5 → 6.3 MB/s). Both are now tests.

## 3. What is true but flattering, and must be said with the caveat

**Structural reach is 97% and diagnostic reach is 0%.** Six of eleven
channels — phases, stages, cycles, decisions, control skips, exclusions,
measurements, completion records — are empty in every log. Verified by
sampling 85 logs: non-empty counts are 0 across all of them. The tool
currently reports *that a log has headings and labeled facts*, not *what
the run did*. Every coverage number before P3/P4 is a statement about
reach, not value.

**`generic_only = 35,766` is the honest half of that.** Those lines are
recorded and citable but not understood as anything in particular. Under
a single `unparsed` number they would have shown up as progress.

## 4. Procedural gaps

**4.1 Predictions have not been pre-registered.** Plan §4 requires each
phase to record "a pre-registered prediction with a tolerance band before
the phase runs", and `PREDICTIONS.md` is a named deliverable in §8.
**This has not been done for P0, P1 or P2.** Every number reported so far
was measured after the fact. That does not make them wrong, but it
removes the check the plan built in: a measurement that can only be
compared to itself cannot surprise me. P3 should start with a written
prediction before any code.

**4.2 Three documentation deliverables are unwritten**:
`EXTRACTOR_REQUIREMENTS.md` v2 (which should carry the corrected §4.7
label measurement and the four contract additions), `DESIGN.md`, and
`MEASUREMENTS.md`. The handoff has absorbed some of their content, which
is the drift the plan warns about — one current-state document is the
rule, but that rule assumed the others existed.

**4.3 The negative controls now sample.** Controls run at a corpus stride
of 5 because the runner would otherwise take nine minutes. Every reported
number is still full-corpus, but a control seeing a fifth of the corpus
is a weaker control, and this was a cost decision made under time
pressure rather than a designed one. Worth a second opinion.

## 5. Risks carried into P3

1. **Performance.** 6.3 MB/s against a 5 MB/s budget. P3 adds three
   parsers. If it does not fit, the honest move is to raise the budget
   with a stated reason, not to let it slide.
2. **The diagnostic channels are where the domain judgment lives.** P0–P2
   were mechanical: a heading is a heading. `Skipping X as Y is set`
   versus `Skipping 9GSD.A - not x-ray` is the distinction the consuming
   project's complaint case depends on, and getting it wrong is not
   visible as a crash.
3. **Test-suite runtime.** 16 s and rising; the scan cache bought one
   cycle of headroom, not more.
4. **`claimed_outside_screen` is 12,289 and unexamined.** It is doing its
   job — form-B titles end in a colon and no screen rule proposes them —
   but nobody has looked at *which* lines are in it. If a later parser
   starts claiming lines the screen never saw, that number is the only
   warning, and an unexamined warning is not a warning.

## 6. What I would do next, in order

1. Sign-off on §2.1 (the screen change), because it is the only item that
   touches a frozen number and it blocks nothing else.
2. Write `PREDICTIONS.md` for P3 **before** writing P3.
3. P3: phases, control skips, item exclusions — with the two-way skip
   distinction as its own corpus invariant from the first commit, not as
   something checked at the end.
