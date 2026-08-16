# Test-tier plan v2 — revised against review

> **RECORD — not edited after the fact.** This is the tiering plan as agreed. Some figures
> here have since moved; the current values are in
> `MEASUREMENTS.md`, which is regenerated from the corpus. A
> pre-registered prediction rewritten afterwards is not a
> prediction, and a review edited to match what was later built is
> not a review.

**Supersedes v1.** The review is adopted in full except where a measurement
changed the answer; those are marked.

> **Tier 1 answers "does every known case still work?"
> Tier 2 answers "does the extractor still behave coherently on diverse
> real logs?"
> Tier 3 answers "are our claims about the full known population still
> true?"**
>
> Three different jobs, not three sizes of the same test.

---

## 0. What the review changed

| # | point | disposition |
|---|---|---|
| 1 | keep Tier 2; the space is worth it | **adopted** |
| 2 | do not claim Tier 2 "preserves what the corpus tier is for" | **adopted, wording replaced** |
| 3 | select the smoke set by a coverage matrix, not channel count | **adopted — and it made the core set SMALLER, §3** |
| 4 | drop the "1-in-10 chance" statistic | **adopted, it was wrong** |
| 5 | fixtures as excerpts with source hash, generated mechanically | **adopted** |
| 6 | keep both the fixture and the invariant that found it | **adopted as the stated philosophy** |
| 7 | classify each test `unit` / `smoke_property` / `full_population` / `full_measurement` | **adopted** |
| 8 | strengthen paired shapes from counts to semantics | **adopted — it was a real defect, §1** |
| 9 | Tier 3 visible but explicitly requested; absent corpus = failure | **adopted** |
| 10 | a machine-readable validation record | **adopted** |
| 11 | full corpus outside the source tree | **open — your call, §8 Q1** |
| 12 | if 21 MB is fine, do the simple thing | **open, but my figure was wrong, §8 Q2** |

---

## 1. Two review points were real defects, not style

**§8 was a genuine bug in the test, and I had documented it wrongly
twice.** The paired-shape invariant built semantic lists and then asserted
only `len(a) == len(b)`. Two logs can each yield 17 items and disagree
about all 17. `EXTRACTOR_ARCHITECTURE.md` and requirements §4.8 both
described it as semantic; it was not. **Fixed:** it now asserts full
semantic identity and names the first divergent item. Verified — all 20
pairs pass strictly.

**§4 was wrong arithmetic on my part.** "A defect in 1 of 253 logs has
roughly a 1-in-10 chance of appearing in the smoke set" holds for a random
sample. The smoke set is deliberately curated, so detection probability is
high for known structural classes and **zero** for an unknown one — D45 is
the proof. Replaced with: *Tier 2 substantially reduces input diversity,
19 logs against 293; there is no useful statistical detection probability
because the set is non-random.*

## 2. What Tier 2 does and does not prove

Adopting the review's wording verbatim, because mine was too strong:

> **Tier 2 is a real-data structural smoke test that runs the same
> invariant machinery on a deliberately diverse subset. It detects broad
> regressions and previously known failure classes. It does not validate
> population-level claims.**

## 3. The smoke set, chosen by set cover — and the core is smaller than v1

I built the coverage vector the review asked for and solved it rather than
hand-picking. **47 distinct features** across all 293 available logs: every
channel (11), every screen rule (4), every refusal rule (11), every source
mode (3), agent / GUI / error shape, identified / abstained, sentinel
cycle, anonymous exclusion, collapsed label, multiple completion events,
named-child and top-level completion, wrapper regions, a section-less
item, CR-vs-newline, quoted agent block, tiny and huge logs.

**A greedy cover needs 8 logs, 434 KB, and hits 47 of 47.**

Adding the flagship acceptance cases (refine trajectory, failed build,
cryo-EM branch, xtriage) and the GUI twin of every agent log — the paired
invariant needs both members — gives:

**19 logs, 1.25 MB, 5.8% of the full corpus.**

Larger than v1's 672 KB, and better justified: v1's was picked by "20
programs, all channels", which the review correctly called necessary but
not sufficient. A manifest ships beside it —

```
log                                why_in_smoke
agent/find_reference_14.log        complaint case; anonymous exclusion; 2 groups
work/ok/fem_407.log                bare CR newlines (D13)
work/ok/sec17-sad__ai_agent_2.log  quoted agent block mid-file; source=unknown
agent/autobuild_4_bromodomain.log  sentinel cycle; control skip; decision
agent_gui/<each of the above>      the paired-shape invariant needs both members
...
```

— so a change to the set is reviewable, and a Tier-3 test asserts the
cover is still complete.

## 4. Test classification

Every test carries one label, so a smoke result can never be read as a
corpus claim:

| label | runs | example |
|---|---|---|
| `unit` | always | every defect regression; parser rules; contract shape |
| `smoke_property` | always | items cite matching text; ledger balances; optimised == naive; paired shapes identical; **every screen rule fires** |
| `full_population` | Tier 3 | identification 0 wrong; every log scans; the 75-label score; coverage gates |
| `full_measurement` | Tier 3 | performance over 17.7 MB; every counted total |

**Adopting review §7's best suggestion:** the smoke set is *constructed*
so that every screen rule fires, which promotes "screen rules all fire"
from a full-corpus coverage claim to a smoke property. The set cover
already guarantees it.

## 5. Fixtures — excerpts, hashed, generated

```python
Fixture(text=..., source="find_reference_14.log",
        start_line=54, end_line=58, source_md5="...")
```

A minimal **real excerpt**, not the triggering line alone, because parsers
depend on neighbouring lines and on parser state. `log_fixtures.py` is
**generated from a manifest plus the corpus** and committed; regeneration
is reproducible and Tier 3 verifies every excerpt against the frozen
source. That removes the human transcription step that produced D10 twice
— an invented fixture asserting the opposite of what it claimed, which
passed a careful reading both times.

## 6. The philosophy, stated because it is the justification

> **Fixtures preserve discovered examples; invariants preserve the general
> property that discovered them.**

So a defect gets both: a Tier-1 fixture pinning the exact case, *and* the
invariant that found it running at Tier 2 every time and Tier 3 before
release. D31 is the model — a summary row reported only by accident,
caught by an accounting identity, now pinned by a fixture *and* still
guarded by the identity.

## 7. Tier 3 — requested, never skipped

```
$ run_all_tests.py
...
  Log Structure Extractor    PASSED (4.1s)     [unit + smoke]
  Full validation: NOT REQUESTED
  Run: phenix.python langchain/log_extraction/validate.py
```

Not a skip and not a failure — explicit scope reporting. Once requested,
**a missing corpus is a hard failure.**

`validate.py` emits a record, so "I think I ran it" and "this revision
passed this frozen corpus" are distinguishable:

```
extractor sha256: ...
corpus id:        9d8689b3861aeeb2c7b32585aad350c2 (wave1) + ...
date:             ...
293/293 scanned · all gates PASS · 47/47 controls moved
```

## 8. What I still need from you

**Q1 — where does the 21 MB corpus live?** The review prefers an external
test-data area over the source tree and I agree in principle. **I do not
know cctbx's conventions and will not invent a mechanism.** If there is a
regression-data location, name it: `validate.py` will look there and fail
loudly with the corpus id if it is absent. If there is none and 21 MB
in-tree is acceptable, say so and I will put it beside the package.

**Q2 — my earlier figure was wrong, which weakens the "just do the simple
thing" option.** I said optimising the corpus tier would take 22 s to
"about 6 s". Measured properly: one shared scan pass over 253 logs is
**3.0 s**, and three tests genuinely need their own passes — purity
(compares two scans), identification-independence (`scan(text)` vs
`scan(text, program_name)`), and the naive-reference equivalence walk.
Realistic floor: **~10 s**, not 6.

So the simple option is ~10 s and 21 MB in-tree on every run; the tiered
option is ~4 s and 1.25 MB, with the full corpus behind one command.
**I now lean tiered**, where before I would have accepted either.

**Q3 — is 1.25 MB acceptable in the tree**, up from the 672 KB you were
asked about? It buys the flagship acceptance cases and complete agent/GUI
pairs. If not, the 434 KB set-cover alone still covers all 47 features and
loses only those two things.
