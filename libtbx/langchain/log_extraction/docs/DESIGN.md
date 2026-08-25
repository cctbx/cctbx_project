# DESIGN — Phenix log extractor

For whoever maintains or extends this. Plan §8 names it a deliverable; it
was written last, in the final review, which is later than it should have
been.

---

## 1. Shape

One module, `log_extractor.py`, standard library only (`bisect`, `re`,
`sys`). About 1,800 lines including the comments that carry the
measurements. Tests in `tests/tst_log_extractor.py`, registered in
`run_all_tests.py`. Negative controls in `run_controls.py`. Corpus
measurements in `run_measurements.py`.

Plan W11 set a threshold: past ~1,200 lines of *code*, or when a parser
needs its own state machine, revisit the package question. Neither has
happened — the file is long because roughly half of it is recorded
measurement, and that is deliberate.

## 2. Three layers, kept apart by construction

**A — generic observation.** Needs no header, no program identity, no
section. Sections, phases, skips, exclusions, stage tables, cycles,
decisions, labeled values, completion records, measurements.

**B — identification.** Reads whatever evidence exists, may return
unknown.

**C — program-specific interpretation.** Out of scope. Named here so
nothing in A or B quietly grows into it.

**The enforcement is the call order**, not a convention: `identify_program`
is the last thing `scan()` does, so no layer-A parser *can* consult a
program name — it does not exist yet. A test asserts that ordering in the
source, and a corpus invariant asserts that `scan(text)` and
`scan(text, program_name=…)` produce identical layer-A items across all
253 logs.

## 3. One pass, one ledger

`scan()` splits the text once (on newlines only — see §5), computes the
wrapper regions, then runs each parser over the line list. Every parser
returns `(items, claims, refusals)`:

- **claims** — line numbers this parser understood
- **refusals** — line numbers it saw and declined, each with a NAMED rule

After the parsers, the ledger walks every line once and files each
screened candidate as `unclaimed`, `generic_only` or `rule_excluded`.

**Specific parsers ignore other parsers' claims; the generic channel
respects them.** A heading that is also a phase is two true facts and both
are emitted. A heading also captured as a nameless "fact" is one fact
counted twice, so `labeled_values` skips claimed lines.

## 4. Adding a parser

1. **Measure the shape first**, on `corpus2/work`, and write the counts
   into the comment above the parser. Every existing parser carries its
   own evidence; a pattern without a measurement is a guess.
2. Write the tests, and run them **failing**, before the implementation.
3. Return `(items, claims, refusals)`. **Refuse by name.** Three separate
   defects in this project were a refusal that was right with a reason
   that was wrong, and a wrong reason never shows up as a failure.
4. Add at least one **corpus-level invariant** — fixture-only tests encode
   the same assumptions as the code.
5. Add a **negative control** in `run_controls.py`: disable the feature,
   confirm the suite fails, report both numbers. If a control changes
   nothing, the feature is untested — that is how D4 was found.
6. If it is an optimisation, prove it **output-identical over the whole
   corpus** before keeping it. An optimisation that changes behaviour is a
   defect wearing a performance argument, and P2 proved that is not
   hypothetical.

## 5. Decisions that look arbitrary and are not

- **Newline-only splitting.** `str.splitlines()` also breaks on CR, VT,
  FF, NEL and U+2028. One corpus log has 100 bare CRs: `splitlines()` gave
  it 539 lines where `grep -n` gives 439, so every cited line after the
  first pointed a human at the wrong place.
- **The frozen candidate screen.** Four shapes, deliberately
  over-selective, and it does not move without a recorded measurement and
  sign-off. It exists so that `unparsed` cannot be reduced by seeing less.
- **`unparsed` has three statuses.** A single number is gameable from both
  ends; `generic_only` is the honest successor and is expected to be large
  (35,089).
- **`source` includes `unknown`.** Marking an ambiguous region `program`
  turns uncertainty into a claim.
- **`section_id=None` is normal.** Otherwise headerless programs become
  second-class inputs and implicit sections get invented.
- **The `end` row of a stage table is not a stage.** It carries the summary
  R-free that hides the trajectory.
- **999.90 is withheld, not reported.** A quotable nonsense R-factor is
  worse than a gap.
- **Identification abstains rather than guesses**, and reports which
  signals fired. A wrong program name is worse than none.

## 6. Using it in other contexts

- `scan(text)` is pure: no file I/O, no global state, no logging, and two
  scans of the same text are equal. Callers decode with
  `errors="replace"`.
- `main()` is the only I/O, and never gates on a `.log` extension.
- `items()` returns every item of every kind in line order, so a consumer
  with no program name and no sections can still group by proximity.
- `reach()` returns the three metrics together, deliberately: structural
  reach without diagnostic reach is a flattering number.
- Runs under `python` and `libtbx.python`: `from __future__` imports, no
  f-strings, dataclasses, walrus, type hints or keyword-only arguments,
  checked by a test. **Never actually executed under `libtbx.python` —
  that is an untested claim and the first thing to verify on a machine
  that has it.**

## 7. What is deliberately absent

No per-program metric vocabulary. No PHIL master or dump reading. No
interpretation, ranking or severity. No failure classification — the
extractor reports completion records observed and never a verdict,
because 27 of 30 failed runs carry no error keyword at all in what the
program itself wrote.
