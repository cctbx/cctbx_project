# MEASUREMENTS

Every number this project asserts, with the corpus it was taken against and
the command that produces it.

**The rule that produced this file:** a claim cited as justification gets
re-run when the code depending on it is written, not only when it is first
taken. D36 — a wrong "phaser writes EXIT STATUS per module" claim — survived
four documents because nothing re-derived it.

**Regenerated after every change that moves a number.** If a figure here
disagrees with another document, this one is right and the other is stale.

## Current-state documents

These must be true now, and are reconciled whenever a number moves: this
file, `EXTRACTOR_ARCHITECTURE.md`, `EXTRACTOR_REQUIREMENTS_v2.md`,
`DESIGN.md`, `HANDOFF.md`, `INSTALL.md`, `UNPARSED_INVENTORY.md`.

The development **records** — pre-registered predictions, the P0–P2 and
flat-log audits, the procedure review, the test-tier and improvement
plans — were kept alongside these during development and were **never
edited after their date**: a prediction rewritten after the fact is not
a prediction, and a review edited to match what was later built is not a
review. They were removed from the repository before release; where one
quoted a figure that has since moved, the current value is here.

**Corpora**

| name | logs | shape | note |
|---|---|---|---|
| `wave1 agent` / `gui` | 20 / 20 | agent-wrapped / stripped | frozen, md5 `9d8689b3861aeeb2c7b32585aad350c2` |
| `corpus2/work` | 253 | GUI | 230 ok + 23 err, the working set |
| `corpus2/holdout` | 93 | GUI | **spent** — opened once at acceptance |
| fresh raw captures | 288 | mixed | sent later; never used to build anything |
| smoke (in tree) | 19 | mixed | 1.25 MB, set-cover over 47 features |

**Reproduce**

    phenix.python langchain/tests/run_all_tests.py            # tiers 1+2
    phenix.python langchain/log_extraction/validate.py        # tier 3
    phenix.python langchain/log_extraction/run_controls.py    # 54 controls
    phenix.python langchain/log_extraction/run_measurements.py

## Reach and identification

| corpus | logs | structural | basic | diagnostic | identified |
|---|---|---|---|---|---|
| corpus2/work ok | 230 | 97% | 96% | **33%** | 56% |
| corpus2/work err | 23 | 96% | 39% | **17%** | 96% |
| fresh raw captures (see note) | 288 | 100% | 100% | **38%** | 24% |
| smoke corpus (in tree) | 19 | 95% | 95% | **53%** | 63% |

*The fresh set is **288 logs**, not the 291 quoted in earlier documents:
Tom's bundle held 305 files, of which 3 are macOS `._` resource forks and
14 duplicate the corpus. The 291 was 288 + the 3 non-logs. Reach figures
are unaffected — the `._` files were always skipped — but the headline
count was wrong and is corrected here.*

*`corpus2/holdout` measured **6% diagnostic** when opened once at
acceptance. It has since been re-read, so any number from it now is
in-sample. **These 288 are the replacement fresh set.***

## Item counts, corpus2/work n=253

- `sections` **13028**
- `phases` **1109**
- `stages` **1012**
- `cycles` **75**
- `decisions` **208**
- `control_skips` **98**
- `exclusions` **49**
- `measurements` **9509**
- `labeled_values` **23747**
- `completion_records` **854**
- `terminal_failures` **22**
- `unparsed` **98879**

## The unparsed ledger, corpus2/work

- unclaimed **54449**
- generic_only **35083**
- rule_excluded **9347**

## Channels

**12**: `sections`, `phases`, `stages`, `cycles`, `decisions`, `control_skips`, `exclusions`, `measurements`, `labeled_values`, `completion_records`, `terminal_failures`, `unparsed`

Diagnostic: `stages`, `cycles`, `decisions`, `control_skips`, `exclusions`, `terminal_failures`

## Fixed numbers of record

- **Suite: 150 tests, 0 failed. Controls: 54, all discriminating.**
- **Wave-1 conformance (suite v4): 20 passed, 0 failed, 0 skipped.**
- **Full validation: 10 of 10 checks**, 148,590 items verified, 151 named /
  **0 wrong**.
- **Holdout, opened once** against predictions hashed `ea3124b77d9a…`:
  0 misidentified · identification 32% · basic reach 93% · diagnostic reach
  **6%** (breached the plan's 15% floor, reported as a failure).
- **75 agent metric labels, scored on GUI-shape: 31 recovered.** The
  agent-shape figure of 73 is circular and must not be quoted.
- **Screen extension declined** (Tom, P8): +64,451 candidates on 41,023 to
  move the label score by ~6.
- **`LABEL_DISTINCT_LIMIT = 50`**, chosen against N=12 (141 groups
  collapsed) and N=20 (54).
- **Identification signals:** the `Starting <program>` banner and the
  `PHENIX <program> <weekday>` header. Across all 657 distinct logs held:
  **44% named, 1 disagreement with the filename — and that one is the
  filename being wrong** (`ligand_identification_68.log` is a GUI job that
  runs LigandFit and says so three times).
- **Terminal failures:** 25 of the 288 fresh logs carry a `Traceback`, 25 a
  `Sorry:`, 19 end `Please try again.` — 27 with some terminal marker (9%),
  against 5 of 253 (2%) in the GUI-shape corpus, where the failure went to
  stderr and was stripped into the agent wrapper.
- **Compound `label: value` pairs:** 9% of key:value lines across 554 logs
  carry two or more.
- **Performance:** a MEASURE, not a gate — ~5.3 MB/s idle, 2.5 MB/s on a
  loaded machine. The gate is a collapse floor of 1.5 MB/s.
- **Data files, not logs:** 7 of 253 working and 2 of 93 holdout have zero
  screened candidates, all FSC two-column tables. Every reach metric carries
  that ~2.5% floor.
