# PHENIX AI Agent - Changelog

## Version 112.79 (Remove directory scanning from graph perceive node)

### Problem

The `perceive()` node in `graph_nodes.py` called
`_discover_companion_files()`, which scanned the directories of ALL
files in `available_files` — including user-supplied input files —
to look for companion outputs (e.g., map coefficients MTZ alongside
a `_data.mtz`).  This meant the agent could pick up unintended files
from user input directories without being told to.

### Root cause

`_discover_companion_files()` was introduced when the session did not
have comprehensive output file discovery.  Since then,
`session.get_available_files()` gained three layers of output
discovery:

1. `_discover_cycle_outputs()` — tries stored paths, then scans the
   expected `sub_NN_program/` output directory
2. `_find_missing_outputs()` — finds companion refine/autobuild files
   from known output paths
3. Step 3 catch-all — scans agent sub-directories for any missed
   output files

These session-side mechanisms are scoped to agent output directories
only, making the graph-side `_discover_companion_files()` both
redundant and overly broad.

### Fix

Removed the `_discover_companion_files()` function from
`graph_nodes.py` entirely.  The `perceive()` node no longer performs
any directory scanning.  It still:
- Injects output files from history entries
- Filters intermediate/temp files via `_filter_intermediate_files()`

### Tests

- Removed 5 companion file discovery tests (tested the deleted
  function)
- Added `test_perceive_no_input_dir_scanning` — verifies the
  perceive pipeline does NOT pick up files from user input
  directories
- Added `test_history_injection_still_works` — verifies output
  files from history entries are still injected

### Documentation

- Updated `ARCHITECTURE.md` companion file discovery layers:
  removed Layer 1 (`_discover_companion_files`), renumbered
  Layers 2-6 to 1-5
- Updated intermediate file filtering description

### Files changed

- `agent/graph_nodes.py` — Removed `_discover_companion_files`
  function and its call from `perceive()`
- `tests/tst_v112_13_fixes.py` — Replaced companion file tests
  with perceive-pipeline tests
- `docs/reference/ARCHITECTURE.md` — Updated layer numbering
- `docs/project/CHANGELOG.md` — This entry

---

## Version 112.78 (GUI mode: map_coeffs_mtz empty after refine; daily usage Sorry; after_program premature stop)

### Problem 1 — map_coeffs_mtz empty after refine

After phenix.refine completes in GUI mode, `best_files["map_coeffs_mtz"]` stays
empty, causing the server to fail when building the ligandfit command
("missing required input: map_coeffs_mtz").

### Root cause

Two compounding bugs:

| # | Location | Bug |
|---|----------|-----|
| 1 | `_execute_command` | Uses `os.getcwd()` as `working_dir`, but the OldStyle runner restores CWD to the parent agent directory. Output files are in `sub_03_refine/` but `extract_output_files` and `scan_directory_for_output_files` look in the parent. Result: `output_files` passed to `record_result` is empty/incomplete, so `best_files` stays `{}`. |
| 2 | Safety net (line ~3233) | Condition `if best_files and not best_files.get("map_coeffs_mtz")` short-circuits when `best_files` is `{}` (empty dict = falsy in Python), skipping the safety net entirely. |

The files DO appear in `active_files` (via `get_available_files` step-3 catch-all
directory scan), but they never pass through `evaluate_file`, so the best_files
tracker has no entries.

### Fix

**Bug 1 — return output_dir from sub-job runner:**
`_execute_sub_job_for_gui` now returns a 4-tuple `(log_text, error_text,
executed_command, gui_output_dir)`.  `_execute_command` uses `gui_output_dir`
as `working_dir` instead of `os.getcwd()` when in GUI mode.  The agent's own
log file is still written to `os.getcwd()` (`log_dir`) — only file scanning
uses the sub-job directory.

**Bug 2 — safety net condition:**
Changed `if best_files and not best_files.get("map_coeffs_mtz")` to
`if best_files is not None and not best_files.get("map_coeffs_mtz")`.

**Bug 3 — `_track_output_files` same wrong directory:**
`_track_output_files` also used `os.getcwd()` for its directory scan.  Now
accepts `working_dir` parameter and uses it when provided.

**Diagnostic enhancement:** The MAP_COEFFS_MTZ diagnostic now prints the
`best_files` entry count.  `0 entries` → `record_result` never processed the
output files (working_dir bug); `>0 entries` → MTZ classification mismatch.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | `_execute_sub_job_for_gui`: return 4-tuple with `gui_output_dir` (3 return paths) |
| `programs/ai_agent.py` | `_execute_command`: unpack 4-tuple, separate `log_dir` from `working_dir` |
| `programs/ai_agent.py` | `_track_output_files`: accept and use `working_dir` param |
| `programs/ai_agent.py` | Safety net: `is not None` instead of truthiness check |
| `programs/ai_agent.py` | Diagnostic: show entry count, differentiate two failure modes |
| `tests/tst_autosol_bugs.py` | +4 tests (Bug 5), total 48 |

### Problem 2 — Daily usage limit not raised as Sorry

When the server returns `daily_usage_reached`, `rest/__init__.py` printed
"Skipping..." and returned `success=False`.  RemoteAgent treated this as
"no response" and returned None.  The cycle loop saw None as "no command
generated" and ended silently, producing a confusing summary that said
xtriage "failed."

Even after `rest/__init__.py` was fixed to raise Sorry, RemoteAgent's
generic `except Exception` around `_send_request()` caught it (Sorry
inherits from Exception), logged it as `"Error calling server: ..."`,
and returned None — swallowing the fatal error.

### Fix — three files

**`rest/__init__.py`:** Raise `Sorry` on `daily_usage_reached` at both
detection points (~lines 784 and 1133) instead of silently returning.
Consistent with how other auth failures are already handled.

**`phenix_ai/remote_agent.py`:** Add `except Sorry: raise` before the
generic `except Exception` handler around `_send_request()`.  This lets
Sorry propagate while still catching transient network errors.

**`programs/ai_agent.py`:** No changes needed — Sorry propagates
naturally up from `decide_next_step()` through the cycle loop.

### Additional files changed

| File | Change |
|------|--------|
| `rest/__init__.py` | Both `daily_usage_reached` paths: raise Sorry (was silent return/break) |
| `phenix_ai/remote_agent.py` | `except Sorry: raise` before generic handler + import Sorry |
| `tests/tst_autosol_bugs.py` | +2 tests (Bug 6), total 50 |

### Problem 3 — `after_program` premature stop on multi-goal requests

When user advice contains multiple goals (e.g., *"improve a map, get symmetry
and map correlation"*), the directive extractor sets `after_program` to a
single program (whichever it recognizes last).  `check_directive_stop()` in
**PERCEIVE** fired immediately when that program completed — before the LLM
ran — short-circuiting the entire cycle and preventing the remaining goals.

Example: mtriage (correlation) ✓ → map_symmetry (symmetry) ✓ → **STOP** ✗
(map improvement never ran).

### Fix

Removed `after_program` hard stop from `check_directive_stop()` in
`agent/perceive_checks.py` (12 lines).  `after_program` is now only used in
the **PLAN** node as a minimum-run guarantee: *"don't auto-stop until at
least this program has run."*  The LLM always runs and decides whether
all goals are met.

`after_cycle` and metrics targets (`r_free_target`, `map_cc_target`) remain
as hard stops — those are explicit numeric limits, not inferred from advice.

**Trade-off:** Single-goal focused tasks ("evaluate these data") now use one
extra LLM call (the LLM runs on the next cycle and decides to stop itself).

### Files changed

| File | Change |
|------|--------|
| `agent/perceive_checks.py` | Removed `after_program` hard stop block + dead `last_command` variable, added explanatory comment |
| `agent/session.py` | Removed parallel `after_program` hard stop from `check_directive_stop_conditions` fallback |
| `agent/directive_extractor.py` | Removed `after_program` hard stop from `check_stop_conditions` (server-side) |
| `agent/graph_nodes.py` | Improved PLAN log message when `after_program` completed + auto-stop |
| `tests/tst_autosol_bugs.py` | +7 tests (Bug 7): functional + regression guards + all 3 parallel implementations + dead code, total 57 |

### Windows compatibility (v112.78)

Four cross-platform issues identified and fixed:

| Issue | File | Fix |
|-------|------|-----|
| Temp-dir filter missed backslash paths | `agent/graph_nodes.py` | Normalize `\\` → `/` before matching `/TEMP`, `/TEMP0/`, `/scratch/` markers |
| Console window flash in GUI | `programs/ai_agent.py` | `CREATE_NO_WINDOW` creationflag on `Popen` when `os.name == 'nt'` |
| Abort detection assumed negative return codes | `programs/ai_agent.py` | Documented that `return_code < 0` is Unix-only; STOPWIZARD is the cross-platform indicator |
| Session JSON used locale encoding | `agent/session.py` | All 3 `open()` calls now specify `encoding='utf-8'` |

+4 Windows tests, total 61.

## Version 112.77 (Autobuild rebuild_in_place stripped by Rule D)

### Problem

The LLM correctly identifies `rebuild_in_place=False` as the fix for a sequence
mismatch error in autobuild, but Rule D in `sanitize_command` strips it.  Rule D
removes bare key=value params not in the program's strategy_flags allowlist.
Autobuild had only 3 flags (quick, nproc, resolution).

### Fix

Added `rebuild_in_place`, `n_cycle_build_max`, and `maps_only` to autobuild's
strategy_flags in programs.yaml.  Added hint about `rebuild_in_place=False`
for sequence mismatch recovery.

### Files changed

| File | Changes |
|------|---------|
| `knowledge/programs.yaml` | +15 lines: 3 new strategy_flags for autobuild + recovery hint |
| `tests/tst_autosol_bugs.py` | +70 lines: 7 new tests (Bug 4), total 44 |

### Broader concern — Rule D is too aggressive

Most programs have only 1–4 strategy_flags but accept dozens of PHIL params.
Rule D strips anything not pre-enumerated, meaning the agent can only use params
that were anticipated at configuration time.  This is safe against hallucination
but blocks legitimate error recovery when the LLM discovers the correct fix.

The v112.76 catch-all blacklist handles the reverse case (bad params that *cause*
errors).  But Rule D has no feedback loop: a stripped correct param never generates
an error, so nothing learns that stripping was wrong.

A future "warn but keep" mode for Rule D could log a warning but preserve unknown
bare params, relying on PHIL validation + catch-all blacklist (N=2 threshold) to
handle bad ones.  Trade-off: at most 2 wasted cycles for a hallucinated param vs
potentially infinite wasted cycles when the agent keeps proposing a correct param
that keeps getting stripped.

Programs most affected: phenix.phaser (2 flags), phenix.autobuild_denmod (2),
phenix.ligandfit (1), phenix.dock_in_map (1).

## Version 112.76 (Catch-all injection blacklist + deterministic atom_type)

Follow-up improvements from the v112.75 autosol bugs analysis.

### Deterministic atom_type — heavier atom wins

**Problem:** LLM directive extraction non-deterministically assigns which
element goes to `atom_type` vs `additional_atom_types`.  Run 109 picked
`atom_type=S` instead of `Se` — the weaker scatterer at 0.9792 Å.

**Design rationale:** Rather than trying to fix the LLM's non-determinism,
apply a deterministic physical rule: the element with higher atomic number (Z)
always provides stronger anomalous signal at typical synchrotron wavelengths.
This is applied as a post-validation swap, not a pre-extraction constraint,
so it catches errors from both LLM and simple-fallback extraction paths.

**Fix:** After the v112.75 dedup check, a new validation compares atomic
numbers (Z) of `atom_type` and `mad_ha_add_list`.  If the primary scatterer
has lower Z, the values are swapped.  This is scientifically correct: heavier
atoms provide stronger anomalous signal in SAD/MAD experiments.

- `_ANOMALOUS_Z` table: 27 common anomalous scatterers (Big 5: S, Se, Zn,
  Fe, Hg, plus Pt, Au, and others from Z=16 to Z=92)
- Handles multi-element `mad_ha_add_list` (e.g., `Se+Zn`)
- Do-no-harm: skips swap when either element is unknown

### Catch-all injection blacklist — supervisor pattern

**Problem:** `bad_inject_params` learning only fires on recognized PHIL error
patterns.  Any unrecognized format causes the agent to loop indefinitely
(9 wasted cycles in run 107).

**Design rationale:** Pattern-based error recognition is inherently fragile
against PHIL's open-ended error space.  Instead of adding patterns one by one,
implement a supervisor that tracks what was injected and correlates it with
failures — no pattern matching required.  The threshold N=2 balances between
transient failures (filesystem/network rarely produce identical fingerprints
twice) and wasted cycles (N=3 would waste an extra cycle on what is almost
certainly deterministic).  The "innocent bystander" trade-off is accepted:
harmless params like `nproc=8` may be blacklisted alongside the real offender,
but breaking the loop outweighs losing a harmless optimization.

**Fix:** Track consecutive same-program failures via error fingerprint.  After
N=2 failures with the same fingerprint, blacklist whatever `inject_user_params`
appended — no pattern matching required.

- `postprocess_command`: new `return_injected=False` kwarg surfaces the list of
  params added by inject_* steps.  Default False for backward compatibility
  with 10+ existing call sites
- `_update_inject_fail_streak()`: streak tracker on ai_agent.py.  Computes
  error fingerprint (first 120 chars, normalized), increments on match, resets
  on mismatch or success, blacklists at threshold
- Recovery retries (`force_retry_program`) excluded — they have their own loop
  guard (v112.74)
- "Innocent bystander" trade-off accepted: harmless params like `nproc=8` may
  be blacklisted alongside the real offender, but breaking the loop outweighs
  losing a harmless optimization

### Files changed

| File | Changes |
|------|---------|
| `agent/command_postprocessor.py` | +120 lines: `_ANOMALOUS_Z` table, `_ensure_primary_scatterer_is_heavier()`, step 2c call, `return_injected` kwarg with injected-list collection |
| `agent/graph_nodes.py` | BUILD + FALLBACK nodes: call `postprocess_command(return_injected=True)`, store `last_injected_params` in graph state |
| `agent/session.py` | +60 lines: `_fix_autosol_atom_type_order()` validates atom_type at directive extraction time (both LLM and simple-fallback paths) |
| `programs/ai_agent.py` | +90 lines: `_update_inject_fail_streak()`, streak call in `_record_command_result`, `last_injected_params` init/extraction in `_get_command_for_cycle` and `_query_agent_for_command` |
| `tests/tst_autosol_bugs.py` | +270 lines: 23 new tests (10 Phase 1 + 9 Phase 2 + 3 directive-level + 1 mock drift), total 37 |

## Version 112.75 (Autosol/autobuild process bugs)

Three process bugs diagnosed from runs 107 and 109, all involving autosol/autobuild
SAD phasing workflows.

### Root cause summary

All three bugs stemmed from gaps in the injection/validation pipeline:

| Bug | Wasted cycles | Root cause category |
|-----|---------------|---------------------|
| 1. wavelength duplication | 9 | Alias-blind duplicate check + unlearnable error pattern |
| 2. atom_type swap | 1 | LLM non-determinism with no post-validation |
| 3. autosol re-run | 1 | `_is_program_already_done` only checked `run_once` strategy |

The `wavelength → autosol.lambda` mapping is uniquely treacherous: it is the only
strategy_flag where the flag leaf (`lambda`) differs completely from the
natural-language key (`wavelength`).  All other mappings have the key as a
substring of the flag.  Bugs 2 and 3 compounded in run 109 — wrong atom type
(S instead of Se) plus unnecessary re-run after both autosol and autobuild had
already succeeded.

### Problem 1 — Duplicate `wavelength=` crashes PHIL (9 wasted cycles)

**Symptom:** Every autosol attempt in run 107 crashes with "One True or False
value expected, autosol.wavelength.added_wavelength=0.9792 found". Repeats 9
times with no self-correction.

**Root cause:** `inject_user_params` extracts `wavelength=0.9792` from user
advice ("wavelength is 0.9792") and appends it to the command even though
`autosol.lambda=0.9792` is already present.  The duplicate check searches for
"wavelength" in the command string but finds only "autosol.lambda" — the
strategy_flags mapping (`wavelength → autosol.lambda={value}`) creates an alias
that the check doesn't know about.  PHIL resolves bare `wavelength=` to
`autosol.wavelength.added_wavelength` (a boolean) → type error.  The error
message doesn't match the `bad_inject_params` learning patterns ("unknown
parameter" / "no such parameter"), so the system never learns to stop.

**Fixes:**
- `command_postprocessor.py`: Build alias map from strategy_flags (`wavelength →
  lambda`) and check alias leaves in duplicate detection before injecting
- `command_postprocessor.py`: Autosol-specific dedup — strip `mad_ha_add_list`
  when identical to `atom_type` (prevents secondary scatterer loss)
- `ai_agent.py`: Add "True or False value expected" to `bad_inject_params`
  error learning patterns, extracting and blacklisting the offending PHIL path
  components

### Problem 2 — `atom_type=S` instead of `Se` (wrong primary scatterer)

**Symptom:** Run 109 extracts `atom_type=S` from "use Se and S as anomalous
atoms" — selenium (the stronger scatterer at 0.9792 Å) is entirely lost from
cycle 2's command (`atom_type=S mad_ha_add_list=S`, both sulfur).

**Root cause:** LLM-based directive extraction is non-deterministic.  Run 107
got `atom_type=Se` (correct); run 109 got `atom_type=S` (wrong).  No validation
catches the case where both parameters are set to the same element.

**Fixes:**
- `command_postprocessor.py`: Post-assembly validation — if `atom_type` and
  `mad_ha_add_list` have the same value, strip the duplicate `mad_ha_add_list`
- `programs.yaml`: Improved hint text explicitly states the heavier element
  (higher Z) should be `atom_type`, and that it must differ from
  `additional_atom_types`

### Problem 3 — Autosol re-runs after autosol+autobuild both succeeded

**Symptom:** Run 109 cycle 4 runs autosol again after both autosol (cycle 2)
and autobuild (cycle 3) succeeded.  Should have proceeded to refinement.

**Root cause:** `_apply_directives` re-adds `phenix.autosol` from
`program_settings` in directives (extracted from user advice mentioning autosol
parameters).  It calls `_is_program_already_done` which only checked `run_once`
programs — autosol uses `set_flag` strategy, so it returned False.  Autosol
gets inserted at position 0 in `valid_programs` and the LLM picks it.  The
belt-and-suspenders filter in `get_valid_programs` runs before
`_apply_directives`, so it can't catch the re-addition.

**Fix:**
- `workflow_engine.py`: Extend `_is_program_already_done` to check non-count
  programs with program-specific done flags (mirrors the existing
  belt-and-suspenders filter).  This prevents `_apply_directives` from
  re-adding any completed non-count program.

### Files changed

| File | Changes |
|------|---------|
| `agent/command_postprocessor.py` | +37 lines: strategy-flag alias map, alias duplicate check, autosol atom_type dedup |
| `agent/workflow_engine.py` | +21 lines: `_is_program_already_done` extended for non-count programs |
| `programs/ai_agent.py` | +35 lines: "True or False" error pattern in `bad_inject_params` learning (both branches) |
| `knowledge/programs.yaml` | Improved `atom_type` and `additional_atom_types` hints |
| `tests/tst_autosol_bugs.py` | +240 lines: 14 regression tests covering all fixes |

## Version 112.74 (Xtriage recovery + ligand misclassification)

### Problem 1 — Xtriage obs_labels recovery never reaches command line

**Symptom:** Xtriage fails with "Multiple equally suitable arrays", recovery
correctly selects I(+) labels, but the retry command is identical to the
failed command — obs_labels is missing.  Agent loops until fallback STOPs.

**Root cause (4-layer failure chain):**
1. `registry.build_command()` silently drops recovery params that don't match
   `strategy_flags` keys in programs.yaml.  Xtriage only has `unit_cell` and
   `space_group` — the fully-qualified `scaling.input.xray_data.obs_labels`
   is ignored.
2. Probe-only sanitizer strips obs_labels (not a file path).
3. Duplicate check consumes `force_retry_program` on identical command.
4. No guard against re-triggering the same failed recovery.

**Fixes:**
- `command_builder.py`: After `_assemble_command`, append any recovery-sourced
  strategy entry missing from the command.
- `command_postprocessor.py`: Whitelist data-label params (`obs_labels`,
  `labels`, `data_labels`, `anomalous_labels`, `r_free_flags_labels`) in
  probe-only sanitizer.
- `ai_agent.py`: Skip duplicate check when `forced_retry` is set.  Guard
  against re-triggering recovery when strategy already exists for file.
- `graph_nodes.py`: Propagate `forced_retry` to top-level BUILD output.

### Problem 2 — Protein PDB files misclassified as ligand

**Symptom:** `1aba.pdb` (729 ATOM + 20 HETATM records) rejected by BUILD:
"LLM file rejected (in excluded category 'ligand'): model=1aba.pdb"

**Root cause:** YAML categorizer patterns for `ligand_pdb` are broad enough
to match protein PDB filenames.  Existing post-processing only validates
`unclassified_pdb` files — files directly placed into `ligand_pdb` were
never content-checked.

**Fix:** `workflow_state.py`: Added post-processing guard in
`_categorize_files()`.  Iterates `ligand_pdb` entries; any PDB file where
`_pdb_is_protein_model()` returns True (>150 atoms, majority ATOM records)
is moved to `unclassified_pdb`/`pdb`/`model`.  Same pattern as existing
half-map validation guard.

### Files modified

| File | Changes |
|---|---|
| `command_builder.py` | Append recovery-sourced strategy entries after build_command |
| `command_postprocessor.py` | Label param exception in probe-only sanitizer |
| `ai_agent.py` | Skip duplicate check for recovery; recovery loop guard |
| `graph_nodes.py` | Propagate forced_retry to top-level BUILD output |
| `workflow_state.py` | Ligand_pdb content validation guard |

---

## Version 112.73 (Output files and best_files lost on restart)

### Problem 1 — Available files empty on restart

After refinement completes and the agent stops then restarts, ligandfit fails
with "missing inputs: map_coeffs_mtz" even though refine output files exist
on disk.  The log shows only 3 available files — refine outputs are absent.

Root cause: `get_available_files()` Step 3 infers `agent_dirs` from tracked
file paths.  With only original inputs tracked (no `sub_*` in paths),
`agent_dirs` is empty and the supplemental scan never runs.

Fix: Seed `agent_dirs` from `_get_session_dir()` — always known from the
session file path.

### Problem 2 — best_files stale after restart (wrong model for refine)

After ligandfit and pdbtools succeed, refine uses the original input model
instead of the `*_with_ligand.pdb` from pdbtools.  Also uses original data
MTZ instead of the R-free-locked MTZ from the first refinement.

Root cause: `_rebuild_best_files_from_cycles()` uses `_discover_cycle_outputs()`
which constructs directory names as `sub_{cycle_num}_{program}*`.  After a
restart, session cycle numbers diverge from GUI directory numbers (e.g.,
session says cycle 2 but GUI created `sub_04_pdbtools`).  The scan finds
nothing, `best_files` stays at originals, and PRIORITY 2 in the command
builder returns the stale model.

Fix: `_discover_cycle_outputs()` now has five strategies:

1. Try stored `output_files` paths as-is
2. Resolve relative paths against session directory
3. Use stored `output_dir` from cycle data (new: set by `record_result`)
4. Scan `sub_{NN}_{program}*` (exact cycle-number match)
5. Scan `sub_*_{program}*` (broad match — handles number mismatch)

### Problem 3 — Directory scan finds files but best_files ignores them

Even after fixes 1-2, polder (and refine) still selects the ligand-free model.
The with_ligand PDBs ARE in `available_files` (13 files), but `best_files`
still shows `model=refine_001_001.pdb`.

Root cause: Architectural gap — `get_available_files()` Step 3 (directory scan)
adds files to the available list, but `_rebuild_best_files_from_cycles()` only
evaluates files from known cycles.  Files discovered only by Step 3 are
invisible to `best_files`.

Fix: New Step 3.5 in `get_available_files()` — after the directory scan,
evaluate all newly discovered files through `best_files.evaluate_file()`.
Uses `highest_cycle` from session data for recency.  Safe because the `seen`
set ensures no double evaluation of files already processed by Steps 1-2.

### Companion fixes (defense-in-depth)

- `_import_mtz_utils()` always returns working functions (workflow_state.py)
- `_should_exclude()` for dual-categorization (command_builder.py)

### Files to deploy

| File | What changed | Role |
|------|-------------|------|
| `agent/session.py` | All 3 fixes above | **Critical** |
| `agent/workflow_state.py` | Self-contained `_import_mtz_utils()` | Defense-in-depth |
| `agent/command_builder.py` | `_should_exclude()` at 4 check points | Defense-in-depth |
| `agent/file_utils.py` | Contains `get_mtz_stage` | Defense-in-depth |

---

## Version 112.72 (Invalid space_group directive injection)

The LLM directive extractor sometimes misinterprets workflow descriptions as
parameter values. User advice mentioning "space group determination" caused
`space_group=determination` to be extracted into `program_settings.default`,
then injected verbatim into `phenix.refine` commands where it caused errors.

### Root cause

The LLM parsed "space group determination" as `space_group = determination`,
treating the workflow description word as a space group symbol. No validation
existed to distinguish real space group symbols (P1, P 21 21 21, C2221) from
English words.

### Fix: Four-layer defense

**Layer 1 — Don't extract it** (`session.py`): New `_sanitize_directives()`
method runs after directive extraction. Validates `space_group` values in all
`program_settings` scopes using `_is_valid_space_group()`. Invalid values are
removed before storage. Logs removal.

**Layer 2 — Don't inject via strategy** (`command_builder.py`): `_build_strategy()`
now validates any `space_group` key in the strategy dict before it reaches
`_assemble_command()`. Invalid values are stripped with a log message.

**Layer 3 — Don't inject via crystal symmetry** (`command_postprocessor.py`):
`inject_crystal_symmetry()` now uses `_is_valid_space_group()` instead of the
weak `_INVALID_SG_PATTERNS` substring check.

**Layer 4 — Strip if already in command** (`command_postprocessor.py`):
`sanitize_command()` Rule B2 validates `space_group=` values in the command
string and strips invalid ones.

### Space group validation (`_is_valid_space_group`)

New validator in `command_postprocessor.py` checks:
- Known placeholder phrases (extended: "determination", "detection", "analysis", "auto", etc.)
- Length limit (>20 chars rejected)
- Space group numbers (1-230 accepted)
- First character must be a crystal system letter (P/I/C/F/R/A/B/H) or digit
- Single English words >6 chars rejected (catches "determination", "monoclinic", etc.)
- Short symbols like "Pnma", "Pbca" correctly accepted

### Files Changed
- `agent/command_postprocessor.py` — New `_is_valid_space_group()`, updated
  `inject_crystal_symmetry()`, new Rule B2 in `sanitize_command()`
- `agent/command_builder.py` — Space group validation in `_build_strategy()`
- `agent/session.py` — New `_sanitize_directives()` method, called after extraction

---

## Version 112.71 (MTZ categorization safety net + dual-categorization fix)

After running `phenix.refine`, the agent could not find `refine_001_001.mtz` as
map coefficients for `phenix.ligandfit`, stopping with `missing inputs:
map_coeffs_mtz`. Three independent categorization gaps caused the failure.

### Bug 1 — Hardcoded categorizer skipped MTZ subcategories

**Root cause:** `_categorize_files_hardcoded()` called `classify_mtz_type()` to
place `refine_001_001.mtz` in the parent `map_coeffs_mtz` category but never
populated the `refine_map_coeffs` subcategory. When the command builder searched
subcategories first via `priority_categories` (denmod → predict_build → refine →
parent), it found nothing in the subcategories and sometimes missed the parent.

**Fix:** After `classify_mtz_type()` returns `map_coeffs_mtz`, now also calls
`get_mtz_stage()` to populate the correct subcategory (`refine_map_coeffs`,
`denmod_map_coeffs`, or `predict_build_map_coeffs`).

### Bug 2 — Dual-categorization caused exclude_categories rejection

**Root cause:** The YAML categorizer's two-step process could place a file in
BOTH `data_mtz` (Step 1: extension match when exclude patterns didn't fire) AND
`map_coeffs_mtz` (Step 2: pattern match → subcategory → bubble-up to parent).
The command builder's `exclude_categories: [data_mtz]` for ligandfit then
rejected the file even though it was correctly in `map_coeffs_mtz`.

**Fix:** Post-categorization safety net now detects files in both `data_mtz` and
`map_coeffs_mtz` and removes them from `data_mtz`.

### Bug 3 — No post-categorization cross-check

**Root cause:** There was no safety net after categorization to catch cases where
YAML patterns and the authoritative regex classifier disagreed. If the server's
`file_categories.yaml` had stale patterns, files could be silently misclassified.

**Fix:** Added a post-categorization safety net at the end of `_categorize_files()`
that cross-checks ALL MTZ files against the authoritative `classify_mtz_type()`
regex. If any file is in the wrong category, it is moved to the correct one with
appropriate subcategory assignment. Logs `WARNING` when corrections are made.

### Diagnostic logging improvements

Added two levels of diagnostic logging to make future occurrences immediately
identifiable:

1. **perceive() node** (`graph_nodes.py`): Logs all MTZ category contents after
   categorization. Emits a WARNING when refinement is in history but
   `map_coeffs_mtz` is empty.

2. **Safety net** (`workflow_state.py`): Logs `WARNING` via Python logger when it
   corrects a misclassification, e.g.:
   `MTZ safety net: moved refine_001_001.mtz from data_mtz to refine_map_coeffs/map_coeffs_mtz`

### Files Changed
- `agent/workflow_state.py` — Added `import logging`, MTZ safety net post-processing
  in `_categorize_files()`, subcategory population in `_categorize_files_hardcoded()`
- `agent/graph_nodes.py` — Added MTZ category diagnostic logging in `perceive()`

---

## Version 112.70 (Ligandfit missing data file + model/ligand swap)

After running `phenix.refine`, the agent tried to run `phenix.ligandfit` with
no data file, swapped model and ligand, and used `atp.pdb` (the small molecule)
as the protein model. Three independent bugs:

### Bug 1 — `refine_001.mtz` classified as `data_mtz` instead of `map_coeffs_mtz`

**Root cause:** `session._rebuild_best_files_from_cycles` had a hardcoded regex
`refine_\d+_001\.mtz$` that only matched two-level serial output
(`refine_001_001.mtz`), NOT the standard single-serial output (`refine_001.mtz`).
So `best_files["map_coeffs_mtz"]` was never populated after refinement, and
ligandfit's `require_best_files_only: true` on `map_coeffs_mtz` meant neither
BUILD nor the safety net could find the data file.

**Fix:** Updated regex to `(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$` — matches
all four patterns: `refine_001.mtz`, `refine_001_001.mtz`, `7qz0_refine_001.mtz`,
`7qz0_refine_001_001.mtz`. Applied in three locations:
- `session._rebuild_best_files_from_cycles` (session load path)
- `session.record_result` (live cycle recording path)
- `file_utils.classify_mtz_type` (shared classifier used by BestFilesTracker)

**Secondary fix:** Added `refine_map_coeffs`, `denmod_map_coeffs`,
`predict_build_map_coeffs`, `original_data_mtz`, and `phased_data_mtz` to
`BestFilesTracker.STAGE_TO_PARENT`. Without these, passing `stage="refine_map_coeffs"`
to `evaluate_file()` caused the category to be inferred from the filename instead
of the stage, creating a silent mismatch.

### Bug 2 — `exclude_patterns` used substring matching

The model slot's `exclude_patterns: [ligand, ...]` used substring matching
(`pat in basename`), so `nsf-d2_noligand.pdb` was incorrectly excluded from the
model slot because "noligand" contains "ligand". This caused the safety net to
skip the protein model and pick `atp.pdb` instead.

**Fix:** New `matches_exclude_pattern()` function in `agent/file_utils.py` uses
word-boundary matching: patterns must appear at the start of the name or after a
separator (`_`, `-`, `.`). Applied in both `CommandBuilder._find_file_for_slot`
(server) and `_inject_missing_required_files._find_candidate_for_slot` (client).
Also applied to `prefer_patterns` for consistency.

### Bug 3 — No content-based small-molecule guard for model slots

Even after fixing exclude_patterns, `atp.pdb` (HETATM-only small molecule) could
still be selected as a "model" because its filename has no ligand-like pattern.

**Fix:** Both `CommandBuilder._find_file_for_slot` and
`_inject_missing_required_files._find_candidate_for_slot` now call
`_pdb_is_small_molecule()` (reads first 8KB, checks for HETATM-only) when filling
model/protein/pdb_file slots. Small-molecule PDB files are rejected from model
slots and left for the ligand slot.

### Bug 4 — Protein model PDB assigned to ligand slot

The LLM assigned `refine_001_001.pdb` (the protein model) to both the `model` and
`ligand` slots of `phenix.ligandfit`, causing it to try fitting the protein as a
ligand instead of the actual ATP file.

**Fix:** Added inverse content-based guard using new `_pdb_is_protein_model()`
function. When filling the `ligand` slot, PDB files positively identified as
protein models (contain ATOM records) are rejected. Applied in three locations:
- `CommandBuilder` LLM file hint validation (rejects LLM's wrong assignment)
- `CommandBuilder._find_file_for_slot` (auto-fill path)
- `_inject_missing_required_files._find_candidate_for_slot` (safety net)

Uses `_pdb_is_protein_model()` (positive protein check) rather than
`not _pdb_is_small_molecule()` because the latter returns False for non-existent
files, which would incorrectly reject valid candidates.

### Bug 5 — `inject_user_params` re-injects stripped bare params

`sanitize_command` (Rule D) correctly stripped hallucinated bare params like
`d_min=2.5` and `elements='Se S'` from autosol/autobuild commands, but
`inject_user_params` re-added them from user advice text.

**Fix:** `inject_user_params` now validates bare (undotted) keys against the
program's strategy_flags allowlist before injection, mirroring Rule D's check.

### Bug 6 — Refinement restraints CIF used as ligand

After the protein PDB guard correctly rejected `refine_001_001.pdb` from the
ligand slot, auto-fill grabbed `refine_001_001.cif` (refinement geometry
restraints) as the ligand instead of the actual ATP file.

**Fix:** Three changes:
1. Added `refine` to `exclude_patterns` for ligandfit ligand slot in programs.yaml.
   Word-boundary matching catches `refine_001_001.cif`, `refine_001.cif`,
   `7qz0_refine_001.cif` but NOT `atp_refined.cif`.
2. LLM-selected files now checked against slot `exclude_patterns` before acceptance
   (previously only auto-fill applied exclude_patterns).
3. New `_pdb_is_protein_model()` in `workflow_state.py` provides positive protein
   detection for ligand slot guard (safe for non-existent files).

### Bug 7 — Agent stops with misleading "all_commands_duplicate" on resume

When resuming after refinement, the agent stopped with `Workflow complete:
all_commands_duplicate` instead of running ligandfit. Two independent issues:

**Root cause A:** `_rebuild_best_files_from_cycles` only evaluated files listed in
the cycle's `output_files`. If the client didn't track the map coefficients MTZ
(e.g. `refine_001.mtz`), then `best_files["map_coeffs_mtz"]` was never populated.
Ligandfit's `require_best_files_only: true` then caused its build to fail.

**Fix A:** Both `_rebuild_best_files_from_cycles` (session load) and `record_result`
(live path) now call `_find_missing_outputs` and evaluate supplemental files through
the best_files tracker. On the live path, discovered files are also appended to the
cycle's `output_files` so `get_available_files` picks them up immediately.
All three file-discovery paths now share the same supplemental logic:
- `get_available_files` → already called `_find_missing_outputs`
- `_rebuild_best_files_from_cycles` → now calls it (session load)
- `record_result` → now calls it (live cycle completion)

**Root cause B:** Fallback reported `all_commands_duplicate` whether the failure was
duplication or inability to build (missing inputs). Misleading stop reason.

**Fix B:** Fallback now distinguishes `cannot_build_any_program`, 
`build_failures_and_duplicates`, and `all_commands_duplicate`. The `abort_message`
includes per-program diagnostics showing which input slots couldn't be filled.
`CommandBuilder.build()` now stores `_last_missing_slots` for diagnostics.

### Bug 8 — Iterative refinement flagged as duplicate command

When the user asked for additional refinement cycles, the duplicate detector
flagged the new `phenix.refine` command as "similar to cycle 1" even though the
input model file changed (e.g. `refine_002.pdb` vs `model.pdb`). The 80%
token-overlap heuristic was comparing basenames of all command tokens, so commands
with many shared params (nproc, macro_cycles, same data file) exceeded the 80%
threshold even when the model file was different.

**Fix:** Before applying the overlap heuristic, extract file tokens (basenames with
crystallographic extensions like .pdb, .mtz, .cif, .sca) from both commands. If the
file token sets differ, the commands are NOT duplicates regardless of overall overlap.
Running `phenix.refine` with a different model is a fundamentally different
computation.

### Bug 9 — `best_files` skipped for slots with specific subcategories

Even after Bugs 1+7 correctly populated `best_files["map_coeffs_mtz"]` from
supplemental file discovery, ligandfit's `map_coeffs_mtz` slot still failed to
build. The `_find_file_for_slot` method skips `best_files` entirely when a slot's
`input_priorities.categories` contains a "specific subcategory" (like
`refine_map_coeffs`). When the category-based lookup also failed (because
`categorized_files` didn't have the file in `refine_map_coeffs`), the extension
fallback was blocked by `require_best_files_only`, and the slot returned None.

**Fix:** Added a `best_files_fallback` path: when both category-based lookup and
extension fallback are exhausted for a specific-subcategory slot, try `best_files`
as a last resort with `exclude_categories` validation. This ensures the
supplemental-discovered map coefficients MTZ is used even when the PERCEIVE node's
file categorizer didn't place it in the specific subcategory.

### Bug 10 — Advice display at verbose level

User advice (`project_advice`) was displayed at `verbose` level in
`_print_iterate_agent_header`, invisible at the default `normal` verbosity.

**Fix:** Changed to `normal` level so users always see their advice in the log.

### Bug 11 — `set_project_info` replaces `original_files` on resume

When resuming a session, `set_project_info(original_files=...)` **replaced** the
existing `original_files` list with whatever files the user supplied on this run.
If the user initially provided `data.mtz model.pdb atp.pdb` but resumed with only
`refine_001_data.mtz refine_001_001.pdb`, the ligand file (`atp.pdb`) was lost.
This caused ligandfit to fail because no ligand file was available.

**Fix:** `set_project_info` now **merges** new files into the existing list on
resume (deduplicating by basename). Original files from the first run are preserved
even if the user doesn't re-supply them.

### Bug 12 — BUILD failure message lacks missing-slot detail

When `CommandBuilder.build()` returned None (missing required inputs), the BUILD
node's `validation_error` was just `"Failed to build command"`. This became the
fallback reasoning shown to the user: `"Fallback: phenix.ligandfit could not be
built (Failed to build command; Failed to build command)"` — no indication of
*which* files were missing.

**Fix:** BUILD node now checks `builder._last_missing_slots` and includes slot
names in the error message: `"Failed to build command (missing: ligand,
map_coeffs_mtz)"`. The fallback reasoning also includes its own `build_failures`
diagnostics (from Bug 7's per-program tracking) when they differ from the BUILD
node errors.

### Bug 13 — Protein guard too aggressive for ligand slot

`_pdb_is_protein_model` returned True if ANY `ATOM` record existed, which
rejected legitimate ligand files (e.g., `atp.pdb` — an ATP nucleotide file
that may include a few ATOM records from extracted protein context). This was
the **root cause** of the `missing inputs: ligand` failure: `atp.pdb` was
present in available files but got excluded by the over-aggressive guard.

**Fix:** `_pdb_is_protein_model` now uses size-based detection:
- Counts total coordinate records (ATOM + HETATM) in first 32 KB
- Small files (≤ 150 coordinate records) → **not protein**, regardless of
  record type (ligands like ATP only have ~31 atoms; smallest crystallographic
  protein has ~500+ atoms)
- Larger files with majority ATOM → protein model (rejected from ligand slot)
- Increased read buffer from 8 KB to 32 KB (~400 lines) for reliable counting
- Also added detailed logging ("Ligand guard passed/excluded") for diagnostics

Additional diagnostic improvements:
- Fallback node logs available .pdb/.cif files and best_files for debugging
- Fallback build logs per program are captured and forwarded to debug_log
- FALLBACK debug logs displayed alongside BUILD/PLAN/PERCEIVE in client

### Bug 14 — `wavelength=` instead of `autosol.lambda=` for autosol

The LLM generated `wavelength=0.9792` (sometimes alongside the correct
`autosol.lambda=0.9792`). Autosol requires `autosol.lambda=`, not bare
`wavelength=`.

**Fix:** Added `wavelength → autosol.lambda` rename in `parameter_fixes.json`.
Also fixed `fix_program_parameters` logic: when the correct parameter already
exists in the command, the wrong parameter is now **removed** (previously it
was just skipped, leaving both in the command).

### Bug 15 — Working directory shows server path in summary

The final session summary displayed the server's temporary working directory
(e.g., `/net/cci-gpu-01/.../rest/c65be5a5-...`) instead of the client's
actual working directory. Both `generate_log_for_summary` and
`_get_working_directory_path` fell back to `os.getcwd()` which is wrong
on the server.

**Fix:** The client now stores `client_working_directory` in session data at
startup. `_get_working_directory_path()` checks this first (before
`input_directory`, session file inference, or `os.getcwd()`).
`generate_log_for_summary()` now calls `_get_working_directory_path()` instead
of raw `os.getcwd()`. `_get_run_name()` also uses `client_working_directory`
as fallback.

### Bug 16 — Advice preprocessor leaks instructions into LLM output

The `advice_preprocessor.py` prompt template contained section headers with
system instructions like `(CRITICAL — translate natural language to PHIL
key=value pairs)`. The LLM echoed these back verbatim in its output, making
them appear in the summary analysis shown to users.

**Fix:** Removed instruction language from section headers. The PHIL translation
guidance remains in the body text but is no longer in a position where the LLM
treats it as content to reproduce.

### Bug 17 — Agent stops prematurely when user provides already-refined model

When the user provides a pre-refined model (e.g., `7qz0.pdb` with R-free=0.204)
and asks for ligand fitting, the agent runs `xtriage` → `model_vs_data` and then
stops, even though ligand fitting was explicitly requested.

**Root causes:**
1. `model_vs_data` as a placement probe set `validation_done=True` (from YAML
   done_tracking), tricking the workflow into thinking validation was complete
2. The model was "at target" (good R-free), so `_is_at_target` returned True,
   causing refinement programs to be removed and STOP to be added
3. The ligandfit YAML condition `refine_count > 0` failed because no refinement
   was done in this session (the model came pre-refined)
4. `_apply_directives` `prefer_programs` only reorders existing valid programs
   — it couldn't add ligandfit since it wasn't in the list

**Fixes:**
- **`_is_at_target` returns False when refinement is a ligandfit prerequisite**:
  When `user_wants_ligandfit=True`, `ligandfit_done=False`, and `refine_count==0`,
  `_is_at_target()` returns False even if R-free is below target. This is the
  primary mechanism: it keeps `phenix.refine` in valid_programs so the YAML
  workflow naturally flows: refine (produces map_coeffs) → ligandfit (YAML
  conditions pass once refine_count > 0).
- **`_apply_directives` defers ligandfit injection until refine runs**: Previously,
  `_apply_directives` injected `phenix.ligandfit` into valid_programs even when
  `refine_count==0`, bypassing YAML conditions. This caused the LLM to pick
  ligandfit prematurely (which would fail in BUILD due to missing map_coeffs_mtz).
  Now, injection is deferred until `refine_count > 0`, letting the YAML conditions
  gate ligandfit correctly.
- **`_get_valid_programs_for_phase` safety net**: Defense-in-depth check that keeps
  refine when ligandfit needs it, in case `_is_at_target` logic changes.
- `_analyze_history`: When `model_vs_data` ran as a placement probe (before any
  refinement), unset `validation_done` to prevent premature workflow completion
- **`require_best_files_only` implemented**: Ligandfit's `map_coeffs_mtz` slot
  has `require_best_files_only: true` in YAML to prevent the raw Fobs MTZ from
  being used. This flag was defined but not honored. Now, slots with this flag
  skip extension-based fallback and only accept files from `best_files`.

### Bug 18 — String metric crash in YAML condition evaluation

Metric values from log analysis can arrive as strings (e.g., `"0.204"` instead of
`0.204`). This caused `'<' not supported between instances of 'str' and 'float'`
in `_check_metric_condition` when evaluating YAML workflow conditions.

**Fixes:**
- `_get_metric()` now casts return values to `float()` — this is the source of all
  metric values in the context, so all downstream comparisons are safe
- `_check_metric_condition()` also casts the value (defense-in-depth for metrics
  that bypass `_get_metric`)

### Bug 19 — `atp.pdb` miscategorized as model instead of ligand

When the user provides a hetcode-named ligand file (e.g., `atp.pdb`) that uses
ATOM records instead of HETATM, the file categorizer classified it as a protein
model. This caused three cascading problems: `has_ligand_file=False`, ligandfit
reported as unavailable ("missing required file: ligand_file"), and `best_files`
selected `atp.pdb` as the model for refinement/validation.

**Root cause:** `_pdb_is_small_molecule()` used a strict HETATM-only test. Some
ligand PDB files (especially those exported from crystal structures) use ATOM
records for all atoms instead of HETATM.

**Fix:** `_pdb_is_small_molecule()` now uses size-based detection aligned with
`_pdb_is_protein_model()`:
- ≤150 total coordinate records → small molecule (regardless of ATOM vs HETATM)
- >150 and HETATM-only → small molecule
- >150 with ATOM records → NOT small molecule (protein)

Both `_categorize_files()` (workflow_state) and `best_files_tracker._is_ligand_file()`
call `_pdb_is_small_molecule()`, so the fix cascades through both code paths.

### Files changed

| File | Change |
|------|--------|
| `agent/session.py` | Fixed regex (2 locations); supplemental file evaluation in `_rebuild_best_files_from_cycles` and `record_result`; duplicate detection respects different input files; `set_project_info` merges original_files on resume; working directory uses `client_working_directory` from session data |
| `agent/file_utils.py` | Fixed `classify_mtz_type` regex; new `matches_exclude_pattern()` |
| `agent/best_files_tracker.py` | Added MTZ/data stage mappings to `STAGE_TO_PARENT` |
| `agent/command_builder.py` | Content guards for model/ligand slots; exclude_patterns on LLM selections; `_last_missing_slots`; `best_files_fallback` for specific-subcategory slots; detailed ligand guard logging; `require_best_files_only` guard |
| `agent/command_postprocessor.py` | `inject_user_params` validates bare keys against strategy_flags allowlist |
| `agent/workflow_state.py` | `_pdb_is_small_molecule()` size-based detection (≤150 atoms = ligand); `_pdb_is_protein_model()` size-based (≤150 atoms = not protein); increased read buffer to 32 KB; placement probe no longer sets `validation_done` |
| `agent/workflow_engine.py` | Ligandfit allowed with pre-refined models; `_apply_directives` defers ligandfit injection until refine runs; refine kept as ligandfit prerequisite; `_get_metric()` casts to float; `_check_metric_condition()` casts to float |
| `agent/advice_preprocessor.py` | Removed instruction leakage from PHIL translation section header |
| `agent/graph_nodes.py` | Fallback diagnostics: per-program build failure tracking, specific stop reasons; BUILD error includes missing slot names; fallback build log capture |
| `agent/planner.py` | `fix_program_parameters` now removes wrong param when target already exists (instead of skipping) |
| `knowledge/parameter_fixes.json` | Added `wavelength → autosol.lambda` for autosol |
| `knowledge/programs.yaml` | Added `refine` to ligandfit ligand slot exclude_patterns |
| `programs/ai_agent.py` | Content guards for model and ligand slots in safety net; advice display at normal level; FALLBACK debug log display; stores `client_working_directory` in session |
| `tests/tst_audit_fixes.py` | 12 new tests (234 total) |

---

## Version 112.69 (Rule D: strip hallucinated bare params + inject_program_defaults)

Two bugs found during integration testing of the `nsf-d2-ligand` scenario.

### Bug 1: `map_type=pre_calculated` not stripped

The LLM hallucinated `map_type=pre_calculated` on a `phenix.refine` command.
The parameter is ambiguous (maps to 5 different PHIL scopes) and causes a
runtime error. The existing Rule C in `sanitize_command` only stripped bare
params for programs with *zero* strategy_flags. `phenix.refine` has 10
strategy_flags, so Rule C didn't fire and bare hallucinated params slipped
through.

**Fix — Rule D:** For programs WITH strategy_flags, strip bare (unscoped, no
dots) `key=value` params not in the program's allowlist. Scoped PHIL params
(e.g., `xray_data.r_free_flags.generate=True`) are preserved because they
contain dots and go through PHIL validation downstream.

### Bug 2: `r_free_flags.generate=True` not reliably present

If the LLM omits `xray_data.r_free_flags.generate=True`, `phenix.refine` runs
without generating R-free flags, which is always wrong for the first refinement
cycle.

**Fix — `inject_program_defaults()`:** New step 4 in the `postprocess_command`
pipeline. Reads `defaults` from `programs.yaml` and appends any missing ones.
If the LLM already includes the parameter, no duplicate is added. Added
`xray_data.r_free_flags.generate: True` as a default for `phenix.refine`.

### Files changed

| File | Change |
|------|--------|
| `agent/command_postprocessor.py` | Added Rule D in `sanitize_command`; added `inject_program_defaults()` function; wired as step 4 in `postprocess_command` |
| `knowledge/programs.yaml` | Added `defaults: {xray_data.r_free_flags.generate: True}` to `phenix.refine` |
| `tests/tst_audit_fixes.py` | Added `test_s5h_rule_d_strips_bare_hallucinated_params` and `test_s5h_inject_program_defaults` (222 total) |

---

## Version 112.68 (Phase 4: Clean execution loop + post-phase cleanup)

Extracted the 200-line post-execution block from `_run_single_cycle` into two
focused methods, and applied cleanup fixes.

### Structural changes

**`_run_single_cycle`** (92 lines, down from 283): Now a clean orchestrator —
start cycle, get command from graph, handle duplicates, record decision, handle
STOP/empty, execute, delegate to `_handle_execution_result`.

**`_handle_execution_result`** (new, 49 lines): Routes to failure or success
path. Success path updates actual program, clears `advice_changed`.

**`_handle_failed_execution`** (new, 150 lines): Full failure pipeline — probe
programs → auto recovery → terminal diagnosis → annotate as recoverable.

### Post-phase cleanup

1. **17 bare `print("[DEBUG ...")` removed** from `ai_agent.py` — development
   debug aids printing unconditionally to stdout on every cycle. Removed from
   `set_defaults` (9 prints), `iterate_agent` (4 prints), `_get_command_for_cycle`
   (4 prints). Replay notification merged into existing `self.vlog.verbose()` call.

2. **SyntaxWarning in session.py** (line 3739): `"_\* Failure..."` had invalid
   escape sequence `\*`. Fixed to `"_\\* Failure..."`.

3. **Redundant `import os as _os_rf`** in `_get_command_for_cycle` → uses
   module-level `os`.

4. **Test path resolution fix** — three tests navigating to
   `command_postprocessor.py` via relative paths broke when
   `_find_ai_agent_path()` resolved via importlib in the PHENIX build tree.
   Fixed to use `_PROJECT_ROOT` directly.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | Extracted `_handle_execution_result` and `_handle_failed_execution` from `_run_single_cycle`; removed 17 debug prints; removed redundant import (5,971 lines) |
| `agent/session.py` | Fixed SyntaxWarning |
| `tests/tst_audit_fixes.py` | Updated 3 structural tests for moved code |

---

## Version 112.67 (Phase 3: Consolidate stop logic & duplicate retries)

The graph's PERCEIVE node is now the sole authority for all stop decisions.
Duplicate retries now route through the normal graph path instead of a parallel
bypass.

### Phase 3a & 3b: Client-side stop checks (already removed)

Both the consecutive-program cap and the directive stop check were already
removed from `_run_single_cycle` in a prior session. The graph's PERCEIVE
already has both:
- `check_directive_stop()` — fires at start of cycle N+1
- `check_consecutive_program_cap()` — fires at start of cycle N+1

**Behavior change:** Directive stop now fires at the START of cycle N+1 instead
of the END of cycle N. This is better — the workflow engine is consulted before
the stop decision, preventing premature stops (e.g. the predict_and_build bug).

### Phase 3c: Duplicate retries through graph

**Problem:** `_retry_duplicate` was a parallel path that rebuilt `session_info`,
gathered files, and called `decide_next_step` directly — duplicating ~80 lines
of `_query_agent_for_command`. Any new field added to `session_info` had to be
added in both places.

**Fix:**

1. Added `duplicate_feedback` parameter to `_query_agent_for_command` — when set,
   feedback string is appended to guidelines before calling `decide_next_step`

2. Refactored `_handle_duplicate_check` — on duplicate detection, builds feedback
   via `_build_duplicate_feedback` (extracted static helper), calls
   `_query_agent_for_command(duplicate_feedback=...)` through the full graph path,
   then applies `_inject_missing_required_files` to the retry command

3. Deleted `_retry_duplicate` method — 82 lines of parallel-path code eliminated

4. Removed `_last_log_content` — was only saved for `_retry_duplicate` to reuse

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | Added `duplicate_feedback` param; refactored `_handle_duplicate_check`; added `_build_duplicate_feedback`; deleted `_retry_duplicate` and `_last_log_content` (6,601 → 6,532 lines) |
| `agent/session.py` | Updated docstring reference from `_retry_duplicate` to `_build_duplicate_feedback` |
| `tests/tst_audit_fixes.py` | Updated `test_s7g` and `test_s8i` for new duplicate retry path |

---

## Version 112.66 (Phases 1–2: Command post-processing migration to BUILD)

Extracted all command post-processing transforms from `ai_agent.py` into
standalone functions in a new `agent/command_postprocessor.py` module, then
migrated them into the graph's BUILD node as the single source of truth. Uses
the Strangler Fig pattern: Phase 1 ran new path alongside old path in shadow
mode; Phase 2 removed old path after validation.

### New module: `agent/command_postprocessor.py`

Four server-safe transforms as standalone functions (no `self`/class
dependencies):

| Function | Purpose |
|----------|---------|
| `sanitize_command()` | Strip placeholder values, blacklisted params, cross-program hallucinations |
| `inject_user_params()` | Append user key=value params missing from command (with scope matching) |
| `inject_crystal_symmetry()` | Append unit_cell/space_group from directives |
| `postprocess_command()` | Single entry point calling all transforms in order |

Each function takes explicit data arguments, making them callable from both the
graph (server-side) and `ai_agent.py` (client-side).

### Phase 2 cutover

Removed all 5 client-side transforms from `_get_command_for_cycle`. Only
`_inject_missing_required_files` remains (needs `os.path.exists()`, client-only).

### Issues found and fixed during implementation

1. **Wrong import path** — `libtbx.langchain.command_postprocessor` →
   `libtbx.langchain.agent.command_postprocessor`

2. **`_retry_duplicate` missing `bad_inject_params`** — retry path built its own
   `session_info` without this field

3. **Transport normalization gap (pre-existing)** — `build_request_v2` allowlist
   was missing `advice_changed`, `unplaced_model_cell`, and `bad_inject_params`.
   All three were silently dropped during transport encoding.

4. **Probe-program file-path stripping (pre-existing)** — `sanitize_command` for
   probe programs stripped ALL `key=value` tokens, including legitimate file
   assignments like `half_map=/path/to/map.mrc`. Fixed by preserving tokens where
   the value contains `/` or has a crystallographic extension.

5. **History replay regression** — replay path bypasses graph; replayed commands
   would get no sanitize/inject. Fixed by calling `postprocess_command` directly
   for replayed commands (detected via `_from_replay` flag).

### Files changed

| File | Change |
|------|--------|
| `agent/command_postprocessor.py` | NEW — 537 lines at Phase 2 cutover |
| `programs/ai_agent.py` | Removed client-side transforms; replay postprocess; `bad_inject_params` in `session_info` |
| `agent/graph_nodes.py` | Added `postprocess_command()` call in `_build_with_new_builder` and `_fallback_with_new_builder` |
| `agent/api_client.py` | Added `bad_inject_params` to `build_session_state`; fixed transport normalization allowlist |
| `agent/graph_state.py` | Added `bad_inject_params: Dict` to `AgentState` and `create_initial_state` |
| `phenix_ai/run_ai_agent.py` | Passes `bad_inject_params` to `create_initial_state` |

---

## Version 112.65 (Systematic audit: file categorization and YAML integrity)

Ran 22 systematic checks across YAML configs and Python code. Found and fixed 5
structural bugs affecting file categorization, category references, and
intermediate file handling.

### Bug 1 — CRITICAL: `parent_category` bubble-up never implemented
**File:** `agent/workflow_state.py`

The YAML defines parent-child relationships (`refined` → `model`,
`refine_map_coeffs` → `map_coeffs_mtz`, `ligand_cif` → `ligand`) but the code
never propagated files from subcategories to their semantic parents. This caused
`model`, `map_coeffs_mtz`, and `ligand` to be permanently empty.

**Fix:** Added bubble-up logic in `_categorize_files_yaml` that propagates files
from child categories to their `parent_category`.

### Bug 2 — `intermediate` category patterns were dead code
**File:** `knowledge/file_categories.yaml`

`intermediate` was marked `is_semantic_parent: true`, which causes both step 1
(extension matching) and step 2 (pattern matching) to skip it entirely.

**Fix:** Removed erroneous `is_semantic_parent` flag.

### Bug 3 — Intermediate files leaked into `model` via `unclassified_pdb`
**Files:** `knowledge/file_categories.yaml`, `agent/workflow_state.py`

Files matching intermediate patterns also match `unclassified_pdb`'s `*`
catch-all, which bubbles up to `model`.

**Fix:** Added missing excludes to `unclassified_pdb` AND added post-processing
that removes any file in `intermediate` from `model`/`search_model`/`pdb`.

### Bug 4 — `phaser_output` lacked extension filter
**File:** `knowledge/file_categories.yaml`

`PHASER*.mtz` matched `phaser_output` (parent: `model`), causing Phaser MTZ
files to bubble up into `model`.

**Fix:** Added `extensions: [".pdb"]` to restrict to PDB files only.

### Bug 5 — Broken category reference
**File:** `knowledge/programs.yaml`

`phenix.phaser.model.exclude_categories` referenced `ligand_fit` which doesn't
exist. Fixed to `ligand_fit_output`.

### Audits that passed clean

All command template `{slot}` placeholders match defined inputs; all workflow
phase transitions reference valid phases; all strategy_flags have proper type
declarations; all invariant `has_strategy` checks have corresponding sources;
all hardcoded `categorized_files.get("xxx")` references in Python match YAML;
no sibling subcategories have overlapping patterns; no category/exclude overlaps
in `input_priorities`; all program output patterns categorize correctly; all
`also_in` and `parent_category` references point to existing categories.

### Files changed

| File | Change |
|------|--------|
| `agent/workflow_state.py` | Implemented `parent_category` bubble-up; intermediate exclusion post-processing |
| `knowledge/file_categories.yaml` | Fixed `intermediate` flag; added `phaser_output` extension filter; added `unclassified_pdb` excludes |
| `knowledge/programs.yaml` | Fixed `ligand_fit` → `ligand_fit_output` reference |
| `tests/tst_audit_fixes.py` | Comprehensive categorization tests |

---

## Version 112.64 (Unit cell / space group reliably propagated to commands)

When the user specifies a unit cell or space group in their advice (e.g.
`"The specified unit cell (116.097, 116.097, 44.175, 90, 90, 120) must be
used for the procedure"`), the value was previously silently ignored — it
never appeared in the command sent to `phenix.model_vs_data` or any other
program.

### Root cause (two independent bugs)

**Bug 1 — LLM returns empty directives despite the correct schema being present.**
The LLM prompt already listed `unit_cell` and `space_group` as extractable
`program_settings` parameters, but when the advice mixes crystal symmetry with
stop-conditions and file preferences, the LLM sometimes returns `{}` for the
entire directive dict (logged as "No actionable directives found"). The
`extract_directives_simple` regex fallback was only called for the ollama
provider, not for the main Google/OpenAI path.

**Bug 2 — Wrong PHIL scope in emitted commands.**
Two code paths (`_inject_crystal_symmetry` in `ai_agent.py`, and the
`KNOWN_PHIL_SHORT_NAMES` PASSTHROUGH in `program_registry.py`) appended bare
`unit_cell="..."` and `space_group=...`, which PHENIX programs may not accept
without the fully-scoped `crystal_symmetry.unit_cell=` form.

### Fix — three layers

**1. Deterministic regex fallback** (`agent/directive_extractor.py`) — New
`_apply_crystal_symmetry_fallback()` function is called inside
`extract_directives()` immediately after `validate_directives()`. It runs the
same regex patterns as `extract_directives_simple` but only fills in fields
that the LLM left empty — it never overwrites a value the LLM extracted
correctly. This means the unit cell is captured even when the LLM returns
`{}`.

**2. Scoped PHIL form in `_inject_crystal_symmetry`** (`programs/ai_agent.py`)
— Changed both append statements from `unit_cell="..."` and `space_group=...`
to `crystal_symmetry.unit_cell="..."` and `crystal_symmetry.space_group=...`.
The fully-scoped form is accepted by all X-ray PHENIX programs.

**3. Scoped PHIL form in `program_registry.py`** — The `KNOWN_PHIL_SHORT_NAMES`
PASSTHROUGH block now maps `unit_cell` and `space_group` to
`crystal_symmetry.unit_cell` / `crystal_symmetry.space_group` before appending
to the command string.

### End-to-end result for the nsf-d2 example

After the fix, even when the LLM returns empty directives, the deterministic
fallback extracts the unit cell and the command becomes:

```
phenix.model_vs_data nsf-d2_noligand.pdb nsf-d2.mtz \
    crystal_symmetry.unit_cell="116.097 116.097 44.175 90 90 120"
```

The same symmetry is also injected into subsequent programs:
`phenix.ligandfit`, `phenix.refine`, `phenix.xtriage`, etc.

### Files changed

| File | Change |
|------|--------|
| `agent/directive_extractor.py` | New `_apply_crystal_symmetry_fallback()` function; called at end of `extract_directives()` after `validate_directives()` |
| `agent/program_registry.py` | `KNOWN_PHIL_SHORT_NAMES` PASSTHROUGH uses `crystal_symmetry.unit_cell=` / `crystal_symmetry.space_group=` scoped form |
| `programs/ai_agent.py` | `_inject_crystal_symmetry()` uses `crystal_symmetry.unit_cell=` and `crystal_symmetry.space_group=` scoped form |
| `tests/tst_audit_fixes.py` | 6 new `test_s4b_*` tests: fallback behavior, non-overwrite, partial fill, source-inspection checks for scoped form |

---

## Version 112.62 (Remove Results page on fatal diagnosis)

When `_diagnose_terminal_failure` fires, `_finalize_session` now skips the
Results summary page entirely.  The diagnosis HTML is already open in the
user's browser — a second Results window saying "1 cycle, 1 failure" adds no
value and risks burying the actionable diagnosis.

### How it works

`_diagnose_terminal_failure` writes `session.data["failure_diagnosis_path"]`
before returning.  `_finalize_session` checks this key:

```python
has_fatal_diagnosis = bool(session.data.get("failure_diagnosis_path"))
if (session.get_num_cycles() > 0
    and not skip_summary
    and not has_fatal_diagnosis   # ← new guard
    and not dry_run):
    self._generate_ai_summary(session)
```

When the key is present, `_generate_ai_summary` (and therefore
`display_results`) is skipped.  The session is still saved to disk and
`self.result` is still populated — only the browser window is suppressed.
A one-line note is printed to the log: `Skipping Results summary — error
diagnosis report already shown. See: <path>`.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | `_finalize_session`: `has_fatal_diagnosis` guard skips summary |
| `tests/tst_audit_fixes.py` | Updated `test_s3a_diagnose_returns_true_no_sorry` and `test_s3a_finalize_runs_after_diagnosis` to verify suppression |

---

## Version 112.61 (Remove Sorry from terminal failure path)

The `_diagnose_terminal_failure` method no longer raises `Sorry`.  Previously
the PHENIX GUI showed a blocking modal error dialog that could end up behind
other windows, preventing the user from proceeding.  The ai_agent job is now
considered successful when it correctly identifies and diagnoses a sub-job
failure.

### Changes

**`_diagnose_terminal_failure`** — Steps 1–4 unchanged.  Step 5 replaced:
- Old: `raise Sorry(short_sorry)` — blocking modal, potentially hidden
- New: `return True` — signals `_run_single_cycle` to break the loop

**`_run_single_cycle`** — Changed from discarding the return value to
propagating it: `return self._diagnose_terminal_failure(...)`.

**`iterate_agent` cycle loop** — Removed the entire `_pending_sorry`
try/except/re-raise scaffold (12 lines).  The loop is now:

```python
for cycle in range(start_cycle, start_cycle + max_cycles):
    should_break = self._run_single_cycle(cycle, session, session_start_time)
    if should_break:
        break
self._finalize_session(session)    # always runs regardless
```

`_finalize_session` remains unconditional.

### User experience

Fatal sub-job failure: diagnosis HTML opens in browser → agent finishes
cleanly → one window, no modals, no buried dialogs.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | `_diagnose_terminal_failure` returns `True`; `_run_single_cycle` propagates it; `iterate_agent` drops `_pending_sorry` |
| `tests/tst_audit_fixes.py` | Replaced `test_s3a_html_write_failure_produces_full_sorry` and `test_s3a_sorry_deferred_until_finalize` with `test_s3a_diagnose_returns_true_no_sorry` and `test_s3a_finalize_runs_after_diagnosis` |

---

## Version 112.60 (Error page and Results page improvements)

Five UX improvements to the terminal error and Results pages.

### Changes

**(1) Heading rename** — `failure_diagnoser.py`: "Terminal Error Diagnosis" →
"Error diagnosis".  Tab title also updated.

**(2) File location on error page** — `build_diagnosis_html` now shows the
full path of the saved HTML file in the footer: `Saved to: <path>`.

**(3) Job context on error page** — `build_diagnosis_html` now accepts
`job_name` and `working_dir`; these appear in the meta bar alongside Program
and Cycle.  `_diagnose_terminal_failure` derives them from `log_dir`.

**(4) Working directory on Results page** — `_extract_summary_data` now
includes the working directory via `_get_working_directory_path()`.
`_format_summary_markdown` shows it as
`**Working directory:** \`/full/path\`` on every run, immediately below the
Session/Cycles line.

**(5) Failure diagnosis reference on Results page** — When a fatal diagnosis
fires, `session.data["failure_diagnosis_path"]` is set and appears in a
"Failure Diagnosis" section of the Results markdown, just above Assessment,
with the full path to `ai_failure_diagnosis.html`.

### Files changed

| File | Change |
|------|--------|
| `agent/failure_diagnoser.py` | `build_diagnosis_html` gains `html_path`, `job_name`, `working_dir` params; heading renamed |
| `programs/ai_agent.py` | `_diagnose_terminal_failure` passes new params; stores `failure_diagnosis_path` in session |
| `agent/session.py` | `_extract_summary_data` adds `working_dir` and `failure_diagnosis_path`; `_format_summary_markdown` renders them |
| `tests/tst_audit_fixes.py` | `test_s3a_build_html_new_fields` verifies heading, path, job name, working dir |

---

## Version 112.59 (Noligand false-positive fix)

Fixes `nsf-d2_noligand.pdb` (a protein structure explicitly lacking a ligand)
being classified as a ligand file due to substring matching: `'ligand.pdb' in
'noligand.pdb'` evaluates to `True`.

### Root cause

`BestFilesTracker._is_ligand_file()` used `any(p in basename for p in
ligand_patterns)` where `ligand_patterns` included `'ligand.pdb'`.  Substring
matching has no concept of word boundaries.

### Fix — word-boundary regex

Replaced substring matching with a regex that requires `ligand` or `lig` to
appear at a word boundary (start of name, or after `_`, `-`, or `.`):

```python
WORD_SEP = r'(?:^|[_\-\.])'
AFTER    = r'(?=[_\-\.]|\.(?:pdb|cif)$|$)'
word_boundary_ligand = re.search(
    r'(?:' + WORD_SEP + r'lig' + AFTER + r'|'
           + WORD_SEP + r'ligand' + AFTER + r')',
    basename
)
```

Matches: `lig.pdb`, `lig_001.pdb`, `ligand.pdb`, `my_ligand.pdb`
Rejects: `nsf-d2_noligand.pdb`, `noligand_model.pdb`

The same pattern was applied to `workflow_state.py`'s hardcoded categorizer
and to `file_categories.yaml`'s `unclassified_pdb` exclude list.

### Files changed

| File | Change |
|------|--------|
| `agent/best_files_tracker.py` | `_is_ligand_file` rewritten with word-boundary regex |
| `agent/workflow_state.py` | `_categorize_files_hardcoded` updated to match |
| `knowledge/file_categories.yaml` | `unclassified_pdb` excludes use specific patterns; `*noligand*` added as positive |
| `tests/tst_audit_fixes.py` | `test_is_ligand_file_noligand_false_positive` (15 cases) |

---

## Version 112.58 (Test path fix)

Fixed `test_s3a_sorry_deferred_until_finalize` using a hard-coded path
(`_PROJECT_ROOT/programs/ai_agent.py`) instead of calling `_find_ai_agent_path()`.
All other tests already used the helper; this was the only outlier.

### Files changed

| File | Change |
|------|--------|
| `tests/tst_audit_fixes.py` | `test_s3a_sorry_deferred_until_finalize`: use `_find_ai_agent_path()` |

---

## Version 112.57 (HETATM-based ligand detection)

Fixes `atp.pdb`, `gdp.pdb`, and similar hetcode-named files being classified as
protein models rather than ligands.  Previously, files named after their PDB
hetcode had no keyword-based ligand signal, so they fell through to the model
category.

### Fix — content-based detection

Added `_pdb_is_small_molecule(path)` which reads the PDB file and returns
`True` when all coordinate records are `HETATM` (no `ATOM` records).  This is
applied as a post-processing pass in three locations:

1. `best_files_tracker.py` `_is_ligand_file()` — content fallback after
   keyword checks
2. `workflow_state.py` YAML categorizer post-processing pass
3. `workflow_state.py` hardcoded categorizer post-processing pass

A file that passes the HETATM test is re-assigned from `model` / `unclassified_pdb`
to `ligand_pdb`, preventing it from being selected as the refinement model.

### Files changed

| File | Change |
|------|--------|
| `agent/workflow_state.py` | `_pdb_is_small_molecule()` helper; post-processing in YAML and hardcoded categorizers |
| `agent/best_files_tracker.py` | HETATM content fallback in `_is_ligand_file()` |
| `tests/tst_audit_fixes.py` | `test_pdb_is_small_molecule_helper`, `test_hetcode_ligand_not_used_as_refine_model` |

---



Fixes the apoferritin AIAgent_165 scenario where the placement probe
(`phenix.map_correlations`) crashed with "model is entirely outside map",
ran a second time with the same result, then caused the agent to quit —
never routing to `dock_in_map`.

### Root cause analysis — three independent failures

**Failure 1 — Tier 1 (cell mismatch) is non-functional in production.**
`_check_cell_mismatch` calls `check_cryoem_cell_mismatch(pdb_path, map_path)`
which immediately calls `os.path.exists(pdb_path)`. On the server, `pdb_path`
is `/Users/terwill/.../1aew_A.pdb` — a client-side path. `os.path.exists`
returns `False`. Function returns `False` (fail-safe). The 4× unit cell
mismatch (184 Å F432 crystal vs 32.5 × 39.65 × 36.4 Å P1 sub-box) is never
detected. **Tier 1 has been silently broken for all RemoteAgent users.**

**Failure 2 — The crash itself carries the answer, but the code discards it.**
The probe runs correctly. `map_correlations` raises:

```
Sorry: Stopping as model is entirely outside map and wrapping=False
```

This is a stronger signal than a low CC — it's categorical proof the model
is not placed. But `_analyze_history` has this at the top of the probe loop:

```python
if _is_failed_result(_result):
    continue   # Ignore failed cycles for probe detection
```

The entire entry is discarded. `placement_probed` stays `False`.

**Failure 3 — The loop.**
With `placement_probed=False`, `placement_uncertain` is still `True` on the
next cycle. The agent runs `map_correlations` again → same crash → same discard
→ eventually the LLM consecutive-failure counter trips and the run stops with
no docking ever attempted.

### Fix 1 — Probe crash detection in `_analyze_history` (workflow_state.py)

Before the `continue`, inspect failed `map_correlations` entries that occurred
before any refine/dock cycle:

```python
if _is_failed_result(_result):
    if "map_correlations" in _ecomb and not _seen_refine_or_dock:
        _rl = (_result or "").lower()
        _outside_signals = [
            "entirely outside map", "outside map", "model is outside",
            "model entirely outside", "stopping as model",
        ]
        if any(s in _rl for s in _outside_signals):
            # Hard evidence: model is not in the map at all
            info["placement_probed"] = True
            info["placement_probe_result"] = "needs_dock"
        elif not info.get("placement_probed"):
            # Unknown failure — prevent infinite probe retry
            info["placement_probed"] = True
            # Leave placement_probe_result as None (inconclusive)
    continue
```

Three outcomes:
- **"outside map" crash** → `placement_probed=True, result="needs_dock"` → routes to `dock_model`
- **Other crash** → `placement_probed=True, result=None` → inconclusive, falls through to `obtain_model`
- **Second failed probe** → `placement_probed` already set, no overwrite; guard on `not info.get("placement_probed")` prevents the inconclusive case from clobbering an earlier definitive result

### Fix 2 — Client-side model cell transport (S2L-b)

The client reads the model's CRYST1 cell (which it has access to) and transmits
it in `session_state["unplaced_model_cell"]`. The server uses this pre-read cell
in `_check_cell_mismatch` instead of trying to open the file:

**Client side** (`programs/ai_agent.py`): Before assembling `session_info`,
read the CRYST1 cell from the first unplaced PDB in `active_files` and add it
to `session_info["unplaced_model_cell"]`. Only populated when placement hasn't
been confirmed by history (no `dock_done`, no `refine_done`).

**Transport** (`agent/api_client.py`): `build_session_state` passes
`unplaced_model_cell` through to `session_state`.

**Server receipt** (`phenix_ai/run_ai_agent.py`): Maps `unplaced_model_cell`
from `session_state` into `session_info`.

**Server use** (`agent/workflow_engine.py`): `build_context` passes
`session_info` down to `_check_cell_mismatch(files, model_cell=...)`.
`_check_cell_mismatch` uses the pre-read cell for comparison against the map
(which the server **can** read — it was just created by `resolve_cryo_em`
on the server). Tier 1 now fires correctly in production:

```
model cell  = (184, 184, 184, 90, 90, 90)   # F432 crystal
map cell    = (32.5, 39.65, 36.4, 90, 90, 90)  # P1 sub-box
→ mismatch > 5% on all three axes → cell_mismatch=True → dock_model
```

### Also fixed: `optimized_full_map` not checked by `_check_cell_mismatch`

After `resolve_cryo_em`, the output `denmod_map.ccp4` is categorized as
`optimized_full_map` (not `full_map`). The original code only checked
`files.get("full_map") or files.get("map")`. Added `optimized_full_map`
as a priority-2 fallback so the freshly-generated density-modified map is
always found.

### Files changed

| File | Change |
|------|--------|
| `agent/workflow_state.py` | `_analyze_history`: probe-crash handler before `continue`; `detect_workflow_state` gains `session_info` parameter |
| `agent/workflow_engine.py` | `_check_cell_mismatch` rewritten with `model_cell` parameter and S2L fast-path; `optimized_full_map` added to map search; `build_context` accepts `session_info`; `get_workflow_state` accepts `session_info` |
| `agent/api_client.py` | `build_session_state` passes `unplaced_model_cell` through |
| `phenix_ai/run_ai_agent.py` | Maps `unplaced_model_cell` from `session_state` → `session_info` |
| `programs/ai_agent.py` | Reads CRYST1 cell client-side, adds to `session_info["unplaced_model_cell"]` |
| `agent/graph_nodes.py` | Passes `session_info=state.get("session_info",{})` to `detect_workflow_state` |
| `tests/tst_audit_fixes.py` | 8 new S2L tests (131 total) |

### Corrected cycle trace for AIAgent_165

| Cycle | Program | Why (after fix) |
|-------|---------|-----------------|
| 1 | `phenix.mtriage` | Map analysis |
| 2 | `phenix.resolve_cryo_em` | Half-maps → denmod_map |
| **3** | **`phenix.dock_in_map`** | **Tier 1: cell_mismatch=True (client cell + server map) → dock_model** |
| 4 | `phenix.real_space_refine` | Model docked, refine |

Without the fix, cycle 3 ran `map_correlations` twice (probe crash not interpreted), then stopped.

### Note on test coverage

Tests 5 and 6 (api_client and workflow_engine) are marked SKIP in non-PHENIX
environments because they require the production libtbx import path. They pass
fully in a PHENIX installation. Tests 1–4, 7–8 pass in all environments.

---

## Version 112.48 (S2k — _inject_user_params Empty Program Guard)

Fixes a subtle Python truth-value bug introduced in v112.47 where the
program-scope filter silently passed all dotted keys through when
`prog_base` was an empty string.

### Problem

```python
# Python: "refinement".startswith("") is True — every string starts with ""
scope_matches = leading_scope.startswith(prog_base) or prog_base.startswith(leading_scope)
```

When `program_name` was not passed to `_inject_user_params`, `prog_base`
defaulted to `""`, so `leading_scope.startswith("")` was always `True` and
every dotted key bypassed the filter — exactly the bug the fix was meant to prevent.

### Fix

Added three guards to the scope-matching expression:

```python
scope_matches = (
    bool(prog_base) and len(prog_base) >= 4 and
    (leading_scope == prog_base or
     (leading_scope.startswith(prog_base) and len(prog_base) >= 4) or
     (prog_base.startswith(leading_scope) and len(leading_scope) >= 4))
)
```

- `bool(prog_base)` — empty string immediately fails; no keys injected
- `len(prog_base) >= 4` — prevents one- or two-character accidental prefix matches
- Prefix matching handles `refine` ↔ `refinement` in both directions

**Decision table:**

| Program | Key | Result |
|---|---|---|
| `phenix.refine` | `refinement.main.number_of_macro_cycles` | **INJECT** (`refine` ⊂ `refinement`) |
| `phenix.ligandfit` | `refinement.main.number_of_macro_cycles` | **SKIP** (no overlap) |
| `phenix.autosol` | `autosol.atom_type` | **INJECT** (exact match) |
| any | `general.nproc` | **INJECT** (universal scope) |
| (empty string) | `refinement.*` | **SKIP** (`bool("")` fails immediately) |

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | `_inject_user_params`: three-part guard on scope-matching expression |
| `tests/tst_audit_fixes.py` | 1 new S2k guard test (123 total) |

---

## Version 112.47 (S2k — _inject_user_params Program-Scoped Filtering)

Extends `_inject_user_params` to accept a `program_name` parameter and filter
dotted PHIL keys to only inject parameters that belong to the target program's
scope.

### Problem

After the server correctly built `phenix.ligandfit ... general.nproc=4`, the
client's `_inject_user_params` scanned the guidelines string, found
`refinement.main.number_of_macro_cycles=2` (from an earlier user directive),
and appended it unconditionally:

```
[inject_user_params] appended: refinement.main.number_of_macro_cycles=2
Final command: phenix.ligandfit ... general.nproc=4 refinement.main.number_of_macro_cycles=2
```

The command was clean when it left the server; the client contaminated it.

### Fix

`_inject_user_params(self, command, guidelines, program_name='')` now
classifies each extracted dotted key before injecting:

```python
_UNIVERSAL_SCOPES = {'general', 'output', 'job', 'data_manager', 'nproc'}
prog_base = program_name.replace('phenix.', '').lower()

for key in extracted_keys:
    if '.' in key:
        leading_scope = key.split('.')[0].lower()
        if leading_scope not in _UNIVERSAL_SCOPES and not scope_matches(leading_scope, prog_base):
            skipped.append(key)
            continue
    inject(key)
```

Universal scopes (`general`, `output`, etc.) are always injected. All other
dotted keys are only injected when their leading scope matches the program name.

**Note:** v112.47 contained the empty-`prog_base` bug fixed in v112.48.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | `_inject_user_params` rewritten with `program_name` parameter and scope filter; call site updated to pass `program_name` |
| `tests/tst_audit_fixes.py` | 1 new S2k test (122 total) |

---

## Version 112.46 (S2k — _inject_user_params STOP Guard)

Prevents `_inject_user_params` from running at all when the command is `STOP`.

### Problem

The `_inject_user_params` call site had no guard for STOP commands. A user
directive containing `refinement.main.number_of_macro_cycles=2` would cause
the function to emit `STOP refinement.main.number_of_macro_cycles=2`, which
then failed command validation.

### Fix

```python
if command and command != 'No command generated.' and \
        command.strip().split()[0] != 'STOP':
    command = self._inject_user_params(command, guidelines, program_name)
```

STOP commands bypass `_inject_user_params` entirely.

### Files changed

| File | Change |
|------|--------|
| `programs/ai_agent.py` | Guard added at `_inject_user_params` call site |

---

## Version 112.45 (S2j — program_registry Passthrough Removal)

Removes the dotted-key passthrough in `program_registry.py` that allowed
LLM strategy parameters to leak across program boundaries.

### Problem

`_build_command_from_registry` had a catchall passthrough: any key in the
LLM strategy dict that contained a dot was appended to the command verbatim.
This meant `refinement.main.number_of_macro_cycles=2` from the LLM's strategy
for one program was silently included in the command for any program.

### Fix

Replaced the unconditional passthrough with an allowlist. A strategy key
is now accepted only when it:

1. Appears in `strategy_flags` for this specific program, or
2. Is in `KNOWN_PHIL_SHORT_NAMES` (`nproc`, `twin_law`, `unit_cell`, etc.), or
3. Contains a literal `=` sign (already a complete PHIL assignment)

Dotted-path keys that pass none of these tests are logged as `DROPPED` and
discarded. User-supplied dotted overrides reach the command through
`_inject_user_params` instead (subject to the program-scope filter added in
v112.47/v112.48).

### Files changed

| File | Change |
|------|--------|
| `agent/program_registry.py` | Passthrough logic replaced with allowlist; dropped keys logged |
| `tests/tst_audit_fixes.py` | 1 new S2j code test (121 total) |

---

## Version 112.44 (S2j — Cross-Program Strategy Contamination: Prompt Fix)

Fixes the LLM prompt so the model only emits strategy parameters that apply to
the program it is currently selecting.

### Root cause

Three places in the prompt instructed the LLM to extract user parameters into
the `strategy` field with no program-specificity qualifier. The LLM obediently
included refinement parameters (`refinement.main.number_of_macro_cycles=2`) in
the strategy for every program it selected — including `phenix.ligandfit` and
even `STOP`.

### Fixes (knowledge/prompts_hybrid.py)

**User advice extraction block** (was "Extract any specific parameters from the
user advice and include them in your strategy field"):

> Extract parameters **ONLY when they apply to the program you are currently
> selecting.** Do NOT include parameters for a different program. For example,
> if selecting `phenix.ligandfit` and the user mentioned
> `refinement.main.number_of_macro_cycles`, do NOT include that parameter —
> it applies to `phenix.refine`, not `phenix.ligandfit`.

**OUTPUT FORMAT strategy field** — Added CRITICAL STRATEGY RULE block:

> Strategy keys must ONLY contain parameters valid for the selected program.
> NEVER include parameters from a different program. When stop=true, strategy
> must be empty {}.

**IMPORTANT RULES** — Added rule 6:

> Strategy is program-specific: never put parameters for program X in a
> strategy for program Y.

### Files changed

| File | Change |
|------|--------|
| `knowledge/prompts_hybrid.py` | User advice extraction, OUTPUT FORMAT strategy, IMPORTANT RULES rule 6 |
| `tests/tst_audit_fixes.py` | 1 new S2j prompt test (120 total) |

---

## Version 112.43 (S2i — STOP Command Normalization, Root Cause Fix)

Fixes the root cause of `STOP refinement.main.number_of_macro_cycles=2` being
generated as a command.

### Root cause analysis

The LLM's OUTPUT FORMAT schema has no `command` field — only `program`,
`strategy`, `files`, `stop`. When the LLM returned:

```json
{"program": "STOP", "strategy": {"refinement.main.number_of_macro_cycles": 2}, "stop": false}
```

The PLAN node normalized `program` to `"STOP"` but never set
`intent["stop"] = True`. The BUILD node then saw `stop=False`, fell through the
STOP guard, and assembled the strategy flags onto the command as normal
arguments, producing `STOP refinement.main.number_of_macro_cycles=2`, which
then failed validation as "not a recognized phenix program."

### Fix 1 — PLAN node (agent/graph_nodes.py)

When `chosen_program == "STOP"`, explicitly set `intent["stop"] = True`:

```python
if chosen_program == "STOP" or (intent.get("stop") and
        chosen_program not in valid_programs):
    chosen_program = "STOP"
    intent["program"] = "STOP"
    intent["stop"] = True        # NEW: always set stop=True when program is STOP
```

### Fix 2 — BUILD node (defence in depth)

Both `_build_with_new_builder` and the legacy `build()` function now short-circuit
when either `intent["stop"]` is True or `intent["program"]` is `"STOP"`:

```python
if intent.get("stop") or intent.get("program") == "STOP":
    state = _log(state, "BUILD: Stop requested, no command needed")
    return {**state, "command": "STOP"}
```

This means the BUILD node never sees strategy flags for STOP, regardless of
how the LLM filled the `stop` boolean.

### Files changed

| File | Change |
|------|--------|
| `agent/graph_nodes.py` | PLAN: `intent["stop"] = True` when program is STOP; BUILD: early-exit guard in both build paths |
| `tests/tst_audit_fixes.py` | 4 new S2i tests (119 total) |

---

## Version 112.42 (S2h — validation_cryoem DataManager PHIL Fix)

Resolves a PHIL argument error in `phenix.validation_cryoem` when run with
DataManager-style file arguments.

### Problem

`phenix.validation_cryoem` expected PHIL-scoped arguments
(`input.pdb.file_name=model.pdb`) but the command builder was emitting
positional DataManager style (`model.pdb map.ccp4`), causing a PHIL parse
error on launch.

### Fix

Added the correct PHIL argument template to `programs.yaml` for
`phenix.validation_cryoem`, matching the argument style the program actually
accepts. The command builder picks this up automatically through the standard
invariants pipeline.

### Files changed

| File | Change |
|------|--------|
| `knowledge/programs.yaml` | `phenix.validation_cryoem`: correct PHIL argument template |
| `tests/tst_audit_fixes.py` | 2 new S2h tests (115 total) |

---

## Version 112.41 (S2g — map_correlations Conditional Resolution)

Fixes `phenix.map_correlations` being handed a `resolution=` argument it does
not accept when running in map-vs-map mode.

### Problem

The command builder always injected `resolution=` from the session context
into `map_correlations` commands. In model-vs-map mode this is harmless; in
map-vs-map mode (used by the placement probe) the program does not accept a
resolution argument and exits with an error.

### Fix

Added a conditional resolution invariant to `map_correlations` in
`programs.yaml`: resolution is only included when a model file is present
in the command's file list.

### Files changed

| File | Change |
|------|--------|
| `knowledge/programs.yaml` | `phenix.map_correlations`: conditional resolution invariant |
| `tests/tst_audit_fixes.py` | 2 new S2g tests (113 total) |

---

## Version 112.40 (S2f — validation_cryoem Auto-Resolution)

Adds automatic resolution filling for `phenix.validation_cryoem`, which
requires an explicit `resolution=` argument and previously had to rely on the
user supplying it.

### Problem

`phenix.validation_cryoem` exited with "resolution not provided" when the
agent did not include a resolution argument. The session always has the
resolution from a prior `mtriage` run, but no invariant was wiring it through.

### Fix

Added a `requires_resolution` invariant to `phenix.validation_cryoem` in
`programs.yaml`. The command builder already handles this invariant for other
programs; the fix was purely a YAML addition.

### Files changed

| File | Change |
|------|--------|
| `knowledge/programs.yaml` | `phenix.validation_cryoem`: `requires_resolution: true` invariant |
| `tests/tst_audit_fixes.py` | 2 new S2f tests (111 total) |

---

## Version 112.36 (S2e — after_program Directive Correctly Suppresses Placement Probe)

Fixes a regression introduced by S2b (v112.33). The S2b fix made
`placement_uncertain` immune to the `model_is_placed` workflow preference
(which can be hallucinated by the LLM). However, it also inadvertently made
it immune to `after_program` directives, which are a reliable signal.

### Problem

`tst_workflow_state.py::test_model_placement_inferred_from_model_vs_data_directive`
failed with `AssertionError: Expected NOT xray_analyzed when model_vs_data requested`.

Scenario: user uploads `1aba.pdb + 1aba.mtz`, requests
`after_program=phenix.model_vs_data`. After S2b, `placement_uncertain=True`
because `has_placed_model_from_history=False` (no dock/refine in history).
This routes to `probe_placement → xray_analyzed` instead of the expected
`refine → xray_refine`.

### Root cause

S2b correctly distinguished two types of "placed" evidence:

| Source | Reliable? | Used by S2b |
|--------|-----------|-------------|
| `model_is_placed: True` (workflow_preferences) | No — LLM hallucination | Ignored by `placement_uncertain` ✓ |
| History/file evidence (dock_done, refined PDB, etc.) | Yes | `has_placed_model_from_history` ✓ |
| `after_program` in programs_requiring_placed | Yes — explicit user request | **Not distinguished — treated as unreliable** ✗ |

### Fix

Added a third placement signal: `has_placed_model_from_after_program`.

New method `_has_placed_model_from_after_program(files, directives)`:
- Returns True when `after_program` is in `programs_requiring_placed`
  (`phenix.refine`, `phenix.model_vs_data`, `phenix.polder`, etc.) AND a
  non-ligand, non-search-model PDB is present
- Returns False for all other cases (no directive, non-placement program, no PDB)

Added to `build_context()` as `context["has_placed_model_from_after_program"]`.

Added to `placement_uncertain` formula:
```python
context["placement_uncertain"] = (
    not context["has_placed_model_from_history"] and
    not context["has_placed_model_from_after_program"] and   # NEW
    not context["cell_mismatch"] and
    ...
)
```

### What is and isn't suppressed

| Directive | Suppresses probe? | Why |
|-----------|-------------------|-----|
| `workflow_preferences: {model_is_placed: True}` | No | LLM hallucination risk (S2b) |
| `stop_conditions: {after_program: phenix.model_vs_data}` | Yes | Explicit program request |
| `stop_conditions: {after_program: phenix.predict_and_build}` | No | Not in programs_requiring_placed |
| History: dock_done / refine_done | Yes | Objective evidence |

### Files changed

| File | Change |
|------|--------|
| `agent/workflow_engine.py` | New `_has_placed_model_from_after_program()` method; `"has_placed_model_from_after_program"` added to `build_context()`; added to `placement_uncertain` formula |
| `tests/tst_audit_fixes.py` | 3 new S2e tests (110 total) |

---

## Version 112.35 (S2d — skip_map_model_overlap_check for real_space_refine)

Adds `skip_map_model_overlap_check=True` as a permanent default for every
`phenix.real_space_refine` run.

### Problem

`phenix.real_space_refine` performs a box/symmetry compatibility check before
refining. When a docked model's CRYST1 record does not match the cryo-EM map
box (common for crystal-structure templates that have been docked), this check
raises a hard error and aborts the run before refinement even begins.

### Fix

Added to `knowledge/programs.yaml` under `phenix.real_space_refine`:

```yaml
defaults:
  skip_map_model_overlap_check: True
```

The `defaults` section maps directly to `key=value` arguments appended to every
command. This flag suppresses the overlap check for all RSR runs — no code
changes required.

### Files changed

| File | Change |
|------|--------|
| `knowledge/programs.yaml` | `defaults: skip_map_model_overlap_check: True` added to `phenix.real_space_refine` |
| `tests/tst_audit_fixes.py` | 1 new S2d test (107 total) |

---

## Version 112.34 (S2c — Crystal PDB to search_model Promotion for Cryo-EM Docking)

Fixes the second stage of the apoferritin AIAgent_104 crash. After S2b (v112.33)
correctly routes unplaced crystal-structure PDBs to `dock_model`, the agent was
immediately stuck because `phenix.dock_in_map` was blocked by a condition that
required AlphaFold post-processing output (`processed_model`). This version
makes the complete path from detection → docking → refinement work.

### Background

The live failure had two stages:

1. **Stage 1 (S2b):** A hallucinated `model_is_placed: True` directive suppressed
   both the cell-mismatch check and the placement-uncertainty probe, allowing
   `phenix.real_space_refine` to run directly against an unplaced model →
   symmetry/box mismatch crash. Fixed in v112.33.

2. **Stage 2 (S2c, this version):** After S2b, detection fires correctly and
   routes to `dock_model`. But `dock_in_map` requires `has_search_model=True`,
   and `1aew_A.pdb` (a plain crystal-structure PDB) is categorised as
   `unclassified_pdb → model` — never `search_model`. Three layered problems
   block docking from ever running.

### Three problems fixed

**Problem 1 — Wrong YAML condition in `dock_model` phase (`knowledge/workflows.yaml`)**

The `dock_model` phase was designed exclusively for the AlphaFold stepwise path
(`predict_and_build → process_predicted_model → dock_in_map`). Its condition
required `processed_model`, which is never set for a plain crystal-structure PDB.

```yaml
# Before:
- has: processed_model

# After:
- has_any: [search_model, processed_model]
```

`has_any` was already supported by `_check_conditions`; no engine changes needed.
The AlphaFold path is unaffected: `processed_predicted` bubbles to `search_model`,
so `has_search_model=True` still holds.

**Problem 2 — File never in `search_model` (`agent/workflow_engine.py`)**

The filename-based categoriser runs at session start with no workflow context. A
plain PDB like `1aew_A.pdb` cannot be distinguished from a positioned model by
name alone, so it lands in `unclassified_pdb → model`. At runtime, once context
confirms docking is needed, the new `_promote_unclassified_for_docking()` method
copies the file into `files["search_model"]`.

Promotion fires when **all** of these hold:
- `experiment_type == "cryoem"` (X-ray unaffected)
- `files["unclassified_pdb"]` non-empty
- `not has_placed_model_from_history` (no dock/refine in session history)
- `not has_search_model` (already populated → no-op)
- Any one of: `placement_uncertain` (Tier 3 pre-probe), `placement_probed AND
  needs_dock` (Tier 3 post-probe), or `cell_mismatch AND not from_history` (Tier 1)

The method never mutates its input dicts. `files["model"]` and
`files["unclassified_pdb"]` are left intact.

**Problem 3 — Promoted files discarded before reaching command builder (`agent/workflow_state.py`)**

`get_workflow_state()` returned the promoted `files` dict, but
`detect_workflow_state()` immediately overwrote `state["categorized_files"]`
with the original pre-promotion `files`. All downstream consumers (PLAN, BUILD,
`CommandContext`) read from `workflow_state["categorized_files"]` and would have
seen an empty `search_model` even after the fix to Problem 2.

```python
# Before (always overwrites with original):
state["categorized_files"] = files

# After (respects promoted files when present; falls back only when key absent):
state["categorized_files"] = state.get("categorized_files", files)
```

Note: `... or files` would be wrong here — Python treats `{}` as falsy and would
overwrite a legitimately empty dict. `.get(key, default)` only falls back when the
key is strictly absent from the dict.

### Data flow after fix

```
detect_workflow_state()
  files = _categorize_files()           # unclassified_pdb=[1aew_A.pdb], search_model=[]
  engine.get_workflow_state()
      build_context(files, ...)          # placement_uncertain=True OR cell_mismatch=True
      _promote_unclassified_for_docking  # files["search_model"]=[1aew_A.pdb]
      detect_phase()                     # dock_model (Tier 1 or Tier 3)
      get_valid_programs()               # dock_in_map: has_any([search_model✓]) + full_map✓
      return {"categorized_files": files, ...}   # promoted files in return dict
  state["categorized_files"] = state.get("categorized_files", files)  # kept ✓
→ PLAN reads categorized_files["search_model"] = [1aew_A.pdb]
→ BUILD / CommandContext: find_in_categories("search_model") → 1aew_A.pdb ✓
→ phenix.dock_in_map 1aew_A.pdb denmod_map.ccp4 resolution=1.9
```

### Complete cycle trace (apoferritin scenario)

| Cycle | Program | How |
|-------|---------|-----|
| 1 | `phenix.mtriage` | Map quality analysis |
| 2 | `phenix.resolve_cryo_em` | Phase 1.5: half-maps → full map |
| 3A | `phenix.dock_in_map` | Tier 1: cell mismatch (production, libtbx available); S2c 5c |
| 3B | `phenix.map_correlations` | Tier 3: no libtbx fallback, probe fires; S2c 5a |
| 4B | `phenix.dock_in_map` | Post-probe needs_dock; S2c 5b |
| 4A/5B | `phenix.real_space_refine` | Model docked, ready to refine ✓ |

### Files changed

| File | Change |
|------|--------|
| `knowledge/workflows.yaml` | `dock_model` `dock_in_map`: `has: processed_model` → `has_any: [search_model, processed_model]` |
| `agent/workflow_engine.py` | New `_promote_unclassified_for_docking()` method; call site after `build_context`; `"categorized_files"` added to return dict |
| `agent/workflow_state.py` | Line 1187: `.get("categorized_files", files)` |
| `agent/graph_nodes.py` | PERCEIVE log for S2c promotion (shows which files were promoted) |
| `tests/tst_audit_fixes.py` | 6 new S2c tests; missing `__main__` runner block restored; SKIP guards added to 32 pre-existing tests that lacked them (no behaviour change, prevents false failures in non-PHENIX environments) |

### Tests

6 new tests in `tst_audit_fixes.py` (106 total, all passing):

| Test | Verifies |
|------|----------|
| `test_s2c_promotion_fires_when_placement_uncertain` | Condition 5a: Tier 3 pre-probe path |
| `test_s2c_promotion_fires_when_probe_says_needs_dock` | Condition 5b: Tier 3 post-probe path |
| `test_s2c_promotion_fires_when_cell_mismatch` | Condition 5c: Tier 1 production path |
| `test_s2c_no_promotion_when_placed_by_history` | Guard: `has_placed_model_from_history` blocks promotion |
| `test_s2c_no_promotion_for_xray` | Guard: X-ray experiment type blocked |
| `test_s2c_categorized_files_propagates_through_get_workflow_state` | Structural fix: promoted files reach command builder |

---

## Version 112.33 (S2b — Directive-Immune Placement Uncertainty)

Fixes the first stage of the apoferritin AIAgent_104 crash. The directive
extractor hallucinated `model_is_placed: True` for "solve the structure",
which triggered a suppression cascade that bypassed both the Tier 1 cell-mismatch
check and the Tier 3 placement probe — allowing `phenix.real_space_refine` to
run directly against an unplaced crystal model in a cryo-EM map.

### Root cause

Two expressions in `workflow_engine.py` used `has_placed_model` (which respects
the `model_is_placed` directive) where they should have used
`has_placed_model_from_history` (directive-immune):

1. **S1 short-circuit** — skips the expensive `phenix.show_map_info` subprocess
   once placement is known. A wrong directive zeroed out `cell_mismatch` before
   Tier 1 could act on it.

2. **`placement_uncertain` formula** — the Tier 3 probe gate. A wrong directive
   set this to `False`, suppressing the probe entirely.

### Fix

Both expressions changed from `has_placed_model` to
`has_placed_model_from_history`. The S1 short-circuit and the probe gate are now
immune to directive mistakes. Only history/file evidence (dock/refine completed,
docked PDB present) can suppress these safety checks.

### Files changed

| File | Change |
|------|--------|
| `agent/workflow_engine.py` | S1 short-circuit: `has_placed_model` → `has_placed_model_from_history`; `placement_uncertain` formula: same change |

### Tests

Covered by existing S2 and R2/R3 tests in `tst_audit_fixes.py`.

---

## Version 112.32 (Probe-Based Model Placement Detection)

Adds a three-tier decision framework that determines whether a supplied atomic
model is already placed in the unit cell / map before choosing between
refinement and MR/docking.  Previously the agent had a blind spot for generic
PDB files with no history and no positioning metadata.

### Background

When a user supplies `model.pdb + data.mtz` (or `model.pdb + map.ccp4`) with
no session history, the agent must decide: is the model already placed (→
refine) or does it need to be placed first (→ MR / docking)?  The old
heuristics relied on file subcategory (`positioned`) or history flags —
neither of which are set for a freshly uploaded PDB.

### Three-tier framework

**Tier 1 — Unit cell comparison (free, instant)**
- Reads CRYST1 from PDB (`read_pdb_unit_cell`) and cell from MTZ/map
  (`read_mtz_unit_cell`, `read_map_unit_cells`)
- Compatible within 5% → falls through to Tier 2
- Incompatible → model cannot be placed here → immediately routes to MR / docking
- Fail-safe: any parse failure returns `False` (no mismatch declared)

**Tier 2 — Existing heuristics (`_has_placed_model`)**
- History flags, file subcategory, user directive
- Clear evidence → done; still ambiguous → Tier 3

**Tier 3 — Diagnostic probe (one program cycle)**
- X-ray: runs `phenix.model_vs_data`; R-free < 0.50 → placed
- Cryo-EM: runs `phenix.map_correlations`; CC > 0.15 → placed
- Probe never repeats: `placement_probed` flag persists in history
- Probe result is detected positionally (first occurrence before any
  refine/dock cycle) with no schema change to history entries

### New module: `agent/placement_checker.py`

Public API:

| Function | Purpose |
|---|---|
| `read_pdb_unit_cell(path)` | Parse CRYST1 line → 6-tuple or None |
| `read_mtz_unit_cell(path)` | iotbx.mtz → 6-tuple, falls back to mtzdump |
| `read_map_unit_cells(path)` | `phenix.show_map_info` → full-map and present-portion cells |
| `cells_are_compatible(a, b, tolerance=0.05)` | Fractional comparison; None → True |
| `check_xray_cell_mismatch(pdb, mtz)` | True only when both readable AND incompatible |
| `check_cryoem_cell_mismatch(pdb, map)` | True only when both readable AND model matches neither map cell |

### Changes to `agent/workflow_engine.py`

- `build_context()` gains four keys: `cell_mismatch`, `placement_probed`,
  `placement_probe_result`, `placement_uncertain`
- `placement_uncertain` is `True` when: model + data present, no history evidence,
  no directive, not predicted, no cell mismatch, probe not yet run
- When `placement_probe_result == "placed"`, `build_context` overrides
  `has_placed_model = True` so normal refine routing takes over
- `_check_cell_mismatch()` private method wires Tier 1 into context building
- `_detect_xray_phase()` / `_detect_cryoem_phase()` each gain three routing
  blocks: Tier 1 mismatch → MR/dock, Tier 3 result → MR/dock or fall-through,
  Tier 3 uncertain → `probe_placement` phase
- `probe_placement` added to both `XRAY_STATE_MAP` and `CRYOEM_STATE_MAP`

### Changes to `agent/workflow_state.py`

- `_analyze_history` initialises `placement_probed=False` and
  `placement_probe_result=None`
- Post-loop pass detects probe by position: `model_vs_data` or
  `map_correlations` found before the first refine/dock cycle
- Only successful cycles contribute (failed runs ignored)
- `cc_volume` used as fallback if `cc_mask` absent

### Changes to `knowledge/workflows.yaml`

- `probe_placement` phase added to both `xray` and `cryoem` workflows
  (positioned between `analyze` and `obtain_model`)
- X-ray probe: `phenix.model_vs_data`, transitions `if_placed → refine`,
  `if_not_placed → molecular_replacement`
- Cryo-EM probe: `phenix.map_correlations`, transitions `if_placed → refine`,
  `if_not_placed → dock_model`

### Changes to `knowledge/programs.yaml`

- `phenix.map_correlations`: added `done_tracking` (was missing entirely —
  pre-existing bug where running it in the validate phase never set
  `validation_done`)
- `phenix.model_vs_data`: added clarifying comment explaining phase-aware
  override (flag unchanged: `validation_done`)

### New tests (34 across all steps, in `tests/tst_audit_fixes.py`)

| Category | Count | What they verify |
|---|---|---|
| R1 | 11 | `placement_checker.py` unit cell parsing and comparison |
| R2 | 11 | `workflow_engine.build_context()` new keys and `placement_uncertain` logic |
| R3 | 12 | Phase routing: probe offered, Tier 1 mismatch bypass, probe results route correctly, probe not re-run |
| R3-extra | 10 | Edge cases: failed probe ignored, CC fallback fields, cryo-EM paths, `build_context` override, `placement_uncertain` clears |

### Directive override protection (S2 fixes — applied after log analysis)

**Root cause identified from runtime log:** the directive extractor LLM set
`model_is_placed: True` from the advice "solve the structure" — a case the
prompt explicitly said should NOT trigger that flag. This cascaded:
`_has_placed_model()` returned True → Tier 1 routing checked
`cell_mismatch AND NOT has_placed_model` → False (model "appeared placed") →
skipped docking → `phenix.real_space_refine` failed with
"Symmetry and/or box (unit cell) dimensions mismatch".

**S2 Fix A — Directive extractor prompt (`agent/directive_extractor.py`)**
- `model_is_placed` is now labelled HIGH-PRECISION: when in doubt, do NOT set it.
- Explicit DO NOT list expanded: "solve the structure", "refine this model",
  "run refinement", "fit a ligand", PDB + cryo-EM map without explicit placement
  confirmation, generic/ambiguous goals.
- Added prominent note: a PDB alongside cryo-EM maps always requires docking first.

**S2 Fix B — `has_placed_model_from_history` context key (`agent/workflow_engine.py`)**
- New `_has_placed_model_from_history()` method — identical logic to
  `_has_placed_model()` but intentionally ignores directives; returns True only
  from history flags (`dock_done`, `phaser_done`, `autobuild_done`,
  `predict_full_done`, `refine_done`) and positioned file subcategories.
- New `has_placed_model_from_history` key in `build_context()` output.
  `has_placed_model` remains the directive-inclusive version for program gating.

**S2 Fix C — Tier 1 routing uses `has_placed_model_from_history`**
- Both `_detect_xray_phase` (MR) and `_detect_cryoem_phase` (dock) changed from
  `not context["has_placed_model"]` to `not context.get("has_placed_model_from_history")`.
- A directive claiming the model is placed cannot override a definitive cell-dimension
  mismatch. Only concrete history evidence (dock ran, phaser ran, etc.) suppresses the check.

**S2 Fix D — S1 short-circuit also uses `has_placed_model_from_history`**
- The subprocess short-circuit (skip `phenix.show_map_info` after placement resolved)
  now checks `has_placed_model_from_history` instead of `has_placed_model`.
  A wrong directive no longer suppresses the check on subsequent cycles either.

### New tests (S2 category, 10 tests in `tests/tst_audit_fixes.py`)

| Test | What it covers |
|---|---|
| `s2_has_placed_model_from_history_method_exists` | Method present on WorkflowEngine |
| `s2_from_history_false_when_only_directive` | Directive cannot fool `_has_placed_model_from_history` |
| `s2_from_history_true_when_dock_done` | `dock_done` history → True |
| `s2_context_has_placed_from_history_key` | Key present in `build_context()` output |
| `s2_directive_model_is_placed_does_not_suppress_cell_mismatch` | Core routing fix: cell_mismatch → dock despite directive |
| `s2_history_placed_does_suppress_cell_mismatch` | `dock_done` history legitimately suppresses re-dock |
| `s2_xray_tier1_uses_from_history` | X-ray MR routing: directive cannot block Tier 1 |
| `s2_short_circuit_uses_from_history_not_directive` | Source-level: short-circuit references correct key |
| `s2_directive_prompt_stronger_do_not_set` | Prompt contains explicit DO NOT cases |
| `s2_full_cryoem_stack_routes_to_dock_not_rsr` | **Regression test for the apoferritin bug** |

### Polish fixes (applied after initial implementation review)

Four issues discovered during code review and corrected before release.

**S1 Fix 1 — YAML validator warnings (`agent/yaml_tools.py`)**
- `if_placed` and `if_not_placed` were not in `valid_transition_fields`, generating
  4 spurious "unknown transition field" warnings every time `_validate_workflows()` ran.
- Added both keys to the set.

**S1 Fix 2 — Redundant import (`agent/workflow_state.py`)**
- The probe detection block used `import re as _re2` inside the loop body.
  `re` is already imported at module level.
- Replaced `_re2.search(...)` with `re.search(...)`.

**S1 Fix 3 — Missing local import fallback (`agent/workflow_engine.py`)**
- `_check_cell_mismatch` only tried the `libtbx.langchain.agent.placement_checker`
  import path; any environment without the libtbx namespace (tests, local dev) would
  silently return `False` without trying the bare `agent.placement_checker` path.
- Added a second `except ImportError` branch matching the pattern used everywhere
  else in the codebase.

**S1 Fix 4 — Subprocess per cycle (`agent/workflow_engine.py`)**
- `_check_cell_mismatch` ran `phenix.show_map_info` as a subprocess on every
  `build_context()` call.  For cryo-EM workflows with many cycles this is wasteful:
  once placement is resolved the check can never change the outcome.
- Added a post-processing short-circuit in `build_context`: after the context dict
  is fully built, if `has_placed_model=True` **or** `placement_probed=True`,
  `cell_mismatch` is forced to `False`.  The check still runs on the first cycle
  when placement is genuinely unknown.  The routing conditions
  (`cell_mismatch AND not has_placed_model`) provide a second safety net.

### New tests (S1 category, 10 tests in `tests/tst_audit_fixes.py`)

| Test | Fix covered |
|---|---|
| `s1_yaml_validator_no_if_placed_warnings` | Fix 1: no warnings from probe_placement transitions |
| `s1_yaml_validator_if_placed_is_in_valid_set` | Fix 1: if_placed/if_not_placed in source |
| `s1_no_redundant_import_re_in_probe_block` | Fix 2: no _re2 alias, module-level re used |
| `s1_probe_re_fallback_still_works` | Fix 2: regex fallback still parses r_free from result text |
| `s1_local_import_fallback_in_check_cell_mismatch` | Fix 3: local path fallback present |
| `s1_placement_checker_importable_locally` | Fix 3: all public functions importable without libtbx |
| `s1_cell_mismatch_short_circuits_when_placed` | Fix 4: False when has_placed_model=True |
| `s1_cell_mismatch_short_circuits_when_probed` | Fix 4: False when placement_probed=True |
| `s1_cell_mismatch_not_short_circuited_first_cycle` | Fix 4: check active on first cycle |
| `s1_short_circuit_order_before_probe_override` | Fix 4: short-circuit before probe-result override |

### Files modified

- `agent/placement_checker.py` — **NEW**
- `agent/workflow_engine.py` — plus S1 fix 3 (import fallback), fix 4 (short-circuit)
- `agent/workflow_state.py` — plus S1 fix 2 (import cleanup)
- `agent/yaml_tools.py` — S1 fix 1 (transition field set)
- `knowledge/workflows.yaml`
- `knowledge/programs.yaml`
- `tests/tst_audit_fixes.py`

---

## Version 112.31 (Session Management, Resume Enhancement, Completed-Workflow Extension)

### P1: Session management keywords populate `self.result` (Fix 27)

**`display_and_stop` / `remove_last_n` left `self.result` unset — GUI calls failed**
- When the agent exited via `_handle_session_management()`, it returned without
  populating `self.result`. Any downstream call to `get_results()` or
  `get_results_as_JSON()` raised `AttributeError`.
- Fix: After session_tools operations complete, `_handle_session_management()`
  now loads the `AgentSession`, calls `_finalize_session(skip_summary=True)`,
  and builds a standard `group_args` result identical to a normal run's result.
  The GUI receives session history, cycle count, and summary with no special cases.
- `_finalize_session` gained a `skip_summary` kwarg (default `False`) that
  suppresses the `_generate_ai_summary()` LLM call — unnecessary since no new
  cycles ran during a display/remove operation.
- Files: `programs/ai_agent.py`

### P3: `get_results()` safe before `run()` (Fix 28)

**`AttributeError: 'Program' object has no attribute 'result'`**
- Any code path that called `get_results()` before `run()` completed — or on
  any early-exit path that bypassed result assignment — raised AttributeError.
- Fix: `run()` assigns `self.result = None` as its very first statement;
  `get_results()` uses `getattr(self, 'result', None)` as a defensive fallback.
- Files: `programs/ai_agent.py`

### P4: `restart_mode` auto-set on session management params (Fix 29)

**`display_and_stop` / `remove_last_n` required explicit `restart_mode=resume`**
- Both session management parameters operate on an existing session directory,
  which requires resume semantics. Forgetting `restart_mode=resume` silently
  cleared the session's log_directory instead of reusing it.
- Fix: At the start of `run()`, before `set_defaults()`, if either parameter
  is set, `restart_mode` is automatically forced to `'resume'`. The
  `display_and_stop` Phil choice default `'None'` (string) is handled
  correctly — the guard checks `!= 'None'` rather than Python truthiness.
- Files: `programs/ai_agent.py`

### Q1: Extending a completed workflow with new `project_advice` (Fix 30)

**Resuming after workflow completion with new advice was silently ignored**

**Scenario:** Agent finishes a ligand-protein complex (xtriage → phaser → refine ×3 →
ligandfit → pdbtools → molprobity). User resumes with
`project_advice="also run polder on the ligand chain B residue 100"`.
Previously the agent replied "workflow complete, nothing to do" and stopped
immediately — the new advice was never acted on.

**Root cause — two walls:**
1. **Wall 1 — AUTO-STOP in PLAN** (already fixed before this version):
   `metrics_trend.should_stop` would terminate before the LLM planned.
   The `advice_changed` flag (set by `_preprocess_user_advice` on hash
   mismatch) suppressed this for one cycle.
2. **Wall 2 — `valid_programs = ['STOP']` in PERCEIVE** (this fix):
   Once the workflow phase was `complete`, `detect_phase` returned the
   terminal phase, and `get_valid_programs` immediately returned `['STOP']`
   before Wall 1's suppression logic ran. The LLM was handed a program
   menu that only said STOP, so it couldn't choose polder even with
   Wall 1 down.

**Fix** (`agent/graph_nodes.py`, PERCEIVE node — 20 lines):
```python
# When advice_changed=True and phase='complete', step back to 'validate'
if (session_info.get("advice_changed") and
        workflow_state.get("phase_info", {}).get("phase") == "complete"):
    _new_valid = engine.get_valid_programs(exp, {"phase": "validate"}, ctx)
    workflow_state["valid_programs"] = _new_valid
    workflow_state["phase_info"] = {"phase": "validate", "reason": "advice_changed"}
```

The `validate` phase contains exactly the right program menu for
post-completion follow-up: `phenix.polder`, `phenix.molprobity`,
`phenix.model_vs_data`, `phenix.map_correlations`, plus `STOP` so the
LLM can still exit if the advice requires no action. After one successful
cycle, `advice_changed` is cleared and normal AUTO-STOP behaviour resumes.

**Complete event flow on resume with new advice:**
```
1. _preprocess_user_advice()
     new hash ≠ stored hash → session.data["advice_changed"] = True

2. PERCEIVE (Q1 fix)
     phase == 'complete' AND advice_changed
     → valid_programs = ['phenix.polder', 'phenix.molprobity',
                         'phenix.model_vs_data', 'phenix.map_correlations',
                         'STOP']
     → phase_info['phase'] set to 'validate'

3. PLAN (Wall 1 fix, pre-existing)
     metrics_trend.should_stop AND advice_changed
     → AUTO-STOP suppressed for this cycle

4. LLM
     sees new advice + validate-phase program menu → chooses phenix.polder

5. Post-cycle cleanup
     advice_changed = False
     → next cycle: normal termination logic resumes
```

**Notable behaviour:** `phenix.polder` intentionally lacks `strategy: run_once`
in `programs.yaml` — different residues and ligands may each need separate
omit maps. `polder_done=True` therefore does NOT block polder from reappearing
in `valid_programs`, allowing additional selections on subsequent resumes.

- Files: `agent/graph_nodes.py`

### Tests added (`tests/tst_audit_fixes.py`)

**P1/P3/P4 tests (20 tests)** covering session_tools functions, real method
extraction via `_build_agent_stub`, `get_results()` safety, and restart_mode
auto-set. See previous session-management context for full list.

**Q1 tests (9 tests):**
- `test_q1_complete_phase_has_only_stop`: baseline — `complete` phase → `['STOP']`
- `test_q1_validate_phase_includes_polder`: validate phase contains `phenix.polder`
- `test_q1_advice_changed_steps_back_to_validate`: core logic — step-back adds polder
- `test_q1_no_step_back_when_advice_unchanged`: unchanged advice keeps `complete` phase
- `test_q1_polder_reruns_allowed_when_already_done`: `polder_done=True` does NOT block re-run
- `test_q1_cryoem_complete_phase_steps_back`: logic is experiment-type agnostic
- `test_q1_graph_nodes_perceive_mutates_state`: end-to-end state mutation test
- `test_q1_advice_cleared_after_one_cycle`: `advice_changed` clears after one cycle
- `test_q1_step_back_does_not_apply_outside_complete`: guard only fires on `complete` phase

**Total tests: 74 (was 65 in v112_30)**

---

## Version 112.14 (Systematic Audit — Categories I, J, E, G, H)

### I1: max_refine_cycles produces bare STOP instead of controlled landing (Fix 21)

**`_apply_directives` returned `["STOP"]` when refinement limit reached**
- When `max_refine_cycles` was reached, the workflow engine stripped
  refinement programs and returned `["STOP"]`, terminating the workflow
  without validation. This left the user with no quality report.
- Fix: After removing refinement programs, `_apply_directives` now
  injects the validate-phase programs appropriate to the experiment type
  (`phenix.molprobity`, `phenix.model_vs_data`, `phenix.map_correlations`
  for X-ray; `phenix.molprobity`, `phenix.validation_cryoem`,
  `phenix.map_correlations` for cryo-EM), then appends STOP so the user
  can still exit immediately if desired.
- Also fixed: cryoem path was reading `context["refine_count"]` to check
  the limit, but cryo-EM refinement is counted in `context["rsr_count"]`.
  Now uses `rsr_count` for cryoem, `refine_count` for xray.
- Design note: `after_program` continues to produce STOP only (it is an
  explicit, unconditional stop). The validate-injection only applies to
  `max_refine_cycles` (a "limit" directive, not a "stop here" directive).
- Files: `agent/workflow_engine.py`

### J2: `_is_failed_result` false-positives on bare ERROR variants (Fix 22)

**Patterns `'ERROR '`, `': ERROR'`, `'ERROR:'` matched non-fatal log text**
- Phenix logs routinely contain strings like "Error model parameter",
  "Expected errors: 0", "No ERROR detected". All of these matched the
  broad `ERROR ` / `ERROR:` / `: ERROR` patterns, causing legitimate
  runs to be classified as failed and their done flags suppressed.
- Priority order for failure detection (per spec J2):
  1. Exit code (handled at the shell layer, before `_is_failed_result`)
  2. Output file check
  3. Log text with specific Phenix terminal phrases
- Fix: Removed the three generic `ERROR` patterns. Retained the seven
  Phenix-specific terminal failure signatures: `FAILED`, `SORRY:`,
  `SORRY `, `*** ERROR`, `FATAL:`, `TRACEBACK`, `EXCEPTION`.
  These cover all real Phenix failure modes without matching non-fatal text.
- Files: `agent/workflow_state.py`

### J5: Zombie state detection — stale done flags block re-execution (Fix 23)

**Missing output files left done flags True, preventing re-run**
- When the agent crashed mid-cycle or the user deleted output files, the
  history record retained `done_flag=True`. The phase detector saw
  `done=True` and skipped the program, but file-based flags
  (`has_full_map`, `has_placed_model`) were False because no file was
  found. The workflow became stuck.
- Fix: `_clear_zombie_done_flags(history_info, available_files)` checks
  each done flag against its expected output file pattern. If the flag
  is True but no matching file exists in `available_files`, it clears
  the flag (and associated file flags) in-memory without modifying
  history. Diagnostic messages are emitted to PERCEIVE.
- Programs covered:
  - `resolve_cryo_em_done` → `denmod_map.ccp4` → also clears `has_full_map`
  - `predict_full_done` → `*_overall_best.pdb` → also clears `has_placed_model`
  - `dock_done` → `*_docked.pdb` → also clears `has_placed_model`
  - `refine_done` → `*_refine_001.pdb` → decrements `refine_count`
  - `rsr_done` → `*_real_space_refined*.pdb` → decrements `rsr_count`
- Files: `agent/workflow_state.py`

### E1: xtriage resolution regex extracts 50.0 instead of 2.3 (Fix 24)

**Dash separator in "50.00 - 2.30" format not handled by skip group**
- The pattern `(?:[0-9.]+\s+)?` was designed to skip the low-resolution
  limit and capture the high-resolution limit. But xtriage formats the
  range as `50.00 - 2.30` (space-dash-space), not `50.00 2.30` (space
  only). With the dash present the optional skip group backtracks, and
  the capture group then matches `50.00` — the wrong value. `pick_min`
  across all lines then returned 50.0 instead of 2.3.
- Fix: Changed `(?:[0-9.]+\s+)?` to `(?:[0-9.]+\s*[-]\s*)?`. The skip
  group now explicitly handles the dash separator. The pattern still
  handles `Resolution: 1.80` (no range, skip group does not fire).
  The negative lookbehind anchors prevent "Completeness in resolution
  range: 1" from matching (J2-era fix retained).
- Files: `knowledge/programs.yaml` (`phenix.xtriage` `resolution` pattern)

### E1: real_space_refine map_cc returns first cycle instead of final (Fix 25)

**`extract: first` (default) captured the initial, worst map_cc value**
- RSR emits one `CC_mask =` line per macro-cycle. With `extract: first`
  the initial (lowest) map_cc was reported rather than the final (best).
  This caused the agent to incorrectly judge model quality as poor and
  over-refine.
- Also fixed: the pattern `CC_mask\s*[=:]\s*([0-9.]+)` was narrower
  than the standard map_cc pattern used by all other programs. Broadened
  to `(?:CC_mask|Map-CC|Model vs map CC)\s*[:=]?\s*([0-9.]+)` for
  consistency.
- Files: `knowledge/programs.yaml` (`phenix.real_space_refine` `map_cc` spec)

### G1: holton_geometry_validation defined but not in any workflow phase (Verified)

`phenix.holton_geometry_validation` is registered in `programs.yaml`
with `done_tracking` but deliberately not in any `workflows.yaml` phase.
Added an audit comment to clarify the intentional status and document
how to activate it (add to the validate phase in `workflows.yaml`).
Files: `knowledge/programs.yaml`

### H1/H3: STATE_MAP comments clarified (Documentation)

- `cryoem_has_model`: The state string assigned to `check_map` and
  `optimize_map` is a legacy misnomer — no model exists at that point.
  No behavioral code gates on this string; `phase_info["phase"]` is
  used for all internal decisions. Comment added.
- `validate` shares state string with `refine` (`xray_refined` /
  `cryoem_refined`): This is intentional for external API compatibility.
  Internal code always uses `phase_info["phase"]` to distinguish them.
  Comment added.
- Files: `agent/workflow_engine.py`

### J5d: Zombie state diagnostics silently discarded (Fix 26)

**`_clear_zombie_done_flags()` return value was ignored**
- `detect_workflow_state()` called `_clear_zombie_done_flags()` but did not
  capture the return value. The done-flag clearing worked correctly (in-memory
  modification), but the diagnostic messages explaining *why* a previously
  "done" program reappeared in `valid_programs` were silently lost.
- Without these diagnostics, users and developers had no way to know that a
  zombie state was detected and resolved — making crash/restart scenarios very
  confusing to debug.
- Fix: captured the return value as `zombie_diagnostics` and added it to the
  state dict under key `"zombie_diagnostics"` (present only when non-empty).
  The PERCEIVE node now logs each diagnostic prefixed with
  `"PERCEIVE: ZOMBIE STATE — "`.
- Files: `agent/workflow_state.py`, `agent/graph_nodes.py`

### Regression tests added: `tests/tst_audit_fixes.py`

27 tests covering all bugs found and fixed in this audit session:

- **J2** (3 tests): `_is_failed_result` true positives, false positive
  elimination, done-flag blocking on failure
- **J5** (5 tests): Zombie detection for resolve_cryo_em, refine, dock;
  preservation when output exists; rsr_count decrement
- **E1/E2** (6 tests): xtriage resolution dash-separator, completeness
  anchor, simple format, multiple ranges; RSR map_cc last-cycle and
  pattern variants
- **I1** (2 tests): X-ray controlled landing (validate + STOP injected);
  cryo-EM rsr_count used (not refine_count)
- **I2** (1 test): after_program → STOP only (no validate injection)
- **I1b** (2 tests): xray and cryoem complete phase → [STOP] after validation_done=True (clean-termination path)
- **J5 zombie surfacing** (2 tests): `zombie_diagnostics` present in state when zombie cleared; absent for clean state
- **YAML spec** (4 tests): xtriage pick_min, RSR map_cc extract:last,
  RSR clashscore extract:last, polder requires_selection invariant

Tests registered in `tests/run_all_tests.py` as "Audit Fix Regressions".



### Explicit program injection overrides multi-step directives (Fix 19)

**`_detect_explicit_program_request` hijacks multi-step workflows**
- When user says "run mtriage, resolve_cryo_em, map_symmetry", the
  directive system correctly parses ordering (start_with_program,
  after_program). But `_detect_explicit_program_request` scans the raw
  text, returns ONE program (e.g., map_symmetry), and injects
  "**IMPORTANT: run phenix.map_symmetry**" into guidelines — overriding
  the directive ordering. The LLM then picks map_symmetry instead of
  resolve_cryo_em (the correct next step).
- This cascades: map_symmetry fails (no full_map, only half maps),
  then resolve_cryo_em runs as a fallback without LLM-guided parameters.
- Fix: Skip explicit program injection when directives contain
  multi-step workflow info (`start_with_program`, `after_program`, or
  `prefer_programs` with 2+ entries). Applied to both the main
  `_query_agent_for_command` path and the retry/duplicate handler.
- Files: `programs/ai_agent.py`

### resolve_cryo_em missing input_priorities (Fix 19b)

**Category-based file selection for half_map slot**
- `phenix.resolve_cryo_em` had empty `input_priorities`, forcing the
  command builder to use extension-only matching for the `half_map`
  slot. This could select sharpened or optimized maps instead of actual
  half maps when multiple .ccp4 files are available.
- Added `input_priorities.half_map` with `categories: [half_map]` and
  `exclude_categories: [full_map, optimized_full_map]`.
- Files: `knowledge/programs.yaml`

### run_once strategy check fix (Fix 17)

**Three broken `tracking.get("run_once")` checks in workflow_engine.py**
- The YAML uses `strategy: "run_once"` but all three checks in
  `workflow_engine.py` used `tracking.get("run_once")` which always
  returned None. This meant run_once programs (map_symmetry, mtriage,
  xtriage) were never filtered from valid_programs.
- Fixed all three locations to check
  `tracking.get("strategy") == "run_once" or tracking.get("run_once")`
  for backward compatibility.
- Additionally, `_apply_directives` could re-add already-done programs
  via `start_with_program`, `program_settings`, and `after_program`
  directives. Added `_is_program_already_done()` helper that checks
  the YAML done_tracking config. All three directive paths now skip
  programs whose done flag is already set.
- Files: `agent/workflow_engine.py`

### Sharpened map mis-categorized as half_map (Fix 18)

**half_map excludes for sharpened/optimized maps**
- `emd_XXXXX_half_map_N_box_sharpened.ccp4` (output of map_sharpening)
  matched `half_map` pattern `*half*`, causing mtriage to use it as
  `half_map=` instead of `full_map=`.
- Added excludes `*sharpened*`, `*sharpen*`, `*optimized*` to the
  `half_map` category in `file_categories.yaml`. Sharpened maps now
  fall through to `full_map` via the map parent category, where
  mtriage's `input_priorities.full_map.categories: [full_map, map]`
  picks them up correctly.
- Files: `knowledge/file_categories.yaml`

### forced_program not enforced in plan node (Fix 19)

**LLM could override multi-step workflow ordering**
- `workflow_engine` sets `forced_program` from `after_program` directive
  (e.g., for "run mtriage, resolve_cryo_em, and map_symmetry"), but the
  plan node never enforced it. The LLM could freely pick any valid
  program, ignoring the directive ordering.
- Added forced_program enforcement block in plan node: when
  `forced_program` is set and valid, the LLM's choice is overridden.
- Also fixed explicit_program injection in perceive: when
  `forced_program` is set from directives, `explicit_program` (from
  `_detect_explicit_program_request` text scanning) is no longer
  injected into valid_programs, preventing conflicting LLM hints.
- Files: `agent/graph_nodes.py`

### map_symmetry offered without full map (Fix 20)

**map_symmetry should not be valid when only half maps exist**
- map_symmetry's input_priorities exclude half_map from the map slot,
  so it always fails to build when only half maps are available. But
  its workflow condition only checked `not_done: map_symmetry`, so it
  appeared in valid_programs and the LLM would pick it.
- Added `has: non_half_map` condition to map_symmetry in workflows.yaml.
- Added composite context key `has_non_half_map` in workflow_engine.py
  that checks `set(map files) - set(half_map files)`. This correctly
  becomes True after map_sharpening produces a sharpened map (which is
  in `map` but not `half_map`), even though the sharpened filename
  contains "half" and doesn't match `full_map`.
- Files: `agent/workflow_engine.py`, `knowledge/workflows.yaml`

### File discovery and filtering (Fixes 14-15)

**Companion file discovery** (REMOVED in v112.79)
*`graph_nodes._discover_companion_files` was removed because it
scanned directories of user-supplied files, picking up unintended
files.  All companion discovery is now handled by the session layer
(`_find_missing_outputs` and `get_available_files` Step 3 scan),
which is scoped to agent output directories only.*
- After `phenix.refine`: discovers map coefficients (`refine_NNN.mtz`) and
  refined model (`refine_NNN.pdb`) from `_data.mtz` prefix. Handles both
  bare (`refine_001.mtz`) and `_001` (`refine_001_001.mtz`) naming.
- After `phenix.autobuild`: discovers `overall_best.pdb` when only
  `overall_best_refine_data.mtz` is tracked by the client.
- After `phenix.pdbtools`: scans sibling `sub_*_pdbtools/` directories in
  the agent directory for `*_with_ligand.pdb` output files.
- Defense-in-depth: also added `session._find_missing_outputs()` as a
  second discovery layer for session-tracked output files.

**Intermediate file filtering** (`graph_nodes._filter_intermediate_files`)
- Filters files from ligandfit's internal `TEMP0/` directories and files
  with `EDITED_` or `TEMP_` prefixes before categorization.
- Prevents intermediate files from being selected as model inputs.
- Also added `EDITED*` and `TEMP*` exclusions to `unclassified_pdb`
  category in `file_categories.yaml`.
- Files: `agent/graph_nodes.py`, `agent/session.py`,
  `knowledge/file_categories.yaml`

### Refinement loop enforcement (Fix 13)

**At-target actively removes refine from valid programs**
- When `_is_at_target` returns True (hopeless R-free > 0.50 after 1+ cycles,
  or hard limit of 3+ cycles), `phenix.refine` and `phenix.real_space_refine`
  are now explicitly **removed** from valid programs in both `validate` and
  `refine` phases. Previously only prevented adding as supplement.
- `STOP` added to valid programs when at target.
- Exception: `needs_post_ligandfit_refine` always allows refinement (model
  changed after ligand fitting, re-refinement is scientifically required).
- Files: `agent/workflow_engine.py`

### Command validation fix (Fix 15)

**Output file arguments excluded from input validation**
- `output.file_name=X.pdb`, `output.prefix=Y`, etc. are now stripped from
  commands before extracting file references for validation. Previously
  `output.file_name=model_with_ligand.pdb` was treated as an input file
  reference and rejected as "not found in available_files".
- Files: `agent/graph_nodes.py`

### best_files excluded category check (Fix 16)

**Prevents ligand fragments from being used as refine model input**
- `best_files["model"]` is now checked against the program's
  `exclude_categories` before being applied as a model override. If the
  best model (e.g., `ligand_fit_1.pdb`) is in an excluded category
  (e.g., `ligand_fit_output` → `ligand`), it is skipped with a log message.
- Applied to both pre-population and LLM override paths in
  `command_builder.py`.
- Files: `agent/command_builder.py`

### MR-SAD phaser condition (Fix 11)

**Composite `has_model_for_mr` context key**
- Added `has_model_for_mr` in `workflow_engine.py` that checks both `model`
  and `search_model` file categories. Phaser condition in `workflows.yaml`
  changed from `has: model` to `has: model_for_mr`.
- Ensures phaser is available when user provides a dedicated search model
  (`search_model.pdb`) that categorizes as `search_model`, not `model`.
- Files: `agent/workflow_engine.py`, `knowledge/workflows.yaml`

### Cryo-EM experiment type inference (Fix 12)

**Advice preprocessor now infers experiment type from file extensions**
- Added experiment type inference rules to the advice preprocessing LLM
  prompt: `.mtz/.sca/.hkl` → X-ray, `.map/.mrc/.ccp4` → cryo-EM,
  half-maps → cryo-EM, `.pdb + .map` → cryo-EM refinement.
- Cosmetic fix only: actual workflow engine already correctly detected
  cryo-EM from file categories. This fixes the user-facing advice text.
- Files: `agent/advice_preprocessor.py`

### Input priority improvements

- `phenix.pdbtools` protein: added `autobuild_output` to
  `prefer_subcategories`, `EDITED` to `exclude_patterns`,
  `overall_best` to `priority_patterns`.
- `phenix.refine` model: `with_ligand` is first in `prefer_subcategories`,
  ensuring the combined protein+ligand model is selected for post-pdbtools
  refinement.
- Files: `knowledge/programs.yaml`

### Test coverage

- 30 tests in `tests/tst_v112_13_fixes.py` covering perceive
  pipeline safety (2, v112.79), intermediate filtering (3),
  file categorization (5), phaser model_for_mr (3), output
  validation (3), program priorities (4), end-to-end
  post-pdbtools selection (2), combine_ligand phase (1),
  sharpened map categorization (2), run_once done-flag
  config (2), map_symmetry condition (1),
  non_half_map context key (2).
- Total: 29/35 passing (6 pre-existing libtbx import failures).

## Version 112.12 (February 2025)

### Done-tracking strategy enum (Fix 10C)

**Replaced `run_once: true` with `strategy` enum; unified detection to one system**

- Added `strategy: "set_flag" | "run_once" | "count"` to done_tracking in
  programs.yaml. `set_flag` is the default (simple done flag), `run_once`
  replaces the old boolean, `count` handles programs that need run counting.
- Added `count_field` and `exclude_markers` to history_detection schema.
  `count_field` specifies the counter name (e.g., "refine_count");
  `exclude_markers` rejects matches (checked BEFORE markers).
- Moved 4 remaining Python-only blocks to YAML: validation (4 programs
  share one flag via markers), phaser (count), refine (count + exclude),
  real_space_refine (count). Only predict_and_build cascade stays in Python.
- Added `ALLOWED_COUNT_FIELDS` whitelist validated at load time — prevents
  typos in YAML count_field from silently creating garbage attributes.
- Removed 3 dead success flags: `refine_success`, `rsr_success`,
  `phaser_success` (set but never read anywhere).
- Unified `_set_simple_done_flags()` → `_set_done_flags()` handling all
  strategies. Removed redundant `detect_programs_in_history()` calls.
- Counter fields initialized dynamically from YAML configs (single source
  of truth) instead of hardcoded in info dict.
- 36 conformance tests passing (was 33): added test_strategy_enum_values,
  test_count_field_validation, test_exclude_markers_prevent_false_matches.
- Files: `knowledge/programs.yaml`, `agent/workflow_state.py`,
  `knowledge/program_registration.py`, `tests/tst_hardcoded_cleanup.py`

## Version 112.11 (February 2025)

### Phase 3: Final hardcoded cleanup (Fixes 6, 8, 10B)

#### Stop-directive patterns → YAML (Fix 8)

**Moved 18 regex patterns from `directive_extractor.py` to `programs.yaml`**
- 9 programs now have `stop_directive_patterns` in YAML (phenix.mtriage,
  xtriage, phaser, ligandfit, refine, autobuild, map_to_model, dock_in_map,
  map_symmetry)
- `_get_stop_directive_patterns()` loads from YAML with length-based sorting
  (longest patterns match first — handles map_to_model vs dock_in_map safely)
- Density modification branching stays in Python (requires experiment-type context)
- Hardcoded fallback with DeprecationWarning if YAML unavailable
- Files: `knowledge/programs.yaml`, `agent/directive_extractor.py`

#### Rules-selector priority lists → YAML (Fix 6)

**Moved 5 priority lists from `rules_selector.py` to `workflows.yaml`**
- `shared/rules_config` section has `default_priority`, `ligand_priority`,
  and `state_aliases`
- Per-phase `rules_priority` lists in workflow phase definitions
- `_load_rules_config()` and `_get_phase_rules_priority()` load from YAML
- r_free/map_cc validation logic stays in Python (behavioral, not config)
- Files: `knowledge/workflows.yaml`, `agent/rules_selector.py`

#### Simple done-flag detection → YAML (Fix 10B)

**Moved 10 simple if/elif blocks from `_analyze_history()` to YAML**
- 10 programs now have `history_detection` in `done_tracking` with markers,
  alt_markers/alt_requires (AND logic), and optional success_flag
- Programs: process_predicted_model, autobuild, autobuild_denmod, autosol,
  ligandfit, pdbtools, dock_in_map, map_to_model, resolve_cryo_em, map_sharpening
- New `_set_simple_done_flags()` replaces ~40 lines of if/elif blocks
- Complex cases stay in Python: validation (3 programs share flag), phaser
  (count + TFZ check), predict_and_build (cascade), refine/rsr (counts)
- Fixed process_predicted_model flag name: YAML now uses `process_predicted_done`
  to match actual usage in workflow_engine.py
- Added design note in programs.yaml: `strategy: "run_once" | "count" | "success_gate"`
  is a cleaner generalization of the current `run_once: true` boolean — noted as
  future consideration, not required now since Fix 10A already nests run_once
  correctly inside done_tracking
- Files: `knowledge/programs.yaml`, `agent/workflow_state.py`

### Test coverage

- 33 conformance tests in tst_hardcoded_cleanup.py (was 30)
- New tests: test_history_detection_coverage, test_history_detection_behavioral,
  test_no_simple_done_flags_in_analyze_history

## Version 112.10 (February 2025)

### Dead code removal in planner.py (Fix 7)

**Removed ~1300 lines of dead code from `agent/planner.py`**
- Investigation confirmed `generate_next_move()`, `construct_command_mechanically()`,
  `get_required_params()`, `extract_clean_command()`, `get_relative_path()`,
  `get_program_keywords()`, and `fix_multiword_parameters()` are never called externally
- All were superseded by the YAML-driven CommandBuilder + rules_selector pipeline
- Retained only `fix_program_parameters()` (called by graph_nodes.py) and
  `extract_output_files()` (called by run_ai_analysis.py)
- Removed heavy imports: langchain_core, phenix_knowledge, validation, memory
- File: `agent/planner.py` (1436 → ~130 lines)

### GUI app_id fallback in programs.yaml (Fix 9)

**Added `gui_app_id` fields to programs.yaml as fallback for headless environments**
- 20 programs now have `gui_app_id` in YAML (3 without GUI windows excluded)
- 3 programs have `gui_app_id_cryoem` for cryo-EM variant windows:
  predict_and_build, map_correlations, map_sharpening
- `_build_program_to_app_id()` falls back to YAML when GUI PHIL is unavailable
- `get_app_id_for_program()` cryo variants now loaded from YAML (with hardcoded baseline)
- Added `gui_app_id`, `gui_app_id_cryoem` to yaml_tools valid fields
- Files: `knowledge/programs.yaml`, `programs/ai_agent.py`, `agent/yaml_tools.py`

### Unified done flag tracking in programs.yaml (Fix 10 Step A)

**Eliminated `_MANUAL_DONE_FLAGS` dict and top-level `run_once` field**
- All done flags now defined in `done_tracking` blocks in programs.yaml
- `run_once` moved from top-level field into `done_tracking.run_once`
- `get_program_done_flag_map()` reads directly from YAML — no hardcoded dict
- Added `done_tracking` to 4 previously missing programs (holton_geometry_validation,
  model_vs_data, validation_cryoem, map_correlations removed stale `run_once: false`)
- `workflow_engine.py` now reads `done_tracking.flag` instead of deriving flag names
- Updated ADDING_PROGRAMS.md, ARCHITECTURE.md, OVERVIEW.md docs
- Files: `knowledge/programs.yaml`, `knowledge/program_registration.py`,
  `agent/workflow_engine.py`, `agent/yaml_tools.py`

## Version 112.8 (February 2025)

### Prerequisite mechanism for resolution-dependent programs

**Programs that need resolution (RSR, dock_in_map, map_to_model) now auto-trigger mtriage**
- When a program's resolution invariant can't be satisfied, the command builder
  detects a `prerequisite: phenix.mtriage` declaration and automatically builds
  the mtriage command instead
- Next cycle: resolution is available from mtriage output → original program
  builds successfully
- Respects `skip_programs`: if user skipped mtriage, returns clear error instead
  of silently failing
- Files: `knowledge/programs.yaml` (prerequisite declarations),
  `agent/command_builder.py` (prerequisite tracking),
  `agent/graph_nodes.py` (prerequisite build logic)

### LLM resolution hallucination guard

**Strip unverified resolution values from LLM strategy**
- LLMs frequently hallucinate resolution values (e.g., "resolution": 3.1)
  which bypassed the prerequisite mechanism
- `_build_strategy()` now checks whether resolution came from the LLM AND
  whether `context.resolution` (verified source) is None — if so, strips it
- Log: `BUILD: Stripped LLM-hallucinated resolution=3.1 (no verified source)`
- File: `agent/command_builder.py`

### Fixed false exclude_pattern matches in file selection

**`_2` pattern in exclude_patterns matched PDB codes like `_23883`**
- Changed mtriage and RSR exclude_patterns from `[half, _1, _2, _a, _b]` to
  `[half_1, half_2, half1, half2, _half]`
- Added `input_priorities` to mtriage for category-based file selection
  (bypasses extension fallback entirely)
- File: `knowledge/programs.yaml`

### Fixed half-map misuse as full_map in mtriage

**Two half maps → one used as full_map → wrong resolution**
- Root cause: half-maps bubbled up to `map` parent category, then selected for
  `full_map` slot via category fallback. Post-selection validation then deleted
  the legitimate half_maps as "redundant"
- Fix 1: Added `exclude_categories: [half_map]` to mtriage's full_map priorities
- Fix 2: Post-selection validation now checks if the "full_map" is actually a
  categorized half-map — if so, removes the mis-selected full_map instead
- Files: `knowledge/programs.yaml`, `agent/command_builder.py`

### Fixed autosol atom_type crash from multi-atom values

**`atom_type="Se, S"` → bare `S` on command line → crash**
- LLMs put multiple atom types in `atom_type` field (e.g., "Se, S") even when
  `additional_atom_types` is correctly set separately
- `_build_strategy()` now sanitizes: splits on comma/space, keeps first atom in
  `atom_type`, moves extras to `additional_atom_types` (if not already set)
- File: `agent/command_builder.py`

### Fixed predict_and_build intermediate file tracking

**Internal `working_model_full_docked.pdb` was tracked as valid output**
- Root cause: `docked` in filename matched `valuable_output_patterns` which
  overrode all intermediate exclusions
- Removed overly broad `(\S*docked\S*\.pdb)` from log parser known_patterns
- Added exclusions for `/local_dock_and_rebuild`, `/local_rebuild` paths
- Added `intermediate_basename_prefixes` for `working_model` (always excluded,
  even if matching a "valuable" pattern)
- Added `working_model*` to `docked` category excludes in file_categories.yaml
- Files: `phenix_ai/log_parsers.py`, `phenix_ai/utilities.py`,
  `programs/ai_agent.py`, `knowledge/file_categories.yaml`

### Fixed .eff file generation for old-style programs

**phenix.refine .eff had `generate=False` despite command saying `generate=True`**
- Root cause: `master_phil.fetch()` requires exact scope paths, but agent uses
  short-form like `xray_data.r_free_flags.generate=True` (full path is
  `refinement.input.xray_data.r_free_flags.generate`)
- Fix: Use `master_phil.command_line_argument_interpreter()` to resolve short
  paths before `fetch()` — same mechanism phenix.refine's own CLI uses
- File: `programs/ai_agent.py`

### Fixed skip_programs causing workflow deadlock

**Skipping xtriage → stuck in "analyze" phase → STOP**
- Root cause: Phase detection checked `xtriage_done` before allowing progression;
  `_apply_directives` removed xtriage from programs but couldn't change the phase
- Fix: Skipped programs are treated as "done" in `build_context()` — their done
  flags are set before phase detection runs
- Done flag mapping now auto-generated via `get_program_done_flag_map()` in
  `program_registration.py` (combines run_once auto-flags with manual mappings)
- Files: `agent/workflow_engine.py`, `knowledge/program_registration.py`

### RSR GUI reload crash fix

**TypeError on `get_output_dir()` after successful RSR execution**
- Status mismatch: native execution returns "complete" but guard checked for
  "success" / "completed" — pkl_path never sent to GUI
- Added "complete" to status check, plus pkl validation before sending
- File: `programs/ai_agent.py`

### Bare command rejection

**`phenix.mtriage` with no arguments hung waiting for input**
- Added explicit bare command check after assembly: commands with fewer than
  2 parts (just the program name) are rejected
- File: `agent/command_builder.py`

## Version 112.3 (February 2025)

### Removed langchain-classic dependency

**Direct implementation replaces deprecated langchain chains/retrievers**
- Replaced `create_stuff_documents_chain` (from `langchain.chains`) with 4-line
  direct implementation: concatenate docs → format prompt → call LLM
- Replaced `ContextualCompressionRetriever` (from `langchain.retrievers`) with
  minimal `_CompressionRetriever(BaseRetriever)` class using `langchain_core`
- Zero `from langchain.` imports remain; all code uses `langchain_core`,
  `langchain_community`, and provider-specific packages only
- Files: `analysis/summarizer.py`, `rag/retriever.py`, `docs/README.md`

### Added phenix.map_correlations support

**New program in YAML registry with multi-mode input support**
- Supports 5 input modes: model+map, model+mtz, map+map, mtz+mtz, map+mtz
- Uses flag-based file assignment (`input_files.model=`, `input_files.map_in_1=`, etc.)
- `map2` and `map_coeffs_2` slots set `auto_fill: false` to prevent the command
  builder from duplicating the same file into both map slots
- Log parsing extracts: `cc_mask`, `cc_volume`, `cc_peaks`, `cc_box`, `map_map_cc`
- Added to both xray and cryoem `validate` phases in `workflows.yaml`
- Added step_metrics entry (`CC_mask: {cc_mask:.3f}`) in `metrics.yaml`
- Added quality_table row with CC_volume detail and assessment in `metrics.yaml`
- Added `cc_mask_assessment` using same thresholds as `map_cc_assessment`
- Files: `knowledge/programs.yaml`, `knowledge/workflows.yaml`,
  `knowledge/metrics.yaml`, `agent/session.py`

### Explicit program request handling

**Hard stop for unregistered `phenix.X` requests, graceful fallback for bare names**
- When user writes `phenix.some_program` explicitly and the program is not in
  the YAML registry, raise Sorry with a clear message
- When a bare name match (e.g., "anomalous signal" matching
  `phenix.anomalous_signal`) refers to an unregistered program, silently ignore
  the match — it's likely a false positive from natural language, not a deliberate
  program request. The agent proceeds with its normal workflow.
- File: `programs/ai_agent.py`

**Explicit program injection into valid_programs**
- Registered explicit programs are injected into `valid_programs` regardless of
  workflow phase, so the LLM can select them even in early phases (e.g.,
  `map_correlations` during `cryoem_initial`)
- File: `agent/graph_nodes.py`

**STOP override for unfulfilled explicit requests**
- If the LLM chooses STOP but the user's explicitly requested program hasn't
  run yet, the plan step overrides STOP and forces the explicit program
- Checks `session_info["explicit_program"]` against history of programs run
- File: `agent/graph_nodes.py`

### Transport pipeline: explicit_program passthrough

**Added `explicit_program` to three transport whitelists**
- `build_session_state()` in `agent/api_client.py`
- `build_request_v2()` normalization in `agent/api_client.py`
- `session_state → session_info` mapping in `phenix_ai/run_ai_agent.py`
- Without these, `explicit_program` was silently dropped during transport
  from client to server, preventing injection and STOP override from working

### Program name resolution fixes

**Bare name ↔ `phenix.` prefix lookup fallback**
- `get_program()` in `yaml_loader.py` now tries `phenix.` + bare_name if
  the initial lookup returns None
- `_resolve_program_patterns()` helper in `metric_patterns.py` does the same
  for metric pattern lookups (used by all 3 lookup functions)
- Root cause: `ai_agent.py` strips the prefix at line 2077
  (`command.split()[0].replace("phenix.", "")`) before passing to metric
  extraction, so lookups against YAML keys (which use full names) failed
- Files: `knowledge/yaml_loader.py`, `knowledge/metric_patterns.py`

### Metric display pattern fixes

**Regex patterns match both raw log and reformatted report formats**
- `format_metrics_report()` transforms `cc_mask` → `Cc Mask` via
  `.replace("_", " ").title()`, so result text stored in cycles uses
  the reformatted form
- Updated all CC patterns in `_extract_final_metrics()` to use
  `CC[_ ]?mask` with `re.IGNORECASE` to match both formats
- Updated `map_cc` pattern in `_extract_metrics_from_result()` similarly
- File: `agent/session.py`

### Minor fixes

**Reasoning truncation increased**
- Main reasoning: 500 → 1000 chars
- Session summary reasoning: 300 → 600 chars
- Files: `programs/ai_agent.py`, `agent/session.py`

**False input_directory warning fixed**
- When user supplies `original_files` directly, no longer warns
  "No input_directory to look for files" if files already present
- File: `programs/ai_agent.py`

**Unused import removed**
- Removed `from langchain_core.documents import Document as LCDocument`
  from `rag/retriever.py`

## Version 112.2 (February 2025)

### Cohere → FlashRank Migration

**Dependency removal - Replace Cohere API with local FlashRank reranker**
- Replaced `CohereRerank` (cloud API) with `FlashrankRerank` (local cross-encoder)
- Model: `ms-marco-MiniLM-L-12-v2` (~34MB, runs on CPU, no API key needed)
- Same `ContextualCompressionRetriever` pattern — callers unchanged
- Removed `COHERE_API_KEY` environment variable requirement
- Removed `CohereApiError` exception handling (local inference has no API errors)
- Updated privacy disclaimers (Cohere no longer contacted)
- Files: `rag/retriever.py`, `analysis/analyzer.py`, `utils/run_utils.py`,
  `programs/ai_agent.py`, `programs/ai_analysis.py`
- Install: `pip install flashrank` (replaces `cohere` + `langchain-cohere`)

**Script cleanup**
- `run_inspect_db.py`: Removed debug print and hardcoded scratch path
- `run_query_docs.py`: Replaced duplicated API key validation with shared
  `validate_api_keys()` from `utils/run_utils.py`

**Usability - Early exit when no inputs provided**
- Agent now stops with a helpful message if launched with no original_files,
  no project_advice, and no README in the input_directory
- Skipped when resuming an existing session (has previous cycles)
- File: `programs/ai_agent.py`

**Fix - Quote multi-word PHIL parameter values in commands**
- `unit_cell=114 114 32.5 90 90 90` was passed unquoted, causing shell/PHIL
  to split the values into separate arguments → crash
- Added `fix_multiword_parameters()` in `agent/planner.py` (LLM command path)
  and inline regex in `agent/program_registry.py` (registry command path)
- Now produces `unit_cell="114 114 32.5 90 90 90"` and `space_group="P 2 21 21"`
- Also handles prefixed forms like `xray_data.unit_cell=...`

**Fix - Strategy passthrough dropping known PHIL short names (HIGH IMPACT)**
- `unit_cell`, `space_group`, etc. from directives were silently dropped
  in `program_registry.py` because the passthrough required a dot in the key
- Server LLM returned `unit_cell=...` (no prefix) vs local returning
  `xray_data.unit_cell=...`, so the bug only manifested on the server
- Added `KNOWN_PHIL_SHORT_NAMES` set to the passthrough check so common
  parameters like `unit_cell`, `space_group`, `resolution`, `nproc`,
  `ncopies`, `twin_law` are accepted without their full PHIL scope prefix
- File: `agent/program_registry.py`

**Fix - Polder has no resolution keyword**
- Polder has no resolution keyword of any kind
- Removed `high_resolution` and `resolution` strategy flags and the
  `auto_fill_rfree_resolution` fix from polder's YAML
- Removed `resolution` (and `high_resolution`, `low_resolution`) from
  `KNOWN_PHIL_SHORT_NAMES` passthrough set in `program_registry.py` —
  these should only go through explicit strategy_flags, not blindly
  passed through to programs that don't support them
- Added polder entry to `parameter_fixes.json` to strip resolution/
  high_resolution/low_resolution/d_min as a safety net
- Added `has_strategy: selection` invariant to block polder if no
  selection is provided, forcing LLM retry
- Added hints telling LLM that selection is required and no resolution
  exists
- Files: `knowledge/programs.yaml`, `agent/program_registry.py`,
  `knowledge/parameter_fixes.json`

**Fix - Anomalous Resolution incorrectly reported as Resolution**
- In xtriage logs with anomalous data, "Anomalous Resolution: 9.80" was
  matched by the generic resolution regex before "Resolution: 2.50"
- Added negative lookbehind `(?<!nomalous )` to all resolution patterns:
  `agent/session.py` (hardcoded fallback), `agent/session_tools.py`,
  `knowledge/patterns.yaml` (centralized pattern)
- Updated xtriage YAML `log_parsing` to use `extract: last` so the summary
  "Resolution: 2.50" at the end of the log is preferred over earlier matches
- Added explicit `anomalous_resolution` pattern to xtriage YAML log_parsing
- Files: `agent/session.py`, `agent/session_tools.py`,
  `knowledge/programs.yaml`, `knowledge/patterns.yaml`

**Fix - AutoSol using Phaser output data instead of original data**
- After Phaser MR, `best_files['data_mtz']` was set to PHASER.1.mtz which
  has lost anomalous signal — useless for SAD phasing
- Added `input_priorities` for `data_mtz` in autosol YAML:
  `categories: [original_data_mtz, data_mtz]`,
  `exclude_categories: [phased_data_mtz]`,
  `skip_rfree_lock: true`
- AutoSol doesn't need rfree handling — added `skip_rfree_lock` support
  to `command_builder.py` PRIORITY 1 so autosol gets original data
- CRITICAL: LLM file choices were bypassing input_priorities entirely.
  Added `exclude_categories` check to the LLM hint acceptance path in
  `_select_files()` — now rejects LLM-chosen files in excluded categories
  (with both full-path and basename matching for robustness)
- Added prompt guidance telling LLM that autosol needs original data
- Files: `knowledge/programs.yaml`, `agent/command_builder.py`,
  `knowledge/prompts_hybrid.py`

**Fix - Wavelength mistakenly extracted as resolution in directives**
- Text like "data collected far from the iron edge (1.1158 Å)" was being
  extracted as `resolution=1.1158` by the directive extractor LLM
- Three fixes applied:
  1. LLM prompt: Added explicit "Do NOT confuse wavelength with resolution"
     warning with guidelines for distinguishing them
  2. Validation: Added post-extraction check that removes resolution values
     < 1.2 Å (wavelength range) or matching any extracted wavelength value
  3. Simple extractor: Removed overly broad patterns like "X Å" standalone;
     added wavelength cross-check; raised minimum from 0.5 to 1.0
- File: `agent/directive_extractor.py`

**Fix - AutoSol getting obs_labels from recovery strategy**
- Recovery strategy from xtriage was applying `autosol.input.xray_data.obs_labels=I(+)`
  to autosol commands — autosol handles labels internally
- Removed autosol from `DATA_LABEL_PARAMETERS` in command_builder.py
- Changed fallback behavior: programs NOT in `DATA_LABEL_PARAMETERS` now skip
  label recovery entirely (instead of using default `obs_labels` parameter)
- Safety net in `parameter_fixes.json` also strips obs_labels from autosol
- Files: `agent/command_builder.py`, `knowledge/parameter_fixes.json`,
  `tests/tst_command_builder.py`

**Fix - MR-SAD workflow skipping phaser and going straight to autosol**
- `after_program=phenix.autosol` directive was forcing autosol immediately,
  even before xtriage or phaser had run
- Root cause TWO-FOLD:
  1. `use_mr_sad` handling only removed autosol from `obtain_model` phase
  2. AutoSol was in YAML `obtain_model` phase with `has: anomalous` condition,
     entering valid_programs through the base path BEFORE `_apply_directives`
- Four fixes:
  1. `get_valid_programs()`: Added MR-SAD guard that removes autosol when
     has_search_model + has_anomalous + not phaser_done — runs BEFORE
     `_apply_directives` so autosol can't leak through the YAML path
  2. `_check_program_prerequisites()`: autosol requires xtriage_done, and
     for implicit MR-SAD (has_search_model + has_anomalous), also phaser_done
  3. `_apply_directives()`: use_mr_sad now removes autosol from ALL phases
     when phaser hasn't run (not just obtain_model)
  4. Directive extraction prompt: Added "CRITICAL: MR-SAD workflow" guidance
- Standalone SAD (no search model) unaffected by the guard
- Files: `agent/workflow_engine.py`, `agent/directive_extractor.py`

### Dependency Cleanup

**Fix - Directive extractor inferring use_experimental_phasing from data**
- LLM was setting `use_experimental_phasing: True` before xtriage ran,
  based on data characteristics (wavelength, atom types) rather than
  explicit user request
- This caused `predict_and_build` to be deprioritized even when the case
  needed AlphaFold to generate a model
- Added CRITICAL guidance to LLM prompt: only set `use_experimental_phasing`
  and `use_mr_sad` when user EXPLICITLY requests SAD/MAD/experimental phasing
- The system already auto-detects anomalous signal via xtriage and adjusts
  the workflow through the `has_anomalous` context flag
- File: `agent/directive_extractor.py`

**Removed langchain-classic dependency**
- In langchain 1.0+, `langchain.chains` and `langchain.retrievers` were
  removed and moved to `langchain_classic`. Rather than depending on the
  legacy package, implemented the functionality directly:
- `analysis/summarizer.py`: Replaced `create_stuff_documents_chain` with
  direct document concatenation + LLM invoke (the function just joins
  document text and passes it through a prompt template)
- `rag/retriever.py`: Replaced `ContextualCompressionRetriever` with
  `_CompressionRetriever`, a minimal `BaseRetriever` subclass that
  retrieves docs then reranks via the compressor — uses `langchain_core`
  which is still maintained
- All `from langchain.` imports eliminated — code now only uses
  `langchain_core`, `langchain_community`, and provider packages
- Removed from README dependencies table
- `langchain-classic` package can now be uninstalled

### Test Infrastructure Fixes

**Fix - Unconditional mock modules breaking PHENIX environment tests**
- Five test files unconditionally overwrote real `libtbx` modules with mocks
- Added `if 'libtbx' not in sys.modules` guards to all mock blocks
- Set `__path__` on mock modules pointing to actual local directories so
  submodule imports (e.g. `libtbx.langchain.agent.utils`) resolve automatically
- Files: `tst_state_serialization.py`, `tst_command_builder.py`,
  `tst_file_categorization.py`, `tst_session_directives.py`

**Fix - `tst_advice_preprocessing.py` import failure in PHENIX**
- `sys.path.insert` was inside `except ImportError` block, so `from tests.tst_utils`
  was unfindable when libtbx imports succeeded
- Moved to top level

**Cleanup - Stale docstring run instructions**
- Updated `Run with: python tests/test_...` → `tst_...` across 21 files

## Version 112.1 (February 2025)

### Cryo-EM State Fix & MR-SAD Workflow Support

**Fix 1 - Cryo-EM workflow stuck at cryoem_initial (HIGH IMPACT)**
- `_detect_cryoem_phase()` gated all progress behind `mtriage_done` check
- Tutorials that skip mtriage (going straight to resolve_cryo_em) got permanently stuck
- Fix: "past analysis" check now also considers `resolve_cryo_em_done`, `dock_done`, `rsr_done`
- Files: `agent/workflow_engine.py`

**Feature 1 - MR-SAD experimental phasing workflow (NEW)**
- Added `experimental_phasing` phase to X-ray workflow for MR-SAD
- Workflow: xtriage → phaser (place model) → autosol (with partpdb_file=PHASER.pdb) → autobuild
- Added `partpdb_file` optional input to autosol (auto-filled from phaser_output category)
- Added `xray_mr_sad` state mapping and detection in workflow engine
- Phase detection: triggered when `phaser_done AND (has_anomalous OR use_mr_sad)`
- Added `use_mr_sad` directive to workflow_preferences
- Directive extractor: "MR-SAD" keywords set `use_mr_sad=true`, do NOT set `after_program=autosol`
- In obtain_model phase with use_mr_sad: phaser prioritized, autosol removed
- Files: `agent/workflow_engine.py`, `knowledge/programs.yaml`, `knowledge/workflows.yaml`,
  `agent/decision_config.json`, `agent/directive_extractor.py`

**Tests: 765 total (+8 new)**
- `test_mr_sad_after_phaser_with_anomalous` - MR-SAD state detection
- `test_mr_sad_not_triggered_without_anomalous` - No false positives
- `test_mr_sad_not_triggered_when_autosol_done` - Skip if already done
- `test_mr_sad_directive_prioritizes_phaser` - Directive removes autosol from obtain_model
- `test_normal_sad_still_works` - Normal SAD pathway unaffected
- `test_autosol_has_partial_model_config` - YAML config correctness
- `test_mr_sad_directive_overrides_no_anomalous` - Directive triggers without has_anomalous
- `test_experimental_phasing_yaml_structure` - YAML phase structure validation

**Documentation updated:**
- ARCHITECTURE.md: State diagram with MR-SAD path, state table
- USER_DIRECTIVES.md: use_mr_sad field, MR-SAD example (Example 5)
- TESTING.md: Test counts (757→765)

## Version 112 (February 2025)

### Summary Quality Improvements: Steps Table Metrics & Consistency Fixes

Seven issues identified and fixed through codebase-wide audit:

**Fix 1 - Steps Performed table shows actual metrics (HIGH IMPACT)**
- `_get_key_metric_for_step()` now uses pre-parsed `cycle["metrics"]` as primary source
- YAML patterns only match raw log format (e.g., `R-free =`) but result text has reformatted
  metrics (e.g., `R Free:`) — caused all YAML pattern extraction to fail silently
- Steps table now shows "R-free: 0.258" instead of fallback text like "Analyzed" or "Built model"
- Files: `agent/session.py`

**Fix 2 - Benign warnings no longer drop metrics (MODERATE)**
- `extract_metrics_from_log()` now called for ALL success paths, not just clean success
- Previously, benign warnings (e.g., "No array of R-free flags found") triggered a code path
  that skipped metrics extraction entirely
- First refinement run (which often recreates R-free flags) could lose all R-factor metrics
- Files: `programs/ai_agent.py`

**Fix 3 - Ligandfit output typed as "ligand" not "model" (MINOR)**
- `_describe_output_file()` now returns `type: "ligand"` for `ligand_fit_*.pdb` files
- Previously typed as "model", so `best_by_type["ligand"]` was never populated in fallback path
- Files: `agent/session.py`

**Fix 4 - Removed dead `ligand_fit_output` from `categories_to_show` (COSMETIC)**
- `ligand_fit_output` is a stage under `"model"` in best_files, not a top-level category key
- Entry could never match anything; ligandfit output found via cycle-scanning (added in v111)
- Files: `agent/session.py`

**Fix 5 - `_is_intermediate_file` case-sensitive patterns fixed (MINOR)**
- `"carryOn"` and `"CarryOn"` compared against `basename_lower` could never match
- Replaced with lowercase `"carryon"`
- Files: `agent/session.py`

**Fix 6 - `detect_program` distinguishes `autobuild_denmod` (MINOR)**
- Added `autobuild_denmod` / `maps_only` detection before generic `autobuild` match
- Also added separate branch in `extract_all_metrics` to prevent autobuild_denmod from
  running autobuild's metrics extractor (which would report misleading R-factors)
- Files: `phenix_ai/log_parsers.py`

**Fix 7 - Added missing YAML `log_parsing` sections (TECH DEBT)**
- Six programs declared `outputs.metrics` but had no `log_parsing` patterns:
  autobuild, autosol, dock_in_map, map_to_model, validation_cryoem, holton_geometry_validation
- Two programs had partial coverage: predict_and_build (missing map_cc),
  real_space_refine (missing clashscore)
- All now have matching YAML patterns, same type of gap that caused v111's predict_and_build issue
- Files: `knowledge/programs.yaml`


## Version 111 (February 2025)

### Summary Output Fixes: Missing R-free Statistics & Ligandfit Files

**Problem 1 - Missing Final Quality Statistics**
- predict_and_build runs internal refinement but R-free/R-work were not extracted from its logs
- `_extract_predict_and_build_metrics()` only extracted pLDDT, map_cc, residues_built
- Added: calls `_extract_refine_metrics()` to also capture R-factors
- Added: `log_parsing` section in `programs.yaml` for predict_and_build
- Added: `autobuild_denmod` entry in `metrics.yaml` step_metrics to prevent it from
  matching generic autobuild config

**Problem 2 - Missing Ligandfit Output in Key Files**
- `ligand_fit_output` is a subcategory under parent `"model"` in best_files_tracker,
  not a top-level key — `best_files.get("ligand_fit_output")` always returned None
- Ligandfit scores 90 vs refined model's 100, so never becomes best "model" entry
- Added: cycle-scanning logic to find ligandfit output files from successful cycles

**Bonus - Fixed fallback cycle filtering**
- `cycle.get("status") != "completed"` checked a field that never exists
- Changed to `"SUCCESS" not in result` to match how success is tracked everywhere else

Files: `phenix_ai/log_parsers.py`, `knowledge/programs.yaml`, `knowledge/metrics.yaml`, `agent/session.py`


## Version 110 (January 2025)

### Major Feature: Experimental Phasing (SAD/MAD) Workflow Support

**Problem**: After xtriage detected anomalous signal, the agent would correctly run `phenix.autosol` on a straight-through workflow. But if the session was restarted (after removing autosol run), the agent would fall back to `phenix.predict_and_build` instead of re-running autosol, even though:
1. The session had `use_experimental_phasing: True` directive
2. The xtriage cycle showed `has_anomalous: True` in its metrics

**Root Cause**: A key mismatch in the transport layer:
- `session.get_history_for_agent()` returned history entries with `"analysis"` key
- `build_request_v2()` looked for `"metrics"` key, so it got empty `{}`
- After transport, history had empty metrics
- `_analyze_history()` looked for `"analysis"` key but found `"metrics"` (now empty)

Result: The anomalous flags (`has_anomalous`, `anomalous_measurability`) were lost during session restart, so autosol's condition `has: anomalous` was not met.

**Solution**:

1. **Fixed key mismatch** (`agent/api_client.py`):
   - `build_request_v2()` now checks both keys: `h.get("analysis", h.get("metrics", {}))`

2. **Fixed history analysis** (`agent/workflow_state.py`):
   - `_analyze_history()` now checks both keys: `entry.get("analysis", entry.get("metrics", {}))`

3. **Added anomalous metrics extraction** (`agent/session.py`):
   - `_extract_metrics_from_result()` now extracts anomalous info from result text
   - Patterns for: "Anomalous Measurability: 0.15", "Has Anomalous: True", etc.

### Code Consolidation: History Analysis

**Problem**: Two overlapping functions were analyzing history:
- `_extract_history_info()` in `graph_nodes.py` - used by rules selector
- `_analyze_history()` in `workflow_state.py` - used by workflow engine

This duplication led to inconsistencies and the key mismatch bug.

**Solution**: Consolidated to single source of truth:

1. **`_analyze_history()` is authoritative** - extracts all metrics from history
2. **`build_context()` includes all info** - passes to workflow_state
3. **`_extract_history_info()` delegates** - now just extracts from pre-computed context

Before (duplicated computation):
```
History → _extract_history_info() → (recompute from history)
History → _analyze_history() → (compute from history)
```

After (single computation):
```
History → _analyze_history() → history_info
                                    ↓
                            build_context() → context
                                    ↓
                            workflow_state (includes context)
                                    ↓
                            _extract_history_info() ← pulls from context
```

### New Context Fields

Added to workflow context for SAD/MAD workflows and refinement strategy:

| Field | Source | Purpose |
|-------|--------|---------|
| `has_anomalous` | xtriage analysis | Enable autosol in obtain_model phase |
| `strong_anomalous` | measurability > 0.10 | Priority boost for autosol |
| `anomalous_measurability` | xtriage | Quantify anomalous signal strength |
| `anomalous_resolution` | xtriage | Anomalous signal extent |
| `has_ncs` | map_symmetry analysis | Enable NCS constraints in refinement |

### New Tests

Added `tests/tst_history_analysis.py` with tests for:
- Transport analysis/metrics key handling
- History analysis from both key types
- NCS extraction
- Weak vs strong anomalous signal
- Session metrics extraction
- Graph nodes context extraction
- Integration: anomalous workflow enablement

### Code Consolidation: Duplicate Code Removal

**Problem**: Several functions were duplicated across the codebase:
1. Session state building - identical code in LocalAgent and RemoteAgent
2. MTZ classification - triplicated logic in workflow_state, template_builder, best_files_tracker
3. Program detection - similar functions in log_parsers and program_registry

**Solution**: Created shared utilities and delegated to canonical implementations:

**1. Session State Building** (`agent/api_client.py`):
```python
def build_session_state(session_info, session_resolution=None):
    """Build session_state dict from session_info. Used by both agents."""
```
- LocalAgent and RemoteAgent now call this shared function
- Reduced ~30 lines of duplicated code

**2. MTZ Classification** (`agent/file_utils.py` - NEW):
```python
def classify_mtz_type(filepath):
    """Classify MTZ as data_mtz or map_coeffs_mtz."""

def get_mtz_stage(filepath, category):
    """Get specific stage like refine_map_coeffs, denmod_map_coeffs."""
```
- Single source of truth for MTZ classification patterns
- workflow_state.py, template_builder.py, best_files_tracker.py now delegate to this
- Also includes is_mtz_file(), is_model_file(), is_map_file(), is_sequence_file()

**3. Program Detection** (`phenix_ai/log_parsers.py`):
- Made detect_program() the canonical implementation
- program_registry._detect_program() now delegates to it
- Added molprobity and polder detection to the canonical function

### New Test Suites

- `tests/tst_file_utils.py`: 12 tests for MTZ classification and file type detection

### Files Changed

- `agent/api_client.py`: Added build_session_state(), fixed analysis/metrics key lookup
- `agent/file_utils.py`: NEW - shared file classification utilities
- `agent/workflow_state.py`: Uses shared classify_mtz_type, fixed key lookup, added has_ncs/strong_anomalous
- `agent/workflow_engine.py`: Added has_ncs to context
- `agent/template_builder.py`: Uses shared classify_mtz_type
- `agent/best_files_tracker.py`: Uses shared classify_mtz_type
- `agent/program_registry.py`: Delegates to log_parsers.detect_program
- `agent/session.py`: Added anomalous metrics extraction
- `agent/graph_nodes.py`: Simplified _extract_history_info to use context
- `phenix_ai/local_agent.py`: Uses build_session_state
- `phenix_ai/remote_agent.py`: Uses build_session_state
- `phenix_ai/log_parsers.py`: Added molprobity and polder detection
- `tests/tst_history_analysis.py`: New test file
- `tests/tst_file_utils.py`: New test file
- `tests/run_all_tests.py`: Added History Analysis and File Utils test suites

---

### Major Feature: Dual MTZ Tracking

**Problem**: The codebase treated all MTZ files as a single category, but there are fundamentally two types:
- **Data MTZ**: Contains measured Fobs and R-free flags (for refinement)
- **Map Coefficients MTZ**: Contains calculated phases (for ligand fitting, visualization)

This caused issues with programs like `phenix.ligandfit` potentially receiving data MTZ instead of map coefficients.

**Solution**: Split MTZ tracking into two explicit categories with different update behaviors:

| Category | Update Rule | Use Case |
|----------|-------------|----------|
| `data_mtz` | First with R-free **locks forever** | Consistent R-free flags across refinement |
| `map_coeffs_mtz` | **Most recent wins** | Maps improve with refinement, use latest |

### Key Changes

**New Categories** (`knowledge/file_categories.yaml`):
- `data_mtz` (parent): original_data_mtz, phased_data_mtz
- `map_coeffs_mtz` (parent): refine_map_coeffs, denmod_map_coeffs, predict_build_map_coeffs

**Updated Programs** (`knowledge/programs.yaml`):
| Program | MTZ Input | Purpose |
|---------|-----------|---------|
| `phenix.xtriage` | `data_mtz` | Analyze Fobs |
| `phenix.phaser` | `data_mtz` | Molecular replacement |
| `phenix.refine` | `data_mtz` | Refine against Fobs |
| `phenix.polder` | `data_mtz` | Calculate omit maps |
| `phenix.ligandfit` | `map_coeffs_mtz` | Fit ligands (needs calculated phases) |

**BestFilesTracker** (`agent/best_files_tracker.py`):
- Added `_evaluate_data_mtz()`: Locks on first R-free flags
- Added `_evaluate_map_coeffs_mtz()`: Prefers most recent cycle
- Added `_classify_mtz_type()`: Auto-classifies by filename pattern

### Additional Fixes in v110

**Stepwise Mode / Automation Path**:
- Added `automation_path` to workflow state ("automated" or "stepwise")
- In stepwise mode, `predict_and_build` stops after prediction
- User then proceeds with `process_predicted_model` → `phaser` → `refine`
- Prevents duplicate predict_and_build runs

**Fallback Program Tracking**:
- Fallback node now correctly sets `program` field in state
- Fixes mismatch where PLAN showed one program but command was different
- Response builder uses `state["program"]` over `intent["program"]` when fallback used

**AutoBuild Scoring**:
- `autobuild_output` stage score increased to 100 (same as `refined`)
- AutoBuild runs internal refinement, so outputs are effectively refined models
- AutoBuild with better R-free now correctly beats earlier refine outputs

**Session Summary Best Files**:
- Removed file existence check in `_get_final_output_files()`
- Files created on client machine don't exist on server
- Now correctly shows best files in markdown summaries

### MTZ Classification Patterns

| Pattern | Category | Stage |
|---------|----------|-------|
| `refine_*_001.mtz` | map_coeffs_mtz | refine_map_coeffs |
| `*map_coeffs*.mtz` | map_coeffs_mtz | varies |
| `*denmod*.mtz` | map_coeffs_mtz | denmod_map_coeffs |
| `*_data.mtz` | data_mtz | original_data_mtz |
| Everything else | data_mtz | data_mtz |

### Files Changed

- `knowledge/file_categories.yaml`: New data_mtz and map_coeffs_mtz hierarchies
- `knowledge/programs.yaml`: All 10 programs updated
- `knowledge/metrics.yaml`: Scoring config for both MTZ types, autobuild_output score
- `knowledge/workflows.yaml`: Updated polder conditions
- `agent/best_files_tracker.py`: New evaluation methods, autobuild scoring
- `agent/workflow_state.py`: Updated parent categories, automation_path
- `agent/workflow_engine.py`: automation_path in context, stepwise mode handling
- `agent/command_builder.py`: Updated slot mappings
- `agent/graph_nodes.py`: Updated sanity context, fallback program tracking
- `agent/rules_selector.py`: Updated file selection
- `agent/session.py`: New get_best_data_mtz(), get_best_map_coeffs_mtz(), fixed best files display
- `agent/template_builder.py`: Updated category detection
- `agent/program_registry.py`: Updated phaser command
- `agent/directive_extractor.py`: Updated file preferences
- `knowledge/prompts_hybrid.py`: Updated recommended files display
- `phenix_ai/run_ai_agent.py`: Use state["program"] for response
- `tests/tst_best_files_tracker.py`: All 48 tests updated + autobuild scoring tests
- `tests/tst_workflow_state.py`: Added stepwise mode tests

---

## Version 109 (January 2025)

### Bug Fix: Empty Directives When Using Ollama Provider

**Problem**: When running with `provider=ollama`, directive extraction returned empty `{}` even when the user's advice clearly specified workflow instructions like "include ligand fitting".

**Root Cause**: Smaller local models (like llama3.2) may not follow complex JSON extraction prompts as reliably as GPT-4 or Gemini. The LLM might return:
- Empty JSON `{}`
- Malformed JSON that fails to parse
- Valid JSON but with content that gets filtered during validation

**Solution**: Multiple improvements for ollama reliability:

1. **Fallback to simple pattern extraction**: When ollama's LLM returns empty or fails, automatically fall back to `extract_directives_simple()` which uses regex patterns

2. **Added prefer_programs patterns**: New patterns to detect workflow preferences like:
   - "include ligand fitting" → `prefer_programs: [phenix.ligandfit]`
   - "fit the ligand" → `prefer_programs: [phenix.ligandfit]`
   - "with ligand fitting" → `prefer_programs: [phenix.ligandfit]`
   - "calculate polder map" → `prefer_programs: [phenix.polder]`

3. **Better logging**: Added debug logging to show:
   - What provider is being used
   - Response length from LLM
   - What sections were parsed
   - Preview of response when parsing fails

**Behavior with ollama**:
```
DIRECTIVES: Got response from ollama (500 chars)
DIRECTIVES: Parsed to empty dict - LLM may not have found actionable directives
DIRECTIVES: Trying simple pattern extraction as ollama fallback
DIRECTIVES: Simple extraction found: ['workflow_preferences']
```

### Files Changed

- `agent/directive_extractor.py`:
  - Added ollama fallback to `extract_directives()` 
  - Added `prefer_program_patterns` to `extract_directives_simple()`
  - Improved logging throughout

---

## Version 108 (January 2025)

### Multiple Summary and Advice Filtering Fixes

**Issue 1: "Predicted model" shown for predict_and_build output**

- Problem: Step metric showed "Predicted model" instead of refinement metrics for full predict_and_build runs
- Fix: Changed `metrics.yaml` step_metrics for predict_and_build to use `r_free` as primary metric with "Built model" fallback

**Issue 2: Key Output Files shows all files instead of best files**

- Problem: Summary listed many files from last cycle instead of the actual best files
- Fix: Modified `_get_final_output_files()` in `session.py` to prioritize `best_files` from session (model, mtz, map, sequence) before falling back to cycle outputs

**Issue 3: Ligandfit filtered out by user advice**

- Problem: User advice like "refinement with ligand fitting" was filtering programs to only `phenix.refine` because it contained "refine"
- Fix: Extended multi-step workflow detection in `_apply_user_advice()`:
  - Added sequencing words: "with", "including", "include", "plus", "also", "workflow", "sequence", "steps", "primary goal", "goal:"
  - Added check for multiple programs mentioned (if 2+ programs mentioned, don't filter)

**Issue 4: Empty directives when using ollama provider**

- Problem: Directive extraction returned `{}` with ollama because the server doesn't support ollama
- Root cause: Even though `run_on_server=False` was set for ollama, the code fell through to server execution when no local RAG database was found
- Fix: Added explicit check in `run_job_on_server_or_locally()` to honor `run_on_server=False` for directive_extraction mode, bypassing the database check since directive extraction only needs the LLM, not the RAG database

### Files Changed

- `knowledge/metrics.yaml` - Fixed predict_and_build step metrics
- `agent/session.py` - Modified `_get_final_output_files()` to use best_files first
- `agent/rules_selector.py` - Extended multi-step detection in `_apply_user_advice()`
- `programs/ai_analysis.py` - Added directive_extraction bypass for local-only mode

---

## Version 107 (January 2025)

### Feature: Graceful Stop on Persistent LLM Failures

**Problem**: When the LLM service was unavailable (rate limited, overloaded, API errors), the agent would silently fall back to rules-only mode without informing the user. This could lead to unexpected behavior.

**Solution**: After 3 consecutive LLM failures, the agent now stops gracefully with a helpful message:

```
The LLM service (google) is currently unavailable after 3 attempts.
Last error: 503 UNAVAILABLE - The model is overloaded

Options:
  1. Wait and try again later
  2. Run with --use_rules_only=True to continue without LLM
  3. Check your API key and network connection
```

**Behavior**:
- Failures 1-2: Fall back to rules-based planning for that cycle, continue workflow
- Failure 3+: Stop gracefully with helpful message
- On success: Reset failure counter

**Implementation**:
- Added `_handle_llm_failure()` function in `graph_nodes.py`
- Tracks `llm_consecutive_failures` in state
- Emits `STOP_DECISION` event with `llm_unavailable=True` flag

### Files Changed

- `agent/graph_nodes.py` - Added `_handle_llm_failure()`, failure tracking, graceful stop

---

## Version 106 (January 2025)

### Bug Fix: Informational Program Mentions Don't Block Workflow

**Problem**: When the LLM's processed advice mentioned programs like `phenix.elbow` or `phenix.ready_set` as suggestions (e.g., "if not provided, generate a CIF restraints file with phenix.elbow"), the directive validator treated these as explicit program requests and blocked the workflow with:

```
DIRECTIVE VALIDATION FAILED
Program 'phenix.elbow' exists in PHENIX but is not available in the AI agent workflow.
```

**Solution**: Program mentions in user advice text are now converted to warnings instead of blocking issues. Only programs explicitly requested in the directives structure (program_settings, stop_conditions.after_program, etc.) will block the workflow.

**Behavior Change**:
- Before: Any mention of unavailable program → VALIDATION FAILED
- After: Text mentions → Warning only, workflow continues

### Files Changed

- `agent/directive_validator.py` - Section 1 (user advice program check) now produces warnings not issues

---

## Version 105 (January 2025)

### Bug Fix: Failed Programs Don't Count as Done

**Problem**: Programs marked with `run_once: true` were being marked as "done" even when they failed, preventing retry attempts.

**Solution**: All `*_done` flags now only get set if the program completed successfully. Failed runs are skipped.

**Failure Detection**: Uses specific patterns to avoid false positives:
- `FAILED`, `SORRY:`, `SORRY `, `ERROR:`, `ERROR `, `: ERROR`, `TRACEBACK`, `EXCEPTION`
- Does NOT match "No ERROR detected" or similar success messages

**Programs Fixed**:

In `knowledge/program_registration.py` (auto-detected `run_once` programs):
- `detect_programs_in_history()` now checks result for failure patterns

In `agent/workflow_state.py` (manually tracked programs):
- Added `_is_failed_result()` helper function
- All `*_done` flag assignments now use this helper:
  - `validation_done`, `phaser_done`, `predict_done`, `predict_full_done`
  - `process_predicted_done`, `autobuild_done`, `autobuild_denmod_done`
  - `autosol_done`, `refine_done`, `rsr_done`, `ligandfit_done`
  - `pdbtools_done`, `dock_done`, `map_to_model_done`
  - `resolve_cryo_em_done`, `map_sharpening_done`

**Behavior Change**:
- Before: `ligandfit` fails → `ligandfit_done=True` → Cannot retry
- After: `ligandfit` fails → `ligandfit_done=False` → Can retry

### Files Changed

- `knowledge/program_registration.py` - Added specific failure patterns to `detect_programs_in_history()`
- `agent/workflow_state.py` - Added `_is_failed_result()` helper, refactored all done flag checks

---

## Version 104 (January 2025)

### Change: Remove Ligandfit Label Specification

**Rationale**: `phenix.ligandfit` can auto-detect map coefficient labels from the MTZ file. Manually specifying labels based on filename patterns was error-prone and could cause failures when the patterns didn't match the actual MTZ contents.

**Changes**:
- Removed `file_info.input_labels` from ligandfit defaults
- Removed label-switching invariants (`denmod_labels`, `predict_and_build_labels`)
- Kept commented-out code in programs.yaml for future restoration if needed

**Before**:
```
phenix.ligandfit model=... data=... ligand=... file_info.input_labels="FP PHIFP" general.nproc=4
```

**After**:
```
phenix.ligandfit model=... data=... ligand=... general.nproc=4
```

### Files Changed

- `knowledge/programs.yaml` - Removed labels from defaults, commented out invariants

---

## Version 103 (January 2025)

### Bug Fix: Intermediate File Handling

Fixed several issues where intermediate/temporary files were incorrectly used as inputs or tracked as best files.

**Issues Fixed:**

1. **Elbow ligand files used instead of fitted ligands**
   - `LIG_lig_ELBOW.*.pdb` files are geometry-optimized ligands, not fitted ligands
   - Added `*ELBOW*` and `*/TEMP*/*` exclusions to `ligand_pdb` category

2. **Superposed predicted models used as best model**
   - `*superposed_predicted_models*` files are alignment intermediates
   - Added exclusion to `predicted` category
   - Added to intermediate patterns in `best_files_tracker.py`

3. **Reference and EDITED files used as final outputs**
   - `*reference*` files are intermediate templates
   - `*EDITED*` files are intermediate edits
   - Added exclusions to `refined` category
   - Added to intermediate patterns in `best_files_tracker.py`

4. **with_ligand models now preserved**
   - Added `with_ligand` to valuable_patterns so combined protein+ligand models are tracked

### Files Changed

- `knowledge/file_categories.yaml`:
  - Added ELBOW/TEMP exclusions to `ligand_pdb`
  - Added `superposed_predicted_models` exclusion to `predicted`
  - Added patterns to `intermediate` category
  - Added exclusions to `refined` category
- `agent/best_files_tracker.py`:
  - Extended `_is_intermediate_file()` with more patterns
  - Added `with_ligand` to valuable patterns

---

## Version 102 (January 2025)

### Bug Fix: Retry on Model Overload (503 UNAVAILABLE)

**Problem**: When Gemini returns a 503 UNAVAILABLE error ("The model is overloaded"), the summarization would fail immediately without retrying.

**Solution**: 
1. Added "503" and "unavailable" to rate limit indicators in `rate_limit_handler.py`
2. Added retry logic to `summarize_log()` with exponential backoff (2s, 4s, 8s) for rate limit and overload errors

### Behavior

When model is overloaded:
```
Summarizing log file (using cheap model)...
Model overloaded/rate limited, waiting 2s before retry...
Summarization retry 2/3...
Model overloaded/rate limited, waiting 4s before retry...
Summarization retry 3/3...
[success or final error]
```

### Files Changed

- `agent/rate_limit_handler.py` - Added "503" and "unavailable" to rate limit indicators
- `phenix_ai/run_ai_analysis.py` - Added retry logic to `summarize_log()` with exponential backoff

---

## Version 101 (January 2025)

### Bug Fix: HTML Summary Table Formatting After Failed Steps

**Problem**: When a step failed, the error message in the "Key Metric" column could contain newlines or pipe characters, breaking the markdown table formatting for all subsequent rows.

Example broken output:
```
| 3 | ligandfit | ✗ | FAILED: Sorry: Sorry LigandFit failed
Please... | | 4 | pdbtools | ✓ | Completed |
```

**Solution**: Sanitize the `key_metric` field before adding to table:
- Replace newlines with spaces
- Replace pipe characters (`|`) with dashes
- Collapse multiple spaces
- Truncate to 60 characters max

### Files Changed

- `agent/session.py` - Sanitize key_metric in `_format_summary_as_markdown()`

---

## Version 100 (January 2025)

### Bug Fix: Ligand File Detection Pattern

**Problem**: `7qz0_ligand.pdb` wasn't detected as a ligand file because patterns only matched `ligand_*.pdb` (ligand at start), not `*_ligand.pdb` (ligand at end).

**Solution**: Added patterns `*_ligand.pdb` and `*_lig.pdb` to `ligand_pdb` category.

### Bug Fix: predict_and_build Map Coefficients for Ligandfit

**Problem**: After `predict_and_build`, ligandfit was using `*_refinement.mtz` which doesn't contain map coefficients. The correct file is `*_map_coeffs.mtz` with labels `FP PHIFP`.

**Solution**:
1. Added new `predict_and_build_mtz` category for `*map_coeffs*.mtz` files
2. Updated ligandfit input priorities to include `predict_and_build_mtz`
3. Added `denmod_mtz` and `predict_and_build_mtz` to `specific_subcategories` in command builder so category-based selection takes priority over `best_files`
4. Added invariant to switch labels to `"FP PHIFP"` when using map_coeffs files

### Files Changed

- `knowledge/file_categories.yaml`:
  - Added `*_ligand.pdb` and `*_lig.pdb` patterns to `ligand_pdb`
  - Added `predict_and_build_mtz` category for map coefficient files
- `knowledge/programs.yaml`:
  - Updated ligandfit mtz priorities: `[denmod_mtz, predict_and_build_mtz, refined_mtz]`
  - Added invariant for predict_and_build labels (`FP PHIFP`)
- `agent/command_builder.py`:
  - Added `denmod_mtz` and `predict_and_build_mtz` to `specific_subcategories`
- `tests/tst_file_categorization.py`:
  - Added `test_ligand_file_patterns`
  - Added `test_predict_and_build_mtz_detection`

---

## Version 99 (January 2025)

### Feature: maximum_automation=False Now Works for X-ray

**Previously**: `maximum_automation=False` (stepwise mode) only affected cryo-EM workflows, forcing `stop_after_predict=True` for `predict_and_build`.

**Now**: Stepwise mode also applies to X-ray workflows. When `maximum_automation=False`:
- `predict_and_build` will use `stop_after_predict=True` in states: `xray_initial`, `xray_placed`
- This gives users more control over the workflow with intermediate checkpoints
- User can then run `process_predicted_model` → `phaser` → `refine` separately

### Usage

```bash
# Stepwise mode - more control with intermediate checkpoints
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"
```

### Workflow Comparison

**Automated (maximum_automation=True, default)**:
```
xray_initial → xtriage → predict_and_build(full) → xray_refined
```

**Stepwise (maximum_automation=False)**:
```
xray_initial → xtriage → predict_and_build(stop_after_predict)
                              ↓
              process_predicted_model → phaser → refine → xray_refined
```

### Files Changed

- `agent/graph_nodes.py` - Extended stepwise mode handling to X-ray states
- `agent/docs_tools.py` - Updated workflow documentation diagrams
- `agent/workflow_state.py` - Updated stepwise hint message
- `docs/README.md` - Added automation modes section and quick start example
- `tests/tst_integration.py` - Added `test_xray_stepwise_forces_stop_after_predict`
- `tests/tst_workflow_state.py` - Added X-ray stepwise tests

---

## Version 98 (January 2025)

### Bug Fix: predict_and_build Counts as Refinement for ligandfit

**Problem**: After running `phenix.predict_and_build`, `phenix.ligandfit` was unavailable because `refine_count=0`:
```
PERCEIVE: phenix.ligandfit unavailable: refine_count=0 does not satisfy condition '> 0'
```

But predict_and_build includes internal refinement cycles and produces the same outputs (map coefficients) that ligandfit needs.

**Solution**: When a full `predict_and_build` run completes (not stopped early), increment `refine_count` so downstream programs like `ligandfit` know there's a refined model with map coefficients.

### Files Changed

- `agent/workflow_state.py` - Increment `refine_count` for successful full predict_and_build runs

---

## Version 97 (January 2025)

### Bug Fix: PredictAndBuild Output Categorized as Model (Not Search Model)

**Problem**: Files like `PredictAndBuild_0_overall_best.pdb` were incorrectly categorized as `search_model` instead of `model`. This happened because:
1. The file contains "predict" in the name
2. Multiple categories had `excludes: ["*predict*"]`
3. The file fell through to the `predicted` category (parent: `search_model`)

This caused the sanity check to fail with:
```
PERCEIVE: RED FLAG [search_model_not_positioned]: Cannot refine: search model found but not yet positioned
```

**Solution**: Added new `predict_and_build_output` category that specifically matches `PredictAndBuild_*_overall_best*.pdb` files and categorizes them as `model` (positioned, ready for refinement).

Also added exclusions to the `predicted` category to ensure these files don't get double-categorized.

### Files Changed

- `knowledge/file_categories.yaml` - Added `predict_and_build_output` category, added exclusions to `predicted`
- `tests/tst_file_categorization.py` - Added test for PredictAndBuild output categorization

---

## Version 96 (January 2025)

### Bug Fix: LLM Slot Alias Mapping for MTZ Files

**Problem**: When the LLM requested an MTZ file using the slot name `data` (e.g., `data=PredictAndBuild_0_overall_best_refinement.mtz`), but the program defined the input slot as `mtz`, the LLM's file choice was ignored and the wrong file was auto-selected.

**Example Debug Output (Before)**:
```
BUILD: LLM requested files: {model=..., data=PredictAndBuild_0_overall_best_refinement.mtz}
BUILD: Skipping best_files for mtz (program needs specific subcategory)
BUILD: Auto-filled mtz=PredictAndBuild_0_refinement_cycle_2.extended_r_free.mtz  # WRONG!
```

**Solution**: Added `SLOT_ALIASES` mapping that translates common LLM slot names to canonical program input names:
- `data` → `mtz`
- `pdb` → `model`
- `seq_file` → `sequence`
- etc.

**Example Debug Output (After)**:
```
BUILD: LLM requested files: {model=..., data=PredictAndBuild_0_overall_best_refinement.mtz}
BUILD: Mapped LLM slot 'data' to 'mtz'
BUILD: Using LLM-selected file for mtz
```

### Files Changed

- `agent/command_builder.py` - Added `SLOT_ALIASES` dict and updated LLM file processing to use aliases
- `tests/tst_command_builder.py` - Added tests for slot alias mapping

---

## Version 94 (January 2025)

### New Feature: Explain Why Programs Are Unavailable

When a program like `phenix.ligandfit` is not available, the debug output now explains WHY:

```
PERCEIVE: Valid programs: phenix.refine, phenix.polder, STOP
PERCEIVE: phenix.ligandfit unavailable: missing required file: ligand_file
```

This helps diagnose issues when expected programs aren't offered.

### Possible Explanations

- `missing required file: ligand_file` - No ligand file (.pdb/.cif) detected
- `missing required file: sequence` - No sequence file (.fa/.seq) detected
- `already completed: ligandfit` - Program already ran (not_done condition failed)
- `r_free=0.40 does not satisfy condition '< 0.35'` - Metric threshold not met
- `refine_count=0 does not satisfy condition '> 0'` - Needs refinement first
- `run_once program already executed` - Program like xtriage already ran

### Files Changed

- `agent/workflow_engine.py` - Added `explain_unavailable_program()` method, added `unavailable_explanations` to workflow state
- `agent/graph_nodes.py` - Added debug logging for unavailable programs

---

## Version 93 (January 2025)

### New Feature: Density Modification Workflow for X-ray

- **Added `phenix.autobuild_denmod` to X-ray workflow**: Before ligand fitting, the agent can now run density modification using `phenix.autobuild maps_only=True`. This creates improved map coefficients (`overall_best_denmod_map_coeffs.mtz` with FWT/PHFWT labels) for better ligand fitting.

### Workflow Changes

The X-ray refine phase now includes:
1. `phenix.refine` (preferred) - standard refinement
2. `phenix.autobuild` - when R-free stuck above threshold
3. **`phenix.autobuild_denmod`** (NEW) - density modification before ligandfit
4. `phenix.ligandfit` - fit ligand when model is good enough
5. `phenix.polder` - omit map calculation

### Technical Changes

- **Prompt clarification**: Added warning that `predict_and_build` is NOT for density modification
- **New file category**: `denmod_mtz` for density-modified MTZ files
- **Ligandfit label switching**: Automatically uses `FWT PHFWT` labels when denmod MTZ is selected
- **Done flag tracking**: Added `autobuild_denmod_done` flag

### Documentation Updates

- **docs/OVERVIEW.md**: Updated workflow example with autobuild_denmod, updated done flags table
- **docs/guides/ADDING_PROGRAMS.md**: Added autobuild_denmod_done to flags table, added note about refine_count/rsr_count
- **docs/guides/TESTING.md**: Updated test table with counts and added "Key Tests for Recent Fixes" section
- **tests/run_all_tests.py**: Updated docstring with current test list and key tests

### Files Changed

- `knowledge/workflows.yaml` - Added autobuild_denmod to refine phase
- `knowledge/programs.yaml` - Added denmod_labels invariant to ligandfit
- `knowledge/file_categories.yaml` - Added denmod_mtz category
- `knowledge/prompts_hybrid.py` - Added autobuild_denmod description, warned about predict_and_build
- `agent/workflow_state.py` - Added autobuild_denmod_done flag and detection
- `agent/command_builder.py` - Added file_matches invariant handling
- `docs/OVERVIEW.md` - Updated workflow examples and done flags
- `docs/guides/ADDING_PROGRAMS.md` - Updated done flags table
- `docs/guides/TESTING.md` - Updated test descriptions
- `tests/run_all_tests.py` - Updated docstring

---

## Version 72 (January 2025)

### Bug Fixes

- **Fixed cycle count to exclude STOP cycles (v92)**: The session summary was reporting "Cycles: 5 (4 successful)" when 4 programs ran and cycle 5 was just STOP. Now STOP cycles are excluded from the count, so it correctly reports "Cycles: 4 (4 successful)".

### Tests Added

- `test_stop_cycle_excluded_from_count` - Verifies STOP cycles are excluded from total_cycles and successful_cycles
- `test_cryoem_done_flags` - Verifies done flags are set for cryo-EM programs (resolve_cryo_em_done, map_sharpening_done, map_to_model_done, dock_done)

### Documentation Updated

- `docs/OVERVIEW.md` - Added "Program Execution Controls" section documenting `not_done` conditions and `run_once` flags with complete table of done flags
- `docs/guides/ADDING_PROGRAMS.md` - Added "Available done flags" table showing all manually-tracked done flags
- `docs/guides/TESTING.md` - Updated test count table with new test files

### Files Changed

- `agent/session.py` - Exclude STOP/None/unknown from total_cycles and successful_cycles
- `agent/session_tools.py` - Exclude STOP cycles from print_session_summary()
- `tests/tst_session_summary.py` - Added test_stop_cycle_excluded_from_count
- `tests/tst_workflow_state.py` - Added test_cryoem_done_flags

---

## Version 71 (January 2025)

### Bug Fixes

- **Fixed predict_and_build running without resolution for X-ray (v91)**: When user provides `program_settings` for predict_and_build (e.g., `rebuilding_strategy=Quick`), the program was being added to valid_programs even before xtriage ran. This caused predict_and_build to run without resolution, forcing `stop_after_predict=True` (prediction only). Now `_check_program_prerequisites` requires xtriage_done (X-ray) or mtriage_done (cryo-EM) before adding predict_and_build from program_settings.

- **Fixed resolution requirement for predict_and_build full workflow**: The command builder now correctly requires resolution for BOTH X-ray and cryo-EM when `stop_after_predict=False`. If resolution is not available, it forces `stop_after_predict=True` with a message suggesting to run xtriage/mtriage first.

### Root Cause

The workflow engine's `_apply_directives` was adding `predict_and_build` to valid_programs whenever the user had `program_settings` for it:

```python
# Before: Always allowed predict_and_build
if program == "phenix.predict_and_build":
    return True  # Always allow - worst case it does prediction-only
```

This bypassed the normal workflow phase ordering (xtriage → obtain_model).

### The Fix

```python
# After: Require xtriage/mtriage to be done first
if program == "phenix.predict_and_build":
    if phase_name in ("obtain_model", "molecular_replacement", "dock_model"):
        return True  # Let the phase conditions handle it
    if context.get("xtriage_done") or context.get("mtriage_done"):
        return True
    return False  # Don't add to early phases
```

### Correct Workflow Now

1. xtriage runs first → extracts resolution
2. predict_and_build runs with resolution → full workflow (prediction + MR + building)

### Files Changed

- `agent/workflow_engine.py` - Fixed `_check_program_prerequisites` to require xtriage/mtriage
- `agent/command_builder.py` - Fixed `_apply_invariants` to require resolution for full workflow

---

## Version 70 (January 2025)

### Bug Fixes

- **Fixed programs running repeatedly without stopping (v90)**: Multiple programs were missing `not_done` conditions in their workflow definitions, allowing the LLM to choose them repeatedly even after successful completion. Added protection to all one-time-run programs.

### Programs Now Protected from Re-runs

**X-ray workflow:**
- `phenix.predict_and_build` - `not_done: predict_full`
- `phenix.phaser` - `not_done: phaser` (in both obtain_model and molecular_replacement phases)

**Cryo-EM workflow:**
- `phenix.predict_and_build` - `not_done: predict`
- `phenix.dock_in_map` - `not_done: dock`
- `phenix.map_to_model` - `not_done: map_to_model`
- `phenix.resolve_cryo_em` - `not_done: resolve_cryo_em` (in all 4 phases where it appears)
- `phenix.map_sharpening` - `not_done: map_sharpening` (in all 4 phases where it appears)

**Previously protected (unchanged):**
- `phenix.autobuild`, `phenix.autosol`, `phenix.ligandfit`, `phenix.map_symmetry`, `phenix.process_predicted_model`

### Programs Intentionally Without Protection (run multiple times)

- `phenix.refine` / `phenix.real_space_refine` - Iterative refinement
- `phenix.molprobity` / `phenix.model_vs_data` / `phenix.validation_cryoem` - Validation after each cycle
- `phenix.polder` - May run for different sites
- `phenix.pdbtools` - May add multiple ligands
- `phenix.xtriage` / `phenix.mtriage` - Already have `run_once: true` in programs.yaml

### Files Changed

- `knowledge/workflows.yaml` - Added `not_done` conditions to 11 program entries
- `agent/workflow_state.py` - Added done flags for `map_to_model`, `resolve_cryo_em`, `map_sharpening`

---

## Version 69 (January 2025)

### Bug Fixes

- **Fixed predict_and_build forcing `stop_after_predict=True` for X-ray (v89)**: The command builder was forcing `stop_after_predict=True` whenever resolution wasn't in the context, even for X-ray where predict_and_build can read resolution directly from the MTZ file. This caused predict_and_build to run in prediction-only mode repeatedly (3 retry cycles due to duplicate detection). Now `stop_after_predict=True` is only forced for cryo-EM without resolution.

### Root Cause

The command builder had this logic:
```python
if program == "phenix.predict_and_build":
    if not context.resolution and "stop_after_predict" not in strategy:
        strategy["stop_after_predict"] = True
```

This applied to BOTH X-ray and cryo-EM, but:
- For X-ray: predict_and_build can determine resolution from the MTZ automatically
- For cryo-EM: mtriage should run first to get resolution

The fix restricts this to cryo-EM only:
```python
if context.experiment_type == "cryoem":
    if not context.resolution and "stop_after_predict" not in strategy:
        strategy["stop_after_predict"] = True
```

### Why 3 LLM Calls Were Happening

1. LLM chose predict_and_build
2. Command builder forced `stop_after_predict=True` (same as previous run)
3. Validate detected duplicate command → `validation_error`
4. Graph looped back to Plan (retry 1)
5. Same result → retry 2
6. Same result → retry 3 → fallback

### Files Changed

- `agent/command_builder.py` - Only force `stop_after_predict=True` for cryo-EM

---

## Version 68 (January 2025)

### Bug Fixes

- **Fixed ligandfit using input MTZ instead of refined MTZ (v88)**: When ligandfit requires `refined_mtz` (an MTZ with map coefficients), and no refined MTZ exists (because refinement failed), the code was falling back to extension-based matching and selecting the input MTZ. Now when a program requires a specific subcategory (like `refined_mtz`), extension-based fallback is disabled.

- **Fixed refine_count incrementing for failed refinements (v88)**: Previously, `refine_count` was incremented regardless of whether refinement succeeded or failed. This caused workflow conditions like `refine_count > 0` (used by ligandfit) to pass even when no successful refinement had occurred. Now only successful refinements increment the count.

- **Fixed rsr_count incrementing for failed RSR (v88)**: Same fix applied to `rsr_count` for `phenix.real_space_refine`.

### Root Cause

The session showed:
```
Categorized files: model=1, search_model=15, mtz=1, sequence=10, ...
```
No `refined_mtz` category because refinement failed. But `refine_count` was 2 (from failed attempts), so ligandfit was allowed. Then command_builder fell back to extension matching and selected the input MTZ.

### Files Changed

- `agent/command_builder.py` - Added check to prevent extension-based fallback when specific subcategory required
- `agent/workflow_state.py` - Added success checking for refine_count and rsr_count in `_analyze_history()`
- `tests/tst_workflow_state.py` - Added `test_failed_refine_not_counted`

---

## Version 67 (January 2025)

### Bug Fixes

- **Fixed stop_after_predict=True being suggested for X-ray (v87)**: The prompt was incorrectly telling the LLM to use `stop_after_predict=True` for ANY stepwise workflow, but this should only apply to cryo-EM stepwise workflows. For X-ray, `predict_and_build` should run the full workflow (prediction → MR → building). Added experiment_type check so the guidance only appears for cryo-EM.

- **Clarified predict_and_build documentation in prompts**: Added explicit notes explaining that by default `predict_and_build` runs the FULL workflow, and `stop_after_predict=True` should only be used for cryo-EM stepwise.

### Root Cause

When user says "stop after PredictAndBuild", the LLM was confusing:
- "Stop the agent workflow after predict_and_build completes" (correct interpretation)
- "Set stop_after_predict=True" (incorrect - this only runs prediction, skipping MR and building)

The prompt was adding `NOTE: Use predict_and_build with strategy: {"stop_after_predict": true}` for stepwise workflows without checking if it's cryo-EM or X-ray.

### Files Changed

- `knowledge/prompts_hybrid.py` - Added experiment_type check for stop_after_predict guidance; clarified predict_and_build documentation

---

## Version 66 (January 2025)

### Critical Bug Fixes

- **Fixed predicted model incorrectly becoming best model after refinement (v85)**: When `phenix.refine` runs, the directory scan picks up ALL files (including pre-existing ones like `PredictAndBuild_0_predicted_model_processed.pdb`). Previously, `record_result()` blindly applied `stage="refined"` to ALL PDB files in `output_files`, causing the predicted model to get a higher score than the actual refined model. Now `record_result()` only applies the program stage to files whose basename matches expected output patterns (e.g., only files containing "refine" get `stage="refined"`).

- **Fixed PHASER models getting inflated scores from refinement metrics**: Similar to above - `PHASER.1.pdb` in refinement's `output_files` was getting `stage="refined"` and metrics from phenix.refine, inflating its score from 70 to 132. Now handled by the same pattern matching fix.

- **Fixed STOP not available after user's workflow completes (v84)**: When user directives specify `after_program` (e.g., "stop after refinement"), STOP is now added to valid_programs after that program completes. Previously STOP was only available after validation, forcing unwanted extra cycles.

### Consistency Fixes (v86)

- **Synchronized pattern matching across all three locations**:
  - `session.py:_rebuild_best_files_from_cycles()` - Rebuild from saved session
  - `session.py:record_result()` - Real-time recording
  - `session_tools.py:rebuild_best_files()` - Manual rebuild tool

- **Added missing program-to-stage mappings**:
  - `session.py:_infer_stage_from_program()` - Added `predict_and_build`, `ligandfit`, `pdbtools`
  - `session_tools.py:infer_stage()` - Added `process_predicted_model`

- **Added missing filename patterns**:
  - `session_tools.py` - Added `processed_predicted` and `autobuild_output` patterns

### Bug Details

**The predicted model bug (v85)**:
- Root cause: `output_files` from directory scanning includes ALL files in the working directory, not just files created by the program
- `record_result()` was giving ALL PDB files the program's stage (e.g., "refined" for phenix.refine)
- This caused `PredictAndBuild_0_predicted_model_processed.pdb` to get `stage="refined"` with the refinement metrics, giving it score 133 vs PHASER.1.pdb's score 132
- Fix: Pattern matching in `record_result()` now matches `_rebuild_best_files_from_cycles()` - only files with matching basenames get the program stage

### Testing

- Added test `test_predicted_model_not_promoted_by_refine` to verify predicted models don't get wrongly promoted when they appear in refine's output_files
- Added test `test_phaser_model_not_promoted_by_refine_metrics` to verify PHASER models don't get refinement metrics
- Added test `test_stop_added_after_after_program_completes` to verify STOP is available after after_program completes

### Files Changed

- `agent/session.py` - Fixed `record_result()` to pattern-match PDB filenames before applying stage; added missing programs to `_infer_stage_from_program()`
- `agent/session_tools.py` - Added `process_predicted_model` to `infer_stage()`; added `processed_predicted` and `autobuild_output` filename patterns
- `agent/workflow_engine.py` - Added after_program completion check in `_apply_directives()`
- `tests/tst_best_files_tracker.py` - Added predicted model and PHASER model promotion tests
- `tests/tst_workflow_state.py` - Added STOP after after_program test

---

## Version 65 (January 2025)

### Major Bug Fixes

- **Fixed map files getting wrong stage in best_files (v79)**: All map files from resolve_cryo_em were incorrectly getting `stage=optimized_full_map` instead of being classified by filename patterns. This caused `initial_map.ccp4` to win over `denmod_map.ccp4`. Now map files are classified based on their basename: `initial*` → intermediate_map, `denmod*` → optimized_full_map, `sharp*` → sharpened, `half*` → half_map.

- **Fixed rebuild function missing MTZ/PDB stage patterns (v80)**: The `_rebuild_best_files_from_cycles()` function was missing:
  - Phased MTZ detection (`*phased*`, `*phases*`, `*solve*` → `phased_mtz`)
  - Generic MTZ fallback (→ `mtz` stage)
  - `processed_predicted` pattern for process_predicted_model outputs
  - `autobuild_output` pattern for autobuild outputs

- **Fixed autobuild_done set even on failure (v77)**: Previously `autobuild_done=True` was set just because autobuild appeared in history, regardless of success/failure. Now checks for "FAIL", "SORRY", "ERROR" in result before marking done, allowing the agent to try alternatives when autobuild fails.

- **Fixed LLM file suggestions with wrong extension (v77)**: Added extension validation when LLM suggests files for input slots. Now rejects files with wrong extension (e.g., `.ccp4` for model slot which expects `.pdb/.cif`).

- **Allow programs from directives even when workflow state is "past" that phase (v81)**: When user has `program_settings` for a program (e.g., `phenix.predict_and_build`), that program is now added to `valid_programs` even if the workflow state thinks we're past that phase. This allows users to explicitly request earlier-phase programs.

### New Features

- **Programs from program_settings added to valid_programs**: If user directives include settings for a specific program, that program is automatically added to the list of valid programs (subject to prerequisite checks). This respects user intent when they've configured a program they want to run.

- **_check_program_prerequisites() helper**: New method in WorkflowEngine that centralizes prerequisite checking for programs being added via directives. Checks:
  - Refinement programs need a model to refine
  - Ligandfit needs prior refinement (for map coefficients)
  - predict_and_build is always allowed (worst case: prediction-only)

### Testing

- Added tests for directive-based program addition in tst_workflow_state.py:
  - `test_program_settings_adds_program_to_valid`
  - `test_program_settings_prioritizes_program`
  - `test_program_settings_respects_prerequisites`
  - `test_default_program_settings_ignored`

### Files Changed

- `agent/workflow_engine.py` - Added `_check_program_prerequisites()`, enhanced `_apply_directives()`
- `agent/workflow_state.py` - Fixed autobuild_done to check for success
- `agent/session.py` - Fixed map file stage assignment, complete rebuild function
- `agent/command_builder.py` - Added LLM file extension validation
- `tests/tst_workflow_state.py` - Added directive program tests

---

## Version 64 (January 2025)

### Major Bug Fixes

- **Fixed best_files rebuild from cycle history (v69)**: Previously, `best_files` was persisted in session.json and could get stale when cycles were removed. Now `Session.load()` always rebuilds `best_files` from the cycle history, ensuring consistency. Removed redundant `best_files.evaluate_file()` call from `_track_output_files()` that was overwriting good evaluations with `metrics=None`.

- **Fixed smart stage assignment in rebuild (v71)**: When rebuilding best_files, stage was blindly applied to all files in a cycle's output_files. Now only applies program-specific stage (e.g., "refined") to files whose basename matches expected patterns (e.g., contains "refine"). This prevents PHASER.1.pdb from incorrectly getting `stage=refined` when it appears in a refinement cycle's output_files.

- **Fixed after_cycle directive for ligand workflows (v72)**: When user says "stop after second refinement" with a ligand workflow, LLM incorrectly extracted `after_cycle: 2`. Extended `_fix_ligand_workflow_conflict()` to also clear `after_cycle <= 4` when ligand constraints are present, since ligand workflows need ~8 cycles minimum.

- **Fixed ligandfit using wrong MTZ file (v73)**: phenix.ligandfit needs an MTZ with map coefficients (2FOFCWT, PH2FOFCWT) from refinement, but was getting the original data MTZ. Added logic to skip generic `best_files["mtz"]` when `input_priorities` specifies a specific subcategory like `refined_mtz`.

- **Fixed cryo-EM dock_in_map using initial_map instead of denmod_map (v74)**: After resolve_cryo_em runs, dock_in_map was selecting `initial_map.ccp4` (intermediate) instead of `denmod_map.ccp4` (density-modified output). Added `optimized_full_map` category with score 100, `intermediate_map` with score 5, and proper pattern matching for `denmod_map`, `density_modified`, etc.

- **Fixed pdbtools output naming**: Added `fixes.output_name` to pdbtools configuration to generate output filenames like `{protein_base}_with_ligand.pdb`, ensuring the combined model is properly categorized for downstream programs.

### New Categories & Scoring

**Map Categories (v74):**
| Stage | Score | Description |
|-------|-------|-------------|
| optimized_full_map | 100 | denmod_map, density_modified, sharpened |
| sharpened | 90 | Sharpened maps |
| full_map | 50 | Regular full reconstructions |
| half_map | 10 | Half-maps for FSC |
| intermediate_map | 5 | initial_map (resolve_cryo_em intermediate) |

### Testing

- Added tests for optimized_full_map scoring and classification
- Added tests for intermediate_map low priority
- Added tests for ligand workflow after_cycle clearing
- Added tests for docked model bubbling to model category

### Files Changed

- `agent/session.py` - Added `_rebuild_best_files_from_cycles()`, always rebuild on load
- `agent/session_tools.py` - Updated `rebuild_best_files()` with smart stage assignment
- `agent/command_builder.py` - Skip best_files for specific subcategories
- `agent/directive_extractor.py` - Extended ligand workflow fix for after_cycle
- `agent/best_files_tracker.py` - Added `intermediate_map` stage and scoring
- `knowledge/file_categories.yaml` - Added `optimized_full_map`, `intermediate_map` categories
- `knowledge/programs.yaml` - Updated resolve_cryo_em outputs, pdbtools output naming
- `programs/ai_agent.py` - Added `placed_model`/`docked` to valuable_output_patterns
- `tests/tst_best_files_tracker.py` - New map scoring tests
- `tests/tst_directive_extractor.py` - New ligand workflow tests
- `tests/tst_file_categorization.py` - New cryo-EM map categorization tests

---

## Version 63 (January 2025)

### Major Bug Fixes

- **Fixed predict_and_build output tracking**: Output files in `CarryOn/` directories were being incorrectly excluded as "intermediate files". Added `valuable_output_patterns` list that overrides intermediate exclusions for important outputs like `*_predicted_model*.pdb`, `*overall_best*.pdb`, `*_processed*.pdb`.

- **Fixed workflow state detection for predicted models**: The `_has_placed_model()` function was incorrectly returning `True` when user directives mentioned `phenix.refine` AND any PDB file existed - even if that PDB was a search_model (unpositioned). Now properly checks that PDB files are not in `search_model`, `predicted`, or `processed_predicted` categories before considering the model "placed".

- **Fixed pdbtools file selection**: Added `input_priorities` to `phenix.pdbtools` in programs.yaml to properly select refined model + ligand instead of predicted model. Added support for `priority_patterns` and `prefer_patterns` in command_builder.py.

- **Fixed ligand workflow directive conflict**: When user wants "refine, fit ligand, refine again", the directive extractor was incorrectly setting `after_program: phenix.refine` which stopped at the first refinement. Added `_fix_ligand_workflow_conflict()` post-processing that clears `after_program` when ligand-related constraints are present.

### New Features

- **Automatic Safety Documentation Generator** (`agent/generate_safety_docs.py`): Script that scans the codebase and generates a comprehensive table of all safety checks, validations, and post-processing corrections. Run with `python agent/generate_safety_docs.py > docs/SAFETY_CHECKS.md`.

- **Simplified Verbosity Levels**: Reduced from 4 levels to 3 levels (quiet/normal/verbose). The `debug` level has been removed; `verbose` now includes all detailed output including file selection, LLM traces, and internal state. For backwards compatibility, `debug` is accepted as input but treated as `verbose`.

- **File existence retry mechanism**: Added retry logic (3 attempts, 0.1s delay) in `resolve_file_path()` to handle race conditions where log mentions a file before it's fully written to disk.

- **Session file fsync**: Added explicit `fsync()` call when saving session files to ensure data is written to disk immediately.

### Documentation

- **New SAFETY_CHECKS.md**: Auto-generated documentation of all 70+ safety checks categorized by type:
  - Sanity Checks (Pre-Execution): 20
  - Directive Validation (Post-LLM): 7
  - File Validation: 4
  - Workflow State Validation: 8
  - Command Building Validation: 3
  - Input Validation: 29
  - Post-Processing Corrections: 4

- **New PROGRAM_CONFIG_ROBUSTNESS.md**: Implementation plan for making program configuration more robust with sensible defaults, validation warnings, and dry_run_file_selection mode.

### Files Changed

- `programs/ai_agent.py` - CarryOn fix, simplified verbosity, fsync
- `agent/workflow_engine.py` - Fixed `_has_placed_model()` to exclude search_models
- `agent/command_builder.py` - Added priority_patterns/prefer_patterns support
- `agent/best_files_tracker.py` - CarryOn fix with valuable_output_patterns
- `agent/directive_extractor.py` - Enhanced LLM prompt, added `_fix_ligand_workflow_conflict()`
- `agent/session.py` - Added fsync to save()
- `phenix_ai/log_parsers.py` - Added file existence retry in resolve_file_path()
- `knowledge/programs.yaml` - Added input_priorities for pdbtools
- `agent/generate_safety_docs.py` - Safety documentation generator
- `docs/SAFETY_CHECKS.md` - New auto-generated safety documentation
- `docs/implementation/PROGRAM_CONFIG_ROBUSTNESS.md` - New implementation plan

---

## Version 42 (January 2025)

### Testing Infrastructure

- **Converted all tests to cctbx-style**: Migrated 8 test files (300+ tests) from `unittest.TestCase` to plain functions with fail-fast behavior
  - Matches PHENIX/cctbx testing conventions
  - First assertion failure stops with full traceback
  - Simpler syntax without class wrappers

- **New test utilities module** (`tests/tst_utils.py`):
  - 20+ assert helper functions (`assert_equal`, `assert_in`, `assert_true`, etc.)
  - `run_tests_with_fail_fast()` for cctbx-style test execution
  - Supports both plain functions and TestCase classes for gradual migration

- **New testing documentation** (`docs/guides/TESTING.md`):
  - Complete guide to writing and running tests
  - Migration guide from unittest to cctbx style
  - Best practices and conventions

### Bug Fixes

- **`sanitize_for_transport()` now handles all types**: Previously converted dicts/lists to string representation. Now recursively sanitizes nested structures while preserving their types.
  - Dicts → sanitized dicts
  - Lists → sanitized lists
  - None/int/float/bool → passed through unchanged
  - Strings → sanitized (control chars, tabs, markers removed)

- **`encode_for_rest()` handles non-string input**: Added type checking to JSON-encode dicts/lists before REST encoding, preventing AttributeError on dict input.

- **`validate_directives()` preserves `file_preferences`**: Added support for boolean preferences (`prefer_anomalous`, `prefer_unmerged`, `prefer_merged`) which were previously dropped.

- **`_fix_program_name()` expanded aliases**: Added mappings for:
  - `sharpen_map`, `auto_sharpen` → `phenix.map_sharpening`
  - `build_model`, `buildmodel` → `phenix.map_to_model`

### New Directive Patterns

- **Atom type extraction**: "use selenium", "Se-SAD", "sulfur SAD" → sets `atom_type` in autosol settings
- **File preferences**: "use anomalous data", "prefer unmerged data" → sets `file_preferences`
- **Workflow preferences**: "skip autobuild", "avoid ligandfit" → adds to `skip_programs` list
- **Stop after refinement**: "stop after refinement" now works (previously required "stop after THE FIRST refinement")

### Files Changed

- `agent/transport.py` - Fixed `sanitize_for_transport()` and `encode_for_rest()`
- `agent/directive_extractor.py` - Added new patterns, fixed `validate_directives()`, expanded `_fix_program_name()`
- `tests/tst_utils.py` - New assert helpers and test runner
- `tests/tst_*.py` - Converted to cctbx-style (8 files)
- `docs/guides/TESTING.md` - New testing documentation
- `docs/README.md` - Updated testing section

---

## Version 41 (January 2025)

### Enhancements
- **Enhanced session summary**: Improved AI-generated summary to include:
  - **Key Output Files**: Now shows file name, type, and descriptive text explaining what each file contains (e.g., "Refined atomic model (X-ray)", "Structure factors with R-free flags and map coefficients")
  - **Key Metrics**: Enhanced prompt requests specific metric values and names; added extraction for Ramachandran outliers, rotamer outliers, and MolProbity score
  - Output files table now formatted with File/Type/Description columns
  - LLM prompt explicitly requests formatted metrics list with values

- **Multi-step workflow support**: Added `start_with_program` directive for handling sequences like "run polder then refine"
  - When user specifies "calculate polder map and then run refinement", the system extracts `start_with_program: phenix.polder`
  - This tells the workflow "run this program first, then continue with normal workflow"
  - Different from `after_program` which means "run this and stop"
  - Cleaner semantics than a `required_programs` list

- **Fixed directive override behavior**: Safer attempt-based strategy
  - First attempt (attempt_number=0): Honor user's directive value (respect explicit request)
  - Retry attempts (attempt_number>0): Trust LLM's interpretation (it may be correcting syntax)
  - This is safer than always trusting one or the other:
    - User's explicit request gets a fair chance first
    - If it fails, LLM can try to fix potential syntax issues
  - Example: User says `selection=solvent molecule MES 88` (invalid syntax)
    - Attempt 0: Uses user's value → fails
    - Attempt 1: LLM interprets as `selection=resname MES and resseq 88` → succeeds

- **Fixed fallback program selection**: Fallback now respects `start_with_program` directive
  - Previously, if LLM failed 3 times, fallback would pick the first valid program (often xtriage)
  - Now fallback prioritizes `start_with_program` if set by directive

### Bug Fixes
- **phenix.polder workflow integration**: Fixed issue where LLM incorrectly assumed polder needs map coefficients from prior refinement
  - Added polder to PROGRAM REFERENCE in LLM system prompt with clear documentation that it works with standard MTZ data (Fobs + R-free flags)
  - Added polder to both `refine` and `validate` phases in workflows.yaml
  - Explicit clarification: "does NOT need pre-calculated map coefficients or phases"

- **Generic PDB file categorization**: Fixed critical bug where unclassified PDB files (e.g., `1aba.pdb`) were being categorized as `search_model` instead of `model`
  - Changed `unclassified_pdb.parent_category` from `search_model` to `model` in file_categories.yaml
  - Added `*search*`, `*sculptor*`, `*chainsaw*` to excludes list to prevent search model files from being miscategorized
  - This ensures generic PDB files are treated as positioned models ready for refinement
  - Files explicitly named as search models (e.g., `search_model.pdb`, `template.pdb`) correctly go to `search_model` category
  - Previously, providing a simple PDB file would cause the workflow to think Phaser/MR was needed
  - Now the workflow correctly recognizes the model is already placed and allows refinement/validation programs

### New Tests
- Added workflow configuration tests for polder (TestPolderWorkflowConfig)
- Added LLM prompt tests for polder (TestPolderLLMPrompt)
- Added file categorization tests (TestUnclassifiedPDBCategorization)
- Tests verify:
  - Polder is in correct workflow phases
  - Prompt clarifies polder doesn't need phases
  - Generic PDB files categorize as `model` not `search_model`

### Files Changed
- `knowledge/prompts_hybrid.py` - Added phenix.polder to VALIDATION PROGRAMS section
- `knowledge/workflows.yaml` - Added phenix.polder to xray refine and validate phases
- `knowledge/file_categories.yaml` - Changed unclassified_pdb parent_category to 'model', added excludes
- `agent/session.py` - Enhanced `_get_final_output_files()` with descriptions, added `_describe_output_file()`, enhanced `_extract_final_metrics()` with more metrics, updated LLM summary prompt
- `agent/directive_extractor.py` - Added `start_with_program` extraction for multi-step workflows
- `agent/directive_validator.py` - Added attempt-based override behavior, `validate_intent()` now accepts `attempt_number`
- `agent/workflow_engine.py` - Added `start_with_program` handling in `_apply_directives()`
- `agent/graph_nodes.py` - Pass `attempt_number` to `validate_intent()`, fallback respects `start_with_program`
- `tests/tst_new_programs.py` - Added TestPolderWorkflowConfig, TestPolderLLMPrompt, TestUnclassifiedPDBCategorization
- `tests/tst_workflow_state.py` - Fixed test_dock_in_map_option to use clear search model filename
- `docs/guides/USER_DIRECTIVES.md` - Added `start_with_program` docs, attempt-based override docs, multi-step example

---

## Version 40 (January 2025)

### New Features
- **USER_REQUEST_INVALID event**: When user requests a program that's not available (e.g., "run phenix.xxx"), the agent now displays a prominent warning explaining why the request can't be fulfilled and what will run instead
- Warning is shown at QUIET verbosity level (always visible)
- Distinguishes between "unknown program" and "wrong workflow state"

### Files Changed
- `agent/event_log.py` - Added USER_REQUEST_INVALID event type
- `agent/event_formatter.py` - Added formatter for prominent warning display
- `agent/graph_nodes.py` - Emit event when user request detected as invalid

---

## Version 39 (January 2025)

### Bug Fixes
- **Event transport plumbing**: Fixed events not flowing through in two edge cases:
  1. Single-shot mode via `run_job_on_server` - events now decoded from server response
  2. API result retrieval via `get_results_as_JSON()` - events now serialized in output_files

### Files Changed
- `programs/ai_agent.py` - Added events serialization in `_build_output_files_from_history`
- `programs/ai_agent.py` - Added events decoding in `run_job_on_server`

---

## Version 38 (January 2025)

### Event System Phase 4: Display Integration
- Added `verbosity` parameter to `phenix.ai_agent` command
- Integrated EventFormatter for consistent output formatting
- Added `_display_cycle_events()` method for event rendering
- Legacy fallback when events not available

### Files Changed
- `programs/ai_agent.py` - Verbosity parameter, EventFormatter integration

---

## Version 37 (January 2025)

### Event System Phase 3: Transport Integration
- Events included in v2 API response schema
- LocalAgent and RemoteAgent parse events from responses
- Events stored in history_record for persistence

### Files Changed
- `phenix_ai/run_ai_agent.py` - Include events in response
- `phenix_ai/local_agent.py` - Parse events from response
- `agent/api_schema.py` - Updated response schema

---

## Version 36 (January 2025)

### Event System Phase 2: Decision Point Instrumentation
- All graph nodes now emit structured events
- Full LLM reasoning captured without truncation
- File selection reasons tracked

### Files Changed
- `agent/graph_nodes.py` - Event emission in perceive, plan, build nodes

---

## Version 34 (January 2025)

### Event System Phase 1: Core Infrastructure
- Created EventLog class for structured logging
- Created EventFormatter for human-readable output
- Defined 17 event types with verbosity levels
- LangGraph state compatibility (list of dicts)

### New Files
- `agent/event_log.py` - EventLog class, EventType constants
- `agent/event_formatter.py` - EventFormatter class

---

## Version 33 (January 2025)

### Cleanup and Production Hardening
- Removed deprecated state.md files
- Removed redundant backup files
- Fixed program registration after import changes
- Updated test suites for new structure

---

## Version 32 (January 2025)

### Pattern Centralization
- Moved all regex patterns to `knowledge/patterns.yaml`
- Created PatternManager for centralized access
- Updated log_parsers.py to use PatternManager

### New Files
- `knowledge/patterns.yaml` - Centralized regex patterns
- `agent/pattern_manager.py` - Pattern loading and compilation

---

## Version 31 (January 2025)

### Unified Command Builder
- Single CommandBuilder class for all programs
- Reads program definitions from YAML
- Consistent file selection across all programs
- Strategy flags and defaults from YAML

### Files Changed
- `agent/command_builder.py` - Complete rewrite

---

## Version 30 (January 2025)

### File Categorization Consolidation
- Centralized file categorization in `file_categorization.py`
- Semantic categories: model vs search_model distinction
- Categories defined in `file_categories.yaml`

### New Files
- `knowledge/file_categories.yaml`
- `agent/file_categorization.py` - Centralized categorization

---

## Version 29 (January 2025)

### BestFilesTracker
- New class to track best file of each type across cycles
- Scores based on metrics (R-free, resolution)
- R-free flag locking after first refinement

### New Files
- `agent/best_files_tracker.py`

---

## Version 28 (January 2025)

### YAML Configuration System
- Programs defined in `programs.yaml`
- Workflows defined in `workflows.yaml`
- Metrics defined in `metrics.yaml`
- Transport rules defined in `transport.yaml`

### New Files
- `knowledge/programs.yaml`
- `knowledge/workflows.yaml`
- `knowledge/metrics.yaml`
- `knowledge/transport.yaml`
- `knowledge/yaml_loader.py`

---

## Version 25-27 (December 2024)

### User Directives System
- Natural language directive parsing
- Stop conditions: "stop after X", "stop when metric < Y"
- Workflow preferences: "skip program", "prefer program"
- Four-layer stop condition checking

### New Files
- `agent/directive_extractor.py`
- `agent/directive_validator.py`
- `docs/guides/USER_DIRECTIVES.md`

---

## Earlier Versions

### Initial Development (2024)
- LangGraph pipeline architecture
- LLM integration (Claude, Gemini)
- Rules-only fallback mode
- Local and remote execution modes
- Session tracking and history
- Sanity checking system
