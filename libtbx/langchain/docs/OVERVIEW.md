# PHENIX AI Agent V2 - Technical Overview

## Introduction

The PHENIX AI Agent is an intelligent automation system for macromolecular
structure determination. It combines:

- **Goal-Directed Planning** - Decomposes structure determination into
  phases with success criteria, retreat logic, and hypothesis testing (v114)
- **LLM Decision Making** - Uses Claude/Gemini to interpret context and
  make intelligent program choices
- **Rules-Based Validation** - Enforces workflow constraints and validates
  all decisions
- **YAML Configuration** - All domain knowledge externalized to editable
  config files
- **Structured Event Logging** - Transparent decision tracking at multiple
  verbosity levels

The agent operates at two levels: a **strategic planner** (v114,
activated by `thinking_level=expert`) that produces a multi-phase
plan at session start and evaluates progress at stage gates, and
a **reactive execution engine** (v112+) that handles per-cycle
program selection, file management, error recovery, and safety
checks. The strategic planner communicates through directives —
the same interface a human user would use — so the reactive
agent's safety checks always apply.

---

## Quick Start

```bash
# X-ray structure determination (expert mode is default)
phenix.ai_agent original_files="data.mtz sequence.fa"

# Cryo-EM structure determination  
phenix.ai_agent original_files="map.mrc sequence.fa"

# With user guidance
phenix.ai_agent original_files="data.mtz model.pdb ligand.pdb" \
    project_advice="Solve the structure and fit the ligand"

# Point at a tutorial directory (auto-discovers files + README)
phenix.ai_agent input_directory=/path/to/tutorial/

# Stop after specific step
phenix.ai_agent original_files="data.mtz seq.fa" \
    project_advice="Stop after one refinement job"

# Rules-only mode (deterministic, no LLM, auto-discovers files)
phenix.ai_agent use_rules_only=True input_directory=/path/to/data/

# Without strategic planning (advanced reasoning only, no plan/gates)
phenix.ai_agent thinking_level=advanced \
    original_files="data.mtz sequence.fa"

# Control output verbosity
phenix.ai_agent verbosity=verbose original_files="data.mtz sequence.fa"

# Run a specific program
phenix.ai_agent original_files="model.pdb map.ccp4" \
    project_advice="Run phenix.map_correlations"

# Tutorial mode (set up tutorial in GUI, then open AI Agent)
# The agent auto-detects the tutorial and reads the README
phenix.ai_agent input_directory=/path/to/tutorial/
```

---

## Active Development

The v119 cluster is the active cycle; see [CHANGELOG.md
v119](CHANGELOG.md) for the per-ship breakdown.

**Scope** — Twenty-one ships stacked over v118 in three phases.
First, an operational-hardening + Phase 2A/B sub-cluster
(`v118 → H1 → H2 → H2.1 → H3 → H3b → H4 → H4.1 → H5 → H5.1 → H5.1.1
→ H6 → H6.1 → H7`): centralized LLM model defaults (H1), a
server-to-client build metadata channel (H2), skip_programs
promotion (H2.1), a startup-canary infrastructure (H3/H3b),
preprocessor telemetry markers with golden-master corpus pinning
(H4/H4.1), a uniform diagnostic-messages relay channel (H5/H5.1)
with two latent-bug cleanups (H5.1.1), the planning-suite
reliability testing framework (H6/H6.1), and Phase 2B's
scanner-first file extraction activation (H7).  Second, a
production-bug sub-cluster (`H7 → H8 → H9 → H10 → H11`) fixing
four bugs surfaced by AIAgent_62 and run_38 batch testing:
template-literal allowlist gap (H8), `predict_and_build` PHIL
scope mismatch (H9), and a paired structural + data fix for
`exclude_patterns` (H10 closes a gap where 7 of 9 selection paths
didn't honor the filter; H11 corrects three YAML patterns that
had been authored with substring semantics in mind despite the
function using word-boundary matching).  Third, three follow-up
ships (`H11 → H12 → H13 → H14`): H12 refactors H10's closure into
a module-level helper and adds a categorizer semantic-pin suite
generalizing the H11 lesson to all of `agent/file_utils.py`; H13
fixes two Ollama-provider bugs surfaced by Tom's `run_39a_ollama`
failure (URL `/v1` suffix not normalized on env-var override,
`OLLAMA_LLM_MODEL` env-var ignored by two consumer sites) and
completes the v118 §3.5 retired-model classification work with
a three-way 404 sub-categorization (RETIRED vs UNAVAILABLE vs
FAILED); H14 fixes three independent regressions surfaced by
`run_39_openai` batch analysis — a phaser false-positive from
goal phrases in `_ACTION_TABLE["solve"]` (the 1029B-sad case),
a duplicate `diagnostic_messages` relay causing `[STEP_1F]`
markers to be emitted twice (60.6% of runs), and a
`space_group` validation gap that let prose phrases through.

Total v119 cluster tests: **189 K-tests + 30 live LLM tests +
14 production-bug K-tests** (4 Bug 8, 4 Bug 9, 5 Bug 10, 1 Bug
11) **+ 8 categorizer semantic-pin K-tests** (H12) **+ 13
provider-error classification K-tests** (H13) **+ 62 H14+H14.1+H14.2
K-tests** (13 solve-keywords, 8 STEP_1F single-emit, 24 space_group
validation including 230/230 official Hermann-Mauguin symbols and
all alternative cell/origin settings, 12 simple-extractor validation
closure, 5 predict_and_build no-NCS config fix).  All ships verified
through act → review → Gemini-critique → ship → production-verify
cadence.

**v119 lessons surfaced** — Seven patterns from the cluster
that generalize:

1. **Sandbox stubs for non-obvious-semantics functions need
   semantic-pin tests.**  H10's sandbox stub for
   `matches_exclude_pattern` did substring matching while the
   real function uses word-boundary matching.  H10's K-suite
   passed in sandbox but the cycle-7 bug remained in production
   because the YAML patterns assumed substring semantics.  H11
   added `test_bug11_matches_exclude_pattern_semantics` — 10
   named cases against the real function — as the template for
   catching such divergences.  H12 then generalized the pattern
   to all of `agent/file_utils.py` via `tst_file_categorizer.py`
   (8 tests, ~75 assertions covering `classify_mtz_type`,
   `is_model_file`, `get_category_for_extension`, and the other
   public functions).  Documented in `DEVELOPER_GUIDE.md` §10.3
   ("Test adequacy").

2. **Configuration grammar must be authored to function
   semantics.**  The `exclude_patterns` YAML grammar requires
   word-boundary-correct patterns: no leading or trailing
   underscores; explicit numeric variants alongside bare patterns
   when alphanumeric suffixes are expected (`half`, `half1`,
   `half2`).  H11 embeds a `DESIGN NOTE` block at the top of
   `programs.yaml` documenting these rules so future YAML authors
   don't repeat the substring-style authoring mistake.

3. **Centralization happens in stages — data first, rules
   second.**  H1 centralized the LLM `DEFAULT_MODELS` table but
   not the env-var precedence rule over it.  Two consumer sites
   that read the table directly didn't honor `OLLAMA_LLM_MODEL`,
   while a third site DID — silently inconsistent.  Tom's
   `run_39a_ollama` exposed this when his `OLLAMA_LLM_MODEL=qwen2.5:72b`
   was ignored.  H13 added `resolve_model_for_provider()` to
   centralize the precedence rule, then routed all three consumer
   sites through it.  Forward policy: when centralizing a data
   table, also centralize the access patterns or expect the
   patterns to diverge across consumers over time.

4. **Error categories should map to operator actions, not just
   error types.**  Pre-H13, all 404s from LLM providers were
   classified as RETIRED, hint "update DEFAULT_MODELS".  Tom's
   actual failure was MODEL_UNAVAILABLE (model not pulled on
   local Ollama) — the right hint is "run `ollama pull X`".
   H13's `_classify_provider_error` splits 404 into three
   operationally-distinct sub-classes (RETIRED, UNAVAILABLE,
   FAILED) plus AUTH_FAILED for 401-class errors, with each
   tag carrying an actionable hint.  Co-occurrence rule
   ("model" word must appear within 80 chars of retirement
   phrase) defends against false-positives from edge-proxy 404
   pages mentioning "endpoint deprecated".

5. **Goal phrases and method requests live in different
   layers.**  H14 surfaced a conflation in
   `_ACTION_TABLE["solve"]`: keyword lists conflated GOAL
   phrasings (`"solve the structure"`) with METHOD requests
   (`"phaser"`, `"molecular replacement"`).  The goal phrase
   matching forced phaser into workflows on SAD datasets where
   autosol is the right method — a regression invisible until
   batch analysis showed 24 phaser-composition failures across
   5 datasets in run_39.  Forward policy: each
   `_ACTION_TABLE` keyword must represent an UNAMBIGUOUS method
   request; goal phrasings belong in `classify_intent` or in
   workflow rules that key on data type, not in method-action
   keyword lists.  The H14 investigation method itself was
   replicable: a controlled three-way log comparison
   (same README ± one line) isolated the trigger to a
   single-line README diff, converting a statistical signal
   from `scan_batch_runs.py` into a mechanistic root cause.

6. **Positive shape checks must validate against the FULL
   domain, not a hand-picked sample.**  H14's `_HM_FORM_RE`
   regex went through two review passes before reaching its
   final form.  The first draft was tested against ~15
   well-known protein space groups and all passed — but
   running it against the full 230 International Tables
   symbols showed only 71 matched (31%).  The missing 159
   included all monoclinic slash forms, all orthorhombic
   mirror/glide groups, and all cubic high-symmetry groups.
   Hand-picked test sets give false confidence; the authoritative
   list is the right benchmark.  Gemini's subsequent external
   review then caught that even the 230/230-matching regex
   rejected alternative cell/origin settings (`R3:H`, `P4/n:2`,
   `P21/c:b`) that appear in real PDB/mmCIF metadata.  Forward
   policy: when authoring a positive shape check for a
   well-defined domain (Hermann-Mauguin symbols, ISO dates,
   semver numbers), enumerate the full canonical set as a
   K-test regression guard, and explicitly consider what
   alternative forms exist beyond the canonical.  The
   `test_hm_form_accepts_all_230_space_groups` test pins
   this with the explicit 230-symbol list inline.

7. **Validators must converge from every producer path.**
   H14.1 surfaced a gap H14 missed: `directive_extractor.py`
   had TWO producers of the directives dict — the LLM-success
   path (validated) and the `extract_directives_simple`
   fallback path (NOT validated).  Tom's 2026-05-26 ollama
   production run exposed the gap when the bogus value
   `space_group=Not explicitly mentio` survived to the agent
   even with H14 installed.  The fix: rather than duplicate
   the validator inline at each producer, make the validator
   the canonical final-sanity step that ALL producers
   converge through — extract → validate, always.  Forward
   policy: when adding a validator for a class of bugs, grep
   for ALL producers of the value being validated; ensure
   each one calls the same validator before returning.
   Concurrent latent-bug audit: check that the validator's
   allow-list contains EVERY key the producers legitimately
   set (H14.1 caught a missing `start_with_program` entry in
   `VALID_STOP_CONDITIONS` this way — pre-existing latent bug
   that the converged-validator pattern would have surfaced
   as a regression in Item 1's behavior).  The lesson
   beneath the technical fix: sandbox K-tests that hit the
   validator directly aren't enough — K-tests must reproduce
   the production entry point with the same provider/failure
   modes (empty LLM response, malformed JSON, `use_rules_only=True`)
   to catch bypass paths.

The v118 cycle is complete; see [CHANGELOG.md
v118](CHANGELOG.md) for the full per-section breakdown.

**v118 scope** — Twelve layers stacked over v117.3 addressing
preprocessor resilience, LLM-shape variability across providers,
operational hardening, and bug-class elimination:
`v117.3 → A → C-prime → B → E → F → 5.1 → G → 6.1–6.7 → 8 → 9 → 10`.
Total sandbox: 204/204 tests passing.

**Categories of work** — Four themes across the twelve layers:

1. **Pipeline diagnostics** (Sections E, C-prime) — adds
   `[DIRECTIVE_EXTRACTION_FAILED]` / `[ADVICE_PREPROCESSING_FAILED]`
   stderr markers and splits the displayed directives diagnostic
   into `directives_user_intent` vs `directives_effective_runtime`.
   These layers proved essential for diagnosing Section 8 (the
   gemini-2.0-flash retirement) — the marker name itself was the
   search key.
2. **Content-shape normalization** (Sections B, 9, 10) — static
   `PHIL_NAMESPACE_TRANSLATIONS` for LLM-emitted PHIL paths;
   `PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE` for cryo-EM-vs-X-ray
   density-modification disambiguation; `_coerce_setting_value`
   for list-vs-scalar JSON shape mismatch between OpenAI and
   Google.  All three follow the same pattern: declarative table
   consumed by a generic validator, with `_corrected_from`-style
   provenance preserved through transformations.
3. **Operational hardening** (Sections A, F, G, 6.7, 8) — file-
   list preservation across preprocessor round-trips; cycle-1
   experiment_type threading + R-free auto-fill in
   `command_builder`; chromadb load resilience against protobuf
   `TypeError`; environment-readiness runtime probe in
   `tst_dependencies.py`; default-model bumps across all three
   hardcoded sites.  Section G was the first verified on both
   Mac and Linux production environments.
4. **Server-side deployment discipline** (Sections F, 8, 9, 10) —
   four server-side sections in v118 (vs zero in v117) require
   file deployment AND server process restart for the fix to
   take effect.

**Principles surfaced** — Three patterns from v118 that
generalize beyond it:

1. **Whole-tree grep before declaring a string-replacement done.**
   Section 8 rev 1 missed `api_client.py`; rev 2 needed a patch
   script; rev 3 shipped the full files.  v119 cadence convention:
   every string-replacement starts with `grep -rn 'STRING'` across
   the full langchain tree.

2. **Stub-module isolation for K-tests.**  Section F established
   the pattern, extended to G and 9/10: K-tests synthesize stub
   modules for `agent.program_registry`,
   `agent.intent_classifier`, etc.  Runs in seconds without the
   full PHENIX conda env; deterministic; no real LLM calls.

3. **Gemini-reviewed plan rev cycles.**  Sections G, 9, and 10
   each cycled through 2–4 plan revs against Gemini critique
   before implementation.  Each revision caught a real concern.
   v118.10 rev 3 in particular caught the `bool([False]) == True`
   truthiness trap that rev 2 would have shipped.

**Architectural watchpoint triggered.**  At v118.10 the project
is at 10 iterations in the preprocessor → extractor → planner →
BUILD pipeline (original threshold was 6).  v119 will pursue
operational hardening (centralized model defaults, server
`/version` endpoint, startup canary, server-stderr-to-client
diagnostic echo) and prompt consolidation; v120 may pursue
Direction A (single structured-extraction LLM call replacing the
three-LLM chain).  See
`v118_next_steps_consolidated_rev4.md` for the full v119+ plan.

### v117 cycle (preceding v118)

The v117 cycle (v117 → v117.1 → v117.2 → v117.3) reworked the
directive-extraction pipeline to make it robust against
preprocessor mangling of user stop intent:

- **v117 Step 1**: `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` —
  extractor LLM now receives raw advice alongside processed
  advice via an AUTHORITY paragraph stating raw is the source
  of truth for intent.  Eliminates the entire class of
  "preprocessor changed `Stop Condition: None` into `Stop
  Condition: None`" bugs.
- **v117.1**: grounding ∧ `stop_after_requested` interaction —
  v116.19a's grounding guardrail conflicted with v117 Step 1
  when both `after_program` and `stop_after_requested=True`
  were set; v117.1 exempts the grounding check in this case.
- **v117.2**: fill `after_program` from raw advice when LLM
  omits it — closes the C1 LLM-test gap where openai sometimes
  failed to set `after_program` even when raw advice was
  explicit.
- **v117.3**: extended stop-intent phrasing recognition —
  adds 5 imperative markers and 2 regex patterns to recognize
  "stop the workflow", "is the last step", etc.

See ARCHITECTURE.md §§9–13 for the v117 architectural changes
and CHANGELOG.md for the per-version details.

### v116.10 cycle (Pre-v117)

The v116.10 cleanup cycle addressed six bugs in advice
parsing, LLM prompting, program selection, and plan
classification, plus drift-detection machinery for the wire
contract and the client's plan-generation logic.  Five files
modified, 89 new tests across 7 new test files plus one
augmented file (`tst_contract_compliance.py`).

**Categories of work** — The fixes operate at four layers:
mechanical filters (Phase 4b strips programs whose inputs aren't
present), prompt engineering (Phase 6a reframes the LLM's
`after_program` directive), state-machine routing (Phase 6b
recognizes sessions that can't analyze), and plan-generation
classification (Phases 1, 3a, 3d fix the user-advice filter and
the `_initialize_plan_inner` standalone-programs classification).
Phase 2 (protocol hygiene) sits alongside and enforces wire-contract
invariants that the existing `tst_contract_compliance.py` suite
couldn't previously catch.

**Principles surfaced** — Three patterns from this cleanup that
generalize beyond v116.10:

1. **Filter before adding to `valid_programs`.** Phase 4b strips
   programs whose inputs aren't present before they reach the LLM,
   which is more robust than relying on the LLM to notice missing
   inputs from the prompt. Any program with hard input requirements
   should be filterable at the workflow_engine layer.

2. **Refactor means no behavior change, period.** Phase 3 was
   initially delivered as a single phase that bundled the
   de-duplication refactor with two new programs added to
   `_STANDALONE_PROGRAMS` (silently changing behavior). Self-review
   caught this. The fix was to split into Phase 3a (pure refactor,
   verified by 51 program×intent trace-equivalence checks) and
   Phase 3d (explicit behavior change, with `tst_dock_and_stop.py`
   including a 48-trace blast-radius regression guard).

3. **Drift catchers belong in the existing test suite.** Phase 2
   added `validate_contract()` to `agent/contract.py`, but its
   drift-detection power came from adding a one-line call to it
   in the existing `tst_contract_compliance.py`. The same pattern
   in Phase 3a: the consistency test imports `_ACTION_TABLE` from
   the source-of-truth module rather than maintaining a duplicate.

**Known gaps** (deferred to future work):

- Phase 3d behavior change is verified by decision-tree traces;
  no tutorial currently exercises "dock and stop with sequence +
  map". Rollback is one line if a regression appears.
- Phase 4b only filters xtriage and mtriage. A systematic
  "every program declares its inputs" audit is unfinished.
- The v116.10 prediction-only allowance in
  `_check_program_prerequisites` is mostly redundant after Phase
  6b but retained as defense in depth.

**Post-Phase-5 additions.** Three user-reported issues, one
internally-discovered fix during integration testing, and two
test suites closing planned verification gaps were patched after
the Phase 5 docs shipped:

- **CC key extraction (S20).** A successful cryo-EM workflow
  displayed the misleading "SESSION STOPPED - INCOMPLETE" banner
  because `_generate_structure_report` looked up the CC metric
  under `map_model_cc` (the form programs print in logs) instead
  of `model_map_cc` (the canonical storage key). The fix is a
  defensive multi-key lookup that accepts both forms. This adds
  a fourth principle to the cycle's collection:
  **when fixing bugs that cross the print/store boundary,
  accept both forms rather than picking the "correct" one.**

- **File encoding (S21).** A Chinese-locale Windows user
  reported a `UnicodeDecodeError` at PHENIX startup because
  Python's `open()` without explicit `encoding=` falls back to
  the system code page. An audit found **305 additional sites**
  across the codebase with the same vulnerability; all were
  patched to specify `encoding='utf-8'` explicitly. The
  regression test uses directory-scan tests that walk the
  production and test trees automatically, so new files are
  covered without test updates. This adds a fifth principle:
  **cross-platform code requires explicit encoding everywhere
  text files are opened; relying on defaults assumes a locale
  the user may not have.**

- **Ligand workflow restart (S22, Phase 6c).** The nsf-d2-ligand
  tutorial restarted its reasoning at cycle 3, with the LLM
  declaring *"As this is the first refinement step"* even after
  successful refinement and ligand fitting had completed. Two
  interacting bugs: `_detect_xray_step` returned "analyze" any
  time `xtriage_done=False` regardless of downstream progress;
  `best_files_tracker` scored `ligand_fit_output` (105) below
  `refined` (≈122) so the agent's "best model" pointer stayed on
  the unliganded refined model. Both fixes mirror existing
  in-codebase patterns: the `past_analysis` check copies
  `_detect_cryoem_step`; the inheritance fix extends the
  existing `with_ligand` machinery to also cover
  `ligand_fit_output`. This adds a sixth principle: **when
  fixing analogous bugs in parallel code paths (X-ray vs cryo-EM,
  with_ligand vs ligand_fit_output), check whether the fix
  pattern already exists nearby.** Both bugs had close analogs
  that simply hadn't been extended to cover the affected case.

- **Tier 1 follow-up tests (S23 + S24).** Two test suites that
  close verification gaps flagged in the v116.10 review:
  `tst_initialize_plan_smoke.py` (9 tests) asserts on
  `_initialize_plan_inner`'s observable side effects — directive
  rewrites, plan-generation gating — for each branch of its
  decision tree, complementing the trace-only tests in
  `tst_dock_and_stop.py`; and
  `tst_phase3d_motivating_tutorial.py` (3 tests) provides
  decision-tree, YAML format, and skip-aware session-based
  verification for the dock-and-stop motivating case. With these
  tests, all four Tier 1 items from the v116.10 follow-up plan
  are now complete (1.1 done in earlier amendment cycles; 1.2
  handled by Tom directly; 1.3 = S23; 1.4 = S24).

- **GUI tutorial override fix.** Discovered during integration
  testing of Phase 6c. When running a tutorial whose README
  mentioned a phenix program not yet supported by the AI Agent,
  the GUI blocked Run with a `Sorry` dialog even when the user
  provided their own files and advice. The fix in
  `wxGUI2/Programs/AIAgent.py`: guard the blocking dialog with
  `not has_files and not has_advice` so it only fires when the
  user is actually relying on the README. The informational
  banner stays. User-confirmed working.

### Prior cycle

The v115 cycle addressed 22 infrastructure items identified from
dual-run evaluation of 21 tutorials plus 5 bug fixes from test
suite review. 20 of 22 infrastructure items completed; I6
(unsupported programs) and I7 (tar.gz input) are deferred pending
workflow engine expansion. See [`PLAN.md`](../PLAN.md) at the
repo root for status.

### Success Metrics

| Metric | Baseline (pre-v115) | Target (v115) | Measured by |
|--------|---------------------|---------------|-------------|
| Wasted program cycles | 41% (46/111 cycles) | < 25% | `analyze_tutorial_runs.py` dual-run |
| Tutorial solve rate (readme mode) | not yet measured | > 80% | `tutorial_expectations.yaml` |
| Terminal failure pivots firing | 0% (no pivot logic) | > 90% | `tst_error_classifier.py` |
| PHIL rejections causing halt | untracked | 0 | `tst_phil_validation.py` |
| Intent misclassification | untracked | < 5% on tutorial set | `tst_intent_classifier.py` |
| "X and stop" sessions | AUTO-STOPs (Bug 1) | preserves valid_programs | `tst_user_advice_filter.py` (v116.10) |
| xtriage offered without .mtz | yes (Bug 3) | filtered | `tst_data_input_filter.py` (v116.10) |
| `_ACTION_TABLE` target drift | manual | automated | `tst_standalone_consistency.py` (v116.10) |
| `CURRENT_PROTOCOL_VERSION` drift | manual | automated | `validate_contract()` + `tst_contract_compliance.py` (v116.10) |
| Cryo-EM runs reported INCOMPLETE | yes (CC bug) | fixed | `tst_cc_key_extraction.py` (v116.10) |
| Non-UTF-8 Windows startup crash | yes | fixed | `tst_file_encoding.py` (v116.10) |
| Text-mode `open()` without `encoding=` | ungoverned | automated | `tst_file_encoding.py` directory-scan tests (v116.10) |
| Workflow state stuck at `xray_initial` post-ligandfit | yes (restart at cycle 3) | advances correctly | `tst_ligand_workflow_restart.py` Section A (v116.10 Phase 6c) |
| `ligand_fit_output` not selected as best model | yes (scored below refined) | inherits metrics, becomes best | `tst_ligand_workflow_restart.py` Section C (v116.10 Phase 6c) |
| `_initialize_plan_inner` side effects untested | only trace-tested | behavioral tests | `tst_initialize_plan_smoke.py` (v116.10 Tier 1 follow-up) |
| Phase 3d dock-and-stop tutorial verification | trace-only | tutorial integration test | `tst_phase3d_motivating_tutorial.py` (v116.10 Tier 1 follow-up) |
| GUI tutorial Run blocked despite user input | yes (bug) | guarded on user input | manual verification + scenario tracing |

---


## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              USER INTERFACE                                  │
│  phenix.ai_agent original_files="data.mtz seq.fa" project_advice="..."      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            programs/ai_agent.py                              │
│  • Parse parameters, initialize session                                      │
│  • Generate plan at session start (v114)                                    │
│  • Choose LocalAgent or RemoteAgent                                          │
│  • After each cycle: gate evaluation + plan update (v114)                   │
│  • Display results with EventFormatter                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                    ┌─────────────────┴─────────────────┐
                    ▼                                   ▼
         ┌──────────────────┐                ┌──────────────────┐
         │   LocalAgent     │                │   RemoteAgent    │
         │ (phenix_ai/)     │                │ (phenix_ai/)     │
         └────────┬─────────┘                └────────┬─────────┘
                  └─────────────────┬─────────────────┘
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       STRATEGIC PLANNER (v114)                                │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌────────────────┐  │
│  │ Plan         │  │ Gate         │  │ Structure    │  │ Explanation    │  │
│  │ Generator    │  │ Evaluator    │  │ Model        │  │ Engine         │  │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘  └──────┬─────────┘  │
│         │                 │                 │                  │            │
│         └─────── directives + advice ───────┘                  │            │
│                           │                                    │            │
│  ┌────────────────────┐   │  ┌───────────────────┐             │            │
│  │ Hypothesis         │───┘  │ Validation        │─────────────┘            │
│  │ Evaluator          │      │ History           │                          │
│  └────────────────────┘      └───────────────────┘                          │
└───────────────────────────────┼─────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         phenix_ai/run_ai_agent.py                            │
│  • Build LangGraph state                                                     │
│  • Execute graph nodes (perceive → think → plan → build → validate → output)│
│  • Retry loop with fallback on persistent validation failures               │
│  • Build response with events                                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            LangGraph Pipeline                                │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐                  │
│  │ perceive │ → │  think   │ → │   plan   │ → │  build   │                  │
│  └────┬─────┘   └──────────┘   └────┬─────┘   └────┬─────┘                  │
│       │                             │   ▲          │                         │
│       │                             │   └──────────│────── (retry < 3)      │
│       │                             │              │                         │
│       │                             │         ┌────┴──────┐                  │
│       │                             │         │ validate  │                  │
│       │                             │         └────┬──────┘                  │
│       │                             │              │                         │
│       │                             │         ┌────┴─────┐                   │
│       │                             │         │ fallback │ (retry >= 3)      │
│       │                             │         └────┬─────┘                   │
│       │                             │              │                         │
│       │                             ▼              ▼                         │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                       ┌──────────┐                                  │    │
│  │                       │  output  │ → END                            │    │
│  │                       └──────────┘                                  │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Core Components

### 1. Graph Nodes (`agent/graph_nodes.py`)

The agent uses a LangGraph pipeline with seven nodes:

| Node | Purpose | Key Actions |
|------|---------|-------------|
| **perceive** | Understand current state | Categorize files, detect workflow state, extract metrics, analyze trends |
| **think** | Expert reasoning (v113+) | Analyze logs with domain expertise; update Structure Model from validation, xtriage, phaser results (v114); generate hypothesis prompts (v114); inject strategic guidance into plan |
| **plan** | Decide next action | Call LLM or rules engine, validate against workflow, check user directives and plan-generated directives (v114) |
| **build** | Generate command | Select files using BestFilesTracker, build command string |
| **validate** | Check output | Sanity checks, red flag detection; retry up to 3× on failure |
| **fallback** | Last resort | If validate fails 3×, use mechanical/rules-based program selection |
| **output** | Format response | Package decision, events, and command into response |

The graph has conditional routing: perceive can skip to output (on red
flag abort), validate can loop back to plan (retry) or route to fallback
(max retries exceeded).

**Between cycles** (in `ai_agent.py`, not in the LangGraph pipeline):
the Gate Evaluator (v114) runs after each cycle to check phase success
criteria, advance/retreat/skip phases, update the plan, and trigger
`advice_changed` on plan revisions. Stage transition summaries are
generated by the Explanation Engine at each gate evaluation.

### 2. Workflow Engine (`agent/workflow_engine.py`)

Interprets `workflows.yaml` to determine:
- **Current phase**: X-ray has 9 phases (analyze, obtain_model, molecular_replacement, experimental_phasing, build_from_phases, refine, combine_ligand, validate, complete); cryo-EM has 9 phases (analyze, obtain_model, dock_model, check_map, optimize_map, ready_to_refine, refine, validate, complete)
- **Valid programs** for current phase
- **Transition conditions**
- **Stop criteria** (target reached, plateau detected)

### 3. Command Builder (`agent/command_builder.py`)

Unified command generation that:
- Reads program definitions from `programs.yaml`
- Selects appropriate input files using BestFilesTracker
- Applies strategy flags and defaults
- Tracks WHY each file was selected (for transparency)
- **Recovery param injection**: Recovery-sourced strategy entries (e.g.,
  `obs_labels` from ambiguous data label recovery) are appended after
  template-based command assembly, since `build_command` only emits params
  matching known `strategy_flags` keys
- **Content-based guards**: rejects small-molecule PDB files from model slots
  (`_pdb_is_small_molecule`) and protein models from ligand slots
  (`_pdb_is_protein_model`)
- **Word-boundary exclude patterns**: `matches_exclude_pattern()` prevents
  false positives (e.g., "noligand" no longer matches "ligand")
- **LLM validation**: LLM file assignments checked against slot `exclude_patterns`
  and content guards before acceptance
- **Diagnostics**: `_last_missing_slots` records which required slots could not
  be filled, enabling specific fallback error messages

### 4. Event System (`agent/event_log.py`, `agent/event_formatter.py`)

Structured logging for transparency:

| Event Type | Verbosity | Description |
|------------|-----------|-------------|
| `cycle_start` | quiet | New cycle beginning |
| `cycle_complete` | quiet | Cycle finished |
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution parsed |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `sanity_check` | normal | Red flag or warning detected |
| `program_selected` | normal | Decision with full reasoning |
| `program_modified` | normal | Program changed by rules/validation |
| `stop_decision` | normal | Whether to continue |
| `directive_applied` | normal | User directive enforced |
| `expert_assessment` | normal | Thinking agent analysis (v113) |
| `user_request_invalid` | quiet | User requested unavailable program |
| `files_selected` | verbose | File selection with reasons |
| `file_scored` | verbose | Individual file scoring detail |
| `command_built` | normal | Final command |
| `thought` | verbose | LLM chain-of-thought/reasoning traces |
| `error` | quiet | Error occurred |
| `warning` | quiet | Non-fatal warning |
| `debug` | verbose | Internal debug information |

### 5. File Categorization (`agent/graph_nodes.py`, `knowledge/file_categories.yaml`)

Semantic file classification using rules from `file_categories.yaml`:

| Category | Description | Examples |
|----------|-------------|----------|
| `data_mtz` | Reflection data (Fobs, R-free) | data.mtz, refine_001_data.mtz |
| `map_coeffs_mtz` | Map coefficients (calculated phases) | refine_001_001.mtz, denmod_map_coeffs.mtz |
| `sequence` | Sequence files | seq.fa, protein.fasta |
| `model` | **Positioned** models | refine_001.pdb, phaser_output.pdb |
| `search_model` | Templates/predictions | template.pdb, alphafold.pdb |
| `full_map` | Cryo-EM maps | map.mrc, emd_1234.ccp4 |
| `half_map` | Half-maps | half1.mrc, map_half_a.mrc |
| `ligand` | Ligand files | ligand.pdb, ATP.cif |

**Semantic distinction**: `model` = already positioned in crystal/map; `search_model` = template not yet placed.

**Post-processing content guards** validate YAML categorizer output:
- **Half-map guard**: Files in `half_map` without "half" in name → reclassified to `full_map`
- **Ligand guard** (v112.74): PDB files in `ligand_pdb` that are actually protein
  models (>150 atoms, majority ATOM records) → rescued to `model`.  Catches
  false positives from broad YAML patterns matching names like `1aba.pdb`.
- **MTZ safety net**: Cross-checks all MTZ files against `classify_mtz_type()` regex
- **Pose file exclusion** (v115.09b): LigandFit individual pose files
  (`ligand_fit_1_pose_5.pdb`) are excluded from `ligand_fit_output` category
  at categorization time. Only the final combined model (`ligand_fit_1.pdb`)
  is categorized. Additionally, `session.get_available_files()` filters
  `_pose_`/`_pose.` from the returned file list — this is the primary
  defense, catching pose files that enter through `_discover_cycle_outputs`
  glob scans.

### 6. Best Files Tracker (`agent/best_files_tracker.py`)

Tracks the best file of each type across cycles:
- Scores files based on metrics (R-free, resolution, cycle number)
- **Dual MTZ tracking**:
  - `data_mtz`: Locks after R-free flags are first generated or supplied
    (consistency for refinement; if the input MTZ already carries an R-free
    array the agent keeps it rather than regenerating — see ARCHITECTURE.md
    "R-free generate guard").  A flagless or unverifiable file later selected
    over the lock (e.g. a `data_mtz` preference for the raw input) is reconciled
    back to it — see ARCHITECTURE.md "R-free Lock Reconciliation — F1 / F2"
  - `map_coeffs_mtz`: Always prefers most recent (maps improve with refinement)
- Provides `best_files["model"]`, `best_files["data_mtz"]`, `best_files["map_coeffs_mtz"]` to CommandBuilder
- **Supplemental file discovery**: Session load (`_rebuild_best_files_from_cycles`)
  and live cycle completion (`record_result`) both call `_find_missing_outputs` to
  discover companion files (e.g., `refine_001.mtz` from `refine_001_data.mtz`) and
  evaluate them through the tracker. This ensures `map_coeffs_mtz` is populated even
  when the client only tracked a subset of output files.
- **MTZ classification**: `file_utils.classify_mtz_type()` uses regex
  `(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$` to correctly identify standard refinement
  output as map coefficients (not raw data)
- **MTZ categorization safety net** (v112.71): After both YAML and hardcoded
  categorization, `_categorize_files()` cross-checks every MTZ file against the
  authoritative `classify_mtz_type()` regex. Catches three failure modes:
  misclassified files (moved to correct category), dual-categorized files
  (removed from `data_mtz` when also in `map_coeffs_mtz`), and missing
  subcategories (added to `refine_map_coeffs` etc.). Logs `WARNING` when
  corrections are made, making future occurrences immediately diagnosable.

### 7. User Directives (`agent/directive_extractor.py`)

Parses natural language guidance from `project_advice`:
- **Stop conditions**: "stop after one refinement", "stop when R-free < 0.25"
- **Workflow preferences**: "skip autobuild", "use SAD phasing"
- **Task focus**: "focus on ligand fitting"

**Stop condition semantics:**

| Condition | Type | Where checked | Behavior |
|-----------|------|---------------|----------|
| `after_cycle` | Hard stop | PERCEIVE | Immediately stops at cycle N |
| `r_free_target` | Hard stop | PERCEIVE | Immediately stops when R-free ≤ target |
| `map_cc_target` | Hard stop | PERCEIVE | Immediately stops when map CC ≥ target |
| `after_program` | Minimum-run guarantee | PLAN | Suppresses auto-stop until target program has run, but the **LLM decides** when to actually stop |

**Why `after_program` is not a hard stop (v112.78):** The directive extractor
can only name one program in `after_program`.  For multi-goal requests like
"improve map, get symmetry and map correlation," it picks one program (e.g.,
`map_symmetry`) and goals beyond that point would be silently dropped if
`after_program` triggered an immediate stop.  By letting the LLM decide, all
goals in the user's advice are honored.

**Ligand-fit exception (v115.09b):** When ligand-fit signals are detected
("fit ATP", "fit ligand", etc.), `after_program` is cleared entirely. Ligand
fitting is a multi-step workflow (refine → ligandfit → pdbtools → refine →
polder → validate) driven by the `refine_placed_ligand` plan template (6
stages). The LLM sets `after_program` to a different program each run
(ligandfit, refine, or polder), and each choice blocks a different
intermediate step. Clearing `after_program` lets the plan gates advance
through all stages naturally.

**General `after_program` resolver (v115.10):** Replaces per-workflow
overlays with a single mechanism.  Uses `_ACTION_TABLE` (14 actions)
to detect program mentions in the user's advice, then applies rules
based on action count and explicit "stop":
- 2+ actions, no stop → clear `after_program` (plan drives)
- 2+ actions + stop → `after_program` = last mentioned
- 1 action + stop → set `after_program` if LLM missed it
- 1 action, no stop → leave as-is

Handles all workflows: MR, cryo-EM density modification, ligand
fitting, predict/build, etc.  Includes negation detection ("don't
build") and predict+build compound rule.

**Parameter requests** ("set resolution to 3", "run 4 cycles") use
a different, three-layer architecture: (1) the LLM extracts the
concept and value, (2) ``strategy_flags`` in programs.yaml maps
concepts to program-specific PHIL paths, (3) the command builder
applies them.  Unlike stop/continue (where the LLM is systematically
unreliable), the LLM handles parameter extraction well.  When a
parameter request fails, the fix is to add a ``strategy_flags``
entry in programs.yaml, not new Python code.

### 8. Safety Checks (`agent/sanity_checker.py`, `agent/directive_extractor.py`)

Multiple layers of validation to prevent errors:

| Layer | Location | Examples |
|-------|----------|----------|
| **Sanity Checks** | Pre-execution | No data for workflow, model not positioned, repeated failures |
| **Directive Validation** | Post-LLM | Invalid program names, conflicting stop conditions |
| **Workflow Validation** | State machine | Invalid stage transitions, wrong experiment type |
| **Post-Processing** | After extraction | Ligand workflow conflict resolution |

Key sanity issues that trigger abort:
- `no_data_for_workflow` - Missing data_mtz (X-ray) or map (cryo-EM)
- `search_model_not_positioned` - Trying to refine before MR/docking
- `no_model_for_refine` - No model available for refinement
- `repeated_failures` - Same error 3+ times

See DEVELOPER_GUIDE.md §6 (Safety Checks) for the complete list.

### 9. RAG Pipeline (`rag/`, `utils/query.py`, `analysis/analyzer.py`)

The agent uses Retrieval-Augmented Generation to ground LLM responses in
PHENIX documentation. The pipeline has three stages:

1. **Retrieval**: A Chroma vector store (built from PHENIX docs via
   `run_build_db.py`) returns the top 20 candidate documents using
   embedding similarity.

2. **Reranking**: FlashRank, a local cross-encoder model
   (`ms-marco-MiniLM-L-12-v2`, ~34MB), reranks the 20 candidates and
   selects the top 8 most relevant. This runs on CPU with no API key —
   the `flashrank` package must be installed where analysis runs locally.

3. **Generation**: The reranked documents are passed as context to the
   LLM (Google or OpenAI) which generates the final response.

The reranking retriever is created in `rag/retriever.py` via
`create_reranking_retriever()` and used by both log analysis
(`analysis/analyzer.py`) and documentation queries (`utils/query.py`).
The reranker itself is `PhenixFlashrankCompressor`, a local
`BaseDocumentCompressor` that calls the `flashrank` package directly (no
`langchain-community` dependency).

### 10. Thinking Levels and the Strategic Planner

The `thinking_level` parameter controls how much intelligence the
agent applies to each cycle. It is a PHIL choice parameter with four
values: `none`, `basic`, `advanced`, and `expert` (default). Each
level is additive — higher levels include everything from lower levels
plus additional capabilities.

A separate parameter, `use_rules_only=True`, is orthogonal to
`thinking_level`. It replaces the LLM call in the PLAN node with
deterministic first-valid-program selection. Everything else — PERCEIVE,
THINK, BUILD, VALIDATE, safety checks, error recovery — runs identically
regardless of `use_rules_only`. The two axes are independent: you can
run `thinking_level=advanced use_rules_only=True` to get structural
validation and expert KB without any LLM calls in PLAN.

#### Capability Matrix

| Capability | `none` | `basic` | `advanced` | `expert` (default) |
|---|---|---|---|---|
| **Inside the LangGraph pipeline (per-cycle)** | | | | |
| PERCEIVE: file categorization, workflow state, metrics, sanity checks | ✓ | ✓ | ✓ | ✓ |
| THINK: LLM call with log analysis | — | ✓ | ✓ | ✓ |
| THINK: strategy memory (persists across cycles) | — | ✓ | ✓ | ✓ |
| THINK: structural validation (Ramachandran, clashscore, rotamers) | — | — | ✓ | ✓ |
| THINK: Structure Model (accumulated structural knowledge) | — | — | ✓ | ✓ |
| THINK: Validation History (per-cycle snapshots, trend analysis) | — | — | ✓ | ✓ |
| THINK: file metadata tracking | — | — | ✓ | ✓ |
| THINK: expert knowledge base (56 rules, IDF-weighted) | — | — | ✓ | ✓ |
| THINK: hypothesis extraction from LLM response | — | — | ✓ | ✓ |
| PLAN: LLM program selection (or rules-only fallback) | ✓ | ✓ | ✓ | ✓ |
| BUILD: deterministic command assembly | ✓ | ✓ | ✓ | ✓ |
| VALIDATE: workflow, file, and duplicate checks | ✓ | ✓ | ✓ | ✓ |
| **Outside the graph (between cycles, in ai_agent.py)** | | | | |
| Plan generation at session start (17 templates) | — | — | — | ✓ |
| Plan-to-directives merging before each cycle | — | — | — | ✓ |
| Gate evaluation after each cycle (advance/retreat/skip/stop) | — | — | — | ✓ |
| Strategy blacklisting on retreat | — | — | — | ✓ |
| Hypothesis lifecycle management (confirm/refute/abandon) | — | — | — | ✓ |
| Cycle commentary (template-based, no LLM) | — | — | — | ✓ |
| Plan enforcement (suppress premature STOP) | — | — | — | ✓ |
| Model placement gate | — | — | — | ✓ |
| **At session end** | | | | |
| Structure report (HTML with SVG trajectory) | — | — | — | ✓ |
| Session summary JSON | — | — | — | ✓ |
| Final/stopped report (template-based, no LLM) | — | — | — | ✓ |
| Failure diagnosis (LLM-generated, on terminal error) | — | — | — | ✓ |

#### Per-Level Behavior

**`none`** — Minimum viable agent. The THINK node is a complete
pass-through. The pipeline is effectively PERCEIVE → PLAN → BUILD →
VALIDATE → OUTPUT. PLAN still calls the LLM for program selection
(unless `use_rules_only=True`), so there is one LLM call per cycle.
No between-cycle operations beyond basic stop checks.

**`basic`** — Adds a second LLM call in the THINK node. The thinking
context includes log sections extracted from the last program's output
(program-specific keyword extraction), current metrics, R-free trend
across history, a brief history summary (last 5 cycles), recent
failures, and strategy memory. The LLM produces an assessment with
an action, confidence, analysis, and guidance. If guidance is produced,
it is prepended to `user_advice` as `[Expert assessment] ...` so the
PLAN node's LLM sees it when selecting a program. Strategy memory
is updated and persisted. The `should_think()` gate applies: THINK
only engages after strategic programs (analysis, model-building,
refinement, ligand), after failures, or when R-free is stalled. It
skips on the first cycle and on routine steps.

**`advanced`** — Adds four subsystems to the THINK node's context:
(A) structural validation via headless `run_validation()` on the
current best model, producing Ramachandran, clashscore, rotamer
analysis, and model contents; (B) a Structure Model that accumulates
cross-cycle knowledge from ground-truth validation results — data
characteristics, model state, R-free trajectory with annotations,
and a strategy blacklist — plus a Validation History with per-cycle
snapshots; (C) file metadata entries for validated models; and
(D) expert knowledge base queries with IDF-weighted tag matching
against 56 crystallographic rules. The LLM call receives all of this
richer context. The between-cycle loop in `ai_agent.py` is unchanged
from `none`/`basic` — no plan, no gate evaluation, no reports.

**`expert`** (default) — Everything from `advanced`, plus the entire
Strategic Planner layer activates in `ai_agent.py`. Inside the graph,
the THINK node runs identically to `advanced` (the value `expert`
is mapped to `advanced` for graph execution, since all planning
operations live outside the graph). The additions are all in the
between-cycle loop, described in Sections 10a–10h below.

#### How `expert` Maps to `advanced` in the Graph

The LangGraph pipeline only knows three thinking levels: `none`,
`basic`, and `advanced`. When the user sets `thinking_level=expert`,
`create_initial_state()` maps it to `advanced` for graph execution.
This is deliberate: the graph nodes (PERCEIVE, THINK, PLAN, BUILD,
VALIDATE) behave identically at `advanced` and `expert`. Everything
that distinguishes `expert` from `advanced` — plan generation, gate
evaluation, hypothesis lifecycle, stage summaries, reports — lives
in the outer loop in `ai_agent.py`, gated on the original
`thinking_level == "expert"` check.

The subsections below (10a–10h) describe the components of the
strategic planner in detail.

Note: The Structure Model (10a) and Validation History are updated
at both `advanced` and `expert` levels inside the THINK node. The
remaining components (10b–10h: Plan Generator, Gate Evaluator,
Hypothesis Engine, Explanation Engine, Model Placement Gate, Plan
Enforcement, Display Data Model, HTML Report) only run at `expert`.

#### 10a. Structure Model (`agent/structure_model.py`)

Maintains a running understanding of the specific structure being solved.
Updated every cycle from validation results (ground truth, not LLM
reasoning). Persists across cycles and session resume.

| Tracking Area | Content | Source |
|---------------|---------|--------|
| Data characteristics | Resolution, space group, twinning, anomalous | xtriage log |
| Model state | Chains, ligands, waters, problem regions | validation_inspector |
| Progress | R-free trajectory with annotations | log analysis + validation |
| Strategy blacklist | Tried-and-failed strategies | gate evaluator retreats |
| Hypotheses | Proposed/active/confirmed/refuted | hypothesis evaluator |

Key methods: `update_from_validation()`, `update_from_xtriage()`,
`update_from_phaser()`, `get_summary(detail_level=)`,
`get_current_problems()`, `blacklist_strategy()`, `is_blacklisted()`.

The `ValidationHistory` (`agent/validation_history.py`) stores per-cycle
validation snapshots with `get_metric_series()` for trend analysis and
`get_phase_start_metrics()` for the monotonic progress gate.

#### 10b. Plan Generator (`agent/plan_generator.py`)

Produces a multi-phase strategy at session start by selecting and
customizing pre-defined templates from `knowledge/plan_templates.yaml`.

Twelve templates cover common scenarios: `mr_refine`,
`mr_refine_ligand`, `mr_refine_lowres`, `mr_refine_highres`,
`mr_refine_twinned`, `predict_refine`, `predict_refine_ligand`,
`mr_sad`, `sad_phasing`, `sad_phasing_ligand`, `cryoem_refine`,
`cryoem_refine_ligand`. Template selection is deterministic
(experiment type + available files + resolution + anomalous
atoms). The LLM only customizes parameters within template
bounds.

Plans are represented as `StructurePlan` objects (`knowledge/plan_schema.py`)
containing a list of `StageDef` phases, each with programs, success
criteria, gate conditions, fallbacks, and skip conditions.

The plan communicates with the reactive agent through:
- `plan_to_directives()` — translates current phase to `prefer_programs`,
  `after_program`, `program_settings`
- `compute_hash()` — Strategy Hash fingerprint; changes trigger
  `advice_changed` in the reactive agent

#### 10c. Gate Evaluator (`agent/gate_evaluator.py`)

After each cycle, evaluates phase progress against the plan's success
criteria. Purely deterministic (no LLM). Returns one of: continue,
advance, retreat, fallback, skip, or stop.

**Success hysteresis**: a 1.5% buffer prevents oscillation on noisy
metrics (e.g., R-free 0.349 → 0.351 → 0.349 doesn't trigger repeated
advance/continue decisions).

**Retreat logic** with five anti-oscillation safeguards:

| Safeguard | Prevention |
|-----------|------------|
| Strategy Blacklist | Never re-try a failed strategy |
| Retreat counter | Max 2 retreats per phase |
| Monotonic progress gate | Only retreat if worse than phase start |
| Retreat cooldown | 2+ cycles between retreats |
| Retreat depth limit | Max 1 phase backwards |

#### 10d. Hypothesis Engine (`agent/hypothesis_evaluator.py`)

Structured hypothesis formation and testing, integrated into the THINK
node. Enforces a single active hypothesis budget (one testing/pending
at a time) to prevent confounded multi-variable experiments.

Lifecycle: proposed → testing → pending → confirmed/refuted/abandoned.
`test_cycles_remaining` provides verification latency (prevents premature
refutation on unstabilized models). Confirmed hypotheses are re-validated
each cycle and demoted if evidence weakens (B-factor drift, RSCC drop).

#### 10e. Explanation Engine (`knowledge/explanation_prompts.py`)

Produces crystallographer-level commentary at three detail levels:

| Function | When | LLM? |
|----------|------|------|
| `generate_cycle_commentary()` | Every cycle | No (template) |
| `generate_stage_summary()` | Stage transitions | No (template) |
| `generate_final_report()` | Session completion | No (template from Structure Model) |
| `generate_stopped_report()` | Early stop | No (template from Structure Model) |

The failure diagnosis (`ai_failure_diagnosis.html`) is separate from the
Explanation Engine — it is LLM-generated and produced by `failure_diagnoser.py`
when the agent detects a terminal error (see ARCHITECTURE.md).

Note: `generate_stage_summary()` no longer includes the "Next:"
line — that information is shown separately in the GUI's stage
transition block ("→ ADVANCING TO: ...").

#### 10f. Model Placement Gate (v114.1, hardened v114.2)

Prevents destructive program selection on models that already fit the
data. This is a single mechanism that addresses the class of bugs where
the LLM runs phaser/autosol on a pre-solved structure, producing a
worse result than simple refinement.

Detection (in `_detect_model_placement`, `programs/ai_agent.py`):
- `model_vs_data` with CC > 0.3 → placed
- `model_vs_data` with R-free < 0.50 from result text → placed (v114.2)
- `refine` with R-free < 0.50 → placed
- `real_space_refine` with map_cc > 0.3 → placed
- `real_space_refine` symmetry mismatch → needs_dock (v114.2)

When placement is confirmed:
1. `session.data["model_is_placed"]` — locked, survives resume
2. `valid_programs` filtered: phaser, autosol, predict_and_build removed
3. Plan fast-forward: MR/phasing stages (and model_rebuilding if
   R-free < 0.35) marked "skipped" (⊘)
4. Conflict warning if user advice mentions MR/phaser
5. Dock keywords ("dock", "docking") in advice cancel "refine" →
   "placed" assumption (v114.2)

#### 10f-b. Plan Enforcement (v114.2)

Prevents premature STOP when the strategic plan has pending stages.
Two enforcement points in `agent/graph_nodes.py`:

1. **AUTO-STOP suppression**: `plan_has_pending_stages` in
   session_info blocks metrics-based auto-stop.
2. **LLM STOP override**: When the LLM returns `program=STOP`
   but the plan has pending stages, the build node selects the
   best program from `prefer_programs`.

Contract field: `plan_has_pending_stages` (bool, default False, v4).

#### 10g. Display Data Model (`agent/display_data_model.py`)

Unified data provider for Results tab, Progress tab, and HTML report.
Eliminates duplication across three display views:

- `outcome_status`: determined / stopped / incomplete
- `outcome_message`: one-line human-readable summary
- `final_metrics`: best R-free/CC from structure model or cycle scan
- `rfree_trajectory` / `cc_trajectory`: for SVG charts
- `timeline`: compact cycle entries
- `format_cycle_compact()`: one-line-per-cycle formatting

#### 10h. Session Reports

At session completion, the agent produces up to four output files.
All reports except the failure diagnosis are template-based with no
LLM call.

**Structure Determination Report** (`structure_determination_report.txt`):
Text report generated by `generate_final_report()` or
`generate_stopped_report()` in `knowledge/explanation_prompts.py`.
Uses data from the Structure Model (ground truth, not LLM reasoning).
Contains: data characteristics (resolution, space group, twinning,
MR/phasing metrics), final model stats (R-free, R-work, clashscore,
Ramachandran, chains, ligands, waters), plan stage outcomes, and
(for stopped sessions) a strategies-that-failed section with metrics
at retreat, outstanding issues, and actionable recommendations.

**HTML Structure Report** (`structure_report.html`):
Template-based HTML report generated by `generate_html_report()` in
`knowledge/html_report_template.py`. Uses the `DisplayDataModel` as
its data source. Includes:
- Title with outcome status icon (✅ determined, ⚠ stopped, ● incomplete)
- Final metrics table (R-free, R-work, CC, clashscore, geometry)
- Inline SVG R-free/CC trajectory chart with retreat markers
- Workflow stage table with status icons (✓ complete, ✗ failed, — skipped)
- Per-cycle timeline (cycle, program, status, key metric)
- Output file paths (model, map)
The GUI provides an "Open Structure Report" button that opens this file.

**Session Summary JSON** (`session_summary.json`):
Machine-readable summary for the evaluation harness and downstream
tools. Generated by `_write_session_summary_json()` in `ai_agent.py`.
Contains: outcome (complete/stopped/incomplete/aborted/reactive),
stop_reason, stop_reason_code, final metrics, data characteristics,
stage outcomes with cycles and metrics, hypothesis outcomes (if any),
strategy blacklist, and R-free/CC metric trajectory.

**Failure Diagnosis** (`ai_failure_diagnosis.html`):
LLM-generated diagnosis produced by `_diagnose_terminal_failure()` in
`ai_agent.py` when the agent detects a terminal program error (e.g.,
crystal symmetry mismatch, missing data columns, PHIL parse failure).
The `DiagnosisDetector` identifies the error type from
`diagnosable_errors.yaml`, reads the failing program's log tail, and
sends both to the LLM (via the server's `failure_diagnosis` analysis
mode). Error types include crystal symmetry mismatch, model outside
map, SHELX not installed, unknown PHIL parameter, polymer special
position, and missing crystal symmetry (v115.09b — covers xtriage
with .sca data lacking unit cell info). The response is structured as three sections: "What went wrong,"
"Most likely cause," and "How to fix it." The HTML report opens
automatically in the user's browser. If the LLM is unavailable
(rules-only mode or API failure), a deterministic fallback uses the
YAML hint text instead. When a failure diagnosis is produced, the
normal structure report is skipped — the diagnosis is the user's
primary output.

---

## Configuration Files

### programs.yaml

Defines each PHENIX program the agent can run:

```yaml
phenix.refine:
  description: "Crystallographic refinement"
  category: refinement
  experiment_types: [xray]
  
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
        priority_patterns: [refine]  # Prefer refined models
      data_mtz:
        extensions: [.mtz, .sca, .hkl, .sdf]
        flag: ""
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
      - pattern: "*_data.mtz"
        type: data_mtz
      - pattern: "*_refine_*.mtz"
        type: map_coeffs_mtz
    metrics:
      - r_free
      - r_work
      - bonds_rmsd
      - angles_rmsd
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-free"
    r_work:
      pattern: 'R-work\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-work"
    bonds_rmsd:
      pattern: '[Bb]onds?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Bonds RMSD"
    angles_rmsd:
      pattern: '[Aa]ngles?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Angles RMSD"
```

**Future: structured results.** The `log_parsing` regex approach
works but is fragile — log format changes silently break extraction.
Newer PHENIX programs built on `ProgramTemplate` expose a
`results_as_json()` method that returns metrics as structured JSON.
As programs adopt this, the agent can read JSON results directly
instead of parsing logs, with `log_parsing` as a fallback for older
programs. See ARCHITECTURE.md "Potential improvements" for the
migration plan.

### workflows.yaml

Defines workflow state machines:

```yaml
xray:
  phases:
    analyze:
      description: "Analyze data quality"
      programs:
        - phenix.xtriage
      transitions:
        on_complete: obtain_model
    
    obtain_model:
      description: "Get initial model"
      programs:
        - program: phenix.predict_and_build
          preferred: true
          conditions:
            - has: sequence
            - not_done: predict_full  # Don't re-run after full workflow completes
        - program: phenix.phaser
          conditions:
            - has: search_model
            - not_done: phaser  # MR should only run once
        - program: phenix.autosol
          conditions:
            - has: sequence
            - has: anomalous        # Requires anomalous signal detected by xtriage
            - not_done: autosol
      transitions:
        on_complete: refine
        if_predict_only: molecular_replacement  # Stepwise mode
    
    # MR-SAD: after phaser places model, use anomalous signal with autosol
    experimental_phasing:
      description: "MR-SAD phasing with placed model"
      programs:
        - program: phenix.autosol
          conditions:
            - not_done: autosol
      transitions:
        on_complete: build_from_phases
    
    refine:
      description: "Improve model"
      programs:
        - program: phenix.refine
        - program: phenix.autobuild
          conditions:
            - r_free: "> autobuild_threshold"
            - has: sequence
            - not_done: autobuild
        - program: phenix.autobuild_denmod
          conditions:
            - has: ligand_file
            - has: sequence
            - not_done: autobuild_denmod
            - refine_count: "> 0"    # Must refine first
          hint: "Run density modification before ligand fitting"
        - program: phenix.ligandfit
          conditions:
            - has: ligand_file
            - not_done: ligandfit
            - r_free: "< 0.35"
            - refine_count: "> 0"    # Must refine first
        - program: phenix.polder
          conditions:
            - has: model
            - has: data_mtz
      repeat:
        max_cycles: 4
        until:
          any:
            - r_free: "< target_r_free"
            - condition: plateau
              cycles: 2
              threshold: 0.005
      transitions:
        on_ligandfit: combine_ligand
        on_target_reached: validate
        on_plateau: validate
        on_max_cycles: validate

  targets:
    r_free:
      default: 0.25
      by_resolution:
        - range: [0, 1.5]
          value: 0.20
        - range: [1.5, 2.5]
          value: 0.25
        - range: [2.5, 3.5]
          value: 0.30
```

#### Program Execution Controls

The workflow uses two mechanisms to prevent programs from running repeatedly:

**1. `not_done` conditions** (in workflows.yaml)

Programs can specify `not_done: <flag>` to prevent re-runs:

| Flag | Program | Meaning |
|------|---------|---------|
| `predict_full` | predict_and_build (X-ray) | Full workflow (prediction+MR+building) completed |
| `predict` | predict_and_build (cryo-EM) | Any prediction completed |
| `process_predicted_model` | process_predicted_model | Model processed for MR |
| `phaser` | phaser | Molecular replacement completed |
| `dock` | dock_in_map | Model docked in map |
| `autobuild` | autobuild | Model building completed |
| `autobuild_denmod` | autobuild_denmod | Density modification completed |
| `autosol` | autosol | Experimental phasing completed |
| `ligandfit` | ligandfit | Ligand fitted |
| `resolve_cryo_em` | resolve_cryo_em | Map optimization completed |
| `map_sharpening` | map_sharpening | Map sharpening completed |
| `map_to_model` | map_to_model | De novo model building completed |
| `map_symmetry` | map_symmetry | Map symmetry analysis completed |

**2. `done_tracking` blocks** (in programs.yaml)

Each program's `done_tracking` block defines its workflow done flag and tracking strategy:

- **`strategy: "set_flag"`** (default) — sets a boolean done flag on success. Most programs use this.
- **`strategy: "run_once"`** — sets the done flag AND filters the program from the valid list after first successful run. Used by `phenix.xtriage`, `phenix.mtriage`, and `phenix.map_symmetry`.
- **`strategy: "count"`** — sets the done flag AND increments a counter (e.g., `refine_count`). Used by `phenix.refine`, `phenix.real_space_refine`, and `phenix.phaser`.

All detection is driven by `history_detection.markers` (substring matching). Programs like `phenix.refine` use `exclude_markers: ["real_space"]` to prevent false matches with `phenix.real_space_refine`. The only program requiring Python-only tracking is `phenix.predict_and_build`, which cascades flags across programs.

**Programs that run multiple times** (intentionally):
- `phenix.refine` / `phenix.real_space_refine` - Iterative refinement
- `phenix.molprobity` - Validation after each cycle
- `phenix.polder` - Can run for different sites

### metrics.yaml

Quality thresholds and display configuration:

```yaml
metrics:
  r_free:
    display_name: "R-free"
    format: "{value:.4f}"
    direction: lower_is_better
    thresholds:
      excellent: 0.20
      good: 0.25
      acceptable: 0.30
    
  map_cc:
    display_name: "Map CC"
    format: "{value:.3f}"
    direction: higher_is_better
    thresholds:
      excellent: 0.80
      good: 0.75
      acceptable: 0.70
```

---

## Data Flow

### Single Cycle

```
1. INPUT
   • available_files: [data.mtz, seq.fa, refine_001.pdb]
   • history: [{cycle: 1, program: "xtriage", ...}]
   • log_text: "R-free = 0.28..."
   • user_advice: "continue refinement"

2. PERCEIVE
   • Categorize files → {mtz: [data.mtz], model: [refine_001.pdb]}
   • Detect state → "refine" phase, valid: [refine, ligandfit, molprobity]
   • Extract metrics → {r_free: 0.28, r_work: 0.24}
   • Analyze trend → "improving, no plateau"
   • Emit: STATE_DETECTED, METRICS_EXTRACTED, METRICS_TREND
   
3. PLAN
   • Check directives → no stop condition met
   • Call LLM or rules → phenix.refine
   • Validate → ✓ in valid_programs
   • Emit: PROGRAM_SELECTED (with full reasoning)
   
4. BUILD
   • Select model: refine_001.pdb (from best_files)
   • Select mtz: data.mtz (rfree_locked)
   • Build: "phenix.refine refine_001.pdb data.mtz output.prefix=refine_002"
   • Emit: FILES_SELECTED, COMMAND_BUILT

5. OUTPUT
   • program: "phenix.refine"
   • command: "phenix.refine ..."
   • events: [{type: "state_detected", ...}, ...]
   • stop: false
```

---

## Execution Modes

The agent's behavior is controlled by two independent axes:
`thinking_level` (how much intelligence per cycle) and
`use_rules_only` (LLM vs deterministic program selection).
The combinations below are the most common configurations.

### Goal-Directed Mode (v114, default)

```bash
phenix.ai_agent original_files="data.mtz seq.fa"
```

- Default mode (`thinking_level=expert`)
- Multi-phase plan generated at session start from templates
- Gate evaluation after each cycle (advance/retreat/skip/stop)
- Model placement gate: detects when model fits data and
  suppresses destructive programs (phaser, autosol)
- Hypothesis testing with verification latency
- Per-cycle expert assessment and stage transition summaries
- Two LLM calls per cycle when THINK engages (expert reasoning + program selection)
- `structure_report.html`, `structure_determination_report.txt`,
  and `session_summary.json` at completion

### Standard Mode (LLM, no planning)

```bash
phenix.ai_agent thinking_level=advanced original_files="data.mtz seq.fa"
```

- `thinking_level=advanced`: THINK node runs with full context
  (structural validation, expert KB, Structure Model) but no
  strategic planning layer
- LLM makes program selection decisions in PLAN
- Two LLM calls per cycle when THINK engages
- No plan generation, gate evaluation, or between-cycle reports
- Useful when the planning overhead is unnecessary or when
  debugging the reactive engine in isolation

### Minimal LLM Mode

```bash
phenix.ai_agent thinking_level=none original_files="data.mtz seq.fa"
```

- THINK node is a pass-through (no expert reasoning)
- One LLM call per cycle (PLAN only)
- Lowest LLM cost; suitable for simple workflows or constrained
  API budgets

### Rules-Only Mode

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz seq.fa"
```

- Deterministic program selection (first valid program from
  workflow engine) — no LLM call in PLAN
- `use_rules_only` only affects PLAN; `thinking_level` still
  controls the THINK node independently. At the default
  `thinking_level=expert`, THINK will attempt its LLM call.
  To suppress all LLM calls, also set `thinking_level=none`.
- Faster, fully reproducible for testing (with `thinking_level=none`)
- Auto-discovers files from input_directory

### Dry-Run Mode

```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic
```

- Simulated execution with predefined outcomes
- For testing workflow logic without running PHENIX programs

---

## Error Handling

### Automatic Error Recovery

Recognized error patterns trigger automatic retry with corrected parameters:
- **Ambiguous data labels**: MTZ with multiple arrays → selects anomalous or merged
  based on workflow context, injects `obs_labels` via recovery param injection
- **Ambiguous experimental phases**: Selects HL coefficients based on context
- **Loop guard** (v112.74): If a recovery strategy already exists for the file,
  re-triggering is skipped to prevent infinite retry loops

See ARCHITECTURE.md "Automatic Error Recovery" for implementation details.

### Parameter Blacklisting (`bad_inject_params`)

When a PHENIX program fails due to an injected parameter, the parameter is
blacklisted so `inject_user_params` never re-injects it. Recognized error patterns:
- **Unknown parameter**: "Unknown command line parameter definition: FOO"
- **No such parameter**: "No such parameter: FOO"
- **Boolean type mismatch** (v112.75): "True or False value expected,
  scope.path.param=value found" — blacklists the full PHIL path and all
  components ≥ 6 characters (catches `wavelength` when PHIL resolves it to
  `autosol.wavelength.added_wavelength`)

**Static (proactive) blacklist (v120):** the above is *reactive* (learned only
after a failure). `AgentSession._STATIC_BAD_INJECT_PARAMS` seeds known-invalid
pairs so they are skipped on the FIRST command — e.g. `ncs_file` for
`phenix.resolve_cryo_em`. `get_all_bad_inject_params()` merges static ∪ learned
and is what the client transmits to the server. See ARCHITECTURE §45.

### User Request Invalid

When user requests an unavailable program:

```
============================================================
  WARNING: Requested program not available
============================================================
  You requested: phenix.refine
  Reason: Not valid in current workflow state 'xray_initial'
  Running instead: phenix.xtriage
  Available programs: phenix.xtriage
  Suggestion: This program requires different conditions.
============================================================
```

The agent:
1. Detects that user mentioned the program in their advice
2. Explains WHY it's not available
3. Suggests what will run instead
4. Always shown (QUIET verbosity level)

### Sanity Checks

Red flags that indicate problems:
- **Experiment type changed** mid-workflow → abort
- **R-free spike** (increased > 0.05) → warning
- **No model** available for refinement → abort
- **Resolution unknown** before refinement → warning

### Server Error Propagation (v112.78, extended v114.2)

Fatal server errors (e.g., daily API usage limit) are raised as `Sorry` in
`rest/__init__.py`.  `RemoteAgent` has a dedicated `except Sorry: raise` handler
before its generic `except Exception` to ensure these propagate cleanly through
to the GUI instead of being silently swallowed as a None result.

**v114.2 extension:** `RemoteAgent._send_request` also re-raises Sorry when
`server_result.success` is False and the `server_message` contains fatal
keywords ("quota", "API key", "rate limit", "authentication", "permission
denied", "cannot continue"). The same check applies to the parsed JSON
response's `error` field. Previously, these errors were logged and returned as
None, causing the agent to silently end with "No command generated".

### Cross-Platform (Windows) Considerations (v112.78)

The agent runs on macOS, Linux, and Windows. Key platform-specific handling:

| Area | Unix/macOS | Windows |
|------|-----------|---------|
| Process tree kill | `psutil` recursive walk → `SIGTERM`; fallback `os.killpg` | `taskkill /F /T /PID` |
| Abort detection | `return_code < 0` (signal) + STOPWIZARD file | STOPWIZARD file only (taskkill returns positive codes) |
| Subprocess GUI | No special handling needed | `CREATE_NO_WINDOW` creationflag prevents console flash |
| Path separators | Forward slash native | Backslash native; normalize to `/` before marker matching |
| Path quoting | `shlex.quote()` (POSIX single-quotes) | Double-quote with escaped inner `"` |
| File encoding | UTF-8 default on modern systems | Locale-dependent; explicit `encoding='utf-8'` on all `open()` |
| `os.path.relpath` | Always works | `ValueError` across drives; caught with `try/except` |
| JSON transport | No issues | v115.05: removed `json_str.replace('\\t', ' ')` from `transport.py` and `api_client.py` — it corrupted paths like `C:\tutorials\test.mtz` (the `\t` matched as tab escape). Tab handling is done correctly by `text_as_simple_string()`. |

---

## Session Management

The agent provides two mechanisms for inspecting and modifying an existing session
**without running new crystallographic cycles**.

### Viewing a Session

```bash
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=basic
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=detailed
```

`basic` prints a one-line-per-cycle summary table (program, R-free, result).
`detailed` prints full reasoning and command for every cycle.
Both modes populate `self.result` identically to a normal run so GUI calls
(`get_results()`, `get_results_as_JSON()`) work without special cases.

`restart_mode=resume` is automatically set when either session management
parameter is active — no manual flag required.

### Removing Cycles

```bash
phenix.ai_agent log_directory=AIAgent_run1 remove_last_n=2
```

Removes the last N cycles from the session, clears the stale AI summary,
rebuilds `active_files.json` and `best_files` from remaining history,
and saves. Useful for pruning a failed run before re-running.

### Extending a Completed Workflow with New Advice

When a workflow has fully completed, resuming with new `project_advice`
triggers follow-up programs via the Q1 mechanism:

```bash
phenix.ai_agent \
    log_directory=AIAgent_run1 \
    restart_mode=resume \
    project_advice="also run polder on chain B residue 100"
```

1. New advice hash detected → `advice_changed=True`
2. PERCEIVE steps `complete` phase back to `validate` (adds polder, molprobity, etc.)
3. PLAN suppresses AUTO-STOP for one cycle
4. LLM acts on the new advice; after success, normal termination resumes

See [USER_GUIDE.md §11](USER_GUIDE.md#11-directives-reference)
and [ARCHITECTURE.md](ARCHITECTURE.md) for details.

---

## Testing

### Test Suites (55+ files, 50+ in runner)

**Standalone (no PHENIX required):**
- API Schema, Best Files Tracker, Transport, State Serialization
- Command Builder, File Categorization, File Utils
- Session Summary, Session Directives, Session Tools, Audit Fix Regressions
- Advice Preprocessing, Directive Extractor, Directive Validator, Directives Integration
- Event System, Metric Patterns, Pattern Manager
- Program Registration, Summary Display, New Programs
- Error Analyzer, Decision Flow, Phaser Multimodel
- History Analysis, Docs Tools, YAML Tools
- Thinking Defense, Strategy Memory, Log Extractor, Thinking Agent (v113)
- Structure Model, Validation History, Plan Generator, Gate Evaluator,
  Hypothesis Evaluator (v114 — 251 tests)

**PHENIX-dependent:**
- Workflow State, YAML Config, Sanity Checker
- Metrics Analyzer, Dry Run, Integration, Directives Integration

**Systematic Testing Framework (v115.08 — 10 phases, ~14s):**

Bottom-up testing of system boundaries where bugs hide. Exercises real
production code paths end-to-end. Designed after code reviews found 5
critical bugs that existing unit tests missed.

| Suite | Phase | Tests | What It Tests |
|-------|-------|-------|---------------|
| S0 | Static Audit | 5 | Parse check, bare except scan, import fallbacks |
| S1 | Contract Gaps | 128 | AST coverage map — 4 modules, boundary gap detection |
| S2 | Path Consistency | 10 | YAML vs hardcoded categorization diff (whitelisted) |
| S3 | Session Round-Trip | 28 | JSON symmetry + AgentSession save/load/pipeline |
| S4 | History Flags | 8 | Flag writer/reader consistency, dead flag audit |
| S5 | Category-Consumer | 14 | input_priorities + fallback_categories alignment |
| S6 | Routing Simulation | 32 | 3-cycle routing through detect_step + get_valid_programs |
| S7 | Command Building | 15 | CommandBuilder.build() with real file combinations |
| S8 | Error Classification | 7 | 3 classifiers × 30+ patterns, overlap + severity |
| S9 | LLM Perturbation | 17 | Filename/program/parameter/truncation/empty resilience |

Phases S6–S8 are skipped in `--quick` mode. All phases produce machine-readable
findings in `findings/`. See `docs/PHASE_REVIEW_REPORT.md` for review details.

**Additional test files (not in run_all_tests.py):**
- tst_template.py (template builder), tst_utils.py (assert helpers)

### Running Tests

```bash
python3 tests/run_all_tests.py        # All tests
python3 tests/run_all_tests.py --quick  # Standalone only (skips S6-S8)
python3 tests/tst_event_system.py    # Single suite
python3 tests/tst_phase7_routing_simulation.py  # Single phase
```

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v120 | **New providers (Portkey + Claude) + `FORCE_NO_AI_SERVER`** (three phases P1–P3): **P1** adds the `FORCE_NO_AI_SERVER=1` environment override in `programs/ai_agent.py::run()` that forces `communication.run_on_server=False` (local execution) without changing the shipped PHIL default — no effect on `predict_and_build`/`predict_model` run separately; **P2** centralizes `SUPPORTED_PROVIDERS` as a single source of truth in `core/llm.py` (was duplicated literals in `graph_nodes.py` and `ai_agent.py` that could drift), with both consumers importing it; **P3** adds **portkey** (the Portkey gateway fronting Azure OpenAI, via the Portkey SDK with `provider="azure-openai"`; env `PORTKEY_AZURE_API_KEY`/`PORTKEY_BASE_URL`/`PHENIX_PORTKEY_MODEL`) and **anthropic** (Claude; env `ANTHROPIC_API_KEY`) as provider choices, wired through every routing point: `core/llm.py` factory + `get_expensive_llm` + new helpers (`portkey_langchain_config`, `sanitize_llm_kwargs`, `_delegate_embeddings_for_nonnative`), `api_client._call_portkey_llm`, `directive_extractor` (`_call_llm` + fallback), `rate_limit_handler.get_portkey_handler`, `graph_nodes.validate_provider` + handler site, `thinking_agent`, the embeddings/RAG DB path (`run_utils`, `summarizer`, `run_query_docs`, `rebuild_ai_database`, `update_ai_database`), PHIL enums, and `install_ai_tools.csh`. `sanitize_llm_kwargs` resolves the asymmetric-args trap (portkey drops `temperature` + remaps to `max_completion_tokens`; anthropic keeps both, which `ChatAnthropic` requires). Embeddings degrade gracefully for non-native providers (delegate or `None`, never raise on construction). Provider travels in request `settings` so local and server run identical code. Derived from a user proof-of-concept patch, hardened for production (no PHIL-default flip, clear `ValueError` over bare `KeyError`, dedicated rate handler). 36 new K-tests across three suites. |
| v119 | **Operational hardening + production-bug cluster** (twenty-eight ships H1–H18.2): H1 centralizes LLM model defaults in `core/llm.DEFAULTS`; H2 adds the server→client `agent_build` metadata channel with version+fingerprint injection; H3/H3b add a startup canary that detects wrong-build deploys and unhealthy LLM environments; H4/H4.1 add the `[STEP_1F]` preprocessing-telemetry marker pinned by a golden-master corpus; H5/H5.1 establish a uniform `diagnostic_messages` relay surfacing server-side stderr markers to the client transparently; H5.1.1 cleans up two latent v118.9 bugs; H6/H6.1 add the planning-suite reliability testing framework; H7 activates Phase 2B scanner-first file extraction; H8–H11 fix four production bugs (template-literal allowlist, predict_and_build PHIL scope, `exclude_patterns` application gap, YAML pattern semantics); H12 refactors the exclude_patterns helper to module level with semantic-pin tests; H13 fixes two Ollama-provider bugs (`OLLAMA_BASE_URL` path construction, `OLLAMA_LLM_MODEL` honored only in `core/llm`) and refines retired-model 404 classification into three actionable sub-categories; H14/H14.1/H14.2 close three independent regressions from `run_39_openai` batch (phaser false-positive in `_ACTION_TABLE["solve"]`, duplicate `[STEP_1F]` emission, space_group `Hermann-Mauguin`-shape validation) and a one-line `programs.yaml` fix (`predict_and_build` does not accept `ncs_spec`); **H15** Tom's bromodomain resume failure (run 144) — `reopen_stages_for_directives` targeted single-stage reopen with O(1) blast radius via the LATEST-stage rule; **H16/H16.1** obs_labels auto-fill for multi-array MTZ via new `agent/mtz_inspector` module + `auto_fill_obs_labels` invariant in `programs.yaml` (closes 88 TIER-1 "Multiple equally suitable arrays" failures across run_25/run_39); **H17/H17.1** `phenix.autobuild` PHIB-required recovery — YAML `recoverable_errors.yaml` entry + new `strip_parameter` resolution in `error_analyzer` + executor-side strip_flags wiring in `programs/ai_agent.py` (Tom's lysozyme-MRSAD cycle-5 reactive recovery); **H18** file-based experiment-type detection — new public helper `agent/file_utils.infer_experiment_type_from_files()` used as PRIMARY signal in BOTH `_apply_experiment_type_program_reprints` AND `_resolve_after_program` (files-win on conflict with text; `[DIRECTIVE_CORRECTION]` marker enriched with `source=files\|text` and `OVERRIDDEN` annotations; Pitfall 1 `[DIRECTIVE_CORRECTION_MIXED]` telemetry for accumulated-files drift; fixes AF_7mjs "density modify and stop" regression where terse advice + cryo-EM files silently produced `phenix.autobuild_denmod` instead of `phenix.resolve_cryo_em`); **H18.1** PHIL declaration deploy-gap hotfix — adds the missing `original_files_for_directives` PHIL declaration to `programs/ai_agent.py`'s `master_params` string (H18 had added it only to `programs/ai_analysis.py`, but `directive_params = copy.deepcopy(self.params)` copies the agent's params which use the agent's schema — production crashed with `AttributeError`, directive extraction silently returned `{}`, and the agent ran the default plan template through to predict_and_build); new `tst_h18_1_phil_roundtrip.py` (6 K-tests) exercises the full PHIL parse → extract → deep-copy → assign path the production code uses, catching the deploy-gap class of bug at sandbox time; **H18.2** third `_resolve_after_program` callsite missed in H18 audit — when the LLM emits `stop_after_requested=True` with NO `after_program` field, a v117.2 fallback at `directive_extractor.py:783` fires to fill in the missing program by parsing raw advice, but pre-H18.2 it called `_resolve_after_program` WITHOUT `original_files` (the H18 audit grepped for sites billed as "experiment-type detection" and missed this site billed as "fill in missing field"); the resolver defaulted `_exp="xray"` via text-only fallback and mapped `denmod` → `phenix.autobuild_denmod` even for cryo-EM file inventories; H18 site 2 (`_apply_workflow_intent_fallback`) couldn't fix it because preprocessed advice's "Stop Condition: None" pushed the resolver into the "leave as-is" branch; H18.2 adds `original_files=original_files` to the v117.2 callsite (one line) plus 3 new K-tests (K21/K22/K23) extending `tst_density_modify_experiment_type.py` to 23 tests including the AF_7mjs production reproduction (K21 mocks the LLM with the exact directives shape Tom's production LLM produced and asserts the full pipeline produces `phenix.resolve_cryo_em`). Total cluster: **249+ K-tests added**. All ships verified through act → review → Gemini-critique → ship cadence. |
| v118 | **Preprocessor resilience + operational hardening** (12 layers stacked over v117.3): Section A file-list preservation (UNION text+context, `_ensure_file_list_in_processed_advice`); Section C-prime diagnostic split (`directives_user_intent` vs `directives_effective_runtime`); Section B PHIL namespace healing (`PHIL_NAMESPACE_TRANSLATIONS` static table); Section E LLM-failure diagnostic markers (`[DIRECTIVE_EXTRACTION_FAILED]`, `[ADVICE_PREPROCESSING_FAILED]` stderr); Section F **server-side** BUILD `experiment_type` threading + R-free auto-fill in `command_builder._select_files` (first v118 section touching server); Section 5.1 PHIL `.help` hyphen/semicolon hotfix; Section G optional dep resilience (`except Exception` not `except ImportError` for chromadb/protobuf chain) + 6.7 environment-readiness runtime probe `tst_dependencies.py` (first v118 verified Mac + Linux); Section 8 model bump `gemini-2.0-flash` → `gemini-2.5-flash-lite` (THREE hardcoded sites in `core/llm.py`, `agent/api_client.py`, `agent/directive_extractor.py` — first v118 requiring server restart for deploy); Section 9 cryo-EM vs X-ray density-mod canonicalization via `PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE` table + `_apply_experiment_type_program_reprints` validator + `[DIRECTIVE_CORRECTION]` log marker + `_corrected_from` sidecar (Tom's "density modify and stop" bug); Section 10 list-to-string coercion via `_coerce_setting_value` helper (single-element unpack type-preserving; multi-element space-join for str-typed fields) applied at TWO call sites in `validate_directives` — fixes Google's `additional_atom_types=["S"]` producing `"['S']"` AND latent silent-drop for list-shaped `after_program`. 9 new test suites, 204/204 sandbox tests pass. Architectural watchpoint triggered (10 iterations); v119 will pursue operational hardening + prompt consolidation. |
| v117.3 | **Extended stop-intent phrasing recognition**: 5 new entries in `_IMPERATIVE_STOP_MARKERS` ("stop the workflow", "immediately after", etc.), 2 new regex patterns in `_POSITIVE_STOP_AFTER_PATTERNS` ("\\bstop\\s+the\\s+workflow\\b", "\\bis\\s+the\\s+last\\s+step\\b"), 4 new phrasings in `DIRECTIVE_EXTRACTION_PROMPT` schema docs. Closes the C1 LLM-test `explicit_stop_after_phaser` failure where openai's extractor produced `after_program=phenix.phaser` but failed to set `stop_after_requested=True`. Asymmetric placement: imperative markers are 300-char window-bounded so safer there than as global regex. Bare addition to existing data structures — no new pipeline steps or state fields. |
| v117.2 | **Fill `after_program` from raw advice when LLM omits it**: in `extract_directives()`, after existing layers run, if `stop_after_requested=True` but `after_program` is missing, scan raw advice for a program name with `_POSITIVE_STOP_AFTER_PATTERNS` and fill. Closes the C1 LLM-test gap where openai sometimes failed to set `after_program` even when raw advice was explicit. |
| v117.1 | **Grounding ∧ `stop_after_requested` interaction fix**: v116.19a's grounding guardrail conflicted with v117 Step 1 — when LLM set BOTH `after_program` and `stop_after_requested=True`, grounding's Failure 2 dropped the directive. v117.1 exempts the grounding check when `stop_after_requested=True` (the user explicitly requested stop; the after_program is the stop target, not a fabricated program assertion). |
| v117 | **Extraction reliability — raw-advice-authoritative directive extraction**: new `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` with AUTHORITY paragraph stating raw advice is the source of truth for intent. Extractor LLM now receives raw advice alongside processed advice. Eliminates the entire class of "preprocessor changed `Stop Condition: None`" bugs that motivated the cycle. Falls back gracefully to single-input form when raw == processed. |
| v116.10 | **Cleanup cycle**: 6 bugs in advice parsing, LLM prompting, program selection, plan classification + drift detection. Phase 4b filter (strip programs whose inputs aren't present), Phase 6a prompt reframe (after_program), Phase 6b state-machine routing (can't-analyze sessions), Phases 1/3a/3d plan classification (`_initialize_plan_inner` standalone-programs), Phase 2 wire-contract drift catcher. S20-S24 post-Phase-5: CC key extraction (`map_model_cc` vs `model_map_cc`); file encoding (305 sites + directory-scan test); ligand workflow restart (`_detect_xray_step` past_analysis check + `with_ligand` extension to `ligand_fit_output`). 89 new tests across 7 new files + augmented `tst_contract_compliance.py`. Three principles: filter before adding to `valid_programs`; refactor means no behavior change; drift catchers belong in the existing test suite. |
| v115.10 | **General `after_program` resolver**: Replaces per-workflow overlays (ligand-fit clearing, denmod stop/clear, continuation_indicators, downstream_tasks, multi_program_patterns) with a single `_ACTION_TABLE`-based mechanism. 14 actions with keywords mapped to xray/cryoem programs. Rules: multiple actions + stop → after_program = last; multiple + no stop → clear; single + stop → set. Features: word-boundary matching, action-specific negation detection ("don't build"), predict+build compound rule, cryo-EM experiment type inference. 33 old regex patterns removed, replaced by keyword-based action detection. 20 unit tests (91 assertions). 1 file modified (`directive_extractor.py`). |
| v115.09b | **GUI Fixes + Ligand Workflow + Production Bugs**: (1) `explicit_program` done-flag guard prevents LLM-driven program loops (bgal_denmod, hipip-refine). (2) Ligand-fitting workflow: `after_program` clearing for multi-step workflows; `combine_ligand` guard forces pdbtools-only; post-ligandfit exemption defers `after_program_done` during combine/refine; `model_is_placed=True` for ligand-fit signals. (3) Pose file exclusion at three levels: `session.get_available_files()` basename filter (primary), `workflow_state.py` categorization, `programs.yaml` `exclude_patterns`. (4) `phaser_sad.atom_type` interception → `additional_atom_types` conversion in autosol. (5) Polder `requires_resolution` invariant with `auto_fill_resolution` + `xray_data.high_resolution` strategy flag. (6) Preprocessing stop override with `_has_explicit_stop` regex. (7) GUI auto-discovery skip for user-selected files. (8) `missing_crystal_symmetry` diagnosable error for xtriage with .sca data. 8 files modified. |
| v115.09 | **Tutorial Routing Fixes**: (1) Cryo-EM `past_analysis` gate: added `map_sharpening_done`, `map_symmetry_done`, `has_optimized_full_map` — unblocks bgal_denmod, apoferritin_denmod, ion_channel_denmod after map_sharpening; regex broadened to match actual output filenames. (2) `.sca-only` data detection in `perceive()` — deferred, needs reactive approach. (3) Validation-only routing: `wants_validation_only` directive extracted via LLM prompt + rules-based fallback + post-LLM overlay → validation shortcut in `_detect_xray_step` → `validate_existing` plan template; PDB scan limit 500→2000 in `_is_valid_file` (3dnd.pdb has 546 header lines); `has_phased_data_mtz` added. (4) MR-SAD routing: `force_mr` flag when `use_mr_sad` + model not categorized as search_model → phaser offered; MR-SAD guard updated. Critical deployment fix: rules-based intent patterns moved to shared `_apply_workflow_intent_fallback()` called as post-LLM overlay. 7 files modified. |
| v115.08 | **Phased File Detection + Systematic Testing Framework**: 4 critical fixes for phased file detection (content-based iotbx+ASCII heuristic replacing filename markers; shared post-processing; category exclusivity; conditional data_mtz removal). 1 additional bug fix (B1: last_program missing from build_context). `[GATE]` diagnostic logging for routing debugging. 10-phase systematic testing framework (S0–S9): static audit, contract gap coverage map, YAML/hardcoded path consistency, AgentSession round-trip, history flag consistency, category-consumer alignment, 32-tutorial routing simulation, command building, error classification, LLM perturbation resilience. 12 new test files, 42 unit tests + ~260 phase checks. Framework reviewed across multiple rounds — 54 issues found and fixed in test scripts (5 can't-fail gates, 3 tautological assertions, 5 wrong best_files keys, etc.). |
| v115.07 | **Run 15b Bug Fixes — Phase 3**: 4 bugs from 371-run analysis. Numeric coercion (`_safe_float` at 6 sites). Two-tier half-map detection. PHIL blocked params for resolve_cryo_em. Terminal diagnosis for unknown chemical elements. Reference model restraint support (hierarchical prefix whitelist, path resolution, strategy rewrites). |
| v115.05 | **Guard Fixes + Polder + Templates**: (1) `_is_at_target` hopeless R-free (> 0.50) now requires `autobuild_done` — prevents premature stop on incomplete models (p9-SAD fix). (2) `_is_at_target` clashscore path requires `refine_count >= 1` (X-ray only; cryo-EM path is theoretical-only gap). (3) Bug F `obtain_model` routing requires `autobuild_done` — gives autobuild a chance before concluding MR is wrong. (4) Negligible-anomalous guard: removes autosol from `valid_programs` when measurability < 0.05 and `has_anomalous=False`. (5) `wants_polder` context flag + polder override fires despite "solve" keyword in advice. (6) `refine_placed_polder` template added (17 templates total). (7) Early rebuild gate: `r_free > 0.50 after 1 cycles → try_rebuilding`; `gate_evaluator` now advances to `model_rebuilding` stage. (8) Phaser copies injection reads `log_analysis["n_copies"]` same-cycle (not 1-cycle delay). (9) `_anti_ligand_patterns` excludes `no_ligand` from ligand file classification; orphaned PDB files promoted to `model` (fixes user-supplied input PDBs like `1aba.pdb` not being recognized). (10) Unregistered `explicit_program` downgrades to warning (no Sorry). (11) `_preprocessing_programs` / `_needs_plan_programs` ensure polder etc. get full plans. (12) Failed programs skip output file tracking. (13) Metrics-based report selection overrides INCOMPLETE status when R-free/CC targets are met. 15 fix-verification tests (`tst_fix_verification.py`). (14) Ligand PDB plan selection: `_build_context()` now checks ligand name hints before setting `has_search_model`, and removes the `len(pdb_files) >= 2` guard — fixes AF_bromodomain_ligand tutorial selecting `mr_refine` (no ligandfit) instead of `predict_refine_ligand`. (15) Windows transport fix: removed `json_str.replace('\\t', ' ')` from `transport.py` and `api_client.py` — the replace corrupted Windows file paths containing `\t` sequences (e.g. `C:\tutorials\test.mtz`), causing "Failed to parse request JSON" on Windows clients. |
| v115 | **Infrastructure Audit + Failure Recovery**: Dual-run evaluation framework across 21 tutorials (41% baseline cycle waste identified). Intent classifier (4-way: solve/solve_constrained/task/tutorial). Tiered error recovery (`error_classifier.py`). PHIL strategy validation. Thinking agent context forwarding. Session bug fixes (resolution contract, intent low-confidence guard, pipeline classification, skip_validation stop trigger). |
| v114.1 | **Model Placement Gate + Display + Evaluation Harness**: Default `thinking_level` changed from `advanced` to `expert`. **Placement gate**: detects when model fits data (model_vs_data CC > 0.3 or refine R-free < 0.50) and locks `model_is_placed` in session — suppresses phaser/autosol/predict_and_build, fast-forwards plan past MR/phasing stages, logs conflict warning when user advice contradicts. **Display**: DisplayDataModel unified data layer for Results/Progress/HTML; HTML structure report with SVG trajectory chart; "Open Structure Report" button in GUI; expert assessments now stored in session JSON; DDM scans all cycles for best metrics. **Templates**: `mr_sad` requires explicit MR-SAD intent (`wants_mr_sad`); predict_and_build in MR stage programs; polder moved to post-ligandfit phase only with `has: ligand_fit` YAML condition; SAD templates require sequence. **Files**: auto-discover from `input_directory`; HETATM ligand detection in input PDBs; ligand PDB filename hints. **GUI**: restart_mode as plain wx.Choice (survives session management reset); stage display "cycle X, up to Y". **Safety**: sanity check threshold 3→4; recent failures injected into THINK prompt; STOP not counted as cycle. **Testing**: 57 scenario tracer tests (PG1-PG5 placement gate, L1-L10 mock LLM, C1-C3 cycle counting); tutorial run analyzer for 5 modes. 31 files modified. |
| v114 | **Goal-Directed Agent** (`thinking_level=expert`): Strategic planner layer with Structure Model (running structural knowledge + strategy blacklist), Plan Generator (12 templates, strategy hash), Gate Evaluator (success hysteresis, 5 anti-oscillation safeguards), Hypothesis Engine (single-budget, verification latency, re-validation), Explanation Engine (cycle/phase/final commentary), GUI stage display (plan header, transition blocks, per-cycle stage context), `session_summary.json` output. 12 new files, 13 modified, 251 tests. Reactive agent unchanged; `advanced` (default) continues to work without planning overhead. |
| v113.10 | **Thinking Level + Validation + Expert KB**: Replaces boolean `use_thinking_agent` with `thinking_level` parameter (`none`/`basic`/`advanced`; v114 adds `expert`). Advanced mode adds structural validation (Ramachandran, clashscore, rotamers, model contents), 56-entry expert knowledge base with IDF-weighted tag matching, file metadata tracking, R-free trend display, and user-facing Expert Assessment block in event formatter. 7 new files, 25 tests. Backward compatible. |
| v113 | **Thinking Agent**: Optional expert crystallographer reasoning node (THINK) between PERCEIVE and PLAN. Second LLM call analyzes program logs with domain expertise, injects strategic guidance via user_advice enrichment. Per-program keyword extraction (xtriage, phaser, autosol, autobuild, refine), priority-ordered sections within character budget. Strategy memory persists across cycles via session_info. GUI checkbox + `[Expert]` display in progress panel. 4 new modules, 4 new test files, 103 thinking-related tests. |
| v112.78 | **GUI mode map_coeffs_mtz + daily usage Sorry + after_program fix + Windows compat**: GUI mode `_record_command_result` and `_track_output_files` used `os.getcwd()` which pointed to parent after CWD restore — now accept explicit `working_dir`; `rest/__init__.py` raises `Sorry` on `daily_usage_reached` and `RemoteAgent` re-raises it; `after_program` changed from hard stop to minimum-run guarantee; Windows: `_filter_intermediate_files` normalizes backslash paths, `Popen` uses `CREATE_NO_WINDOW`, session JSON uses explicit UTF-8 encoding |
| v112.77 | **Autobuild rebuild_in_place**: Rule D stripped `rebuild_in_place=False` because it wasn't in strategy_flags; added `rebuild_in_place`, `n_cycle_build_max`, `maps_only` to autobuild allowlist; recovery hint for sequence mismatch errors |
| v112.76 | **Catch-all injection blacklist + deterministic atom_type**: heavier-atom-wins rule swaps `atom_type`/`mad_ha_add_list` when primary has lower Z (27-element table); catch-all streak tracker blacklists injected params after 2 consecutive same-error failures (`return_injected` kwarg, `_update_inject_fail_streak`); recovery retries excluded |
| v112.75 | **Autosol/autobuild process bugs**: strategy-flag alias awareness in `inject_user_params` (wavelength→lambda dedup); `bad_inject_params` learning expanded to PHIL boolean-type errors; autosol atom_type/mad_ha_add_list same-value dedup; `_is_program_already_done` extended to non-count programs (prevents `_apply_directives` re-adding completed autosol from program_settings); improved atom_type hint in programs.yaml |
| v112.74 | **Xtriage recovery + ligand misclassification**: recovery param injection survives command builder and probe-only sanitizer; ligand-as-model misclassification guard; obs_labels error recovery loop guard |
| v112.70 | **Ligandfit file selection**: fixed refine MTZ classification regex (3 locations); word-boundary `exclude_patterns`; content-based PDB guards for model/ligand slots; protein-in-ligand-slot rejection; refinement CIF exclusion; `inject_user_params` bare-key validation; supplemental file discovery on session load and live path; fallback diagnostics (per-program missing slots); duplicate detection respects different input files |
| v112.31 | **Session management**: `display_and_stop` / `remove_last_n` populate `self.result`; `get_results()` safe before `run()`; `restart_mode` auto-set; **Q1**: resuming with new advice after workflow completion steps back from `complete` to `validate` phase, enabling follow-up programs (polder etc.) |
| v112 | **Steps table metrics**: cycle metrics as primary source; benign warning metrics extraction; ligand typing fix; case-sensitive pattern fix; autobuild_denmod detection; YAML log_parsing for 8 programs |
| v111 | **Summary output fixes**: predict_and_build R-free extraction; ligandfit output in final file list; fallback cycle status check fix |
| v110 | **Stepwise mode**: automation_path controls predict_and_build behavior; fallback program tracking; autobuild scoring equals refined; best files in summary fix |
| v40 | Fixes 12-21: Ligandfit MTZ exclusion, stop condition on failed runs, summary display, predict_and_build resolution handling; USER_REQUEST_INVALID event |
| v39 | Event system plumbing fixes for single-shot mode |
| v38 | Event system Phase 4: display integration |
| v36-37 | Event system Phases 2-3: instrumentation and transport |
| v34 | Event system Phase 1: EventLog and EventFormatter |
| v30-33 | YAML centralization, BestFilesTracker, CommandBuilder unification |

---

## Automation Modes

The agent supports two automation modes controlled by `maximum_automation`:

### Automated Mode (default)

```bash
phenix.ai_agent maximum_automation=True original_files="data.mtz sequence.fa"
```

- `predict_and_build` runs the complete workflow (prediction → MR → building)
- Fewer checkpoints, faster end-to-end processing
- Best for well-understood datasets

### Stepwise Mode

```bash
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"
```

- `predict_and_build` stops after prediction only (`stop_after_predict=True`)
- User can inspect predicted model before proceeding
- Workflow continues: `process_predicted_model` → `phaser` → `refine`
- Best for troubleshooting or when intermediate inspection is needed

The `automation_path` is set in workflow_state and propagated to all decision-making components to ensure consistent behavior throughout the pipeline.

---

## Dependencies

The AI agent requires the following Python packages beyond the standard PHENIX
installation. Install via `phenix.python -m pip install <package>` or via the
`install_ai_tools.csh` script.

| Package | Purpose | Required on |
|---------|---------|-------------|
| `langchain-core` | LLM orchestration core | Server and local |
| `langchain-google-genai` | Google Gemini LLM provider | Server and local |
| `langchain-openai` | OpenAI LLM provider | Server and local |
| `langchain-anthropic`, `anthropic` | Anthropic (Claude) LLM provider | Server and local |
| `portkey-ai` | Portkey gateway SDK (Azure-OpenAI upstream); optional — a manual fallback path works without it | Server and local |
| `langchain-chroma` | Chroma vector store for document retrieval | Server and local |
| `flashrank` | Local cross-encoder reranking (no API key needed), wrapped by `PhenixFlashrankCompressor` | Server, or local if `run_on_server=False` |
| `pypdf`, `beautifulsoup4` | Local document loaders (PDF and HTML) for building the docs DB | Where the docs DB is built |
| `markdown-it-py` | HTML rendering of analysis output | Server and local |

**Note:** `flashrank` downloads its model (~34MB) automatically on first use.
Reranking runs entirely locally via `PhenixFlashrankCompressor` (in
`rag/retriever.py`), which calls `flashrank` directly — there is no dependency on
`langchain-community` (the loaders and reranker were localized to remove it).

### LLM Providers and Environment Variables

The `communication.provider` PHIL parameter selects the LLM backend.
Supported values: `ollama` (local), `google` (default), `openai`,
`anthropic` (Claude), `portkey` (Portkey gateway fronting Azure OpenAI).
Each provider reads its credentials from the environment at call time
(only the selected provider's keys are checked):

| Provider | Required environment variables | Optional |
|----------|--------------------------------|----------|
| `google` | `GOOGLE_API_KEY` | — |
| `openai` | `OPENAI_API_KEY` | — |
| `anthropic` | `ANTHROPIC_API_KEY` | `ANTHROPIC_LLM_MODEL` (model override) |
| `ollama` | (reachable Ollama server) | `OLLAMA_BASE_URL`, `OLLAMA_LLM_MODEL`, `OLLAMA_EMBED_MODEL` |
| `portkey` | `PORTKEY_AZURE_API_KEY`, `PORTKEY_BASE_URL` | `PHENIX_PORTKEY_MODEL` (upstream model, default `gpt-5`) |

**`FORCE_NO_AI_SERVER`** — set to `1` to force local execution
(`communication.run_on_server=False`) regardless of the PHIL default.
Useful for running entirely on one machine (e.g. the Portkey local-only
deployment).  It applies to both `phenix.ai_agent` and `phenix.ai_analysis`
(all analysis_modes); it does not affect `predict_and_build` / `predict_model`
when those are run separately.  The flag is **absolute**: it never contacts
the Phenix server.  If local execution is impossible — `standard` analysis
mode needs a local RAG database and none exists for the chosen provider — the
job errors out with actionable guidance (build the database with
`phenix.rebuild_ai_database <provider>`, use an LLM-only mode, or unset the
flag) rather than silently falling back to the server.  The LLM-only modes
(directive_extraction, advice_preprocessing, failure_diagnosis, agent_session)
never need the database, so they always run local under the flag.  Behaviour
is otherwise identical to server mode; only the locus of execution changes.

**Embeddings note:** `anthropic` has no native embeddings endpoint and
`portkey` requires an embeddings-capable deployment.  Chat works for
both regardless; the standalone `phenix.ai_analysis` RAG tool
(`standard` mode) is the only path that needs a working embeddings
backend.  For `anthropic`, embeddings transparently delegate to OpenAI
or Google when a key is available.  `anthropic` is intentionally not
offered by the database-build tools (`rebuild_ai_database` /
`update_ai_database`) since it cannot produce embeddings.

**LLM-unavailable notice:** if the requested provider cannot be reached, the
agent falls back to RULES-BASED program selection (deterministic workflow
rules, no LLM) rather than stopping outright.  This used to be silent.  It is
now surfaced clearly, in both CLI and GUI and on both local and server runs:
- **Per cycle**, each cycle that falls back prints a `NOTICE: LLM PROVIDER
  UNAVAILABLE` line naming the provider and stating that rules-based selection
  is in use for that cycle.  (Shown every affected cycle, so a long run never
  *looks* like the LLM quietly recovered.)
- **At the end of the run**, an `IMPORTANT: THIS RUN DID NOT USE THE LLM`
  banner is printed, the Results-tab summary is prefixed with the same warning,
  and the result object carries an `llm_ever_unavailable` flag.

This is observability only — it does not change which program the agent
selects, only how visibly the rules-fallback is reported.  (Setting
`use_rules_only=True` explicitly is a separate, already-announced mode; the
notice covers the *unintended* fallback when a provider is down.)

**Note:** The `langchain-classic` package is **not required**. The agent implements
document chain and compression retriever functionality directly using `langchain-core`
base classes, avoiding the deprecated `langchain.chains` and `langchain.retrievers`
modules.

---


---

## Command-Line Tools

### Configuration Management

```bash
# List all YAML configuration files
python3 agent/yaml_tools.py list

# Validate all YAML files for syntax and structural errors
python3 agent/yaml_tools.py validate

# Validate a specific file
python3 agent/yaml_tools.py validate programs.yaml

# Display formatted contents of a YAML file
python3 agent/yaml_tools.py display programs

# Compare two YAML files or directories
python3 agent/yaml_tools.py compare programs.yaml programs_backup.yaml

# Show overview of all configuration
python3 agent/yaml_tools.py summary

# Show all defined terms in the configuration system
python3 agent/yaml_tools.py terms
python3 agent/yaml_tools.py terms --detail full  # With cross-references
```

### Session Management

```bash
# Show current session summary
python3 agent/session_tools.py --show

# Show detailed session info (files, metrics, reasoning for each cycle)
python3 agent/session_tools.py --show --detailed

# Remove last N cycles from session
python3 agent/session_tools.py --remove-last 2

# Reset entire session
python3 agent/session_tools.py --reset

# Dry-run (show what would be done without saving)
python3 agent/session_tools.py --remove-last 3 --dry-run

# Use a specific session directory
python3 agent/session_tools.py --dir /path/to/session --show
```

### Documentation Generation

```bash
# Generate safety checks documentation
python3 agent/generate_safety_docs.py  # output reviewed in DEVELOPER_GUIDE.md §6
```

### RAG Documentation Database

The agent uses a Retrieval-Augmented Generation pipeline to ground LLM
responses in PHENIX documentation. Documents are stored in a Chroma vector
database and retrieved with FlashRank cross-encoder reranking (top 8 of 20
candidates). See OVERVIEW.md §9 for architecture details.

```bash
# Build vector database from PHENIX documentation
python3 run_build_db.py

# Inspect database contents
python3 run_inspect_db.py

# Query the documentation
python3 run_query_docs.py "How do I set up SAD phasing in phenix.autosol?"
```

### Program Validation

```bash
# Validate a specific program's configuration completeness
python3 agent/program_validator.py phenix.polder

# Validate all programs
python3 agent/program_validator.py --all

# List all configured programs
python3 agent/program_validator.py --list
```

### Pattern Management

```bash
# Validate all metric extraction patterns and run tests
python3 agent/pattern_manager.py
```

### Directive Validation

```bash
# Run directive validator self-test (checks program availability detection)
python3 agent/directive_validator.py
```

### Testing

```bash
# Run all standalone tests (no PHENIX required)
python3 tests/run_all_tests.py --quick

# Run all tests (including PHENIX-dependent)
python3 tests/run_all_tests.py

# Verbose output
python3 tests/run_all_tests.py --verbose

# Individual test file
python3 tests/tst_file_utils.py

# Tests matching a pattern
python3 tests/run_all_tests.py --pattern "directive"
```

---


---

## Directory Structure

```
improved_agent_v2/              # AI Agent system (libtbx/langchain/ + phenix/phenix/ entry points)
├── __init__.py                 # Package marker
├── VERSION                     # Package version string
├── agent/                      # Core agent logic
│   ├── advice_preprocessor.py  # README discovery, advice processing
│   ├── api_client.py           # V2 API request/response building
│   ├── best_files_tracker.py   # Track best files per type
│   ├── command_builder.py      # Unified command generation (with content guards)
│   ├── command_postprocessor.py # Server-safe command transforms (sanitize, inject)
│   ├── command_templates.json  # Program command templates and file slots
│   ├── config_loader.py        # Load decision_config.json thresholds
│   ├── contract.py             # session_info field registry + protocol version
│   ├── decision_config.json    # Tiered decision rules and thresholds
│   ├── directive_extractor.py  # Parse user directives from advice
│   ├── directive_validator.py  # Pre-validate user requests
│   ├── display_data_model.py   # Unified display data (progress/results tabs + HTML report)
│   ├── docs_tools.py           # CLI: Documentation generation
│   ├── dry_run_manager.py      # Testing: dry-run workflow simulation
│   ├── error_analyzer.py       # YAML-driven automatic error recovery
│   ├── error_classifier.py     # PERCEIVE error classes (TERMINAL/PHIL/LABEL/RETRYABLE)
│   ├── event_formatter.py      # Output formatting (verbosity levels)
│   ├── event_log.py            # Structured event logging
│   ├── failure_diagnoser.py    # LLM diagnosis of terminal errors
│   ├── file_metadata.py        # Output-file metadata (contents, quality, provenance)
│   ├── file_utils.py           # Shared file classification (MTZ type, exclude patterns)
│   ├── format_validation.py    # Format validation results as THINK-context text
│   ├── gate_evaluator.py       # Strategic planner: gate/stop evaluation
│   ├── generate_logic_doc.py   # CLI: Decision logic documentation
│   ├── generate_safety_docs.py # CLI: Safety checks documentation
│   ├── graph.py                # LangGraph state machine definition
│   ├── graph_nodes.py          # LangGraph node implementations
│   ├── graph_state.py          # Agent state type definitions
│   ├── hypothesis_evaluator.py # Strategic planner: hypothesis evaluation
│   ├── intent_classifier.py    # Classify advice intent (solve/task/tutorial) -> stop behavior
│   ├── kb_tags.py              # Derive KB lookup tags (category + tags) from state
│   ├── log_section_extractor.py # Extract informative log sections within a char budget
│   ├── memory.py               # Persistent learned syntax tips
│   ├── metric_evaluator.py     # Metric quality evaluation
│   ├── metrics_analyzer.py     # Metric trends and convergence
│   ├── mtz_inspector.py        # Client-side MTZ inspection (R-free columns)
│   ├── nl_to_phil.py           # Natural-language-to-PHIL conversion
│   ├── pattern_manager.py      # Regex pattern management
│   ├── perceive_checks.py      # PERCEIVE-node stop-condition checks (pure helpers)
│   ├── phenix_utils.py         # REST encoding, standalone PHENIX utilities
│   ├── phil_validator.py       # PHIL parameter validation
│   ├── placement_checker.py    # Unit cell comparison for model placement
│   ├── plan_generator.py       # Strategic planner: plan generation
│   ├── planner.py              # Agent planning and next-move generation
│   ├── program_registry.py     # YAML program registry
│   ├── program_validator.py    # CLI: Program config validation
│   ├── rate_limit_handler.py   # LLM rate limiting with backoff
│   ├── raw_advice_scanner.py   # Scan raw advice for crystallography filenames (LLM-free)
│   ├── rules_selector.py       # Rules-only program selection
│   ├── sanity_checker.py       # Red flag detection
│   ├── session.py              # Persistent session tracking (AgentSession, dup detection)
│   ├── session_tools.py        # CLI: Session management
│   ├── strategy_memory.py      # Persistent strategy memory across sessions
│   ├── structure_model.py      # Structure Model state (xtriage/phaser-derived)
│   ├── template_builder.py     # YAML-driven command templates
│   ├── thinking_agent.py       # THINK node: should_think()/run_think_node()
│   ├── transport.py            # Sanitization and encoding
│   ├── tst_agent_enhancements.py # In-tree tests for agent enhancements
│   ├── utils.py                # General utility functions
│   ├── validation_history.py   # Per-cycle validation history tracking
│   ├── validation_inspector.py # Headless structural validation (GUI-free) -> dict
│   ├── workflow_engine.py      # YAML workflow interpreter
│   ├── workflow_state.py       # State detection, PDB content guards
│   └── yaml_tools.py           # CLI: YAML validation and inspection
├── knowledge/                  # Configuration & domain knowledge
│   ├── api_schema.py           # V2 API schema definitions
│   ├── diagnosable_errors.yaml # Terminal errors needing LLM diagnosis
│   ├── errors.yaml             # THINK stop-reason codes (STOP_REASON_CODES source of truth)
│   ├── expert_knowledge_base_v2.yaml # Expert crystallography knowledge base
│   ├── explanation_prompts.py  # Explanation-engine prompts
│   ├── file_categories.yaml    # File categorization
│   ├── html_report_template.py # Template HTML report (DisplayDataModel; inline SVG charts)
│   ├── kb_loader.py            # Expert KB runtime loader (category/tag filtering)
│   ├── metric_patterns.py      # YAML-driven metric extraction
│   ├── metrics.yaml            # Quality thresholds
│   ├── parameter_fixes.json    # Wrong->correct parameter name mappings
│   ├── patterns.yaml           # Regex patterns
│   ├── phenix_programs.py      # Program discovery and introspection
│   ├── plan_schema.py          # Plan schema definitions/validation
│   ├── plan_template_loader.py # Load plan_templates.yaml (resolve inheritance) -> StructurePlan
│   ├── plan_templates.yaml     # Expert plan skeletons for standard workflows
│   ├── program_registration.py # Program detection from logs
│   ├── programs.yaml           # Program definitions
│   ├── prompts.py              # LLM prompts for planning and commands
│   ├── prompts_hybrid.py       # Hybrid planning prompts (rules + LLM)
│   ├── recoverable_errors.yaml # Error recovery patterns
│   ├── solve_run_info.yaml     # Minimal solve-mode hints (project_advice) per tutorial
│   ├── summary_display.py      # Quality table formatting
│   ├── thinking_prompts.py     # THINK-node prompts
│   ├── transport.yaml          # Sanitization rules
│   ├── tutorial_expectations.yaml # Expected programs/method/results per tutorial (dual-run)
│   ├── workflows.yaml          # Workflow state machines
│   └── yaml_loader.py          # Configuration loading
├── phenix_ai/                  # Runtime entry points (in phenix/phenix/)
│   ├── local_agent.py          # Local execution (same process)
│   ├── remote_agent.py         # Server execution (REST API)
│   ├── run_ai_agent.py         # Decision engine (graph execution)
│   ├── run_ai_analysis.py      # Log analysis (standalone)
│   ├── log_parsers.py          # Log metric extraction
│   ├── utilities.py            # Shared runtime utilities
│   ├── agent_interface.py      # BaseAgent ABC (decide_next_step interface)
│   ├── install_ai_tools.csh    # C-shell installer for the AI tools
│   ├── cohere_api_key.html     # Cohere API-key setup help page
│   ├── google_api_key.html     # Google API-key setup help page
│   ├── openai_api_key.html     # OpenAI API-key setup help page
│   └── README                  # Package README
├── programs/                   # PHENIX program integration (AI entry points, in phenix/phenix/)
│   ├── ai_agent.py             # Main PHENIX entry point
│   └── ai_analysis.py          # Log analysis entry point
├── analysis/                   # Log analysis and post-run analysis
│   ├── analyzer.py             # RAG-based log analysis
│   ├── log_info.py             # High-level log info extraction
│   ├── state_extractor.py      # Project state extraction from logs
│   ├── summarizer.py           # Map-reduce log summarization
│   └── agent_session_analyzer.py # Session performance analysis
├── core/                       # LLM integration
│   ├── llm.py                  # Provider abstraction (Google, OpenAI, etc.)
│   ├── types.py                # Core data types (AgentPlan, etc.)
│   ├── _build_info.py          # Auto-generated build info
│   └── _version.py             # Auto-generated version string
├── commands/                   # Command building framework
│   └── base.py                 # Abstract base class for command builders
├── strategies/                 # Planning strategy framework
│   └── base.py                 # Abstract base class for planning strategies
├── examples/                   # Example scripts (package marker only)
├── rag/                        # RAG (Retrieval-Augmented Generation)
│   ├── document_loader.py      # Document loading and chunking
│   ├── retriever.py            # Vector DB retrieval with FlashRank reranking
│   ├── vector_store.py         # Chroma vector store creation
│   └── _chroma_resilience.py   # Lazy chromadb import/probe (resilient optional dependency)
├── utils/                      # General utilities
│   ├── query.py                # Documentation query interface
│   ├── run_utils.py            # Log parsing and HTML output
│   └── text_processing.py      # Text block extraction helpers
├── phenix_knowledge.py         # Phenix program allow-list and syntax hints
├── phenix_learned_info/        # Persistent learned knowledge
│   └── phenix_learned_memory.json # Learned syntax tips per program
├── run_build_db.py             # CLI: Build RAG documentation database
├── run_inspect_db.py           # CLI: Inspect RAG database contents
├── run_query_docs.py           # CLI: Query PHENIX documentation
├── findings/                   # Generated static-audit reports
│   ├── phase_0_static_audit.yaml                 # Static audit
│   ├── phase_1_contract_gaps.yaml                # Contract gaps
│   ├── phase_2_path_divergences.yaml             # Path divergences
│   ├── phase_3_roundtrip_failures.yaml           # Round-trip failures
│   ├── phase_4_flag_mismatches.yaml              # Flag mismatches
│   ├── phase_5_classification_disagreements.yaml # Classifier disagreements
│   ├── phase_6_missing_inputs.yaml               # Missing inputs
│   ├── phase_8_command_failures.yaml             # Command failures
│   └── phase_9_llm_perturbation.yaml             # LLM perturbation
├── tests/                      # Test suites (abbreviated; many tst_*.py)
│   ├── run_all_tests.py        # Test runner (registers all suites)
│   ├── tst_langchain_tools.py  # Legacy-module unit tests
│   ├── tst_hardcoded_cleanup.py # Conformance guards
│   ├── tst_utils.py            # Assert helpers (cctbx-style)
│   └── scenarios/              # Dry-run test scenarios
└── docs/                       # Documentation
    ├── OVERVIEW.md             # Technical overview (this file)
    ├── ARCHITECTURE.md         # Component deep-dive
    ├── DEVELOPER_GUIDE.md      # How-to guides (adding programs, testing, safety)
    ├── USER_GUIDE.md           # User guide + directives reference
    ├── CHANGELOG.md            # Version history
    ├── CHANGELOG_ARCHIVE.md    # Archived older version history
    ├── CCTBX_LLM_PROGRAMMING_GUIDELINES.md     # cctbx LLM-assisted dev guide
    ├── AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md  # Agent-specific supplement
    └── prompts/                # LLM workflow prompts & handoff templates
        ├── HOW_TO_PROGRAM_WITH_AN_LLM.md       # Methodology guide
        ├── WORKFLOW.md                         # Plan/Review/Execute workflow
        ├── WORKFLOW_PROMPT.txt                 # Workflow-stage prompt
        ├── PLAN_PROMPT.txt                     # Plan-stage prompt
        ├── REVIEW.txt                          # Review-stage prompt
        ├── CONTINUE_PROMPT.txt                 # Continuation prompt
        ├── HANDOFF.json                        # Session handoff checkpoint
        ├── CCTBX_LLM_PROGRAMMING_GUIDELINES.md     # (prompt copy)
        └── AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md  # (prompt copy)
```

---

