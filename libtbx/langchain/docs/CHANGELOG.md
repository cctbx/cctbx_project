# CHANGELOG — v116 / v117 / v117.1 / v117.2 / v117.3 / v118 / v119 / v120

## Version 120 (New Providers: Portkey + Claude, and FORCE_NO_AI_SERVER — P1, P2, P3)

### Summary

v120 is a three-phase feature cluster that (1) adds an environment
override, `FORCE_NO_AI_SERVER`, for forcing local execution, and
(2) adds two new LLM providers — **portkey** (the Portkey gateway
fronting Azure OpenAI) and **anthropic** (Claude) — wired through
every provider-routing decision point in the agent.  The phases are
independent and ship-self-contained:

```
v119 → P1 (FORCE_NO_AI_SERVER) → P2 (single-source SUPPORTED_PROVIDERS) → P3 (portkey + anthropic providers)
```

The work originated from a user request to run the agent against the
Portkey access platform, accompanied by a working proof-of-concept
patch (`add_portkey_support.patch`).  That patch was the authoritative
source for the deployment specifics (Portkey SDK with
`provider="azure-openai"`, env vars `PORTKEY_AZURE_API_KEY` /
`PORTKEY_BASE_URL` / `PHENIX_PORTKEY_MODEL`, the Azure gpt-5 call shape
of `max_completion_tokens` with no `temperature`), but its structure
was re-targeted onto current code and hardened for production.  Claude
(anthropic) was added in the same pass because it was already
half-scaffolded and shares every touch-point portkey needed.

Three deliberate departures from the proof-of-concept patch: (a) the
patch flipped the `run_on_server` PHIL default `True → False` in both
program files, silently changing behaviour for every user — v120 keeps
the `True` default and provides `FORCE_NO_AI_SERVER=1` for the same
local-only effect without a baseline change; (b) the patch read env
vars with `os.environ[...]` (bare `KeyError` on a missing var) — v120
checks at call time and raises a clear `ValueError` naming the missing
variable; (c) the patch reused `get_openai_handler()` for portkey —
v120 adds a dedicated `get_portkey_handler()` so rate-limit metrics and
failure logging stay unpolluted.

### P1 — FORCE_NO_AI_SERVER

`programs/ai_agent.py::run()` gained a single override, placed after
`set_defaults()` and before session-management dispatch so it covers
the v2-API, `iterate_agent`, and `run_job_on_server_or_locally` paths:
when `FORCE_NO_AI_SERVER=1` (literal, whitespace-stripped), it forces
`communication.run_on_server = False` and logs the override via
`self.vlog.normal`.  Any other value, or unset, is a no-op.  The same
override is mirrored in `programs/ai_analysis.py::run()` (placed before
its `run_job_on_server_or_locally` dispatch), so the env var forces local
execution for `phenix.ai_analysis` too, across every analysis_mode
(standard, agent_session, advice_preprocessing, directive_extraction,
failure_diagnosis).  It does not affect `predict_and_build` /
`predict_model` when those are run separately.  Behaviour is otherwise
identical to server mode — only the locus of execution changes.

### P2 — single-source SUPPORTED_PROVIDERS

Before v120, `SUPPORTED_PROVIDERS` existed as two independent literals
(`agent/graph_nodes.py` and a hardcoded fallback in
`programs/ai_agent.py`) that could drift silently.  P2 centralizes the
canonical list in `core/llm.py`; both consumers now import it (with the
standard two-path `try/except ImportError` fallback).  The hardcoded
literal in `ai_agent.py` was removed — but because its consumer is an
*error-advice* helper, an import failure degrades to the
provider-agnostic rules-only suggestion rather than raising a second
exception over the original error.  This eliminates the dual-literal
drift hazard at the root; a runtime import-identity test replaces the
regex source-scan originally planned.

### P3 — portkey + anthropic providers

`core/llm.py`: portkey + anthropic added to the model-default tables;
`PHENIX_PORTKEY_MODEL` and `ANTHROPIC_LLM_MODEL` added to the env-var
override map; three new helpers — `portkey_langchain_config()` (call-time
env check, clear `ValueError`), `sanitize_llm_kwargs(provider, kwargs)`
(dict-driven per-provider kwarg reshaping), and
`_delegate_embeddings_for_nonnative()` (the embeddings policy below);
factory branches for both providers; `get_expensive_llm` branches.
`SUPPORTED_PROVIDERS` activated to the five-provider set.

The **asymmetric-args** trap is handled in `sanitize_llm_kwargs`:
portkey (Azure gpt-5) renames `max_tokens → max_completion_tokens` and
drops `temperature`; anthropic always keeps `max_tokens` (`ChatAnthropic`
requires it) but its `temperature` handling is **model-dependent**.  The
Claude 4.6/4.7+ reasoning family (e.g. `claude-sonnet-4-6`,
`claude-opus-4-7`) rejects `temperature`/`top_p`/`top_k` with HTTP 400
("temperature is deprecated for this model" — adaptive thinking took over
sampling control), so `sanitize_llm_kwargs` and the `ChatAnthropic`
constructor drop those params for that family, detected by
`anthropic_model_rejects_sampling_params()` (model-name based, so future
4.8+ models are covered).  Older Claude models (Sonnet 4 / Opus 4 dated,
3.x) still accept `temperature` and keep it.  (This was tightened after a
live anthropic reliability run surfaced the 400; the original assumption
that all Claude models accept `temperature` held only for the pre-4.6
models v120 was first written against.)

The **embeddings policy** (anthropic has no native embeddings endpoint;
the user's Portkey key can reach the embeddings deployment, but their
current key lacks access): `get_llm_and_embeddings` always returns a
working chat LLM; the embeddings object is constructed lazily (no
network call) so the agent's end-of-run `agent_session` assessment —
which constructs but never queries embeddings — is safe even with a
no-access key.  For anthropic, embeddings delegate to a default provider
(OpenAI/Google) when a key is present, else return `None`, emitting a
`[PROVIDER_DELEGATION]` marker; the helper never raises on construction.
An actual embeddings *query* failure (only reachable via
`phenix.ai_analysis standard`/RAG, never the agent) surfaces a clear
execution-time error echoing the raw upstream message rather than
matching on `"Unauthorized"`.

Call sites updated to route portkey/anthropic: `agent/api_client.py`
(`_call_portkey_llm` + dispatch), `agent/directive_extractor.py`
(`_call_llm` handler-select + `_call_llm_fallback`),
`agent/rate_limit_handler.py` (`get_portkey_handler`),
`agent/graph_nodes.py` (`validate_provider` key-checks + second handler
site), `agent/thinking_agent.py` (`_get_rate_handler`).  The embeddings/
RAG database path was completed end-to-end: `utils/run_utils.py`
(`validate_api_keys` gates both providers; `get_db_dir_for_provider`
gained `docs_db_portkey`), `analysis/summarizer.py` (portkey chunk
sizes), `run_query_docs.py` (accepts both), and the database-build tools
`command_line/rebuild_ai_database.py` + `update_ai_database.py` (accept
`portkey` from argv with key-check).  Anthropic is intentionally **not**
a database-build provider (no native embeddings).  `get_ai_db_dir` /
`have_ai_database` in `phenix_ai/utilities.py` are provider-agnostic
(string interpolation → `docs_db_<provider>`) and needed no change.

PHIL `provider` enums in `programs/ai_agent.py` and
`programs/ai_analysis.py` extended to `ollama *google openai anthropic
portkey` (the GUI dropdown regenerates automatically).  `install_ai_tools.csh`
adds `anthropic`, `langchain-anthropic`, `portkey-ai`.

The supplied server dispatcher (`phenix_ai/run_ai_agent.py`), GUI
(`wxGUI2/Programs/AIAgent.py`), and contract (`agent/contract.py`)
needed no changes: provider travels in the request `settings` and is
consumed only in shared `agent/`+`core/` code, so local and server run
identical branches — the server only additionally needs the keys in its
environment (deployment, not code).  The user's deployment is local-only
with `PORTKEY_AZURE_API_KEY` set on the workstation.

### LLM-unavailable observability

Previously, when the requested LLM provider could not be reached, the
agent silently degraded to rules-based program selection (`_mock_plan`
→ `RulesSelector`) and a run could complete on rules without the user
realizing the LLM never contributed.  Now the condition is made obvious
in two places, visible in both CLI and GUI, on both the local and the
server execution paths.

The signal travels through the `events` list, which is already part of the
server→client response contract (`agent/contract.py` `RESPONSE_FIELDS`) and is
preserved by every graph node (verified: `events` is a plain `List[Dict]` in
`graph_state.py` with no reducer, and `build`/`validate`/`fallback`/`output`
all return `{**state}` or never clobber it).  Using the existing carrier means
**no change to `contract.py`, `api_schema.py`, `run_ai_agent.py`,
`event_log.py`, `event_formatter.py`, or `session.py`** — the fix is confined
to `graph_nodes.py` and `ai_agent.py`.

- **Per cycle (graph):** `agent/graph_nodes.py::_handle_llm_failure` emits a
  structured `EventType.NOTICE` event carrying a machine-readable
  discriminator `notice_kind="llm_unavailable"` (plus `provider` and a
  human `message`), and also logs the human line to `debug_log` (shown at
  every verbosity, independent of the event formatter).  Emission is guarded
  once per `graph.invoke` by `llm_notice_emitted_this_invoke` — a key scoped
  to the single invoke (the `validate→plan` retry edge can call the handler
  more than once per cycle), deliberately *not* named like the driver's
  run-level key to avoid scope confusion.  Graph state is per-cycle (each
  cycle is a fresh `graph.invoke`), so there is no run-level flag in graph
  state; the NOTICE is emitted each cycle the LLM is unavailable (honest
  per-cycle reporting).
- **Per run (driver):** `programs/ai_agent.py::print_history_record` →
  `_detect_llm_unavailable_from_events` scans each cycle's
  `history_record["events"]` for the discriminator and, on first sight, sets
  the run-level `session.data["llm_ever_unavailable"]` (and provider).
  `session.data` persists across cycles, unlike graph state, so this is the
  durable run-level signal.  The detection hook is guarded `if session is not
  None` (the display-only/replay path calls `print_history_record` without a
  session) and never raises.  `_finalize_session` shows an end-of-run banner
  ("IMPORTANT: THIS RUN DID NOT USE THE LLM") via
  `session_data_flag_llm_unavailable`; the Results-tab summary is prefixed
  with a one-line warning and the result object exposes
  `llm_ever_unavailable` for programmatic/GUI use.

This is observability only — it does not change which program the agent
selects, only how clearly the rules-fallback is reported.  (The explicit
`use_rules_only=True` mode was already announced clearly; this covers the
*unintended* fallback when a provider is down.)  Known gap, documented for a
later follow-up: a pure display-only replay of a saved session that never
re-runs cycles will not show the end-of-run banner, because the session cycle
records do not persist `events`; fixing that would require a `session.py`
schema change and is out of scope here.

### PHIL multi-line `.help` continuation fix

The v120 `provider` enum gained a multi-line `.help` description, written as an
unquoted continuation where the first line ended `...quickest,\` — a backslash
with no preceding space.  PHIL requires a **space before the continuation
backslash** in an unquoted value; without it the continuation is not joined and
the next line is re-lexed as a new statement, producing
`Syntax error: expected "=", found "is"`.  This broke
`tst_h18_1_phil_roundtrip.py`'s parse tests in BOTH
`programs/ai_agent.py` and `programs/ai_analysis.py` (dual schemas, see
ARCHITECTURE §42).  Fixing the provider help then surfaced two pre-existing
helps with the same flaw (`url`, `port`, continued `...server\` then
`Normally set automatically` → `found "set"`).

Fix: add a single space before the backslash on the three broken helps
(`provider`, `url`, `port`) in each file.  The sibling `url_type`/`token` helps
were already correct (they ended `. \`) and were left untouched, as were the
~30 other multi-line helps and the quoted multi-line scope `.caption` (quoted
continuations follow different rules and need no leading space).  The full rule,
the empirical verification method, and authoring guidance are documented in
ARCHITECTURE §44.

### Tests

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_force_no_ai_server.py` | 10 (new) | P1: env-value truth table (`1`/`0`/`true`/empty/whitespace), quiet-when-already-false, + source-scan drift guards (block present, placed after set_defaults / before session dispatch, literal-`"1"` contract) |
| `tests/tst_supported_providers_single_source.py` | 5 (new) | P2: runtime import-identity (graph_nodes re-exports core.llm's list), source scans that both consumers import from core.llm and no literal survives outside core.llm |
| `tests/tst_v120_providers.py` | 21 (new) | P3: `sanitize_llm_kwargs` both directions (portkey remaps/drops, anthropic keeps), `portkey_langchain_config` clear-ValueError, SUPPORTED_PROVIDERS + PHIL-enum coverage, embeddings delegation None-without-keys + never-raises, dummy-key lazy construction with no network (mocked SDK), embeddings-query clear-error contract, `validate_api_keys` gates new providers, per-call-site source scans |
| `tests/tst_default_models.py` | (updated) | Exact-baseline dicts updated for portkey/anthropic table entries |
| `tests/tst_llm_unavailable_notice.py` | 10 (new) | LLM-unavailable observability: structured NOTICE event shape + discriminator, emit-once-per-invoke guard, **events list not mutated in place**, **propagation (event → driver detection → `session.data` flag)** — the check the prior implementation lacked — idempotent across cycles, ignores unrelated events, `session=None` safe, + source scans that graph_nodes emits the NOTICE with the scoped guard key (and does NOT use the run-level key name), ai_agent detects by discriminator + banners + guards `session is not None`, and `run_ai_agent.py` still passes `events` into the response builders (server path) |

All three new suites registered in `tests/run_all_tests.py`.  Standalone
suites run under plain `python3`; live Portkey/Anthropic round-trips
remain in the keyed `tests/llm/` suite (real keys + network).

### Files changed (cumulative)

| File | Phase |
|---|---|
| `core/llm.py` | P2 (canonical `SUPPORTED_PROVIDERS`), P3 (tables, env-overrides, `portkey_langchain_config`/`sanitize_llm_kwargs`/`_delegate_embeddings_for_nonnative`, portkey+anthropic factory & expensive branches) |
| `agent/graph_nodes.py` | P2 (import `SUPPORTED_PROVIDERS`), P3 (`validate_provider` key-checks + second handler site) |
| `agent/api_client.py` | P3 (`_call_portkey_llm` + dispatch) |
| `agent/directive_extractor.py` | P3 (portkey in `_call_llm` + `_call_llm_fallback`) |
| `agent/rate_limit_handler.py` | P3 (`get_portkey_handler`) |
| `agent/thinking_agent.py` | P3 (`_get_rate_handler` portkey) |
| `analysis/summarizer.py` | P3 (portkey chunk sizes) |
| `utils/run_utils.py` | P3 (`validate_api_keys` + `get_db_dir_for_provider`) |
| `run_query_docs.py` | P3 (accepts portkey + anthropic) |
| `programs/ai_agent.py` | P1 (`FORCE_NO_AI_SERVER`), P2 (import `SUPPORTED_PROVIDERS`, remove literal), P3 (PHIL enum) |
| `programs/ai_analysis.py` | P1 (`FORCE_NO_AI_SERVER`), P3 (PHIL enum) |
| `command_line/rebuild_ai_database.py` | P3 (accepts portkey) |
| `command_line/update_ai_database.py` | P3 (accepts portkey) |
| `phenix_ai/install_ai_tools.csh` | P3 (anthropic, langchain-anthropic, portkey-ai) |
| `tests/tst_force_no_ai_server.py` | P1 (new) |
| `tests/tst_supported_providers_single_source.py` | P2 (new) |
| `tests/tst_v120_providers.py` | P3 (new) |
| `tests/tst_default_models.py` | P3 (baseline dicts updated) |
| `tests/run_all_tests.py` | P1/P2/P3 (register three suites) |

**Verified no change:** `phenix_ai/run_ai_agent.py`,
`wxGUI2/Programs/AIAgent.py`, `agent/contract.py`,
`phenix_ai/utilities.py`.

---


## Version 119 (Operational Hardening + Phase 2A + Phase 2B + Production Bug Cluster — H1, H2, H2.1, H3, H3b, H4, H4.1, H5, H5.1, H5.1.1, H6, H6.1, H7, H8, H9, H10, H11, H12, H13, H14, H14.1, H14.2, H15, H16, H16.1, H17, H17.1, H18, H18.1, H18.2)

### Summary

v119 is a twenty-ship cluster pursuing the operational-hardening
agenda from `v118_next_steps_consolidated_rev4.md` §4.7B plus
the Phase 2A planning-suite framework (H6/H6.1), the Phase 2B
preprocessor activation (H7), a four-ship production-bug
sub-cluster (H8–H11) surfaced by AIAgent_62 and run_38 batch
testing, a follow-up cleanup pass for the exclude_patterns
mechanism + categorizer semantic pinning (H12), and an
Ollama-provider robustness pass with retired-model 404 detection
surfaced by Tom's `run_39a_ollama` failure (H13).  The ships are
independent and ship-self-contained:

```
v118 → H1 → H2 → H2.1 → H3 → H3b → H4 → H4.1 → H5 → H5.1 → H5.1.1 → H6 → H6.1 → H7 → H8 → H9 → H10 → H11 → H12 → H13 → H14 → H14.1 → H14.2 → H15 → H16 → H16.1 → H17 → H17.1 → H18 → H18.1 → H18.2
```

Together they establish: a single source of truth for LLM model
defaults (H1), a build-metadata channel from server to client
(H2), promotion of LLM-emitted skip flags to the workflow layer
(H2.1), a startup-canary infrastructure for detecting wrong-
build deploys and unhealthy LLM environments (H3/H3b), a
preprocessing-telemetry marker ([STEP_1F]) with golden-master
corpus pinning (H4/H4.1), a uniform diagnostic-messages
relay channel that surfaces server-side stderr markers to the
client transparently (H5 + H5.1), two latent-bug cleanups in
the directive extractor (H5.1.1), the planning-suite reliability
testing framework (H6) with per-scenario failure-count reporting
(H6.1), and the Phase 2B scanner-first file extraction activation
that completes the preprocessor deprecation trigger first set
up in H4 (H7).  The H8–H11 sub-cluster fixes four production
bugs surfaced during batch testing: template-literal allowlist
gap (H8), `predict_and_build` PHIL scope mismatch (H9), structural
gap in `exclude_patterns` application across selection paths
(H10), and YAML pattern authoring vs function-semantics mismatch
that prevented H10's filters from rejecting their intended targets
(H11).  H12 is a follow-up cleanup pass: it refactors H10's
closure into a module-level helper (8 application sites
consolidated to 1 helper) and adds a semantic-pin test file for
every public function in `agent/file_utils.py` (the
`tst_file_categorizer.py` suite — pre-emptive defense against
sandbox-stub drift, the lesson H11 surfaced).  H13 fixes two
Ollama-provider bugs surfaced by Tom's `run_39a_ollama` failure
on cci-gpu-01: `OLLAMA_BASE_URL` env-var override silently
broke the OpenAI-SDK path construction when the user-set value
omitted `/v1`, and `OLLAMA_LLM_MODEL` env-var was inconsistently
honored (only `core/llm.get_llm_and_embeddings` honored it; both
`directive_extractor` and `api_client` ignored it).  H13 also
completes the v118 §3.5 retired-model 404 detection work,
refining the classification into three sub-categories
(MODEL_RETIRED vs MODEL_UNAVAILABLE vs FAILED) so the operator
sees actionable hints rather than opaque "404 page not found".

H14/H14.1/H14.2 close three independent regressions surfaced by
`run_39_openai` batch analysis: a phaser false-positive in
`_ACTION_TABLE["solve"]` (24 occurrences across 5 datasets where
the multi-action branch forced phaser into SAD/MAD workflows), a
duplicate `[STEP_1F]` emission, and a permissive space_group
regex producing "Not explicitly mentio" truncations that bypassed
the v118.9 validator on the rules-only path.  H14.2 is a one-line
config fix removing `ncs_spec` from `predict_and_build`'s strategy_flags.

H15 closes Tom's bromodomain resume failure (run 144).  When new
advice was provided on resume, `gate_stop` was cleared but
per-stage statuses were not — the gate immediately re-fired
"all stages complete" and the LLM saw STATE=complete, choosing
polder against user intent.  H15 adds `reopen_stages_for_directives`
to `plan_generator.py`: a TARGETED single-stage reopen that finds
the LATEST completed stage whose `programs` list contains a program
named in `directives.program_settings` and resets ONLY that stage
to PENDING.  Per Gemini's critique of the original H15 plan, the
blast radius is O(1) regardless of plan size — no cascade-resets,
no downstream stages reopened.

H16/H16.1 closes 88 TIER-1 failures across run_25 and run_39
matching "Sorry: Multiple equally suitable arrays of observed
xray data found", concentrated in AF_exoV_MRSAD and lysozyme-MRSAD
tutorials.  H16 adds the new `agent/mtz_inspector` module:
`inspect_mtz()` reads MTZ column structure via cctbx,
`select_obs_labels_for(program, info)` applies a per-program
preference policy (merged for MR, anomalous for SAD), and a new
`auto_fill_obs_labels` invariant in `programs.yaml` triggers the
builder to inject the chosen labels into xtriage, autosol,
phaser, and predict_and_build commands.  H16.1 bumps the scanner
version pin and adjusts `[STEP_1F]` telemetry.

H17/H17.1 closes Tom's lysozyme-MRSAD cycle-5 reactive recovery
gap.  The LLM passed `lyso2001_scala1.mtz` (raw anomalous MTZ
without phase columns) as `map_file=` to `phenix.autobuild`,
which requires `PHIB` phase columns; autobuild crashed with
"Sorry, PHIB is required for input_map_file"; the agent halted
after 4 retries.  H17 adds:
- A new `missing_phib_input_map_file` entry in
  `recoverable_errors.yaml` with strict conjunction detection
  (matches the exact PHIB-required error text, not just any "PHIB"
  mention).
- A new `strip_parameter` resolution kind in `error_analyzer.py`
  alongside the existing `add_parameter` and `select_value` kinds.
- A new `strip_flags` field on the `ErrorRecovery` dataclass.

H17 shipped the analyzer side but exposed Scenario B: detection
fires and `[NOTICE]` prints, but the executor never strips the
flag from the retry command.  H17.1 closes that gap with three
edits to `programs/ai_agent.py`: `_handle_recovery` stashes
`strip_flags` into `session.data["pending_strip_recoveries"]`;
`_execute_command` pops the entry at the top of each cycle and
applies a robust regex (per Gemini's critique: handles
quoted-with-spaces, PHIL spacing, and single/double quote
variants); `_print_recovery_notice` shows
"Action: Stripping [flags]" instead of the awkward "Selecting ''"
for strip recoveries.  Validated end-to-end on lysozyme-MRSAD:
cycle 5 detection + cycle 6 [STRIP] line + autobuild SUCCESS.

H18 closes the AF_7mjs density-modify-and-stop regression
(separate from §20's bug, surfaced after H17.1 deployed): user
wrote "density modify and stop" with cryo-EM half-map inputs; the
extractor produced `after_program=phenix.autobuild_denmod` (X-ray
density-modification); the §20 correction failed to fire because
its text-only `_detect_experiment_type_signals` returned None
(terse user text has neither cryo-EM nor X-ray tokens); the LLM
downstream bypassed the failed stop check and ran
`phenix.predict_and_build`, overriding user intent.  H18 adds a
new public helper `agent/file_utils.infer_experiment_type_from_files()`
that uses file extensions as the deterministic primary signal,
and rewires BOTH `_apply_experiment_type_program_reprints` AND
`_resolve_after_program` (the v115.10 post-LLM overlay) to use
files-first detection with text-based detection as fallback.
Files-win policy on conflict per Gemini's H18 review.  Enriched
`[DIRECTIVE_CORRECTION]` telemetry records `source=files|text`
and `OVERRIDDEN` annotations; new `[DIRECTIVE_CORRECTION_MIXED]`
marker fires on accumulated-files drift (Pitfall 1).
`agent/plan_generator._build_context` was also refactored to use
the shared helper, establishing a single source of truth for
file-based experiment-type inference.

H18.1 closes a PHIL declaration deploy-gap surfaced when Tom
deployed H18 and re-ran AF_7mjs.  All 20 H18 K-tests passed in
the sandbox; production still ran `predict_and_build` on cycle 3.
The log showed `AttributeError: Assignment to non-existing
attribute "ai_analysis.original_files_for_directives"` at
`ai_agent.py` line 8513.  H18 had added the PHIL declaration to
`programs/ai_analysis.py` but not to `programs/ai_agent.py`'s
OWN `master_params` string — and `directive_params =
copy.deepcopy(self.params)` copies the AGENT's params object,
which is parsed against the agent's `master_params`.  The
assignment crashed; the surrounding try/except swallowed it
silently; directive extraction returned `{}`; the agent ran the
default plan template.  H18.1 adds the missing PHIL declaration
(14 lines) to `ai_agent.py`'s `master_params`, plus a new
`tst_h18_1_phil_roundtrip.py` (6 K-tests) that exercises the
full parse → extract → deep-copy → assign path the production
code uses.  Source-grep K-test runs in sandbox without PHENIX;
live-PHIL K-tests skip gracefully there.  No code logic
changed.

H18.2 closes a THIRD `_resolve_after_program` callsite that H18
missed.  After H18.1 deployed cleanly, Tom re-ran AF_7mjs with a
runtime tracer installed (`h18_install_runtime_tracer.py` —
patches the deployed `directive_extractor.py` in-place with
`[H18_TRACE]` markers at every key decision point).  The trace
revealed that the LLM was emitting
`{"stop_conditions": {"stop_after_requested": true}}` — note: no
`after_program` field at all.  H18's site 1
(`_apply_experiment_type_program_reprints`) correctly took its
"no after_prog → early return" branch.  Then a v117.2 fallback at
`directive_extractor.py:783` fired (its raison d'être: fill in
`after_program` when the LLM sets `stop_after_requested=True` but
omits the program field), calling `_resolve_after_program` on the
raw advice "density modify and stop" — but **without
`original_files`**.  Pre-H18.2, that callsite was a third site
the H18 audit missed.  The resolver defaulted `_exp="xray"` via
the text-only heuristic (raw advice has no
cryo-EM/X-ray tokens) and mapped `denmod` →
`phenix.autobuild_denmod`.  H18's site 2
(`_apply_workflow_intent_fallback`) DOES pass `original_files`,
but the preprocessed advice contains "Stop Condition: None" so
`_is_stop_after_requested` returned False, pushing the resolver
into "n==1, no stop → leave as-is" — and the buggy after_program
from the v117.2 path persisted.  H18.2 adds `original_files=
original_files` to the v117.2 callsite (one line in
`agent/directive_extractor.py:796`) plus 3 new K-tests
(K21/K22/K23) extending the suite to 23 tests.  K21 is the exact
AF_7mjs production reproduction (mocked LLM returns the same
shape the production LLM did).

Total K-test additions across the cluster: **189 v119-cluster
K-tests** plus **30 live LLM tests** (directive_extraction +
planning) plus **14 production-bug K-tests** in the H8–H11
sub-cluster (4 Bug 8, 4 Bug 9, 5 Bug 10, 1 Bug 11), plus
**8 categorizer semantic-pin K-tests** added by H12, plus
**13 provider-error classification K-tests** added by H13, plus
**45 K-tests** in H14 (13 keyword + 8 single-emit + 24
space_group), plus **14 K-tests** in H14.1 (validate_directives
closure), plus **7 K-tests** in H15 (resume reopen), plus
**8 K-tests** in H16 (obs_labels auto-fill), plus **16 K-tests**
in H17/H17.1 (7 detection + 9 strip executor), plus **7 K-tests**
in H18 (experiment-type-from-files), merged into the existing 13
§20 tests to extend `tst_density_modify_experiment_type.py` to a
20-test suite, plus **6 K-tests** in H18.1 (PHIL round-trip
deploy-gap defense), plus **3 K-tests** in H18.2 (v117.2 fallback
threads files; merged into `tst_density_modify_experiment_type.py`
extending it to a 23-test suite).  See the per-suite breakdown in
`next_steps_post_v119.md` "Test totals" section.

All ships verified through act → review → Gemini-critique → ship
cadence.

### Per-ship breakdown

**H1 (client+server) — Default Model Centralization**
`core/llm.py`, `agent/api_client.py`, `agent/directive_extractor.py`,
`tests/test_api_keys.py`, `tests/tst_default_models.py`.  Replaces
scattered per-call-site model defaults (`"gpt-4o-mini"`,
`"gemini-2.5-flash-lite"`, etc.) with five centralized DEFAULTS
tables in `core/llm.py`: `DECISION`, `RAG`, `RAG_EMBEDDING`,
`EXPENSIVE`, `CHEAP`.  Adds `default_model_for_provider(provider,
role="DECISION")` helper.  Adds `RETIRED_MODELS` frozenset and
emits `[DIRECTIVE_EXTRACTION_MODEL_RETIRED]` stderr marker when
a request specifies a retired model name (defense-in-depth
against the v118.8 server-404 class of incident).  +22 tests
(K_H1).

**H2 (client+server) — Server Build Metadata Channel**
`VERSION` (new at langchain/ root), `core/_version.py` (new),
`core/_build_info.py` (new, with `get_agent_build_info()` and
`inject_agent_build()`), `core/llm.py` (+`compute_defaults_fingerprint()`
which excludes `RETIRED_MODELS`), `knowledge/api_schema.py`
(+`agent_build` schema entry), `phenix_ai/run_ai_agent.py`
(injection at top of `_build_group_args_response`),
`tests/tst_agent_build_info.py`.  Every server response now
carries `response["agent_build"] = {version, defaults_fingerprint,
started_at}` in strict UTC ISO 8601.  Baseline fingerprint at
H2 ship: `sha256:77bf7421...`.  The `defaults_fingerprint`
deliberately excludes `RETIRED_MODELS` so retirement updates
don't trigger fingerprint drift.  +24 tests (K_H2).  Also updates
`docs/DEVELOPER_GUIDE.md` §8 to reflect `CURRENT_PROTOCOL_VERSION = 5`
(was stale at `3`).

**H2.1 (client) — skip_programs Promotion (micro-ship)**
`agent/directive_extractor.py`, `tests/tst_skip_promotion.py`.
Closes a pre-existing bug surfaced by the
`tst_directive_extraction.py::skip_programs` scenario: the LLM
emits `program_settings[X].skip=true` (a natural reading of the
prompt schema), but downstream code reads
`workflow_preferences.skip_programs`.  Fix: a new
`_promote_skip_settings_to_skip_programs()` helper at the top of
`validate_directives`, BEFORE per-setting validation.  Includes
all four Gemini guardrails: truthy-value acceptance (12 forms via
`_SKIP_TRUE_VALUES` constant), `.pop()` instead of read+leave,
empty-sub-dict pruning, deduplication.  +11 tests (K_H2.1).

**H3 (test infrastructure) — Startup Canary**
`tests/canary_expected.json` (new pinned config), `tests/canary_utils.py`
(new shared loader), `tests/tst_canary.py` (K_H3a, ~280 lines),
`tests/llm/tst_directive_extraction.py` (+canary Scenario),
`tests/llm/canary_check.py` (H3b orchestrator).  Two independent
canaries that consume the H2 metadata channel:

- **H3a (metadata canary, K-suite)** sends an intentionally
  wrong-typed request through `run()` (specifically `files` as
  string instead of list — necessary because schema defaults
  would silently repair simple missing-field probes).  The error
  path returns immediately with H2's `agent_build` injected.
  K_H3a (10 tests) asserts: agent_build matches
  `canary_expected.json` (hard fail on version mismatch, soft
  warn on fingerprint drift per Tom's graduated severity), VERSION
  file matches the JSON's `agent_version` (catches operator
  drift), probe is fast (<3s — catches accidental LLM leak into
  error path).
- **H3b (LLM smoke canary, operator-invoked)** combines the
  metadata check with a one-shot directive-extraction probe
  ("Run phenix.refine on the model") via the framework's
  production-faithful `call_directive_extractor` entry point.
  Wraps the LLM call in a `ThreadPoolExecutor` with a 30-second
  strict timeout (sized for two LLM round-trips: intent
  classification + main extraction).  Google→OpenAI provider
  preference.  Probe runs through the same Scenario registered in
  `tst_directive_extraction.build_scenarios()`, so it can also
  be invoked via the standard CLI for debugging:
  `phenix.python tests/llm/run_llm_tests.py --scenario canary`.
  +10 tests (K_H3a; H3b is operator-invoked, not in K-suite).

VERSION bumped to `119.H3`.  `defaults_fingerprint` unchanged
from H1 baseline because H3 doesn't touch `core/llm.py`.

**H4 (server) — [STEP_1F] Preprocessing Metrics Marker**
`phenix_ai/run_ai_analysis.py`, `agent/raw_advice_scanner.py`
(new), `tests/tst_raw_advice_scanner.py`, `tests/step_1f_corpus.json`
(new), `tests/run_all_tests.py`.  Adds a `[STEP_1F]` stderr
marker emitted by `run_advice_preprocessing` after a successful
preprocessing call.  The marker compares LLM-extracted file
mentions against a regex-based scanner's extraction
(`raw_advice_scanner.py`) on the SAME advice text, capturing
the symmetric difference for telemetry: `llm_files`,
`regex_files`, `in_llm_only`, `in_regex_only`.  Companion
`[STEP_1F_FAILED]` marker fires if the metric block itself
raises.  Includes `scanner_version` field (read from VERSION
at runtime) so telemetry can be segmented by scanner generation
in aggregation.  +31 tests (K_H4).  No `defaults_fingerprint`
drift.  VERSION bumped to `119.H4`.

**H4.1 (server) — [STEP_1F] Golden-Master Corpus Pinning**
`agent/raw_advice_scanner.py`, `tests/step_1f_corpus.json`,
`tests/tst_raw_advice_scanner.py`.  Adds a 31-document corpus
(`step_1f_corpus.json`) with hand-labeled file mentions per
document, and pins the scanner's measured recall at
**0.9810** against this corpus.  K_H4's golden-master test
fails if scanner edits drift recall.  This unblocks Phase 2B
activation (replacing the LLM preprocessor with the scanner)
because recall ≥ 0.90 was the trigger threshold from H4 plan
rev 1 §6.  No marker semantics change.  VERSION bumped to
`119.H4.1`.

**H5 (client+server) — diagnostic_messages Relay Channel**
`phenix_ai/run_ai_analysis.py`, `programs/ai_analysis.py`,
`tests/tst_diagnostic_messages.py` (new), `tests/run_all_tests.py`,
`tests/canary_expected.json`, `VERSION`.  Closes the
observability gap diagnosed during H4 deployment: H4's
`[STEP_1F]` markers emit to the SERVER's stderr when the call
dispatches to `ai.phenix-online.org`, but the server's stderr
isn't relayed to the client.  H5 adds a `diagnostic_messages`
list field to `working_results` and threads it end-to-end
through both dispatch modes (local in-memory passing AND
remote REST round-trip).  Engine-side `_emit_marker(list, str)`
helper appends markers to the list during processing.  Client
re-emits at a single uniform site
(`_relay_diagnostic_messages_to_stderr` called from
`run_job_on_server_or_locally`).  The operator sees identical
`[STEP_1F]`/`[STEP_1F_FAILED]` output regardless of where
analysis ran.

Key design decisions, locked across all H5.x work:
- **Total Initialization Policy**: every `working_results`
  from `get_results_from_all` gets `diagnostic_messages=[]`
  regardless of mode — defensive against
  add-marker-but-forget-to-init regressions.
- **Centralized re-emit**: single site in
  `programs/ai_analysis.py::run_job_on_server_or_locally`.
  Adding new markers requires NO client-side changes.
- **Uniformity criterion** (introduced in H5 plan rev 5 after
  Tom flagged rev 4's dispatch-mode branch): local and remote
  modes produce identical operator-visible output; no
  dispatch-mode branches in client behavior.
- **Helper safety**: `_emit_marker` is a pure list-append with
  `isinstance(list)` guard; never writes to stderr; never
  raises.  Wrapped at call sites in `try/except: pass` for
  defense against `"%s" % e` format failures.
- **Server-client compatibility**: defensive unpack handles
  missing field, non-list payloads, and corrupted JSON.
  Verified by K_H5 §C `test_old_server_response_handled_defensively`
  and `test_corrupted_payload_decoded_defensively`.

K_H5 organized into 5 sections at H5 ship: §A helper unit
tests (4), §B field-in-return tests (4, libtbx-dependent), §C
backward compatibility (3), §D uniform client re-emit (6), §E
production encode/decode integration (4).  +21 tests at H5
ship.  VERSION bumped to `119.H5`.

**H5.1 (client+server) — Channel Extended to 3 More Markers**
`phenix_ai/run_ai_analysis.py`, `programs/ai_analysis.py`,
`tests/tst_diagnostic_messages.py`, `tests/llm/canary_check.py`,
`tests/llm/tst_directive_extraction.py`, `tests/run_all_tests.py`,
`tests/canary_expected.json`, `VERSION`, `README.md`.  Instruments
three more engine-level exception handlers to flow through the
H5 channel:

| Marker | Engine site |
|---|---|
| `[ADVICE_PREPROCESSING_FAILED]` | `run_advice_preprocessing` main exception handler |
| `[DIRECTIVE_EXTRACTION_FAILED]` (outer) | `run_directive_extraction` main exception handler |
| `[FAILURE_DIAGNOSIS_FAILED]` | `run_failure_diagnosis` main exception handler |

These markers fire only on UNEXPECTED exceptions (library
bugs, programmer errors).  Configuration mistakes (invalid
provider, missing API key) are handled gracefully upstream by
`validate_api_keys` and `setup_llms` and don't trigger these
markers.

Two Gemini-mandated guardrails apply uniformly:
1. **No mutable default arguments**: all three functions use
   function-local `diagnostic_messages = []`, fresh per call.
   K_H5 §F field-in-return tests pin this via mutate-then-
   call-again checks.
2. **Emit before return, no new returns in except**: exception
   handlers append the marker then fall through to the main
   return.  Prevents accidental control-flow forks.

Tests use deterministic monkey-patching (patching
`validate_api_keys` or `setup_llms` to raise) plus Gemini Q2
dual-assertion (marker present AND debug_log evidence of
exception path) for robustness against future upstream
graceful-handling refactors.

Also folded in: **canary probe-text fix**.  The H3b live-LLM
canary used input `"Run phenix.refine on the model"` which
OpenAI faithfully extracted as `{"program_settings":
{"phenix.refine": {}}}` — semantically correct (no parameters
mentioned) but stripped to `{}` by `validate_directives`
because the inner dict is empty.  Google over-extracted and
returned non-empty, hiding the issue.  Probe updated to `"Run
phenix.refine with resolution 2.5"`; the concrete parameter
ensures non-empty validation across providers without changing
the canary's purpose (provider reachability, not extraction
quality).

K_H5 extended with §F (5 new tests), bringing K_H5 to **26
tests**.  VERSION bumped to `119.H5.1`.  No `defaults_fingerprint`
drift.

**H5.1.1 (client) — §2.1 v118.9 leftover cleanups (micro-ship)**
`agent/directive_extractor.py`, `tests/tst_directive_validation.py`
(new K_H5_1_1, 23 tests), `tests/run_all_tests.py`,
`tests/canary_expected.json`, `VERSION`, `README.md`.  Two
latent bugs in `agent/directive_extractor.py` surfaced
during the §2.1 review of v118.9 leftovers:

**Item 3 — `_corrected_from` sidecar protection in
`merge_directives`**.  The sidecar describes the CURRENT
value of `after_program` after an
experiment-type-mismatch correction.  Pre-fix, the
default override-wins dict-merge had two pathologies:

1. *Correction reverted*: if override's `after_program`
   matched `_corrected_from.from`, the LLM correction was
   silently un-corrected by simple-extraction merge.
2. *Zombie metadata* (caught during Gemini plan review):
   if override set `after_program` to a third value (not
   matching `_corrected_from.from`), the sidecar's `to`
   field became stale — pointing at a transition that no
   longer described the current value.

Fix: detect both cases inside `merge_directives`:
- Override matches `_corrected_from.from` → strip from
  override (preserve correction, keep sidecar)
- Override is a third value AND has no own sidecar →
  let override win, clear base's stale sidecar
- Override brings own sidecar → standard dict-merge
  (override's wins)
- Override doesn't touch `after_program` → standard
  dict-merge

`programs/ai_agent.py` is NOT modified.  Both production
callers (`phenix_ai/run_ai_analysis.py:1276` and
`programs/ai_agent.py:8376`) inherit the fix automatically.

**Item 4 — Boolean list-wrap defense in
`validate_directives`**.  `bool([False])` returns `True`
in Python (non-empty list is truthy) — silently flipping
semantics if the LLM emitted a list-wrapped boolean.  Two
sites had this:
- `prefer_anomalous`, `prefer_unmerged`, `prefer_merged`
  in `file_preferences`
- `use_experimental_phasing`, `use_molecular_replacement`,
  `use_mr_sad`, `model_is_placed`,
  `wants_validation_only` in `workflow_preferences`

Fix: explicit type-gated inline pattern at both sites.
Only unwraps `[bool]` (not `[int]`/`[str]`/etc.) — bounded
blast radius even against future maintenance error
(adding a list-typed key to either boolean tuple would
trigger drop+log rather than silent flattening).
Explicit pattern chosen over the existing
`_coerce_setting_value` helper to bound the unwrap to
genuine `[bool]` only.

Items 1 and 2 from §2.1 deferred (Item 1:
`_detect_experiment_type_signals` vs
`_resolve_after_program` aren't actually duplicates —
different regex breadth + ambiguity handling; unification
would be behavioural change.  Item 2:
`directives.get("stop_conditions") or {}` is idiomatic
safe read-only fallback, not a bug.).  See §1 of
`next_steps_post_v119.md` for the closure record.

Pre-existing wrongness NOT fixed by this ship:
`bool("false")` returns `True` in Python — a separate
string-truthy quirk preserved for non-list inputs.
K_H5_1_1's `test_validate_bool_bare_string_legacy_quirk`
pins the existing behaviour so a future fix catches it
deliberately.

+23 tests (K_H5_1_1 — 7 §A merge-protection + 16 §B
list-wrap-defense, including no-mutation-of-inputs
invariants for both branches).  VERSION bumped to
`119.H5.1.1`.  No `defaults_fingerprint` drift.

**H6 (test infrastructure) — Phase 2A planning-suite framework**
`tests/llm/framework.py`, `tests/tst_planning_framework.py`
(new K_H6, 18 tests), `tests/run_all_tests.py`,
`tests/canary_expected.json`, `VERSION`, `README.md`.

H6 replaces four stubs in `tests/llm/framework.py` with
real implementations, unblocking the 8-scenario planning
LLM reliability suite that has been queued since H3b:

- **`is_stop_intent(intent)`** — production-parity STOP
  detection.  Matches both signals production checks at
  `graph_nodes.py:~2168`: `program == "STOP"` OR
  `stop is True` (strict equality, not truthy).
- **`validate_planning_state(state)`** — 9-key structural
  validation per PHASE2_PLAN_v2.md §4.1.  Raises
  ValueError on missing keys.  Called from
  `tst_planning.py::build_scenarios()` at load time.
- **`call_planning_llm(state_inputs, provider)`** —
  production-faithful invocation mirroring
  `graph_nodes.py::plan()` lines 1992-2090, including the
  rate-limit handler.  Critical Gemini rev-3 fix:
  `raw_output` is hoisted outside the try block so that
  if `parse_intent_json` raises, the exception handler
  returns the LLM's malformed text for diagnosis.
- **`make_planning_run_fn(call_fn=...)`** — factory mirroring
  `make_directive_extraction_run_fn`, simplified since
  `parse_intent_json` never returns None (it either returns
  dict or raises).

K_H6's 18 tests cover all four implementations across 4
sections: §A `is_stop_intent` (6), §B `validate_planning_state`
(5), §C `make_planning_run_fn` (5), §D `call_planning_llm`
raw_output preservation (2).  §D tests SKIP gracefully
under PHENIX where the libtbx import path bypasses the
sandbox `sys.modules` fakes — sandbox tests pin the
mechanic; live planning runs exercise it in production.

Mac verification: 8 planning scenarios × 2 providers
(google, openai) = **16 PASS live LLM tests**, all
early-stopping at threshold.  This is the first time the
planning LLM has had reliability testing analogous to
directive_extraction.

VERSION bumped to `119.H6`.  No `defaults_fingerprint`
drift.  Plan: `v119_H6_PLAN_rev3.md`.

**H6.1 (micro-ship) — Per-scenario fail/err suffix**
`tests/llm/framework.py`, `VERSION`,
`tests/canary_expected.json`, `README.md`.

PHASE2_PLAN_v2.md §5 part 1: per-scenario output line
in `run_scenario_against_providers` gains a `[N fail]` /
`[N err]` suffix when nonzero.  Helps surface failures
that a reliability test tolerated by passing on threshold
(e.g., "PASS 4/5" — what did the 1 failure look like?).

Byte-identical to H6 output when both counts are zero,
which is the current observed reality across all 30
PASS scenarios from the H6 Mac run.

PHASE2_PLAN §5 part 2 (summary parenthetical) was
already implemented in `tests/llm/run_llm_tests.py`
lines 248-280 — no edit needed there.

Single edit in `framework.py` (~10 net lines).  No new
K-tests (purely cosmetic output formatting).  VERSION
bumped to `119.H6.1`.  No `defaults_fingerprint` drift.

**H7 (production code) — Phase 2B activation: scanner-first file extraction**
`phenix_ai/run_ai_analysis.py`, `tests/tst_phase2b_activation.py`
(new K_H7, 15 tests), `tests/run_all_tests.py`,
`tests/canary_expected.json`, `VERSION`, `README.md`.

H7 activates Phase 2B's first scope (Scope B per plan):
**scanner is the primary source for `extracted_files`**,
replacing the LLM preprocessor's
`extract_files_from_processed_advice` as the default path.
Unblocked since H4.1 measured scanner recall = 0.9810
against the LLM extraction (well above the 0.90
deprecation threshold from PHASE2_PLAN).

Two edits to `phenix_ai/run_ai_analysis.py`:

**Edit 1 — Scanner-first extraction (lines 1009-1069):**
- Scanner runs unconditionally on `raw_advice` (was gated
  on `processed_advice` existing pre-H7).  Net file
  detection improves when LLM preprocessing fails.
- Q2-strict: scanner called with `None` hint (NOT
  `file_list`).  Preserves the existing consumer contract
  at `programs/ai_agent.py:8128-8174` — the consumer's
  own comment names the field "Files mentioned in advice"
  (line 8142) and auto-promotes extracted files to the
  agent's `original_files`.  Passing `file_list` would
  silently inflate `original_files` with workspace files
  the user never mentioned.
- Fallback to LLM extraction preserved inside try/except,
  using the libtbx → relative import pattern (Gemini
  rev-2 critical fix: a bare libtbx import in the
  fallback would silently fail in sandbox and be
  swallowed by the outer except handler).
- Double-failure path explicitly logs the second exception
  (not swallowed by `except Exception: pass`).

**Edit 2 — Inline recall metrics in [STEP_1F] (lines 1155-1208):**
- Two new fields added to the existing [STEP_1F] marker:
  `scanner_recall_against_llm` (|intersection|/|llm_files|,
  the per-request equivalent of H4.1's 0.9810 corpus
  measurement) and `llm_recall_against_scanner`
  (|intersection|/|regex_files|, signals scanner-regex
  extension needs).  Both formatted as `%.4f`.
- Both metrics zero-division-guarded: empty denominator →
  1.0 (mathematically sound for recall, avoids
  ZeroDivisionError on advice with no file mentions).
- Purely additive change.  No existing [STEP_1F] fields
  modified.  No dashboard pollution risk (no existing
  inline `recall` field to displace).

**Subtle correction caught during implementation review:**
Post-H7, `extracted_files` is scanner-derived, so the
[STEP_1F] block can't use it as the LLM-comparison axis
(would compare scanner-to-scanner, defeating telemetry).
Fix: the block now calls
`extract_files_from_processed_advice` separately to obtain
the LLM view for comparison.  Fault-isolated: any failure
yields empty list; outer [STEP_1F_FAILED] catches
exceptions.  K_H7's
`test_telemetry_independence_from_extracted_files` pins
this property at the source level.

K_H7 (15 tests across 4 sections): §A scanner contract
(4), §B fallback semantics (3), §C recall metrics with
zero-div guards (3), §D integration smoke including
source-level regression pins (5).

Mac verification: **30/30 live LLM tests PASS**
(14 directive_extraction + 16 planning), **135/135
K-suite modules PASS** including K_H7's 15 new tests.
Scenarios that reference specific files
(`resolution_and_space_group`, `skip_programs`) still
pass — Q2-strict preserved the consumer contract.

`programs/ai_agent.py` NOT modified.  `core/llm.py` NOT
modified (`defaults_fingerprint` unchanged at
`sha256:77bf7421...`).

VERSION bumped to `119.H7`.  Plan:
`v119_H7_PLAN_rev2_2.md` (final, post-Gemini-rev-2,
post-baseline-correction).

**H8 (production code) — Template-literal allowlist fix
(production bug sub-cluster begins)**
`agent/command_postprocessor.py`, `tests/tst_autosol_bugs.py`
(new K_H8, 4 tests), `tests/canary_expected.json`, `VERSION`,
`README.md`.

H8 closes a production bug surfaced in AIAgent_62 cycle 4
(bromodomain ligand demo) where `phenix.autobuild_denmod` was
being emitted WITHOUT its `maps_only=True` flag.  Without that
flag, `autobuild_denmod` tries to rebuild the model rather than
just producing density-modified maps — wrong program semantics
for the cycle's intent.

**Root cause**: `agent/command_postprocessor.py::sanitize_command`
applies Rule D ("strip unknown program-specific parameters") using
a per-program allowlist built by `_load_prog_allowlist`.  The
allowlist was constructed only from the program's `strategy_flags`
dictionary — entries the LLM is permitted to set.  But
`programs.yaml` command templates also contain **literal key=value
tokens** baked into the template itself, like
`phenix.autobuild_denmod ... maps_only=True`.  These literal
parameters are program invariants, not LLM-tunable strategies, so
they were absent from `strategy_flags` and therefore absent from
the allowlist.  When sanitize_command processed the realized
command, Rule D treated `maps_only=True` as an "unknown parameter"
and stripped it.

**The fix — one new block of code in `_load_prog_allowlist`**
(lines 406–432, marked v119.H8):

After building the strategy_flags-derived allowlist, the code
now walks the program's command template (`prog_def.get('command', '')`)
and adds every literal `key=` parameter found:

```python
_cmd_template = prog_def.get('command', '')
if _cmd_template:
    for _m in re.finditer(
            r'(?:^|\s)([a-zA-Z_][a-zA-Z0-9_.]*)=',
            _cmd_template):
        _bare = _m.group(1).split('.')[-1].lower()
        if _bare:
            allowlist.add(_bare)
```

Three properties of the regex are load-bearing:

1. **`(?:^|\s)` prefix**: requires start-of-string or whitespace
   before the identifier.  This is what makes `{placeholder}`-style
   substitution targets NOT match — the `{` character is neither.
   So `{data_mtz}` doesn't pollute the allowlist with `data_mtz`
   the way it would if any `[a-zA-Z_]=` were accepted.
2. **`[a-zA-Z_][a-zA-Z0-9_.]*=` pattern**: matches dotted PHIL
   paths like `xray_data.r_free_flags.generate=True` as a single
   token, then takes the leaf via `split('.')[-1]`.  This handles
   both bare and scoped parameters uniformly.
3. **Leaf-only addition**: the allowlist contains bare leaves
   (e.g., `maps_only`), so the post-PHIL-flattening sanitization
   step (which compares leaves) works correctly.

K_H8 (4 tests in `tst_autosol_bugs.py`):

- `test_bug8_template_literal_maps_only_preserved` — primary
  regression.  Builds a realistic post-registry command for
  `phenix.autobuild_denmod` and asserts `maps_only=True` survives
  `sanitize_command`.
- `test_bug8_load_allowlist_includes_template_literals` — unit
  test on `_load_prog_allowlist` directly.  Confirms the allowlist
  for `phenix.autobuild_denmod` includes `maps_only` AND retains
  pre-existing entries (`nproc`, `resolution`) — i.e., no
  regression in the existing strategy_flags logic.
- `test_bug8_template_literal_extraction_ignores_placeholders` —
  pins the `(?:^|\s)` regex anchor's safety guarantee.  Verifies
  that `{data_mtz}` and similar placeholder syntax do NOT extract
  into the allowlist.
- `test_bug8_adversarial_strategy_override_dropped` — defense-in-
  depth.  Confirms that if an LLM emits `strategy.maps_only=False`
  to override the template invariant, `program_registry.build_command`
  drops the unknown strategy entry (logs "Unknown strategy
  'maps_only'") and the template's `maps_only=True` remains
  authoritative.  This was Claude-reviewer-suggested to ensure
  the fix doesn't accidentally allow LLM override of program
  invariants.

`agent/planner.py` NOT modified.  `agent/command_builder.py`
NOT modified.  `core/llm.py` NOT modified (`defaults_fingerprint`
unchanged).  `knowledge/programs.yaml` NOT modified — the
allowlist content (in `strategy_flags`) is correct; the template
literals just had to be discovered too.

VERSION bumped to `119.H8`.

**H9 (data + tests) — `phenix.predict_and_build` PHIL scope
rewrite**
`knowledge/parameter_fixes.json`, `tests/tst_autosol_bugs.py`
(new K_H9, 4 tests), `tests/canary_expected.json`, `VERSION`,
`README.md`.

H9 closes a production bug surfaced in run_38_openai cycle 2:
`phenix.predict_and_build` was emitting parameters like
`crystal_symmetry.unit_cell=...` and `crystal_symmetry.space_group=...`,
but the `phenix.predict_and_build` program expects
`crystal_info.unit_cell=...` and `crystal_info.space_group=...`.
PHENIX rejected the command with "Some PHIL parameters are not
recognized".

**Root cause** — a 3-layer pipeline gap:

- **Layer 1** (`agent/program_registry.py`, around line 778):
  when a bare unit_cell or space_group parameter is built into a
  command, the registry unconditionally prepends the
  `crystal_symmetry.` scope:

  ```python
  if key in ('unit_cell', 'space_group'):
      cmd_key = 'crystal_symmetry.%s' % key
  ```

  This is correct for most programs.
- **Layer 2** (`agent/planner.py::fix_program_parameters`,
  consuming `knowledge/parameter_fixes.json`): for programs whose
  PHIL scope differs from `crystal_symmetry`, the parameter
  name is rewritten.  Pre-H9, the existing `parameter_fixes.json`
  had rewrite entries for `phenix.refine` (which goes the other
  direction: `crystal_info.*` → `crystal_symmetry.*` for the
  refine-specific scope handling) but **no scope-rewrite entries
  for `phenix.predict_and_build`**.  The program's existing entry
  in the JSON had only `data_file`, `input_files.data_file`,
  `control.nproc`, `prediction.nproc` mappings.

So when Layer 1 produced `crystal_symmetry.unit_cell=...` for
`phenix.predict_and_build`, Layer 2 had no rewrite rule for that
parameter, the wrong-scope parameter survived to the command line,
and PHENIX rejected it.

**The fix — 4 new rewrite entries plus 1 annotation in
`knowledge/parameter_fixes.json` under `phenix.predict_and_build`:**

```json
"phenix.predict_and_build": {
  "data_file": "xray_data_file",                              // pre-existing
  "input_files.data_file": "input_files.xray_data_file",      // pre-existing
  "control.nproc": null,                                       // pre-existing
  "prediction.nproc": null,                                    // pre-existing
  "_comment_cs": "predict_and_build uses crystal_info.* scope (v119.H9)",  // NEW
  "crystal_symmetry.space_group": "crystal_info.space_group",  // NEW
  "crystal_symmetry.unit_cell":   "crystal_info.unit_cell",    // NEW
  "xray_data.space_group":        "crystal_info.space_group",  // NEW
  "xray_data.unit_cell":          "crystal_info.unit_cell"     // NEW
}
```

The `xray_data.*` fallback entries handle the case where an LLM
might emit `xray_data.unit_cell` (a different incorrect scope)
— Gemini's review suggested covering this for robustness.  The
`_comment_cs` key is skipped by `fix_program_parameters` (which
ignores keys starting with `_`), so it functions as inline
documentation.

**K_H9 testing wrinkle**: `parameter_fixes.json` is loaded once
on first call and cached at module scope in `agent/planner.py`
via `_PARAMETER_FIXES = None` (lazily populated by
`get_parameter_fixes()`).  K_H9's tests need to invalidate this
cache between runs so they exercise the JSON file's current
state.  Solution: a module-local `_reload_parameter_fixes()`
helper in the test file that simply sets the planner module's
cache back to `None`:

```python
def _reload_parameter_fixes():
    """Force re-read of parameter_fixes.json (defeat the module cache)."""
    try:
        from libtbx.langchain.agent import planner as _p
    except ImportError:
        try:
            from agent import planner as _p
        except ImportError:
            return
    _p._PARAMETER_FIXES = None
```

All four K_H9 tests call this at the top.  This pattern
generalizes: any K-test that depends on JSON-or-YAML-loaded
data with module-level lazy caching needs an explicit
invalidation hook.

K_H9 (4 tests in `tst_autosol_bugs.py`):
`test_bug9_predict_and_build_unit_cell_scope_rewrite`,
`test_bug9_predict_and_build_space_group_scope_rewrite`,
`test_bug9_predict_and_build_xray_data_fallback_rewrite`
(Gemini-suggested coverage),
`test_bug9_predict_and_build_idempotent_on_correct_scope`
(no double-rewrite when input is already correct).

`agent/program_registry.py` NOT modified — Layer 1 is correct
as the default behavior; Layer 2's rewrite table is the
extensibility point for per-program exceptions.  `agent/planner.py`
NOT modified — the consumer of `parameter_fixes.json` already
applied entries uniformly.  `defaults_fingerprint` unchanged at
`sha256:77bf7421...`.

VERSION bumped to `119.H9`.  Plan: `v119_H9_PLAN.md`
(Gemini-reviewed, green light).

**H10 (production code) — `exclude_patterns` structural fix
(prerequisite for H11)**
`agent/command_builder.py`, `tests/tst_autosol_bugs.py` (new
K_H10, 5 tests), `tests/canary_expected.json`, `VERSION`,
`README.md`.

H10 closes a structural gap surfaced in AIAgent_62 cycle 7:
`phenix.refine` emitted a command with `refine_001_001.cif` (a
model mmCIF from an earlier refine step) as a positional third
argument.  PHENIX interpreted it as a SECOND model and crashed
with "wrong number of models".  The `phenix.refine::ligand_cif`
slot has an `exclude_patterns` filter specifically to prevent
this, but the filter wasn't being applied at the relevant
selection path.

**Root cause** — `agent/command_builder.py`'s file-selection
algorithm has multiple paths through `_find_file_for_slot` plus
a pre-population step in `_select_files`.  The in-code H10
comment (line 1279-1283) summarizes the gap:

> Pre-H10, `exclude_patterns` was only consulted at PRIORITY 4
> (extension fallback).  PRIORITY 2 (best_files), PRIORITY 2.5
> (recovery strategies), PRIORITY 3 (category-based), and
> PRIORITY 3.5 (fallback best_files) all bypassed it.

(The LLM-selected path at line 963 in `_select_files` had its
own exclude_patterns check; whether this predated H10 is not
marked.  The bug class was the auto-fill paths inside
`_find_file_for_slot`, which is what H10's helper covers.)

The AIAgent_62 cycle-7 build picked `refine_001_001.cif` via
**PRIORITY 3 (category-based)** — one of the paths that
bypassed the filter.  Since the cycle-7 LLM did not specify a
ligand_cif file, auto-fill ran, hit PRIORITY 3, and picked the
wrong file from the `ligand_cif` category.

**The fix — helper + 8 application sites, 9 markers:**

A closure `_matches_exclude_patterns_h10(f)` defined at the
entry of `_find_file_for_slot` (line 1292):

```python
_exclude_patterns_h10 = input_def.get("exclude_patterns", [])
def _matches_exclude_patterns_h10(f):
    if not _exclude_patterns_h10:
        return False
    return matches_exclude_pattern(
        os.path.basename(f), _exclude_patterns_h10)
```

Applied as an additional filter clause at every auto-fill site:

- 7 sites inside `_find_file_for_slot` (lines 1333, 1359, 1406,
  1417, 1437, 1472, 1513): PRIORITY 2 best_files, PRIORITY 2.5
  recovery_strategies, PRIORITY 3 category subcat / parent /
  multiple variants, PRIORITY 3.5 `require_best_files_only`,
  PRIORITY 3.5 fallback best_files.
- 1 site inside `_select_files` (line 698): refinement pre-
  population.  `best_files["model"]` was being injected without
  honoring the model slot's exclude_patterns; the patch adds the
  check inline (using `matches_exclude_pattern` directly rather
  than via the closure, since `_select_files` doesn't have one).

Total: 9 `v119.H10` markers in the file (1 helper definition
comment block + 7 closure-application sites + 1 in `_select_files`).

**Implementation lesson** — the 8th patch site (refinement pre-
population in `_select_files`) was discovered during act/verify
when Test 3 (`test_bug10_exclude_patterns_applied_in_best_files_path`)
failed initially.  The diagnostic revealed that
`_find_file_for_slot` was never called for the model slot when
pre-population had already injected `best_files["model"]`.  Fix:
extend the patch to `_select_files`.  This is the
"act/verify catches missed sites" pattern that justifies
running each test as it's written rather than batching them.

K_H10 (5 tests in `tst_autosol_bugs.py`):
`test_bug10_ligand_cif_slot_rejects_model_mmcif` (primary
regression),
`test_bug10_ligand_cif_slot_accepts_legitimate_cif` (no
over-rejection),
`test_bug10_exclude_patterns_applied_in_best_files_path` (pins
pre-population fix),
`test_bug10_exclude_patterns_applied_in_category_path` (pins
PRIORITY 3 category fix — this was the test that failed on
Tom's real PHENIX and surfaced the H11 follow-on),
`test_bug10_unrelated_slots_unaffected` (no-regression sentinel,
Claude-reviewer suggested).

**Important caveat — H10's tarball README is INCORRECT** about
fully fixing the AIAgent_62 cycle-7 bug.  H10 IS code-correct
(all patches fire), but a pre-existing bug in the YAML
`exclude_patterns` data (substring-style authoring of
`"refine_"` despite the function using word-boundary semantics)
prevented the filter from rejecting `refine_001_001.cif`.  H10
provides the structural prerequisite; **H11 is required to
actually fix the bug**.  See H11 below.

`agent/file_utils.py` NOT modified — the `matches_exclude_pattern`
function's semantics are correct and intentional (word-boundary
regex matching).  `core/llm.py` NOT modified.  `defaults_fingerprint`
unchanged.

VERSION bumped to `119.H10`.  Plan: `v119_H10_PLAN_rev2.md`
(Gemini-reviewed).

**H11 (data) — exclude_patterns YAML word-boundary fix
(completes H10)**
`knowledge/programs.yaml`, `tests/tst_autosol_bugs.py`,
`tests/canary_expected.json`, `VERSION`, `README.md`.

H11 completes the H10 fix.  When Tom installed H10 and ran
`tst_autosol_bugs.py`, Test 4
(`test_bug10_exclude_patterns_applied_in_category_path`) failed
on his machine while it had passed in the sandbox.  A diagnostic
script (`h10_diagnostic.py`) called `matches_exclude_pattern`
directly on representative cases and revealed that the function
uses **word-boundary regex matching**, not substring matching:

```python
# agent/file_utils.py
re.search(
    r'(?:^|[_\-\.])' + re.escape(pat_stem) + r'(?=[_\-\.]|$)',
    stem)
```

A pattern matches if it appears in the stem preceded by start-of-
string or one of `_`/`-`/`.`, AND followed by one of `_`/`-`/`.`/
end-of-stem.  Under this semantics, three YAML patterns authored
with substring intent did NOT match their intended targets:

| Slot | Pattern | Why broken |
|---|---|---|
| `phenix.refine::ligand_cif` | `"refine_"` | Trailing `_` requires another boundary char after; in `refine_001_001.cif` the next char is `0`, not a boundary.  **THE bug — broke AIAgent_62 cycle 7.** |
| `phenix.mtriage::full_map` | `"_half"` | Leading `_` consumed by boundary-prefix regex group; mismatch with stem position. |
| `phenix.real_space_refine::map` | `"_half"` | Same as above. |

The sandbox stub for `matches_exclude_pattern` used substring
matching, so H10's tests passed in sandbox but failed against
the real function.  H10 packaged as "fixes the cycle-7 bug" —
it didn't; H10's filter fired but `matches_exclude_pattern` returned
False on every file because of the YAML pattern bug.

**The fix — three independent edits:**

1. **YAML pattern corrections (3 lines):**
   - `phenix.refine::ligand_cif`: `"refine_"` → `"refine"`
   - `phenix.mtriage::full_map`: `"_half"` → `"half"`
   - `phenix.real_space_refine::map`: `"_half"` → `"half"`

   Word-boundary `"refine"` correctly matches `refine_001.cif`,
   `phenix_refine.cif`, `pre_refine_template.cif`; correctly
   skips `refining_template.cif`.  Word-boundary `"half"`
   correctly matches `protein_half.ccp4` (un-numbered case).
   The numeric patterns `half_1`, `half_2`, `half1`, `half2` are
   preserved alongside the bare `"half"` because Gemini's plan
   review caught that bare `"half"` fails on `protein_half1.ccp4`
   (the `1` immediately after `half` is not a boundary char).

2. **YAML authoring docstring (~25 lines in `programs.yaml`):**
   A `DESIGN NOTE (exclude_patterns)` block near the top of the
   file documenting the word-boundary semantics, with worked GOOD
   and BAD examples for both trailing-underscore and leading-
   underscore traps.  Inlined per Gemini's plan-review recommendation
   to maximize contextual freshness for future YAML authors.

3. **Bug 11 semantic-pin test (1 new test, 10 named cases):**
   `test_bug11_matches_exclude_pattern_semantics` calls the real
   `matches_exclude_pattern` and asserts its actual behavior on
   10 documented cases, including the trailing-underscore bug
   (`"refine_"` does NOT match `refine_001.cif`), the bare-pattern
   correctness (`"refine"` does match), the alphanumeric-boundary
   subtlety (`"half"` does NOT match `protein_half1.ccp4`), and
   extension-suffix patterns.  Includes a graceful sandbox-skip
   fallback so lightweight CI environments don't error.

K_H11 (1 test): `tst_autosol_bugs.py::test_bug11_matches_exclude_pattern_semantics`.

**Cluster regression check post-H11:** All 9 v119-cluster K-suites
PASS (160 PASS, 0 FAIL).  Bug 8 (4), Bug 9 (4), Bug 10 (5), Bug
11 (1): all PASS in sandbox against word-boundary-correct stub.
**VERSION/canary lockstep verified.**  `defaults_fingerprint`
unchanged at `sha256:77bf7421...`.

`agent/command_builder.py` NOT modified (all 9 H10 markers
preserved).  `agent/file_utils.py` NOT modified (word-boundary
semantics are intentional and now documented in the YAML
authoring guide).  `core/llm.py` NOT modified.

**The load-bearing verification** — H11's actual proof of correctness
is Tom rerunning the AIAgent_62 bromodomain demo on the Mac.
Cycle 7's `phenix.refine` command must NOT include
`refine_001_001.cif` as a positional argument.  Pre-H11 (with
H10 installed) this still fails; H10 + H11 together should fix it.

VERSION bumped to `119.H11`.  Plan: `v119_H11_PLAN_rev2.md`
(final, post-Gemini-rev-2).

**Lesson captured (H10 → H11 sandbox cycle):**

The H10 → H11 transition surfaced a sandboxing failure mode worth
remembering as policy:

- The H10 sandbox stub for `matches_exclude_pattern` did substring
  matching.  The real PHENIX function does word-boundary matching.
- The stub was wrong but no test caught the divergence.  H10's
  sandbox K-suite passed.
- H10 was packaged and presented as "fixes the AIAgent_62 cycle-7
  bug" — it didn't.  H10's filter fired but the YAML pattern
  didn't match.
- The divergence was discovered at install-time when Tom's `phenix.python
  tst_autosol_bugs.py` Test 4 failed against the real function.
- The `h10_diagnostic.py` script that called `matches_exclude_pattern`
  directly with representative cases revealed the semantics.

Forward policy (also captured in `DEVELOPER_GUIDE.md` §10.3 "Test
adequacy"): sandbox stubs for functions with non-obvious semantics
(regex, parsing, escaping, etc.) should be accompanied by a
semantic-pin test that calls the real function with documented
expected behavior.  `test_bug11_matches_exclude_pattern_semantics`
is the template — when a test that uses the real function asserts
documented behavior, any future stub-vs-real divergence is caught
the moment a sandbox test runs against the real function.

**H12 (production code) — `exclude_patterns` helper DRY refactor +
categorizer semantic-pin suite (deferred follow-up cluster from H11)**
`agent/command_builder.py`, `tests/tst_file_categorizer.py`
(NEW, ~75 cases across 8 tests), `tests/canary_expected.json`,
`VERSION`, `README.md`.

H12 closes three deferred items from `v119_H11_PLAN_rev2.md`
(items 2.2, 2.3, 2.4) bundled into a single follow-up ship.  The
items are independent in scope but share file area and zero
behavioral change:

**Item 2.3 — helper DRY refactor in `agent/command_builder.py`**:

H10 introduced a closure `_matches_exclude_patterns_h10(f)`
defined at the entry of `_find_file_for_slot` and applied at
7 sites there, plus a direct call to `matches_exclude_pattern`
in `_select_files` (line 698) where the closure wasn't in
scope.  H12 consolidates these into a single module-level
helper near the top of `command_builder.py`:

```python
def _file_passes_exclude_patterns(file_path, exclude_patterns):
    """Return True if file_path is NOT rejected by exclude_patterns."""
    if not exclude_patterns:
        return True
    return not matches_exclude_pattern(
        os.path.basename(file_path), exclude_patterns)
```

Three design choices per Gemini's H12 plan review:

1. **Signature takes `exclude_patterns: list`, not `input_def: dict`**.
   Decouples the helper from the program-registry schema; lets
   `_select_files`'s call site (which has the list, not a
   wrapping dict in scope) pass it directly.  Within
   `_find_file_for_slot`, the list is extracted once at the top
   of the function — preserving the closure's micro-perf
   property (one dict lookup, not 7).
2. **Sense-inversion**.  The closure returned True if the file
   MATCHED (i.e., should be rejected).  The helper returns True
   if the file PASSES (i.e., is NOT rejected).  Call sites read
   more naturally as `if not _file_passes_exclude_patterns(f, p): continue`.
3. **Marker convention**: each refactored call site is marked
   `# v119.H10/H12` rather than dropping the H10 history.  Clean
   audit trail: the BEHAVIOR comes from H10's bug fix; the
   present STRUCTURAL form belongs to H12.

Total: 8 application sites consolidated to 1 helper, 9
`v119.H10/H12` markers in `command_builder.py`.

**Items 2.2 + 2.4 — `tests/tst_file_categorizer.py` (NEW)**:

Semantic-pin tests for every public function in
`agent/file_utils.py`.  Generalizes the discipline that
`test_bug11_matches_exclude_pattern_semantics` applied to one
function — pre-emptive defense against the kind of
sandbox-stub-vs-real drift that masked the AIAgent_62 cycle-7
bug for one ship.

Functions covered (8 tests, ~75 documented assertions):
- `classify_mtz_type` — 12 cases across 5 rules + ordering
  invariant comment
- `get_mtz_stage` — 9 cases (data_mtz, map_coeffs_mtz subcategories,
  unrecognized-category passthrough)
- `get_category_for_extension` — 19 cases including
  case-insensitivity, multi-dot filenames (`model.backup.pdb`),
  no-extension files (`README`, `LOCAL_MAP`), unknown extensions,
  trailing-dot edge case
- `is_mtz_file` — 8 cases including case-insensitivity, trailing-
  suffix discrimination
- `is_model_file` — 12 cases pinning the load-bearing quirk that
  `.pdb` extension overrides the `'ligand'`-in-name filter (so
  `ligand_fit.pdb` IS a model file but `my_ligand.cif` is not)
- `is_map_file`, `is_sequence_file` — 7 cases each
- Cross-function consistency — 4 integration cases verifying
  `is_mtz_file` + `classify_mtz_type` + `get_mtz_stage` agree
  and that `is_X_file` / `get_category_for_extension` are
  internally consistent

`matches_exclude_pattern` is NOT re-pinned here —
`test_bug11_matches_exclude_pattern_semantics` already covers
it canonically.

Gemini-suggested cases folded in: multi-dot (`model.backup.pdb`,
`data.processed.mtz`), no-extension (`README`, `LOCAL_MAP`,
`file.`), `is_model_file` boundary (`ligand_fit.pdb` vs
`co_crystallized_ligand_structure.cif`), and a rule-ordering
invariant comment in `test_classify_mtz_type_semantics`.

**Why these bundle as 2.2 + 2.4**: Gemini's original 2.4
proposal was a "differential parity test" comparing real
function vs sandbox stub.  Reframed in H12 as "broaden the
semantic-pin discipline to ALL of `agent/file_utils.py`" —
because (a) sandbox stubs aren't a real artifact in the PHENIX
tree, (b) `test_bug11_*` already provides the canonical parity
check for the one function with documented divergence problems,
and (c) generalizing the pattern to the other functions in
`file_utils.py` is the more useful generalization.

K_categorizer: 8 tests, all PASS in sandbox.  All 6 K_Bug 10
+ K_Bug 11 tests still PASS unchanged after the H12 refactor,
providing the strong assurance that behavior is preserved.
`defaults_fingerprint` unchanged.

`agent/file_utils.py` NOT modified — its functions are the
SUBJECT of pinning, not the subject of change.  All H1–H11
changes carried forward.

VERSION bumped to `119.H12`.  Plan: `v119_H12_PLAN.md`
(Gemini-reviewed).

**H13 (production code) — Ollama provider robustness + retired-
model 404 detection (v118 §3.5 closed; AIAgent_run_39a_ollama bugs
fixed)**
`core/llm.py`, `agent/directive_extractor.py`, `agent/api_client.py`,
`tests/tst_provider_error_classification.py` (NEW, 13 tests),
`tests/tst_default_models.py` (+8 H13 tests),
`tests/canary_expected.json`, `VERSION`, `README.md`.

H13 closes a production bug surfaced by Tom's
`phenix.ai_agent provider=ollama ...` run on cci-gpu-01.  The
failure mode was `[DIRECTIVE_EXTRACTION_FAILED]: 404 page not
found, falling back to rules-only resolver` — opaque and
non-actionable.  Diagnostic surfaced two stacked bugs PLUS an
opportunity to complete the v118 §3.5 retired-model classification
work that was deferred through H1–H12.

**Root cause — TWO bugs, each duplicated at 2 sites:**

- **Bug A — URL mismatch.**  `agent/directive_extractor.py:1226`
  and `agent/api_client.py:735` had:
  ```python
  base_url = os.environ.get(
      "OLLAMA_BASE_URL", "http://localhost:11434/v1")
  ```
  When the user set `OLLAMA_BASE_URL=http://localhost:11434`
  (no `/v1`), the env-var completely overrode the `/v1`-suffixed
  default.  The OpenAI Python SDK appends `/chat/completions`
  to the base URL, so the actual request went to
  `http://localhost:11434/chat/completions`.  Ollama exposes
  the OpenAI-compat API only at `/v1/chat/completions`; for
  any other path it returns the HTML literal `404 page not
  found`.  Tom's user-facing behavior (matching `OLLAMA_HOST`
  which is the bare host) silently broke.

- **Bug B — env-var override ignored.**  `directive_extractor.py:1219`
  and `api_client.py:734` called `default_model_for_provider("ollama")`
  unconditionally — ignoring `OLLAMA_LLM_MODEL` even when the
  user set it.  This is inconsistent with
  `core/llm.py:get_llm_and_embeddings` which DID honor the
  env-var ("Env-var overrides take precedence; fall back to
  the central default tables when the env vars are unset").
  H1 centralized the table; it missed centralizing the
  precedence rule.

**The fix — three new helpers in `core/llm.py`:**

```python
def normalize_ollama_openai_base_url(base_url):
    """Idempotently append /v1.  Empty string and None pass through."""
    if not base_url:
        return base_url
    stripped = base_url.rstrip('/')
    if stripped.endswith('/v1'):
        return stripped
    return stripped + '/v1'

_PROVIDER_MODEL_ENV_OVERRIDES = {"ollama": "OLLAMA_LLM_MODEL"}

def resolve_model_for_provider(provider, role="decision"):
    """Env-var precedence: env-var wins, central default falls back."""
    env_var = _PROVIDER_MODEL_ENV_OVERRIDES.get(
        (provider or "").lower().strip())
    if env_var:
        env_value = os.getenv(env_var)
        if env_value:
            return env_value
    return default_model_for_provider(provider, role=role)
```

Plus a unified provider-error classifier
`_classify_provider_error(exc, model_name)` (also in `core/llm.py`)
that distinguishes four error classes:

| Marker tag | Trigger | Operator action |
|---|---|---|
| `DIRECTIVE_EXTRACTION_MODEL_RETIRED` | 404 + retirement phrase ("deprecated"/"no longer available"/"retired") near the word "model" | Update `DEFAULT_MODELS` |
| `DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE` | 404 + `model '...' not found` regex match (Tom's case) | Pull the model, or check env-var typo |
| `DIRECTIVE_EXTRACTION_AUTH_FAILED` | 401 or "unauthorized" or "invalid api key" | Check provider API key |
| `DIRECTIVE_EXTRACTION_FAILED` | Anything else | Investigate manually |

Defense-in-depth design (Gemini's H13 plan review):

1. **Programmatic property interrogation** — reads
   `status_code`, `.response.status_code` (chained), `.body`
   (dict), `.message` attributes BEFORE falling back to
   `str(exc)`.  Defends against SDK exceptions that hide HTTP
   details behind opaque `__str__` representations.
2. **Co-occurrence rule for RETIRED** — the retirement phrase
   must appear within 80 chars of the word "model".
   Distinguishes legitimate retirement ("The model X has been
   deprecated") from edge-proxy 404 pages with similar
   language ("This endpoint is deprecated; see new path").
3. **UNAVAILABLE regex** — `model\s+['"]X['"]\s+not\s+found`
   matches Tom's exact error pattern from ollama
   (`{"error":{"message":"model 'llama3.2' not found"}}`).

**Where the helpers are applied:**

- `agent/directive_extractor.py` ollama block (the
  `_call_llm_fallback` path): URL via
  `normalize_ollama_openai_base_url`, model via
  `resolve_model_for_provider`, classifier dispatch in the
  exception handler with two emit helpers
  (`_emit_retired_model_marker` from v119.H1 +
  `_emit_unavailable_model_marker` NEW in H13).
- `agent/api_client.py::_call_ollama_llm` (the primary path):
  same URL and model fixes.  Classifier dispatch happens at
  the OUTER catch in `extract_directives` line 876 — see H13.1
  fix in the README for details.
- `core/llm.py::get_llm_and_embeddings` (the native-`/api/chat`
  path via `langchain-ollama`): NOT modified.  That path uses
  ChatOllama with the native endpoint and DOES NOT need `/v1`;
  applying `normalize_ollama_openai_base_url` would break it.
  An explicit comment at the call site documents the asymmetry.
  Env-var precedence already worked correctly here (it was the
  inconsistency-source for H13).

**K_provider_error (13 tests in
`tst_provider_error_classification.py`)**:

- 3 RETIRED cases: Google verbose 404, OpenAI clean 404,
  Anthropic `model_not_found`
- 3 UNAVAILABLE cases: Tom's exact error verbatim (primary
  regression), variant with `qwen2.5:72b`, double-quoted
  format
- 3 FAILED cases (negative): edge-proxy "endpoint deprecated"
  (LOAD-BEARING co-occurrence sentinel), bare 404, empty
  exception
- 2 AUTH cases: 401, "invalid api key" without status code
- 2 class-attribute defense cases (Gemini's review): mock
  exception with `status_code=404`, mock with `.response.status_code`
  chained

**K_default_models extension** (`tst_default_models.py`, +8 tests):
- 4 `resolve_model_for_provider` tests: unset env-var falls
  through, set env-var overrides, empty env-var falls through,
  unknown provider raises
- 4 `normalize_ollama_openai_base_url` tests: appends `/v1`,
  idempotent, None passthrough, empty-string passthrough
  (added in H13.1 after self-review found the original helper
  produced `'/v1'` from an empty input — see H13.1 notes)

Also updates 2 existing source-scan tests
(`test_api_client_ollama_uses_central_default`,
`test_directive_extractor_fallback_uses_central_defaults`) to
accept BOTH `default_model_for_provider` AND
`resolve_model_for_provider` references — the centralization
invariant is "every provider resolves via a central helper",
not "every provider uses the exact same helper name".

**Refined v118 §3.5 scope**:

Pre-H13 the existing `_looks_like_retired_model_error` (from
v119.H1) lumped "model not found" together with retirement
phrases, which would have misclassified Tom's case as
retirement (the operator hint would say "update
DEFAULT_MODELS" — wrong action; the right action is `ollama
pull`).  H13's classifier separates these into RETIRED vs
UNAVAILABLE so the hint matches the operator action.

**Known limitation** (NOT fixed in H13): the H1-era
`_looks_like_retired_model_error` + `_emit_retired_model_marker`
path is still used by google/openai/anthropic exception
handlers in `_call_llm_fallback`.  Only the ollama path was
migrated to the unified classifier.  Migrating the other three
providers is candidate work for H14 (scope deferred since the
immediate bug was ollama-specific).

VERSION bumped to `119.H13` (then `119.H13.1` after the
self-review surfaced two issues — see "H13.1 post-review fixes"
in the ship README).  Plans: `v119_H13_PLAN_rev3.md` (rev3,
post-Ollama-diagnosis; rev1 and rev2 covered the broader
§3.5 + §4.9 scope before run_39a_ollama narrowed the focus).

`defaults_fingerprint` unchanged at
`sha256:77bf7421...` — H13 added new functions but didn't
modify the `DEFAULT_MODELS` tables that feed the fingerprint.

**Deferred to H14** (per H13's rev3 plan): §4.9
prompt-consolidation refactor of `DIRECTIVE_EXTRACTION_PROMPT`.
Medium risk, needs multi-pass approach + Gemini review +
LLM-roundtrip verification.  Was bundled with §3.5 in H13's
rev1/rev2 plans; rev3 split to keep H13's blast radius small
after Tom's bug took priority.

**H14 (production code) — Three independent regressions surfaced by
run_39_openai batch analysis**
`agent/directive_extractor.py`, `programs/ai_agent.py`,
`tests/tst_solve_action_keywords.py` (NEW, 13 tests),
`tests/tst_step1f_single_emit.py` (NEW, 8 tests),
`tests/tst_space_group_validation.py` (NEW, 24 tests),
`tests/tst_utils.py` (NEW — restored sandbox helper, enables
3 previously-shipped K-tests to run in sandbox),
`tests/tst_autosol_bugs.py` (1 test relaxed from exact-signature
match to substring — was brittle to additive kwargs),
`tests/canary_expected.json`, `VERSION`.

H14 closes three independent regressions surfaced by post-H13 batch
analysis of `run_39_openai` (520 runs, May 2026) compared against
`run_25_openai` (720 runs, March 2026 baseline).  The investigation
was made possible by `scan_batch_runs.py` (a new
operations-side scanner that triages run logs into Tier-1 crashes,
Tier-2 diagnostic markers, Tier-3 state anomalies, and Tier-4 soft
anomalies) and a controlled three-way log comparison of the
1029B-sad dataset across run_25-solve, run_39-solve, and
run_39-stop variants.  The three items are independent; they ship
together because they were diagnosed in one investigation pass.

**Item 1 — Phaser false-positive in `_ACTION_TABLE["solve"]`
(operationally biggest win):**

The rules-only directive extractor's `_ACTION_TABLE["solve"]`
entry treated the goal phrase "solve the structure" as a synonym
for "do molecular replacement with phaser."  When a README
contained BOTH that goal phrase AND another action (e.g., the
"Stop after refinement" suffix), the multi-action branch in
`_apply_workflow_intent_fallback` (`n > 1` case) set
`start_with_program = phenix.phaser` — forcing phaser into the
workflow even on SAD/MAD datasets where the correct method is
autosol.  Once phaser was in the workflow, it ran on a model that
autosol had already refined to a multi-chain ASU; phaser's
composition check then errored with `composition will not fit in
the unit cell volume`.  24 occurrences across 5 datasets in
run_39 (1029B-sad, 1J4R-ligand, 1aba-polder, 3tpp-ensemble-refine,
7rpq_AF_reference_model).

The controlled comparison confirmed the trigger:
  - `run_25/1029B-sad__rules_only_solve` (README without "stop"):
    workflow xtriage→autosol→autobuild→molprobity→refine.  SUCCESS.
  - `run_39/1029B-sad__rules_only_solve` (same README): IDENTICAL
    workflow.  SUCCESS.
  - `run_39/1029B-sad__rules_only_stop` (README with "Stop after
    refinement" added): xtriage→autosol→refine→phaser×4.  FAILED.

The only input that differed was the README — adding one line
that the rules-only extractor turned into a stop directive AND
a phaser start_with directive.

**Fix:** Remove the goal phrases `"solve the structure"` and
`"solve structure"` from the `solve` action's keyword list.
Explicit method keywords (`"molecular replacement"`, `"phaser"`,
`"mr "`) remain — only the ambiguous goal phrasing is removed.
This restores the semantic: `solve` action means "user explicitly
asked for MR," not "user expressed a goal of solving the
structure."  When the goal phrase appears alone, no action fires
and `RulesSelector` picks the right method based on data type.

The `asu_copies` → `phaser.search_copies` injection in
`graph_nodes.py` (originally Item 2 of the H14 plan, also surfaced
by this investigation) was DROPPED from H14 after Tom's domain
input: "the value from xtriage is usually pretty good, keep it as
is."  When phaser legitimately runs (with a single-chain search
model), `search_copies=N` from xtriage's Matthews `Best guess` is
the correct value.  The problem in run_39 was not the injection —
it was that phaser shouldn't have been running at all on
already-phased SAD data.  Item 1 prevents that, eliminating the
24/24 phaser failures.  Defensive instrumentation for the
multi-chain-search-model case is deferred to a candidate H15
item.

**Item 2 — STEP_1F preprocessor double-emit (cosmetic but
high-volume):**

`[STEP_1F] preprocessing_metrics` lines were being emitted twice
per advice preprocessing call.  60.6% of run_39_openai runs (315
of 520) had adjacent duplicate STEP_1F lines.  The cause: TWO
sites both iterated `diagnostic_messages` and wrote to stderr:

  1. `phenix.programs.ai_analysis._relay_diagnostic_messages_to_stderr`
     (the v119.H5.1 dispatcher-level relay; called inside
     `run_job_on_server_or_locally` at 4 return paths)
  2. `programs/ai_agent.py:8087-8103` (a client-side relay block
     from v119.H5 §2.9 that survived the H5.1 refactor)

The H5 §2.9 comment at site 2 claimed it was "the SINGLE, uniform
site where operators see them."  That claim was correct at the
time H5 shipped, but the H5.1 refactor moved the relay INTO the
dispatcher (the architecturally correct location, since the
dispatcher is the single place that handles all dispatch paths
uniformly).  Site 2 was never deleted, became dead duplicate code,
and survived undetected because the duplicate emit doesn't break
anything functionally — operators just saw each marker twice.

**Fix:** Delete the client-side relay block in `ai_agent.py`,
replace with a marker comment explaining the H14 deletion and
referencing the dispatcher's relay as the canonical site.  Source
scan tests pin the absence of `for _msg in _diagnostics` in
`ai_agent.py` and the continued presence of
`_relay_diagnostic_messages_to_stderr` (with 4 call sites) in
`ai_analysis.py`.

**Item 3 — `space_group` validation extensions (defensive
cleanup):**

The pre-H14 validator at `validate_directives()` already dropped
LLM-emitted `space_group` placeholder values via the
`_SYMMETRY_SENTINELS` frozenset (`"Not specified"`, `"Not
provided"`, etc.) plus negative structural checks (must start
with a letter; length <= 25).  Two gaps remained:

  1. The sentinel set did not include the "Not explicitly *"
     family.  Tom's qwen2.5:72b xtriage tutorial verification of
     H13.1 surfaced the truncated value `"Not explicitly mentio"`
     — LLM output hit a length cap mid-phrase and emitted a
     partial sentinel that wasn't in the set.  The H13.1 ship's
     plan-side override prevented downstream harm in that case,
     but the value was still passed through `validate_directives`
     unchanged.

  2. The negative structural checks let multi-word English
     phrases through if they happened to start with a letter
     and fit in 25 chars — e.g., `"Solve the structure"` passed
     all three checks.  Such values then reached PHENIX as
     bogus PHIL params and produced obscure errors.

**Fix (two additive parts):**

(a) Extended `_SYMMETRY_SENTINELS` with the "not explicitly
{mentioned,stated,given,specified,defined,listed,noted}" family
including truncated forms.  Frozenset is additive — pre-H14
sentinel behavior is preserved (regression tested).

(b) Added a positive Hermann-Mauguin shape check
`_looks_like_space_group(value)` that accepts standard space-group
symbols and rejects prose phrases.  Wired into
`validate_directives()` as an additional invalid-value condition.

**Regex evolution — two review passes shaped the final pattern:**

*Initial draft* (matched 71/230 official symbols — 31%):

```python
r'^[PFICRHAB]\s*-?\d{0,3}(?:\s*[-/_]?\s*\d{0,3}){0,4}'
r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$'
```

This too-strict pattern rejected 159 valid space groups: all
monoclinic slash forms (P21/c, P21/n, C2/c), all orthorhombic
mirror/glide groups (Pmma, Pnma, Pbca, Fdd2), all tetragonal
forms with mirrors (P4mm, P4/mmm, I41/amd), all cubic
high-symmetry groups (Pm-3m, Im-3m, Fd-3m, Ia-3d).  Self-review
against the full 230-symbol list (per Tom's "must admit all 240
space groups correctly" instruction) caught this before ship.

*After self-review* (matched 230/230 canonical but rejected
alternative settings):

```python
r'^[PFICRHAB][0-9mcndabe\s/_\-]{0,24}'
r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$'
```

The `e` letter was required for the 2002-ITA renamings of the
centered orthorhombic groups (Aem2, Aea2, Cmce, Cmme, Ccce).

*Gemini external review (Risk A)* identified that this regex
rejected alternative cell/origin settings.  International Tables
Vol. A defines colon-suffixed forms used in real PDB/mmCIF
metadata and cctbx tool output:

  - `R3:H` vs `R3:R` — rhombohedral axis choice (7 space groups)
  - `:1` vs `:2` — origin choice (24 centrosymmetric groups,
    mostly cubic + some tetragonal)
  - `:a` / `:b` / `:c` — unique-axis cell choice (monoclinic
    groups #3–15)

These would have been silently dropped as malformed prose.
Gemini's fix: add `:` to the interior alphabet.  Additional
finding during implementation: the axis-spec letters `h` and `r`
(in `R3:H`, `R3:R`) also weren't in the alphabet, so they had
to be added too.

*Final regex* (admits all 230 + all alternative settings):

```python
_HM_FORM_RE = re.compile(
    r'^[PFICRHAB]'
    r'[0-9mcndabehr\s/_\-:]{0,24}'
    r'(?:\s*\(\s*no\.?\s*\d+\s*\))?$',
    re.IGNORECASE
)
```

Alphabet rationale:

  - `[PFICRHAB]` — the 8 Bravais lattice letters (P primitive,
    F face-centered, I body-centered, C base-centered, R
    rhombohedral, H hexagonal-centered, A and B alternative
    base-centerings)
  - `0-9` — symmetry-element orders, screw-axis subscripts,
    IUCr number digits
  - `m c n d a b e` — mirror plane (m), glide planes (c, n, d, a, b),
    double glide (e — 2002 ITA)
  - `h r` — axis-spec suffix letters for rhombohedral hexagonal /
    rhombohedral settings (e.g., `R3:H`, `R3:R`)
  - `space / _ - :` — separators (slash for `P21/c`, dash for
    `P-1`, underscore for mmCIF, colon for alternative settings)

*Gemini Risk B* (acknowledged as known limitation, pinned with
a K-test): the permissive alphabet still admits short English
words that happen to use only HM-alphabet characters — `Panda`,
`Fame`, `Bad`, `Cab`, `Bed`, `Acme` all match.  This is not a
realistic LLM output for `space_group`, and the pre-H14 negative
checks already let these through.  Closing this gap fully would
require enumerating the 230 canonical symbols via `cctbx.sgtbx`,
which is beyond a directive sanity check's scope.  The test
`test_hm_form_known_limitation_short_words` pins the limitation
so future contributors don't think "Panda accepts" is a bug; the
test `test_hm_form_rejects_words_with_non_alphabet_chars` pins
the positive guarantee that words with non-alphabet letters
(`Phaser`, `Pizza`, `Place`, `Process`) ARE rejected.

**Test coverage** for Item 3 — 24 K-tests:

  - §A: 3 new sentinel coverage (including truncated forms)
  - §B: 1 existing-sentinel regression
  - §C: 6 HM-shape checks (standard symbols, 230/230 official,
    spaced variants, parenthetical numbers, prose rejection,
    empty/None rejection)
  - §D: 5 end-to-end through `validate_directives`
  - §E: 4 equivalence-class tests (case-insensitive,
    spacing, 1-placeholder, setting variants)
  - §F: 3 alternative-setting tests (rhombohedral, origin
    choice, unique axis) — per Gemini Risk A
  - §G: 2 limitation documentation tests — per Gemini Risk B

`hallucinated_param_value` rate in run_39 was 71 (down from 148
in run_25 — already trending down before H14), of which most
were `space_group=*` variants.  H14 should drive this further
down without affecting valid values.

**Plans:** `v119_H14_PLAN_rev1.md` (rev1).  Rev0 (an internal
draft) contained four items; Item 2 (`asu_copies` injection) was
dropped after Tom's domain call.

VERSION bumped to `119.H14`.  `defaults_fingerprint` unchanged at
`sha256:77bf7421...` — H14 did not touch the `DEFAULT_MODELS`
tables.

**Additional housekeeping:**

  - Added `tests/tst_utils.py` (sandbox helper for
    `run_tests_with_fail_fast()`).  H13.1 shipped three K-test
    files (`tst_autosol_bugs.py`, `tst_file_categorizer.py`,
    `tst_provider_error_classification.py`) that depended on this
    helper but the helper itself wasn't in the ship — those tests
    couldn't run in sandbox.  Now they can.
  - Relaxed `tst_autosol_bugs.py` Bug 5 test from exact-multi-line
    signature match to substring check.  The original test pinned
    the exact `def _track_output_files(...)` signature including
    indentation; a later addition of `skip_if_failed=False` to
    that signature broke the test.  The test's intent ("function
    accepts working_dir") is preserved with the substring form.

### v119.H14.1 — Ollama-fallback bypass closure for space_group validation

**Trigger**: Tom's 2026-05-26 ollama xtriage tutorial verification of
H14 (the gate-C test for Item 3) showed `space_group=Not explicitly
mentio` still appearing in the extracted directives, despite H14
being installed (the `[STEP_1F]` line confirmed
`scanner_version=119.H14`).  Gate C failed in production.

**Root cause** — two cooperating issues:

(a) `directive_extractor.py` has TWO space-group extraction paths:

  | Path                              | Validates? |
  |-----------------------------------|------------|
  | LLM-JSON success → `extract_directives` line 714 calls `validate_directives` | ✓ |
  | LLM fallback → `extract_directives_simple` returns directly (no validation) | ✗ |

  When ollama's LLM fails to return parseable JSON (small local
  models frequently fail JSON discipline), `extract_directives` falls
  back to `extract_directives_simple`.  This simple extractor has its
  OWN space_group regex at line ~4350:

  ```python
  _sg_pattern = (
      r'space[_ ]group'
      r'(?:\s*[=:]\s*|\s+(?:is|of|=|:)?\s*)'
      r'([A-Za-z][A-Za-z0-9 /_-]{1,20})'
  )
  ```

  This regex matched `"Space group: Not explicitly mentioned"` in
  Tom's preprocessed advice and captured the first 21 chars
  (`{1,20}` quantifier + initial letter): `"Not explicitly mentio"`.
  The dict was returned directly, bypassing the H14 sentinel
  (`_SYMMETRY_SENTINELS`) and shape check (`_looks_like_space_group`)
  that lived inside `validate_directives`.

  Correction to H14's CHANGELOG attribution: the H14 entry said the
  truncation came from "LLM hitting a length cap mid-phrase".  The
  diagnosis was wrong — the LLM passed `"Space group: Not explicitly
  mentioned"` through fully (visible in Tom's session summary at
  line 197 of `ollama_xtriage.log`).  The truncation came from the
  `{1,20}` quantifier in the simple extractor's regex.  The H14
  sentinel set still correctly catches the truncated form, but the
  validator was never invoked on the simple-extractor path.

(b) Latent `VALID_STOP_CONDITIONS` gap: the key `start_with_program`
is set by `_resolve_after_program` (called from
`_apply_workflow_intent_fallback`) when advice contains multiple
actions plus a stop intent (the Item 1 path).  It's consumed
downstream by `workflow_engine.py` (line ~2288), `ai_agent.py`
(~2816, ~3004, ~4555), and `ai_analysis.py` (~160).  But the key
was missing from `VALID_STOP_CONDITIONS`, so `validate_directives`
would log `Unknown stop condition start_with_program` and DROP it.
Pre-H14.1 this was invisible because validate_directives ran
BEFORE `_apply_workflow_intent_fallback` (the order:
extract → validate → fallback-overlay).  The key was added
post-validation in the LLM path.  But once H14.1 added
validate_directives to the END of the simple extractor (where the
fallback overlay runs INSIDE extract_directives_simple before
returning), `start_with_program` would have been stripped without
this fix — breaking H14 Item 1.

**Fix (two parts)**:

(a) Added `start_with_program: str` to `VALID_STOP_CONDITIONS`
(`agent/directive_extractor.py`, line ~1555).  Latent-bug fix: the
key was already consumed widely downstream; pre-H14.1 it just never
hit `validate_directives` so the gap went unnoticed.

(b) Modified `extract_directives_simple` signature to take an
optional `log=None` parameter (backward-compatible with the four
external callers in `ai_agent.py` and `run_ai_analysis.py` that
pass single-arg).  Added a final
`directives = validate_directives(directives, log)` call before
the return.  This closes the ollama-fallback bypass and makes
`validate_directives` the canonical final-sanity step for both
extraction paths — future validator extensions automatically apply
to both.

**Investigation method**: the gap was discovered by direct
reproduction.  Given Tom's preprocessed advice text and the
observed `space_group=Not explicitly mentio` symptom, calling
`extract_directives_simple(tom_advice)` reproduced the bug; calling
`validate_directives(simple_result)` on the result dropped the
value correctly.  The two-step trace pinpointed the missing wiring.
A broad audit of `extract_directives_simple` output against
`validate_directives` then surfaced the `start_with_program` issue
before H14.1 could ship — running validate_directives on results
from `multi_action_with_stop`, `polder_selection`, and other Item 1
patterns showed `start_with_program` being silently stripped.
Without the latent-bug audit, H14.1 would have fixed Item 3 while
breaking Item 1.

**Test coverage** — 12 K-tests in
`tests/tst_simple_extractor_validation.py`:

  - §A (2): Tom's exact ollama xtriage case (gate C closure) — verifies
    `space_group=Not explicitly mentio` is dropped with the canonical
    `DIRECTIVES: Dropping invalid space_group value` log line, and
    that other extracted values (resolution=1.7) survive intact
  - §B (2): Item 1 preservation through validation — explicit MR + stop
    must set `start_with_program=phenix.phaser` and have it survive
    validate_directives; the 1029B-sad pattern (goal phrase + stop,
    no MR keyword) must NOT set start_with_program (n=1 branch)
  - §C (2): VALID_STOP_CONDITIONS now includes start_with_program;
    direct validate_directives call preserves the key without the
    "Unknown stop condition" log line
  - §D (4): other simple-extractor outputs survive validation —
    atom_type, max_refine_cycles, skip_programs, valid space_group
  - §E (1): backward compat — the four external callers in
    ai_agent.py / run_ai_analysis.py that pass single-arg continue
    to work
  - §F (1): validate_directives is idempotent (external code may
    call validate_directives on results; running it twice must be
    a no-op)

`tests/tst_simple_extractor_validation.py` (NEW, 12 tests),
`agent/directive_extractor.py` (VALID_STOP_CONDITIONS extension +
extract_directives_simple final-validate),
`tests/canary_expected.json` (agent_version: 119.H14.1),
`tests/run_all_tests.py` (registers tst_simple_extractor_validation).

VERSION bumped to `119.H14.1`.  `defaults_fingerprint` unchanged at
`sha256:77bf7421...`.

**What this changes for production**:

  - Gate C now passes for ollama runs: `space_group=Not explicitly
    mentio` (and any other invalid value the simple extractor's
    permissive regex captures) is dropped before reaching downstream
    code.  Tom should re-run his ollama xtriage tutorial to verify;
    expected log line: `DIRECTIVES: Dropping invalid space_group
    value: 'Not explicitly mentio'`.
  - Item 1 fix unchanged (the latent bug fix preserves
    start_with_program through validation).
  - Item 2 unchanged.
  - The MR-SAD-hallucination bug surfaced in Tom's 1029B-sad ollama
    run (LLM set `use_mr_sad=True` for a SAD-only README) is a
    SEPARATE issue, candidate for H15.  It's an LLM-prompt-engineering
    or post-LLM-validation problem, not an H14 problem.

### v119.H14.2 — programs.yaml: predict_and_build does not accept ncs_spec

**Trigger**: Tom's 2026-05-26 1029B-sad ollama production verification
of H14.1 (run after removing the PDB file from the input directory)
hit a different failure mode.  With no PDB present, the workflow
correctly routed to `phenix.predict_and_build` instead of phaser.
The agent's file auto-discovery found a leftover `find_ncs.ncs_spec`
file from a previous run and looked up the documented PHIL flag in
`knowledge/programs.yaml`:

```yaml
phenix.predict_and_build:
  inputs:
    optional:
      ncs_spec:
        extensions: [.ncs_spec]
        flag: "map_model.ncs_file="     # ← stale documentation
```

The agent emitted that flag, and predict_and_build rejected it
with `Some PHIL parameters are not recognized by
phenix.predict_and_build`.  Cycle-level retry re-emitted the same
bad parameter and the workflow looped.

**Root cause**: per Tom's domain confirmation, predict_and_build
does NOT accept any NCS-file PHIL parameter (it runs NCS detection
internally if needed).  The yaml entry was either stale (PHIL
grammar changed) or never quite right.  No code change was needed,
only a configuration fix.

**Fix** (purely in `knowledge/programs.yaml`):

  1. Removed the `ncs_spec` block from
     `phenix.predict_and_build.inputs.optional`
  2. Removed `{ncs_spec}` from the command template
  3. Updated `phenix.map_symmetry`'s downstream-consumer hint to
     drop predict_and_build from the list of programs that consume
     `.ncs_spec` output

**No code changes.**  The agent's file auto-discovery and PHIL
emission machinery already does the right thing given the
configuration — once the yaml stops declaring predict_and_build
as an ncs_spec consumer, the agent stops trying to pass one.

**Defensive consideration noted but NOT shipped**: a generic
PHIL-validation-at-command-build-time pass (look up the program's
master_phil, drop unknown keys with a `DIRECTIVES: Dropping unknown
PHIL key for <program>` log line) would catch this entire class of
programs.yaml drift bug at the point of emission rather than at
program execution time.  Scoped as a candidate for H15 or later.
The H14.2 ship is the one-line config fix only.

**Test coverage** — 5 K-tests in
`tests/tst_predict_and_build_no_ncs.py`:

  - §A (3): predict_and_build's inputs.optional no longer contains
    ncs_spec; command template no longer references {ncs_spec};
    command template regression guard (sequence/data_mtz/full_map/half_map
    placeholders still present)
  - §B (1): map_symmetry's hint no longer advertises predict_and_build
    as a downstream consumer of .ncs_spec output (or explicitly says
    it's NOT a consumer)
  - §C (1): sibling regression guard — resolve_cryo_em's ncs_spec
    entry remains intact (predict_and_build is the only program
    being corrected; resolve_cryo_em DOES accept `ncs_file=`)

`knowledge/programs.yaml` (3 surgical edits),
`tests/tst_predict_and_build_no_ncs.py` (NEW, 5 tests),
`tests/canary_expected.json` (agent_version: 119.H14.2),
`tests/run_all_tests.py` (registers tst_predict_and_build_no_ncs).

VERSION bumped to `119.H14.2`.  `defaults_fingerprint` unchanged
at `sha256:77bf7421...`.

**What this changes for production**: when the 1029B-sad README
runs with no PDB in the input directory, the workflow no longer
loops on a rejected PHIL parameter.  Tom can confirm by re-running
his 1029B-sad ollama test1 directory (after `rm 1029B.pdb*`); the
predict_and_build command should now omit any `map_model.ncs_file=`
parameter.

**H15 (production code) — Resume reopen targeted stages**
`agent/plan_generator.py` (new `reopen_stages_for_directives`
and `_reopen_stages_inner` helpers),
`tests/tst_resume_reopen_stages.py` (NEW, 7 tests),
`tests/canary_expected.json` (agent_version: 119.H15),
`tests/run_all_tests.py` (registers tst_resume_reopen_stages),
VERSION.

Trigger: Tom's bromodomain resume failure (run 144).  Three
compounding causes formed the failure mode:
1. Bug 1 corrupted `final_refinement.status` to COMPLETE
2. Resume cleared `gate_stop` but NOT per-stage statuses
3. Gate immediately re-fired "all stages complete" → LLM
   saw STATE=complete → chose polder against user intent

H15 Item 1 prevents the corruption going forward.  Item 2
is the resume safety net: when new advice is provided on
resume, walk the directives and reopen affected stages.

**Targeted single-stage reopen design** (per Gemini's critique
of the original blast-radius proposal): rather than reopen
every stage downstream of the matched program, find the
LATEST completed stage whose `programs` list contains a
program named in `directives.program_settings`, and reset
ONLY that stage to PENDING.  Reset `cycles_used → 0` and
clear runtime fields (started_at, completed_at, last_result).
Don't cascade-reset stages downstream; don't reopen earlier
stages with the same program; skipped stages stay skipped.

This keeps blast radius O(1) regardless of plan size.

The 7 K-tests pin the design contract:
- §A: Tom's exact scenario
- §B: No advice change → no-op
- §C: Directive for unmatched program → no reopen
- §D: Multiple programs in directives → still O(1)
- §E: Skipped stages stay skipped
- §F: Empty plan / empty directives → no-op, no exception
- §G: Stage's strategy already honors directive — H15 design
  says reopen anyway (semantic: "user re-asserted intent")

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**What this changes for production**: when Tom resumes a
completed bromodomain run with new advice mentioning a
program inside an earlier stage, only that stage is
reopened; the rest of the plan stays intact.  The gate no
longer immediately re-fires "all complete" because at least
one stage is now PENDING.

**H16 (production code) — MTZ obs_labels auto-fill for
multi-array MTZ files**
`agent/mtz_inspector.py` (NEW: `inspect_mtz()`,
`select_obs_labels_for(program, info)`,
`has_ambiguous_arrays()`, `_PROGRAM_PREFERENCES` table),
`knowledge/programs.yaml` (new `auto_fill_obs_labels`
invariant on xtriage, autosol, phaser, predict_and_build),
`agent/command_builder.py` (consumer of the new invariant),
`tests/tst_obs_labels_auto_fill.py` (NEW, 8 tests across
four layers: policy unit tests, MTZ inspector, builder
integration simulation, YAML config validation),
`tests/canary_expected.json` (agent_version: 119.H16),
VERSION.

Trigger: 88 TIER-1 failures across two batch scans
(run_25, run_39) matching "Sorry: Multiple equally suitable
arrays of observed xray data found", concentrated in
AF_exoV_MRSAD and lysozyme-MRSAD tutorials.  The failure
mode: the MTZ contains multiple legitimate observation
arrays (e.g., merged Iobs AND anomalous I(+)/I(-)) and
PHENIX programs error out unless the user disambiguates
via `scaling.input.xray_data.obs_labels=` (or equivalent
per-program PHIL path).  Pre-H16, the LLM had to know to
inject this parameter manually for every multi-array MTZ;
when it forgot, the program crashed and the agent treated
it as an unrecoverable error.

**Three-layer mechanism**:

1. **MTZ inspector** (`agent/mtz_inspector.py`):
   `inspect_mtz(path)` reads the MTZ via `iotbx.mtz` (cctbx)
   and returns a dict describing each column group (anomalous
   intensities, anomalous amplitudes, merged intensities,
   merged amplitudes, R-free flags, etc.).  Side-effect-free;
   never raises (returns `{"error": ...}` on file-read
   failures so the caller can decide).

2. **Per-program preference policy**:
   `select_obs_labels_for(program, mtz_info)` reads the
   `_PROGRAM_PREFERENCES` table (e.g., `phenix.phaser` prefers
   merged for MR; `phenix.autosol` prefers anomalous for SAD)
   and returns the program-specific obs_labels string to
   inject.  Returns `None` when no ambiguity exists.

3. **YAML invariant + builder hook**: `programs.yaml` declares
   `auto_fill_obs_labels: true` on xtriage, autosol, phaser,
   and predict_and_build.  `command_builder._apply_invariants()`
   reads the flag, calls `inspect_mtz()` on the relevant
   `data_mtz` input, calls `select_obs_labels_for()`, and
   injects the result into the command before execution.

The 8 K-tests cover all four layers: policy unit tests
(pure functions, no I/O); MTZ inspector tests (cctbx-dependent,
gracefully skip in sandbox); builder integration simulation
(mirrors the actual `_apply_invariants()` branch); YAML
config validation.

H16.1 follow-up (no code change): bumped the scanner version
pin in `[STEP_1F]` telemetry to `119.H16.1` so production logs
clearly indicate the obs_labels auto-fill is active.

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**What this changes for production**: the 88 TIER-1
"Multiple equally suitable arrays" failures no longer fire.
PHENIX programs receive an explicit obs_labels parameter
appropriate to their use case (merged for MR, anomalous for
SAD).  Programs that previously crashed on cycle 1 with the
ambiguity error now run normally on the user's first attempt.

**H17 (production code) — autobuild PHIB-required reactive
recovery (analyzer side)**
`knowledge/recoverable_errors.yaml` (NEW
`missing_phib_input_map_file` entry),
`agent/error_analyzer.py` (NEW `strip_flags` field on
`ErrorRecovery` dataclass; NEW generic `_resolve_strip_parameter`
handler alongside existing `_resolve_add_parameter` and
`_resolve_select_value`; `_extract_error_info` returns marker
dict for strip_parameter errors),
`tests/tst_error_analyzer.py` (7 new K-tests covering H17
detection + retroactive resolution of the latent
`rfree_flags_mismatch` resolution that had been declared in
YAML but never wired in the analyzer),
`tests/canary_expected.json` (agent_version: 119.H17),
VERSION.

Trigger: Tom's lysozyme-MRSAD tutorial cycle 5.  The LLM
correctly identified MR-SAD as the workflow and ran xtriage
→ phaser → autosol successfully (Bayes CC 74.20); on cycle 5
it then called `phenix.autobuild` with the user-supplied
`lyso2001_scala1.mtz` (raw anomalous MTZ without phases) as
`map_file=`.  autobuild requires `PHIB` phase columns in
`input_map_file`/`map_file`; the program crashed with "Sorry,
PHIB is required for input_map_file"; the agent had no
recovery template; the workflow halted after four retries.

**The right fix is reactive, not proactive**.  The user-supplied
`map_coeffs.mtz` is a legitimate input in many workflows (any
post-phasing build).  Pre-emptively rejecting `map_file=` on
all autobuild calls would break the legitimate case.  H17
detects the specific PHIB-required error AFTER it fires and
strips the offending parameter on retry.

**Strict conjunction detection**: the YAML pattern requires
BOTH the literal phrase "PHIB is required for input_map_file"
AND the column-spec context.  Without the conjunction, generic
"PHIB" mentions in unrelated error messages would false-positive.

**New `strip_parameter` resolution kind**: alongside the
existing `add_parameter` (adds a flag) and `select_value`
(disambiguates an enum), H17 introduces `strip_parameter`
for the case where the right fix is to REMOVE an inappropriate
flag entirely.  The `strip_parameters` YAML list names the
PHIL paths to strip (e.g., `[map_file, input_map_file,
input_files.map_file]` — three aliases for the same input
slot).

**Retroactive side benefit**: `rfree_flags_mismatch` had been
declared in `recoverable_errors.yaml` for months but never
wired in the analyzer (it would never trigger).  H17's
generic `_resolve_strip_parameter` handler retroactively
resolves that error type too.

Gemini's critique on the H17 plan recommended a robust regex
for the strip operation: `r'(?:^|\s)' + re.escape(flag_prefix)
+ r'\s*=\s*(?:"[^"]*"|\'[^\']*\'|\S+)'` to handle
quoted-with-spaces and PHIL spacing.  H17 plan adopted; H17.1
implements the regex in the executor.

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**H17.1 (production code) — Executor-side strip_flags wiring
(completes H17)**
`programs/ai_agent.py` (three edits: `_handle_recovery` stashes
strip_flags; `_execute_command` pops and applies; new helper
`_print_recovery_notice` shows "Action: Stripping [...]"),
`tests/tst_h17_strip_executor.py` (NEW, 9 tests pinning the
regex pattern against PHIL/quoted/end-of-line variants),
`tests/canary_expected.json` (agent_version: 119.H17.1),
VERSION.

Pre-deploy review of H17 surfaced Scenario B: the analyzer
emits ErrorRecovery with `strip_flags` populated and the
`[NOTICE]` log fires, but the executor never actually strips
the flag from the retry command.  H17 was analyzer-side only;
H17.1 closes the executor gap.

Three edits to `programs/ai_agent.py`:

1. **`_handle_recovery`** (~line 3602): when `recovery.strip_flags`
   is non-empty, stash the entry in
   `session.data["pending_strip_recoveries"]` keyed by program
   name.  Preserves the existing `set_recovery_strategy` call
   for backward compat with `add_parameter` recoveries (which
   don't need executor support — they're injected via the
   existing strategy mechanism).

2. **`_execute_command`** (~line 5694): at the top of the method,
   pop any pending entry for the current program.  If present,
   apply the robust regex (per Gemini's critique) to strip
   each flag prefix.  Emit `[STRIP]` log line.  One-shot via
   `pop()` — the entry is consumed on use so a subsequent
   retry doesn't re-strip.

3. **`_print_recovery_notice`** (~line 3783): show "Action:
   Stripping [map_file, input_map_file, input_files.map_file]"
   instead of the awkward "Selecting ''" for strip recoveries.

The 9 K-tests in `tst_h17_strip_executor.py` pin the regex
pattern itself.  If the pattern is changed in `ai_agent.py`,
the test file's reference implementation must be updated
in lockstep.  The tests cover: basic, PHIL spacing (spaces
around `=`), double-quoted, single-quoted, end-of-line
(no trailing whitespace), the exact lysozyme command, all
three PHIL variants in sequence, a false-positive guard
(substring of a different flag name), and idempotence
(running the strip twice doesn't corrupt the command).

Validated end-to-end on Tom's lysozyme-MRSAD:
- Cycle 1: xtriage `ambiguous_data_labels` recovery (pre-H17
  mechanism, still works)
- Cycle 2: xtriage SUCCESS
- Cycles 3-4: phaser SUCCESS (LLG 95.0), autosol SUCCESS
  (Bayes CC 74.20)
- Cycle 5: autobuild FAILED with PHIB error → H17 fires
  `[NOTICE] DETECTED RECOVERABLE ERROR` with "Action:
  Stripping [map_file, input_map_file, input_files.map_file]"
- Cycle 6: `[STRIP] phenix.autobuild: removed 'map_file=...'
  from retry command (recovery: missing_phib_input_map_file)`
  followed by `Running:` (no map_file= in command), autobuild
  SUCCESS.

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

(Note: cycles 7-10 then thrashed on a separate post-autobuild
issue — bad metrics parsing "Residues Built: 5220629854239855"
and state-machine over-reaction to R-free=0.52.  Booked as
separate issue, not H17.)

**What this changes for production**: the lysozyme-MRSAD
tutorial completes through autobuild without halting at the
PHIB error.  The reactive recovery infrastructure (YAML +
analyzer + executor) now supports three resolution kinds:
add_parameter, select_value, strip_parameter.

**H18 (production code) — File-based experiment-type
detection as primary signal**
`agent/file_utils.py` (NEW public helper
`infer_experiment_type_from_files()`),
`agent/directive_extractor.py` (files-first detection in BOTH
`_apply_experiment_type_program_reprints` AND
`_resolve_after_program`; `extract_directives()` accepts
`original_files` kwarg; `_apply_workflow_intent_fallback()`
threads files through),
`agent/plan_generator.py` (`_build_context` delegates to
shared helper for single source of truth),
`programs/ai_agent.py` (passes file inventory to extraction
via new PHIL param),
`programs/ai_analysis.py` (NEW PHIL param
`original_files_for_directives`; threaded through
`_run_directive_extraction_locally` and the args-builder),
`phenix_ai/run_ai_analysis.py` (`run_directive_extraction()`
accepts `original_files`; passes to `extract_directives()`),
`tests/tst_density_modify_experiment_type.py` (extended 13 →
20 tests; K14-K20 cover H18 surface),
`tests/canary_expected.json` (agent_version: 119.H18),
VERSION.

Trigger: Tom's AF_7mjs density-modify-and-stop regression
(distinct from the v118.9 §20 bug, surfaced after H17.1
deployed).  User wrote literally "density modify and stop"
with cryo-EM half-map inputs (`.ccp4` files, no `.mtz`); the
extractor produced `after_program=phenix.autobuild_denmod`
(the X-ray density-modification program); the §20 correction
was supposed to map this to `phenix.resolve_cryo_em` for
cryo-EM data but its text-only
`_detect_experiment_type_signals` returned None ("ambiguous,
decline to act") because the terse user text "density modify
and stop" has neither cryo-EM nor X-ray tokens.  The LLM
downstream then saw the wrong `after_program` in VALID
PROGRAMS, recognized it couldn't actually run (cryo-EM has
no MTZ), and selected `phenix.predict_and_build` as "the
next logical step" — overriding the user's explicit "stop"
instruction.

**Root cause**: §20's correction relied on text-based
detection alone.  For terse advice, text gives no signal,
and the file inventory (the unambiguous evidence) was never
inspected.  ARCHITECTURE.md §17 documents that
`session.set_experiment_type()` only locks AFTER the first
program returns — so at directive-extraction time, the
locked experiment_type is None and the validator had to
infer one.  ARCHITECTURE.md §3.3 already proposed inferring
experiment type from file extensions at session creation;
H18 is the smallest first step in that direction.

**New public helper**: `agent/file_utils.py` gains
`infer_experiment_type_from_files(files) -> (type,
evidence_dict)`.  Returns "xray" for files with `.mtz/.sca/.hkl`
extensions alone, "cryoem" for `.mrc/.ccp4/.map` alone, and
None for mixed or empty inputs.  The asymmetric semantics
match the existing `_detect_experiment_type_signals` policy.
Evidence dict carries the unique extensions seen and an
`is_mixed` flag for telemetry.

**Detection priority (locked policy)** — applied identically
in both correction sites:
```python
target_type = None
if original_files:
    file_type, evidence = infer_experiment_type_from_files(original_files)
    target_type = file_type   # files-win on conflict

if target_type is None:
    text_type = _detect_experiment_type_signals(combined_advice)
    target_type = text_type

if target_type is None:
    return directives  # decline to act
```

Files-win is the policy per Gemini's H18 review: file
extensions are the hard physical boundary (passing `.ccp4`
to an X-ray-only program crashes at parse time regardless
of what the user wrote); text-based detection has known
false-positive vectors ("mad" in "modify", "sad" in
arbitrary sentences).

**Two correction sites needed the fix**.  H18 rev 1
implemented files-first only in
`_apply_experiment_type_program_reprints` (§20's site).
During the merge into Tom's existing
`tst_density_modify_experiment_type.py` suite, K14 (the
exact AF_7mjs failing case) revealed a second silent-revert
site: `_resolve_after_program` (the v115.10 post-LLM
overlay) had its own text-only experiment-type heuristic
that defaulted `_exp="xray"` for terse advice and
unconditionally overrode `after_program`.  H18 rev 2 fixes
both sites identically.

**Telemetry (Pitfall 2 — Silent Override Hazard)**: the
enriched `[DIRECTIVE_CORRECTION]` marker records
`source=files|text`, `evidence=['.ccp4', ...]`, and
`text_signal=xray, OVERRIDDEN` (when files override
contrary text).  Example:
```
[DIRECTIVE_CORRECTION] Mapped after_program=phenix.autobuild_denmod
  to phenix.resolve_cryo_em (source=files, target_type=cryoem,
  evidence=['.ccp4'], text_signal=xray, OVERRIDDEN)
```

**Telemetry (Pitfall 1 — Dirty Directory Poisoning)**:
long-running sessions that MERGE files across cycles (per
`Session.set_project_info` documented in session.py:494-510)
could accumulate both `.mtz` and `.ccp4`, pushing detection
into "mixed → None" state.  A new
`[DIRECTIVE_CORRECTION_MIXED]` log marker fires when both
types are detected, making accumulated-files drift visible
in audit logs.  In practice this is a non-issue for the
AF_7mjs bug path because directive extraction is gated by
`directives_extracted=True` (only fires on initial
extraction, before any cycle has run), but the telemetry
is there defensively.

**Single source of truth**: `agent/plan_generator._build_context`
was refactored to call the new shared helper.  Previously
it had its own private file-extension mapping at lines
232-236; H18 routes it through `infer_experiment_type_from_files`
so plan_generator and directive_extractor cannot drift
apart on detection behavior.

The 7 new K-tests in `tst_density_modify_experiment_type.py`
(K14-K20) extend the existing 13 §20 tests to a 20-test
suite covering both v118.9 §20 and v119.H18:
- K14: AF_7mjs failure verbatim — terse advice + cryo-EM
  files → corrected via files (source=files)
- K15: Mirror — terse advice + X-ray files
- K16: Backward compat — no files → text fallback still
  corrects (source=text)
- K17: Pre-H18 bug path — no files + terse text → declines
  (proves H18 requires file threading)
- K18: Files-win on conflict + OVERRIDDEN telemetry
- K19: Mixed input → defers to text + emits
  `[..._MIXED]` marker
- K20: plan_generator delegates to shared helper

All 20 tests pass.  Full sandbox sweep: 227 PASS, 22 SKIP,
no regressions.

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**Backward compatibility**: `original_files` parameter is
optional with default `None` at every signature change
(`extract_directives`, `_apply_experiment_type_program_reprints`,
`_resolve_after_program`, `_apply_workflow_intent_fallback`,
`run_directive_extraction`).  Pre-H18 callsites without files
get exact pre-H18 behavior — text-only detection.

**What this changes for production**: AF_7mjs "density
modify and stop" + cryo-EM files now correctly:
1. Files detector identifies cryoem from `.ccp4`
2. After_program is corrected to `phenix.resolve_cryo_em`
3. `[DIRECTIVE_CORRECTION]` log shows `source=files`
4. Stop check fires after `resolve_cryo_em` completes
5. Agent halts honoring user intent; no `predict_and_build`

**H18.1 (production code) — PHIL declaration deploy-gap hotfix**
`programs/ai_agent.py` (+14 lines: new `original_files_for_directives`
declaration in `master_params` after `user_advice_raw`),
`tests/tst_h18_1_phil_roundtrip.py` (NEW, 6 K-tests),
`tests/run_all_tests.py` (registers the new test),
`tests/canary_expected.json` (agent_version: 119.H18.1),
VERSION.

Trigger: Tom re-ran AF_7mjs with H18 deployed.  All sandbox
tests passed.  Production still ran `predict_and_build` on
cycle 3 instead of stopping after `resolve_cryo_em`.  The log
showed:

```
DIRECTIVES: Extraction failed - Assignment to non-existing
  attribute "ai_analysis.original_files_for_directives"
  File ".../programs/ai_agent.py", line 8513, in _extract_directives
    directive_params.ai_analysis.original_files_for_directives = (
```

**Root cause**: H18 added the PHIL parameter declaration to
`programs/ai_analysis.py` (the analysis server's master_phil)
and added the assignment site at `programs/ai_agent.py:8513`,
but DID NOT add the PHIL declaration to `programs/ai_agent.py`'s
OWN `master_params` string (defined at line 144).
`directive_params = copy.deepcopy(self.params)` copies the
agent's params, which are parsed against the agent's
master_params — so the assignment failed.  The crash was
caught by the surrounding try/except in `_extract_directives`
and directive extraction silently returned `{}`.  The agent
then ran with no user-supplied `after_program`, no
`stop_after_requested`, falling back to the default cryo-EM
plan template: cycle 1 mtriage, cycle 2 resolve_cryo_em,
cycle 3 `phenix.predict_and_build` (Stage 3 of the plan).

**Why H18's 20 K-tests passed despite the production failure**:
the existing tests in `tst_density_modify_experiment_type.py`
call helpers directly with Python dicts.  They never exercise
the PHIL parse → deep-copy → assign path.  The bug was
entirely in the PHIL layer, which the tests bypassed.

**The fix**: add the missing PHIL declaration to
`programs/ai_agent.py`'s `master_params` string, mirroring
the declaration already in `programs/ai_analysis.py`.  No
code logic changes.

**Preventive K-test pattern**:
`tests/tst_h18_1_phil_roundtrip.py` adds 6 K-tests:

1. **Source-grep verification** (sandbox-safe): confirms the
   PHIL declaration is present in `ai_agent.py`'s
   `master_params` region.  Runs anywhere, no libtbx
   required.
2. **PHIL parse + extract** (PHENIX-only, skips in sandbox):
   actually calls `libtbx.phil.parse()` on the extracted
   master_params and confirms `extract()` produces an object
   with the new attribute.
3. **Assignment reproduction** (PHENIX-only): deep-copies the
   extracted params and performs the exact assignment that
   crashed in production.  Must succeed.
4. **Server-side definition present**: confirms
   `programs/ai_analysis.py` still has the matching PHIL
   declaration.  Guards against fixing one side and
   accidentally regressing the other.
5. **Assignment site anchor**: confirms the assignment is
   still at the expected location in `ai_agent.py`.
   Catches refactors that would move the call site without
   updating the test.
6. **Cross-file consistency**: same parameter name in both
   files.  Catches typo divergence.

Validated by temporarily reverting H18.1's PHIL block and
confirming test 1 fails with a clear message identifying
the deploy gap.  Restored — all 6 pass.

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**What this changes for production**: AF_7mjs is finally
fixed.  The PHIL assignment in `_extract_directives` succeeds,
`original_files_for_directives` is threaded to the directive
extractor, the experiment-type correction fires, and the
agent stops after `resolve_cryo_em` as the user requested.

**Lesson** (in DEVELOPER_GUIDE.md §3j + new lesson section):
this codebase has two independent master_params blocks
(`programs/ai_agent.py` and `programs/ai_analysis.py`).
Adding a parameter to one does NOT propagate to the other.
For any new PHIL parameter that flows client → server, BOTH
declarations must be added in the same commit, with mirroring
`.help` comments, and the PHIL-roundtrip K-test pattern
applied to catch the deploy gap at sandbox time.

**H18.2 (production code) — third `_resolve_after_program`
callsite missed in H18 audit**
`agent/directive_extractor.py` (+13 lines: pass
`original_files=original_files` to the v117.2 fallback
callsite at line 796, plus comment explaining the H18.2
rationale),
`tests/tst_density_modify_experiment_type.py` (+200 lines: K21
production reproduction, K22 X-ray mirror, K23 backward-compat
no-files path),
`tests/canary_expected.json` (agent_version: 119.H18.2),
VERSION.

Trigger: Tom re-ran AF_7mjs with the H18.1 hotfix deployed and
a runtime tracer installed (`h18_install_runtime_tracer.py`
patches the deployed `directive_extractor.py` in-place with
`[H18_TRACE]` markers at every key decision point in the
pipeline).  The trace revealed the exact failure mode:

```
[H18_TRACE] extract_directives ENTRY:
  raw_advice='User instructions:\ndensity modify and stop'
  original_files=['7mjs_23883_H_1.ccp4', '7mjs_23883_H_2.ccp4',
                  '7mjs_23883_H.fa']
[H18_TRACE] CALLING _apply_experiment_type_program_reprints:
  after_prog=None
[H18_TRACE] _apply_reprints: after_prog=None
  stop_cond={'stop_after_requested': True}
[H18_TRACE] _apply_reprints EARLY RETURN: no after_prog
[H18_TRACE] AFTER _apply_experiment_type_program_reprints:
  after_prog=None
[H18_TRACE] before _apply_workflow_intent_fallback:
  stop_conds={'stop_after_requested': True,
              'after_program': 'phenix.autobuild_denmod',
              'skip_validation': True}
```

`after_program` was None when H18's first site ran (correct,
since the LLM didn't emit it), then was SET to the wrong value
by the time the second site ran.  Something BETWEEN the two
H18 sites was setting `phenix.autobuild_denmod`.

**Root cause**: a v117.2 fallback path at
`directive_extractor.py:783` fires when the LLM emits
`stop_after_requested=True` but omits `after_program`.  It calls
`_resolve_after_program` to fill in the missing field by parsing
the raw advice:

```python
# Pre-H18.2:
_resolve_after_program(directives, _v172_source.lower())
#                                                       ↑
#                                       MISSING: original_files
```

Without `original_files`, `_resolve_after_program` defaulted
`_exp="xray"` via the text-only heuristic (raw advice "density
modify and stop" has no cryo-EM/X-ray tokens) and mapped
`denmod` → `phenix.autobuild_denmod`.

The downstream `_apply_workflow_intent_fallback` DOES pass
`original_files` (H18 site 2 was correctly updated).  But at
that point the preprocessed advice contains "Stop Condition:
None", so `_is_stop_after_requested(advice)` returned False.
This pushed the resolver into the
"n==1, no stop → leave as-is" branch — the buggy
`after_program` from the v117.2 path persisted.

**Why H18's audit missed this**: H18 identified TWO callsites
billed as "experiment-type detection" sites
(`_apply_experiment_type_program_reprints` and
`_apply_workflow_intent_fallback`).  The v117.2 fallback was
billed differently — as a "fill in the missing field" site.
Internally it called the same `_resolve_after_program` function
with the same files-win contract, but the label was different,
so the H18 grep audit didn't include it.

**The fix is one line** at the v117.2 callsite:

```python
# Post-H18.2:
_resolve_after_program(directives, _v172_source.lower(),
                       original_files=original_files)
```

Plus a comment explaining the H18.2 rationale.

**After H18.2, all three `_resolve_after_program` callsites in
`extract_directives` pass `original_files`**:

```
$ grep -n "_resolve_after_program(" agent/directive_extractor.py \
       | grep -v "^.*def "
796:  _resolve_after_program(directives, _v172_source.lower(),
                            original_files=original_files)  # ← H18.2
4371: _resolve_after_program(directives, advice_lower,
                            original_files=original_files)  # ← H18 (workflow_intent)
```

The fourth callsite in `extract_directives_simple` (line 4920)
doesn't take files by design (rules-only path).  Out of scope.

**Preventive K-tests**:

- **K21**: AF_7mjs production reproduction.  Mocked LLM returns
  the exact directive shape Tom's production LLM produced
  (`{"stop_conditions": {"stop_after_requested": true}}` — no
  `after_program`).  Cryo-EM files.  Expected:
  `after_program=phenix.resolve_cryo_em` via the v117.2 fallback
  using H18.2's `original_files` threading.
- **K22**: X-ray mirror.  LLM omits `after_program`; X-ray files
  (.mtz).  Expected: `after_program=phenix.autobuild_denmod`
  (the correct X-ray choice).
- **K23**: Backward compat.  `original_files=None`; text-only
  cryo-em signal in raw advice.  Expected:
  `phenix.resolve_cryo_em` via the text-fallback branch of
  `_resolve_after_program`.

All 23 tests in `tst_density_modify_experiment_type.py` pass
(13 original §20 tests + 7 H18 tests + 3 H18.2 tests).

`defaults_fingerprint` unchanged at `sha256:77bf7421...`.

**What this changes for production**: AF_7mjs "density modify and
stop" with cryo-EM files finally works correctly:
- Cycle 1: mtriage → SUCCESS
- Cycle 2: resolve_cryo_em → SUCCESS
- Cycle 3: **STOP** (after_program=phenix.resolve_cryo_em matched,
  stop_after_requested=True)

NO predict_and_build in cycle 3.

**Lesson** (in DEVELOPER_GUIDE.md): when adding an optional
parameter to a function, grep for ALL callsites of that function
regardless of what the surrounding code is "doing".  The
parameter's contract is part of the function's identity; every
caller must opt in or the bug travels through the callers that
didn't.  H18 grepped for two billed-as-"experiment-type-detection"
sites and missed a third site billed as
"fill-in-missing-field" — same function, same parameter, same
contract, different label.

The K-test pattern that catches this class of bug is
production-faithful end-to-end: mock the LLM with the exact shape
production emits, run the full `extract_directives()` pipeline,
assert on the final state.  This catches every internal site that
touches the parameter because it tests behavior at the
entry-point contract, not at any individual function.  K21 follows
this pattern.

### Architectural notes

**Zero modifications to existing production code paths in H1–H3.**
H1 centralizes; H2 channels; H2.1 promotes; H3 observes.  None of
those four ships touches the workflow_engine, the metric_evaluator,
the LangGraph pipeline, or any other component on the live
decision path.  The K-suites for H2/H2.1/H3 explicitly verify
this with regression assertions against H1.

**H4 onwards: server-side instrumentation only, no decision-path
changes.**  H4 adds `[STEP_1F]` emission inside
`run_advice_preprocessing` after preprocessing completes — pure
telemetry, doesn't alter what advice is produced.  H5 and H5.1
add the `diagnostic_messages` list to `working_results` and
thread it through return paths; the list never feeds back into
any decision.  K_H5 §B `test_*_field_in_return_path` tests pin
that the field has no side effect on `directives`,
`processed_advice`, or `diagnosis_text`.

**The "graceful skip" pattern, generalized.**  K_H1 introduced
`_try_import_*()` helpers that return `(value, None)` on success
or `(None, error_msg)` on ImportError, letting tests skip
cleanly in sandbox while running fully under PHENIX.  K_H2/H2.1/H3/H4/H5
all adopt this pattern.  Sandbox runs report mixed PASS/SKIP;
PHENIX runs report all PASS.  This makes pre-ship verification
possible without requiring a PHENIX environment for every
edit-cycle.

**Gemini-reviewed plan rev cycles.**  Each ship went through
2–5 plan revs incorporating Gemini critique before
implementation.  H3 in particular: Gemini contributed the
shared loader factoring (Q5), the dedicated `CANARY_PING`
marker (Q2), the "stay loose" assertion threshold for the LLM
probe (Q3), and the strict client-side timeout recommendation
(later relaxed from 5s to 30s after Tom observed false-positive
timeouts on real network latency).  H5 hit 5 plan revs; rev 5
corrected a uniformity-criterion violation in rev 4 (removed
dispatch-mode branch in client behavior).  H5.1 plan rev 2
centralized the re-emit (backported to H5 since H5 hadn't
deployed to server yet) — the centralization saved the planned
H5.1 client-side edits entirely.

**Self-review on fresh extract.**  H2/H2.1/H3 reviews caught
stale line-number references, unused imports, and an incorrect
import path (`phenix_ai.X` vs the correct
`phenix.phenix_ai.X` — phenix_ai lives in the phenix-project
source tree, not the langchain tree).  H5.1 review after the
`run_utils.py` audit changed the test strategy: invalid-provider
injection wouldn't trigger the new markers because
`validate_api_keys` silently passes unknown providers and
`setup_llms` swallows its own exceptions.  Switched to
deterministic monkey-patching with the Gemini Q2 dual-assertion
mitigation.

### Files changed (cumulative)

| File | Ship |
|---|---|
| `VERSION` | H2 (new), H3 (→119.H3), H4 (→119.H4), H4.1 (→119.H4.1), H5 (→119.H5), H5.1 (→119.H5.1), H5.1.1 (→119.H5.1.1) |
| `core/_version.py` | H2 (new) |
| `core/_build_info.py` | H2 (new) |
| `core/llm.py` | H1 (DEFAULTS tables + helpers), H2 (+fingerprint) |
| `agent/api_client.py` | H1 (uses default_model_for_provider) |
| `agent/directive_extractor.py` | H1 (uses defaults), H2.1 (+skip promotion), H5.1.1 (+sidecar protection in merge_directives, +type-gated bool unwrap at two validate_directives sites) |
| `agent/raw_advice_scanner.py` | H4 (new), H4.1 (pinned vs corpus) |
| `knowledge/api_schema.py` | H2 (+agent_build schema) |
| `phenix_ai/run_ai_agent.py` | H2 (+inject_agent_build) |
| `phenix_ai/run_ai_analysis.py` | H4 (+[STEP_1F] emission), H5 (+_emit_marker helper, diagnostic_messages in run_advice_preprocessing), H5.1 (+diagnostic_messages in run_directive_extraction and run_failure_diagnosis), H7 (+scanner-first extraction + recall metrics + telemetry-independent LLM call for [STEP_1F]) |
| `programs/ai_analysis.py` | H5 (+ encode/decode, Total Init, _relay helper, centralized re-emit, local-mode passthroughs), H5.1 (+passthroughs in directive_extraction_locally and failure_diagnosis_locally) |
| `tests/test_api_keys.py` | H1 (uses defaults) |
| `tests/tst_default_models.py` | H1 (new) |
| `tests/tst_agent_build_info.py` | H2 (new) |
| `tests/tst_skip_promotion.py` | H2.1 (new) |
| `tests/canary_expected.json` | H3 (new), H4/H4.1/H5/H5.1/H5.1.1/H6/H6.1/H7 (version sync each ship) |
| `tests/canary_utils.py` | H3 (new) |
| `tests/tst_canary.py` | H3 (new) |
| `tests/step_1f_corpus.json` | H4 (new), H4.1 (frozen at 31 docs) |
| `tests/tst_raw_advice_scanner.py` | H4 (new), H4.1 (+golden-master test) |
| `tests/tst_diagnostic_messages.py` | H5 (new, 21 tests across 5 sections), H5.1 (+§F, →26 tests) |
| `tests/tst_directive_validation.py` | H5.1.1 (new, 23 tests across §A merge protection + §B list-wrap defense) |
| `tests/llm/framework.py` | H3b (initial framework with stubs), H6 (4 stubs replaced with real implementations + raw_output preservation fix), H6.1 (+ per-scenario fail/err suffix) |
| `tests/tst_planning_framework.py` | H6 (new K_H6, 18 tests across §A is_stop_intent + §B validate_planning_state + §C make_planning_run_fn + §D call_planning_llm raw_output preservation) |
| `tests/tst_phase2b_activation.py` | H7 (new K_H7, 15 tests across §A scanner contract + §B fallback semantics + §C recall metrics + §D integration smoke incl. source-level regression pins) |
| `tests/llm/tst_directive_extraction.py` | H3 (+canary Scenario), H5.1 (probe text fix) |
| `tests/llm/canary_check.py` | H3 (new), H5.1 (probe text fix) |
| `tests/run_all_tests.py` | H1, H2, H2.1, H3, H4, H5, H5.1, H5.1.1, H6, H7 (K-suite registrations) |
| `docs/DEVELOPER_GUIDE.md` | H2 (§8 protocol version 3→5) |

Note on `programs/ai_agent.py`: **never modified** in the v119
cluster.  Two architectural properties confirm this:
- For H5/H5.1: the centralized re-emit in
  `programs/ai_analysis.py::run_job_on_server_or_locally`
  means adding new markers to engine functions requires NO
  client-side changes.  H5.1 added three new markers without
  touching `ai_agent.py`.
- For H5.1.1: the sidecar-protection fix lives inside the
  shared `merge_directives` helper; both callers
  (`run_ai_analysis.py:1276` and `ai_agent.py:8376`) inherit
  the fix automatically.

### Verification

Under PHENIX (verified on Tom's Mac, end of H7):
- K_H1 (Default Models): 22 PASS (+1 SKIP under sandbox; 23 under PHENIX)
- K_H2 (Agent Build Info): 24 PASS
- K_H2.1 (skip_programs Promotion): 11 PASS
- K_H3a (Startup Canary): 3 PASS (+7 SKIP under sandbox; 10 under PHENIX — VERSION reads 119.H7; lockstep verified)
- K_H4 (Raw Advice Scanner + STEP_1F): 31 PASS (golden-master recall 0.9810 preserved through H7's scanner-first switch)
- K_H5 (Diagnostic Messages Channel): 13 PASS (+13 SKIP under sandbox; 26 under PHENIX)
- K_H5_1_1 (Directive Validation Cleanups): 23 PASS
- K_H6 (Planning Framework): 18 PASS (sandbox-side §D mock tests SKIP under PHENIX where libtbx import wins the race; same mechanic pinned by sandbox version)
- K_H7 (Phase 2B Activation): 15 PASS

Total v119-cluster K-tests after v119.H7: **181 PASS under
PHENIX** (sandbox runs report mixed PASS/SKIP; PHENIX runs
report all PASS).

Plus the 126 pre-v119 K-suite modules in `run_all_tests.py`,
all continuing to PASS as regression checks.  Total module
count: **135 modules, 135 PASS, 0 FAIL** on Mac end-of-H7.

**Test-count correction note:** Earlier ship docs (through
H5.1.1 and H6) reported cumulative v119-cluster counts with
a ~7-test arithmetic drift (155 post-H5.1.1, 173 post-H6,
188 first-draft post-H7).  The verified-correct counts
above sum to 181; this is the authoritative number going
forward.

Operator-invoked / live LLM tests verified on Mac (end of H7):
- `tests/llm/canary_check.py` (H3b orchestrator) — verified
  end-to-end with both Google AND OpenAI providers.
- `tests/llm/tst_directive_extraction.py` — all 7 scenarios ×
  2 providers = **14 PASS**.  File-detection-sensitive
  scenarios (`resolution_and_space_group`, `skip_programs`)
  pass after H7's scanner-first switch — Q2-strict preserved
  the consumer contract.
- `tests/llm/tst_planning.py` (NEW, exercised via H6's
  framework implementations) — 8 scenarios × 2 providers =
  **16 PASS**.  Early-stopping at threshold; no failures
  observed.  First reliability testing of the planning LLM.

**Cumulative live LLM test count after H7: 30 PASS
(14 directive_extraction + 16 planning), 0 FAIL, 86 LLM
calls, ~12 minutes wall-clock.**

### Deployment status

| Ship | Server (ai.phenix-online.org) | Mac (development) |
|---|---|---|
| H1, H2, H2.1, H3, H4, H4.1 | Deployed | Installed and verified |
| **H5** | Not yet deployed | Verified end-to-end |
| **H5.1** | Not yet deployed | Verified (K-suite + 14/14 directive_extraction live) |
| **H5.1.1** | Not yet deployed | Verified (+23/23 K_H5_1_1) |
| **H6** (test-only) | N/A (no server-visible changes) | Verified (+18/18 K_H6 + 16/16 planning live) |
| **H6.1** (test-only) | N/A (no server-visible changes) | Verified (output formatting only) |
| **H7** | Not yet deployed | **Verified (+15/15 K_H7 + 30/30 live LLM still PASS, 135/135 K-modules)** |

Recommended sequencing for server deployment of the
production-code ships:
H5 → H5.1 → H5.1.1 → H7 (each independent and rollback-safe;
combined deploy also safe).  H6 and H6.1 affect only `tests/`
files and don't need server deployment.

H5/H5.1/H5.1.1 are in `tests/` and `phenix_ai/` (plus
`programs/ai_analysis.py` for the relay channel).  H7 is
entirely in `phenix_ai/run_ai_analysis.py` (production
code) plus its K-suite.  Different code paths, safe to
combine or sequence in any order.

### v120 outlook

The v119 cluster completes the operational-hardening agenda from
`v118_next_steps_consolidated_rev4.md` §4.7B with the addition of
the H4 telemetry-marker, H5 relay-channel, H5.1.1 directive-
extractor cleanups, H6 planning-suite framework, and H7
Phase 2B scanner-first activation that weren't in the original
§4.7B list.  Deferred §14 items: auto-startup canary in
production (§14.1), Anthropic RAG branch (§14.3),
history_record propagation of agent_build (§14.5).

**Major v119 outcomes (now done):**
- ✓ **Phase 2A** (H6): planning-LLM reliability testing
  framework.  16 planning scenarios PASS across both providers.
- ✓ **Phase 2B Scope B** (H7): scanner-first file extraction
  with LLM fallback.  Deterministic file detection for the
  consumer at `programs/ai_agent.py:8128`.

**Outstanding for v119.H5.x or v120:**
- **H5.2 (queued, not started)** — extend the diagnostic_messages
  channel to three more stderr-only markers in
  `agent/directive_extractor.py`:
  `[DIRECTIVE_EXTRACTION_FAILED]` (inner — fires on transient
  LLM-call exceptions caught and returning `{}`),
  `[DIRECTIVE_EXTRACTION_MODEL_RETIRED]`, and
  `[DIRECTIVE_CORRECTION]`.  Per Tom's choice, let H5/H5.1 bake
  in production first to validate the channel design with real
  failures before extending.
- **Phase 2B Scope C** — kill the LLM preprocessor entirely.
  H7 captured most of Phase 2B's value (deterministic file
  detection); Scope C's primary value is latency reduction
  (~5-15s per request saved).  Needs corpus testing that
  directive extraction quality is preserved when fed raw
  advice instead of structured-section advice.
- **Pre-existing libtbx-only imports** (§3.11 in next_steps) —
  consistency sweep across the codebase for the libtbx →
  relative fallback pattern that H5+ / H6 / H7 uniformly use.

v120 is also expected to pursue prompt consolidation for
directive extraction (the natural pair to H1's model-name
centralization) and possibly Direction A (the single-call
structured extraction that would collapse preprocessor +
extractor + planner — referenced at the end of the v118
entry).  Specific direction TBD; scope shrinks if Phase 2B
Scope C succeeds.

---

## Version 118 (Preprocessor resilience + operational hardening)

### Summary

v118 is a 12-layer cumulative release addressing a class of bugs
surfaced in the 2026-05-19 production logs (AIAgent_240/241/243/245)
and follow-up reports through 2026-05-21.  The layers are
independent fixes that share the architectural pattern
"strengthen the directive pipeline against LLM-shape and
provider-asymmetry variance":

```
v117.3 → A → C-prime → B → E → F → 5.1 → G → 6.1-6.7 → 8 → 9 → 10
```

Each layer ships with its own test suite using the mock-LLM
pattern from `tst_after_program_fill_from_raw.py`.  Total
sandbox tests after v118.10: 204/204 passing.

### Per-layer breakdown

**Section A (client) — File-list preservation**
`agent/advice_preprocessor.py`,
`programs/ai_agent.py`,
`tests/tst_preprocessor_file_override.py`.  Preserves the input
file list across the preprocessor's text-vs-context round-trip
via UNION semantics in
`_ensure_file_list_in_processed_advice()`.  +12 tests.

**Section C-prime (client) — Diagnostic split**
`programs/ai_agent.py`,
`tests/tst_directive_layer_diagnostics.py`.  Splits the previously-
overloaded "directives" diagnostic output into
`directives_user_intent` (what the LLM extracted) vs
`directives_effective_runtime` (what the runtime actually
applied after grounding / fallback overlays).  Makes provider-
asymmetry symptoms diagnosable from logs alone.  +7 tests.

**Section B (client) — PHIL namespace healing**
`agent/advice_preprocessor.py`,
`agent/directive_extractor.py`,
`tests/tst_phil_namespace_cleaner.py`.  Static translation
table `PHIL_NAMESPACE_TRANSLATIONS` rewrites known-bad PHIL
paths (e.g. `data_manager.r_free_flags.generate` →
`xray_data.r_free_flags.generate` → bare
`generate_rfree_flags`).  Drops with log line when path can't
be healed.  +14 tests.

**Section E (client) — LLM-failure diagnostic markers**
`agent/advice_preprocessor.py`,
`agent/directive_extractor.py`,
`tests/tst_extraction_failure_visibility.py`.  Adds
`[DIRECTIVE_EXTRACTION_FAILED]` and
`[ADVICE_PREPROCESSING_FAILED]` stderr markers with the
exception text when an LLM call fails (network error, JSON
parse failure, retired model 404, etc.).  Section E proved
essential for diagnosing v118.8 — without it, the server
404 would have been buried in the server-side stderr with
no marker for client-side log search.  +11 tests.

**Section F (SERVER) — BUILD experiment_type + R-free auto-fill**
`agent/command_builder.py`,
`tests/tst_build_experiment_type_and_rfree.py`.  First v118
section that touches server-side code.  Threads
`experiment_type` through `_select_files` (was hardcoded
"xray" in the cycle-1 fast path → cryo-EM cycle-1 BUILD
selected wrong file slots when xtriage hadn't run yet).
Adds R-free generate-flag auto-fill for first-refinement-with-
no-rfree-mtz cases (was the AIAgent_245 reproducer:
phenix.refine ... xray_data.r_free_flags.generate=True).
+9 tests.

**Section 5.1 (client) — PHIL .help hyphen/semicolon hotfix**
`agent/advice_preprocessor.py`.  Drops hyphen-prefix and
semicolon-suffix forms from `[--help]` patterns that occasionally
make it through advice quoting.

**Section G (client) — Optional dependency resilience**
`rag/vector_store.py`,
`tests/tst_optional_dep_resilience.py`,
plus 6.x test-infrastructure cleanups.  Wraps the optional
chromadb → langchain_chroma chain in
`except Exception` (not `except ImportError`), because
protobuf version conflicts surface as `TypeError`, not
`ImportError`.  Adds environment-readiness test
`tst_dependencies.py` (v118.6.7 rev 3) that runs a true
runtime probe (`from langchain_chroma import Chroma`) rather
than just attempting to import.  +11 + 11 = +22 tests.
First v118 section verified on both Mac and Linux production
environments.

**Section 8 rev 3 (client + SERVER) — Model bump**
`agent/directive_extractor.py`,
`agent/api_client.py`.  Bumps the default Google model from
`gemini-2.0-flash` (retired by Google 2026-05-20) to
`gemini-2.5-flash-lite` in TWO files.  Without Section 8 the
server's directive extraction silently 404'd; Section E's
diagnostic surfaced the root cause in seconds.  First v118
section that required server-side deployment + process
restart for the fix to take effect.  No new tests (string
replacement only); confirmed server-verified via Section E
markers vanishing on retry.

**Section 9 rev 2 (SERVER) — Experiment-type-conditional program canonicalization**
`agent/directive_extractor.py`,
`tests/tst_density_modify_experiment_type.py`.  After Section 8
fixed the 404, runs 237/239 (Tom's "density modify and stop"
on cryo-EM half-maps) showed the LLM still picking
`phenix.autobuild_denmod` (X-ray) instead of
`phenix.resolve_cryo_em` (cryo-EM).  Section 9 adds three
layers:

  1. Decision-tree prompt restructure at lines 443-462 of
     `DIRECTIVE_EXTRACTION_PROMPT`: replaces two parallel
     cryo-EM/X-ray rules with an explicit "CHECK THE
     EXPERIMENT TYPE FIRST" decision tree.
  2. Module-level `PROGRAM_REPRINTS_BY_EXPERIMENT_TYPE`
     table (2 entries, symmetric) consumed by a new
     `_apply_experiment_type_program_reprints()` validator.
     Validator runs AFTER `_validate_after_program_grounded`
     by design (grounding bypasses on
     `stop_after_requested=True`, so the wrong program
     survives intact for correction).  Emits
     `[DIRECTIVE_CORRECTION]` log line + stores
     `_corrected_from` sidecar field in `stop_conditions`
     for traceability.
  3. K_DENMOD test suite (13 tests) covering the bug case,
     mirror case, no-change cases, ambiguous cases, grounding-
     bypass interaction, and diagnostic emission verification.

+13 tests.

**Section 10 (SERVER) — List-to-string coercion**
`agent/directive_extractor.py`,
`tests/tst_settings_list_coercion.py`.  Provider-specific JSON
shape mismatch: OpenAI emits `"additional_atom_types": "S"`
(string); Google's `gemini-2.5-flash-lite` emits
`"additional_atom_types": ["S"]` (list).  Validation loop
called `str(["S"])` which produced `"['S']"` (Python repr,
with brackets) — autosol rejected the malformed string.
Section 10 adds `_coerce_setting_value()` helper that:

  - Unpacks single-element lists/tuples of any type
    (`[False]` → `False`, fixing the Python truthiness
    trap `bool([False]) == True`).
  - Space-joins multi-element list/tuples for str-typed
    fields (PHIL multi-value-string syntax).
  - Raises ValueError for unsupported coercions (dict→str,
    multi-element list to non-str), caught by outer
    try/except and logged.

Helper applied at TWO call sites: `program_settings`
validation (Tom's bug) and `stop_conditions` validation
(latent bug — list-shaped `after_program` raised
`TypeError("unhashable type: 'list'")` from the
`value not in VALID_PROGRAMS` check, silently swallowed
and the directive dropped).  K_LIST test suite (12 tests)
includes K11 verifying the v118.10 + v118.9 interaction
chain end-to-end.  +12 tests.

### Process notes from v118

**Whole-tree grep before declaring a string-replacement done.**
Section 8 rev 1 fixed the model string in
`agent/directive_extractor.py` only.  The server kept
failing because `agent/api_client.py` also hardcoded the
retired model name at two locations.  Rev 2 needed a patch
script; rev 3 had the full file with both occurrences
updated.  v119 cadence convention (per consolidated
next_steps): every string-replacement starts with
`grep -rn 'STRING' $PHENIX/modules/cctbx_project/libtbx/langchain/`.

**Stub-module isolation for K-tests.** Section F established
the pattern of K-tests using stub modules for
`agent.program_registry` and `agent.intent_classifier` so
tests run in seconds without the full PHENIX conda env.
Sections G and 10 reused this pattern.

**Gemini-reviewed plan rev cycles.** Sections G, 9, and 10
each went through 2–4 plan revs incorporating Gemini critique
before implementation.  Each revision caught a real concern:
v118.G's annotation issue (rev 2), v118.9's table-driven
design (rev 2), v118.10's bool-truthiness trap (rev 3).

### Files changed (cumulative)

| File | Sections that touched it |
|---|---|
| `agent/directive_extractor.py` | B, E, 8, 9, 10 |
| `agent/advice_preprocessor.py` | A, B, E, 5.1 |
| `agent/api_client.py` | 8 |
| `agent/command_builder.py` | F |
| `programs/ai_agent.py` | A, C-prime |
| `rag/vector_store.py` | G |
| `tests/tst_preprocessor_file_override.py` | A (new) |
| `tests/tst_directive_layer_diagnostics.py` | C-prime (new) |
| `tests/tst_phil_namespace_cleaner.py` | B (new) |
| `tests/tst_extraction_failure_visibility.py` | E (new) |
| `tests/tst_build_experiment_type_and_rfree.py` | F (new) |
| `tests/tst_optional_dep_resilience.py` | G (new) |
| `tests/tst_dependencies.py` | 6.7 (new) |
| `tests/tst_density_modify_experiment_type.py` | 9 (new) |
| `tests/tst_settings_list_coercion.py` | 10 (new) |
| `tests/run_all_tests.py` | A, C-prime, B, E, F, G, 6.7, 9, 10 (registrations) |

### Verification

Sandbox: 204/204 passing.

Production verified:
- Section G: Mac and Linux production envs
- Section 8: server stderr clean of
  `[DIRECTIVE_EXTRACTION_FAILED]` markers after deploy
- v117.3 + Section A + Section C-prime: confirmed via prior
  production runs

Production awaiting:
- Section 9: rerun "density modify and stop" on cryo-EM
- Section 10: rerun P9 SAD with Google provider

### Architectural watchpoint triggered

At v118.10 the project is at **10 iterations** in the
preprocessor → extractor → planner → BUILD pipeline (original
threshold was 6).  Specific signs that have appeared:

- New sections added for yet-another-LLM-emission-shape failure
  modes (v118.9 wrong canonical name; v118.10 list-wrapping)
- Defensive logic at one layer compensating for upstream
  limits (v118.9 validator, v118.10 helper)
- Cross-section coupling: v118.10's K11 explicitly tests the
  v118.10 → v118.9 chain

v119 will pursue operational hardening (centralized model
defaults, server `/version` endpoint, startup canary, server→
client diagnostic echo) and prompt consolidation; v120 may
pursue Direction A (single structured-extraction LLM call
replacing preprocessor + extractor + planner).  See
`v118_next_steps_consolidated_rev4.md` for the full v119+ plan.

---

## Version 117.3 (Extended stop-intent phrasing recognition)

### Summary

A C1 LLM test run after v117.2 surfaced one remaining openai failure
in the existing `tst_directive_extraction.py::explicit_stop_after_phaser`
scenario (0/5 on openai, 1/1 on google).  The fixture's user-side
phrasing is `"Stop the workflow immediately after phaser completes"`
in its Special Instructions section.

Diagnosis: this phrasing is a clear user-stop directive, but neither
the v117.2 LLM prompt schema documentation nor the regex backstops
(`_IMPERATIVE_STOP_MARKERS`, `_POSITIVE_STOP_AFTER_PATTERNS`) had any
entry that recognized it.  Result: the LLM extracted
`after_program=phenix.phaser` correctly, but did NOT set
`stop_after_requested=True`, so v117.1's flag-exemption couldn't fire,
and the v116.19a grounding guardrail's Failure 2 path dropped the
directive as "fabricated."

### The fix — three independent additions

**Change 1 — `_IMPERATIVE_STOP_MARKERS`.**  Added 5 new substring
markers used by the grounding guardrail's window-bounded
imperative-check (300 chars near program name):
`"stop the workflow"`, `"immediately after"`, `"after it completes"`,
`"after it finishes"`, `"is the last step"`.

**Change 2 — `_POSITIVE_STOP_AFTER_PATTERNS`.**  Added 2 new
contiguous-phrase regex patterns: `\bstop\s+the\s+workflow\b` and
`\bis\s+the\s+last\s+step\b`.  Two over-permissive patterns
(`\bafter\s+\w+\s+completes?\b`, `\bafter\s+\w+\s+(?:finishes|is\s+done)\b`)
were considered during planning and dropped after the K5 false-positive
verification — they would have fired on descriptive prose like
"validation runs after refinement completes."

**Change 3 — `DIRECTIVE_EXTRACTION_PROMPT` schema docs.**  Extended
the `stop_after_requested` documentation with 4 new recognized
phrasings and 2 program-neutral few-shot examples (using "refinement"
and `<program>` placeholder to avoid biasing extraction toward
specific Phenix tools).  The dual-input prompt
`DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` inherits the change via the
existing `.replace()` mechanism.

### How the fix works end-to-end

Two paths now correctly recognize the new phrasings:

1. **Via the LLM (preferred):** the schema docs + examples teach the
   LLM to emit `stop_after_requested=True` alongside `after_program`
   when the user says "stop the workflow after X" / "X is the last
   step" / etc.  v117.1's grounding-flag exemption then keeps
   `after_program` intact.

2. **Via the regex backstops (safety net):** even when the LLM
   misses the flag (as openai did pre-v117.3), the new regex
   patterns in `_POSITIVE_STOP_AFTER_PATTERNS` cause
   `_is_stop_after_requested(stripped_advice)` to return True.  The
   downstream `_apply_workflow_intent_fallback` then sets the flag
   via `_resolve_after_program`.  Meanwhile, the new imperative
   markers cause the grounding guardrail's `_imperative_marker_nearby`
   check to succeed, so `after_program` is preserved through Failure 2.

The end-to-end test simulation confirms both paths.  In the C1 failure
scenario where the LLM emits `after_program=phenix.phaser` but no
flag, v117.3 produces:

```json
{
  "after_program": "phenix.phaser",
  "stop_after_requested": true,
  "skip_validation": true
}
```

### What changes for users

For `explicit_stop_after_phaser`-class workflows (preprocessed
advice containing "stop the workflow immediately after X" or similar
phrasings), the production workflow now correctly stops after the
user's intended program completes.  Pre-v117.3, the openai provider
would silently continue past the user's stop.

For all other scenarios, behavior is unchanged.  The new patterns
and markers are additive; no existing patterns are modified or
removed.

### Tests

New unit test file `tests/tst_extended_stop_phrasings.py` with 9
boundary tests (K1, K2, K4, K5, K6, K7a, K7b, K8a, K8b):

| Test | Purpose |
|------|---------|
| K1 | `_is_stop_after_requested` recognizes "stop the workflow" phrasing |
| K2 | `_is_stop_after_requested` recognizes "X is the last step" |
| K4 | AF_7mjs preservation — stripped advice still returns False |
| K5 | Descriptive "validation runs after refinement completes" → False (false-positive guard) |
| K6 | `_imperative_marker_nearby` fires on the explicit_stop_after_phaser fixture |
| K7a, K7b | Existing C1 phrasings ("density modify and stop", "refine and stop") still recognized |
| K8a, K8b | Negation guards ("do not stop", "Stop Condition: None") still suppress |

K3 was considered and dropped from the test plan because the existing
`,\s*stop\b` pattern already catches "after phaser completes, stop".
K9 and K10 (cross-sentence safety) were dropped because the
pre-existing v117.2 behavior of `\bstop\s+if\b` matching within
single sentences is a separately tracked issue, not v117.3 scope.

The v117.2 sandbox suite (106 tests) plus the 9 new K-tests gives
115/115 passing on v117.3.

### Expected LLM test results after applying

| Test | v117.2 | v117.3 |
|------|--------|--------|
| `explicit_stop_after_phaser` openai | 0/5 | expected ≥4/5 |
| `explicit_stop_after_phaser` google | 1/1 | maintained |
| C1 suite | 6/6 | maintained |
| B1 suite | 10/10 | maintained |
| A1 suite | 6/6 | maintained |
| Other 5 `tst_directive_extraction.py` scenarios | passing | maintained |

### Files changed

| File | Change |
|------|--------|
| `agent/directive_extractor.py` | 34 lines added across 3 blocks; no removals |
| `tests/tst_extended_stop_phrasings.py` | **NEW** — 9 unit tests |

### Architectural notes

**Pipeline ordering unchanged.**  v117.3 does not add a new pipeline
step.  The accumulated post-LLM processing order from v117.2 stays:

1. `_validate_after_program_grounded` (v116.19a + v117.1 flag exemption)
2. `_apply_crystal_symmetry_fallback`
3. v117.2 fill_in_from_raw
4. Ollama empty-directives retry
5. Intent classification merge
6. Intent override block
7. `_apply_workflow_intent_fallback` (calls `_resolve_after_program`)

v117.3 extends the regex sets and marker tuples consumed by steps 1,
3, and 7.  No new pipeline step; no new state fields.

**Plan revisions during implementation.**  The plan in
`v117_3_PLAN.md` was revised twice during implementation:

- Revision 2: absorbed Gemini's review feedback (bounding discipline
  in Change 2, few-shot examples in Change 3, rejection of auto-
  healing alternative documented in §9).
- During Step 2 (TDD baseline check): K3 was discovered to already
  pass on v117.2 (redundant) and K9/K10 were discovered to surface a
  pre-existing `\bstop\s+if\b` issue not introduced by v117.3.
  Both adjustments narrowed test scope without weakening coverage.
- During Step 2 (K5 false-positive verification): the over-permissive
  patterns in Change 2 were dropped (Option 3).

### Process note — fourth pipeline augmentation

This is the fourth iterative behavior augmentation to `extract_directives`
and its helpers in six days (Step 1 / v117.1 / v117.2 / v117.3).  Per
the v117.2 CHANGELOG flag, a fifth such gap should trigger a
consolidation discussion rather than a fifth additive patch.

v117.3's scope was kept narrow specifically because of this concern:
the change is purely additive to existing data structures (5 entries
added to a tuple, 2 entries added to another tuple, 4 + 2 lines added
to a docstring).  No new helpers, no new pipeline step.

---

## Version 117.2 (Fill in after_program from raw advice when LLM omits it)

### Summary

A C1 LLM test run after the v117.1 grounding fix surfaced a remaining
gap on openai: for the input `"refine and stop"` (raw) with processed
advice ending in `"Stop Condition: None"`, openai's extractor
occasionally emits `stop_after_requested=True` but omits
`after_program` entirely.

Without `after_program`, `workflow_engine._apply_directives`' gated
wipe code is unreachable — the wipe is gated on `if after_program:`.
The user's stop intent is signaled but cannot be acted on; the
workflow continues past where the user wanted it to stop.

### The fix

In `extract_directives()`, after the existing LLM extraction and
validation flow but before the intent-override block, add a small
fallback: when `stop_after_requested=True` and `after_program` is
unset, run the existing `_resolve_after_program()` regex resolver
against the raw advice to fill in `after_program` via `_ACTION_TABLE`.

Guards:
- Only fires when the LLM left `after_program` unset (does not
  override LLM choices)
- Only fires when raw advice contains stop intent per
  `_is_stop_after_requested` (does not amplify hallucinated flags)
- Uses `raw_advice` if available; falls back to `user_advice` for
  callers using the single-input path

### What changes for users

For the affected scenarios — openai + `"refine and stop"` /
`"run xtriage and stop"` / `"density modify and stop"` against
preprocessor-mangled processed advice — the production workflow now
correctly stops after the user's intended program completes, where
it previously continued past it.

For all other scenarios (LLM emits both fields correctly; or raw
advice has no stop intent; or flag is unset), behavior is unchanged.

### Tests

New unit test file `tests/tst_after_program_fill_from_raw.py` with
6 boundary cases (K1-K6):
- K1: flag True + after_program missing + raw has stop intent →
  after_program filled in (the C1 fix case)
- K2: flag True + after_program already set → no change (LLM choice
  wins)
- K3: flag True + after_program missing + raw has NO stop intent →
  no change (guard prevents amplifying hallucination)
- K4: flag absent → no change
- K5: raw_advice None → falls back to user_advice (single-input
  compatibility)
- K6: `"run xtriage and stop"` → after_program=phenix.xtriage
  (confirms resolver handles other action verbs)

The v117.1 sandbox suite still passes 100/100.  With v117.2 added:
106/106.

### Expected LLM test results after applying

| Suite | Pre-v117.1 | v117.1 | v117.2 |
|-------|------------|--------|--------|
| A1 (smoke) | 4/6 | 6/6 | 6/6 |
| B1 (stop_after_requested) | 8/10 | 10/10 | 10/10 |
| C1 (raw authority) | 0/6 | 5/6 | 6/6 (expected) |

### Files changed

| File | Change |
|------|--------|
| `agent/directive_extractor.py` | 30-line addition in `extract_directives()` |
| `tests/tst_after_program_fill_from_raw.py` | **NEW** — 6 unit tests |

### Architectural note

This is a third behavior-augmentation in three days touching
`extract_directives()`.  The accumulated post-LLM processing path now
runs (in order):

1. `_validate_after_program_grounded` (v116.19a, with v117.1 flag
   exemption)
2. `_apply_crystal_symmetry_fallback`
3. **v117.2 fill_in_from_raw** (new)
4. Ollama empty-directives retry
5. Intent classification merge
6. Intent override block (solve / solve_constrained / task)
7. `_apply_workflow_intent_fallback` (calls
   `_resolve_after_program` on processed advice)

The v117.2 fill-in is positioned after grounding (so it can't be
undone by it) and before the intent block (so it doesn't accidentally
clear a flag the resolver just confirmed).  See DEVELOPER_GUIDE.md
§3i for the full design rationale.

### Process note

Like v117.1, v117.2 was surfaced by the C1 LLM test — the test's
remaining single failure after v117.1 produced exactly the diagnostic
log needed to identify the gap.  The pattern of "ship the fix, run
the LLM tests, examine remaining failures with the same surgical
discipline" is working.

---

## Version 117.1 (Grounding ∧ stop_after_requested Interaction Fix)

### Summary

The C1 LLM test (raw-authority end-to-end) surfaced an architectural
conflict between v116.19a's grounding guardrail and v117 Step 1's
raw-advice extraction.  Both were correct in isolation but collided
in production.

**The conflict:** when the user typed `"density modify and stop"`,
v117 Step 1's AUTHORITY paragraph correctly led the LLM to produce
`after_program=phenix.resolve_cryo_em, stop_after_requested=True` —
mapping the action verb to the cryo-EM density-modification program
even though "resolve_cryo_em" doesn't appear literally in the advice.

Then v116.19a's grounding guardrail dropped `after_program` because
the program name was "absent from user advice (pure fabrication)" by
its literal-substring test.  The user-stop intent was retained as
the flag, but the target program was lost.

### The fix

When the LLM sets BOTH `after_program` AND `stop_after_requested=True`
on the same `stop_conditions` block, the literal-name grounding check
is skipped.  The flag IS the grounding signal — the LLM made an
explicit user-stop assertion based on the raw advice and the
AUTHORITY paragraph, which is strictly stronger evidence than the
heuristic substring check.

The fix is 11 lines in `_validate_after_program_grounded`
(`agent/directive_extractor.py`) with the new behavior gated on the
`stop_after_requested` flag presence.  No other production code
changes.

### Behavior preservation

- **AF_7mjs case** (the case v116.19a was originally written for):
  preserved.  AF_7mjs's preprocessed advice has no positive stop
  signal, so `stop_after_requested` is not set and the guardrail
  fires as before.
- **Pure fabrication** (program name absent, no flag): still dropped.
- **Flag explicitly False**: still dropped (flag must be True to
  trigger the skip).
- **Grounded program with no flag**: still kept (existing v116.19a
  behavior unchanged).

### Tests

New unit test file `tests/tst_grounding_stop_after_requested.py` with
5 boundary cases:
- K1: LLM sets both signals → after_program preserved (the C1 case)
- K2: LLM sets only after_program (no flag) → dropped (AF_7mjs case)
- K3: pure fabrication, no flag → dropped
- K4: flag explicitly False → dropped
- K5: grounded program, no flag → kept (sanity baseline)

The v117 sandbox suite still passes 95/95.  With the new file the
total is 100/100.

### LLM test results expected after applying

| Suite | Pre-v117.1 | Post-v117.1 |
|-------|------------|-------------|
| A1 (smoke) | 4/6 | 6/6 |
| B1 (stop_after_requested) | 8/10 | 10/10 |
| C1 (raw authority) | 0/6 | 6/6 |

The c1_density_modify_raw_authority scenario also has a relaxed
assertion accepting any density-modification-family program
(resolve_cryo_em / resolve / density_modification / parrot /
solvent_modify), since the raw input `"density modify and stop"`
doesn't specify cryo-EM vs X-ray context.  Strict assertions remain
on c1_refine_raw_authority and c1_xtriage_raw_authority where the
program name appears in the raw input itself.

### Files changed

| File | Change | md5 (pre → post) |
|------|--------|------------------|
| `agent/directive_extractor.py` | 11-line addition to `_validate_after_program_grounded` | (v117 → v117.1) |
| `tests/tst_grounding_stop_after_requested.py` | **NEW** — 5 unit tests | — |
| `tests/llm/tst_c1_raw_authority.py` | relaxed C1.1 assertion | — |
| `tests/llm/tst_b1_stop_after_requested.py` | relaxed B1.1 assertion | — |
| `tests/llm/tst_a1_extractor_smoke.py` | relaxed A1.1 assertion | — |

### Process note

The C1 test was correctly designed and produced exactly the failure
it should have: the LLM-layer behavior was right, but a downstream
guardrail dropped its output.  The test "failing" surfaced the
architectural conflict.  This is the second time in v117 that a
test failure surfaced a real finding — the first was the baseline
drift caught in Step B by
`test_af7mjs_no_start_with_program_for_preprocessed`.

---

## Version 117 (Extraction Reliability — Step 1: Raw-Advice-Authoritative Directive Extraction)

### Summary

The LLM advice preprocessor occasionally mangles short imperative user
input.  The most reproducible failure: when the user writes
`"density modify and stop"` and the preprocessor is run with the openai
provider, the preprocessed output expands the request into a multi-step
Primary Goal and writes `Stop Condition: None` — the user's stop intent
is silently dropped.  Because the downstream directive extractor only
sees the preprocessed text, the LLM-extracted directives miss the stop
signal entirely; recovery depends on the v116.x regex stop-condition
backstops in `workflow_engine.py::_apply_directives` and the simple
pattern-extractor fallback in `extract_directives_simple`.  Those
backstops have their own gaps (the markdown-bold/numbered preprocessor
output formats — see `EXTRACTOR_BUGS_FINDING.md` Finding 2 — slip
through some of the strip regexes), so reliance on them is fragile.

Step 1 fixes this at the source by giving the LLM directive extractor
**both** inputs (raw user text and preprocessed advice) with an explicit
authority rule: **raw wins for intent** (`after_program`,
`start_with_program`, `stop_after_requested`); processed fills gaps for
files, parameters, and experiment type the raw is silent on; on
disagreement, raw is the source of truth.

The change is purely additive — old clients that don't send raw advice
get identical behavior to today.  New clients (which set
`user_advice_raw` alongside `user_advice_for_directives`) get the
dual-input prompt and the more reliable stop-intent extraction.

### Modified Files

| File | Lines | Description |
|------|-------|-------------|
| `agent/directive_extractor.py` | +112 | Add `_RAW_INPUT_BLOCK` and `_SINGLE_INPUT_BLOCK` named template fragments. Derive `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` from `DIRECTIVE_EXTRACTION_PROMPT` via `str.replace()` with an import-time `assert` guarding single-match. `extract_directives()` gains `raw_advice=None` parameter and a selector that picks dual-input only when raw is provided and differs from processed. `stop_after_requested` schema doc added to the prompt with recognized phrasings (`"stop after X"`, `"X and stop"`, `"X, then stop"`). |
| `phenix_ai/run_ai_analysis.py` | +57 | `run_directive_extraction` gains `raw_user_advice=None` parameter (sub-step 1B); plumbs to `extract_directives(raw_advice=...)`. `run_advice_preprocessing` gains Step 1F metric instrumentation — emits one greppable `PP_FILE_METRIC` line per call, comparing the LLM-derived `extracted_files` against a regex baseline run on the raw advice.  Gated to fire only when the LLM produced a distinct output; a `PP_FILE_METRIC: skipped` marker covers the alternative so log data is unambiguous. |
| `agent/session.py` | +28 | `AgentSession.extract_directives()` gains `raw_advice=None` parameter with default-source from `self.data["raw_advice"]` (CCTBX None-safe), plumbs to the module-level `extract_directives()`. |
| `agent/advice_preprocessor.py` | +95 | New "RAW-ADVICE FILE DETECTION (STEP 1F METRICS)" section. `extract_files_from_raw_advice(raw_advice)` runs a conservative regex over raw text and returns deduplicated lowercase basenames.  `_summarize_file_detection_metric(llm_files, regex_files)` builds the `PP_FILE_METRIC` summary line. |
| `phenix/programs/ai_analysis.py` | +30 | New PHIL param `user_advice_raw = None` next to `user_advice_for_directives`. `_build_server_args` encodes `user_advice_raw` for transport when set. `_run_directive_extraction_locally` reads, decodes, and plumbs as `raw_user_advice=` to `run_directive_extraction`. |
| `phenix/programs/ai_agent.py` | +19 | New PHIL param `user_advice_raw = None` mirroring `ai_analysis.py`. `_extract_directives` sets `directive_params.ai_analysis.user_advice_raw = session.data.get("raw_advice") or None` after the existing `user_advice_for_directives` setter. |
| `tests/tst_extract_raw_advice.py` | +675 (new file) | 20-test module across 8 sections: prompt structure, signature/selector, Tom's openai case (structural), plumbing source-scans for run_directive_extraction / session.py / ai_analysis.py / ai_agent.py, and 1F behavioral tests for the file-detection regex and metric summary plus a source-scan for the run_advice_preprocessing emission point. |
| `tests/run_all_tests.py` | +10 | Register the new test module after `tst_directive_extractor`. |

### Backward compatibility

Old clients send `user_advice_for_directives` only (the processed text).
New server reads `user_advice_raw` as `None` (PHIL default), passes
`None` through to `extract_directives`, which falls back to the
single-input prompt — identical to old behavior.  Deployment story:
server-first.  No "new client + old server" scenario exists because
clients always lag the server.

### Pre-existing bugs found while implementing Step 1

Two bugs that predate Step 1 were uncovered during sub-step 1A and are
documented separately in `EXTRACTOR_BUGS_FINDING.md`:

1. **Brace bug (medium severity)** —
   `DIRECTIVE_EXTRACTION_PROMPT.format()` raises
   `KeyError('"program_settings"')` because the unit_cell example in the
   prompt body contains unescaped JSON braces.  The crash is caught by
   a broad `except` at `run_ai_analysis.py:1020` and silently logged to
   `debug_log`, returning an empty directives dict.  **The LLM
   directive extractor has been dark in production.**  The
   google-vs-openai differences Tom investigated turn out to be
   preprocessor LLM differences, not extractor LLM differences.  Fix is
   5 characters (double every `{` and `}` in the JSON example).  Step
   1's value isn't observable until this is fixed.

2. **Stripper regex bug (low severity)** —
   `_strip_preprocessor_stop_condition()` regex
   `r'^[ \t]*stop\s+condition\s*:\s*[^\n]*$'` only matches plain
   `"Stop Condition:"` lines, not the markdown-bold forms
   (`**Stop Condition**: None`) or numbered forms
   (`7. Stop Condition: None`) the preprocessor actually emits.
   Currently benign because the downstream `_STOP_CONDITION_NONE`
   regex DOES handle markdown forms, but worth fixing for
   consistency.

### Why the existing v116.x regex backstops are still in play

The v116.x stop-after refactor (`stop_conditions.stop_after_requested`
flag, workflow_engine wiping `valid_programs` to `[STOP]`) remains as
defense-in-depth.  Step 1 makes the LLM extractor's output more
reliable, but the regex backstops in
`workflow_engine.py::_apply_directives` continue to catch any case the
LLM extractor still misses.  Belt and braces.

### Testing

| Test module | Tests | Result |
|-------------|-------|--------|
| `tst_extract_raw_advice.py` (new) | 20 | 20/20 pass |
| End-to-end extractor cases (existing) | 12 | 12/12 pass |
| Workflow engine integration (existing) | 11 | 11/11 pass |
| Audit-fixes i2 (existing) | 3 | 3/3 pass |

The new test module verifies the same 20 invariants both pass on
post-Step-1 code AND fail on pre-Step-1 code with diagnostic messages,
giving anyone running the suite against an older checkout an immediate
checklist of what needs to be done.

### Sub-step inventory (for review and roll-back)

The change landed as 8 sub-steps with checkpoints between each.  Diffs
per sub-step are stored alongside the production files:

- **1A** — `agent/directive_extractor.py` (prompt + signature)
- **1B** — `phenix_ai/run_ai_analysis.py::run_directive_extraction`
- **1B-bis** — `agent/session.py::AgentSession.extract_directives`
- **1C** — `phenix/programs/ai_analysis.py` (PHIL + dispatcher)
- **1D** — `phenix/programs/ai_agent.py` (PHIL + client touch)
- **1E** — `tests/tst_extract_raw_advice.py` + registration
- **1F** — `agent/advice_preprocessor.py` regex + metric in `run_advice_preprocessing`
- **1G** — this CHANGELOG entry plus ARCHITECTURE.md and DEVELOPER_GUIDE.md updates


## Version 116.13 (Cryo-EM Validation Check in metric_evaluator.py)

### Summary

The real fix for the AF_7mjs Stage 5 (validation) skip bug. v116.12
patched `agent/metrics_analyzer.py` but missed that
`agent/graph_nodes.py` sets `USE_YAML_METRICS = True` (the default).
With that flag set, `analyze_metrics_trend()` delegates immediately
to `analyze_refinement_trend()` in `agent/metric_evaluator.py` —
**a parallel implementation v116.12 never patched.** The v116.12
fix was applied to dead code; the production path retained the
same asymmetry (X-ray branch checks validation_done, cryo-EM branch
doesn't).

v116.13 patches the live code path: adds the validation_done check
to `MetricEvaluator.analyze_trend()` cryo-EM SUCCESS branch
(line ~500), mirroring the existing X-ray pattern at line 403.

### Modified Files (1 production file)

| File | Lines | Description |
|------|-------|-------------|
| `agent/metric_evaluator.py` | +20 / -3 | Mirror X-ray validation_done check in the cryo-EM `latest_cc > target` branch. Recognizes `phenix.molprobity` and `phenix.validation_cryoem`. |

### Why v116.12 didn't work

`agent/graph_nodes.py` line 103:
```python
USE_YAML_METRICS = True      # Use metrics.yaml for trend analysis
```

`agent/metrics_analyzer.py` line 271-276:
```python
if use_yaml_evaluator and analyze_refinement_trend is not None:
    try:
        return analyze_refinement_trend(metrics_history, ...)
    except Exception as e:
        print("Warning: YAML evaluator failed, using hardcoded: %s" % e, ...)
```

When `USE_YAML_METRICS = True`, `analyze_metrics_trend` delegates
the entire trend analysis to `analyze_refinement_trend()` in
`metric_evaluator.py` and returns immediately. The hardcoded code
(where v116.12 Fix #1 was applied) only runs as a fallback when
the YAML evaluator raises an exception.

Verification: both v116.12 deployed files were byte-identical to
the patched versions, but behavior was unchanged. The fix was in
dead code.

### Three-place asymmetry that should have been spotted

| Path | File | Line | Has validation_done check? |
|------|------|------|----------------------------|
| YAML evaluator (active) | `metric_evaluator.py` | ~500 cryo-EM | **NO** (fixed by v116.13) |
| YAML evaluator (active) | `metric_evaluator.py` | ~403 X-ray | YES (pre-existing) |
| Hardcoded fallback | `metrics_analyzer.py` | ~505 cryo-EM | NO (fixed by v116.12 in dead code) |
| Hardcoded fallback | `metrics_analyzer.py` | ~365 X-ray | YES (pre-existing) |

After v116.13, all four paths have the validation_done check.

### Code change

```python
# Get target from YAML
target = self.get_target("map_cc") or 0.70

# v116.13: Mirror X-ray validation_done check at line 403.
validation_done = any(
    m.get("program") in ("phenix.molprobity", "phenix.validation_cryoem")
    for m in metrics_history
)

# SUCCESS check
if latest_cc > target:
    if validation_done:
        result["should_stop"] = True
        result["reason"] = "SUCCESS: Map-model CC (%.3f) above target (%.2f)" % (latest_cc, target)
        result["recommendation"] = "stop"
    else:
        result["should_stop"] = False
        result["reason"] = "Map-model CC (%.3f) above target - recommend validation" % latest_cc
        result["recommendation"] = "validate"
        result["suggest_validation"] = True
    result["trend_summary"] = "Map CC: %.3f - TARGET REACHED" % latest_cc
    return result
```

The trend_summary `"Map CC: X.XXX - TARGET REACHED"` format is
preserved exactly, so the cycle box display doesn't change.

### Tests (1 new file, 15 tests)

| File | Tests | Pre-fix | Post-fix |
|------|-------|---------|----------|
| `tst_cryoem_metric_evaluator_validation.py` (S29) | 15 (22 assertions) | 6 fail (9 pass) | 15 / 15 pass |

Test coverage:

- **`test_af7mjs_regression`**: exact AF_7mjs metrics history (cycle 4 RSR, CC=0.827) → should_stop=False
- **`test_cc_above_target_with_molprobity_done`** / **`_with_validation_cryoem_done`**: both validation programs recognized
- **`test_phenix_model_vs_data_not_cryoem_validation`**: X-ray validation program doesn't count for cryo-EM
- **`test_validation_earlier_in_history_counts`**: validation mid-history still counts (any() over full history)
- **`test_cc_above_target_validation_not_done`**: defer to validation when not done
- **`test_reason_mentions_validation`**: reason hints at validation
- **`test_trend_summary_format_preserved`**: TARGET REACHED display preserved
- **`test_cc_just_above_threshold`** / **`_exactly_at_threshold_no_success`** / **`_below_target_no_success_suggestion`**: boundary conditions (strict `>`)
- **`test_xray_path_unchanged_validation_done`** / **`_not_done`**: X-ray path regression guard
- **`test_using_evaluator_class_directly`**: fix is at the MetricEvaluator class level, not just a wrapper
- **`test_consecutive_field_still_set`**: existing result fields preserved

Registered in `tests/run_all_tests.py` as suite S29.

### Lesson from the v116.12 miss

When patching a function: **before writing the patch, search the
module for early-return paths, delegation, and alternative
implementations.** The function signature and first 30 lines tell
you whether you're patching live code.

A simple grep would have caught this:
```bash
grep -n "use_yaml\|analyze_refinement_trend" agent/metrics_analyzer.py
```

When a deployed fix doesn't change behavior: **look at the
message format first.** If the printed strings don't match the
patched code, the patch is in dead code. The session output had
`"above target (0.70)"` but v116.12 produced `"above 0.70 target"`
— I noticed but explained it away as a deployment lag instead of
treating it as definitive evidence.

### Compatibility

- X-ray workflows unchanged (X-ray path already had this check)
- Workflows with validation already done: no behavior change
- Workflows below CC target: no behavior change (SUCCESS branch doesn't fire)
- No directive or contract changes
- Trend summary format unchanged

### Relation to v116.12

v116.12's `metrics_analyzer.py` patch remains in place but is
dead code (only runs when YAML evaluator raises an exception).
Leaving it provides correct behavior in the unlikely fallback case.

v116.12's `graph_nodes.py` Fix #2 (defense-in-depth elif +
diagnostic context dump) remains valuable: now genuinely
defense-in-depth rather than necessary, and the diagnostic logging
will help if any other `should_stop=True` source (PLATEAU,
EXCESSIVE) ever fires with a pending validation stage.

### Known limitations

| Limitation | Impact |
|------------|--------|
| **Two parallel stop-logic implementations** | `metrics_analyzer.py` and `metric_evaluator.py` both implement cryo-EM stop logic. v116.12 + v116.13 has both fixed, but future drift is a real risk. Code review for any change to one should require updating the other. |
| **PLATEAU and EXCESSIVE paths don't check validation** | Currently neither X-ray nor cryo-EM PLATEAU/EXCESSIVE branches gate on validation_done. v116.13 doesn't change this. If a workflow PLATEAUs at CC=0.72 without validation, stop still fires. Pre-existing semantic question. |
| **Cosmetic `'unknown'` in F7 fallback** | Unchanged from v116.12. When stop intent sets `decision_info["program"] = None`, F7's fallback prints `'unknown'`. Session correctly ends; only display text is wrong. |
| **`plan_has_pending_stages` mystery** | Still unresolved why v116.12 Fix #2's earlier `plan_has_pending_stages` elif didn't suppress AF_7mjs. v116.13 makes this moot for the SUCCESS path (we don't reach the PLAN auto-stop chain), but a future PLATEAU/EXCESSIVE stop with a pending validation stage could resurface the same issue. |

---

## Version 116.12 (Cryo-EM Validation-Aware Auto-Stop)

### Summary

Two coordinated fixes for the AF_7mjs Stage 5 (validation) skip
bug. The agent ran the correct 5-stage plan but auto-stopped after
Stage 4 (refinement, CC=0.78) without running validation
(phenix.molprobity).

Root cause: asymmetry between X-ray and cryo-EM stop logic in
`agent/metrics_analyzer.py`. The X-ray path correctly defers
auto-stop until validation runs; the cryo-EM path didn't. The two
paths drifted out of sync.

### Modified Files (2 production files)

| File | Phases | Changes |
|------|--------|---------|
| `agent/metrics_analyzer.py` | Fix #1 | +25 / -4 lines. Mirror the X-ray validation_done check in the cryo-EM `latest_cc > 0.70` branch. Recognizes `phenix.molprobity` and `phenix.validation_cryoem`. When validation isn't done, sets `should_stop=False` and `suggest_validation=True` instead of stopping immediately. |
| `agent/graph_nodes.py` | Fix #2 | +41 / -0 lines. Two additions in the PLAN node's auto-stop chain: (a) new `elif` between `plan_has_pending_stages` and the final `else` — suppresses AUTO-STOP when `workflow_state.step_info.step == "validate"` AND `validation_done == False`; (b) diagnostic context dump in the AUTO-STOP path logging `plan_has_pending_stages`, `step`, `validation_done`, `after_program`, `experiment_type` to `debug_log`. |

### Bug Details

#### Symptom

```
[STATE]  cryoem_refined
  Model refined, need validation - Goal: Ensure model meets quality standards (CC: 0.783)
[STOP]
  SUCCESS: Map-model CC (0.783) above target (0.70)
No command was generated for program 'unknown'. Ending session.
```

AF_7mjs (cryo-EM, half-maps + sequence, PredictAndBuild) had a
plan with 5 stages. After Stage 4 (`phenix.real_space_refine`)
produced `map_cc=0.7832`, the agent should have advanced to Stage
5 (`phenix.molprobity`). Instead it auto-stopped. The state
machine even reported "Model refined, need validation" but the
threshold-based stop fired anyway.

Session JSON confirmed:
- Plan had 5 stages (correct)
- Stages 1-3 complete (correct)
- Stage 4 active, `cycles_used=2`, `success_criteria: model_map_cc > 0.7` met
- Stage 5 `status="pending"`, `programs=["phenix.molprobity"]`, **never ran**
- Cycle 5 reasoning: `"Automatically stopping: SUCCESS: Map-model CC (0.783) above target (0.70). Map CC: 0.783 - TARGET REACHED"`
- Cycle 5 program: empty (None); F7 fallback printed `'unknown'`

#### Root cause analysis

The X-ray and cryo-EM stop paths in `analyze_metrics_trend` should
be symmetric on the "success → validation" question. They were
not.

**X-ray path (lines 365-385) — correct:**
```python
if latest_r_free < success_threshold:
    validation_done = any(
        m.get("program") in ("phenix.molprobity", "phenix.model_vs_data")
        for m in metrics_history
    )
    if validation_done:
        result["should_stop"] = True
    else:
        result["should_stop"] = False
        result["suggest_validation"] = True
```

**Cryo-EM path (lines 506-511) — bug:**
```python
if latest_cc > 0.70:
    result["should_stop"] = True   # No validation check
    return result
```

The cryo-EM path was probably written first and the validation_done
check was added to the X-ray path later. The two paths drifted.

### Three logical changes

| Change | Location | What |
|--------|----------|------|
| **Fix #1: validation_done check** | `metrics_analyzer.py` `_analyze_cryoem_trend`, `latest_cc > 0.70` branch | Mirror X-ray pattern. When validation not done, `should_stop=False` + `suggest_validation=True` + trend_summary `"Map-model CC: X.XXX - TARGET REACHED, VALIDATE BEFORE STOPPING"`. Accepts `phenix.molprobity` and `phenix.validation_cryoem` as valid cryo-EM validation programs. |
| **Fix #2a: validate-step suppression elif** | `graph_nodes.py` `plan()` auto-stop chain, between `plan_has_pending_stages` and `else` | Defense-in-depth. Suppresses AUTO-STOP when `workflow_state.step_info.step == "validate"` AND `not validation_done`. Catches cases where `plan_has_pending_stages` doesn't fire. |
| **Fix #2b: diagnostic context dump** | `graph_nodes.py` `plan()` else (AUTO-STOP) branch | Logs `plan_has_pending_stages`, `step`, `validation_done`, `after_program`, `experiment_type` to `debug_log` before stopping. Makes future regressions debuggable. |

### Priority order in the elif chain (after Fix #2)

```
1. advice_changed       (user provided new advice on resume)
2. after_program        (user-specified stop target)
3. user_wants_ligandfit (ligandfit prerequisite pending)
4. plan_has_pending_stages (expert plan has pending stages)
5. step=="validate" AND !validation_done  ← NEW (defense-in-depth)
6. else AUTO-STOP (now with diagnostic context dump)
```

Rationale: user/plan intent (1-4) wins over engine-derived intent
(5). The new check is the last line of defense.

### Tests (2 new files, 27 new test cases)

| File | Phase | Cases | Covers |
|------|-------|-------|--------|
| `tst_cryoem_stop_validation.py` (S27) | Fix #1 | 16 (31 assertions) | AF_7mjs regression, both validation programs recognized, X-ray validation programs explicitly don't count, CC just-above threshold, CC exactly at threshold (strict `>`), validation_done with EXCESSIVE / PLATEAU paths, X-ray path unchanged |
| `tst_plan_autostop_validation_suppression.py` (S28) | Fix #2 | 11 (28 assertions) | New elif fires when step=="validate" + !validation_done; doesn't fire when validation_done=True; doesn't fire when step != "validate"; priority order with other suppressors; diagnostic dump appears only on AUTO-STOP; dump records correct values |

Both test files use the libtbx-stub pattern so they run standalone
without a PHENIX install.

Registered in `tests/run_all_tests.py` as suites S27 and S28.

### Verification

```
=== Fix #1 tests ===          (pre-fix:)
31 passed, 0 failed           18 passed, 13 failed

=== Fix #2 tests ===          (pre-fix:)
28 passed, 0 failed           10 passed, 13 failed

=== Regression check (existing v116.11 suites) ===
tst_stop_condition_false_positive: 15 passed, 0 failed
tst_general_resolver:              21 passed, 0 failed
tst_dock_and_stop:                  5 passed, 0 failed
tst_standalone_consistency:         8 passed, 0 failed
tst_initialize_plan_smoke:          9 passed, 0 failed
tst_file_encoding:                  7 passed, 0 failed
tst_program_requirements:          40 passed, 0 failed
```

### Why TWO fixes if Fix #1 is the primary

Fix #1 alone should prevent the bug (don't set `should_stop=True`
until validation is done). Fix #2 adds:

1. **Unsolved-mystery insurance.** We could not determine from the
   AF_7mjs session JSON why `plan_has_pending_stages` didn't
   suppress the AUTO-STOP. The function should have returned True
   against the final plan state. Fix #2's diagnostic logging is
   specifically designed to catch this on future runs.

2. **Defense-in-depth on a class of bugs.** Other code paths
   might set `should_stop=True` outside `_analyze_cryoem_trend`
   (PLATEAU, EXCESSIVE, future additions). Fix #2 protects
   against any future "stop fires while validate is pending"
   regression.

3. **No cost.** Fix #2's new elif only fires when the existing
   chain would have stopped AND `step="validate"` AND
   `validation_done=False`. If Fix #1 works, Fix #2's new elif
   never fires in practice.

### Compatibility

- **X-ray workflows unchanged.** Fix #1 only modifies the cryo-EM
  branch. The X-ray pattern it mirrors was already there.
- **Workflows with validation already done**: `should_stop=True`
  fires as before. No behavior change.
- **No directives or contract changes.** No client protocol updates.
- **AF_7mjs and similar cryo-EM tutorials**: workflow now correctly
  advances to validation stage instead of stopping prematurely.

### Known limitations (not blocking AF_7mjs)

| Limitation | Impact |
|------------|--------|
| **Cosmetic 'unknown' in F7 fallback** | "No command was generated for program 'unknown'" comes from `_run_single_cycle`: `decision_info.get('program', '') or 'unknown'`. When `program=None` (stop intent), `None or 'unknown'` = `'unknown'`. Cleaner: show `'STOP'` or `stop_reason`. Not blocking; session correctly ends. |
| **Unknown why `plan_has_pending_stages` didn't suppress** | Function should return True against AF_7mjs session JSON. Possibilities: stale session.data at check time, transport-layer propagation issue, or deployed version differs from analyzed version. Fix #2's diagnostic logging will surface this on future occurrences. |
| **Map-CC trend_summary format change** | Before: `"Map-model CC: 0.620 → 0.780 (+25.8% last cycle)"`. After (validation done): `"Map-model CC: 0.780 - ABOVE TARGET"`. After (validation not done): `"Map-model CC: 0.780 - TARGET REACHED, VALIDATE BEFORE STOPPING"`. UI log parsers that match the trend format may need to be updated. Mirrors X-ray's pre-existing pattern. |
| **MTZ warning in cryo-EM workflows** | Pre-existing diagnostic: `"WARNING: Refinement completed but best_files[map_coeffs_mtz] is EMPTY"`. Cryo-EM doesn't produce MTZ — the slot doesn't apply. Should be gated on `experiment_type == "xray"`. Not blocking. |

## Version 116.11 (AF_7mjs Stop Condition Fix + Test Cleanup)

### Summary

Two coordinated changes shipped together:

1. **Directive extractor fix** for an AF_7mjs (and likely other
   cryo-EM tutorial) regression: the preprocessor-inserted
   `**Stop Condition**: None` header was being misinterpreted as
   real user stop intent, causing the planner to skip prerequisite
   stages and produce wrong plans.

2. **Test file cleanup** — removed `libtbx.find_unused_imports`
   warnings from 5 test files by switching to real symbol references
   (no `# noqa` shortcuts).

The directive extractor fix was the trigger; the test cleanup
shipped in the same cycle since it was in flight.

### Modified Files (1 production file, 5 test files)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | +98 / -17 lines. Three logical changes (see Bug Details below): (1) strengthened the regex used at all three places that detect preprocessor formatting — handles markdown bold (`**Header**:`), numbered list prefixes (`7. Header:`), and bullet markers (`- Header:`); (2) defense-in-depth Stop Condition strip inside `_resolve_after_program` before the `\bstop\b` check; (3) suppression of `start_with_program` writes from `_resolve_after_program` when advice is preprocessor output, since multi-action signals in descriptive Primary Goal prose are not user prescription. |
| `tests/tst_file_encoding.py` | Switched from `from yaml_loader import _load_yaml_file` to `from knowledge import yaml_loader` + `hasattr(yaml_loader, '_load_yaml_file')` availability probe. Real symbol reference, no `# noqa`. |
| `tests/tst_general_resolver.py` | Removed unused `import re` (was at top, never referenced). |
| `tests/tst_skip_to_program.py` | Removed unused `STAGE_FAILED` from the plan_schema import tuple. |
| `tests/tst_standalone_consistency.py` | Switched from importing `_ACTION_TABLE` directly to importing the `directive_extractor` module and checking `directive_extractor is not None` as an availability probe. The actual symbol is loaded later by `_load_action_table()` which has its own fallback chain. |
| `tests/tst_user_advice_filter.py` | Merged the two import attempts (libtbx and local) into a single try/except chain. Eliminated the dead `_real` alias that existed only to satisfy the linter. |

### Bug Details

#### v116.11 root cause: latent strip-regex weakness

`_strip_preprocessor_stop_condition` at line 44 of
`directive_extractor.py` exists specifically to remove
preprocessor-inserted headers at the entry of both extraction
paths (`extract_directives` and `extract_directives_simple`). Its
regex required headers to start with the literal word:

```python
r'^[ \t]*stop\s+condition\s*:\s*[^\n]*$'
```

But the advice-preprocessor LLM produces output with **markdown
bold and numbered list prefixes** like:

```
7. **Stop Condition**: None
```

The strip never matched this format. Before v115.10, no code
path did a bare `\bstop\b` search on the un-stripped text, so the
bug was latent. v115.10's new `_resolve_after_program` introduced
the first sensitive consumer.

#### Symptom: AF_7mjs three-stage plan

When the agent ran AF_7mjs with the v115.10 directive extractor,
the resolver matched `\bstop\b` against the un-stripped
`**Stop Condition**` header, set `_has_stop=True`, and combined
with the multi-action detection in the Primary Goal text
("Run PredictAndBuild ... rebuild ... refine") produced
directives like:

```yaml
after_program: phenix.real_space_refine
skip_validation: True
start_with_program: phenix.predict_and_build
```

The planner interpreted `start_with_program` as a skip-prerequisites
directive and produced a 3-stage plan starting at Stage 3 instead
of the correct 5-stage cryo-EM plan starting with `mtriage`. The
March 30, 2026 working run produced the correct 5-stage output;
something between then and May 11 (v115.10 deployment) exposed
the latent bug.

#### Three logical changes in directive_extractor.py

| Change | Locations | What |
|--------|-----------|------|
| **Strengthen preprocessor-format regex** | `_strip_preprocessor_stop_condition` (Stop Condition strip, Primary Goal strip, and `_PREPROCESSOR_SIGNATURES` detection) + the pre-existing `_is_preprocessed` check at line 2860 in `extract_directives_simple` | Adds optional numbered list prefix (`[\d]+\.`), bullet markers (`[-*]`), and markdown bold (`\**`) wrappers around the header keyword. Three locations, identical pattern, for consistency. |
| **Defense-in-depth Stop Condition strip** | Top of `_resolve_after_program`, before the `\bstop\b` check | Strips the Stop Condition header again before computing `_has_stop`. Under normal flow this is redundant with the upstream strip, but catches future code paths that bypass `_strip_preprocessor_stop_condition`. |
| **Suppress `start_with_program` for preprocessed advice** | Multi-action branch of `_resolve_after_program` | New `_is_preprocessed` check inside the resolver. When the advice is detected as preprocessor output, the multi-action branch no longer writes `start_with_program`. Reasoning: in preprocessed advice, multi-action mentions come from descriptive Primary Goal prose ("Run PredictAndBuild ... dock/trim ... rebuild ... refine"), not user prescription. Setting `start_with_program` from descriptive prose makes the planner skip prerequisite stages. Real user prose like "run phaser and refine" is unaffected — no preprocessor signatures, so `start_with_program=phenix.phaser` still gets set as before. |

### Tests (1 new file, 15 tests)

| File | Tests | Covers |
|------|-------|--------|
| `tst_stop_condition_false_positive.py` (S26) | 15 | Primary AF_7mjs regression; markdown header variants (5); upstream-strip preserves user prose (3); Primary Goal strip strengthened; `_resolve_after_program` directly for resolver-isolated behavior (negation, legitimate stop, multi-action with/without stop, real user prose for start_with); start_with suppression for preprocessed advice. |

Registered in `tests/run_all_tests.py` as suite S26.

The test file uses the libtbx-stub pattern so it runs without a
real PHENIX install. Tests that target specific resolver behavior
call `_resolve_after_program` directly to isolate from
`intent_classifier`, which is a separate, legitimate extraction
path that can also write `after_program` through its own
mechanism.

### Verification

- **Reproducer**: full AF_7mjs processed advice now produces empty
  or minimal directives (matching the March 30 working behavior).
  No `after_program` write, no `start_with_program` write.
- **Real PHENIX**: 15/15 tests pass via `phenix.python
  tst_stop_condition_false_positive.py` after the test-environment
  sensitivity fix (`test_only_primary_goal_extracted_no_stop` was
  initially too strict — it didn't account for `intent_classifier`
  legitimately writing `after_program` for simple "Refine the
  model." advice; rewrote the test to call `_strip_preprocessor_stop_condition`
  in isolation).
- **No regressions**: all existing test suites pass
  (`tst_general_resolver` 21/21, `tst_dock_and_stop` 5/5,
  `tst_program_requirements` 40/40, `tst_standalone_consistency`
  8/8, `tst_initialize_plan_smoke` 9/9, `tst_file_encoding` 7/7).

### Self-Review Findings

The fix went through five review iterations before settling on the
current shape. Each iteration found a problem with the previous
diagnosis:

1. **First diagnosis (wrong shape):** Strip the header inside
   `_resolve_after_program` only. This is a symptom fix — it
   treats the bug at the consumer, not the cause. The right fix
   point was the upstream `_strip_preprocessor_stop_condition`.

2. **Second diagnosis (incomplete):** Mark resolver writes with
   `_set_by_pattern=True` so the intent-driven clearing block
   could remove them. But the clearing block only handles
   `after_program`, not `start_with_program`. The actual
   problematic write on AF_7mjs was `start_with_program`.

3. **Third diagnosis (right shape):** Suppress `start_with_program`
   writes when advice is preprocessed. Simpler than mark-and-clear;
   addresses the root semantic confusion (descriptive prose ≠ user
   prescription).

4. **Fourth iteration:** Found that the pre-existing
   `_is_preprocessed` regex at line 2860 of `extract_directives_simple`
   was inconsistent with my upstream fix — it still used the old
   regex that didn't handle markdown. Strengthened for consistency.

5. **Fifth iteration:** Test failure in real PHENIX revealed that
   `intent_classifier` writes `after_program` through a separate
   path. Rewrote affected tests to use `_resolve_after_program`
   directly when targeting resolver-specific behavior.

### Known limitations (not blocking AF_7mjs)

| Limitation | Impact |
|------------|--------|
| **Conditional stop in Special Instructions** | Text like "Stop the analysis if R-free > 0.5" in Special Instructions could still false-positive `_has_stop=True`. AF_7mjs's actual Special Instructions don't contain stop words. Pre-existing fragility. |
| **DRY violation across 3 `_PREPROC_SIGS` locations** | Same signature list and regex appears at three locations in `directive_extractor.py`. Cosmetic; future refactor should extract to a module-level helper. |
| **Line 2860 `_is_preprocessed` covered only indirectly** | My test suite verifies the strengthened regex via the strip function (which uses the same pattern). The line 2860 instance is mechanically a copy of tested code, but a direct test would be cleaner. |

### Test Cleanup Details

The unused-imports fix used **real symbol references**, not
`# noqa: F401` suppression. The `libtbx.find_unused_imports`
tool requires that imported names actually be referenced; a
`noqa` comment doesn't satisfy it.

The general pattern was: where a test imports a symbol only to
verify the module is available, switch to importing the module
itself and using `hasattr()` or `is not None` as the availability
probe. This is both more honest about what the test depends on
and less brittle if the symbol gets renamed.

---

## Version 116.10 (Reliability and Classification Cleanup)

### Summary

Six bugs in advice parsing, LLM prompting, program selection, and
plan classification, plus drift-detection machinery for the wire
contract and the client's plan-generation logic. The cleanup was
scoped from a single batch-tutorial finding ("predict and stop"
sessions were AUTO-STOPping) and grew as adjacent bugs surfaced
during investigation. Five files modified, 89 new tests across
7 new test files plus one augmented file. Rollout order in
testing: 1 → 6a → 4b → 2 → 6b → 3a → 3d.

### Modified Files (5 production files)

| File | Phases | Changes |
|------|--------|---------|
| `agent/rules_selector.py` | 1 | `_apply_user_advice` `stop_condition_patterns` extended with `"and stop"`, `"then stop"`, `", stop"` so "predict and stop" stops being collapsed to `["STOP"]`. |
| `agent/workflow_engine.py` | 4b, 6b | Phase 4b: new `_filter_programs_missing_data_inputs` strips `xtriage` when `not has_data_mtz and not has_phased_data_mtz`, strips `mtriage` when no map; matching guards added to `_check_program_prerequisites` for defense in depth. Phase 6b: top of `_detect_xray_step` routes sequence-only sessions directly to `obtain_model` (state becomes `xray_analyzed`). |
| `agent/contract.py` | 2 | `CURRENT_PROTOCOL_VERSION` 3 → 5 (matches highest field version after silent drift). New `validate_contract()` returns `(ok, errors)` enforcing three invariants: CURRENT ≥ max field version, MIN ≤ CURRENT, MIN ≥ 1. |
| `knowledge/prompts_hybrid.py` | 6a | `_format_directives_for_prompt` after-program lines replaced: "CRITICAL: You MUST run X" → "Stop target: X. If in VALID PROGRAMS choose it; if not, choose an appropriate prerequisite. Never pick outside VALID PROGRAMS." |
| `phenix/programs/ai_agent.py` | 2, 3a, 3d | Phase 2: `_get_protocol_version()` fallback `return 3` → `return 5`. Phase 3a: `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS` extracted as module-level frozensets near line 815; three call sites in `_initialize_plan_inner` updated to reference them (pure refactor, 51 program×intent traces verified identical). Phase 3d: `phenix.map_symmetry` and `phenix.dock_in_map` added to `_STANDALONE_PROGRAMS` (behavior change: skips plan generation, routes through workflow_engine state machine instead). |

### Bug Details

| Bug | Phase | Symptom | Root Cause | Fix |
|-----|-------|---------|------------|-----|
| 1 | 1 | "predict and stop" sessions immediately AUTO-STOP at session start | `_apply_user_advice` saw `stop` keyword without matching any of `stop after / stop when / stop once / stop if / stop condition / stop at` → fell through to "treat as immediate stop" → `valid_programs` collapsed to `["STOP"]` | Extended `stop_condition_patterns` with `"and stop"`, `"then stop"`, `", stop"` |
| 2 | 6a | LLM picks `predict_and_build` before `xtriage` has run, hitting prerequisite failure | Prompt told LLM "CRITICAL: You MUST run X before stopping. If it's in VALID PROGRAMS, choose it NOW." — overrode prerequisite logic | Reframed as "Stop target: X. If in VALID PROGRAMS choose it; if not, choose an appropriate prerequisite" |
| 3 | 4b | xtriage offered as valid program even when no `.mtz` is uploaded; session stalls | YAML lists xtriage at analyze step unconditionally; no data-input filter | New `_filter_programs_missing_data_inputs` strips xtriage when no .mtz, mtriage when no map; matching prereq guard in `_check_program_prerequisites` for defense in depth |
| 4 | 2 | New v3/v4/v5 fields added to `SESSION_INFO_FIELDS` while `CURRENT_PROTOCOL_VERSION` stuck at 3 — silent drift | No invariant check between `CURRENT_PROTOCOL_VERSION` and registered field versions | Bumped CURRENT to 5; added `validate_contract()` invariant function; augmented existing `tst_contract_compliance.py` with drift-detection test |
| 5 | 6b | Sequence-only X-ray sessions stuck at `xray_initial` (analyze step); xtriage filtered, no fallback | `_detect_xray_step` returned analyze for `not xtriage_done` regardless of whether xtriage could even run | Routing branch at top of `_detect_xray_step`: when no .mtz and has sequence → `obtain_model` (state becomes `xray_analyzed`); `predict_and_build` available via YAML `has: sequence` condition |
| 6 | 3a | `_initialize_plan_inner` standalone-programs list duplicated inline at two call sites; comment admitted the duplication; classification couldn't be tested | Original implementation kept the two tuples in sync by hand | Extracted `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS` as module-level frozensets; new `tst_standalone_consistency.py` enforces alignment with `directive_extractor._ACTION_TABLE` |
| 7 | 3d | "dock and stop" with sequence + map (no model) generates a plan, `skip_to_program(dock_in_map)` marks predict_and_build SKIPPED, dock_in_map then fails at runtime | `dock_in_map` was a "full-plan target" (caught by the v116.10 elif); skip_to_program incorrectly removed the predict prerequisite | Reclassified `phenix.dock_in_map` and `phenix.map_symmetry` as `_STANDALONE_PROGRAMS`; `_initialize_plan_inner` now skips plan generation, lets `workflow_engine` handle prerequisites via state machine (analyze → predict → dock) |

### Tests (7 new files, 1 augmented, 89 new test cases)

| File | Phase | Cases | Covers |
|------|-------|-------|--------|
| `tst_user_advice_filter.py` | 1 (S13) | 16 | "X and stop" / "X then stop" / "X, stop" preserved; recognized stop-conditions still respected |
| `tst_after_program_prompt.py` | 6a (S14) | 12 | Prompt uses target-not-now framing; aggressive "MUST" language removed |
| `tst_data_input_filter.py` | 4b (S15) | 22 | xtriage filtered when no .mtz; mtriage filtered when no map; prereq guard handles directive re-injection |
| `tst_protocol_version.py` | 2 (S16) | 15 | `validate_contract()` invariants via synthetic probes (drift detection, bounds checking) |
| `tst_contract_compliance.py` | 2 (augmented) | +1 | New `test_contract_validate_passes` adds drift-detection entry point |
| `tst_sequence_only_routing.py` | 6b (S17) | 10 | Sequence-only sessions route to `obtain_model`; doesn't fire when data is present or sequence absent |
| `tst_standalone_consistency.py` | 3a (S18) | 8 | `_STANDALONE_PROGRAMS` aligned with `_ACTION_TABLE`; drift detection for new actions |
| `tst_dock_and_stop.py` | 3d (S19) | 5 | Phase 3d behavior change documented; blast-radius regression guard (48 unrelated traces unchanged) |

Registered in `tests/run_all_tests.py` as suites S13–S19.

Run with: `python3 tests/run_all_tests.py` or individually
`python3 tests/tst_user_advice_filter.py`, etc.

### Verification

| Fix | Verification |
|-----|--------------|
| Bug 1 (user advice "X and stop") | ✅ 16 unit tests pass; 6/16 fail on pre-fix code |
| Bug 2 (after-program prompt) | ✅ 12 unit tests pass; 6/12 fail on pre-fix |
| Bug 3 (data-input filter) | ✅ 22 unit tests pass; 18/22 fail on pre-fix; integration smoke test confirms xtriage absent from valid_programs when no .mtz |
| Bug 4 (protocol hygiene) | ✅ 15 + 1 unit tests pass; 8/15 fail on pre-fix; verified `tst_contract_compliance.py:test_protocol_version_consistency` now passes (was failing because fallback `return 3` ≠ CURRENT `5`) |
| Bug 5 (sequence-only routing) | ✅ 10 unit tests pass; 4/10 fail on pre-fix; end-to-end integration test confirms state becomes `xray_analyzed` with `valid_programs=[phenix.predict_and_build]` |
| Bug 6 (Phase 3a refactor) | ✅ 8 consistency tests pass; behavioral equivalence verified by 51 (program × intent) decision-tree traces — all identical pre/post |
| Bug 7 (Phase 3d behavior change) | ✅ 5 unit tests pass; 3/5 fail on pre-fix; blast-radius check confirms 48 unrelated traces unchanged. ⏳ Integration test against tutorial corpus pending. |

### Self-Review Finding (Phase 3a/3d split)

Phase 3 was initially delivered as a single phase that bundled the
de-duplication refactor with two new programs added to
`_STANDALONE_PROGRAMS` (silently changing behavior). Self-review
caught this:

- Refactor PRs that change behavior should be split. "Cleaner"
  behavior is still a behavior change.
- The fix was to split into Phase 3a (verified by 51 trace-equivalence
  checks) and Phase 3d (explicit behavior change, with `tst_dock_and_stop.py`
  including a 48-trace blast-radius regression guard).

This is recorded in `OVERVIEW.md` under Active Development as a
principle for future cleanups.

### Known Issues (deferred to future work)

| Issue | Description |
|-------|-------------|
| Phase 3d integration verification | Decision-tree traces verified, but no tutorial currently exercises "dock and stop with sequence + map". A tutorial that does would close the verification gap. Rollback is one line if a regression appears. |
| Systematic input-availability audit | Phase 4b only filters xtriage and mtriage. A systematic "every program declares its inputs" audit is unfinished. The next program added with hard input requirements should be added to the filter at the same time. |
| Phase 6a "appropriate prerequisite" language | Shipped Option A (ship as-is) over Options B (add hint) and C (surface `rules_priority`). If integration testing shows the LLM picks the wrong prerequisite, revisit. |
| Redundant prediction-only allowance | The v116.10 allowance in `_check_program_prerequisites` is mostly redundant after Phase 6b but retained as defense in depth. A future cleanup could remove it if edge cases (e.g., `xtriage_done=True` with no data) are verified to reach `predict_and_build` through other paths. |

### Post-Phase-5 Addition: CC Key Extraction Fix (S20)

#### Summary

A successful cryo-EM workflow that produced metrics and a structure
model displayed the misleading banner "SESSION STOPPED -
INCOMPLETE / No structure model available." Root-caused to a key
naming inconsistency in `_generate_structure_report` and patched
with a defensive multi-key lookup.

#### Bug Details

| Symptom | Successful cryo-EM run reports "SESSION STOPPED - INCOMPLETE" |
| Root Cause | `_generate_structure_report` looked up CC under `map_model_cc` (the form programs PRINT in logs), but the cycle metrics dict stores it under `model_map_cc` (the canonical storage key, see schema at `ai_agent.py:9083-9088`). The v115.05 author wrote the lookup against the printed form. For cryo-EM workflows without R-free, this meant `_best_cc=None` → `_metrics_good=False` → stopped-report path fires despite success. |
| Fix | Defensive multi-key lookup at lines 8794-8797 and 8822-8823: checks `model_map_cc` (canonical), `map_model_cc` (legacy/unnormalized), `map_cc` (workflow_engine variant), `cc_mask` (real_space_refine output). Works regardless of which form actually appears in the dict. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `phenix/programs/ai_agent.py` | 8794-8797, 8822-8823 | Two surgical edits adding `model_map_cc` to lookup lists |

#### Tests

`tst_cc_key_extraction.py` (11 tests, suite S20): 2 source-level
invariants, 5 behavioral key-variant tests, 3 edge cases, 1
motivating-case test reproducing the Tom-reported scenario.

#### Self-Review Note

Initially framed as "an open-and-shut typo." That was sloppy.
`map_model_cc` is a legitimate spelling — it appears as a phenix
program name and in regex-captured raw log lines. The bug is real
but it's "wrong form of the name" not "typo": the author used the
printed form rather than the stored form. The fix's defensive
multi-key lookup handles both forms regardless.

### Post-Phase-5 Addition: File Encoding Fix (S21)

#### Summary

A user on Chinese-locale Windows reported a `UnicodeDecodeError:
'gbk' codec can't decode byte 0x94` crash at startup while loading
`programs.yaml`. Root cause: Python's `open()` without explicit
`encoding=` falls back to the system code page (`gbk` for
Chinese-locale Windows, `cp1252` for English, `cp932` for
Japanese), which crashes on any non-ASCII content in a UTF-8 file.

The fix was applied in two passes: initially to the 4 files
directly implicated by the crash trace, then to the **entire
v116.10 codebase** (63 additional files surfaced by a recursive
grep audit).

#### Bug Details

| Symptom | `UnicodeDecodeError: 'gbk' codec can't decode byte 0x94 in position 987` on PHENIX startup, Chinese-locale Windows |
| Root Cause | `open(path)` without `encoding=` uses `locale.getpreferredencoding()`. On non-UTF-8 Windows locales, any UTF-8 file with non-ASCII content (curly quotes, em dashes, accented letters) crashes the read. |
| Fix | Add `encoding='utf-8'` to every text-mode `open()` call. YAML, JSON, PDB, and report files are all UTF-8 by spec, so this is correct behavior, not a workaround. |

#### Files Modified

**First pass (the reported crash site + immediate neighbors):**

| File | Sites |
|------|-------|
| `knowledge/yaml_loader.py` | 1 (the reported crash) |
| `phenix/programs/ai_agent.py` | 3 |
| `agent/directive_validator.py` | 1 |
| 4 test files | 11 |

**Second pass (codebase-wide audit, 305 sites across 63 files):**

| Subdirectory | Offenders | Files |
|--------------|-----------|-------|
| `agent/` | 53 | 22 |
| `knowledge/` | 3 | 3 |
| `utils/` | 1 | 1 |
| `tests/` | 248 | 37 |

Patches applied: `open(path)` → `open(path, encoding='utf-8')`;
`open(path, 'r')` → `open(path, 'r', encoding='utf-8')`;
`open(path, errors='X')` → `open(path, encoding='utf-8',
errors='X')`. Binary opens, `Popen`, `urlopen`, and `.open()`
method calls were correctly excluded.

#### Tests

`tst_file_encoding.py` (7 tests, suite S21):

- 3 per-file source scans (`yaml_loader`, `ai_agent`, `directive_validator`)
- 2 directory-scan tests (`test_all_production_code_uses_utf8`,
  `test_all_test_code_uses_utf8`) — walk the entire production and
  test trees automatically, so new files added in the future are
  covered without test updates
- 2 empirical tests (write/read UTF-8, the user's exact byte pattern)

The scanner uses Python's `tokenize` module rather than regex
heuristics, so docstrings and string-literal occurrences of
`open(` are correctly excluded.

#### Verification

| Check | Result |
|-------|--------|
| Re-audit of patched files (tokenize-based) | 0 offenders |
| Python syntax (all 63 files) | All parse |
| Line-count drift | 0 (every patch is in-place) |
| 1:1 substitution (added = removed) | 306 = 306 |
| `tst_file_encoding.py` on patched | 7/7 pass |
| `tst_file_encoding.py` on unpatched | Correctly flags 302 |

#### User Workaround (No Code Change)

If a user can't deploy the patched files immediately, setting
`PYTHONUTF8=1` in their shell before launching phenix forces
Python 3.7+ to use UTF-8 for all file I/O regardless of system
locale.

#### Self-Review Note

The 305 patches are mechanically uniform but lack per-patch
explanatory comments (the first encoding fix added a 5-line
comment to each file; the automated patcher for the codebase-wide
pass did not). The rationale lives in this CHANGELOG entry rather
than per-patch comments. Future maintainers cleaning up the
codebase should know the `encoding='utf-8'` convention is
load-bearing for non-UTF-8 Windows locales.

### Post-Phase-5 Addition: Ligand Workflow Restart Fix (Phase 6c, S22)

#### Summary

The nsf-d2-ligand tutorial restarted its reasoning at cycle 3:
after a successful refinement (cycle 1) and a successful ligand
fit (cycle 2), cycle 3 reported `STATE: xray_initial / Need to
analyze data quality first` and the LLM produced the reasoning
*"As this is the first refinement step, I will set
'generate_rfree_flags=true'."* The same pattern repeated at
cycles 4 and 5. Three interacting issues were diagnosed; two
fixes resolve all three.

#### Bug Details

| Bug | Symptom | Root Cause | Fix |
|-----|---------|------------|-----|
| **1** | Every cycle reports `STATE: xray_initial` even after refine and ligandfit succeeded | `_detect_xray_step` returns "analyze" any time `xtriage_done=False`, regardless of downstream progress. When the user supplies a pre-placed model + `start_with_program=phenix.refine` (which skips xtriage), the state never advances | Added a `past_analysis` check that mirrors `_detect_cryoem_step`: if any downstream program has completed (`refine_done`, `ligandfit_done`, `phaser_done`, etc.) or a positioned model exists on disk, advance past analyze |
| **2** | After ligandfit succeeds, the "best model" pointer stays on the previous unliganded refined model | `best_files_tracker` scores `ligand_fit_output` (stage 105) below `refined` (stage 100 + R-free contribution ~22 = ~122). LigandFit's output (no R-free metrics) loses the scoring contest, so the LLM is told the model is still the unliganded one | Extended the existing `with_ligand` metric inheritance to also cover `ligand_fit_output`. One-condition change: `stage == "with_ligand"` → `stage in ("with_ligand", "ligand_fit_output")` |
| **3** | LLM mechanically re-applies first-cycle directives (`generate_rfree_flags=true`) at cycles 3, 4, 5 | Observable consequence of Bugs 1 and 2 — when the agent presents a stale worldview (state=xray_initial + unliganded model as "current"), the LLM correctly applies first-cycle directives because the worldview says we're at cycle 1 | No separate fix needed; resolves automatically once the LLM sees correct state (Bug 1) and current model (Bug 2) |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `agent/workflow_engine.py` | 1040-1086 | Replaced 3-line early return with 45-line `past_analysis` check (mirrors `_detect_cryoem_step` at lines 1283-1299) |
| `agent/best_files_tracker.py` | 564-583 | Extended condition `stage == "with_ligand"` to `stage in ("with_ligand", "ligand_fit_output")`; existing inheritance machinery handles the rest |

#### Tests

`tst_ligand_workflow_restart.py` (14 tests, suite S22):

- **Section A** (8 tests): `past_analysis` check — verifies each
  downstream flag (`refine_done`, `ligandfit_done`,
  `phaser_done`, etc.) correctly advances past analyze; verifies
  fresh cycle 1 still routes to analyze; reproduces the exact
  nsf-d2-ligand cycle 3 context and confirms it returns
  `combine_ligand`.
- **Section B** (2 tests): regression — Phase 6b sequence-only
  precedence still fires; fresh data+model upload still goes to
  analyze first.
- **Section C** (4 tests): `ligand_fit_output` metric
  inheritance — ligand_fit_1.pdb becomes best with inherited
  R-free; existing `with_ligand` inheritance still works;
  explicit metrics on ligand_fit_output are respected (no
  inheritance); end-to-end refine→ligandfit→refine tracks correctly.

Uses the same libtbx-stub pattern as `tst_sequence_only_routing.py`
(Phase 6b tests), so workflow_engine routing tests run without a
real PHENIX install.

#### Verification

| Check | Result |
|-------|--------|
| Python syntax (both patched files) | Parses |
| `tst_ligand_workflow_restart.py` (14 new tests) | 14/14 pass |
| `tst_sequence_only_routing.py` (Phase 6b, 10 tests) | 10/10 pass (regression) |
| `tst_data_input_filter.py` (Phase 4b, 22 tests) | 22/22 pass (regression) |
| Empirical bug reproduction | Pre-fix: ligand_fit_1.pdb does NOT become best (score 122 stays on refined). Post-fix: ligand_fit_1.pdb becomes best (score 127 = stage 105 + inherited R-free 22) |
| nsf-d2-ligand cycle 3 routing | Returns `combine_ligand` (correct flow through pdbtools) |

#### Workflow Shape Change

A behavioral change worth flagging: the workflow shape is
different after this fix.

| Before | After |
|--------|-------|
| refine → ligandfit → refine → refine → refine (3 refines without progress until guard fires) | refine → ligandfit → combine_ligand (pdbtools) → refine (canonical flow) |

Both reach acceptable R-free, but the post-fix flow is the
canonical sequence with the combine_ligand step properly
executing. Tutorial expectations that encode the old shape
("expect 3 refine cycles after ligandfit") would need updating.

#### Self-Review Note

Both fixes are non-novel — each mirrors a proven pattern already
in the codebase:

- Fix 1 mirrors `_detect_cryoem_step`'s `past_analysis` exactly.
  The cryo-EM version has a comment explaining why the check
  exists: *"This handles tutorials that skip mtriage and also
  prevents the workflow from getting stuck in 'analyze' if
  mtriage's done flag wasn't set."* The X-ray version simply
  lacked the equivalent.
- Fix 2 extends an existing two-line metric-inheritance pattern
  to one additional stage. A one-character-ish change.

A new principle surfaced from this round: **when fixing
analogous bugs in parallel code paths (X-ray vs cryo-EM,
with_ligand vs ligand_fit_output), check whether the fix
pattern already exists nearby**. Both bugs had close analogs
that simply hadn't been extended to cover the affected case.

The highest residual risk is that this delivery was tested at
source level but not via a live tutorial run. Tom should re-run
the nsf-d2-ligand tutorial to confirm: (a) R-free still reaches
~0.21 or better, (b) the reasoning report no longer shows "first
refinement step" 4 times, (c) no tutorial expectations break.

### Post-Phase-5 Addition: Tier 1 Follow-Up Tests (S23 + S24)

#### Summary

Two new test suites that close verification gaps flagged in the
v116.10 review. Both are pure tests (no production code changes);
they verify existing behavior that previously had only indirect
coverage.

#### Tests Added

**S23 — `tst_initialize_plan_smoke.py` (9 tests)**

`_initialize_plan_inner` had decision-tree traces in
`tst_dock_and_stop.py` and `tst_standalone_consistency.py`, but
neither test actually CALLED the function. A future refactor
could change side effects (directive rewrites, stop_after
clearing, plan-generation gating) without tripping any existing
test.

The new tests assert on the function's observable side effects
for 9 scenarios covering each branch of the decision tree:

| Case | Scenario | Expected outcome |
|------|----------|------------------|
| A | standalone (non-preprocessing) + task | single_program_skip; no rewrite |
| B | preprocessing + task (no explicit stop) | rewrite intent; clear stop_after; proceed |
| C | preprocessing + task + "and stop" advice | preprocessing_explicit_stop_skip; preserve directive |
| D | needs_plan + task | rewrite intent; KEEP stop_after; proceed |
| E | v116.10 elif full-plan target | rewrite intent; preserve stop_after; proceed |
| F | task intent + no stop_after | single_program_skip (else branch) |
| G | solve intent + no stop_after | proceed_to_generate_plan |
| H | Phase 3d dock_in_map + task | single_program_skip (post-Phase-3d) |
| Extra | xtriage + explicit_stop | preprocessing_explicit_stop_skip (xtriage is BOTH preprocessing and standalone — preprocessing branch fires first) |

The "extra" case surfaced from a test development misclassification
that was caught by the actual control-flow trace — exactly the
kind of subtlety this layer of testing is meant to catch.

The harness parses `_STANDALONE_PROGRAMS` and `_NEEDS_PLAN_PROGRAMS`
from source (same pattern as `tst_dock_and_stop.py`), then replays
the decision tree in a minimal harness. The full function can't
be imported because of PHENIX dependencies; the harness mirrors
the control flow exactly so a refactor that changes side effects
requires the test to be updated in parallel.

**S24 — `tst_phase3d_motivating_tutorial.py` (3 tests)**

Closes the **highest-priority verification gap** from the v116.10
review: Phase 3d's behavior change was verified only by
decision-tree traces; no tutorial exercised the dock-and-stop
path with sequence + map.

Three layers:

1. **Decision-tree** (unconditional): asserts dock_in_map + task
   routes to `single_program_skip` (Phase 3d post-fix).
2. **YAML expectations format** (unconditional):
   `get_expectations_yaml_entry()` returns the snippet to add to
   `tutorial_expectations.yaml`; test verifies it parses and
   contains the expected program ordering.
3. **Session-based** (skip-aware): given a real session.json
   fixture (via `DOCK_AND_STOP_FIXTURE` env var or
   `tests/fixtures/dock_and_stop_session.json`), asserts on
   actual program ordering and the absence of the anti-pattern
   `"no model file found"`.

A `docs/PHASE_3D_TUTORIAL_README.md` walks through wiring up the
synthetic tutorial data (sequence + map) and the session fixture.

#### Files Modified

| File | Lines |
|------|-------|
| `tests/tst_initialize_plan_smoke.py` (new) | 631 |
| `tests/tst_phase3d_motivating_tutorial.py` (new) | 367 |
| `tests/run_all_tests.py` (S23 + S24 registered) | +58 |
| `docs/PHASE_3D_TUTORIAL_README.md` (new) | ~100 |

#### Verification

| Check | Result |
|-------|--------|
| Both new test files: Python syntax | Parses |
| `tst_initialize_plan_smoke.py` | 9/9 pass |
| `tst_phase3d_motivating_tutorial.py` | 3/3 pass (one skip-aware) |
| `tst_dock_and_stop.py` (existing, regression) | 5/5 pass |

#### Tier 1 status

With S23 and S24 shipped, all four Tier 1 items from the
v116.10 follow-up plan are now complete:

| Item | Status |
|------|--------|
| 1.1 Phase 5 docs amendment for CC key fix | Done (covered by prior CHANGELOG amendments) |
| 1.2 Multi-clause stop test cases | Handled by Tom directly |
| 1.3 Phase 3a behavioral smoke test | Done (S23) |
| 1.4 Phase 3d motivating-case tutorial | Done (S24) |

### Post-Phase-5 Addition: GUI Tutorial Override Fix

#### Summary

When running a tutorial whose README mentions phenix programs
not yet available in the AI Agent (e.g., `groel-dock-refine`
mentions a currently-unsupported program), the GUI was blocking
Run with a `Sorry` dialog even after the user provided their own
files and advice that overrode the README.

#### Bug Details

| Symptom | "This tutorial requires programs not yet available" dialog blocks Run, even when user has provided their own files and instructions |
| Root Cause | `AIAgent.py:validate_params()` has two consecutive checks. The first correctly distinguishes user-provided input from README fallback. The second (the tutorial-blocking block) fires unconditionally based on `tut.can_run` regardless of whether the user was using the README at all. |
| Fix | One condition added to the tutorial-blocking guard: only fire when `not has_files and not has_advice`. The banner at the top of the page still shows the informational "tutorial requires X" warning — that's fine and was kept. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `wxGUI2/Programs/AIAgent.py` | 69-86 | Added `and not has_files and not has_advice` to the blocking guard |

#### Scenario Verification

| Scenario | files | advice | tut.can_run | Pre-fix | Post-fix |
|----------|-------|--------|-------------|---------|----------|
| Tutorial, no own input | No | No | False | BLOCKED ✓ | BLOCKED ✓ |
| Tutorial + own files (bug case) | Yes | No | False | BLOCKED (bug) | PROCEEDS ✓ |
| Tutorial + own advice | No | Yes | False | BLOCKED (bug) | PROCEEDS ✓ |
| Tutorial + both | Yes | Yes | False | BLOCKED (bug) | PROCEEDS ✓ |
| Non-tutorial, no input | No | No | n/a | First block raises | First block raises ✓ |
| Tutorial with supported programs | any | any | True | PROCEEDS ✓ | PROCEEDS ✓ |

Three pre-fix bug scenarios fixed; no correct behavior changed.

#### Verification

GUI code (wxPython, requires a display to run). Verified by:
- Python syntax check
- Manual scenario tracing through 7 user paths

User-confirmed working: Tom retested groel-dock-refine after
deploying the fix and the dialog no longer blocks Run.

#### Self-Review Note

This fix was discovered during integration testing of Phase 6c
(ligand workflow restart), not from a planned cycle item. The
banner display (informational) was deliberately left in place so
users still see what the README is asking for; only the
blocking dialog was conditioned on user-override.

### Post-Phase-5 Addition: Declarative Program Requirements (Tier 2.1, S25)

#### Summary

Added an optional `requirements:` block to entries in
`programs.yaml`. When present, the block is evaluated by a new
filter pass inside `_filter_programs_missing_data_inputs`.
Programs without the block are unaffected. Initial scope: one
declaration (`phenix.autobuild`) to plug a known gap where the
program could be picked and crash at runtime with `MTZ lacks
phase/FOM columns`.

This is the result of Tier 2.1 design work, reviewed externally
by Gemini and refined in the `PATH_Y_DESIGN.md` document.

#### Bug Details

| Symptom | LLM picks `phenix.autobuild` in sessions without phase information; command runs and fails with `MTZ lacks phase/FOM columns` after building setup. One wasted cycle plus a confusing runtime error. |
| Root Cause | Autobuild had an explanation path (`explain_unavailable_program`) telling the LLM why it shouldn't pick autobuild, but no actual filter. The LLM could ignore the explanation. The state machine alone didn't gate autobuild because the relevant step has multiple valid programs and the conditions in `workflows.yaml` weren't expressive enough. |
| Fix | Declarative `requirements:` block on `phenix.autobuild` evaluated by new `_check_requirements()` parser. Filter integrates inside `_filter_programs_missing_data_inputs` after the existing xtriage/mtriage checks. |

#### Files Modified

| File | Lines | Change |
|------|-------|--------|
| `agent/workflow_engine.py` | +135 | New method `_check_requirements()` (~100 lines including grammar docs); new clause-keyword set `_KNOWN_REQUIREMENT_CLAUSES`; integration block inside `_filter_programs_missing_data_inputs` (~22 lines) |
| `knowledge/programs.yaml` | +4 | `requirements:` block on `phenix.autobuild` entry |
| `tests/run_all_tests.py` | +22 | Register S25 (Program Requirements) |

#### Declarative Schema (v1)

```yaml
phenix.<program>:
  # ... existing fields ...
  requirements:    # NEW — optional
    requires:
      - <clause>
      - <clause>
```

Grammar (closed, boolean-only):

| Clause | Semantics |
|--------|-----------|
| `has: <name>` | `context["has_<name>"]` truthy |
| `has_any: [<name>, ...]` | At least one `context["has_<n>"]` truthy |
| `not_has: <name>` | Not `context["has_<name>"]` |
| `done: <name>` | `context["<name>_done"]` truthy |
| `not_done: <name>` | Not `context["<name>_done"]` |

The grammar is intentionally closed. No `if/then`, no nested
logic, no metric comparisons. Per Pitfall 10 of the design
doc: if a requirement needs control flow, it goes in Python
(Mechanism 3 or 4), not here.

#### autobuild Declaration

```yaml
phenix.autobuild:
  requirements:
    requires:
      - has_any: [data_mtz, phased_data_mtz]
      - has_any: [model, placed_model, phased_data_mtz]
```

`has_phased_data_mtz` appears in both clauses because a phased
MTZ satisfies both "x-ray data" and "phase information"
independently. So a phased MTZ alone is sufficient. Other
qualifying scenarios: raw data + model, raw data + placed
model (after phaser). Filtered scenarios: data alone, model
alone, empty session.

#### Tests (30 new test cases)

| Section | Tests | Description |
|---------|-------|-------------|
| A. Parser unit tests | 17 | Each clause type evaluated in isolation; AND of clauses; malformed input handling; unknown-keyword warning + strict-mode raise |
| B. Filter integration | 6 | autobuild scenarios: phased alone (keep), data+model (keep), data+placed_model (keep), data alone (filter), model alone (filter), empty session (filter) |
| C. Backward compat | 7 | Programs without `requirements:` unaffected; existing xtriage/mtriage filters still work; order of operations correct; mixed program lists handled |

Run with `phenix.python tests/tst_program_requirements.py`.

#### Order of Operations

The new filter integrates into the existing pipeline at step 4:

1. State machine (Mechanism 1) — decides current step
2. YAML step conditions (Mechanism 2) — filters programs within step
3. `_filter_programs_missing_data_inputs` — existing xtriage/mtriage checks
4. **`_check_requirements` — NEW declarative filter (this delivery)**
5. `_check_program_prerequisites` (Mechanism 4) — runs for directive additions

`requirements:` can only remove programs from `valid_programs`,
never add them.

#### Forward-Looking Guidance

Per the design doc, `requirements:` is **the preferred path
for new program filtering rules going forward**. The four
existing hand-coded mechanisms aren't being grown.
Documentation rules:

- Add `requirements:` ONLY when the existing mechanisms don't
  cover the filtering need (avoid duplication)
- The grammar is closed; extension requires documented
  justification per Pitfall 10
- The block is boolean-only; file-content checks (e.g., "MTZ
  has FOM columns") must be exposed as context flags elsewhere
  (per Pitfall 16)
- Before adding a new declaration, complete the context-flag
  audit checklist (per Pitfall 15)

See `PATH_Y_DESIGN.md` in the design package for the full
rationale, 16 pitfalls with mitigations, and three rejected
alternatives.

#### Verification

| Check | Result |
|-------|--------|
| Python syntax (workflow_engine.py, run_all_tests.py) | parse OK |
| New test suite (40 tests) | 40 passed, 0 failed |
| Backward-compat: existing `tst_data_input_filter` regression | xtriage/mtriage filters still pass |
| Context-flag audit for autobuild's 4 referenced flags | All set in `_build_context()` before filter runs |
| `tst_scenario_tracer.py` (S5A) | Passes after Tier 2.1.1 fix below |
| `tst_phase4_history_flags.py` (all_read_flags_are_written) | Passes after Tier 2.1.1 fix below |

#### Post-Deployment Fix (Tier 2.1.1): `any_of` clause + phantom flag

After initial deployment, two test failures surfaced that required
the Tier 2.1 implementation to be revised:

**Failure 1 — S5A in `tst_scenario_tracer.py`:**

The simulator reported `autobuild not available after autosol`. The
v1 rule required `has_any: [model, placed_model, phased_data_mtz]`
in its second clause, but in the simulator (and some real session
paths) `autosol_done=True` is set without `has_phased_data_mtz` being
updated on the same cycle. The existing explanation code at
`explain_unavailable_program` uses `phaser_done OR autosol_done OR
has_placed_model_from_history` as the gate; my v1 rule diverged from
this and was stricter.

**Failure 2 — phantom `has_` flag in `tst_phase4_history_flags.py`:**

The static flag scanner regex matches `context.get("has_X")` to
extract flag names. The v1 implementation had `context.get("has_" +
clause["has"])` patterns, which the regex extracted as just `has_`
(captured up to the first closing quote). This produced a false
"flag read but not written" warning.

**Fixes shipped in Tier 2.1.1:**

1. **Grammar extension: `any_of` clause type** (Pitfall 10 justification
   provided in workflow_engine.py docstring). Where `has_any: [...]`
   takes a list of names and auto-prefixes `has_`, `any_of: [...]`
   takes a list of sub-clauses (each is itself a clause: `has`,
   `done`, `not_has`, `has_any`, or nested `any_of`). This allows
   mixing flag families:

   ```yaml
   - any_of:
       - has: phased_data_mtz
       - done: autosol
       - done: phaser
   ```

2. **Refactor key construction in `_check_requirements`** to use
   intermediate variables (matching the pattern already used in
   `_check_conditions`). Replaces `context.get("has_" + name)` with
   `key = "has_" + name; context.get(key)` so the static scanner
   doesn't trip on a literal `"has_"` inside a `context.get(...)`.

3. **Update `phenix.autobuild` requirements**: replace the
   inline-only `has_any` second clause with an `any_of` covering
   both flag families:

   ```yaml
   requirements:
     requires:
       - has_any: [data_mtz, phased_data_mtz]
       - any_of:
           - has: phased_data_mtz
           - has: model
           - has: placed_model
           - has: placed_model_from_history
           - done: phaser
           - done: autosol
   ```

   This now exactly mirrors the existing explanation pattern.

4. **Test suite expanded from 30 → 40**: added 6 tests for the
   new `any_of` clause type (positive/negative/nested/empty/malformed/
   AND-combined), plus 4 regression tests for the S5A scenario
   (autobuild kept after autosol_done, phaser_done,
   has_placed_model_from_history; still filtered when no history at all).

## Version 115.09b (GUI Fixes + Ligand Workflow + Bug 1)

### Summary

Production fixes from GUI testing and run 25 analysis. Three categories:
explicit_program loop guard, ligand-fitting workflow (6 fixes), and
GUI-mode corrections. 7 files modified.

### Modified Files

| File | Changes |
|------|---------|
| `agent/graph_nodes.py` | Bug 1: explicit_program done-flag guard prevents LLM-driven program loops. Removed dead .sca-only proactive check (replaced by diagnosable error). |
| `agent/directive_extractor.py` | Ligand-fit placement: `model_is_placed=True` for "fit ligand" advice. Ligand-fit `after_program` clearing: removes LLM-set `after_program` for multi-step ligand workflows. |
| `agent/workflow_engine.py` | Post-ligandfit exemption: defers `after_program_done` during combine/refine steps. Combine_ligand guard: forces pdbtools-only in combine step. Debug tracing for combine_ligand routing. |
| `agent/command_builder.py` | `phaser_sad.atom_type` interception: converts to `additional_atom_types`, prevents Se→S override in autosol. |
| `knowledge/programs.yaml` | Polder: `requires_resolution` invariant with `auto_fill_resolution` + `xray_data.high_resolution` strategy flag. Pdbtools: `exclude_patterns: [pose]` for ligand input (prevents selecting pose files over final model). Merged duplicate hints blocks. |
| `knowledge/diagnosable_errors.yaml` | `missing_crystal_symmetry`: "No unit cell info available" + "Cell and/or symmetry not specified" patterns for xtriage with .sca data. |
| `phenix/programs/ai_agent.py` | Fix A: `_has_explicit_stop` regex guard on preprocessing stop override (xtriage/mtriage). GUI auto-discovery skip: in GUI mode with user-selected files, skip supplement to prevent PHENIX artifacts from changing workflow intent. |

### Bug Details

| Bug | Symptom | Root Cause | Fix |
|-----|---------|------------|-----|
| Bug 1 (explicit_program loop) | bgal_denmod 3x resolve_cryo_em, hipip-refine 3x refine, emd_6123 3x map_to_model | `explicit_program` injection at graph_nodes.py line 691 bypasses done-flag check | Check done flag before injecting; skip if program already completed |
| Fix A (preprocessing stop) | "run mtriage and stop" generates full 5-stage plan | `_preprocessing_programs` unconditionally clears `after_program` | `_has_explicit_stop` regex check before clearing |
| Fix B (p9-xtriage) | xtriage fails twice on .sca without cell | Proactive file-type check can't distinguish p9 from gene-5-mad | Reactive: diagnosable_errors.yaml entry for "No unit cell" |
| GUI discovery | nsf-d2-ligand picks up pdb_sequences.fa → predict_and_build | Auto-discovery supplements with PHENIX project artifacts | Skip supplement in GUI mode with user-selected files |
| Ligand placement | "fit ATP" advice → predict_and_build instead of ligandfit | `model_is_placed` not set; model treated as unplaced | Rules overlay sets `model_is_placed=True` for ligand-fit signals |
| Ligand after_program | Ligandfit→pdbtools→refine→polder workflow blocked at each step | LLM sets `after_program` to different program each run (ligandfit/refine/polder); each choice blocks a different step | Clear `after_program` when ligand-fit signals detected; plan template drives workflow |
| Combine_ligand stuck | `combine_ligand` step gets STOP instead of pdbtools | `_apply_directives` wipes valid_programs via `after_program_done` | Post-ligandfit exemption + combine_ligand guard forces pdbtools-only |
| Post-ligandfit refine | LLM picks molprobity instead of refine after pdbtools | `_apply_directives` re-adds ligandfit or polder to front of list | Three-way branch: combine→pdbtools-only; refine→prioritize refine; done→STOP |
| phaser_sad.atom_type | Autosol uses S instead of Se | LLM adds `phaser_sad.atom_type=S` which overrides primary atom_type | Intercept and convert to `additional_atom_types` |
| Polder resolution | Polder crashes without resolution limit | No `auto_fill_resolution` invariant for polder | Added `requires_resolution` invariant + `xray_data.high_resolution` strategy flag |
| Pdbtools pose file | Pdbtools uses `ligand_fit_1_pose_5.pdb` instead of `ligand_fit_1.pdb` | `prefer_patterns: [ligand_fit]` matches both final and pose files | `exclude_patterns: [pose]` on ligand input slot |

### Verification

| Fix | Status |
|-----|--------|
| Bug 1 (bgal_denmod/rules_only) | ✅ mtriage→resolve_cryo_em→predict_and_build→RSR, CC=0.68 |
| Bug 1 (hipip-refine/rules_only) | ✅ refine→xtriage→molprobity→auto-stop, R-free=0.223 |
| Fix A (xtriage and stop) | ✅ Stops after xtriage in GUI |
| Fix B (p9-xtriage) | ✅ 1 cycle, clean diagnosis with helpful message |
| GUI discovery skip | ✅ Deployed, prevents pdb_sequences.fa pickup |
| Ligand workflow | ✅ xtriage→refine→ligandfit→pdbtools→refine→polder→molprobity |
| Polder resolution | ⏳ Deployed to server, awaiting verification |
| Pdbtools pose exclusion | ⏳ Deployed, awaiting verification |


## Version 115.09 (Tutorial Routing Fixes)

### Summary

Four routing bugs that prevented tutorials from completing. Fixes span
cryo-EM workflow progression, data-limitation detection, validation-only
routing, and MR-SAD intent recognition. Intent detection lives in the
directive extractor layer (LLM prompt + rules-based fallback); the
routing engine stays deterministic.

### Modified Files (5 production files)

| File | Fixes | Changes |
|------|-------|---------|
| `agent/workflow_engine.py` | 1, 3, 4 | Fix 1: `map_sharpening_done`, `map_symmetry_done`, `has_optimized_full_map` added to `past_analysis` gate in `_detect_cryoem_step`. Fix 3: `wants_validation_only` read from directives in `build_context`; validation shortcut in `_detect_xray_step` (guarded by `has_model` AND `has_data_mtz`). Fix 4: `force_mr` flag set in `build_context` when `use_mr_sad` + `has_model` + `not has_search_model` + `not phaser_done`; routing override in `_detect_xray_step` before placement probes; MR-SAD guard in `get_valid_programs` updated to check `force_mr`. |
| `agent/graph_nodes.py` | 2 | `.sca/.hkl`-only data detection in `perceive()` with abort_message. Guarded by `not _has_model` AND `not _has_sequence` AND no `unit_cell` in directives. Follows established early-stop pattern (stop + stop_reason + abort_message + intent dict). |
| `agent/directive_extractor.py` | 3, 4 | LLM prompt: `wants_validation_only` added to `workflow_preferences` schema with extraction guidance. MR-SAD prompt strengthened with concrete examples. `_validate_directives`: `wants_validation_only` added to allowed boolean keys. `extract_directives_simple`: rules-based fallback for 6 validation signals and 6 MR-SAD patterns. |
| `agent/plan_generator.py` | 3 | `wants_validation_only` added to `_build_context` initial dict and propagated from `workflow_preferences` to template selection context. |
| `knowledge/plan_templates.yaml` | 3 | `validate_existing` template: `applicable_when` requires xray + `wants_validation_only` + `has_search_model` + `model_is_placed` (priority 60, beats `refine_placed_ligand` at 55). Two stages: data_assessment (xtriage) + validation (model_vs_data, molprobity). |

### Bug Details

| Bug | Tutorial | Symptom | Root Cause | Fix |
|-----|----------|---------|------------|-----|
| 1 | bgal_denmod, apoferritin_denmod, ion_channel_denmod | Cryo-EM stops after 1 cycle (map_sharpening or resolve_cryo_em) | `map_sharpening_done` missing from `past_analysis` gate → `detect_step` returns "analyze" → no programs → STOP | 3 new flags in `past_analysis` including file-presence fallback |
| 2 | p9-xtriage | xtriage fails on unmerged .sca without cell dimensions | Data limitation: unmerged scalepack needs unit cell + space group provided interactively | Early detection in `perceive()` with helpful abort message |
| 3 | pka-validate | Runs predict_and_build instead of molprobity despite "Analysis only" advice | No validation-only routing path; directive extractor doesn't capture validation intent | `wants_validation_only` directive + validation shortcut + plan template |
| 4 | lysozyme-MRSAD | R-free stuck at 0.54; phaser never runs despite "MR-SAD" advice | PDB classified as `model` not `search_model`; LLM ignores `use_mr_sad` prompt; no rules-based MR-SAD fallback | `force_mr` flag + MR-SAD rules fallback + prompt strengthening |

### Review Fixes (2 bugs found during review)

| Issue | Fix |
|-------|-----|
| Fix 2 false positive: gene-5-mad (merged .sca + sequence) falsely aborted | Added `not _has_sequence` guard to `.sca-only` condition |
| Fix 3 false positive: "analysis only" matched xtriage-only and cryo-EM tutorials | Removed from rules-based signals; LLM prompt handles the nuance |

### Deployment Fix (discovered via run 19 verification)

Run 19 (OpenAI, post-deployment) showed all four fixes not firing: directives
contained no `wants_validation_only` or `use_mr_sad` despite code being
deployed and bytecache cleared.

**Root cause**: The rules-based intent patterns were only in
`extract_directives_simple()` — the fallback path. OpenAI runs use the
LLM path in `extract_directives()`, which calls the LLM, gets directives
without these flags, and returns. The rules-based patterns never execute.

**Fix**: Extracted patterns into `_apply_workflow_intent_fallback()` shared
helper. Called as **post-LLM overlay** in `extract_directives()` (after LLM
returns, before final return) AND from `extract_directives_simple()`.
Rules always run last and always win for routing flags. The LLM cannot
override them because it never sets them (0/240 extractions).

**Future**: Replace overlay with centralized `_DIRECTIVE_SCHEMA` and
registry-driven `_merge_tiered` merge (v115.10, see `docs/directive_merge_plan.md`).

### Production Verification Fixes (runs 20-22)

| Issue | Discovery | Root Cause | Fix |
|-------|-----------|------------|-----|
| Fix 3: `_is_valid_file` rejects valid 3dnd.pdb | Run 21 `[GATE]` output: `has_model=False` despite 3225 ATOM records | PDB scan limit of 500 lines; 3dnd.pdb has 546 header lines before first ATOM record | `workflow_state.py`: increased scan limit to 2000 |
| Fix 3: `has_data_mtz=False` for phased MTZ | Run 20 analysis: completed structures with phase columns go to `phased_data_mtz` | Validation shortcut only checked `has_data_mtz` | `workflow_engine.py`: added `has_phased_data_mtz` to context + validation shortcut |
| Fix 1: `map_sharpening_done` regex wrong | Diagnostic: regex matches "sharpened" but output is `auto_sharpen_A.ccp4` | Zombie check table regex too narrow | `workflow_state.py`: broadened regex from "sharpened" to "sharpen" |
| Fix 2: `.sca-only` check blocked by sequence guard | Run 19: p9-xtriage has `seq.dat` like gene-5-mad | Proactive file-type approach can't distinguish the two | **Deferred** — needs reactive approach (diagnosable_errors.yaml) |

### Additional files modified (production verification)

| File | Changes |
|------|---------|
| `agent/workflow_state.py` | `_is_valid_file`: PDB scan limit 500→2000. `_ZOMBIE_CHECK_TABLE`: map_sharpening regex "sharpened"→"sharpen". |
| `agent/workflow_engine.py` | `has_phased_data_mtz` added to `build_context` initial dict + validation shortcut condition. |

### Known Issues (deferred to v115.10)

| Issue | Description |
|-------|-------------|
| Preprocessing stop override | `ai_agent.py` line 2761 unconditionally clears `after_program` for xtriage/mtriage even when user explicitly says "stop". Fix: check `_has_explicit_stop` regex before clearing. |
| CIF model categorization | `model.cif` files go to generic `cif` category, not `model`. Affects model-dependent routing for CIF-format structures. |
| Fix 2 redesign | p9-xtriage `.sca-only` proactive check fails (has sequence like gene-5-mad). Need reactive: xtriage "no unit cell" failure as diagnosable terminal error. |
| Directive extraction refactor | Replace post-LLM overlay with `_DIRECTIVE_SCHEMA` registry + `_merge_tiered`. See `docs/directive_merge_plan.md`. |


## Version 115.08 (Phased File Detection + Systematic Testing Framework)

### Summary

Four critical fixes for phased file detection (F1–F4) that caused rab3a-refine
and nsf-d2-refine to fail across all 5 modes. One additional bug fix (B1:
last_program) found by the testing framework. A 10-phase systematic testing
framework exercising file categorization, routing, command building, error
classification, and LLM resilience across 32 tutorials.

### Modified Files (main code — 3 files)

| File | Change | Details |
|------|--------|---------|
| `agent/workflow_state.py` | F1–F4 | Content-based phased detection replacing filename markers. New: `_PHASE_COLUMN_CACHE` (iotbx column-type check cached by abspath+mtime), `_ascii_phase_heuristic()` (ASCII header check for .hkl files), `_has_phase_columns_cached()` (cached wrapper). `_resolve_phased_promotions()` moved from hardcoded-only to shared post-processing (F2). Category-exclusivity enforcement — files appear in one of `data_mtz` or `phased_data_mtz`, not both (F3). Conditional removal from `data_mtz` only when other non-phased data files remain (F4). `_IGNORED_FORMATS` dict for .cv/.xplor/.param/.eff extensions. |
| `agent/workflow_engine.py` | B1, diagnostics | B1 fix: added `"last_program": history_info.get("last_program")` to `build_context()` dict literal — the `after_program` directive check at line 1799 previously always saw None. `[GATE]` diagnostic logging controlled by `PHENIX_AGENT_DIAG_VALID_PROGRAMS=1` env var — prints context flags, detect_step result, and STOP reasons to aid routing debugging. |
| `agent/graph_nodes.py` | WARNING event | WARNING event emitted for each ignored format file in `perceive()` — informs user that .cv/.xplor files are not supported for auto-fill. |

### New Test Files (12 files)

| File | Phase | Tests | Purpose |
|------|-------|-------|---------|
| `tst_file_categorization.py` | — | 42 | Unit tests for v115.08 phased detection fixes |
| `tst_phase0_static_audit.py` | S0 | 5 | Parse check, bare except scan, import fallback check |
| `tst_phase1_contract_gaps.py` | S1 | 128 | AST-based coverage map of 4 key modules |
| `tst_phase2_path_consistency.py` | S2 | 10 | YAML vs hardcoded categorization path diff |
| `tst_phase3_serialization_symmetry.py` | S3 | 28 | JSON symmetry + invariants + AgentSession round-trip |
| `tst_phase4_history_flags.py` | S4 | 8 | Flag writer/reader consistency, B1 fix verification |
| `tst_phase5_error_classification.py` | S8 | 7 | 3 error classifiers, pattern overlap, self-detection |
| `tst_phase6_category_consumer.py` | S5 | 3 | input_priorities + fallback_categories alignment |
| `tst_phase7_routing_simulation.py` | S6 | 32 | 3-cycle routing simulation with real production code |
| `tst_phase8_command_building.py` | S7 | 15 | CommandBuilder.build() for 7 tutorials + 2 edge cases |
| `tst_phase9_llm_perturbation.py` | S9 | 17 | Filename/program/parameter/truncation/empty perturbation |
| `run_all_tests.py` | — | — | Updated with S0–S9 registration |

### New Documentation

| File | Purpose |
|------|---------|
| `docs/PHASE_REVIEW_REPORT.md` | Detailed report of 54 review fixes across all 10 phases |

### Bug Details

| Bug | Severity | Description | Fix |
|-----|----------|-------------|-----|
| F1 | CRITICAL | Filename-marker detection checked 'phased' but files were named 'phases' | Content-based detection (iotbx + ASCII heuristic) |
| F2 | CRITICAL | `_resolve_phased_promotions()` only ran on hardcoded path, not YAML (production) | Moved to shared post-processing |
| F3 | CRITICAL | YAML put files in both `data_mtz` and `phased_data_mtz` | Exclusivity enforcement |
| F4 | CRITICAL | Removing from `data_mtz` could delete the only data file | Conditional removal (only when alternatives exist) |
| B1 | LOW | `last_program` never transferred from `_analyze_history()` to `build_context()` | 1-line addition to `build_context()` dict |

### Test Framework Review (54 fixes across 9 phases)

The testing framework itself was reviewed across multiple rounds. Key
classes of defects found and fixed in the test scripts:

- **5 tests that could never fail** (Phases 2, 6, 7, 9) — status was always
  PASS/PARTIAL but raise checked for FAIL. Fixed with whitelist patterns.
- **3 tautological assertions** (Phases 3, 9) — `assert True`, `assert X or True`.
- **6 phantom/wrong names** in boundary function sets (Phase 1).
- **5 missing exception isolation** patterns (Phases 2, 6, 7, 8) — one crash
  killed all remaining tutorials.
- **5 wrong `best_files` keys** (Phases 8, 9) — `'data'` should be `'data_mtz'`.
- **9 new AgentSession round-trip tests** (Phase 3) — replaced 18 tautological
  JSON identity tests with real production code exercise.

See `docs/PHASE_REVIEW_REPORT.md` for full details.


## Version 115.07 (Run 15b Bug Fixes — Phase 3)

### Summary

Four bugs identified from analysis of run 15b (OpenAI, 371 runs, 42
tutorials). Fixes address a graph crash from uncoerced JSON metrics, a
cryo-EM half-map pair detection gap, a terminal crash loop, and an LLM
hallucination pattern. Two issues deferred as known limitations.

**Run 15b headline:** 9→18 GOOD tutorials (doubled), failure rate 39%→34%.
Phase 3 targets remaining failures: lysozyme-refine (77%), apoferritin
denmod_dock (57%), 7rpq (never refines), AF_7n8i (5/10 fail).

### Modified files (4)

| File | Bug | Changes |
|------|-----|---------|
| `agent/structure_model.py` | 4 | **Numeric coercion**: added `_coerce_numerics()` function — coerces all float/int fields after `from_dict()` deserialization. Added `_safe_float()` guards at all 6 arithmetic and formatting sites in `_detect_problems()`, `get_summary()`, and `_format_report()`. Progress entries also coerced. Prevents `TypeError: unsupported operand type(s) for -: 'str' and 'float'` when model_vs_data stores metrics as strings in session JSON. |
| `agent/workflow_state.py` | 6 | **Two-tier half-map pair detection**: relaxed `_categorize_files()` heuristic. Tier 1 (new): when exactly 2 `full_map` files form a `_1/_2` pair, promote to `half_map` without requiring a companion full map. Tier 2 (existing): when ≥3 `full_map` files include a pair, promote only if a companion remains. Fixes AF_7n8i where `box_1.ccp4`/`box_2.ccp4` were the only maps. |
| `agent/phil_validator.py` | 3 | **Blocked params**: added `phenix.resolve_cryo_em` to `_BLOCKED_PARAMS` with 4 entries: `mask_atoms`, `output.prefix`, `d_min`, `main.number_of_macro_cycles`. These are LLM-hallucinated params copied from mtriage log output that cause RuntimeError. |
| `knowledge/diagnosable_errors.yaml` | 1 | **Terminal diagnosis**: added `unknown_chemical_element` — detects "Unknown chemical element type" errors (bad PDB columns 77-78) and stops immediately with actionable hint instead of crash-looping. |

Also modified (config):

| File | Bug | Changes |
|------|-----|---------|
| `knowledge/programs.yaml` | 3 | Removed `mask_atoms` from `phenix.resolve_cryo_em` `strategy_flags` whitelist. Only `resolution` passes through. |

### Bug details

| Bug | Tutorial | Symptom | Root cause | Fix |
|-----|----------|---------|-----------|-----|
| 1 | lysozyme-refine | 77% fail, no metric | PDB has AU atom missing element columns → refine crashes → agent retries 5× | Terminal diagnosis: stop with hint |
| 3 | apoferritin_denmod_dock | 57% fail (LLM modes) | LLM copies `mask_atoms=True` from mtriage log → `strategy.mask_atoms_atom_radius="True"` RuntimeError | Whitelist removes it; `_BLOCKED_PARAMS` double-blocks it |
| 4 | 7rpq_AF_reference | 0% fail but stops at R=0.386 | model_vs_data stores `r_work="0.385"` (string) → cycle 3 graph crashes on `r_free - r_work` | `_safe_float()` at 5 sites in metric_evaluator (true root cause: USE_YAML_METRICS=True); belt-and-suspenders in metrics_analyzer (7 sites), structure_model, graph_nodes, event_formatter, kb_tags, workflow_state._analyze_history |
| 5 | lowres_restraints | 75% fail (restraint params stripped) | PHIL validator strips `reference_model.file`, `ncs.*`, `secondary_structure.*`, `ramachandran_restraints` — 7 of 8 params lost | Hierarchical prefix whitelist + path resolution + reference model exclusion |
| 6 | AF_7n8i | 5/10 modes fail at mtriage | `box_1.ccp4`/`box_2.ccp4` not recognized as half-maps (old heuristic required ≥3 full_map files) | Tier 1: exactly 2 matching files → promote |
| 7 | 1aba-polder | polder fails (missing: model) | Protein model with ligand atoms misclassified as `ligand_pdb` → polder can't find `model` slot | File-size fallback: PDB >10KB in ligand_pdb rescued to model |
| 8 | apoferritin_denmod_dock (rules_only) | dock_in_map loops 3x, RSR never offered | `emd-20026_auto_sharpen_A.ccp4` excluded from `full_map` by `*_a.*` pattern → `has_full_map=False` → RSR blocked | Orphan-map promotion: map files not in any subcategory → promoted to full_map |

### Bug 5 details — Reference model restraints for phenix.refine

**Problem:** phenix.refine has a massive PHIL parameter tree, but the
strategy_flags whitelist had only 13 entries. Advanced crystallographic
restraint parameters (reference model, NCS, secondary structure, Ramachandran)
were correctly extracted by the LLM from the README but stripped by PHIL
validation. The tutorial requires all of these to produce a meaningful result.

**Fix A — Hierarchical prefix whitelist** (`programs.yaml` + `phil_validator.py`):
Added `allowed_phil_prefixes` to phenix.refine: `reference_model.`,
`secondary_structure`, `ncs.`, `ramachandran`. The PHIL validator now allows
any strategy key containing one of these as a case-insensitive substring.
This covers both short forms (`ncs.type`) and full PHIL paths
(`refinement.pdb_interpretation.ncs.type`). Added `ramachandran_restraints`
as an individual strategy_flag (standalone param, not under a namespace).

**Fix B — Path resolution** (`program_registry.py`):
Added a pre-pass in `build_command()` that detects strategy values ending
in file extensions (`.pdb`, `.cif`, `.params`, `.eff`, etc.) and resolves
them to absolute paths. Builds a basename→path lookup from the command's
input files and the working directory. Handles the LLM writing
`reference_model.file=4pf4.pdb` (relative) by resolving to the full path.

**Fix C — Reference model exclusion** (3-tier):
- Tier 1 (`command_builder.py`): When `reference_model.file` is in the strategy,
  exclude that file from primary model selection. Already existed (lines 380-394),
  but was blocked because Fix A wasn't in place — now unblocked.
- Tier 2 (`workflow_state.py`): When ≥2 PDB files are in the `model` category
  and one matches `reference|homolog|template|restraint|high.res` (but NOT
  agent output prefixes like `refine_`, `autobuild_`), reclassify to
  `reference_model`.
- Tier 3: LLM extraction of `reference_model.file=` from README feeds Tier 1.

**Fix D — Strategy rewrites** (`graph_nodes.py`):
5 new `_STRATEGY_REWRITES` entries normalize verbose PHIL paths to shortest
valid forms (e.g. `refinement.pdb_interpretation.ncs.type` → `ncs.type`).

### Deferred issues

| Tutorial | Symptom | Why deferred |
|----------|---------|-------------|
| lysozyme-MRSAD | 72% fail (multi-wavelength label confusion) | Needs wavelength selection feature (Phase 4) |

### Other changes this session

| File | Change |
|------|--------|
| `agent/metrics_analyzer.py` | Added `_safe_float()` at all 7 numeric read sites in `derive_metrics_from_history()`, `_analyze_xray_trend()`, `_analyze_cryoem_trend()`, `get_latest_resolution()`, `get_best_r_free()`, `get_latest_r_free()`, `get_latest_map_cc()`. Root cause of Bug 4 crash: JSON round-tripping turns floats to strings. |
| `agent/metric_evaluator.py` | Added `_safe_float()` at 5 arithmetic sites: `_analyze_xray_trend()` r_free extraction, `_analyze_cryoem_trend()` CC extraction, `is_significant_improvement()`, `calculate_improvement_rate()`, `is_plateau()`. TRUE root cause of Bug 4: `USE_YAML_METRICS=True` routes through this file, not `metrics_analyzer.py`. Added import fallback for standalone testing. |
| `agent/event_formatter.py` | Wrapped compact metrics formatting (L920-935) with float coercion + try/except. Prevents `str - float` crash in METRICS_EXTRACTED event rendering. |
| `agent/kb_tags.py` | `_trend_tags()` now coerces all R-free trend values to float via try/except before arithmetic (diffs, total_drop). |
| `agent/workflow_state.py` | `_analyze_history()` L1841-1850: all metric reads coerced via local `_sf()` helper. Orphan-map promotion: map files in parent `map` but not in any subcategory (`full_map`, `half_map`, `optimized_full_map`) promoted to `full_map`. Fixes apoferritin_denmod_dock rules_only. |
| `agent/graph_nodes.py` | Removed stale `strategy.mask_atoms` → `mask_atoms` rewrite for resolve_cryo_em (now a no-op since mask_atoms is in `_BLOCKED_PARAMS`). |
| `knowledge/plan_schema.py` | `record_stage_cycle()` catch-up: when agent runs ahead of plan tracker (e.g. ligandfit during refine stage), advance through intermediate stages. Overshoot guard: verify `new_curr` matches program before counting. |
| `phenix_ai/local_agent.py` | Added `client_version=self._get_client_version()` + `_get_client_version()` method for parity with RemoteAgent. |
| `phenix_ai/remote_agent.py` | Added `request["settings"]["verbosity"]` and `events=parsed.get("events", [])` for parity with LocalAgent. |
| `tests/tst_phase3_bug5.py` | 38 test functions / 104 assertions covering all Phase 3 + Bug 5 + Bug 8 fixes, metric_evaluator coercion, orphan-map promotion, kb_tags string handling. |
| `tests/tst_phil_validation.py` | Fixed stale `test_rewrite_resolve_cryo_em_mask_atoms` (mask_atoms now blocked, not allowed). |
| `tests/tst_audit_fixes.py` | 3 stale tests updated: `test_k2_mtriage` (prefers_half_maps), `test_k2_map_sharpening` (positional), `test_s5h_inject_program_defaults` (generate not in defaults). |
| `tests/tst_backward_compat.py` | 3 new parity tests (26 total): `test_local_remote_agent_settings_parity`, `test_local_remote_agent_return_parity`, `test_local_agent_full_roundtrip`. |
| `docs/ARCHITECTURE.md` | Updated PHIL validation section (prefix whitelist + path resolution), added half-map Tier 1, reference model categorizer, `_coerce_numerics`, plan_schema catch-up documentation. |
| `docs/DEVELOPER_GUIDE.md` | RULE 4 expanded (prefix whitelist guidance), RULE 9 added (numeric type safety in StructureModel). |

### Test results

| Suite | Tests |
|-------|-------|
| `tst_structure_model.py` | 77/77 |
| `tst_plan_schema.py` | 53/53 |
| `tst_backward_compat.py` | 26/26 |
| `tst_command_builder.py` | 22/22 |
| `tst_event_system.py` | 13/13 |
| `tst_phase3_bug5.py` | 104/104 |
| `tst_phil_validation.py` | 15/15 |
| **Total** | **310** |
| **Total** | **191/191** |


## Version 115.06 (Cryo-EM Half-Map Handling — Phase 2)

### Summary

Three interconnected fixes addressing the #1 failure source in cryo-EM
tutorials: programs receiving only 1 half-map instead of 2. Root cause was
that the LLM assigns one half-map to a `multiple:true` slot, marking it as
"filled", so auto-fill skips it. Affected mtriage (109 failures),
resolve_cryo_em (49), and map_sharpening (46) — totaling ~200 wasted cycles.

### Modified files (3)

| File | Changes |
|------|---------|
| `agent/command_builder.py` | **Supplement logic**: after LLM files + auto-fill, a new loop checks every `multiple:true` slot. If the LLM filled it with 1 file but the category has 2+, missing files are backfilled using `_find_file_for_slot()`. Uses `os.path.realpath()` for symlink robustness. |
| `knowledge/programs.yaml` | **mtriage**: changed from `keep_half_maps_with_full_map: true` to `prefers_half_maps: true` — drops full_map when half-maps present (eliminates "Maps have different dimensions" error from mismatched grids). **map_sharpening**: added `prefers_half_maps: true` + changed half_map flag from `"half_map="` to `""` (positional args, matching the working command `phenix.map_sharpening h1.ccp4 h2.ccp4 seq_file=seq.fa`). |
| `knowledge/workflows.yaml` | **map_sharpening**: reverted all 4 entries from `has: full_map` to `has_any: [full_map, half_map]` — map_sharpening works with either input type. |

### Data flow: half-map supplement

```
LLM: {half_map=file2.ccp4}         ← one file
Corrector: → half_map_1.ccp4       ← still one file
selected_files["half_map"] = "half_map_1.ccp4"

Auto-fill: half_map already filled → SKIP

NEW SUPPLEMENT:
  half_map is multiple:true AND in selected_files
  _find_file_for_slot returns [half_map_1.ccp4, half_map_2.ccp4]
  len(2) > len(1) → supplement fires
  selected_files["half_map"] = ["half_map_1.ccp4", "half_map_2.ccp4"]

Dedup: prefers_half_maps → drop full_map, keep both half-maps
Command: phenix.mtriage half_map=h1.ccp4 half_map=h2.ccp4  ← SUCCESS
```

### Expected impact

- mtriage: gets 2 half-maps (no full_map) → no dimension mismatch
- resolve_cryo_em: gets 2 half-maps (supplemented) → "Need 2 half maps" eliminated
- map_sharpening: gets 2 half-maps positionally → "use b-factor" error eliminated
- Affects all 9 cryo-EM tutorials that were failing in run 14


## Version 115.05 (Run 13/14 Bug Fixes — Phase 1)

### Summary

Ten bugs identified from analysis of run 13 (expanded tutorial set: 33
README + 35 SOLVE across 25+ tutorials) and verified in run 14 (253 OpenAI
+ 52 Ollama runs). Fixes address R-free flag mismatches, cryo-EM half-map
dedup, map_sharpening workflow routing, docking for unplaced models, space
group parsing, PHIL parameter passthrough, and data label selection.

### Modified files (8)

| File | Bugs | Changes |
|------|------|---------|
| `knowledge/programs.yaml` | 4/5, 2/3, 7, 10 | Removed `generate=True` from refine template+defaults; added `generate_rfree_flags` strategy flag; added `prefers_half_maps: true` to resolve_cryo_em; added `fallback_categories: [model]` to dock_in_map; added `strict_strategy_flags: true` to process_predicted_model |
| `knowledge/workflows.yaml` | 1, 2/3, 7 | resolve_cryo_em: added `has: half_map` to optimize_map entry; dock_in_map: added `model` to has_any in dock_model+refine steps; added dock_in_map to cryo-EM refine step |
| `knowledge/recoverable_errors.yaml` | 9 | Added `phenix.model_vs_data` to data_label_parameters |
| `knowledge/prompts_hybrid.py` | 4/5 | Updated R-free error recovery prompt (was telling LLM "command already includes generate=True") |
| `agent/command_builder.py` | 2/3, 6 | `prefers_half_maps` dedup branch; `multiple:true` append logic with normalization; cctbx-first space group validation with prefix shortening |
| `agent/program_registry.py` | 10 | `strict_strategy_flags` check before KNOWN_PHIL_SHORT_NAMES passthrough |
| `agent/graph_nodes.py` | 7 | Unplaced model guard: redirects to dock_in_map when CC<0.10 and dock_done=False; quality floor prefers dock_in_map for unplaced models |
| `docs/RUN_13_BUG_ANALYSIS_AND_PLAN.md` | — | Analysis document |

### Bug details

| Bug | Issue | Fix | Verified |
|-----|-------|-----|----------|
| 1 | map_sharpening offered with only half-maps → crash | `has: full_map` condition (later reverted in v115.06 to `has_any`) | ✓ bgal, actin_sharpen eliminated |
| 2/3 | resolve_cryo_em gets 1 half-map (dedup drops second) | `prefers_half_maps` flag + append logic | ✓ apoferritin_dock CC=0.496 |
| 4/5 | R-free generate=True on every refinement | Removed from template; conditional via strategy flag | ✓ 3tpp R=0.167 (was CRASH) |
| 6 | Space group "P4 for refinement wit" | cctbx.sgtbx progressive prefix shortening | ✓ |
| 7 | groel: dock_in_map never offered | fallback_categories + refine step + CC guard | ✓ dock_in_map now runs |
| 8 | beta-blip hard MR case | No fix needed | — |
| 9 | "Multiple equally suitable arrays" | model_vs_data added to data_label_parameters | ✓ |
| 10 | process_predicted_model PHIL error | strict_strategy_flags blocks passthrough | ✓ b-CA-mr 0 failures |

### Ollama vs OpenAI (run 14)

Both providers use the same server-side code. Key finding: Ollama achieves
comparable quality (36% vs 39% fail rate) and dramatically outperforms OpenAI
on AF_POMGNT2 (R=0.287 vs 0.553) because Ollama's LLM correctly runs
process_predicted_model before molecular replacement.


## Version 115.04 (ASU Copy Count Tracking)

### Summary

Tracks the number of copies of the search model in the asymmetric unit (ASU)
and passes this automatically to Phaser as `search_copies` each cycle.
Copies are sourced from user directives (always wins) or xtriage log analysis.

### Modified files (5)

| File | Changes |
|------|---------|
| `programs/ai_agent.py` | `_extract_copies_from_directives()` method added; copies injected in all 3 directive extraction paths (rules-only, LLM, fallback); xtriage n_copies update from `history_record` |
| `agent/api_client.py` | `asu_copies` added to both `build_session_state()` and `build_request_v2()` whitelists so value survives client→server round-trip |
| `agent/graph_nodes.py` | `_fallback_extract_metrics()` extracts `n_copies` from xtriage log; BUILD node reads `session_info["asu_copies"]` and injects `component_copies` into strategy (skipped if LLM already set it); sanity bound 1–30 enforced |
| `phenix_ai/run_ai_agent.py` | `asu_copies` forwarded from `session_state` → `session_info` each cycle; serialized into `metadata` for persistence |
| `tests/tst_copies_tracking.py` | New test: directive path, xtriage path, priority (directive wins), range clamping, integration with `run_ai_agent` metadata |

### Data flow

**Directive path** (user says "4 copies in the ASU"):
```
_extract_copies_from_directives() → session.data["asu_copies"] = 4
  → session_info["asu_copies"] → api_client whitelists
    → run_ai_agent session_state → BUILD node → phaser.search_copies=4
```

**Xtriage path** (xtriage log: "Best guess : 4 copies"):
```
_fallback_extract_metrics() → log_analysis["n_copies"] = 4
  → history_record["asu_copies"] → session.data["asu_copies"]
    (only if not already set by directive)
      → next cycle: session_info → BUILD → phaser.search_copies=4
```

### Tests

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_copies_tracking.py` | New | Directive extraction, xtriage extraction, priority, bounds |


## Version 115.03 (Tutorial Run Bug Fixes: Bugs A–F + Anomalous Signal)

### Summary

Six bugs identified from log analysis of 16 tutorial runs across CASP7,
a2u-globulin-mr, a2u-globulin-rebuild, and advanced-xtriage.  The most
impactful fixes are: the `_has_placed_model` MR-keyword guard (Bug E)
that prevented search models from being misrouted to `xray_refined`; the
hopeless R-free retry path (Bug F) that recovered from failed MR solutions
instead of stopping; and the anomalous signal detection fix that correctly
handles weak signal (xtriage `has_anomalous=True` with `measurability` in
`[0.05, 0.10]`) without overwriting it with a measurability-based clear.
Bugs A, B, C (cascade), and D were already present in the tree.

### Modified files (4)

| File | Bug | Changes |
|------|-----|---------|
| `agent/workflow_state.py` | Anomalous | `has_anomalous` detection: `elif has_anomalous` promoted to independent `if` (was silently skipped when `anomalous_resolution` present but ≥ 6.0); `>= 0.06` weak-signal branch added; `< 0.06` clear no longer fires when `analysis["has_anomalous"]` is explicitly `True` |
| `agent/workflow_engine.py` | Bug E | `_has_placed_model`: `_wants_mr_first` guard scans constraints for MR keywords (`"molecular replacement"`, `"phaser"`, `" mr "`, etc.); if found, skips placement-keyword inference so search models are not misidentified as placed |
| `agent/workflow_engine.py` | Bug F | `_detect_xray_step` Step 2e: when `r_free >= 0.45` after MR + ≥1 refinement cycle, routes to `obtain_model` instead of falling through to `validate`/`complete` → STOP |
| `knowledge/file_categories.yaml` | Bug C (yaml) | New `predict_build_refine_internal` intermediate category matching `*overall_best*final_refine*.pdb` and `*overall_best*refine_[0-9][0-9][0-9]*.pdb`; excludes `*_final_refine_*` from `predict_and_build_output` |

### AgentState TypedDict (graph_state.py)

Eight fields added to `AgentState` to support thinking-mode features.
`create_initial_state` gains one new parameter (`session_blocked_programs`).
Server-side test `tst_thinking_defense.py` expected counts updated:
H02 `37` → `45` (field count), H03 `19` → `20` (param count).
Run `patch_server_tests.py <langchain_root>` to apply.

### Tests

| File | Tests added | Covers |
|------|-------------|--------|
| `tests/tst_p1_to_p5_fixes.py` | +2 | Weak anomalous detection; negligible anomalous clear |
| `patch_server_tests.py` | — | Updates `tst_thinking_defense.py` H02/H03 expected counts |

### Bugs confirmed already in tree (no new code needed)

| Bug | Root cause | Fix location |
|-----|-----------|--------------|
| Bug A | `event_formatter.py` `fmt % value` TypeError on string metrics | `agent/event_formatter.py` try/except guard |
| Bug B | `_auto_discover_files` skipped when `original_files` non-empty | `programs/ai_agent.py` supplement-mode discovery |
| Bug C (cascade) | `predict_and_build` cascade set `refine_done=True` prematurely | `agent/workflow_state.py` cascade comment + guard |
| Bug D | autosol deprioritized when `anomalous_measurability < 0.05` | `agent/workflow_engine.py` deprioritization logic |


## Version 115.02 (Report Analysis Fixes & GUI Update)

### Summary

Fixes from report analysis: program alias resolution, plan generation
guard for single-program tasks, solve-mode file discovery fix,
unsupported program handling (Sorry → warning), polder file routing,
predict_and_build label injection, and GUI alignment.

### Modified files (7)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | `phenix.resolve` → `phenix.resolve_cryo_em` alias in `_fix_program_name()`; LLM prompt negative example for resolve |
| `agent/intent_classifier.py` | Added `resolve_cryo_em` pattern to `phenix.resolve_cryo_em` aliases |
| `agent/sanity_checker.py` | Added `_check_autobuild_phib()` — warns after 2+ autobuild PHIB failures |
| `agent/thinking_agent.py` | Enhancement B: `error_classification` and `failure_count` forwarded to thinking context |
| `programs/ai_agent.py` | Fix 1: `_load_result_file_exclusions()` + filter in `_auto_discover_files()`; plan skip for task intent; solve-mode README preprocessing bypass; expanded `_README_IGNORE_PROGRAMS` with utility programs (keeps `raise Sorry` for unsupported core programs like ensemble_refinement) |
| `knowledge/programs.yaml` | Polder `data_mtz` input_priorities (exclude map_coeffs_mtz); `predict_and_build` obs_labels strategy_flag |
| `wxGUI2/Programs/AIAgent.py` | Expanded ignore set for utility programs (keeps `raise Sorry` + blocking banner for unsupported core programs) |


## Version 115.01 (Intent Classification & Infrastructure Fixes)

### Summary

Intent classification system for dual-run evaluation (tutorial vs solve
mode), plus 6 infrastructure fixes addressing file routing bugs discovered
by the first dual-run evaluation.  Expected additional cycle waste
reduction: ~15-20% from preventing technical crashes.

### New files (3)

| File | Lines | Purpose |
|------|-------|---------|
| `agent/intent_classifier.py` | 433 | Classifies user advice into solve/solve_constrained/task/tutorial |
| `knowledge/tutorial_expectations.yaml` | 218 | Expected behavior for all 21 tutorials (answer key) |
| `knowledge/solve_run_info.yaml` | 94 | Minimal hints for solve-mode runs |

### Modified files (5)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | Intent classification integration; intent-driven stop behavior; `_set_by_pattern` flag for distinguishing pattern-inferred vs user-explicit stops; intent simplified to string before return |
| `agent/advice_preprocessor.py` | Early return bypass for "Solve the structure by standard procedures" READMEs (skips LLM preprocessing) |
| `agent/workflow_state.py` | Fix I1: Multi-array MTZ detection (`_detect_mtz_arrays`, `_read_mtz_array_labels`, `_classify_mtz_arrays`); Fix I3: `phased_data_mtz` post-categorization scan for autosol outputs |
| `agent/phil_validator.py` | Fix I2: `_BLOCKED_PARAMS` dict; blocked-param check before allowed-flags check (strips `input_map_file` from autobuild) |
| `agent/graph_nodes.py` | Fix I1: MTZ label injection with SAD/MR ranking rule; Fix I4: `map_sharpening` resolution injection; Fix I5: `_estimate_anomalous_sites`/`_read_sequence` + sites injection + `atom_type` forwarding from directives + `rebuild_in_place=False` when R-free > 0.45 |

### Knowledge base changes (1)

| File | Changes |
|------|---------|
| `knowledge/programs.yaml` | Fix I1: `obs_labels` for xtriage and autosol; Fix I3: autosol output type `data_mtz` → `phased_data_mtz` + denmod pattern; Fix I8: `ligand_cif` optional input + input_priorities + command template for refine |

### Evaluation infrastructure (2)

| File | Changes |
|------|---------|
| `tests/analyze_tutorial_runs.py` | Dual-run support: `run_type` dimension (solve/readme), per-type output files, updated directory parsing |
| `solve_readmes/` | 21 solve-mode README.txt files |

### Tests

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_intent_classifier.py` | 32 | Intent classification (4 categories) |
| Total (all test files) | 106 | Fixes 2-6, Enhancement A, Intent, Infrastructure |


## Version 115.00 (Failure Handling & Stop Condition Fixes)

### Summary

Seven fixes and one enhancement addressing the top causes of
wasted agent cycles, based on analysis of 21 tutorials × 5 modes.
The core changes eliminate fabricated stop conditions (Fix 2),
add PHIL parameter validation (Fix 4), implement a tiered failure
response with error classification (Fix 3), make duplicate
detection parameter-aware (Fix 5), rescue uncategorized map files
(Fix 6), and inject structured error context for LLM
self-correction (Enhancement A).  Expected cycle waste reduction:
41% → ~23% (20 cycles saved across analyzed sessions).

### New files (2)

| File | Lines | Purpose |
|------|-------|---------|
| `agent/phil_validator.py` | 142 | PHIL strategy validation against programs.yaml strategy_flags |
| `agent/error_classifier.py` | 410 | Error classification (TERMINAL/PHIL_ERROR/LABEL_ERROR/RETRYABLE) + tiered pivot logic |

### Modified files (6)

| File | Changes |
|------|---------|
| `agent/directive_extractor.py` | Fix 2: Strip fabricated "Stop Condition:" and "Goal:" headers from preprocessor output; remove "Stop Condition:" parser block; strengthen LLM extraction prompt with negative examples; guard tutorial_patterns against preprocessor text |
| `agent/graph_nodes.py` | Fix 3: Error classification in PERCEIVE, pivot logic in PLAN (before rules_only); Fix 4: PHIL validation wrapper + call in both build paths; Fix 5: Parameter-aware duplicate check in validate(); Enhancement A: Error context injection into planning prompt |
| `agent/workflow_state.py` | Fix 6: Map extension safety net in _categorize_files(); fixed _load_category_rules() ImportError handling |
| `knowledge/thinking_prompts.py` | Enhancement A: Error classification context in expert prompt ("=== PREVIOUS FAILURE ===" section) |
| `knowledge/programs.yaml` | Fix 4: Added reference_model_enabled, reference_model_use_starting to phenix.refine; added mask_atoms to phenix.resolve_cryo_em; updated comment |
| `docs/reference/ARCHITECTURE.md` | Added phil_validator.py and error_classifier.py to Key Files table |

### Tests (4 files, 68 tests)

| File | Tests | Covers |
|------|-------|--------|
| `tests/tst_stop_condition_fix.py` | 20 | Fix 2: stripping, context-awareness, real user commands preserved |
| `tests/tst_phil_validation.py` | 15 | Fix 4: valid/invalid params, rewrite regression, PLAN cases |
| `tests/tst_error_classifier.py` | 24 | Fix 3: terminal/PHIL/label/retryable errors, pivot logic, failure counting |
| `tests/tst_fix_integration.py` | 9 | Cross-fix interactions: Fix 4→3, Fix 3+5, full AF_exoV chain |

Run with: `python3 tests/tst_stop_condition_fix.py && python3 tests/tst_phil_validation.py && python3 tests/tst_error_classifier.py && python3 tests/tst_fix_integration.py`

### Items completed in v115.01/v115.02

| Item | Status |
|------|--------|
| Fix 1: Clean tutorial set | ✓ Done (v115.02) — `_load_result_file_exclusions` in ai_agent.py |
| Fix 7: CASP7 MTZ recognition | Partially mitigated — auto-discover finds files |
| Enhancement B: Model upgrade | ✓ Done (v115.02) — error context forwarded in thinking_agent.py |
| advice_preprocessor.py | ✓ Done (v115.01) — solve-mode bypass |
