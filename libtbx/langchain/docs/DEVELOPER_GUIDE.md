# PHENIX AI Agent — Developer Guide

> Quick orientation: **§1** for the big picture, **§2** for how a cycle works,
> **§4** to add a program, **§5** to run the tests.
> Full architecture: [ARCHITECTURE.md](ARCHITECTURE.md)

---

---

## 1. Architecture at a Glance

The PHENIX AI Agent is a two-layer system for automated
structure determination:

```
┌────────────────────────────────────────────────┐
│ STRATEGIC PLANNER (v114, thinking_level=expert) │
│                                                │
│  Plan Generator → Gate Evaluator               │
│  Structure Model → Hypothesis Engine           │
│  Explanation Engine                             │
│        │                                       │
│        │ directives (prefer_programs,           │
│        │   after_program, strategy settings)    │
│        ▼                                       │
├────────────────────────────────────────────────┤
│ REACTIVE EXECUTION ENGINE (v112+)              │
│                                                │
│  LangGraph Pipeline:                           │
│    PERCEIVE → THINK → PLAN → BUILD →           │
│    VALIDATE → OUTPUT                           │
│                                                │
│  Workflow Engine (YAML-driven)                 │
│  Command Builder (programs.yaml)               │
│  Error Recovery (error_analyzer)               │
│  Safety Checks (sanity_checker)                │
│  Session Management (session.py)               │
└────────────────────────────────────────────────┘
```

**Key design principle:** The strategic planner is
additive — it communicates with the reactive engine
through the same directive system that a human user
would use. The reactive engine is unchanged by the
planning layer. Every safety check still applies.

When `thinking_level` is `advanced` (the default),
the planning layer is inactive and the reactive
engine operates alone.

### Import dependency rules

The codebase is split across directories with strict
import rules. Violating these breaks either the
standalone test environment or the PHENIX build:

```
agent/       → knowledge/    OK (shared code)
agent/       → libtbx.*      OK (cctbx core)
agent/       → phenix.*      NEVER
agent/       → phenix_ai/*   NEVER

knowledge/   → (no imports from agent/ or phenix/)

programs/    → agent/, libtbx.*, phenix.*    OK

phenix_ai/   → agent/, knowledge/           OK
phenix_ai/   → langchain/*, google.genai/*  OK

wxGUI2/      → agent/, libtbx.*, wx.*       OK

tests/       → agent/, knowledge/           OK
tests/       → programs/*, wxGUI2/*         With
               phenix_ai/*                  os.path.isfile
                                            guard
```

The `agent/` directory is **shared code** — it ships
with both the PHENIX client and the REST server. This
is why `agent/` modules use try/except ImportError
fallbacks for `libtbx.langchain` imports: they must
work both inside PHENIX and standalone.

### File inventory

| Directory | Files | Role |
|-----------|-------|------|
| `agent/` | ~100 .py | Shared client/server logic |
| `knowledge/` | ~15 .py + 10 .yaml | Domain knowledge, configs |
| `programs/` | 4 .py | PHENIX program driver |
| `phenix_ai/` | 14 .py | Server, agents, parsers |
| `wxGUI2/` | 6 .py | GUI panels |
| `tests/` | 55 .py | Test suites (~1100 tests) |
| `docs/` | 15+ .md | Documentation |

---

## 2. How a Cycle Works

This section traces a single cycle from start to
finish through the actual code path. Use this to
orient yourself when debugging or extending the agent.

### Session startup (once)

```
ai_agent.py: iterate_agent()
  │
  ├─ Load/create AgentSession
  ├─ Clean stale STOPWIZARD files
  ├─ Read thinking_level from PHIL
  │   _use_planning = (thinking_level == "expert")
  │
  ├─ Preprocess user advice
  │   _preprocess_user_advice() → LLM or passthrough
  │
  ├─ Extract directives from advice
  │   _extract_directives() → session.data["directives"]
  │
  ├─ [expert] Generate plan
  │   _initialize_plan() → session.data["plan"]
  │   Logs plan header with ✓/●/○ indicators
  │
  └─ Enter cycle loop: for cycle in 1..max_cycles
```

### Per-cycle flow

```
_run_single_cycle(cycle, session)
  │
  ├─ Check STOPWIZARD file (GUI abort)
  ├─ Check gate_stop flag (expert mode)
  │
  ├─ _query_agent_for_command()
  │   │
  │   ├─ [expert] Merge plan directives
  │   │   plan_to_directives() → prefer_programs,
  │   │   after_program, strategy settings
  │   │
  │   ├─ Build session_info dict
  │   │   Includes: history, files, directives,
  │   │   metrics, strategy_memory, structure_model
  │   │
  │   ├─ Call the graph (local or remote)
  │   │   run_ai_agent.run(request_json=...)
  │   │   │
  │   │   ├─ create_initial_state(...)
  │   │   │   Maps thinking_level="expert" →
  │   │   │   "advanced" for graph
  │   │   │
  │   │   └─ compiled_graph.invoke(state)
  │   │       │
  │   │       ├─ perceive(state)
  │   │       │   Workflow detection, valid programs,
  │   │       │   file categorization, loop detection
  │   │       │
  │   │       ├─ think(state)
  │   │       │   [advanced/expert] Structural
  │   │       │   validation, KB lookup, LLM
  │   │       │   reasoning, Structure Model update
  │   │       │
  │   │       ├─ plan(state)
  │   │       │   LLM selects program + strategy,
  │   │       │   or rules-based fallback
  │   │       │
  │   │       ├─ build(state)
  │   │       │   Command assembly from programs.yaml,
  │   │       │   file selection, parameter injection
  │   │       │
  │   │       ├─ validate(state)
  │   │       │   Directive validation, safety checks
  │   │       │
  │   │       └─ output_node(state)
  │   │           Package response for client
  │   │
  │   ├─ Extract command from response
  │   ├─ Persist structure_model, validation_history,
  │   │   strategy_memory from response
  │   │
  │   └─ Return (command, decision_info)
  │
  ├─ Record decision in session
  ├─ Send agent_cycle callback to GUI
  │
  ├─ Execute command
  │   _execute_command() → subprocess via easy_run
  │
  ├─ Handle result
  │   _handle_execution_result()
  │   Success: update session, clear advice_changed
  │   Failure: probe → auto-recovery → diagnosis
  │
  ├─ [expert] Gate evaluation
  │   GateEvaluator.evaluate() →
  │     advance / retreat / skip / stop / continue
  │   Send agent_gate_transition callback
  │
  ├─ [expert] Hypothesis evaluation
  │   evaluate_hypotheses() →
  │     confirmed / refuted / countdown
  │
  ├─ [expert] Cycle commentary
  │   generate_cycle_commentary() → log message
  │
  └─ Return: should_break?
```

### Session finalization (once)

```
_finalize_session(session)
  │
  ├─ [expert] Generate structure report
  │   _generate_structure_report()
  │   → structure_determination_report.txt
  │
  ├─ [expert] Write session summary JSON
  │   _write_session_summary_json()
  │   → session_summary.json
  │
  ├─ Generate AI summary (LLM or skip)
  │   _generate_ai_summary()
  │
  └─ Save session, set self.result
```

### Key files in the cycle path

| Step | File | Function |
|------|------|----------|
| Session setup | `programs/ai_agent.py` | `iterate_agent()` |
| Plan generation | `agent/plan_generator.py` | `generate_plan()` |
| State creation | `agent/graph_state.py` | `create_initial_state()` |
| Workflow detection | `agent/workflow_state.py` | `detect_workflow_state()` |
| Valid programs | `agent/workflow_engine.py` | `get_valid_programs()` |
| Expert reasoning | `agent/thinking_agent.py` | `expert_think()` |
| Program selection | `agent/graph_nodes.py` | `plan()` |
| Command assembly | `agent/command_builder.py` | `build_command()` |
| File tracking | `agent/best_files_tracker.py` | `BestFilesTracker` |
| Gate evaluation | `agent/gate_evaluator.py` | `GateEvaluator.evaluate()` |
| Session storage | `agent/session.py` | `AgentSession` |

---

## 3. Key Subsystems

### 3a. Workflow Engine

**Files:** `knowledge/workflows.yaml`,
`agent/workflow_engine.py`

The workflow engine determines which PHENIX programs
are valid to run at any given point. It reads
`workflows.yaml`, which defines a state machine for
each experiment type (X-ray and cryo-EM).

**How it works:**

1. `detect_workflow_state()` examines the session
   history and available files to determine the
   current workflow step (e.g., "analyze", "refine",
   "complete").

2. `get_valid_programs()` reads the programs listed
   for that phase in `workflows.yaml`, then applies
   filters:
   - Remove programs whose `done_tracking` flag is
     already set (prevents re-running)
   - Remove programs whose YAML conditions aren't met
     (e.g., "has: model" when no model exists)
   - Apply user directives: `skip_programs` removes
     programs, `prefer_programs` reorders them
   - Apply plan directives (Expert mode):
     `prefer_programs` from the current plan stage

**`after_program` semantics (v112.78, refactored v116.x):**
`after_program` in `stop_conditions` is a **minimum-run
guarantee**, not a hard stop.  When the named target
program has not yet run, it is added to the front of
`valid_programs` to ensure it runs.

Once the target program HAS run, behavior depends on
`stop_conditions.stop_after_requested`:

- `stop_after_requested=True` (user/README explicitly
  requested a stop-after) → `valid_programs` is wiped to
  `[STOP]`.  The user said stop; the workflow_engine
  enforces it.
- `stop_after_requested=False` (typically a plan-progression
  hint from `plan_to_directives`) → no action.  The workflow
  advances naturally to the next step / plan stage.

The flag is set by the directive extractor
(`agent/directive_extractor.py`) via `_is_stop_after_requested()`
when the raw user advice or README contains explicit stop-after
phrasing.  Plan-injected `after_program` does NOT set the flag.
`prompts_hybrid.py` also gates on the flag and emits "Stop
target: X" prompt directives only when it is True.

See `ARCHITECTURE.md` "Stop-after directive routing (v116.x)"
for the full discussion, including how this replaces the
previous "append-STOP" approach and why a `skip_validation`
carve-out was rejected.

3. The resulting list is passed to the LLM (or
   rules-based selector) which picks one.

**Stage transitions** happen when a program sets its
done flag and the phase's transition conditions are
met. For example, completing xtriage in the "analyze"
stage transitions to "obtain_model".

**Key configuration:**

```yaml
# workflows.yaml structure
xray:
  phases:
    analyze:
      programs:
        - program: phenix.xtriage
          conditions: [has: data]
      transitions:
        on_complete: obtain_model
    obtain_model:
      programs:
        - program: phenix.phaser
          conditions: [has: search_model, has: data]
        - program: phenix.predict_and_build
          conditions: [has: sequence]
      transitions:
        on_complete: refine
```

#### 3a.1. Stop-condition processing (v116.x)

Stop conditions flow through three layers:

1. **Detection** (`agent/directive_extractor.py`):
   `_is_stop_after_requested(advice)` returns True when the raw
   user advice contains explicit stop-after phrasing.  Positive
   patterns: `"stop after X"`, `"X and stop"`, `"only run X"`,
   `"just (do|run) X"`, `"stop when"`, `"stop once"`, `"stop if"`,
   `"stop at"`, `"Stop Condition: <real value>"`.  Negative
   patterns: `"don't stop"`, `"do not stop"`, `"never stop"`.
   Special-cased to False: `"Stop Condition: None"`, `"Stop
   Condition: not specified"`, heading-only `"Stop Conditions:"`,
   absence of any stop signal.

2. **Persistence** (`agent/directive_extractor.py`,
   `VALID_STOP_CONDITIONS`): `stop_after_requested: bool` is part
   of the directives schema and survives `validate_directives()`.
   The directive extractor sets it alongside `after_program` and
   `skip_validation` at every site where explicit user intent is
   detected (resolver `+stop` branches, intent="task",
   `tutorial_patterns`, `denmod_patterns`, `"stop after
   refinement"` regex).  It also pops the flag wherever
   `after_program` is removed (intent="solve" / "solve_constrained"
   carve-outs, ligand-workflow conflict, downstream-work conflict,
   skip-program conflict, resolver multi-action-no-stop branch).

3. **Consumption** (two consumers gate on the flag):
   - `agent/workflow_engine.py::_apply_directives` — when
     `after_program` is set AND target has run AND
     `stop_after_requested=True`, wipes `valid_programs` to
     `[STOP]`.  Otherwise (no flag, plan progression), takes no
     action.
   - `knowledge/prompts_hybrid.py::_format_directives_for_prompt`
     — emits "Stop target: X" prompt directives only when
     `stop_after_requested=True`.  Without this gating, plan-
     injected `after_program` would produce misleading prompts
     telling the LLM to stop mid-plan.  As a side effect, the
     `**Stop Conditions:**` section header is now emitted only
     when there is non-empty body content under it.

**Why this design (vs simpler alternatives):**

- A simpler `if after_program is set, treat as hard stop` was tried
  and reverted (v112.78): too many cases where `after_program` is a
  plan-progression hint, not a user stop.
- A pure min-run `if after_program is set, ensure it runs, never
  stop` was tried and reverted (v112.78 restoration patch): the
  workflow_engine no longer enforces user stops, leaving the LLM
  as the only mechanism — works most of the time but no structural
  guarantee.
- The flag distinguishes the two cases at the source, so the
  consumers can act unambiguously without subtle prompt-engineering
  workarounds or `skip_validation` carve-outs.

**Adding a new stop-condition trigger:**

If you want a new phrase to count as a user-explicit stop:

1. Add the pattern to `_POSITIVE_STOP_AFTER_PATTERNS` in
   `directive_extractor.py` (a compiled `re.Pattern`).
2. Add a test case to `tst_audit_fixes.py` covering both the
   positive case (sets the flag) and a negation case (doesn't).
3. No changes to consumers — they read the flag, they don't
   re-implement detection.

If you want to suppress a phrase that's currently triggering a
false positive:

1. Either add to `_NEGATIVE_STOP_PATTERNS` (if the phrase always
   means "don't stop") or to the `_STOP_CONDITION_NONE` regex (if
   it's a no-value marker like `"None"` / `"not specified"`).
2. Add a test case.

**Adding a new code path that emits `after_program`:**

If you add a new place where directives gain an `after_program`
based on user intent, you MUST also set
`directives["stop_conditions"]["stop_after_requested"] = True`
in the same code path.  If you add a code path that REMOVES
`after_program` (e.g. due to a new conflict), you MUST also
`pop("stop_after_requested", None)` from the same dict.  The
flag is a parallel signal to `after_program`; the two must stay
consistent.

If your code path emits `after_program` from a PLAN context
(not user intent), do NOT set the flag.  Plan-injected
`after_program` is by design a min-run hint only.

### 3b. Command Builder

**Files:** `knowledge/programs.yaml`,
`agent/command_builder.py`,
`agent/best_files_tracker.py`

The command builder translates a program name into an
executable PHENIX command with the correct files and
parameters.

**How it works:**

1. `programs.yaml` defines each program's inputs
   (required/optional files with extensions and
   flags), outputs, command template, and metric
   extraction rules.

2. `BestFilesTracker` maintains a catalog of all
   files produced during the session. When the command
   builder needs a "model" file, the tracker provides
   the best one (most recent, from the most reliable
   source).

3. `build_command()` resolves the template:
   - Selects files from BestFilesTracker
   - Applies strategy settings from directives
   - Handles special cases (crystal symmetry
     injection, resolution limits)
   - Content-based guards prevent misuse (e.g.,
     small-molecule PDB used as protein model)

4. `postprocess_command()` sanitizes the final
   command (removes invalid parameters, strips
   unverified resolution claims).

**programs.yaml key sections:**

```yaml
phenix.refine:
  inputs:
    required:
      model: {extensions: [.pdb], flag: ""}
      data:  {extensions: [.mtz], flag: ""}
    optional:
      ligand_cif: {extensions: [.cif], flag: ""}
  command: "phenix.refine {model} {data} {ligand_cif}"
  strategy_flags:
    ordered_solvent: "ordered_solvent=True"
    simulated_annealing: "simulated_annealing=True"
  log_parsing:
    r_free: {pattern: "R-free\\s*=\\s*([\\d.]+)"}
    r_work: {pattern: "R-work\\s*=\\s*([\\d.]+)"}
  done_tracking:
    flag: "refine_done"
    strategy: "count"
```

### 3c. THINK Node (Expert Reasoning)

**Files:** `agent/thinking_agent.py`,
`agent/thinking_prompts.py`,
`agent/validation_inspector.py`,
`knowledge/expert_knowledge_base_v2.yaml`,
`agent/kb_tags.py`, `agent/kb_loader.py`

The THINK node adds expert-level reasoning between
PERCEIVE and PLAN. It is controlled by
`thinking_level`:

| Level | Behavior |
|-------|----------|
| `none` | Pass-through (no THINK) |
| `basic` | Log analysis + LLM reasoning |
| `advanced` | Validation + KB + metadata + LLM |
| `expert` | Same as `advanced` in graph (planning is gated in ai_agent.py) |

**`advanced` mode adds:**

- **Structural validation** via
  `validation_inspector.py`: Ramachandran, rotamer
  outliers, clashscore, bonds/angles, model contents
  (chains, ligands, waters, ions), difference peaks
  (F2 with KD-tree), ligand RSCC (F1 with Z-score).

- **Expert knowledge base**: 56 entries in
  `expert_knowledge_base_v2.yaml`, matched by
  IDF-weighted tags against the current session
  context. Provides domain-specific guidance for
  common situations.

- **Structure Model** (`agent/structure_model.py`):
  A running understanding of the specific structure,
  updated each cycle from validation results. Tracks
  data characteristics, model state, progress
  trajectory, problems, hypotheses, and strategy
  blacklist.

- **File metadata tracking**: Content-aware queries
  like "find the best model" or "find the latest
  model with a ligand."

### 3d. Strategic Planner (Expert Mode)

**Files:** `knowledge/plan_templates.yaml`,
`knowledge/plan_template_loader.py`,
`agent/plan_generator.py`,
`agent/gate_evaluator.py`,
`agent/hypothesis_evaluator.py`,
`knowledge/explanation_prompts.py`,
`knowledge/plan_schema.py`

Activated by `thinking_level=expert`. The planning
layer generates a multi-phase plan at session start
and evaluates progress after each cycle.

**Plan templates** (`plan_templates.yaml`): 12
pre-defined plan skeletons covering MR, SAD, MR-SAD,
predict+build, cryo-EM, and variants with ligands,
twinning, and resolution extremes. Template selection
is deterministic (rule-based on experiment type,
available files, resolution, anomalous atoms).

**Template anatomy:**

```yaml
sad_phasing:
  description: "SAD/MAD phasing (X-ray)"
  applicable_when:
    experiment_type: xray       # string match
    has_anomalous_atoms: true   # required boolean

  priority: 30  # higher wins ties

  phases:
    - id: data_assessment
      programs: [phenix.xtriage]
      max_cycles: 1
      # SUCCESS CRITERIA: gate checks these after
      # each cycle. When ALL met → "advance".
      success:
        xtriage_completed: "true"

    - id: initial_refinement
      programs: [phenix.refine]
      max_cycles: 3
      success:
        r_free: "<0.30"        # advance when met
      # GATE CONDITIONS: checked when success NOT
      # met. Trigger retreats or stops.
      gate:
        - if: "r_free > 0.45 after 2 cycles"
          action: "retreat_to model_rebuilding"
```

Key distinction: `success` criteria tell the gate
when to **advance** (things going well). `gate`
conditions tell it when to **retreat or stop** (things
going badly). Success is checked first.

**Gate evaluator** (`gate_evaluator.py`): Evaluates
after each cycle. The evaluation order is:

1. Success with hysteresis (1.5% buffer prevents
   oscillation around threshold)
2. Phase exhaustion (max_cycles reached)
3. Gate conditions (retreat triggers)
4. Anti-oscillation safeguards:
   - Strategy blacklist (don't repeat failed
     approaches)
   - Retreat counter (max 3 retreats)
   - Monotonic progress gate
   - Cooldown (2 cycles between retreats)
   - Depth limit (max 2 levels of retreat)

**Hypothesis engine** (`hypothesis_evaluator.py`):
Single active hypothesis budget with verification
latency (wait N cycles for test to complete before
evaluating). Re-validates confirmed hypotheses if
evidence decays.

**Explanation engine** (`explanation_prompts.py`):
Three levels of commentary:
- `generate_cycle_commentary()`: template-based,
  no LLM, every cycle
- `generate_stage_summary()`: template-based,
  no LLM, at phase transitions
- `generate_final_report()` /
  `generate_stopped_report()`: template-based, at
  session end (data from Structure Model, no LLM)

The only LLM-generated report is the failure diagnosis
(`ai_failure_diagnosis.html`), produced by
`failure_diagnoser.py` when a terminal error is
detected.

### 3e. Model Placement Gate (v114.1, hardened v114.2)

**Files:** `programs/ai_agent.py`
(`_detect_model_placement`,
`_skip_plan_stages_for_placement`),
`agent/workflow_engine.py`
(`model_is_placed_confirmed` in context),
`agent/workflow_state.py` (symmetry mismatch
detection)

Prevents destructive program selection on models
that already fit the data. This is the fix for the
class of bugs where LLM modes underperform
rules_only by running phaser/autosol on pre-solved
structures.

**Data flow:**

```
_detect_model_placement (post-result)
  → session.data["model_is_placed"] = True
    → session_info["model_is_placed"]
      → api_client.build_session_state
        → run_ai_agent.create_session_info
          → workflow_engine context
            → get_valid_programs filters out
              phaser, autosol, predict_and_build
```

**Detection sources** (v114.2 additions in bold):

- model_vs_data: CC > 0.3 in metrics dict
- **model_vs_data: R-free < 0.50 parsed from result text**
  (fallback when metrics dict is empty)
- **model_vs_data: R-work < 0.45 parsed from result text**
- refine: R-free < 0.50 or map_cc > 0.3
- **RSR symmetry mismatch → `needs_dock`** (failed
  real_space_refine with "dimensions mismatch" sets
  `placement_probe_result` in workflow_state)

**Plan stage skipping:** When placement is detected,
pending `molecular_replacement` and
`experimental_phasing` stages are marked "skipped".
**Additionally** (v114.2), when the placement evidence
shows R-free < 0.35 or R-work < 0.30 or CC > 0.5,
`model_rebuilding` is also skipped (prevents
autobuild from damaging a well-refined model).

**Dock keyword handling** (v114.2): Words like "dock",
"docking", "dock_in_map", "place the model" in user
advice (including README text) cancel the "refine" →
"model is placed" assumption that suppresses the
placement probe. Both `workflow_engine.py` and
`plan_generator.py` share this keyword list.

**HETATM detection:** At startup,
`_scan_input_models_for_ligands` scans input PDB
files for non-water HETATM records, sets
`session.data["input_has_ligand"]`. This enables
polder without a ligandfit step when the model
already contains ligands.

### 3f. Plan Enforcement (v114.2)

**Files:** `agent/graph_nodes.py` (plan_node,
build), `phenix/programs/ai_agent.py`
(`_plan_has_pending_stages`),
`agent/contract.py` (v4 field)

Prevents premature STOP when the strategic plan
has pending stages.  Without this, both LLMs
(Gemini and OpenAI) frequently stop after the first
stage completes (e.g., after resolve_cryo_em in a
4-stage cryo-EM plan).

**Two enforcement points:**

1. **AUTO-STOP suppression** (plan_node): When
   `metrics_trend.should_stop` fires but
   `session_info.plan_has_pending_stages` is True,
   `should_stop` is cleared and the LLM/rules
   selector runs normally.

2. **LLM STOP override** (build): When the LLM
   returns `program=STOP` but
   `plan_has_pending_stages` is True, the node
   selects the best available program from
   `prefer_programs` in the plan's current stage
   directives.

**Contract field:** `plan_has_pending_stages`
(bool, default False, v4).

### 3f-b. Dotted PHIL Passthrough (v114.2)

**Files:** `agent/program_registry.py`,
`agent/command_postprocessor.py`

The sanitizer (Rule D) intentionally passes dotted
PHIL paths through for runtime validation. But
`program_registry.build_command` used to drop them
as "unknown dotted keys". This broke every valid
PHIL parameter the LLM generated that wasn't
pre-registered in `strategy_flags`.

**Fix:** Dotted PHIL keys now pass through to the
command line with a PASSTHROUGH log. PHENIX's PHIL
interpreter validates at runtime; wrong paths produce
clear errors.

**`strict_strategy_flags: true`:** Programs that crash
on unexpected PHIL parameters (e.g., `process_predicted_model`
does not accept `crystal_symmetry.space_group`) can set
this flag in `programs.yaml`. When set, only defined
`strategy_flags` are allowed — `KNOWN_PHIL_SHORT_NAMES`
(space_group, unit_cell, nproc, etc.) are NOT passed
through, and dotted PHIL paths are skipped.

**Rewrite system:** `_PARAM_REWRITES` in
`command_postprocessor.py` expands short names to
full PHIL paths (e.g., `mask_atoms` →
`strategy.mask_atoms` for resolve_cryo_em).
`_STRATEGY_REWRITES` in `graph_nodes.py` maps
LLM's dotted keys to canonical strategy_flag names
before the builder processes them.

### 3f-c. Client-Server Contract

**Files:** `agent/contract.py`,
DEVELOPER_GUIDE.md §8 (Backward Compatibility & Contract)

The agent has a client-server architecture: the
client (`programs/ai_agent.py`) runs on the user's
machine, and the agent's "brain" (the LangGraph
pipeline) may run on a remote server. When the server
is updated, old clients still connect to it.

**The contract** (`agent/contract.py`) is the single
source of truth for the client-server interface. It
defines every `session_info` field, its default value,
and the protocol version it was introduced in.

**Rules for maintaining compatibility:**

- Never read `session_info["key"]` — always use
  `session_info.get("key", default)`.
- Every new field must be registered in
  `contract.py` with a default value.
- The `normalize_session_info()` function fills in
  missing fields with defaults before the graph runs.
- Protocol version is checked at the PERCEIVE node;
  old clients with missing fields are handled
  gracefully.

**To add a new session_info field:**

1. Add it to the registry in `contract.py`
2. Provide a sensible default value
3. Bump the protocol version if the field is required
   (old clients won't send it)
4. Add a test in `tests/tst_contract_compliance.py`
   verifying the field is registered

See DEVELOPER_GUIDE.md §8 (Backward Compatibility & Contract) for the
full contract system.

**Critical: the transport triple.** Every new
`session_info` field that must survive across cycles
also requires updates to three whitelists — missing
any one silently drops the value between cycles:

1. `build_session_state()` return dict in
   `agent/api_client.py`
2. `build_request_v2()` normalization allowlist in
   `agent/api_client.py`
3. `session_state → session_info` mapping in
   `phenix_ai/run_ai_agent.py`

`normalize_session_info()` in `contract.py` is safe —
it only fills defaults for explicitly defined fields;
unknown fields pass through unchanged, so a field that
survives transport won't be zeroed out by normalization.

Historical additions (each required all three sites):
`strategy_memory` (v113), `bad_inject_params`,
`unplaced_model_cell`, `explicit_program`,
`session_blocked_programs`, `asu_copies` (v115).

### 3f-d. ASU Copy Count Tracking (v115)

**Files:** `agent/directive_extractor.py`,
`agent/graph_nodes.py`, `programs/ai_agent.py`,
`phenix_ai/run_ai_agent.py`, `agent/api_client.py`

Tracks the number of copies of the search model in
the asymmetric unit (ASU) and passes this to Phaser
as `search_copies` each cycle.

**Data flow — directive path:**
```
User advice: "4 copies in the ASU"
  → _extract_copies_from_directives()
    → session.data["asu_copies"] = 4  (all 3 extraction paths)
      → session_info["asu_copies"] = 4  (each cycle, ai_agent.py)
        → api_client.build_session_state  → build_request_v2
          → run_ai_agent session_state["asu_copies"]
            → BUILD: strategy["component_copies"] = 4
              → phaser.search_copies=4
```

**Data flow — xtriage path:**
```
xtriage log: "Best guess : 4 copies in the ASU"
  → _fallback_extract_metrics(): log_analysis["n_copies"] = 4
    → _attach_thinking_metadata() 3-tier fallback:
        session_info → log_analysis["n_copies"] → raw log_text regex
          → metadata["asu_copies"] = 4
            → history_record["asu_copies"] = 4
              → session.data["asu_copies"] = 4  (only if not set by directive)
                → next cycle: session_info → BUILD → phaser.search_copies=4
```

**Priority:** Directive always wins — directive path
sets `session.data["asu_copies"]` unconditionally;
xtriage path is guarded by `if not session.data.get("asu_copies")`.

**BUILD guard:** Injection is skipped if the LLM
already set `component_copies` in the strategy dict
(prevents overwriting an explicit LLM decision).

**Sanity bound:** All extraction points enforce
`1 ≤ n ≤ 30`; values outside this range are ignored.

**Transport:** `asu_copies` is included in both
whitelists in `api_client.py` (`build_session_state`
and `build_request_v2`), so it survives the
client→server round-trip each cycle.

**v115.05 same-cycle fix:** The BUILD node now also
reads `state["log_analysis"]["n_copies"]` directly,
so when Phaser runs and reports the copy count in the
same cycle, the copies are available immediately —
eliminating the one-cycle delay where the first
post-Phaser autobuild would miss the copies count.

### 3f. Error Recovery

**Files:** `agent/error_analyzer.py`, `agent/error_classifier.py`,
`agent/failure_diagnoser.py`, `programs/ai_agent.py`

When a PHENIX program fails, the agent attempts
automatic recovery through three systems that execute
in sequence:

1. **ErrorAnalyzer** (YAML-driven, `recoverable_errors.yaml`):
   Detects errors the agent can auto-fix. Three resolution
   kinds are currently supported:
   - **`select_value`**: disambiguates an enum the LLM
     left ambiguous (e.g., `ambiguous_data_labels` for
     multi-array MTZ).  Picks a default per a YAML-defined
     preference policy.
   - **`add_parameter`**: injects a missing flag (e.g.,
     `xray_data.r_free_flags.generate=True` when the user's
     MTZ lacks R-free flags).
   - **`strip_parameter`** (v119.H17): removes an
     inappropriate flag entirely on retry.  Used when the
     LLM passed a flag whose value type doesn't match what
     the program accepts (e.g., raw `data.mtz` as
     autobuild's `map_file=` parameter, which requires
     PHIB/FOM phase columns).  Executor support in
     `programs/ai_agent.py` (H17.1) applies a robust regex
     to strip the flag from the retry command before
     execution, preserving PHIL spacing, quoted-with-spaces,
     and single-quoted values.

   ErrorAnalyzer produces file-level recovery flags and
   forces a retry.  See `tests/tst_error_analyzer.py` and
   `tests/tst_h17_strip_executor.py` for the contract.

2. **DiagnosisDetector** (YAML-driven, `diagnosable_errors.yaml`):
   Detects terminal errors needing human intervention
   (crystal symmetry mismatch, model outside map, SHELX
   not installed, unknown PHIL parameter, polymer on
   special position). Stops the run and produces an
   LLM-generated HTML diagnosis page.

3. **Error classifier** (`error_classifier.py`, v115):
   Classifies failures for the graph's THINK and PLAN
   nodes. Five categories: TERMINAL, PHIL_ERROR,
   AMBIGUOUS_PHIL, LABEL_ERROR, RETRYABLE. Extracts
   bad parameter names and provides context for the
   `should_pivot()` function, which excludes the failed
   program from valid_programs on the next cycle.

4. **`_classify_error`** (hardcoded, `ai_agent.py`):
   The oldest classifier. Maps errors to INPUT_ERROR
   (agent's fault, don't count) or REAL_FAILURE (real
   problem). Used only for cycle history bookkeeping.

**Known design tension:** These three systems evolved
independently and have overlapping patterns. See
ARCHITECTURE.md "Design tensions" section for a
detailed analysis and recommended consolidation path.

The **catch-all injection blacklist** prevents the
agent from injecting the same recovery parameter
twice (which would cause an infinite retry loop).

### 3g. Session Management

**Files:** `agent/session.py`

`AgentSession` manages all persistent state:

- Cycle history (program, command, result, metrics)
- Directives (extracted from user advice)
- Structure Model (v114)
- Validation History (v114)
- Plan data (v114, expert mode)
- Best files catalog
- Strategy memory

State is serialized to `session.json` after every
cycle. On resume, the session is restored from this
file. The `advice_changed` flag triggers re-evaluation
of the plan and directives when the user provides
new advice.

### 3h. GUI Integration

**Files:** `wxGUI2/Programs/AIAgent.py`

The GUI communicates with the agent through a callback
system. The agent (running in a child process) sends
callbacks via `libtbx.call_back`; the GUI
(`RunNotebook`) intercepts them in the main thread.

**Callback types:**

| Message | When | Data |
|---------|------|------|
| `agent_status` | Session start | Status text |
| `agent_directives` | After extraction | Directives, files, advice |
| `agent_plan` | Plan generated | Plan dict, header text |
| `agent_cycle` | Each decision/result | Cycle, program, reasoning, plan context |
| `agent_gate_transition` | Stage change | Action, stages, reason |
| `sub_job_complete` | Program finished | Job info for project DB |

**Display panels:**

- `ConfigPanel`: Setup (files, advice, settings).
  Includes tutorial banner (`_detect_tutorial`,
  `_display_tutorial_banner`) when working directory
  is a Phenix tutorial.
- `RunNotebook` → `_progress_panel`: Live progress
  (text append only). Stage headers appear once per
  stage (tracked via `_last_stage_id`). Cycle results
  display inline as `Running prog ... [OK]`.
- `ResultPanel`: Structured results after completion.

**Tutorial detection** (v114.2): `_detect_tutorial()`
loads `tutorials.phil` from `phenix_examples`, matches
CWD basename against tutorial names (longest first to
avoid prefix collisions), requires README file to
confirm. `_find_unsupported_programs()` scans the
README for `phenix.xxx` tokens not in the agent's
program registry.

All new callback fields must use
`getattr(data, 'field', default)` on the GUI side
for backward compatibility. All callback sends must
be wrapped in `try/except Exception: pass` (non-
critical, must not crash the agent).

> **YAML is the source of truth for all subsystem behavior.** The workflow
> engine, command builder, and program registry are all driven by YAML files
> in `knowledge/`. When reading the subsystem descriptions above:
>
> - To understand what a YAML field does → §13 Configuration Reference
> - To add or modify a program's behavior → §4 Adding a New Program
> - To edit YAML safely → see CCTBX_LLM_PROGRAMMING_GUIDELINES.md §5


---

### 3i. Stop-after Routing and Raw-Advice Extraction (v117 / v117.1 / v117.2 / v117.3)

The v117 family introduced two related changes to directive
extraction and stop-after handling.  This subsection covers what a
new developer needs to know to safely modify either area.

#### The two kinds of `after_program`

`stop_conditions.after_program` carries two distinct intents
depending on whether `stop_conditions.stop_after_requested` is set:

```
stop_after_requested=True   → user-explicit stop
                              wipe valid_programs to ["STOP"] when done
                              prompt emits "Stop target: X" 4-line block

stop_after_requested=False  → plan-progression hint
or absent                     no wipe; plan progresses naturally
                              prompt emits nothing for after_program
```

The flag is set in two ways:

1. **Directive extractor regex helper** `_is_stop_after_requested`
   detects 11 positive patterns: `"X and stop"`, `"stop after X"`,
   `"only run X"`, `"just run X"`, `"Stop Condition: <real value>"`,
   the period-stop pattern `"\. stop"`, and others.  Run automatically
   on the user advice; flag set if any pattern matches.

2. **LLM directly** — the JSON schema documents `stop_after_requested`
   as a boolean field.  The dual-input prompt (see below) instructs
   the LLM to set it based on the raw user instruction.

The downstream code that respects this flag:
- `workflow_engine._apply_directives` — gates the wipe
- `prompts_hybrid._format_directives_for_prompt` — gates the 4-line block

#### Raw-advice dual-input extraction

`extract_directives()` accepts an optional `raw_advice` parameter.
When supplied and different from `user_advice`, the extractor uses
`DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` (dual-input + AUTHORITY
paragraph) instead of the single-input prompt.  When equal or None,
behavior is unchanged (backward compat).

The AUTHORITY paragraph tells the LLM that raw is the source of truth
for intent (after_program, start_with_program, stop_after_requested)
and processed fills gaps for files, parameters, and experiment type.

Plumbing path: client sets `session.data["raw_advice"]` →
`programs/ai_agent.py::_extract_directives` reads it →
`directive_params.ai_analysis.user_advice_raw` →
`programs/ai_analysis.py` PHIL parameter →
`phenix_ai/run_ai_analysis.py::run_directive_extraction` →
`extract_directives(..., raw_advice=...)`.

#### The grounding ∧ flag interaction (v117.1)

The v116.19a `_validate_after_program_grounded` guardrail drops
`after_program` when the program name doesn't appear literally in the
advice — protecting against fabrication.  But when v117 Step 1's
AUTHORITY paragraph correctly leads the LLM from a verb in the raw
advice (e.g. "density modify") to a learned mapping (e.g.
`phenix.resolve_cryo_em`), the literal-substring check fires
incorrectly.

v117.1 resolves this: when the LLM sets BOTH `after_program` AND
`stop_after_requested=True`, the literal grounding check is skipped.
The flag is the LLM's explicit user-stop assertion — strictly
stronger evidence than the heuristic check.

**Rule for modifying this area:** before changing either the
grounding guardrail or the flag handling, ensure all 5 cases in
`tests/tst_grounding_stop_after_requested.py` (K1-K5) still pass.
They lock the interaction at the boundary.

#### Filling `after_program` from raw advice (v117.2)

After v117.1 fixed the grounding guardrail, a separate gap remained:
the LLM (specifically openai, observed on the C1 test) occasionally
emits `stop_after_requested=True` but omits `after_program` for
ambiguous inputs.  Without `after_program`, the gated wipe in
`workflow_engine._apply_directives` is unreachable (the wipe path is
gated on `if after_program:`).

v117.2 closes the gap in `extract_directives()`: when the flag is set
but the program field is unset, call `_resolve_after_program()` on
the raw advice to fill in `after_program` via `_ACTION_TABLE`.

The fix's three guards:

```python
if stop_after_requested is True and not after_program:
    if raw_advice and _is_stop_after_requested(raw_advice):
        _resolve_after_program(directives, raw_advice.lower())
```

- The `not after_program` guard means the fix never overrides an LLM
  choice; it only fills in when the LLM omitted the field.
- The `_is_stop_after_requested(raw_advice)` guard means a
  hallucinated flag (LLM sets True against advice with no stop
  signal) does not get amplified — the resolver isn't called.
- The `raw_advice or user_advice` fallback maintains backward
  compatibility with callers that pass only `user_advice` (single-
  input path; the fallback gives them the same benefit).

The fix lives at line ~700 in `agent/directive_extractor.py`,
between `_apply_crystal_symmetry_fallback` and the ollama-empty
retry, after grounding has already validated/dropped any LLM-set
`after_program`.

**Rule for modifying this area**: before changing the fill-in
logic, ensure all 6 cases in
`tests/tst_after_program_fill_from_raw.py` (K1-K6) still pass.
K1 (the C1 fix) and K3 (no stop intent in raw → no fill) are the
most important — they together define the lock between "do useful
work" and "don't amplify bad data."

#### Extended stop-intent phrasings (v117.3)

After v117.2 fixed the C1 "missing after_program" gap, openai's
extractor still failed on the pre-existing
`tst_directive_extraction.py::explicit_stop_after_phaser` test
(0/5).  The fixture's user-side phrasing in the Special Instructions
section is `"Stop the workflow immediately after phaser completes"`
— a clear user-stop directive that:

- v117.2's regex helpers don't recognize (no pattern matches it)
- The LLM (openai) doesn't translate into `stop_after_requested=True`
  (the schema docs gave no example of this phrasing class)

So openai produced `after_program=phenix.phaser` without the flag,
v117.1's exemption couldn't fire, and Failure 2 in
`_validate_after_program_grounded` dropped the directive.

v117.3 closes the gap with three small additive changes:

1. **`_IMPERATIVE_STOP_MARKERS`** (in `agent/directive_extractor.py`
   around line 1283) gains 5 new entries:
   `"stop the workflow"`, `"immediately after"`,
   `"after it completes"`, `"after it finishes"`,
   `"is the last step"`.

2. **`_POSITIVE_STOP_AFTER_PATTERNS`** (around line 2812) gains 2
   new regex patterns:
   `r'\bstop\s+the\s+workflow\b'` and `r'\bis\s+the\s+last\s+step\b'`.

3. **`DIRECTIVE_EXTRACTION_PROMPT`** schema docs for
   `stop_after_requested` (around line 244) gain 4 new recognized
   phrasings and 2 program-neutral few-shot examples.

**Why some markers are NOT in the regex list**: this is a deliberate
asymmetry.  Marker matching is window-bounded (300 chars near program
name), so `"immediately after"` is safe there but unsafe as a global
regex (would fire on "the workflow continues immediately after
prediction").  The over-permissive `\bafter\s+\w+\s+completes?\b`
and `\bafter\s+\w+\s+(?:finishes|is\s+done)\b` patterns were
considered for the regex list and **dropped** during planning after
K5 verified the false positive on "validation runs after refinement
completes."

**Rule for modifying this area**: before changing any of the three
v117.3 additions, ensure all 9 cases in
`tests/tst_extended_stop_phrasings.py` (K1, K2, K4, K5, K6, K7a/b,
K8a/b) still pass.  K4 (AF_7mjs preservation) and K5 (descriptive
false-positive guard) are the most important: they lock the
preservation properties that the rejected over-permissive patterns
would have broken.

**End-to-end flow with v117.3 applied** for the
`explicit_stop_after_phaser` failure mode:

- LLM emits `after_program=phenix.phaser` (correct) but no
  `stop_after_requested` flag (the openai gap)
- Grounding guardrail's `_imperative_marker_nearby` fires because
  the new marker `"stop the workflow"` matches within 300 chars of
  the phaser occurrence → `after_program` preserved
- Downstream `_apply_workflow_intent_fallback` calls
  `_resolve_after_program`; `_is_stop_after_requested` on the
  stripped advice returns True because of the new regex pattern →
  `stop_after_requested=True` set automatically

Result: full correct directives reach `workflow_engine._apply_directives`.

#### Modifying these areas — checklist

- [ ] If changing `_is_stop_after_requested`, also check the 11
      positive patterns vs the unit tests in
      `tst_extract_raw_advice.py::test_is_stop_after_requested_helper_basic_cases`
- [ ] If changing `_apply_directives` after_program logic, verify
      M1-M7 in `tst_validate_step_after_program_guard.py`
- [ ] If changing `_format_directives_for_prompt` after_program
      handling, verify K1-K3 + N1-N5 in `tst_after_program_prompt.py`.
      N5 is an anti-regression test for the old "CRITICAL: You MUST
      run X" wording — keep the new wording.
- [ ] If changing `_validate_after_program_grounded`, verify K1-K5
      in `tst_grounding_stop_after_requested.py`.  Specifically, K2
      (AF_7mjs case still dropped without flag) and K1 (C1 case
      preserved with flag).
- [ ] If changing the v117.2 fill-in logic in `extract_directives`,
      verify K1-K6 in `tst_after_program_fill_from_raw.py`.  K1 (fix
      fires) and K3 (guard against hallucinated flag) are the
      essential pair.
- [ ] If changing v117.3 markers or patterns in
      `_IMPERATIVE_STOP_MARKERS` / `_POSITIVE_STOP_AFTER_PATTERNS`,
      verify all 9 cases in `tst_extended_stop_phrasings.py`.  K4
      (AF_7mjs preservation) and K5 (descriptive false-positive
      guard) are essential — they prevent reintroduction of the
      over-permissive patterns rejected during plan review.
- [ ] If changing the extractor LLM prompt — both
      `DIRECTIVE_EXTRACTION_PROMPT` and
      `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` — verify the schema
      documentation for `stop_after_requested` and the AUTHORITY
      paragraph remain consistent.

#### Anti-patterns to avoid

- **Re-introducing v112.78's "wipe everywhere".**  Wiping
  `valid_programs` to `["STOP"]` whenever `after_program_done` would
  regress AF_7mjs and any future per-stage min-run hints.  The wipe
  must be gated on `stop_after_requested=True`.
- **Tightening grounding to drop after_program when the flag is
  True.**  This would re-resurface the C1 failure.  The flag is the
  grounding signal in that case.
- **Adding `after_program` patterns that don't also pop
  `stop_after_requested`** (or vice versa).  Pattern-based
  assignments are tentative — the intent classifier may override
  them.  All pop sites must clear both.  See
  `agent/directive_extractor.py`: 8 pop sites including
  `_set_by_pattern` cleanup.
- **Reading `stop_after_requested` from anywhere except
  `stop_conditions`.**  The flag lives only there; downstream readers
  should always destructure from `directives["stop_conditions"]`.
- **Removing the `_is_stop_after_requested(raw_advice)` guard in
  v117.2's fill-in code.**  This guard is what prevents a
  hallucinated `stop_after_requested=True` from getting amplified
  into a hallucinated `after_program`.  If the LLM sets the flag
  against advice that genuinely has no stop intent, the resolver
  must not be invoked — leave the bad state alone rather than
  compound it.
- **Reintroducing the rejected `\bafter\s+\w+\s+completes?\b`-class
  patterns into `_POSITIVE_STOP_AFTER_PATTERNS`.**  These were considered
  during v117.3 planning and dropped after K5 verified the false
  positive on descriptive prose like "validation runs after refinement
  completes."  If the failure mode reappears in a new LLM test
  fixture, address it through `_IMPERATIVE_STOP_MARKERS` (window-
  bounded) instead, or via prompt schema docs.  The bare regex
  patterns are too permissive.

#### Why the docs were updated late

v117 originally shipped without full doc reconciliation (see
`next_steps.md` §4.3) because the stop_refactor delivery's docs were
built against a baseline that predated some v116.14+ content.
Wholesale replacement would have regressed that content.  This
section was hand-merged in v117.1.

### 3j. PHIL Parameter Authoring

When adding a new PHIL parameter to a program (`programs/<file>.py`),
the parameter block goes inside a PHIL scope and looks like:

```
my_param = None
  .type = str
  .short_caption = Brief one-line caption
  .help = Longer description that may span \
          multiple lines using backslash-newline continuation.
  .style = hidden
```

**Syntax constraints that may not be obvious:**

- **No semicolons in `.help` strings.**  The PHIL parser does not
  accept semicolons inside `.help` values.  A `.help` value like
  `"Optional; when set the server uses it"` will cause the parser
  to fail or misread the field.  Use periods, commas, or em-dashes
  instead.  Observed during v117 Step 1 when the original `.help`
  for `user_advice_raw` used a semicolon and broke PHIL parsing.
- **Hyphens in `.help` strings are risky too.**  Hyphenated words
  like `backward-compatible` or `multi-step` can cause the PHIL
  parser to interpret the hyphenated token as a malformed
  definition name and crash with
  `RuntimeError: Syntax error: improper definition name`.
  The failure mode depends on tokenizer context: hyphens inside
  parentheses are often tolerated, but hyphens after sentence-
  ending punctuation are not.  Observed in v118 when
  `user_advice_raw`'s `.help` text contained `backward-compatible`
  after `Optional;` — and notably, those two rule violations
  were on the same line, so it's not always obvious which one is
  the actual trigger.  When the error message names a hyphenated
  token, check the surrounding text for a semicolon that may have
  prematurely terminated the `.help` value and put the tokenizer
  in the wrong state.  Safest practice: write `backward compatible`
  or `multi step` (space-separated), or keep the help to a single
  line where the entire value-after-`=` is taken verbatim.
- **Backslash-continuation for multi-line `.help`.**  PHIL strings
  end at the next unescaped newline.  Use `\` at the end of each
  continued line.
- **Caption is one line.**  `.short_caption` cannot span lines and
  cannot contain newlines or backslash-continuations.
- **`.style = hidden` keeps the parameter out of GUI auto-rendering.**
  Use for server-internal params (e.g., the `user_advice_raw` field
  added in v117 Step 1) that aren't meant to be set by GUI users.
- **PHIL parameter names are lowercase with underscores.**  No
  camelCase, no hyphens.  Sub-scopes use dotted notation
  (`my_scope.my_param`).

**Mirroring PHIL on client and server:**

For parameters that flow client → server through PHIL serialization
(see §3i for the v117 raw-advice path as an example), the PHIL
parameter MUST be defined identically in both
`programs/ai_agent.py` (client) and `programs/ai_analysis.py`
(server).

**These are independent schemas.**  `ai_agent.py` parses its OWN
`master_params` string (defined at line ~144) when the user runs
`phenix.ai_agent`.  `ai_analysis.py` parses its OWN `master_params`
when the analysis server is invoked.  Adding a parameter to ONE
does NOT propagate to the other — they're parsed independently.

**Failure mode if you forget to mirror:**

The client does `directive_params = copy.deepcopy(self.params)`
and then assigns a new attribute.  If that attribute is NOT
declared in the CLIENT's master_params, the assignment raises:

```
AttributeError: Assignment to non-existing attribute
  "ai_analysis.your_new_field"
```

In production this is caught by the surrounding try/except and
silently swallowed — directive extraction returns `{}`, the
agent ignores the parameter, and behavior drifts in ways that
may take a production run to notice.  This was the H18 → H18.1
production failure mode.

**Defensive test pattern** (`tests/tst_h18_1_phil_roundtrip.py`):

When adding a new PHIL parameter, add a K-test that exercises
the full parse + extract + deep-copy + assign path:

```python
def test_phil_roundtrip_for_new_param():
    # Source-grep variant — runs in sandbox without libtbx
    with open("programs/ai_agent.py") as f:
        content = f.read()
    start = content.find('master_params = """')
    end = content.find('master_phil = libtbx.phil.parse')
    master_region = content[start:end]
    assert "your_new_field" in master_region

    # Live-PHIL variant — needs PHENIX, skips gracefully
    try:
        import libtbx.phil, copy
    except ImportError:
        return
    master_phil = libtbx.phil.parse(master_region.split('"""', 1)[1])
    params = master_phil.extract()
    directive_params = copy.deepcopy(params)
    directive_params.ai_analysis.your_new_field = "test"
    assert directive_params.ai_analysis.your_new_field == "test"
```

This catches the entire class of "I forgot to declare the param
in master_params" bugs at sandbox time.

**Checklist before shipping a new PHIL parameter that flows
client → server:**

1. Add declaration to `programs/ai_agent.py:master_params` AND
   `programs/ai_analysis.py:master_params`.
2. Add a `.help` comment in each declaration noting the
   cross-file mirror requirement.
3. Add a `tst_*_phil_roundtrip.py` style K-test that exercises
   the PHIL parse → assign path for the new field.
4. Run sandbox tests; the source-grep check should pass.
5. Deploy to PHENIX and run the live-PHIL test; should also pass.

When adding such a parameter, add it to both files in the same
commit, verify the parse on both sides with
`libtbx.python -c 'from programs.ai_agent import master_phil_str; ...'`,
and run the cross-tree tests (`tst_audit_fixes.py` exercises the
end-to-end PHIL flow when run in the PHENIX tree).

### Optional dependency handling

`libtbx.langchain` integrates several optional Python packages
(`langchain_openai`, `langchain_anthropic`, `langchain_google_genai`,
`langchain_chroma`, `chromadb`, etc.).  These packages can fail to
import in user environments for many reasons beyond "not installed":

1. **Not installed** — `ImportError`.
2. **Version conflicts among transitive deps** — surface as
   `TypeError`, `RuntimeError`, etc. raised from inside the
   package's own import-time code.  (Observed in v118 when
   `chromadb` → `opentelemetry` → `protobuf` raised
   `TypeError: Descriptors cannot be created directly` because
   the protobuf version was too new for chromadb's bundled
   `_pb2.py` files.)
3. **Missing native binaries** — `OSError`.
4. **API breakage between versions** — `AttributeError` at
   import time.

To keep cctbx modules importable in any environment, use this
pattern for optional dependencies.

**Don't (eager at module top):**
```python
from langchain_chroma import Chroma  # at module top

def use_chroma(embeddings):
    return Chroma(embedding_function=embeddings)
```

This fails to import the entire module — including unrelated
helper functions — in any env where the dep is unusable.

**Do (lazy with probe + shared helper when 2+ files share a dep):**

Create a private helper module (e.g.
`rag/_chroma_resilience.py`):
```python
_PROBE_STATE = None
_IMPORT_ERROR = None
_MODULE = None
_CLASS = None


def ensure_chroma():
    global _PROBE_STATE, _IMPORT_ERROR, _MODULE, _CLASS
    if _PROBE_STATE is not None:
        return _MODULE, _CLASS
    try:
        import chromadb as _chromadb_mod
        from langchain_chroma import Chroma as _Chroma
        _MODULE = _chromadb_mod
        _CLASS = _Chroma
        _PROBE_STATE = True
    except Exception as e:        # NOT just ImportError
        _IMPORT_ERROR = "%s: %s" % (type(e).__name__, e)
        _PROBE_STATE = False
    return _MODULE, _CLASS


def is_chroma_available():
    return ensure_chroma()[0] is not None


def chroma_unavailable_error():
    err = _IMPORT_ERROR or "unknown"
    return RuntimeError(
        "chromadb is required for this feature but not usable: %s\n"
        "Install: pip install langchain_chroma chromadb" % err)
```

Then in each consumer:
```python
from libtbx.langchain.rag._chroma_resilience import (
    ensure_chroma, chroma_unavailable_error)


def use_chroma(embeddings):
    _, Chroma = ensure_chroma()
    if Chroma is None:
        raise chroma_unavailable_error()
    return Chroma(embedding_function=embeddings)
```

**Critical rules:**

- **Catch `Exception` broadly, not just `ImportError`.**  Version
  conflicts in transitive deps surface as other exception types.
- **Don't make module-load conditional on the dep being usable.**
  Always import the module successfully.  Failure happens at
  use-time, not import-time.
- **Provide an informative error at use-time** with install hint.
- **Extract to a shared helper** when 2+ files in the same
  package need the same probe.  Single source of truth for the
  cache and error message.
- **For type annotations referring to the optional class**, use
  `Optional[Any]` from `typing` rather than the optional class
  itself.  Both `from __future__ import annotations` (PEP 563)
  and string literals like `-> "Chroma"` defer evaluation but
  create a hidden `NameError` trap if anything calls
  `typing.get_type_hints()` on the module while the dep is
  unavailable.  `Optional[Any]` has no such trap and preserves
  the "may return None" semantic.  Document the specific class
  in the docstring instead.

**For tests** that genuinely need the optional dep:

```python
_RAG_SKIP_REASON = None
try:
    import langchain_chroma  # noqa: F401
    import chromadb  # noqa: F401
except Exception as e:
    _RAG_SKIP_REASON = "%s: %s" % (type(e).__name__, e)


def _skip_if_rag_missing(cls):
    if _RAG_SKIP_REASON:
        for attr_name in list(dir(cls)):
            if attr_name.startswith("test_"):
                setattr(cls, attr_name,
                        unittest.skip(_RAG_SKIP_REASON)(
                            getattr(cls, attr_name)))
    return cls


@_skip_if_rag_missing
class TestRAGFeatures(unittest.TestCase):
    def test_chroma_query(self):
        ...
```

The class-level decorator skips an entire group of tests when the
probe detects the dep is unusable.  Use `Exception` (not just
`ImportError`) in the probe.

**Architectural note — LLM-visible tools:**

If a future architecture exposes RAG functionality as a tool the
LLM can choose (e.g., a `QueryDocumentationTool` in the planner's
tool list), gate the tool registration on `is_chroma_available()`:

```python
from libtbx.langchain.rag._chroma_resilience import is_chroma_available

if is_chroma_available():
    register_tool(QueryDocumentationTool())
else:
    log("Optional dependency 'chromadb' unavailable. "
        "Documentation search tool disabled.")
```

This prevents the LLM from seeing a tool it cannot use.  v118 does
NOT expose RAG as an LLM tool (the planner chooses among PHENIX
programs configured in `knowledge/programs.yaml`, none of which
depend on chromadb), so no tool-registration gating is currently
needed.

**Why "catch Exception" — concrete example:**

The chromadb failure that prompted this pattern manifests as:
```
TypeError: Descriptors cannot be created directly.
If this call came from a _pb2.py file, your generated code is
out of date and must be regenerated with protoc >= 3.19.0.
```

This is raised from inside `opentelemetry`'s `_pb2.py` during
chromadb's import.  `except ImportError` would not catch it.  A
defensive probe must use `except Exception`.

## 4. Adding a New Program

This guide documents the process for adding a new PHENIX program to the AI Agent. YAML definitions get the program running in the pipeline with no code changes to the core pipeline. In practice, testing will reveal edge cases that require Python guards — but the process always starts with YAML.

### Quick Overview

| File | What to Add | Required? |
|------|-------------|-----------|
| `knowledge/programs.yaml` | Complete program definition with metrics | ✅ Yes |
| `knowledge/workflows.yaml` | Add to appropriate workflow phase(s) | ✅ Yes |
| `knowledge/file_categories.yaml` | New output file types (if any) | If new file types |
| `knowledge/metrics.yaml` | Custom summary display (quality table, step metrics) | Usually not needed |
| `agent/directive_extractor.py` | Tutorial/procedure patterns | Optional |
| `phenix_ai/log_parsers.py` | Complex metric extraction (tables, multi-line) | Only if YAML can't handle it |

**To get the program running**, you only need `programs.yaml` and `workflows.yaml`. This integrates the program into BUILD, VALIDATE, the sanitizer, and the workflow engine automatically.

The system automatically handles:
- ✅ Metric extraction from logs (via `log_parsing` in programs.yaml)
- ✅ `<program>_done` tracking flags (via `done_tracking` in programs.yaml)
- ✅ Session summary display (default formatting works for most metrics)
- ✅ Workflow context building (automatic from programs.yaml)
- ✅ Prerequisite auto-execution (e.g., mtriage for resolution)
- ✅ skip_programs → workflow phase advancement (done flags auto-set)
- ✅ LLM hallucination sanitization (resolution stripped if unverified)
- ✅ .eff file generation with short-form PHIL resolution

---

### Minimal Example (2 files)

For a simple analysis program that runs once and extracts one metric:

**programs.yaml:**
```yaml
phenix.new_analysis:
  description: "Analyze something in the data"
  category: analysis
  experiment_types: [xray, cryoem]
  done_tracking:
    flag: "new_analysis_done"
    strategy: "run_once"
  
  inputs:
    required:
      data:
        extensions: [.mtz, .mrc]
        flag: ""
  
  command: "phenix.new_analysis {data}"
  
  log_parsing:
    quality_score:
      pattern: 'Quality Score[:\s]+([0-9.]+)'
      type: float
```

**workflows.yaml:**
```yaml
xray:
  phases:
    analyze:
      programs:
        - program: phenix.new_analysis
          conditions:
            - not_done: new_analysis
```

That's it! The agent will now:
- Run `phenix.new_analysis` during the analyze phase
- Extract `quality_score` from logs
- Track that it's been run (won't run again)
- Display the metric in session summaries

---

### Step 1: Define the Program in `programs.yaml`

Location: `knowledge/programs.yaml`

This is the **single source of truth** for program configuration.

```yaml
phenix.new_program:
  description: "Brief description of what the program does"
  category: analysis|phasing|building|refinement|validation
  experiment_types: [xray, cryoem]  # or just one
  
  # Done tracking: controls workflow done flags and strategy behavior
  done_tracking:
    flag: "new_program_done"   # Flag name set in workflow context
    strategy: "run_once"       # "set_flag" (default) | "run_once" | "count"
    # count_field: "new_count" # Required when strategy: "count"
    history_detection:         # Auto-set done flag from history
      markers: ["new_program"] # Strings to match (substring, OR logic)
      # exclude_markers: ["other_thing"]  # Reject match if present (checked FIRST)
      # alt_markers: ["alt_name"]  # Alt match (with alt_requires, AND logic)
      # alt_requires: ["special_flag"]

  # GUI app_id: maps to wxGUI2 window for project History (optional)
  gui_app_id: "NewProgram"           # From wxGUI2/Programs/__init__.params
  # gui_app_id_cryoem: "NewProgramCryoEM"  # If cryo-EM uses a different window

  # Stop directive patterns: regex for user stop conditions (optional)
  # stop_directive_patterns:
  #   - 'new_program'

  inputs:
    required:
      input_name:
        extensions: [.ext1, .ext2]
        flag: "input_flag="  # or "" for positional
    optional:
      optional_input:
        extensions: [.ext]
        flag: "optional_flag="
        # Set auto_fill: false to prevent the command builder from
        # automatically filling this slot with a matching file.
        # Use when a slot should only be filled if the LLM or user
        # explicitly assigns a file to it (e.g., a second map for
        # map-map correlation — without this, both map slots would
        # get the same file).
        auto_fill: false  # Default is true

  # Control which file categories are preferred for each input
  # IMPORTANT: Without this, the agent may select wrong files!
  input_priorities:
    input_name:
      categories: [preferred_category, fallback_category]
      prefer_subcategories: [refined, with_ligand]  # Within category, prefer these
      exclude_categories: [excluded_category]
      # Filename pattern filters (word-boundary matching):
      exclude_patterns: [ligand, refine]  # Reject files with these words in name
      prefer_patterns: [with_ligand]      # Prefer files matching these words
      # Content-based guards (automatic, no config needed):
      # - Model slots: HETATM-only PDB files rejected (small molecules)
      # - Ligand slot: PDB files with ATOM records rejected (protein models)
      require_best_files_only: true  # Only use best_files, not file scan fallback

  outputs:
    files:
      - pattern: "*.output_ext"
        type: output_type  # Used for file categorization
    metrics:
      - metric_name  # Metrics this program produces

  command: "phenix.new_program {input_name}"

  # Optional parameters that can be set via strategy
  strategy_flags:
    parameter_name:
      flag: "parameter={value}"
      type: float|int|string|bool
      hint: "Description of parameter"

  # IMPORTANT: Define metric extraction patterns here
  # These are used automatically by log_parsers.py and session.py
  log_parsing:
    metric_name:
      pattern: 'Regex pattern with (capture group)'
      type: float|int|string|bool
      display_name: "Human Readable Name"    # For summary display
      summary_format: "{value:.2f}"          # Format for display
      # Optional: Handle "no result" cases
      no_match_pattern: 'No.*found'
      no_match_value: "None"

  hints:
    - "Hint for when/how to use this program"

  # Keywords that trigger this program from user advice
  user_advice_keywords:
    - "keyword1"
    - "keyword2"
```

### Key Fields for Automatic Handling

| Field | Effect |
|-------|--------|
| `done_tracking.flag` | Names the workflow done flag (e.g., `xtriage_done`) |
| `done_tracking.strategy` | `"set_flag"` (default), `"run_once"` (filter after first run), or `"count"` (increment counter) |
| `done_tracking.count_field` | Counter name for `strategy: "count"` (e.g., `"refine_count"`) — must be in ALLOWED_COUNT_FIELDS |
| `done_tracking.history_detection` | YAML-driven done-flag setting (markers, exclude_markers, alt_markers, success_flag) |
| `gui_app_id` | wxGUI2 app_id for project History (fallback when PHIL unavailable) |
| `gui_app_id_cryoem` | Cryo-EM variant app_id if the program uses a different GUI window |
| `stop_directive_patterns` | Regex patterns for matching user stop-condition directives |
| `log_parsing:` | Patterns used by log_parsers.py AND session.py for metric extraction |
| `log_parsing.*.display_name` | Used in session summary display |
| `log_parsing.*.summary_format` | Format string for summary display |
| `exclude_patterns` | Reject files with matching words (word-boundary matching: `ligand` matches `atp_ligand.pdb` but NOT `noligand.pdb`) |
| `prefer_patterns` | Prefer files matching these words (same word-boundary semantics) |
| `require_best_files_only` | Only fill this slot from `best_files`, skip extension-scan fallback (e.g., ligandfit's `map_coeffs_mtz`) |
| `auto_fill` | Set to `false` to prevent auto-fill; slot only filled if LLM/user assigns a file |

---

### Step 2: Add to Workflow in `workflows.yaml`

Location: `knowledge/workflows.yaml`

Find the appropriate phase and add your program:

```yaml
xray:  # or cryoem
  phases:
    phase_name:
      programs:
        # Simple entry (always available)
        - phenix.new_program
        
        # Or with conditions
        - program: phenix.new_program
          conditions:
            - has: required_file_type  # e.g., has: sequence
            - not_done: new_program    # Don't run if already done
          hint: "When to use this program"
```

**Common conditions:**
- `has: <file_type>` - Requires file (sequence, model, full_map, half_map, ligand_file, anomalous, etc.)
- `has_any: [type1, type2]` - Requires at least one of these files
- `not_done: <program>` - Program hasn't been run yet (flag name from `done_tracking.flag`)
- `not_has: <file_type>` - File type is NOT present
- `r_free: "< 0.35"` - R-free below threshold (also supports `"> autobuild_threshold"`, `"< target_r_free"`)
- `refine_count: "> 0"` - At least one successful refinement completed (also `rsr_count`)

**Available done flags:**

All done flags are defined in programs.yaml under `done_tracking.flag`. Common flags:

| Flag | Set When | Used By |
|------|----------|---------|
| `predict_done` | predict_and_build runs (any mode) | cryo-EM workflow |
| `predict_full_done` | predict_and_build completes full workflow | X-ray workflow |
| `process_predicted_done` | process_predicted_model succeeds | cryo-EM workflow |
| `phaser_done` | phaser completes MR | X-ray workflow |
| `phaser_count` | Count of successful phaser runs | Internal tracking |
| `dock_done` | dock_in_map completes | cryo-EM workflow |
| `autobuild_done` | autobuild succeeds | X-ray workflow |
| `autobuild_denmod_done` | autobuild maps_only completes | X-ray workflow (before ligandfit) |
| `autosol_done` | autosol succeeds | X-ray workflow |
| `ligandfit_done` | ligandfit completes | Both workflows |
| `pdbtools_done` | pdbtools completes | X-ray workflow (combine_ligand) |
| `refine_done` | refine completes (any) | Internal tracking |
| `refine_count` | Count of successful refinements | Workflow conditions |
| `rsr_done` | real_space_refine completes | Internal tracking |
| `rsr_count` | Count of successful RSR runs | Workflow conditions |
| `resolve_cryo_em_done` | resolve_cryo_em completes | cryo-EM workflow |
| `map_sharpening_done` | map_sharpening completes | cryo-EM workflow |
| `map_to_model_done` | map_to_model completes | cryo-EM workflow |
| `validation_done` | molprobity/validation completes | Internal tracking |

**Note:** `refine_count` and `rsr_count` only count SUCCESSFUL runs (not failed ones).

---

### Step 2.5: Add Done Tracking (if needed)

Location: `knowledge/programs.yaml` (same file as Step 1)

**All done flags are defined in the `done_tracking` block** in programs.yaml.
No Python file edits are needed — `get_program_done_flag_map()` reads directly
from YAML.

```yaml
phenix.new_program:
  # ... other fields ...
  done_tracking:
    flag: "new_program_done"    # Flag name in workflow context
    strategy: "run_once"        # "set_flag" | "run_once" | "count"
    history_detection:          # YAML-driven done-flag detection
      markers: ["new_program"]  # Match any of these in history (OR logic)
      # exclude_markers: [...]  # Reject match if present (checked FIRST)
      # alt_markers: ["other_name"]  # Alternative markers (with alt_requires)
      # alt_requires: ["some_flag"]  # ALL must match alongside alt_markers
      # success_flag: "new_program_success"  # Extra flag set on success
```

**Strategies:**
- `"set_flag"` (default): Set done flag on success. Most programs use this.
- `"run_once"`: Set done flag + filter program from valid list after first run.
  Use for analysis programs that should never re-run (xtriage, mtriage).
- `"count"`: Set done flag + increment a counter. Requires `count_field`.
  Use for programs that can run multiple times (refine, rsr).

**`history_detection`** replaces if/elif blocks in `_analyze_history()`.
Use it when your program only needs a boolean done flag set on success.
Do NOT use it if your program needs cascading flags (like
predict_and_build → refine_done), or exclusion logic (like refine excluding
real_space).

**`stop_directive_patterns`** (optional, top-level): If users may reference your
program in stop conditions (e.g., "stop after running new_program"), add regex
patterns:

```yaml
phenix.new_program:
  stop_directive_patterns:
    - 'new_program'
    - 'new\\s+program'
```

Patterns are sorted by length at load time (longest first) so more specific
patterns match before shorter ones.

When the user's advice matches one of these patterns AND contains an
explicit stop-after phrasing recognised by `_is_stop_after_requested()`
(v116.x — `"stop after X"`, `"X and stop"`, `"only run X"`, etc.), the
directive extractor sets both `stop_conditions.after_program=phenix.new_program`
AND `stop_conditions.stop_after_requested=True`.  Without the flag,
`after_program` is treated as a min-run hint only — see section 3a
"`after_program` semantics".  If you want your program to support
user-explicit stops, the helper's standard positive patterns cover the
common phrasings without per-program work.  Adding patterns specific to
your program (e.g., `"run_only_new_program"`) is only needed if users
might use a phrasing the helper doesn't recognise.

**When to add `done_tracking`:**
- Your program gates a workflow phase (e.g., xtriage gates "analyze")
- Users might skip it via `skip_programs` (the flag must be set for phase advancement)
- You want `strategy: "run_once"` behavior (program filtered from valid list after completion)
- You need to track how many times a program ran (`strategy: "count"`)

**When to skip:** Validation programs like `phenix.map_correlations` that don't
gate any phase and can run freely don't need `done_tracking`.

**Also ensure** the done flag is set in `build_context()` in `agent/workflow_engine.py`
(the flags are populated from `history_info`):

```python
context = {
    # ... existing flags ...
    "new_program_done": history_info.get("new_program_done", False),
}
```

And in `_analyze_history()` in `agent/workflow_state.py`, detect when the program
has been run and set the flag in `history_info`.

---

### Step 3: Add Output File Types (if needed)

Location: `knowledge/file_categories.yaml`

If your program produces a **new type of output file**, define it here:

```yaml
new_file_type:
  description: "Description of what this file contains"
  extensions: [.ext]
  patterns:
    - "*.ext"
    - "specific_name.ext"
  notes: "Which programs produce/consume this file type"
```

**Skip this step if** your program only produces standard file types (`.pdb`, `.mtz`, `.mrc`, etc.) that are already defined.

**Example** - map_symmetry produces `.ncs_spec` files:
```yaml
ncs_spec:
  description: "NCS symmetry operators file"
  extensions: [.ncs_spec]
  patterns:
    - "*.ncs_spec"
  notes: "From map_symmetry, used by resolve_cryo_em, predict_and_build"
```

---

### Step 4: Add Summary Display (usually not needed)

Location: `knowledge/metrics.yaml`

**Most programs don't need this step.** The default formatting works for standard metrics. Only add custom display if you need:
- Special formatting in the Final Quality table
- Custom step metric display in the workflow steps table

Location: `knowledge/metrics.yaml`

If your program produces metrics that should appear in the session summary, add them to the `summary_display` section:

### For the Final Quality Table:

```yaml
summary_display:
  quality_table:
    # Add your metric
    - metrics: [new_metric]
      label: "New Metric"
      format: "{new_metric:.2f}"
      experiment_types: [xray]  # Optional filter
      assessment: true          # Show Good/Acceptable/Needs Improvement
```

### For Per-Step Display:

```yaml
summary_display:
  step_metrics:
    phenix.new_program:
      format: "Result: {new_metric:.2f}"
      metrics: [new_metric]
      fallback_format: "Completed"  # When metrics unavailable
```

---

### Step 5: Add Tutorial Detection (Optional)

Location: `agent/directive_extractor.py`

If users might ask to run this specific program as a tutorial or procedure:

```python
tutorial_patterns = [
    # ... existing patterns ...
    (r'(?:run|analyze).*new_program', 'phenix.new_program'),
    (r'keyword.*pattern', 'phenix.new_program'),
]
```

---

### Step 6: Add Tests

Verify your program works correctly:

```bash
# Validate program configuration
python agent/program_validator.py phenix.new_program

# Run all tests — including the integration suite
python tests/run_all_tests.py --quick
```

The validator checks:
- ✅ Program defined in programs.yaml
- ✅ Added to workflows.yaml
- ✅ Log parsing patterns defined
- ✅ Summary display configured (if applicable)

### Run the full suite, not just your new tests

`workflow_state.py` is the integration boundary. When adding new behavior there
(zombie checks, file validation, diagnostic surfacing), it should be treated as a
refactor of the public API — requiring a full suite run, not just the feature's own
tests. Any PR touching `detect_workflow_state()` must show passing results for
`tst_workflow_state.py` and `tst_integration.py` before merging.

The reason this matters: the pre-existing tests in those files use hypothetical
filename strings (e.g. `"data.mtz"`, `"refined_model.pdb"`) without creating real
files on disk, and use generic output names rather than phenix naming conventions.
New features that inspect files or match output patterns by name can break these
tests invisibly — passing their own unit tests while silently breaking all callers.

---

### Step 6.5: Zombie State Checker (if your program produces a model or map)

If your program produces a **model file** (`.pdb`, `.cif`) or **map file**
(`.mrc`, `.ccp4`) and uses `strategy: "set_flag"` or `strategy: "count"`, you
should add an entry to `_ZOMBIE_CHECK_TABLE` in `agent/workflow_state.py`.

### What is a zombie state?

A zombie state occurs when the agent crashes mid-cycle or output files are
deleted after the fact. History records `done_flag=True`, but the output file
is no longer on disk. The phase detector sees `done=True` and skips the
program — the workflow becomes stuck.

The zombie checker (`_clear_zombie_done_flags`) runs every PERCEIVE cycle. It
checks each done flag against an expected output file pattern. If the output is
missing, it clears the done flag so the program can re-run.

### When to add a zombie check entry

Add an entry when **all** of the following are true:
- The program produces an output file that later steps depend on
- The done flag is used to gate a workflow phase (e.g., `has_placed_model`,
  `has_refined_model`, or a count like `refine_count`)
- A crash or cleanup that removes the output would leave the workflow stuck

You **do not** need an entry for:
- Analysis programs (xtriage, mtriage) — their outputs are metrics, not files
- Programs with `strategy: "run_once"` for analysis — they don't gate file-dependent phases
- Programs where the done flag is advisory only (e.g., `validation_done`)

### How to add the entry

In `agent/workflow_state.py`, find `_ZOMBIE_CHECK_TABLE` and add a tuple:

```python
_ZOMBIE_CHECK_TABLE = [
    # ... existing entries ...

    # Your new program
    ("new_program_done",
     _re.compile(r"<output_pattern>.*\.pdb$", _re.IGNORECASE),
     "has_placed_model",   # context flag to also clear (or None)
     True),                # accept_any_pdb (see below)
]
```

**Fields:**
| Field | Type | Description |
|-------|------|-------------|
| `done_flag` | str | Key in `history_info` to check (e.g., `"dock_done"`) |
| `output_pattern` | regex | Matched against basenames of available files |
| `file_flag_to_clear` | str or None | Context flag also cleared when zombie detected (e.g., `"has_placed_model"`) |
| `accept_any_pdb` | bool | If True, any `.pdb` on disk suppresses zombie detection as a fallback |

### Pattern alignment with file_categories.yaml

The `output_pattern` regex should be consistent with the patterns in
`knowledge/file_categories.yaml` for the same output type. If the file category
uses `*dock*map*` as a pattern, the zombie regex should also match `dock.*map.*\.pdb`.

**Always check:** does the zombie regex match the filenames phenix actually produces,
as well as what the file category accepts? They must agree.

```bash
# Quick check: test your pattern against real phenix output names
python3 -c "
import re
pat = re.compile(r'your_pattern_here', re.IGNORECASE)
names = ['expected_output_001.pdb', 'alternate_name.pdb']
for n in names: print(n, bool(pat.search(n)))
"
```

### The `accept_any_pdb` flag

| Value | Use when |
|-------|----------|
| `True` | Program produces a **model** file (`.pdb`/`.cif`). Any PDB on disk is plausible evidence of success, preventing false zombie clears when files are renamed or tests use generic names. |
| `False` | Program produces a **specifically-named** file that downstream steps reference by name (e.g., `denmod_map.ccp4` from `resolve_cryo_em`). Pattern match must be exact. |

**Trade-off:** `accept_any_pdb=True` weakens zombie detection in the edge case
where an older model file from a *different* program is on disk when the newer
program's output was deleted. This is an acceptable heuristic — zombie detection
is best-effort, not a guarantee. Document this choice in a comment.

### Example: dock_in_map

```python
# dock_in_map produces *_docked.pdb, *dock*map*.pdb, placed_model*.pdb, or *_placed*.pdb
# Pattern mirrors the "docked" file category in file_categories.yaml.
# accept_any_pdb=True: any PDB is acceptable fallback (model renamed or test proxy)
("dock_done",
 _re.compile(r"docked.*\.pdb$|.*_docked.*\.pdb$|dock.*map.*\.pdb$|placed_model.*\.pdb$|.*_placed.*\.pdb$",
             _re.IGNORECASE),
 "has_placed_model", True),
```

### Count-based programs

For programs with `strategy: "count"` (refine, rsr), the zombie checker
automatically decrements the count when it clears the done flag:

```python
("refine_done",
 _re.compile(r"refine_\d+\.pdb$", _re.IGNORECASE),
 None,   # No separate context flag — refine_count is decremented automatically
 True),
```

The decrement logic is built into `_clear_zombie_done_flags`. You don't need to
add code for it — just set `file_flag_to_clear=None` and the function handles
`refine_count` and `rsr_count` automatically.

---

### Example: Adding phenix.map_symmetry

Here's the complete configuration for map_symmetry:

### programs.yaml

```yaml
phenix.map_symmetry:
  description: "Find point-group symmetry in cryo-EM map"
  category: analysis
  experiment_types: [cryoem]
  done_tracking:
    flag: "map_symmetry_done"
    strategy: "run_once"
    history_detection:
      markers: ["map_symmetry"]

  inputs:
    required:
      map:
        extensions: [.mrc, .ccp4, .map]
        flag: ""

  input_priorities:
    map:
      categories: [full_map, optimized_full_map, map]
      exclude_categories: [half_map]

  outputs:
    files:
      - pattern: "*.ncs_spec"
        type: ncs_spec
    metrics:
      - symmetry_type
      - ncs_copies
      - ncs_cc

  command: "phenix.map_symmetry {map}"

  # Metric extraction - used by log_parsers.py AND session.py
  log_parsing:
    symmetry_type:
      pattern: 'NCS\s+type[:\s]+([A-Z]\d+(?:\s*\([a-z]\))?)'
      type: string
      display_name: "Symmetry Type"
      summary_format: "{value}"
      no_match_pattern: 'No suitable symmetry found'
      no_match_value: "None"
    ncs_copies:
      pattern: 'Number of NCS copies[:\s]+(\d+)'
      type: int
      display_name: "NCS Copies"
      summary_format: "{value}"
    ncs_cc:
      pattern: 'NCS correlation[:\s]+([0-9.]+)'
      type: float
      display_name: "NCS CC"
      summary_format: "{value:.2f}"

  hints:
    - "Run to detect point-group symmetry in cryo-EM maps"
```

### workflows.yaml

```yaml
cryoem:
  phases:
    analyze:
      programs:
        - phenix.mtriage
        - program: phenix.map_symmetry
          conditions:
            - not_done: map_symmetry
          hint: "Run to detect map symmetry"
```

### metrics.yaml (summary_display section)

```yaml
summary_display:
  quality_table:
    - metrics: [symmetry_type]
      label: "Symmetry"
      format: "{symmetry_type}"
      optional_format: " ({ncs_copies} copies)"
      optional_metrics: [ncs_copies]
      detail_format: "CC: {ncs_cc:.2f}"
      detail_metrics: [ncs_cc]
      experiment_types: [cryoem]

  step_metrics:
    phenix.map_symmetry:
      format: "Symmetry: {symmetry_type}"
      metrics: [symmetry_type]
      optional_format: " ({ncs_copies} copies)"
      optional_metrics: [ncs_copies]
      fallback_format: "Analyzed"
```

This gets the program running — 3 YAML files, no Python changes.

### Reality check: what testing reveals

The YAML definitions above are the starting point, not the finish
line. When you run the new program through the tutorial suite, you
will likely encounter edge cases that need Python guards. As of
v115.05, 17 of 23 registered programs have at least some
program-specific Python code (~60 individual guard blocks). Only
6 programs are pure YAML.

Common edge cases that surface during testing:

| Category | Where | Example |
|----------|-------|---------|
| File categorization guards | `workflow_state.py` | Ligand vs protein PDB, half-map vs full-map |
| Content-based workflow checks | `workflow_state.py` | Se heavy-atom files misidentified as ligands |
| Error recovery patterns | `recoverable_errors.yaml` | Ambiguous data labels for this program |
| Exclude patterns | `command_postprocessor.py` | Program strips certain params from other programs |
| Special-case command building | `graph_nodes.py` BUILD | `rebuild_in_place=False` for autobuild |
| Done-tracking edge cases | `graph_nodes.py` | Program that should run more than once |

**Process**: Add the YAML → run tutorials → fix the edge cases that
appear → re-run. Expect 2-3 iterations. The core pipeline
(`PERCEIVE → THINK → PLAN → BUILD → VALIDATE → OUTPUT`) never
changes — all guards are in the surrounding code.

---

### Example: Adding phenix.map_correlations (Multi-Mode Inputs)

This example demonstrates a program with multiple input modes (model+map, map+map,
mtz+mtz) where `auto_fill: false` prevents the command builder from duplicating files.

### programs.yaml

```yaml
phenix.map_correlations:
  description: "Compute map-model or map-map correlation coefficients"
  category: validation
  experiment_types: [xray, cryoem]

  inputs:
    # All inputs are optional — different modes use different combinations
    optional:
      model:
        extensions: [.pdb, .cif]
        description: "Model file"
        flag: "input_files.model="
      full_map:
        extensions: [.ccp4, .mrc, .map]
        description: "Map 1"
        flag: "input_files.map_in_1="
      map2:
        extensions: [.ccp4, .mrc, .map]
        description: "Map 2 for map-map CC"
        flag: "input_files.map_in_2="
        auto_fill: false  # CRITICAL: prevents same map filling both slots
      data_mtz:
        extensions: [.mtz]
        description: "Map coefficients 1"
        flag: "input_files.map_coeffs_1="
      map_coeffs_2:
        extensions: [.mtz]
        description: "Map coefficients 2"
        flag: "input_files.map_coeffs_2="
        auto_fill: false  # Same reason as map2

  command: "phenix.map_correlations"  # Files appended via flags

  log_parsing:
    cc_mask:
      pattern: "CC_mask\\s*[=:]\\s*([0-9.]+)"
      type: float
      display_name: "CC_mask"
      summary_format: "{value:.3f}"
    cc_volume:
      pattern: "CC_volume\\s*[=:]\\s*([0-9.]+)"
      type: float
      display_name: "CC_volume"
      summary_format: "{value:.3f}"
    # ... additional metrics
```

**Key design decisions:**

1. **All inputs optional:** Different modes need different combinations. The command
   builder appends only the slots that have files assigned.

2. **`auto_fill: false` on secondary slots:** Without this, if the user provides one
   `.ccp4` map, the builder would auto-fill both `full_map` AND `map2` with the same
   file, producing a meaningless self-correlation command.

3. **Flag-based append:** No `{placeholder}` in the command template. Files are
   appended as `flag + filepath` (e.g., `input_files.model=/path/to/model.pdb`).

### metrics.yaml additions

```yaml
# In summary_display.step_metrics:
phenix.map_correlations:
  format: "CC_mask: {cc_mask:.3f}"
  metrics: [cc_mask]
  fallback_format: "Computed correlations"

# In summary_display.quality_table:
- metrics: [cc_mask]
  label: "CC_mask"
  format: "{cc_mask:.3f}"
  detail_format: "CC_volume: {cc_volume:.3f}"
  detail_metrics: [cc_volume]
  assessment: true
```

### Important: Pattern matching for reformatted metrics

The `format_metrics_report()` function transforms metric keys via
`.replace("_", " ").title()`, so `cc_mask` becomes `Cc Mask` in the result text
stored in cycle records. Any hardcoded regex patterns in `session.py` must match
both forms:

```python
# Matches both "CC_mask  : 0.849" (raw log) and "Cc Mask: 0.849" (report)
r'CC[_ ]?mask\s*[=:]\s*([0-9.]+)'  # with re.IGNORECASE
```

---

### Key Design Patterns

### Prerequisites (Auto-Run Dependencies)

If your program requires a value that must come from another program (e.g.,
resolution from mtriage), declare a `prerequisite` on the invariant:

```yaml
phenix.new_program:
  invariants:
    resolution:
      required_key: resolution
      source: session
      prerequisite: phenix.mtriage   # ← Auto-runs mtriage if resolution missing
```

When the command builder detects the invariant can't be satisfied, it
automatically builds and runs the prerequisite program first. The original
program runs on the next cycle with the newly-available value.

**Important:** The prerequisite program must be able to produce the required
value. The build node checks `skip_programs` directives — if the user skipped
the prerequisite, a clear error is returned instead of silently failing.

### File Selection: exclude_categories for Map Inputs

When a program accepts maps, always exclude half-maps from the full_map slot.
Half-maps bubble up to the parent `map` category, so without exclusion, a
half-map can be selected as a full_map via the category fallback:

```yaml
input_priorities:
  full_map:
    categories: [full_map, optimized_full_map, map]
    exclude_categories: [half_map]   # ← Prevents half-map selection
```

### Half-Map Handling Flags

Three flags control how the dedup logic handles programs with both `full_map`
and `half_map` slots:

| Flag | Effect | Programs |
|------|--------|----------|
| `prefers_half_maps: true` | Drop full_map, keep half-maps | resolve_cryo_em, mtriage, map_sharpening |
| `keep_half_maps_with_full_map: true` | Keep both | map_to_model |
| `requires_full_map: true` | Strip any half-maps | dock_in_map, real_space_refine, validation_cryoem |
| (none) | Drop half-maps, keep full_map | Default behavior |

Use `prefers_half_maps` when the program works best with half-maps (e.g., FSC
resolution from half-map pairs) or when full maps may have different grid
dimensions from post-processing. Use `keep_half_maps_with_full_map` only when
the program genuinely needs both simultaneously.

### Supplement Logic for `multiple:true` Slots

LLMs often assign only one file to a `multiple:true` slot (e.g.,
`half_map=file2.ccp4`). Auto-fill skips the slot because it's already
"filled". A supplement loop in `CommandBuilder._select_files()` checks each
`multiple:true` slot after auto-fill and backfills missing files from the
category using `_find_file_for_slot()`. This ensures programs like
resolve_cryo_em always receive both half-maps.

The supplement respects `auto_fill: false` — if the LLM never assigned a file
to the slot, it remains empty (the supplement only completes partially-filled
slots, not empty ones).

### Intermediate File Exclusions

If your program writes intermediate files to subdirectories, ensure they don't
get tracked as outputs. Add exclusions in three places:

1. **`phenix_ai/log_parsers.py`** — `extract_output_files()`: Add path/basename
   patterns to the exclusion checks
2. **`phenix_ai/utilities.py`** — `scan_log_for_files()`: Same exclusion patterns
3. **`programs/ai_agent.py`** — `_track_output_files()`: Add to
   `intermediate_patterns` list and/or `intermediate_basename_prefixes`

Example for predict_and_build's internal docking intermediates:
```python
intermediate_patterns = [
    '/local_dock_and_rebuild',   # predict_and_build intermediates
    '/local_rebuild',
]
intermediate_basename_prefixes = [
    'working_model',              # Never a valid output
]
```

Also update `knowledge/file_categories.yaml` to exclude intermediates from
semantic categories (e.g., add `working_model*` to the `docked` category's
`excludes` list).

### LLM Strategy Sanitization

If your program has parameters where the LLM commonly hallucinates values,
add sanitization in `_build_strategy()` in `agent/command_builder.py`.

Current sanitizations:
- **resolution**: Stripped from LLM strategy when no verified source exists
  (prevents hallucinated resolution values from bypassing the mtriage prerequisite)
- **atom_type** (autosol): Multi-atom values like "Se, S" are split — first atom
  stays in `atom_type`, extras move to `additional_atom_types`

### .eff File Generation

Programs that use old-style master_phil (like phenix.refine) have their .eff
files generated via `command_line_argument_interpreter()` which resolves
short-form PHIL paths to full scope paths. This is handled in
`_generate_eff_file()` in `programs/ai_agent.py` (Attempt 3).

---

### Architecture Notes

### How It Works

1. **Metric Extraction** (`knowledge/metric_patterns.py`)
   - Reads `log_parsing` from programs.yaml
   - Compiles patterns once, caches them
   - Used by both `log_parsers.py` and `session.py`

2. **Program Registration** (`knowledge/program_registration.py`)
   - Reads `done_tracking` blocks from programs.yaml
   - Provides `get_program_done_flag_map()` (all programs → done flags)
   - Provides `get_trackable_programs()` (strategy: "run_once" programs only)

3. **Summary Display** (`knowledge/summary_display.py`)
   - Reads `summary_display` from metrics.yaml
   - Formats quality table and step metrics
   - Used by `session.py` for summary generation

### Adding Complex Parsing

For programs that need complex parsing (tables, multi-line, etc.), you can still add hardcoded extractors in `log_parsers.py`. The YAML patterns are tried first, then hardcoded extractors fill in any missing metrics.

**Future alternative: `results_as_json()`.** If the program you're
adding is built on `ProgramTemplate` and supports
`results_as_json()`, it may be possible to read structured JSON
results directly instead of writing `log_parsing` patterns. This
is not yet implemented in the agent's PERCEIVE node, but is planned
(see ARCHITECTURE.md "Potential improvements"). For now, add
`log_parsing` patterns as described above — they'll serve as the
fallback path once JSON results are supported.

```python
# In log_parsers.py
def _extract_new_program_metrics(log_text):
    """Handle complex parsing that YAML can't express."""
    metrics = {}
    # Complex table parsing, multi-line context, etc.
    return metrics
```

Add the call in `extract_all_metrics()`:

```python
elif "new_program" in prog_lower:
    hardcoded_metrics = _extract_new_program_metrics(log_text)
    for k, v in hardcoded_metrics.items():
        if k not in metrics:  # Don't overwrite YAML results
            metrics[k] = v
```

---

### Validation

Use the program validator to check configuration:

```bash
# Validate a specific program
python agent/program_validator.py phenix.map_symmetry

# Validate all programs
python agent/program_validator.py --all

# List all defined programs
python agent/program_validator.py --list
```

---

### Debugging Tips

1. **Program not appearing in valid_programs?**
   - Check conditions in `workflows.yaml`
   - Verify file categorization is correct
   - Check if `_done` flag is blocking it
   - If `after_program` is set and target has run AND
     `stop_after_requested=True` (v116.x), valid_programs is
     wiped to `[STOP]`.  Check the directive extraction logs for
     `"DIRECTIVES: stop_after_requested"`.  If False / absent,
     the workflow advances naturally — see section 3a
     "`after_program` semantics".

2. **Metrics not extracted?**
   - Test regex pattern against actual log output
   - Check pattern in `programs.yaml` `log_parsing` section
   - Run: `python -c "from knowledge.metric_patterns import extract_metrics_for_program; print(extract_metrics_for_program(open('log.txt').read(), 'phenix.new_program'))"`

3. **Metrics not in summary?**
   - Check `summary_display` section in `metrics.yaml`
   - Verify `experiment_types` filter matches your data
   - Check required metrics are being extracted

4. **Command built incorrectly?**
   - Check `programs.yaml` input definitions
   - Verify `input_priorities` categories match file categorization

5. **Wrong file selected for input?**
   - Add explicit `input_priorities` with `exclude_categories`
   - Common mistake: protein input picks ligand file → add `exclude_categories: [ligand, ligand_fit_output]`
   - Common mistake: protein input picks search model → add `exclude_categories: [search_model, predicted]`

---

### Common Mistakes to Avoid

| Mistake | Symptom | Fix |
|---------|---------|-----|
| Missing `input_priorities` | Wrong file selected | Add explicit priorities with excludes |
| No `exclude_categories` for protein/model input | Picks ligand or predicted model | Add `exclude_categories: [ligand, search_model]` |
| Forgetting `done_tracking` with `strategy: "run_once"` | Program runs multiple times | Add `done_tracking: {flag: ..., strategy: "run_once"}` |
| Wrong `experiment_types` | Program offered for wrong data | Set `[xray]`, `[cryoem]`, or both |
| Missing workflow entry | Program never offered | Add to appropriate phase in `workflows.yaml` |

### Example: Proper input_priorities for a refinement program

```yaml
phenix.my_refine:
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
      data_mtz:
        extensions: [.mtz]
        flag: ""
  
  # CRITICAL: Prevent selecting wrong files
  input_priorities:
    model:
      categories: [model]
      prefer_subcategories: [refined, with_ligand, phaser_output]
      exclude_categories: [ligand, search_model, ligand_fit_output, predicted]
    data_mtz:
      categories: [data_mtz]
```

---

## 5. Testing

> **Prompt templates:** When writing a new graph node that calls an LLM, the
> prompt must live in `agent/prompts/` — see
> [CCTBX_LLM_PROGRAMMING_GUIDELINES.md](CCTBX_LLM_PROGRAMMING_GUIDELINES.md) §3.

This document describes the testing infrastructure for the PHENIX AI Agent.

### Overview

The test suite uses a **cctbx-style testing approach** with plain functions and fail-fast behavior, meaning the first assertion failure stops execution with a full traceback. This matches the testing conventions used throughout the PHENIX/cctbx ecosystem.

### Quick Start

```bash
# Run all tests
cd improved_agent_v2v1
python tests/run_all_tests.py

# Run quick mode (standalone tests only, no PHENIX required)
python tests/run_all_tests.py --quick

# Run a single test file
python tests/tst_directive_extractor.py

# Run tests matching a pattern
python tests/run_all_tests.py --pattern "directive"
```

### Test Architecture

### Fail-Fast Behavior

Unlike Python's standard `unittest` which collects all failures, our tests crash immediately on the first failure:

```python
def test_something():
    assert_equal(result, expected)  # If this fails, execution stops here
    assert_true(condition)          # This line never runs if above fails
```

This provides immediate feedback with full tracebacks, making debugging faster.

### Assert Helper Functions

All tests use helper functions from `tests/tst_utils.py` instead of `unittest.TestCase` assertions:

| Function | Purpose |
|----------|---------|
| `assert_equal(a, b)` | Check `a == b` |
| `assert_not_equal(a, b)` | Check `a != b` |
| `assert_true(x)` | Check `x` is truthy |
| `assert_false(x)` | Check `x` is falsy |
| `assert_none(x)` | Check `x is None` |
| `assert_not_none(x)` | Check `x is not None` |
| `assert_in(item, container)` | Check `item in container` |
| `assert_not_in(item, container)` | Check `item not in container` |
| `assert_is(a, b)` | Check `a is b` |
| `assert_is_not(a, b)` | Check `a is not b` |
| `assert_is_instance(obj, cls)` | Check `isinstance(obj, cls)` |
| `assert_greater(a, b)` | Check `a > b` |
| `assert_greater_equal(a, b)` | Check `a >= b` |
| `assert_less(a, b)` | Check `a < b` |
| `assert_less_equal(a, b)` | Check `a <= b` |
| `assert_almost_equal(a, b, places)` | Check floats equal to N decimal places |
| `assert_raises(exc, func, *args)` | Check function raises exception |
| `assert_startswith(s, prefix)` | Check string starts with prefix |
| `assert_endswith(s, suffix)` | Check string ends with suffix |
| `assert_contains(text, substring)` | Check substring in text |

### Test File Structure

Each test file follows this structure:

```python
"""
Tests for module_name.

Run with: python tests/tst_module_name.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in,
    run_tests_with_fail_fast
)

from agent.module_name import function_to_test


# =============================================================================
# SECTION NAME
# =============================================================================

def test_basic_functionality():
    """Description of what this tests."""
    result = function_to_test("input")
    assert_equal(result, "expected")


def test_edge_case():
    """Test edge case handling."""
    result = function_to_test(None)
    assert_none(result)


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
```

### Test Discovery

The `run_tests_with_fail_fast()` function automatically discovers tests:

1. Finds all functions starting with `test_` in the calling module
2. Also finds `unittest.TestCase` classes (for backward compatibility)
3. Runs them in sorted order
4. Stops at first failure with full traceback

### Test Files

| File | Tests | Description |
|------|-------|-------------|
| `tst_directive_extractor.py` | 87 | Directive extraction and validation; crystal symmetry (unit_cell/space_group) extraction and normalisation |
| `tst_transport.py` | 66 | REST transport sanitization |
| `tst_workflow_state.py` | 66 | Workflow state detection, done flags, refine counting, MR-SAD |
| `tst_best_files_tracker.py` | 50 | File tracking, scoring, model stage detection |
| `tst_directive_validator.py` | 48 | Intent validation and modification |
| `tst_new_programs.py` | 44 | YAML config for new programs (polder, map_sharpening, autobuild_denmod) |
| `tst_error_analyzer.py` | 43 | Error recovery system |
| `tst_langchain_tools.py` | 34 | Legacy module tests (core, analysis, RAG, validation, prompts, memory) |
| `tst_sanity_checker.py` | 29 | Sanity check logic |
| `tst_yaml_config.py` | 29 | YAML configuration validation |
| `tst_file_categorization.py` | 26 | File type detection (including denmod_mtz) |
| `tst_api_schema.py` | 26 | API request/response validation |
| `tst_decision_flow.py` | 21 | Directive flow architecture |
| `tst_metrics_analyzer.py` | 19 | Metric extraction and trends |
| `tst_pattern_manager.py` | 17 | Pattern management |
| `tst_directives_integration.py` | 16 | End-to-end directive tests |
| `tst_command_builder.py` | 22 | Unified command generation, MR-SAD partpdb_file |
| `tst_rfree_generate_guard.py` | 13 | phenix.refine never regenerates an existing R-free test set (tri-state client-fact detection; server-undetermined→strip) |
| `tst_rfree_lock_reconciliation.py` | 10 | F1/F2 (v120.3–v120.4): a flagless or unverifiable data file selected over the locked R-free MTZ is reconciled back to it; flagged files and the lock-as-slot case are untouched |
| `tst_input_mtz_has_rfree_plumbing.py` | 4 | input_mtz_has_rfree tri-state survives the client→server wire (False not dropped) |
| `tst_event_system.py` | 13 | Event logging system |
| `tst_program_registration.py` | 13 | Program registry tests |
| `tst_integration.py` | 13 | End-to-end workflow tests |
| `tst_file_utils.py` | 12 | Shared file classification utilities |
| `tst_session_directives.py` | 12 | Session-level directive handling |
| `tst_state_serialization.py` | 12 | State packaging/unpackaging |
| `tst_summary_display.py` | 12 | Summary formatting |
| `tst_advice_preprocessing.py` | 11 | README discovery, advice processing |
| `tst_metric_patterns.py` | 11 | Log parsing patterns |
| `tst_history_analysis.py` | 10 | History analysis, anomalous workflow support |
| `tst_session_summary.py` | 12 | Session summary generation, STOP cycle exclusion, working directory in output, failure diagnosis reference |
| `tst_yaml_tools.py` | 9 | YAML validation and inspection |
| `tst_docs_tools.py` | 8 | Documentation generation |
| `tst_dry_run.py` | 8 | Dry run manager functionality |
| `tst_session_tools.py` | 7 | Session management utilities |
| `tst_template.py` | 5 | Template builder |
| `tst_phaser_multimodel.py` | 3 | Phaser multi-model handling |
| `tst_utils.py` | 2 | Assert helpers |
| `tst_v112_13_fixes.py` | — | Companion files, intermediate filtering, file categorisation (v112.13) |
| `tst_audit_fixes.py` | 230 | Audit regressions: `_is_failed_result`, zombie state, xtriage resolution, RSR map_cc, max_refine_cycles landing (v112.14); session management keywords, `get_results()` safety, restart_mode auto-set (v112.31 P1/P3/P4); completed-workflow extension via `advice_changed` phase step-back (v112.31 Q1); S2L probe crash detection (v112.52); diagnosable terminal errors — detection, HTML report, no-Sorry UX (v112.55–56); HETATM-based ligand detection, noligand false-positive fix (v112.57–59); Results page working directory + failure diagnosis reference, fatal-diagnosis skip (v112.60–62); crystal symmetry regex fallback, scoped PHIL form (v112.64); ligandfit file selection: refine MTZ classification, word-boundary patterns, content-based PDB guards, supplemental file discovery, duplicate detection (v112.70) |
| `tst_hardcoded_cleanup.py` | 36 | Hardcoded categorizer cleanup tests |

Total: **1080+ tests across 38 files**

### Key Tests for Recent Fixes

| Test | File | What It Verifies |
|------|------|------------------|
| `test_summary_includes_working_directory` | tst_session_summary.py | `_extract_summary_data` includes `working_dir`; markdown renders `**Working directory:**` |
| `test_failure_diagnosis_path_in_summary` | tst_session_summary.py | "Failure Diagnosis" section appears when `failure_diagnosis_path` set; absent on normal runs |
| `test_normalize_unit_cell_parenthesised` | tst_directive_extractor.py | `_normalize_unit_cell` converts `(a, b, c, al, be, ga)` to space-separated string |
| `test_validate_directives_normalises_unit_cell` | tst_directive_extractor.py | `validate_directives` normalises parenthesised unit_cell in program_settings |
| `test_extract_crystal_symmetry_simple_parenthesised` | tst_directive_extractor.py | Rules extractor captures unit cell from prose like "unit cell (116.097, …)" |
| `test_extract_directives_simple_unit_cell_end_to_end` | tst_directive_extractor.py | Full extraction from real nsf-d2 preprocessed advice (the AIAgent_35 failure case) |
| `test_unit_cell_in_valid_settings` | tst_directive_extractor.py | `unit_cell` and `space_group` present in `VALID_SETTINGS` with `str` type |
| `test_s3a_detect_crystal_symmetry_mismatch` | tst_audit_fixes.py | DiagnosisDetector matches crystal symmetry mismatch pattern |
| `test_s3a_build_html_new_fields` | tst_audit_fixes.py | Diagnosis HTML: heading is "Error diagnosis", file path, job name, and working dir are shown |
| `test_s3a_diagnose_returns_true_no_sorry` | tst_audit_fixes.py | `_diagnose_terminal_failure` returns True (no Sorry); `_finalize_session` skips Results page on fatal diagnosis |
| `test_s3a_finalize_runs_after_diagnosis` | tst_audit_fixes.py | `_finalize_session` runs unconditionally even when diagnosis fires; Results summary suppressed |
| `test_pdb_is_small_molecule_helper` | tst_audit_fixes.py | `_pdb_is_small_molecule` correctly identifies HETATM-only files as small molecules |
| `test_hetcode_ligand_not_used_as_refine_model` | tst_audit_fixes.py | atp.pdb / gdp.pdb (HETATM-only) not selected as refinement model |
| `test_is_ligand_file_noligand_false_positive` | tst_audit_fixes.py | `nsf-d2_noligand.pdb` not classified as ligand; word-boundary regex works correctly |
| `test_s4b_fallback_populates_unit_cell_from_empty_directives` | tst_audit_fixes.py | `_apply_crystal_symmetry_fallback` extracts unit_cell when LLM returns {} |
| `test_s4b_fallback_does_not_overwrite_llm_unit_cell` | tst_audit_fixes.py | Fallback preserves an existing LLM-extracted unit_cell, never overwrites |
| `test_s4b_fallback_populates_space_group_only_when_missing` | tst_audit_fixes.py | Partial fill: space_group added without disturbing existing unit_cell |
| `test_s4b_program_registry_uses_crystal_symmetry_scope` | tst_audit_fixes.py | `program_registry` PASSTHROUGH emits `crystal_symmetry.unit_cell=` form |
| `test_s4b_inject_crystal_symmetry_uses_scoped_form` | tst_audit_fixes.py | `_inject_crystal_symmetry` appends `crystal_symmetry.unit_cell=` not bare form |
| `test_s4b_fallback_called_in_extract_directives` | tst_audit_fixes.py | `_apply_crystal_symmetry_fallback` called after `validate_directives` in source |
| `tst_automation_path_in_workflow_state` | tst_workflow_state.py | automation_path correctly set in state |
| `tst_stepwise_mode_blocks_predict_and_build_after_prediction` | tst_workflow_state.py | predict_and_build blocked in stepwise mode after prediction |
| `tst_autobuild_beats_earlier_refine_with_better_metrics` | tst_best_files_tracker.py | Autobuild with better R-free becomes best model |
| `tst_yaml_stage_scores` | tst_best_files_tracker.py | autobuild_output has same score as refined (100) |
| `tst_stop_cycle_excluded_from_count` | tst_session_summary.py | STOP cycles not counted in total_cycles |
| `tst_cryoem_done_flags` | tst_workflow_state.py | cryo-EM done flags (resolve_cryo_em, map_sharpening, etc.) |
| `tst_failed_refine_not_counted` | tst_workflow_state.py | Failed refinements don't increment refine_count |
| `tst_mr_sad_after_phaser_with_anomalous` | tst_workflow_state.py | MR-SAD state after phaser + anomalous data |
| `tst_mr_sad_not_triggered_without_anomalous` | tst_workflow_state.py | MR-SAD not triggered without anomalous signal |
| `tst_mr_sad_not_triggered_when_autosol_done` | tst_workflow_state.py | MR-SAD skipped if autosol already ran |
| `tst_mr_sad_directive_prioritizes_phaser` | tst_workflow_state.py | use_mr_sad prioritizes phaser, removes autosol from obtain_model |
| `tst_normal_sad_still_works` | tst_workflow_state.py | Normal SAD workflow unaffected by MR-SAD changes |
| `tst_autosol_has_partial_model_config` | tst_workflow_state.py | autosol YAML has partpdb_file optional input |
| `tst_mr_sad_directive_overrides_no_anomalous` | tst_workflow_state.py | use_mr_sad triggers MR-SAD even without has_anomalous |
| `tst_experimental_phasing_yaml_structure` | tst_workflow_state.py | workflows.yaml experimental_phasing phase structure |
| `tst_predict_and_build_blocked_after_full_completion` | tst_workflow_state.py | predict_and_build not re-offered after full run, even with directives |
| `tst_autosol_partpdb_file_in_command` | tst_command_builder.py | autosol MR-SAD command includes partpdb_file=PHASER.1.pdb |
| `tst_model_scoring_prefers_refined` | tst_best_files_tracker.py | Refined models scored higher than MR output |
| `tst_predicted_model_exclusion` | tst_best_files_tracker.py | Predicted models excluded from model category |
| `test_j2_is_failed_result_false_positives_eliminated` | tst_audit_fixes.py | "No ERROR detected" and similar must not suppress done flags |
| `test_j5_zombie_cleared_when_output_missing` | tst_audit_fixes.py | resolve_cryo_em_done cleared when denmod_map.ccp4 missing |
| `test_j5_refine_zombie_decrements_count` | tst_audit_fixes.py | refine_done zombie decrements refine_count (not just clears flag) |
| `test_e1_e2_xtriage_resolution_dash_separator` | tst_audit_fixes.py | "50.00 - 2.30" format extracts 2.30 (not 50.0) |
| `test_e2_xtriage_resolution_anchor_blocks_completeness_line` | tst_audit_fixes.py | "Completeness in resolution range: 1" must not match resolution pattern |
| `test_e1_rsr_map_cc_uses_last_cycle` | tst_audit_fixes.py | RSR map_cc returns final macro-cycle value |
| `test_i1_max_refine_cycles_xray_controlled_landing` | tst_audit_fixes.py | Hitting max_refine_cycles injects validate + STOP, not bare STOP |
| `test_i1_max_refine_cycles_cryoem_uses_rsr_count` | tst_audit_fixes.py | Cryo-EM limit uses rsr_count, not refine_count |
| `test_i2_after_program_beats_quality_gate` | tst_audit_fixes.py | At step=refine with at_target=True and after_program done, refine programs are removed and STOP added.  Docstring updated (v116.x) to reflect that this is upheld by the at_target gate; the explicit-wipe-on-stop only fires when `stop_after_requested=True`. Test passes unchanged. |
| `test_i2_v112_78_after_program_min_run_not_hard_stop` | tst_audit_fixes.py | v116.x regression test for Tom's nsf-d2-ligand ligandfit case: after_program done at refine step with NO `stop_after_requested` (plan progression hint) must leave ligandfit available, NOT add STOP, NOT wipe to [STOP]. |
| `test_i2_v112_78_user_explicit_stop_still_offers_STOP` | tst_audit_fixes.py | v116.x complementary test: user-explicit "X and stop" (`stop_after_requested=True`) wipes `valid_programs` to `[STOP]` once X has run. |
| `test_p1_handle_session_management_display_sets_result` | tst_audit_fixes.py | `display_and_stop` populates `self.result` via `_finalize_session(skip_summary=True)` |
| `test_p2_remove_last_n_saves_and_populates_result` | tst_audit_fixes.py | `remove_last_n` saves to disk AND updates `result.session_data` |
| `test_p3_get_results_safe_before_run` | tst_audit_fixes.py | `get_results()` on fresh instance returns None, not AttributeError |
| `test_p4_restart_mode_set_to_resume_for_display` | tst_audit_fixes.py | `display_and_stop` auto-sets `restart_mode=resume` before `set_defaults()` |
| `test_q1_complete_phase_has_only_stop` | tst_audit_fixes.py | Baseline: completed workflow → `valid_programs=['STOP']` (no Q1 fix applied) |
| `test_q1_advice_changed_steps_back_to_validate` | tst_audit_fixes.py | Core logic: `advice_changed` + `phase=complete` → validate-phase programs incl. polder |
| `test_q1_no_step_back_when_advice_unchanged` | tst_audit_fixes.py | Unchanged advice leaves `complete` phase as-is (no spurious step-back) |
| `test_q1_cryoem_complete_phase_steps_back` | tst_audit_fixes.py | Same step-back logic works for cryo-EM `complete` phase |
| `test_q1_polder_reruns_allowed_when_already_done` | tst_audit_fixes.py | `polder_done=True` does NOT block polder re-run (no `run_once` strategy) |
| `test_q1_advice_cleared_after_one_cycle` | tst_audit_fixes.py | `advice_changed` flag cleared after one successful cycle; next cycle terminates normally |
| `test_q1_step_back_does_not_apply_outside_complete` | tst_audit_fixes.py | Guard fires only on `complete` phase; `refine`, `validate`, `obtain_model` untouched |
| `test_s5j_refine_mtz_classified_as_map_coeffs` | tst_audit_fixes.py | `refine_001.mtz` classified as `map_coeffs_mtz` (not `data_mtz`); regex + `classify_mtz_type` |
| `test_s5j_matches_exclude_pattern_word_boundary` | tst_audit_fixes.py | Word-boundary matching: `ligand` matches `atp_ligand.pdb` but NOT `nsf-d2_noligand.pdb` |
| `test_s5j_content_guard_ligand_and_model_slots` | tst_audit_fixes.py | HETATM-only PDB rejected from model slot; ATOM-bearing PDB rejected from ligand slot |
| `test_s5j_refine_cif_excluded_from_ligand_slot` | tst_audit_fixes.py | `refine_001.cif` (geometry restraints) excluded from ligandfit ligand slot |
| `test_s5j_session_rebuild_map_coeffs` | tst_audit_fixes.py | `_rebuild_best_files_from_cycles` populates `map_coeffs_mtz` from `refine_001.mtz` |
| `test_s5j_session_rebuild_supplemental_map_coeffs` | tst_audit_fixes.py | Supplemental file discovery on session load finds `refine_001.mtz` from `refine_001_data.mtz` |
| `test_s5j_record_result_discovers_supplemental_map_coeffs` | tst_audit_fixes.py | `record_result` discovers supplemental map coefficients on live path |
| `test_s5j_duplicate_detection_different_model_not_duplicate` | tst_audit_fixes.py | `phenix.refine` with different model file NOT flagged as duplicate despite high param overlap |

### Writing New Tests

### 1. Create Test File

```bash
touch tests/tst_my_feature.py
```

### 2. Add Boilerplate

```python
"""
Tests for my_feature module.

Run with: python tests/tst_my_feature.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_in,
    run_tests_with_fail_fast
)

from agent.my_feature import my_function


# =============================================================================
# MY FEATURE TESTS
# =============================================================================

def test_my_function_basic():
    """Test basic functionality."""
    result = my_function("input")
    assert_equal(result, "expected")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
```

### 3. Run Tests

```bash
python tests/tst_my_feature.py
```

### Handling Optional Dependencies

Some tests require PHENIX/libtbx. Use runtime checks to skip gracefully:

```python
# Try to import - may fail without libtbx
try:
    from agent.workflow_engine import WorkflowEngine
    HAS_WORKFLOW_ENGINE = True
except ImportError:
    HAS_WORKFLOW_ENGINE = False
    WorkflowEngine = None


def test_workflow_engine_feature():
    """Test that requires workflow engine."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return
    
    engine = WorkflowEngine()
    result = engine.some_method()
    assert_equal(result, expected)
```

### Test Categories

### Unit Tests (No External Dependencies)

These test individual functions in isolation:
- `tst_directive_extractor.py`
- `tst_transport.py`
- `tst_best_files_tracker.py`

### Integration Tests (May Need PHENIX)

These test component interactions:
- `tst_directives_integration.py`
- `tst_decision_flow.py`
- `tst_session_directives.py`

### Scenario Tests (Full Workflows)

Located in `tests/scenarios/`, these test complete workflows:
- `dry_run_xray_basic.py`
- `dry_run_cryoem_basic.py`

### Continuous Integration

The test suite is designed for CI environments:

```bash
# Quick CI run (no PHENIX required)
python tests/run_all_tests.py --quick

# Full CI run (with PHENIX)
python tests/run_all_tests.py
```

Exit codes:
- `0`: All tests passed
- `1`: One or more tests failed

### Debugging Test Failures

When a test fails, you get a full traceback:

```
  test_something ... FAIL
Traceback (most recent call last):
  File "tests/tst_module.py", line 42, in <module>
    run_all_tests()
  File "tests/tst_utils.py", line 357, in run_tests_with_fail_fast
    func()
  File "tests/tst_module.py", line 25, in test_something
    assert_equal(result, expected)
  File "tests/tst_utils.py", line 47, in assert_equal
    raise AssertionError(msg)
AssertionError: expected 'foo', got 'bar'
```

To debug interactively:

```python
# Add this before the failing assertion
import pdb; pdb.set_trace()
```

### Best Practices

1. **One assertion focus per test** - Test one specific behavior
2. **Descriptive names** - `test_extract_resolution_from_angstrom_notation`
3. **Clear docstrings** - Explain what the test verifies
4. **Minimal setup** - Keep tests fast and independent
5. **Section headers** - Group related tests with `# ===` comments
6. **Edge cases** - Test None, empty, boundary conditions

### Migration from unittest

The codebase previously used `unittest.TestCase`. If you encounter old-style tests:

**Before (unittest style):**
```python
class TestSomething(unittest.TestCase):
    def setUp(self):
        self.obj = create_object()
    
    def test_feature(self):
        self.assertEqual(self.obj.method(), "expected")
```

**After (cctbx style):**
```python
def test_feature():
    """Test the feature."""
    obj = create_object()
    assert_equal(obj.method(), "expected")
```

The test runner supports both styles during migration, but new tests should use plain functions.

---

## 6. Safety Checks

This document is auto-generated by `agent/generate_safety_docs.py`.

### Overview

The PHENIX AI Agent includes multiple layers of safety checks:

| Category | Count |
|----------|-------|
| Sanity Checks (Pre-Execution) | 20 |
| Model Placement Gate (v114.1) | 1 |
| Directive Validation (Post-LLM) | 7 |
| File Validation | 4 |
| Workflow State Validation | 8 |
| Command Building Validation | 3 |
| Input Validation | 29 |
| Post-Processing Corrections | 4 |

### Model Placement Gate (v114.1, hardened v114.2)

When the agent detects that a model already fits the data
(via `model_vs_data` CC > 0.3 or `refine` R-free < 0.50),
it locks `model_is_placed` in the session and suppresses
destructive programs (phaser, autosol, predict_and_build)
from `valid_programs`. This prevents the LLM from running
unnecessary molecular replacement on pre-solved structures.

Pending MR/phasing plan stages are marked "skipped" (⊘).
When the model is well-refined (R-free < 0.35), the
`model_rebuilding` stage is also skipped.
A conflict warning is logged if user advice mentions MR/phaser.

**v114.2 additions:**
- Result text parsing fallback (R-free/R-work from log output
  when metrics dict is empty)
- RSR symmetry mismatch → `needs_dock` (failed
  real_space_refine routes to dock_in_map)
- Dock keywords in advice cancel placement assumption
- Plan enforcement: `plan_has_pending_stages` suppresses
  premature STOP when plan has more work to do

### Sanity Check Thresholds

The `repeated_failures` check triggers after **4** consecutive
identical failures of the same program (raised from 3 in v114.1
to allow reference model refinement and similar tutorials that
need a first failure to learn correct parameters).

### Sanity Check Issues

These are the specific issues that can be detected by the SanityChecker:

| Code | Severity | Message |
|------|----------|---------|
| `no_data_for_workflow` | critical | No reflection data file found for X-ray workflow |
| `no_data_for_workflow` | critical | No map file found for cryo-EM workflow |
| `no_model_for_refine` | critical | Refinement requested but no model file available |
| `search_model_not_positioned` | critical | Cannot refine: search model found but not yet positioned |
| `resolution_unknown` | warning | Entering refinement phase without established resolution |

### Sanity Check Functions

### `agent/directive_validator.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `ValidationResult` | 58 | class | Result of directive validation. |

### `agent/sanity_checker.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `SanityChecker` | 67 | class | Checks for impossible or nonsensical agent states.  Loads check definitions from... |
| `SanityIssue` | 32 | class | A single sanity check failure. |
| `SanityResult` | 46 | class | Result of sanity checking. |
| `_check_data_exists` | 286 | function | Check that required experimental data exists for the workflow. |
| `_check_experiment_type_stable` | 185 | function | Check that experiment type hasn't changed unexpectedly. |
| `_check_metric_anomalies` | 392 | function | Check for dramatic metric changes that indicate problems. |
| `_check_model_for_refine` | 203 | function | Check that a POSITIONED model exists when entering refine state.  With semantic ... |
| `_check_repeated_failures` | 315 | function | Check for repeated identical failures. |
| `_check_resolution_established` | 377 | function | Warn if entering refinement without established resolution. |
| `check` | 99 | function | (No docstring) |
| `experiment_type_changed` | 192 | sanity_issue | Severity: critical |
| `map_cc_drop` | 425 | sanity_issue | Severity: warning |
| `no_data_for_workflow` | 296 | sanity_issue | Severity: critical |
| `no_data_for_workflow` | 305 | sanity_issue | Severity: critical |
| `no_model_for_refine` | 274 | sanity_issue | Severity: critical |
| `r_free_spike` | 409 | sanity_issue | Severity: warning |
| `repeated_failures` | 345 | sanity_issue | Severity: critical |
| `resolution_unknown` | 383 | sanity_issue | Severity: warning |
| `search_model_not_positioned` | 260 | sanity_issue | Severity: critical |

### Directive Validation (Post-LLM)

### `agent/directive_extractor.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `check_stop_conditions` | 939 | function | Check if any stop conditions are met.  Args: directives: Directives dict cycle_n... |
| `validate_directives` | 547 | function | Validate and clean extracted directives.  Removes invalid entries, converts type... |

### `agent/directive_validator.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_max_program_cycles` | 1023 | function | Check if max cycles for a specific program have been reached.  Args: directives:... |
| `_check_stop_conditions` | 955 | function | Check all stop conditions from directives.  Args: directives: Directives dict cy... |
| `check_program_available` | 708 | function | Quick check if a specific program is available. |
| `validate_directives` | 505 | function | (No docstring) |
| `validate_intent` | 757 | function | Validate and modify an LLM intent against user directives.  This applies program... |

### File Validation

### `agent/best_files_tracker.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_is_intermediate_file` | 1120 | function | Check if a file is an intermediate that shouldn't be tracked.  Args: path: File ... |

### `agent/yaml_tools.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_validate_file_categories` | 293 | function | Validate file_categories.yaml structure. |
| `validate_files` | 1295 | function | Validate YAML files and report issues. |

### `knowledge/yaml_loader.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate_yaml_files` | 673 | function | Validate that all YAML files are present and well-formed.  Returns: tuple: (is_v... |

### Workflow State Validation

### `agent/program_validator.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `check_workflow_engine` | 208 | function | Check if program is in workflow engine context (if needed). |
| `check_workflow_state` | 180 | function | Check if program has history tracking. |
| `check_workflows_yaml` | 99 | function | Check if program is included in workflows.yaml. |

### `agent/workflow_engine.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_conditions` | 852 | function | Check if program conditions are met. |
| `_check_metric_condition` | 878 | function | Check a metric condition like "> 0.35" or "< target_r_free". |
| `_check_priority_condition` | 1047 | function | Check if a priority_when condition is satisfied.  Args: condition: Condition str... |

### `agent/workflow_state.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate_program_choice` | 683 | function | Validate that a program choice is allowed in the current state.  Args: chosen_pr... |

### `agent/yaml_tools.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_validate_workflows` | 594 | function | Validate workflows.yaml structure. |

### Command Building Validation

### `agent/template_builder.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_invariant` | 708 | function | Check if an invariant condition is satisfied.  Supports: - has_file: [list of fi... |
| `validate_and_fix` | 649 | function | Validate program invariants and apply fixes if needed.  This is the single place... |
| `validate_intent` | 532 | function | Validate that an intent's files exist in the available files list.  Args: progra... |

### Input Validation

### `agent/error_analyzer.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_retry_limits` | 531 | function | Check if we've exceeded retry limits for this error type.  Returns (can_retry, r... |

### `agent/event_log.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `has_errors` | 268 | function | Check if any error events were recorded. |

### `agent/graph_nodes.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate` | 1984 | function | Node D: Validate command before execution.  Checks: 1. Command is not empty 2. P... |
| `validate_provider` | 152 | function | Validate that the provider is supported and properly configured.  Returns: tuple... |

### `agent/pattern_manager.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate_all` | 454 | function | Validate all patterns against their test cases.  Returns: dict: {pattern_name: (... |
| `validate_pattern` | 389 | function | Validate a pattern against its test cases.  Args: name: Pattern name  Returns: t... |

### `agent/program_validator.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `check_directive_extractor` | 261 | function | Check if program has tutorial/directive patterns. |
| `check_log_parsers` | 134 | function | Check if program has log parsing support. |
| `check_programs_yaml` | 55 | function | Check if program is defined in programs.yaml. |
| `check_session_summary` | 236 | function | Check if program metrics appear in session summary. |
| `validate_all` | 360 | function | Validate all programs. |
| `validate_program` | 286 | function | Validate that a program is fully configured. |

### `agent/rate_limit_handler.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_decay` | 101 | function | Check if enough time has passed to decay the delay. |

### `agent/session.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_validate_directives` | 302 | function | Validate directives against available capabilities.  Checks if the user is reque... |
| `check_directive_stop_conditions` | 393 | function | Check if any directive stop conditions are met.  Args: cycle_number: Current cyc... |
| `check_max_program_cycles` | 504 | function | Check if max cycles for a program have been reached per directives.  Args: progr... |

### `agent/transport.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `verify_roundtrip` | 865 | function | Verify that a request survives the encode/decode roundtrip.  This is useful for ... |

### `agent/yaml_tools.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_validate_check_recursive` | 578 | function | Recursively validate check structure (handles any_of, all_of). |
| `_validate_metrics` | 703 | function | Validate metrics.yaml structure. |
| `_validate_programs` | 359 | function | Validate programs.yaml structure. |
| `validate_yaml_structure` | 260 | function | Validate the structure of a YAML file based on its type.  Returns: list: List of... |

### `knowledge/api_schema.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate_request` | 304 | function | Validate a v2 request.  Args: request: Dict containing request data strict: If T... |
| `validate_response` | 350 | function | Validate a v2 response.  Args: response: Dict containing response data strict: I... |

### `programs/ai_agent.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_check_benign_error` | 2164 | function | Check if error is benign. |
| `_check_for_fatal_error` | 1902 | function | Check if error text contains fatal error markers. |
| `_validate_command` | 1814 | function | Validate that a command is safe to execute.  Only allows commands that start wit... |
| `_validate_directives` | 2645 | function | Validate directives against available capabilities.  If validation fails (user r... |
| `validate` | 610 | function | (No docstring) |

### `programs/ai_analysis.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `validate` | 321 | function | (No docstring) |

### Post-Processing Corrections

### `agent/directive_extractor.py`

| Name | Line | Type | Description |
|------|------|------|-------------|
| `_fix_ligand_workflow_conflict` | 700 | function | Fix conflicting directives when ligand fitting workflow is detected.  If user wa... |
| `_fix_ligand_workflow_conflict` | 700 | function | Fix conflicting directives when ligand fitting workflow is detected.  If user wa... |
| `_fix_program_name` | 738 | function | Try to fix common variations in program names.  Args: name: Potentially incorrec... |
| `_fix_program_name` | 738 | function | Try to fix common variations in program names.  Args: name: Potentially incorrec... |

---

## 7. Validation Subsystems

### Overview

The PHENIX AI Agent includes a sanity-checking system that detects when the agent
is in an impossible or nonsensical state. When critical issues are detected, the
agent aborts gracefully with a clear explanation and suggestions for resolution.

### Design Principles

1. **Abort only when truly necessary** - Red flags indicate real problems, not edge cases
2. **User control** - Configurable abort behavior
3. **Clear, actionable messages** - Explain what went wrong and how to fix it
4. **Resumable sessions** - After fixing issues, users can continue the session

### Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `abort_on_red_flags` | `True` | Abort on critical issues (experiment type change, impossible states) |
| `abort_on_warnings` | `False` | Also abort on warnings (metric anomalies, missing resolution) |

### Red Flag Categories

### Critical Issues (abort_on_red_flags)

| Code | Check | Trigger |
|------|-------|---------|
| `experiment_type_changed` | Experiment type stability | Type differs from initial (xray↔cryoem) |
| `no_model_for_refine` | Model required for refinement | Refine state reached with no PDB files |
| `no_data_for_workflow` | Data required for workflow | No MTZ (xray) or map (cryoem) |
| `repeated_failures` | Same error repeated | Same program fails 3+ times with identical error |

### Warnings (abort_on_warnings)

| Code | Check | Trigger |
|------|-------|---------|
| `resolution_unknown` | Resolution should be established | Entering refine phase without resolution |
| `r_free_spike` | R-free anomaly | R-free increases by >0.15 in one cycle |
| `map_cc_drop` | Map CC anomaly | Map CC drops by >0.3 in one cycle |

### Validation Gates

Before the workflow can stop with success, validation must be completed:

1. **Quality Check Required**: Must run `phenix.molprobity` before stopping
2. **Metrics Must Be Acceptable**: R-free or map CC must meet thresholds
3. **No Red Flags**: No critical issues detected

### Validation Gate Logic

```
if workflow_should_stop:
    if not validation_completed:
        return "Run molprobity first"
    if red_flags:
        return "Cannot stop: " + red_flag_message
    return "STOP"
```

### Implementation

### SanityChecker Class

Located in `agent/sanity_checker.py`:

```python
from agent.sanity_checker import SanityChecker

checker = SanityChecker()

# Run all checks
result = checker.check(
    context={
        "experiment_type": "xray",
        "state": "xray_refined",
        "has_model": True,
        "has_data_mtz": True,
        "history": [...],
        "metrics_history": [...],
    },
    session_info={...},
    abort_on_red_flags=True,
    abort_on_warnings=False
)

if result.red_flags:
    print("Critical issues:", result.red_flags)
if result.warnings:
    print("Warnings:", result.warnings)
```

### YAML Configuration

Sanity checks are defined in `knowledge/workflows.yaml`:

```yaml
sanity_checks:
  critical:
    - name: experiment_type_stable
      message: "Experiment type changed from {initial} to {current}"
      suggestion: "Start a new session for different experiment types"

    - name: model_required_for_refine
      message: "No model file found but refinement is required"
      suggestion: "Check that model building succeeded"

  warnings:
    - name: resolution_established
      message: "Resolution not established before refinement"
      suggestion: "Run xtriage or mtriage first"
```

### Response Format

When validation issues are detected, they appear in the response:

```json
{
  "api_version": "2.0",
  "stop": true,
  "stop_reason": "red_flag",
  "metadata": {
    "red_flags": [
      {
        "code": "experiment_type_changed",
        "message": "Experiment type changed from xray to cryoem",
        "suggestion": "Start a new session for different experiment types"
      }
    ],
    "warnings": []
  }
}
```

### Testing

Run sanity checker tests:

```bash
python tests/tst_sanity_checker.py
```

Test scenarios include:
- Experiment type change detection
- Missing model detection
- Repeated failure detection
- Metric spike detection
- Validation gate enforcement

---

## 8. Backward Compatibility & Contract

### Golden Rule

**A user with yesterday's PHENIX install must be able to connect to
today's server and get correct results.** Their client will send the same
`session_info` fields it sent yesterday. The server must accept that
input, produce valid commands, and return a well-formed response. If you
cannot maintain this guarantee for a particular change, you MUST bump the
protocol version and handle old clients gracefully (Rule 6 below).

This applies to every kind of change:

| Change type | Safe? | Why |
|-------------|-------|-----|
| `programs.yaml` — add program, change conditions, edit flags | **Yes** | Server-side knowledge; client never reads it |
| `workflows.yaml` — reorder steps, change conditions | **Yes** | Server-side workflow; client never reads it |
| `command_builder.py` — change file selection, dedup logic | **Yes** | Runs on server only; client receives the final command string |
| `graph_nodes.py` — change planning, validation, guards | **Yes** | Runs on server only |
| `recoverable_errors.yaml` — add/change error recovery | **Yes** | Server-side knowledge |
| `prompts_hybrid.py` — change LLM prompts | **Yes** | Server-side only |
| `session_info` — add a new field server reads | **Safe IF** default provided | Old clients won't send it; server must `.get(field, default)` |
| `session_info` — remove a field server reads | **UNSAFE** | Old clients still send it; server ignoring it changes behavior |
| `session_info` — change meaning of existing field | **UNSAFE** | Old clients send the old meaning |
| Response — remove a field client reads | **UNSAFE** | Old clients will KeyError or get wrong defaults |
| Response — add a new field | **Yes** | Old clients ignore unknown fields |
| `agent/` shared code — add new function | **Risky** | Runs on both client and server; old client doesn't have it |
| `agent/` shared code — change function signature | **UNSAFE** | Old client calls with old signature |

**UNSAFE changes require a protocol version bump:** increment
`CURRENT_PROTOCOL_VERSION` in `agent/contract.py`, add
`MIN_SUPPORTED_PROTOCOL_VERSION` handling for the transition period, and
log a deprecation warning for old clients via `get_deprecation_warnings()`.

### The Problem

Users run the PHENIX AI Agent client (GUI/CLI) on their local machines. The agent's "brain" runs on a remote REST server. When you update the server, some users will still be running older client versions. A server change that assumes a new field in `session_info` will silently break for every old client that doesn't send it.

This document describes the backward compatibility system that is now **implemented and enforced by automated tests**.

### Architecture

```
CLIENT (user's PHENIX)                  SERVER (your REST endpoint)
─────────────────────                   ──────────────────────────
programs/ai_agent.py                    agent/graph_nodes.py
  └── agent/session.py                    ├── perceive()    ← normalization + version gate
  └── agent/best_files_tracker.py         ├── plan()        ← stop guard
                                          ├── build()       ← stop guard
                                          │   └── agent/command_builder.py
                                          │   └── agent/workflow_state.py
                                          │   └── agent/workflow_engine.py
                                          ├── validate()    ← stop guard
                                          ├── fallback()    ← stop guard
                                          └── output_node() ← response field guarantees

  ── session_info, files, history ──→
  ←── history_record (command, etc) ──
```

The `agent/` directory is **shared code** — it ships with both the client (in PHENIX) and the server. When you update the server, the `agent/` code on the server is newer than the `agent/` code on the client. This is the primary source of compatibility risk.

### The Contract (`agent/contract.py`)

The contract registry is the **single source of truth** for the client↔server interface. It defines every field, its default value, and the protocol version it was introduced in.

### Client → Server (Request)

The client calls `decide_next_step()` with these arguments:

| Argument | Type | Description |
|----------|------|-------------|
| `log_content` | str | Log text from the last program run |
| `history` | list[dict] | Cycle records with `program`, `command`, `result`, `output_files`, `analysis` |
| `files` | list[str] | Available file paths (validated on client) |
| `guidelines` | str | User advice / project instructions |
| `session_resolution` | float or None | Resolution in Å |
| `session_info` | dict | **The main extensibility surface** — see below |
| `abort_on_red_flags` | bool | Whether to abort on sanity check failures |
| `abort_on_warnings` | bool | Whether to abort on warnings |

**`session_info` fields** (registered in `agent/contract.py :: SESSION_INFO_FIELDS`):

| Field | Default | Version | Description |
|-------|---------|---------|-------------|
| `experiment_type` | `""` | v1 | `"xray"` or `"cryoem"` |
| `best_files` | `{}` | v1 | `{category: path}` from BestFilesTracker |
| `rfree_mtz` | `None` | v1 | Locked R-free MTZ path |
| `directives` | `{}` | v1 | Parsed user directives |
| `rfree_resolution` | `None` | v2 | Resolution limit of R-free flags |
| `force_retry_program` | `None` | v2 | One-shot forced program for error recovery |
| `recovery_strategies` | `{}` | v2 | Error recovery hints |
| `explicit_program` | `None` | v2 | User-requested program from advice |
| `advice_changed` | `False` | v2 | Whether advice changed on resume |
| `bad_inject_params` | `{}` | v2 | Parameter injection blacklist |
| `unplaced_model_cell` | `None` | v3 | Pre-extracted CRYST1 cell `[a,b,c,α,β,γ]` |
| `model_hetatm_residues` | `None` | v3 | Pre-extracted HETATM data `[[chain,resseq,resname],...]` |
| `client_protocol_version` | `1` | v3 | Protocol version of the sending client |

### Server → Client (Response)

The server returns a `history_record` dict. `output_node()` guarantees these fields always exist with safe defaults:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `next_move` | dict | — | `{command, program, explanation, process_log}` |
| `debug_log` | list[str] | `[]` | Diagnostic messages from all nodes |
| `events` | list[dict] | `[]` | Structured events for display |
| `experiment_type` | str | — | Detected experiment type |
| `stop_reason` | str or None | `None` | Why the agent stopped |
| `abort_message` | str or None | `None` | Red flag abort explanation |
| `red_flag_issues` | list | `[]` | Issues that triggered abort |
| `warnings` | list[str] | `[]` | Deprecation/advisory messages |

### The `warnings` Deprecation Channel

The `warnings` field provides a server→client communication channel for non-fatal advisories. Currently it's used for protocol version deprecation:

- **Server side** (`perceive()` via `get_deprecation_warnings()`): Automatically appends a warning when `client_protocol_version < CURRENT_PROTOCOL_VERSION`.
- **Client side** (`ai_agent.py`): After receiving the response, prints each warning via `self.vlog.quiet("[AI Server Warning] ...")`.
- **Old clients** that don't check `warnings` simply ignore the field — they already use `.get()` / `getattr()` for all response fields.

### Protocol Version Constants

Defined in `agent/contract.py`:

| Constant | Value | Purpose |
|----------|-------|---------|
| `CURRENT_PROTOCOL_VERSION` | `5` | What the latest client sends |
| `MIN_SUPPORTED_PROTOCOL_VERSION` | `1` | Oldest client the server accepts |

The client imports `CURRENT_PROTOCOL_VERSION` from contract.py via `_get_protocol_version()` in `ai_agent.py`, so bumping the constant updates both server and client automatically.

### The Rules

### RULE 1: New `session_info` fields MUST have safe defaults

Every server-side `.get()` on a `session_info` field MUST supply an explicit default that matches the contract:

```python
# CORRECT
hetatm = session_info.get("model_hetatm_residues", None)
strategies = session_info.get("recovery_strategies", {})
changed = session_info.get("advice_changed", False)

# WRONG — crashes if old client doesn't send it
hetatm = session_info["model_hetatm_residues"]

# FRAGILE — returns None instead of contract default {}
strategies = session_info.get("recovery_strategies")
```

**Enforced by**: `tst_contract_compliance.py :: test_no_bare_session_info_bracket_reads` and `test_get_without_default_on_non_none_fields`. Both are hard assertions — the test suite fails if either is violated.

### RULE 2: Never REMOVE a response field the client reads

Old clients read specific fields from the response. If the server stops sending a field, old clients break. Fields can be added freely.

**Enforced by**: `output_node()` in `graph_nodes.py` guarantees `warnings`, `debug_log`, `events`, `stop_reason`, `abort_message`, and `red_flag_issues` always exist in the final state.

### RULE 3: Never change the SEMANTICS of existing fields

If `best_files.model` means "best positioned model path", don't repurpose it. Add a new field instead.

### RULE 4: New programs.yaml entries are safe; removing entries is not

Adding a new program to `programs.yaml` is safe. Removing a program that old clients might request via `explicit_program` or directives could cause confusing behavior.

When adding PHIL parameters for a program:
- **Individual params**: Add to `strategy_flags` with a short name, flag template, type, and hint. Use for standalone params (e.g., `ramachandran_restraints`).
- **Parameter namespaces**: Add to `allowed_phil_prefixes` instead of listing every sub-param. Use for PHIL scopes with many sub-parameters (e.g., `ncs.` covers `ncs.type`, `ncs.constraints`, `ncs.find_ncs`, etc.). The prefix match is case-insensitive and uses substring matching, so `"ncs."` in the prefix list matches both `ncs.type` and `refinement.pdb_interpretation.ncs.type`.
- **Blocked params**: Add to `_BLOCKED_PARAMS` in `phil_validator.py` for params that cause crashes. Blocked params are stripped even if they match a strategy_flag or prefix.

### RULE 5: File list is the client's responsibility

The server must never assume files exist on its filesystem. All file discovery, validation, and content reading must happen on the client, with results passed through `session_info`.

### RULE 6: Hard rejection for unsupported clients

When a client is too old to serve correctly, the server returns a clean STOP with a human-readable error message instead of silently producing wrong results.

**Implemented in**: `perceive()` calls `check_client_version()` from `contract.py`. If `client_protocol_version < MIN_SUPPORTED_PROTOCOL_VERSION`, the pipeline returns immediately with `stop=True`, `stop_reason="unsupported_client"`, and an `abort_message` telling the user to update.

Today `MIN_SUPPORTED_PROTOCOL_VERSION = 1` (accept everything). When you bump it, the `client_protocol_version` logging data tells you how many users would be affected.

### RULE 7: The `agent/` shared code trap

Code in `agent/` runs on both client and server. If the server adds a new function to `agent/workflow_engine.py` that old clients don't have, any code path that calls it on the client will crash.

Mitigations:
- Guard new utility functions with `try/except ImportError` when callable from client paths
- Prefer adding server-only logic in `graph_nodes.py` rather than modifying shared files
- Verify new dependencies exist in the PHENIX distribution

**Enforced by**: `tst_shared_code_imports.py :: test_no_forbidden_imports_in_shared` (blocks LLM SDK imports in shared code) and `test_agent_imports_reference_known_modules` (blocks imports of unknown modules).

### RULE 8: LocalAgent and RemoteAgent must be identical

`LocalAgent` (`phenix_ai/local_agent.py`) and `RemoteAgent` (`phenix_ai/remote_agent.py`) must build identical requests and return identical response formats. The only difference is transport: LocalAgent decodes the request locally, RemoteAgent sends it over HTTP.

When you add or change:
- A new `session_info` field → both agents receive it via `_query_agent_for_command()` (shared code, no action needed)
- A new `request["settings"]` entry → add it in **both** `decide_next_step()` methods
- A new field in the return `group_args` → add it in **both** return blocks

LocalAgent intentionally performs the full encode/decode roundtrip (`prepare_request_for_transport` → `process_request_from_transport`) even though it could pass the dict directly. This ensures local mode exercises the same transport code path as remote mode, catching serialization bugs before they reach production.

**Enforced by**: `tst_backward_compat.py` verifies structural invariants. Manual review is required for settings parity (no automated diff exists yet between the two `decide_next_step` methods).

### RULE 9: Numeric fields must survive JSON roundtrip everywhere

JSON round-tripping between client and server can turn floats into strings
(e.g. `0.385` → `"0.385"`). This affects StructureModel, metrics_history,
event data, and history analysis. Any code that performs arithmetic on
metric values (`r_free - r_work`, `previous - latest_r_free`) or formats
them (`"%.3f" % r_free`) must use `_safe_float()` or equivalent coercion.

**Files that require numeric coercion (all fixed as of v115.07+):**

| File | What to guard | Method |
|------|--------------|--------|
| `structure_model.py` | `model_state` fields (r_free, r_work, etc.) | `_coerce_numerics()` in `from_dict()` + `_safe_float()` at read sites |
| `metric_evaluator.py` | `metrics_history` r_free/CC extraction + arithmetic | `_safe_float()` at 5 sites |
| `metrics_analyzer.py` | `derive_metrics_from_history` + trend analysis | `_safe_float()` at 7 sites |
| `graph_nodes.py` | METRICS_EXTRACTED event emission (L788/791) | `float()` with try/except |
| `event_formatter.py` | Compact metrics `r_free - r_free_prev` (L925) | `float()` with try/except |
| `kb_tags.py` | `_trend_tags()` diffs and total_drop | `float()` with try/except |
| `workflow_state.py` | `_analyze_history()` metric reads (L1841-1850) | Local `_sf()` helper |

**When adding new numeric reads from history or metrics_history:**
1. Use `_safe_float()` or `float()` with try/except at every read site
2. Add a test in `tst_phase3_bug5.py` with a string-valued input
3. Never assume values from JSON dicts are the correct Python type

**Enforced by**: `tst_phase3_bug5.py :: test_bug4_*` tests (8 tests
covering string coercion, garbage handling, crash prevention, gap
detection, None preservation, metric_evaluator, cryo-EM CC, and kb_tags).

### RULE 10: When adding a new LLM-driven extractor, plumb raw advice alongside processed

**Background.** The advice preprocessor sits between the user's raw
text and any downstream LLM that consumes that text.  The preprocessor
LLM occasionally mangles short imperative input — the openai-provider
"density modify and stop" failure case is documented in
`CHANGELOG.md` under Step 1 / v117.  When a downstream LLM
extractor sees only the preprocessed text, the user's original intent
can be silently lost.

**The rule.** When you write a new LLM-driven extractor that consumes
user advice (any new prompt that takes `user_advice=` in
`directive_extractor`-shape), give it BOTH raw and processed inputs
with an explicit authority rule for the fields you care about.  Use
the Step 1 plumbing as a template:

1. The extractor function accepts `raw_advice=None` alongside
   `user_advice=` and selects between single-input and dual-input
   prompts based on `raw_advice is not None and raw_advice.strip()
   and raw_advice != user_advice`.
2. The dual-input prompt is derived from the single-input prompt via
   `str.replace()` of a named block, guarded by an import-time
   `assert ... .count(_SINGLE_INPUT_BLOCK) == 1` so future prompt edits
   that break the derivation fail at import.
3. Add the corresponding PHIL field
   (`user_advice_raw`-shape) to both `ai_analysis.py` (server) and
   `ai_agent.py` (client).  Wire it through `_build_server_args` and
   the local dispatcher.  Backward-compatible: defaults to `None` and
   old clients fall through to single-input.
4. Set the client-side PHIL value from `session.data["raw_advice"]`
   in the function that builds `directive_params`.  Use
   `session.data.get("raw_advice") or None` (None-safe per
   guideline §3) so an unset key or stored-as-empty doesn't put a
   bogus empty string on the wire.

**Files that establish the template:**

| Role | File | Function/symbol |
|------|------|-----------------|
| Prompt + selector | `agent/directive_extractor.py` | `extract_directives()`, `_RAW_INPUT_BLOCK`, `DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW` |
| Module wrapper | `phenix_ai/run_ai_analysis.py` | `run_directive_extraction(raw_user_advice=...)` |
| Session wrapper | `agent/session.py` | `AgentSession.extract_directives(raw_advice=...)` |
| Server PHIL + dispatch | `phenix/programs/ai_analysis.py` | `user_advice_raw` PHIL, `_build_server_args`, `_run_directive_extraction_locally` |
| Client touch + PHIL | `phenix/programs/ai_agent.py` | `user_advice_raw` PHIL, `_extract_directives` setter |

**Why not just read raw from somewhere downstream?** Because the
agent's preprocessing flow is what produces `session.data["raw_advice"]`,
and that data is client-side state at the time we build PHIL.  The
clean way to get it to the server is over the wire as a PHIL field;
adding it to `session_info` (which goes elsewhere) would conflate two
different transport channels.

**Enforced by**: `tst_extract_raw_advice.py` source-scans every link in
the chain.  Run on pre-Step-1 code, every test in that module fails
with a diagnostic message pointing at the missing piece.

### What's Implemented

### Runtime Protection (server side)

All in `agent/graph_nodes.py`:

1. **`perceive()` — Normalization and version gate**: At the top of the pipeline, before any `session_info` access:
   - Calls `normalize_session_info()` — fills in safe defaults for every registered field (mutable defaults are copied to avoid cross-request contamination)
   - Calls `check_client_version()` — returns a clean STOP if client is below `MIN_SUPPORTED`
   - Calls `get_deprecation_warnings()` — appends advisory messages to `state["warnings"]`
   - Wrapped in `try/except ImportError` — graceful no-op if contract.py is somehow missing

2. **Stop guards on all 5 pipeline nodes**: Every node after `perceive()` begins with:
   ```python
   if state.get("stop"):
       return state
   ```
   This ensures a stop from any upstream node (version rejection, red flag, etc.) flows cleanly through without any node accidentally overwriting it or trying to build a command.

3. **`output_node()` — Response field guarantees**: Before returning the final state, ensures `warnings`, `debug_log`, `events`, `stop_reason`, `abort_message`, and `red_flag_issues` all exist with safe defaults.

### Client-Side Changes (`programs/ai_agent.py`)

1. **`client_protocol_version`**: Sent in `session_info`, imported from `contract.py` via `_get_protocol_version()` helper (falls back to `3` if import fails).

2. **`warnings` handling**: After receiving the server response, checks `history_record.get('warnings')` (with fallback to `result_record.warnings`) and prints each as `[AI Server Warning] ...`.

### Contract Registry (`agent/contract.py`)

Single file containing:
- `CURRENT_PROTOCOL_VERSION` and `MIN_SUPPORTED_PROTOCOL_VERSION` constants
- `SESSION_INFO_FIELDS` — every field with its default, version, and description
- `RESPONSE_FIELDS` and `NEXT_MOVE_FIELDS` — response shape documentation
- `HISTORY_ENTRY_FIELDS` — history list entry shape
- `STATE_PROMOTED_FIELDS` — fields copied from session_info to top-level state
- `normalize_session_info()` — runtime helper
- `check_client_version()` — version gate
- `get_deprecation_warnings()` — advisory message generator

Zero external dependencies — safe to import on both client and server.

### Test Suite

### `tests/tst_contract_compliance.py` — 10 tests

Static analysis and contract verification:

| Test | What it catches |
|------|----------------|
| `test_no_bare_session_info_bracket_reads` | `session_info["X"]` without `.get()` — would KeyError on old clients |
| `test_all_accessed_fields_in_contract` | Fields used but not registered in contract |
| `test_contract_defaults_consistency` | `.get("X", wrong_default)` — default doesn't match contract |
| `test_get_without_default_on_non_none_fields` | `.get("X")` where contract default is `{}`, `False`, `""` — fragile without normalization |
| `test_protocol_version_consistency` | Client sends a version that doesn't match `CURRENT_PROTOCOL_VERSION` |
| `test_version_bounds` | `MIN_SUPPORTED > CURRENT` (impossible state) |
| `test_warnings_in_response_contract` | `warnings` missing from `RESPONSE_FIELDS` |
| `test_client_handles_warnings` | Client code doesn't check for server warnings |
| `test_normalize_covers_all_fields` | `normalize_session_info()` misses a registered field |
| `test_normalize_mutable_isolation` | Mutable defaults leak between calls |

### `tests/tst_old_client_compat.py` — 6 tests

Frozen fixture tests against real old-client payloads:

| Test | What it catches |
|------|----------------|
| `test_fixtures_valid_structure` | Fixture JSON is well-formed |
| `test_normalize_fills_missing_fields` | Old client payloads get all fields filled in |
| `test_junk_fields_preserved` | Unknown fields aren't stripped (forward compat) |
| `test_version_check_accepts_all` | Current fixtures aren't rejected |
| `test_deprecation_warnings_for_old_clients` | Old clients get warnings, current clients don't |
| `test_full_pipeline_if_available` | Full perceive→output pipeline with old payloads (in PHENIX env) |

**Fixtures** (in `tests/fixtures/`, NEVER updated after creation):

| Fixture | Protocol | Fields | Scenario |
|---------|----------|--------|----------|
| `client_v1/xray_ligandfit.json` | v1 | 5 (missing 8) + junk field | X-ray after 2 refines, needs ligandfit |
| `client_v1/cryoem_initial.json` | v1 | 4 (missing 9) | Cryo-EM first cycle with half-maps |
| `client_v3/xray_ligandfit.json` | v3 | 13 (complete) | Same X-ray scenario, full fields |

### `tests/tst_shared_code_imports.py` — 4 tests

Import safety for shared `agent/` modules:

| Test | What it catches |
|------|----------------|
| `test_no_forbidden_imports_in_shared` | LLM SDKs (langchain, openai, anthropic) in shared code |
| `test_agent_imports_reference_known_modules` | Import of `agent.X` where X doesn't exist in PHENIX |
| `test_server_only_has_llm_imports` | Sanity check: `graph_nodes.py` should have LLM imports |
| `test_shared_imports_are_guarded` | Top-level unguarded intra-agent imports (warning-level) |

### Running All Tests

```bash
cd /path/to/agent_debug
PYTHONPATH=. python tests/tst_contract_compliance.py      # 10 tests
PYTHONPATH=. python tests/tst_old_client_compat.py        # 6 tests
PYTHONPATH=. python tests/tst_shared_code_imports.py      # 4 tests
PYTHONPATH=. python tests/tst_audit_fixes.py              # 31 tests (pre-existing)
```

### Version Lifecycle

The protocol version system provides a deprecation lifecycle:

```
 v1 shipped ──→ v2 shipped ──→ v3 shipped ──→ MIN bumped to v2
                                                │
                   soft warning ◄────────────────┘
                   for 6 months
                        │
                   hard rejection
                   after deadline
```

**To make a breaking change**:

1. Ship the new client (e.g., v4) with the new field/behavior
2. Server automatically warns clients below v4 (via `get_deprecation_warnings()`)
3. Monitor `client_protocol_version` logs to see who's still on old versions
4. Set a deadline and add it to the warning message
5. After the deadline, bump `MIN_SUPPORTED_PROTOCOL_VERSION` to v4
6. Old clients get a clear "please update" error, not silent breakage

### Checklist for Every Server Update

### Code Review

- [ ] No new `session_info["X"]` bracket reads — always `.get("X", default)`
- [ ] Every `.get("X", default)` uses the contract-registered default
- [ ] New session_info fields added to `contract.py :: SESSION_INFO_FIELDS`
- [ ] No removed response fields — old clients depend on them
- [ ] No changed field semantics — add new fields instead
- [ ] New `programs.yaml` entries are additive — no removals
- [ ] No new `os.path.exists()` guards on file paths — use `_file_is_available()`
- [ ] No new `open()` calls without fallback — server can't read client files
- [ ] New features degrade gracefully when session_info field is missing
- [ ] Shared `agent/` code: no new imports that don't exist in shipped PHENIX
- [ ] Version-gated features: if it needs a new client field, gate on `client_protocol_version`

### Tests

- [ ] `tst_contract_compliance.py` passes (10 tests)
- [ ] `tst_old_client_compat.py` passes (6 tests)
- [ ] `tst_shared_code_imports.py` passes (4 tests)
- [ ] `tst_audit_fixes.py` passes (31 tests)
- [ ] Manual test: one cycle with current PHENIX client against new server

### Deployment

1. Merge changes to server code
2. Run full test suite including compatibility tests
3. Deploy to staging server
4. Run current PHENIX release against staging (smoke test)
5. Deploy to production
6. Monitor for errors from clients (log `client_protocol_version`)

### How to Add a New session_info Field

1. Add the entry to `agent/contract.py :: SESSION_INFO_FIELDS` with the next version number and a safe default
2. In `programs/ai_agent.py`, add the field to the `session_info` dict in `_query_agent()`
3. On the server, access it **only** via `session_info.get("field", default)` where the default matches the contract
4. Run `tst_contract_compliance.py` to confirm compliance
5. Create a new frozen fixture in `tests/fixtures/client_vN/` if the change represents a meaningful new scenario

### Files

| File | Role |
|------|------|
| `agent/contract.py` | Contract registry, version constants, runtime helpers |
| `agent/graph_nodes.py` | Pipeline nodes with normalization, version gate, stop guards, response guarantees |
| `programs/ai_agent.py` | Client: sends `client_protocol_version`, handles `warnings` |
| `tests/tst_contract_compliance.py` | 10 static analysis + contract tests |
| `tests/tst_old_client_compat.py` | 6 frozen fixture tests |
| `tests/tst_shared_code_imports.py` | 4 import safety tests |
| `tests/tst_utils.py` | Test infrastructure (assert helpers, test runner) |
| `tests/fixtures/client_v1/` | 2 frozen old-client payloads |
| `tests/fixtures/client_v3/` | 1 frozen current-client payload |

---

## 9. Event System

### Status: IMPLEMENTED ✓

The event system provides structured, transparent logging of the agent's decision-making process. Events flow from graph nodes through the API to the display layer, ensuring consistent output regardless of local or remote execution.

**Current version:** v40 (January 2025)

### Implementation Files

| File | Purpose |
|------|---------|
| `agent/event_log.py` | EventType constants, Verbosity levels, EventLog class |
| `agent/event_formatter.py` | EventFormatter for human-readable output |
| `agent/graph_nodes.py` | Event emission via `_emit()` helper |
| `phenix_ai/run_ai_agent.py` | Events in API response |
| `programs/ai_agent.py` | EventFormatter integration, verbosity parameter |

### Goals

1. **Transparency** - Show what decisions are made and why at each step
2. **Consistency** - Same output format for local and remote execution
3. **Controllability** - Support quiet, normal, and verbose verbosity levels
4. **Full Reasoning** - Show complete LLM reasoning without truncation
5. **User Feedback** - Prominently warn when user requests can't be fulfilled

---

### Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          Decision Making                                 │
│  ┌─────────────┐   ┌─────────────┐   ┌─────────────┐   ┌─────────────┐  │
│  │  perceive   │ → │    plan     │ → │    build    │ → │  validate   │  │
│  └──────┬──────┘   └──────┬──────┘   └──────┬──────┘   └──────┬──────┘  │
│         │                 │                 │                 │          │
│         ▼                 ▼                 ▼                 ▼          │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │                   state["events"] (list of dicts)               │    │
│  └─────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     run_ai_agent.py                                      │
│  response = create_response(..., events=state["events"])                │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                    ┌───────────────┴───────────────┐
                    ▼                               ▼
            ┌─────────────┐                 ┌─────────────┐
            │ LocalAgent  │                 │ RemoteAgent │
            └──────┬──────┘                 └──────┬──────┘
                   └───────────────┬───────────────┘
                                   ▼
            ┌─────────────────────────────────────────────┐
            │              ai_agent.py                     │
            │  EventFormatter(verbosity).format_cycle()   │
            │  print(output, file=self.logger)            │
            └─────────────────────────────────────────────┘
```

---

### Event Types

| Event Type | Verbosity | Description | Key Fields |
|------------|-----------|-------------|------------|
| `cycle_start` | quiet | New cycle beginning | `cycle_number` |
| `cycle_complete` | quiet | Cycle finished | `program`, `success` |
| `state_detected` | normal | Workflow state determined | `workflow_state`, `experiment_type`, `reason`, `valid_programs` |
| `metrics_extracted` | normal | Metrics parsed from log | `r_free`, `r_work`, `resolution`, `map_cc`, `*_prev` |
| `metrics_trend` | normal | Trend analysis result | `improving`, `plateau`, `cycles_analyzed` |
| `sanity_check` | normal | Sanity check result | `passed`, `red_flags`, `warnings` |
| `program_selected` | normal | Program decision | `program`, `reasoning` (full), `source`, `provider` |
| `program_modified` | normal | LLM choice overridden | `original`, `selected`, `reason` |
| `user_request_invalid` | **quiet** | User requested unavailable program | `requested_program`, `reason`, `selected_instead`, `valid_programs`, `suggestion` |
| `files_selected` | verbose | Input file selection | `selections` dict with per-input details |
| `command_built` | normal | Final command | `command`, `program` |
| `stop_decision` | normal | Whether to continue | `stop`, `reason` |
| `directive_applied` | normal | User directive triggered | `directive`, `action` |
| `error` | quiet | Error occurred | `message`, `details` |
| `warning` | normal | Warning issued | `message` |
| `thought` | verbose | LLM chain-of-thought/reasoning traces | `message` |
| `file_scored` | verbose | Individual file scoring | `file`, `score`, `reason` |
| `debug` | verbose | Debug trace | `message` |

---

### Verbosity Levels

```python
class Verbosity:
    QUIET = "quiet"    # Errors, warnings, and cycle summaries
    NORMAL = "normal"  # Key decisions and metrics (default)
    VERBOSE = "verbose"  # Full detail including file selection, internal state
```

Note: `debug` is an event type (not a verbosity level). Events of type `debug` are shown at `verbose` verbosity.

**PHIL Parameter:**
```phil
verbosity = normal
  .type = choice(quiet, normal, verbose)
  .help = Controls how much detail is shown in agent output.
```

---

### Example Output

### Normal Verbosity

```
================================================================================
 CYCLE 3
================================================================================

State: xray_refined
  Experiment type: xray
  Reason: Has refined model and reflection data

Metrics:
  R-free: 0.2600 → 0.2400 (↓ improved)
  R-work: 0.2200 → 0.2100 (↓ improved)
  Resolution: 2.50 Å

Trend Analysis:
  Improving (no plateau detected)
  Cycles analyzed: 3

Sanity Check: PASSED

Decision: phenix.refine
  Source: LLM (google)
  Reasoning: R-free continues to improve (0.26 → 0.24). No plateau detected
  after 3 cycles. Continue refinement to optimize geometry and reduce R-factors.
  The model shows good density fit with clashscore of 5.2.

Command:
  phenix.refine refine_002_001.pdb data.mtz output.prefix=refine_003

--------------------------------------------------------------------------------
```

### User Request Invalid Warning

When a user requests an unavailable program, this warning is **always shown** (even at quiet verbosity):

```
============================================================
  WARNING: Requested program not available
============================================================
  You requested: phenix.refine
  Reason: Not valid in current workflow state 'xray_initial'
  Running instead: phenix.xtriage
  Available programs: phenix.xtriage
  Suggestion: This program requires different conditions. Need to analyze data first.
============================================================
```

### Verbose (adds file selection)

```
...
File Selection:
  Model: refine_002_001.pdb
    Reason: best_files
  Mtz: data.mtz
    Reason: rfree_locked
  Sequence: sequence.fa
    Reason: auto_selected
...
```

### Quiet

```
CYCLE 3: phenix.refine
```

---

### Implementation Files

| File | Purpose |
|------|---------|
| `agent/event_log.py` | EventType constants, Verbosity levels, EventLog class |
| `agent/event_formatter.py` | EventFormatter class for human-readable output |
| `agent/graph_nodes.py` | Event emission at decision points |
| `agent/command_builder.py` | File selection tracking with reasons |
| `knowledge/api_schema.py` | `events` field in response schema |
| `phenix_ai/run_ai_agent.py` | Events included in response building |
| `programs/ai_agent.py` | Verbosity parameter, display integration |

---

### Usage

### Emitting Events in Graph Nodes

```python
from agent.event_log import EventType

def perceive(state):
    # Initialize events if not present
    if "events" not in state:
        state = {**state, "events": []}
    
    # Use _emit helper
    state = _emit(state, EventType.STATE_DETECTED,
        workflow_state="xray_refining",
        experiment_type="xray",
        reason="Has refined model")
    
    return {**state, ...}
```

### Formatting Events for Display

```python
from agent.event_formatter import EventFormatter
from agent.event_log import Verbosity

formatter = EventFormatter(verbosity=Verbosity.NORMAL)
output = formatter.format_cycle(events, cycle_number=3)
print(output, file=self.logger)
```

---

### File Selection Tracking

The CommandBuilder records WHY each file was selected:

| Reason | Description |
|--------|-------------|
| `best_files` | From BestFilesTracker (highest scored) |
| `rfree_locked` | Locked R-free MTZ from first refinement |
| `llm_selected` | LLM explicitly chose this file |
| `best_files_override` | LLM choice overridden by best_files |
| `auto_selected` | Fallback based on category/extension |

---

### Testing

```bash
# Run event system tests
python3 tests/tst_event_system.py

# Test coverage includes:
# - EventType and Verbosity constants
# - EventLog emit and filtering
# - EventFormatter at all verbosity levels
# - USER_REQUEST_INVALID formatting
# - Full reasoning preservation (no truncation)
```

---

### Change History

| Version | Date | Changes |
|---------|------|---------|
| v40 | 2025-01 | Added USER_REQUEST_INVALID event type |
| v38 | 2025-01 | Phase 4: Display integration, verbosity parameter |
| v37 | 2025-01 | Phase 3: Transport integration (API schema, response building) |
| v36 | 2025-01 | Phase 2: Instrumented graph_nodes.py and command_builder.py |
| v34 | 2025-01 | Phase 1: Created EventLog, EventFormatter, 17 event types |

---

## 10. LLM-Assisted Development

This section describes the workflow that PHENIX developers
follow when using LLMs (Claude, Gemini, ChatGPT, etc.)
to write or modify code. The goal is not to restrict LLM
use — it is to ensure that LLM-generated code meets the
same standards as hand-written code: correct, tested,
portable, and maintainable.

Two companion documents contain the coding standards that
the LLM should follow:

- **`CCTBX_LLM_PROGRAMMING_GUIDELINES.md`** — general
  cctbx/PHENIX coding standards, style rules, pitfalls,
  and checklist. Applies to all PHENIX code.
- **`AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md`** — additional
  patterns specific to the AI Agent codebase: parameter
  verification against programs.yaml, logging conventions,
  import fallbacks, analysis mode routing, the three error
  classification systems, client-server code path
  awareness, and session state persistence.

---

### 1. Session Setup

Before asking the LLM to write any code, provide it with
context. The quality of LLM output is directly proportional
to the quality of context you provide.

**Required context (always attach):**

- **`CCTBX_LLM_PROGRAMMING_GUIDELINES.md`** — attaching
  this at the start of every session is the single
  highest-leverage step you can take. It grounds the LLM
  in the coding standards, patterns, and pitfalls it needs
  to follow.
- **`AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md`** — attach
  this as well when working on AI Agent code (`agent/`,
  `knowledge/`, `programs/ai_agent.py`). It covers
  agent-specific patterns not in the general guide.
- **The existing source file(s)** to be modified, in full.
  LLMs produce better code when they can see the
  surrounding style, imports, and conventions.

**Recommended context (attach when relevant):**

- `OVERVIEW.md` or `ARCHITECTURE.md` — grounds the LLM
  in the system design so it doesn't reinvent existing
  infrastructure.
- The test file(s) that cover the code being modified —
  helps the LLM understand what's already tested and write
  compatible new tests.
- Error messages or log output that motivated the change —
  gives the LLM concrete evidence of the problem.
- `programs.yaml` — essential when the work involves
  command building, strategy flags, or file selection.
- Run logs from tutorial runs — when debugging agent
  behavior, real log output is far more useful than
  descriptions of the problem.

**Never rely on the LLM's training data for PHENIX
specifics.** LLMs hallucinate PHIL parameters, API
signatures, and command-line flags that don't exist. Always
provide the actual source code or YAML definitions as
context rather than asking the LLM to recall them.

---

### 2. Plan-Verify-Execute Workflow

All LLM-assisted development follows a three-phase cycle.
Do not skip phases, even for "simple" changes — the bugs
that escape review are always the ones that looked simple.

### Phase 1: Plan

Ask the LLM to produce a markdown plan before writing any
code. The plan must include:

- **Problem statement**: What is broken or missing, with
  specific evidence (error messages, failing test output,
  user report).
- **Approach**: How the fix or feature will work, at a
  level of detail sufficient for another developer to
  review. For scientific code, this includes the
  mathematical formulas or algorithms being implemented.
- **Dependencies and edge cases**: What existing code is
  affected, what inputs could be unusual (None values,
  empty files, Windows paths, resumed sessions), and how
  each is handled.
- **Implementation steps**: A numbered sequence of small,
  testable increments.
- **Test plan**: Which new tests will be written, which
  existing tests might be affected, and how to verify
  correctness.

**Review the plan before proceeding.** Read it yourself.
Optionally, show it to a second LLM and ask: *"What could
go wrong with this approach?"* This cross-check catches
design-level issues before any code is written.

### Phase 2: Implement and verify

The LLM implements the plan one step at a time. After each
step:

1. **The LLM self-reviews** (this is built into the
   guidelines document — the LLM is instructed to check
   for edge cases, scoping issues, and style before
   presenting each step).
2. **You ask**: *"Are there any other fixes or additions
   you would like to make to this change before we move
   on?"* This prompt is the single most effective quality
   gate. It triggers the LLM to re-examine its own work
   and consistently catches bugs that the initial
   implementation misses.
3. **The LLM runs tests** and confirms they pass before
   presenting the step.
4. **You review** the code. Check that it matches the plan,
   follows the coding standards, and makes sense
   scientifically.

Do not batch multiple steps and review at the end. The
cost of finding a bug in step 2 while you're reviewing
step 5 is much higher than catching it immediately.

### Phase 3: Final audit

After all steps are complete, ask the LLM to run a full
audit:

- Parse-check every changed file
- Run all relevant test suites
- Check for line-length violations in changed regions
- Verify the checklists in both guidelines documents

Then ask: *"Check every changed file for the patterns in
the Common Pitfalls section."* This catches the classes of
bugs that LLMs produce most frequently (None-safety,
path separators, floating point comparison, parameter
hallucination).

---

### 3. What the Programmer Is Responsible For

The LLM is a tool. The programmer is responsible for the
code that gets committed.

### Scientific correctness

LLM-generated code is a hypothesis, not a proof. The
programmer must verify:

- **Algorithms match the literature.** If the code
  implements a crystallographic calculation (structure
  factor computation, map coefficient generation, peak
  searching), check the math against the relevant
  reference, not just the LLM's assertion that it's
  correct.
- **Units are consistent.** Distances in Ångströms, angles
  in degrees, R-factors as fractions (0.25, not 25),
  B-factors in Å². An LLM that has seen both conventions
  in its training data may use either one.
- **Coordinate-space operations are cache-aware.** Spatial
  data structures (KD-trees, neighbor lists) must be
  invalidated when the unit cell or symmetry changes.
  LLMs build correct trees but frequently omit the
  invalidation logic.
- **Numerical stability.** Division by zero, log of zero,
  very small denominators in correlation coefficients.
  LLMs rarely add these guards unprompted.

### Platform correctness

PHENIX runs on macOS, Linux, and Windows. The programmer
must verify:

- No hardcoded path separators (`/` or `\`) — only
  `os.path.join()`.
- Explicit `encoding='utf-8'` on file operations where
  the content may contain non-ASCII characters.
- No reliance on Unix-specific behavior (signal handling,
  symlinks, case-sensitive filesystems).

### Test adequacy

The programmer must verify that the LLM's tests actually
test what they claim to:

- **Coverage of edge cases.** Does the test exercise None
  inputs, empty lists, missing keys, and boundary values?
  LLMs tend to write tests for the happy path.
- **Meaningful assertions.** `assert result is not None`
  doesn't verify correctness — it verifies existence.
  Tests should check specific values, ideally against
  known crystallographic results.
- **Floating point comparisons.** Tests must use
  `approx_equal()` (from `libtbx.test_utils`), never
  `==`, for any value derived from a calculation.
- **Sandbox stubs lie if not pinned to real-function
  semantics.** When a test imports a function from
  `agent.X` and the sandbox provides a stub for `agent.X`,
  the test exercises stub semantics, not real semantics.
  An LLM may write a stub that "looks right" but does
  something subtly different from the real function
  (substring match vs word-boundary regex match, eager
  vs lazy evaluation, normalization steps that change
  the answer).  See "Semantic-pin tests" below for the
  policy that prevents this.

#### Semantic-pin tests (the v119.H10 → H11 lesson)

The H10 → H11 ship cycle in May 2026 exposed a sandbox
failure mode worth capturing as policy.

**What happened**: H10 added an `exclude_patterns` filter
at 8 file-selection paths in `command_builder.py` to
reject model mmCIF files (`refine_*.cif`) from being
injected positionally as a second model into
`phenix.refine`.  The sandbox stub for
`agent.file_utils.matches_exclude_pattern` did
substring matching.  H10's K-suite passed (15+ tests).
H10 was packaged as "fixes AIAgent_62 cycle-7 crash."

When Tom installed H10 and ran the test suite against
the real PHENIX environment, Test 4 failed.  A diagnostic
script that called `matches_exclude_pattern` directly
revealed the real function uses **word-boundary regex
matching**, not substring matching:

```python
# agent/file_utils.py — real semantics
re.search(
    r'(?:^|[_\-\.])' + re.escape(pat_stem) + r'(?=[_\-\.]|$)',
    stem)
```

The YAML pattern `"refine_"` (with trailing underscore)
matches `"refine_001_001.cif"` under substring semantics
but NOT under word-boundary semantics (the regex requires
another boundary char after the trailing `_`, but `0` is
not a boundary).  H10's filter fired but rejected
nothing.

**The policy**: every function with non-obvious semantics
(regex matching, parsing, escaping, normalization,
canonicalization) that's stubbed for sandbox testing
must have a paired **semantic-pin test** in the test
suite.  The semantic-pin test:

1. Calls the REAL function (not the stub), via
   try/except import to handle sandbox-without-real
   environments.
2. Asserts behavior on a fixture of named cases
   including at least one case where the stub-vs-real
   semantics would diverge.
3. Has a graceful sandbox-skip fallback when the
   real function isn't importable, so lightweight CI
   environments don't error.

**Template**: `tests/tst_autosol_bugs.py::test_bug11_matches_exclude_pattern_semantics`.
10 named cases pinning the real function's behavior,
including the trailing-underscore case that distinguishes
substring from word-boundary semantics.  Read the test
and copy its shape when introducing any similar
sandbox-stubbed function.

**Why this works as a parity check**: when a sandbox
test imports `matches_exclude_pattern` and the import
resolution finds a stub, the same test that pins the
real function's behavior will FAIL against the stub IF
the stub semantics diverge from real semantics.  The
divergence is caught at the moment a sandbox test runs.
No separate "differential parity test" is needed unless
you also want to enforce stub-real equivalence across
fixtures that don't already appear in the semantic pin.

**Forward-looking corollary**: if you write a stub, write
it to match the real function's semantics, not to a
plausibly-correct intuition.  When in doubt, read the
real function's source.  When the real function lives
in a module you can't import in sandbox, document the
stub's intended semantics in a docstring and add a
semantic-pin test that exercises both the stub and the
real function via try/import.

**H12 extension — broaden the pinning beyond the one bug**:
the H10/H11 incident pinned `matches_exclude_pattern`
specifically because that was the function with documented
divergence.  But the same risk applies to every function in
`agent/file_utils.py` — small, non-trivial logic that's easy
to stub badly.  H12 generalized the pattern by adding
`tests/tst_file_categorizer.py` with 8 tests / ~75
documented assertions covering EVERY public function in
`file_utils.py`: `classify_mtz_type` (rule ordering invariant),
`get_mtz_stage`, `get_category_for_extension` (multi-dot,
no-extension edge cases), `is_mtz_file`, `is_model_file`
(load-bearing `.pdb`-overrides-`'ligand'`-filter quirk),
`is_map_file`, `is_sequence_file`, plus cross-function
consistency.

**Pattern**: when one function in a module exhibits the
sandbox-stub divergence problem, the other functions in the
same module are candidates for the same defense.  The cost
of a comprehensive semantic-pin suite is low (one new test
file, no production code change); the benefit is pre-emptive
detection of future divergences.

**The highest-value pinning targets are non-obvious quirks**.
A test that pins the obvious behavior is less useful than
one that pins a documented exception.  When you find
yourself writing "BTW, this case is the exception to the
rule" in a docstring, that case should have a test.
`is_model_file('ligand_fit.pdb') == True` (because `.pdb`
extension overrides the `'ligand'`-in-name filter) is the
template — a load-bearing quirk that's easy to "simplify"
away in a refactor.

#### Centralizing access patterns, not just data tables (v119.H1 → H13 lesson)

The H13 ship (May 2026) surfaced a different category of
inconsistency: H1 centralized the LLM `DEFAULT_MODELS` table
in `core/llm.py`, but each consumer of the table acquired
its own ad-hoc access pattern over time.  Three of the four
ollama consumer sites read the table directly via
`default_model_for_provider("ollama")`.  The fourth site
(in `get_llm_and_embeddings`) honored `OLLAMA_LLM_MODEL`
env-var as an override.  Tom set `OLLAMA_LLM_MODEL=qwen2.5:72b`
expecting the documented override behavior; two of the four
consumers silently ignored it.

**The policy**: centralization happens in stages.  Stage 1
extracts the data table (H1 did this).  Stage 2 extracts the
access patterns or rules over the table.  H13's
`resolve_model_for_provider()` is the Stage 2 helper: it
wraps the precedence rule "env-var wins, central default
falls back" and is consumed by every site that needs a
model name.

When introducing centralization:

1. **Identify both the data AND the rules**.  Don't ship the
   table without auditing the consumers' access patterns —
   if they're inconsistent, the table's central status is
   illusory.
2. **Extract the rule before the rule diverges further**.  A
   new consumer added after Stage 1 but before Stage 2 will
   add its own access pattern, increasing the eventual
   migration cost.
3. **Empty-string handling for env-var overrides**: treat
   empty string as unset, matching shell intuition.
   `setenv OLLAMA_LLM_MODEL ""` should fall through to the
   default, not silently break the consumer.  Pinned in
   `tst_default_models.py::test_resolve_model_for_provider_empty_env_falls_through`.
4. **Source-scan tests need to accept multiple helper names**
   during a transition.  `tst_default_models.py`'s per-site
   tests now accept BOTH `default_model_for_provider` AND
   `resolve_model_for_provider` references — the invariant is
   "every consumer uses a central helper", not "every consumer
   uses the same helper name."

#### Error categories should map to operator actions (v119.H13 lesson)

H13's `_classify_provider_error` splits 404 errors into
three sub-categories: `MODEL_RETIRED` (provider deprecated
the model), `MODEL_UNAVAILABLE` (model unknown to this
server), and `FAILED` (generic).  Each has a distinct
operator-action hint.  Pre-H13, all 404s were classified
as RETIRED with hint "update DEFAULT_MODELS" — wrong action
for Tom's actual case (`ollama pull X` was the right action).

**The policy**: when adding a diagnostic marker for an error
class, ask "what should the operator do when they see this?"
If multiple actions are possible depending on context, the
marker is too coarse — split it.

Design constraints when splitting:

1. **Co-occurrence rules** to avoid false-positives.  H13's
   RETIRED classifier requires "model" word within 80 chars
   of retirement phrase, distinguishing legitimate model
   retirement from edge-proxy 404 pages mentioning
   "deprecated endpoint."  A naive single-phrase match
   would mis-route the latter.
2. **Class-attribute interrogation** for SDK robustness.
   Read `exc.status_code` and `exc.response.status_code`
   programmatically BEFORE falling back to string matching
   on `str(exc)`.  Some SDK exceptions hide details behind
   opaque `__str__` representations.
3. **Load-bearing test sentinels**.  Pin the load-bearing
   distinctions with explicit test cases — e.g.,
   `test_classify_edge_proxy_404_not_retired`.  If a future
   refactor weakens the co-occurrence rule, the sentinel
   fires loudly.
4. **Backward-compatible emission**.  When adding richer
   markers, ALSO emit the legacy marker (in addition, not
   in place of).  Log-grep pipelines and monitoring
   dashboards depend on the legacy markers; H13 emits both
   `[DIRECTIVE_EXTRACTION_FAILED]` (legacy) and the richer
   tag (new) so existing infrastructure continues to work
   while operators get the actionable hint.

#### Goal phrases vs method requests in keyword lists (v119.H14 lesson)

H14 fixed a long-standing bug in `_ACTION_TABLE["solve"]`: the
keyword list conflated GOAL phrases (`"solve the structure"`)
with METHOD requests (`"molecular replacement"`, `"phaser"`,
`"mr "`).  Adding the goal phrase made the extractor inject
`start_with_program = phenix.phaser` on any README that
expressed the goal — even on SAD/MAD datasets where phaser is
the wrong method.

**Rule for keyword-list authoring**: each entry in
`_ACTION_TABLE` represents a METHOD request — "the user
explicitly asked for THIS tool."  Goal phrases ("solve",
"refine to convergence", "validate the model") belong in
intent classification (`classify_intent` in
`agent/intent_classifier.py`), not in method-action
keyword lists.  The two layers serve different purposes:

| Layer | Inputs | Outputs | Use case |
|---|---|---|---|
| `_ACTION_TABLE` | exact keyword match | program name + workflow directive | "user asked for phaser → run phaser" |
| `classify_intent` | full advice text | intent classification (`solve`/`refine`/`validate`/etc.) | "user wants to solve → adjust workflow rules accordingly" |

**Check before adding a keyword to `_ACTION_TABLE`**:

1. Is this keyword UNAMBIGUOUSLY a method name?  If a user
   says it, do they specifically mean "use this tool"?
   ("phaser" → yes.  "solve" → no, the method depends on
   data type.)
2. Could this keyword combine with another action's keyword
   to trigger the `n > 1` branch in
   `_apply_workflow_intent_fallback`?  If yes, walk through
   the resulting `start_with_program` + `after_program` and
   verify the workflow is correct.
3. Pin the decision with a K-test: assert the keyword
   triggers the action AND assert that goal-only phrasings
   do NOT trigger it.

#### Single-emit invariant for output channels (v119.H14 lesson)

H14 removed a duplicate `diagnostic_messages` relay that
survived a v119.H5 → H5.1 refactor.  Both the dispatcher and
the client iterated the list and wrote to stderr, causing
every marker to be emitted twice.  The bug was invisible
functionally — operators just saw duplicate lines — and
survived until run_39 batch analysis surfaced the
60.6%-of-runs pattern.

**Rule for adding output channels**: any structured output
channel (markers, telemetry, log lines) must have ONE
canonical relay site, documented as such.  When refactoring
to introduce a new central relay, audit ALL pre-existing
relay sites and remove them in the SAME ship — not later,
not "as cleanup."  Surviving duplicate relay sites cause
silent over-emission that's easy to miss and hard to detect
without batch-level statistical analysis.

**Check before adding a marker emit**:

1. Is there already a relay function for this channel
   (`_relay_diagnostic_messages_to_stderr`, etc.)?  If so,
   USE it, don't write to stderr inline.
2. If you're adding a relay function, GREP for any existing
   inline stderr writes that emit markers from the same
   channel.  Delete them in the same change set.
3. Document the relay function's role in its docstring:
   "this is the SINGLE site that writes marker X to
   stderr."  If a later refactor invalidates the claim,
   the docstring becomes the audit trail that the next
   investigator follows.

H14 also enforces this with a K-test that source-scans
`programs/ai_agent.py` for `for _msg in _diagnostics` — if
the duplicate ever returns, the test catches it.

#### Positive shape checks complement sentinel sets (v119.H14 lesson)

H14 strengthened `validate_directives` for `space_group`.
Pre-H14 the validator had a sentinel set
(`_SYMMETRY_SENTINELS`) and negative structural checks
(starts-with-letter, length <= 25).  Those caught most cases
but missed prose phrases that happened to look like real
values (`"Solve the structure"` passes all three checks).
H14 added a positive Hermann-Mauguin shape check.

**Rule for validating LLM-emitted parameter values**: use a
LAYERED check.

| Layer | What it catches | When to use |
|---|---|---|
| Sentinel set (frozenset) | Known placeholder phrases | LLMs frequently emit these; cheap match |
| Negative structural checks | Obvious shape violations (empty, too long, starts with digit) | First-pass quick rejection |
| Positive shape check (regex) | Anything that ISN'T the expected form | Final guarantee against prose / hallucination |

The positive check is the SOUNDER guarantee but is the most
costly to author (need to know the value's grammar).  Use it
when the value space is well-defined (Hermann-Mauguin
symbols, ISO dates, semver numbers) and the cost of
passing a bad value through is high (PHENIX crashes with
obscure errors, data corruption, etc.).

**Anti-pattern to avoid**: relying ONLY on a sentinel set
when the value space has structural rules.  Sentinel sets
are necessarily incomplete because LLMs can emit arbitrarily
new placeholder phrasings.  A positive shape check closes
that open-ended gap.

#### Validate positive shape checks against the FULL domain (v119.H14 lesson)

H14's positive HM shape check went through two review passes
before reaching its final form.  Both reviews caught real bugs
that would have shipped silent false-negatives — valid user
inputs being dropped as malformed.

**Lesson 1 — verify against an authoritative list, not your
hand-picked sample.**  The initial H14 regex was tested against
a hand-picked set of ~15 well-known protein space groups and
all passed.  But running it against the full 230 International
Tables symbols showed only 71 matched (31%).  The missing 159
included all monoclinic slash forms (`P21/c`, `P21/n`), all
orthorhombic mirror/glide groups (`Pmma`, `Pnma`, `Pbca`), all
tetragonal mirror forms, and all cubic high-symmetry groups
(`Pm-3m`, `Im-3m`, `Fd-3m`).  Hand-picked test sets give false
confidence; the full domain enumeration is the right
benchmark.

Concrete K-test pattern when the domain is enumerable:

```python
def test_pattern_accepts_full_domain():
    """Pin the regex behavior against the AUTHORITATIVE list,
    not a hand-picked subset."""
    all_canonical = [...]  # complete enumeration, ideally from
                            # a library (cctbx.sgtbx, etc.)
    assert len(all_canonical) == EXPECTED_TOTAL, "..."
    failures = [v for v in all_canonical if not pattern.match(v)]
    assert not failures, "Domain values rejected: %r" % failures
```

**Lesson 2 — consider alternative forms, not just canonical
ones.**  After the regex was fixed to match all 230 canonical
symbols, Gemini's external review (Risk A) caught that
alternative cell/origin settings were still rejected:

- `R3:H` vs `R3:R` (rhombohedral axis choice)
- `P4/n:1` vs `P4/n:2` (origin choice)
- `P21/c:b` (unique-axis cell choice)

These aren't in the 230-symbol canonical list — they're
*derived* forms that appear in real PDB/mmCIF metadata.  The
fix required adding `:`, `h`, `r` to the alphabet.  When
validating a domain with multiple notation conventions,
enumerate the conventions explicitly:

| Question | For the HM regex |
|---|---|
| What's the canonical form? | 230 short symbols (ITA Vol A) |
| What alternative forms exist? | Colon-suffixed settings (`:H`/`:R`, `:1`/`:2`, `:a`/`:b`/`:c`) |
| What separators / whitespace are accepted? | Space-separated full forms, slash, dash, underscore, colon |
| What case variants are accepted? | Lower and upper (PHENIX is case-insensitive) |

**Lesson 3 — document acceptable residual leaks instead of
chasing perfection.**  Even the final regex admits short English
words that happen to use only HM-alphabet characters (`Panda`,
`Bed`, `Cab` — Gemini Risk B).  The right response is NOT to
keep tightening the regex.  No realistic LLM emits "Panda" as a
space-group value, and closing this gap fully would require
enumerating the 230 + alternative-setting symbols via
`cctbx.sgtbx` — a substantive dependency.  Instead, pin the
limitation with a K-test that documents the leak and the
acceptable upstream filters:

```python
def test_pattern_known_limitation():
    """DOCUMENTED LIMITATION: ...
    The acceptable upstream filters are: ...
    The future path forward (if production surfaces this class
    of hallucination) is to enumerate the domain via <library>
    rather than tighten the regex."""
    # Pin the current behavior, not assert it goes away
```

This makes the limitation discoverable (future contributors
won't think it's a bug) and points to the right next step
(enumeration via library, not regex iteration).

#### Validate on every path, not just the happy path (v119.H14.1 lesson)

H14 added a positive Hermann-Mauguin shape check on
`space_group` inside `validate_directives`.  Local sandbox
tests proved the check worked.  Tom's 2026-05-26 production
verification on ollama showed the bug was NOT closed in
production — the bogus `Not explicitly mentio` value still
reached the displayed directives even though
`scanner_version=119.H14` confirmed H14 was installed.

**Root cause**: `directive_extractor.py` had TWO paths that
produce a directives dict:

1. LLM-success path → `extract_directives` calls
   `validate_directives` on the result
2. LLM-failure / `use_rules_only` path → falls back to
   `extract_directives_simple` and returns the result
   DIRECTLY (no validation)

The simple extractor had its own permissive space_group regex
that captured `"Not explicitly mentio"` from the
preprocessed advice (the truncation was the regex's
`{1,20}` quantifier, not an LLM length cap — the
H14 CHANGELOG misattributed this).  The dict went straight
to the agent.  H14's sentinel check was bypassed.

**Lesson**: when a validator is added to fix a class of bugs,
trace EVERY code path that produces the value being validated
— not just the path you originally identified.  Two
diagnostic moves to make this concrete:

1. **Grep for all producers**.  For `directive_extractor.py`,
   the validator was `validate_directives` and the producers
   were `extract_directives` (with validation) and
   `extract_directives_simple` (without).  A single grep for
   `def extract_` would have surfaced both producers and
   prompted the question "which of these call
   `validate_directives`?".

2. **End-to-end production verification, not just sandbox
   K-tests**.  The H14 K-tests called `validate_directives`
   directly and passed.  Tom's production run was what
   exposed the gap.  K-tests that reproduce the production
   entry point (`extract_directives(advice, provider="ollama")`)
   with empty LLM responses would have caught the bypass.

**Architectural fix pattern: converged validators**.  Rather
than duplicate the validator inline at each producer site,
make the validator the canonical final-sanity step that ALL
producers converge through:

```python
def extract_thing_simple(advice, log=None):
    # ... pattern-matching extraction ...
    return validate_thing(result, log)   # ← final-sanity

def extract_thing_llm(advice, ...):
    # ... LLM extraction ...
    return validate_thing(result, log)   # ← same validator
```

Now: adding a new extractor adds one line (call the
validator before return); adding a new validator rule changes
one site (the validator function); neither requires touching
the other.

**Prerequisite check before applying converged-validator
pattern**: audit whether `validate_thing` strips anything
that the producers legitimately set.  H14.1 surfaced a latent
bug during this audit — `start_with_program` was set by the
simple extractor's `_resolve_after_program` overlay but
wasn't in `VALID_STOP_CONDITIONS`.  Pre-H14.1, the bug was
invisible because the LLM path called validate BEFORE the
overlay added the key.  Adding the converged-validator
pattern would have stripped the key, silently breaking H14
Item 1.  The audit caught this in time:

```bash
# For each key the simple extractor sets:
$ grep -n 'directives\["X"\]' agent/directive_extractor.py

# Check whether that key is in the validator's allow-list:
$ grep -A 10 'VALID_X_CONDITIONS' agent/directive_extractor.py
```

Any key in the producer but not in the allow-list is a
latent bug surfaced by the converged-validator pattern.  Fix
it BEFORE applying the converged pattern, not after.

#### Targeted reset blast radius (v119.H15 lesson)

H15 fixed Tom's bromodomain resume failure (run 144).  The
original H15 plan proposed reopening every stage downstream
of any matched program — "if the user mentions phaser, reset
phaser AND every refine/build stage after it."  Gemini's
critique flagged this as O(N) blast radius: on a 12-stage
cryo-EM pipeline, mentioning an early program in resume
advice would reset 10+ stages of completed work.

**Lesson**: when a re-trigger mechanism needs to invalidate
prior state, prefer the SMALLEST blast radius consistent with
the user's stated intent.  Three principles:

1. **Single-stage default**.  Pick ONE stage to reset based
   on a deterministic rule.  H15's rule: "the LATEST
   completed stage whose `programs` list contains the named
   program."  This matches the strongest semantic
   interpretation: "the user wants to re-do this program
   in its most recent context."  Reopening earlier stages
   would discard work done after them.

2. **No cascade**.  Don't reset stages downstream of the
   reopened one — they'll re-fire naturally when the
   reopened stage completes, IF they need to.  Cascade
   resets assume downstream stages depend on the reopened
   one, which is often but not always true.

3. **Respect skip semantics**.  Skipped stages stay
   skipped.  The user can't un-skip a stage by mentioning
   its program; that requires explicit directives.

The result is O(1) blast radius regardless of plan size.
The 7 K-tests in `tests/tst_resume_reopen_stages.py` pin
this contract — specifically test §G ("stage's strategy
already honors directive") proves that the design reopens
the stage anyway, treating user re-assertion as the
authoritative signal.

**Generalize**: any "user re-asserts X" mechanism should
identify the SMALLEST unit of state affected by X and reset
only that.  Cascade resets are an anti-pattern for resume
flows — they pretend to be safe but discard legitimate
prior work.

#### Reactive recovery requires executor support (v119.H17 → H17.1 lesson)

H17 added a new `strip_parameter` resolution kind to the
recoverable-errors system: when the LLM passes an
inappropriate flag (e.g., `map_file=` to autobuild without
phases), strip it from the retry command rather than
trying to add a fix-up flag.  The implementation looked
complete: YAML entry declared the pattern and the
`strip_parameters` list; `error_analyzer.py` gained a
`strip_flags` field on `ErrorRecovery`; the pre-deploy
review concluded "ready to ship."

In production, H17 fired correctly — the `[NOTICE]
DETECTED RECOVERABLE ERROR` log appeared with "Action:
Stripping [...]" — but the retry command still contained
the offending flag.  H17 was analyzer-side only.  The
executor in `programs/ai_agent.py` had no code to actually
strip the flag from the command before running it.

**Lesson**: when adding a new resolution kind to a
detect-and-recover system, trace the END-TO-END path
through the executor.  Three concrete checks:

1. **Detection emits a SIGNAL the executor can act on.**
   `ErrorRecovery.strip_flags` is the signal; the
   executor must read it.
2. **The executor consumes the signal.**  In ai_agent.py,
   `_handle_recovery` was the obvious site for stashing
   the strip directive; `_execute_command` was the obvious
   site for applying it.  Pre-H17.1 neither was wired.
3. **The on-disk command differs after the strip.**
   K-tests should assert the command STRING differs
   between pre-strip and post-strip — not just that the
   recovery struct contains the strip directive.

The H17 → H17.1 pattern: ship the detection side first
(YAML + analyzer) because that exercises the detection
logic against real production errors.  Once detection is
confirmed firing, the executor wire-up is mechanical.
But don't conflate "detection works" with "fix works" —
the K-tests for H17.1 (`tst_h17_strip_executor.py`)
specifically test the executor regex against PHIL spacing,
quoted-with-spaces, single quotes, and end-of-line — the
forms a naive `\S+` pattern would corrupt.

**Generalize**: any add-a-new-resolution-kind change to a
reactive system needs THREE patches (declaration,
detection wiring, executor wiring) and corresponding
K-tests for each layer.  Missing the executor layer is
the most common gap because the analyzer-side change
"looks complete" without it.

#### Multi-site detection — fix every text-based site (v119.H18 lesson)

H18 fixed the AF_7mjs density-modify-and-stop regression
by switching `_apply_experiment_type_program_reprints`
from text-based detection to file-based detection (with
text as fallback).  The patch looked complete and the
6 new K-tests passed in the standalone H18 test file.

During the merge into the existing
`tst_density_modify_experiment_type.py` suite, K14 (the
exact AF_7mjs production input) failed: the correction
fired correctly (logs showed `[DIRECTIVE_CORRECTION]
source=files`), but the final `after_program` reverted to
`phenix.autobuild_denmod` anyway.

The cause: `_resolve_after_program` (the v115.10 post-LLM
overlay called from `_apply_workflow_intent_fallback`) had
its OWN text-only experiment-type heuristic.  For terse
advice ("density modify and stop"), its heuristic returned
no signal, defaulted `_exp="xray"`, then looked up
`_ACTION_TABLE["density_modify"]["xray"]` =
`phenix.autobuild_denmod` and unconditionally overwrote the
H18 correction.

**Lesson**: when fixing a categorical bug (like "wrong
experiment-type inference"), grep for every site that does
the same kind of inference.  Don't trust that fixing the
"obvious" site is sufficient.

Concrete grep commands that would have caught this earlier:

```bash
# Find every function that infers experiment_type from text:
grep -n 'is_cryoem\|is_xray\|cryo-?em\|x-?ray' agent/directive_extractor.py

# Find every site that defaults to a specific experiment type:
grep -n '_exp\s*=\s*"xray"\|_exp\s*=\s*"cryoem"' agent/directive_extractor.py

# Find every consumer of _ACTION_TABLE:
grep -n '_ACTION_TABLE\[' agent/directive_extractor.py
```

Either of the first two would have surfaced
`_resolve_after_program` as a second site with the same
flaw.

**The merge-to-existing-suite caught this; standalone
tests would not.**  K1-K13 (the pre-H18 §20 tests) all
used preprocessor-shaped advice with EXPLICIT
"Experiment Type: cryo-EM" / "Experiment Type: X-ray"
text markers.  The text-based heuristic in
`_resolve_after_program` got the right answer for those
inputs because the text contained the signal.  Only K14's
TERSE advice ("Primary Goal: Perform density modification
and stop", with NO experiment-type marker — matching the
real AF_7mjs preprocessor output) exposed the second-site
revert.

**Generalize**: when shipping a fix for a categorical bug,
add at least one K-test that exercises the SAME entry
point as the production failure (not the validator in
isolation).  And merge new K-tests into the existing
suite for that subsystem — separate test files miss
interactions between old and new code.

**See also (v119.H18.2)**: this lesson applies to function
CALLSITES too, not just internal heuristic sites.  H18 fixed
two callsites of `_resolve_after_program` but missed a third
that called the same function with a different surrounding
comment.  See the "Grep for ALL callsites when adding an
optional parameter" lesson below.

#### PHIL parameter declarations are per-master_params, not transitive (v119.H18.1 lesson)

H18 introduced a new PHIL parameter
`ai_analysis.original_files_for_directives` for threading the
input file inventory through to directive extraction.  The
declaration was added to `programs/ai_analysis.py` (the
analysis server's master_phil).  The USAGE site at
`programs/ai_agent.py:8513` then did:

```python
directive_params = copy.deepcopy(self.params)
directive_params.ai_analysis.original_files_for_directives = (
    ",".join(_basenames))
```

The 20 K-tests in `tst_density_modify_experiment_type.py`
passed because they called the helpers directly (with Python
dicts as inputs) — they never exercised the PHIL layer.

In production, every directive-extraction attempt crashed
with:

```
AttributeError: Assignment to non-existing attribute
  "ai_analysis.original_files_for_directives"
```

Directive extraction returned `{}` and the agent ran the
default cryo-EM plan stages through to predict_and_build
because no user-supplied stop directive ever made it into
the session.

**Root cause**: this codebase has TWO independent
master_params blocks:

- `programs/ai_agent.py` line 144 — the agent's PHIL.  This
  is what gets parsed into `self.params` when the user runs
  `phenix.ai_agent`.
- `programs/ai_analysis.py` line ~120 — the analysis server's
  PHIL.  This is what gets parsed when the analysis server
  is invoked as a sub-process or via the remote-agent path.

These two declarations are NOT transitive.  Adding a parameter
to `ai_analysis.py`'s master_params does NOT make it available
on `ai_agent.py`'s params object.  The agent calls
`copy.deepcopy(self.params)` to build `directive_params`,
which preserves the agent's PHIL schema.  Any assignment to
an attribute that schema doesn't declare raises
`AttributeError`.

**Lesson — PHIL declarations live with the consumer, not the
producer**.  When adding a new field:

1. **Identify EVERY entry point that produces a params object
   matching this scope.**  In this codebase, the question
   "where is `params.ai_analysis.<field>` ever extracted or
   assigned?" has TWO answers: ai_agent.py and ai_analysis.py.
2. **Add the PHIL declaration to EVERY one of them.**  Both
   must declare the field — they're independent schemas.
3. **Cross-reference the declarations in comments.**  A
   one-line `.help` note like "This PHIL definition must mirror
   the one in programs/<other>.py" prevents the gap from
   reopening during future refactors.

**The preventive K-test pattern** (`tst_h18_1_phil_roundtrip.py`):

```python
def test_master_params_parses_with_new_attribute():
    """Source-grep + actual PHIL parse + deep-copy + assignment.
    Catches the deploy gap at sandbox time."""
    # 1. Extract master_params from ai_agent.py
    with open("programs/ai_agent.py") as f:
        content = f.read()
    start = content.find('master_params = """') + 18
    end = content.find('"""', start)
    master_params_str = content[start:end]

    # 2. Parse with libtbx.phil (skip gracefully in sandbox
    #    without PHENIX)
    try:
        import libtbx.phil
    except ImportError:
        return  # sandbox-safe skip
    master_phil = libtbx.phil.parse(master_params_str)
    params = master_phil.extract()

    # 3. Deep-copy (mirroring ai_agent.py:8483)
    import copy
    directive_params = copy.deepcopy(params)

    # 4. Reproduce the production assignment
    directive_params.ai_analysis.new_field = "test_value"

    # 5. Verify it took
    assert directive_params.ai_analysis.new_field == "test_value"
```

Plus a source-grep variant for sandbox use (no PHENIX
required):

```python
def test_master_params_string_contains_param():
    """Sandbox-safe: confirm the PHIL declaration is present
    in the master_params region."""
    with open("programs/ai_agent.py") as f:
        content = f.read()
    start = content.find('master_params = """')
    end = content.find('master_phil = libtbx.phil.parse')
    master_region = content[start:end]
    assert "new_field" in master_region, (
        "PHIL declaration missing — deploy gap will crash "
        "production")
```

The source-grep variant runs anywhere — no libtbx required —
and is the load-bearing test for catching the gap before
deploy.

**Generalize: schema declarations don't propagate.**  This
lesson extends beyond PHIL to any place where data objects
are constructed against a schema:

- TypedDict subclasses (`AgentState` in `graph_state.py`)
- Dataclasses with explicit fields (`ErrorRecovery` in
  `error_analyzer.py`)
- Pydantic models (none currently in this codebase, but the
  pattern is the same)

For each, adding a new field requires updating the schema in
EVERY producer module.  When in doubt: grep for the type name
and audit each construction site.

#### Grep for ALL callsites when adding an optional parameter (v119.H18.2 lesson)

H18 added an optional `original_files` parameter to
`_resolve_after_program()` in `agent/directive_extractor.py`,
along with files-first experiment-type detection inside it.
The H18 audit identified TWO callsites that needed to pass
`original_files`:

1. `_apply_experiment_type_program_reprints` (the §20 site)
2. `_apply_workflow_intent_fallback` (the post-LLM overlay)

Both were updated.  H18's 20 K-tests passed.  The H18.1 PHIL
deploy gap was then fixed.  Tom re-ran AF_7mjs.  Production
STILL failed: the agent ran `phenix.predict_and_build` on
cycle 3 instead of stopping after `phenix.resolve_cryo_em`.

The runtime tracer (`h18_install_runtime_tracer.py`) revealed
the exact sequence:

```
[H18_TRACE] _apply_reprints: after_prog=None
[H18_TRACE] _apply_reprints EARLY RETURN: no after_prog
[H18_TRACE] AFTER _apply_experiment_type_program_reprints: after_prog=None
[H18_TRACE] before _apply_workflow_intent_fallback:
  stop_conds={..., 'after_program': 'phenix.autobuild_denmod', ...}
```

`after_program` was None when H18's first site ran, then set
to the WRONG value by the time the second site ran.  Something
between the two sites was writing `phenix.autobuild_denmod`.

That something was a THIRD callsite I missed.  At
`directive_extractor.py:783` lives a v117.2 fallback path that
fires when the LLM emits `stop_after_requested=True` but
omits `after_program`.  It calls `_resolve_after_program` to
fill in the missing field by parsing the raw advice:

```python
# Pre-H18.2 code at line 783:
_v172_source = raw_advice if raw_advice else user_advice
if _v172_source and _is_stop_after_requested(_v172_source):
    _resolve_after_program(directives, _v172_source.lower())
    #                                                       ↑
    #                                       MISSING: original_files
```

Without `original_files`, the resolver defaulted `_exp="xray"`
via the text-only heuristic ("density modify and stop" has no
cryo-EM/X-ray tokens) and mapped `denmod` →
`phenix.autobuild_denmod`.

The downstream `_apply_workflow_intent_fallback` DOES pass
`original_files`, but at THAT point the preprocessed advice
contains "Stop Condition: None", so
`_is_stop_after_requested(advice)` returns False, pushing the
resolver into the `n==1, no stop → leave as-is` branch.  The
buggy `after_program` from the v117.2 path persisted.

**Why H18's audit missed this**: the H18 plan correctly
identified the `_apply_experiment_type_program_reprints` site
and the `_apply_workflow_intent_fallback` site — both billed
as "experiment-type detection" sites in the codebase.  The
v117.2 fallback was billed differently: as a "fill in the
missing field" site.  Internally it called the same
`_resolve_after_program` function with the same files-win
contract — but the LABEL was different, so I never grepped
for it.

**The fix is one line** at the v117.2 callsite:

```python
# Post-H18.2:
_resolve_after_program(directives, _v172_source.lower(),
                       original_files=original_files)
```

**Lesson — when adding an optional parameter to a function,
grep for ALL callsites of that function regardless of what
the surrounding code is "doing"**.  The optional parameter's
contract is part of the function's identity; every caller
must opt in to the new contract or the bug travels through
the callers that didn't.

Concrete grep step that would have caught H18.2 at H18 time:

```bash
# Step 1 in every "add optional parameter to function X" PR:
grep -n "X(" *.py | grep -v "^.*def X"
```

Audit each line.  If the new parameter changes function
behavior in a way that matters at THAT callsite, pass it.  If
not, document why.  Either way, don't leave callsites
silently passing the legacy default.

**The K-test pattern that catches this class of bug** is the
production-faithful end-to-end test
(`tst_density_modify_experiment_type.py::K21`):

- Mock the LLM to return the EXACT shape Tom saw in
  production (`{"stop_conditions": {"stop_after_requested":
  true}}` — note: no `after_program`).
- Run the full `extract_directives()` pipeline.
- Assert on the final state, not on intermediate function
  outputs.

This pattern catches every internal site that touches
`after_program` because it tests behavior at the entry-point
contract, not at any individual function.

**Generalize the K-test pattern**: when a fix changes behavior
inside a complex pipeline, the regression test for that fix
should drive through the pipeline ENTRY POINT — not the
function the fix is in.  Internal-function tests catch
internal-function regressions; entry-point tests catch
multi-site regressions.  Both are necessary; entry-point is
the load-bearing one.

#### Sanity-check "missing input" — observed pattern, deliberately not auto-handled (v119 design note)

A `run_40_openai` batch scan revealed 15 of 21 logs with the
same failure shape: `EMPTY_RUN_NO_PROGRAM` on every variant
of the HNP3-mr-rosetta tutorial — across rules-only, LLM, and
think extraction modes.  Identical failure across LLM and
rules-only paths is a tell: the bug lives BELOW the directive-
extraction layer.

Root cause: the tutorial ships its data inside a single
`mr_rosetta.tar.gz` archive (README documents this).  File
discovery finds the tarball but has no visibility into its
contents.  PERCEIVE correctly notes "no MTZ, no model, no
map" and sanity-aborts.

**Both behaviors are individually correct**.  File discovery
correctly enumerates the directory.  PERCEIVE correctly
sanity-aborts when there's no data.  The case is interesting
because the `.tar.gz` is INPUT DATA from the user's
perspective, INVISIBLE from the discovery layer's
perspective, and ABSENT from PERCEIVE's perspective.

**Design decision: don't auto-extract.**  We considered
adding an auto-extraction hook between `_auto_discover_files`
and `set_project_info` that would unpack the archive in place
when no usable data was present alongside it.  Prototyped it,
wrote a bounded-contract module with path-escape defense,
symlink refusal, no-overwrite, and 15 K-tests.  Then decided
not to ship.

The reasoning:

1. **Real users don't ship tarballs to the agent.**  In
   production, users have a directory with `.mtz` and `.pdb`
   files already laid out.  The HNP3 case is an artifact of
   how the PHENIX tutorial corpus packages its examples — a
   tutorial-distribution decision, not an agent-design
   constraint.
2. **Production code with no production caller is a
   liability.**  Auto-extracting archives carries real risk
   surface — path escapes, symlink attacks, clobbered
   files, runaway extraction on giant archives — that
   would need to be maintained forever for a problem only
   the batch test suite hits.
3. **The right fix is at the test-harness level.**  Either
   pre-extract HNP3's tarball before invoking the agent in
   the batch suite, or exclude HNP3-mr-rosetta from the
   batch corpus with a comment explaining the
   incompatibility.  Both keep the agent codebase clean.
4. **Failure mode is already legible.**  The agent's
   sanity-abort message names exactly what's missing
   ("No reflection data file found for X-ray workflow")
   and points at the directory.  A human looking at the
   log can see `mr_rosetta.tar.gz` and infer the fix in
   seconds.

**What to do when this pattern recurs:**

If a real user reports the same failure (data in a tarball
the agent can't see), the right response is:

1. First, ask: is this a real workflow or a tutorial
   artifact?  Real workflows almost never look like this.
2. If real: tell the user to extract the archive and re-run.
   The agent will pick up the files on the second pass.
3. Only if the same case recurs frequently from genuine
   user workflows should auto-extraction be revisited.  At
   that point, the prototype design (bounded trigger,
   no-overwrite, path-escape defense, symlink refusal,
   idempotent, never-raises) is the starting point.

**The general principle (worth keeping as a design lens):**
when a sanity check fires on missing inputs, ask TWO
questions:

1. Is the input genuinely absent?  (User forgot to provide
   it.)
2. Is the input present but invisible?  (Compressed in a
   tarball, in a subdirectory, behind a symlink to an
   unusual path, etc.)

If the second is plausible AND it's a real user pattern (not
a test-harness artifact), the discovery layer should either
expand visibility or surface an actionable error that
distinguishes between the two.  Don't conflate "I don't see
it" with "it's not there."

**Generalize beyond data ingest**: the "absent vs invisible"
distinction applies anywhere a pipeline reports "X missing".
Authentication systems should distinguish "no credentials"
from "credentials expired."  Compilers should distinguish
"header not on path" from "header exists but include guard
didn't expand."  Error messages that say "X missing" without
distinguishing the two confuse users.

But — and this is the key takeaway from the HNP3 case — the
principle doesn't always require implementation.  Sometimes
the right move is to recognize the pattern, document it,
fix the test harness, and keep production code lean.

### Commit quality

Before committing LLM-generated code:

- Run the full test suite, not just the tests the LLM
  wrote. Changes in shared modules can break downstream
  code that the LLM never saw.
- Review the diff, not just the final file. LLMs
  occasionally make unrelated "cleanup" changes that
  alter behavior.
- Write a commit message that describes the change, not
  the process. "Fix None crash in space group lookup"
  is a good message. "Changes suggested by Claude" is
  not.

---

### 4. Coding Standards (Summary)

The full coding standards are in the two companion
guidelines documents. The programmer should be familiar
with these and enforce them during review:

**Style**: 2-space indentation, 80-character line width,
descriptive names, no trailing whitespace, no unused
imports.

**None-safety**: `(d.get("key") or "").lower()`, never
`d.get("key", "").lower()` — the latter crashes when the
key exists with value `None`.

**Error handling**: No bare `except:` (catches
`SystemExit`). Every `except Exception: pass` gets a
comment explaining why. Non-critical functions use the
"Never raises" pattern with logging.

**Paths**: `os.path.join()` always. Never string
concatenation with `/`.

**Tests**: Fail-fast style, `approx_equal` for floats,
explicit edge-case coverage.

**Serialization**: `to_dict()` / `from_dict()` must
round-trip. `from_dict()` must tolerate missing keys.

**Agent-specific** (see `AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md`):
Verify strategy flags against `programs.yaml` before use.
Check all three error classification systems when adding
error patterns. Identify client vs server code path for
changes to `programs/ai_agent.py`. Ensure session state
fields appear in `session.data`, `create_initial_state()`,
and `contract.py`.

---

### 5. When NOT to Use LLMs

LLMs are effective for:

- Implementing well-defined algorithms with clear
  specifications
- Writing test suites for existing code
- Refactoring for consistency (renaming, reformatting,
  extracting functions)
- Debugging with specific error messages as input
- Writing documentation for existing code
- Boilerplate (PHIL definitions, serialization, CLI
  wrappers)
- Analyzing large codebases to find patterns, inconsistencies,
  or undocumented behavior

LLMs are unreliable for:

- **Novel scientific algorithms** where correctness
  depends on domain knowledge the LLM may not have.
  Use LLMs to implement the algorithm you've designed,
  not to design it.
- **Performance-critical inner loops** where the LLM
  can't profile or benchmark. Write the inner loop
  yourself; use the LLM for the surrounding
  infrastructure.
- **Security-sensitive code** (authentication, file
  permissions, subprocess sandboxing). LLMs produce
  code that works but may have subtle vulnerabilities.
- **Guessing at APIs they haven't seen.** If the LLM
  doesn't have the source code in context, it will
  invent plausible-looking function signatures that
  don't exist.

In these cases, write the critical code yourself and use
the LLM for the non-critical surrounding code (tests,
documentation, boilerplate).

---

### 6. Prompt Reference Card

Keep these prompts handy during LLM sessions:

| When | Prompt |
|------|--------|
| Start | Attach both guidelines + source files |
| Before coding | *"Create a plan as a markdown file"* |
| Plan review | *"What could go wrong with this?"* |
| Each step | *"Any other fixes before we move on?"* |
| Each step | *"Run all relevant tests"* |
| All done | *"Run a full audit of all changes"* |
| Final | *"Check changed files for Common Pitfalls"* |

---

### 7. Checklist for Every Pull Request

Before committing LLM-assisted code, the programmer
confirms:

- [ ] A plan was created and reviewed before coding began
- [ ] Each step was reviewed individually, not batched
- [ ] "Any other fixes?" was asked after each step
- [ ] All changed files parse cleanly
- [ ] The full test suite passes (not just new tests)
- [ ] The diff contains only intentional changes
- [ ] No hardcoded path separators
- [ ] No bare `except:` statements
- [ ] No `d.get("key", "").lower()` on potentially-None
      values
- [ ] `approx_equal()` used for all floating-point
      comparisons in tests
- [ ] Scientific units are correct (Å, fractions, degrees)
- [ ] Agent-specific: strategy flags verified against
      `programs.yaml`
- [ ] Agent-specific: new error patterns checked against
      all three classification systems
- [ ] Agent-specific: session state fields in all three
      locations (session.data, create_initial_state,
      contract.py)
- [ ] Documentation updated for user-visible changes
- [ ] Commit message describes the change, not the tool

---

## 11. Worked Example: Standard MR

### Worked Example: Standard X-ray MR (Scenario 1)

This traces the core decision logic for a standard X-ray molecular
replacement run. For full cycle-by-cycle detail see
[`docs/reference/traces/`](reference/traces/).

**Input:** `raster.mtz`, `raster.fa`, `raster.pdb` (search model)
**Mode:** `maximum_automation=True`

**Key decision points:**

1. **Data analysis** — `phenix.xtriage` runs first because `xtriage_done=False`.
   No anomalous signal → MR pathway. Twinning check passes.

2. **Model placement** — `phenix.phaser` (MR). `_has_placed_model()` uses
   unit cell comparison (Tier 1) — placed model detected when cell matches.
   Subsequent cycles gate refinement on placement status.

3. **Refinement loop** — `phenix.refine` cycles until R-free plateau or
   `max_refine_cycles` limit. `should_pivot()` fires after 3+ RETRYABLE
   failures of the same program.

4. **Validation** — `phenix.molprobity` + `phenix.model_vs_data` when
   R-free < 0.35 and `validation_done=False`.

5. **Stop** — either convergence (R-free stable, validation passed) or
   safety stop (R-free stuck above 0.45 after 8+ cycles).

**Files involved:** `workflow_engine.py` (step detection),
`graph_nodes.py` (PLAN/BUILD nodes), `directive_validator.py` (stop logic).

### Key Decision Points

### 1. Experiment Type Detection

Based on file categories:
- Has `.mtz` + no maps → **X-ray**
- Has `.mrc`/`.ccp4` maps → **Cryo-EM**

### 2. Dual MTZ Tracking (v110)

| Category | Update Rule | Used By |
|----------|-------------|---------|
| `data_mtz` | First with R-free **locks forever** | phenix.refine, phenix.phaser |
| `map_coeffs_mtz` | **Most recent wins** | phenix.ligandfit |

This prevents refinement from using its own output maps as input data.

### 3. Automation Path (v110)

| Mode | `maximum_automation` | predict_and_build Behavior |
|------|---------------------|---------------------------|
| Automated | True (default) | Runs full workflow |
| Stepwise | False | Stops after prediction |

### 4. Anomalous Workflow Triggering

`phenix.autosol` becomes available when:
- `has_anomalous` = True (from xtriage or history)
- `has_sequence` = True
- `autosol_done` = False

The `use_experimental_phasing` directive prioritizes autosol over predict_and_build.

### 5. Plateau Detection

From `workflows.yaml`:
```yaml
repeat:
  max_cycles: 4
  until:
    any:
      - r_free: "< target_r_free"
      - condition: plateau
        cycles: 2
        threshold: 0.005
```

Agent stops refinement when:
- Target R-free reached, OR
- 2 consecutive cycles with < 0.005 improvement

### 6. File Selection Logic

**Best model selection:**
1. Check `best_files["model"]` from BestFilesTracker
2. Prefer by stage score: refined (100) > autobuild_output (100) > predict_build_output (90)
3. Among equal scores, prefer better R-free

**MTZ selection:**
- For refinement: Use locked `data_mtz` (original with R-free flags)
- For ligandfit: Use `map_coeffs_mtz` (latest calculated phases)

**R-free generate guard (v120):** `command_builder.py`'s phenix.refine invariant
adds `xray_data.r_free_flags.generate=True` ONLY on positive evidence there is no
test set, and never overwrites an existing one.  Because the server cannot read
client file paths (`files_local=False`), the input-MTZ R-free check is
**client-extracted**: `ai_agent.py` inspects the original data MTZ (from
`session.data["original_files"]`, via `classify_mtz_type`/`inspect_mtz`), caches the
result once, and ships a tri-state fact `input_mtz_has_rfree`
(True/False/None) through the transport (`build_session_state` →
`build_request_v2` whitelist → `run_ai_agent` map-back, each with an explicit
`is not None` guard so a confirmed-`False` is not dropped).
`CommandBuilder._input_mtz_rfree_state(files, context)` reads, in order: the client
fact (authoritative, works server-side); `context.mtz_inspection` positive-only; and
`inspect_mtz` of the selected data MTZ **only when `files_local`** (so a server read
can't masquerade as a confident False).  The decision is three-valued: generate only
when state is `False`; strip/suppress when `rfree_mtz` is locked, state is `True`, or
state is `None` (undetermined → conservative, since silently overwriting a test set
is worse than a recoverable `rfree_flags_missing` error).  Field registered in
`contract.py` (v7).  Pinned by `tst_rfree_generate_guard.py` and
`tst_input_mtz_has_rfree_plumbing.py`.

---

### Edge Cases

### What if predict_and_build fails?

The agent would:
1. See error in log
2. If anomalous signal detected: try `phenix.autosol`
3. If template available: try `phenix.phaser`
4. Eventually STOP with error if no path forward

### What if user provides a pre-built model?

If `model.pdb` is provided alongside `data.mtz`:
- File categorization checks filename patterns
- If looks refined: skip to refinement phase
- If looks like template: use phaser to place it

### What if ambiguous data arrays?

Error recovery system:
1. Detects "Multiple equally suitable arrays" error
2. Selects appropriate array based on context (merged vs anomalous)
3. Retries program with explicit label selection
4. Max 3 retries before giving up

### What if ligand file is named ambiguously?

From `file_categories.yaml`:
```yaml
ligand:
  extensions: [.pdb, .cif]
  prefer_patterns: [lig, ligand]
```

Files matching `lig*` or `ligand*` are categorized as ligands, not protein models.

**Full scenario traces** (stepwise, ligand, SAD, restart):
[`docs/reference/traces/`](reference/traces/)

---

## 12. Common Development Tasks

### 4a. Adding a New PHENIX Program

This is the most common development task. For most
programs, you only need to edit two files:

1. **`knowledge/programs.yaml`**: Define the program's
   inputs, outputs, command template, metric
   extraction, and done_tracking.
2. **`knowledge/workflows.yaml`**: Add the program to
   the appropriate workflow step(s) with conditions.

Optionally:

3. **`knowledge/plan_templates.yaml`**: Add the
   program to relevant plan template phases (if it
   should be part of the expert-mode strategy).
4. **`knowledge/file_categories.yaml`**: Define new
   output file types if the program produces files
   the agent hasn't seen before.

The system automatically handles metric extraction,
done flags, session summary display, and LLM
hallucination sanitization.

See §4 (Adding a New Program) for the
complete guide with examples.

### 4b. Adding a Plan Template

**When you need a new template:**
- The agent has a blank plan for a combination of
  inputs you want to support (check by running
  `generate_plan()` with your test inputs)
- An existing template's phase structure is wrong
  for a specific scenario

**When to modify an existing template instead:**
- You just need to adjust thresholds or add a program
  to an existing phase
- You want to add a gate condition to an existing
  phase

**Template anatomy** (see Section 3d for the
annotated YAML). Key points:

- `applicable_when` conditions are matched against
  the session context. Higher score (more matching
  conditions) wins. Ties broken by `priority`.
- `success` criteria → gate says "advance"
- `gate` conditions → gate says "retreat" or "stop"
- `extends` lets a template inherit from another and
  override specific phases

**Testing template selection:**

```python
from knowledge import plan_template_loader
plan_template_loader._templates_cache = None
from agent.plan_generator import generate_plan

plan = generate_plan(
    available_files=['data.sca', 'seq.dat'],
    user_advice='SAD with selenium',
    directives={'program_settings': {
        'phenix.autosol': {'atom_type': 'Se'}}})
print(plan.template_id)  # → "sad_phasing"
```

### 4c. Modifying Decision Logic

**Where program selection happens:**

1. `workflow_engine.py:get_valid_programs()` →
   determines which programs CAN run (rule-based)
2. `graph_nodes.py:plan()` → selects which program
   WILL run (LLM or rules-based fallback)
3. `graph_nodes.py:validate()` → confirms the
   selection is safe

To change which programs are available in a given
state, edit `workflows.yaml`. To change how the agent
chooses between available programs, modify the LLM
prompt in `knowledge/prompts.py` or the rules-based
selector in `graph_nodes.py`.

**Adding a new workflow step:**

1. Add the phase to `workflows.yaml` under the
   appropriate experiment type
2. Define the phase's programs, conditions, and
   transitions
3. Add transition rules from existing phases to
   the new one
4. Test with `detect_workflow_state()` to verify
   the phase is reached

### 4d. Adding Expert KB Entries

**File:** `knowledge/expert_knowledge_base_v2.yaml`

Each entry has tags (matched against session context
using IDF-weighted scoring) and a guidance text
(injected into the THINK prompt when matched).

```yaml
- tags: [r_free_stuck, refinement, convergence]
  guidance: |
    R-free stalling during refinement often indicates
    a problem with the model rather than the data.
    Consider running autobuild to rebuild problematic
    regions before continuing refinement.
```

Tags are matched against the current session state
(programs run, metrics, workflow step). IDF weighting
means rare tags (like "twinning") get more weight
than common ones (like "refinement").

Test that entries fire by running the agent with
`verbosity=verbose` and looking for `[KB rules
consulted]` in the Expert Assessment output.

### 4e. Debugging Agent Behavior

**Tools available:**

| Tool | What it shows |
|------|--------------|
| `verbosity=verbose` | All decisions, file selection, valid programs, directive modifications |
| `display_and_stop=detailed` | Full session history without running new cycles |
| `dry_run=True` | Simulated execution with predefined outcomes |
| `session.json` | Complete session state (open in a JSON viewer) |
| Event log (verbose) | PERCEIVE/PLAN/BUILD debug lines |

**Common debugging patterns:**

- **"Why did the agent choose program X?"** → Run
  with `verbosity=verbose` and look for `PLAN:`
  lines. Check `valid_programs` in PERCEIVE output.
- **"Why is program X not available?"** → Check
  `workflows.yaml` conditions and the program's
  `done_tracking` flag in session data.
- **"Why did the gate evaluate to retreat?"** →
  Look for `[GATE]` lines in the log. Check the
  phase's success criteria and gate conditions.
- **"Why is the agent stuck in a loop?"** → Check
  for `LOOP WARNING` in PERCEIVE output. Look at
  the valid_programs list across consecutive cycles.

---

## 13. Configuration Reference

### Programs (`knowledge/programs.yaml`)

| Program | Category | Experiment Types |
|---------|----------|-----------------|
| phenix.xtriage | analysis | xray |
| phenix.mtriage | analysis | cryoem |
| phenix.phaser | model_building | xray |
| phenix.autosol | model_building | xray |
| phenix.autobuild | model_building | xray |
| phenix.predict_and_build | model_building | xray, cryoem |
| phenix.process_predicted_model | model_building | xray |
| phenix.dock_in_map | model_building | cryoem |
| phenix.map_to_model | model_building | cryoem |
| phenix.refine | refinement | xray |
| phenix.real_space_refine | refinement | cryoem |
| phenix.ligandfit | ligand | xray, cryoem |
| phenix.autobuild_denmod | map_improvement | xray |
| phenix.resolve_cryo_em | map_optimization | cryoem |
| phenix.map_sharpening | map_optimization | cryoem |
| phenix.molprobity | validation | xray, cryoem |
| phenix.model_vs_data | validation | xray |
| phenix.validation_cryoem | validation | cryoem |
| phenix.holton_geometry_validation | validation | xray, cryoem |
| phenix.map_correlations | validation | xray, cryoem |
| phenix.map_symmetry | analysis | cryoem |
| phenix.polder | map_analysis | xray |
| phenix.pdbtools | utility | xray, cryoem |

*This table is derived from `programs.yaml`. Edit the
YAML, not this table.*

### Plan Templates (`knowledge/plan_templates.yaml`)

| Template | Applicable When | Priority |
|----------|----------------|----------|
| mr_refine | xray + search model | 10 |
| mr_refine_ligand | xray + search model + ligand | 20 |
| mr_refine_lowres | xray + search model + >3.0Å | 15 |
| mr_refine_highres | xray + search model + <1.5Å | 15 |
| mr_refine_twinned | xray + search model + twinned | 25 |
| predict_refine | xray + sequence, no model, no anomalous | 10 |
| predict_refine_ligand | xray + sequence + ligand, no model | 20 |
| mr_sad | xray + search model + anomalous | 40 |
| sad_phasing | xray + anomalous (no search model) | 30 |
| sad_phasing_ligand | xray + anomalous + ligand | 35 |
| cryoem_refine | cryo-EM | 10 |
| cryoem_refine_ligand | cryo-EM + ligand | 20 |

*This table is derived from `plan_templates.yaml`.*

---

## 14. Testing (Legacy Notes)

The test suite uses cctbx-style fail-fast behavior
with plain functions. See §5 (Testing)
for the complete testing guide.

**Quick reference:**

```bash
# Run all tests
python tests/run_all_tests.py

# Quick mode (standalone only, no PHENIX required)
python tests/run_all_tests.py --quick

# Single file
python tests/tst_plan_generator.py

# Pattern match
python tests/run_all_tests.py --pattern "plan"
```

**Key test suites for the planning layer:**

| Suite | Tests | What it covers |
|-------|-------|---------------|
| tst_plan_schema | 53 | StageDef, StructurePlan serialization |
| tst_plan_generator | 39 | Template selection, context building |
| tst_gate_evaluator | 54 | Success/retreat/advance logic |
| tst_structure_model | 77 | Model state tracking, round-trips |
| tst_validation_history | 56 | Per-cycle metric recording |
| tst_thinking_defense | 55 | Source code contracts |
| tst_display_data_model | 41 | DDM properties, HTML report generation |
| tst_scenario_tracer | 57 | End-to-end scenario testing (see below) |

**Scenario tracer** (`tests/tst_scenario_tracer.py`):

Calls real agent functions with mock data at each
decision point: `generate_plan`, `GateEvaluator`,
`build_thinking_prompt`, `parse_intent_json`,
`parse_assessment`. 57 scenarios covering:

| Group | Count | What |
|-------|-------|------|
| S1-S15 | 27 | Crystallographic scenarios (template + gate + prompt) |
| C1-C3 | 3 | Cycle counting edge cases |
| G1-G4 | 4 | Gate numeric logic |
| M1-M2 | 2 | Multi-cycle plan progression |
| P1-P6 | 6 | PHENIX workflow state (needs libtbx) |
| L1-L10 | 10 | Mock LLM output parsing and validation |
| PG1-PG5 | 5 | Model placement gate (v114.1) |

```bash
python tests/tst_scenario_tracer.py            # All
python tests/tst_scenario_tracer.py S1A        # One
python tests/tst_scenario_tracer.py --list     # List
python tests/tst_scenario_tracer.py --failures # Failures only
```

**Tutorial run analyzer**
(`tests/analyze_tutorial_runs.py`):

Compares agent performance across modes. Reads
`agent_session.json` from tutorial directories named
`{tutorial}__{mode}` (double underscore), produces
CSV, Markdown, and JSON reports.

```bash
python tests/analyze_tutorial_runs.py /path/to/results/
```

Supports 5 modes: `rules_only`, `llm`, `llm_think`,
`llm_think_advanced`, `llm_think_expert`.

**When adding tests:**

1. Create `tests/tst_your_module.py`
2. Use assert helpers from `tests/tst_utils.py`
3. Add a `run_tests()` or `run_all_tests()` function
4. Register in `tests/run_all_tests.py`
5. Guard any file reads outside `agent/`/`knowledge/`
   with `os.path.isfile()`

---

## 15. Future Directions

### Evaluation Harness

Tutorial run analyzer is operational. Next steps:
expand to 30+ tutorials, establish baseline metrics,
run regression tests on code changes.

### GUI Enhancements

- **Pinned plan header**: Replace the bare
  `wx.TextCtrl` progress panel with a composite
  widget showing live ✓/●/○ phase indicators.
- **Documentation tooltips**: "?" icons next to Expert
  Assessment metrics linking to the User Guide.
- **Interactive plan editing**: Let the user inspect
  and modify the plan before execution.

### Template Coverage

- MR-SAD + ligand combined template
- Neutron diffraction support
- Multi-dataset processing
- Phasing with non-crystallographic symmetry

### Agent Intelligence

- **Hypothesis engine expansion**: More types (wrong
  space group, register shift errors, twinning
  detection from metrics).
- **Learning from failures**: Staging unverified rules
  learned from failed runs for human review before
  promotion to the KB.
- **Multi-crystal support**: Handle multiple crystals
  or datasets in a single session with shared
  decision context.
- **Structured results via `results_as_json()`**:
  Newer PHENIX programs built on `ProgramTemplate`
  expose a `results_as_json()` method returning
  metrics as structured JSON. Migrating the agent
  from regex log parsing to JSON results would
  eliminate extraction fragility, provide richer
  data (per-residue validation, per-chain stats),
  and reduce per-program integration cost. The
  migration is incremental — programs that support
  it can switch one at a time with log parsing as
  fallback. See ARCHITECTURE.md "Potential
  improvements" for details.

---

*This guide covers PHENIX AI Agent v114. For the
latest information, see the PHENIX documentation
and the in-tree docs (ARCHITECTURE.md, CHANGELOG.md,
guides/).*

*Reference tables in Section 5 are derived from YAML
configuration files. When YAML changes, regenerate
these tables — the YAML is the single source of
truth.*
