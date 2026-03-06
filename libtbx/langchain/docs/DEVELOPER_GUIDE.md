# PHENIX AI Agent — Developer Guide

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
   current workflow phase (e.g., "analyze", "refine",
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
     `prefer_programs` from the current plan phase

3. The resulting list is passed to the LLM (or
   rules-based selector) which picks one.

**Phase transitions** happen when a program sets its
done flag and the phase's transition conditions are
met. For example, completing xtriage in the "analyze"
phase transitions to "obtain_model".

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
- `generate_phase_summary()`: LLM-based, at phase
  transitions
- `generate_final_report()` /
  `generate_stopped_report()`: template-based, at
  session end

### 3e. Model Placement Gate (v114.1)

**Files:** `programs/ai_agent.py`
(`_detect_model_placement`,
`_skip_plan_phases_for_placement`),
`agent/workflow_engine.py`
(`model_is_placed_confirmed` in context)

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

**Detection thresholds:** CC > 0.3 for
model_vs_data, R-free < 0.50 for refine,
map_cc > 0.3 for real_space_refine. These are
deliberately low — the question is "does the model
correspond to the data at all?" not "is the model
good?"

**Plan fast-forward:** When placement is detected,
pending `molecular_replacement` and
`experimental_phasing` phases are marked "skipped"
and the plan advances to the next actionable phase.

**HETATM detection:** At startup,
`_scan_input_models_for_ligands` scans input PDB
files for non-water HETATM records, sets
`session.data["input_has_ligand"]`. This enables
polder without a ligandfit step when the model
already contains ligands.

### 3f. Client-Server Contract

**Files:** `agent/contract.py`,
`docs/guides/BACKWARD_COMPATIBILITY.md`

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

See `docs/guides/BACKWARD_COMPATIBILITY.md` for the
full contract system.

### 3f. Error Recovery

**Files:** `agent/error_analyzer.py`,
`programs/ai_agent.py` (`_handle_execution_result`)

When a PHENIX program fails, the agent attempts
automatic recovery:

1. **Probe programs** (`phenix.model_vs_data`): Run
   a quick diagnostic to gather information.
2. **Auto-recovery**: Match the error against known
   patterns and inject corrective parameters (e.g.,
   reduce macro-cycles, change strategy).
3. **Terminal diagnosis**: If the error is
   unrecoverable, generate an HTML diagnosis page
   explaining what went wrong and how to fix it.
4. **Continue**: If none of the above applies,
   annotate the failure and let the agent try a
   different program on the next cycle.

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
| `agent_gate_transition` | Phase change | Action, phases, reason |
| `sub_job_complete` | Program finished | Job info for project DB |

**Display panels:**

- `ConfigPanel`: Setup (files, advice, settings)
- `RunNotebook` → `_progress_panel`: Live progress
  (text append only)
- `ResultPanel`: Structured results after completion

All new callback fields must use
`getattr(data, 'field', default)` on the GUI side
for backward compatibility. All callback sends must
be wrapped in `try/except Exception: pass` (non-
critical, must not crash the agent).

---

## 4. Common Development Tasks

### 4a. Adding a New PHENIX Program

This is the most common development task. For most
programs, you only need to edit two files:

1. **`knowledge/programs.yaml`**: Define the program's
   inputs, outputs, command template, metric
   extraction, and done_tracking.
2. **`knowledge/workflows.yaml`**: Add the program to
   the appropriate workflow phase(s) with conditions.

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

See `docs/guides/ADDING_PROGRAMS.md` for the
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

**Adding a new workflow phase:**

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
(programs run, metrics, workflow phase). IDF weighting
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

## 5. Configuration Reference

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

## 6. Testing

The test suite uses cctbx-style fail-fast behavior
with plain functions. See `docs/guides/TESTING.md`
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
| tst_plan_schema | 53 | PhaseDef, StructurePlan serialization |
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

## 7. Future Directions

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

---

*This guide covers PHENIX AI Agent v114. For the
latest information, see the PHENIX documentation
and the in-tree docs (ARCHITECTURE.md, CHANGELOG.md,
guides/).*

*Reference tables in Section 5 are derived from YAML
configuration files. When YAML changes, regenerate
these tables — the YAML is the single source of
truth.*
