# User Directives System

## Overview

The user directives system provides a "trust but verify" approach for user advice in the PHENIX AI Agent. The system uses an LLM to extract structured directives from natural language user advice, stores them persistently across cycles, and applies them through a clean layered architecture.

## Goals

1. Allow users to give complex, nuanced instructions in natural language
2. Ensure instructions are followed consistently across all cycles
3. Catch and correct LLM mistakes that contradict user intent
4. Respect user stop conditions even when they conflict with validation requirements

## Architecture

The directive system follows a layered decision flow:

```
User Advice (natural language)
            │
            ▼
┌─────────────────────────────────┐
│  PHASE 1: DIRECTIVE EXTRACTION  │
│  (Once at session start, or     │
│   when advice changes)          │
│                                 │
│  LLM parses user advice into    │
│  structured JSON directives     │
│  directive_extractor.py         │
└─────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────┐
│  STORED DIRECTIVES              │
│  (Persistent in session.json)   │
│                                 │
│  - program_settings             │
│  - stop_conditions              │
│  - file_preferences             │
│  - workflow_preferences         │
│  - constraints                  │
└─────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────┐
│  PHASE 2: WORKFLOW ENGINE       │
│  (Each cycle)                   │
│                                 │
│  _apply_directives() modifies   │
│  valid_programs list:           │
│  - Add STOP if skip_validation  │
│  - Add after_program target     │
│  - Apply skip/prefer programs   │
│  workflow_engine.py             │
└─────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────┐
│  PHASE 3: LLM PLANNING          │
│  (Each cycle)                   │
│                                 │
│  LLM sees raw user advice +     │
│  extracted directives + context │
│  and outputs intent             │
│                                 │
│  directive_validator applies    │
│  program settings to strategy   │
│  graph_nodes.py                 │
└─────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────┐
│  PHASE 4: POST-EXECUTION CHECK  │
│  (After each program runs)      │
│                                 │
│  ONLY place stop conditions     │
│  are evaluated:                 │
│  - Check after_program match    │
│  - Check after_cycle reached    │
│  - Check metric targets         │
│  ai_agent._run_single_cycle()   │
└─────────────────────────────────┘
            │
            ▼
        Final Decision
        (continue or stop)
```

## Directive Schema

```json
{
  "program_settings": {
    "phenix.autosol": {
      "resolution": 3.0,
      "atom_type": "Se"
    },
    "phenix.refine": {
      "anisotropic_adp": true
    },
    "default": {
      "resolution": 2.5
    }
  },

  "stop_conditions": {
    "after_program": "phenix.refine",
    "start_with_program": "phenix.polder",
    "after_cycle": 4,
    "max_refine_cycles": 2,
    "skip_validation": true,
    "r_free_target": 0.25,
    "map_cc_target": 0.70
  },

  "file_preferences": {
    "model": "beta.pdb",
    "sequence": "beta.seq",
    "exclude": ["old_model.pdb"]
  },

  "workflow_preferences": {
    "skip_programs": ["phenix.autobuild"],
    "prefer_programs": ["phenix.phaser"],
    "use_experimental_phasing": true,
    "use_molecular_replacement": false,
    "use_mr_sad": false
  },

  "constraints": [
    "Do not add waters until R-free < 0.30",
    "Use TLS refinement after cycle 3"
  ]
}
```

### Stop Condition Fields

| Field | Description |
|-------|-------------|
| `after_program` | Stop after this program completes (final goal, unconditional) |
| `start_with_program` | Run this program first, then continue workflow |
| `after_cycle` | Stop after this cycle number |
| `max_refine_cycles` | Limit refinement to N cycles, then proceed to validation (controlled landing) |
| `skip_validation` | Skip molprobity validation at end |
| `r_free_target` | Stop when R-free reaches this value |
| `map_cc_target` | Stop when map CC reaches this value |

**`max_refine_cycles` vs `after_program`** — these are fundamentally different:

- **`max_refine_cycles: N`** — "Run at most N refinement cycles, then
  proceed to validation." When the limit is reached the engine removes
  refinement from the valid-program list and injects the validate-phase
  programs (`phenix.molprobity`, `phenix.model_vs_data` for X-ray;
  `phenix.molprobity`, `phenix.validation_cryoem` for cryo-EM) plus
  `STOP`. The user gets a quality report before finishing.

- **`after_program: "phenix.refine"`** — "Run until this program
  completes, then stop unconditionally." No validation programs are
  injected; `STOP` is the only option after the named program finishes.
  Use when you want to check MR worked before committing to a full run.

### Program Settings Enable Programs

When you specify `program_settings` for a program, that program is automatically added to the list of valid programs (if prerequisites are met). This allows you to run earlier-phase programs even when the workflow state has progressed past that phase.

**Example**: If you have a model file but want to run predict_and_build anyway:

```json
{
  "program_settings": {
    "phenix.predict_and_build": {"rebuilding_strategy": "Quick"}
  }
}
```

The agent will add `phenix.predict_and_build` to valid programs even if the workflow state is `xray_refined` (because a model already exists).

**Prerequisites are still checked:**
- `phenix.refine` / `phenix.real_space_refine`: Requires a model to refine
- `phenix.ligandfit`: Requires prior refinement (for map coefficients)
- `phenix.predict_and_build`: Always allowed (worst case: prediction-only)

## Files

### Core Directive System Files

| File | Description |
|------|-------------|
| `agent/directive_extractor.py` | LLM-based directive extraction, `check_stop_conditions()` |
| `agent/directive_validator.py` | Validates and augments LLM intent with program settings |
| `agent/workflow_engine.py` | `_apply_directives()` modifies valid_programs list |
| `agent/session.py` | Directive storage and retrieval methods |

### Test Files

| File | Description |
|------|-------------|
| `tests/tst_directive_extractor.py` | 73 unit tests for extraction |
| `tests/tst_directive_validator.py` | 53 unit tests for validation |
| `tests/tst_session_directives.py` | 12 unit tests for session methods |
| `tests/tst_directives_integration.py` | 16 integration tests |
| `tests/tst_decision_flow.py` | Decision flow architecture tests |

### Modified Files

| File | Changes |
|------|---------|
| `agent/graph_nodes.py` | Simple validation gate in `plan()`, no stop logic |
| `agent/graph_state.py` | Added `directives: Dict` to AgentState TypedDict and `create_initial_state()` |
| `agent/workflow_state.py` | Pass directives through to engine |
| `agent/api_client.py` | Include directives in REST request/response |
| `knowledge/api_schema.py` | Added directives to session_state schema |
| `knowledge/prompts_hybrid.py` | Show extracted directives in LLM prompt |
| `programs/ai_agent.py` | Post-execution stop check in `_run_single_cycle()` |
| `phenix_ai/local_agent.py` | Include directives in session_state |
| `phenix_ai/remote_agent.py` | Include directives in session_state |
| `phenix_ai/run_ai_agent.py` | Extract and pass directives to graph |

> **Important**: The `directives` field must be defined in the `AgentState` TypedDict 
> (`agent/graph_state.py`) for LangGraph to properly pass it between nodes. Missing 
> TypedDict fields are silently dropped by LangGraph's StateGraph.

## Data Flow

### Complete Flow (Local and Remote Identical)

```
Session.extract_directives()
    ↓
session.data["directives"]
    ↓
ai_agent.py: session.get_directives() → session_info["directives"]
    ↓
local_agent.py / remote_agent.py: session_info["directives"] → session_state["directives"]
    ↓
build_request_v2(): session_state["directives"] → normalized_session_state["directives"]
    ↓
create_request(): session_state passed through
    ↓
[JSON over REST if remote]
    ↓
run_ai_agent.py: session_state["directives"] → directives variable
    ↓
create_initial_state(directives=...) → state["directives"]
    ↓
workflow_engine._apply_directives(): Modifies valid_programs
graph_nodes.py plan(): 
  - directive_validator.validate_intent() applies program settings
  - Simple validation: is choice in valid_programs?
    ↓
graph_nodes.py build() + ai_agent.py execute()
    ↓
ai_agent._run_single_cycle(): check_stop_conditions() → stop if matched
```

### Key Implementation Files

| File | Role |
|------|------|
| `agent/directive_extractor.py` | Extract directives from advice, `check_stop_conditions()` |
| `agent/directive_validator.py` | Apply program settings to intent (no stop logic) |
| `agent/workflow_engine.py` | `_apply_directives()` modifies valid_programs |
| `agent/graph_nodes.py` | Simple validation gate in `plan()` |
| `programs/ai_agent.py` | Post-execution stop check in `_run_single_cycle()` |
| `agent/session.py` | Store/retrieve directives |

## Usage Examples

### Example 1: Program-Specific Resolution

**User Advice:**
```
Use resolution 3.0 in autosol (anomalous signal is weak at high resolution)
but refine to 2.5 Angstroms.
```

**Extracted Directives:**
```json
{
  "program_settings": {
    "phenix.autosol": {"resolution": 3.0},
    "default": {"resolution": 2.5}
  }
}
```

**Behavior:**
- Cycle 1 (autosol): LLM sees directive, uses resolution=3.0 ✓
- Cycle 2 (refine): LLM forgets → Validator adds resolution=2.5 ✓
- Cycle 3 (refine): LLM uses 3.0 → Validator overrides to 2.5, logs warning ✓

### Directive Override Behavior (Attempt-Based)

When LLM and directive values differ, the system uses an attempt-based strategy:

| Attempt | Behavior | Rationale |
|---------|----------|-----------|
| 0 (first) | Use directive value | Honor user's explicit request |
| 1+ (retry) | Use LLM's value | LLM may be correcting syntax errors |

**Example:** User provides invalid selection syntax

```
User advice: selection=solvent molecule MES 88  (invalid Phenix syntax)

Attempt 0:
  - Directive says: selection=solvent molecule MES 88
  - LLM interprets: selection=resname MES and resseq 88
  - Decision: Use directive value (first attempt)
  - Result: FAILS (invalid syntax)

Attempt 1:
  - Directive says: selection=solvent molecule MES 88
  - LLM interprets: selection=resname MES and resseq 88
  - Decision: Use LLM value (retry - let it fix syntax)
  - Result: SUCCEEDS
```

This approach:
1. Respects user's explicit request first
2. Allows recovery from user syntax mistakes
3. Logs warnings so user knows what happened

### Example 2: Stop After First Refinement

**User Advice:**
```
Just run phaser and one round of refinement to check if MR worked.
Don't bother with validation.
```

**Extracted Directives:**
```json
{
  "stop_conditions": {
    "after_program": "phenix.refine",
    "max_refine_cycles": 1,
    "skip_validation": true
  }
}
```

**Behavior:**
- Cycle 1 (phaser): Runs normally
- Cycle 2 (refine): Runs normally
- Cycle 3: `max_refine_cycles=1` limit fires → refinement removed from
  valid programs. Because `skip_validation: true` is also set, the
  engine adds STOP only (no molprobity injected). Workflow ends.
- Without `skip_validation`, molprobity would be injected for a
  controlled landing before STOP.

### Example 3: Complex SAD Workflow

**User Advice:**
```
This is a SAD dataset with selenium. Use resolution 3.0 for autosol
(the anomalous signal is weak at high resolution) but refine to 2.0.
Stop after getting an initial model - I'll do ligand fitting manually.
Use anisotropic refinement for the final cycles.
```

**Extracted Directives:**
```json
{
  "program_settings": {
    "phenix.autosol": {"resolution": 3.0, "atom_type": "Se"},
    "phenix.refine": {"resolution": 2.0, "anisotropic_adp": true}
  },
  "stop_conditions": {
    "after_program": "phenix.autobuild",
    "skip_validation": true
  },
  "constraints": ["User will do ligand fitting manually"]
}
```

### Example 4: Workflow Preferences

**User Advice:**
```
Use experimental phasing, not molecular replacement.
Skip autobuild - I want to build the model manually.
```

**Extracted Directives:**
```json
{
  "workflow_preferences": {
    "use_experimental_phasing": true,
    "skip_programs": ["phenix.autobuild", "phenix.phaser"]
  }
}
```

**Behavior:**
- phenix.phaser removed from valid programs
- phenix.autobuild removed from valid programs
- Agent prefers phenix.autosol pathway
- **Phase advancement:** Skipped programs are treated as "done" for workflow
  phase detection. Without this, skipping a phase-gating program (like xtriage
  which gates the "analyze" phase) would deadlock the workflow — the phase
  would only offer the skipped program, leaving no valid programs.

### Example 5: MR-SAD Workflow

**User Advice:**
```
Run MR-SAD. Use the AlphaFold model for molecular replacement first,
then use the placed model for SAD phasing to locate Fe sites.
atom_type=Fe, sites=4, wavelength=1.1158
```

**Extracted Directives:**
```json
{
  "program_settings": {
    "phenix.autosol": {
      "atom_type": "Fe",
      "sites": 4,
      "wavelength": 1.1158
    }
  },
  "workflow_preferences": {
    "use_mr_sad": true,
    "use_experimental_phasing": true
  },
  "stop_conditions": {
    "after_program": "phenix.autosol",
    "skip_validation": true
  }
}
```

**Behavior:**
- In `obtain_model` phase: phaser prioritized, autosol removed (phaser must run first)
- After phaser: workflow enters `xray_mr_sad` state (experimental_phasing phase)
- autosol runs with `partpdb_file=PHASER.1.pdb` (auto-filled from phaser output)
- Workflow: xtriage → phaser → autosol (with partpdb_file) → STOP

## Testing

Run all directive-related tests:

```bash
cd improved_agent_v2

# Unit tests
python tests/tst_directive_extractor.py      # 73 tests
python tests/tst_directive_validator.py      # 53 tests
python tests/tst_session_directives.py       # 12 tests

# Integration tests
python tests/tst_directives_integration.py   # 16 tests

# All tests
python -m pytest tests/tst_directive*.py tests/tst_session_directives.py -v
```

## Rollback Plan

If issues arise, the directive system can be disabled by:

1. Setting `directives = {}` in session (disables all directive logic)
2. The validation step becomes a no-op when directives is empty
3. Original LLM-only behavior is preserved

## Future Enhancements

See [FUTURE_PLANS.md](../project/FUTURE_PLANS.md) for planned features:

1. **Interactive directive refinement** - Ask user to confirm extracted directives
2. **Directive learning** - Log which directives work well for future suggestions
3. **Directive templates** - Pre-built directive sets for common scenarios (SAD, MR, etc.)
4. **Directive conflict detection** - Warn when directives may conflict with each other

## Implemented Features

1. **Mid-session directive updates** - Users can now update advice when restarting the agent.
   The agent detects advice changes via hash comparison and re-extracts directives automatically.
   See [ARCHITECTURE.md](ARCHITECTURE.md#advice-change-detection) for details.

2. **Tutorial/Procedure Detection** - When the agent detects that user advice describes a
   specific limited procedure (like a tutorial), it automatically adds appropriate stop conditions.

### Example 5: Tutorial - Twinning Analysis

**README Content:**
```
Run Xtriage on the reflection data to analyze for twinning.
Check for evidence of twinning in the porin dataset.
```

**Extracted Directives:**
```json
{
  "stop_conditions": {
    "after_program": "phenix.xtriage",
    "skip_validation": true
  },
  "constraints": [
    "Run Xtriage on the reflection data to analyze for twinning.",
    "Check for evidence of twinning."
  ]
}
```

**Behavior:**
- Agent runs phenix.xtriage
- Agent stops after xtriage completes (doesn't continue to phasing/refinement)
- No molprobity validation required

### Automatic Stop Detection Patterns

The directive extractor recognizes these tutorial patterns and adds stop conditions:

| Pattern | Stop After |
|---------|------------|
| "run xtriage", "check for twinning", "analyze data quality" | phenix.xtriage |
| "run phaser", "test MR", "try molecular replacement" | phenix.phaser |
| "run mtriage", "analyze map" | phenix.mtriage |
| "run one refinement", "quick refinement test" | phenix.refine (1 cycle) |
| "run polder", "polder map", "omit map" | phenix.polder |
| "dock in map", "fit model to map" | phenix.dock_in_map |
| "map to model", "build model into map" | phenix.map_to_model |

### Example 6: Multi-Step Workflow (Polder then Refine)

**User Advice:**
```
Calculate a polder map for ligand MES 88 and then run refinement.
```

**Extracted Directives:**
```json
{
  "program_settings": {
    "phenix.polder": {
      "selection": "resname MES and resseq 88"
    }
  },
  "stop_conditions": {
    "start_with_program": "phenix.polder"
  }
}
```

**Behavior:**
- `start_with_program` is added to valid_programs list
- Cycle 1: Agent runs phenix.polder first
- Cycle 2+: Normal workflow continues (refinement, etc.)

**Note:** When user says "X then Y", the extractor sets:
- `start_with_program: X` (run this first)
- Does NOT set `after_program` (allows workflow to continue)

---

## Extending a Completed Workflow

Sometimes you want to run an additional program after the agent has already
declared the workflow complete (R-free at target, molprobity done). Rather
than starting over, you can resume the session with new advice.

### How It Works

When the agent resumes with advice whose hash differs from the stored value,
it sets an internal `advice_changed` flag. On the next cycle this flag:

1. **Steps the workflow back from `complete` to `validate` phase** (Q1 fix),
   giving the LLM access to post-refinement programs:
   `phenix.polder`, `phenix.molprobity`, `phenix.model_vs_data`,
   `phenix.map_correlations`
2. **Suppresses the AUTO-STOP check** that would otherwise terminate
   before the LLM even gets to plan

After one successful cycle the flag is cleared and the normal completion
logic resumes.

### Example 7: Polder Map on a Ligand After Completion

**Scenario:** A ligand-protein structure has been fully solved and validated.
You want an omit map to confirm the ligand density.

**Command:**
```bash
phenix.ai_agent \
    log_directory=AIAgent_42 \
    restart_mode=resume \
    project_advice="also run polder on the MES ligand in chain B residue 100"
```

**What happens:**
- Agent detects new advice hash → `advice_changed = True`
- PERCEIVE: steps back from `complete` to `validate` phase →
  `phenix.polder` appears in the program menu
- PLAN: AUTO-STOP suppressed for this cycle
- LLM sees new advice + polder in menu → selects `phenix.polder`
  with `selection=resname MES and resseq 100`
- Polder runs; omit map coefficients saved in session
- `advice_changed` cleared; next cycle would stop normally

**Note:** You do not need to specify `selection` in the advice — the LLM
will infer it. For better precision, include the PHENIX selection syntax:

```bash
project_advice="run polder selection='chain B and resname MES and resseq 100'"
```

### Example 8: Additional Validation After Completion

```bash
phenix.ai_agent \
    log_directory=AIAgent_42 \
    restart_mode=resume \
    project_advice="run model_vs_data to get a comprehensive validation report"
```

The LLM will select `phenix.model_vs_data` from the validate-phase menu.

### Programs Available for Follow-Up

The validate phase program menu available after step-back:

| Program | Purpose |
|---------|---------|
| `phenix.polder` | Ligand/residue omit maps (OMIT electron density) |
| `phenix.molprobity` | Geometry validation (re-run if you refined more) |
| `phenix.model_vs_data` | Comprehensive model-vs-data statistics |
| `phenix.map_correlations` | Map-model correlation analysis |

### One Step at a Time

Each resume with new advice triggers **one additional cycle**. For multiple
follow-up programs, run multiple resumes (each with different advice) or
use a multi-step directive:

```bash
# Two programs in sequence: polder then map_correlations
phenix.ai_agent \
    log_directory=AIAgent_42 \
    restart_mode=resume \
    project_advice="run polder on the MES ligand, then run map_correlations"
```

The directive extractor will recognise the multi-step request and set
`start_with_program: phenix.polder` with normal workflow continuation.
