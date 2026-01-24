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
    "use_molecular_replacement": false
  },

  "constraints": [
    "Do not add waters until R-free < 0.30",
    "Use TLS refinement after cycle 3"
  ]
}
```

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
| `tests/test_directive_extractor.py` | 31 unit tests for extraction |
| `tests/test_directive_validator.py` | 26 unit tests for validation |
| `tests/test_session_directives.py` | 12 unit tests for session methods |
| `tests/test_directives_integration.py` | 16 integration tests |
| `tests/test_decision_flow.py` | Decision flow architecture tests |

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
- Cycle 3: Validator detects stop condition → STOP (no molprobity required)

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

## Testing

Run all directive-related tests:

```bash
cd improved_agent_v2

# Unit tests
python tests/test_directive_extractor.py      # 31 tests
python tests/test_directive_validator.py      # 26 tests
python tests/test_session_directives.py       # 12 tests

# Integration tests
python tests/test_directives_integration.py   # 16 tests

# All tests
python -m pytest tests/test_directive*.py tests/test_session_directives.py -v
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
