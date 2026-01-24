# Unified Agent Transparency & Logging System

## Status: IMPLEMENTED ✓

The event system provides structured, transparent logging of the agent's decision-making process. Events flow from graph nodes through the API to the display layer, ensuring consistent output regardless of local or remote execution.

**Current version:** v40 (January 2025)

## Implementation Files

| File | Purpose |
|------|---------|
| `agent/event_log.py` | EventType constants, Verbosity levels, EventLog class |
| `agent/event_formatter.py` | EventFormatter for human-readable output |
| `agent/graph_nodes.py` | Event emission via `_emit()` helper |
| `phenix_ai/run_ai_agent.py` | Events in API response |
| `programs/ai_agent.py` | EventFormatter integration, verbosity parameter |

## Goals

1. **Transparency** - Show what decisions are made and why at each step
2. **Consistency** - Same output format for local and remote execution
3. **Controllability** - Support quiet, normal, verbose, and debug verbosity levels
4. **Full Reasoning** - Show complete LLM reasoning without truncation
5. **User Feedback** - Prominently warn when user requests can't be fulfilled

---

## Architecture

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

## Event Types

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
| `debug` | debug | Debug trace | `message` |

---

## Verbosity Levels

```python
class Verbosity:
    QUIET = "quiet"    # Errors, warnings, and cycle summaries
    NORMAL = "normal"  # Key decisions and metrics (default)
    VERBOSE = "verbose"  # Includes file selection details
    DEBUG = "debug"    # Full internal state and traces
```

**PHIL Parameter:**
```phil
verbosity = normal
  .type = choice(quiet, normal, verbose, debug)
  .help = Controls how much detail is shown in agent output.
```

---

## Example Output

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

## Implementation Files

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

## Usage

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

## File Selection Tracking

The CommandBuilder records WHY each file was selected:

| Reason | Description |
|--------|-------------|
| `best_files` | From BestFilesTracker (highest scored) |
| `rfree_locked` | Locked R-free MTZ from first refinement |
| `llm_selected` | LLM explicitly chose this file |
| `best_files_override` | LLM choice overridden by best_files |
| `auto_selected` | Fallback based on category/extension |

---

## Testing

```bash
# Run event system tests
python3 tests/test_event_system.py

# Test coverage includes:
# - EventType and Verbosity constants
# - EventLog emit and filtering
# - EventFormatter at all verbosity levels
# - USER_REQUEST_INVALID formatting
# - Full reasoning preservation (no truncation)
```

---

## Change History

| Version | Date | Changes |
|---------|------|---------|
| v40 | 2025-01 | Added USER_REQUEST_INVALID event type |
| v38 | 2025-01 | Phase 4: Display integration, verbosity parameter |
| v37 | 2025-01 | Phase 3: Transport integration (API schema, response building) |
| v36 | 2025-01 | Phase 2: Instrumented graph_nodes.py and command_builder.py |
| v34 | 2025-01 | Phase 1: Created EventLog, EventFormatter, 17 event types |
