# PHENIX AI Agent Architecture

## System Overview

The PHENIX AI Agent is an automated crystallographic workflow system that:
1. Analyzes experimental data and logs
2. Decides which program to run next
3. Executes the program
4. Tracks results and quality metrics
5. Repeats until structure is solved or workflow completes

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                              CLIENT                                      │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                         ai_agent.py                                 ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │   Session   │  │    Agent    │  │   Command   │  │    Log     │ ││
│  │  │   Tracker   │  │  Interface  │  │  Executor   │  │   Parser   │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  │        │                │                │                │        ││
│  │        ▼                ▼                ▼                ▼        ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │ BestFiles   │  │ Local/Remote│  │  easy_run   │  │ Metrics    │ ││
│  │  │ Tracker     │  │   Agent     │  │  subprocess │  │ Extraction │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  └─────────────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    │ v2 JSON Request
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         DECISION ENGINE                                  │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                       run_ai_agent.py                               ││
│  │  ┌─────────────────────────────────────────────────────────────────┐││
│  │  │                      LangGraph                                  │││
│  │  │  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌───────┐│││
│  │  │  │PERCEIVE │─▶│  PLAN   │─▶│  BUILD  │─▶│VALIDATE │─▶│OUTPUT ││││
│  │  │  └─────────┘  └────┬────┘  └─────────┘  └────┬────┘  └───────┘│││
│  │  │                    ▲                          │                │││
│  │  │                    └──────── retry < 3 ───────┤                │││
│  │  │                                               ▼                │││
│  │  │                                          ┌──────────┐          │││
│  │  │                                          │ FALLBACK │──▶OUTPUT │││
│  │  │                                          └──────────┘          │││
│  │  └─────────────────────────────────────────────────────────────────┘││
│  └─────────────────────────────────────────────────────────────────────┘│
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────────────┐│
│  │                        Knowledge Layer                              ││
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ ││
│  │  │   YAML      │  │   Rules     │  │  Templates  │  │  Workflow  │ ││
│  │  │  Programs   │  │  Selector   │  │   Builder   │  │   State    │ ││
│  │  └─────────────┘  └─────────────┘  └─────────────┘  └────────────┘ ││
│  └─────────────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────────┘
```

## Unified Request Flow

Both LocalAgent and RemoteAgent use the same v2 JSON API and **identical transport encoding**:

```
┌─────────────────────────────────────────────────────────────────┐
│                         Client Side                              │
├─────────────────────────────────────────────────────────────────┤
│  LocalAgent / RemoteAgent                                        │
│       │                                                          │
│       ▼                                                          │
│  build_request_v2()                                              │
│       │                                                          │
│       ▼                                                          │
│  prepare_request_for_transport()                                 │
│       │                                                          │
│       ├── 1. sanitize_request() [YAML-driven]                    │
│       │       ├── Remove ZZxxZZ markers                          │
│       │       ├── Truncate long quoted strings                   │
│       │       ├── Replace tabs with spaces                       │
│       │       └── Remove control characters                      │
│       │                                                          │
│       ├── 2. json.dumps()                                        │
│       │                                                          │
│       └── 3. encode_for_rest() [text_as_simple_string]           │
│                    │                                             │
│     ┌──────────────┴──────────────┐                              │
│     │                             │                              │
│  LocalAgent:                 RemoteAgent:                        │
│  process_request_from_      send encoded to                      │
│  transport() locally        REST server                          │
│     │                             │                              │
│     └──────────────┬──────────────┘                              │
└────────────────────│────────────────────────────────────────────┘
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                      run_ai_agent.py                             │
├─────────────────────────────────────────────────────────────────┤
│  process_request_from_transport()                                │
│       │                                                          │
│       ├── 1. decode_from_rest() [simple_string_as_text]          │
│       └── 2. json.loads()                                        │
│                    │                                             │
│                    ▼                                             │
│              LangGraph Execution                                 │
│                    │                                             │
│                    ▼                                             │
│  prepare_response_for_transport()                                │
│       │                                                          │
│       ├── 1. sanitize_response()                                 │
│       ├── 2. json.dumps()                                        │
│       └── 3. encode_for_rest()                                   │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Decision**: LocalAgent performs the full encode/decode roundtrip locally, even though it could skip encoding. This ensures that transport bugs are caught during local testing rather than only appearing in production server scenarios.

## Component Responsibilities

### Client Components

#### ai_agent.py
Main entry point and orchestrator:
- Manages the iterative workflow loop
- Coordinates all other components
- Handles command execution
- Maintains session state

#### Session Tracker (session.py)
Persists workflow state across cycles:
- Experiment type (X-ray/cryo-EM)
- Resolution
- R-free MTZ path
- Cycle history
- **User Directives** (extracted structured instructions)

#### Directive System
Extracts and enforces user instructions:
- `directive_extractor.py`: Parses natural language → structured JSON
- `directive_validator.py`: Validates LLM decisions against directives
- Integrated into graph nodes for consistent enforcement
- See [USER_DIRECTIVES.md](../guides/USER_DIRECTIVES.md) for details

#### BestFilesTracker (best_files_tracker.py)
Tracks highest-quality files:
- Scores files using YAML-based criteria
- Maintains current best by category
- Provides files for server decisions

#### Agent Interface (LocalAgent/RemoteAgent)
Both agents use identical interface, v2 JSON format, **and transport encoding**:
- `LocalAgent`: Full encode/decode roundtrip, then calls `run_ai_agent.run()`
- `RemoteAgent`: Encodes and sends to REST server
- Same `prepare_request_for_transport()` for encoding
- Same `history_record` response format
- Transport module: `agent/transport.py`
- Configuration: `knowledge/transport.yaml`

#### Log Parsers (phenix_ai/log_parsers.py)
Extracts metrics and output files from program log output:

**Key functions:**
- `extract_all_metrics(log_text, program)`: Main entry point for metric extraction
- `detect_program(log_text)`: Identifies which program generated a log
- `extract_output_files(log_text, working_dir)`: Finds output file paths
- `format_metrics_report(metrics)`: Generates "FINAL QUALITY METRICS REPORT"

**Note**: Metric extraction patterns are defined in `programs.yaml` and loaded via `metric_patterns.py`. Hardcoded extractors in log_parsers.py only handle complex cases (tables, multi-line context) that YAML patterns can't express.

**YAML `log_parsing` pattern conventions** (see `programs.yaml`):

| Field | Values | Usage |
|-------|--------|-------|
| `extract` | `first` (default), `last` | `last` required when a program emits one metric line per cycle (RSR, autobuild, phaser) |
| `pick_min` | `true` | After collecting all matches, take the minimum value. Used for xtriage `resolution` (multiple lines; smallest value = highest resolution) |
| `type` | `float`, `int`, `string`, `boolean` | Controls return type |

**Important audited patterns** (fixes applied in v112.14):
- `phenix.xtriage` `resolution`: uses `pick_min: true`. Pattern anchored with `^\s*` and negative lookbehind to prevent "Completeness in resolution range: 1" from matching. Skip group `(?:[0-9.]+\s*[-]\s*)?` handles "50.00 - 2.30" format.
- `phenix.real_space_refine` `map_cc`: uses `extract: last`. RSR emits one CC_mask line per macro-cycle; `last` gives the final value.

**Adding new program support:**
See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide.

### Centralized Configuration Modules

The agent uses a centralized YAML-driven architecture that reduces program configuration from 7 files to 2-3 files.

#### Metric Patterns (knowledge/metric_patterns.py)
Centralizes metric extraction patterns from programs.yaml:

```
programs.yaml          metric_patterns.py
     │                        │
     ▼                        ▼
log_parsing: ────────▶ extract_metrics_for_program()
  pattern:                    │
  type:              ┌────────┴────────┐
  display_name:      ▼                 ▼
                log_parsers.py    session.py
```

**Key functions:**
- `get_all_metric_patterns()`: Load all patterns from YAML (cached)
- `extract_metrics_for_program(log_text, program)`: Extract metrics using YAML patterns
- `format_metric_value(program, metric, value)`: Format using YAML config

#### Program Registration (knowledge/program_registration.py)
Reads `done_tracking` blocks from programs.yaml for workflow done flags:

```
programs.yaml             program_registration.py
     │                            │
     ▼                            ▼
done_tracking:           get_program_done_flag_map()  ← ALL programs
  flag: "xxx_done"       get_trackable_programs()     ← strategy: "run_once"
  strategy: "run_once"            │
  # or "count"           ┌────────┴────────┐
  # or "set_flag"        ▼                 ▼
                 workflow_state.py   workflow_engine.py
                 (_set_done_flags)   (context building)
```

**Key functions:**
- `get_program_done_flag_map()`: All programs → done flag names (from YAML)
- `get_trackable_programs()`: Programs with `strategy: "run_once"` (filtered from valid list after completion)

#### Summary Display (knowledge/summary_display.py)
Configures session summary display from metrics.yaml:

```
metrics.yaml           summary_display.py
     │                        │
     ▼                        ▼
summary_display: ─────▶ format_quality_table_rows()
  quality_table:              │
  step_metrics:        ┌──────┴──────┐
                       ▼             ▼
                session.py     session.py
             (Final Quality)  (Steps table)
```

**Key functions:**
- `get_quality_table_config()`: Load quality table row configs
- `get_step_metrics_config()`: Load per-program step metric configs
- `format_quality_table_rows()`: Format Final Quality table
- `format_step_metric()`: Format a step's key metric

#### Program Validator (agent/program_validator.py)
Validates that programs are fully configured:
```bash
python agent/program_validator.py phenix.map_symmetry  # Check one program
python agent/program_validator.py --all                 # Check all programs
```

### Decision Engine Components

#### run_ai_agent.py
Single entry point for all requests:
- Parses v2 JSON requests
- Validates request schema
- Creates graph initial state
- Invokes LangGraph
- Returns v2 JSON response

#### LangGraph Pipeline
Decision-making workflow:

```
PERCEIVE ──┬──▶ PLAN ─▶ BUILD ─▶ VALIDATE ──┬──▶ OUTPUT ─▶ END
           │     │        │         │    │    │
           │     │        │         │    │    │
        Analyze Select   Build   Validate │  Format
        inputs  program  command response │  decision
           │                        │     │
           │               retry <3 │     │ retry >=3
           │                   ┌────┘     │
           │                   ▼          ▼
           │                  PLAN     FALLBACK ─▶ OUTPUT ─▶ END
           │
           └──▶ OUTPUT ─▶ END  (if red flag abort)
```

#### Decision Flow Architecture

The decision flow follows a clean, layered architecture where each component has a single responsibility:

```
┌─────────────────────────────────────────────────────────────────┐
│                     SESSION START                                │
│  directive_extractor.py extracts structured directives          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 1: WORKFLOW ENGINE                         │
│  workflow_engine.py determines valid_programs                    │
│                                                                  │
│  Input: phase, context, directives                              │
│  Output: valid_programs list (includes STOP if appropriate)     │
│                                                                  │
│  Responsibilities:                                               │
│  - Get programs from workflow phase definition                   │
│  - Add directive-required programs (after_program target)        │
│  - Add start_with_program to front of list (multi-step flows)   │
│  - Add STOP if skip_validation=true                             │
│  - Apply workflow preferences (skip/prefer)                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 2: LLM DECISION (PLAN)                     │
│  graph_nodes.py plan() function                                  │
│                                                                  │
│  Input: valid_programs, context, history                        │
│  Output: chosen program (must be in valid_programs)             │
│                                                                  │
│  Responsibilities:                                               │
│  - LLM selects best program from valid_programs                 │
│  - Simple validation: is choice in valid_programs?              │
│  - If invalid choice, use first valid program                   │
│                                                                  │
│  NOT responsible for:                                            │
│  - Directive logic (already handled in Layer 1)                 │
│  - Stop condition evaluation (handled in Layer 4)               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 3: BUILD & EXECUTE                         │
│  graph_nodes.py build(), ai_agent.py execute                    │
│                                                                  │
│  Input: chosen program, files, context                          │
│  Output: command result, output files                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 LAYER 4: POST-EXECUTION CHECK                    │
│  ai_agent.py _run_single_cycle() - SINGLE PLACE                 │
│                                                                  │
│  Input: completed program, cycle number, directives             │
│  Output: should_stop (bool), stop_reason                        │
│                                                                  │
│  Responsibilities:                                               │
│  - Check if after_program condition met                         │
│  - Check if after_cycle condition met                           │
│  - Check if metric targets met (r_free, map_cc)                 │
│  - This is the ONLY place stop conditions are evaluated         │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Principles:**

1. **Single Responsibility**: Each layer handles one aspect of decision-making
2. **Clear Data Flow**: Directives flow from extraction → workflow engine → validated programs
3. **Post-Execution Stop**: Stop conditions are only checked after program execution, not during planning
4. **Simple Validation Gate**: The plan() function validates that LLM choice is in valid_programs list

#### Knowledge Layer

**YAML Programs** (`knowledge/programs.yaml`):
- Program definitions
- Input/output specifications
- Invariants (rules that must be satisfied)
- Scoring rules

**Rules Selector** (`rules_selector.py`):
- Matches workflow state to applicable programs
- Filters by experiment type
- Applies priority rules

**Command Builder** (`command_builder.py`):
- **Single entry point** for all command generation
- 4-stage pipeline:
  1. `_select_files()` - Priority-based file selection
  2. `_build_strategy()` - Strategy with output_prefix
  3. `_apply_invariants()` - Auto-fill resolution, R-free flags
  4. `_assemble_command()` - Final command string
- Uses `CommandContext` dataclass for all parameters
- Replaces fragmented logic from template_builder and graph_nodes

**Template Builder** (`template_builder.py`):
- Legacy interface (delegates to CommandBuilder)
- YAML-to-slot mapping
- Invariant checking

**Workflow State** (`workflow_state.py`):
- Determines current workflow position
- Categorizes available files
- Detects experiment type

## LLM Provider Architecture

The agent supports three LLM providers: **google**, **openai**, and **ollama**. The
provider is set via `params.communication.provider` (or the `LLM_PROVIDER` env var)
and flows consistently through every LLM call point — there is no cross-provider
fallback.

### LLM Call Points

The system uses three distinct LLM roles, each with provider-specific model defaults:

| Activity | Function | ollama | google | openai |
|---|---|---|---|---|
| **Agent decisions** (PLAN node) | `get_planning_llm()` → `get_expensive_llm()` | qwen3:32b (json_mode) | gemini-2.5-pro | gpt-5 |
| **Log summarization** | `get_cheap_llm()` | qwen2.5:7b | gemini-2.5-flash-lite | gpt-5-nano |
| **RAG analysis** | `get_expensive_llm()` | qwen3:32b | gemini-2.5-pro | gpt-5 |
| **Directive extraction** | `call_llm_simple()` | direct ollama HTTP | gemini-2.0-flash | (provider default) |

### Provider Flow

```
params.communication.provider = "<provider>"
         │
         ├─► LocalAgent.decide_next_step()
         │       └─► request["settings"]["provider"] = "<provider>"
         │            └─► run_ai_agent.py: create_initial_state(provider=...)
         │                 └─► graph_nodes.py PLAN: get_planning_llm(provider)
         │                      └─► get_expensive_llm(provider)
         │
         ├─► _extract_directives()
         │       ├─► ollama: forces run_on_server=False (local ollama not on server)
         │       └─► run_directive_extraction(provider=...)
         │            └─► call_llm_simple(provider=...)
         │
         └─► run_ai_analysis.run() [summarization + RAG]
                 └─► setup_llms(provider=...)
                      ├─► expensive_llm = get_expensive_llm(provider)
                      ├─► cheap_llm = get_cheap_llm(provider)
                      └─► embeddings = provider-specific embeddings
```

### Key Files

| File | Role |
|------|------|
| `core/llm.py` | `get_llm_and_embeddings()`, `get_expensive_llm()`, `get_cheap_llm()` — model creation |
| `agent/graph_nodes.py` | `get_planning_llm()` — cached LLM for PLAN node decisions |
| `agent/api_client.py` | `call_llm_simple()` — lightweight direct calls (directive extraction) |
| `utils/run_utils.py` | `setup_llms()` — creates all three LLM roles for analysis |
| `phenix_ai/local_agent.py` | Reads provider from `params.communication.provider`, passes to request |
| `phenix_ai/run_ai_agent.py` | Reads provider from settings, passes to `create_initial_state()` |

### Ollama-Specific Behavior

- **Directive extraction** forces `run_on_server=False` because the Phenix REST
  server doesn't have access to the user's local ollama service
- **Planning LLM** uses `json_mode=True` (sets `format="json"` in ChatOllama)
  because ollama models need explicit JSON formatting; Google and OpenAI handle
  structured output without this flag
- **GPU requirement**: `rest_server.requires_gpu = True` is set for ollama to
  ensure server dispatch targets GPU-capable machines

### Adding a New Provider

1. Add provider case to `get_llm_and_embeddings()` in `core/llm.py`
2. Add validation check in `validate_provider()` in `agent/graph_nodes.py`
3. Add `_call_<provider>_llm()` function in `agent/api_client.py`
4. Add rate limit handler in `agent/rate_limit_handler.py` (optional)

## Data Flow

### Request Building (Client)

```python
# Both agents build requests the same way
session_state = {
    "resolution": 2.5,
    "experiment_type": "xray",
    "rfree_mtz": "/path/to/data.mtz",
    "best_files": {"model": "/path/to/model.pdb"}
}

request = build_request_v2(
    files=available_files,
    cycle_number=5,
    log_content=last_log,
    history=cycle_history,
    session_state=session_state,
    user_advice=guidelines,
    ...\n)
```

**`best_files` value type contract:** Most entries are plain strings (`category → path`). However, multi-file entries such as `half_map` may be stored as a list:

```python
"best_files": {
    "model":    "/path/to/model.pdb",       # single file → str
    "half_map": ["/path/map1.mrc", "/path/map2.mrc"],  # two files → list
}
```

All server-side code that reads `best_files` values must handle both types. Use `CommandBuilder._best_path(value)` to safely extract a single path from either a string or a list. Do **not** pass `best_files` values directly to `os.path` functions — this caused a cycle=2 crash (`TypeError: expected str, bytes or os.PathLike, not list`) that only surfaced when `session_state` was re-sent from a client that stored `half_map` as a list.

**Programs that need both `full_map` and `half_map` simultaneously:** Several programs (`phenix.mtriage`, `phenix.predict_and_build`, `phenix.map_to_model`) accept both a full map AND half maps at once — the half maps are not redundant, they enable FSC-based resolution calculation or density modification. These programs are marked `keep_half_maps_with_full_map: true` in `programs.yaml`. The post-selection validation in `CommandBuilder._select_files()` checks this flag before removing half maps. When adding a new program that needs both simultaneously, add this flag to its YAML entry.

### Request Processing (Server)

```python
# run_ai_agent.py processes all requests uniformly
def run(request_json):
    request = _parse_request(request_json)
    response = _process_request(request)
    return _build_group_args_response(response)
```

### Response Handling (Client)

```python
# Both agents receive the same response format
result = agent.decide_next_step(...)
history_record = result  # Contains program, command, reasoning, etc.
command = history_record["command"]
```

## File Organization

See the [README.md](../README.md#directory-structure) for the complete directory tree.

Key directories:
- `agent/` — Core agent logic (session, graph nodes, command builder, workflow engine, etc.)
- `knowledge/` — YAML configuration files and supporting Python modules
- `phenix_ai/` — Runtime entry points (local/remote agent, log parsers)
- `programs/` — PHENIX program integration (main entry point `ai_agent.py`)
- `analysis/` — Post-run log analysis and session evaluation
- `core/` — LLM provider abstraction
- `validation/` — Command validation framework
- `tests/` — 34 test files with 747+ tests

## Key Design Decisions

### 1. Unified v2 JSON API

**Rationale**: Consistency between local and remote execution.

**Implementation**:
- Both agents use `build_request_v2()` to create requests
- Single `run_ai_agent.run(request_json)` entry point
- Identical response format for both paths

### 2. Client-Server Separation

**Rationale**: Allows server-side updates without client reinstallation.

**Implementation**:
- Client only gathers inputs and executes commands
- Server makes all decisions
- Clean API boundary via JSON protocol

### 3. Configuration Over Code

**Rationale**: Domain experts can modify behavior without code changes.

**Implementation**:
- Program definitions in YAML
- Scoring rules in YAML
- Workflow transitions in YAML

### 4. Session State Tracking

**Rationale**: Ensures consistency across workflow cycles.

**Implementation**:
- R-free MTZ locked after first refinement
- Experiment type locked after detection
- Best files tracked with quality scores
- All state transmitted in `session_state`

### 5. LangGraph Pipeline

**Rationale**: Modular, testable decision-making.

**Implementation**:
- Discrete nodes (PERCEIVE, PLAN, BUILD, VALIDATE, OUTPUT)
- State passed between nodes
- Easy to add/modify nodes

### 6. Automation Modes

**Rationale**: Different users need different levels of control over the workflow.

**Implementation**:
Two modes controlled by `maximum_automation` parameter:

| Mode | Setting | `predict_and_build` Behavior |
|------|---------|------------------------------|
| **Automated** | `maximum_automation=True` (default) | Runs complete workflow (prediction → MR → building) |
| **Stepwise** | `maximum_automation=False` | Stops after prediction only |

The `automation_path` is set in workflow state and propagated to all decision points:

```python
# workflow_engine.py
context["automation_path"] = "automated" if maximum_automation else "stepwise"

# workflow_engine.py (_check_program_prerequisites)
# Block after full workflow has completed (both automated and stepwise)
if context.get("predict_full_done"):
    return False
# Block in stepwise mode after prediction-only step
if automation_path == "stepwise" and context.get("predict_done"):
    return False
```

**Stepwise workflow example**:
```
xtriage → predict_and_build(stop_after_predict) → process_predicted_model → phaser → refine
```

### 7. Fallback Program Tracking

**Rationale**: When LLM planning fails 3 times, the fallback mechanism must correctly report which program it selected.

**Implementation**:
- Fallback node now sets `state["program"]` when returning a command
- Response builder uses `state["program"]` if available, falls back to `intent["program"]`
- Prevents mismatch where PLAN shows one program but command is different

## Workflow States

### X-ray Crystallography Workflow

```
xray_initial → xtriage → xray_analyzed
                              │
              ┌───────────────┼───────────────┐
              ↓               ↓               ↓
       predict_and_build    phaser         autosol†
              ↓               │               ↓
    xray_has_prediction       │        xray_has_phases
              ↓               │               ↓
  process_predicted_model     │          autobuild
              ↓               │               │
    xray_model_processed ─────┘               │
              ↓                               │
           phaser ────────────────────────────┘
              │
              ├── [anomalous data] → xray_mr_sad
              │                         ↓
              │                 autosol (partpdb_file=PHASER.pdb)
              │                         ↓
              │                  xray_has_phases
              │                         ↓
              │                    autobuild ──┐
              ↓                               │
       xray_has_model                         │
              ↓                               │
         refine (loop) ←──────────────────────┘
              ↓
       xray_refined
         ↓    ↓    ↓
    molprobity refine STOP
              │
              ↓
       [if ligand] ligandfit → pdbtools → refine
```

† Prioritized when strong anomalous signal detected (measurability > 0.10)
‡ MR-SAD: phaser places model first, then autosol uses it as partpdb_file

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `xray_initial` | Starting point | xtriage |
| `xray_analyzed` | After data analysis | predict_and_build, phaser, autosol |
| `xray_has_prediction` | Have AlphaFold model | process_predicted_model |
| `xray_mr_sad` | After phaser + anomalous data (MR-SAD) | autosol (with partpdb_file) |
| `xray_has_phases` | After experimental phasing | autobuild |
| `xray_has_model` | Have placed model | refine |
| `xray_refined` | After refinement **or** validation (`refine` and `validate` phases share this external state; internal `phase_info["phase"]` distinguishes them) | refine, molprobity, autobuild, STOP |

### Cryo-EM Workflow

```
cryoem_initial → mtriage → cryoem_analyzed
                                │
                ┌───────────────┼───────────────┐
                ↓               ↓               ↓
      predict_and_build    dock_in_map    [half-maps only]
                │               │               ↓
                └───────┬───────┘         resolve_cryo_em
                        ↓                       ↓
                cryoem_has_model          map_sharpening
                        ↓                       │
              real_space_refine (loop) ←────────┘
                        ↓
                 cryoem_refined
                   ↓    ↓    ↓
              molprobity RSR  STOP
```

| State | Description | Valid Programs |
|-------|-------------|----------------|
| `cryoem_initial` | Starting point | mtriage |
| `cryoem_analyzed` | After map analysis | predict_and_build, dock_in_map |
| `cryoem_has_model` | Half-map optimisation / model check (legacy name — no model yet; `check_map` and `optimize_map` phases) | resolve_cryo_em, map_sharpening |
| `cryoem_docked` | Model docked, ready for first real-space refinement (`ready_to_refine` phase) | real_space_refine |
| `cryoem_refined` | After refinement **or** validation (`refine` and `validate` phases share this external state; internal `phase_info["phase"]` distinguishes them) | real_space_refine, molprobity, STOP |

## Workflow History and Done Flags

`workflow_state.py` maintains two complementary systems for tracking
completed programs:

### _analyze_history: done flag extraction

`_analyze_history(history)` reads the history list and sets done flags
(`refine_done`, `resolve_cryo_em_done`, …) based on marker strings in
each record's `result` field. Key invariants:

- **Failure blocks flags**: `_is_failed_result(result)` is checked
  first. If the result is a failure, all done-flag extraction for that
  record is skipped. Failure signals (in priority order):
  1. Shell exit code (caught at the subprocess layer before this function)
  2. Specific Phenix terminal phrases: `FAILED`, `SORRY:`, `SORRY `,
     `*** ERROR`, `FATAL:`, `TRACEBACK`, `EXCEPTION`
  - Generic `ERROR` words without those prefixes are **not** failure
    signals — Phenix logs contain "Error model parameter", "Expected
    errors: 0", etc. and these must not suppress done flags.

- **predict_and_build cascade**: `predict_full_done` additionally sets
  `refine_done=True` and increments `refine_count` (Python-only; no
  YAML `history_detection` entry) because a full `predict_and_build`
  run includes refinement internally.

### _clear_zombie_done_flags: crash recovery

After `_analyze_history`, `detect_workflow_state` calls
`_clear_zombie_done_flags(history_info, available_files)`.

A **zombie state** occurs when the agent crashed mid-cycle or the user
deleted output files: history records `done_flag=True`, but the output
file is missing from disk. The phase detector sees `done=True` and
skips the program; the workflow becomes stuck.

The function checks each crashable done flag against its expected
output filename pattern:

| Done flag | Expected output | Also clears | Any-PDB fallback |
|-----------|----------------|-------------|-----------------|
| `resolve_cryo_em_done` | `denmod*.ccp4/mrc` | `has_full_map` | No |
| `predict_full_done` | `*_overall_best.pdb` | `has_placed_model` | Yes |
| `dock_done` | `*docked*.pdb`, `*dock*map*.pdb`, `placed_model*.pdb`, `*_placed*.pdb` | `has_placed_model` | Yes |
| `refine_done` | `*_refine_NNN.pdb` | decrements `refine_count` | Yes |
| `rsr_done` | `*_real_space_refined*.pdb` | decrements `rsr_count` | Yes |

**Any-PDB fallback**: when `True`, the presence of any `.pdb` file on disk suppresses
zombie detection even if the pattern doesn't match. Used for model-producing programs
where output may be renamed. `resolve_cryo_em_done` uses `False` because `denmod_map.ccp4`
is referenced by name downstream.

Flags are cleared **in-memory only** — history is never rewritten.
The function returns a `zombie_diagnostics` list which `detect_workflow_state`
attaches to the state dict under the key `"zombie_diagnostics"`. The PERCEIVE
node logs each entry prefixed `"PERCEIVE: ZOMBIE STATE — "`. This makes
crash/restart re-runs self-explaining: the user can see exactly which done flag
was cleared and why the program appears again.

## Best Files Tracking

The agent maintains a **Best Files Tracker** that identifies and tracks the highest-quality file of each type throughout a session, ensuring programs always receive optimal inputs.

### Categories Tracked

| Category | Description | Example |
|----------|-------------|---------|
| `model` | Best atomic model (PDB/mmCIF) | `refine_002_001.pdb` |
| `map` | Best full cryo-EM map | `denmod_map.ccp4` |
| `mtz` | Best reflection data | `data_with_rfree.mtz` |
| `map_coefficients` | Best map coefficients | `refine_001_001.mtz` |
| `sequence` | Sequence file | `sequence.fa` |
| `ligand_cif` | Ligand restraints | `LIG.cif` |

### Scoring System

Each file is scored based on **processing stage** (0-100 points) and **quality metrics** (0-100 points). All scoring parameters are configurable in `knowledge/metrics.yaml` under `best_files_scoring`.

**Model stage scores:** refined=100, rsr_output=100, autobuild_output=80, docked=60, processed_predicted=50, predicted=40, phaser_output=30, pdb=10

**Model metric scores:** R-free (40 pts, linear_inverse), Map CC (30 pts, linear), Clashscore (30 pts, linear_inverse)

**Map stage scores:** optimized_full_map=100, sharpened=90, density_modified=80, full_map=50, half_map=10. Plus resolution bonus (0-30 pts).

### R-free MTZ Locking (X-ray Only)

Once an MTZ with R-free flags is identified, it is **locked for the entire session**. This ensures consistent cross-validation statistics throughout refinement. The locked MTZ has absolute priority (Priority 0) over all other file selection mechanisms.

### File Selection Priority

When building commands, files are selected in this order:

0. **Locked R-free MTZ** — For MTZ inputs in X-ray refinement (highest priority)
1. **Best Files** — From the BestFilesTracker scoring system
2. **Categorized Files** — From workflow_state file categorization
3. **Extension Matching** — Search available_files by extension

### Key Files

| File | Role |
|------|------|
| `agent/best_files_tracker.py` | Core tracker class with scoring |
| `agent/session.py` | Integration with session persistence |
| `agent/template_builder.py` | Uses best_files for command building |

### Companion File Discovery

Some clients only track a subset of program output files. The agent discovers
missing companion files in two layers:

**Layer 1: `graph_nodes._discover_companion_files()`** — Runs in the perceive
node before file categorization. Triggered by file patterns in available_files:

| Trigger Pattern | Companions Discovered |
|----------------|-----------------------|
| `refine_NNN_data.mtz` | `refine_NNN.mtz` (map coefficients), `refine_NNN.pdb` (model) |
| `overall_best_*.mtz` | `overall_best.pdb` (autobuild model) |
| Files in `sub_NN_*/` dirs | Scans `sub_*_pdbtools/` for `*_with_ligand.pdb` |

All discovered files are checked with `os.path.exists()` and deduplicated.

**Layer 2: `session._find_missing_outputs()`** — Runs in `get_available_files()`
to supplement cycle output_files from session data.

### Intermediate File Filtering

`graph_nodes._filter_intermediate_files()` removes temporary/intermediate files
before categorization. Runs after companion discovery in the perceive node.

**Filtered patterns:**
- Files in `/TEMP`, `/temp`, `/TEMP0/`, `/scratch/` directories
- Files with `EDITED_` or `TEMP_` prefixes

These are internal working files from programs like ligandfit that should never
be used as inputs to other programs.

### Best Files Exclusion Check

When `best_files["model"]` is applied as a model override for refinement
programs, it is first checked against the program's `exclude_categories`. If
the best model is in an excluded category (e.g., `ligand_fit_1.pdb` in
`ligand_fit_output` → parent `ligand`), it is skipped and category-based
selection runs instead. This prevents ligand fragments from being used as the
protein model for refinement.

Applied in `command_builder.py` at two points: (1) pre-population of the model
slot from best_files, and (2) LLM override where best_files would normally
take precedence over the LLM's model choice.
| `knowledge/metrics.yaml` | Scoring configuration (best_files_scoring section) |

## Metrics and Stop Conditions

### Quality Assessment

Metrics are evaluated against YAML-defined thresholds:

| Quality | R-free (2.0Å) | Map CC | Anomalous Measurability |
|---------|---------------|--------|-------------------------|
| Good | < 0.25 | > 0.80 | > 0.10 |
| Acceptable | < 0.30 | > 0.70 | > 0.05 |
| Poor | ≥ 0.30 | ≤ 0.70 | ≤ 0.05 |

### Stop Conditions

The workflow stops when:

1. **Success**: Quality metrics reach targets AND validation done
2. **Plateau**: < 0.5% improvement for 2+ consecutive cycles
3. **Hopeless bailout**: R-free > 0.50 after 1+ refinement cycle
4. **Hard limit**: 3+ refinement cycles regardless of R-free
5. **Validation Gate**: Must run molprobity before stopping if R-free is good

### Refinement Loop Enforcement

When `_is_at_target()` returns True (conditions 2-4 above),
`get_valid_programs()` actively **removes** `phenix.refine` and
`phenix.real_space_refine` from valid programs in both `validate` and
`refine` phases, and adds `STOP`. This prevents the LLM from selecting
refinement even when it appears as a phase-preferred program.

**Exception:** `needs_post_ligandfit_refine` always allows refinement.
After ligand fitting changes the model, re-refinement is scientifically
required regardless of the current cycle count or R-free value.

The validation gate prevents stopping without validation: if R-free is below the success threshold or 3+ refinement cycles have completed, STOP is removed from valid_programs until validation runs.

## Extension Points

### Adding a New Program

See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide. In summary:

1. Add program definition to `knowledge/programs.yaml` (inputs, outputs, log_parsing)
2. Add to appropriate workflow phase in `knowledge/workflows.yaml`
3. Add file categories to `knowledge/file_categories.yaml` (if new file types)
4. Add hardcoded extractor to `phenix_ai/log_parsers.py` (only if YAML patterns insufficient)

### Adding a New Workflow

1. Define states in `workflow_state.py`
2. Add transitions in workflow YAML
3. Update rules selector if needed

### Adding New Session State

1. Add to `session_state` in request schema
2. Map to `session_info` in `_process_request()`
3. Use in graph nodes as needed

## Performance Considerations

### Graph Caching

```python
_cached_graph = None

def _get_graph():
    global _cached_graph
    if _cached_graph is None:
        _cached_graph = build_agent_graph()
    return _cached_graph
```

### Request Size

- Log content can be large (truncate if needed)
- History grows with cycles (summarize old cycles)
- Files list should only include relevant files

### Response Time

- Graph invocation: ~1-5 seconds
- LLM calls: ~2-10 seconds (if used)
- Rules-only mode: <1 second

## Command Builder Architecture

The CommandBuilder provides a unified, single-entry-point system for generating
PHENIX commands. It replaces the previous fragmented approach spread across
multiple modules.

### Pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         CommandBuilder.build()                           │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────┐   ┌──────────┐ │
│  │ _select_     │──▶│ _build_      │──▶│ _apply_      │──▶│_assemble_│ │
│  │   files()    │   │  strategy()  │   │ invariants() │   │command() │ │
│  └──────────────┘   └──────────────┘   └──────────────┘   └──────────┘ │
│        │                  │                   │                 │       │
│        ▼                  ▼                   ▼                 ▼       │
│  Priority order:     Auto-fill:          Auto-fill:       Final cmd    │
│  1. LLM hints       - output_prefix     - resolution      string       │
│  2. Locked rfree    - from history      - R-free flags                 │
│  3. Best files                          - twin_law                     │
│  4. Categories                                                          │
│  5. Extensions                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### CommandContext

All parameters are passed via a `CommandContext` dataclass:

```python
@dataclass
class CommandContext:
    cycle_number: int = 1
    experiment_type: str = ""
    resolution: float = None
    best_files: Dict[str, str] = None
    rfree_mtz: str = None
    categorized_files: Dict[str, List[str]] = None
    workflow_state: str = ""
    history: List[Dict] = None
    llm_files: Dict[str, str] = None      # Optional LLM hints
    llm_strategy: Dict[str, Any] = None   # Optional LLM hints
    log: Callable = None
```

### File Selection Priority

1. **LLM hints** - If LLM specified files and they exist
2. **Locked R-free MTZ** - For MTZ slots in X-ray refinement
3. **Best files** - From BestFilesTracker (quality-scored)
4. **Categories** - From `input_priorities` in YAML
5. **Extension fallback** - Most recent file with matching extension

**Best files exclusion:** Before applying `best_files["model"]` as a model
override, the system checks if the file is in a program's `exclude_categories`.
E.g., `ligand_fit_1.pdb` (in `ligand_fit_output` → `ligand`) is excluded from
refine's model slot, falling through to category-based selection.

### Command Validation

The `validate` graph node checks that all file references in a generated command
exist in `available_files`. Output arguments (`output.file_name=X.pdb`,
`output.prefix=Y`) are stripped before extraction since they reference files
that don't exist yet. Duplicate commands and empty commands are also rejected.

### Integration with Graph Nodes

```python
# In graph_nodes.py build() function:
if USE_NEW_COMMAND_BUILDER:
    return _build_with_new_builder(state)

# _build_with_new_builder creates CommandContext from state
# and delegates to CommandBuilder.build()
```

## Automatic Error Recovery

When a PHENIX program fails with a recognized error pattern, the agent can automatically
recover and retry with corrected parameters. This is handled by `agent/error_analyzer.py`
with patterns defined in `knowledge/recoverable_errors.yaml`.

### How It Works

1. **Detection**: After a program fails, the error output is matched against
   `detection_patterns` in `recoverable_errors.yaml`
2. **Extraction**: Regex patterns extract the parameter name and available choices
   from the error message
3. **Resolution**: A strategy selects the best choice based on program type and
   workflow context (e.g., anomalous vs merged data)
4. **Retry**: The command is rebuilt with the corrected parameter and re-executed

### Currently Handled Errors

- **`ambiguous_data_labels`**: MTZ file contains multiple data arrays (e.g., both
  merged intensities and anomalous pairs). The agent picks the right array based on
  whether the program needs anomalous or merged data.
- **`ambiguous_experimental_phases`**: MTZ file contains multiple phase arrays.
  The agent selects standard vs anomalous HL coefficients based on context.

### Configuration

`knowledge/recoverable_errors.yaml` defines:
- **`errors`**: Detection patterns, extraction regexes, resolution strategies, max retries
- **`label_patterns`**: How to classify data labels as anomalous or merged
- **`program_data_preferences`**: Which programs need anomalous vs merged data
- **`context_keywords`**: Workflow hints from user advice (e.g., "SAD" → anomalous)

### Adding a New Recoverable Error

1. Add an entry under `errors:` in `recoverable_errors.yaml` with detection patterns
   and extraction regexes
2. Add a resolution handler in `error_analyzer.py` (`resolve_error()` method)
3. The fallback node in the graph will automatically use the new pattern

### Key Files

- `knowledge/recoverable_errors.yaml` — Error patterns and resolution config
- `agent/error_analyzer.py` — Detection, extraction, and resolution logic
- `agent/graph_nodes.py` — Fallback node triggers error analysis on failure

## Session Summary Generation

The agent generates structured summaries of completed sessions with optional LLM assessment.

### Summary Components

1. **Input Section**: Files, user advice, experiment type, resolution
2. **Input Data Quality**: Metrics from xtriage/mtriage (resolution, completeness, twinning)
3. **Workflow Path**: High-level description of strategy taken
4. **Steps Performed**: Table of all cycles with key metrics
5. **Final Quality**: Final metrics with quality assessments
6. **Output Files**: Key output files from the session
7. **Assessment** (optional): LLM-generated evaluation

### Quality Assessments

Metrics are automatically assessed based on resolution-dependent thresholds:

| Metric | Good | Acceptable | Needs Improvement |
|--------|------|------------|-------------------|
| R-free (2.5Å) | ≤0.25 | ≤0.30 | >0.30 |
| Map CC | ≥0.80 | ≥0.70 | <0.70 |
| Clashscore | ≤5 | ≤10 | >20 |

### LLM Assessment Prompt

When enabled, the LLM evaluates:
1. **Input Data Quality**: Resolution, completeness, issues
2. **Goal and Strategy**: User's goal and agent's approach
3. **Strategy Assessment**: Was it appropriate? Goal achieved?
4. **Current Status**: Ready for deposition?
5. **Next Steps**: Recommendations

### Usage

```python
# Generate summary without LLM
result = session.generate_agent_session_summary(include_llm_assessment=False)
print(result["markdown"])

# Get concise summary for LLM assessment
llm_input = session.get_summary_for_llm_assessment()
```

### Files

- `agent/session.py`: `generate_agent_session_summary()`, `get_summary_for_llm_assessment()`
- `analysis/agent_session_analyzer.py`: LLM assessment integration
- `phenix_ai/run_ai_analysis.py`: `run_agent_session_analysis()`
- `knowledge/prompts_hybrid.py`: Assessment prompt template

---

## Advice Preprocessing

The agent preprocesses user input to create structured, actionable guidance, with
built-in protection against prompt injection attacks and automatic change detection
for mid-session updates.

### Components

| Component | File | Function |
|-----------|------|----------|
| README discovery | `agent/advice_preprocessor.py` | `find_readme_file()` |
| README reading | `agent/advice_preprocessor.py` | `read_readme_file()` |
| Advice gathering | `agent/advice_preprocessor.py` | `gather_raw_advice()` |
| Input sanitization | `agent/advice_preprocessor.py` | `sanitize_advice()` |
| LLM preprocessing | `agent/advice_preprocessor.py` | `preprocess_advice()` |
| Change detection | `programs/ai_agent.py` | `_preprocess_user_advice()` |
| Server routing | `phenix_ai/run_ai_analysis.py` | `run_advice_preprocessing()` |

### Sources

1. **Direct advice**: `project_advice` PHIL parameter
2. **README files**: Found in `input_directory` if specified

### Process

1. **Gather**: Collect advice from all sources
2. **Hash**: Compute MD5 hash of raw advice for change detection
3. **Compare**: Check if advice changed from previous run
4. **Sanitize**: Remove potential prompt injection patterns
5. **Combine**: Merge with clear labeling
6. **Preprocess**: Use LLM to structure (optional)
7. **Store**: Save raw, processed, and hash in session

### Advice Change Detection

The agent detects when advice has changed between runs and automatically reprocesses:

- **Same hash** → reuse cached processed advice
- **Different hash** → reprocess advice AND re-extract directives
- **No new advice** → reuse cached version

This enables mid-session advice updates - users can provide new instructions
(like "stop" or "focus on the ligand") and have them take effect immediately.

### Advice Change Detection

The agent detects when advice has changed between runs and automatically reprocesses:

- **Same hash** → reuse cached processed advice
- **Different hash** → reprocess advice AND re-extract directives
- **No new advice** → reuse cached version

This enables mid-session advice updates - users can provide new instructions
(like "stop" or "focus on the ligand") and have them take effect immediately.

#### Extending a Completed Workflow (Q1)

A special case arises when the workflow has already reached the `complete` phase
and the user resumes with new advice. Two independent mechanisms both need to be
active to allow follow-up work:

**Wall 1 — AUTO-STOP in PLAN (pre-existing):**
When `metrics_trend.should_stop` is True, the PLAN node normally terminates
immediately. The `advice_changed` flag suppresses this for exactly one cycle,
letting the LLM plan before reverting to normal termination.

**Wall 2 — `valid_programs = ['STOP']` in PERCEIVE (Q1 fix):**
The `complete` phase returns `['STOP']` from `get_valid_programs`. Even with
Wall 1 down, the LLM is presented with a menu of only STOP and cannot select
any follow-up program. Fix: when `advice_changed=True` AND `phase=='complete'`,
PERCEIVE steps back to the `validate` phase, which contains:

```
['phenix.polder', 'phenix.molprobity', 'phenix.model_vs_data',
 'phenix.map_correlations', 'STOP']
```

After one successful cycle, `advice_changed` is cleared in the post-execution
check and normal termination logic resumes on the next cycle.

**Complete flow on resume with new advice:**

```
_preprocess_user_advice()
  ↓ new hash ≠ stored hash → session.data["advice_changed"] = True

PERCEIVE (graph_nodes.py)
  ↓ phase == 'complete' AND advice_changed
  ↓ valid_programs ← get_valid_programs(exp, {"phase":"validate"}, ctx)
  ↓ phase_info ← {"phase": "validate", "reason": "advice_changed: stepped back"}

PLAN (graph_nodes.py)
  ↓ metrics_trend.should_stop AND advice_changed → skip AUTO-STOP this cycle
  ↓ LLM: sees new advice + validate-phase program menu → chooses phenix.polder

BUILD / VALIDATE / OUTPUT → polder command executed

Post-execution (iterate_agent.py)
  ↓ session.data["advice_changed"] = False
  ↓ next cycle: normal AUTO-STOP and complete-phase logic resume
```

**Usage:**
```bash
# Completed workflow: xtriage → phaser → refine ×3 → ligandfit → molprobity
phenix.ai_agent \
    log_directory=AIAgent_run1 \
    restart_mode=resume \
    project_advice="also run polder on chain B residue 100"
```

**Key design note:** `phenix.polder` intentionally has no `run_once` strategy
in `programs.yaml`. Different residues and ligands may each require a separate
omit map, so `polder_done=True` does not gate polder out of `valid_programs`.

### Session Management Keywords

Two Phil parameters allow inspecting and modifying an existing session without
running new crystallographic cycles.

#### `display_and_stop`

```bash
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=basic
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=detailed
```

Prints the session history (`basic`: one line per cycle; `detailed`: full
reasoning + command), then exits. Populates `self.result` via
`_finalize_session(skip_summary=True)` so `get_results()` / `get_results_as_JSON()`
work identically to a normal run.

#### `remove_last_n`

```bash
phenix.ai_agent log_directory=AIAgent_run1 remove_last_n=2
```

Removes the last N cycles from the session JSON, clears the stale AI summary
(since the cycle set has changed), rebuilds `active_files.json` and `best_files`
from the remaining history, then saves and exits. Useful for pruning a failed
refinement run before re-running.

#### Auto-set `restart_mode=resume`

Both parameters operate on an existing session directory and therefore always
require resume semantics. `run()` automatically forces `restart_mode='resume'`
before `set_defaults()` when either parameter is set, so the user does not need
to remember the flag.

#### `get_results()` safety

`run()` sets `self.result = None` as its very first statement. `get_results()`
uses `getattr(self, 'result', None)` as a defensive fallback. This prevents
`AttributeError` on any early-exit path (session management, red-flag abort,
bad parameters) that would otherwise bypass the normal result-assignment code.

### Input Sanitization

The `sanitize_advice()` function removes potentially malicious patterns:

- Instruction override attempts ("ignore all previous instructions")
- System prompt manipulation (`<system>`, `[system]`)
- Role manipulation ("you are now a...", "act as a...")
- Hidden text (null bytes, control characters)
- Excessive repetition (potential buffer overflow)

```python
from agent.advice_preprocessor import sanitize_advice, is_suspicious

# Check if input contains suspicious patterns
if is_suspicious(user_input):
    clean_input = sanitize_advice(user_input)
```

### README Discovery

Searches for these files (case-insensitive):
- `README`, `README.txt`, `README.dat`, `README.md`
- `notes.txt`, `NOTES.txt`

Files are truncated at `max_readme_chars` (default: 5000).

### LLM Preprocessing

When `preprocess_advice=True` (default), the LLM extracts:

1. **Input Files Found**: Data files mentioned in the text
2. **Experiment Type**: SAD, MAD, MR, cryo-EM, etc.
3. **Primary Goal**: What the user wants to accomplish
4. **Key Parameters**: Wavelength, resolution, sites, heavy atom type
5. **Special Instructions**: Ligands, quality targets, etc.

### PHIL Parameters

```phil
input_directory = None
  .type = path
  .help = Directory containing input files and optional README

preprocess_advice = True
  .type = bool
  .help = Use LLM to preprocess and clarify user advice

readme_file_patterns = README README.txt README.dat README.md notes.txt
  .type = strings
  .help = Filenames to look for when extracting advice

max_readme_chars = 5000
  .type = int
  .help = Truncate README files longer than this
```

### Usage

```bash
# With direct advice
phenix.ai_agent data.mtz seq.fa project_advice="Solve by MR, R-free < 0.25"

# With README in directory
phenix.ai_agent input_directory=/data/project/

# Both combined
phenix.ai_agent input_directory=/data/project/ project_advice="Prioritize geometry"

# Update advice mid-session (on restart)
phenix.ai_agent session_file=session.json project_advice="stop"
```

### Testing

```bash
python tests/tst_advice_preprocessing.py
```

Tests cover sanitization, README discovery, advice combination, file extraction,
and change detection logic.

---

## Event System

The agent uses a structured event system for transparent decision logging.

### Event Flow

```
┌────────────────────────────────────────────────────────────────┐
│                        Graph Nodes                              │
│  perceive() → plan() → build() → validate()                    │
│       │          │         │          │                        │
│       │          │         │     ┌────┴─────┐                  │
│       │          │         │     │fallback()│                  │
│       │          │         │     └────┬─────┘                  │
│       │          │         │          │                        │
│       │          │         │     output_node()                  │
│       ▼          ▼         ▼          ▼                        │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │               state["events"] (list of dicts)            │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────┐
│              Response Building (run_ai_agent.py)               │
│  response["events"] = state["events"]                          │
│  response["events_as_simple_string"] = json.dumps(events)      │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────┐
│                    Display (ai_agent.py)                        │
│  formatter = EventFormatter(verbosity=params.verbosity)        │
│  output = formatter.format_cycle(events, cycle_number)         │
│  print(output, file=self.logger)                               │
└────────────────────────────────────────────────────────────────┘
```

### Event Types

| Type | Level | Description |
|------|-------|-------------|
| `cycle_start` | quiet | Cycle beginning |
| `cycle_complete` | quiet | Cycle finished |
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `sanity_check` | normal | Red flag or warning detected |
| `program_selected` | normal | Decision with reasoning |
| `program_modified` | normal | Program changed by rules/validation |
| `stop_decision` | normal | Whether to continue |
| `directive_applied` | normal | User directive enforced |
| `user_request_invalid` | quiet | User request unavailable |
| `files_selected` | verbose | File selection details |
| `file_scored` | verbose | Individual file scoring detail |
| `command_built` | normal | Final command |
| `thought` | verbose | LLM chain-of-thought/reasoning |
| `error` | quiet | Error occurred |
| `warning` | quiet | Non-fatal warning |
| `debug` | verbose | Internal debug information |

### Verbosity Levels

```phil
verbosity = normal
  .type = choice(quiet, normal, verbose)
```

- **quiet**: Errors, warnings, cycle summaries only
- **normal**: Key decisions and metrics (default)
- **verbose**: Full details including file selection, LLM traces, debug info

Note: `debug` is accepted as an alias for `verbose` (3 levels total).

### Implementation Files

- `agent/event_log.py` - EventType, Verbosity, EventLog class
- `agent/event_formatter.py` - EventFormatter class
- `agent/graph_nodes.py` - Event emission (`_emit()` helper)

See [TRANSPARENCY_LOGGING.md](../project/TRANSPARENCY_LOGGING.md) for full details.
