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
│  │  │  └─────────┘  └─────────┘  └─────────┘  └─────────┘  └───────┘│││
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

#### Log Parsers (log_parsers.py)
Extracts metrics from program output:
- R-factors, resolution, clashscore
- Output file paths
- Success/failure status

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
PERCEIVE ─▶ PLAN ─▶ BUILD ─▶ VALIDATE ─▶ OUTPUT
    │         │        │         │          │
    │         │        │         │          │
 Analyze   Select   Build    Validate    Format
 inputs    program  command  response   decision
```

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
    ...
)
```

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

```
improved_agent_v2/
├── ai_agent.py              # Main client orchestrator
├── local_agent.py           # Local agent (142 lines)
├── remote_agent.py          # Remote agent (211 lines)
├── run_ai_agent.py          # Decision engine entry point (328 lines)
├── session.py               # Session state tracking
├── utilities.py             # Shared utilities
├── log_parsers.py           # Log parsing functions
│
├── agent/
│   ├── api_client.py        # v2 request/response adapters
│   ├── graph.py             # LangGraph definition
│   ├── graph_state.py       # Graph state schema
│   ├── graph_nodes.py       # Graph node implementations
│   ├── rules_selector.py    # Program selection rules
│   ├── template_builder.py  # Command template builder
│   └── workflow_state.py    # Workflow state detection
│
├── knowledge/
│   ├── api_schema.py        # v2 API schema definitions
│   ├── best_files_config.yaml
│   ├── programs/            # YAML program definitions
│   └── workflows/           # Workflow configurations
│
├── tests/
│   ├── run_all_tests.py     # Unified test runner
│   ├── test_api_schema.py   # API tests (25 tests)
│   └── test_best_files_tracker.py  # Tracker tests (37 tests)
│
└── docs/
    ├── README.md
    ├── API_DOCUMENTATION.md
    └── ARCHITECTURE.md
```

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

## Extension Points

### Adding a New Program

1. Create YAML definition in `knowledge/programs/`
2. Add to workflow transitions if needed
3. Add log parser if output format is unique

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

### Integration with Graph Nodes

```python
# In graph_nodes.py build() function:
if USE_NEW_COMMAND_BUILDER:
    return _build_with_new_builder(state)

# _build_with_new_builder creates CommandContext from state
# and delegates to CommandBuilder.build()
```
