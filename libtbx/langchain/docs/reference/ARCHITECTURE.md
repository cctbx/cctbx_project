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
Auto-generates tracking flags from `run_once: true` programs:

```
programs.yaml           program_registration.py
     │                          │
     ▼                          ▼
run_once: true ────────▶ get_trackable_programs()
                                │
                     ┌──────────┴──────────┐
                     ▼                     ▼
            workflow_state.py      workflow_engine.py
            (done flags)           (context building)
```

**Key functions:**
- `get_trackable_programs()`: Get programs with `run_once: true`
- `get_initial_history_flags()`: Get initial `{program}_done: False` flags
- `detect_programs_in_history()`: Detect completed programs in history

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
PERCEIVE ─▶ PLAN ─▶ BUILD ─▶ VALIDATE ─▶ OUTPUT
    │         │        │         │          │
    │         │        │         │          │
 Analyze   Select   Build    Validate    Format
 inputs    program  command  response   decision
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
python tests/test_advice_preprocessing.py
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
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution |
| `program_selected` | normal | Decision with reasoning |
| `user_request_invalid` | quiet | User request unavailable |
| `files_selected` | verbose | File selection details |
| `command_built` | normal | Final command |
| `error` | quiet | Error occurred |

### Verbosity Levels

```phil
verbosity = normal
  .type = choice(quiet, normal, verbose, debug)
```

- **quiet**: Errors, warnings, cycle summaries only
- **normal**: Key decisions and metrics (default)
- **verbose**: File selection details
- **debug**: Full internal state

### Implementation Files

- `agent/event_log.py` - EventType, Verbosity, EventLog class
- `agent/event_formatter.py` - EventFormatter class
- `agent/graph_nodes.py` - Event emission (`_emit()` helper)

See [TRANSPARENCY_LOGGING.md](../project/TRANSPARENCY_LOGGING.md) for full details.
