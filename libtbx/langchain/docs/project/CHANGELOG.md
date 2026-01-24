# PHENIX AI Agent - Changelog

## Version 40 (January 2025)

### New Features
- **USER_REQUEST_INVALID event**: When user requests a program that's not available (e.g., "run phenix.xxx"), the agent now displays a prominent warning explaining why the request can't be fulfilled and what will run instead
- Warning is shown at QUIET verbosity level (always visible)
- Distinguishes between "unknown program" and "wrong workflow state"

### Files Changed
- `agent/event_log.py` - Added USER_REQUEST_INVALID event type
- `agent/event_formatter.py` - Added formatter for prominent warning display
- `agent/graph_nodes.py` - Emit event when user request detected as invalid

---

## Version 39 (January 2025)

### Bug Fixes
- **Event transport plumbing**: Fixed events not flowing through in two edge cases:
  1. Single-shot mode via `run_job_on_server` - events now decoded from server response
  2. API result retrieval via `get_results_as_JSON()` - events now serialized in output_files

### Files Changed
- `programs/ai_agent.py` - Added events serialization in `_build_output_files_from_history`
- `programs/ai_agent.py` - Added events decoding in `run_job_on_server`

---

## Version 38 (January 2025)

### Event System Phase 4: Display Integration
- Added `verbosity` parameter to `phenix.ai_agent` command
- Integrated EventFormatter for consistent output formatting
- Added `_display_cycle_events()` method for event rendering
- Legacy fallback when events not available

### Files Changed
- `programs/ai_agent.py` - Verbosity parameter, EventFormatter integration

---

## Version 37 (January 2025)

### Event System Phase 3: Transport Integration
- Events included in v2 API response schema
- LocalAgent and RemoteAgent parse events from responses
- Events stored in history_record for persistence

### Files Changed
- `phenix_ai/run_ai_agent.py` - Include events in response
- `phenix_ai/local_agent.py` - Parse events from response
- `agent/api_schema.py` - Updated response schema

---

## Version 36 (January 2025)

### Event System Phase 2: Decision Point Instrumentation
- All graph nodes now emit structured events
- Full LLM reasoning captured without truncation
- File selection reasons tracked

### Files Changed
- `agent/graph_nodes.py` - Event emission in perceive, plan, build nodes

---

## Version 34 (January 2025)

### Event System Phase 1: Core Infrastructure
- Created EventLog class for structured logging
- Created EventFormatter for human-readable output
- Defined 17 event types with verbosity levels
- LangGraph state compatibility (list of dicts)

### New Files
- `agent/event_log.py` - EventLog class, EventType constants
- `agent/event_formatter.py` - EventFormatter class

---

## Version 33 (January 2025)

### Cleanup and Production Hardening
- Removed deprecated state.md files
- Removed redundant backup files
- Fixed program registration after import changes
- Updated test suites for new structure

---

## Version 32 (January 2025)

### Pattern Centralization
- Moved all regex patterns to `knowledge/patterns.yaml`
- Created PatternManager for centralized access
- Updated log_parsers.py to use PatternManager

### New Files
- `knowledge/patterns.yaml` - Centralized regex patterns
- `agent/pattern_manager.py` - Pattern loading and compilation

---

## Version 31 (January 2025)

### Unified Command Builder
- Single CommandBuilder class for all programs
- Reads program definitions from YAML
- Consistent file selection across all programs
- Strategy flags and defaults from YAML

### Files Changed
- `agent/command_builder.py` - Complete rewrite

---

## Version 30 (January 2025)

### File Categorization Consolidation
- Centralized file categorization in `file_categorization.py`
- Semantic categories: model vs search_model distinction
- Categories defined in `file_categories.yaml`

### New Files
- `knowledge/file_categories.yaml`
- `agent/file_categorization.py` - Centralized categorization

---

## Version 29 (January 2025)

### BestFilesTracker
- New class to track best file of each type across cycles
- Scores based on metrics (R-free, resolution)
- R-free flag locking after first refinement

### New Files
- `agent/best_files_tracker.py`

---

## Version 28 (January 2025)

### YAML Configuration System
- Programs defined in `programs.yaml`
- Workflows defined in `workflows.yaml`
- Metrics defined in `metrics.yaml`
- Transport rules defined in `transport.yaml`

### New Files
- `knowledge/programs.yaml`
- `knowledge/workflows.yaml`
- `knowledge/metrics.yaml`
- `knowledge/transport.yaml`
- `knowledge/yaml_loader.py`

---

## Version 25-27 (December 2024)

### User Directives System
- Natural language directive parsing
- Stop conditions: "stop after X", "stop when metric < Y"
- Workflow preferences: "skip program", "prefer program"
- Four-layer stop condition checking

### New Files
- `agent/directive_extractor.py`
- `agent/directive_validator.py`
- `docs/guides/USER_DIRECTIVES.md`

---

## Earlier Versions

### Initial Development (2024)
- LangGraph pipeline architecture
- LLM integration (Claude, Gemini)
- Rules-only fallback mode
- Local and remote execution modes
- Session tracking and history
- Sanity checking system
