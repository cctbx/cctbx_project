# API Versioning Plan - COMPLETED

## Summary

This document describes the completed migration from the legacy v1 PHIL-based API to the unified v2 JSON API.

## Problem Statement

The original system had two different communication paths:
- **LocalAgent**: Passed parameters directly to `run_ai_agent()` with full `session_info`
- **RemoteAgent**: Used PHIL-based args via `run_job_on_server()` with limited state

This caused inconsistency - LocalAgent could use `rfree_mtz` and `best_files`, but RemoteAgent could not.

## Solution

Unified both agents on a single v2 JSON API:

```
LocalAgent                              RemoteAgent
    │                                        │
    │ build_request_v2()                     │ build_request_v2()
    ▼                                        ▼
┌─────────────────┐                  ┌─────────────────┐
│  JSON Request   │                  │  JSON Request   │
└────────┬────────┘                  └────────┬────────┘
         │                                    │
         │ (local call)                       │ (REST API)
         └──────────────┬─────────────────────┘
                        ▼
              ┌─────────────────────┐
              │   run_ai_agent.py   │
              │   run(request_json) │
              └─────────────────────┘
```

## Implementation Phases (All Complete)

### Phase 1: Schema Definition ✅

Created `knowledge/api_schema.py`:
- Request/response schema definitions
- Validation functions
- Default value handling
- Version helpers

Created `agent/api_client.py`:
- `build_request_v2()` - Build requests
- `parse_response_v2()` - Parse responses
- Serialization functions

### Phase 2: Integration ✅

Updated `RemoteAgent` to use v2 JSON:
- Builds request with `build_request_v2()`
- Sends JSON to server
- Parses JSON response

Updated `run_ai_agent.py` to accept v2 JSON:
- Single `run(request_json)` entry point
- Validates and processes requests
- Returns JSON response

### Phase 3: Session State ✅

Full session state now transmitted:
- `session_state.resolution`
- `session_state.experiment_type`
- `session_state.rfree_mtz`
- `session_state.best_files`

Server uses all state for better decisions.

### Phase 4: Unified Architecture ✅

Removed all v1 code:
- No more PHIL-based parameter passing
- No more separate code paths
- Both agents use identical JSON format

Updated `LocalAgent` to use v2:
- Builds JSON request like RemoteAgent
- Calls `run_ai_agent.run(request_json=...)`
- Same response format

## Final Architecture

### Files Changed

| File | Before | After | Change |
|------|--------|-------|--------|
| `local_agent.py` | 68 lines (params) | 142 lines (JSON) | Unified |
| `remote_agent.py` | 403 lines (v1+v2) | 211 lines (v2 only) | -192 lines |
| `run_ai_agent.py` | 656 lines (v1+v2) | 328 lines (v2 only) | -328 lines |
| `agent/api_client.py` | 518 lines | 450 lines | -68 lines |

### Total Reduction

- Before: ~1,645 lines across 4 files
- After: ~1,131 lines across 4 files
- **Removed: ~514 lines of duplicate/legacy code**

## Request Schema

```json
{
  "api_version": "2.0",
  "files": [...],
  "cycle_number": 5,
  "log_content": "...",
  "history": [...],
  "session_state": {
    "resolution": 2.5,
    "experiment_type": "xray",
    "rfree_mtz": "/path/to/data.mtz",
    "best_files": {"model": "/path/to/model.pdb"}
  },
  "user_advice": "...",
  "settings": {
    "provider": "google",
    "abort_on_red_flags": true
  }
}
```

## Response Schema

```json
{
  "api_version": "2.0",
  "decision": {
    "program": "phenix.refine",
    "command": "phenix.refine model.pdb data.mtz",
    "reasoning": "..."
  },
  "stop": false,
  "metadata": {
    "experiment_type": "xray",
    "workflow_state": "xray_refined"
  },
  "debug": {
    "timing_ms": 1234
  }
}
```

## Benefits Achieved

1. **Consistency**: Both agents use identical communication
2. **Full State**: Server receives all session state (rfree_mtz, best_files)
3. **Simplicity**: Single code path, no version switching
4. **Maintainability**: 514 fewer lines of code
5. **Testability**: One format to test

## Tests

All tests pass:
- `test_api_schema.py`: 25 tests
- `test_best_files_tracker.py`: 37 tests

Run with:
```bash
python tests/run_all_tests.py
```
