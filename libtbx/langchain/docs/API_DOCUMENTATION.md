# PHENIX AI Agent API Documentation

## Overview

The PHENIX AI Agent uses a client-server architecture where:
- **Client**: User's PHENIX installation (runs `ai_agent.py`)
- **Server**: Decision-making engine (runs `run_ai_agent.py`)

Both LocalAgent (same process) and RemoteAgent (REST server) use identical v2 JSON API **and transport encoding** for consistency.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Client Side                              │
├─────────────────────────────────────────────────────────────────┤
│  LocalAgent / RemoteAgent                                        │
│       │                                                          │
│       ▼                                                          │
│  build_request_v2() + prepare_request_for_transport()            │
│       │                                                          │
│       ├── sanitize (remove ZZxxZZ, truncate quotes, etc.)        │
│       ├── json.dumps()                                           │
│       └── encode_for_rest() [ZZxxZZ markers]                     │
│                    │                                             │
│  LocalAgent: decode locally    RemoteAgent: send to server       │
└─────────────────────────────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                      run_ai_agent.py                             │
│       │                                                          │
│       ▼                                                          │
│  process_request_from_transport()                                │
│       │                                                          │
│       ├── decode_from_rest()                                     │
│       └── json.loads()                                           │
│                    │                                             │
│                    ▼                                             │
│              Graph execution → JSON Response                     │
└─────────────────────────────────────────────────────────────────┘
```

## Transport Layer

The transport layer (`agent/transport.py`) ensures reliable communication:

### Sanitization (before encoding)
1. Remove ZZxxZZ markers from log content (prevents double-encoding)
2. Truncate long quoted strings (e.g., pdb70_text='...')
3. Replace tabs with spaces
4. Remove control characters (except newline)

### Configuration
Transport settings are defined in `knowledge/transport.yaml`:
- Field-specific max lengths
- Quote truncation settings
- Sanitization patterns

### Usage
```python
from agent.transport import (
    prepare_request_for_transport,
    process_request_from_transport,
    verify_roundtrip,
)

# Encode request
encoded, original = prepare_request_for_transport(request, do_encode=True)

# Decode request
decoded = process_request_from_transport(encoded, was_encoded=True)

# Verify integrity
success, msg, details = verify_roundtrip(request)
```

## Request Schema

```json
{
  "api_version": "2.0",
  "client_version": "1.22.0",
  
  "log_content": "string - log text from last command",
  
  "files": [
    "/absolute/path/to/file1.mtz",
    "/absolute/path/to/file2.pdb"
  ],
  
  "history": [
    {
      "cycle": 1,
      "program": "phenix.xtriage",
      "command": "phenix.xtriage data.mtz",
      "result": "SUCCESS",
      "output_files": ["/path/to/output.log"]
    }
  ],
  
  "session_state": {
    "resolution": 2.5,
    "experiment_type": "xray",
    "rfree_mtz": "/path/to/refine_001_data.mtz",
    "best_files": {
      "model": "/path/to/best.pdb",
      "mtz": "/path/to/best.mtz"
    }
  },
  
  "user_advice": "Solve the structure using data to 3 A",
  
  "settings": {
    "provider": "google",
    "abort_on_red_flags": true,
    "abort_on_warnings": false,
    "max_cycles": 20,
    "use_rules_only": false
  },
  
  "cycle_number": 5
}
```

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `api_version` | string | Must be "2.0" |
| `files` | array | List of available file paths |
| `cycle_number` | integer | Current cycle number |

### Optional Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `client_version` | string | null | Client software version |
| `log_content` | string | "" | Log text from previous command |
| `history` | array | [] | Previous cycle records |
| `session_state` | object | {} | Session state (see below) |
| `user_advice` | string | "" | User instructions |
| `settings` | object | {} | Execution settings |

### Session State Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `resolution` | float | null | Data resolution in Angstroms |
| `experiment_type` | string | null | "xray" or "cryoem" |
| `rfree_mtz` | string | null | Path to locked R-free MTZ |
| `best_files` | object | {} | Best files by category |

### Settings Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `provider` | string | "google" | LLM provider |
| `abort_on_red_flags` | boolean | true | Abort on critical issues |
| `abort_on_warnings` | boolean | false | Abort on warnings |
| `max_cycles` | integer | 20 | Maximum cycles |
| `use_rules_only` | boolean | false | Skip LLM, use rules only |

## Response Schema

```json
{
  "api_version": "2.0",
  "server_version": "1.23.0",
  
  "decision": {
    "program": "phenix.refine",
    "command": "phenix.refine model.pdb data.mtz output.prefix=refine_001",
    "reasoning": "The model needs refinement to improve R-free...",
    "strategy": {
      "output_prefix": "refine_001",
      "resolution": 3.0
    },
    "confidence": "high"
  },
  
  "stop": false,
  "stop_reason": null,
  
  "metadata": {
    "experiment_type": "xray",
    "workflow_state": "xray_refined",
    "warnings": [],
    "red_flags": []
  },
  
  "debug": {
    "log": ["PERCEIVE: ...", "PLAN: ...", "BUILD: ..."],
    "timing_ms": 1234
  },
  
  "error": null
}
```

### Decision Fields

| Field | Type | Description |
|-------|------|-------------|
| `program` | string | Selected program name |
| `command` | string | Full command to execute |
| `reasoning` | string | Explanation of decision |
| `strategy` | object | Strategy options used |
| `confidence` | string | "high", "medium", "low", "unknown" |

### Stop Response

When workflow should stop:

```json
{
  "decision": {
    "program": "STOP",
    "command": "STOP",
    "reasoning": "R-free converged at 0.22"
  },
  "stop": true,
  "stop_reason": "converged"
}
```

Stop reasons:
- `converged` - Quality metrics indicate completion
- `red_flag` - Critical issue detected
- `max_cycles` - Reached cycle limit
- `error` - Processing error

## Agent Classes

### LocalAgent

Runs the decision-making graph locally (same process).

```python
from local_agent import LocalAgent

agent = LocalAgent(params, logger=sys.stdout)
result = agent.decide_next_step(
    log_content="...",
    history=[...],
    files=[...],
    guidelines="...",
    session_resolution=2.5,
    session_info={
        "experiment_type": "xray",
        "rfree_mtz": "/path/to/data.mtz",
        "best_files": {"model": "/path/to/model.pdb"}
    }
)
```

### RemoteAgent

Sends requests to PHENIX REST server for decision-making.

```python
from remote_agent import RemoteAgent

agent = RemoteAgent(params, logger=sys.stdout)
result = agent.decide_next_step(
    log_content="...",
    history=[...],
    files=[...],
    guidelines="...",
    session_resolution=2.5,
    session_info={
        "experiment_type": "xray",
        "rfree_mtz": "/path/to/data.mtz",
        "best_files": {"model": "/path/to/model.pdb"}
    }
)
```

### Common Interface

Both agents implement the same interface:

```python
def decide_next_step(
    self,
    log_content,              # Log text from previous command
    history,                  # List of previous cycle records
    files,                    # List of available file paths
    guidelines="",            # User instructions
    session_resolution=None,  # Resolution in Angstroms
    session_info=None,        # Session state dict
    abort_on_red_flags=True,  # Abort on critical issues
    abort_on_warnings=False   # Abort on warnings
) -> dict:
    """Returns history_record with decision."""
```

Both agents:
1. Build identical v2 JSON requests using `build_request_v2()`
2. Pass them to `run_ai_agent.run(request_json=...)`
3. Receive identical `history_record` responses

## Session Info Structure

The `session_info` dict carries tracked state from the client:

```python
session_info = {
    # Locked experiment type (prevents switching mid-workflow)
    "experiment_type": "xray",  # or "cryoem"
    
    # Locked R-free MTZ (ensures consistent R-free flags)
    "rfree_mtz": "/path/to/refine_001_data.mtz",
    
    # Best files by category (from BestFilesTracker)
    "best_files": {
        "model": "/path/to/best_model.pdb",
        "mtz": "/path/to/best_data.mtz",
        "map": "/path/to/best_map.ccp4"
    }
}
```

## Error Handling

### Client-Side

```python
result = agent.decide_next_step(...)
if result is None:
    # Communication failed
    handle_error()
elif result.get("error"):
    # Server returned error
    print(f"Error: {result['error']}")
elif result.get("stop"):
    # Workflow should stop
    reason = result.get("stop_reason")
    if reason == "red_flag":
        print(f"Aborted: {result.get('abort_message')}")
```

### Server-Side

```python
# Validation errors
if not is_valid:
    return _build_error_response(f"Invalid request: {errors}")

# Processing errors
try:
    final_state = graph.invoke(initial_state)
except Exception as e:
    return _build_error_response(f"Graph execution failed: {e}")
```

## Extending the API

### Adding New Request Fields

1. Add field to `REQUEST_V2_SCHEMA` in `api_schema.py` with default value
2. Update `_process_request()` in `run_ai_agent.py` to use the field
3. Update `build_request_v2()` in `api_client.py` to accept the field
4. Existing code automatically uses default value

### Adding New Response Fields

1. Add field to `RESPONSE_V2_SCHEMA` in `api_schema.py`
2. Update `create_response()` to include the field
3. Update `parse_response_v2()` to extract the field
4. Existing code ignores unknown fields

## Files Reference

| File | Lines | Description |
|------|-------|-------------|
| `local_agent.py` | 142 | Builds v2 request, calls run_ai_agent locally |
| `remote_agent.py` | 211 | Builds v2 request, sends to REST server |
| `run_ai_agent.py` | 328 | Entry point, parses JSON, runs graph |
| `knowledge/api_schema.py` | ~450 | Schema definitions, validation, defaults |
| `agent/api_client.py` | ~450 | Request building, response parsing |
| `tests/test_api_schema.py` | ~760 | API tests |
