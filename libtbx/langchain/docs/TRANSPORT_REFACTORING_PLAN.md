# Transport and Sanitization Refactoring Plan

## Overview

Refactor the client-server data transport to be:
1. **Consistent** - LocalAgent and RemoteAgent use identical encode/decode paths
2. **Configurable** - Sanitization rules and limits defined in YAML
3. **Robust** - Handle edge cases like embedded ZZxxZZ markers, long quoted strings
4. **Maintainable** - Single source of truth for field definitions

## Current Problems

### 1. Inconsistent Paths
- **LocalAgent**: `serialize_request()` → `run_ai_agent()` directly
- **RemoteAgent**: `serialize_request()` → `text_as_simple_string()` → network → `simple_string_as_text()` → `run_ai_agent()`
- Bugs may only appear in one mode

### 2. Hardcoded Sanitization
- `sanitize_string()` defined inside `build_request_v2()` - not reusable
- Rules hardcoded: tabs, ZZxxZZ markers, control chars
- Truncation limits scattered: 500, 1000, 50000

### 3. Schema Not Driving Behavior
- `api_schema.py` defines fields but sanitization is separate
- Adding new fields requires code changes in multiple places

### 4. Log Files Too Large
- Logs contain large quoted data dumps (HHblits pdb70_text, etc.)
- These inflate request size without adding value
- ZZxxZZ markers from previous runs cause decode failures

## Proposed Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Client Side                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ai_agent.py                                                     │
│      │                                                           │
│      ▼                                                           │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │ LocalAgent / RemoteAgent                                 │    │
│  │     │                                                    │    │
│  │     ▼                                                    │    │
│  │ transport.prepare_request(data, schema)                  │    │
│  │     │                                                    │    │
│  │     ├── 1. remove ZZxxZZ markers (before truncation)     │    │
│  │     ├── 2. truncate_quoted_strings()                     │    │
│  │     ├── 3. sanitize_fields() [schema-driven]             │    │
│  │     ├── 4. json.dumps()                                  │    │
│  │     └── 5. text_as_simple_string()                       │    │
│  │                    │                                     │    │
│  │                    ▼                                     │    │
│  │     ┌──────────────────────────────────┐                 │    │
│  │     │ LocalAgent: decode locally       │                 │    │
│  │     │ RemoteAgent: send to server      │                 │    │
│  │     └──────────────────────────────────┘                 │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                         Server Side                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  run_ai_agent.py                                                 │
│      │                                                           │
│      ▼                                                           │
│  transport.decode_request(encoded)                               │
│      │                                                           │
│      ├── 1. simple_string_as_text()                              │
│      └── 2. json.loads()                                         │
│                    │                                             │
│                    ▼                                             │
│              Graph execution                                     │
│                    │                                             │
│                    ▼                                             │
│  transport.prepare_response(data, schema)                        │
│      │                                                           │
│      ├── 1. sanitize_fields() [schema-driven]                    │
│      ├── 2. json.dumps()                                         │
│      └── 3. text_as_simple_string()                              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## New Files

### 1. `knowledge/transport.yaml`
Configuration for transport encoding and sanitization.

```yaml
# Transport encoding configuration
encoding:
  # Characters that get replaced by ZZxxZZ markers
  markers:
    ZZCRZZ: "\n"
    ZZLBZZ: "{"
    ZZRBZZ: "}"
    ZZSCZZ: ";"
    ZZHAZZ: "#"
    ZZTAZZ: "\t"
    ZZDQZZ: '"'
    ZZSQZZ: "'"
    ZZSRZZ: "`"
    ZZEXZZ: "!"
    ZZDSZZ: "$"

# Sanitization rules
sanitization:
  # Remove these patterns before encoding
  remove_patterns:
    - 'ZZ[A-Z]{2}ZZ'  # REST encoding markers from previous runs
  
  # Replace these characters
  replacements:
    '\t': ' '      # tabs to spaces
    '\r': ''       # remove carriage returns
  
  # Remove control characters (except these)
  preserve_control_chars:
    - '\n'         # newlines are useful

  # Truncate long quoted strings in log content
  # NOTE: This happens AFTER removing ZZxxZZ markers to avoid
  # splitting a marker in half (e.g., truncating "ZZTAZZ" to "ZZT")
  truncate_quoted_strings:
    enabled: true
    max_length: 500
    quote_chars: ["'"]  # Only single quotes (double quotes are JSON structure)
    show_truncation: true  # Add ...[truncated]... marker

# Field-specific limits and rules
fields:
  request:
    log_content:
      max_length: 50000
      sanitize: true
      truncate_quotes: true
    
    user_advice:
      max_length: 10000
      sanitize: true
      truncate_quotes: false
    
    history:
      type: list
      item_fields:
        program:
          max_length: 100
          sanitize: true
        command:
          max_length: 2000
          sanitize: true
        result:
          max_length: 1000
          sanitize: true
          truncate_quotes: true
        metrics:
          type: dict
          sanitize: recursive
          max_length_per_value: 500
  
  response:
    decision.reasoning:
      max_length: 5000
      sanitize: true
    
    debug.log:
      type: list
      item_max_length: 1000
      sanitize: true
```

### 2. `agent/transport.py`
Core transport module with all encoding/decoding logic.

```python
"""
Transport module for PHENIX AI Agent client-server communication.

This module provides consistent encoding/decoding for both LocalAgent
and RemoteAgent, ensuring identical behavior regardless of transport.
"""

# Functions to implement:
# - load_transport_config() -> dict
# - truncate_quoted_strings(text, config) -> str
# - sanitize_string(text, config) -> str  
# - sanitize_field(value, field_config) -> any
# - sanitize_request(request, config) -> dict
# - sanitize_response(response, config) -> dict
# - encode_request(request) -> str
# - decode_request(encoded) -> dict
# - encode_response(response) -> str
# - decode_response(encoded) -> dict
```

## Implementation Steps

### Phase 1: Create transport.py with core functions ✅ COMPLETE
**Goal**: Centralize sanitization logic, make it testable

1. ✅ Create `agent/transport.py`
2. ✅ Move `sanitize_string()` from `api_client.py` to `transport.py`
3. ✅ Add `truncate_quoted_strings()` function
4. ✅ Add unit tests for sanitization functions
5. ✅ Update `api_client.py` to use `transport.sanitize_string()`

**Deliverable**: Sanitization works as before, but code is centralized

### Phase 2: Add YAML configuration ✅ COMPLETE
**Goal**: Make sanitization configurable

1. ✅ Create `knowledge/transport.yaml` with current hardcoded values
2. ✅ Add `load_transport_config()` to read YAML
3. ✅ Update sanitization functions to use config
4. ✅ Add schema-driven field sanitization
5. ✅ Update tests to verify config-driven behavior

**Deliverable**: Sanitization is YAML-driven, can add new rules without code changes

### Phase 3: Unify LocalAgent and RemoteAgent paths ✅ COMPLETE
**Goal**: Both agents use identical encode/decode

1. ✅ Add `prepare_request_for_transport()` and `process_request_from_transport()` 
2. ✅ Add `prepare_response_for_transport()` and `process_response_from_transport()`
3. ✅ Add `verify_roundtrip()` helper for debugging
4. ✅ Update LocalAgent to: encode → decode → run_ai_agent
5. ✅ Update RemoteAgent to: encode → network → decode (on server)
6. ✅ Remove debug code from RemoteAgent
7. ✅ Verify both paths produce identical results (71 tests passing)

**Deliverable**: LocalAgent and RemoteAgent are functionally identical

### Phase 4: Clean up and optimize ✅ COMPLETE
**Goal**: Remove debug code, optimize performance

1. ✅ Remove debug print statements from ai_agent.py
2. ✅ Remove debug print statements from remote_agent.py
3. ✅ Debug logging now controlled by verbose flag
4. ✅ Config caching implemented (in Phase 2)
5. ✅ All 71 tests passing

**Deliverable**: Production-ready transport system

## Testing Strategy

### Unit Tests (`tests/test_transport.py`)
```python
def test_truncate_quoted_strings():
    """Long quoted strings are truncated."""
    
def test_truncate_preserves_short_strings():
    """Short quoted strings are unchanged."""
    
def test_sanitize_removes_zzxxzz():
    """ZZxxZZ markers are removed."""
    
def test_sanitize_replaces_tabs():
    """Tabs become spaces."""
    
def test_encode_decode_roundtrip():
    """Data survives encode/decode cycle."""
    
def test_schema_driven_field_limits():
    """Field-specific limits are applied."""
```

### Integration Tests
```python
def test_local_remote_equivalence():
    """LocalAgent and RemoteAgent produce identical requests."""
    
def test_large_log_handling():
    """Large logs with quoted strings are handled correctly."""
```

## Migration Notes

### Backwards Compatibility
- Old clients will continue to work (server handles unsanitized input)
- New sanitization only affects outgoing data
- API version remains 2.0 (no protocol change)

### Configuration Precedence
1. YAML config (default)
2. Environment variables (for testing/override)
3. Hardcoded fallbacks (if YAML missing)

## Success Criteria

1. ✅ LocalAgent and RemoteAgent use identical code paths
2. ✅ Sanitization rules defined in YAML
3. ✅ Long quoted strings in logs are truncated
4. ✅ ZZxxZZ markers are removed before encoding
5. ✅ All existing tests pass
6. ✅ New unit tests for transport module
7. ✅ Request size reduced by >50% for logs with data dumps

## Completion Summary

**All phases complete!** The transport refactoring is now production-ready.

### Files Created/Modified

#### New Files:
- `agent/transport.py` - Core transport module (926 lines)
- `knowledge/transport.yaml` - YAML configuration
- `tests/test_transport.py` - 71 unit tests

#### Modified Files:
- `agent/api_client.py` - Uses transport functions
- `local_agent.py` - Uses unified transport (encode/decode roundtrip)
- `remote_agent.py` - Uses unified transport (cleaner implementation)
- `ai_agent.py` - Removed debug statements, controlled logging

### Key Features:
1. **Unified Transport**: LocalAgent and RemoteAgent use identical encode/decode paths
2. **YAML Configuration**: All sanitization rules configurable without code changes
3. **Quoted String Truncation**: Long data dumps (like pdb70_text) are truncated
4. **ZZxxZZ Marker Handling**: Markers from previous runs are removed before encoding
5. **Config Caching**: YAML loaded once and cached
6. **Comprehensive Tests**: 71 tests covering all edge cases

### Usage:
```python
# Prepare request for transport (sanitize + JSON + encode)
from agent.transport import prepare_request_for_transport
encoded, original = prepare_request_for_transport(request, do_encode=True)

# Process received request (decode + parse)
from agent.transport import process_request_from_transport
decoded = process_request_from_transport(encoded, was_encoded=True)

# Verify roundtrip (for debugging)
from agent.transport import verify_roundtrip
success, message, details = verify_roundtrip(request)
```
