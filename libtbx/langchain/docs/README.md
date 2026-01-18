# PHENIX AI Agent - Documentation Index

This directory contains detailed technical documentation for the PHENIX AI Agent.

## Documents

| Document | Description |
|----------|-------------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | System architecture, component design, data flow |
| [API_DOCUMENTATION.md](API_DOCUMENTATION.md) | V2 API request/response schemas, transport layer |
| [TRANSPORT_REFACTORING_PLAN.md](TRANSPORT_REFACTORING_PLAN.md) | Transport layer design and implementation |
| [RED_FLAG_DETECTION_PLAN.md](RED_FLAG_DETECTION_PLAN.md) | Sanity checking, abort logic, validation gates |
| [YAML_SCORING_PLAN.md](YAML_SCORING_PLAN.md) | YAML-based program scoring system |
| [API_VERSIONING_PLAN.md](API_VERSIONING_PLAN.md) | API versioning strategy |

## Quick Start

See the main [README.md](../README.md) in the project root for:
- Overview and key features
- Directory structure
- YAML configuration guide
- Workflow diagrams (X-ray and cryo-EM)
- Operating modes (standard, rules-only, dry-run)
- How to modify agent behavior

## Key Components

### Transport Layer (`agent/transport.py`)
Handles client-server communication with:
- Request/response sanitization
- ZZxxZZ REST encoding/decoding
- Quoted string truncation
- YAML-driven configuration

### API Client (`agent/api_client.py`)
Builds V2 API requests and parses responses:
- `build_request_v2()` - Create request from client data
- `parse_response_v2()` - Extract decision from response

### Agents
- `local_agent.py` - Local execution with full encode/decode roundtrip
- `remote_agent.py` - Remote server execution via REST

## Testing

```bash
# Run all tests
python run_all_tests.py

# Run transport tests (71 tests)
python tests/test_transport.py -v

# Run specific test file
python tests/test_yaml_config.py
```
