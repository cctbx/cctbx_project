# Implementation Plan: Reducing Program Configuration Points

## Status: ✅ COMPLETED (January 2026)

This document describes the implementation of a centralized YAML-driven architecture that reduced the number of files needed when adding a new PHENIX program from 7 files to 2-3 files.

---

## Summary of Changes

### Before (7 files required)
1. `knowledge/programs.yaml` - Program definition
2. `knowledge/workflows.yaml` - Workflow phase
3. `phenix_ai/log_parsers.py` - Metric extraction function
4. `agent/workflow_state.py` - `<program>_done` flag
5. `agent/workflow_engine.py` - Context building
6. `agent/session.py` - Summary display
7. `agent/directive_extractor.py` - Tutorial patterns (optional)

### After (2-3 files required)
1. `knowledge/programs.yaml` - Complete definition including metrics ✅
2. `knowledge/workflows.yaml` - Workflow phase ✅
3. `agent/directive_extractor.py` - Tutorial patterns (optional)

---

## Implementation Phases

### Phase 1: Centralized Metric Patterns ✅ COMPLETED

**Goal**: Define metric extraction patterns once in `programs.yaml`, use them everywhere.

**Files Created:**
- `knowledge/metric_patterns.py` - Centralized pattern loader
- `tests/test_metric_patterns.py` - 11 test cases

**Files Modified:**
- `knowledge/programs.yaml` - Extended `log_parsing` sections with `display_name`, `summary_format`, `no_match_pattern`, `no_match_value`
- `phenix_ai/log_parsers.py` - Uses YAML patterns via `extract_metrics_for_program()`
- `agent/session.py` - Uses YAML patterns for metric extraction

**How it works:**
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

---

### Phase 2: Auto-Registration of Programs ✅ COMPLETED

**Goal**: Automatically generate `<program>_done` tracking flags from `programs.yaml`.

**Files Created:**
- `knowledge/program_registration.py` - Auto-registration module
- `tests/test_program_registration.py` - 13 test cases

**Files Modified:**
- `agent/workflow_state.py` - Uses `get_initial_history_flags()` and `detect_programs_in_history()`
- `agent/workflow_engine.py` - Uses `get_all_done_flags()` for context building

**How it works:**
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

---

### Phase 3: Declarative Summary Display ✅ COMPLETED

**Goal**: Define how metrics appear in session summaries via YAML configuration.

**Files Created:**
- `knowledge/summary_display.py` - YAML-driven formatting
- `tests/test_summary_display.py` - 12 test cases

**Files Modified:**
- `knowledge/metrics.yaml` - Added `summary_display` section with `quality_table` and `step_metrics`
- `agent/session.py` - Uses `format_quality_table_rows()` and `format_step_metric()`

**How it works:**
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

---

### Phase 4: Integration and Cleanup ✅ COMPLETED

**Changes:**
- Removed backward compatibility fallbacks from all modules
- Updated `docs/ADDING_PROGRAMS.md` to reflect simplified 2-3 file process
- Updated `docs/ARCHITECTURE.md` with new module descriptions
- Verified all tests pass

---

## New Modules

### knowledge/metric_patterns.py

| Function | Purpose |
|----------|---------|
| `get_all_metric_patterns()` | Load and compile all patterns from YAML |
| `extract_metrics_for_program()` | Extract metrics using YAML patterns |
| `get_metric_display_config()` | Get display configuration for a metric |
| `format_metric_value()` | Format using YAML format strings |

### knowledge/program_registration.py

| Function | Purpose |
|----------|---------|
| `get_trackable_programs()` | Get all programs with `run_once: true` |
| `get_initial_history_flags()` | Get initial `{program}_done: False` dict |
| `detect_programs_in_history()` | Detect completed programs in history |
| `get_all_done_flags()` | Get list of all flag names |

### knowledge/summary_display.py

| Function | Purpose |
|----------|---------|
| `get_quality_table_config()` | Load quality table row configs |
| `get_step_metrics_config()` | Load per-program step metric configs |
| `format_quality_table_rows()` | Format Final Quality table |
| `format_step_metric()` | Format a step's key metric |

---

## Test Coverage

| Test File | Test Count | Status |
|-----------|------------|--------|
| `test_metric_patterns.py` | 11 | ✅ PASS |
| `test_program_registration.py` | 13 | ✅ PASS |
| `test_summary_display.py` | 12 | ✅ PASS |

**Total new tests: 36**

---

## Benefits

1. **Single Source of Truth**: Metric patterns defined once in `programs.yaml`
2. **Automatic Flag Generation**: `run_once: true` creates tracking flags automatically
3. **Declarative Display**: Summary display configured in YAML, not code
4. **Reduced Maintenance**: Adding programs requires 2-3 files instead of 7
5. **Consistency**: Same patterns used everywhere (extraction, display)
6. **Testability**: Each module has comprehensive tests

---

## How to Add a New Program

See [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) for the complete guide.

**Quick Summary:**

1. **`programs.yaml`**: Define program with `log_parsing` patterns and optionally `run_once: true`
2. **`workflows.yaml`**: Add to appropriate phase
3. **`metrics.yaml`** (optional): Add to `summary_display` if custom display needed
4. **`directive_extractor.py`** (optional): Add tutorial patterns

That's it - metric extraction, history tracking, and summary display are handled automatically!
