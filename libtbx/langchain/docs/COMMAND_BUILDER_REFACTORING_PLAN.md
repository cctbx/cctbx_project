# Command Builder Refactoring Plan

## Status

| Phase | Status | Notes |
|-------|--------|-------|
| Phase 1: Create New Module | ✅ COMPLETE | `command_builder.py` created with CommandBuilder class |
| Phase 2: Add Compatibility Layer | ✅ COMPLETE | TemplateBuilder.USE_NEW_BUILDER flag added |
| Phase 3: Update Graph Nodes | ✅ COMPLETE | graph_nodes.USE_NEW_COMMAND_BUILDER flag added |
| Phase 4: Remove Old Code | ⏸️ DEFERRED | Will remove after production testing |
| Phase 5: Testing & Documentation | 🔄 IN PROGRESS | Basic tests done, more needed |

### How to Enable New Builder

```python
# Option 1: Enable globally in graph_nodes.py
from agent import graph_nodes
graph_nodes.USE_NEW_COMMAND_BUILDER = True

# Option 2: Enable in template_builder.py (for fallback path)
from agent.template_builder import TemplateBuilder
TemplateBuilder.USE_NEW_BUILDER = True
```

---

## Executive Summary

Consolidate the current fragmented command generation system (4 entry points across 2 overlapping modules) into a single, clean `CommandBuilder` class with one entry point and a clear pipeline.

---

## Current State Analysis

### Current Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        programs.yaml                             │
│  command templates, inputs, strategy_flags, defaults, invariants │
└─────────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┴─────────────────────┐
        ▼                                           ▼
┌───────────────────────┐                 ┌───────────────────────┐
│  program_registry.py  │                 │  template_builder.py  │
│  - build_command()    │◄────────────────│  - build_command()    │
│  - get_program()      │    delegates    │  - build_command_     │
│  - get_strategy_flags │                 │    for_program()      │
│  - get_defaults()     │                 │  - validate_and_fix() │
└───────────────────────┘                 └───────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        graph_nodes.py                            │
│  - build_node(): 2 paths (LLM with/without files)               │
│  - fallback_build_node(): rules-only path                        │
│  - scattered file selection/override logic                       │
└─────────────────────────────────────────────────────────────────┘
```

### Current Entry Points

| Entry Point | Location | When Used |
|-------------|----------|-----------|
| `builder.build_command(prog, files, strategy)` | graph_nodes.py:1275 | LLM provided files |
| `builder.build_command_for_program(prog, files, ...)` | graph_nodes.py:1151 | LLM chose program only |
| `builder.build_command_for_program(prog, files, context, ...)` | graph_nodes.py:1487 | Rules-only fallback |
| `registry.build_command(prog, files, strategy)` | program_registry.py:722 | Direct utility call |

### Problems

1. **Fragmentation**: 4 entry points, unclear which to use when
2. **Overlapping modules**: `template_builder.py` often just wraps `program_registry.py`
3. **Inconsistent invariant application**: Sometimes applied, sometimes not
4. **Scattered file selection**: Logic spread across graph_nodes.py (lines 1046-1240)
5. **Ad-hoc context passing**: `best_files`, `rfree_mtz`, `context` as separate params
6. **Hard to test**: No single point to unit test command generation

---

## Proposed Architecture

### New Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        programs.yaml                             │
│                      (format unchanged)                          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                     command_builder.py                           │
│                                                                  │
│  CommandContext (dataclass)                                      │
│  ├── cycle_number: int                                          │
│  ├── experiment_type: str                                       │
│  ├── resolution: float                                          │
│  ├── best_files: dict                                           │
│  ├── rfree_mtz: str                                             │
│  ├── categorized_files: dict                                    │
│  ├── workflow_state: str                                        │
│  └── llm_hints: dict (files, strategy from LLM)                 │
│                                                                  │
│  CommandBuilder                                                  │
│  └── build(program, available_files, context) ◄── SINGLE ENTRY  │
│      │                                                          │
│      ├── 1. _select_files()                                     │
│      ├── 2. _build_strategy()                                   │
│      ├── 3. _apply_invariants()                                 │
│      └── 4. _assemble_command()                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        graph_nodes.py                            │
│  - build_node(): creates CommandContext, calls builder.build()  │
│  - fallback_build_node(): same pattern                          │
│  - NO file selection logic here                                 │
└─────────────────────────────────────────────────────────────────┘
```

### The Pipeline

```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│   SELECT    │    │    BUILD    │    │    APPLY    │    │  ASSEMBLE   │
│   FILES     │───►│  STRATEGY   │───►│ INVARIANTS  │───►│  COMMAND    │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
      │                  │                  │                  │
      ▼                  ▼                  ▼                  ▼
 Use priorities    From LLM hints     Auto-fill res      Template +
 best_files        or defaults        Auto-fill prefix   files +
 rfree_mtz lock                       Validate files     strategy +
 categories                                              defaults
```

---

## Implementation Plan

### Phase 1: Create New Module (Non-Breaking)

**Goal**: Create `command_builder.py` alongside existing modules.

#### Step 1.1: Create CommandContext dataclass

```python
# agent/command_builder.py

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any

@dataclass
class CommandContext:
    """All context needed for command generation."""
    
    # Cycle info
    cycle_number: int = 1
    
    # Experiment info
    experiment_type: str = ""  # 'xray' or 'cryoem'
    resolution: Optional[float] = None
    
    # File tracking
    best_files: Dict[str, str] = field(default_factory=dict)
    rfree_mtz: Optional[str] = None  # Locked R-free MTZ (X-ray only)
    categorized_files: Dict[str, List[str]] = field(default_factory=dict)
    
    # Workflow info
    workflow_state: str = ""
    
    # LLM suggestions (optional - used when LLM provides hints)
    llm_files: Optional[Dict[str, Any]] = None  # {slot: path} from LLM
    llm_strategy: Optional[Dict[str, Any]] = None  # {flag: value} from LLM
    
    @classmethod
    def from_state(cls, state: dict) -> 'CommandContext':
        """Create context from graph state dict."""
        session_info = state.get("session_info", {})
        workflow_state = state.get("workflow_state", {})
        
        return cls(
            cycle_number=state.get("cycle_number", 1),
            experiment_type=session_info.get("experiment_type", ""),
            resolution=state.get("resolution") or workflow_state.get("resolution"),
            best_files=session_info.get("best_files", {}),
            rfree_mtz=session_info.get("rfree_mtz"),
            categorized_files=workflow_state.get("categorized_files", {}),
            workflow_state=workflow_state.get("current_state", ""),
            llm_files=state.get("llm_files"),
            llm_strategy=state.get("llm_strategy"),
        )
```

#### Step 1.2: Create CommandBuilder class skeleton

```python
class CommandBuilder:
    """
    Single entry point for all command generation.
    
    Usage:
        builder = CommandBuilder()
        context = CommandContext.from_state(state)
        command = builder.build("phenix.refine", available_files, context)
    """
    
    def __init__(self):
        # Reuse existing registry for YAML access
        from agent.program_registry import ProgramRegistry
        self._registry = ProgramRegistry()
    
    def build(self, program: str, available_files: List[str], 
              context: CommandContext) -> Optional[str]:
        """
        Build a command for the given program.
        
        This is the SINGLE ENTRY POINT for all command generation.
        
        Args:
            program: Program name (e.g., "phenix.refine")
            available_files: List of available file paths
            context: CommandContext with all needed information
            
        Returns:
            Command string, or None if command cannot be built
        """
        # 1. Select files
        files = self._select_files(program, available_files, context)
        if files is None:
            return None
        
        # 2. Build strategy
        strategy = self._build_strategy(program, context)
        
        # 3. Apply invariants (may modify files and strategy)
        files, strategy = self._apply_invariants(program, files, strategy, context)
        
        # 4. Assemble final command
        command = self._assemble_command(program, files, strategy)
        
        return command
```

#### Step 1.3: Implement _select_files()

Move file selection logic from:
- `template_builder.py` lines 175-340 (`build_command_for_program`)
- `graph_nodes.py` lines 1046-1140 (model override logic)
- `graph_nodes.py` lines 1168-1240 (auto-fill logic)

Into single method with clear priority order:
1. LLM hints (if provided and valid)
2. Locked rfree_mtz (for MTZ slots in X-ray)
3. Best files (from BestFilesTracker)
4. Category-based selection (input_priorities)
5. Extension-based fallback

#### Step 1.4: Implement _build_strategy()

```python
def _build_strategy(self, program: str, context: CommandContext) -> Dict[str, Any]:
    """Build initial strategy from context and LLM hints."""
    strategy = {}
    
    # Start with LLM suggestions if provided
    if context.llm_strategy:
        strategy.update(context.llm_strategy)
    
    return strategy
```

#### Step 1.5: Implement _apply_invariants()

Move from `template_builder.py` `validate_and_fix()` method.

```python
def _apply_invariants(self, program: str, files: dict, strategy: dict,
                      context: CommandContext) -> Tuple[dict, dict]:
    """Apply invariants - may modify files and strategy."""
    
    invariants = self._registry.get_invariants(program)
    
    for inv in invariants:
        files, strategy, msg = self._apply_single_invariant(
            inv, program, files, strategy, context
        )
        if msg:
            self._log("INVARIANT: %s" % msg)
    
    return files, strategy
```

#### Step 1.6: Implement _assemble_command()

Use existing `program_registry.build_command()` logic but cleaner:

```python
def _assemble_command(self, program: str, files: dict, 
                      strategy: dict) -> Optional[str]:
    """Assemble final command from program, files, and strategy."""
    
    # Get command template
    prog_def = self._registry.get_program(program)
    template = prog_def.get("command", program)
    
    # Substitute file placeholders
    cmd = self._substitute_files(template, files, prog_def)
    
    # Add strategy flags
    cmd = self._add_strategy_flags(cmd, strategy, prog_def)
    
    # Add defaults
    cmd = self._add_defaults(cmd, strategy, prog_def)
    
    return cmd
```

### Phase 2: Add Compatibility Layer

**Goal**: Make new builder work with existing code without breaking anything.

#### Step 2.1: Add wrapper in template_builder.py

```python
# template_builder.py

def build_command_for_program(self, program, available_files, ...):
    """DEPRECATED: Use CommandBuilder.build() instead."""
    
    # Create context from parameters
    context = CommandContext(
        cycle_number=context_dict.get("cycle_number", 1) if context_dict else 1,
        resolution=context_dict.get("session_resolution") if context_dict else None,
        best_files=best_files or {},
        rfree_mtz=rfree_mtz,
        categorized_files=categorized_files or {},
    )
    
    # Delegate to new builder
    from agent.command_builder import CommandBuilder
    builder = CommandBuilder()
    return builder.build(program, available_files, context)
```

#### Step 2.2: Add tests comparing old vs new

```python
# tests/test_command_builder.py

def test_compatibility_refine():
    """Ensure new builder produces same output as old."""
    files = ["model.pdb", "data.mtz"]
    
    # Old way
    old_cmd = template_builder.build_command_for_program("phenix.refine", files, ...)
    
    # New way
    context = CommandContext(cycle_number=1, ...)
    new_cmd = command_builder.build("phenix.refine", files, context)
    
    assert old_cmd == new_cmd
```

### Phase 3: Update Graph Nodes

**Goal**: Simplify graph_nodes.py to use new builder.

#### Step 3.1: Simplify build_node()

```python
# graph_nodes.py

def build_node(state):
    """Build command using LLM decision."""
    
    program = state.get("program")
    available_files = state.get("available_files", [])
    
    # Create context from state (includes LLM hints)
    context = CommandContext.from_state(state)
    context.llm_files = state.get("corrected_files")  # From LLM
    context.llm_strategy = state.get("strategy")  # From LLM
    
    # Single call to builder
    builder = get_command_builder()
    command = builder.build(program, available_files, context)
    
    if not command:
        return {**state, "command": "", "validation_error": "Failed to build command"}
    
    return {**state, "command": command, "validation_error": None}
```

#### Step 3.2: Simplify fallback_build_node()

```python
def fallback_build_node(state):
    """Build command using rules only (no LLM)."""
    
    available_files = state.get("available_files", [])
    context = CommandContext.from_state(state)
    # No llm_files or llm_strategy - pure rule-based
    
    builder = get_command_builder()
    
    # Try each valid program
    for program in get_valid_programs(context.workflow_state):
        command = builder.build(program, available_files, context)
        if command and not is_duplicate(command, state):
            return {**state, "command": command, "fallback_used": True}
    
    return {**state, "command": "", "validation_error": "No valid command found"}
```

### Phase 4: Remove Old Code (DEFERRED)

**Status**: DEFERRED until new CommandBuilder is tested in production.

**Rationale**: The feature flag approach (`USE_NEW_COMMAND_BUILDER`) allows safe testing
of the new builder alongside existing code. Old code will be removed once:
1. New builder has been tested with real workflows
2. No regressions are found
3. Feature flag has been enabled in production for a reasonable period

**Goal**: Delete deprecated modules and clean up.

#### Step 4.1: Remove template_builder.py methods (DEFERRED)

Methods marked as deprecated but NOT yet deleted:
- `build_command()` - DEPRECATED, use CommandBuilder.build()
- `build_command_for_program()` - DEPRECATED, use CommandBuilder.build()
- `validate_and_fix()` - Keep for backward compatibility (used by tests)

#### Step 4.2: Simplify program_registry.py (DEFERRED)

Keep all methods for now. Will remove `build_command()` after migration complete.

#### Step 4.3: Update imports (DEFERRED)

Will update imports after old code is removed.

### Phase 5: Testing & Documentation

#### Step 5.1: Comprehensive tests

```
tests/
├── test_command_builder.py
│   ├── test_select_files_*  (file selection tests)
│   ├── test_build_strategy_*  (strategy building tests)
│   ├── test_apply_invariants_*  (invariant tests)
│   ├── test_assemble_command_*  (assembly tests)
│   └── test_build_*  (end-to-end tests)
```

#### Step 5.2: Update documentation

- Update ARCHITECTURE.md with new command building flow
- Add docstrings to all public methods
- Add examples in command_builder.py header

---

## File Changes Summary

### New Files
| File | Description |
|------|-------------|
| `agent/command_builder.py` | New unified command builder |
| `tests/test_command_builder.py` | Tests for new builder |

### Modified Files
| File | Changes |
|------|---------|
| `agent/graph_nodes.py` | Simplify build_node, fallback_build_node |
| `agent/template_builder.py` | Phase 2: compatibility wrapper, Phase 4: delete most code |
| `agent/program_registry.py` | Phase 4: remove build_command |
| `docs/ARCHITECTURE.md` | Update with new flow |

### Deleted Files (Phase 4)
| File | Reason |
|------|--------|
| None | We simplify rather than delete files |

---

## Rollback Plan

If issues arise:

1. **Phase 1-2**: New code is additive, just delete `command_builder.py`
2. **Phase 3**: Revert graph_nodes.py changes
3. **Phase 4**: Don't do until Phase 3 is stable in production

---

## Success Criteria

1. **Single entry point**: All command generation goes through `CommandBuilder.build()`
2. **Testable**: Can unit test each pipeline stage independently
3. **Understandable**: New developer can understand flow in 5 minutes
4. **No regressions**: All existing tests pass
5. **Same output**: Commands generated are identical to current system

---

## Timeline Estimate

| Phase | Effort | Risk |
|-------|--------|------|
| Phase 1: Create new module | 2-3 hours | Low |
| Phase 2: Compatibility layer | 1 hour | Low |
| Phase 3: Update graph_nodes | 1-2 hours | Medium |
| Phase 4: Remove old code | 1 hour | Medium |
| Phase 5: Testing & docs | 1-2 hours | Low |

**Total**: 6-9 hours

---

## Open Questions

1. **YAML format changes?** Current plan keeps YAML unchanged. Future improvement could simplify it.

2. **Logging strategy?** Should CommandBuilder have its own logger or use callback?

3. **Error handling?** Return None vs raise exceptions vs return Result object?

4. **Caching?** Should CommandBuilder cache program definitions?
