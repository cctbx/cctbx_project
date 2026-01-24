# Future Plans and Roadmap

This document describes planned features that have been designed but not yet implemented.

For documentation of implemented features, see:
- [ARCHITECTURE.md](../reference/ARCHITECTURE.md) - System design and components
- [USER_DIRECTIVES.md](../guides/USER_DIRECTIVES.md) - User directive system
- [VALIDATION.md](../reference/VALIDATION.md) - Sanity checking and validation gates
- [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) - Centralized YAML architecture (completed)
- [ADDING_PROGRAMS.md](../guides/ADDING_PROGRAMS.md) - How to add new programs

---

## Recently Completed Features ✅

The following features have been implemented:

| Feature | Status | Description |
|---------|--------|-------------|
| Centralized Metric Patterns | ✅ Done | Patterns in programs.yaml used by log_parsers.py and session.py |
| Auto-Registration of Programs | ✅ Done | `run_once: true` auto-generates tracking flags |
| Declarative Summary Display | ✅ Done | Summary display configured in metrics.yaml |
| Mid-Session Directive Updates | ✅ Done | Advice changes detected and re-processed |

See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for implementation details.

---

## 1. YAML-Based Best Files Scoring

**Status:** Planned  
**Priority:** Medium  
**Effort:** Medium

### Overview

Move the best files scoring configuration from hardcoded Python to YAML,
consistent with our "configuration over code" principle.

### Current State

Scoring is mostly defined in `knowledge/metrics.yaml` under `best_files_scoring`,
but some stage scores are still in `agent/best_files_tracker.py`.

### Target State

Complete migration to YAML:

```yaml
best_files_scoring:
  model:
    stage_scores:
      refined: 100
      autobuild_output: 80
      docked: 60
      phaser_output: 30
      pdb: 10
    
    metric_scores:
      r_free:
        max_points: 40
        formula: linear_inverse
        best_value: 0.20
        worst_value: 0.40
```

### Implementation Steps

1. Verify all scoring is in `knowledge/metrics.yaml`
2. Update `BestFilesTracker` to load exclusively from YAML
3. Remove any remaining hardcoded scores
4. Add tests for YAML-based scoring

---

## 2. Interactive Directive Refinement

**Status:** Concept  
**Priority:** Low  
**Effort:** Medium

### Overview

After extracting directives from user advice, present them to the user for confirmation:

```
Extracted Directives:
  Program Settings:
    phenix.autosol: resolution=3.0
    phenix.refine: anisotropic_adp=True
  
  Stop Conditions:
    after_program: phenix.refine
    skip_validation: True

Is this correct? [Y/n/edit]
```

### Benefits

- Catches misinterpretations early
- Allows users to fine-tune without rewriting advice
- Builds trust in the directive system

### Considerations

- Requires interactive mode support
- Need good default behavior for non-interactive runs
- UI design for editing directives

---

## 3. Directive Learning

**Status:** Concept  
**Priority:** Low  
**Effort:** High

### Overview

Log which directives lead to successful outcomes, building a knowledge base:

```json
{
  "scenario": "SAD_Se_high_resolution",
  "success_rate": 0.85,
  "common_directives": {
    "program_settings": {
      "phenix.autosol": {"resolution": "data_resolution * 1.2"}
    }
  }
}
```

### Use Cases

- Suggest directives for similar datasets
- Warn when directives conflict with successful patterns
- Generate scenario-specific templates

### Implementation Considerations

- Privacy: User data handling
- Storage: Where to keep learned patterns
- Feedback loop: How to measure "success"

---

## 4. Directive Conflict Detection

**Status:** Concept  
**Priority:** Low  
**Effort:** Medium

### Overview

Warn when directives may conflict:

```
Warning: Potential directive conflict detected:
  - stop_conditions.after_cycle: 3
  - stop_conditions.r_free_target: 0.20

R-free target of 0.20 is unlikely to be reached in only 3 cycles.
Proceed anyway? [y/N]
```

### Conflict Types to Detect

- Unrealistic metric targets given cycle limits
- Contradictory workflow preferences (e.g., skip autobuild but require model)
- Resolution mismatches between programs
- Stop conditions that may never trigger

---

## 5. Self-Improving Tips System

**Status:** Concept  
**Priority:** Medium  
**Effort:** High

### Overview

Automatically generate tips and rules from failed runs. When the agent
encounters errors or suboptimal outcomes, extract lessons learned:

```yaml
# Auto-generated tip
tips:
  - condition: "autosol failed with 'insufficient anomalous signal'"
    action: "Check anomalous_measurability > 0.05 before running autosol"
    source: "auto-generated from session abc123"
```

### Implementation Considerations

- Error pattern recognition
- Tip validation (avoid false correlations)
- Storage and retrieval system
- Integration with decision-making

---

## Implementation Priority Summary

| Feature | Priority | Effort | Status |
|---------|----------|--------|--------|
| YAML-based best files scoring | Medium | Medium | Planned |
| Interactive directive refinement | Low | Medium | Concept |
| Directive learning | Low | High | Concept |
| Directive conflict detection | Low | Medium | Concept |
| Self-improving tips system | Medium | High | Concept |
