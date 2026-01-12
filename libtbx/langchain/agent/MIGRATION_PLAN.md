# Migration Plan: YAML-Driven Agent Architecture

## Vision

Transform the PHENIX AI Agent from a code-centric system to a configuration-driven system where:
- **Crystallography knowledge is data, not code**
- **A crystallographer can modify behavior by editing YAML files**
- **Python code is a generic workflow engine with no domain-specific hardcoding**

## Current Architecture Problems

### Knowledge is Scattered
To understand or modify "what happens when", you must look at:
- `workflow_state.py` (800+ lines) - state detection, validation logic
- `graph_nodes.py` (1000+ lines) - LLM calls, command building
- `decision_config.json` - program rankings, some conditions
- `command_templates.json` - command patterns
- `prompts_hybrid.py` - program descriptions, LLM instructions
- `log_parsers.py` - metric extraction patterns

### Code Mixes Infrastructure and Domain Knowledge
Functions like `_detect_xray_state()` combine:
- Generic logic (check what files exist)
- Crystallography knowledge (xtriage should run before phaser)

### Hard to Test and Teach
- Can't easily run "what if" scenarios
- LLM makes decisions opaque
- No clear separation of "what's possible" vs "what's preferred"

## Target Architecture

### Three YAML Files Define Everything

```
/knowledge/
  programs.yaml     # What tools exist, their inputs/outputs
  workflows.yaml    # What should happen when (state machine)
  metrics.yaml      # What metrics exist, how to evaluate them
```

### Python Code Becomes Generic Engine

```
/agent/
  engine.py         # Workflow interpreter (reads YAML, executes)
  file_manager.py   # Track available files
  metric_tracker.py # Extract and track metrics
  llm_selector.py   # Optional LLM for choosing among valid programs
  command_builder.py # Build commands from templates
```

## File Specifications

### 1. programs.yaml

Defines every program the agent can use:

```yaml
# programs.yaml
phenix.xtriage:
  description: "Analyze X-ray data quality, detect twinning and anomalous signal"
  category: analysis
  experiment_types: [xray]
  
  inputs:
    required:
      - type: mtz
    optional:
      - type: sequence
  
  outputs:
    files:
      - pattern: "xtriage_*.log"
        type: log
    metrics:
      - resolution
      - completeness
      - twinning
      - anomalous_signal
      - space_group
  
  command: "phenix.xtriage {mtz}"
  
  log_parsing:
    resolution:
      pattern: 'High resolution limit\s*:\s*([0-9.]+)'
      type: float
    twinning:
      pattern: 'Twin fraction\s*:\s*([0-9.]+)'
      type: float
    anomalous_signal:
      pattern: 'Anomalous signal.*?:\s*([0-9.]+)'
      type: float

phenix.refine:
  description: "Refine atomic model against X-ray data"
  category: refinement
  experiment_types: [xray]
  
  inputs:
    required:
      - type: mtz
      - type: model
    optional:
      - type: cif  # ligand restraints
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
      - pattern: "*_refine_*.mtz"
        type: mtz
      - pattern: "*_refine_*.log"
        type: log
    metrics:
      - r_work
      - r_free
      - bonds_rmsd
      - angles_rmsd
  
  command: "phenix.refine {model} {mtz} output.prefix={prefix}"
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
    r_work:
      pattern: 'R-work\s*[=:]\s*([0-9.]+)'
      type: float

phenix.real_space_refine:
  description: "Refine atomic model against cryo-EM map"
  category: refinement
  experiment_types: [cryoem]
  
  inputs:
    required:
      - type: map
      - type: model
    optional:
      - type: half_maps
      - type: resolution  # from session
  
  outputs:
    files:
      - pattern: "*_real_space_refined*.pdb"
        type: model
    metrics:
      - map_cc
      - clashscore
  
  command: "phenix.real_space_refine {model} {map} resolution={resolution}"

# ... more programs
```

### 2. workflows.yaml

Defines the state machine for each experiment type:

```yaml
# workflows.yaml
xray:
  description: "X-ray crystallography structure determination"
  
  phases:
    analyze:
      description: "Analyze data quality"
      goal: "Understand data characteristics before proceeding"
      programs:
        - xtriage
      transitions:
        default: obtain_model
    
    obtain_model:
      description: "Get an initial model"
      goal: "Obtain atomic coordinates to refine"
      programs:
        - predict_and_build:
            preferred: true
            conditions:
              - has: sequence
              - has: mtz
        - phaser:
            conditions:
              - has: search_model
              - has: mtz
      transitions:
        default: refine
        if_program: predict_and_build
        and_condition: stop_after_predict
        then: molecular_replacement
    
    molecular_replacement:
      description: "Place predicted model in unit cell"
      goal: "Position model using diffraction data"
      programs:
        - phaser:
            conditions:
              - has: prediction
              - has: mtz
      transitions:
        default: refine
    
    refine:
      description: "Improve model against data"
      goal: "Reduce R-free to target"
      programs:
        - refine
      repeat:
        max_cycles: 5
        until:
          any:
            - metric: r_free
              operator: "<"
              value: target_r_free
            - condition: no_improvement
              cycles: 2
      transitions:
        default: validate
    
    validate:
      description: "Check model quality"
      goal: "Ensure model meets quality standards"
      programs:
        - molprobity
      transitions:
        if_quality: acceptable
        then: complete
        else: refine
    
    complete:
      description: "Workflow finished"
      stop: true
      reason: "Structure determination complete"
  
  # Quality targets (can be overridden by user)
  targets:
    r_free:
      default: 0.25
      by_resolution:
        - range: [0, 1.5]
          value: 0.20
        - range: [1.5, 2.5]
          value: 0.25
        - range: [2.5, 3.5]
          value: 0.30
        - range: [3.5, 999]
          value: 0.35
    clashscore:
      default: 5
    ramachandran_outliers:
      default: 0.5


cryoem:
  description: "Cryo-EM structure determination"
  
  phases:
    analyze:
      description: "Analyze map quality"
      programs:
        - mtriage
      transitions:
        default: obtain_model
    
    obtain_model:
      description: "Build or dock initial model"
      programs:
        - predict_and_build:
            preferred: true
            conditions:
              - has: sequence
              - has: map
        - dock_in_map:
            conditions:
              - has: search_model
              - has: map
      transitions:
        default: refine
    
    optimize_map:
      description: "Improve map quality"
      optional: true
      programs:
        - resolve_cryo_em
        - map_sharpening
      transitions:
        default: refine
    
    refine:
      description: "Refine model against map"
      programs:
        - real_space_refine:
            requires_resolution: true
      repeat:
        max_cycles: 5
        until:
          any:
            - metric: map_cc
              operator: ">"
              value: 0.8
            - condition: no_improvement
              cycles: 2
      transitions:
        default: validate
    
    validate:
      description: "Check model quality"
      programs:
        - molprobity
      transitions:
        if_quality: acceptable
        then: complete
        else: refine
    
    complete:
      stop: true
      reason: "Cryo-EM structure determination complete"
  
  targets:
    map_cc:
      default: 0.75
    clashscore:
      default: 5
```

### 3. metrics.yaml

Defines metrics, their meaning, and how to evaluate them:

```yaml
# metrics.yaml
metrics:
  r_free:
    description: "Cross-validation R-factor (lower is better)"
    type: float
    direction: minimize
    good: 0.25
    acceptable: 0.30
    poor: 0.35
    format: "{:.4f}"
  
  r_work:
    description: "Working R-factor (lower is better)"
    type: float
    direction: minimize
    good: 0.20
    acceptable: 0.25
  
  resolution:
    description: "High resolution limit in Angstroms"
    type: float
    direction: minimize  # lower is better resolution
    unit: "Å"
    format: "{:.2f} Å"
  
  map_cc:
    description: "Map-model correlation coefficient"
    type: float
    direction: maximize
    good: 0.80
    acceptable: 0.70
    poor: 0.60
  
  clashscore:
    description: "Steric clashes per 1000 atoms (lower is better)"
    type: float
    direction: minimize
    good: 2
    acceptable: 5
    poor: 10
  
  ramachandran_outliers:
    description: "Percentage of residues in disallowed regions"
    type: float
    direction: minimize
    unit: "%"
    good: 0.2
    acceptable: 0.5
    poor: 1.0
  
  completeness:
    description: "Data completeness"
    type: float
    direction: maximize
    unit: "%"
    good: 99
    acceptable: 95
  
  twinning:
    description: "Twin fraction (0 = no twinning)"
    type: float
    direction: minimize
    threshold: 0.1  # above this suggests twinning

# Composite quality assessment
quality_rules:
  acceptable:
    all:
      - metric: clashscore
        operator: "<="
        value: acceptable
      - metric: ramachandran_outliers
        operator: "<="
        value: acceptable
    any:
      - metric: r_free
        operator: "<="
        value: acceptable
      - metric: map_cc
        operator: ">="
        value: acceptable

# Improvement detection
improvement:
  significant_change:
    r_free: 0.005
    map_cc: 0.01
    clashscore: 1.0
```

## Implementation Plan

### Phase 1: Create YAML Infrastructure (Foundation)

**Goal:** Create the YAML files and loaders without changing existing behavior.

**Tasks:**
1. Create `knowledge/programs.yaml` with all current programs
2. Create `knowledge/workflows.yaml` with current X-ray and cryo-EM workflows
3. Create `knowledge/metrics.yaml` with all tracked metrics
4. Create `agent/yaml_loader.py` to read and validate YAML files
5. Write tests that verify YAML matches current hardcoded behavior

**Deliverables:**
- Three complete YAML files
- Loader with validation
- Test suite comparing YAML to current behavior

**Estimated effort:** Medium

### Phase 2: Migrate Program Definitions

**Goal:** Programs are defined only in `programs.yaml`.

**Tasks:**
1. Create `agent/program_registry.py` that reads from YAML
2. Migrate command building to use YAML templates
3. Migrate log parsing patterns to use YAML definitions
4. Remove hardcoded program lists from `workflow_state.py`
5. Remove `command_templates.json` (absorbed into `programs.yaml`)
6. Update `prompts_hybrid.py` to generate program descriptions from YAML

**Deliverables:**
- Single source of truth for programs
- `command_templates.json` deleted
- Tests pass

**Estimated effort:** Medium-High

### Phase 3: Migrate Workflow State Machine

**Goal:** Workflow logic comes from `workflows.yaml`, not Python.

**Tasks:**
1. Create `agent/workflow_engine.py` - generic workflow interpreter
2. Replace `_detect_xray_state()` with workflow engine
3. Replace `_detect_cryoem_state()` with workflow engine
4. Replace validation logic with YAML-driven rules
5. Implement phase transitions from YAML
6. Implement repeat/until logic from YAML

**Deliverables:**
- `workflow_state.py` reduced to thin wrapper or deleted
- All workflow logic in YAML
- Same behavior, different implementation

**Estimated effort:** High

### Phase 4: Migrate Metrics System

**Goal:** Metrics defined only in `metrics.yaml`.

**Tasks:**
1. Create `agent/metric_evaluator.py` reading from YAML
2. Migrate quality assessment to use YAML rules
3. Migrate improvement detection to use YAML thresholds
4. Migrate target calculations (resolution-dependent R-free)
5. Remove hardcoded metric definitions from `metrics_analyzer.py`

**Deliverables:**
- Single source of truth for metrics
- Dynamic quality assessment from YAML
- `metrics_analyzer.py` simplified

**Estimated effort:** Medium

### Phase 5: Simplify LLM Integration

**Goal:** LLM is optional selector, not decision-maker.

**Tasks:**
1. Refactor so workflow engine determines valid programs
2. LLM only chooses among valid options (or uses rules)
3. Add "rules-only" mode that doesn't need LLM
4. Generate LLM prompts from YAML (program descriptions, current phase)
5. Make decision logging transparent ("Phase: X, Valid: [A,B], Selected: A because Y")

**Deliverables:**
- Agent works without LLM (deterministic mode)
- Clear separation of "what's valid" vs "what's preferred"
- Better decision transparency

**Estimated effort:** Medium

### Phase 6: Polish and Documentation

**Goal:** System is teachable and maintainable.

**Tasks:**
1. Write user guide for modifying YAML files
2. Create example modifications (add program, change threshold)
3. Add YAML validation with helpful error messages
4. Create visualization of workflow (generate Mermaid diagram from YAML)
5. Add simulation/dry-run mode
6. Create test scenarios with expected outcomes

**Deliverables:**
- Complete documentation
- Teaching materials
- Robust validation
- Optional visualization

**Estimated effort:** Medium

## Migration Strategy

### Parallel Implementation
- Keep existing code working throughout
- New YAML-driven code alongside old code
- Feature flag to switch between implementations
- Gradual cutover as each phase completes

### Testing Strategy
- Create "behavior snapshots" of current system
- After each phase, verify behavior matches snapshots
- Add new tests for YAML-specific features

### Rollback Plan
- Each phase is independent
- Can stop migration at any phase
- Old code remains until new code is proven

## Success Criteria

After migration, these tasks should require only YAML edits:

1. ✅ Add a new program to the workflow
2. ✅ Change when the workflow stops (different R-free target)
3. ✅ Add a new metric to track
4. ✅ Change program preferences (prefer phaser over predict_and_build)
5. ✅ Add a condition (if twinned, suggest twin refinement)
6. ✅ Create a simplified workflow (analysis only)
7. ✅ Adjust quality thresholds

And these should NOT require touching YAML:

1. ✅ Core graph execution
2. ✅ File tracking
3. ✅ Session management
4. ✅ Server communication
5. ✅ LLM API calls

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| YAML becomes too complex | Keep hierarchical, add validation, good defaults |
| Performance regression | Profile after each phase, optimize YAML loading |
| Edge cases not covered by YAML | Allow "escape hatch" Python plugins |
| Migration takes too long | Each phase delivers value independently |
| Breaking existing users | Feature flag, extensive testing |

## Timeline Estimate

| Phase | Effort | Dependencies |
|-------|--------|--------------|
| Phase 1: YAML Infrastructure | 1-2 days | None |
| Phase 2: Program Definitions | 2-3 days | Phase 1 |
| Phase 3: Workflow Engine | 3-4 days | Phase 1, 2 |
| Phase 4: Metrics System | 1-2 days | Phase 1 |
| Phase 5: LLM Simplification | 2-3 days | Phase 3 |
| Phase 6: Polish | 2-3 days | All above |

**Total: ~12-17 days of focused work**

## Next Steps

1. Review and approve this plan
2. Start Phase 1: Create initial YAML files
3. Iterate based on what we learn

---

*This document is a living plan. Update as implementation progresses.*
