# YAML-Based Best Files Scoring - Implementation Plan

## Overview

Move the best files scoring configuration from hardcoded Python to YAML,
consistent with our "configuration over code" principle.

## Current State

Scoring is hardcoded in `agent/best_files_tracker.py`:
- `MODEL_STAGE_SCORES` dict
- `MAP_STAGE_SCORES` dict  
- `MTZ_STAGE_SCORES` dict
- `_score_model()`, `_score_map()`, `_score_mtz()` methods with formulas

## Target State

Scoring configuration in `knowledge/metrics.yaml`:
- Stage scores per category
- Metric formulas with parameters
- Special rules (MTZ locking)
- Python loads from YAML with hardcoded fallback

---

## Implementation Plan

### Phase 1: Add YAML Configuration

**File: `knowledge/metrics.yaml`**

Add new `best_files_scoring` section with:
- Categories: model, map, mtz, map_coefficients, sequence, ligand_cif
- Stage scores for each category
- Metric scoring formulas
- Special rules

**Formula types to support:**
- `linear` - Higher value is better (e.g., map_cc)
- `linear_inverse` - Lower value is better (e.g., r_free, clashscore)
- `boolean` - True/false flag (e.g., has_rfree_flags)

### Phase 2: Add YAML Loading to BestFilesTracker

**File: `agent/best_files_tracker.py`**

Changes:
1. Add `_load_scoring_config()` method
2. Add `_get_default_scoring()` method (current hardcoded values)
3. Add `_apply_formula()` method for generic formula evaluation
4. Modify `_calculate_score()` to use loaded config
5. Modify `_score_model()`, `_score_map()`, `_score_mtz()` to use generic scoring
6. Keep hardcoded fallback for backward compatibility

### Phase 3: Add Tests

**File: `tests/test_best_files_tracker.py`**

New tests:
1. `test_yaml_scoring_loads` - Verify YAML config loads
2. `test_yaml_stage_scores` - Verify stage scores from YAML
3. `test_yaml_metric_formulas` - Verify each formula type
4. `test_yaml_fallback_to_defaults` - Verify fallback works
5. `test_yaml_custom_category` - Verify extensibility

### Phase 4: Update Documentation

**File: `AGENT_LOGIC.md`**

Update Best Files Tracking section to reference YAML configuration.

### Phase 5: Remove Hardcoded Values

**File: `agent/best_files_tracker.py`**

After testing:
1. Remove `MODEL_STAGE_SCORES`, `MAP_STAGE_SCORES`, `MTZ_STAGE_SCORES` class variables
2. Keep `_get_default_scoring()` as fallback only
3. Simplify scoring methods to use generic `_calculate_score()`

---

## YAML Schema Design

```yaml
best_files_scoring:
  # Each category has stage_scores, metric_scores, and optional special_rules
  
  model:
    stage_scores:
      refined: 100
      rsr_output: 100
      autobuild_output: 80
      docked: 60
      processed_predicted: 50
      predicted: 40
      phaser_output: 30
      pdb: 10
      _default: 10  # Fallback for unknown stages
    
    metric_scores:
      r_free:
        max_points: 40
        formula: linear_inverse
        best_value: 0.20
        worst_value: 0.40
      
      map_cc:
        max_points: 30
        formula: linear
        best_value: 1.0
        worst_value: 0.0
      
      clashscore:
        max_points: 30
        formula: linear_inverse
        best_value: 0
        worst_value: 20

  map:
    stage_scores:
      optimized_full_map: 100
      sharpened: 90
      density_modified: 80
      full_map: 50
      map: 40
      half_map: 10
      _default: 40
    
    metric_scores:
      resolution:
        max_points: 30
        formula: linear_inverse
        best_value: 1.0
        worst_value: 4.0

  mtz:
    stage_scores:
      refined_mtz: 70
      original: 50
      mtz: 40
      _default: 40
    
    metric_scores:
      has_rfree_flags:
        max_points: 30
        formula: boolean
    
    special_rules:
      lock_on_rfree: true

  map_coefficients:
    stage_scores:
      refined: 80
      _default: 50

  sequence:
    stage_scores:
      _default: 50

  ligand_cif:
    stage_scores:
      _default: 50
```

---

## Formula Implementations

### linear
Higher value is better.
```
score = max_points × (value - worst_value) / (best_value - worst_value)
score = max(0, min(max_points, score))  # Clamp to [0, max_points]
```

Example: map_cc=0.75, best=1.0, worst=0.0, max=30
```
score = 30 × (0.75 - 0.0) / (1.0 - 0.0) = 22.5
```

### linear_inverse
Lower value is better.
```
score = max_points × (worst_value - value) / (worst_value - best_value)
score = max(0, min(max_points, score))
```

Example: r_free=0.25, best=0.20, worst=0.40, max=40
```
score = 40 × (0.40 - 0.25) / (0.40 - 0.20) = 40 × 0.15 / 0.20 = 30
```

### boolean
True gives max points, false gives 0.
```
score = max_points if value else 0
```

Example: has_rfree_flags=True, max=30
```
score = 30
```

---

## Test Cases

### test_yaml_scoring_loads
```python
def test_yaml_scoring_loads():
    tracker = BestFilesTracker()
    assert tracker.scoring is not None
    assert "model" in tracker.scoring
    assert "stage_scores" in tracker.scoring["model"]
```

### test_yaml_stage_scores
```python
def test_yaml_stage_scores():
    tracker = BestFilesTracker()
    # Refined model should get 100 points from stage
    tracker.evaluate_file("/path/refined.pdb", cycle=1, stage="refined")
    entry = tracker.get_best("model")
    assert entry.score == 100  # Stage only, no metrics
```

### test_yaml_metric_formulas
```python
def test_yaml_linear_formula():
    tracker = BestFilesTracker()
    # map_cc=0.75 should give ~22.5 points
    tracker.evaluate_file("/path/refined.pdb", cycle=1, stage="refined",
                         metrics={"map_cc": 0.75})
    entry = tracker.get_best("model")
    assert 122 < entry.score < 123  # 100 + 22.5

def test_yaml_linear_inverse_formula():
    tracker = BestFilesTracker()
    # r_free=0.25 should give 30 points
    tracker.evaluate_file("/path/refined.pdb", cycle=1, stage="refined",
                         metrics={"r_free": 0.25})
    entry = tracker.get_best("model")
    assert 129 < entry.score < 131  # 100 + 30

def test_yaml_boolean_formula():
    tracker = BestFilesTracker()
    # has_rfree_flags=True should give 30 points
    tracker.evaluate_file("/path/data.mtz", cycle=1, stage="original",
                         metrics={"has_rfree_flags": True})
    entry = tracker.get_best("mtz")
    assert entry.score == 80  # 50 + 30
```

### test_yaml_fallback_to_defaults
```python
def test_yaml_fallback_to_defaults():
    # Temporarily break YAML loading
    tracker = BestFilesTracker()
    tracker.scoring = None
    tracker._ensure_scoring_loaded()  # Should load defaults
    assert tracker.scoring is not None
```

---

## File Changes Summary

| File | Changes |
|------|---------|
| `knowledge/metrics.yaml` | Add `best_files_scoring` section |
| `agent/best_files_tracker.py` | Load from YAML, generic formula evaluation |
| `tests/test_best_files_tracker.py` | Add YAML-specific tests |
| `AGENT_LOGIC.md` | Update documentation |

---

## Implementation Order

1. **Add YAML config** to `knowledge/metrics.yaml`
2. **Add loading code** to `best_files_tracker.py`:
   - `_load_scoring_config()`
   - `_get_default_scoring()`
   - `_apply_formula()`
3. **Modify scoring methods** to use YAML config
4. **Run existing tests** - Should all pass (same values)
5. **Add new YAML-specific tests**
6. **Update documentation**
7. **Remove hardcoded class variables** (optional cleanup)

---

## Success Criteria

1. All existing tests pass (backward compatible)
2. New YAML-specific tests pass
3. Scoring results identical to hardcoded version
4. YAML config is human-readable and self-documenting
5. Fallback to defaults works if YAML missing/invalid
