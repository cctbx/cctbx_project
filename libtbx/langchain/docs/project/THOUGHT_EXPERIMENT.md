# Agent Decision Thought Experiment

This document traces through how the PHENIX AI Agent would handle standard crystallography cases, based on the YAML configurations and workflow engine logic.

**Updated for v110** - Includes dual MTZ tracking, stepwise mode, and error recovery. Code-level improvements in v111-v112 (metrics extraction, summary display) do not change workflow traces.

---

## Scenario 1: X-ray Structure Solution from Scratch (Automated Mode)

**Input:**
- `raster.mtz` - reflection data (classified as `data_mtz`)
- `raster.fa` - sequence file
- User directive: "solve the structure"
- Mode: `maximum_automation=True` (default)

**Assumptions:**
- Normal metrics throughout (no anomalous issues, no twinning)
- R-free improves steadily: ~0.45 → 0.35 → 0.28 → 0.25 → plateau at ~0.24
- Resolution: ~2.0 Å

### Cycle-by-Cycle Trace

#### Cycle 1: Data Analysis

**State Detection:**
- `_categorize_files()` → `{data_mtz: [raster.mtz], sequence: [raster.fa]}`
- `_detect_experiment_type()` → `xray` (has MTZ, no maps)
- `detect_phase()`:
  - `xtriage_done` = False (no history)
  - → **Phase: `analyze`**
- **State: `xray_initial`**
- **Valid programs: `[phenix.xtriage]`**

**Program Selection (LLM or rules):**
- Only one valid program → **`phenix.xtriage`**

**Command Built:**
```
phenix.xtriage raster.mtz
```

**Output Files:** (none - analysis only)  
**Metrics Extracted:** resolution=2.0, no twinning detected, no anomalous signal

**Best Files Tracker:**
- `data_mtz` locked to `raster.mtz` (first MTZ with R-free flags)

---

#### Cycle 2: Model Building with AlphaFold

**State Detection:**
- `xtriage_done` = True (from history)
- `has_placed_model` = False
- `has_predicted_model` = False
- `automation_path` = "automated" (from `maximum_automation=True`)
- → **Phase: `obtain_model`**
- **State: `xray_analyzed`**
- **Valid programs:** `[phenix.predict_and_build, phenix.phaser, phenix.autosol]`

**Conditions:**
- `phenix.predict_and_build`: requires sequence ✓
- `phenix.phaser`: requires search_model ✗
- `phenix.autosol`: requires sequence + has_anomalous ✗ (no anomalous signal)

**Program Selection:**
- `phenix.predict_and_build` is preferred (first in list, has sequence)
- **`phenix.predict_and_build`** selected

**Strategy from workflows.yaml:**
- `stop_after_predict: false` (automated mode runs full workflow)

**Command Built:**
```
phenix.predict_and_build input_files.seq_file=raster.fa \
    input_files.data_file=raster.mtz \
    crystal_info.resolution=2.0
```

**Output Files:** 
- `overall_best_model.pdb` (placed model)
- `overall_best_map_coeffs.mtz` (map coefficients)

**Metrics Extracted:** r_free=0.45, r_work=0.40

**Best Files Tracker:**
- `model` = `overall_best_model.pdb` (score: 100, predict_build_output)
- `map_coeffs_mtz` = `overall_best_map_coeffs.mtz` (most recent wins)

---

#### Cycle 3: First Refinement

**State Detection:**
- `has_placed_model` = True (overall_best_model.pdb categorized as `model`)
- `has_refined_model` = False
- → **Phase: `improve_model`**
- **Valid programs:** `[phenix.refine, phenix.autobuild, phenix.ligandfit, phenix.polder]`

**Conditions:**
- `phenix.autobuild`: r_free > 0.35 ✓, has sequence ✓, not_done ✓
- `phenix.ligandfit`: has ligand_file ✗
- `phenix.polder`: has ligand_file ✗

**Program Selection:**
- R-free is 0.45 (high), autobuild is available
- But `phenix.refine` is preferred for initial refinement
- **`phenix.refine`** selected

**File Selection (dual MTZ tracking):**
- Model: `overall_best_model.pdb` (best model)
- Data: `raster.mtz` (locked `data_mtz` - original with R-free flags)
- NOT `overall_best_map_coeffs.mtz` (that's map coefficients, not data)

**Command Built:**
```
phenix.refine overall_best_model.pdb raster.mtz output.prefix=refine_001
```

**Output Files:** 
- `refine_001_001.pdb` (refined model)
- `refine_001_001.mtz` (map coefficients)
- `refine_001_data.mtz` (copy of input data)

**Metrics Extracted:** r_free=0.35, r_work=0.30

**Best Files Tracker:**
- `model` = `refine_001_001.pdb` (score: 100, refined)
- `map_coeffs_mtz` = `refine_001_001.mtz` (most recent, better maps)

---

#### Cycles 4-6: Continue Refinement

Each cycle:
- Uses best model from previous cycle
- Uses locked `data_mtz` (raster.mtz) - NOT refined MTZ outputs
- Updates `map_coeffs_mtz` to most recent refine output

**Metrics progression:**
- Cycle 4: r_free=0.28
- Cycle 5: r_free=0.25
- Cycle 6: r_free=0.24 (plateauing)

---

#### Cycle 7: Validation

**State Detection:**
- R-free = 0.24 (at or below target for 2.0 Å resolution)
- `validation_done` = False
- → **Phase: `validate`**
- **Valid programs:** `[phenix.molprobity]`

**Program Selection:**
- **`phenix.molprobity`** selected

**Output:** clashscore=3.2, Ramachandran outliers=0.3%

---

#### Cycle 8: Complete

**State Detection:**
- R-free at target
- Validation done
- Quality acceptable
- → **Phase: `complete`**
- **STOP** with reason: "R-free target achieved, model quality acceptable"

---

### Summary: Scenario 1 (Automated)

| Cycle | Phase | Program | R-free | Best Model |
|-------|-------|---------|--------|------------|
| 1 | analyze | phenix.xtriage | - | - |
| 2 | obtain_model | phenix.predict_and_build | 0.45 | overall_best_model.pdb |
| 3 | improve_model | phenix.refine | 0.35 | refine_001_001.pdb |
| 4 | improve_model | phenix.refine | 0.28 | refine_002_001.pdb |
| 5 | improve_model | phenix.refine | 0.25 | refine_003_001.pdb |
| 6 | improve_model | phenix.refine | 0.24 | refine_004_001.pdb |
| 7 | validate | phenix.molprobity | - | - |
| 8 | complete | STOP | - | - |

**Total: 8 cycles**

---

## Scenario 2: X-ray with Stepwise Mode

**Input:** Same as Scenario 1
- Mode: `maximum_automation=False` (stepwise)

**Key Difference:** `predict_and_build` stops after prediction, giving user a checkpoint.

### Cycle-by-Cycle Trace

#### Cycle 1: Data Analysis (same as Scenario 1)

#### Cycle 2: Prediction Only

**State Detection:**
- `automation_path` = "stepwise"
- → Strategy includes `stop_after_predict: true`

**Command Built:**
```
phenix.predict_and_build input_files.seq_file=raster.fa \
    input_files.data_file=raster.mtz \
    stop_after_predict=true
```

**Output Files:**
- `predicted_model.pdb` (AlphaFold prediction, NOT placed)

**Note:** Model is predicted but not yet placed into density.

---

#### Cycle 3: Process Predicted Model

**State Detection:**
- `predict_done` = True
- `predict_full_done` = False (stopped after predict)
- `process_predicted_done` = False
- → **Valid programs:** `[phenix.process_predicted_model]`

**Command Built:**
```
phenix.process_predicted_model predicted_model.pdb
```

**Output:** Processed model ready for molecular replacement

---

#### Cycle 4: Molecular Replacement

**State Detection:**
- `process_predicted_done` = True
- `phaser_done` = False
- → **Valid programs:** `[phenix.phaser]`

**Command Built:**
```
phenix.phaser processed_model.pdb raster.mtz
```

**Output:** Placed model with initial phases

---

#### Cycles 5+: Refinement (same as Scenario 1)

### Summary: Scenario 2 (Stepwise)

| Cycle | Phase | Program | Notes |
|-------|-------|---------|-------|
| 1 | analyze | phenix.xtriage | |
| 2 | obtain_model | phenix.predict_and_build | stop_after_predict=true |
| 3 | obtain_model | phenix.process_predicted_model | Prepare for MR |
| 4 | obtain_model | phenix.phaser | Place model |
| 5+ | improve_model | phenix.refine | Multiple cycles |
| N | validate | phenix.molprobity | |
| N+1 | complete | STOP | |

**Total: ~10-12 cycles** (more granular control)

---

## Scenario 3: X-ray with Ligand Fitting

**Input:**
- `raster.mtz` - reflection data
- `raster.fa` - sequence file
- `m29lig.pdb` - ligand structure
- User directive: "solve the structure and fit the ligand"

### Key Decision Points

#### When is ligandfit triggered?

From `workflows.yaml`:
```yaml
phenix.ligandfit:
  conditions:
    has: ligand
    not_done: ligandfit
    metric: r_free
    less_than: 0.35
```

Ligand fitting waits until:
- R-free < 0.35 (model quality threshold)
- Ligand file is present
- ligandfit hasn't been run yet

This prevents fitting ligand into poor density.

#### File Selection for ligandfit

**Critical:** `phenix.ligandfit` needs map coefficients (phases), not data!

From `programs.yaml`:
```yaml
phenix.ligandfit:
  input_files:
    required:
      - name: map_coeffs_mtz   # Uses map_coeffs_mtz category
        slot: map_coeffs_file
      - name: model
        slot: model_file
      - name: ligand
        slot: ligand_file
```

**Best Files Tracker provides:**
- `model` = `refine_002_001.pdb` (best refined model)
- `map_coeffs_mtz` = `refine_002_001.mtz` (latest map coefficients)
- `ligand` = `m29lig.pdb`

**Command Built:**
```
phenix.ligandfit model=refine_002_001.pdb \
    map_coeffs_file=refine_002_001.mtz \
    ligand=m29lig.pdb
```

### Summary: Scenario 3

| Cycle | Phase | Program | Notes |
|-------|-------|---------|-------|
| 1 | analyze | phenix.xtriage | |
| 2 | obtain_model | phenix.predict_and_build | |
| 3 | improve_model | phenix.refine | R-free 0.45→0.35 |
| 4 | improve_model | phenix.refine | R-free 0.35→0.28 |
| 5 | improve_model | phenix.ligandfit | R-free < 0.35, fit ligand |
| 6 | improve_model | phenix.pdbtools | Combine protein + ligand |
| 7 | improve_model | phenix.refine | Refine combined |
| 8+ | improve_model | phenix.refine | Continue to convergence |
| N | validate | phenix.molprobity | |
| N+1 | complete | STOP | |

---

## Scenario 4: SAD Phasing with Anomalous Data

**Input:**
- `sad_data.mtz` - reflection data with anomalous signal
- `sequence.fa` - sequence file
- User directive: "use experimental phasing" OR anomalous signal detected

### Cycle-by-Cycle Trace

#### Cycle 1: Data Analysis

**Metrics Extracted from xtriage:**
```
anomalous_measurability: 0.15
anomalous_resolution: 3.5
has_anomalous: true
```

**Best Files Tracker:**
- Stores anomalous metrics in session

---

#### Cycle 2: Experimental Phasing

**State Detection:**
- `has_anomalous` = True (from xtriage metrics)
- `has_sequence` = True
- `autosol_done` = False
- → **Valid programs include:** `phenix.autosol`

**Directive Priority:**
If user said "use experimental phasing":
- `autosol` moved to front of valid programs list
- `predict_and_build` moved to end

**Program Selection:**
- **`phenix.autosol`** selected (anomalous data + directive)

**Command Built:**
```
phenix.autosol data=sad_data.mtz seq_file=sequence.fa
```

---

#### Potential Error Recovery

If MTZ has both merged and anomalous arrays:

**Cycle 2 Attempt 1:** autosol fails
```
Multiple equally suitable arrays of observed xray data found.
Choices: IMEAN, I(+)/I(-)
```

**Error Recovery:**
1. ErrorAnalyzer detects ambiguous data error
2. Context: program=autosol, has_anomalous=true
3. Selection: `I(+)/I(-)` (anomalous needed for SAD)
4. Saves recovery strategy to session

**Cycle 3:** autosol retries with fix
```
phenix.autosol data=sad_data.mtz seq_file=sequence.fa \
    scaling.input.xray_data.obs_labels="I(+)"
```
→ SUCCESS

---

## Scenario 5: Session Restart After Partial Workflow

**Situation:** User ran xtriage, detected anomalous signal, then stopped. Now restarting.

### Key: Metrics Preserved Across Restart

When session loads:
1. `get_history_for_agent()` returns previous cycles
2. Each cycle includes `analysis` dict with stored metrics
3. `_analyze_history()` extracts `has_anomalous=True` from xtriage cycle

**State on Restart:**
```python
history_info = {
    "xtriage_done": True,
    "has_anomalous": True,
    "anomalous_measurability": 0.15,
    "strong_anomalous": True,  # measurability > 0.10
}
```

**Valid Programs:**
- `phenix.autosol` is available (has_anomalous=True, has_sequence=True)
- Even though xtriage isn't re-run, its metrics persist

---

## Key Decision Points

### 1. Experiment Type Detection

Based on file categories:
- Has `.mtz` + no maps → **X-ray**
- Has `.mrc`/`.ccp4` maps → **Cryo-EM**

### 2. Dual MTZ Tracking (v110)

| Category | Update Rule | Used By |
|----------|-------------|---------|
| `data_mtz` | First with R-free **locks forever** | phenix.refine, phenix.phaser |
| `map_coeffs_mtz` | **Most recent wins** | phenix.ligandfit |

This prevents refinement from using its own output maps as input data.

### 3. Automation Path (v110)

| Mode | `maximum_automation` | predict_and_build Behavior |
|------|---------------------|---------------------------|
| Automated | True (default) | Runs full workflow |
| Stepwise | False | Stops after prediction |

### 4. Anomalous Workflow Triggering

`phenix.autosol` becomes available when:
- `has_anomalous` = True (from xtriage or history)
- `has_sequence` = True
- `autosol_done` = False

The `use_experimental_phasing` directive prioritizes autosol over predict_and_build.

### 5. Plateau Detection

From `workflows.yaml`:
```yaml
repeat:
  max_cycles: 4
  until:
    any:
      - r_free: "< target_r_free"
      - condition: plateau
        cycles: 2
        threshold: 0.005
```

Agent stops refinement when:
- Target R-free reached, OR
- 2 consecutive cycles with < 0.005 improvement

### 6. File Selection Logic

**Best model selection:**
1. Check `best_files["model"]` from BestFilesTracker
2. Prefer by stage score: refined (100) > autobuild_output (100) > predict_build_output (90)
3. Among equal scores, prefer better R-free

**MTZ selection:**
- For refinement: Use locked `data_mtz` (original with R-free flags)
- For ligandfit: Use `map_coeffs_mtz` (latest calculated phases)

---

## Edge Cases

### What if predict_and_build fails?

The agent would:
1. See error in log
2. If anomalous signal detected: try `phenix.autosol`
3. If template available: try `phenix.phaser`
4. Eventually STOP with error if no path forward

### What if user provides a pre-built model?

If `model.pdb` is provided alongside `data.mtz`:
- File categorization checks filename patterns
- If looks refined: skip to refinement phase
- If looks like template: use phaser to place it

### What if ambiguous data arrays?

Error recovery system:
1. Detects "Multiple equally suitable arrays" error
2. Selects appropriate array based on context (merged vs anomalous)
3. Retries program with explicit label selection
4. Max 3 retries before giving up

### What if ligand file is named ambiguously?

From `file_categories.yaml`:
```yaml
ligand:
  extensions: [.pdb, .cif]
  prefer_patterns: [lig, ligand]
```

Files matching `lig*` or `ligand*` are categorized as ligands, not protein models.
