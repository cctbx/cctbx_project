# Agent Decision Thought Experiment

This document traces through how the PHENIX AI Agent would handle standard crystallography cases, based on the YAML configurations and workflow engine logic.

---

## Scenario 1: X-ray Structure Solution from Scratch

**Input:**
- `raster.mtz` - reflection data
- `raster.fa` - sequence file
- User directive: "solve the structure"

**Assumptions:**
- Normal metrics throughout (no anomalous issues, no twinning)
- R-free improves steadily: ~0.45 → 0.35 → 0.28 → 0.25 → plateau at ~0.24
- Resolution: ~2.0 Å

### Cycle-by-Cycle Trace

#### Cycle 1: Data Analysis

**State Detection:**
- `_categorize_files()` → `{mtz: [raster.mtz], sequence: [raster.fa]}`
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
**Metrics Extracted:** resolution=2.0, no twinning detected

---

#### Cycle 2: Model Building with AlphaFold

**State Detection:**
- `xtriage_done` = True (from history)
- `has_placed_model` = False
- `has_predicted_model` = False
- → **Phase: `obtain_model`**
- **State: `xray_analyzed`**
- **Valid programs:** `[phenix.predict_and_build, phenix.phaser, phenix.autosol]`

**Conditions:**
- `phenix.predict_and_build`: requires sequence ✓
- `phenix.phaser`: requires search_model ✗
- `phenix.autosol`: requires sequence + anomalous ✗

**Program Selection:**
- `phenix.predict_and_build` is preferred (first in list)
- Has sequence, has MTZ data
- **`phenix.predict_and_build`** selected

**Strategy from workflows.yaml:**
- `stop_after_predict: false` (default for X-ray automated path)

**Command Built:**
```
phenix.predict_and_build input_files.seq_file=raster.fa input_files.xray_data_file=raster.mtz crystal_info.resolution=2.0
```

**Output Files:** `predict_and_build_001_overall_best_model.pdb`  
**Metrics Extracted:** r_free=0.45, r_work=0.40

---

#### Cycle 3: First Refinement

**State Detection:**
- `has_placed_model` = True (overall_best_model.pdb categorized as `model`)
- `has_refined_model` = False
- → **Phase: `refine`**
- **State: `xray_refined`** (misleading name - means "refining state")
- **Valid programs:** `[phenix.refine, phenix.autobuild, phenix.ligandfit]`

**Conditions:**
- `phenix.autobuild`: r_free > 0.35 ✓, has sequence ✓, not_done ✓
- `phenix.ligandfit`: has ligand_file ✗

**Program Selection:**
- R-free is 0.45 (high), could try autobuild
- But `phenix.refine` is preferred
- **`phenix.refine`** selected (standard approach first)

**Command Built:**
```
phenix.refine predict_and_build_001_overall_best_model.pdb raster.mtz output.prefix=refine_001
```

**Output Files:** `refine_001_001.pdb`, `refine_001_001.mtz`  
**Metrics Extracted:** r_free=0.35, r_work=0.30

---

#### Cycle 4: Second Refinement

**State Detection:**
- `has_refined_model` = True (refine_001_001.pdb exists)
- `r_free` = 0.35
- Not at target (target ~0.25 for 2.0 Å)
- → **Phase: `refine`**
- **Valid programs:** `[phenix.refine, phenix.autobuild, phenix.ligandfit]`

**Program Selection:**
- R-free still above autobuild threshold (0.35)
- Could try autobuild, but let's say LLM chooses to continue refine
- **`phenix.refine`** selected

**Command Built:**
```
phenix.refine refine_001_001.pdb raster.mtz output.prefix=refine_002
```

**Output Files:** `refine_002_001.pdb`  
**Metrics Extracted:** r_free=0.28, r_work=0.24

---

#### Cycle 5: Third Refinement

**State Detection:**
- R-free = 0.28, improving
- Not at target yet
- → **Phase: `refine`**

**Program Selection:**
- Autobuild now excluded (r_free < 0.35)
- **`phenix.refine`** selected

**Metrics after:** r_free=0.25

---

#### Cycle 6: Fourth Refinement

**Metrics after:** r_free=0.24 (plateauing)

---

#### Cycle 7: Target Reached → Validation

**State Detection:**
- R-free = 0.24 (at or below target for 2.0 Å resolution)
- `validation_done` = False
- → **Phase: `validate`**
- **Valid programs:** `[phenix.molprobity, phenix.model_vs_data]`

**Program Selection:**
- `phenix.molprobity` is required
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

### Summary: Scenario 1

| Cycle | State | Program | R-free |
|-------|-------|---------|--------|
| 1 | xray_initial | phenix.xtriage | - |
| 2 | xray_analyzed | phenix.predict_and_build | 0.45 |
| 3 | xray_refined | phenix.refine | 0.35 |
| 4 | xray_refined | phenix.refine | 0.28 |
| 5 | xray_refined | phenix.refine | 0.25 |
| 6 | xray_refined | phenix.refine | 0.24 |
| 7 | xray_refined | phenix.molprobity | - |
| 8 | complete | STOP | - |

**Total: 8 cycles**

---

## Scenario 2: X-ray with Ligand Fitting

**Input:**
- `raster.mtz` - reflection data
- `raster.fa` - sequence file
- `m29lig.pdb` - ligand structure
- User directive: "solve the structure and fit the ligand"

**Key Difference:** Agent needs to recognize the ligand file and fit it after refinement reaches acceptable quality.

### Cycle-by-Cycle Trace

#### Cycles 1-5: Same as Scenario 1

The workflow proceeds identically until the model is reasonably refined (R-free ~0.28).

#### Cycle 6: Ligand Fitting (when R-free < 0.35)

**State Detection:**
- `has_refined_model` = True
- `has_ligand_file` = True (m29lig.pdb categorized as `ligand`)
- R-free = 0.28 (below 0.35 threshold)
- `ligandfit_done` = False
- → **Phase: `refine`**
- **Valid programs:** `[phenix.refine, phenix.ligandfit]`

**Conditions:**
- `phenix.ligandfit`: has ligand_file ✓, not_done ✓, r_free < 0.35 ✓

**Program Selection:**
- LLM sees user advice "fit the ligand"
- `phenix.ligandfit` conditions are met
- **`phenix.ligandfit`** selected

**File Selection:**
- `model`: refine_002_001.pdb (best refined model)
- `mtz`: raster.mtz (original data with R-free flags)
- `ligand`: m29lig.pdb

**Command Built:**
```
phenix.ligandfit model=refine_002_001.pdb data=raster.mtz ligand=m29lig.pdb
```

**Output Files:** `ligand_fit_001.pdb`

---

#### Cycle 7: Combine Ligand with Protein

**State Detection:**
- `has_ligand_fit` = True (ligand_fit_001.pdb exists)
- `pdbtools_done` = False
- → **Phase: `combine_ligand`**
- **Valid programs:** `[phenix.pdbtools]`

**Program Selection:**
- **`phenix.pdbtools`** selected

**File Selection:**
- `protein`: refine_002_001.pdb (refined protein, NOT ligand file)
- `ligand`: ligand_fit_001.pdb (fitted ligand coordinates)

**Command Built:**
```
phenix.pdbtools refine_002_001.pdb ligand_fit_001.pdb output.file_name=protein_with_ligand.pdb
```

**Output Files:** `protein_with_ligand.pdb`

---

#### Cycle 8: Refine with Ligand

**State Detection:**
- `pdbtools_done` = True
- Back to refine phase with combined model
- → **Phase: `refine`**

**Program Selection:**
- **`phenix.refine`** selected

**File Selection:**
- `model`: protein_with_ligand.pdb (combined model)
- `mtz`: raster.mtz

**Command Built:**
```
phenix.refine protein_with_ligand.pdb raster.mtz output.prefix=refine_003
```

---

#### Cycles 9-10: Continue Refinement

Continue until R-free plateaus (~0.22 with ligand)

---

#### Cycle 11: Validation

**Program:** `phenix.molprobity`

---

#### Cycle 12: Complete

**STOP** with reason: "R-free target achieved, ligand fitted, model quality acceptable"

---

### Summary: Scenario 2

| Cycle | State | Program | Notes |
|-------|-------|---------|-------|
| 1 | xray_initial | phenix.xtriage | Data analysis |
| 2 | xray_analyzed | phenix.predict_and_build | AlphaFold + build |
| 3 | xray_refined | phenix.refine | R-free 0.45→0.35 |
| 4 | xray_refined | phenix.refine | R-free 0.35→0.28 |
| 5 | xray_refined | phenix.refine | R-free 0.28→0.26 |
| 6 | xray_refined | phenix.ligandfit | Fit ligand |
| 7 | xray_combined | phenix.pdbtools | Combine P+L |
| 8 | xray_refined | phenix.refine | Refine combined |
| 9 | xray_refined | phenix.refine | R-free 0.24→0.22 |
| 10 | xray_refined | phenix.refine | Plateau |
| 11 | xray_refined | phenix.molprobity | Validation |
| 12 | complete | STOP | Done |

**Total: 12 cycles**

---

## Key Decision Points

### 1. Experiment Type Detection

Based on file categories:
- Has `.mtz` + no maps → **X-ray**
- Has `.mrc`/`.ccp4` maps → **Cryo-EM**

### 2. Model Building Strategy

For X-ray with sequence only (no search model):
- **`phenix.predict_and_build`** is the default choice
- AlphaFold predicts structure, then places it via internal phaser

For X-ray with template PDB:
- **`phenix.phaser`** for molecular replacement

### 3. Ligand Fitting Timing

The agent waits until:
- R-free < 0.35 (model quality threshold)
- Ligand file is present
- ligandfit hasn't been run yet

This prevents fitting ligand into poor density.

### 4. Plateau Detection

From `workflows.yaml`:
```yaml
repeat:
  max_cycles: 4
  until:
    any:
      - r_free: "< target_r_free"
      - condition: plateau
        cycles: 2
        threshold: 0.005  # Less than 0.5% improvement
```

Agent stops refinement when:
- Target R-free reached, OR
- 2 consecutive cycles with < 0.005 improvement

### 5. File Selection Logic

From `command_builder.py` and `file_categorization.py`:

**Best model selection:**
1. Check `best_files["model"]` from BestFilesTracker
2. Fall back to category-based selection (prefer `refined` > `model` > `pdb`)
3. Exclude ligand files when selecting protein

**MTZ selection for refinement:**
1. Use locked R-free MTZ if set
2. Otherwise use original data MTZ (not refined outputs)

---

## Edge Cases

### What if predict_and_build fails?

The agent would:
1. See error in log
2. Try alternative: `phenix.phaser` if a template becomes available
3. Or `phenix.autosol` if anomalous signal is strong
4. Eventually STOP with error if no path forward

### What if user provides a pre-built model?

If `model.pdb` is provided alongside `data.mtz`:
- File categorization puts it in `search_model` (if looks like template) or `model` (if looks refined)
- If `model` category: skip to refinement phase
- If `search_model` category: use phaser to place it

### What if ligand file is named ambiguously?

The agent uses patterns from `programs.yaml`:
```yaml
ligand:
  extensions: [.pdb, .cif]
  prefer_patterns: [lig, ligand]
```

Files matching `lig*` or `ligand*` are categorized as ligands, not protein models.
