# Implementation Plan: PDB File Category Refactoring

## Executive Summary

The current PHENIX AI agent treats all PDB files as potential "models" with a scoring hierarchy. This conflates several semantically distinct concepts:

- **Models**: Coordinate files that explain experimental data (placed in the correct unit cell/box, with compatible symmetry)
- **Search models**: Coordinate files used as templates for molecular replacement or docking (not yet positioned)
- **Ligands**: Small molecule coordinates for fitting
- **Intermediate files**: Temporary byproducts that should never be used

This document proposes a refactoring to create distinct top-level categories with clear semantics about what each can be used for.

---

## Problem Statement

### Current Issues

1. **Semantic confusion**: A "predicted model" isn't really a model - it's a search template that hasn't been placed in the crystal/map yet.

2. **Invalid input selection**: The `best_files["model"]` tracking can select files that are invalid for refinement (e.g., predicted models that aren't in the correct unit cell).

3. **Scoring complexity**: We're using a scoring system (predicted=40, processed=50, phaser=70, refined=100) to work around a fundamental categorization problem.

4. **Brittle validation**: We have to manually exclude categories like `predicted` and `intermediate_mr` from refinement inputs.

### What Is a "Model"?

A **model** in crystallography/cryo-EM is a coordinate file that:
1. Is positioned in the correct reference frame (unit cell for X-ray, map box for cryo-EM)
2. Has compatible symmetry operators applied
3. Represents an explanation of the observed experimental data
4. Can be directly used for refinement

Things that are NOT models:
- **Predicted structures** (AlphaFold, ESMFold) - not positioned, no symmetry
- **Processed predicted structures** - trimmed but still not positioned
- **Ligand coordinates** - small molecules, not the macromolecular model
- **Search models** from the PDB - templates for MR, not positioned

---

## Proposed Architecture

### New Top-Level Categories

```yaml
# =============================================================================
# PDB FILE CATEGORIES (by semantic purpose)
# =============================================================================

# 1. MODELS - Coordinate files that explain experimental data
model:
  description: "Atomic coordinates positioned in experimental reference frame"
  extensions: [.pdb, .cif]
  semantic: "Can be directly used for refinement"
  subcategories:
    - refined          # Output from phenix.refine
    - rsr_output       # Output from real_space_refine
    - phaser_output    # Positioned by molecular replacement
    - autobuild_output # Built into density
    - docked           # Docked into cryo-EM map

# 2. SEARCH_MODEL - Templates for molecular replacement/docking
search_model:
  description: "Coordinate templates not yet positioned in experimental frame"
  extensions: [.pdb, .cif]
  semantic: "Must be placed/docked before refinement"
  subcategories:
    - predicted           # Raw AlphaFold/ESMFold output
    - processed_predicted # Trimmed for MR
    - pdb_template        # From PDB database

# 3. LIGAND - Small molecule coordinates
ligand:
  description: "Small molecule coordinates for fitting"
  extensions: [.pdb, .cif]
  semantic: "For ligandfit, not macromolecular refinement"
  subcategories:
    - ligand_pdb      # Ligand coordinates
    - ligand_cif      # Ligand with restraints

# 4. INTERMEDIATE - Temporary files (never use)
intermediate:
  description: "Intermediate processing files"
  extensions: [.pdb]
  semantic: "Never use as program input"
  subcategories:
    - intermediate_mr  # MR intermediates
    - autobuild_temp   # AutoBuild temporaries
```

### Category Relationships

```
Before (flat hierarchy with scoring):
┌─────────────────────────────────────────────────────────────┐
│ pdb (all PDB files)                                         │
│   ├── refined (100)                                         │
│   ├── rsr_output (100)                                      │
│   ├── autobuild_output (80)                                 │
│   ├── phaser_output (70)                                    │
│   ├── docked (60)                                           │
│   ├── processed_predicted (50)                              │
│   ├── predicted (40)                                        │
│   ├── ligand_pdb (excluded)                                 │
│   └── intermediate_mr (excluded)                            │
└─────────────────────────────────────────────────────────────┘

After (semantic hierarchy):
┌─────────────────────────────────────────────────────────────┐
│ model (ready for refinement)                                │
│   ├── refined                                               │
│   ├── rsr_output                                            │
│   ├── autobuild_output                                      │
│   ├── phaser_output                                         │
│   └── docked                                                │
├─────────────────────────────────────────────────────────────┤
│ search_model (needs placement first)                        │
│   ├── predicted                                             │
│   ├── processed_predicted                                   │
│   └── pdb_template                                          │
├─────────────────────────────────────────────────────────────┤
│ ligand (small molecules)                                    │
│   ├── ligand_pdb                                            │
│   └── ligand_cif                                            │
├─────────────────────────────────────────────────────────────┤
│ intermediate (never use)                                    │
│   └── intermediate_mr                                       │
└─────────────────────────────────────────────────────────────┘
```

---

## Implementation Steps

### Phase 1: Update YAML Configuration (Low Risk)

**Goal**: Define new categories without breaking existing code.

#### Step 1.1: Update `file_categories.yaml`

Add new top-level categories while keeping backward compatibility:

```yaml
# NEW: Top-level semantic categories
model:
  description: "Atomic coordinates positioned in experimental reference frame"
  extensions: [.pdb, .cif]
  semantic: "Can be directly used for refinement"
  is_parent: true  # This is a parent category, not a file matcher
  valid_for: [refinement]
  
search_model:
  description: "Coordinate templates not yet positioned"
  extensions: [.pdb, .cif]
  semantic: "Must be placed/docked before refinement"
  is_parent: true
  valid_for: [molecular_replacement, docking, prediction]
  invalid_for: [refinement]  # EXPLICIT: cannot be used for refinement
  
ligand:
  description: "Small molecule coordinates"
  extensions: [.pdb, .cif]
  semantic: "For ligand fitting"
  is_parent: true
  valid_for: [ligand_fitting]
  invalid_for: [refinement]

# Update existing subcategories to reference parents
phaser_output:
  description: "Models from molecular replacement (Phaser)"
  parent_category: model  # NEW: explicit parent
  patterns: [...]
  
predicted:
  description: "AI-predicted structures"
  parent_category: search_model  # NEW: changes from pdb to search_model
  patterns: [...]
```

#### Step 1.2: Update `metrics.yaml` Scoring

Remove scoring for non-model categories since they won't compete:

```yaml
best_files_scoring:
  model:
    stage_scores:
      refined: 100
      rsr_output: 100
      autobuild_output: 80
      phaser_output: 70
      docked: 60
      _default: 50
    # No predicted/processed_predicted - they're not models
    
  # NEW: Separate tracking for search models
  search_model:
    stage_scores:
      processed_predicted: 60  # Better for MR
      predicted: 50
      pdb_template: 40
      _default: 40
```

#### Step 1.3: Update `programs.yaml` Input Priorities

Make refinement programs explicitly require `model` category:

```yaml
phenix.refine:
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
  input_priorities:
    model:
      categories: [model]  # SIMPLIFIED: just "model"
      # No need for exclude_categories - search_model isn't model
      
phenix.phaser:
  inputs:
    required:
      search_model:  # RENAMED from model
        extensions: [.pdb]
        flag: "model="
  input_priorities:
    search_model:
      categories: [search_model, model]  # Can use both
```

**Testing checkpoint**: Run existing tests to verify backward compatibility.

---

### Phase 2: Update File Categorization Logic (Medium Risk)

**Goal**: Update the categorization code to use parent categories.

#### Step 2.1: Update `yaml_tools.py` Categorization

Modify `categorize_files()` to return both specific and parent categories:

```python
def categorize_files(files):
    """
    Categorize files with both specific and parent categories.
    
    Returns:
        dict: {
            'mtz': [...],
            'model': [...],           # Parent category
            'refined': [...],         # Specific category (also in model)
            'phaser_output': [...],   # Specific category (also in model)
            'search_model': [...],    # Parent category
            'predicted': [...],       # Specific category (also in search_model)
            ...
        }
    """
```

Key changes:
- When a file matches `predicted`, add it to BOTH `predicted` AND `search_model`
- When a file matches `refined`, add it to BOTH `refined` AND `model`
- This maintains backward compatibility while enabling parent-level queries

#### Step 2.2: Update `BestFilesTracker`

Track `model` and `search_model` as separate categories:

```python
CATEGORIES = [
    "model",         # Best positioned model for refinement
    "search_model",  # Best template for MR/docking
    "map",
    "mtz",
    "sequence",
    "ligand_cif",
]
```

Update `_classify_category()`:

```python
def _classify_category(self, path):
    """Classify into semantic categories, not just extension-based."""
    lower = path.lower()
    basename = os.path.basename(lower)
    
    if lower.endswith('.pdb') or lower.endswith('.cif'):
        # Check if it's a model (positioned) or search_model (template)
        if self._is_model(basename):
            return "model"
        elif self._is_search_model(basename):
            return "search_model"
        elif self._is_ligand(basename):
            return "ligand"
        # Default for unknown PDB files
        return "model"  # Assume model if can't determine
    
    # ... other categories

def _is_model(self, basename):
    """Check if file is a positioned model."""
    model_patterns = ['refine', 'phaser', 'placed', 'dock', 'autobuild', 
                      'overall_best', 'rsr_', '_rsr']
    return any(p in basename for p in model_patterns)

def _is_search_model(self, basename):
    """Check if file is a search template."""
    search_patterns = ['predict', 'alphafold', 'colabfold', 'esmfold', 'template']
    processed_patterns = ['processed']
    return (any(p in basename for p in search_patterns) and 
            not any(p in basename for p in ['phaser', 'refine']))
```

**Testing checkpoint**: Run `test_best_files_tracker.py` with new categories.

---

### Phase 3: Update Command Builder (Medium Risk)

**Goal**: Use semantic categories in file selection.

#### Step 3.1: Update `command_builder.py`

Simplify file selection by using parent categories:

```python
def _find_file_for_slot(self, program, input_name, input_def, 
                         available_files, context, basename_to_path):
    """
    Find appropriate file for an input slot.
    
    For 'model' inputs: use context.best_files["model"]
    For 'search_model' inputs: use context.best_files["search_model"]
    """
    # Get category priorities from program definition
    priorities = self._get_input_priorities(program, input_name)
    
    # Check best_files first
    if "model" in priorities:
        best_model = context.best_files.get("model")
        if best_model and os.path.exists(best_model):
            return best_model
    
    if "search_model" in priorities:
        best_search = context.best_files.get("search_model")
        if best_search and os.path.exists(best_search):
            return best_search
    
    # Fall back to categorized files
    ...
```

#### Step 3.2: Update Validation

Add explicit validation that refinement inputs are models:

```python
def _validate_refinement_inputs(self, program, selected_files, context):
    """Ensure refinement programs get actual models, not search templates."""
    if program not in ["phenix.refine", "phenix.real_space_refine"]:
        return []
    
    issues = []
    model_file = selected_files.get("model")
    if model_file:
        category = self._get_file_category(model_file)
        if category in ["predicted", "processed_predicted", "search_model"]:
            issues.append(SanityIssue(
                level="critical",
                message=f"Cannot refine {category} - must place in unit cell first",
                suggestion="Run phaser or dock_in_map before refinement"
            ))
    return issues
```

**Testing checkpoint**: Run `test_command_builder.py` with validation tests.

---

### Phase 4: Update Workflow Logic (Low-Medium Risk)

**Goal**: Update workflow state detection to use semantic categories.

#### Step 4.1: Update `workflow_state.py`

Use semantic categories for state detection:

```python
def detect_workflow_state(categorized_files, history, metrics, experiment_type):
    """
    Detect workflow state using semantic categories.
    """
    has_model = bool(categorized_files.get("model"))
    has_search_model = bool(categorized_files.get("search_model"))
    
    # X-ray workflow states
    if experiment_type == "xray":
        if has_model and _has_rfree_stats(metrics):
            return "xray_refined"
        elif has_model:
            return "xray_has_model"  # Ready for refinement
        elif has_search_model:
            return "xray_has_search_model"  # Ready for MR
        ...
```

#### Step 4.2: Update Valid Programs Logic

Programs are valid based on what categories are available:

```python
def get_valid_programs(categorized_files, workflow_state):
    """
    Get valid programs based on available categories.
    """
    has_model = bool(categorized_files.get("model"))
    has_search_model = bool(categorized_files.get("search_model"))
    has_mtz = bool(categorized_files.get("mtz"))
    
    valid = []
    
    # Refinement needs a model
    if has_model and has_mtz:
        valid.append("phenix.refine")
    
    # MR needs a search model
    if has_search_model and has_mtz:
        valid.append("phenix.phaser")
    
    # Processing needs a predicted model
    if categorized_files.get("predicted"):
        valid.append("phenix.process_predicted_model")
    
    ...
```

**Testing checkpoint**: Run `test_workflow_state.py` with new logic.

---

### Phase 5: Update Sanity Checker (Low Risk)

**Goal**: Add explicit checks for category misuse.

#### Step 5.1: Add Category Validation Checks

```python
class SanityChecker:
    def _check_refinement_inputs(self, context, session_info):
        """Ensure refinement uses models, not search templates."""
        issues = []
        
        planned_program = context.get("planned_program")
        if planned_program not in ["phenix.refine", "phenix.real_space_refine"]:
            return issues
        
        # Check if best_model is actually a model
        best_model = session_info.get("best_files", {}).get("model")
        if best_model:
            basename = os.path.basename(best_model).lower()
            if 'predict' in basename and 'phaser' not in basename:
                issues.append(SanityIssue(
                    level="critical",
                    message="Best model is a prediction, not a positioned model",
                    suggestion="Run phaser first to position the model"
                ))
        
        return issues
```

**Testing checkpoint**: Run `test_sanity_checker.py` with new checks.

---

### Phase 6: Migration and Cleanup (Low Risk)

**Goal**: Clean up deprecated patterns and update tests.

#### Step 6.1: Deprecation Warnings

Add warnings for old patterns:

```python
# In yaml_tools.py
def categorize_files(files):
    result = _categorize_files_new(files)
    
    # Backward compatibility: pdb includes all PDB-like categories
    # but warn if code uses it
    if 'pdb' in result:
        warnings.warn(
            "Category 'pdb' is deprecated. Use 'model', 'search_model', or 'ligand'",
            DeprecationWarning
        )
    
    return result
```

#### Step 6.2: Update Tests

Add tests for new categories:

```python
class TestSemanticCategories(unittest.TestCase):
    def test_predicted_is_search_model(self):
        files = ["predicted_model.pdb"]
        cats = categorize_files(files)
        self.assertIn(files[0], cats.get("search_model", []))
        self.assertNotIn(files[0], cats.get("model", []))
    
    def test_phaser_output_is_model(self):
        files = ["PHASER.1.pdb"]
        cats = categorize_files(files)
        self.assertIn(files[0], cats.get("model", []))
        self.assertNotIn(files[0], cats.get("search_model", []))
    
    def test_refinement_requires_model(self):
        """Refinement should reject search_model inputs."""
        ...
```

#### Step 6.3: Update Documentation

Update README and ARCHITECTURE docs to explain new categories.

---

## Backward Compatibility

### Preserved Behaviors

1. **Subcategory names unchanged**: `refined`, `predicted`, etc. still exist
2. **Programs still work**: Just with cleaner input selection
3. **Scoring still works**: But only within semantic categories
4. **Existing tests pass**: With minor updates for new categories

### Breaking Changes (Intentional)

1. **`pdb` category deprecated**: Use `model`, `search_model`, or `ligand`
2. **Predicted models can't be "best model"**: They're "best search_model"
3. **Stricter validation**: Won't allow refinement of search templates

---

## Risk Assessment

| Phase | Risk | Mitigation |
|-------|------|------------|
| 1. YAML Config | Low | Backward-compatible additions |
| 2. Categorization | Medium | Maintain dual categories |
| 3. Command Builder | Medium | Add validation, not removal |
| 4. Workflow Logic | Low-Medium | Feature flags for rollback |
| 5. Sanity Checker | Low | Additive checks only |
| 6. Migration | Low | Gradual deprecation |

---

## Testing Strategy

### Unit Tests
- `test_file_categorization.py`: Add semantic category tests
- `test_best_files_tracker.py`: Test separate model/search_model tracking
- `test_command_builder.py`: Test category-based file selection
- `test_workflow_state.py`: Test state detection with new categories

### Integration Tests
- Run full X-ray workflow: xtriage → predict → process → phaser → refine
- Verify PHASER output is selected for refinement (not predicted)
- Verify predicted model is selected for MR (not refined)

### Regression Tests
- Run all existing tests
- Compare outputs on standard test cases

---

## Success Criteria

1. **Clear semantics**: Each PDB file is in exactly one top-level category
2. **Correct selection**: Refinement always gets models, MR always gets templates
3. **No scoring hacks**: Don't need score manipulation to get correct behavior
4. **Backward compatible**: Existing workflows still work
5. **Well-documented**: Clear guidance for adding new categories

---

## Timeline Estimate

| Phase | Estimated Effort |
|-------|------------------|
| Phase 1: YAML Config | 2-3 hours |
| Phase 2: Categorization | 4-6 hours |
| Phase 3: Command Builder | 3-4 hours |
| Phase 4: Workflow Logic | 2-3 hours |
| Phase 5: Sanity Checker | 1-2 hours |
| Phase 6: Migration | 2-3 hours |
| **Total** | **14-21 hours** |

---

## Appendix: Current vs Proposed Category Mapping

| Current Category | Current Parent | Proposed Parent | Rationale |
|-----------------|----------------|-----------------|-----------|
| `refined` | `pdb` | `model` | Refined = positioned & explains data |
| `rsr_output` | `pdb` | `model` | RSR output = positioned in map |
| `phaser_output` | `pdb` | `model` | PHASER output = positioned in cell |
| `autobuild_output` | `pdb` | `model` | Built into density = positioned |
| `docked` | `pdb` | `model` | Docked = positioned in map |
| `predicted` | `pdb` | `search_model` | Not positioned |
| `processed_predicted` | `pdb` | `search_model` | Still not positioned |
| `ligand_pdb` | `pdb` | `ligand` | Small molecule |
| `intermediate_mr` | `pdb` | `intermediate` | Never use |

---

## Design Decisions (Resolved)

### 1. Should `pdb_template` be a new subcategory of `search_model`?

**Decision: Yes.**

**Rationale:** A `pdb_template` (e.g., a homolog downloaded from the PDB) is functionally distinct from a prediction. It might need chain renaming or pruning using `sculptor`, which is different from how AlphaFold models are processed with `process_predicted_model`.

**Implementation:**
```yaml
pdb_template:
  description: "Homologous structure from PDB for molecular replacement"
  parent_category: search_model
  patterns:
    - "*template*"
    - "*homolog*"
    - "*pdb[0-9][a-z0-9][a-z0-9][a-z0-9]*"  # PDB ID patterns
  notes: "Consider running Sculptor for chain editing before MR"
```

This allows specific advice like: "If using a raw PDB template, consider running Sculptor first."

---

### 2. Should there be a `model_with_ligand` subcategory distinct from `refined`?

**Decision: Yes, but as a "tag" or secondary trait, not a replacement.**

**Rationale:** A `model_with_ligand` is still a `refined` model (usually). You don't want to lose the "refined" status just because you added a ligand. The ligand presence is additional information, not a replacement classification.

**Implementation:** Allow files to belong to multiple subcategories using `also_in`:
```yaml
with_ligand:
  description: "Models that include fitted ligands"
  parent_category: model
  patterns:
    - "*with_ligand*"
    - "*_liganded*"
    - "*plus_lig*"
  also_in: [refined]  # Don't lose refined status
  notes: "Combined protein + ligand models"
```

**Example:** A file named `refine_005_with_ligand.pdb` should match BOTH:
- `refined` (so it scores 100 points as a model)
- `with_ligand` (so the agent knows the ligand is present)

The categorization code should add the file to both categories.

---

### 3. How should we handle mmCIF files (`.cif`) that could be models or restraints?

**Decision: The `ligand` vs `model` top-level split handles this, with content detection as fallback.**

**Logic:**
- If the `.cif` contains restraint data (e.g., `data_comp_`, `_chem_comp`), it's a `ligand_cif`
- If the `.cif` contains atomic coordinates (e.g., `_atom_site.`), it's a `model_cif`

**Implementation:** Use naming conventions first, content detection second:
```python
def _classify_cif_file(self, path):
    """Classify CIF file as ligand restraints or model coordinates."""
    basename = os.path.basename(path).lower()
    
    # Naming convention checks first (fast)
    if any(p in basename for p in ['restraint', 'constraint', 'lig_', 'ligand_']):
        return "ligand_cif"
    if any(p in basename for p in ['refine', 'model', 'coord']):
        return "model_cif"
    
    # Content detection fallback (slower but accurate)
    try:
        with open(path, 'r') as f:
            content = f.read(4096)  # Read first 4KB
            if '_atom_site.' in content or '_atom_site_' in content:
                return "model_cif"
            if '_chem_comp.' in content or 'data_comp_' in content:
                return "ligand_cif"
    except:
        pass
    
    # Default based on file size (restraints are usually small)
    if os.path.getsize(path) < 50000:  # < 50KB
        return "ligand_cif"
    return "model_cif"
```

---

### 4. Should the agent "promote" a `search_model` to `model`?

**Decision: No. The PROCESS promotes it, not the agent.**

**Rationale:** You don't just "re-label" a file. A `search_model` becomes a `model` only by passing through a program like Phaser or dock_in_map. The file type transformation IS the proof that the scientific work was done.

**The Flow:**
```
1. Input:   prediction.pdb        (Type: search_model/predicted)
2. Process: phenix.phaser
3. Output:  PHASER.1.pdb          (Type: model/phaser_output)
```

**Implementation:** No special "promotion" logic needed. The output file from Phaser has a different name (`PHASER.*.pdb`), which naturally classifies it as `model/phaser_output`. The transformation is implicit in the program execution.

**Key Insight:** This is why the category is determined by filename patterns - the program's output naming convention encodes the semantic transformation that occurred.

---

## Concrete Implementation Details

This section provides the exact file changes needed for each phase.

### Phase 1A: Update `file_categories.yaml`

Replace the current flat hierarchy with semantic parent categories:

```yaml
# =============================================================================
# TOP-LEVEL SEMANTIC CATEGORIES
# =============================================================================
# These are parent categories that define WHAT a file can be used for.
# Files are classified into exactly ONE of these based on their semantic purpose.

model:
  description: "Atomic coordinates positioned in experimental reference frame"
  extensions: [.pdb, .cif]
  is_semantic_parent: true
  valid_for:
    - refinement
    - validation
    - ligand_fitting
    - analysis
  notes: "Can be directly used for refinement - positioned in unit cell (X-ray) or map box (cryo-EM)"

search_model:
  description: "Coordinate templates not yet positioned in experimental frame"
  extensions: [.pdb, .cif]
  is_semantic_parent: true
  valid_for:
    - molecular_replacement
    - docking
    - model_processing
  invalid_for:
    - refinement  # CRITICAL: Cannot refine an unpositioned model
  notes: "Must be placed via Phaser (X-ray) or dock_in_map (cryo-EM) before refinement"

ligand:
  description: "Small molecule coordinates and restraints"
  extensions: [.pdb, .cif]
  is_semantic_parent: true
  valid_for:
    - ligand_fitting
    - restraint_generation
  invalid_for:
    - refinement  # Not a macromolecular model
  notes: "For ligandfit, not macromolecular refinement"

intermediate:
  description: "Intermediate processing files - never use as input"
  extensions: [.pdb]
  is_semantic_parent: true
  valid_for: []  # Nothing!
  invalid_for:
    - all
  notes: "Temporary byproducts that should never be used as program inputs"

# =============================================================================
# MODEL SUBCATEGORIES (children of 'model')
# =============================================================================

refined:
  description: "X-ray refined models"
  parent_category: model
  patterns:
    - "*refine*"
    - "refine_[0-9][0-9][0-9]_*.pdb"
  excludes:
    - "*real_space*"
    - "*rsr*"
  notes: "Output from phenix.refine"

rsr_output:
  description: "Real-space refined models (cryo-EM)"
  parent_category: model
  patterns:
    - "*real_space_refined*"
    - "*rsr_*"
    - "*_rsr*"
  notes: "Output from phenix.real_space_refine"

phaser_output:
  description: "Models from molecular replacement (Phaser)"
  parent_category: model
  patterns:
    - "phaser*"
    - "PHASER*"
    - "*phaser*"
  notes: "Positioned in unit cell via MR, ready for refinement"

autobuild_output:
  description: "Models from automated building"
  parent_category: model
  patterns:
    - "*autobuild*"
    - "*overall_best*"
    - "*built*"
  excludes:
    - "*predict*"
  notes: "Output from phenix.autobuild - built into density"

docked:
  description: "Models docked into cryo-EM maps"
  parent_category: model
  patterns:
    - "*dock*map*"
    - "*docked*"
    - "placed_model*"
  notes: "Output from phenix.dock_in_map"

with_ligand:
  description: "Models that include fitted ligands"
  parent_category: model
  patterns:
    - "*with_ligand*"
    - "*_liganded*"
  also_in: [refined]  # Preserve refined status
  notes: "Combined protein + ligand models"

# =============================================================================
# SEARCH_MODEL SUBCATEGORIES (children of 'search_model')
# =============================================================================

predicted:
  description: "AI-predicted models (AlphaFold, ESMFold, etc.)"
  parent_category: search_model
  patterns:
    - "*predict*"
    - "*alphafold*"
    - "*colabfold*"
    - "*AF-*"
    - "*ESMFold*"
  excludes:
    - "*phaser*"  # PHASER output contains 'predict' in path but IS a model
    - "*processed*"
  notes: "Raw predictions - need processing before MR"

processed_predicted:
  description: "Processed predicted models (trimmed, B-factors adjusted)"
  parent_category: search_model
  patterns:
    - "*processed*predict*"
    - "*processed*model*"
  notes: "Ready for MR but NOT for refinement (not yet positioned)"

pdb_template:
  description: "Homologous structure from PDB for molecular replacement"
  parent_category: search_model
  patterns:
    - "*template*"
    - "*homolog*"
    - "*sculptor*"
  notes: "Consider running Sculptor for chain editing before MR"

# =============================================================================
# LIGAND SUBCATEGORIES (children of 'ligand')
# =============================================================================

ligand_pdb:
  description: "Standalone ligand coordinate files"
  parent_category: ligand
  patterns:
    - "lig*.pdb"
    - "ligand*.pdb"
    - "LIG*.pdb"
  excludes:
    - "*ligand_fit*"
    - "*with_ligand*"
  max_basename_length: 20
  notes: "Small molecule coordinates for fitting"

ligand_cif:
  description: "Ligand restraint/geometry files"
  parent_category: ligand
  extensions: [.cif]
  patterns:
    - "*restraint*"
    - "*constraint*"
    - "lig_*.cif"
  excludes:
    - "*refine*"
  notes: "Chemical restraints for ligand refinement"

# =============================================================================
# INTERMEDIATE SUBCATEGORIES (children of 'intermediate')
# =============================================================================

intermediate_mr:
  description: "Intermediate MR files from dock_in_map"
  parent_category: intermediate
  patterns:
    - "run_mr*"
    - "*mr.[0-9]*"
    - "*_mr_*temp*"
  notes: "NEVER use - temporary files from internal MR"

autobuild_temp:
  description: "Temporary AutoBuild files"
  parent_category: intermediate
  patterns:
    - "*AutoBuild_run_*"
    - "*_temp_*"
  notes: "NEVER use - intermediate AutoBuild files"

carryover_temp:
  description: "CarryOn directory intermediate files"
  parent_category: intermediate
  patterns:
    - "*CarryOn/*"
    - "*_CarryOn/*"
  notes: "NEVER use - predict_and_build intermediate files"
```

### Phase 1B: Update `metrics.yaml` Scoring

Separate scoring for `model` and `search_model` categories:

```yaml
best_files_scoring:
  # -------------------------------------------------------------------------
  # MODEL SCORING (for refinement input selection)
  # -------------------------------------------------------------------------
  # Only files in the 'model' parent category compete here.
  # Score determines which positioned model is used for refinement.
  
  model:
    stage_scores:
      refined: 100           # Best - already refined
      rsr_output: 100        # Best for cryo-EM
      autobuild_output: 80   # Built into density
      phaser_output: 70      # Positioned via MR
      docked: 60             # Positioned via docking
      with_ligand: 100       # Same as refined (it's a tag, not replacement)
      _default: 50           # Unknown model type
    
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

  # -------------------------------------------------------------------------
  # SEARCH_MODEL SCORING (for MR/docking input selection)
  # -------------------------------------------------------------------------
  # Only files in the 'search_model' parent category compete here.
  # Score determines which template is used for molecular replacement.
  
  search_model:
    stage_scores:
      processed_predicted: 70  # Best - trimmed and ready for MR
      pdb_template: 60         # Good - may need sculptor
      predicted: 50            # Raw prediction - needs processing
      _default: 40
    
    metric_scores:
      # Could add pLDDT scoring for AlphaFold models
      plddt_mean:
        max_points: 30
        formula: linear
        best_value: 90
        worst_value: 50

  # -------------------------------------------------------------------------
  # LIGAND SCORING
  # -------------------------------------------------------------------------
  ligand:
    stage_scores:
      ligand_cif: 60   # Prefer restraints over bare coordinates
      ligand_pdb: 50
      _default: 40
```

### Phase 2: Update `BestFilesTracker`

Key changes to `agent/best_files_tracker.py`:

```python
class BestFilesTracker:
    """
    Tracks the best file of each SEMANTIC category throughout a session.
    
    Categories tracked:
        - model: Best positioned model for refinement
        - search_model: Best template for MR/docking
        - map: Best full cryo-EM map
        - mtz: Best reflection data (with R-free flags)
        - map_coefficients: Best map coefficients
        - sequence: Sequence file
        - ligand: Best ligand file (coordinates or restraints)
    """
    
    # Semantic categories we track
    CATEGORIES = [
        "model",         # Positioned models for refinement
        "search_model",  # Templates for MR/docking
        "map",
        "mtz",
        "map_coefficients",
        "sequence",
        "ligand",
    ]
    
    # Mapping from subcategory to parent category
    PARENT_CATEGORIES = {
        # model subcategories
        "refined": "model",
        "rsr_output": "model",
        "phaser_output": "model",
        "autobuild_output": "model",
        "docked": "model",
        "with_ligand": "model",
        # search_model subcategories
        "predicted": "search_model",
        "processed_predicted": "search_model",
        "pdb_template": "search_model",
        # ligand subcategories
        "ligand_pdb": "ligand",
        "ligand_cif": "ligand",
        # intermediate - NOT tracked
        "intermediate_mr": None,
        "autobuild_temp": None,
        "carryover_temp": None,
    }
    
    def _classify_category(self, path):
        """
        Classify file into semantic parent category.
        
        Returns the PARENT category (model, search_model, ligand, etc.),
        not the subcategory.
        """
        if not path:
            return None
        
        lower = path.lower()
        basename = os.path.basename(lower)
        
        # Check for intermediate files FIRST - never track these
        if self._is_intermediate_file(path):
            return None
        
        # PDB/CIF files need semantic classification
        if lower.endswith('.pdb') or (lower.endswith('.cif') and self._is_model_cif(path)):
            # Determine subcategory first
            subcategory = self._classify_pdb_subcategory(basename)
            # Map to parent category
            return self.PARENT_CATEGORIES.get(subcategory, "model")
        
        # CIF files that are restraints
        if lower.endswith('.cif'):
            return "ligand"
        
        # Other file types (unchanged)
        if lower.endswith('.mtz'):
            return "mtz"
        if lower.endswith(('.mrc', '.ccp4', '.map')):
            return "map"
        if lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            return "sequence"
        
        return None
    
    def _classify_pdb_subcategory(self, basename):
        """
        Classify PDB file into specific subcategory.
        
        Order matters! More specific patterns should come first.
        """
        # MODEL subcategories (positioned, ready for refinement)
        if 'refine' in basename and 'real_space' not in basename:
            return "refined"
        if 'real_space_refined' in basename or 'rsr_' in basename:
            return "rsr_output"
        if 'phaser' in basename:
            return "phaser_output"
        if 'overall_best' in basename or 'autobuild' in basename:
            return "autobuild_output"
        if 'placed' in basename or 'dock' in basename:
            return "docked"
        if 'with_ligand' in basename or '_liganded' in basename:
            return "with_ligand"
        
        # SEARCH_MODEL subcategories (templates, not positioned)
        if 'processed' in basename:
            return "processed_predicted"
        if 'predict' in basename or 'alphafold' in basename:
            return "predicted"
        if 'template' in basename or 'homolog' in basename:
            return "pdb_template"
        
        # LIGAND subcategories
        if self._is_ligand_pdb(basename):
            return "ligand_pdb"
        
        # Default: assume it's a model (conservative)
        return "refined"  # Default to model category
    
    def _is_intermediate_file(self, path):
        """Check if file is intermediate/temporary."""
        intermediate_patterns = [
            '/run_mr/',
            'run_mr.',
            '_mr.',
            '/AutoBuild_run_',
            '/CarryOn/',
            '_CarryOn/',
            '/temp/',
            '/tmp/',
            '.tmp.',
            'mask.ccp4',
        ]
        path_lower = path.lower()
        return any(p.lower() in path_lower for p in intermediate_patterns)
    
    def _is_ligand_pdb(self, basename):
        """Check if PDB file is a ligand (small molecule)."""
        ligand_patterns = ['lig.pdb', 'ligand', 'lig_', 'lig-']
        is_ligand = any(p in basename for p in ligand_patterns)
        is_not_fit = 'ligand_fit' not in basename and 'ligandfit' not in basename
        is_small_name = len(basename) < 20
        return is_ligand and is_not_fit and is_small_name
```

### Phase 3: Update `command_builder.py`

Simplify file selection using semantic categories:

```python
def _find_file_for_slot(self, program, input_name, input_def,
                        available_files, context, basename_to_path):
    """
    Find appropriate file for an input slot using semantic categories.
    """
    # Get the required category for this input
    required_category = self._get_required_category(program, input_name)
    
    # For model inputs: only accept files from 'model' category
    if required_category == "model":
        best_model = context.best_files.get("model")
        if best_model and os.path.exists(best_model):
            self._log(context, f"BUILD: Using best model: {os.path.basename(best_model)}")
            return best_model
        
        # Fall back to categorized files in model category
        for cat in ["refined", "rsr_output", "phaser_output", "autobuild_output", "docked"]:
            if cat in context.categorized_files:
                files = context.categorized_files[cat]
                if files:
                    self._log(context, f"BUILD: Using {cat} file: {os.path.basename(files[0])}")
                    return files[0]
        
        # CRITICAL: Do NOT fall back to search_model!
        self._log(context, "BUILD: WARNING - No model found, cannot proceed with refinement")
        return None
    
    # For search_model inputs (MR, docking): accept search_model OR model
    if required_category == "search_model":
        # Prefer search_model (that's what MR wants)
        best_search = context.best_files.get("search_model")
        if best_search and os.path.exists(best_search):
            self._log(context, f"BUILD: Using best search_model: {os.path.basename(best_search)}")
            return best_search
        
        # Can also use a model as template (less common but valid)
        best_model = context.best_files.get("model")
        if best_model and os.path.exists(best_model):
            self._log(context, f"BUILD: Using model as MR template: {os.path.basename(best_model)}")
            return best_model
        
        return None
    
    # ... handle other categories ...

def _get_required_category(self, program, input_name):
    """
    Get the semantic category required for a program input.
    """
    # Refinement programs require 'model'
    if program in ["phenix.refine", "phenix.real_space_refine"]:
        if input_name in ["model", "pdb", "pdb_file"]:
            return "model"
    
    # MR/docking programs want 'search_model'
    if program in ["phenix.phaser", "phenix.dock_in_map"]:
        if input_name in ["model", "pdb", "search_model"]:
            return "search_model"
    
    # Processing programs want raw predictions
    if program == "phenix.process_predicted_model":
        if input_name == "model":
            return "search_model"  # Specifically predicted subcategory
    
    return None  # No specific requirement
```

### Phase 4: Add Validation in `sanity_checker.py`

```python
def _check_refinement_model_category(self, context, session_info):
    """
    Ensure refinement programs get actual models, not search templates.
    """
    issues = []
    
    planned_program = context.get("planned_program")
    if planned_program not in ["phenix.refine", "phenix.real_space_refine"]:
        return issues
    
    # Check what's in best_files
    best_files = session_info.get("best_files", {})
    
    # If there's no model but there IS a search_model, that's wrong
    if not best_files.get("model") and best_files.get("search_model"):
        issues.append(SanityIssue(
            level="critical",
            code="REFINE_NEEDS_MODEL",
            message="Cannot refine: only search templates available, no positioned model",
            suggestion="Run phenix.phaser (X-ray) or phenix.dock_in_map (cryo-EM) first to position the model",
            context={
                "has_search_model": True,
                "search_model": os.path.basename(best_files.get("search_model", "")),
            }
        ))
    
    return issues
```
