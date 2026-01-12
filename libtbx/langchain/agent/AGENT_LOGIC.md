# PHENIX AI Agent Logic Documentation

This document describes all decision logic used by the AI agent.
Generated automatically from source code and configuration files.

## Table of Contents
1. [Tiered Decision Architecture](#tiered-decision-architecture)
2. [Workflow States](#workflow-states)
3. [Program Templates](#program-templates)
4. [Stop Conditions](#stop-conditions)
5. [File Categorization](#file-categorization)
6. [Key Thresholds](#key-thresholds)

---
## Tiered Decision Architecture

The agent uses a three-tier decision system:

| Tier | Type | Description | Override? |
|------|------|-------------|-----------|
| 1 | Hard Constraints | Must be followed | No |
| 2 | Strong Defaults | Applied automatically | Yes (with warning) |
| 3 | Soft Guidance | Suggestions only | Yes (no warning) |

### Tier 1: Hard Constraints

Cannot be overridden under any circumstances

| Rule | Description | Applies To |
|------|-------------|------------|
| xtriage_before_mr | Must run xtriage before molecular replacement for X-ray | xray |
| mtriage_before_refine | Must run mtriage before refinement for cryo-EM | cryoem |
| phaser_requires_model | Phaser requires a search model (.pdb file) | xray |
| autobuild_requires_phases | Autobuild requires phased MTZ (from phaser or refine output) | xray |
| autobuild_requires_sequence | Autobuild requires sequence file | xray |
| refine_requires_model_and_data | Refine requires both model (.pdb/.cif) and data (.mtz) | xray |
| rsr_requires_model_and_map | Real-space refine requires both model and map | cryoem |
| no_duplicate_commands | Cannot run identical command twice in succession | all |
| ligandfit_requires_ligand | Ligandfit requires ligand coordinates file | all |

### Tier 2: Strong Defaults

Applied automatically; LLM can override but warning is logged

| Default | Condition | Value | Override Warning |
|---------|-----------|-------|------------------|
| generate_rfree_flags | refine_count == 0 | True | Overriding R-free flag generat... |
| use_existing_rfree_flags | refine_count > 0 | False | Generating new R-free flags af... |
| water_building | {'refine_count': {'operator': '>=', 'val... | True | Skipping water building when c... |
| twinning | {'has_twinning': {'operator': '==', 'val... |  | Ignoring detected twinning |
| autobuild_resolution_limit | resolution < 2.0 | 2.0 | Running autobuild at full reso... |
| autobuild_include_model | has_refined_model | auto_select_best_refined | Running autobuild without star... |

### Tier 3: Soft Guidance

Suggestions only; LLM decides freely, no warnings logged

- **When** `r_free < good_model_threshold`: Consider stopping - R-free {r_free:.3f} is below target {good_model_threshold:.3f}
- **When** `resolution < 1.5`: Consider anisotropic ADP refinement for high resolution ({resolution:.1f}Å) data
- **When** `resolution < 1.2`: Consider adding riding hydrogens for very high resolution ({resolution:.1f}Å) data
- **When** `r_free < good_model_threshold and not validation_done`: Consider running validation to confirm model quality
- **When** `refine_count >= 3 and not validation_done`: Consider running validation after {refine_count} refinement cycles
- **When** `has_useful_anomalous and not has_model`: Strong anomalous signal detected - consider experimental phasing with autosol

### Program Rankings by State

The graph recommends programs in priority order:

**xray_initial** - Before any analysis

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.xtriage | always |

**xray_analyzed** - After xtriage, need to get model

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.autosol | has_useful_anomalous and has_sequence |
| 2 | phenix.predict_and_build | has_sequence |
| 3 | phenix.autosol | has_sequence and not has_useful_anomalous |
| 4 | phenix.phaser | has_search_model |

**xray_has_prediction** - Have AlphaFold prediction, need to process

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.process_predicted_model | always |

**xray_model_processed** - Prediction processed, need MR

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.phaser | always |

**xray_has_phases** - After experimental phasing, need to build

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.autobuild | always |

**xray_has_model** - Have placed model, need refinement

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.refine | always |

**xray_refined** - After refinement, multiple options

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.ligandfit | has_ligand and not has_ligand_fit and model_is_goo... |
| 2 | phenix.autobuild | r_free_unknown_or_high and has_sequence and has_ph... |
| 3 | phenix.molprobity | suggest_validation |
| 4 | phenix.refine | r_free_unknown_or_above_target |
| 5 | phenix.ligandfit | has_ligand and not has_ligand_fit and not model_is... |
| 6 | STOP | always |

**xray_has_ligand** - Ligand fitted, need to combine

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.pdbtools | always |

**cryoem_initial** - Before any analysis

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.mtriage | always |

**cryoem_analyzed** - After mtriage, need model

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.predict_and_build | has_sequence |
| 2 | phenix.dock_in_map | has_model |

**cryoem_has_model** - Have model in map, need refinement

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.real_space_refine | always |

**cryoem_refined** - After refinement, multiple options

| Priority | Program | Condition |
|----------|---------|-----------|
| 1 | phenix.ligandfit | has_ligand and not has_ligand_fit |
| 2 | phenix.real_space_refine | map_cc < success_cc |
| 3 | phenix.validation_cryoem | suggest_validation |
| 4 | STOP | map_cc >= success_cc |

---
## Workflow States

### X-ray Crystallography Workflow

```
xray_initial → xtriage → xray_analyzed
                              ↓
              [if sequence] predict_and_build → xray_has_prediction
              [if model]    phaser ─────────────────┐
                                                    ↓
                     xray_has_prediction → process_predicted_model
                                                    ↓
                                    xray_model_processed → phaser
                                                              ↓
                                                    xray_has_model → refine
                                                              ↓
                                                       xray_refined
                                                     ↓    ↓    ↓
                                              refine  autobuild  STOP
```

| State | Valid Programs | Description |
|-------|----------------|-------------|
| xray_initial | phenix.xtriage | Need to analyze data quality with xtriage first |
| xray_has_prediction | phenix.process_predicted_model | Have AlphaFold prediction |
| xray_model_processed | phenix.phaser | Predicted model processed |
| xray_analyzed | STOP | STUCK: No sequence file for AlphaFold and no searc |
| xray_analyzed | phenix.autobuild | Experimental phasing complete |
| xray_has_model | phenix.refine | Have placed model |
| xray_has_ligand | phenix.pdbtools | Ligand fitted |

### Cryo-EM Workflow

**Automated path** (maximum_automation=True):
```
cryoem_initial → mtriage → cryoem_analyzed → predict_and_build(full)
                                                      ↓
                                            cryoem_has_model → real_space_refine
                                                      ↓
                                               cryoem_refined → STOP
```

**Stepwise path** (maximum_automation=False):
```
cryoem_initial → mtriage → cryoem_analyzed → predict(stop_after_predict)
                                                      ↓
                                          process_predicted_model → dock_in_map
                                                      ↓
                                            cryoem_has_model → real_space_refine
```

| State | Valid Programs | Description |
|-------|----------------|-------------|
| cryoem_initial | phenix.mtriage | Need to analyze map quality with mtriage first |
| cryoem_analyzed | STOP | STUCK: No sequence for AlphaFold and no model to d |
| cryoem_has_model | phenix.real_space_refine | Have model in map |

---
## Program Templates

### phenix.xtriage

**Description:** Data quality analysis for X-ray crystallography.

**Required Files:**
- `data`: .mtz, .sca, .hkl, .sdf (required)
- `sequence`: .fa, .fasta, .seq, .dat (optional)

### phenix.predict_and_build

**Description:** AlphaFold prediction and model building - use when you have sequence but no model!

**Required Files:**
- `sequence`: .fa, .fasta, .seq, .dat (required)
- `data`: .mtz, .sca, .hkl, .sdf (optional)

**Default Flags:**
- `control.nproc=4`

**Strategy Options:**
- `stop_after_predict`: true → `predict_and_build.stop_after_predict=True`
- `quick`: true → `quick=True`
- `nproc`: control.nproc={}
- `resolution`: crystal_info.resolution={}

**Usage Hints:**
- If stop_after_predict is False (default), you MUST provide resolution (e.g., 2.5)
- If you only want the AlphaFold model without building, set stop_after_predict: true
- Resolution can be found in xtriage output or MTZ file header

### phenix.process_predicted_model

**Description:** Process AlphaFold model for molecular replacement.

**Required Files:**
- `model`: .pdb (required)

### phenix.phaser

**Description:** Molecular replacement using a search model.

**Required Files:**
- `data`: .mtz, .sca, .hkl, .sdf (required)
- `model`: .pdb (required)
- `sequence`: .fa, .fasta, .seq, .dat (optional)

**Default Flags:**
- `phaser.mode=MR_AUTO`

**Strategy Options:**
- `component_copies`: phaser.ensemble.copies={}

### phenix.refine

**Description:** Automated refinement of atomic models against X-ray data.

**Required Files:**
- `model`: .pdb, .cif (required)
- `data`: .mtz, .sca, .hkl, .sdf (required)

**Strategy Options:**
- `use_simulated_annealing`: true → `main.simulated_annealing=True`
- `generate_rfree_flags`: true → `xray_data.r_free_flags.generate=True`
- `optimize_xyz`: true → `refine.strategy=individual_sites`
- `optimize_adp`: true → `refine.strategy=individual_adp`
- `add_waters`: true → `ordered_solvent=True`
- `twin_law`: refinement.twinning.twin_law={}
- `anisotropic_adp`: true → `refine.adp.individual.anisotropic="not (water or element H)"`
- `resolution_limit`: xray_data.high_resolution={}
- `nproc`: refinement.main.nproc={}

**Usage Hints:**
- FIRST refinement after phaser: use generate_rfree_flags=true
- SUBSEQUENT refinements: do NOT generate new R-free flags (use existing from MTZ)
- Add waters (ordered_solvent=True) only after 2+ refine cycles AND R-free < 0.35 AND resolution <= 3.0Å
- Do NOT use simulated_annealing for initial refinement
- Only use simulated_annealing=True if R-free > 0.40 AND previous refinement didn't improve
- If TWINNING detected, include twin_law (e.g., twin_law='-h,-k,l')
- At very high resolution (<1.5Å), consider anisotropic_adp=true

### phenix.ligandfit

**Description:** Fit a ligand into difference density.

**Required Files:**
- `model`: .pdb (required)
- `data`: .mtz, .sca, .hkl, .sdf (required)
- `ligand`: .pdb, .cif (required)

**Default Flags:**
- `file_info.input_labels="2FOFCWT PH2FOFCWT"`
- `general.nproc=4`

**Strategy Options:**
- `nproc`: general.nproc={}

### phenix.pdbtools

**Description:** Combine PDB files (e.g., protein + ligand).

**Required Files:**
- `protein`: .pdb (required)
- `ligand`: .pdb (required)

**Strategy Options:**
- `output_name`: output.file_name={}

### phenix.real_space_refine

**Description:** Refinement of atomic models against Cryo-EM maps.

**Required Files:**
- `model`: .pdb, .cif (required)
- `map`: .mrc, .ccp4, .map (required)

**Default Flags:**
- `run=minimization`

**Strategy Options:**
- `nproc`: nproc={}

### phenix.autobuild

**Description:** Iterative model building and refinement.

**Required Files:**
- `data`: .mtz, .sca, .hkl, .sdf (required)
- `model`: .pdb (optional)
- `sequence`: .fa, .fasta, .seq, .dat (required)

**Default Flags:**
- `nproc=4`

**Strategy Options:**
- `quick`: true → `quick=True`
- `nproc`: nproc={}
- `resolution`: resolution={}

**Usage Hints:**
- REQUIRES phased MTZ (from phaser or refine output) - will not work with raw F/SIGF data
- If resolution is better than 2.0Å, set resolution=2.0 to speed up building
- Use the refined MTZ (with map coefficients) not the original data MTZ
- IMPORTANT: If a refined model exists, include it with model= for better results

### phenix.mtriage

**Description:** Map quality analysis for cryo-EM.

**Required Files:**
- `map`: .mrc, .ccp4, .map (required)
- `half_map_1`: .mrc, .ccp4 (optional)
- `half_map_2`: .mrc, .ccp4 (optional)

### phenix.dock_in_map

**Description:** Dock model into cryo-EM map.

**Required Files:**
- `model`: .pdb (required)
- `map`: .mrc, .ccp4, .map (required)
- `sequence`: .fa, .fasta, .seq, .dat (optional)

### phenix.autosol

**Description:** Experimental phasing (SAD/MAD) - use when MTZ has anomalous data and no model exists.

**Required Files:**
- `data`: .mtz, .sca, .hkl, .sdf (required)
- `sequence`: .fa, .fasta, .seq, .dat (required)

**Default Flags:**
- `nproc=4`

**Strategy Options:**
- `nproc`: nproc={}
- `quick`: true → `quick=True`
- `atom_type`: autosol.atom_type={}
- `wavelength`: autosol.wavelength={}
- `sites`: autosol.sites={}

**Usage Hints:**
- Use for SAD/MAD experimental phasing when data has anomalous signal
- Requires anomalous data (F+/F- or I+/I-) in the data file
- Set atom_type to the anomalous scatterer (e.g., 'Se', 'S', 'Zn', 'Fe')
- Alternative to MR when no homologous model exists

### phenix.map_to_model

**Description:** Build atomic model directly into cryo-EM map.

**Required Files:**
- `map`: .mrc, .ccp4, .map (required)
- `sequence`: .fa, .fasta, .seq, .dat (required)
- `model`: .pdb (optional)

**Default Flags:**
- `nproc=4`

**Strategy Options:**
- `nproc`: nproc={}
- `resolution`: resolution={}

**Usage Hints:**
- Use when you have a cryo-EM map and sequence but no model
- Alternative to docking a predicted model

### phenix.molprobity

**Description:** Comprehensive model validation (geometry, clashes, Ramachandran).

**Required Files:**
- `model`: .pdb, .cif (required)

**Usage Hints:**
- Run after refinement to check model quality
- Good values: Ramachandran outliers <0.5%, rotamer outliers <1%, clashscore <4, MolProbity score <1.5
- Poor values: any metric >2x the good threshold indicates problems

### phenix.model_vs_data

**Description:** Validate model against X-ray data (R-factors, map quality).

**Required Files:**
- `model`: .pdb, .cif (required)
- `data`: .mtz, .sca, .hkl, .sdf (required)

**Usage Hints:**
- Run after X-ray refinement to validate model vs data
- Reports R-work, R-free, and map-model correlation

### phenix.validation_cryoem

**Description:** Validate model against cryo-EM map.

**Required Files:**
- `model`: .pdb, .cif (required)
- `map`: .mrc, .ccp4, .map (required)

**Strategy Options:**
- `resolution`: resolution={}

**Usage Hints:**
- Run after cryo-EM refinement to validate model
- Reports map-model correlation and geometry statistics

### phenix.holton_geometry_validation

**Description:** Detailed geometry validation with energy-based scoring.

**Required Files:**
- `model`: .pdb, .cif (required)

**Usage Hints:**
- Provides detailed geometry analysis with energy scores
- Overall geometry energy should be close to expected value (~67 for ideal)
- Reports worst deviations in each category (bonds, angles, clashes, etc.)

---
## Stop Conditions

The agent will recommend stopping when any of these conditions are met:

### X-ray Success
- **Description:** R-free below dynamic target (resolution/10, bounded 0.20-0.30)
- **Formula:** `R-free < max(0.20, min(0.30, resolution/10)) - 0.02`

### Plateau Detection
- **Description:** Less than 0.5% improvement for 2+ consecutive cycles
- **Formula:** `improvement < 0.5% for last 2 cycles`

### Excessive Refinement
- **Description:** Too many consecutive refinement cycles
- **Formula:** `consecutive_refines >= 5`

### Cryo-EM Success
- **Description:** Map-model correlation above threshold
- **Formula:** `map_cc > 0.75`

---
## File Categorization

Files are categorized by extension and naming patterns:

| Category | Extensions/Pattern | Description |
|----------|-------------------|-------------|
| MTZ files | .mtz | X-ray data |
| Sequence files | .fa, .fasta, .seq | Protein sequence |
| Map files | .mrc, .ccp4, .map | Cryo-EM maps |
| PDB files | .pdb | Model coordinates |
| Phaser output | 'phaser' in filename | MR solution |
| Refined model | 'refine' in filename (not 'real_space') | X-ray refined |
| Predicted model | 'predict' or 'alphafold' in filename | AlphaFold prediction |
| Autobuild output | 'autobuild', 'buccaneer', 'build' in filename | Model building output |
| Ligand CIF | .cif without 'refine' | Ligand restraints |

---
## Key Thresholds

Important decision thresholds used by the agent:

| Threshold | Value | Description |
|-----------|-------|-------------|
| High Resolution (< 1.5Å) | autobuild=0.25, good=0.20, ligandfit=0.28 | Resolution-dependent R-free thresholds |
| Medium Resolution (1.5-2.5Å) | autobuild=0.35, good=0.25, ligandfit=0.35 | Resolution-dependent R-free thresholds |
| Low Resolution (> 2.5Å) | autobuild=0.38, good=0.30, ligandfit=0.38 | Resolution-dependent R-free thresholds |
| Plateau detection | < 0.3% for 3 cycles | Trigger plateau warning |
| Excessive refinement | >= 8 cycles | Too many consecutive refinement cycles |
| Cryo-EM success | CC > 0.70 | Map-model correlation target |

---

## Modifying Agent Behavior

To change agent behavior, edit these files:

| What to Change | File | Location |
|----------------|------|----------|
| **Thresholds & Defaults** | `agent/decision_config.json` | tiers, thresholds, program_rankings |
| Program commands/flags | `agent/command_templates.json` | file_slots, defaults, strategy_flags |
| Workflow state transitions | `agent/workflow_state.py` | _detect_xray_state(), _detect_cryoem_state() |
| Stop conditions | `agent/metrics_analyzer.py` | analyze_metrics_trend() |
| LLM prompts/guidance | `knowledge/prompts_hybrid.py` | SYSTEM_PROMPT, get_planning_prompt() |
| File categorization | `agent/workflow_state.py` | _categorize_files() |

### Configuration Priority

1. `decision_config.json` - Centralized thresholds and tier definitions (preferred)
2. `workflow_state.py` - State machine logic (uses config_loader)
3. `metrics_analyzer.py` - Stop condition detection

