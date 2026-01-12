# PHENIX AI Agent Decision-Making Architecture

## Overview

The PHENIX AI Agent is a LangGraph-based system that automates crystallographic and cryo-EM structure determination workflows. It operates as a cyclic graph that perceives the current state, plans the next action, builds commands, validates them, and executes them until the workflow completes.

## Architecture Components

### 1. Graph Structure

The agent operates as a state machine with the following nodes:

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│    ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌─────────┐ │
│ ──►│ PERCEIVE │───►│   PLAN   │───►│  BUILD   │───►│VALIDATE │ │
│    └──────────┘    └──────────┘    └──────────┘    └────┬────┘ │
│         ▲                                               │      │
│         │              ┌──────────┐                     │      │
│         │         ┌───►│ FALLBACK │◄────────────────────┤      │
│         │         │    └──────────┘   (on 3 failures)   │      │
│         │         │                                     │      │
│         │         ▼                                     ▼      │
│    ┌──────────┐◄──────────────────────────────────┌─────────┐  │
│    │  OUTPUT  │                                   │ EXECUTE │  │
│    └──────────┘◄──────────────────────────────────└─────────┘  │
│         │                                                      │
│         ▼                                                      │
│      [END or LOOP]                                             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 2. Node Responsibilities

| Node | Function |
|------|----------|
| **PERCEIVE** | Analyzes history, categorizes files, detects workflow state, extracts metrics |
| **PLAN** | Uses LLM to decide next program based on workflow state and recommendations |
| **BUILD** | Constructs command with correct files and parameters, applies Tier 2 defaults |
| **VALIDATE** | Checks program validity, file existence, and duplicate prevention |
| **EXECUTE** | Runs the PHENIX command (handled externally) |
| **FALLBACK** | Mechanical backup when LLM fails 3 times |
| **OUTPUT** | Formats results and determines if workflow should continue |

## Decision-Making Layers

### Layer 1: Workflow State Detection

The system first determines the current state based on:
- **History analysis**: What programs have run, their outcomes
- **File categorization**: Available PDBs, MTZs, maps, sequences
- **Metrics extraction**: R-free, R-work, resolution, map CC

#### X-ray Crystallography States

| State | Description | Typical Valid Programs |
|-------|-------------|----------------------|
| `xray_initial` | Starting point | xtriage |
| `xray_analyzed` | After xtriage | predict_and_build, phaser |
| `xray_has_prediction` | Have AlphaFold model | process_predicted_model |
| `xray_model_processed` | Prediction trimmed | phaser |
| `xray_has_phases` | After experimental phasing | autobuild |
| `xray_has_model` | Have placed model | refine |
| `xray_refined` | After refinement | refine, autobuild, molprobity, ligandfit, STOP |
| `xray_has_ligand` | Ligand fitted | pdbtools |
| `xray_combined` | Ligand + protein combined | refine, molprobity, STOP |

#### Cryo-EM States

| State | Description | Typical Valid Programs |
|-------|-------------|----------------------|
| `cryoem_initial` | Starting point | mtriage |
| `cryoem_analyzed` | After mtriage | predict_and_build, dock_in_map |
| `cryoem_has_model` | Model in map | real_space_refine |
| `cryoem_refined` | After RSR | real_space_refine, molprobity, STOP |

### Layer 2: Tiered Decision System

Decisions are governed by a three-tier configuration system defined in `decision_config.json`:

#### Tier 1: Hard Constraints (Cannot Override)

These are absolute rules that must be followed:

| Constraint | Description |
|------------|-------------|
| `xtriage_before_mr` | Must run xtriage before molecular replacement |
| `phaser_requires_model` | Phaser needs a search model |
| `autobuild_requires_phases` | Autobuild needs phased MTZ |
| `refine_requires_model_and_data` | Refine needs both PDB and MTZ |
| `no_duplicate_commands` | Cannot run identical command twice |

#### Tier 2: Strong Defaults (Override with Warning)

Applied automatically; LLM can override but a warning is logged:

| Default | Condition | Value | Override Warning |
|---------|-----------|-------|------------------|
| `generate_rfree_flags` | First refinement | True | May compromise R-free statistics |
| `use_existing_rfree_flags` | Subsequent refinements | True | Generating new flags loses continuity |
| `add_waters` | R-free < 0.35, resolution ≤ 3.0Å | True | Skipping water building prematurely |
| `twin_law` | Twinning detected | Include twin law | Ignoring detected twinning |

#### Tier 3: Soft Guidance (Suggestions Only)

Presented to LLM but no enforcement or warnings:

- Consider anisotropic ADP for resolution < 1.5Å
- Consider riding hydrogens for resolution < 1.2Å
- Consider stopping when R-free below target
- Suggest validation after 3+ refinement cycles

### Layer 3: Program Rankings

For each state, programs are ranked by priority with conditions:

```json
"xray_refined": {
  "rankings": [
    {"program": "phenix.ligandfit", "priority": 1, "condition": "has_ligand and model_is_good"},
    {"program": "phenix.autobuild", "priority": 2, "condition": "r_free_high and has_sequence"},
    {"program": "phenix.molprobity", "priority": 3, "condition": "suggest_validation"},
    {"program": "phenix.refine", "priority": 4, "condition": "r_free_above_target"},
    {"program": "STOP", "priority": 5, "condition": "always"}
  ]
}
```

### Layer 4: Validation Gate

A critical safety mechanism that **prevents stopping without validation**:

**Conditions for requiring validation before STOP:**
1. R-free is below success threshold (target - 0.02), OR
2. 3+ refinement cycles have been completed, OR
3. R-free is below the good model threshold

**When validation is required:**
- `STOP` is removed from `valid_programs`
- `phenix.molprobity` is added at highest priority
- Reason includes "VALIDATION REQUIRED BEFORE STOPPING"

**After validation completes:**
- `STOP` is restored to `valid_programs` (at position 0 if model is good)
- Workflow can terminate

### Layer 5: LLM Planning

The LLM receives a structured prompt containing:
1. Current workflow state and valid programs
2. Ranked program recommendations with reasons
3. Tier 2 defaults that will be applied
4. Hard constraints that cannot be violated
5. Soft suggestions to consider
6. Recent history and metrics trend

The LLM outputs a JSON decision:
```json
{
  "program": "phenix.refine",
  "reasoning": "Model quality improving, continue refinement...",
  "files": {"data": "data.mtz", "model": "model.pdb"},
  "strategy": {"ordered_solvent": true},
  "stop": false
}
```

### Layer 6: Command Building

The BUILD node:
1. Validates LLM's program choice against `valid_programs`
2. Overrides if LLM tries to STOP when STOP isn't allowed
3. Corrects file paths (maps basenames to full paths)
4. Applies Tier 2 defaults if LLM didn't specify
5. Logs any overrides with warnings
6. Constructs the final command string

### Layer 7: Validation & Fallback

**Validation checks:**
1. Program is in valid_programs for current state
2. All referenced files exist
3. Command is not a duplicate of previous cycle

**On validation failure:**
- Increment attempt counter
- Return to PLAN for retry
- After 3 failures, go to FALLBACK

**FALLBACK behavior:**
1. Try each valid program in order (skip duplicates)
2. If molprobity is valid (validation required), run it
3. Last resort: STOP with reason "all_commands_duplicate"

## Metrics and Stop Conditions

### Resolution-Dependent Thresholds

| Resolution | Autobuild R-free | Good Model R-free | Success Threshold |
|------------|------------------|-------------------|-------------------|
| < 1.5Å | 0.30 | 0.20 | 0.18 |
| 1.5-2.5Å | 0.35 | 0.25 | 0.23 |
| 2.5-3.5Å | 0.38 | 0.28 | 0.26 |
| > 3.5Å | 0.38 | 0.30 | 0.28 |

### Stop Conditions

The workflow stops when:
1. **SUCCESS**: R-free below success threshold AND validation done
2. **PLATEAU**: < 0.3% improvement for 3+ consecutive cycles
3. **EXCESSIVE**: 8+ consecutive refinement cycles
4. **USER REQUEST**: LLM chooses STOP (only when allowed)
5. **ALL DUPLICATES**: No non-duplicate commands available

## Data Flow

```
┌─────────────┐
│   History   │────┐
│  (cycles)   │    │
└─────────────┘    │
                   ▼
┌─────────────┐  ┌─────────────────┐  ┌──────────────┐
│  Available  │─►│  WORKFLOW STATE │─►│ valid_programs│
│    Files    │  │    DETECTION    │  │ ranked_programs│
└─────────────┘  └─────────────────┘  │ recommendations│
                   │                   └──────────────┘
                   ▼                          │
              ┌─────────────┐                 │
              │   METRICS   │                 │
              │  ANALYZER   │                 │
              └─────────────┘                 │
                   │                          │
                   ▼                          ▼
              ┌─────────────┐          ┌─────────────┐
              │ should_stop │          │   PROMPT    │
              │   reason    │─────────►│ CONSTRUCTION│
              │   trend     │          └─────────────┘
              └─────────────┘                 │
                                             ▼
                                       ┌─────────────┐
                                       │     LLM     │
                                       │   DECISION  │
                                       └─────────────┘
                                             │
                                             ▼
                                       ┌─────────────┐
                                       │   COMMAND   │
                                       │   BUILDER   │
                                       └─────────────┘
```

## Configuration Files

| File | Purpose |
|------|---------|
| `decision_config.json` | Tiered rules, thresholds, program rankings |
| `command_templates.json` | Program file slots, defaults, strategy flags |
| `prompts_hybrid.py` | LLM prompt templates and formatting |

---

## Limitations

### 1. Resolution Value Extraction

The system sometimes extracts incorrect resolution values (e.g., 47.31Å instead of 2.1Å). This appears to be a parsing issue where the low-resolution limit is captured instead of the high-resolution limit. This affects threshold calculations and decision-making.

### 2. LLM Reliability

- The LLM may choose STOP even when told not to (requires override logic)
- The LLM may generate duplicate commands, requiring retry cycles
- Different LLM providers (Google, OpenAI) may behave differently
- LLM reasoning doesn't always match its actual decision

### 3. File Selection

- The LLM sometimes selects wrong files (e.g., original MTZ instead of refined MTZ with R-free flags)
- File path correction relies on basename matching, which can fail with duplicate names
- No semantic understanding of which file is "best" for a given purpose

### 4. State Detection Edge Cases

- Intermediate states during ligand fitting workflow can be ambiguous
- The system may not correctly detect when autobuild has been run if output files have non-standard names
- Cryo-EM workflows are less thoroughly tested than X-ray

### 5. Metrics Parsing

- Autobuild table format requires special parsing
- Some programs don't output metrics in a parseable format
- The "Final Quality Metrics Report" may show duplicate/outdated values

### 6. Validation Limitations

- Molprobity is the only validation tool currently integrated
- No integration with wwPDB validation server
- Validation metrics aren't used to guide further refinement decisions

### 7. Single-Threaded Execution

- Commands run sequentially, not in parallel
- Long-running jobs (autobuild, predict_and_build) block the workflow
- No checkpoint/resume capability for interrupted long jobs

### 8. Error Recovery

- Some errors (e.g., R-free flag mismatch) require user intervention
- The system may get stuck in retry loops
- No automatic rollback to previous good state

---

## Future Improvements

### Short-Term (Implementation Ready)

1. **Fix Resolution Parsing**
   - Add specific regex for high-resolution limit
   - Cross-validate resolution from multiple sources (xtriage, refine)
   - Reject obviously wrong values (> 10Å for typical macromolecular data)

2. **Improve File Selection**
   - Track file provenance (which program created which file)
   - Prefer files from most recent successful cycle
   - Implement explicit MTZ label checking for R-free flags

3. **Better Duplicate Handling**
   - Vary strategy parameters (e.g., cycles, optimization target) to avoid duplicates
   - Track what strategies have been tried
   - Smarter differentiation of "same command" vs "same intent"

4. **Enhanced Metrics Display**
   - De-duplicate metrics in final report
   - Show metrics progression graph
   - Highlight improvements/regressions

### Medium-Term (Design Required)

5. **Multi-Model Support**
   - Track multiple alternative models
   - Compare models and select best
   - Support ensemble refinement

6. **Ligand Workflow Improvements**
   - Better ligand file detection and handling
   - Integration with eLBOW for restraint generation
   - Support for multiple ligands

7. **Advanced Validation**
   - Integrate wwPDB OneDep validation
   - Use validation results to guide refinement strategy
   - Set validation-based stop conditions

8. **Checkpoint/Resume**
   - Save state after each successful cycle
   - Resume from last checkpoint after interruption
   - Support for long-running jobs

### Long-Term (Research Required)

9. **Adaptive Strategy Selection**
   - Learn from successful/failed strategies
   - Per-resolution and per-space-group optimization
   - Integration with crystallographic knowledge bases

10. **Parallel Execution**
    - Run independent tasks in parallel
    - Distributed computing support
    - Job queue management

11. **Interactive Mode**
    - Allow user intervention at decision points
    - Support for "suggest but don't execute" mode
    - Better explanation of decisions

12. **Expanded Workflow Support**
    - SAD/MAD phasing workflows
    - Neutron crystallography
    - Time-resolved crystallography
    - Microcrystal electron diffraction (MicroED)

13. **Quality Prediction**
    - Predict expected R-free improvement before running
    - Estimate time to convergence
    - Early stopping based on predicted diminishing returns

---

## Summary

The PHENIX AI Agent implements a sophisticated multi-layer decision system:

1. **Perception** → Understand current state from history and files
2. **Rules** → Apply hard constraints and resolution-dependent thresholds
3. **Recommendations** → Rank programs by priority with conditions
4. **Gates** → Enforce validation before stopping
5. **LLM Planning** → Intelligent decision-making with full context
6. **Defaults** → Apply best practices automatically
7. **Validation** → Catch errors before execution
8. **Fallback** → Mechanical backup when LLM fails

This architecture balances flexibility (LLM can adapt to unusual situations) with safety (hard constraints and validation gates prevent bad decisions). The tiered configuration system makes it easy to adjust thresholds and rules without code changes.

The system successfully automates routine structure determination while maintaining the quality standards expected in crystallography.
