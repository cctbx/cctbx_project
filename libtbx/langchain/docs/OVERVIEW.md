# PHENIX AI Agent V2 - Technical Overview

## Introduction

The PHENIX AI Agent is an intelligent automation system for macromolecular
structure determination. It combines:

- **Goal-Directed Planning** - Decomposes structure determination into
  phases with success criteria, retreat logic, and hypothesis testing (v114)
- **LLM Decision Making** - Uses Claude/Gemini to interpret context and
  make intelligent program choices
- **Rules-Based Validation** - Enforces workflow constraints and validates
  all decisions
- **YAML Configuration** - All domain knowledge externalized to editable
  config files
- **Structured Event Logging** - Transparent decision tracking at multiple
  verbosity levels

The agent operates at two levels: a **strategic planner** (v114,
activated by `thinking_level=expert`) that produces a multi-phase
plan at session start and evaluates progress at stage gates, and
a **reactive execution engine** (v112+) that handles per-cycle
program selection, file management, error recovery, and safety
checks. The strategic planner communicates through directives —
the same interface a human user would use — so the reactive
agent's safety checks always apply.

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              USER INTERFACE                                  │
│  phenix.ai_agent original_files="data.mtz seq.fa" project_advice="..."      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            programs/ai_agent.py                              │
│  • Parse parameters, initialize session                                      │
│  • Generate plan at session start (v114)                                    │
│  • Choose LocalAgent or RemoteAgent                                          │
│  • After each cycle: gate evaluation + plan update (v114)                   │
│  • Display results with EventFormatter                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                    ┌─────────────────┴─────────────────┐
                    ▼                                   ▼
         ┌──────────────────┐                ┌──────────────────┐
         │   LocalAgent     │                │   RemoteAgent    │
         │ (phenix_ai/)     │                │ (phenix_ai/)     │
         └────────┬─────────┘                └────────┬─────────┘
                  └─────────────────┬─────────────────┘
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       STRATEGIC PLANNER (v114)                                │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌────────────────┐  │
│  │ Plan         │  │ Gate         │  │ Structure    │  │ Explanation    │  │
│  │ Generator    │  │ Evaluator    │  │ Model        │  │ Engine         │  │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘  └──────┬─────────┘  │
│         │                 │                 │                  │            │
│         └─────── directives + advice ───────┘                  │            │
│                           │                                    │            │
│  ┌────────────────────┐   │  ┌───────────────────┐             │            │
│  │ Hypothesis         │───┘  │ Validation        │─────────────┘            │
│  │ Evaluator          │      │ History           │                          │
│  └────────────────────┘      └───────────────────┘                          │
└───────────────────────────────┼─────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         phenix_ai/run_ai_agent.py                            │
│  • Build LangGraph state                                                     │
│  • Execute graph nodes (perceive → think → plan → build → validate → output)│
│  • Retry loop with fallback on persistent validation failures               │
│  • Build response with events                                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            LangGraph Pipeline                                │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐                  │
│  │ perceive │ → │  think   │ → │   plan   │ → │  build   │                  │
│  └────┬─────┘   └──────────┘   └────┬─────┘   └────┬─────┘                  │
│       │                             │   ▲          │                         │
│       │                             │   └──────────│────── (retry < 3)      │
│       │                             │              │                         │
│       │                             │         ┌────┴──────┐                  │
│       │                             │         │ validate  │                  │
│       │                             │         └────┬──────┘                  │
│       │                             │              │                         │
│       │                             │         ┌────┴─────┐                   │
│       │                             │         │ fallback │ (retry >= 3)      │
│       │                             │         └────┬─────┘                   │
│       │                             │              │                         │
│       │                             ▼              ▼                         │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                       ┌──────────┐                                  │    │
│  │                       │  output  │ → END                            │    │
│  │                       └──────────┘                                  │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Core Components

### 1. Graph Nodes (`agent/graph_nodes.py`)

The agent uses a LangGraph pipeline with seven nodes:

| Node | Purpose | Key Actions |
|------|---------|-------------|
| **perceive** | Understand current state | Categorize files, detect workflow state, extract metrics, analyze trends |
| **think** | Expert reasoning (v113+) | Analyze logs with domain expertise; update Structure Model from validation, xtriage, phaser results (v114); generate hypothesis prompts (v114); inject strategic guidance into plan |
| **plan** | Decide next action | Call LLM or rules engine, validate against workflow, check user directives and plan-generated directives (v114) |
| **build** | Generate command | Select files using BestFilesTracker, build command string |
| **validate** | Check output | Sanity checks, red flag detection; retry up to 3× on failure |
| **fallback** | Last resort | If validate fails 3×, use mechanical/rules-based program selection |
| **output** | Format response | Package decision, events, and command into response |

The graph has conditional routing: perceive can skip to output (on red
flag abort), validate can loop back to plan (retry) or route to fallback
(max retries exceeded).

**Between cycles** (in `ai_agent.py`, not in the LangGraph pipeline):
the Gate Evaluator (v114) runs after each cycle to check phase success
criteria, advance/retreat/skip phases, update the plan, and trigger
`advice_changed` on plan revisions. Stage transition summaries are
generated by the Explanation Engine at each gate evaluation.

### 2. Workflow Engine (`agent/workflow_engine.py`)

Interprets `workflows.yaml` to determine:
- **Current phase**: X-ray has 9 phases (analyze, obtain_model, molecular_replacement, experimental_phasing, build_from_phases, refine, combine_ligand, validate, complete); cryo-EM has 9 phases (analyze, obtain_model, dock_model, check_map, optimize_map, ready_to_refine, refine, validate, complete)
- **Valid programs** for current phase
- **Transition conditions**
- **Stop criteria** (target reached, plateau detected)

### 3. Command Builder (`agent/command_builder.py`)

Unified command generation that:
- Reads program definitions from `programs.yaml`
- Selects appropriate input files using BestFilesTracker
- Applies strategy flags and defaults
- Tracks WHY each file was selected (for transparency)
- **Recovery param injection**: Recovery-sourced strategy entries (e.g.,
  `obs_labels` from ambiguous data label recovery) are appended after
  template-based command assembly, since `build_command` only emits params
  matching known `strategy_flags` keys
- **Content-based guards**: rejects small-molecule PDB files from model slots
  (`_pdb_is_small_molecule`) and protein models from ligand slots
  (`_pdb_is_protein_model`)
- **Word-boundary exclude patterns**: `matches_exclude_pattern()` prevents
  false positives (e.g., "noligand" no longer matches "ligand")
- **LLM validation**: LLM file assignments checked against slot `exclude_patterns`
  and content guards before acceptance
- **Diagnostics**: `_last_missing_slots` records which required slots could not
  be filled, enabling specific fallback error messages

### 4. Event System (`agent/event_log.py`, `agent/event_formatter.py`)

Structured logging for transparency:

| Event Type | Verbosity | Description |
|------------|-----------|-------------|
| `cycle_start` | quiet | New cycle beginning |
| `cycle_complete` | quiet | Cycle finished |
| `state_detected` | normal | Workflow state determined |
| `metrics_extracted` | normal | R-free, CC, resolution parsed |
| `metrics_trend` | normal | Improvement/plateau analysis |
| `sanity_check` | normal | Red flag or warning detected |
| `program_selected` | normal | Decision with full reasoning |
| `program_modified` | normal | Program changed by rules/validation |
| `stop_decision` | normal | Whether to continue |
| `directive_applied` | normal | User directive enforced |
| `expert_assessment` | normal | Thinking agent analysis (v113) |
| `user_request_invalid` | quiet | User requested unavailable program |
| `files_selected` | verbose | File selection with reasons |
| `file_scored` | verbose | Individual file scoring detail |
| `command_built` | normal | Final command |
| `thought` | verbose | LLM chain-of-thought/reasoning traces |
| `error` | quiet | Error occurred |
| `warning` | quiet | Non-fatal warning |
| `debug` | verbose | Internal debug information |

### 5. File Categorization (`agent/graph_nodes.py`, `knowledge/file_categories.yaml`)

Semantic file classification using rules from `file_categories.yaml`:

| Category | Description | Examples |
|----------|-------------|----------|
| `data_mtz` | Reflection data (Fobs, R-free) | data.mtz, refine_001_data.mtz |
| `map_coeffs_mtz` | Map coefficients (calculated phases) | refine_001_001.mtz, denmod_map_coeffs.mtz |
| `sequence` | Sequence files | seq.fa, protein.fasta |
| `model` | **Positioned** models | refine_001.pdb, phaser_output.pdb |
| `search_model` | Templates/predictions | template.pdb, alphafold.pdb |
| `full_map` | Cryo-EM maps | map.mrc, emd_1234.ccp4 |
| `half_map` | Half-maps | half1.mrc, map_half_a.mrc |
| `ligand` | Ligand files | ligand.pdb, ATP.cif |

**Semantic distinction**: `model` = already positioned in crystal/map; `search_model` = template not yet placed.

**Post-processing content guards** validate YAML categorizer output:
- **Half-map guard**: Files in `half_map` without "half" in name → reclassified to `full_map`
- **Ligand guard** (v112.74): PDB files in `ligand_pdb` that are actually protein
  models (>150 atoms, majority ATOM records) → rescued to `model`.  Catches
  false positives from broad YAML patterns matching names like `1aba.pdb`.
- **MTZ safety net**: Cross-checks all MTZ files against `classify_mtz_type()` regex

### 6. Best Files Tracker (`agent/best_files_tracker.py`)

Tracks the best file of each type across cycles:
- Scores files based on metrics (R-free, resolution, cycle number)
- **Dual MTZ tracking**:
  - `data_mtz`: Locks after first R-free flags generated (consistency for refinement)
  - `map_coeffs_mtz`: Always prefers most recent (maps improve with refinement)
- Provides `best_files["model"]`, `best_files["data_mtz"]`, `best_files["map_coeffs_mtz"]` to CommandBuilder
- **Supplemental file discovery**: Session load (`_rebuild_best_files_from_cycles`)
  and live cycle completion (`record_result`) both call `_find_missing_outputs` to
  discover companion files (e.g., `refine_001.mtz` from `refine_001_data.mtz`) and
  evaluate them through the tracker. This ensures `map_coeffs_mtz` is populated even
  when the client only tracked a subset of output files.
- **MTZ classification**: `file_utils.classify_mtz_type()` uses regex
  `(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$` to correctly identify standard refinement
  output as map coefficients (not raw data)
- **MTZ categorization safety net** (v112.71): After both YAML and hardcoded
  categorization, `_categorize_files()` cross-checks every MTZ file against the
  authoritative `classify_mtz_type()` regex. Catches three failure modes:
  misclassified files (moved to correct category), dual-categorized files
  (removed from `data_mtz` when also in `map_coeffs_mtz`), and missing
  subcategories (added to `refine_map_coeffs` etc.). Logs `WARNING` when
  corrections are made, making future occurrences immediately diagnosable.

### 7. User Directives (`agent/directive_extractor.py`)

Parses natural language guidance from `project_advice`:
- **Stop conditions**: "stop after one refinement", "stop when R-free < 0.25"
- **Workflow preferences**: "skip autobuild", "use SAD phasing"
- **Task focus**: "focus on ligand fitting"

**Stop condition semantics:**

| Condition | Type | Where checked | Behavior |
|-----------|------|---------------|----------|
| `after_cycle` | Hard stop | PERCEIVE | Immediately stops at cycle N |
| `r_free_target` | Hard stop | PERCEIVE | Immediately stops when R-free ≤ target |
| `map_cc_target` | Hard stop | PERCEIVE | Immediately stops when map CC ≥ target |
| `after_program` | Minimum-run guarantee | PLAN | Suppresses auto-stop until target program has run, but the **LLM decides** when to actually stop |

**Why `after_program` is not a hard stop (v112.78):** The directive extractor
can only name one program in `after_program`.  For multi-goal requests like
"improve map, get symmetry and map correlation," it picks one program (e.g.,
`map_symmetry`) and goals beyond that point would be silently dropped if
`after_program` triggered an immediate stop.  By letting the LLM decide, all
goals in the user's advice are honored.

### 8. Safety Checks (`agent/sanity_checker.py`, `agent/directive_extractor.py`)

Multiple layers of validation to prevent errors:

| Layer | Location | Examples |
|-------|----------|----------|
| **Sanity Checks** | Pre-execution | No data for workflow, model not positioned, repeated failures |
| **Directive Validation** | Post-LLM | Invalid program names, conflicting stop conditions |
| **Workflow Validation** | State machine | Invalid stage transitions, wrong experiment type |
| **Post-Processing** | After extraction | Ligand workflow conflict resolution |

Key sanity issues that trigger abort:
- `no_data_for_workflow` - Missing data_mtz (X-ray) or map (cryo-EM)
- `search_model_not_positioned` - Trying to refine before MR/docking
- `no_model_for_refine` - No model available for refinement
- `repeated_failures` - Same error 3+ times

See `docs/SAFETY_CHECKS.md` for the complete auto-generated list.

### 9. RAG Pipeline (`rag/`, `utils/query.py`, `analysis/analyzer.py`)

The agent uses Retrieval-Augmented Generation to ground LLM responses in
PHENIX documentation. The pipeline has three stages:

1. **Retrieval**: A Chroma vector store (built from PHENIX docs via
   `run_build_db.py`) returns the top 20 candidate documents using
   embedding similarity.

2. **Reranking**: FlashRank, a local cross-encoder model
   (`ms-marco-MiniLM-L-12-v2`, ~34MB), reranks the 20 candidates and
   selects the top 8 most relevant. This runs on CPU with no API key —
   the `flashrank` package must be installed where analysis runs locally.

3. **Generation**: The reranked documents are passed as context to the
   LLM (Google or OpenAI) which generates the final response.

The reranking retriever is created in `rag/retriever.py` via
`create_reranking_retriever()` and used by both log analysis
(`analysis/analyzer.py`) and documentation queries (`utils/query.py`).

### 10. Strategic Planner (v114)

The goal-directed layer sits above the reactive agent and communicates
through the directives system — the same interface a human user would
use. The reactive agent is unchanged; every safety check still applies.

Activated by `thinking_level=expert`. The four thinking levels are:

| Level | THINK node | Structure Model | Planning layer |
|-------|-----------|-----------------|---------------|
| `none` | off | off | off |
| `basic` | log analysis + LLM | off | off |
| `advanced` (default) | validation + KB + LLM | updated | off |
| `expert` | validation + KB + LLM | updated | **on** |

The Structure Model is updated at both `advanced` and `expert`
levels (it feeds the expert assessment display). The planning
layer (plan generation, gate evaluation, hypothesis testing,
reports) only activates at `expert`.

#### 10a. Structure Model (`agent/structure_model.py`)

Maintains a running understanding of the specific structure being solved.
Updated every cycle from validation results (ground truth, not LLM
reasoning). Persists across cycles and session resume.

| Tracking Area | Content | Source |
|---------------|---------|--------|
| Data characteristics | Resolution, space group, twinning, anomalous | xtriage log |
| Model state | Chains, ligands, waters, problem regions | validation_inspector |
| Progress | R-free trajectory with annotations | log analysis + validation |
| Strategy blacklist | Tried-and-failed strategies | gate evaluator retreats |
| Hypotheses | Proposed/active/confirmed/refuted | hypothesis evaluator |

Key methods: `update_from_validation()`, `update_from_xtriage()`,
`update_from_phaser()`, `get_summary(detail_level=)`,
`get_current_problems()`, `blacklist_strategy()`, `is_blacklisted()`.

The `ValidationHistory` (`agent/validation_history.py`) stores per-cycle
validation snapshots with `get_metric_series()` for trend analysis and
`get_phase_start_metrics()` for the monotonic progress gate.

#### 10b. Plan Generator (`agent/plan_generator.py`)

Produces a multi-phase strategy at session start by selecting and
customizing pre-defined templates from `knowledge/plan_templates.yaml`.

Twelve templates cover common scenarios: `mr_refine`,
`mr_refine_ligand`, `mr_refine_lowres`, `mr_refine_highres`,
`mr_refine_twinned`, `predict_refine`, `predict_refine_ligand`,
`mr_sad`, `sad_phasing`, `sad_phasing_ligand`, `cryoem_refine`,
`cryoem_refine_ligand`. Template selection is deterministic
(experiment type + available files + resolution + anomalous
atoms). The LLM only customizes parameters within template
bounds.

Plans are represented as `StructurePlan` objects (`knowledge/plan_schema.py`)
containing a list of `StageDef` phases, each with programs, success
criteria, gate conditions, fallbacks, and skip conditions.

The plan communicates with the reactive agent through:
- `plan_to_directives()` — translates current phase to `prefer_programs`,
  `after_program`, `program_settings`
- `compute_hash()` — Strategy Hash fingerprint; changes trigger
  `advice_changed` in the reactive agent

#### 10c. Gate Evaluator (`agent/gate_evaluator.py`)

After each cycle, evaluates phase progress against the plan's success
criteria. Purely deterministic (no LLM). Returns one of: continue,
advance, retreat, fallback, skip, or stop.

**Success hysteresis**: a 1.5% buffer prevents oscillation on noisy
metrics (e.g., R-free 0.349 → 0.351 → 0.349 doesn't trigger repeated
advance/continue decisions).

**Retreat logic** with five anti-oscillation safeguards:

| Safeguard | Prevention |
|-----------|------------|
| Strategy Blacklist | Never re-try a failed strategy |
| Retreat counter | Max 2 retreats per phase |
| Monotonic progress gate | Only retreat if worse than phase start |
| Retreat cooldown | 2+ cycles between retreats |
| Retreat depth limit | Max 1 phase backwards |

#### 10d. Hypothesis Engine (`agent/hypothesis_evaluator.py`)

Structured hypothesis formation and testing, integrated into the THINK
node. Enforces a single active hypothesis budget (one testing/pending
at a time) to prevent confounded multi-variable experiments.

Lifecycle: proposed → testing → pending → confirmed/refuted/abandoned.
`test_cycles_remaining` provides verification latency (prevents premature
refutation on unstabilized models). Confirmed hypotheses are re-validated
each cycle and demoted if evidence weakens (B-factor drift, RSCC drop).

#### 10e. Explanation Engine (`knowledge/explanation_prompts.py`)

Produces crystallographer-level commentary at three detail levels:

| Function | When | LLM? |
|----------|------|------|
| `generate_cycle_commentary()` | Every cycle | No (template) |
| `generate_stage_summary()` | Stage transitions | Yes |
| `generate_final_report()` | Session completion | Yes |
| `generate_stopped_report()` | Early stop | Yes |

#### 10f. Model Placement Gate (v114.1)

Prevents destructive program selection on models that already fit the
data. This is a single mechanism that addresses the class of bugs where
the LLM runs phaser/autosol on a pre-solved structure, producing a
worse result than simple refinement.

Detection (in `_detect_model_placement`, `programs/ai_agent.py`):
- `model_vs_data` with CC > 0.3 → placed
- `refine` with R-free < 0.50 → placed
- `real_space_refine` with map_cc > 0.3 → placed

When placement is confirmed:
1. `session.data["model_is_placed"] = True` — locked, survives resume
2. `valid_programs` filtered: phaser, autosol, predict_and_build removed
3. Plan fast-forward: MR/phasing stages marked "skipped" (⊘)
4. Conflict warning if user advice mentions MR/phaser

#### 10g. Display Data Model (`agent/display_data_model.py`)

Unified data provider for Results tab, Progress tab, and HTML report.
Eliminates duplication across three display views:

- `outcome_status`: determined / stopped / incomplete
- `outcome_message`: one-line human-readable summary
- `final_metrics`: best R-free/CC from structure model or cycle scan
- `rfree_trajectory` / `cc_trajectory`: for SVG charts
- `timeline`: compact cycle entries
- `format_cycle_compact()`: one-line-per-cycle formatting

#### 10h. HTML Report Generator (`knowledge/html_report_template.py`)

Template-based HTML report using DisplayDataModel. No LLM call.
Includes inline SVG R-free/CC trajectory chart with retreat markers.
Generated at session completion for all modes (not just expert).

---

## Configuration Files

### programs.yaml

Defines each PHENIX program the agent can run:

```yaml
phenix.refine:
  description: "Crystallographic refinement"
  category: refinement
  experiment_types: [xray]
  
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
        priority_patterns: [refine]  # Prefer refined models
      data_mtz:
        extensions: [.mtz, .sca, .hkl, .sdf]
        flag: ""
  
  outputs:
    files:
      - pattern: "*_refine_*.pdb"
        type: model
      - pattern: "*_data.mtz"
        type: data_mtz
      - pattern: "*_refine_*.mtz"
        type: map_coeffs_mtz
    metrics:
      - r_free
      - r_work
      - bonds_rmsd
      - angles_rmsd
  
  log_parsing:
    r_free:
      pattern: 'R-free\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-free"
    r_work:
      pattern: 'R-work\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "R-work"
    bonds_rmsd:
      pattern: '[Bb]onds?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Bonds RMSD"
    angles_rmsd:
      pattern: '[Aa]ngles?\s*(?:RMSD)?\s*[=:]\s*([0-9.]+)'
      type: float
      display_name: "Angles RMSD"
```

### workflows.yaml

Defines workflow state machines:

```yaml
xray:
  phases:
    analyze:
      description: "Analyze data quality"
      programs:
        - phenix.xtriage
      transitions:
        on_complete: obtain_model
    
    obtain_model:
      description: "Get initial model"
      programs:
        - program: phenix.predict_and_build
          preferred: true
          conditions:
            - has: sequence
            - not_done: predict_full  # Don't re-run after full workflow completes
        - program: phenix.phaser
          conditions:
            - has: search_model
            - not_done: phaser  # MR should only run once
        - program: phenix.autosol
          conditions:
            - has: sequence
            - has: anomalous        # Requires anomalous signal detected by xtriage
            - not_done: autosol
      transitions:
        on_complete: refine
        if_predict_only: molecular_replacement  # Stepwise mode
    
    # MR-SAD: after phaser places model, use anomalous signal with autosol
    experimental_phasing:
      description: "MR-SAD phasing with placed model"
      programs:
        - program: phenix.autosol
          conditions:
            - not_done: autosol
      transitions:
        on_complete: build_from_phases
    
    refine:
      description: "Improve model"
      programs:
        - program: phenix.refine
        - program: phenix.autobuild
          conditions:
            - r_free: "> autobuild_threshold"
            - has: sequence
            - not_done: autobuild
        - program: phenix.autobuild_denmod
          conditions:
            - has: ligand_file
            - has: sequence
            - not_done: autobuild_denmod
            - refine_count: "> 0"    # Must refine first
          hint: "Run density modification before ligand fitting"
        - program: phenix.ligandfit
          conditions:
            - has: ligand_file
            - not_done: ligandfit
            - r_free: "< 0.35"
            - refine_count: "> 0"    # Must refine first
        - program: phenix.polder
          conditions:
            - has: model
            - has: data_mtz
      repeat:
        max_cycles: 4
        until:
          any:
            - r_free: "< target_r_free"
            - condition: plateau
              cycles: 2
              threshold: 0.005
      transitions:
        on_ligandfit: combine_ligand
        on_target_reached: validate
        on_plateau: validate
        on_max_cycles: validate

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
```

#### Program Execution Controls

The workflow uses two mechanisms to prevent programs from running repeatedly:

**1. `not_done` conditions** (in workflows.yaml)

Programs can specify `not_done: <flag>` to prevent re-runs:

| Flag | Program | Meaning |
|------|---------|---------|
| `predict_full` | predict_and_build (X-ray) | Full workflow (prediction+MR+building) completed |
| `predict` | predict_and_build (cryo-EM) | Any prediction completed |
| `process_predicted_model` | process_predicted_model | Model processed for MR |
| `phaser` | phaser | Molecular replacement completed |
| `dock` | dock_in_map | Model docked in map |
| `autobuild` | autobuild | Model building completed |
| `autobuild_denmod` | autobuild_denmod | Density modification completed |
| `autosol` | autosol | Experimental phasing completed |
| `ligandfit` | ligandfit | Ligand fitted |
| `resolve_cryo_em` | resolve_cryo_em | Map optimization completed |
| `map_sharpening` | map_sharpening | Map sharpening completed |
| `map_to_model` | map_to_model | De novo model building completed |
| `map_symmetry` | map_symmetry | Map symmetry analysis completed |

**2. `done_tracking` blocks** (in programs.yaml)

Each program's `done_tracking` block defines its workflow done flag and tracking strategy:

- **`strategy: "set_flag"`** (default) — sets a boolean done flag on success. Most programs use this.
- **`strategy: "run_once"`** — sets the done flag AND filters the program from the valid list after first successful run. Used by `phenix.xtriage`, `phenix.mtriage`, and `phenix.map_symmetry`.
- **`strategy: "count"`** — sets the done flag AND increments a counter (e.g., `refine_count`). Used by `phenix.refine`, `phenix.real_space_refine`, and `phenix.phaser`.

All detection is driven by `history_detection.markers` (substring matching). Programs like `phenix.refine` use `exclude_markers: ["real_space"]` to prevent false matches with `phenix.real_space_refine`. The only program requiring Python-only tracking is `phenix.predict_and_build`, which cascades flags across programs.

**Programs that run multiple times** (intentionally):
- `phenix.refine` / `phenix.real_space_refine` - Iterative refinement
- `phenix.molprobity` - Validation after each cycle
- `phenix.polder` - Can run for different sites

### metrics.yaml

Quality thresholds and display configuration:

```yaml
metrics:
  r_free:
    display_name: "R-free"
    format: "{value:.4f}"
    direction: lower_is_better
    thresholds:
      excellent: 0.20
      good: 0.25
      acceptable: 0.30
    
  map_cc:
    display_name: "Map CC"
    format: "{value:.3f}"
    direction: higher_is_better
    thresholds:
      excellent: 0.80
      good: 0.75
      acceptable: 0.70
```

---

## Data Flow

### Single Cycle

```
1. INPUT
   • available_files: [data.mtz, seq.fa, refine_001.pdb]
   • history: [{cycle: 1, program: "xtriage", ...}]
   • log_text: "R-free = 0.28..."
   • user_advice: "continue refinement"

2. PERCEIVE
   • Categorize files → {mtz: [data.mtz], model: [refine_001.pdb]}
   • Detect state → "refine" phase, valid: [refine, ligandfit, molprobity]
   • Extract metrics → {r_free: 0.28, r_work: 0.24}
   • Analyze trend → "improving, no plateau"
   • Emit: STATE_DETECTED, METRICS_EXTRACTED, METRICS_TREND
   
3. PLAN
   • Check directives → no stop condition met
   • Call LLM or rules → phenix.refine
   • Validate → ✓ in valid_programs
   • Emit: PROGRAM_SELECTED (with full reasoning)
   
4. BUILD
   • Select model: refine_001.pdb (from best_files)
   • Select mtz: data.mtz (rfree_locked)
   • Build: "phenix.refine refine_001.pdb data.mtz output.prefix=refine_002"
   • Emit: FILES_SELECTED, COMMAND_BUILT

5. OUTPUT
   • program: "phenix.refine"
   • command: "phenix.refine ..."
   • events: [{type: "state_detected", ...}, ...]
   • stop: false
```

---

## Execution Modes

### Goal-Directed Mode (v114, default)

```bash
phenix.ai_agent original_files="data.mtz seq.fa"
```

- Default mode (`thinking_level=expert`)
- Multi-phase plan generated at session start from templates
- Gate evaluation after each cycle (advance/retreat/skip/stop)
- Model placement gate: detects when model fits data and
  suppresses destructive programs (phaser, autosol)
- Hypothesis testing with verification latency
- Per-cycle expert assessment and stage transition summaries
- `structure_report.html`, `structure_determination_report.txt`,
  and `session_summary.json` at completion

### Standard Mode (LLM, no planning)

```bash
phenix.ai_agent thinking_level=advanced original_files="data.mtz seq.fa"
```

- LLM makes program selection decisions
- Structural validation and expert KB
- Rules engine validates all choices
- Full reasoning captured in events
- No strategic plan or gate evaluation

### Rules-Only Mode

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz seq.fa"
```

- Deterministic program selection (first valid program)
- No LLM calls
- Faster, reproducible for testing
- Auto-discovers files from input_directory (no README
  parsing needed)

### Dry-Run Mode

```bash
phenix.ai_agent dry_run=True dry_run_scenario=xray_basic
```

- Simulated execution with predefined outcomes
- For testing workflow logic without running PHENIX programs

---

## Error Handling

### Automatic Error Recovery

Recognized error patterns trigger automatic retry with corrected parameters:
- **Ambiguous data labels**: MTZ with multiple arrays → selects anomalous or merged
  based on workflow context, injects `obs_labels` via recovery param injection
- **Ambiguous experimental phases**: Selects HL coefficients based on context
- **Loop guard** (v112.74): If a recovery strategy already exists for the file,
  re-triggering is skipped to prevent infinite retry loops

See ARCHITECTURE.md "Automatic Error Recovery" for implementation details.

### Parameter Blacklisting (`bad_inject_params`)

When a PHENIX program fails due to an injected parameter, the parameter is
blacklisted so `inject_user_params` never re-injects it. Recognized error patterns:
- **Unknown parameter**: "Unknown command line parameter definition: FOO"
- **No such parameter**: "No such parameter: FOO"
- **Boolean type mismatch** (v112.75): "True or False value expected,
  scope.path.param=value found" — blacklists the full PHIL path and all
  components ≥ 6 characters (catches `wavelength` when PHIL resolves it to
  `autosol.wavelength.added_wavelength`)

### User Request Invalid

When user requests an unavailable program:

```
============================================================
  WARNING: Requested program not available
============================================================
  You requested: phenix.refine
  Reason: Not valid in current workflow state 'xray_initial'
  Running instead: phenix.xtriage
  Available programs: phenix.xtriage
  Suggestion: This program requires different conditions.
============================================================
```

The agent:
1. Detects that user mentioned the program in their advice
2. Explains WHY it's not available
3. Suggests what will run instead
4. Always shown (QUIET verbosity level)

### Sanity Checks

Red flags that indicate problems:
- **Experiment type changed** mid-workflow → abort
- **R-free spike** (increased > 0.05) → warning
- **No model** available for refinement → abort
- **Resolution unknown** before refinement → warning

### Server Error Propagation (v112.78)

Fatal server errors (e.g., daily API usage limit) are raised as `Sorry` in
`rest/__init__.py`.  `RemoteAgent` has a dedicated `except Sorry: raise` handler
before its generic `except Exception` to ensure these propagate cleanly through
to the GUI instead of being silently swallowed as a None result.

### Cross-Platform (Windows) Considerations (v112.78)

The agent runs on macOS, Linux, and Windows. Key platform-specific handling:

| Area | Unix/macOS | Windows |
|------|-----------|---------|
| Process tree kill | `psutil` recursive walk → `SIGTERM`; fallback `os.killpg` | `taskkill /F /T /PID` |
| Abort detection | `return_code < 0` (signal) + STOPWIZARD file | STOPWIZARD file only (taskkill returns positive codes) |
| Subprocess GUI | No special handling needed | `CREATE_NO_WINDOW` creationflag prevents console flash |
| Path separators | Forward slash native | Backslash native; normalize to `/` before marker matching |
| Path quoting | `shlex.quote()` (POSIX single-quotes) | Double-quote with escaped inner `"` |
| File encoding | UTF-8 default on modern systems | Locale-dependent; explicit `encoding='utf-8'` on all `open()` |
| `os.path.relpath` | Always works | `ValueError` across drives; caught with `try/except` |

---

## Session Management

The agent provides two mechanisms for inspecting and modifying an existing session
**without running new crystallographic cycles**.

### Viewing a Session

```bash
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=basic
phenix.ai_agent log_directory=AIAgent_run1 display_and_stop=detailed
```

`basic` prints a one-line-per-cycle summary table (program, R-free, result).
`detailed` prints full reasoning and command for every cycle.
Both modes populate `self.result` identically to a normal run so GUI calls
(`get_results()`, `get_results_as_JSON()`) work without special cases.

`restart_mode=resume` is automatically set when either session management
parameter is active — no manual flag required.

### Removing Cycles

```bash
phenix.ai_agent log_directory=AIAgent_run1 remove_last_n=2
```

Removes the last N cycles from the session, clears the stale AI summary,
rebuilds `active_files.json` and `best_files` from remaining history,
and saves. Useful for pruning a failed run before re-running.

### Extending a Completed Workflow with New Advice

When a workflow has fully completed, resuming with new `project_advice`
triggers follow-up programs via the Q1 mechanism:

```bash
phenix.ai_agent \
    log_directory=AIAgent_run1 \
    restart_mode=resume \
    project_advice="also run polder on chain B residue 100"
```

1. New advice hash detected → `advice_changed=True`
2. PERCEIVE steps `complete` phase back to `validate` (adds polder, molprobity, etc.)
3. PLAN suppresses AUTO-STOP for one cycle
4. LLM acts on the new advice; after success, normal termination resumes

See [USER_DIRECTIVES.md](guides/USER_DIRECTIVES.md#extending-a-completed-workflow)
and [ARCHITECTURE.md](reference/ARCHITECTURE.md#advice-change-detection) for details.

---

## Testing

### Test Suites (42+ files, 40+ in runner)

**Standalone (no PHENIX required):**
- API Schema, Best Files Tracker, Transport, State Serialization
- Command Builder, File Categorization, File Utils
- Session Summary, Session Directives, Session Tools, Audit Fix Regressions
- Advice Preprocessing, Directive Extractor, Directive Validator, Directives Integration
- Event System, Metric Patterns, Pattern Manager
- Program Registration, Summary Display, New Programs
- Error Analyzer, Decision Flow, Phaser Multimodel
- History Analysis, Docs Tools, YAML Tools
- Thinking Defense, Strategy Memory, Log Extractor, Thinking Agent (v113)
- Structure Model, Validation History, Plan Generator, Gate Evaluator,
  Hypothesis Evaluator (v114 — 251 tests)

**PHENIX-dependent:**
- Workflow State, YAML Config, Sanity Checker
- Metrics Analyzer, Dry Run, Integration, Directives Integration

**Additional test files (not in run_all_tests.py):**
- tst_template.py (template builder), tst_utils.py (assert helpers)

### Running Tests

```bash
python3 tests/run_all_tests.py        # All tests
python3 tests/run_all_tests.py --quick  # Standalone only
python3 tests/tst_event_system.py    # Single suite
```

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v114.1 | **Model Placement Gate + Display + Evaluation Harness**: Default `thinking_level` changed from `advanced` to `expert`. **Placement gate**: detects when model fits data (model_vs_data CC > 0.3 or refine R-free < 0.50) and locks `model_is_placed` in session — suppresses phaser/autosol/predict_and_build, fast-forwards plan past MR/phasing stages, logs conflict warning when user advice contradicts. **Display**: DisplayDataModel unified data layer for Results/Progress/HTML; HTML structure report with SVG trajectory chart; "Open Structure Report" button in GUI; expert assessments now stored in session JSON; DDM scans all cycles for best metrics. **Templates**: `mr_sad` requires explicit MR-SAD intent (`wants_mr_sad`); predict_and_build in MR stage programs; polder moved to post-ligandfit phase only with `has: ligand_fit` YAML condition; SAD templates require sequence. **Files**: auto-discover from `input_directory`; HETATM ligand detection in input PDBs; ligand PDB filename hints. **GUI**: restart_mode as plain wx.Choice (survives session management reset); stage display "cycle X, up to Y". **Safety**: sanity check threshold 3→4; recent failures injected into THINK prompt; STOP not counted as cycle. **Testing**: 57 scenario tracer tests (PG1-PG5 placement gate, L1-L10 mock LLM, C1-C3 cycle counting); tutorial run analyzer for 5 modes. 31 files modified. |
| v114 | **Goal-Directed Agent** (`thinking_level=expert`): Strategic planner layer with Structure Model (running structural knowledge + strategy blacklist), Plan Generator (12 templates, strategy hash), Gate Evaluator (success hysteresis, 5 anti-oscillation safeguards), Hypothesis Engine (single-budget, verification latency, re-validation), Explanation Engine (cycle/phase/final commentary), GUI stage display (plan header, transition blocks, per-cycle stage context), `session_summary.json` output. 12 new files, 13 modified, 251 tests. Reactive agent unchanged; `advanced` (default) continues to work without planning overhead. |
| v113.10 | **Thinking Level + Validation + Expert KB**: Replaces boolean `use_thinking_agent` with `thinking_level` parameter (`none`/`basic`/`advanced`; v114 adds `expert`). Advanced mode adds structural validation (Ramachandran, clashscore, rotamers, model contents), 56-entry expert knowledge base with IDF-weighted tag matching, file metadata tracking, R-free trend display, and user-facing Expert Assessment block in event formatter. 7 new files, 25 tests. Backward compatible. |
| v113 | **Thinking Agent**: Optional expert crystallographer reasoning node (THINK) between PERCEIVE and PLAN. Second LLM call analyzes program logs with domain expertise, injects strategic guidance via user_advice enrichment. Per-program keyword extraction (xtriage, phaser, autosol, autobuild, refine), priority-ordered sections within character budget. Strategy memory persists across cycles via session_info. GUI checkbox + `[Expert]` display in progress panel. 4 new modules, 4 new test files, 103 thinking-related tests. |
| v112.78 | **GUI mode map_coeffs_mtz + daily usage Sorry + after_program fix + Windows compat**: GUI mode `_record_command_result` and `_track_output_files` used `os.getcwd()` which pointed to parent after CWD restore — now accept explicit `working_dir`; `rest/__init__.py` raises `Sorry` on `daily_usage_reached` and `RemoteAgent` re-raises it; `after_program` changed from hard stop to minimum-run guarantee; Windows: `_filter_intermediate_files` normalizes backslash paths, `Popen` uses `CREATE_NO_WINDOW`, session JSON uses explicit UTF-8 encoding |
| v112.77 | **Autobuild rebuild_in_place**: Rule D stripped `rebuild_in_place=False` because it wasn't in strategy_flags; added `rebuild_in_place`, `n_cycle_build_max`, `maps_only` to autobuild allowlist; recovery hint for sequence mismatch errors |
| v112.76 | **Catch-all injection blacklist + deterministic atom_type**: heavier-atom-wins rule swaps `atom_type`/`mad_ha_add_list` when primary has lower Z (27-element table); catch-all streak tracker blacklists injected params after 2 consecutive same-error failures (`return_injected` kwarg, `_update_inject_fail_streak`); recovery retries excluded |
| v112.75 | **Autosol/autobuild process bugs**: strategy-flag alias awareness in `inject_user_params` (wavelength→lambda dedup); `bad_inject_params` learning expanded to PHIL boolean-type errors; autosol atom_type/mad_ha_add_list same-value dedup; `_is_program_already_done` extended to non-count programs (prevents `_apply_directives` re-adding completed autosol from program_settings); improved atom_type hint in programs.yaml |
| v112.74 | **Xtriage recovery + ligand misclassification**: recovery param injection survives command builder and probe-only sanitizer; ligand-as-model misclassification guard; obs_labels error recovery loop guard |
| v112.70 | **Ligandfit file selection**: fixed refine MTZ classification regex (3 locations); word-boundary `exclude_patterns`; content-based PDB guards for model/ligand slots; protein-in-ligand-slot rejection; refinement CIF exclusion; `inject_user_params` bare-key validation; supplemental file discovery on session load and live path; fallback diagnostics (per-program missing slots); duplicate detection respects different input files |
| v112.31 | **Session management**: `display_and_stop` / `remove_last_n` populate `self.result`; `get_results()` safe before `run()`; `restart_mode` auto-set; **Q1**: resuming with new advice after workflow completion steps back from `complete` to `validate` phase, enabling follow-up programs (polder etc.) |
| v112 | **Steps table metrics**: cycle metrics as primary source; benign warning metrics extraction; ligand typing fix; case-sensitive pattern fix; autobuild_denmod detection; YAML log_parsing for 8 programs |
| v111 | **Summary output fixes**: predict_and_build R-free extraction; ligandfit output in final file list; fallback cycle status check fix |
| v110 | **Stepwise mode**: automation_path controls predict_and_build behavior; fallback program tracking; autobuild scoring equals refined; best files in summary fix |
| v40 | Fixes 12-21: Ligandfit MTZ exclusion, stop condition on failed runs, summary display, predict_and_build resolution handling; USER_REQUEST_INVALID event |
| v39 | Event system plumbing fixes for single-shot mode |
| v38 | Event system Phase 4: display integration |
| v36-37 | Event system Phases 2-3: instrumentation and transport |
| v34 | Event system Phase 1: EventLog and EventFormatter |
| v30-33 | YAML centralization, BestFilesTracker, CommandBuilder unification |

---

## Automation Modes

The agent supports two automation modes controlled by `maximum_automation`:

### Automated Mode (default)

```bash
phenix.ai_agent maximum_automation=True original_files="data.mtz sequence.fa"
```

- `predict_and_build` runs the complete workflow (prediction → MR → building)
- Fewer checkpoints, faster end-to-end processing
- Best for well-understood datasets

### Stepwise Mode

```bash
phenix.ai_agent maximum_automation=False original_files="data.mtz sequence.fa"
```

- `predict_and_build` stops after prediction only (`stop_after_predict=True`)
- User can inspect predicted model before proceeding
- Workflow continues: `process_predicted_model` → `phaser` → `refine`
- Best for troubleshooting or when intermediate inspection is needed

The `automation_path` is set in workflow_state and propagated to all decision-making components to ensure consistent behavior throughout the pipeline.
