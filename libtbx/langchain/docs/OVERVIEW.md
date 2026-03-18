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

## Quick Start

```bash
# X-ray structure determination (expert mode is default)
phenix.ai_agent original_files="data.mtz sequence.fa"

# Cryo-EM structure determination  
phenix.ai_agent original_files="map.mrc sequence.fa"

# With user guidance
phenix.ai_agent original_files="data.mtz model.pdb ligand.pdb" \
    project_advice="Solve the structure and fit the ligand"

# Point at a tutorial directory (auto-discovers files + README)
phenix.ai_agent input_directory=/path/to/tutorial/

# Stop after specific step
phenix.ai_agent original_files="data.mtz seq.fa" \
    project_advice="Stop after one refinement job"

# Rules-only mode (deterministic, no LLM, auto-discovers files)
phenix.ai_agent use_rules_only=True input_directory=/path/to/data/

# Without strategic planning (advanced reasoning only, no plan/gates)
phenix.ai_agent thinking_level=advanced \
    original_files="data.mtz sequence.fa"

# Control output verbosity
phenix.ai_agent verbosity=verbose original_files="data.mtz sequence.fa"

# Run a specific program
phenix.ai_agent original_files="model.pdb map.ccp4" \
    project_advice="Run phenix.map_correlations"

# Tutorial mode (set up tutorial in GUI, then open AI Agent)
# The agent auto-detects the tutorial and reads the README
phenix.ai_agent input_directory=/path/to/tutorial/
```

---

## Active Development

The v115 cycle addresses two categories of work.

**Infrastructure fixes** — 22 items identified from dual-run evaluation
of 21 tutorials, tracked in [`PLAN.md`](../PLAN.md) at the repo root.
20 of 22 are complete; I6 (unsupported programs) and I7 (tar.gz input)
are deferred pending workflow engine expansion.

**Bug fixes** — five issues found during v115 test suite review:
see [CHANGELOG.md v115.03](CHANGELOG.md).

Contributors should read `PLAN.md` before starting new v115 work to
avoid duplicating completed items or conflicting with deferred ones.

### Success Metrics

| Metric | Baseline (pre-v115) | Target (v115) | Measured by |
|--------|---------------------|---------------|-------------|
| Wasted program cycles | 41% (46/111 cycles) | < 25% | `analyze_tutorial_runs.py` dual-run |
| Tutorial solve rate (readme mode) | not yet measured | > 80% | `tutorial_expectations.yaml` |
| Terminal failure pivots firing | 0% (no pivot logic) | > 90% | `tst_error_classifier.py` |
| PHIL rejections causing halt | untracked | 0 | `tst_phil_validation.py` |
| Intent misclassification | untracked | < 5% on tutorial set | `tst_intent_classifier.py` |

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

See DEVELOPER_GUIDE.md §6 (Safety Checks) for the complete list.

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

### 10. Thinking Levels and the Strategic Planner

The `thinking_level` parameter controls how much intelligence the
agent applies to each cycle. It is a PHIL choice parameter with four
values: `none`, `basic`, `advanced`, and `expert` (default). Each
level is additive — higher levels include everything from lower levels
plus additional capabilities.

A separate parameter, `use_rules_only=True`, is orthogonal to
`thinking_level`. It replaces the LLM call in the PLAN node with
deterministic first-valid-program selection. Everything else — PERCEIVE,
THINK, BUILD, VALIDATE, safety checks, error recovery — runs identically
regardless of `use_rules_only`. The two axes are independent: you can
run `thinking_level=advanced use_rules_only=True` to get structural
validation and expert KB without any LLM calls in PLAN.

#### Capability Matrix

| Capability | `none` | `basic` | `advanced` | `expert` (default) |
|---|---|---|---|---|
| **Inside the LangGraph pipeline (per-cycle)** | | | | |
| PERCEIVE: file categorization, workflow state, metrics, sanity checks | ✓ | ✓ | ✓ | ✓ |
| THINK: LLM call with log analysis | — | ✓ | ✓ | ✓ |
| THINK: strategy memory (persists across cycles) | — | ✓ | ✓ | ✓ |
| THINK: structural validation (Ramachandran, clashscore, rotamers) | — | — | ✓ | ✓ |
| THINK: Structure Model (accumulated structural knowledge) | — | — | ✓ | ✓ |
| THINK: Validation History (per-cycle snapshots, trend analysis) | — | — | ✓ | ✓ |
| THINK: file metadata tracking | — | — | ✓ | ✓ |
| THINK: expert knowledge base (56 rules, IDF-weighted) | — | — | ✓ | ✓ |
| THINK: hypothesis extraction from LLM response | — | — | ✓ | ✓ |
| PLAN: LLM program selection (or rules-only fallback) | ✓ | ✓ | ✓ | ✓ |
| BUILD: deterministic command assembly | ✓ | ✓ | ✓ | ✓ |
| VALIDATE: workflow, file, and duplicate checks | ✓ | ✓ | ✓ | ✓ |
| **Outside the graph (between cycles, in ai_agent.py)** | | | | |
| Plan generation at session start (17 templates) | — | — | — | ✓ |
| Plan-to-directives merging before each cycle | — | — | — | ✓ |
| Gate evaluation after each cycle (advance/retreat/skip/stop) | — | — | — | ✓ |
| Strategy blacklisting on retreat | — | — | — | ✓ |
| Hypothesis lifecycle management (confirm/refute/abandon) | — | — | — | ✓ |
| Cycle commentary (template-based, no LLM) | — | — | — | ✓ |
| Plan enforcement (suppress premature STOP) | — | — | — | ✓ |
| Model placement gate | — | — | — | ✓ |
| **At session end** | | | | |
| Structure report (HTML with SVG trajectory) | — | — | — | ✓ |
| Session summary JSON | — | — | — | ✓ |
| Final/stopped report (template-based, no LLM) | — | — | — | ✓ |
| Failure diagnosis (LLM-generated, on terminal error) | — | — | — | ✓ |

#### Per-Level Behavior

**`none`** — Minimum viable agent. The THINK node is a complete
pass-through. The pipeline is effectively PERCEIVE → PLAN → BUILD →
VALIDATE → OUTPUT. PLAN still calls the LLM for program selection
(unless `use_rules_only=True`), so there is one LLM call per cycle.
No between-cycle operations beyond basic stop checks.

**`basic`** — Adds a second LLM call in the THINK node. The thinking
context includes log sections extracted from the last program's output
(program-specific keyword extraction), current metrics, R-free trend
across history, a brief history summary (last 5 cycles), recent
failures, and strategy memory. The LLM produces an assessment with
an action, confidence, analysis, and guidance. If guidance is produced,
it is prepended to `user_advice` as `[Expert assessment] ...` so the
PLAN node's LLM sees it when selecting a program. Strategy memory
is updated and persisted. The `should_think()` gate applies: THINK
only engages after strategic programs (analysis, model-building,
refinement, ligand), after failures, or when R-free is stalled. It
skips on the first cycle and on routine steps.

**`advanced`** — Adds four subsystems to the THINK node's context:
(A) structural validation via headless `run_validation()` on the
current best model, producing Ramachandran, clashscore, rotamer
analysis, and model contents; (B) a Structure Model that accumulates
cross-cycle knowledge from ground-truth validation results — data
characteristics, model state, R-free trajectory with annotations,
and a strategy blacklist — plus a Validation History with per-cycle
snapshots; (C) file metadata entries for validated models; and
(D) expert knowledge base queries with IDF-weighted tag matching
against 56 crystallographic rules. The LLM call receives all of this
richer context. The between-cycle loop in `ai_agent.py` is unchanged
from `none`/`basic` — no plan, no gate evaluation, no reports.

**`expert`** (default) — Everything from `advanced`, plus the entire
Strategic Planner layer activates in `ai_agent.py`. Inside the graph,
the THINK node runs identically to `advanced` (the value `expert`
is mapped to `advanced` for graph execution, since all planning
operations live outside the graph). The additions are all in the
between-cycle loop, described in Sections 10a–10h below.

#### How `expert` Maps to `advanced` in the Graph

The LangGraph pipeline only knows three thinking levels: `none`,
`basic`, and `advanced`. When the user sets `thinking_level=expert`,
`create_initial_state()` maps it to `advanced` for graph execution.
This is deliberate: the graph nodes (PERCEIVE, THINK, PLAN, BUILD,
VALIDATE) behave identically at `advanced` and `expert`. Everything
that distinguishes `expert` from `advanced` — plan generation, gate
evaluation, hypothesis lifecycle, stage summaries, reports — lives
in the outer loop in `ai_agent.py`, gated on the original
`thinking_level == "expert"` check.

The subsections below (10a–10h) describe the components of the
strategic planner in detail.

Note: The Structure Model (10a) and Validation History are updated
at both `advanced` and `expert` levels inside the THINK node. The
remaining components (10b–10h: Plan Generator, Gate Evaluator,
Hypothesis Engine, Explanation Engine, Model Placement Gate, Plan
Enforcement, Display Data Model, HTML Report) only run at `expert`.

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
| `generate_stage_summary()` | Stage transitions | No (template) |
| `generate_final_report()` | Session completion | No (template from Structure Model) |
| `generate_stopped_report()` | Early stop | No (template from Structure Model) |

The failure diagnosis (`ai_failure_diagnosis.html`) is separate from the
Explanation Engine — it is LLM-generated and produced by `failure_diagnoser.py`
when the agent detects a terminal error (see ARCHITECTURE.md).

Note: `generate_stage_summary()` no longer includes the "Next:"
line — that information is shown separately in the GUI's stage
transition block ("→ ADVANCING TO: ...").

#### 10f. Model Placement Gate (v114.1, hardened v114.2)

Prevents destructive program selection on models that already fit the
data. This is a single mechanism that addresses the class of bugs where
the LLM runs phaser/autosol on a pre-solved structure, producing a
worse result than simple refinement.

Detection (in `_detect_model_placement`, `programs/ai_agent.py`):
- `model_vs_data` with CC > 0.3 → placed
- `model_vs_data` with R-free < 0.50 from result text → placed (v114.2)
- `refine` with R-free < 0.50 → placed
- `real_space_refine` with map_cc > 0.3 → placed
- `real_space_refine` symmetry mismatch → needs_dock (v114.2)

When placement is confirmed:
1. `session.data["model_is_placed"]` — locked, survives resume
2. `valid_programs` filtered: phaser, autosol, predict_and_build removed
3. Plan fast-forward: MR/phasing stages (and model_rebuilding if
   R-free < 0.35) marked "skipped" (⊘)
4. Conflict warning if user advice mentions MR/phaser
5. Dock keywords ("dock", "docking") in advice cancel "refine" →
   "placed" assumption (v114.2)

#### 10f-b. Plan Enforcement (v114.2)

Prevents premature STOP when the strategic plan has pending stages.
Two enforcement points in `agent/graph_nodes.py`:

1. **AUTO-STOP suppression**: `plan_has_pending_stages` in
   session_info blocks metrics-based auto-stop.
2. **LLM STOP override**: When the LLM returns `program=STOP`
   but the plan has pending stages, the build node selects the
   best program from `prefer_programs`.

Contract field: `plan_has_pending_stages` (bool, default False, v4).

#### 10g. Display Data Model (`agent/display_data_model.py`)

Unified data provider for Results tab, Progress tab, and HTML report.
Eliminates duplication across three display views:

- `outcome_status`: determined / stopped / incomplete
- `outcome_message`: one-line human-readable summary
- `final_metrics`: best R-free/CC from structure model or cycle scan
- `rfree_trajectory` / `cc_trajectory`: for SVG charts
- `timeline`: compact cycle entries
- `format_cycle_compact()`: one-line-per-cycle formatting

#### 10h. Session Reports

At session completion, the agent produces up to four output files.
All reports except the failure diagnosis are template-based with no
LLM call.

**Structure Determination Report** (`structure_determination_report.txt`):
Text report generated by `generate_final_report()` or
`generate_stopped_report()` in `knowledge/explanation_prompts.py`.
Uses data from the Structure Model (ground truth, not LLM reasoning).
Contains: data characteristics (resolution, space group, twinning,
MR/phasing metrics), final model stats (R-free, R-work, clashscore,
Ramachandran, chains, ligands, waters), plan stage outcomes, and
(for stopped sessions) a strategies-that-failed section with metrics
at retreat, outstanding issues, and actionable recommendations.

**HTML Structure Report** (`structure_report.html`):
Template-based HTML report generated by `generate_html_report()` in
`knowledge/html_report_template.py`. Uses the `DisplayDataModel` as
its data source. Includes:
- Title with outcome status icon (✅ determined, ⚠ stopped, ● incomplete)
- Final metrics table (R-free, R-work, CC, clashscore, geometry)
- Inline SVG R-free/CC trajectory chart with retreat markers
- Workflow stage table with status icons (✓ complete, ✗ failed, — skipped)
- Per-cycle timeline (cycle, program, status, key metric)
- Output file paths (model, map)
The GUI provides an "Open Structure Report" button that opens this file.

**Session Summary JSON** (`session_summary.json`):
Machine-readable summary for the evaluation harness and downstream
tools. Generated by `_write_session_summary_json()` in `ai_agent.py`.
Contains: outcome (complete/stopped/incomplete/aborted/reactive),
stop_reason, stop_reason_code, final metrics, data characteristics,
stage outcomes with cycles and metrics, hypothesis outcomes (if any),
strategy blacklist, and R-free/CC metric trajectory.

**Failure Diagnosis** (`ai_failure_diagnosis.html`):
LLM-generated diagnosis produced by `_diagnose_terminal_failure()` in
`ai_agent.py` when the agent detects a terminal program error (e.g.,
crystal symmetry mismatch, missing data columns, PHIL parse failure).
The `DiagnosisDetector` identifies the error type from
`diagnosable_errors.yaml`, reads the failing program's log tail, and
sends both to the LLM (via the server's `failure_diagnosis` analysis
mode). The response is structured as three sections: "What went wrong,"
"Most likely cause," and "How to fix it." The HTML report opens
automatically in the user's browser. If the LLM is unavailable
(rules-only mode or API failure), a deterministic fallback uses the
YAML hint text instead. When a failure diagnosis is produced, the
normal structure report is skipped — the diagnosis is the user's
primary output.

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

**Future: structured results.** The `log_parsing` regex approach
works but is fragile — log format changes silently break extraction.
Newer PHENIX programs built on `ProgramTemplate` expose a
`results_as_json()` method that returns metrics as structured JSON.
As programs adopt this, the agent can read JSON results directly
instead of parsing logs, with `log_parsing` as a fallback for older
programs. See ARCHITECTURE.md "Potential improvements" for the
migration plan.

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

The agent's behavior is controlled by two independent axes:
`thinking_level` (how much intelligence per cycle) and
`use_rules_only` (LLM vs deterministic program selection).
The combinations below are the most common configurations.

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
- Two LLM calls per cycle when THINK engages (expert reasoning + program selection)
- `structure_report.html`, `structure_determination_report.txt`,
  and `session_summary.json` at completion

### Standard Mode (LLM, no planning)

```bash
phenix.ai_agent thinking_level=advanced original_files="data.mtz seq.fa"
```

- `thinking_level=advanced`: THINK node runs with full context
  (structural validation, expert KB, Structure Model) but no
  strategic planning layer
- LLM makes program selection decisions in PLAN
- Two LLM calls per cycle when THINK engages
- No plan generation, gate evaluation, or between-cycle reports
- Useful when the planning overhead is unnecessary or when
  debugging the reactive engine in isolation

### Minimal LLM Mode

```bash
phenix.ai_agent thinking_level=none original_files="data.mtz seq.fa"
```

- THINK node is a pass-through (no expert reasoning)
- One LLM call per cycle (PLAN only)
- Lowest LLM cost; suitable for simple workflows or constrained
  API budgets

### Rules-Only Mode

```bash
phenix.ai_agent use_rules_only=True original_files="data.mtz seq.fa"
```

- Deterministic program selection (first valid program from
  workflow engine) — no LLM call in PLAN
- `use_rules_only` only affects PLAN; `thinking_level` still
  controls the THINK node independently. At the default
  `thinking_level=expert`, THINK will attempt its LLM call.
  To suppress all LLM calls, also set `thinking_level=none`.
- Faster, fully reproducible for testing (with `thinking_level=none`)
- Auto-discovers files from input_directory

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

### Server Error Propagation (v112.78, extended v114.2)

Fatal server errors (e.g., daily API usage limit) are raised as `Sorry` in
`rest/__init__.py`.  `RemoteAgent` has a dedicated `except Sorry: raise` handler
before its generic `except Exception` to ensure these propagate cleanly through
to the GUI instead of being silently swallowed as a None result.

**v114.2 extension:** `RemoteAgent._send_request` also re-raises Sorry when
`server_result.success` is False and the `server_message` contains fatal
keywords ("quota", "API key", "rate limit", "authentication", "permission
denied", "cannot continue"). The same check applies to the parsed JSON
response's `error` field. Previously, these errors were logged and returned as
None, causing the agent to silently end with "No command generated".

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
| JSON transport | No issues | v115.05: removed `json_str.replace('\\t', ' ')` from `transport.py` and `api_client.py` — it corrupted paths like `C:\tutorials\test.mtz` (the `\t` matched as tab escape). Tab handling is done correctly by `text_as_simple_string()`. |

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

See [USER_GUIDE.md §11](USER_GUIDE.md#11-directives-reference)
and [ARCHITECTURE.md](ARCHITECTURE.md) for details.

---

## Testing

### Test Suites (55+ files, 50+ in runner)

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

**Systematic Testing Framework (v115.08 — 10 phases, ~14s):**

Bottom-up testing of system boundaries where bugs hide. Exercises real
production code paths end-to-end. Designed after code reviews found 5
critical bugs that existing unit tests missed.

| Suite | Phase | Tests | What It Tests |
|-------|-------|-------|---------------|
| S0 | Static Audit | 5 | Parse check, bare except scan, import fallbacks |
| S1 | Contract Gaps | 128 | AST coverage map — 4 modules, boundary gap detection |
| S2 | Path Consistency | 10 | YAML vs hardcoded categorization diff (whitelisted) |
| S3 | Session Round-Trip | 28 | JSON symmetry + AgentSession save/load/pipeline |
| S4 | History Flags | 8 | Flag writer/reader consistency, dead flag audit |
| S5 | Category-Consumer | 14 | input_priorities + fallback_categories alignment |
| S6 | Routing Simulation | 32 | 3-cycle routing through detect_step + get_valid_programs |
| S7 | Command Building | 15 | CommandBuilder.build() with real file combinations |
| S8 | Error Classification | 7 | 3 classifiers × 30+ patterns, overlap + severity |
| S9 | LLM Perturbation | 17 | Filename/program/parameter/truncation/empty resilience |

Phases S6–S8 are skipped in `--quick` mode. All phases produce machine-readable
findings in `findings/`. See `docs/PHASE_REVIEW_REPORT.md` for review details.

**Additional test files (not in run_all_tests.py):**
- tst_template.py (template builder), tst_utils.py (assert helpers)

### Running Tests

```bash
python3 tests/run_all_tests.py        # All tests
python3 tests/run_all_tests.py --quick  # Standalone only (skips S6-S8)
python3 tests/tst_event_system.py    # Single suite
python3 tests/tst_phase7_routing_simulation.py  # Single phase
```

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v115.09 | **Tutorial Routing Fixes**: (1) Cryo-EM `past_analysis` gate: added `map_sharpening_done`, `map_symmetry_done`, `has_optimized_full_map` — unblocks bgal_denmod, apoferritin_denmod, ion_channel_denmod after map_sharpening; regex broadened to match actual output filenames. (2) `.sca-only` data detection in `perceive()` — deferred, needs reactive approach. (3) Validation-only routing: `wants_validation_only` directive extracted via LLM prompt + rules-based fallback + post-LLM overlay → validation shortcut in `_detect_xray_step` → `validate_existing` plan template; PDB scan limit 500→2000 in `_is_valid_file` (3dnd.pdb has 546 header lines); `has_phased_data_mtz` added. (4) MR-SAD routing: `force_mr` flag when `use_mr_sad` + model not categorized as search_model → phaser offered; MR-SAD guard updated. Critical deployment fix: rules-based intent patterns moved to shared `_apply_workflow_intent_fallback()` called as post-LLM overlay. 7 files modified. |
| v115.08 | **Phased File Detection + Systematic Testing Framework**: 4 critical fixes for phased file detection (content-based iotbx+ASCII heuristic replacing filename markers; shared post-processing; category exclusivity; conditional data_mtz removal). 1 additional bug fix (B1: last_program missing from build_context). `[GATE]` diagnostic logging for routing debugging. 10-phase systematic testing framework (S0–S9): static audit, contract gap coverage map, YAML/hardcoded path consistency, AgentSession round-trip, history flag consistency, category-consumer alignment, 32-tutorial routing simulation, command building, error classification, LLM perturbation resilience. 12 new test files, 42 unit tests + ~260 phase checks. Framework reviewed across multiple rounds — 54 issues found and fixed in test scripts (5 can't-fail gates, 3 tautological assertions, 5 wrong best_files keys, etc.). |
| v115.07 | **Run 15b Bug Fixes — Phase 3**: 4 bugs from 371-run analysis. Numeric coercion (`_safe_float` at 6 sites). Two-tier half-map detection. PHIL blocked params for resolve_cryo_em. Terminal diagnosis for unknown chemical elements. Reference model restraint support (hierarchical prefix whitelist, path resolution, strategy rewrites). |
| v115.05 | **Guard Fixes + Polder + Templates**: (1) `_is_at_target` hopeless R-free (> 0.50) now requires `autobuild_done` — prevents premature stop on incomplete models (p9-SAD fix). (2) `_is_at_target` clashscore path requires `refine_count >= 1` (X-ray only; cryo-EM path is theoretical-only gap). (3) Bug F `obtain_model` routing requires `autobuild_done` — gives autobuild a chance before concluding MR is wrong. (4) Negligible-anomalous guard: removes autosol from `valid_programs` when measurability < 0.05 and `has_anomalous=False`. (5) `wants_polder` context flag + polder override fires despite "solve" keyword in advice. (6) `refine_placed_polder` template added (17 templates total). (7) Early rebuild gate: `r_free > 0.50 after 1 cycles → try_rebuilding`; `gate_evaluator` now advances to `model_rebuilding` stage. (8) Phaser copies injection reads `log_analysis["n_copies"]` same-cycle (not 1-cycle delay). (9) `_anti_ligand_patterns` excludes `no_ligand` from ligand file classification; orphaned PDB files promoted to `model` (fixes user-supplied input PDBs like `1aba.pdb` not being recognized). (10) Unregistered `explicit_program` downgrades to warning (no Sorry). (11) `_preprocessing_programs` / `_needs_plan_programs` ensure polder etc. get full plans. (12) Failed programs skip output file tracking. (13) Metrics-based report selection overrides INCOMPLETE status when R-free/CC targets are met. 15 fix-verification tests (`tst_fix_verification.py`). (14) Ligand PDB plan selection: `_build_context()` now checks ligand name hints before setting `has_search_model`, and removes the `len(pdb_files) >= 2` guard — fixes AF_bromodomain_ligand tutorial selecting `mr_refine` (no ligandfit) instead of `predict_refine_ligand`. (15) Windows transport fix: removed `json_str.replace('\\t', ' ')` from `transport.py` and `api_client.py` — the replace corrupted Windows file paths containing `\t` sequences (e.g. `C:\tutorials\test.mtz`), causing "Failed to parse request JSON" on Windows clients. |
| v115 | **Infrastructure Audit + Failure Recovery**: Dual-run evaluation framework across 21 tutorials (41% baseline cycle waste identified). Intent classifier (4-way: solve/solve_constrained/task/tutorial). Tiered error recovery (`error_classifier.py`). PHIL strategy validation. Thinking agent context forwarding. Session bug fixes (resolution contract, intent low-confidence guard, pipeline classification, skip_validation stop trigger). |
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

---

## Dependencies

The AI agent requires the following Python packages beyond the standard PHENIX
installation. Install via `phenix.python -m pip install <package>` or via the
`install_ai_tools.csh` script.

| Package | Purpose | Required on |
|---------|---------|-------------|
| `langchain-core`, `langchain-community` | LLM orchestration core and community integrations | Server and local |
| `langchain-google-genai` | Google Gemini LLM provider | Server and local |
| `langchain-openai` | OpenAI LLM provider | Server and local |
| `langchain-chroma` | Chroma vector store for document retrieval | Server and local |
| `flashrank` | Local cross-encoder reranking (no API key needed) | Server, or local if `run_on_server=False` |
| `markdown-it-py` | HTML rendering of analysis output | Server and local |

**Note:** `flashrank` downloads its model (~34MB) automatically on first use.
No Cohere API key is required — reranking runs entirely locally.

**Note:** The `langchain-classic` package is **not required**. The agent implements
document chain and compression retriever functionality directly using `langchain-core`
base classes, avoiding the deprecated `langchain.chains` and `langchain.retrievers`
modules.

---


---

## Command-Line Tools

### Configuration Management

```bash
# List all YAML configuration files
python3 agent/yaml_tools.py list

# Validate all YAML files for syntax and structural errors
python3 agent/yaml_tools.py validate

# Validate a specific file
python3 agent/yaml_tools.py validate programs.yaml

# Display formatted contents of a YAML file
python3 agent/yaml_tools.py display programs

# Compare two YAML files or directories
python3 agent/yaml_tools.py compare programs.yaml programs_backup.yaml

# Show overview of all configuration
python3 agent/yaml_tools.py summary

# Show all defined terms in the configuration system
python3 agent/yaml_tools.py terms
python3 agent/yaml_tools.py terms --detail full  # With cross-references
```

### Session Management

```bash
# Show current session summary
python3 agent/session_tools.py --show

# Show detailed session info (files, metrics, reasoning for each cycle)
python3 agent/session_tools.py --show --detailed

# Remove last N cycles from session
python3 agent/session_tools.py --remove-last 2

# Reset entire session
python3 agent/session_tools.py --reset

# Dry-run (show what would be done without saving)
python3 agent/session_tools.py --remove-last 3 --dry-run

# Use a specific session directory
python3 agent/session_tools.py --dir /path/to/session --show
```

### Documentation Generation

```bash
# Generate safety checks documentation
python3 agent/generate_safety_docs.py  # output reviewed in DEVELOPER_GUIDE.md §6
```

### RAG Documentation Database

The agent uses a Retrieval-Augmented Generation pipeline to ground LLM
responses in PHENIX documentation. Documents are stored in a Chroma vector
database and retrieved with FlashRank cross-encoder reranking (top 8 of 20
candidates). See OVERVIEW.md §9 for architecture details.

```bash
# Build vector database from PHENIX documentation
python3 run_build_db.py

# Inspect database contents
python3 run_inspect_db.py

# Query the documentation
python3 run_query_docs.py "How do I set up SAD phasing in phenix.autosol?"
```

### Program Validation

```bash
# Validate a specific program's configuration completeness
python3 agent/program_validator.py phenix.polder

# Validate all programs
python3 agent/program_validator.py --all

# List all configured programs
python3 agent/program_validator.py --list
```

### Pattern Management

```bash
# Validate all metric extraction patterns and run tests
python3 agent/pattern_manager.py
```

### Directive Validation

```bash
# Run directive validator self-test (checks program availability detection)
python3 agent/directive_validator.py
```

### Testing

```bash
# Run all standalone tests (no PHENIX required)
python3 tests/run_all_tests.py --quick

# Run all tests (including PHENIX-dependent)
python3 tests/run_all_tests.py

# Verbose output
python3 tests/run_all_tests.py --verbose

# Individual test file
python3 tests/tst_file_utils.py

# Tests matching a pattern
python3 tests/run_all_tests.py --pattern "directive"
```

---


---

## Directory Structure

```
improved_agent_v2/
├── agent/                      # Core agent logic
│   ├── graph.py                # LangGraph state machine definition
│   ├── graph_state.py          # Agent state type definitions
│   ├── graph_nodes.py          # LangGraph node implementations
│   ├── planner.py              # Agent planning and next-move generation
│   ├── workflow_engine.py      # YAML workflow interpreter
│   ├── workflow_state.py       # State detection, PDB content guards
│   ├── command_builder.py      # Unified command generation (with content guards)
│   ├── command_postprocessor.py # Server-safe command transforms (sanitize, inject)
│   ├── template_builder.py     # YAML-driven command templates
│   ├── file_utils.py           # Shared file classification (MTZ type, exclude patterns)
│   ├── best_files_tracker.py   # Track best files per type
│   ├── placement_checker.py    # Unit cell comparison for model placement
│   ├── error_analyzer.py       # Automatic error recovery
│   ├── advice_preprocessor.py  # README discovery, advice processing
│   ├── directive_extractor.py  # Parse user directives from advice
│   ├── directive_validator.py  # Pre-validate user requests
│   ├── sanity_checker.py       # Red flag detection
│   ├── rules_selector.py       # Rules-only program selection
│   ├── config_loader.py        # Load decision_config.json thresholds
│   ├── memory.py               # Persistent learned syntax tips
│   ├── metrics_analyzer.py     # Metric trends and convergence
│   ├── metric_evaluator.py     # Metric quality evaluation
│   ├── session.py              # Persistent session tracking (AgentSession, duplicate detection, supplemental file discovery)
│   ├── api_client.py           # V2 API request/response building
│   ├── phenix_utils.py         # REST encoding, standalone PHENIX utilities
│   ├── transport.py            # Sanitization and encoding
│   ├── rate_limit_handler.py   # LLM rate limiting with backoff
│   ├── program_registry.py     # YAML program registry
│   ├── pattern_manager.py      # Regex pattern management
│   ├── event_log.py            # Structured event logging
│   ├── event_formatter.py      # Output formatting (verbosity levels)
│   ├── dry_run_manager.py      # Testing: dry-run workflow simulation
│   ├── utils.py                # General utility functions
│   ├── command_templates.json  # Program command templates and file slots
│   ├── decision_config.json    # Tiered decision rules and thresholds
│   ├── parameter_fixes.json    # Wrong→correct parameter name mappings
│   ├── yaml_tools.py           # CLI: YAML validation and inspection
│   ├── session_tools.py        # CLI: Session management
│   ├── docs_tools.py           # CLI: Documentation generation
│   ├── generate_safety_docs.py # CLI: Safety checks documentation
│   ├── generate_logic_doc.py   # CLI: Decision logic documentation
│   └── program_validator.py    # CLI: Program config validation
├── knowledge/                  # Configuration & domain knowledge
│   ├── programs.yaml           # Program definitions
│   ├── workflows.yaml          # Workflow state machines
│   ├── metrics.yaml            # Quality thresholds
│   ├── file_categories.yaml    # File categorization
│   ├── patterns.yaml           # Regex patterns
│   ├── recoverable_errors.yaml # Error recovery patterns
│   ├── transport.yaml          # Sanitization rules
│   ├── api_schema.py           # V2 API schema definitions
│   ├── yaml_loader.py          # Configuration loading
│   ├── metric_patterns.py      # YAML-driven metric extraction
│   ├── phenix_programs.py      # Program discovery and introspection
│   ├── program_registration.py # Program detection from logs
│   ├── prompts.py              # LLM prompts for planning and commands
│   ├── prompts_hybrid.py       # Hybrid planning prompts (rules + LLM)
│   └── summary_display.py      # Quality table formatting
├── phenix_ai/                  # Runtime entry points
│   ├── local_agent.py          # Local execution (same process)
│   ├── remote_agent.py         # Server execution (REST API)
│   ├── run_ai_agent.py         # Decision engine (graph execution)
│   ├── run_ai_analysis.py      # Log analysis (standalone)
│   ├── log_parsers.py          # Log metric extraction
│   └── utilities.py            # Shared runtime utilities
├── programs/                   # PHENIX program integration
│   ├── ai_agent.py             # Main PHENIX entry point
│   └── ai_analysis.py          # Log analysis entry point
├── analysis/                   # Log analysis and post-run analysis
│   ├── analyzer.py             # RAG-based log analysis
│   ├── log_info.py             # High-level log info extraction
│   ├── state_extractor.py      # Project state extraction from logs
│   ├── summarizer.py           # Map-reduce log summarization
│   └── agent_session_analyzer.py # Session performance analysis
├── core/                       # LLM integration
│   ├── llm.py                  # Provider abstraction (Google, OpenAI, etc.)
│   └── types.py                # Core data types (AgentPlan, etc.)
├── commands/                   # Command building framework
│   └── base.py                 # Abstract base class for command builders
├── strategies/                 # Planning strategy framework
│   └── base.py                 # Abstract base class for planning strategies
├── validation/                 # Command validation framework
│   ├── base.py                 # Abstract base class for validators
│   ├── core_validator.py       # Syntax validation and LLM-based fixing
│   ├── mtz_utils.py            # MTZ R-free flag checking
│   ├── phenix_refine.py        # phenix.refine-specific validation
│   └── registry.py             # Validator auto-discovery and lookup
├── rag/                        # RAG (Retrieval-Augmented Generation)
│   ├── document_loader.py      # Document loading and chunking
│   ├── retriever.py            # Vector DB retrieval with FlashRank reranking
│   └── vector_store.py         # Chroma vector store creation
├── utils/                      # General utilities
│   ├── query.py                # Documentation query interface
│   ├── run_utils.py            # Log parsing and HTML output
│   └── text_processing.py      # Text block extraction helpers
├── phenix_knowledge.py         # Phenix program allow-list and syntax hints
├── phenix_learned_info/        # Persistent learned knowledge
│   └── phenix_learned_memory.json # Learned syntax tips per program
├── run_build_db.py             # CLI: Build RAG documentation database
├── run_inspect_db.py           # CLI: Inspect RAG database contents
├── run_query_docs.py           # CLI: Query PHENIX documentation
├── tst_langchain_tools.py      # Unit tests for external modules
├── tests/                      # Test suites (~745+ tests)
│   ├── run_all_tests.py        # Test runner (32 registered suites)
│   ├── tst_utils.py           # Assert helpers (cctbx-style)
│   └── scenarios/              # Dry-run test scenarios
└── docs/                       # Documentation
    ├── README.md               # This file
    ├── OVERVIEW.md             # Technical overview
    ├── DEVELOPER_GUIDE.md      # How-to guides (adding programs, testing, safety)
    ├── USER_GUIDE.md           # User guide + directives reference
    ├── reference/              # API and logic reference
    │   ├── ARCHITECTURE.md     # Component deep-dive
    │   ├── API_REFERENCE.md    # V2 JSON API spec
    │   └── VALIDATION.md       # Validation system
    ├── project/                # Design history and changelog
    │   ├── CHANGELOG.md        # Version history
    │   ├── THOUGHT_EXPERIMENT.md # Example workflow traces
    │   └── TRANSPARENCY_LOGGING.md # Event system design
```

---

