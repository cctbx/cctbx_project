# PHENIX AI Agent — Project Status

**Date:** March 16, 2026
**Based on:** Run 17 (partial, 22 of ~47 tutorials) with all Phase 3 + Bug 5 fixes deployed

---

## Overall Metrics (Run 15 → Run 17, 21 common tutorials)

|                    | Run 15 (baseline) | Run 17 (current) |
|--------------------|-------------------|------------------|
| Failure rate       | 27%               | 24%              |
| GOOD tutorials     | 10                | 7*               |
| Total cycles       | 789               | 464              |

\* 3 of the 4 apparent "degradations" are artifacts of the partial run (missing
modes), not real regressions. See details below.

---

## Fix Verification

| Bug | Tutorial | Status | Evidence |
|-----|----------|--------|----------|
| 1 | lysozyme-refine | ⏳ Not in Run 17 yet | — |
| 3 | apoferritin_denmod_dock | ✅ **Verified** | 57% → 0% fail; CC=0.688; resolve_cryo_em succeeds |
| 4 | 7rpq_AF_reference | ✅ **Crash fixed** | 0% fail (was crashing in R16 at cycle 3). Agent doesn't advance past model_vs_data — separate workflow issue, not a crash |
| 5 | lowres_restraints | ⏳ Not in Run 17 yet | — |
| 6 | AF_7n8i | ✅ **Verified** | 26% → 3% fail; CC=0.843; half-maps correctly detected |
| 7 | 1aba-polder | ✅ **Verified** | Polder runs in 7/10 modes (was 0 in R16). ligand_pdb rescue worked. NOMET is correct — polder produces map coefficients, not R-free |

### Regression Guards (all stable)

| Tutorial | Run 15 | Run 17 | Status |
|----------|--------|--------|--------|
| 3tpj-ensemble-refine | R=0.154 | R=0.154 | ✅ Identical |
| AF_7mjs_H | CC=0.807 | CC=0.807 | ✅ Identical |
| a2u-globulin-ligand | R=0.247 | R=0.248 | ✅ Stable |
| 1J4R-ligand | R=0.204 | R=0.204 | ✅ Identical |
| AF_bromodomain_ligand | R=0.238 | R=0.259 | ✅ Stable (within normal LLM variance) |

### Apparent Degradations (all explained)

| Tutorial | R15 → R17 | Explanation |
|----------|-----------|-------------|
| 1aba-polder | GOOD → NOMET | **Fix working.** Polder runs now. NOMET because polder produces maps, not R-free. Was GOOD in R15 because solve modes ran refine. |
| AF_POMGNT2 | GOOD → POOR | **Partial run.** Only 1 mode (rules_only, R=0.553) present. In R15, LLM modes reached R=0.288. |
| AF_exoV_MRSAD | GOOD → NOMET | **Partial run.** Pipeline runs (phaser→autosol→autobuild) but autobuild fails. Same pattern as R15 — the GOOD came from modes not yet run. |
| advanced-xtriage | OK → NOMET | **Partial run.** Only xtriage modes present. In R15, solve modes ran refine → R=0.300. |

---

## What the Agent Can Do

### X-ray Crystallography
- **Data triage:** phenix.xtriage (space group, twinning, anomalous signal)
- **Molecular replacement:** phenix.phaser (multi-ensemble, twin-aware)
- **Experimental phasing:** phenix.autosol (SAD, MAD, MR-SAD)
- **Model building:** phenix.autobuild, phenix.predict_and_build
- **Refinement:** phenix.refine (NCS, reference model, restraints, twin laws)
- **Ligand fitting:** phenix.ligandfit (residue code or ligand PDB)
- **Validation:** phenix.molprobity, phenix.model_vs_data
- **Special tools:** phenix.polder (omit maps), phenix.process_predicted_model (AlphaFold)

### Cryo-EM
- **Map analysis:** phenix.mtriage (resolution, FSC)
- **Density modification:** phenix.resolve_cryo_em (half-map denoising)
- **Map operations:** phenix.map_sharpening, phenix.map_symmetry
- **Model building:** phenix.predict_and_build, phenix.map_to_model
- **Model docking:** phenix.dock_in_map
- **Refinement:** phenix.real_space_refine

### Agent Capabilities
- **5 thinking levels:** rules_only (no LLM), llm, llm_think, llm_think_advanced, llm_think_expert
- **Workflow state machine:** Tracks experiment progress, ensures correct program ordering
- **Error recovery:** Diagnosable errors stop cleanly; recoverable errors retry with fixes
- **Strategic planning** (expert mode): Multi-stage plans with gate evaluation
- **File categorization:** Automatic detection of model, data, map, sequence, ligand files
- **PHIL validation:** Whitelist + prefix + blocked params system prevents invalid parameters

---

## Architecture

```
CLIENT (ai_agent.py)                    SERVER (run_ai_agent.py)
┌─────────────────────┐                ┌──────────────────────────────────┐
│ File discovery       │                │  PERCEIVE (no LLM)               │
│ Directive extraction │  session_info  │    File categorization           │
│ Plan generation      │ ──────────►   │    Workflow state detection       │
│ Plan advancement     │                │    Metrics extraction            │
│ Gate evaluation      │                │    Sanity checks                 │
│                      │                │                                  │
│                      │  ◄──────────  │  THINK (optional LLM)            │
│ Program execution    │   command      │    Expert reasoning              │
│ Log collection       │                │    Structure model update        │
│ Output tracking      │                │                                  │
│                      │                │  PLAN (LLM or rules)             │
│                      │                │    Program selection             │
│                      │                │    File + strategy assignment    │
│                      │                │                                  │
│                      │                │  BUILD (deterministic)           │
│                      │                │    Command construction          │
│                      │                │    PHIL validation               │
│                      │                │    Path resolution               │
│                      │                │                                  │
│                      │                │  VALIDATE → FALLBACK → OUTPUT    │
└─────────────────────┘                └──────────────────────────────────┘
```

---

## Test Suite

| Suite | Tests | What it covers |
|-------|-------|----------------|
| tst_structure_model | 77 | Numeric coercion, metric tracking, problem detection |
| tst_plan_schema | 53 | Plan stages, advancement, directives |
| tst_phase3_bug5 | 85 | All Phase 3 + Bug 5 fixes (31 test functions) |
| tst_backward_compat | 26 | API compatibility, request/response format |
| tst_command_builder | 22 | Command generation, file selection |
| tst_phil_validation | 15 | PHIL whitelist, prefix matching, blocked params |
| tst_event_system | 13 | Event emission and formatting |
| **Total** | **291** | |

---

## Known Issues

### Will fix with full Run 17 data

| Issue | Impact | Status |
|-------|--------|--------|
| 7rpq doesn't advance past model_vs_data | R stays at 0.386 instead of refining to ~0.329 | Crash fixed; workflow advancement is separate issue |

### Deferred to Phase 4

| Issue | Impact | Complexity |
|-------|--------|------------|
| lysozyme-MRSAD multi-wavelength labels | 72% fail | Large — needs wavelength selection feature |
| hipip-refine [4Fe-4S] crash loop | 51% fail | Medium — refine crashes on Fe-S cluster, recovers to R=0.22 |
| a2u-globulin-rebuild autobuild crash | 63% fail | Medium — autobuild fails, retries waste cycles |

### Not bugs (correct behavior)

| Tutorial | Behavior | Why it's correct |
|----------|----------|-----------------|
| 3tpj/3tpp README 0 cycles | "requires ensemble_refinement" | Program not supported — correct rejection |
| HNP3-mr-rosetta, if5a-textal | "Cannot build" (missing data) | Test data files not present |
| apoferritin_chimerax, douse, etc. | 0 cycles | Programs not supported by agent |
| advanced-xtriage NOMET | Runs xtriage only | Correct — tutorial only asks for triage |

---

## Files Changed This Session

9 code files, 3 docs, 3 test files across 7 bug fixes:

| File | Changes |
|------|---------|
| `diagnosable_errors.yaml` | Bug 1: terminal diagnosis for bad element types |
| `programs.yaml` | Bug 3: mask_atoms removed; Bug 5A: prefix whitelist |
| `phil_validator.py` | Bug 3: blocked params; Bug 5A: prefix matching |
| `structure_model.py` | Bug 4: `_coerce_numerics` belt-and-suspenders |
| `metrics_analyzer.py` | Bug 4: `_safe_float` at 7 sites (root cause) |
| `workflow_state.py` | Bug 6: half-map Tier 1; Bug 5C: ref model categorizer; Bug 7: ligand_pdb rescue |
| `graph_nodes.py` | Bug 5D: strategy rewrites; cleanup |
| `program_registry.py` | Bug 5B: path resolution |
| `command_builder.py` | Bug 5C Tier 1: pre-existing, unblocked by Fix A |

---

## Next Steps

1. **Complete Run 17** — 25 tutorials still running. Key ones to watch: lysozyme-refine (Bug 1), lowres_restraints (Bug 5), hipip-refine (regression guard).

2. **7rpq workflow advancement** — The crash is fixed but the agent stops after model_vs_data. Need to investigate why the workflow state machine doesn't offer phenix.refine after a successful placement probe.

3. **Phase 4 sprint** — Multi-wavelength label selection (lysozyme-MRSAD) is the last major unsupported data type.
