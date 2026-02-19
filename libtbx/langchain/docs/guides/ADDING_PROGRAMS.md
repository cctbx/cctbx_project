# Adding New Programs to PHENIX AI Agent

This guide documents the process for adding a new PHENIX program to the AI Agent. The configuration has been centralized so that **adding a new program only requires editing 2-3 files** for most cases.

## Quick Overview

| File | What to Add | Required? |
|------|-------------|-----------|
| `knowledge/programs.yaml` | Complete program definition with metrics | ✅ Yes |
| `knowledge/workflows.yaml` | Add to appropriate workflow phase(s) | ✅ Yes |
| `knowledge/file_categories.yaml` | New output file types (if any) | If new file types |
| `knowledge/metrics.yaml` | Custom summary display (quality table, step metrics) | Usually not needed |
| `agent/directive_extractor.py` | Tutorial/procedure patterns | Optional |
| `phenix_ai/log_parsers.py` | Complex metric extraction (tables, multi-line) | Only if YAML can't handle it |

**For most programs**, you only need to edit **programs.yaml** and **workflows.yaml**.

The system automatically handles:
- ✅ Metric extraction from logs (via `log_parsing` in programs.yaml)
- ✅ `<program>_done` tracking flags (via `done_tracking` in programs.yaml)
- ✅ Session summary display (default formatting works for most metrics)
- ✅ Workflow context building (automatic from programs.yaml)
- ✅ Prerequisite auto-execution (e.g., mtriage for resolution)
- ✅ skip_programs → workflow phase advancement (done flags auto-set)
- ✅ LLM hallucination sanitization (resolution stripped if unverified)
- ✅ .eff file generation with short-form PHIL resolution

---

## Minimal Example (2 files)

For a simple analysis program that runs once and extracts one metric:

**programs.yaml:**
```yaml
phenix.new_analysis:
  description: "Analyze something in the data"
  category: analysis
  experiment_types: [xray, cryoem]
  done_tracking:
    flag: "new_analysis_done"
    strategy: "run_once"
  
  inputs:
    required:
      data:
        extensions: [.mtz, .mrc]
        flag: ""
  
  command: "phenix.new_analysis {data}"
  
  log_parsing:
    quality_score:
      pattern: 'Quality Score[:\s]+([0-9.]+)'
      type: float
```

**workflows.yaml:**
```yaml
xray:
  phases:
    analyze:
      programs:
        - program: phenix.new_analysis
          conditions:
            - not_done: new_analysis
```

That's it! The agent will now:
- Run `phenix.new_analysis` during the analyze phase
- Extract `quality_score` from logs
- Track that it's been run (won't run again)
- Display the metric in session summaries

---

## Step 1: Define the Program in `programs.yaml`

Location: `knowledge/programs.yaml`

This is the **single source of truth** for program configuration.

```yaml
phenix.new_program:
  description: "Brief description of what the program does"
  category: analysis|phasing|building|refinement|validation
  experiment_types: [xray, cryoem]  # or just one
  
  # Done tracking: controls workflow done flags and strategy behavior
  done_tracking:
    flag: "new_program_done"   # Flag name set in workflow context
    strategy: "run_once"       # "set_flag" (default) | "run_once" | "count"
    # count_field: "new_count" # Required when strategy: "count"
    history_detection:         # Auto-set done flag from history
      markers: ["new_program"] # Strings to match (substring, OR logic)
      # exclude_markers: ["other_thing"]  # Reject match if present (checked FIRST)
      # alt_markers: ["alt_name"]  # Alt match (with alt_requires, AND logic)
      # alt_requires: ["special_flag"]

  # GUI app_id: maps to wxGUI2 window for project History (optional)
  gui_app_id: "NewProgram"           # From wxGUI2/Programs/__init__.params
  # gui_app_id_cryoem: "NewProgramCryoEM"  # If cryo-EM uses a different window

  # Stop directive patterns: regex for user stop conditions (optional)
  # stop_directive_patterns:
  #   - 'new_program'

  inputs:
    required:
      input_name:
        extensions: [.ext1, .ext2]
        flag: "input_flag="  # or "" for positional
    optional:
      optional_input:
        extensions: [.ext]
        flag: "optional_flag="
        # Set auto_fill: false to prevent the command builder from
        # automatically filling this slot with a matching file.
        # Use when a slot should only be filled if the LLM or user
        # explicitly assigns a file to it (e.g., a second map for
        # map-map correlation — without this, both map slots would
        # get the same file).
        auto_fill: false  # Default is true

  # Control which file categories are preferred for each input
  # IMPORTANT: Without this, the agent may select wrong files!
  input_priorities:
    input_name:
      categories: [preferred_category, fallback_category]
      prefer_subcategories: [refined, with_ligand]  # Within category, prefer these
      exclude_categories: [excluded_category]

  outputs:
    files:
      - pattern: "*.output_ext"
        type: output_type  # Used for file categorization
    metrics:
      - metric_name  # Metrics this program produces

  command: "phenix.new_program {input_name}"

  # Optional parameters that can be set via strategy
  strategy_flags:
    parameter_name:
      flag: "parameter={value}"
      type: float|int|string|bool
      hint: "Description of parameter"

  # IMPORTANT: Define metric extraction patterns here
  # These are used automatically by log_parsers.py and session.py
  log_parsing:
    metric_name:
      pattern: 'Regex pattern with (capture group)'
      type: float|int|string|bool
      display_name: "Human Readable Name"    # For summary display
      summary_format: "{value:.2f}"          # Format for display
      # Optional: Handle "no result" cases
      no_match_pattern: 'No.*found'
      no_match_value: "None"

  hints:
    - "Hint for when/how to use this program"

  # Keywords that trigger this program from user advice
  user_advice_keywords:
    - "keyword1"
    - "keyword2"
```

### Key Fields for Automatic Handling

| Field | Effect |
|-------|--------|
| `done_tracking.flag` | Names the workflow done flag (e.g., `xtriage_done`) |
| `done_tracking.strategy` | `"set_flag"` (default), `"run_once"` (filter after first run), or `"count"` (increment counter) |
| `done_tracking.count_field` | Counter name for `strategy: "count"` (e.g., `"refine_count"`) — must be in ALLOWED_COUNT_FIELDS |
| `done_tracking.history_detection` | YAML-driven done-flag setting (markers, exclude_markers, alt_markers, success_flag) |
| `gui_app_id` | wxGUI2 app_id for project History (fallback when PHIL unavailable) |
| `gui_app_id_cryoem` | Cryo-EM variant app_id if the program uses a different GUI window |
| `stop_directive_patterns` | Regex patterns for matching user stop-condition directives |
| `log_parsing:` | Patterns used by log_parsers.py AND session.py for metric extraction |
| `log_parsing.*.display_name` | Used in session summary display |
| `log_parsing.*.summary_format` | Format string for summary display |

---

## Step 2: Add to Workflow in `workflows.yaml`

Location: `knowledge/workflows.yaml`

Find the appropriate phase and add your program:

```yaml
xray:  # or cryoem
  phases:
    phase_name:
      programs:
        # Simple entry (always available)
        - phenix.new_program
        
        # Or with conditions
        - program: phenix.new_program
          conditions:
            - has: required_file_type  # e.g., has: sequence
            - not_done: new_program    # Don't run if already done
          hint: "When to use this program"
```

**Common conditions:**
- `has: <file_type>` - Requires file (sequence, model, full_map, half_map, ligand_file, anomalous, etc.)
- `has_any: [type1, type2]` - Requires at least one of these files
- `not_done: <program>` - Program hasn't been run yet (flag name from `done_tracking.flag`)
- `not_has: <file_type>` - File type is NOT present
- `r_free: "< 0.35"` - R-free below threshold (also supports `"> autobuild_threshold"`, `"< target_r_free"`)
- `refine_count: "> 0"` - At least one successful refinement completed (also `rsr_count`)

**Available done flags:**

All done flags are defined in programs.yaml under `done_tracking.flag`. Common flags:

| Flag | Set When | Used By |
|------|----------|---------|
| `predict_done` | predict_and_build runs (any mode) | cryo-EM workflow |
| `predict_full_done` | predict_and_build completes full workflow | X-ray workflow |
| `process_predicted_done` | process_predicted_model succeeds | cryo-EM workflow |
| `phaser_done` | phaser completes MR | X-ray workflow |
| `phaser_count` | Count of successful phaser runs | Internal tracking |
| `dock_done` | dock_in_map completes | cryo-EM workflow |
| `autobuild_done` | autobuild succeeds | X-ray workflow |
| `autobuild_denmod_done` | autobuild maps_only completes | X-ray workflow (before ligandfit) |
| `autosol_done` | autosol succeeds | X-ray workflow |
| `ligandfit_done` | ligandfit completes | Both workflows |
| `pdbtools_done` | pdbtools completes | X-ray workflow (combine_ligand) |
| `refine_done` | refine completes (any) | Internal tracking |
| `refine_count` | Count of successful refinements | Workflow conditions |
| `rsr_done` | real_space_refine completes | Internal tracking |
| `rsr_count` | Count of successful RSR runs | Workflow conditions |
| `resolve_cryo_em_done` | resolve_cryo_em completes | cryo-EM workflow |
| `map_sharpening_done` | map_sharpening completes | cryo-EM workflow |
| `map_to_model_done` | map_to_model completes | cryo-EM workflow |
| `validation_done` | molprobity/validation completes | Internal tracking |

**Note:** `refine_count` and `rsr_count` only count SUCCESSFUL runs (not failed ones).

---

## Step 2.5: Add Done Tracking (if needed)

Location: `knowledge/programs.yaml` (same file as Step 1)

**All done flags are defined in the `done_tracking` block** in programs.yaml.
No Python file edits are needed — `get_program_done_flag_map()` reads directly
from YAML.

```yaml
phenix.new_program:
  # ... other fields ...
  done_tracking:
    flag: "new_program_done"    # Flag name in workflow context
    strategy: "run_once"        # "set_flag" | "run_once" | "count"
    history_detection:          # YAML-driven done-flag detection
      markers: ["new_program"]  # Match any of these in history (OR logic)
      # exclude_markers: [...]  # Reject match if present (checked FIRST)
      # alt_markers: ["other_name"]  # Alternative markers (with alt_requires)
      # alt_requires: ["some_flag"]  # ALL must match alongside alt_markers
      # success_flag: "new_program_success"  # Extra flag set on success
```

**Strategies:**
- `"set_flag"` (default): Set done flag on success. Most programs use this.
- `"run_once"`: Set done flag + filter program from valid list after first run.
  Use for analysis programs that should never re-run (xtriage, mtriage).
- `"count"`: Set done flag + increment a counter. Requires `count_field`.
  Use for programs that can run multiple times (refine, rsr).

**`history_detection`** replaces if/elif blocks in `_analyze_history()`.
Use it when your program only needs a boolean done flag set on success.
Do NOT use it if your program needs cascading flags (like
predict_and_build → refine_done), or exclusion logic (like refine excluding
real_space).

**`stop_directive_patterns`** (optional, top-level): If users may reference your
program in stop conditions (e.g., "stop after running new_program"), add regex
patterns:

```yaml
phenix.new_program:
  stop_directive_patterns:
    - 'new_program'
    - 'new\\s+program'
```

Patterns are sorted by length at load time (longest first) so more specific
patterns match before shorter ones.

**When to add `done_tracking`:**
- Your program gates a workflow phase (e.g., xtriage gates "analyze")
- Users might skip it via `skip_programs` (the flag must be set for phase advancement)
- You want `strategy: "run_once"` behavior (program filtered from valid list after completion)
- You need to track how many times a program ran (`strategy: "count"`)

**When to skip:** Validation programs like `phenix.map_correlations` that don't
gate any phase and can run freely don't need `done_tracking`.

**Also ensure** the done flag is set in `build_context()` in `agent/workflow_engine.py`
(the flags are populated from `history_info`):

```python
context = {
    # ... existing flags ...
    "new_program_done": history_info.get("new_program_done", False),
}
```

And in `_analyze_history()` in `agent/workflow_state.py`, detect when the program
has been run and set the flag in `history_info`.

---

## Step 3: Add Output File Types (if needed)

Location: `knowledge/file_categories.yaml`

If your program produces a **new type of output file**, define it here:

```yaml
new_file_type:
  description: "Description of what this file contains"
  extensions: [.ext]
  patterns:
    - "*.ext"
    - "specific_name.ext"
  notes: "Which programs produce/consume this file type"
```

**Skip this step if** your program only produces standard file types (`.pdb`, `.mtz`, `.mrc`, etc.) that are already defined.

**Example** - map_symmetry produces `.ncs_spec` files:
```yaml
ncs_spec:
  description: "NCS symmetry operators file"
  extensions: [.ncs_spec]
  patterns:
    - "*.ncs_spec"
  notes: "From map_symmetry, used by resolve_cryo_em, predict_and_build"
```

---

## Step 4: Add Summary Display (usually not needed)

Location: `knowledge/metrics.yaml`

**Most programs don't need this step.** The default formatting works for standard metrics. Only add custom display if you need:
- Special formatting in the Final Quality table
- Custom step metric display in the workflow steps table

Location: `knowledge/metrics.yaml`

If your program produces metrics that should appear in the session summary, add them to the `summary_display` section:

### For the Final Quality Table:

```yaml
summary_display:
  quality_table:
    # Add your metric
    - metrics: [new_metric]
      label: "New Metric"
      format: "{new_metric:.2f}"
      experiment_types: [xray]  # Optional filter
      assessment: true          # Show Good/Acceptable/Needs Improvement
```

### For Per-Step Display:

```yaml
summary_display:
  step_metrics:
    phenix.new_program:
      format: "Result: {new_metric:.2f}"
      metrics: [new_metric]
      fallback_format: "Completed"  # When metrics unavailable
```

---

## Step 5: Add Tutorial Detection (Optional)

Location: `agent/directive_extractor.py`

If users might ask to run this specific program as a tutorial or procedure:

```python
tutorial_patterns = [
    # ... existing patterns ...
    (r'(?:run|analyze).*new_program', 'phenix.new_program'),
    (r'keyword.*pattern', 'phenix.new_program'),
]
```

---

## Step 6: Add Tests

Verify your program works correctly:

```bash
# Validate program configuration
python agent/program_validator.py phenix.new_program

# Run all tests — including the integration suite
python tests/run_all_tests.py --quick
```

The validator checks:
- ✅ Program defined in programs.yaml
- ✅ Added to workflows.yaml
- ✅ Log parsing patterns defined
- ✅ Summary display configured (if applicable)

### Run the full suite, not just your new tests

`workflow_state.py` is the integration boundary. When adding new behavior there
(zombie checks, file validation, diagnostic surfacing), it should be treated as a
refactor of the public API — requiring a full suite run, not just the feature's own
tests. Any PR touching `detect_workflow_state()` must show passing results for
`tst_workflow_state.py` and `tst_integration.py` before merging.

The reason this matters: the pre-existing tests in those files use hypothetical
filename strings (e.g. `"data.mtz"`, `"refined_model.pdb"`) without creating real
files on disk, and use generic output names rather than phenix naming conventions.
New features that inspect files or match output patterns by name can break these
tests invisibly — passing their own unit tests while silently breaking all callers.

---

## Step 6.5: Zombie State Checker (if your program produces a model or map)

If your program produces a **model file** (`.pdb`, `.cif`) or **map file**
(`.mrc`, `.ccp4`) and uses `strategy: "set_flag"` or `strategy: "count"`, you
should add an entry to `_ZOMBIE_CHECK_TABLE` in `agent/workflow_state.py`.

### What is a zombie state?

A zombie state occurs when the agent crashes mid-cycle or output files are
deleted after the fact. History records `done_flag=True`, but the output file
is no longer on disk. The phase detector sees `done=True` and skips the
program — the workflow becomes stuck.

The zombie checker (`_clear_zombie_done_flags`) runs every PERCEIVE cycle. It
checks each done flag against an expected output file pattern. If the output is
missing, it clears the done flag so the program can re-run.

### When to add a zombie check entry

Add an entry when **all** of the following are true:
- The program produces an output file that later steps depend on
- The done flag is used to gate a workflow phase (e.g., `has_placed_model`,
  `has_refined_model`, or a count like `refine_count`)
- A crash or cleanup that removes the output would leave the workflow stuck

You **do not** need an entry for:
- Analysis programs (xtriage, mtriage) — their outputs are metrics, not files
- Programs with `strategy: "run_once"` for analysis — they don't gate file-dependent phases
- Programs where the done flag is advisory only (e.g., `validation_done`)

### How to add the entry

In `agent/workflow_state.py`, find `_ZOMBIE_CHECK_TABLE` and add a tuple:

```python
_ZOMBIE_CHECK_TABLE = [
    # ... existing entries ...

    # Your new program
    ("new_program_done",
     _re.compile(r"<output_pattern>.*\.pdb$", _re.IGNORECASE),
     "has_placed_model",   # context flag to also clear (or None)
     True),                # accept_any_pdb (see below)
]
```

**Fields:**
| Field | Type | Description |
|-------|------|-------------|
| `done_flag` | str | Key in `history_info` to check (e.g., `"dock_done"`) |
| `output_pattern` | regex | Matched against basenames of available files |
| `file_flag_to_clear` | str or None | Context flag also cleared when zombie detected (e.g., `"has_placed_model"`) |
| `accept_any_pdb` | bool | If True, any `.pdb` on disk suppresses zombie detection as a fallback |

### Pattern alignment with file_categories.yaml

The `output_pattern` regex should be consistent with the patterns in
`knowledge/file_categories.yaml` for the same output type. If the file category
uses `*dock*map*` as a pattern, the zombie regex should also match `dock.*map.*\.pdb`.

**Always check:** does the zombie regex match the filenames phenix actually produces,
as well as what the file category accepts? They must agree.

```bash
# Quick check: test your pattern against real phenix output names
python3 -c "
import re
pat = re.compile(r'your_pattern_here', re.IGNORECASE)
names = ['expected_output_001.pdb', 'alternate_name.pdb']
for n in names: print(n, bool(pat.search(n)))
"
```

### The `accept_any_pdb` flag

| Value | Use when |
|-------|----------|
| `True` | Program produces a **model** file (`.pdb`/`.cif`). Any PDB on disk is plausible evidence of success, preventing false zombie clears when files are renamed or tests use generic names. |
| `False` | Program produces a **specifically-named** file that downstream steps reference by name (e.g., `denmod_map.ccp4` from `resolve_cryo_em`). Pattern match must be exact. |

**Trade-off:** `accept_any_pdb=True` weakens zombie detection in the edge case
where an older model file from a *different* program is on disk when the newer
program's output was deleted. This is an acceptable heuristic — zombie detection
is best-effort, not a guarantee. Document this choice in a comment.

### Example: dock_in_map

```python
# dock_in_map produces *_docked.pdb, *dock*map*.pdb, placed_model*.pdb, or *_placed*.pdb
# Pattern mirrors the "docked" file category in file_categories.yaml.
# accept_any_pdb=True: any PDB is acceptable fallback (model renamed or test proxy)
("dock_done",
 _re.compile(r"docked.*\.pdb$|.*_docked.*\.pdb$|dock.*map.*\.pdb$|placed_model.*\.pdb$|.*_placed.*\.pdb$",
             _re.IGNORECASE),
 "has_placed_model", True),
```

### Count-based programs

For programs with `strategy: "count"` (refine, rsr), the zombie checker
automatically decrements the count when it clears the done flag:

```python
("refine_done",
 _re.compile(r"refine_\d+\.pdb$", _re.IGNORECASE),
 None,   # No separate context flag — refine_count is decremented automatically
 True),
```

The decrement logic is built into `_clear_zombie_done_flags`. You don't need to
add code for it — just set `file_flag_to_clear=None` and the function handles
`refine_count` and `rsr_count` automatically.

---

## Example: Adding phenix.map_symmetry

Here's the complete configuration for map_symmetry:

### programs.yaml

```yaml
phenix.map_symmetry:
  description: "Find point-group symmetry in cryo-EM map"
  category: analysis
  experiment_types: [cryoem]
  done_tracking:
    flag: "map_symmetry_done"
    strategy: "run_once"
    history_detection:
      markers: ["map_symmetry"]

  inputs:
    required:
      map:
        extensions: [.mrc, .ccp4, .map]
        flag: ""

  input_priorities:
    map:
      categories: [full_map, optimized_full_map, map]
      exclude_categories: [half_map]

  outputs:
    files:
      - pattern: "*.ncs_spec"
        type: ncs_spec
    metrics:
      - symmetry_type
      - ncs_copies
      - ncs_cc

  command: "phenix.map_symmetry {map}"

  # Metric extraction - used by log_parsers.py AND session.py
  log_parsing:
    symmetry_type:
      pattern: 'NCS\s+type[:\s]+([A-Z]\d+(?:\s*\([a-z]\))?)'
      type: string
      display_name: "Symmetry Type"
      summary_format: "{value}"
      no_match_pattern: 'No suitable symmetry found'
      no_match_value: "None"
    ncs_copies:
      pattern: 'Number of NCS copies[:\s]+(\d+)'
      type: int
      display_name: "NCS Copies"
      summary_format: "{value}"
    ncs_cc:
      pattern: 'NCS correlation[:\s]+([0-9.]+)'
      type: float
      display_name: "NCS CC"
      summary_format: "{value:.2f}"

  hints:
    - "Run to detect point-group symmetry in cryo-EM maps"
```

### workflows.yaml

```yaml
cryoem:
  phases:
    analyze:
      programs:
        - phenix.mtriage
        - program: phenix.map_symmetry
          conditions:
            - not_done: map_symmetry
          hint: "Run to detect map symmetry"
```

### metrics.yaml (summary_display section)

```yaml
summary_display:
  quality_table:
    - metrics: [symmetry_type]
      label: "Symmetry"
      format: "{symmetry_type}"
      optional_format: " ({ncs_copies} copies)"
      optional_metrics: [ncs_copies]
      detail_format: "CC: {ncs_cc:.2f}"
      detail_metrics: [ncs_cc]
      experiment_types: [cryoem]

  step_metrics:
    phenix.map_symmetry:
      format: "Symmetry: {symmetry_type}"
      metrics: [symmetry_type]
      optional_format: " ({ncs_copies} copies)"
      optional_metrics: [ncs_copies]
      fallback_format: "Analyzed"
```

That's it - **only 3 files modified**, and most of the content is in programs.yaml!

---

## Example: Adding phenix.map_correlations (Multi-Mode Inputs)

This example demonstrates a program with multiple input modes (model+map, map+map,
mtz+mtz) where `auto_fill: false` prevents the command builder from duplicating files.

### programs.yaml

```yaml
phenix.map_correlations:
  description: "Compute map-model or map-map correlation coefficients"
  category: validation
  experiment_types: [xray, cryoem]

  inputs:
    # All inputs are optional — different modes use different combinations
    optional:
      model:
        extensions: [.pdb, .cif]
        description: "Model file"
        flag: "input_files.model="
      full_map:
        extensions: [.ccp4, .mrc, .map]
        description: "Map 1"
        flag: "input_files.map_in_1="
      map2:
        extensions: [.ccp4, .mrc, .map]
        description: "Map 2 for map-map CC"
        flag: "input_files.map_in_2="
        auto_fill: false  # CRITICAL: prevents same map filling both slots
      data_mtz:
        extensions: [.mtz]
        description: "Map coefficients 1"
        flag: "input_files.map_coeffs_1="
      map_coeffs_2:
        extensions: [.mtz]
        description: "Map coefficients 2"
        flag: "input_files.map_coeffs_2="
        auto_fill: false  # Same reason as map2

  command: "phenix.map_correlations"  # Files appended via flags

  log_parsing:
    cc_mask:
      pattern: "CC_mask\\s*[=:]\\s*([0-9.]+)"
      type: float
      display_name: "CC_mask"
      summary_format: "{value:.3f}"
    cc_volume:
      pattern: "CC_volume\\s*[=:]\\s*([0-9.]+)"
      type: float
      display_name: "CC_volume"
      summary_format: "{value:.3f}"
    # ... additional metrics
```

**Key design decisions:**

1. **All inputs optional:** Different modes need different combinations. The command
   builder appends only the slots that have files assigned.

2. **`auto_fill: false` on secondary slots:** Without this, if the user provides one
   `.ccp4` map, the builder would auto-fill both `full_map` AND `map2` with the same
   file, producing a meaningless self-correlation command.

3. **Flag-based append:** No `{placeholder}` in the command template. Files are
   appended as `flag + filepath` (e.g., `input_files.model=/path/to/model.pdb`).

### metrics.yaml additions

```yaml
# In summary_display.step_metrics:
phenix.map_correlations:
  format: "CC_mask: {cc_mask:.3f}"
  metrics: [cc_mask]
  fallback_format: "Computed correlations"

# In summary_display.quality_table:
- metrics: [cc_mask]
  label: "CC_mask"
  format: "{cc_mask:.3f}"
  detail_format: "CC_volume: {cc_volume:.3f}"
  detail_metrics: [cc_volume]
  assessment: true
```

### Important: Pattern matching for reformatted metrics

The `format_metrics_report()` function transforms metric keys via
`.replace("_", " ").title()`, so `cc_mask` becomes `Cc Mask` in the result text
stored in cycle records. Any hardcoded regex patterns in `session.py` must match
both forms:

```python
# Matches both "CC_mask  : 0.849" (raw log) and "Cc Mask: 0.849" (report)
r'CC[_ ]?mask\s*[=:]\s*([0-9.]+)'  # with re.IGNORECASE
```

---

## Key Design Patterns

### Prerequisites (Auto-Run Dependencies)

If your program requires a value that must come from another program (e.g.,
resolution from mtriage), declare a `prerequisite` on the invariant:

```yaml
phenix.new_program:
  invariants:
    resolution:
      required_key: resolution
      source: session
      prerequisite: phenix.mtriage   # ← Auto-runs mtriage if resolution missing
```

When the command builder detects the invariant can't be satisfied, it
automatically builds and runs the prerequisite program first. The original
program runs on the next cycle with the newly-available value.

**Important:** The prerequisite program must be able to produce the required
value. The build node checks `skip_programs` directives — if the user skipped
the prerequisite, a clear error is returned instead of silently failing.

### File Selection: exclude_categories for Map Inputs

When a program accepts maps, always exclude half-maps from the full_map slot.
Half-maps bubble up to the parent `map` category, so without exclusion, a
half-map can be selected as a full_map via the category fallback:

```yaml
input_priorities:
  full_map:
    categories: [full_map, optimized_full_map, map]
    exclude_categories: [half_map]   # ← Prevents half-map selection
```

### Intermediate File Exclusions

If your program writes intermediate files to subdirectories, ensure they don't
get tracked as outputs. Add exclusions in three places:

1. **`phenix_ai/log_parsers.py`** — `extract_output_files()`: Add path/basename
   patterns to the exclusion checks
2. **`phenix_ai/utilities.py`** — `scan_log_for_files()`: Same exclusion patterns
3. **`programs/ai_agent.py`** — `_track_output_files()`: Add to
   `intermediate_patterns` list and/or `intermediate_basename_prefixes`

Example for predict_and_build's internal docking intermediates:
```python
intermediate_patterns = [
    '/local_dock_and_rebuild',   # predict_and_build intermediates
    '/local_rebuild',
]
intermediate_basename_prefixes = [
    'working_model',              # Never a valid output
]
```

Also update `knowledge/file_categories.yaml` to exclude intermediates from
semantic categories (e.g., add `working_model*` to the `docked` category's
`excludes` list).

### LLM Strategy Sanitization

If your program has parameters where the LLM commonly hallucinates values,
add sanitization in `_build_strategy()` in `agent/command_builder.py`.

Current sanitizations:
- **resolution**: Stripped from LLM strategy when no verified source exists
  (prevents hallucinated resolution values from bypassing the mtriage prerequisite)
- **atom_type** (autosol): Multi-atom values like "Se, S" are split — first atom
  stays in `atom_type`, extras move to `additional_atom_types`

### .eff File Generation

Programs that use old-style master_phil (like phenix.refine) have their .eff
files generated via `command_line_argument_interpreter()` which resolves
short-form PHIL paths to full scope paths. This is handled in
`_generate_eff_file()` in `programs/ai_agent.py` (Attempt 3).

---

## Architecture Notes

### How It Works

1. **Metric Extraction** (`knowledge/metric_patterns.py`)
   - Reads `log_parsing` from programs.yaml
   - Compiles patterns once, caches them
   - Used by both `log_parsers.py` and `session.py`

2. **Program Registration** (`knowledge/program_registration.py`)
   - Reads `done_tracking` blocks from programs.yaml
   - Provides `get_program_done_flag_map()` (all programs → done flags)
   - Provides `get_trackable_programs()` (strategy: "run_once" programs only)

3. **Summary Display** (`knowledge/summary_display.py`)
   - Reads `summary_display` from metrics.yaml
   - Formats quality table and step metrics
   - Used by `session.py` for summary generation

### Adding Complex Parsing

For programs that need complex parsing (tables, multi-line, etc.), you can still add hardcoded extractors in `log_parsers.py`. The YAML patterns are tried first, then hardcoded extractors fill in any missing metrics.

```python
# In log_parsers.py
def _extract_new_program_metrics(log_text):
    """Handle complex parsing that YAML can't express."""
    metrics = {}
    # Complex table parsing, multi-line context, etc.
    return metrics
```

Add the call in `extract_all_metrics()`:

```python
elif "new_program" in prog_lower:
    hardcoded_metrics = _extract_new_program_metrics(log_text)
    for k, v in hardcoded_metrics.items():
        if k not in metrics:  # Don't overwrite YAML results
            metrics[k] = v
```

---

## Validation

Use the program validator to check configuration:

```bash
# Validate a specific program
python agent/program_validator.py phenix.map_symmetry

# Validate all programs
python agent/program_validator.py --all

# List all defined programs
python agent/program_validator.py --list
```

---

## Debugging Tips

1. **Program not appearing in valid_programs?**
   - Check conditions in `workflows.yaml`
   - Verify file categorization is correct
   - Check if `_done` flag is blocking it

2. **Metrics not extracted?**
   - Test regex pattern against actual log output
   - Check pattern in `programs.yaml` `log_parsing` section
   - Run: `python -c "from knowledge.metric_patterns import extract_metrics_for_program; print(extract_metrics_for_program(open('log.txt').read(), 'phenix.new_program'))"`

3. **Metrics not in summary?**
   - Check `summary_display` section in `metrics.yaml`
   - Verify `experiment_types` filter matches your data
   - Check required metrics are being extracted

4. **Command built incorrectly?**
   - Check `programs.yaml` input definitions
   - Verify `input_priorities` categories match file categorization

5. **Wrong file selected for input?**
   - Add explicit `input_priorities` with `exclude_categories`
   - Common mistake: protein input picks ligand file → add `exclude_categories: [ligand, ligand_fit_output]`
   - Common mistake: protein input picks search model → add `exclude_categories: [search_model, predicted]`

---

## Common Mistakes to Avoid

| Mistake | Symptom | Fix |
|---------|---------|-----|
| Missing `input_priorities` | Wrong file selected | Add explicit priorities with excludes |
| No `exclude_categories` for protein/model input | Picks ligand or predicted model | Add `exclude_categories: [ligand, search_model]` |
| Forgetting `done_tracking` with `strategy: "run_once"` | Program runs multiple times | Add `done_tracking: {flag: ..., strategy: "run_once"}` |
| Wrong `experiment_types` | Program offered for wrong data | Set `[xray]`, `[cryoem]`, or both |
| Missing workflow entry | Program never offered | Add to appropriate phase in `workflows.yaml` |

### Example: Proper input_priorities for a refinement program

```yaml
phenix.my_refine:
  inputs:
    required:
      model:
        extensions: [.pdb, .cif]
        flag: ""
      data_mtz:
        extensions: [.mtz]
        flag: ""
  
  # CRITICAL: Prevent selecting wrong files
  input_priorities:
    model:
      categories: [model]
      prefer_subcategories: [refined, with_ligand, phaser_output]
      exclude_categories: [ligand, search_model, ligand_fit_output, predicted]
    data_mtz:
      categories: [data_mtz]
```
