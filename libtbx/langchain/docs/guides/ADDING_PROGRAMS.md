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
- ✅ `<program>_done` tracking flags (via `run_once: true`)
- ✅ Session summary display (default formatting works for most metrics)
- ✅ Workflow context building (automatic from programs.yaml)

---

## Minimal Example (2 files)

For a simple analysis program that runs once and extracts one metric:

**programs.yaml:**
```yaml
phenix.new_analysis:
  description: "Analyze something in the data"
  category: analysis
  experiment_types: [xray, cryoem]
  run_once: true  # Auto-creates new_analysis_done flag
  
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
  
  # Set run_once: true for programs that should only run once per session
  # This automatically creates <program>_done tracking flags
  run_once: true

  inputs:
    required:
      input_name:
        extensions: [.ext1, .ext2]
        flag: "input_flag="  # or "" for positional
    optional:
      optional_input:
        extensions: [.ext]
        flag: "optional_flag="

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
| `run_once: true` | Automatically creates `<program>_done` flag in workflow_state.py |
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
- `has: <file_type>` - Requires file (sequence, model, full_map, etc.)
- `has_any: [type1, type2]` - Requires at least one of these files
- `not_done: <program>` - Program hasn't been run yet (auto-generated for `run_once: true` programs)
- `not_has: <file_type>` - File type is NOT present

**Available done flags:**

For programs with `run_once: true`, flags are auto-generated (e.g., `xtriage_done`, `mtriage_done`).

For other programs, these flags are manually tracked in `workflow_state.py`:

| Flag | Set When | Used By |
|------|----------|---------|
| `predict_done` | predict_and_build runs (any mode) | cryo-EM workflow |
| `predict_full_done` | predict_and_build completes full workflow | X-ray workflow |
| `phaser_done` | phaser completes MR | X-ray workflow |
| `dock_done` | dock_in_map completes | cryo-EM workflow |
| `autobuild_done` | autobuild succeeds | X-ray workflow |
| `autobuild_denmod_done` | autobuild maps_only completes | X-ray workflow (before ligandfit) |
| `autosol_done` | autosol succeeds | X-ray workflow |
| `ligandfit_done` | ligandfit completes | Both workflows |
| `refine_done` | refine completes (any) | Internal tracking |
| `refine_count` | Count of successful refinements | Workflow conditions |
| `rsr_done` | real_space_refine completes | Internal tracking |
| `rsr_count` | Count of successful RSR runs | Workflow conditions |
| `resolve_cryo_em_done` | resolve_cryo_em completes | cryo-EM workflow |
| `map_sharpening_done` | map_sharpening completes | cryo-EM workflow |
| `map_to_model_done` | map_to_model completes | cryo-EM workflow |

**Note:** `refine_count` and `rsr_count` only count SUCCESSFUL runs (not failed ones).

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

# Run all tests
python tests/run_all_tests.py --quick
```

The validator checks:
- ✅ Program defined in programs.yaml
- ✅ Added to workflows.yaml
- ✅ Log parsing patterns defined
- ✅ Summary display configured (if applicable)

---

## Example: Adding phenix.map_symmetry

Here's the complete configuration for map_symmetry:

### programs.yaml

```yaml
phenix.map_symmetry:
  description: "Find point-group symmetry in cryo-EM map"
  category: analysis
  experiment_types: [cryoem]
  run_once: true  # Creates map_symmetry_done flag automatically

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

## Architecture Notes

### How It Works

1. **Metric Extraction** (`knowledge/metric_patterns.py`)
   - Reads `log_parsing` from programs.yaml
   - Compiles patterns once, caches them
   - Used by both `log_parsers.py` and `session.py`

2. **Program Registration** (`knowledge/program_registration.py`)
   - Reads programs with `run_once: true`
   - Auto-generates `<program>_done` flags
   - Updates `workflow_state.py` and `workflow_engine.py` automatically

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
| Forgetting `run_once: true` | Program runs multiple times | Add `run_once: true` |
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
      mtz:
        extensions: [.mtz]
        flag: ""
  
  # CRITICAL: Prevent selecting wrong files
  input_priorities:
    model:
      categories: [model]
      prefer_subcategories: [refined, with_ligand, phaser_output]
      exclude_categories: [ligand, search_model, ligand_fit_output, predicted]
    mtz:
      categories: [mtz]
      prefer_subcategories: [refined_mtz]
```
