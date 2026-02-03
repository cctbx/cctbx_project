# Error Recovery System

The PHENIX AI Agent includes automatic error recovery for certain well-defined error conditions. When a program fails with a recoverable error, the agent can automatically determine the correct fix and retry.

## Overview

Some PHENIX programs fail with structured, prescriptive error messages that can be automatically resolved. For example, when an MTZ file contains multiple data arrays, the program tells you exactly which parameter to use and what choices are available.

The error recovery system:
1. **Detects** known error patterns in program output
2. **Extracts** structured information (choices, keywords, affected files)
3. **Resolves** the error using context-aware logic
4. **Retries** the program with the fix applied

## Supported Errors

### 1. Ambiguous Data Labels

**When it occurs:** MTZ file contains multiple suitable data arrays (e.g., both merged intensities and anomalous pairs)

**Error message:**
```
Multiple equally suitable arrays of observed xray data found.

Possible choices:
  /path/to/data.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
  /path/to/data.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-),merged

Please use scaling.input.xray_data.obs_labels
to specify an unambiguous substring of the target label.
```

**How it's resolved:**
- For SAD/MAD phasing programs (autosol): Select anomalous data `I(+)/I(-)`
- For standard refinement: Select merged intensities `IMEAN`
- User advice like "MRSAD" or "anomalous" influences the selection

### 2. Ambiguous Experimental Phases

**When it occurs:** MTZ file contains both standard and anomalous HL coefficients

**Error message:**
```
Multiple equally suitable arrays of experimental phases found.

Possible choices:
  /path/to/data.mtz:HLAM,HLBM,HLCM,HLDM
  /path/to/data.mtz:HLanomA,HLanomB,HLanomC,HLanomD
```

**How it's resolved:**
- For refinement: Select standard HL coefficients
- For SAD phasing: Select anomalous HL coefficients

## How It Works

### Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     Program Execution                            │
│                            │                                     │
│                    [Program FAILED]                              │
│                            │                                     │
│                            ▼                                     │
│              ┌─────────────────────────┐                        │
│              │     ErrorAnalyzer       │                        │
│              │  analyze(log_text,      │                        │
│              │          program,       │                        │
│              │          context)       │                        │
│              └───────────┬─────────────┘                        │
│                          │                                       │
│           ┌──────────────┼──────────────┐                       │
│           ▼              ▼              ▼                       │
│    [Not Recoverable] [Recoverable]  [Max Retries]               │
│           │              │              │                       │
│           ▼              ▼              ▼                       │
│      Continue        Save to        Report &                    │
│      (normal)        Session        Continue                    │
│                          │                                       │
│                          │  • recovery_strategies[file] = fix   │
│                          │  • force_retry_program = program     │
└──────────────────────────┼───────────────────────────────────────┘
                           │
                           ▼
                    [Next Cycle: Retry with fix applied]
```

### Key Components

| Component | File | Purpose |
|-----------|------|---------|
| ErrorAnalyzer | `agent/error_analyzer.py` | Detect and resolve errors |
| Error Patterns | `knowledge/recoverable_errors.yaml` | YAML-defined error signatures |
| Session Storage | `agent/session.py` | Store recovery strategies per file |
| Command Builder | `agent/command_builder.py` | Apply fixes to commands |

### YAML Configuration

Error patterns are defined in `knowledge/recoverable_errors.yaml`:

```yaml
errors:
  ambiguous_data_labels:
    description: "MTZ file contains multiple suitable data arrays"
    max_retries: 3
    
    detection_patterns:
      - "Multiple equally suitable arrays"
      - "Please use.*to specify an unambiguous substring"
    
    # Extract parameter name from error message
    keyword_extraction: 'Please use\s+(\S+)\s+to specify'
    
    # Extract available choices
    choice_extraction: '^\s*(\S+\.mtz):(.+)$'
    
    resolution: context_based_label_selection

label_patterns:
  anomalous_indicators:
    - '\(\+\)'      # I(+), F(+)
    - '\(-\)'       # I(-), F(-)
    - 'anom'        # Fanom, Ianom
```

## Usage

### Automatic (Default)

Error recovery is enabled by default. No configuration needed:

```bash
phenix.ai_agent original_files="data.mtz sequence.fa"
```

If a recoverable error occurs, the agent will:
1. Detect the error type
2. Determine the appropriate fix
3. Retry the program automatically

### Disable Recovery

For debugging or when you want manual control:

```bash
phenix.ai_agent auto_recovery=False original_files="data.mtz sequence.fa"
```

### View Recovery Actions

The session summary shows what recoveries occurred:

```
============================================================
SESSION SUMMARY
============================================================
Cycle 2: phenix.autosol - FAILED
         [Recovery: Selected anomalous data I_CuKa(+) for phasing]
Cycle 3: phenix.autosol - SUCCESS (retry with fix)
```

## Example Walkthrough

### Scenario: SAD Phasing with Ambiguous Data

**Cycle 1:** User starts workflow
```bash
phenix.ai_agent original_files="lyso.mtz lyso.fa" \
    project_advice="Solve using SAD phasing"
```

**Cycle 2:** Agent runs xtriage (succeeds, detects anomalous signal)

**Cycle 3:** Agent runs autosol, which fails:
```
phenix.autosol lyso.mtz lyso.fa
→ FAILED: Multiple equally suitable arrays of observed xray data found.
   Choices: IMEAN_CuKa, I_CuKa(+)...
```

**Error Analysis:**
1. ErrorAnalyzer detects "Multiple equally suitable arrays"
2. Extracts keyword: `scaling.input.xray_data.obs_labels`
3. Extracts choices: `["IMEAN_CuKa", "I_CuKa(+)"]`
4. Context: program=autosol, advice="SAD phasing"
5. Selection: `I_CuKa(+)` (anomalous data needed for SAD)

**Cycle 4:** Agent retries autosol with fix:
```
phenix.autosol lyso.mtz lyso.fa scaling.input.xray_data.obs_labels="I_CuKa(+)"
→ SUCCESS
```

## Design Principles

1. **Extract from error message** - Don't hardcode parameter names; parse them from the error
2. **Context-aware resolution** - Use program type and user advice to choose correctly
3. **Max 3 retries** - Prevent infinite loops on unrecoverable situations
4. **Per-file tracking** - Store recovery strategies keyed by filename
5. **Strategy injection** - Use existing parameter system, not command string hacking
6. **YAML-driven** - Error patterns externalized for easy maintenance
7. **Opt-out available** - `auto_recovery=False` for debugging

## Adding New Recoverable Errors

To add support for a new error type:

1. **Add pattern to YAML** (`knowledge/recoverable_errors.yaml`):
```yaml
errors:
  new_error_type:
    description: "Description of the error"
    max_retries: 3
    detection_patterns:
      - "pattern to match in log"
    keyword_extraction: 'regex to extract parameter name'
    resolution: resolution_function_name
```

2. **Add resolution logic** (`agent/error_analyzer.py`):
```python
def _resolve_new_error_type(self, error_info, context):
    # Determine the correct fix based on context
    return {
        "parameter_name": "correct_value",
    }
```

3. **Add tests** (`tests/test_error_analyzer.py`)

## Future Enhancements

Potential additional recoverable errors:

| Error | Detection | Potential Recovery |
|-------|-----------|-------------------|
| Missing resolution | "Please set the resolution" | Run mtriage/xtriage first |
| Space group mismatch | "Space group.*does not match" | May need user input |
| Missing sequence | "No sequence provided" | Check for .fa/.fasta files |
| Memory error | "Out of memory" | Reduce nproc, add memory flags |

## See Also

- [OVERVIEW.md](../OVERVIEW.md) - Technical architecture
- [USER_DIRECTIVES.md](USER_DIRECTIVES.md) - How user advice influences recovery decisions
- `knowledge/recoverable_errors.yaml` - Error pattern definitions
- `agent/error_analyzer.py` - Implementation
