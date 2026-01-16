"""
Prompt templates for agent planning and command generation.

This module contains all prompts used by the agent for:
- Strategic planning (deciding what to do next)
- Command writing (constructing valid commands)
- Keyword extraction (getting program parameters)
- Documentation queries

Usage:
    from libtbx.langchain.knowledge.prompts import get_strategic_planning_prompt

    prompt = get_strategic_planning_prompt()
"""
from __future__ import absolute_import, division, print_function

from langchain_core.prompts import PromptTemplate


def get_strategic_planning_prompt() -> PromptTemplate:
    """
    Returns the prompt for strategic planning - deciding what program to run next.
    """
    template = """You are a Lead Crystallographer supervising a structure solution project.
You have reports from {num_runs} previous jobs. Your goal is to determine the SINGLE best next step.
I will provide project advice enclosed in <project_advice> tags,
the project state enclosed in <project_state> tags,
the original input files enclosed in <original_files> tags,
and the project history enclosed in <project_history> tags.

**User's Goal / Advice:**
<project_advice>
{project_advice}
</project_advice>

**Project Database (Current Scientific State):**
<project_state>
{project_state}
</project_state>

**Original Input Files:**
<original_files>
{original_files}
</original_files>

**Project History:**
<history>
{history}
</history>

================================================================================
WORKFLOW RULES
================================================================================

**STANDARD X-RAY WORKFLOW (follow this order):**

1. phenix.xtriage (data quality analysis)
2. phenix.predict_and_build stop_after_predict=True (generate AlphaFold model)
3. phenix.process_predicted_model (prepare model for MR)
4. phenix.phaser (molecular replacement - places model in unit cell)
5. phenix.refine (refinement)
6. phenix.ligandfit (if ligand file provided)
7. phenix.pdbtools (combine protein + ligand)
8. phenix.refine (final refinement of complex)

**STANDARD CRYO-EM WORKFLOW:**

1. phenix.mtriage (map quality analysis)
2. phenix.predict_and_build stop_after_predict=True (generate AlphaFold model)
3. phenix.dock_in_map or phenix.em_placement (dock model into map)
4. phenix.real_space_refine (refinement)
5. phenix.ligandfit (if ligand file provided)

**FILE PROGRESSION (which files to use at each step):**

| After This Program | Use These Output Files |
|--------------------|------------------------|
| predict_and_build (stop_after_predict=True) | *_predicted_model.pdb |
| process_predicted_model | *_processed_*.pdb |
| phaser | PHASER.1.pdb, PHASER.1.mtz |
| refine | *_refine_001.pdb, *_refine_001.mtz |
| ligandfit | ligand_fit_1.pdb (combine with protein) |
| pdbtools (combine) | *_with_ligand.pdb or complex.pdb |

================================================================================
CRITICAL RULES
================================================================================

**PLACEMENT RULE:**
Models from `stop_after_predict=True` are in P1 (arbitrary orientation).
You MUST run phaser (X-ray) or dock_in_map (cryo-EM) BEFORE refinement.
NEVER run phenix.refine directly on *_predicted_model.pdb or *_rebuilt.pdb.

**ANTI-LOOPING RULE:**
- Never run the same program twice on the same input files
- Never go backwards in the workflow
- After successful refinement, proceed to ligandfit (if ligand exists) or STOP

**LIGAND RULE:**
- Never pass ligand files to predict_and_build (causes crashes)
- Never use ligand files as phaser ensemble models
- Ligands are fitted AFTER protein refinement using phenix.ligandfit

**SUCCESS CRITERIA:**
- X-ray: R-free < 0.30 is good, < 0.25 is excellent
- Cryo-EM: Map CC > 0.7 is good

================================================================================
FORBIDDEN PROGRAMS
================================================================================

Do NOT use these programs (use alternatives instead):
- phenix.cosym → use phenix.xtriage
- phenix.reindex → use phenix.xtriage
- phenix.predict_model → use phenix.predict_and_build with stop_after_predict=True
- phenix.model_vs_data → R-values from phenix.refine are sufficient

================================================================================
ERROR RECOVERY
================================================================================

If the last job FAILED, read the error message carefully:
- "R-free flags" error → add xray_data.r_free_flags.generate=True
- "No search procedure" → model file is missing or invalid
- "space group mismatch" → model needs placement (run phaser first)
- "input labels" error → add appropriate input_labels parameter

================================================================================
OUTPUT FORMAT
================================================================================

Respond with ONLY a valid JSON object (no markdown, no code blocks):

{{"long_term_plan": "Step 1: [Done], Step 2: [Current], Step 3: [Goal]", "reasoning": "Your explanation", "selected_program": "phenix.program_name", "strategy_details": "key=value pairs", "input_files": "file1.pdb file2.mtz"}}

If the project is complete or stuck, use: "selected_program": "STOP: Success" or "STOP: Manual Intervention"
"""

    return PromptTemplate(
        template=template,
        input_variables=["num_runs", "history", "project_advice", "original_files", "project_state"]
    )


def get_command_writer_prompt() -> PromptTemplate:
    """
    Returns the prompt for constructing valid Phenix command lines.
    """
    template = """You are a Phenix Command-Line Expert.
Construct a valid command line for "{program}".
I will provide project advice enclosed in <project_advice> tags,
the strategy_details enclosed in <strategy_details> tags,
suggested input files enclosed in <input_files> tags.
and all available input files enclosed in <available_files> tags.

**User Instructions:** <project_advice>{project_advice}</project_advice>
**Strategy:** <strategy_details>{strategy_details}</strategy_details>
**Input Files:** <input_files>{input_files}</input_files>
**Available Files:** <available_files>{original_files}</available_files>

================================================================================
REQUIRED PARAMETERS (MUST INCLUDE EXACTLY AS SHOWN)
================================================================================

{required_params}

If parameters are listed above, include each one exactly as shown.

================================================================================
VALID PARAMETERS FOR {program}
================================================================================

{valid_keywords}

IMPORTANT: Only use parameters from this list. Invalid parameters cause crashes.

================================================================================
LEARNED TIPS (from previous errors)
================================================================================

{learned_tips}

================================================================================
PROGRAM-SPECIFIC SYNTAX
================================================================================

**phenix.refine:**
```
phenix.refine model.pdb data.mtz [options]
```
- Model and data are positional (no prefixes)
- Add xray_data.r_free_flags.generate=True if needed

**phenix.phaser:**
```
phenix.phaser data.mtz model.pdb sequence.fa phaser.mode=MR_AUTO
```
- Requires a PDB model (not just sequence)
- Do NOT use search.copies=N (causes errors)

**phenix.ligandfit:**
```
phenix.ligandfit model=refined.pdb data=refined.mtz ligand=lig.pdb file_info.input_labels="2FOFCWT PH2FOFCWT"
```
- model= is REQUIRED (command fails without it)
- Use the REFINED model (*_refine_001.pdb), not PHASER.1.pdb
- Use general.nproc= for processors

**phenix.predict_and_build:**
```
phenix.predict_and_build input_files.seq_file=seq.fa input_files.xray_data_file=data.mtz predict_and_build.stop_after_predict=True
```
- Use input_files.xray_data_file= (not input_files.data_file=)
- NEVER include ligand files

**phenix.pdbtools (combine files):**
```
phenix.pdbtools protein.pdb ligand.pdb output.file_name=combined.pdb
```

**phenix.xtriage:**
```
phenix.xtriage data.mtz
```

**phenix.mtriage:**
```
phenix.mtriage map.ccp4
```

================================================================================
RULES
================================================================================

1. Use ACTUAL filenames from Input Files or Available Files
2. Never use placeholders like /path/to/, model.pdb, data.mtz
3. After ligand combination, refine the COMBINED file (*_with_ligand.pdb), not the old model
4. If a parameter isn't in the valid parameters list, omit it entirely

================================================================================
OUTPUT
================================================================================

Output ONLY the command starting with "phenix." - no explanations, no markdown.

Example: phenix.refine PHASER.1.pdb PHASER.1.mtz xray_data.r_free_flags.generate=True
"""

    return PromptTemplate(
        template=template,
        input_variables=[
            "program", "strategy_details", "input_files",
            "original_files", "valid_keywords", "learned_tips",
            "project_advice", "required_params"
        ]
    )


def get_keywords_prompt() -> PromptTemplate:
    """
    Returns the prompt for extracting keywords and usage patterns from documentation.
    """
    template = """Extract command-line information for "{program_name}" from the documentation below.

  I will provide documentation between <documentation_context> tags.

---BEGIN DOCUMENTATION---
<documentation_context>
{context}
</documentation_context>
---END DOCUMENTATION---

Provide TWO sections:

**SECTION 1: SYNTAX & EXAMPLES**
List any command-line examples found (e.g., `phenix.refine model.pdb data.mtz`)

**SECTION 2: VALID PARAMETERS**
List parameters as `key=value` format (e.g., `xray_data.high_resolution=2.0`)
"""
    return PromptTemplate(template=template, input_variables=["context", "program_name"])


def get_docs_query_prompt() -> PromptTemplate:
    """
    Returns the general-purpose prompt for querying the documentation RAG.
    """
    template = """Answer the question based on the Phenix documentation below.

I will provide context between <context> tags and the question between <input> tags

**Context:**
<context>{context}</context>

**Question:**
<input>{input}</input>

Focus on Phenix tools. Name specific programs with their inputs, outputs, and purpose.
If the information isn't in the context, say so.
"""
    return PromptTemplate(template=template, input_variables=["context", "input"])

