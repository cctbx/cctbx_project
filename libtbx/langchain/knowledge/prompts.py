"""
Prompt templates for agent planning and command generation.

This module contains all prompts used by the agent for:
- Strategic planning (deciding what to do next)
- Command writing (constructing valid commands)
- Keyword extraction (getting program parameters)
- Documentation queries

Usage:
    from libtbx.langchain.knowledge import get_strategic_planning_prompt

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

    **User's Goal / Advice:**
    {project_advice}

    **Project Database (Current Scientific State):**
    {project_state}

    **Original Input Files (Ground Truth):**
    {original_files}

    **Project History:**
    {history}

    **FORBIDDEN PROGRAMS:**
    1. **Do NOT use `phenix.cosym`**: Use `phenix.xtriage`.
    2. **Do NOT use `pointless`**: Use `phenix.xtriage`.
    3. **Do NOT use `phenix.reindex`**: Use `phenix.xtriage` to diagnose first.
    4. **Do NOT use `phenix.predict_model` or `phenix.predict_chain`**: Use `phenix.predict_and_build` with the keyword `stop_after_predict=True`.
    5. **Do NOT try to include ligands in `predict_and_build`**: Fit ligands separately after initial model building using `phenix.ligandfit`.

    **LIGAND TOOLS (When user provides a ligand file):**
    - **For fitting a known ligand into density:** Use `phenix.ligandfit` with:
      - `model=<refined_model.pdb>` - the protein model
      - `data=<map_coeffs.mtz>` - MTZ file with map coefficients
      - `ligand=<ligand.pdb>` - the ligand coordinates
      - `input_labels="2FOFCWT PH2FOFCWT"` - REQUIRED for map coefficients from refinement
    - **Example:** `phenix.ligandfit model=refined.pdb data=refined.mtz ligand=lig.pdb input_labels="2FOFCWT PH2FOFCWT"`
    - **After ligandfit:** The output `ligand_fit_1.pdb` contains ONLY the fitted ligand. Before refinement, you MUST combine it with the protein model using `phenix.pdbtools`.
    - **For refinement with ligand:** Use `phenix.refine` with a SINGLE PDB file containing both protein AND ligand. Do NOT pass multiple PDB files.
    - **Ligand workflow:** Model building -> Initial refinement -> `phenix.ligandfit` -> Combine protein+ligand -> Final refinement with combined model.

    **ANALYSIS PROTOCOL (Follow in Order):**

    1. **TERMINATION CHECK (Stop cleanly):**
       - **SUCCESS:** If R-free < 0.25 (X-ray) or map CC is high (CryoEM) AND validation is good:
         -> **DECISION:** Output `STOP: Success`.
       - **MANUAL:** If R-factors are stalled/high and automation is failing:
         -> **DECISION:** Output `STOP: Manual Intervention`.

    2. **STATUS CHECK (Fix Crashes First):**
       - **Did the last job fail?** (e.g. "Sorry: ...", "Error", "Exception").
         -> **ACTION:** You MUST retry the *same* program (or a direct alternative) to fix the error.
         -> **CRITICAL:** Read the FULL error message carefully. It often tells you EXACTLY what parameters to add or change.
         -> **CRITICAL:** If error mentions "input labels" or "map type", add the suggested `input_labels=` or `lig_map_type=` parameters.
         -> **CRITICAL:** If error is about "R-free flags", add `xray_data.r_free_flags.generate=True`.
         -> **CRITICAL:** If error is "No search procedure", add `ncopies=1`.
       - **Did it succeed?**
         -> **ACTION:** Move to Protocol #3.

    3. **DEADLOCK CHECK (Symmetry):**
       - **Did a previous attempt to change the space group FAIL?** (e.g. "Incompatible unit cell").
       - **OR has Xtriage already been run multiple times?**
         -> **HARD STOP:** The lattice cannot support the higher symmetry.
         -> **REQUIRED ACTION:** Abandon re-indexing. Proceed using the **ORIGINAL (Input)** space group.
         -> **STRATEGY:** Explicitly handle this as a **TWINNING** case in the lower symmetry.

    4. **PIPELINE ORDER:**
       - **Standard:** Data Analysis (Xtriage) -> Phasing (Phaser/predict_and_build) -> Refinement -> Validation.
       - **With Ligand:** Data Analysis -> Phasing -> Initial Refinement -> **Ligand Fitting** -> Final Refinement -> Validation.
       - Do not skip to Refinement if Phasing failed or hasn't run.
       - **MODEL GENERATION:** If you have a sequence but no model, run `phenix.predict_and_build` (`stop_after_predict=True`).
       - **LIGAND FITTING:** If Original Files include a ligand file (e.g., `*_ligand.pdb`, `*_ligand.cif`, or user mentions "ligand"), run `phenix.ligandfit` AFTER initial refinement to place the ligand into density.

    **FILE REALITY PROTOCOL (Strict):**
    1. **Existing Files Only:** You can ONLY use files explicitly listed in "Original Input Files" OR the "Project History".
    2. **Missing Files:** If you need a file but it is NOT in the history:
       - **ACTION:** Do NOT stop. Identify the tool that *creates* that file and run THAT tool instead.

    **ANTI-LOOPING (Crucial):**
    - You CANNOT run the exact same command twice. If you retry a program, you MUST change at least one parameter.

    **OUTPUT RULES:**
    - **No Placeholders:** Provide exact values.
    - **Strategy Details:** Clean `key=value` pairs.

    **Output Format (Strict JSON):**
    {{
        "long_term_plan": "Step 1: [Done], Step 2: [Current], Step 3: [Goal]",
        "reasoning": "Explain decision. If fixing an error, be explicit.",
        "selected_program": "phenix.program_name",
        "strategy_details": "Exact parameters.",
        "input_files": "List of input files."
    }}
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
    Your task is to construct a valid command line for "{program}".

    **USER'S SPECIFIC INSTRUCTIONS (HIGHEST PRIORITY):**
    {project_advice}

    **The Strategy:**
    {strategy_details}

    **Input Files:**
    {input_files}

    **Original Input Files (Look here for paths):**
    {original_files}

    **Reference Documentation:**
    {valid_keywords}

    **LEARNED HISTORY (Mistakes from previous runs):**
    {learned_tips}

    **CRITICAL RULES:**
    1. **USER INSTRUCTIONS FIRST:** If the user gave specific instructions above (e.g., "use nproc=4", "use quick mode"), you MUST include those parameters.
    2. **PATH RESOLUTION (Crucial):** If a file in "Input Files" (e.g. `data.mtz`) matches a file in "Original Input Files" (e.g. `subdir/data.mtz`), you MUST use the FULL PATH from the Original list. Do NOT use bare filenames if a path is known.
    3. **Prioritize History:** If "LEARNED HISTORY" warns about a specific syntax error, you MUST follow that advice.
    4. **Follow the Reference:** Use the syntax rules found in the "Reference Documentation".
    5. **Syntax Priority:**
       - **Default:** Standard Phenix uses `key=value`.
       - **Exception:** If Reference says to use double-dashes `--flag`), use that.
    6. **Completeness:** Ensure every defined object has an action.

    **Output:**
    Provide ONLY the command string. No markdown.
    """
    return PromptTemplate(
        template=template,
        input_variables=[
            "program", "strategy_details", "input_files",
            "original_files", "valid_keywords", "learned_tips",
            "project_advice"
        ]
    )


def get_keywords_prompt() -> PromptTemplate:
    """
    Returns the prompt for extracting keywords and usage patterns from documentation.
    """
    template = """You are a Phenix documentation expert.
    Your task is to provide the critical information needed to construct a valid command for "{program_name}".

    **Goal:**
    The user needs to run this program for a specific task.
    You must extract:
    1. **Valid Parameters:** (e.g. `start_temperature=5000`)
    2. **Usage Examples:** (e.g. `phenix.refine data.mtz model.pdb annealing=True`) - THIS IS CRITICAL.

    **Instructions:**
    1. Scan the text for "Usage" or "Examples" sections.
    2. Look for patterns: how are input files specified? Do they use flags (e.g. `hklin=...`) or positional arguments?
    3. Look for the specific parameter syntax (e.g. is it `search.copies=1` or `ncopies=1`?).
    4. If broken lines appear (e.g. `param \\n = value`), reconstruct them to `param=value`.

    ---BEGIN CONTEXT---
    {context}
    ---END CONTEXT---

    **Output:**
    Provide a concise summary in two sections:

    SECTION 1: SYNTAX RULES & EXAMPLES
    (Paste any relevant command-line examples found in the text here. If none, describe the standard usage pattern.)

    SECTION 2: RELEVANT KEYWORDS
    (List of clean `key=value` strings found.)
    """
    return PromptTemplate(template=template, input_variables=["context", "program_name"])


def get_docs_query_prompt() -> PromptTemplate:
    """
    Returns the general-purpose prompt for querying the documentation RAG.
    """
    template = """You are an expert assistant for the Phenix software suite.
Answer the user's question based only on the following context from documentation and papers.
Provide detailed, helpful answers.
Do not discuss limitations of the sources.
Do not make up any answers.
You can say that you do not have enough information to reply.
Do not include any information on CryoFit, ShelxD, or CCP4 tools unless there
is a specific question about them.
Consider this question in the context of the process of structure determinaion in Phenix.
Focus on using Phenix tools, but include the use of Coot or Isolde if appropriate.
Name the tools that are to be used, along with their inputs and outputs and what they do.

Provide your answer based on the context and question below.

---BEGIN CONTEXT AND QUESTION---
Context:
{context}

Question:
{input}
"""
    return PromptTemplate(template=template, input_variables=["context", "input"])

