"""
Log file summarization using map-reduce approach.

This module handles:
- Chunking log files for processing
- Map phase: summarize each chunk
- Reduce phase: combine chunk summaries into final report

"""
from __future__ import absolute_import, division, print_function

import asyncio
from typing import List, Iterable
import os

from langchain_core.documents import Document
from langchain_core.prompts import PromptTemplate
from langchain_classic.chains.combine_documents import create_stuff_documents_chain
from langchain_text_splitters import RecursiveCharacterTextSplitter

from libtbx import group_args

# =============================================================================
# Debug flag - set to True to enable verbose output
# =============================================================================
DEBUG_SUMMARIZER = True #os.getenv("DEBUG_SUMMARIZER", "false").lower() == "true"

def debug_print(msg):
    """Print debug message if DEBUG_SUMMARIZER is enabled."""
    if DEBUG_SUMMARIZER:
        print(f"[SUMMARIZER DEBUG] {msg}")

# =============================================================================
# Program-Specific Summary Templates
# =============================================================================

PROGRAM_SPECIFIC_GUIDANCE = {
    "phenix.xtriage": """
**PROGRAM: phenix.xtriage - DATA QUALITY ANALYSIS**

This is a DATA ANALYSIS tool, NOT refinement. It does NOT produce Rwork/Rfree.

EXTRACT THESE METRICS:
- Completeness (%)
- Resolution range (e.g., "47.3 - 2.1 Å")
- I/sigma (overall signal-to-noise)
- Wilson B-factor
- Twin fraction estimates (from L-test, H-test, Britton analysis)
- Twin laws identified (e.g., "-h,-k,l")
- Space group suggestions

DO NOT EXTRACT OR REPORT:
- R-values, Rwork, Rfree (xtriage does NOT do refinement)
- The "R" values you may see (R_merge, R_sym, R_abs_twin = 0.527) are DATA QUALITY
  metrics measuring symmetry agreement, NOT refinement R-factors

OUTPUT FILES: None (xtriage is analysis only, produces no output files)
""",

    "phenix.phaser": """
**PROGRAM: phenix.phaser - MOLECULAR REPLACEMENT**

This places a search model into the unit cell using Patterson methods.

EXTRACT THESE METRICS:
- TFZ (Translation Function Z-score) - SUCCESS if > 8
- RFZ (Rotation Function Z-score)
- LLG (Log-Likelihood Gain) - higher is better
- eLLG (expected LLG)
- Number of copies placed
- Space group used/determined
- VRMS (coordinate error estimate)

DO NOT EXTRACT OR REPORT:
- Rwork/Rfree (phaser does NOT do refinement)

OUTPUT FILES: Look for PHASER.1.pdb, PHASER.1.mtz (numbered outputs)
""",

    "phenix.refine": """
**PROGRAM: phenix.refine - CRYSTALLOGRAPHIC REFINEMENT**

This refines atomic coordinates and B-factors against diffraction data.

EXTRACT THESE METRICS:
- Rwork and Rfree (THESE are the R-values to report)
- Bonds RMSD (target ~0.02 Å)
- Angles RMSD (target ~2.0°)
- Clashscore (lower is better)
- Ramachandran favored/outliers (%)
- Resolution used
- Number of atoms/residues

OUTPUT FILES: Look for *_001.pdb, *_001.mtz, *_data.mtz
""",

    "phenix.predict_and_build": """
**PROGRAM: phenix.predict_and_build - ALPHAFOLD PREDICTION & MODEL BUILDING**

This predicts structure with AlphaFold and optionally builds into density.

EXTRACT THESE METRICS:
- pLDDT scores (model confidence, 0-100, higher is better)
- Number of residues predicted/built
- Sequence coverage
- Map correlation (if building was performed)

DO NOT EXTRACT OR REPORT:
- Rwork/Rfree (unless full refinement was included in the run)

OUTPUT FILES: Look for *_predicted_model.pdb, *_predicted_*.pdb in output directories
""",

    "phenix.autobuild": """
**PROGRAM: phenix.autobuild - AUTOMATED MODEL BUILDING**

This automatically builds a model into electron density.

EXTRACT THESE METRICS:
- R-factor and R-free (if refinement included)
- Number of residues built
- Map correlation
- Completeness of model
- Number of fragments/chains

OUTPUT FILES: Look for overall_best.pdb, overall_best_denmod_map_coeffs.mtz
""",

    "phenix.ligandfit": """
**PROGRAM: phenix.ligandfit - LIGAND FITTING**

This fits a ligand into difference density.

EXTRACT THESE METRICS:
- Correlation coefficient (CC) of ligand fit
- Number of atoms placed
- Real-space R-factor for ligand

OUTPUT FILES: Look for ligand_fit_*.pdb, *_fitted.pdb
""",

    "phenix.molprobity": """
**PROGRAM: phenix.molprobity - STRUCTURE VALIDATION**

This validates model geometry and quality.

EXTRACT THESE METRICS:
- Clashscore
- Ramachandran favored/allowed/outliers (%)
- Rotamer outliers (%)
- C-beta deviations
- Bad bonds/angles count
- MolProbity score (overall quality, lower is better)

OUTPUT FILES: None (validation produces reports, not models)
""",

    "phenix.ready_set": """
**PROGRAM: phenix.ready_set - MODEL PREPARATION**

This prepares a model for refinement (adds hydrogens, restraints).

EXTRACT THESE METRICS:
- Number of atoms added
- Ligand codes processed
- Restraint files generated

OUTPUT FILES: Look for *.updated.pdb, *.cif (restraint files)
""",

    "phenix.mtriage": """
**PROGRAM: phenix.mtriage - CRYO-EM MAP ANALYSIS**

This analyzes cryo-EM map quality.

EXTRACT THESE METRICS:
- Resolution estimates (d99, d_model, d_FSC)
- Map-model correlation
- Half-map correlation (FSC)

DO NOT EXTRACT OR REPORT:
- Rwork/Rfree (this is map analysis, not refinement)

OUTPUT FILES: None (analysis only)
""",
}


def get_program_specific_guidance(program_name: str) -> str:
    """
    Returns program-specific extraction guidance for the given program.

    Args:
        program_name: Name of the Phenix program (e.g., "phenix.xtriage")

    Returns:
        str: Specific guidance for that program, or generic guidance
    """
    if not program_name:
        return ""

    program_lower = program_name.lower().strip()

    # Direct match
    if program_lower in PROGRAM_SPECIFIC_GUIDANCE:
        return PROGRAM_SPECIFIC_GUIDANCE[program_lower]

    # Partial match (e.g., "xtriage" matches "phenix.xtriage")
    for key, guidance in PROGRAM_SPECIFIC_GUIDANCE.items():
        if key.replace("phenix.", "") in program_lower or program_lower in key:
            return guidance

    # Generic fallback
    return """
**PROGRAM: Unknown/Generic**

Extract standard crystallographic metrics if present:
- Resolution, completeness, R-values (only if this is a refinement program)
- Space group information
- Any warnings or errors
- Output files explicitly mentioned as written

Be conservative - only report metrics you're confident are present.
"""


def detect_program_from_text(text: str) -> str:
    """
    Detects the Phenix program name from log text.
    
    Checks multiple patterns:
    1. Client-added header "COMMAND THAT WAS RUN: phenix.xxx ..."
    2. Standard Phenix log header "Starting phenix.xxx on ..."
    3. Program-specific content patterns
    
    Args:
        text: Log file content
        
    Returns:
        Program name (e.g., "phenix.refine") or empty string
    """
    if not text:
        return ""
    
    # Check first 5000 chars for efficiency
    text_start = text[:5000]
    text_start_lower = text_start.lower()
    
    # Pattern 1: Client-added header "COMMAND THAT WAS RUN: phenix.xxx ..."
    for line in text_start.splitlines():
        if "COMMAND THAT WAS RUN:" in line:
            parts = line.split()
            for p in parts:
                if p.startswith("phenix."):
                    debug_print(f"Detected program from COMMAND header: {p}")
                    return p
    
    # Pattern 2: Standard Phenix log header "Starting phenix.xxx on ..."
    for line in text_start.splitlines():
        line_stripped = line.strip()
        if line_stripped.startswith("Starting phenix."):
            parts = line_stripped.split()
            for p in parts:
                if p.startswith("phenix."):
                    debug_print(f"Detected program from 'Starting' header: {p}")
                    return p
    
    # Pattern 3: Program-specific content patterns
    # Order matters - check more specific patterns first
    program_patterns = [
        ("phenix.predict_and_build", ["predict_and_build", "phenix.predict_and_build", "alphafold prediction"]),
        ("phenix.refine", ["phenix.refine", "refinement.input", "macro_cycle", "r_work", "r_free"]),
        ("phenix.phaser", ["phenix.phaser", "phaser.", "translation function z-score", "molecular replacement", "tfz="]),
        ("phenix.xtriage", ["phenix.xtriage", "xtriage", "scaling.input", "wilson statistics", "twinning analysis"]),
        ("phenix.autobuild", ["phenix.autobuild", "autobuild", "autobuild_cycle"]),
        ("phenix.ligandfit", ["phenix.ligandfit", "ligandfit", "ligand fitting", "fitting ligand"]),
        ("phenix.molprobity", ["phenix.molprobity", "molprobity", "clashscore", "ramachandran"]),
        ("phenix.ready_set", ["phenix.ready_set", "ready_set", "readyset"]),
        ("phenix.mtriage", ["phenix.mtriage", "mtriage", "half-map", "d_fsc"]),
        ("phenix.real_space_refine", ["phenix.real_space_refine", "real_space_refine"]),
        ("phenix.dock_in_map", ["phenix.dock_in_map", "dock_in_map"]),
        ("phenix.autosol", ["phenix.autosol", "autosol", "sad phasing", "experimental phasing"]),
    ]
    
    for program_name, patterns in program_patterns:
        for pattern in patterns:
            if pattern in text_start_lower:
                debug_print(f"Detected program from content pattern '{pattern}': {program_name}")
                return program_name
    
    debug_print("Could not detect program from text")
    return ""

# =============================================================================
# Prompt Templates
# =============================================================================

def get_log_map_prompt(program_name: str = "") -> PromptTemplate:
    """
    Returns the prompt for summarizing a single chunk of a log file.
    
    Args:
        program_name: Optional program name for specific guidance
        
    Returns:
        PromptTemplate configured for the program type
    """
    program_guidance = get_program_specific_guidance(program_name)
    
    template = f"""You are a crystallography data extraction expert analyzing a Phenix log file.

{program_guidance}

**CRITICAL INSTRUCTIONS:**
1. Analyze the ENTIRE log text below - not just the beginning or end
2. Extract ONLY information that is explicitly present in the log
3. Do NOT invent or assume information not shown
4. For output files, only list files explicitly shown as "written to" or "Writing:" or "Output:"

**EXTRACT THE FOLLOWING (if present):**

1. PROGRAM NAME: The specific phenix program (e.g., phenix.xtriage, phenix.phaser)
2. INPUT FILES: File paths for .mtz, .pdb, .fa, .seq files
3. OUTPUT FILES: ONLY files explicitly written (look for "Writing:", "written to", "Output:")
4. SPACE GROUP: Look for "Space group:" or "SPACEGROUP"
5. KEY METRICS: Based on the program-specific guidance above
6. WARNINGS/ERRORS: Lines with "Warning:", "Sorry:", "Error:", or tracebacks

**FORMAT:** Respond with a bulleted list using these exact section headers.

---
LOG CHUNK TO ANALYZE:
"{{text}}"
---

Extracted information:
"""
    return PromptTemplate.from_template(template)

def get_log_combine_prompt(program_name: str = "") -> PromptTemplate:
    """
    Returns the prompt for combining/summarizing log chunks.

    Args:
        program_name: Optional program name for specific guidance

    Returns:
        PromptTemplate for the reduce phase
    """
    program_guidance = get_program_specific_guidance(program_name)

    template = f"""You are an expert crystallographer creating a final summary report.

{program_guidance}

**CRITICAL RULES:**

1. **OUTPUT FILES**: ONLY list files that appear with "Writing:", "written to", "Output:",
   or "HKLOUT" in the summaries. If none are explicitly mentioned, write "None".
   NEVER invent filenames like "refined_001.mtz" unless explicitly shown.

2. **METRICS**: Only report metrics appropriate for this program type (see guidance above).
   Do NOT report Rwork/Rfree for programs that don't do refinement.

3. **R-VALUES CONFUSION**:
   - For phenix.xtriage: R_merge, R_sym, R_abs_twin are DATA QUALITY metrics, NOT refinement R-factors
   - Only phenix.refine produces Rwork/Rfree
   - Do NOT report "R-values" for xtriage, phaser, or other non-refinement programs

---

Combine these intermediate summaries into ONE final report:

{{context}}

---

**FINAL REPORT FORMAT:**

1. **Input Files:** (List filenames with type: X-ray, Cryo-EM, Sequence, Model)

2. **Program Name:** (The specific Phenix program)

3. **Key Steps:** (2-4 bullet points of major operations)

4. **Key Metrics:** (Table format with metrics appropriate for this program)
   | Metric | Value |
   |--------|-------|
   | ...    | ...   |

5. **Key Space Group:** (Space group and any recommendations)

6. **Warnings/Errors:** (List any issues, or "None")

7. **Key Output Files:** (ONLY explicitly written files, or "None")

Be concise and accurate. Do not invent information.
"""
    return PromptTemplate(template=template, input_variables=["context"])

# =============================================================================
# Helper Functions
# =============================================================================

def get_chunk_size(provider: str = None):
    """
    Returns appropriate chunk size based on provider.

    Args:
        provider: 'google' or 'openai' or 'ollama'

    Returns:
        tuple: (chunk_size, chunk_overlap)
    """
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    provider = provider.lower()
    if provider == "openai":
        chunk_size = 100000
        chunk_overlap = 10000
        debug_print(f"OpenAI: chunk_size={chunk_size}, overlap={chunk_overlap}")
    elif provider == "google":
        chunk_size = 750000
        chunk_overlap = 50000
        debug_print(f"Google: chunk_size={chunk_size}, overlap={chunk_overlap}")
    elif provider == "ollama":
        chunk_size = 80000
        chunk_overlap = 8000
        debug_print(f"Ollama: chunk_size={chunk_size}, overlap={chunk_overlap}")
    else:
        chunk_size = 100000
        chunk_overlap = 10000
        debug_print(f"Unknown provider '{provider}': using default chunk_size={chunk_size}")
    return chunk_size, chunk_overlap


def _iter_batches(seq: List[Document], size: int) -> Iterable[List[Document]]:
    """Yield successive batches from sequence."""
    for i in range(0, len(seq), size):
        yield seq[i:i+size]


def _custom_log_chunker(log_text: str, provider: str = "google") -> List[Document]:
    """
    Chunks a log text with a special rule for the 'Files are in the directory' section.

    Args:
        log_text: The full log file text
        provider: 'google' or 'openai' or 'ollama'

    Returns:
        List of Document objects, each containing a chunk
    """
    chunk_size, chunk_overlap = get_chunk_size(provider)

    debug_print(f"Input log_text length: {len(log_text)} chars")

    trigger_phrase = "Files are in the directory"
    end_phrase = "Citations"

    final_chunks = []
    start_index = log_text.lower().find(trigger_phrase.lower())

    if start_index == -1:
        debug_print("No 'Files are in the directory' section found - using standard splitting")
        documents = [Document(page_content=log_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        chunks = standard_splitter.split_documents(documents)
        debug_print(f"Standard splitting produced {len(chunks)} chunks")
        for i, chunk in enumerate(chunks):
            debug_print(f"  Chunk {i+1}: {len(chunk.page_content)} chars, starts with: {chunk.page_content[:100]!r}...")
        return chunks

    debug_print(f"Found 'Files are in the directory' at position {start_index}")

    before_text = log_text[:start_index].strip()
    if before_text:
        debug_print(f"Before section: {len(before_text)} chars")
        before_doc = [Document(page_content=before_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        before_chunks = standard_splitter.split_documents(before_doc)
        debug_print(f"Before section split into {len(before_chunks)} chunks")
        final_chunks.extend(before_chunks)

    special_section_text = log_text[start_index:]
    end_index = special_section_text.find(end_phrase)

    if end_index != -1:
        special_chunk_content = special_section_text[:end_index].strip()
        debug_print(f"Special section (before '{end_phrase}'): {len(special_chunk_content)} chars")
    else:
        special_chunk_content = special_section_text.strip()
        debug_print(f"Special section (to end): {len(special_chunk_content)} chars")

    if special_chunk_content:
        final_chunks.append(Document(page_content=special_chunk_content))

    debug_print(f"Total chunks after custom splitting: {len(final_chunks)}")
    for i, chunk in enumerate(final_chunks):
        debug_print(f"  Chunk {i+1}: {len(chunk.page_content)} chars")
        debug_print(f"    Starts: {chunk.page_content[:80]!r}...")
        debug_print(f"    Ends: ...{chunk.page_content[-80:]!r}")

    return final_chunks


# =============================================================================
# Main Summarization Function
# =============================================================================

async def summarize_log_text(
    text: str,
    llm,
    timeout: int = 120,
    batch_size: int = 3,
    pause_between_batches: int = 1,
    use_throttling: bool = True,
    provider: str = None,
    program_name: str = None,
):

    """
    Performs a map-reduce summarization with batching to respect API rate limits.

    Args:
        text: The log file content to summarize
        llm: Language model to use for summarization
        timeout: Timeout in seconds for each batch
        batch_size: Number of chunks to process in parallel
        pause_between_batches: Seconds to pause between batches
        use_throttling: Whether to pause between batches
        provider: 'google' or 'openai' or 'ollama'
        program_name: Optional program name for specific guidance (e.g., 'phenix.xtriage').
                      If not provided, will be auto-detected from log content.

    Returns:
        group_args with:
            - group_args_type: 'log_summary'
            - log_summary: The final summary text (or None if failed)
            - error: Error message (or None if successful)
    """

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    min_log_length = 500  # Minimum chars for a real log file

    # Check for initialization markers that indicate no real log
    initialization_markers = [
        "Project initialized",
        "Starting analysis of provided files",
        "Startup Advice",
        "No log file",
    ]

    # Also check that the text does NOT contain real log content markers
    real_log_markers = [
        "phenix.",
        "Working directory:",
        "COMMAND THAT WAS RUN:",
        "Resolution range:",
        "Space group:",
        "R-work",
        "R-free",
    ]

    text_lower = text.lower()
    text_start = text_lower[:1000]

    is_initialization = any(marker.lower() in text_start for marker in initialization_markers)
    has_real_content = any(marker.lower() in text_lower for marker in real_log_markers)
    is_too_short = len(text.strip()) < min_log_length

    # Skip if it looks like initialization AND doesn't have real log content
    if (is_initialization or is_too_short) and not has_real_content:
        debug_print(f"Skipping summarization - initialization text or too short ({len(text)} chars)")
        debug_print(f"  is_initialization={is_initialization}, has_real_content={has_real_content}, is_too_short={is_too_short}")
        return group_args(
            group_args_type='log_summary',
            log_summary="No log file to analyze. This is the project initialization phase.",
            error=None
        )

    # NEW: Auto-detect program if not provided
    if not program_name:
        program_name = detect_program_from_text(text)
        debug_print(f"Auto-detected program: {program_name or 'unknown'}")
    else:
        debug_print(f"Using EXPLICIT program_name: {program_name}")  # Make it clear

    debug_print("=" * 60)
    debug_print(f"SUMMARIZE_LOG_TEXT STARTING")
    debug_print(f"  Provider: {provider}")
    debug_print(f"  LLM model: {getattr(llm, 'model', getattr(llm, 'model_name', 'unknown'))}")
    debug_print(f"  Input text length: {len(text)} chars")
    debug_print(f"  Timeout: {timeout}s, Batch size: {batch_size}")
    debug_print("=" * 60)

    docs = _custom_log_chunker(text, provider=provider)
    if not docs:
        debug_print("ERROR: No chunks produced!")
        return group_args(
            group_args_type='log_summary',
            log_summary=None,
            error="Log file produced no content to summarize."
        )

    # NEW: Use program-specific map prompt
    map_prompt = get_log_map_prompt(program_name)
    map_chain = map_prompt | llm

    debug_print(f"\n--- MAP PHASE: Summarizing {len(docs)} chunks ---")

    all_intermediate_summaries = []

    num_batches = (len(docs) + batch_size - 1) // batch_size
    print(f"Summarizing {len(docs)} chunks in {num_batches} batches to respect rate limits...")

    for i, batch in enumerate(_iter_batches(docs, batch_size)):
        print(f"  - Processing batch {i + 1} of {num_batches}...")
        debug_print(f"\n  BATCH {i+1}/{num_batches}: {len(batch)} chunks")

        for j, doc in enumerate(batch):
            debug_print(f"    Chunk {j+1} input ({len(doc.page_content)} chars):")
            debug_print(f"      First 200 chars: {doc.page_content[:200]!r}")
            debug_print(f"      Last 200 chars: ...{doc.page_content[-200:]!r}")

        tasks = [map_chain.ainvoke({"text": doc.page_content}) for doc in batch]

        try:
            map_results = await asyncio.wait_for(asyncio.gather(*tasks), timeout=timeout)
            
            for j, result in enumerate(map_results):
                summary_text = result.content
                all_intermediate_summaries.append(summary_text)
                debug_print(f"\n    Chunk {j+1} MAP OUTPUT ({len(summary_text)} chars):")
                debug_print("-" * 40)
                # Print full intermediate summary for debugging
                for line in summary_text.split('\n')[:30]:  # First 30 lines
                    debug_print(f"      {line}")
                if summary_text.count('\n') > 30:
                    debug_print(f"      ... ({summary_text.count(chr(10)) - 30} more lines)")
                debug_print("-" * 40)

        except asyncio.TimeoutError:
            print(f"Batch {i + 1} timed out after {timeout} seconds. Proceeding with partial results.")
            debug_print(f"  BATCH {i+1} TIMEOUT!")
            continue

        if use_throttling and (i < num_batches - 1):
            print(f"  - Pausing {pause_between_batches} seconds ...")
            await asyncio.sleep(pause_between_batches)

    if not all_intermediate_summaries:
        debug_print("ERROR: No intermediate summaries produced!")
        return group_args(
            group_args_type='log_summary',
            log_summary=None,
            error="Log summarization failed to produce any results."
        )

    debug_print(f"\n--- REDUCE PHASE: Combining {len(all_intermediate_summaries)} summaries ---")

    summary_docs = [Document(page_content=s) for s in all_intermediate_summaries]

    # Debug: show what's being combined
    debug_print("Intermediate summaries being combined:")
    for i, doc in enumerate(summary_docs):
        debug_print(f"\n  === INTERMEDIATE SUMMARY {i+1} ({len(doc.page_content)} chars) ===")
        debug_print(doc.page_content)
        debug_print("  === END ===\n")

    # NEW: Use program-specific combine prompt
    combine_prompt = get_log_combine_prompt(program_name)

    
    # Debug: show the combine prompt
    debug_print("COMBINE PROMPT TEMPLATE:")
    debug_print(combine_prompt.template[:500] + "...")

    reduce_chain = create_stuff_documents_chain(llm, combine_prompt)
    final_output = reduce_chain.invoke({"context": summary_docs})

    debug_print(f"\n--- FINAL OUTPUT ({len(final_output)} chars) ---")
    debug_print("=" * 60)
    debug_print(final_output)
    debug_print("=" * 60)

    return group_args(
        group_args_type='log_summary',
        log_summary=final_output,
        error=None
    )

