"""
A library of modular functions for building and querying LangChain RAG pipelines.
Required packages: google-generativeai langchain-cohere cohere langchain langchain-google-genai langchain-chroma pypdf unstructured beautifulsoup4 tqdm markdown-it-py langchain_community linkify-it-py

This module contains methods for setting up a RAG database from files
in a directory, querying the database, summarizing a log file, and
analyzing a log file summary in the context of the database.
"""

import nest_asyncio
nest_asyncio.apply()
import chromadb
import cohere
from cohere.core.api_error import ApiError as CohereApiError
from langchain_google_genai._common import GoogleGenerativeAIError
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
import openai # For handling OpenAI-specific exceptions

from bs4 import BeautifulSoup
from langchain_community.document_loaders.base import BaseLoader

import json
import os
import time
import asyncio
import glob
from concurrent.futures import TimeoutError
from libtbx import group_args
from google.api_core import exceptions as google_exceptions
from langchain_cohere import CohereRerank
from langchain.retrievers import ContextualCompressionRetriever
from langchain_core.runnables import RunnablePassthrough
from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAIEmbeddings

from langchain_community.document_loaders import TextLoader, PyPDFLoader, UnstructuredHTMLLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.chains.combine_documents import create_stuff_documents_chain
from langchain_core.prompts import PromptTemplate
from langchain_core.documents import Document
from langchain_chroma import Chroma
from typing import Iterable, List
from libtbx.langchain.run_analyze_log import save_as_html

try:
    from libtbx.langchain import phenix_knowledge as pk
except ImportError:
    import phenix_knowledge as pk

# --- GLOBAL VARIABLE to track the last query time ---
_last_query_time = 0

os.environ['GRPC_ENABLE_FORK_SUPPORT'] = "false"   # suppress a warning message

# --- LEARNING / MEMORY FUNCTIONS ---

def mtz_has_rfree_flags(mtz_file: str) -> bool:
    """
    Check if an MTZ file contains R-free flags.
    Returns True if R-free flags are present, False otherwise.
    """
    from libtbx import easy_run

    try:
        # Use phenix.mtz.dump to check contents
        result = easy_run.fully_buffered(f"phenix.mtz.dump {mtz_file}")
        output = "\n".join(result.stdout_lines)

        # Look for R-free-flags column
        # Common names: R-free-flags, FreeR_flag, FREE, RFREE
        rfree_indicators = [
            'R-free-flags',
            'FreeR_flag',
            'FREE',
            'RFREE',
            'R_FREE_FLAG',
            'test_flag'
        ]

        output_lower = output.lower()
        for indicator in rfree_indicators:
            if indicator.lower() in output_lower:
                return True

        return False

    except Exception as e:
        # If we can't check, assume no flags (safer)
        return False


def add_rfree_generation_if_needed(command: str, program: str) -> str:
    """
    For phenix.refine, check if the MTZ file has R-free flags.
    If not, automatically add xray_data.r_free_flags.generate=True
    """
    if 'phenix.refine' not in program:
        return command

    # Check if command already has r_free_flags parameter
    if 'r_free_flags.generate' in command:
        return command  # Already handled

    # Extract MTZ file from command
    parts = command.split()
    mtz_file = None
    for part in parts:
        if part.endswith('.mtz') and not part.startswith('-'):
            mtz_file = part
            break

    if not mtz_file:
        return command  # No MTZ file found

    # Check if MTZ has R-free flags
    if not mtz_has_rfree_flags(mtz_file):
        # Add the parameter
        command += " xray_data.r_free_flags.generate=True"

    return command

def validate_phenix_command(command: str) -> tuple[bool, str]:
    """
    Validates a Phenix command by running with --dry-run flag.
    This validates all parameters without actually running the job.
    """
    from libtbx import easy_run

    parts = command.split()
    if not parts:
        return False, "Empty command"

    program = parts[0]

    # Extract parameters (anything with '=')
    params_to_test = []
    for part in parts[1:]:
        if '=' in part:
            params_to_test.append(part)

    if not params_to_test:
        return True, ""

    # Build test command with --dry-run
    # This will validate parameters without running the actual job
    test_command = f"{program} --dry-run {' '.join(params_to_test)}"

    try:
        result = easy_run.fully_buffered(test_command)

        # Check both stderr and stdout for errors
        stderr_text = "\n".join(result.stderr_lines)
        stdout_text = "\n".join(result.stdout_lines)

        # Common error patterns
        error_patterns = [
            "Sorry:",
            "Error:",
            "ERROR:",
            "Ambiguous parameter",
            "not recognized",
            "Unknown parameter",
            "Invalid parameter",
            "Unrecognized PHIL"
        ]

        for pattern in error_patterns:
            if pattern in stderr_text or pattern in stdout_text:
                # Extract error context
                error_lines = []
                all_lines = result.stderr_lines + result.stdout_lines
                found_error = False
                for line in all_lines:
                    if pattern in line:
                        found_error = True
                    if found_error:
                        error_lines.append(line)
                        if len(error_lines) >= 15:  # Get more context
                            break

                error_msg = "\n".join(error_lines)
                return False, error_msg

        return True, ""

    except Exception as e:
        # If validation itself fails, be permissive
        return True, ""


def fix_command_syntax(command: str, error_message: str, llm) -> str:
    """
    Ask the LLM to fix a command based on the specific error message.
    Uses the SAME llm that's already loaded (expensive model).
    """
    fix_prompt = f"""You are a Phenix command-line expert. A command has a syntax error.

**Failed Command:**
{command}

**Error Message:**
{error_message}

**Your Task:**
Fix ONLY the syntax error. Do not change the logic or files.
Return ONLY the corrected command string, nothing else. No markdown, no explanation.

**Common Fixes:**
- "Ambiguous parameter" → Use full PHIL path (e.g., `xray_data.r_free_flags.generate=True` not `r_free_flags.generate=True`)
- "Ambiguous parameter" for twin_law → Use `refinement.main.twin_law=...`
- "not recognized" → Check spelling and PHIL hierarchy
- If the error mentions specific valid options, use one of those

**Critical:** Return ONLY the command, one line, no formatting.
"""

    try:
        response = llm.invoke(fix_prompt)
        fixed_command = response.content.strip()

        # Clean up any markdown formatting
        if "```" in fixed_command:
            lines = fixed_command.split('\n')
            for line in lines:
                if line and not line.startswith('```') and not line.startswith('#'):
                    fixed_command = line.strip()
                    break

        return fixed_command

    except Exception as e:
        return command  # Return original if fix fails


def get_state_update_prompt() -> PromptTemplate:
    template = """You are a Research Data Manager.
    Your task is to update the 'Project State' based on the log from a recent software run.

    **Current Project State:**
    {current_state}

    **Recent Job Summary:**
    {log_summary}

    **Instructions:**
    1. Identify if this run generated new data files or determined new parameters.
    2. **General Fields:** sequence_file.
    3. **Crystallography Fields:** original_data, data_with_free_r_flags, space_group, unit_cell, resolution, twin_law, twin_fraction.
    4. **Cryo-EM Fields:** half_map_1, half_map_2, full_map, resolution.
    5. **Models:** If a new model file (pdb/cif) was generated, ALWAYS add it to 'model_list'. Include metrics (R-work, R-free, pLDDT) if available, otherwise just list the filename.
    6. Return a JSON object containing ONLY the fields that should be UPDATED or ADDED.

    **Output Format (JSON):**
    {{
      "sequence_file": "seq.fa",
      "crystallography": {{ "space_group": "P 21 21 21" }},
      "model_list": [ {{ "file": "refine_1.pdb", "r_free": 0.22 }} ]
    }}
    """
    return PromptTemplate(template=template, input_variables=["current_state", "log_summary"])


def get_phenix_program_list() -> list:
    """
    Returns a sorted list of all available Phenix programs.
    Scans the bin directory associated with the current environment.
    """
    import sys

    # Attempt to find the bin directory
    import libtbx.load_env
    bin_dir = abs(libtbx.env.bin_path)
    if not os.path.isdir(bin_dir):
        # Fallback for some environments
        return []

    programs = []
    for filename in os.listdir(bin_dir):
        # Filter for likely Phenix/CCTBX programs
        if filename.find("development") > -1:
          continue
        if filename.startswith("phenix."):
            programs.append(filename)

    return sorted(programs)

def get_memory_file_path(db_dir):
    """
    Determines the path for the learned memory file.
    Creates the directory 'phenix_learned_info' in the parent of db_dir.
    """
    if not db_dir:
        return "phenix_learned_memory.json"

    try:
        abs_db_dir = os.path.abspath(db_dir)
        parent_dir = os.path.dirname(abs_db_dir)
        info_dir = os.path.join(parent_dir, "phenix_learned_info")

        if not os.path.exists(info_dir):
             os.makedirs(info_dir)

        return os.path.join(info_dir, "phenix_learned_memory.json")
    except Exception as e:
        print(f"Warning: Could not create/access memory dir: {e}")
        return "phenix_learned_memory.json"

def load_learned_memory(memory_file="phenix_learned_memory.json"):
    """Loads the database of learned syntax tips."""
    if os.path.exists(memory_file):
        try:
            with open(memory_file, 'r') as f:
                return json.load(f)
        except:
            pass
    return {}

def save_learned_memory(memory, memory_file="phenix_learned_memory.json"):
    """Saves the database of learned syntax tips."""
    with open(memory_file, 'w') as f:
        json.dump(memory, f, indent=2)

async def learn_from_history(run_history, llm, memory_file="phenix_learned_memory.json", logger=print):
    """
    Scans history for [Fail] -> [Success] patterns and extracts a lesson.
    Uses LLM to consolidate and deduplicate lessons.
    """
    import json
    import re

    if len(run_history) < 2:
        return

    last_run = run_history[-1]
    prev_run = run_history[-2]

    prev_summary = str(prev_run.get('summary', '')).lower()
    curr_summary = str(last_run.get('summary', '')).lower()

    prev_failed = "error" in prev_summary or "failed" in prev_summary or "exception" in prev_summary
    curr_success = "error" not in curr_summary and "failed" not in curr_summary and "exception" not in curr_summary

    prog_prev = prev_run.get('program', 'Unknown')
    prog_curr = last_run.get('program', 'Unknown')

    # Ignore empty or unknown program names
    if not prog_curr or prog_curr.lower() == 'unknown':
        return

    # Only learn if we retried the SAME program
    if prev_failed and curr_success and prog_curr == prog_prev:
        logger(f"\n[Learning] Detected a successful fix for {prog_curr}. Extracting lesson...")

        # 1. Extract the raw lesson from the logs
        extract_prompt = f"""
        Compare these two attempts to run {prog_curr}.

        FAILED ATTEMPT SUMMARY:
        {prev_run.get('summary', '')[:2000]}

        SUCCESSFUL ATTEMPT SUMMARY:
        {last_run.get('summary', '')[:2000]}

        Identify the SPECIFIC syntax or parameter change that fixed the error.
        Be concise. Format as a rule. (e.g. "Do not use 'hklin=', use bare filename.")
        """

        try:
            lesson = llm.invoke(extract_prompt).content.strip()

            # 2. Load Existing Memory
            memory = load_learned_memory(memory_file)
            if prog_curr not in memory:
                memory[prog_curr] = []

            existing_tips = memory[prog_curr]

            # 3. Consolidate: Ask LLM to merge the new tip into the list and remove duplicates
            cleanup_prompt = f"""
            You are maintaining a list of troubleshooting tips for the software '{prog_curr}'.

            Current Tips:
            {json.dumps(existing_tips)}

            New Candidate Tip:
            "{lesson}"

            INSTRUCTIONS:
            1. If the "New Candidate Tip" is redundant with Current Tips, ignore it.
            2. If it is new, add it.
            3. Merge any duplicate or very similar tips in the list into concise, unique rules.
            4. Return ONLY the final JSON list of strings (e.g. ["Tip 1", "Tip 2"]).
            """

            try:
                response = llm.invoke(cleanup_prompt).content

                # Parse JSON (handling potential Markdown wrappers)
                if "```" in response:
                    match = re.search(r'\[.*\]', response, re.DOTALL)
                    if match:
                         response = match.group(0)
                    else:
                         response = response.replace("```json", "").replace("```", "").strip()

                cleaned_list = json.loads(response)

                if isinstance(cleaned_list, list) and all(isinstance(s, str) for s in cleaned_list):
                    memory[prog_curr] = cleaned_list
                    save_learned_memory(memory, memory_file)
                    logger(f"[Learning] Consolidated tips for {prog_curr}. (Total: {len(cleaned_list)})")
                else:
                    raise ValueError("LLM did not return a valid list of strings")

            except Exception as e:
                # Fallback: Just append if consolidation fails
                logger(f"[Learning] Consolidation failed ({e}). Appending raw tip.")
                if lesson not in existing_tips:
                    memory[prog_curr].append(lesson)
                    save_learned_memory(memory, memory_file)

        except Exception as e:
            logger(f"[Learning] Failed to extract lesson: {e}")

def get_run_history(log_directory, max_history=5):
    """
    Reads JSON files from the log_directory, sorts them by time,
    and returns the last N entries.
    """
    if not os.path.exists(log_directory):
        print(f"Error: Log directory '{log_directory}' does not exist.")
        return []

    # Find all json files matching the pattern job_*.json
    search_path = os.path.join(log_directory, "job_*.json")
    files = glob.glob(search_path)

    if not files:
        return []

    # Sort files by modification time (Oldest -> Newest)
    files.sort(key=os.path.getmtime)

    # Load the data
    history = []
    for fpath in files:
        try:
            with open(fpath, 'r') as f:
                data = json.load(f)
                history.append(data)
        except Exception as e:
            print(f"Skipping corrupt file {fpath}: {e}")

    # Slice to keep only the last N
    if len(history) > max_history:
        history = history[-max_history:]

    return history



def get_keywords_as_phil_string(program: str) -> str:
    from libtbx import easy_run

    # Use the external Knowledge Base
    hint_text = ""
    if hasattr(pk, 'USAGE_HINTS') and program in pk.USAGE_HINTS:
        hint_text = f"\n{'-'*40}\n{pk.USAGE_HINTS[program]}\n{'-'*40}\n\n"

    print(f"Introspecting {program} to get valid keywords...")

    # Strategy 1: PHIL Defaults
    cmd_defaults = f"{program} --show_defaults 3"
    result = easy_run.fully_buffered(cmd_defaults)

    if result.return_code == 0 and len(result.stdout_lines) > 10 and \
       "unrecognized argument" not in "\n".join(result.stderr_lines):
        return hint_text + "\n".join(result.stdout_lines)

    # Strategy 2: Help Text
    print(f"  --show_defaults failed. Trying --help...")
    cmd_help = f"{program} --help"
    result_help = easy_run.fully_buffered(cmd_help)
    help_text = "\n".join(result_help.stdout_lines + result_help.stderr_lines)
    if len(help_text) > 100:
        return hint_text + help_text

    # Strategy 3: Naked Call
    print(f"  --help failed. Trying run with no arguments...")
    result_none = easy_run.fully_buffered(f"{program}")
    none_text = "\n".join(result_none.stdout_lines + result_none.stderr_lines)
    if len(none_text) > 100:
        return hint_text + none_text

    return f"ERROR: Could not extract keywords for {program}. It might not be in the path."

# --- Prompt Management Functions ---

def get_log_map_prompt() -> PromptTemplate:
    """Returns the prompt for summarizing a single chunk of a log file."""
    template = """
    You are a highly skilled data extraction bot. Your task is to scan a chunk of a Phenix log file and pull out only the most critical information.

    **Instructions:**
    1.  **Identify Key Steps:** Look for lines that indicate the start or end of a major computational step.
    2.  **Extract File Names:** List any input/output files mentioned. Ignore .dat files. Capture 'overall_best' files.
    3.  **Capture Metrics:** Record specific numbers (CC, resolution, R-values, etc.). Report final/refined values.
    4.  **Identify Data Type:** Note if data is X-ray or CryoEM.
    5.  **Note Warnings/Errors:** Capture critical warnings.
    6.  **Note Space Group:** Record space group and any recommendations.
    7.  **Identify Program:** Name the program if clear.
    8.  **Be Concise:** Brief, bulleted list.

    Log Chunk:
    "{text}"

    Concise, structured summary of this chunk:
    """
    return PromptTemplate.from_template(template)

def get_log_combine_prompt() -> PromptTemplate:
    """Returns the prompt for synthesizing log chunk summaries into a final report."""
    template = (
        "You are an expert research scientist. Synthesize the following summaries "
        "from a long log file into a structured report. Your report "
        "must contain the key information requested below.\n\n"
        "Final Report Requirements:\n"
        "1. **Input Files:** List key input files and data type (X-ray or cryo-EM).\n"
        "2. **Program Name:** The name of the Phenix program that was run.\n"
        "3. **Key Steps:** A high-level, bulleted list of the main steps carried out.\n"
        "4. **Key Metrics:** A concise table of the *final* key metrics (e.g., CC, R-values, resolution).\n"
        "5. **Key Space Group:** Specify the space group used and any recommendations.\n"
        "6. **Warnings/Errors:** A list of any critical warnings or errors.\n"
        "7. **Key Output Files:** A comma-separated list of the most important output filenames generated by this run (e.g. model.pdb, data.mtz).\n\n"

        "**IMPORTANT:** Be structured and clear. Do NOT provide 'Detailed descriptions' or 'Full reproductions' of tables. "
        "Focus on the most critical, final information. "
        "Do not add conversational text or offer help.\n\n"

        "Now, synthesize the following summaries:\n"
        "{context}"
    )
    return PromptTemplate(template=template, input_variables=["context"])




def get_strategic_planning_prompt() -> PromptTemplate:
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

    **ANALYSIS PROTOCOL (Follow in Order):**

    1. **TERMINATION CHECK (Stop cleanly):**
       - **SUCCESS:** If R-free < 0.25 (X-ray) or map CC is high (CryoEM) AND validation is good:
         -> **DECISION:** Output `STOP: Success`.
       - **MANUAL:** If R-factors are stalled/high and automation is failing:
         -> **DECISION:** Output `STOP: Manual Intervention`.

    2. **STATUS CHECK (Fix Crashes First):**
       - **Did the last job fail?** (e.g. "Sorry: ...", "Error", "Exception").
         -> **ACTION:** You MUST retry the *same* program (or a direct alternative) to fix the error.
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
       - **1. Data (Xtriage)** -> **2. Phasing (Phaser)** -> **3. Refinement**.
       - Do not skip to Refinement if Phasing failed or hasn't run.
       - **MODEL GENERATION:** If you have a sequence but no model, run `phenix.predict_and_build` (`stop_after_predict=True`).

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

    return PromptTemplate(template=template, input_variables=["num_runs", "history", "project_advice", "original_files", "project_state"])


def get_command_writer_prompt() -> PromptTemplate:
    template = """You are a Phenix Command-Line Expert.
    Your task is to construct a valid command line for "{program}".

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
    1. **PATH RESOLUTION (Crucial):** If a file in "Input Files" (e.g. `data.mtz`) matches a file in "Original Input Files" (e.g. `subdir/data.mtz`), you MUST use the FULL PATH from the Original list. Do NOT use bare filenames if a path is known.
    2. **Prioritize History:** If "LEARNED HISTORY" warns about a specific syntax error, you MUST follow that advice.
    3. **Follow the Reference:** Use the syntax rules found in the "Reference Documentation".
    4. **Syntax Priority:**
       - **Default:** Standard Phenix uses `key=value`.
       - **Exception:** If Reference says to use double-dashes (`--flag`), use that.
    5. **Completeness:** Ensure every defined object has an action.

    **Output:**
    Provide ONLY the command string. No markdown.
    """
    # Added "original_files" to input_variables
    return PromptTemplate(template=template, input_variables=["program", "strategy_details", "input_files", "valid_keywords", "learned_tips", "original_files"])

def get_log_analysis_prompt() -> PromptTemplate:
    """Returns the prompt for analyzing a log summary..."""
    template = (
        "You are expert in crystallography and cryo-EM."
        "You are a Phenix power-user. Your task is to analyze "
        "a program summary in the context of the provided documentation "
        "and research papers.\n\n"
        "Based on the data provided at the end, please perform the following "
        " analysis. Consider the events of the log summary in the broader "
        "context of a "
        "typical Phenix structure determination workflow as described in "
        "the documentation and papers, but do not describe this context. "
        "The analysis should include the following for the run "
        "described in the log summary:\n\n"

        "1. Evaluate whether the run described in the summary was useful. "
          "Provide a short summary of the usefulness of the run."
          "Then List reported metrics and expected values of these metrics and "
          "consider the goals of the program. Note any warnings, errors,"
          "or advisories obtained. "
          "If the results of the run indicate a low confidence solution,"
          "multiple solutions, or no solution, clearly state this observation.\n\n"

        "2. Considering whether the data are from crystallography or"
        "cryo-EM and considering the normal sequence of Phenix tool use"
        "for that type of data, suggest "
        "three concrete next steps in structure determination using Phenix "
        "or graphical tools such as Coot or Isolde that I should take, "
        "justifying each suggestion with information from the provided "
        "documentation and papers. "
        "If the results of the run indicate a low confidence solution,"
        "multiple solutions, or no solution, then instead suggest concrete"
        "steps for backtracking and figuring out what went wrong."
        "Do not include any information on CryoFit, ShelxD, Parrot, "
        "or CCP4 tools "
        "unless there is a specific question about them. "
        "If appropriate, include validation as one of your next steps. "
        "Name the tools that are to be used, along with their inputs "
        " and outputs and what they do. "
        "Do not suggest depositing the model. "
        "Do not suggest analyzing the biological relevance. \n\n"

        "3.List the inputs and briefly describe what was done."
        "Report whether the data are from crystallography (X-ray or neutron)"
        " or from cryo-EM\n\n"

        "4. List the key output files from this run, along with the values "
        "of any available metrics describing their utilities. "
        "If no metrics are available, do not provide any. \n\n"

        "**Please note: no offers of help*** Do not offer to help the "
        "user with additional analyses and do not mention that you are"
        "not to offer to help."
        "\n\n---BEGIN DATA FOR ANALYSIS---\n"
        "Documentation Context:\n{context}\n\n"
        "Log File Summary:\n{log_summary}"
    )
    return PromptTemplate(template=template, input_variables=["context", "log_summary"])


def get_keywords_prompt() -> PromptTemplate:
    """
    Returns the prompt for extracting BOTH parameters AND usage examples.
    """
    template = """You are a Phenix command-line expert.
    Your task is to provide the critical information needed to construct a valid command for "{program_name}".

    **Goal:**
    The user needs to run this program for a specific task (e.g., "{program_name} simulated annealing").
    You must extract:
    1. **Valid Parameters:** (e.g. `start_temperature=5000`)
    2. **Usage Examples:** (e.g. `phenix.refine data.mtz model.pdb annealing=True`) - THIS IS CRITICAL.

    **Instructions:**
    1. Scan the text for "Usage" or "Examples" sections.
    2. Look for patterns: how are input files specified? Do they use flags (e.g. `hklin=...`) or positional arguments?
    3. Look for the specific parameter syntax (e.g. is it `search.copies=1` or `ncopies=1`?).
    4. If broken lines appear (e.g. `param \n = value`), reconstruct them to `param=value`.

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
    """Returns the general-purpose prompt for querying the documentation RAG."""
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
    return PromptTemplate(
        template=template, input_variables=["context", "input"]
    )

# --- Core Functionality ---

def _custom_chunker(docs: List[Document], keyword_phrase: str = "List of all available keywords") -> List[Document]:
    """
    Chunks documents with UNIVERSAL context injection.
    Prepend the source filename to EVERY chunk.
    """
    final_chunks = []
    import re

    standard_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
    keyword_splitter = RecursiveCharacterTextSplitter(chunk_size=2000, chunk_overlap=200)

    for doc in docs:
        content = doc.page_content
        metadata = doc.metadata

        source_path = metadata.get('source', '')
        filename = os.path.basename(source_path)
        context_header = f"** DOCUMENT: {filename} **\n"

        match = re.search(r"list\s+of\s+all\s+available\s+keywords", content, re.IGNORECASE)
        start_index = match.start() if match else -1

        current_doc_chunks = []

        if start_index != -1:
            # Part 1: Intro
            before_content = content[:start_index]
            if before_content.strip():
                before_chunks = standard_splitter.create_documents(
                    [before_content], metadatas=[metadata])
                current_doc_chunks.extend(before_chunks)

            # Part 2: Keywords
            keyword_chunk_content = content[start_index:]
            kw_chunks = keyword_splitter.create_documents(
                [keyword_chunk_content], metadatas=[metadata])
            current_doc_chunks.extend(kw_chunks)

        else:
            # Part 1: All content
            current_doc_chunks = standard_splitter.split_documents([doc])

        # --- UNIVERSAL STAMPING ---
        for chunk in current_doc_chunks:
            chunk.page_content = context_header + chunk.page_content
            final_chunks.append(chunk)

    return final_chunks


class PhenixHTMLLoader(BaseLoader):
    """
    Custom loader that parses HTML and ensures newlines are inserted
    between tags (like list items), preventing 'one long line' issues.
    """
    def __init__(self, file_path: str):
        self.file_path = file_path

    def load(self) -> List[Document]:
        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            soup = BeautifulSoup(f, "html.parser")

        for script in soup(["script", "style"]):
            script.extract()

        text = soup.get_text(separator="\n")
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        cleaned_text = "\n".join(lines)

        metadata = {"source": self.file_path}
        return [Document(page_content=cleaned_text, metadata=metadata)]

def load_all_docs_from_folder(folder_path: str,
         excluded_dirs: list[str] = None) -> (List[Document], List[str]):
    print(f"Loading documents from {folder_path}...")

    if excluded_dirs is None:
        excluded_dirs = set()
    else:
        excluded_dirs = set(excluded_dirs)

    all_docs = []
    processed_files = []
    loaders = {
    '.txt': TextLoader,
    '.pdf': PyPDFLoader,
    '.html': PhenixHTMLLoader,
    }

    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=True):
        dirnames[:] = [d for d in dirnames if d not in excluded_dirs]
        for filename in filenames:
            if filename.startswith('.'):
                continue

            file_path = os.path.join(dirpath, filename)
            file_ext = os.path.splitext(filename)[1].lower()

            if file_ext in loaders:
                loader_cls = loaders[file_ext]
                try:
                    loader = loader_cls(file_path)
                    all_docs.extend(loader.load())
                    processed_files.append(file_path)
                    print("LOADING: %s (%s)" %(file_path, len(all_docs)))
                except Exception as e:
                    print(f"Error loading file {file_path}: {e}")

    print(f"Loaded {len(all_docs)} documents.")
    return all_docs, processed_files

def load_specific_docs(file_path_list: List[str]) -> List[Document]:
    all_docs = []
    loaders = {
        '.txt': TextLoader,
        '.pdf': PyPDFLoader,
        '.html': UnstructuredHTMLLoader
    }

    print(f"Loading {len(file_path_list)} specific files...")

    for file_path in file_path_list:
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext in loaders:
            loader_cls = loaders[file_ext]
            try:
                loader = loader_cls(file_path)
                all_docs.extend(loader.load())
                print(f"  - Successfully loaded: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"  - WARNING: Error loading file {file_path}: {e}")
        else:
            print(f"  - WARNING: No loader available for file type: {file_path}")

    print(f"Successfully loaded content from {len(all_docs)} documents.")
    return all_docs

def get_llm_and_embeddings(
    provider: str = "google",
    llm_model_name: str = None,
    embedding_model_name: str = None,
    temperature: float = 0.1,
    timeout: int = 60,
    batch_size: int = 100
):
    provider = provider.lower()
    if provider == "google":
        if llm_model_name is None:
            llm_model_name = "gemini-2.5-flash-lite"
        if embedding_model_name is None:
            embedding_model_name = "models/embedding-001"

        llm = ChatGoogleGenerativeAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=0
        )
        embeddings = GoogleGenerativeAIEmbeddings(
            model=embedding_model_name,
            timeout=timeout,
            batch_size=batch_size,
            google_api_key=os.getenv("GOOGLE_API_KEY")
        )
        print("Using Google Gemini models.")

    elif provider == "openai":
        if llm_model_name is None:
            llm_model_name = "gpt-5-nano"
        if embedding_model_name is None:
            embedding_model_name = "text-embedding-3-small"

        llm = ChatOpenAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=2
        )
        embeddings = OpenAIEmbeddings(
            model=embedding_model_name,
            chunk_size=batch_size
        )
        print("Using OpenAI models.")

    else:
        raise ValueError("Unsupported provider. Choose 'google' or 'openai'.")

    return llm, embeddings

async def get_log_info(text, llm, embeddings, timeout: int = 120,
    provider: str = 'google'):
    try:
        log_summary_info = await summarize_log_text(
            text, llm, timeout=timeout, provider = provider)
        processed_log_dict = get_processed_log_dict(
           log_summary_info.log_summary)

        if log_summary_info.error:
           return group_args(
               group_args_type = 'error', error=log_summary_info.error)

        return group_args(group_args_type = 'log summary',
          summary = log_summary_info.log_summary,
          summary_as_html = save_as_html(log_summary_info.log_summary),
          processed_log_dict = processed_log_dict,
          error = None)

    except asyncio.TimeoutError:
        error_message = "Analysis timed out, try increasing timeout."
        print(error_message)
        return group_args(group_args_type = 'error', error=error_message)

    except openai.AuthenticationError:
        error_message = "OPENAI API key is invalid"
        return group_args(group_args_type='error', error=error_message)

    except Exception as e:
        # Inspect the exception's original cause to find the specific error
        original_cause = getattr(e, '__cause__', None)

        if isinstance(original_cause, google_exceptions.PermissionDenied):
            error_message = "Google AI API key is invalid"
        elif isinstance(original_cause, google_exceptions.ResourceExhausted):
            error_message = "Google AI API quota exceeded"
        else:
            # Handle any other general exception
            msg = str(e)
            if "API_KEY_INVALID" in msg or "API_KEY_IP_ADDRESS_BLOCKED" in msg:
               error_message = "Google API key is invalid"
            elif "free_tier_requests, limit: 0" in msg:
               error_message = "Google API key is not activated"
            elif "limit: 0" in msg:
                  error_message = "Google AI API key has a zero quota"
            elif "request timed out" in msg:
                  error_message = "Summarizing timed out."
            else:
              error_message = f"An unexpected error occurred during summarization: {e}"
            print("Summarize log failed")

        print(error_message)
        return group_args(group_args_type='error', error=error_message)

async def get_program_keywords(program_name: str, llm, embeddings,
   db_dir: str = "./docs_db",
   timeout: int = 60):
  """
  Retrieves valid keywords by RUNNING the program (Introspection),
  not by querying the database.
  """
  try:
    # 1. Get Ground Truth from the executable
    raw_help_text = get_keywords_as_phil_string(program_name)

    if "ERROR:" in raw_help_text:
        return group_args(group_args_type='error', error=raw_help_text)

    # 2. Limit the size (Context Window Safety)
    max_chars = 500000
    if len(raw_help_text) > max_chars:
        print(f"Warning: Truncating parameter list from {len(raw_help_text)} to {max_chars} chars.")
        raw_help_text = raw_help_text[:max_chars] + "\n... (truncated)"

    return group_args(group_args_type = 'keywords',
      keywords = raw_help_text,
      error = None)

  except Exception as e:
      return group_args(group_args_type='error', error=str(e))

async def analyze_log_summary(log_info, llm, embeddings,
   db_dir: str = "./docs_db",
   timeout: int = 60):
  try:
    vectorstore = load_persistent_db(embeddings, db_dir = db_dir)
    analysis_prompt = get_log_analysis_prompt()
    retriever = create_reranking_retriever(vectorstore, llm,
      timeout = timeout)
    analysis_rag_chain = create_log_analysis_chain(retriever, llm,
        analysis_prompt)
    retriever_query = ("Here is a summary of the %s log file:\n\n " %(
       log_info.processed_log_dict['phenix_program']) +
       log_info.processed_log_dict['summary'] +
       "\n\nConsidering whether the input data are from crystallography "+
       "(X-ray or neutron) or from cryo-EM, and considering the "+
       "the normal procedure for structure determination "
       "in Phenix, what are the next steps that I should carry out?" +
       "Consider this question in the context of the process of structure "+
       "determination in Phenix. Focus on using Phenix tools, but include "+
       "the use of Coot or Isolde if appropriate. Name the tools that are "+
       "to be used, along with their inputs and outputs and what they do."
      )

    final_analysis = analysis_rag_chain.invoke({
        "input": retriever_query,
        "log_summary": log_info.summary,
    })

    return group_args(group_args_type = 'answer',
      analysis = final_analysis.content,
      error = None)

  # --- EXCEPTION HANDLERS ---
  except (TimeoutError, google_exceptions.DeadlineExceeded) as e:
      error_message = (
          "Network timeout with Google,"
          "You might try increasing the timeout in Preferences or in "
          f"AnalyzeLog (currently {timeout} sec)."
      )
      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except google_exceptions.ResourceExhausted as e:
      error_message = (
          "ERROR: Google AI API quota exceeded."
          f"Details: {e}"
      )
      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except CohereApiError as e:
      if hasattr(e, 'http_status') and e.http_status == 401:
            error_message = "Invalid Cohere API key.\n"
      elif str(e).find("invalid api token") > -1:
          error_message = "Cohere API key is invalid"
      else:
            error_message = f"ERROR: A Cohere API error occurred. Details: {e}"

      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)

  except Exception as e:
      error_message = "Reranking failed - try again in a couple minutes..."+str(e)
      print(error_message)
      return group_args(group_args_type = 'answer',
        analysis = None,
        error = error_message)


def _iter_batches(seq: List[Document], size: int) -> Iterable[List[Document]]:
    for i in range(0, len(seq), size):
        yield seq[i:i+size]

def create_and_persist_db(
    docs: List[Document],
    embeddings,
    db_dir: str = "./docs_db",
    *,
    add_batch_size: int = 200,
    pause_between_batches: float = 2.0,
    max_attempts: int = 6,
    max_backoff: float = 60.0
) -> Chroma:
    """
    Build a Chroma vector store with batching to avoid SQLite limits and API rate limits.
    """

    docs_chunks = _custom_chunker(docs)

    seen = set()
    unique_docs: List[Document] = []
    for d in docs_chunks:
        key = (d.page_content.strip(), tuple(sorted((d.metadata or {}).items())))
        if key not in seen:
            seen.add(key)
            unique_docs.append(d)

    client = chromadb.PersistentClient(path=db_dir)

    vectorstore = Chroma(
        client=client,
        collection_name="docs",
        embedding_function=embeddings,
    )

    for batch in _iter_batches(unique_docs, add_batch_size):
        backoff = 2.0
        for attempt in range(1, max_attempts + 1):
            try:
                vectorstore.add_documents(batch)
                break
            except Exception as e:
                msg = str(e).lower()
                if "api_key_ip_address_blocked" in msg:
                     print("Google AI API key does not allow access from this server")
                     return None
                elif "you exceeded your current quota" in msg:
                   if "limit: 0" in msg:
                     print("Google AI API key has a zero quota")
                     return None
                   else:
                     print("Google AI API quota exceeded")
                     return None
                elif "429" in msg or "rate" in msg:
                    time.sleep(backoff)
                    backoff = min(backoff * 2.0, max_backoff)
                    if attempt == max_attempts:
                        raise RuntimeError("Failed to use Google API (Rate Limit)")
                else:
                    raise RuntimeError(f"Failed to add batch: {e}")
        time.sleep(pause_between_batches)

    print(f"Database stored in: {db_dir}")
    return vectorstore

def load_persistent_db(embeddings: GoogleGenerativeAIEmbeddings,
      db_dir: str = "./docs_db") -> Chroma:
    """Loads a persisted Chroma vector store from disk."""
    if not os.path.exists(db_dir):
        raise FileNotFoundError(f"Database directory not found at '{db_dir}'.")
    return Chroma(persist_directory=db_dir, embedding_function=embeddings)

def create_reranking_retriever(vectorstore: Chroma,
         llm: ChatGoogleGenerativeAI,
         timeout: int = 60,
         top_n: int = 8):

    base_retriever = vectorstore.as_retriever(search_kwargs={"k": 20})

    cohere_client = cohere.ClientV2(
        api_key=os.getenv("COHERE_API_KEY"),
        timeout=timeout
    )

    reranker = CohereRerank(client=cohere_client, model="rerank-english-v3.0", top_n=top_n)

    compression_retriever = ContextualCompressionRetriever(
        base_compressor=reranker, base_retriever=base_retriever
    )

    return compression_retriever

def create_reranking_rag_chain(retriever, llm: ChatGoogleGenerativeAI,
        prompt: PromptTemplate):
    def format_docs(docs):
        return "\n\n".join(doc.page_content for doc in docs)

    rag_chain = (
        {"context": retriever | format_docs, "input": RunnablePassthrough()}
        | prompt
        | llm
    )
    return rag_chain

def create_keyword_lookup_chain(retriever, llm, prompt: PromptTemplate):
    inputs = {
        "context": lambda x: retriever.invoke(f"command line examples and parameters for {x['program_name']}"),
        "program_name": lambda x: x["program_name"]
    }

    chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return chain

def create_log_analysis_chain(retriever, llm: ChatGoogleGenerativeAI, prompt: PromptTemplate):
    inputs = {
        "context": lambda x: retriever.invoke(x["input"]),
        "log_summary": lambda x: x["log_summary"]
    }

    analysis_rag_chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return analysis_rag_chain

def get_processed_log_dict(log_text: str,
   summary_text_block_start = "**8",
   program_text_block_start = "**3") -> dict:
  summary = ""
  phenix_program = find_text_block(
    log_text, program_text_block_start, end_text = "**")
  if not phenix_program:
    phenix_program = find_text_block(
      log_text, program_text_block_start.replace("*","#"), end_text = "##")

  dd = {'phenix_program':phenix_program,
        'summary':summary}
  return dd

def find_text_block(log_text: str, target_text: str, end_text: str = "**"):
  text_lines = []
  started = False
  if log_text is None:
    log_text = ""
  for line in log_text.splitlines():
    if (not started) and line.strip().replace(" ","").startswith(target_text):
       started = True
    elif (started) and line.strip().replace(" ","").startswith(end_text):
       break
    if started:
      text_lines.append(line)
  return " ".join(text_lines)

def query_docs(query_text, llm=None, embeddings=None,
        db_dir: str = "./docs_db",
        timeout: int = 60,
        max_attempts: int = 5,
        use_throttling: bool = False,
        provider: str = "google"):
    """Query Phenix docs with query_text,
    with automatic retries on rate limits."""

    global _last_query_time

    if use_throttling:
        min_interval = 4
        elapsed_time = time.time() - _last_query_time
        if elapsed_time < min_interval:
            wait_time = min_interval - elapsed_time
            print(f"Throttling: Waiting {wait_time:.1f} seconds ...")
            time.sleep(wait_time)
        _last_query_time = time.time()

    if not llm:
      try:
        llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout)
      except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM with %s" %(provider))

    query_text += """Consider this question in the context of the process
        of structure determinaion in Phenix. Focus on using Phenix tools,
        but include the use of Coot or Isolde if appropriate. Name the tools
        that are to be used, along with their inputs and outputs and what
        they do."""

    vectorstore = load_persistent_db(embeddings, db_dir=db_dir)
    retriever = create_reranking_retriever(vectorstore, llm, timeout=timeout)
    prompt = get_docs_query_prompt()
    rag_chain = create_reranking_rag_chain(retriever, llm, prompt)

    backoff_time = 2.0
    for attempt in range(max_attempts):
        try:
            response = rag_chain.invoke(query_text)
            return response.content

        except openai.RateLimitError:
            if attempt < max_attempts - 1:
                print(f"OpenAI rate limit exceeded. Waiting {backoff_time:.1f} sec...")
                time.sleep(backoff_time)
                backoff_time *= 2
            else:
                print("OpenAI API rate limit exceeded after multiple retries.")
                return None
        except openai.AuthenticationError:
            print("OpenAI API key is invalid")
            return None

        except (google_exceptions.ResourceExhausted, GoogleGenerativeAIError) as e:
            error_text = str(e).lower()
            if "you exceeded your current quota" in error_text:
                if "limit: 0" in error_text:
                  print("Google AI API has a zero quota. Check plan.")
                  return None
                else:
                  print("Google AI API quota exceeded. Check plan.")
                  return None

            if "429" in error_text or "rate limit" in error_text:
                if attempt < max_attempts - 1:
                    print(f"Rate limit exceeded. Waiting for {backoff_time:.1f} sec...")
                    time.sleep(backoff_time)
                    backoff_time *= 2
                else:
                  print("Google AI API quota exceeded after retries.")
                  return None
            else:
                print(f"An unexpected Google AI error occurred: {e}")
                return None

        except Exception as e:
            print(f"An unexpected error occurred during query: {e}")
            return None

    print("Query docs failed...")
    return None


def get_chunk_size(provider: str = 'google'):
    provider = provider.lower()
    if provider == "openai":
        chunk_size = 100000
        chunk_overlap = 10000
        print("Using smaller chunk size for OpenAI.")
    elif provider == "google":
        chunk_size = 750000
        chunk_overlap = 50000
        print("Using larger chunk size for Google Gemini.")
    else:
        chunk_size = 100000
        chunk_overlap = 10000
        print(f"Warning: Unknown provider '{provider}'. Defaulting to smaller chunk size.")
    return chunk_size, chunk_overlap

def _custom_log_chunker(log_text: str,
     provider: str = "google") -> List[Document]:
    """
    Chunks a log text with a special rule for the 'Files are in the directory' section.
    """

    chunk_size, chunk_overlap = get_chunk_size(provider)

    trigger_phrase = "Files are in the directory"
    end_phrase = "Citations"

    final_chunks = []
    start_index = log_text.lower().find(trigger_phrase.lower())

    if start_index == -1:
        documents = [Document(page_content=log_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        return standard_splitter.split_documents(documents)

    before_text = log_text[:start_index].strip()
    if before_text:
        before_doc = [Document(page_content=before_text)]
        standard_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size, chunk_overlap=chunk_overlap
        )
        final_chunks.extend(standard_splitter.split_documents(before_doc))

    special_section_text = log_text[start_index:]
    end_index = special_section_text.find(end_phrase)

    if end_index != -1:
        special_chunk_content = special_section_text[:end_index].strip()
    else:
        special_chunk_content = special_section_text.strip()

    if special_chunk_content:
        final_chunks.append(Document(page_content=special_chunk_content))

    return final_chunks

async def summarize_log_text(
    text: str,
    llm: ChatGoogleGenerativeAI,
    timeout: int = 120,
    batch_size: int = 3,
    pause_between_batches: int = 1,
    use_throttling: bool = True,
    provider: str = 'google',
) -> str:
    """
    Performs a map-reduce summarization with batching to respect API rate limits.
    """
    docs = _custom_log_chunker(text, provider = provider)
    if not docs:
      return group_args(group_args_type = 'log_summary',
        log_summary = None,
        error = "Log file produced no content to summarize.")

    map_prompt = get_log_map_prompt()
    map_chain = map_prompt | llm

    all_intermediate_summaries = []

    num_batches = (len(docs) + batch_size - 1) // batch_size
    print(f"Summarizing {len(docs)} chunks in {num_batches} batches to respect rate limits...")

    for i, batch in enumerate(_iter_batches(docs, batch_size)):
        print(f"  - Processing batch {i + 1} of {num_batches}...")

        tasks = [map_chain.ainvoke({"text": doc.page_content}) for doc in batch]

        try:
            map_results = await asyncio.wait_for(asyncio.gather(*tasks), timeout=timeout)
            all_intermediate_summaries.extend([result.content for result in map_results])
        except asyncio.TimeoutError:
            print(f"Batch {i + 1} timed out after {timeout} seconds. Proceeding with partial results.")
            continue

        if use_throttling and (i < num_batches - 1):
            print(f"  - Pausing {pause_between_batches} seconds ...")
            await asyncio.sleep(pause_between_batches)

    if not all_intermediate_summaries:
      return group_args(group_args_type = 'log_summary',
        log_summary = None,
        error = "Log summarization failed to produce any results.")

    summary_docs = [Document(page_content=s) for s in all_intermediate_summaries]

    combine_prompt = get_log_combine_prompt()
    reduce_chain = create_stuff_documents_chain(llm, combine_prompt)
    final_output = reduce_chain.invoke({"context": summary_docs})
    return group_args(group_args_type = 'log_summary',
      log_summary = final_output,
      error = None)

async def extract_project_state_updates(log_summary, current_state, llm):
    """Uses LLM to identify updates to the persistent project database."""
    try:
        if not log_summary: return {}

        prompt = get_state_update_prompt()
        chain = prompt | llm

        # Serialize current state for context
        state_str = json.dumps(current_state, indent=2) if current_state else "Empty"

        response = await chain.ainvoke({
            "current_state": state_str,
            "log_summary": log_summary
        })

        # content is likely JSON string, possibly wrapped in markdown
        text = response.content.strip()
        if "```" in text:
            import re
            match = re.search(r'\{.*\}', text, re.DOTALL)
            if match: text = match.group(0)

        return json.loads(text)
    except Exception as e:
        print(f"State extraction failed: {e}")
        return {}

def extract_output_files(summary_text):
    """Parses the 'Key Output Files' section from a summary."""
    if not summary_text: return []
    files = []
    try:
        if "**Key Output Files:**" in summary_text:
            # Grab text after header
            post = summary_text.split("**Key Output Files:**")[1]
            # Stop at next header or end
            content = post.split("**")[0].strip()
            # Split by commas or newlines and clean
            raw_items = content.replace('\n', ',').split(',')
            for item in raw_items:
                clean = item.strip()
                # Basic validation: must look like a file
                if clean and '.' in clean and len(clean) > 2:
                    files.append(clean)
    except:
        pass
    return files

async def generate_next_move(
    run_history: list,
    llm,
    embeddings,
    db_dir: str = "./docs_db",
    timeout: int = 60,
    project_advice: str = None,
    original_files: list = None,
    project_state: dict = None
):
    from langchain_core.output_parsers import JsonOutputParser

    # --- LOGGING CAPTURE ---
    process_log = []
    def log(msg):
        print(msg)
        process_log.append(str(msg))

    # Get cheap llm
    try:
        cheap_llm, _ = get_llm_and_embeddings(
            provider='google',
            timeout=timeout,
        )
    except Exception as e:
        log(f"Warning: Could not load cheap LLM for learning: {e}")
        cheap_llm = llm  # Fallback to expensive if cheap fails

    # --- 1. LEARNING PHASE ---
    memory_file_path = get_memory_file_path(db_dir)
    await learn_from_history(run_history, cheap_llm, memory_file=memory_file_path, logger=log)

    log(f"thinking... (Analyzing {len(run_history)} past runs)")

    # --- BUILD COMMAND HISTORY ---
    # We collect all commands run previously to prevent duplicates
    past_commands = set()
    for run in run_history:
        nm = run.get('next_move')
        if nm and isinstance(nm, dict) and 'command' in nm:
            cmd = nm['command'].strip()
            # We normalize whitespace to catch "cmd arg" vs "cmd  arg"
            if cmd:
                norm_cmd = " ".join(cmd.split())
                past_commands.add(norm_cmd)

    # Format Inputs
    advice_text = project_advice if project_advice else "No specific advice provided. Optimize for best structure."

    if original_files:
        if isinstance(original_files, str):
            original_files_list = original_files.split()
        else:
            original_files_list = []
            for item in original_files:
                if item:
                    original_files_list.extend(str(item).split())

        orig_files_text = ", ".join(original_files_list)
    else:
        orig_files_text = "None provided."
        original_files_list = []

    history_text = ""
    for i, run in enumerate(run_history):
        history_text += f"\n--- JOB {i+1} ---\nSummary: {run.get('summary', 'N/A')}\nAnalysis: {run.get('analysis', 'N/A')}\n"

    strategy_prompt = get_strategic_planning_prompt()
    strategy_chain = strategy_prompt | llm | JsonOutputParser()

    plan = {}
    program = "Unknown"
    long_term_plan = "No plan generated."
    command = "No command generated."

    max_retries = 5
    for attempt in range(max_retries):
        try:
            # Format state for display
            state_text = json.dumps(project_state, indent=2) if project_state else "No database available."

            plan = strategy_chain.invoke({
                "num_runs": len(run_history),
                "history": history_text,
                "project_advice": advice_text,
                "original_files": orig_files_text,
                "project_state": state_text
            })


            program = plan.get('selected_program', '').strip()
            long_term_plan = plan.get('long_term_plan', 'No plan provided.')

            # --- VALIDATION CHECKS ---

            # 0. Handle Explicit Stop (Success, Manual, or Missing)
            if "STOP" in program.upper() or "MISSING" in program.upper():

                # Determine the header based on the type of stop
                if "SUCCESS" in program.upper():
                    header = "**PROJECT COMPLETE (SUCCESS):**"
                elif "MANUAL" in program.upper():
                    header = "**AUTOMATION STOPPED (MANUAL STEPS REQUIRED):**"
                else:
                    header = "**MISSING INPUT / ERROR:**"

                return group_args(
                    group_args_type='next_move',
                    command="No command generated.",
                    explanation=f"{header}\n{plan.get('reasoning', '')}",
                    program="STOP",
                    strategy=plan.get('strategy_details', ""),
                    process_log="\n".join(process_log),
                    error=None
                )

            # 1. Program Name Format
            if not program.startswith("phenix."):
                 log(f"Correction: '{program}' is not a Phenix tool. Retrying...")
                 history_text += f"\n\nSYSTEM NOTE: '{program}' is not a valid Phenix command. You must select a tool starting with 'phenix.' (e.g. phenix.refine).\n"
                 continue

            # 2. Auto-Correct Known Hallucinations
            if "alphafold" in program:
                log(f"Notice: Agent selected 'phenix.alphafold'. Auto-switching to 'phenix.predict_model'.")
                program = "phenix.predict_model"
                plan['selected_program'] = program

            if "cosym" in program:
                log(f"Notice: Agent selected 'phenix.cosym'. Auto-switching to 'phenix.xtriage'.")
                program = "phenix.xtriage"
                plan['selected_program'] = program

            if "reindex" in program:
                log(f"Notice: Agent selected 'phenix.reindex'. Auto-switching to 'phenix.xtriage'.")
                program = "phenix.xtriage"
                plan['selected_program'] = program

            # 3. Check Allow-List
            # Prefer the Knowledge Base, fall back to dynamic list
            allowed_programs = getattr(pk, 'VALID_PHENIX_PROGRAMS', None)
            if allowed_programs is None:
                 allowed_programs = get_phenix_program_list()

            if program not in allowed_programs:
                 log(f"Correction: '{program}' is not in the allowed list. Retrying...")
                 history_text += f"\n\nSYSTEM NOTE: '{program}' is not supported. Choose a standard tool.\n"
                 continue

            # 4. Robust Strategy Formatting
            strategy_details = plan.get('strategy_details', "")
            if isinstance(strategy_details, dict):
                log("Notice: Agent returned dictionary for strategy. Converting to string.")
                strategy_details = ", ".join([f"{k}={v}" for k, v in strategy_details.items()])
            elif isinstance(strategy_details, list):
                log("Notice: Agent returned list for strategy. Converting to string.")
                strategy_details = ", ".join([str(s) for s in strategy_details])
            plan['strategy_details'] = str(strategy_details)

            # 5. Placeholder Check
            strategy_upper = plan['strategy_details'].upper()
            placeholders = ["HIGHER_", "FROM_", "INSERT_", "TODO", "SELECT_", "CORRECT_", "APPROPRIATE_"]
            if any(p in strategy_upper for p in placeholders):
                 log(f"Correction: Strategy contains prohibited placeholder. Retrying...")
                 history_text += "\n\nSYSTEM NOTE: You used a placeholder. If missing values, run a diagnostic tool first.\n"
                 continue

            # 6. Anti-Looping Check
            if len(run_history) > 0:
                last_run = run_history[-1]
                last_program = last_run.get('program', '').strip()
                last_summary = last_run.get('summary', '').lower()
                last_run_failed = " error " in last_summary or "failed" in last_summary or "exception" in last_summary

                if program == last_program and not last_run_failed:
                    diagnostic_tools = ["phenix.xtriage", "phenix.mtz.dump", "phenix.reflection_file_converter", "phenix.explore_metric_symmetry"]
                    if any(tool in program for tool in diagnostic_tools):
                         log(f"Correction: Preventing diagnostic loop for '{program}'. Retrying...")
                         history_text += f"\n\nSYSTEM NOTE: You just ran '{program}' successfully. Do not run it again.\n"
                         continue

            # 7. Unit Cell Validation
            if "unit_cell" in plan['strategy_details'].lower():
                uc_val = ""
                parts = plan['strategy_details'].split(",")
                for p in parts:
                    if "unit_cell" in p:
                        uc_val = p
                        break
                if " a " in uc_val or " b " in uc_val or "alpha" in uc_val:
                     log(f"Correction: Unit cell contains variables. Retrying...")
                     history_text += "\n\nSYSTEM NOTE: 'unit_cell' must be numbers. Provide values or omit.\n"
                     continue


            # --- 8. FILE AVAILABILITY CHECK ---
            input_files_raw = plan.get('input_files', [])
            files_to_check = []
            if isinstance(input_files_raw, str):
                files_to_check = [f.strip() for f in input_files_raw.replace(',', ' ').split() if '.' in f]
            elif isinstance(input_files_raw, list):
                files_to_check = [str(f).strip() for f in input_files_raw]

            # A. Build List of Available Files
            available_files = set()

            # 1. Add Original Files (Ground Truth)
            if original_files_list:
                for f in original_files_list:
                    available_files.add(os.path.basename(f))

            # 2. Add Files from History (Using Regex)
            for run in run_history:
                outs = run.get('output_files', [])
                if not outs:
                    outs = extract_output_files(run.get('summary', ''))
                for f in outs:
                    available_files.add(os.path.basename(f))

            # --- DEBUG: Print what the Agent sees ---
            # This fulfills your request to see available files during the thought process
            log(f"DEBUG: Available Files: {list(available_files)}")

            # B. Check Files
            missing_files = []
            for f in files_to_check:
                if len(f) < 3: continue
                f_base = os.path.basename(f)

                if hasattr(pk, 'INVALID_FILENAMES') and f_base in pk.INVALID_FILENAMES:
                    missing_files.append(f)
                    continue

                if f_base not in available_files:
                    missing_files.append(f)

                if missing_files:
                 log(f"Notice: Agent tried to use missing files: {missing_files}. Retrying with valid file list...")

                 # Create a list of valid files to help the Agent correct itself
                 valid_list_str = ", ".join(sorted(list(available_files)))

                 history_text += (
                     f"\n\nSYSTEM NOTE: The file(s) {missing_files} do NOT exist. "
                     f"You cannot use them.\n"
                     f"Here is the list of VALID files you can use: [{valid_list_str}].\n"
                     "Please correct the filename in your strategy. "
                     "If the file you need is truly missing and not in this list, output 'DECISION: STOP'."
                 )
                 continue


            # --- STEP 2A: FETCH KEYWORDS (Moved Inside Loop) ---
            log(f"fetching keywords for {program}...")
            kw_result = await get_program_keywords(program, llm, embeddings, db_dir=db_dir, timeout=timeout)

            if kw_result.error or not kw_result.keywords:
                # If keywords fail, retry (Agent might have picked bad program name)
                log(f"Correction: Could not find keywords for '{program}'. Retrying...")
                history_text += f"\n\nSYSTEM NOTE: Could not find keywords/help for '{program}'. It might be invalid. Try a different tool.\n"
                continue

            # --- STEP 3A: CONSTRUCT COMMAND (Moved Inside Loop) ---
            log("constructing command...")

            memory = load_learned_memory(memory_file_path)
            relevant_tips = ""
            if program in memory:
                tips_list = memory[program]
                relevant_tips = "\n".join([f"- {tip}" for tip in tips_list])
                log(f"  -> Applying {len(tips_list)} learned tips for {program}")
            else:
                relevant_tips = "No specific history for this program."

            cmd_prompt = get_command_writer_prompt()
            cmd_chain = cmd_prompt | llm

            # Generate the command object
            command_obj = cmd_chain.invoke({
                "program": program,
                "strategy_details": plan.get('strategy_details', ""),
                "input_files": plan.get('input_files', ""),
                "valid_keywords": kw_result.keywords,
                "learned_tips": relevant_tips,
                "original_files": orig_files_text,
            })
            # Extract string content
            command = command_obj.content.strip()


            # --- 9. DUPLICATE COMMAND CHECK ---
            norm_cmd = " ".join(command.split())
            if norm_cmd in past_commands:
                log(f"Correction: Agent generated exact duplicate command. Retrying...")
                history_text += f"\n\nSYSTEM NOTE: You just generated the command:\n'{command}'\nThis command was ALREADY RUN in a previous job and failed/insufficient. You CANNOT run the exact same command again. You MUST change at least one parameter (e.g. add `ncopies`, change resolution, add `twinning=True`) OR run a different command.\n"
                continue
            # ----------------------------------------


            break

        except Exception as e:
            if attempt == max_retries - 1:
                log(f"ERROR: Strategy generation failed: {e}")
                # --- FIX: Return complete object to prevent AttributeError ---
                return group_args(
                    group_args_type='error',
                    error=f"Strategy generation failed: {e}",
                    program="CRASH",
                    explanation="Strategy generation failed after retries.",
                    strategy="",
                    command="No command generated.",
                    process_log="\n".join(process_log)
                )
                # -------------------------------------------------------------
            log(f"Strategy generation error: {e}. Retrying...")


    if not command or command == "No command generated.":
      log(f"ERROR: No valid command was generated after {max_retries} attempts")
      return group_args(
        group_args_type='error',
        error="Failed to generate valid command",
        program="UNKNOWN",
        explanation="Command generation failed",
        strategy="",
        command="No command generated.",
        process_log="\n".join(process_log)
       )


    # --- Pre-Flight Validation ---
    log("\n--- Pre-Flight Validation ---")
    log(f"Initial command: {command}")

    is_valid, error = validate_phenix_command(command)

    fix_attempts = 0
    max_fix_attempts = 3

    while not is_valid and fix_attempts < max_fix_attempts:
        fix_attempts += 1
        log(f"Validation failed (attempt {fix_attempts}/{max_fix_attempts})")
        log(f"Error: {error[:200]}...")

        # Ask LLM to fix (using the same expensive llm already loaded)
        fixed_command = fix_command_syntax(command, error, llm)

        log(f"Fixed command: {fixed_command}")

        # Validate the fix
        command = fixed_command
        is_valid, error = validate_phenix_command(command)

    if not is_valid:
        log(f" WARNING: Command still invalid after {max_fix_attempts} attempts")
        log(f" Will return anyway. Last error: {error[:200]}")
    else:
        log(f"✓ Command validation passed")

    # Preemptive fix for common issues
    if 'phenix.refine' in program:
        original_command = command
        command = add_rfree_generation_if_needed(command, program)
        if command != original_command:
            log(f" Auto-added R-free flag generation (MTZ lacks test set)")
            log(f" Updated command: {command}")

    log("-----------------------------\n")
    # --- End Validation ---

    # Now log the final decision and return
    log(f"Decision: Run {program}")
    log(f"Reasoning: {plan.get('reasoning', 'No reasoning provided')}")

    return group_args(
        group_args_type='next_move',
        command=command,
        explanation=f"PLAN: {long_term_plan}\n\nREASONING: {plan.get('reasoning', '')}",
        program=program,
        strategy=plan.get('strategy_details', ""),
        process_log="\n".join(process_log),
        error=None
    )
