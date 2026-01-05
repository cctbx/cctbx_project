"""
Agent memory and learning functions.

This module handles:
- Persistent storage of learned syntax tips
- Learning from success/failure patterns
- Run history management

Usage:
    from libtbx.langchain.agent import load_learned_memory, learn_from_history

    memory = load_learned_memory(memory_file)
    await learn_from_history(run_history, llm, memory_file)
"""
from __future__ import absolute_import, division, print_function
from langchain_core.prompts import PromptTemplate
import os
import json
import glob
import re


# =============================================================================
# Memory File Management
# =============================================================================

def get_memory_file_path(db_dir):
    """
    Determines the path for the learned memory file.
    Creates the directory 'phenix_learned_info' in the parent of db_dir.

    Args:
        db_dir: Path to the database directory

    Returns:
        str: Path to the memory file
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
    """
    Loads the database of learned syntax tips.

    Args:
        memory_file: Path to the memory file

    Returns:
        dict: Program name -> list of tips
    """
    if os.path.exists(memory_file):
        try:
            with open(memory_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            pass
    return {}


def save_learned_memory(memory, memory_file="phenix_learned_memory.json"):
    """
    Saves the database of learned syntax tips.

    Args:
        memory: Dictionary of program -> tips
        memory_file: Path to save to
    """
    with open(memory_file, 'w') as f:
        json.dump(memory, f, indent=2)


# =============================================================================
# Learning from History
# =============================================================================

async def learn_from_history(run_history, llm, memory_file="phenix_learned_memory.json", logger=print):
    """
    Scans history for [Fail] -> [Success] patterns and extracts a lesson.
    Uses LLM to consolidate and deduplicate lessons.

    Args:
        run_history: List of run dictionaries
        llm: Language model for extracting lessons
        memory_file: Path to memory file
        logger: Logging function (default: print)
    """
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

        # 1. Extract the raw lesson from the logs using the Safe Pattern
        template = """Compare these two attempts to run {prog_curr}.

I will provide the failed attempt summary in <failed_attempt> tags,
and the successful attempt summary in <successful_attempt> tags.

FAILED ATTEMPT SUMMARY:
<failed_attempt>
{failed_summary}
</failed_attempt>

SUCCESSFUL ATTEMPT SUMMARY:
<successful_attempt>
{success_summary}
</successful_attempt>

Identify the SPECIFIC syntax or parameter change that fixed the error.
Be concise. Format as a rule. (e.g. "Do not use 'hklin=', use bare filename.")
"""

        # Create the safe prompt object
        prompt_template = PromptTemplate.from_template(template)

        # Generate the string safely
        extract_prompt = prompt_template.format(
            prog_curr=prog_curr,
            failed_summary=prev_run.get('summary', '')[:2000],
            success_summary=last_run.get('summary', '')[:2000]
        )


        try:
            lesson = llm.invoke(extract_prompt).content.strip()

            # 2. Load Existing Memory
            memory = load_learned_memory(memory_file)
            if prog_curr not in memory:
                memory[prog_curr] = []

            existing_tips = memory[prog_curr]

            # 3. Consolidate: Ask LLM to merge the new tip using the Safe Pattern
            template = """You are maintaining a list of troubleshooting tips for the software '{prog_curr}'.

I will provide the current list of tips in <current_tips> tags,
and a new candidate tip in <new_candidate_tip> tags.

Current Tips:
<current_tips>
{current_tips}
</current_tips>

New Candidate Tip:
<new_candidate_tip>
{new_tip}
</new_candidate_tip>

INSTRUCTIONS:
1. If the "New Candidate Tip" is redundant with Current Tips, ignore it.
2. If it is new, add it.
3. Merge any duplicate or very similar tips in the list into concise, unique rules.
4. Return ONLY the final JSON list of strings (e.g. ["Tip 1", "Tip 2"]).
"""

            # Create the safe prompt object
            prompt_template = PromptTemplate.from_template(template)

            # Generate the string safely
            # Note: json.dumps happens HERE, not inside the template string
            cleanup_prompt = prompt_template.format(
                prog_curr=prog_curr,
                current_tips=json.dumps(existing_tips),
                new_tip=lesson
            )

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


# =============================================================================
# Run History Management
# =============================================================================

def get_run_history(log_directory, max_history=5):
    """
    Reads JSON files from the log_directory, sorts them by time,
    and returns the last N entries.

    Args:
        log_directory: Path to directory containing job_*.json files
        max_history: Maximum number of entries to return

    Returns:
        list: List of run dictionaries, oldest to newest
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
