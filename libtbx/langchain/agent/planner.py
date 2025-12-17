"""
Agent planning and next-move generation.

This module contains the main agent logic for:
- Deciding what program to run next
- Constructing valid commands
- Validating and fixing commands

Usage:
    from libtbx.langchain.agent import generate_next_move

    result = await generate_next_move(run_history, llm, embeddings, db_dir)
"""
from __future__ import absolute_import, division, print_function

import os
import json

from langchain_core.output_parsers import JsonOutputParser

from libtbx import group_args
from libtbx.langchain.core import get_llm_and_embeddings
from libtbx.langchain.knowledge import (
    get_phenix_program_list,
    get_keywords_as_phil_string,
    get_strategic_planning_prompt,
    get_command_writer_prompt,
)
from libtbx.langchain.validation import (
    validate_phenix_command,
    fix_command_syntax,
    add_rfree_generation_if_needed,
)
from libtbx.langchain.agent.memory import (
    get_memory_file_path,
    load_learned_memory,
    learn_from_history,
)

# Import phenix_knowledge for constants
try:
    from libtbx.langchain import phenix_knowledge as pk
except ImportError:
    pk = None


# =============================================================================
# Helper Functions
# =============================================================================

def extract_output_files(summary_text):
    """Parses the 'Key Output Files' section from a summary."""
    if not summary_text:
        return []
    files = []
    try:
        if "**Key Output Files:**" in summary_text:
            # Grab text after header
            post = summary_text.split("**Key Output Files:**")[1]
            # Stop at next header or end
            content = post.split("**")[0].strip()
            # Extract filenames (simple pattern)
            import re
            files = re.findall(r'[\w\-\.]+\.\w{2,4}', content)
    except:
        pass
    return files


async def get_program_keywords(program_name: str, llm, embeddings,
                               db_dir: str = "./docs_db",
                               timeout: int = 60):
    """
    Retrieves valid keywords by RUNNING the program (Introspection),
    not by querying the database.

    Args:
        program_name: Name of the Phenix program
        llm: Language model (not used, kept for API compatibility)
        embeddings: Embeddings model (not used, kept for API compatibility)
        db_dir: Database directory (not used, kept for API compatibility)
        timeout: Timeout in seconds

    Returns:
        group_args with 'keywords' or 'error'
    """
    try:
        # 1. Get Ground Truth from the executable
        raw_help_text = get_keywords_as_phil_string(program_name)

        if "ERROR:" in raw_help_text:
            return group_args(group_args_type='error', error=raw_help_text, keywords=None)

        # 2. Limit the size (Context Window Safety)
        max_chars = 500000
        if len(raw_help_text) > max_chars:
            print(f"Warning: Truncating parameter list from {len(raw_help_text)} to {max_chars} chars.")
            raw_help_text = raw_help_text[:max_chars] + "\n... (truncated)"

        return group_args(
            group_args_type='keywords',
            keywords=raw_help_text,
            error=None
        )

    except Exception as e:
        return group_args(group_args_type='error', error=str(e), keywords=None)


# =============================================================================
# Main Planning Function
# =============================================================================

async def generate_next_move(
    run_history: list,
    llm,
    embeddings,
    db_dir: str = "./docs_db",
    timeout: int = 60,
    project_advice: str = None,
    original_files: list = None,
    project_state: dict = None,
    file_list: list = None,
):

    """
    Generates the next move for the structure determination agent.

    This is the main agent function that:
    1. Learns from past failures/successes
    2. Analyzes project history
    3. Selects the next program to run
    4. Constructs and validates the command

    Args:
        run_history: List of past run dictionaries
        llm: Language model for planning
        embeddings: Embeddings model
        db_dir: Path to documentation database
        timeout: Timeout in seconds
        project_advice: User's goal or advice
        original_files: List of original input files
        project_state: Current project state dictionary
        file_list: List of all available files (from active_file_list.json)
    Returns:
        group_args with:
            - command: The command to run
            - explanation: Explanation of the decision
            - program: Selected program name
            - strategy: Strategy details
            - process_log: Debug log
            - error: Error message if any
    """
    # --- LOGGING CAPTURE ---
    process_log = []
    def log(msg):
        print(msg)
        process_log.append(str(msg))

    # Get cheap llm for learning
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
    past_commands = set()
    for run in run_history:
        nm = run.get('next_move')
        if nm and isinstance(nm, dict) and 'command' in nm:
            cmd = nm['command'].strip()
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

            # 0. Handle Explicit Stop
            if "STOP" in program.upper() or "MISSING" in program.upper():
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
            allowed_programs = getattr(pk, 'VALID_PHENIX_PROGRAMS', None) if pk else None
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
                raw_files = [f.strip() for f in input_files_raw.replace(',', ' ').split() if '.' in f]
                # Strip parameter prefixes like "data=", "model=", "ligand=" etc.
                for f in raw_files:
                    if '=' in f:
                        f = f.split('=', 1)[-1]  # Take everything after the '='
                    # Also strip quotes
                    f = f.strip('"\'')
                    if f:
                        files_to_check.append(f)
            elif isinstance(input_files_raw, list):
                for f in input_files_raw:
                    f_str = str(f).strip()
                    if '=' in f_str:
                        f_str = f_str.split('=', 1)[-1]
                    f_str = f_str.strip('"\'')
                    if f_str:
                        files_to_check.append(f_str)

            # A. Build List of Available Files
            available_files = set()

            # 1. Add Original Files (Ground Truth)
            if original_files_list:
                for f in original_files_list:
                    available_files.add(os.path.basename(f))

            # 1b. Add Files from file_list (active files)
            if file_list:
                for f in file_list:
                    available_files.add(os.path.basename(f))


            # 2. Add Files from History
            for run in run_history:
                outs = run.get('output_files', [])
                if not outs:
                    outs = extract_output_files(run.get('summary', ''))
                for f in outs:
                    available_files.add(os.path.basename(f))

            log(f"DEBUG: Available Files: {list(available_files)}")

            # B. Check Files
            missing_files = []
            for f in files_to_check:
                if len(f) < 3:
                    continue
                f_base = os.path.basename(f)

                if pk and hasattr(pk, 'INVALID_FILENAMES') and f_base in pk.INVALID_FILENAMES:
                    missing_files.append(f)
                    continue

                if f_base not in available_files:
                    missing_files.append(f)

            if missing_files:
                log(f"Notice: Agent tried to use missing files: {missing_files}. Retrying with valid file list...")
                valid_list_str = ", ".join(sorted(list(available_files)))
                history_text += (
                    f"\n\nSYSTEM NOTE: The file(s) {missing_files} do NOT exist. "
                    f"You cannot use them.\n"
                    f"Here is the list of VALID files you can use: [{valid_list_str}].\n"
                    "Please correct the filename in your strategy. "
                    "If the file you need is truly missing and not in this list, output 'DECISION: STOP'."
                )
                continue

            # --- STEP 2A: FETCH KEYWORDS ---
            log(f"fetching keywords for {program}...")
            kw_result = await get_program_keywords(program, llm, embeddings, db_dir=db_dir, timeout=timeout)

            if kw_result.error or not kw_result.keywords:
                log(f"Correction: Could not find keywords for '{program}'. Retrying...")
                history_text += f"\n\nSYSTEM NOTE: Could not find keywords/help for '{program}'. It might be invalid. Try a different tool.\n"
                continue

            # --- STEP 3A: CONSTRUCT COMMAND ---
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

            command_obj = cmd_chain.invoke({
                "program": program,
                "strategy_details": plan.get('strategy_details', ""),
                "input_files": plan.get('input_files', ""),
                "valid_keywords": kw_result.keywords,
                "learned_tips": relevant_tips,
                "original_files": orig_files_text,
                "project_advice": advice_text,
            })

            command = command_obj.content.strip()

            # --- 9. DUPLICATE COMMAND CHECK ---
            norm_cmd = " ".join(command.split())
            if norm_cmd in past_commands:
                log(f"Correction: Agent generated exact duplicate command. Retrying...")
                history_text += f"\n\nSYSTEM NOTE: You just generated the command:\n'{command}'\nThis command was ALREADY RUN in a previous job and failed/insufficient. You CANNOT run the exact same command again. You MUST change at least one parameter (e.g. add `ncopies`, change resolution, add `twinning=True`) OR run a different command.\n"
                continue

            break

        except Exception as e:
            if attempt == max_retries - 1:
                log(f"ERROR: Strategy generation failed: {e}")
                return group_args(
                    group_args_type='error',
                    error=f"Strategy generation failed: {e}",
                    program="CRASH",
                    explanation="Strategy generation failed after retries.",
                    strategy="",
                    command="No command generated.",
                    process_log="\n".join(process_log)
                )
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

        fixed_command = fix_command_syntax(command, error, llm)
        log(f"Fixed command: {fixed_command}")

        command = fixed_command
        is_valid, error = validate_phenix_command(command)

    # Mechanical fallback: strip unrecognized parameters if LLM fix failed
    if not is_valid and ("not recognized" in error.lower() or
                          "unknown" in error.lower() or
                          "ambiguous" in error.lower()):
        import re
        # Try to extract and remove the bad parameter
        patterns = [
            r"parameter.*?:\s*(\w+)\s*=",
            r"definition:\s*(\w+)\s*=",
        ]
        for pattern in patterns:
            match = re.search(pattern, error, re.IGNORECASE)
            if match:
                bad_param = match.group(1)
                stripped = re.sub(rf'\s*\b{bad_param}\s*=\s*\S+', '', command)
                stripped = ' '.join(stripped.split())
                if stripped != command:
                    log(f" Stripping invalid parameter '{bad_param}'...")
                    is_valid, error = validate_phenix_command(stripped)
                    if is_valid:
                        command = stripped
                        log(f" ✓ Parameter strip succeeded")
                    break

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

