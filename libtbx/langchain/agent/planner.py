"""
Agent planning and next-move generation.

This module contains the main agent logic for:
- Deciding what program to run next
- Constructing valid commands
- Validating and fixing commands

Usage:
    from libtbx.langchain.agent.planner import generate_next_move

    result = await generate_next_move(run_history, llm, embeddings, db_dir = xx, cheap_llm = xxx)
"""
from __future__ import absolute_import, division, print_function

import re
import os
import json

from langchain_core.output_parsers import JsonOutputParser

from libtbx import group_args
from libtbx.langchain.knowledge.phenix_programs import (
    get_phenix_program_list,
    get_keywords_as_phil_string,
)
from libtbx.langchain.knowledge.prompts import (
    get_strategic_planning_prompt,
    get_command_writer_prompt,
)
from libtbx.langchain.validation.core_validator import (
    validate_phenix_command,
    fix_command_syntax,
)
from libtbx.langchain.validation.mtz_utils import add_rfree_generation_if_needed
from libtbx.langchain.agent.memory import (
    get_memory_file_path,
    load_learned_memory,
    learn_from_history,
)

# Import phenix_knowledge for constants
try:
    from libtbx.langchain import phenix_knowledge as pk
except Exception as e:
    pk = None

# =============================================================================
# Helper Functions
# =============================================================================
_PARAMETER_FIXES = None

def construct_command_mechanically(program, plan, original_files_list, required_params, log):
    """
    Construct commands mechanically for key programs.

    Args:
        program: The program name (e.g., 'phenix.refine')
        plan: Plan dict from planner (contains 'input_files', 'strategy_details')
        original_files_list: List of original files provided by user
        required_params: Dict of required parameters (e.g., {'model': 'model=xxx.pdb'})
        log: Logging function

    Returns:
        str: Constructed command, or None if no mechanical construction available
    """
    command = None

    input_files_list = plan.get('input_files', '').split()
    strategy_details = plan.get('strategy_details', '').lower()

    def extract_filename(param_value):
        """Extract just the filename from 'key=filename' or 'filename'."""
        if not param_value:
            return None
        if '=' in param_value:
            return param_value.split('=', 1)[1]
        return param_value

    if program == "phenix.refine":
        # Find model (prefer *_with_ligand.pdb, then *_refine_*.pdb, then any .pdb)
        model_file = None
        data_file = None

        for f in input_files_list:
            if f.endswith('.pdb'):
                if 'with_ligand' in f or 'complex' in f:
                    model_file = f  # Highest priority - combined model
                elif model_file is None:
                    model_file = f
                elif '_refine_' in f and 'with_ligand' not in (model_file or ''):
                    model_file = f  # Prefer refined over unrefined
            elif f.endswith('.mtz'):
                data_file = f

        # Override with required_params if set (extract just filename for positional args)
        if required_params.get('model'):
            model_file = extract_filename(required_params['model'])

        if model_file and data_file:
            command = f"phenix.refine {model_file} {data_file} xray_data.r_free_flags.generate=True"
            log(f"DEBUG: Mechanically constructed refine command: {command}")

    elif program == "phenix.ligandfit":
        # required_params already has 'model': 'model=filename.pdb' format
        model_param = required_params.get('model', '')
        data_param = required_params.get('data', '')
        ligand_param = required_params.get('ligand', '')

        # Fallback: find files from input_files_list (need to add key= prefix)
        if not model_param:
            for f in input_files_list:
                if '_refine_' in f and f.endswith('.pdb') and 'ligand' not in f.lower():
                    model_param = f"model={f}"
                    break
        if not data_param:
            for f in input_files_list:
                if '_refine_' in f and f.endswith('.mtz'):
                    data_param = f"data={f}"
                    break
        if not ligand_param:
            for f in input_files_list:
                if f.endswith('.pdb') and ('lig' in f.lower() or 'ligand' in f.lower()):
                    if '_refine_' not in f:
                        ligand_param = f"ligand={f}"
                        break

        if model_param and data_param and ligand_param:
            command = f"phenix.ligandfit {model_param} {data_param} {ligand_param} file_info.input_labels=\"2FOFCWT PH2FOFCWT\" general.nproc=4"
            log(f"DEBUG: Mechanically constructed ligandfit command: {command}")

    elif program == "phenix.pdbtools":
        # Look for protein model and ligand to combine
        protein_file = None
        ligand_file = None

        for f in input_files_list:
            if f.endswith('.pdb'):
                if 'ligand_fit' in f.lower():
                    ligand_file = f
                elif f.lower().startswith('lig') and '_refine_' not in f:
                    ligand_file = f
                elif '_refine_' in f and 'with_ligand' not in f:
                    protein_file = f
                elif protein_file is None and 'with_ligand' not in f and 'ligand_fit' not in f.lower():
                    protein_file = f

        # Also search original_files for ligand if not found
        if not ligand_file:
            for f in original_files_list:
                if f.endswith('.pdb') and ('lig' in f.lower() or 'ligand' in f.lower()):
                    if '_refine_' not in f and 'with_ligand' not in f:
                        ligand_file = f
                        break

        if protein_file and ligand_file:
            # Generate output name
            base = protein_file.replace('.pdb', '')
            output_name = f"{base}_with_ligand.pdb"
            command = f"phenix.pdbtools {protein_file} {ligand_file} output.file_name={output_name}"
            log(f"DEBUG: Mechanically constructed pdbtools command: {command}")

    elif program == "phenix.phaser":
        model_file = None
        data_file = None
        seq_file = None

        for f in input_files_list:
            if f.endswith('.pdb'):
                model_file = f
            elif f.endswith('.mtz'):
                data_file = f
            elif f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                seq_file = f

        # Fallback to original_files
        if not data_file:
            for f in original_files_list:
                if f.endswith('.mtz'):
                    data_file = f
                    break
        if not seq_file:
            for f in original_files_list:
                if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                    seq_file = f
                    break

        if model_file and data_file:
            command = f"phenix.phaser {data_file} {model_file}"
            if seq_file:
                command += f" {seq_file}"
            command += " phaser.mode=MR_AUTO"
            log(f"DEBUG: Mechanically constructed phaser command: {command}")

    elif program == "phenix.xtriage":
        data_file = None
        seq_file = None

        for f in input_files_list:
            if f.endswith('.mtz'):
                data_file = f
            elif f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                seq_file = f

        # Fallback to original_files
        if not data_file:
            for f in original_files_list:
                if f.endswith('.mtz'):
                    data_file = f
                    break
        if not seq_file:
            for f in original_files_list:
                if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                    seq_file = f
                    break

        if data_file:
            command = f"phenix.xtriage {data_file}"
            if seq_file:
                command += f" scaling.input.asu_contents.sequence_file={seq_file}"
            log(f"DEBUG: Mechanically constructed xtriage command: {command}")

    elif program == "phenix.predict_and_build":
        data_file = None
        seq_file = None

        for f in input_files_list:
            if f.endswith('.mtz'):
                data_file = f
            elif f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                seq_file = f

        # Fallback to original_files
        if not data_file:
            for f in original_files_list:
                if f.endswith('.mtz'):
                    data_file = f
                    break
        if not seq_file:
            for f in original_files_list:
                if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                    seq_file = f
                    break

        if seq_file:
            command = f"phenix.predict_and_build input_files.seq_file={seq_file}"
            if data_file:
                command += f" input_files.xray_data_file={data_file}"
            # Check if stop_after_predict was in strategy
            if 'stop_after_predict' in strategy_details:
                command += " predict_and_build.stop_after_predict=True"
            log(f"DEBUG: Mechanically constructed predict_and_build command: {command}")

    elif program == "phenix.mtriage":
        map_file = None
        half_map_1 = None
        half_map_2 = None

        for f in input_files_list:
            if f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                if 'half' in f.lower() or '_1' in f or '_2' in f:
                    if half_map_1 is None:
                        half_map_1 = f
                    else:
                        half_map_2 = f
                else:
                    map_file = f

        # Fallback to original_files
        if not map_file and not half_map_1:
            for f in original_files_list:
                if f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                    map_file = f
                    break

        if half_map_1 and half_map_2:
            command = f"phenix.mtriage half_map={half_map_1} half_map={half_map_2}"
            if map_file:
                command = f"phenix.mtriage map={map_file} half_map={half_map_1} half_map={half_map_2}"
            log(f"DEBUG: Mechanically constructed mtriage command: {command}")
        elif map_file:
            command = f"phenix.mtriage {map_file}"
            log(f"DEBUG: Mechanically constructed mtriage command: {command}")

    elif program == "phenix.process_predicted_model":
        model_file = None

        for f in input_files_list:
            if f.endswith('.pdb') and 'predicted' in f.lower():
                model_file = f
                break

        # Fallback: any pdb file
        if not model_file:
            for f in input_files_list:
                if f.endswith('.pdb'):
                    model_file = f
                    break

        if model_file:
            command = f"phenix.process_predicted_model {model_file}"
            log(f"DEBUG: Mechanically constructed process_predicted_model command: {command}")

    elif program == "phenix.dock_in_map":
        model_file = None
        map_file = None
        seq_file = None

        for f in input_files_list:
            if f.endswith('.pdb'):
                model_file = f
            elif f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                map_file = f
            elif f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                seq_file = f

        # Fallback to original_files
        if not map_file:
            for f in original_files_list:
                if f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                    map_file = f
                    break
        if not seq_file:
            for f in original_files_list:
                if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.seq'):
                    seq_file = f
                    break

        if model_file and map_file:
            command = f"phenix.dock_in_map {model_file} {map_file}"
            if seq_file:
                command += f" {seq_file}"
            log(f"DEBUG: Mechanically constructed dock_in_map command: {command}")

    elif program == "phenix.real_space_refine":
        model_file = None
        map_file = None

        for f in input_files_list:
            if f.endswith('.pdb'):
                model_file = f
            elif f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                map_file = f

        # Fallback to original_files
        if not map_file:
            for f in original_files_list:
                if f.endswith('.ccp4') or f.endswith('.mrc') or f.endswith('.map'):
                    map_file = f
                    break

        if model_file and map_file:
            command = f"phenix.real_space_refine {model_file} {map_file}"
            log(f"DEBUG: Mechanically constructed real_space_refine command: {command}")

    return command


def get_required_params(program, input_files, original_files):
    """
    Pre-compute required parameters for specific programs.
    Returns a dict of param_name: value that MUST be included in the command.
    """
    required = {}

    # Combine all available files
    all_files = input_files.split() if input_files else []
    if original_files:
        all_files.extend(original_files.split())

    if program == "phenix.ligandfit":
        for f in all_files:
            # Find refined model (prefer *_refine_*.pdb)
            if "_refine_" in f and f.endswith(".pdb") and "ligand" not in f.lower():
                if "model" not in required:
                    required["model"] = f
            # Find refined MTZ
            if "_refine_" in f and f.endswith(".mtz"):
                if "data" not in required:
                    required["data"] = f
            # Find ligand file
            if f.endswith(".pdb") and ("lig" in f.lower() or "ligand" in f.lower()):
                if "_refine_" not in f:
                    required["ligand"] = f

    elif program == "phenix.refine":
        # After ligand combination, use the combined model
        for f in all_files:
            if "with_ligand" in f and f.endswith(".pdb"):
                required["model"] = f
                break
            if "complex" in f and f.endswith(".pdb"):
                required["model"] = f
                break

    return required

# Load parameter fixes from JSON
def load_parameter_fixes():
    """Load program-specific parameter fixes from JSON."""
    fixes = {}

    # JSON file is in same directory as this script
    path = os.path.join(os.path.dirname(__file__), 'parameter_fixes.json')

    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                fixes = json.load(f)
            print(f"Loaded parameter fixes from {path}")
        except Exception as e:
            print(f"Warning: Could not load parameter fixes from {path}: {e}")
    else:
        print(f"Warning: parameter_fixes.json not found at {path}")

    return fixes


def get_parameter_fixes():
    """Get parameter fixes, loading from file if needed."""
    global _PARAMETER_FIXES
    if _PARAMETER_FIXES is None:
        _PARAMETER_FIXES = load_parameter_fixes()
    return _PARAMETER_FIXES

def fix_program_parameters(command, program):
    """
    Fix common parameter mistakes for a specific program.

    Uses a JSON config that maps wrong parameter names to correct ones.

    Args:
        command: The generated command string
        program: The program name (e.g., 'phenix.ligandfit')

    Returns:
        Fixed command string
    """
    fixes = get_parameter_fixes()

    if program not in fixes:
        return command

    program_fixes = fixes[program]

    for wrong_param, right_param in program_fixes.items():
        # Skip comment keys
        if wrong_param.startswith('_'):
            continue

        # Skip if right_param is already in the command (avoid double-fixing)
        if right_param and right_param + '=' in command:
            continue

        # Pattern to match parameter=value (handles quoted values too)
        # Use word boundary and negative lookbehind to avoid matching if already prefixed
        # e.g., don't match "input_labels" if "file_info.input_labels" exists
        pattern = rf'(?<![.\w]){re.escape(wrong_param)}=("(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'|\S+)'

        if right_param == "":
            # Remove the parameter name but keep the value
            command = re.sub(pattern, r'\1', command)
        elif right_param is None:
            # Remove entirely
            command = re.sub(pattern, '', command)
        else:
            # Replace with correct parameter name
            command = re.sub(pattern, rf'{right_param}=\1', command)

    # Clean up any double spaces
    command = ' '.join(command.split())

    return command

def extract_clean_command(raw_output):
    """
    Extract the actual Phenix command from model output that may contain
    thinking tokens, chain-of-thought, or other garbage.
    """

    # Remove <think>...</think> blocks
    raw_output = re.sub(r'<think>.*?</think>', '', raw_output, flags=re.DOTALL)

    # Remove </think> if present without opening tag
    raw_output = raw_output.replace('</think>', '')

    # Look for a line that starts with phenix.
    lines = raw_output.strip().split('\n')
    for line in lines:
        line = line.strip()
        if line.startswith('phenix.'):
            return line

    # If no phenix. command found, try to find any line that looks like a command
    for line in lines:
        line = line.strip()
        # Skip empty lines and lines that look like natural language
        if not line:
            continue
        if len(line.split()) > 20:  # Too long to be a command
            continue
        if line[0].islower() and ' ' in line and not line.startswith('phenix'):
            continue  # Looks like natural language
        # Check if it starts with a word that could be a program name
        first_word = line.split()[0] if line.split() else ''
        if first_word and '.' not in first_word and not first_word.startswith('phenix'):
            continue
        return line

    # Fallback: return first non-empty line (will likely fail, but at least it's clean)
    for line in lines:
        if line.strip():
            return line.strip()

    return raw_output.strip()


def get_relative_path(filepath, working_dir=None):
    """
    Extract a usable relative path from a filepath.

    Since the server may have a different working directory than the client,
    we look for known subdirectory patterns to preserve the path structure
    needed for Phenix commands.

    Args:
        filepath: The file path (absolute or relative)
        working_dir: The working directory (optional, may not match client's)

    Returns:
        Relative path string suitable for use in Phenix commands
    """
    if not filepath:
        return filepath

    # If already a simple relative path (no absolute indicators), use as-is
    if not os.path.isabs(filepath) and not filepath.startswith('/'):
        return filepath

    # Normalize path separators
    filepath = filepath.replace('\\', '/')

    # Known output subdirectory patterns from Phenix programs
    # These patterns identify directories that contain output files
    # Order matters - check more specific patterns first
    known_subdir_patterns = [
        '_CarryOn/',           # e.g., PredictAndBuild_0_CarryOn/
        '_CarryOn',            # At end of path
        'local_prediction_',   # e.g., PredictAndBuild_0/local_prediction_1_0/
        'PredictAndBuild_',    # e.g., PredictAndBuild_0/
        'PHASER_',
        'Refine_',
        'AutoBuild_',
        'LigandFit_',
        'DockAndRebuild_',
    ]

    # First, check if path contains a _CarryOn directory - this is the most specific
    if '_CarryOn/' in filepath or filepath.endswith('_CarryOn'):
        # Find the parent directory of _CarryOn (e.g., PredictAndBuild_0_CarryOn)
        # We want to return "PredictAndBuild_0_CarryOn/filename.pdb"
        parts = filepath.split('/')
        for i, part in enumerate(parts):
            if '_CarryOn' in part:
                # Return from this directory onwards
                return '/'.join(parts[i:])

    # Check for local_prediction subdirectories
    if 'local_prediction_' in filepath:
        parts = filepath.split('/')
        for i, part in enumerate(parts):
            if part.startswith('PredictAndBuild_') or part.startswith('local_prediction_'):
                return '/'.join(parts[i:])

    # Check for other known patterns
    for pattern in known_subdir_patterns:
        if pattern in filepath:
            parts = filepath.split('/')
            for i, part in enumerate(parts):
                if pattern.rstrip('/') in part:
                    return '/'.join(parts[i:])

    # No known subdirectory found - return basename
    return os.path.basename(filepath)

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
            # Extract filenames (pattern includes paths with /)
            files = re.findall(r'[\w\-\.\/]+\.\w{2,4}', content)
    except Exception as e:
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
    cheap_llm = None,
    command_llm = None,
    db_dir: str = None,
    provider: str = None,
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
            # IMPORTANT: The server cannot validate files on disk.
            # We trust the file_list from the client (already validated) and
            # original_files. We do NOT trust files extracted from summaries
            # because the LLM may have hallucinated them.

            available_files = set()
            # Also keep a mapping from basename to full relative path
            basename_to_relpath = {}

            # Known hallucinated file patterns to reject
            hallucinated_patterns = [
                "PHENIX.1.pdb",  # Common hallucination
                "HKLOUT",  # Placeholder name, not a real file
                "file1.", "file2.",  # Generic hallucinations
                "/path/to/",  # Path placeholder
                "/nonexistent/",  # Obvious placeholder
            ]

            def is_likely_hallucinated(filename):
                """Check if filename looks like a hallucination."""
                if not filename:
                    return True
                basename = os.path.basename(filename)
                # Check against known patterns
                for pattern in hallucinated_patterns:
                    if pattern in basename or pattern in filename:
                        return True
                # Generic placeholder check
                if basename.startswith("file") and basename[4:5].isdigit():
                    return True
                return False

            # 1. Add Original Files (Ground Truth from client)
            if original_files_list:
                for f in original_files_list:
                    rel_path = get_relative_path(f)
                    basename = os.path.basename(rel_path)

                    if not is_likely_hallucinated(basename):
                        available_files.add(rel_path)
                        # Map basename to relative path for lookup
                        if basename not in basename_to_relpath:
                            basename_to_relpath[basename] = rel_path

            # 2. Add Files from file_list (validated by client)
            # This is the AUTHORITATIVE source - client validated these exist
            if file_list:
                for f in file_list:
                    rel_path = get_relative_path(f)
                    basename = os.path.basename(rel_path)

                    if not is_likely_hallucinated(basename):
                        available_files.add(rel_path)
                        # Map basename to relative path for lookup
                        # Prefer paths with subdirectories over plain basenames
                        if basename not in basename_to_relpath or '/' in rel_path:
                            basename_to_relpath[basename] = rel_path

            # 3. Add Files from History - BE VERY SKEPTICAL
            # Only add if the file was mentioned as a RESULT of a SUCCESSFUL run
            # AND matches expected output patterns for that program
            expected_output_patterns = {
                'phenix.phaser': ['PHASER.', '.sol', '_MR.'],
                'phenix.refine': ['_refine_', '_refined', '_001.pdb', '_001.mtz'],
                'phenix.autobuild': ['_autobuild_', '_rebuilt'],
                'phenix.ligandfit': ['ligand_fit_', '_lig.pdb'],
                'phenix.predict_and_build': ['_rebuilt', '_predicted', '_overall_best'],
                'phenix.dock_and_rebuild': ['_rebuilt', '_docked'],
            }

            for run in run_history:
                # Only trust output files from successful runs
                result = run.get('result', '')
                if not (isinstance(result, str) and result.startswith('SUCCESS')):
                    continue

                run_program = run.get('program', '')
                outs = run.get('output_files', [])
                if not outs:
                    outs = extract_output_files(run.get('summary', ''))

                for f in outs:
                    # Preserve path structure if present
                    rel_path = get_relative_path(f)
                    basename = os.path.basename(rel_path)

                    # Skip obvious hallucinations
                    if is_likely_hallucinated(basename):
                        log(f"WARNING: Rejecting likely hallucinated file: {basename}")
                        continue

                    # If we have expected patterns for this program, validate
                    patterns = expected_output_patterns.get(run_program, [])
                    if patterns:
                        matches_pattern = any(p in basename for p in patterns)
                        if not matches_pattern:
                            log(f"WARNING: File {basename} doesn't match expected output for {run_program}")
                            continue

                    # If already in file_list (validated by client), use that path
                    if file_list:
                        for fl in file_list:
                            if os.path.basename(fl) == basename:
                                rel_path = get_relative_path(fl)
                                available_files.add(rel_path)
                                if basename not in basename_to_relpath or '/' in rel_path:
                                    basename_to_relpath[basename] = rel_path
                                break
                        else:
                            # Not in file_list, be more skeptical
                            if basename.endswith(('.pdb', '.mtz', '.cif')):
                                cycle_num = run.get('job_id', '0')
                                try:
                                    if int(str(cycle_num).replace('job_', '')) > 1:
                                        available_files.add(rel_path)
                                        if basename not in basename_to_relpath:
                                            basename_to_relpath[basename] = rel_path
                                except Exception as e:
                                    pass
                    else:
                        # No file_list, use what we have
                        if basename.endswith(('.pdb', '.mtz', '.cif')):
                            available_files.add(rel_path)
                            if basename not in basename_to_relpath:
                                basename_to_relpath[basename] = rel_path

            log(f"DEBUG: Available Files (filtered): {sorted(list(available_files))}")
            log(f"DEBUG: Basename to path mapping: {basename_to_relpath}")

            # B. Check Files
            missing_files = []
            for f in files_to_check:
                if len(f) < 3:
                    continue
                f_base = os.path.basename(f)

                if pk and hasattr(pk, 'INVALID_FILENAMES') and f_base in pk.INVALID_FILENAMES:
                    missing_files.append(f)
                    continue

                # Check if file exists in available_files (either as full path or basename)
                file_found = False
                if f in available_files:
                    file_found = True
                elif f_base in basename_to_relpath:
                    file_found = True
                else:
                    # Check if any available file ends with this path
                    for avail in available_files:
                        if avail.endswith(f) or os.path.basename(avail) == f_base:
                            file_found = True
                            break

                if not file_found:
                    missing_files.append(f)

            if missing_files:
                log(f"Notice: Agent tried to use missing files: {missing_files}. Retrying with valid file list...")
                valid_list_str = ", ".join(sorted(list(available_files)))
                history_text += (
                    f"\n\nSYSTEM NOTE: The file(s) {missing_files} do NOT exist. "
                    f"You cannot use them.\n"
                    f"Here is the list of VALID files you can use: [{valid_list_str}].\n"
                    "Please correct the filename in your strategy. Use the FULL PATH shown above including any subdirectories. "
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
                log(f"  -> No learned tips applied for {program}")


            cmd_prompt = get_command_writer_prompt()
            cmd_chain = cmd_prompt | (command_llm if command_llm else llm)

            # Pass available files with full paths to command writer
            available_files_text = "\n".join(sorted(list(available_files)))

            # Get required parameters for this program
            required_params = get_required_params(
                program,
                plan.get('input_files', ""),
                orig_files_text
            )
            required_params_text = ""
            if required_params:
                required_params_text = "\n".join([f"{k}={v}" for k, v in required_params.items()])

            log(f"DEBUG program to run: {program}")
            log(f"DEBUG Strategy details: {plan.get('strategy_details', 'No details provided')}")
            log(f"DEBUG input files: {plan.get('input_files', 'No files provided')}")
            log(f"DEBUG original files: {orig_files_text}")
            log(f"DEBUG tips: {relevant_tips}")
            log(f"DEBUG advice: {advice_text}")
            log(f"DEBUG required_params_text: {required_params_text}")

            command_obj = cmd_chain.invoke({
                "program": program,
                "strategy_details": plan.get('strategy_details', ""),
                "input_files": plan.get('input_files', ""),
                "valid_keywords": kw_result.keywords,
                "learned_tips": relevant_tips,
                "original_files": f"{orig_files_text}\n\nALL AVAILABLE FILES (use exact paths including subdirectories):\n{available_files_text}",
                "project_advice": advice_text,
                "required_params": required_params_text,
            })

            command = command_obj.content.strip()
            log(f"DEBUG: Original command: {command}")
            # Strip thinking tokens and extract just the command
            command = extract_clean_command(command)

            # Save the agent's original command
            agent_command = command

            # Generate a backup mechanical command
            mechanical_command = construct_command_mechanically(
                program,
                plan,
                original_files_list,
                required_params,
                log
            )

            if mechanical_command:
                log(f"DEBUG: Agent command: {agent_command}")
                log(f"DEBUG: Backup mechanical command: {mechanical_command}")

            # Use agent's command first (don't override)
            # mechanical_command will be used as fallback if agent's command fails

            # NEW: Fix common parameter mistakes for this program
            original_command = command
            command = fix_program_parameters(command, program)
            if command != original_command:
                log(f"DEBUG: Fixed parameters: {original_command} -> {command}")


            # --- FIX BASENAME-ONLY REFERENCES ---
            # If the command uses just a basename but we have a full path, substitute it
            log(f"\nDEBUG: Command before basename fix: {command}")
            command_parts = command.split()
            fixed_parts = []
            for part in command_parts:
                # Check if this part is a file reference (contains . and looks like a filename)
                if '=' in part:
                    # Handle key=value pairs
                    key, value = part.split('=', 1)
                    value = value.strip('"\'')
                    value_basename = os.path.basename(value)
                    # If value is just a basename and we have a full path, use the full path
                    if value_basename == value and value_basename in basename_to_relpath:
                        fixed_parts.append(f"{key}={basename_to_relpath[value_basename]}")
                    else:
                        fixed_parts.append(part)
                elif '.' in part and not part.startswith('phenix.') and not part.startswith('-'):
                    # Could be a filename
                    part_clean = part.strip('"\'')
                    basename = os.path.basename(part_clean)
                    # Only fix if it's just the basename (not already a path)
                    if basename == part_clean and basename in basename_to_relpath:
                        fixed_parts.append(basename_to_relpath[basename])
                    else:
                        fixed_parts.append(part)
                else:
                    fixed_parts.append(part)

            command = ' '.join(fixed_parts)
            log(f"DEBUG: Command after basename fix: {command}")

            # Hard block - remove docked_model_file if it references a ligand
            LIGAND_FILE_PATTERNS = ['lig.pdb', 'ligand.pdb', '_lig.pdb', '_ligand.pdb', 'lig.cif', 'ligand.cif']

            if 'docked_model_file' in command:
                match = re.search(r'docked_model_file[=\s]+["\']?([^\s"\']+)', command)
                if match:
                    referenced_file = match.group(1).lower()
                    # Check against known ligand patterns
                    is_ligand = any(pattern in referenced_file for pattern in LIGAND_FILE_PATTERNS)
                    # Also check if filename contains 'lig' anywhere
                    is_ligand = is_ligand or 'lig' in os.path.basename(referenced_file).lower()

                    if is_ligand:
                        log(f"BLOCKING: Removing docked_model_file={match.group(1)} - ligand files cannot be used with this parameter")
                        command = re.sub(r'\s*input_files\.docked_model_file[=\s]+["\']?[^\s"\']+["\']?', '', command)
                        command = re.sub(r'\s*docked_model_file[=\s]+["\']?[^\s"\']+["\']?', '', command)
                        command = ' '.join(command.split())


            # Hard block - reject commands with placeholder paths
            placeholder_patterns = [
                '/path/to/',
                '/dev/null',
                '/input/',
                '/output/',
            ]

            # These should only match as standalone filenames, not substrings
            standalone_placeholder_files = [
                'model.pdb',
                'data.mtz',
                'input.pdb',
                'input.mtz',
            ]

            has_placeholder = False

            # Check path patterns (substring match is fine)
            for pattern in placeholder_patterns:
                if pattern in command:
                    has_placeholder = True
                    log(f"ERROR: Command contains placeholder '{pattern}'")
                    break

            # Check standalone file placeholders (must be exact filename match)
            if not has_placeholder:
                import re
                for placeholder in standalone_placeholder_files:
                    # Match only if it's a standalone filename (not part of a longer name)
                    # Looks for the placeholder preceded by space, =, or start, and followed by space or end
                    pattern = rf'(?:^|[\s=])({re.escape(placeholder)})(?:\s|$)'
                    if re.search(pattern, command):
                        has_placeholder = True
                        log(f"ERROR: Command contains placeholder '{placeholder}'")
                        break

            if has_placeholder:
                # Get actual available files for the error message
                pdb_files = [f for f in available_files if f.endswith('.pdb')]
                mtz_files = [f for f in available_files if f.endswith('.mtz')]
                history_text += (
                    f"\n\nSYSTEM ERROR: You used placeholder paths like '/path/to/' or 'model.pdb'. "
                    f"You MUST use ACTUAL filenames from the available files.\n"
                    f"Available PDB files: {pdb_files}\n"
                    f"Available MTZ files: {mtz_files}\n"
                    f"Example correct command: phenix.refine PHASER.1.pdb PHASER.1.mtz\n"
                )
                continue  # Retry

            # Hard block - remove ensemble.pdb if it references a ligand
            if 'ensemble.pdb' in command or 'ensemble_pdb' in command:
                match = re.search(r'ensemble[._]pdb[=\s]+["\']?([^\s"\']+)', command)
                if match:
                    referenced_file = match.group(1).lower()
                    is_ligand = any(pattern in referenced_file for pattern in LIGAND_FILE_PATTERNS)
                    is_ligand = is_ligand or 'lig' in os.path.basename(referenced_file).lower()

                    if is_ligand:
                        log(f"BLOCKING: Removing ensemble.pdb={match.group(1)} - ligand files cannot be used as phaser ensemble")
                        command = re.sub(r'\s*phaser\.ensemble\.pdb[=\s]+["\']?[^\s"\']+["\']?', '', command)
                        command = re.sub(r'\s*ensemble\.pdb[=\s]+["\']?[^\s"\']+["\']?', '', command)
                        command = ' '.join(command.split())


            # Validate it's actually a Phenix command
            from libtbx.langchain.validation.core_validator import validate_is_phenix_command

            is_valid, error_msg = validate_is_phenix_command(command)
            if not is_valid:
              log(f"ERROR: Invalid command generated: {error_msg}")
              log(f"Raw output was: {command_obj.content[:200]}...")
              # Retry with feedback
              history_text += f"\n\nSYSTEM NOTE: Your previous output was not a valid Phenix command. {error_msg}. You MUST output ONLY a valid phenix command starting with 'phenix.' - no explanations, no thinking, just the command.\n"
              continue  # This continues the retry loop


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
        mechanical_command=mechanical_command,
        explanation=f"PLAN: {long_term_plan}\n\nREASONING: {plan.get('reasoning', '')}",
        program=program,
        strategy=plan.get('strategy_details', ""),
        process_log="\n".join(process_log),
        error=None
    )

