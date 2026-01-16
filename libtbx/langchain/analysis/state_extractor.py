"""
Project state extraction from log summaries.

This module uses LLM to identify updates to the persistent project database
based on log summaries. It extracts:
- New data files
- Crystallographic parameters (space group, resolution, etc.)
- Cryo-EM parameters
- Model files with metrics

Usage:
    from libtbx.langchain.analysis.state_extractor import extract_project_state_updates

    updates = await extract_project_state_updates(log_summary, current_state, llm)
    # updates is a dict like: {"crystallography": {"space_group": "P 21 21 21"}}
"""
from __future__ import absolute_import, division, print_function

import json
import re

from langchain_core.prompts import PromptTemplate


# =============================================================================
# Prompt Templates
# =============================================================================

def get_state_update_prompt() -> PromptTemplate:
    """Returns the prompt for extracting state updates from a log summary."""
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


# =============================================================================
# Main Extraction Function
# =============================================================================

async def extract_project_state_updates(log_summary, current_state, llm):
    """
    Uses LLM to identify updates to the persistent project database.

    Includes validation to prevent hallucinated file names from corrupting state.
    """
    try:
        if not log_summary:
            return {}

        # Don't extract state from initialization summaries
        initialization_indicators = [
            "No log file to analyze",
            "project initialization phase",
            "No log content",
            "This is the project initialization",
        ]

        summary_lower = log_summary.lower()
        for indicator in initialization_indicators:
            if indicator.lower() in summary_lower:
                print("Skipping state extraction - initialization phase")
                return {}

        prompt = get_state_update_prompt()
        chain = prompt | llm

        # Serialize current state for context
        state_str = json.dumps(current_state, indent=2) if current_state else "Empty"

        response = await chain.ainvoke({
            "current_state": state_str,
            "log_summary": log_summary
        })

        # Content is likely JSON string, possibly wrapped in markdown
        text = response.content.strip()
        if "```" in text:
            match = re.search(r'\{.*\}', text, re.DOTALL)
            if match:
                text = match.group(0)

        updates = json.loads(text)

        # === VALIDATION PHASE ===
        # Detect which program was run from the summary
        program_name = _detect_program_from_summary(log_summary)

        # Validate and filter the updates
        validated_updates = _validate_state_updates(updates, program_name, log_summary)

        return validated_updates

    except Exception as e:
        print(f"State extraction failed: {e}")
        return {}


def _detect_program_from_summary(log_summary):
    """Detect which Phenix program produced this summary."""
    summary_lower = log_summary.lower()

    if 'xtriage' in summary_lower:
        return 'xtriage'
    elif 'phaser' in summary_lower:
        return 'phaser'
    elif 'refine' in summary_lower:
        return 'refine'
    elif 'autobuild' in summary_lower:
        return 'autobuild'
    elif 'ligandfit' in summary_lower or 'ligand_fit' in summary_lower:
        return 'ligandfit'
    elif 'predict_and_build' in summary_lower:
        return 'predict_and_build'
    elif 'mtriage' in summary_lower:
        return 'mtriage'
    else:
        return 'unknown'


def _validate_state_updates(updates, program_name, log_summary):
    """
    Filter out updates that couldn't have come from the detected program.
    This prevents hallucinated file names from corrupting the state.
    """
    if not updates:
        return {}

    validated = {}

    # Programs that CAN produce model files
    model_producing_programs = ['phaser', 'refine', 'autobuild', 'ligandfit',
                                 'predict_and_build', 'dock_and_rebuild']

    # Programs that CAN produce MTZ with R-free flags
    rfree_mtz_programs = ['refine', 'autobuild', 'reflection_file_converter']

    # Diagnostic programs that produce NO structural files
    diagnostic_programs = ['xtriage', 'mtriage', 'mtz.dump', 'explore_metric_symmetry']

    # Hallucination markers - generic placeholder names
    hallucination_markers = [
        'refined_001.mtz',
        'model.pdb',
        'data.mtz',
        'file1.', 'file2.',
        '/path/to/',
    ]

    for key, value in updates.items():
        if key == 'model_list':
            # Only accept model_list updates from programs that create models
            if program_name in diagnostic_programs:
                print(f"WARNING: Rejecting model_list from {program_name} (diagnostic tool)")
                continue

            # Filter individual models for hallucinations
            if isinstance(value, list):
                valid_models = []
                for model in value:
                    if isinstance(model, dict):
                        filename = model.get('file', '')
                        # Check for hallucination markers
                        is_hallucinated = any(marker in filename for marker in hallucination_markers)
                        # Check if filename is explicitly mentioned in log with "written" context
                        is_in_log = filename.lower() in log_summary.lower()
                        written_context = any(
                            phrase in log_summary.lower()
                            for phrase in ['written to', 'output:', 'saved to', 'solution #']
                        )

                        if is_hallucinated:
                            print(f"WARNING: Rejecting hallucinated model: {filename}")
                        elif is_in_log and (written_context or program_name in model_producing_programs):
                            valid_models.append(model)
                        elif not is_in_log:
                            print(f"WARNING: Rejecting model not in log: {filename}")
                        else:
                            valid_models.append(model)

                if valid_models:
                    validated['model_list'] = valid_models

        elif key == 'crystallography':
            if isinstance(value, dict):
                valid_cryst = {}
                for cryst_key, cryst_val in value.items():
                    # Check for hallucinated mtz files from wrong programs
                    if cryst_key in ['data_with_free_r_flags', 'data_file']:
                        if program_name in diagnostic_programs:
                            print(f"WARNING: Rejecting {cryst_key} from {program_name}")
                            continue
                        # Check for hallucination markers
                        if any(marker in str(cryst_val) for marker in hallucination_markers):
                            print(f"WARNING: Rejecting hallucinated {cryst_key}: {cryst_val}")
                            continue

                    # Space group and resolution are OK from diagnostic tools
                    valid_cryst[cryst_key] = cryst_val

                if valid_cryst:
                    validated['crystallography'] = valid_cryst

        else:
            # Other keys pass through with basic hallucination check
            str_val = str(value)
            if not any(marker in str_val for marker in hallucination_markers):
                validated[key] = value

    return validated

