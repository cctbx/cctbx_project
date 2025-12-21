"""
Project state extraction from log summaries.

This module uses LLM to identify updates to the persistent project database
based on log summaries. It extracts:
- New data files
- Crystallographic parameters (space group, resolution, etc.)
- Cryo-EM parameters
- Model files with metrics

Usage:
    from libtbx.langchain.analysis import extract_project_state_updates

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

    Args:
        log_summary: Summary text from a recent job
        current_state: Current project state dict
        llm: Language model to use for extraction

    Returns:
        dict: Updates to apply to project state, or empty dict if extraction fails

    Example:
        updates = await extract_project_state_updates(
            log_summary="Refinement complete. R-free = 0.22. Output: refine_1.pdb",
            current_state={"crystallography": {"space_group": "P 21 21 21"}},
            llm=llm
        )
        # Returns: {"model_list": [{"file": "refine_1.pdb", "r_free": 0.22}]}
    """
    try:
        if not log_summary:
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

        return json.loads(text)

    except Exception as e:
        print(f"State extraction failed: {e}")
        return {}
