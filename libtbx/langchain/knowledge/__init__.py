"""
Knowledge module - Phenix program information and prompts.

This module contains:
- phenix_programs.py: Program discovery and keyword introspection
- prompts.py: All prompt templates for the agent
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.knowledge.phenix_programs import (
    get_phenix_program_list,
    get_keywords_as_phil_string,
)

from libtbx.langchain.knowledge.prompts import (
    get_strategic_planning_prompt,
    get_command_writer_prompt,
    get_keywords_prompt,
    get_docs_query_prompt,
)

__all__ = [
    # Program discovery
    'get_phenix_program_list',
    'get_keywords_as_phil_string',
    # Prompts
    'get_strategic_planning_prompt',
    'get_command_writer_prompt',
    'get_keywords_prompt',
    'get_docs_query_prompt',
]
