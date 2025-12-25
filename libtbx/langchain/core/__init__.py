"""
Core module - fundamental types and infrastructure.

This module contains:
- Data types used throughout the system (types.py)
- LLM setup and configuration (llm.py)
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.core.types import (
    AgentPlan,
    ValidationResult,
    CommandResult,
    ProjectState,
    NextMove,
)

from libtbx.langchain.core.llm import (
    get_llm_and_embeddings,
    get_expensive_llm,
    get_cheap_llm,
    get_ollama_llm,
    get_ollama_fast_llm,
    get_ollama_smart_llm,
)

__all__ = [
    # Types
    'AgentPlan',
    'ValidationResult',
    'CommandResult',
    'ProjectState',
    'NextMove',
    # LLM
    'get_llm_and_embeddings',
    'get_expensive_llm',
    'get_cheap_llm',
    'get_ollama_llm',
    'get_ollama_fast_llm',
    'get_ollama_smart_llm',
]

