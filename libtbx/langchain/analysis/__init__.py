"""
Analysis module - log summarization and analysis.

This module contains:
- summarizer.py: Log file summarization (map-reduce)
- state_extractor.py: Project state extraction from log summaries
- analyzer.py: Log analysis with RAG
- log_info.py: High-level log information extraction
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.analysis.summarizer import (
    summarize_log_text,
    get_log_map_prompt,
    get_log_combine_prompt,
    get_chunk_size,
)

from libtbx.langchain.analysis.state_extractor import (
    extract_project_state_updates,
    get_state_update_prompt,
)

from libtbx.langchain.analysis.analyzer import (
    analyze_log_summary,
    get_log_analysis_prompt,
)

from libtbx.langchain.analysis.log_info import (
    get_log_info,
)

__all__ = [
    # Summarizer
    'summarize_log_text',
    'get_log_map_prompt',
    'get_log_combine_prompt',
    'get_chunk_size',
    # State extractor
    'extract_project_state_updates',
    'get_state_update_prompt',
    # Analyzer
    'analyze_log_summary',
    'get_log_analysis_prompt',
    # Log info
    'get_log_info',
]
