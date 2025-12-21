"""
Utility functions module.

This module contains:
- text_processing.py: Text extraction and processing
- query.py: Documentation query utilities
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.utils.text_processing import (
    find_text_block,
    get_processed_log_dict,
)

from libtbx.langchain.utils.query import (
    query_docs,
)

__all__ = [
    'find_text_block',
    'get_processed_log_dict',
    'query_docs',
]
