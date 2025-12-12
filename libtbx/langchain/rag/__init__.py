"""
RAG (Retrieval-Augmented Generation) module.

This module contains:
- document_loader.py: Document loading and chunking
- vector_store.py: Vector store creation and persistence
- retriever.py: Retriever and chain creation
"""
from __future__ import absolute_import, division, print_function

from libtbx.langchain.rag.document_loader import (
    load_all_docs_from_folder,
    load_specific_docs,
    PhenixHTMLLoader,
)

from libtbx.langchain.rag.vector_store import (
    create_and_persist_db,
)

from libtbx.langchain.rag.retriever import (
    load_persistent_db,
    create_reranking_retriever,
    create_reranking_rag_chain,
    create_log_analysis_chain,
    create_keyword_lookup_chain,
)

__all__ = [
    # Document loading
    'load_all_docs_from_folder',
    'load_specific_docs',
    'PhenixHTMLLoader',
    # Vector store
    'create_and_persist_db',
    # Retriever
    'load_persistent_db',
    'create_reranking_retriever',
    'create_reranking_rag_chain',
    'create_log_analysis_chain',
    'create_keyword_lookup_chain',
]

