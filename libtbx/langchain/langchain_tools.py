"""
Backward compatibility layer for langchain_tools.

This module re-exports all functions from the new modular structure
to maintain backward compatibility with existing code.

New code should import directly from submodules:
    from libtbx.langchain.core import get_llm_and_embeddings
    from libtbx.langchain.analysis import summarize_log_text
    from libtbx.langchain.agent import generate_next_move

Legacy code can continue to use:
    from libtbx.langchain import langchain_tools as lct
    lct.get_llm_and_embeddings(...)
"""
from __future__ import division

# =============================================================================
# Required setup (must be at top)
# =============================================================================
import nest_asyncio
nest_asyncio.apply()

import os
os.environ['GRPC_ENABLE_FORK_SUPPORT'] = "false"

# =============================================================================
# Re-export from core module
# =============================================================================
from libtbx.langchain.core import get_llm_and_embeddings
from libtbx.langchain.core import get_expensive_llm
from libtbx.langchain.core import get_cheap_llm

# =============================================================================
# Re-export from analysis module
# =============================================================================
from libtbx.langchain.analysis import (
    summarize_log_text,
    analyze_log_summary,
    extract_project_state_updates,
    get_log_info,
    get_log_map_prompt,
    get_log_combine_prompt,
    get_log_analysis_prompt,
    get_state_update_prompt,
    get_chunk_size,
)

# Internal helpers from analysis
from libtbx.langchain.analysis.summarizer import (
    _custom_log_chunker,
    _iter_batches,
)

# =============================================================================
# Re-export from rag module
# =============================================================================
from libtbx.langchain.rag import (
    load_all_docs_from_folder,
    load_specific_docs,
    create_and_persist_db,
    load_persistent_db,
    create_reranking_retriever,
    create_reranking_rag_chain,
    create_log_analysis_chain,
    create_keyword_lookup_chain,
    PhenixHTMLLoader,
)

# Internal helpers from rag
from libtbx.langchain.rag.document_loader import (
    _custom_chunker,
)

# =============================================================================
# Re-export from validation module
# =============================================================================
from libtbx.langchain.validation import (
    validate_phenix_command,
    fix_command_syntax,
    mtz_has_rfree_flags,
    add_rfree_generation_if_needed,
)

# =============================================================================
# Re-export from knowledge module
# =============================================================================
from libtbx.langchain.knowledge import (
    get_phenix_program_list,
    get_keywords_as_phil_string,
    get_strategic_planning_prompt,
    get_command_writer_prompt,
    get_keywords_prompt,
    get_docs_query_prompt,
)

# =============================================================================
# Re-export from agent module
# =============================================================================
from libtbx.langchain.agent import (
    generate_next_move,
    get_run_history,
    get_memory_file_path,
    load_learned_memory,
    save_learned_memory,
    learn_from_history,
    get_program_keywords,
    extract_output_files,
    AgentSession,
)

# =============================================================================
# Re-export from utils module
# =============================================================================
from libtbx.langchain.utils import (
    find_text_block,
    get_processed_log_dict,
    query_docs,
)

# =============================================================================
# Module-level state (for backward compatibility)
# =============================================================================
_last_query_time = 0

# =============================================================================
# Verification
# =============================================================================
if __name__ == "__main__":
    print("langchain_tools compatibility layer loaded successfully")
    print(f"get_llm_and_embeddings: {get_llm_and_embeddings}")
    print(f"generate_next_move: {generate_next_move}")
    print(f"summarize_log_text: {summarize_log_text}")
