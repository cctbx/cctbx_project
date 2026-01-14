"""
Vector store creation and persistence.

This module handles:
- Creating Chroma vector stores from documents
- Batched embedding with rate limit handling
- Persisting vector stores to disk

Usage:
    from libtbx.langchain.rag.vector_store import create_and_persist_db

    vectorstore = create_and_persist_db(docs, embeddings, './docs_db')
"""
from __future__ import absolute_import, division, print_function

import time
from typing import List, Iterable
assert Iterable is not None

import chromadb
from langchain_core.documents import Document
from langchain_chroma import Chroma

from libtbx.langchain.rag.document_loader import _custom_chunker


# =============================================================================
# Helper Functions
# =============================================================================

def _iter_batches(seq: List[Document], size: int) -> Iterable[List[Document]]:
    """Yield successive batches from sequence."""
    for i in range(0, len(seq), size):
        yield seq[i:i+size]


# =============================================================================
# Vector Store Creation
# =============================================================================

def create_and_persist_db(
    docs: List[Document],
    embeddings,
    db_dir: str = "./docs_db",
    *,
    add_batch_size: int = 200,
    pause_between_batches: float = 2.0,
    max_attempts: int = 6,
    max_backoff: float = 60.0
) -> Chroma:
    """
    Build a Chroma vector store with batching to avoid SQLite limits and API rate limits.

    Args:
        docs: List of documents to embed and store
        embeddings: Embedding model to use
        db_dir: Directory to persist the database
        add_batch_size: Number of documents per batch
        pause_between_batches: Seconds to wait between batches
        max_attempts: Maximum retry attempts per batch
        max_backoff: Maximum backoff time in seconds

    Returns:
        Chroma: The created vector store, or None if creation failed

    Example:
        from libtbx.langchain.core.llm import get_llm_and_embeddings
        from libtbx.langchain.rag.document_loader import load_all_docs_from_folder
        from libtbx.langchain.rag.vector_store import create_and_persist_db

        _, embeddings = get_llm_and_embeddings()
        docs, _ = load_all_docs_from_folder('./docs/')
        vectorstore = create_and_persist_db(docs, embeddings, './docs_db')
    """
    docs_chunks = _custom_chunker(docs)

    # Deduplicate
    seen = set()
    unique_docs: List[Document] = []
    for d in docs_chunks:
        key = (d.page_content.strip(), tuple(sorted((d.metadata or {}).items())))
        if key not in seen:
            seen.add(key)
            unique_docs.append(d)

    client = chromadb.PersistentClient(path=db_dir)

    vectorstore = Chroma(
        client=client,
        collection_name="docs",
        embedding_function=embeddings,
    )

    for batch in _iter_batches(unique_docs, add_batch_size):
        backoff = 2.0
        for attempt in range(1, max_attempts + 1):
            try:
                vectorstore.add_documents(batch)
                break
            except Exception as e:
                msg = str(e).lower()
                if "api_key_ip_address_blocked" in msg:
                    print("Google AI API key does not allow access from this server")
                    return None
                elif "you exceeded your current quota" in msg:
                    if "limit: 0" in msg:
                        print("Google AI API key has a zero quota")
                        return None
                    else:
                        print("Google AI API quota exceeded")
                        return None
                elif "429" in msg or "rate" in msg:
                    time.sleep(backoff)
                    backoff = min(backoff * 2.0, max_backoff)
                    if attempt == max_attempts:
                        raise RuntimeError("Failed to use Google API (Rate Limit)")
                else:
                    raise RuntimeError(f"Failed to add batch: {e}")
        time.sleep(pause_between_batches)

    print(f"Database stored in: {db_dir}")
    return vectorstore
