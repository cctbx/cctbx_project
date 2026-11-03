"""
RAG (Retrieval-Augmented Generation) retriever and chain functions.

This module handles:
- Loading persisted vector databases
- Creating reranking retrievers (using FlashRank, a local cross-encoder)
- Building RAG chains for various purposes

Usage:
    from libtbx.langchain.rag.retriever import load_persistent_db, create_reranking_retriever

    vectorstore = load_persistent_db(embeddings, db_dir='./docs_db')
    retriever = create_reranking_retriever(vectorstore, llm)
"""
from __future__ import absolute_import, division, print_function

import os

from langchain_core.retrievers import BaseRetriever
from langchain_core.runnables import RunnablePassthrough
from langchain_core.prompts import PromptTemplate
from langchain_core.documents import Document, BaseDocumentCompressor
# v118.G1: chromadb has known protobuf version conflicts in some envs
# (TypeError: Descriptors cannot be created directly).  Lazy-import via
# shared helper so this module imports cleanly in any env; failure is
# deferred to function-call time with a clear install hint.
# See docs/DEVELOPER_GUIDE.md "Optional dependency handling".
from libtbx.langchain.rag._chroma_resilience import (
    ensure_chroma,
    chroma_unavailable_error,
)
assert PromptTemplate is not None


class _CompressionRetriever(BaseRetriever):
    """Minimal replacement for langchain ContextualCompressionRetriever.

    Retrieves documents from a base retriever then reranks/compresses them
    using a document compressor. Implements BaseRetriever so it works in
    LCEL chains (e.g., retriever | format_docs | llm).

    This replaces langchain.retrievers.ContextualCompressionRetriever which
    was removed in langchain 1.0 (moved to langchain_classic).
    """
    base_compressor: object  # BaseDocumentCompressor
    base_retriever: object   # BaseRetriever

    class Config:
        arbitrary_types_allowed = True

    def _get_relevant_documents(self, query, *, run_manager=None):
        docs = self.base_retriever.invoke(query)
        compressed = self.base_compressor.compress_documents(docs, query)
        return list(compressed)


# Module-level Ranker cache: building a flashrank Ranker may download/load model
# weights, so we build each (model_name) once and reuse it.  Kept OUTSIDE the
# pydantic model below to avoid private-attribute machinery on
# BaseDocumentCompressor (a pydantic v2 model in langchain-core 1.x).
_RANKER_CACHE = {}


def _get_flashrank_ranker(model_name):
    """Return a cached flashrank Ranker for model_name, building it on first use."""
    ranker = _RANKER_CACHE.get(model_name)
    if ranker is None:
        from flashrank import Ranker
        ranker = Ranker(model_name=model_name)
        _RANKER_CACHE[model_name] = ranker
    return ranker


class PhenixFlashrankCompressor(BaseDocumentCompressor):
    """Local FlashRank reranker -- drop-in replacement for the deprecated
    langchain-community FlashrankRerank, with no langchain-community dependency.

    Reranks retrieved documents with a local cross-encoder (default
    ms-marco-MiniLM-L-12-v2, ~34MB, CPU, no API key) via the `flashrank` package
    directly.  Subclasses BaseDocumentCompressor and implements compress_documents
    so it plugs into _CompressionRetriever unchanged.

    Behaviour matches the old langchain-community FlashrankRerank defaults: results
    are filtered by score_threshold (default 0.0, i.e. drop negative-score docs) and
    then truncated to top_n.
    """
    # All fields are standard scalars, so no `arbitrary_types_allowed` / Config is
    # needed (the flashrank Ranker is held in the module-level cache, NOT as a
    # field).  BaseDocumentCompressor is a pydantic v2 model and accepts these
    # defaulted fields directly, exactly as the community FlashrankRerank did.
    model_name: str = "ms-marco-MiniLM-L-12-v2"
    top_n: int = 8
    score_threshold: float = 0.0   # matches community FlashrankRerank default

    # Only compress_documents is abstract in BaseDocumentCompressor; the async
    # acompress_documents has a concrete default that delegates to this method via
    # run_in_executor (verified in langchain_core source), so we don't implement
    # the async variant.
    def compress_documents(self, documents, query, callbacks=None):
        documents = list(documents)
        if not documents:
            return []
        from flashrank import RerankRequest
        # Tag each doc with its index so we can map flashrank's result back to the
        # exact Document object (preserving page_content + metadata).
        passages = [
            {"id": i, "text": doc.page_content, "meta": doc.metadata or {}}
            for i, doc in enumerate(documents)
        ]
        ranker = _get_flashrank_ranker(self.model_name)
        ranked = ranker.rerank(RerankRequest(query=query, passages=passages))
        # flashrank returns ALL passages, already sorted best-first.  Apply the
        # score_threshold filter (community default 0.0) THEN truncate to top_n --
        # the order the community wrapper used.
        out = []
        for item in ranked:
            score = float(item.get("score", 0.0))   # raw is np.float32
            if score < self.score_threshold:
                continue
            src = documents[item["id"]]
            md = dict(src.metadata or {})
            md["relevance_score"] = score
            out.append(Document(page_content=src.page_content, metadata=md))
            if len(out) >= self.top_n:
                break
        return out


# =============================================================================
# Database Loading
# =============================================================================

# In libtbx/langchain/rag/retriever.py

def load_persistent_db(embeddings, db_dir: str = "./docs_db", collection_name: str = "docs"):
    """
    Loads a persisted Chroma vector store from disk.
    Now defaults to collection_name="docs" to match your build script.
    """
    # v118.G1: probe chromadb lazily.  When the chromadb stack is
    # unavailable (version conflicts, not installed, etc.), raise a
    # clear RuntimeError with install hint instead of letting a
    # confusing TypeError bubble up from inside opentelemetry.
    _chromadb_mod, Chroma = ensure_chroma()
    if Chroma is None:
        raise chroma_unavailable_error()

    if not os.path.exists(db_dir):
        raise FileNotFoundError(f"Database directory not found at '{db_dir}'.")

    print(f"Loading Vector Store from '{db_dir}' (Collection: '{collection_name}')...")

    return Chroma(
        persist_directory=db_dir,
        embedding_function=embeddings,
        collection_name=collection_name  # <--- THIS IS THE CRITICAL FIX
    )

# =============================================================================
# Retriever Creation
# =============================================================================

def create_reranking_retriever(vectorstore, llm, timeout=60, top_n=8):
    """
    Creates a retriever with FlashRank reranking for improved relevance.

    Uses a local cross-encoder model (ms-marco-MiniLM-L-12-v2, ~34MB) to
    rerank documents. Runs on CPU, no API key required.

    Args:
        vectorstore: Chroma vector store to retrieve from
        llm: Language model (not directly used, kept for API compatibility)
        timeout: Unused, kept for backward compatibility
        top_n: Number of top results to return after reranking

    Returns:
        BaseRetriever: Retriever with reranking

    Example:
        retriever = create_reranking_retriever(vectorstore, llm, top_n=5)
        docs = retriever.invoke("How do I refine a model?")
    """
    base_retriever = vectorstore.as_retriever(search_kwargs={"k": 20})

    # Local FlashRank compressor (no langchain-community dependency).
    reranker = PhenixFlashrankCompressor(
        model_name="ms-marco-MiniLM-L-12-v2",
        top_n=top_n,
    )

    compression_retriever = _CompressionRetriever(
        base_compressor=reranker, base_retriever=base_retriever
    )

    return compression_retriever


# =============================================================================
# Chain Creation
# =============================================================================

def create_reranking_rag_chain(retriever, llm, prompt: PromptTemplate):
    """
    Creates a RAG chain with the given retriever, LLM, and prompt.

    Args:
        retriever: Retriever to use for fetching context
        llm: Language model for generation
        prompt: Prompt template with 'context' and 'input' variables

    Returns:
        Runnable chain
    """
    def format_docs(docs):
        return "\n\n".join(doc.page_content for doc in docs)

    rag_chain = (
        {"context": retriever | format_docs, "input": RunnablePassthrough()}
        | prompt
        | llm
    )
    return rag_chain


def create_log_analysis_chain_debug(retriever, llm, prompt: PromptTemplate):
    """
    Creates a chain specifically for log analysis with RAG context.
    """

    # --- DEFINE THE DEBUG FUNCTION INSIDE HERE ---
    def retrieve_and_print(x):
        query = x["input"]

        # Now it can see 'retriever' because it's in the same scope
        docs = retriever.invoke(query)

        # Debug print
        print(f"\n[DEBUG] RAG Query: '{query}'")
        print(f"[DEBUG] Found {len(docs)} docs.")
        for i, doc in enumerate(docs):
            src = doc.metadata.get('source', 'unknown')
            snippet = doc.page_content[:100].replace('\n', ' ')
            print(f"   {i+1}. [{os.path.basename(src)}] {snippet}...")
        print("-" * 40)

        return docs
    # ---------------------------------------------

    inputs = {
        # Use the local function here
        "context": retrieve_and_print,
        "log_summary": lambda x: x["log_summary"]
    }

    analysis_rag_chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return analysis_rag_chain

def create_log_analysis_chain(retriever, llm, prompt: PromptTemplate):
    """
    Creates a chain specifically for log analysis with RAG context.

    Args:
        retriever: Retriever to use for fetching documentation context
        llm: Language model for analysis
        prompt: Prompt template with 'context', 'input', and 'log_summary' variables

    Returns:
        Runnable chain for log analysis
    """
    inputs = {
        "context": lambda x: retriever.invoke(x["input"]),
        "log_summary": lambda x: x["log_summary"]
    }

    analysis_rag_chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return analysis_rag_chain


def create_keyword_lookup_chain(retriever, llm, prompt: PromptTemplate):
    """
    Creates a chain for looking up program keywords/parameters.

    Args:
        retriever: Retriever to use
        llm: Language model
        prompt: Prompt template with 'context' and 'program_name' variables

    Returns:
        Runnable chain for keyword lookup
    """
    inputs = {
        "context": lambda x: retriever.invoke(f"command line examples and parameters for {x['program_name']}"),
        "program_name": lambda x: x["program_name"]
    }

    chain = (
        RunnablePassthrough.assign(**inputs)
        | prompt
        | llm
    )
    return chain
