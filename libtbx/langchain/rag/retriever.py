"""
RAG (Retrieval-Augmented Generation) retriever and chain functions.

This module handles:
- Loading persisted vector databases
- Creating reranking retrievers (using Cohere)
- Building RAG chains for various purposes

Usage:
    from libtbx.langchain.rag import load_persistent_db, create_reranking_retriever

    vectorstore = load_persistent_db(embeddings, db_dir='./docs_db')
    retriever = create_reranking_retriever(vectorstore, llm)
"""
from __future__ import absolute_import, division, print_function

import os

import cohere
from langchain_cohere import CohereRerank
from langchain.retrievers import ContextualCompressionRetriever
from langchain_core.runnables import RunnablePassthrough
from langchain_core.prompts import PromptTemplate
from langchain_chroma import Chroma


# =============================================================================
# Database Loading
# =============================================================================

def load_persistent_db(embeddings, db_dir: str = "./docs_db"):
    """
    Loads a persisted Chroma vector store from disk.

    Args:
        embeddings: Embedding model to use
        db_dir: Path to the persisted database directory

    Returns:
        Chroma: The loaded vector store

    Raises:
        FileNotFoundError: If database directory doesn't exist

    Example:
        from libtbx.langchain.core import get_llm_and_embeddings
        _, embeddings = get_llm_and_embeddings()
        vectorstore = load_persistent_db(embeddings, './docs_db')
    """
    if not os.path.exists(db_dir):
        raise FileNotFoundError(f"Database directory not found at '{db_dir}'.")
    return Chroma(persist_directory=db_dir, embedding_function=embeddings)


# =============================================================================
# Retriever Creation
# =============================================================================

def create_reranking_retriever(vectorstore, llm, timeout: int = 60, top_n: int = 8):
    """
    Creates a retriever with Cohere reranking for improved relevance.

    Args:
        vectorstore: Chroma vector store to retrieve from
        llm: Language model (not directly used, but kept for API compatibility)
        timeout: Timeout in seconds for Cohere API
        top_n: Number of top results to return after reranking

    Returns:
        ContextualCompressionRetriever: Retriever with reranking

    Example:
        retriever = create_reranking_retriever(vectorstore, llm, top_n=5)
        docs = retriever.invoke("How do I refine a model?")
    """
    base_retriever = vectorstore.as_retriever(search_kwargs={"k": 20})

    cohere_client = cohere.ClientV2(
        api_key=os.getenv("COHERE_API_KEY"),
        timeout=timeout
    )

    reranker = CohereRerank(client=cohere_client, model="rerank-english-v3.0", top_n=top_n)

    compression_retriever = ContextualCompressionRetriever(
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
