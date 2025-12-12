"""
LLM and Embeddings setup for the Phenix Crystallography Agent.

This module handles initialization of Language Models and Embedding models
from different providers (Google, OpenAI).

Usage:
    from libtbx.langchain.core import get_llm_and_embeddings

    llm, embeddings = get_llm_and_embeddings(provider='google')

    # Or with specific models:
    llm, embeddings = get_llm_and_embeddings(
        provider='google',
        llm_model_name='gemini-2.5-pro',
        temperature=0.2
    )
"""
from __future__ import absolute_import, division, print_function

import os


def get_llm_and_embeddings(
    provider: str = "google",
    llm_model_name: str = None,
    embedding_model_name: str = None,
    temperature: float = 0.1,
    timeout: int = 60,
    batch_size: int = 100
):
    """
    Initialize LLM and Embeddings from the specified provider.

    Args:
        provider: Which AI provider to use ('google' or 'openai')
        llm_model_name: Specific model name, or None for default
        embedding_model_name: Specific embedding model, or None for default
        temperature: Sampling temperature (0.0 = deterministic, 1.0 = creative)
        timeout: Request timeout in seconds
        batch_size: Batch size for embedding requests

    Returns:
        tuple: (llm, embeddings) - initialized model objects

    Raises:
        ValueError: If provider is not supported

    Example:
        # Use defaults (Google, fast model)
        llm, embeddings = get_llm_and_embeddings()

        # Use expensive model for complex reasoning
        llm, embeddings = get_llm_and_embeddings(
            provider='google',
            llm_model_name='gemini-2.5-pro'
        )

        # Use OpenAI
        llm, embeddings = get_llm_and_embeddings(provider='openai')
    """
    # Import here to avoid import errors if packages not installed
    from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAIEmbeddings
    from langchain_openai import ChatOpenAI, OpenAIEmbeddings

    provider = provider.lower()

    if provider == "google":
        if llm_model_name is None:
            llm_model_name = "gemini-2.5-flash-lite"
        if embedding_model_name is None:
            embedding_model_name = "models/embedding-001"

        llm = ChatGoogleGenerativeAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=0
        )
        embeddings = GoogleGenerativeAIEmbeddings(
            model=embedding_model_name,
            timeout=timeout,
            batch_size=batch_size,
            google_api_key=os.getenv("GOOGLE_API_KEY")
        )
        print("Using Google Gemini models.")

    elif provider == "openai":
        if llm_model_name is None:
            llm_model_name = "gpt-5-nano"
        if embedding_model_name is None:
            embedding_model_name = "text-embedding-3-small"

        llm = ChatOpenAI(
            model=llm_model_name,
            temperature=temperature,
            timeout=timeout,
            max_retries=2
        )
        embeddings = OpenAIEmbeddings(
            model=embedding_model_name,
            chunk_size=batch_size
        )
        print("Using OpenAI models.")

    else:
        raise ValueError(f"Unsupported provider: '{provider}'. Choose 'google' or 'openai'.")

    return llm, embeddings

