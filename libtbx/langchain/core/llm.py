"""
LLM and Embeddings setup for the Phenix Crystallography Agent.

This module handles initialization of Language Models and Embedding models
from different providers (Ollama, Google, OpenAI).

Usage:
    from libtbx.langchain.core import get_llm_and_embeddings

    # Use local Ollama (recommended)
    llm, embeddings = get_llm_and_embeddings(provider='ollama')

    # Or with Google
    llm, embeddings = get_llm_and_embeddings(provider='google')
"""
from __future__ import absolute_import, division, print_function

import os


def get_llm_and_embeddings(
    provider: str = None,
    llm_model_name: str = None,
    embedding_model_name: str = None,
    temperature: float = 0.1,
    timeout: int = 120,
    batch_size: int = 100,
    ollama_base_url: str = None,
    num_ctx: int = 8192,
):
    """
    Initialize LLM and Embeddings from the specified provider.

    Args:
        provider: Which AI provider to use ('ollama', 'google', or 'openai').
                  If None, reads from LLM_PROVIDER env var, defaults to 'ollama'.
        llm_model_name: Specific model name, or None for default
        embedding_model_name: Specific embedding model, or None for default
        temperature: Sampling temperature (0.0 = deterministic, 1.0 = creative)
        timeout: Request timeout in seconds
        batch_size: Batch size for embedding requests
        ollama_base_url: Base URL for Ollama server (default from env or localhost)
        num_ctx: Context window size for Ollama models (default: 8192)

    Returns:
        tuple: (llm, embeddings) - initialized model objects

    Raises:
        ValueError: If provider is not supported
    """
    # Get provider from argument or environment variable
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")
    provider = provider.lower()

    if provider == "ollama":
        from langchain_ollama import ChatOllama, OllamaEmbeddings

        if llm_model_name is None:
            llm_model_name = os.getenv("OLLAMA_LLM_MODEL", "llama3.1:70b")
        if embedding_model_name is None:
            embedding_model_name = os.getenv("OLLAMA_EMBED_MODEL", "nomic-embed-text")
        if ollama_base_url is None:
            ollama_base_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

        llm = ChatOllama(
            model=llm_model_name,
            base_url=ollama_base_url,
            temperature=temperature,
            num_ctx=num_ctx,
            num_predict=4096,
        )
        embeddings = OllamaEmbeddings(
            model=embedding_model_name,
            base_url=ollama_base_url,
        )
        print(f"Using Ollama at {ollama_base_url}")
        print(f"  LLM: {llm_model_name}, Embeddings: {embedding_model_name}")

    elif provider == "google":
        from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAIEmbeddings

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
        print(f"Using Google Gemini: {llm_model_name}")

    elif provider == "openai":
        from langchain_openai import ChatOpenAI, OpenAIEmbeddings

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
        print(f"Using OpenAI: {llm_model_name}")

    else:
        raise ValueError(f"Unsupported provider: '{provider}'. Choose 'ollama', 'google', or 'openai'.")

    return llm, embeddings

# Get expensive or cheap LLMs

def get_expensive_llm(provider = None, timeout = None):
      try:
        if provider == "google":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='gemini-2.5-pro')
          print(f"Using expensive model for analysis: {expensive_llm.model}")
        elif provider == "openai":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='gpt-5')
          print(f"Using expensive model for analysis: {expensive_llm.model_name}")
        elif provider == "ollama":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='llama3.1:70b')
          print(f"Using expensive model for analysis: {expensive_llm.model}")
        else:
          raise ValueError("Sorry, unable to set up LLM. Check llm provider (%s)" %provider)
      except Exception as e:
          raise ValueError("Sorry, unable to set up LLM. API keys for provider (%s)" %provider)
      return expensive_llm, embeddings

def get_cheap_llm(provider = None, timeout = None):
  if provider is None:
    provider = os.getenv("LLM_PROVIDER", "ollama")

  try:
    if provider in ["google", "openai"]:
      cheap_llm, _ = get_llm_and_embeddings(provider='google', timeout=timeout)
    else:
      cheap_llm, _ = get_llm_and_embeddings(
       provider=provider, timeout=timeout, llm_model_name='llama3.1:8b')
    # Handle different attribute names
    model_name = getattr(cheap_llm, 'model', None) or getattr(cheap_llm, 'model_name', 'unknown')
    print(f"Using cheap/fast model for summarization: {model_name}")

  except ValueError as e:
    print(e)
    raise ValueError("Sorry, unable to set up LLM (%s). Check API keys." %(
         provider))
  return cheap_llm


# Convenience functions for Ollama
def get_ollama_llm(
    model_name: str = None,
    base_url: str = None,
    temperature: float = 0.1,
    num_ctx: int = 8192,
):
    """Get just the Ollama LLM without embeddings."""
    from langchain_ollama import ChatOllama

    if model_name is None:
        model_name = os.getenv("OLLAMA_LLM_MODEL", "llama3.1:70b")
    if base_url is None:
        base_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

    return ChatOllama(
        model=model_name,
        base_url=base_url,
        temperature=temperature,
        num_ctx=num_ctx,
        num_predict=4096,
    )


def get_ollama_fast_llm(base_url: str = None):
    """Get a fast Ollama LLM for simple tasks."""
    return get_ollama_llm(
        model_name="llama3.1:8b",
        base_url=base_url,
        temperature=0.0,
        num_ctx=4096,
    )


def get_ollama_smart_llm(base_url: str = None):
    """Get a powerful Ollama LLM for complex reasoning."""
    return get_ollama_llm(
        model_name="llama3.1:70b",
        base_url=base_url,
        temperature=0.1,
        num_ctx=8192,
    )
