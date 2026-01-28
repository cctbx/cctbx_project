"""
LLM and Embeddings setup for the Phenix Crystallography Agent.

This module handles initialization of Language Models and Embedding models
from different providers (Ollama, Google, OpenAI).

Usage:
    from libtbx.langchain.core.llm import get_llm_and_embeddings

    # Use local Ollama (recommended)
    llm, embeddings = get_llm_and_embeddings(provider='ollama')

    # Or with Google
    llm, embeddings = get_llm_and_embeddings(provider='google')
"""
from __future__ import absolute_import, division, print_function

import os


# Global verbosity setting for LLM module
# Can be set via set_llm_verbosity() or LLM_VERBOSITY environment variable
_llm_verbosity = None


def set_llm_verbosity(level):
    """Set the verbosity level for LLM-related messages.

    Args:
        level: One of 'quiet', 'normal', 'verbose', 'debug'
    """
    global _llm_verbosity
    _llm_verbosity = level


def get_llm_verbosity():
    """Get the current verbosity level."""
    global _llm_verbosity
    if _llm_verbosity is not None:
        return _llm_verbosity
    return os.getenv("LLM_VERBOSITY", "normal")


def _llm_log(msg, level='normal'):
    """Print message if verbosity level permits."""
    levels = ['quiet', 'normal', 'verbose', 'debug']
    current = get_llm_verbosity()
    current_idx = levels.index(current) if current in levels else 1
    msg_idx = levels.index(level) if level in levels else 1
    if msg_idx <= current_idx:
        print(msg)


def get_llm_and_embeddings(
    provider: str = None,
    llm_model_name: str = None,
    embedding_model_name: str = None,
    temperature: float = 0.0,
    timeout: int = 120,
    batch_size: int = 100,
    ollama_base_url: str = None,
    json_mode: bool = False,
    num_ctx: int = 8192,
    seed: int = 42,
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
        json_mode: If True, force JSON output (for structured responses).
                   Only applies to Ollama provider.
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

        # Build kwargs conditionally
        llm_kwargs = {
            "model": llm_model_name,
            "base_url": ollama_base_url,
            "temperature": temperature,
            "num_ctx": num_ctx,
            "num_predict": 4096,
            "seed": seed,
            "think": False,
        }
        if json_mode:
            llm_kwargs["format"] = "json"

        llm = ChatOllama(**llm_kwargs)
        embeddings = OllamaEmbeddings(
            model=embedding_model_name,
            base_url=ollama_base_url,
        )
        mode_str = " (JSON mode)" if json_mode else ""
        _llm_log(f"Using Ollama at {ollama_base_url}{mode_str}", level='verbose')
        _llm_log(f"  LLM: {llm_model_name}, Embeddings: {embedding_model_name}", level='verbose')


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
        _llm_log(f"Using Google Gemini: {llm_model_name}", level='verbose')

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
        _llm_log(f"Using OpenAI: {llm_model_name}", level='verbose')

    else:
        raise ValueError(f"Unsupported provider: '{provider}'. Choose 'ollama', 'google', or 'openai'.")

    return llm, embeddings

# Get expensive or cheap LLMs

def get_expensive_llm(provider = None, timeout = None, json_mode=False):
      try:
        if provider == "google":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='gemini-2.5-pro')
          _llm_log(f"Using expensive model for analysis: {expensive_llm.model}", level='verbose')
        elif provider == "openai":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='gpt-5')
          _llm_log(f"Using expensive model for analysis: {expensive_llm.model_name}", level='verbose')
        elif provider == "ollama":
          expensive_llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout, llm_model_name='qwen3:32b', #'llama3.1:70b',# 'qwen2.5:72b',# 'llama3.1:405b', #'qwen2.5:72b', #llm_model_name='llama3.1:70b',
             json_mode=json_mode)
          _llm_log(f"Using expensive model for analysis: {expensive_llm.model}", level='verbose')
        else:
          raise ValueError("Sorry, unable to set up LLM. Check llm provider (%s)" %provider)
      except Exception as e:
          raise ValueError("Sorry, unable to set up LLM. API keys for provider (%s)" %provider)
      return expensive_llm, embeddings

def get_cheap_llm(provider=None, timeout=None):
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    try:
        if provider in ["google", "openai"]:
            cheap_llm, _ = get_llm_and_embeddings(provider=provider, timeout=timeout)
        else:
            cheap_llm, _ = get_llm_and_embeddings(
                provider=provider, timeout=timeout, llm_model_name='qwen2.5:7b') #'qwen2.5:14b') #'qwen2.5:32b') #'qwen3:32b')# 'qwen2.5:7b') # qwen2.5:72b') #'llama3.1:8b')

        # Handle different attribute names across providers
        model_name = getattr(cheap_llm, 'model', None) or getattr(cheap_llm, 'model_name', 'unknown')
        _llm_log(f"Using cheap/fast model for summarization: {model_name}", level='verbose')

        return cheap_llm

    except ValueError as e:
        _llm_log(str(e), level='quiet')
        raise ValueError("Sorry, unable to set up LLM (%s). Check API keys." % provider)
