"""
Documentation query utilities.

This module provides the high-level query_docs function
for querying the Phenix documentation RAG.
"""
from __future__ import absolute_import, division, print_function

import time
import os

import openai
from google.api_core import exceptions as google_exceptions
from langchain_google_genai._common import GoogleGenerativeAIError

from libtbx.langchain.core import get_llm_and_embeddings
from libtbx.langchain.rag import (
    load_persistent_db,
    create_reranking_retriever,
    create_reranking_rag_chain,
)
from libtbx.langchain.knowledge import get_docs_query_prompt


# Global variable to track the last query time for throttling
_last_query_time = 0


def query_docs(query_text, llm=None, embeddings=None,
               db_dir: str = "./docs_db",
               timeout: int = 60,
               max_attempts: int = 5,
               use_throttling: bool = False,
               provider: str = None):
    """
    Query Phenix docs with query_text, with automatic retries on rate limits.

    Args:
        query_text: The question to ask
        llm: Language model (optional, will create if not provided)
        embeddings: Embeddings model (optional, will create if not provided)
        db_dir: Path to the documentation database
        timeout: Timeout in seconds
        max_attempts: Maximum retry attempts
        use_throttling: Whether to throttle queries
        provider: 'google' or 'openai' or 'ollama'

    Returns:
        str: The answer, or None if query failed

    Example:
        answer = query_docs("How do I run molecular replacement?")
        print(answer)
    """
    global _last_query_time

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")
    if use_throttling:
        min_interval = 4
        elapsed_time = time.time() - _last_query_time
        if elapsed_time < min_interval:
            wait_time = min_interval - elapsed_time
            print(f"Throttling: Waiting {wait_time:.1f} seconds ...")
            time.sleep(wait_time)
        _last_query_time = time.time()

    if not llm:
        try:
            llm, embeddings = get_llm_and_embeddings(
                provider=provider, timeout=timeout)
        except ValueError as e:
            print(e)
            raise ValueError("Sorry, unable to set up LLM with %s" % (provider))

    query_text += """Consider this question in the context of the process
        of structure determinaion in Phenix. Focus on using Phenix tools,
        but include the use of Coot or Isolde if appropriate. Name the tools
        that are to be used, along with their inputs and outputs and what
        they do."""

    vectorstore = load_persistent_db(embeddings, db_dir=db_dir)
    retriever = create_reranking_retriever(vectorstore, llm, timeout=timeout)
    prompt = get_docs_query_prompt()
    rag_chain = create_reranking_rag_chain(retriever, llm, prompt)

    backoff_time = 2.0
    for attempt in range(max_attempts):
        try:
            response = rag_chain.invoke(query_text)
            return response.content

        except openai.RateLimitError:
            if attempt < max_attempts - 1:
                print(f"OpenAI rate limit exceeded. Waiting {backoff_time:.1f} sec...")
                time.sleep(backoff_time)
                backoff_time *= 2
            else:
                print("OpenAI API rate limit exceeded after multiple retries.")
                return None
        except openai.AuthenticationError:
            print("OpenAI API key is invalid")
            return None

        except (google_exceptions.ResourceExhausted, GoogleGenerativeAIError) as e:
            error_text = str(e).lower()
            if "you exceeded your current quota" in error_text:
                if "limit: 0" in error_text:
                    print("Google AI API has a zero quota. Check plan.")
                    return None
                else:
                    print("Google AI API quota exceeded. Check plan.")
                    return None

            if "429" in error_text or "rate limit" in error_text:
                if attempt < max_attempts - 1:
                    print(f"Rate limit exceeded. Waiting for {backoff_time:.1f} sec...")
                    time.sleep(backoff_time)
                    backoff_time *= 2
                else:
                    print("Google AI API quota exceeded after retries.")
                    return None
            else:
                print(f"An unexpected Google AI error occurred: {e}")
                return None

        except Exception as e:
            print(f"An unexpected error occurred during query: {e}")
            return None

    print("Query docs failed...")
    return None
