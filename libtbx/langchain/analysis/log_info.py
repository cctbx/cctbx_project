"""
High-level log information extraction.

This module provides the get_log_info function that combines
summarization and processing into a single operation.
"""
from __future__ import absolute_import, division, print_function

import os
import asyncio

import openai
from google.api_core import exceptions as google_exceptions

from libtbx import group_args
from libtbx.langchain.analysis.summarizer import summarize_log_text
from libtbx.langchain.utils.text_processing import get_processed_log_dict


async def get_log_info(text, llm, embeddings, timeout: int = 120,
                       provider: str = None, program_name: str = None):
    """
    Summarizes a log file and extracts key information.

    This is the main entry point for log analysis, combining:
    - Log summarization (map-reduce)
    - Key information extraction (program name, etc.)

    Args:
        text: The log file content
        llm: Language model for summarization
        embeddings: Embeddings model (not currently used, kept for API compatibility)
        timeout: Timeout in seconds
        provider: 'google' or 'openai' or 'ollama'
        program_name: Explicit program name (e.g., 'phenix.refine') to ensure
                      correct summary template is used. If None, will be auto-detected.

    Returns:
        group_args with:
            - summary: The log summary text
            - summary_as_html: HTML version of the summary
            - processed_log_dict: Dict with extracted information
            - error: Error message if any

    Example:
        result = await get_log_info(log_text, llm, embeddings, program_name='phenix.refine')
        if result.error:
            print(f"Error: {result.error}")
        else:
            print(result.summary)
    """
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    # Debug: Log the program name being used
    if program_name:
        print(f"[LOG_INFO] Using explicit program_name: {program_name}")
    else:
        print(f"[LOG_INFO] No program_name provided - will auto-detect")

    try:
        log_summary_info = await summarize_log_text(
            text, llm, timeout=timeout, provider=provider,
            program_name=program_name)
        processed_log_dict = get_processed_log_dict(
            log_summary_info.log_summary)

        if log_summary_info.error:
            return group_args(
                group_args_type='error', error=log_summary_info.error)

        # Import save_as_html here to avoid circular imports
        from libtbx.langchain.run_analyze_log import save_as_html

        return group_args(
            group_args_type='log summary',
            summary=log_summary_info.log_summary,
            summary_as_html=save_as_html(log_summary_info.log_summary),
            processed_log_dict=processed_log_dict,
            error=None
        )

    except asyncio.TimeoutError:
        error_message = "Analysis timed out, try increasing timeout."
        print(error_message)
        return group_args(group_args_type='error', error=error_message)

    except openai.AuthenticationError:
        error_message = "OPENAI API key is invalid"
        return group_args(group_args_type='error', error=error_message)

    except Exception as e:
        # Inspect the exception's original cause to find the specific error
        original_cause = getattr(e, '__cause__', None)

        if isinstance(original_cause, google_exceptions.PermissionDenied):
            error_message = "Google AI API key is invalid"
        elif isinstance(original_cause, google_exceptions.ResourceExhausted):
            error_message = "Google AI API quota exceeded"
        else:
            # Handle any other general exception
            msg = str(e)
            if "API_KEY_INVALID" in msg or "API_KEY_IP_ADDRESS_BLOCKED" in msg:
                error_message = "Google API key is invalid"
            elif "free_tier_requests, limit: 0" in msg:
                error_message = "Google API key is not activated"
            elif "limit: 0" in msg:
                error_message = "Google AI API key has a zero quota"
            elif "request timed out" in msg:
                error_message = "Summarizing timed out."
            else:
                error_message = f"An unexpected error occurred during summarization: {e}"
            print("Summarize log failed")

        print(error_message)
        return group_args(group_args_type='error', error=error_message)

