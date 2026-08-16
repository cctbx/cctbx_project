"""
Log analysis using RAG (Retrieval-Augmented Generation).

This module analyzes log summaries in the context of Phenix documentation
to provide insights and suggest next steps.

Usage:
  from libtbx.langchain.analysis.analyzer import analyze_log_summary

  result = await analyze_log_summary(log_info, llm, embeddings, db_dir)
  if result.error:
    print(f"Error: {result.error}")
  else:
    print(result.analysis)
"""
from __future__ import absolute_import, division, print_function

from concurrent.futures import TimeoutError

from langchain_core.prompts import PromptTemplate
from google.api_core import exceptions as google_exceptions

from libtbx import group_args
# v118.G2: retriever functions moved to function-body scope below.
# rag.retriever imports langchain_chroma which has known protobuf
# version conflicts in some envs.  Eager top-level import here would
# crash this module entirely (including prompt-only functions like
# get_log_analysis_prompt), so the imports are deferred to
# analyze_log_summary() where they are actually used.
# See docs/DEVELOPER_GUIDE.md "Optional dependency handling".


# =============================================================================
# Prompt Templates
# =============================================================================

def get_log_analysis_prompt() -> PromptTemplate:
  """Returns the prompt for analyzing a log summary."""
  template = (
    "You are expert in crystallography and cryo-EM."
    "You are a Phenix power-user. Your task is to analyze "
    "a program summary in the context of the provided documentation "
    "and research papers.\n\n"
    "Based on the data provided at the end, please perform the following "
    " analysis. Consider the events of the log summary in the broader "
    "context of a "
    "typical Phenix structure determination workflow as described in "
    "the documentation and papers, but do not describe this context. "
    "The analysis should include the following for the run "
    "described in the log summary:\n\n"

    "1. Evaluate whether the run described in the summary was useful. "
      "Provide a short summary of the usefulness of the run."
      "Then List reported metrics and expected values of these metrics and "
      "consider the goals of the program. Note any warnings, errors,"
      "or advisories obtained. "
      "If the results of the run indicate a low confidence solution,"
      "multiple solutions, or no solution, clearly state this observation.\n\n"

    "2. Considering whether the data are from crystallography or"
    "cryo-EM and considering the normal sequence of Phenix tool use"
    "for that type of data, suggest "
    "three concrete next steps in structure determination using Phenix "
    "or graphical tools such as Coot or Isolde that I should take, "
    "justifying each suggestion with information from the provided "
    "documentation and papers. "
    "If the results of the run indicate a low confidence solution,"
    "multiple solutions, or no solution, then instead suggest concrete"
    "steps for backtracking and figuring out what went wrong."
    "Do not include any information on CryoFit, ShelxD, Parrot, "
    "or CCP4 tools "
    "unless there is a specific question about them. "
    "If appropriate, include validation as one of your next steps. "
    "Name the tools that are to be used, along with their inputs "
    " and outputs and what they do. "
    "Do not suggest depositing the model. "
    "Do not suggest analyzing the biological relevance. \n\n"

    "3.List the inputs and briefly describe what was done."
    "Report whether the data are from crystallography (X-ray or neutron)"
    " or from cryo-EM\n\n"

    "4. List the key output files from this run, along with the values "
    "of any available metrics describing their utilities. "
    "If no metrics are available, do not provide any. \n\n"

    "**Please note: no offers of help*** Do not offer to help the "
    "user with additional analyses and do not mention that you are"
    "not to offer to help."
    "\n\n---BEGIN DATA FOR ANALYSIS---\n"
    "Documentation Context:\n{context}\n\n"
    "Log File Summary:\n{log_summary}"
  )
  return PromptTemplate(
    template=template,
    input_variables=["context", "log_summary"])


# =============================================================================
# Main Analysis Function
# =============================================================================

async def analyze_log_summary(log_info, llm, embeddings,
                db_dir: str = "./docs_db",
                timeout: int = 60):
  """
  Analyzes a log summary using RAG with Phenix documentation.

  Args:
    log_info: Object with 'summary' and 'processed_log_dict' attributes
    llm: Language model for analysis
    embeddings: Embedding model for retrieval
    db_dir: Path to the documentation vector database
    timeout: Timeout in seconds

  Returns:
    group_args with:
      - group_args_type: 'answer'
      - analysis: The analysis text (or None if failed)
      - error: Error message (or None if successful)

  Example:
    result = await analyze_log_summary(log_info, llm, embeddings)
    if result.error:
      print(f"Failed: {result.error}")
    else:
      print(result.analysis)
  """
  try:
    # Retrieval removed.  retriever.py:275 is
    #     "context": lambda x: retriever.invoke(x["input"])
    # so retrieval keyed on the query string ALONE.  That string was
    # built from processed_log_dict['phenix_program'], which the
    # summary scrape returns as '' on the format the model actually
    # emits, and ['summary'], which is assigned "" and never
    # reassigned.  Measured across four captured runs covering four
    # programs, the query was
    #     'Here is a summary of the  log file:\n\n ' + fixed boilerplate
    # byte-identical for every log.  Retrieval selected documentation
    # against a constant string; removing it costs nothing because it
    # contributed nothing.
    #
    # `embeddings` and `db_dir` are now unused here and kept only so
    # callers do not have to change.
    payload = getattr(log_info, 'analysis_payload', None)
    if not payload:
      return group_args(
        group_args_type='answer',
        analysis=None,
        error="no analysis payload was built for this log",
      )

    final_analysis = llm.invoke(payload)
    content = getattr(final_analysis, 'content', final_analysis)

    # Verifier, SHADOW MODE: log only, surface nothing.
    #
    # Phase 0a found zero fabricated figures in twenty reports, so this
    # guards a problem we do not have and its whole value is rarity.
    # Measured false-positive rate: 10% in-sample on the 40 study
    # reports, 23% held-out on 13 shipped reports; detection unaffected
    # (one fabricated figure injected into each of twenty reports is
    # caught 20 of 20).  The entire residual is one category -- the
    # R-free/R-work gap, a difference of two log values.  Do not expose
    # flags to users until that is decided.
    try:
      from libtbx.langchain.analysis.report_verifier import (
        check_numbers, check_program_names, summarise)
      from libtbx.langchain.analysis.program_identity import load_registry
      flags = check_numbers(content, getattr(log_info, 'log_text', '') or '')
      blocked = check_program_names(content, set(load_registry()))
      note = summarise(flags, blocked)
      getattr(log_info, 'debug_log', []).append("VERIFIER (shadow): %s"
                                                % note)
      for flag in flags[:10]:
        getattr(log_info, 'debug_log', []).append("  flag: %s" % flag)
    except Exception as verifier_error:              # noqa: BLE001
      # The verifier must never break a report it is only observing.
      getattr(log_info, 'debug_log', []).append(
        "VERIFIER unavailable: %s" % verifier_error)

    return group_args(
      group_args_type='answer',
      analysis=content,
      error=None
    )

  # --- EXCEPTION HANDLERS ---
  except (TimeoutError, google_exceptions.DeadlineExceeded) as e:
    error_message = (
      "Network timeout with Google. "
      "You might try increasing the timeout in Preferences or in "
      f"AnalyzeLog (currently {timeout} sec)."
    )
    print(error_message)
    return group_args(
      group_args_type='answer',
      analysis=None,
      error=error_message
    )

  except google_exceptions.ResourceExhausted as e:
    # Log full details for debugging, return clean user message.
    # No period before the action text — contains_sorry() in rest/__init__.py
    # truncates at the first '.' or '\n' and appends "Please wait a minute"
    # unless "please" already appears in the first chunk.
    print(f"Google API quota exceeded. Full error: {e}")
    error_message = (
      "Google API quota exceeded, please try another provider "
      "(eg provider=openai) or wait for quota reset"
    )
    print(error_message)
    return group_args(
      group_args_type='answer',
      analysis=None,
      error=error_message
    )

  except Exception as e:
    e_str = str(e)
    e_lower = e_str.lower()
    # Langchain wraps google quota errors as generic exceptions with
    # "RESOURCE_EXHAUSTED" in the message — the ResourceExhausted branch
    # above only catches the raw google exception.
    _is_quota = (
      "resource_exhausted" in e_lower or
      ("quota" in e_lower and "exceeded" in e_lower)
    )
    if _is_quota:
      # Log full details for debugging, return clean user message.
      print(f"Google API quota exceeded. Full error: {e_str}")
      error_message = (
        "Google API quota exceeded, please try another provider "
        "(eg provider=openai) or wait for quota reset"
      )
    else:
      error_message = (
        "Reranking failed - try again in a"
        " couple minutes..." + e_str)
    print(error_message)
    return group_args(
      group_args_type='answer',
      analysis=None,
      error=error_message
    )
