"""
Log analysis using RAG (Retrieval-Augmented Generation).

This module analyzes log summaries in the context of Phenix documentation
to provide insights and suggest next steps.

Usage:
    from libtbx.langchain.analysis import analyze_log_summary

    result = await analyze_log_summary(log_info, llm, embeddings, db_dir)
    if result.error:
        print(f"Error: {result.error}")
    else:
        print(result.analysis)
"""
from __future__ import absolute_import, division, print_function

from concurrent.futures import TimeoutError

from langchain_core.prompts import PromptTemplate
from cohere.core.api_error import ApiError as CohereApiError
from google.api_core import exceptions as google_exceptions

from libtbx import group_args
from libtbx.langchain.rag import (
    load_persistent_db,
    create_reranking_retriever,
    create_log_analysis_chain,
)


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
    return PromptTemplate(template=template, input_variables=["context", "log_summary"])


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
        vectorstore = load_persistent_db(embeddings, db_dir=db_dir)
        analysis_prompt = get_log_analysis_prompt()
        retriever = create_reranking_retriever(vectorstore, llm, timeout=timeout)
        analysis_rag_chain = create_log_analysis_chain(retriever, llm, analysis_prompt)

        retriever_query = (
            "Here is a summary of the %s log file:\n\n " % (
                log_info.processed_log_dict['phenix_program']) +
            log_info.processed_log_dict['summary'] +
            "\n\nConsidering whether the input data are from crystallography "
            "(X-ray or neutron) or from cryo-EM, and considering the "
            "the normal procedure for structure determination "
            "in Phenix, what are the next steps that I should carry out?" +
            "Consider this question in the context of the process of structure "
            "determination in Phenix. Focus on using Phenix tools, but include "
            "the use of Coot or Isolde if appropriate. Name the tools that are "
            "to be used, along with their inputs and outputs and what they do."
        )

        final_analysis = analysis_rag_chain.invoke({
            "input": retriever_query,
            "log_summary": log_info.summary,
        })

        return group_args(
            group_args_type='answer',
            analysis=final_analysis.content,
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
        error_message = (
            "ERROR: Google AI API quota exceeded. "
            f"Details: {e}"
        )
        print(error_message)
        return group_args(
            group_args_type='answer',
            analysis=None,
            error=error_message
        )

    except CohereApiError as e:
        if hasattr(e, 'http_status') and e.http_status == 401:
            error_message = "Invalid Cohere API key.\n"
        elif str(e).find("invalid api token") > -1:
            error_message = "Cohere API key is invalid"
        else:
            error_message = f"ERROR: A Cohere API error occurred. Details: {e}"

        print(error_message)
        return group_args(
            group_args_type='answer',
            analysis=None,
            error=error_message
        )

    except Exception as e:
        error_message = "Reranking failed - try again in a couple minutes..." + str(e)
        print(error_message)
        return group_args(
            group_args_type='answer',
            analysis=None,
            error=error_message
        )

