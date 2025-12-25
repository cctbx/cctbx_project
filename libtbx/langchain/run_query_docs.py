"""
Ask a direct question to the Phenix documentation RAG.
Uses the advanced reranking retrieval chain.
"""
from __future__ import division
import sys
import time
import os

from libtbx.langchain.core import get_llm_and_embeddings
from libtbx.langchain.utils import query_docs
from libtbx.langchain.run_analyze_log import save_as_html


def run(query_text=None, output_file_path=None, db_dir=None,
        timeout=60, provider=None):
    """
    Loads the reranking RAG and queries it with a user-provided question.
    """

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    if provider == 'google':
        if not os.getenv("GOOGLE_API_KEY"):
            raise ValueError("GOOGLE_API_KEY environment variable not set.")
    elif provider == 'openai':
        if not os.getenv("OPENAI_API_KEY"):
            raise ValueError("OPENAI_API_KEY environment variable not set.")
    elif provider == 'ollama':
        pass # ok
    else:
        raise ValueError(f"Unknown provider: {provider}")

    if len(sys.argv) < 2:
        print("Usage: python3 run_query_docs.py \"<Your question here>\"")
        sys.exit(1)

    # Set up the LLM
    try:
        llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout)
    except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM with %s" % (provider))

    print("\nQuerying the Phenix documentation with '%s'..." % (query_text))
    # Answer the question
    answer = query_docs(query_text=query_text, db_dir=db_dir,
                        timeout=timeout, provider=provider)

    if not answer:  # no result
        print("No answer obtained")
        return answer

    # Put it in an html window

    # Decide where to write files
    if not output_file_path:
        import tempfile
        dd = tempfile.TemporaryDirectory()
        output_file_path = dd.name

    fn = os.path.join(output_file_path, 'query.html')
    save_as_html(answer, file_name=fn,
                 title='Answer to: %s' % (query_text))
    print("Loading answer summary at %s" % (fn))
    try:
        from phenix.command_line.doc import load_url
        load_url(fn)
    except Exception as e:
        # phenix is not available or no viewer. Just skip
        print("Unable to load viewer")
    time.sleep(1)

    return answer


if __name__ == "__main__":
    query_text = " ".join(sys.argv[1:])
    run(query_text=query_text)
