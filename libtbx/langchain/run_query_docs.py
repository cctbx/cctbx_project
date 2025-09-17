"""
Ask a direct question to the Phenix documentation RAG.
Uses the advanced reranking retrieval chain.
"""
from __future__ import division
from libtbx.langchain import langchain_tools as lct
from libtbx.langchain.run_analyze_log import save_as_html

import sys
import os

def run(query_text = None, output_file_path = None, db_dir = None,
     timeout= 60):
    """
    Loads the reranking RAG and queries it with a user-provided question.
    """
    if not os.getenv("GOOGLE_API_KEY"):
        raise ValueError("GOOGLE_API_KEY environment variable not set.")

    if len(sys.argv) < 2:
        print("Usage: python3 run_query_docs.py \"<Your question here>\"")
        sys.exit(1)


    # Set up the LLM
    llm = lct.get_llm()
    embeddings = lct.get_embeddings()

    print("\nQuerying the Phenix documentation with '%s'..." %(query_text))
    # Answer the question
    answer = lct.query_docs(query_text = query_text, db_dir = db_dir,
      timeout = timeout)

    # Put it in an html window

    # Decide where to write files
    if not output_file_path:
      output_file_path = '/var/tmp'

    fn = os.path.join(output_file_path,'query.html')
    save_as_html(answer, file_name = fn,
       title = 'Answer to: %s' %(query_text))
    print("Loading answer summary at %s" %(fn))
    try:
      from phenix.command_line.doc import load_url
      load_url(fn)
    except Exception as e:
      # phenix is not available or no viewer.  Just skip
      print("Unable to load viewer...see text in the file: '%s'" %(fn))

    return answer


if __name__ == "__main__":
    query_text = " ".join(sys.argv[1:])
    run(query_text = query_text)
