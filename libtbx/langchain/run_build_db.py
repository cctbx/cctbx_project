"""
Script to build and persist the documentation vector database.
"""
from __future__ import division
from libtbx.langchain import langchain_tools as lct
import os

# In run_build_db.py

def run(docs_folder_path_list = ["./data_docs/"], db_dir = "./docs_db",
   excluded_dirs = None, provider = 'google', timeout = 300):
    """
    Initializes models and builds the documentation vector database.
    Returns a list of the file paths that were processed.
    """
    if provider == 'google':
      if not os.getenv("GOOGLE_API_KEY"):
        raise ValueError("GOOGLE_API_KEY environment variable not set.")
    elif provider == 'openai':
      if not os.getenv("OPENAI_API_KEY"):
        raise ValueError("OPENAI_API_KEY environment variable not set.")


    if not docs_folder_path_list:
      docs_folder_path_list = ["./data_docs/"]

    # Initialize a list for all processed file paths ---
    all_processed_files = []
    all_docs = []
    for docs_folder_path in docs_folder_path_list:
      # Unpack the two return values ---
      new_docs, new_files = lct.load_all_docs_from_folder(
          docs_folder_path,
          excluded_dirs = excluded_dirs
      )
      all_docs.extend(new_docs)
      all_processed_files.extend(new_files)

    try:
        llm, embeddings = lct.get_llm_and_embeddings(
            provider=provider, timeout=timeout)
    except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM with %s" %(provider))

    lct.create_and_persist_db(all_docs, embeddings, db_dir)

    # Return the list of processed files ---
    return all_processed_files

if __name__ == "__main__":
    run()
