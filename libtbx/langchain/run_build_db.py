"""
Script to build and persist the documentation vector database.
"""
from __future__ import division
from libtbx.langchain import langchain_tools as lct
import os

def run(docs_folder_path_list = "./data_docs/", db_dir = "./docs_db",
   excluded_dirs = None):
    """
    Initializes models and builds the documentation vector database.
    """

    # Ensure GOOGLE_API_KEY is set before proceeding
    if not os.getenv("GOOGLE_API_KEY"):
        raise ValueError("GOOGLE_API_KEY environment variable not set.")

    if not docs_folder_path_list:
      docs_folder_path_list = ["./data_docs/"]

    # Load all documentation files
    all_docs = []
    for docs_folder_path in docs_folder_path_list:
      all_docs += lct.load_all_docs_from_folder(docs_folder_path,
        excluded_dirs = excluded_dirs)

    # Initialize the embedding model
    embeddings = lct.get_embeddings()

    # Create and persist the database
    lct.create_and_persist_db(all_docs, embeddings, db_dir)

if __name__ == "__main__":
    run()
