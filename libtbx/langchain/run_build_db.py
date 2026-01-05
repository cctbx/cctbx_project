"""
Script to build and persist the documentation vector database.
"""
from __future__ import division
import os
import shutil  # <--- Added for cleanup

from libtbx.langchain.core import get_llm_and_embeddings
# Import the chunker explicitly
from libtbx.langchain.rag.document_loader import load_all_docs_from_folder, _custom_chunker
from libtbx.langchain.rag import create_and_persist_db

def run(docs_folder_path_list=["./data_docs/"], db_dir=None,
        excluded_dirs=None, provider=None, timeout=300):


    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "ollama")

    # 1. CLEANUP: Delete existing DB to prevent duplicates
    if db_dir and os.path.exists(db_dir):
        print(f"Removing existing database at {db_dir} to rebuild fresh...")
        shutil.rmtree(db_dir)

    if not docs_folder_path_list:
        docs_folder_path_list = ["./data_docs/"]

    all_processed_files = []
    raw_docs = []

    # 1. Create a set to track unique file paths
    seen_files = set()

    for docs_folder_path in docs_folder_path_list:
        new_docs, new_files = load_all_docs_from_folder(
            docs_folder_path,
            excluded_dirs=excluded_dirs
        )

        # 2. Filter out duplicates
        unique_new_docs = []
        for doc in new_docs:
            source = doc.metadata.get('source')
            # Only add if we haven't seen this file path before
            if source and source not in seen_files:
                unique_new_docs.append(doc)
                seen_files.add(source)
            elif not source:
                 # Safety for docs without source metadata
                 unique_new_docs.append(doc)

        raw_docs.extend(unique_new_docs)

        # Same for the file list
        for f in new_files:
            if f not in all_processed_files:
                all_processed_files.append(f)

    print(f"Total unique documents loaded: {len(raw_docs)}")

    # 2. CHUNKING: Apply your smart chunker logic here
    print(f"Chunking {len(raw_docs)} documents...")
    chunked_docs = _custom_chunker(raw_docs)
    print(f"Created {len(chunked_docs)} chunks.")

    try:
        llm, embeddings = get_llm_and_embeddings(
            provider=provider, timeout=timeout)
    except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM with %s" % (provider))

    # Pass the CHUNKED docs, not the raw docs
    create_and_persist_db(chunked_docs, embeddings, db_dir)

    return all_processed_files

if __name__ == "__main__":
    run()

