# run_inspect_db.py
"""
Runner script to inspect the contents of the persisted vector database.
"""
from __future__ import division
try:
    from libtbx.langchain import langchain_tools as lct
except ImportError:
    import langchain_tools as lct
from collections import Counter
import sys

def run(db_dir = "./docs_db", provider = "google"):
    """
    Loads the database and prints a summary of its contents.
    """
    try:
        # FIX: Use get_llm_and_embeddings instead of get_embeddings
        # We ignore the LLM (_) as we only need embeddings here
        print(f"Loading embeddings for provider: {provider}...")
        _, embeddings = lct.get_llm_and_embeddings(provider=provider)
        
        vectorstore = lct.load_persistent_db(embeddings, db_dir)

        # Directly access the collection to get metadata
        collection_data = vectorstore._collection.get(include=["metadatas"])

        if not collection_data or not collection_data.get('ids'):
            print("The database is empty or could not be loaded correctly.")
            return

        total_chunks = len(collection_data['ids'])
        print(f"\nDatabase contains a total of {total_chunks} document chunks.")

        print("\n--- Source Files in Database ---")
        # Handle cases where metadata might be None
        metadatas = collection_data.get('metadatas', [])
        if not metadatas:
             print("No metadata found in collection.")
             return

        source_paths = [m.get('source', 'Unknown') for m in metadatas if m]

        if not source_paths:
            print("No source file information found in the database metadata.")
            return

        # Count occurrences of each source file
        source_counts = Counter(source_paths)

        for source_file, count in sorted(source_counts.items()):
            print(f"  - {source_file} ({count} chunks)")

        print("--------------------------------\n")

    except FileNotFoundError as e:
        print(f"ERROR: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Usage: phenix.python run_inspect_db.py [db_dir] [provider]
    
    db_dir = "/net/cci-gpu-01/raid1/scratch1/terwill/build_dir/modules/phenix/phenix/phenix_ai/docs_db_google"
    provider = "google"
    
    if len(sys.argv) > 1:
        db_dir = sys.argv[1]
    if len(sys.argv) > 2:
        provider = sys.argv[2]
        
    run(db_dir=db_dir, provider=provider)

