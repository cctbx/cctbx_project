# run_inspect_db.py
"""
Runner script to inspect the contents of the persisted vector database.
"""
import langchain_tools as lct
from collections import Counter

def run(db_dir = "./docs_db"):
    """
    Loads the database and prints a summary of its contents.
    """
    try:
        embeddings = lct.get_embeddings()
        vectorstore = lct.load_persistent_db(embeddings, db_dir)

        # Directly access the collection to get metadata
        collection_data = vectorstore._collection.get(include=["metadatas"])

        if not collection_data or not collection_data.get('ids'):
            print("The database is empty or could not be loaded correctly.")
            return

        total_chunks = len(collection_data['ids'])
        print(f"\nDatabase contains a total of {total_chunks} document chunks.")

        print("\n--- Source Files in Database ---")
        source_paths = [metadata['source'] for metadata in collection_data['metadatas'] if 'source' in metadata]

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
    import sys
    if len(sys.argv) > 1:
      run(db_dir=sys.argv[1])
    else:
      run()
