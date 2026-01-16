from __future__ import division
import asyncio
import sys
import os
import glob
import json
import time

# Force using local library if needed
sys.path.insert(0, ".") 
import langchain_tools as lct

def get_run_history(log_directory, max_history=5):
    """
    Reads JSON files from the log_directory, sorts them by time,
    and returns the last N entries.
    """
    if not os.path.exists(log_directory):
        print(f"Error: Log directory '{log_directory}' does not exist.")
        return []

    # Find all json files matching the pattern job_*.json
    search_path = os.path.join(log_directory, "job_*.json")
    files = glob.glob(search_path)
    
    if not files:
        print(f"No history files found in {log_directory}")
        return []

    # Sort files by modification time (Oldest -> Newest)
    # This ensures the Agent sees the timeline correctly.
    files.sort(key=os.path.getmtime)

    # Load the data
    history = []
    print(f"Found {len(files)} history files.")
    
    for fpath in files:
        try:
            with open(fpath, 'r') as f:
                data = json.load(f)
                history.append(data)
        except Exception as e:
            print(f"Skipping corrupt file {fpath}: {e}")

    # Slice to keep only the last N
    if len(history) > max_history:
        print(f"Trimming history to last {max_history} runs (from {len(history)} total).")
        history = history[-max_history:]
    
    return history
def run(log_directory=".", max_history=5, db_dir = './docs_db'):
    # --- 1. Setup ---
    # Increase timeout to 180 seconds (3 minutes) to prevent 504 errors
    TIMEOUT = 180 
    
    try:
        print("Initializing Agent (Super-Analyst)...")
        # Use the PRO model for high-level reasoning
        llm, embeddings = lct.get_llm_and_embeddings(
            provider='google',
            llm_model_name='gemini-2.5-pro',
            timeout=TIMEOUT  # <--- UPDATE 1: Set LLM timeout here
        )
    except Exception as e:
        print(f"Setup failed: {e}")
        return

    # --- 2. Load Real History ---
    print(f"Reading history from: {log_directory}")
    history = get_run_history(log_directory, max_history)

    if not history:
        print("Cannot proceed without history.")
        return

    # --- 3. Run Super-Analysis ---
    print("\n--- Starting Super-Analysis ---")
    
    # We pass db_dir explicitly to ensure it finds your new database
    result = asyncio.run(lct.generate_next_move(
        run_history=history,
        llm=llm,
        embeddings=embeddings,
        db_dir=db_dir,
        timeout=TIMEOUT # <--- UPDATE 2: Pass timeout to the function
    ))

    # --- 4. Output Results ---
    if result.error:
        print(f"\nERROR: {result.error}")
    else:
        print("\n" + "="*40)
        print("RECOMMENDED ACTION")
        print("="*40)
        print(f"**Program:** {result.program}")
        print(f"**Explanation:**\n{result.explanation}\n")
        print(f"**Strategy Details:**\n{result.strategy}\n")
        print("-" * 20)
        print("GENERATED COMMAND:")
        print(f"{result.command}")
        print("="*40)

if __name__ == "__main__":
    # Simple CLI parsing
    import argparse
    parser = argparse.ArgumentParser(description="Phenix Agent")
    parser.add_argument("--log_dir", required=True, help="Directory containing job_*.json files")
    parser.add_argument("--max_history", type=int, default=5, help="Number of past jobs to consider")
    
    parser.add_argument("--db_dir", type=str, default='./docs_db', help="database_directory")
    
    args = parser.parse_args()
    
    run(log_directory=args.log_dir, max_history=args.max_history,
            db_dir = args.db_dir)
