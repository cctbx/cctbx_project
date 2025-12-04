"""
Summarizes a log file (usually Phenix) and suggests next steps in the
context of the (Phenix) documentation.  Uses a two-stage pipeline with
an advanced reranking RAG.
"""
from __future__ import division
import sys, os, json
import time

def run(file_name = None,
      log_as_text = None,
      output_file_path = None,
      db_dir: str = "./docs_db",
      timeout: int = 60,
      existing_summary: str  = None,
      existing_analysis: str  = None,
      display_results: bool = True,
      text_to_append_to_summary: str = None,
      text_to_append_to_analysis: str = None,
      summary_html_file_name = None,  # write to these files if supplied
      analysis_html_file_name = None, # write to these files if supplied
      provider: str = 'google',
      max_analyze_log_tries = 3,
      log_directory: str = None, # Where to save the history
      job_id: str = None,        # The ID for this specific run (e.g., '23')
      run_agent: bool = False,   # are we going to run the agent
      history_files: list = None,  # text files here
      history_simple_string: str = None # content from client
      ):

    # What we are going to return (Backwards Compatible Object)
    from libtbx import group_args
    log_info = group_args(group_args_type = 'log summary',
      summary = existing_summary,
      analysis = existing_analysis,
      log_text = None,
      processed_log_dict = None,
      error = None,
      next_move = None,     # NEW field
      history_record = None # NEW field
      )

    def create_history_record_dict(log_info, next_move=None):
        """Helper to create the dict for JSON serialization"""
        rec = {
            "job_id": job_id,
            "timestamp": time.time(),
            "file_name": file_name,
            "program": log_info.processed_log_dict.get('phenix_program', 'Unknown') if log_info.processed_log_dict else 'Unknown',
            "summary": log_info.summary,
            "analysis": log_info.analysis,
            "error": log_info.error,
        }
        if next_move:
            rec['next_move'] = next_move
        return rec

    # Get the text for the log file
    if not log_as_text:
      if not os.path.isfile(file_name):
        raise ValueError("Sorry, the file %s is missing" %(file_name))
      log_as_text = open(
         file_name, 'r', encoding='utf-8', errors='ignore').read()
      print("Summarizing the log file '%s'..." %(file_name))
    elif file_name:
      print("Summarizing the log file '%s'..." %(file_name))
    else:
      print("Summarizing supplied text as a log file")
      file_name = "Supplied text"

    # Decide where to write files
    if not output_file_path:
      import tempfile
      dd = tempfile.TemporaryDirectory()
      output_file_path = dd.name
      # Note: will be deleted when execution ends

    if (not log_info.summary) or (not log_info.analysis):
      # need to get the info
      if provider == 'google':
        if not os.getenv("GOOGLE_API_KEY"):
           log_info.error = "GOOGLE_API_KEY environment variable not set."
           return log_info
      elif provider == 'openai':
        if not os.getenv("OPENAI_API_KEY"):
           log_info.error = "OPENAI_API_KEY environment variable not set."
           return log_info
      if not os.getenv("COHERE_API_KEY"):
         log_info.error = "COHERE_API_KEY environment variable not set."
         return log_info

      # --- NEW: Automatically select DB directory based on provider ---
      if db_dir == "./docs_db":
        if provider == 'google':
            db_dir = "./docs_db_google"
            print(f"Provider is 'google', using default database: {db_dir}")
        elif provider == 'openai':
            db_dir = "./docs_db_openai"
            print(f"Provider is 'openai', using default database: {db_dir}")

      if not os.path.exists(db_dir):
          log_info.error = (
              f"Database not found for provider '{provider}'. "
              f"Expected at: {db_dir}\n"
              "Please run the database build script for this provider."
          )
          return log_info
      # ---------------------------------------------------------------

      try:
        from libtbx.langchain import langchain_tools as lct
        import asyncio
      except Exception as e:
        raise ValueError("Sorry, unable to analyze the file %s " %(file_name))

      # --- MODIFIED SECTION: Set up two LLMs ---
      print("Setting up LLMs...")
      try:
        # 1. the EXPENSIVE model ('gpt-5') for final analysis or gemini-2.5-pro
        if provider == "google":
          expensive_llm, embeddings = lct.get_llm_and_embeddings(
            provider=provider,
            timeout=timeout,
            llm_model_name='gemini-2.5-pro' # Your requested powerful model
          )
          print(f"Using expensive model for analysis: {expensive_llm.model}")
        else:
          expensive_llm, embeddings = lct.get_llm_and_embeddings(
            provider=provider,
            timeout=timeout,
            llm_model_name='gpt-5' # Your requested powerful model
          )
          print(f"Using expensive model for analysis: {expensive_llm.model_name}")

        # 2. Get the CHEAP & FAST Google model for summarization
        cheap_llm, _ = lct.get_llm_and_embeddings(
            provider='google',
            timeout=timeout
        )
        print(f"Using cheap/fast model for summarization: {cheap_llm.model}")

      except ValueError as e:
        print(e)
        raise ValueError("Sorry, unable to set up LLM. Check both GOOGLE_API_KEY and OPENAI_API_KEY.")
      # --- END OF MODIFIED SECTION ---


      # Summarize the log file (using the CHEAP & FAST Google LLM)
      print("Summarizing log file (using cheap Google model)...")
      result = asyncio.run(lct.get_log_info(
          log_as_text,
          cheap_llm,  # <-- Pass Google's cheap_llm
          embeddings, # <-- Pass OpenAI's embeddings (summarize doesn't use them)
          timeout = timeout,
          provider = 'google' # Tell the function which provider we're using
      ))

      if result.error or not result.summary: # failed
        print("Log file summary failed")
        log_info.error = result.error
        return log_info # Return object

      log_info = result
      log_info.analysis = None

    log_info.log_text = log_as_text

    # Put log summary in an html window
    if display_results:
      fn = summary_html_file_name
      if not fn:
        fn = os.path.join(output_file_path,'log_summary.html')
      text = log_info.summary
      if text_to_append_to_summary:
         text += text_to_append_to_summary
      save_as_html(text, file_name = fn,
       title = 'Summary of %s' %(file_name))
      print("Loading log summary at %s" %(fn))
      try:
        from phenix.command_line.doc import load_url
        load_url(fn)
      except Exception as e:
        # phenix is not available or no viewer.  Just skip
        print("Unable to load viewer")

    # Analyze the log summary in the context of the docs
    print("\nAnalyzing summary in context of documentation (using expensive model)...")
    if (db_dir and (not log_info.analysis)):
      ok = False
      from copy import deepcopy
      log_info_std = deepcopy(log_info)

      # --- MODIFIED RETRY LOGIC ---
      last_error = "" # Keep track of the last error
      for i in range(max_analyze_log_tries):
        log_info = deepcopy(log_info_std)
        result = asyncio.run(
          lct.analyze_log_summary(log_info, expensive_llm, embeddings,
          db_dir = db_dir, timeout = timeout))

        if (not result.error) and (result.analysis): # Success
          log_info.analysis = result.analysis
          ok = True
          break # Exit the loop on success
        else:
          # Failure, store the error and try again
          last_error = result.error or "Unknown analysis error"
          print(f"Analysis failed (Attempt {i+1}/{max_analyze_log_tries}): {last_error}")
          time.sleep(1) # Optional: short pause before retry

      if (not ok):
        # If we get here, all retries failed
        log_info.error = f"Analysis failed after {max_analyze_log_tries} attempts. Last error: {last_error}"
        log_info.analysis = ""
        print("Unable to carry out analysis of log file.")
        return log_info # Return object
      # --- END OF MODIFIED LOGIC ---

    # Put it in an html window
    if display_results:
      text = log_info.analysis
      if text_to_append_to_analysis:
         text += text_to_append_to_analysis
      fn = analysis_html_file_name
      if not fn:
        fn = os.path.join(output_file_path,'analysis.html')
      save_as_html(text, file_name = fn,
       title = 'Analysis of %s' %(file_name))
      print("Loading analysis at %s" %(fn))
      try:
        from phenix.command_line.doc import load_url
        load_url(fn)
      except Exception as e:
        # phenix is not available or no viewer.  Just skip
        print("Unable to load viewer")
    # Make sure viewer has enough time to load
    time.sleep(0.5)

    # Save History to JSON ---
    # Convert to dict for saving, but we keep log_info as object for return
    history_record_dict = create_history_record_dict(log_info)
    log_info.history_record = history_record_dict # Attach to object

    if log_directory and \
          log_info.summary and log_info.analysis:
        try:
            if not os.path.exists(log_directory):
                os.makedirs(log_directory)
                print(f"Created log directory: {log_directory}")

            # Create a filename. Use job_id if provided, else timestamp.
            if job_id:
                json_filename = f"job_{job_id}.json"
            else:
                json_filename = f"job_{int(time.time())}.json"

            json_path = os.path.join(log_directory, json_filename)

            with open(json_path, 'w') as f:
                json.dump(history_record_dict, f, indent=2)

            print(f"Saved run history to: {json_path}")

        except Exception as e:
            print(f"Warning: Failed to save history JSON: {e}")
    # ---------------------------------

    # --- AGENT LOGIC (The New Brain) ---
    if run_agent:
        print("\n--- Running Agent for Next Move ---")
        try:
            # 1. Ensure LLM is loaded (Fixes the "Missing Brain" crash)
            if 'expensive_llm' not in locals() or 'embeddings' not in locals():
                from libtbx.langchain import langchain_tools as lct
                if provider == "google":
                    expensive_llm, embeddings = lct.get_llm_and_embeddings(
                        provider=provider, timeout=timeout, llm_model_name='gemini-2.5-pro')
                else:
                    expensive_llm, embeddings = lct.get_llm_and_embeddings(
                        provider=provider, timeout=timeout, llm_model_name='gpt-5')

            # 2. Load Past History
            past_history = []

            # Method A: String from Client (Robust for Remote)
            if history_simple_string:
                from phenix.rest import simple_string_as_text
                try:
                    json_str = simple_string_as_text(history_simple_string)
                    past_history = json.loads(json_str)
                except Exception as e:
                    print(f"Error decoding history string: {e}")

            # Method B: Files (Local backup)
            elif history_files:
                print(f"Loading {len(history_files)} past history files...")
                for fpath in history_files:
                    if os.path.isfile(fpath):
                        try:
                            with open(fpath, 'r') as f:
                                past_history.append(json.load(f))
                        except: pass

            # 3. Append Current Run
            full_history = past_history + [history_record_dict]

            # 4. Run Agent
            from libtbx.langchain import langchain_tools as lct
            agent_result = asyncio.run(lct.generate_next_move(
                run_history=full_history,
                llm=expensive_llm,
                embeddings=embeddings,
                db_dir=db_dir,
                timeout=timeout*3
            ))


            if agent_result:
                next_move_dict = {
                    "program": agent_result.program,
                    "explanation": agent_result.explanation,
                    "strategy": agent_result.strategy,
                    "command": agent_result.command,
                    "process_log": getattr(agent_result, 'process_log', ""),
                    "error": agent_result.error
                }
                # Attach to object AND dict
                log_info.next_move = next_move_dict
                history_record_dict['next_move'] = next_move_dict
                log_info.history_record = history_record_dict

                if log_directory and json_path:
                     with open(json_path, 'w') as f:
                        json.dump(history_record_dict, f, indent=2)

                print("\n=== AGENT RECOMMENDATION ===")
                print(agent_result.command)
                print("============================\n")

        except Exception as e:
            print(f"Agent failed: {e}")

    return log_info # RETURN OBJECT (Backwards Compatible)

def save_as_html(markdown_string: str,
     title: str = "Summary", file_name: str = None):
    """
    Converts a Markdown string to an HTML file.
    """
    html_content = ""
    try:
        from markdown_it import MarkdownIt
        md = MarkdownIt("gfm-like")
        html_content = md.render(markdown_string)
        print("Using 'markdown-it-py' for rich HTML rendering.")
    except ImportError:
        print("Warning: 'markdown-it-py' not found. Falling back to plain text rendering.")
        html_content = f"<pre>{markdown_string}</pre>"

    html_with_style = f"""
    <html>
    <head>
        <title>{title}</title>
        <style>
            body {{ font-family: sans-serif; line-height: 1.6; padding: 2em; max-width: 800px; margin: auto; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
            th, td {{ border: 1px solid #dddddd; text-align: left; padding: 8px; }}
            th {{ background-color: #f2f2f2; }}
            code {{ background-color: #eee; padding: 2px 4px; border-radius: 3px; }}
            pre {{ background-color: #f6f8fa; padding: 16px; border-radius: 6px; overflow: auto; white-space: pre-wrap; word-wrap: break-word;}}
        </style>
    </head>
    <body>
        <h2></b>{title}</b></h2>
        {html_content}
    </body>
    </html>
    """

    if file_name:
      with open(file_name, "w", encoding='utf-8') as f:
        f.write(html_with_style)
      print(f"\nSaved formatted output to: {file_name}")

    return html_with_style

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze Phenix Log or Error")

    # Standard file input
    parser.add_argument("file_name", nargs="?", help="Path to log file")

    # NEW: Allow raw error text input
    parser.add_argument("--error", type=str, help="Paste the error message string here")

    parser.add_argument("--db_dir", default="./docs_db", help="Database directory")
    parser.add_argument("--log_dir", default=None, help="Directory to save history JSON")
    parser.add_argument("--job_id", default=None, help="Job ID for history")
    parser.add_argument("--provider", default="google", help="LLM Provider")

    args = parser.parse_args()

    log_text = None
    file_label = "Error Message"

    # Logic: Did user provide a file or an error string?
    if args.error:
        log_text = f"COMMAND FAILURE REPORT:\n{args.error}"
        print(f"Analyzing error message...")
    elif args.file_name:
        file_label = args.file_name
        # The run function handles file reading, so we leave log_text None
    else:
        print("Usage: phenix.python run_analyze_log.py <file> OR --error 'Error text'")
        sys.exit(1)

    answer = run(
        file_name=file_label if not log_text else None, # Pass None filename if using text
        log_as_text=log_text,                           # Pass the raw error text
        db_dir=args.db_dir,
        log_directory=args.log_dir,
        job_id=args.job_id,
        provider=args.provider
    )

