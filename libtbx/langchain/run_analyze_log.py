"""
Summarizes a log file (usually Phenix) and suggests next steps in the
context of the (Phenix) documentation.  Uses a two-stage pipeline with
an advanced reranking RAG.
"""
from __future__ import division
import sys, os, json
import time
import asyncio

def run(file_name = None,
      log_as_text = None,
      output_file_path = None,
      db_dir: str = None,
      timeout: int = 60,
      existing_summary: str  = None,
      existing_analysis: str  = None,
      display_results: bool = True,
      text_to_append_to_summary: str = None,
      text_to_append_to_analysis: str = None,
      summary_html_file_name = None,  # write to these files if supplied
      analysis_html_file_name = None, # write to these files if supplied
      provider: str = None, 
      max_analyze_log_tries = 3,
      log_directory: str = None, # Where to save the history
      job_id: str = None,        # The ID for this specific run (e.g., '23')
      run_agent: bool = False,   # are we going to run the agent
      history_files: list = None,  # text files here
      history_simple_string: str = None, # content from client
      project_advice: str = None,  # User guidance
      original_files: list = None,  # List of ground-truth files
      project_state_json: str = None, # Database
      file_list_as_simple_string: str = None,
      ):

    if provider is None:
      provider = os.getenv("LLM_PROVIDER", "ollama")
    debug_log = []

    # --- Parse file list ---
    file_list = None
    if file_list_as_simple_string:
        from phenix.rest import simple_string_as_text
        try:
            file_list_text = simple_string_as_text(file_list_as_simple_string)
            file_list = [f.strip() for f in file_list_text.split('\n') if f.strip()]
        except Exception as e:
            print(f"Error parsing file list: {e}")
    # -----------------------


    # --- FIX: Handle string "None" from command line/server args ---
    if str(log_directory).lower() == 'none' or log_directory == '':
        log_directory = None
    # ---------------------------------------------------------------


    # What we are going to return (Backwards Compatible Object)
    from libtbx import group_args
    log_info = group_args(group_args_type = 'log summary',
      summary = existing_summary,
      analysis = existing_analysis,
      log_text = None,
      processed_log_dict = None,
      error = None,
      next_move = None,
      history_record = None,
      debug_log = debug_log,
      )

    # Ensure tools are imported to extract files. Not needed when run locally
    try:
        from libtbx.langchain import langchain_tools as lct
    except ImportError:
        lct = None
    except Exception as e:
        lct = None

    def create_history_record_dict(log_info, next_move=None, debug_log = None):
        """Helper to create the dict for JSON serialization"""
        if debug_log is None:
             debug_log = []
        program_name = 'Unknown'
        if log_info.processed_log_dict:
             program_name = log_info.processed_log_dict.get('phenix_program', 'Unknown')

        # --- Fallback Program Name Extraction (For Crashed Jobs) ---
        # If the summary didn't find the name (because of a crash), scan the raw log.
        if (program_name == 'Unknown' or not program_name) and log_info.log_text:
            for line in log_info.log_text.splitlines():
                # Look for the header your client script writes
                if "COMMAND THAT WAS RUN:" in line:
                    # Line format: "COMMAND THAT WAS RUN: phenix.refine ..."
                    parts = line.strip().split()
                    for p in parts:
                        if p.startswith("phenix."):
                            program_name = p
                            break
                    break

        # --- Extract Output Files ---

        state_updates = getattr(log_info, 'state_updates', {})

        # CAPTURE DEBUG
        debug_log.append(f"DEBUG CREATE_HIST: state_updates retrieved = {state_updates}")
        debug_log.append(f"DEBUG CREATE_HIST: Type = {type(state_updates)}")

        out_files = []
        if lct and log_info.summary:
           out_files = lct.extract_output_files(log_info.summary)

        rec = {
            "job_id": job_id,
            "timestamp": time.time(),
            "file_name": file_name,
            "program": program_name,
            "summary": log_info.summary,
            "analysis": log_info.analysis,
            "error": log_info.error,
            "output_files": out_files,
            "state_updates": state_updates,
        }
        debug_log.append(f"DEBUG CREATE_HIST: rec['state_updates'] = {rec.get('state_updates')}")

        if next_move:
            rec['next_move'] = next_move
        return rec

    # Get the text for the log file
    something_to_analyze = True
    if file_name and (not log_as_text):
      if not os.path.isfile(file_name):
        raise ValueError("Sorry, the file %s is missing" %(file_name))
      log_as_text = open(
         file_name, 'r', encoding='utf-8', errors='ignore').read()
      print("Summarizing the log file '%s'..." %(file_name))
    elif file_name and log_as_text:
      print("Summarizing the log file '%s'..." %(file_name))
    elif log_as_text:
      print("Summarizing supplied text as a log file")
      file_name = "Supplied text"
    elif log_info.summary or log_info.analysis:
      print("Summarizing existing info")
      file_name = "Existing info"
    else:
      print("No text to summarize")
      file_name = "No supplied text"
      something_to_analyze = False
    # Decide where to write files
    if not output_file_path:
      import tempfile
      dd = tempfile.TemporaryDirectory()
      output_file_path = dd.name
      # Note: will be deleted when execution ends


    if something_to_analyze and (
       (not log_info.summary) or (not log_info.analysis)):

      # --- 1. DETECT CRASH (Construct Summary) ---
      if log_as_text:
          lines = log_as_text.splitlines()
          for i, line in enumerate(lines):
              if line.strip().startswith("Sorry:"):
                  # Capture the Sorry line plus the next 20 lines (context/suggestions)
                  context = "\n".join(lines[i:i+20])
                  log_info.error = context
                  print(f"NOTE: Found fatal error in log: {line.strip()}")

                  # Set summary so we skip expensive summarization but run Analysis on the error
                  log_info.summary = f"The job failed with the following error and context:\n{context}"
                  break
      # -------------------------------------------
      # --- 2. SETUP (Always run) ---
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

      if db_dir is None:
        if provider == 'google':
            db_dir = "./docs_db_google"
            print(f"Provider is 'google', using default database: {db_dir}")
        elif provider == 'openai':
            db_dir = "./docs_db_openai"
            print(f"Provider is 'openai', using default database: {db_dir}")
        elif provider == 'ollama':
            db_dir = "./docs_db_ollama"
            print(f"Provider is 'ollama', using default database: {db_dir}")

      if not os.path.exists(db_dir):
          log_info.error = (
              f"Database not found for provider '{provider}'. "
              f"Expected at: {db_dir}\n"
              "Please run the database build script for this provider."
          )
          return log_info


      print("Setting up LLMs...")
      expensive_llm, embeddings = lct.get_expensive_llm(
          provider=provider, timeout=timeout)
      cheap_llm = lct.get_cheap_llm(
          provider=provider, timeout=timeout)
         
      # ------ DONE WITH SETUP --------


      # --- 3. SUMMARIZE (Skip if we found a crash) ---
      if lct and (not log_info.summary):
          print("Summarizing log file (using cheap Google model)...")
          result = asyncio.run(lct.get_log_info(
              log_as_text, cheap_llm, embeddings, timeout=timeout, provider='google'))

          if result.error or not result.summary:
            print("Log file summary failed")
            log_info.error = result.error
            return log_info

          log_info = result
          log_info.analysis = None


    elif something_to_analyze:
      print(f"NOTE: Analysis already prepared")
    else:  # nothing to analyze
      print(f"NOTE: No text analyzed")
      log_info.summary = None
      log_info.analysis = None
      # -------------------------------------------
    log_info.log_text = log_as_text

    # Put log summary in an html window
    if display_results and log_info.summary:
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
    if log_info.summary and not log_info.error:
        print("\nAnalyzing summary in context of documentation (using expensive model)...")
        if lct and db_dir and (not log_info.analysis):
          ok = False
          from copy import deepcopy
          log_info_std = deepcopy(log_info)

          # --- RETRY LOGIC ---
          last_error = ""
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
              last_error = result.error or "Unknown analysis error"
              print(f"Analysis failed (Attempt {i+1}/{max_analyze_log_tries}): {last_error}")
              time.sleep(1)

          if (not ok):
            log_info.error = f"Analysis failed after {max_analyze_log_tries} attempts. Last error: {last_error}"
            log_info.analysis = ""
            print("Unable to carry out analysis of log file.")
            return log_info # Return object

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


    # --- Parse Project State ---
    project_state = {}
    if project_state_json:
        try:
             from phenix.rest import simple_string_as_text
             # Decode if simple string, or just load if json
             try:
                 s = simple_string_as_text(project_state_json)
                 project_state = json.loads(s)
             except Exception as e:
                 project_state = json.loads(project_state_json)
        except Exception as e: pass

    # --- Extract Updates ---
    state_updates = {}
    # Ensure we have something to analyze and the agent is enabled
    if (log_info.summary or log_info.error) and run_agent:
        try:
             state_updates = asyncio.run(lct.extract_project_state_updates(
                 log_info.summary, project_state, expensive_llm))
             # CAPTURE DEBUG
             debug_log.append(f"DEBUG EXTRACT: state_updates = {state_updates}")
             debug_log.append(f"DEBUG EXTRACT: Type = {type(state_updates)}, Empty = {len(state_updates) == 0}")

        except Exception as e:
             print(f"Warning: Could not extract state updates: {e}")
             debug_log.append(f"DEBUG EXTRACT: Exception = {e}")

    # --- PACK FOR TRANSPORT (Critical Fix) ---
    if state_updates:
        try:
             from phenix.rest import text_as_simple_string
             updates_json = json.dumps(state_updates)
             # Attach to log_info so server sends it back
             log_info.state_updates_as_simple_string = text_as_simple_string(updates_json)
             log_info.state_updates = state_updates
             # CAPTURE DEBUG
             debug_log.append(f"DEBUG PACK: Successfully packed state_updates")
             debug_log.append(f"DEBUG PACK: Has state_updates attr = {hasattr(log_info, 'state_updates')}")

        except Exception as e:
             print(f"Warning: Failed to encode state_updates: {e}")
             debug_log.append(f"DEBUG PACK: Encoding exception = {e}")
    else:
        debug_log.append(f"DEBUG PACK: state_updates is empty, NOT packing")

    # -----------------------------------------

    # Save History to JSON ---
    if something_to_analyze:
      history_record_dict = create_history_record_dict(log_info, debug_log=debug_log)
    else:
      history_record_dict = {}

    # --- ADD UPDATES TO RECORD ---
    if state_updates:
        history_record_dict['state_updates'] = state_updates
    # -----------------------------

    log_info.history_record = history_record_dict
    # ADD: Attach debug_log to history record for transport
    if debug_log:
        history_record_dict['debug_log'] = debug_log
        log_info.debug_log = debug_log


    # --- UPDATED CONDITION: Save if Summary OR Error exists ---
    if log_directory and (
          (log_info.summary and log_info.analysis) or log_info.error or
          (not something_to_analyze)
       ):
        try:
            if not os.path.exists(log_directory):
                os.makedirs(log_directory)

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
    else:
      json_path = None
    # ---------------------------------

    # --- AGENT LOGIC ---
    if run_agent:
        print("\n--- DATA AVAILABILITY CHECK ---")
        if original_files:
            print(f"Original Files (User Provided): {original_files}")
        else:
            print("Original Files: None provided")

    # Run if we have a summary OR if we caught a fatal error (to fix it) OR
    #   there was no log file at all
    from libtbx.utils import Sorry
    # Run if we have a summary OR if we caught a fatal error OR if we are just starting/resuming (no text)
    if run_agent and (
      log_info.summary or log_info.error or (not something_to_analyze)):
        print("\n--- Running Agent for Next Move ---")

        try:
            # 1. Ensure LLM is loaded (Fix for reused summaries)
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
                        except Exception as e: pass

            # 3. Append Current Run
            if something_to_analyze:
              full_history = past_history + [history_record_dict]
            else:
              full_history = past_history

            # 4. Run Agent
            from libtbx.langchain import langchain_tools as lct
            agent_result = asyncio.run(lct.generate_next_move(
                run_history=full_history,
                llm=expensive_llm,
                embeddings=embeddings,
                db_dir=db_dir,
                cheap_llm=cheap_llm,
                timeout=timeout*3,
                project_advice=project_advice, # Pass User Context
                original_files=original_files,  # Pass File Context
                project_state=project_state,
                file_list=file_list,
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
            # --- FIX: Report crash to client ---
            next_move_dict = {
                "program": "CRASH",
                "explanation": "The Agent crashed while generating a move.",
                "strategy": "Debug",
                "command": "No command generated.",
                "process_log": f"AGENT CRASH TRACEBACK:\n{e}",
                "error": str(e)
            }
            # Attach to object AND dict
            log_info.next_move = next_move_dict
            history_record_dict['next_move'] = next_move_dict
            log_info.history_record = history_record_dict
            # -----------------------------------


    # Debug check at end of run()
    if hasattr(log_info, 'next_move') and log_info.next_move:
        print("DEBUG SERVER: Returning log_info with next_move attached.")
    else:
        print("DEBUG SERVER: log_info has NO next_move.")

    time.sleep(1.0)  # give document loader time to load before deleting files
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

    # Allow raw error text input
    parser.add_argument("--error", type=str, help="Paste the error message string here")

    parser.add_argument("--db_dir", default="./docs_db", help="Database directory")
    parser.add_argument("--log_dir", default=None, help="Directory to save history JSON")
    parser.add_argument("--job_id", default=None, help="Job ID for history")
    parser.add_argument("--provider", default=None, help="LLM Provider (ollama, google, openai)")


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
        provider = args.provider or os.getenv("LLM_PROVIDER", "ollama")
    )
