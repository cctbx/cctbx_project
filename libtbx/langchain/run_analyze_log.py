"""
Summarizes a log file (usually Phenix) and suggests next steps in the
context of the (Phenix) documentation.  Uses a two-stage pipeline with
an advanced reranking RAG.
"""
from __future__ import division
import sys, os
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
      provider: str = 'openai',
      max_analyze_log_tries = 3,
      ):

    # What we are going to return
    from libtbx import group_args
    log_info = group_args(group_args_type = 'log summary',
      summary = existing_summary,
      analysis = existing_analysis,
      log_text = None,
      processed_log_dict = None,
      error = None)

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

      if not db_dir:
        raise ValueError("Sorry, the database is missing")

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
        #    This model is fast, cheap, and built for this kind of work.
        #    We let it use its default: gemini-2.5-flash-lite
        cheap_llm, _ = lct.get_llm_and_embeddings(
            provider='google',
            timeout=timeout
        )
        # Note: We discard Google's embeddings (_) since we only need its LLM.
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
        return log_info
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

# In run_analyze_log.py, in the "run" function

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
        return log_info
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
    return log_info

# In run_analyze_log.py

def save_as_html(markdown_string: str,
     title: str = "Summary", file_name: str = None):
    """
    Converts a Markdown string to an HTML file.
    If markdown-it-py is not installed, it falls back to a simple
    pre-formatted text block.
    """
    html_content = ""
    try:
        # Try to import the library
        from markdown_it import MarkdownIt

        # If successful, render the full HTML
        md = MarkdownIt("gfm-like")
        html_content = md.render(markdown_string)
        print("Using 'markdown-it-py' for rich HTML rendering.")

    except ImportError:
        # If the import fails, use the fallback
        print("Warning: 'markdown-it-py' not found. Falling back to plain text rendering.")
        # Wrap the raw text in <pre> tags to preserve formatting
        html_content = f"<pre>{markdown_string}</pre>"

    # The rest of the function remains the same, wrapping the content
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
    import sys
    if len(sys.argv) < 2:
        print("Usage: phenix.python run_analyze_log.py <path_to_log_file.txt>")
        sys.exit(1)
    fn = sys.argv[1]
    answer = run(file_name = fn)
