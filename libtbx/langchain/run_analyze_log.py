"""
Summarizes a log file (usually Phenix) and suggests next steps in the
context of the (Phenix) documentation.  Uses a two-stage pipeline with
an advanced reranking RAG.
"""
from __future__ import division
from libtbx.langchain import langchain_tools as lct
import sys, os
import asyncio

def run(file_name = None, text = None, output_file_path = None,
      db_dir: str = "./docs_db",
      timeout: int = 60,
      existing_summary: str  = None,
      existing_analysis: str  = None,
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
    if not text:
      if not os.path.isfile(file_name):
        raise ValueError("Sorry, the file %s is missing" %(file_name))
      text = open(file_name).read()
      print("Summarizing the log file '%s'..." %(file_name))
    else:
      print("Summarizing supplied text as a log file")
      file_name = "Supplied text"

    # Decide where to write files
    if not output_file_path:
      output_file_path = '/var/tmp'

    if (not log_info.summary) or (not log_info.analysis):
      # need to get the info
      if not os.getenv("GOOGLE_API_KEY"):
         log_info.error = "GOOGLE_API_KEY environment variable not set."
         return log_info
      if not os.getenv("COHERE_API_KEY"):
         log_info.error = "COHERE_API_KEY environment variable not set."
         return log_info


      # Set up the LLM (same one for both analyzing log file and processing)
      print("Setting up LLM")
      llm = lct.get_llm()
      embeddings = lct.get_embeddings()

      # Summarize the log file
      print("Summarizing log file")
      result = asyncio.run(lct.get_log_info(text, llm, embeddings,
        timeout = timeout))
      if result.error: # failed
        print("Log file summary failed")
        log_info.error = result.error
        return log_info
      log_info = result
      log_info.analysis = None

    log_info.log_text = text

    # Put log summary in an html window
    fn = os.path.join(output_file_path,'log_summary.html')
    lct.save_as_html(log_info.summary, file_name = fn,
       title = 'Summary of %s' %(file_name))
    print("Loading log summary at %s" %(fn))
    try:
      from phenix.command_line.doc import load_url
      load_url(fn)
    except Exception as e:
      # phenix is not available or no viewer.  Just skip
      print("Unable to load viewer...see text in the file: '%s'" %(fn))

    # Analyze the log summary in the context of the docs
    print("\nAnalyzing summary in context of documentation...")
    if (db_dir and (not log_info.analysis)):
      result = asyncio.run(
        lct.analyze_log_summary(log_info, llm, embeddings,
        db_dir = db_dir, timeout = timeout))
      if (result.error): # failed
        log_info.error = result.error
        print("Unable to carry out analysis of log file")
        return log_info
      else:
        log_info.analysis = result.analysis

    # Put it in an html window
    fn = os.path.join(output_file_path,'analysis.html')
    lct.save_as_html(log_info.analysis, file_name = fn,
       title = 'Analysis of %s' %(file_name))
    print("Loading analysis at %s" %(fn))
    try:
      from phenix.command_line.doc import load_url
      load_url(fn)
    except Exception as e:
      # phenix is not available or no viewer.  Just skip
      print("Unable to load viewer...see text in the file: '%s'" %(fn))

    return log_info

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: phenix.python run_analyze_log.py <path_to_log_file.txt>")
        sys.exit(1)
    fn = sys.argv[1]
    answer = run(file_name = fn)
