"""Run analyze_log to analyze a log file"""
from __future__ import division, print_function
import os
import sys
import glob
import json

from phenix.program_template import ProgramTemplate

from libtbx import group_args
from libtbx.utils import null_out
from libtbx.utils import Sorry
from phenix.command_line.analyze_log import master_params
from phenix.rest import simple_string_as_text, text_as_simple_string

import libtbx.callbacks
import time

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_ai_dir():
  import phenix.phenix_ai
  from pathlib import Path
  return Path(phenix.phenix_ai.__file__).parent

def have_ai_database(provider):
  ai_dir = get_ai_dir()
  if (not ai_dir):
    return False
  ai_db_dir = os.path.join(ai_dir,"docs_db_%s" %(provider))
  if not os.path.isdir(ai_db_dir):
    return False
  return True

def get_recent_history_files(log_directory, max_history=5):
    """
    Locally scans the log_directory for job_*.json files and returns the
    paths to the N most recent ones.
    """
    if not log_directory or not os.path.isdir(log_directory):
        return []

    search_path = os.path.join(log_directory, "job_*.json")
    files = glob.glob(search_path)

    if not files:
        return []

    # Sort files by modification time (Oldest -> Newest)
    files.sort(key=os.path.getmtime)

    # Keep only the last N
    if len(files) > max_history:
        files = files[-max_history:]

    return files

def save_history_locally(history_record, log_directory, out=sys.stdout):
    """
    Saves the history record returned by the server to the local disk.
    """
    if not log_directory or not history_record:
        return

    try:
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)

        # Determine filename
        job_id = history_record.get('job_id')
        if job_id:
            json_filename = f"job_{job_id}.json"
        else:
            json_filename = f"job_{int(time.time())}.json"
        
        json_path = os.path.join(log_directory, json_filename)

        with open(json_path, 'w') as f:
            json.dump(history_record, f, indent=2)
        
        print(f"Saved run history locally to: {json_path}", file=out)

    except Exception as e:
        print(f"Warning: Failed to save local history JSON: {e}", file=out)

def text_to_append_to_summary():
  return """

**For more help:**

You can open up a [Phenix chat](https://phenix-online.org/chatbot) and ask questions or paste in a part of the summary above to get more details and context.

Have a look at the Analysis of this run as well.
"""

def text_to_append_to_analysis():
  return """

**For more help:**

You can open up a [Phenix chat](https://phenix-online.org/chatbot) and
ask questions or paste in a part of the analysis to get more details
and context.

Have a look at the Summary of this run as well.
"""

def set_api_keys(params):
    original_api_key_info = group_args(
      group_args_type = 'original_api_key_info',
      original_google_key = os.environ.get("GOOGLE_API_KEY",None),
      original_openai_key = os.environ.get("OPENAI_API_KEY",None),
      original_cohere_key = os.environ.get("COHERE_API_KEY",None),
)
    if params.communication.openai_api_key:
       os.environ["OPENAI_API_KEY"] = params.communication.google_api_key
    if params.communication.google_api_key:
       os.environ["GOOGLE_API_KEY"] = params.communication.google_api_key
    if params.communication.cohere_api_key:
       os.environ["COHERE_API_KEY"] = params.communication.cohere_api_key
    return original_api_key_info

def unset_api_keys(original_api_info):
    if original_api_info.original_google_key:
       os.environ["GOOGLE_API_KEY"] = original_api_info.original_google_key
    elif os.environ.get("GOOGLE_API_KEY",None):
       del os.environ["GOOGLE_API_KEY"]

    if original_api_info.original_openai_key:
       os.environ["OPENAI_API_KEY"] = original_api_info.original_openai_key
    elif os.environ.get("OPENAI_API_KEY",None):
       del os.environ["OPENAI_API_KEY"]

    if original_api_info.original_cohere_key:
       os.environ["COHERE_API_KEY"] = original_api_info.original_cohere_key
    elif os.environ.get("COHERE_API_KEY",None):
       del os.environ["COHERE_API_KEY"]

def get_files_from_return_value(return_value):
    summary = None
    summary_file_name = None
    analysis = None
    analysis_file_name = None
    for (fn, t) in return_value:
        if t == "Log summary":
          summary_file_name = fn
          if fn and os.path.isfile(fn):
            summary = open(fn).read()
          else:
            summary_file_name = None
            summary = None
        if t == "Log analysis":
          analysis_file_name = fn
          if analysis_file_name and os.path.isfile(analysis_file_name):
            analysis = open(fn).read()
          else:
            analysis_file_name = None
            analysis = None
    return summary, summary_file_name, analysis, analysis_file_name

def markdown_fn_to_html_fn(fn):
  if not fn:
    return ""
  base, ext = os.path.splitext(fn)
  return "%s_%s.html" %(base, ext.replace(".",""))

def get_program_name_from_log_file_name(log_file):
  if not log_file:
    return None
  base, ext = os.path.splitext(os.path.basename(log_file))
  if ext != ".log":
    return None
  if base.find("_") < -1:
    return None
  program_name = "_".join(base.split("_")[:-1])
  return program_name

def get_info_from_display_text_and_file_list(
   params = None, out = sys.stdout):
  """ return info from display_text_as_simple_string or file_list_as_simple_string
  """
  text = ""
  if not params:
    return text

  if params.analyze_log.display_text_as_simple_string:
     display_text = simple_string_as_text(
        params.analyze_log.display_text_as_simple_string)
     if display_text:
       text += """
    Key results displayed in GUI:
%s
       """  %(display_text)

  if params.analyze_log.file_list_as_simple_string:
     file_list = simple_string_as_text(
        params.analyze_log.file_list_as_simple_string)
     text += "\n  Key output files displayed in GUI: \n"
     for fn in file_list.split():
       text += "%s\n" %(fn)

  return text

def get_results_from_all(params = None,
      return_value = None,
      history_record = None,
      analysis = None,
      summary = None):
  if history_record:
    # Check for object or dict
    if isinstance(history_record, dict):
        analysis = history_record.get('analysis')
        summary = history_record.get('summary')
    else:
        analysis = getattr(history_record, 'analysis', None)
        summary = getattr(history_record, 'summary', None)

  if return_value:
    existing_summary, summary_file_name, \
         existing_analysis, analysis_file_name = \
      get_files_from_return_value(return_value)
  elif analysis and summary:
    existing_summary = summary
    existing_analysis = analysis
  else:
    return_value = []
    existing_summary = simple_string_as_text(
       params.analyze_log.summary_as_simple_string)
    existing_analysis = simple_string_as_text(
       params.analyze_log.analysis_as_simple_string)

  if existing_summary and existing_analysis:

    if params.analyze_log.write_files and \
       params.analyze_log.summary_file_name:
      f = open(params.analyze_log.summary_file_name, 'w')
      print(existing_summary, file = f)
      f.close()
    if params.analyze_log.write_files and \
       params.analyze_log.analysis_file_name:
      f = open(params.analyze_log.analysis_file_name, 'w')
      print(existing_analysis, file = f)
      f.close()

  elif params.analyze_log.summary_file_name and (
      os.path.isfile(params.analyze_log.summary_file_name)) and \
    params.analyze_log.analysis_file_name and (
      os.path.isfile(params.analyze_log.analysis_file_name)) and (
      params.analyze_log.load_existing_analysis):

    existing_summary = open(params.analyze_log.summary_file_name).read()
    existing_analysis = open(params.analyze_log.analysis_file_name).read()

  if (not return_value) and existing_summary and existing_analysis:
    return_value = []
    return_value.append(
      (params.analyze_log.summary_file_name, "Log summary"))
    return_value.append(
      (params.analyze_log.analysis_file_name, "Log analysis"))

  if not history_record:  # construct one
    history_record = {
        "job_id": None,
        "timestamp": time.time(),
        "file_name": None,
        "program": 'Unknown',
        "summary": existing_summary,
        "analysis": existing_analysis,
        "error": None,
    }
  return group_args(group_args_type = 'working result',
     return_value = return_value,
     summary = existing_summary,
     analysis = existing_analysis,
     summary_file_name = params.analyze_log.summary_file_name,
     analysis_file_name = params.analyze_log.analysis_file_name,
     history_record = history_record,
     next_move = getattr(history_record, 'get', lambda x: None)('next_move') if isinstance(history_record, dict) else getattr(history_record, 'next_move', None)
     )

def display_results(params = None, working_results = None):
    if not params.analyze_log.display_results:
       return # nothing to do
    from libtbx.langchain.run_analyze_log import run as run_analyze_log
    # This run call is just for display
    run_analyze_log(
       log_as_text = simple_string_as_text(
         params.analyze_log.log_as_simple_string),
       file_name = params.analyze_log.log_file,
       existing_summary = working_results.summary,
       existing_analysis = working_results.analysis,
       timeout = params.analyze_log.timeout,
       text_to_append_to_summary = text_to_append_to_summary(),
       text_to_append_to_analysis = text_to_append_to_analysis(),
       summary_html_file_name = markdown_fn_to_html_fn(
        params.analyze_log.summary_file_name),
       analysis_html_file_name = markdown_fn_to_html_fn(
        params.analyze_log.analysis_file_name),
       )

def run_job_on_server(params, out = sys.stdout):
    """ Run analyze_log job on server"""

    print("Running analyze_log job on server")
    log_as_simple_string = params.analyze_log.log_as_simple_string
    if not log_as_simple_string:
       log_as_simple_string = text_as_simple_string(
         open(params.analyze_log.log_file).read())

    # Set the server
    params.rest_server.url_type = 'ai'

    args = [
         "log_file=%s" %(params.analyze_log.log_file),
         "log_as_simple_string=%s" %(log_as_simple_string),
         "display_text_as_simple_string=%s" %(
           params.analyze_log.display_text_as_simple_string),
         "file_list_as_simple_string=%s" %(
           params.analyze_log.file_list_as_simple_string),
         "display_results=False",
         "analysis_file_name=None",
         "summary_file_name=None",
         "write_files=False",
         "analyze_log.timeout=%s" %(params.analyze_log.timeout),
         "provider=%s" %(params.communication.provider),
         "log_directory=%s" %(params.analyze_log.log_directory),
         "job_id=%s" %(params.analyze_log.job_id),
         "run_agent=%s" %(params.analyze_log.run_agent)
        ]

    # --- Send User Inputs ---
    if params.analyze_log.project_advice:
         args.append("analyze_log.project_advice=\"%s\"" %(params.analyze_log.project_advice))
    
    if params.analyze_log.original_files:
         files_str = " ".join(params.analyze_log.original_files)
         args.append("analyze_log.original_files=\"%s\"" %(files_str))

    # --- Send history files (Content, not Paths) ---
    if params.analyze_log.history_file:
        import json
        history_content = []
        for hf in params.analyze_log.history_file:
             if os.path.isfile(hf):
                 try:
                     with open(hf, 'r') as f:
                         history_content.append(json.load(f))
                 except: pass

        if history_content:
            # Serialize and encode safely for transport
            json_str = json.dumps(history_content)
            simple_str = text_as_simple_string(json_str)
            args.append("analyze_log.history_simple_string=%s" %(simple_str))
    # -------------------------------------------------

    if params.communication.openai_api_key:
         args.append("communication.openai_api_key=%s" %(
           params.communication.openai_api_key))
    elif params.communication.google_api_key:
         args.append("communication.google_api_key=%s" %(
           params.communication.google_api_key))
    if params.communication.cohere_api_key:
         args.append("communication.cohere_api_key=%s" %(
           params.communication.cohere_api_key))
    program = 'analyze_log'

    from phenix.rest import run_on_server
    server_result = run_on_server(
       args = args,
       program = program,
       rest_server_params = params.rest_server,
       max_wait_time = params.communication.max_wait_time,
       cycle = params.communication.cycle,
       update_wait_time_if_busy =
          params.communication.update_wait_time_if_busy,
       stop_if_internet_not_available =
          params.communication.stop_if_internet_not_available,
       update_wait_time_if_down =
          params.communication.update_wait_time_if_down,
       verbose = params.communication.verbose,
       wait_for_server = params.communication.wait_for_server,
       max_server_tries = params.communication.max_server_tries,
       out = out)

    if server_result.success:
      print("Successful run", file = out)
      
      # --- DECODE NEXT MOVE ---
      next_move = None
      next_move_str = server_result.output_files.get('next_move_as_simple_string', None)
      if next_move_str:
          try:
             import json
             next_move_json = simple_string_as_text(next_move_str)
             next_move = json.loads(next_move_json)
          except Exception as e:
             print(f"Warning: Failed to decode next_move: {e}")
      # ------------------------

      summary = simple_string_as_text(
        server_result.output_files.get('summary_as_simple_string',''))
      analysis = simple_string_as_text(
        server_result.output_files.get('analysis_as_simple_string',''))

      history_record = {
        "job_id": server_result.output_files.get('job_id',''),
        "timestamp": server_result.output_files.get('timestamp',''),
        "file_name": server_result.output_files.get('file_id',''),
        "program": server_result.output_files.get('program',''),
        "summary": summary,
        "analysis": analysis,
        "error": server_result.output_files.get('error',''),
        "next_move": next_move # Add to local record
    }
      
      # --- SAVE LOCALLY (New) ---
      save_history_locally(history_record, params.analyze_log.log_directory, out=out)
      # --------------------------

      working_results = get_results_from_all(params = params,
        history_record = history_record)
      return working_results

    else:
        print("Run failed", file = out)
        print("Message: %s" %(server_result.server_message), file = out)
        raise Sorry("%s" %(server_result.server_message))


def run_job_locally(params, out = sys.stdout):

  """ Run analyze_log job locally

  Analyze log file specified by user based on database in phenix_ai/docs_db

  Args:
    params : Params file

  """
  print("Running analyze_log job locally")
  log_file = params.analyze_log.log_file
  log_as_text = simple_string_as_text(params.analyze_log.log_as_simple_string)
  timeout = params.analyze_log.timeout

  if (not log_file) and (not log_as_text):
    print("Please name a log file to analyze")
    return None

  if not log_as_text:
    if not os.path.isfile(log_file):
      raise Sorry("The file %s is missing" %(log_file))
    log_as_text = open(log_file).read()

  program_name_as_text = get_program_name_from_log_file_name(log_file)
  summary_info_as_text = get_info_from_display_text_and_file_list(params,
     out = out)

  if program_name_as_text:
    log_as_text = "Log file for phenix.%s\n\n%s" %(
        program_name_as_text, log_as_text)
  if summary_info_as_text:
    log_as_text = "%s\n\n%s" %(log_as_text, summary_info_as_text)

  # Get directory with database
  ai_dir = get_ai_dir()
  if (not ai_dir):
    print("Missing phenix_ai/ directory...cannot use database")
    return None
  ai_db_dir = os.path.join(ai_dir,"docs_db_%s" %(params.communication.provider))
  if not os.path.isdir(ai_db_dir):
    raise Sorry("Please run `phenix.install_ai_tools` to install database.")
  
  # Query the database or just put up results if present
  from libtbx.langchain.run_analyze_log import run as run_analyze_log
  
  # --- NEW: Extract history files AND pass new params locally ---
  history_files = params.analyze_log.history_file
  if params.analyze_log.log_directory and not history_files:
      files = get_recent_history_files(
          params.analyze_log.log_directory,
          params.analyze_log.max_history
      )
  
  log_info = run_analyze_log(file_name = log_file,
     log_as_text = log_as_text,
     db_dir = ai_db_dir,
     timeout = timeout,
     display_results = params.analyze_log.display_results,
     text_to_append_to_summary = text_to_append_to_summary(),
     text_to_append_to_analysis = text_to_append_to_analysis(),
     provider = params.communication.provider,
     summary_html_file_name = markdown_fn_to_html_fn(
       params.analyze_log.summary_file_name),
     analysis_html_file_name = markdown_fn_to_html_fn(
       params.analyze_log.analysis_file_name),
     log_directory = params.analyze_log.log_directory,
     job_id = params.analyze_log.job_id,
     run_agent = params.analyze_log.run_agent, 
     history_files = history_files,
     # --- Pass new params ---
     project_advice = params.analyze_log.project_advice,
     original_files = params.analyze_log.original_files
     )

  # --- CRITICAL CHANGE: Do NOT raise Sorry here ---
  # Check if there was a fatal error returned in the object
  err = getattr(log_info, 'error', None)
  if err:
      print(f"NOTE: Log analysis reported an error: {err}. Proceeding to return history record.", file=out)
  # ------------------------------------------------

  # --- UNWRAP FOR LOCAL RUNS ---
  history_record = getattr(log_info, 'history_record', None)
  
  # --- SAVE LOCALLY (New) ---
  save_history_locally(history_record, params.analyze_log.log_directory, out=out)

  working_results = get_results_from_all(params = params,
    history_record = history_record)
  return working_results

# =============================================================================
# CLASS DEFINITION
# =============================================================================

class Program(ProgramTemplate):
  program_name = 'phenix.analyze_log'
  description = '''
  Analyzes a specified log file. This program takes a log file as input
  and processes it using an external analysis function.
  '''
  datatypes = ['phil']

  master_phil_str = master_params

  def validate(self):
    print('Validating inputs for phenix.analyze_log...', file=self.logger)
    if (not self.params.analyze_log.log_file) and (
        not self.params.analyze_log.log_as_simple_string):
      raise Sorry("An input log file must be specified using "+
         "'analyze_log.log_file='.")
    if self.params.analyze_log.log_file and (
        not self.params.analyze_log.log_as_simple_string) and (
        not os.path.isfile(self.params.analyze_log.log_file)):
        raise Sorry(
     f"The log file does not exist: {self.params.analyze_log.log_file}")

    print('Inputs validated successfully.', file=self.logger)

  def print_version_info(self):
    import time
    print ("\n"+60*"*"+"\n"+10*" "+"PHENIX %s" %(self.program_name) +\
           " "+str(time.asctime())+"\n"+60*"*"+"\n",file=self.logger)
    print ("Working directory: ",os.getcwd(),"\n",file=self.logger)
    print ("PHENIX VERSION: ",os.environ.get('PHENIX_VERSION','svn'),"\n",
           file=self.logger)

  def run(self):
    self.print_version_info()
    print('Starting log file analysis with phenix.analyze_log...',
       file=self.logger)
    self.result = run_job_on_server_or_locally(params = self.params,
      out = self.logger)

    # --- Print the Agent Results to console (if available) ---
    if hasattr(self.result, 'history_record') and self.result.history_record:
        hr = self.result.history_record
        nm = None
        if isinstance(hr, dict):
            nm = hr.get('next_move')
        else:
            nm = getattr(hr, 'next_move', None)

        if nm and isinstance(nm, dict):
            # 1. Print Thought Process (Log)
            if nm.get('process_log'):
                 print("\n" + "="*40, file=self.logger)
                 print("AGENT THOUGHT PROCESS:", file=self.logger)
                 print(nm['process_log'], file=self.logger)
                 print("="*40 + "\n", file=self.logger)

            # 2. Print Recommendation & Command (The Final Output)
            print("\n" + "="*40, file=self.logger)
            print("RECOMMENDED ACTION", file=self.logger)
            print("="*40, file=self.logger)
            print(f"**Program:** {nm.get('program', 'Unknown')}", file=self.logger)
            print(f"**Explanation:**\n{nm.get('explanation', '')}\n", file=self.logger)
            print(f"**Strategy Details:**\n{nm.get('strategy', '')}\n", file=self.logger)
            print("-" * 20, file=self.logger)
            print("GENERATED COMMAND:", file=self.logger)
            print(f"{nm.get('command', 'No command generated.')}", file=self.logger)
            print("="*40, file=self.logger)
    # -----------------------------------------------------

    return self.result.return_value

  def get_results(self):
    return self.result

  def get_results_as_JSON(self):
    # Handle result being an object (Standard) or a dict (from server logic sometimes)
    def get_val(obj, key):
        if isinstance(obj, dict): return obj.get(key)
        return getattr(obj, key, None)

    # Use the stored history_record (dict) if available, otherwise use object attrs
    hist = getattr(self.result, 'history_record', None)
    
    output_files = {}
    
    if hist:
        output_files = {
          'summary_as_simple_string': text_as_simple_string(get_val(hist, 'summary')),
          'analysis_as_simple_string': text_as_simple_string(get_val(hist, 'analysis')),
          'job_id': get_val(hist, 'job_id'),
          'timestamp': get_val(hist, 'timestamp'),
          'file_name': get_val(hist, 'file_name'),
          'program': get_val(hist, 'program'),
          'error': get_val(hist, 'error'),
        }
        
        # --- SERIALIZE NEXT_MOVE SAFELY ---
        next_move = get_val(hist, 'next_move')
        if next_move:
            import json
            try:
                next_move_json = json.dumps(next_move)
                output_files['next_move_as_simple_string'] = text_as_simple_string(next_move_json)
            except Exception as e:
                print(f"Error serializing next_move: {e}")
        # ----------------------------------
        
    else:
        # Fallback for legacy behavior
        output_files = {
          'summary_as_simple_string': text_as_simple_string(get_val(self.result, 'summary')),
          'analysis_as_simple_string': text_as_simple_string(get_val(self.result, 'analysis')),
          'error': get_val(self.result, 'error')
        }

    json_results = {
      'output_files': output_files
    }
    return json_results

