"""
Utility functions for log analysis and LLM operations.
"""
from __future__ import division
import os
import sys
import time

# =============================================================================
# STRING/INPUT PARSING
# =============================================================================

def parse_simple_string_to_list(simple_string, debug_log=None):
  """
  Parse a simple_string (encoded text) into a list of strings.
  Returns empty list on failure.
  """
  if not simple_string:
    return []
  if debug_log is None:
    debug_log = []
  try:
    from phenix.rest import simple_string_as_text
    text = simple_string_as_text(simple_string)
    return [f.strip() for f in text.split('\n') if f.strip()]
  except Exception as e:
    debug_log.append(f"Error parsing simple string to list: {e}")
    return []


def normalize_none_string(value):
  """
  Handle string "None" from command line/server args.
  Returns None if value is 'none', 'None', or empty string.
  """
  if value is None:
    return None
  if str(value).lower() == 'none' or value == '':
    return None
  return value


def load_log_text(file_name, debug_log=None, out=sys.stdout):
  """
  Load log text from a file.
  Returns (log_text, file_name, something_to_analyze) tuple.
  """
  if debug_log is None:
    debug_log = []

  if file_name:
    if not os.path.isfile(file_name):
      raise ValueError(f"Sorry, the file {file_name} is missing")
    log_text = open(file_name, 'r', encoding='utf-8', errors='ignore').read()
    print(f"Summarizing the log file '{file_name}'...", file=out)
    debug_log.append(f"Summarizing the log file '{file_name}'...")
    return log_text, file_name, True

  return None, None, False


def determine_input_source(file_name, log_as_text, existing_summary, existing_analysis,
                           debug_log=None):
  """
  Determine the input source and whether there's something to analyze.
  Returns (log_text, file_name, something_to_analyze) tuple.
  """
  if debug_log is None:
    debug_log = []

  if file_name and not log_as_text:
    if not os.path.isfile(file_name):
      raise ValueError(f"Sorry, the file {file_name} is missing")
    log_text = open(file_name, 'r', encoding='utf-8', errors='ignore').read()
    debug_log.append(f"Summarizing the log file '{file_name}'...")
    return log_text, file_name, True

  if file_name and log_as_text:
    debug_log.append(f"Summarizing the log file '{file_name}'...")
    return log_as_text, file_name, True

  if log_as_text:
    debug_log.append("Summarizing supplied text as a log file")
    return log_as_text, "Supplied text", True

  if existing_summary or existing_analysis:
    debug_log.append("Summarizing existing info")
    return None, "Existing info", True

  debug_log.append("No text to summarize")
  return None, "No supplied text", False


def get_temp_output_path():
  """Create and return a temporary directory path for output files."""
  import tempfile
  dd = tempfile.TemporaryDirectory()
  return dd.name, dd  # Return both path and object to prevent cleanup

# =============================================================================
# API KEY VALIDATION
# =============================================================================

def validate_api_keys(provider, debug_log=None):
  """
  Validate that required API keys are set for the provider.
  Returns error message string if validation fails, None if OK.
  """
  if debug_log is None:
    debug_log = []

  if provider == 'google':
    if not os.getenv("GOOGLE_API_KEY"):
      return "GOOGLE_API_KEY environment variable not set."
  elif provider == 'openai':
    if not os.getenv("OPENAI_API_KEY"):
      return "OPENAI_API_KEY environment variable not set."

  return None

# =============================================================================
# DATABASE DIRECTORY
# =============================================================================

def get_db_dir_for_provider(provider, db_dir=None, debug_log=None):
  """
  Get the database directory for the given provider.
  Uses provided db_dir if set, otherwise returns default for provider.
  Returns (db_dir, error_message) tuple.
  """
  if debug_log is None:
    debug_log = []

  if db_dir is not None:
    if not os.path.exists(db_dir):
      error = (
        f"Database not found for provider '{provider}'. "
        f"Expected at: {db_dir}\n"
        "Please run the database build script for this provider."
      )
      return None, error
    return db_dir, None

  # Default directories based on provider
  defaults = {
    'google': "./docs_db_google",
    'openai': "./docs_db_openai",
    'ollama': "./docs_db_ollama",
  }

  db_dir = defaults.get(provider, "./docs_db_ollama")
  debug_log.append(f"Provider is '{provider}', using default database: {db_dir}")

  if not os.path.exists(db_dir):
    error = (
      f"Database not found for provider '{provider}'. "
      f"Expected at: {db_dir}\n"
      "Please run the database build script for this provider."
    )
    return None, error

  return db_dir, None

# =============================================================================
# LLM SETUP
# =============================================================================

def setup_llms(provider, timeout, debug_log=None):
  """
  Setup LLMs for the given provider.
  Returns (expensive_llm, cheap_llm, planning_llm, embeddings) tuple.
  Returns (None, None, None, None) on failure.
  """
  if debug_log is None:
    debug_log = []

  try:
    from libtbx.langchain.core.llm import get_expensive_llm
    from libtbx.langchain.core.llm import get_cheap_llm


    debug_log.append("Setting up LLMs...")

    expensive_llm, embeddings = get_expensive_llm(
      provider=provider, timeout=timeout, json_mode=False
    )

    # Only create separate planning LLM for Ollama (needs JSON mode)
    if provider == 'ollama':
      planning_llm, _ = get_expensive_llm(
        provider=provider, timeout=timeout, json_mode=True
      )
    else:
      planning_llm = expensive_llm

    cheap_llm = get_cheap_llm(provider=provider, timeout=timeout)

    return expensive_llm, cheap_llm, planning_llm, embeddings

  except Exception as e:
    debug_log.append(f"Error setting up LLMs: {e}")
    return None, None, None, None

# =============================================================================
# SUMMARY VALIDATION
# =============================================================================

HALLUCINATION_MARKERS = [
  "/path/to/file",
  "/nonexistent/",
  "file1.mtz",
  "file2.pdb",
  "file1.pdb",
  "file2.mtz",
]


def validate_summary_for_hallucinations(summary_text, debug_log=None):
  """
  Check if summary contains hallucination markers.
  Returns (is_valid, cleaned_summary) tuple.
  If hallucinations detected, returns replacement text.
  """
  if debug_log is None:
    debug_log = []

  if not summary_text:
    return True, summary_text

  has_hallucination = any(marker in summary_text for marker in HALLUCINATION_MARKERS)

  if has_hallucination:
    debug_log.append("DEBUG SUMMARY: HALLUCINATION DETECTED - replacing summary")
    return False, "Log analysis produced unreliable results. Please re-run the analysis."

  return True, summary_text

# =============================================================================
# HTML OUTPUT
# =============================================================================

def save_as_html(markdown_string, title="Summary", file_name=None, out=sys.stdout):
  """
  Converts a Markdown string to an HTML file.
  Returns the HTML content string.
  """
  html_content = ""
  try:
    from markdown_it import MarkdownIt
    md = MarkdownIt("gfm-like")
    html_content = md.render(markdown_string)
  except ImportError:
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
    print(f"\nSaved formatted output to: {file_name}", file=out)

  return html_with_style

# =============================================================================
# LOG INFO OBJECT
# =============================================================================

def create_log_info_object(existing_summary=None, existing_analysis=None, debug_log=None):
  """
  Create the standard log_info return object.
  """
  if debug_log is None:
    debug_log = []

  from libtbx import group_args
  return group_args(
    group_args_type='log summary',
    summary=existing_summary,
    analysis=existing_analysis,
    log_text=None,
    processed_log_dict=None,
    error=None,
    next_move=None,
    history_record=None,
    debug_log=debug_log,
  )

# =============================================================================
# HISTORY FILE OPERATIONS
# =============================================================================

def save_history_to_json(history_record, log_directory, job_id=None, debug_log=None):
  """
  Save history record to JSON file.
  Returns the path to saved file, or None on failure.
  """
  if debug_log is None:
    debug_log = []

  if not log_directory:
    return None

  try:
    if not os.path.exists(log_directory):
      os.makedirs(log_directory)

    if job_id:
      json_filename = f"job_{job_id}.json"
    else:
      json_filename = f"job_{int(time.time())}.json"

    json_path = os.path.join(log_directory, json_filename)

    import json
    with open(json_path, 'w') as f:
      json.dump(history_record, f, indent=2)

    debug_log.append(f"Saved run history to: {json_path}")
    return json_path

  except Exception as e:
    debug_log.append(f"Warning: Failed to save history JSON: {e}")
    return None

# =============================================================================
# PROGRAM NAME DETECTION
# =============================================================================

def detect_program_from_log(log_text, debug_log=None):
  """
  Try to detect program name from log text.
  Returns program name string or None.
  """
  if debug_log is None:
    debug_log = []

  if not log_text:
    return None

  try:
    from libtbx.langchain.analysis.summarizer import detect_program_from_text
    program_name = detect_program_from_text(log_text)
    if program_name:
      debug_log.append(f"[RUN] Auto-detected program from log: {program_name}")
      return program_name
    else:
      debug_log.append("[RUN] Could not detect program from log")
      return None
  except ImportError:
    debug_log.append("[RUN] Could not import detect_program_from_text")
    return None
