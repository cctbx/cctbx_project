"""Tool to recursively add docstrings to python files in a directory"""

from __future__ import division

try:
  import google.generativeai as genai
  from google.api_core.exceptions import ResourceExhausted # Add this line
except Exception as e:
  genai = None

import os
import time
import argparse

# --- Configuration ---
CHUNK_SIZE_LIMIT = 150  # Process files larger than this in chunks
CHUNK_SIZE_LIMIT_MIN = 100  # Minimum size to try after first run
CHUNK_SIZE_LIMIT_MAX = 300  # Maximum size to try after first run
MAX_SIZE_LIMITS = 4 # maximum tries with different size limits
MAX_LOCAL_TRIES = 10 #  maximum tries if gemini does not return all sent text

# --- Prompts ---
PROMPT_FIRST_CHUNK = """
        Please add docstrings to the uploaded Python code.
        - Add a module-level docstring at the start of the file.
        - Add a docstring for each method and each class.
        - Include arguments and return values and their types in the docstrings.
        - Make sure that each docstring you add starts with a triple quote
           and ends with a triple quote.
        - VERY IMPORTANT: Do not add, change or remove any existing lines of code.
        - VERY IMPORTANT: Do not change any indentation of any code lines.
        - VERY IMPORTANT: Do not change any spacing or formatting of code lines.
        Then provide the fully updated file content as a raw code block.
"""

PROMPT_SUBSEQUENT_CHUNK = """
        Please add docstrings to the uploaded Python code.
        - Add a docstring for each method and each class.
        - Include arguments and return values and their types in the docstrings.
        - Make sure that each docstring you add starts with a triple quote
           and ends with a triple quote.
        - VERY IMPORTANT: Do not change any indentation of any code lines.
        - VERY IMPORTANT: Do not change any spacing or formatting of code lines.
        - VERY IMPORTANT: Do not add, change or remove any existing lines of code.
        Then provide the fully updated file content as a raw code block.
"""


# --- Global variable for the model ---
model = None

def initialize_api():
    """Initializes the Gemini API connection."""
    global model
    try:
        api_key = os.environ.get("GOOGLE_API_KEY")
        if not api_key:
         raise ValueError("""
The GOOGLE_API_KEY environment variable is not set.

   How to Get Your API Key 🔑

Go to Google AI Studio: Navigate to makersuite.google.com/app/apikey. You may need to sign in with your Google account.

Create a New Key: Click the button labeled "Create API key in new project".

Copy Your Key: A pop-up window will appear displaying your new API key. Click the copy icon next to the key to copy it to your clipboard.

   How to use Your API Key 🔑

Set the GOOGLE_API_KEY environment variable:

csh:   setenv GOOGLE_API_KEY <your-api-key>
bash:  set GOOGLE_API_KEY=<your-api-key>

""")

        if not genai:
         raise ValueError("""
          This tool requires google-generativeai. Install with:

          phenix.python -m pip install -q -U google-generativeai
          """)

        genai.configure(api_key=api_key)
        model = genai.GenerativeModel('gemini-1.5-pro-latest')
        print("✅ Gemini API initialized successfully.")
    except Exception as e:
        print(f"❌ Error initializing Gemini API: {e}")
        exit(1)

def strip_markdown_fences(text):
    """Removes markdown code fences from a string."""
    lines = text.splitlines()
    # Check if the first and last lines are fences and remove them
    if lines and lines[0].strip().startswith("```"):
        lines.pop(0)
    if lines and lines[-1].strip() == "```":
        lines.pop()
    return "\n".join(lines)

def file_needs_processing(file_path):
    """Checks if a file contains at least one class or method definition."""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.strip().startswith(('def ', 'class ')):
                    return True
    except Exception as e:
        print(f"  Could not read {file_path}: {e}")
    return False

def get_code_line_ranges(original_lines, doc_string_line_ranges = None):
  code_line_ranges = []
  code_line_range = None
  for i in range(len(original_lines)):
     if original_lines[i].strip() != '' and (
         not is_in_doc_string_line_ranges(i, doc_string_line_ranges)): # coding
       if code_line_range:
          code_line_range[1] = i
       else:
          code_line_range = [i,i]
          code_line_ranges.append(code_line_range)
     else:
       if code_line_range: # end coding
         code_line_range = None
       else:
         pass
  return code_line_ranges

def get_total_lines(doc_string_line_ranges):
  n = 0
  for i1, i2 in doc_string_line_ranges:
    n += i2 + 1 - i1
  return n

def get_included_lines(doc_string_line_ranges, all_lines):
  included_lines = []
  for i1, i2 in doc_string_line_ranges:
    included_lines += all_lines[i1:i2+1]
  return included_lines


def all_words_in_doc_string_are_present(doc_string_lines, lines):
  if not lines:
    return False
  doc_string_text = " ".join(doc_string_lines)
  lines_text = " ".join(lines)
  for x in ['"""',"'''",'"',"'","#","."]:
    doc_string_text = doc_string_text.replace(x," ")
    lines_text = lines_text.replace(x," ")
  doc_words = doc_string_text.split()
  lines_words = lines_text.split()
  if doc_words == lines_words:
    return True
  if doc_words and not lines_words:
    return False
  if lines_words and not doc_words:
    return False
  if len(doc_words) > len(lines_words):
    return False
  if not doc_words[0] in lines_words:
    return False

  ii = lines_words.index(doc_words[0])

  test_words = lines_words[ii:ii+len(doc_words)]
  if test_words == doc_words:
    return True
  else:
    return False

def get_doc_string_block_lines(i,original_lines, orig_doc_string_line_ranges):
  first_line = None
  last_line = None
  for i1, i2 in orig_doc_string_line_ranges:
    if (first_line is None):
      if i == i1: # found range where i appears
        first_line = i1
        last_line = i2
      else:  # not yet there
        continue
    elif i1 > last_line + 1: # done
      return original_lines[first_line:last_line+1]
    else:  # keep going, contiguous lines
      last_line = i2
  if first_line is not None:
    return original_lines[first_line:last_line+1]
  else:
    return []

def get_code_line_ranges_from_doc_line_ranges(lines, doc_string_ranges):
  """ Everything not in doc_string ranges is a code line range"""
  code_line_ranges = []
  last_doc_string_range = None
  if not doc_string_ranges:  # take all
    if not lines:
      return code_line_ranges
    else:
      code_line_range = [0,len(lines)-1]
      code_line_ranges.append(code_line_range)
      return code_line_ranges
  for i1, i2 in doc_string_ranges:
    if last_doc_string_range is None:
      if i1 > 0:
        code_line_range = [0, i1 -1]
        code_line_ranges.append(code_line_range)
      last_doc_string_range = [i1, i2]
    if i1 > last_doc_string_range[1] + 1:
        code_line_range = [last_doc_string_range[1]+1, i1 -1]
        code_line_ranges.append(code_line_range)
    last_doc_string_range = [i1, i2]
  if last_doc_string_range and (len(lines) > last_doc_string_range[1]+1):
        code_line_range = [last_doc_string_range[1]+1, len(lines)-1]
        code_line_ranges.append(code_line_range)
  return code_line_ranges

def clean_chunk_content(original_lines, gemini_lines,
      return_analysis = False):
    """
    Match coding lines in original with gemini, then insert gemini docstrings
    """

    orig_doc_string_line_ranges = get_line_ranges_with_doc_strings(
       original_lines, get_comment_and_empty_lines_too = True)
    gemini_doc_string_line_ranges = get_line_ranges_with_doc_strings(
       gemini_lines, get_comment_and_empty_lines_too = True)
    orig_code_line_ranges = get_code_line_ranges_from_doc_line_ranges(
         original_lines, orig_doc_string_line_ranges)
    gemini_code_line_ranges = get_code_line_ranges_from_doc_line_ranges(
         gemini_lines, gemini_doc_string_line_ranges)


    orig_code_lines = get_included_lines(orig_code_line_ranges, original_lines)
    gemini_code_lines = get_included_lines(gemini_code_line_ranges, gemini_lines)

    if get_total_lines(orig_code_line_ranges) != get_total_lines(gemini_code_line_ranges):
      print("Gemini code does not match original...")
      orig_code = remove_doc_string_line_ranges(original_lines,
             get_comment_and_empty_lines_too = True)
      gemini_code = remove_doc_string_line_ranges(gemini_lines,
             get_comment_and_empty_lines_too = True)
      return None

    total_gemini_lines = len(gemini_lines)
    total_code_lines_in_common = get_total_lines(orig_code_line_ranges)
    total_code_lines = 0
    total_gemini_doc_lines_included = 0
    total_original_doc_lines_included = 0

    # step through original_lines...insert gemini doc strings
    new_lines = []
    gemini_position = get_first_position(gemini_code_line_ranges)
    prev_gemini_position = 0
    gemini_lines_to_add = []
    in_doc_string_block = False
    for i in range(len(original_lines)):  # go through original code + docstr
      if gemini_position is not None and (
          gemini_position > prev_gemini_position + 1):
        # there are doc strings before next gemini code line
        gemini_lines_to_add = get_lines_to_add(
             gemini_lines,first_line = prev_gemini_position + 1,
             last_line = gemini_position - 1)
        new_lines += gemini_lines_to_add
        total_gemini_doc_lines_included += len(gemini_lines_to_add)
        prev_gemini_position = gemini_position
        in_doc_string_block = False

      if (not is_in_doc_string_line_ranges(i,orig_doc_string_line_ranges)):
        # original code
        new_lines.append(original_lines[i])
        total_code_lines += 1
        prev_gemini_position = gemini_position
        gemini_position = get_next_position(gemini_code_line_ranges,
           gemini_position)
        gemini_lines_to_add = []
        in_doc_string_block = False
      else:  # a doc string line. Insert if not all doc strings in this
             # block are in gemini_lines_to_add
        if (not in_doc_string_block):  # start doc_string_block
          in_doc_string_block = True
          doc_string_block_lines = get_doc_string_block_lines(i,original_lines,
             orig_doc_string_line_ranges)
          if len(doc_string_block_lines) == 1 and (
             not doc_string_block_lines[0].strip()):  # nothing there (blank)
            include_doc_string_block = True
          elif all_words_in_doc_string_are_present(
              doc_string_block_lines, gemini_lines_to_add):
            include_doc_string_block = False
          else:
            include_doc_string_block = True

        if include_doc_string_block:  # include it
          new_lines.append(original_lines[i])
          total_original_doc_lines_included += 1

    if gemini_position is not None and (
           gemini_position > prev_gemini_position + 1):
        # there are doc strings before next gemini code line
        gemini_lines_to_add = get_lines_to_add(
             gemini_lines,first_line = prev_gemini_position + 1,
             last_line = gemini_position - 1)
        new_lines += gemini_lines_to_add
        total_gemini_doc_lines_included += len(gemini_lines_to_add)


    if return_analysis:
      from libtbx import group_args

      return group_args(group_args_type = "doc string analysis",
        total_code_lines = total_code_lines,
        total_gemini_doc_lines_included = total_gemini_doc_lines_included,
        total_original_doc_lines_included = total_original_doc_lines_included,
        total_gemini_lines = total_gemini_lines,
        total_code_lines_in_common = total_code_lines_in_common,
        total_output_lines = len(new_lines))

    else: # usual
      return new_lines

def get_lines_to_add( gemini_lines,first_line = None, last_line = None):
  new_lines = []
  for i in range(first_line, last_line + 1):
    new_lines.append(gemini_lines[i])
  return new_lines

def get_first_position(orig_code_line_ranges):
  if orig_code_line_ranges:
    return orig_code_line_ranges[0][0]
  else:
    return None

def get_next_position(orig_code_line_ranges, orig_position):
  # Find the next position in orig_code_line_ranges

  if orig_position is None:
    return None

  for (i1, i2), (next_i1, next_i2) in zip(orig_code_line_ranges,
       orig_code_line_ranges[1:] + [(None, None)]):
    if orig_position >= i1 and orig_position <= i2: # found it
      if orig_position < i2:  # next line is it
        return orig_position + 1
      elif (next_i1 is not None):
        return next_i1
      else:
        return None
  return None



def get_new_size_limit(size_limits_tried, min_size = None,
     max_size = None):
  import random
  for i in range(100):
    i = random.randint(min_size, max_size)
    if not i in size_limits_tried:
      return i
  return None

def process_file(input_path, output_path, size_limits_tried = None,
      remove_existing = False):
    """
    Manages the chunking and processing for a single file.
    If size_limits_tried already...try a new one
    """


    try:
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
            original_lines = f.read().splitlines()
    except Exception as e:
        print(f"  ❌ Error reading file: {e}")
        return

    if original_lines and remove_existing:  # remove all docstring lines
       original_lines = remove_doc_string_line_ranges(original_lines,
          get_comment_and_empty_lines_too = False)

    if original_lines and original_lines[0].startswith('"""'):
      print("The file %s starts with a docstring" %(input_path)+
           "...skipping.  " +
           "\nUse 'remove_existing' to continue and remove existing docstrings")
      return

    if size_limits_tried is None:
      size_limit = min(len(original_lines),CHUNK_SIZE_LIMIT)
      size_limits_tried = [size_limit]
      first_try = True
    else:
      size_limit = get_new_size_limit(size_limits_tried,
          min_size = min(len(original_lines)//3,CHUNK_SIZE_LIMIT_MIN),
          max_size = min(len(original_lines),CHUNK_SIZE_LIMIT_MAX),
          )
      size_limits_tried.append(size_limit)
      first_try = False
      if size_limit is None or len(size_limits_tried) > MAX_SIZE_LIMITS:
        print("\nGiving up on %s" %(input_path))
        return

    print(f"▶️  Processing: {input_path} with size limit {size_limit}")

    # Convert doc strings to comments
    original_lines = doc_to_comment(original_lines)

    # --- Flexible Chunking Logic ---
    chunks = []
    if len(original_lines) <= size_limit:
        chunks.append(original_lines)
    else:
        print(f"  File is large ({len(original_lines)} lines), splitting into chunks...")

        split_indices = [i for i, line in enumerate(original_lines) if
            line.lstrip().startswith(('def ', 'class '))]

        start_index = 0
        while start_index < len(original_lines):
            # If there are no more logical split points, take the rest of the file as one chunk.
            if not any(i > start_index for i in split_indices):
                chunks.append(original_lines[start_index:])
                break

            # Find the next split point that is GREATER than the chunk size limit.
            # This ensures each chunk is at least size_limit lines long.
            threshold = start_index + size_limit
            next_split = -1
            for split_pos in split_indices:
                if split_pos > threshold:
                    next_split = split_pos
                    break

            # If a split point beyond the threshold is found, use it.
            # Otherwise, it means the remaining code is smaller than the threshold, so take it all.
            final_end_index = next_split if next_split != -1 else len(original_lines)

            chunks.append(original_lines[start_index:final_end_index])
            start_index = final_end_index


    # --- Process each chunk ---
    master_final_lines = []
    temp_filename = "temp_chunk_for_api.txt"

    last_class_line = None
    for i, chunk_lines in enumerate(chunks):
        is_first_chunk = (i == 0)
        print(f"    - Processing chunk {i + 1} of {len(chunks)} (lines: {len(chunk_lines)})...")

        if not chunk_lines:
            continue


        if (is_first_chunk):
          chunk_lines_to_write = chunk_lines
          chunk_lines_at_start = []
        else:  # add a class statement at beginning if nec
          chunk_lines_at_start = add_class_at_start_if_nec(
            chunk_lines, last_class_line = last_class_line)
          chunk_lines_to_write = chunk_lines_at_start + chunk_lines


        with open(temp_filename, "w", encoding='utf-8') as f:
            f.write("\n".join(chunk_lines_to_write))

        uploaded_file = None

        prompt = PROMPT_FIRST_CHUNK if \
            is_first_chunk else PROMPT_SUBSEQUENT_CHUNK

        gemini_lines = None
        resource_tries = 1
        for j in range(MAX_LOCAL_TRIES):
          try:
            uploaded_file = genai.upload_file(path=temp_filename,
                display_name=temp_filename)
            response = model.generate_content([prompt, uploaded_file],
                request_options={'timeout': 300})

            cleaned_text = strip_markdown_fences(response.text)
            gemini_lines = cleaned_text.splitlines()
          except ResourceExhausted as e:
            # This is the 429 error, so we wait and retry
            wait_time = 2 * (2 ** resource_tries)
            resource_tries += 1
            print(f"    Rate limit exceeded. Waiting {wait_time} sec...")
            time.sleep(wait_time)
            gemini_lines = None # Ensure we retry

          except Exception as e:
            print(f"    ❌ Error processing chunk: {e}")
            gemini_lines = None
          if gemini_lines:  # make sure that we have enough

            # remove lines before original chunk lines
            gemini_lines = remove_lines_before_first_chunk_line(
              gemini_lines, chunk_lines)
            # Remove any comment characters inside doc strings
            gemini_lines = remove_comment_inside_doc_strings(gemini_lines)
            # make sure only changes in doc strings
            gemini_lines = check_for_only_changes_in_doc_strings(
              chunk_lines_to_write, gemini_lines)

          if gemini_lines: # got it
            print("      Got chunk %s on try %s with %s lines" %(
              i+1, j+1, len(gemini_lines)))
            break
          else:
            time.sleep(1)

        if gemini_lines:

          gemini_lines = clean_chunk_content(chunk_lines_to_write,
              gemini_lines)

        if gemini_lines:

          # remove lines before original chunk lines
          cleaned_lines = remove_lines_before_first_chunk_line(
              gemini_lines, chunk_lines)

          master_final_lines.extend(cleaned_lines)

          if uploaded_file:
              genai.delete_file(uploaded_file.name)
          if os.path.exists(temp_filename):
              os.remove(temp_filename)
          if len(chunks) > 1:
              time.sleep(1)

        else:
          print("Failed to get chunk %s. Using original lines" %(i+1))
          master_final_lines.extend(chunk_lines_to_write)

    # Save last class line in chunk
    last_class_line = get_last_class_line(chunk_lines)

    # Postprocess to add comment at top of file if not present
    post_master_final_lines = postprocess(master_final_lines)
    if post_master_final_lines:
      master_final_lines = post_master_final_lines
    else:
      print("Failed in post processing")
      if master_final_lines:
        print("Wrote failed text to %s" %(f.name))
      else:
        print("Skipped with no text")
      master_final_lines = None

    # Make sure that the only thing that is changed is the addition of
    # docstrings in master_final_lines
    if master_final_lines:
      master_final_lines = check_for_only_changes_in_doc_strings(
        original_lines, master_final_lines)

    if master_final_lines:
      test_original_lines = check_for_only_changes_in_doc_strings(
        master_final_lines, original_lines)
      if not test_original_lines:
        print("Final text had extra code lines...skipping")
        master_final_lines = None


    if not master_final_lines:  # failed somewhere
      print("\nFailed to process...trying different size limit\n")
      return process_file(input_path, output_path,
          size_limits_tried = size_limits_tried)

    # Summarize changes from initial to final file
    result = clean_chunk_content(original_lines, master_final_lines,
      return_analysis = True)
    if result:
      print("\nSummary for %s:" %(output_path) +
       " Code lines: %s Doc lines (original): %s (AI): %s" %(
         result.total_code_lines,
         result.total_original_doc_lines_included,
         result.total_gemini_doc_lines_included,))
    else:
      print("\nUnable to obtain summary for %s:" %(output_path))


    # Create the destination directory and save the final file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding='utf-8') as f:
        f.write("\n".join(master_final_lines) + "\n")
    print(f"  ✅ Successfully created: {output_path}\n")

def remove_lines_before_first_chunk_line(cleaned_lines, chunk_lines):
  if (not chunk_lines) or (not cleaned_lines):
    return cleaned_lines

  new_lines = []
  found = False
  for line in cleaned_lines:
    if found or (line.strip()[:50] == chunk_lines[0].strip()[:50]):
      found = True
    if found:
      new_lines.append(line)
  return new_lines

def remove_comment_inside_doc_strings(lines):
  # Remove special text:  triple quote (""" or ''')
  special_texts = ["Checks if a line contains a triple quote",
     "Checks if a line contains exactly two triple quotes"]

  new_lines = []
  for line in lines:
    for special_text in special_texts:
      if line.strip().startswith(special_text):
        line = special_text
    new_lines.append(line)
  lines = new_lines
  line_ranges = get_line_ranges_with_doc_strings(lines)
  # Remove any comment characters in the doc string
  for i1, i2 in line_ranges:
    for ii in range(i1, i2+1):
      lines[ii] = lines[ii].replace("#","")
  return lines

def doc_to_comment(lines):
  line_ranges = get_line_ranges_with_doc_strings(lines)
  # Now convert doc strings to comments
  for i1, i2 in line_ranges:
    indentation = find_next_indentation(lines, i2)
    if not indentation:
      indentation = get_indentation(lines[i1])

    for ii in range(i1, i2+1):
      lines[ii] =  " " * indentation + "#" + lines[ii].replace('"""','').replace("'''",'')
  return lines

def add_class_at_start_if_nec(chunk_lines, last_class_line = None):
  if not chunk_lines or not last_class_line:
    return []
  elif chunk_lines[0].startswith("class ") or chunk_lines[0].startswith("def "):
    return []
  else:
    return [last_class_line]

def get_last_class_line(lines):
  last_class_line = """class dummy_class:"""
  for line in lines:
    if line.startswith("class"):
      last_class_line = line
  return last_class_line

def check_for_only_changes_in_doc_strings(original_lines, master_final_lines):
   orig_lines_no_doc = remove_doc_string_line_ranges(original_lines)
   new_lines_no_doc = remove_doc_string_line_ranges(master_final_lines)

   if "\n".join(orig_lines_no_doc).strip() == \
      "\n".join(new_lines_no_doc).strip():  # definitely ok
     return master_final_lines

   elif len(orig_lines_no_doc) != len(new_lines_no_doc): # definitely not ok
     return None

   else: # What about code block indentation
     indent_orig = None
     indent_new = None
     in_method = False
     for orig, new in zip (orig_lines_no_doc, new_lines_no_doc):
       if orig.lstrip() != new.lstrip(): # definitely not ok
         return None
       elif orig.lstrip().startswith("def ") or \
           orig.lstrip().startswith("class"): # restart
         indent_orig = None
         indent_new = None
         in_method = False
       elif (not in_method) and orig.split('#')[0].strip().endswith(':'):
         in_method = True

       elif in_method and (indent_orig is None):  # set indentation
         indent_orig = orig.find(orig.lstrip())
         indent_new = new.find(new.lstrip())
       elif in_method and (indent_orig is not None):
         work_indent_orig = orig.find(orig.lstrip())
         work_indent_new = new.find(new.lstrip())
         if (indent_orig-indent_new) != (work_indent_orig - work_indent_new):
           # indentation changed
           return None
     # all differ by just indentation
     return master_final_lines

def is_in_doc_string_line_ranges(i, doc_string_line_ranges):
   for i1, i2 in doc_string_line_ranges:
     if i >=i1 and i <= i2:
       return True
   return False

def postprocess(master_final_lines):

    # Adjust indentation if necessary
    master_final_lines = fix_indentation(master_final_lines)

    # Add a starter class statement if necessary
    if master_final_lines and master_final_lines[0].startswith("  ") and \
       master_final_lines[0].strip().startswith("def "):
      extra_line = "class dummy_class(dummy1):"
      master_final_lines = [extra_line] + master_final_lines
    else:
      extra_line = None
    from cctbx_website.command_line.comment_to_docstring import \
       convert_comments_to_docstrings
    original_code = "\n".join(master_final_lines)
    try:
      modified_code = convert_comments_to_docstrings(original_code,
         raise_on_errors = True)
    except Exception as e:
      return None


    lines = modified_code.splitlines()
    if extra_line:
      new_lines = []
      for line in lines:
        if not line.strip() == extra_line.strip():
          new_lines.append(line)
      lines = new_lines

    return lines

def fix_indentation(lines):
  line_ranges = get_line_ranges_with_doc_strings(lines)
  # Now fix indentation
  for i1, i2 in line_ranges:
    indentation = find_next_indentation(lines, i2)
    if indentation and indentation != get_indentation(lines[i1]):
      for ii in range(i1, i2+1):
        old_line = lines[ii]
        lines[ii] = set_indentation(lines[ii], indentation)
        #print("SET",get_indentation(old_line),'\n',old_line,'\n',lines[ii])
  return lines



def set_indentation(line, indentation):

  return " " * indentation + line.lstrip()
def find_next_indentation(lines, i):
  # Find indentation of next line after i that is not a comment
  indentation = None
  for ii in range(i+1, len(lines)):
    if lines[ii] and lines[ii].strip() and (
        not lines[ii].strip().startswith("#")):
      return get_indentation(lines[ii])

def get_indentation(line):
  chars = line.strip()
  indentation = line.find(chars)
  return indentation


def remove_doc_string_line_ranges(original_lines,
        get_comment_and_empty_lines_too = True):
   doc_string_line_ranges = get_line_ranges_with_doc_strings(original_lines,
     get_comment_and_empty_lines_too = get_comment_and_empty_lines_too)
   new_lines = []
   if not original_lines:
     return new_lines
   for i in range(len(original_lines)):
     if not is_in_doc_string_line_ranges(i, doc_string_line_ranges):
       new_lines.append(original_lines[i])
   return new_lines

def line_contains_two_triple_quotes(line):
  """Return True if line contains two triple quotes of either kind.
   Ignore anything after hash symbol"""
  line = line.split("#")[0]
  if len(line.strip().split('"""')) == 3:
    return True
  elif len(line.strip().split("'''")) == 3:
    return True
  else:
    return False

def line_contains_triple_quote(line, start = False):
  """Return True if line contains triple quote of either kind.
   Ignore quoted triple string.
   Ignore anything after hash symbol"""
  if start:
    if (line.strip().startswith('"""')) or (line.strip().startswith("'''")):
      return True
    else:
      return False

  else: # usual, but check for quoted triple quote and ignore
    line = line.split("#")[0]
    triple = '"""'
    other_triple = "'''"
    quoted_triple = "'%s'" %(triple)
    quoted_other_triple = '"%s"' %(other_triple)
    if ((line.find(triple)>-1) or
        (line.find(other_triple)> -1)) and (
         line.find(quoted_triple) < 0) and (
         line.find(quoted_other_triple) < 0):
      return True
    else:
      return False

def get_line_ranges_with_doc_strings(lines,
       get_comment_and_empty_lines_too = False):
  # Find all triple-quoted doc strings.  Previous line not ending with \\
  #  and not quoted triple-quoted strings
  in_quote = False
  quote_is_bare = None
  line_range = None
  line_ranges = []
  if not lines:
    return line_ranges
  for i in range(len(lines)):
    line = lines[i]
    if line and line.strip() and \
        (not in_quote) and get_comment_and_empty_lines_too:
      if line.strip()[0] == "#":
        line_ranges.append([i,i])
        continue # go on to next line
    elif (not line) and (not in_quote) and get_comment_and_empty_lines_too:
        line_ranges.append([i,i])
        continue # go on to next line
    line = line.split("#")[0]
    if (not line_contains_triple_quote(line)):
       continue
    # we have a """ or ''':
    if in_quote:
      if quote_is_bare:
        line_range.append(i)
        line_ranges.append(line_range)
      in_quote = False
      line_range = None
      quote_is_bare = None
      if line_contains_two_triple_quotes(line):
        in_quote = True
        quote_is_bare = False

    elif i > 0 and lines[i-1].endswith("\\"): # an assignment
      in_quote = True
      quote_is_bare = False
      if line_contains_two_triple_quotes(line):
        in_quote = False
        quote_is_bare = None
    elif line_contains_triple_quote(line, start = True):
      quote_is_bare = True
      in_quote = True
      line_range = [i]
      if line_contains_two_triple_quotes(line):
        line_range.append(i)
        line_ranges.append(line_range)
        in_quote = False
        quote_is_bare = None
    else:
      in_quote = True
      quote_is_bare = False
      if line_contains_two_triple_quotes(line):
        in_quote = False
        quote_is_bare = None
  return line_ranges


# --- Main execution ---
if __name__ == "__main__":
    import sys
    if 'remove_existing' in sys.argv:
      sys.argv.remove('remove_existing')
      remove_existing = True
      print("Removing existing docstrings")
    else:
      remove_existing = False
    parser = argparse.ArgumentParser(
        description="Recursively find Python files, add docstrings using Gemini, and save to a new directory tree.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_dir", help="The source directory to scan for .py files.")
    parser.add_argument("output_dir", help="The destination directory to save modified files.")
    args = parser.parse_args()
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory not found at '{args.input_dir}'")
        exit(1)

    initialize_api()

    print(f"\nScanning for processable .py files in '{args.input_dir}'...")
    for root, dirs, files in os.walk(args.input_dir):
        dirs[:] = [d for d in dirs if not d.startswith(('.', '__')) and d not in ('venv', '.venv')]
        for filename in files:
            if filename.endswith(".py"):
                input_filepath = os.path.join(root, filename)
                if file_needs_processing(input_filepath):
                    relative_path = os.path.relpath(input_filepath, args.input_dir)
                    output_filepath = os.path.join(args.output_dir, relative_path)
                    if os.path.isfile(output_filepath):
                      print("Skipping %s which already exists" %(
                        output_filepath))
                    else:
                      process_file(input_filepath, output_filepath,
                        remove_existing = remove_existing)
                      time.sleep(1)
                else:
                    print(f"⏩ Skipping (no classes/methods found): {input_filepath}")

    if os.path.exists("temp_chunk_for_api.txt"):
        os.remove("temp_chunk_for_api.txt")
    print("\n🎉 All done!")

