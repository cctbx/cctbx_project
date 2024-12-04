from __future__ import absolute_import, division, print_function
import os, sys


from libtbx.utils import display_context


''' Tool to find instances in code where PDB formatting is assumed and
    should be made general to include PDB or mmcif formatting
'''

# =============================================================================
def find_pdb_mmcif_problems(path, n_context = 7, overall_exclude = None,
    mark_files = None, unmark_files = None):

  # Run on all paths supplied in file if a .list file is supplied
  if os.path.isfile(path) and path.endswith(".list"):
    print("Working on all paths listed in '%s'" %(path))
    all_problems = []
    total_files = 0
    for full_fn in open(path).readlines():
      full_fn = full_fn.strip()
      print("Working on '%s'" %(full_fn))
      if os.path.isdir(full_fn) or full_fn.endswith('.py'):
        problems, number_of_files_examined = find_pdb_mmcif_problems(full_fn,
           n_context = n_context, overall_exclude = overall_exclude,
            mark_files = mark_files, unmark_files = unmark_files)
        all_problems += problems
        total_files += number_of_files_examined
    return all_problems, total_files

  # Run recursively if a directory is supplied
  if os.path.isdir(path):
    print("Working recursively on the path '%s'" %(path))
    all_problems = []
    total_files = 0
    for fn in os.listdir(path):
      if fn.startswith("."): continue
      if fn.startswith("tst_"): continue
      full_fn = os.path.join(path, fn)
      if os.path.isdir(full_fn) or full_fn.endswith('.py'):
        problems, number_of_files_examined = find_pdb_mmcif_problems(full_fn,
           n_context = n_context, overall_exclude = overall_exclude,
           mark_files = mark_files, unmark_files = unmark_files)
        all_problems += problems
        total_files += number_of_files_examined
    return all_problems, total_files

  # Normal run here on one file in path

  if not os.path.isfile(path):
    print("The path %s is not a file" %path)
    return [],0
  if (not path.endswith('.py')):
    print("The file %s is not a python file" %path)
    return [], 0
  if (os.path.split(path)[-1].startswith('tst_')):
    print("The file %s is a test file" %path)
    return [], 0

  # Run on one file
  all_problems = []
  text = open(path).read()
  for problem in (
     write_model_file_without_assignment,
     pdb_write_statements,
     pdb_format_interpretation,
     raw_records,
     pdb_file_name,
   ):
    all_problems += problem(text, file_name = path, n_context = n_context,
      overall_exclude = overall_exclude)

  if mark_files: # mark all problems in the files themselves
    mark_problems_in_file(text,path, all_problems, n_context = n_context)
  if unmark_files: # unmark all problems in the files
    unmark_problems_in_file(text,path)

  return all_problems, 1

def unmark_problems_in_file(text,path):
  if text.find(" # XXX CHECK PDB:") < 0:
     return # nothing to do

  lines = text.splitlines()
  for i in range(len(lines)):
    line = lines[i]
    index = line.find(" # XXX CHECK PDB:")
    if index > -1:
      lines[i] = line[:index]
  f = open(path,'w')
  print("\n".join(lines), file = f)
  f.close()
  print("Removed CHECK PDB text from %s" %(path))


def mark_problems_in_file(text, path, all_problems,
     n_context = 5):
  if not all_problems:
     return  # nothing to do
  lines = text.splitlines()
  for p in all_problems:
    sw = [p.search_word, p.required_word]
    if None in sw: sw.remove(None)
    line = lines[p.line_number-1].rstrip()
    next_ending = None
    for w in sw:
      if line.find(w):
        next_ending = w
        break
    if (not next_ending):
       next_ending = " ".join(sw)
    new_text = " # XXX CHECK PDB: %s" %(next_ending)
    for i in range(p.line_number-1, p.line_number + n_context - 1):
      line = lines[i]
      if (not line.endswith("\\")) and (not line.find("XXX CHECK PDB") > -1):
        line = line.rstrip()
        blanks = max(0, 80 - len(line) - len(new_text))
        if blanks > 0:
          new_text = blanks *" " + new_text
        line +=  new_text
        lines[i] = line
        break
  f = open(path,'w')
  print("\n".join(lines), file = f)
  f.close()
  print("Marked problems with CHECK PDB text in %s" %(path))

def run(args, n_context = 7, overall_exclude = None):
  if len(args) < 1:
    print("phenix.python find_pdb_mmcif_problems.py <file_name or directory> <n_context> <mark_files> <unmark_files>")
    return

  mark_files = False
  unmark_files = False
  if 'mark_files' in args:
    args.remove('mark_files')
    mark_files = True
  if 'unmark_files' in args:
    args.remove('unmark_files')
    unmark_files = True

  path = args[0]
  if not os.path.exists(path):
    print("The path '%s' does not exist?" %path)
    return []

  if not overall_exclude:
    overall_exclude = ['get_refine_file_stem',
    'DEBUG',        # debugging code
    'PDB OK',       # marked as ok
    'PDB REQUIRED', # marked as ok
    'must be PDB',  # marked as ok
    'MUST BE PDB',  # marked as ok
    'XXX OK',  # marked as ok
    'XXX ok',  # marked as ok
    'set_model_ext_and_target_output_format',  # function to set the extension
    'forward_compatible', # using pdb deliberately
    'transfer_ext', # function to set the extension
    'get_cif_or_pdb_file_if_present', # function to set the extension
    '.pdb_or_mmcif_string_info', # function to set the extension
    'mmcif', # the author has considered cif
    '.cif',  # the author has considered cif
    'write_model(',  # the write_model function is cif-aware
    'write_output_file(',  # the write_output_file function is cif-aware
    '.type = ',  # This is in parameters
    ]

  if len(args) > 1:
    n_context = int(args[1])
    print("Set n_context to ",n_context)
  all_problems, total_files  = find_pdb_mmcif_problems(
     path, n_context = n_context, overall_exclude = overall_exclude,
     mark_files = mark_files, unmark_files = unmark_files)

  # Summarize results

  all_file_names = get_all_unique(all_problems, key = 'file_name')
  all_categories = get_all_unique(all_problems, key = 'category')
  all_search_words = get_all_unique(all_problems, key = 'search_word')

  print("All unique categories: ",all_categories)
  print("All unique file_names: ",all_file_names)
  print("All unique search words: ",all_search_words)

  print("\nAll problems:")
  for file_name in all_file_names:
    print("\n"+70*"=")
    print("FILE NAME: %s" %(file_name))
    print(70*"=")
    problem_list = select_problems(all_problems, key = 'file_name',
      value = file_name, n_context = n_context)
    for p in problem_list:
      display_problem(p)

  print("\nSummary of problems by file:")
  for file_name in all_file_names:
    problem_list = select_problems(all_problems, key = 'file_name',
      value = file_name, n_context = n_context)
    total_problems = 0
    from six.moves import StringIO
    f = StringIO()
    for sw in all_search_words:
       search_words = [sw]
       p_list = select_problems(problem_list, key = 'search_word',
         value = sw, n_context = 0)
       if p_list:
         for p in p_list:
           if p.required_word and (not p.required_word in search_words):
              search_words.append(p.required_word)
         print(" TEXT FOUND:   %s  NUMBER OF TIMES: %s" %(
           search_words, len(p_list)), file = f)
         total_problems += len(p_list)
    print("\n%s: %s problems" %(file_name,total_problems))
    print(f.getvalue())
  total_files_with_problems = len(all_file_names)
  print("\nTotal files with problems: %s (of %s)" %(
    total_files_with_problems, total_files))

  print("\nSummary of problems by code:\n")
  total_problems = 0
  for sw in all_search_words:
     search_words = [sw]
     p_list = select_problems(all_problems, key = 'search_word',
       value = sw, n_context = 0)
     if p_list:
       for p in p_list:
         if p.required_word and (not p.required_word in search_words):
            search_words.append(p.required_word)
       print(" TEXT FOUND:   %s  NUMBER OF TIMES: %s" %(
         search_words, len(p_list)))
       total_problems += len(p_list)

  print("\nTotal problems found in %s: %s" %(path,total_problems))


def display_problem(p, out = sys.stdout):
  print("%s  Line: %s  :: %s\n" %(
     p.file_name,p.line_number,p.category)+ 70*"-", file = out)
  sw = [p.search_word, p.required_word]
  if None in sw: sw.remove(None)
  next_ending = None
  for line in p.text_block.splitlines():
    line = line.rstrip()
    if line.startswith("  **"):
      for w in sw:
        if line.find(w):
          next_ending = w
          break
      if (not next_ending):
          next_ending = " ".join(sw)
    if next_ending and (not line.endswith("\\")):
      line = line.rstrip()
      new_text = " # XXX CHECK PDB: %s" %(next_ending)
      blanks = max(0, 80 - len(line) - len(new_text))
      if blanks > 0:
        new_text = blanks *" " + new_text
      line +=  new_text
      next_ending = None
    print(line, file = out)
  print(70*"-", file = out)

def select_problems(all_problems, key = 'file_name', value = None,
    n_context = 7):
  new_problems = []
  for p in all_problems:
    if getattr(p,key,None) == value:
      new_problems.append(p)
  new_problems = sorted(new_problems, key = lambda p: p.line_number)
  sieved_problems = []
  last_line_number = None
  for p in new_problems:
    if (last_line_number is None) or (
      p.line_number >= last_line_number + n_context - 1):
      sieved_problems.append(p)
      last_line_number = p.line_number
  return sieved_problems

def get_all_unique(all_problems, key = 'file_name'):
  all_unique = []
  for p in all_problems:
    x = getattr(p,key,None)
    if (x is not None) and not x in all_unique:
      all_unique.append(x)
  return all_unique

def write_model_file_without_assignment(text, file_name = None, n_context = 7,
    overall_exclude = None):
  all_problems = []
  for search_word in [ "dm.write_model_file(",
     "data_manager.write_model_file(", ]:
    all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context, search_word = search_word,quiet = True,
       excluded_words = [
         '=dm.write_model_file',
         '= dm.write_model_file',
         '=data_manager.write_model_file',
         '= data_manager.write_model_file',
         '=self.dm.write_model_file',
         '= self.dm.write_model_file',
         '=self.data_manager.write_model_file',
         '= self.data_manager.write_model_file',
            ] + overall_exclude,
       category = 'write_model_file_without_assignment')
  return all_problems

def pdb_write_statements(text, file_name = None, n_context = 7,
    overall_exclude = None):
  all_problems = []
  for search_word in [".model_as_pdb(", ".as_pdb_string(", # PDB OK
      ".write_pdb_file("]:
    all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context, search_word = search_word,quiet = True,
       excluded_words = overall_exclude,
       category = 'pdb_write_statements')
  return all_problems

def pdb_format_interpretation(text, file_name = None, n_context = 7,
      overall_exclude = None):
  all_problems = []
  for search_word in ['.open(']:
    for required_word in ['.pdb"',".pdb'",'pdb_file','model_file']: # PDB OK
      all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context,
       search_word = search_word,
       required_word=required_word,
       excluded_words= overall_exclude,
       quiet = True,
       category = 'pdb_format_interpretation')
  all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context,
       search_word = ".splitlines(",
       required_word="pdb",
       excluded_words=["cmds.","iotbx.pdb",
          'pdb_interpretation.process',
          'traceback','phil_string',] + overall_exclude,
       quiet = True,
       category = 'pdb_format_interpretation')
  for search_word in ['HETATM','"ATOM',"'ATOM",'"TER',"'TER",
         '"BREAK',"'BREAK",' CA ',' N ']: # PDB OK
    for required_word in ['startswith','find(','re.search(']:
      all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context,
       search_word = search_word,
       required_word=required_word,
       excluded_words = ['group_PDB'] + overall_exclude,
       quiet = True,
       category = 'pdb_format_interpretation')
  return all_problems

def raw_records(text, file_name = None, n_context = 7, overall_exclude = None):
  all_problems = []
  for search_word in ['raw_records=','raw_records =']:
    for required_word in ['pdb']:
      all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context,
       search_word = search_word,
       required_word=required_word,
       excluded_words=['pdb.input','iotbx.pdb','input.pdb',
         '.pdb_or_mmcif_string_info','get_cif_or_pdb_file_if_present'
          'pdb_interpretation.process'] +
         overall_exclude,
       quiet = True,
       category = 'raw_records')
  return all_problems

def pdb_file_name(text, file_name = None,
      n_context = 7, overall_exclude = None):
  all_problems = []
  for search_word in ['file_name',' fn=',' fn = ','filename']:
    for required_word in ['.pdb"',".pdb'"]:
      all_problems += display_context(file_name = file_name, text = text,
       n_context = n_context,
       search_word = search_word,
       required_word=required_word,
       excluded_words = ['.pdb_or_mmcif_string_info',
         'get_cif_or_pdb_file_if_present','write_model',
          'pdb.input','iotbx.pdb','.cif','input.pdb',
          'map_file_name','restraint_file',] + overall_exclude,
       quiet = True,
       category = 'pdb_file_name')
  return all_problems


# =============================================================================
if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

# =============================================================================
# end
