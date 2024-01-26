from __future__ import absolute_import, division, print_function
import os, sys
from libtbx import group_args

from libtbx.utils import display_context


''' Tool to find instances in code where PDB formatting is assumed and
    should be made general to include PDB or mmcif formatting
'''

# =============================================================================
def find_pdb_mmcif_problems(path, n_context = 5):

  # Run recursively if a directory is supplied
  if os.path.isdir(path):
    print("Working recursively on the path '%s'" %(path))
    all_problems = []
    for fn in os.listdir(path):
      if fn.startswith("."): continue
      full_fn = os.path.join(path, fn)
      if os.path.isdir(full_fn) or full_fn.endswith('.py'):
        all_problems += find_pdb_mmcif_problems(full_fn,
           n_context = n_context)
    return all_problems

  if not os.path.isfile(path):
    print("The path %s is not a file" %path)
    return []
  if (not path.endswith('.py')):
    print("The file %s is not a python file" %path)
    return []

  # Run on one file
  all_problems = []
  text = open(path).read()
  for problem in (
     pdb_write_statements,
     pdb_format_interpretation,
     raw_records,
     pdb_file_name,
   ):
    all_problems += problem(text, header = path, n_context = n_context)

  return all_problems

def run(args, n_context = 5):
  if len(args) < 1:
    print("phenix.python find_pdb_mmcif_problems.py <file_name or directory> <n_context>")
    return
  path = args[0]
  if not os.path.exists(path):
    print("The path '%s' does not exist?" %path)
    return []

  if len(args) > 1:
    n_context = int(args[1])
  all_problems = find_pdb_mmcif_problems(path, n_context = n_context)

  # Summarize results

  all_headers = get_all_unique(all_problems, key = 'header')
  all_categories = get_all_unique(all_problems, key = 'category')
  all_search_words = get_all_unique(all_problems, key = 'search_word')

  print("All unique categories: ",all_categories)
  print("All unique headers: ",all_headers)
  print("All unique search words: ",all_search_words)

  print("\nAll problems:")
  for header in all_headers:
    print("\n"+70*"=")
    print("FILE NAME: %s" %(header)) 
    print(70*"=")
    problem_list = select_problems(all_problems, key = 'header',
      value = header, n_context = n_context)
    for p in problem_list:
      display_problem(p)

  print("\nSummary of problems by file:")
  for header in all_headers:
    problem_list = select_problems(all_problems, key = 'header',
      value = header, n_context = n_context)
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
    print("\n%s: %s problems" %(header,total_problems))
    print(f.getvalue())

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
  print(70*"-", file = out)
  sw = [p.search_word, p.required_word]
  if None in sw: sw.remove(None)
  print("Search words: %s \nExcluded words: %s\n" %(
     sw, p.excluded_words,
         ), file = out)
  print("%s  Line: %s" %(p.category, p.line_number), file = out)
  print(p.text_block, file = out)
  print(70*"-", file = out)

def select_problems(all_problems, key = 'header', value = None,
    n_context = 5):
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

def get_all_unique(all_problems, key = 'header'):
  all_unique = []
  for p in all_problems:
    x = getattr(p,key,None)
    if (x is not None) and not x in all_unique:
      all_unique.append(x)
  return all_unique

def pdb_write_statements(text, header = None, n_context = 5):
  all_problems = []
  for search_word in [".model_as_pdb(", ".as_pdb_string(",
     ".write_model_file("]:
    all_problems += display_context(header = header, text = text,
       n_context = n_context, search_word = search_word,quiet = True,
       category = 'pdb_write_statements')
  return all_problems

def pdb_format_interpretation(text, header = None, n_context = 5):
  all_problems = []
  for search_word in ['.open(']:
    for required_word in ['.pdb','pdb_file','model_file']:
      all_problems += display_context(header = header, text = text,
       n_context = n_context, 
       search_word = search_word,
       required_word=required_word,
       excluded_words=['transfer_ext','forward_compatible','iotbx.pdb'],
       quiet = True,
       category = 'pdb_format_interpretation')
  all_problems += display_context(header = header, text = text,
       n_context = n_context, 
       search_word = ".splitlines(", 
       required_word=" in ",
       excluded_words=["cmds.","iotbx.pdb"], # if this is present, it is ok
       quiet = True,
       category = 'pdb_format_interpretation')
  for search_word in ['HETATM','ATOM','TER','BREAK',' CA ',' N ']:
    for required_word in ['startswith','find(']:
      all_problems += display_context(header = header, text = text,
       n_context = n_context, 
       search_word = search_word,
       required_word=required_word,
       quiet = True,
       category = 'pdb_format_interpretation')
  return all_problems

def raw_records(text, header = None, n_context = 5):
  all_problems = []
  for search_word in ['raw_records=','lines=','raw_records =','lines =']:
    for required_word in ['.pdb','open(']:
      all_problems += display_context(header = header, text = text,
       n_context = n_context, 
       search_word = search_word,
       required_word=required_word,
       excluded_words=['pdb.input','iotbx.pdb','input.pdb'],
       quiet = True,
       category = 'pdb_format_interpretation')
  return all_problems

def pdb_file_name(text, header = None, n_context = 5):
  all_problems = []
  for search_word in ['file_name',' fn ',' fn=','filename']:
    for required_word in ['.pdb']:
      all_problems += display_context(header = header, text = text,
       n_context = n_context, 
       search_word = search_word,
       required_word=required_word,
       excluded_words = ['.pdb_or_mmcif_string_info',
          'pdb.input','iotbx.pdb','.cif','input.pdb'],
       quiet = True,
       category = 'pdb_format_interpretation')
  return all_problems


# =============================================================================
if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

# =============================================================================
# end
