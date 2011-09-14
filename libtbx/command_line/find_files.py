from libtbx.path import walk_source_tree
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from libtbx.option_parser import option_parser
from fnmatch import fnmatch
import re
import sys, os

def read_lines_if_possible(file_path):
  try: f = open(file_path, "r")
  except IOError: return []
  return f.read().splitlines()

def run(args, command_name="libtbx.find_files"):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s [options] pattern ..." % command_name,
    description="Recursively finds all files matching patterns,\n"
      "excluding CVS and .svn directories and .pyc files.")
    .option("-t", "--top",
      action="append",
      type="string",
      metavar="PATH",
      help="top-level directory where search starts"
           " (default is current working directory)")
    .option("-g", "--grep",
      action="append",
      type="string",
      metavar="PATTERN",
      help="find regular expression pattern in each file (multiple"
           " -g/--grep options can be given)")
    .option("-i", "--ignore_case",
      action="store_true",
      default=False,
      help="with -g/--grep: case-insensitive match")
    .option("-f", "--file_names_only",
      action="store_true",
      default=False,
      help="with -g/--grep: show file names only, not the matching lines")
    .option("-q", "--quote",
      action="store_true",
      default=False,
      help="quote file names")
  ).process(args=args)
  fn_patterns = command_line.args
  co = command_line.options
  grep_flags = 0
  if (co.ignore_case):
    grep_flags |= re.IGNORECASE
  if (len(fn_patterns) == 0):
    fn_patterns = ["*"]
  tops = co.top
  if (tops is None):
    tops = ["."]
  for top in tops:
    if (not os.path.isdir(top)):
      raise Sorry("Not a directory: %s" % show_string(top))
    for file_path in walk_source_tree(top=top):
      file_name = os.path.basename(file_path)
      for fn_pattern in fn_patterns:
        if (fnmatch(file_name, fn_pattern)):
          if (co.quote): fp = show_string(file_path)
          else: fp = file_path
          if (co.grep is None):
            print fp
          else:
            is_binary_file = co.file_names_only
            for line in read_lines_if_possible(file_path=file_path):
              if (not is_binary_file):
                is_binary_file = "\0" in line
              def line_matches_all_grep_patterns():
                for grep_pattern in co.grep:
                  if (re.search(
                        pattern=grep_pattern,
                        string=line,
                        flags=grep_flags) is None):
                    return False
                return True
              if (line_matches_all_grep_patterns()):
                if (co.file_names_only):
                  print fp
                  break
                elif (is_binary_file):
                  print "%s: match in binary file" % fp
                  break
                else:
                  print "%s: %s" % (fp, line)

if (__name__ == "__main__"):
  run(sys.argv[1:])
