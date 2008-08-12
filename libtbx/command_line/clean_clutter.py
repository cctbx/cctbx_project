import sys
import os.path
from fileinput import input, isfirstline, filename, isstdin
from libtbx import subversion
from libtbx.option_parser import option_parser

def clean_clutter_in(files, tabsize=8):
  n_empty = 0
  for line in input([ f for f in files if not os.path.isdir(f) ], inplace=1):
    if (isfirstline()):
      if (not isstdin()):
        print >> sys.__stdout__, filename() + ':'
      n_empty = 0
    clean_line = line.expandtabs(tabsize).rstrip()
    if (len(clean_line) == 0):
      n_empty += 1
    else:
      for i in xrange(n_empty): sys.stdout.write("\n")
      n_empty = 0
      sys.stdout.write(clean_line)
      sys.stdout.write("\n")

def run():
  opt_parser = (option_parser(
    usage="""
clean_clutter [-t n | --tabsize==n] file1 file2 ...
clean_clutter [-t n | --tabsize==n] --committing|-c""",
    description="""The first form cleans the specified files whereas the second
form cleans those files which would be committed by running svn commit.""")
    .option("-t", "--tabsize",
      action="store",
      type="int",
      default=8,
      help="the number of spaces a tab is to be replaced by",
      metavar="INT")
     .option("-c", "--committing",
      action="store_true",
      default=False,
      help="whether to clean the files which are to be committed")
  )
  command_line = opt_parser.process(args=sys.argv[1:])
  co = command_line.options
  files = command_line.args
  if co.committing and files:
      opt_parser.show_help()
      exit(1)
  if co.committing:
    try:
      files = list(subversion.marked_for_commit())
    except RuntimeError, err:
      print err
      exit(1)
  elif not files:
    opt_parser.show_help()
    sys.exit(1)
  clean_clutter_in(files, tabsize=co.tabsize)

if (__name__ == "__main__"):
  import sys
  run()
