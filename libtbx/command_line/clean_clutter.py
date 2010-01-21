import sys
import os.path
from fileinput import input, isfirstline, filename, isstdin
from libtbx import subversion
from libtbx.option_parser import option_parser
import libtbx.command_line.file_clutter

def clean_clutter_in(files, tabsize=8):
  if not files: return
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
clean_clutter [-t n | --tabsize=n] file1 file2 ...
clean_clutter [-t n | --tabsize=n] [directory]
clean_clutter [-t n | --tabsize=n] [--committing|-c]""",
    description="""The first form cleans the specified files whereas the second
form cleans all files in the hierarchy rooted in the given directory or
the current directory is none is given.
The  -c options restricts cleaning to those files which would be committed
by running svn commit.""")
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
  else:
    if len(files) <= 1:
      if not files: dir = '.'
      else: dir = files[0]
      files = [ c.path for c in libtbx.command_line.file_clutter.gather([dir])
                if c.is_cluttered(flag_x=False) ]
  clean_clutter_in(files, tabsize=co.tabsize)

if (__name__ == "__main__"):
  import sys
  run()
