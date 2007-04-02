import sys
from fileinput import input, isfirstline, filename, isstdin
from string import expandtabs, rstrip
from libtbx import subversion

def clean_clutter_in(files):
  n_empty = 0
  for line in input(files, inplace=1):
    if (isfirstline()):
      if (not isstdin()):
        print >> sys.__stdout__, filename() + ':'
      n_empty = 0
    clean_line = line.expandtabs().rstrip()
    if (len(clean_line) == 0):
      n_empty += 1
    else:
      for i in xrange(n_empty): sys.stdout.write("\n")
      n_empty = 0
      sys.stdout.write(clean_line)
      sys.stdout.write("\n")

def run():
  files = []
  for arg in sys.argv[1:]:
    if arg == '--committing' or arg == '-c':
      if files:
        print_usage()
        exit(1)
      else:
        files = list(subversion.marked_for_commit())
    else:
      files.append(arg)
  if not files:
    print_usage()
    exit(1)
  clean_clutter_in(files)

def print_usage():
  print """
Usage:
  clean_clutter file1 file2 ...
  clean_clutter --committing|-c

The first form cleans the specified files whereas the second form cleans those
files which would be committed by running "svn commit"
"""

if (__name__ == "__main__"):
  run()
