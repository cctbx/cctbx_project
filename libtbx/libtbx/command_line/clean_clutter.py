import sys
from fileinput import input, isfirstline, filename, isstdin
from string import expandtabs, rstrip
n_empty = 0
for line in input(inplace=1):
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
