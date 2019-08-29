from __future__ import absolute_import, division, print_function
import sys, os
from six.moves import zip

def run():
  results = []
  for file in sys.argv[1:]:
    p = os.popen("grep u+s " + file + " | cut -d: -f3", "r")
    results.append([float(x) for x in p.readlines()])
    p.close()
  assert len(results) == 2
  os.system("hostname")
  for x, y in zip(results[0], results[1]):
    print(x, y, x/y)

if (__name__ == "__main__"):
  run()
