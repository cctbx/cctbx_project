from __future__ import absolute_import, division, print_function

import os
import sys


def recurse(parentformat, filename):
  for subformat in parentformat._children:
    understood = subformat.understand(filename)
    print("%s: %s" % (subformat.__name__, understood))
    if understood:
      recurse(subformat, filename)

def show_matching_formats(files):
  from dxtbx.format.Registry import Registry

  for filename in files:
    print('\n=== %s ===' % filename)
    if os.path.exists(filename):
      for fmt in Registry.get():
        understood = fmt.understand(filename)
        print("%s: %s" % (fmt.__name__, understood))
        if understood:
          recurse(fmt, filename)
    else:
      print("File not found.")

if __name__ == '__main__':
  show_matching_formats(sys.argv[1:])
