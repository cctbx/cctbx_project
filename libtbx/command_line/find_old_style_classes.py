from __future__ import absolute_import, division, print_function
import os, fnmatch
from libtbx import python_code_parsing

def run(args):
  if not args: args = [ '.' ]
  work = set()
  arg_filenames = []
  for arg in args:
    if os.path.isdir(arg):
      for dirpath, dirnames, filenames in os.walk(arg):
        work.update( os.path.join(dirpath, f)
                     for f in fnmatch.filter(filenames, '*.py') )
    else:
      arg_filenames.append(arg)
  work.update(fnmatch.filter(arg_filenames, '*.py'))
  for filename in work:
    try:
      old_style = python_code_parsing.find_old_style_classes(
        python_source_filename=filename)
    except Exception as e:
      import traceback
      print(filename)
      print(traceback.format_exc())
      continue
    if old_style:
      print('In file %s:' % filename)
      print(old_style)
      print()


if __name__ == '__main__':
  import sys
  run(args=sys.argv[1:])
