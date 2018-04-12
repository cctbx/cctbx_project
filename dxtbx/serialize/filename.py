from __future__ import absolute_import, division, print_function
import os

def load_path(path, directory=None):
  ''' Load a filename from a JSON file.

  First expand any environment and user variables. Then create the absolute path
  from the current directory (i.e. the directory in which the JSON file is
  located.

  Params:
    path The path to the file.

  '''
  if not path:
    return ""
  path = os.path.expanduser(os.path.expandvars(path))
  if directory and not os.path.isabs(path):
    path = os.path.join(directory, path)
  return os.path.abspath(path)
