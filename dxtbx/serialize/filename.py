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
  if path is None or path == "":
    return ""
  if directory is not None and not os.path.isabs(path):
    path = os.path.join(directory, path)
  return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
