from __future__ import absolute_import, division, print_function
from contextlib import contextmanager
import os
import threading

temp_chdir_lock = threading.RLock()

@contextmanager
def temp_chdir(path):
  ''' A context manager to temporarily change the current directory, perform a
  task and then change back to the previous working directory. '''
  with temp_chdir_lock:
    cwd = os.getcwd()
    try:
      os.chdir(path)
      yield
    finally:
      os.chdir(cwd)

def load_path(path):
  ''' Load a filename from a JSON file.

  First expand any environment and user variables. Then create the absolute path
  from the current directory (i.e. the directory in which the JSON file is
  located.

  Params:
    path The path to the file.

  '''
  if path is None or path == "":
    return ""
  return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
