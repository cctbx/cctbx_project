"""
Tools to safely add and remove files
"""
from __future__ import absolute_import, division, print_function
import random
import stat
import os
import os.path as op
from shutil import rmtree
from six.moves import range

def make_paths_writable_if_possible(paths):
  """Safely make a list of paths writable (if possible)"""
  os_stat = os.stat
  os_chmod = os.chmod
  irwusr = stat.S_IRUSR | stat.S_IWUSR
  irwxusr = irwusr | stat.S_IXUSR
  for path in paths:
    try:
      mode = os_stat(path).st_mode
    except KeyboardInterrupt: raise
    except Exception:
      pass
    else:
      if (op.isdir(path)):
        new_mode = mode | irwxusr
      else:
        new_mode = mode | irwusr
      if (new_mode != mode):
        try:
          os_chmod(path, new_mode)
        except KeyboardInterrupt: raise
        except Exception:
          pass

def remove_directories_if_possible(paths):
  """Safely remove a list of directories (if possible)"""
  remaining = []
  for path in paths:
    if (op.isdir(path)):
      try:
        rmtree(path, ignore_errors=True)
      except KeyboardInterrupt: raise
      except Exception:
        pass
    if (op.isdir(path)):
      remaining.append(path)
  return remaining

def remove_files_if_possible(paths):
  """Safely remove a list of files (if possible)"""
  remaining = []
  for path in paths:
    if (op.exists(path)):
      try:
        os.remove(path)
      except KeyboardInterrupt: raise
      except Exception:
        pass
    if (op.exists(path)):
      remaining.append(path)
  return remaining

def rename_files_and_directories_if_possible(paths):
  """Safely rename a list of files or directories and append _OBSOLETE_xxx
    (if possible)"""
  remaining = []
  for path in paths:
    if (op.exists(path)):
      path = op.normpath(path) # strips trailing slash
      for i_trial in range(1000):
        new_path = path + "_OBSOLETE_%8.8X" % random.randrange(2**32)
        if (not op.exists(new_path)):
          try:
            os.rename(path, new_path)
          except KeyboardInterrupt: raise
          except Exception:
            pass
          break
    if (op.exists(path)):
      remaining.append(path)
  return remaining

def remove_or_rename_files_and_directories_if_possible(paths):
  """Safely remove a list of files and directories,
    and if not successful, rename them with the extension _OBSOLETE_xxx
    (if possible)"""
  make_paths_writable_if_possible(paths=paths)
  remaining_dirs = remove_directories_if_possible(paths=paths)
  remaining_files = remove_files_if_possible(paths=paths)
  remaining = []
  for remaining_paths in [remaining_dirs, remaining_files]:
    remaining.extend(
      rename_files_and_directories_if_possible(paths=remaining_paths))
  return remaining
