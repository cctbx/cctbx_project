from libtbx.str_utils import show_string
import os

op = os.path

def norm_join(*args):
  return op.normpath(op.join(*args))

def abs_real_norm(path):
  return op.normpath(op.realpath(op.abspath(path)))

def abs_norm(path):
  return op.normpath(op.abspath(path))

def create_target_dir(target_file):
  target_dir = op.split(target_file)[0]
  if (not op.isdir(target_dir)):
    os.makedirs(target_dir)

def canonical_path(file_name, effective_current_working_directory=None):
  if (not op.isabs(file_name)):
    if (effective_current_working_directory is None):
      effective_current_working_directory = os.getcwd()
    file_name = op.join(effective_current_working_directory, file_name)
  return op.normpath(file_name)

def is_same_canoncial_file(file_names):
  assert len(file_names) == 2
  if (file_names[0] == file_names[1]): return True
  if (hasattr(op, "samefile")):
    return op.samefile(file_names[0], file_names[1])
  return False

def is_same_file(file_names, effective_current_working_directory=None):
  return is_same_canoncial_file(
    [canonical_path(file_name, effective_current_working_directory)
      for file_name in file_names])

def full_command_path(command, search_first=[], search_last=[]):
  dirs = search_first + os.environ["PATH"].split(os.pathsep) + search_last
  for path in dirs:
    path_command = op.join(path, command)
    if (op.exists(path_command)):
      return abs_norm(path=path_command)
  return None

class directory(object):

  def __init__(self, path):
    self.path = path

  def get(self, name, must_exist=True):
    assert name is not None
    result = op.join(self.path, name)
    if (must_exist and not op.exists(result)):
      raise RuntimeError("No such file or directory: %s" % show_string(result))
    return result

  def sub_directory(self, name, must_exist=True):
    result = directory(self.get(name))
    if (must_exist and not op.isdir(result.path)):
      raise RuntimeError("Not a directory: %s" % show_string(result.path))
    return result

def walk_source_tree(top, arg=None):
  def visitor(result, dirname, names):
    names_keep = []
    for name in names:
      path = op.join(dirname, name)
      if (not op.isdir(path)):
        if (not name.endswith(".pyc")):
          result.append(path)
        continue
      def is_file_in_subdir(name):
        return op.isfile(op.join(path, name))
      if (   (name == "CVS" and is_file_in_subdir("Entries"))
          or (name == ".svn"
                and (   is_file_in_subdir("README.txt")
                     or is_file_in_subdir("entries")))):
        continue
      names_keep.append(name)
    if (len(names_keep) != len(names)):
      del names[:]
      names.extend(names_keep)
  result = []
  op.walk(top, visitor, result)
  return result
