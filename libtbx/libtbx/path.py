import os

def norm(path):
  return os.path.normcase(os.path.normpath(path))

def norm_join(*args):
  return norm(apply(os.path.join, args))

def create_target_dir(target_file):
  target_dir = os.path.split(target_file)[0]
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)

def canonical_path(file_name, effective_current_working_directory=None):
  if (not os.path.isabs(file_name)):
    if (effective_current_working_directory is None):
      effective_current_working_directory = os.getcwd()
    file_name = os.path.join(effective_current_working_directory, file_name)
  return norm(file_name)

def is_same_canoncial_file(file_names):
  assert len(file_names) == 2
  if (file_names[0] == file_names[1]): return True
  if (hasattr(os.path, "samefile")):
    return os.path.samefile(file_names[0], file_names[1])
  return False

def is_same_file(file_names, effective_current_working_directory=None):
  return is_same_canoncial_file(
    [canonical_path(file_name, effective_current_working_directory)
      for file_name in file_names])

def full_command_path(command, search_first=[], search_last=[]):
  dirs = search_first + os.environ["PATH"].split(os.pathsep) + search_last
  for path in dirs:
    path_command = os.path.join(path, command)
    if (os.path.exists(path_command)):
      return norm(os.path.abspath(path_command))
  return None
