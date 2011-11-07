import os
op = os.path

def norm_join(*args):
  return op.normpath(op.join(*args))

def abs_real_norm(path):
  return op.normpath(op.realpath(op.abspath(path)))

def abs_norm(path):
  return op.normpath(op.abspath(path))

def tail_levels(path, number_of_levels):
  return op.join(*path.split(op.sep)[-number_of_levels:])

def create_target_dir(target_file):
  target_dir = op.split(target_file)[0]
  if (not op.isdir(target_dir)):
    os.makedirs(target_dir)

def move_old(path, serial_sep="_", serial_fmt="%03d"):
  if (not op.exists(path)): return
  bns = op.basename(path) + serial_sep
  dn = op.dirname(op.abspath(path))
  max_i = 0
  for ex in os.listdir(dn):
    if (ex.startswith(bns)):
      s = ex[len(bns):]
      try: i = int(s)
      except ValueError: pass
      else: max_i = max(max_i, i)
  nn = op.join(dn, bns + serial_fmt % (max_i+1))
  os.rename(path, nn)

def move_old_create_new_directory(path, serial_sep="_", serial_fmt="%03d"):
  move_old(path=path, serial_sep=serial_sep, serial_fmt=serial_fmt)
  os.makedirs(path)

def canonical_path(file_name, effective_current_working_directory=None):
  if not isinstance(file_name, (str, unicode)):
    file_name = abs(file_name)
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
      from libtbx.str_utils import show_string
      raise RuntimeError("No such file or directory: %s" % show_string(result))
    return result

  def sub_directory(self, name, must_exist=True):
    result = directory(self.get(name))
    if (must_exist and not op.isdir(result.path)):
      from libtbx.str_utils import show_string
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

def random_new_directory_name(prefix="tmp_dir_", number_of_hex_code_digits=8):
  from libtbx.utils import random_hex_code
  for i_trial in xrange(10**6):
    name = prefix + random_hex_code(number_of_digits=number_of_hex_code_digits)
    if (not op.exists(name)):
      return name
  else:
    raise AssertionError

def makedirs_race(
      path=None,
      max_trials=10,
      delay_if_exists=0.001,
      delay_after_exception=0.5):
  if (path is None):
    path = random_new_directory_name()
  import time
  for i_trial in xrange(max_trials):
    if (op.exists(path)):
      if (delay_if_exists is not None):
        # in case the OS needs time to finalize makedirs from another process
        time.sleep(delay_if_exists)
      break
    try:
      os.makedirs(path)
    except Exception:
      if (delay_after_exception is not None):
        time.sleep(delay_after_exception)
  if (not op.isdir(path)):
    raise RuntimeError("makedirs_race(%s) failure." % path)
  return path


class path_mixin(object):

  def __truediv__(self, path):
    return self.__div__(path)

  def isdir(self):
    return os.path.isdir(abs(self))

  def isfile(self):
    return os.path.isfile(abs(self))

  def exists(self):
    return os.path.exists(abs(self))

  def open(self, *args, **kwds):
    return open(abs(self), *args, **kwds)

  def makedirs(self):
    os.makedirs(abs(self))

  def remove(self):
    if self.exists(): os.remove(abs(self))
    assert not self.exists()

  def remove_tree(self):
    from distutils.dir_util import remove_tree
    if self.isdir():
      remove_tree(abs(self))
    else:
      self.remove()

  def listdir(self):
    return os.listdir(abs(self))

  def chmod(self, *args, **kwds):
    return os.chmod(abs(self), *args, **kwds)

  def access(self, *args, **kwds):
    return os.access(abs(self), *args, **kwds)

  def basename(self):
    return os.path.basename(abs(self))

  def ext(self):
    return os.path.splitext(self.basename())[1]

  def split(self):
    return (self.dirname(), self.basename())

  def samefile(self, other):
    if isinstance(other, str) or isinstance(other, unicode):
      return os.path.samefile(abs(self), other)
    else:
      return os.path.samefile(abs(self), abs(other))


class absolute_path(path_mixin):

  def __init__(self, path, case_sensitive=False):
    assert os.path.isabs(path)
    if not case_sensitive:
      path = os.path.normcase(path)
    path = os.path.normpath(path)
    self._path = path

  def __div__(self, other):
    return absolute_path(os.path.join(self._path, other))

  def __abs__(self):
    return self._path

  def __add__(self, ext):
    return absolute_path(self._path + ext)

  def __repr__(self):
    return 'absolute_path("%s")' % self._path

  def dirname(self):
    return os.path.dirname(self._path)


class relocatable_path(path_mixin):

  def __init__(self, rooted, relocatable):
    self._rooted = rooted
    if os.path.isabs(relocatable):
      relocatable = os.path.relpath(relocatable, rooted.root_path)
    self.relocatable = relocatable

  def root(self):
    return self._rooted.root_path
  root = property(root)

  def __div__(self, path):
    return relocatable_path(self._rooted,
                            os.path.join(self.relocatable, path))

  def __idiv__(self, path):
    self.relocatable = os.path.join(self.relocatable, path)
    return self

  def __add__(self, ext):
    return relocatable_path(self._rooted, self.relocatable + ext)

  def self_or_abs_if(self, flag):
    if flag:
      return abs(self)
    else:
      return self

  def __abs__(self):
    return os.path.abspath(os.path.join(self.root, self.relocatable))

  def __repr__(self):
    return 'relocatable_path(root="%s", relocatable="%s")' % (self.root,
                                                              self.relocatable)

  def dirname(self):
    assert self.relocatable
    return relocatable_path(self._rooted, os.path.dirname(self.relocatable))

  def basename(self):
    return os.path.basename(self.relocatable)

  def normcase(self):
    return relocatable_path(self._rooted, os.path.normcase(self.relocatable))

  def __eq__(self, other):
    return (self._rooted == other._rooted
            and self.relocatable == other.relocatable)
