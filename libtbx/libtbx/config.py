import sys, os, pickle
from libtbx.path import norm_join

class UserError(Exception): pass

def full_path(command, search_first=[], search_last=[]):
  dirs = search_first + os.environ["PATH"].split(os.pathsep) + search_last
  for path in dirs:
    path_command = os.path.join(path, command)
    if (os.path.exists(path_command)):
      return os.path.normpath(os.path.abspath(path_command))
  return None

def get_hostname():
  try: import socket
  except: return None
  try: return socket.gethostname()
  except: return None

def python_include_path(must_exist=0001):
  if (sys.platform == "win32"):
    include_path = sys.prefix + r"\include"
  else:
    include_path = sys.prefix + "/include/python" + sys.version[0:3]
  include_path = norm_join(include_path)
  if (must_exist and not os.path.isdir(include_path)):
    raise RuntimeError("Cannot locate Python's include directory: %s"
      % include_path)
  return include_path

def python_api_from_process(must_exist=0001):
  try: return str(sys.api_version) # Python 2.3 or higher
  except AttributeError: pass
  include_path = python_include_path(must_exist)
  if (not os.path.isdir(include_path)): return None
  modsupport_h = open(norm_join(include_path, "modsupport.h")).readlines()
  python_api_version = None
  for line in modsupport_h:
    if (line.startswith("#define")):
      flds = line.split()
      if (len(flds) == 3 and flds[1] == "PYTHON_API_VERSION"):
        python_api_version = flds[2]
        break
  assert python_api_version is not None
  return python_api_version

def python_api_version_file_name(libtbx_build):
  return norm_join(libtbx_build, "libtbx", "PYTHON_API_VERSION")

adaptor_toolbox_suffix = "_adaptbx"

class package_pair:

  def __init__(self, name):
    if (name.lower().endswith(adaptor_toolbox_suffix)):
      self.primary = name[:-len(adaptor_toolbox_suffix)]
      self.adaptbx = name
    else:
      self.primary = name
      self.adaptbx = name + adaptor_toolbox_suffix

  def primary_first(self):
    return (self.primary, self.adaptbx)

  def adaptbx_first(self):
    return (self.adaptbx, self.primary)

  def zip_primary_first(self, other):
    return ((self.primary, other.primary), (self.adaptbx, other.adaptbx))

  def zip_adaptbx_first(self, other):
    return ((self.adaptbx, other.adaptbx), (self.primary, other.primary))

class env:

  def __init__(self):
    libtbx_build = os.environ["LIBTBX_BUILD"]
    file_name = os.path.join(libtbx_build, "libtbx_env")
    try:
      f = open(file_name)
    except:
      raise UserError("Cannot open libtbx environment file: " + file_name)
    try:
      e = pickle.load(f)
    except:
      raise UserError("Corrupt libtbx environment file: " + file_name)
    self.__dict__.update(e)
    assert self.LIBTBX_BUILD == libtbx_build
    if (not self.python_version_major_minor == sys.version_info[:2]):
      raise UserError("Python version incompatible with this build."
       + " Version used to configure: %d.%d." % self.python_version_major_minor
       + " Version used now: %d.%d." % sys.version_info[:2])
    self._dispatcher_front_end_exe = None

  def write_api_file(self):
    api_file_name = python_api_version_file_name(self.LIBTBX_BUILD)
    api_from_process = python_api_from_process()
    print >> open(api_file_name, "w"), api_from_process

  def dist_name(self, package_name):
    return package_name.upper() + "_DIST"

  def has_dist(self, package_name):
    return self.dist_name(package_name) in self.dist_paths

  def dist_path(self, package_name):
    return self.dist_paths[self.dist_name(package_name)]

  def current_working_directory_is_libtbx_build(self):
    return os.path.normpath(os.getcwd()) == self.LIBTBX_BUILD

  def dispatcher_front_end_exe(self):
    if (    os.name == "nt"
        and self._dispatcher_front_end_exe is None):
      self._dispatcher_front_end_exe = open(os.path.join(
        self.dist_path("libtbx"), "dispatcher_front_end.exe"), "rb").read()
    return self._dispatcher_front_end_exe

class include_registry:

  def __init__(self):
    self.scan_boost()
    self._had_message = {}

  def scan_boost(self, flag=00000):
    self._scan_boost = flag
    return self

  def scan_flag(self, path):
    if (not self._scan_boost and path.lower().endswith("boost")):
      if (not path in self._had_message):
        print "libtbx.scons: implicit dependency scan disabled for directory",
        print path
        self._had_message[path] = 1
      return 00000
    return 0001

  def prepend_include_switch(self, env, path):
    assert isinstance(path, str)
    return env["INCPREFIX"] + path

  def append(self, env, paths):
    assert isinstance(paths, list)
    for path in paths:
      if (self.scan_flag(path)):
        env.Append(CPPPATH=[path])
      else:
        ipath = self.prepend_include_switch(env, path)
        env.Append(CXXFLAGS=[ipath])
        env.Append(SHCXXFLAGS=[ipath])

  def prepend(self, env, paths):
    assert isinstance(paths, list)
    paths.reverse()
    for path in paths:
      if (self.scan_flag(path)):
        env.Prepend(CPPPATH=[path])
      else:
        ipath = self.prepend_include_switch(env, path)
        env.Prepend(CXXFLAGS=[ipath])
        env.Prepend(SHCXXFLAGS=[ipath])

def get_mipspro_version():
  version = os.popen4("CC -version", "r")[1].read().strip().split()
  if (version[:3] == "MIPSpro Compilers: Version ".split()):
    if (version[3].startswith("7.3")):
      return "73"
    elif (version[3].startswith("7.4")):
      return "74"
  sys.tracebacklimit = 0
  raise RuntimeError("Unknown MIPSpro compiler (CC -version).")

class _build_options:

  def set(self, **kw):
    assert env().current_working_directory_is_libtbx_build()
    self.__dict__.update(kw)

build_options = _build_options()
