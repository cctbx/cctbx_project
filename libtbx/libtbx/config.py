from libtbx.path import norm_join
import sys, os, pickle
from os.path import normpath, join, abspath, dirname, isdir, isfile
norm = normpath

class UserError(Exception): pass

def full_path(command, search_first=[], search_last=[]):
  dirs = search_first + os.environ["PATH"].split(os.pathsep) + search_last
  for path in dirs:
    path_command = join(path, command)
    if (os.path.exists(path_command)):
      return norm(abspath(path_command))
  return None

def get_hostname():
  try: import socket
  except: return None
  try: return socket.gethostname()
  except: return None

def python_include_path(must_exist=True):
  if (sys.platform == "win32"):
    include_path = sys.prefix + r"\include"
  else:
    include_path = sys.prefix + "/include/python" + sys.version[0:3]
  include_path = norm_join(include_path)
  if (must_exist and not isdir(include_path)):
    raise RuntimeError("Cannot locate Python's include directory: %s"
      % include_path)
  return include_path

def python_api_from_process(must_exist=True):
  try: return str(sys.api_version) # Python 2.3 or higher
  except AttributeError: pass
  include_path = python_include_path(must_exist)
  if (not isdir(include_path)): return None
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

class resolve_redirection:

  def __init__(self, dist_root, name):
    self.effective_root = dist_root
    self.dist_path = norm_join(self.effective_root, name)
    if (isfile(self.dist_path)):
      try:
        redirection = open(self.dist_path).readlines()[0][:-1]
        assert len(redirection) > 0
      except:
        raise UserError("Error reading redirection file: %s" % self.dist_path)
      # first attempt: interpret the redirection as absolute path
      self.effective_root = norm(abspath(redirection))
      self.dist_path = norm_join(self.effective_root, name)
      if (not isdir(self.dist_path)):
        # second attempt: interpret the redirection as relative path
        self.effective_root = norm(abspath(join(dist_root, redirection)))
        self.dist_path = norm(join(self.effective_root, name))
        if (not isdir(self.dist_path)):
          raise UserError('Redirection to non-existing directory: "%s"'
            % self.effective_root)

adaptor_toolbox_suffix = "_adaptbx"

class package_pair:

  def __init__(self, name, dist_root=None):
    if (name.lower().endswith(adaptor_toolbox_suffix)):
      self.primary = name[:-len(adaptor_toolbox_suffix)]
      self.adaptbx = name
    else:
      self.primary = name
      self.adaptbx = name + adaptor_toolbox_suffix
    if (dist_root is not None):
      self.primary = resolve_redirection(
        dist_root=dist_root, name=self.primary).dist_path
      self.adaptbx = resolve_redirection(
        dist_root=dist_root, name=self.adaptbx).dist_path

  def primary_first(self):
    return (self.primary, self.adaptbx)

  def adaptbx_first(self):
    return (self.adaptbx, self.primary)

  def zip_primary_first(self, other):
    return ((self.primary, other.primary), (self.adaptbx, other.adaptbx))

  def zip_adaptbx_first(self, other):
    return ((self.adaptbx, other.adaptbx), (self.primary, other.primary))

def read_libtbx_config(path):
  try:
    f = open(path)
  except:
    raise UserError("Cannot open configuration file: " + path)
  try:
    result = eval(" ".join(f.readlines()), {}, {})
  except:
    raise UserError("Corrupt file: " + path)
  f.close()
  return result

class env:

  def __init__(self, application_prefix="LIBTBX"):
    self.application_prefix = application_prefix
    libtbx_build = os.environ["LIBTBX_BUILD"] # XXX XXX XXX XXX XXX XXX XXX XXX
    file_name = join(libtbx_build, "libtbx_env")
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

  def effective_root(self, package_name):
    return dirname(self.dist_path(package_name))

  def current_working_directory_is_libtbx_build(self):
    return norm(os.getcwd()) == self.LIBTBX_BUILD

  def dispatcher_front_end_exe(self,
       application_prefix_placeholder="3A7P0P0L4I8C3A5T5I6O3N490P0R7E6F8I6X4"):
    if (    os.name == "nt"
        and self._dispatcher_front_end_exe is None):
      self._dispatcher_front_end_exe = open(join(
        self.dist_path("libtbx"), "dispatcher_front_end.exe"), "rb").read()
      application_prefix_index = self._dispatcher_front_end_exe.find(
        application_prefix_placeholder)
      assert application_prefix_index >= 0
      assert len(self.application_prefix) <= len(application_prefix_placeholder)
      self._dispatcher_front_end_exe \
        = self._dispatcher_front_end_exe[:application_prefix_index] \
        + self.application_prefix + "\0" \
        + self._dispatcher_front_end_exe[
            application_prefix_index+len(self.application_prefix)+1:]
    return self._dispatcher_front_end_exe

class include_registry:

  def __init__(self):
    self.scan_boost()
    self._had_message = {}

  def scan_boost(self, flag=False):
    self._scan_boost = flag
    return self

  def scan_flag(self, path):
    if (not self._scan_boost and path.lower().endswith("boost")):
      if (not path in self._had_message):
        print "libtbx.scons: implicit dependency scan disabled for directory",
        print path
        self._had_message[path] = 1
      return False
    return True

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
