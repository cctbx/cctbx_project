from libtbx.path import norm_join
from libtbx.utils import UserError
import sys, os, pickle
from os.path import normpath, join, abspath, dirname, isdir, isfile
norm = normpath

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
  return norm_join(libtbx_build, "lib", "PYTHON_API_VERSION")

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

def _windows_pathext():
  result = os.environ.get("PATHEXT", "").lower().split(os.pathsep)
  for ext in [".py", ".exe", ".bat"]:
    if (ext not in result):
      result.insert(0, ext)
  return result

if (os.name == "nt"):
  windows_pathext = _windows_pathext()

class env:

  def __init__(self):
    libtbx_build = os.environ["LIBTBX_BUILD"]
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
    self._dispatcher_precall_commands = None
    self.partially_customized_windows_dispatcher = None
    self.windows_dispatcher_unique_pattern = None

  def write_api_file(self):
    api_file_name = python_api_version_file_name(self.LIBTBX_BUILD)
    if (not isdir(os.path.dirname(api_file_name))):
      os.makedirs(os.path.dirname(api_file_name))
    api_from_process = python_api_from_process()
    print >> open(api_file_name, "w"), api_from_process

  def dist_name(self, package_name):
    return package_name.upper() + "_DIST"

  def has_dist(self, package_name):
    return self.dist_name(package_name) in self.dist_paths

  def dist_path(self, package_name, default=KeyError):
    if (default is KeyError):
      return self.dist_paths[self.dist_name(package_name)]
    return self.dist_paths.get(self.dist_name(package_name), default)

  def effective_root(self, package_name):
    return dirname(self.dist_path(package_name))

  def current_working_directory_is_libtbx_build(self):
    return norm(os.getcwd()) == self.LIBTBX_BUILD

  def dispatcher_precall_commands(self):
    if (self._dispatcher_precall_commands is None):
      lines = []
      if (    self.python_version_major_minor == (2,2)
          and sys.platform == "linux2"
          and os.path.isfile("/etc/redhat-release")):
        try: red_hat_linux_release = open("/etc/redhat-release").readline()
        except: pass
        else:
          if (    red_hat_linux_release.startswith("Red Hat Linux release")
              and red_hat_linux_release.split()[4] == "9"):
            lines.extend([
              'if [ ! -n "$LD_ASSUME_KERNEL" ]; then',
              '  LD_ASSUME_KERNEL=2.4.1',
              '  export LD_ASSUME_KERNEL',
              'fi'])
      if (os.name == "posix" and self.compiler == "icc"):
        addl_lines = self.create_posix_icc_ld_preload()
        if (addl_lines is None):
          raise UserError("Cannot determine LD_PRELOAD for icc.")
        lines.extend(addl_lines)
      self._dispatcher_precall_commands = lines
    return self._dispatcher_precall_commands

  def create_posix_icc_ld_preload(self):
    path_icc = full_path("icc")
    if (path_icc is None): return None
    path_lib = os.sep.join(path_icc.split(os.sep)[:-2] + ["lib"])
    if (not isdir(path_lib)): return None
    ld_preload = []
    path_libirc_a = join(path_lib, "libirc.a")
    path_libirc_so = join(path_lib, "libirc.so")
    if (isfile(path_libirc_so)):
      ld_preload.append(path_libirc_so)
    else:
      if (isfile(path_libirc_a)):
        path_libirc_so = join(self.LIBTBX_BUILD, "lib/libirc.so")
        if (not isdir(os.path.dirname(path_libirc_so))):
          os.makedirs(os.path.dirname(path_libirc_so))
        if (not isfile(path_libirc_so)):
          cmd = "%s -shared -o %s %s" % (
            path_icc, path_libirc_so, path_libirc_a)
          print cmd
          sys.stdout.flush()
          os.system(cmd)
        ld_preload.append(path_libirc_so)
    path_libunwind_so = None
    best_version = None
    for file_name in os.listdir(path_lib):
      if (file_name.startswith("libunwind.so.")):
        try: version = int(file_name.split(".")[2])
        except: version = None
        if (version is not None):
          if (best_version is None or version > best_version):
            path_libunwind_so = join(path_lib, file_name)
            best_version = version
    if (path_libunwind_so is not None):
      ld_preload.append(path_libunwind_so)
    if (len(ld_preload) == 0): return None
    return [
      'LD_PRELOAD="%s"' % os.pathsep.join(ld_preload),
      'export LD_PRELOAD']

  def dispatcher_include(self):
    if (not hasattr(self, "_dispatcher_include")):
      file_name = norm(join(self.LIBTBX_BUILD, "dispatcher_include.sh"))
      if (isfile(file_name)):
        try: lines = open(file_name).read().splitlines()
        except IOError, e: raise UserError(str(e))
        lines.insert(0, "# included from %s" % file_name)
        m = max([len(line) for line in lines])
        lines.insert(0, "# " + "-"*(m-2))
        lines.append(lines[0])
        self._dispatcher_include = lines
      else:
        self._dispatcher_include = []
    return self._dispatcher_include

  def create_bin_sh_dispatcher(self, source_file, target_file):
    f = open(target_file, "w")
    print >> f, '#! /bin/sh'
    print >> f, '# LIBTBX_DISPATCHER DO NOT EDIT'
    print >> f, 'unset PYTHONHOME'
    print >> f, 'LIBTBX_BUILD="%s"' % self.LIBTBX_BUILD
    print >> f, 'export LIBTBX_BUILD'
    essentials = [("PYTHONPATH", self.PYTHONPATH)]
    if (sys.platform.startswith("darwin")):
      ld_library_path = "DYLD_LIBRARY_PATH"
    else:
      ld_library_path = "LD_LIBRARY_PATH"
    essentials.append((ld_library_path, self.LD_LIBRARY_PATH))
    essentials.append(("PATH", self.PATH))
    for n,v in essentials:
      v = ":".join(v)
      print >> f, 'if [ ! -n "$%s" ]; then' % n
      print >> f, '  %s="%s"' % (n, v)
      print >> f, 'else'
      print >> f, '  %s="%s:$%s"' % (n, v, n)
      print >> f, 'fi'
      print >> f, 'export %s' % n
    precall_commands = self.dispatcher_precall_commands()
    if (precall_commands is not None):
      for line in precall_commands:
        print >> f, line
    for line in self.dispatcher_include():
      print >> f, line
    cmd = ""
    if (source_file.lower().endswith(".py")):
      cmd += " '"+self.LIBTBX_PYTHON_EXE+"'"
    cmd += " '"+source_file+"'"
    print >> f, 'if [ -n "$LIBTBX__VALGRIND_FLAG__" ]; then'
    print >> f, "  exec $LIBTBX_VALGRIND"+cmd, '"$@"'
    print >> f, "elif [ $# -eq 0 ]; then"
    print >> f, "  exec"+cmd
    print >> f, "else"
    print >> f, "  exec"+cmd, '"$@"'
    print >> f, "fi"
    f.close()
    os.chmod(target_file, 0755)

  def windows_dispatcher(self, command_path,
        unique_pattern="0W6I0N6D0O2W8S5_0D0I8S1P4A3T6C4H9E4R7",
        libtbx_build="3L0I2B2T9B4X2_8B5U5I5L2D4",
        python_executable="5P2Y5T7H2O5N8_0E7X9E7C8U6T4A9B9L5E3",
        pythonpath="2P0Y1T7H3O2N7P7A2T5H8",
        main_path="1M5A1I0N4_8P7A0T9H9",
        target_command="5T4A3R7G8E3T7_6C5O0M0M3A8N8D2",
        dispatcher_exe_file_name="windows_dispatcher.exe"):
    if (os.name == "nt"
        and self.partially_customized_windows_dispatcher is None):
      try:
        self.partially_customized_windows_dispatcher = open(join(
          self.dist_path("libtbx"), dispatcher_exe_file_name), "rb").read()
      except IOError, e:
        raise RuntimeError(str(e))
      if (self.partially_customized_windows_dispatcher.find(unique_pattern)<0):
        raise RuntimeError('Unique pattern "%s" not found in file %s' % (
          unique_pattern, dispatcher_exe_file_name))
      self.windows_dispatcher_unique_pattern = unique_pattern
      for place_holder,actual_value in [
           (libtbx_build, self.LIBTBX_BUILD),
           (python_executable, self.LIBTBX_PYTHON_EXE),
           (pythonpath, os.pathsep.join(self.PYTHONPATH)),
           (main_path, os.pathsep.join([
                         norm(join(self.LIBTBX_BUILD, "bin")),
                         norm(join(self.LIBTBX_BUILD, "lib"))]))]:
       self.partially_customized_windows_dispatcher = patch_windows_dispatcher(
         dispatcher_exe_file_name=dispatcher_exe_file_name,
         binary_string=self.partially_customized_windows_dispatcher,
         place_holder=place_holder,
         actual_value=actual_value)
    if (command_path is None):
      return self.partially_customized_windows_dispatcher
    return patch_windows_dispatcher(
      dispatcher_exe_file_name=dispatcher_exe_file_name,
      binary_string=self.partially_customized_windows_dispatcher,
      place_holder=target_command,
      actual_value=command_path)

  def create_win32_dispatcher(self, source_file, target_file):
    open(target_file+".exe", "wb").write(
      self.windows_dispatcher(command_path=source_file))

  def create_dispatcher(self, source_file, target_file):
    if (os.name == "nt"):
      action = self.create_win32_dispatcher
      ext = ".exe"
    else:
      action = self.create_bin_sh_dispatcher
      ext = ""
      try: os.chmod(source_file, 0755)
      except OSError: pass
    target_file_ext = target_file + ext
    try: os.remove(target_file_ext)
    except OSError:
      try: os.remove(target_file_ext+".old")
      except OSError: pass
      try: os.rename(target_file_ext, target_file_ext+".old")
      except OSError: pass
    try: action(source_file, target_file)
    except IOError, e: print "  Ignored:", e

  def create_dispatcher_in_bin(self, source_file, target_file):
    self.create_dispatcher(
      source_file=source_file,
      target_file=norm(join(self.LIBTBX_BUILD, "bin", target_file)))

def patch_windows_dispatcher(
      dispatcher_exe_file_name,
      binary_string,
      place_holder,
      actual_value):
  place_holder_start = binary_string.find(place_holder)
  if (place_holder_start < 0):
    raise RuntimeError('Place holder "%s" not found in file %s' % (
      place_holder, dispatcher_exe_file_name))
  place_holder_end = binary_string.find("\0", place_holder_start)
  assert len(actual_value) <= place_holder_end - place_holder_start
  return binary_string[:place_holder_start] \
       + actual_value + "\0" \
       + binary_string[place_holder_start+len(actual_value)+1:]

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
