from __future__ import generators
import libtbx.path
from libtbx.optparse_wrapper import option_parser
from libtbx import introspection
from libtbx.utils import UserError
import pickle
from cStringIO import StringIO
import sys, os

def get_hostname():
  try: import socket
  except: return None
  try: return socket.gethostname()
  except: return None

def get_mipspro_version():
  version = os.popen4("CC -version", "r")[1].read().strip().split()
  if (version[:3] == "MIPSpro Compilers: Version ".split()):
    if (version[3].startswith("7.3")):
      return "73"
    elif (version[3].startswith("7.4")):
      return "74"
  sys.tracebacklimit = 0
  raise RuntimeError("Unknown MIPSpro compiler (CC -version).")

def python_include_path(must_exist=True):
  if (sys.platform == "win32"):
    include_path = sys.prefix + r"\include"
  else:
    include_path = sys.prefix + "/include/python" + sys.version[0:3]
  include_path = libtbx.path.norm_join(include_path)
  if (must_exist and not os.path.isdir(include_path)):
    raise RuntimeError("Cannot locate Python's include directory: %s"
      % include_path)
  return include_path

def python_api_from_process(include_must_exist=True):
  try: return str(sys.api_version) # Python 2.3 or higher
  except AttributeError: pass
  include_path = python_include_path(must_exist=include_must_exist)
  if (not os.path.isdir(include_path)): return "UNKNOWN"
  modsupport_h = open(libtbx.path.norm_join(
    include_path, "modsupport.h")).readlines()
  python_api_version = None
  for line in modsupport_h:
    if (line.startswith("#define")):
      flds = line.split()
      if (len(flds) == 3 and flds[1] == "PYTHON_API_VERSION"):
        python_api_version = flds[2]
        break
  assert python_api_version is not None
  return python_api_version

def ld_library_path_var_name():
  if (os.name == "nt"):
    return "PATH"
  if (sys.platform.startswith("darwin")):
    return "DYLD_LIBRARY_PATH"
  else:
    return "LD_LIBRARY_PATH"

def highlight_dispatcher_include_lines(lines):
  m = max([len(line) for line in lines])
  lines.insert(0, "# " + "-"*(m-2))
  lines.append(lines[0])

def source_specific_dispatcher_include(pattern, source_file):
  try: source_lines = open(source_file).read().splitlines()
  except IOError: return []
  lines = ["# lines marked " + pattern]
  for line in source_lines:
    pattern_begin = line.find(pattern)
    if (pattern_begin >= 0):
      pattern_end = pattern_begin + len(pattern)
      to_include = line[pattern_end:]
      if (not to_include[:1].isalnum()):
        if (to_include.startswith(" ")):
          to_include = to_include[1:]
        lines.append(to_include)
  highlight_dispatcher_include_lines(lines)
  return lines

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

def write_incomplete_libtbx_environment(f):
  message = [
    "*******************************************",
    "Fatal Error: Incomplete libtbx environment!",
    "*******************************************",
    "Please re-run the libtbx/configure.py command."]
  if (os.name != "nt"):
    for line in message: print >> f, '  echo "%s"' % line
    print >> f, '  echo ""'
  else:
    for line in message: print >> f, '  echo %s' % line
    print >> f, '  echo.'

class common_setpaths:

  def __init__(self, env, shell, suffix):
    self.env = env
    self.shell = shell
    self.suffix = suffix
    self.s = open_info(env.under_build("setpaths%s.%s" % (suffix, shell)))
    if (suffix == "_debug"):
      self.u = open_info(env.under_build("unsetpaths.%s" % shell))
    else:
      self.u = StringIO() # /dev/null equivalent

  def all_and_debug(self):
    if (self.suffix != ""):
      self.setenv("LIBTBX_BUILD", self.env.build_path)
      for module in self.env.module_list:
        for name,path in module.name_and_dist_path_pairs():
          self.setenv(name.upper()+"_DIST", path)
    if (self.suffix == "_debug"):
      self.update_path("PYTHONPATH", self.env.pythonpath)
      self.update_path(ld_library_path_var_name(), [self.env.lib_path])

class unix_setpaths(common_setpaths):

  def __init__(self, env, shell, suffix):
    assert shell in ["sh", "csh"]
    common_setpaths.__init__(self, env, shell, suffix)
    if (self.shell == "sh"):
      self._setenv = "%s="
    else:
      self._setenv = "setenv %s "

  def setenv(self, var_name, val):
    if (self.shell == "sh"):
      print >> self.s, '  %s="%s"' % (var_name, val)
      print >> self.s, '  export %s' % var_name
      print >> self.u, '  unset %s' % var_name
    else:
      print >> self.s, '  setenv %s "%s"' % (var_name, val)
      print >> self.u, '  unsetenv %s' % var_name

  def update_path(self, var_name, val):
    val = os.pathsep.join(val)
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      print >> f, '''  %s"`'%s' '%s' %s %s '%s'`"''' % (
        self._setenv % var_name,
        self.env.python_exe,
        self.env.path_utility,
        action,
        var_name,
        val)
      if (f is self.s and self.shell == "sh"):
        print >> f, '  export %s' % var_name
      if (self.shell == "sh"):
        print >> f, '  if [ "$%s" == "E_M_P_T_Y" ]; then unset %s; fi' % (
          var_name, var_name)
      else:
        print >> f, '  if ("$%s" == "E_M_P_T_Y") unsetenv %s' % (
          var_name, var_name)

class windows_setpaths(common_setpaths):

  def __init__(self, env, suffix):
    common_setpaths.__init__(self, env, "bat", suffix)

  def setenv(self, var_name, val):
    print >> self.s, '  set %s=%s' % (var_name, val)
    print >> self.u, '  set %s=' % var_name

  def update_path(self, var_name, val):
    val = os.pathsep.join(val)
    fmt = '''\
  for /F "delims=" %%%%i in ('%s "%s" %s %s "%s"') do set %s=%%%%i'''
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      print >> f, fmt % (
        self.env.python_exe,
        self.env.path_utility,
        action,
        var_name,
        val,
        var_name)
      print >> f, '  if "%%%s%%" == "E_M_P_T_Y" set %s=' % (
        var_name, var_name)

def _windows_pathext():
  result = os.environ.get("PATHEXT", "").lower().split(os.pathsep)
  for ext in [".py", ".exe", ".bat"]:
    if (ext not in result):
      result.insert(0, ext)
  return result

if (os.name == "nt"):
  windows_pathext = _windows_pathext()

def open_info(path, mode="w", info="Creating:"):
  print info, path
  try: return open(path, mode)
  except IOError, e:
    raise UserError(str(e))

def remove_or_rename(path):
  try: os.remove(path)
  except OSError:
    try: os.remove(path+".old")
    except OSError: pass
    try: os.rename(path, path+".old")
    except OSError: pass

class environment:

  def __init__(self, build_path=None):
    self.python_version_major_minor = sys.version_info[:2]
    self.build_path = build_path
    self._shortpath_bat = None
    self.build_path = self.abs_path_clean(build_path)
    self.manage_python_api_version()
    self.reset_dispatcher_support()
    self.set_derived_paths()
    if (os.name == "nt"):
      self.python_exe = self.abs_path_clean(os.path.join(
        sys.prefix, "python.exe"))
    else:
      self.python_exe = self.abs_path_clean(os.path.join(
        sys.prefix, "bin/python"))
    self.read_command_version_suffix()
    self.build_options = None
    self.repository_paths = []
    self.reset_module_registry()
    self.scons_dist_path = None
    self.pythonpath = []

  def manage_python_api_version(self):
    self.python_api_version = python_api_from_process(include_must_exist=False)
    path = os.path.join(self.build_path, "lib")
    if (not os.path.isdir(path)):
      os.makedirs(path)
    path = os.path.join(path, "PYTHON_API_VERSION")
    if (not os.path.isfile(path)):
      prev_python_api_version = "UNKNOWN"
    else:
      prev_python_api_version = open(path).read().strip()
      if (prev_python_api_version != self.python_api_version):
        if (prev_python_api_version != "UNKNOWN"):
          raise UserError(
            "Incompatible Python API's\n"
            "  Current version:        %s\n"
            "  Used to build binaries: %s" % (
              self.python_api_version,
              prev_python_api_version))
    if (prev_python_api_version != self.python_api_version):
      print >> open(path, "w"), self.python_api_version

  def reset_dispatcher_support(self):
    self._shortpath_bat = None
    self._dispatcher_precall_commands = None
    self.partially_customized_windows_dispatcher = None
    self.windows_dispatcher_unique_pattern = None

  def reset_module_registry(self):
    self.module_list = []
    self.module_dict = {}
    self.module_dist_paths = {}
    self.missing_for_build = {}
    self.missing_for_use = {}
    self.missing_optional = {}

  def abs_path_short(self, abs_path):
    if (os.name != "nt"): return abs_path
    if (self._shortpath_bat is None):
      self._shortpath_bat = self.under_build("shortpath.bat")
      assert os.path.exists(self._shortpath_bat)
    return os.popen('call "%s" "%s"' %
      (self._shortpath_bat, abs_path), "r").readline().rstrip()

  def abs_path_clean(self, path):
    abs_path = os.path.normpath(os.path.abspath(path))
    if (os.name != "nt" or abs_path.find(" ") < 0): return abs_path
    short = self.abs_path_short(abs_path).split(os.sep)
    orig = abs_path.split(os.sep)
    clean = []
    for o,s in zip(orig, short):
      if (o.find(" ") < 0):
        clean.append(o)
      else:
        clean.append(s)
    return os.sep.join(clean)

  def set_derived_paths(self):
    if (self.build_path is None):
      self.bin_path = None
      self.exe_path = None
      self.lib_path = None
      self.include_path = None
    else:
      self.bin_path = self.under_build("bin")
      self.exe_path = self.under_build("exe")
      self.lib_path = self.under_build("lib")
      self.include_path = self.under_build("include")

  def under_build(self, path):
    return libtbx.path.norm_join(self.build_path, path)

  def under_dist(self, module_name, path):
    return libtbx.path.norm_join(self.module_dist_paths[module_name], path)

  def dist_path(self, module_name, default=KeyError):
    if (default is KeyError):
      return self.module_dist_paths[module_name]
    return self.module_dist_paths.get(module_name, default)

  def dist_paths(self):
    for module in self.module_list:
      for dist_path in module.dist_paths_active():
        yield dist_path

  def clear_bin_directory(self):
    if (not os.path.isdir(self.bin_path)): return
    for file_name in os.listdir(self.bin_path):
      path = os.path.join(self.bin_path, file_name)
      if (os.path.isfile(path)):
        remove_or_rename(path)

  def write_command_version_suffix(self):
    assert self.command_version_suffix is not None
    path = self.under_build("command_version_suffix")
    try: f = open(path, "w")
    except IOError:
      raise UserError("Cannot write command_version_suffix file: " + path)
    print >> f, self.command_version_suffix

  def read_command_version_suffix(self):
    path = self.under_build("command_version_suffix")
    if (not os.path.isfile(path)):
        self.command_version_suffix = None
    else:
      try:
        self.command_version_suffix = open(path).read().strip()
      except IOError:
        raise UserError("Cannot read command_version_suffix file: " + path)

  def register_module(self, dependent_module, module):
    if (dependent_module is None):
      self.module_list.append(module)
    else:
      for dependent_index,registered_module in enumerate(self.module_list):
        if (registered_module is dependent_module): break
      else: raise RuntimeError(
        "Internal error: dependent_module not in module_list.")
      self.module_list.insert(dependent_index, module)
    self.module_dict[module.name] = module
    for name,path in module.name_and_dist_path_pairs():
      self.module_dist_paths[name] = path

  def find_dist_path(self, module_name, optional=False):
    for path in self.repository_paths:
      dist_path = self.abs_path_clean(libtbx.path.norm_join(path, module_name))
      if (os.path.isdir(dist_path)):
        return dist_path
    if (not optional):
      msg = ["Module not found: %s" % module_name,
             "  Repository directories searched:"]
      for path in self.repository_paths:
        msg.append("    " + path)
      raise UserError("\n".join(msg))
    return None

  def process_module(self, dependent_module, module_name, optional):
    dist_path = self.find_dist_path(module_name, optional=optional)
    if (dist_path is None): return False
    new_module = module(env=self, name=module_name, dist_path=dist_path)
    new_name_normcase = os.path.normcase(new_module.name)
    for name in self.module_dict:
      if (os.path.normcase(name) == new_name_normcase): return True
    new_module.find_mate()
    new_module.process_libtbx_config()
    self.register_module(dependent_module=dependent_module, module=new_module)
    new_module.process_dependencies()
    return True

  def add_repository(self, path):
    path = self.abs_path_clean(path)
    path_normcase = os.path.normcase(path)
    for repository_path in self.repository_paths:
      if (os.path.normcase(repository_path) == path_normcase):
        break
    else:
      self.repository_paths.append(path)

  def process_args(self, args, default_repository=None):
    cold_start = default_repository is not None
    if (len(args) == 0): args = ["--help"]
    if (default_repository is None):
      command_name = "libtbx.configure"
    else:
      command_name = "libtbx/configure.py"
    parser = option_parser(usage="%s [options] module_name ..." % command_name)
    if (not cold_start):
      if ("--only" in args): self.reset_module_registry()
      parser.option(None, "--only",
        action="store_true",
        help="disable previously configured modules")
    else:
      parser.option("-r", "--repository",
        action="callback",
        type="string",
        callback=self.option_repository,
        help="path to source code repository",
        metavar="DIRECTORY")
      parser.option(None, "--build",
        choices=("release", "quick", "debug"),
        default="release",
        help="build mode",
        metavar="release|quick|debug")
      parser.option(None, "--compiler",
        action="store",
        type="string",
        default="default",
        help="select non-standard compiler",
        metavar="STRING")
      parser.option(None, "--static_libraries",
        action="store_true",
        default=False,
        help="build all libraries statically")
      parser.option(None, "--static_exe",
        action="store_true",
        default=False,
        help="link all executables statically (implies --static_libraries)")
      parser.option(None, "--scan_boost",
        action="store_true",
        default=False,
        help="enable implicit dependency scan")
      parser.option(None, "--command_version_suffix",
        action="store",
        type="string",
        default=None,
        help="version suffix for commands in bin directory",
        metavar="STRING")
    command_line = parser.process(args=args)
    if (default_repository is not None):
      self.add_repository(default_repository)
    if (cold_start):
      self.build_options = build_options(
        compiler=command_line.options.compiler,
        mode=command_line.options.build,
        static_libraries=command_line.options.static_libraries,
        static_exe=command_line.options.static_exe,
        scan_boost=command_line.options.scan_boost)
      if (command_line.options.command_version_suffix is not None):
        self.command_version_suffix = \
          command_line.options.command_version_suffix
        self.write_command_version_suffix()
    module_names = list(command_line.args)
    module_names.append("libtbx")
    module_names.reverse()
    for module_name in module_names:
      self.process_module(
        dependent_module=None, module_name=module_name, optional=False)
    self.scons_dist_path = self.find_dist_path("scons", optional=True)
    self.path_utility = self.under_dist(
      "libtbx", "libtbx/command_line/path_utility.py")

  def option_repository(self, option, opt, value, parser):
    if (not os.path.isdir(value)):
      raise UserError("Not a directory: --repository %s" % value)
    self.add_repository(value)

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
      if (os.name == "posix" and self.build_options.compiler == "icc"):
        addl_lines = self.create_posix_icc_ld_preload()
        if (addl_lines is None):
          raise UserError("Cannot determine LD_PRELOAD for icc.")
        lines.extend(addl_lines)
      self._dispatcher_precall_commands = lines
    return self._dispatcher_precall_commands

  def create_posix_icc_ld_preload(self):
    path_icc = libtbx.path.full_command_path("icc")
    if (path_icc is None): return None
    path_lib = os.sep.join(path_icc.split(os.sep)[:-2] + ["lib"])
    if (not os.path.isdir(path_lib)): return None
    ld_preload = []
    path_libirc_a = os.path.join(path_lib, "libirc.a")
    path_libirc_so = os.path.join(path_lib, "libirc.so")
    if (os.path.isfile(path_libirc_so)):
      ld_preload.append(path_libirc_so)
    else:
      if (os.path.isfile(path_libirc_a)):
        path_libirc_so = self.under_build("lib/libirc.so")
        if (not os.path.isdir(os.path.dirname(path_libirc_so))):
          os.makedirs(os.path.dirname(path_libirc_so))
        if (not os.path.isfile(path_libirc_so)):
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
            path_libunwind_so = os.path.join(path_lib, file_name)
            best_version = version
    if (path_libunwind_so is not None):
      ld_preload.append(path_libunwind_so)
    if (len(ld_preload) == 0): return None
    return [
      'LD_PRELOAD="%s"' % os.pathsep.join(ld_preload),
      'export LD_PRELOAD']

  def dispatcher_include(self):
    if (not hasattr(self, "_dispatcher_include")):
      self._dispatcher_include = []
      for file_name in os.listdir(self.build_path):
        if (not os.path.isfile(file_name)): continue
        if (    file_name.startswith("dispatcher_include")
            and file_name.endswith(".sh")):
          try: lines = open(file_name).read().splitlines()
          except IOError, e: raise UserError(str(e))
          lines.insert(0, "# included from %s" % file_name)
          highlight_dispatcher_include_lines(lines)
          self._dispatcher_include.extend(lines)
    return self._dispatcher_include

  def write_bin_sh_dispatcher(self, source_file, target_file):
    f = open(target_file, "w")
    print >> f, '#! /bin/sh'
    print >> f, '# LIBTBX_DISPATCHER DO NOT EDIT'
    print >> f, 'unset PYTHONHOME'
    print >> f, 'LIBTBX_BUILD="%s"' % self.build_path
    print >> f, 'export LIBTBX_BUILD'
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((ld_library_path_var_name(), [self.lib_path]))
    essentials.append(("PATH", [self.bin_path]))
    for n,v in essentials:
      if (len(v) == 0): continue
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
    for line in source_specific_dispatcher_include(
                  pattern="LIBTBX_PRE_DISPATCHER_INCLUDE_SH",
                  source_file=source_file):
      print >> f, line
    for line in self.dispatcher_include():
      print >> f, line
    for line in source_specific_dispatcher_include(
                  pattern="LIBTBX_POST_DISPATCHER_INCLUDE_SH",
                  source_file=source_file):
      print >> f, line
    cmd = ""
    if (source_file.lower().endswith(".py")):
      cmd += " '"+self.python_exe+"'"
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
        self.partially_customized_windows_dispatcher = open(self.under_dist(
          "libtbx", dispatcher_exe_file_name), "rb").read()
      except IOError, e:
        raise RuntimeError(str(e))
      if (self.partially_customized_windows_dispatcher.find(unique_pattern)<0):
        raise RuntimeError('Unique pattern "%s" not found in file %s' % (
          unique_pattern, dispatcher_exe_file_name))
      self.windows_dispatcher_unique_pattern = unique_pattern
      for place_holder,actual_value in [
           (libtbx_build, self.build_path),
           (python_executable, self.python_exe),
           (pythonpath, os.pathsep.join(self.pythonpath)),
           (main_path, os.pathsep.join([self.bin_path, self.lib_path]))]:
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

  def write_win32_dispatcher(self, source_file, target_file):
    open(target_file+".exe", "wb").write(
      self.windows_dispatcher(command_path=source_file))

  def write_dispatcher(self, source_file, target_file):
    if (os.name == "nt"):
      action = self.write_win32_dispatcher
      ext = ".exe"
    else:
      action = self.write_bin_sh_dispatcher
      ext = ""
      try: os.chmod(source_file, 0755)
      except OSError: pass
    target_file_ext = target_file + ext
    remove_or_rename(target_file_ext)
    try: action(source_file, target_file)
    except IOError, e: print "  Ignored:", e

  def write_dispatcher_in_bin(self, source_file, target_file):
    self.write_dispatcher(
      source_file=source_file,
      target_file=self.under_build("bin/"+target_file))

  def write_setpaths_sh(self, suffix):
    setpaths = unix_setpaths(self, "sh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, 'libtbx_pythonhome_save="$PYTHONHOME"'
      print >> f, 'unset PYTHONHOME'
      print >> f, '"%s" -V > /dev/null 2>&1' % self.python_exe
      print >> f, 'if [ $? -ne 0 -o ! -f "%s" ]; then' % self.path_utility
      write_incomplete_libtbx_environment(f)
      print >> f, 'else'
    setpaths.update_path("PATH", [self.bin_path])
    for command in ["setpaths_all", "unsetpaths"]:
      print >> s, """  alias libtbx.%s='. "%s/%s.sh"'""" % (
       command, self.build_path, command)
    print >> u, '  unalias libtbx.unsetpaths > /dev/null 2>&1'
    setpaths.all_and_debug()
    for f in s, u:
      print >> f, 'fi'
      print >> f, 'if [ -n "$libtbx_pythonhome_save" ]; then'
      print >> f, '  PYTHONHOME="$libtbx_pythonhome_save"'
      print >> f, '  export PYTHONHOME'
      print >> f, 'fi'
      print >> f, 'unset libtbx_pythonhome_save'

  def write_setpaths_csh(self, suffix):
    setpaths = unix_setpaths(self, "csh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, 'if ($?PYTHONHOME) then'
      print >> f, '  set libtbx_pythonhome_save="$PYTHONHOME"'
      print >> f, '  unsetenv PYTHONHOME'
      print >> f, 'endif'
      print >> f, '"%s" -V >& /dev/null' % self.python_exe
      print >> f, 'if ($status != 0 || ! -f "%s") then' % self.path_utility
      write_incomplete_libtbx_environment(f)
      print >> f, 'else'
    setpaths.update_path("PATH", [self.bin_path])
    for command in ["setpaths_all", "unsetpaths"]:
      print >> s, """  alias libtbx.%s 'source "%s/%s.csh"'""" % (
       command, self.build_path, command)
    print >> u, '  unalias libtbx.unsetpaths'
    setpaths.all_and_debug()
    for f in s, u:
      print >> f, 'endif'
      print >> f, 'if ($?libtbx_pythonhome_save) then'
      print >> f, '  setenv PYTHONHOME "$libtbx_pythonhome_save"'
      print >> f, '  unset libtbx_pythonhome_save'
      print >> f, 'endif'

  def write_setpaths_bat(self, suffix):
    setpaths = windows_setpaths(self, suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, '@ECHO off'
      print >> f, 'if not exist "%s" goto fatal_error' % self.python_exe
      print >> f, 'if exist "%s" goto update_path' % self.path_utility
      print >> f, ':fatal_error'
      write_incomplete_libtbx_environment(f)
      print >> f, '  goto end_of_script'
      print >> f, ':update_path'
      print >> f, '  set PYTHONHOME='
    setpaths.update_path("PATH", [self.bin_path])
    for command in ["setpaths_all", "unsetpaths"]:
      print >> s, '  doskey libtbx.%s=%s\\%s.bat $*' % (
        command, self.build_path, command)
    print >> u, '  doskey libtbx.unsetpaths='
    setpaths.all_and_debug()
    if (suffix == "_debug"):
      print >> s, '  set PYTHONCASEOK=1' # no unset
    for f in s, u:
      print >> f, ':end_of_script'

  def write_SConstruct(self):
    SConstruct_path = self.under_build("SConstruct")
    f = open_info(SConstruct_path)
    print >> f, 'SConsignFile()'
    for path in self.repository_paths:
      print >> f, 'Repository(r"%s")' % path
    for module in self.module_list:
      name,path  = list(module.name_and_dist_path_pairs())[-1]
      if (os.path.isfile(os.path.join(path, "SConscript"))):
        print >> f, 'SConscript("%s/SConscript")' % name
    f.close()

  def collect_test_scripts(self):
    result = []
    for module in self.module_list:
      result.extend(module.collect_test_scripts())
    return result

  def write_run_tests_csh(self):
    test_scripts = self.collect_test_scripts()
    if (len(test_scripts) > 0):
      path = self.under_build("run_tests.csh")
      f = open_info(path)
      print >> f, "#! /bin/csh -f"
      print >> f, "set noglob"
      print >> f, "set verbose"
      for file_name in test_scripts:
        print >> f, 'libtbx.python "%s" $*' % file_name
      f.close()
      os.chmod(path, 0755)

  def pickle(self):
    self.reset_dispatcher_support()
    file_name = libtbx.path.norm_join(self.build_path, "libtbx_env")
    pickle.dump(self, open(file_name, "wb"), 0)

  def write_setpath_files(self):
    print "Python:", sys.version.split()[0], sys.executable
    if (len(self.missing_for_build) == 0):
      self.build_options.report()
    print "command_version_suffix:", self.command_version_suffix
    print "Top-down list of all packages involved:"
    top_down_module_list = list(self.module_list)
    top_down_module_list.reverse()
    for module in top_down_module_list:
      print " ", "+".join(module.names_active())
    if (len(self.missing_for_use) > 0):
      raise UserError("Missing modules: "
        + " ".join(self.missing_for_use.keys()))
    if (len(self.missing_for_build) != 0):
      if (self.scons_dist_path is not None):
        print "***********************************"
        print "Warning: modules missing for build:"
        for module_name in self.missing_for_build.keys():
          print " ", module_name
        print "***********************************"
    for suffix in ["", "_all", "_debug"]:
      if (hasattr(os, "symlink")):
        self.write_setpaths_sh(suffix)
        self.write_setpaths_csh(suffix)
      else:
        self.write_setpaths_bat(suffix)
    self.pickle()

  def write_python_and_show_path_duplicates(self):
    module_names = {}
    for file_name in os.listdir(self.bin_path):
      if (file_name.startswith(".")): continue
      file_name_lower = file_name.lower()
      if (file_name_lower.startswith("libtbx.")): continue
      if (   file_name_lower == "python"
          or file_name_lower.startswith("python.")): continue
      module_names[file_name.split(".")[0]] = None
    module_names = module_names.keys()
    for module_name in module_names:
      self.write_dispatcher_in_bin(
        source_file=self.python_exe,
        target_file=module_name+".python")
    for command in ["show_build_path", "show_dist_paths"]:
      source_file = self.under_dist(
        "libtbx", "libtbx/command_line/"+command+".py")
      for module_name in module_names:
        self.write_dispatcher_in_bin(
          source_file=source_file,
          target_file=module_name+"."+command)

  def write_command_version_duplicates(self):
    if (self.command_version_suffix is None): return
    suffix = "_" + self.command_version_suffix
    for file_name in os.listdir(self.bin_path):
      if (file_name.startswith(".")): continue
      source_file = os.path.join(self.bin_path, file_name)
      if (os.name == "nt" and file_name.lower().endswith(".exe")):
        target_file = source_file[:-4] + suffix + source_file[-4:]
      else:
        target_file = source_file + suffix
      remove_or_rename(target_file)
      try: open(target_file, "wb").write(open(source_file, "rb").read())
      except IOError: pass

  def assemble_pythonpath(self):
    pythonpath = [self.lib_path]
    for module in self.module_list:
      pythonpath.extend(module.assemble_pythonpath())
    pythonpath.reverse()
    hash = {}
    self.pythonpath = []
    for path in pythonpath:
      path_normcase = os.path.normcase(path)
      if (path_normcase in hash): continue
      hash[path_normcase] = None
      self.pythonpath.append(path)

  def refresh(self):
    self.assemble_pythonpath()
    self.write_setpath_files()
    if (len(self.missing_for_build) == 0):
      self.write_SConstruct()
    if (os.name != "nt"):
      self.write_run_tests_csh()
    self.clear_bin_directory()
    if (not os.path.isdir(self.bin_path)):
      os.makedirs(self.bin_path)
    for file_name in ("libtbx.python", "python"):
      self.write_dispatcher_in_bin(
        source_file=self.python_exe,
        target_file=file_name)
    for module in self.module_list:
      module.process_command_line_directories()
    if (os.path.isdir(self.exe_path)):
      print "Processing:", self.exe_path
      for file_name in os.listdir(self.exe_path):
        if (file_name[0] == "."): continue
        self.write_dispatcher_in_bin(
          source_file=libtbx.path.norm_join(self.exe_path, file_name),
          target_file=file_name)
    self.write_python_and_show_path_duplicates()
    self.write_command_version_duplicates()

class module:

  def __init__(self, env, name, dist_path, mate_suffix="_adaptbx"):
    self.env = env
    self.mate_suffix = mate_suffix
    if (os.path.normcase(name).endswith(os.path.normcase(self.mate_suffix))):
      self.name = name[:-len(self.mate_suffix)]
      self.names = [self.name, name]
      self.dist_paths = [None, dist_path]
    else:
      self.name = name
      self.names = [name, name + self.mate_suffix]
      self.dist_paths = [dist_path, None]

  def names_active(self):
    for name,path in zip(self.names, self.dist_paths):
      if (path is not None): yield name

  def dist_paths_active(self):
    for path in self.dist_paths:
      if (path is not None): yield path

  def name_and_dist_path_pairs(self, all=False):
    for name,path in zip(self.names, self.dist_paths):
      if (all or path is not None): yield (name,path)

  def find_mate(self):
    new_dist_paths = []
    for name,path in self.name_and_dist_path_pairs(all=True):
      if (path is None):
        path = self.env.find_dist_path(module_name=name, optional=True)
      new_dist_paths.append(path)
    self.dist_paths = new_dist_paths

  def process_libtbx_config(self):
    self.python_paths = []
    self.required_for_build = []
    self.required_for_use = []
    self.optional = []
    dist_paths = []
    for dist_path in self.dist_paths:
      if (dist_path is not None):
        while True:
          path = libtbx.path.norm_join(dist_path, "libtbx_config")
          if (not os.path.isfile(path)):
            config = None
            break
          try: f = open(path)
          except IOError: raise UserError(
            "Cannot open configuration file: " + path)
          try: config = eval(" ".join(f.readlines()), {}, {})
          except KeyboardInterrupt: raise
          except: raise UserError("Corrupt configuration file: " + path)
          f.close()
          redirection = config.get("redirection", None)
          if (redirection is None):
            break
          if (not isinstance(redirection, str)):
            raise UserError(
              "Corrupt configuration file:\n"
              "  file = %s\n"
              "  redirection must be a Python string" % path)
          new_dist_path = os.path.expandvars(redirection)
          if (not os.path.isabs(new_dist_path)):
            new_dist_path = libtbx.path.norm_join(dist_path, new_dist_path)
          new_dist_path = self.env.abs_path_clean(new_dist_path)
          if (not os.path.isdir(new_dist_path)):
            raise UserError(
              "Invalid redirection:\n"
              "  file = %s\n"
              "  redirection = %s\n"
              "  resulting target = %s" % (path, redirection, new_dist_path))
          dist_path = new_dist_path
        if (config is not None):
          self.required_for_build.extend(config.get(
            "modules_required_for_build", []))
          self.required_for_use.extend(config.get(
            "modules_required_for_use", []))
          self.optional.extend(config.get(
            "optional_modules", []))
      dist_paths.append(dist_path)
    self.dist_paths = dist_paths

  def process_dependencies(self):
    for module_name in self.required_for_build:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_for_build[module_name] = None
    for module_name in self.required_for_use:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_for_use[module_name] = None
    for module_name in self.optional:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_optional[module_name] = None

  def assemble_pythonpath(self):
    result = []
    for dist_path in self.dist_paths_active():
      for sub_dir in ["", "/"+self.name]:
        sub_file = sub_dir + "/__init__.py"
        path = os.path.normpath(dist_path + sub_file)
        if (os.path.isfile(path)):
          result.append(os.path.dirname(os.path.dirname(path)))
          break
    return result

  def write_dispatcher(self, source_dir, file_name):
    source_file = libtbx.path.norm_join(source_dir, file_name)
    if (not os.path.isfile(source_file)): return
    file_name_lower = file_name.lower()
    if (file_name_lower.startswith("__init__.py")): return
    if (file_name_lower.endswith(".pyc")): return
    if (file_name_lower.endswith(".pyo")): return
    if (file_name[0] == "."): return
    if (os.name == "nt"):
      ext = os.path.splitext(file_name_lower)[1]
      if (ext not in windows_pathext): return
    target_file = self.name
    if (file_name_lower != "main.py"):
      target_file += "." + os.path.splitext(file_name)[0]
    self.env.write_dispatcher_in_bin(
      source_file=source_file,
      target_file=target_file)

  def process_command_line_directories(self):
    for dist_path in self.dist_paths_active():
      for source_dir in [
            libtbx.path.norm_join(dist_path, "command_line"),
            libtbx.path.norm_join(dist_path, self.name, "command_line")]:
        if (not os.path.isdir(source_dir)): continue
        print "Processing:", source_dir
        for file_name in os.listdir(source_dir):
          self.write_dispatcher(source_dir=source_dir, file_name=file_name)

  def collect_test_scripts(self,
        file_names=["run_tests.py", "run_examples.py"]):
    result = []
    for dist_path in self.dist_paths_active():
      for file_name in file_names:
        path = libtbx.path.norm_join(dist_path, file_name)
        if (os.path.isfile(path)): result.append(path)
    return result

class build_options:

  def __init__(self, compiler, mode, static_libraries, static_exe, scan_boost):
    introspection.adopt_init_args()
    assert self.mode in ["release", "quick", "debug"]
    self.optimization = (self.mode == "release")
    self.debug_symbols = (self.mode == "debug")

  def report(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Compiler:", self.compiler
    print >> f, "Build mode:", self.mode
    print >> f, "Static libraries:", self.static_libraries
    print >> f, "Static exe:", self.static_exe
    print >> f, "Scan Boost headers:", self.scan_boost

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

def cold_start(args):
  env = environment(build_path=os.path.normpath(os.path.abspath(os.getcwd())))
  env.process_args(
    args=args[1:],
    default_repository=os.path.dirname(os.path.dirname(args[0])))
  env.refresh()

def unpickle():
  file_name = libtbx.path.norm_join(os.environ["LIBTBX_BUILD"], "libtbx_env")
  env = pickle.load(open(file_name, "rb"))
  if (env.python_version_major_minor != sys.version_info[:2]):
    raise UserError("Python version incompatible with this build.\n"
     + "  Version used to configure: %d.%d\n" % env.python_version_major_minor
     + "  Version used now: %d.%d" % sys.version_info[:2])
  return env

def warm_start(args):
  env = unpickle()
  env.process_args(args=args[1:])
  env.refresh()
