import libtbx.path
from libtbx.path import relocatable_path, absolute_path
from libtbx.str_utils import show_string
from libtbx.utils import detect_binary_file
from libtbx import adopt_init_args
import platform
import shutil
import pickle
from cStringIO import StringIO
import re
import sys, os
op = os.path

if os.environ.get('LIBTBX_WINGIDE_DEBUG'):
  import wingdbstub # special import

# XXX backward compatibility 2011-03-29
qnew = 2
if (qnew == 1 and sys.version_info[:2] < (2,7)): qnew = 0
qnew = ["", " -Qwarn", " -Qnew"][qnew]

if (os.name == "nt"):
  exe_suffix = ".exe"
else:
  exe_suffix = ""

def bool_literal(b):
  if (b): return "True"
  return "False"

default_write_full_flex_fwd_h = sys.platform.startswith("irix")
default_msvc_arch_flag = ["None", "SSE2"][int(os.name == "nt")]
default_build_boost_python_extensions = True
default_enable_boost_threads = False
default_enable_openmp_if_possible = False
default_enable_cuda = False
default_opt_resources = False

def is_64bit_architecture():
  return (platform.architecture()[0] == "64bit")

def unique_paths(paths):
  hash = set()
  result = []
  for path in paths:
    try: path_normcase = abs(path.normcase())
    except AttributeError: path_normcase = op.normcase(path)
    if (path_normcase in hash): continue
    hash.add(path_normcase)
    result.append(path)
  return result

def darwin_shlinkcom(env_etc, env, lo, dylib):
  if env_etc.compiler.startswith('darwin_'):
    if (env_etc.mac_cpu == "powerpc" or env_etc.compiler == "darwin_gcc"):
      dylib1 = "-ldylib1.o"
    else :
      dylib1 = " ".join(env_etc.shlinkflags)
    if (    env_etc.gcc_version is not None
        and env_etc.gcc_version >= 40201):
      opt_m = ""
    else:
      opt_m = " -m"
    shlinkcom = [
      "ld -dynamic%s -r -d -bind_at_load -o %s $SOURCES" % (opt_m, lo),
      "$SHLINK -nostartfiles -undefined dynamic_lookup -Wl,-dylib"
        " %s -o %s %s" % (dylib1, dylib, lo)]
    if (env_etc.mac_os_use_dsymutil):
      shlinkcom.append('dsymutil "%s"' % dylib)
    env.Replace(SHLINKCOM=shlinkcom)

def get_darwin_gcc_build_number(gcc='gcc'):
  from libtbx import easy_run
  gcc_version = (easy_run.fully_buffered(command='%s --version' % gcc)
                         .raise_if_errors()
                         .stdout_lines[0].strip())
  m = re.search(r"Apple Inc. build (\d+)\)", gcc_version)
  if m is None: return None
  try:
    return int(m.group(1))
  except ValueError:
    return None

def is_llvm_compiler(gcc='gcc') :
  from libtbx import easy_run
  try :
    gcc_version = (easy_run.fully_buffered(command='%s --version' % gcc)
                           .raise_if_errors()
                           .stdout_lines[0].strip())
    return ("llvm" in gcc_version)
  except Exception, e :
    return False

def get_gcc_version(command_name="gcc"):
  from libtbx import easy_run
  buffer = easy_run.fully_buffered(
    command="%s -dumpversion" % command_name)
  if (len(buffer.stderr_lines) != 0):
    return None
  if (len(buffer.stdout_lines) != 1):
    return None
  major_minor_patchlevel = buffer.stdout_lines[0].split(".")
  if (len(major_minor_patchlevel) not in [2,3]):
    return None
  num = []
  for fld in major_minor_patchlevel:
    try: i = int(fld)
    except ValueError:
      return None
    num.append(i)
  if (len(num) == 2): num.append(0) # substitute missing patchlevel
  return ((num[0]*100)+num[1])*100+num[2]

def get_hostname():
  try: import socket
  except KeyboardInterrupt: raise
  except Exception: return None
  try: return socket.gethostname()
  except KeyboardInterrupt: raise
  except Exception: return None

def get_ldd_output(target=None):
  if (target is None): target = sys.executable
  from libtbx import easy_run
  return easy_run.go(command="ldd '%s'" % target).stdout_lines

def get_hp_ux_acc_version():
  from libtbx import easy_run
  run_out = easy_run.go(command="aCC -V").stdout_lines
  version = run_out
  if (len(version) > 0):
    version = version[0].strip().split()
  # aCC: HP aC++/ANSI C B3910B A.06.01 [Jan 05 2005]
  if (len(version) >= 6 and version[5].startswith("A.")):
    return version[5]
  raise RuntimeError(
    "\n  ".join(["Unknown C++ compiler (aCC -V):"] + run_out))

def get_mipspro_version():
  from libtbx import easy_run
  run_out = easy_run.go(command="CC -version").stdout_lines
  version = run_out
  if (len(version) > 0):
    version = version[0].strip().split()
  if (version[:3] == "MIPSpro Compilers: Version ".split()):
    if (version[3].startswith("7.3")):
      return "73"
    elif (version[3].startswith("7.4")):
      return "74"
  raise RuntimeError(
    "\n  ".join(["Unknown MIPSpro compiler (CC -version):"] + run_out))

def python_include_path(must_exist=True):
  if (sys.platform == "win32"):
    include_path = sys.prefix + r"\include"
  else:
    include_path = sys.prefix + "/include/python%d.%d" % sys.version_info[:2]
  include_path = libtbx.path.norm_join(include_path)
  if (must_exist and not op.isdir(include_path)):
    raise RuntimeError("Cannot locate Python's include directory: %s"
      % include_path)
  return include_path

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
  try: source_lines = source_file.open().read().splitlines()
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

def write_do_not_edit(f, win_bat=False):
  if (win_bat): s = "rem"
  else:         s = "#"
  print >> f, s+' THIS IS AN AUTOMATICALLY GENERATED FILE.'
  print >> f, s+' DO NOT EDIT! CHANGES WILL BE LOST.'

def open_info(path, mode="w", info="   "):
  print info, path.basename()
  try: return path.open(mode)
  except IOError, e:
    raise RuntimeError(str(e))

class common_setpaths(object):

  def __init__(self, env, shell, suffix):
    self.env = env
    self.shell = shell
    self.suffix = suffix
    self.s = open_info(env.under_build("setpaths%s.%s" % (suffix, shell),
                                       return_relocatable_path=True))
    if (suffix == "_debug"):
      self.u = open_info(env.under_build("unsetpaths.%s" % shell,
                                         return_relocatable_path=True))
    else:
      self.u = StringIO() # /dev/null equivalent

  def all_and_debug(self):
    if (self.suffix == "_debug"):
      self.update_path("PYTHONPATH",
        os.pathsep.join([
          self.path_script_value(_) for _ in self.env.pythonpath]))
      self.update_path(
        ld_library_path_var_name(), self.path_script_value(self.env.lib_path))

  def set_unset_vars(self):
    if (self.suffix != ""):
      for var_name,path in self.env.var_name_and_build_or_dist_path_pairs():
        self.setenv(var_name=var_name, val=self.path_script_value(path))

class unix_setpaths(common_setpaths):

  def __init__(self, env, shell, suffix):
    assert shell in ["sh", "csh"]
    common_setpaths.__init__(self, env, shell, suffix)
    if (self.shell == "sh"):
      self._setenv = "%s="
    else:
      self._setenv = "setenv %s "

  def path_script_value(self, path_obj):
    return path_obj.sh_value()

  def setenv(self, var_name, val):
    if (self.shell == "sh"):
      print >> self.s, '%s="%s"' % (var_name, val)
      print >> self.s, 'export %s' % var_name
      print >> self.u, 'unset %s' % var_name
    else:
      print >> self.s, 'setenv %s "%s"' % (var_name, val)
      print >> self.u, 'unsetenv %s' % var_name

  def update_path(self, var_name, val, var_name_in=None):
    if (var_name_in is None): var_name_in = var_name
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      if (self.shell == "sh"):
        print >> f, '''%s`libtbx.path_utility %s %s "%s" < /dev/null`''' % (
          self._setenv % var_name, action, var_name_in, val)
      else:
        print >> f, '''%s"`libtbx.path_utility %s %s '%s' < /dev/null`"''' % (
          self._setenv % var_name, action, var_name_in, val)
      if (f is self.s and self.shell == "sh"):
        print >> f, 'export %s' % var_name
      if (self.shell == "sh"):
        print >> f, \
          'if [ "$%s" = "L_I_B_T_B_X_E_M_P_T_Y" ]; then unset %s; fi' % (
          var_name, var_name)
      else:
        print >> f, 'if ("$%s" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv %s' % (
          var_name, var_name)

class windows_setpaths(common_setpaths):

  def __init__(self, env, suffix):
    common_setpaths.__init__(self, env, "bat", suffix)

  def path_script_value(self, path_obj):
    return path_obj.bat_value()

  def setenv(self, var_name, val):
    print >> self.s, '@set %s=%s' % (var_name, val)
    print >> self.u, '@set %s=' % var_name

  def update_path(self, var_name, val, var_name_in=None):
    if (var_name_in is None): var_name_in = var_name
    fmt = '''\
@for /F "delims=" %%%%i in ('libtbx.path_utility %s %s "%s"') do @set %s=%%%%i'''
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      print >> f, fmt % (
        action,
        var_name_in,
        val,
        var_name)
      print >> f, '@if "%%%s%%" == "L_I_B_T_B_X_E_M_P_T_Y" @set %s=' % (
        var_name, var_name)

def _windows_pathext():
  result = os.environ.get("PATHEXT", "").lower().split(os.pathsep)
  for ext in [".py", ".exe", ".bat"]:
    if (ext not in result):
      result.insert(0, ext)
  return result

if (os.name == "nt"):
  windows_pathext = _windows_pathext()

def remove_or_rename(path):
  assert path is not None
  if (not isinstance(path, str)):
    path = abs(path)
  if (op.isfile(path)):
    try: os.remove(path)
    except OSError:
      try: os.remove(path+".old")
      except OSError: pass
      try: os.rename(path, path+".old")
      except OSError: pass

class environment:

  def __init__(self, build_path):
    self.python_version_major_minor = sys.version_info[:2]
    self.set_build_path(build_path)
    self.manage_python_version_major_minor()
    self.reset_dispatcher_support()
    self.set_derived_paths()
    self.python_exe = relocatable_path(self, sys.executable)
    if self.python_exe.relocatable.startswith('..'):
      self.python_exe = absolute_path(sys.executable)
    # sanity checks
    assert self.python_exe.isfile()
    assert self.python_exe.access(os.X_OK)

    self.read_command_version_suffix()
    self.build_options = None
    self.repository_paths = []
    self.command_line_redirections = {}
    self.reset_module_registry()
    self.scons_dist_path = None
    self.pythonpath = []
    self.relocatable = True

  def raise_python_version_incompatible(self, prev_pvmm=None):
    if (prev_pvmm is None):
      prev_pvmm = "%d.%d" % self.python_version_major_minor
    raise RuntimeError("Python version incompatible with this build:\n"
      + "  Build directory: %s\n" % show_string(self.build_path)
      + "  Python version used initially: %s\n" % prev_pvmm
      + "  Python version in use now:     %d.%d" % sys.version_info[:2])

  def manage_python_version_major_minor(self):
    path = self.build_path / "lib"
    if not path.isdir():
      path.makedirs()
    path /= "PYTHON_VERSION_MAJOR_MINOR"
    pvmm = "%d.%d" % self.python_version_major_minor
    if not path.isfile():
      path.open("w").write("""\
# DO NOT EDIT THIS FILE UNDER ANY CIRCUMSTANCE!
# The version number below is purely to insure against accidental use
# of another Python version when re-configuring an existing build.
# To use a different Python version, initialize a new build directory
# with the libtbx/configure.py command.
%s
""" % pvmm)
    else:
      prev_pvmm = path.open().read().splitlines()
      if (len(prev_pvmm) == 0):
        prev_pvmm = None
      else:
        prev_pvmm = prev_pvmm[-1].strip()
      if (prev_pvmm != pvmm):
        self.raise_python_version_incompatible(prev_pvmm=prev_pvmm)

  def reset_dispatcher_support(self):
    self._dispatcher_precall_commands = None

  def reset_module_registry(self):
    self.module_list = []
    self.module_dict = {}
    self.module_dist_paths = {}
    self.missing_for_build = set()
    self.missing_for_use = set()
    self.missing_optional = set()

  def is_ready_for_build(self):
    return (len(self.missing_for_build) == 0)

  def as_relocatable_path(self, path):
    if isinstance(path, libtbx.path.path_mixin): return path
    return relocatable_path(self, path)

  def set_build_path(self, build_path):
    build_path = op.realpath(op.normcase(op.normpath(build_path)))
    d, b = op.split(build_path)
    self.root_path = d
    self.build_path = relocatable_path(self, b)

  def set_derived_paths(self):
    self.bin_path     = self.build_path / 'bin'
    self.exe_path     = self.build_path / 'exe'
    self.lib_path     = self.build_path / 'lib'
    self.include_path = self.build_path / 'include'

  def under_build(self, path, return_relocatable_path=False):
    result = self.build_path / path
    if return_relocatable_path:
      return result
    else:
      return abs(result)

  def under_dist(self, module_name, path, default=KeyError, test=None,
                 return_relocatable_path=False):
    if (default is KeyError):
      result = self.module_dist_paths[module_name] / path
    else:
      mdp = self.module_dist_paths.get(module_name)
      if (mdp is None): return default
      result = mdp / path
    if (test is None or test(abs(result))):
      if return_relocatable_path:
        return result
      else:
        return abs(result)
    return None

  def dist_path(self, module_name, default=KeyError,
                return_relocatable_path=False):
    if (default is KeyError):
      result = self.module_dist_paths[module_name]
    else:
      result = self.module_dist_paths.get(module_name, default)
    if (isinstance(result, relocatable_path) and not return_relocatable_path):
      result = abs(result)
    return result

  def has_module(self, name):
    return self.module_dist_paths.has_key(name)

  def require_module(self, name, error=RuntimeError):
    if (not self.has_module(name)):
      build_path = self.build_path
      from libtbx import introspection
      nproc = introspection.number_of_processors()
      raise error("""\
The %(name)s module is needed but not configured.
Run:
  cd %(build_path)s
  libtbx.configure %(name)s
  libtbx.scons -j %(nproc)d
Wait for the command to finish, then try again.""" % vars())

  def dist_paths(self, return_relocatable_path=False):
    for module in self.module_list:
      for dist_path in module.dist_paths_active():
        if return_relocatable_path:
          yield dist_path
        else:
          yield abs(dist_path)

  def var_name_and_build_or_dist_path_pairs(self):
    yield ("LIBTBX_BUILD", self.build_path)
    for module in self.module_list:
      for name,path in module.name_and_dist_path_pairs():
        yield (name.upper()+"_DIST", path)

  def set_os_environ_all_dist(self):
    for module in self.module_list:
      for name,path in module.name_and_dist_path_pairs():
        os.environ[name.upper()+"_DIST"] = abs(path)

  def clear_bin_directory(self):
    if not self.bin_path.isdir(): return
    buffer = []
    have_libtbx_command = False
    for file_name in self.bin_path.listdir():
      if (    not have_libtbx_command
          and file_name.lower().startswith("libtbx.")):
        have_libtbx_command = True
      path = self.bin_path / file_name
      if path.isfile():
        buffer.append(path)
    if (len(buffer) != 0):
      if (not have_libtbx_command):
        raise RuntimeError("""Existing "bin" sub-directory safety-guard:
  A "bin" sub-directory exists already in the current working directory,
  but it does not contain any "libtbx." commands. Therefore the current
  working directory does not appear to be an existing libtbx-managed
  build directory. To resolve this problem:
    - If this command was accidentally run in the wrong directory,
      change to the correct directory.
    - Remove the "bin" subdirectory, then run this command again.""")
      #
      for path in buffer:
        remove_or_rename(path)

  def write_command_version_suffix(self):
    assert self.command_version_suffix is not None
    path = self.under_build("command_version_suffix")
    try: f = open(path, "w")
    except IOError:
      raise RuntimeError(
        'Cannot write command_version_suffix file: "%s"' % path)
    print >> f, self.command_version_suffix

  def read_command_version_suffix(self):
    path = self.build_path / "command_version_suffix"
    if not path.isfile():
        self.command_version_suffix = None
    else:
      try:
        self.command_version_suffix = path.open().read().strip()
      except IOError:
        raise RuntimeError(
          'Cannot read command_version_suffix file: "%s"' % path)

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

  def raise_not_found_in_repositories(self, message):
    if (isinstance(message, str)): message = [message]
    else: message = list(message)
    message.append("  Repository directories searched:")
    for path in self.repository_paths:
      message.append("    %s" % show_string(abs(path)))
    raise RuntimeError("\n".join(message))

  def listdir_in_repositories(self, test=None):
    for path in self.repository_paths:
      for name in path.listdir():
        if (test is None or test(abs(path / name))):
          yield path, name

  def match_in_repositories(self,
        relative_path_pattern,
        test=op.isdir,
        optional=True,
        must_be_unique=True):
    compiled_pattern = re.compile(relative_path_pattern)
    all_matches = []
    for path,name in self.listdir_in_repositories(test=test):
      if (compiled_pattern.match(name) is not None):
        if (len(all_matches) > 0 and all_matches[-1][0] != path):
          break
        all_matches.append((path, name))
    if (len(all_matches) == 0):
      if (not optional):
        self.raise_not_found_in_repositories(
          message="Cannot locate: %s" % show_string(relative_path_pattern))
      return None
    all_matches.sort()
    if (len(all_matches) != 1):
      if (must_be_unique):
        message = ["Multiple matches for search pattern %s:" %
          show_string(relative_path_pattern)]
        message.append("  Repository directory: %s" %
          show_string(all_matches[0][0]))
        for path,name in all_matches:
          message.append('    %s' % show_string(name))
        raise RuntimeError("\n".join(message))
    return libtbx.path.norm_join(*all_matches[0])

  def find_in_repositories(self,
        relative_path,
        test=op.isdir,
        optional=True,
        return_relocatable_path=False):
    assert len(relative_path) != 0
    for path in self.repository_paths:
      result = path / relative_path
      if test is None or test(abs(result)):
        if return_relocatable_path:
          return result
        else:
          return abs(result)
    if (not optional):
      self.raise_not_found_in_repositories(
        message="Cannot locate: %s" % show_string(relative_path))
    return None

  def find_dist_path(self, module_name, optional=False,
                     return_relocatable_path=False):
    dist_path = self.command_line_redirections.get(module_name, None)
    if (dist_path is not None):
      return dist_path.self_or_abs_if(return_relocatable_path)
    dist_path = self.find_in_repositories(relative_path=module_name,
                                          return_relocatable_path=True)
    if (dist_path is not None):
      return dist_path.self_or_abs_if(return_relocatable_path)
    trial_module = module(env=self, name=module_name)
    mate_name = trial_module.names[1]
    if (mate_name != module_name):
      if (self.find_in_repositories(relative_path=mate_name,
                                    return_relocatable_path=True) is not None):
        dist_path = self.match_in_repositories(
          relative_path_pattern="%s(?!_%s)" % (
            module_name, trial_module.mate_suffix))
        if (dist_path is not None):
          return dist_path.self_or_abs_if(return_relocatable_path)
    if (not optional):
      self.raise_not_found_in_repositories(
        message="Module not found: %s" % module_name)
    return None

  def process_module(self, dependent_module, module_name, optional):
    dist_path = self.find_dist_path(module_name, optional=optional)
    if (dist_path is None): return False
    new_module = module(env=self, name=module_name, dist_path=dist_path)
    new_name_normcase = op.normcase(new_module.name)
    for name in self.module_dict:
      if (op.normcase(name) == new_name_normcase): return True
    new_module.find_mate()
    new_module.process_libtbx_config()
    self.register_module(dependent_module=dependent_module, module=new_module)
    new_module.process_dependencies()
    return True

  def add_repository(self, path):
    path = relocatable_path(self, path)
    if path not in self.repository_paths:
      self.repository_paths.append(path)

  def process_args(self, pre_processed_args):
    command_line = pre_processed_args.command_line
    for path in pre_processed_args.repository_paths:
      if not op.isabs(path):
        path = abs(self.build_path / path)
      self.add_repository(path=path)
    module_names = []
    for module_name in command_line.args:
      if (len(module_name) == 0): continue # ignore arguments like ""
      if (module_name == ".."):
        raise RuntimeError('Invalid module name: ".."')
      if (module_name == "."): module_name = "libtbx"
      if (module_name in self.command_line_redirections):
        del self.command_line_redirections[module_name]
      elif (module_name.count("=") == 1
            and self.find_dist_path(module_name, optional=True) is None):
        module_name, redirection = module_name.split("=")
        dist_path = relocatable_path(self, op.expandvars(redirection))
        if not dist_path.isdir():
          raise RuntimeError(
            'Invalid command line redirection:\n'
            '  module name = "%s"\n'
            '  redirection = "%s"\n'
            '  resulting target = "%s"' % (
              module_name, redirection, dist_path))
        self.command_line_redirections[module_name] = dist_path
      module_names.append(module_name)
    if (pre_processed_args.warm_start):
      if (not command_line.options.only):
        for module in self.module_list:
          module_names.append(module.name)
    else:
      self.build_options = build_options(
        compiler=command_line.options.compiler,
        mode=command_line.options.build,
        warning_level=command_line.options.warning_level,
        precompile_headers=command_line.options.precompile_headers,
        static_libraries=command_line.options.static_libraries,
        static_exe=command_line.options.static_exe,
        scan_boost=command_line.options.scan_boost,
        write_full_flex_fwd_h=command_line.options.write_full_flex_fwd_h,
        boost_python_no_py_signatures
          =command_line.options.boost_python_no_py_signatures,
        boost_python_bool_int_strict
          =command_line.options.boost_python_bool_int_strict,
        enable_boost_threads=command_line.options.enable_boost_threads,
        enable_openmp_if_possible
          =command_line.options.enable_openmp_if_possible,
        enable_cuda=command_line.options.enable_cuda,
        opt_resources=command_line.options.opt_resources,
        use_environment_flags=command_line.options.use_environment_flags,
        force_32bit=command_line.options.force_32bit,
        msvc_arch_flag=command_line.options.msvc_arch_flag)
      self.build_options.get_flags_from_environment()
      if (command_line.options.command_version_suffix is not None):
        self.command_version_suffix = \
          command_line.options.command_version_suffix
        self.write_command_version_suffix()
    if (command_line.options.build_boost_python_extensions is not None):
      self.build_options.build_boost_python_extensions \
        = command_line.options.build_boost_python_extensions
    self.reset_module_registry()
    module_names.append("libtbx")
    module_names.reverse()
    for module_name in module_names:
      self.process_module(
        dependent_module=None, module_name=module_name, optional=False)
    self.scons_dist_path = self.find_dist_path("scons", optional=True)
    self.path_utility = self.under_dist(
      "libtbx", "command_line/path_utility.py",
      return_relocatable_path=True)
    assert self.path_utility.isfile()

  def dispatcher_precall_commands(self):
    if (self._dispatcher_precall_commands is None):
      lines = []
      if (    self.python_version_major_minor == (2,2)
          and sys.platform.startswith("linux")
          and op.isfile("/etc/redhat-release")):
        try: red_hat_linux_release = open("/etc/redhat-release").readline()
        except KeyboardInterrupt: raise
        except Exception: pass
        else:
          if (    red_hat_linux_release.startswith("Red Hat Linux release")
              and red_hat_linux_release.split()[4] == "9"):
            lines.extend([
              'if [ ! -n "$LD_ASSUME_KERNEL" ]; then',
              '  LD_ASSUME_KERNEL=2.4.1',
              '  export LD_ASSUME_KERNEL',
              'fi'])
      self._dispatcher_precall_commands = lines
    return self._dispatcher_precall_commands

  def write_dispatcher_include_template(self):
    if (os.name == "nt"): return
    print "    dispatcher_include_template.sh"
    f = self.under_build("dispatcher_include_template.sh",
                         return_relocatable_path=True).open("w")
    print >> f, "# include at start"
    print >> f, "#   Commands to be executed at the start of the"
    print >> f, "#   auto-generated dispatchers in bin."
    print >> f, "#"
    print >> f, "# include before command"
    print >> f, "#   Commands to be executed before the target command"
    print >> f, "#   is called by the auto-generated dispatchers in bin."
    print >> f, "#"
    print >> f, "# To see how the dispatchers work, look at an example:"
    print >> f, "#   %s" % show_string(self.under_build("bin/libtbx.help"))
    print >> f, "#"
    f.close()

  def reset_dispatcher_bookkeeping(self):
    self._dispatcher_registry = {}
    self._dispatcher_include_at_start = []
    self._dispatcher_include_before_command = []
    include_files = []
    for file_name in self.build_path.listdir():
      path = self.under_build(file_name, return_relocatable_path=True)
      if not path.isfile(): continue
      if (    file_name.startswith("dispatcher_include")
          and file_name.endswith(".sh")
          and file_name != "dispatcher_include_template.sh"):
        include_files.append(path)
    include_files.sort()
    for path in include_files:
      print "Processing: %s" % show_string(abs(path))
      lines = path.open().read().splitlines()
      lines_at_start = []
      lines_before_command = []
      buffer = lines_before_command
      for line in lines:
        l = " ".join(line.split()).lower()
        if (l.startswith("#include ")):
          l = "# " + l[1:]
        if   (l == "# include at start"):
          buffer = lines_at_start
        elif (l == "# include before command"):
          buffer = lines_before_command
        else:
          buffer.append(line)
      for buffer,target in [(lines_at_start,
                             self._dispatcher_include_at_start),
                            (lines_before_command,
                             self._dispatcher_include_before_command)]:
        if (len(buffer) == 0): continue
        buffer.insert(0, "# included from %s" % abs(path))
        highlight_dispatcher_include_lines(buffer)
        target.extend(buffer)

  def dispatcher_include(self, where):
    assert where in ["at_start", "before_command"]
    assert hasattr(self, "_dispatcher_include_at_start")
    if (where == "at_start"):
      return self._dispatcher_include_at_start
    return self._dispatcher_include_before_command

  def opt_resources_ld_preload(self):
    if (self.build_options.compiler in ["icc", "icpc"]):
      raise RuntimeError(
        "--opt-resources not supported in combination with --compiler=icc")
    def get_libs_dir():
      if (sys.platform.startswith("linux")) :
        libs = ["libimf.so", "libirc.so"]
        if (is_64bit_architecture()):
          return "linux64", libs
        return "linux32", libs
      if (sys.platform == "darwin"):
        libs = ["libimf.dylib", "libirc.dylib"]
        if (is_64bit_architecture()):
          return "darwin64", libs
        if (    platform.mac_ver()[0].startswith("10.6")
            and platform.machine() != "Power Macintosh"):
          return None, None # some extensions hang
        return "darwin32", libs
      return None, None
    libs_dir, libs = get_libs_dir()
    if (libs_dir is None):
      return None
    d = self.find_in_repositories(
      relative_path=op.join("opt_resources", libs_dir),
      optional=False)
    result = []
    for l in libs:
      p = d / l
      if not p.isfile():
        raise RuntimeError(
          "Missing file: %s" % show_string(abs(p)))
      if (p.find(":") >= 0):
        raise RuntimeError(
          "File name with embedded colon not supported: %s"
          % show_string(abs(p)))
      result.append(p)
    return ":".join(result)

  def write_bin_sh_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):
    f = target_file.open("w")
    if (source_file is not None):
      print >> f, '#! /bin/sh'
      print >> f, '# LIBTBX_DISPATCHER DO NOT EDIT'
    else:
      print >> f, '# LIBTBX_DISPATCHER_HEAD DO NOT EDIT'
      print >> f, '#'
      print >> f, '# This file is intended to be sourced from other scripts.'
      print >> f, '# It is like the dispatcher scripts in the bin directory,'
      print >> f, '# but only sets up the environment without calling a'
      print >> f, '# command at the end.'
    print >> f, '#'
    write_do_not_edit(f=f)
    print >> f, '# To customize this auto-generated script create'
    print >> f, '#'
    print >> f, '#   dispatcher_include*.sh'
    print >> f, '#'
    print >> f, '# files in %s and run' % show_string(abs(self.build_path))
    print >> f, '#'
    print >> f, '#   libtbx.refresh'
    print >> f, '#'
    print >> f, '# to re-generate the dispatchers (libtbx.refresh is a subset'
    print >> f, '# of the functionality of the libtbx/configure.py command).'
    print >> f, '#'
    print >> f, '# See also:'
    print >> f, '#   %s' \
      % show_string(self.under_build("dispatcher_include_template.sh"))
    print >> f, '#'
    print >> f, 'unset PYTHONHOME'
    print >> f, 'LC_ALL=C'
    print >> f, 'export LC_ALL'
    print >> f, 'LIBTBX_BUILD="$(cd "$(dirname "$0")" && cd .. && pwd -P)"'
    print >> f, 'export LIBTBX_BUILD'
    print >> f, 'LIBTBX_ROOT="$(dirname "$LIBTBX_BUILD")"'
    print >> f, 'export LIBTBX_ROOT'
    print >> f, 'LIBTBX_PYEXE_BASENAME="%s"' % self.python_exe.basename()
    print >> f, 'export LIBTBX_PYEXE_BASENAME'
    source_is_py = False
    if (source_file is not None):
      dispatcher_name = target_file.basename()
      if (dispatcher_name.find('"') >= 0):
        raise RuntimeError(
          "Dispatcher target file name contains double-quote: %s\n"
            % dispatcher_name
          + "  source file: %s" % source_file)
      print >> f, 'LIBTBX_DISPATCHER_NAME="%s"' % target_file.basename()
      print >> f, 'export LIBTBX_DISPATCHER_NAME'
      if source_file.ext().lower() == ".py":
        source_is_py = True
    for line in self.dispatcher_include(where="at_start"):
      print >> f, line
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((ld_library_path_var_name(), [self.lib_path]))
    essentials.append(("PATH", [self.bin_path]))
    for n,v in essentials:
      if (len(v) == 0): continue
      v = ":".join([ op.join('$LIBTBX_ROOT', p.relocatable) for p in v ])
      print >> f, 'if [ -n "$%s" ]; then' % n
      print >> f, '  %s="%s:$%s"' % (n, v, n)
      print >> f, '  export %s' % n
      print >> f, 'else'
      print >> f, '  %s="%s"' % (n, v)
      print >> f, '  export %s' % n
      print >> f, 'fi'
    precall_commands = self.dispatcher_precall_commands()
    if (precall_commands is not None):
      for line in precall_commands:
        print >> f, line
    if (source_is_py):
      scan_for_dispatcher_includes = True
    elif source_file is None or not source_file.isfile():
      scan_for_dispatcher_includes = False
    else:
      scan_for_dispatcher_includes = not detect_binary_file.from_initial_block(
        file_name=abs(source_file))
    if (scan_for_dispatcher_includes):
      for line in source_specific_dispatcher_include(
                    pattern="LIBTBX_PRE_DISPATCHER_INCLUDE_SH",
                    source_file=source_file):
        print >> f, line
    for line in self.dispatcher_include(where="before_command"):
      print >> f, line
    if (scan_for_dispatcher_includes):
      for line in source_specific_dispatcher_include(
                    pattern="LIBTBX_POST_DISPATCHER_INCLUDE_SH",
                    source_file=source_file):
        print >> f, line
    if (self.build_options.opt_resources):
      ldpl = self.opt_resources_ld_preload()
      if (ldpl is not None):
        print >> f, 'if [ "${LIBTBX_NO_LD_PRELOAD-UNSET}" == UNSET ]; then'
        print >> f, '  LD_PRELOAD="%s"' % ldpl
        print >> f, '  export LD_PRELOAD'
        print >> f, 'fi'
    print >> f, 'LIBTBX_PYEXE="%s"' % (
      self.python_exe.dirname() / "$LIBTBX_PYEXE_BASENAME").sh_value()
    print >> f, 'export LIBTBX_PYEXE'
    if (source_file is not None):
      def pre_cmd():
        return ['', '/usr/bin/arch -i386 '][self.build_options.force_32bit]
      cmd = ""
      if (source_is_py or source_is_python_exe):
        cmd += ' %s"$LIBTBX_PYEXE"%s' % (pre_cmd(), qnew)
      start_python = False
      if (source_is_py):
        if (len(source_specific_dispatcher_include(
                  pattern="LIBTBX_START_PYTHON",
                  source_file=source_file)) > 3):
          start_python = True
      if (not start_python and not source_is_python_exe):
        cmd += ' "%s"' % source_file.sh_value()
      print >> f, 'if [ -n "$LIBTBX__VALGRIND_FLAG__" ]; then'
      print >> f, "  exec $LIBTBX_VALGRIND"+cmd, '"$@"'
      print >> f, "elif [ $# -eq 0 ]; then"
      print >> f, "  exec"+cmd
      print >> f, "else"
      print >> f, "  exec"+cmd, '"$@"'
      print >> f, "fi"
    f.close()
    target_file.chmod(0755)

  def write_win32_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):
    f = target_file.open('w')
    # By default, changes to environment variables are  permanent on Windows,
    # i.e. it is as if export VAR was added after each set VAR=...
    # As a result, e.g. set PYTHONPATH=...; %PYTHONPATH% results in growing
    # PYTHONPATH each time a dispatcher script is run.
    # Thus setlocal essential (endlocal is implied)
    print >>f, '@setlocal'
    print >>f, '@set LIBTBX_BUILD=%~dp0'
    print >>f, '@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%'
    print >>f, r'@for %%F in ("%LIBTBX_BUILD%") do @set LIBTBX_BUILD=%%~dpF'
    print >>f, '@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%'
    print >>f, r'@for %%F in ("%LIBTBX_BUILD%") do @set LIBTBX_ROOT=%%~dpF'
    print >>f, '@set LIBTBX_ROOT=%LIBTBX_ROOT:~0,-1%'
    print >>f, '@set LIBTBX_DISPATCHER_NAME=%~nx0'
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((ld_library_path_var_name(), [self.lib_path]))
    essentials.append(("PATH", [self.bin_path]))
    for n,v in essentials:
      if (len(v) == 0): continue
      v = ';'.join([ op.join('%LIBTBX_ROOT%', p.relocatable) for p in v ])
      print >>f, '@set %s=%s;%%%s%%' % (n, v, n)
    print >>f, '@set LIBTBX_PYEXE=%s' % self.python_exe.bat_value()
    if source_file.ext().lower() == '.py':
      print >>f, '@"%%LIBTBX_PYEXE%%"%s "%s" %%*' % (
        qnew, source_file.bat_value())
    elif source_file.basename().lower() == 'python.exe':
      print >>f, '@"%%LIBTBX_PYEXE%%"%s %%*' % qnew
    else:
      print >>f, '@"%s" %%*' % source_file.bat_value()
    f.close()

  def write_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):
    source_file = self.as_relocatable_path(source_file)
    target_file = self.as_relocatable_path(target_file)
    reg = self._dispatcher_registry.setdefault(target_file, source_file)
    if (reg != source_file):
      if not reg.isfile():
        self._dispatcher_registry[target_file] = source_file
      elif (source_file.isfile()
            and (   not hasattr(os.path, "samefile")
                 or not reg.samefile(source_file))
            and  reg != source_file):
        raise RuntimeError("Multiple sources for dispatcher:\n"
          + "  target file:\n"
          + "    %s\n" % show_string(abs(target_file))
          + "  source files:\n"
          + "    %s\n" % show_string(abs(reg))
          + "    %s" % show_string(abs(source_file)))
    if (os.name == "nt"):
      action = self.write_win32_dispatcher
    else:
      action = self.write_bin_sh_dispatcher
    target_file_ext = target_file
    if os.name == 'nt':
      target_file_ext += '.bat'
    remove_or_rename(target_file_ext)
    try: action(source_file, target_file_ext, source_is_python_exe)
    except IOError, e: print "  Ignored:", e

  def _write_dispatcher_in_bin(self,
        source_file, target_file, source_is_python_exe=False):
    self.write_dispatcher(
      source_file=source_file,
      target_file=self.under_build("bin/"+target_file,
                                   return_relocatable_path=True),
      source_is_python_exe=source_is_python_exe)

  def write_dispatcher_in_bin(self, source_file, target_file):
    self._write_dispatcher_in_bin(
      source_file=source_file,
      target_file=target_file)

  def write_lib_dispatcher_head(self, target_file="dispatcher_head.sh"):
    if (os.name == "nt"): return
    print "   ", target_file
    self.write_bin_sh_dispatcher(
      source_file=None,
      target_file=self.under_build(target_file, return_relocatable_path=True))

  def write_setpaths_sh(self, suffix):
    setpaths = unix_setpaths(self, "sh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      write_do_not_edit(f=f)
      f.write("""\
ocwd="`pwd`"
cd "%s"
LIBTBX_BUILD=`pwd -P`
export LIBTBX_BUILD
cd ..
LIBTBX_ROOT=`pwd -P`
export LIBTBX_ROOT
LIBTBX_OPATH="$PATH"
export LIBTBX_OPATH
PATH="$LIBTBX_BUILD/bin:$PATH"
export PATH
cd "$ocwd"
ocwd=
""" % abs(self.build_path))
    s.write("""\
alias libtbx.setpaths_all=". '$LIBTBX_BUILD/setpaths_all.sh'"
alias libtbx.unsetpaths=". '$LIBTBX_BUILD/unsetpaths.sh'"
""")
    print >> u, 'unalias libtbx.unsetpaths > /dev/null 2>&1'
    if (self.is_development_environment()):
      print >> s, '''alias cdlibtbxbuild="cd '$LIBTBX_BUILD'"'''
      print >> u, 'unalias cdlibtbxbuild > /dev/null 2>&1'
    setpaths.all_and_debug()
    setpaths.set_unset_vars()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.sh_value(),
      var_name_in="LIBTBX_OPATH")
    for f in s, u:
      print >> f, 'LIBTBX_OPATH='
      if (suffix != "_all"):
        print >> f, 'LIBTBX_ROOT='
      if (suffix == ""):
        print >> f, 'LIBTBX_BUILD='

  def write_setpaths_csh(self, suffix):
    setpaths = unix_setpaths(self, "csh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, "#! /bin/tcsh"
      write_do_not_edit(f=f)
      f.write("""\
set ocwd="$cwd"
cd "%s"
setenv LIBTBX_BUILD "`/bin/sh -c 'pwd -P'`"
cd ..
setenv LIBTBX_ROOT "`/bin/sh -c 'pwd -P'`"
setenv LIBTBX_OPATH "$PATH"
setenv PATH "$LIBTBX_BUILD/bin:$PATH"
cd "$ocwd"
unset ocwd
""" % abs(self.build_path))
    s.write("""\
alias libtbx.setpaths_all "source '$LIBTBX_BUILD/setpaths_all.csh'"
alias libtbx.unsetpaths "source '$LIBTBX_BUILD/unsetpaths.csh'"
""")
    print >> u, 'unalias libtbx.unsetpaths'
    if (self.is_development_environment()):
      print >> s, '''alias cdlibtbxbuild "cd '$LIBTBX_BUILD'"'''
      print >> u, 'unalias cdlibtbxbuild'
    setpaths.all_and_debug()
    setpaths.set_unset_vars()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.sh_value(),
      var_name_in="LIBTBX_OPATH")
    for f in s, u:
      print >> f, 'unsetenv LIBTBX_OPATH'
      if (suffix != "_all"):
        print >> f, 'unsetenv LIBTBX_ROOT'
      if (suffix == ""):
        print >> f, 'unsetenv LIBTBX_BUILD'

  def write_setpaths_bat(self, suffix):
    setpaths = windows_setpaths(self, suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      write_do_not_edit(f=f, win_bat=True)
      print >> f, r'''@set LIBTBX_BUILD=%~dp0
@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%
@for %%F in ("%LIBTBX_BUILD%") do @set LIBTBX_ROOT=%%~dpF
@set LIBTBX_ROOT=%LIBTBX_ROOT:~0,-1%
@set LIBTBX_OPATH=%PATH%'''
      print >> f, '@set PATH=%s;%%PATH%%' % self.bin_path.bat_value()
    setpaths.all_and_debug()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.bat_value(),
      var_name_in="LIBTBX_OPATH")
    for command in ["setpaths_all", "unsetpaths"]:
      print >> s, '@doskey libtbx.%s="%s\\%s.bat"' % (
        command, "%LIBTBX_BUILD%", command)
    print >> u, '@doskey libtbx.unsetpaths='
    if (self.is_development_environment()):
      print >> s, '@doskey cdlibtbxbuild=cd "%LIBTBX_BUILD%"'
      print >> u, '@doskey cdlibtbxbuild='
    if (suffix == "_debug"):
      print >> s, '@set PYTHONCASEOK=1' # no unset
    setpaths.set_unset_vars()
    for f in s, u:
      print >> f, '@set LIBTBX_OPATH='
      if (suffix != "_all"):
        print >> f, '@set LIBTBX_ROOT='
      if (suffix == ""):
        print >> f, '@set LIBTBX_BUILD='

  def write_SConstruct(self):
    f = open_info(self.under_build("SConstruct",
                                   return_relocatable_path=True))
    write_do_not_edit(f=f)
    print >> f, 'SConsignFile()'
    for path in self.repository_paths:
      print >> f, 'Repository(r"%s")' % abs(path)
    for module in self.module_list:
      name,path  = list(module.name_and_dist_path_pairs())[-1]
      for script_name in ["libtbx_SConscript", "SConscript"]:
        if (path / script_name).isfile():
          print >> f, 'SConscript("%s/%s")' % (name, script_name)
          break
    f.close()

  def write_Makefile(self):
    if (op.isfile("Makefile")):
      os.rename("Makefile", "Makefile.old")
        # make cja seems to get confused if the file is simply overwritten
    f = open_info(self.under_build("Makefile",
                                   return_relocatable_path=True))
    lsj = './bin/libtbx.scons -j "`./bin/libtbx.show_number_of_processors`"'
    f.write("""\
# DO NOT EDIT THIS FILE!
# This file will be overwritten by the next libtbx/configure.py,
# libtbx.configure, or libtbx.refresh.

default:
\t%(lsj)s

nostop:
\t%(lsj)s -k

bp:
\t%(lsj)s -k boost_python_tests=1

reconf:
\t./bin/libtbx.configure .
\t%(lsj)s

redo:
\t./bin/libtbx.configure . --clear-scons-memory
\t%(lsj)s

# example
selfx:
\trm -rf selfx_tmp
\tmkdir selfx_tmp
\tcd selfx_tmp ; \\
\t\tlibtbx.start_binary_bundle example boost ; \\
\t\tlibtbx.bundle_as_selfx example build_id ; \\
\t\tmv example_build_id.selfx .. ; \\
\t\tcd .. ; \\
\t\tls -l example_build_id.selfx
""" % vars())
    f.close()

  def collect_test_scripts(self):
    result = []
    for module in self.module_list:
      result.extend(module.collect_test_scripts())
    return result

  def write_run_tests_csh(self):
    test_scripts = self.collect_test_scripts()
    if (len(test_scripts) > 0):
      path = self.under_build("run_tests.csh", return_relocatable_path=True)
      f = open_info(path)
      print >> f, "#! /bin/csh -f"
      print >> f, "set noglob"
      print >> f, "set verbose"
      for file_name in test_scripts:
        print >> f, 'libtbx.python "%s" $*' % abs(file_name)
      f.close()
      path.chmod(0755)

  def pickle(self):
    self.reset_dispatcher_support()
    file_name = self.build_path / "libtbx_env"
    pickle.dump(self, file_name.open("wb"), 0)

  def show_module_listing(self):
    print "Rooted at: %s" % self.root_path
    print "Top-down list of all modules involved:"
    top_down_module_list = list(self.module_list)
    top_down_module_list.reverse()
    labels = [module.names_for_module_listing()
      for module in top_down_module_list]
    if (len(labels) == 0): return
    fmt = "  %%-%ds  %%s" % max([len(label) for label in labels])
    for label,module in zip(labels,top_down_module_list):
      for dist_path in module.dist_paths_active():
        print fmt % (label, show_string(abs(dist_path)))
        label = ""

  def show_build_options_and_module_listing(self):
    print "Python: %s %s" % (
      sys.version.split()[0], show_string(sys.executable))
    if (self.is_ready_for_build()):
      self.build_options.report()
    print "command_version_suffix:", self.command_version_suffix
    self.show_module_listing()
    if (len(self.missing_for_use) > 0):
      raise RuntimeError("Missing modules: "
        + " ".join(sorted(self.missing_for_use)))
    if (not self.is_ready_for_build()):
      if (self.scons_dist_path is not None):
        print "***********************************"
        print "Warning: modules missing for build:"
        for module_name in sorted(self.missing_for_build):
          print " ", module_name
        print "***********************************"
      remove_or_rename(self.under_build("SConstruct"))

  def write_setpath_files(self):
    for suffix in ["", "_all", "_debug"]:
      if (hasattr(os, "symlink")):
        self.write_setpaths_sh(suffix)
        self.write_setpaths_csh(suffix)
      else:
        self.write_setpaths_bat(suffix)

  def process_exe(self):
    for path in [self.exe_path,
                 self.under_build("exe_dev", return_relocatable_path=True)]:
      if path.isdir():
        print "Processing: %s" % show_string(abs(path))
        for file_name in path.listdir():
          if (file_name.startswith(".")): continue
          target_file = file_name
          if (os.name == "nt"):
            fnl = file_name.lower()
            if (fnl.endswith(".exe.manifest")): continue
            if (fnl.endswith(".exe")): target_file = file_name[:-4]
          self._write_dispatcher_in_bin(source_file=path / file_name,
                                        target_file=target_file)

  def write_python_and_show_path_duplicates(self):
    module_names = []
    for module in self.module_list:
      if (   len(module.command_line_directory_paths()) != 0
          or len(module.assemble_pythonpath()) != 0):
        module_names.append(module.name.lower())
    for module_name in module_names:
      self._write_dispatcher_in_bin(
        source_file=self.python_exe,
        target_file=module_name+".python",
        source_is_python_exe=True)
    d, b = self.python_exe.split()
    pythonw_exe = d / b.replace("python", "pythonw")
    if pythonw_exe.isfile():
      for module_name in module_names:
        self._write_dispatcher_in_bin(
          source_file=pythonw_exe,
          target_file=module_name+".pythonw",
          source_is_python_exe=True)
    def have_ipython():
      for file_name in self.bin_path.listdir():
        file_name_lower = file_name.lower()
        if (file_name_lower.startswith("libtbx.")):
          if (   file_name_lower == "libtbx.ipython"
              or file_name_lower.startswith("libtbx.ipython.")):
            return
    commands = ["show_build_path", "show_dist_paths"]
    if (have_ipython()): commands.append("ipython")
    for command in commands:
      source_file = self.under_dist(
        "libtbx", "command_line/"+command+".py")
      for module_name in module_names:
        self._write_dispatcher_in_bin(
          source_file=source_file,
          target_file=module_name+"."+command)

  def write_command_version_duplicates(self):
    if (self.command_version_suffix is None): return
    suffix = "_" + self.command_version_suffix
    for file_name in self.bin_path.listdir():
      if (file_name.startswith(".")): continue
      source_file = abs(self.bin_path / file_name)
      if (os.name == "nt" and file_name.lower().endswith(".exe")):
        target_file = source_file[:-4] + suffix + source_file[-4:]
      else:
        target_file = source_file + suffix
      remove_or_rename(target_file)
      try: shutil.copy(source_file, target_file)
      except IOError: pass

  def assemble_pythonpath(self):
    pythonpath = [self.lib_path]
    for module in self.module_list:
      pythonpath.extend(module.assemble_pythonpath())
    pythonpath.reverse()
    self.pythonpath = unique_paths(paths=pythonpath)

  def is_development_environment(self):
    libtbx_cvs_root = self.under_dist("libtbx", "CVS/Root",
                                      return_relocatable_path=True)
    if libtbx_cvs_root.isfile():
      try: libtbx_cvs_root = open(libtbx_cvs_root).read()
      except IOError: pass
      else:
        if (libtbx_cvs_root.lower().find("ccp") >= 0): return False
    for module in self.module_list:
      if (module.is_version_controlled()):
        return True
    return False

  def clear_scons_memory(self):
    (self.build_path / ".sconsign.dblite").remove()
    (self.build_path / ".sconf_temp").remove_tree()

  def refresh(self):
    completed_file_name = (self.build_path / "libtbx_refresh_is_completed")
    completed_file_name.remove()
    self.assemble_pythonpath()
    self.show_build_options_and_module_listing()
    self.reset_dispatcher_bookkeeping()
    print "Creating files in build directory: %s" \
      % show_string(abs(self.build_path))
    self.write_dispatcher_include_template()
    self.write_lib_dispatcher_head()
    self.write_setpath_files()
    self.pickle()
    os.environ["LIBTBX_BUILD"] = abs(self.build_path) # to support libtbx.load_env
    if (self.is_ready_for_build()):
      self.write_SConstruct()
      if (os.name != "nt"):
        self.write_Makefile()
    if (os.name != "nt"):
      self.write_run_tests_csh()
    self.clear_bin_directory()
    if not self.bin_path.isdir():
      self.bin_path.makedirs()
    python_dispatchers = ["libtbx.python"]
    if (self.is_development_environment()):
      python_dispatchers.append("python")
    for file_name in python_dispatchers:
      self._write_dispatcher_in_bin(
        source_file=self.python_exe,
        target_file=file_name,
        source_is_python_exe=True)
    for module in self.module_list:
      module.process_command_line_directories()
    for path in self.pythonpath:
      sys.path.insert(0, abs(path))
    for module in self.module_list:
      module.process_libtbx_refresh_py()
    self.write_python_and_show_path_duplicates()
    self.process_exe()
    self.write_command_version_duplicates()
    self.pickle()
    print >> completed_file_name.open("w"), "libtbx_refresh_is_completed"

  def get_module(self, name, must_exist=True):
    result = self.module_dict.get(name, None)
    if (result is None and must_exist):
      raise RuntimeError("libtbx.env.get_module(name=%s): unknown module" % (
        show_string(name)))
    return result

class module:

  def __init__(self, env, name, dist_path=None, mate_suffix="adaptbx"):
    self.env = env
    self.mate_suffix = mate_suffix
    mate_suffix = "_" + mate_suffix
    if (op.normcase(name).endswith(op.normcase(mate_suffix))):
      self.name = name[:-len(mate_suffix)]
      self.names = [self.name, name]
      if (dist_path is not None):
        self.dist_paths = [None, dist_path]
    else:
      self.name = name
      self.names = [name, name + mate_suffix]
      if (dist_path is not None):
        self.dist_paths = [dist_path, None]

  def names_active(self):
    for name,path in zip(self.names, self.dist_paths):
      if (path is not None): yield name

  def names_for_module_listing(self):
    names_active = list(self.names_active())
    if (len(names_active) == 1): return names_active[0]
    return "+".join([self.name, self.mate_suffix])

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
    self.exclude_from_binary_bundle = []
    dist_paths = []
    for dist_path in self.dist_paths:
      if (dist_path is not None):
        while True:
          path = dist_path / "libtbx_config"
          if not path.isfile():
            config = None
            break
          try: f = path.open()
          except IOError: raise RuntimeError(
            'Cannot open configuration file: "%s"' % abs(path))
          try: config = eval(" ".join(f.readlines()), {}, {})
          except KeyboardInterrupt: raise
          except Exception: raise RuntimeError(
            'Corrupt configuration file: "%s"' % abs(path))
          f.close()
          redirection = config.get("redirection", None)
          if (redirection is None):
            break
          if (not isinstance(redirection, str)):
            raise RuntimeError(
              'Corrupt configuration file:\n'
              '  file = "%s"\n'
              '  redirection must be a Python string' % path)
          new_dist_path = op.expandvars(redirection)
          if (not op.isabs(new_dist_path)):
            new_dist_path = libtbx.path.norm_join(dist_path, new_dist_path)
          if (not op.isdir(new_dist_path)):
            raise RuntimeError(
              'Invalid redirection:\n'
              '  file = "%s"\n'
              '  redirection = "%s"\n'
              '  resulting target = "%s"' % (path, redirection, new_dist_path))
          dist_path = new_dist_path
        if (config is not None):
          self.required_for_build.extend(config.get(
            "modules_required_for_build", []))
          self.required_for_use.extend(config.get(
            "modules_required_for_use", []))
          self.optional.extend(config.get(
            "optional_modules", []))
          sep = os.sep
          for re_pattern in config.get("exclude_from_binary_bundle", []):
            if (sep != "/"):
              re_pattern = re_pattern.replace("/", "\\"+sep)
            self.exclude_from_binary_bundle.append(re_pattern)
      dist_paths.append(dist_path)
    self.dist_paths = dist_paths

  def process_dependencies(self):
    for module_name in self.required_for_build:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_for_build.add(module_name)
    for module_name in self.required_for_use:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_for_use.add(module_name)
    for module_name in self.optional:
      if (not self.env.process_module(
           dependent_module=self, module_name=module_name, optional=True)):
        self.env.missing_optional.add(module_name)

  def assemble_pythonpath(self):
    result = []
    for dist_path in self.dist_paths_active():
      path = dist_path / "pythonpath"
      if path.isdir():
        result.append(path)
      for sub_dir in ["", self.name]:
        path = dist_path / sub_dir / "__init__.py"
        if path.isfile():
          result.append(path.dirname().dirname())
    return result

  def has_top_level_directory(self, directory_name):
    for dist_path in self.dist_paths_active():
      if (dist_path / directory_name).isdir():
        return True
    return False

  def is_version_controlled(self):
    for directory_name in ["CVS", ".svn", ".git"]:
      if (self.has_top_level_directory(directory_name=directory_name)):
        return True
    return False

  def write_dispatcher(self,
        source_dir,
        file_name,
        suppress_warning,
        target_file_name_infix="",
        scan_for_libtbx_set_dispatcher_name=False):
    assert target_file_name_infix == "" or not scan_for_libtbx_set_dispatcher_name
    source_dir = self.env.as_relocatable_path(source_dir)
    if (len(file_name) == 0): return
    source_file = source_dir / file_name
    if not source_file.isfile(): return
    file_name_lower = file_name.lower()
    if (file_name_lower.startswith("__init__.py")): return
    if (file_name_lower.endswith(".pyc")): return
    if (file_name_lower.endswith(".pyo")): return
    if (file_name.startswith(".")): return
    if (file_name.endswith("~")): return # ignore emacs backup files
    if (file_name == "ipython_shell_start.py" and self.name == "libtbx"):
      try: import IPython
      except ImportError: return
    ext = op.splitext(file_name_lower)[1]
    if (scan_for_libtbx_set_dispatcher_name):
      read_size = 1000
    else:
      read_size = 0
    check_for_hash_bang = False
    if (ext == ".launch"):
      assert read_size != 0
    elif (os.name == "nt"):
      if (ext not in windows_pathext): return
    elif (ext == ".bat"):
      return
    elif (ext not in [".sh", ".py"]):
      read_size = max(2, read_size)
      check_for_hash_bang = True
    target_files = []
    if (read_size != 0):
      try: source_text = source_file.open().read(read_size)
      except IOError:
        raise RuntimeError('Cannot read file: "%s"' % source_file)
      if (check_for_hash_bang and not source_text.startswith("#!")):
        if (not suppress_warning):
          msg = 'WARNING: Ignoring file "%s" due to missing "#!"' % (
            source_file)
          print "*"*len(msg)
          print msg
          print "*"*len(msg)
        return
      for line in source_text.splitlines():
        flds = line.split()
        try:
          i = flds.index("LIBTBX_SET_DISPATCHER_NAME")
        except ValueError:
          pass
        else:
          if (i+1 < len(flds)):
            target_files.append(flds[i+1])
        if (ext == ".launch" and "LIBTBX_LAUNCH_EXE" in flds):
          source_file = self.env.under_build(
            op.join(self.name, "exe", file_name[:-len(ext)]),
            return_relocatable_path=True)
          source_file.relocatable += exe_suffix
    if (len(target_files) == 0):
      target_file = self.name.lower() + target_file_name_infix
      if (not file_name_lower.startswith("main.")
           or file_name_lower.count(".") != 1):
        target_file += "." + op.splitext(file_name)[0]
      target_files.append(target_file)
    for target_file in target_files:
      self.env._write_dispatcher_in_bin(
        source_file=source_file,
        target_file=target_file)

  def command_line_directory_paths(self):
    result = []
    for dist_path in self.dist_paths_active():
      for sub_dir in ["command_line", self.name+"/command_line"]:
        path = dist_path / sub_dir
        if path.isdir():
          result.append(path)
    return result

  def process_command_line_directories(self):
    for source_dir in self.command_line_directory_paths():
      print "Processing: %s" % show_string(abs(source_dir))
      def is_py_sh(file_name):
        return file_name.endswith(".sh") \
            or file_name.endswith(".py")
      nodes = source_dir.listdir()
      py_sh_dict = {}
      for file_name in nodes:
        if (is_py_sh(file_name)):
          py_sh_dict.setdefault(file_name[:-3], []).append(file_name[-2:])
      for file_name in nodes:
        if (    is_py_sh(file_name)
            and len(py_sh_dict[file_name[:-3]]) == 2):
          if (os.name == "nt"): skip = ".sh"
          else:                 skip = ".py"
          if (file_name.endswith(skip)):
            continue
        self.write_dispatcher(
          source_dir=source_dir,
          file_name=file_name,
          suppress_warning=False,
          scan_for_libtbx_set_dispatcher_name=True)

  def process_python_command_line_scripts(self,
        source_dir,
        print_prefix="  ",
        target_file_name_infix="",
        scan_for_libtbx_set_dispatcher_name=False):
    source_dir = self.env.as_relocatable_path(source_dir)
    print print_prefix+"Processing: %s" % show_string(abs(source_dir))
    for file_name in source_dir.listdir():
      if (not file_name.endswith(".py")): continue
      self.write_dispatcher(
        source_dir=source_dir,
        file_name=file_name,
        suppress_warning=False,
        target_file_name_infix=target_file_name_infix,
        scan_for_libtbx_set_dispatcher_name
          =scan_for_libtbx_set_dispatcher_name)

  def process_libtbx_refresh_py(self):
    for dist_path in self.dist_paths_active():
      custom_refresh = dist_path / "libtbx_refresh.py"
      if custom_refresh.isfile():
        print "Processing: %s" % show_string(abs(custom_refresh))
        execfile(abs(custom_refresh), {}, {"self": self})

  def collect_test_scripts(self,
        file_names=["run_tests.py", "run_examples.py"]):
    result = []
    for dist_path in self.dist_paths_active():
      for file_name in file_names:
        path = dist_path / file_name
        if path.isfile(): result.append(path)
    return result

  def remove_obsolete_pyc_if_possible(self, pyc_file_names):
    for file_name in pyc_file_names:
      for dist_path in self.dist_paths_active():
        path = dist_path / file_name
        if path.isfile(): path.remove()

class build_options:

  supported_modes = [
    "release",
    "max_optimized",
    "quick",
    "debug",
    "debug_optimized",
    "profile"]

  def __init__(self,
        compiler,
        mode,
        warning_level,
        static_libraries,
        static_exe,
        scan_boost,
        write_full_flex_fwd_h=default_write_full_flex_fwd_h,
        build_boost_python_extensions=default_build_boost_python_extensions,
        boost_python_no_py_signatures=False,
        boost_python_bool_int_strict=True,
        enable_boost_threads=default_enable_boost_threads,
        enable_openmp_if_possible=default_enable_openmp_if_possible,
        enable_cuda=default_enable_cuda,
        opt_resources=default_opt_resources,
        precompile_headers=False,
        use_environment_flags=False,
        force_32bit=False,
        msvc_arch_flag=default_msvc_arch_flag):
    adopt_init_args(self, locals())
    assert self.mode in build_options.supported_modes
    assert self.warning_level >= 0
    self.optimization = (self.mode in [
      "release", "max_optimized", "debug_optimized", "profile"])
    self.max_optimized = (self.mode in [
      "max_optimized", "debug_optimized", "profile"])
    self.debug_symbols = (self.mode in [
      "debug", "debug_optimized", "profile"])
    if (self.static_exe):
      self.static_libraries = True
    if (self.msvc_arch_flag == "None"): self.msvc_arch_flag = None

  def get_flags_from_environment(self):
    if (self.use_environment_flags ):
      # get compiler flags from environment vars
      # they will be stored at configure in the file libtbx_env
      # and used during build
      self.env_cxxflags = ""
      self.env_cflags = ""
      self.env_cppflags = ""
      flg = os.environ.get("CXXFLAGS")
      if flg is not None:
        self.env_cxxflags = flg
      flg = os.environ.get("CFLAGS")
      if flg is not None:
        self.env_cflags = flg
      flg = os.environ.get("CPPFLAGS")
      if flg is not None:
        self.env_cppflags = flg

  def report(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Compiler:", self.compiler
    print >> f, "Build mode:", self.mode
    print >> f, "Warning level:", self.warning_level
    print >> f, "Precompiled Headers:", self.precompile_headers
    print >> f, "Static libraries:", self.static_libraries
    print >> f, "Static exe:", self.static_exe
    print >> f, "Scan Boost headers:", self.scan_boost
    print >> f, "Write full flex_fwd.h files:", self.write_full_flex_fwd_h
    print >> f, "Build Boost.Python extensions:", \
      self.build_boost_python_extensions
    print >> f, "Define BOOST_PYTHON_NO_PY_SIGNATURES:", \
      self.boost_python_no_py_signatures
    print >> f, "Define BOOST_PYTHON_BOOL_INT_STRICT:", \
      self.boost_python_bool_int_strict
    print >> f, "Boost threads enabled:", self.enable_boost_threads
    print >> f, "Enable OpenMP if possible:", self.enable_openmp_if_possible
    print >> f, "Enable CUDA:", self.enable_cuda
    print >> f, "Use opt_resources if available:", self.opt_resources
    print >> f, "Use environment flags:", self.use_environment_flags
    if( self.use_environment_flags ):
      print >>f, "  CXXFLAGS = ", self.env_cxxflags
      print >>f, "  CFLAGS = ", self.env_cflags
      print >>f, "  CPPFLAGS = ", self.env_cppflags

class include_registry:

  def __init__(self):
    self._boost_dir_name = "boost"
    self._had_message = {}
    self.scan_boost()

  def scan_boost(self, flag=False):
    self._scan_boost = flag
    return self

  def set_boost_dir_name(self, path):
    self._boost_dir_name = op.basename(path).lower()

  def scan_flag(self, path):
    if (not self._scan_boost
        and op.basename(path).lower() == self._boost_dir_name):
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
    paths = unique_paths(paths=paths)
    for path in paths:
      if (self.scan_flag(path)):
        env.Append(CPPPATH=[path])
      else:
        ipath = self.prepend_include_switch(env, path)
        env.Append(CXXFLAGS=[ipath])
        env.Append(SHCXXFLAGS=[ipath])

  def prepend(self, env, paths):
    assert isinstance(paths, list)
    paths = unique_paths(paths=paths)
    paths.reverse()
    for path in paths:
      if (self.scan_flag(path)):
        env.Prepend(CPPPATH=[path])
      else:
        ipath = self.prepend_include_switch(env, path)
        env.Prepend(CXXFLAGS=[ipath])
        env.Prepend(SHCXXFLAGS=[ipath])

class pre_process_args:

  def __init__(self, args, default_repositories=None):
    self.repository_paths = []
    if (len(args) == 0): args = ["--help"]
    if (default_repositories is None):
      command_name = "libtbx.configure"
      self.warm_start = True
    else:
      command_name = "libtbx/configure.py"
      self.warm_start = False
    from libtbx.option_parser import option_parser
    parser = option_parser(
      usage="%s [options] module_name[=redirection_path] ..." % command_name)
    if (self.warm_start):
      parser.option(None, "--only",
        action="store_true",
        default=False,
        help="disable previously configured modules")
    else:
      parser.option("-r", "--repository",
        action="callback",
        type="string",
        callback=self.option_repository,
        help="path to source code repository"
          " (may be specified multiple times;"
          " paths are searched in the order given)",
        metavar="DIRECTORY")
      if (hasattr(os.path, "samefile")):
        parser.option(None, "--current_working_directory",
          action="store",
          type="string",
          default=None,
          help="preferred spelling of current working directory"
            " (to resolve ambiguities due to soft links)",
          metavar="DIRECTORY")
      parser.option(None, "--build",
        choices=build_options.supported_modes,
        default="release",
        help="build mode (default: release)",
        metavar="|".join(build_options.supported_modes))
      parser.option(None, "--compiler",
        action="store",
        type="string",
        default="default",
        help="select non-standard compiler (platform dependent)",
        metavar="STRING")
      parser.option(None, "--warning_level",
        action="store",
        type="int",
        default=0,
        help="manipulate warning options (platform dependent)")
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
      parser.option(None, "--write_full_flex_fwd_h",
        action="store_true",
        default=default_write_full_flex_fwd_h,
        help="create full flex_fwd.h files to work around platform-specific"
             " problems (see comments in cctbx.source_generators.flex_fwd_h)")
      parser.option(None, "--command_version_suffix",
        action="store",
        type="string",
        default=None,
        help="version suffix for commands in bin directory",
        metavar="STRING")
      parser.option(None, "--use_environment_flags",
        action="store_true",
        default=False,
        help="add compiler flags from environment variables: CXXFLAGS, CFLAGS,"
             " CPPFLAGS")
      parser.option(None, "--force_32bit",
        action="store_true",
        default=False,
        help="Force 32-bit compilation on Mac OS 10.6 (Snow Leopard)\n"
             "Not compatible with /usr/bin/python: please run configure\n"
             "with /System/Library/Frameworks/Python.framework/"
             "Versions/2.x/bin/python")
      msvc_arch_flag_choices = ("None", "SSE", "SSE2")
      parser.option(None, "--msvc_arch_flag",
        choices=msvc_arch_flag_choices,
        default=default_msvc_arch_flag,
        help="choose MSVC CPU architecture instruction set"
             " for optimized builds",
        metavar="|".join(msvc_arch_flag_choices))
    parser.option(None, "--build_boost_python_extensions",
      action="store",
      type="bool",
      default=default_build_boost_python_extensions,
      help="build Boost.Python extension modules (default: %s)"
        % bool_literal(default_build_boost_python_extensions),
      metavar="True|False")
    parser.option(None, "--enable_boost_threads",
      action="store_true",
      default=False,
      help="enable threads in Boost")
    parser.option(None, "--enable_openmp_if_possible",
      action="store",
      type="bool",
      default=default_enable_openmp_if_possible,
      help="use OpenMP if available and known to work (default: %s)"
        % bool_literal(default_enable_openmp_if_possible),
      metavar="True|False")
    parser.option(None, "--enable_cuda",
      action="store_true",
      default=default_enable_cuda,
      help="Use optimized CUDA routines for certain calculations.  Requires at least one NVIDIA GPU with compute capability of 2.0 or higher, and CUDA Toolkit 4.0 or higher (default: don't)")
    parser.option(None, "--opt_resources",
      action="store",
      type="bool",
      default=default_opt_resources,
      help="use opt_resources if available (default: %s)"
        % bool_literal(default_opt_resources),
      metavar="True|False")
    parser.option(None, "--precompile_headers",
      action="store_true",
      default=False,
      help="Precompile headers, especially Boost Python ones (default: don't)")
    if (not self.warm_start):
      parser.option(None, "--boost_python_no_py_signatures",
        action="store_true",
        default=False,
        help="disable Boost.Python docstring Python signatures")
      parser.option(None, "--boost_python_bool_int_strict",
        action="store_false",
        default=True,
        help="disable Boost.Python implicit bool<->int conversions")
    parser.option(None, "--clear_scons_memory",
      action="store_true",
      default=False,
      help="remove scons build signatures and config cache")
    self.command_line = parser.process(args=args)
    if (len(self.command_line.args) == 0):
      raise RuntimeError(
        "At least one module name is required"
        " (use --help to obtain more information).")
    if (not hasattr(os.path, "samefile")):
      self.command_line.options.current_working_directory = None
    if (default_repositories is not None):
      self.repository_paths.extend(default_repositories)
    if (not self.warm_start):
      if (self.command_line.options.force_32bit):
        if (sys.platform != "darwin"):
          raise RuntimeError(
            "The --force_32bit option is only valid on Mac OS systems.")
        if (sys.maxint > 2**31-1):
          raise RuntimeError(
            'The --force_32bit option can only be used with 32-bit Python.\n'
            '  See also: "man python"')
        from libtbx import easy_run
        buffers = easy_run.fully_buffered(
          command="/usr/bin/arch -i386 /bin/ls /")
        if (   len(buffers.stderr_lines) != 0
            or len(buffers.stdout_lines) == 0):
          raise RuntimeError(
            "The --force_32bit option is not valid for this platform.")
      if (    self.command_line.options.msvc_arch_flag != "None"
          and os.name != "nt"):
        raise RuntimeError(
          "The --msvc_arch_flag option is not valid for this platform.")

  def option_repository(self, option, opt, value, parser):
    if (not op.isdir(value)):
      raise RuntimeError(
        "Not a directory: --repository %s" % show_string(value))
    self.repository_paths.append(value)

def set_preferred_sys_prefix_and_sys_executable(build_path):
  if (not hasattr(os.path, "samefile")): return
  dirname, basename = op.split(sys.prefix)
  if (op.samefile(build_path, dirname)):
    new_prefix = op.join(build_path, basename)
    if (op.samefile(new_prefix, sys.prefix)):
      p = op.normpath(op.normcase(sys.prefix))
      if (sys.prefix != new_prefix):
        sys.prefix = new_prefix
      e = op.normpath(op.normcase(sys.executable))
      if (e.startswith(p)):
        new_executable = sys.prefix + e[len(p):]
        if (op.samefile(new_executable, sys.executable)):
          if (sys.executable != new_executable):
            sys.executable = new_executable

def cold_start(args):
  cwd_was_empty_at_start = True
  for file_name in os.listdir("."):
    if (not file_name.startswith(".")):
      cwd_was_empty_at_start = False
      break
  default_repositories = []
  r = op.dirname(op.dirname(args[0]))
  b = op.basename(r)
  if (b.lower().startswith("cctbx_project")):
    default_repositories.append(op.dirname(r))
  default_repositories.append(r)
  pre_processed_args = pre_process_args(
    args=args[1:],
    default_repositories=default_repositories)
  build_path=pre_processed_args.command_line.options.current_working_directory
  if (build_path is None):
    build_path = os.getcwd()
  else:
    if (not op.isabs(build_path)):
      raise RuntimeError("Not an absolute path name:"
        " --current_working_directory %s" % show_string(build_path))
    if (not op.isdir(build_path)):
      raise RuntimeError("Not a directory:"
        " --current_working_directory %s" % show_string(build_path))
    if (not op.samefile(build_path, os.getcwd())):
      raise RuntimeError("Not equivalent to the current working directory:"
        " --current_working_directory %s" % show_string(build_path))
    n = len(os.sep)
    while (len(build_path) > n and build_path.endswith(os.sep)):
      build_path = build_path[:-n]
  set_preferred_sys_prefix_and_sys_executable(build_path=build_path)
  env = environment(build_path=build_path)
  env.process_args(pre_processed_args=pre_processed_args)
  if (   pre_processed_args.command_line.options.clear_scons_memory
      or cwd_was_empty_at_start):
    env.clear_scons_memory()
  env.refresh()

def unpickle():
  build_path = os.environ["LIBTBX_BUILD"]
  set_preferred_sys_prefix_and_sys_executable(build_path=build_path)
  libtbx_env = open(op.join(build_path, "libtbx_env"), "rb")
  env = pickle.load(libtbx_env)
  if (env.python_version_major_minor != sys.version_info[:2]):
    env.raise_python_version_incompatible()
  # XXX backward compatibility 2009-04-06
  if (not hasattr(env.build_options, "boost_python_bool_int_strict")):
    env.build_options.boost_python_bool_int_strict = True
  # XXX backward compatibility 2009-04-27
  if( not hasattr(env.build_options, "use_environment_flags") ):
    env.build_options.use_environment_flags = False
    env.build_options.env_cxxflags = ""
    env.build_options.env_cflags = ""
    env.build_options.env_cppflags = ""
  # XXX backward compatibility 2009-10-11
  if (not hasattr(env.build_options, "force_32bit")) :
    env.build_options.force_32bit = False
  # XXX backward compatibility 2009-10-13
  if (not hasattr(env.build_options, "msvc_arch_flag")) :
    env.build_options.msvc_arch_flag = default_msvc_arch_flag
  # XXX backward compatibility 2010-05-28
  if (not hasattr(env.build_options, "precompile_headers")):
    env.build_options.precompile_headers = False
  # XXX backward compatibility 2011-04-01
  if (not hasattr(env.build_options, "opt_resources")):
    env.build_options.opt_resources = False
  # XXX backward compatibility 2011-07-05
  if (not hasattr(env.build_options, "enable_cuda")):
    env.build_options.enable_cuda = False
  # XXX backward incompatibility 2011-10
  if not hasattr(env, 'relocatable'):
    print ("Please re-configure from scratch your cctbx_build:"
           "cd cctbx_build; "
           "/your/path/to/python ../cctbx_project/libtbx/configure.py [options] [modules]")
    sys.exit(1)
  env.set_build_path(build_path)
  return env

def warm_start(args):
  pre_processed_args = pre_process_args(args=args[1:])
  env = unpickle()
  env.process_args(pre_processed_args=pre_processed_args)
  if (pre_processed_args.command_line.options.clear_scons_memory):
    env.clear_scons_memory()
  env.refresh()

if (__name__ == "__main__"):
  if (len(sys.argv) == 2 and sys.argv[1] == "__libtbx_refresh__"):
    unpickle().refresh()
  else:
    warm_start(sys.argv)
  print "Done."
