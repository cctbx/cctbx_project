from __future__ import absolute_import, division, print_function
import libtbx.path
from libtbx.auto_build import regenerate_module_files
from libtbx.auto_build.installer_utils import call
from libtbx.path import relocatable_path, absolute_path
from libtbx.str_utils import show_string
from libtbx.utils import detect_binary_file, to_str
from libtbx import adopt_init_args
import platform
import shutil
from six.moves import zip, map
from six.moves import cPickle as pickle

import io
import os
import re
import site
import sys
import sysconfig
op = os.path

if os.environ.get('LIBTBX_WINGIDE_DEBUG'):
  import wingdbstub # special import

qnew = " -Qnew"

if (os.name == "nt"):
  exe_suffix = ".exe"
else:
  exe_suffix = ""

default_write_full_flex_fwd_h = sys.platform.startswith("irix")
default_msvc_arch_flag = ["None", "SSE2"][int(os.name == "nt")]
default_build_boost_python_extensions = True
default_enable_openmp_if_possible = False
default_enable_boost_threads = True
default_enable_cuda = False
default_enable_kokkos = False
default_opt_resources = False
default_enable_cxx11 = False
default_cxxstd = None
default_use_conda = False

def is_64bit_architecture():
  # this appears to be most compatible (hat tip: James Stroud)
  # http://stackoverflow.com/questions/1405913/how-do-i-determine-if-my-python-shell-is-executing-in-32bit-or-64bit-mode-on-os
  import struct
  nbits = 8 * struct.calcsize("P")
  return (nbits == 64)

def using_conda_python():
  '''
  Return True if Python is from conda, False otherwise.
  This check is independent of being in an active conda environment.
  https://stackoverflow.com/questions/47608532/how-to-detect-from-within-python-whether-packages-are-managed-with-conda?noredirect=1&lq=1
  '''
  conda_prefix = sys.prefix
  if (sys.platform == 'darwin'):
    if ('python.app' in conda_prefix):
      conda_prefix = conda_prefix.split('python.app')[0]
  return os.path.exists(os.path.join(conda_prefix, 'conda-meta'))

def get_conda_prefix():
  '''
  Return the root directory of the conda environment. This function will
  try to figure out the root directory if the environment is not active.
  A special case exists for macOS where the framework package (python.app)
  is used for GUI programs.

  A RuntimeError is raised if the root directory of the conda environment
  cannot be determined.
  '''
  conda_prefix = sys.prefix
  if (using_conda_python()):
    if (sys.platform == 'darwin'):  # case where python.app is used
      if ('python.app' in conda_prefix):
        conda_prefix = conda_prefix.split('python.app')[0]
  return conda_prefix

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
      "ld -dynamic%s -headerpad_max_install_names -r -bind_at_load -o %s $SOURCES" %
        (opt_m, lo),
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

def is_llvm_compiler(gcc='gcc'):
  from libtbx import easy_run
  try:
    gcc_version = (easy_run.fully_buffered(command='%s --version' % gcc)
                           .raise_if_errors()
                           .stdout_lines[0].strip())
    return "llvm" in gcc_version
  except Exception:
    return False

def get_gcc_version(command_name="gcc"):
  # run command in shell subprocess
  from libtbx import easy_run
  buffer = easy_run.fully_buffered(
    command="%s -dumpversion" % command_name)
  # if something went wrong `buffer.stdout_lines` might have len=0
  if len(buffer.stdout_lines) < 1:
      return None

  major_minor_patchlevel = buffer.stdout_lines[0].split(".")
  # output is not a valid version format:
  if (len(major_minor_patchlevel) not in [1,2,3]):
    return None
  # parse version number
  num = []
  for fld in major_minor_patchlevel:
    try: i = int(fld)
    except ValueError:
      return None
    num.append(i)
  # substitute missing minor and patchlevel for gcc7 on Ubuntu18
  if (len(num) == 1): num.append(0); num.append(0)
  if (len(num) == 2): num.append(0) # substitute missing patchlevel

  # format version numbers
  return ((num[0]*100)+num[1])*100+num[2]

def get_hostname():
  try:
    import socket
    return socket.gethostname()
  except Exception:
    return None

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

def python_include_path():
  if (sys.platform == "win32"):
    include_path = sys.prefix + r"\include"
  else:
    include_path = sys.prefix + "/include/python%d.%d" % sys.version_info[:2]
  include_path = libtbx.path.norm_join(include_path)

  # I believe the above code to be unnecessary.
  # If this flags up no problems then remove the code above and
  # the assertion test directly below in a month or so.
  if op.isdir(include_path):
    sysconfig_path = sysconfig.get_paths()['include']
    if sys.platform == "win32":
      include_path = include_path.lower()
      sysconfig_path = sysconfig_path.lower()
    assert include_path == sysconfig_path, \
        "%s != %s" % (include_path, sysconfig_path)

  include_path = sysconfig.get_paths()['include']
  if not op.isdir(include_path):
    try:  # conda environment
      conda_prefix = get_conda_prefix()
      include_path = os.path.join(conda_prefix, 'include',
                                  'python%d.%dm' % sys.version_info[:2])
      if not op.isdir(include_path):
        include_path = include_path[:-1]
    except RuntimeError:
      pass
  if not op.isdir(include_path):
    raise RuntimeError("Cannot locate Python's include directory: %s" % include_path)
  return include_path

def highlight_dispatcher_include_lines(lines):
  m = max([len(line) for line in lines])
  if (os.name == "nt"):
    lines.insert(0, "@REM " + "-"*(m-5))
  else :
    lines.insert(0, "# " + "-"*(m-2))
  lines.append(lines[0])

def source_specific_dispatcher_include(pattern, source_file):
  try:
    with io.open(abs(source_file), encoding='utf-8', errors='ignore') as fh:
      source_lines = to_str(fh.read()).splitlines()
  except IOError: return []
  if (os.name == "nt"):
    lines = ["@REM lines marked " + pattern]
  else :
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
  if (win_bat): s = "@rem"
  else:         s = "#"
  print(s+' THIS IS AN AUTOMATICALLY GENERATED FILE.', file=f)
  print(s+' DO NOT EDIT! CHANGES WILL BE LOST.', file=f)

def open_info(path, mode="w", info="   "):
  print(info, path.basename())
  try: return path.open(mode)
  except IOError as e:
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
      self.u = open(os.devnull, 'w')

  def all_and_debug(self):
    if (self.suffix == "_debug"):
      self.update_path("PYTHONPATH",
        os.pathsep.join([
          self.path_script_value(_) for _ in self.env.pythonpath]))
      self.update_path(
        self.env.ld_library_path_var_name(),
        os.pathsep.join([self.path_script_value(_)
          for _ in self.env.ld_library_path_additions()]))

  def set_unset_vars(self):
    if (self.suffix != ""):
      for var_name,path in self.env.var_name_and_build_or_dist_path_pairs():
        if (var_name == "LIBTBX_BUILD"): continue
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
      print('%s="%s"' % (var_name, val), file=self.s)
      print('export %s' % var_name, file=self.s)
      print('unset %s' % var_name, file=self.u)
    else:
      print('setenv %s "%s"' % (var_name, val), file=self.s)
      print('unsetenv %s' % var_name, file=self.u)

  def update_path(self, var_name, val, var_name_in=None):
    if (var_name_in is None): var_name_in = var_name
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      if (self.shell == "sh"):
        print('if [ -n "$%s" ]; then' % var_name_in, file=f)
        print('  LIBTBX_TMPVAL="$%s"' % var_name_in, file=f)
        print('else', file=f)
        print('  LIBTBX_TMPVAL=', file=f)
        print('fi', file=f)
        print('export LIBTBX_TMPVAL', file=f)
        fmt = \
          '''%s`libtbx.path_utility %s LIBTBX_TMPVAL "%s" < /dev/null`'''
      else:
        print('if ($?%s) then' % var_name_in, file=f)
        print('  setenv LIBTBX_TMPVAL "$%s"' % var_name_in, file=f)
        print('else', file=f)
        print('  unsetenv LIBTBX_TMPVAL', file=f)
        print('endif', file=f)
        fmt = \
          '''%s"`libtbx.path_utility %s LIBTBX_TMPVAL '%s' < /dev/null`"'''
      print(fmt % (self._setenv % var_name, action, val), file=f)
      if (self.shell == "sh"):
        if (f is self.s):
          print('export %s' % var_name, file=f)
        print('if [ "$%s" = "L_I_B_T_B_X_E_M_P_T_Y" ]; then unset %s; fi' % (
          var_name, var_name), file=f)
      else:
        print('if ("$%s" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv %s' % (
          var_name, var_name), file=f)

class windows_setpaths(common_setpaths):

  def __init__(self, env, suffix):
    common_setpaths.__init__(self, env, "bat", suffix)

  def path_script_value(self, path_obj):
    return path_obj.bat_value()

  def setenv(self, var_name, val):
    print('@set %s=%s' % (var_name, val), file=self.s)
    print('@set %s=' % var_name, file=self.u)

  def update_path(self, var_name, val, var_name_in=None):
    if (var_name_in is None): var_name_in = var_name
    fmt = '''\
@for /F "delims=" %%%%i in ('libtbx.path_utility %s %s "%s"') do @set %s=%%%%i'''
    for f,action in [(self.s, "prepend"), (self.u, "delete")]:
      print(fmt % (
        action,
        var_name_in,
        val,
        var_name), file=f)
      print('@if "%%%s%%" == "L_I_B_T_B_X_E_M_P_T_Y" @set %s=' % (
        var_name, var_name), file=f)

def _windows_pathext():
  result = os.environ.get("PATHEXT", "").lower().split(os.pathsep)
  for ext in [".py", ".exe", ".bat"]:
    if (ext not in result):
      result.insert(0, ext)
  return result

if os.name == "nt":
  windows_pathext = _windows_pathext()

def remove_or_rename(path):
  assert path is not None
  if not isinstance(path, str):
    path = abs(path)
  if op.isfile(path):
    try: os.remove(path)
    except OSError:
      try: os.remove(path+".old")
      except OSError: pass
      try: os.rename(path, path+".old")
      except OSError: pass

# Constant used in class environment.write_bin_sh_dispatcher -----------------

_SHELLREALPATH_CODE = '''
# ----------------------------------------------------------------------------
# The shellrealpath function resolves an absolute physical path of its
# first argument and stores it in a global shell variable RESULT.
# The function returns nonzero for unreadable or invalid symlinks
# and resets the RESULT to an empty string.

shellrealpath() {
    local ORGDIR="$PWD"
    local TARGET="$1"
    RESULT=""
    # This test fails for a symlink loop.  We can do without resolution
    # of symlinks that point to existing unreadable files.
    [ -r "$TARGET" ] || return $?
    # Check if the readlink command exists.
    type readlink >/dev/null || return $?
    while true; do
        cd "$(dirname "$TARGET")"
        TARGET="$(basename "$TARGET")"
        if [ -L "$TARGET" ]; then
            TARGET="$(readlink "$TARGET")"
            continue
        fi
        RESULT="$(pwd -P)/$TARGET"
        break
    done
    cd "$ORGDIR"
}
# ----------------------------------------------------------------------------
'''

# ----------------------------------------------------------------------------

class environment:

  def __init__(self, build_path):
    self.python_version_major_minor = sys.version_info[:2]
    self.build_path = absolute_path(build_path)
    self.manage_python_version_major_minor()
    self.reset_dispatcher_support()
    self.set_derived_paths()
    self.python_exe = self.as_relocatable_path(sys.executable)
    self.installed = False
    self.installed_modules = []  # used to track installed modules in local env
    self.installed_order = []
    # sanity checks
    assert self.python_exe.isfile()
    assert self.python_exe.access(os.X_OK)

    self.read_command_version_suffix()
    # the following are modules explicitly specified on the command-line
    # --only reset the list and --exclude prunes from it
    self.explicitly_requested_modules = set()
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
      + "  Build directory: %s\n" % show_string(abs(self.build_path))
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
    return True
    #return (len(self.missing_for_build) == 0)

  def as_relocatable_path(self, path):
    if isinstance(path, libtbx.path.path_mixin): return path
    return relocatable_path(self.build_path, path, resolve_symlinks=False)

  def set_derived_paths(self):
    r = self.as_relocatable_path
    self.bin_path     = r("bin")
    self.exe_path     = r("exe")
    self.lib_path     = r("lib")
    self.include_path = r("include")

  def check_installed_env(self, function, *args, **kwargs):
    """
    Loads the installed environment, if available, and runs the function
    in that environment. This is used in the local environment to check
    the installed environment.
    """
    installed_env = get_installed_env()
    if installed_env is not None and not self.installed:
      return getattr(installed_env, function)(*args, **kwargs)
    return None

  def check_local_env(self, function, *args, **kwargs):
    """
    Loads the local environment, if available, and runs the function
    in that environment. This is used in the installed environment to
    check the local environment.
    """
    local_env = get_local_env()
    if local_env is not None and self.installed:
      return getattr(local_env, function)(*args, **kwargs)
    return None

  def get_installed_module_path(self, name):
    """
    Returns the path to the installed path of a module. The path may not
    always be stored in the installed environment. None is returned if
    the module is not installed.
    """
    module_path = None
    installed_env = get_installed_env()
    if installed_env is not None:
      for p in installed_env.repository_paths:
        installed_path = os.path.join(abs(p), name)
        if os.path.isdir(installed_path):
          module_path = installed_path
          break
    return module_path

  def module_is_installed(self, name):
    """
    Returns True if module is in installed environment, False otherwise
    """
    module_path = self.get_installed_module_path(name)
    return module_path is not None

  def under_root(self, path, return_relocatable_path=False):
    return abs(self.build_path / '..' / path)

  def under_base(self, path, return_relocatable_path=False):
    result = self.build_path / '..' / 'conda_base' / path
    if self.installed:
      result = self.bin_path.anchor / path
    elif not self.build_options.use_conda:
      result = self.build_path / '..' / 'base' / path
    if not return_relocatable_path:
      result = abs(result)
    return result

  def under_build(self, path, return_relocatable_path=False):
    result = self.as_relocatable_path(path)
    if return_relocatable_path:
      return result
    else:
      return abs(result)

  def under_dist(self, module_name, path, default=KeyError, test=None,
                 return_relocatable_path=False):
    # check installed environment first
    result = self.check_installed_env('_under_dist', module_name, path, None, test, return_relocatable_path)
    if result is None:
      result = self.get_installed_module_path(module_name)
      if result is not None:
        result = self.as_relocatable_path(result) / path
        if not return_relocatable_path:
          result = abs(result)
    if result is not None:
      return result
    # check current environment
    result = self._under_dist(module_name, path, None, test, return_relocatable_path)
    if result is not None:
      return result
    # then check local environment
    result = self.check_local_env('_under_dist', module_name, path, default, test, return_relocatable_path)
    return result

  def _under_dist(self, module_name, path, default=KeyError, test=None,
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
    # check installed environment first
    result = self.check_installed_env('_dist_path', module_name, None, return_relocatable_path)
    if result is None:
      result = self.get_installed_module_path(module_name)
      if result is not None:
        result = self.as_relocatable_path(result)
        if not return_relocatable_path:
          result = abs(result)
    if result is not None:
      return result
    # check current environment
    result = self._dist_path(module_name, None, return_relocatable_path)
    if result is not None:
      return result
    # then check local environment
    result = self.check_local_env('_dist_path', module_name, default, return_relocatable_path)
    return result

  def _dist_path(self, module_name, default=KeyError,
                 return_relocatable_path=False):
    if (default is KeyError):
      result = self.module_dist_paths[module_name]
    else:
      result = self.module_dist_paths.get(module_name, default)
    if (isinstance(result, relocatable_path) and not return_relocatable_path):
      result = abs(result)
    return result

  def has_module(self, name):
    # check installed environment first
    result = self.check_installed_env('_has_module', name)
    if result:
      return result
    # check current environment
    result = self._has_module(name)
    if result:
      return result
    # then check local environment
    result = self.check_local_env('_has_module', name)
    return result

  def _has_module(self, name):
    return name in self.module_dist_paths

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
    print(self.command_version_suffix, file=f)

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
    # check installed environment first
    result = self.check_installed_env(
      '_find_in_repositories', relative_path, test, True, return_relocatable_path)
    if result is not None:
      return result
    # check current environment
    result = self._find_in_repositories(relative_path, test, True, return_relocatable_path)
    if result is not None:
      return result
    # then check local environment
    result = self.check_local_env(
      '_find_in_repositories', relative_path, test, optional, return_relocatable_path)
    return result

  def _find_in_repositories(self,
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
    result = None

    if module_name=='amber': return result # because amber_adaptbx is not an adapter

    # check installed environment first
    result = self.check_installed_env('_find_dist_path', module_name, True, return_relocatable_path)
    if result is not None:
      return result
    # check current environment
    result = self._find_dist_path(module_name, True, return_relocatable_path)
    if result is not None:
      return result
    # then check local environment
    result = self.check_local_env('_find_dist_path', module_name, optional, return_relocatable_path)
    return result

  def _find_dist_path(self, module_name, optional=False,
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
    if self.module_is_installed(module_name):
      print("{module_name} is already installed".format(module_name=module_name))
      installed_env = get_installed_env()
      if module_name in installed_env.module_dict:
        self.installed_modules.append(installed_env.module_dict[module_name])
      else:
        installed_module = module(env=self, name=module_name, dist_path=dist_path)
        self.installed_modules.append(installed_module)
      return True
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
    path = self.as_relocatable_path(path)
    if path not in self.repository_paths:
      self.repository_paths.append(path)

  def process_args(self, pre_processed_args):
    command_line = pre_processed_args.command_line
    self.no_bin_python = command_line.options.no_bin_python
    for path in pre_processed_args.repository_paths:
      if not op.isabs(path):
        path = abs(self.build_path / path)
      self.add_repository(path=path)
    module_names = []
    installed_env = get_installed_env()
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
        dist_path = relocatable_path(
          self.build_path, op.expandvars(redirection))
        if not dist_path.isdir():
          raise RuntimeError(
            'Invalid command line redirection:\n'
            '  module name = "%s"\n'
            '  redirection = "%s"\n'
            '  resulting target = "%s"' % (
              module_name, redirection, dist_path))
        self.command_line_redirections[module_name] = dist_path
      # remove any trailing dir separators
      module_name = module_name.rstrip("/\\")
      module_names.append(module_name)
    if (pre_processed_args.warm_start):
      self.explicitly_requested_modules |= set(module_names)
      if (not command_line.options.only):
        excludes = []
        if (command_line.options.exclude):
          excludes = command_line.options.exclude.split(",")
        self.explicitly_requested_modules -= set(excludes)
        for module in self.module_list:
          if (not module.name in excludes):
            module_names.append(module.name)
      else:
        self.explicitly_requested_modules = set(module_names)
    else:
      self.explicitly_requested_modules = set(module_names)
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
        enable_openmp_if_possible
          =command_line.options.enable_openmp_if_possible,
        enable_boost_threads
          =command_line.options.enable_boost_threads,
        enable_cuda=command_line.options.enable_cuda,
        enable_kokkos=command_line.options.enable_kokkos,
        use_conda=command_line.options.use_conda,
        opt_resources=command_line.options.opt_resources,
        use_environment_flags=command_line.options.use_environment_flags,
        force_32bit=command_line.options.force_32bit,
        msvc_arch_flag=command_line.options.msvc_arch_flag,
        enable_cxx11=command_line.options.enable_cxx11,
        cxxstd=command_line.options.cxxstd,
        skip_phenix_dispatchers=command_line.options.skip_phenix_dispatchers)
      self.build_options.get_flags_from_environment()
      # if an installed environment exists, override with build_options
      # from installed environment
      if installed_env is not None:
        self.build_options = installed_env.build_options
      if (self.build_options.use_conda):
        get_conda_prefix()
      if (command_line.options.command_version_suffix is not None):
        self.command_version_suffix = \
          command_line.options.command_version_suffix
        self.write_command_version_suffix()
    if (command_line.options.build_boost_python_extensions is not None):
      self.build_options.build_boost_python_extensions \
        = command_line.options.build_boost_python_extensions
    self.reset_module_registry()
    module_names.insert(0, "libtbx")

    # check for installed environment and remove installed modules
    if installed_env is not None and not self.installed:
      module_set = set(module_names)
      installed_module_names = set([module.name for module in self.installed_modules])
      for module_name in module_names:
        dist_path = installed_env.get_installed_module_path(module_name)
        if module_name in installed_env.module_dict \
          or dist_path is not None:
          try:
            module_set.remove(module_name)
            if module_name not in installed_module_names \
              and module_name != 'boost':
              from libtbx.env_config import module
              if dist_path is not None:
                dist_path = installed_env.as_relocatable_path(dist_path)
              installed_module = module(env=self, name=module_name, dist_path=dist_path)
              self.installed_modules.append(installed_module)
              installed_module_names.add(module_name)
          except (KeyError, ValueError):
            pass
      module_names = list(module_set)
      if len(module_names) == 0:
        print("All modules have already been installed. No new configuration is necessary.")
        sys.exit()

    for module_name in module_names:
      self.process_module(
        dependent_module=None, module_name=module_name, optional=False)
    self.scons_dist_path = self.find_dist_path("scons", optional=True)
    if self.scons_dist_path is None:
      # try to use SCons in Python installation
      try:
        import SCons
      except ImportError:
        pass
      else:
        self.scons_dist_path = self.as_relocatable_path(sys.prefix)
    if installed_env is not None and not self.installed:
      self.path_utility = installed_env.under_dist(
        "libtbx", "command_line/path_utility.py",
        return_relocatable_path=True)
    else:
      self.path_utility = self.under_dist(
        "libtbx", "command_line/path_utility.py",
        return_relocatable_path=True)
    assert self.path_utility.isfile()

  def dispatcher_precall_commands(self):
    if (self._dispatcher_precall_commands is None):
      lines = []
      self._dispatcher_precall_commands = lines
    return self._dispatcher_precall_commands

  def write_dispatcher_include_template(self):
    if (os.name == "nt"):
      print("    dispatcher_include_template.bat")
      f = self.under_build("dispatcher_include_template.bat",
                           return_relocatable_path=True).open("w")
      cp = "@REM"
    else :
      print("    dispatcher_include_template.sh")
      f = self.under_build("dispatcher_include_template.sh",
                           return_relocatable_path=True).open("w")
      cp = "#"
    print(cp+" include at start", file=f)
    print(cp+"   Commands to be executed at the start of the", file=f)
    print(cp+"   auto-generated dispatchers in bin.", file=f)
    print(cp+"", file=f)
    print(cp+" include before command", file=f)
    print(cp+"   Commands to be executed before the target command", file=f)
    print(cp+"   is called by the auto-generated dispatchers in bin.", file=f)
    print(cp+"", file=f)
    print(cp+" To see how the dispatchers work, look at an example:", file=f)
    print(cp+"   %s" % show_string(self.under_build("bin/libtbx.help")), file=f)
    print(cp+"", file=f)
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
          and ((file_name.endswith(".sh") and os.name != "nt") or
               (file_name.endswith(".bat") and os.name == "nt"))
          and file_name not in ["dispatcher_include_template.sh",
                                "dispatcher_include_template.bat"]):
        include_files.append(path)
    include_files.sort()
    for path in include_files:
      print("Processing: %s" % show_string(abs(path)))
      lines = path.open().read().splitlines()
      lines_at_start = []
      lines_before_command = []
      buffer = lines_before_command
      for line in lines:
        l = " ".join(line.split()).lower()
        if (l.startswith("#include ") or l.startswith("REM include")):
          l = "# " + l[1:]
        if   (l in ["# include at start", "REM include at start"]):
          buffer = lines_at_start
        elif (l in ["# include before command", "REM include before command"]):
          buffer = lines_before_command
        else:
          buffer.append(line)
      for buffer,target in [(lines_at_start,
                             self._dispatcher_include_at_start),
                            (lines_before_command,
                             self._dispatcher_include_before_command)]:
        if (len(buffer) == 0): continue
        if (os.name == "nt"):
          buffer.insert(0, "@REM included from %s" % abs(path))
        else :
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
        "--opt_resources not supported in combination with --compiler=icc")
    def get_libs_dir():
      if (sys.platform.startswith("linux")):
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

  def ld_library_path_var_name(self):
    if (os.name == "nt"):
      return "PATH"
    if (sys.platform.startswith("darwin")):
      if (self.build_options.use_conda):
        return "DYLD_FALLBACK_LIBRARY_PATH"
      return "DYLD_LIBRARY_PATH"
    else:
      return "LD_LIBRARY_PATH"

  def ld_library_path_additions(self):
    result = [self.lib_path]
    dirs = []
    if (is_64bit_architecture()):
      dirs.append("lib64")
    dirs.append("lib")
    for d in dirs:
      p = op.join(sys.prefix, d)
      if (op.isdir(p)):
        result.append(self.as_relocatable_path(p))
    return result

  def write_conda_dispatcher(self, source_file, target_file,
                             source_is_python_exe=False):
    '''
    Simplified dispatcher for conda package since many of the environment
    variables are no longer necessary.
    '''
    with target_file.open('w') as f:
      if sys.platform == 'win32':  # windows
        print('@setlocal', file=f)
        print('@set LIBTBX_PREFIX=%~dp0', file=f)
        print('@set LIBTBX_PREFIX=%LIBTBX_PREFIX:~0,-1%', file=f)
        print(r'@for %%F in ("%LIBTBX_PREFIX%") do @set LIBTBX_PREFIX=%%~dpF', file=f)
        print('@set LIBTBX_PREFIX=%LIBTBX_PREFIX:~0,-1%', file=f)
        print('@set LIBTBX_DISPATCHER_NAME=%~nx0', file=f)
        print('@set PATH=%LIBTBX_PREFIX%\\..;%LIBTBX_PREFIX%\\..\\mingw-w64\\bin;%LIBTBX_PREFIX%\\..\\bin;%LIBTBX_PREFIX%\\..\\..\\Scripts;%PATH%', file=f)
        def write_dispatcher_include(where):
          for line in self.dispatcher_include(where=where):
            if (line.startswith("@")):
              print(line, file=f)
            else :
              print("@" + line, file=f)
        write_dispatcher_include(where="at_start")
        print('@set LIBTBX_PYEXE=%s' % self.python_exe.bat_value(anchor_var='LIBTBX_PREFIX'), file=f)
        write_dispatcher_include(where="before_command")
        qnew_tmp = qnew
        if self.python_version_major_minor[0] == 3:
          qnew_tmp = '' # -Q is gone in Python3.
        if source_file.ext().lower() == '.py':
          print('@"%%LIBTBX_PYEXE%%"%s "%s" %%*' % (
            qnew_tmp, source_file.bat_value(anchor_var='LIBTBX_PREFIX')), file=f)
        elif source_file.basename().lower() == 'python.exe':
          print('@"%%LIBTBX_PYEXE%%"%s %%*' % qnew_tmp, file=f)
        else:
          print('@"%s" %%*' % source_file.bat_value(anchor_var='LIBTBX_PREFIX'), file=f)
      else:  # linux and macOS
        if (source_file is not None):
          print('#! /bin/sh', file=f)
          print('# LIBTBX_DISPATCHER DO NOT EDIT', file=f)
        else:
          print('# LIBTBX_DISPATCHER_HEAD DO NOT EDIT', file=f)
          print('#', file=f)
          print('# This file is intended to be sourced from other scripts.', file=f)
          print('# It is like the dispatcher scripts in the bin directory,', file=f)
          print('# but only sets up the environment without calling a', file=f)
          print('# command at the end.', file=f)
        print('#', file=f)
        write_do_not_edit(f=f)
        print('# To customize this auto-generated script create', file=f)
        print('#', file=f)
        print('#   dispatcher_include*.sh', file=f)
        print('#', file=f)
        print('# files in %s and run' % show_string(abs(self.build_path)), file=f)
        print('#', file=f)
        print('#   libtbx.refresh', file=f)
        print('#', file=f)
        print('# to re-generate the dispatchers (libtbx.refresh is a subset', file=f)
        print('# of the functionality of the libtbx/configure.py command).', file=f)
        print('#', file=f)
        print('# See also:', file=f)
        print('#   %s' \
          % show_string(self.under_build("dispatcher_include_template.sh")), file=f)
        print('#', file=f)
        print(_SHELLREALPATH_CODE, file=f)
        print('LIBTBX_PREFIX="$(shellrealpath "$0" && cd "$(dirname "$RESULT")/.." && pwd)"', file=f)
        print('export LIBTBX_PREFIX', file=f)
        print('LIBTBX_PYEXE_BASENAME="%s"' % self.python_exe.basename(), file=f)
        print('export LIBTBX_PYEXE_BASENAME', file=f)
        print('# Set the CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES environment variable', file=f)
        print('# if you want python to use the following environment variables.', file=f)
        print('# Otherwise, this environment takes priority at runtime.', file=f)
        print('if [ -z "${CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES}" ]; then', file=f)
        print('  unset PYTHONHOME', file=f)
        print('  unset PYTHONPATH', file=f)
        print('  unset LD_LIBRARY_PATH', file=f)
        print('  unset DYLD_LIBRARY_PATH', file=f)
        print('  unset DYLD_FALLBACK_LIBRARY_PATH', file=f)
        print('  export PATH="${LIBTBX_PREFIX}/bin:${PATH}"', file=f)
        print('fi', file=f)
        source_is_py = False
        if (source_file is not None):
          dispatcher_name = target_file.basename()
          if (dispatcher_name.find('"') >= 0):
            raise RuntimeError(
              "Dispatcher target file name contains double-quote: %s\n"
                % dispatcher_name
              + "  source file: %s" % source_file)
          print('LIBTBX_DISPATCHER_NAME="%s"' % target_file.basename(), file=f)
          print('export LIBTBX_DISPATCHER_NAME', file=f)
          if source_file.ext().lower() == ".py":
            source_is_py = True
          else:
            with open(abs(source_file), 'rb') as fh:
              first_line = fh.readline()
            if first_line.startswith(b'#!') and b'python' in first_line.lower():
              source_is_py = True
        for line in self.dispatcher_include(where="at_start"):
          print(line, file=f)

        precall_commands = self.dispatcher_precall_commands()
        if (precall_commands is not None):
          for line in precall_commands:
            print(line, file=f)
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
            print(line, file=f)
        for line in self.dispatcher_include(where="before_command"):
          print(line, file=f)
        if (scan_for_dispatcher_includes):
          for line in source_specific_dispatcher_include(
                        pattern="LIBTBX_POST_DISPATCHER_INCLUDE_SH",
                        source_file=source_file):
            print(line, file=f)
        if (self.build_options.opt_resources):
          ldpl = self.opt_resources_ld_preload()
          if (ldpl is not None):
            print('if [ "${LIBTBX_NO_LD_PRELOAD-UNSET}" == UNSET ]; then', file=f)
            print('  LD_PRELOAD="%s"' % ldpl, file=f)
            print('  export LD_PRELOAD', file=f)
            print('fi', file=f)
        print('LIBTBX_PYEXE="%s"' % (
          self.python_exe.dirname() / "$LIBTBX_PYEXE_BASENAME").sh_value(anchor_var='LIBTBX_PREFIX'), file=f)
        print('export LIBTBX_PYEXE', file=f)

        if (source_file is not None):
          cmd = ""
          if (source_is_py or source_is_python_exe):
            qnew_tmp = qnew
            if self.python_version_major_minor[0] == 3:
              qnew_tmp = '' # -Q is gone in Python3.
            cmd += ' "$LIBTBX_PYEXE"%s' % qnew_tmp
          start_python = False
          if (source_is_py):
            if (len(source_specific_dispatcher_include(
                      pattern="LIBTBX_START_PYTHON",
                      source_file=source_file)) > 3):
              start_python = True
          if (not start_python and not source_is_python_exe):
            cmd += ' "%s"' % source_file.sh_value(anchor_var='LIBTBX_PREFIX')
          print('if [ -n "$LIBTBX__VALGRIND_FLAG__" ]; then', file=f)
          print("  exec $LIBTBX_VALGRIND"+cmd, '"$@"', file=f)
          tmp_reloc = os.path.basename(source_file.relocatable)
          if tmp_reloc.endswith('.py') and cmd.find('-Qnew')>-1:
            print('elif [ -n "$LIBTBX__CPROFILE_FLAG__" ]; then', file=f)
            print('  exec %s "$@"' % cmd.replace(
              '-Qnew',
              '-Qnew -m cProfile -o %s.profile' % os.path.basename(target_file.relocatable),
              ), file=f)
          print("elif [ $# -eq 0 ]; then", file=f)
          print("  exec"+cmd, file=f)
          print("else", file=f)
          print("  exec"+cmd, '"$@"', file=f)
          print("fi", file=f)
        target_file.chmod(0o755)

  def write_bin_sh_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):

    # set SSL_CERT_FILE if certifi is available
    cert_file = None
    try:
      import certifi
      cert_file = self.as_relocatable_path(certifi.where())
    except ImportError:
      pass

    # set OPENBLAS_NUM_THREADS
    # This prevents the binary pip installation of numpy from spawning threads
    # when flex is imported ("from scitbx.array_family import flex")
    # Problem seems to only appear on linux, but should not hurt other operating
    # systems.
    # According to https://github.com/xianyi/OpenBLAS#usage , there are 3
    # environment variables that control the spawning of threads
    #   OPENBLAS_NUM_THREADS
    #   GOTO_NUM_THREADS
    #   OMP_NUM_THREADS
    # where OPENBLAS_NUM_THREADS > GOTO_NUM_THREADS > OMP_NUM_THREADS in terms
    # of precedence. Setting OPENBLAS_NUM_THREADS should allow OMP_NUM_THREADS
    # to be used for OpenMP sections. Using multiple threads for numpy functions
    # that depend on OpenBLAS will require changing the OPENBLAS_NUM_THREADS
    # environment variable in the code that wants that functionality.
    openblas_num_threads = 1

    # determine LC_ALL from environment (Python UTF-8 compatibility in Linux)
    LC_ALL = os.environ.get('LC_ALL')     # user setting
    if (LC_ALL is not None):
      if ( ('UTF-8' not in LC_ALL) and ('utf8' not in LC_ALL) ):
        LC_ALL = None
    if (LC_ALL is None):
      LC_ALL = 'en_US.UTF-8'              # default

    f = target_file.open("w")
    if (source_file is not None):
      print('#! /bin/sh', file=f)
      print('# LIBTBX_DISPATCHER DO NOT EDIT', file=f)
    else:
      print('# LIBTBX_DISPATCHER_HEAD DO NOT EDIT', file=f)
      print('#', file=f)
      print('# This file is intended to be sourced from other scripts.', file=f)
      print('# It is like the dispatcher scripts in the bin directory,', file=f)
      print('# but only sets up the environment without calling a', file=f)
      print('# command at the end.', file=f)
    print('#', file=f)
    write_do_not_edit(f=f)
    print('# To customize this auto-generated script create', file=f)
    print('#', file=f)
    print('#   dispatcher_include*.sh', file=f)
    print('#', file=f)
    print('# files in %s and run' % show_string(abs(self.build_path)), file=f)
    print('#', file=f)
    print('#   libtbx.refresh', file=f)
    print('#', file=f)
    print('# to re-generate the dispatchers (libtbx.refresh is a subset', file=f)
    print('# of the functionality of the libtbx/configure.py command).', file=f)
    print('#', file=f)
    print('# See also:', file=f)
    print('#   %s' \
      % show_string(self.under_build("dispatcher_include_template.sh")), file=f)
    print('#', file=f)
    print(_SHELLREALPATH_CODE, file=f)
    print('unset PYTHONHOME', file=f)
    print('LC_ALL=' + LC_ALL, file=f)
    print('export LC_ALL', file=f)
    print('LIBTBX_BUILD="$(shellrealpath "$0" && cd "$(dirname "$RESULT")/.." && pwd)"', file=f)
    print('export LIBTBX_BUILD', file=f)
    print('LIBTBX_PYEXE_BASENAME="%s"' % self.python_exe.basename(), file=f)
    print('export LIBTBX_PYEXE_BASENAME', file=f)
    source_is_py = False
    if (source_file is not None):
      dispatcher_name = target_file.basename()
      if (dispatcher_name.find('"') >= 0):
        raise RuntimeError(
          "Dispatcher target file name contains double-quote: %s\n"
            % dispatcher_name
          + "  source file: %s" % source_file)
      print('LIBTBX_DISPATCHER_NAME="%s"' % target_file.basename(), file=f)
      print('export LIBTBX_DISPATCHER_NAME', file=f)
      if source_file.ext().lower() == ".py":
        source_is_py = True
      else:
        with open(abs(source_file), 'rb') as fh:
          first_line = fh.readline()
        if first_line.startswith(b'#!') and b'python' in first_line.lower():
          source_is_py = True
    for line in self.dispatcher_include(where="at_start"):
      print(line, file=f)
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((
      self.ld_library_path_var_name(),
      self.ld_library_path_additions()))
    essentials.append(("PATH", [self.bin_path]))

    if (cert_file is not None):
      print('SSL_CERT_FILE="%s"' % cert_file.sh_value(), file=f)
      print('export SSL_CERT_FILE', file=f)

    print('OPENBLAS_NUM_THREADS="%s"' % openblas_num_threads, file=f)
    print('export OPENBLAS_NUM_THREADS', file=f)

    if sys.version_info[0] == 2:
      pangorc = abs(self.build_path / '..' / 'base' / 'etc' / 'pango' / 'pangorc')
      if os.path.exists(pangorc):
        print('PANGO_RC_FILE="$LIBTBX_BUILD/../base/etc/pango/pangorc"', file=f)
        print('export PANGO_RC_FILE', file=f)
      fontconfig = abs(self.build_path / '..' / 'base' / 'etc' / 'fonts')
      if os.path.exists(fontconfig):
        print('FONTCONFIG_PATH="$LIBTBX_BUILD/../base/etc/fonts"', file=f)
        print('export FONTCONFIG_PATH', file=f)

      # set paths for fontconfig and gdk-pixbuf due to gtk2
      # checks the location of the conda environment
      if self.build_options.use_conda and sys.platform.startswith('linux'):
        fontconfig_path = '{conda_base}/etc/fonts'
        fontconfig_file = '$FONTCONFIG_PATH/fonts.conf'
        gdk_pixbuf_module_file = '{conda_base}/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache'
        gtk_path = '{conda_base}/lib/gtk-2.0/2.10.0'
        conda_base = get_conda_prefix()
        if (os.path.exists(fontconfig_path.format(conda_base=conda_base)) and
            os.path.exists(gdk_pixbuf_module_file.format(conda_base=conda_base))):
          if conda_base == abs(self.build_path / '..' / 'conda_base'):
            conda_base = '$LIBTBX_BUILD/../conda_base'
          print('unset GTK_MODULES', file=f)
          print('unset GTK2_RC_FILES', file=f)
          print('unset GTK_RC_FILES', file=f)
          print('export FONTCONFIG_PATH=' +
                fontconfig_path.format(conda_base=conda_base), file=f)
          print('export FONTCONFIG_FILE=' + fontconfig_file, file=f)
          print('export GDK_PIXBUF_MODULE_FILE=' +
                gdk_pixbuf_module_file.format(conda_base=conda_base), file=f)
          print('export GTK_PATH=' + gtk_path.format(conda_base=conda_base), file=f)

    for n,v in essentials:
      if (len(v) == 0): continue
      v = ":".join([p.sh_value() for p in v])
      print('if [ -n "$%s" ]; then' % n, file=f)
      print('  %s="%s:$%s"' % (n, v, n), file=f)
      print('  export %s' % n, file=f)
      print('else', file=f)
      print('  %s="%s"' % (n, v), file=f)
      print('  export %s' % n, file=f)
      print('fi', file=f)
    precall_commands = self.dispatcher_precall_commands()
    if (precall_commands is not None):
      for line in precall_commands:
        print(line, file=f)
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
        print(line, file=f)
    for line in self.dispatcher_include(where="before_command"):
      print(line, file=f)
    if (scan_for_dispatcher_includes):
      for line in source_specific_dispatcher_include(
                    pattern="LIBTBX_POST_DISPATCHER_INCLUDE_SH",
                    source_file=source_file):
        print(line, file=f)
    if (self.build_options.opt_resources):
      ldpl = self.opt_resources_ld_preload()
      if (ldpl is not None):
        print('if [ "${LIBTBX_NO_LD_PRELOAD-UNSET}" == UNSET ]; then', file=f)
        print('  LD_PRELOAD="%s"' % ldpl, file=f)
        print('  export LD_PRELOAD', file=f)
        print('fi', file=f)
    print('LIBTBX_PYEXE="%s"' % (
      self.python_exe.dirname() / "$LIBTBX_PYEXE_BASENAME").sh_value(), file=f)
    print('export LIBTBX_PYEXE', file=f)

    # Since El Capitan, Apple Python does not allow relative rpath in shared
    # libraries. Thus any cctbx-based script will fail with an import error
    # because each and every .so references lib/libboost_python.dylib.
    # There are two ways to fix this issue:
    # - we could put the absolute path at compile-time
    # - or we can fix it at runtime
    # The former would make the cctbx build directory non-relocatable. Thus
    # we prefer the latter.
    if sys.platform == "darwin" and abs(self.python_exe).startswith("/usr/bin/"):
      print("""/usr/bin/perl <<'FIXRPATH'
        for $lib(<$ENV{LIBTBX_BUILD}/lib/*.so>) {
            open OTOOL, "-|", "otool", "-L", $lib;
            while(<OTOOL>) {
                m{^\\s+(lib\\S+)} and $libs{$lib}{$1}++;
            }
        }
        while(($so, $relative_libs) = each %libs) {
            for $lib(keys %$relative_libs) {
                system "install_name_tool",
                       "-change", $lib, "\\@loader_path/../$lib", $so;
            }
        }
      """, file=f)
      print("FIXRPATH", file=f)
    if (source_file is not None):
      def pre_cmd():
        return ['', '/usr/bin/arch -i386 '][self.build_options.force_32bit]
      cmd = ""
      if (source_is_py or source_is_python_exe):
        qnew_tmp = qnew
        if self.python_version_major_minor[0] == 3:
          qnew_tmp = '' # -Q is gone in Python3.
        cmd += ' %s"$LIBTBX_PYEXE"%s' % (pre_cmd(), qnew_tmp)
      start_python = False
      if (source_is_py):
        if (len(source_specific_dispatcher_include(
                  pattern="LIBTBX_START_PYTHON",
                  source_file=source_file)) > 3):
          start_python = True
      if (not start_python and not source_is_python_exe):
        cmd += ' "%s"' % source_file.sh_value()
      print('if [ -n "$LIBTBX__VALGRIND_FLAG__" ]; then', file=f)
      print("  exec $LIBTBX_VALGRIND"+cmd, '"$@"', file=f)
      tmp_reloc = os.path.basename(source_file.relocatable)
      if tmp_reloc.endswith('.py') and cmd.find('-Qnew')>-1:
        print('elif [ -n "$LIBTBX__CPROFILE_FLAG__" ]; then', file=f)
        print('  exec %s "$@"' % cmd.replace(
          '-Qnew',
          '-Qnew -m cProfile -o %s.profile' % os.path.basename(target_file.relocatable),
          ), file=f)
      print("elif [ $# -eq 0 ]; then", file=f)
      print("  exec"+cmd, file=f)
      print("else", file=f)
      print("  exec"+cmd, '"$@"', file=f)
      print("fi", file=f)
    f.close()
    target_file.chmod(0o755)

  def write_win32_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):
    f = target_file.open('w')
    # By default, changes to environment variables are  permanent on Windows,
    # i.e. it is as if export VAR was added after each set VAR=...
    # As a result, e.g. set PYTHONPATH=...; %PYTHONPATH% results in growing
    # PYTHONPATH each time a dispatcher script is run.
    # Thus setlocal essential (endlocal is implied)
    print('@setlocal', file=f)
    print('@set LIBTBX_BUILD=%~dp0', file=f)
    print('@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%', file=f)
    print(r'@for %%F in ("%LIBTBX_BUILD%") do @set LIBTBX_BUILD=%%~dpF', file=f)
    print('@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%', file=f)
    print('@set LIBTBX_DISPATCHER_NAME=%~nx0', file=f)
    def write_dispatcher_include(where):
      for line in self.dispatcher_include(where=where):
        if (line.startswith("@")):
          print(line, file=f)
        else :
          print("@" + line, file=f)
    write_dispatcher_include(where="at_start")
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((self.ld_library_path_var_name(), [self.lib_path]))
    bin_path = [self.bin_path]
    if self.build_options.use_conda:
      bin_path.append(self.as_relocatable_path(get_conda_prefix()) / 'Library' / 'bin')
      bin_path.append(self.as_relocatable_path(get_conda_prefix()) / 'Library' / 'mingw-w64' / 'bin')
    essentials.append(("PATH", bin_path))
    for n,v in essentials:
      if (len(v) == 0): continue
      v = ';'.join([ op.join('%LIBTBX_BUILD%', p.relocatable) for p in v ])
      print('@set %s=%s;%%%s%%' % (n, v, n), file=f)
    print('@set LIBTBX_PYEXE=%s' % self.python_exe.bat_value(), file=f)
    write_dispatcher_include(where="before_command")
    qnew_tmp = qnew
    if self.python_version_major_minor[0] == 3:
      qnew_tmp = '' # -Q is gone in Python3.
    if source_file.ext().lower() == '.py':
      print('@"%%LIBTBX_PYEXE%%"%s "%s" %%*' % (
        qnew_tmp, source_file.bat_value()), file=f)
    elif source_file.basename().lower() == 'python.exe':
      print('@"%%LIBTBX_PYEXE%%"%s %%*' % qnew_tmp, file=f)
    else:
      print('@"%s" %%*' % source_file.bat_value(), file=f)
    f.close()

  def write_dispatcher(self,
        source_file, target_file, source_is_python_exe=False):
    source_file = self.as_relocatable_path(source_file)
    target_file = self.as_relocatable_path(target_file)
    # always skip amber.python
    if target_file.basename() in ['amber.python',
                                  #'rosetta.python',
                                  #'afitt.python',
                                 ]:
      return
    if "phenix" not in self.module_dict and self.build_options.skip_phenix_dispatchers and "phenix" in target_file.basename():
      return
    reg = self._dispatcher_registry.setdefault(target_file, source_file)
    if reg != source_file:
      if not reg.isfile():
        self._dispatcher_registry[target_file] = source_file
      elif (source_file.isfile()
            and (   not hasattr(os.path, "samefile")
                 or not reg.samefile(source_file))):
        raise RuntimeError("Multiple sources for dispatcher:\n"
          + "  target file:\n"
          + "    %s\n" % show_string(abs(target_file))
          + "  source files:\n"
          + "    %s\n" % show_string(abs(reg))
          + "   =%s\n" % reg
          + "    %s\n" % show_string(abs(source_file))
          + "   =%s" % source_file)
    if abs(self.build_path) == os.path.abspath(get_conda_prefix()) or \
      (os.name == "nt" and abs(self.build_path).lower().endswith('library')):
      action = self.write_conda_dispatcher
    elif (os.name == "nt"):
      action = self.write_win32_dispatcher
    else:
      action = self.write_bin_sh_dispatcher
    target_file_ext = target_file
    if os.name == 'nt':
      target_file_ext += '.bat'
    remove_or_rename(target_file_ext)
    try: action(source_file, target_file_ext, source_is_python_exe)
    except IOError as e: print("  Ignored:", e)

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
    print("   ", target_file)
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
if [ -n "$LIBTBX_BUILD_RELOCATION_HINT" ]; then
  cd "$LIBTBX_BUILD_RELOCATION_HINT"
  LIBTBX_BUILD_RELOCATION_HINT=
  export LIBTBX_BUILD_RELOCATION_HINT
elif [ -n "$BASH_SOURCE" ]; then
  LIBTBX_BUILD=`dirname "$BASH_SOURCE[0]"`
  cd "$LIBTBX_BUILD"
else
  cd "%s"
fi
LIBTBX_BUILD=`pwd -P`
export LIBTBX_BUILD
LIBTBX_OPATH="$PATH"
export LIBTBX_OPATH
PATH="$LIBTBX_BUILD/bin:$PATH"
export PATH
cd "$ocwd"
ocwd=
""" % abs(self.build_path))
    s.write("""\
alias libtbx.setpaths_all=". \\"$LIBTBX_BUILD/setpaths_all.sh\\""
alias libtbx.unsetpaths=". \\"$LIBTBX_BUILD/unsetpaths.sh\\""
""")
    print('unalias libtbx.unsetpaths > /dev/null 2>&1', file=u)
    if (self.is_development_environment()):
      print('''alias cdlibtbxbuild="cd \\"$LIBTBX_BUILD\\""''', file=s)
      print('unalias cdlibtbxbuild > /dev/null 2>&1', file=u)
    setpaths.all_and_debug()
    setpaths.set_unset_vars()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.sh_value(),
      var_name_in="LIBTBX_OPATH")
    for f in s, u:
      print('LIBTBX_TMPVAL=', file=f)
      print('LIBTBX_OPATH=', file=f)
      if (suffix == ""):
        print('LIBTBX_BUILD=', file=f)

  def write_setpaths_csh(self, suffix):
    setpaths = unix_setpaths(self, "csh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      write_do_not_edit(f=f)
      f.write("""\
set ocwd="$cwd"
if ($?LIBTBX_BUILD_RELOCATION_HINT) then
  cd "$LIBTBX_BUILD_RELOCATION_HINT"
  unsetenv LIBTBX_BUILD_RELOCATION_HINT
else
  cd "%s"
endif
setenv LIBTBX_BUILD "`/bin/sh -c 'pwd -P'`"
setenv LIBTBX_OPATH "$PATH"
setenv PATH "$LIBTBX_BUILD/bin:$PATH"
cd "$ocwd"
unset ocwd
""" % abs(self.build_path))
    s.write("""\
alias libtbx.setpaths_all "source '$LIBTBX_BUILD/setpaths_all.csh'"
alias libtbx.unsetpaths "source '$LIBTBX_BUILD/unsetpaths.csh'"
""")
    print('unalias libtbx.unsetpaths', file=u)
    if (self.is_development_environment()):
      print('''alias cdlibtbxbuild "cd '$LIBTBX_BUILD'"''', file=s)
      print('unalias cdlibtbxbuild', file=u)
    setpaths.all_and_debug()
    setpaths.set_unset_vars()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.sh_value(),
      var_name_in="LIBTBX_OPATH")
    for f in s, u:
      print('unsetenv LIBTBX_TMPVAL', file=f)
      print('unsetenv LIBTBX_OPATH', file=f)
      if (suffix == ""):
        print('unsetenv LIBTBX_BUILD', file=f)

  def write_setpaths_bat(self, suffix):
    setpaths = windows_setpaths(self, suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      write_do_not_edit(f=f, win_bat=True)
      print(r'''@set LIBTBX_BUILD=%~dp0
@set LIBTBX_BUILD=%LIBTBX_BUILD:~0,-1%
@set LIBTBX_OPATH=%PATH%''', file=f)
      print('@set PATH=%s;%%PATH%%' % self.bin_path.bat_value(), file=f)
    setpaths.all_and_debug()
    setpaths.update_path(
      var_name="PATH",
      val=self.bin_path.bat_value(),
      var_name_in="LIBTBX_OPATH")
    for command in ["setpaths_all", "unsetpaths"]:
      print('@doskey libtbx.%s="%s\\%s.bat"' % (
        command, "%LIBTBX_BUILD%", command), file=s)
    print('@doskey libtbx.unsetpaths=', file=u)
    if (self.is_development_environment()):
      print('@doskey cdlibtbxbuild=cd "%LIBTBX_BUILD%"', file=s)
      print('@doskey cdlibtbxbuild=', file=u)
    if (suffix == "_debug"):
      print('@set PYTHONCASEOK=1', file=s) # no unset
    setpaths.set_unset_vars()
    for f in s, u:
      print('@set LIBTBX_OPATH=', file=f)
      if (suffix == ""):
        print('@set LIBTBX_BUILD=', file=f)

  def write_SConstruct(self):
    f = open_info(self.under_build("SConstruct",
                                   return_relocatable_path=True))
    write_do_not_edit(f=f)
    print('SConsignFile()', file=f)
    if os.getenv("SCONS_CACHE"):
      cache_dir = os.getenv("SCONS_CACHE")
      print("Using build cache in %s" % cache_dir)
      print("CacheDir(%r)" % cache_dir, file=f)
      assert os.path.exists(cache_dir) and os.path.isdir(cache_dir), \
        "Specified build cache dir does not exist"

    repository_paths = self.repository_paths
    repository_names = set([abs(p) for p in repository_paths])
    module_list = self.module_list
    module_names = set([module.name for module in module_list])
    # insert repositories and modules from installed environment
    installed_env = get_installed_env()
    if installed_env is not None:
      for repository_path in installed_env.repository_paths:
        if abs(repository_path) in repository_names:
          for p in repository_paths:
            if abs(p) == abs(repository_path):
              try:
                repository_paths.remove(p)
              except ValueError:
                pass
      repository_paths = installed_env.repository_paths + repository_paths
      # collect modules
      all_modules = installed_env.module_list
      local_env = get_local_env()
      if local_env is not None:
        all_modules += local_env.module_list
      all_modules += self.module_list
      # reorder modules
      module_list = []
      for ordered_name in installed_env.installed_order:
        for module in all_modules:
          if ordered_name == module.name:
            module_list.append(module)
            all_modules.remove(module)
            break
      # add remaining modules
      module_list_names = set([module.name for module in module_list])
      for module in all_modules:
        if module.name not in module_list_names:
          module_list.append(module)
          module_list_names.add(module.name)

    for path in repository_paths:
      print('Repository(r"%s")' % abs(path), file=f)
    for module in module_list:
      name,path  = list(module.name_and_dist_path_pairs())[-1]
      for script_name in ["libtbx_SConscript", "SConscript"]:
        if (path / script_name).isfile():
          print('SConscript("%s/%s")' % (name, script_name), file=f)
          break
    f.close()

  def write_Makefile(self):
    if (op.isfile("Makefile")):
      os.rename("Makefile", "Makefile.old")
        # make cja seems to get confused if the file is simply overwritten
    f = open_info(self.under_build("Makefile",
                                   return_relocatable_path=True))
    current_path = os.environ.get('PATH')
    lsj = './bin/libtbx.scons -j "`./bin/libtbx.show_number_of_processors`"'
    lc = './bin/libtbx.configure'
    if get_installed_env() is not None:
      lsj = 'libtbx.scons -j "`libtbx.show_number_of_processors`"'
      lc = 'libtbx.configure'
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
\t%(lc)s .
\t%(lsj)s

redo:
\t%(lc)s . --clear_scons_memory
\t%(lsj)s

clean:
\t%(lsj)s -c

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
      print("#! /bin/csh -f", file=f)
      print("set noglob", file=f)
      print("set verbose", file=f)
      for file_name in test_scripts:
        print('libtbx.python "%s" $*' % abs(file_name), file=f)
      f.close()
      path.chmod(0o755)

  def pickle(self):
    self.reset_dispatcher_support()
    file_name = self.build_path / "libtbx_env"
    with file_name.open("wb") as f:
      pickle.dump(self, f, 0)

  def show_module_listing(self):
    print("Relocatable paths anchored at: %s" % abs(self.build_path))
    print("Top-down list of all modules involved:")
    top_down_module_list = list(self.module_list)
    top_down_module_list.reverse()
    labels = [module.names_for_module_listing()
      for module in top_down_module_list]
    if (len(labels) == 0): return
    fmt = "  %%-%ds  %%s" % max([len(label) for label in labels])
    for label,module in zip(labels,top_down_module_list):
      for dist_path in module.dist_paths_active():
        print(fmt % (label, show_string(abs(dist_path))))
        label = ""

  def show_build_options_and_module_listing(self):
    print("Python: %s %s" % (
      sys.version.split()[0], show_string(sys.executable)))
    if (self.is_ready_for_build()):
      self.build_options.report()
    print("command_version_suffix:", self.command_version_suffix)
    self.show_module_listing()
    if (len(self.missing_for_use) > 0):
      raise RuntimeError("Missing modules: "
        + " ".join(sorted(self.missing_for_use)))
    if (not self.is_ready_for_build()):
      if (self.scons_dist_path is not None):
        print("***********************************")
        print("Warning: modules missing for build:")
        for module_name in sorted(self.missing_for_build):
          print(" ", module_name)
        print("***********************************")
      remove_or_rename(self.under_build("SConstruct"))

  def write_setpath_files(self):
    for suffix in ["", "_all", "_debug"]:
      if (sys.platform == "win32"):
        self.write_setpaths_bat(suffix)
      else:
        self.write_setpaths_sh(suffix)
        self.write_setpaths_csh(suffix)

  def process_exe(self):
    for path in [self.exe_path,
                 self.under_build("exe_dev", return_relocatable_path=True)]:
      if path.isdir():
        print("Processing: %s" % show_string(abs(path)))
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

  @staticmethod
  def get_setuptools_script_dir():
    '''
    Find the location of python entry point console_scripts, ie. things like
    'pip', 'pytest', ...
    This is different from simple /base/bin, eg. on MacOS.

    https://stackoverflow.com/questions/25066084/get-entry-point-script-file-location-in-setuputils-package
    '''
    from setuptools import Distribution
    from setuptools.command.install import install
    class OnlyGetScriptPath(install):
      def run(self):
        # does not call install.run() by design
        self.distribution.install_scripts = self.install_scripts
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary, but to be safe
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts

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
    # explicitly add site-packages path to avoid conflicts with per user site-packages
    # https://www.python.org/dev/peps/pep-0370/
    if (hasattr(site, 'getsitepackages')):
      pythonpath = [self.as_relocatable_path(p) for p in site.getsitepackages()]
    else:
      # fallback in case of virtualenv
      # https://github.com/pypa/virtualenv/issues/355
      from distutils.sysconfig import get_python_lib
      pythonpath = [self.as_relocatable_path(get_python_lib())]
    pythonpath.append(self.lib_path)
    for module in self.module_list:
      pythonpath.extend(module.assemble_pythonpath())
    pythonpath.reverse()
    self.pythonpath = unique_paths(paths=pythonpath)

  def is_development_environment(self):
    for module in self.module_list:
      if (module.is_version_controlled()):
        return True
    return False

  def relocate_python_paths_if_necessary(self):
    base_directory = abs(self.build_path / '..' / 'base')
    def is_in_local_base_directory(path):
      return os.path.commonprefix([os.path.abspath(path), base_directory]) == base_directory

    if not is_in_local_base_directory(sys.executable):
      return # no relocation required when using non-local (presumably: system) python
    site_packages = site.getsitepackages()
    if any(map(is_in_local_base_directory, site_packages)):
      return # There is a site directory inside the local python path.
             # This means relocation is not required.

    print("libtbx.configure determined that a python relocation is required.")
    old_base_path = os.path.commonprefix(site_packages)
    while old_base_path:
      old_base_path, directory = os.path.split(old_base_path)
      if directory == 'base': break
    if directory != 'base':
      print("WARNING: Relocation not possible, could not determine original base path from '{}'".format(os.path.commonprefix(site_packages)))
      return
    old_base_path = os.path.join(old_base_path, 'base')
    print("Attempting relocate from: {}\n           relocate to  : {}".format(old_base_path, base_directory))

    import libtbx.auto_build.rpath
    libtbx.auto_build.rpath.run(['--otherroot', old_base_path, base_directory])

  def clear_scons_memory(self):
    (self.build_path / ".sconsign.dblite").remove()
    (self.build_path / ".sconf_temp").remove_tree()

  def refresh(self):
    # check if self is modifiable
    if self.installed:
      print("Installed environments cannot be refreshed.")
      try:
        env = unpickle(os.getcwd())
      except FileNotFoundError:
        env = None
      if env is not None:
        print("Updating environment in {}".format(os.getcwd()))
        env.refresh()
        env.pickle()  # need a separate call
    # continue normally
    else:
      completed_file_name = (self.build_path / "libtbx_refresh_is_completed")
      completed_file_name.remove()
      self.assemble_pythonpath()
      self.show_build_options_and_module_listing()
      self.reset_dispatcher_bookkeeping()
      print("Creating files in build directory: %s" \
        % show_string(abs(self.build_path)))
      self.write_dispatcher_include_template()
      self.write_lib_dispatcher_head()
      self.write_setpath_files()
      self.pickle()
      libtbx.env = self
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
      self.relocate_python_paths_if_necessary()
      python_dispatchers = ["libtbx.python"]
      if (self.is_development_environment() and not self.no_bin_python):
        python_dispatchers.append("python")
      for file_name in python_dispatchers:
        self._write_dispatcher_in_bin(
          source_file=self.python_exe,
          target_file=file_name,
          source_is_python_exe=True)
      for module in self.module_list:
        module.process_command_line_directories()
        # Reload the libtbx_config in case dependencies have changed
        module.process_libtbx_config()

      for path in self.pythonpath:
        sys.path.insert(0, abs(path))

      # Some libtbx_refresh scripts do a `pip install --editable`. The default
      # behavior was changed in setuptools 64 and currently we expect the old
      # behavior. See: https://github.com/pypa/setuptools/pull/3265
      sef_previous = os.environ.get('SETUPTOOLS_ENABLE_FEATURES')
      os.environ['SETUPTOOLS_ENABLE_FEATURES'] = 'legacy_editable'

      for module in self.module_list:
        module.process_libtbx_refresh_py()

      if sef_previous is not None:
          os.environ['SETUPTOOLS_ENABLE_FEATURES'] = sef_previous
      else:
          os.environ.pop('SETUPTOOLS_ENABLE_FEATURES')

      self.write_python_and_show_path_duplicates()
      self.process_exe()
      self.write_command_version_duplicates()
      if (os.name != "nt"):     # LD_LIBRARY_PATH for dependencies
        os.environ[self.ld_library_path_var_name()] = ":".join(
          [abs(p) for p in self.ld_library_path_additions()])
      if sys.version_info[0] == 2:
        if self.build_options.use_conda:
          # refresh loaders.cache for gdk-pixbuf on linux due to gtk2
          if sys.platform.startswith("linux"):
            conda_base = get_conda_prefix()
            command = "{conda_base}/bin/gdk-pixbuf-query-loaders"
            loaders = "{conda_base}/lib/gdk-pixbuf-2.0/2.10.0/loaders/*.so"
            cache = "{conda_base}/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache"
            if os.path.isfile(command.format(conda_base=conda_base)):
              command = command + " " + loaders + " > " + cache
              command = command.format(conda_base=conda_base)
              call(command)
        else:
          regenerate_module_files.run(libtbx.env.under_base('.'), only_if_needed=True)
      self.pickle()
      print("libtbx_refresh_is_completed", file=completed_file_name.open("w"))

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
    self.extra_command_line_locations = []
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
          self.extra_command_line_locations.extend(config.get(
            "extra_command_line_locations", []))
          self.required_for_build.extend(config.get(
            "modules_required_for_build", []))
          self.required_for_use.extend(config.get(
            "modules_required_for_use", []))
          self.optional.extend(config.get(
            "optional_modules", []))
          self.optional.extend(
            set(config.get("optional_modules_only_if_explicit_request", [])) &
            self.env.explicitly_requested_modules)
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
      path = dist_path / "src"
      if (path / self.name).isdir() or (path / "__init__.py").isfile():
        result.append(path)
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
    if (file_name_lower.endswith(".md")): return # ignore markdown files
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
      try:
        with io.open(abs(source_file), encoding='utf-8', errors='ignore') as fh:
          source_text = to_str(fh.read(read_size))
      except IOError:
        raise RuntimeError('Cannot read file: "%s"' % source_file)
      if (check_for_hash_bang and not source_text.startswith("#!")):
        if (not suppress_warning):
          msg = 'WARNING: Ignoring file "%s" due to missing "#!"' % (
            source_file)
          print("*"*len(msg))
          print(msg)
          print("*"*len(msg))
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
            op.join(self.name, "exe", file_name[:-len(ext)]+exe_suffix),
            return_relocatable_path=True)
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
      for sub_dir in ["command_line", os.path.join(self.name, "command_line") ]+ \
        [os.path.join(a,"command_line") for a in getattr(self,"extra_command_line_locations",[])]:
        path = dist_path / sub_dir
        if path.isdir():
          result.append(path)
    return result

  def program_directory_paths(self):
    '''
    Returns the directory for the program templates in a module
    Generally, this "programs" directory is at the same level as the
    "command_line" directory
    '''
    result = []
    for dist_path in self.dist_paths_active():
      for sub_dir in ['programs', os.path.join(self.name, 'programs')] + \
        [os.path.join(a,'programs') for a in getattr(self,"extra_command_line_locations",[])]:
        path = dist_path / sub_dir
        if path.isdir():
          result.append(path)
    return result

  def process_command_line_directories(self):
    for source_dir in self.command_line_directory_paths():
      print("Processing: %s" % show_string(abs(source_dir)))
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

  def process_libtbx_refresh_py(self):
    for dist_path in self.dist_paths_active():
      custom_refresh = dist_path / "libtbx_refresh.py"
      if custom_refresh.isfile():
        print("Processing: %s" % show_string(abs(custom_refresh)))
        global_vars = globals()
        global_vars["__name__"] = dist_path.basename() + ".libtbx_refresh"
        global_vars["self"] = self
        with io.open(abs(custom_refresh), encoding='utf-8', errors='ignore') as fh:
          exec(to_str(fh.read()), global_vars)

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
        enable_openmp_if_possible=default_enable_openmp_if_possible,
        enable_boost_threads=True,
        enable_cuda=default_enable_cuda,
        enable_kokkos=default_enable_kokkos,
        use_conda=default_use_conda,
        opt_resources=default_opt_resources,
        precompile_headers=False,
        use_environment_flags=False,
        force_32bit=False,
        msvc_arch_flag=default_msvc_arch_flag,
        enable_cxx11=default_enable_cxx11,
        cxxstd=default_cxxstd,
        skip_phenix_dispatchers=False):

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
      self.env_ldflags = ""
      flg = os.environ.get("CXXFLAGS")
      if flg is not None:
        self.env_cxxflags = flg
      flg = os.environ.get("CFLAGS")
      if flg is not None:
        self.env_cflags = flg
      flg = os.environ.get("CPPFLAGS")
      if flg is not None:
        self.env_cppflags = flg
      flg = os.environ.get("LDFLAGS")
      if flg is not None:
        self.env_ldflags = flg

  def report(self, f=None):
    if (f is None): f = sys.stdout
    print("Compiler:", self.compiler, file=f)
    print("Build mode:", self.mode, file=f)
    print("Warning level:", self.warning_level, file=f)
    print("Precompiled Headers:", self.precompile_headers, file=f)
    print("Static libraries:", self.static_libraries, file=f)
    print("Static exe:", self.static_exe, file=f)
    print("Scan Boost headers:", self.scan_boost, file=f)
    print("Write full flex_fwd.h files:", self.write_full_flex_fwd_h, file=f)
    print("Build Boost.Python extensions:", \
      self.build_boost_python_extensions, file=f)
    print("Define BOOST_PYTHON_NO_PY_SIGNATURES:", \
      self.boost_python_no_py_signatures, file=f)
    print("Define BOOST_PYTHON_BOOL_INT_STRICT:", \
      self.boost_python_bool_int_strict, file=f)
    print("Enable OpenMP if possible:", self.enable_openmp_if_possible, file=f)
    print("Boost threads enabled:", self.enable_boost_threads, file=f)
    print("Enable CUDA:", self.enable_cuda, file=f)
    print("Enable KOKKOS:", self.enable_kokkos, file=f)
    print("Use conda:", self.use_conda, file=f)
    print("Use opt_resources if available:", self.opt_resources, file=f)
    print("Use environment flags:", self.use_environment_flags, file=f)
    print("Enable C++11:", self.enable_cxx11, file=f)
    if( self.use_environment_flags ):
      print("  CXXFLAGS = ", self.env_cxxflags, file=f)
      print("  CFLAGS = ", self.env_cflags, file=f)
      print("  CPPFLAGS = ", self.env_cppflags, file=f)
      print("  LDFLAGS = ", self.env_ldflags, file=f)

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
        print("libtbx.scons: implicit dependency scan disabled for directory", end=' ')
        print(path)
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
      parser.option(None, "--exclude",
        action="store",
        type="string",
        default=None,
        help="Modules to leave out from configuration")
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
             " CPPFLAGS, LDFLAGS")
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
        % default_build_boost_python_extensions,
      metavar="True|False")
    parser.option(None, "--enable_openmp_if_possible",
      action="store",
      type="bool",
      default=default_enable_openmp_if_possible,
      help="use OpenMP if available and known to work (default: %s)"
        % default_enable_openmp_if_possible,
      metavar="True|False")
    parser.option(None, "--enable_boost_threads",
      action="store",
      type="bool",
      default=default_enable_boost_threads,
      help="make Boost.Threads available")
    parser.option(None, "--enable_cuda",
      action="store_true",
      default=default_enable_cuda,
      help="Use optimized CUDA routines for certain calculations.  Requires at least one NVIDIA GPU with compute capability of 2.0 or higher, and CUDA Toolkit 4.0 or higher (default: %s)"
        % default_enable_cuda)
    parser.option(None, "--enable_kokkos",
      action="store_true",
      default=default_enable_kokkos,
      help="Use optimized KOKKOS routines for certain calculations. Which backend (CUDA/HIP/OpenMP) is used, depends on the system (default: %s)"
        % default_enable_kokkos)
    parser.option(None, "--use_conda",
      action="store_true",
      default=default_use_conda,
      help="Use conda as the source for Python and dependencies (default: %s)"
        % default_use_conda)
    parser.option(None, "--opt_resources",
      action="store",
      type="bool",
      default=default_opt_resources,
      help="use opt_resources if available (default: %s)"
        % default_opt_resources,
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
    parser.option(None, "--no_bin_python",
                  action="store_true",
                  default=False,
                  help="do not create <build directory>/bin/python even in a development "
                       "environment")
    parser.option(None, "--enable_cxx11",
      action="store_true",
      default=default_enable_cxx11,
      help="use C++11 standard")
    parser.option(None, "--cxxstd",
      action="store",
      type="choice",
      default=default_cxxstd,
      choices=['c++11', 'c++14'], # this should just be the argument to the -std flag
      help="Set the C++ standard. This cannot be set along with --enable_cxx11")
    parser.option("--skip_phenix_dispatchers",
      action="store_true",
      default=False,
      help="Skip all dispatchers with 'phenix' in the title")
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
        if (sys.maxsize > 2**31-1):
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

      # check that only enable_cxx11 or cxxstd is set
      if self.command_line.options.enable_cxx11 \
        and self.command_line.options.cxxstd is not None:
        raise RuntimeError('''
Both --enable_cxx11 and --cxxstd have been set. Please only set one of
these options.
      ''')

  def option_repository(self, option, opt, value, parser):
    if (not op.isdir(value)):
      raise RuntimeError(
        "Not a directory: --repository %s" % show_string(value))
    self.repository_paths.append(value)

def set_preferred_sys_prefix_and_sys_executable(build_path):
  if (not hasattr(os.path, "samefile")): return
  dirname, basename = op.split(sys.prefix)
  if (dirname and op.samefile(build_path, dirname)):
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

def raise_if_source_directory_suspected():
  likely_sources = []
  for level1 in os.listdir("."):
    if (not op.isdir(level1)): continue
    if (level1 == "cctbx_project"):
      likely_sources.append(level1)
    else:
      for file_name in [
            "libtbx_config",
            "libtbx_refresh.py",
            "libtbx_SConscript"]:
        level2 = op.join(level1, file_name)
        if (op.isfile(level2)):
          likely_sources.append(level1)
  if likely_sources:
    likely_sources = sorted(set(likely_sources))
    if (len(likely_sources) == 1): t = "this"; s = "y"
    else:                          t = "these"; s = "ies"
    msg = ["Safety guard:",
      "  The current working directory appears to be a source directory",
      "  as determined by the presence of %s sub-director%s:" % (t, s)]
    msg.extend(["    %s" % show_string(d) for d in likely_sources])
    msg.extend([
      "  To resolve this problem:",
      "  - If this command was accidentally run in the wrong directory,",
      "    change to the correct build directory.",
      "  - Remove or rename the director%s listed above." %s])
    raise RuntimeError("\n".join(msg))

def cold_start(args):
  raise_if_source_directory_suspected()
  cwd_was_empty_at_start = True
  for file_name in os.listdir("."):
    if not file_name.startswith("."):
      cwd_was_empty_at_start = False
      break
  default_repositories = []
  r = op.dirname(op.dirname(args[0]))
  b = op.basename(r)
  if b.lower().startswith("cctbx_project"):
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

def unpickle(build_path=None, env_name="libtbx_env"):
  '''
  Function for loading a libtbx_env file. The build_path and env_name
  parameters are for checking additional environment files.

  Parameters
  ----------
    build_path: str
      The directory containing the environment file
    env_name: str
      The filename for the environment file, default is "libtbx_env"

  Returns
  -------
    env: environment object
  '''
  # try default location in build directory
  if build_path is None:
    build_path = os.getenv("LIBTBX_BUILD")
  # try default installed location
  if not build_path:
    build_path = get_installed_path()
  set_preferred_sys_prefix_and_sys_executable(build_path=build_path)
  with open(op.join(build_path, env_name), "rb") as libtbx_env:
    env = pickle.load(libtbx_env)
  if (env.python_version_major_minor != sys.version_info[:2]):
    env.raise_python_version_incompatible()
  if (op.realpath(build_path) != op.realpath(abs(env.build_path))):
    env.build_path.reset(build_path)
  # XXX backward compatibility 2018-12-10
  if not hasattr(env.build_options, "use_conda"):
    env.build_options.use_conda = False
  # XXX backward compatibility 2020-08-21
  # for installed copies of cctbx, the installed environment is not modifiable
  if not hasattr(env, "installed"):
    env.installed = False
  if not hasattr(env, "installed_modules"):
    env.installed_modules = []
  if not hasattr(env, "installed_order"):
    env.installed_order = []
  # XXX backward compatibility 2022-01-19
  if not hasattr(env.build_options, "enable_kokkos"):
    env.build_options.enable_kokkos = False
  # XXX backward compatibility 2022-12-07
  if not hasattr(env.build_options, "cxxstd"):
    env.build_options.cxxstd = None
  # update installed location
  if env.installed:
    sys_prefix = get_conda_prefix()
    if sys.platform == 'win32':
      sys_prefix = op.join(sys_prefix, 'library')
    sys_prefix = absolute_path(sys_prefix)
    for i in range(len(env.repository_paths)):
      env.repository_paths[i]._anchor = sys_prefix
    env.bin_path._anchor = sys_prefix
    env.exe_path._anchor = sys_prefix
    env.include_path._anchor = sys_prefix
    env.lib_path._anchor = sys_prefix
    env.path_utility._anchor = sys_prefix
    if sys.platform == 'win32':
      for module in env.module_list:
        for dist_path in module.dist_paths:
          if dist_path is not None:
            dist_path._anchor = sys_prefix
  return env

def warm_start(args):
  env = unpickle()
  # ---------------------------------------------------------------------------
  def _warm_start(env, args):
    pre_processed_args = pre_process_args(args=args[1:])
    env.process_args(pre_processed_args=pre_processed_args)
    if (pre_processed_args.command_line.options.clear_scons_memory):
      env.clear_scons_memory()
    env.refresh()
  # ---------------------------------------------------------------------------
  # check if the loaded environment is modifiable
  if env.installed:
    # fix classes
    from libtbx.env_config import build_options as _build_options
    from libtbx.env_config import environment as _environment
    from libtbx.env_config import module as _module

    globals()['build_options'] = _build_options
    globals()['environment'] = _environment
    globals()['module'] = _module

    # load an existing, modifiable environment from the current directory
    if op.exists(op.join(os.getcwd(), "libtbx_env")):
      env = unpickle(build_path=os.getcwd())
      _warm_start(env, args)
    # or create new environment in current directory if it does not exist
    # currently can only run configuration in the "build" directory, so
    # "modules" is one level up.
    else:
      repository_paths = ['-r', os.path.join('..', 'modules')]
      cctbx_project = os.path.join('..', 'modules', 'cctbx_project')
      if os.path.isdir(cctbx_project):
        repository_paths += ['-r', cctbx_project]
      cold_start(args + repository_paths + ['--no_bin_python'])
      env = unpickle(build_path=os.getcwd())
      env.pickle()  # need separate call
  # continue normally
  else:
    _warm_start(env, args)

def get_installed_path():
  """
  Returns the default location of an installed environment
  """
  if sys.platform == 'win32':
    installed_path = os.path.join(sys.prefix, 'Library', 'share', 'cctbx')
  else:
    installed_path = os.path.join(get_conda_prefix(), 'share', 'cctbx')
  if not os.path.isdir(installed_path):
    # try libtbx/core
    import sysconfig
    paths = sysconfig.get_paths()
    for key in ['purelib', 'platlib']:
      site_packages = paths[key]
      test_path = os.path.join(site_packages, 'libtbx', 'core', 'share', 'cctbx')
      if os.path.isdir(test_path):
        installed_path = test_path
      break
  return installed_path

def _get_env(build_path, env_name='libtbx_env'):
  current_env = None
  if op.isfile(op.join(build_path, env_name)):
    current_env = unpickle(build_path, env_name)
  return current_env

def get_installed_env():
  """
  Returns the installed environment in the default location or None if
  it does not exist.
  """
  return _get_env(get_installed_path())

def get_local_env(build_dir=None):
  """
  Returns the local environment from build_dir or None
  if it does not exist. By default, build_dir is set to the current
  working directory.

  Parameters
  ----------
    build_dir: str
      The directory with the local environment. Defaults to the current
      working directory

  Returns
  -------
    env: libtbx.env_config.environment or None
      The environment loaded from build_dir or None if it does not
  """
  if build_dir is None:
    build_dir = os.getcwd()
  env = None
  for p in [build_dir] + sys.path:
    e = _get_env(p)
    if e is not None:
      env = e
      break
  return env

def get_boost_library_with_python_version(name, libpath):
  """
  Standard Boost.Python libraries may have the Python version appended
  as a suffix. This function returns the name with the current Python
  version if there is a file with that name in LIBPATH. Otherwise, it
  will return the original name. For example, libboost_python.so may be
  named libboost_python27.so for Python 2.7.

  Parameters
  ----------
  name: str
    The base name for a library (e.g. "boost_python")
  libpath: list
    The paths to search for this library

  Returns
  -------
  name: str
    The input name modified with the current Python version, if available
  """

  version = str(sys.version_info.major) + str(sys.version_info.minor)
  for p in libpath:
    name_version = name + version
    if sys.platform == 'win32':
      full_names = [os.path.join(p, name_version + '.dll'),
                    os.path.join(p, name_version + '.lib')]
    else:
      full_name = os.path.join(p, 'lib' + name_version)
      if sys.platform == 'darwin':
        full_name += '.dylib'
      else:
        full_name += '.so'
      full_names = [full_name]
    for full_name in full_names:
      if os.path.isfile(full_name):
        return name_version
  return name

if (__name__ == "__main__"):
  if (len(sys.argv) == 2 and sys.argv[1] == "__libtbx_refresh__"):
    unpickle().refresh()
  else:
    warm_start(sys.argv)
  print("Done.")
