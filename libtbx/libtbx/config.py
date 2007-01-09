from __future__ import generators
import libtbx.path
from libtbx.option_parser import option_parser
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx import easy_run
import shutil
import pickle
from cStringIO import StringIO
import re
import sys, os

def get_hostname():
  try: import socket
  except KeyboardInterrupt: raise
  except: return None
  try: return socket.gethostname()
  except KeyboardInterrupt: raise
  except: return None

def get_ldd_output(target=None):
  if (target is None): target = sys.executable
  return easy_run.go(command="ldd '%s'" % target).stdout_lines

def get_hp_ux_acc_version():
  run_out = easy_run.go(command="aCC -V").stdout_lines
  version = run_out
  if (len(version) > 0):
    version = version[0].strip().split()
  # aCC: HP aC++/ANSI C B3910B A.06.01 [Jan 05 2005]
  if (len(version) >= 6 and version[5].startswith("A.")):
    return version[5]
  raise Sorry(
    "\n  ".join(["Unknown C++ compiler (aCC -V):"] + run_out))

def get_mipspro_version():
  run_out = easy_run.go(command="CC -version").stdout_lines
  version = run_out
  if (len(version) > 0):
    version = version[0].strip().split()
  if (version[:3] == "MIPSpro Compilers: Version ".split()):
    if (version[3].startswith("7.3")):
      return "73"
    elif (version[3].startswith("7.4")):
      return "74"
  raise Sorry(
    "\n  ".join(["Unknown MIPSpro compiler (CC -version):"] + run_out))

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
    "Fatal Error: Incomplete libtbx environment\\!",
    "*******************************************",
    "Please re-run the libtbx/configure.py command."]
  if (os.name != "nt"):
    for line in message: print >> f, '  echo "%s"' % line
    print >> f, '  echo ""'
  else:
    for line in message: print >> f, '  echo %s' % line
    print >> f, '  echo.'

def open_info(path, mode="w", info="   "):
  print info, os.path.basename(path)
  try: return open(path, mode)
  except IOError, e:
    raise Sorry(str(e))

class common_setpaths(object):

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
      for var_name,path in self.env.var_name_and_build_or_dist_path_pairs():
        self.setenv(var_name=var_name, val=path)
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
        print >> f, \
          '  if [ "$%s" = "L_I_B_T_B_X_E_M_P_T_Y" ]; then unset %s; fi' % (
          var_name, var_name)
      else:
        print >> f, '  if ("$%s" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv %s' % (
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
      print >> f, '  if "%%%s%%" == "L_I_B_T_B_X_E_M_P_T_Y" set %s=' % (
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
  if (os.path.isfile(path)):
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
    self.manage_python_version_major_minor()
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
    self.command_line_redirections = {}
    self.reset_module_registry()
    self.scons_dist_path = None
    self.pythonpath = []

  def raise_python_version_incompatible(self, prev_pvmm=None):
    if (prev_pvmm is None):
      prev_pvmm = "%d.%d" % self.python_version_major_minor
    raise Sorry("Python version incompatible with this build:\n"
      + "  Build directory: %s\n" % show_string(self.build_path)
      + "  Python version used initially: %s\n" % prev_pvmm
      + "  Python version in use now:     %d.%d" % sys.version_info[:2])

  def manage_python_version_major_minor(self):
    path = os.path.join(self.build_path, "lib")
    if (not os.path.isdir(path)):
      os.makedirs(path)
    path = os.path.join(path, "PYTHON_VERSION_MAJOR_MINOR")
    pvmm = "%d.%d" % self.python_version_major_minor
    if (not os.path.isfile(path)):
      print >> open(path, "w"), pvmm
    else:
      prev_pvmm = open(path).read().strip()
      if (prev_pvmm != pvmm):
        self.raise_python_version_incompatible(prev_pvmm=prev_pvmm)

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

  def is_ready_for_build(self):
    return (len(self.missing_for_build) == 0)

  def abs_path_short(self, abs_path):
    if (os.name != "nt"): return abs_path
    if (self._shortpath_bat is None):
      self._shortpath_bat = self.under_build("shortpath.bat")
      assert os.path.exists(self._shortpath_bat)
    return easy_run.fully_buffered(
      command='call "%s" "%s"' % (self._shortpath_bat, abs_path)) \
        .raise_if_errors() \
        .stdout_lines[0].rstrip()

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

  def under_dist(self, module_name, path, default=KeyError, test=None):
    if (default is KeyError):
      result = libtbx.path.norm_join(self.module_dist_paths[module_name], path)
    else:
      mdp = self.module_dist_paths.get(module_name)
      if (mdp is None): return default
      result = libtbx.path.norm_join(mdp, path)
    if (test is None or test(result)):
      return result
    return None

  def dist_path(self, module_name, default=KeyError):
    if (default is KeyError):
      return self.module_dist_paths[module_name]
    return self.module_dist_paths.get(module_name, default)

  def dist_paths(self):
    for module in self.module_list:
      for dist_path in module.dist_paths_active():
        yield dist_path

  def var_name_and_build_or_dist_path_pairs(self):
    yield ("LIBTBX_BUILD", self.build_path)
    for module in self.module_list:
      for name,path in module.name_and_dist_path_pairs():
        yield (name.upper()+"_DIST", path)

  def set_os_environ_all_dist(self):
    for module in self.module_list:
      for name,path in module.name_and_dist_path_pairs():
        os.environ[name.upper()+"_DIST"] = path

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
      raise Sorry('Cannot write command_version_suffix file: "%s"' % path)
    print >> f, self.command_version_suffix

  def read_command_version_suffix(self):
    path = self.under_build("command_version_suffix")
    if (not os.path.isfile(path)):
        self.command_version_suffix = None
    else:
      try:
        self.command_version_suffix = open(path).read().strip()
      except IOError:
        raise Sorry('Cannot read command_version_suffix file: "%s"' % path)

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

  def raise_sorry_not_found_in_repositories(self, message):
    if (isinstance(message, str)): message = [message]
    else: message = list(message)
    message.append("  Repository directories searched:")
    for path in self.repository_paths:
      message.append("    %s" % show_string(path))
    raise Sorry("\n".join(message))

  def listdir_in_repositories(self, test=None):
    for path in self.repository_paths:
      for name in os.listdir(path):
        if (test is None or test(os.path.join(path, name))):
          yield path, name

  def match_in_repositories(self,
        relative_path_pattern,
        test=os.path.isdir,
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
        self.raise_sorry_not_found_in_repositories(
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
        raise Sorry("\n".join(message))
    return libtbx.path.norm_join(*all_matches[0])

  def find_in_repositories(self,
        relative_path,
        test=os.path.isdir,
        optional=True):
    for path in self.repository_paths:
      result = self.abs_path_clean(
        libtbx.path.norm_join(path, relative_path))
      if (test is None or test(result)):
        return result
    if (not optional):
      self.raise_sorry_not_found_in_repositories(
        message="Cannot locate: %s" % show_string(relative_path))
    return None

  def find_dist_path(self, module_name, optional=False):
    dist_path = self.command_line_redirections.get(module_name, None)
    if (dist_path is not None): return dist_path
    dist_path = self.find_in_repositories(relative_path=module_name)
    if (dist_path is not None): return dist_path
    trial_module = module(env=self, name=module_name)
    mate_name = trial_module.names[1]
    if (mate_name != module_name):
      if (self.find_in_repositories(relative_path=mate_name) is not None):
        dist_path = self.match_in_repositories(
          relative_path_pattern="%s(?!_%s)" % (
            module_name, trial_module.mate_suffix))
        if (dist_path is not None): return dist_path
    if (not optional):
      self.raise_sorry_not_found_in_repositories(
        message="Module not found: %s" % module_name)
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

  def process_args(self, pre_processed_args):
    command_line = pre_processed_args.command_line
    for path in pre_processed_args.repository_paths:
      self.add_repository(path=path)
    module_names = []
    for module_name in command_line.args:
      if (module_name in self.command_line_redirections):
        del self.command_line_redirections[module_name]
      elif (module_name.count("=") == 1
            and self.find_dist_path(module_name, optional=True) is None):
        module_name, redirection = module_name.split("=")
        dist_path = self.abs_path_clean(os.path.expandvars(redirection))
        if (not os.path.isdir(dist_path)):
          raise Sorry(
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
        static_libraries=command_line.options.static_libraries,
        static_exe=command_line.options.static_exe,
        scan_boost=command_line.options.scan_boost)
      if (command_line.options.command_version_suffix is not None):
        self.command_version_suffix = \
          command_line.options.command_version_suffix
        self.write_command_version_suffix()
    self.reset_module_registry()
    module_names.append("libtbx")
    module_names.reverse()
    for module_name in module_names:
      self.process_module(
        dependent_module=None, module_name=module_name, optional=False)
    self.scons_dist_path = self.find_dist_path("scons", optional=True)
    self.path_utility = self.under_dist(
      "libtbx", "libtbx/command_line/path_utility.py")
    assert os.path.isfile(self.path_utility)

  def dispatcher_precall_commands(self):
    if (self._dispatcher_precall_commands is None):
      lines = []
      if (    self.python_version_major_minor == (2,2)
          and sys.platform == "linux2"
          and os.path.isfile("/etc/redhat-release")):
        try: red_hat_linux_release = open("/etc/redhat-release").readline()
        except KeyboardInterrupt: raise
        except: pass
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
    f = open(self.under_build("dispatcher_include_template.sh"), "w")
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

  def reset_dispatcher_include_cache(self):
    self._dispatcher_include_at_start = []
    self._dispatcher_include_before_command = []
    include_files = []
    for file_name in os.listdir(self.build_path):
      path = self.under_build(file_name)
      if (not os.path.isfile(path)): continue
      if (    file_name.startswith("dispatcher_include")
          and file_name.endswith(".sh")
          and file_name != "dispatcher_include_template.sh"):
        include_files.append(path)
    include_files.sort()
    for path in include_files:
      print "Processing: %s" % show_string(path)
      try: lines = open(path).read().splitlines()
      except IOError, e: raise Sorry(str(e))
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
        buffer.insert(0, "# included from %s" % path)
        highlight_dispatcher_include_lines(buffer)
        target.extend(buffer)

  def dispatcher_include(self, where):
    assert where in ["at_start", "before_command"]
    assert hasattr(self, "_dispatcher_include_at_start")
    if (where == "at_start"):
      return self._dispatcher_include_at_start
    return self._dispatcher_include_before_command

  def write_bin_sh_dispatcher(self, source_file, target_file):
    f = open(target_file, "w")
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
    print >> f, '# To customize this auto-generated script create'
    print >> f, '#'
    print >> f, '#   dispatcher_include*.sh'
    print >> f, '#'
    print >> f, '# files in %s and run' % show_string(self.build_path)
    print >> f, '#'
    print >> f, '#   libtbx.refresh'
    print >> f, '#'
    print >> f, '# to re-generate the dispatchers.'
    print >> f, '#'
    print >> f, '# See also:'
    print >> f, '#   %s' \
      % show_string(self.under_build("dispatcher_include_template.sh"))
    print >> f, '#'
    print >> f, 'unset PYTHONHOME'
    print >> f, 'LIBTBX_BUILD="%s"' % self.build_path
    print >> f, 'export LIBTBX_BUILD'
    if (source_file is not None):
      dispatcher_name = os.path.basename(target_file)
      assert dispatcher_name.find('"') < 0
      print >> f, 'LIBTBX_DISPATCHER_NAME="%s"' % os.path.basename(target_file)
      print >> f, 'export LIBTBX_DISPATCHER_NAME'
    for line in self.dispatcher_include(where="at_start"):
      print >> f, line
    essentials = [("PYTHONPATH", self.pythonpath)]
    essentials.append((ld_library_path_var_name(), [self.lib_path]))
    essentials.append(("PATH", [self.bin_path]))
    for n,v in essentials:
      if (len(v) == 0): continue
      v = ":".join(v)
      if (sys.platform == "irix6" and n == "LD_LIBRARY_PATH"):
        n32 = "LD_LIBRARYN32_PATH"
        print >> f, 'if [ -n "$%s" ]; then' % n32
        print >> f, '  %s="%s:$%s"' % (n32, v, n32)
        print >> f, '  export %s' % n32
        print >> f, 'elif [ -n "$%s" ]; then' % n
      else:
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
    if (source_file is not None):
      for line in source_specific_dispatcher_include(
                    pattern="LIBTBX_PRE_DISPATCHER_INCLUDE_SH",
                    source_file=source_file):
        print >> f, line
    for line in self.dispatcher_include(where="before_command"):
      print >> f, line
    if (source_file is not None):
      for line in source_specific_dispatcher_include(
                    pattern="LIBTBX_POST_DISPATCHER_INCLUDE_SH",
                    source_file=source_file):
        print >> f, line
      start_python = False
      cmd = ""
      if (source_file.lower().endswith(".py")):
        cmd += " '"+self.python_exe+"'"
        if (len(source_specific_dispatcher_include(
                  pattern="LIBTBX_START_PYTHON",
                  source_file=source_file)) > 3):
          start_python = True
      if (not start_python):
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

  def windows_dispatcher(self, command_path, dispatcher_name,
        unique_pattern="0W6I0N6D0O2W8S5_0D0I8S1P4A3T6C4H9E4R7",
        libtbx_build="3L0I2B2T9B4X2_8B5U5I5L2D4",
        libtbx_dispatcher_name="6L6I7B2T3B2X5_6D8I7S0P2A0T5C8H1E8R3_1N0A9M9E",
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
           (libtbx_dispatcher_name, dispatcher_name),
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
    open(target_file, "wb").write(
      self.windows_dispatcher(
        command_path=source_file,
        dispatcher_name=os.path.splitext(os.path.basename(target_file))[0]))

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
    try: action(source_file, target_file_ext)
    except IOError, e: print "  Ignored:", e

  def write_dispatcher_in_bin(self, source_file, target_file):
    self.write_dispatcher(
      source_file=source_file,
      target_file=self.under_build("bin/"+target_file))

  def write_lib_dispatcher_head(self, target_file="dispatcher_head.sh"):
    if (os.name == "nt"): return
    print "   ", target_file
    self.write_bin_sh_dispatcher(
      source_file=None,
      target_file=self.under_build(target_file))

  def write_setpaths_sh(self, suffix):
    setpaths = unix_setpaths(self, "sh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, 'libtbx_pyhome_save="$PYTHONHOME"'
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
    if (self.is_development_environment()):
      print >> s, """  alias cdlibtbxbuild='cd "%s"'""" % self.build_path
      print >> u, '  unalias cdlibtbxbuild > /dev/null 2>&1'
    setpaths.all_and_debug()
    for f in s, u:
      print >> f, 'fi'
      print >> f, 'if [ -n "$libtbx_pyhome_save" ]; then'
      print >> f, '  PYTHONHOME="$libtbx_pyhome_save"'
      print >> f, '  export PYTHONHOME'
      print >> f, 'fi'
      print >> f, 'unset libtbx_pyhome_save'

  def write_setpaths_csh(self, suffix):
    setpaths = unix_setpaths(self, "csh", suffix)
    s, u = setpaths.s, setpaths.u
    for f in s, u:
      print >> f, 'if ($?PYTHONHOME) then'
      print >> f, '  set libtbx_pyhome_save="$PYTHONHOME"'
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
    if (self.is_development_environment()):
      print >> s, """  alias cdlibtbxbuild 'cd "%s"'""" % self.build_path
      print >> u, '  unalias cdlibtbxbuild'
    setpaths.all_and_debug()
    for f in s, u:
      print >> f, 'endif'
      print >> f, 'if ($?libtbx_pyhome_save) then'
      print >> f, '  setenv PYTHONHOME "$libtbx_pyhome_save"'
      print >> f, '  unset libtbx_pyhome_save'
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
      print >> s, '  doskey libtbx.%s="%s\\%s.bat"' % (
        command, self.build_path, command)
    print >> u, '  doskey libtbx.unsetpaths='
    if (self.is_development_environment()):
      print >> s, '  doskey cdlibtbxbuild=cd "%s"' % self.build_path
      print >> u, '  doskey cdlibtbxbuild='
    setpaths.all_and_debug()
    if (suffix == "_debug"):
      print >> s, '  set PYTHONCASEOK=1' # no unset
    for f in s, u:
      print >> f, ':end_of_script'

  def write_SConstruct(self):
    f = open_info(self.under_build("SConstruct"))
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

  def show_module_listing(self):
    print "Top-down list of all modules involved:"
    top_down_module_list = list(self.module_list)
    top_down_module_list.reverse()
    labels = [module.names_for_module_listing()
      for module in top_down_module_list]
    if (len(labels) == 0): return
    fmt = "  %%-%ds  %%s" % max([len(label) for label in labels])
    for label,module in zip(labels,top_down_module_list):
      for dist_path in module.dist_paths_active():
        print fmt % (label, dist_path)
        label = ""

  def show_build_options_and_module_listing(self):
    print 'Python: %s "%s"' % (sys.version.split()[0], sys.executable)
    if (self.is_ready_for_build()):
      self.build_options.report()
    print "command_version_suffix:", self.command_version_suffix
    self.show_module_listing()
    if (len(self.missing_for_use) > 0):
      raise Sorry("Missing modules: "
        + " ".join(self.missing_for_use.keys()))
    if (not self.is_ready_for_build()):
      if (self.scons_dist_path is not None):
        print "***********************************"
        print "Warning: modules missing for build:"
        for module_name in self.missing_for_build.keys():
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
    for path in [self.exe_path, self.under_build("exe_dev")]:
      if (os.path.isdir(path)):
        print 'Processing: "%s"' % path
        for file_name in os.listdir(path):
          if (file_name[0] == "."): continue
          self.write_dispatcher_in_bin(
            source_file=libtbx.path.norm_join(path, file_name),
            target_file=file_name)

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
      try: shutil.copy(source_file, target_file)
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

  def is_development_environment(self):
    libtbx_cvs_root = self.under_dist("libtbx", "CVS/Root")
    if (os.path.isfile(libtbx_cvs_root)):
      try: libtbx_cvs_root = open(libtbx_cvs_root).read()
      except IOError: pass
      else:
        if (libtbx_cvs_root.lower().find("ccp") >= 0): return False
    for module in self.module_list:
      if (module.has_top_level_cvs_directory()):
        return True
    return False

  def refresh(self):
    self.assemble_pythonpath()
    self.show_build_options_and_module_listing()
    self.reset_dispatcher_include_cache()
    print 'Creating files in build directory:\n  "%s"' % self.build_path
    self.write_dispatcher_include_template()
    self.write_lib_dispatcher_head()
    self.write_setpath_files()
    self.pickle()
    if (self.is_ready_for_build()):
      self.write_SConstruct()
    if (os.name != "nt"):
      self.write_run_tests_csh()
    self.clear_bin_directory()
    if (not os.path.isdir(self.bin_path)):
      os.makedirs(self.bin_path)
    python_dispatchers = ["libtbx.python"]
    if (self.is_development_environment()):
      python_dispatchers.append("python")
    for file_name in python_dispatchers:
      self.write_dispatcher_in_bin(
        source_file=self.python_exe,
        target_file=file_name)
    for module in self.module_list:
      module.process_command_line_directories()
    os.environ["LIBTBX_BUILD"] = self.build_path # to support libtbx.load_env
    for path in self.pythonpath:
      sys.path.insert(0, path)
    for module in self.module_list:
      module.process_libtbx_refresh_py()
    self.write_python_and_show_path_duplicates()
    self.process_exe()
    self.write_command_version_duplicates()

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
    if (os.path.normcase(name).endswith(os.path.normcase(mate_suffix))):
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
    dist_paths = []
    for dist_path in self.dist_paths:
      if (dist_path is not None):
        while True:
          path = libtbx.path.norm_join(dist_path, "libtbx_config")
          if (not os.path.isfile(path)):
            config = None
            break
          try: f = open(path)
          except IOError: raise Sorry(
            'Cannot open configuration file: "%s"' % path)
          try: config = eval(" ".join(f.readlines()), {}, {})
          except KeyboardInterrupt: raise
          except: raise Sorry('Corrupt configuration file: "%s"' % path)
          f.close()
          redirection = config.get("redirection", None)
          if (redirection is None):
            break
          if (not isinstance(redirection, str)):
            raise Sorry(
              'Corrupt configuration file:\n'
              '  file = "%s"\n'
              '  redirection must be a Python string' % path)
          new_dist_path = os.path.expandvars(redirection)
          if (not os.path.isabs(new_dist_path)):
            new_dist_path = libtbx.path.norm_join(dist_path, new_dist_path)
          new_dist_path = self.env.abs_path_clean(new_dist_path)
          if (not os.path.isdir(new_dist_path)):
            raise Sorry(
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

  def has_top_level_cvs_directory(self):
    for dist_path in self.dist_paths_active():
      if (os.path.isdir(os.path.join(dist_path, "CVS"))):
        return True
    return False

  def write_dispatcher(self,
        source_dir,
        file_name,
        suppress_warning,
        target_file_name_infix=""):
    source_file = libtbx.path.norm_join(source_dir, file_name)
    if (not os.path.isfile(source_file)): return
    file_name_lower = file_name.lower()
    if (file_name_lower.startswith("__init__.py")): return
    if (file_name_lower.endswith(".pyc")): return
    if (file_name_lower.endswith(".pyo")): return
    if (file_name[0] == "."): return
    ext = os.path.splitext(file_name_lower)[1]
    if (os.name == "nt"):
      if (ext not in windows_pathext): return
    elif (ext != ".sh" and ext != ".py"):
      try: hash_bang = open(source_file).read(2)
      except IOError:
        raise Sorry('Cannot read file: "%s"' % source_file)
      if (hash_bang != "#!"):
        if (ext != ".bat" and not suppress_warning):
          msg = 'WARNING: Ignoring file "%s" due to missing "#!"' % (
            source_file)
          print "*"*len(msg)
          print msg
          print "*"*len(msg)
        return
    if (file_name == "ipython.py" and self.name == "libtbx"):
      try: import IPython
      except ImportError: return
    target_file = self.name + target_file_name_infix
    if (not file_name_lower.startswith("main.")
         or file_name_lower.count(".") != 1):
      target_file += "." + os.path.splitext(file_name)[0]
    self.env.write_dispatcher_in_bin(
      source_file=source_file,
      target_file=target_file)

  def process_command_line_directories(self):
    for dist_path in self.dist_paths_active():
      for source_dir,suppress_warning in ([
            ("command_line", False),
            (self.name+"/command_line", False)]):
        source_dir = libtbx.path.norm_join(dist_path, source_dir)
        if (not os.path.isdir(source_dir)): continue
        print 'Processing: "%s"' % source_dir
        for file_name in os.listdir(source_dir):
          self.write_dispatcher(
            source_dir=source_dir,
            file_name=file_name,
            suppress_warning=suppress_warning)

  def process_python_command_line_scripts(self,
        source_dir,
        print_prefix="  ",
        target_file_name_infix=""):
    print print_prefix+'Processing: %s' % show_string(source_dir)
    for file_name in os.listdir(source_dir):
      if (not file_name.endswith(".py")): continue
      self.write_dispatcher(
        source_dir=source_dir,
        file_name=file_name,
        suppress_warning=False,
        target_file_name_infix=target_file_name_infix)

  def process_libtbx_refresh_py(self):
    for dist_path in self.dist_paths_active():
      custom_refresh = libtbx.path.norm_join(dist_path, "libtbx_refresh.py")
      if (os.path.isfile(custom_refresh)):
        print 'Processing: "%s"' % custom_refresh
        execfile(custom_refresh, {}, {"self": self})

  def collect_test_scripts(self,
        file_names=["run_tests.py", "run_examples.py"]):
    result = []
    for dist_path in self.dist_paths_active():
      for file_name in file_names:
        path = libtbx.path.norm_join(dist_path, file_name)
        if (os.path.isfile(path)): result.append(path)
    return result

class build_options:

  def __init__(self,
        compiler,
        mode,
        warning_level,
        static_libraries,
        static_exe,
        scan_boost):
    adopt_init_args(self, locals())
    assert self.mode in ["release", "quick", "debug", "debug_optimized"]
    assert self.warning_level >= 0
    self.optimization = (self.mode in ["release", "debug_optimized"])
    self.debug_symbols = (self.mode in ["debug", "debug_optimized"])
    if (self.static_exe):
      self.static_libraries = True

  def report(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Compiler:", self.compiler
    print >> f, "Build mode:", self.mode
    print >> f, "Warning level:", self.warning_level
    print >> f, "Static libraries:", self.static_libraries
    print >> f, "Static exe:", self.static_exe
    print >> f, "Scan Boost headers:", self.scan_boost

class include_registry:

  def __init__(self):
    self._boost_dir_name = "boost"
    self._had_message = {}
    self.scan_boost()

  def scan_boost(self, flag=False):
    self._scan_boost = flag
    return self

  def set_boost_dir_name(self, path):
    self._boost_dir_name = os.path.basename(path).lower()

  def scan_flag(self, path):
    if (not self._scan_boost
        and os.path.basename(path).lower() == self._boost_dir_name):
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

class pre_process_args:

  def __init__(self, args, default_repository=None):
    self.repository_paths = []
    if (len(args) == 0): args = ["--help"]
    if (default_repository is None):
      command_name = "libtbx.configure"
      self.warm_start = True
    else:
      command_name = "libtbx/configure.py"
      self.warm_start = False
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
        help="path to source code repository",
        metavar="DIRECTORY")
      if (hasattr(os.path, "samefile")):
        parser.option(None, "--current_working_directory",
          action="store",
          type="string",
          default=None,
          help="preferred spelling of current working directory",
          metavar="DIRECTORY")
      parser.option(None, "--build",
        choices=("release", "quick", "debug", "debug_optimized"),
        default="release",
        help="build mode",
        metavar="release|quick|debug|debug_optimized")
      parser.option(None, "--compiler",
        action="store",
        type="string",
        default="default",
        help="select non-standard compiler",
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
      parser.option(None, "--command_version_suffix",
        action="store",
        type="string",
        default=None,
        help="version suffix for commands in bin directory",
        metavar="STRING")
    self.command_line = parser.process(args=args)
    if (not hasattr(os.path, "samefile")):
      self.command_line.options.current_working_directory = None
    if (default_repository is not None):
      self.repository_paths.append(default_repository)

  def option_repository(self, option, opt, value, parser):
    if (not os.path.isdir(value)):
      raise Sorry("Not a directory: --repository %s" % show_string(value))
    self.repository_paths.append(value)

def set_preferred_sys_prefix_and_sys_executable(build_path):
  if (not hasattr(os.path, "samefile")): return
  dirname, basename = os.path.split(sys.prefix)
  if (os.path.samefile(build_path, dirname)):
    new_prefix = os.path.join(build_path, basename)
    if (os.path.samefile(new_prefix, sys.prefix)):
      p = os.path.normpath(os.path.normcase(sys.prefix))
      if (sys.prefix != new_prefix):
        sys.prefix = new_prefix
      e = os.path.normpath(os.path.normcase(sys.executable))
      if (e.startswith(p)):
        new_executable = sys.prefix + e[len(p):]
        if (os.path.samefile(new_executable, sys.executable)):
          if (sys.executable != new_executable):
            sys.executable = new_executable

def cold_start(args):
  pre_processed_args = pre_process_args(
    args=args[1:],
    default_repository=os.path.dirname(os.path.dirname(args[0])))
  build_path=pre_processed_args.command_line.options.current_working_directory
  if (build_path is None):
    build_path = os.getcwd()
  else:
    if (not os.path.isabs(build_path)):
      raise Sorry("Not an absolute path name:"
        " --current_working_directory %s" % show_string(build_path))
    if (not os.path.isdir(build_path)):
      raise Sorry("Not a directory:"
        " --current_working_directory %s" % show_string(build_path))
    if (not os.path.samefile(build_path, os.getcwd())):
      raise Sorry("Not equivalent to the current working directory:"
        " --current_working_directory %s" % show_string(build_path))
    n = len(os.sep)
    while (len(build_path) > n and build_path.endswith(os.sep)):
      build_path = build_path[:-n]
  set_preferred_sys_prefix_and_sys_executable(build_path=build_path)
  env = environment(build_path=build_path)
  env.process_args(pre_processed_args=pre_processed_args)
  env.refresh()

def unpickle():
  build_path = os.environ["LIBTBX_BUILD"]
  set_preferred_sys_prefix_and_sys_executable(build_path=build_path)
  libtbx_env = open(os.path.join(build_path, "libtbx_env"), "rb")
  env = pickle.load(libtbx_env)
  if (env.python_version_major_minor != sys.version_info[:2]):
    env.raise_python_version_incompatible()
  return env

def warm_start(args):
  pre_processed_args = pre_process_args(args=args[1:])
  env = unpickle()
  env.process_args(pre_processed_args=pre_processed_args)
  env.refresh()
