#! /usr/bin/env python

import sys, os, copy, pickle
from os.path import normpath, join, abspath, dirname, isdir, isfile
norm = normpath

import libtbx.config
import libtbx.path
from libtbx.utils import UserError

class registry:

  def __init__(self):
    self.dict = {}
    self.list = []
    self.missing_for_build = []

  def insert(self, position, key, value):
    self.dict[key] = value
    self.list.insert(position, key)

  def merge(self, other):
    i = len(self.list)
    other_reverse = other.list[:]
    other_reverse.reverse()
    for key in other_reverse:
      if (key in self.dict):
        i = self.list.index(key)
      else:
        self.insert(i, key, other.dict[key])
    other_reverse = other.missing_for_build[:]
    for key in other_reverse:
      if (key in self.missing_for_build):
        i = self.missing_for_build.index(key)
      else:
        self.missing_for_build.insert(i, key)

class package:

  def __init__(self, dist_root, name, must_exist=1):
    self.dist_root = dist_root
    self.name = name
    paths = libtbx.config.resolve_redirection(dist_root=dist_root, name=name)
    self.effective_root = paths.effective_root
    self.dist_path = paths.dist_path
    self.config = None
    self.SConscript_path = None
    self.needs_adaptbx = False
    self.effective_root_adaptbx = None
    self.python_path = None
    if (not isdir(self.dist_path)):
      if (must_exist):
        raise UserError("Not a directory: " + self.dist_path)
      self.dist_path = None
    else:
      ppn = libtbx.config.package_pair(name=self.name)
      ppd = libtbx.config.package_pair(name=self.name, dist_root=dist_root)
      for dist_path_suf in ppd.adaptbx_first():
        path = norm(join(dist_path_suf, "libtbx_config"))
        if (isfile(path)):
          self.config = libtbx.config.read_libtbx_config(path=path)
          if (not self.needs_adaptbx):
            self.needs_adaptbx = (    self.name == ppn.primary
                                  and dist_path_suf is ppd.adaptbx)
          break
      for dist_path_suf,name_suf in ppd.zip_adaptbx_first(ppn):
        if (isfile(norm(join(dist_path_suf, "SConscript")))):
          self.SConscript_path = norm(join(name_suf, "SConscript"))
          if (not self.needs_adaptbx):
            self.needs_adaptbx = (    self.name == ppn.primary
                                  and dist_path_suf is ppd.adaptbx)
          break
      for name_suf in ppn.adaptbx_first():
        for dist_path_suf in ppd.adaptbx_first():
          if (isfile(norm(join(dist_path_suf, name_suf, "__init__.py")))):
            self.python_path = norm(dist_path_suf)
            if (not self.needs_adaptbx):
              self.needs_adaptbx = (    self.name == ppn.primary
                                    and dist_path_suf is ppd.adaptbx)
            break
      if (self.needs_adaptbx):
        self.effective_root_adaptbx = dirname(ppd.adaptbx)
      self.test_scripts = []
      for dist_path_suf in ppd.adaptbx_first():
        self.test_scripts.extend(find_test_scripts(directory=dist_path_suf))
      self._build_dependency_registry()

  def _build_dependency_registry(self):
    self.dependency_registry = registry()
    self.dependency_registry.insert(0, self.name, self)
    if (self.config is not None):
      for package_name in self.config.get("packages_required_for_use", []):
        p = package(self.dist_root, package_name)
        self.dependency_registry.merge(p.dependency_registry)
      for package_name in self.config.get("packages_required_for_build", []):
        p = package(self.dist_root, package_name, must_exist=0)
        if (p.dist_path is None):
          self.dependency_registry.missing_for_build.append(package_name)
        else:
          self.dependency_registry.merge(p.dependency_registry)
      for package_name in self.config.get("optional_packages", []):
        try:
          p = package(self.dist_root, package_name)
        except UserError, e:
          e.reset_tracebacklimit()
        else:
          if (p.dist_path is not None):
            self.dependency_registry.merge(p.dependency_registry)

def insert_normed_path(path_list, addl_path):
  addl_path = norm(addl_path)
  if (not addl_path in path_list):
    i = 0
    if (path_list[:1] == ["."]): i = 1
    path_list.insert(i, addl_path)

def open_info(path, mode="w", info="Creating:"):
  print info, path
  try: return open(path, mode)
  except IOError, e:
    raise UserError(str(e))

class libtbx_env:

  def __init__(self, cwd, libtbx_dist):
    self.LIBTBX_BUILD = norm(abspath(os.getcwd()))
    self._shortpath_bat = None
    self.LIBTBX_BUILD = self.abs_path_clean(self.LIBTBX_BUILD)
    libtbx_dist = self.abs_path_clean(libtbx_dist)
    self.python_version_major_minor = sys.version_info[:2]
    if (os.name == "nt"):
      self.LIBTBX_PYTHON_EXE = self.abs_path_clean(join(
        sys.prefix, "python.exe"))
    else:
      self.LIBTBX_PYTHON_EXE = self.abs_path_clean(join(
        sys.prefix, "bin/python"))
    self.LD_LIBRARY_PATH = [norm(join(self.LIBTBX_BUILD, "libtbx"))]
    self.PYTHONPATH = [norm(join(self.LIBTBX_BUILD, "libtbx")), libtbx_dist]
    self.PATH = [norm(join(self.LIBTBX_BUILD, "libtbx/bin"))]
    self.libtbx_dist_root = dirname(libtbx_dist)
    self.libtbx_python = norm(join(
      self.LIBTBX_BUILD,"libtbx","bin","libtbx.python"))
    self.libtbx_path_utility = norm(join(
      libtbx_dist, "libtbx/command_line/path_utility.py"))
    self.package_list = []
    self.effective_roots = []
    self.dist_paths = {"LIBTBX_DIST": libtbx_dist}
    self.scons_in_dist_root = False
    if (os.path.isdir(join(self.libtbx_dist_root, "scons"))):
      self.scons_in_dist_root = True

  def add_package(self, package, explicit_adaptbx):
    self.package_list.insert(0, package.name)
    if (package.effective_root not in self.effective_roots):
      self.effective_roots.append(package.effective_root)
    if (package.effective_root_adaptbx is not None
        and package.effective_root_adaptbx not in self.effective_roots):
      self.effective_roots.append(package.effective_root_adaptbx)
    self.dist_paths[package.name.upper() + "_DIST"] = self.abs_path_clean(
      package.dist_path)
    if (not explicit_adaptbx):
      ppn = libtbx.config.package_pair(name=package.name)
      if (package.name == ppn.primary):
        ppd = libtbx.config.package_pair(
          name=package.name,
          dist_root=package.dist_root)
        self.dist_paths[ppn.adaptbx.upper() + "_DIST"] = self.abs_path_clean(
          ppd.adaptbx)
    if (package.python_path is not None):
      insert_normed_path(self.PYTHONPATH, package.python_path)

  def items(self):
    return self.__dict__.items() + self.dist_paths.items()

  def pickle_dict(self):
    file_name = norm(join(self.LIBTBX_BUILD, "libtbx_env"))
    f = open_info(file_name)
    pickle.dump(self.__dict__, f)
    f.close()

  def python_api_from_libtbx_build_libtbx(self):
    api_file_name = libtbx.config.python_api_version_file_name(
      self.LIBTBX_BUILD)
    if (isfile(api_file_name)):
      return open(api_file_name).readline().strip()
    return None

  def check_python_api(self):
    api_from_process = libtbx.config.python_api_from_process(must_exist=False)
    if (api_from_process is None): return
    api_from_build = self.python_api_from_libtbx_build_libtbx()
    if (api_from_build is None): return
    if (api_from_process != api_from_build):
      raise UserError(("Incompatible Python API's: current version: %s,"
        + " used to build binaries: %s") % (api_from_process, api_from_build))

  def abs_path_short(self, abs_path):
    if (os.name != "nt"): return abs_path
    if (self._shortpath_bat is None):
      self._shortpath_bat = "shortpath.bat"
      libtbx_build = self.LIBTBX_BUILD
      if (libtbx_build is not None):
        self._shortpath_bat = os.path.join(libtbx_build, self._shortpath_bat)
    assert os.path.exists(self._shortpath_bat)
    return os.popen('call "%s" "%s"' %
      (self._shortpath_bat, abs_path), "r").readline().rstrip()

  def abs_path_clean(self, path):
    abs_path = norm(abspath(path))
    if (os.name != "nt"): return abs_path
    short = self.abs_path_short(abs_path).split(os.sep)
    orig = abs_path.split(os.sep)
    clean = []
    for o,s in zip(orig, short):
      if (o.find(" ") < 0):
        clean.append(o)
      else:
        clean.append(s)
    return os.sep.join(clean)

def emit_dispatcher_include_sh(env):
  file_name = norm(join(env.LIBTBX_BUILD, "dispatcher_include.sh"))
  if (not isfile(file_name)):
    f = open_info(file_name)
    print >> f, \
      "# sh commands to be included in the auto-generated dispatchers"
    f.close()

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

class unix_update_path:

  def __init__(self, env, shell, s, u):
    assert shell in ["sh", "csh"]
    self.shell = shell
    self.env = env
    self.s = s
    self.u = u
    if (self.shell == "sh"):
      self.setenv = "%s="
    else:
      self.setenv = "setenv %s "

  def write(self, prefixes, var_name, val):
    val = os.pathsep.join(val)
    for f,prefix,action in [(self.s, prefixes[0], "prepend"),
                            (self.u, prefixes[1], "delete")]:
      print >> f, '''%s%s"`'%s' '%s' %s %s '%s'`"''' % (
        prefix,
        self.setenv % var_name,
        self.env.LIBTBX_PYTHON_EXE,
        self.env.libtbx_path_utility,
        action,
        var_name,
        val)
      if (f is self.s and self.shell == "sh"):
        print >> f, '%sexport %s' % (prefixes[0], var_name)

class windows_update_path:

  def __init__(self, env, s, u):
    self.env = env
    self.s = s
    self.u = u

  def write(self, prefixes, var_name, val):
    val = os.pathsep.join(val)
    fmt = '''\
%sfor /F "delims=" %%%%i in ('%s "%s" %s %s "%s"') do set %s=%%%%i'''
    for f,prefix,action in [(self.s, prefixes[0], "prepend"),
                            (self.u, prefixes[1], "delete")]:
      print >> f, fmt % (
        prefix,
        self.env.LIBTBX_PYTHON_EXE,
        self.env.libtbx_path_utility,
        action,
        var_name,
        val,
        var_name)
      print >> f, '%sif "%%%s%%" == "E_M_P_T_Y" set %s=' % (
        prefix, var_name, var_name)

def emit_setpaths_sh(env):
  setpaths_path = norm(join(env.LIBTBX_BUILD, "setpaths.sh"))
  unsetpaths_path = norm(join(env.LIBTBX_BUILD, "unsetpaths.sh"))
  s = open_info(setpaths_path)
  u = open_info(unsetpaths_path)
  update_path = unix_update_path(env, "sh", s, u)
  for f in s, u:
    print >> f, '"%s" -V > /dev/null 2>&1' % env.libtbx_python
    print >> f, \
      'if [ $? -ne 0 -o ! -f "%s" ]; then' % env.libtbx_path_utility
    write_incomplete_libtbx_environment(f)
    print >> f, 'else'
  update_path.write(("  ", "  "), "PATH", env.PATH)
  for un in ["", "un"]:
    print >> s, \
      """  alias libtbx.%ssetpaths='source "%s/%ssetpaths.sh"'""" % (
        un, env.LIBTBX_BUILD, un)
  print >> u, '  unalias libtbx.unsetpaths > /dev/null 2>&1'
  print >> s, '  if [ $# -ne 0 ]; then'
  print >> s, \
    '    if [ $# -ne 1 -o \( "$1" != "all" -a "$1" != "debug" \) ]; then'
  print >> s, '      echo "usage: source setpaths.sh [all|debug]"'
  print >> s, '    else'
  print >> s, '      %s="%s"' % ("LIBTBX_BUILD", env.LIBTBX_BUILD)
  print >> s, '      export %s' % "LIBTBX_BUILD"
  print >> u, '  unset %s' % "LIBTBX_BUILD"
  for var_name, values in env.items():
    if (not var_name.endswith("_DIST")): continue
    print >> s, '      %s="%s"' % (var_name, values)
    print >> s, '      export %s' % var_name
    print >> u, '  unset %s' % var_name
  print >> s, '      if [ "$1" == "debug" ]; then'
  update_path.write(("        ", "  "), "PYTHONPATH", env.PYTHONPATH)
  if (sys.platform.startswith("darwin")):
    ld_library_path = "DYLD_LIBRARY_PATH"
  else:
    ld_library_path = "LD_LIBRARY_PATH"
  update_path.write(("        ", "  "), ld_library_path, env.LD_LIBRARY_PATH)
  print >> s, '      fi'
  print >> s, '    fi'
  print >> s, '  fi'
  for f in s, u:
    print >> f, 'fi'
  s.close()
  u.close()

def emit_env_run_sh(env):
  env_run_sh = norm(join(env.LIBTBX_BUILD, "env_run.sh"))
  f = open_info(env_run_sh)
  print >> f, '. "%s" all' % (os.path.join(env.LIBTBX_BUILD, "setpaths.sh"))
  print >> f, '"%s" "%s" "$@"' % (
    "$LIBTBX_BUILD/libtbx/bin/python",
    "$LIBTBX_DIST/libtbx/command_line/env_run.py")
  f.close()
  os.chmod(env_run_sh, 0755)

def emit_setpaths_csh(env, all):
  if (all): s = "_all"
  else:     s = ""
  setpaths_csh_path = norm(join(env.LIBTBX_BUILD, "setpaths%s.csh"%s))
  unsetpaths_csh_path = norm(join(env.LIBTBX_BUILD, "unsetpaths.csh"))
  s = open_info(setpaths_csh_path)
  u = open_info(unsetpaths_csh_path)
  update_path = unix_update_path(env, "csh", s, u)
  for f in s, u:
    print >> f, '"%s" -V >& /dev/null' % env.libtbx_python
    print >> f, \
      'if ($status != 0 || ! -f "%s") then' % env.libtbx_path_utility
    write_incomplete_libtbx_environment(f)
    print >> f, 'else'
  update_path.write(("  ", "  "), "PATH", env.PATH)
  for un in ["", "un"]:
    print >> s, \
      """  alias libtbx.%ssetpaths 'source "%s/%ssetpaths.csh"'""" % (
        un, env.LIBTBX_BUILD, un)
  print >> u, '  unalias libtbx.unsetpaths'
  if (all): c = "#"
  else:     c = ""
  print >> s, '  %sif ($#argv != 0) then' % c
  print >> s, \
    '    %sif ($#argv != 1 || ("$1" != "all" && "$1" != "debug")) then' % c
  print >> s, '    %s  echo "usage: source setpaths.csh [all|debug]"' % c
  print >> s, '    %selse' % c
  print >> s, '      setenv %s "%s"' % ("LIBTBX_BUILD", env.LIBTBX_BUILD)
  print >> u, '  unsetenv %s' % "LIBTBX_BUILD"
  for var_name, values in env.items():
    if (not var_name.endswith("_DIST")): continue
    print >> s, '      setenv %s "%s"' % (var_name, values)
    print >> u, '  unsetenv %s' % var_name
  print >> s, '      if ("$1" == "debug") then'
  update_path.write(("        ", "  "), "PYTHONPATH", env.PYTHONPATH)
  if (sys.platform.startswith("darwin")):
    ld_library_path = "DYLD_LIBRARY_PATH"
  else:
    ld_library_path = "LD_LIBRARY_PATH"
  update_path.write(("        ", "  "), ld_library_path, env.LD_LIBRARY_PATH)
  print >> s, '      endif'
  print >> s, '    %sendif' % c
  print >> s, '  %sendif' % c
  for f in s, u:
    print >> f, 'endif'
  s.close()
  u.close()

def emit_setpaths_bat(env):
  setpaths_path = norm(join(env.LIBTBX_BUILD, "setpaths.bat"))
  unsetpaths_path = norm(join(env.LIBTBX_BUILD, "unsetpaths.bat"))
  s = open_info(setpaths_path)
  u = open_info(unsetpaths_path)
  update_path = windows_update_path(env, s, u)
  for f in s, u:
    print >> f, '@ECHO off'
    print >> f, 'if not exist "%s" goto fatal_error' % env.LIBTBX_PYTHON_EXE
    print >> f, 'if not exist "%s" goto fatal_error' % env.libtbx_path_utility
    print >> f, 'if exist "%s.exe" goto update_path' % env.libtbx_python
    write_incomplete_libtbx_environment(f)
    print >> f, '  goto end_of_script'
    print >> f, ':update_path'
  update_path.write(("  ", "  "), "PATH", env.PATH)
  for un in ["", "un"]:
    print >> s, '  doskey libtbx.%ssetpaths="%s\\%ssetpaths.bat" $*' % (
      un, env.LIBTBX_BUILD, un)
  print >> u, '  doskey libtbx.unsetpaths='
  print >> s, '  if "%1" == "" goto end_of_script'
  print >> s, '  if not "%2" == "" goto show_usage'
  print >> s, '  if "%1" == "all" goto set_all'
  print >> s, '  if "%1" == "debug" goto set_all'
  print >> s, ':show_usage'
  print >> s, '  echo usage: setpaths [all^|debug]'
  print >> s, '  goto end_of_script'
  print >> s, ':set_all'
  print >> s, '  set %s="%s"' % ("LIBTBX_BUILD", env.LIBTBX_BUILD)
  print >> u, '  set %s=' % "LIBTBX_BUILD"
  for var_name, values in env.items():
    if (not var_name.endswith("_DIST")): continue
    print >> s, '  set %s="%s"' % (var_name, values)
    print >> u, '  set %s=' % var_name
  print >> s, '  if not "%1" == "debug" goto end_of_script'
  update_path.write(("    ", "  "), "PYTHONPATH", env.PYTHONPATH)
  update_path.write(("    ", "  "), "PATH", env.LD_LIBRARY_PATH)
  for f in s, u:
    print >> f, ':end_of_script'
  s.close()
  u.close()

class build_options_t:

  def __init__(self):
    self.compiler = "default"
    self.mode = "release"
    self.static_libraries = False
    self.static_exe = False
    self.scan_boost = False

  def report(self):
    print "Compiler:", self.compiler
    print "Build mode:", self.mode
    print "Static libraries:", self.static_libraries
    print "Static exe:", self.static_exe
    print "Scan Boost headers:", self.scan_boost

  def add_to_libtbx_env(self, env):
    env.build_options_compiler = self.compiler
    env.build_options_mode = self.mode
    env.build_options_static_libraries = self.static_libraries
    env.build_options_static_exe = self.static_exe
    env.build_options_scan_boost = self.scan_boost

  def get_from_libtbx_env(self, env):
    if (hasattr(env, "build_options_compiler")):
      self.compiler = env.build_options_compiler
    if (hasattr(env, "build_options_mode")):
      self.mode = env.build_options_mode
    if (hasattr(env, "build_options_static_libraries")):
      self.static_libraries = env.build_options_static_libraries
    if (hasattr(env, "build_options_static_exe")):
      self.static_exe = env.build_options_static_exe
    if (hasattr(env, "build_options_scan_boost")):
      self.scan_boost = env.build_options_scan_boost
    if (self.static_exe):
      self.static_libraries = True

def emit_SConstruct(env, build_options, packages_dict):
  SConstruct_path = norm(join(env.LIBTBX_BUILD, "SConstruct"))
  f = open_info(SConstruct_path)
  print >> f, 'import libtbx.config'
  print >> f, 'libtbx.config.build_options.set('
  print >> f, '  compiler="%s",' % build_options.compiler
  print >> f, '  optimization=%d,' % int(build_options.mode == "release")
  print >> f, '  debug_symbols=%d,' % int(build_options.mode == "debug")
  print >> f, '  static_libraries=%d,' % int(build_options.static_libraries)
  print >> f, '  static_exe=%d,' % int(build_options.static_exe)
  print >> f, '  scan_boost=%d)' % int(build_options.scan_boost)
  print >> f
  print >> f, 'SConsignFile()'
  for effective_root in env.effective_roots:
    print >> f, 'Repository(r"%s")' % effective_root
  print >> f, 'SConscript("libtbx/SConscript")'
  done = {}
  for package_name in env.package_list:
    p = packages_dict[package_name].SConscript_path
    if (p is not None and not done.has_key(p)):
      print >> f, 'SConscript("%s")' % p
      done[p] = 0
  f.close()

def find_test_scripts(
      directory,
      file_names=["run_tests.py", "run_examples.py"]):
  result = []
  for file_name in file_names:
    path = norm(join(directory, file_name))
    if (isfile(path)):
      result.append(path)
  return result

def collect_test_scripts(env, packages):
  result = find_test_scripts(directory=env.dist_paths["LIBTBX_DIST"])
  package_names = list(packages.list)
  package_names.reverse()
  for package_name in package_names:
    for file_name in packages.dict[package_name].test_scripts:
      if (not file_name in result):
        result.append(file_name)
  return result

def emit_run_tests_csh(env, packages):
  test_scripts = collect_test_scripts(env, packages)
  if (len(test_scripts) > 0):
    path = norm(join(env.LIBTBX_BUILD, "run_tests.csh"))
    f = open_info(path)
    print >> f, "#! /bin/csh -f"
    print >> f, "set noglob"
    print >> f, "set verbose"
    for file_name in test_scripts:
      print >> f, "python %s $*" % file_name
    f.close()
    os.chmod(path, 0755)

def show_help(old_env):
  if (old_env):
    cmd = "libtbx.configure"
  else:
    cmd = "[path_to/]python [path_to/]libtbx/configure.py".replace("/", os.sep)
  print
  print   cmd, "[options] package [...]"
  print
  print   "options:"
  print   "  -h, --help                      show this help message and exit"
  if (old_env):
    print "  --only                          disable previously configured"
    print "                                    packages"
  print   "  --build=[quick|release|debug]   select build mode"
  print   "  --compiler=COMPILER             select non-standard compiler"
  print   "  --static_libraries              build all libraries statically"
  print   "  --static_exe                    link all executables statically"
  print   "                                    (implies --static_libraries)"
  print   "  --scan_boost                    enable implicit dependency scan"
  print   "                                    of boost header files"
  print

def run(libtbx_dist, args, old_env=None):
  env = libtbx_env(os.getcwd(), libtbx_dist)
  build_options = build_options_t()
  if (old_env is not None):
    build_options.get_from_libtbx_env(env=old_env)
  remaining_args = []
  option_only = False
  for arg in args:
    if (arg in ["-h", "--help"]):
      show_help(old_env)
      sys.exit(0)
    elif (arg == "--only"):
      option_only = True
    elif (arg.startswith("--build=")):
      build_options.mode = arg.split("=", 1)[1].strip().lower()
      assert build_options.mode in ("quick", "release", "debug")
    elif (arg.startswith("--compiler=")):
      build_options.compiler = arg.split("=", 1)[1].strip().lower()
    elif (arg == "--static_libraries"):
      build_options.static_libraries = True
    elif (arg == "--static_exe"):
      build_options.static_libraries = True
      build_options.static_exe = True
    elif (arg == "--scan_boost"):
      build_options.scan_boost = True
    elif (arg.startswith("--")):
      show_help(old_env)
      raise UserError("Unknown option: " + arg)
    else:
      remaining_args.append(arg)
  env.compiler = build_options.compiler
  args = remaining_args
  if (len(args) > 0 and not option_only and old_env is not None):
    args.extend(old_env.package_list)
  packages = registry()
  for arg in args:
    packages.merge(package(env.libtbx_dist_root, arg).dependency_registry)
  if (len(packages.list) == 0):
    show_help(old_env)
    raise UserError("At least one package must be specified.")
  print "Python:", sys.version.split()[0], sys.executable
  if (len(packages.missing_for_build) == 0):
    build_options.report()
  print "Top-down list of all packages involved:"
  for package_name in packages.list:
    p = packages.dict[package_name]
    pp = libtbx.config.package_pair(name=package_name)
    if (p.needs_adaptbx and pp.adaptbx not in packages.dict):
      explicit_adaptbx = False
      print "  %s+%s" % pp.primary_first()
    else:
      explicit_adaptbx = True
      print " ", package_name
    env.add_package(packages.dict[package_name], explicit_adaptbx)
  if (len(packages.missing_for_build) != 0):
    if (env.scons_in_dist_root):
      print "************************************"
      print "Warning: packages missing for build:"
      for package_name in packages.missing_for_build:
        print " ", package_name
      print "************************************"
    env.check_python_api()
  build_options.add_to_libtbx_env(env=env)
  env.pickle_dict()
  if (hasattr(os, "symlink")):
    emit_dispatcher_include_sh(env)
    emit_setpaths_sh(env)
    emit_env_run_sh(env)
    for all in [False, True]:
      emit_setpaths_csh(env, all)
  else:
    emit_setpaths_bat(env)
  if (len(packages.missing_for_build) == 0):
    emit_SConstruct(env, build_options, packages.dict)
  emit_run_tests_csh(env, packages)
  return env

def cold_start(args):
  try:
    env = run(
      libtbx_dist=norm(dirname(norm(abspath(args[0])))),
      args=args[1:])
  except UserError, e:
    print "Error:", e
  else:
    os.environ["LIBTBX_BUILD"] = env.LIBTBX_BUILD
    from libtbx.command_line import refresh
    refresh.run()

def warm_start(args):
  from libtbx import config
  old_env = config.env()
  try:
    if (old_env.LIBTBX_BUILD != norm(abspath(os.getcwd()))):
      raise UserError(
        "Current working directory must be: " + old_env.LIBTBX_BUILD)
    run(
      libtbx_dist=old_env.dist_path("libtbx"),
      args=args[1:],
      old_env=old_env)
  except UserError, e:
    print "Error:", e
  else:
    from libtbx.command_line import refresh
    refresh.run()

if (__name__ == "__main__"):
  from libtbx.command_line import configure
  configure.warm_start(sys.argv)
