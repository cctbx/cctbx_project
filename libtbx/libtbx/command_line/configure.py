#! /usr/bin/env python

import sys, os, pickle
from os.path import normpath, join, abspath, dirname, isdir, isfile
norm = normpath

import libtbx.config
from libtbx.config import UserError

class registry:

  def __init__(self):
    self.dict = {}
    self.list = []
    self.missing_for_build = []

  def append(self, key, value):
    self.dict[key] = value
    self.list.append(key)

  def insert(self, position, key, value):
    self.dict[key] = value
    self.list.insert(position, key)

  def merge(self, other):
    i = 0
    for key in other.list:
      if (not key in self.dict):
        self.insert(i, key, other.dict[key])
        i += 1
    i = 0
    for key in other.missing_for_build:
      if (not key in self.missing_for_build):
        self.missing_for_build.insert(i, key)
        i += 1

class package:

  def __init__(self, dist_root, name, must_exist=1):
    self.dist_root = dist_root
    self.name = name
    self.dist_path = norm(join(dist_root, name))
    self.config = None
    self.SConscript_path = None
    self.needs_adaptbx = 00000
    self.python_path = None
    if (not isdir(self.dist_path)):
      if (must_exist):
        raise UserError("Not a directory: " + self.dist_path)
      self.dist_path = None
    else:
      ppn = libtbx.config.package_pair(self.name)
      ppd = libtbx.config.package_pair(self.dist_path)
      for dist_path_suf in ppd.adaptbx_first():
        path = norm(join(dist_path_suf, "libtbx_config"))
        if (isfile(path)):
          try:
            f = open(path)
          except:
            raise UserError("Cannot open configuration file: " + path)
          try:
            self.config = eval(" ".join(f.readlines()))
          except:
            raise UserError("Corrupt file: " + path)
          f.close()
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
      self._build_dependency_registry()

  def _build_dependency_registry(self):
    self.dependency_registry = registry()
    self._resolve_dependencies(self.dependency_registry)

  def _resolve_dependencies(self, registry):
    if (self.name in registry.dict):
      raise UserError("Dependency cycle detected: "
        + str(registry.list) + " + " + self.name)
    registry.append(self.name, self)
    if (self.config is not None):
      try:
        required_packages = self.config["packages_required_for_use"]
      except:
        pass
      else:
        for required_package_name in required_packages:
          package(self.dist_root, required_package_name)._resolve_dependencies(
            registry)
      try:
        required_packages = self.config["packages_required_for_build"]
      except:
        pass
      else:
        for required_package_name in required_packages:
          p = package(self.dist_root, required_package_name, must_exist=0)
          if (p.dist_path is None):
            registry.missing_for_build.append(required_package_name)
          else:
            p._resolve_dependencies(registry)

def insert_normed_path(path_list, addl_path):
  addl_path = norm(addl_path)
  if (not addl_path in path_list):
    i = 0
    if (path_list[:1] == ["."]): i = 1
    path_list.insert(i, addl_path)

def open_info(path, mode="w", info="Creating:"):
  print info, path
  return open(path, mode)

class libtbx_env:

  def __init__(self, cwd, libtbx_dist):
    self.python_version_major_minor = sys.version_info[:2]
    if (os.name == "nt"):
      self.LIBTBX_PYTHON_EXE = norm(abspath(join(sys.prefix, "python.exe")))
    else:
      self.LIBTBX_PYTHON_EXE = norm(abspath(join(sys.prefix, "bin", "python")))
    self.LIBTBX_BUILD = norm(abspath(os.getcwd()))
    self.LIBTBX_DIST_ROOT = norm(dirname(libtbx_dist))
    self.LD_LIBRARY_PATH = [norm(join(self.LIBTBX_BUILD, "libtbx"))]
    self.PYTHONPATH = [norm(join(self.LIBTBX_BUILD, "libtbx")), libtbx_dist]
    self.PATH = [norm(join(self.LIBTBX_BUILD, "libtbx/bin"))]
    self.package_list = []
    self.dist_paths = {"LIBTBX_DIST": libtbx_dist}
    self.scons_in_dist_root = 00000
    if (os.path.isdir(join(self.LIBTBX_DIST_ROOT, "scons"))):
      self.scons_in_dist_root = 0001

  def add_package(self, package, explicit_adaptbx):
    self.package_list.insert(0, package.name)
    self.dist_paths[package.name.upper() + "_DIST"] = package.dist_path
    if (not explicit_adaptbx):
      pp = libtbx.config.package_pair(package.name)
      if (package.name == pp.primary):
        self.dist_paths[pp.adaptbx.upper() + "_DIST"] = \
          package.dist_path + libtbx.config.adaptor_toolbox_suffix
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
    if (os.path.isfile(api_file_name)):
      return open(api_file_name).readline().strip()
    return None

  def check_python_api(self):
    api_from_include = libtbx.config.python_api_from_include(must_exist=00001)
    if (api_from_include is None): return
    api_from_build = self.python_api_from_libtbx_build_libtbx()
    if (api_from_build is None): return
    if (api_from_include != api_from_build):
      raise UserError(("Incompatible Python API's: current version: %s,"
        + " used to build binaries: %s") % (api_from_include, api_from_build))

def emit_setpaths_sh(env):
  for file_name in ("setpaths.sh", "env_run.sh"):
    full_path = norm(join(env.LIBTBX_BUILD, file_name))
    f = open_info(full_path)
    for var_name, values in env.items():
      if (var_name.upper() != var_name): continue
      if (var_name == "LD_LIBRARY_PATH" and sys.platform.startswith("darwin")):
        var_name = "DYLD_LIBRARY_PATH"
      if (type(values) == type([])):
        val = os.pathsep.join(values)
        print >> f, 'if [ ! -n "$%s" ]; then' % (var_name,)
        print >> f, '  %s=""' % (var_name,)
        print >> f, 'fi'
        print >> f, '%s="%s%s$%s"' % (var_name, val, os.pathsep, var_name)
      else:
        print >> f, '%s="%s"' % (var_name, values)
      print >> f, 'export %s' % (var_name,)
    if (file_name == "env_run.sh"):
      print >> f, '"%s" "%s" "$@"' % (
        "$LIBTBX_BUILD/libtbx/bin/libtbx.python",
        "$LIBTBX_DIST/libtbx/command_line/env_run.py")
    f.close()
    os.chmod(full_path, 0755)

def emit_setpaths_csh(env):
  libtbx_python = norm(join(
    env.LIBTBX_BUILD, "libtbx","bin","libtbx.python"))
  libtbx_path_utility = norm(join(
    env.LIBTBX_DIST_ROOT, "libtbx","libtbx","command_line","path_utility.py"))
  setpaths_csh_path = norm(join(env.LIBTBX_BUILD, "setpaths.csh"))
  unsetpaths_csh_path = norm(join(env.LIBTBX_BUILD, "unsetpaths.csh"))
  s = open_info(setpaths_csh_path)
  u = open_info(unsetpaths_csh_path)
  for f in s, u:
    print >> f, '%s -V >& /dev/null' % libtbx_python
    print >> f, 'if ($status != 0 || ! -f "%s") then' % libtbx_path_utility
    print >> f, '  echo "*******************************************"'
    print >> f, '  echo "Fatal Error: Incomplete libtbx environment!"'
    print >> f, '  echo "*******************************************"'
    print >> f, '  echo "Please re-run the libtbx/configure.py command."'
    print >> f, '  echo ""'
    print >> f, 'else'
  for var_name, values in env.items():
    if (var_name.upper() != var_name): continue
    if (var_name == "LD_LIBRARY_PATH" and sys.platform.startswith("darwin")):
      var_name = "DYLD_LIBRARY_PATH"
    if (type(values) == type([])):
      val = os.pathsep.join(values)
      fmt_args = (var_name, libtbx_python, libtbx_path_utility, var_name, val)
      print >> s, '''  setenv %s "`%s %s prepend %s '%s'`"''' % fmt_args
      print >> u, '''  setenv %s "`%s %s delete %s '%s'`"''' % fmt_args
    else:
      print >> s, '  setenv %s "%s"' % (var_name, values)
      print >> u, '  unsetenv %s' % var_name
  for f in s, u:
    print >> f, 'endif'
  s.close()
  u.close()

def join_path_ld_library_path(env):
  joined_path = env.PATH
  for path in env.LD_LIBRARY_PATH:
    if (not path in joined_path):
      joined_path.append(path)
  return joined_path

def emit_setpaths_bat(env):
  setpaths_bat_path = norm(join(env.LIBTBX_BUILD, "setpaths.bat"))
  f = open_info(setpaths_bat_path)
  print >> f, '@ECHO off'
  for var_name, values in env.items():
    if (var_name.upper() != var_name): continue
    if (type(values) == type([])):
      if (var_name == "LD_LIBRARY_PATH"): continue
      if (var_name == "PATH"):
        values = join_path_ld_library_path(env)
      val = os.pathsep.join(values)
      print >> f, 'if not defined %s set %s=' % (var_name, var_name)
      print >> f, 'set %s=%s%s%%%s%%' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'set %s=%s' % (var_name, values)
  print >> f, 'if not defined PATHEXT set PATHEXT='
  print >> f, 'set PATHEXT=.PX;.PY;%PATHEXT%'
  print >> f, \
   'call "%LIBTBX_PYTHON_EXE%" "%LIBTBX_DIST%\\libtbx\\assoc_ftype.py"'
  f.close()

class build_options_t:

  def __init__(self):
    self.compiler = "default"
    self.mode = "release"
    self.static_libraries = 00000
    self.static_exe = 00000
    self.scan_boost = 00000

  def report(self):
    print "Compiler:", self.compiler
    print "Build mode:", self.mode
    print "Static libraries:", self.static_libraries
    print "Static exe:", self.static_exe
    print "Scan Boost headers:", self.scan_boost

def emit_SConstruct(env, build_options, packages_dict):
  SConstruct_path = norm(join(env.LIBTBX_BUILD, "SConstruct"))
  f = open_info(SConstruct_path)
  print >> f, 'import libtbx.config'
  print >> f, 'import libtbx.dblite'
  print >> f, 'libtbx.config.build_options.set('
  print >> f, '  compiler="%s",' % build_options.compiler
  print >> f, '  optimization=%d,' % int(build_options.mode == "release")
  print >> f, '  debug_symbols=%d,' % int(build_options.mode == "debug")
  print >> f, '  static_libraries=%d,' % int(build_options.static_libraries)
  print >> f, '  static_exe=%d,' % int(build_options.static_exe)
  print >> f, '  scan_boost=%d)' % int(build_options.scan_boost)
  print >> f, 'if ("SConsignFile" in dir()):'
  print >> f, '  SConsignFile(name=".sconsign", dbm_module=libtbx.dblite)'
  print >> f
  print >> f, 'Repository(r"%s")' % (env.LIBTBX_DIST_ROOT,)
  print >> f, 'SConscript("libtbx/SConscript")'
  done = {}
  for package_name in env.package_list:
    p = packages_dict[package_name].SConscript_path
    if (p is not None and not done.has_key(p)):
      print >> f, 'SConscript("%s")' % p
      done[p] = 0
  f.close()

def run(libtbx_dist, args):
  env = libtbx_env(os.getcwd(), libtbx_dist)
  build_options = build_options_t()
  remaining_args = []
  for arg in args:
    if (arg.startswith("--compiler=")):
      build_options.compiler = arg.split("=", 1)[1].strip().lower()
    elif (arg.startswith("--build=")):
      build_options.mode = arg.split("=", 1)[1].strip().lower()
      assert build_options.mode in ("quick", "release", "debug")
    elif (arg == "--static_libraries"):
      build_options.static_libraries = 0001
    elif (arg == "--static_exe"):
      build_options.static_libraries = 0001
      build_options.static_exe = 0001
    elif (arg == "--scan_boost"):
      build_options.scan_boost = 0001
    elif (arg.startswith("--")):
      raise UserError("Unknown option: " + arg)
    else:
      remaining_args.append(arg)
  env.compiler = build_options.compiler
  args = remaining_args
  packages = registry()
  for arg in args:
    packages.merge(package(env.LIBTBX_DIST_ROOT, arg).dependency_registry)
  if (len(packages.list) == 0):
    raise UserError("At least one package must be specified.")
  if (len(packages.missing_for_build) == 0):
    build_options.report()
  print "Top-down list of all packages involved:"
  for package_name in packages.list:
    p = packages.dict[package_name]
    pp = libtbx.config.package_pair(package_name)
    if (p.needs_adaptbx and pp.adaptbx not in packages.dict):
      explicit_adaptbx = 00000
      print "  %s+%s" % pp.primary_first()
    else:
      explicit_adaptbx = 0001
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
  env.pickle_dict()
  if (hasattr(os, "symlink")):
    emit_setpaths_sh(env)
    emit_setpaths_csh(env)
  else:
    emit_setpaths_bat(env)
  if (len(packages.missing_for_build) == 0):
    emit_SConstruct(env, build_options, packages.dict)
  return env

def cold_start(args):
  try:
    env = run(libtbx_dist=norm(dirname(norm(abspath(args[0])))), args=args[1:])
  except UserError, e:
    print "Error:", e
  else:
    sys.path.insert(0, norm(join(env.dist_paths["LIBTBX_DIST"], "libtbx")))
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
    run(libtbx_dist=old_env.dist_paths["LIBTBX_DIST"], args=args[1:])
  except UserError, e:
    print "Error:", e
  else:
    from libtbx.command_line import refresh
    refresh.run()

if (__name__ == "__main__"):
  from libtbx.command_line import configure
  configure.warm_start(sys.argv)
