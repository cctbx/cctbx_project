import sys, os, pickle
from os.path import normpath, join, abspath, dirname, isdir
norm = normpath

class UserError(Exception): pass

class registry:

  def __init__(self):
    self.dict = {}
    self.list = []

  def append(self, key, value):
    self.dict[key] = value
    self.list.append(key)

  def insert(self, position, key, value):
    self.dict[key] = value
    self.list.insert(position, key)

  def merge(self, other):
    i = 0
    for name in other.list:
      if (not name in self.dict):
        self.insert(i, name, other.dict[name])
        i += 1

class package:

  def __init__(self, dist_root, name):
    self.dist_root = dist_root
    self.name = name
    self.dist_path = norm(join(dist_root, name))
    if (not isdir(self.dist_path)):
      raise UserError("Not a package directory: " + self.dist_path)
    config_path = norm(join(self.dist_path, "libtbx_config"))
    try:
      f = open(config_path)
    except:
      self.config = None
    else:
      try:
        self.config = eval(" ".join(f.readlines()))
      except:
        raise UserError("Corrupt file: " + config_path)
    self._build_dependency_registry()

  def _build_dependency_registry(self):
    self.dependency_registry = registry()
    self._resolve_dependencies(self.dependency_registry)

  def _resolve_dependencies(self, registry):
    if (self.name in registry.dict):
      raise UserError("Dependency cycle detected: "
        + str(registry.list) + " + " + self.name)
    registry.append(self.name, self)
    if (self.config != None):
      for required_package_name in self.config["required_packages"]:
        package(self.dist_root, required_package_name)._resolve_dependencies(
          registry)

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
    self.LIBTBX_PYTHON_EXE = norm(abspath(sys.executable))
    self.LIBTBX_BUILD = norm(abspath(os.getcwd()))
    self.LIBTBX_DIST_ROOT = norm(dirname(libtbx_dist))
    self.LD_LIBRARY_PATH = [norm(join(self.LIBTBX_BUILD, "libtbx"))]
    self.PYTHONPATH = [norm(join(self.LIBTBX_BUILD, "libtbx")), libtbx_dist]
    self.PATH = [norm(join(self.LIBTBX_BUILD, "libtbx/bin"))]
    self.package_list = []
    self.dist_paths = {"LIBTBX_DIST": libtbx_dist}

  def add_package(self, package):
    self.package_list.insert(0, package.name)
    self.dist_paths[package.name.upper() + "_DIST"] = package.dist_path
    insert_normed_path(self.PYTHONPATH, package.dist_path)

  def items(self):
    return self.__dict__.items() + self.dist_paths.items()

  def pickle_dict(self):
    file_name = norm(join(self.LIBTBX_BUILD, "libtbx_env"))
    f = open_info(file_name)
    pickle.dump(self.__dict__, f)
    f.close()

def emit_env_run_sh(env):
  env_run_sh_path = norm(join(env.LIBTBX_BUILD, "env_run.sh"))
  f = open_info(env_run_sh_path)
  for var_name, values in env.items():
    if (var_name.upper() != var_name): continue
    if (type(values) == type([])):
      val = os.pathsep.join(values)
      print >> f, 'if [ ! -n "$%s" ]; then' % (var_name,)
      print >> f, '  %s=""' % (var_name,)
      print >> f, 'fi'
      print >> f, '%s="%s%s$%s"' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, '%s="%s"' % (var_name, values)
    print >> f, 'export %s' % (var_name,)
  print >> f, '$LIBTBX_PYTHON_EXE "$LIBTBX_DIST/libtbx/command_line/env_run.py" $*'
  f.close()
  os.chmod(env_run_sh_path, 0755)

def emit_setpaths_csh(env):
  setpaths_csh_path = norm(join(env.LIBTBX_BUILD, "setpaths.csh"))
  f = open_info(setpaths_csh_path)
  for var_name, values in env.items():
    if (var_name.upper() != var_name): continue
    if (var_name == "LD_LIBRARY_PATH" and sys.platform.startswith("darwin")):
      var_name = "DYLD_LIBRARY_PATH"
    if (type(values) == type([])):
      val = os.pathsep.join(values)
      print >> f, 'if (! $?%s) then' % (var_name,)
      print >> f, '  setenv %s ""' % (var_name,)
      print >> f, 'endif'
      print >> f, 'setenv %s "%s%s$%s"' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'setenv %s "%s"' % (var_name, values)
  f.close()

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
  print >> f, 'set PATHEXT=.PY;%PATHEXT%'
  f.close()

def emit_SConstruct(env):
  SConstruct_path = norm(join(env.LIBTBX_BUILD, "SConstruct"))
  f = open_info(SConstruct_path)
  print >> f, 'import os, os.path'
  print >> f, 'norm = os.path.normpath'
  print >> f, 'assert norm(os.getcwd()) == norm(os.environ["LIBTBX_BUILD"])'
  print >> f, 'Repository(r"%s")' % (env.LIBTBX_DIST_ROOT,)
  print >> f, 'try:'
  print >> f, '  CScanSetFlags('
  print >> f, '    python=0,'
  for package_name in env.package_list:
    flag = 1
    if (package_name == "boost"): flag = 0
    print >> f, '    %s=%d,' % (package_name, flag)
  print >> f, '  )'
  print >> f, 'except:'
  print >> f, '  pass'
  print >> f, '#SetContentSignatureType("timestamp")'
  print >> f, 'SConscript("libtbx/SConscript")'
  print >> f, '''\

def use_SConscript_if_present(package_name):
  dist = os.environ[package_name.upper() + "_DIST"]
  if (os.path.isfile(dist + "/SConscript")):
    SConscript(package_name + "/SConscript")
'''
  for package_name in env.package_list:
    print >> f, 'use_SConscript_if_present("%s")' % package_name
  f.close()

def run(libtbx_dist, args):
  env = libtbx_env(os.getcwd(), libtbx_dist)
  packages = registry()
  for arg in args:
    packages.merge(package(env.LIBTBX_DIST_ROOT, arg).dependency_registry)
  if (len(packages.list) == 0):
    print "Error: At least one package must be specified."
    return
  print "Top-down list of all packages involved:"
  for package_name in packages.list:
    print " ", package_name
    env.add_package(packages.dict[package_name])
  env.pickle_dict()
  if (hasattr(os, "symlink")):
    emit_env_run_sh(env)
    emit_setpaths_csh(env)
  else:
    emit_setpaths_bat(env)
  emit_SConstruct(env)
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
  warm_start(sys.argv)
