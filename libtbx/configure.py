import sys, os, os.path, pprint
norm = os.path.normpath
join = os.path.join

class empty: pass

def open_info(path, mode="w", info="Creating:"):
  print info, path
  return open(path, mode)

def insert_normed_path(path_list, addl_path):
  addl_path = norm(addl_path)
  if (not addl_path in path_list):
    i = 0
    if (path_list[:1] == ["."]): i = 1
    path_list.insert(i, addl_path)

def update_libtbx_info(env, package_name, package_dist):
  libtbx_info_path = norm(join(env.libtbx_build, "libtbx_info"))
  package_dist_varname = package_name.upper() + "_DIST"
  libtbx_info = {
    "LD_LIBRARY_PATH": [],
    "PYTHONPATH": [norm(join(env.libtbx_build, "libtbx")), env.libtbx_dist],
    "PATH": [norm(join(env.libtbx_dist, "libtbx/command_line"))],
    "LIBTBX_PYTHON_BIN": env.python_bin,
    "LIBTBX_DIST": env.libtbx_dist,
    "LIBTBX_SCONS": env.libtbx_scons,
    "LIBTBX_BUILD": env.libtbx_build,
    package_dist_varname: package_dist,
    "package_names": [],
  }
  try:
    f = open(libtbx_info_path, "r")
  except:
    libtbx_info_write_action = "Creating:"
  else:
    libtbx_info_write_action = "Updating:"
    libtbx_info.update(eval(" ".join(f.readlines())))
    f.close()
  if (libtbx_info["LIBTBX_PYTHON_BIN"] != env.python_bin):
    print >> sys.stderr
    print >> sys.stderr, "#" * 78
    print >> sys.stderr, "FATAL ERROR: python executable mismatch:"
    print >> sys.stderr, "  previously:", libtbx_info["LIBTBX_PYTHON_BIN"]
    print >> sys.stderr, "         now:", sys.executable
    print >> sys.stderr, "#" * 78
    print >> sys.stderr
    sys.exit(1)
  assert libtbx_info["LIBTBX_DIST"] == env.libtbx_dist
  assert libtbx_info["LIBTBX_BUILD"] == env.libtbx_build
  assert libtbx_info[package_dist_varname] == package_dist
  insert_normed_path(libtbx_info["LD_LIBRARY_PATH"],
    join(env.libtbx_build, "libtbx"))
  insert_normed_path(libtbx_info["PYTHONPATH"],
    package_dist)
  insert_normed_path(libtbx_info["PATH"],
    join(package_dist, package_name, "command_line"))
  package_names = libtbx_info["package_names"]
  if (not package_name in package_names):
    package_names.append(package_name)
  f = open_info(libtbx_info_path, "w", libtbx_info_write_action)
  pprint.pprint(libtbx_info, f)
  f.close()
  return libtbx_info

def emit_env_run_sh(libtbx_build, libtbx_info):
  env_run_sh_path = norm(join(libtbx_build, "env_run.sh"))
  f = open_info(env_run_sh_path)
  for var_name, values in libtbx_info.items():
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
  print >> f, 'PATH="$LIBTBX_PYTHON_BIN%s$PATH"' % (os.pathsep,)
  print >> f, 'python "$LIBTBX_DIST/libtbx/command_line/env_run.py" $*'
  f.close()
  os.chmod(env_run_sh_path, 0755)

def emit_setpaths_csh(libtbx_build, libtbx_info):
  setpaths_csh_path = norm(join(libtbx_build, "setpaths.csh"))
  f = open_info(setpaths_csh_path)
  for var_name, values in libtbx_info.items():
    if (var_name.upper() != var_name): continue
    if (type(values) == type([])):
      val = os.pathsep.join(values)
      print >> f, 'if (! $?%s) then' % (var_name,)
      print >> f, '  setenv %s ""' % (var_name,)
      print >> f, 'endif'
      print >> f, 'setenv %s "%s%s$%s"' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'setenv %s "%s"' % (var_name, values)
  f.close()

def join_path_ld_library_path(libtbx_info):
  joined_path = libtbx_info["PATH"]
  for path in libtbx_info["LD_LIBRARY_PATH"]:
    if (not path in joined_path):
      joined_path.append(path)
  return joined_path

def emit_setpaths_bat(libtbx_build, libtbx_info):
  setpaths_bat_path = norm(join(libtbx_build, "setpaths.bat"))
  f = open_info(setpaths_bat_path)
  print >> f, '@ECHO off'
  for var_name, values in libtbx_info.items():
    if (var_name.upper() != var_name): continue
    if (type(values) == type([])):
      if (var_name == "LD_LIBRARY_PATH"): continue
      if (var_name == "PATH"):
        values = join_path_ld_library_path(libtbx_info)
      val = os.pathsep.join(values)
      print >> f, 'if not defined %s set %s=' % (var_name, var_name)
      print >> f, 'set %s=%s%s%%%s%%' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'set %s=%s' % (var_name, values)
  print >> f, 'if not defined PATHEXT set PATHEXT='
  print >> f, 'set PATHEXT=.PY;%PATHEXT%'
  f.close()

def find_scons(env):
  from_env = 0
  try: env.libtbx_scons = os.environ["LIBTBX_SCONS"]
  except: env.libtbx_scons = norm(join(env.dist_root, "scons/src/engine"))
  else: from_env = 1
  sys.path.insert(0, env.libtbx_scons)
  try: import SCons
  except: del sys.path[0]
  else: return
  env.libtbx_scons = None
  if (not from_env):
    try: import SCons
    except: pass
    else: env.libtbx_scons = "default"
  if (env.libtbx_scons == None):
    print >> sys.stderr, "Cannot find SCons (Software Construction Tool)"
    print >> sys.stderr, "Please refer to file: XXX"
    sys.exit(1)

def emit_SConstruct(env, libtbx_info):
  SConstruct_path = norm(join(env.libtbx_build, "SConstruct"))
  f = open_info(SConstruct_path)
  print >> f, 'Repository(r"%s")' % (env.dist_root,)
  print >> f, 'SConscript("libtbx/SConscript")'
  for package_name in libtbx_info["package_names"]:
    print >> f, 'SConscript("%s/SConscript")' % (package_name,)
  f.close()

def run():
  env = empty()
  env.python_bin = norm(os.path.split(sys.executable)[0])
  env.libtbx_build = norm(os.path.abspath(os.getcwd()))
  env.libtbx_dist = norm(os.path.dirname(norm(os.path.abspath(sys.argv[0]))))
  env.dist_root = norm(os.path.dirname(env.libtbx_dist))
  find_scons(env)
  libtbx_info = None
  for arg in sys.argv[1:]:
    package_dist = norm(join(env.dist_root, arg))
    if (not os.path.isdir(package_dist)):
      print "Error: Not a directory or platform identifier:", arg
      sys.exit(1)
    else:
      libtbx_info = update_libtbx_info(env, arg, package_dist)
      if (hasattr(os, "symlink")):
        emit_env_run_sh(env.libtbx_build, libtbx_info)
        emit_setpaths_csh(env.libtbx_build, libtbx_info)
      else:
        emit_setpaths_bat(env.libtbx_build, libtbx_info)
  if (libtbx_info == None):
    print "Specify at least one package name."
  else:
    emit_SConstruct(env, libtbx_info)

if (__name__ == "__main__"):
  run()
