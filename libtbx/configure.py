import sys, os, os.path, pprint
norm = os.path.normpath
join = os.path.join

class empty: pass

def insert_normed_path(path_list, addl_path):
  addl_path = norm(addl_path)
  if (not addl_path in path_list):
    i = 0
    if (path_list[:1] == ["."]): i = 1
    path_list.insert(i, addl_path)

def update_paths_raw(env, package_name, package_dist):
  paths_raw_path = norm(join(env.libtbx_build, "paths_raw"))
  package_dist_varname = package_name.upper() + "_DIST"
  paths_raw_dict = {
    "LD_LIBRARY_PATH": [],
    "PYTHONPATH": [],
    "PATH": [],
    "LIBTBX_PYTHON_BIN": env.python_bin,
    "LIBTBX_DIST": env.libtbx_dist,
    "LIBTBX_BUILD": env.libtbx_build,
    "LIBTBX_PACKAGES": [],
    package_dist_varname: package_dist,
  }
  try:
    f = open(paths_raw_path, "r")
  except:
    pass
  else:
    paths_raw_dict.update(eval(" ".join(f.readlines())))
    f.close()
  if (paths_raw_dict["LIBTBX_PYTHON_BIN"] != env.python_bin):
    print >> sys.stderr
    print >> sys.stderr, "#" * 78
    print >> sys.stderr, "FATAL ERROR: python executable mismatch:"
    print >> sys.stderr, "  previously:", paths_raw_dict["LIBTBX_PYTHON_BIN"]
    print >> sys.stderr, "         now:", sys.executable
    print >> sys.stderr, "#" * 78
    print >> sys.stderr
    sys.exit(1)
  assert paths_raw_dict["LIBTBX_DIST"] == env.libtbx_dist
  assert paths_raw_dict["LIBTBX_BUILD"] == env.libtbx_build
  assert paths_raw_dict[package_dist_varname] == package_dist
  insert_normed_path(paths_raw_dict["LD_LIBRARY_PATH"],
    join(env.libtbx_build, "libtbx"))
  insert_normed_path(paths_raw_dict["PYTHONPATH"],
    join(env.libtbx_build, "libtbx"))
  insert_normed_path(paths_raw_dict["PYTHONPATH"],
    package_dist)
  insert_normed_path(paths_raw_dict["PATH"],
    join(package_dist, package_name, "command_line"))
  package_names = paths_raw_dict["LIBTBX_PACKAGES"]
  if (not package_name in package_names):
    package_names.append(package_name)
  f = open(paths_raw_path, "w")
  pprint.pprint(paths_raw_dict, f)
  f.close()
  return paths_raw_dict

def emit_env_run_sh(libtbx_build, paths_raw_dict):
  env_run_sh_path = norm(join(libtbx_build, "env_run.sh"))
  print "Updating:", env_run_sh_path
  f = open(env_run_sh_path, "w")
  for var_name, values in paths_raw_dict.items():
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
  print >> f, 'python "$LIBTBX_DIST/env_run.py" $*'
  os.chmod(env_run_sh_path, 0755)
  f.close()

def emit_setpaths_csh(libtbx_build, paths_raw_dict):
  setpaths_csh_path = norm(join(libtbx_build, "setpaths.csh"))
  print "Updating:", setpaths_csh_path
  f = open(setpaths_csh_path, "w")
  for var_name, values in paths_raw_dict.items():
    if (type(values) == type([])):
      val = os.pathsep.join(values)
      print >> f, 'if (! $?%s) then' % (var_name,)
      print >> f, '  setenv %s ""' % (var_name,)
      print >> f, 'endif'
      print >> f, 'setenv %s "%s%s$%s"' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'setenv %s "%s"' % (var_name, values)
  f.close()

def emit_setpaths_bat(libtbx_build, paths_raw_dict):
  setpaths_bat_path = norm(join(libtbx_build, "setpaths.bat"))
  print "Updating:", setpaths_bat_path
  f = open(setpaths_bat_path, "w")
  for var_name, values in paths_raw_dict.items():
    if (type(values) == type([])):
      if (var_name == "LD_LIBRARY_PATH"): continue
      val = os.pathsep.join(values)
      print >> f, 'if not defined %s set %s=' % (var_name, var_name)
      print >> f, 'set %s=%s%s%%%s%%' % (var_name, val, os.pathsep, var_name)
    else:
      print >> f, 'set %s=%s' % (var_name, values)
  print >> f, 'if not defined PATHEXT set PATHEXT='
  print >> f, 'set PATHEXT=.PY;%PATHEXT%'
  f.close()

def emit_SConstruct(env, paths_raw_dict):
  SConstruct_path = norm(join(env.libtbx_build, "SConstruct"))
  print "Updating:", SConstruct_path
  f = open(SConstruct_path, "w")
  print >> f, 'Repository(r"%s")' % (env.dist_root,)
  print >> f, 'SConscript("libtbx/build/SConscript")'
  print >> f, 'SConscript("libtbx/SConscript")'
  for package_name in paths_raw_dict["LIBTBX_PACKAGES"]:
    print >> f, 'SConscript("libtbx/%s_boost/SConscript")' % (package_name,)
    print >> f, 'SConscript("%s/SConscript")' % (package_name,)
  f.close()

def run():
  env = empty()
  env.python_bin = norm(os.path.split(sys.executable)[0])
  env.libtbx_build = norm(os.path.abspath(os.getcwd()))
  env.libtbx_dist = norm(os.path.dirname(norm(os.path.abspath(sys.argv[0]))))
  env.dist_root = norm(os.path.dirname(env.libtbx_dist))
  for arg in sys.argv[1:]:
    package_dist = norm(join(env.dist_root, arg))
    if (not os.path.isdir(package_dist)):
      print "Error: Not a directory or platform identifier:", arg
      sys.exit(1)
    else:
      paths_raw_dict = update_paths_raw(env, arg, package_dist)
      if (hasattr(os, "symlink")):
        emit_env_run_sh(env.libtbx_build, paths_raw_dict)
        emit_setpaths_csh(env.libtbx_build, paths_raw_dict)
      else:
        emit_setpaths_bat(env.libtbx_build, paths_raw_dict)
  emit_SConstruct(env, paths_raw_dict)

if (__name__ == "__main__"):
  run()
