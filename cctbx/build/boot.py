import sys, os, shutil, pprint
norm              = os.path.normpath
join              = os.path.join
pyexe             = sys.executable
this_script       = norm(os.path.abspath(sys.argv[0]))
build_system_dir  = norm(os.path.dirname(this_script))

def create_makefile(pkg_dir, configuration, subdir, package):
  print "Creating Makefile in " + subdir
  exe_g = {}
  exe_l = {}
  execfile(pkg_dir + "/" + subdir + "/Dependencies.py", exe_g, exe_l)
  try:
    os.makedirs(subdir)
  except OSError:
    pass
  h = exe_l['write_makefiles'](subdir, configuration, package)
  f = open(subdir + "/Makefile", "wb")
  h.write(f)
  f.close()
  shutil.copy(subdir + "/Makefile", subdir + "/Makefile.nodepend")

def create_lib_python_init_py(path_name):
  f = open(path_name + "/__init__.py", "w")
  print >> f, """import sys
if (sys.platform == "linux2"):
  if (hasattr(sys, "setdlopenflags")):
    sys.setdlopenflags(0x100|0x2)
"""
  f.close()

def create_lib_dir_and_lib_python_dir(pkg):
  try: os.makedirs("../lib")
  except OSError: pass
  lib_python_dir = "../lib_python/%s_boost"%(pkg.name)
  try: os.makedirs(lib_python_dir)
  except OSError: pass
  create_lib_python_init_py(lib_python_dir)
  for subdir in pkg.tbx_subpkg:
    try: os.makedirs(lib_python_dir + "/" + subdir)
    except OSError: pass
    create_lib_python_init_py(lib_python_dir + "/" + subdir)

def insert_normed_path(path_list, addl_path):
  addl_path = norm(addl_path)
  if (not addl_path in path_list):
    i = 0
    if (path_list[:1] == ["."]): i = 1
    path_list.insert(i, addl_path)

def update_paths_raw(build_root_dir, package_name, structured_package_dir):
  paths_raw_path = norm(join(build_root_dir, "paths_raw"))
  python_bin = norm(os.path.split(sys.executable)[0])
  package_dist = package_name.upper() + "_DIST"
  paths_raw_dict = {
    "LD_LIBRARY_PATH": [],
    "PYTHONPATH": [],
    "PATH": [],
    "PHENIX_PYTHON_BIN": python_bin,
    "PHENIX_BUILD": build_root_dir,
    package_dist: norm(structured_package_dir)
  }
  try:
    f = open(paths_raw_path, "r")
  except:
    pass
  else:
    paths_raw_dict.update(eval(" ".join(f.readlines())))
    f.close()
  if (paths_raw_dict["PHENIX_PYTHON_BIN"] != python_bin):
    print >> sys.stderr
    print >> sys.stderr, "#" * 78
    print >> sys.stderr, "FATAL ERROR: python executable mismatch:"
    print >> sys.stderr, "  previously:", paths_raw_dict["PHENIX_PYTHON_BIN"]
    print >> sys.stderr, "         now:", sys.executable
    print >> sys.stderr, "#" * 78
    print >> sys.stderr
    sys.exit(1)
  assert paths_raw_dict["PHENIX_BUILD"] == build_root_dir
  assert paths_raw_dict[package_dist] == norm(structured_package_dir)
  insert_normed_path(paths_raw_dict["LD_LIBRARY_PATH"],
    join(build_root_dir, "lib"))
  insert_normed_path(paths_raw_dict["PYTHONPATH"],
    join(build_root_dir, "lib_python"))
  insert_normed_path(paths_raw_dict["PYTHONPATH"],
    structured_package_dir)
  insert_normed_path(paths_raw_dict["PATH"],
    join(structured_package_dir, package_name, "command_line"))
  f = open(paths_raw_path, "w")
  pprint.pprint(paths_raw_dict, f)
  f.close()
  return paths_raw_dict

def emit_env_run_sh(build_root_dir, paths_raw_dict):
  env_run_sh_path = norm(join(build_root_dir, "env_run.sh"))
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
  print >> f, 'PATH="$PHENIX_PYTHON_BIN%s$PATH"' % (os.pathsep,)
  print >> f, 'python "$CCTBX_DIST/cctbx/command_line/env_run.py" $*'
  os.chmod(env_run_sh_path, 0755)
  f.close()

def emit_setpaths_csh(build_root_dir, paths_raw_dict):
  setpaths_csh_path = norm(join(build_root_dir, "setpaths.csh"))
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

def emit_setpaths_bat(build_root_dir, paths_raw_dict):
  setpaths_bat_path = norm(join(build_root_dir, "setpaths.bat"))
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

def run():
  sys.path.insert(0,os.getcwd())
  from PackageStructure import PackageStructure as package
  print "***Booting package: ", package.name.upper()
  pks = (package.name,) + package.supporting

  custom_subdir = None
  if len(sys.argv) == 2:
    structured_package_dir = sys.argv[1]
  elif os.path.isfile("PackageSource"):
    f = open("PackageSource","r")
    structured_package_dir = f.readline().strip()
    f.close()
    if (len(sys.argv) == 3):
      assert sys.argv[1] == "--custom"
      custom_subdir = sys.argv[2]
    else:
      assert len(sys.argv) == 1

  assert os.path.isdir(structured_package_dir)

  sys.path.insert(0, build_system_dir)
  from make import read_configuration

  cf = read_configuration(path_root = os.path.dirname(structured_package_dir),
                          supporting = pks)
  platform = cf[0]
  if (    platform in ("vc60", "mingw32", "win32_mwcc")
      and hasattr(os, "symlink")):
    print "Error: Must run under Windows!"
    sys.exit(1)

  if (custom_subdir):
    create_makefile(structured_package_dir, cf, custom_subdir, package)
    return

  for subdir in package.externals + package.toolboxes + package.examples:
    create_makefile(structured_package_dir, cf, subdir, package)

  if (hasattr(os, "symlink")):
    try: os.symlink(
      structured_package_dir + "/examples/python", "examples/python")
    except: pass
  create_lib_dir_and_lib_python_dir(package)

  build_root_dir = norm(join(os.getcwd(), ".."))
  paths_raw_dict = update_paths_raw(
    build_root_dir, package.name, structured_package_dir)
  if (hasattr(os, "symlink")):
    emit_env_run_sh(build_root_dir, paths_raw_dict)
    emit_setpaths_csh(build_root_dir, paths_raw_dict)
  else:
    emit_setpaths_bat(build_root_dir, paths_raw_dict)

  if (hasattr(os, "symlink")):
    for file in ("Makefile", "test_imports.py"):
      print "Linking:", file
      try: os.symlink(structured_package_dir + "/build/" + file, file)
      except: pass
  else:
    for file in ( "test_imports.py",):
      print "Copying:", file
      shutil.copy(structured_package_dir + "/build/" + file, file)

if (__name__ == "__main__"):
  run()
