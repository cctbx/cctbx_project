import sys, os, shutil
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

  # emit the pythonpath command file
  if (hasattr(os, "symlink")):
    for file in ("Makefile", "test_imports.py"):
      print "Linking:", file
      try: os.symlink(structured_package_dir + "/build/" + file, file)
      except: pass
    file = norm("../setpythonpath.csh")
    proc = """
if ( ! $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH ""
endif
setenv LD_LIBRARY_PATH "%s:${LD_LIBRARY_PATH}"
if ( ! $?PYTHONPATH ) then
  setenv PYTHONPATH ""
endif
setenv PYTHONPATH ".:%s:%s:${PYTHONPATH}"
set path = (%s/command_line $path)
    """ % ( norm(join(os.getcwd(),"../lib")),
            structured_package_dir,
            norm(join(os.getcwd(),"../lib_python")),
            structured_package_dir)

  else:
    for file in ( "test_imports.py",):
      print "Copying:", file
      shutil.copy(structured_package_dir + "/build/" + file, file)
    file = norm("..\setpythonpath.bat")
    proc = r"""
if not defined PYTHONPATH set PYTHONPATH=
if not defined PATHEXT set PATHEXT=
set PYTHONPATH=.;%s;%s;%%PYTHONPATH%%
set PATH=%s\command_line;%%PATH%%
set PATHEXT=.PY;%%PATHEXT%%
    """ % ( structured_package_dir,
            norm(join(os.getcwd(),"../lib_python")),
            structured_package_dir)
  print "Updating:", file
  f = open(file, "a")
  f.write(proc)
  f.close()

if (__name__ == "__main__"):
  run()
