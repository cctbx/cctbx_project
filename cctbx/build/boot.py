import sys, os, shutil

def create_makefile(path_cctbx, configuration, subdir):
  print "Creating Makefile in " + subdir
  exe_g = {}
  exe_l = {}
  execfile(path_cctbx + "/" + subdir + "/Dependencies.py", exe_g, exe_l)
  try:
    os.makedirs(subdir)
  except OSError:
    pass
  h = exe_l['write_makefiles'](subdir, configuration)
  f = open(subdir + "/Makefile", "wb")
  h.write(f)
  f.close()
  shutil.copy(subdir + "/Makefile", subdir + "/Makefile.nodepend")

if (__name__ == "__main__"):
  path_cctbx = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
  sys.path.insert(0, path_cctbx + "/build")
  from make import read_configuration
  cf = read_configuration(path_root = os.path.dirname(path_cctbx))
  platform = cf[0]
  if (    platform in ("vc60", "mingw32", "win32_mwcc")
      and hasattr(os, "symlink")):
    print "Error: Must run under Windows!"
    sys.exit(1)
  if (len(sys.argv) == 1):
    for subdir in ("external/boost_python",
                   "eltbx",
                   "uctbx",
                   "sgtbx",
                   "adptbx",
                   "sftbx",
                   "fftbx",
                   "examples/cpp"):
      create_makefile(path_cctbx, cf, subdir)
  else:
    for subdir in sys.argv[1:]:
      create_makefile(path_cctbx, cf, subdir)
  if (hasattr(os, "symlink")):
    try: os.symlink(path_cctbx + "/examples/python", "examples/python")
    except: pass
  if (hasattr(os, "symlink")):
    for file in ("Makefile", "make.py", "test_imports.py"):
      print "Linking:", file
      try: os.symlink(path_cctbx + "/build/" + file, file)
      except: pass
    file = "setpythonpath.csh"
    set = "setenv PYTHONPATH '.:%s/lib_python'\n"
  else:
    for file in ("make.py", "test_imports.py"):
      print "Copying:", file
      shutil.copy(path_cctbx + "/build/" + file, file)
    file = "setpythonpath.bat"
    set = "set PYTHONPATH=.;%s\\lib_python\n"
  if (not os.path.exists(file)):
    print "Creating:", file
    f = open(file, "w")
    f.write(set % (os.getcwd(),))
    f.close()
