import sys, os, shutil

def create_makefiles(path_cctbx, configuration):
  for subdir in ("eltbx", "sgtbx", "uctbx", "examples/cpp"):
    print "Creating Makefile in " + subdir
    exe_g = {}
    exe_l = {}
    execfile(path_cctbx + "/" + subdir + "/Dependencies.py", exe_g, exe_l)
    try:
      os.makedirs(subdir)
    except OSError:
      pass
    h = exe_l['write_makefiles'](configuration)
    f = open(subdir + "/Makefile", "wb")
    h.write(f)
    f.close()
    shutil.copy(subdir + "/Makefile", subdir + "/Makefile.nodepend")

if (__name__ == "__main__"):
  path_cctbx = os.path.dirname(os.path.dirname(sys.argv[0]))
  sys.path.insert(0, path_cctbx + "/build")
  from make import read_configuration
  cf = read_configuration(os.path.dirname(path_cctbx))
  platform = cf[0]
  if (platform in ("vc60", "mingw32") and hasattr(os, "symlink")):
    print "Error: Must run under Windows!"
    sys.exit(1)
  create_makefiles(path_cctbx, cf)
  if (hasattr(os, "symlink")):
    for file in ("Makefile", "make.py", "test.py"):
      print "Linking:", file
      try: os.symlink(path_cctbx + "/build/" + file, file)
      except: pass
    file = "setpythonpath.csh"
    set = "setenv PYTHONPATH '%s/lib/python'\n"
  else:
    for file in ("make.py", "test.py"):
      print "Copying:", file
      shutil.copy(path_cctbx + "/build/" + file, file)
    file = "setpythonpath.bat"
    set = "set PYTHONPATH=%s\\lib\\python\n"
  if (not os.path.exists(file)):
    print "Creating:", file
    f = open(file, "w")
    f.write(set % (os.getcwd(),))
    f.close()
