import sys, os, shutil, string

def expand_cf(cf, pattern, substitute, normpath = 0):
  for i in xrange(len(cf)):
    c = string.find(cf[i], pattern)
    if (c >= 0):
      cf[i] = cf[i][:c] + substitute + cf[i][c + len(pattern):]
      if (normpath):
        e = string.find(cf[i], "=")
        if (e >= 0):
          e = e + 1
          cf[i] = cf[i][:e] + os.path.normpath(cf[i][e:])

def read_configuration(path_root = None):
  try:
    f = open("configuration", "r")
  except:
    print "Cannot open configuration file."
    sys.exit(1)
  cf = f.readlines()
  f.close()
  for i in xrange(len(cf)): cf[i] = cf[i][:-1] # remove new-line
  while (len(cf)):
    cf[0] = string.strip(cf[0])
    if (len(cf[0]) != 0 and cf[0][0] != "#"): break
    del cf[0] # remove leading comment lines
  if (path_root):
    path_root = os.path.normpath(path_root)
    while (path_root[-1] == os.sep): path_root = path_root[:-1]
    expand_cf(cf, "@(ROOT)", path_root, normpath = 1)
  expand_cf(cf, "@(CWD)", os.getcwd(), normpath = 1)
  return cf

def system_verbose(command):
  print command
  return os.system(command)

def copy_verbose(src, dst):
  print "copy %s %s" % (src, dst)
  return shutil.copy(src, dst)

def run_in_subdirs(subdirs, command_line, verbose = 0):
  for s in subdirs:
    if (verbose): print s + ": " + command_line
    cwd = os.getcwd()
    os.chdir(s)
    os.system(command_line)
    os.chdir(cwd)

def make_libdir(platform):
  try: os.makedirs("lib")
  except OSError: pass
  if (platform == "vc60"):
    libext = ".lib"
  else:
    libext = ".a"
  copy_verbose("eltbx/libeltbx" + libext, "lib")
  copy_verbose("sgtbx/libsgtbx" + libext, "lib")
  copy_verbose("uctbx/libuctbx" + libext, "lib")

def make_libpythondir(platform):
  try: os.makedirs("lib/python/eltbx")
  except OSError: pass
  open("lib/python/eltbx/__init__.py", "a+").close()
  if (platform in ("vc60", "mingw32")):
    libpyd = ".pyd"
  else:
    libpyd = ".so"
  if (hasattr(os, "symlink")):
    system_verbose("cp eltbx/*%s lib/python/eltbx" % (libpyd,))
    system_verbose("cp sgtbx/*%s lib/python" % (libpyd,))
    system_verbose("cp uctbx/*%s lib/python" % (libpyd,))
    system_verbose("cp adptbx/*%s lib/python" % (libpyd,))
  else:
    system_verbose(r"copy eltbx\*.pyd lib\python\eltbx")
    system_verbose(r"copy sgtbx\*.pyd lib\python")
    system_verbose(r"copy uctbx\*.pyd lib\python")
    system_verbose(r"copy adptbx\*.pyd lib\python")

if (__name__ == "__main__"):
  cf = read_configuration()
  platform = cf[0]
  if (platform == "vc60" and not hasattr(os, "symlink")):
    make = "nmake"
  else:
    make = "make"
  compile_dev = "compile_dev" in sys.argv
  compile_all = "compile_all" in sys.argv or compile_dev

  toolboxes = ("eltbx", "sgtbx", "uctbx", "adptbx")
  examples = ("examples/cpp",)

  if (hasattr(os, "symlink")):
    if ("softlinks" in sys.argv):
      run_in_subdirs(toolboxes + examples, make + " softlinks")
    if ("unlink" in sys.argv):
      run_in_subdirs(toolboxes + examples, make + " unlink")
    if ("rm" in sys.argv):
      run_in_subdirs(toolboxes + examples, make + " rm")
  else:
    if ("copy" in sys.argv):
      run_in_subdirs(toolboxes + examples, make + " copy")
    if ("del" in sys.argv):
      run_in_subdirs(toolboxes + examples, make + " del")

  if ("clean" in sys.argv):
    run_in_subdirs(toolboxes, make + " clean")

  if (platform != "vc60"):
    if ("depend" in sys.argv or compile_dev):
      run_in_subdirs(toolboxes + examples,
        make + " -f Makefile.nodepend depend > Makefile", verbose = 1)

  if ("compile" in sys.argv or compile_all):
    run_in_subdirs(toolboxes, make + " compile")

  if ("libdir" in sys.argv or compile_all):
    make_libdir(platform)

  if ("libpythondir" in sys.argv or compile_all):
    make_libpythondir(platform)

  if ("compile_examples" in sys.argv or compile_all):
    run_in_subdirs(examples, make + " compile")
