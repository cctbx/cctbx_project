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

def read_configuration(config_file = "configuration", path_root = None):
  try:
    f = open(config_file, "r")
  except:
    print "Cannot open configuration file:", config_file
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
  python_executable = sys.executable
  if (cf[0] in ("vc60", "win32_mwcc")):
    python_include = sys.prefix + r"\include"
    python_lib = sys.prefix + r"\libs\python%s%s.lib" % (
      sys.version[0], sys.version[2])
  elif (cf[0] == "mingw32"):
    python_include = r"$(MINGW32_USR)\include\python" + sys.version[0:3]
    python_lib = r"$(MINGW32_USR)\lib\libpython%s%s.a" % (
      sys.version[0], sys.version[2])
  else:
    python_include = sys.prefix + "/include/python" + sys.version[0:3]
    python_lib = "%s/lib/python%s/config/libpython%s.a" % (
      sys.prefix, sys.version[0:3], sys.version[0:3])
  expand_cf(cf, "@(python_executable)", python_executable, normpath = 1)
  expand_cf(cf, "@(python_include)", python_include, normpath = 1)
  expand_cf(cf, "@(python_lib)", python_lib, normpath = 1)
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

if (__name__ == "__main__"):
  cf = read_configuration()
  platform = cf[0]
  if (platform in ("vc60", "win32_mwcc") and not hasattr(os, "symlink")):
    make = "nmake"
  else:
    make = "make"
  make_all = "all" in sys.argv

  externals = ("external/boost_python",)
  toolboxes = ("uctbx", "sgtbx", "arraytbx",
               "adptbx", "eltbx", "sftbx", "fftbx", "lbfgs")
  examples = ("examples/cpp",)
  all_targets = externals + toolboxes + examples

  if (hasattr(os, "symlink")):
    if ("softlinks" in sys.argv):
      run_in_subdirs(all_targets, make + " softlinks")
    if ("unlink" in sys.argv):
      run_in_subdirs(all_targets, make + " unlink")
    if ("rm" in sys.argv):
      run_in_subdirs(all_targets, make + " rm")
  else:
    if ("copy" in sys.argv):
      run_in_subdirs(all_targets, make + " copy")
    if ("del" in sys.argv):
      run_in_subdirs(all_targets, make + " del")

  if ("clean" in sys.argv):
    run_in_subdirs(all_targets, make + " clean")

  if (platform != "vc60"):
    if ("depend" in sys.argv):
      run_in_subdirs(all_targets,
        make + " -f Makefile.nodepend depend > Makefile", verbose = 1)

  if ("compile" in sys.argv or make_all):
    run_in_subdirs(all_targets, make + " compile")
