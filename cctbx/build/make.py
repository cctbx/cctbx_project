import sys, os, shutil, string
join = os.path.join
norm = os.path.normpath

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

def read_configuration(config_file = "configuration", path_root = None, supporting=None):
  # Get the configuration template
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

  simple_replacements = {}
  if (path_root):
    path_root = os.path.normpath(path_root)
    while (path_root[-1] == os.sep): path_root = path_root[:-1]
    simple_replacements["@(ROOT)"] = path_root
  simple_replacements["@(CWD)"] = os.getcwd()
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
  simple_replacements["@(python_executable)"] = python_executable
  simple_replacements["@(python_include)"] = python_include
  simple_replacements["@(python_lib)"] = python_lib
  libdir = os.path.abspath(norm(join(os.getcwd(),"../lib")))
  libpyth= os.path.abspath(norm(join(os.getcwd(),"../lib_python")))
  trans_to_win32= string.maketrans("/", "\\")
  trans_to_unix = string.maketrans("\\", "/")
  simple_replacements["@(SP_LIBDIR_UNIX)"] = libdir.translate(trans_to_unix)
  simple_replacements["@(SP_LIBDIR_WIN)"] = libdir.translate(trans_to_win32)
  simple_replacements["@(SP_LIBPYTHON_ROOT_UNIX)"] =libpyth.translate(trans_to_unix)
  simple_replacements["@(SP_LIBPYTHON_ROOT_WIN)"] =libpyth.translate(trans_to_win32)

  if supporting: # List of supporting structured packages which must appear in makefile
                 # These string expansions must be done before the simple_replacements
    newcf = []
    for line in cf:
      if line.find("@(STRUCTPACK)") < 0:
        newcf.append(line)
      else:
        for package in supporting:
          split_line = line.split("@(STRUCTPACK)")
          assert len(split_line) >= 2 # depend on this structure; otherwise re-implement
          templine = split_line[0] + package.upper() + split_line[1]
          if len(split_line) > 2:
            for elem in split_line[2:]:
              templine += package + elem
          newcf.append(templine)
    cf = newcf

  for key in simple_replacements.keys():
    expand_cf(cf, key, simple_replacements[key], normpath = 1 )
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

def get_package_structure():
  sys.path.insert(0,os.getcwd())
  from PackageStructure import PackageStructure as package
  return package

if (__name__ == "__main__"):
  package = get_package_structure()
  cf = read_configuration()
  platform = cf[0]
  if (platform in ("vc60", "win32_mwcc") and not hasattr(os, "symlink")):
    make = "nmake"
  else:
    make = "make"
  make_all = "all" in sys.argv

  externals = package.externals
  toolboxes = package.toolboxes
  examples  = package.examples
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
