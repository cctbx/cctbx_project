import sys, os, shutil, traceback, string
norm              = os.path.normpath
join              = os.path.join
pyexe             = sys.executable
default_package   = "../../cctbx"

this_script       = norm(os.path.abspath(sys.argv[0]))
build_system_dir  = norm(os.path.dirname(this_script))
this_target       = os.getcwd()

all_build_files   = os.listdir(build_system_dir)
allowed_platforms = []
for filename in all_build_files:
  if filename.find("configuration_") == 0:
    allowed_platforms.append(filename.split("configuration_")[1])

__doc__="""Generic build system for C++ & Boost.Python modules
Usage:

argv[0] name of this script ("%s") located in %s
argv[1] platform for this build,
        <%s>,
        or path name of user-defined configuration file "configuration_xxx"
argv[2:] list of directories containing structured packages, default "%s"

"""%(os.path.basename(this_script),
     build_system_dir,
     string.join(allowed_platforms,"|"),
     os.path.basename(default_package))

def assert_dir(x): assert os.path.isdir(x); return 1

try:
  if sys.argv[1] in allowed_platforms:
    platform_identifier = sys.argv[1]
    configuration_file  = join(build_system_dir,
                              "configuration_%s"%(platform_identifier))
    assert (os.path.isfile(configuration_file))
  else:
    configuration_file  = norm(os.path.abspath(sys.argv[1]))
    assert (os.path.isfile(configuration_file))
    bases = os.path.basename(configuration_file).split("configuration_")
    assert bases[0]==""
    assert bases[1]!=None
    platform_identifier = bases[1]

  if len(sys.argv) == 2:
    structured_package_dirs = [join(build_system_dir,default_package)]
  else:
    structured_package_dirs = sys.argv[2:]

  structured_package_dirs = map(os.path.abspath,structured_package_dirs)
  [assert_dir(x) for x in structured_package_dirs]

except:
  print "Incorrect arguments for fast track build"
  traceback.print_exc()
  print __doc__
  sys.exit(1)

print "valid build call for",structured_package_dirs
sys.stdout.flush()

sys.path.insert(0, build_system_dir)
from make import read_configuration
cf = read_configuration(configuration_file)
platform = cf[0]
assert platform == platform_identifier

if (platform == "vc60"):
  make = "nmake"
else:
  make = "make"

if (hasattr(os, "symlink")):
  os.environ["LD_LIBRARY_PATH"] = os.path.abspath("./lib")
  link_cmd = "softlinks"
else:
  link_cmd = "copy"

def system_checked(cmd):
  if (os.system(cmd) != 0):
    sys.exit(1)
  sys.stdout.flush()
  sys.stderr.flush()

for sp in structured_package_dirs:
  os.chdir(this_target)
  try:
    structured_package = os.path.basename(sp)

    if (not os.path.isdir(structured_package)):
      os.mkdir(structured_package)
    os.chdir(structured_package)
    shutil.copy(configuration_file, "configuration")

    structure = join(sp,"PackageStructure.py")
    assert os.path.isfile(structure)
    shutil.copy(structure, "PackageStructure.py")

    f = open("PackageSource","w")
    f.write(sp)
    f.close()

    system_checked(pyexe + " " + norm(build_system_dir + "/boot.py") + " " + sp)
    system_checked(pyexe + " " + build_system_dir + "/make.py " + link_cmd)
    system_checked(pyexe + " " + build_system_dir + "/make.py all")
    system_checked(pyexe + " test_imports.py")
    if os.path.isfile(norm("examples/cpp/getting_started")):
      system_checked(norm("examples/cpp/getting_started"))
      os.environ["PYTHONPATH"] = norm(join(os.getcwd(),"../lib_python"))
      system_checked(
        pyexe + " " + norm(sp + "/examples/python/getting_started.py"))
    print "Done with %s." % (structured_package,)
  except:
    # for now, failure to compile sp doesn't kill entire list
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]
