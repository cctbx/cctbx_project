import sys, os, shutil
norm = os.path.normpath

this_script = norm(os.path.abspath(sys.argv[0]))
build_dir = norm(os.path.dirname(this_script))
cctbx_dir = norm(os.path.dirname(build_dir))

if (len(sys.argv) != 2):
  print "usages:"
  print " ", os.path.basename(this_script), "<platform-identifier>"
  print " ", os.path.basename(this_script), "<cctbx-configuration-file>"
  sys.exit(1)

config_file = sys.argv[1]
if (os.path.isdir(config_file) or not os.path.exists(config_file)):
  config_file = norm(build_dir + "/configuration_" + sys.argv[1])
  if (not os.path.exists(config_file)):
    print "Files not found:"
    print " ", sys.argv[1]
    print " ", config_file
    sys.exit(1)
config_file = norm(os.path.abspath(config_file))

sys.path.insert(0, build_dir)
from make import read_configuration
cf = read_configuration(config_file)
platform = cf[0]

if (platform == "vc60"):
  make = "nmake"
else:
  make = "make"

if (hasattr(os, "symlink")):
  link_cmd = "softlinks"
else:
  link_cmd = "copy"

if (os.path.isdir("cctbx")):
  print "Warning: directory cctbx already exists"
else:
  os.mkdir("cctbx")
os.chdir("cctbx")
shutil.copy(config_file, "configuration")

def system_checked(cmd):
  if (os.system(cmd) != 0):
    sys.exit(1)

system_checked("python " + norm(build_dir + "/boot.py"))
system_checked("python make.py " + link_cmd)
system_checked("python make.py all")
system_checked("python test_imports.py")
system_checked(norm("examples/cpp/getting_started"))
os.environ["PYTHONPATH"] = norm(os.getcwd() + "/lib_python")
system_checked(
  "python " + norm(cctbx_dir + "/examples/python/getting_started.py"))
print "Done."
