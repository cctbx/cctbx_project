import sys, os, shutil
norm = os.path.normpath
PACKAGES = os.path.abspath(sys.argv[0])
if (len(sys.argv) != 4):
  print "usage: " + os.path.basename(PACKAGES) \
        + " <platform> <boost-makefile> <cctbx-configuration-file>"
  sys.exit(1)
for i in xrange(2):
  if (not os.path.exists(sys.argv[i + 2])):
    print "File not found: " + sys.argv[i + 2]
    sys.exit(1)
PACKAGES = norm(PACKAGES)
for i in xrange(3): PACKAGES = os.path.dirname(PACKAGES)
while (PACKAGES[-1] == os.sep): PACKAGES = PACKAGES[:-1]
platform = sys.argv[1]
boost_mak = os.path.abspath(sys.argv[2])
cctbx_cf = os.path.abspath(sys.argv[3])
if (hasattr(os, "symlink")):
  link_cmd = "softlinks"
else:
  link_cmd = "copy"
if (platform == "vc60"):
  make = "nmake"
else:
  make = "make"
if (os.path.isdir(platform)):
  print "directory " + platform + " already exists"
else:
  os.mkdir(platform)
os.chdir(platform)
# Run the cctbx boot command before installing Boost.Python.
# The boost script will issue a warning if the platform identifier
# is unsupported.
if (os.path.isdir("cctbx")):
  print "directory cctbx already exists"
else:
  os.mkdir("cctbx")
os.chdir("cctbx")
shutil.copy(cctbx_cf, "configuration")
if (os.system("python " + norm(PACKAGES + "/cctbx/build/boot.py")) != 0):
  sys.exit(1)
os.chdir("..")
# Now install Boost.Python.
if (os.path.isdir("boost")):
  print "directory boost already exists"
else:
  os.mkdir("boost")
os.chdir("boost")
shutil.copy(boost_mak, "Makefile")
os.system(make + ' ROOT="%s" %s' % (PACKAGES, link_cmd))
os.system(make + ' ROOT="%s"' % (PACKAGES,))
os.system(make + " test")
os.chdir("..")
# Continue the installation of the cctbx.
os.chdir("cctbx")
os.system("python make.py " + link_cmd)
os.system("python make.py compile_all")
os.system("python test.py")
os.system(norm("examples/cpp/getting_started"))
os.environ["PYTHONPATH"] = norm(os.getcwd() + "/lib/python")
os.system("python "
          + norm(PACKAGES + "/cctbx/examples/python/getting_started.py"))
print "Done."
