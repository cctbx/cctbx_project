"""
This module is deprecated and no longer used as of Phenix Version 1.10.1-2155
where create_windows_installer.py supersedes it to invoke creating a NullSoft setup installer
"""


# XXX like the Mac equivalent, this module is designed to be run independently
# of the rest of CCTBX if necessary, although it will use installed resources
# if found

from __future__ import absolute_import, division, print_function
try :
  import libtbx.load_env
  libtbx_env = libtbx.env
except ImportError:
  libtbx_env = None
import optparse
import shutil
import stat
import os
import sys

def run(args, out=sys.stdout):
  if (sys.platform != "win32"):
    print("This application will only run on Windows systems.", file=out)
    return 1
  parser = optparse.OptionParser(
    description="Utility for creating an iconified Windows launcher for the specified command, which must be present in %LIBTBX_BUILD%\\bin.")
  bin_path = icns_path = None
  if (libtbx_env is not None):
    bin_path = os.path.join(abs(libtbx_env.build_path), "bin")
    ico_path = libtbx_env.find_in_repositories(
      relative_path="gui_resources/icons/custom/WinPhenix.ico",
      test=os.path.exists)
  else :
    bin_path = os.getcwd()
  parser.add_option("--bin_dir", dest="bin_dir", action="store",
    help="Directory containing target executable or batch script.",
    default=bin_path)
  parser.add_option("--exe_name", dest="exe_name", action="store",
    help="Name of iconified program", default=None)
  parser.add_option("--icon", dest="icon", action="store",
    help="Path to .ico file", default=ico_path)
  parser.add_option("--dest", dest="dest", action="store",
    help="Destination path", default=os.getcwd())
  parser.add_option("--bundle_all", dest="bundle_all", action="store_true",
    help="Bundle Python interpreter, etc. into .exe", default=False)
  options, args = parser.parse_args(args)
  if (len(args) == 0):
    return parser.error("Executable name not specified.")
  if (options.bin_dir is None):
    return parser.error("Executables directory not specified.")
  program_name = args[-1]
  bin_dir = options.bin_dir
  bin_files = os.listdir(bin_dir)
  program_cmd_file = None
  for file_name in bin_files :
    base, ext = os.path.splitext(file_name)
    if (base == program_name):
      if (ext == ".bat") : # preferred
        program_cmd_file = file_name
        break
      elif (ext == ".exe"):
        program_cmd_file = file_name
  if (program_cmd_file is None):
    print("No program named '%s' found in %s." % (program_name,
      bin_dir), file=out)
    return 1
  exe_name = program_name
  if (options.exe_name is not None):
    exe_name = options.exe_name
  if (os.path.isdir("py2exe_tmp")):
    shutil.rmtree("py2exe_tmp")
  os.mkdir("py2exe_tmp")
  os.chdir("py2exe_tmp")
  f = open("%s.py" % exe_name, "w")
  # XXX for reasons unknown to me, the method used on Mac (os.spawnv) will
  # not work for windows, but subprocess.call appears to do the job nicely,
  # with the minor complaint that it leaves the phenix.exe command running
  # (and visible in the taskbar) until the actual app closes.
  f.write("""
import subprocess
subprocess.call(r"%s")
""" % os.path.join(bin_dir, program_cmd_file))
  f.close()
  bundle_files = 3
  #zip_file = "'%s.zip'" %
  if (options.bundle_all):
    bundle_files = 1 # won't work on win64
    zip_file = None
  icon_rsrc = ""
  if (options.icon is not None):
    icon_rsrc = "'icon_resources':[(0,r'%s')]," % options.icon
  f2 = open("setup.py", "w")
  f2.write("""
from distutils.core import setup
import py2exe
setup(
  console=[
    {
      'script':'%s.py',
      %s
    },
  ],
  zipfile=None,
  options={
    'py2exe': {
       'includes': ['subprocess'],
       'dll_excludes' : ['w9xpopen.exe'],
       'bundle_files': %d,
    },
  })
""" % (exe_name, icon_rsrc, bundle_files))
  f2.close()
  # XXX use sys.executable to avoid -Qnew behavior (crashes py2exe)
  import subprocess
  rc = subprocess.call([sys.executable, "setup.py", "py2exe"])
  if (rc != 0):
    return rc
  dist_path = os.path.join(os.getcwd(), "dist")
  dist_files = os.listdir(dist_path)
  exe_file = os.path.join(dist_path, "%s.exe" % exe_name)
  assert (os.path.isfile(exe_file))
  os.chdir(options.dest)
  print("", file=out)
  for file_name in dist_files :
    if os.path.exists(file_name):
      print("WARNING: %s already exists" % file_name, file=out)
      continue
      # XXX Even for Windows, this is incredibly broken
      #full_path = os.path.join(options.dest, file_name)
      #print "removing %s" % file_name
      #subprocess.call("del /F /Q %s" % os.path.join(os.getcwd(), file_name))
      #os.chmod(full_path, stat.S_IWRITE)
      #os.unlink(full_path)
    print("moving %s..." % file_name, file=out)
    shutil.move(os.path.join(dist_path, file_name), os.getcwd())
    os.chmod(file_name, stat.S_IWRITE)
  return 0

if (__name__ == "__main__"):
  sys.exit(run(sys.argv[1:]))
