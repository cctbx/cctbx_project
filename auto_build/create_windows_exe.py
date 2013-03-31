
# XXX like the Mac equivalent, this module is designed to be run independently
# of the rest of CCTBX if necessary, although it will use installed resources
# if found

from __future__ import division
try :
  import libtbx.load_env as libtbx_env
except ImportError, e :
  libtbx_env = None
import argparse
import shutil
import os
import sys

def run (args, out=sys.stdout) :
  if (sys.platform != "win32") :
    print >> out, "This application will only run on Windows systems."
    return 1
  parser = argparse.ArgumentParser(
    description="Utility for creating an iconified Windows launcher for the specified command, which must be present in %LIBTBX_BUILD%\\bin.")
  bin_path = icns_path = None
  if (libtbx_env is not None) :
    bin_path = os.path.join(abs(libtbx_env.build_path), "bin")
    ico_path = libtbx_env.find_in_repositories(
      relative_path="gui_resources/icons/custom/WinPhenix.ico",
      test=os.path.exists)
  else :
    bin_path = os.getcwd()
  parser.add_argument("--bin_dir", dest="bin_dir", action="store",
    help="Directory containing target executable or batch script.",
    default=bin_path)
  parser.add_argument("--exe_name", dest="exe_name", action="store",
    help="Name of iconified program", default=None)
  parser.add_argument("--icon", dest="icon", action="store",
    help="Path to .ico file", default=ico_path)
  parser.add_argument("--dest", dest="dest", action="store",
    help="Destination path", default=os.getcwd())
  options, args = parser.parse_known_args(args)
  if (len(args) == 0) :
    return parser.error("Executable name not specified.")
  if (options.bin_dir is None) :
    return parser.error("Executables directory not specified.")
  program_name = args[-1]
  bin_dir = options.bin_dir
  bin_files = os.listdir(bin_dir)
  program_cmd_file = None
  for file_name in bin_files :
    base, ext = os.path.splitext(file_name)
    if (base == program_name) :
      if (ext == ".bat") : # preferred
        program_cmd_file = file_name
        break
      elif (ext == ".exe") :
        program_cmd_file = file_name
  if (program_cmd_file is None) :
    print >> out, "No program named '%s' found in %s." % (program_name,
      bin_dir)
    return 1
  exe_name = program_name
  if (options.exe_name is not None) :
    exe_name = options.exe_name
  if (os.path.isdir("py2exe_tmp")) :
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
""" % os.path.join(bin_dir, program_cmd_name))
  f.close()
  bundle_files = 3
  if (options.bundle_all) :
    bundle_files = 1 # won't work on win64
  icon_rsrc = ""
  if (options.icon is not None) :
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
       'bundle_files': %d,
    }
  })
""" % (exe_name, icon_rsrc, bundle_files))
  f2.close()
  # XXX use sys.executable to avoid -Qnew behavior (crashes py2exe)
  import subprocess
  rc = subprocess.call([sys.executable, "setup.py", "py2exe"])
  if (rc != 0) :
    return rc
  exe_file = os.path.join("dist", "%s.exe" % exe_name)
  assert (os.path.isfile(exe_file))
  os.chdir(options.dest)
  if (os.path.exists("%s.exe" % exe_name)) :
    shutil.rmtree("%s.exe" % exe_name)
  shutil.move(exe_file, os.getcwd())
  print >> out, "Created %s" % os.path.join(os.getcwd(), "%s.exe" % exe_name)
  return 0

if (__name__ == "__main__") :
  sys.exit(run(sys.argv[1:]))
