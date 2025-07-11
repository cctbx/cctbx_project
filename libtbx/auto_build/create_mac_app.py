
# XXX this module is designed to be run independently of the rest of CCTBX if
# necessary, although it will use installed resources if found

from __future__ import absolute_import, division, print_function

import optparse
import os
import re
import shutil
import sys

try :
  import libtbx.load_env
# can be either ImportError or SyntaxError, depending on whether we're using
# the bundled Python version or something more ancient
except Exception:
  libtbx_env = None
else :
  libtbx_env = libtbx.env

def run(args, out=sys.stdout):
  if (sys.platform != "darwin"):
    print("This application will only run on Mac systems.", file=out)
    return 1
  parser = optparse.OptionParser(
    description="Utility for creating an iconified Mac launcher for the specified command, which must be present in $LIBTBX_BUILD/bin.")
  bin_path = icns_path = None
  if (libtbx_env is not None):
    bin_path = os.path.join(abs(libtbx_env.build_path), "bin")
    icns_path = libtbx_env.find_in_repositories(
      relative_path="gui_resources/icons/custom/phenix.icns",
      test=os.path.exists)
  parser.add_option("--bin_dir", dest="bin_dir", action="store",
    help="Directory containing target executable.", default=bin_path)
  parser.add_option("--app_name", dest="app_name", action="store",
    help="Name of iconified program", default=None)
  parser.add_option("--icon", dest="icon", action="store",
    help="Path to .icns file", default=icns_path)
  parser.add_option("--dest", dest="dest", action="store",
    help="Destination path", default=os.getcwd())
  parser.add_option("--alias_build", dest="alias_build", action="store_true",
    help="Generate alias build without Python interpreter", default=False)
  parser.add_option("--python_interpreter", dest="python_interpreter",
    action="store", help="Python interpreter to use for final app",
    default=None)
  parser.add_option("--extra_lines", dest="extra_lines", action="store",
    help="Extra line(s) to be added to the Python script that is run",
    default="")
  options, args = parser.parse_args(args)
  if (len(args) == 0):
    return parser.error("Executable name not specified.")
  if (options.bin_dir is None):
    return parser.error("Executables directory not specified.")
  program_name = args[-1]
  build_dir = abs(libtbx.env.build_path)
  bin_dir = os.path.join(build_dir, "bin")
  if libtbx.env.installed:
    bin_dir = abs(libtbx.env.bin_path)
  if (not program_name in os.listdir(bin_dir)):
    print("No program named '%s' found in %s." % (program_name,
      bin_dir), file=out)
    return 1
  try :
    import py2app.script_py2applet
  except ImportError:
    print("py2app not installed.", file=out)
    return 1
  app_name = program_name
  if (options.app_name is not None):
    app_name = options.app_name
  py2app_bin_dir = "py2app_tmp"
  if (options.alias_build):
    py2app_bin_dir = os.path.join("py2app_bin", app_name)
  if (os.path.isdir(py2app_bin_dir)):
    shutil.rmtree(py2app_bin_dir)
  os.makedirs(py2app_bin_dir)
  os.chdir(py2app_bin_dir)
  f = open("%s.py" % app_name, "w")
  f.write("""
import os
import sys
os.environ["PYTHONPATH"] = ""
%s
os.spawnv(os.P_NOWAIT, "%s", ["%s"])
""" % (options.extra_lines, os.path.join(bin_dir, program_name), app_name))
  f.close()
  script_name = re.sub(".pyc$", ".py", py2app.script_py2applet.__file__)
  import subprocess
  executable = sys.executable
  if (options.python_interpreter is not None):
    executable = options.python_interpreter
  elif (libtbx_env is not None):
    executable = abs(libtbx.env.python_exe)
  try :
    args = [executable, "-c", "'import py2app'"]
    rc = subprocess.call(args)
    if (rc != 0) : raise RuntimeError("oops!")
  except RuntimeError :
    print("py2app not available, aborting .app creation.", file=out)
    return 1
  args = [executable, script_name, "--make-setup", "%s.py" % app_name]
  if (options.icon is not None):
    args.append(options.icon)
  rc = subprocess.call(args)
  if (rc != 0):
    return rc
  args = [executable, "setup.py", "py2app"]
  if (options.alias_build):
    args.append("-A")
  rc = subprocess.call(args)
  if (rc != 0):
    return rc
  app_path = os.path.abspath(os.path.join("dist", "%s.app" % app_name))
  assert os.path.isdir(app_path), app_path
  os.chdir(options.dest)
  if (os.path.exists("%s.app" % app_name)):
    shutil.rmtree("%s.app" % app_name)
  shutil.move(app_path, os.getcwd())
  print("Created %s" % os.path.join(os.getcwd(), "%s.app" % app_name), file=out)
  return 0

if (__name__ == "__main__"):
  sys.exit(run(sys.argv[1:]))
