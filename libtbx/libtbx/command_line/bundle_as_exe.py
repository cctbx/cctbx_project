from libtbx.command_line import create_unzipsfx
import libtbx.config
import sys, os

def create_autorun(bundle_prefix):
  return """\
$AUTORUN$>.\\%(bundle_prefix)s_install_script.bat
.
""" % vars()

def run(args):
  "usage: libtbx.bundle_up bundle_prefix platform_string [addl_files...]"
  if (len(args) < 2):
    print run.__doc__
    return
  path_zip = libtbx.config.full_path(command="zip.exe", search_first=["."])
  if (path_zip is None):
    raise RuntimeError("Fatal: zip executable not found.")
  bundle_prefix = args[0]
  platform_string = args[1]
  addl_files = args[2:]
  zip_file_name = "%(bundle_prefix)s_%(platform_string)s.zip" % vars()
  open("autorun", "w").write(create_autorun(bundle_prefix))
  cmd = "%(path_zip)s -q -r -z %(zip_file_name)s" \
      + " %(bundle_prefix)s_sources" \
      + " %(bundle_prefix)s_build" \
      + " %(bundle_prefix)s_install_script.bat" % vars()
  for addl in addl_files:
    cmd += " " + addl
  cmd += " < autorun"
  print cmd
  os.system(cmd)
  create_unzipsfx.create(zip_file_name=zip_file_name)

if (__name__ == "__main__"):
  run(sys.argv[1:])
