import sys, os

def create_autorun(bundle_prefix, single_dir=False):
  if (single_dir) :
    return """\
$AUTORUN$>.\\%(bundle_prefix)s_install_script.bat
.
""" % vars()
  else :
    return """\
$AUTORUN$>.\\%(bundle_prefix)s\\%(bundle_prefix)s_install_script.bat
.
""" % vars()

def run(args):
  no_unzipsfx = (len(args) > 0 and args[0] == "--no-unzipsfx")
  if (no_unzipsfx):
    args = args[1:]
  single_dir = (len(args) > 0 and args[0] == "--single_directory")
  if (single_dir) :
    args = args[1:]
  if (len(args) < 2):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage(
      "%s [--no-unzipsfx] [--single_directory] bundle_prefix platform_string [addl_files...]"
        % libtbx.env.dispatcher_name)
  if (os.name == "nt"):
    exe_suffix = ".exe"
  else:
    exe_suffix = ""
  import libtbx.path
  path_zip = libtbx.path.full_command_path(
    command="zip"+exe_suffix, search_first=["."])
  if (path_zip is None):
    raise RuntimeError("Fatal: zip executable not found.")
  bundle_prefix = args[0]
  if (single_dir) and (not os.path.isdir(bundle_prefix)) :
    from libtbx.utils import Sorry
    raise Sorry("%s does not exist or is not a directory." % bundle_prefix)
  platform_string = args[1]
  addl_files = args[2:]
  zip_file_name = "%(bundle_prefix)s_%(platform_string)s.zip" % vars()
  open("autorun", "w").write(create_autorun(bundle_prefix, single_dir))
  if (single_directory) :
    cmd = ("\"%(path_zip)s\" -q -r -z %(zip_file_name)s"
        + " %(bundle_prefix)s") % vars()
  else :
    cmd = ("\"%(path_zip)s\" -q -r -z %(zip_file_name)s"
        + " %(bundle_prefix)s_sources"
        + " %(bundle_prefix)s_build"
        + " %(bundle_prefix)s_install_script.bat") % vars()
  for addl in addl_files:
    cmd += " " + addl
  cmd += " < autorun"
  print cmd
  from libtbx import easy_run
  easy_run.fully_buffered(command=cmd).raise_if_errors().show_stdout()
  if (not no_unzipsfx):
    from libtbx.command_line import create_unzipsfx
    create_unzipsfx.create(zip_file_name=zip_file_name)

if (__name__ == "__main__"):
  run(sys.argv[1:])
