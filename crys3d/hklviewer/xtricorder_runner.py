from __future__ import absolute_import, division, print_function
from libtbx import group_args
from libtbx.utils import Sorry
import os
import os.path
import sys

def redirect_stdout_to_file(file_path):
    """
    Redirects stdout to the given file using os.dup2.
    Works for Python 3.8+ and affects all C-level stdout writes.
    """
    try:
        # Open the target file for writing (create if not exists, truncate if exists)
        fd_target = os.open(file_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)

        # Duplicate the file descriptor for stdout (fd=1)
        os.dup2(fd_target, sys.stdout.fileno())

        # Close the original file descriptor (no longer needed)
        os.close(fd_target)

    except OSError as e:
        sys.stderr.write(f"Error redirecting stdout: {e}\n")
        sys.exit(1)



def external_cmd(parent, master_phil, firstpart):
  from phasertng.scripts import xtricorder
  # Provide a temp directory for xtricorder in current working directory and replace any
  # backslashes on Windows with forwardslashes for the sake of phasertng. Append a random
  # number to tempdir to avoid race conditions if another instance of xtricorder is running
  import random
  from pathlib import PurePath

  logfname = firstpart + "_xtricorder.log"
  original_stdout_fd = os.dup(sys.stdout.fileno())
  redirect_stdout_to_file( logfname )

  tempdir = PurePath(os.path.join( os.getcwd(), "HKLviewerXtricorder")).as_posix() + str(random.randrange(100000))
  os.mkdir(tempdir)
  parent.hklin =  firstpart + "_xtricorder.mtz"
  tabname = "Xtricorder"
  (retobj) = xtricorder.xtricorder(
  r'''phasertng {
              hklin.filename = "%s"
              reflections.wavelength = 1.0
              suite.store = logfile
              suite.level = logfile
              mtzout = "%s"
              suite.database = "%s"
            }
  ''' %(parent.loaded_file_name, parent.hklin, tempdir)
  )
  # reset stdout from logfname
  os.dup2(original_stdout_fd, sys.stdout.fileno())
  os.close(original_stdout_fd)

  retval = retobj.exit_code()
  errormsg = retobj.error_type() + " error, " + retobj.error_message()
  import shutil
  mtzs = retobj.get_filenames(["mtz"])
  if len(mtzs):
    parent.update_from_philstr("openfilename=" + parent.hklin) # resets all PHIL parameters
    parent.params.external_cmd = "runXtricorder" # allow this PHIL parameter to be shown
    parent.currentphil = master_phil.format(python_object=parent.params)

  else:
    raise Sorry("Could not find the mtz file from running Xtricorder")

  shutil.rmtree(tempdir)
  return group_args(tabname=tabname, logfname=logfname, retval=retval, errormsg=errormsg)
