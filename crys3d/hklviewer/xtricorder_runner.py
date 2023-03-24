from __future__ import absolute_import, division, print_function
from libtbx import group_args
import os, os.path, glob, threading

def external_cmd(parent, master_phil, firstpart, tempdir):
  from phasertng.scripts import xtricorder

  tabname = "Xtricorder"
  (retobj) = xtricorder.xtricorder(
  r'''phasertng {
              hklin.filename = "%s"
              reflections.wavelength = 1.0
              suite.mute = True
              suite.store = logfile
              suite.level = logfile
              suite.database = "%s"
              commensurate.patterson.percent = 10
            }
  ''' %(parent.loaded_file_name, tempdir)
  )
  retval = retobj.exit_code()
  errormsg = retobj.error_type() + " error, " + retobj.error_message()
  import shutil, glob
  mtzs = retobj.get_filenames(["mtz"])
  if len(mtzs):
    xtricordermtz = mtzs[-1]
    parent.hklin =  firstpart + "_xtricorder.mtz"
    shutil.copyfile( xtricordermtz, parent.hklin ) # copy the last file only
    parent.update_from_philstr("openfilename=" + parent.hklin)

  logs = glob.glob(tempdir + "/**/*.logfile.log", recursive=True)
  timesortedlogs = sorted( [ (p, os.path.getmtime(p) )   for p in logs ], key=lambda e: e[1] )
  mstr = ''
  for fname, t in timesortedlogs:
    with open(fname, 'r') as f:
      mstr += f.read() + '\\n'
  # The name of logfile and tab should be present in ldic after running exec().
  # cctbx.python sends this back to HKLviewer from HKLViewFrame.run_external_cmd()
  logfname = firstpart + "_xtricorder.log"
  with open(logfname, 'w') as f:
    f.write(mstr)
  if len(mtzs) == 0:
    raise Sorry("Could not find the mtz file from running Xtricorder")

  shutil.rmtree(tempdir)
  return group_args(tabname=tabname, logfname=logfname, retval=retval, errormsg=errormsg)
