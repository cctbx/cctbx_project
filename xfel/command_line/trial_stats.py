from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.trial_stats
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from cxi_xdr_xes.cftbx.cspad_ana import db as db
import libtbx.phil
from libtbx.utils import Usage, Sorry
import sys

master_phil = libtbx.phil.parse("""
  trial_id = None
    .type = int
  hit_cutoff = 16
    .type = int
""")

def run (args) :
  user_phil = []
  # TODO: replace this stuff with iotbx.phil.process_command_line_with_files
  # as soon as I can safely modify it
  for arg in args :
    #if (os.path.isdir(arg)) :
      #user_phil.append(libtbx.phil.parse("""status_dir=\"%s\"""" % arg))
    #elif (not "=" in arg) :
    if (not "=" in arg) :
      try :
        user_phil.append(libtbx.phil.parse("""trial_id=%d""" % int(arg)))
      except ValueError, e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = master_phil.fetch(sources=user_phil).extract()
  if (params.trial_id is None) :
    master_phil.show()
    raise Usage("trial_id must be defined (either trial_id=XXX, or the integer "+
      "ID alone).")
  assert (params.hit_cutoff is not None) and (params.hit_cutoff > 0)

  dbobj = db.dbconnect()


  cursor = dbobj.cursor()
  cmd = "SELECT DISTINCT(run) FROM %s WHERE trial = %s ORDER BY run"
  cursor.execute(cmd, (db.table_name,params.trial_id)


  for runId in cursor.fetchall():
    run = int(runId[0])
    cursor.execute("SELECT id, eventstamp, hitcount, distance, sifoil, wavelength FROM %s \
        WHERE trial = %s AND run = %s"%(db.table_name,params.trial_id,run))

    numframes = numhits = 0
    for id, eventstamp, hitcount, distance, sifoil, wavelength in cursor.fetchall():
      numframes +=1
      if hitcount >= params.hit_cutoff:
        numhits += 1

    print "Run: %3d, number of hits: %6d, number of frames: %6d, hitrate: %3.1f%%"%(run,numhits,numframes,100*numhits/numframes)

  dbobj.close()

if (__name__ == "__main__") :
  run(sys.argv[1:])
