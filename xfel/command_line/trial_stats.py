from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.trial_stats
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import libtbx.phil
from libtbx.utils import Usage, Sorry
import sys

master_phil = libtbx.phil.parse("""
  trial_id = None
    .type = int
  hit_cutoff = 16
    .type = int
  run_start = None
    .type = int
  run_end = None
    .type = int
  db {
    host = None
      .type = str
    name = None
      .type = str
    table_name = None
      .type = str
    user = None
      .type = str
    password = None
      .type = str
    tags = None
      .type = str
  }
""")

def run (args) :
  try:
    from cxi_xdr_xes.cftbx.cspad_ana import db as db
  except ImportError:
    raise Sorry("Trial logging not supported for this installation. Conact the developers for access.")

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
      except ValueError as e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = master_phil.fetch(sources=user_phil).extract()
  if (params.trial_id is None) :
    master_phil.show()
    raise Usage("trial_id must be defined (either trial_id=XXX, or the integer "+
      "ID alone).")
  assert (params.hit_cutoff is not None) and (params.hit_cutoff > 0)

  extra_cmd = ""
  if params.run_start is not None:
    extra_cmd += "AND run >= %d" % params.run_start
  if params.run_end is not None:
    extra_cmd += "AND run <= %d" % params.run_end

  dbobj = db.dbconnect(host=params.db.host, db=params.db.name, username=params.db.user, password=params.db.password)


  cursor = dbobj.cursor()
  cmd = "SELECT DISTINCT(run) FROM %s WHERE trial = %%s %s ORDER BY run"%(params.db.table_name, extra_cmd)
  cursor.execute(cmd, params.trial_id)

  frames_total = 0
  hits_total = 0
  indexed_total = 0

  for runId in cursor.fetchall():
    run = int(runId[0])
    cmd = "SELECT id, eventstamp, hitcount, distance, sifoil, wavelength, indexed FROM %s \
        WHERE trial = %s AND run = %s"
    if params.db.tags is not None:
      for tag in params.db.tags.split(','):
        cmd += """ AND tags LIKE "%%{0}%%" """.format(tag)
    cursor.execute(cmd%(params.db.table_name,params.trial_id,run))

    numframes = numhits = numindexed = 0
    for id, eventstamp, hitcount, distance, sifoil, wavelength, indexed in cursor.fetchall():
      numframes +=1
      if hitcount >= params.hit_cutoff:
        numhits += 1
      if indexed:
        numindexed += 1

    if numhits == 0:
      hitrate = 0
    else:
      hitrate = 100*numhits/numframes
    if numindexed == 0:
      indexingrate = 0
    else:
      indexingrate = 100*numindexed/numframes

    print("Run: %3d, number of hits: %6d, number of frames: %6d, hitrate: %4.1f%%. Number indexed: %6d (%4.1f%%)"%(run,numhits,numframes,hitrate,numindexed,indexingrate))
    frames_total += numframes
    hits_total += numhits
    indexed_total += numindexed

  if hits_total == 0:
    hitrate = 0
  else:
    hitrate = 100*hits_total/frames_total
  if indexed_total == 0:
    indexingrate = 0
  else:
    indexingrate = 100*indexed_total/frames_total

  print("Totals: frames: %d, hits: %d (%4.1f%%), indexed: %d (%4.1f%%)"%(frames_total,hits_total,hitrate,indexed_total,indexingrate))
  dbobj.close()

if (__name__ == "__main__") :
  run(sys.argv[1:])
