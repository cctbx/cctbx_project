from __future__ import division

#-----------------------------------------------------------------------
# Monitor a directory for streams and submit them
#-----------------------------------------------------------------------

import libtbx.phil
#from libtbx.utils import Usage, Sorry
from cxi_xdr_xes.cftbx.cspad_ana import db as db
import os
import sys
import time
import libtbx

master_phil = libtbx.phil.parse("""
  xtc_dir = None
    .type = path
  output_dir = None
    .type = path
  trial_id = None
    .type = int
  stream_count = 5
    .type = int
  num_procs = 12
    .type = int
  queue = 'psfehmpiq'
    .type = str
  start_run = None
    .type = int
  config_file = 'onlyhitfind.cfg'
    .type = str
  experiment = 'cxi78513'
    .type = str
  submit_as_group = True
    .type = bool
""")

submitted_runs = []

class _run:
  def __init__(self, id):
    self.id = id
    self.files = []

  def __eq__(self, other): return other == self.id
  def __ne__(self, other): return not __eq__(self,other)

def match_runs(dir):
  runs = []
  files = os.listdir(dir)
  for file in files:
    #print file.split('.')[-1] != 'xtc'
    if file.split('.')[-1] != 'xtc':
      continue
    r = s = c = None
    #print file.split("-")
    for str in file.split("-"):
      if 'inprogress' in str:
        break
      if 'c' in str:
        #print "yarp: %s"%str.split(".")[0].strip('c')
        c = int(str.split(".")[0].strip('c'))
      if 'r' in str:
        r = int(str.strip('r'))
      elif 's' in str:
        s = int(str.strip('s'))
    #print r, s, c
    if r is not None and s is not None and c is not None and c == 0:
      foundIt = False
      for rn in runs:
        if rn.id == r:
          foundIt = True
          if not file in rn.files:
            rn.files.append(file)
      if not foundIt:
        rn = _run(r)
        rn.files.append(file)
        runs.append(rn)
  return runs

def run (args) :
  path = os.getcwd()
  if os.path.basename(path) != "myrelease":
    raise Sorry("You must run this script from within your pyana myrelease directory.")

  user_phil = []
  # TODO: replace this stuff with iotbx.phil.process_command_line_with_files
  # as soon as I can safely modify it
  for arg in args :
    if (os.path.isdir(arg)) :
      user_phil.append(libtbx.phil.parse("""xtc_dir=\"%s\"""" % arg))
    elif (not "=" in arg) :
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
  if (params.xtc_dir is None) :
     master_phil.show()
     raise Usage("xtc_dir must be defined!")
  elif (not os.path.isdir(params.xtc_dir)) :
     raise Sorry("%s does not exist or is not a directory!" % params.xtc_dir)

  if (params.output_dir is None) :
     master_phil.show()
     raise Usage("output_dir must be defined!")
  elif (not os.path.isdir(params.output_dir)) :
     raise Sorry("%s does not exist or is not a directory!" % params.output_dir)

  if (params.config_file is None) :
     master_phil.show()
     raise Usage("output_dir must be defined!")
  elif (not os.path.isfile(params.config_file)) :
     raise Sorry("%s does not exist or is not a file!" % params.output_dir)

  assert (params.stream_count is not None) and (params.stream_count > 0)
  assert (params.num_procs is not None) and (params.stream_count > 0)
  assert (params.queue is not None)
  assert (params.experiment is not None)
  assert (params.submit_as_group is not None)

  submitted_runs = []
  submitted_files = [] # used in single stream submit mode

  print "Note, it is not recommended that you run this program while you have new jobs pending."
  print "Starting indefinite loop to scan directory '%s'"%params.xtc_dir
  #rs = match_runs(params.xtc_dir)
  #return

  try:
    while(True):
      rs = match_runs(params.xtc_dir)
      add_runs = []
      for r in rs:
        if not (params.start_run is not None and r.id < params.start_run):
          if not db.run_in_trial(r.id, params.trial_id):
            doit = True
            for test in submitted_runs:
              if test.id == r.id:
                doit = False
                break
            if doit: add_runs.append(r)

      submitted_a_run = False
      if len(add_runs) > 0:
        for r in add_runs:
	  if params.submit_as_group:
            if len(r.files) != params.stream_count:
              print "Waiting to queue run %s.  %s/%s streams ready."% \
                (r.id,len(r.files),params.stream_count)
              continue

            print "Preparing to queue run %s into trial %s"%(r.id,params.trial_id)
            cmd = "./lsf.sh -c %s -p %s -x %s -o %s -t %s -r %s -q %s"%(params.config_file,params.num_procs,params.experiment,
              params.output_dir,params.trial_id,r.id,params.queue)

            print "Command to execute: %s"%cmd
            os.system(cmd)
            print "Run %s queued."%r.id

            submitted_a_run = True
            submitted_runs.append(r)
	  else:
	    for f in r.files:
	      if not f in submitted_files:
	        s = None
                for str in f.split('-'):
		  if 's' in str:
                    s = int(str.strip('s'))
		if s is None:
		  print "Couldn't get the stream number from file %s"%os.path.basename(f)

		print "Preparing to queue stream %s into trial %s"%(os.path.basename(f),params.trial_id)
                cmd = "./single_lsf.sh -c %s -p %s -x %s -o %s -t %s -r %s -q %s -s %s"%(params.config_file,params.num_procs,params.experiment,
                  params.output_dir,params.trial_id,r.id,params.queue,s)

                print "Command to execute: %s"%cmd
                os.system(cmd)
                print "Stream %s queued."%r.id

		submitted_a_run = True
		submitted_files.append(f)

      if not submitted_a_run:
        print "No new data... sleepy..."

      time.sleep(10)
  except KeyboardInterrupt:
    pass



if (__name__ == "__main__") :
  run(sys.argv[1:])
