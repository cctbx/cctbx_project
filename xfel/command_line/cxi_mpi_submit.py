from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.mpi_submit
#
# Submit a cxi.xtc_project to the cluster using bsub, after first creating an
# appropiate trial directory in the requested location and copying and modifying
# the given psana config files and cctbx.xfel phil files
#
import os, sys, shutil
from libtbx.utils import Sorry
from libtbx.phil import parse
from libtbx import easy_run

help_str = """
Use cxi.mpi_submit to submit a cxi.xtc_process job to the cluster for analyzing
diffraction data in XTC streams.  cxi.mpi_submit will create a directory for
you in the location specified by output.output_dir using the run and trial
numbers specified.  If no trial number is specified, one will be chosen auto-
matically by examining the output directory.  The output directory will be
created if it doesn't exist.

Examples:

cxi.mpi_submit input.cfg=cxi49812/thermo.cfg input.experiment=cxi49812 \\
  input.run_num=25 output.output_dir=\\
  /reg/d/psdm/cxi/cxi49812/ftc/username/results mp.nproc=100 mp.queue=psanaq

This will submit run 25 of experiment cxi49812 to the psana queue for
processing using the modules specified by cxi49812/thermo.cfg.

cxi.mpi_submit submit.phil

Use the phil parameters in submit.phil to control job submission.

Before the job is submitted, the following occurs:

Directory /reg/d/psdm/cxi/cxi49812/ftc/username/results is created if it
doesn't exist.

Directory /reg/d/psdm/cxi/cxi49812/ftc/username/results/r0025 is created if it
doesn't exist.

The run directory is searched for the next available trial ID, say 003.  That
directory is created.

Under /reg/d/psdm/cxi/cxi49812/ftc/username/results/r0025/003, the following is
created:

Directory stdout.  The log for the experiment will go here.
File psana_orig.cfg.  This is a verbatim copy of the input config file.
File psana.cfg.  This is a copy of the input config file, re-written to respect
the new paths for this trial.  Any parameters with _dirname are re-written
to include this trial directory and the associated directories are created. Any
phil files found are copied to the trial directory, renamed, and updated.  This
includes phil files included by copied phil files.

The final bsub command is saved in submit.sh for future use.
"""


phil_scope = parse('''
  input {
    cfg = None
      .type = str
      .help = "Path to psana config file"
    experiment = None
      .type = str
      .help = "Experiment identifier, e.g. cxi84914"
    run_num = None
      .type = int
      .help = "Run number or run range to process"
    trial = None
      .type = int
      .help = "Trial number for this job.  Leave blank to auto determine."
  }
  output {
    output_dir = "."
      .type = str
      .help = "Directory for output files"
  }
  mp {
    nproc = 1
      .type = int
      .help = "Number of processes"
    queue = "psanacsq"
      .type = str
      .help = "Queue to submit MPI job to"
  }
''', process_includes=True)

def copy_target(target, dest_dir, root_name):
  """ Copy a phil file to a directory, and all of its included phil files. Recurively
  copies the files and changes the includes statements to match the new file names.
  @param target Phil file to copy
  @param dest_dir Full path to directory to copy file to
  @param root_name Name to give to copied file. Included phil files will be given
  this name + _N where N is incremented for each file included in the phil file
  """
  # Each included phil file will be named root_name_N.phil where N is this number,
  # incremented for each included phil file
  num_sub_targets = 1
  f = open(os.path.join(dest_dir, root_name + ".phil"), 'w')
  for line in open(target).readlines():
    if "include" in line:
      inc_str, include_type, sub_target = line.strip().split() # Example: include file cxi-8.2.phil
      if include_type != 'file':
        raise Sorry("Include isn't a file") # FIXME look up what other values are possible here
      sub_target_root_name = "%s_%d"%(root_name, num_sub_targets)
      line = " ".join([inc_str, include_type, sub_target_root_name  + ".phil\n"])
      # recursive call to check for other included files
      copy_target(os.path.join(os.path.dirname(target), sub_target), dest_dir, sub_target_root_name)
      num_sub_targets += 1
    f.write(line)

class Script(object):
  """ Script to submit XFEL data at LCLS for processing"""
  def __init__(self):
    pass

  def run(self):
    """ Set up run folder and submit the job. """
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv or "-c" in sys.argv:
      print help_str
      print "Showing phil parameters:"
      print phil_scope.as_str(attributes_level = 2)
      return

    user_phil = []
    for arg in sys.argv[1:]:
      if (os.path.isfile(arg)):
        user_phil.append(parse(file_name=arg))
      else:
        try:
          user_phil.append(parse(arg))
        except RuntimeError, e:
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

    params = phil_scope.fetch(sources=user_phil).extract()

    assert params.input.cfg is not None
    assert params.input.experiment is not None
    assert params.input.run_num is not None

    if not os.path.exists(params.input.cfg):
      raise Sorry("Config file not found: %s"%params.input.cfg)

    print "Submitting run %d of experiment %s using config file %s"%(params.input.run_num, params.input.experiment, params.input.cfg)

    if not os.path.exists(params.output.output_dir):
       os.makedirs(params.output.output_dir)

    rundir = os.path.join(params.output.output_dir, "r%04d"%params.input.run_num)
    if not os.path.exists(rundir):
      os.mkdir(rundir)

    # If a trial number wasn't included, find the next available, up to 999 trials
    if params.input.trial is None:
      found_one = False
      for i in xrange(1000):
        trialdir = os.path.join(rundir, "%03d"%i)
        if not os.path.exists(trialdir):
          found_one = True
          break
      if found_one:
        params.input.trial = i
      else:
        raise Sorry("All trial numbers in use")
    else:
      trialdir = os.path.join(rundir, "%03d"%params.input.trial)
      if os.path.exists(trialdir):
        raise Sorry("Trial %d already in use"%params.input.trial)

    print "Using trial", params.input.trial
    os.mkdir(trialdir)

    # log file will live here
    stdoutdir = os.path.join(trialdir, "stdout")
    os.mkdir(stdoutdir)

    # make a copy of the original cfg file
    shutil.copy(params.input.cfg, os.path.join(trialdir, "psana_orig.cfg"))

    config_path = os.path.join(trialdir, "psana.cfg")

    # Re-write the config file, changing paths to be relative to the trial directory.
    # Also copy and re-write included phil files, while updating the config file to
    # the new paths
    # Note, some legacy rewrites done by cxi.lsf are not included here.
    target_num = 1
    f = open(config_path, 'w')
    for line in open(params.input.cfg).readlines():
      if "[pyana]" in line:
        raise Sorry("Pyana not supported. Check your config file.")
      if "RUN_NO" in line:
        line = line.replace("RUN_NO", str(params.input.run_num))
      if "trial_id" in line:
        key, val = line.split("=")
        line = "%s= %d\n"%(key,params.input.trial)
      elif "_dirname" in line:
        key, val = line.split("=")
        val = os.path.join(trialdir, os.path.basename(val.strip()))
        if not os.path.exists(val):
          os.mkdir(val)
        line = "%s= %s\n"%(key,val)
      elif "xtal_target" in line:
        key, val = line.split("=")
        val = val.strip()
        if not os.path.exists(val):
          raise Sorry("One of the xtal_target files in the cfg file doesn't exist: %s"%val)
        new_target = "params_%d"%target_num
        copy_target(val, trialdir, new_target)
        target_num += 1
        line = "%s= %s.phil\n"%(key,os.path.join(trialdir, new_target))
      f.write(line)
    f.close()

    command = "bsub -a mympi -n %d -o %s -q %s cxi.xtc_process input.cfg=%s input.experiment=%s input.run_num=%d"%( \
      params.mp.nproc, os.path.join(stdoutdir, "log.out"), params.mp.queue, config_path,
      params.input.experiment, params.input.run_num)

    # Write out a script for submitting this job exactly (even though it isn't used here)
    submit_path = os.path.join(trialdir, "submit.sh")
    f = open(submit_path, 'w')
    f.write("#! /bin/sh\n")
    f.write("\n")
    f.write("%s\n"%command)
    f.close()

    try:
      result = easy_run.fully_buffered(command=command).raise_if_errors()
    except Exception, e:
      if not "Warning: job being submitted without an AFS token." in str(e):
        raise e

    print "Job submitted.  Output in", trialdir

if __name__ == "__main__":
  script = Script()
  script.run()
