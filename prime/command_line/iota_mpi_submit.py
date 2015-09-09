from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME prime.iota.mpi_submit
#
# Submit a iota processing job to a cluster using bsub or qsub, after first
# creating an appropiate trial directory in the requested location and copying
# and modifying the given iota config files and cctbx.xfel phil files
#

import os, sys, shutil
from libtbx.utils import Sorry
from libtbx.phil import parse
from libtbx import easy_run

help_str = """
Use prime.iota.mpi_submit to submit a prime.linear_iota job to the cluster for
analyzing stills diffraction data in single-file format.  prime.iota.mpi_submit
will create a directory for you in the location specified by output.output_dir
using the run and trial numbers specified.  If no trial number is specified, one
will be chosen automatically by examining the output directory.  The output
directory will be created if it doesn't exist.

Examples:

prime.iota.mpi_submit input.param=iota.param output.output_dir=\\
  /reg/d/psdm/cxi/cxi49812/ftc/username/results mp.nproc=100 mp.queue=psanaq

This will submit the data indicated in the givin iota.param file to the psana
queue for processing using iota. The user may optionally provide input.run_num=X
if the input data is organized into runs, say from an XFEL experiment.

prime.iota.mpi_submit submit.phil

Use the phil parameters in submit.phil to control job submission.

Before the job is submitted, the following occurs:

Directory /reg/d/psdm/cxi/cxi49812/ftc/username/results is created if it
doesn't exist.

If a run number is provided, say run 25, directory\\
/reg/d/psdm/cxi/cxi49812/ftc/username/results/r0025 is created if it
doesn't exist.

The run directory is searched for the next available trial ID, say 003.  That
directory is created.

Under /reg/d/psdm/cxi/cxi49812/ftc/username/results/r0025/003, the following is
created:

Directory stdout.  The log for the experiment will go here.
File iota_orig.param.  This is a verbatim copy of the input param file.
File iota.param.  This is a copy of the input config file, re-written to respect
the new paths for this trial.  Any phil files found are copied to the trial
directory, renamed, and updated.  This includes phil files included by copied phil
files.

The final bsub or qsub command is saved in submit.sh for future use.
"""


phil_scope = parse('''
  dry_run = False
    .type = bool
    .help = If True, the program will create the trial directory but not submit the job, \
            and will show the command that would have been executed.
  input {
    param = None
      .type = str
      .help = Path to iota param file. User can include RUN_NO, RUN_STR and NPROC, which \
              will be replaced with the appropiate values when the param file is copied \
              to the trial directory.
    run_num = None
      .type = int
      .help = Run number to process
    trial = None
      .type = int
      .help = Trial number for this job.  Leave blank to auto determine.
    dispatcher = prime.linear_iota
      .type = str
      .help = Which program to run. prime.linear_iota is the mpi enabled iota
  }
  output {
    output_dir = "."
      .type = str
      .help = Directory for output files
    split_logs = False
      .type = bool
      .help = Option to split error and log files into separate per process
  }
  mp {
    method = *mpi sge pbs custom
      .type = choice
      .help = Muliprocessing method
    nproc = 1
      .type = int
      .help = Number of processes
    queue = "psanaq"
      .type = str
      .help = Queue to submit MPI job to
    sge {
      memory = 4g
        .type = str
        .help = How much memory to request for an SGE job
    }
    pbs {
      env_script = None
        .type = str
        .help = Path to bash script with extra environment settings. Ran before PSDM \
                (if available) and cctbx are sourced.
      walltime = "00:00:30"
        .type = str
        .help = Maximum time this job will be scheduled to run for.
    }
    custom {
      extra_args = ""
        .type = str
        .help = Extra arguments to qsub as needed
      submit_template = None
        .type = str
        .help = Submission script for qsub. The script will be copied to the trial \
                directory and modified. There should be one instance of the string \
                <command> (including the <> brackets) in the template script, which \
                will be replaced with the processing command. <queue> and <nproc> \
                will similarly be replaced.
    }
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
  """ Script to submit iota data for processing"""
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

    assert params.input.param is not None

    if not os.path.exists(params.input.param):
      raise Sorry("Param file not found: %s"%params.input.param)

    if params.input.run_num is None:
      print "Submitting iota job using param file %s"%params.input.param
    else:
      print "Submitting run %d as iota job using param file %s"%(params.input.run_num, params.input.param)

    if not os.path.exists(params.output.output_dir):
       os.makedirs(params.output.output_dir)

    if params.input.run_num is None:
      rundir = params.output.output_dir
    else:
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
    if params.output.split_logs:# test parameter for split_log then open and close log file and loop over nprocs
      for i in xrange(params.mp.nproc):
        error_files = os.path.join(stdoutdir,"error_rank%04d.out"%i)
        log_files = os.path.join(stdoutdir,"log_rank%04d.out"%i)
        open(log_files,'a').close()
        open(error_files,'a').close()
    # make a copy of the original params file
    shutil.copy(params.input.param, os.path.join(trialdir, "iota_orig.param"))

    param_path = os.path.join(trialdir, "iota.param")

    # Re-write the config file, changing paths to be relative to the trial directory.
    # Also copy and re-write included phil files, while updating the config file to
    # the new paths
    # Note, some legacy rewrites done by cxi.lsf are not included here.
    target_num = 1
    f = open(param_path, 'w')
    for line in open(params.input.param).readlines():
      if "RUN_NO" in line:
        line = line.replace("RUN_NO", str(params.input.run_num))
      if "RUN_STR" in line:
        line = line.replace("RUN_STR", "r%04d"%(params.input.run_num))
      if "NPROC" in line:
        line = line.replace("NPROC", "%d"%(params.mp.nproc))
      if "trial_id" in line:
        key, val = line.split("=")
        line = "%s= %d\n"%(key,params.input.trial)
      elif "target" in line and line.split('=')[0].strip() == "target":
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

    # Write out a script for submitting this job and submit it
    submit_path = os.path.join(trialdir, "submit.sh")

    if params.mp.method == "mpi":
      command_process = "bsub -a mympi -n %d -o %s -q %s %s %s --mpi process"%(
        params.mp.nproc, os.path.join(stdoutdir, "log.out"), params.mp.queue, params.input.dispatcher, param_path)

      f = open(submit_path, 'w')
      f.write("#! /bin/sh\n")
      f.write("\n")
      f.write("%s\n"%command_process)
      f.close()
    elif params.mp.method == "sge":
      raise NotImplementedError()
    elif params.mp.method == "pbs":
      # Write out a script for submitting this job and submit it
      submit_path = os.path.join(trialdir, "submit.csh")

      command_process = "qsub -o %s -d %s %s"%(os.path.join(stdoutdir, "log.out"), trialdir, submit_path)

      import libtbx.load_env
      setpaths_path = os.path.join(abs(libtbx.env.build_path), "setpaths.sh")
      psdmenv_path = os.path.join(os.environ.get("SIT_ROOT", ""), "etc", "ana_env.sh")
      assert os.path.exists(setpaths_path)

      f = open(submit_path, 'w')
      f.write("#!/bin/bash\n")
      f.write("#PBS -q %s\n"%params.mp.queue)
      f.write("#PBS -l mppwidth=%d\n"%params.mp.nproc)
      f.write("#PBS -l walltime=%s\n"%params.mp.pbs.walltime)
      f.write("#PBS -N run%d\n"%params.input.run_num)
      f.write("#PBS -j oe\n")
      f.write("\n")
      f.write("cd $PBS_O_WORKDIR\n")
      f.write("\n")
      if params.mp.pbs.env_script is not None:
        f.write("source %s\n"%params.mp.pbs.env_script)
      f.write("\n")
      if os.path.exists(psdmenv_path):
        f.write("source %s\n"%psdmenv_path)
      f.write("source %s\n"%setpaths_path)
      if os.path.exists(psdmenv_path):
        f.write("sit_setup\n")

      f.write("aprun -n %d %s %s --mpi process\n"%( \
        params.mp.nproc, params.input.dispatcher, param_path))
      f.close()
    elif params.mp.method == "custom":
      if not os.path.exists(params.mp.custom.submit_template):
        raise Sorry("Custom submission template file not found: %s"%params.mp.custom.submit_template)

      command_process = "qsub -o %s %s %s"%(os.path.join(stdoutdir, "log.out"), params.mp.custom.extra_args, submit_path)

      processing_command = "%s %s --mpi process\n"%( \
         params.input.dispatcher, param_path)

      f = open(submit_path, 'w')
      for line in open(params.mp.custom.submit_template).readlines():
        if "<command>" in line:
          line = line.replace("<command>", processing_command)
        if "<queue>" in line:
          line = line.replace("<queue>", params.mp.queue)
        if "<nproc>" in line:
          line = line.replace("<nproc>", str(params.mp.nproc))

        f.write(line)
      f.close()

    else:
      raise Sorry("Multiprocessing method %s not recognized"%params.mp.method)

    print "Performing iota initialization..."
    cwd = os.getcwd()
    os.chdir(trialdir)
    command_init = "%s %s --mpi init"%(params.input.dispatcher, param_path)
    result = easy_run.fully_buffered(command=command_init).raise_if_errors()
    result.show_stdout()

    if params.dry_run:
      print "Dry run: job not submitted. Trial directory created here:", trialdir
      print "Execute this command to submit the job:"
      print command_process
    else:
      try:
        result = easy_run.fully_buffered(command=command_process).raise_if_errors()
      except Exception, e:
        if not "Warning: job being submitted without an AFS token." in str(e):
          raise e

      print "Job submitted.  Output in", trialdir

    os.chdir(cwd)
    print "When job has concluded, analyze with these commands:"
    print "cd %s"%trialdir
    print "%s %s --mpi analyze"%(params.input.dispatcher, param_path)

if __name__ == "__main__":
  script = Script()
  script.run()
