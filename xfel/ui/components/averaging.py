from __future__ import absolute_import, division, print_function

from xfel.util.mp import get_lsf_submit_command
import os


class AveragingCommand(get_lsf_submit_command):
  def __init__(self, run, db_params, raw=False):
    # command setup
    super(AveragingCommand, self).__init__("cxi.mpi_average", None, None, db_params.mp, job_name="average")
    # run parameters
    run_no = int(run.run)
    self.args.append("-r %d" % run_no)
    # current settings and experiment parameters
    nproc = int(db_params.mp.nproc)
    self.command = "-n %d %s" % (nproc, self.command) # to make sure nproc is passed to mpirun
    expt = str(db_params.facility.lcls.experiment)
    self.args.append("-x %s" % expt)
    outdir = os.path.join(str(db_params.output_folder), "averages")
    if not os.path.exists(outdir):
      os.mkdir(outdir)
    self.args.append("-o %s" % outdir)
    logdir = os.path.join(outdir, "logs")
    if not os.path.exists(logdir):
      os.mkdir(logdir)
    self.stdoutdir = logdir
    self.submit_path = os.path.join(logdir, "r%04d.sh" % run_no)
    self.log_name = "r%04d.log" % run_no
    ffb = db_params.facility.lcls.use_ffb
    if ffb:
      self.args.append("-f")
    if raw:
      self.args.append("-R")
    # rungroup parameters (using the most recent, active rungroup
    rungroups = run.get_rungroups()
    rungroup_ids = [rg.id for rg in rungroups]
    rungroup_obj = rungroups[rungroup_ids.index(max(rungroup_ids))]
    address = rungroup_obj.detector_address
    self.args.append("-a %s" % address)
    beamx = rungroup_obj.beamx
    if beamx:
      self.args.append("-X %.1f" % beamx)
    beamy = rungroup_obj.beamy
    if beamy:
      self.args.append("-Y %.1f" % beamy)
    distance = rungroup_obj.detz_parameter
    self.args.append("-d %.3f" % distance)
    config = rungroup_obj.config_str
    if config:
      self.args.append("-c %s" % config)
    pickle = rungroup_obj.format == "pickle"
    if pickle:
      self.args.append("-p")
    binning = rungroup_obj.binning
    if binning:
      self.args.append("-b %d" % binning)
    self.args.append("-v")
  def customize_for_method(self):
    self.submit_head = "bsub"
    self.command = "mpirun %s" % self.command
