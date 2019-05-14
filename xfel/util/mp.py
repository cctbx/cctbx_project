from __future__ import absolute_import, division, print_function
#
# Handle multiprocessing with any of the implemented methods so that this step
# is abstracted away from the use case (e.g. cxi_mpi_submit).
#
from libtbx.utils import Sorry
import os

mp_phil_str = '''
  mp {
    method = *lsf sge pbs shifter custom
      .type = choice
      .help = Computing environment
    use_mpi = True
      .type = bool
      .help = Use mpi multiprocessing
    nproc = 1
      .type = int
      .help = Number of processes total (== nnodes x nproc_per_node). \
              If two of the three params (nproc, nnodes, nproc_per_node) are \
              specified, the last will be determined by modular arithmetic. \
              If all three are specified, nnodes is ignored. nproc alone is \
              sufficient for most methods.
    nnodes = 1
      .type = int
      .help = Number of nodes to request
    nproc_per_node = 1
      .type = int
      .help = Number of processes to allocate per node
    queue = None
      .type = str
      .help = Queue to submit multiprocessing job to (optional for some methods)
    memory = None
      .type = int
      .help = Memory (in MB) to allocate for a job (optional)
    wall_time = None
      .type = int
      .help = Wall time limit (in minutes) to impose (optional)
    extra_options = None
      .type = str
      .multiple = True
      .help = Any other options to be included in the job submission command
    extra_args = None
      .type = str
      .multiple = True
      .help = Any other arguments to the main command
    env_script = None
      .type = str
      .multiple = True
      .help = Path to script sourcing a particular environment (optional)
    shifter {
      submit_command = "sbatch "
        .type = str
        .help = Command used to run the zero-level script sbatch.sh.
      sbatch_script_template = None
        .type = path
        .help = Script template to be run with sbatch. The script will be copied \
                to the trial directory as sbatch.sh and modified. Must contain \
                exactly one srun command of the format srun [...] <srun_script> \
                (including the <> brackets). May contain additional srun commands. \
                May also contain substitutables <walltime>, <partition>, <nnodes> \
                and <nproc>.
      srun_script_template = None
        .type = path
        .help = Script template to be run with srun. The script will be copied \
                to the trial directory as srun.sh and modified. Must contain \
                exactly one instance of the string <command> (including the <> \
                brackets) after setting up the necessary environment.
      partition = regular
        .type = str
        .help = Partition to run jobs on, e.g. regular or debug.
      jobname = LCLS_EXP
        .type = str
        .help = Jobname
      nnodes = 1
        .type = int
        .help = Number of nodes to request with sbatch -N <nnodes>.
      nproc = 32
        .type = int
        .help = Number of processors (total) to request with srun -n <nproc>.
    }
    custom {
      submit_command_template = None
        .type = str
        .help = Job submission command template. There should be one instance of the \
                string <script> (including the <> brackets) in the command template. \
                Fields <queue>, <nproc>, <memory>, <walltime>, <outfile>, <errfile>, \
                <envscripts>, and <args> will similarly be replaced if present.
      submit_script_template = None
        .type = path
        .help = Submission script template. The script will be copied to the trial \
                directory and modified. There should be one instance of the string \
                <command> (including the <> brackets) in the template script, which \
                will be replaced with the processing command. <queue>, <nproc>, \
                <memory>, <walltime>, <outfile>, <errfile>, <envscripts>, and <args> \
                will similarly be replaced if present.
      wall_time_string = None
        .type = str
        .help = Fully formatted wall time limit (e.g. 00:30:00). For custom computing \
                environments, mp.wall_time is ignored because the format is unknown.
    }
    encapsulate_submit_script = True
      .type = bool
      .help = Encapsulate the submission command itself in another script containing \
              the qsub or bsub job submission command.
  }
'''

class get_submit_command(object):
  def __init__(self, command, submit_path, stdoutdir, params,
               log_name="log.out", err_name="log.err", job_name=None):
    """ Get a submit command for the various compute environments
    @param command Any command line program and its arguments
    @param submit_path Submit script will be written here
    @param stdoutdir Log file will be created in this directory
    @param params Multiprocessing phil params (see mp_phil_scope)
    @param log_name Filename for stdout (optional).
    @param err_name Filename for stderr (if None, combined with the stdout).
    @param job_name For applicable queueing systems, identifier for the job (optional).
    """
    self.shell_path = "/bin/sh"
    self.source_env_scripts = []
    self.options_inside_submit_script = []
    self.submit_head = "qsub"
    self.submit_path = submit_path
    self.stdoutdir = stdoutdir
    self.log_name = log_name
    self.err_name = err_name
    self.params = params
    self.job_name = job_name
    self.command = command
    self.options = []
    self.args = []

  def customize_for_method(self):
    pass

  def eval_params(self):
    pass

  def substitute(self, template, marker, value):
    if marker in template:
      if value is None:
        raise Sorry("No value found for %s" % marker)
      return template.replace(marker, value)
    else:
      return template

  def make_executable(self, file):
    import stat
    st = os.stat(file)
    os.chmod(file, st.st_mode | stat.S_IXUSR)

  def write_script(self):
    command_str = " ".join([self.command] + self.args)
    f = open(self.submit_path, 'wb')
    f.write("#! %s\n" % self.shell_path)
    for line in self.options_inside_submit_script:
      f.write("%s\n" % line)
    for line in self.source_env_scripts:
      f.write("%s\n" % line)
    f.write("\n")
    f.write("%s\n" % command_str)
    f.close()
    self.make_executable(self.submit_path)

  def generate_submit_command(self):
    return " ".join([self.submit_head] + self.options + [self.submit_path])

  def encapsulate_submit(self):
    path, ext = os.path.splitext(self.submit_path)
    encapsulate_path = path + "_submit" + ext
    f = open(encapsulate_path, 'wb')
    f.write("#! /bin/%s\n\n" % ext[1:])
    f.write(self.generate_submit_command())
    f.write("\n")

  def __call__(self):
    self.customize_for_method()
    self.eval_params()
    self.write_script()
    if self.params.encapsulate_submit_script:
      self.encapsulate_submit()
    return self.generate_submit_command()

class get_lsf_submit_command(get_submit_command):

  def customize_for_method(self):
    self.submit_head = "bsub"
    if self.params.use_mpi:
      self.command = "mpirun %s mp.method=mpi" % self.command

  def eval_params(self):
    # -n <nproc>
    nproc_str = "-n %d" % self.params.nproc
    self.options.append(nproc_str)

    # -o <outfile>
    out_str = "-o %s" % os.path.join(self.stdoutdir, self.log_name)
    self.options.append(out_str)

    # -e <errfile> (optional)
    if self.err_name is not None:
      err_str = "-e %s" % os.path.join(self.stdoutdir, self.err_name)
      self.options.append(err_str)

    # -q <queue> (optional on LSF)
    if self.params.queue is not None:
      queue_str = "-q %s" % self.params.queue
      self.options.append(queue_str)

    # -W <wall_time_limit> (optional)
    if self.params.wall_time is not None:
      hours = self.params.wall_time // 60
      minutes = self.params.wall_time % 60
      wt_str = "-W %2d:%02d" % (hours, minutes)
      self.options.append(wt_str)

    # -R "rusage[mem=<memory_requested>]" (optional)
    if self.params.memory is not None:
      memory_str = "-R \"rusage[mem=%d]\"" % self.params.memory
      self.options.append(memory_str)

    # <extra_options> (optional, preceding the command)
    for cmd in self.params.extra_options:
      self.options.append(cmd)

    # source </path/to/env.sh> (optional)
    for env in self.params.env_script:
      env_str = "source %s\n" % env
      self.source_env_scripts.append(env_str)

    # <args> (optional, following the command)
    for arg in self.params.extra_args:
      self.args.append(arg)

class get_sge_submit_command(get_submit_command):

  def customize_for_method(self):
    self.shell_path += " -q"
    self.options.append("-cwd")
    self.options.append("mp.method=sge")

  def eval_params(self):
    # -t 1-<nproc>
    if self.params.nproc > 1:
      nproc_str = "-t 1-%d" % self.params.nproc
      self.options.append(nproc_str)

    # -o <outfile>
    out_str = "-o %s" % os.path.join(self.stdoutdir, self.log_name)
    self.options.append(out_str)

    # -j [y/n] -e <errfile> (optional)
    if self.err_name is not None:
      err_str = "-j n -e %s" % os.path.join(self.stdoutdir, self.err_name)
      self.options.append(err_str)
    else:
      self.options.append("-j y")

    # -q <queue>
    if self.params.queue is None:
      raise Sorry("Queue not specified.")
    queue_str = "-q %s" % self.params.queue
    self.options.append(queue_str)

    # -l h_rt=<wall_time_limit> (optional)
    if self.params.wall_time is not None:
      hours = self.params.wall_time // 60
      minutes = self.params.wall_time % 60
      wt_str = "-l h_rt=%02d:%02d:00" % (hours, minutes)
      self.options.append(wt_str)

    # -l mem_free=<memory_requested> (optional)
    if self.params.memory is not None:
      memory_str = "-l mem_free=%dM" % self.params.memory
      self.options.append(memory_str)

    # -N <job_name>
    if self.job_name is not None:
      name_str = "-N %s" % self.job_name
      self.options.append(name_str)

    # <extra_options> (optional, preceding the command)
    for cmd in self.params.extra_options:
      self.options.append(cmd)

    # source </path/to/env.sh> (optional)
    for env in self.params.env_script:
      env_str = "source %s\n" % env
      self.source_env_scripts.append(env_str)

    # <args> (optional, following the command)
    for arg in self.params.extra_args:
      self.args.append(arg)

class get_pbs_submit_command(get_submit_command):

  def customize_for_method(self):
    if (self.params.nnodes > 1) or (self.params.nproc_per_node > 1):
      self.params.nproc = self.params.nnodes * self.params.nproc_per_node
    if self.params.use_mpi:
      self.command = "mpiexec --hostfile $PBS_NODEFILE %s mp.method=mpi" % (self.command)

  def eval_params(self):

    # # -t 1-<nproc> # deprecated
    # if self.params.nproc > 1:
    #   nproc_str = "#PBS -l mppwidth=%d" % self.params.nproc
    #   self.options_inside_submit_script.append(nproc_str)

    # -l nodes=<nnodes>:ppn=<procs_per_node>
    if max(self.params.nproc, self.params.nproc_per_node, self.params.nnodes) > 1:
      # If specified, nproc overrides procs_per_node and procs_per_node overrides
      # nnodes. One process per node is requested if only nproc is specified.
      if self.params.nproc > 1:
        import math
        if self.params.nproc <= self.params.nproc_per_node:
          procs_per_node = self.params.nproc
          nnodes = 1
        elif self.params.nproc_per_node > 1:
          procs_per_node = self.params.nproc_per_node
          nnodes = int(math.ceil(self.params.nproc/procs_per_node))
        elif self.params.nnodes > 1:
          procs_per_node = int(math.ceil(self.params.nproc/self.params.nnodes))
          nnodes = self.params.nnodes
        else: # insufficient information; allocate 1 proc per node
          procs_per_node = 1
          nnodes = self.params.nproc
      else:
        procs_per_node = self.params.nproc_per_node
        nnodes = self.params.nnodes
      nproc_str = "#PBS -l nodes=%d:ppn=%d" % (nnodes, procs_per_node)
      self.options_inside_submit_script.append(nproc_str)

    # -o <outfile>
    out_str = "#PBS -o %s" % os.path.join(self.stdoutdir, self.log_name)
    self.options_inside_submit_script.append(out_str)

    # [-j oe/-e <errfile>] (optional)
    if self.err_name is not None:
      err_str = "#PBS -e %s" % os.path.join(self.stdoutdir, self.err_name)
      self.options_inside_submit_script.append(err_str)
    else:
      self.options_inside_submit_script.append("#PBS -j oe")

    # -q <queue>
    if self.params.queue is None:
      raise Sorry("Queue not specified.")
    queue_str = "#PBS -q %s" % self.params.queue
    self.options_inside_submit_script.append(queue_str)

    # -l walltime=<wall_time_limit> (optional)
    if self.params.wall_time is not None:
      hours = self.params.wall_time // 60
      minutes = self.params.wall_time % 60
      wt_str = "#PBS -l walltime=%2d:%02d:00" % (hours, minutes)
      self.options_inside_submit_script.append(wt_str)

    # -l mem_free=<memory_requested> (optional)
    if self.params.memory is not None:
      memory_str = "#PBS -l mem=%dmb" % self.params.memory
      self.options_inside_submit_script.append(memory_str)

    # -N <job_name>
    if self.job_name is not None:
      name_str = "#PBS -N %s" % self.job_name
      self.options_inside_submit_script.append(name_str)

    # <extra_options> (optional, preceding the command)
    for cmd in self.params.extra_options:
      cmd_str = "#PBS %s" % cmd
      self.options_inside_submit_script.append(cmd_str)

    # source </path/to/env.sh> (optional)
    for env in self.params.env_script:
      env_str = "source %s\n" % env
      self.source_env_scripts.append(env_str)

    # <args> (optional, following the command)
    for arg in self.params.extra_args:
      self.args.append(arg)

class get_shifter_submit_command(get_submit_command):

  def customize_for_method(self):
    # template for sbatch.sh
    self.sbatch_template = self.params.shifter.sbatch_script_template
    if self.sbatch_template is None:
      raise Sorry("sbatch script template required for shifter")
    sb = open(self.sbatch_template, "rb")
    self.sbatch_contents = sb.read()
    self.destination = os.path.dirname(self.submit_path)
    sb.close()
    self.sbatch_path = os.path.join(self.destination, "sbatch.sh")

    # template for srun.sh
    self.srun_template = self.params.shifter.srun_script_template
    if self.srun_template is None:
      raise Sorry("srun script template required for shifter")
    sr = open(self.srun_template, "rb")
    self.srun_contents = sr.read()
    sr.close()
    self.srun_path = os.path.join(self.destination, "srun.sh")

#    self.destination = os.path.dirname(self.submit_path)

  def eval_params(self):
    # -N <nnodes> (optional)
    self.sbatch_contents = self.substitute(self.sbatch_contents, "<nnodes>",
      str(self.params.shifter.nnodes))

    # -n <nproc> (optional)
    self.sbatch_contents = self.substitute(self.sbatch_contents, "<nproc>",
      str(self.params.shifter.nproc))

    # -W <walltime> (optional)
    if self.params.wall_time is not None:
      hours = self.params.wall_time // 60
      minutes = self.params.wall_time % 60
      wt_str = "%02d:%02d:00" % (hours, minutes)
      self.sbatch_contents = self.substitute(self.sbatch_contents, "<walltime>",
        wt_str)

    # -p <partition>
    self.sbatch_contents = self.substitute(self.sbatch_contents, "<partition>",
      self.params.shifter.partition)
    # --job-name
    self.sbatch_contents = self.substitute(self.sbatch_contents, "<jobname>",
      self.params.shifter.jobname)

    # <srun_script>
    self.sbatch_contents = self.substitute(self.sbatch_contents, "<srun_script>",
      self.srun_path)

    # <command> and any extra args
    if len(self.params.extra_args) > 0:
      self.srun_contents = self.substitute(self.srun_contents, "<command>",
        "<command> %s" % " ".join(self.params.extra_args))
    self.srun_contents = self.substitute(self.srun_contents, "<command>",
      self.command)
  def generate_submit_command(self):
    return self.params.shifter.submit_command + self.sbatch_path

  def encapsulate_submit(self):
    pass

  def generate_sbatch_script(self):
    sb = open(self.sbatch_path, "wb")
    sb.write(self.sbatch_contents)
    sb.write("\n")
    sb.close()
    self.make_executable(self.sbatch_path)

  def generate_srun_script(self):
    sr = open(self.srun_path, "wb")
    sr.write(self.srun_contents)
    sr.write("\n")
    sr.close()
    self.make_executable(self.srun_path)

  def write_script(self):
    self.generate_sbatch_script()
    self.generate_srun_script()

class get_custom_submit_command(get_submit_command):

  def customize_for_method(self):
    # template for the script to be submitted, beginning with #!
    self.script_template = self.params.custom.submit_script_template
    if not os.path.exists(self.template):
      raise Sorry("Custom submission template file not found: %s" % self.template)

    # template for the submission command itself, e.g. qsub -n <nproc> -q <queue> script.sh
    self.command_template = self.params.custom.submit_command_template
    if self.command_template is None:
      raise Sorry("Custom submit command must be specified for custom environments.")

  def eval_params(self):
    # any changes to the script to be submitted
    script = open(self.script_template, "rb")
    self.script_contents = script.read()
    script.close()

    # <command> and any <args>
    if len(self.params.extra_args) > 0:
      self.script_contents = self.script_contents.replace("<command>",
        "<command> %s" % " ".join(self.params.extra_args))
    self.script_contents = self.script_contents.replace("<command>", self.command)

    # other changes to the contents of the script
    for marker, value in [
      ("<queue>", self.params.queue),
      ("<nproc>", self.params.nproc),
      ("<memory>", self.params.memory),
      ("<walltime>", self.params.custom.wall_time_string),
      ("<outfile>", os.path.join(self.stdoutdir, self.log_name)),
      ("<errfile>", os.path.join(self.stdoutdir, self.err_name)),
      ("<envscripts>", self.params.env_script)]:
      self.script_contents = self.substitute(self.script_contents, marker, value)

    # any changes to the submission command
    # <script> and any extra <options>
    if len(self.params.extra_options) > 0:
      self.submit_command_contents = self.params.custom.submit_command_template.replace("<script>",
        "%s <script>" % " ".join(self.params.extra_options))
    self.submit_command_contents = self.submit_command_contents.replace("<script>",
      self.submit_path)

    # other changes to the submission command
    for marker, value in [
      ("<queue>", self.params.queue),
      ("<nproc>", self.params.nproc),
      ("<memory>", self.params.memory),
      ("<walltime>", self.params.custom.wall_time_string),
      ("<outfile>", os.path.join(self.stdoutdir, self.log_name)),
      ("<errfile>", os.path.join(self.stdoutdir, self.err_name))]:
      self.submit_command_contents = self.substitute(self.submit_command_contents, marker, value)

  def write_script(self):
    f = open(self.submit_path, "wb")
    f.write(self.script_contents)
    f.write("\n")
    f.close()

  def generate_submit_command(self):
    return self.submit_command_contents

def get_submit_command_chooser(command, submit_path, stdoutdir, params,
                               log_name="log.out", err_name="log.err", job_name=None):
  if params.method == "lsf":
    choice = get_lsf_submit_command
  elif params.method == "sge":
    choice = get_sge_submit_command
  elif params.method == "pbs":
    choice = get_pbs_submit_command
  elif params.method == "shifter":
    choice = get_shifter_submit_command
  elif params.method == "custom":
    choice = get_custom_submit_command
  else:
    raise Sorry("Multiprocessing method %s not recognized" % params.method)
  command_generator = choice(command, submit_path, stdoutdir, params,
                             log_name=log_name, err_name=err_name, job_name=job_name)
  return command_generator()
