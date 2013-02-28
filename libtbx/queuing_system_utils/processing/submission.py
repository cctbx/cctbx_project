from __future__ import division
from __future__ import with_statement

import subprocess

class Submission(object):
  """
  Handles job submissions
  """

  def __init__(self, command, name, out, err, extra):

    self.command = command
    self.name = name
    self.out = out
    self.err = err
    self.extra = extra


  def subprocess_command_list(self, name):

    outfile = "%s.out" % name
    errfile = "%s.err" % name
    commands = ( [ self.command ] + self.extra
      + [ self.name, name, self.out, outfile, self.err, errfile ] )

    return ( commands, outfile, errfile )


class Synchronous(Submission):
  """
  Submits jobs synchronously
  """

  SCRIPT = \
"""\
source %s
%s << EOF
%s
EOF
"""

  def __call__(self, name, executable, script, include, extra):

    ( commands, outfile, errfile ) = self.subprocess_command_list( name = name )
    process = build_submission_for( commands = commands )
    process.stdin.write( self.SCRIPT % ( include, executable, script ) )
    process.stdin.close()

    from libtbx.queuing_system_utils.processing import status

    return status.Synchronous(
      outfile = outfile,
      errfile = errfile,
      extra = extra,
      process = process,
      )


  @classmethod
  def SGE(cls, command = "qsub", extra = []):

    return cls(
      command = command,
      name = "-N",
      out = "-o",
      err = "-e",
      extra = extra + [ "-S", "/bin/sh", "-cwd", "-sync", "y" ],
      )


  @classmethod
  def LSF(cls, command = "bsub", extra = []):

    return cls(
      command = command,
      name_switch = "-J",
      out_switch = "-o",
      err_switch = "-e",
      extra = extra + [ "-K" ],
      )


class AsynchronousCmdLine(Submission):
  """
  Submits jobs asynchronously
  """

  def __init__(self, command, name, out, err, extra, extract, poller, handler):

    super( AsynchronousCmdLine, self ).__init__(
      command = command,
      name = name,
      out = out,
      err = err,
      extra = extra,
      )
    self.extract = extract
    self.poller = poller
    self.handler = handler


  def __call__(self, name, executable, script, include, extra):

    ( commands, outfile, errfile ) = self.subprocess_command_list( name = name )
    process = build_submission_for( commands = commands )

    ( out, err ) = process.communicate(
      input = self.handler.script( include = include, executable = executable, script = script )
      )

    if err:
      raise RuntimeError, err

    return self.handler(
      jobid = self.extract( output = out ),
      poller = self.poller,
      outfile = outfile,
      errfile = errfile,
      extra = extra,
      )


  @classmethod
  def SGE(cls, poller, handler, command = "qsub", extra = []):

    return cls(
      command = command,
      name = "-N",
      out = "-o",
      err = "-e",
      extra = extra + [ "-S", "/bin/sh", "-cwd", "-terse" ],
      extract = cls.generic_jobid_extract,
      poller = poller,
      handler = handler,
      )


  @classmethod
  def LSF(cls, poller, command = "bsub", extra = []):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      command = command,
      name = "-J",
      out = "-o",
      err = "-e",
      extra = extra,
      extract = cls.lsf_jobid_extract,
      poller = poller,
      handler = status.StdStreamStrategy,
      )


  @classmethod
  def PBS(cls, poller, command = "qsub", extra = []):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      command = command,
      name = "-N",
      out = "-o",
      err = "-e",
      extra = extra + [ "-d", "." ],
      extract = cls.generic_jobid_extract,
      poller = poller,
      handler = status.StdStreamStrategy,
      )


  @staticmethod
  def lsf_jobid_extract(output):

    match = get_lsf_jobid_extract_regex().search( output )

    if not match:
      raise RuntimeError, "Unexpected response from queuing system"

    return match.group(1)


  @staticmethod
  def generic_jobid_extract(output):

    return output.strip()


class AsynchronousScript(object):
  """
  Submits jobs asynchronously. Options are passed using a script
  """

  CONDOR_SCRIPT = \
"""
Executable     = %s
Universe       = vanilla
initialdir = .

Output  = %s
Error   = %s
Log = %s
Log_xml = True
Notification = Never
Queue
"""

  def __init__(self, command, script, extra, extract, poller, handler):

    self.command = command
    self.script = script
    self.extra = extra
    self.extract = extract
    self.poller = poller
    self.handler = handler


  def __call__(self, name, executable, script, include, extra):

    outfile = "%s.out" % name
    errfile = "%s.err" % name
    scriptfile = "%s.condor.script" % name
    logfile = "%s.condor.xml" % name

    with open( scriptfile, "w" ) as ofile:
      ofile.write(
        self.handler.script( include = include, executable = executable, script = script )
        )

    import os
    import stat
    os.chmod( scriptfile, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR )

    process = build_submission_for( commands = [ self.command ]  + self.extra )

    ( out, err ) = process.communicate(
      input = self.script % ( scriptfile, outfile, errfile, logfile )
      )

    if err:
      raise RuntimeError, err

    return self.handler(
      jobid = self.extract( output = out ),
      poller = self.poller,
      outfile = outfile,
      errfile = errfile,
      logfile = logfile,
      extra = extra + [ scriptfile ],
      )

  @classmethod
  def Condor(cls, poller, command = "condor_submit", script = CONDOR_SCRIPT, extra = []):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      command = command,
      script = script,
      extra = extra,
      extract = cls.condor_jobid_extract,
      poller = poller,
      handler = status.LogfileStrategy,
      )


  @staticmethod
  def condor_jobid_extract(output):

    match = get_condor_jobid_extract_regex().search( output )

    if not match:
      raise RuntimeError, "Unexpected response from queuing system"

    return match.group(1)


# Helpers
def build_submission_for(commands):

  try:
    process = subprocess.Popen(
      commands,
      stdin = subprocess.PIPE,
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE,
      )

  except OSError, e:
    raise RuntimeError, "Error: '%s': %s" % ( " ".join( commands ), e )

  return process

# Regex caching
LSF_JOBID_EXTRACT_REGEX = None

def get_lsf_jobid_extract_regex():

  global LSF_JOBID_EXTRACT_REGEX

  if LSF_JOBID_EXTRACT_REGEX is None:
    import re
    LSF_JOBID_EXTRACT_REGEX = re.compile( r"Job <(\d+)> is submitted" )

  return LSF_JOBID_EXTRACT_REGEX


CONDOR_JOBID_EXTRACT_REGEX = None

def get_condor_jobid_extract_regex():

  global CONDOR_JOBID_EXTRACT_REGEX

  if CONDOR_JOBID_EXTRACT_REGEX is None:
    import re
    CONDOR_JOBID_EXTRACT_REGEX = re.compile(
      r"\d+ job\(s\) submitted to cluster (\d+)\."
      )

  return CONDOR_JOBID_EXTRACT_REGEX

