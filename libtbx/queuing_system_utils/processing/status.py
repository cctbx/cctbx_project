from __future__ import division

import os
import subprocess

from libtbx.queuing_system_utils.processing import errors


class JobStatus(object):
  """
  Status class
  """

  def __init__(self, outfile, errfile, additional):

    self.outfile = outfile
    self.errfile = errfile
    self.files = additional + [ self.outfile, self.errfile ]


  def get_stdout(self):

    import time

    if self.outfile is not None:
      while True:
        try:
          infile = open( self.outfile )

        except IOError:
          time.sleep( 1 )

        else:
          break

      stdout = infile.read().strip()
      infile.close()

    else:
      stdout = None

    return stdout


  def get_stderr(self):

    import time

    if self.errfile is not None:
      while True:
        try:
          infile = open( self.errfile )

        except IOError:
          time.sleep( 1 )

        else:
          break

      stderr = infile.read()
      infile.close()

    else:
      stderr = None

    return stderr


  def cleanup(self):

    for fname in self.files:
      if os.path.exists( fname ):
        os.remove( fname )

    self.outfile = None
    self.errfile = None


class Synchronous(JobStatus):
  """
  Determines job status for synchronous jobs
  """

  def __init__(self, process, outfile, errfile, additional):

    super( Synchronous, self ).__init__(
      outfile = outfile,
      errfile = errfile,
      additional = additional,
      )
    self.process = process


  def is_finished(self):

    return self.process.poll() is not None


  def terminate(self):

    self.process.terminate()


  def results(self):

    assert self.is_finished()

    try:
      stdout = self.get_stdout()
      stderr = self.get_stderr()

    except IOError:
      stdout = self.process.stdout.read()
      stderr = self.process.stderr.read()
      raise errors.AbnormalExitError, "Queue error:%s\n%s" % ( stdout, stderr )

    return ( stdout, stderr, self.process.poll() )


class SlurmStream(object):
  """
  Specialized handler for Slurm, because file handles of the shell are redirected

  Note that this can fill the buffers, resulting in the job never completing
  """

  def __init__(self, process, additional):

    self.process = process
    self.additional = additional


  def is_finished(self):

    return self.process.poll() is not None


  def terminate(self):

    self.process.terminate()


  def results(self):

    assert self.is_finished()

    stdout = self.process.stdout.read()
    stderr = self.process.stderr.read()
    return ( stdout, stderr, self.process.poll() )


  def cleanup(self):

    for fname in self.additional:
      if os.path.exists( fname ):
        os.remove( fname )


class Asynchronous(JobStatus):
  """
  Handler for asynchronous jobs
  """

  def __init__(self, jobid, poller, outfile, errfile, additional, qdel):

    super( Asynchronous, self ).__init__(
      outfile = outfile,
      errfile = errfile,
      additional = additional,
      )
    self.jobid = jobid
    self.poller = poller
    self.qdel = qdel

    # Allow poller to update without actual refresh
    self.poller.new_job_submitted( jobid = self.jobid )


  def is_finished(self):

    try:
      result = self.poller.is_finished( jobid = self.jobid )

    except ValueError:
      result = True

    return result


  def terminate(self):

    try:
      process = subprocess.Popen( self.qdel + [ self.jobid ] )

    except OSError, e:
      raise errors.ExecutableError, "'%s %s': %s" % ( self.qdel, self.jobid, e )

    process.communicate()


  @classmethod
  def script(cls, include, executable, script, cwd = "."):

    return cls.SCRIPT % ( cwd, include, executable, script )


class StdStreamStrategy(Asynchronous):
  """
  Establishes exit code through stderr
  """

  SCRIPT = \
"""\
cd %s
source %s
%s 2>&1 << EOF
%s
EOF

echo exit_status $? 1>&2
"""

  def results(self):

    # Error file is used to return exit code
    stderr = self.get_stderr()

    if stderr is not None:
      try:
        exit_code = extract_exit_code_text( output = stderr )

      except errors.ExtractionError:
        exit_code = 1 # Assume that job died/cancelled

    else:
      exit_code = 0

    if exit_code != 0:
      return ( "", self.get_stdout(), exit_code )

    else:
      return ( self.get_stdout(), "", exit_code )


class SlurmStdStreamStrategy(StdStreamStrategy):
  """
  Stdstream strategy specifically for Slurm (different script)
  """

  SCRIPT = \
"""\
#!/bin/sh
cd %s
source %s
srun %s 2>&1 << EOF
%s
EOF

echo exit_status $? 1>&2
"""


class AccountingStrategy(Asynchronous):
  """
  Establishes exit code through accounting
  """

  SCRIPT = \
"""\
cd %s
source %s
%s << EOF
%s
EOF
"""

  def results(self):

    # Use accounting file to return exit code
    try:
      process = subprocess.Popen(
        [ "qacct", "-j", self.jobid ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    except OSError, e:
      raise errors.ExecutableError, "'qacct -j %s': %s" % ( self.jobid, e )

    ( out, err ) = process.communicate()

    if process.poll():
      raise errors.AbnormalExitError, "SGE accounting error:\n%s" % err

    exit_code = extract_exit_code_text( output = out )

    return ( self.get_stdout(), self.get_stderr(), exit_code )


# Helper method
def extract_exit_code_text(output):

  for line in output.splitlines():
    tokens = line.split()

    if len( tokens ) == 2 and tokens[0] == "exit_status":
      exit_code = int( tokens[1] ) # this should work
      return exit_code

  else:
    raise errors.ExtractionError, "Unexpected output: %r" % output


class LogfileStrategy(Asynchronous):
  """
  Establishes exit code through queuing system logfile
  """

  SCRIPT = \
"""\
cd %s
#!/bin/sh
source %s
%s << EOF
%s
EOF
"""

  def __init__(self, jobid, poller, outfile, errfile, logfile, additional):

    self.logfile = logfile
    super( LogfileStrategy, self ).__init__(
      jobid = jobid,
      poller = poller,
      outfile = outfile,
      errfile = errfile,
      additional = additional + [ self.logfile ],
      )


  def get_stdlog(self):

    if self.logfile is not None:
      assert os.path.exists( self.logfile )
      stdlog = open( self.logfile ).read()

    else:
      stdlog = None

    return stdlog


  def results(self):

    # Use logfile to return exit code
    stdlog = self.get_stdlog()

    if stdlog is not None:
      import xml.etree.ElementTree as ET
      logxml = "<events>%s</events>" % stdlog.strip()
      root = ET.fromstring( logxml )
      res = [ e for e in root.iter( "a" ) if e.attrib.get( "n" ) == "ReturnValue" ]

      if not res:
        exit_code = 1 # assume an error has happened

      else:
        exit_code = int( res[-1].find( "i" ).text )

    else:
      exit_code = 0

    return ( self.get_stdout(), self.get_stderr(), exit_code )


  def cleanup(self):

    super( LogfileStrategy, self ).cleanup()
    self.logfile = None

