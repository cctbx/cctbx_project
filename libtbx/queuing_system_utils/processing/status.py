from __future__ import division

import os
import subprocess


class JobStatus(object):
  """
  Status class
  """

  def __init__(self, outfile, errfile, additional):

    self.outfile = outfile
    self.errfile = errfile
    self.files = additional + [ self.outfile, self.errfile ]


  def get_stdout(self):

    if self.outfile is not None:
      stdout = open( self.outfile ).read().strip()

    else:
      stdout = None

    return stdout


  def get_stderr(self):

    if self.errfile is not None:
      stderr = open( self.errfile ).read()

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


  def results(self):

    assert self.is_finished()

    try:
      stdout = self.get_stdout()

    except IOError:
      stdout = self.process.stdout.read()

    try:
      stderr = self.get_stderr()

    except IOError:
      stderr = self.process.stderr.read()

    return ( stdout, stderr, self.process.poll() )


class Asynchronous(JobStatus):
  """
  Handler for asynchronous jobs
  """

  def __init__(self, jobid, poller, outfile, errfile, additional):

    super( Asynchronous, self ).__init__(
      outfile = outfile,
      errfile = errfile,
      additional = additional,
      )
    self.jobid = jobid
    self.poller = poller

    # Allow poller to update without actual refresh
    self.poller.new_job_submitted( jobid = self.jobid )


  def is_finished(self):

    try:
      result = self.poller.is_finished( jobid = self.jobid )

    except ValueError:
      result = True

    return result


  @classmethod
  def script(cls, include, executable, script):

    return cls.SCRIPT % ( include, executable, script )


class StdStreamStrategy(Asynchronous):
  """
  Establishes exit code through stderr
  """

  SCRIPT = \
"""\
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

      except RuntimeError:
        exit_code = 1 # Assume that job died/cancelled

    else:
      exit_code = 0

    if exit_code != 0:
      return ( "", self.get_stdout(), exit_code )

    else:
      return ( self.get_stdout(), "", exit_code )


class AccountingStrategy(Asynchronous):
  """
  Establishes exit code through accounting
  """

  SCRIPT = \
"""\
source %s
%s << EOF
%s
EOF
"""

  def results(self):

    # Use accounting file to return exit code
    process = subprocess.Popen(
      [ "qacct", "-j", self.jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise RuntimeError, "SGE error:\n%s" % err

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
    raise RuntimeError, "Unexpected output:\n%s" % output


class LogfileStrategy(Asynchronous):
  """
  Establishes exit code through queuing system logfile
  """

  SCRIPT = \
"""\
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

