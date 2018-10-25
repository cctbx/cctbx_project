from __future__ import absolute_import, division, print_function

import os
import subprocess
import time

from libtbx.queuing_system_utils.processing import errors

# Timeouts
class TimedTimeout(object):
  """
  Timeout after given time
  """

  def __init__(self, max_delay):

    self.max_delay = max_delay


  def __call__(self, waittime):

    if 0 < self.max_delay:
      waittime = min( self.max_delay, waittime )
      self.max_delay -= waittime
      time.sleep( waittime )
      return False

    return True


def NoTimeout(waittime):
  """
  No timeout
  """

  time.sleep( waittime )
  return False


# Status
class NotSubmitted(object):
  """
  Signal that the job has not been submitted yet
  Note: this is a singleton
  """

  @staticmethod
  def is_finished():

    raise RuntimeError("job has not been submitted yet")


  @staticmethod
  def is_submitted():

    return False


  @staticmethod
  def terminate():

    raise RuntimeError("job has not been submitted yet")


  @staticmethod
  def exit_code():

    raise RuntimeError("job has not been submitted yet")


  @staticmethod
  def outcome(timeout, polltime):

    raise RuntimeError("job has not been submitted yet")


  @staticmethod
  def cleanup():

    raise RuntimeError("job has not been submitted yet")


class Finished(object):
  """
  Signals that the job is complete
  """

  def __init__(self, stdout, stderr, exitcode):

    self.stdout = stdout
    self.stderr = stderr
    self.exitcode = exitcode


  @staticmethod
  def is_finished():

    return True


  @staticmethod
  def is_submitted():

    return True


  @staticmethod
  def terminate():

    pass


  def exit_code(self):

    return self.exitcode


  def outcome(self, timeout, polltime):

    return self


  @staticmethod
  def cleanup():

    pass


  # Extra method only available for this class
  def strip(self):

    self.stdout = None
    self.stderr = None


class JobStatus(object):
  """
  Status class
  """

  def __init__(self, outfile, errfile, additional):

    self.outfile = outfile
    self.errfile = errfile
    self.files = [ self.outfile, self.errfile ]
    self.files.extend( additional )


  @staticmethod
  def is_submitted():

    return True


  def cleanup(self):

    for fname in self.files:
      if os.path.exists( fname ):
        os.remove( fname )

    self.outfile = None
    self.errfile = None


  # Protected methods
  def get_stdout(self, timeout, polltime):

    return self.get_file_content(
      filename = self.outfile,
      timeout = timeout,
      polltime = polltime,
      )


  def get_stderr(self, timeout, polltime):

    return self.get_file_content(
      filename = self.errfile,
      timeout = timeout,
      polltime = polltime,
      )


  def get_file_content(self, filename, timeout, polltime):

    if filename is not None:
      timeoutobj = self.get_timeout( timeout = timeout )

      while True:
        try:
          with open( filename ) as infile:
            return infile.read().lstrip()

        except IOError:
          pass

        if timeoutobj( waittime = polltime ):
          return None

    else:
      return None


  # Internal
  def get_timeout(self, timeout):

    if timeout is not None:
      return TimedTimeout( max_delay = timeout )

    else:
      return NoTimeout


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


  def exit_code(self):

    self.process.poll()


  def outcome(self, timeout, polltime):

    exitcode = self.exit_code()
    assert exitcode is not None

    stdout = self.get_stdout( timeout = timeout, polltime = polltime )
    stderr = self.get_stderr( timeout = timeout, polltime = polltime )

    if stdout is None or stderr is None:
      stdout = self.process.stdout.read()
      stderr = self.process.stderr.read()
      raise errors.AbnormalExitError("Queue error:%s\n%s" % ( stdout, stderr ))

    return Finished( stdout = stdout, stderr = stderr, exitcode = exitcode )


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


  @staticmethod
  def is_submitted():

    return True


  def terminate(self):

    self.process.terminate()


  def exit_code(self):

    self.process.poll()


  def outcome(self, timeout, polltime):

    exitcode = self.exit_code()
    assert exitcode is not None

    stdout = self.process.stdout.read()
    stderr = self.process.stderr.read()

    return Finished( stdout = stdout, stderr = stderr, exitcode = exitcode )


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

    except OSError as e:
      raise errors.ExecutableError("'%s %s': %s" % ( self.qdel, self.jobid, e ))

    process.communicate()


  def exit_code(self):

    return None


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

  def outcome(self, timeout, polltime):

    # Error file is used to return exit code
    stderr = self.get_stderr( timeout = timeout, polltime = polltime )

    if stderr is not None:
      try:
        exit_code = extract_exit_code_text( output = stderr )

      except errors.ExtractionError:
        exit_code = 1 # Assume that job died/cancelled

    else:
      exit_code = 0

    stdout = self.get_stdout( timeout = timeout, polltime = polltime )

    if exit_code != 0:
      return Finished( stdout = "", stderr = stdout, exitcode = exit_code )

    else:
      return Finished( stdout = stdout, stderr = "", exitcode = exit_code )


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

  def outcome(self, timeout, polltime):

    # Use accounting file to return exit code
    try:
      process = subprocess.Popen(
        [ "qacct", "-j", self.jobid ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    except OSError as e:
      raise errors.ExecutableError("'qacct -j %s': %s" % ( self.jobid, e ))

    ( out, err ) = process.communicate()

    if process.poll():
      raise errors.AbnormalExitError("Accounting error:\n%s" % err)

    exit_code = extract_exit_code_text( output = out )

    return Finished(
      stdout = self.get_stdout( timeout = timeout, polltime = polltime ),
      stderr = self.get_stderr( timeout = timeout, polltime = polltime ),
      exitcode = exit_code,
      )


# Helper method
def extract_exit_code_text(output):

  for line in output.splitlines():
    tokens = line.split()

    if len( tokens ) == 2 and tokens[0] == "exit_status":
      exit_code = int( tokens[1] ) # this should work
      return exit_code

  else:
    raise errors.ExtractionError("Unexpected output: %r" % output)


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

  def __init__(self, jobid, poller, outfile, errfile, logfile, additional, qdel):

    self.logfile = logfile
    super( LogfileStrategy, self ).__init__(
      jobid = jobid,
      poller = poller,
      outfile = outfile,
      errfile = errfile,
      additional = additional + [ self.logfile ],
      qdel = qdel,
      )


  def get_stdlog(self, timeout, polltime):

    return self.get_file_content(
      filename = self.logfile,
      timeout = timeout,
      polltime = polltime,
      )


  def outcome(self, timeout, polltime):

    # Use logfile to return exit code
    stdlog = self.get_stdlog( timeout = timeout, polltime = polltime )

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

    return Finished(
      stdout = self.get_stdout( timeout = timeout, polltime = polltime ),
      stderr = self.get_stderr( timeout = timeout, polltime = polltime ),
      exitcode = exit_code,
      )


  def cleanup(self):

    super( LogfileStrategy, self ).cleanup()
    self.logfile = None

