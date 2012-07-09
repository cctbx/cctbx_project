"""
Generic module to provide parallel job execution on queuing systems

Provides drop-in replacement classes to those defined in the multiprocessing
module (Queue and Process), with certain restrictions placed by the pickle
module
"""

import cPickle as pickle
import subprocess
import os
import time
import itertools
import glob
import re
from Queue import Empty as QueueEmptyException

class InstantTimeout(object):
  """
  Timeout immediately
  """

  def delay(self, waittime):

    raise QueueEmptyException, "No data found in queue"


class TimedTimeout(object):
  """
  Timeout after given time
  """

  def __init__(self, max_delay):

      self.max_delay = max_delay


  def delay(self, waittime):

    if 0 < self.max_delay:
      waittime = min( self.max_delay, waittime )
      self.max_delay -= waittime
      time.sleep( waittime )

    else:
      raise QueueEmptyException, "No data found in queue within timeout"


class NoTimeout(object):
  """
  No timeout
  """

  def delay(self, waittime):

    time.sleep( waittime )


class Queue(object):
  """
  Queue object to receive data from jobs running on remote machines

  Data transfer is achieved via files. It is safe to use any number of
  Queue objects in the same directory, even with a matching identifier
  """

  def __init__(self, identifier, waittime = 1):

    self.waittime = waittime
    self.root = "%s_%d_%d" % ( identifier, os.getpid(), id( self ) )
    self.count = itertools.count()


  def put(self, obj):

    index = self.count.next()
    # Writing a tempfile and renaming it may prevent reading incomplete files
    tmp_name = "tmp_%s.%d" % ( self.root, index )
    assert not os.path.exists( tmp_name )
    ofile = open( tmp_name, "wb" )
    pickle.dump( obj, ofile )
    ofile.close()
    target_name = "%s.%d" % ( self.root, index )
    assert not os.path.exists( target_name )
    os.rename( tmp_name, target_name )


  def get(self, block = True, timeout = None):

    if not block:
      predicate = InstantTimeout()

    else:
      if timeout is not None:
        predicate = TimedTimeout( max_delay = timeout )

      else:
        predicate = NoTimeout()

    while True:
      try:
        data = self.next()

      except StopIteration:
        predicate.delay( waittime = self.waittime )
        continue

      break

    return data


  def next(self):

    waiting = []

    for fname in glob.glob( "%s.*" % self.root ):
      ( root, ext ) = os.path.splitext( fname )

      try:
        number = int( ext[1:] )

      except ValueError:
        continue

      waiting.append( ( number, fname ) )

    if not waiting:
      raise StopIteration

    selected = min( waiting, key = lambda p: p[0] )[1]
    data = pickle.load( open( selected, "rb" ) )
    os.remove( selected )
    return data


class TemporaryFileInput(object):
  """
  Sends data by writing out a temporary file
  """

  SCRIPT= "( target, args, kwargs ) = pickle.load( open( \"%s\" ) )"

  def __init__(self, name, target, args, kwargs):

    self.target = "%s.target" % name
    ifile = open( self.target, "wb" )
    pickle.dump( ( target, args, kwargs ), ifile )
    ifile.close()


  def script(self):

    return self.SCRIPT % self.target


  def files(self):

    return [ self.target ]


class StdinInput(object):
  """
  Sends data by sending along the command file as a pickled string
  """

  SCRIPT = "( target, args, kwargs ) = pickle.loads( \"%s\" ) )"

  def __init__(self, name, target, args, kwargs):

    self.target = target
    self.args = args
    self.kwargs = kwargs


  def script(self):

    data = pickle.dumps( ( target, args, kwargs ) )
    return self.SCRIPT % data


  def files(self):

    return []


class Job(object):
  """
  Job object to execute function calls on remote machines accessible via
  queuing systems

  Restrictions: target has to be pickleable
  """

  SCRIPT = \
"""\
source %s
libtbx.python << EOF
import pickle
%s
target( *args, **kwargs )
EOF
"""

  def __init__(self, qinterface, target, args = (), kwargs = {}):

    self.qinterface = qinterface
    self.name = "%s_%d" % ( self.qinterface.root, id( self ) )
    self.data = self.qinterface.input(
        name = self.name,
        target = target,
        args = args,
        kwargs = kwargs,
        )
    self.status = None


  def start(self):

    if self.status is not None:
      raise RuntimeError, "start called second time"

    self.status = self.qinterface.submit(
      name = self.name,
      out = self.out_file(),
      err = self.err_file(),
      input = self.SCRIPT % ( self.qinterface.include, self.data.script() ),
      )


  def is_alive(self):

    if self.status is None:
      raise RuntimeError, "job has not been submitted yet"

    return not self.status.is_finished()


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )

    if os.path.exists( self.err_file() ):
      error = open( self.err_file() ).read()

      if error:
        raise RuntimeError, error

    for fname in [ self.out_file(), self.err_file() ] + self.data.files():
      if os.path.exists( fname ):
        os.remove( fname )


  def out_file(self):

    return "%s.out" % self.name


  def err_file(self):

    return "%s.err" % self.name


class SynchronousJobStatus(object):
  """
  Determines job status for synchronous jobs
  """

  def __init__(self, process):

    self.process = process


  def is_finished(self):

    return self.process.poll() is not None


class AsynchronousJobStatus(object):
  """
  Determines job status for asynchronous jobs
  """

  def __init__(self, jobid, poller):

    self.jobid = jobid
    self.poller = poller


  def is_finished(self):

    return self.poller.is_finished( jobid = jobid )


# These are not very efficient
class SGEPoller(object):
  """
  Polls job status for SGE
  """

  MISSING = re.compile( r"Following jobs do not exist" )

  def is_finished(self, jobid):

    process = subprocess.Popen(
      [ "qstat", "-j", jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      if self.MISSING.search( err ):
        raise ValueError, "Unknown job id"

      else:
        raise RuntimeError, "SGE error:\n%s" % err

    else:
      return False


class LSFPoller(object):
  """
  Polls job status for LSF
  """

  def is_finished(self, jobid):

    process = subprocess.Popen(
      [ "bjobs", jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise ValueError, "Unknown job id"

    else:
      return False


class PBSPoller(object):
  """
  Polls job status for PBS
  """

  STATE = re.compile( r"job_state\s*=\s*(\w+)" )
  MISSING = re.compile( r"Unknown Job Id" )

  def is_finished(self, jobid):

    process = subprocess.Popen(
      [ "qstat", "-f", jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      if self.MISSING.search( err ):
        raise ValueError, "Unknown job id"

      else:
        raise RuntimeError, "PBS error:\n%s" % err

    state = self.STATE.search( out )

    if not state:
      raise RuntimeError, "Unexpected response from queue:\n%s" % out

    return state.group(1) == "C"


class Submission(object):
  """
  Handles job submissions
  """

  def __init__(self, switches):

    self.switches = switches


  def get_process_object(self, command_list):

    try:
      process = subprocess.Popen(
        command_list + self.switches,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        )

    except OSError, e:
      raise RuntimeError, "Error while executing: '%s': %s" % (
        " ".join( command_list ),
        e,
        )

    return process


class SynchronousSubmission(Submission):
  """
  Submits jobs synchronously
  """

  def __call__(self, command_list, input):

    process = self.get_process_object( command_list = command_list )
    process.stdin.write( input )
    process.stdin.close()
    return SynchronousJobStatus( process = process )


  @classmethod
  def SGE(cls):

    return cls( switches = [ "-sync", "y" ] )


  @classmethod
  def LSF(cls):

    return cls( switches = [ "-K" ] )


class AsynchronousSubmission(Submission):
  """
  Submits jobs asynchronously
  """

  def __init__(self, switches, extract, poller):

    super( AsynchronousSubmission, self ).__init__( switches = switches )
    self.extract = extract
    self.poller = poller


  def __call__(self, command_list, input):

    process = self.get_process_object( command_list = command_list )
    ( out, err ) = process.communicate( input = input )

    if err:
      raise RuntimeError, err

    return AsynchronousJobStatus(
      jobid = self.extract( output = out ),
      poller = self.poller,
      )


  @classmethod
  def SGE(cls, poller):

    return cls(
        switches = [ "-terse" ],
        extract = cls.generic_jobid_extract,
        poller = poller,
        )


  @classmethod
  def LSF(cls, poller):

    return cls(
      switches = [],
      extract = cls.lsf_jobid_extract,
      poller = poller,
      )


  @classmethod
  def PBS(cls, poller):

    return cls(
        switches = [],
        extract = cls.generic_jobid_extract,
        poller = poller,
        )


  @staticmethod
  def lsf_jobid_extract(output):

    import re
    regex = re.compile( r"Job <(\d+)> is submitted" )

    match = regex.search( output )

    if not match:
      raise RuntimeError, "Unexpected response from queuing system"

    return match.group(1)


  @staticmethod
  def generic_jobid_extract(output):

    return output.strip()


class QueueHandler(object):
  """
  Handles submission requests for common queuing systems
  """

  def __init__(
    self,
    command,
    name_switch,
    out_switch,
    err_switch,
    extra_switches,
    submission,
    input,
    include,
    root = None,
    ):

    self.command = command
    self.name_switch = name_switch
    self.out_switch = out_switch
    self.err_switch = err_switch
    self.extra_switches = extra_switches
    self.submission = submission
    self.root = "%s%s" % ( ( root + "_" ) if root else "", os.getpid() )
    self.input = input
    self.include = include


  def command_line(self, name, out, err):

    return ( [ self.command ] + self.extra_switches
      + [ self.name_switch, name, self.out_switch, out, self.err_switch, err ] )


  def submit(self, name, out, err, input):

    return self.submission(
      command_list = self.command_line( name = name, out = out, err = err ),
      input = input,
      )


  def Job(self, target, args = (), kwargs = {}):

    return Job(
      qinterface = self,
      target = target,
      args = args,
      kwargs = kwargs,
      )


  @classmethod
  def SGE(cls, command, switches, submission, input, include, root):

    return cls(
      command = command,
      name_switch = "-N",
      out_switch = "-o",
      err_switch = "-e",
      extra_switches = switches,
      submission = submission,
      input = input,
      include = include,
      root = root,
      )


  @classmethod
  def LSF(cls, command, switches, submission, input, include, root):

    return cls(
      command = command,
      name_switch = "-J",
      out_switch = "-o",
      err_switch = "-e",
      extra_switches = switches,
      submission = submission,
      input = input,
      include = include,
      root = root,
      )


  @classmethod
  def PBS(cls, command, switches, submission, input, include, root):

    return cls(
      command = command,
      name_switch = "-N",
      out_switch = "-o",
      err_switch = "-e",
      extra_switches = switches,
      submission = submission,
      input = input,
      include = include,
      root = root,
      )


def SGE(
  name = "libtbx_python",
  command = None,
  asynchronous = False,
  switches = [ "-S", "/bin/sh", "-cwd" ],
  input = TemporaryFileInput,
  include = None,
  ):

  if not command:
    ( command, extra ) = ( "qsub", [] )

  else:
    ( command, extra ) = lex_command_line( command = command )

  if asynchronous:
    submission = AsynchronousSubmission.SGE( poller = SGEPoller() )

  else:
    submission = SynchronousSubmission.SGE()

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler.SGE(
    command = command,
    switches = extra + switches,
    submission =  submission,
    input = input,
    include = include,
    root = name,
    )


def LSF(
  name = "libtbx_python",
  command = None,
  asynchronous = False,
  switches = [],
  input = TemporaryFileInput,
  include = None,
  ):

  if not command:
    ( command, extra ) = ( "bsub", [] )

  else:
    ( command, extra ) = lex_command_line( command = command )

  if asynchronous:
    submission = AsynchronousSubmission.LSF( poller = LSFPoller() )

  else:
    submission = SynchronousSubmission.LSF()

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler.LSF(
    command = command,
    switches = extra + switches,
    submission =  submission,
    input = input,
    include = include,
    root = name,
    )


def PBS(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  switches = [ "-d", "." ],
  input = TemporaryFileInput,
  include = None,
  ):

  if not command:
    ( command, extra ) = ( "qsub", [] )

  else:
    ( command, extra ) = lex_command_line( command = command )

  if asynchronous:
    submission = AsynchronousSubmission.PBS( poller = PBSPoller() )

  else:
    raise RuntimeError, "PBS does not support synchronous submission"

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler.PBS(
    command = command,
    switches = extra + switches,
    submission =  submission,
    input = input,
    include = include,
    root = name,
    )


def lex_command_line(command):

  assert command
  import shlex
  lexed = shlex.split( command )
  assert lexed
  return ( lexed[0], lexed[1:] )


def get_libtbx_env_setpaths():

  import libtbx.load_env
  return libtbx.env.under_build( "setpaths.sh" )


INTERFACE_FOR = {
  "sge": SGE,
  "lsf": LSF,
  "pbs": PBS,
  }

def qsub (
  target,
  name="libtbx_python",
  platform="sge",
  command=None,
  asynchronous=True
  ):

  assert hasattr(target, "__call__")

  if platform not in INTERFACE_FOR:
    raise RuntimeError, "Unknown platform: %s" % platform

  qinterace = INTERFACE_FOR[ platform ](
      name = name,
      command = command,
      asynchronous = asynchronous,
      )
  return qinterace.Job( target = target )

