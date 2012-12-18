"""
Generic module to provide parallel job execution on queuing systems

Provides drop-in replacement classes to those defined in the multiprocessing
module (Queue and Process), with certain restrictions placed by the pickle
module
"""
from __future__ import division
from __future__ import with_statement

import cPickle as pickle
import os
import time
import itertools
import glob
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

    with open( tmp_name, "wb" ) as ofile:
      pickle.dump( obj, ofile )

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


class Job(object):
  """
  Job object to execute function calls on remote machines accessible via
  queuing systems

  Restrictions: target, args and kwargs has to be pickleable
  """

  def __init__(self, qinterface, target, args = (), kwargs = {}):

    self.qinterface = qinterface

    self.target = target
    self.args = args
    self.kwargs = kwargs

    self.status = None


  @property
  def name(self):

    return "%s_%d" % ( self.qinterface.root, id( self ) )

  @property
  def jobid (self) :
    return getattr(self.status, "jobid", None)

  def start(self):

    if self.status is not None:
      raise RuntimeError, "start called second time"

    data = self.qinterface.input(
      name = self.name,
      target = self.target,
      args = self.args,
      kwargs = self.kwargs,
      )

    self.status = self.qinterface.submitter(
      name = self.name,
      executable = self.qinterface.executable,
      script = self.qinterface.script % data.script(),
      include = self.qinterface.include,
      extra = data.files(),
      )


  def is_alive(self):

    if self.status is None:
      raise RuntimeError, "job has not been submitted yet"

    return not self.status.is_finished()


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )

    ( stdout, stderr, exit_code ) = self.status.results()

    if stdout:
      print stdout

    if stderr:
      print stderr

    self.status.cleanup()

    if exit_code != 0:
      raise RuntimeError, "%s: exit code = %s" % ( self, exit_code )


  def __str__(self):

    return "%s(name = '%s')" % ( self.__class__.__name__, self.name )


class QueueHandler(object):
  """
  Handles submission requests for common queuing systems
  """
  SCRIPT = \
"""\
import cPickle as pickle
%s
target( *args, **kwargs )
"""

  def __init__(
    self,
    submitter,
    input,
    include,
    root,
    executable = "libtbx.python",
    script = SCRIPT,
    ):

    self.submitter = submitter
    self.root = "%s%s" % ( root, os.getpid() )
    self.input = input
    self.include = include

    self.executable = executable
    self.script = script


  def Job(self, target, args = (), kwargs = {}):

    return Job(
      qinterface = self,
      target = target,
      args = args,
      kwargs = kwargs,
      )


def get_libtbx_env_setpaths():

  import libtbx.load_env
  return libtbx.env.under_build( "setpaths.sh" )


def SGE(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  ):

  if not command:
    ( command, switches ) = ( "qsub", [] )

  else:
    ( command, switches ) = lex_command_line( command = command )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if handler is None:
      from libtbx.queuing_system_utils.processing import status
      handler = status.StdStreamStrategy

    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.SGECentralPoller()

    submitter = submission.AsynchronousCmdLine.SGE(
      poller = poller,
      handler = handler,
      command = command,
      extra = switches,
      )

  else:
    submitter = submission.Synchronous.SGE( command = command, extra = switches )

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    )


def LSF(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  ):

  if not command:
    ( command, switches ) = ( "bsub", [] )

  else:
    ( command, switches ) = lex_command_line( command = command )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if handler is None:
      from libtbx.queuing_system_utils.processing import status
      handler = status.StdStreamStrategy

    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.LSFPoller()

    submitter = submission.AsynchronousCmdLine.LSF(
      poller = poller,
      handler = handler,
      command = command,
      extra = switches,
      )

  else:
    submitter = submission.Synchronous.LSF( command = command, extra = switches )

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    )


def PBS(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  ):

  if not command:
    ( command, switches) = ( "qsub", [] )

  else:
    ( command, switches ) = lex_command_line( command = command )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if handler is None:
      from libtbx.queuing_system_utils.processing import status
      handler = status.StdStreamStrategy

    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.PBSPoller()

    submitter = submission.AsynchronousCmdLine.PBS(
      poller = poller,
      handler = handler,
      command = command,
      extra = switches,
      )

  else:
    raise RuntimeError, "PBS does not support synchronous submission"

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    )


def Condor(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  ):

  if not command:
    ( command, switches ) = ( "condor_submit", [] )

  else:
    ( command, switches ) = lex_command_line( command = command )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CondorCentralPoller()

    submitter = submission.AsynchronousScript.Condor(
      poller = poller,
      command = command,
      extra = switches,
      )

  else:
    raise RuntimeError, "Condor does not support synchronous submission"

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    )


def lex_command_line(command):
  assert isinstance(command, str)
  import shlex
  lexed = shlex.split( command )
  assert lexed
  return ( lexed[0], lexed[1:] )


INTERFACE_FOR = {
  "sge": SGE,
  "lsf": LSF,
  "pbs": PBS,
  "condor": Condor,
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

  qinterface = INTERFACE_FOR[ platform ](
      name = name,
      command = command,
      asynchronous = asynchronous,
      )
  return qinterface.Job( target = target )

