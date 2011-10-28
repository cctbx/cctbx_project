"""
Generic module to provide parallel job execution on queuing systems

Provides drop-in replacement classes to those defined in the multiprocessing
module (Queue and Job), with certain restrictions placed by the pickle module
"""

import cPickle as pickle
import subprocess
import os
import time
import itertools
import glob
import re
from Queue import Empty as QueueEmptyException

import libtbx.load_env


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

    if waittime <= self.max_delay:
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

  WAITTIME = 0.1

  def __init__(self, identifier):

    self.root = "%s_%d_%d" % ( identifier, os.getpid(), id( self ) )
    self.count = itertools.count()
    self.waiting = []


  def put(self, obj):

    index = self.count.next()
    # Writing a tempfile and renaming it may prevent reading incomplete files
    tmp_name = "%s.%d.tmp" % ( self.root, index )
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
      fname = self.next()

      if fname is not None:
        break

      predicate.delay( waittime = self.WAITTIME )

    assert fname is not None
    data = pickle.load( open( fname, "rb" ) )
    os.remove( fname )
    return data


  def next(self):

    if not self.waiting:
      self.read_waiting()

    try:
      return self.waiting.pop( 0 )

    except IndexError:
      return None


  def read_waiting(self):

    waiting = []

    for fname in glob.glob( "%s.*" % self.root ):
      if fname[-4] == ".tmp":
        continue

      ( root, ext ) = os.path.splitext( fname )

      try:
        number = int( ext[1:] )

      except ValueError:
        continue

      waiting.append( ( number, fname ) )

    waiting.sort( key = lambda p: p[0] )
    self.waiting = [ p[1] for p in waiting ]


class Job(object):
  """
  Job object to execute function calls on remote machines accessible via
  queuing systems

  Data transfer is achieved via files. It is safe to use any number of
  Job objects in the same directory, even with a matching identifier

  Restrictions: target has to be pickleable
  """

  SCRIPT = \
"""\
source %s
libtbx.python << EOF
import pickle
( target, args, kwargs ) = pickle.load( open( "%s.target" ) )
target( *args, **kwargs )
EOF
"""
  SETPATHS = libtbx.env.under_build("setpaths.sh")

  def __init__(self, name, target, qinterface, args = (), kwargs = {}):

    self.name = "%s_%d_%d" % ( name, os.getpid(), id( self ) )
    self.target = target
    self.args = args
    self.kwargs = kwargs
    self.qinterface = qinterface
    self.process = None


  def start(self):

    if self.process is not None:
      raise RuntimeError, "start called second time"

    self.write_input_data()

    cmd = self.qinterface(
      name = self.name,
      out = self.out_file(),
      err = self.err_file()
      )

    try:
        self.process = subprocess.Popen(
          cmd,
          stdin = subprocess.PIPE,
          stdout = subprocess.PIPE,
          stderr = subprocess.STDOUT
          )

    except OSError, e:
        raise RuntimeError, "Error while executing: '%s': %s" % (
            " ".join( cmd ),
            e,
            )

    self.process.stdin.write( self.SCRIPT % ( self.SETPATHS, self.name ) )
    self.process.stdin.close()


  def is_alive(self):

    if self.process is None:
      raise RuntimeError, "job has not been submitted yet"

    return self.process.poll() is None


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )

    if os.path.exists( self.err_file() ):
      error = open( self.err_file() ).read()

      if error:
        raise RuntimeError, error

    for fname in [ self.target_file(), self.out_file(), self.err_file() ]:
      if os.path.exists( fname ):
        os.remove( fname )


  def write_input_data(self):

    ifile = open( self.target_file(), "wb" )
    pickle.dump( ( self.target, self.args, self.kwargs ), ifile )
    ifile.close()


  def target_file(self):

    return "%s.target" % self.name


  def out_file(self):

    return "%s.out" % self.name


  def err_file(self):

    return "%s.err" % self.name


class PBSJob(Job):
  """
  Job object to execute function calls on remote machines accessible via
  Portable Batch System (PBS)

  Data transfer is achieved via files. It is safe to use any number of
  Job objects in the same directory, even with a matching identifier

  If PBS support for 'wait for job to finish' is available, this class can be
  merged into Job

  Restrictions: target has to be picklable
  """

  REGEX = re.compile( r"job_state\s*=\s*(\w+)" )

  def __init__(self, name, target, qinterface, args = (), kwargs = {}):

    self.name = "%s_%d_%d" % ( name, os.getpid(), id( self ) )
    self.target = target
    self.args = args
    self.kwargs = kwargs
    self.qinterface = qinterface
    self.jobid = None


  def start(self):

    if self.jobid is not None:
      raise RuntimeError, "start called second time"

    self.write_input_data()

    cmd = self.qinterface(
      name = self.name,
      out = self.out_file(),
      err = self.err_file()
      )

    try:
        process = subprocess.Popen(
          cmd,
          stdin = subprocess.PIPE,
          stdout = subprocess.PIPE,
          stderr = subprocess.PIPE
          )

    except OSError, e:
        raise RuntimeError, "Error while executing: '%s': %s" % (
            " ".join( cmd ),
            e,
            )

    ( out, err ) = process.communicate(
      input = self.SCRIPT % ( self.SETPATHS, self.name )
      )

    if err:
      raise RuntimeError, err

    assert out is not None
    self.jobid = out.strip()


  def is_alive(self):

    return self.job_status() != "C"


  def job_status(self):

    if self.jobid is None:
      raise RuntimeError, "job has not been submitted yet"

    process = subprocess.Popen(
      ( "qstat", "-f", self.jobid ),
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      # may be better to indicate finish, in case job has already been deleted
      raise RuntimeError, "Jobid: %s\nPBS error: %s" % ( self.jobid, err )

    m = self.REGEX.search( out )

    if not m:
      raise RuntimeError, "Incorrect qstat output: %s" % out

    return m.group( 1 )


class sge_interface(object):
  """
  Interface to Sun Grid Engine (SGE)
  """

  COMMAND = "qsub"

  def __init__(self, command):

    if command is None:
      self.command = ( self.COMMAND, )

    else:
      self.command = command


  def __call__(self, name, out, err):

    return self.command + ( "-S", "/bin/sh", "-cwd", "-N", name, "-sync", "y",
      "-o", out, "-e", err )


class lsf_interface(object):
  """
  Interface to Load Sharing Facility (LSF)
  """

  COMMAND = "bsub"

  def __init__(self, command):

    if command is None:
      self.command = ( self.COMMAND, )

    else:
      self.command = command


  def __call__(self, name, out, err):

    return self.command + ( "-K", "-J", name,
      "-o", out, "-e", err )


class pbs_interface(object):
  """
  Interface to Portable Batch System (PBS)
  """

  COMMAND = "qsub"

  def __init__(self, command):

    if command is None:
      self.command = ( self.COMMAND, )

    else:
      self.command = command


  def __call__(self, name, out, err):

    return self.command + ( "-d", ".", "-N", name,
      "-o", out, "-e", err )

