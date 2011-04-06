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

import libtbx.load_env

class Queue(object):
  """
  Queue object to receive data from jobs running on remote machines

  Data transfer is achieved via files. It is safe to use any number of
  Queue objects in the same directory, even with a matching identifier
  """

  def __init__(self, identifier):

    self.root = "%s_%d_%d" % ( identifier, os.getpid(), id( self ) )
    self.count = itertools.count()


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


  def get(self):

    while True:
      fname = self.next()

      if fname is not None:
        break

      time.sleep( 0.1 )

    data = pickle.load( open( fname, "rb" ) )
    os.remove( fname )
    return data


  def next(self):

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

    if not waiting:
      return None

    else:
      return waiting[0][1]


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
    self.process = subprocess.Popen(
      cmd,
      stdin = subprocess.PIPE,
      stdout = subprocess.PIPE,
      stderr = subprocess.STDOUT
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

  def __init__(self, name, target, args = (), kwargs = {}):

    self.name = "%s_%d_%d" % ( name, os.getpid(), id( self ) )
    self.target = target
    self.args = args
    self.kwargs = kwargs
    self.jobid = None


  def start(self):

    if self.jobid is not None:
      raise RuntimeError, "start called second time"

    self.write_input_data()

    cmd = pbs_interface(
      name = self.name,
      out = self.out_file(),
      err = self.err_file()
      )
    process = subprocess.Popen(
      cmd,
      stdin = subprocess.PIPE,
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
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


def sge_interface(name, out, err):
  """
  Interface to Sun Grid Engine (SGE)
  """

  return ( "qsub", "-S", "/bin/sh", "-cwd", "-N", name, "-sync", "y",
      "-o", out, "-e", err )


def lsf_interface(name, out, err):
  """
  Interface to Load Sharing Facility (LSF)
  """

  return ( "bsub", "-K", "-J", name,
      "-o", out, "-e", err )


def pbs_interface(name, out, err):
  """
  Interface to Portable Batch System (PBS)
  """

  return ( "qsub", "-d", ".", "-N", name,
      "-o", out, "-e", err )

