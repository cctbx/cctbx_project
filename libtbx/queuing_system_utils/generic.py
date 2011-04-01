import cPickle as pickle
import subprocess
import os
import time
import itertools
import glob

import libtbx.load_env

class Queue(object):
  """
  Queuing system queue, single use only
  """

  def __init__(self, identifier, work_dir = "."):

    assert os.path.isdir( work_dir )
    self.root = os.path.join(
        work_dir,
        "%s_%d_%d" % ( identifier, os.getpid(), id( self ) )
        )
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
  Generic-purpose queueing system Job
  """

  SCRIPT = \
"""\
source %s
cd %s
libtbx.python << EOF
import pickle
( target, args, kwargs ) = pickle.load( open( "%s.target" ) )
target( *args, **kwargs )
EOF
"""
  SETPATHS = libtbx.env.under_build("setpaths.sh")

  def __init__(self, name, target, qinterface, args = (), kwargs = {}, work_dir = "."):

    assert os.path.isdir( work_dir )
    self.name = "%s_%d_%d" % ( name, os.getpid(), id( self ) )
    self.target = target
    self.args = args
    self.kwargs = kwargs
    self.qinterface = qinterface
    self.work_dir = work_dir
    self.process = None


  def start(self):

    if self.process is not None:
      raise RuntimeError, "start called second time"

    ifile = open( os.path.join( self.work_dir, "%s.target" % self.name ), "wb" )
    pickle.dump( ( self.target, self.args, self.kwargs ), ifile )
    ifile.close()

    cmd = self.qinterface( name = self.name, work_dir = self.work_dir )
    self.process = subprocess.Popen(
        cmd,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.STDOUT

        )
    self.process.stdin.write(
        self.SCRIPT % ( self.SETPATHS, self.work_dir, self.name )
        )
    self.process.stdin.close()


  def is_alive(self):

    if self.process is None:
      raise RuntimeError, "job has not been submitted yet"

    return self.process.poll() is None


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )

    root = os.path.join( self.work_dir, self.name )
    cleanups = [ "%s.target" % root, "%s.out" % root, "%s.err" % root ]

    if os.path.exists( cleanups[-1] ):
      error = open( cleanups[-1] ).read()

      if error:
        raise RuntimeError, error

    for fname in cleanups:
      if os.path.exists( fname ):
        os.remove( fname )


def sge_interface(name, work_dir):

  root = os.path.join( work_dir, name )
  return ( "qsub", "-S", "/bin/sh", "-cwd", "-N", name, "-sync", "y",
      "-o", "%s.out" % root, "-e", "%s.err" % root )


def lsf_interface(name, work_dir):

  root = os.path.join( work_dir, name )
  return ( "bsub", "-K", "-J", name,
      "-o", "%s.out" % root, "-e", "%s.err" % root )

