from __future__ import division

import subprocess
import time

from libtbx.object_oriented_patterns import lazy_initialization

# These are not very efficient
def get_regex(pattern):

  import re
  return re.compile( pattern )


SGE_REGEX = lazy_initialization(
  calculation = lambda: get_regex( pattern = r"Following jobs do not exist" ),
  )


def sge_single_evaluate(out, err):

  if err:
    if SGE_REGEX().search( err ):
      return True

    else:
      raise RuntimeError, "SGE error:\n%s" % err

  else:
    return False


def lsf_single_evaluate(out, err):

  if err:
    return True

  else:
    return False


PBS_SEARCH_REGEX = lazy_initialization(
  calculation = lambda: get_regex( pattern = r"Unknown Job Id" )
  )
PBS_EVAL_REGEX = lazy_initialization(
  calculation = lambda: get_regex( pattern = r"job_state\s*=\s*(\w+)" )
  )


def pbs_single_evaluate(out, err):

  if err:
    if PBS_SEARCH_REGEX.search( err ):
      return True

    else:
      raise RuntimeError, "PBS error:\n%s" % err

  state = PBS_EVAL_REGEX.search( out )

  if not state:
    raise RuntimeError, "Unexpected response from queue:\n%s" % out

  return state.group(1) == "C"


class SinglePoller(object):
  """
  Polls status for single jobs
  """

  def __init__(self, cmdline, evaluator):

    self.cmdline = cmdline
    self.evaluator = evaluator


  def is_finished(self, jobid):

    process = subprocess.Popen(
      self.cmdline + ( jobid, ),
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    return self.evaluator( out = out, err = err )


  def new_job_submitted(self, jobid):

    pass


class CentralPoller(object):
  """
  Polls job status in batch
  """

  def __init__(self, waittime = 5):

    self.waittime = waittime
    self.running = set()
    self.update()
    self.polltime = time.time()


  def is_finished(self, jobid):

    now = time.time()

    if self.waittime < ( now - self.polltime ) or now < self.polltime:
      self.update()
      self.polltime = now

    if jobid in self.running:
      return False

    else:
      raise ValueError, "Unknown job id"


  def new_job_submitted(self, jobid):

    self.running.add( jobid )


class SGECentralPoller(CentralPoller):
  """
  Polls job status for SGE in batch
  """

  def update(self):

    process = subprocess.Popen(
      [ "qstat", "-xml" ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise RuntimeError, "SGE error:\n%s" % err

    import xml.etree.ElementTree as ET
    root = ET.fromstring( out )
    self.running = set()

    for n in root.iter( "job_list" ):
      status_node = n.find( "JB_job_number" )
      assert status_node is not None
      self.running.add( status_node.text )


class CondorCentralPoller(CentralPoller):
  """
  Polls job status for Condor in batch
  """

  def __init__(self, waittime = 5):

    import re
    self.regex = re.compile( "^.*?(?=<\?xml)", re.DOTALL )
    super( CondorCentralPoller, self ).__init__( waittime = waittime )


  def update(self):

    process = subprocess.Popen(
      [ "condor_q", "-xml" ], # "-attribute", "ClusterId"
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise RuntimeError, "Condor error:\n%s" % err

    # Necessary to fix broken XML from Condor 7.2
    xml = self.regex.sub( "", out )

    import xml.etree.ElementTree as ET
    root = ET.fromstring( xml )
    self.running = set()

    for job_node in root.iter( "c" ):
      res = [ e for e in job_node.iter( "a" ) if e.attrib.get( "n" ) == "ClusterId" ]
      assert res
      number_node = res[-1].find( "i" )
      assert number_node is not None
      self.running.add( number_node.text )


class PBSCentralPoller(object):
  """
  Polls job status for PBS in batch
  """

  def __init__(self, waittime = 5):

    self.waittime = waittime
    self.running = set()
    self.completed = set()
    self.update()
    self.polltime = time.time()


  def is_finished(self, jobid):

    now = time.time()

    if self.waittime < ( now - self.polltime ) or now < self.polltime:
      self.update()
      self.polltime = now

    if jobid in self.running:
      return False

    elif jobid in self.completed:
      return True

    else:
      raise ValueError, "Unknown job id"


  def new_job_submitted(self, jobid):

    self.running.add( jobid )


  def update(self):

    process = subprocess.Popen(
      [ "qstat", "-x" ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise RuntimeError, "PBS error:\n%s" % err

    self.running = set()
    self.completed = set()

    if not out:
      return

    import xml.etree.ElementTree as ET
    root = ET.fromstring( out )
    for n in root.iter( "Job" ):
      status_node = n.find( "job_state" )
      assert status_node is not None

      if status_node.text == "C":
        self.completed.add( n.text )

      else:
        self.running.add( n.text )

