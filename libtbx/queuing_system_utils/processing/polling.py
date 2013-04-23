from __future__ import division

import subprocess
import time

from libtbx.queuing_system_utils.processing import util

# These are not very efficient
SGE_REGEX = util.get_lazy_initialized_regex( pattern = r"Following jobs do not exist" )


def sge_single_evaluate(out, err):
  "Evaluates SGE text output in single mode"

  if err:
    if SGE_REGEX().search( err ):
      return True

    else:
      raise RuntimeError, "SGE error:\n%s" % err

  else:
    return False


def lsf_single_evaluate(out, err):
  "Evaluates LSF text output in single mode"

  if err:
    return True

  else:
    return False


PBS_SEARCH_REGEX = util.get_lazy_initialized_regex( pattern = r"Unknown Job Id" )
PBS_EVAL_REGEX = util.get_lazy_initialized_regex( pattern = r"job_state\s*=\s*(\w+)" )


def pbs_single_evaluate(out, err):
  "Evaluates PBS text output in single mode"

  if err:
    if PBS_SEARCH_REGEX().search( err ):
      return True

    else:
      raise RuntimeError, "PBS error:\n%s" % err

  state = PBS_EVAL_REGEX().search( out )

  if not state:
    raise RuntimeError, "Unexpected response from queue:\n%s" % out

  return state.group(1) == "C"


class SinglePoller(object):
  """
  Polls status for single jobs
  """

  def __init__(self, command, evaluate):

    self.command = command
    self.evaluate = evaluate


  def is_finished(self, jobid):

    process = subprocess.Popen(
      self.command + [ jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    return self.evaluate( out = out, err = err )


  def new_job_submitted(self, jobid):

    pass


  @classmethod
  def SGE(cls):

    return cls( command = [ "qstat", "-j" ], evaluate = sge_single_evaluate )


  @classmethod
  def LSF(cls):

    return cls( command = [ "bjobs" ], evaluate = lsf_single_evaluate )


  @classmethod
  def PBS(cls):

    return cls( command = [ "qstat", "-f" ], evaluate = pbs_single_evaluate )


class CentralPoller(object):
  """
  Polls job status in batch
  """

  def __init__(self, command, evaluate, waittime = 5):

    self.waittime = waittime
    self.command = command
    self.evaluate = evaluate
    self.running = set()
    self.completed = set()
    self.update()


  def update(self):

    process = subprocess.Popen(
      self.command,
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      raise RuntimeError, "Polling error:\n%s" % err

    self.running.clear()
    self.completed.clear()

    self.evaluate( out = out, running = self.running, completed = self.completed )
    self.polltime = time.time()


  def is_finished(self, jobid):

    now = time.time()

    if self.waittime < ( now - self.polltime ) or now < self.polltime:
      self.update()

    if jobid in self.running:
      return False

    elif jobid in self.completed:
      return True

    else:
      raise ValueError, "Unknown job id"


  def new_job_submitted(self, jobid):

    self.running.add( jobid )


  @classmethod
  def SGE(cls):

    return cls( command = [ "qstat", "-xml" ], evaluate = sge_xml_evaluate )


  @classmethod
  def LSF(cls):

    return cls( command = [ "bjobs" ], evaluate = lsf_text_evaluate )


  @classmethod
  def PBS(cls):

    return cls( command = [ "qstat", "-x" ], evaluate = pbs_xml_evaluate )


  @classmethod
  def Condor(cls):

    return cls( command = [ "condor_q", "-xml" ], evaluate = condor_xml_evaluate )


def sge_xml_evaluate(out, running, completed):
  "Parses SGE xml output"

  import xml.etree.ElementTree as ET
  root = ET.fromstring( out )

  for n in root.iter( "job_list" ):
    status_node = n.find( "JB_job_number" )
    assert status_node is not None
    running.add( status_node.text )


CONDOR_XML_OUTPUT_REGEX = util.get_lazy_initialized_regex(
  pattern = "^.*?(?=<\?xml)",
  flags = [ "DOTALL" ],
  )

def condor_xml_evaluate(out, running, completed):
  "Parses Condor xml output"

  # Necessary to fix broken XML from Condor 7.2
  xml = CONDOR_XML_OUTPUT_REGEX().sub( "", out )

  import xml.etree.ElementTree as ET
  root = ET.fromstring( xml )

  for job_node in root.iter( "c" ):
    res = [ e for e in job_node.iter( "a" ) if e.attrib.get( "n" ) == "ClusterId" ]
    assert res
    number_node = res[-1].find( "i" )
    assert number_node is not None
    running.add( number_node.text )


def pbs_xml_evaluate(out, running, completed):
  "Parses PBS xml output"

  if not out:
    return

  import xml.etree.ElementTree as ET
  root = ET.fromstring( out )
  for n in root.iter( "Job" ):
    status_node = n.find( "job_state" )
    assert status_node is not None

    if status_node.text == "C":
      completed.add( n.text )

    else:
      running.add( n.text )

LSF_CENTRAL_NO_RESULTS_REGEX = util.get_lazy_initialized_regex(
  pattern = r"No unfinished job found"
  )
LSF_CENTRAL_HEADER_REGEX = util.get_lazy_initialized_regex(
  pattern = r"JOBID\s+USER\s+STAT\s+QUEUE\s+FROM_HOST\s+EXEC_HOST\s+JOB_NAME\s+SUBMIT_TIME"
  )
LSF_CENTRAL_JOBID_REGEX = util.get_lazy_initialized_regex(
  pattern = r"\s*(\S+)\s.+"
  )

def lsf_text_evaluate(out, running, completed):
  "Parses LSF bjobs text output"

  if LSF_CENTRAL_NO_RESULTS_REGEX().search( out ):
    return

  if not LSF_CENTRAL_HEADER_REGEX().match( out ):
    raise RuntimeError, "Unexpected response from queue"

  regex = LSF_CENTRAL_JOBID_REGEX()

  for line in out.splitlines()[1:]:
    match = regex.match( line )

    if not match:
      raise RuntimeError, "Unexpected response from queue"

    jobid = match.group( 1 )
    running.add( jobid )

