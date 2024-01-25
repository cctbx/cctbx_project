from __future__ import absolute_import, division, print_function

import subprocess
import time

from libtbx.queuing_system_utils.processing import util
from libtbx.queuing_system_utils.processing import errors

# These are not very efficient
SGE_REGEX = util.get_lazy_initialized_regex( pattern = r"Following jobs do not exist" )


def sge_single_evaluate(out, err):
  "Evaluates SGE text output in single mode"

  if err:
    if SGE_REGEX().search( err ):
      return True

    else:
      raise errors.AbnormalExitError("SGE error:\n%s" % err)

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
      raise errors.AbnormalExitError("PBS error:\n%s" % err)

  state = PBS_EVAL_REGEX().search( out )

  if not state:
    raise errors.ExtractionError("Unexpected response from PBS: %r" % out)

  return state.group(1) == "C"


SLURM_SEARCH_REGEX = util.get_lazy_initialized_regex(
  pattern = r"slurm_load_jobs error: Invalid job id specified",
  )
SLURM_CODES = set( [ "PD", "R", "CA", "CF", "CG", "CD", "F", "TO", "NF" ] )

def slurm_single_evaluate(out, err):
  "Evaluates Slurm text output in single mode"

  if err:
    if SLURM_SEARCH_REGEX().search( err ):
      return True

    else:
      raise errors.AbnormalExitError("Slurm error:\n%s" % err)

  state = out.strip()

  if state not in SLURM_CODES:
    raise errors.ExtractionError("Unexpected response from Slurm: %r" % out)

  return state == "CD"


class SinglePoller(object):
  """
  Polls status for single jobs
  """

  def __init__(self, command, evaluate):

    self.command = command
    self.evaluate = evaluate


  def is_finished(self, jobid):

    try:
      process = subprocess.Popen(
        self.command + [ jobid ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    except OSError as e:
      cmdline = " ".join( self.command + [ jobid ] )
      raise errors.ExecutableError("'%s': %s" % ( cmdline, e ))

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


  @classmethod
  def Slurm(cls):

    return cls(
      command = [ "squeue", "-o", "%t", "-h", "-j" ],
      evaluate = slurm_single_evaluate,
      )


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

    try:
      process = subprocess.Popen(
        self.command,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    except OSError as e:
      raise errors.ExecutableError("'%s': %s" % ( " ".join( self.command ), e ))

    ( out, err ) = process.communicate()

    if process.poll():
      message = "Poll error: '%s' exited abnormally (code %s, message %s)" % (
        " ".join( self.command ),
        process.poll(),
        err,
        )
      raise errors.AbnormalExitError(message)

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
      raise ValueError("Unknown job id")


  def new_job_submitted(self, jobid):

    self.running.add( jobid )


  @classmethod
  def SGE(cls, waittime = 5):

    return cls(
      command = [ "qstat", "-xml" ],
      evaluate = sge_xml_evaluate,
      waittime = waittime,
      )


  @classmethod
  def LSF(cls, waittime = 5):

    return cls(
      command = [ "bjobs" ],
      evaluate = lsf_text_evaluate,
      waittime = waittime,
      )


  @classmethod
  def PBS(cls, waittime = 5):

    return cls(
      command = [ "qstat", "-x" ],
      evaluate = pbs_xml_evaluate,
      waittime = waittime,
      )


  @classmethod
  def PBSPro(cls, waittime = 5):

    return cls(
      command = [ "qstat" ],
      evaluate = pbspro_text_evaluate,
      waittime = waittime,
      )


  @classmethod
  def Condor(cls, waittime = 5):

    return cls(
      command = [ "condor_q", "-xml" ],
      evaluate = condor_xml_evaluate,
      waittime = waittime,
      )


  @classmethod
  def Slurm(cls, waittime = 0.5):

    return cls(
      command = [ "squeue", "-h", "-o", "%i %t" ],
      evaluate = slurm_text_evaluate,
      waittime = waittime,
      )


def sge_xml_evaluate(out, running, completed):
  "Parses SGE xml output"

  import xml.etree.ElementTree as ET
  root = ET.fromstring( out )

  for n in root.iter( "job_list" ):
    status_node = n.find( "JB_job_number" )
    assert status_node is not None
    running.add( status_node.text )


CONDOR_XML_OUTPUT_REGEX = util.get_lazy_initialized_regex(
  pattern = r"^.*?(?=<\?xml)",
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
    raise errors.ExtractionError("Unexpected response from LSF: %r" % out)

  regex = LSF_CENTRAL_JOBID_REGEX()

  for line in out.splitlines()[1:]:
    match = regex.match( line )

    if not match:
      raise errors.ExtractionError("Unexpected response from LSF: %r" % out)

    jobid = match.group( 1 )
    running.add( jobid )


PBSPRO_CENTRAL_HEADER_REGEX = util.get_lazy_initialized_regex(
  pattern = r"Job id\s+Name\s+User\s+Time Use\s+S\s+Queue\s*\n[ -]*"
  )
PBSPRO_CENTRAL_JOBID_REGEX = LSF_CENTRAL_JOBID_REGEX

def pbspro_text_evaluate(out, running, completed):
  "Parses PBSPro qstat text output"

  if not out.strip():
    return

  if not PBSPRO_CENTRAL_HEADER_REGEX().match( out ):
    raise errors.ExtractionError("Unexpected response from PBSPro: %r" % out)

  regex = PBSPRO_CENTRAL_JOBID_REGEX()

  for line in out.splitlines()[2:]:
    match = regex.match( line )

    if not match:
      raise errors.ExtractionError("Unexpected response from PBSPro: %r" % out)

    jobid = match.group( 1 )
    running.add( jobid )


def slurm_text_evaluate(out, running, completed):
  "Parses Slurm squeue text output"

  for line in out.splitlines():
    pieces = line.split()

    if len( pieces ) != 2:
      raise errors.ExtractionError("Unexpected response from Slurm: %r" % line)

    ( jobid, status ) = pieces

    assert status in SLURM_CODES

    if status == "CD":
      completed.add( jobid )

    else:
      running.add( jobid )
