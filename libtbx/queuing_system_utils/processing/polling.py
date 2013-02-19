from __future__ import division

import subprocess
import time

# These are not very efficient
class SGEPoller(object):
  """
  Polls job status for SGE
  """

  def __init__(self):

    import re
    self.regex = re.compile( r"Following jobs do not exist" )


  def is_finished(self, jobid):

    process = subprocess.Popen(
      [ "qstat", "-j", jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      if self.regex.search( err ):
        return True

      else:
        raise RuntimeError, "SGE error:\n%s" % err

    else:
      return False


  def new_job_submitted(self, jobid):

    pass


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
      return True

    else:
      return False


  def new_job_submitted(self, jobid):

    pass


class PBSPoller(object):
  """
  Polls job status for PBS
  """

  def __init__(self):

    import re
    self.state = re.compile( r"job_state\s*=\s*(\w+)" )
    self.missing = re.compile( r"Unknown Job Id" )


  def is_finished(self, jobid):

    process = subprocess.Popen(
      [ "qstat", "-f", jobid ],
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE
      )
    ( out, err ) = process.communicate()

    if err:
      if self.missing.search( err ):
        return True

      else:
        raise RuntimeError, "PBS error:\n%s" % err

    state = self.state.search( out )

    if not state:
      raise RuntimeError, "Unexpected response from queue:\n%s" % out

    return state.group(1) == "C"


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

