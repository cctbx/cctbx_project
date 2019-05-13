from __future__ import division
from __future__ import print_function

from libtbx import easy_run

class JobStopper(object):
  def __init__(self, queueing_system):
    self.queueing_system = queueing_system
    if self.queueing_system in ["mpi", "lsf"]:
      self.command = "bkill %s"
    elif self.queueing_system == 'pbs':
      self.command = "qdel %s"
    else:
      raise NotImplementedError("job stopper not implemented for %s queueing system" \
      % self.queueing_system)

  def stop_job(self, submission_id):
    result = easy_run.fully_buffered(command=self.command%submission_id)
    status = "\n".join(result.stdout_lines)
    error = "\n".join(result.stderr_lines)
    print(status)
    print(error)

class QueueInterrogator(object):
  """A queue monitor that returns the status of a given queued job, or ERR if the job cannot
  be found in the queue."""
  def __init__(self, queueing_system):
    self.queueing_system = queueing_system
    if self.queueing_system in ["mpi", "lsf"]:
      self.command = "bjobs %s | grep %s | awk '{ print $3 }'"
    elif self.queueing_system == 'pbs':
      self.command = "qstat -H %s | tail -n 1 | awk '{ print $10 }'"
    else:
      raise NotImplementedError(
      "queue interrogator not implemented for %s queueing system"%self.queueing_system)

  def query(self, submission_id):
    if self.queueing_system in ["mpi", "lsf"]:
      result = easy_run.fully_buffered(command=self.command % \
        (submission_id, submission_id))
    elif self.queueing_system == 'pbs':
      result = easy_run.fully_buffered(command=self.command%submission_id)
    status = "\n".join(result.stdout_lines)
    error = "\n".join(result.stderr_lines)
    if error != "" and not "Warning: job being submitted without an AFS token." in error:
      if "not found" in error:
        return "ERR"
      else:
        return error
    else:
      return status

class LogReader(object):
  """A log reader that distinguishes between expected results of successful and unsuccessful
  log file termination, and returns an error message if the log file cannot be found or read."""
  def __init__(self, queueing_system):
    self.queueing_system = queueing_system
    if self.queueing_system in ["mpi", "lsf", "pbs"]:
      self.command = "tail -17 %s | head -1"
    else:
      raise NotImplementedError(
      "queue interrogator not implemented for %s queueing system"%self.queueing_system)

  def read_result(self, log_path):
    result = easy_run.fully_buffered(command=self.command % log_path)
    status = "\n".join(result.stdout_lines)
    error = "\n".join(result.stderr_lines)
    if error != "" and not "Warning: job being submitted without an AFS token." in error:
      return "Error reading log file."
    else:
      return status

class SubmissionTracker(object):
  """An object that uses the QueueInterrogator and LogReader to query a queueing system and log
  file to determine the status of a queued job."""
  def __init__(self, params):
    self.queueing_system = params.mp.method
    self.interrogator = QueueInterrogator(self.queueing_system)
    self.reader = LogReader(self.queueing_system)
  def track(self, submission_id, log_path):
    raise NotImplementedError("Override me!")

class LSFSubmissionTracker(SubmissionTracker):
  def track(self, submission_id, log_path):
    from xfel.ui.db.job import known_job_statuses
    if submission_id is None:
      return "UNKWN"
    status = self.interrogator.query(submission_id)
    if status == "ERR":
      log_status = self.reader.read_result(log_path)
      if log_status == "Successfully completed.":
        return "DONE"
      else:
        return "ERR" # error querying the queueing system
    elif status in known_job_statuses:
      return status
    else:
      print("Found an unknown status", status)

class PBSSubmissionTracker(SubmissionTracker):
  def track(self, submission_id, log_path):
    if submission_id is None:
      return "UNKWN"
    status = self.interrogator.query(submission_id)
    if status == "F":
      return "DONE"
    elif status in ["Q", "H", "S"]:
      return "PEND"
    elif status in ["R", "T", "W", "E"]:
      return "RUN"
    else:
      print("Found an unknown status", status)

class TrackerFactory(object):
  @staticmethod
  def from_params(params):
    if params.mp.method in ['mpi', 'lsf']:
      return LSFSubmissionTracker(params)
    elif params.mp.method == 'pbs':
      return PBSSubmissionTracker(params)
