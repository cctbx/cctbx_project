from __future__ import absolute_import, division, print_function

from libtbx import easy_run

# Mapping from raw Slurm (sacct) job states to cctbx job statuses.
SLURM_STATUS_MAP = {'COMPLETED': 'DONE',
                    'COMPLETING': 'RUN',
                    'FAILED': 'EXIT',
                    'PENDING': 'PEND',
                    'PREEMPTED': 'SUSP',
                    'RUNNING': 'RUN',
                    'SUSPENDED': 'SUSP',
                    'STOPPED': 'SUSP',
                    'CANCELLED': 'EXIT',
                    'TIMEOUT': 'TIMEOUT',
                    'OUT_OF_ME': 'EXIT',      # truncated 'OUT_OF_MEMORY'
                    'OUT_OF_MEMORY': 'EXIT',
                   }

def _map_slurm_state(raw):
  """Normalize a raw sacct state string and map it to a cctbx job status.

  Handles both truncated forms (e.g. 'CANCELLED+', 'OUT_OF_ME') and full forms
  with trailing detail (e.g. 'CANCELLED by 12345')."""
  status = raw.split()[0].rstrip('+') if raw else ''
  if status not in SLURM_STATUS_MAP:
    print('Unknown job status', status)
  return SLURM_STATUS_MAP.get(status, 'UNKWN')

# Per-chunk statuses that are terminal and stable, so they never need to be
# re-queried. 'UNKWN' is deliberately excluded (it is transient -- a job not yet
# registered with the scheduler reads as UNKWN until it appears).
CACHEABLE_TERMINAL = frozenset(["DONE", "EXIT", "DELETED", "ERR", "SUBMIT_FAIL", "TIMEOUT"])

def _aggregate_ensemble(statuses):
  """Collapse the statuses of a job's chunks (comma-separated submission ids,
  e.g. ensemble refinement) into one overall status. Active states win so the
  job keeps being tracked; once nothing is active, a single failed chunk marks
  the whole job EXIT, and only an all-clear set is DONE."""
  s = [st for st in statuses if st is not None]
  if not s:
    return "UNKWN"
  if "RUN" in s:
    return "RUN"
  if "PEND" in s or "SUBMITTED" in s:
    return "PEND"
  if "SUSP" in s:
    return "SUSP"
  if any(st in ("EXIT", "ERR", "TIMEOUT") for st in s):
    return "EXIT"
  if "UNKWN" in s:
    return "UNKWN"
  return "DONE"

class JobStopper(object):
  def __init__(self, queueing_system):
    self.queueing_system = queueing_system
    if self.queueing_system in ["mpi", "lsf"]:
      self.command = "bkill %s"
    elif self.queueing_system in ["pbs", "sge"]:
      self.command = "qdel %s"
    elif self.queueing_system == 'local':
      pass
    elif self.queueing_system == 'slurm' or self.queueing_system == 'shifter':
      # The current implementation of the shifter mp method assumes that we're
      # running on NERSC's systems => jobs should be tracked using the _slurm_
      # submission tracker.
      self.command = "scancel %s"
    elif self.queueing_system == 'htcondor':
      self.command = "condor_rm %s"
    else:
      raise NotImplementedError("job stopper not implemented for %s queueing system" \
      % self.queueing_system)

  def stop_job(self, submission_id):
    if not submission_id: return
    for sid in submission_id.split(','):
      if self.queueing_system == 'local':
        import psutil
        try:
          process = psutil.Process(int(sid))
        except psutil.NoSuchProcess:
          return
        for child in process.children(recursive=True): child.kill()
        process.kill()
      else:
        result = easy_run.fully_buffered(command=self.command%sid)
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
    elif self.queueing_system == 'sge':
      self.command = "qstat |grep %s | awk '{print $5}'" # previous one "qstat -j %s | awk '/job_state/ {print $3}'"
    elif self.queueing_system == 'local':
      pass
    elif self.queueing_system == 'slurm' or  self.queueing_system == 'shifter':
      # The current implementation of the shifter mp method assumes that we're
      # running on NERSC's systems => jobs should be tracked using the _slurm_
      # submission tracker.
      self.command = "sacct --job %s --format state --noheader"
    elif self.queueing_system == 'htcondor':
      self.command1 = "condor_q %s -nobatch | grep -A 1 OWNER | tail -n 1"
      self.command2 = "condor_history %s | grep -A 1 OWNER | tail -n 1"
    else:
      raise NotImplementedError(
      "queue interrogator not implemented for %s queueing system"%self.queueing_system)

  def query(self, submission_id):
    if self.queueing_system in ["mpi", "lsf"]:
      result = easy_run.fully_buffered(command=self.command % \
        (submission_id, submission_id))
    elif self.queueing_system == 'pbs':
      result = easy_run.fully_buffered(command=self.command%submission_id)
    elif self.queueing_system == 'sge':
      result = easy_run.fully_buffered(command=self.command%submission_id)
    elif self.queueing_system == 'local':
      import psutil
      try:
        process = psutil.Process(int(submission_id))
      except psutil.NoSuchProcess:
        return "DONE"
      statuses = [p.status() for p in [process] + process.children(recursive=True)]
      if 'running' in statuses: return "RUN"
      if 'sleeping' in statuses: return "RUN"
      # zombie jobs can be left because the GUI process that forked them is still running
      if len(statuses) == 1 and statuses[0] == 'zombie': return "DONE"
      return ", ".join(statuses)
    elif self.queueing_system == 'slurm' or self.queueing_system == "shifter":
      # The current implementation of the shifter mp method assumes that we're
      # running on NERSC's systems => jobs should be tracked using the _slurm_
      # submission tracker.
      result = easy_run.fully_buffered(command=self.command%submission_id)
      if len(result.stdout_lines) == 0: return 'UNKWN'
      return _map_slurm_state(result.stdout_lines[0].strip())
    elif self.queueing_system == 'htcondor':
      # (copied from the man page)
      # H = on hold, R = running, I = idle (waiting for a machine to execute on), C = completed,
      # X = removed, S = suspended (execution of a running job temporarily suspended on execute node),
      # < = transferring input (or queued to do so), and > = transferring output (or queued to do so).
      statuses = {'H': 'HOLD',
                  'R': 'RUN',
                  'I': 'PEND',
                  'C': 'DONE',
                  'X': 'DONE',
                  'S': 'SUSP',
                  '<': 'RUN',
                  '>': 'RUN'}
      for c in [self.command1, self.command2]:
        result = easy_run.fully_buffered(command=c%submission_id)
        if len(result.stdout_lines) != 1 or len(result.stdout_lines[0]) == 0: continue
        status = result.stdout_lines[0].split()[5]
        return statuses[status] if status in statuses else 'UNKWN'
      return 'ERR'
    status = "\n".join(result.stdout_lines)
    error = "\n".join(result.stderr_lines)
    if error != "" and not "Warning: job being submitted without an AFS token." in error:
      if "not found" in error:
        return "ERR"
      else:
        return error
    else:
      return status

  def query_many(self, submission_ids):
    """Query many slurm/shifter jobs with minimal load on the accounting DB.

    Live state comes from the controller (slurmctld) via a single `squeue` call,
    which is cheap and in-memory. Only jobs that have already left the queue are
    looked up in the accounting DB (slurmdbd) via a single batched `sacct`, to
    classify their terminal state once. Returns {submission_id: cctbx status};
    ids unknown to both map to 'UNKWN'. Slurm/shifter only."""
    assert self.queueing_system in ('slurm', 'shifter')
    ids = [str(sid) for sid in submission_ids if sid]
    result = {}
    if not ids:
      return result
    # 1) Controller state -- avoids slurmdbd entirely. --states=all is needed so
    #    that jobs still held in the controller in a finished state are reported.
    squeue_command = "squeue --jobs=%s --noheader --format='%%i|%%T' --states=all" \
      % ",".join(ids)
    squeue_result = easy_run.fully_buffered(command=squeue_command)
    for line in squeue_result.stdout_lines:
      fields = line.split('|')
      if len(fields) < 2:
        continue
      jobid, state = fields[0].strip(), fields[1].strip()
      if '.' in jobid:  # array/step sub-entries
        continue
      result[jobid] = _map_slurm_state(state)
    # 2) Jobs no longer known to the controller have finished and left the queue;
    #    resolve their terminal state once from the accounting DB, batched.
    missing = [sid for sid in ids if sid not in result]
    if missing:
      sacct_command = "sacct --jobs=%s --format=JobID,State --noheader -P" \
        % ",".join(missing)
      sacct_result = easy_run.fully_buffered(command=sacct_command)
      for line in sacct_result.stdout_lines:
        fields = line.split('|')
        if len(fields) < 2:
          continue
        jobid, state = fields[0].strip(), fields[1].strip()
        # Skip step rows (e.g. '12345.batch', '12345.extern', '12345.0').
        if '.' in jobid:
          continue
        if jobid in missing:
          result[jobid] = _map_slurm_state(state)
    # 3) Anything still unresolved (e.g. very recently submitted, not yet
    #    registered with either the controller or the accounting DB).
    for sid in ids:
      result.setdefault(sid, 'UNKWN')
    return result

  def get_mysql_server_hostname(self, submission_id):
    if self.queueing_system in ["mpi", "lsf"]:
      print("method to obtain hostname running MySQL server not implemented for ", self.queueing_system)
      pass
    elif self.queueing_system == 'pbs':
      print("method to obtain hostname running MySQL server not implemented for ", self.queueing_system)
      pass
    elif self.queueing_system == 'sge':
      print("method to obtain hostname running MySQL server not implemented for ", self.queueing_system)
      pass
    elif self.queueing_system == 'slurm' or self.queueing_system == "shifter":
      hostname_command = "sacct --job %s -o NODELIST --noheader | tail -n 1"
      result = easy_run.fully_buffered(command=hostname_command%submission_id).raise_if_errors()
      if result.stdout_lines:
        return result.stdout_lines[0].strip()
      else:
        return ""
    elif self.queueing_system == 'htcondor':
      print("method to obtain hostname running MySQL server not implemented for ", self.queueing_system)
      pass

class LogReader(object):
  """A log reader that distinguishes between expected results of successful and unsuccessful
  log file termination, and returns an error message if the log file cannot be found or read."""
  def __init__(self, queueing_system):
    self.queueing_system = queueing_system
    if self.queueing_system in ["mpi", "lsf", "pbs", "local", "sge"]:
      self.command = "tail -23 %s | head -1"
    elif self.queueing_system in ["slurm", "shifter", "htcondor"]:
      pass # no log reader used
    else:
      raise NotImplementedError(
      "queue interrogator not implemented for %s queueing system"%self.queueing_system)

  def read_result(self, log_path):
    if not log_path:
      return "Error reading log file, no path to log"
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
    if submission_id is None:
      return "UNKWN"
    all_statuses = [self._track(sid, log_path) for sid in submission_id.split(',')]
    # Collapse a job's (comma-separated) submission ids into one status. All ids
    # must agree, otherwise the job is considered UNKWN (e.g. an ensemble job
    # whose chunks are still in mixed states).
    if all_statuses and all(all_statuses[0] == s for s in all_statuses[1:]):
      return all_statuses[0]
    return "UNKWN"

  def track_many(self, submission_ids, log_paths):
    """Return a list of statuses, one per job. Default implementation simply loops
    over track(); SlurmSubmissionTracker overrides this to batch into one query."""
    return [self.track(sid, lp) for sid, lp in zip(submission_ids, log_paths)]

  def _track(self, submission_id, log_path):
    raise NotImplementedError("Override me!")

class LSFSubmissionTracker(SubmissionTracker):
  def _track(self, submission_id, log_path):
    from xfel.ui.db.job import known_job_statuses
    status = self.interrogator.query(submission_id)
    if status == "ERR":
      log_status = self.reader.read_result(log_path)
      if log_status == "Successfully completed.":
        return "DONE"
      elif "exit" in log_status.lower():
        return "EXIT"
      else:
        return "ERR" # error querying the queueing system
    elif status in known_job_statuses:
      return status
    else:
      print("Found an unknown status", status)

class PBSSubmissionTracker(SubmissionTracker):
  def _track(self, submission_id, log_path):
    status = self.interrogator.query(submission_id)
    if status == "F":
      return "DONE"
    elif status in ["Q", "H", "S"]:
      return "PEND"
    elif status in ["R", "T", "W", "E"]:
      return "RUN"
    else:
      print("Found an unknown status", status)

class SGESubmissionTracker(SubmissionTracker):
  def _track(self, submission_id, log_path):
    status = self.interrogator.query(submission_id)
    if status == "":
      return "DONE"
    elif status in ["s", "t", "S"]:
      return "PEND"
    elif status == "r":
      return "RUN"
    elif status == "qw":
      return "WAITING"
    else:
      print("Found an unknown status", status)

class SlurmSubmissionTracker(SubmissionTracker):
  def __init__(self, params):
    super(SlurmSubmissionTracker, self).__init__(params)
    # Cache of terminal per-chunk statuses, so finished chunks (e.g. of an
    # ensemble job) are never re-queried for the life of this tracker.
    self._chunk_cache = {}

  def _track(self, submission_id, log_path):
    return self.interrogator.query(submission_id)

  def track_many(self, submission_ids, log_paths):
    # Collect every not-yet-terminal chunk id across all jobs, query them in a
    # single squeue (+ at most one sacct) call, and aggregate per job. Chunks
    # already known to be terminal are served from the cache and not re-queried.
    to_query = set()
    for submission_id in submission_ids:
      if submission_id is None:
        continue
      for cid in submission_id.split(','):
        if cid and cid not in self._chunk_cache:
          to_query.add(cid)
    fresh = {}
    if to_query:
      fresh = self.interrogator.query_many(sorted(to_query))
      for cid, status in fresh.items():
        if status in CACHEABLE_TERMINAL:
          self._chunk_cache[cid] = status
    def chunk_status(cid):
      return self._chunk_cache.get(cid) or fresh.get(cid, "UNKWN")
    results = []
    for submission_id in submission_ids:
      if submission_id is None:
        results.append("UNKWN")
        continue
      statuses = [chunk_status(cid) for cid in submission_id.split(',') if cid]
      # Single-chunk jobs report their exact status (preserving e.g. TIMEOUT);
      # multi-chunk (ensemble) jobs are collapsed by precedence.
      if len(statuses) == 1:
        results.append(statuses[0])
      else:
        results.append(_aggregate_ensemble(statuses))
    return results

class HTCondorSubmissionTracker(SubmissionTracker):
  def _track(self, submission_id, log_path):
    return self.interrogator.query(submission_id)

class LocalSubmissionTracker(SubmissionTracker):
  def _track(self, submission_id, log_path):
    return self.interrogator.query(submission_id)

class TrackerFactory(object):
  @staticmethod
  def from_params(params):
    if params.mp.method in ['mpi', 'lsf']:
      return LSFSubmissionTracker(params)
    elif params.mp.method == 'pbs':
      return PBSSubmissionTracker(params)
    elif params.mp.method == 'sge' :
      return SGESubmissionTracker(params)
    elif params.mp.method == 'local':
      return LocalSubmissionTracker(params)
    elif params.mp.method == 'slurm':
      return SlurmSubmissionTracker(params)
    elif params.mp.method == 'shifter':
      # The current implementation of the shifter mp method assumes that we're
      # running on NERSC's systems => jobs should be tracked using the _slurm_
      # submission tracker.
      return SlurmSubmissionTracker(params)
    elif params.mp.method == 'htcondor':
      return HTCondorSubmissionTracker(params)
