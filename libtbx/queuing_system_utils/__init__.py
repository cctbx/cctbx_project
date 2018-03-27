from __future__ import absolute_import, division, print_function

import sys

from libtbx import Auto

class chunk_manager(object):

  def __init__(self, n, i):
    assert n > 0
    assert i >= 0
    assert i < n
    self.n = n
    self.i = i
    self.queuing_system_info = None

  def easy_all(self, log_format=Auto, out=Auto):
    self.queuing_system_overrides_chunk()
    self.redirect_chunk_stdout_and_stderr(log_format=log_format, out=out)
    return self

  def skip_iteration(self, i):
    return (i % self.n != self.i)

  def queuing_system_overrides_chunk(self):
    from libtbx.queuing_system_utils import pbs_utils, sge_utils
    pbs_info = pbs_utils.chunk_info()
    sge_info = sge_utils.info()
    assert [pbs_info, sge_info].count(None) <= 1
    if pbs_info.have_array():
      self.queuing_system_info = pbs_info
      n, i = pbs_info.as_n_i_pair()
      self.n = max(self.n, n)
      self.i = i
    elif sge_info.have_array():
      self.queuing_system_info = sge_info
      self.n = max(self.n, sge_info.last)
      self.i = sge_info.id - 1
    return self

  def redirect_chunk_stdout_and_stderr(self,
        log_format=Auto,
        out=Auto,
        have_array=False):
    if self.n == 1: return
    log_name = None
    if not have_array:
      i = self.queuing_system_info
      if i is not None and i.have_array():
        have_array = True
    if have_array:
      if log_format is Auto: log_format="log%%0%dd"
      fmt = log_format % max(3, len("%d" % (self.n-1)))
      log_name = fmt % self.i
      log = open(log_name, "w")
      sys.stdout = log
      sys.stderr = log
    from libtbx.utils import host_and_user
    if out is Auto: out = sys.stdout
    if out is not None:
      host_and_user().show(out=out)
      print("chunk.n:", self.n, file=out)
      print("chunk.i:", self.i, file=out)
      if log_name:
        print("log_name:", log_name, file=out)
      print(file=out)
    return self

# XXX tested on SGE only so far (2012-12-19)
def qdel(job_id, platform):
  """
  Stop a queue job.  Supports the same platforms as 'processing' sub-module,
  but primarily used by the Phenix GUI.
  """
  from libtbx import easy_run
  assert platform in ("sge", "lsf", "pbs", "condor", "pbspro", "slurm")
  cmd = None
  if platform in ("sge", "pbs", "pbspro", "slurm"):
    cmd = "qdel %s" % job_id
  elif platform == "lsf":
    cmd = "bkill %s" % job_id
  elif platform == "condor":
    cmd = "condor_rm %s" % job_id
  assert cmd
  qdel_out = easy_run.fully_buffered(
    command=cmd).raise_if_errors().stdout_lines
  print("\n".join(qdel_out))
  # XXX this is specific to SGE - need error handling for other systems too
  for line in qdel_out:
    if "denied" in line:
      if "does not exist" in line:  # SGE job does not exist
        pass
      else:
        raise RuntimeError("\n".join(qdel_out))
  return True
