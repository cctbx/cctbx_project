"Portable Batch System (PBS) utilities"
from __future__ import absolute_import, division, print_function

import sys
import os

def eval_env(variable_name):
  return eval(os.environ.get(variable_name, "None"))

def raise_if_none(value, variable_name):
  if (value is None):
    raise RuntimeError(
      "Required environment variable not defined: %s" % variable_name)

class chunk_info(object):

  def __init__(self):
    self.pbs_array_size = eval_env("MY_PBS_ARRAY_SIZE")
    self.pbs_arrayid_offset = eval_env("MY_PBS_ARRAYID_OFFSET")
    self.pbs_arrayid = eval_env("PBS_ARRAYID")
    self.mpi_num_procs = eval_env("OMPI_MCA_ns_nds_num_procs")
    self.mpi_vpid = eval_env("OMPI_MCA_ns_nds_vpid")

  def show(self, out=None, prefix="", even_if_none=False):
    if (out is None): out = sys.stdout
    if (self.pbs_array_size is not None or even_if_none):
      print(prefix+"MY_PBS_ARRAY_SIZE =", self.pbs_array_size, file=out)
    if (self.pbs_arrayid_offset is not None or even_if_none):
      print(prefix+"MY_PBS_ARRAYID_OFFSET =", self.pbs_arrayid_offset, file=out)
    if (self.pbs_arrayid is not None or even_if_none):
      print(prefix+"PBS_ARRAYID =", self.pbs_arrayid, file=out)
    if (self.mpi_num_procs is not None or even_if_none):
      print(prefix+"OMPI_MCA_ns_nds_num_procs =", self.mpi_num_procs, file=out)
    if (self.mpi_vpid is not None or even_if_none):
      print(prefix+"OMPI_MCA_ns_nds_vpid =", self.mpi_vpid, file=out)
    return self

  def have_array(self):
    return self.pbs_array_size is not None

  def as_n_i_pair(self):
    if (not self.have_array()): return 1, 0
    raise_if_none(self.pbs_arrayid_offset, "MY_PBS_ARRAYID_OFFSET")
    raise_if_none(self.pbs_arrayid, "PBS_ARRAYID")
    if (self.mpi_num_procs is None and self.mpi_vpid is None):
      n = self.pbs_array_size
      i = self.pbs_arrayid-1 + self.pbs_arrayid_offset
    else:
      raise_if_none(self.mpi_num_procs, "OMPI_MCA_ns_nds_num_procs")
      raise_if_none(self.mpi_vpid, "OMPI_MCA_ns_nds_vpid")
      n = self.pbs_array_size * self.mpi_num_procs
      i = (self.pbs_arrayid-1 + self.pbs_arrayid_offset) * self.mpi_num_procs + self.mpi_vpid
    return n, i

# XXX this is probably only going to work with newer versions, e.g. Torque
def qstat_parse():
  from libtbx.queuing_system_utils.sge_utils import qstat_items
  from libtbx import easy_run
  from xml.dom import minidom
  qstat_out = easy_run.fully_buffered(
    command="qstat -x").raise_if_errors().stdout_lines
  result = []
  if (len(qstat_out) == 0):
    return result
  xml = minidom.parseString("\n".join(qstat_out))
  jobs = xml.getElementsByTagName("Job")
  def get_tag_content(node, tag_name):
    node = node.getElementsByTagName(tag_name)[0].childNodes[0]
    assert (node.nodeType == node.TEXT_NODE)
    return node.data
  for job in jobs :
    id = get_tag_content(job, "Job_Id")
    name = get_tag_content(job, "Job_Name")
    user = get_tag_content(job, "Job_Owner").split("@")[0]
    state = get_tag_content(job, "job_state")
    prior = get_tag_content(job, "Priority")
    submit = get_tag_content(job, "start_time")
    queue = get_tag_content(job, "queue")
    nodect = get_tag_content(job, "nodect")
    result.append(qstat_items(
      job_id=id,
      prior=prior,
      name=name,
      user=user,
      state=state,
      submit=submit,
      queue=queue,
      slots=nodect,
      ja_task_id="",
      qtype="pbs")) # XXX is there any equivalent for this?
  return result

if (__name__ == "__main__"):
  n,i = chunk_info().show(prefix="*** ", even_if_none=True).as_n_i_pair()
  print("n,i:", n,i)
  print("OK")
