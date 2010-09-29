"Portable Batch System (PBS) utilities"

import sys, os

def eval_env(variable_name):
  return eval(os.environ.get(variable_name, "None"))

def raise_if_none(value, variable_name):
  if (value is None):
    raise RuntimeError(
      "Required environment variable not defined: %s" % variable_name)

class chunk_info(object):

  def __init__(O):
    O.pbs_array_size = eval_env("MY_PBS_ARRAY_SIZE")
    O.pbs_arrayid_offset = eval_env("MY_PBS_ARRAYID_OFFSET")
    O.pbs_arrayid = eval_env("PBS_ARRAYID")
    O.mpi_num_procs = eval_env("OMPI_MCA_ns_nds_num_procs")
    O.mpi_vpid = eval_env("OMPI_MCA_ns_nds_vpid")

  def show(O, out=None, prefix="", even_if_none=False):
    if (out is None): out = sys.stdout
    if (O.pbs_array_size is not None or even_if_none):
      print >> out, prefix+"MY_PBS_ARRAY_SIZE =", O.pbs_array_size
    if (O.pbs_arrayid_offset is not None or even_if_none):
      print >> out, prefix+"MY_PBS_ARRAYID_OFFSET =", O.pbs_arrayid_offset
    if (O.pbs_arrayid is not None or even_if_none):
      print >> out, prefix+"PBS_ARRAYID =", O.pbs_arrayid
    if (O.mpi_num_procs is not None or even_if_none):
      print >> out, prefix+"OMPI_MCA_ns_nds_num_procs =", O.mpi_num_procs
    if (O.mpi_vpid is not None or even_if_none):
      print >> out, prefix+"OMPI_MCA_ns_nds_vpid =", O.mpi_vpid
    return O

  def have_array(O):
    return O.pbs_array_size is not None

  def as_n_i_pair(O):
    if (not O.have_array()): return 1, 0
    raise_if_none(O.pbs_arrayid_offset, "MY_PBS_ARRAYID_OFFSET")
    raise_if_none(O.pbs_arrayid, "PBS_ARRAYID")
    if (O.mpi_num_procs is None and O.mpi_vpid is None):
      n = O.pbs_array_size
      i = O.pbs_arrayid-1 + O.pbs_arrayid_offset
    else:
      raise_if_none(O.mpi_num_procs, "OMPI_MCA_ns_nds_num_procs")
      raise_if_none(O.mpi_vpid, "OMPI_MCA_ns_nds_vpid")
      n = O.pbs_array_size * O.mpi_num_procs
      i = (O.pbs_arrayid-1 + O.pbs_arrayid_offset) * O.mpi_num_procs + O.mpi_vpid
    return n, i

if (__name__ == "__main__"):
  n,i = chunk_info().show(prefix="*** ", even_if_none=True).as_n_i_pair()
  print "n,i:", n,i
  print "OK"
