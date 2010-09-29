from libtbx import Auto
import sys

class chunk_manager(object):

  def __init__(O, n, i):
    assert n > 0
    assert i >= 0
    assert i < n
    O.n = n
    O.i = i
    O.queuing_system_info = None

  def easy_all(O, log_format=Auto, out=Auto):
    O.queuing_system_overrides_chunk()
    O.redirect_chunk_stdout_and_stderr(log_format=log_format, out=out)
    return O

  def skip_iteration(O, i):
    return (i % O.n != O.i)

  def queuing_system_overrides_chunk(O):
    from libtbx.queuing_system_utils import pbs_utils, sge_utils
    pbs_info = pbs_utils.chunk_info()
    sge_info = sge_utils.info()
    assert [pbs_info, sge_info].count(None) <= 1
    if (pbs_info.have_array()):
      O.queuing_system_info = pbs_info
      n, i = pbs_info.as_n_i_pair()
      O.n = max(O.n, n)
      O.i = i
    elif (sge_info.have_array()):
      O.queuing_system_info = sge_info
      O.n = max(O.n, sge_info.last)
      O.i = sge_info.id - 1
    return O

  def redirect_chunk_stdout_and_stderr(O,
        log_format=Auto,
        out=Auto,
        have_array=False):
    if (O.n == 1): return
    log_name = None
    if (not have_array):
      i = O.queuing_system_info
      if (i is not None and i.have_array()):
        have_array = True
    if (have_array):
      if (log_format is Auto): log_format="log%%0%dd"
      fmt = log_format % max(3, len("%d" % (O.n-1)))
      log_name = fmt % O.i
      log = open(log_name, "w")
      sys.stdout = log
      sys.stderr = log
    from libtbx.utils import host_and_user
    if (out is Auto): out = sys.stdout
    if (out is not None):
      host_and_user().show(out=out)
      print >> out, "chunk.n:", O.n
      print >> out, "chunk.i:", O.i
      if (log_name is not None):
        print >> out, "log_name:", log_name
      print >> out
    return O
