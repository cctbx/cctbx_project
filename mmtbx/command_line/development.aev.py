from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.aev
import sys
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import mmtbx.atomic_environment_vectors as aev

def main(filename, precision):
  t0 = time.time()
  pdb_inp = iotbx.pdb.input(file_name = filename)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    build_grm     = True,
    log           = null_out())
  a = aev.AEV(model = model)
  b = aev.compare(a)
  print(b)
  recs = aev.format_HELIX_records_from_AEV(b, float(precision))
  print("\n".join(recs))
  print('time', time.time()-t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
