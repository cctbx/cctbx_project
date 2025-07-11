from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.aev
import sys
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import __init__ as aev

def main(filename):
  t0 = time.time()
  pdb_inp = iotbx.pdb.input(file_name = filename)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    log           = null_out())
  sel = model.selection(string="protein")
  model = model.select(selection=sel)
  model.crystal_symmetry()
  a = aev.AEV(model = model)
  CC_value = aev.compare(a)
  print(CC_value)
  recs = aev.format_HELIX_records_from_AEV(CC_value)
  print("\n".join(recs))
  print('time', time.time()-t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
