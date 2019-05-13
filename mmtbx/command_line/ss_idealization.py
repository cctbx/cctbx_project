from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.ss_idealization

from mmtbx.secondary_structure import build as ssb
import mmtbx.model
import iotbx.pdb
import os, sys



def run(args):
  pdb_inp = iotbx.pdb.input(source_info=None,
    file_name=args[0])
  model = mmtbx.model.manager(
      model_input=pdb_inp)
  params = ssb.master_phil.extract()
  params.ss_idealization.file_name_before_regularization="before_reg.pdb"
  params.ss_idealization.enabled=True
  rm = ssb.substitute_ss(
    model = model,
    params = params,
    log=sys.stdout)
  out_fname = "%s_ss_ideal.pdb" % os.path.basename(args[0])
  txt = model.model_as_pdb()
  with open(out_fname, 'w') as f:
    f.write(txt)
  print("File saved: %s" % out_fname)
  print("All done.")

if __name__ == "__main__" :
  run(sys.argv[1:])
