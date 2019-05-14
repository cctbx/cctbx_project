
from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
import os.path

def exercise():
  from mmtbx.wwpdb import utils as wwpdb_utils
  for fn in ["3qgg.pdb", "3gqq.mtz", "3gqq-sf.cif"] :
    if (os.path.isfile(fn)) : os.remove(fn)
  pdb_id = "3gqq"
  wwpdb_utils.fetch_pdb_data("3gqq")
  program, program_full = wwpdb_utils.get_program("3gqq.pdb")
  assert (program == "PHENIX.REFINE")
  data = wwpdb_utils.find_data_arrays("3gqq.mtz", log=null_out())
  filter = wwpdb_utils.filter_pdb_file("3gqq.pdb", log=null_out())
  assert (filter.n_semet == 24) and (filter.n_unknown == 30)
  print("OK")

if (__name__ == "__main__"):
  exercise()
