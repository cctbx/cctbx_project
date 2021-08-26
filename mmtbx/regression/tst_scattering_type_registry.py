from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out


# Testing behavior of scattering_type_registry

def check_scattering_type_registry():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log = null_out())
  model.process(make_restraints=True)
  xrs1 = model.get_xray_structure()
  xrs2 = model.get_hierarchy().extract_xray_structure(
        crystal_symmetry=model.crystal_symmetry())
  xrs1.scattering_type_registry(table = "electron")
  xrs2.scattering_type_registry(table = "electron")
  xrs1.scattering_type_registry().show()
  xrs2.scattering_type_registry().show()
  # TODO: Assert to the same value once added
  assert (xrs1.scattering_type_registry().gaussian("O1-") is None)
  assert (xrs2.scattering_type_registry().gaussian("O1-") is None)

pdb_str = """
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
ATOM      1  O   ARG M 284     148.445 116.581 116.082  1.00 56.57           O1-
"""

if (__name__ == "__main__"):
  t0 = time.time()
  check_scattering_type_registry()
  print("Finished. Time:", round(time.time()-t0, 2))
  print("OK")
