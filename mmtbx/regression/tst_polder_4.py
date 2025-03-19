from __future__ import division, print_function
import iotbx.pdb
import time, os
from libtbx import easy_run
import iotbx.pdb

def exercise_00(prefix="tst_phenix_regression_maps_polder"):
  """
  Make sure EXACTLY THIS comand below runs (= does not crash):

  phenix.polder model.pdb data.mtz scattering_table=electron selectioon="resseq 2"
  """
  pdb_str= """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 10.00           O
HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O
END
"""
  pdb_inp = iotbx.pdb.input(lines = pdb_str, source_info=None)
  pdb_inp.write_pdb_file(file_name="%s.pdb"%prefix)
  xrs = pdb_inp.xray_structure_simple()
  fo = abs(xrs.structure_factors(d_min=1).f_calc())
  mtz_dataset = fo.as_mtz_dataset(column_root_label = "FOBS")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "%s.mtz"%prefix)
  cmd = " ".join([
    "phenix.polder",
    "%s.mtz"%prefix,
    "%s.pdb"%prefix,
    "selection='resseq 1' ",
    "scattering_table=electron",
    "> %s.zlog"%prefix
    ])
  #print(cmd)
  assert not easy_run.call(cmd)
  assert os.path.isfile("tst_phenix_regression_maps_polder_polder_map_coeffs.mtz")

  os.remove('tst_phenix_regression_maps_polder.mtz')
  os.remove('tst_phenix_regression_maps_polder.pdb')
  os.remove('tst_phenix_regression_maps_polder.zlog')
  os.remove('tst_phenix_regression_maps_polder_polder_map_coeffs.mtz')

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  print("Time: %6.4f" % (time.time()-t0))
