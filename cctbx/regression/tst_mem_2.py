from __future__ import absolute_import, division, print_function
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run
from iotbx import mtz

def run():
  pdb_str="""
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  5.00           C
END
"""
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  fc = xrs.structure_factors(d_min = 1, algorithm = "direct",
    anomalous_flag=True).f_calc()
  mtz_dataset = fc.as_mtz_dataset(column_root_label="F")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "tst_mem_2.mtz")
  assert fc.anomalous_flag()
  cmd = "phenix.maximum_entropy_map tst_mem_2.mtz > tst_mem_2.log"
  assert not easy_run.call(cmd)
  o = mtz.object(file_name="tst_mem_2_mem.mtz")
  mas = o.as_miller_arrays()
  assert len(mas) == 2
  for ma in mas:
    assert ma.anomalous_flag() is False

if __name__ == "__main__" :
  if libtbx.env.has_module("phenix"):
    run()
    print("OK")
  else:
    print("Skipping test: phenix not available")
