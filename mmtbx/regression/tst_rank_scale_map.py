from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.mrcfile
from scitbx.array_family import flex
from libtbx import easy_run
from libtbx.test_utils import approx_equal

pdb_str = """
CRYST1   10.240   10.010   10.990  90.00  90.00  90.00 P1
ATOM      1  N   ASP L   1       5.347   5.804   5.380  1.00 34.60           N
"""

def exercise_00(prefix="tst_rank_scale_map"):
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  fc = xrs.structure_factors(d_min=2).f_calc()
  cs = fc.crystal_symmetry()
  fft_map = fc.fft_map(resolution_factor=0.25)
  m = fft_map.real_map_unpadded()
  iotbx.mrcfile.write_ccp4_map(
    file_name="%s.ccp4"%prefix,
    unit_cell=cs.unit_cell(),
    space_group=cs.space_group(),
    #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
    #gridding_last=n_real,  # This causes a bug (map gets shifted)
    map_data=m,
    labels=flex.std_string([""]))
  assert not easy_run.call("phenix.rank_scale_map %s.ccp4"%prefix)
  #
  ccp4_map = iotbx.mrcfile.map_reader(
    file_name="%s.ccp4_rank_scaled.ccp4"%prefix)
  m = ccp4_map.data.as_double()
  mmm = m.as_1d().min_max_mean().as_tuple()
  assert approx_equal(mmm, [0,1,0.5], 0.01)

if (__name__ == "__main__"):
  exercise_00()
  print("OK")
