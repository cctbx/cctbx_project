from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx
import time

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

pdb_str = """
remark CRYST1   20.000   25.000   30.000  90.00  90.00  90.00 P 1
CRYST1   20.000   25.000   30.000  96.06  99.11  91.14 P 1
ATOM      1  N   HIS A 109      10.000  14.000   9.000  1.10 30.10           N
ATOM      2  C   HIS A 109      13.000  10.000  10.000  0.90 20.00           C
ATOM      3  O   HIS A 109      12.000  12.000  11.000  0.50 10.00           O
TER
END
"""

def run(debug, use_exp_table, table, d_min=2):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  cs = pdb_inp.crystal_symmetry()
  pdb_inp.write_pdb_file(file_name="m.pdb")
  xrs = pdb_inp.xray_structure_simple()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = d_min/2.
    )
  n_real = crystal_gridding.n_real()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(
    table = table,
    d_min = d_min,
    types_without_a_scattering_contribution=["?"])
  # FFT
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(crystal_gridding = crystal_gridding)
  mFFT = fft_map.real_map_unpadded()
  # VRM
  o = qmap.compute(xray_structure=xrs,
                   n_real=n_real, resolution=d_min, resolutions=None,
                   use_exp_table=use_exp_table,
                   debug=debug)
  mVRM = o.map_data()
  #
  cc = flex.linear_correlation(
    x=mFFT.as_1d(), y=mVRM.as_1d()).coefficient()
  assert cc > 0.99, cc
  mmFFT = iotbx.map_manager.map_manager(
    map_data                   = mFFT,
    unit_cell_grid             = mFFT.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mmFFT.write_map("m_fft.ccp4")
  mmVRM = iotbx.map_manager.map_manager(
    map_data                   = mVRM,
    unit_cell_grid             = mVRM.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mmVRM.write_map("m_bcr.ccp4")
  cc = flex.linear_correlation(
    x=mmFFT.map_data().as_1d(),
    y=mmVRM.map_data().as_1d()).coefficient()
  assert cc > 0.99

if (__name__ == "__main__"):
  start = time.perf_counter()
  for table in ["electron", "wk1995"]:
    for it in [[True, False], [False, True]]:
      debug, use_exp_table = it
      run(debug=debug, use_exp_table=use_exp_table, table=table)
  print("Time:", time.perf_counter()-start)
  print("OK")
