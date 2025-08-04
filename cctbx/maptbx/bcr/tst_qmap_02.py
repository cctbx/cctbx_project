from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.map_manager
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

def run(d_min=2):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  cs = pdb_inp.crystal_symmetry()
  pdb_inp.write_pdb_file(file_name="m.pdb")
  xrs = pdb_inp.xray_structure_simple()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.3
    #pre_determined_n_real = (32, 40, 48)
    )
  n_real = crystal_gridding.n_real()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(
        table = "wk1995",
        d_min = 2,
        types_without_a_scattering_contribution=["?"])

  fc = xrs.structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(crystal_gridding = crystal_gridding)
  m1 = fft_map.real_map_unpadded()
  OmegaMap, _, _, _ = qmap.compute(xray_structure=xrs,
                     n_real=n_real, resolution=d_min, resolutions=None,
                     debug=True)
  cc = flex.linear_correlation(
    x=m1.as_1d(), y=OmegaMap.as_1d()).coefficient()
  print(cc)
  #assert cc > 0.99, cc
  mm1 = iotbx.map_manager.map_manager(
    map_data                   = m1,
    unit_cell_grid             = m1.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mm1.write_map("m_fft.ccp4")
  mm2 = iotbx.map_manager.map_manager(
    map_data                   = OmegaMap,
    unit_cell_grid             = OmegaMap.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mm2.write_map("m_bcr.ccp4")
  cc = flex.linear_correlation(
    x=mm1.map_data().as_1d(),
    y=mm2.map_data().as_1d()).coefficient()
  print(cc)
  #assert cc > 0.99

if (__name__ == "__main__"):
  start = time.perf_counter()
  run()
  print("Time:", time.perf_counter()-start)
  print("OK")
