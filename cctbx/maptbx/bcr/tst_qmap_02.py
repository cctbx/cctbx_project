from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.map_manager
import iotbx.pdb
from cctbx.maptbx.bcr import qmap

from cctbx import maptbx

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

#pdb_str = """
#CRYST1   20.000   22.000   24.000  90.00  90.00  90.00 P 1
#ATOM      1  N   HIS A 109       9.500  14.000   9.000  1.10 60.10           N
#ATOM      2  C   HIS A 109      20.000   0.000   0.000  0.90 20.00           C
#ATOM      3  O   HIS A 109       0.000   0.000   0.000  0.50 40.00           O
#TER
#END
#"""

pdb_str = """
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A 109      10.000  14.000   9.000  1.10 30.10           N
ATOM      2  C   HIS A 109      13.000  10.000  10.000  0.90 20.00           C
ATOM      3  O   HIS A 109      12.000  12.000  11.000  0.50 10.00           O
TER
END
"""

def run(d_min=2, prefix="12"):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  cs = pdb_inp.crystal_symmetry()
  hierarchy = pdb_inp.construct_hierarchy()
  pdb_inp.write_pdb_file(file_name="m%s.pdb"%prefix)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.7
    #pre_determined_n_real = (32, 40, 48)
    )
  n_real = crystal_gridding.n_real()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(
        table = "it1992",
        d_min = 2,
        types_without_a_scattering_contribution=["?"])

  fc = xrs.structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(crystal_gridding = crystal_gridding)
  m1 = fft_map.real_map_unpadded()
  OmegaMap = qmap.compute(hierarchy=hierarchy, unit_cell=cs.unit_cell(),
                     n_real=n_real, resolution=d_min, resolutions=None,
                     debug=True)
  print(flex.linear_correlation(x=m1.as_1d(), y=OmegaMap.as_1d()).coefficient())
  mm1 = iotbx.map_manager.map_manager(
    map_data                   = m1,
    unit_cell_grid             = m1.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mm1.write_map("m_fft_%s.ccp4"%prefix)
  mm2 = iotbx.map_manager.map_manager(
    map_data                   = OmegaMap,
    unit_cell_grid             = OmegaMap.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mm2.write_map("m_bcr_%s.ccp4"%prefix)

if (__name__ == "__main__"):
  run()
