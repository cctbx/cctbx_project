from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx
import time
import mmtbx.model
import cctbx.maptbx.bcr
from scitbx.array_family import flex
import mmtbx.refinement.real_space.adp
from libtbx.test_utils import approx_equal, not_approx_equal
import scitbx.minimizers
from cctbx import adptbx

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

pdb_str_good = """
CRYST1   20.000   25.000   30.000  96.06  99.11  91.14 P 1
remark CRYST1   20.000   25.000   30.000  90.00  90.00  90.00 P 1
ATOM      0  N   HIS A 109      10.000  14.000   9.000  1.10 10.00           N
ATOM      1  C   HIS A 109      13.000  10.000  14.000  0.90 15.00           C
ATOM      2  O   HIS A 109      12.000  10.000  11.000  0.50 20.00           O
TER
END
"""

pdb_str_poor = """
CRYST1   20.000   25.000   30.000  96.06  99.11  91.14 P 1
remark CRYST1   20.000   25.000   30.000  90.00  90.00  90.00 P 1
ATOM      0  N   HIS A 109      10.000  14.000   9.000  1.10  1.10           N
ATOM      1  C   HIS A 109      13.000  10.000  14.000  0.90  2.00           C
ATOM      2  O   HIS A 109      12.000  10.000  11.000  0.50  3.00           O
TER
END
"""

def get_map_data_and_crystal_gridding(resolution, table):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_good)
  model = mmtbx.model.manager(model_input=pdb_inp)
  model.setup_scattering_dictionaries(scattering_table = table)
  xrs = model.get_xray_structure()
  cs = model.crystal_symmetry()
  uc = cs.unit_cell()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.5)
  n_real = crystal_gridding.n_real()
  #
  bcr_scatterers = cctbx.maptbx.bcr.scatterers(
    xray_structure = xrs,
    resolution     = resolution,
    resolutions    = None,
    RadFact        = 2.0,
    RadAdd         = 0.5)
  #
  o = qmap.compute(
    xray_structure = xrs,
    n_real         = n_real,
    resolution     = resolution,
    resolutions    = None,
    use_exp_table  = False,
    debug          = False,
    refinement     = True)
  #
  map_data_inv = o.map_data_inversed()
  o.gradients(map_data = o.map_data_inversed())
  assert approx_equal(o.target(), 0)
  # EXAMPLE how to convert VRM map to Phenix convention
  o = qmap.compute(
    xray_structure = xrs,
    n_real         = n_real,
    resolution     = resolution,
    resolutions    = None,
    use_exp_table  = False,
    debug          = False,
    refinement     = False)
  map_data = o.map_data()
  assert map_data.all()     == (45, 54, 64)
  assert map_data_inv.all() == (64, 54, 45)
  assert not_approx_equal(map_data, map_data_inv)
  nx,ny,nz = map_data.all() # phenix
  new_map = flex.double(flex.grid([nz,ny,nx])) # AU's VRM convention
  for iz in range(0, nz):
    for iy in range(0, ny):
      for ix in range(0, nx):
        new_map[iz,iy,ix] = map_data[ix,iy,iz]
  assert approx_equal(new_map, map_data_inv)
  #
  return o.map_data_inversed(), crystal_gridding, \
    xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)

def get_bcr_scatterers(resolution, table):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_poor)
  model = mmtbx.model.manager(model_input=pdb_inp)
  model.setup_scattering_dictionaries(scattering_table = table)
  xrs = model.get_xray_structure()
  #
  bcr_scatterers = cctbx.maptbx.bcr.scatterers(
    xray_structure = xrs,
    resolution     = resolution,
    resolutions    = None,
    RadFact        = 2.0,
    RadAdd         = 0.5)
  return bcr_scatterers, xrs

def run(resolution=2.0, table="wk1995"):
  # Generate "experimental map"
  map_data, crystal_gridding, answer = get_map_data_and_crystal_gridding(
    resolution = resolution, table = table)
  # Read in model to refine
  bcr_scatterers, xrs = get_bcr_scatterers(
    resolution = resolution, table = table)
  vrm = qmap.compute(
    xray_structure = xrs,
    bcr_scatterers = bcr_scatterers,
    n_real         = crystal_gridding.n_real(),
    resolution     = resolution,
    resolutions    = None,
    use_exp_table  = False,
    debug          = False,
    refinement     = True)
  #
  for it in [1,]:
    x = xrs.extract_u_iso_or_u_equiv()
    lower = flex.double(x.size(), 0)
    upper = flex.double(x.size(), 3) # 1 = 78.95683520871486
    calculator = mmtbx.refinement.real_space.adp.tg_vrm(
      vrm               = vrm,
      map_data          = map_data,
      x                 = x,
      restraints_weight = None,
      bound_flags       = flex.int(x.size(), 2),
      lower_bound       = lower,
      upper_bound       = upper)
    m = scitbx.minimizers.lbfgs(
              mode='lbfgsb', max_iterations=100, calculator=calculator)
  #
  result = xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
  assert approx_equal(answer, result)

if (__name__ == "__main__"):
  start = time.perf_counter()
  run()
  print("Time:", time.perf_counter()-start)
  print("OK")
