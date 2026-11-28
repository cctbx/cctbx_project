from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import iotbx.phil
import mmtbx.model
import mmtbx.refinement.real_space.occupancy
from iotbx import map_manager, map_model_manager
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

pdb_str = """\
CRYST1   15.000   16.000   17.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  SER A   1       6.200   5.800   5.000  1.00 10.00           C
ATOM      3  C   SER A   1       7.400   4.900   5.200  1.00 10.00           C
ATOM      4  O   SER A   1       7.300   3.700   5.000  1.00 10.00           O
ATOM      5  CB  SER A   1       6.300   6.800   6.200  1.00 10.00           C
ATOM      6  OG ASER A   1       5.200   7.700   6.200  %4.2f 10.00           O
ATOM      7  OG BSER A   1       7.500   7.500   6.400  %4.2f 10.00           O
ATOM      8  N   GLY A   2       8.500   5.500   5.600  1.00 10.00           N
ATOM      9  CA  GLY A   2       9.700   4.700   5.800  1.00 10.00           C
ATOM     10  C   GLY A   2      10.900   5.600   6.000  1.00 10.00           C
ATOM     11  O   GLY A   2      10.800   6.800   5.800  %4.2f 10.00           O
TER
END
"""

# same definition as in phenix.real_space_refine parameters
occupancies_phil_str = """
individual = None
  .type = str
  .multiple = True
remove_selection = None
  .type = str
  .multiple = True
constrained_group
  .multiple = True
{
  selection = None
    .type = str
    .multiple = True
}
"""

def get_model(occ_a, occ_b, occ_o, table="wk1995", b_add=0):
  pdb_inp = iotbx.pdb.input(
    source_info=None, lines=pdb_str%(occ_a, occ_b, occ_o))
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.setup_scattering_dictionaries(scattering_table=table, d_min=1.0)
  if(b_add != 0): model.set_b_iso(values = model.get_b_iso() + b_add)
  return model

def run_one(params, d_min=1.0, table="wk1995", b_add=0):
  # Map: altloc occupancies 0.7/0.3, carbonyl O of GLY 2 at 0.6
  fc = get_model(0.7, 0.3, 0.6, table, b_add).get_xray_structure(
    ).structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  # Start: 0.5/0.5 and 1
  model = get_model(0.5, 0.5, 1.0, table, b_add)
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = True)
  mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)
  mmtbx.refinement.real_space.occupancy.ncs_aware_refinement(
    map_model_manager = mmm,
    d_min             = d_min,
    atom_radius       = 3.0,
    occupancy_params  = params,
    log               = null_out())
  return mmm.model().get_hierarchy().atoms().extract_occ()

def exercise():
  """
  Real-space occupancy refinement. Default: only alternative conformations are
  refined, constrained to sum to one. User selections: an individual atom is
  refined too; an atom listed in remove_selection is not refined.
  """
  master = iotbx.phil.parse(occupancies_phil_str)
  # default
  occ = run_one(params = master.extract())
  print("default:   ", " ".join(["%4.2f"%o for o in occ]))
  assert approx_equal(occ[5], 0.7, 0.05)
  assert approx_equal(occ[6], 0.3, 0.05)
  assert approx_equal(occ[5]+occ[6], 1.0, 1.e-3)
  sel = flex.bool([True]*5+[False,False]+[True]*4)
  assert occ.select(sel).all_eq(1.0)
  # individual
  params = master.fetch(iotbx.phil.parse(
    "individual = chain A and resseq 2 and name O")).extract()
  occ = run_one(params = params)
  print("individual:", " ".join(["%4.2f"%o for o in occ]))
  assert approx_equal(occ[10], 0.6, 0.07)
  assert approx_equal(occ[5], 0.7, 0.05)
  assert approx_equal(occ[5]+occ[6], 1.0, 1.e-3)
  assert occ[:5].all_eq(1.0) and occ[7:10].all_eq(1.0)
  # remove_selection
  params = master.fetch(iotbx.phil.parse(
    "remove_selection = chain A and resseq 1")).extract()
  occ = run_one(params = params)
  print("remove:    ", " ".join(["%4.2f"%o for o in occ]))
  assert approx_equal(occ[5], 0.5) and approx_equal(occ[6], 0.5)
  assert occ[:5].all_eq(1.0) and occ[7:].all_eq(1.0)

def exercise_low_resolution():
  """
  Electron scattering table, 2 and 3 A, larger B: occupancies are still
  recovered. They were not when the model lost its scattering table on the
  way and the mask around atoms cut off a part of their density (the latter
  is what adp.fill_outside_mask=False does).
  """
  import mmtbx.refinement.real_space.adp as adp
  master = iotbx.phil.parse(occupancies_phil_str)
  params = master.fetch(iotbx.phil.parse(
    "individual = chain A and resseq 2 and name O")).extract()
  for d_min, b_add, eps in [(2., 20, 0.04), (3., 40, 0.07)]:
    occ = run_one(params = params, d_min = d_min, table = "electron",
      b_add = b_add)
    print("d_min %3.1f:"%d_min, " ".join(["%4.2f"%o for o in occ]))
    assert approx_equal(occ[5], 0.7, eps)
    assert approx_equal(occ[5]+occ[6], 1.0, 1.e-3)
    assert approx_equal(occ[10], 0.6, eps)
  adp.fill_outside_mask = False
  try:
    occ = run_one(params = params, d_min = 3., table = "electron", b_add = 40)
  finally:
    adp.fill_outside_mask = True
  print("d_min 3.0, no fill:", " ".join(["%4.2f"%o for o in occ]))
  assert abs(occ[10] - 0.6) > 0.15, occ[10]
  assert abs(occ[5] - 0.7) > 0.07, occ[5]

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  exercise_low_resolution()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
