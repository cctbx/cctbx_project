from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
from mmtbx.refinement.real_space import individual_sites
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import maptbx
from cctbx import uctbx
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

"""
Real-space refinement of coordinates as shifts from the start position with
L-BFGS-B, each shift bounded (individual_sites.simple with max_shift, the way
whole models are refined), compared with refinement of coordinates
themselves with L-BFGS (max_shift=None).
"""

def get_inputs(d_min = 2.0, shake = 0.3):
  """
  Helix with geometry restraints, map calculated from it, shaken coordinates.
  """
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWY")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 5)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.process(make_restraints = True)
  xrs = model.get_xray_structure()
  fft_map = xrs.structure_factors(d_min = d_min).f_calc().fft_map(
    resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  xrs_shaken = xrs.deep_copy_scatterers()
  flex.set_random_seed(0)
  xrs_shaken.shake_sites_in_place(mean_distance = shake)
  return xrs, xrs_shaken, fft_map.real_map_unpadded(), \
    model.get_restraints_manager().geometry

def refine(xrs, map_data, grm, max_shift, selection=None, weight=20,
           max_iterations=150, gradients_method="tricubic"):
  if(selection is None): selection = flex.bool(xrs.scatterers().size(), True)
  o = individual_sites.simple(
    target_map                  = map_data,
    selection                   = selection,
    geometry_restraints_manager = grm,
    gradients_method            = gradients_method,
    max_iterations              = max_iterations,
    max_shift                   = max_shift)
  o.refine(weight = weight, xray_structure = xrs.deep_copy_scatterers())
  return o

def exercise_default():
  """
  Whole-model refinement (easy, weight.run) refines bounded shifts, with a
  bound of 10 A; simple (used for fragments too) does not unless asked.
  """
  import inspect
  import mmtbx.refinement.real_space.weight
  assert individual_sites.default_max_shift == 10.0
  def default(f):
    s = inspect.signature(f)
    return s.parameters["max_shift"].default
  assert default(individual_sites.easy.__init__) == 10.0
  assert default(mmtbx.refinement.real_space.weight.run.__init__) == 10.0
  assert default(individual_sites.simple.__init__) is None

def exercise_convergence():
  """
  Both ways bring shaken model back to the answer, to about the same
  distance from it.
  """
  xrs, xrs_shaken, map_data, grm = get_inputs()
  d_start = flex.mean(xrs.distances(xrs_shaken))
  assert d_start > 0.25
  d = {}
  for max_shift in [None, 5.0]:
    o = refine(xrs_shaken, map_data, grm, max_shift = max_shift)
    d[max_shift] = flex.mean(flex.sqrt((o.sites_cart()-xrs.sites_cart()).dot()))
    b, a = o.rmsds()
    print("max_shift=%s: distance to answer %6.3f -> %6.3f bonds %6.4f "
      "angles %5.2f largest shift %s"%(str(max_shift), d_start,
      d[max_shift], b, a, str(o.largest_shift)))
    assert d[max_shift] < 0.15, d
    assert b < 0.02 and a < 2.5
  assert abs(d[None]-d[5.0]) < 0.03, d
  # largest shift is reported and is far below the bound
  assert 0.3 < o.largest_shift < 2, o.largest_shift

def exercise_bound():
  """
  No coordinate moves by more than the bound; with a bound smaller than the
  shifts needed some coordinates end at the bound.
  """
  xrs, xrs_shaken, map_data, grm = get_inputs()
  start = xrs_shaken.sites_cart().as_double()
  for max_shift in [0.05, 0.2]:
    o = refine(xrs_shaken, map_data, grm, max_shift = max_shift)
    shifts = flex.abs(o.sites_cart().as_double() - start)
    assert flex.max(shifts) <= max_shift + 1.e-9, [max_shift, flex.max(shifts)]
    assert approx_equal(o.largest_shift, flex.max(shifts), 1.e-9)
    n_at_bound = (shifts > max_shift - 1.e-6).count(True)
    print("max_shift=%4.2f: coordinates at the bound: %d of %d"%(
      max_shift, n_at_bound, shifts.size()))
    assert n_at_bound > 5, n_at_bound
  # with the unbounded minimizer the same coordinates move further
  o = refine(xrs_shaken, map_data, grm, max_shift = None)
  shifts = flex.abs(o.sites_cart().as_double() - start)
  assert flex.max(shifts) > 0.3

def exercise_selection():
  """
  Atoms that are not selected for refinement do not move, and bounds apply
  to those that are.
  """
  xrs, xrs_shaken, map_data, grm = get_inputs()
  n = xrs.scatterers().size()
  selection = flex.bool([i%2 == 0 for i in range(n)])
  start = xrs_shaken.sites_cart()
  o = refine(xrs_shaken, map_data, grm, max_shift = 0.2,
    selection = selection)
  d = flex.sqrt((o.sites_cart() - start).dot())
  assert flex.max(d.select(~selection)) < 1.e-9
  assert flex.min(d.select(selection)) > 1.e-4
  shifts = flex.abs((o.sites_cart() - start).as_double())
  assert flex.max(shifts) <= 0.2 + 1.e-9

def exercise_gradients_methods():
  """
  Every gradients method works with bounded shifts.
  """
  xrs, xrs_shaken, map_data, grm = get_inputs()
  for method in ["fd", "linear", "tricubic"]:
    o = refine(xrs_shaken, map_data, grm, max_shift = 5.0,
      gradients_method = method)
    d = flex.mean(flex.sqrt((o.sites_cart()-xrs.sites_cart()).dot()))
    assert d < 0.2, [method, d]

def exercise_low_level():
  """
  maptbx.real_space_refinement_simple.lbfgs with shifts_bound: no
  restraints, map only; atoms move uphill in the map and stay within bounds.
  """
  xrs, xrs_shaken, map_data, grm = get_inputs()
  sites_start = xrs_shaken.sites_cart()
  import scitbx.lbfgs
  for bound in [0.1, 5.0]:
    o = maptbx.real_space_refinement_simple.lbfgs(
      sites_cart                 = sites_start.deep_copy(),
      density_map                = map_data,
      gradients_method           = "tricubic",
      unit_cell                  = xrs.unit_cell(),
      real_space_gradients_delta = 0.25,
      lbfgs_termination_params   = scitbx.lbfgs.termination_parameters(
        max_iterations = 50),
      shifts_bound               = bound)
    assert o.f_final < o.f_start
    shifts = flex.abs((o.sites_cart_variable - sites_start).as_double())
    assert flex.max(shifts) <= bound + 1.e-9
    assert flex.max(shifts) > 0.05
    assert approx_equal(o.max_shift, flex.max(shifts), 1.e-9)
  assert sites_start.size() == xrs.scatterers().size()

class restraints_manager(object):
  """ What easy() needs: an object with .geometry """
  def __init__(self, geometry): self.geometry = geometry

class states(object):
  def __init__(self): self.n = 0
  def add(self, sites_cart): self.n += 1

def exercise_easy_special_cases():
  """
  easy(): nothing selected for refinement - nothing moves (it used to crash
  the interpreter with bounded shifts); states are collected also when the
  weight is searched for.
  """
  import iotbx.pdb
  xrs, xrs_shaken, map_data, grm = get_inputs()
  n = xrs.scatterers().size()
  ph = iotbx.pdb.input(source_info = None, lines = flex.split_lines(
    xrs_shaken.as_pdb_file())).construct_hierarchy()
  for w in [None, 5.]:
    o = individual_sites.easy(
      map_data                    = map_data,
      xray_structure              = xrs_shaken.deep_copy_scatterers(),
      pdb_hierarchy               = ph.deep_copy(),
      geometry_restraints_manager = restraints_manager(grm),
      selection                   = flex.bool(n, False),
      w                           = w)
    d = flex.sqrt((o.xray_structure.sites_cart()-xrs_shaken.sites_cart()).dot())
    assert flex.max(d) < 1.e-9, [w, flex.max(d)]
  o = individual_sites.simple(
    target_map = map_data, selection = flex.bool(n, False),
    geometry_restraints_manager = grm, max_shift = 10.)
  o.refine(weight = 5, xray_structure = xrs_shaken.deep_copy_scatterers())
  assert approx_equal(o.sites_cart(), xrs_shaken.sites_cart(), 1.e-9, out=None)
  for w in [None, 5.]:
    s = states()
    o = individual_sites.easy(
      map_data                    = map_data,
      xray_structure              = xrs_shaken.deep_copy_scatterers(),
      pdb_hierarchy               = ph.deep_copy(),
      geometry_restraints_manager = restraints_manager(grm),
      selection                   = flex.bool([i%2 == 0 for i in range(n)]),
      states_accumulator          = s,
      max_iterations              = 50,
      w                           = w)
    assert s.n > 20, [w, s.n]

if(__name__ == "__main__"):
  t0 = time.time()
  exercise_easy_special_cases()
  exercise_default()
  exercise_convergence()
  exercise_bound()
  exercise_selection()
  exercise_gradients_methods()
  exercise_low_level()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
