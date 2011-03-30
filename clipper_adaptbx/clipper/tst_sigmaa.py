import clipper
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys

if (1):
  flex.set_random_seed(2)

def exercise_SFweight_spline_core(structure, d_min, verbose=0):
  structure.scattering_type_registry(d_min=d_min)
  f_obs = abs(structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc())
  if (0 or verbose):
    f_obs.show_summary()
  f_obs = miller.array(
    miller_set=f_obs,
    data=f_obs.data(),
    sigmas=flex.sqrt(f_obs.data()))
  partial_structure = xray.structure(
    crystal_symmetry=structure,
    scatterers=structure.scatterers()[:-2])
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=partial_structure).f_calc()
  test_set_flags = (flex.random_double(size=f_obs.indices().size()) < 0.1)
  sfweight = clipper.SFweight_spline_interface(
    unit_cell=f_obs.unit_cell(),
    space_group=f_obs.space_group(),
    miller_indices=f_obs.indices(),
    anomalous_flag=f_obs.anomalous_flag(),
    f_obs_data=f_obs.data(),
    f_obs_sigmas=f_obs.sigmas(),
    f_calc=f_calc.data(),
    test_set_flags=test_set_flags,
    n_refln=f_obs.indices().size()//10,
    n_param=20)
  if (0 or verbose):
    print "number_of_spline_parameters:",sfweight.number_of_spline_parameters()
    print "mean fb: %.8g" % flex.mean(flex.abs(sfweight.fb()))
    print "mean fd: %.8g" % flex.mean(flex.abs(sfweight.fd()))
    print "mean phi: %.8g" % flex.mean(sfweight.centroid_phases())
    print "mean fom: %.8g" % flex.mean(sfweight.figures_of_merit())
  return sfweight

def exercise_with_fixed_structure():
  structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(46.7058, 46.7058, 79.3998, 90, 90, 120),
      space_group_symbol="P 31"),
    scatterers=flex.xray_scatterer(
      [xray.scatterer(scattering_type="const", site=site) for site in [
        (0.0169, 0.8953, 0.1115),
        (0.9395, 0.1282, 0.1780),
        (0.2998, 0.3497, 0.0593),
        (0.8220, 0.8814, 0.1601),
        (0.6478, 0.4879, 0.3141)]]))
  sfweight = exercise_SFweight_spline_core(
    structure=structure, d_min=5, verbose="--Verbose" in sys.argv[1:])
  assert approx_equal(flex.mean(flex.abs(sfweight.fb())), 1.7545459)
  assert approx_equal(flex.mean(flex.abs(sfweight.fd())), 1.8437204)
  assert approx_equal(flex.mean(sfweight.centroid_phases()), -0.033979132)
  assert approx_equal(flex.mean(sfweight.figures_of_merit()), 0.018943642)

def exercise_SFweight_spline(
      space_group_info,
      n_scatterers=10,
      d_min=4,
      verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["C"]*n_scatterers,
    volume_per_atom=1000)
  sfweight = exercise_SFweight_spline_core(
    structure=structure, d_min=d_min, verbose=verbose)

def run_call_back(flags, space_group_info):
  exercise_SFweight_spline(space_group_info, verbose=flags.Verbose)

def run():
  structure = exercise_with_fixed_structure()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
