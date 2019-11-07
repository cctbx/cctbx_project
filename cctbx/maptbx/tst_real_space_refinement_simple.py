from __future__ import absolute_import, division, print_function
from cctbx.maptbx import real_space_refinement_simple
import cctbx.geometry_restraints.manager
from cctbx import xray
from cctbx.array_family import flex
from scitbx.graph import test_cases_tardy_pdb
import scitbx.math
from scitbx import matrix
from libtbx.utils import null_out, format_cpu_times
import random
import sys
from six.moves import zip

def exercise_lbfgs(test_case, use_geo, out, d_min=2):
  sites_cart, geo_manager = cctbx.geometry_restraints.manager.\
    construct_non_crystallographic_conserving_bonds_and_angles(
      sites_cart=flex.vec3_double(test_case.sites),
      edge_list_bonds=test_case.bonds,
      edge_list_angles=test_case.angles())
  scatterers = flex.xray_scatterer(
    sites_cart.size(), xray.scatterer(scattering_type="C", b=20))
  for sc,lbl in zip(scatterers, test_case.labels):
    sc.label = lbl
  structure = xray.structure(
    crystal_symmetry=geo_manager.crystal_symmetry,
    scatterers=scatterers)
  structure.set_sites_cart(sites_cart=sites_cart)
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  fft_map = f_calc.fft_map()
  fft_map.apply_sigma_scaling()
  if (use_geo):
    axis = matrix.col(flex.random_double_point_on_sphere())
    rot = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis=axis, angle=25, deg=True)
    trans = matrix.col(flex.random_double_point_on_sphere()) * 1.0
    structure.apply_rigid_body_shift(rot=rot, trans=trans)
    geo_manager.energies_sites(sites_cart=structure.sites_cart()).show(f=out)
    minimized = real_space_refinement_simple.lbfgs(
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      geometry_restraints_manager=geo_manager,
      gradients_method="fd",
      real_space_target_weight=1,
      real_space_gradients_delta=d_min/3)
    geo_manager.energies_sites(sites_cart=minimized.sites_cart).show(f=out)
  else:
    minimized = real_space_refinement_simple.lbfgs(
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      unit_cell=structure.unit_cell(),
      gradients_method="fd",
      real_space_gradients_delta=d_min/3)
  rmsd_start = sites_cart.rms_difference(structure.sites_cart())
  rmsd_final = sites_cart.rms_difference(minimized.sites_cart)
  print("RMSD start, final:", rmsd_start, rmsd_final, file=out)
  if (use_geo):
    assert rmsd_start >= 1-1e-6
    assert rmsd_final < 0.2
  def show_f_g(label, f, g):
    print(label, "f, |g|:", f, flex.mean_sq(g)**0.5, file=out)
  show_f_g(label="start", f=minimized.f_start, g=minimized.g_start)
  show_f_g(label="final", f=minimized.f_final, g=minimized.g_final)
  assert minimized.f_final <= minimized.f_start
  return minimized

def run(args):
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  out = null_out()
  remaining_args = []
  for arg in args:
    if (arg == "--verbose"): out = sys.stdout
    else: remaining_args.append(arg)
  test_cases = test_cases_tardy_pdb.select_test_cases(
    tags_or_indices=remaining_args)
  for test_case in test_cases:
    print("test case %d: %s" % (test_case.index, test_case.tag), file=out)
    minimized = []
    for use_geo in [False, True]:
      minimized.append(exercise_lbfgs(test_case, use_geo=use_geo, out=out))
    m0, m1 = minimized
    assert m0.real_space_target_weight == 1
    assert m1.real_space_target_weight == 1
    assert m1.f_final < m0.f_start * 0.98
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
