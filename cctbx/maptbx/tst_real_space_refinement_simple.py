from __future__ import division
from cctbx.maptbx import real_space_refinement_simple
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
from cctbx import xray
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx import uctbx
from cctbx.array_family import flex
from scitbx.graph import tst_tardy_pdb
import scitbx.math
from scitbx import matrix
from libtbx.utils import null_out, format_cpu_times
import random
import sys

def construct_geometry_restraints_manager(test_case):
  sites = test_case.sites
  bond_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=None)
  for edge_list,weight in [(test_case.bonds, 100), (test_case.angles(), 50)]:
    for i,j in edge_list:
      distance = abs(sites[i] - sites[j])
      bond_proxies.process(geometry_restraints.bond_simple_proxy(
        i_seqs=(i,j), distance_ideal=distance, weight=weight))
  bond_params_table = geometry_restraints.extract_bond_params(
    n_seq=len(sites),
    bond_simple_proxies=bond_proxies.simple)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=flex.vec3_double(sites),
    buffer_layer=5)
  asu_mappings = box.crystal_symmetry().special_position_settings() \
    .asu_mappings(
      buffer_thickness=5,
      sites_cart=box.sites_cart)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  geometry_restraints.add_pairs(bond_asu_table, bond_proxies.simple)
  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
    pair_asu_table=bond_asu_table,
    max_shell=3)
  shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
    for shell_asu_table in shell_asu_tables]
  nonbonded_types = flex.std_string(bond_params_table.size(), "Default")
  nonbonded_params = geometry_restraints.nonbonded_params()
  nonbonded_params.distance_table.setdefault(
    "Default")["Default"] = 1.2
  return box.sites_cart, geometry_restraints.manager.manager(
    crystal_symmetry=box.crystal_symmetry(),
    site_symmetry_table=asu_mappings.site_symmetry_table(),
    bond_params_table=bond_params_table,
    shell_sym_tables=shell_sym_tables,
    nonbonded_params=nonbonded_params,
    nonbonded_types=nonbonded_types,
    nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
    max_reasonable_bond_distance=10)

def exercise_lbfgs(test_case, use_geo, out, d_min=2):
  sites_cart, geo_manager = construct_geometry_restraints_manager(
    test_case=test_case)
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
      gradients_delta=d_min/3,
      geometry_restraints_manager=geo_manager,
      real_space_weight=1)
    geo_manager.energies_sites(sites_cart=minimized.sites_cart).show(f=out)
  else:
    minimized = real_space_refinement_simple.lbfgs(
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      gradients_delta=d_min/3,
      unit_cell=structure.unit_cell())
  rmsd_start = sites_cart.rms_difference(structure.sites_cart())
  rmsd_final = sites_cart.rms_difference(minimized.sites_cart)
  print >> out, "RMSD start, final:", rmsd_start, rmsd_final
  if (use_geo):
    assert rmsd_start >= 1-1e-6
    assert rmsd_final < 0.2
  def show_f_g(f, g):
    print >> out, "start f, |g|:", f, flex.mean_sq(g)**0.5
  show_f_g(f=minimized.f_start, g=minimized.g_start)
  show_f_g(f=minimized.f_final, g=minimized.g_final)
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
  test_cases = tst_tardy_pdb.select_test_cases(tags_or_indices=remaining_args)
  for test_case in test_cases:
    print >> out, "test case %d: %s" % (test_case.index, test_case.tag)
    minimized = []
    for use_geo in [False, True]:
      minimized.append(exercise_lbfgs(test_case, use_geo=use_geo, out=out))
    m0, m1 = minimized
    assert m0.rs_weight is None
    assert m1.rs_weight == 1
    assert m1.f_final < m0.f_start * 0.99
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
