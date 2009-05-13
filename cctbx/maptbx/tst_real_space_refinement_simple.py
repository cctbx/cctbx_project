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
import sys
from cctbx.development import random_structure
from cctbx import sgtbx
from cctbx import adptbx
import random

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

def run(args):
  d_min = 3
  test_cases = tst_tardy_pdb.select_test_cases(tags_or_indices=args)
  for test_case in test_cases:
    print "test case %d: %s" % (test_case.index, test_case.tag)
    sites_cart, geo_manager = construct_geometry_restraints_manager(
      test_case=test_case)
    scatterers = flex.xray_scatterer(
      sites_cart.size(), xray.scatterer(scattering_type="C"))
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
    structure.shake_sites_in_place(rms_difference=0.1)
    geo_manager.energies_sites(sites_cart=structure.sites_cart()).show()
    minimized = real_space_refinement_simple.lbfgs(
      geometry_restraints_manager=geo_manager,
      unit_cell=structure.unit_cell(),
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      gradients_delta=d_min/4,
      real_space_weight_scale=0.01)
    structure.set_sites_cart(sites_cart=minimized.sites_cart)
    geo_manager.energies_sites(sites_cart=structure.sites_cart()).show()
    print "f, |g|, w:", \
      minimized.f, flex.mean_sq(minimized.g)**0.5, minimized.rs_weight
    print
  print "OK"

def run_minimization():
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    unit_cell        = (10, 20, 30, 70, 80, 120),
    elements         =(("O","N","C")*10),
    volume_per_atom  = 50,
    min_distance     = 2,
    u_iso            = adptbx.b_as_u(10.),
    use_u_iso        = True)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.5)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  print "Start:", start_error
  map_current, miller_array, crystal_gridding = compute_map(xray_structure = xrs_sh)
  steps1 = [1.0,0.5,0.25,0.1,0.05,0.01,0.001]
  steps1.reverse()
  steps2 = [1.0,0.5,0.25,0.1,0.05,0.01,0.001]
  for step in steps1+steps2:
    minimized = real_space_refinement_simple.minimization(
      xray_structure   = xrs_sh,
      miller_array     = miller_array,
      crystal_gridding = crystal_gridding,
      map_target       = map_target,
      max_iterations   = 500,
      min_iterations   = 25,
      step             = step)
    xrs_sh = minimized.xray_structure
    map_current = minimized.map_current
    final_error = flex.mean(xrs.distances(other = minimized.xray_structure))
    print "Final:", final_error
  assert start_error >= 0.5
  assert final_error < 0.021

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
  run_minimization()
