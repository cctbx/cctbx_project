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

def exercise_lbfgs(test_case, use_geo, d_min=2):
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
  if (use_geo):
    structure.shake_sites_in_place(rms_difference=0.1)
    geo_manager.energies_sites(sites_cart=structure.sites_cart()).show()
    minimized = real_space_refinement_simple.lbfgs(
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      gradients_delta=d_min/3,
      geometry_restraints_manager=geo_manager,
      real_space_weight=0.01)
    structure.set_sites_cart(sites_cart=minimized.sites_cart)
    geo_manager.energies_sites(sites_cart=structure.sites_cart()).show()
  else:
    minimized = real_space_refinement_simple.lbfgs(
      sites_cart=structure.sites_cart(),
      density_map=fft_map.real_map(),
      gradients_delta=d_min/3,
      unit_cell=structure.unit_cell())
  print "f, |g|, w:", \
    minimized.f, flex.mean_sq(minimized.g)**0.5, minimized.rs_weight
  print
  return minimized

def run(args):
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  test_cases = tst_tardy_pdb.select_test_cases(tags_or_indices=args)
  for test_case in test_cases:
    print "test case %d: %s" % (test_case.index, test_case.tag)
    minimized = []
    for use_geo in [False, True]:
      minimized.append(exercise_lbfgs(test_case, use_geo=use_geo))
    f0 =  minimized[0].f
    f1 =  minimized[1].f / minimized[1].rs_weight
    assert f1 > f0 * 1.1
    assert f1 < f0 * 0.7
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
