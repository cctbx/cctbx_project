import iotbx.pdb
from iotbx.pymol import pml_stick, pml_write
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
import cctbx.geometry_restraints.lbfgs
from cctbx import xray
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx import sgtbx
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx import matrix
from libtbx.itertbx import count
from stdlib import math
from cStringIO import StringIO
import sys

def exercise(verbose=0):
  distance_ideal = 1.8
  default_vdw_distance = 3.6
  vdw_1_4_factor = 3.5/3.6
  sites_cart_manual = flex.vec3_double([
    (1,3,0), (2,3,0), (3,2,0), (3,1,0), (4,1,0), (3,4,0), (4,3,0), (5,3,0),
    (6,2,0), (7,2,0), (8,3,0), (7,4,0), (6,4,0), (7,5,0), (6,6,0), (8,6,0)])
  bond_proxies = geometry_restraints.bond_sorted_asu_proxies(asu_mappings=None)
  for i_seqs in [(0,1),(1,2),(2,3),(3,4),(1,5),(2,6),(5,6),
                 (6,7),(7,8),(8,9),(9,10),(10,11),(11,12),
                 (12,7),(11,13),(13,14),(14,15),(15,13)]:
    bond_proxies.process(geometry_restraints.bond_simple_proxy(
      i_seqs=i_seqs, distance_ideal=distance_ideal, weight=100))
  angle_proxies = geometry_restraints.shared_angle_proxy()
  for i_seqs,angle_ideal in [[(0,1,2),135],
                             [(0,1,5),135],
                             [(1,2,3),135],
                             [(3,2,6),135],
                             [(2,3,4),120],
                             [(1,2,6),90],
                             [(2,6,5),90],
                             [(6,5,1),90],
                             [(5,1,2),90],
                             [(2,6,7),135],
                             [(5,6,7),135],
                             [(6,7,8),120],
                             [(6,7,12),120],
                             [(7,8,9),120],
                             [(8,9,10),120],
                             [(9,10,11),120],
                             [(10,11,12),120],
                             [(11,12,7),120],
                             [(12,7,8),120],
                             [(10,11,13),120],
                             [(12,11,13),120],
                             [(11,13,15),150],
                             [(11,13,14),150],
                             [(13,15,14),60],
                             [(15,14,13),60],
                             [(14,13,15),60]]:
    angle_proxies.append(geometry_restraints.angle_proxy(
      i_seqs=i_seqs, angle_ideal=angle_ideal, weight=1))
  if (0 or verbose):
    f = open("manual.pdb", "w")
    for serial,site in zip(count(1), sites_cart_manual):
      print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
    print >> f, "END"
    f.close()
  sites_cart = sites_cart_manual.deep_copy()
  assert bond_proxies.asu.size() == 0
  bond_params_table = geometry_restraints.extract_bond_params(
    n_seq=sites_cart.size(),
    bond_simple_proxies=bond_proxies.simple)
  manager = geometry_restraints.manager.manager(
    bond_params_table=bond_params_table,
    angle_proxies=angle_proxies)
  minimized = geometry_restraints.lbfgs.lbfgs(
    sites_cart=sites_cart,
    geometry_restraints_manager=manager,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations=1000))
  sites_cart_minimized_1 = sites_cart.deep_copy()
  s = StringIO()
  manager.show_interactions(f=s)
  if (0 or verbose):
    f = open("minimized_1.pdb", "w")
    for serial,site in zip(count(1), sites_cart_minimized_1):
      print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
    print >> f, "END"
    f.close()
  bond_deltas = geometry_restraints.bond_deltas(
    sites_cart=sites_cart_minimized_1,
    proxies=bond_proxies.simple)
  angle_deltas = geometry_restraints.angle_deltas(
    sites_cart=sites_cart_minimized_1,
    proxies=angle_proxies)
  if (0 or verbose):
    for proxy,delta in zip(bond_proxies.simple, bond_deltas):
      print "bond:", proxy.i_seqs, delta
    for proxy,delta in zip(angle_proxies, angle_deltas):
      print "angle:", proxy.i_seqs, delta
  assert flex.max(flex.abs(bond_deltas)) < 1.e-6
  assert flex.max(flex.abs(angle_deltas)) < 1.e-6
  sites_cart += matrix.col((1,1,0)) - matrix.col(sites_cart.min())
  unit_cell_lengths = list(  matrix.col(sites_cart.max())
                           + matrix.col((1,-1.2,4)))
  unit_cell_lengths[1] *= 2
  xray_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=unit_cell_lengths,
      space_group_symbol="P112"))
  for serial,site in zip(count(1), sites_cart):
    xray_structure.add_scatterer(xray.scatterer(
      label="C%02d"%serial,
      site=xray_structure.unit_cell().fractionalize(site)))
  if (0 or verbose):
    xray_structure.show_summary().show_scatterers()
  p1_structure = (xray_structure
    .apply_shift((-.5,-.5,0))
    .expand_to_p1()
    .apply_shift((.5,.5,0)))
  for shift in [(1,0,0), (0,1,0), (0,0,1)]:
    p1_structure.add_scatterers(p1_structure.apply_shift(shift).scatterers())
  if (0 or verbose):
    open("p1_structure.pdb", "w").write(p1_structure.as_pdb_file())
  nonbonded_cutoff = 7
  asu_mappings = xray_structure.asu_mappings(
    buffer_thickness=nonbonded_cutoff)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  geometry_restraints.add_pairs(bond_asu_table, bond_proxies.simple)
  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
    pair_asu_table=bond_asu_table,
    max_shell=3)
  shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
    for shell_asu_table in shell_asu_tables]
  bond_params_table = geometry_restraints.extract_bond_params(
    n_seq=sites_cart.size(),
    bond_simple_proxies=bond_proxies.simple)
  atom_energy_types = flex.std_string(sites_cart.size(), "Default")
  nonbonded_params = geometry_restraints.nonbonded_params(
    factor_1_4_interactions=vdw_1_4_factor,
    const_shrink_1_4_interactions=0,
    default_distance=default_vdw_distance)
  nonbonded_params.distance_table.setdefault(
    "Default")["Default"] = default_vdw_distance
  pair_proxies = geometry_restraints.pair_proxies(
    nonbonded_params = nonbonded_params,
    nonbonded_types=atom_energy_types,
    bond_params_table=bond_params_table,
    shell_asu_tables=shell_asu_tables,
    bonded_distance_cutoff=0,
    nonbonded_distance_cutoff=nonbonded_cutoff,
    nonbonded_buffer=0)
  if (0 or verbose):
    print "pair_proxies.bond_proxies.n_total():", \
           pair_proxies.bond_proxies.n_total(),
    print "simple:", pair_proxies.bond_proxies.simple.size(),
    print "sym:", pair_proxies.bond_proxies.asu.size()
    print "pair_proxies.nonbonded_proxies.n_total():", \
           pair_proxies.nonbonded_proxies.n_total(),
    print "simple:", pair_proxies.nonbonded_proxies.simple.size(),
    print "sym:", pair_proxies.nonbonded_proxies.asu.size()
    print "pair_proxies.n_nonbonded:", pair_proxies.n_nonbonded
    print "pair_proxies.n_1_4:      ", pair_proxies.n_1_4
    print "min_distance_nonbonded: %.2f" % flex.min(
      geometry_restraints.nonbonded_deltas(
        sites_cart=sites_cart,
        sorted_asu_proxies=pair_proxies.nonbonded_proxies,
        function=geometry_restraints.prolsq_repulsion_function()))
  vdw_1_sticks = []
  vdw_2_sticks = []
  for proxy in pair_proxies.nonbonded_proxies.simple:
    if (proxy.vdw_distance == default_vdw_distance):
      vdw_1_sticks.append(pml_stick(
        begin=sites_cart[proxy.i_seqs[0]],
        end=sites_cart[proxy.i_seqs[1]]))
    else:
      vdw_2_sticks.append(pml_stick(
        begin=sites_cart[proxy.i_seqs[0]],
        end=sites_cart[proxy.i_seqs[1]]))
  mps = asu_mappings.mappings()
  for proxy in pair_proxies.nonbonded_proxies.asu:
    if (proxy.vdw_distance == default_vdw_distance):
      vdw_1_sticks.append(pml_stick(
        begin=mps[proxy.i_seq][0].mapped_site(),
        end=mps[proxy.j_seq][proxy.j_sym].mapped_site()))
    else:
      vdw_2_sticks.append(pml_stick(
        begin=mps[proxy.i_seq][0].mapped_site(),
        end=mps[proxy.j_seq][proxy.j_sym].mapped_site()))
  if (0 or verbose):
    pml_write(f=open("vdw_1.pml", "w"), label="vdw_1", sticks=vdw_1_sticks)
    pml_write(f=open("vdw_2.pml", "w"), label="vdw_2", sticks=vdw_2_sticks)
  #
  i_pdb = count(2)
  for use_crystal_symmetry in [False, True]:
    if (not use_crystal_symmetry):
      crystal_symmetry = None
      site_symmetry_table = None
    else:
      crystal_symmetry = xray_structure
      site_symmetry_table = xray_structure.site_symmetry_table()
    for sites_cart in [sites_cart_manual.deep_copy(),
                       sites_cart_minimized_1.deep_copy()]:
      manager = geometry_restraints.manager.manager(
        crystal_symmetry=crystal_symmetry,
        site_symmetry_table=site_symmetry_table,
        nonbonded_params=nonbonded_params,
        nonbonded_types=atom_energy_types,
        nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
        bond_params_table=bond_params_table,
        shell_sym_tables=shell_sym_tables,
        nonbonded_distance_cutoff=nonbonded_cutoff,
        nonbonded_buffer=1,
        angle_proxies=angle_proxies)
      pair_proxies = manager.pair_proxies(sites_cart=sites_cart)
      if (0 or verbose):
        print "len(vdw_1):", pair_proxies.n_nonbonded
        print "len(vdw_2):", pair_proxies.n_1_4
      minimized = geometry_restraints.lbfgs.lbfgs(
        sites_cart=sites_cart,
        geometry_restraints_manager=manager,
        lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
          max_iterations=1000))
      s = StringIO()
      manager.show_interactions(f=s)
      if (0 or verbose):
        minimized.final_target_result.show()
        print "number of function evaluations:", minimized.minimizer.nfun()
        print "n_updates_pair_proxies:", manager.n_updates_pair_proxies
      if (not use_crystal_symmetry):
        assert minimized.final_target_result.bond_residual_sum < 1.e-3
        assert minimized.final_target_result.nonbonded_residual_sum < 0.1
      else:
        assert minimized.final_target_result.bond_residual_sum < 1.e-2
        assert minimized.final_target_result.nonbonded_residual_sum < 0.1
      assert minimized.final_target_result.angle_residual_sum < 1.e-3
      if (0 or verbose):
        pdb_file_name = "minimized_%d.pdb" % i_pdb.next()
        print "Writing file:", pdb_file_name
        f = open(pdb_file_name, "w")
        for serial,site in zip(count(1), sites_cart):
          print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
        print >> f, "END"
        f.close()
  print "OK"

if (__name__ == "__main__"):
  exercise(verbose=("--verbose" in sys.argv[1:]))
