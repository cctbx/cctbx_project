from iotbx.pdb.tst_pdb import dump_pdb
from iotbx.pymol import pml_stick, pml_write
from cctbx import geometry_restraints
from cctbx import adp_restraints
import cctbx.geometry_restraints.manager
import cctbx.geometry_restraints.lbfgs
from cctbx import xray
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx import sgtbx
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx import matrix
from libtbx.test_utils import approx_equal, is_below_limit, show_diff
from itertools import count
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
    dump_pdb(file_name="manual.pdb", sites_cart=sites_cart_manual)
  for traditional_convergence_test in [True,False]:
    for sites_cart_selection in [True, False]:
      sites_cart = sites_cart_manual.deep_copy()
      if sites_cart_selection:
        sites_cart_selection = flex.bool(sites_cart.size(), True)
        sites_cart_selection[1] = False
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
          traditional_convergence_test=traditional_convergence_test,
          drop_convergence_test_max_drop_eps=1.e-20,
          drop_convergence_test_iteration_coefficient=1,
          max_iterations=1000),
        sites_cart_selection=sites_cart_selection,
        )
      assert minimized.minimizer.iter() > 100
      sites_cart_minimized_1 = sites_cart.deep_copy()
      s = StringIO()
      manager.show_interactions(f=s)
      if (0 or verbose):
        dump_pdb(
          file_name="minimized_1.pdb", sites_cart=sites_cart_minimized_1)
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
      assert is_below_limit(
        value=flex.max(flex.abs(bond_deltas)), limit=0, eps=1.e-6)
      assert is_below_limit(
        value=flex.max(flex.abs(angle_deltas)), limit=0, eps=2.e-6)
  sites_cart += matrix.col((1,1,0)) - matrix.col(sites_cart.min())
  unit_cell_lengths = list(  matrix.col(sites_cart.max())
                           + matrix.col((1,-1.2,4)))
  unit_cell_lengths[1] *= 2
  unit_cell_lengths[2] *= 2
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
  nonbonded_cutoff = 6.5
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
    bond_params_table=bond_params_table,
    shell_asu_tables=shell_asu_tables,
    model_indices=None,
    conformer_indices=None,
    nonbonded_params=nonbonded_params,
    nonbonded_types=atom_energy_types,
    nonbonded_distance_cutoff_plus_buffer=nonbonded_cutoff)
  if (0 or verbose):
    print "pair_proxies.bond_proxies.n_total():", \
           pair_proxies.bond_proxies.n_total(),
    print "simple:", pair_proxies.bond_proxies.simple.size(),
    print "sym:", pair_proxies.bond_proxies.asu.size()
    print "pair_proxies.nonbonded_proxies.n_total():", \
           pair_proxies.nonbonded_proxies.n_total(),
    print "simple:", pair_proxies.nonbonded_proxies.simple.size(),
    print "sym:", pair_proxies.nonbonded_proxies.asu.size()
    print "min_distance_nonbonded: %.2f" % flex.min(
      geometry_restraints.nonbonded_deltas(
        sites_cart=sites_cart,
        sorted_asu_proxies=pair_proxies.nonbonded_proxies))
  s = StringIO()
  pair_proxies.bond_proxies.show_histogram_of_model_distances(
    sites_cart=sites_cart,
    f=s,
    prefix="[]")
  assert s.getvalue().splitlines()[0] == "[]Histogram of bond lengths:"
  assert s.getvalue().splitlines()[5].startswith("[]      1.80 -     1.80:")
  s = StringIO()
  pair_proxies.bond_proxies.show_histogram_of_deltas(
    sites_cart=sites_cart,
    f=s,
    prefix="][")
  assert s.getvalue().splitlines()[0] == "][Histogram of bond deltas:"
  assert s.getvalue().splitlines()[5].startswith("][     0.000 -    0.000:")
  s = StringIO()
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    max_items=3,
    f=s,
    prefix=":;")
  l = s.getvalue().splitlines()
  assert l[0] == ":;Bond restraints: 18"
  assert l[1] == ":;Sorted by residual:"
  assert l[2].startswith(":;bond ")
  assert l[3].startswith(":;     ")
  assert l[4] == ":;  ideal  model  delta    sigma   weight residual"
  for i in [5,-2]:
    assert l[i].startswith(":;  1.800  1.800 ")
  assert l[-1] == ":;... (remaining 15 not shown)"
  s = StringIO()
  pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
    sites_cart=sites_cart,
    f=s,
    prefix="]^")
  assert not show_diff(s.getvalue(), """\
]^Histogram of nonbonded interaction distances:
]^      2.16 -     3.03: 3
]^      3.03 -     3.89: 12
]^      3.89 -     4.75: 28
]^      4.75 -     5.61: 44
]^      5.61 -     6.48: 54
""")
  s = StringIO()
  pair_proxies.nonbonded_proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    max_items=7,
    f=s,
    prefix=">,")
  assert not show_diff(s.getvalue(), """\
>,Nonbonded interactions: 141
>,Sorted by model distance:
>,nonbonded 15
>,          15
>,   model   vdw sym.op.
>,   2.164 3.600 -x+2,-y+1,z
...
>,nonbonded 4
>,          8
>,   model   vdw
>,   3.414 3.600
>,... (remaining 134 not shown)
""",
    selections=[range(6), range(-5,0)])
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
        angle_proxies=angle_proxies,
        plain_pairs_radius=5)
      manager = manager.select(selection=flex.bool(sites_cart.size(), True))
      manager = manager.select(
        iselection=flex.size_t_range(stop=sites_cart.size()))
      pair_proxies = manager.pair_proxies(sites_cart=sites_cart)
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
        dump_pdb(file_name=pdb_file_name, sites_cart=sites_cart)
      if (manager.site_symmetry_table is None):
        additional_site_symmetry_table = None
      else:
        additional_site_symmetry_table = sgtbx.site_symmetry_table()
      assert manager.new_including_isolated_sites(
        n_additional_sites=0,
        site_symmetry_table=additional_site_symmetry_table,
        nonbonded_types=flex.std_string()).plain_pairs_radius \
          == manager.plain_pairs_radius
      if (crystal_symmetry is not None):
        assert len(manager.plain_pair_sym_table) == 16
        if (0 or verbose):
          manager.plain_pair_sym_table.show()
  #
  xray_structure.set_u_iso(values=flex.double([
    0.77599982480241358, 0.38745781137212021, 0.20667558236418682,
    0.99759840171302094, 0.8917287406687805, 0.64780251325379845,
    0.24878590382983534, 0.59480621182194615, 0.58695637792905142,
    0.33997130213653637, 0.51258699130743735, 0.79760289141276675,
    0.39996577657875021, 0.4329328819341467, 0.70422156561726479,
    0.87260110626999332]))
  class parameters: pass
  parameters.sphere_radius = 5
  parameters.distance_power = 0.7
  parameters.average_power = 0.9
  parameters.wilson_b_weight = 1.3952
  parameters.wilson_b_weight_auto = False
  adp_energies = adp_restraints.energies_iso(
    geometry_restraints_manager=manager,
    xray_structure=xray_structure,
    parameters=parameters,
    wilson_b=None,
    use_hd=False,
    use_u_local_only = False,
    compute_gradients=False,
    gradients=None,
    normalization=False,
    collect=True)
  assert adp_energies.number_of_restraints == 69
  assert approx_equal(adp_energies.residual_sum, 6.24865382467)
  assert adp_energies.gradients is None
  assert adp_energies.u_i.size() == adp_energies.number_of_restraints
  assert adp_energies.u_j.size() == adp_energies.number_of_restraints
  assert adp_energies.r_ij.size() == adp_energies.number_of_restraints
  for wilson_b in [None, 10, 100]:
    finite_difference_gradients = flex.double()
    eps = 1.e-6
    for i_scatterer in xrange(xray_structure.scatterers().size()):
      rs = []
      for signed_eps in [eps, -eps]:
        xray_structure_eps = xray_structure.deep_copy_scatterers()
        xray_structure_eps.scatterers()[i_scatterer].u_iso += signed_eps
        adp_energies = adp_restraints.energies_iso(
          geometry_restraints_manager=manager,
          xray_structure=xray_structure_eps,
          parameters=parameters,
          wilson_b=wilson_b,
          use_u_local_only = False,
          use_hd=False,
          compute_gradients=True,
          gradients=None,
          normalization=False,
          collect=False)
        rs.append(adp_energies.residual_sum)
        assert adp_energies.gradients.size() \
            == xray_structure.scatterers().size()
        assert adp_energies.u_i == None
        assert adp_energies.u_j == None
        assert adp_energies.r_ij == None
      finite_difference_gradients.append((rs[0]-rs[1])/(2*eps))
    sel = flex.bool(xray_structure.scatterers().size(), True)
    xray_structure.scatterers().flags_set_grad_u_iso(sel.iselection())
    adp_energies = adp_restraints.energies_iso(
      geometry_restraints_manager=manager,
      xray_structure=xray_structure,
      parameters=parameters,
      wilson_b=wilson_b,
      use_u_local_only = False,
      use_hd=False,
      compute_gradients=True,
      gradients=None,
      normalization=False,
      collect=False)
    assert approx_equal(adp_energies.gradients, finite_difference_gradients)
  print "OK"

if (__name__ == "__main__"):
  exercise(verbose=("--verbose" in sys.argv[1:]))
