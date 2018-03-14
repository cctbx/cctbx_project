from __future__ import division
import mmtbx.refinement.tardy
from mmtbx.monomer_library import pdb_interpretation
import iotbx.pdb
import iotbx.phil
import cctbx.geometry_restraints
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
from libtbx.utils import Sorry, format_cpu_times
from libtbx import Auto, group_args
import random
import sys, os

class potential_object(object):

  allowed_origin_shifts_need_to_be_suppressed = False

  def __init__(O,
        density_map,
        geo_manager,
        reduced_geo_manager,
        prolsq_repulsion_function_changes,
        real_space_gradients_delta,
        real_space_target_weight,
        ideal_sites_cart,
        site_labels=None,
        orca_experiments=False):
    O.density_map = density_map
    O.geo_manager = geo_manager
    O.reduced_geo_manager = reduced_geo_manager
    assert isinstance(
      geo_manager.nonbonded_function,
      cctbx.geometry_restraints.prolsq_repulsion_function)
    c = prolsq_repulsion_function_changes
    O.custom_nonbonded_function = geo_manager.nonbonded_function \
      .customized_copy(
        c_rep=c.c_rep,
        k_rep=c.k_rep,
        irexp=c.irexp,
        rexp=c.rexp)
    O.real_space_gradients_delta = real_space_gradients_delta
    O.real_space_target_weight = real_space_target_weight
    O.ideal_sites_cart = ideal_sites_cart
    O.site_labels = site_labels
    O.orca_experiments = orca_experiments
    #
    O.last_sites_moved = None
    O.f = None
    O.g = None
    O.last_grms = None

  def e_pot(O, sites_moved):
    if (   O.last_sites_moved is None
        or O.last_sites_moved.id() is not sites_moved.id()):
      O.last_sites_moved = sites_moved
      sites_cart = sites_moved
      #
      if (O.reduced_geo_manager is None):
        flags = None
        if (O.orca_experiments):
          flags = cctbx.geometry_restraints.flags.flags(
            nonbonded=False, default=True)
      else:
        # computing nonbonded interactions only with geo_manager,
        # contributions from other restraints below with reduced_geo_manager
        flags = cctbx.geometry_restraints.flags.flags(
          nonbonded=True, default=False)
      geo_energies = O.geo_manager.energies_sites(
        sites_cart=sites_cart,
        flags=flags,
        custom_nonbonded_function=O.custom_nonbonded_function,
        compute_gradients=True)
      if (0): # XXX orca_experiments
        print "geo_energies:"
        geo_energies.show()
      if (0): # XXX orca_experiments
        O.geo_manager.show_sorted(site_labels=O.site_labels)
      O.f = geo_energies.target
      O.g = geo_energies.gradients
      if (O.reduced_geo_manager is not None):
        reduced_geo_energies = O.reduced_geo_manager.energies_sites(
          sites_cart=sites_cart,
          compute_gradients=True)
        O.f += reduced_geo_energies.target
        O.g += reduced_geo_energies.gradients
      O.last_grms = group_args(geo=flex.mean_sq(O.g.as_double())**0.5)
      #
      if (O.density_map is not None):
        rs_f = maptbx.real_space_target_simple(
          unit_cell=O.geo_manager.crystal_symmetry.unit_cell(),
          density_map=O.density_map,
          sites_cart=sites_cart,
          selection=flex.bool(sites_cart.size(),True))
        rs_g = maptbx.real_space_gradients_simple(
          unit_cell=O.geo_manager.crystal_symmetry.unit_cell(),
          density_map=O.density_map,
          sites_cart=sites_cart,
          delta=O.real_space_gradients_delta,
          selection=flex.bool(sites_cart.size(),True))
        rs_f *= -O.real_space_target_weight
        rs_g *= -O.real_space_target_weight
        O.f += rs_f
        O.g += rs_g
        O.last_grms.real = flex.mean_sq(rs_g.as_double())**0.5
        O.last_grms.real_or_xray = "real"
      #
      O.last_grms.total = flex.mean_sq(O.g.as_double())**0.5
    return O.f

  def d_e_pot_d_sites(O, sites_moved):
    O.e_pot(sites_moved=sites_moved)
    return O.g

def cartesian_random_displacements(sites_cart, target_rmsd, max_trials=10):
  assert sites_cart.size() != 0
  for i in xrange(max_trials):
    shifts_cart = flex.vec3_double(flex.random_double(
      size=sites_cart.size()*3, factor=2) - 1)
    rmsd = (flex.sum(shifts_cart.dot()) / sites_cart.size()) ** 0.5
    if (rmsd > 1.e-6): break # to avoid numerical problems
  else:
    raise RuntimeError
  shifts_cart *= (target_rmsd / rmsd)
  return sites_cart + shifts_cart

def run_test(params, pdb_files, other_files, callback=None, log=None):
  if (log is None): log = sys.stdout
  if (params.random_seed is not None):
    random.seed(params.random_seed)
    flex.set_random_seed(value=params.random_seed)
  #
  if (len(pdb_files) != 0):
    print >> log, "PDB files:"
    for file_name in pdb_files:
      print >> log, " ", file_name
    print >> log
  if (len(other_files) != 0):
    print >> log, "Other files:"
    for file_name in other_files:
      print >> log, " ", file_name
    print >> log
  #
  assert len(pdb_files) in [1, 2]
  #
  pdb_interpretation_params = pdb_interpretation.master_params.extract()
  pdb_interpretation_params.dihedral_function_type \
    = params.dihedral_function_type
  processed_pdb_files = pdb_interpretation.run(
    args=pdb_files[-1:]+other_files,
    params=pdb_interpretation_params,
    strict_conflict_handling=False,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    return_all_processed_pdb_files=True,
    log=log)
  assert len(processed_pdb_files) == 1
  print >> log
  #
  xs = processed_pdb_files[0].xray_structure()
  geo_manager = processed_pdb_files[0].geometry_restraints_manager()
  labels = [sc.label for sc in xs.scatterers()]
  ideal_sites_cart = xs.sites_cart()
  sites = ideal_sites_cart
  masses = xs.atomic_weights()
  tardy_tree_simple_connectivity = geo_manager.construct_tardy_tree(sites=sites)
  rmsd_calculator = tardy_tree_simple_connectivity.rmsd_calculator()
  #
  if (params.tardy_displacements is not None):
    def get_tardy_model_no_potential():
      return scitbx.rigid_body.tardy_model(
        labels=labels,
        sites=sites,
        masses=masses,
        tardy_tree=tardy_tree_simple_connectivity,
        potential_obj=None)
    def get_tardy_model_no_density():
      tardy_tree = scitbx.graph.tardy_tree.construct(
        n_vertices=len(sites), edge_list=[])
      tardy_tree.build_tree()
      potential_obj = potential_object(
        density_map=None,
        geo_manager=geo_manager,
        reduced_geo_manager=None,
        prolsq_repulsion_function_changes=
          params.prolsq_repulsion_function_changes,
        real_space_gradients_delta=None,
        real_space_target_weight=None,
        ideal_sites_cart=None)
      return scitbx.rigid_body.tardy_model(
        labels=labels,
        sites=sites,
        masses=masses,
        tardy_tree=tardy_tree,
        potential_obj=potential_obj)
    if (params.tardy_displacements is Auto):
      auto_params = params.tardy_displacements_auto
      target_rmsd = \
          params.structure_factors_high_resolution \
        * auto_params.rmsd_vs_high_resolution_factor
      target_rmsd_tol = \
          params.structure_factors_high_resolution \
        * auto_params.rmsd_tolerance
      assert target_rmsd > 0
      assert target_rmsd_tol > 0
      print >> log, "Random displacements (%s):" \
        % params.tardy_displacements_auto.parameterization
      print >> log, "  high resolution: %.6g" \
        % params.structure_factors_high_resolution
      print >> log, "  target rmsd: %.6g" % target_rmsd
      print >> log, "  target rmsd tolerance: %.6g" % target_rmsd_tol
      log.flush()
      def raise_max_steps_exceeded(var_name, rmsd_history):
        msg = [
          "tardy_displacements_auto.max_steps exceeded:",
          "  %        -13s  rmsd" % var_name]
        for var_rmsd in rmsd_history:
          msg.append("  %13.6e  %13.6e" % var_rmsd)
        raise Sorry("\n".join(msg))
      if (params.tardy_displacements_auto.parameterization == "cartesian"):
        multiplier = 1.5
        rmsd_history = []
        for i_step in xrange(auto_params.max_steps):
          sites = cartesian_random_displacements(
            sites_cart=ideal_sites_cart,
            target_rmsd=target_rmsd*multiplier)
          tardy_model = get_tardy_model_no_density()
          tardy_model.minimization(max_iterations=20)
          sites = tardy_model.sites_moved()
          sites_moved = sites
          rmsd = rmsd_calculator(sites_moved, ideal_sites_cart)
          rmsd_history.append((multiplier, rmsd))
          print >> log, "    multiplier, rmsd: %13.6e, %13.6e" \
            % rmsd_history[-1]
          log.flush()
          if (rmsd < target_rmsd - target_rmsd_tol):
            if (rmsd != 0):
              multiplier = min(
                multiplier*2, max(
                  multiplier*1.2,
                    target_rmsd / rmsd))
            else:
              multiplier *= 2
          else:
            if (rmsd <= target_rmsd + target_rmsd_tol):
              tardy_model.minimization(max_iterations=500)
              sites = tardy_model.sites_moved()
              sites_moved = sites
              rmsd = rmsd_calculator(sites_moved, ideal_sites_cart)
              rmsd_history.append((0, rmsd))
              print >> log, "    multiplier, rmsd: %13.6e, %13.6e" \
                % rmsd_history[-1]
              log.flush()
              break
            multiplier *= max(0.5, target_rmsd/rmsd)
        else:
          raise_max_steps_exceeded(
            var_name="multiplier", rmsd_history=rmsd_history)
        del rmsd_history
        print >> log, "  actual rmsd: %.6g" % rmsd
        print >> log
      elif (params.tardy_displacements_auto.parameterization == "constrained"):
        tardy_model = get_tardy_model_no_potential()
        tardy_model.assign_random_velocities()
        delta_t = auto_params.first_delta_t
        rmsd_history = []
        assert auto_params.max_steps > 0
        for i_step in xrange(auto_params.max_steps):
          prev_q = tardy_model.pack_q()
          prev_qd = tardy_model.pack_qd()
          tardy_model.dynamics_step(delta_t=delta_t)
          sites_moved = tardy_model.sites_moved()
          rmsd = rmsd_calculator(sites_moved, ideal_sites_cart)
          rmsd_history.append((delta_t, rmsd))
          if (rmsd < target_rmsd - target_rmsd_tol):
            delta_t *= 2 - rmsd / target_rmsd
          else:
            if (rmsd <= target_rmsd + target_rmsd_tol):
              break
            tardy_model.unpack_q(q_packed=prev_q)
            tardy_model.unpack_qd(qd_packed=prev_qd)
            delta_t *= 0.5
          prev_q = None
          prev_qd = None
        else:
          raise_max_steps_exceeded(
            var_name="delta_t", rmsd_history=rmsd_history)
        del rmsd_history
        print >> log, "  actual rmsd: %.6g" % rmsd
        print >> log, "  tardy_displacements=%s" % ",".join(
          ["%.6g" % v for v in tardy_model.pack_q()])
        print >> log
        sites = tardy_model.sites_moved()
      else:
        raise AssertionError
    else:
      tardy_model = get_tardy_model_no_potential()
      q = tardy_model.pack_q()
      if (len(params.tardy_displacements) != len(q)):
        print >> log, "tardy_displacements:", params.tardy_displacements
        hinge_edges = tardy_model.tardy_tree.cluster_manager.hinge_edges
        assert len(hinge_edges) == tardy_model.bodies_size()
        dofej = tardy_model.degrees_of_freedom_each_joint()
        qsej = tardy_model.q_size_each_joint()
        for ib,(i,j) in enumerate(hinge_edges):
          if (i == -1): si = "root"
          else: si = tardy_model.labels[i]
          sj = tardy_model.labels[j]
          print >> log, "%21s - %-21s: %d dof, %d q_size" % (
            si, sj, dofej[ib], qsej[ib])
        print >> log, "Zero displacements:"
        print >> log, "  tardy_displacements=%s" % ",".join(
          [str(v) for v in q])
        raise Sorry("Incompatible tardy_displacements.")
      tardy_model.unpack_q(q_packed=flex.double(params.tardy_displacements))
      sites = tardy_model.sites_moved()
  #
  if (params.emulate_cartesian):
    tardy_tree = scitbx.graph.tardy_tree.construct(
      n_vertices=len(sites), edge_list=[])
    tardy_tree.build_tree()
  else:
    tardy_tree = tardy_tree_simple_connectivity
  print >> log, "tardy_tree summary:"
  tardy_tree.show_summary(vertex_labels=labels, out=log, prefix="  ")
  print >> log
  #
  if (len(pdb_files) == 2):
    ideal_pdb_inp = iotbx.pdb.input(file_name=pdb_files[0])
    ideal_pdb_hierarchy = ideal_pdb_inp.construct_hierarchy()
    assert ideal_pdb_hierarchy.is_similar_hierarchy(
      processed_pdb_files[0].all_chain_proxies.pdb_hierarchy)
    ideal_sites_cart = ideal_pdb_hierarchy.atoms().extract_xyz()
    xs.set_sites_cart(sites_cart=ideal_sites_cart)
  fft_map = xs.structure_factors(
    d_min=params.structure_factors_high_resolution).f_calc().fft_map()
  fft_map.apply_sigma_scaling()
  #
  assert not params.orca_experiments or not params.emulate_cartesian
  if (params.orca_experiments):
    from mmtbx.refinement import orca
    x = orca.expand(
      labels=labels,
      sites_cart=sites,
      masses=masses,
      geo_manager=geo_manager)
    labels = x.labels
    sites = x.sites_cart
    masses = x.masses
    geo_manager = x.geo_manager
    tardy_tree = x.tardy_tree
    x_ideal_sites_cart = flex.vec3_double()
    for i_orc,i_seq in x.indices:
      x_ideal_sites_cart.append(ideal_sites_cart[i_seq])
    ideal_sites_cart = x_ideal_sites_cart
    rmsd_calculator = x.rmsd_calculator(
      tardy_tree_rmsd_calculator=rmsd_calculator)
  #
  if (params.emulate_cartesian or params.keep_all_restraints):
    reduced_geo_manager = None
  else:
    reduced_geo_manager = geo_manager.reduce_for_tardy(tardy_tree=tardy_tree)
  real_space_gradients_delta = \
      params.structure_factors_high_resolution \
    * params.real_space_gradients_delta_resolution_factor
  potential_obj = potential_object(
    density_map=fft_map.real_map(),
    geo_manager=geo_manager,
    reduced_geo_manager=reduced_geo_manager,
    prolsq_repulsion_function_changes=params.prolsq_repulsion_function_changes,
    real_space_gradients_delta=real_space_gradients_delta,
    real_space_target_weight=params.real_space_target_weight,
    ideal_sites_cart=ideal_sites_cart,
    site_labels=labels,
    orca_experiments=params.orca_experiments)
  tardy_model = scitbx.rigid_body.tardy_model(
    labels=labels,
    sites=sites,
    masses=masses,
    tardy_tree=tardy_tree,
    potential_obj=potential_obj)
  mmtbx.refinement.tardy.action(
    tardy_model=tardy_model,
    params=params,
    rmsd_calculator=rmsd_calculator,
    callback=callback,
    log=log)
  print >> log

def get_master_phil():
  return iotbx.phil.parse(
    input_string=mmtbx.refinement.tardy.master_phil_str + """\
structure_factors_high_resolution = 3
  .type = float
real_space_target_weight = 100
  .type = float
real_space_gradients_delta_resolution_factor = 1/3
  .type = float
orca_experiments = False
  .type = bool
keep_all_restraints = False
  .type = bool
random_seed = None
  .type = int
tardy_displacements = None
  .type = floats
tardy_displacements_auto {
  parameterization = *constrained cartesian
    .type = choice
    .optional = False
  rmsd_vs_high_resolution_factor = 1/3
    .type = float
  rmsd_tolerance = 0.1
    .type = float
  first_delta_t = 0.001
    .type = float
  max_steps = 100
    .type = int
}
%(dihedral_function_type_params_str)s
""" % pdb_interpretation.__dict__)

def run(args, callback=None):
  master_phil = get_master_phil()
  argument_interpreter = master_phil.command_line_argument_interpreter()
  master_phil_str_overrides = """\
start_temperature_kelvin = 2500
final_temperature_kelvin = 300
number_of_cooling_steps = 44
number_of_time_steps = 1
time_step_pico_seconds = 0.001
minimization_max_iterations = None
"""
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  pdb_files = []
  other_files = []
  for arg in args:
    if (os.path.isfile(arg)):
      if (iotbx.pdb.is_pdb_file(file_name=arg)):
        pdb_files.append(arg)
        arg = None
      else:
        try: phil_obj = iotbx.phil.parse(file_name=arg)
        except KeyboardInterrupt: raise
        except Exception:
          other_files.append(arg)
        else:
          phil_objects.append(phil_obj)
        arg = None
    if (arg is not None):
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Command-line argument not recognized: %s" % arg)
      else: phil_objects.append(command_line_params)
  params = master_phil.fetch(sources=phil_objects).extract()
  master_phil.format(params).show()
  print
  run_test(
    params=params,
    pdb_files=pdb_files,
    other_files=other_files,
    callback=callback)
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
