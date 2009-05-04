from mmtbx.dynamics import \
  kinetic_energy_as_temperature, \
  temperature_as_kinetic_energy
from mmtbx.dynamics.constants import akma_time_as_pico_seconds
from cctbx import xray
import cctbx.geometry_restraints
from cctbx.array_family import flex
from scitbx.rigid_body.essence import tst_molecules
from scitbx.graph import tardy_tree
from scitbx import matrix
from libtbx.utils import Sorry
from libtbx.str_utils import format_value
from libtbx import group_args

master_phil_str = """\
  start_temperature_kelvin = 5000
    .type = float
  final_temperature_kelvin = 300
    .type = float
  velocity_scaling = True
    .type = bool
  temperature_cap_factor = 1.5
    .type = float
  excessive_temperature_factor = 5
    .type = float
  number_of_cooling_steps = 20
    .type = int
  number_of_time_steps = 50
    .type = int
  time_step_pico_seconds = 0.004
    .type = float
  minimization_max_iterations = 25
    .type = int
  nonbonded_attenuation_factor = 0.5
    .type = float
  omit_bonds_with_slack_greater_than = 0
    .type = float
  constrain_dihedrals_with_sigma_less_than = 10
    .type = float
"""

class potential_object(object):

  def __init__(O,
        nonbonded_attenuation_factor,
        fmodels,
        model,
        target_weights,
        reduced_geo_manager):
    O.fmodels = fmodels
    O.model = model
    O.weights = target_weights.xyz_weights_result
    O.reduced_geo_manager = reduced_geo_manager
    if (nonbonded_attenuation_factor is None):
      O.custom_nonbonded_function = None
    else:
      assert nonbonded_attenuation_factor > 0
      assert nonbonded_attenuation_factor <= 1
      nonbonded_function = model.restraints_manager.geometry.nonbonded_function
      assert isinstance(
        nonbonded_function,
        cctbx.geometry_restraints.prolsq_repulsion_function)
      O.custom_nonbonded_function = nonbonded_function.customized_copy(
        k_rep=nonbonded_attenuation_factor)
    O.fmodels.create_target_functors()
    O.fmodels.prepare_target_functors_for_minimization()
    O.last_sites_moved = None
    O.f = None
    O.g = None
    O.last_grms = None

  def e_pot(O, sites_moved):
    if (O.last_sites_moved is not sites_moved):
      O.last_sites_moved = sites_moved
      xs = O.fmodels.fmodel_xray().xray_structure
      assert len(sites_moved) == xs.scatterers().size()
      sites_cart = flex.vec3_double(sites_moved)
      xs.set_sites_cart(sites_cart=sites_cart)
      O.fmodels.update_xray_structure(update_f_calc=True)
      xs.scatterers().flags_set_grads(state=False)
      xs.scatterers().flags_set_grad_site(
        iselection=flex.size_t_range(xs.scatterers().size()))
      tg = O.fmodels.target_and_gradients(
        weights=O.weights,
        compute_gradients=True)
      O.f = tg.target()
      O.g = tg.gradients()
      assert O.g.size() == len(sites_moved) * 3
      reduced_geo_energies = O.reduced_geo_manager.energies_sites(
        sites_cart=sites_cart,
        compute_gradients=True)
      other_energies = O.model.restraints_manager.energies_sites(
        sites_cart=sites_cart,
        geometry_flags=cctbx.geometry_restraints.flags.flags(nonbonded=True),
        custom_nonbonded_function=O.custom_nonbonded_function,
        compute_gradients=True)
      nfw = other_energies.normalization_factor * O.weights.w
      O.f += other_energies.target * O.weights.w
      O.f += reduced_geo_energies.target * nfw
      gg = other_energies.gradients * O.weights.w
      gg += reduced_geo_energies.gradients * nfw
      assert nfw != 0
      scale = 1 / nfw
      O.last_grms = group_args(
        geo=scale*flex.mean_sq(gg.as_double())**0.5,
        xray=scale*flex.mean_sq(O.g)**0.5)
      xray.minimization.add_gradients(
        scatterers=xs.scatterers(),
        xray_gradients=O.g,
        site_gradients=gg)
      O.f *= scale
      O.g *= scale
      O.last_grms.total = flex.mean_sq(O.g)**0.5
    return O.f

  def d_e_pot_d_sites(O, sites_moved):
    O.e_pot(sites_moved=sites_moved)
    return matrix.col_list(flex.vec3_double(O.g))

# XXX move to geometry_restraints_manager.py
def reduce_geo_manager(
      geo_manager,
      cluster_manager,
      omit_slack_greater_than=0):
  from cctbx import crystal
  from cctbx import sgtbx
  from scitbx.graph.utils import construct_edge_sets
  #
  cluster_indices = cluster_manager.cluster_indices
  n_vertices = len(cluster_indices)
  loop_edge_sets = construct_edge_sets(
    n_vertices=n_vertices,
    edge_list=cluster_manager.loop_edges)
  #
  def get():
    result = crystal.pair_sym_table(n_vertices)
    bond_params_table = geo_manager.bond_params_table
    for i,pair_sym_dict in enumerate(geo_manager.shell_sym_tables[0]):
      reduced_pair_sym_dict = result[i]
      for j,sym_ops in pair_sym_dict.items():
        reduced_sym_ops = sgtbx.stl_vector_rt_mx()
        for sym_op in sym_ops:
          if (sym_op.is_unit_mx()):
            if (bond_params_table[i][j].slack > omit_slack_greater_than):
              continue
            if (j not in loop_edge_sets[i]):
              continue
          reduced_sym_ops.append(sym_op)
        if (reduced_sym_ops.size() != 0):
          reduced_pair_sym_dict[j] = reduced_sym_ops
    return [result]
  reduced_shell_sym_tables = get()
  #
  def get():
    result = cctbx.geometry_restraints.shared_angle_proxy()
    for proxy in geo_manager.angle_proxies:
      i,j,k = proxy.i_seqs
      if (   j in loop_edge_sets[i]
          or j in loop_edge_sets[k]):
        result.append(proxy)
    if (result.size() == 0):
      return None
    return result
  reduced_angle_proxies = get()
  #
  def get():
    result = cctbx.geometry_restraints.shared_dihedral_proxy()
    for proxy in geo_manager.dihedral_proxies:
      if (len(set([cluster_indices[i] for i in proxy.i_seqs])) != 1):
        result.append(proxy)
    if (result.size() == 0):
      return None
    return result
  reduced_dihedral_proxies = get()
  #
  return cctbx.geometry_restraints.manager.manager(
    crystal_symmetry=geo_manager.crystal_symmetry,
    site_symmetry_table=geo_manager.site_symmetry_table,
    bond_params_table=geo_manager.bond_params_table,
    shell_sym_tables=reduced_shell_sym_tables,
    angle_proxies=reduced_angle_proxies,
    dihedral_proxies=reduced_dihedral_proxies)

def run(fmodels, model, target_weights, params, log):
  assert fmodels.fmodel_neutron() is None # not implemented
  assert model.ias_selection is None # tardy+ias is not a useful combination
  xs = fmodels.fmodel_xray().xray_structure
  sites_cart_start = xs.sites_cart()
  sites = matrix.col_list(sites_cart_start)
  tt = tardy_tree.construct(
    sites=sites,
    edge_list=model.restraints_manager.geometry.simple_edge_list(
      omit_slack_greater_than=params.omit_bonds_with_slack_greater_than),
    external_clusters=model.restraints_manager.geometry
      .rigid_clusters_due_to_dihedrals_and_planes(
        constrain_dihedrals_with_sigma_less_than
          =params.constrain_dihedrals_with_sigma_less_than))
  tt.finalize()
  for i,j in tt.collinear_bonds_edge_list:
    s = xs.scatterers()
    print >> log, "tardy collinear bond:", s[i].label
    print >> log, "                     ", s[j].label
  log.flush()
  potential_obj = potential_object(
    nonbonded_attenuation_factor=params.nonbonded_attenuation_factor,
    fmodels=fmodels,
    model=model,
    target_weights=target_weights,
    reduced_geo_manager=reduce_geo_manager(
      geo_manager=model.restraints_manager.geometry,
      cluster_manager=tt.cluster_manager))
  sim = tst_molecules.simulation(
    labels=[sc.label for sc in xs.scatterers()],
    sites=sites,
    bonds=tt.edge_list,
    cluster_manager=tt.cluster_manager,
    potential_obj=potential_obj,
    bodies=tst_molecules.construct_bodies(
      sites=sites,
      masses=xs.atomic_weights(),
      cluster_manager=tt.cluster_manager))
  del sites
  qd_e_kin_scales = sim.assign_random_velocities()
  def e_as_t(e):
    return kinetic_energy_as_temperature(dof=sim.degrees_of_freedom, e=e)
  def t_as_e(t):
    return temperature_as_kinetic_energy(dof=sim.degrees_of_freedom, t=t)
  time_step_akma = params.time_step_pico_seconds / akma_time_as_pico_seconds
  print >> log, "tardy dynamics:"
  print >> log, "  kinetic energy sensitivity to generalized velocities:"
  qd_e_kin_scales.min_max_mean().show(out=log, prefix="    ")
  print >> log, "  time step: %7.5f pico seconds" % (
    params.time_step_pico_seconds)
  print >> log, "  velocity scaling:", params.velocity_scaling
  log.flush()
  n_time_steps = 0
  for i_cool_step in xrange(params.number_of_cooling_steps+1):
    t_target = params.start_temperature_kelvin \
             - i_cool_step * (  params.start_temperature_kelvin
                              - params.final_temperature_kelvin) \
                           / params.number_of_cooling_steps
    e_kin_target = t_as_e(t=t_target)
    def reset_e_kin(msg):
      sim.reset_e_kin(e_kin_target=e_kin_target)
      print >> log, "  %s temperature: %8.2f K" % (
        msg, e_as_t(e=sim.e_kin()))
      log.flush()
      return True
    show_column_headings = reset_e_kin("new target")
    for i_time_step in xrange(params.number_of_time_steps):
      assert params.temperature_cap_factor > 1.0
      assert params.excessive_temperature_factor > params.temperature_cap_factor
      if (sim.e_kin() > e_kin_target * params.temperature_cap_factor):
        print >> log, "  system temperature is too high:"
        if (sim.e_kin() > e_kin_target * params.excessive_temperature_factor):
          print >> log, "    excessive_temperature_factor: %.6g" % \
            params.excessive_temperature_factor
          print >> log, "    excessive temperature limit: %.2f K" % (
            t_target * params.excessive_temperature_factor)
          print >> log, "    time_step_pico_seconds: %.6g" % (
            params.time_step_pico_seconds)
          log.flush()
          raise Sorry(
            "Excessive system temperature in torsion angle dynamics:\n"
            "  Please try again with a smaller time_step_pico_seconds.")
        print >> log, "     temperature_cap_factor: %.6g" % \
          params.temperature_cap_factor
        print >> log, "     temperature cap: %.2f K" % (
          t_target * params.temperature_cap_factor)
        show_column_headings = reset_e_kin("resetting")
      e_kin_before, e_tot_before = sim.e_kin(), sim.e_tot()
      sim.dynamics_step(delta_t=time_step_akma)
      e_kin_after, e_tot_after = sim.e_kin(), sim.e_tot()
      if (e_tot_before > e_tot_after):
        fluct_e_tot = -e_tot_after / e_tot_before
      elif (e_tot_after > e_tot_before):
        fluct_e_tot = e_tot_before / e_tot_after
      else:
        fluct_e_tot = None
      if (params.velocity_scaling):
        sim.reset_e_kin(e_kin_target=e_kin_target)
      n_time_steps += 1
      if (show_column_headings):
        show_column_headings = False
        log.write("""\
          coordinate                   fluctuations        gradient rms
    step      rmsd        temp        temp   e_total     geo    xray   total
""")
      print >> log, "    %4d  %8.4f A  %8.2f K  %8.2f K  %s" \
        "  %6.2f  %6.2f  %6.2f" % (
          n_time_steps,
          xs.sites_cart().rms_difference(sites_cart_start),
          e_as_t(e=sim.e_kin()),
          e_as_t(e=e_kin_after-e_kin_before),
          format_value(format="%6.3f", value=fluct_e_tot),
          sim.potential_obj.last_grms.geo,
          sim.potential_obj.last_grms.xray,
          sim.potential_obj.last_grms.total)
      log.flush()
  if (params.minimization_max_iterations > 0):
    print >> log, "tardy gradient-driven minimization:"
    log.flush()
    def show_rms(minimizer=None):
      print >> log, "  coor. rmsd: %8.4f" % (
        xs.sites_cart().rms_difference(sites_cart_start))
      log.flush()
    sim.minimization(
      max_iterations=params.minimization_max_iterations,
      callback_after_step=show_rms)
    print >> log, "After tardy minimization:"
    show_rms()
