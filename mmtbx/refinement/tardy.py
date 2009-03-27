from mmtbx.dynamics.constants import \
  akma_time_as_pico_seconds, \
  boltzmann_constant_akma
from cctbx import xray
import cctbx.geometry_restraints
from cctbx.array_family import flex
from scitbx.rigid_body.essence import tst_molecules
from scitbx.graph import tardy_tree
from scitbx import matrix

master_phil_str = """\
  number_of_time_steps = 10
    .type = int
  time_step_pico_seconds = 0.004
    .type = float
  temperature_cap_kelvin = 300
    .type = float
  minimization_max_iterations = 10
    .type = int
  nonbonded_attenuation_factor = 0.5
    .type = float
  omit_bonds_with_slack_greater_than = 0
    .type = float
"""

class potential_object(object):

  def __init__(O,
        nonbonded_attenuation_factor,
        fmodels,
        model,
        target_weights):
    O.fmodels = fmodels
    O.model = model
    O.weights = target_weights.xyz_weights_result
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
    O.e_pot_factor = None

  def e_pot_and_normalization_factor(O, sites_moved):
    if (O.last_sites_moved is not sites_moved):
      O.last_sites_moved = sites_moved
      xs = O.fmodels.fmodel_xray().xray_structure
      assert len(sites_moved) == xs.scatterers().size()
      xs.set_sites_cart(sites_cart=flex.vec3_double(sites_moved))
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
      stereochemistry_residuals = O.model.restraints_manager_energies_sites(
        custom_nonbonded_function=O.custom_nonbonded_function,
        compute_gradients=True)
      O.f += stereochemistry_residuals.target * O.weights.w
      xray.minimization.add_gradients(
        scatterers=xs.scatterers(),
        xray_gradients=O.g,
        site_gradients=stereochemistry_residuals.gradients*O.weights.w)
      O.e_pot_normalization_factor = \
          stereochemistry_residuals.normalization_factor \
        * O.weights.w
    return O.f, O.e_pot_normalization_factor

  def d_e_pot_d_sites(O, sites_moved):
    O.e_pot_and_normalization_factor(sites_moved=sites_moved)
    return matrix.col_list(flex.vec3_double(O.g))

def run(fmodels, model, target_weights, params, log):
  assert fmodels.fmodel_neutron() is None # not implemented
  assert model.ias_selection is None # tardy+ias is not a useful combination
  xs = fmodels.fmodel_xray().xray_structure
  sites_cart_start = xs.sites_cart()
  sites = matrix.col_list(sites_cart_start)
  tt = tardy_tree.construct(
    sites=sites,
    edge_list=model.restraints_manager.geometry.simple_edge_list(
      omit_slack_greater_than=params.omit_bonds_with_slack_greater_than))
  tt.finalize()
  for i,j in tt.collinear_bonds_edge_list:
    s = xs.scatterers()
    print >> log, "tardy collinear bond:", s[i].label
    print >> log, "                     ", s[j].label
  potential_obj = potential_object(
    nonbonded_attenuation_factor=params.nonbonded_attenuation_factor,
    fmodels=fmodels,
    model=model,
    target_weights=target_weights)
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
  temp_as_e_kin = 0.5 * sim.degrees_of_freedom * boltzmann_constant_akma
  e_kin_cap = params.temperature_cap_kelvin * temp_as_e_kin
  time_step_akma = params.time_step_pico_seconds / akma_time_as_pico_seconds
  print >> log, "tardy dynamics:"
  for i_time_step in xrange(params.number_of_time_steps):
    sim.dynamics_step(delta_t=time_step_akma, e_kin_cap=e_kin_cap)
    print >> log, "  time step: %3d * %7.5f pico seconds," % (
      i_time_step+1, params.time_step_pico_seconds),
    print >> log, "rms: %8.4f, temperature: %7.2f K" % (
      xs.sites_cart().rms_difference(sites_cart_start),
      sim.e_kin/temp_as_e_kin),
    ekbvs = sim.e_kin_before_velocity_scaling
    if (ekbvs is not None and ekbvs != sim.e_kin):
      print >> log, "(before velocity scaling: %7.2f K)" % (
        ekbvs/temp_as_e_kin),
    print >> log
  if (params.minimization_max_iterations > 0):
    print >> log, "tardy gradient-driven minimization:"
    def show_rms(minimizer=None):
      print >> log, "  rms: %8.4f" % (
        xs.sites_cart().rms_difference(sites_cart_start))
    sim.minimization(
      max_iterations=params.minimization_max_iterations,
      callback_after_step=show_rms)
  print >> log, "After tardy minimization:"
  show_rms()
