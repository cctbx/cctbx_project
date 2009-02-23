from cctbx import xray
from cctbx.array_family import flex
from scitbx.rigid_body.essence import tst_molecules
from scitbx.graph import tardy_tree
from scitbx import matrix

master_phil_str = """\
  number_of_time_steps = 10
    .type = int
  time_step = 0.001
    .type = float
  minimization_max_iterations = 10
    .type = int
"""

class potential_object(object):

  def __init__(O, fmodels, model, target_weights):
    O.fmodels = fmodels
    O.model = model
    O.weights = target_weights.xyz_weights_result
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

def run(fmodels, model, target_weights, params):
  assert fmodels.fmodel_neutron() is None # not implemented
  assert model.ias_selection is None # tardy+ias is not a useful combination
  sst = model.restraints_manager.geometry.shell_sym_tables[0]
  tt = tardy_tree.construct(
    edge_list=sst.simple_edge_list(), n_vertices=sst.size()).finalize()
  xs = fmodels.fmodel_xray().xray_structure
  potential_obj = potential_object(
    fmodels=fmodels, model=model, target_weights=target_weights)
  sites_cart_start = xs.sites_cart()
  sites = matrix.col_list(sites_cart_start)
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
  def show_rms(minimizer=None):
    print "rms:", xs.sites_cart().rms_difference(sites_cart_start)
  for i_time_step in xrange(params.number_of_time_steps):
    print "tardy time step:", i_time_step
    sim.dynamics_step(delta_t=params.time_step)
    show_rms()
  if (params.minimization_max_iterations > 0):
    sim.minimization(
      max_iterations=params.minimization_max_iterations,
      callback_after_step=show_rms)
  print "After minimization:"
  show_rms()
