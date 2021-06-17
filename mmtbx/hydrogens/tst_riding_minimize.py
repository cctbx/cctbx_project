from __future__ import absolute_import, division, print_function
import time
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx import geometry_restraints
import scitbx.lbfgs
import mmtbx.utils
# For lbfgs class:
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
import math
# for riding H:
from mmtbx.hydrogens import riding

class lbfgs(object):

  def __init__(self,
               xray_structure,
               geometry_restraints,
               states,
               riding_h_manager,
               use_riding,
               verbose,
               max_iterations = 1000,
               min_iterations = 0,
               correct_special_position_tolerance = 1.0):
    adopt_init_args(self, locals())
    # f = refinement target
    self.f = None
    self.correct_special_position_tolerance = correct_special_position_tolerance
    # self.x = coordinate shifts
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self.h_parameterization = self.riding_h_manager.h_parameterization
    self._scatterers_start = self.xray_structure.scatterers()
    # -----------------------------------------------------
    if self.use_riding:
      self.hd_selection = self.xray_structure.hd_selection()
      self.x = flex.double(
        self.xray_structure.select(~self.hd_selection).n_parameters(), 0)
      self._scatterers_start = self.xray_structure.scatterers().select(
        ~self.hd_selection)
    # -----------------------------------------------------
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations,
      min_iterations = min_iterations)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    self.compute_target(compute_gradients = False)

  def apply_shifts(self):
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers_shifted[i_seq].site = crystal.correct_special_position(
        crystal_symmetry = self.xray_structure,
        special_op       = site_symmetry_table.get(i_seq).special_op(),
        site_frac        = scatterers_shifted[i_seq].site,
        site_label       = scatterers_shifted[i_seq].label,
        tolerance        = self.correct_special_position_tolerance)
    # -----------------------------------------------------
    if self.use_riding:
      new_sites = self.xray_structure.sites_frac().set_selected(
        ~self.hd_selection, scatterers_shifted.extract_sites())
      self.xray_structure.set_sites_frac(new_sites)
      self.riding_h_manager.idealize_riding_h_positions(
          xray_structure = self.xray_structure)
    else:
      self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)
    self.states.add(sites_cart = self.xray_structure.sites_cart())

  def compute_target(self, compute_gradients):
    target_and_grads = self.geometry_restraints.energies_sites(
      sites_cart = self.xray_structure.sites_cart(),
      compute_gradients = True)
    self.f = target_and_grads.target
    if(compute_gradients):
      self.grads = target_and_grads.gradients
    # -----------------------------------------------------
    # if use_riding hydrogen, modify the gradients
      if self.use_riding:
        self.grads = self.riding_h_manager.gradients_reduced_cpp(
          gradients    = self.grads,
          sites_cart   = self.xray_structure.sites_cart(),
          hd_selection = self.hd_selection)
    # -----------------------------------------------------
      self.g = self.grads.as_double()

  def callback_after_step(self, minimizer):
    if(self.verbose > 0):
      print("refinement.minimization step: f,iter,nfun:", self.f,minimizer.iter(),minimizer.nfun())

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    self.compute_target(compute_gradients = True)
    if(self.verbose > 1):
      print ("xray.minimization line search: f,rms(g):",
        self.f, math.sqrt(flex.mean_sq(self.g)))
    return self.f, self.g

# distorted Tyrosine residue
pdb_str = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      1  N   TYR A   7       9.035   7.089   5.637  1.00 15.00           N
ATOM      2  CA  TYR A   7       9.887   7.460   6.606  1.00 15.00           C
ATOM      3  C   TYR A   7      11.285   7.169   6.250  1.00 15.00           C
ATOM      4  O   TYR A   7      11.563   6.972   4.911  1.00 15.00           O
ATOM      5  CB  TYR A   7       9.830   6.306   7.775  1.00 15.00           C
ATOM      6  CG  TYR A   7       8.528   6.445   8.572  1.00 15.00           C
ATOM      7  CD1 TYR A   7       7.544   5.558   7.939  1.00 15.00           C
ATOM      8  CD2 TYR A   7       8.348   7.014   9.545  1.00 15.00           C
ATOM      9  CE1 TYR A   7       6.338   5.396   8.586  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.068   7.273  10.230  1.00 15.00           C
ATOM     11  CZ  TYR A   7       6.090   6.501   9.765  1.00 15.00           C
ATOM     12  OH  TYR A   7       4.857   6.463  10.553  1.00 15.00           O
ATOM     13  HA  TYR A   7       9.831   8.340   7.271  1.00 15.00           H
ATOM     14  HB2 TYR A   7       9.943   5.616   7.588  1.00 15.00           H
ATOM     15  HB3 TYR A   7      10.644   6.563   8.354  1.00 15.00           H
ATOM     16  HD1 TYR A   7       7.598   4.974   7.330  1.00 15.00           H
ATOM     17  HD2 TYR A   7       9.064   7.637  10.081  1.00 15.00           H
ATOM     18  HE1 TYR A   7       5.756   5.132   8.319  1.00 15.00           H
ATOM     19  HE2 TYR A   7       7.206   7.495  11.005  1.00 15.00           H
ATOM     20  HH  TYR A   7       5.208   7.016  11.256  1.00 15.00           H
TER
"""

def run():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = None,
    raw_records    = pdb_str,
    force_symmetry = True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xray_structure = processed_pdb_file.xray_structure()
  #pdb_hierarchy.write_pdb_file(
  #  file_name        = "distorted.pdb",
  #  crystal_symmetry = xray_structure.crystal_symmetry())
  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)
  #geometry_restraints.write_geo_file(file_name='start.geo')
  states = mmtbx.utils.states(
    pdb_hierarchy  = pdb_hierarchy)
  states.add(sites_cart = xray_structure.sites_cart())

  riding_h_manager = riding.manager(
    pdb_hierarchy       = pdb_hierarchy,
    geometry_restraints = geometry_restraints)

  states.add(sites_cart = xray_structure.sites_cart())
  #xray_structure.scatterers().flags_set_grads(state=False)
  xray_structure.scatterers().flags_set_grad_site(
    iselection = xray_structure.all_selection().iselection())

  use_riding = True
  minimized = lbfgs(
    xray_structure      = xray_structure,
    states              = states,
    geometry_restraints = geometry_restraints,
    riding_h_manager    = riding_h_manager,
    use_riding          = use_riding,
    verbose             = 0)
  #minimized.states.write(
  #  file_name        = "minimized_all_states.pdb",
  #  crystal_symmetry = xray_structure.crystal_symmetry())
  pdb_hierarchy.adopt_xray_structure(minimized.xray_structure)
  #pdb_hierarchy.write_pdb_file(
  #  file_name        = "minimized.pdb",
  #  crystal_symmetry = xray_structure.crystal_symmetry())

  target_final = geometry_restraints.energies_sites(
    sites_cart = xray_structure.sites_cart(),
    compute_gradients = True).target
  print('final target', target_final)
  assert (target_final < 1.5), 'Target of final riding model is too large'

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time:", round(time.time()-t0, 2), "seconds")
