from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx import geometry_restraints
import scitbx.lbfgs
from cctbx import xray
import mmtbx.utils
from scitbx.array_family import flex
import random

# For lbfgs class
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
import math

if (1): # fixed random seed
  random.seed(1)
  flex.set_random_seed(1)

class lbfgs(object):

  def __init__(self,
               xray_structure,
               geometry_restraints,
               states,
               max_iterations = 100,
               min_iterations = 0,
               verbose = 0,
               correct_special_position_tolerance = 1.0):
    adopt_init_args(self, locals())
    self.f=None
    self.correct_special_position_tolerance = correct_special_position_tolerance
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
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
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)
    self.states.add(sites_cart = self.xray_structure.sites_cart())

  def compute_target(self, compute_gradients):
    target_and_grads = self.geometry_restraints.energies_sites(
      sites_cart = self.xray_structure.sites_cart(),
      compute_gradients = True)
    self.f = target_and_grads.target
    if(compute_gradients):
      self.g = target_and_grads.gradients.as_double()

  def callback_after_step(self, minimizer):
    if(self.verbose > 0):
      print("refinement.minimization step: f,iter,nfun:", end=' ')
      print(self.f,minimizer.iter(),minimizer.nfun())

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    self.compute_target(compute_gradients = True)
    if(self.verbose > 1):
      print("xray.minimization line search: f,rms(g):", end=' ')
      print(self.f, math.sqrt(flex.mean_sq(self.g)))
    return self.f, self.g

pdb_str = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A   7       9.035   7.190   5.709  1.00 15.00           N
ATOM      2  CA  TYR A   7      10.069   7.387   6.719  1.00 15.00           C
ATOM      3  C   TYR A   7      11.456   7.352   6.086  1.00 15.00           C
ATOM      4  O   TYR A   7      11.630   6.849   4.976  1.00 15.00           O
ATOM      5  CB  TYR A   7       9.962   6.321   7.810  1.00 15.00           C
ATOM      6  CG  TYR A   7       8.638   6.326   8.541  1.00 15.00           C
ATOM      7  CD1 TYR A   7       7.575   5.550   8.098  1.00 15.00           C
ATOM      8  CD2 TYR A   7       8.451   7.106   9.674  1.00 15.00           C
ATOM      9  CE1 TYR A   7       6.363   5.551   8.762  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.243   7.114  10.345  1.00 15.00           C
ATOM     11  CZ  TYR A   7       6.203   6.335   9.885  1.00 15.00           C
ATOM     12  OH  TYR A   7       4.998   6.339  10.550  1.00 15.00           O
ATOM     13  HA  TYR A   7       9.948   8.256   7.132  1.00 15.00           H
ATOM     14  HB2 TYR A   7      10.073   5.446   7.405  1.00 15.00           H
ATOM     15  HB3 TYR A   7      10.662   6.472   8.464  1.00 15.00           H
ATOM     16  HD1 TYR A   7       7.680   5.020   7.340  1.00 15.00           H
ATOM     17  HD2 TYR A   7       9.151   7.632   9.987  1.00 15.00           H
ATOM     18  HE1 TYR A   7       5.660   5.026   8.454  1.00 15.00           H
ATOM     19  HE2 TYR A   7       7.133   7.641  11.103  1.00 15.00           H
ATOM     20  HH  TYR A   7       5.037   6.856  11.211  1.00 15.00           H
TER
END
"""

def run():
  # Read and process PDB file
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
  #
  xrs_dc = xray_structure.deep_copy_scatterers()
  pdb_hierarchy.write_pdb_file(
    file_name        = "start.pdb",
    crystal_symmetry = xray_structure.crystal_symmetry())
  # Create grm
  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)
  states = mmtbx.utils.states(
    xray_structure = xray_structure,
    pdb_hierarchy  = pdb_hierarchy)
  states.add(sites_cart = xray_structure.sites_cart())
  # shake coordinates
  xray_structure.shake_sites_in_place(rms_difference=5.0)
  pdb_hierarchy.adopt_xray_structure(xray_structure)
  pdb_hierarchy.write_pdb_file(
    file_name        = "distorted.pdb",
    crystal_symmetry = xray_structure.crystal_symmetry())
  states.add(sites_cart = xray_structure.sites_cart())
  #
  xray_structure.scatterers().flags_set_grads(state=False)
  xray_structure.scatterers().flags_set_grad_site(
    iselection = xray_structure.all_selection().iselection())
  minimized = lbfgs(
    xray_structure      = xray_structure,
    states              = states,
    geometry_restraints = geometry_restraints,
    verbose             = 0)
  minimized.states.write(
    file_name        = "minimized_all_states.pdb",
    crystal_symmetry = xray_structure.crystal_symmetry())
  pdb_hierarchy.adopt_xray_structure(minimized.xray_structure)
  pdb_hierarchy.write_pdb_file(
    file_name        = "minimized.pdb",
    crystal_symmetry = xray_structure.crystal_symmetry())

if (__name__ == "__main__"):
  run()
