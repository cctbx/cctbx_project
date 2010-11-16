# LIBTBX_SET_DISPATCHER_NAME phenix.geometry_minimization

import iotbx.phil
from cctbx import geometry_restraints
import cctbx.geometry_restraints.lbfgs
import scitbx.lbfgs
import sys

master_params = iotbx.phil.parse("""\
  alternate_nonbonded_off_on=False
    .type = bool
  max_iterations=500
    .type = int
  macro_cycles=1
    .type = int
  show_geometry_restraints=False
    .type = bool
""")

class lbfgs(geometry_restraints.lbfgs.lbfgs):

  def __init__(self,
        sites_cart,
        geometry_restraints_manager,
        geometry_restraints_flags,
        lbfgs_termination_params,
        sites_cart_selection=None,
        lbfgs_exception_handling_params=None,
        rmsd_bonds_termination_cutoff=0,
        rmsd_angles_termination_cutoff=0):
    self.rmsd_bonds_termination_cutoff = rmsd_bonds_termination_cutoff
    self.rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff
    geometry_restraints.lbfgs.lbfgs.__init__(self,
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=lbfgs_termination_params,
      sites_cart_selection=sites_cart_selection,
      lbfgs_exception_handling_params=lbfgs_exception_handling_params)

  def callback_after_step(self, minimizer):
    self.apply_shifts()
    if([self.rmsd_angles, self.rmsd_bonds].count(None) == 0):
      if(self.rmsd_angles < self.rmsd_angles_termination_cutoff and
         self.rmsd_bonds < self.rmsd_bonds_termination_cutoff):
        return True

def run(processed_pdb_file, params = master_params.extract(), log =sys.stdout):
  co = params
  geometry_restraints_flags = geometry_restraints.flags.flags(default=True)
  all_chain_proxies = processed_pdb_file.all_chain_proxies
  geometry_restraints_manager = processed_pdb_file.\
    geometry_restraints_manager(show_energies = False)
  special_position_settings = all_chain_proxies.special_position_settings
  sites_cart = all_chain_proxies.sites_cart_exact().deep_copy()
  atom_labels = [atom.id_str() for atom in all_chain_proxies.pdb_atoms]
  geometry_restraints_manager.site_symmetry_table \
    .show_special_position_shifts(
      special_position_settings=special_position_settings,
      site_labels=atom_labels,
      sites_cart_original=all_chain_proxies.sites_cart,
      sites_cart_exact=sites_cart,
      out=log,
      prefix="  ")
  if (co.show_geometry_restraints):
    geometry_restraints_manager.show_sorted(
      flags=geometry_restraints_flags,
      sites_cart=sites_cart,
      site_labels=atom_labels)
  pair_proxies =  geometry_restraints_manager.pair_proxies(
    sites_cart=all_chain_proxies.sites_cart,
    flags=geometry_restraints_flags)
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=atom_labels,
    f=log,
    max_items=10)
  if (pair_proxies.nonbonded_proxies is not None):
    pair_proxies.nonbonded_proxies.show_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      site_labels=atom_labels,
      f=log,
      max_items=10)
  del pair_proxies
  print
  log.flush()
  if (co.alternate_nonbonded_off_on and co.macro_cycles % 2 != 0):
    co.macro_cycles += 1
    print "INFO: Number of macro cycles increased by one to ensure use of"
    print "      nonbonded interactions in last macro cycle."
    print
  for i_macro_cycle in xrange(co.macro_cycles):
    if (co.alternate_nonbonded_off_on):
      geometry_restraints_flags.nonbonded = bool(i_macro_cycle % 2)
      print "Use nonbonded interactions this macro cycle:", \
        geometry_restraints_flags.nonbonded
    minimized = lbfgs(
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=co.max_iterations))
    print "Energies at start of minimization:"
    minimized.first_target_result.show()
    print
    print "Number of minimization iterations:", minimized.minimizer.iter()
    print "Root-mean-square coordinate difference: %.3f" % (
      all_chain_proxies.sites_cart.rms_difference(sites_cart))
    print
    print "Energies at end of minimization:"
    minimized.final_target_result.show()
    print
    geometry_restraints_manager.pair_proxies(
      sites_cart=sites_cart,
      flags=geometry_restraints_flags) \
        .bond_proxies.show_sorted(
          by_value="residual",
          sites_cart=sites_cart,
          site_labels=atom_labels,
          f=log,
          max_items=10)
    print
  assert geometry_restraints_flags.nonbonded
  return sites_cart

if (__name__ == "__main__"):
  # XXX temporary message
  print "\n***Command chage: use phenix.pdbtools for geometry_minimization.***\n"
