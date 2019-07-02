from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.real_space.rigid_body
import scitbx.lbfgs
import iotbx.pdb

def run(pdb_hierarchy,
        target_map,
        unit_cell,
        real_space_gradients_delta,
        max_allowed_shift = 1.5,
        max_iterations = 50,
        log = None):
  lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
    max_iterations = max_iterations)
  get_class = iotbx.pdb.common_residue_names_get_class
  def target(target_map, sites_cart, unit_cell):
    sites_frac = unit_cell.fractionalize(sites_cart)
    result = 0
    for site_frac in sites_frac:
      result += target_map.eight_point_interpolation(site_frac)
    return result
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            atoms = residue.atoms()
            if(get_class(name=residue.resname)=="common_water" and len(atoms)>1):
              if(log is not None):
                print("chain %s resname %s resseq %s"%(
                  chain.id, residue.resname, residue.resseq), file=log)
              sites_cart_start = atoms.extract_xyz()
              target_start = target(target_map, sites_cart_start, unit_cell)
              if(log is not None):
                print("  target_start: %6.4f"%target_start, file=log)
              target_current = target_start
              sites_cart_best = sites_cart_start.deep_copy()
              shift_range = [-0.3,0,0.3]
              for x_shift in shift_range:
                for y_shift in shift_range:
                  for z_shift in shift_range:
                    shift = flex.vec3_double(
                      [(x_shift, y_shift, z_shift)]*sites_cart_start.size())
                    sites_cart = sites_cart_start + shift
                    residue.atoms().set_xyz(sites_cart)
                    minimized = mmtbx.refinement.real_space.rigid_body.refine(
                      residue                     = residue,
                      density_map                 = target_map,
                      geometry_restraints_manager = None,
                      real_space_target_weight    = 1,
                      real_space_gradients_delta  = real_space_gradients_delta,
                      lbfgs_termination_params    = lbfgs_termination_params,
                      unit_cell                   = unit_cell)
                    sites_cart = minimized.sites_cart_residue
                    distance_moved = flex.mean(
                      flex.sqrt((sites_cart - sites_cart_start).dot()))
                    t = target(target_map, sites_cart, unit_cell)
                    if(t>=target_current and distance_moved<max_allowed_shift):
                      sites_cart_best = sites_cart.deep_copy()
                      target_current = t
              residue.atoms().set_xyz(sites_cart_best)
              target_final = target(target_map, sites_cart_best, unit_cell)
              distance_moved = flex.mean(
                flex.sqrt((sites_cart_best - sites_cart_start).dot()))
              if(log is not None):
                print("  target_final: %6.4f" % target_final, file=log)
                print("  dist. moved : %6.4f" % distance_moved, file=log)
