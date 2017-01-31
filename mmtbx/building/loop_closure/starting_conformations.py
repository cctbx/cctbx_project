from __future__ import division

from mmtbx.building.loop_closure import utils
from mmtbx.validation import ramalyze
import itertools
from libtbx.utils import null_out

import boost.python
ext = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

def set_rama_angles(moving_h, angles, direction_forward=True):
  """
  angles = [(phi, psi), (phi, psi), ... (phi, psi)]
  phi or psi == None means we don't change this angle
  returns deep-copied hierarchy with new angles. Change occurs from first to
  last angle so starting point would be in the same place.
  This function should produce up to all possible favored conformations.
  This function doesn't change moving_h
  direction_forward==True - set from beginning to end - the end residue moves
  direction_forward==False - set from end to beginning, the first residue moves
  """
  # print "angles", angles
  # STOP()
  result_h = moving_h.deep_copy()
  result_h.reset_atom_i_seqs()
  phi_psi_atoms = utils.get_phi_psi_atoms(moving_h)
  assert len(phi_psi_atoms) == len(angles)
  if not direction_forward:
    phi_psi_atoms.reverse()
    angles.reverse()
  for ps_atoms, target_angle_pair in zip(phi_psi_atoms, angles):
    phi_psi_pair = ps_atoms[0]
    phi_psi_angles = utils.get_pair_angles(phi_psi_pair)
    # print "ps_atoms, target_angle_pair", phi_psi_angles, target_angle_pair
    # phi
    if target_angle_pair[0] is not None:
      rotation_angle = -phi_psi_angles[0]+target_angle_pair[0]
      # if not direction_forward:
      #   rotation_angle = -rotation_angle
      utils.rotate_atoms_around_bond(
          result_h,
          phi_psi_pair[0][1],
          phi_psi_pair[0][2],
          angle=rotation_angle,
          direction_forward=direction_forward)
    # psi
    if target_angle_pair[1] is not None:
      rotation_angle = -phi_psi_angles[1]+target_angle_pair[1]
      # if not direction_forward:
      #   rotation_angle = -rotation_angle
      utils.rotate_atoms_around_bond(
          result_h,
          phi_psi_pair[1][1],
          phi_psi_pair[1][2],
          angle=rotation_angle,
          direction_forward=direction_forward)
  # print utils.list_rama_outliers_h(result_h)
  # result_h.write_pdb_file(file_name="variant_%s.pdb" % direction_forward)
  # STOP()
  return result_h

def is_not_none_combination(comb):
  for pair in comb:
    if pair != (None, None):
      return True
  return False

def get_sampled_rama_favored_angles(rama_key, r=None, step=20):
  if r is None:
    r = rama_eval()
  result = []
  for i in range(-180, 180, step):
    for j in range(-180, 180, step):
      score = r.evaluate_angles(ramalyze.res_types[rama_key], i,j)
      r_ev = ramalyze.ramalyze.evalScore(ramalyze.res_types[rama_key], score)
      if r_ev == ramalyze.RAMALYZE_FAVORED:
        result.append((i,j))
  return result

# Refactoring idea: combine these two functions
def get_all_starting_conformations(moving_h, change_radius,
    direction_forward=True, cutoff=50, log=null_out()):
  variants = []
  r = rama_eval()
  phi_psi_atoms = utils.get_phi_psi_atoms(moving_h)
  n_rama = len(phi_psi_atoms)
  change_angles = range((n_rama)//2-change_radius, (n_rama)//2+change_radius+1)
  # print "  change_angles", change_angles
  for i, (phi_psi_pair, rama_key) in enumerate(phi_psi_atoms):
    if i in change_angles or (utils.rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_OUTLIER):
      if utils.rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_OUTLIER:
        vs = get_sampled_rama_favored_angles(rama_key, r)
      else:
        vs = ramalyze.get_favored_regions(rama_key)
      variants.append(vs)
      # variants.append(ramalyze.get_favored_regions(rama_key))
    else:
      variants.append([(None, None)])
  print >> log, "variants", variants
  all_angles_combination = list(itertools.product(*variants))
  result = []
  i = 0
  n_added = 0
  n_all_combination = len(all_angles_combination)
  i_max = min(cutoff, n_all_combination)
  while n_added < i_max:
    comb = all_angles_combination[i]
    if is_not_none_combination(comb):
      result.append(set_rama_angles(moving_h, list(comb),direction_forward=direction_forward))
      print >> log, "Model %d, angles:" % i, comb
      n_added += 1
    i += 1
  # STOP()
  return result

def get_starting_conformations(moving_h, direction_forward=True, cutoff=50, log=null_out()):
  """
  modify only ramachandran outliers.
  """
  variants = []
  r = rama_eval()
  phi_psi_atoms = utils.get_phi_psi_atoms(moving_h)
  for phi_psi_pair, rama_key in phi_psi_atoms:
    if (utils.rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_OUTLIER):
      vs = get_sampled_rama_favored_angles(rama_key, r)
      # print len(vs)
      # print vs
      # STOP()
      variants.append(vs)
      # variants.append(ramalyze.get_favored_regions(rama_key))
    else:
      variants.append([(None, None)])
  result = []
  print >> log, "variants", variants
  if variants.count([(None, None)]) == len(variants):
    print "Nothing to CCD"
    return result
  all_angles_combination = list(itertools.product(*variants))
  i = 0
  n_added = 0
  n_all_combination = len(all_angles_combination)
  i_max = min(cutoff, n_all_combination)
  while n_added < i_max:
    comb = all_angles_combination[i]
    if is_not_none_combination(comb):
      result.append(set_rama_angles(moving_h, list(comb),direction_forward=direction_forward))
      print >> log, "Model %d, angles:" % i, comb, direction_forward
      n_added += 1
      # result[-1].write_pdb_file(file_name="start_h_%s_%d.pdb" % (direction_forward, i))
    i += 1
  # STOP()
  return result
