from __future__ import absolute_import, division, print_function

from mmtbx.building.loop_closure import utils
from mmtbx.validation import ramalyze
import itertools
import numpy
import random
from libtbx.utils import null_out

import boost_adaptbx.boost.python as bp
from six.moves import zip
from six.moves import range
ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

from six.moves import cStringIO as StringIO

def set_rama_angles(moving_h, angles, direction_forward=True, check_omega=False):
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
  fixed_omega = False
  phi_psi_atoms = utils.get_phi_psi_atoms(moving_h, omega=True)
  assert len(phi_psi_atoms) == len(angles), "%d != %d" % (len(phi_psi_atoms), len(angles))
  if not direction_forward:
    phi_psi_atoms.reverse()
    angles.reverse()
  for ps_atoms, target_angle_pair in zip(phi_psi_atoms, angles):
    phi_psi_pair = ps_atoms[0]
    # print "phi_psi_pair", phi_psi_pair
    omega = ps_atoms[2]
    phi_psi_angles = utils.get_pair_angles(phi_psi_pair)
    # print "ps_atoms, target_angle_pair", phi_psi_angles, target_angle_pair
    # phi
    if target_angle_pair[0] is not None and phi_psi_angles[0] is not None:
      rotation_angle = -phi_psi_angles[0]+target_angle_pair[0]
      # print "rot angle", rotation_angle
      # if not direction_forward:
      #   rotation_angle = -rotation_angle
      utils.rotate_atoms_around_bond(
          result_h,
          phi_psi_pair[0][1],
          phi_psi_pair[0][2],
          angle=rotation_angle,
          direction_forward=direction_forward)
    # psi
    if target_angle_pair[1] is not None and phi_psi_angles[1] is not None:
      rotation_angle = -phi_psi_angles[1]+target_angle_pair[1]
      # print "rot angle", rotation_angle
      # if not direction_forward:
      #   rotation_angle = -rotation_angle
      utils.rotate_atoms_around_bond(
          result_h,
          phi_psi_pair[1][1],
          phi_psi_pair[1][2],
          angle=rotation_angle,
          direction_forward=direction_forward)
    # omega
    if omega is not None and abs(abs(omega)-180) > 10 and check_omega:
      rotation_angle= -omega+180
      # print "Omega rotation:", omega, rotation_angle
      utils.rotate_atoms_around_bond(
          result_h,
          phi_psi_pair[0][0],
          phi_psi_pair[0][1],
          angle=rotation_angle,
          direction_forward=direction_forward)
      fixed_omega = True
  # print utils.list_rama_outliers_h(result_h)
  # result_h.write_pdb_file(file_name="variant_%s.pdb" % direction_forward)
  # STOP()
  return result_h, fixed_omega

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

def get_all_starting_conformations(moving_h, change_radius,
    include_allowed, n_outliers,
    direction_forward=True, cutoff=50, change_all=True, log=null_out(), check_omega=False):
  if log is None:
    log = StringIO()
  variants = []
  result = []
  r = rama_eval()
  phi_psi_atoms = utils.get_phi_psi_atoms(moving_h, omega=True)
  # print "N residue groups in h", [x.resseq for x in moving_h.residue_groups()]
  if len(phi_psi_atoms) == 0:
    print("Strange input to starting conformations!!!", file=log)
    return result
  n_rama = len(phi_psi_atoms)
  # print "n_rama", n_rama
  change_angles = [None]
  if change_all:
    change_angles = range((n_rama)//2-change_radius-n_outliers//2, (n_rama)//2+change_radius+1+n_outliers//2)
    # if change_angles[0] < 0:
    #   change_angles = range(change_angles[-1]-change_angles[0])
  has_twisted = False
  if check_omega:
    omegas = [x[2] for x in phi_psi_atoms]
    for o in omegas:
      if o is not None and abs(abs(o)-180) > 30:
        has_twisted = True
  print("n_outliers", n_outliers, file=log)
  for i, (phi_psi_pair, rama_key, omega) in enumerate(phi_psi_atoms):
    angle_is_outlier = utils.rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_OUTLIER
    angle_is_outlier = angle_is_outlier or (include_allowed and utils.rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_ALLOWED)
    twisted = omega is not None and ((abs(abs(omega)-180) > 30) and check_omega)
    print("in cycle, N, outlier?, change?, twisted?", i, angle_is_outlier, i in change_angles, twisted, file=log)
    if angle_is_outlier and n_outliers < 3:
      vs = get_sampled_rama_favored_angles(rama_key, r)
    elif (i in change_angles) or angle_is_outlier or has_twisted:
      # vs = get_sampled_rama_favored_angles(rama_key, r)
      vs = ramalyze.get_favored_regions(rama_key)
    else:
      vs = [(None, None)]
    variants.append(vs)
  print("variants", variants, file=log)

  # Filtering them, since could be
  # [len(x) for x in variants] = [129, 129, 4, 129, 129]
  # resulting in 1107691524 all_angles_combination
  n_comb = numpy.prod([len(x) for x in variants])
  if n_comb > cutoff:
    # still aiming for ~1000
    n_in_each = int(1000 ** (1/len(variants)))
    variants = [random.sample(x, n_in_each) if len(x)>n_in_each else x for x in variants]
  all_angles_combination = list(itertools.product(*variants))
  # filter none combinations
  # print "len(all_angles_combination)", len(all_angles_combination)
  all_angles_combination_f = []
  for comb in all_angles_combination:
    if is_not_none_combination(comb):
      all_angles_combination_f.append(comb)
  print("len(all_angles_combination_f)", len(all_angles_combination_f), file=log)
  return all_angles_combination_f
  # if len(all_angles_combination_f) == 0:
  #   print "In starting conformations - outlier was fixed?"
  #   return result
  # n_added = 0
  # n_all_combination = len(all_angles_combination_f)
  # i_max = min(cutoff, n_all_combination)
  # assert i_max > 0
  # step = float(n_all_combination-1)/float(i_max-1)
  # if step < 1:
  #   step = 1
  # for i in range(i_max):
  #   comb = all_angles_combination_f[int(round(step*i))]
  #   result.append(set_rama_angles(moving_h, list(comb),direction_forward=direction_forward))
  #   print >> log, "Model %d, angles:" % i, comb
  # return result
