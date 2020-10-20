from __future__ import absolute_import, division, print_function

from libtbx import adopt_init_args
import mmtbx.utils
from mmtbx.building.loop_closure import utils
from mmtbx.validation.ramalyze import ramalyze, RAMALYZE_OUTLIER, \
    RAMALYZE_ALLOWED, RAMALYZE_FAVORED
from scitbx.matrix import project_point_on_axis
import math
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

import boost_adaptbx.boost.python as bp
from six.moves import zip,range
ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

ext2 = bp.import_ext("mmtbx_building_loop_closure_ext")
from mmtbx_building_loop_closure_ext import ccd_cpp

@bp.inject_into(ccd_cpp)
class _():

  def run(self, direction_forward=True, save_states=False, avoid_allowed_region=False):
    if save_states:
      self.states = mmtbx.utils.states(pdb_hierarchy=self.moving_h)
      self.states.add(sites_cart=self.moving_h.atoms().extract_xyz())

    phi_psi_atoms = utils.get_phi_psi_atoms(self.moving_h)
    if not direction_forward:
      phi_psi_atoms.reverse()

    # here we can start a ccd cycle
    self.n_iter = 0
    rmsd_good = 1000
    previous_rmsd = 1000
    self.early_exit = False
    # self.moving_h.write_pdb_file(file_name="start_ccd.pdb")
    while (rmsd_good > self.needed_rmsd and
        self.n_iter <= self.max_number_of_iterations and not self.early_exit):
      # print_rama_stats(phi_psi_atoms, r)
      # for phi_psi_pair in phi_psi_atoms[:-1]:

      # check rama again separately before the cycle
      # list_rama_outliers(phi_psi_atoms, r)

      for phi_psi_pair, rama_key in phi_psi_atoms:
        before_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
        rama_score = before_rama_score
        # print "rama score:", rama_score, "->",
        for i, atoms in enumerate(phi_psi_pair):
          # current phi-psi angles:
          # find the optimal angle
          if atoms is None:
            continue
          if direction_forward:
            ccd_angle = self._find_angle(atoms[1].xyz, atoms[2].xyz)
          else:
            ccd_angle = self._find_angle(atoms[2].xyz, atoms[1].xyz)
          # print "phi_psi_angles", phi_psi_angles
          # rama_score = r.evaluate("general", phi_psi_angles)
          # print "rama_score", rama_score


          angle_modified = self._modify_angle(ccd_angle)
          # angle_modified = ccd_angle
          # print ("  ccd_angle", ccd_angle, angle_modified)
          # angle_modified = - angle_modified

          phi_psi_angles = utils.get_pair_angles(phi_psi_pair)
          # print "phi_psi_angles", phi_psi_angles
          before_rotation_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
          if (ramalyze.evalScore(rama_key, before_rotation_rama_score) == RAMALYZE_OUTLIER
            or (avoid_allowed_region and ramalyze.evalScore(rama_key, before_rotation_rama_score) == RAMALYZE_ALLOWED)):
            # assert i == 0
            if i != 0:
              # this is a error, we should spot rama outliers on the first angle
              print("i", i)
              print(pair_info(phi_psi_pair))
              print("rama_key", rama_key)
              print("before_rotation_rama_score", before_rotation_rama_score, end=' ')
              print(ramalyze.evalScore(rama_key, before_rotation_rama_score))
              break

            # correct it to the nearest non-outlier region
            target_phi_psi = utils.find_nearest_non_outlier_region(phi_psi_pair, self.r, rama_key)
            # print "For outlier:", phi_psi_angles, target_phi_psi
            # here we want to correct outlier regardless the target function
            # outcome and proceed to the next phi-psi pair
            now_psi_angle0 = utils.get_dihedral_angle(phi_psi_pair[1])
            utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
                angle=-phi_psi_angles[0]+target_phi_psi[0],direction_forward=direction_forward)
            # now psi angle
            now_psi_angle = utils.get_dihedral_angle(phi_psi_pair[1])

            # print "psi angles:", now_psi_angle0, now_psi_angle
            angles_ok = (approx_equal(now_psi_angle0-now_psi_angle, 0, eps=1e-4, out=null_out()) or
                approx_equal(now_psi_angle0-now_psi_angle, 360, eps=1e-4, out=null_out()) or
                approx_equal(now_psi_angle0-now_psi_angle, -360, eps=1e-4, out=null_out()))
            if not angles_ok:
                approx_equal(now_psi_angle0-now_psi_angle, 0, eps=1e-4)
                approx_equal(now_psi_angle0-now_psi_angle, 360, eps=1e-4)
                approx_equal(now_psi_angle0-now_psi_angle, -360, eps=1e-4)
            assert angles_ok
            # approx_equal(now_psi_angle0, now_psi_angle)
            # assert now_psi_angle0 == now_psi_angle
            utils.rotate_atoms_around_bond(self.moving_h, atoms[2], atoms[3],
                angle=-now_psi_angle+target_phi_psi[1], direction_forward=direction_forward)

            # approx_equal(utils.get_dihedral_angle(phi_psi_pair[0]), target_phi_psi[0])
            # approx_equal(utils.get_dihedral_angle(phi_psi_pair[1]), target_phi_psi[1])
            resulting_rama_ev = utils.rama_evaluate(phi_psi_pair, self.r, rama_key)
            # print "evaluation:", resulting_rama_ev, RAMALYZE_FAVORED
            assert resulting_rama_ev == RAMALYZE_FAVORED, resulting_rama_ev
            break # we are done with this phi_psi_pair
          # rotate the whole thing around
          utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
              angle=angle_modified, direction_forward=direction_forward)
          after_rotation_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
          # print "before/after rotation rama:", before_rotation_rama_score, after_rotation_rama_score
          # if before_rotation_rama_score > after_rotation_rama_score:
          eval_score_after_rotation = ramalyze.evalScore(rama_key, after_rotation_rama_score)
          if eval_score_after_rotation == RAMALYZE_OUTLIER or \
             eval_score_after_rotation == RAMALYZE_ALLOWED:
            # rotate back!!! / not always
            # print "  rotate back"
            if True: # always
              utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
                  angle=-angle_modified,direction_forward=direction_forward)
          s = utils.get_rama_score(phi_psi_pair, self.r, rama_key,round_coords=True)
          assert utils.rama_score_evaluate(rama_key, s) != RAMALYZE_OUTLIER, s
          if avoid_allowed_region:
            assert utils.rama_score_evaluate(rama_key, s) != RAMALYZE_ALLOWED, "%s, %s" % (s, after_rotation_rama_score)

        # new rama score:
        after_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key)
        if after_rama_score + 1e-7 < before_rama_score:
          pass
          # print "before, after", before_rama_score, after_rama_score
          # STOP()

      rmsd_good = utils.get_rmsd_xyz_fixed(
          self.fixed_ref_atoms,
          [self.moving_h.atoms()[x] for x in self.moving_ref_atoms_iseqs])
      self.resulting_rmsd = rmsd_good
      # print "n_iter, rmsd:", self.n_iter, rmsd_good
      # print get_main_chain_rmsd(moving_h, original_h)

      if save_states:
        self.states.add(sites_cart=self.moving_h.atoms().extract_xyz())
      # if n_iter % 100 == 0:
      #   moving_h.write_pdb_file(file_name="int_%d.pdb" % n_iter)
      self.n_iter += 1
      self.early_exit = abs(previous_rmsd - rmsd_good) < self.convergence_diff
      # if self.early_exit:
      #   print "  Early exit:", self.early_exit, previous_rmsd - rmsd_good
      previous_rmsd = rmsd_good
    # print "number of iterations:", n_iter
    # print_rama_stats(phi_psi_atoms, r)
    # moving_h.write_pdb_file(file_name="int_%d.pdb" % n_iter)
    # states.write(file_name="all_states.pdb")

    # return rmsd_good, states, n_iter


class ccd_python():
  def __init__(self, fixed_ref_atoms, moving_h, moving_ref_atoms_iseqs,
      max_number_of_iterations=500, needed_rmsd=0.1):
    """
    fixed_ref_atoms - list of 3 atom objects, actually, only xyz's are needed
    moving_ref_atoms_iseqs - list of 3 indeces matching atoms in
      moving_h.atoms()[<here!>].
    moving_h - hierarchy to make closure. Atom positions in it will be changed!

    """
    assert len(fixed_ref_atoms) == 3
    assert len(moving_ref_atoms_iseqs) == 3
    assert moving_h is not None
    assert moving_h.atoms_size() > 10 # arbitrary
    # adopt_init_args(self, locals())
    self.moving_h = moving_h
    self.fixed_ref_atoms = fixed_ref_atoms
    self.moving_ref_atoms_iseqs = moving_ref_atoms_iseqs
    self.max_number_of_iterations = max_number_of_iterations
    self.needed_rmsd = needed_rmsd
    self.set_modify_angle_procedure(self._modify_angle)
    self.r = rama_eval()
    # self.states = mmtbx.utils.states(pdb_hierarchy=moving_h)
    self.convergence_diff = 1e-5
    # will be bool, True if converged before max_number_of_iterations reached
    self.early_exit = None
    self.resulting_rmsd = None


  def set_modify_angle_procedure(self, procedure):
    """
    can be used to set external procedure for angle modification.
    the only argument should be angle, should return new angle in degrees
    """
    self.modify_angle_procedure = procedure

  def _modify_angle(self,angle):
    """
    change angle found by minimization. Primary use - to avoid huge turns
    in first phi-psi angles.
    """
    threshold = 1
    if abs(angle) > threshold:
      if angle > 0:
        return threshold
      else:
        return -threshold
    else:
      return angle

  @staticmethod
  def _get_f_r_s(axis_point_1,axis_point_2, moving_coor, fixed_coor):
    fc_proj = project_point_on_axis(axis_point_1, axis_point_2, fixed_coor)
    mc_proj = project_point_on_axis(axis_point_1, axis_point_2, moving_coor)
    f = (fixed_coor[0]-fc_proj[0],fixed_coor[1]-fc_proj[1],fixed_coor[2]-fc_proj[2])
    r = (moving_coor[0]-mc_proj[0],moving_coor[1]-mc_proj[1],moving_coor[2]-mc_proj[2])
    ap_21 = (axis_point_2[0]-axis_point_1[0], axis_point_2[1]-axis_point_1[1], axis_point_2[2]-axis_point_1[2])
    r_norm = math.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
    r_home = flex.vec3_double([(r[0]/r_norm, r[1]/r_norm, r[2]/r_norm)])
    ap_21_norm = math.sqrt(ap_21[0]*ap_21[0]+ap_21[1]*ap_21[1]+ap_21[2]*ap_21[2])
    theta_home = flex.vec3_double([(ap_21[0]/ap_21_norm, ap_21[1]/ap_21_norm, ap_21[2]/ap_21_norm)])
    tt = theta_home.cross(r_home)
    s_home = tt*(1/tt.norm())
    return flex.vec3_double([f]), s_home, r_norm, r_home

  def _find_angle(self, axis_point_1, axis_point_2):
    f_all = []
    s_home_all = []
    r_all = []
    r_home_all = []
    for fixed_xyz, moving_xyz in zip([x.xyz for x in self.fixed_ref_atoms],
        [self.moving_h.atoms()[x].xyz for x in self.moving_ref_atoms_iseqs]):
      f, s_home, r_norm, r_home = ccd_python._get_f_r_s(
          axis_point_1, axis_point_2, moving_xyz, fixed_xyz)
      f_all.append(f)
      s_home_all.append(s_home)
      r_all.append(r_norm)
      r_home_all.append(r_home)
    # calculating
    b = 0
    c = 0
    for i in range(3):
      b += list(2*r_all[i]*(f_all[i].dot(r_home_all[i])))[0]
      c += list(2*r_all[i]*(f_all[i].dot(s_home_all[i])))[0]
    znam = math.sqrt(b*b+c*c)
    sin_alpha = c/znam
    cos_alpha = b/znam
    alpha = math.atan2(sin_alpha, cos_alpha)
    # print "ver3 alpha:", math.degrees(alpha)
    return math.degrees(alpha)

  def run(self):
    # self.states.add(sites_cart=self.moving_h.atoms().extract_xyz())

    phi_psi_atoms = utils.get_phi_psi_atoms(self.moving_h)

    # here we can start a ccd cycle
    self.n_iter = 0
    rmsd_good = 1000
    previous_rmsd = 1000
    self.early_exit = False
    while (rmsd_good > self.needed_rmsd and
        self.n_iter <= self.max_number_of_iterations and not self.early_exit):
      # print_rama_stats(phi_psi_atoms, r)
      # for phi_psi_pair in phi_psi_atoms[:-1]:

      # check rama again separately before the cycle
      # list_rama_outliers(phi_psi_atoms, r)

      for phi_psi_pair, rama_key in phi_psi_atoms:
        before_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
        rama_score = before_rama_score
        # print "rama score:", rama_score, "->",
        for i, atoms in enumerate(phi_psi_pair):
          # current phi-psi angles:
          # find the optimal angle
          ccd_angle = self._find_angle(atoms[1].xyz, atoms[2].xyz)
          # print "phi_psi_angles", phi_psi_angles
          # rama_score = r.evaluate("general", phi_psi_angles)
          # print "rama_score", rama_score
          angle_modified = self.modify_angle_procedure(ccd_angle)

          phi_psi_angles = utils.get_pair_angles(phi_psi_pair)
          before_rotation_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
          if (ramalyze.evalScore(rama_key, before_rotation_rama_score) == RAMALYZE_OUTLIER):
              # or ramalyze.evalScore(rama_key, before_rotation_rama_score) == RAMALYZE_ALLOWED):
            # assert i == 0
            if i != 0:
              # this is a error, we should spot rama outliers on the first angle
              print("i", i)
              print(pair_info(phi_psi_pair))
              print("rama_key", rama_key)
              print("before_rotation_rama_score", before_rotation_rama_score, end=' ')
              print(ramalyze.evalScore(rama_key, before_rotation_rama_score))
              break

            # correct it to the nearest non-outlier region
            target_phi_psi = utils.find_nearest_non_outlier_region(phi_psi_pair, self.r, rama_key)
            # print "For outlier:", phi_psi_angles, target_phi_psi
            # here we want to correct outlier regardless the target function
            # outcome and proceed to the next phi-psi pair
            now_psi_angle0 = utils.get_dihedral_angle(phi_psi_pair[1])
            utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
                angle=-phi_psi_angles[0]+target_phi_psi[0])
            # now psi angle
            now_psi_angle = utils.get_dihedral_angle(phi_psi_pair[1])

            # print "psi angles:", now_psi_angle0, now_psi_angle
            angles_ok = (approx_equal(now_psi_angle0-now_psi_angle, 0) or
                approx_equal(now_psi_angle0-now_psi_angle, 360) or
                approx_equal(now_psi_angle0-now_psi_angle, -360))

            assert angles_ok
            # approx_equal(now_psi_angle0, now_psi_angle)
            # assert now_psi_angle0 == now_psi_angle
            utils.rotate_atoms_around_bond(self.moving_h, atoms[2], atoms[3],
                angle=-now_psi_angle+target_phi_psi[1])

            approx_equal(utils.get_dihedral_angle(phi_psi_pair[0]), target_phi_psi[0])
            approx_equal(utils.get_dihedral_angle(phi_psi_pair[1]), target_phi_psi[1])
            resulting_rama_ev = utils.rama_evaluate(phi_psi_pair, self.r, rama_key)
            assert resulting_rama_ev == RAMALYZE_FAVORED, resulting_rama_ev
            break # we are done with this phi_psi_pair
          # rotate the whole thing around
          utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
              angle=angle_modified)
          after_rotation_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key, round_coords=True)
          # print "before/after rotation rama:", before_rotation_rama_score, after_rotation_rama_score
          # if before_rotation_rama_score > after_rotation_rama_score:
          if ramalyze.evalScore(rama_key, after_rotation_rama_score) == RAMALYZE_OUTLIER:
            # rotate back!!! / not always
            # print "  rotate back"
            if True: # always
              utils.rotate_atoms_around_bond(self.moving_h, atoms[1], atoms[2],
                  angle=-angle_modified)
          s = utils.get_rama_score(phi_psi_pair, self.r, rama_key,round_coords=True)
          assert utils.rama_score_evaluate(rama_key, s) != RAMALYZE_OUTLIER, s

        # new rama score:
        after_rama_score = utils.get_rama_score(phi_psi_pair, self.r, rama_key)
        if after_rama_score + 1e-7 < before_rama_score:
          pass
          # print "before, after", before_rama_score, after_rama_score
          # STOP()

      rmsd_good = utils.get_rmsd(
          self.fixed_ref_atoms,
          [self.moving_h.atoms()[x].xyz for x in self.moving_ref_atoms_iseqs])
      self.resulting_rmsd = rmsd_good
      # print "n_iter, rmsd:", n_iter, rmsd_good,
      # print get_main_chain_rmsd(moving_h, original_h)

      # self.states.add(sites_cart=self.moving_h.atoms().extract_xyz())
      # if n_iter % 100 == 0:
      #   moving_h.write_pdb_file(file_name="int_%d.pdb" % n_iter)
      self.n_iter += 1
      self.early_exit = previous_rmsd - rmsd_good < self.convergence_diff
      previous_rmsd = rmsd_good
    # print "number of iterations:", n_iter
    # print_rama_stats(phi_psi_atoms, r)
    # moving_h.write_pdb_file(file_name="int_%d.pdb" % n_iter)
    # states.write(file_name="all_states.pdb")

    # return rmsd_good, states, n_iter
