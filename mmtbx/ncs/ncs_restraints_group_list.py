from __future__ import division
from StringIO import StringIO
from scitbx.math import superpose
from scitbx.array_family import flex
from scitbx import matrix
import mmtbx.ncs.ncs_utils as nu
import scitbx.rigid_body


class NCS_copy():
  def __init__(self,copy_iselection, rot, tran):
    """
    used for NCS groups list copies

    Attributes:
      iselection (flex.size_t): NCS copy selection
      r (matrix obj): rotation matrix from master to this copy
      t (matrix obj): translation vector from master to this copy
    """
    self.iselection = copy_iselection
    self.r = rot
    self.t = tran

  def deep_copy(self):
    res = NCS_copy(self.iselection.deep_copy(), self.r, self.t)
    return res

  def select(self, selection):
    self.iselection = iselection_select(self.iselection, selection)

def iselection_select(isel, sel):
  x = flex.bool(sel.size(), False)
  x.set_selected(isel, True)
  res = x.select(sel).iselection()
  return res

class NCS_restraint_group(object):

  def __init__(self,master_iselection):
    """
    used for NCS groups list

    Attributes:
      master_iselection (flex.size_t): NCS group master copy selection
      copies (list): list of NCS_copy objects
    """
    self.master_iselection = master_iselection
    self.copies = []

  def get_iselections_list(self):
    """
    Returns all iselections in the group in one list
    """
    return [self.master_iselection] + [c.iselection for c in self.copies]

  def get_number_of_copies(self):
    return len(self.copies)

  def deep_copy(self):
    result = NCS_restraint_group(self.master_iselection.deep_copy())
    for c in self.copies:
      result.copies.append(c.deep_copy())
    return result

  def select(self, selection):
    assert isinstance(selection, flex.bool)
    self.master_iselection = iselection_select(self.master_iselection, selection)
    for c in self.copies:
      c.select(selection)

  def whole_group_iselection(self):
    isel = self.master_iselection.deep_copy()
    for cp in self.copies:
      isel.extend(cp.iselection)
    # make sure sequential order of selection indices
    return flex.sorted(isel)

  def make_nth_copy_master(self, n):
    """
    n - index of the copy to become master, starting with 0
    master becomes nth copy
    """
    # switch master and copy selection
    assert n < self.get_number_of_copies()
    self.master_iselection, self.copies[n].iselection = \
      self.copies[n].iselection, self.master_iselection
    # Adjust rotation and translation for the new master
    r = self.copies[n].r = (self.copies[n].r.transpose())
    t = self.copies[n].t = -(self.copies[n].r * self.copies[n].t)
    # change all other rotations and translations to the new master
    for i, c in enumerate(self.copies):
      if i == n: continue
      # change translation before rotation
      c.t = (c.r * t + c.t)
      c.r = (c.r * r)

class class_ncs_restraints_group_list(list):
  def __init__(self, *args):
    super(class_ncs_restraints_group_list, self).__init__(*args)

  def get_n_groups(self):
    return len(self)

  def deep_copy(self):
    result = class_ncs_restraints_group_list()
    for gr in self:
      result.append(gr.deep_copy())
    return result

  def select(self, selection):
    assert isinstance(selection, flex.bool)
    result = self.deep_copy()
    for gr in result:
      gr.select(selection)
    return result

  def filter_ncs_restraints_group_list(self, whole_h):
    """ Remove ncs groups where master or copy does not cover whole chain
    (some atoms are left behind).
    Reason for this - when big moves are likely (e.g. in real-space refine or
    model idealization), the chain can get a big gap in place where NCS ends.
    This leads to undesired artefacts in refinement.
    """
    def whole_chain_in_ncs(whole_h, master_iselection):
      m_c = whole_h.select(master_iselection)
      m_c_id = m_c.only_model().chains()[0].id
      for chain in whole_h.only_model().chains():
        if chain.id == m_c_id:
          n_non_h_atoms = 0
          for a in chain.atoms():
            # print "'%s'" % a.element
            if not a.element_is_hydrogen():
              n_non_h_atoms += 1
          # print "n_non_h_atoms, master_iselection.size()", n_non_h_atoms, master_iselection.size()
          if n_non_h_atoms <= master_iselection.size():
            return True
          else:
            return False
    n_gr_to_remove = []
    for i, ncs_gr in enumerate(self):
      if not whole_chain_in_ncs(whole_h, ncs_gr.master_iselection):
        n_gr_to_remove.append(i)
        continue
      for c in ncs_gr.copies:
        if not whole_chain_in_ncs(whole_h, c.iselection):
          n_gr_to_remove.append(i)
          break
    result = self.deep_copy()
    for i in reversed(n_gr_to_remove):
      del result[i]
    return result

  def recalculate_ncs_transforms(self, asu_site_cart):
    """
    Re-evaluate the rotation and translation in the ncs groups list, base on
    the ncs groups selection and the atoms location.
    Updates self.

    Args:
      asu_site_cart (flex.vec_3): the complete ASU sites cart (coordinates)
    """
    for gr in self:
      m_sel = gr.master_iselection
      for cp in gr.copies:
        c_sel = cp.iselection
        # other_sites are the master, reference_sites are the copies
        lsq_fit_obj = superpose.least_squares_fit(
            reference_sites = asu_site_cart.select(c_sel),
            other_sites     = asu_site_cart.select(m_sel))
        cp.r = lsq_fit_obj.r
        cp.t = lsq_fit_obj.t

  def check_for_max_rmsd(self,
      sites_cart,
      chain_max_rmsd=10.0,
      log=StringIO()):
    """
    Check that all copies relate correctly to master via the transforms

    Args:
      sites_cart: Atom coordinates
      chain_max_rmsd (float): maximum allowed rmsd between coordinates copies

    Returns:
      nrgl_ok (bool): True when ncs_restraints_group_list is OK
    """
    nrgl_ok = True
    for i,gr in enumerate(self):
      master_xyz = sites_cart.select(gr.master_iselection)
      for j,cp in enumerate(gr.copies):
        copy_xyz = sites_cart.select(cp.iselection)
        xyz = cp.r.elems * master_xyz + cp.t
        rmsd = copy_xyz.rms_difference(xyz)
        nrgl_ok &= (rmsd <= chain_max_rmsd)
        if (rmsd > chain_max_rmsd):
          print >>log,'Allowed rmsd : {}, rmsd: {}'.format(chain_max_rmsd,rmsd)
    return nrgl_ok

  def shift_translation_to_center(self, shifts):
    """
    Add shifts to the translation component of ncs_restraints_group_list
    towards the center of coordinates

    Args:
      shifts (list): [mu_1, mu_1, mu_2...] where the mu stands
        for the shift of the master copy to the coordinate center mu is (dx,dy,dz)

    Returns:
      new ncs_restraints_group_list
    """
    new_list = []
    if bool(shifts):
      new_list = self.deep_copy()
      i = 0
      for nrg in new_list:
        for ncs_copy in nrg.copies:
          mu = shifts[i]
          i += 1
          # Only the translation is changing
          t = ncs_copy.r.elems * mu + ncs_copy.t - mu
          ncs_copy.t = matrix.col(t[0])
    return new_list

  def shift_translation_back_to_place(self, shifts):
    """
    shifts to the translation component of ncs_restraints_group_list from the
    center of coordinates back to place

    Args:
      shifts (list): [mu_1, mu_1, mu_2...] where the mu stands
        for the shift of the master copy to the coordinate center mu is (dx,dy,dz)

    Returns:
      new ncs_restraints_group_list
    """
    if bool(shifts):
      i = 0
      new_list = self.deep_copy()
      for nrg in new_list:
        for ncs_copy in nrg.copies:
          mu = shifts[i]
          i += 1
          # Only the translation is changing
          t = mu - ncs_copy.r.elems * mu + ncs_copy.t
          ncs_copy.t = matrix.col(t[0])
    else:
      new_list = ncs_restraints_group_list
    return new_list

  def get_ncs_groups_centers(self, sites_cart):
    """
    calculate the center of coordinate for the master of each ncs copy

    Args:
      sites_cart

    Returns:
      shifts (list): [mu_1, mu_1, mu_2...] where the mu stands
      for the shift of the master copy to the coordinate center mu is (dx,dy,dz)
    """
    shifts = []
    for nrg in self:
      master_ncs_selection = nrg.master_iselection
      master_xyz = sites_cart.select(master_ncs_selection)
      mu_m = flex.vec3_double([master_xyz.mean()])
      # add a copy of the master coordinate center for each copy
      for ncs_copy in nrg.copies:
        shifts.append(mu_m)
    return shifts

  def get_extended_ncs_selection(self, refine_selection):
    """
    Args:
      refine_selection (flex.siz_t): of all ncs related copies and
        non ncs related parts to be included in selection (to be refined)

    Returns:
      (flex.siz_t): selection of all ncs groups master ncs selection and
        non ncs related portions that are being refined (exclude NCS copies)
    """
    if not refine_selection:
      refine_selection = []
    refine_selection = set(refine_selection)
    total_master_ncs_selection = set()
    total_ncs_related_selection = set()
    for nrg in self:
      master_ncs_selection = nrg.master_iselection
      total_master_ncs_selection.update(set(master_ncs_selection))
      for ncs_copy in nrg.copies:
        asu_selection = ncs_copy.iselection
        total_ncs_related_selection.update(set(asu_selection))
    if refine_selection:
      # make sure all ncs related parts are in refine_selection
      all_ncs = total_master_ncs_selection | total_ncs_related_selection
      not_all_ncs_related_atoms_selected = bool(all_ncs - refine_selection)
      if not_all_ncs_related_atoms_selected:
        msg = 'refine_selection does not contain all ncs related atoms'
        raise Sorry(msg)
      #
      extended_ncs_selection = refine_selection - total_ncs_related_selection
      return flex.size_t(list(extended_ncs_selection))
    else:
      # if refine_selection is None
      return flex.size_t(list(total_master_ncs_selection))

  def concatenate_rot_tran(self):
    """
    Concatenate rotation angles, corresponding to the rotation
    matrices and scaled translation vectors to a single long flex.double object

    Returns:
      flex.double : [(alpha_1,beta_1,gamma_1,Tx_1,Ty_1,Tz_1)...]
    """
    x = []
    for gr in self:
      for tr in gr.copies:
        x.extend(list(nu.rotation_to_angles(rotation=tr.r.elems))
                 + list(tr.t.elems))
    return flex.double(x)

  def update_rot_tran(self, x):
    """
    Convert the refinable parameters, rotations angles and
    scaled translations, back to rotation matrices and translation vectors and
    updates the transforms_obj (ncs_restraints_group_list)
    !!! IN PLACE !!!

    Args:
      x : a flex.double of the form (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
        theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s). where n is the number of
        transformations.

    Returns:
      Nothing
    """
    i = 0
    for gr in self:
      copies = []
      for tr in gr.copies:
        the,psi,phi =x[i*6:i*6+3]
        rot = scitbx.rigid_body.rb_mat_xyz(
          the=the, psi=psi, phi=phi, deg=False)
        tran = matrix.rec(x[i*6+3:i*6+6],(3,1))
        tr.r = (rot.rot_mat())
        tr.t = tran
        copies.append(tr)
        i += 1
      gr.copies = copies

  def get_rotation_translation_as_list(self):
    """
    XXX
    XXX Consider deletion. Used only in tests tst_minimization_ncs_constraints_real_space.py,
    XXX tst_ncs_utils.py
    XXX

    Get rotations and translations vectors from ncs_restraints_group_list or
    transforms_obj

    Returns:
      r (list): list of rotation matrices
      t (list): list of translation vectors
    """
    r = []
    t = []
    for nrg in self:
      for tr in nrg.copies:
        r.append(tr.r)
        t.append(tr.t)
    return r,t
