from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
from scitbx.math import superpose
from scitbx.array_family import flex
from scitbx import matrix
import mmtbx.ncs.ncs_utils as nu
import scitbx.rigid_body
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
import iotbx.cif.model
from six.moves import zip
import sys

class NCS_copy():
  def __init__(self,copy_iselection, rot, tran, str_selection=None, rmsd=999):
    """
    used for NCS groups list copies

    Attributes:
      iselection (flex.size_t): NCS copy selection
      r (matrix obj): rotation matrix from master to this copy
      t (matrix obj): translation vector from master to this copy
    """
    self.iselection = copy_iselection
    self.str_selection = str_selection
    self.r = rot
    self.t = tran
    self.rmsd = rmsd

  def __eq__(self, other):
    return approx_equal(self.r, other.r, out=null_out()) and approx_equal(self.t, other.t, out=null_out())

  def deep_copy(self):
    res = NCS_copy(
        self.iselection.deep_copy(),
        matrix.sqr(self.r),
        matrix.col(self.t),
        self.str_selection,
        self.rmsd)
    return res

  def select(self, selection):
    self.iselection = iselection_select(self.iselection, selection)
    self.str_selection = None # it is not valid anymore
    self.rmsd = 999

def iselection_select(isel, sel):
  x = flex.bool(sel.size(), False)
  x.set_selected(isel, True)
  res = x.select(sel).iselection()
  return res

class NCS_restraint_group(object):

  def __init__(self,master_iselection, str_selection=None):
    """
    used for NCS groups list

    Attributes:
      master_iselection (flex.size_t): NCS group master copy selection
      copies (list): list of NCS_copy objects
    """
    self.master_iselection = master_iselection
    self.master_str_selection = str_selection
    self.copies = []

  def __eq__(self, other):
    result = True
    for sc, oc in zip(self.copies, other.copies):
      result &= sc == oc
    return result

  def setup_selection_set(self):
    self.set_master_iselection = set(self.master_iselection)
    self.list_master_iselection = list(self.master_iselection)

  def get_iselections_list(self):
    """
    Returns all iselections in the group in one list
    """
    return [self.master_iselection] + [c.iselection for c in self.copies]

  def update_i_seqs(self, old_i_seqs):
    """
    correct iseqs using supplied list
    """
    if old_i_seqs is None:
      return
    for n,i in enumerate(self.master_iselection):
      self.master_iselection[n] = old_i_seqs[i]
    for c in self.copies:
      for n, i in enumerate(c.iselection):
        c.iselection[n] = old_i_seqs[i]

  def split_by_chains(self, hierarchy):
    # actually splitting. Looking for chains only in master. If some corner
    # case arise when it is not enough, will have to do something.
    # For example, chains in master and copy do not match with each other by
    # number of atoms. E.g. atoms belong to these chains:
    # master: AAAAABBBB
    #   copy: CCCDDDDDD
    #
    # Note, that I don't recalculate rotation/translation matrices!
    #
    # Looks like there's no tests for this functionality
    #
    from mmtbx.ncs.ncs_search import get_bool_selection_to_keep

    if len(hierarchy.select(self.master_iselection).only_model().chains()) == 1:
      return [self.deep_copy()]
    result = []
    selected_h = hierarchy.select(self.master_iselection)
    assert self.master_iselection.size() == selected_h.atoms_size()
    # shortcut (same chain ids)
    c_ids = [c.id for c in selected_h.only_model().chains()]
    if len(set(c_ids)) == 1:
      # assuming the configuration of copies is the same
      return [self.deep_copy()]

    for chain in selected_h.only_model().chains():
      c_iseqs = chain.atoms().extract_i_seq()
      to_keep = get_bool_selection_to_keep(
          big_selection=self.master_iselection,
          small_selection=c_iseqs)
      new_group = NCS_restraint_group(
          master_iselection = self.master_iselection.select(to_keep),
          str_selection = None)
      for old_copy in self.copies:
        new_copy = NCS_copy(
            copy_iselection=old_copy.iselection.select(to_keep),
            rot=matrix.sqr(old_copy.r),
            tran=matrix.col(old_copy.t),
            str_selection=None)
        new_group.append_copy(new_copy)
      result.append(new_group)
    return result

  def append_copy(self, copy):
    self.copies.append(copy)

  def get_number_of_copies(self):
    return len(self.copies)

  def deep_copy(self):
    result = NCS_restraint_group(
        master_iselection=self.master_iselection.deep_copy(),
        str_selection=self.master_str_selection)
    for c in self.copies:
      result.copies.append(c.deep_copy())
    return result

  def select(self, selection):
    """
    Modifies the selections of master and copies according the "selection"
    - Keep the order of selected atoms
    - Keep only atoms that appear in master and ALL copies
    Also modify "selection" to include ncs related atoms only if selected in
    both master and ALL ncs copies (The modified selection is not returned in
    current version)

    Args:
      selection (flex.bool): atom selection
    """
    from mmtbx.ncs.ncs_utils import selected_positions, remove_items_from_selection
    assert isinstance(selection, flex.bool)


    iselection = selection.iselection(True)
    sel_set = set(iselection)
    m = set(self.master_iselection)
    m_list = [(pos,indx) for pos,indx in enumerate(list(self.master_iselection))]
    m_in_sel = m.intersection(sel_set)
    common_selection_pos = {pos for (pos,indx) in m_list if indx in m_in_sel}
    for ncs in self.copies:
      c = set(ncs.iselection)
      c_list = [(pos,indx) for pos,indx in enumerate(list(ncs.iselection))]
      copy_in_sel = c.intersection(sel_set)
      include_set = {pos for (pos,indx) in c_list if indx in copy_in_sel}
      common_selection_pos.intersection_update(include_set)
      if not bool(common_selection_pos): break
    # use the common_selection_pos to update all selections
    self.master_iselection, not_included = selected_positions(
      self.master_iselection,common_selection_pos)
    iselection = remove_items_from_selection(iselection,not_included)
    for ncs in self.copies:
      ncs.iselection, not_included = selected_positions(
        ncs.iselection,common_selection_pos)
      iselection = remove_items_from_selection(iselection,not_included)
    for c in self.copies:
      assert self.master_iselection.size() == c.iselection.size(), "%s\n%s" % (
          list(self.master_iselection), list(c.iselection))

    # This is to handle renumbering properly
    self.master_iselection = iselection_select(self.master_iselection, selection)
    self.master_str_selection = None
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
    t_sel = self.copies[n].str_selection
    self.copies[n].str_selection = self.master_str_selection
    self.master_str_selection = t_sel
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
    for g in self:
      assert isinstance(g, NCS_restraint_group)

  def __eq__(self, other):
    result = self.get_n_groups() == other.get_n_groups()
    for sg, og in zip(self, other):
      result &= (sg == og)
    return result

  def setup_sets(self):
    for g in self:
      g.setup_selection_set()

  def get_copy_iseqs(self, iseqs):
    """get iseqs from copies for proxy. E.g. for bond:
    iseqs = [1,2]

    Args:
        iseqs (iterable): iseqs of original proxy

    Returns:
        [[3,4], [5,6], [7,8]]
    """
    result = []
    # self.setup_sets()
    # print("iseqs in get_copy_iseqs:", iseqs)
    for gr in self:
      if iseqs[0] in gr.set_master_iselection:
        # check the rest are in:
        for iseq in iseqs[1:]:
          assert iseq in gr.set_master_iselection
        # now iterate over input iseqs and populate the result
        for i in range(gr.get_number_of_copies()):
          result.append([])
        for in_iseq in iseqs:
          # find the index:
          iseq_idex = gr.list_master_iselection.index(in_iseq)
          for i, c in enumerate(gr.copies):
            result[i].append(c.iselection[iseq_idex])
        return result
    return result

  def get_n_groups(self):
    return len(self)

  def update_str_selections_if_needed(
      self, hierarchy, asc=None, chains_info=None):
    from mmtbx.ncs.ncs_search import get_chains_info
    from iotbx.pdb.atom_selection import selection_string_from_selection
    if asc is None:
      asc = hierarchy.atom_selection_cache()
    if chains_info is None:
      chains_info = get_chains_info(hierarchy)
    for gr in self:
      if gr.master_str_selection is None:
        gr.master_str_selection = selection_string_from_selection(
            hierarchy,
            gr.master_iselection,
            chains_info,
            asc)
        for c in gr.copies:
          if c.str_selection is None:
            c.str_selection = selection_string_from_selection(
                hierarchy,
                c.iselection,
                chains_info,
                asc)

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

  def split_by_chains(self, hierarchy):
    new_groups = class_ncs_restraints_group_list()
    for g in self:
      new_gs = g.split_by_chains(hierarchy)
      new_groups += new_gs
    return new_groups

  def filter_out_small_groups(self, min_n_atoms=3):
    new_groups = class_ncs_restraints_group_list()
    for g in self:
      if g.master_iselection.size() >= min_n_atoms:
        new_groups.append(g.deep_copy())
    return new_groups

  def update_i_seqs(self, old_i_seqs):
    """
    correct iseqs using supplied list
    """
    for group in self:
      group.update_i_seqs(old_i_seqs)

  def filter_ncs_restraints_group_list(self, whole_h, ncs_obj):
    """ Remove ncs groups where master or copy does not cover whole chain
    (some atoms are left behind).
    Reason for this - when big moves are likely (e.g. in real-space refine or
    model idealization), the chain can get a big gap in place where NCS ends.
    This leads to undesired artefacts in refinement.
    """
    def whole_chain_in_ncs(whole_h, master_iselection):
      m_c_id = whole_h.atoms()[master_iselection[0]].parent().parent().parent().id
      for chain in ncs_obj.truncated_hierarchy.only_model().chains():
        if chain.id == m_c_id:
          if chain.atoms_size() <= master_iselection.size():
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
        cp.rmsd = asu_site_cart.select(c_sel).rms_difference(lsq_fit_obj.other_sites_best_fit())

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
          print('Allowed rmsd : {}, rmsd: {}'.format(chain_max_rmsd,rmsd), file=log)
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

  def _show(self, hierarchy=None, brief=True, out=sys.stdout):
    """
    For debugging
    """
    #print("debugging output of ncs_restraints_group_list")
    for group in self:
      print("Master str selection:", group.master_str_selection, file=out)
      if not brief:
        print(list(group.master_iselection), file=out)
      if hierarchy is not None:
        print(hierarchy.select(group.master_iselection).as_pdb_string(), file=out) # PDB OK - debugging output
      for c in group.copies:
        print("Copy str selection:", c.str_selection, file=out)
        if not brief:
          print(list(c.iselection), file=out)
        # print "rot", list(c.r)
        # print "tran", list(c.t)
        if hierarchy is not None:
          print(hierarchy.select(c.iselection).as_pdb_string(), file=out) # PDB OK - debugging output
      print("="*30, file=out)
    #print("end debugging output of ncs_restraints_group_list")


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

  def get_ncs_groups_shifts(self, sites_cart, current):
    shifts = []
    for nrg in self:
      master_ncs_selection = nrg.master_iselection
      master_xyz = sites_cart.select(master_ncs_selection)
      mu_m = flex.vec3_double([master_xyz.mean()])
      shifts.append(mu_m-current.mean())
      for ncs_copy in nrg.copies:
        asu_selection = ncs_copy.iselection
        asu_xyz = sites_cart.select(asu_selection)
        mu_i = flex.vec3_double([asu_xyz.mean()])
        shifts.append(mu_i-current.mean())
    return shifts

  def get_all_copies_selection(self):
    result = flex.size_t()
    for nrg in self:
      for c in nrg.copies:
        result.extend(c.iselection)
    return flex.sorted(result)

  def get_extended_ncs_selection(self, refine_selection):
    """
    Args:
      refine_selection (flex.size_t): of all ncs related copies and
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

  def get_array_of_str_selections(self):
    """
    Returns array of phil selection strings e.g. for the exapmle above in
    print_ncs_phil_param:
    [['(Chain A)','(chain C)','(chain E)'],['(chain B)','(chain D)','(chain F)']]
    """
    result = []
    for gr in self:
      group = [gr.master_str_selection]
      for c in gr.copies:
        group.append(c.str_selection)
      result.append(group)
    return result

  def unique_with_biomt(self, hierarchy):
    # first we need to check if it is possible.
    # Criteria:
    # - all model is covered by NCS
    # - every chain is fully covered by NCS
    # self.filter_ncs_restraints_group_list
    # also limiting to 1 NCS group. Not clear how to describe multiple groups
    # where operations need to performed on a different chains

    assert len(self) == 1

    cif_block = iotbx.cif.model.block()


    resulting_hierarchy = hierarchy.select(self[0].master_iselection)
    master_label_asym_ids = []
    for c in resulting_hierarchy.only_model().chains():
      lai = resulting_hierarchy.get_label_asym_id(c.residue_groups()[0])
      master_label_asym_ids.append(lai)

    pdbx_struct_assembly_gen_loop = iotbx.cif.model.loop(header=(
        '_pdbx_struct_assembly_gen.assembly_id',
        '_pdbx_struct_assembly_gen.oper_expression',
        '_pdbx_struct_assembly_gen.asym_id_list',))

    pdbx_struct_assembly_loop = iotbx.cif.model.loop(header=(
        '_pdbx_struct_assembly.id',
        '_pdbx_struct_assembly.details',
        '_pdbx_struct_assembly.method_details',
        '_pdbx_struct_assembly.oligomeric_details',
        '_pdbx_struct_assembly.oligomeric_count',))

    pdbx_struct_oper_list_loop = iotbx.cif.model.loop(header=(
       '_pdbx_struct_oper_list.id',
       '_pdbx_struct_oper_list.type',
       '_pdbx_struct_oper_list.name',
       '_pdbx_struct_oper_list.symmetry_operation',
       '_pdbx_struct_oper_list.matrix[1][1]',
       '_pdbx_struct_oper_list.matrix[1][2]',
       '_pdbx_struct_oper_list.matrix[1][3]',
       '_pdbx_struct_oper_list.matrix[2][1]',
       '_pdbx_struct_oper_list.matrix[2][2]',
       '_pdbx_struct_oper_list.matrix[2][3]',
       '_pdbx_struct_oper_list.matrix[3][1]',
       '_pdbx_struct_oper_list.matrix[3][2]',
       '_pdbx_struct_oper_list.matrix[3][3]',
       '_pdbx_struct_oper_list.vector[1]',
       '_pdbx_struct_oper_list.vector[2]',
       '_pdbx_struct_oper_list.vector[3]',))

    master_asym_id = ','.join(master_label_asym_ids)
    pdbx_struct_assembly_gen_loop.add_row({
        '_pdbx_struct_assembly_gen.assembly_id':'1',
        '_pdbx_struct_assembly_gen.oper_expression':'(1-%d)' % (len(self[0].copies)+1),
        '_pdbx_struct_assembly_gen.asym_id_list': master_asym_id,
      })

    pdbx_struct_assembly_loop.add_row({
        '_pdbx_struct_assembly.id': '1',
        '_pdbx_struct_assembly.details': 'Symmetry assembly',
        '_pdbx_struct_assembly.method_details': '?',
        '_pdbx_struct_assembly.oligomeric_details': '?',
        '_pdbx_struct_assembly.oligomeric_count': '?',
      })

    # put in identity transform
    oper_id = 1
    row = [oper_id, 'point symmetry operation', '?', '?']
    row.extend([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
    pdbx_struct_oper_list_loop.add_row(row)

    for i_ncs_copy, ncs_copy in enumerate(self[0].copies):
      oper_id = i_ncs_copy + 2
      row = [oper_id, 'point symmetry operation', '?', '?']
      row.extend(ncs_copy.r)
      row.extend(ncs_copy.t)
      pdbx_struct_oper_list_loop.add_row(row)

    cif_block.add_loop(pdbx_struct_assembly_gen_loop)
    cif_block.add_loop(pdbx_struct_assembly_loop)
    cif_block.add_loop(pdbx_struct_oper_list_loop)

    return self[0].master_iselection, cif_block


  def as_cif_block(self, cif_block, hierarchy, scattering_type, ncs_type):
    """
    Let me lay out what I found and please correct me if I am wrong in any detail.
    There are some questions along as well. Then I will implement the correct NCS
    output in Phenix.
    We look at NCS groups in the following way:
    ncs_group {
      reference = chain 'C'
      selection = chain 'E'
      selection = chain 'G'
      selection = chain 'I'
    }
    ncs_group {
      reference = chain 'D'
      selection = chain 'F'
      selection = chain 'H'
      selection = chain 'J'
    }
    In this example we have 2 NCS groups with 4 chains in each. For simplicity,
    let's assume that whole chains are included. This means that chains C,E,G,I
    are NCS-related and (almost) identical. The same is for chains D,F,H,J.
    Chains in different NCS groups does not match each other. We also know relation
    (rotation/translation matrices) between chains C<--E, C<--G and C<--I.
    Same goes to the second group.

    Now I will try to translate this into mmCIF terminology:
    Ensemble - each of ncs_group.
    Domain - each of the reference/selection chains.
    If so, then mmCIF should have loops as following:
     _struct_ncs_ens.id
     _struct_ncs_ens.details
    en1 ?
    en2 ?
    We almost never know any useful details, so I'm putting ? here.
    Or will it be better to put something else very generic?
    Then we list 'domains' (we will put valid Phenix atom selection syntax in details):
    loop_
        _struct_ncs_dom.id
        _struct_ncs_dom.pdbx_ens_id
        _struct_ncs_dom.details
         d1 en1  'Chains C'
         d2 en1  'Chains E'
         d3 en1  'Chains G'
         d4 en1  'Chains I'
         d1 en2  'Chains D'
         d2 en2  'Chains F'
         d3 en2  'Chains H'
         d4 en2  'Chains J'

    And to relate everything to matrices:
    _struct_ncs_ens_gen.dom_id_1
    _struct_ncs_ens_gen.dom_id_2
    _struct_ncs_ens_gen.ens_id
    _struct_ncs_ens_gen.oper_id
    d1 d2 en1 op1
    d1 d3 en1 op2
    d1 d4 en1 op3
    d1 d2 en2 op4
    d1 d3 en2 op5
    d1 d4 en2 op6

    Question: Here is an important question: should we output matrices defining
    how to align dom_id_2 onto dom_id_1 or another way around?

    Answer (John Berrisford): "NCS operators are defined by struct_ncs_ens_gen.
    This defines which domain moves onto the other. Dom_id_1 doesn't move,
    dom_id_2 is transformed by the operator."
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ncs_ens_gen.html


    Content of _struct_ncs_oper
    (matrices themselves) seems trivial, so omitting from this example.

    Now I'm guessing _struct_ncs_dom_lim was designed with the ability to
    provide more than one residue-level interval for each 'domain' and
    pdbx_component_id is for numbering these intervals. In our simple case there
    will be only one row for each domain:
    loop_
        _struct_ncs_dom_lim.dom_id
        _struct_ncs_dom_lim.pdbx_ens_id
        _struct_ncs_dom_lim.pdbx_component_id
        _struct_ncs_dom_lim.beg_label_alt_id
        _struct_ncs_dom_lim.beg_label_asym_id
        _struct_ncs_dom_lim.beg_label_comp_id
        _struct_ncs_dom_lim.beg_label_seq_id
        _struct_ncs_dom_lim.end_label_alt_id
        _struct_ncs_dom_lim.end_label_asym_id
        _struct_ncs_dom_lim.end_label_comp_id
        _struct_ncs_dom_lim.end_label_seq_id
         d1 en1 1 .  C PRO  1  . C GLY  29
         d2 en1 1 .  E PRO  1  . E GLY  29
         d3 en1 1 .  G PRO  1  . G GLY  29
         d4 en1 1 .  I PRO  1  . I GLY  29
         d1 en2 1 .  D PRO  31 . D GLY  59
         d2 en2 1 .  F PRO  31 . F GLY  59
         d3 en2 1 .  H PRO  31 . H GLY  59
         d4 en2 1 .  J PRO  31 . J GLY  59
    Can we omit _struct_ncs_dom_lim loop? We will put valid Phenix atom selection
    syntax in selection_details in this loop.

    JB: "Please do not omit _struct_ncs_dom_lim - this provides
    the exact details of which residues from which chain have to be transformed
    with the NCS operator.  This is the bit read by users such as pdb_redo when
    they are re-refining the structures. Note in struct_ncs_dom_lim there are
    fields for auth_asym_id (chain ID in PDB speak) and
    label_asym_id (mmCIF only chain identifiers). Please do not mix these two."

    Then, after refinement, we can construct _refine_ls_restr_ncs where
    pdbx_ordinal - is just ordinal number in this loop
    dom_id - domain id (== _struct_ncs_ens_gen.dom_id_1 == _struct_ncs_dom.id)
    pdbx_ens_id - ensemble id (== _struct_ncs_ens_gen.ens_id == _struct_ncs_ens.id )
    pdbx_asym_id - chain id from the first component (_struct_ncs_dom_lim.pdbx_component_id)?
    There is an example of domain with 3 components:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ncs_dom_lim.html
    and chain A is chosen for this asym_id:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/refine_ls_restr_ncs.html
    Do I understand correctly that you expect information on only 3 domains from
    each ensemble (marked 'selection' in Phenix notation) and it should show
    how this particular domain relates to the first one ('reference' in Phenix notation)?

    JB: "yes, refine_ls_restr_ncs describes the refinement result and gives the
    details of moving domain 2 onto domain 1 - i.e. how well do they overlay if
    the transform is applied. So domain 1 doesn't need to be described here - it
    would fit perfectly to itself (we hope!)"

    _struct_ncs_oper.code is 'given' since in Phenix we never omit NCS-related
    parts of model.
    """
    #
    def _consecutive_ranges(isel):
      '''
      Split selection into consecutive ranges that can be represented
      by ncs_dom_lim. Maybe useful elsewhere.
      '''
      result = []
      start_index = 0
      cur_index = 1
      atoms = hierarchy.atoms()
      # shortcut: if all atoms form consecutive range and in one chain:
      if ( (isel[-1] -isel[0] == len(isel)-1)
          and (hierarchy.get_label_asym_id_iseq(isel[0]) == \
                hierarchy.get_label_asym_id_iseq(isel[-1]))):
        return [isel]

      seq_id = hierarchy.get_label_seq_id_iseq(isel[0])
      seq_id = int(seq_id) if seq_id !='.' else 0
      asym_id = hierarchy.get_label_asym_id_iseq(isel[0])
      while cur_index < len(isel):
        if cur_index >= len(isel):
          break
        seq_id_1 = seq_id
        asym_id_1 = asym_id
        #
        seq_id = hierarchy.get_label_seq_id_iseq(isel[cur_index])
        seq_id = int(seq_id) if seq_id !='.' else 0
        asym_id = hierarchy.get_label_asym_id_iseq(isel[cur_index])
        # print ('  DEBUG: seq_id, asym_id,', seq_id, asym_id, seq_id_1, asym_id_1)
        if (
            # consecutive indices and no chain break
            (isel[cur_index]-isel[cur_index-1] == 1
                and asym_id == asym_id_1)
            # same chain, same or next resid, no indices check (missing atom, etc)
            or (asym_id == asym_id_1
                and abs(seq_id - seq_id_1) <= 1 )
            ):
          cur_index += 1
        else: # end range
          # assert start_index != cur_index-1 # maybe needed for 1 atom chains, like ions
          result.append(isel[start_index:cur_index])
          start_index = cur_index
          cur_index += 1
      if start_index != cur_index-1:
        result.append(isel[start_index:cur_index])
      # print ('Debug result', [list(x) for x in result])
      return result

    def _get_struct_ncs_dom_lim_row(dom_id, comp_id, r):
      return {
          "_struct_ncs_dom_lim.pdbx_ens_id": struct_ncs_ens_id,
          "_struct_ncs_dom_lim.dom_id": dom_id,
          "_struct_ncs_dom_lim.pdbx_component_id": comp_id,
          "_struct_ncs_dom_lim.beg_label_alt_id": hierarchy.get_label_alt_id_iseq(r[0]),
          "_struct_ncs_dom_lim.beg_label_asym_id": hierarchy.get_label_asym_id_iseq(r[0]),
          "_struct_ncs_dom_lim.beg_label_comp_id": hierarchy.atoms()[r[0]].parent().resname.strip(),
          "_struct_ncs_dom_lim.beg_label_seq_id": hierarchy.get_label_seq_id_iseq(r[0]),
          "_struct_ncs_dom_lim.end_label_alt_id": hierarchy.get_label_alt_id_iseq(r[-1]),
          "_struct_ncs_dom_lim.end_label_asym_id": hierarchy.get_label_asym_id_iseq(r[-1]),
          "_struct_ncs_dom_lim.end_label_comp_id": hierarchy.atoms()[r[-1]].parent().resname.strip(),
          "_struct_ncs_dom_lim.end_label_seq_id": hierarchy.get_label_seq_id_iseq(r[-1])}

    ncs_ens_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens.id",
      "_struct_ncs_ens.details"))
    ncs_dom_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom.pdbx_ens_id",
      "_struct_ncs_dom.id",
      "_struct_ncs_dom.details"))
    ncs_dom_lim_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom_lim.pdbx_ens_id",
      "_struct_ncs_dom_lim.dom_id",
      "_struct_ncs_dom_lim.pdbx_component_id",
      #"_struct_ncs_dom_lim.pdbx_refine_code", # no need, from CCP4
      "_struct_ncs_dom_lim.beg_label_alt_id",
      "_struct_ncs_dom_lim.beg_label_asym_id",
      "_struct_ncs_dom_lim.beg_label_comp_id",
      "_struct_ncs_dom_lim.beg_label_seq_id",
      "_struct_ncs_dom_lim.end_label_alt_id",
      "_struct_ncs_dom_lim.end_label_asym_id",
      "_struct_ncs_dom_lim.end_label_comp_id",
      "_struct_ncs_dom_lim.end_label_seq_id",))

    ncs_oper_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_oper.id",
      "_struct_ncs_oper.code",
      "_struct_ncs_oper.matrix[1][1]",
      "_struct_ncs_oper.matrix[1][2]",
      "_struct_ncs_oper.matrix[1][3]",
      "_struct_ncs_oper.matrix[2][1]",
      "_struct_ncs_oper.matrix[2][2]",
      "_struct_ncs_oper.matrix[2][3]",
      "_struct_ncs_oper.matrix[3][1]",
      "_struct_ncs_oper.matrix[3][2]",
      "_struct_ncs_oper.matrix[3][3]",
      "_struct_ncs_oper.vector[1]",
      "_struct_ncs_oper.vector[2]",
      "_struct_ncs_oper.vector[3]",
      "_struct_ncs_oper.details"))

    ncs_ens_gen_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens_gen.ens_id",
      "_struct_ncs_ens_gen.dom_id_1",
      "_struct_ncs_ens_gen.dom_id_2",
      "_struct_ncs_ens_gen.oper_id"))

    refine_ls_restr_ncs_loop = iotbx.cif.model.loop(header=(
      "_refine_ls_restr_ncs.pdbx_ordinal",
      "_refine_ls_restr_ncs.pdbx_ens_id",
      "_refine_ls_restr_ncs.dom_id",
      "_refine_ls_restr_ncs.pdbx_refine_id",
      "_refine_ls_restr_ncs.pdbx_asym_id",
      "_refine_ls_restr_ncs.pdbx_type",
      "_refine_ls_restr_ncs.weight_position",
      "_refine_ls_restr_ncs.weight_B_iso",
      "_refine_ls_restr_ncs.rms_dev_position",
      "_refine_ls_restr_ncs.rms_dev_B_iso",
      "_refine_ls_restr_ncs.ncs_model_details"))

    self.update_str_selections_if_needed(hierarchy)
    self.recalculate_ncs_transforms(hierarchy.atoms().extract_xyz())
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    if len(self) == 0:
      return cif_block
    refine_ls_restr_ncs_pdbx_ordinal = 1
    n_oper = 1
    for i_group, group in enumerate(self):
      struct_ncs_ens_id = 'ens_%d' % (i_group+1)
      ncs_ens_loop.add_row({"_struct_ncs_ens.id" : struct_ncs_ens_id})
      # master loops
      master_dom_id = 'd_1'
      ncs_dom_loop.add_row({
          "_struct_ncs_dom.pdbx_ens_id":struct_ncs_ens_id,
          "_struct_ncs_dom.id":master_dom_id,
          "_struct_ncs_dom.details":"%s" % group.master_str_selection.replace("'", '"')})
      ranges = _consecutive_ranges(group.master_iselection)
      for i_r, r in enumerate(ranges):
        ncs_dom_lim_loop.add_row(_get_struct_ncs_dom_lim_row(master_dom_id, i_r+1, r))
      # now get to copies
      for i_ncs_copy, ncs_copy in enumerate(group.copies):
        copy_dom_id = 'd_%d' % (i_ncs_copy+2) # no 0, 1-master
        ncs_dom_loop.add_row({
            "_struct_ncs_dom.pdbx_ens_id":struct_ncs_ens_id,
            "_struct_ncs_dom.id":copy_dom_id,
            "_struct_ncs_dom.details":'%s' % ncs_copy.str_selection.replace("'", '"')})

        ranges = _consecutive_ranges(ncs_copy.iselection)
        for i_r, r in enumerate(ranges):
          ncs_dom_lim_loop.add_row(_get_struct_ncs_dom_lim_row(copy_dom_id, i_r+1, r))
        oper_id = 'op_%d' % n_oper
        n_oper += 1
        ncs_ens_gen_loop.add_row({
            "_struct_ncs_ens_gen.dom_id_1":copy_dom_id,
            "_struct_ncs_ens_gen.dom_id_2":master_dom_id,
            "_struct_ncs_ens_gen.ens_id":struct_ncs_ens_id,
            "_struct_ncs_ens_gen.oper_id":oper_id})
        row = [oper_id, 'given']
        row.extend(ncs_copy.r)
        row.extend(ncs_copy.t)
        row.append('?')
        ncs_oper_loop.add_row(row)

        refine_ls_restr_ncs_loop.add_row({
            "_refine_ls_restr_ncs.pdbx_ordinal": refine_ls_restr_ncs_pdbx_ordinal,
            "_refine_ls_restr_ncs.dom_id": copy_dom_id,
            "_refine_ls_restr_ncs.pdbx_refine_id": scattering_type,
            "_refine_ls_restr_ncs.pdbx_ens_id": struct_ncs_ens_id,
            "_refine_ls_restr_ncs.pdbx_asym_id": "'%s'" % hierarchy.get_label_asym_id_iseq(group.master_iselection[0]),
            "_refine_ls_restr_ncs.pdbx_type": ncs_type,
            "_refine_ls_restr_ncs.weight_position": '?', # weight_position
            "_refine_ls_restr_ncs.weight_B_iso": '?', # weight_B_iso
            "_refine_ls_restr_ncs.rms_dev_position": ncs_copy.rmsd,
            "_refine_ls_restr_ncs.rms_dev_B_iso": '?', # rms_dev_B_iso
            "_refine_ls_restr_ncs.ncs_model_details": '?', # model_details
            })
        refine_ls_restr_ncs_pdbx_ordinal += 1
      cif_block.add_loop(ncs_ens_loop)
      cif_block.add_loop(ncs_dom_loop)
      cif_block.add_loop(ncs_dom_lim_loop)
      cif_block.add_loop(ncs_oper_loop)
      cif_block.add_loop(ncs_ens_gen_loop)
      cif_block.add_loop(refine_ls_restr_ncs_loop)
    return cif_block
