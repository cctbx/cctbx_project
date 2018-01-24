from __future__ import division
from scitbx.array_family import flex
from scitbx.math import superpose
from libtbx.utils import Sorry
from scitbx import matrix
import scitbx.rigid_body
from cctbx import xray
import string
import random
import math
import sys
from copy import deepcopy
import mmtbx.monomer_library.server
from mmtbx.refinement.flip_peptide_side_chain import should_be_flipped, \
    flip_residue


__author__ = 'Youval, massively rewritten by Oleg'

def flip_atoms_in_ncs_groups(hierarchy, ncs_restraints_group_list, mon_lib_srv=None):
  """
  XXX
  XXX not used, only tested. May have some value.
  XXX
  This function will actually modify hierarchy by making necessary flips
  in ncs-related residues. Flip will be made by exchanging atom coordinates.
  Will make all copies consistent with master.
  """
  if mon_lib_srv is None:
    mon_lib_srv = mmtbx.monomer_library.server.server()
  for ncs_gr in ncs_restraints_group_list:
    master_isel = ncs_gr.master_iselection
    chains_master = hierarchy.select(master_isel).only_model().chains()
    for copy in ncs_gr.copies:
      copy_isel = copy.iselection
      chains_copy = hierarchy.select(copy_isel).only_model().chains()
      for ch_m, ch_c in zip(chains_master, chains_copy):
        for r_m, r_c in zip(ch_m.residues(), ch_c.residues()):
          # print "working on ", r_m.id_str(), r_c.id_str()
          if should_be_flipped(r_m, r_c):
            flip_residue(r_c, mon_lib_srv)


def concatenate_rot_tran(transforms_obj=None,
                         ncs_restraints_group_list=None):
  """
  XXX
  XXX Remove transforms_obj parameter
  XXX Make it method of class_ncs_restraints_group_list
  XXX

  Concatenate rotation angles, corresponding to the rotation
  matrices and scaled translation vectors to a single long flex.double object

  Args:
    transforms_obj : (mmtbx.refinement.minimization_ncs_constraints
      ncs_group_object) containing information on Rotation matrices (lists of
      objects matrix.rec) and Translation vectors (lists of objects matrix.rec)
    ncs_restraints_group_list : a list of ncs_restraint_group objects

  Returns:
    flex.double : [(alpha_1,beta_1,gamma_1,Tx_1,Ty_1,Tz_1)...]
  """
  x = []
  if (not ncs_restraints_group_list) and transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  if ncs_restraints_group_list:
    for gr in ncs_restraints_group_list:
      for tr in gr.copies:
        x.extend(list(rotation_to_angles(rotation=tr.r.elems))
                 + list(tr.t.elems))
  return flex.double(x)

def get_rotation_translation_as_list(transforms_obj=None,
                                     ncs_restraints_group_list=None):
  """
  XXX
  XXX Consider deletion. Used only in tests tst_minimization_ncs_constraints_real_space.py,
  XXX tst_ncs_utils.py
  XXX Alternatively, make it method of class_ncs_restraints_group_list
  XXX

  Get rotations and translations vectors from ncs_restraints_group_list or
  transforms_obj

  Returns:
    r (list): list of rotation matrices
    t (list): list of translation vectors
  """
  r = []
  t = []
  if (not ncs_restraints_group_list) and transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  if ncs_restraints_group_list:
    for nrg in ncs_restraints_group_list:
      for tr in nrg.copies:
        r.append(tr.r)
        t.append(tr.t)
  return r,t

def update_transforms(transforms_obj,rm,tv):
  """
  XXX
  XXX Consider removing it. Only used in
  XXX mmtbx/regression/ncs/tst_minimization_ncs_constraints_real_space.py
  XXX

  Update transforms_obj with the rotation matrices (rm) and translation
  vectors (tv) """
  assert len(transforms_obj.transform_order) == len(rm)
  assert len(rm) == len(tv)
  for tr,r,t in zip(transforms_obj.transform_order,rm,tv):
    transforms_obj.ncs_transform[tr].r = r
    transforms_obj.ncs_transform[tr].t = t
  return transforms_obj

# MARKED_FOR_DELETION_OLEG
# REASON: not used not tested
def update_ncs_restraints_group_list(ncs_restraints_group_list,rm,tv):
  """
  Update ncs_restraints_group_list with the rotation matrices (rm) and
  translation vectors (tv) """
  assert len(rm) == len(tv)
  new_list = []
  for gr in ncs_restraints_group_list:
    for tr in gr.copies:
      tr.r = rm.pop(0)
      tr.t = tv.pop(0)
    new_list.append(gr)
  return new_list
# END_MARKED_FOR_DELETION_OLEG

def update_rot_tran(x,transforms_obj=None,ncs_restraints_group_list=None):
  """
  XXX
  XXX Make it method of class_ncs_restraints_group_list
  XXX Remove transforms_obj parameter.
  XXX

  Convert the refinable parameters, rotations angles and
  scaled translations, back to rotation matrices and translation vectors and
  updates the transforms_obj (ncs_restraints_group_list)

  Args:
    x : a flex.double of the form (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
      theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s). where n is the number of
      transformations.
    transforms_obj : (ncs_group_object) containing information on Rotation
      matrices, Translation vectors and NCS
  ncs_restraints_group_list : a list of ncs_restraint_group objects

  Returns:
    The same type of input object with converted transforms
  """
  assert bool(transforms_obj) == (not bool(ncs_restraints_group_list))
  if transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  if ncs_restraints_group_list:
    i = 0
    for gr in ncs_restraints_group_list:
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
    if transforms_obj:
      transforms_obj.update_using_ncs_restraints_group_list(
        ncs_restraints_group_list)
      return transforms_obj
    else:
      return ncs_restraints_group_list

def rotation_to_angles(rotation, deg=False):
  """
  Get the rotation angles around the axis x,y,z for rotation r
  Such that r = Rx*Ry*Rz
  Those angles are the Tait-Bryan angles form of Euler angles

  Note that typically there are two solutions, and this function will return
  only one. In the case that cos(beta) == 0 there are infinite number of
  solutions, the function returns the one where gamma = 0

  Args:
    r : (flex.double) of the form (Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz)
    deg : When False use radians, when True use degrees

  Returns:
    angles (flex.double): (alpha, beta, gamma) rotation angles around the x,y,z
  """
  # make sure the rotation data type is flex.double
  if not isinstance(rotation,type(flex.double([]))):
    rotation = flex.double(rotation)
  (Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz) = rotation.round(8)
  if Rxz not in [1,-1]:
    beta = math.asin(Rxz)
    # beta2 = math.pi - beta
    # using atan2 to take into account the possible different angles and signs
    alpha = math.atan2(-Ryz/math.cos(beta),Rzz/math.cos(beta))
    gamma = math.atan2(-Rxy/math.cos(beta),Rxx/math.cos(beta))
    # alpha2 = math.atan2(-Ryz/math.cos(beta2),Rzz/math.cos(beta2))
    # gamma2 = math.atan2(-Rxy/math.cos(beta2),Rxx/math.cos(beta2))
  elif Rxz == 1:
    beta = math.pi/2
    alpha = math.atan2(Ryx,Ryy)
    gamma = 0.0
  elif Rxz == -1:
    beta = -math.pi/2
    alpha = math.atan2(-Ryx,Ryy)
    gamma = 0.0
  else:
    raise ArithmeticError("Can't calculate rotation angles")

  angles = flex.double((alpha,beta,gamma))
  # angles2 = flex.double((alpha2,beta2,gamma2))

  if deg:
    # Convert to degrees
    angles = 180*angles/math.pi
    angles = angles.round(5)
    # angles2 = 180*angles2/math.pi
  return angles

def shake_transformations(x,
                          shake_angles_sigma      = 0.035,
                          shake_translation_sigma = 0.5):
  """
  XXX
  XXX Used here in get_weight().
  XXX Not clear what relation MTRIX have to this function at all...
  XXX

  Shake rotation matrices and translation vectors of a rotation matrices and
  translation vectors from the MTRIX records in a PDB file.

  Args:
    x (flex.double): [(alpha_1,beta_1,gamma_1,Tx_1/s,Ty_1/s,Tz_1/s)...]
    shake_angles_sigma (float): the sigma (in radians) of the random gaussian
      shaking of the rotation angles
    shake_translation_sigma (float): the sigma (in angstrom) of the random
      gaussian shaking of the translation

  Return:
    new_x (flex.double): The shaken x
  """
  new_x = flex.double([])
  for i in xrange(0,len(x),6):
    new_x.append(random.gauss(x[i+0],shake_angles_sigma))
    new_x.append(random.gauss(x[i+1],shake_angles_sigma))
    new_x.append(random.gauss(x[i+2],shake_angles_sigma))
    new_x.append(random.gauss(x[i+3],shake_translation_sigma))
    new_x.append(random.gauss(x[i+4],shake_translation_sigma))
    new_x.append(random.gauss(x[i+5],shake_translation_sigma))
  return new_x

def compute_transform_grad(grad_wrt_xyz,
                           xyz_asu,
                           x,
                           ncs_restraints_group_list=None,
                           transforms_obj=None):
  """
  XXX
  XXX Remove transforms_obj from parameter list.
  XXX Consider making it method of class_ncs_restraints_group_list
  XXX

  Compute gradient in respect to the rotation angles and the translation
  vectors. R = Rx(the)Ry(psi)Rz(phi)

  Args:
    grad_wrt_xyz (flex.double): gradients with respect to xyz.
    ncs_restraints_group_list: list containing ncs_restraint_group objects
    transforms_obj (ncs_group_object): containing information in rotation
      matrices and to which chains they apply
    xyz_asu (flex.vec3): The coordinates sites cart of the complete ASU
    x (flex double): The angles, in the form
      (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
      theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s)

  Returns:
    g (flex.double): the gradient
  """
  assert bool(transforms_obj) == (not bool(ncs_restraints_group_list))
  if transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  g = []
  grad_wrt_xyz = flex.vec3_double(grad_wrt_xyz)
  i = 0
  for nrg in ncs_restraints_group_list:
    xyz_ncs_transform = xyz_asu.select(nrg.master_iselection)
    xyz_len = xyz_ncs_transform.size()
    # calc the coordinates of the master NCS at its coordinates center system
    mu_c = flex.vec3_double([xyz_ncs_transform.sum()]) * (1/xyz_len)
    xyz_cm = xyz_ncs_transform - flex.vec3_double(list(mu_c) * xyz_len)
    for nrg_copy in nrg.copies:
      grad_ncs_wrt_xyz = grad_wrt_xyz.select(nrg_copy.iselection)
      assert xyz_len == grad_ncs_wrt_xyz.size()
      grad_wrt_t = list(grad_ncs_wrt_xyz.sum())
      # Sum angles gradient over the coordinates
      # Use the coordinate center for rotation
      m = grad_ncs_wrt_xyz.transpose_multiply(xyz_cm)
      m = matrix.sqr(m)
      # Calculate gradient with respect to the rotation angles
      the,psi,phi = x[i*6:i*6+3]
      rot = scitbx.rigid_body.rb_mat_xyz(
        the=the, psi=psi, phi=phi, deg=False)
      g_the = (m * rot.r_the().transpose()).trace()
      g_psi = (m * rot.r_psi().transpose()).trace()
      g_phi = (m * rot.r_phi().transpose()).trace()
      g.extend([g_the, g_psi, g_phi])
      g.extend(grad_wrt_t)
      i += 1
  return flex.double(g)

def get_ncs_sites_cart(ncs_obj=None,
                       fmodel=None,
                       xray_structure=None,
                       sites_cart=None,
                       extended_ncs_selection=None):
  """
  XXX The same name of the function as in cctbx/maptbx/segment_and_split_map.py
  XXX The only call is from minimization_ncs_constraints.py, passing whole
  XXX lbfgs class as 'ncs_obj'!
  XXX TODO: remove ncs_obj parameter, pass everything needed explicitly
  XXX not yet clear the goal of this func
  XXX

  Args::
    ncs_obj: an object that contains fmodel, sites_cart or xray_structure
      and an atom selection flags for a single NCS copy.

  Returns:
    (flex.vec3): coordinate sites cart of the single NCS copy
  """
  if ncs_obj:
    if hasattr(ncs_obj, 'extended_ncs_selection'):
      extended_ncs_selection = ncs_obj.extended_ncs_selection
    else:
      assert extended_ncs_selection
    if hasattr(ncs_obj, 'sites_cart'):
      return ncs_obj.sites_cart().select(extended_ncs_selection)
    elif hasattr(ncs_obj, 'fmodel'):
      xrs_one_ncs = ncs_obj.fmodel.xray_structure.select(extended_ncs_selection)
      return xrs_one_ncs.sites_cart()
    elif  hasattr(ncs_obj, 'xray_structure') or xray_structure:
      xrs_one_ncs = ncs_obj.xray_structure.sites_cart()
      return xrs_one_ncs.select(extended_ncs_selection)
  else:
    if sites_cart:
      return sites_cart().select(extended_ncs_selection)
    elif fmodel:
      xrs_one_ncs = fmodel.xray_structure.select(extended_ncs_selection)
      return xrs_one_ncs.sites_cart()
    elif  xray_structure:
      xrs_one_ncs = xray_structure.sites_cart()
      return xrs_one_ncs.select(extended_ncs_selection)


def get_weight(fmodel=None,
               restraints_manager=None,
               sites=None,
               transformations=None,
               u_iso=None,
               ncs_restraints_group_list=None,
               refine_selection=None,
               minimized_obj=None):
  """
  Calculates weights for refinements by slightly shaking the minimized
  parameters and taking the ratio:
  (restraint manager grad norm / parameters gradient norm)


  When calling this function during refinement macro cycle, the minimized
  object, "minimized_obj" , may contains fmodel, restraints_manager,
  refinement type info (sites, transformations, u_iso) and
  ncs_restraints_group_list.


  Args:
    fmodel : F-model object
    restraints_manager: Restraints manager object
    sites (bool): Refine by sites
    u_iso (bool): Refine using u_iso
    transformations (bool): Refine using transformations
      rotations, translations (matrix objects):
    ncs_restraints_group_list: list of ncs_restraint_group objects
    refine_selection (flex.size_t): selection of all ncs related copies and
      non ncs related parts to be included in selection (to be refined)
    minimized_obj:Minimization object containing all the other
      parameters above

  Returns:
    weight (int):

  Example:
  >>>get_weight(minimized_obj=minimized_obj)

  or

  >>>get_weight(fmodel=fmodel,
                restraints_manager=grm,
                sites=sites,
                transformations=transformations,
                u_iso=u_iso,
                ncs_restraints_group_list=ncs_restraints_group_list,
                refine_selection=refine_selection)
  """
  grm  = restraints_manager
  extended_ncs_selection = None
  # extract parameters from minimized_obj
  if minimized_obj:
    mo = minimized_obj
    fmodel = mo.fmodel
    grm = mo.grm
    sites = mo.sites
    transformations = mo.transformations
    u_iso = mo.u_iso
    ncs_restraints_group_list = mo.ncs_restraints_group_list
  # del: consider deleting the code below
  #   if hasattr(mo,'extended_ncs_selection'):
  #     extended_ncs_selection = mo.extended_ncs_selection
  # if not extended_ncs_selection:
  #   extended_ncs_selection = get_extended_ncs_selection(
  #     ncs_restraints_group_list=ncs_restraints_group_list,
  #     refine_selection=refine_selection)

  # make sure sufficient input is provided
  assert bool(fmodel), 'F-model is not provided'
  assert bool(grm), 'F-restraints_manager is not provided'
  assert [sites,transformations,u_iso].count(True)==1, 'Refinement type Error'
  #
  have_transforms = ncs_restraints_group_list != []
  fmdc = fmodel.deep_copy()
  if sites:
    fmdc.xray_structure.shake_sites_in_place(mean_distance=0.3)
  elif u_iso:
    fmdc.xray_structure.shake_adp()
  elif transformations and have_transforms:
    x = concatenate_rot_tran(
      ncs_restraints_group_list = ncs_restraints_group_list)
    x = shake_transformations(
      x = x,
      shake_angles_sigma=0.035,
      shake_translation_sigma=0.5)
  fmdc.update_xray_structure(xray_structure = fmdc.xray_structure,
    update_f_calc=True)
  fmdc.xray_structure.scatterers().flags_set_grads(state=False)
  if sites or transformations:
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # restraints manager, energy sites gradients
    gc = grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
  elif u_iso:
    # Create energies_site gradient, to create geometry_restraints_manager
    # plain_pair_sym_table needed for the energies_adp_iso
    import mmtbx.refinement.adp_refinement
    temp = mmtbx.refinement.adp_refinement.adp_restraints_master_params
    iso_restraints = temp.extract().iso
    gc = grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      u_iso      = True)
    # fmodel gradients
    gxc = fmdc.one_time_gradients_wrt_atomic_parameters(
      u_iso = True).as_double()
    # manager restraints, energy sites gradients
    gc = grm.energies_adp_iso(
      xray_structure    = fmdc.xray_structure,
      parameters        = iso_restraints,
      use_u_local_only  = iso_restraints.use_u_local_only,
      use_hd            = False,
      compute_gradients = True).gradients
  if transformations and have_transforms:
    # Apply NCS relations to gradients
    gxc = compute_transform_grad(
      grad_wrt_xyz      = gxc.as_double(),
      ncs_restraints_group_list = ncs_restraints_group_list,
      xyz_asu           = fmdc.xray_structure.sites_cart(),
      x                 = x)
    gc = compute_transform_grad(
      grad_wrt_xyz      = gc.as_double(),
      ncs_restraints_group_list = ncs_restraints_group_list,
      xyz_asu           = fmdc.xray_structure.sites_cart(),
      x                 = x)

  weight = 1.
  gc_norm  = gc.norm()
  gxc_norm = gxc.norm()
  if(gxc_norm != 0.0):
    weight = gc_norm / gxc_norm

  weight =min(weight,1e6) # limit the weight max value
  return weight

def apply_transforms(ncs_coordinates,
                     ncs_restraints_group_list,
                     total_asu_length,
                     extended_ncs_selection,
                     round_coordinates = True,
                     center_of_coordinates = None):
  """
  Apply transformation to ncs_coordinates,
  and round the results if round_coordinates is True

  Args:
    ncs_coordinates (flex.vec3): master ncs coordinates
    ncs_restraints_group_list: list of ncs_restraint_group objects
    total_asu_length (int): Complete ASU length
      extended_ncs_selection (flex.size_t): master ncs and non-ncs related parts
    center_of_coordinates : when not None, contains the center of coordinate of
      the master for each ncs copy

  Returns:
    (flex.vec3_double): Asymmetric or biological unit parts that are related via
      ncs operations
  """
  asu_xyz = flex.vec3_double([(0,0,0)]*total_asu_length)
  asu_xyz.set_selected(extended_ncs_selection,ncs_coordinates)

  # get the rotation and translation for the native coordinate system
  if bool(center_of_coordinates):
    ncs_restraints_group_list = ncs_restraints_group_list.shift_translation_back_to_place(
        shifts = center_of_coordinates)
  for nrg in ncs_restraints_group_list:
    master_ncs_selection = flex.bool(total_asu_length,nrg.master_iselection)
    for ncs_copy in nrg.copies:
      copy_selection = flex.bool(total_asu_length,ncs_copy.iselection)
      ncs_xyz = asu_xyz.select(master_ncs_selection)
      new_sites = ncs_copy.r.elems * ncs_xyz + ncs_copy.t
      asu_xyz.set_selected(copy_selection,new_sites)
  if round_coordinates:
    return flex.vec3_double(asu_xyz).round(3)
  else:
    return flex.vec3_double(asu_xyz)

def ncs_groups_selection(ncs_restraints_group_list,selection):
  """
  XXX
  XXX Move to iotbx/ncs/__init__.py, member func of class_ncs_restraints_group_list()
  XXX Figure out difference between this and def select() there.
  XXX

  Modifies the selections of master and copies according the "selection"
  - Keep the order of selected atoms
  - Keep only atoms that appear in master and ALL copies
  Also modify "selection" to include ncs related atoms only if selected in
  both master and ALL ncs copies (The modified selection is not returned in
  current version)

  Args:
    ncs_restraints_group_list (list): list of ncs_restraints_group objects
    selection (flex.bool or flex.size_t): atom selection

  Returns:
    new_nrg_list (list): list of modified ncs_restraints_group objects
  """
  if isinstance(selection,flex.bool): selection = selection.iselection(True)
  sel_set = set(selection)
  new_nrg_list = ncs_restraints_group_list.deep_copy()
  # check what are the selection that shows in both master and all copies
  for nrg in new_nrg_list:
    m = set(nrg.master_iselection)
    m_list = [(pos,indx) for pos,indx in enumerate(list(nrg.master_iselection))]
    m_in_sel = m.intersection(sel_set)
    common_selection_pos = {pos for (pos,indx) in m_list if indx in m_in_sel}
    for ncs in nrg.copies:
      c = set(ncs.iselection)
      c_list = [(pos,indx) for pos,indx in enumerate(list(ncs.iselection))]
      copy_in_sel = c.intersection(sel_set)
      include_set = {pos for (pos,indx) in c_list if indx in copy_in_sel}
      common_selection_pos.intersection_update(include_set)
      if not bool(common_selection_pos): break
    # use the common_selection_pos to update all selections
    nrg.master_iselection, not_included = selected_positions(
      nrg.master_iselection,common_selection_pos)
    selection = remove_items_from_selection(selection,not_included)
    for ncs in nrg.copies:
      ncs.iselection, not_included = selected_positions(
        ncs.iselection,common_selection_pos)
      selection = remove_items_from_selection(selection,not_included)

  return new_nrg_list


def selected_positions(selection,positions):
  """
  Returns only the selected indices in the positions specified in "positions"
  keeping the order

  Args:
    selection (flex.size_t): Atoms selection
    positions (set or list): the allowed positions in the selections

  Returns:
    (flex.size_t, flex.size_t): (selected atoms, atoms, not selected)

  Examples::
    >>>a = flex.size_t([1,2,5,6,4])
    >>>pos = {0,3,4}
    >>>s,d = selected_positions(a,pos)
    >>>list(s)
    [1,6,4]
    >>>list(d)
    [2,5]
  """
  assert isinstance(selection,flex.size_t)
  if isinstance(positions,set): positions = flex.size_t(list(positions))
  if isinstance(positions,list): positions = flex.size_t(positions)
  include = flex.bool(selection.size(),positions)
  not_include = ~include
  return selection.select(include), selection.select(not_include)

def remove_items_from_selection(selection,remove):
  """
  Remove a set of atoms from "selection"

  Args:
    selection (flex.size_t): atom selection
    remove (flex.size_t): atoms to remove from selection

  Returns:
    (flex.size_t): modified atom selection

  Examples::
    >>>a = flex.size_t([1,2,5,6,4])
    >>>r = flex.size_t([2,5])
    >>>s = remove_items_from_selection(a,r,10)
    >>>list(s)
    [1,6,4]
  """
  selection = list(selection)
  remove = set(remove)
  new_selection = [x for x in selection if not (x in remove)]
  return flex.size_t(new_selection)

def get_list_of_best_ncs_copy_map_correlation(
      ncs_groups,
      xray_structure=None,
      fmodel=None,
      map_data=None,
      d_min=None):
  """
  Finds the copy with best map correlation in each ncs group

  Returns:
    best_list (list of int): list of the copy with the best map correlation.
      (the master copy is 0)
  """
  assert [fmodel,d_min].count(None) in [0,1]
  assert [fmodel,map_data].count(None)==1
  assert [d_min,map_data].count(None) in [0,2]
  assert [xray_structure, fmodel].count(None)==1
  import mmtbx.maps.correlation
  best_list = []
  if(fmodel is None):
    mp = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      xray_structure = xray_structure,
      fmodel         = fmodel,
      map_data       = map_data,
      d_min          = d_min)
  else:
    mp = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      fmodel = fmodel)
  for nrg in ncs_groups:
    selections = [nrg.master_iselection]
    for ncs in nrg.copies:
      selections.append(ncs.iselection)
    cc = mp.cc(selections=selections)
    i_seq = cc.index(max(cc)) # best matching copy
    if(i_seq == 0): continue
    #
    c_i = i_seq-1
    # switch master and copy selection
    nrg.master_iselection, nrg.copies[c_i].iselection = \
      nrg.copies[c_i].iselection, nrg.master_iselection
    # Adjust rotation and translation for the new master
    r = nrg.copies[c_i].r = (nrg.copies[c_i].r.transpose())
    t = nrg.copies[c_i].t = -(nrg.copies[c_i].r * nrg.copies[c_i].t)
    # change all other rotations and translations to the new master
    for i in xrange(len(nrg.copies)):
      if i == c_i: continue
      # change translation before rotation
      nrg.copies[i].t = (nrg.copies[i].r * t + nrg.copies[i].t)
      nrg.copies[i].r = (nrg.copies[i].r * r)

def get_refine_selection(refine_selection=None,number_of_atoms=None):
  """ populate refine_selection with all atoms if no selection is given  """
  if not bool(refine_selection):
      # select to refine all atoms
      assert bool(number_of_atoms)
      selection_list = range(number_of_atoms)
      refine_selection = flex.size_t(selection_list)
  return refine_selection

def make_unique_chain_names(unique_chain_names,number_of_names=1):
  """
  XXX
  XXX Consider removing, unifying with iotbx.pdb.utils.all_chain_ids
  XXX

  Produce a sorted list of new unique chain names.
  Chain names are strings of one or two characters long.

  Args:
    unique_chain_names (set): Current names
    number_of_names (int): number of new names

  Returns:
    new_names_list (list): sorted list on new names
  """
  # check availability of one letter chain names
  chr_list1 = list(set(string.ascii_uppercase) - set(unique_chain_names))
  chr_list2 = list(set(string.ascii_lowercase) - set(unique_chain_names))
  chr_list1.sort()
  chr_list2.sort()
  new_names_list = chr_list1 + chr_list2
  if len(new_names_list) < number_of_names:
    # calc how many more chain names we need
    n_names =  number_of_names - len(new_names_list)
    # the number of character needed to produce new names
    chr_number = int(math.sqrt(n_names)) + 1
    # build character list
    chr_list = list(string.ascii_uppercase) + \
               list(string.ascii_lowercase) + \
               list(string.digits)
    # take only as many characters as needed
    chr_list = chr_list[:chr_number]
    extra_names = set([ x+y for x in chr_list for y in chr_list])
    # make sure not using existing names
    extra_names = list(extra_names - set(unique_chain_names))
    extra_names.sort()
    new_names_list.extend(extra_names)
  return new_names_list[:number_of_names]
