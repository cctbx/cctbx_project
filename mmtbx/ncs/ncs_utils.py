from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx import matrix
import scitbx.rigid_body
from cctbx import xray
import random
import math
import mmtbx.monomer_library.server
from mmtbx.refinement.flip_peptide_side_chain import should_be_flipped, \
    flip_residue
from six.moves import zip
from six.moves import range


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
  for i in range(0,len(x),6):
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
                           ncs_restraints_group_list):
  """
  XXX
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
    x = ncs_restraints_group_list.concatenate_rot_tran()
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
    selections = nrg.get_iselections_list()
    cc = mp.cc(selections=selections)
    i_seq = cc.index(max(cc)) # best matching copy
    if(i_seq == 0): continue
    #
    c_i = i_seq-1
    nrg.make_nth_copy_master(c_i)

def get_refine_selection(refine_selection=None,number_of_atoms=None):
  """ populate refine_selection with all atoms if no selection is given  """
  if not bool(refine_selection):
      # select to refine all atoms
      assert bool(number_of_atoms)
      selection_list = range(number_of_atoms)
      refine_selection = flex.size_t(selection_list)
  return refine_selection
