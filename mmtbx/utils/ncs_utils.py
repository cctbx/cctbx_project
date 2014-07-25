from __future__ import division
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from mmtbx import monomer_library
from scitbx import matrix
import scitbx.rigid_body
from cctbx import xray
import random
import math


def concatenate_rot_tran(transforms_obj=None,
                         ncs_restraints_group_list=None,
                         deg=True, s=1):
  """
  Concatenate rotation angles, corresponding to the rotation
  matrices and scaled translation vectors to a single long flex.double object

  Arguments:
  transforms_obj : (mmtbx.refinement.minimization_ncs_constraints
                   ncs_group_object) containing information on Rotation
                   matrices(lists of objects matrix.rec) and Translation vectors
                   (lists of objects matrix.rec)
  ncs_restraints_group_list : a list of ncs_restraint_group objects
  s : (float) scaling factor, scale translation to be of the
      same order ot the translation

  Return:
  flex.double object of the form
  [(alpha_1,beta_1,gamma_1,Tx_1/s,Ty_1/s,Tz_1/s)...]
  """
  x = []
  if (not ncs_restraints_group_list) and transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  if ncs_restraints_group_list:
    for gr in ncs_restraints_group_list:
      for tr in gr.copies:
        x.extend(list(rotation_to_angles(rotation=tr.r.elems,deg=deg))
                 + list((tr.t/s).elems))
  return flex.double(x)

def get_rotation_translation_as_list(transforms_obj=None,
                                     ncs_restraints_group_list=None):
  """
  Collect and returns rotations matrices and translations vectors to two lists
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
  """ Update of the rotation matrices (rm) and translation vectors (vt) """
  assert len(transforms_obj.transform_order) == len(rm)
  assert len(rm) == len(tv)
  for tr,r,t in zip(transforms_obj.transform_order,rm,tv):
    transforms_obj.ncs_transform[tr].r = r
    transforms_obj.ncs_transform[tr].t = t
  return transforms_obj

def update_ncs_restraints_group_list(ncs_restraints_group_list,rm,tv):
  """ Update of the rotation matrices (rm) and translation vectors (tv) """
  assert len(rm) == len(tv)
  new_list = []
  for gr in ncs_restraints_group_list:
    for tr in gr.copies:
      tr.r = rm.pop(0)
      tr.t = tv.pop(0)
    new_list.append(gr)
  return new_list

def update_rot_tran(x,s=1.0,
                    transforms_obj=None,
                    ncs_restraints_group_list=None,
                    deg=True):
  """
  Convert the refinemable parameters, rotations angles and
  scaled translations, back to rotation matrices and translation vectors and
  updates the transforms_obj (ncs_restraints_group_list)

  Arguments:
  x : a flex.double of the form (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
      theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s). where n is the number of
      transformations.
  transforms_obj : (ncs_group_object) containing information on Rotation
                   matrices, Translation vectors and NCS
  ncs_restraints_group_list : a list of ncs_restraint_group objects
  s : (float) scaling factor, scale translation to be of the
      same order ot the translation

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
          the=the, psi=psi, phi=phi, deg=deg)
        tran = matrix.rec(x[i*6+3:i*6+6],(3,1))*s
        tr.r = (rot.rot_mat()).round(8)
        tr.t = tran.round(8)
        copies.append(tr)
        i += 1
      gr.copies = copies
    if transforms_obj:
      transforms_obj.update_using_ncs_restraints_group_list(
        ncs_restraints_group_list)
      return transforms_obj
    else:
      return ncs_restraints_group_list

def rotation_to_angles(rotation, deg=True):
  """
  Get the rotation angles around the axis x,y,x for rotation r
  Such that r = Rx*Ry*Rz

  Note that typically there are two solutions, and this function will return
  only one. In the case that cos(beta) == 0 there are infinite number of
  solutions, the function returns the one where gamma = 0

  Arguments:
  r : (flex.double) of the form (Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz)
  deg : When False use radians, when True use degrees

  Return:
  angles: (flex.double) containing rotation angles in  the form
          (alpha, beta, gamma)
  """
  # make sure the rotation data type is flex.double
  if not isinstance(rotation,type(flex.double())):
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
    gamma = 0
  elif Rxz == -1:
    beta = -math.pi/2
    alpha = math.atan2(-Ryx,Ryy)
    gamma = 0
  else:
    raise ArithmeticError("Can't calculate rotation angles")

  angles = flex.double((alpha,beta,gamma))
  # angles2 = flex.double((alpha2,beta2,gamma2))

  if deg:
    # Convert to degrees
    angles = 180*angles/math.pi
    # angles2 = 180*angles2/math.pi
  return angles.round(5)

def angles_to_rotation(angles_xyz, deg=True, rotation_is_tuple=False):
  """
  Calculate rotation matrix R, such that R = Rx(alpha)*Ry(beta)*Rz(gamma)

  Arguments:
  angles_xyz : (flex.double) (alpha,beta,gamma)
  deg : (bool) When False use radians, when True degrees
  rotation_is_tuple : (bool) when False, return flxe.double object,
                      when True return tuple

  Returns:
  R : (tuple or flex.double) the components of a rotation matrix
  """
  assert len(angles_xyz) == 3
  alpha,beta,gamma = angles_xyz
  rot = scitbx.rigid_body.rb_mat_xyz(the=alpha, psi=beta, phi=gamma, deg=deg)
  R = rot.rot_mat()
  if rotation_is_tuple:
    return R.round(6).elems
  else:
    return flex.double(R.round(6))

def shake_transformations(x,
                          shake_angles_sigma      = 0.035,
                          shake_translation_sigma = 0.5):
  """
  Shake rotation matrices and translation vectors of a rotation matrices and
  translation vectors from the MTRIX records in a PDB file.

  Argument:
  x: flex.double object of the form
     [(alpha_1,beta_1,gamma_1,Tx_1/s,Ty_1/s,Tz_1/s)...]
  shake_angles_sigma: (float) the sigma (in radians) of the random gaussian
                      shaking of the rotation angles
  shake_translation_sigma: (float) the sigma (in angstrom) of the random
                           gaussian shaking of the translation

  Return:
  Shaken x, rotation matrices and translation vectors
  """
  new_x = flex.double()
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
                           transforms_obj=None,
                           deg=True):
  """
  Compute gradient in respect to the rotation angles and the translation
  vectors. R = Rx(the)Ry(psi)Rz(phi)

  Arguments:
  grad_wrt_xyz : (flex.double) gradients with respect to xyz.
  ncs_restraints_group_list: list containing ncs_restraint_group objects
  transforms_obj : (ncs_group_object) containing information in rotation
                    matrices and to which chains they apply
  xyz_asu: (flex.vec3) The coordinates sites cart of the complete ASU
  x: (flex double) The angles, in the form (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
      theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s)

  Returns:
  g : (flex.double) a gradient
  """
  assert bool(transforms_obj) == (not bool(ncs_restraints_group_list))
  if transforms_obj:
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
  g = []
  grad_wrt_xyz = flex.vec3_double(grad_wrt_xyz)
  i = 0
  for nrg in ncs_restraints_group_list:
    xyz_ncs_transform = xyz_asu.select(nrg.master_iselection)
    for nrg_copy in nrg.copies:
      grad_ncs_wrt_xyz = grad_wrt_xyz.select(nrg_copy.copy_iselection)
      xyz_len = xyz_ncs_transform.size()
      assert xyz_len == grad_ncs_wrt_xyz.size()
      grad_wrt_t = list(grad_ncs_wrt_xyz.sum())
      # Use the coordinate center for rotation
      mu_c = flex.vec3_double([xyz_ncs_transform.sum()]) * (1/xyz_len)
      xyz_cm = xyz_ncs_transform - flex.vec3_double(list(mu_c) * xyz_len)
      # Sum angles gradient over the coordinates
      m = grad_ncs_wrt_xyz.transpose_multiply(xyz_cm)
      m = matrix.sqr(m)
      # Calculate gradient with respect to the rotation angles
      the,psi,phi = x[i*6:i*6+3]
      rot = scitbx.rigid_body.rb_mat_xyz(
        the=the, psi=psi, phi=phi, deg=deg)
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
  Argument:
  ncs_obj: an object that contains fmodel, sites_cart or xray_structure
           and an atom selection flags for a single NCS copy.

  Return: (flex.vec3) coordinate sites cart of the single NCS copy
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


def get_weight(minimization_obj=None,
               fmodel=None,
               grm=None,
               sites=None,
               transformations=None,
               u_iso=None,
               ncs_restraints_group_list=None,
               refine_selection=None):
  """
  Calculates weights for refinements

  Minimization object must contain the following methods and attributes:
  - fmodel
  - sites, u_iso, transformations: (bool)
  - rotations, translations: (matrix objects)
  - grm: Restraints manager
  - iso_restraints: (libtbx.phil.scope_extract)
                    object used for u_iso refinement parameters
  - ncs_atom_selection: (flex bool) for a single ncs atom selection

  Return:
  weight: (int)
  """
  if minimization_obj:
    mo = minimization_obj
    fmodel = mo.fmodel
    grm = mo.grm
    sites = mo.sites
    transformations = mo.transformations
    u_iso = mo.u_iso
    ncs_restraints_group_list = mo.ncs_restraints_group_list
    extended_ncs_selection = mo.extended_ncs_selection
  else:
    extended_ncs_selection = get_extended_ncs_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      refine_selection=refine_selection)

  have_transforms = ncs_restraints_group_list != []
  assert [sites,u_iso,transformations].count(True)==1
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
  if sites:
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # manager restraints, energy sites gradients
    gc = grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
  elif u_iso:
    # Create energies_site gradient, to create
    # geometry_restraints_manager.plain_pair_sym_table
    # needed for the energies_adp_iso
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
  elif transformations and have_transforms:
    xyz_ncs = get_ncs_sites_cart(
      fmodel=fmodel, extended_ncs_selection=extended_ncs_selection)
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc_xyz = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # manager restraints, energy sites gradients
    gc_xyz = grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
    gxc = compute_transform_grad(
      grad_wrt_xyz      = gxc_xyz.as_double(),
      ncs_restraints_group_list = ncs_restraints_group_list,
      xyz_asu           = fmdc.xray_structure.sites_cart(),
      x                 = x)
    gc = compute_transform_grad(
      grad_wrt_xyz      = gc_xyz.as_double(),
      ncs_restraints_group_list = ncs_restraints_group_list,
      xyz_asu           = fmdc.xray_structure.sites_cart(),
      x                 = x)

  weight = 1.
  gc_norm  = gc.norm()
  gxc_norm = gxc.norm()
  if(gxc_norm != 0.0):
    weight = gc_norm / gxc_norm

  weight =min(weight,1e6)
  return weight

def get_restraints_manager(pdb_file_name=None,pdb_string=None):
  """  Generate restraint manager from a PDB file or a PDB string  """
  assert [pdb_file_name,pdb_string].count(None)==1
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if pdb_string: pdb_lines = pdb_string.splitlines()
  else: pdb_lines = None
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = pdb_file_name,
    raw_records    = pdb_lines,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = False, plain_pairs_radius = 5.0)
  return mmtbx.restraints.manager(
    geometry = geometry, normalization = False)

def apply_transforms(ncs_coordinates,
                     ncs_restraints_group_list,
                     total_asu_length,
                     extended_ncs_selection,
                     round_coordinates = True):
  """
  Apply transformation to ncs_coordinates,
  and round the results if round_coordinates is True

  Argument:
  ncs_coordinates: (flex.vec3) master ncs coordinates
  ncs_restraints_group_list: list of ncs_restraint_group objects
  total_asu_length: (int) Complete ASU length
  extended_ncs_selection: (flex.size_t) master ncs and non-ncs related parts

  Returns:
  Asymmetric or biological unit parts that are related via ncs operations
  """
  asu_xyz = flex.vec3_double([(0,0,0)]*total_asu_length)
  asu_xyz.set_selected(extended_ncs_selection,ncs_coordinates)

  for nrg in ncs_restraints_group_list:
    master_ncs_selection = nrg.master_iselection
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.copy_iselection
      ncs_xyz = asu_xyz.select(master_ncs_selection)
      new_sites = ncs_copy.r.elems * ncs_xyz + ncs_copy.t
      asu_xyz.set_selected(asu_selection,new_sites)
  if round_coordinates:
    return flex.vec3_double(asu_xyz).round(3)
  else:
    return flex.vec3_double(asu_xyz)

def get_extended_ncs_selection(ncs_restraints_group_list,refine_selection):
  """
  :param ncs_restraints_group_list: list of ncs_restraint_group objects
  total_asu_length: (int) Complete ASU length
  :param refine_selection: (flex.siz_t) of all ncs related copies and
  non ncs related parts to be included in selection
  :return: (flex.siz_t) selection of all ncs groups master ncs selection and
  non ncs related portions that are being refined
  """
  refine_selection = set(refine_selection)
  total_master_ncs_selection = set()
  total_ncs_related_selection = set()
  for nrg in ncs_restraints_group_list:
    master_ncs_selection = nrg.master_iselection
    total_master_ncs_selection.update(set(master_ncs_selection))
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.copy_iselection
      total_ncs_related_selection.update(set(asu_selection))
  # make sure all ncs related parts are in refine_selection
  all_ncs = total_master_ncs_selection | total_ncs_related_selection
  msg = 'refine_selection does not contain all ncs related atoms'
  assert not bool(all_ncs - refine_selection), msg
  #
  extended_ncs_selection = refine_selection - total_ncs_related_selection
  return flex.size_t(list(extended_ncs_selection))

def get_ncs_related_selection(ncs_restraints_group_list,asu_size):
  """
  :param ncs_restraints_group_list: list of ncs_restraint_group objects
  total_asu_length: (int) Complete ASU length
  :param asu_size: (int) the total size of the ASU
  :return: (flex.bool) selection of all ncs related atom in the ASU
  """
  total_master_ncs_selection = set()
  total_ncs_related_selection = set()
  for nrg in ncs_restraints_group_list:
    master_ncs_selection = nrg.master_iselection
    total_master_ncs_selection.update(set(master_ncs_selection))
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.copy_iselection
      total_ncs_related_selection.update(set(asu_selection))
  #
  total_ncs_related_selection.update(total_master_ncs_selection)
  ts = flex.size_t(list(total_ncs_related_selection))
  selection = flex.bool(asu_size, ts)
  return selection

