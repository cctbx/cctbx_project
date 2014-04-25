from __future__ import division
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from mmtbx import monomer_library
from scitbx import matrix
import scitbx.rigid_body
from cctbx import xray
import random
import math


def concatenate_rot_tran(rot,tran,s=1):
    """
    Concatenate rotation angles, corresponding to the rotation
    matrices and scaled translation vectors to a single long flex.double object

    Arguments:
    rot : (lists of objects matrix.rec) Rotation matrices
    tran : (lists of objects matrix.rec) Tanslation vectors
    s : (float) scaling factor, scale translation to be of the
        same order ot the translation

    Return:
    flex.double object of the form
    [(alpha_1,beta_1,gamma_1,Tx_1/s,Ty_1/s,Tz_1/s)...]
    """
    x = []
    [x.extend(list(rotation_to_angles(r.elems)) + list((t/s).elems))
     for (r,t) in zip(rot,tran)]
    assert len(x) == 6*(len(rot))
    assert len(x) != 0
    return flex.double(x)

def separate_rot_tran(x,ncs_copies,s=1.0):
  """
  Convert the refinemable parameters, rotations angles and
  scaled translations, back to rotation matrices and translation vectors

  Arguments:
  x : a flex.double of the form (theta_1,psi_1,phi_1,tx_1,ty_1,tz_1,..
      theta_n,psi_n,phi_n,tx_n/s,ty_n/s,tz_n/s). where n is the number of
      transformations.
  ncs_copies : (int) expected number of ncs copies
  s : (float) scaling factor, scale translation to be of the
      same order ot the translation

  Returns:
  rot : matrix.rec type, size (3,3)
  tran: matrix.rec type, size (3,1)
  """
  rot = []
  tran = []

  for i in range(ncs_copies):
    the,psi,phi =x[i*6:i*6+3]
    rot_obj = scitbx.rigid_body.rb_mat_xyz(
      the=the, psi=psi, phi=phi, deg=False)
    rot.append(rot_obj.rot_mat())
    tran.append(matrix.rec(x[i*6+3:i*6+6],(3,1))*s)

  assert len(rot) == ncs_copies
  assert len(tran) == ncs_copies
  return rot,tran

def rotation_to_angles(rotation, deg=False):
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
    # using atan2 and dividing by cos(beta) to take into account the possible
    # different angles and signs
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
  return angles

def angles_to_rotation(angles_xyz, deg=False, rotation_is_tuple=False):
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
    return R.round(8).elems
  else:
    return flex.double(R.round(8))

def shake_transformations(rotation_matrices,
                          translation_vectors,
                          shake_angles_sigma      = 0.035,
                          shake_translation_sigma = 0.5,
                          return_all_transforms   = False):
  """
  Shake rotation matrices and translation vectors of a rotation matrices and
  translation vectors from the MTRIX records in a PDB file.

  Argument:
  rotation_matrices: (lists of objects matrix.rec)
  translation_vectors:(lists of objects matrix.rec)
  shake_angles_sigma: (float) the sigma (in radians) of the random gaussian
                      shaking of the rotation angles
  shake_translation_sigma: (float) the sigma (in angstrom) of the random
                           gaussian shaking of the translation
  return_all_transforms: (bool) When False, returns only shaken transformations

  Return:
  Shaken rotation matrices and translation vectors
  """
  rt = []
  tv = []
  iter_obj = enumerate(zip(rotation_matrices,translation_vectors))
  for i,(r,t) in iter_obj:
    if not (r.is_r3_identity_matrix() and t.elems == (0,0,0)):
      # shake rotations angles (alpha,beta,gamma)
      angles = rotation_to_angles(r.elems,deg=False)
      new_angles = [random.gauss(x,shake_angles_sigma) for x in angles]
      new_r_elems = angles_to_rotation(
        angles_xyz=new_angles,rotation_is_tuple=True)
      rt.append(matrix.sqr(new_r_elems))
      rotation_matrices[i] = matrix.sqr(new_r_elems)
      # Shake translation vectors
      new_t_elems = [random.gauss(x,shake_translation_sigma) for x in t.elems]
      tv.append(matrix.rec(new_t_elems,(3,1)))
      translation_vectors[i] = matrix.rec(new_t_elems,(3,1))

  if return_all_transforms:
    return rotation_matrices,translation_vectors
  else:
    return rt,tv

def compute_transform_grad(grad_wrt_xyz,
                           rotation_matrices,
                           xyz_ncs,
                           x):
  """
  Compute gradient in respect to the rotation angles and the translation
  vectors. R = Rx(the)Ry(psi)Rz(phi)

  Arguments:
  grad_wrt_xyz : (flex.double) gradients with respect to xyz.
  rotation_matrices: (list of matrix objects) rotation matrices
  xyz_ncs: (flex.vec3) The coordinates sites cart of the single NCS copy
  x: (flex double) Of the form (x1,y1,z1,...,xn,yn,zn)

  Returns:
  g : (flex.double) a gradient
  """
  g = []
  number_of_ncs_copies = len(rotation_matrices) + 1
  n_grad_in_ncs = len(grad_wrt_xyz)//number_of_ncs_copies
  assert len(xyz_ncs.as_double()) == n_grad_in_ncs
  # collect the derivative of all transformations
  for i in range(number_of_ncs_copies - 1):
    # get the portion of the gradient corresponding to the n'th NCS copy
    # not looking at the gradient of the original NCS copy
    grad_ncs_wrt_xyz = grad_wrt_xyz[(i+1)*n_grad_in_ncs:(i+2)*n_grad_in_ncs]
    # Translation derivatives are the same for all transformations
    grad_wrt_t = list(flex.vec3_double(grad_ncs_wrt_xyz).sum())
    # Sum angles gradient over the coordinates
    m = flex.vec3_double(grad_ncs_wrt_xyz).transpose_multiply(xyz_ncs)
    m = matrix.sqr(m)
    # Calculate gradient with respect to the rotation angles
    the,psi,phi = x[i*6:i*6+3]
    rot = scitbx.rigid_body.rb_mat_xyz(
      the=the, psi=psi, phi=phi, deg=False)
    g_the = (m*rot.r_the().transpose()).trace()
    g_psi = (m*rot.r_psi().transpose()).trace()
    g_phi = (m*rot.r_phi().transpose()).trace()
    g.extend([g_the, g_psi, g_phi])
    g.extend(grad_wrt_t)
  assert len(g) == 6*(number_of_ncs_copies - 1)
  return flex.double(g)

def get_ncs_sites_cart(ncs_obj):
  """
  Argument:
  ncs_obj: an object that contains fmodel and atom selection of a single NCS
  copy

  Return: (flex.vec3) coordinate sites cart
  """
  xrs_one_ncs = ncs_obj.fmodel.xray_structure.select(ncs_obj.ncs_atom_selection)
  return xrs_one_ncs.sites_cart()

def get_weight(minimization_obj):
  """
  Calculates weights for refinements

  Minimization object must contain the following methods and attributes:
  fmodel
  sites, u_iso, transformations: (bool)
  rotations, translations: (matrix objects)
  grm: Restraints manager
  iso_restraints: (libtbx.phil.scope_extract)
                  object used for u_iso refinement parameters
  ncs_atom_selection: (flex bool) for a single ncs atom selection

  Return:
  weight: (int)
  """
  mo = minimization_obj
  assert [mo.sites,mo.u_iso,mo.transformations].count(True)==1
  fmdc = mo.fmodel.deep_copy()
  if mo.sites:
    fmdc.xray_structure.shake_sites_in_place(mean_distance=0.3)
  elif mo.u_iso:
    fmdc.xray_structure.shake_adp()
  elif mo.transformations:
    rotation_matrices,translation_vectors = shake_transformations(
      rotation_matrices = mo.rotations,
      translation_vectors = mo.translations,
      shake_angles_sigma=0.01,
      shake_translation_sigma=0.1)
    x = concatenate_rot_tran(
      rotation_matrices,translation_vectors)
  fmdc.update_xray_structure(xray_structure = fmdc.xray_structure,
    update_f_calc=True)
  fmdc.xray_structure.scatterers().flags_set_grads(state=False)
  if mo.sites:
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # manager restraints, energy sites gradients
    gc = mo.grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
  elif mo.u_iso:
    # Create energies_site gradient, to create
    # geometry_restraints_manager.plain_pair_sym_table
    # needed for the energies_adp_iso
    gc = mo.grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      u_iso      = True)
     # fmodel gradients
    gxc = fmdc.one_time_gradients_wrt_atomic_parameters(
      u_iso = True).as_double()
    # manager restraints, energy sites gradients
    gc = mo.grm.energies_adp_iso(
      xray_structure    = fmdc.xray_structure,
      parameters        = mo.iso_restraints,
      use_u_local_only  = mo.iso_restraints.use_u_local_only,
      use_hd            = False,
      compute_gradients = True).gradients
  elif mo.transformations:
    xyz_ncs = get_ncs_sites_cart(mo)
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc_xyz = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # manager restraints, energy sites gradients
    gc_xyz = mo.grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
    gxc = compute_transform_grad(
      grad_wrt_xyz      = gxc_xyz.as_double(),
      rotation_matrices = rotation_matrices,
      xyz_ncs           = xyz_ncs,
      x                 = x)
    gc = compute_transform_grad(
      grad_wrt_xyz      = gc_xyz.as_double(),
      rotation_matrices = rotation_matrices,
      xyz_ncs           = xyz_ncs,
      x                 = x)

  weight = 1.
  gc_norm  = gc.norm()
  gxc_norm = gxc.norm()
  if(gxc_norm != 0.0):
    weight = gc_norm / gxc_norm

  weight =min(weight,1e6)
  return weight

def get_restraints_manager(pdb_file_name=None,pdb_string=None):
  """
  Generate restraint manager from a PDB file or a PDB string
  """
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
