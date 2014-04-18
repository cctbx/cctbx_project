from __future__ import division
from scitbx.array_family import flex
import scitbx.rigid_body
import  math

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
