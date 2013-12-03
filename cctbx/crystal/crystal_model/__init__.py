from __future__ import division
from scitbx import matrix
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import space_group as SG
from cctbx.sgtbx import space_group_symbols
from cctbx.crystal_orientation import crystal_orientation

class crystal_model(object):
  '''Simple model for the crystal lattice geometry and symmetry

  A crystal is initialised from the elements of its real space axes
  a, b, and c. Space group information must also be provided, either
  in the form of a symbol, or an existing
  cctbx.sgtbx.space_group object. If space_group_symbol is provided,
  it is passed to the cctbx.sgtbx.space_group_symbols constructor.
  This accepts either extended Hermann Mauguin format, or Hall format
  with the prefix 'Hall:'. E.g.

  space_group_symbol = "P b a n:1"
      or
  space_group_symbol = "Hall:P 2 2 -1ab"

  Optionally the crystal mosaicity value may be set, with the deg
  parameter controlling whether this value is treated as being an
  angle in degrees or radians.'''

  def __init__(self, real_space_a, real_space_b, real_space_c,
               space_group_symbol=None, space_group=None,
               mosaicity=None, deg=True):

    # Set the space group
    assert [space_group_symbol, space_group].count(None) == 1
    if space_group_symbol:
      self._sg = SG(space_group_symbols(space_group_symbol))
    else: self._sg = space_group

    # Set the mosaicity
    if mosaicity is not None:
      self.set_mosaicity(mosaicity, deg=deg)
    else:
      self._mosaicity = 0.0

    # setting matrix at initialisation
    real_space_a = matrix.col(real_space_a)
    real_space_b = matrix.col(real_space_b)
    real_space_c = matrix.col(real_space_c)
    A = matrix.sqr(real_space_a.elems +  real_space_b.elems + \
                   real_space_c.elems).inverse()

    # unit cell
    self.set_unit_cell(real_space_a, real_space_b, real_space_c)

    # reciprocal space orthogonalisation matrix (is the transpose of the
    # real space fractionalisation matrix, see http://goo.gl/H3p1s)
    self._update_B()

    # initial orientation matrix
    self._U = A * self._B.inverse()

  def __str__(self):
    uc = self.get_unit_cell().parameters()
    sg = str(self.get_space_group().info())
    umat = self.get_U().mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    bmat = self.get_B().mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    amat = (self.get_U() * self.get_B()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()

    msg =  "Crystal:\n"
    msg += "    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc + "\n"
    msg += "    Space group: " + sg + "\n"
    msg += "    U matrix:  " + umat[0] + "\n"
    msg += "               " + umat[1] + "\n"
    msg += "               " + umat[2] + "\n"
    msg += "    B matrix:  " + bmat[0] + "\n"
    msg += "               " + bmat[1] + "\n"
    msg += "               " + bmat[2] + "\n"
    msg += "    A = UB:    " + amat[0] + "\n"
    msg += "               " + amat[1] + "\n"
    msg += "               " + amat[2] + "\n"
    return msg

  def set_unit_cell(self, real_space_a, real_space_b, real_space_c):
    cell = (real_space_a.length(),
            real_space_b.length(),
            real_space_c.length(),
            real_space_b.angle(real_space_c, deg = True),
            real_space_c.angle(real_space_a, deg = True),
            real_space_a.angle(real_space_b, deg = True))
    self._uc = unit_cell(cell)
    self._update_B()

  def _update_B(self):
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

  def set_U(self, U):

    # check U is a rotation matrix.
    assert(U.is_r3_rotation_matrix())
    self._U = U

  def get_U(self):
    return self._U

  def get_B(self):
    return self._B

  def set_B(self, B):

    # also set the unit cell
    co = crystal_orientation(B,True)
    self._uc = co.unit_cell()
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

  def get_unit_cell(self):
    return self._uc

  def get_space_group(self):
    return self._sg

  def get_mosaicity(self, deg=True):
    from math import pi
    if deg == True:
      return self._mosaicity * 180.0 / pi

    return self._mosaicity

  def set_mosaicity(self, mosaicity, deg=True):
    from math import pi
    if deg == True:
      self._mosaicity = mosaicity * pi / 180.0
    else:
      self._mosaicity = mosaicity

  def get_A(self):
    return self._U * self._B

  def __eq__(self, other, eps=1e-7):
    if isinstance(other, crystal_model):
      d_mosaicity = abs(self._mosaicity - other._mosaicity)
      d_U = sum([abs(u1 - u2) for u1, u2 in zip(self._U, other._U)])
      d_B = sum([abs(b1 - b2) for b1, b2 in zip(self._B, other._B)])
      return (d_mosaicity <= eps and
              d_U <= eps and
              d_B <= eps and
              self._sg == other._sg)
    return NotImplemented

  def get_real_space_vectors(self):
    A_inv = self.get_A().inverse()
    return (matrix.col(A_inv[:3]),
            matrix.col(A_inv[3:6]),
            matrix.col(A_inv[6:9]))

  def change_basis(self, change_of_basis_op):
    direct_matrix = self.get_A().inverse()
    M = matrix.sqr(change_of_basis_op.c_inv().r().transpose().as_double())
    # equation 2.19 of Giacovazzo
    new_direct_matrix = M * direct_matrix
    real_space_a = new_direct_matrix[:3]
    real_space_b = new_direct_matrix[3:6]
    real_space_c = new_direct_matrix[6:9]
    other = crystal_model(real_space_a,
                          real_space_b,
                          real_space_c,
                          space_group=self.get_space_group().change_basis(
                            change_of_basis_op),
                          mosaicity=self.get_mosaicity())
    return other

def crystal_model_from_mosflm_matrix(mosflm_A_matrix,
                                     unit_cell=None,
                                     wavelength=None,
                                     space_group=None):
  '''Create a crystal_model from a Mosflm A matrix (a*, b*, c*); N.B. assumes
  the mosflm coordinate frame:

                                                   /!
                      Y-axis                      / !
                        ^                        /  !
                        !                       /   !
                        !                      /    !
                        !   /                 /  Xd !
                        !  /                 / * ^  !
                        ! /                  ! 3 !  !
                        !/      X-ray beam   !   !  !
                        /------------------------/--!---->X-axis
                       /                     !  / *1!
                    <-/-                     ! /    !
                     /  \+ve phi             ! Yd  /
                    /   /                    ! 2  /
                   /                         ! * /
                  Z-axis                  Ys ^ _/
                Rotation                     ! /| Xs
                 axis                        !/
                                             O

  Also assume that the mosaic spread is 0. If space_group is None spacegroup
  will be assigned as P1.
  '''

  if not space_group:
    space_group = SG('P1')

  A_star = matrix.sqr(mosflm_A_matrix)
  A = A_star.inverse()

  if unit_cell:
    unit_cell_constants = unit_cell.parameters()
    a = matrix.col(A.elems[0:3])
    wavelength = unit_cell_constants[0] / a.length()
    A *= wavelength
  elif wavelength:
    A *= wavelength
  else:
    # assume user has pre-scaled
    pass

  a = A.elems[0:3]
  b = A.elems[3:6]
  c = A.elems[6:9]
  rotate_mosflm_to_imgCIF = matrix.sqr((0, 0, 1, 0, 1, 0, -1, 0, 0))
  _a = rotate_mosflm_to_imgCIF * a
  _b = rotate_mosflm_to_imgCIF * b
  _c = rotate_mosflm_to_imgCIF * c

  return crystal_model(_a, _b, _c, space_group=space_group)
