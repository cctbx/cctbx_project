from __future__ import absolute_import, division, print_function

from scitbx import matrix
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import space_group as SG
from cctbx.sgtbx import space_group_symbols
from cctbx.crystal_orientation import crystal_orientation

class crystal_model_old(object):
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

    # set up attributes for scan-varying model
    self.reset_scan_points()

    # set up attributes for unit cell errors
    self.reset_unit_cell_errors()

    return

  @property
  def num_scan_points(self):
    return self._num_scan_points

  def show(self, show_scan_varying=False, out=None):
    if out is None:
      import sys
      out = sys.stdout
    uc = self.get_unit_cell().parameters()
    sg = str(self.get_space_group().info())
    umat = self.get_U().mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    bmat = self.get_B().mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    amat = (self.get_U() * self.get_B()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()

    msg =  ["Crystal:"]
    msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
    msg.append("    Space group: " + sg)
    msg.append("    U matrix:  " + umat[0])
    msg.append("               " + umat[1])
    msg.append("               " + umat[2])
    msg.append("    B matrix:  " + bmat[0])
    msg.append("               " + bmat[1])
    msg.append("               " + bmat[2])
    msg.append("    A = UB:    " + amat[0])
    msg.append("               " + amat[1])
    msg.append("               " + amat[2])
    if self.num_scan_points > 0:
      msg.append("    A sampled at " + str(self.num_scan_points) \
                 + " scan points")
      if show_scan_varying:
        for i in range(self.num_scan_points):
          A = self.get_A_at_scan_point(i)
          B = self.get_B_at_scan_point(i)
          U = self.get_U_at_scan_point(i)
          uc = self.get_unit_cell_at_scan_point(i).parameters()
          umat = U.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          bmat = B.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          amat = A.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          msg.append("  Scan point #%i:" %(i+1))
          msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
          msg.append("    U matrix:  " + umat[0])
          msg.append("               " + umat[1])
          msg.append("               " + umat[2])
          msg.append("    B matrix:  " + bmat[0])
          msg.append("               " + bmat[1])
          msg.append("               " + bmat[2])
          msg.append("    A = UB:    " + amat[0])
          msg.append("               " + amat[1])
          msg.append("               " + amat[2])
    print("\n".join(msg), file=out)

  def __str__(self):
    from cStringIO import StringIO
    s = StringIO()
    msg = self.show(out=s)
    s.seek(0)
    return s.read()

  def set_unit_cell(self, real_space_a, real_space_b, real_space_c):
    cell = (real_space_a.length(),
            real_space_b.length(),
            real_space_c.length(),
            real_space_b.angle(real_space_c, deg = True),
            real_space_c.angle(real_space_a, deg = True),
            real_space_a.angle(real_space_b, deg = True))
    self._uc = unit_cell(cell)
    self._update_B()

  def get_unit_cell_at_scan_point(self, t):
    B = self.get_B_at_scan_point(t)
    co = crystal_orientation(B)
    return co.unit_cell()

  def _update_B(self):
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

  def set_U(self, U):

    # check U is a rotation matrix.
    assert(U.is_r3_rotation_matrix())
    self._U = U

    # reset scan-varying data, if the static U has changed
    self.reset_scan_points()

  def get_U(self):
    return self._U

  def get_B(self):
    return self._B

  def set_B(self, B):

    # also set the unit cell
    co = crystal_orientation(B,True)
    self._uc = co.unit_cell()
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

    # reset scan-varying data, if the static B has changed
    self.reset_scan_points()

    # reset unit cell errors
    self.reset_unit_cell_errors()

  def get_B_covariance(self):
    '''Return the 9*9 covariance matrix of elements of B, if available,
    otherwise None. The order of elements of this matrix are determined by
    'flattening' B to form a 1D vector using row major ordering.

    :rtype: :py:class:`scitbx.matrix.sqr`
    '''
    return self._cov_B

  def set_B_covariance(self, cov):

    cov = matrix.sqr(cov)

    # check cov is of the right size. No other checks made
    assert cov.n == (9, 9)

    self._cov_B = cov

    return

  def _calc_cell_parameter_sd(self):

    from scitbx.math.lefebvre import matrix_inverse_error_propagation
    from scitbx.math import angle_derivative_wrt_vectors
    from math import pi, acos, sqrt

    # self._cov_B is the covariance matrix of elements of the B matrix. We
    # need to construct the covariance matrix of elements of the
    # transpose of B. The vector of elements of B is related to the
    # vector of elements of its transpose by a permutation, P.
    P = matrix.sqr((1,0,0,0,0,0,0,0,0,
                    0,0,0,1,0,0,0,0,0,
                    0,0,0,0,0,0,1,0,0,
                    0,1,0,0,0,0,0,0,0,
                    0,0,0,0,1,0,0,0,0,
                    0,0,0,0,0,0,0,1,0,
                    0,0,1,0,0,0,0,0,0,
                    0,0,0,0,0,1,0,0,0,
                    0,0,0,0,0,0,0,0,1))

    # We can use P to replace var_cov with the covariance matrix of
    # elements of the transpose of B.
    var_cov = P*self._cov_B*P

    # From B = (O^-1)^T we can convert this
    # to the covariance matrix of the real space orthogonalisation matrix
    Bt = self._B.transpose()
    O = Bt.inverse()
    cov_O = matrix_inverse_error_propagation(Bt, var_cov)

    # The real space unit cell vectors are given by
    vec_a = (O * matrix.col((1,0,0)))
    vec_b = (O * matrix.col((0,1,0)))
    vec_c = (O * matrix.col((0,0,1)))

    # So the unit cell parameters are
    a, b, c = vec_a.length(), vec_b.length(), vec_c.length()
    alpha = acos(vec_b.dot(vec_c) / (b*c))
    beta =  acos(vec_a.dot(vec_c) / (a*c))
    gamma = acos(vec_a.dot(vec_b) / (a*b))

    # The estimated errors are calculated by error propagation from cov_O. In
    # each case we define a function F(O) that converts the matrix O into the
    # unit cell parameter of interest. To do error propagation to get the
    # variance of that cell parameter we need the Jacobian of the function.
    # This is a 1*9 matrix of derivatives in the order of the elements of O
    #
    #     / dF   dF   dF   dF   dF   dF   dF   dF   dF  \
    # J = | ---, ---, ---, ---, ---, ---, ---, ---, --- |
    #     \ da1  db1  dc1  da2  db2  dc2  da3  db3  dc3 /
    #

    # For cell length |a|, F = sqrt(a1^2 + a2^2 + a3^2)
    jacobian = matrix.rec(
      (vec_a[0]/a, 0, 0, vec_a[1]/a, 0, 0, vec_a[2]/a, 0, 0), (1, 9))
    var_a = (jacobian * cov_O * jacobian.transpose())[0]

    # For cell length |b|, F = sqrt(b1^2 + b2^2 + b3^2)
    jacobian = matrix.rec(
      (0, vec_b[0]/b, 0, 0, vec_b[1]/b, 0, 0, vec_b[2]/b, 0), (1, 9))
    var_b = (jacobian * cov_O * jacobian.transpose())[0]

    # For cell length |c|, F = sqrt(c1^2 + c2^2 + c3^2)
    jacobian = matrix.rec(
      (0, 0, vec_c[0]/c, 0, 0, vec_c[1]/c, 0, 0, vec_c[2]/c), (1, 9))
    jacobian_t = jacobian.transpose()
    var_c = (jacobian * cov_O * jacobian.transpose())[0]

    # For cell volume (a X b).c,
    # F = c1(a2*b3 - b2*a3) + c2(a3*b1 - b3*a1) + c3(a1*b2 - b1*a2)
    a1, a2, a3 = vec_a
    b1, b2, b3 = vec_b
    c1, c2, c3 = vec_c
    jacobian = matrix.rec((c3*b2 - c2*b3,
                           c1*b3 - c3*b1,
                           c2*b1 - c1*b2,
                           c2*a3 - c3*a2,
                           c3*a1 - c1*a3,
                           c1*a2 - c2*a1,
                           a2*b3 - b2*a3,
                           a3*b1 - b3*a1,
                           a1*b2 - b1*a2), (1,9))
    jacobian_t = jacobian.transpose()
    var_V = (jacobian * cov_O * jacobian.transpose())[0]
    self._cell_volume_sd = sqrt(var_V)

    # For the unit cell angles we need to calculate derivatives of the angles
    # with respect to the elements of O
    dalpha_db, dalpha_dc = angle_derivative_wrt_vectors(vec_b, vec_c)
    dbeta_da, dbeta_dc = angle_derivative_wrt_vectors(vec_a, vec_c)
    dgamma_da, dgamma_db = angle_derivative_wrt_vectors(vec_a, vec_b)

    # For angle alpha, F = acos( b.c / |b||c|)
    jacobian = matrix.rec(
      (0, dalpha_db[0], dalpha_dc[0], 0, dalpha_db[1], dalpha_dc[1], 0,
      dalpha_db[2], dalpha_dc[2]), (1, 9))
    var_alpha = (jacobian * cov_O * jacobian.transpose())[0]

    # For angle beta, F = acos( a.c / |a||c|)
    jacobian = matrix.rec(
      (dbeta_da[0], 0, dbeta_dc[0], dbeta_da[1], 0, dbeta_dc[1], dbeta_da[2],
      0, dbeta_dc[2]), (1, 9))
    var_beta = (jacobian * cov_O * jacobian.transpose())[0]

    # For angle gamma, F = acos( a.b / |a||b|)
    jacobian = matrix.rec(
      (dgamma_da[0], dgamma_db[0], 0, dgamma_da[1], dgamma_db[1], 0,
      dgamma_da[2], dgamma_db[2], 0), (1, 9))
    var_gamma = (jacobian * cov_O * jacobian.transpose())[0]

    # Symmetry constraints may mean variances of the angles should be zero.
    # Floating point error could lead to negative variances. Ensure these are
    # caught before taking their square root!
    var_alpha = max(0, var_alpha)
    var_beta = max(0, var_beta)
    var_gamma = max(0, var_gamma)

    from math import pi
    rad2deg = 180. / pi
    self._cell_sd = (sqrt(var_a),
                     sqrt(var_b),
                     sqrt(var_c),
                     sqrt(var_alpha) * rad2deg,
                     sqrt(var_beta) * rad2deg,
                     sqrt(var_gamma) * rad2deg)
    return

  def get_cell_parameter_sd(self):
    '''Return the estimated standard deviations of unit cell parameters in
    units of Angstroms and degrees, if available, otherwise None.

    :rtype: :py:class:`tuple`
    '''
    if self._cov_B is None: return None
    if self._cell_sd is not None: return self._cell_sd

    # propagate covariance of B to errors in the unit cell parameters
    self._calc_cell_parameter_sd()

    return self._cell_sd

  def get_cell_volume_sd(self):
    '''Return the estimated standard deviation of the unit cell volume in
    units of Angstroms^3, if available, otherwise None.

    :rtype: :py:class:`tuple`
    '''
    if self._cov_B is None: return None
    if self._cell_volume_sd is not None: return self._cell_volume_sd

    # propagate covariance of B to errors in the unit cell parameters and
    # volume
    self._calc_cell_parameter_sd()

    return self._cell_volume_sd

  def set_A(self, A):
    uc = unit_cell(orthogonalization_matrix=A.transpose().inverse())
    B = matrix.sqr(uc.fractionalization_matrix()).transpose()
    U = A * B.inverse()
    self.set_B(B)
    self.set_U(U)

  def set_A_at_scan_points(self, A_list):
    '''Set the setting matrix A at a series of checkpoints within a rotation
    scan. There would typically be n+1 points, where n is the number of images
    in the scan. The first point is the state at the beginning of the rotation
    scan. The final point is the state at the end of the rotation scan.
    Intervening points give the state at the start of the rotation at the 2nd
    to the nth image.

    This data is held separately from the 'static' U and B because per-image
    setting matrices may be refined whilst restraining to a previously
    determined best-fit static UB. The values will be reset if any changes are
    made to the static U and B matrices.
    '''

    self._A_at_scan_points = [matrix.sqr(e) for e in A_list]
    self._num_scan_points = len(A_list)

  def get_A_at_scan_point(self, t):
    '''Return the setting matrix with index t.

    This will typically have been
    set with reference to a particular scan, such that it equals the UB matrix
    appropriate at the start of the rotation for the image with array index t

    :param t: index into the list of scan points.

    :returns: The A matrix at index t
    :rtype: :py:class:`scitbx.matrix.sqr`
    '''

    return self._A_at_scan_points[t]

  def get_U_at_scan_point(self, t):
    '''
    Return orientation matrix with index t.

    :param t: index into the list of scan points.

    :returns: The U matrix at index t
    :rtype: :py:class:`scitbx.matrix.sqr`
    '''

    Bt = self.get_B_at_scan_point(t)
    At = self._A_at_scan_points[t]

    return At * Bt.inverse()

  def get_B_at_scan_point(self, t):
    '''
    Return orthogonalisation matrix with index t.

    :param t: index into the list of scan points.

    :returns: The B matrix at index t
    :rtype: :py:class:`scitbx.matrix.sqr`
    '''

    At = self._A_at_scan_points[t]
    uc = unit_cell(orthogonalization_matrix=At.transpose().inverse())

    return matrix.sqr(uc.fractionalization_matrix()).transpose()

  def get_unit_cell_at_scan_point(self, t):
    '''
    Return unit cell with index t.

    :param t: index into the list of scan points.
    :type t: int

    :returns: The unit cell at index t
    :rtype: :py:class:`cctbx.uctbx.unit_cell`
    '''

    At = self._A_at_scan_points[t]
    uc = unit_cell(orthogonalization_matrix=At.transpose().inverse())

    return uc

  def reset_scan_points(self):
    self._num_scan_points = 0
    self._A_at_scan_points = None

  def reset_unit_cell_errors(self):
    self._cov_B = None
    self._cell_sd = None
    self._cell_volume_sd = None

  def get_unit_cell(self):
    '''
    :returns: The unit cell
    :rtype: :py:class:`cctbx.uctbx.unit_cell`
    '''
    return self._uc

  def set_space_group(self, space_group):
    '''
    :param space_group:
    :type space_group: cctbx.sgtbx.space_group
    '''
    self._sg = space_group

  def get_space_group(self):
    '''
    :returns: The current space group
    :rtype: :py:class:`cctbx.sgtbx.space_group`
    '''
    return self._sg

  def get_mosaicity(self, deg=True):
    '''
    Get the current value of the mosaicity.

    :param deg: If True, return the mosaicity in degrees, else return the
                mosaicity in radians.
    :type deg: bool
    :returns: The mosaicity
    :rtype: float or None
    '''
    from math import pi
    if deg == True:
      return self._mosaicity * 180.0 / pi

    return self._mosaicity

  def set_mosaicity(self, mosaicity, deg=True):
    '''
    Update the value of the mosaicity.

    :param mosaicity: The mosaicity
    :type mosaicity: float
    :param deg: If True, assume the mosaicity is given in degrees, else
                assume it is given in radians.
    :type deg: bool
    :returns: The mosaicity
    :rtype: float or None
    '''
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
      if self.num_scan_points > 0:
        if other.num_scan_points != self.num_scan_points: return False
        for i in range(self.num_scan_points):
          A1, A2 = self.get_A_at_scan_point(i), other.get_A_at_scan_point(i)
          d_A = sum([abs(a1 - a2) for a1, a2 in zip(A1, A2)])
          if d_A > eps: return False
      return (d_mosaicity <= eps and
              d_U <= eps and
              d_B <= eps and
              self._sg == other._sg)
    return NotImplemented

  def is_similar_to(self, other, angle_tolerance=0.01,
                    mosaicity_tolerance=0.8, uc_rel_length_tolerance=0.01,
                    uc_abs_angle_tolerance=1.):
    '''
    Test similarity of this to another crystal model

    :param other: the crystal model to test against
    :type other: crystal_model
    :param angle_tolerance: maximum tolerated orientation difference in degrees
    :type angle_tolerance: float
    :param mosaicity_tolerance: minimum tolerated fraction of the larger
            mosaicity for the smaller mosaicity to retain similarity
    :type mosaicity_tolerance:float
    :param uc_rel_length_tolerance: relative length tolerance to pass to
            unit_cell.is_similar_to
    :type uc_rel_length_tolerance: float
    :param uc_abs_angle_tolerance: absolute angle tolerance to pass to
            unit_cell.is_similar_to
    :type uc_abs_angle_tolerance: float
    :returns: Whether the other crystal model is similar to this
    :rtype: bool
    '''

    # space group test
    if self.get_space_group() != other.get_space_group(): return False

    # mosaicity test
    m_a, m_b = self.get_mosaicity(deg=True), other.get_mosaicity(deg=True)
    min_m, max_m = min(m_a, m_b), max(m_a, m_b)
    if min_m <= 0.0:
      if max_m > 0.0: return False
    elif min_m / max_m < mosaicity_tolerance: return False

    # static orientation test
    U_a, U_b = self.get_U(), other.get_U()
    assert U_a.is_r3_rotation_matrix()
    assert U_b.is_r3_rotation_matrix()
    R_ab = U_b * U_a.transpose()
    uq = R_ab.r3_rotation_matrix_as_unit_quaternion()
    angle = uq.unit_quaternion_as_axis_and_angle(deg=True)[0]
    if abs(angle) > angle_tolerance: return False

    # static unit cell test
    uc_a, uc_b = self.get_unit_cell(), other.get_unit_cell()
    if not uc_a.is_similar_to(uc_b,
      relative_length_tolerance=uc_rel_length_tolerance,
      absolute_angle_tolerance=uc_abs_angle_tolerance): return False

    # scan varying tests
    if self.num_scan_points > 0:
      if other.num_scan_points != self.num_scan_points: return False
      for i in range(self.num_scan_points):
        U_a, U_b = self.get_U_at_scan_point(i), other.get_U_at_scan_point(i)
        assert U_a.is_r3_rotation_matrix()
        assert U_b.is_r3_rotation_matrix()
        R_ab = U_b * U_a.transpose()
        uq = R_ab.r3_rotation_matrix_as_unit_quaternion()
        angle = uq.unit_quaternion_as_axis_and_angle(deg=True)[0]
        if abs(angle) > angle_tolerance: return False

        B_a, B_b = self.get_B_at_scan_point(i), other.get_B_at_scan_point(i)
        uc_a = crystal_orientation(B_a,True).unit_cell()
        uc_b = crystal_orientation(B_b,True).unit_cell()
        if not uc_a.is_similar_to(uc_b,
          relative_length_tolerance=uc_rel_length_tolerance,
          absolute_angle_tolerance=uc_abs_angle_tolerance): return False

    return True

  def get_real_space_vectors(self):
    '''
    Get the real space unit cell basis vectors.

    :returns: Real space unit cell basis vectors
    :rtype: list of :py:class:`scitbx.matrix.col`
    '''
    A_inv = self.get_A().inverse()
    return (matrix.col(A_inv[:3]),
            matrix.col(A_inv[3:6]),
            matrix.col(A_inv[6:9]))

  def change_basis(self, change_of_basis_op):
    '''
    Returns a copy of the current crystal model transformed by the given
    change of basis operator to the new basis.

    :param change_of_basis_op: The change of basis operator.
    :type cctbx.sgtbx.change_of_basis_op:
    :returns: The crystal model transformed to the new basis.
    :rtype: :py:class:`crystal_model`
    '''
    # cctbx change of basis matrices and those Giacovazzo are related by
    # inverse and transpose, i.e. Giacovazzo's "M" is related to the cctbx
    # cb_op as follows:
    #   M = cb_op.c_inv().r().transpose()
    #   M_inverse = cb_op_to_minimum.c().r().transpose()

    # (Giacovazzo calls the direct matrix "A",
    #  we call the reciprocal matrix "A")
    # Therefore, from equation 2.19 in Giacovazzo:
    #   A' = M A

    # and:
    #   (A')^-1 = (M A)^-1
    #   (A')^-1 = A^-1 M^-1

    #reciprocal_matrix = self.get_A()
    #rm_cb = reciprocal_matrix * M.inverse()
    #dm_cb = rm_cb.inverse()
    #from libtbx.test_utils import approx_equal
    #assert approx_equal(dm_cb.elems, new_direct_matrix.elems)

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
    if self.num_scan_points > 0:
      M_inv = M.inverse()
      other.set_A_at_scan_points(
        [At * M_inv for At in self._A_at_scan_points])
      assert other.num_scan_points == self.num_scan_points
    return other

  def update(self, other):
    '''
    Update the current crystal model such that self == other.

    :param other: A crystal model to update self with.
    :type other: crystal_model
    '''
    self._U = other._U
    self._B = other._B
    self._A_at_scan_points = other._A_at_scan_points
    self._num_scan_points = other._num_scan_points
    self._sg = other._sg
    self._uc = other._uc
    self._mosaicity = other._mosaicity

  def rotate_around_origin(self, axis, angle, deg=True):
    '''
    Rotate the model around an axis and angle

    :param axis: The axis to rotate around
    :param angle: The angle to rotate around
    :param deg: Degrees or radians

    '''
    from scitbx import matrix

    # Compute the matrix
    R = matrix.col(axis).axis_and_angle_as_r3_rotation_matrix(angle, deg=deg)

    # Update U
    self._U = R * self._U

    # Update A at scan points
    if self._A_at_scan_points is not None:
      for i, At in enumerate(self._A_at_scan_points):
        uc = unit_cell(orthogonalization_matrix=At.transpose().inverse())
        Bt = matrix.sqr(uc.fractionalization_matrix()).transpose()
        Ut = At * Bt.inverse()
        Ut = R * Ut
        At = Ut * Bt
        self._A_at_scan_points[i] = At

