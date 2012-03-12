from scitbx.array_family import flex
import scitbx.matrix


def savitzky_golay_coefficients(n_left, n_right, degree, derivative=0):
  """
    Compute the convolution coefficients to be used for smoothing data by the
    method of Savitzky and Golay [1].
    Parameters
    ----------
    n_left and n_right are the number of past and future data points used. The
    total number of data points used is n_left + n_right + 1.
    degree is the degree or order of the smoothing polynomial.
    Notes
    -----
    See also:
    http://www.scipy.org/Cookbook/SavitzkyGolay
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
  # the highest conserved moment is equal to the degree of the smoothing polynomial
  assert derivative <= degree
  a = flex.double([[k**i for i in range(degree+1)] for k in range(-n_left, n_right+1)])
  b = flex.double(degree+1, 0)
  b[derivative] = 1 # choose which derivative we want
  lu = a.matrix_transpose_multiply(a)
  pivot_indices = lu.matrix_lu_decomposition_in_place()
  x = lu.matrix_lu_back_substitution(pivot_indices=pivot_indices, b=b)
  coefficients = flex.double(n_left+n_right+1, 0)
  for i, k in enumerate(range(-n_left, n_right+1)):
    c = scitbx.matrix.col([k**n for n in range(degree+1)]).dot(
      scitbx.matrix.row(x))
    coefficients[i] = c
  return coefficients



def exercise_savitzky_golay_coefficients():
  from libtbx.test_utils import approx_equal
  coeffs = savitzky_golay_coefficients(5, 5, 4)
  assert approx_equal(
    coeffs,
    (0.042, -0.105, -0.023, 0.140, 0.280, 0.333, 0.280, 0.140, -0.023, -0.105, 0.042), eps=1e-3)

  coeffs = savitzky_golay_coefficients(4, 4, 4)
  assert approx_equal(
    coeffs,
    (0.035, -0.128, 0.070, 0.315, 0.417, 0.315, 0.070, -0.128, 0.035), eps=1e-3)

  coeffs = savitzky_golay_coefficients(4, 0, 2)
  assert approx_equal(
    coeffs,
    (0.086, -0.143, -0.086, 0.257, 0.886), eps=1e-3)


if __name__ == '__main__':
  exercise_savitzky_golay_coefficients()
  print "OK"
