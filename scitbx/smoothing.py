from scitbx.array_family import flex
import scitbx.matrix


def savitzky_golay_coefficients(n_left, n_right, degree, derivative=0, wraparound=True):
  """
    Compute the convolution coefficients to be used for smoothing data by the
    method of Savitzky and Golay [1].
    Parameters
    ----------
    n_left and n_right are the number of past and future data points used. The
    total number of data points used is n_left + n_right + 1.
    degree is the degree or order of the smoothing polynomial.
    wraparound: if True then the coefficients are returned in "wraparound" order
    i.e. point 0 is at index[0], 1 at index[1], etc. and point -1 is at index[-1],
    point -2 at index[-2], etc.
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
    np = n_left+n_right+1
    if wraparound:
      coefficients[(np - k) %np] = c
    else:
      coefficients[i] = c
  return coefficients


def savitzky_golay_filter(x, y, half_window, degree, derivative=0):
  from scitbx import smoothing
  #y = y.as_numpy_array()
  # pad the signal at the extremes with
  # values taken from the signal itself
  firstvals = y[1:half_window+1].reversed()
  lastvals = y[-half_window-1:-1].reversed()
  firstvals.extend(y)
  firstvals.extend(lastvals)
  y = firstvals
  # discrete convolution
  coeffs = smoothing.savitzky_golay_coefficients(
    half_window, half_window, degree, derivative=derivative, wraparound=False)
  x, y = x, smoothing.convolve(y, coeffs)[half_window:-half_window]
  y = y[half_window:-half_window]
  return x, y


def convolve(x, y, mode="full"):
  assert mode in ("full", "same", "valid")
  P, Q, N = len(x), len(y), len(x)+len(y)-1
  z = []
  for k in range(N):
    t, lower, upper = 0, max(0, k-(Q-1)), min(P-1, k)
    for i in range(lower, upper+1):
      t = t + x[i] * y[k-i]
    z.append(t)
  z = flex.double(z)
  if mode == "full":
    return flex.double(z)
  elif mode == "same":
    padding = (N - P)//2
    if (N - P) % 2 == 0:
      return z[padding:-padding]
    else:
      return z[padding:-padding-1]
  elif mode == "valid":
    padding = N - P
    return z[padding:-padding]
