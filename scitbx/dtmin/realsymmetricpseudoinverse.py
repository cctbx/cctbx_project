from __future__ import division, print_function
from scitbx.array_family import flex
from scitbx.linalg import eigensystem
from sys import float_info

def RealSymmetricPseudoInverse(A_, filtered, min_to_filter=0):
  """Compute the pseudo-inverse of a real symmetric matrix

  filters (i.e. makes 0) small and negative eigen values before
  calculating the pseudo-inverse.
  """

  A = A_.deep_copy() # make a copy to avoid changing the original when scaling

  #A is a flex array, size N by N
  (M, N) = A.all()
  assert(M == N)
  A_scale = []
  for i in range(N):
    assert(A[i,i] > 0)
    A_scale.append(A[i,i] ** -0.5)
  for i in range(N):
    for j in range(N):
      A[i,j] = A_scale[i]*A_scale[j]*A[i,j]

  # Call the scitbx routine to find eigenvalues and eigenvectors
  eigen = eigensystem.real_symmetric(A)

  # Set condition number for range of eigenvalues.
  evCondition = min(1.e9, 0.01/float_info.epsilon)
  evMax = eigen.values()[0]       # (Eigen values are stored in descending order.)
  evMin = evMax/evCondition
  ZERO = 0.
  lambdaInv = flex.double(N*N, ZERO)
  lambdaInv.reshape(flex.grid(N,N)) # Matrix of zeros

  filtered = 0
  for i in range(N-1,-1, -1):
    if (eigen.values()[i] < evMin) or (filtered < min_to_filter):
      lambdaInv[i,i] = ZERO
      filtered += 1
    else:
      lambdaInv[i,i] = 1.0 / eigen.values()[i]

  # Compute the result that we are after - the pseudo inverse
  # The formula from wikipedia is (A^-1) = (Q)(lam^-1)(Q^T)
  # where Q has eigen vectors as columns
  # eigen.vectors() has eigenvectors as rows
  # let V be eigen.vectors(). The formula is now: (A^-1) = (V^T)(lam^-1)(V)

  # get a transposed matrix
  vectors_transposed = eigen.vectors().deep_copy()
  vectors_transposed.matrix_transpose_in_place()
  inverse_matrix = vectors_transposed.matrix_multiply( lambdaInv.matrix_multiply(eigen.vectors()) )

  # Rescale the matrix
  for i in range(N):
    for j in range(N):
      inverse_matrix[i,j] = A_scale[i]*A_scale[j]*inverse_matrix[i,j]

  return (inverse_matrix, filtered)

# nb here the c++ version has a debug section turned on by a #define.
# this does not.
