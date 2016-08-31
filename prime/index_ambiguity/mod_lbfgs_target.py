from __future__ import division
from cctbx.array_family import flex
import numpy as np
"""
lbfgs_target_handler:
Calculate indexing ambiguity resolving target function according to Brehm & Diederics, 2014.
"""
class lbfgs_target_handler(object):
  def __init__(self):
    """
    Intialitze parameters
    """

  def func(self, x, r_matrix):
    # reshape x
    x_set = np.array(x).reshape((len(x)/2,2))
    n_x = len(x_set)
    # calculate sum_n-1(sum_n((r_ij - dot(x_i, x_j))^2))
    error = flex.double([sum([abs(r_matrix[i,j] - np.dot(x_set[i], x_set[j])) for j in range(i+1,len(x_set))]) for i in range(len(x_set)-1)])
    return error
