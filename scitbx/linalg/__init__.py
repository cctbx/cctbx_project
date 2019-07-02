from __future__ import absolute_import, division, print_function
from scitbx.linalg.ext import *
from scitbx.linalg.householder import *
from scitbx.array_family import flex
import scitbx.matrix

class random_normal_matrix_generator(ext.random_normal_matrix_generator):

  def __init__(self, m, n):
    super(random_normal_matrix_generator, self).__init__(m, n)
    self.state = flex.random_generator.getstate()


def matrix_equality_ratio(a, b, eps=None):
  if isinstance(a, scitbx.matrix.rec):
    a = a.as_flex_double_matrix()
  if isinstance(b, scitbx.matrix.rec):
    b = b.as_flex_double_matrix()
  if eps is not None:
    return ext.matrix_equality_ratio(a, b, eps)
  else:
    return ext.matrix_equality_ratio(a, b)
