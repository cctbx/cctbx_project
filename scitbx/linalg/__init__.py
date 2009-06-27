from scitbx.linalg.ext import *
from scitbx.array_family import flex

class random_normal_matrix_generator(ext.random_normal_matrix_generator):

  def __init__(self, m, n):
    super(random_normal_matrix_generator, self).__init__(m, n)
    self.state = flex.random_generator.getstate()
