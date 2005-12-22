from cctbx import sgtbx

class group_search(object):

  def __init__(self, modulus=2):
    self.ext = sgtbx.lattice_symmetry_group_search(modulus)

  def n_potential_axes(self):
    return self.ext.n_potential_axes()

  def __call__(self,
        minimum_cell,
        max_delta=3.,
        enforce_max_delta_for_generated_two_folds=True):
    return self.ext(
      minimum_cell, max_delta, enforce_max_delta_for_generated_two_folds)

group = group_search()

def find_max_delta(minimum_cell, group):
  return sgtbx.lattice_symmetry_find_max_delta(minimum_cell, group)
