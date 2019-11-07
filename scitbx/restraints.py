
"""
Handling for normalized target and gradients.
"""

from __future__ import absolute_import, division, print_function
from scitbx.stdlib import math
from libtbx.utils import Sorry

class energies(object):
  """
  Convenient wrapper class for encapsulating the output of any arbitrary
  restraining target function, such as the geometry restraints used elsewhere
  in CCTBX.  Tracks the number of restraints to facilitate normalized targets.
  """
  def __init__(O,
        compute_gradients,
        gradients,
        gradients_size,
        gradients_factory,
        normalization):
    O.number_of_restraints = 0
    O.residual_sum = 0
    if (not compute_gradients):
      O.gradients = None
    else:
      if (gradients is None):
        assert gradients_size is not None
        O.gradients = gradients_factory(gradients_size)
      else:
        if (gradients_size is not None):
          assert gradients.size() == gradients_size
        O.gradients = gradients
    O.normalization = normalization
    O.normalization_factor = None

  def __iadd__(O, rhs):
    O.number_of_restraints += rhs.number_of_restraints
    O.residual_sum += rhs.residual_sum
    if (O.gradients is not None):
      assert rhs.gradients is not None
      if (O.gradients.id() != rhs.gradients.id()):
        O.gradients += rhs.gradients
    return O

  def __setattr1__(self, attr, value):
    if attr=="target":
      if math.isnan(value):
        raise Sorry('target value (energy) is "nan"')
    object.__setattr__(self, attr, value)

  def finalize_target_and_gradients(O):
    O.target = O.residual_sum
    if (O.normalization):
      O.normalization_factor = 1.0 / max(1, O.number_of_restraints)
      O.target *= O.normalization_factor
      if (O.gradients is not None):
        O.gradients *= O.normalization_factor
