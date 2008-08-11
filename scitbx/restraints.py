class energies(object):

  def __init__(self,
        compute_gradients,
        gradients,
        gradients_size,
        gradients_factory,
        normalization):
    self.number_of_restraints = 0
    self.residual_sum = 0
    if (not compute_gradients):
      self.gradients = None
    else:
      if (gradients is None):
        assert gradients_size is not None
        self.gradients = gradients_factory(gradients_size)
      else:
        if (gradients_size is not None):
          assert gradients.size() == gradients_size
        self.gradients = gradients
    self.normalization = normalization

  def __iadd__(self, rhs):
    self.number_of_restraints += rhs.number_of_restraints
    self.residual_sum += rhs.residual_sum
    if (self.gradients is not None):
      assert rhs.gradients is not None
      if (self.gradients.id() != rhs.gradients.id()):
        self.gradients += rhs.gradients
    return self

  def finalize_target_and_gradients(self):
    self.target = self.residual_sum
    if (self.normalization):
      normalization_factor = 1.0 / max(1, self.number_of_restraints)
      self.target *= normalization_factor
      if (self.gradients is not None):
        self.gradients *= normalization_factor
