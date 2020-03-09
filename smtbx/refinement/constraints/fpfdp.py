from __future__ import absolute_import, division, print_function
import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
from six.moves import range

class shared_fp(object):
  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  @property
  def constrained_parameters(self):
    return tuple((idx, "fp") for idx in self.indices[1:])

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    ref = reparametrisation.add_new_fp_parameter(self.indices[0])
    for i in range(1, len(self.indices)):
      param = reparametrisation.add(
        _.shared_fp,
        reference=ref,
        scatterer = scatterers[self.indices[i]])
      reparametrisation.shared_fps[self.indices[i]] = ref
      reparametrisation.asu_scatterer_parameters[self.indices[i]].fp = param
    self.value = ref

class shared_fdp(object):
  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  @property
  def constrained_parameters(self):
    return tuple((idx, "fdp") for idx in self.indices[1:])

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    ref = reparametrisation.add_new_fdp_parameter(self.indices[0])
    for i in range(1, len(self.indices)):
      param = reparametrisation.add(
        _.shared_fdp,
        reference=ref,
        scatterer = scatterers[self.indices[i]])
      reparametrisation.shared_fps[self.indices[i]] = ref
      reparametrisation.asu_scatterer_parameters[self.indices[i]].fdp = param
    self.value = ref
