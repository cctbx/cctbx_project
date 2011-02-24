import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint

class dependent_occupancy(object):
  """ occupancy of a site depend on the occupancy of the other site
  """

  def __init__(self, var_refs, var_minus_one_refs):
    if (len(var_refs) + len(var_minus_one_refs)) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.as_var = var_refs
    self.as_one_minus_var = var_minus_one_refs

  def __eq__(self, other):
    if self.indices != other.indices:  return False
    return True

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    if len(self.as_var) != 0:
      sc_idx = self.as_var[0][0]
      ratio = 1./self.as_var[0][1]
    else:
      sc_idx = self.as_one_minus_var[0][0]
      ratio = 1./self.as_one_minus_var[0][1]
    occupancy = reparametrisation.add_new_occupancy_parameter(sc_idx)
    for sc in self.as_var:
      if sc[0] == sc_idx:  continue
      param = reparametrisation.add(
        _.dependent_occupancy,
        occupancy = occupancy,
        multiplier = sc[1]*ratio,
        as_one = True,
        scatterer = reparametrisation.structure.scatterers()[sc[0]])
      reparametrisation.asu_scatterer_parameters[sc[0]].occupancy = param
      reparametrisation.shared_occupancies[sc[0]] = occupancy
    as_one = len(self.as_var) == 0  # only if both lists are not empty
    for sc in self.as_one_minus_var:
      if sc[0] == sc_idx:  continue
      param = reparametrisation.add(
        _.dependent_occupancy,
        occupancy = occupancy,
        multiplier = sc[1]*ratio,
        as_one = as_one,
        scatterer = reparametrisation.structure.scatterers()[sc[0]])
      reparametrisation.asu_scatterer_parameters[sc[0]].occupancy = param
      reparametrisation.shared_occupancies[sc[0]] = occupancy
    self.value = occupancy
