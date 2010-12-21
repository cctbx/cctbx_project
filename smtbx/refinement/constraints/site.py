import smtbx.refinement.constraints as _

class shared_site(object):
  """ a site shared by two or more scatterers
  """

  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  def __eq__(self, other):
    if self.indices != other.indices:  return False
    return True

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    site = reparametrisation.add_new_site_parameter(self.indices[0])
    for i in xrange(1, len(self.indices)):
      param = reparametrisation.add(
        _.shared_site,
        site=site,
        scatterer = reparametrisation.structure.scatterers()[
          self.indices[i]])
      reparametrisation.asu_scatterer_parameters[self.indices[i]].site = param
      reparametrisation.shared_sites[self.indices[i]] = u_c

