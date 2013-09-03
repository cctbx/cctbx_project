from __future__ import division
import smtbx.refinement.constraints as _

class shared_site(object):
  """ a site shared by two or more scatterers
  """

  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  def get_parameter_set(self, reparametrisation):
    rv_l = []
    for s in self.indices[1:]: rv_l.append("%s_xyz" %s)
    rv = set(rv_l)
    if len(rv_l) != len(rv) or len(reparametrisation.constrained_parameters&rv) != 0:
      print("Redundant atoms in %s - '%s' skipping" %(
        self.__class__.__name__,
        reparametrisation.format_scatter_list(self.indices)))
      return None
    return rv

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    site = reparametrisation.add_new_site_parameter(self.indices[0])
    for i in xrange(1, len(self.indices)):
      param = reparametrisation.add(
        _.shared_site,
        reference=site,
        scatterer = reparametrisation.structure.scatterers()[
          self.indices[i]])
      reparametrisation.asu_scatterer_parameters[self.indices[i]].site = param
      reparametrisation.shared_sites[self.indices[i]] = site

