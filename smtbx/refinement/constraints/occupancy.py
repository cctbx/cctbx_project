from __future__ import absolute_import, division, print_function
import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
from scitbx.array_family import flex
import itertools
from six.moves import range

class occupancy_affine_constraint(object):
  """ Constraint a0 occ0 + a1 occ1 +... == b """

  def __init__(self, scatterer_indices, a, b):
    self.a = []
    self.scatterer_indices = []
    for i in range(len(a)):
      if(a[i]!=0.0):
        self.a += [a[i]]
        self.scatterer_indices += [scatterer_indices[i]]
    self.a = flex.double(self.a)
    self.b = b

  @property
  def constrained_parameters(self):
    return ((self.scatterer_indices[0], 'occupancy'),)

  def add_to(self, reparametrisation):
    sc = reparametrisation.structure.scatterers()
    dependees = [reparametrisation.add_new_occupancy_parameter(i)\
      for i in self.scatterer_indices[1:]]
    sidx = self.scatterer_indices[0]
    param = reparametrisation.add(_.affine_asu_occupancy_parameter,
                                  dependees=dependees, a=-1.0*self.a[1:]/(self.a[0]), b=self.b/(self.a[0]),
                                  scatterer=sc[sidx])
    reparametrisation.asu_scatterer_parameters[sidx].occupancy = param
    for idx,i in enumerate(self.scatterer_indices[1:]):
      reparametrisation.shared_occupancies[i] = dependees[idx]
    reparametrisation.shared_occupancies[sidx] = param
    self.value = param

class occupancy_pair_affine_constraint(object):
  """ Constraint a0 occ0 + a1 occ1 == b """

  def __init__(self, scatterer_indices, linear_form):
    self.scatterer_indices = scatterer_indices
    self.linear_form = linear_form

  @property
  def constrained_parameters(self):
    return ((self.scatterer_indices[1], 'occupancy'),)

  def add_to(self, reparametrisation):
    sc = reparametrisation.structure.scatterers()
    (a0, a1), b = self.linear_form
    i0, i1 = self.scatterer_indices
    occ0 = reparametrisation.add_new_occupancy_parameter(i0)
    reparametrisation.shared_occupancies[i0] = occ0
    param = reparametrisation.add(_.affine_asu_occupancy_parameter,
                                  dependee=occ0, a=-a0/a1, b=b/a1,
                                  scatterer=sc[i1])
    reparametrisation.asu_scatterer_parameters[i1].occupancy = param
    reparametrisation.shared_occupancies[i1] = param


class dependent_occupancy(object):
  """ occupancy of a site depend on the occupancy of the other site
  """

  def __init__(self, var_refs, var_minus_one_refs):
    if (len(var_refs) + len(var_minus_one_refs)) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.as_var = var_refs
    self.as_one_minus_var = var_minus_one_refs

  @property
  def constrained_parameters(self):
    ids = [(sc[0], 'occupancy')
            for sc in itertools.chain(self.as_var, self.as_one_minus_var)]
    return tuple(ids[1:])

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    if len(self.as_var) != 0:
      sc_idx = self.as_var[0][0]
      original_mult = self.as_var[0][1]
    else:
      sc_idx = self.as_one_minus_var[0][0]
      original_mult = self.as_one_minus_var[0][1]
    self.value = reparametrisation.add_new_occupancy_parameter(sc_idx)
    if sc_idx not in reparametrisation.shared_occupancies:
      reparametrisation.shared_occupancies[sc_idx] = self.value
    for sc in self.as_var:
      if sc[0] == sc_idx:  continue
      param = reparametrisation.add(
        _.dependent_occupancy,
        occupancy = self.value,
        original_multiplier = original_mult,
        multiplier = sc[1],
        as_one = True,
        scatterer = scatterers[sc[0]])
      reparametrisation.asu_scatterer_parameters[sc[0]].occupancy = param
      reparametrisation.shared_occupancies[sc[0]] = param
    as_one = len(self.as_var) == 0  # only if both lists are not empty
    for sc in self.as_one_minus_var:
      if sc[0] == sc_idx:  continue
      param = reparametrisation.add(
        _.dependent_occupancy,
        occupancy = self.value,
        original_multiplier = original_mult,
        multiplier = sc[1],
        as_one = as_one,
        scatterer = scatterers[sc[0]])
      reparametrisation.asu_scatterer_parameters[sc[0]].occupancy = param
      reparametrisation.shared_occupancies[sc[0]] = param
