""" All X-H bond lengths are in Angstrom and their values are taken from
ShelXL documentation (p. 4-3) """

import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint, geometrical
from scitbx.matrix import col

class hydrogens(geometrical.any):

  need_pivot_neighbour_substituent = False

  def add_to(self, reparametrisation):
    i_pivot = self.pivot
    scatterers = reparametrisation.structure.scatterers()
    conformer_indices = reparametrisation.connectivity_table.conformer_indices
    if conformer_indices is not None:
      constrained_site_conformer = conformer_indices[
        self.constrained_site_indices[0]]
      for i in self.constrained_site_indices:
        assert conformer_indices[i] == constrained_site_conformer
    else: constrained_site_conformer = 0
    pivot_site = scatterers[i_pivot].site
    pivot_site_param = reparametrisation.add_new_site_parameter(i_pivot)
    pivot_neighbour_sites = ()
    pivot_neighbour_site_params = ()
    pivot_neighbour_substituent_site_param = None
    for j, ops in reparametrisation.pair_sym_table[i_pivot].items():
      if j in self.constrained_site_indices: continue
      for op in ops:
        if (conformer_indices is None or
            conformer_indices[j] == 0 or
            constrained_site_conformer == 0 or
            (conformer_indices[j] == constrained_site_conformer)):
          s = reparametrisation.add_new_site_parameter(j, op)
          pivot_neighbour_site_params += (s,)
          pivot_neighbour_sites += (op*scatterers[j].site,)
          if (self.need_pivot_neighbour_substituent
              and pivot_neighbour_substituent_site_param is None):
            for k, ops_k in reparametrisation.pair_sym_table[j].items():
              if k != i_pivot:
                pivot_neighbour_substituent_site_param = \
                  reparametrisation.add_new_site_parameter(k, ops_k[0])
                break

    bond_length = reparametrisation.add(
      _.independent_scalar_parameter,
      value=self.ideal_bond_length(scatterers[i_pivot],
                                   reparametrisation.temperature),
      variable=self.stretching)

    hydrogens = tuple(
      [ scatterers[i_sc] for i_sc in self.constrained_site_indices ])

    param = self.add_hydrogen_to(
      reparametrisation=reparametrisation,
      bond_length=bond_length,
      pivot_site=pivot_site,
      pivot_neighbour_sites=pivot_neighbour_sites,
      pivot_site_param=pivot_site_param,
      pivot_neighbour_site_params=pivot_neighbour_site_params,
      pivot_neighbour_substituent_site_param=
        pivot_neighbour_substituent_site_param,
      hydrogens=hydrogens)
    for i_sc in self.constrained_site_indices:
      reparametrisation.asu_scatterer_parameters[i_sc].site = param

  def ideal_bond_length(self, pivot, temperature):
    pivot_element = pivot.scattering_type
    d = self.room_temperature_bond_length.get(pivot_element)
    if d is None:
      raise InvalidConstraint(
        "Invalid %s constraint involving %s:"
        " ideal bond length not defined to atom type %s" %(
          self.__class__.__name__, pivot.label, pivot_element))
    if temperature is not None:
      if   temperature < -70: d += 0.02
      elif temperature < -20: d += 0.01
    return d


class terminal_tetrahedral_xhn_site(hydrogens):

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    if len(pivot_neighbour_site_params) != 1:
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    azimuth = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotating)
    uc = reparametrisation.structure.unit_cell()
    return reparametrisation.add(
      getattr(_, self.__class__.__name__),
      pivot=pivot_site_param,
      pivot_neighbour=pivot_neighbour_site_params[0],
      length=bond_length,
      azimuth=azimuth,
      e_zero_azimuth=uc.orthogonalize(
        col(hydrogens[0].site) - col(pivot_site)),
      hydrogen=hydrogens)


class terminal_tetrahedral_xh_site(terminal_tetrahedral_xhn_site):
  n_constrained_sites = 1
  room_temperature_bond_length = { 'O' : 0.82,
                                   }

class terminal_tetrahedral_xh3_sites(terminal_tetrahedral_xhn_site):
  n_constrained_sites = 3
  room_temperature_bond_length = { 'C' : 0.96,
                                   'N' : 0.89,
                                   }


class tertiary_xh_site(hydrogens):

  n_constrained_sites = 1
  room_temperature_bond_length = { 'C' : 0.98,
                                   'N' : 0.91,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    if len(pivot_neighbour_site_params) != 3:
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    return reparametrisation.add(
      _.tertiary_xh_site,
      pivot=pivot_site_param,
      pivot_neighbour_0=pivot_neighbour_site_params[0],
      pivot_neighbour_1=pivot_neighbour_site_params[1],
      pivot_neighbour_2=pivot_neighbour_site_params[2],
      length=bond_length,
      hydrogen=hydrogens[0])


class secondary_xh2_sites(hydrogens):

  n_constrained_sites = 2
  room_temperature_bond_length = { 'C' : 0.97,
                                   'N' : 0.90,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    if len(pivot_neighbour_site_params) != 2:
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    x_h = [ col(h.site) for h in hydrogens ]
    x_p = col(pivot_site)
    uc = reparametrisation.structure.unit_cell()
    theta = col(uc.orthogonalize(x_h[0] - x_p)).angle(
      col(uc.orthogonalize(x_h[1] - x_p)))
    flapping = reparametrisation.add(_.independent_scalar_parameter,
                                     value=theta,
                                     variable=True)
    return reparametrisation.add(
      _.secondary_xh2_sites,
      pivot=pivot_site_param,
      pivot_neighbour_0=pivot_neighbour_site_params[0],
      pivot_neighbour_1=pivot_neighbour_site_params[1],
      length=bond_length,
      h_c_h_angle=flapping,
      hydrogen_0=hydrogens[0],
      hydrogen_1=hydrogens[1])


class secondary_planar_xh_site(hydrogens):

  n_constrained_sites = 1
  room_temperature_bond_length = { 'C' : 0.93,
                                   'N' : 0.86,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    # e.g. Carbon atoms in Cyclopentadienyl complexes will have
    #      3 pivot neighbours
    if len(pivot_neighbour_site_params) not in (2, 3):
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    return reparametrisation.add(
      _.secondary_planar_xh_site,
      pivot=pivot_site_param,
      pivot_neighbour_0=pivot_neighbour_site_params[0],
      pivot_neighbour_1=pivot_neighbour_site_params[1],
      length=bond_length,
      hydrogen=hydrogens[0])


class terminal_planar_xh2_sites(hydrogens):

  n_constrained_sites = 2
  need_pivot_neighbour_substituent = True

  room_temperature_bond_length = \
    secondary_planar_xh_site.room_temperature_bond_length

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      pivot_neighbour_substituent_site_param,
                      hydrogens, **kwds):
    if len(pivot_neighbour_site_params) != 1:
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    return reparametrisation.add(
      _.terminal_planar_xh2_sites,
      pivot=pivot_site_param,
      pivot_neighbour=pivot_neighbour_site_params[0],
      pivot_neighbour_substituent=pivot_neighbour_substituent_site_param,
      length=bond_length,
      hydrogen_0=hydrogens[0],
      hydrogen_1=hydrogens[1])


class terminal_linear_ch_site(hydrogens):

  n_constrained_sites = 1
  room_temperature_bond_length = { 'C' : 0.93,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    if len(pivot_neighbour_site_params) != 1:
      raise InvalidConstraint(_.bad_connectivity_msg %(
        self.__class__.__name__, pivot_site_param.scatterers[0].label))
    return reparametrisation.add(
      _.terminal_linear_ch_site,
      pivot=pivot_site_param,
      pivot_neighbour=pivot_neighbour_site_params[0],
      length=bond_length,
      hydrogen=hydrogens[0])

class staggered_terminal_tetrahedral_xh3_sites(hydrogens):

  n_constrained_sites = 3
  staggered = True
  room_temperature_bond_length = \
    terminal_tetrahedral_xh3_sites.room_temperature_bond_length

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    raise NotImplementedError


class staggered_terminal_tetrahedral_xh_site(hydrogens):

  n_constrained_sites = 1
  staggered = True
  room_temperature_bond_length = \
    terminal_tetrahedral_xh_site.room_temperature_bond_length

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    raise NotImplementedError


class polyhedral_bh_site(hydrogens):

  n_constrained_sites = 5
  room_temperature_bond_length = { 'B': 1.10, }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens, **kwds):
    raise NotImplementedError
