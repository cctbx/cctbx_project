""" All X-H bond lengths are in Angstrom and their values are taken from
ShelXL documentation (p. 4-3) """

import iotbx.constraints.geometrical as _input
import smtbx.refinement.constraints as _
from scitbx.matrix import col

class geometrical_hydrogens_mixin(object):

  def add_to(self, reparametrisation):
    i_pivot = self.pivot
    scatterers = reparametrisation.structure.scatterers()
    pivot_site = scatterers[i_pivot].site
    pivot_site_param = reparametrisation.add_new_site_parameter(i_pivot)
    pivot_neighbour_sites = ()
    pivot_neighbour_site_params = ()
    for j, ops in reparametrisation.pair_sym_table[i_pivot].items():
      if j in self.constrained_site_indices: continue
      for op in ops:
        s = reparametrisation.add_new_site_parameter(j, op)
        pivot_neighbour_site_params += (s,)
        pivot_neighbour_sites += (op*scatterers[j].site,)

    bond_length = reparametrisation.add(
      _.independent_scalar_parameter,
      value=self.ideal_bond_length(scatterers[i_pivot].scattering_type,
                                   reparametrisation.temperature),
      variable=self.stretching)

    hydrogens = tuple(
      [ scatterers[i_sc] for i_sc in self.constrained_site_indices ])

    param = self.add_hydrogen_to(reparametrisation, bond_length,
                                 pivot_site      , pivot_neighbour_sites,
                                 pivot_site_param, pivot_neighbour_site_params,
                                 hydrogens)
    for i_sc in self.constrained_site_indices:
      reparametrisation.asu_scatterer_parameters[i_sc].site = param

  def ideal_bond_length(self, pivot_element, temperature):
    d = self.room_temperature_bond_length[pivot_element]
    if temperature is not None:
      if   temperature < -20: d += 0.1
      elif temperature < -70: d += 0.2
    return d


class terminal_tetrahedral_xhn_site_mixin(geometrical_hydrogens_mixin):

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens):
    assert len(pivot_neighbour_site_params) == 1
    azimuth = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=True)
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


class terminal_tetrahedral_xh_site(_input.terminal_tetrahedral_xh_site,
                                   terminal_tetrahedral_xhn_site_mixin):

  room_temperature_bond_length = { 'O' : 0.82,
                                   }


class tertiary_ch_site(_input.tertiary_ch_site,
                       geometrical_hydrogens_mixin):

  room_temperature_bond_length = { 'C' : 0.98,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens):
    assert len(pivot_neighbour_site_params) == 3
    return reparametrisation.add(
      _.tertiary_ch_site,
      pivot=pivot_site_param,
      pivot_neighbour_0=pivot_neighbour_site_params[0],
      pivot_neighbour_1=pivot_neighbour_site_params[1],
      pivot_neighbour_2=pivot_neighbour_site_params[2],
      length=bond_length,
      hydrogen=hydrogens[0])


class secondary_ch2_sites(_input.secondary_ch2_sites,
                          geometrical_hydrogens_mixin):

  room_temperature_bond_length = { 'C' : 0.97,
                                   }

  def add_hydrogen_to(self, reparametrisation, bond_length,
                      pivot_site      , pivot_neighbour_sites,
                      pivot_site_param, pivot_neighbour_site_params,
                      hydrogens):
    assert len(pivot_neighbour_site_params) == 2
    flapping = reparametrisation.add(_.angle_starting_tetrahedral,
                                     variable=True)
    return reparametrisation.add(
      _.secondary_ch2_sites,
      pivot=pivot_site_param,
      pivot_neighbour_0=pivot_neighbour_site_params[0],
      pivot_neighbour_1=pivot_neighbour_site_params[1],
      length=bond_length,
      h_c_h_angle=flapping,
      hydrogen_0=hydrogens[0],
      hydrogen_1=hydrogens[1])
