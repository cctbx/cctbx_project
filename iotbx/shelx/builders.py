from cctbx import crystal
from cctbx import xray
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager

from iotbx.shelx import util

class crystal_symmetry_builder(object):

  def make_crystal_symmetry(self, unit_cell, space_group):
    self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                             space_group=space_group)


class crystal_structure_builder(crystal_symmetry_builder,
                                util.behaviour_of_variable):

  def __init__(self,
               set_grad_flags=True,
               min_distance_sym_equiv=0.5):
    super(crystal_structure_builder, self).__init__()
    self.set_grad_flags = set_grad_flags
    self.min_distance_sym_equiv = min_distance_sym_equiv

  def make_structure(self):
    self.structure = xray.structure(
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=self.crystal_symmetry,
        min_distance_sym_equiv=self.min_distance_sym_equiv))

  def add_scatterer(self, scatterer, behaviour_of_variable):
    """ If the parameter set_grad_flags passed to the constructor was True,
        the scatterer.flags.grad_xxx() will be set to True
        if the corresponding variables have been found to be refined
        by the parser using this builder.
    """
    if self.set_grad_flags:
      f = scatterer.flags
      if behaviour_of_variable[0:3].count(self.fixed) != 3:
        f.set_grad_site(True)
      if behaviour_of_variable[3] != self.fixed:
        f.set_grad_occupancy(True)
      if f.use_u_iso():
        if behaviour_of_variable[4] != self.fixed:
          f.set_grad_u_iso(True)
      else:
        if behaviour_of_variable[-6:].count(self.fixed) != 3:
          f.set_grad_u_aniso(True)
    self.structure.add_scatterer(scatterer)


class afixed_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, *args, **kwds):
    super(afixed_crystal_structure_builder, self).__init__(*args, **kwds)
    self.afixed = []
    self.afix = None

  def start_afix(self, constraint_type, kwds):
    self.afix = (constraint_type, kwds,
                 len(self.structure.scatterers()))

  def end_afix(self):
    last = len(self.structure.scatterers())
    constraint_type, kwds, first = self.afix
    kwds['constrained_scatterer_indices'] = tuple(xrange(first, last))
    self.afixed.append((constraint_type, kwds))

  def finish(self):
    pass


class restrained_crystal_structure_builder(afixed_crystal_structure_builder):
  def __init__(self, *args, **kwds):
    super(restrained_crystal_structure_builder, self).__init__(*args, **kwds)
    self.bond_sym_proxies = []
    self.angle_proxies = geometry_restraints.shared_angle_proxy()
    self.planarity_proxies = geometry_restraints.shared_planarity_proxy()
    self.dihedral_proxies = geometry_restraints.shared_dihedral_proxy()

  def add_restraint(self, restraint_type, kwds):
    restraint=restraint_type(**kwds)
    if 'bond_sym' in restraint_type.__name__:
      self.bond_sym_proxies.append(restraint)
    elif 'planarity' in restraint_type.__name__:
      self.planarity_proxies.append(restraint)
    elif 'angle' in restraint_type.__name__:
      self.angle_proxies.append(restraint)
    elif 'dihedral' in restraint_type.__name__:
      self.shared_proxies[restraint_type].append(restraint)

  def finish_restraints(self):
    max_bond_distance = 2
    bond_params_table = geometry_restraints.bond_params_table(
      self.structure.scatterers().size())
    asu_mappings = self.structure.asu_mappings(buffer_thickness=max_bond_distance*3)
    bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    for proxy in self.bond_sym_proxies:
      i_seq, j_seq = proxy.i_seqs
      bond_params_table.update(
        i_seq=i_seq,
        j_seq=j_seq,
        params=proxy)
      bond_asu_table.add_pair(
          i_seq=i_seq,
          j_seq=j_seq,
          rt_mx_ji=proxy.rt_mx_ji)
    self.geometry_restraints_manager = geometry_restraints.manager.manager(
      crystal_symmetry=self.crystal_symmetry,
      site_symmetry_table=self.structure.site_symmetry_table(),
      bond_params_table=bond_params_table,
      shell_sym_tables=[bond_asu_table.extract_pair_sym_table()],
      angle_proxies=self.angle_proxies,
      planarity_proxies= self.planarity_proxies,
      dihedral_proxies=self.dihedral_proxies)
