from __future__ import division
from cctbx import crystal
from cctbx import xray
from cctbx import sgtbx
from cctbx import adp_restraints, geometry_restraints
import scitbx.math

import libtbx.load_env
if (libtbx.env.has_module(name="smtbx")):
  from smtbx.refinement.restraints import adp_restraints as smtbx_adp_restraints
else:
  smtbx_adp_restraints = None

import iotbx.constraints.commonplace
import iotbx.constraints.geometrical

class crystal_symmetry_builder(object):

  def make_crystal_symmetry(self, unit_cell, space_group):
    self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                             space_group=space_group)


class crystal_structure_builder(crystal_symmetry_builder):

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

  def add_scatterer(self, scatterer, behaviour_of_variable,
                    occupancy_includes_symmetry_factor):
    """ If the parameter set_grad_flags passed to the constructor was True,
        the scatterer.flags.grad_xxx() will be set to True
        if the corresponding variables have been found to be refined
        by the parser using this builder.
    """
    _ = iotbx.constraints.commonplace
    if self.set_grad_flags:
      f = scatterer.flags
      if behaviour_of_variable[0:3].count(_.constant_parameter) != 3:
        f.set_grad_site(True)
      if behaviour_of_variable[3] != _.constant_parameter:
        f.set_grad_occupancy(True)
      if f.use_u_iso():
        if behaviour_of_variable[4] != _.constant_parameter:
          f.set_grad_u_iso(True)
      else:
        if behaviour_of_variable[-6:].count(_.constant_parameter) != 3:
          f.set_grad_u_aniso(True)
    self.structure.add_scatterer(scatterer)

    if occupancy_includes_symmetry_factor:
      sc = self.structure.scatterers()[-1]
      sc.occupancy /= sc.weight_without_occupancy()
      occ = scitbx.math.continued_fraction.from_real(sc.occupancy, eps=1e-5)
      r_occ = occ.as_rational()
      sc.occupancy = round(r_occ.numerator() / r_occ.denominator(), 5)


class constrained_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, constraint_factory=iotbx.constraints.geometrical,
               *args, **kwds):
    super(constrained_crystal_structure_builder, self).__init__(*args, **kwds)
    self.constraint_factory = constraint_factory
    self.geometrical_constraints = []

  def start_geometrical_constraint(self, type_,
                                   bond_length, rotating, stretching,
                                   pivot_relative_pos):
    self.first = len(self.structure.scatterers())

    self.current = type_(rotating=rotating,
                         stretching=stretching,
                         bond_length=bond_length,
                         pivot=self.first + pivot_relative_pos)

  def end_geometrical_constraint(self):
    last = len(self.structure.scatterers())
    self.current.finalise(self.first, last)
    self.geometrical_constraints.append(self.current)

  def finish(self):
    pass

class restrained_crystal_structure_builder(crystal_structure_builder):

  geometry_restraint_types = {
    'bond': geometry_restraints.bond_simple_proxy,
    'angle': geometry_restraints.angle_proxy,
    'dihedral': geometry_restraints.dihedral_proxy,
    'planarity': geometry_restraints.planarity_proxy,
    'chirality': geometry_restraints.chirality_proxy,
    'bond_similarity': geometry_restraints.bond_similarity_proxy,
  }
  if smtbx_adp_restraints is not None:
    adp_restraint_types = {
      'adp_similarity': smtbx_adp_restraints.adp_similarity_restraints,
      'rigid_bond': smtbx_adp_restraints.rigid_bond_restraints,
      'isotropic_adp': smtbx_adp_restraints.isotropic_adp_restraints
    }
  else:
    adp_restraint_types = {}

  def __init__(self, *args, **kwds):
    super(restrained_crystal_structure_builder, self).__init__(*args, **kwds)
    geom = geometry_restraints
    adp = adp_restraints
    self._proxies = {}

    self._proxies = {
      'bond': geometry_restraints.shared_bond_simple_proxy(),
      'angle': geometry_restraints.shared_angle_proxy(),
      'dihedral': geometry_restraints.shared_dihedral_proxy(),
      'planarity': geometry_restraints.shared_planarity_proxy(),
      'chirality': geometry_restraints.shared_chirality_proxy(),
      'bond_similarity': geometry_restraints.shared_bond_similarity_proxy(),
      'adp_similarity': adp_restraints.shared_adp_similarity_proxy(),
      'rigid_bond': adp_restraints.shared_rigid_bond_proxy(),
      'isotropic_adp': adp_restraints.shared_isotropic_adp_proxy(),
    }

  def add_proxy(self, restraint_type, *args, **kwds):
    if restraint_type in self.adp_restraint_types:
      kwds['proxies'] = self._proxies[restraint_type]
      kwds['xray_structure'] = self.structure
      self.adp_restraint_types[restraint_type](**kwds)
    else:
      proxy_type = self.geometry_restraint_types[restraint_type]
      proxy=proxy_type(**kwds)
      self._proxies[restraint_type].append(proxy)

  def process_restraint(self, restraint_type, *args, **kwds):
    def replace_None_with_unit_matrix(sym_ops):
      for i, sym_op in enumerate(sym_ops):
        if sym_op is None:
          sym_ops[i] = sgtbx.rt_mx()
      return sym_ops
    if 'sym_ops' in kwds:
      sym_ops = kwds['sym_ops']
      if restraint_type == 'bond_similarity':
        for i, sym_op in enumerate(sym_ops):
          if sym_op is not None and not isinstance(sym_op, sgtbx.rt_mx):
            if len(sym_op) == 2:
              sym_op_pair = sym_op
              for j, sym_op in enumerate(sym_op_pair):
                if sym_op is not None and not isinstance(sym_op, sgtbx.rt_mx):
                  sym_op_pair[j] = sgtbx.rt_mx(sym_op)
              sym_ops[i] = sym_op_pair
            else:
              sym_ops[i] = sgtbx.rt_mx(sym_op)
      else:
        for i, sym_op in enumerate(sym_ops):
          if sym_op is not None and not isinstance(sym_op, sgtbx.rt_mx):
            sym_ops[i] = sgtbx.rt_mx(sym_op)
      if sym_ops.count(None) == len(sym_ops):
        del kwds['sym_ops']
      elif restraint_type == 'bond':
        # map to asu if necessary
        if sym_ops.count(None) == 2:
          rt_mx_ji = None
        elif sym_ops.count(None) == 1:
          rt_mx_ji = sym_ops[1]
          if rt_mx_ji is None:
            rt_mx_ji = sym_ops[0]
            kwds['i_seqs'].reverse()
        else:
          rt_mx_ji_1 = sym_ops[0]
          rt_mx_ji_2 = sym_ops[1]
          rt_mx_ji_inv = rt_mx_ji_1.inverse()
          rt_mx_ji = rt_mx_ji_inv.multiply(rt_mx_ji_2)
        kwds['rt_mx_ji'] = rt_mx_ji
        del kwds['sym_ops']
      elif restraint_type == 'bond_similarity':
        kwds['sym_ops'] = []
        for i, sym_ops_ in enumerate(sym_ops):
          if sym_ops_ is None or isinstance(sym_ops_, sgtbx.rt_mx): continue
          if sym_ops_.count(None) == 2:
            rt_mx_ji = None
          elif sym_ops_.count(None) == 1:
            rt_mx_ji = sym_ops_[1]
            if rt_mx_ji is None:
              rt_mx_ji = sym_ops_[0]
              kwds['i_seqs'][i].reverse()
          else:
            rt_mx_ji_1 = sym_ops_[0]
            rt_mx_ji_2 = sym_ops_[1]
            rt_mx_ji_inv = rt_mx_ji_1.inverse()
            rt_mx_ji = rt_mx_ji_inv.multiply(rt_mx_ji_2)
          kwds['sym_ops'].append(rt_mx_ji)
        if kwds['sym_ops'].count(None) == len(sym_ops):
          del kwds['sym_ops']
        else:
          kwds['sym_ops'] = replace_None_with_unit_matrix(kwds['sym_ops'])
      else:
        kwds['sym_ops'] = replace_None_with_unit_matrix(sym_ops)
    elif 'rt_mx_ji' in kwds and kwds['rt_mx_ji'] is None:
      del kwds['rt_mx_ji']
    self.add_proxy(restraint_type, **kwds)

  def finish(self):
    pass

  def proxies(self):
    return dict([(proxy_type, proxies) for proxy_type, proxies in self._proxies.iteritems()
                 if len(proxies) != 0])
