from __future__ import division
from cctbx import crystal
from cctbx import xray
from cctbx import sgtbx
from cctbx import geometry_restraints
from cctbx import adp_restraints
import scitbx.math

from iotbx.shelx import util
import iotbx.constraints

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

    ## scatterer.occupancy has been filled with the chemical occupancy times
    ## the special position occupancy by the parsers, and this  needs to be
    ## corrected for.
    sc = self.structure.scatterers()[-1]
    sc.occupancy /= sc.weight_without_occupancy()
    occ = scitbx.math.continued_fraction.from_real(sc.occupancy, eps=1e-5)
    r_occ = occ.as_rational()
    sc.occupancy = round(r_occ.numerator() / r_occ.denominator(), 5)


class constrained_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, *args, **kwds):
    super(constrained_crystal_structure_builder, self).__init__(*args, **kwds)
    self.geometrical_constraints = []

  def start_geometrical_constraint(self, type,
                                   bond_length, rotating, use_pivot):
    self.first = len(self.structure.scatterers())
    if use_pivot: pivot = self.first - 1
    else: pivot = None
    self.current = type(bond_length=bond_length,
                        rotating=rotating,
                        pivot=pivot)

  def end_geometrical_constraint(self):
    last = len(self.structure.scatterers())
    self.current.constrained_site_indices = tuple(
      xrange(self.first, last))
    self.geometrical_constraints.append(self.current)

  def finish(self):
    pass


class restrained_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, *args, **kwds):
    super(restrained_crystal_structure_builder, self).__init__(*args, **kwds)
    geom = geometry_restraints
    adp = adp_restraints
    self._proxies = {
      geom.bond_simple_proxy: geom.shared_bond_simple_proxy(),
      geom.angle_proxy: geom.shared_angle_proxy(),
      geom.dihedral_proxy: geom.shared_dihedral_proxy(),
      geom.chirality_proxy: geom.shared_chirality_proxy(),
      geom.planarity_proxy: geom.shared_planarity_proxy(),
      geom.bond_similarity_proxy: geom.shared_bond_similarity_proxy(),
      adp.adp_similarity_proxy: adp.shared_adp_similarity_proxy(),
      adp.isotropic_adp_proxy: adp.shared_isotropic_adp_proxy(),
      adp.rigid_bond_proxy: adp.shared_rigid_bond_proxy(),
    }
    import libtbx.load_env
    if (not libtbx.env.has_module(name="smtbx")):
      raise RuntimeError("smtbx module is not available.")
    from smtbx.refinement.restraints import \
      adp_restraints as smtbx_adp_restraints
    self.adp_proxy_builders = {
      adp.adp_similarity_proxy: smtbx_adp_restraints.adp_similarity_restraints,
      adp.rigid_bond_proxy: smtbx_adp_restraints.rigid_bond_restraints,
      adp.isotropic_adp_proxy: smtbx_adp_restraints.isotropic_adp_restraints
    }
    self.shelx_cmd_to_proxy_type_mappings = {
      'DFIX': geometry_restraints.bond_simple_proxy,
      'DANG': geometry_restraints.bond_sym_proxy,
      'FLAT': geometry_restraints.planarity_proxy,
      'CHIV': geometry_restraints.chirality_proxy,
      'SADI': geometry_restraints.bond_similarity_proxy,
      'SIMU': adp_restraints.adp_similarity_proxy,
      'DELU': adp_restraints.rigid_bond_proxy,
      'ISOR': adp_restraints.isotropic_adp_proxy,
    }

  def add_proxy(self, cmd, kwds):
    proxy_type = self.shelx_cmd_to_proxy_type_mappings[cmd]
    if proxy_type in self.adp_proxy_builders:
      kwds['proxies'] = self._proxies[proxy_type]
      kwds['xray_structure'] = self.structure
      self.adp_proxy_builders[proxy_type](**kwds)
    else:
      proxy=proxy_type(**kwds)
      self._proxies[proxy_type].append(proxy)

  def process_restraint(self, cmd, shelx_kwds):
    def replace_None_with_unit_matrix(sym_ops):
      for i, sym_op in enumerate(sym_ops):
        if sym_op is None:
          sym_ops[i] = sgtbx.rt_mx()
      return sym_ops
    kwds = {}
    i_seqs = shelx_kwds['i_seqs']
    kwds['i_seqs'] = i_seqs
    if cmd in ('DFIX','DANG'):
      sym_ops = shelx_kwds['sym ops']
      if sym_ops.count(None) == 2:
        rt_mx_ji = None # unit matrix
      elif sym_ops.count(None) == 1:
        sym_op = sym_ops[1]
        if sym_op is None:
          sym_op = sym_ops[0]
          kwds['i_seqs'].reverse()
        rt_mx_ji = sgtbx.rt_mx(sym_op)
      else:
        rt_mx_ji_1 = sgtbx.rt_mx(sym_ops[0])
        rt_mx_ji_2 = sgtbx.rt_mx(sym_ops[1])
        rt_mx_ji_inv = rt_mx_ji_1.inverse()
        rt_mx_ji = rt_mx_ji_inv.multiply(rt_mx_ji_2)
      if rt_mx_ji is not None:
        kwds['rt_mx_ji'] = rt_mx_ji
      kwds['distance_ideal'] = shelx_kwds['d']
      esd = shelx_kwds['s']
      kwds['weight'] = 1/(esd**2)
      self.add_proxy(cmd, kwds)
    elif cmd == 'FLAT':
      esd = shelx_kwds['s']
      kwds['weights'] = [1/(esd**2)]*len(i_seqs)
      self.add_proxy(cmd, kwds)
    elif cmd == 'SADI':
      sym_ops = []
      for i in range(len(i_seqs)):
        symmetry_operations = shelx_kwds['sym ops'][i]
        if symmetry_operations.count(None) == 2:
          rt_mx_ji = None
        elif symmetry_operations.count(None) == 1:
          sym_op = symmetry_operations[1]
          if sym_op is None:
            sym_op = symmetry_operations[0]
            kwds['i_seqs'][i].reverse()
          rt_mx_ji = sgtbx.rt_mx(sym_op)
        else:
          rt_mx_ji_1 = sgtbx.rt_mx(symmetry_operations[0])
          rt_mx_ji_2 = sgtbx.rt_mx(symmetry_operations[1])
          rt_mx_ji_inv = rt_mx_ji_1.inverse()
          rt_mx_ji = rt_mx_ji_inv.multiply(rt_mx_ji_2)
        sym_ops.append(rt_mx_ji)
      esd = shelx_kwds['s']
      kwds['weights'] = [1/(esd**2)]*len(i_seqs)
      if sym_ops.count(None) != len(sym_ops):
        kwds['sym_ops'] = replace_None_with_unit_matrix(sym_ops)
      self.add_proxy(cmd, kwds)
    elif cmd in ('DELU', 'ISOR', 'SIMU'):
      self.add_proxy(cmd, shelx_kwds)

  def finish(self):
    pass

  def proxies(self):
    return dict([(proxy_type, proxies) for proxy_type, proxies in self._proxies.iteritems()
                 if len(proxies) != 0])
