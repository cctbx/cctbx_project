from __future__ import division

import scitbx.math
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import xray
from cctbx import sgtbx
from cctbx import adp_restraints, geometry_restraints
import iotbx.constrained_parameters

import libtbx.load_env
if (libtbx.env.has_module(name="smtbx")):
  from smtbx.refinement.restraints import adp_restraints as smtbx_adp_restraints
else:
  smtbx_adp_restraints = None


def mixin_builder_class(mixin_name, *mixed_builders):
  """ This function will make on-the-fly a builder class with the given
  name, that inherits from the given builders """
  return type.__new__(type, mixin_name, mixed_builders, {})


class crystal_symmetry_builder(object):

  def make_crystal_symmetry(self, unit_cell, space_group):
    self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                             space_group=space_group)


class electron_density_peak(object):

  __slots__ = ('site', 'height')

  def __init__(self, site, height):
    self.site = site
    self.height = height


class crystal_structure_builder(crystal_symmetry_builder):

  def __init__(self,
               set_grad_flags=True,
               min_distance_sym_equiv=0.5):
    super(crystal_structure_builder, self).__init__()
    self.set_grad_flags = set_grad_flags
    self.min_distance_sym_equiv = min_distance_sym_equiv
    self.conformer_indices = None
    self.sym_excl_indices = None
    self.electron_density_peaks = []

  def make_structure(self):
    self.structure = xray.structure(
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=self.crystal_symmetry,
        min_distance_sym_equiv=self.min_distance_sym_equiv))

  def add_scatterer(self, scatterer, behaviour_of_variable,
                    occupancy_includes_symmetry_factor,
                    conformer_index=None,
                    sym_excl_index=None):
    """ If the parameter set_grad_flags passed to the constructor was True,
        the scatterer.flags.grad_xxx() will be set to True
        if the corresponding variables have been found to be refined
        by the parser using this builder.
    """
    _ = iotbx.constrained_parameters
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
        if behaviour_of_variable[-6:].count(_.constant_parameter) != 6:
          f.set_grad_u_aniso(True)
    self.structure.add_scatterer(scatterer)

    if occupancy_includes_symmetry_factor:
      sc = self.structure.scatterers()[-1]
      sc.occupancy /= sc.weight_without_occupancy()
      occ = scitbx.math.continued_fraction.from_real(sc.occupancy, eps=1e-5)
      r_occ = occ.as_rational()
      sc.occupancy = round(r_occ.numerator() / r_occ.denominator(), 5)

    if conformer_index is not None:
      if self.conformer_indices is None: self.conformer_indices = flex.size_t()
      self.conformer_indices.append(conformer_index)
    if sym_excl_index is not None:
      if self.sym_excl_indices is None: self.sym_excl_indices = flex.size_t()
      self.sym_excl_indices.append(sym_excl_index)

  def add_electron_density_peak(self, site, height):
    self.electron_density_peaks.append(electron_density_peak(site, height))


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
      for i, sym_op in enumerate(sym_ops):
        if sym_op is not None and not isinstance(sym_op, sgtbx.rt_mx):
          if len(sym_op) == 2:
            assert restraint_type == 'bond_similarity'
            sym_op_pair = sym_op
            for j, sym_op in enumerate(sym_op_pair):
              if sym_op is not None and not isinstance(sym_op, sgtbx.rt_mx):
                sym_op_pair[j] = sgtbx.rt_mx(sym_op)
            sym_ops[i] = sym_op_pair
          else:
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

  class restraints_manager(libtbx.property):
    def fget(self):
      from smtbx.refinement.restraints import manager
      kwds = dict([ ("%s_proxies" % name, value)
                    for name, value in self.proxies().iteritems() ])
      return manager(**kwds)

  def proxies(self):
    return dict([
      (proxy_type, proxies) for proxy_type, proxies in self._proxies.iteritems()
      if len(proxies) != 0])

class reflection_data_source_builder(object):
  """ A builder which is passed the information about the reflections
  corresponding to the parsed

  Attributes:

    - reflection_file_format: a string to identify the format, compatible with
              iotbx.reflection_file_reader.any_reflection_file
    - data_change_of_basis_op: the instance of sgtbx.change_of_basis_op that
                          shall be applied to the Miller indices to match
                          the parsed model
    - data_scale: the scale of the data and their standard deviations
  """

  def create_shelx_reflection_data_source(self,
                                          format,
                                          indices_transform=None,
                                          change_of_basis_op=None,
                                          data_scale=1):
    """ format is one of 3, 4, 5, etc.
    data_scale scales the data and their standard deviations
    """
    assert [indices_transform, change_of_basis_op].count(None) == 1
    if change_of_basis_op is None:
      if indices_transform.is_unit_mx():
        change_of_basis_op = sgtbx.change_of_basis_op()
      else:
        r = sgtbx.rt_mx(indices_transform.new_denominator(24).transpose())
        change_of_basis_op = sgtbx.change_of_basis_op(r).inverse()
    self.reflection_file_format = "hklf%i" % format
    self.data_change_of_basis_op = change_of_basis_op
    self.data_scale = data_scale

class twinning_builder(object):
  """ Construct twin components
      They come as a tuple of instances of xray.twin_component.
      If R is the twin law, the i-th component corresponds to the domain
      where Miller indices are to be transformed by R^(i+1). The domain
      with untransformed Miller indices is omitted as its fraction is just
      1 minus the sum of the listed fractions.
      As a result, the 0th component holds the twin law.
  """

  twin_components = ()
  non_merohedral_twin_components_with_transformed_hkl = ()

  def make_merohedral_twinning(self, twin_law, fractions):
    self.twin_components = []
    if isinstance(twin_law, sgtbx.rt_mx): twin_law = twin_law.r()
    t = twin_law
    for f in fractions:
      self.twin_components.append(xray.twin_component(twin_law=t,
                                                      twin_fraction=f,
                                                      grad_twin_fraction=True))
      t = t.multiply(twin_law)

  def make_non_merohedral_twinning_with_transformed_hkl(self, fractions):
    self.non_merohedral_twin_components_with_transformed_hkl = fractions


# **********************************************
# Conditional import of smtbx-dependent builders
# **********************************************
if (libtbx.env.has_module(name="smtbx")):
  from iotbx.builders_depending_on_smtbx import *
