import boost.python
ext = boost.python.import_ext("smtbx_refinement_constraints_ext")
from smtbx_refinement_constraints_ext import *

import scitbx.sparse # import dependency
from scitbx.array_family import flex # import dependency
import libtbx.utils

class InvalidConstraint(libtbx.utils.Sorry):
  __module__ = Exception.__module__

bad_connectivity_msg = "Invalid %s constraint involving %s: bad connectivity"

class _(boost.python.injector, ext.parameter):

  def arguments(self):
    """ An iterator over its arguments """
    for i in xrange(self.n_arguments):
      yield self.argument(i)

  def __str__(self):
    """ String representation using the graphviz DOT language """
    try:
      scatt = ', '.join([ sc.label for sc in self.scatterers ])
      scatt = " (%s)" % scatt
    except AttributeError:
      scatt = ""
    info = []
    if not self.is_variable: info.append('cst')
    if info: info = ' [%s]' % ', '.join(info)
    else: info = ''
    lbl = '%i [label="%s%s%s #%s"]' % (
      self.index,
      self.__class__.__name__,
      info,
      scatt, self.index)
    return lbl

class _(boost.python.injector, ext.reparametrisation):

  def __str__(self):
    """ String representation using the graphviz DOT language """
    self.finalise()
    bits = []
    for p in self.parameters():
      for q in p.arguments():
        bits.append("%i -> %i" % (p.index, q.index))
    dsu_bits = []
    for p in self.parameters():
      dsu_bits.append((p.index, str(p)))
    dsu_bits.sort()
    bits.extend([ p for i,p in dsu_bits ])
    return "digraph dependencies {\n%s\n}" % ';\n'.join(bits)


# The order in which constraints are added MAKES a difference, shared site, U and/or
# occupancy constraints must be added first for proper bookkeeping

class reparametrisation(ext.reparametrisation):
  """ Enhance the C++ level reparametrisation class for ease of use """

  temperature = 20 # Celsius
  twin_fractions = None

  def __init__(self,
               structure,
               constraints,
               connectivity_table,
               **kwds):
    """ Construct for the given instance of xray.structure subject to the
    given sequence of constraints. Each constraint instance shall understand:
    constraint.add_to(self). That method shall perform 2 tasks:

      - add to self the parameters relevant to the reparametrisation
        associated with that constraint;

      - update self.asu_scatterer_parameters.

    The latter is an array containing one instance of scatterer_parameters
    for each scatterer in the a.s.u.
    C.f. module geometrical_hydrogens in this package for a typical example
    """
    super(reparametrisation, self).__init__(structure.unit_cell())
    #association of scatterer_idx:parameter, bookkeeping
    self.shared_Us = {}
    self.shared_sites = {}
    self.shared_occupancies = {}

    self.structure = xs = structure
    self.connectivity_table = connectivity_table
    self.pair_sym_table = \
        connectivity_table.pair_asu_table.extract_pair_sym_table(
          skip_j_seq_less_than_i_seq=False,
          all_interactions_from_inside_asu=True)
    scatterers = xs.scatterers()
    self.site_symmetry_table_ = self.structure.site_symmetry_table()
    libtbx.adopt_optional_init_args(self, kwds)
    self.asu_scatterer_parameters = shared_scatterer_parameters(xs.scatterers())
    self.independent_scalar_parameters = shared_independent_shared_parameters()

    for constraint in constraints:
      constraint.add_to(self)

    for i_sc in xrange(len(self.asu_scatterer_parameters)):
      self.add_new_site_parameter(i_sc)
      self.add_new_thermal_displacement_parameter(i_sc)
      self.add_new_occupancy_parameter(i_sc)
    if self.twin_fractions is not None:
      for twin_fraction in self.twin_fractions:
        self.add_new_independent_scalar_parameter(twin_fraction)
    self.finalise()

  def finalise(self):
    super(reparametrisation, self).finalise()
    self.mapping_to_grad_fc = \
        self.asu_scatterer_parameters.mapping_to_grad_fc()
    self.mapping_to_grad_fc_independent_scalars = \
        self.independent_scalar_parameters.mapping_to_grad_fc()
    self.mapping_to_grad_fc_all = self.mapping_to_grad_fc.deep_copy()
    self.mapping_to_grad_fc_all.extend(self.mapping_to_grad_fc_independent_scalars)

  def apply_shifts(self, shifts):
    ext.reparametrisation.apply_shifts(self, shifts)
    if self.twin_fractions is not None:
      self.twin_fractions += shifts[-self.twin_fractions.size():]

  class component_annotations(libtbx.property):
    def fget(self):
      return self.__dict__.setdefault(
        "_component_annotations",
        self.asu_scatterer_parameters.component_annotations().split(',')[:-1])

  class n_independent_params(libtbx.property):
    def fget(self):
      return self.jacobian_transpose.n_rows

  def jacobian_transpose_matching_grad_fc(self):
    """ The columns of self.jacobian_transpose corresponding to crystallographic
    parameters for the scatterers, in the same order as the derivatives in
    grad Fc. In this class, the latter is assumed to follow the convention of
    smtbx.structure_factors
    """
    return self.jacobian_transpose_matching(self.mapping_to_grad_fc)

  def jacobian_transpose_matching(self, mapping):
    """ The columns of self.jacobian_transpose corresponding to the given
    crystallographic parameters.
    """
    return self.jacobian_transpose.select_columns(mapping)

  def add_new_occupancy_parameter(self, i_sc):
    if self.shared_occupancies.has_key(i_sc):
      return self.shared_occupancies[i_sc]
    occ = self.asu_scatterer_parameters[i_sc].occupancy
    if occ is None:
      sc = self.structure.scatterers()[i_sc]
      occ = self.add(independent_occupancy_parameter, sc)
      self.asu_scatterer_parameters[i_sc].occupancy = occ
    return occ

  def add_new_site_parameter(self, i_scatterer, symm_op=None):
    if self.shared_sites.has_key(i_scatterer):
      return self.shared_sites[i_scatterer]
    s = self.asu_scatterer_parameters[i_scatterer].site
    if s is None:
      site_symm = self.site_symmetry_table_.get(i_scatterer)
      sc = self.structure.scatterers()[i_scatterer]
      if site_symm.is_point_group_1():
        s = self.add(independent_site_parameter, sc)
      else:
        s = self.add(special_position_site_parameter, site_symm, sc)
      self.asu_scatterer_parameters[i_scatterer].site = s
    if symm_op is not None and not symm_op.is_unit_mx():
      s = self.add(symmetry_equivalent_site_parameter, s, symm_op)
    return s

  def add_new_site_proxy_parameter(self, param, i, i_sc):
    if self.shared_sites.has_key(i_sc):
      return self.shared_sites[i_sc]
    sc = self.structure.scatterers()[i_sc]
    self.shared_sites[i_sc] = self.add(rigid_site_proxy, param, i)

  def add_new_thermal_displacement_parameter(self, i_scatterer):
    if self.shared_Us.has_key(i_scatterer):
      return self.shared_Us[i_scatterer]
    u = self.asu_scatterer_parameters[i_scatterer].u
    if u is None:
      sc = self.structure.scatterers()[i_scatterer]
      assert not (sc.flags.use_u_iso() and sc.flags.use_u_aniso())
      if sc.flags.use_u_iso():
        u = self.add(independent_u_iso_parameter, sc)
      else:
        site_symm = self.site_symmetry_table_.get(i_scatterer)
        if site_symm.is_point_group_1():
          u = self.add(independent_u_star_parameter, sc)
        else:
          u = self.add(special_position_u_star_parameter,
                       site_symm,
                       sc)
      self.asu_scatterer_parameters[i_scatterer].u = u
    return u

  def add_new_independent_scalar_parameter(self, value, variable=True):
    p = self.add(independent_scalar_parameter, value=value, variable=variable)
    self.independent_scalar_parameters.append(p)
    return p
