import boost.python
boost.python.import_ext("smtbx_refinement_constraints_ext")
import smtbx_refinement_constraints_ext as ext
from smtbx_refinement_constraints_ext import *

import libtbx
import scitbx.sparse
from cctbx import crystal
from cctbx.eltbx import covalent_radii

class _parameter(boost.python.injector, ext.parameter):

  def arguments(self):
    for i in xrange(self.n_arguments):
      yield self.argument(i)


class reparametrisation(ext.reparametrisation):

  temperature = 20 # Celsius
  covalent_bond_tolerance = 0.5 # Angstrom

  def __init__(self,
               structure,
               geometrical_constraints,
               **kwds):
    super(reparametrisation, self).__init__(structure.unit_cell())
    self.structure = xs = structure
    scatterers = xs.scatterers()
    self.site_symmetry_table_ = self.structure.site_symmetry_table()
    libtbx.adopt_optional_init_args(self, kwds)
    self.asu_scatterer_parameters = shared_scatterer_parameters(
      len(xs.scatterers()))

    radii = [
      covalent_radii.table(elt).radius()
      for elt in xs.scattering_type_registry().type_index_pairs_as_dict() ]
    self.buffer_thickness = max(radii) + self.covalent_bond_tolerance

    asu_mappings = xs.asu_mappings(buffer_thickness=self.buffer_thickness)
    bond_table = crystal.pair_asu_table(asu_mappings)
    bond_table.add_covalent_pairs(xs.scattering_types(),
                                  tolerance=self.covalent_bond_tolerance)
    self.pair_sym_table = bond_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq=False)

    for constraint in geometrical_constraints:
      constraint.add_to(self)

    for i_sc in xrange(len(self.asu_scatterer_parameters)):
      if self.asu_scatterer_parameters[i_sc].site is None:
        self.add_new_site_parameter(i_sc)
      if self.asu_scatterer_parameters[i_sc].u is None:
        self.add_new_thermal_displacement_parameter(i_sc)
      if self.asu_scatterer_parameters[i_sc].occupancy is None:
        sc = scatterers[i_sc]
        occ = self.add(independent_scalar_parameter,
                       value=sc.occupancy,
                       variable=sc.flags.grad_occupancy())
        self.asu_scatterer_parameters[i_sc].occupancy = occ
    self.finalise()

  def add_new_site_parameter(self, i_scatterer, symm_op=None):
    s = self.asu_scatterer_parameters[i_scatterer].site
    if s is None:
      site_symm = self.site_symmetry_table_.get(i_scatterer)
      sc = self.structure.scatterers()[i_scatterer]
      if site_symm.is_point_group_1():
        s = self.add(independent_site_parameter, sc)
      else:
        s = self.add(special_position_site, site_symm, sc)
      if symm_op is not None and not symm_op.is_unit_mx():
        s = self.add(symmetry_equivalent_site_parameter, s)
      self.asu_scatterer_parameters[i_scatterer].site = s
    return s

  def add_new_thermal_displacement_parameter(self, i_scatterer):
    u = self.asu_scatterer_parameters[i_scatterer].u
    if u is None:
      sc = self.structure.scatterers()[i_scatterer]
      assert sc.flags.use_u_iso() ^ sc.flags.use_u_aniso()
      if sc.flags.use_u_iso():
        u = self.add(independent_scalar_parameter,
                     value=sc.u_iso,
                     variable=sc.flags.grad_u_iso())
      else:
        site_symm = self.site_symmetry_table_.get(i_scatterer)
        if site_symm.is_point_group_1():
          u = self.add(independent_cartesian_adp, sc)
        else:
          u = self.add(special_position_cartesian_adp,
                       site_symm,
                       self.structure.unit_cell(),
                       sc)
      self.asu_scatterer_parameters[i_scatterer].u = u
    return u
