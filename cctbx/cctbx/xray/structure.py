from cctbx.xray import ext
from cctbx.xray import structure_factors
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import eltbx
from cctbx import matrix
from cctbx.array_family import flex
import types
import sys

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings, scatterers=None,
                     scattering_dict=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    if (scatterers is None):
      self.erase_scatterers()
    else:
      self._scatterers = scatterers.deep_copy()
      self.all_apply_symmetry()
    self._scattering_dict = scattering_dict
    self._scattering_dict_is_out_of_date = 0001

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._special_position_indices = other._special_position_indices
    self._scattering_dict = other._scattering_dict
    self._scattering_dict_is_out_of_date=other._scattering_dict_is_out_of_date

  def erase_scatterers(self):
    self._scatterers = flex.xray_scatterer()
    self._special_position_indices = flex.size_t()
    self._scattering_dict_is_out_of_date = 0001

  def deep_copy_scatterers(self):
    cp = structure(self, scattering_dict=self._scattering_dict)
    cp._scatterers = self._scatterers.deep_copy()
    cp._special_position_indices = self._special_position_indices.deep_copy()
    return cp

  def scatterers(self):
    return self._scatterers

  def special_position_indices(self):
    return self._special_position_indices

  def scattering_dict(self, custom_dict=None, d_min=None, d_min_it1992=1):
    if (   self._scattering_dict_is_out_of_date
        or custom_dict is not None
        or d_min is not None):
      new_dict = {"const": eltbx.caasf.custom(1) }
      if (self._scattering_dict is not None and d_min is None):
        for k,v in self._scattering_dict.dict().items():
          new_dict[k] = v.coefficients
      if (custom_dict is not None):
        new_dict.update(custom_dict)
      if (d_min is None): d_min = 0
      self._scattering_dict = ext.scattering_dictionary(self.scatterers())
      for key_undef in self._scattering_dict.find_undefined():
        if (new_dict.has_key(key_undef)):
          val = new_dict[key_undef]
        else:
          if (d_min >= d_min_it1992):
            val = eltbx.caasf.it1992(key_undef, 1).as_custom()
          else:
            val = eltbx.caasf.wk1995(key_undef, 1).as_custom()
        self._scattering_dict.assign(key_undef, val)
      self._scattering_dict_is_out_of_date = 00000
    return self._scattering_dict

  def __getitem__(self, slice_object):
    assert type(slice_object) == types.SliceType
    assert self.scatterers() is not None
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().__getitem__(slice_object),
      scattering_dict=self._scattering_dict)

  def add_scatterer(self, scatterer):
    i = self.scatterers().size()
    self._scatterers.append(scatterer)
    self._scattering_dict_is_out_of_date = 0001
    self.apply_symmetry(i)

  def apply_symmetry(self, i, update_special_position_indices=0001):
    site_symmetry = self._scatterers[i].apply_symmetry(
      self.unit_cell(),
      self.space_group(),
      self.min_distance_sym_equiv(),
      self.u_star_tolerance(),
      self.assert_is_positive_definite(),
      self.assert_min_distance_sym_equiv())
    if (update_special_position_indices
        and not site_symmetry.is_point_group_1()):
      self.special_position_indices().append(i)

  def add_scatterers(self, scatterers):
    n = self.scatterers().size()
    special_position_indices = self._all_apply_symmetry(scatterers) + n
    self._scatterers.append(scatterers)
    self._scattering_dict_is_out_of_date = 0001
    self._special_position_indices.append(special_position_indices)

  def all_apply_symmetry(self):
    self._special_position_indices =self._all_apply_symmetry(self.scatterers())

  def _all_apply_symmetry(self, scatterers):
    return ext.apply_symmetry(
      self.unit_cell(),
      self.space_group(),
      scatterers,
      self.min_distance_sym_equiv(),
      self.u_star_tolerance(),
      self.assert_is_positive_definite(),
      self.assert_min_distance_sym_equiv())

  def replace_scatterers(self, scatterers):
    self.erase_scatterers()
    self.add_scatterers(scatterers)

  def set_occupancy(self, i, value):
    self.scatterers()[i].occupancy = value
    self.apply_symmetry(i, update_special_position_indices=00000)

  def shift_occupancy(self, i, delta):
    self.scatterers()[i].occupancy += delta
    self.apply_symmetry(i, update_special_position_indices=00000)

  def structure_factors(self, anomalous_flag=None, d_min=None,
                              algorithm=None,
                              cos_sin_table=00000,
                              quality_factor=None,
                              u_extra=None,
                              b_extra=None):
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors.from_scatterers(
      crystal_symmetry=self,
      d_min=d_min,
      cos_sin_table=cos_sin_table,
      quality_factor=quality_factor,
      u_extra=u_extra,
      b_extra=b_extra)(
        xray_structure=self,
        miller_set=miller_set,
        algorithm=algorithm)

  def show_summary(self, f=sys.stdout):
    print >> f, "Number of scatterers:", self.scatterers().size()
    print >> f, "At special positions:", self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Label  M  Coordinates            Occ  Uiso or Ucart"
    for scatterer in self.scatterers():
      scatterer.show(f=f, unit_cell=self.unit_cell())
    return self

  def apply_special_position_ops_d_target_d_site(self, d_target_d_site):
    for i in self.special_position_indices():
      site_symmetry = self.site_symmetry(self.scatterers()[i].site)
      r = matrix.sqr(site_symmetry.special_op().r().as_double())
      d_target_d_site[i] = (matrix.row(d_target_d_site[i]) * r).elems

  def asymmetric_unit_in_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_dict=self._scattering_dict)
    new_structure._scatterers = self.scatterers().deep_copy()
    new_structure._special_position_indices = \
      self.special_position_indices().deep_copy()
    return new_structure

  def expand_to_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_dict=self._scattering_dict)
    for scatterer in self.scatterers():
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      site_symmetry = self.site_symmetry(scatterer.site)
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      new_scatterer = scatterer.copy()
      for site in equiv_sites.coordinates():
        new_scatterer.site = site
        new_structure.add_scatterer(new_scatterer)
    return new_structure

  def change_basis(self, cb_op):
    new_structure = structure(
      crystal.special_position_settings.change_basis(self, cb_op),
      scattering_dict=self._scattering_dict)
    for scatterer in self.scatterers():
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      new_structure.add_scatterer(scatterer.copy(site=cb_op(scatterer.site)))
    return new_structure

  def change_hand(self):
    ch_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(ch_op)

  def apply_shift(self, shift):
    shifted_scatterers = self.scatterers().deep_copy()
    shifted_scatterers.set_sites(
      shifted_scatterers.extract_sites() + shift)
    return structure(
      special_position_settings=self,
      scatterers=shifted_scatterers,
      scattering_dict=self._scattering_dict)

  def sort(self, by_value="occupancy", reverse=00000):
    assert by_value in ("occupancy",)
    assert reverse in (00000, 0001)
    p = flex.sort_permutation(
      self.scatterers().extract_occupancies(),
      reverse)
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().select(p),
      scattering_dict=self._scattering_dict)

  def as_emma_model(self):
    from cctbx import euclidean_model_matching as emma
    positions = []
    for scatterer in self.scatterers():
      positions.append(emma.position(scatterer.label, scatterer.site))
    return emma.model(self, positions)

  def atomic_weights(self):
    from cctbx.eltbx import tiny_pse
    result = flex.double()
    for scatterer in self.scatterers():
      scattering_type = scatterer.scattering_type
      try:
        label = eltbx.caasf.wk1995(scattering_type, 1).label()
      except RuntimeError:
        raise RuntimeError("Unknown atomic weight: " + scattering_type)
      result.append(tiny_pse.table(label).weight())
    return result

  def center_of_mass(self, atomic_weights=None):
    if (atomic_weights is None):
      atomic_weights = self.atomic_weights()
    sites_cart = (self.unit_cell().orthogonalization_matrix()
                  * self.scatterers().extract_sites())
    sum_w = 0
    sum_wc = matrix.col((0,0,0))
    for i,site_cart in sites_cart.items():
      w = atomic_weights[i]
      sum_w += w
      sum_wc += matrix.col(site_cart) * w
    if (sum_w == 0): return sum_wc
    return sum_wc / sum_w

  def n_parameters(self, gradient_flags):
    n_scatterers = self.scatterers().size()
    n_anisotropic = self.scatterers().count_anisotropic()
    n_isotropic = n_scatterers - n_anisotropic
    result = 0
    if (gradient_flags.site): result += n_scatterers * 3
    if (gradient_flags.u_iso): result += n_isotropic
    if (gradient_flags.u_aniso): result += n_anisotropic * 6
    if (gradient_flags.occupancy): result += n_scatterers
    if (gradient_flags.fp): result += n_scatterers
    if (gradient_flags.fdp): result += n_scatterers
    return result
