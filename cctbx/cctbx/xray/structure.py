from cctbx.xray import ext
from cctbx.xray import structure_factors
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import matrix
from cctbx.array_family import flex
import types
import sys

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings, scatterers=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    if (scatterers is None):
      self._scatterers = flex.xray_scatterer()
      self._special_position_indices = flex.size_t()
    else:
      self._scatterers = scatterers.deep_copy()
      self.all_apply_symmetry()

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._special_position_indices = other._special_position_indices

  def deep_copy_scatterers(self):
    cp = structure(self)
    cp._scatterers = self._scatterers.deep_copy()
    cp._special_position_indices = self._special_position_indices.deep_copy()
    return cp

  def scatterers(self):
    return self._scatterers

  def special_position_indices(self):
    return self._special_position_indices

  def __getitem__(self, slice_object):
    assert type(slice_object) == types.SliceType
    assert self.scatterers() is not None
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().__getitem__(slice_object))

  def add_scatterer(self, scatterer):
    i = self.scatterers().size()
    self._scatterers.append(scatterer)
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

  def set_occupancy(self, i, value):
    self.scatterers()[i].occupancy = value
    self.apply_symmetry(i, update_special_position_indices=00000)

  def shift_occupancy(self, i, delta):
    self.scatterers()[i].occupancy += delta
    self.apply_symmetry(i, update_special_position_indices=00000)

  def structure_factors(self, anomalous_flag=None, d_min=None,
                              direct=00000, fft=00000):
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors.from_scatterers(
      crystal_symmetry=self,
      d_min=d_min)(
        xray_structure=self,
        miller_set=miller_set,
        direct=direct,
        fft=fft)

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

  def expand_to_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)))
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
      crystal.special_position_settings.change_basis(self, cb_op))
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
      scatterers=shifted_scatterers)

  def sort(self, by_value="occupancy", reverse=00000):
    assert by_value in ("occupancy",)
    assert reverse in (00000, 0001)
    p = flex.sort_permutation(
      self.scatterers().extract_occupancies(),
      reverse)
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().select(p))

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
      label = scatterer.caasf.label()
      assert label != "const", "Unknown atomic weight: "+scatterer.label
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
