from cctbx.xray import ext
from cctbx.xray import structure_factors
from cctbx import miller
from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import sgtbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.itertbx import count
from stdlib import math
import types
import sys

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings=None, scatterers=None,
                     scattering_dict=None, crystal_symmetry=None):
    assert [special_position_settings, crystal_symmetry].count(None) == 1
    if (special_position_settings is None):
      special_position_settings = crystal.special_position_settings(
        crystal_symmetry=crystal_symmetry)
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
    if (getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    return cp

  def scatterers(self):
    return self._scatterers

  def sites_frac(self):
    return self.scatterers().extract_sites()

  def sites_cart(self):
    return self.unit_cell().orthogonalization_matrix() * self.sites_frac()

  def special_position_indices(self):
    return self._special_position_indices

  def scattering_dict(self, custom_dict=None, d_min=None, table=None):
    assert table in [None, "n_gaussian", "it1992", "wk1995"]
    if (table == "it1992"): assert d_min in [0,None] or d_min >= 1/4.
    if (table == "wk1995"): assert d_min in [0,None] or d_min >= 1/12.
    if (   self._scattering_dict_is_out_of_date
        or custom_dict is not None
        or d_min is not None
        or table is not None):
      new_dict = {"const": eltbx.xray_scattering.gaussian(1) }
      if (    self._scattering_dict is not None
          and d_min is None
          and table is None):
        for k,v in self._scattering_dict.dict().items():
          new_dict[k] = v.gaussian
      if (custom_dict is not None):
        new_dict.update(custom_dict)
      if (d_min is None): d_min = 0
      self._scattering_dict = ext.scattering_dictionary(self.scatterers())
      for key_undef in self._scattering_dict.find_undefined():
        if (new_dict.has_key(key_undef)):
          val = new_dict[key_undef]
        else:
          if (table == "it1992"):
            val = eltbx.xray_scattering.it1992(key_undef, 1).fetch()
          elif (table == "wk1995"):
            val = eltbx.xray_scattering.wk1995(key_undef, 1).fetch()
          else:
            val = eltbx.xray_scattering.n_gaussian_table_entry(
              key_undef, d_min, 0).gaussian()
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
    self._scatterers.extend(scatterers)
    self._scattering_dict_is_out_of_date = 0001
    self._special_position_indices.extend(special_position_indices)

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
                              u_base=None,
                              b_base=None):
    if (anomalous_flag is None):
      if (self.scatterers().count_anomalous() != 0):
        anomalous_flag = 0001
      else:
        anomalous_flag = 00000
    elif (not anomalous_flag):
      if (self.scatterers().count_anomalous() != 0):
        raise RuntimeError(
            "xray.structure with anomalous scatterers"
          + " but miller.array is non-anomalous.")
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors.from_scatterers(
      crystal_symmetry=self,
      d_min=d_min,
      cos_sin_table=cos_sin_table,
      quality_factor=quality_factor,
      u_base=u_base,
      b_base=b_base)(
        xray_structure=self,
        miller_set=miller_set,
        algorithm=algorithm)

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Number of scatterers:", self.scatterers().size()
    print >> f, "At special positions:", self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso"
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

  def expand_to_p1(self, append_number_to_labels=00000):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_dict=self._scattering_dict)
    for scatterer in self.scatterers():
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      site_symmetry = self.site_symmetry(scatterer.site)
      if (append_number_to_labels):
        if (site_symmetry.multiplicity() >= 100):
          fmt = "_%03d"
        elif (site_symmetry.multiplicity() >= 10):
          fmt = "_%02d"
        else:
          fmt = "_%d"
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      new_scatterer = scatterer.copy()
      for i,site in zip(count(), equiv_sites.coordinates()):
        if (append_number_to_labels):
          new_scatterer.label = scatterer.label + fmt % i
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

  def random_shift_sites(self, max_shift_cart=0.2):
    shifts = flex.vec3_double(
      (flex.random_double(self.scatterers().size()*3)*2-1) * max_shift_cart)
    return self.apply_shift(self.unit_cell().fractionalization_matrix()*shifts)

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
        label = eltbx.xray_scattering.wk1995(scattering_type, 1).label()
      except RuntimeError:
        raise RuntimeError("Unknown atomic weight: " + scattering_type)
      result.append(tiny_pse.table(label).weight())
    return result

  def center_of_mass(self, atomic_weights=None):
    if (atomic_weights is None):
      atomic_weights = self.atomic_weights()
    sites_cart = self.sites_cart()
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

  def asu_mappings(self, buffer_thickness,
                         is_inside_epsilon=None,
                         sym_equiv_epsilon=1.e-6):
    result = crystal.direct_space_asu.asu_mappings(
      space_group=self.space_group(),
      asu=self.direct_space_asu().as_float_asu(
        is_inside_epsilon=is_inside_epsilon),
      buffer_thickness=buffer_thickness,
      sym_equiv_epsilon=sym_equiv_epsilon)
    ext.asu_mappings_process(
      asu_mappings=result,
      scatterers=self.scatterers())
    return result

  def difference_vectors_cart(self, other):
    return other.sites_cart() - self.sites_cart()

  def rms(self, other):
    d = self.difference_vectors_cart(other)
    return math.sqrt(flex.mean(d.dot(d)))
