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

  def __init__(self,
        special_position_settings=None,
        scatterers=None,
        site_symmetry_table=None,
        scattering_dict=None,
        crystal_symmetry=None):
    assert [special_position_settings, crystal_symmetry].count(None) == 1
    assert scatterers is not None or site_symmetry_table is None
    if (special_position_settings is None):
      special_position_settings = crystal.special_position_settings(
        crystal_symmetry=crystal_symmetry)
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self.erase_scatterers()
    self._scattering_dict = scattering_dict
    if (scatterers is not None):
      self.add_scatterers(
        scatterers=scatterers,
        site_symmetry_table=site_symmetry_table)

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._site_symmetry_table = other._site_symmetry_table
    self._scattering_dict = other._scattering_dict
    self._scattering_dict_is_out_of_date=other._scattering_dict_is_out_of_date

  def erase_scatterers(self):
    self._scatterers = flex.xray_scatterer()
    self._site_symmetry_table = sgtbx.site_symmetry_table()
    self._scattering_dict_is_out_of_date = 0001

  def deep_copy_scatterers(self):
    cp = structure(self, scattering_dict=self._scattering_dict)
    cp._scatterers = self._scatterers.deep_copy()
    cp._site_symmetry_table = self._site_symmetry_table.deep_copy()
    if (getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    return cp

  def scatterers(self):
    return self._scatterers

  def sites_frac(self):
    return self.scatterers().extract_sites()

  def set_sites_frac(self, sites_frac):
    assert sites_frac.size() == self._scatterers.size()
    self._scatterers.set_sites(sites_frac)

  def sites_cart(self):
    return self.unit_cell().orthogonalization_matrix() * self.sites_frac()

  def set_sites_cart(self, sites_cart):
    self.set_sites_frac(
      self.unit_cell().fractionalization_matrix() * sites_cart)

  def site_symmetry_table(self):
    return self._site_symmetry_table

  def special_position_indices(self):
    return self._site_symmetry_table.special_position_indices()

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
    sel = flex.slice_indices(
      array_size=self._scatterers.size(),
      python_slice=slice_object)
    return structure(
      special_position_settings=self,
      scatterers=self._scatterers.select(sel),
      site_symmetry_table=self._site_symmetry_table.select(sel),
      scattering_dict=self._scattering_dict)

  def add_scatterer(self, scatterer, site_symmetry_ops=None):
    self._scatterers.append(scatterer)
    if (site_symmetry_ops is None):
      site_symmetry_ops = self._scatterers[-1].apply_symmetry(
        unit_cell=self.unit_cell(),
        space_group=self.space_group(),
        min_distance_sym_equiv=self.min_distance_sym_equiv(),
        u_star_tolerance=self.u_star_tolerance(),
        assert_is_positive_definite=self.assert_is_positive_definite(),
        assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    elif (not site_symmetry_ops.is_point_group_1()):
      self._scatterers[-1].apply_symmetry_site(
        site_symmetry_ops=site_symmetry_ops)
      self._scatterers[-1].apply_symmetry_u_star(
        unit_cell=self.unit_cell(),
        site_symmetry_ops=site_symmetry_ops,
        u_star_tolerance=self.u_star_tolerance(),
        assert_is_positive_definite=self.assert_is_positive_definite(),
        assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    self._site_symmetry_table.process(site_symmetry_ops)
    self._scattering_dict_is_out_of_date = 0001

  def add_scatterers(self, scatterers, site_symmetry_table=None):
    if (site_symmetry_table is None):
      site_symmetry_table = sgtbx.site_symmetry_table()
    else:
      assert site_symmetry_table.indices().size() == scatterers.size()
    self._scatterers.extend(scatterers)
    ext.add_scatterers_ext(
      unit_cell=self.unit_cell(),
      space_group=self.space_group(),
      scatterers=self._scatterers,
      site_symmetry_table=self._site_symmetry_table,
      site_symmetry_table_for_new=site_symmetry_table,
      min_distance_sym_equiv=self.min_distance_sym_equiv(),
      u_star_tolerance=self.u_star_tolerance(),
      assert_is_positive_definite=self.assert_is_positive_definite(),
      assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    self._scattering_dict_is_out_of_date = 0001

  def replace_scatterers(self, scatterers, site_symmetry_table=None):
    self.erase_scatterers()
    self.add_scatterers(
      scatterers=scatterers,
      site_symmetry_table=site_symmetry_table)

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
      special_op = self._site_symmetry_table.get(i).special_op()
      r = matrix.sqr(special_op.r().as_double())
      d_target_d_site[i] = (matrix.row(d_target_d_site[i]) * r).elems

  def asymmetric_unit_in_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_dict=self._scattering_dict)
    new_structure._scatterers = self.scatterers().deep_copy()
    new_structure._site_symmetry_table = self.site_symmetry_table().deep_copy()
    return new_structure

  def expand_to_p1(self, append_number_to_labels=00000):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_dict=self._scattering_dict)
    for i_seq,scatterer in enumerate(self.scatterers()):
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      if (append_number_to_labels):
        if (scatterer.multiplicity() >= 100):
          fmt = "_%03d"
        elif (scatterer.multiplicity() >= 10):
          fmt = "_%02d"
        else:
          fmt = "_%d"
      equiv_sites = sgtbx.sym_equiv_sites(
        unit_cell=self.unit_cell(),
        space_group=self.space_group(),
        original_site=scatterer.site,
        site_symmetry_ops=self._site_symmetry_table.get(i_seq))
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
    for i_seq,scatterer in enumerate(self.scatterers()):
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      new_structure.add_scatterer(
        scatterer.copy(site=cb_op(scatterer.site)),
        site_symmetry_ops=
          self._site_symmetry_table.get(i_seq).change_basis(cb_op))
    return new_structure

  def change_hand(self):
    ch_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(ch_op)

  def apply_shift(self, shift, recompute_site_symmetries=00000):
    shifted_scatterers = self.scatterers().deep_copy()
    shifted_scatterers.set_sites(shifted_scatterers.extract_sites() + shift)
    if (recompute_site_symmetries):
      site_symmetry_table = None
    else:
      site_symmetry_table = self._site_symmetry_table
    return structure(
      special_position_settings=self,
      scatterers=shifted_scatterers,
      site_symmetry_table=site_symmetry_table,
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
      scatterers=self._scatterers.select(p),
      site_symmetry_table=self._site_symmetry_table.select(p),
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

  def asu_mappings(self, buffer_thickness, is_inside_epsilon=None):
    result = crystal.direct_space_asu.asu_mappings(
      space_group=self.space_group(),
      asu=self.direct_space_asu().as_float_asu(
        is_inside_epsilon=is_inside_epsilon),
      buffer_thickness=buffer_thickness)
    ext.asu_mappings_process(
      asu_mappings=result,
      scatterers=self._scatterers,
      site_symmetry_table=self._site_symmetry_table)
    return result

  def difference_vectors_cart(self, other):
    return other.sites_cart() - self.sites_cart()

  def rms(self, other):
    d = self.difference_vectors_cart(other)
    return math.sqrt(flex.mean(d.dot(d)))
