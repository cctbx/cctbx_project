from cctbx.xray import ext
from cctbx.xray import structure_factors
from cctbx import miller
from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import sgtbx
import cctbx.eltbx.xray_scattering
from cctbx import adptbx
from cctbx import eltbx
from cctbx.array_family import flex
from scitbx import matrix
from stdlib import math
import types
import sys
import random

class structure(crystal.special_position_settings):

  def __init__(self,
        special_position_settings=None,
        scatterers=None,
        site_symmetry_table=None,
        scattering_type_registry=None,
        crystal_symmetry=None):
    assert [special_position_settings, crystal_symmetry].count(None) == 1
    assert scatterers is not None or site_symmetry_table is None
    if (special_position_settings is None):
      special_position_settings = crystal.special_position_settings(
        crystal_symmetry=crystal_symmetry)
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self.erase_scatterers()
    self._scattering_type_registry = scattering_type_registry
    if (scatterers is not None):
      self.add_scatterers(
        scatterers=scatterers,
        site_symmetry_table=site_symmetry_table)
    self.u_cart_group = None

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._site_symmetry_table = other._site_symmetry_table
    self._scattering_type_registry = other._scattering_type_registry
    self._scattering_type_registry_is_out_of_date \
      = other._scattering_type_registry_is_out_of_date

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell = self.unit_cell(),
      space_group_info = self.space_group_info())

  def erase_scatterers(self):
    self._scatterers = flex.xray_scatterer()
    self._site_symmetry_table = sgtbx.site_symmetry_table()
    self._scattering_type_registry_is_out_of_date = True

  def deep_copy_scatterers(self):
    cp = structure(self,
      scattering_type_registry=self._scattering_type_registry)
    cp._scatterers = self._scatterers.deep_copy()
    cp._site_symmetry_table = self._site_symmetry_table.deep_copy()
    if (getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    return cp

  def scatterers(self):
    return self._scatterers

  def approx_equal(self, other):
    assert self.scatterers().size() == other.scatterers().size()
    atom_atom_distances = flex.sqrt(self.difference_vectors_cart(other).dot())
    assert flex.mean_default(atom_atom_distances,0) < 1.e-6
    u_iso_1 = self.scatterers().extract_u_iso()
    u_iso_2 = other.scatterers().extract_u_iso()
    assert flex.mean(flex.abs(u_iso_1-u_iso_2)) < 1.e-6

  def set_u_iso(self, value=None, values=None, allow_mixed=False):
    assert [value, values].count(None) == 1
    s = self._scatterers
    if (not allow_mixed and s.count_anisotropic() > 0):
      raise RuntimeError("set_u_iso: all scatterers must be isotropic.")
    if (value is not None):
      s.set_u_iso(flex.double(s.size(), value))
    else:
      assert values.size() == s.size()
      s.set_u_iso(values)
    return self

  def set_b_iso(self, value=None, values=None, allow_mixed=False):
    assert [value, values].count(None) == 1
    s = self._scatterers
    if (not allow_mixed and s.count_anisotropic() > 0):
      raise RuntimeError("set_b_iso: all scatterers must be isotropic.")
    if (value is not None):
      s.set_u_iso(flex.double(s.size(), adptbx.b_as_u(value)))
    else:
      assert values.size() == s.size()
      b_iso = values
      u_iso_values = b_iso*adptbx.b_as_u(1)
      s.set_u_iso(u_iso_values)

  def random_remove_sites_selection(self, fraction):
    scatterers_size = self._scatterers.size()
    if(abs(fraction-0.0) < 1.e-3):
       return flex.bool(scatterers_size, True)
    if(fraction < 0.01 or fraction > 0.99):
       raise RuntimeError("fraction must be between 0.01 and 0.99.")
    tol = 999.
    selection = None
    l = max(fraction - 0.05, 0.0)
    r = min(fraction + 0.05, 1.0)
    for i in xrange(5):
      while l <= r:
        arr = flex.random_double(scatterers_size)-l
        sel = arr > 0.0
        deleted = float((scatterers_size - sel.count(True))) / scatterers_size
        if abs(fraction - deleted) < tol:
           tol = abs(fraction - deleted)
           selection = sel
        l += 0.0001
    return selection

  def replace_sites_cart(self, new_sites):
    cp = structure(self,
      scattering_type_registry=self._scattering_type_registry)
    new_scatterers = self._scatterers.deep_copy()
    sites_frac_new = self.unit_cell().fractionalization_matrix()*new_sites
    new_scatterers.set_sites(sites_frac_new)
    cp._scatterers = new_scatterers
    cp._site_symmetry_table = self._site_symmetry_table.deep_copy()
    if(getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    return cp

  def translate(self, x=0, y=0, z=0):
    sites_cart = self.sites_cart()
    sites_cart_size = sites_cart.size()
    shift_vector = flex.vec3_double(sites_cart_size,[x,y,z])
    sites_cart_new = sites_cart + shift_vector
    cp = structure(self,
      scattering_type_registry=self._scattering_type_registry)
    new_scatterers = self._scatterers.deep_copy()
    sites_frac_new = self.unit_cell().fractionalization_matrix()*sites_cart_new
    new_scatterers.set_sites(sites_frac_new)
    cp._scatterers = new_scatterers
    cp._site_symmetry_table = self._site_symmetry_table.deep_copy()
    if(getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    return cp

  def shake_sites(self, mean_error):
    tolerance = 0.00005
    sites_cart = self.sites_cart()
    sites_cart_size = sites_cart.size()
    current_mean_error = 0.0
    tolerance_scale = 1./5
    if(mean_error >= 0.1 and mean_error < 1.0):
       left  = mean_error - mean_error*0.1
       right = mean_error + mean_error*0.3
    elif(mean_error >= 1.0 and mean_error <= 3.0):
       left  = mean_error - mean_error*0.1
       right = mean_error + mean_error*1.0
    elif(abs(mean_error-0.0) < 1.e-3):
       return self
    else:
       raise RuntimeError("mean_error requested is too big or too small")
    while abs(mean_error - current_mean_error) > tolerance:
      two_left = 2.0 * left
      shift_xyz = (flex.random_double(sites_cart_size*3) - 0.5) * two_left
      sites_cart_new = sites_cart + flex.vec3_double(shift_xyz)
      left += tolerance * tolerance_scale
      #current_mean_error = sites_cart.rms_difference(sites_cart_new)
      # not the same in my definition
      current_mean_error = \
                      flex.mean(flex.sqrt((sites_cart - sites_cart_new).dot()))
      if(left >= right):
        raise RuntimeError("mean_error is not achieved within specified range")
    cp = structure(self,
      scattering_type_registry=self._scattering_type_registry)
    new_scatterers = self._scatterers.deep_copy()
    sites_frac_new = self.unit_cell().fractionalization_matrix()*sites_cart_new
    new_scatterers.set_sites(sites_frac_new)
    cp._scatterers = new_scatterers
    cp._site_symmetry_table = self._site_symmetry_table.deep_copy()
    if(getattr(self, "scatterer_pdb_records", None) is not None):
      cp.scatterer_pdb_records = self.scatterer_pdb_records
    #assert abs(sites_cart.rms_difference(cp.sites_cart())-mean_error) <= \
    #                                                                  tolerance
    assert abs(flex.mean(flex.sqrt((sites_cart - sites_cart_new).dot()))- \
                                                       mean_error) <= tolerance
    return cp

  def mean_distance(self, other):
    s1 = self.sites_cart()
    s2 = other.sites_cart()
    if(s1.size() != s2.size()):
       raise RuntimeError("models must be exactly aligned and of equal size.")
    return flex.mean(flex.sqrt((s1 - s2).dot()))

  def max_distance(self, other):
    s1 = self.sites_cart()
    s2 = other.sites_cart()
    if(s1.size() != s2.size()):
       raise RuntimeError("models must be exactly aligned and of equal size.")
    return flex.max(flex.sqrt((s1 - s2).dot()))

  def min_distance(self, other):
    s1 = self.sites_cart()
    s2 = other.sites_cart()
    if(s1.size() != s2.size()):
       raise RuntimeError("models must be exactly aligned and of equal size.")
    return flex.min(flex.sqrt((s1 - s2).dot()))

  def randomize_adp(self, random_u_iso_scale=1.0, random_u_cart_scale=1.0):
    for sc in self._scatterers:
        if(sc.flags.use()):
           if(sc.flags.use_u_iso()):
              u_iso = random.random() * random_u_iso_scale
              sc.u_iso = u_iso
           if(sc.flags.use_u_aniso()):
              site_symmetry = sc.apply_symmetry(self.unit_cell(),
                                                self.space_group(),
                                                self.min_distance_sym_equiv())
              run_away_counter = 0
              while 1:
                 run_away_counter += 1
                 assert run_away_counter < 100
                 u_cart = adptbx.random_u_cart(u_scale = random_u_cart_scale)
                 sc.u_star = site_symmetry.average_u_star(
                             adptbx.u_cart_as_u_star(self.unit_cell(), u_cart))
                 u_cart = adptbx.u_star_as_u_cart(self.unit_cell(), sc.u_star)
                 eigenvalues = adptbx.eigenvalues(u_cart)
                 if(min(eigenvalues) > 0.001): break

  def set_b_iso_random(self, allow_mixed=False):
    s = self._scatterers
    if (not allow_mixed and s.count_anisotropic() > 0):
      raise RuntimeError("set_b_iso_random: all scatterers must be isotropic.")
    b_iso_new = flex.random_double(s.size())*100.
    self.set_b_iso(values = b_iso_new)

  def shake_b_iso(self, deviation, allow_mixed=False):
    assert deviation >= 0.0 and deviation <= 100.0
    s = self._scatterers
    if (not allow_mixed and s.count_anisotropic() > 0):
      raise RuntimeError("shake_b_iso: all scatterers must be isotropic.")
    b_isos = s.extract_u_iso()/adptbx.b_as_u(1)
    shift_abs = b_isos * deviation/100.
    r_set = flex.random_double(s.size())-0.5
    b_iso_shifted = flex.double()
    for i in xrange(s.size()):
      if(r_set[i] >= 0): sign = 1.0
      if(r_set[i] <  0): sign =-1.0
      shift = shift_abs[i] * sign
      new_value = b_isos[i] + shift
      if(new_value > 0.0):
        b_iso_shifted.append(new_value)
      else:
        b_iso_shifted.append(b_isos[i])
    self.set_b_iso(values = b_iso_shifted)

  def n_undefined_multiplicities(self):
    return ext.n_undefined_multiplicities(self._scatterers)

  def sites_frac(self):
    return self._scatterers.extract_sites()

  def set_sites_frac(self, sites_frac):
    assert sites_frac.size() == self._scatterers.size()
    self._scatterers.set_sites(sites_frac)

  def sites_cart(self):
    return self.unit_cell().orthogonalization_matrix() * self.sites_frac()

  def set_sites_cart(self, sites_cart):
    self.set_sites_frac(
      self.unit_cell().fractionalization_matrix() * sites_cart)

  def extract_u_iso_or_u_equiv(self):
    return self._scatterers.extract_u_iso_or_u_equiv(
      unit_cell=self.unit_cell())

  def convert_to_isotropic(self):
    self._scatterers.convert_to_isotropic(unit_cell=self.unit_cell())

  def convert_to_anisotropic(self):
    self._scatterers.convert_to_anisotropic(unit_cell=self.unit_cell())

  def show_u_statistics(self, text="", out=None):
    if(out is None): out = sys.stdout
    size = self._scatterers.size()
    n_anisotropic = self._scatterers.count_anisotropic()
    n_isotropic = size - n_anisotropic
    ipd = self.is_positive_definite_u()
    npd = ipd.count(True)
    nnpd = ipd.count(False)
    bisos = self.extract_u_iso_or_u_equiv() * 8*math.pi**2
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print >> out, part1 + "-"*n + part2
    part1 = "| Total number of ADPs = %-d%-s"%(size,":")
    n = 79 - len(part1 + "|")
    print >> out, part1+" "*n+"|"
    print >> out, "|   iso = %-5d   aniso = %-5d   pos. def. = %-5d   "\
                "non-pos. def. = %-5d   |"%(n_isotropic,n_anisotropic,npd,nnpd)
    print >> out, "| Equivalent B: min = %-6.2f   max = %-6.2f   mean ="\
        " %-6.2f"%(flex.min(bisos),flex.max(bisos),flex.mean(bisos))+" "*19+"|"
    print >> out, "|" +"-"*77+"|"

  def site_symmetry_table(self):
    return self._site_symmetry_table

  def special_position_indices(self):
    return self._site_symmetry_table.special_position_indices()

  def scattering_type_registry(self,
        custom_dict=None,
        d_min=None,
        table=None,
        types_without_a_scattering_contribution=None):
    assert table in [None, "n_gaussian", "it1992", "wk1995"]
    if (table == "it1992"): assert d_min in [0,None] or d_min >= 1/4.
    if (table == "wk1995"): assert d_min in [0,None] or d_min >= 1/12.
    if (   self._scattering_type_registry_is_out_of_date
        or custom_dict is not None
        or d_min is not None
        or table is not None
        or types_without_a_scattering_contribution is not None):
      new_dict = {"const": eltbx.xray_scattering.gaussian(1) }
      old_dict = {}
      if (self._scattering_type_registry is not None):
        ugs = self._scattering_type_registry.unique_gaussians_as_list()
        tip = self._scattering_type_registry.type_index_pairs_as_dict()
        for t,i in tip.items():
          if (ugs[i] is not None): old_dict[t] = ugs[i]
        if (d_min is None and table is None):
          new_dict.update(old_dict)
      if (types_without_a_scattering_contribution is not None):
        for t in types_without_a_scattering_contribution:
          new_dict[t] = eltbx.xray_scattering.gaussian(0)
      if (custom_dict is not None):
        new_dict.update(custom_dict)
      if (d_min is None): d_min = 0
      self._scattering_type_registry = ext.scattering_type_registry()
      self._scattering_type_registry.process(self._scatterers)
      for t_undef in self._scattering_type_registry.unassigned_types():
        val = new_dict.get(t_undef, None)
        if (val is None):
          try:
            if (table == "it1992"):
              val = eltbx.xray_scattering.it1992(t_undef, True).fetch()
            elif (table == "wk1995"):
              val = eltbx.xray_scattering.wk1995(t_undef, True).fetch()
            else:
              if (t_undef == "D"): t_undef = "H"
              val = eltbx.xray_scattering.n_gaussian_table_entry(
                t_undef, d_min, 0).gaussian()
          except RuntimeError:
            val = old_dict.get(t_undef, None)
            if (val is None): raise
        self._scattering_type_registry.assign(t_undef, val)
      self._scattering_type_registry_is_out_of_date = False
    return self._scattering_type_registry

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
      scattering_type_registry=self._scattering_type_registry)

  def select(self, selection, negate=False):
    assert self.scatterers() is not None
    if (negate): selection = ~selection
    return structure(
      special_position_settings=self,
      scatterers=self._scatterers.select(selection),
      site_symmetry_table=self._site_symmetry_table.select(selection),
      scattering_type_registry=self._scattering_type_registry)


  def select_inplace(self, selection):
    assert self.scatterers() is not None
    self._scatterers          = self._scatterers.select(selection)
    self._site_symmetry_table = self._site_symmetry_table.select(selection)
    self._scattering_type_registry_is_out_of_date = True

  def add_scatterer(self, scatterer, site_symmetry_ops=None):
    self._scatterers.append(scatterer)
    if (site_symmetry_ops is None):
      site_symmetry_ops = self._scatterers[-1].apply_symmetry(
        unit_cell=self.unit_cell(),
        space_group=self.space_group(),
        min_distance_sym_equiv=self.min_distance_sym_equiv(),
        u_star_tolerance=self.u_star_tolerance(),
        assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    else:
      self._scatterers[-1].apply_symmetry(
        site_symmetry_ops=site_symmetry_ops,
        u_star_tolerance=self.u_star_tolerance())
    self._site_symmetry_table.process(site_symmetry_ops)
    self._scattering_type_registry_is_out_of_date = True

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
      assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    self._scattering_type_registry_is_out_of_date = True

  def concatenate(self, other):
    result = self.deep_copy_scatterers()
    result.add_scatterers(
      scatterers=other._scatterers,
      site_symmetry_table=other._site_symmetry_table)
    return result

  def replace_scatterers(self, scatterers, site_symmetry_table="existing"):
    if (site_symmetry_table == "existing"):
      site_symmetry_table = self._site_symmetry_table
    if (site_symmetry_table is not None):
      assert site_symmetry_table.indices().size() == self.scatterers().size()
    self.erase_scatterers()
    self.add_scatterers(
      scatterers=scatterers,
      site_symmetry_table=site_symmetry_table)

  def structure_factors(self, anomalous_flag=None, d_min=None,
                              algorithm=None,
                              cos_sin_table=False,
                              quality_factor=None,
                              u_base=None,
                              b_base=None,
                              wing_cutoff=None):
    if (anomalous_flag is None):
      if (self.scatterers().count_anomalous() != 0):
        anomalous_flag = True
      else:
        anomalous_flag = False
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
      b_base=b_base,
      wing_cutoff=wing_cutoff)(
        xray_structure=self,
        miller_set=miller_set,
        algorithm=algorithm)

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix + "Number of scatterers:", \
      self.scatterers().size()
    print >> f, prefix + "At special positions:", \
      self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f=f, prefix=prefix)
    return self

  def show_scatterers(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso"
    for scatterer in self.scatterers():
      scatterer.show(f=f, unit_cell=self.unit_cell())
    return self

  def show_special_position_shifts(self,
        sites_frac_original=None,
        sites_cart_original=None,
        out=None,
        prefix=""):
    self._site_symmetry_table.show_special_position_shifts(
      special_position_settings=self,
      site_labels=self.scatterers().extract_labels(),
      sites_frac_original=sites_frac_original,
      sites_cart_original=sites_cart_original,
      sites_frac_exact=self.scatterers().extract_sites(),
      out=out,
      prefix=prefix)

  def is_positive_definite_u(self, u_cart_tolerance=None):
    if (u_cart_tolerance is None):
      return ext.is_positive_definite_u(
        scatterers=self._scatterers,
        unit_cell=self.unit_cell())
    else:
      return ext.is_positive_definite_u(
        scatterers=self._scatterers,
        unit_cell=self.unit_cell(),
        u_cart_tolerance=u_cart_tolerance)

  def tidy_us(self, u_min):
    ext.tidy_us(
      scatterers=self._scatterers,
      unit_cell=self.unit_cell(),
      site_symmetry_table=self._site_symmetry_table,
      u_min=u_min)

  def shift_us(self, u_shift=None, b_shift=None):
    assert [u_shift, b_shift].count(None) == 1
    if (u_shift is None):
      u_shift = adptbx.b_as_u(b_shift)
    ext.shift_us(
      scatterers=self._scatterers,
      unit_cell=self.unit_cell(),
      u_shift=u_shift)

  def apply_symmetry_sites(self):
    ext.apply_symmetry_sites(
      site_symmetry_table=self._site_symmetry_table,
      scatterers=self._scatterers)

  def apply_symmetry_u_stars(self):
    ext.apply_symmetry_u_stars(
      site_symmetry_table=self._site_symmetry_table,
      scatterers=self._scatterers,
      u_star_tolerance=self.u_star_tolerance())

  def re_apply_symmetry(self, i_scatterer):
    self._scatterers[i_scatterer].apply_symmetry(
      site_symmetry_ops=self._site_symmetry_table.get(i_scatterer),
      u_star_tolerance=self.u_star_tolerance())

  def apply_special_position_ops_d_target_d_site(self, d_target_d_site):
    for i in self.special_position_indices():
      special_op = self._site_symmetry_table.get(i).special_op()
      r = matrix.sqr(special_op.r().as_double())
      d_target_d_site[i] = (matrix.row(d_target_d_site[i]) * r).elems

  def asymmetric_unit_in_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)),
      scattering_type_registry=self._scattering_type_registry)
    new_structure._scatterers = self.scatterers().deep_copy()
    new_structure._site_symmetry_table = self.site_symmetry_table().deep_copy()
    return new_structure

  def change_basis(self, cb_op):
    return structure(
      special_position_settings
        =crystal.special_position_settings.change_basis(self, cb_op),
      scatterers=ext.change_basis(scatterers=self._scatterers, cb_op=cb_op),
      site_symmetry_table=self._site_symmetry_table.change_basis(cb_op=cb_op),
      scattering_type_registry=self._scattering_type_registry)

  def change_hand(self):
    ch_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(ch_op)

  def expand_to_p1(self, append_number_to_labels=False):
    return structure(
      special_position_settings
        =crystal.special_position_settings(
          crystal.symmetry.cell_equivalent_p1(self)),
      scatterers=ext.expand_to_p1(
        unit_cell=self.unit_cell(),
        space_group=self.space_group(),
        scatterers=self._scatterers,
        site_symmetry_table=self._site_symmetry_table,
        append_number_to_labels=append_number_to_labels),
      scattering_type_registry=self._scattering_type_registry)

  def apply_shift(self, shift, recompute_site_symmetries=False):
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
      scattering_type_registry=self._scattering_type_registry)

  def random_shift_sites(self, max_shift_cart=0.2):
    shifts = flex.vec3_double(
      (flex.random_double(self.scatterers().size()*3)*2-1) * max_shift_cart)
    return self.apply_shift(self.unit_cell().fractionalization_matrix()*shifts)

  def sort(self, by_value="occupancy", reverse=False):
    assert by_value in ("occupancy",)
    assert reverse in (False, True)
    p = flex.sort_permutation(
      data=self.scatterers().extract_occupancies(),
      reverse=reverse)
    return structure(
      special_position_settings=self,
      scatterers=self._scatterers.select(p),
      site_symmetry_table=self._site_symmetry_table.select(p),
      scattering_type_registry=self._scattering_type_registry)

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
    return self.sites_cart().mean_weighted(weights=atomic_weights)


  def n_parameters(self):
    #XXX move to C++ (after anisotropic_flag is gone)
    result_ = 0
    for sc in self.scatterers():
        if(sc.flags.grad_site()     ): result_ +=3
        if(sc.flags.grad_u_iso()   and sc.anisotropic_flag==False): result_ +=1
        if(sc.flags.grad_u_aniso() and sc.anisotropic_flag==True): result_ +=6
        if(sc.flags.grad_occupancy()): result_ +=1
        if(sc.flags.grad_fp()       ): result_ +=1
        if(sc.flags.grad_fdp()      ): result_ +=1
    return result_

  def n_parameters_XXX(self):
    #XXX move to C++ (after anisotropic_flag is gone)
    result_ = 0
    for sc in self.scatterers():
        if(sc.flags.grad_site()): result_ +=3
        if(sc.flags.grad_u_iso() and sc.flags.use_u_iso()): result_ +=1
        if(sc.flags.grad_u_aniso() and sc.flags.use_u_aniso()): result_ +=6
        if(sc.flags.grad_occupancy()): result_ +=1
        if(sc.flags.grad_fp()       ): result_ +=1
        if(sc.flags.grad_fdp()      ): result_ +=1
    return result_

  def asu_mappings(self, buffer_thickness, is_inside_epsilon=None):
    result = crystal.symmetry.asu_mappings(self,
      buffer_thickness=buffer_thickness,
      is_inside_epsilon=is_inside_epsilon)
    ext.asu_mappings_process(
      asu_mappings=result,
      scatterers=self._scatterers,
      site_symmetry_table=self._site_symmetry_table)
    return result

  def pair_asu_table(self,
        distance_cutoff=None,
        asu_mappings_buffer_thickness=None,
        asu_mappings_is_inside_epsilon=None):
    assert distance_cutoff is not None or asu_mappings_buffer_thickness is not None
    if (asu_mappings_buffer_thickness is None):
      asu_mappings_buffer_thickness = distance_cutoff
    asu_mappings = self.asu_mappings(
      buffer_thickness=asu_mappings_buffer_thickness)
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    if (distance_cutoff is not None):
      pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
    return pair_asu_table

  def show_distances(self,
        distance_cutoff=None,
        asu_mappings_buffer_thickness=None,
        asu_mappings_is_inside_epsilon=None,
        pair_asu_table=None,
        show_cartesian=False,
        keep_pair_asu_table=False,
        out=None):
    assert [distance_cutoff, pair_asu_table].count(None) == 1
    if (pair_asu_table is None):
      pair_asu_table = self.pair_asu_table(
        distance_cutoff=distance_cutoff,
        asu_mappings_buffer_thickness=asu_mappings_buffer_thickness,
        asu_mappings_is_inside_epsilon=asu_mappings_is_inside_epsilon)
    return pair_asu_table.show_distances(
      site_labels=self.scatterers().extract_labels(),
      sites_frac=self.sites_frac(),
      show_cartesian=show_cartesian,
      keep_pair_asu_table=keep_pair_asu_table,
      out=out)

  def difference_vectors_cart(self, other):
    return other.sites_cart() - self.sites_cart()

  def rms_difference(self, other):
    return self.sites_cart().rms_difference(other.sites_cart())
