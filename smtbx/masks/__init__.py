import cctbx.masks
from cctbx import maptbx
from cctbx import miller
from cctbx import xray
from cctbx.array_family import flex

class mask(object):
  def __init__(self, xray_structure, observations):
    self.xray_structure = xray_structure
    self.observations = observations
    self.mask = None
    self._f_mask = None
    self.masked_diff_map = None
    self.f_000 = None
    self.f_000_s = None
    self.f_000_cell = None

  def compute(self,
              solvent_radius,
              shrink_truncation_radius,
              ignore_hydrogen_atoms=False,
              crystal_gridding=None,
              resolution_factor=1./4,
              atom_radii=None):
    xs = self.xray_structure
    if crystal_gridding is None:
      self.crystal_gridding = maptbx.crystal_gridding(
        unit_cell=xs.unit_cell(),
        space_group_info=self.xray_structure.space_group_info(),
        d_min=self.observations.d_min(),
        resolution_factor=resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry)
    else:
      self.crystal_gridding = crystal_gridding
    if atom_radii is None:
      atom_radii = cctbx.masks.vdw_radii_from_xray_structure(xs)
    xs_p1 = self.xray_structure.expand_to_p1()
    self.mask = cctbx.masks.around_atoms(
      unit_cell=xs_p1.unit_cell(),
      space_group_order_z=xs_p1.space_group().order_z(),
      sites_frac=xs_p1.sites_frac(),
      atom_radii=cctbx.masks.vdw_radii_from_xray_structure(xs_p1),
      gridding_n_real=self.crystal_gridding.n_real(),
      solvent_radius=solvent_radius,
      shrink_truncation_radius=shrink_truncation_radius)
    self.solvent_accessible_volume = float(
      self.mask.data.count(1))/self.mask.data.size() * self.xray_structure.unit_cell().volume()
    print "Solvent accessible volume = %.1f [%.1f%%]" %(
      self.solvent_accessible_volume, 100.*
      self.solvent_accessible_volume/xs.unit_cell().volume())

  def structure_factors(self, max_cycles, scale_factor=None):
    """P. van der Sluis and A. L. Spek, Acta Cryst. (1990). A46, 194-201."""
    assert self.mask is not None
    f_obs = self.observations.as_amplitude_array()
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_obs,
      cos_sin_table=True)
    self.f_calc = sf(self.xray_structure, f_obs).f_calc()
    if scale_factor is None:
      self.scale_factor = f_obs.quick_scale_factor_approximation(
        self.f_calc, cutoff_factor=0.8)
    f_obs_minus_f_calc = f_obs.f_obs_minus_f_calc(
      1./self.scale_factor, self.f_calc)
    self.fft_scale = self.xray_structure.unit_cell().volume()\
              / self.crystal_gridding.n_grid_points()
    for i in range(max_cycles):
      diff_map = miller.fft_map(self.crystal_gridding, f_obs_minus_f_calc, self.f_000_s)
      diff_map.apply_volume_scaling()
      stats = diff_map.statistics()
      masked_diff_map = diff_map.real_map_unpadded().set_selected(
        self.mask.data.as_double() == 0, 0)
      self.f_000_cell = flex.sum(diff_map.real_map_unpadded()) * self.fft_scale
      self.f_000 = flex.sum(masked_diff_map) * self.fft_scale
      self.f_000_s = self.f_000 * (
        float(masked_diff_map.size())/(masked_diff_map.size() - self.mask.data.count(1)))
      #print "F000 cell: %.1f" %self.f_000_cell
      #print "F000 void: %.1f" %self.f_000_s
      masked_diff_map += self.f_000_s/self.xray_structure.unit_cell().volume()
      self._f_mask = f_obs.structure_factors_from_map(map=masked_diff_map)
      self._f_mask *= self.fft_scale
      f_model = self.f_model()
      f_obs_minus_f_calc = f_obs.phase_transfer(f_model).f_obs_minus_f_calc(
        1./self.scale_factor, self.f_calc)
    print "F000 cell: %.1f" %self.f_000_cell
    print "F000 void: %.1f" %self.f_000_s
    self.masked_diff_map = masked_diff_map
    return self._f_mask

  def f_mask(self):
    assert self._f_mask is not None
    return self._f_mask

  def f_model(self):
    assert self._f_mask is not None
    data = self.f_calc.data()+self.f_mask().data()
    return miller.array(miller_set=self.f_calc, data=data)

  def modified_structure_factors(self):
    """Subtracts the solvent contribution from the observed structure
    factors to obtain modified structure factors, suitable for refinement
    with other refinement programs such as ShelXL"""
    f_obs = self.observations.as_amplitude_array()
    f_model = self.f_model()
    f_obs = f_obs.phase_transfer(phase_source=f_model)
    assert self._f_mask is not None
    f_mask = self.f_mask()
    modified_f_obs = miller.array(
      miller_set=f_obs,
      data=(f_obs.data() - f_mask.data()*self.scale_factor),
      sigmas=f_obs.sigmas())
    modified_fo_sq = modified_f_obs.as_intensity_array()
    return modified_fo_sq
