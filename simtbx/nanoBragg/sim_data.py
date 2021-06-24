from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import shapetype, nanoBragg
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.nanoBragg_beam import NBbeam
from copy import deepcopy


def Amatrix_dials2nanoBragg(crystal):
  """
  returns the A matrix from a cctbx crystal object
  in nanoBragg format
  :param crystal: cctbx crystal
  :return: Amatrix as a tuple
  """
  sgi = crystal.get_space_group().info()
  if sgi.type().lookup_symbol().startswith('C'):
    raise ValueError('You need to convert your crystal model to its primitive setting first')
  Amatrix = sqr(crystal.get_A()).transpose()
  return Amatrix


def determine_spot_scale(beam_size_mm, crystal_thick_mm, mosaic_vol_mm3):
  """
  :param beam_size_mm:  diameter of beam focus (millimeter)
  :param crystal_thick_mm: thickness of crystal (millimeter)
  :param mosaic_vol_mm3:  volume of a mosaic block in crystal (cubic mm)
  :return: roughly the number of exposed mosaic blocks
  """
  if beam_size_mm <= crystal_thick_mm:
    illum_xtal_vol = crystal_thick_mm * beam_size_mm ** 2
  else:
    illum_xtal_vol = crystal_thick_mm ** 3
  return illum_xtal_vol / mosaic_vol_mm3


class SimData:

  def __init__(self, default_crystal=True):
    self.detector = SimData.simple_detector(180, 0.1, (512, 512))
    self.seed = 1
    self.crystal = NBcrystal(default=default_crystal)
    self.add_air = False
    self.add_water = True
    self.water_path_mm = 0.005
    self.air_path_mm = 0
    nbBeam = NBbeam()
    nbBeam.unit_s0 = (0, 0, -1)
    self.beam = nbBeam
    self.using_cuda = False
    self.using_omp = False
    self.rois = None
    self.readout_noise = 3
    self.gain = 1
    self.psf_fwhm = 0
    self.include_noise = True
    self.background_raw_pixels = None  # background raw pixels, should be a 2D flex double array
    self.backrground_scale = 1  # scale factor to apply to background raw pixels
    self.functionals = []
    self.mosaic_seeds = 777, 777
    self.D = None # nanoBragg instance

  @property
  def background_raw_pixels(self):
    return self._background_raw_pixels

  @background_raw_pixels.setter
  def background_raw_pixels(self, val):
    self._background_raw_pixels = val

  @property
  def gain(self):
    return self._gain

  @gain.setter
  def gain(self, val):
    self._gain = val

  @property
  def air_path_mm(self):
    return self._air_path_mm

  @air_path_mm.setter
  def air_path_mm(self, val):
    self._air_path_mm = val

  @property
  def water_path_mm(self):
    return self._water_path_mm

  @water_path_mm.setter
  def water_path_mm(self, val):
    self._water_path_mm = val

  @property
  def crystal(self):
    return self._crystal

  @crystal.setter
  def crystal(self, val):
    self._crystal = val

  @property
  def beam(self):
    return self._beam

  @beam.setter
  def beam(self, val):
    self._beam = val

  @property
  def seed(self):
    return self._seed

  @seed.setter
  def seed(self, val):
    self._seed = val

  @staticmethod
  def Umats(mos_spread_deg, n_mos_doms, isotropic=True, seed=777, norm_dist_seed=777):
    import scitbx
    import scitbx.math
    import math
    UMAT_nm = flex.mat3_double()
    mersenne_twister = flex.mersenne_twister(seed=seed)
    scitbx.random.set_random_seed(norm_dist_seed)
    rand_norm = scitbx.random.normal_distribution(mean=0, sigma=(mos_spread_deg * math.pi / 180.0))
    g = scitbx.random.variate(rand_norm)
    mosaic_rotation = g(n_mos_doms)
    sites, angles = [], []
    for m in mosaic_rotation:
      site = mersenne_twister.random_double_point_on_sphere()
      if mos_spread_deg > 0:
        sites.append(site)
        angles.append(m)
      else:
        sites.append(site)
        angles.append(0)
      if isotropic and mos_spread_deg > 0:
        sites.append(site)
        angles.append(-m)
    UMAT_nm = scitbx.math.r3_rotation_axis_and_angle_as_matrix(sites, angles)

    return UMAT_nm

  @property
  def psf_fwhm(self):
    return self._psf_fwhm

  @psf_fwhm.setter
  def psf_fwhm(self, val):
    self._psf_fwhm = val

  @property
  def readout_noise(self):
    return self._readout_noise

  @readout_noise.setter
  def readout_noise(self, val):
    self._readout_noise = val

  @property
  def add_air(self):
    return self._add_air

  @add_air.setter
  def add_air(self, val):
    self._add_air = val

  @property
  def add_water(self):
    return self._add_water

  @add_water.setter
  def add_water(self, val):
    self._add_water = val

  @property
  def detector(self):
    return self._detector

  @detector.setter
  def detector(self, val):
    self._detector = val

  @property
  def rois(self):
    return self._rois

  @rois.setter
  def rois(self, val):
    self._rois = val

  @property
  def using_omp(self):
    return self._using_omp

  @using_omp.setter
  def using_omp(self, val):
    assert val in (True, False)
    self._using_omp = val

  @property
  def using_cuda(self):
    return self._using_cuda

  @using_cuda.setter
  def using_cuda(self, val):
    assert val in (True, False)
    self._using_cuda = val

  @property
  def include_noise(self):
    return self._include_noise

  @include_noise.setter
  def include_noise(self, val):
    self._include_noise = val

  def update_Fhkl_tuple(self):
    if self.crystal.miller_array is not None:
      self.D.Fhkl_tuple = (
        self.crystal.miller_array.indices(), self.crystal.miller_array.data())

  def _crystal_properties(self):
    if self.crystal is not None:
      self.D.xtal_shape = self.crystal.xtal_shape
      self.update_Fhkl_tuple()
      self.D.Amatrix = Amatrix_dials2nanoBragg(self.crystal.dxtbx_crystal)
      self.D.Ncells_abc = self.crystal.Ncells_abc
      self.D.mosaic_spread_deg = self.crystal.mos_spread_deg
      self.D.mosaic_domains = self.crystal.n_mos_domains
      self.D.set_mosaic_blocks(SimData.Umats(self.crystal.mos_spread_deg, self.crystal.n_mos_domains,
                                             seed=self.mosaic_seeds[0], norm_dist_seed=self.mosaic_seeds[1]) )

  def _beam_properties(self):
    self.D.xray_beams = self.beam.xray_beams
    self.D.beamsize_mm = self.beam.size_mm

  def _seedlings(self):
    self.D.seed = self.seed
    self.D.calib_seed = self.seed
    self.D.mosaic_seed = self.seed

  def determine_spot_scale(self):
    if self.crystal is None:
      return 1
    if self.beam.size_mm <= self.crystal.thick_mm:
      illum_xtal_vol = self.crystal.thick_mm * self.beam.size_mm ** 2
    else:
      illum_xtal_vol = self.crystal.thick_mm ** 3
    mosaic_vol = self.D.xtal_size_mm[0] * self.D.xtal_size_mm[1] * self.D.xtal_size_mm[2]
    return illum_xtal_vol / mosaic_vol

  def update_nanoBragg_instance(self, parameter, value):
    setattr(self.D, parameter, value)

  @property
  def panel_id(self):
    return self._panel_id

  @panel_id.setter
  def panel_id(self, val):
    if val >= len(self.detector):
      raise ValueError("panel id cannot be larger than the number of panels in detector (%d)" % len(self.detector))
    if val <0:
      raise ValueError("panel id cannot be negative!")
    if self.D is None:
      self._panel_id = 0
    else:
      self.D.set_dxtbx_detector_panel(self.detector[int(val)], self.beam.nanoBragg_constructor_beam.get_s0())
      self._panel_id = int(val)

  def instantiate_nanoBragg(self, verbose=0, oversample=0, device_Id=0, adc_offset=0, default_F=1000.0, interpolate=0,
                            pid=0):
    self.D = nanoBragg(self.detector, self.beam.nanoBragg_constructor_beam, verbose=verbose,
                       panel_id=int(pid))
    self._seedlings()
    self.D.interpolate = interpolate
    self._crystal_properties()
    self._beam_properties()
    self.D.spot_scale = self.determine_spot_scale()
    self.D.adc_offset_adu = adc_offset
    self.D.default_F = default_F
    if oversample > 0:
      self.D.oversample = oversample
    if self.using_cuda:
      self.D.device_Id = device_Id
    self._full_roi = self.D.region_of_interest

  def generate_simulated_image(self, instantiate=False):
    if instantiate:
      self.instantiate_nanoBragg()
    self._add_nanoBragg_spots()
    self.D.raw_pixels /= len(self.beam.xray_beams)
    self._add_background()
    if self.include_noise:
      self._add_noise()
    return self.D.raw_pixels.as_numpy_array()

  def _add_nanoBragg_spots(self):
    rois = self.rois
    if rois is None:
      rois = [self._full_roi]
    _rawpix = None # cuda_add_spots doesnt add spots, it resets each time.. hence we need this
    for roi in rois:
      self.D.region_of_interest = roi
      if self.using_cuda:
        self.D.add_nanoBragg_spots_cuda()
        if _rawpix is None and len(rois) > 1:
          _rawpix = deepcopy(self.D.raw_pixels)
        elif _rawpix is not None:
          _rawpix += self.D.raw_pixels

      elif self.using_omp:
        from boost_adaptbx.boost.python import streambuf  # will deposit printout into dummy StringIO as side effect
        from six.moves import StringIO
        self.D.add_nanoBragg_spots_nks(streambuf(StringIO()))
      else:
        self.D.add_nanoBragg_spots()

    if self.using_cuda and _rawpix is not None:
      self.D.raw_pixels = _rawpix

  def _add_background(self):
    if self.background_raw_pixels is not None:
      self.D.raw_pixels += self.background_raw_pixels
    else:
      if self.add_water:
        print('add water %f mm' % self.water_path_mm)
        water_scatter = flex.vec2_double([
          (0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.18, 7.32), (0.2, 6.75),
          (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
        self.D.Fbg_vs_stol = water_scatter
        self.D.amorphous_sample_thick_mm = self.water_path_mm
        self.D.amorphous_density_gcm3 = 1
        self.D.amorphous_molecular_weight_Da = 18
        self.D.add_background(1, 0)
      if self.add_air:
        print('add air %f mm' % self.air_path_mm)
        air_scatter = flex.vec2_double([(0, 14.1), (0.045, 13.5), (0.174, 8.35), (0.35, 4.78), (0.5, 4.22)])
        self.D.Fbg_vs_stol = air_scatter
        self.D.amorphous_sample_thick_mm = self.air_path_mm
        self.D.amorphous_density_gcm3 = 0.0012
        self.D.amorphous_sample_molecular_weight_Da = 28
        self.D.add_background(1, 0)

  def _add_noise(self):
    self.D.detector_psf_kernel_radius_pixels = 5
    self.D.detector_psf_type = shapetype.Unknown
    self.D.detector_psf_fwhm_mm = self.psf_fwhm
    self.D.readout_noise = self.readout_noise
    self.D.quantum_gain = self.gain
    self.D.add_noise()

  @staticmethod
  def simple_detector(detector_distance_mm, pixelsize_mm, image_shape, fast=(1, 0, 0), slow=(0, -1, 0)):
    from dxtbx.model.detector import DetectorFactory
    import numpy as np
    trusted_range = (0, 200000000000000.0)
    detsize_s = image_shape[0] * pixelsize_mm
    detsize_f = image_shape[1] * pixelsize_mm
    cent_s = (detsize_s + pixelsize_mm * 2) / 2.0
    cent_f = (detsize_f + pixelsize_mm * 2) / 2.0
    beam_axis = np.cross(fast, slow)
    origin = -np.array(fast) * cent_f - np.array(slow) * cent_s + beam_axis * detector_distance_mm
    return DetectorFactory.make_detector('', fast, slow, origin, (
      pixelsize_mm, pixelsize_mm), image_shape, trusted_range)


if __name__ == '__main__':
  S = SimData()
  img = S.generate_simulated_image(instantiate=True)
  print('Maximum pixel value: %.3g' % img.max())
  print('Minimum pixel value: %.3g' % img.min())
