from __future__ import absolute_import, division, print_function
import json, h5py, numpy as np
from scipy.stats import binned_statistic
from scipy import signal
from scipy import constants
from scipy.signal import argrelmax, argrelmin, savgol_filter
import time
from scitbx.array_family import flex
from simtbx.nanoBragg import nanoBragg
from simtbx.nanoBragg.nanoBragg_beam import NBbeam
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from scitbx.matrix import sqr

ENERGY_CONV = 10000000000.0 * constants.c * constants.h / constants.electron_volt


def ensure_p1(Crystal, Famp):
  high_symm_symbol = Famp.space_group_info().type().lookup_symbol()
  cb_op = Famp.space_group_info().change_of_basis_op_to_primitive_setting()
  dtrm = sqr(cb_op.c().r().as_double()).determinant()
  if not dtrm == 1:
    Crystal = Crystal.change_basis(cb_op)
    Famp = Famp.change_basis(cb_op)
  return Crystal, Famp


def flexBeam_sim_colors(CRYSTAL, DETECTOR, BEAM, Famp, energies, fluxes,
                        pids=None, cuda=False, oversample=0, Ncells_abc=(50, 50, 50),
                        mos_dom=1, mos_spread=0, beamsize_mm=0.001, device_Id=0, omp=False,
                        show_params=True, crystal_size_mm=0.01, printout_pix=None, time_panels=True,
                        verbose=0, default_F=0, interpolate=0, recenter=True, profile="gauss",
                        spot_scale_override=None, background_raw_pixels=None, include_noise=False,
                        add_water = False, add_air=False, water_path_mm=0.005, air_path_mm=0, rois_perpanel=None,
                        adc_offset=0, readout_noise=3, psf_fwhm=0, gain=1, mosaicity_random_seeds=None, nopolar=False):
  """
  :param CRYSTAL: dxtbx Crystal model
  :param DETECTOR: dxtbx detector model
  :param BEAM: dxtbx beam model
  :param Famp: cctbx miller array (amplitudes)
  :param energies: list of energies to simulate the scattering
  :param fluxes:  list of pulse fluences per energy (same length as energies)
  :param pids: panel ids to simulate on (None means all panels)
  :param cuda: whether to use GPU (only works for nvidia builds)
  :param oversample: pixel oversample factor (0 means nanoBragg will decide)
  :param Ncells_abc: number of unit cells along each crystal direction in the mosaic block
  :param mos_dom: number of mosaic domains in used to sample mosaic spread (texture)
  :param mos_spread: mosaicity in degrees (spherical cap width)
  :param beamsize_mm: focal size of the beam
  :param device_Id: cuda device id (ignore if cuda=False)
  :param omp: whether to use open mp (required open MP build configuration)
  :param show_params: show the nanoBragg parameters
  :param crystal_size_mm: size of the crystal (increases the intensity of the spots)
  :param printout_pix: debug pixel position : tuple of (pixel_fast_coord, pixel_slow_coord)
  :param time_panels: show timing info
  :param verbose: verbosity level for nanoBragg (0-10), 0 is quiet
  :param default_F: default amplitude value for nanoBragg
  :param interpolate: whether to interpolate for small mosaic domains
  :param recenter: recenter for tilted cameras, deprecated
  :param profile: profile shape, can be : gauss, round, square, or tophat
  :param spot_scale_override: scale the simulated scattering bythis amounth (overrides value based on crystal thickness)
  :param background_raw_pixels: dictionary of {panel_id: raw_pixels}, add these background pixels to the simulated Bragg
  :param include_noise: add noise to simulated pattern
  :param add_water: add water to similated pattern
  :param add_air: add ait to simulated pattern
  :param water_path_mm: length of water the beam travels through
  :param air_path_mm: length of air the beam travels through
  :param rois_perpanel: regions of intererest on each panel
  :param adc_offset: add this value to each pixel in simulated pattern
  :param readout_noise: readout noise level (usually 3-5 ADU)
  :param psf_fwhm: point spread kernel FWHM
  :param gain: photon gain
  :param mosaicity_random_seeds: random seeds to simulating mosaic texture
  :param nopolar: switch of polarization
  :return: list of [(panel_id0,simulated pattern0), (panel_id1, simulated_pattern1), ...]
  """

  CRYSTAL, Famp = ensure_p1(CRYSTAL, Famp)

  if pids is None:
    pids = range(len(DETECTOR))

  if background_raw_pixels is None:
    background_raw_pixels = {pid: None for pid in pids}

  if rois_perpanel is None:
    rois_perpanel = {pid: None for pid in pids}

  nbBeam = NBbeam()
  nbBeam.size_mm = beamsize_mm
  nbBeam.unit_s0 = BEAM.get_unit_s0()
  wavelengths = ENERGY_CONV / np.array(energies)
  nbBeam.spectrum = list(zip(wavelengths, fluxes))

  nbCrystal = NBcrystal()
  nbCrystal.dxtbx_crystal = CRYSTAL
  nbCrystal.miller_array = Famp
  nbCrystal.Ncells_abc = Ncells_abc
  nbCrystal.symbol = CRYSTAL.get_space_group().info().type().lookup_symbol()
  nbCrystal.thick_mm = crystal_size_mm
  nbCrystal.xtal_shape = profile
  nbCrystal.n_mos_domains = mos_dom
  nbCrystal.mos_spread_deg = mos_spread

  panel_images = []

  tinit = time.time()
  S = SimData()
  S.detector = DETECTOR
  S.beam = nbBeam
  S.crystal = nbCrystal
  S.using_cuda = cuda
  S.using_omp = omp
  S.add_air = add_air
  S.air_path_mm = air_path_mm
  S.add_water = add_water
  S.water_path_mm = water_path_mm
  S.readout_noise = readout_noise
  S.gain = gain
  S.psf_fwhm = psf_fwhm
  S.include_noise = include_noise

  if mosaicity_random_seeds is not None:
    S.mosaic_seeds = mosaicity_random_seeds

  S.instantiate_nanoBragg(verbose=verbose, oversample=oversample, interpolate=interpolate, device_Id=device_Id,
                            default_F=default_F, adc_offset=adc_offset)

  if printout_pix is not None:
    S.update_nanoBragg_instance("printout_pixel_fastslow", printout_pix)
  if spot_scale_override is not None:
    S.update_nanoBragg_instance("spot_scale", spot_scale_override)
  S.update_nanoBragg_instance("nopolar", nopolar)

  if show_params:
    S.D.show_params()

  for pid in pids:
    t_panel = time.time()
    S.background_raw_pixels = background_raw_pixels[pid]
    S.panel_id = pid
    S.rois = rois_perpanel[pid]

    S.generate_simulated_image()

    if show_params:
      S.D.show_params()
      print('spot scale: %2.7g' % S.D.spot_scale)
    panel_image = S.D.raw_pixels.as_numpy_array()
    panel_images.append([pid, panel_image])
    S.D.raw_pixels*=0
    if time_panels:
      tdone = time.time() - tinit
      t_panel = time.time() - t_panel
      print('Panel %d took %.4f seconds (Total sim time = %.4f seconds)' % (pid,t_panel, tdone))

  S.D.free_all()

  return panel_images


def sim_background(DETECTOR, BEAM, wavelengths, wavelength_weights, total_flux, pidx=0, beam_size_mm=0.001,
                   Fbg_vs_stol=None, sample_thick_mm=100, density_gcm3=1, molecular_weight=18):
  """
  :param DETECTOR:
  :param BEAM: see sim_spots
  :param wavelengths: see sim_spots
  :param wavelength_weights: see sim_spots
  :param total_flux: see sim_spots
  :param pidx: see sim_spots
  :param beam_size_mm: see sim_spots
  :param Fbg_vs_stol: list of tuples where each tuple is (Fbg, sin theta over lambda)
  :param sample_thick_mm: path length of background that is exposed by the beam
  :param density_gcm3: density of background  (defaults to water)
  :param molecular_weight: molecular weight of background (defaults to water)
  :return: raw_pixels as flex array, these can be passed to sim_spots function below
  """
  wavelength_weights = np.array(wavelength_weights)
  weights = wavelength_weights / wavelength_weights.sum() * total_flux
  spectrum = list(zip(wavelengths, weights))
  xray_beams = get_xray_beams(spectrum, BEAM)
  SIM = nanoBragg(DETECTOR, BEAM, panel_id=(int(pidx)))
  SIM.beamsize_mm = beam_size_mm
  SIM.xray_beams = xray_beams
  if Fbg_vs_stol is None:
    Fbg_vs_stol = flex.vec2_double([
      (0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.18, 7.32), (0.2, 6.75),
      (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
  SIM.flux = sum(weights)
  SIM.Fbg_vs_stol = Fbg_vs_stol
  SIM.amorphous_sample_thick_mm = sample_thick_mm
  SIM.amorphous_density_gcm3 = density_gcm3
  SIM.amorphous_molecular_weight_Da = molecular_weight
  SIM.progress_meter = False
  SIM.add_background()
  background_raw_pixels = SIM.raw_pixels.deep_copy()
  SIM.free_all()
  del SIM
  return background_raw_pixels


def get_xray_beams(spectrum, beam):
  """
  :param spectrum:  list of tuples where one tuple is (wavelength_Angstrom, flux)
  :param beam: beam where we derive the s0 vector and polarization and divergence
  :return: flex_Beam array to be set as a nanoBragg property
  """
  nb = NBbeam()
  nb.spectrum = spectrum
  nb.polarization_fraction = beam.get_polarization_fraction()
  nb.divergence = beam.get_divergence()
  nb.unit_s0 = beam.get_unit_s0()
  return nb.xray_beams


class H5AttributeGeomWriter:

  def __init__(self, filename, image_shape, num_images, detector, beam, dtype=None, compression_args=None, detector_and_beam_are_dicts=False):
    """
    Simple class for writing dxtbx compatible HDF5 files

    :param filename:  input file path
    :param image_shape: shape of a single image (Npanel x Nfast x Nslow)
    :param num_images: how many images will you be writing to the file
    :param detector: dxtbx detector model
    :param beam: dxtbx beam model
    :param dtype: datatype for storage
    :param compression_args: compression arguments for h5py, lzf is performant and simple
        if you only plan to read file in python
        Examples:
          compression_args={"compression": "lzf"}  # Python only
          comression_args = {"compression": "gzip", "compression_opts":9}
    :param detector_and_beam_are_dicts:
    """
    if compression_args is None:
      compression_args = {}
    self.file_handle = h5py.File(filename, 'w')
    self.beam = beam
    self.detector = detector
    self.detector_and_beam_are_dicts = detector_and_beam_are_dicts
    if dtype is None:
      dtype = np.float64
    dset_shape = (
                   num_images,) + tuple(image_shape)
    self.image_dset = (self.file_handle.create_dataset)('images', shape=dset_shape, dtype=dtype, **compression_args)
    self._write_geom()
    self._counter = 0

  def add_image(self, image):
    """
    :param image: a single image as numpy image, same shape as used to instantiate the class
    """
    if self._counter >= self.image_dset.shape[0]:
      raise IndexError('Maximum number of images is %d' % self.image_dset.shape[0])
    self.image_dset[self._counter] = image
    self._counter += 1

  def _write_geom(self):
    beam = self.beam
    det = self.detector
    if not self.detector_and_beam_are_dicts:
      beam = beam.to_dict()
      det = det.to_dict()
    self.image_dset.attrs['dxtbx_beam_string'] = json.dumps(beam)
    self.image_dset.attrs['dxtbx_detector_string'] = json.dumps(det)

  def __exit__(self, exc_type, exc_val, exc_tb):
    self.file_handle.close()

  def __enter__(self):
    return self

  def close_file(self):
    """
    close the file handle (if instantiated using `with`, then this is done automatically)
    """
    self.file_handle.close()


def fcalc_from_pdb(resolution, algorithm=None, wavelength=0.9, anom=True, ucell=None, symbol=None, as_amplitudes=True):
  pdb_lines = 'HEADER TEST\nCRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1\nATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O\nATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O\nATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O\nATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O\nATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O\nATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE\nEND\n'
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None, lines=pdb_lines)
  xray_structure = pdb_inp.xray_structure_simple()
  if ucell is not None:
    assert symbol is not None
    from cctbx.xray import structure
    from cctbx import crystal
    crystal_sym = crystal.symmetry(unit_cell=ucell, space_group_symbol=symbol)
    xray_structure = structure(scatterers=(xray_structure.scatterers()), crystal_symmetry=crystal_sym)
  scatterers = xray_structure.scatterers()
  if anom:
    from cctbx.eltbx import henke
    for sc in scatterers:
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()

  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  fcalc = P1_primitive_xray_structure.structure_factors(d_min=resolution,
                                                        anomalous_flag=anom,
                                                        algorithm=algorithm).f_calc()
  if as_amplitudes:
    fcalc = fcalc.amplitudes().set_observation_type_xray_amplitude()
  return fcalc


def downsample_spectrum(energies, fluences, total_flux=1e12, nbins=100, method=0, ev_width=1.5, baseline_sigma=3.5,
                        method2_param=None):
  """
  :param energies:
  :param fluences:
  :param total_flux:
  :param nbins:
  :param method:
  :param ev_width:
  :param baseline_sigma:
  :return:
  """
  if method2_param is None:
    method2_param = {"filt_freq": 0.07, "filt_order": 3, "tail": 50, "delta_en": 1}
  if method == 0:
    energy_bins = np.linspace(energies.min() - 1e-6, energies.max() + 1e-6, nbins + 1)
    fluences = np.histogram(energies, bins=energy_bins, weights=fluences)[0]
    energies = .5 * (energy_bins[:-1] + energy_bins[1:])

    # only simulate if significantly above the baselein (TODO make more accurate)
    cutoff = np.median(fluences) * 0.8
    is_finite = fluences > cutoff
    fluences = fluences[is_finite]
    energies = energies[is_finite]
  elif method==1:
    w = fluences
    med = np.median(np.hstack((w[:100] ,w[-100:])))
    sigma = np.std(np.hstack((w[:100] ,w[-100:])))
    baseline = med + baseline_sigma*sigma
    width = ev_width/((energies[-1] - energies[0]) / len(energies))
    idx_min = argrelmin(savgol_filter(w,21, 11),order=int(width/3.))[0]
    idx_max = argrelmax(savgol_filter(w,21, 11),order=int(width/3.))[0]
    idx = sorted(np.hstack((idx_min, idx_max)))
    kept_idx = [i for i in idx if w[i] > baseline]
    energies = energies[kept_idx]
    fluences = fluences[kept_idx]
  elif method==2:
    delta_en = method2_param["delta_en"]
    tail = method2_param["tail"]
    filt_order = method2_param["filt_order"]
    filt_freq = method2_param["filt_freq"]
    xdata = np.hstack((energies[:tail], energies[-tail:]))
    ydata = np.hstack((fluences[:tail], fluences[-tail:]))
    pfit = np.polyfit(xdata, ydata, deg=1)
    baseline = np.polyval(pfit, energies)
    denz = signal.filtfilt(*signal.butter(filt_order, filt_freq,'low'), fluences-baseline)
    enbin = np.arange(energies.min(), energies.max(), delta_en)
    fluences, _, _ = binned_statistic(energies, denz, bins=enbin)
    energies = (enbin[1:] + enbin[:-1])*0.5
    fluences[fluences < 0] = 0  # probably ok

  fluences /= fluences.sum()
  fluences *= total_flux
  return energies, fluences
