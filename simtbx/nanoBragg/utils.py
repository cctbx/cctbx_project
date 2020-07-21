from __future__ import absolute_import, division, print_function
import json, h5py, numpy as np
from scipy import constants
import time
from scitbx.array_family import flex
from dxtbx_model_ext import flex_Beam
from dxtbx.model import BeamFactory
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from simtbx.nanoBragg.nanoBragg_beam import NBbeam
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData, Amatrix_dials2nanoBragg

ENERGY_CONV = 10000000000.0 * constants.c * constants.h / constants.electron_volt


def flexBeam_sim_colors(CRYSTAL, DETECTOR, BEAM, Famp, energies, fluxes, pids=None, cuda=False, oversample=0, Ncells_abc=(50, 50, 50), mos_dom=1, mos_spread=0, exposure_s=1, beamsize_mm=0.001, device_Id=0, show_params=False, crystal_size_mm=0.01, printout_pix=None, time_panels=True, verbose=0, default_F=0, interpolate=0, seed=3142019, recenter=True):
    nbBeam = NBbeam()
    nbBeam.size_mm = beamsize_mm
    nbBeam.unit_s0 = BEAM.get_unit_s0()
    wavelengths = ENERGY_CONV / np.array(energies)
    nbBeam.spectrum = list(zip(wavelengths, fluxes))
    nbCrystal = NBcrystal()
    nbCrystal.miller_array = Famp
    nbCrystal.thick_mm = crystal_size_mm
    mos_spread_matrices = SimData.Umats(mos_spread, mos_dom)
    if pids is None:
        pids = range(len(DETECTOR))
    panel_images = []
    for pid in pids:
        tinit = time.time()
        D = nanoBragg(DETECTOR, BEAM, verbose=verbose,
          panel_id=(int(pid)))
        if recenter:
            D.beam_center_mm = DETECTOR[int(pid)].get_beam_centre(BEAM.get_s0())
        if printout_pix is not None:
            D.printout_pixel_fastslow = printout_pix
        D.seed = seed
        D.calib_seed = seed
        D.mosaic_seed = seed
        D.exposure_s = exposure_s
        D.interpolate = interpolate
        D.Fhkl_tuple = (
         nbCrystal.miller_array.indices(), nbCrystal.miller_array.data())
        D.Amatrix = Amatrix_dials2nanoBragg(CRYSTAL)
        D.Ncells_abc = Ncells_abc
        D.mosaic_spread_deg = mos_spread
        D.mosaic_domains = mos_dom
        D.set_mosaic_blocks(mos_spread_matrices)
        D.xray_beams = nbBeam.xray_beams
        D.beamsize_mm = nbBeam.size_mm
        D.flux = sum(fluxes)
        mosaic_vol = D.xtal_size_mm[0] * D.xtal_size_mm[1] * D.xtal_size_mm[2]
        D.spot_scale = determine_spot_scale(nbBeam, nbCrystal, mosaic_vol)
        D.default_F = default_F
        if cuda:
            D.device_Id = int(device_Id)
        D.xtal_shape = nbCrystal.xtal_shape
        if oversample > 0:
            D.oversample = oversample
        if cuda:
            D.add_nanoBragg_spots_cuda()
        else:
            D.add_nanoBragg_spots()
        if show_params:
            D.show_params()
            print('Mosaic domain volume: %2.7g (mm^3)' % mosaic_vol)
            print('spot scale: %2.7g' % D.spot_scale)
        panel_image = D.raw_pixels.as_numpy_array()
        if cuda:
            panel_image /= len(energies)
        panel_images.append(panel_image)
        D.free_all()
        if time_panels:
            tdone = time.time() - tinit
            print('Panel %d took %.4f seconds' % (pid, tdone))
        del D

    return np.array(panel_images)


def sim_background(DETECTOR, BEAM, wavelengths, wavelength_weights, total_flux, pidx=0, beam_size_mm=0.001, Fbg_vs_stol=None, sample_thick_mm=100, density_gcm3=1, molecular_weight=18):
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
         (0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.2, 6.75), (0.18, 7.32),
         (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
    SIM.flux = total_flux
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


def sim_spots(CRYSTAL, DETECTOR, BEAM, Famp, wavelengths, wavelength_weights, total_flux, pidx=0, cuda=False, oversample=0, mosaic_vol_A3=27000000000, mos_dom=1, mos_spread=0, profile=None, crystal_size_mm=0.01, beam_size_mm=0.001, device_Id=0, show_params=False, printout_pix=None, time_panels=True, verbose=0, default_F=0, interpolate=0, noise_seed=None, calib_seed=None, mosaic_seed=None, recenter=True, spot_scale_override=None, add_noise=True, adc_offset=10, readout_noise_adu=3, gain=1, background_raw_pixels=None, background_scale=1):
    """
    :param CRYSTAL:  dxtbx crystal model (e.g. from dials.index or dials.stills_process)
    :param DETECTOR: dxtbx detector model
    :param BEAM: dxtbx beam model
    :param Famp: cctbx miller array (.data() attribute should point to SF amplitudes)
    :param wavelengths: wavelengths in Angstrom
    :param wavelength_weights: weights for each wavelength in Angstroms (for pink beam simulation)
    :param total_flux: total flux per pulse
    :param pidx: panel index
    :param cuda: whether to use cuda, if so provide correct device Id (use nvidia-smi to check)
    :param oversample: oversample the pixel by this factor. If 0, it will be auto determined by nanoBragg
    :param mosaic_vol_A3: volume of mosaic domain in crystal
    :param mos_dom: for mosacic spread calculations, use this many mosaic domains to sample the spread
    :param mos_spread: width of mosaic spread distribution (degrees) for simuting spherical cap effect
    :param profile: profile of RELP (can be "round", "tophat", "square", or "gauss"). default is "square"
    :param crystal_size_mm: diameter of crystal (in mm)
    :param beam_size_mm: diameter of beam focus (in mm)
    :param device_Id: for running GPU, provide the device Id (usualy 0 for hosts with 1 GPU)
    :param show_params: display the nanoBragg prams at the start
    :param printout_pix: fast scan/ slow scan coordinate of a pixel to printout low level nanoBragg parameters (useful for debug)
    :param time_panels: time the simulation
    :param verbose: verbosity level for nanobragg (0=silent 10=blarg)
    :param default_F: default amplitude for miller indices not included in Famp
    :param interpolate: use interpolation for small mosaic domains (set to 0 to disable)
    :param noise_seed: seed for generating random poisson and readout noise
    :param calib_seed: seed for generating gain calibration noise
    :param mosaic_seed: seed for generating mosaic spread
    :param recenter: recenter the panel after instantiating (do this if the detector and incident beam are not squred up, e.g. for tilted geometries)
    :param spot_scale_override: override the spot scale parameter which is currently determined by taking ration of crystal volume and mosaic domain volume
    :param add_noise: whether to add noise
    :param adc_offset: offset to pixel values
    :param readout_noise_adu: readout noise level (default is 3)
    :param gain: quantum gain (default is 1)
    :param background_raw_pixels: flex array of background pixels (output from sim_background function)
    :param background_scale: option to boost ot decrease the background level
    :return: simulated pixels as a numpy array, that can then be written to an hdf5 file
    """
    if not len(wavelengths) == len(wavelength_weights):
        raise AssertionError
    else:
        wavelength_weights = np.array(wavelength_weights)
        weights = wavelength_weights / wavelength_weights.sum() * total_flux
        spectrum = list(zip(wavelengths, weights))
        xray_beams = get_xray_beams(spectrum, BEAM)
        tinit = time.time()
        SIM = nanoBragg(DETECTOR, BEAM, verbose=verbose,
          panel_id=(int(pidx)))
        if profile is not None:
            if profile == 'gauss':
                SIM.xtal_shape = shapetype.Gauss
            else:
                if profile == 'tophat':
                    SIM.xtal_shape = shapetype.Tophat
                else:
                    if profile == 'round':
                        SIM.xtal_shape = shapetype.Round
                    elif profile == 'square':
                        SIM.xtal_shape = shapetype.Square
        if recenter:
            SIM.beam_center_mm = DETECTOR[int(pidx)].get_beam_centre(BEAM.get_s0())
        if printout_pix is not None:
            SIM.printout_pixel_fastslow = printout_pix
        if noise_seed is not None:
            SIM.seed = noise_seed
        if calib_seed is not None:
            SIM.calib_seed = calib_seed
        if mosaic_seed is not None:
            SIM.mosaic_seed = mosaic_seed
        SIM.exposure_s = 1
        SIM.interpolate = interpolate
        Nunit_cell = mosaic_vol_A3 / CRYSTAL.get_unit_cell().volume()
        N = np.power(Nunit_cell, 0.3333333333333333)
        SIM.Ncells_abc = (int(N), int(N), int(N))
        SIM.mosaic_spread_deg = mos_spread
        SIM.mosaic_domains = mos_dom
        if isinstance(Famp, float):
            SIM.default_F = Famp
        else:
            SIM.Fhkl = Famp
        SIM.Amatrix = Amatrix_dials2nanoBragg(CRYSTAL)
        SIM.spot_scale = determine_spot_scale(beam_size_mm, crystal_size_mm, mosaic_vol_A3)
        if spot_scale_override is not None:
            SIM.spot_scale = spot_scale_override
        SIM.flux = total_flux
        SIM.beamsize_mm = beam_size_mm
        SIM.xray_beams = xray_beams
        SIM.default_F = default_F
        if cuda:
            SIM.device_Id = int(device_Id)
        if oversample > 0:
            SIM.oversample = oversample
        if cuda:
            SIM.add_nanoBragg_spots_cuda()
        else:
            SIM.add_nanoBragg_spots()
        if show_params:
            SIM.show_params()
            print('Mosaic domain volume: %2.7g (mm^3)' % mosaic_vol_A3)
            print('spot scale: %2.7g' % SIM.spot_scale)
        SIM.raw_pixels /= len(wavelengths)
        if background_raw_pixels is not None:
            SIM.raw_pixels += background_raw_pixels * background_scale
        if add_noise:
            SIM.adc_offset_adu = adc_offset
            SIM.detector_psf_fwhm_mm = 0
            SIM.quantum_gain = gain
            SIM.readout_noise_adu = readout_noise_adu
            SIM.add_noise()
        raw_pixels = SIM.raw_pixels.as_numpy_array()
        SIM.free_all()
        if time_panels:
            tdone = time.time() - tinit
            print('Panel %d took %.4f seconds' % (pidx, tdone))
    del SIM
    return raw_pixels


def determine_spot_scale(beam_size_mm, crystal_thick_mm, mosaic_vol_A3):
    """

    :param beam_size_mm:  diameter of beam focus (millimeter)
    :param crystal_thick_mm: thickness of crystal (millimeter)
    :param mosaic_vol_A3:  volume of a mosaic block in crystal (cubic angstrom)
    :return: roughly the number of exposed mosaic blocks
    """
    if beam_size_mm <= crystal_thick_mm:
        illum_xtal_vol = crystal_thick_mm * beam_size_mm ** 2
    else:
        illum_xtal_vol = crystal_thick_mm ** 3
    return illum_xtal_vol / mosaic_vol_A3 * 1e+21


def get_xray_beams(spectrum, beam_originator):
    """

    :param spectrum:  list of tuples where one tuple is (wavelength_Angstrom, flux)
    :param beam_originator: beam where we derive the s0 vector and polarization and divergence
    :return: flex_Beam array to be set as a nanoBragg property
    """
    xray_beams = flex_Beam()
    for wavelen, flux in spectrum:
        beam = BeamFactory.simple(wavelen * 1e-10)
        beam.set_flux(flux)
        beam.set_unit_s0(beam_originator.get_unit_s0())
        beam.set_polarization_fraction(beam_originator.get_polarization_fraction())
        beam.set_divergence(beam_originator.get_divergence())
        xray_beams.append(beam)

    return xray_beams


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
# okay decompiling utils.cpython-36.pyc
