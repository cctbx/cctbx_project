
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.beam import BeamFactory
from scitbx.array_family import flex
from simtbx.nanoBragg.tst_nanoBragg_basic import fcalc_from_pdb
from simtbx.nanoBragg import shapetype
from simtbx.diffBragg import diffBragg


class SimData:
    def __init__(self):
        pass

    @property
    def seed(self):
        return 1

    @property
    def detector_distance_mm(self):
        return 180

    @property
    def pixelsize_mm(self):
        return 0.1

    @property
    def wavelength_A(self):
        return 1.24

    @property
    def image_shape(self):
        return (1024,1024)

    @property
    def detector(self):

        det_descr = {'panels':
                         [{'fast_axis': (-1.0, 0.0, 0.0),
                           'gain': 1.0,
                           'identifier': '',
                           'image_size':self.image_shape,
                           'mask': [],
                           'material': '',
                           'mu': 0.0,
                           'name': 'Panel',
                           'origin': (self.image_shape[0]/2*self.pixelsize_mm, -self.image_shape[0]/2*self.pixelsize_mm,
                                      -self.detector_distance_mm),
                           'pedestal': 0.0,
                           'pixel_size': (self.pixelsize_mm, self.pixelsize_mm),
                           'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                           'raw_image_offset': (0, 0),
                           'slow_axis': (0.0, 1.0, 0.0),
                           'thickness': 0.0,
                           'trusted_range': (0.0, 2**14),
                           'type': ''}]}

        return DetectorFactory.from_dict(det_descr)

    @property
    def beam(self):

        beam_descr = {'direction': (0.0, 0.0, 1.0),
                      'divergence': 0.0,
                      'flux': 1,
                      'polarization_fraction': 1.,
                      'polarization_normal': (0.0, 1.0, 0.0),
                      'sigma_divergence': 0.0,
                      'transmission': 1.0,
                      'wavelength': self.wavelength_A}

        return BeamFactory.from_dict(beam_descr)

    @property
    def flux(self):
        return 2e12

    @property
    def Fhkl(self):
        return fcalc_from_pdb(resolution=4, algorithm="fft", wavelength=self.beam.get_wavelength())

    @property
    def mosaic_domains(self):
        return 10

    @property
    def mosaic_spread_deg(self):
        return 0.01

    @property
    def Ncells_abc(self):
        return (15,15,15)

    @property
    def xtal_shape(self):
        return shapetype.Gauss

    @property
    def beamsize_mm(self):
        return 5e-3

    @property
    def spot_scale(self):
        return 1e6

    @property
    def progress_meter(self):
        return False

    @property
    def verbose(self):
        return 0

    def instantiate_diffBragg(self):

        self.D = diffBragg(self.detector, self.beam, verbose=0, panel_id=0)
        self.D.xtal_shape = self.xtal_shape
        self.D.Ncells_abc = self.Ncells_abc
        self.D.wavelength_A = self.beam.get_wavelength()
        self.D.flux = self.flux
        self.D.mosaic_spread_deg = self.mosaic_spread_deg
        self.D.mosaic_domains = self.mosaic_domains
        self.D.Fhkl = self.Fhkl
        self.D.spot_scale = self.spot_scale
        self.D.beamsize_mm = self.beamsize_mm
        self.D.progress_meter = self.progress_meter
        self.D.verbose = self.verbose
        self.D.mosaic_seed = self.seed
        self.D.seed = self.seed
        self.D.calib_seed = self.seed
        self.D.vectorize_umats()

    def generate_simulated_image(self):
        self.instantiate_diffBragg()
        print "add spots"
        self._add_diffBragg_spots()
        print "add background"
        self._add_background()
        print("add noise")
        self._add_noise()
        print("Done!")
        return self.D.raw_pixels.as_numpy_array()

    def _add_diffBragg_spots(self):
        self.D.add_diffBragg_spots()

    def _add_background(self):
        water_scatter = flex.vec2_double(
            [(0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.2, 6.75), (0.18, 7.32),
             (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
        self.D.Fbg_vs_stol = water_scatter
        self.D.amorphous_sample_thick_mm = 0.005  # typical Gself.DVN jet thickness but could vary
        self.D.amorphous_density_gcm3 = 1
        self.D.amorphous_molecular_weight_Da = 18
        self.D.add_background()

    def _add_noise(self):
        self.D.adc_offset_adu = 10
        self.D.detector_psf_kernel_radius_pixels = 5
        self.D.detector_psf_type = shapetype.Unknown  # for CSPAself.D
        self.D.detector_psf_fwhm_mm = 0
        self.D.quantum_gain = 1
        self.D.add_noise()


if __name__=="__main__":
    S = SimData()
    img = S.generate_simulated_image()
    print "Maximum pixel value: %.3g" % img.max()
    print "Minimum pixel value: %.3g" % img.min()
