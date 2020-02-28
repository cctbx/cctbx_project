
from scitbx.array_family import flex
from simtbx.nanoBragg import shapetype, nanoBragg
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.nanoBragg_beam import nanoBragg_beam
from simtbx.diffBragg import diffBragg


class SimData:
    def __init__(self):
        self.detector = SimData.simple_detector(180, 0.1, (512, 512))
        self.seed = 1
        self.crystal = nanoBragg_crystal()
        self.add_air = False
        self.add_water = True
        self.water_path_mm = 0.005
        self.air_path_mm = 0
        nbBeam = nanoBragg_beam()
        nbBeam.unit_s0 = (0, 0, -1)
        self.beam = nbBeam
        self.using_cuda = False
        self.panel_id = 0
        self.include_noise = True
        self.using_diffBragg_spots = False
        self.functionals = []

    @property
    def air_path(self):
        return self._air_path

    @air_path.setter
    def air_path(self, val):
        self._air_path = val

    @property
    def water_path(self):
        return self._water_path

    @water_path.setter
    def water_path(self, val):
        self._water_path = val

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
    def Umats(mos_spread_deg, n_mos_doms, isotropic=True,
              seed=777, norm_dist_seed=777):
        import scitbx
        from scitbx.matrix import col
        import math
        UMAT_nm = flex.mat3_double()
        mersenne_twister = flex.mersenne_twister(seed=seed)
        scitbx.random.set_random_seed(norm_dist_seed)
        rand_norm = scitbx.random.normal_distribution(mean=0, sigma=mos_spread_deg * math.pi / 180.)
        g = scitbx.random.variate(rand_norm)
        mosaic_rotation = g(n_mos_doms)
        for m in mosaic_rotation:
            site = col(mersenne_twister.random_double_point_on_sphere())
            if mos_spread_deg > 0:
                UMAT_nm.append(site.axis_and_angle_as_r3_rotation_matrix(m, deg=False))
            else:
                UMAT_nm.append(site.axis_and_angle_as_r3_rotation_matrix(0, deg=False))
            if isotropic and mos_spread_deg > 0:
                UMAT_nm.append(site.axis_and_angle_as_r3_rotation_matrix(-m, deg=False))
        return UMAT_nm

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
    def panel_id(self):
        return self._panel_id

    @panel_id.setter
    def panel_id(self, val):
        self._panel_id = val

    @property
    def detector(self):
        return self._detector

    @detector.setter
    def detector(self, val):
        self._detector = val

    @property
    def using_diffBragg_spots(self):
        return self._using_diffBragg_spots

    @using_diffBragg_spots.setter
    def using_diffBragg_spots(self, val):
        assert(val in [True, False])
        self._using_diffBragg_spots = val

    @property
    def using_cuda(self):
        return self._using_cuda

    @using_cuda.setter
    def using_cuda(self, val):
        assert(val in [True, False])
        self._using_cuda = val

    @property
    def include_noise(self):
        return self._include_noise

    @include_noise.setter
    def include_noise(self, val):
        self._include_noise = val

    def _crystal_properties(self):
        self.D.xtal_shape = self.crystal.xtal_shape

        if self.crystal.miller_array is not None:
            self.D.Fhkl_tuple = self.crystal.miller_array.indices(), self.crystal.miller_array.data()

        ## TODO: am I unnecessary?
        #self.D.unit_cell_tuple = self.crystal.dxtbx_crystal.get_unit_cell().parameters()

        self.D.Omatrix = self.crystal.Omatrix
        self.D.Bmatrix = self.crystal.dxtbx_crystal.get_B() #
        self.D.Umatrix = self.crystal.dxtbx_crystal.get_U()

        self.D.Ncells_abc = self.crystal.Ncells_abc[0]
        self.D.mosaic_spread_deg = self.crystal.mos_spread_deg
        self.D.mosaic_domains = self.crystal.n_mos_domains
        self.D.set_mosaic_blocks(SimData.Umats(
                self.crystal.mos_spread_deg, self.crystal.n_mos_domains))

    def _beam_properties(self):
        self.D.xray_beams = self.beam.xray_beams
        # TODO: make me a circular size ?
        self.D.beamsize_mm = self.beam.size_mm

    def _seedlings(self):
        self.D.seed = self.seed
        self.D.calib_seed = self.seed
        self.D.mosaic_seed = self.seed

    def determine_spot_scale(self):
        if self.beam.size_mm <= self.crystal.thick_mm:
            illum_xtal_vol = self.crystal.thick_mm * self.beam.size_mm**2
        else:
            illum_xtal_vol = self.crystal.thick_mm**3
        mosaic_vol = self.D.xtal_size_mm[0]*self.D.xtal_size_mm[1]*self.D.xtal_size_mm[2]
        return illum_xtal_vol / mosaic_vol

    def update_nanoBragg_instance(self, parameter, value):
        setattr(self.D, parameter, value)

    def instantiate_diffBragg(self, verbose=0, oversample=0, device_Id=0,
                              adc_offset=0, default_F=1e3, interpolate=0):

        self.D = diffBragg(self.detector, self.beam.nanoBragg_constructor_beam,
                           verbose=verbose, panel_id=self.panel_id)
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

        self.D.vectorize_umats()  # NOTE: dont forget me!

    def generate_simulated_image(self, instantiate=False):
        if instantiate:
            self.instantiate_diffBragg()
        print ("add spots")
        self._add_nanoBragg_spots()
        self._add_background()
        if self.include_noise:
            print("add noise")
            self._add_noise()
        print("Done!")
        return self.D.raw_pixels.as_numpy_array()

    def _add_nanoBragg_spots(self):
        if self.using_cuda:
            self.D.add_nanoBragg_spots_cuda()
        else:
            if self.using_diffBragg_spots:
                self.D.add_diffBragg_spots()
            else:
                self.D.add_nanoBragg_spots()

    def _add_background(self):
        if self.add_water:
            print("add water %f mm" % self.water_path_mm)
            water_scatter = flex.vec2_double(
                [(0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.2, 6.75), (0.18, 7.32),
                 (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
            self.D.Fbg_vs_stol = water_scatter
            self.D.amorphous_sample_thick_mm = self.water_path_mm
            self.D.amorphous_density_gcm3 = 1
            self.D.amorphous_molecular_weight_Da = 18
            self.D.add_background(1, 0)

        if self.add_air:
            print("add air %f mm" % self.air_path_mm)
            air_scatter = flex.vec2_double([(0, 14.1), (0.045, 13.5), (0.174, 8.35), (0.35, 4.78), (0.5, 4.22)])
            self.D.Fbg_vs_stol = air_scatter
            self.D.amorphous_sample_thick_mm = self.air_path_mm
            self.D.amorphous_density_gcm3 = 1.2e-3
            self.D.amorphous_sample_molecular_weight_Da = 28  # nitrogen = N2
            self.D.add_background(1, 0)

    def _add_noise(self):
        self.D.detector_psf_kernel_radius_pixels = 5
        self.D.detector_psf_type = shapetype.Unknown  # for CSPAself.D
        self.D.detector_psf_fwhm_mm = 0
        #self.D.quantum_gain = 1 #self.gain
        self.D.add_noise()

    @staticmethod
    def simple_detector(detector_distance_mm, pixelsize_mm, image_shape,
                        fast=(1, 0, 0), slow=(0, -1, 0)):
        from dxtbx.model.detector import DetectorFactory
        import numpy as np
        trusted_range = 0, 2e14
        detsize_s = image_shape[0]*pixelsize_mm
        detsize_f = image_shape[1]*pixelsize_mm
        cent_s = (detsize_s + pixelsize_mm*2)/2.
        cent_f = (detsize_f + pixelsize_mm*2)/2.
        beam_axis = np.cross(fast, slow)
        origin = -np.array(fast)*cent_f - np.array(slow)*cent_s + beam_axis*detector_distance_mm

        return DetectorFactory.make_detector("", fast, slow, origin,
                                             (pixelsize_mm, pixelsize_mm), image_shape, trusted_range)


if __name__ == "__main__":
    S = SimData()
    img = S.generate_simulated_image()
    print "Maximum pixel value: %.3g" % img.max()
    print "Minimum pixel value: %.3g" % img.min()
