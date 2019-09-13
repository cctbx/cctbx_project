
import scitbx
from scitbx.array_family import flex
import numpy as np
from simtbx.diffBragg.sim_data2 import SimData
from IPython import embed
import pylab as plt


class RefineLS49Rot:

    def __init__(self):

        # get information from the simulated data image
        from simtbx.diffBragg.load_ls49 import process_ls49_image
        data = process_ls49_image()
        self.dxcry = data["dxcrystal"]
        self.dxdet = data["dxdetector"]
        self.dxbeam = data["dxbeam"]
        self.spot_hkl = data["mill_idx"]
        self.spot_rois = data["bboxes_x1x2y1y2"]
        self.abc_init = data["tilt_abc"]
        self.spectrum = data["spectrum"]
        self.Fhkl = data["sfall"]
        self.img = data["data_img"]

        self._make_nanoBragg_crystal()
        self._make_nanoBragg_beam()
        # NOTE spot_rois are x1,x2,y1,y2 where x=fast, y=slow

        self.n_spots = len(self.spot_rois)
        self.n_background_params = 3*self.n_spots
        self.n_rot_params = 3
        self.n_gain_params = 1
        self.n = self.n_background_params + self.n_rot_params + self.n_gain_params
        self.x = flex.double(self.n)
        self.x[-1] = 1

        self._cache_roi_arrays()
        self._move_abc_init_to_x()
        self._set_diffBragg_instance()

        #for roi in self.nanoBragg_rois:
        #    self.D.region_of_interest = roi
        #    #self.D.init_raw_pixels_roi()
        #    self.D.add_diffBragg_spots()
        #    ((x1,x2),(y1,y2)) = roi
        #    #img0 = self.D.get_raw_pixels_roi().as_numpy_array()
        #    img1 = self.D.raw_pixels.as_numpy_array()[y1:y2, x1:x2]
        #    embed()
        #    #assert np.allclose(img0, img1 , atol=1e-4)
        #exit()

        self.terminator = scitbx.lbfgs.termination_parameters(
            traditional_convergence_test=False,
            traditional_convergence_test_eps=0.05, #1.e-3 significantly (4x) quicker than 1.e-4
            drop_convergence_test_max_drop_eps=1.e-5, # using drop convergence test (since traditional = False), drastic effect
            #min_iterations=min_iterations,
            #max_iterations = None,
            max_calls=100)
        self.minimizer = scitbx.lbfgs.run(
            target_evaluator=self,
            termination_params=self.terminator)

    def _make_nanoBragg_beam(self):
        from simtbx.diffBragg import nanoBragg_beam
        self.nbbeam = nanoBragg_beam.nanoBragg_beam()
        self.nbbeam.spectrum = self.spectrum
        self.nbbeam.size_mm = 0.003
        self.nbbeam.unit_s0 = self.dxbeam.get_unit_s0()

    def _make_nanoBragg_crystal(self):
        from simtbx.diffBragg import nanoBragg_crystal
        self.nbcry = nanoBragg_crystal.nanoBragg_crystal()
        self.nbcry.dxtbx_crystal = self.dxcry
        self.nbcry.mos_spread_deg = 0.05
        self.nbcry.n_mos_domains = 25
        self.nbcry.Ncells_abc = 72, 72, 72
        self.nbcry.miller_array = self.Fhkl

    def _move_abc_init_to_x(self):
        for i in range(self.n_spots):
            self.x[i] = self.abc_init[i, 0]
            self.x[self.n_spots+i] = self.abc_init[i, 1]
            self.x[2*self.n_spots+i] = self.abc_init[i, 2]

    def _cache_roi_arrays(self):
        self.nanoBragg_rois = []  # special nanoBragg format
        self.xrel, self.yrel, self.roi_img = [], [], []
        for x1, x2, y1, y2 in self.spot_rois:

            self.nanoBragg_rois.append(((x1, x2), (y1, y2)))
            yr, xr = np.indices((y2-y1+1, x2-x1+1))
            self.xrel.append(xr)
            self.yrel.append(yr)
            self.roi_img.append(self.img[y1:y2+1, x1:x2+1])

    def _set_diffBragg_instance(self):
        self.S = SimData()
        self.S.using_diffBragg_spots = True
        self.S.crystal = self.nbcry
        self.S.beam = self.nbbeam
        self.S.detector = self.dxdet
        self.S.seed = 1
        self.S.instantiate_diffBragg(adc_offset=10,
                                     default_F=0,
                                     verbose=0,
                                     oversample=1)
        self.D = self.S.D
        self.D.spot_scale = 93.75
        self.D.refine(0)
        self.D.refine(1)
        self.D.refine(2)
        self.D.initialize_managers()

    def _run_diffBragg_current(self, i_spot):
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.set_value(0, self.thetaX)
        self.D.set_value(1, self.thetaY)
        self.D.set_value(2, self.thetaZ)
        self.D.add_diffBragg_spots()

    def _set_background_plane(self, i_spot):
        xr = self.xrel[i_spot]
        yr = self.yrel[i_spot]
        self.tilt_plane = xr*self.a + yr*self.b + self.c

    def _extract_pixel_data(self, i_spot):
        #x1, x2, y1, y2 = self.spot_rois[i_spot]
        self.dRotX = self.D.get_derivative_pixels(0)
        self.dRotX = self.dRotX.as_numpy_array()

        self.dRotY = self.D.get_derivative_pixels(1)
        self.dRotY = self.dRotY.as_numpy_array()

        self.dRotZ = self.D.get_derivative_pixels(2)
        self.dRotZ = self.dRotZ.as_numpy_array()

        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _evaluate_averageI(self):
        """model_Lambda means expected intensity in the pixel"""
        self.model_Lambda = self.tilt_plane + self.scale_fac * self.model_bragg_spots

    def _unpack_params(self, i_spot):
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots*2 + i_spot]
        self.thetaX = self.x[self.n_spots*3]
        self.thetaY = self.x[self.n_spots*3+1]
        self.thetaZ = self.x[self.n_spots*3+2]
        self.scale_fac = self.x[-1]

    def _evaluate_log_averageI(self):
        # fix log(x<=0)
        self.log_Lambda = np.log(self.model_Lambda)
        self.log_Lambda[self.model_Lambda <= 0] = 0

    def compute_functional_and_gradients(self):
        self._set_diffBragg_instance()
        f = 0
        g = flex.double(len(self.x))
        for i_spot in range(self.n_spots):

            self._unpack_params(i_spot)
            self._run_diffBragg_current(i_spot)
            self._set_background_plane(i_spot)
            self._extract_pixel_data(i_spot)
            self._evaluate_averageI()
            self._evaluate_log_averageI()

            Imeas = self.roi_img[i_spot]

            f += (self.model_Lambda - Imeas*self.log_Lambda).sum()

            one_minus_k_over_Lambda = (1. - Imeas / self.model_Lambda)

            # compute gradients for background plane constants a,b,c
            xr = self.xrel[i_spot]  # fast scan pixels
            yr = self.yrel[i_spot]  # slow scan pixels
            g[i_spot] += (xr * one_minus_k_over_Lambda).sum()  # from handwritten notes
            g[self.n_spots + i_spot] += (yr * one_minus_k_over_Lambda).sum()
            g[self.n_spots*2 + i_spot] += one_minus_k_over_Lambda.sum()
            plt.cla()
            plt.subplot(121)
            plt.imshow(self.model_Lambda)
            plt.subplot(122)
            plt.imshow(Imeas)
            plt.suptitle("Spot %d / %d" % (i_spot+1, self.n_spots))
            plt.draw()
            plt.pause(0.7)

            # rotation derivative
            g[self.n_spots*3] += (one_minus_k_over_Lambda * (self.dRotX)).sum()
            g[self.n_spots*3+1] += (one_minus_k_over_Lambda * (self.dRotY)).sum()
            g[self.n_spots*3+2] += (one_minus_k_over_Lambda * (self.dRotZ)).sum()

            # scale factor derivative
            g[-1] += ((self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

        plt.cla()
        plt.title("f=%g"%f)
        plt.imshow(self.D.raw_pixels.as_numpy_array(), vmax=200)
        plt.draw()
        plt.pause(2.2)

        #self.D.raw_pixels *= 0
        #self.D.initialize_managers()
        self.print_step("LBFGS stp", f)
        return f, g

    def print_step(self, message, target):
        print ("%s %10.4f" % (message, target),
               "[", " ".join(["%9.6f" % a for a in self.x]), "]")


if __name__=="__main__":
    RefineLS49Rot()
