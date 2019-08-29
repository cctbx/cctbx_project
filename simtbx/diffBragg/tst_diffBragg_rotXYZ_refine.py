
import scitbx
from scitbx.array_family import flex
import numpy as np
from simtbx.diffBragg.process_simdata import process_simdata
from simtbx.diffBragg.sim_data_for_tests import SimData
from IPython import embed


class RefineRot:

    def __init__(self):

        # get information from the simulated data image
        self.spot_rois, self.spot_hkl, self.Amat_init, self.Amat_known, self.abc_init, self.img = \
            process_simdata(plot=False)

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
        embed()

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

    def _move_abc_init_to_x(self):
        for i in range(self.n_spots):
            self.x[i] = self.abc_init[i,0]
            self.x[self.n_spots+i] = self.abc_init[i,1]
            self.x[2*self.n_spots+i] = self.abc_init[i,2]

    def _cache_roi_arrays(self):
        self.nanoBragg_rois = []  # special nanoBragg format
        self.xrel, self.yrel, self.roi_img = [], [],[]
        for x1,x2,y1,y2 in self.spot_rois:

            self.nanoBragg_rois.append( ((x1,x2),(y1,y2)))
            yr, xr = np.indices( (y2-y1, x2-x1))
            self.xrel.append(xr)
            self.yrel.append(yr)
            self.roi_img.append(self.img[y1:y2, x1:x2])

    def _set_diffBragg_instance(self):
        self.S = SimData()
        self.S.instantiate_diffBragg()
        self.D = self.S.D
        self.D.Amatrix = self.Amat_init.inverse()  # modify the Amatrix!
        #self.D.vectorize_umats()

    def _run_diffBragg_current(self, i_spot):

        xr = self.xrel[i_spot]
        yr = self.yrel[i_spot]
        tilt_plane = xr*self.a + yr*self.b + self.c
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.set_value(0, self.thetaX)
        self.D.set_value(1, self.thetaY)
        self.D.set_value(2, self.thetaZ)
        self.D.initialize_managers()  # TODO verify zero-ing the raw_pixel images
        self.D.add_diffBragg_spots()

        self.dRotX = self.D.get_derivative_pixels(0).as_numpy_array()
        self.dRotY = self.D.get_derivative_pixels(1).as_numpy_array()
        self.dRotZ = self.D.get_derivative_pixels(2).as_numpy_array()

        # TODO : convert to numpy  only ROI pixel
        self.model_lambda = tilt_plane + self.D.raw_pixels.as_numpy_array()
        #self.model_lambda[ self.model_lambda < 0] = 1e-9  # TODO does this even happen ?

        self.D.raw_pixels *=0  # TODO, only 0 the ROI

    def _unpack_params(self, i_spot):
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots*2 +i_spot]
        self.thetaX = self.x[self.n_spots*3]
        self.thetaY = self.x[self.n_spots*3+1]
        self.thetaZ = self.x[self.n_spots*3+2]
        self.scale_fac = self.x[self.n_spots*3+3]

    def compute_functional_and_gradients(self):

        f = 0
        g = flex.double(3)
        for i_spot in range( self.n_spots):
            xr = self.xrel[i_spot]
            yr = self.yrel[i_spot]
            self._unpack_params(i_spot)
            self._run_diffBragg_current(i_spot)

            datapts = self.roi_img[i_spot]

            f += (self.model_lambda - np.log(self.model_lambda)*datapts).sum()
            one_minus_k_over_lambda = (1. - datapts / self.model_lambda)

            g[i_spot ] += (xr * one_minus_k_over_lambda).sum()  # from handwritten notes
            g[self.n_spots +i_spot] += (yr * one_minus_k_over_lambda).sum()
            g[self.n_spots*2 + i_spot] += one_minus_k_over_lambda.sum()

            g[self.n_spots*3] += (one_minus_k_over_lambda * self.dRotX).sum()
            g[self.n_spots*3+1] += (one_minus_k_over_lambda * self.dRotY).sum()
            g[self.n_spots*3+2] += (one_minus_k_over_lambda * self.dRotZ).sum()
            g[-1] += 0  #(datapts * one_minus_k_over_lambda).sum()

        self.print_step("LBFGS stp", f)
        return f, g

    def print_step(self, message, target):
        print ("%s %10.4f" % (message, target),
               "[", " ".join(["%9.3f" % a for a in self.x]), "]")

        #f = 0.
        #g = flex.double(self.n)
        #for ispot in range(self.n_spots):
        #    F = self.sb_data[ispot].focus()
        #    for x in range(F[1]):
        #        for y in range(F[2]):
        #            model_lambda = self.a[3 * ispot + 0] * x + self.a[3 * ispot + 1] * y + self.a[3 * ispot + 2] + \
        #                           self.a[-1] * self.roi_model_pixels[ispot][x, y]
        #            datapt = self.sb_data[ispot][0, x, y]  # not sure the right datapt when model_lambda<0
        #            if model_lambda <= 0:
        #                f += model_lambda  # complete kludge, guard against math domain error
        #            else:
        #                f += model_lambda - datapt * math.log(model_lambda)
        #            g[3 * ispot + 0] += x * one_minus_k_over_lambda  # from handwritten notes
        #            g[3 * ispot + 1] += y * one_minus_k_over_lambda
        #            g[3 * ispot + 2] += one_minus_k_over_lambda
        #            # for this paper, never refine rotX as it is parallel to beam and well-determined
        #            g[-4] += one_minus_k_over_lambda * self.a[
        #                -1] * 0.  # always fix rotx.  self.roi_dxyz["Amat_dx"][ispot][x,y]
        #            g[-3] += one_minus_k_over_lambda * self.a[-1] * self.roi_dxyz["Amat_dy"][ispot][x, y]
        #            g[-2] += one_minus_k_over_lambda * self.a[-1] * self.roi_dxyz["Amat_dz"][ispot][x, y]
        #            g[-1] += self.roi_model_pixels[ispot][x, y] * one_minus_k_over_lambda

if __name__=="__main__":
    RefineRot()
