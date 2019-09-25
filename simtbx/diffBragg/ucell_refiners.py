
from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
import pylab as plt


class RefineUcell(RefineRot):

    def __init__(self, ucell_manager, *args, **kwargs):

        RefineRot.__init__(self, *args, **kwargs)
        self.ucell_manager = ucell_manager
        self.S.crystal.dxtbx_crystal.show()
        self.n_ucell_param = len(self.ucell_manager.variables)

    def _setup(self):
        self._setup_lbfgs_x_array()
        self._cache_roi_arrays()
        self._move_abc_init_to_x()
        self._set_diffBragg_instance()

    def _setup_lbfgs_x_array(self):
        self.n_spots = len(self.spot_rois)
        self.n_background_params = 3*self.n_spots
        self.n_gain_params = 1
        self.n = self.n_background_params + self.n_ucell_param + self.n_gain_params
        self.x = flex.double(self.n)
        for i in range(self.n_ucell_param):
            self.x[3*self.n_spots + i] = self.ucell_manager.variables[i]
        self.x[-1] = 1

    def _set_diffBragg_instance(self):
        self.D = self.S.D
        for i in range(self.n_ucell_param):
            self.D.refine(i+3)
        self.D.initialize_managers()

    def _send_gradients_to_derivative_managers(self):
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i+3,
                self.ucell_manager.derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self._send_gradients_to_derivative_managers()
        self.D.Bmatrix = self.ucell_manager.B_recipspace
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        self.ucell_derivatives = []
        for i in range(self.n_ucell_param):
            self.ucell_derivatives.append(self.D.get_derivative_pixels(3+i).as_numpy_array())
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots*2 + i_spot]

    def _update_orientation(self):
        self.ucell_manager.variables = list(self.x[self.n_spots*3:self.n_spots*3+self.n_ucell_param])
        self._send_gradients_to_derivative_managers()

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(len(self.x))
        self._update_orientation()
        self.scale_fac = self.x[-1]
        for i_spot in range(self.n_spots):
            self._run_diffBragg_current(i_spot)
            self._unpack_bgplane_params(i_spot)
            self._set_background_plane(i_spot)
            self._extract_pixel_data()
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
            if self.plot_images:
                plt.cla()
                plt.subplot(121)
                im = plt.imshow(self.model_Lambda)
                plt.subplot(122)
                im2 = plt.imshow(Imeas)
                im.set_clim(im2.get_clim())
                plt.suptitle("Spot %d / %d" % (i_spot+1, self.n_spots))
                plt.draw()
                plt.pause(.02)

            # unit cell derivative
            for i_ucell_p in range(self.n_ucell_param):
                g[self.n_spots*3+i_ucell_p] += (one_minus_k_over_Lambda * (self.ucell_derivatives[i_ucell_p])).sum()
            #g[self.n_spots*3+1] += (one_minus_k_over_Lambda * (self.ucell_derivatives[1])).sum()
            #g[self.n_spots*3+2] += (one_minus_k_over_Lambda * (self.ucell_derivatives[2])).sum()
            #g[self.n_spots*3+3] += (one_minus_k_over_Lambda * (self.ucell_derivatives[3])).sum()

            # scale factor derivative
            g[-1] += (self.model_bragg_spots * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        return f, g

    #def print_step(self, message, target):
    #    #names = self.ucell_manager.nam
    #    print ("unit cell a= %2.7g, c=%2.7g .. scale = %2.7g" % (self.x[-5], self.x[-4], self.x[-3], self.x[-2], self.x[-1]))
    #    #print ("unit cell a= %2.7g, c=%2.7g .. scale = %2.7g" % (self.x[-5], self.x[-4], self.x[-3], self.x[-2], self.x[-1]))
    #    #print ("unit cell a= %2.7g, c=%2.7g .. scale = %2.7g" % (self.x[-3], self.x[-2], self.x[-1]))
    def print_step(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        labels = []
        for n, v in zip(names, vals):
            #if "beta" in n or "gamma" in n or "alpha" in n:
            #    v = v * 180 / np.pi
            labels.append('%s=%2.7g' % (n, v))
        print ", ".join(labels)
