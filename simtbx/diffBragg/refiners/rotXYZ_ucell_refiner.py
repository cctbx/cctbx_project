from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
import pylab as plt
import numpy as np


class RefineMissetAndUcell(RefineRot):

    def __init__(self, ucell_manager, rotXYZ_refine=(True, True, True), init_gain=1, init_scale=1, *args, **kwargs):
        """
        :param ucell_manager:
        :param rotXYZ_refine:
        :param args:
        :param kwargs:
        """

        RefineRot.__init__(self, *args, **kwargs)

        self.refine_rotX = rotXYZ_refine[0]
        self.refine_rotY = rotXYZ_refine[1]
        self.refine_rotZ = rotXYZ_refine[2]
        self.n_rot_param = 3  # self.refine_rotX + self.refine_rotY + self.refine_rotZ
        self.ucell_manager = ucell_manager
        self.n_ucell_param = len(self.ucell_manager.variables)
        self._init_scale = init_scale
        self._init_gain = init_gain
        self.best_image = np.zeros_like(self.img)
        self.use_curvatures = False

    def _setup(self):
        # total number of refinement parameters
        if self.refine_background_planes:
            raise NotImplementedError("No planes do I care about:w")
        else:
            n_bg = 0

        n_spotscale = 2
        self.n = n_bg + self.n_rot_param + self.n_ucell_param + n_spotscale
        self.x = flex.double(self.n)

        self.rotX_xpos = n_bg
        self.rotY_xpos = n_bg + 1
        self.rotZ_xpos = n_bg + 2

        # populate the x-array with initial values
        if self.refine_background_planes:
            self._move_abc_init_to_x()
        self.x[self.rotX_xpos] = 0
        self.x[self.rotY_xpos] = 0
        self.x[self.rotZ_xpos] = 0
        self.ucell_xstart = n_bg + self.n_rot_param
        for i in range(self.n_ucell_param):
            self.x[self.ucell_xstart + i] = self.ucell_manager.variables[i]
        self.x[-2] = self._init_gain  # initial gain for experiment
        self.x[-1] = self._init_scale  # initial scale factor

        # setup the diffBragg instance
        self.D = self.S.D
        self.D.refine(0)  # rotX
        self.D.refine(1)  # rotY
        self.D.refine(2)  # rotZ
        for i in range(self.n_ucell_param):
            self.D.refine(i + 3)  # unit cell params
        self.D.initialize_managers()

    @property
    def x(self):
        """LBFGS parameter array"""
        return self._x

    @x.setter
    def x(self, val):
        self._x = val

    def _send_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i+3,
                self.ucell_manager.derivative_matrices[i])
            if self.use_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i+3, self.ucell_manager.second_derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.set_value(0, self.x[self.rotX_xpos])
        self.D.set_value(1, self.x[self.rotY_xpos])
        self.D.set_value(2, self.x[self.rotZ_xpos])
        self.D.Bmatrix = self.ucell_manager.B_recipspace
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        self.rot_second_deriv = [0, 0, 0]
        if self.refine_rotX:
            self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
            if self.use_curvatures:
                self.rot_second_deriv[0] = self.D.get_second_derivative_pixels(0).as_numpy_array()
        if self.refine_rotY:
            self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
            if self.use_curvatures:
                self.rot_second_deriv[1] = self.D.get_second_derivative_pixels(1).as_numpy_array()
        if self.refine_rotZ:
            self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()
            if self.use_curvatures:
                self.rot_second_deriv[2] = self.D.get_second_derivative_pixels(2).as_numpy_array()

        self.ucell_derivatives = []
        self.ucell_second_derivatives = []
        for i in range(self.n_ucell_param):
            self.ucell_derivatives.append(self.D.get_derivative_pixels(3 + i).as_numpy_array())
            if self.use_curvatures:
                self.ucell_second_derivatives.append(self.D.get_second_derivative_pixels(3 + i).as_numpy_array())
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        if self.refine_background_planes:
            raise NotImplementedError("No planish")
            self.a = self.x[i_spot]
            self.b = self.x[self.n_spots + i_spot]
            self.c = self.x[self.n_spots * 2 + i_spot]
        else:
            self.a, self.b, self.c = self.abc_init[i_spot]

    def _update_ucell(self):
        s = slice(self.ucell_xstart, self.ucell_xstart + self.n_ucell_param, 1)
        self.ucell_manager.variables = list(self.x[s])
        self._send_gradients_to_derivative_managers()

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(self.n)
        if self.use_curvatures:
            self.curv = flex.double(self.n)
        self.gain_fac = self.x[-2]
        self.scale_fac = self.x[-1]
        G2 = self.gain_fac**2
        S2 = self.scale_fac**2
        self._update_ucell()
        for i_spot in range(self.n_spots):
            self._run_diffBragg_current(i_spot)
            self._unpack_bgplane_params(i_spot)
            self._set_background_plane(i_spot)
            self._extract_pixel_data()
            self._evaluate_averageI()
            self._evaluate_log_averageI()

            Imeas = self.roi_img[i_spot]
            f += (self.model_Lambda - Imeas * self.log_Lambda).sum()
            one_over_Lambda = 1./self.model_Lambda
            one_minus_k_over_Lambda = (1. - Imeas * one_over_Lambda)
            if self.use_curvatures:
                k_over_squared_Lambda = Imeas * one_over_Lambda * one_over_Lambda

            # compute gradients for background plane constants a,b,c
            if self.refine_background_planes:
                xr = self.xrel[i_spot]  # fast scan pixels
                yr = self.yrel[i_spot]  # slow scan pixels
                raise NotImplementedError("Stopped caring about background planes..")
                g[i_spot] += (xr * G2*one_minus_k_over_Lambda).sum()  # from handwritten notes
                g[self.n_spots + i_spot] += (yr*G2*one_minus_k_over_Lambda).sum()
                g[self.n_spots * 2 + i_spot] += (G2*one_minus_k_over_Lambda).sum()
            if self.plot_images:
                m = Imeas[Imeas > 1e-9].mean()
                s = Imeas[Imeas > 1e-9].std()
                vmax = m+5*s
                vmin = m-s
                m2 = self.model_Lambda.mean()
                s2 = self.model_Lambda.std()
                vmax2 = m2+5*s2
                vmin2 = m2-s2
                self.ax1.images[0].set_data(self.model_Lambda)
                #self.ax1.images[0].set_clim(vmin2, vmax2)
                self.ax1.images[0].set_clim(vmin, vmax)
                #self.ax1.images[0].set_data(self.model_bragg_spots)
                self.ax2.images[0].set_data(Imeas)
                self.ax2.images[0].set_clim(vmin, vmax)
                plt.suptitle("Iterations = %d, image %d / %d"
                             % (self.iterations, i_spot+1, self.n_spots))
                self.fig.canvas.draw()
                plt.pause(.02)

            if self.refine_Amatrix:
                for ii, xpos in enumerate([self.rotX_xpos, self.rotY_xpos, self.rotZ_xpos]):
                    d = self.rot_deriv[ii]
                    g[xpos] += (one_minus_k_over_Lambda * S2*G2*d).sum()
                    if self.use_curvatures:
                        d2 = self.rot_second_deriv[ii]
                        cc = S2*G2*(d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda)
                        #if cc.sum() < 0:
                        #    from IPython import embed
                        #    embed()
                        self.curv[xpos] += cc.sum()

                # unit cell derivative
                for i_ucell_p in range(self.n_ucell_param):
                    xpos = self.ucell_xstart + i_ucell_p
                    d = self.ucell_derivatives[i_ucell_p]
                    g[xpos] += (S2 * G2 * one_minus_k_over_Lambda * d).sum()

                    if self.use_curvatures:
                        d2 = self.ucell_second_derivatives[i_ucell_p]
                        cc = S2 * G2 * (d2 * one_minus_k_over_Lambda + d * d * k_over_squared_Lambda)
                        #if cc.sum() < 0:
                        #    from IPython import embed
                        #    embed()

                        self.curv[xpos] += cc.sum()

            #if self.refine_gain_fac:
            #    g[-2] += (2*self.gain_fac*(self.tilt_plane + S2*self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

            #if self.refine_crystal_scale:
            #    # scale factor derivative
            #    g[-1] += (G2*2*self.scale_fac*self.model_bragg_spots * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        self.iterations += 1
        return f, g

    def print_step(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))

        rotX = self.x[self.rotX_xpos]
        rotY = self.x[self.rotY_xpos]
        rotZ = self.x[self.rotZ_xpos]
        rot_labels = ["rotX=%+2.7g" % rotX, "rotY=%+2.7g" % rotY, "rotZ=%+2.7g" % rotZ]
        print ("Ucell: %s *** Missets: %s ** s=%2.7g, g=%2.7g" %
               (", ".join(ucell_labels),
                ", ".join(rot_labels),
                self.x[-1], self.x[-2]))

    def curvatures(self):
        return self.curv

    def get_refined_Bmatrix(self):
        return self.ucell_manager.B_recipspace
