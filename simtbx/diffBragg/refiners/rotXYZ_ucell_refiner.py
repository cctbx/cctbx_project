from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
import pylab as plt


class RefineMissetAndUcell(RefineRot):

    def __init__(self, ucell_manager, rotXYZ_refine=(True, True, True), init_scale=1, *args, **kwargs):
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

    def _setup(self):
        # total number of refinement parameters
        n_bg = 3 * self.n_spots
        n_spotscale = 1
        n_params = n_bg + self.n_rot_param + self.n_ucell_param + n_spotscale
        self.x = flex.double(n_params)

        self.rotX_xpos = n_bg
        self.rotY_xpos = n_bg + 1
        self.rotZ_xpos = n_bg + 2

        # populate the x-array with initial values
        self._move_abc_init_to_x()
        self.x[self.rotX_xpos] = 0
        self.x[self.rotY_xpos] = 0
        self.x[self.rotZ_xpos] = 0
        for i in range(self.n_ucell_param):
            self.x[3 * self.n_spots + self.n_rot_param + i] = self.ucell_manager.variables[i]
        self.x[-1] = self._init_scale  # initial scale factor

        # setup the diffBragg instance
        self.D = self.S.D
        if self.refine_rotX:
            self.D.refine(0)  # rotX
        if self.refine_rotY:
            self.D.refine(1)  # rotY
        if self.refine_rotZ:
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
                i + 3,
                self.ucell_manager.derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        if self.refine_rotX:
            self.D.set_value(0, self.x[self.rotX_xpos])
        if self.refine_rotY:
            self.D.set_value(1, self.x[self.rotY_xpos])
        if self.refine_rotZ:
            self.D.set_value(2, self.x[self.rotZ_xpos])
        self.D.Bmatrix = self.ucell_manager.B_recipspace
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        if self.refine_rotX:
            self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
        if self.refine_rotY:
            self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
        if self.refine_rotZ:
            self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()

        self.ucell_deriv = []
        for i in range(self.n_ucell_param):
            self.ucell_deriv.append(self.D.get_derivative_pixels(3 + i).as_numpy_array())
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots * 2 + i_spot]

    def _update_ucell(self):
        _s = slice(self.n_spots * 3 + self.n_rot_param, self.n_spots * 3 + self.n_rot_param + self.n_ucell_param, 1)
        self.ucell_manager.variables = list(self.x[_s])
        self._send_gradients_to_derivative_managers()

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(len(self.x))
        self.scale_fac = self.x[-1]
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
            one_minus_k_over_Lambda = (1. - Imeas / self.model_Lambda)

            # compute gradients for background plane constants a,b,c
            xr = self.xrel[i_spot]  # fast scan pixels
            yr = self.yrel[i_spot]  # slow scan pixels
            g[i_spot] += (xr * one_minus_k_over_Lambda).sum()  # from handwritten notes
            g[self.n_spots + i_spot] += (yr * one_minus_k_over_Lambda).sum()
            g[self.n_spots * 2 + i_spot] += one_minus_k_over_Lambda.sum()
            if self.plot_images:
                plt.cla()
                plt.subplot(121)
                im = plt.imshow(self.model_Lambda)
                plt.subplot(122)
                im2 = plt.imshow(Imeas)
                im.set_clim(im2.get_clim())
                plt.suptitle("Spot %d / %d" % (i_spot + 1, self.n_spots))
                plt.draw()
                plt.pause(.02)

            # might not be necessary to have if statements # TODO : try taking out the if
            if self.refine_rotX:
                g[self.rotX_xpos] += (one_minus_k_over_Lambda * (self.rot_deriv[0])).sum()
            if self.refine_rotY:
                g[self.rotY_xpos] += (one_minus_k_over_Lambda * (self.rot_deriv[1])).sum()
            if self.refine_rotZ:
                g[self.rotZ_xpos] += (one_minus_k_over_Lambda * (self.rot_deriv[2])).sum()

            # unit cell derivative
            for i_ucell_p in range(self.n_ucell_param):
                g[self.n_spots * 3 + self.n_rot_param + i_ucell_p] += (
                            one_minus_k_over_Lambda * (self.ucell_deriv[i_ucell_p])).sum()

            # scale factor derivative
            g[-1] += (self.model_bragg_spots * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
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
        print ("Ucell: %s *** Missets: %s" %
               (", ".join(ucell_labels),
                ", ".join(rot_labels)))

    def get_refined_Bmatrix(self):
        return self.ucell_manager.B_recipspace
