from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
import pylab as plt
import numpy as np
from dxtbx.model import Panel


class RefineAll(RefineRot):

    def __init__(self, ucell_manager, rotXYZ_refine=(True, True, True), init_gain=1, init_scale=1,
                 panel_id=0, *args, **kwargs):
        """
        :param ucell_manager:
        :param rotXYZ_refine:
        :param args:
        :param kwargs:
        """

        RefineRot.__init__(self, *args, **kwargs)
        self.calc_func = True
        self.f_vals =[]

        self.refine_rotX = rotXYZ_refine[0]
        self.refine_rotY = rotXYZ_refine[1]
        self.refine_rotZ = rotXYZ_refine[2]
        self.n_rot_param = 3  # self.refine_rotX + self.refine_rotY + self.refine_rotZ
        self.ucell_manager = ucell_manager
        self._ncells_id = 9
        self._originZ_id = 10
        self.n_ucell_param = len(self.ucell_manager.variables)
        self._init_scale = init_scale
        self._init_gain = init_gain
        self._panel_id = panel_id
        self.best_image = np.zeros_like(self.img)

    def _setup(self):
        # total number of refinement parameters
        n_bg = 0 #3 * self.n_spots
        n_spotscale = 2
        n_origin_params = 1
        n_ncells_params = 1
        self.n = n_bg + self.n_rot_param + self.n_ucell_param + n_ncells_params + n_origin_params + n_spotscale
        self.n = self.n_rot_param + self.n_ucell_param + n_ncells_params + n_origin_params + n_spotscale
        self.x = flex.double(self.n)

        self.rotX_xpos = n_bg
        self.rotY_xpos = n_bg + 1
        self.rotZ_xpos = n_bg + 2

        self.ucell_xstart = n_bg + self.n_rot_param
        # populate the x-array with initial values
        self._move_abc_init_to_x()
        self.x[self.rotX_xpos] = 0
        self.x[self.rotY_xpos] = 0
        self.x[self.rotZ_xpos] = 0
        for i_uc in range(self.n_ucell_param):
            self.x[self.ucell_xstart + i_uc] = self.ucell_manager.variables[i_uc]
        self.x[-4] = self.S.crystal.Ncells_abc[0]
        self.x[-3] = self.S.detector[self._panel_id].get_origin()[2]
        self.x[-2] = self._init_gain  # initial gain for experiment
        self.x[-1] = self._init_scale  # initial scale factor

        # setup the diffBragg instance
        self.D = self.S.D
        if self.refine_Amatrix:
            if self.refine_rotX:
                self.D.refine(0)  # rotX
            if self.refine_rotY:
                self.D.refine(1)  # rotY
            if self.refine_rotZ:
                self.D.refine(2)  # rotZ
            # TODO if refine_ucell:
            for i in range(self.n_ucell_param):
                self.D.refine(i + 3)  # unit cell params
        if self.refine_ncells:
            self.D.refine(self._ncells_id)
        if self.refine_detdist:
            self.D.refine(self._originZ_id)
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

    def _update_ncells(self):
        self.D.set_value(self._ncells_id, self.x[-4])

    def _update_dxtbx_detector(self):
        det = self.S.detector
        node = det[self._panel_id]
        node_d = node.to_dict()
        new_originZ = self.x[-3]
        node_d["origin"] = node_d["origin"][0], node_d["origin"][1], new_originZ
        det[self._panel_id] = Panel.from_dict(node_d)
        self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)

    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        if self.refine_Amatrix:
            if self.refine_rotX:
                self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
            if self.refine_rotY:
                self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
            if self.refine_rotZ:
                self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()

            self.ucell_deriv = []
            for i in range(self.n_ucell_param):
                self.ucell_deriv.append(self.D.get_derivative_pixels(3 + i).as_numpy_array())

        self.ncells_deriv = self.detdist_deriv = 0
        if self.refine_ncells:
            self.ncells_deriv = self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()
        if self.refine_detdist:
            self.detdist_deriv = self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()

        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _update_best_image(self, i_spot):
        x1, x2, y1, y2 = self.spot_rois[i_spot]
        self.best_image[y1:y2 + 1, x1:x2 + 1] = self.model_Lambda

    def _unpack_bgplane_params(self, i_spot):
        self.a, self.b, self.c = self.abc_init[i_spot]
        #self.a = self.x[i_spot]
        #self.b = self.x[self.n_spots + i_spot]
        #self.c = self.x[self.n_spots * 2 + i_spot]

    def _update_ucell(self):
        _s = slice(self.ucell_xstart, self.ucell_xstart + self.n_ucell_param, 1)
        pars = list(self.x[_s])
        if pars[0] == 0:
            from IPython import embed
            embed()
        self.ucell_manager.variables = list(self.x[_s])
        self._send_gradients_to_derivative_managers()

    def compute_functional_and_gradients(self):

        if self.calc_func:
            f = 0
            g = flex.double(self.n)
            self.best_image = np.zeros_like(self.img)
            self.scale_fac = self.x[-1]
            self.gain_fac = self.x[-2]
            G2 = self.gain_fac**2
            S2 = self.scale_fac**2
            self._update_dxtbx_detector()
            self._update_ucell()
            self._update_ncells()
            for i_spot in range(self.n_spots):
                self._run_diffBragg_current(i_spot)
                self._unpack_bgplane_params(i_spot)
                self._set_background_plane(i_spot)
                self._extract_pixel_data()
                self._evaluate_averageI()
                self._evaluate_log_averageI()
                self._update_best_image(i_spot)

                Imeas = self.roi_img[i_spot]
                f += (self.model_Lambda - Imeas * self.log_Lambda).sum()
                one_minus_k_over_Lambda = (1. - Imeas / self.model_Lambda)

                # compute gradients for background plane constants a,b,c
                xr = self.xrel[i_spot]  # fast scan pixels
                yr = self.yrel[i_spot]  # slow scan pixels
                #if self.refine_background_planes:
                #    g[i_spot] += (xr * G2*one_minus_k_over_Lambda).sum()  # from handwritten notes
                #    g[self.n_spots + i_spot] += (yr*G2*one_minus_k_over_Lambda).sum()
                #    g[self.n_spots * 2 + i_spot] += (G2*one_minus_k_over_Lambda).sum()

                if self.plot_images and self.iterations %5==0:
                    if self.plot_residuals:
                        self.ax.clear()
                        residual = self.model_Lambda - Imeas
                        if i_spot == 0:
                            x = residual.max()
                        else:
                            x = np.mean([x, residual.max()])

                        self.ax.plot_surface(xr, yr, residual, rstride=2, cstride=2, alpha=0.3, cmap='coolwarm')
                        self.ax.contour(xr, yr, residual, zdir='z', offset=-x, cmap='coolwarm')
                        self.ax.set_yticks(range(yr.min(), yr.max()))
                        self.ax.set_xticks(range(xr.min(), xr.max()))
                        self.ax.set_xticklabels([])
                        self.ax.set_yticklabels([])
                        self.ax.set_zlim(-x, x)
                        self.ax.set_title("residual (photons)")
                    else:
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
                    # might not be necessary to have if statements # TODO : try taking out the if
                    if self.refine_rotX:
                        g[self.rotX_xpos] += (one_minus_k_over_Lambda * S2*G2*self.rot_deriv[0]).sum()
                    if self.refine_rotY:
                        g[self.rotY_xpos] += (one_minus_k_over_Lambda * S2*G2*self.rot_deriv[1]).sum()
                    if self.refine_rotZ:
                        g[self.rotZ_xpos] += (one_minus_k_over_Lambda * S2*G2*self.rot_deriv[2]).sum()

                    # unit cell derivative
                    for i_ucell_p in range(self.n_ucell_param):
                        g[self.ucell_xstart + i_ucell_p] += (
                                    one_minus_k_over_Lambda * S2*G2*self.ucell_deriv[i_ucell_p]).sum()

                if self.refine_ncells:
                    # NOTE: we are refining u=log(ncells), hence multiply be deriviative of ncells w.r.t. u, e.g. exp^u
                    g[-4] += (S2*G2*self.ncells_deriv*one_minus_k_over_Lambda).sum()

                if self.refine_detdist:
                    g[-3] += (S2*G2*self.detdist_deriv*one_minus_k_over_Lambda).sum()

                if self.refine_gain_fac:
                    g[-2] += (2*self.gain_fac*(self.tilt_plane + S2*self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

                if self.refine_crystal_scale:
                    # scale factor derivative
                    g[-1] += (G2*2*self.scale_fac*self.model_bragg_spots * one_minus_k_over_Lambda).sum()


            self._f = f
            self._g = g
            self.D.raw_pixels *= 0
            self.print_step("LBFGS stp", f)
            self.print_step_grads("LBFGS GRADS", f)
            self.iterations += 1
            self.f_vals.append(f)

        return self._f, self._g

    def print_step(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))
        rotX = self.x[self.rotX_xpos]
        rotY = self.x[self.rotY_xpos]
        rotZ = self.x[self.rotZ_xpos]
        rot_labels = ["rotX=%+3.7g" % rotX, "rotY=%+3.7g" % rotY, "rotZ=%+3.4g" % rotZ]
        print("%s: residual=%3.8g, ncells=%f, detdist=%3.8g, gain=%3.4g, scale=%4.8g"
              % (message, target, self.x[-4], self.x[-3], self.x[-2]**2, self.x[-1]**2))
        print ("Ucell: %s *** Missets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def print_step_grads(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for i,(n, v) in enumerate(zip(names, vals)):
            grad = self._g[self.ucell_xstart+i]
            ucell_labels.append('%s=%+2.7g' % (n, grad))
        rotX = self._g[self.rotX_xpos]
        rotY = self._g[self.rotY_xpos]
        rotZ = self._g[self.rotZ_xpos]
        rot_labels = ["GrotX=%+3.7g" % rotX, "GrotY=%+3.7g" % rotY, "GrotZ=%+3.4g" % rotZ]
        print("%s: Gresidual=%3.8g, Gncells=%f, Gdetdist=%3.8g, Ggain=%3.4g, Gscale=%4.8g"
              % (message, target, self._g[-4], self._g[-3], self._g[-2], self._g[-1]))
        print ("GUcell: %s *** GMissets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def get_refined_Bmatrix(self):
        return self.ucell_manager.B_recipspace
