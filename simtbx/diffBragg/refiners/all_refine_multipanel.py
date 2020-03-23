from __future__ import print_function
from simtbx.diffBragg.refiners import RefineRot, BreakToUseCurvatures
from scitbx.array_family import flex
import pylab as plt
import numpy as np
import sys
from IPython import embed


class RefineAllMultiPanel(RefineRot):

    def __init__(self, ucell_manager, rotXYZ_refine=(True, True, True), init_gain=1, init_scale=1, *args, **kwargs):
        """
        :param ucell_manager:
        :param rotXYZ_refine:
        :param args:
        :param kwargs:
        """

        RefineRot.__init__(self, *args, **kwargs)
        self.calc_func = True
        self.multi_panel = True
        self.f_vals = []

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
        self.best_image = np.zeros_like(self.img)
        self.num_positive_curvatures = 0
        self._panel_id = None

        self.ave_ucell = [78.95, 38.12]
        self.sig_ucell = [0.025, 0.025]
        self.sig_rot = 0.01

    def _setup(self):
        # total number of refinement parameters
        n_bg = 0  #3 * self.n_spots
        n_spotscale = 2  # a crystal scsale and an overall scale (e.g. gain)
        n_origin_params = 1  # Currently just the Z-component of the origin is being refined..
        n_ncells_params = 1

        # make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {i: pid for i, pid in enumerate(np.unique(self.panel_ids))}
        self.idx_from_pid = {pid: i for i, pid in enumerate(np.unique(self.panel_ids))}
        n_panels = len(self.pid_from_idx)

        self.n = n_bg + self.n_rot_param + self.n_ucell_param + n_ncells_params + n_origin_params*n_panels + n_spotscale
        self.x = flex.double(self.n)

        self.rotX_xpos = n_bg
        self.rotY_xpos = n_bg + 1
        self.rotZ_xpos = n_bg + 2
        self.x[self.rotX_xpos] = 0
        self.x[self.rotY_xpos] = 0
        self.x[self.rotZ_xpos] = 0

        self._move_abc_init_to_x()

        self.ucell_xstart = self.rotZ_xpos + 1
        # populate the x-array with initial values
        for i_uc in range(self.n_ucell_param):
            self.x[self.ucell_xstart + i_uc] = self.ucell_manager.variables[i_uc]

        # put in Ncells abc estimate
        self.ncells_xpos = self.ucell_xstart + self.n_ucell_param
        self.x[self.ncells_xpos] = self.S.crystal.Ncells_abc[0]

        # put in estimates for origin vectors
        self.origin_xstart = self.ncells_xpos + 1
        for i_pan in range(n_panels):
            pid = self.pid_from_idx[i_pan]
            self.x[self.origin_xstart + i_pan] = self.S.detector[pid].get_local_origin()[2]

        # lastly, the scale factors
        self.x[-2] = self._init_gain  # initial gain for experiment
        self.x[-1] = self._init_scale  # initial scale factor

        # setup the diffBragg instance
        self.D = self.S.D

        if self.refine_Umatrix:
            if self.refine_rotX:
                self.D.refine(0)  # rotX
            if self.refine_rotY:
                self.D.refine(1)  # rotY
            if self.refine_rotZ:
                self.D.refine(2)  # rotZ
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.D.refine(i + 3)  # unit cell params

        if self.refine_ncells:
            self.D.refine(self._ncells_id)
        if self.refine_detdist:
            self.D.refine(self._originZ_id)
        self.D.initialize_managers()
        self.x_init = self.x.as_numpy_array()

    @property
    def x(self):
        """LBFGS parameter array"""
        return self._x

    @x.setter
    def x(self, val):
        self._x = val

    #@profile
    def _send_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i + 3,
                self.ucell_manager.derivative_matrices[i])
            if self.calc_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i + 3, self.ucell_manager.second_derivative_matrices[i])

    #@profile
    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        print("a")
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        print("b")
        if self.iterations==18:
            from IPython import embed
            embed()
        self.D.add_diffBragg_spots()
        print("c")

    #@profile
    def _update_rotXYZ(self):
        if self.refine_rotX:
            self.D.set_value(0, self.x[self.rotX_xpos])
        if self.refine_rotY:
            self.D.set_value(1, self.x[self.rotY_xpos])
        if self.refine_rotZ:
            self.D.set_value(2, self.x[self.rotZ_xpos])

    #@profile
    def _update_ncells(self):
        self.D.set_value(self._ncells_id, self.x[self.ncells_xpos])

    #@profile
    def _update_dxtbx_detector(self):

        det = self.S.detector
        self.S.panel_id = self._panel_id
        # TODO: select hierarchy level at this point
        # NOTE: what does fast-axis and slow-axis mean
        # for the different hierarchy levels?
        node = det[self._panel_id]
        orig = node.get_local_origin()
        xpos = self.origin_xstart + self.idx_from_pid[self._panel_id]
        new_originZ = self.x[xpos]
        new_origin = orig[0], orig[1], new_originZ
        node.set_local_frame(node.get_local_fast_axis(),
                             node.get_local_slow_axis(),
                             new_origin)
        self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)

    #@profile
    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        self.rot_second_deriv = [0, 0, 0]
        if self.refine_Umatrix:
            if self.refine_rotX:
                self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[0] = self.D.get_second_derivative_pixels(0).as_numpy_array()
            if self.refine_rotY:
                self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[1] = self.D.get_second_derivative_pixels(1).as_numpy_array()
            if self.refine_rotZ:
                self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[2] = self.D.get_second_derivative_pixels(2).as_numpy_array()

        self.ucell_derivatives = [0]*self.n_ucell_param
        self.ucell_second_derivatives = [0]*self.n_ucell_param
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.ucell_derivatives[i] = self.D.get_derivative_pixels(3 + i).as_numpy_array()
                if self.calc_curvatures:
                    self.ucell_second_derivatives[i] = self.D.get_second_derivative_pixels(3 + i).as_numpy_array()

        self.ncells_deriv = self.detdist_deriv = 0
        self.ncells_second_deriv = self.detdist_second_deriv = 0
        if self.refine_ncells:
            self.ncells_deriv = self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()
            if self.calc_curvatures:
                self.ncells_second_deriv = self.D.get_second_derivative_pixels(self._ncells_id).as_numpy_array()
        if self.refine_detdist:
            self.detdist_deriv = self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()
            if self.calc_curvatures:
                self.detdist_second_deriv = self.D.get_second_derivative_pixels(self._originZ_id).as_numpy_array()

        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    #@profile
    def _update_best_image(self, i_spot):
        pid = self.panel_ids[i_spot]
        x1, x2, y1, y2 = self.spot_rois[i_spot]
        self.best_image[pid, y1:y2 + 1, x1:x2 + 1] = self.model_Lambda

    #@profile
    def _unpack_bgplane_params(self, i_spot):
        self.a, self.b, self.c = self.abc_init[i_spot]
        #self.a = self.x[i_spot]
        #self.b = self.x[self.n_spots + i_spot]
        #self.c = self.x[self.n_spots * 2 + i_spot]

    #@profile
    def _update_ucell(self):
        _s = slice(self.ucell_xstart, self.ucell_xstart + self.n_ucell_param, 1)
        pars = list(self.x[_s])
        self.ucell_manager.variables = pars
        self._send_gradients_to_derivative_managers()
        self.D.Bmatrix = self.ucell_manager.B_recipspace

    #@profile
    def compute_functional_and_gradients(self):
        if self.calc_func:
            if self.verbose:
                if self.use_curvatures:
                    print("Compute functional and gradients Iter %d (Using Curvatures)\n<><><><><><><><><><><><><>"
                          % (self.iterations+1))
                else:
                    print("Compute functional and gradients Iter %d PosCurva %d\n<><><><><><><><><><><><><>"
                          % (self.iterations+1, self.num_positive_curvatures))
            f = 0
            g = flex.double(self.n)
            if self.calc_curvatures:
                self.curv = flex.double(self.n)
            self.best_image *= 0
            self.gain_fac = self.x[-2]
            self.scale_fac = self.x[-1]
            G2 = self.gain_fac**2
            S2 = self.scale_fac**2
            self._update_ucell()
            self._update_ncells()
            self._update_rotXYZ()
            for i_spot in range(self.n_spots):
                self._panel_id = self.panel_ids[i_spot]
                if self.verbose:
                    print( "\rRunning diffBragg over spot %d/%d on panel %d " % (i_spot+1, self.n_spots, self._panel_id), flush=True)

                self._update_dxtbx_detector()  # TODO : chek that I wrk
                self._run_diffBragg_current(i_spot)
                self._unpack_bgplane_params(i_spot)
                self._set_background_plane(i_spot)
                self._extract_pixel_data()
                self._evaluate_averageI()
                if self.poisson_only:
                    self._evaluate_log_averageI()
                else:
                    self._evaluate_log_averageI_plus_sigma_readout()
                self._update_best_image(i_spot)
                self.Imeas = self.roi_img[i_spot]

                # helper terms for doing derivatives
                one_over_Lambda = 1. / self.model_Lambda
                self.one_minus_k_over_Lambda = (1. - self.Imeas * one_over_Lambda)
                self.k_over_squared_Lambda = self.Imeas * one_over_Lambda * one_over_Lambda

                self.u = self.Imeas - self.model_Lambda
                self.one_over_v = 1./(self.model_Lambda + self.sigma_r**2)
                self.one_minus_2u_minus_u_squared_over_v = 1-2*self.u-self.u*self.u*self.one_over_v

                f += self._target_accumulate()

                # compute gradients for background plane constants a,b,c
                xr = self.xrel[i_spot]  # fast scan pixels
                yr = self.yrel[i_spot]  # slow scan pixels

                if self.plot_images and self.iterations % self.plot_stride == 0:
                    if self.plot_residuals:
                        self.ax.clear()
                        residual = self.model_Lambda - self.Imeas
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
                        m = self.Imeas[self.Imeas > 1e-9].mean()
                        s = self.Imeas[self.Imeas > 1e-9].std()
                        vmax = m+5*s
                        vmin = m-s
                        m2 = self.model_Lambda.mean()
                        s2 = self.model_Lambda.std()
                        self.ax1.images[0].set_data(self.model_Lambda)
                        self.ax1.images[0].set_clim(vmin, vmax)
                        self.ax2.images[0].set_data(self.Imeas)
                        self.ax2.images[0].set_clim(vmin, vmax)
                    plt.suptitle("Iterations = %d, image %d / %d"
                                 % (self.iterations, i_spot+1, self.n_spots))
                    self.fig.canvas.draw()
                    plt.pause(0.5)

                if self.refine_Umatrix:
                    for ii, xpos in enumerate([self.rotX_xpos, self.rotY_xpos, self.rotZ_xpos]):
                        d = S2*G2*self.rot_deriv[ii]
                        g[xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = S2*G2*self.rot_second_deriv[ii]
                            self.curv[xpos] += self._curv_accumulate(d, d2)

                if self.refine_Bmatrix:
                    # unit cell derivative
                    for i_ucell_p in range(self.n_ucell_param):
                        xpos = self.ucell_xstart + i_ucell_p
                        d = S2*G2*self.ucell_derivatives[i_ucell_p]
                        g[xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = S2*G2*self.ucell_second_derivatives[i_ucell_p]
                            self.curv[xpos] += self._curv_accumulate(d, d2)

                if self.refine_ncells:
                    d = S2*G2*self.ncells_deriv
                    g[self.ncells_xpos] += self._grad_accumulate(d)
                    if self.calc_curvatures:
                        d2 = S2*G2*self.ncells_second_deriv
                        self.curv[self.ncells_xpos] += self._curv_accumulate(d, d2)

                if self.refine_detdist:
                    raise NotImplementedError("Cannot refined detdist (yet...)")
                    #if self.calc_curvatures:
                    #    raise NotImplementedError("Cannot use curvatures and refine detdist (yet...)")
                    #origin_xpos = self.origin_xstart + self.idx_from_pid[self._panel_id]
                    #g[origin_xpos] += (S2*G2*self.detdist_deriv*one_minus_k_over_Lambda).sum()

                if self.refine_gain_fac:
                    d = 2*self.gain_fac*(self.tilt_plane + S2*self.model_bragg_spots)
                    g[-2] += self._grad_accumulate(d)
                    if self.calc_curvatures:
                        d2 = d / self.gain_fac
                        self.curv[-2] += self._curv_accumulate(d, d2)

                if self.refine_crystal_scale:
                    d = G2*2*self.scale_fac*self.model_bragg_spots
                    g[-1] += self._grad_accumulate(d)
                    if self.calc_curvatures:
                        d2 = d / self.scale_fac
                        self.curv[-1] += self._curv_accumulate(d, d2)
            if self.calc_curvatures and not self.use_curvatures:
                if np.all(self.curv.as_numpy_array() >= 0):
                    self.num_positive_curvatures += 1
                else:
                    self.num_positive_curvatures = 0

            # TODO add in the priors:
            if self.use_ucell_priors and self.refine_Bmatrix:
                for jj in range(self.n_ucell_param):
                    xpos = self.ucell_xstart + jj
                    ucell_p = self.x[xpos]
                    sig_square = self.sig_ucell[jj]**2
                    f += (ucell_p - self.ave_ucell[jj])**2 / 2 / sig_square
                    g[xpos] += (ucell_p-self.ave_ucell[jj])/sig_square
                    if self.calc_curvatures:
                        self.curv[xpos] += 1/sig_square

            if self.use_rot_priors and self.refine_Umatrix:
                x_positions = [self.rotX_xpos,
                               self.rotY_xpos,
                               self.rotZ_xpos]
                for xpos in x_positions:
                    rot_p = self.x[xpos]
                    sig_square = self.sig_rot**2
                    f += rot_p**2 / 2 / sig_square
                    g[xpos] += rot_p / sig_square
                    if self.calc_curvatures:
                        self.curv[xpos] += 1/sig_square

            self._f = f
            self._g = g
            self.D.raw_pixels *= 0
            gnorm = np.linalg.norm(g)
            if self.verbose:
                self.print_step("LBFGS stp", f)
                self.print_step_grads("LBFGS GRADS", gnorm )
            self.iterations += 1
            self.f_vals.append(f)

            if self.calc_curvatures and not self.use_curvatures:
                if self.num_positive_curvatures == self.use_curvatures_threshold:
                    raise BreakToUseCurvatures

        return self._f, self._g

    def _poisson_target(self):
        fterm = (self.model_Lambda - self.Imeas * self.log_Lambda).sum()
        return fterm

    def _poisson_d(self, d):
        gterm = (d*self.one_minus_k_over_Lambda).sum()
        return gterm

    def _poisson_d2(self, d, d2):
        cterm = d2*self.one_minus_k_over_Lambda + d*d*self.k_over_squared_Lambda
        return cterm.sum()

    def _gaussian_target(self):
        fterm = .5 * (self.log2pi + self.log_Lambda_plus_sigma_readout + self.u * self.u * self.one_over_v).sum()
        return fterm

    def _gaussian_d(self, d):
        gterm = .5*(d*self.one_over_v*self.one_minus_2u_minus_u_squared_over_v).sum()
        return gterm

    def _gaussian_d2(self, d, d2):
        cterm = self.one_over_v*(d2*self.one_minus_2u_minus_u_squared_over_v -
                                 d*d*(self.one_over_v*self.one_minus_2u_minus_u_squared_over_v -
                                      (2+2*self.u*self.one_over_v + self.u*self.u*self.one_over_v*self.one_over_v)))
        cterm = .5*cterm.sum()
        return cterm

    def _evaluate_log_averageI_plus_sigma_readout(self):
        L = self.model_Lambda + self.sigma_r**2
        is_neg = (L <= 0).ravel()
        if np.any(is_neg):
            print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda_plus_sigma_readout = np.log(L)
        self.log_Lambda_plus_sigma_readout[L <= 0] = 0

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
              % (message, target, self.x[self.ncells_xpos], self.x[self.origin_xstart], self.x[-2]**2, self.x[-1]**2))
        print ("Ucell: %s *** Missets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def print_step_grads(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for i,(n, v) in enumerate(zip(names, vals)):
            grad = self._g[self.ucell_xstart+i]
            ucell_labels.append('G%s=%+2.7g' % (n, grad))
        rotX = self._g[self.rotX_xpos]
        rotY = self._g[self.rotY_xpos]
        rotZ = self._g[self.rotZ_xpos]
        rot_labels = ["GrotX=%+3.7g" % rotX, "GrotY=%+3.7g" % rotY, "GrotZ=%+3.4g" % rotZ]
        xnorm = np.linalg.norm(self.x)
        print("%s: |x|=%f, |g|=%f, Gncells=%f, Gdetdist=%3.8g, Ggain=%3.4g, Gscale=%4.8g"
              % (message, xnorm, target, self._g[self.ncells_xpos], self._g[self.origin_xstart], self._g[-2], self._g[-1]))
        print ("GUcell: %s *** GMissets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def get_refined_Bmatrix(self):
        return self.ucell_manager.B_recipspace

    def curvatures(self):
        return self.curv
