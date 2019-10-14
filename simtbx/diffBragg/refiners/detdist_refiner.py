from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
from dxtbx.model import Panel
import pylab as plt
import numpy as np


class RefineDetdist(RefineRot):

    def __init__(self, panel_id=0, init_gain=1, init_scale=1, *args, **kwargs):
        """
        :param ucell_manager:
        :param rotXYZ_refine:
        :param args:
        :param kwargs:
        """
        self.f_vals = []

        RefineRot.__init__(self, *args, **kwargs)
        self._init_scale = init_scale
        self._init_gain = init_gain
        self._panel_id = panel_id
        self._originZ_id = 10  #TODO make this less clunky, this is the derivative manager ID for detector distance
        self.best_image = np.zeros_like(self.img)
        self.refine_detdist = True

    def _setup(self):
        # total number of refinement parameters
        n_bg = 0 #3 * self.n_spots
        n_spotscale = 2
        n_origin = 1  # number of dxtbx detector model origin parameters (1 for det distance)
        self.n = n_bg + n_origin + n_spotscale
        self.x = flex.double(self.n)

        # populate the x-array with initial values
        #self._move_abc_init_to_x()
        from copy import deepcopy
        self.x[-3] = deepcopy(self.S.detector[self._panel_id].get_origin()[2])
        self.x[-2] = self._init_gain  # initial gain for experiment
        self.x[-1] = self._init_scale  # initial scale factor

        # setup the diffBragg instance
        self.D = self.S.D
        self.D.refine(self._originZ_id)
        self.D.initialize_managers()

    def _update_dxtbx_detector(self):
        from copy import deepcopy
        self._mod_det = deepcopy(self.S.detector)
        node = self._mod_det[self._panel_id]
        node_d = node.to_dict()
        new_originZ = self.x[-3]
        node_d["origin"] = node_d["origin"][0], node_d["origin"][1], new_originZ
        self._mod_det[self._panel_id] = Panel.from_dict(node_d)
        #self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(self._mod_det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)

    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        self.detdist_deriv = self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a, self.b, self.c = self.abc_init[i_spot]
        #self.a = self.x[i_spot]
        #self.b = self.x[self.n_spots + i_spot]
        #self.c = self.x[self.n_spots * 2 + i_spot]

    def _update_best_image(self, i_spot):
        x1, x2, y1, y2 = self.spot_rois[i_spot]
        self.best_image[y1:y2 + 1, x1:x2 + 1] = self.model_Lambda

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(self.n)
        self.scale_fac = self.x[-1]
        self.gain_fac = self.x[-2]
        G2 = self.gain_fac**2
        S2 = self.scale_fac**2
        self._update_dxtbx_detector()
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

            if self.plot_images:
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
                    #self.ax.set_zticklabels([])
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

            if self.refine_detdist:
                g[-3] += (S2*G2*self.detdist_deriv*one_minus_k_over_Lambda).sum()

            if self.refine_gain_fac:
                g[-2] += (2*self.gain_fac*(self.tilt_plane + S2*self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

            if self.refine_crystal_scale:
                # scale factor derivative
                g[-1] += (G2*2*self.scale_fac*self.model_bragg_spots * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        self.iterations += 1
        self.f_vals.append(f)
        return f, g

    def print_step(self, message, target):
        print("%s: residual=%2.7g, detdist=%2.7g, gain=%2.7g, scale=%2.7g"
              % (message, target, self.x[-3], self.x[-2], self.x[-1]))
