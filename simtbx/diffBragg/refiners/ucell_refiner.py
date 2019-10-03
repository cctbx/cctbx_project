
from simtbx.diffBragg.refiners import RefineRot
from scitbx.array_family import flex
import pylab as plt
import numpy as np


class RefineUcell(RefineRot):

    def __init__(self, ucell_manager, *args, **kwargs):

        RefineRot.__init__(self, *args, **kwargs)
        self.ucell_manager = ucell_manager
        self.n_ucell_param = len(self.ucell_manager.variables)

    def _setup(self):
        self._setup_lbfgs_x_array()
        self._move_abc_init_to_x()
        self._set_diffBragg_instance()

    def _setup_lbfgs_x_array(self):
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
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i+3,
                self.ucell_manager.derivative_matrices[i])
            if self.use_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i+3, self.ucell_manager.second_derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.Bmatrix = self.ucell_manager.B_recipspace
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        self.ucell_derivatives = []
        self.ucell_second_derivatives = []
        for i in range(self.n_ucell_param):
            self.ucell_derivatives.append(self.D.get_derivative_pixels(3+i).as_numpy_array())
            if self.use_curvatures:
                self.ucell_second_derivatives.append(self.D.get_second_derivative_pixels(3+i).as_numpy_array())
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        #self.a = np.exp(self.x[i_spot])
        #self.b = np.exp(self.x[self.n_spots + i_spot])
        #self.c = np.exp(self.x[self.n_spots*2 + i_spot])
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots*2 + i_spot]


    def _update_orientation(self):
        self.ucell_manager.variables = list(self.x[self.n_spots*3:self.n_spots*3+self.n_ucell_param])
        self._send_gradients_to_derivative_managers()

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(len(self.x))
        if self.use_curvatures:
            self.curv = flex.double(len(self.x))
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
            self.II = Imeas
            f += (self.model_Lambda - Imeas*self.log_Lambda).sum()
            one_over_Lambda = 1./self.model_Lambda
            one_minus_k_over_Lambda = (1. - Imeas * one_over_Lambda)
            if self.use_curvatures:
                k_over_squared_Lambda = Imeas * one_over_Lambda * one_over_Lambda

            # compute gradients for background plane constants a,b,c
            xr = self.xrel[i_spot]  # fast scan pixels
            yr = self.yrel[i_spot]  # slow scan pixels
            #da = self.a*xr
            #db = self.b*yr
            #dc = self.c
            da = xr
            db = yr
            dc = 1
            if self.refine_background_planes:
                g[i_spot] += (da*one_minus_k_over_Lambda).sum()
                g[self.n_spots + i_spot] += (db*one_minus_k_over_Lambda).sum()
                g[self.n_spots*2 + i_spot] += (dc*one_minus_k_over_Lambda).sum()
            if self.use_curvatures and self.refine_background_planes:
                #da2 = da
                #db2 = db
                #dc2 = dc
                da2 = 0
                db2 = 0
                dc2 = 0
                self.curv[i_spot] += \
                    (da2*one_minus_k_over_Lambda + (da2**2) * k_over_squared_Lambda).sum()
                self.curv[self.n_spots + i_spot] += \
                    (db2*one_minus_k_over_Lambda + (db2**2) * k_over_squared_Lambda).sum()
                self.curv[self.n_spots*2 + i_spot] += \
                    (dc2*one_minus_k_over_Lambda + (dc2**2) * k_over_squared_Lambda).sum()

            if self.plot_images:
                m = Imeas[Imeas > 1e-9].mean()
                s = Imeas[Imeas > 1e-9].std()
                vmax = m+5*s
                vmin = m-s
                self.ax1.images[0].set_data(self.model_Lambda)
                self.ax1.images[0].set_clim(vmin, vmax)
                self.ax2.images[0].set_data(Imeas)
                self.ax2.images[0].set_clim(vmin, vmax)
                plt.suptitle("Iterations = %d, image %d / %d"
                             % (self.iterations, i_spot+1, self.n_spots))
                self.fig.canvas.draw()
                plt.pause(.02)

            # unit cell derivative
            for i_ucell_p in range(self.n_ucell_param):
                x_idx = self.n_spots * 3 + i_ucell_p
                d = self.ucell_derivatives[i_ucell_p]
                g[x_idx] += (one_minus_k_over_Lambda * d).sum()

                if self.use_curvatures:
                    d2 = self.ucell_second_derivatives[i_ucell_p]
                    self.curv[x_idx] += (d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda).sum()

            # scale factor derivative
            g[-1] += (2*self.scale_fac*self.model_bragg_spots * one_minus_k_over_Lambda).sum()
            if self.use_curvatures:
                self.curv[-1] += (self.model_bragg_spots * self.model_bragg_spots * k_over_squared_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        self.iterations += 1
        return f, g

    def print_step(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        labels = []
        for n, v in zip(names, vals):
            labels.append('%s=%2.7g' % (n, v))
        print ", ".join(labels)

    def curvatures(self):
        return self.curv


# # FIXME ? Maybe not ?
# # THIS IS BROKEN, CANT GET RSTBX parameter reduction TO WORKOUT T_T
#from cctbx import sgtbx, crystal_orientation
#from dxtbx.model.crystal import Crystal
#from rstbx.symmetry.constraints import parameter_reduction
#from scitbx.matrix import sqr
#from simtbx.diffBragg import utils
#import numpy as np
#class RefineUnitCellParameterReduction(RefineRot):
#
#    def __init__(self, symbol, *args, **kwargs):
#
#        RefineRot.__init__(self, *args, **kwargs)
#
#        self.symbol = symbol
#        self.point_group = sgtbx.space_group_info(symbol=self.symbol).group().build_derived_point_group()
#        self.param_ucell_tool = parameter_reduction.symmetrize_reduce_enlarge(self.point_group)
#        # TODO: get the a_real, b_real, c_real of the aligned crystal
#        self.a_real, self.b_real, self.c_real = \
#            sqr(
#                self.S.crystal.dxtbx_crystal.get_unit_cell().orthogonalization_matrix()).transpose().as_list_of_lists()
#        self._set_orientation()
#        self.init_param = self.param_ucell_tool.forward_independent_parameters()
#        self.n_ucell_param = len(self.init_param)
#        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1, ncols=2)
#
#    def get_da_dG(self):
#        matrices = self.param_ucell_tool.forward_gradients()
#        da_dG = []
#        for i in range(self.n_ucell_param):
#            da_dG.append(utils.lower_triangle(matrices[i]))
#        return da_dG
#
#    def _setup(self):
#        self._setup_lbfgs_x_array()
#        self._cache_roi_arrays()
#        self._move_abc_init_to_x()
#        self._set_diffBragg_instance()
#
#    def _set_orientation(self):
#        CO = crystal_orientation.crystal_orientation(
#            self.a_real + self.b_real + self.c_real,
#            crystal_orientation.basis_type.direct)
#        self.param_ucell_tool.set_orientation(CO)
#        # _ = self.param_ucell_tool.forward_independent_parameters()  # WHY?!
#
#    def _setup_lbfgs_x_array(self):
#        self.n_spots = len(self.spot_rois)
#        self.n_background_params = 3 * self.n_spots
#        self.n_gain_params = 1
#        self.n = self.n_background_params + self.n_ucell_param + self.n_gain_params
#        self.x = flex.double(self.n)
#        for i in range(self.n_ucell_param):
#            # NOTE: using exp of parameter
#            self.x[3 * self.n_spots + i] = np.log(self.init_param[i])
#        self.x[-1] = 1
#
#    def _set_diffBragg_instance(self):
#        self.D = self.S.D
#        self.D.initialize_managers()
#
#    def _run_diffBragg_current(self, i_spot):
#        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
#        self.D.Bmatrix = self.Brecip
#        self.D.add_diffBragg_spots()
#
#    def _extract_pixel_data(self):
#        # get the derivatives of the target w.r.t. the 6 unit cell parameters
#        self.dChi_da = []
#        for i in range(6):
#            self.dChi_da.append(self.D.get_derivative_pixels(3 + i).as_numpy_array())
#        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()
#
#    def _unpack_bgplane_params(self, i_spot):
#        self.a = self.x[i_spot]
#        self.b = self.x[self.n_spots + i_spot]
#        self.c = self.x[self.n_spots * 2 + i_spot]
#
#    def _update_orientation(self):
#        # self.param_ucell_tool = parameter_reduction.symmetrize_reduce_enlarge(self.point_group)
#        self.nicks_special_uc_params = [np.exp(i) for i in
#                                        self.x[self.n_spots * 3: self.n_spots * 3 + self.n_ucell_param]]
#        B = self.param_ucell_tool.backward_orientation(independent=self.nicks_special_uc_params)
#        self.Breal = B.direct_matrix()
#
#        # self.Brecip = B.reciprocal_matrix()
#        # NOTE: rstbx B matrix are all in lower diagonal format
#        self.a_real = self.Breal[0], self.Breal[1], self.Breal[2]
#        self.b_real = self.Breal[3], self.Breal[4], self.Breal[5]
#        self.c_real = self.Breal[6], self.Breal[7], self.Breal[8]
#        C = Crystal(self.a_real, self.b_real, self.c_real, self.symbol)
#        self.Brecip = C.get_B()
#        self._set_orientation()
#
#    def compute_functional_and_gradients(self):
#        f = 0
#        g = flex.double(len(self.x))
#        self._update_orientation()
#        self.scale_fac = self.x[-1]
#        self.da_dG = self.get_da_dG()
#        for i_spot in range(self.n_spots):
#            self._run_diffBragg_current(i_spot)
#            self._unpack_bgplane_params(i_spot)
#            self._set_background_plane(i_spot)
#            self._extract_pixel_data()
#            self._evaluate_averageI()
#            self._evaluate_log_averageI()
#
#            Imeas = self.roi_img[i_spot]
#            f += (self.model_Lambda - Imeas * self.log_Lambda).sum()
#            one_minus_k_over_Lambda = (1. - Imeas / self.model_Lambda)
#
#            # compute gradients for background plane constants a,b,c
#            xr = self.xrel[i_spot]  # fast scan pixels
#            yr = self.yrel[i_spot]  # slow scan pixels
#            g[i_spot] += (xr * one_minus_k_over_Lambda).sum()  # from handwritten notes
#            g[self.n_spots + i_spot] += (yr * one_minus_k_over_Lambda).sum()
#            g[self.n_spots * 2 + i_spot] += one_minus_k_over_Lambda.sum()
#            if self.plot_images:
#                self.ax1.clear()
#                self.ax2.clear()
#                self.ax1.imshow(self.model_Lambda)
#                self.ax2.imshow(Imeas)
#                self.ax1.images[0].set_clim(self.ax2.images[0].get_clim())
#                plt.suptitle("Spot %d / %d" % (i_spot + 1, self.n_spots))
#                plt.draw()
#                plt.pause(.2)
#
#            # do the unit cell derivatives, apply the  chain rule
#            # to get the derivative w.r.t. the symmetrized unit cell parameter G
#            for i_ucell_p in range(self.n_ucell_param):
#                # compute dChi_dG here using chain rule
#                dChi_dG = np.zeros_like(self.dChi_da[0])
#                for i in range(6):
#                    dChi_dG += (self.dChi_da[i] * self.da_dG[i_ucell_p][i])
#                # additional factor because using exponential of the G-tensor parameters
#                dG_dp = np.exp(self.x[self.n_spots * 3 + i_ucell_p])
#                g[self.n_spots * 3 + i_ucell_p] += (one_minus_k_over_Lambda * dChi_dG * dG_dp).sum()
#
#            # scale factor derivative
#            g[-1] += ((self.model_bragg_spots) * one_minus_k_over_Lambda).sum()
#
#        self.D.raw_pixels *= 0
#        self.print_step("LBFGS stp", f)
#        return f, g
#
#    def print_step(self, message, target):
#        ucell_t = self.D.unit_cell_tuple
#        dat = ucell_t + (self.x[-1],)
#        print ("unit cell %2.7g, %2.7g, %2.7g, %2.7g, %2.7g, %2.7g,  scale = %2.7g" % dat)
#        # print ("%2.7g %2.7g %2.7g %2.7g" % tuple(self.nicks_special_uc_params))
#        # print()
