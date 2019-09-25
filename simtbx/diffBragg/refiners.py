
from dxtbx.model.crystal import Crystal
import scitbx
from scitbx.array_family import flex
from scitbx.matrix import col
import numpy as np
import pylab as plt
from cctbx import sgtbx
from scitbx.matrix import sqr
from rstbx.symmetry.constraints import parameter_reduction
from cctbx import crystal_orientation
from simtbx.diffBragg import utils
from IPython import embed
import time


class RefineRot(object):

    def __init__(self, spot_rois, abc_init, img, SimData_instance,
                 plot_images=False):
        """

        :param spot_rois:
        :param abc_init:
        :param img:
        :param SimData_instance:
        :param plot_images:
        """

        self.plot_images = plot_images
        self.spot_rois = spot_rois
        self.abc_init = abc_init
        self.img = img
        self.S = SimData_instance
        self.trad_conv = False
        self.trad_conv_eps = 0.05
        self.drop_conv_max_eps = 1e-5
        self.mn_iter = None
        self.mx_iter = None
        self.max_calls = 1000
        self.ignore_line_search = False

    def _setup(self):
        # NOTE spot_rois are x1,x2,y1,y2 where x=fast, y=slow
        self._setup_lbfgs_x_array()
        self._cache_roi_arrays()
        self._move_abc_init_to_x()
        self._set_diffBragg_instance()

    def run(self):
        self._setup()
        self.terminator = scitbx.lbfgs.termination_parameters(
            traditional_convergence_test=self.trad_conv,
            traditional_convergence_test_eps=self.trad_conv_eps,
            drop_convergence_test_max_drop_eps=self.drop_conv_max_eps,
            min_iterations=self.mx_iter,
            max_iterations=self.mn_iter,
            max_calls=self.max_calls)

        self.handlers = scitbx.lbfgs.exception_handling_parameters(
            ignore_line_search_failed_step_at_lower_bound=self.ignore_line_search
        )

        self.minimizer = scitbx.lbfgs.run(
            target_evaluator=self,
            exception_handling_params=self.handlers,
            termination_params=self.terminator)

    def _setup_lbfgs_x_array(self):
        self.n_spots = len(self.spot_rois)
        self.n_background_params = 3*self.n_spots
        self.n_rot_params = 3
        self.n_gain_params = 1
        self.n = self.n_background_params + self.n_rot_params + self.n_gain_params
        self.x = flex.double(self.n)
        self.x[-1] = 1

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
        self.D = self.S.D
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

    def _extract_pixel_data(self):
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
        f = 0
        g = flex.double(len(self.x))
        for i_spot in range(self.n_spots):

            self._unpack_params(i_spot)
            self._run_diffBragg_current(i_spot)
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
                #plt.imshow(self.model_bragg_spots > 1e-6)
                plt.subplot(122)
                im2 = plt.imshow(Imeas)
                im.set_clim(im2.get_clim())
                plt.suptitle("Spot %d / %d" % (i_spot+1, self.n_spots))
                #plt.draw()
                plt.show()
                #plt.pause(.2)

            # rotation derivative
            g[self.n_spots*3] += (one_minus_k_over_Lambda * (self.dRotX)).sum()
            g[self.n_spots*3+1] += (one_minus_k_over_Lambda * (self.dRotY)).sum()
            g[self.n_spots*3+2] += (one_minus_k_over_Lambda * (self.dRotZ)).sum()

            # scale factor derivative
            g[-1] += ((self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        return f, g

    def print_step(self, message, target):
        print ("%s %10.4f" % (message, target),
               "[", " ".join(["%9.6f" % a for a in self.x]), "]")

    def get_correction_misset(self, as_axis_angle_deg=False):
        """
        return the current state of the perturbation matrix
        :return: scitbx.matrix sqr
        """
        angles = self.x[-4:-1]
        x = col((-1, 0, 0))
        y = col((0, -1, 0))
        z = col((0, 0, -1))
        RX = x.axis_and_angle_as_r3_rotation_matrix(angles[0], deg=False)
        RY = y.axis_and_angle_as_r3_rotation_matrix(angles[1], deg=False)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(angles[2], deg=False)
        M = RX*RY*RZ
        if as_axis_angle_deg:
            q = M.r3_rotation_matrix_as_unit_quaternion()
            rot_ang, rot_ax = q.unit_quaternion_as_axis_and_angle(deg=True)
            return rot_ang, rot_ax
        else:
            return M


# CLASS TO REFINE UNIT CELL PARAMETERS!
class RefineUnitCell(RefineRot):

    def __init__(self, symbol, *args, **kwargs):

        RefineRot.__init__(self, *args, **kwargs)

        self.symbol = symbol
        self.point_group = sgtbx.space_group_info(symbol=self.symbol).group().build_derived_point_group()
        self.param_ucell_tool = parameter_reduction.symmetrize_reduce_enlarge(self.point_group)
        # TODO: get the a_real, b_real, c_real of the aligned crystal
        self.a_real, self.b_real, self.c_real = \
            sqr(self.S.crystal.dxtbx_crystal.get_unit_cell().orthogonalization_matrix()).transpose().as_list_of_lists()
        self._set_orientation()
        self.init_param = self.param_ucell_tool.forward_independent_parameters()
        self.n_ucell_param = len(self.init_param)
        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1, ncols=2)

    def get_da_dG(self):
        matrices = self.param_ucell_tool.forward_gradients()
        da_dG = []
        for i in range(self.n_ucell_param):
            da_dG.append(utils.lower_triangle(matrices[i]))
        from IPython import embed
        embed()
        return da_dG

    def _setup(self):
        self._setup_lbfgs_x_array()
        self._cache_roi_arrays()
        self._move_abc_init_to_x()
        self._set_diffBragg_instance()

    def _set_orientation(self):
        CO = crystal_orientation.crystal_orientation(
            self.a_real + self.b_real + self.c_real,
            crystal_orientation.basis_type.direct)
        self.param_ucell_tool.set_orientation(CO)
        #_ = self.param_ucell_tool.forward_independent_parameters()  # WHY?!

    def _setup_lbfgs_x_array(self):
        self.n_spots = len(self.spot_rois)
        self.n_background_params = 3*self.n_spots
        self.n_gain_params = 1
        self.n = self.n_background_params + self.n_ucell_param + self.n_gain_params
        self.x = flex.double(self.n)
        for i in range(self.n_ucell_param):
            # NOTE: using exp of parameter
            self.x[3*self.n_spots + i] = np.log(self.init_param[i])
        self.x[-1] = 1

    def _set_diffBragg_instance(self):
        self.D = self.S.D
        self.D.initialize_managers()

    def _run_diffBragg_current(self, i_spot):
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self.D.Bmatrix = self.Brecip
        self.D.add_diffBragg_spots()

    def _extract_pixel_data(self):
        # get the derivatives of the target w.r.t. the 6 unit cell parameters
        self.dChi_da = []
        for i in range(6):
            self.dChi_da.append(self.D.get_derivative_pixels(3+i).as_numpy_array())
        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a = self.x[i_spot]
        self.b = self.x[self.n_spots + i_spot]
        self.c = self.x[self.n_spots*2 + i_spot]

    def _update_orientation(self):
        #self.param_ucell_tool = parameter_reduction.symmetrize_reduce_enlarge(self.point_group)
        self.nicks_special_uc_params = [np.exp(i) for i in self.x[self.n_spots*3: self.n_spots*3+self.n_ucell_param]]
        B = self.param_ucell_tool.backward_orientation(independent=self.nicks_special_uc_params)
        self.Breal = B.direct_matrix()

        #self.Brecip = B.reciprocal_matrix()
        # NOTE: rstbx B matrix are all in lower diagonal format
        self.a_real = self.Breal[0], self.Breal[1], self.Breal[2]
        self.b_real = self.Breal[3], self.Breal[4], self.Breal[5]
        self.c_real = self.Breal[6], self.Breal[7], self.Breal[8]
        C = Crystal(self.a_real, self.b_real, self.c_real, self.symbol)
        self.Brecip = C.get_B()
        self._set_orientation()

    def compute_functional_and_gradients(self):
        f = 0
        g = flex.double(len(self.x))
        self._update_orientation()
        self.scale_fac = self.x[-1]
        self.da_dG = self.get_da_dG()
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
                self.ax1.clear()
                self.ax2.clear()
                self.ax1.imshow(self.model_Lambda)
                self.ax2.imshow(Imeas)
                self.ax1.images[0].set_clim(self.ax2.images[0].get_clim())
                plt.suptitle("Spot %d / %d" % (i_spot+1, self.n_spots))
                plt.draw()
                plt.pause(.2)

            # do the unit cell derivatives, apply the  chain rule
            # to get the derivative w.r.t. the symmetrized unit cell parameter G
            for i_ucell_p in range(self.n_ucell_param):
                # compute dChi_dG here using chain rule
                dChi_dG = np.zeros_like(self.dChi_da[0])
                for i in range(6):
                    dChi_dG += (self.dChi_da[i] * self.da_dG[i_ucell_p][i])
                # additional factor because using exponential of the G-tensor parameters
                embed()
                dG_dp = np.exp(self.x[self.n_spots*3+i_ucell_p])
                g[self.n_spots*3+i_ucell_p] += (one_minus_k_over_Lambda * dChi_dG*dG_dp).sum()

            # scale factor derivative
            g[-1] += ((self.model_bragg_spots) * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        self.print_step("LBFGS stp", f)
        return f, g

    def print_step(self, message, target):
        ucell_t = self.D.unit_cell_tuple
        dat = ucell_t + (self.x[-1],)
        print ("unit cell %2.7g, %2.7g, %2.7g, %2.7g, %2.7g, %2.7g,  scale = %2.7g" % dat)
        #print ("%2.7g %2.7g %2.7g %2.7g" % tuple(self.nicks_special_uc_params))
        #print()


# CLASS TO REFINE UNIT CELL PARAMETERS!
class RefineTetragonal(RefineRot):

    def __init__(self, *args, **kwargs):

        RefineRot.__init__(self, *args, **kwargs)
        self.a_real, self.b_real, self.c_real = \
            sqr(self.S.crystal.dxtbx_crystal.get_unit_cell().orthogonalization_matrix()).transpose().as_list_of_lists()
        self.Tet = utils.TetragonalManager()
        self.Tet.a = self.a_real[0]
        self.Tet.c = self.c_real[2]
        self.S.crystal.dxtbx_crystal.show()
        self.n_ucell_param = 2

        # for plotting:
        self.iterations = 0
        self.all_a = []
        self.all_c = []
        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1,ncols=2)
        self.ax1.set_xlabel("iteration number")
        self.ax2.set_xlabel("iteration number")
        self.ax1.set_title("a $\AA$")
        self.ax2.set_title("c $\AA$")
        self.ax1.plot([0], [0], color='C1')
        self.ax2.plot([0], [0], color='C0')

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
        self.x[3*self.n_spots] = self.Tet.a
        self.x[3*self.n_spots+1] = self.Tet.c
        self.x[-1] = 0

    def _set_diffBragg_instance(self):
        self.D = self.S.D
        self.D.refine(3)
        self.D.refine(4)
        self.D.initialize_managers()

    def _send_gradients_to_derivative_managers(self):
        self.D.set_ucell_derivative_matrix(3, self.Tet.dB_da_real)
        self.D.set_ucell_derivative_matrix(4, self.Tet.dB_dc_real)

    def _run_diffBragg_current(self, i_spot):
        self.D.region_of_interest = self.nanoBragg_rois[i_spot]
        self._send_gradients_to_derivative_managers()
        self.D.Bmatrix = self.Tet.B_recipspace
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
        self.Tet.a = self.x[-3]
        self.Tet.c = self.x[-2]
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
            g[self.n_spots*3] += (one_minus_k_over_Lambda * (self.ucell_derivatives[0])).sum()
            g[self.n_spots*3+1] += (one_minus_k_over_Lambda * (self.ucell_derivatives[1])).sum()

            # scale factor derivative
            g[-1] += (self.model_bragg_spots * one_minus_k_over_Lambda).sum()

        self.D.raw_pixels *= 0
        #self.D.raw_pixels_roi *= 0
        self.print_step("LBFGS stp", f)
        return f, g

    def print_step(self, message, target):
        self.iterations += 1
        self.all_a.append(self.x[-3])
        self.all_c.append(self.x[-2])
        self.ax1.lines[0].set_data(range(self.iterations), self.all_a)
        self.ax2.lines[0].set_data(range(self.iterations), self.all_c)
        self.ax1.set_xlim(0, self.iterations+1)
        self.ax2.set_xlim(0, self.iterations+1)
        self.ax1.set_ylim(53, 58)
        self.ax2.set_ylim(74, 80)
        plt.draw()
        plt.pause(0.1)

        print ("unit cell a= %2.7g, c=%2.7g .. scale = %2.7g" % (self.x[-3], self.x[-2], self.x[-1]))
