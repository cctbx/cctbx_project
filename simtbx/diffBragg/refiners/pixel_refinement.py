
import scitbx
import numpy as np
from abc import ABCMeta, abstractproperty, abstractmethod
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in
from scitbx.lbfgs.tst_mpi_split_evaluator import mpi_split_evaluator_run


# used in pixel refinement
class BreakToUseCurvatures(Exception):
    pass


class PixelRefinement(lbfgs_with_curvatures_mix_in):
    """
    This is the base class for pixel refinement based on
    ROI sub-images, where its understood that each sub-image
    contains Bragg scattering
    """

    __metaclass__ = ABCMeta

    run_on_init = False

    def __init__(self):
        self.use_curvatures = False  # whether to use the curvatures
        self.refine_background_planes = True  # whether to refine the background planes
        self.refine_gain_fac = False  # whether to refine the gain factor
        self.multi_panel = False  # whether the camera is multi panel or single panel
        self.split_evaluation = False  # whether to use split evaluation run method
        self.refine_ncells = False  # whether to refine Ncells abc
        self.hit_break_to_use_curvatures = False  # internal flag if calculating curvatures
        self.refine_detdist = False  # whether to refine the detdist
        self.refine_Amatrix = True  # whether to refine the  Amatrix (deprecated)
        self.refine_Bmatrix = True  # whether to refine the Bmatrx
        self.has_pre_cached_roi_data = False  # only for use in global refinement mode
        self.use_curvatures_threshold = 7   # how many positive curvature iterations required before breaking
        self.curv = None  # curvatures array used internally
        self.refine_Umatrix = True  # whether to refine the Umatrix
        self.verbose = True  # whether to print during iterations
        self.refine_crystal_scale = True  # whether to refine the crystal scale factor
        self.plot_images = False  # whether to plot images
        self.refine_rotX = True  # whether to refine the X rotation
        self.iterations = 0  # iteration counter , used internally
        self.refine_rotY = True  # whether to refine Y rotations
        self.refine_rotZ = True  # whether to refine Z rotations
        self.plot_residuals = False  # whether to plot residuals
        self.trad_conv = False  # traditional convergenve
        self.calc_curvatures = False  # whether to calc curvatures until a region of positive curvature is reached
        self.trad_conv_eps = 0.05  # converges whern |g| <= max(|x|,1) * trad_conv_eps
        self.plot_statistics = False  # whether to plot stats (global refinement mode)
        self.drop_conv_max_eps = 1e-5  # not sure, used in the other scitbx lbfgs convergence test
        self.mn_iter = None  # not sure used in lbfgs
        self.mx_iter = None  # not sure used in lbfgs
        self.max_calls = 1000  # how many overall iterations
        self.plot_stride = 10  # update plots after this many iterations
        self.ignore_line_search = False  # leave False, lbfgs
        self.panel_ids = None  # list of panel_ids (same length as roi images, spot_rois, tilt_abc etc)
        self.update_curvatures_every = 3  # every 3 consecutive all positive curvatures we will update them
        self.shot_idx = 0  # place holder because global refinement is across multiple shots
        self.shot_ids = None  # for global refinement ,

        self.poisson_only = True  # use strictly Poissonian statistics
        self.sigma_r = 3
        self.log2pi = np.log(np.pi*2)

        self.request_diag_once = False  # property of the parent class
        lbfgs_with_curvatures_mix_in.__init__(self, run_on_init=False)


    @property
    def _grad_accumulate(self):
        if self.poisson_only:
            return self._poisson_d
        else:
            return self._gaussian_d

    @property
    def _curv_accumulate(self):
        if self.poisson_only:
            return self._poisson_d2
        else:
            return self._gaussian_d2

    @property
    def _target_accumulate(self):
        if self.poisson_only:
            return self._poisson_target
        else:
            return self._gaussian_target

    def _poisson_target(self):
        pass

    def _gaussian_target(self):
        pass

    def _poisson_d(self, d):
        pass

    def _poisson_d2(self, d, d2):
        pass

    def _gaussian_d(self, d):
        pass

    def _gaussian_d2(self, d, d2):
        pass

    @abstractmethod
    def compute_functional_and_gradients(self):
        pass

    @abstractproperty
    def x(self):
        pass

    @abstractproperty
    def n(self):
        pass

    @property
    def plot_images(self):
        return self._plot_images

    @plot_images.setter
    def plot_images(self, val):
        if not isinstance(val, bool):
            raise ValueError("plot_images must be True or False")
        self._plot_images = val

    @property
    def spot_rois(self):
        return self._spot_rois

    @spot_rois.setter
    def spot_rois(self, val):
        if not len(val[0]) == 4:
            raise TypeError("spot_rois should be a list of 4-tuples")
        self._spot_rois = val

    @property
    def n_spots(self):
        try:
            return len(self._spot_rois)
        except AttributeError:
            return None

    @property
    def abc_init(self):
        return self._abc_init

    @abc_init.setter
    def abc_init(self, val):
        if not len(val[0]) == 3:
            raise TypeError("abc_init should be a list of 3-tuples")
        self._abc_init = val

    @property
    def img(self):
        return self._img

    @img.setter
    def img(self, val):
        self._img = val

    @property
    def S(self):
        """An instance of simtbx.diffBragg.sim_data2.SimData, the simulation workhorse"""
        return self._S

    @S.setter
    def S(self, val):
        if not hasattr(val, 'D'):
            print("S should be an instance of SimData after running SimData.instantiate_diffBragg()")
        self._S = val

    @property
    def trad_conv(self):
        return self._trad_conv_test

    @trad_conv.setter
    def trad_conv(self, val):
        self._trad_conv_test = val

    @property
    def trad_conv_eps(self):
        return self._trad_conv_eps

    @trad_conv_eps.setter
    def trad_conv_eps(self, val):
        self._trad_conv_eps = val

    @property
    def drop_conv_max_eps(self):
        return self._drop_conv_max_eps

    @drop_conv_max_eps.setter
    def drop_conv_max_eps(self, val):
        self._drop_conv_max_eps = val

    @property
    def mn_iter(self):
        return self._mn_iter

    @mn_iter.setter
    def mn_iter(self, val):
        self._mn_iter = val

    @property
    def mx_iter(self):
        return self._mx_iter

    @mx_iter.setter
    def mx_iter(self, val):
        self._mx_iter = val

    @property
    def max_calls(self):
        return self._max_calls

    @max_calls.setter
    def max_calls(self, val):
        self._max_calls = val

    @property
    def _terminator(self):
        return scitbx.lbfgs.termination_parameters(
                    traditional_convergence_test=self.trad_conv,
                    traditional_convergence_test_eps=self.trad_conv_eps,
                    drop_convergence_test_max_drop_eps=self.drop_conv_max_eps,
                    min_iterations=self.mx_iter,
                    max_iterations=self.mn_iter,
                    max_calls=self.max_calls)

    @property
    def _handler(self):
        return scitbx.lbfgs.exception_handling_parameters(
            ignore_line_search_failed_step_at_lower_bound=False)

    def _setup(self):
        """
        Optional place holder for class organization
        This is called just before running the minimizer
        Typically this involves populating the x array
        with initial values and configuring the diffBragg
        instance
        """
        pass

    def run(self, curvature_min_verbose=False, setup=True, cache_roi=True):
        """runs the LBFGS minimizer"""

        if setup:
            self._setup()
            self._cache_roi_arrays()

        # if not working in MPI mode (MPI used in global refinement)
        if not self.split_evaluation:
            if self.use_curvatures:
                self.minimizer = self.lbfgs_run(
                    target_evaluator=self,
                    min_iterations=self.mn_iter,
                    max_iterations=self.mx_iter,
                    traditional_convergence_test=self.trad_conv,
                    traditional_convergence_test_eps=self.trad_conv_eps,
                    use_curvatures=True,
                    verbose=curvature_min_verbose)

            else:
                try:
                    from scitbx.lbfgs import core_parameters
                    C = core_parameters()
                    C.gtol = 1
                    self.minimizer = scitbx.lbfgs.run(
                        target_evaluator=self,
                        #core_params=C,
                        exception_handling_params=self._handler,
                        termination_params=self._terminator)
                except BreakToUseCurvatures:
                    self.hit_break_to_use_curvatures = True
                    pass

        else:
            if self.use_curvatures:
                self.diag_mode = "always"
                self.minimizer = mpi_split_evaluator_run(
                    target_evaluator=self,
                    termination_params=self._terminator,
                    core_params=None,
                    exception_handling_params=None, log=None)
            else:
                try:
                    self.diag_mode = None
                    self.minimizer = mpi_split_evaluator_run(
                        target_evaluator=self,
                        termination_params=self._terminator,
                        core_params=None,
                        exception_handling_params=None, log=None)
                    # NOTE: best to leave log=None, not sure what would happen to log in MPI mode
                except BreakToUseCurvatures:
                    self.hit_break_to_use_curvatures = True
                    pass

    def _filter_spot_rois(self):
        """
        This is important to handle the edge case where an ROI occurs along
        the boundary of an image. This arises because
        NanoBragg assumes inclusive ROI bounds, but an exclusive raw_image
        """
        if self.multi_panel:
            nslow, nfast = self.img[0].shape
        else:
            nslow, nfast = self.img.shape
        for i, (_, x2, _, y2) in enumerate(self.spot_rois):
            if x2 == nfast:
                self.spot_rois[i][1] = x2-1  # update roi_xmax
            if y2 == nslow:
                self.spot_rois[i][3] = y2-1  # update roi_ymax

    def _cache_roi_arrays(self):
        """useful cache for iterative LBFGS step"""
        if self.has_pre_cached_roi_data:
            return
        nanoBragg_rois = []  # special nanoBragg format
        xrel, yrel, roi_img = [], [], []
        self._filter_spot_rois()
        for i_roi, (x1, x2, y1, y2) in enumerate(self.spot_rois):
            nanoBragg_rois.append(((x1, x2), (y1, y2)))
            yr, xr = np.indices((y2-y1+1, x2-x1+1))
            xrel.append(xr)
            yrel.append(yr)
            if self.multi_panel:
                pid = self.panel_ids[i_roi]
                roi_img.append(self.img[pid, y1:y2 + 1, x1:x2 + 1])
            else:
                roi_img.append(self.img[y1:y2+1, x1:x2+1])

        self.nanoBragg_rois = nanoBragg_rois
        self.roi_img = roi_img
        self.xrel = xrel
        self.yrel = yrel

    @property
    def roi_img(self):
        """the pixel value in each ROI (2D array)"""
        return self._roi_img

    @roi_img.setter
    def roi_img(self, val):
        self._roi_img = val

    @property
    def nanoBragg_rois(self):
        """rois in the nanoBragg format e.g. ((x1,x2), (y1,y2)"""
        return self._nanoBragg_rois

    @nanoBragg_rois.setter
    def nanoBragg_rois(self, val):
        self._nanoBragg_rois = val

    @property
    def xrel(self):
        """this is the fast-scan coordinate of each pixel in an ROI as 2D array"""
        return self._xrel

    @xrel.setter
    def xrel(self, val):
        self._xrel = val

    @property
    def yrel(self):
        """this is the slow-scan coordinate of each pixel in an ROI as 2D array"""
        return self._yrel

    @yrel.setter
    def yrel(self, val):
        self._yrel = val

    @property
    def use_curvatures(self):
        return self._use_curvatures

    @use_curvatures.setter
    def use_curvatures(self, val):
        if not isinstance(val, bool):
            raise ValueError("use_curvatures should be boolean")
        if val:
            self.calc_curvatures = True
        self._use_curvatures = val

    @property
    def refine_background_planes(self):
        return self._refine_background_planes

    @refine_background_planes.setter
    def refine_background_planes(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_background_planes = val

    @property
    def refine_Amatrix(self):
        return self._refine_Amatrix

    @refine_Amatrix.setter
    def refine_Amatrix(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_Amatrix = val

    @property
    def refine_Umatrix(self):
        return self._refine_Umatrix

    @refine_Umatrix.setter
    def refine_Umatrix(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_Umatrix = val

    @property
    def refine_Bmatrix(self):
        return self._refine_Bmatrix

    @refine_Bmatrix.setter
    def refine_Bmatrix(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_Bmatrix = val

    @property
    def refine_crystal_scale(self):
        return self._refine_crystal_scale

    @refine_crystal_scale.setter
    def refine_crystal_scale(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_crystal_scale = val

    @property
    def refine_gain_fac(self):
        return self._refine_gain_fac

    @refine_gain_fac.setter
    def refine_gain_fac(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_gain_fac = val

    @property
    def refine_detdist(self):
        return self._refine_detdist

    @refine_detdist.setter
    def refine_detdist(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_detdist = val

    @property
    def refine_ncells(self):
        return self._refine_ncells

    @refine_ncells.setter
    def refine_ncells(self, val):
        if not isinstance(val, bool):
            raise ValueError("refine background planes should be a boolean")
        self._refine_ncells = val


