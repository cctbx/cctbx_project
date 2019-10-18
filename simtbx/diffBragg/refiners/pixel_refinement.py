
import scitbx
import numpy as np
from abc import ABCMeta, abstractproperty, abstractmethod
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in


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
        self.use_curvatures = False
        self.refine_background_planes = True
        self.refine_gain_fac = False
        self.multi_panel = False
        self.refine_ncells = False
        self.hit_break_to_use_curvatures = False
        self.refine_detdist = True
        self.refine_Amatrix = True
        self.refine_Bmatrix = True
        self.use_curvatures_threshold = 7
        self.curv = None  # curvatures array
        self.refine_Umatrix = True
        self.verbose = True
        self.refine_crystal_scale = True
        self.plot_images = False
        self.trad_conv = False
        self.calc_curvatures = False
        self.trad_conv_eps = 0.05
        self.drop_conv_max_eps = 1e-5
        self.mn_iter = None
        self.mx_iter = None
        self.max_calls = 1000
        self.plot_stride = 10
        self.ignore_line_search = False
        self.panel_ids = None
        lbfgs_with_curvatures_mix_in.__init__(self, run_on_init=False)

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
        return len(self._spot_rois)

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

    def run(self, curvature_min_verbose=False, setup=True):
        """runs the LBFGS minimizer"""
        if setup:
            self._setup()
            self._cache_roi_arrays()
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

    def _cache_roi_arrays(self):
        """useful cache for iterative LBFGS step"""
        nanoBragg_rois = []  # special nanoBragg format
        xrel, yrel, roi_img = [], [], []
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


