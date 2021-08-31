
from __future__ import division
import scitbx
from abc import ABCMeta, abstractproperty, abstractmethod
from scitbx.array_family import flex
from scitbx.lbfgs import core_parameters


# used in pixel refinement
class BreakToUseCurvatures(Exception):
    pass

class BreakBecauseSignal(Exception):
    pass

class ReachedMaxIterations(Exception):
    pass


class BaseRefiner:
    """
    This is the base class for pixel refinement
    """

    __metaclass__ = ABCMeta

    run_on_init = False

    def __init__(self):
        self.hit_break_signal = False  # internal flag in case of hitting the beak signal
        self.S = None   # instance of simtbx.nanoBragg.sim_data.SimData
        self.ignore_line_search_failed_step_at_lower_bound = False  # TODO: why is this sometimes necessary ?
        self.refine_lambda0 = False  # spectrea offset
        self.refine_lambda1 = False  # spectra scale factor
        self.m = 5  # LBFGS default core parameters
        self.maxfev = 20  # LBFGS default core parameters
        self.gtol = 0.9  # LBFGS default core parameters
        self.xtol = 1.e-16  # LBFGS default core parameters
        self.stpmin = 1.e-20  # LBFGS default core parameters
        self.stpmax = 1.e20  # LBFGS default core parameters
        self.trad_conv_eps = 0.05  # LBFGS terminator param converges whern |g| <= max(|x|,1) * trad_conv_eps
        self.drop_conv_max_eps = 1e-5  # LBFGS terminator param not sure, used in the other scitbx lbfgs convergence test
        self.mn_iter = 0  # LBFGS terminator param not sure used in lbfgs
        self.mx_iter = None  # LBFGS terminator param not sure used in lbfgs
        self.max_calls = 100000  # LBFGS terminator param how many overall iterations
        self.diag_mode = "always"  # LBFGS refiner property, whether to update curvatures at each iteration
        self.d = None   # place holder for a second derivative diagonal
        self.binner_dmin = 2  # if Fref is not None, then this defines R-factor and CC resolution binner
        self.binner_dmax = 999  # if Fref is not None, then this defines R-factor and CC resolution binner
        self.binner_nbin = 10  # if Fref is not None, then this defines R-factor and CC resolution binner
        self.Fref = None  # place holder for Fhkl reference (for computing R factors during refinement)
        self.use_curvatures = False  # whether to use the curvatures
        self.multi_panel = False  # whether the camera is multi panel or single panel
        self.hit_break_to_use_curvatures = False  # internal flag if calculating curvatures
        self.has_pre_cached_roi_data = False  # only for use in global refinement mode
        self.curv = None  # curvatures array used internally
        self.gradient_only = False  # parameter for LBFGS run method (internal to the Fortran code, see scitbx.lbfgs.__init__.py method run_c_plus_plus
        self.trad_conv = False  # traditional convergenve
        self.calc_curvatures = False  # whether to calc curvatures until a region of positive curvature is reached
        self.panel_ids = None  # list of panel_ids (same length as roi images, spot_rois, tilt_abc etc)
        self.poisson_only = False  # use strictly Poissonian statistics
        self._refinement_millers = None  # flex array of refinement miller indices (computed by GlobalRefiner _setup method)

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
    def S(self):
        """An instance of simtbx.nanoBragg.sim_data.SimData, the simulation workhorse"""
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
                    min_iterations=self.mn_iter,
                    max_iterations=self.mx_iter,
                    max_calls=self.max_calls)

    @property
    def _core_param(self):
        core_param = core_parameters()
        core_param.gtol = self.gtol
        core_param.xtol = self.xtol
        core_param.stpmin = self.stpmin
        core_param.stpmax = self.stpmax
        core_param.maxfev = self.maxfev
        core_param.m = self.m
        return core_param

    @property
    def _handler(self):
        return scitbx.lbfgs.exception_handling_parameters(
            ignore_line_search_failed_step_at_lower_bound=\
                self.ignore_line_search_failed_step_at_lower_bound)

    def _setup(self):
        """
        Optional place holder for class organization
        This is called just before running the minimizer
        Typically this involves populating the x array
        with initial values and configuring the diffBragg
        instance
        """
        pass

    def run(self, setup=True, setup_only=False):
        """runs the LBFGS minimizer"""

        if setup:
            self._setup()
            if setup_only:
                return

        if self.use_curvatures:
            try:
                self.minimizer = scitbx.lbfgs.run(
                    target_evaluator=self,
                    core_params=self._core_param,
                    exception_handling_params=self._handler,
                    termination_params=self._terminator,
                    gradient_only=self.gradient_only)
            except BreakToUseCurvatures:
                self.hit_break_to_use_curvatures = True
            except BreakBecauseSignal:
                self.hit_break_signal = True
                pass

        else:
            try:
                self.diag_mode = None
                self.minimizer = scitbx.lbfgs.run(
                    target_evaluator=self,
                    core_params=self._core_param,
                    exception_handling_params=self._handler,
                    termination_params=self._terminator,
                    gradient_only=self.gradient_only)
            except BreakToUseCurvatures:
                self.hit_break_to_use_curvatures = True
                pass
            except BreakBecauseSignal:
                self.hit_break_signal = True
                pass

    @property
    def use_curvatures(self):
        return self._use_curvatures

    @use_curvatures.setter
    def use_curvatures(self, val):
        if val:
            self.calc_curvatures = True
        self._use_curvatures = val

    def compute_functional_gradients_diag(self):
        self.f, self.g = self.compute_functional_and_gradients()
        self.d = flex.double(self.curv.as_numpy_array())
        self._verify_diag()
        return self.f, self.g, self.d

    def _verify_diag(self):
        sel = (self.g != 0)
        self.d.set_selected(~sel, 1000)
        assert self.d.select(sel).all_gt(0)
        self.d = 1 / self.d
