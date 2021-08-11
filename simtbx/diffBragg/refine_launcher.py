from __future__ import absolute_import, division, print_function
from simtbx.diffBragg.utils import map_hkl_list
from simtbx.diffBragg.refiners.local_refiner import LocalRefiner
import numpy as np
from dials.array_family import flex
from simtbx.diffBragg import utils
from copy import deepcopy
import os
import h5py
from dxtbx.model.experiment_list import ExperimentListFactory
import logging

LOGGER = logging.getLogger("main")


def local_refiner_from_parameters(refls, expt, params, miller_data=None):
    launcher = LocalRefinerLauncher(params)
    return launcher.launch_refiner(refls, expt, miller_data=miller_data)


class LocalRefinerLauncher:

    def __init__(self, params):
        self.params = self.check_parameter_integrity(params)

        self.rotXYZ_inits = {0: [0, 0, 0]}
        self.Modelers = {}

        self.panel_groups_refined = None
        self.panel_group_from_id = None
        self.panel_reference_from_id = None
        self.n_panel_groups = None
        self.verbose = False

        self.asu_from_idx = {}
        self.idx_from_asu = {}
        self.Hi = {}
        self.Hi_asu = {}
        self.num_hkl_global = 0
        self.Fref = None

        self.RUC = None
        self.SIM = None
        self.NCELLS_MASK = None
        self.NPIX_TO_ALLOC = -1
        self.DEVICE_ID = 0
        self.WATCH_MISORIENTATION = None
        self.SCALE_INIT_PER_SHOT = None
        self.NCELLS_INIT_PER_SHOT = None
        self.sausages_init = [0,0,0,1]

        self.n_local_unknowns = None
        self.n_global_unknowns = None
        self.local_idx_start = 0
        self.global_idx_start = None

        self.symbol = None

        self._rationalize_ncells_refinement_protocol()

    @property
    def NPIX_TO_ALLOC(self):
        return self._NPIX_TO_ALLOC

    @NPIX_TO_ALLOC.setter
    def NPIX_TO_ALLOC(self, val):
        assert val> 0 or val == -1
        self._NPIX_TO_ALLOC = int(val)

    def _alias_refiner(self):
        self._Refiner = LocalRefiner

    @staticmethod
    def check_parameter_integrity(params):
        if params.refiner.max_calls is None or len(params.refiner.max_calls) == 0:
            raise ValueError("Cannot refine because params.refiner.max_calls is empty")

        if os.environ.get("DIFFBRAGG_CUDA") is not None:
            params.refiner.use_cuda = True

        return params

    @staticmethod
    def _check_experiment_integrity(expt):
        for model in ["crystal", "detector", "beam", "imageset"]:
            if not hasattr(expt, model):
                raise ValueError("No %s in experiment, exiting. " % model)

    def launch_refiner(self, refls, expt, miller_data=None):
        self._alias_refiner()
        self._check_experiment_integrity(expt)
        self.create_cache_dir()

        self.verbose = self.params.refiner.verbose > 0
        self.DEVICE_ID = self.params.simulator.device_id

        shot_data = self.load_roi_data(refls, expt)
        if "rlp" in refls:
            shot_data.reso = {}

        if shot_data is None:
            raise ValueError("Cannot refine!")

        self._init_panel_group_information(expt.detector)

        self._init_simulator(expt, miller_data)
        UcellMan = utils.manager_from_crystal(self.SIM.crystal.dxtbx_crystal)

        if "rlp" in refls:
            self.shot_reso = {0: 1/np.linalg.norm(refls["rlp"], axis=1)}
        self.shot_ucell_managers = {0: UcellMan}
        self.shot_rois = {0: shot_data.rois}
        self.shot_nanoBragg_rois = {0: shot_data.nanoBragg_rois}
        self.shot_roi_imgs = {0: shot_data.roi_imgs}
        self.shot_roi_darkRMS = {0: shot_data.roi_darkRMS} if shot_data.roi_darkRMS else None
        self.shot_spectra = {0: self.SIM.beam.spectrum}
        self.shot_crystal_models = {0: expt.crystal}
        self.shot_crystal_model_refs = {0: deepcopy(expt.crystal)}
        self.shot_xrel = {0: shot_data.xrel}
        self.shot_yrel = {0: shot_data.yrel}
        self.shot_abc_inits = {0: shot_data.tilt_abc}
        self.shot_panel_ids = {0: shot_data.pids}
        self.shot_originZ_init = {0: 0}
        self.shot_selection_flags = {0: shot_data.selection_flags}
        self.shot_background = {0: shot_data.background}
        # TODO , get the experiment name in here, such that self.shot_expernames = {0: expername}
        self.symbol = expt.crystal.get_space_group().type().lookup_symbol()

        # <><><><><><>><><><><><><><><><>
        # determine number of parameters:
        # <><><><><><>><><><><><><><><><><>
        if not any(self.NCELLS_MASK):
            n_ncells_param = 3
        elif all(self.NCELLS_MASK):
            n_ncells_param = 1
        else:
            n_ncells_param = 2

        n_ncells_def_param = 3

        nrot_params = 3
        n_unitcell_params = len(UcellMan.variables)
        n_spotscale_params = 1
        n_originZ_params = 1
        n_eta_params = 3
        n_tilt_params = 3 * len(shot_data.nanoBragg_rois)
        n_sausage_params = 4*self.params.simulator.crystal.num_sausages
        n_perspot_param = len(shot_data.nanoBragg_rois)
        n_local_unknowns = nrot_params + n_unitcell_params + n_ncells_param + n_spotscale_params + n_originZ_params \
                           + n_tilt_params + n_eta_params + n_sausage_params + n_ncells_def_param + n_perspot_param
        n_local_unknowns += len(shot_data.nanoBragg_rois)
        self.panel_groups_refined = self.determine_refined_panel_groups(shot_data.pids, shot_data.selection_flags)

        #n_spectra_params = 0
        #if self.params.refiner.refine_spectra is not None and any(self.params.refiner.refine_spectra):
        n_spectra_params = 2
        n_panelRot_params = 3*self.n_panel_groups
        n_panelXYZ_params = 3*self.n_panel_groups
        n_global_unknowns = n_spectra_params + n_panelRot_params + n_panelXYZ_params
        self.n_ncells_param = n_ncells_param
        self.n_spectra_params = n_spectra_params
        self._launch(n_local_unknowns, n_global_unknowns, local_idx_start=0, global_idx_start=n_local_unknowns)
        return self.RUC

    def prep_for_pershot_Fcell_refinement(self, refls):
        observed_Hi_p1, _ = map(set, self.SIM.D.Fhkl_tuple)
        hkl_observed = observed_Hi_p1.intersection(set(refls["miller_index"]))
        missing_h = []
        for h in refls["miller_index"]:
            if h not in hkl_observed:
                missing_h.append(h)
        if missing_h:
            print(missing_h)
            raise Exception("The above indices dont have amplitudes in the input miller array")

        Hi = list(refls["miller_index"])
        self.Hi[0] = Hi
        self.Hi_asu[0] = map_hkl_list(Hi, True, self.symbol)

    def _launch(self, n_local_unknowns, n_global_unknowns, local_idx_start, global_idx_start):
        """
        Usually this method should be modified when new features are added to refinement
        """
        # TODO return None or refiner instance
        LOGGER.info("begin _launch")
        x_init = None
        nmacro = self.params.refiner.num_macro_cycles
        n_trials = len(self.params.refiner.max_calls)
        for i_trial in range(n_trials*nmacro):

            self.RUC = self._init_refiner(n_local_unknowns, n_global_unknowns, local_idx_start, global_idx_start)

            self.RUC.FNAMES = None # self.shot_expernames if self.shot_expernames else None

            self.RUC.print_end = self.params.refiner.print_end

            if not self.params.refiner.only_predict_model:
                if self.will_refine(self.params.refiner.refine_Bmatrix):
                    self.RUC.refine_Bmatrix = (self.params.refiner.refine_Bmatrix*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_Umatrix):
                    self.RUC.refine_Umatrix = (self.params.refiner.refine_Umatrix*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_ncells):
                    self.RUC.refine_ncells = (self.params.refiner.refine_ncells*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_ncells_def):
                    self.RUC.refine_ncells_def = (self.params.refiner.refine_ncells_def*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_bg):
                    self.RUC.refine_background_planes = (self.params.refiner.refine_bg*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_spot_scale):
                    self.RUC.refine_crystal_scale = (self.params.refiner.refine_spot_scale*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_spectra):
                    self.RUC.refine_spectra = (self.params.refiner.refine_spectra*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_detdist):
                    self.RUC.refine_detdist = (self.params.refiner.refine_detdist*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_panelZ):
                    self.RUC.refine_panelZ = (self.params.refiner.refine_panelZ*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_panelRotO):
                    self.RUC.refine_panelRotO = (self.params.refiner.refine_panelRotO*nmacro)[i_trial]
                if self.will_refine(self.params.refiner.refine_panelRotF):
                    self.RUC.refine_panelRotF = (self.params.refiner.refine_panelRotF*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_panelRotS):
                    self.RUC.refine_panelRotS = (self.params.refiner.refine_panelRotS*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_panelXY):
                    self.RUC.refine_panelXY = (self.params.refiner.refine_panelXY*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_per_spot_scale):
                    self.RUC.refine_per_spot_scale = (self.params.refiner.refine_per_spot_scale*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_eta):
                    self.RUC.refine_eta = (self.params.refiner.refine_eta*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_blueSausages):
                    self.RUC.refine_blueSausages = (self.params.refiner.refine_blueSausages*nmacro)[i_trial]

                if self.will_refine(self.params.refiner.refine_Fcell):
                    self.RUC.refine_Fcell = (self.params.refiner.refine_Fcell*nmacro)[i_trial]

            if self.RUC.refine_detdist and self.RUC.refine_panelZ:
                raise ValueError("Cannot refine panelZ and detdist simultaneously")

            if self.params.refiner.io.output_dir is not None:
                self.RUC.output_dir = self.params.refiner.io.output_dir

            self.RUC.parameter_hdf5_path = self.get_parameter_hdf5_path()

            self.RUC.break_signal = self.params.refiner.break_signal

            self.RUC.panel_group_from_id = self.panel_group_from_id

            self.RUC.panel_reference_from_id = self.panel_reference_from_id
            self.RUC.panel_groups_being_refined = self.panel_groups_refined
            self.RUC.panelRot_sigma = self.params.refiner.sensitivity.panelRotOFS
            self.RUC.panelX_sigma, self.RUC.panelY_sigma = self.params.refiner.sensitivity.panelXY
            self.RUC.panelZ_sigma = self.params.refiner.sensitivity.panelZ
            self.RUC.panelX_range = self.params.refiner.ranges.panel_X
            self.RUC.panelY_range = self.params.refiner.ranges.panel_Y
            self.RUC.panelZ_range = self.params.refiner.ranges.panel_Z
            self.RUC.panelRot_range = [[ang*np.pi/180 for ang in self.params.refiner.ranges.panel_rotO],
                                  [ang*np.pi/180 for ang in self.params.refiner.ranges.panel_rotF],
                                  [ang*np.pi/180 for ang in self.params.refiner.ranges.panel_rotS]]

            # TODO verify not refining Fcell in case of local refiner
            self.RUC.max_calls = (self.params.refiner.max_calls*nmacro)[i_trial]
            if self.params.refiner.refine_eta is not None and any(self.params.refiner.refine_eta):
                self.RUC.update_eta = any(self.params.refiner.refine_eta)

            self.RUC.x_init = x_init
            self.RUC.only_pass_refined_x_to_lbfgs = False
            self.RUC.bg_extracted = False
            self.RUC.save_model = self.params.refiner.save_models

            if self.params.refiner.refine_spectra is not None:
                self.RUC.update_spectra_during_refinement = any(self.params.refiner.refine_spectra)

            self.RUC.n_ncells_param = self.n_ncells_param
            self.RUC.recenter = True
            self.RUC.rescale_params = True

            self.RUC.ignore_line_search_failed_step_at_lower_bound = True

            self._initialize_some_refinement_parameters()

            # SIGMA VALUES
            self.RUC.rotX_sigma, self.RUC.rotY_sigma, self.RUC.rotZ_sigma = self.params.refiner.sensitivity.rotXYZ
            self.RUC.detector_distance_sigma = self.params.refiner.sensitivity.originZ
            self.RUC.ucell_sigmas = utils.unitcell_sigmas(self.Modelers[list(self.Modelers.keys())[0]].ucell_man,
                                                          self.params.refiner.sensitivity.unitcell)
            self.RUC.m_sigma = self.params.refiner.sensitivity.ncells_abc
            self.RUC.ncells_def_sigma = self.params.refiner.sensitivity.ncells_def
            self.RUC.spot_scale_sigma = self.params.refiner.sensitivity.spot_scale
            self.RUC.a_sigma, self.RUC.b_sigma, self.RUC.c_sigma = self.params.refiner.sensitivity.tilt_abc
            self.RUC.fcell_sigma_scale = self.params.refiner.sensitivity.fcell
            self.RUC.eta_sigma = self.params.refiner.sensitivity.eta
            self.RUC.eta_min = self.params.refiner.ranges.eta[0]
            self.RUC.eta_max = self.params.refiner.ranges.eta[1]
            if self.RUC.eta_min < 0:
                raise ValueError("Eta min should always be >= 0")

            self.RUC.print_all_missets = True

            self.RUC.n_spectra_param = self.n_spectra_params
            self.RUC.spectra_coefficients_sigma = self.params.refiner.sensitivity.spectra_coefficients  # .01, .01
            self.RUC.spectra_coefficients_init = self.params.refiner.init.spectra_coefficients  # 0, 1
            self.RUC.lambda_coef_ranges = [self.params.refiner.ranges.spectra0, self.params.refiner.ranges.spectra1]
            self.RUC.detector_distance_range = self.params.refiner.ranges.originZ

            Zvalues = [self.Modelers[i_exp].originZ_init for i_exp in self.Modelers]
            self.RUC.pershot_detdist_shifts = np.any(list(Zvalues))
            self.RUC.update_detector_during_refinement = False  # TODO is this ok ? maybe some tests will fail

            self.RUC.compute_image_model_correlation = self.params.refiner.compute_image_model_correlation

            # plot things
            self.RUC.sigma_r_scalar = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon
            self.RUC.trial_id = i_trial
            self.RUC.refine_rotZ = not self.params.refiner.fix_rotZ
            self.RUC.plot_images = self.params.refiner.plot.display
            self.RUC.plot_residuals = self.params.refiner.plot.as_residuals
            self.RUC.plot_stride = self.params.refiner.plot.iteration_stride

            # Fcell stuff
            self.RUC.rescale_fcell_by_resolution = self.params.refiner.rescale_fcell_by_resolution
            self.RUC.Fref = self.Fref
            self.RUC.merge_stat_frequency = self.params.refiner.stage_two.merge_stat_freq
            self.RUC.min_multiplicity = self.params.refiner.stage_two.min_multiplicity
            self.RUC.print_resolution_bins = self.params.refiner.stage_two.print_reso_bins
            if self.RUC.print_resolution_bins and self.RUC.Fref is None:
                raise ValueError("Fref cannot be None when printing resolution bins")
            if self.Hi:  # if miller indices stored per-shot
                self.RUC.Hi = self.Hi
            self.RUC.fcell_resolution_bin_Id = None
            self.RUC.log_fcells = True
            self.RUC.idx_from_asu = self.idx_from_asu
            self.RUC.asu_from_idx = self.asu_from_idx
            self.RUC.scale_r1 = True
            self.RUC.binner_dmax = self.params.refiner.stage_two.d_max
            self.RUC.binner_dmin = self.params.refiner.stage_two.d_min
            self.RUC.binner_nbin = self.params.refiner.stage_two.n_bin
            # end Fcell stuff
            self.RUC.rotXYZ_inits = self.rotXYZ_inits

            self.RUC.big_dump = self.params.refiner.big_dump
            self.RUC.request_diag_once = False
            self.RUC.S = self.SIM
            if self.params.refiner.mask is not None:
                self.RUC.MASK = utils.load_mask(self.params.refiner.mask)
            #if any() : self.RUC.refine_blueSausages is not None and any(self.RUC.refine_blueSausages):
            #TODO move to init_some_params method
            if self.params.refiner.refine_blueSausages is not None and any(self.params.refiner.refine_blueSausages):
                self.RUC.num_sausages = self.params.simulator.crystal.num_sausages
                self.RUC.sausages_init = self.sausages_init #{0: self.sausages_init}
                print("SAUSAGES!: ", self.sausages_init)
                #self.RUC.S.update_nanoBragg_instance("num_sausages", self.RUC.num_sausages)
            self.RUC.restart_file = self.params.refiner.io.restart_file
            self.RUC.has_pre_cached_roi_data = True
            self.RUC.trad_conv = True
            self.RUC.S.update_nanoBragg_instance('update_oversample_during_refinement',
                                                 self.params.refiner.update_oversample_during_refinement)
            self.RUC.S.update_nanoBragg_instance('use_cuda', self.params.refiner.use_cuda)
            self.RUC.S.update_nanoBragg_instance("Npix_to_allocate", self.NPIX_TO_ALLOC)
            if self.params.refiner.use_cuda:
                #TODO ensemble
                if self.params.refiner.randomize_devices:
                    self.RUC.S.update_nanoBragg_instance('device_Id', np.random.choice(self.params.refiner.num_devices))
                else:
                    self.RUC.S.update_nanoBragg_instance('device_Id', self.DEVICE_ID)
            self.RUC.refine_gain_fac = False
            self.RUC.use_curvatures_threshold = self.params.refiner.use_curvatures_threshold
            if not self.params.refiner.curvatures:
                self.RUC.S.update_nanoBragg_instance('compute_curvatures', False)
                self.RUC.S.update_nanoBragg_instance('verbose', self.params.refiner.verbose)
            self.RUC.calc_curvatures = self.params.refiner.curvatures
            self.RUC.poisson_only = self.params.refiner.poissononly
            self.RUC.trad_conv_eps = self.params.refiner.tradeps
            self.RUC.verbose = self.verbose
            if self.params.refiner.quiet:
                self.RUC.verbose = False
            # TODO optional properties.. make this obvious
            self.RUC.PROC_FNAMES = None
            self.RUC.PROC_IDX = None
            self.RUC.BBOX_IDX = None
            self.RUC.output_dir = self.params.refiner.io.output_dir
            #self.RUC.iterations=0
            LOGGER.info("_launch run setup")
            self.RUC.run(setup_only=True)
            LOGGER.info("_launch done run setup")
            # for debug purposes:
            #if not self.params.refiner.quiet:
            #    print("\n<><><><><><><><>TRIAL %d refinement status:" % i_trial)
            #    self.RUC.S.D.print_if_refining()

            self.RUC.num_positive_curvatures = 0
            self.RUC.use_curvatures = self.params.refiner.start_with_curvatures
            self.RUC.hit_break_to_use_curvatures = False

            # selection flags set here:
            #self.RUC.selection_flags = self.shot_selection_flags
            #if self.params.refiner.res_ranges is not None:
            #    assert self.shot_reso is not None, "cant set reso flags is rlp is not in refl tables"
            #    nshots = len(self.shot_selection_flags)
            #    more_sel_flags = {}
            #    res_ranges = utils.parse_reso_string(self.params.refiner.res_ranges)
            #    for i_shot in range(nshots):
            #        rhigh, rlow = (res_ranges*nmacro)[i_trial]
            #        sel_flags = self.shot_selection_flags[i_shot]
            #        res_flags = [rhigh < r < rlow for r in self.shot_reso[i_shot]]
            #        more_sel_flags[i_shot] = [flag1 and flag2 for flag1,flag2 in zip(sel_flags, res_flags)]
            #    self.RUC.selection_flags = more_sel_flags

            self.RUC.record_model_predictions = self.params.refiner.record_xy_calc

            #if True: #self.params.refiner.tryscipy:
            #    self.RUC.calc_curvatures = False
            #    self.RUC._setup()
            #    self.RUC.calc_func = True
            #    self.RUC.compute_functional_and_gradients()

            #    from scitbx.array_family import flex
            #    def func(x, RUC):
            #        RUC.calc_func = True
            #        RUC.x = flex.double(x)
            #        f, g = RUC.compute_functional_and_gradients()
            #        return f

            #    def fprime(x, RUC):
            #        RUC.calc_func = False
            #        RUC.x = flex.double(x)
            #        RUC.x = flex.double(x)
            #        f, g = RUC.compute_functional_and_gradients()
            #        return g.as_numpy_array()

            #    from scipy.optimize import fmin_l_bfgs_b
            #    out = fmin_l_bfgs_b(func=func, x0=np.array(self.RUC.x),
            #                        fprime=fprime, args=[self.RUC], factr=1e7)# args.scipyfactr)

            #else:
            LOGGER.info("_launcher runno setup")
            self.RUC.run(setup=False)
            LOGGER.info("_launcher done runno setup")
            if self.RUC.hit_break_to_use_curvatures:
                self.RUC.fix_params_with_negative_curvature = False
                self.RUC.num_positive_curvatures = 0
                self.RUC.use_curvatures = True
                self.RUC.run(setup=False)

            if self.RUC.hit_break_signal:
                if self.params.profile:
                    self.RUC.S.D.show_timings(self.RUC.rank)
                self.RUC._MPI_barrier()
                break

            if self.params.refiner.debug_pixel_panelfastslow is not None:
                utils.show_diffBragg_state(self.RUC.S.D, self.params.refiner.debug_pixel_panelfastslow)
                s = self.RUC._get_spot_scale(0)
                print("refiner spot scale=%f" % (s**2))

            x_init = self.RUC.x
            if self.params.refiner.only_predict_model:
                if self.RUC.gnorm > 0:
                    raise ValueError("Only predciting model, but the gradient is finite! This means the model changed, somethings wrong!")

            if self.params.profile:
                self.RUC.S.D.show_timings(self.RUC.rank)
            if self.params.refiner.use_cuda:
                self.RUC.S.D.gpu_free()


    def will_refine(self, param):
        return param is not None and any(param)

    def _initialize_some_refinement_parameters(self):
        # TODO verify spot scale init is squared or what
        self.RUC.spot_scale_init = {0: np.sqrt(self.params.refiner.init.spot_scale)}
        # eta_init in the refiner is a 3-tuple (6-tuple not yet supported)
        if self.params.simulator.crystal.anisotropic_mosaicity is None:
            eta_init = [self.params.simulator.crystal.mosaicity, 0, 0]
        else:
            eta_init = self.params.simulator.crystal.anisotropic_mosaicity
            if len(eta_init) == 6:
                raise NotImplementedError("No support for 6 parameter mosaic model")
            if self.params.simulator.crystal.crystal_for_anisotropic_mosaicity is not None:
                El = ExperimentListFactory.from_json_file(\
                    self.params.simulator.crystal.crystal_for_anisotropic_mosaicity, check_format=False)
                assert len(El)==1
                assert El[0].crystal is not None
                self.RUC.crystal_for_mosaicity_model = El[0].crystal
            else:
                self.RUC.crystal_for_mosaicity_model = deepcopy(self.SIM.crystal.dxtbx_crystal)
        self.RUC.eta_init = {0: eta_init}
        # MOSAICBLOCK
        m_init = self.params.simulator.crystal.ncells_abc
        if self.n_ncells_param == 2:
            INVERTED_NCELLS_MASK = [int(not mask_val) for mask_val in self.NCELLS_MASK]
            self.RUC.ncells_mask = INVERTED_NCELLS_MASK
            m_init = [m_init[i_ncell] for i_ncell in sorted(set(INVERTED_NCELLS_MASK))]
        self.RUC.m_init = {0: m_init}
        if self.params.refiner.ranges.ncells_abc is not None:
            self.RUC.m_range = self.params.refiner.ranges.ncells_abc
            self.RUC.use_Ncells_range = True

        ncells_def_init = self.params.simulator.crystal.ncells_def
        self.RUC.ncells_def_init = {0: ncells_def_init}
        self.RUC.ncells_def_range = self.params.refiner.ranges.ncells_def

        # UNITCELL
        self.RUC.ucell_inits = {0: self.shot_ucell_managers[0].variables}
        if self.params.refiner.ranges.ucell_edge_percentage is not None:
            names = self.shot_ucell_managers[0].variable_names
            maxs, mins = [], []
            for i_n, n in enumerate(names):
                val = self.shot_ucell_managers[0].variables[i_n]
                if "Ang" in n:
                    perc = self.params.refiner.ranges.ucell_edge_percentage * 0.01
                    valmin = val - val * perc
                    valmax = val + val * perc
                else:
                    deviation = self.params.refiner.ranges.ucell_angle_deviation
                    valmin = val - deviation / 2.
                    valmax = val + deviation / 2.
                mins.append(valmin)
                maxs.append(valmax)
            self.RUC.ucell_mins = {0: mins}
            self.RUC.ucell_maxs = {0: maxs}
            self.RUC.use_ucell_ranges = True


    def load_roi_data(self, refls, expt):
        background_mask = utils.load_mask(self.params.roi.background_mask)
        hotpix_mask = utils.load_mask(self.params.roi.hotpixel_mask)
        if hotpix_mask is not None:
            hotpix_mask = ~hotpix_mask

        img_data = utils.image_data_from_expt(expt)
        darkRMS_data = None
        if self.params.refiner.use_perpixel_dark_rms:
            darkRMS_data = utils.get_pedestalRMS_from_jungfrau(expt)
        if self.params.refiner.adu_per_photon is not None:
            img_data /= self.params.refiner.adu_per_photon
            if darkRMS_data is not None:
                darkRMS_data /= self.params.refiner.adu_per_photon
        img_data = img_data.astype("float32")

        # prepare ROI information
        panel_selection = [True]*len(refls)
        if self.params.roi.panels is not None:
            keeper_panels = utils.parse_panel_input_string(self.params.roi.panels)
            for i_pid, pid in enumerate(refls['panel']):
                panel_selection[i_pid] = pid in keeper_panels

        if self.params.roi.cachefile_dir is not None and not self.params.roi.make_cache_dir:
            rois, pids, tilt_abc, selection_flags, background = self.extract_roi_data_from_cachefile(expt)
        else:
            if darkRMS_data is not None:
                sigma_rdout = darkRMS_data
            elif self.params.refiner.sigma_r is not None:
                sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon
            else:
                sigma_rdout = 0
            # TODO unweighted fit
            roi_packet = utils.get_roi_background_and_selection_flags(
                refls, img_data, shoebox_sz=self.params.roi.shoebox_size,
                reject_edge_reflections=self.params.roi.reject_edge_reflections,
                reject_roi_with_hotpix=self.params.roi.reject_roi_with_hotpix,
                background_mask=background_mask, hotpix_mask=hotpix_mask,
                bg_thresh=self.params.roi.background_threshold,
                use_robust_estimation=not self.params.roi.fit_tilt,
                set_negative_bg_to_zero=self.params.roi.force_negative_background_to_zero,
                pad_for_background_estimation=self.params.roi.pad_shoebox_for_background_estimation,
                sigma_rdout=sigma_rdout, deltaQ=self.params.roi.deltaQ, experiment=expt,
                weighted_fit=self.params.roi.fit_tilt_using_weights)
            if self.params.roi.cachefile_dir is not None and self.params.roi.make_cache_dir:
                self.save_roi_data_to_file(expt, roi_packet)

            rois, pids, tilt_abc, selection_flags, background = roi_packet

        selection_flags = [sel1 and sel2 for sel1, sel2 in zip(selection_flags, panel_selection)]
        if not np.any(selection_flags):
            print("No spots for refinement")
            return None

        #nanoBragg_rois, xrel, yrel, roi_imgs, roi_gainmaps = [], [], [], [], []
        nanoBragg_rois, xrel, yrel, roi_imgs, roi_darkRMS = [], [], [], [], []
        for i_roi, (x1, x2, y1, y2) in enumerate(rois):
            nanoBragg_rois.append(((int(x1), int(x2)), (int(y1), int(y2))))
            yr, xr = np.indices((y2 - y1, x2 - x1))
            xrel.append(xr)
            yrel.append(yr)
            panel_id = pids[i_roi]
            roi_imgs.append(img_data[panel_id, y1:y2, x1:x2])
            if darkRMS_data is not None:
                roi_darkRMS.append(darkRMS_data[panel_id, y1:y2, x1:x2])

        shot_data = ShotData()
        shot_data.nanoBragg_rois = nanoBragg_rois
        shot_data.background = background
        shot_data.roi_imgs = roi_imgs
        shot_data.roi_darkRMS = roi_darkRMS if darkRMS_data is not None else None
        shot_data.pids = pids
        shot_data.rois = rois
        shot_data.xrel = xrel
        shot_data.yrel = yrel
        shot_data.tilt_abc = tilt_abc
        shot_data.selection_flags = selection_flags

        img_data = None
        del img_data

        return shot_data

    def _rationalize_ncells_refinement_protocol(self):
        self.NCELLS_MASK = utils.get_ncells_mask_from_string(self.params.refiner.ncells_mask)
        if all(self.NCELLS_MASK):
            if not self.params.simulator.crystal.has_isotropic_ncells:
                print("WARNING: NCELLS mask specifies isotropic ncells, but simulator.crystal.has_isotropic_ncells is set to False")
            self.params.simulator.crystal.has_isotropic_ncells = True
        else:
            if self.params.simulator.crystal.has_isotropic_ncells:
                print("WARNING: NCELLS mask specifies anisotropic ncells, but simulator.crystal.has_isotropic_ncells is set to True")
            self.params.simulator.crystal.has_isotropic_ncells = False

    def _prep_blue_sausages(self):
        if self.params.refiner.refine_blueSausages is not None and any(self.params.refiner.refine_blueSausages):
            from scitbx.array_family import flex
            init = self.params.refiner.init.sausages
            if init is None:
                init = []
                scale_fac = 0.1
                #if self.params.refiner.n_sausages_one_scale:
                #    scale_fac = 0   # NOTE if all sausage scales start at the same value, then the dervatives will all be equal
                #for i_sausage in range(self.params.simulator.crystal.num_sausages):
                #    init += [0, 0, 0, 1+scale_fac*i_sausage]
                init += [0, 0, 0, 1] * self.params.simulator.crystal.num_sausages
            else:
                assert len(init) % 4 == 0  # 4 parameters per mosaic texture block (sausage)

            self.SIM.D.update_number_of_sausages(self.params.simulator.crystal.num_sausages)
            x = flex.double(init[0::4])
            y = flex.double(init[1::4])
            z = flex.double(init[2::4])
            scale = flex.double(init[3::4])
            self.sausages_init = {0:init}
            self.SIM.D.set_sausages(x, y, z, scale)

    def _init_simulator(self, expt, miller_data):
        self.SIM = utils.simulator_from_expt_and_params(expt, self.params)
        # note self.SIM.D is a now diffBragg instance
        # include mosaic texture ?
        self._prep_blue_sausages()

        # update the miller data ?
        if miller_data is not None:
            self.SIM.crystal.miller_array = miller_data.as_amplitude_array()
            self.SIM.update_Fhkl_tuple()

    def _init_panel_group_information(self, detector):
        # default is one group per panel:
        self.panel_group_from_id = {pid: 0 for pid in range(len(detector))}
        # if specified in a file, then overwrite:
        if self.params.refiner.panel_group_file is not None:
            self.panel_group_from_id = utils.load_panel_group_file(self.params.refiner.panel_group_file)
            if not self.panel_group_from_id:
                raise ValueError("Loading panel group file %s  produced an panel group dict!"
                                 % self.params.refiner.panel_group_file)

        # how many unique panel groups :
        panel_groups = set(self.panel_group_from_id.values())
        self.n_panel_groups = len(panel_groups)

        # make dict where the key is the panel group id, and the value is a list of the panel ids within the group
        panels_per_group = {group_id: [] for group_id in panel_groups}
        for pid in self.panel_group_from_id:
            group_id = self.panel_group_from_id[pid]
            panels_per_group[group_id].append(pid)

        # we should rotate each panel in a group about the same reference point
        # Make a dict where key is panel id, and the  value is a reference origin
        self.panel_reference_from_id = {}
        for pid in self.panel_group_from_id:
            group_id = self.panel_group_from_id[pid]

            # take as reference, the origin of the first panel in the group
            reference_panel = detector[panels_per_group[group_id][0]]
            self.panel_reference_from_id[pid] = reference_panel.get_origin()

    def determine_refined_panel_groups(self, pids, selection_flags):
        refined_groups = []
        for i, pid in enumerate(pids):
            if selection_flags[i] and self.panel_group_from_id[pid] not in refined_groups:
                refined_groups.append(self.panel_group_from_id[pid])
        return refined_groups

    def _init_refiner(self, n_local_unknowns, n_global_unknowns, local_idx_start, global_idx_start):
        ref_crystals = None
        if self.WATCH_MISORIENTATION:
            ref_crystals = self.shot_crystal_model_refs
        RUC = self._Refiner(self.Modelers, self.symbol, self.params)
        return RUC

    def extract_roi_data_from_cachefile(self, expt):
        filename = self.roi_filename(expt)
        with h5py.File(filename, "r") as h5:
            rois = h5["rois"].value
            pids = h5["pids"].value
            tilt_abc = h5["tilt_abc"].value
            selection_flags = h5["selection_flags"].value
            background = h5["background"].value
        return rois, pids, tilt_abc, selection_flags, background


    def save_roi_data_to_file(self, expt, roi_packet):
        if self.params.roi.cachefile_dir is None:
            raise ValueError("Need to supply a cahcefile directory in order to store the roi data for faster re-loads")
        filename = self.roi_filename(expt)
        rois, pids, tilt_abc, selection_flags, background = roi_packet
        with h5py.File(filename, "w") as h5:
            h5.create_dataset("rois", data=np.array(rois), dtype=np.int32)
            h5.create_dataset("pids", data=pids, dtype=np.int32)
            h5.create_dataset("tilt_abc", data=np.array(tilt_abc), dtype=np.float32)
            h5.create_dataset("selection_flags", data=selection_flags, dtype=np.bool)
            h5.create_dataset("background", data=background, dtype=np.float32)

    def roi_filename(self, expt):
        basename = expt.imageset.get_path(0).replace("/", "+_+")
        filename = "_roi_data_%s-%d.h5" % (basename, expt.imageset.indices()[0])
        filename = os.path.join(self.params.roi.cachefile_dir, filename+".h5")
        return filename

    def create_cache_dir(self):
        if self.params.roi.cachefile_dir is not None and self.params.roi.make_cache_dir:
            if not os.path.exists(self.params.roi.cachefile_dir):
                #TODO multi rank make cache dir
                os.makedirs(self.params.roi.cachefile_dir)

    def get_parameter_hdf5_path(self):
        return self.params.refiner.parameter_hdf5_path


class ShotData:
    """for organizing per-shot information"""
    def __init__(self):
        self.nanoBragg_rois = None  # list of [((x1,x2),(y1,y2)),  ..]  # regions of interest in nanoBragg format
        self.background = None  # background data for this particular image
        self.roi_imgs = None  # list of image data regions of interest [ roi_1, roi_2 , ...]
        self.roi_darkRMS = None
        self.rois = None
        self.tilt_abc = None
        self.xrel = None
        self.yrel = None
        self.pids = None  # list of panel id
        self.selection_flags = None
        self.reso = None
