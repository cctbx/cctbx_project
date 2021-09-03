from __future__ import absolute_import, division, print_function

from libtbx.mpi4py import MPI

COMM = MPI.COMM_WORLD
if not hasattr(COMM, "rank"):
    COMM.rank=0
    COMM.size=1
import time
import warnings
import signal
import logging

LOGGER = logging.getLogger("main")
warnings.filterwarnings("ignore")


class SignalHandler:
    def __init__(self):
        self.t = time.time()

    def handle(self, signum, frame):
        t = time.time()-self.t
        print("Recived signal ",signum," after program running for %f sec" % t)
        raise BreakBecauseSignal


SIGHAND = SignalHandler()


class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

try:
    import pandas
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# TODO : consider PEP-8 ing these numpy imports, but do a NERSC massively MPI time-test first...
# for now, if it aint broke, dont fix it ...
import numpy as np
import os

from simtbx.diffBragg.refiners import BreakBecauseSignal, BreakToUseCurvatures
from dials.array_family import flex
from simtbx.diffBragg.refiners import BaseRefiner
from collections import Counter
from cctbx import miller, sgtbx


class LocalRefiner(BaseRefiner):

    def __init__(self, shot_modelers, sgsymbol, params):
        BaseRefiner.__init__(self)

        self.params = params
        self.min_multiplicity = self.params.refiner.stage_two.min_multiplicity
        self.trad_conv_eps = self.params.refiner.tradeps
        self.calc_curvatures = self.params.refiner.curvatures
        self.poisson_only = self.params.refiner.poissononly
        self.break_signal = self.params.refiner.break_signal


        self.save_model_freq = 10  # save pixel model values after this many iterations
        self.saveZ_freq = 5  # save Zscore data every N iterations
        self.break_signal = None  # check for this signal during refinement, and break refinement if signal is received (see python signal module)
        self.print_end = "\n"  # value appended to the end of each printed string
        self.refine_blueSausages = False  # refine multiple crystals per image (e.g. James Holtons blue sausage plot)
        self.refine_eta = False  # refine the mosaic spread angle
        self.refine_per_spot_scale = False  # experimental, refine a per spot scale factor for each ROI
        self.save_model = False  # whether to save the model
        self.idx_from_asu = {}  # maps global fcell index to asu hkl
        self.asu_from_idx = {}  # maps asu hkl to global fcell index
        self.rescale_params = True  # whether to rescale parameters during refinement  # TODO this will always be true, so remove the ability to disable
        self.fcell_sigma_scale = 0.005  # sensitivity for Fcell during refinement
        self.pause_after_iteration = 0.001  # pause for this long after each iteration (not used currently)
        self.request_diag_once = False  # LBFGS refiner property
        self.output_dir = self.params.refiner.io.output_dir  # directory to dump progress files, these can be used to restart simulation later
        self.min_multiplicity = 1  # only refine a spots Fhkl if multiplicity greater than this number
        self.restart_file = None  # output file from previous run refinement
        self.trial_id = 0  # trial id in case multiple trials are run in sequence
        self.x_init = None  # used to restart the refiner (e.g. self.x gets updated with this)
        self.log_fcells = True  # to refine Fcell using logarithms to avoid negative Fcells
        self.refine_background_planes = False  # whether to refine the background planes
        self.refine_ncells = False  # whether to refine Ncells abc
        self.refine_ncells_def = False  # whether to refine Ncells abc
        self.refine_detdist = False  # whether to refine the detdist
        self.refine_panelXY = False  # whether to refine panel origin X and Y components
        self.refine_panelZ = False  # whether to refine panel origin X and Y components
        self.refine_panelRotO = False  # whether to refine the panel rotation
        self.refine_panelRotF = False  # whether to refine the panel rotation
        self.refine_panelRotS = False  # whether to refine the panel rotation
        self.refine_Umatrix = False  # whether to refine the Umatrix
        self.refine_Bmatrix = False  # whether to refine the Bmatrx
        self.refine_crystal_scale = False  # whether to refine the crystal scale factor
        self.refine_Fcell = False  # whether to refine Fhkl for each shoebox ROI
        self.use_curvatures_threshold = 7  # how many positive curvature iterations required before breaking, after which simulation can be restart with use_curvatures=True
        self.verbose = True  # whether to print during iterations
        self.iterations = 0  # iteration counter , used internally
        self.shot_ids = None  # for global refinement ,
        self.sigma_r = 3.  # readout noise mean in ADU
        self.log2pi = np.log(np.pi*2)

        self._sig_hand = None  # method for handling the break_signal, e.g. SIGHAND.handle defined above (theres an MPI version in global_refiner that overwrites this in stage 2)
        self._is_trusted = None  # used during refinement, 1-D array or trusted pixels corresponding to the pixels in the ROI

        self.rank = COMM.rank
        self.save_model_freq = self.params.refiner.stage_two.save_model_freq
        self.use_nominal_h = self.params.refiner.stage_two.use_nominal_hkl
        self.init_image_corr = None

        self.Modelers = shot_modelers
        self.shot_ids = sorted(self.Modelers.keys())
        self.n_shots = len(shot_modelers)
        self.n_shots_total = COMM.bcast(COMM.reduce(self.n_shots))
        LOGGER.debug("Loaded %d shots across all ranks" % self.n_shots_total)
        self.f_vals = []  # store the functional over time

        self._ncells_id = 9  # diffBragg internal index for Ncells derivative manager
        self._detector_distance_id = 10  # diffBragg internal index for detector_distance derivative manager
        self._panelRotO_id = 14  # diffBragg internal index for derivative manager
        self._panelRotF_id = 17  # diffBragg internal index for derivative manager
        self._panelRotS_id = 18  # diffBragg internal index for derivative manager
        self._panelX_id = 15  # diffBragg internal index for  derivative manager
        self._panelY_id = 16  # diffBragg internal index for  derivative manager
        self._fcell_id = 11  # diffBragg internal index for Fcell derivative manager
        self._eta_id = 19  # diffBragg internal index for eta derivative manager
        self._lambda0_id = 12  # diffBragg interneal index for lambda derivatives
        self._lambda1_id = 13  # diffBragg interneal index for lambda derivatives
        self._sausage_id = 20
        self._ncells_def_id = 21

        self.symbol = sgsymbol
        self.space_group = sgtbx.space_group(sgtbx.space_group_info(symbol=self.symbol).type().hall_symbol())
        self.I_AM_ROOT = COMM.rank==0

    def print(self, s, *args, **kwargs):
        """cheap logger"""
        if self.verbose:
            if isinstance(s, str):
                for line in s.split("\n"):
                    print(line, *args, **kwargs, end=self.print_end)
            else:
                print(s, *args, **kwargs, end=self.print_end)

    def __call__(self, *args, **kwargs):
        _, _ = self.compute_functional_and_gradients()
        return self.x, self._f, self._g, self.d

    @property
    def n(self):
        """LBFGS property"""
        return len(self.x)

    @property
    def n_global_fcell(self):
        return len(self.idx_from_asu)

    @property
    def image_shape(self):
        panelXdim, panelYdim = self.S.detector[0].get_image_size()
        Npanels = len(self.S.detector)
        return Npanels, panelYdim, panelXdim

    @property
    def x(self):
        """LBFGS parameter array"""
        return self._x

    @x.setter
    def x(self, val):
        self._x = val

    def _check_keys(self, shot_dict):
        """checks that the dictionary keys are the same"""
        if not sorted(shot_dict.keys()) == self.shot_ids:
            raise KeyError("input data funky, check GlobalRefiner inputs")
        return shot_dict

    def _evaluate_averageI(self):
        """model_Lambda means expected intensity in the pixel"""
        self.model_Lambda = self.Modelers[self._i_shot].all_background + self.model_bragg_spots

    def make_output_dir(self):
        if self.I_AM_ROOT and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.Zdir = os.path.join(self.output_dir, "Z")
        self.model_dir = os.path.join(self.output_dir, "model")
        for dirname in (self.Zdir, self.model_dir):
            if self.I_AM_ROOT and not os.path.exists(dirname):
                os.makedirs(dirname)
        COMM.barrier()

    def _setup(self):
        # Here we go!  https://youtu.be/7VvkXA6xpqI
        LOGGER.info("Setup begins!")
        if self.refine_Fcell and not self.asu_from_idx:
            raise ValueError("Need to supply a non empty asu from idx map")
        if self.refine_Fcell and not self.idx_from_asu:  # # TODO just derive from its inverse
            raise ValueError("Need to supply a non empty idx from asu map")

        self.make_output_dir()

        self.shot_mapping = self._get_shot_mapping()
        self.n_total_shots = len(self.shot_mapping)

        test_shot = self.shot_ids[0]
        N_PARAM_PER_SHOT = 2
        self.n_ucell_param = len(self.Modelers[test_shot].PAR.ucell_man.variables)
        self.n_total_params = self.n_total_shots*N_PARAM_PER_SHOT + self.n_global_fcell

        self.spot_scale_xpos = {}
        self.Bfactor_xpos = {}
        for shot_id in self.shot_ids:
            self.spot_scale_xpos[shot_id] = self.shot_mapping[shot_id]*N_PARAM_PER_SHOT
            self.Bfactor_xpos[shot_id] = self.shot_mapping[shot_id]*N_PARAM_PER_SHOT + 1

        LOGGER.info("--0 create an Fcell mapping")
        if self.refine_Fcell:
            #idx, data = self.S.D.Fhkl_tuple
            #self.idx_from_p1 = {h: i for i, h in enumerate(idx)}
            self._make_p1_equiv_mapping()
            # self.p1_from_idx = {i: h for i, h in zip(idx, data)}

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        self.x = flex.double(np.ones(self.n_total_params))
        LOGGER.info("--Setting up per shot parameters")

        self.fcell_xstart = self.n_total_shots*N_PARAM_PER_SHOT

        self.hkl_totals = []
        if self.refine_Fcell:
            for i_shot in self.shot_ids:
                for i_h, h in enumerate(self.Modelers[i_shot].Hi_asu):
                    self.hkl_totals.append(self.idx_from_asu[h])
            self.hkl_totals = self._MPI_reduce_broadcast(self.hkl_totals)

        self._MPI_setup_global_params()
        self._MPI_sync_fcell_parameters()
        # reduce then broadcast fcell
        LOGGER.info("--combining parameters across ranks")
        self._MPI_sync_hkl_freq()

        if self.x_init is not None:
            LOGGER.info("Initializing with provided x_init array")
            self.x = self.x_init
        elif self.restart_file is not None:
            LOGGER.info("Restarting from parameter file %s" % self.restart_file)
            self.x = flex.double(np.load(self.restart_file)["x"])

        # setup the diffBragg instance
        self.D = self.S.D

        self.D.refine(self._fcell_id)
        self.D.initialize_managers()

        for sid in self.shot_ids:
            Modeler = self.Modelers[sid]
            Modeler.all_fcell_global_idx = np.array([self.idx_from_asu[h] for h in Modeler.hi_asu_perpix])
            Modeler.unique_i_fcell = set(Modeler.all_fcell_global_idx)
            Modeler.i_fcell_slices = self._get_i_fcell_slices(Modeler)
            self.Modelers[sid] = Modeler  # TODO: VERIFY IF THIS IS NECESSARY ?

        self._MPI_barrier()
        LOGGER.info("Setup ends!")

    def _get_i_fcell_slices(self, Modeler):
        """finds the boundaries for each fcell in the 1-D array of per-shot data"""
        # TODO move this to Data Modeler class ?
        splitter = np.where(np.diff(Modeler.all_fcell_global_idx) != 0)[0]+1
        npix = len(Modeler.all_fcell_global_idx)
        slices = [slice(V[0], V[-1]+1, 1) for V in np.split(np.arange(npix), splitter)]
        i_fcells = [V[0] for V in np.split(Modeler.all_fcell_global_idx, splitter)]
        i_fcell_slices = {}
        for i_fcell, slc in zip(i_fcells, slices):
            if i_fcell not in i_fcell_slices:
                i_fcell_slices[i_fcell] = [slc]
            else:
                i_fcell_slices[i_fcell].append(slc)
        return i_fcell_slices

    def _get_shot_mapping(self):
        """each modeled shot maps to an integer along interval [0,Nshots) """
        all_shot_ids = COMM.gather(self.shot_ids)
        shot_mapping = None
        if COMM.rank == 0:
            unique_shot_ids = set([sid for shot_ids in all_shot_ids for sid in shot_ids])
            shot_mapping = {shot_id: i_shot for i_shot, shot_id in enumerate(unique_shot_ids)}
        shot_mapping = COMM.bcast(shot_mapping)
        return shot_mapping

    def _make_p1_equiv_mapping(self):
        self.num_equivs_for_i_fcell = {}
        self.update_indices = []
        for i_fcell in range(self.n_global_fcell):
            hkl_asu = self.asu_from_idx[i_fcell]

            equivs = [i.h() for i in miller.sym_equiv_indices(self.space_group, hkl_asu).indices()]
            self.num_equivs_for_i_fcell[i_fcell] = len(equivs)
            self.update_indices += equivs
        self.update_indices = flex.miller_index(self.update_indices)

    def _MPI_setup_global_params(self):
        if self.I_AM_ROOT:
            self.print("--2 Setting up global parameters")
            if self.output_dir is not None:
                np.save(os.path.join(self.output_dir, "f_asu_map"), self.asu_from_idx)

            self._setup_fcell_params()

    def _setup_fcell_params(self):
        if self.refine_Fcell:
            self.print("----loading fcell data")
            # this is the number of observations of hkl (accessed like a dictionary via global_fcell_index)
            self.print("---- -- counting hkl totes")
            LOGGER.info("compute HKL multiplicity")
            self.hkl_frequency = Counter(self.hkl_totals)
            LOGGER.info("save HKL multiplicity")
            np.save(os.path.join(self.output_dir, "f_asu_multi"), self.hkl_frequency)
            LOGGER.info("Done ")

            LOGGER.info("local refiner symbol=%s ; nanoBragg crystal symbol: %s" % (self.symbol, self.S.crystal.symbol))
            ma = self.S.crystal.miller_array_high_symmetry.map_to_asu()
            LOGGER.info("make an Fhkl map")
            ma_map = {h: d for h,d in zip(ma.indices(), ma.data())}
            LOGGER.info("make fcell_init")
            self.fcell_init_from_i_fcell = np.array([ma_map[self.asu_from_idx[i_fcell]] for i_fcell in range(self.n_global_fcell)])
            self.fcell_sigmas_from_i_fcell = 1
            LOGGER.info("DONE make fcell_init")

    def _get_sausage_parameters(self, i_shot):
        pass

    def _get_rotXYZ(self, i_shot):
        vals = [self.Modelers[i_shot].RotXYZ[i_rot].init for i_rot in range(3)]
        return vals

    def _get_rotX(self, i_shot):
        pass

    def _get_rotY(self, i_shot):
        pass

    def _get_rotZ(self, i_shot):
        pass

    def _get_spectra_coefficients(self):
        pass

    def _get_ucell_vars(self, i_shot):
        vars = []
        for i in range(self.n_ucell_param):
            var = self.Modelers[i_shot].PAR.ucell[i].init
            vars.append(var)
        return vars

    def _get_panelRot_val(self, panel_id):
        pass

    def _get_panelXYZ_val(self, panel_id, i_shot=0):
        pass

    def _get_detector_distance_val(self, i_shot):
        return self.Modelers[i_shot].PAR.detz_shift.init

    def _get_ncells_def_vals(self, i_shot):
        pass

    def _get_m_val(self, i_shot):
        vals = [self.Modelers[i_shot].PAR.Nabc[i_N].init for i_N in range(3)]
        return vals

    def _get_eta(self, i_shot):
        pass

    def _get_spot_scale(self, i_shot):
        xval = self.x[self.spot_scale_xpos[i_shot]]
        PAR = self.Modelers[i_shot].PAR
        sig = PAR.Scale.sigma
        init = PAR.Scale.init
        val = sig*(xval-1) + init
        return val

    def _get_bfactor(self, i_shot):
        xval = self.x[self.Bfactor_xpos[i_shot]]
        PAR = self.Modelers[i_shot].PAR
        sig = PAR.B.sigma
        init = PAR.B.init
        val = sig*(xval-1) + init
        return val

    def _get_bg_vals(self, i_shot, i_spot):
        pass

    def _send_ucell_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        pass

    def _run_diffBragg_current(self):
        LOGGER.info("run diffBragg for shot %d" % self._i_shot)
        pfs = self.Modelers[self._i_shot].pan_fast_slow
        if self.use_nominal_h:
            nom_h = self.Modelers[self._i_shot].all_nominal_hkl
            self.D.add_diffBragg_spots(pfs, nom_h)
        else:
            self.D.add_diffBragg_spots(pfs)
        LOGGER.info("finished diffBragg for shot %d" % self._i_shot)

    def _store_updated_Fcell(self):
        if not self.refine_Fcell:
            return
        xvals = self.x[self.fcell_xstart: self.fcell_xstart+self.n_global_fcell]
        if self.rescale_params and self.log_fcells:
            sigs = self.fcell_sigmas_from_i_fcell
            inits = self.fcell_init_from_i_fcell
            if self.log_fcells:
                vals = np.exp(sigs*(xvals - 1))*inits
            else:
                vals = sigs*(xvals - 1) + inits
                vals[vals < 0] = 0
        else:
            if self.log_fcells:
                vals = np.exp(xvals)
            else:
                vals = xvals
                vals [vals < 0] = 0
        self._fcell_at_i_fcell = vals

    def _update_Fcell(self):
        if not self.refine_Fcell:
            return
        update_amps = []
        for i_fcell in range(self.n_global_fcell):
            new_Fcell_amplitude = self._fcell_at_i_fcell[i_fcell]
            update_amps += [new_Fcell_amplitude] * self.num_equivs_for_i_fcell[i_fcell]

        update_amps = flex.double(update_amps)
        self.S.D.quick_Fhkl_update((self.update_indices, update_amps))

    def _update_spectra_coefficients(self):
        pass

    def _update_eta(self):
        pass

    def _set_background_plane(self):
        self.tilt_plane = self.Modelers[self._i_shot].all_background[self.roi_sel]

    def _update_sausages(self):
        pass

    def _update_rotXYZ(self):
        pass

    def _update_ncells(self):
        vals = self._get_m_val(self._i_shot)
        self.D.set_ncells_values(tuple(vals))

    def _update_ncells_def(self):
        pass

    def _update_dxtbx_detector(self):
        shiftZ = self._get_detector_distance_val(self._i_shot)
        self.S.D.shift_origin_z(self.S.detector,  shiftZ)

    def _extract_spectra_coefficient_derivatives(self):
        pass

    def _pre_extract_deriv_arrays(self):
        npix = len(self.Modelers[self._i_shot].all_data)
        self._model_pix = self.D.raw_pixels_roi[:npix].as_numpy_array()

        if self.refine_Fcell:
            dF = self.D.get_derivative_pixels(self._fcell_id)
            self._extracted_fcell_deriv = dF[:npix].as_numpy_array()
            if self.calc_curvatures:
                d2F = self.D.get_second_derivative_pixels(self._fcell_id)
                self._extracted_fcell_second_deriv = d2F[:npix].as_numpy_array()

    def _extract_sausage_derivs(self):
        pass

    def _extract_Umatrix_derivative_pixels(self):
        pass

    def _extract_Bmatrix_derivative_pixels(self):
        pass

    def _extract_ncells_def_derivative_pixels(self):
        pass

    def _extract_mosaic_parameter_m_derivative_pixels(self):
        pass

    def _extract_detector_distance_derivative_pixels(self):
        pass

    def _extract_panelRot_derivative_pixels(self):
        pass

    def _extract_panelXYZ_derivative_pixels(self):
        pass

    def _extract_Fcell_derivative_pixels(self):
        # TODO pre-extract
        self.fcell_deriv = self.fcell_second_deriv = 0
        if self.refine_Fcell:
            SG = self.scale_fac
            self.fcell_deriv = SG*(self._extracted_fcell_deriv)
            # handles Nan's when Fcell is 0 for whatever reason
            if self.calc_curvatures:
                self.fcell_second_deriv = SG*self._extracted_fcell_second_deriv

    def _get_per_spot_scale(self, i_shot, i_spot):
        pass

    def _extract_pixel_data(self):
        #Mod = self.Modelers[self._i_shot]
        #self.Bfactor_qterm = Mod.all_q_perpix**2 / 4.
        #self._expBq = np.exp(-self.b_fac**2 * self.Bfactor_qterm)
        #self.model_bragg_spots = self._expBq*self.scale_fac*(self._model_pix)
        self.model_bragg_spots = self.scale_fac*self._model_pix
        self._extract_Fcell_derivative_pixels()

    def _update_ucell(self):
        self.D.Bmatrix = self.Modelers[self._i_shot].PAR.Bmatrix

    def _update_umatrix(self):
        self.D.Umatrix = self.Modelers[self._i_shot].PAR.Umatrix

    def _update_beams(self):
        # sim_data instance has a nanoBragg beam object, which takes spectra and converts to nanoBragg xray_beams
        self.S.beam.spectrum = self.Modelers[self._i_shot].spectra
        self.D.xray_beams = self.S.beam.xray_beams

    def _get_panels_fasts_slows(self):
        pass

    def compute_functional_gradients_diag(self):
        self.compute_functional_and_gradients()
        return self._f, self._g, self.d

    def compute_functional_and_gradients(self):
        t = time.time()
        out = self._compute_functional_and_gradients()
        t = time.time()-t
        self.print("TOok %.4f sec to compute functional and grad" % t)
        return out

    def _compute_functional_and_gradients(self):
        LOGGER.info("BEGIN FUNC GRAD ; iteration %d" % self.iterations)
        #if self.verbose:
        #    self._print_iteration_header()

        self.target_functional = 0

        self.grad = flex.double(self.n_total_params)
        if self.calc_curvatures:
            self.curv = flex.double(self.n_total_params)

        LOGGER.info("start update Fcell")
        self._store_updated_Fcell()
        self._update_Fcell()  # update the structure factor with the new x
        LOGGER.info("done update Fcell")
        self._MPI_save_state_of_refiner()
        self._update_spectra_coefficients()  # updates the diffBragg lambda coefficients if refinining spectra

        tshots = time.time()

        LOGGER.info("Iterate over %d shots" % len(self.shot_ids))
        self._shot_Zscores = []
        save_model = self.save_model_freq is not None and self.iterations % self.save_model_freq == 0
        if save_model:
            self._save_model_dir = os.path.join(self.model_dir, "iter%d" % self.iterations)

            if COMM.rank == 0 and not os.path.exists(self._save_model_dir):
                os.makedirs(self._save_model_dir)
            COMM.barrier()

        for self._i_shot in self.shot_ids:
            self.scale_fac = self._get_spot_scale(self._i_shot)**2
            self.b_fac = self._get_bfactor(self._i_shot)

            # TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
            #LOGGER.info("update models shot %d " % self._i_shot)
            self._update_beams()
            self._update_umatrix()
            self._update_ucell()
            self._update_ncells()
            self._update_ncells_def()
            self._update_rotXYZ()
            self._update_eta()  # mosaic spread
            self._update_dxtbx_detector()
            self._update_sausages()

            self._run_diffBragg_current()

            # CHECK FOR SIGNAL INTERRUPT HERE
            if self.break_signal is not None:
                signal.signal(self.break_signal, self._sig_hand.handle)
                self._MPI_check_for_break_signal()

            # TODO pre-extractions for all parameters
            self._pre_extract_deriv_arrays()
            self._extract_pixel_data()
            self._evaluate_averageI()
            if self.poisson_only:
                self._evaluate_log_averageI()
            else:
                self._evaluate_log_averageI_plus_sigma_readout()

            self._derivative_convenience_factors()

            if self.iterations % self.saveZ_freq == 0:
                MOD = self.Modelers[self._i_shot]
                self._spot_Zscores = []
                for i_fcell in MOD.unique_i_fcell:
                    for slc in MOD.i_fcell_slices[i_fcell]:
                        sigZ = self._Zscore[slc]
                        trus = MOD.all_trusted[slc]
                        sigZ = sigZ[trus].std()
                        #sigZ = np.abs(sigZ[trus]).mean()
                        self._spot_Zscores.append((i_fcell, sigZ))
                self._shot_Zscores.append(self._spot_Zscores)

            if save_model:
                MOD = self.Modelers[self._i_shot]
                P = MOD.all_pid
                F = MOD.all_fast
                S = MOD.all_slow
                M = self.model_Lambda
                B = MOD.all_background
                D = MOD.all_data
                C = self.model_bragg_spots
                Z = self._Zscore
                iF = MOD.all_fcell_global_idx
                iROI = MOD.roi_id
                trust = MOD.all_trusted

                model_info = {"p": P, "f": F, "s": S, "model": M,
                        "background": B, "data": D, "bragg": C,
                        "Zscore": Z, "i_fcell": iF, "trust": trust,
                        "i_roi": iROI}
                self._save_model(model_info)
            self._shot_Zscores.append(self._spot_Zscores)
            self._is_trusted = self.Modelers[self._i_shot].all_trusted
            self.target_functional += self._target_accumulate()
            self._spot_scale_derivatives()
            #self._Bfactor_derivatives()
            self._Fcell_derivatives()
            self._shot_Zscores.append(self._spot_Zscores)
        tshots = time.time()-tshots
        LOGGER.info("Time rank worked on shots=%.4f" % tshots)
        tmpi = time.time()
        LOGGER.info("MPI aggregation of func and grad")
        self._mpi_aggregation()
        tmpi = time.time() - tmpi
        LOGGER.info("Time for MPIaggregation=%.4f" % tmpi)

        LOGGER.info("Aliases")
        self._f = self.target_functional
        self._g = self.g = self.grad
        self.d = self.curv
        LOGGER.info("curvature analysis")
        self._curvature_analysis()

        # reset ROI pixels TODO: is this necessary
        LOGGER.info("Zero pixels")
        self.D.raw_pixels_roi *= 0
        self.gnorm = -1

        tsave = time.time()
        LOGGER.info("DUMP param and Zscore data")
        self._save_Zscore_data()
        tsave = time.time()-tsave
        LOGGER.info("Time to dump param and Zscore data: %.4f" % tsave)

        self.iterations += 1
        self.f_vals.append(self.target_functional)
        #time.sleep(self.pause_after_iteration)

        if self.calc_curvatures and not self.use_curvatures:
            if self.num_positive_curvatures == self.use_curvatures_threshold:
                raise BreakToUseCurvatures

        LOGGER.info("DONE WITH FUNC GRAD")
        return self._f, self._g

    def _save_model(self, model_info):
        LOGGER.info("SAVING MODEL FOR SHOT %d" % self._i_shot)
        df = pandas.DataFrame(model_info)
        df["shot_id"] = self._i_shot
        outdir = self._save_model_dir
        outname = os.path.join(outdir, "rank%d_shot%d_ITER%d.pkl" % (COMM.rank, self._i_shot, self.iterations))
        df.to_pickle(outname)

    def _save_Zscore_data(self):
        if not self.iterations % self.saveZ_freq == 0:
            return
        outdir = os.path.join(self.Zdir, "rank%d_Zscore" % self.rank)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        fname = os.path.join(outdir, "sigZ_iter%d_rank%d" % (self.iterations, self.rank))
        np.save(fname, self._shot_Zscores)

    def _sanity_check_grad(self):
        pass

    def _Fcell_derivatives(self):
        if not self.refine_Fcell:
            return
        MOD = self.Modelers[self._i_shot]
        for i_fcell in MOD.unique_i_fcell:

            multi = self.hkl_frequency[i_fcell]
            if multi < self.min_multiplicity:
                continue

            xpos = self.fcell_xstart + i_fcell
            Famp = self._fcell_at_i_fcell[i_fcell]
            #resolution_id = self.res_group_id_from_fcell_index[i_fcell]
            sig = 1  # self.sigma_for_res_id[resolution_id] * self.fcell_sigma_scale
            for slc in MOD.i_fcell_slices[i_fcell]:
                self.fcell_dI_dtheta = self.fcell_deriv[slc]

                if self.log_fcells:
                    # case 2 rescaling
                    sig_times_fcell = sig*Famp
                    d = sig_times_fcell*self.fcell_dI_dtheta
                else:
                    # case 1 rescaling
                    d = sig*self.fcell_dI_dtheta

                gterm = self.common_grad_term[slc]
                g_accum = d*gterm
                trust = MOD.all_trusted[slc]
                self.grad[xpos] += (g_accum[trust].sum())*.5

    def _spot_scale_derivatives(self, return_derivatives=False):
        if not self.refine_crystal_scale:
            return
        S = np.sqrt(self.scale_fac)
        dI_dtheta = (2./S)*self.model_bragg_spots
        d2I_dtheta2 = (2./S/S)*self.model_bragg_spots
        # second derivative is 0 with respect to scale factor
        sig = self.Modelers[self._i_shot].PAR.Scale.sigma
        d = dI_dtheta*sig
        d2 = d2I_dtheta2 *(sig**2)

        xpos = self.spot_scale_xpos[self._i_shot]
        self.grad[xpos] += self._grad_accumulate(d)
        if self.calc_curvatures:
            self.curv[xpos] += self._curv_accumulate(d, d2)

        if return_derivatives:
            return d, d2

    def _Bfactor_derivatives(self):
        LOGGER.info("derivatives of Bfactors for shot %d: current B=%e Ang^2" % (self._i_shot, self.b_fac**2))
        if self.params.fix.B:
            return
        dI_dtheta = -.5*self.model_bragg_spots*self.Bfactor_qterm * self.b_fac
        d2I_dtheta2 = 0 #-.5*self.model_bragg_spots*self.Bfactor_qterm
        # second derivative is 0 with respect to scale factor
        sig = self.Modelers[self._i_shot].PAR.B.sigma
        d = dI_dtheta*sig
        d2 = d2I_dtheta2*(sig**2)
        #from IPython import embed
        #embed()

        xpos = self.Bfactor_xpos[self._i_shot]
        self.grad[xpos] += self._grad_accumulate(d)
        if self.calc_curvatures:
            self.curv[xpos] += self._curv_accumulate(d, d2)

    def _mpi_aggregation(self):
        # reduce the broadcast summed results:
        LOGGER.info("aggregate barrier")
        self._MPI_barrier()
        LOGGER.info("Functional")
        self.target_functional = self._MPI_reduce_broadcast(self.target_functional)
        LOGGER.info("gradients")
        self.grad = self._MPI_reduce_broadcast(self.grad)
        if self.calc_curvatures:
            self.curv = self._MPI_reduce_broadcast(self.curv)

    def _curvature_analysis(self):
        self.tot_neg_curv = 0
        self.neg_curv_shots = []
        if self.calc_curvatures:
            self.is_negative_curvature = self.curv.as_numpy_array() < 0
            self.tot_neg_curv = sum(self.is_negative_curvature)

        if self.calc_curvatures and not self.use_curvatures:
            if self.tot_neg_curv == 0:
                self.num_positive_curvatures += 1
                self.d = self.curv
                self._verify_diag()
            else:
                self.num_positive_curvatures = 0
                self.d = None

        if self.use_curvatures:
            assert self.tot_neg_curv == 0
            self.request_diag_once = False
            self.diag_mode = "always"  # TODO is this proper place to set ?
            self.d = self.curv
            self._verify_diag()
        else:
            self.d = None

    def _get_refinement_string_label(self):
        refine_str = "refining "
        if self.refine_Fcell:
            refine_str += "fcell, "
        if self.refine_ncells:
            refine_str += "Ncells, "
        if self.refine_ncells_def:
            refine_str += "Ncells_def, "
        if self.refine_Bmatrix:
            refine_str += "Bmat, "
        if self.refine_Umatrix:
            refine_str += "Umat, "
        if self.refine_crystal_scale:
            refine_str += "scale, "
        if self.refine_background_planes:
            refine_str += "bkgrnd, "
        if self.refine_detdist:
            refine_str += "detector_distance, "
        if self.refine_panelRotO:
            refine_str += "panelRotO, "
        if self.refine_panelRotF:
            refine_str += "panelRotF, "
        if self.refine_panelRotS:
            refine_str += "panelRotS, "
        if self.refine_panelXY:
            refine_str += "panelXY, "
        if self.refine_panelZ:
            refine_str += "panelZ, "
        if self.refine_lambda0:
            refine_str += "Lambda0 (offset), "
        if self.refine_lambda1:
            refine_str += "Lambda1 (scale), "
        if self.refine_per_spot_scale:
            refine_str += "Per-spot scales, "
        if self.refine_eta:
            refine_str += "Eta, "
        if self.refine_blueSausages:
            refine_str += "Mosaic texture, "
        return refine_str

    def _print_iteration_header(self):
        refine_str = self._get_refinement_string_label()
        border = "<><><><><><><><><><><><><><><><>"
        if self.use_curvatures:

            self.print(
                "%s%s%s%s\nTrial%d (%s): Compute functional and gradients Iter %d %s(Using Curvatures)%s\n%s%s%s%s"
                % (Bcolors.HEADER, border,border,border, self.trial_id + 1, refine_str, self.iterations + 1, Bcolors.OKGREEN, Bcolors.HEADER, border,border,border, Bcolors.ENDC))
        else:
            self.print("%s%s%s%s\n, Trial%d (%s): Compute functional and gradients Iter %d PosCurva %d\n%s%s%s%s"
                  % (Bcolors.HEADER, border, border, border, self.trial_id + 1, refine_str, self.iterations + 1, self.num_positive_curvatures, border, border,border, Bcolors.ENDC))

    def _MPI_save_state_of_refiner(self):
        if self.I_AM_ROOT and self.output_dir is not None and self.refine_Fcell:
            outf = os.path.join(self.output_dir, "_fcell_trial%d_iter%d" % (self.trial_id, self.iterations))
            np.savez(outf, fvals=self._fcell_at_i_fcell)

    def _show_plots(self, i_spot, n_spots):
        pass

    def _poisson_target(self):
        fterm = self.model_Lambda - self.Imeas * self.log_Lambda
        if self._is_trusted is not None:
            fterm = fterm[self._is_trusted]
        fterm = fterm.sum()
        return fterm

    def _poisson_d(self, d):
        gterm = d * self.one_minus_k_over_Lambda
        if self._is_trusted is not None:
            gterm = gterm[self._is_trusted]
        gterm = gterm.sum()
        return gterm

    def _poisson_d2(self, d, d2):
        cterm = d2 * self.one_minus_k_over_Lambda + d * d * self.k_over_squared_Lambda
        if self._is_trusted is not None:
            cterm = cterm[self._is_trusted]
        return cterm.sum()

    def _gaussian_target(self):
        fterm = self.log2pi + self.log_v + self.u*self.u*self.one_over_v
        if self._is_trusted is not None:
            fterm = fterm[self._is_trusted]
        fterm = 0.5*fterm.sum()
        return fterm

    def _gaussian_d(self, d):
        gterm = d * self.one_over_v * self.one_minus_2u_minus_u_squared_over_v
        if self._is_trusted is not None:
            gterm = gterm[self._is_trusted]
        gterm = 0.5*gterm.sum()
        return gterm

    def _gaussian_d2(self, d, d2):
        cterm = self.one_over_v * (d2*self.one_minus_2u_minus_u_squared_over_v -
                                   d*d*(self.one_over_v_times_one_minus_2u_minus_u_squared_over_v -
                                        (2 + 2*self.u_times_one_over_v + self.u_u_one_over_v*self.one_over_v)))
        if self._is_trusted is not None:
            cterm = cterm[self._is_trusted]
        cterm = .5 * (cterm.sum())
        return cterm

    def _derivative_convenience_factors(self):
        self.Imeas = self.Modelers[self._i_shot].all_data
        #one_over_Lambda = 1. / self.model_Lambda
        #self.one_minus_k_over_Lambda = (1. - self.Imeas * one_over_Lambda)
        #self.k_over_squared_Lambda = self.Imeas * one_over_Lambda * one_over_Lambda


        self.u = self.Imeas - self.model_Lambda
        self.one_over_v = 1. / (self.model_Lambda + self.sigma_r ** 2)
        self.one_minus_2u_minus_u_squared_over_v = 1 - 2 * self.u - self.u * self.u * self.one_over_v
        #self.u_times_one_over_v = self.u*self.one_over_v
        #self.u_u_one_over_v = self.u*self.u_times_one_over_v
        #self.one_over_v_times_one_minus_2u_minus_u_squared_over_v = self.one_over_v*self.one_minus_2u_minus_u_squared_over_v
        self.common_grad_term = self.one_over_v * self.one_minus_2u_minus_u_squared_over_v
        #if self.compute_Z:
        self._Zscore = self.u*np.sqrt(self.one_over_v)
        #self.one_over_v_data = 1. / (self.Imeas + self.sigma_r ** 2)
        #self._Zscore = self.u*np.sqrt(self.one_over_v_data)

    def _evaluate_log_averageI(self):  # for Poisson only stats
        try:
            self.log_Lambda = np.log(self.model_Lambda)
        except FloatingPointError:
            pass
        if any((self.model_Lambda <= 0).ravel()):
            is_bad = self.model_Lambda <= 0
            self.log_Lambda[is_bad] = 1e-6
            self.print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL (negative_models=%d)!!!!!!!!!\n<><><><><><><><><>\n" % self.num_negative_model)
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda[self.model_Lambda <= 0] = 0

    def _evaluate_log_averageI_plus_sigma_readout(self):
        v = self.model_Lambda + self.sigma_r ** 2
        v_is_neg = (v <= 0).ravel()
        if any(v_is_neg):
            self.print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_v = np.log(v)
        self.log_v[v <= 0] = 0  # but will I ever negative_model ?

    def get_refined_Bmatrix(self, i_shot):
        return self.Modelers[i_shot].PAR.ucell_man.B_recipspace

    def curvatures(self):
        return self.curv

    def _MPI_sync_hkl_freq(self):
            if self.refine_Fcell:
                if self.rank != 0:
                    self.hkl_frequency = None
                self.hkl_frequency = COMM.bcast(self.hkl_frequency)

    def _MPI_sync_fcell_parameters(self):
        if not self.I_AM_ROOT:
            self.sigma_for_res_id = None
            self.res_group_id_from_fcell_index = None
            self.resolution_ids_from_i_fcell = self.fcell_sigmas_from_i_fcell = self.fcell_init_from_i_fcell = None

        if self.rescale_params:
            if self.refine_Fcell:
                self.fcell_sigmas_from_i_fcell = COMM.bcast(self.fcell_sigmas_from_i_fcell)
                self.fcell_init_from_i_fcell = COMM.bcast(self.fcell_init_from_i_fcell)

    def _MPI_sync_panel_params(self):
        if not self.I_AM_ROOT:
            self.panelRot_params = None
            self.panelX_params = None
            self.panelY_params = None
            self.panelZ_params = None
        self.panelRot_params = COMM.bcast(self.panelRot_params)
        self.panelX_params = COMM.bcast(self.panelX_params)
        self.panelY_params = COMM.bcast(self.panelY_params)
        self.panelZ_params = COMM.bcast(self.panelZ_params)

    def _MPI_reduce_broadcast(self, var):
        var = COMM.reduce(var, MPI.SUM, root=0)
        var = COMM.bcast(var, root=0)
        return var

    def _MPI_barrier(self):
        COMM.barrier()
