from __future__ import absolute_import, division, print_function
from simtbx.diffBragg.utils import nearest_non_zero

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

import h5py

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
from copy import deepcopy
from cctbx import miller, sgtbx


class LocalRefiner(BaseRefiner):

    def __init__(self,
                 shot_modelers, sgsymbol):
        BaseRefiner.__init__(self)
        self.rank = COMM.rank

        self.Modelers = shot_modelers
        self.shot_ids = sorted(self.Modelers.keys())
        self.n_shots = len(shot_modelers)
        self.n_shots_total = COMM.bcast(COMM.reduce(self.n_shots))
        LOGGER.debug("Loaded %d shots across all ranks" % self.n_shots_total)
        self.init_image_corr = None
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
        self.model_Lambda = self.tilt_plane + self.model_bragg_spots

    def make_output_dir(self):
        if self.I_AM_ROOT and self.output_dir is not None and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def _dump_parameters_to_hdf5(self):
        if self.parameter_hdf5_path is not None and self.iterations % self.parameter_hdf5_write_freq == 0:
            with h5py.File(self.parameter_hdf5_path, 'w') as h5:
                for exp_name in self.parameters.keys:
                    h5.create_dataset("Ncells_abc/%s" % exp_name, data=self.parameters.Ncells_abc[exp_name])
                    h5.create_dataset("Ncells_def/%s" % exp_name, data=self.parameters.Ncells_def[exp_name])
                    h5.create_dataset("RotXYZ/%s" % exp_name, data=self.parameters.rotXYZ[exp_name])
                    h5.create_dataset("Bmat/%s" % exp_name, data=self.parameters.Bmatrix[exp_name])
                    h5.create_dataset("eta/%s" % exp_name, data=self.parameters.eta[exp_name])
                    h5.create_dataset("spot_scale/%s" % exp_name, data=self.parameters.spot_scale[exp_name])
                    h5.create_dataset("wavelen_offset/%s" % exp_name, data=self.parameters.wavelen_offset[exp_name])
                    h5.create_dataset("wavelen_scale/%s" % exp_name, data=self.parameters.wavelen_scale[exp_name])

                # add global parameters
                if self.I_AM_ROOT:
                    h5.create_dataset("panelX", data=self.parameters.panelX)
                    h5.create_dataset("panelY", data=self.parameters.panelY)
                    h5.create_dataset("panelZ", data=self.parameters.panelZ)
                    h5.create_dataset("panelO", data=self.parameters.panelO)
                    h5.create_dataset("panelF", data=self.parameters.panelF)
                    h5.create_dataset("panelS", data=self.parameters.panelS)
                    h5.create_dataset("panelOrig", data=self.parameters.panelOrig)
                    h5.create_dataset("panelFast", data=self.parameters.panelFast)
                    h5.create_dataset("panelSlow", data=self.parameters.panelSlow)
                #if self.refine_detdist:
                #    panZ = [self._get_detector_distance_val(i_shot)] * len(self.S.detector)

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
        N_PARAM_PER_SHOT = 1
        self.n_ucell_param = len(self.Modelers[test_shot].PAR.ucell_man.variables)
        #N_PARAM_PER_SHOT = len(self.Modelers[test_shot].PARLIST)
        self.n_total_params = self.n_total_shots*N_PARAM_PER_SHOT + self.n_global_fcell

        self.spot_scale_xpos = {}
        for shot_id in self.shot_ids:
            self.spot_scale_xpos[shot_id] = self.shot_mapping[shot_id]*N_PARAM_PER_SHOT

        LOGGER.info("--0 create an Fcell mapping")
        if self.refine_Fcell:
            idx, data = self.S.D.Fhkl_tuple
            self.idx_from_p1 = {h: i for i, h in enumerate(idx)}
            self._make_p1_equiv_mapping()
            # self.p1_from_idx = {i: h for i, h in zip(idx, data)}

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        self.x = flex.double(np.ones(self.n_total_params))
        LOGGER.info("--1 Setting up per shot parameters")

        self.fcell_xstart = self.n_total_shots*N_PARAM_PER_SHOT

        LOGGER.info("REduction of global data layout")
        self.hkl_totals = []
        if self.refine_Fcell:
            for i_shot in self.shot_ids:
                for i_h, h in enumerate(self.Modelers[i_shot].Hi_asu):
                    self.hkl_totals.append(self.idx_from_asu[h])
            self.hkl_totals = self._MPI_reduce_broadcast(self.hkl_totals)

        self._MPI_setup_global_params()
        self._MPI_sync_fcell_parameters()
        # reduce then broadcast fcell
        LOGGER.info("--3 combining parameters across ranks")
        self._MPI_sync_hkl_freq()

        if self.x_init is not None:
            print("Initializing with provided x_init array")
            self.x = self.x_init
        elif self.restart_file is not None:
            print("Restarting from parameter file %s" % self.restart_file)
            self.x = flex.double(np.load(self.restart_file)["x"])

        # setup the diffBragg instance
        self.D = self.S.D

        self.D.refine(self._fcell_id)
        self.D.initialize_managers()
        self._MPI_barrier()

    def _get_shot_mapping(self):

        all_shot_ids = COMM.gather(self.shot_ids)
        shot_mapping = None
        if COMM.rank == 0:
            unique_shot_ids = set([sid for shot_ids in all_shot_ids for sid in shot_ids])
            shot_mapping = {shot_id: i_shot for i_shot, shot_id in enumerate(unique_shot_ids)}
        shot_mapping = COMM.bcast(shot_mapping)
        return shot_mapping

    def _make_p1_equiv_mapping(self):

        self.p1_indices_from_i_fcell = {}
        for i_fcell in range(self.n_global_fcell):
            hkl_asu = self.asu_from_idx[i_fcell]

            self.p1_indices_from_i_fcell[i_fcell] = []
            equivs = [i.h() for i in miller.sym_equiv_indices(self.space_group, hkl_asu).indices()]
            for h_equiv in equivs:
                # get the nanoBragg p1 miller table index corresponding to this hkl equivalent
                try:
                    p1_idx = self.idx_from_p1[h_equiv]  # T
                except KeyError:
                    self.print("Whoops, missing index", h_equiv)
                    continue
                self.p1_indices_from_i_fcell[i_fcell].append(p1_idx)

    def _MPI_setup_global_params(self):
        if self.I_AM_ROOT:
            self.print("--2 Setting up global parameters")
            if self.output_dir is not None:
                # np.save(os.path.join(self.output_dir, "f_truth"), self.f_truth)  #FIXME by adding in the correct truth from Fref
                np.save(os.path.join(self.output_dir, "f_asu_map"), self.asu_from_idx)

            self._setup_fcell_params()

    def _setup_fcell_params(self):
        if self.refine_Fcell:
            self.print("----loading fcell data")
            # this is the number of observations of hkl (accessed like a dictionary via global_fcell_index)
            self.print("---- -- counting hkl totes")
            self.hkl_frequency = Counter(self.hkl_totals)
            np.save(os.path.join(self.output_dir, "f_asu_multi"), self.hkl_frequency)

            # initialize the Fhkl global values
            self.print("--- --- --- inserting the Fhkl array in the parameter array... ")
            asu_idx = [self.asu_from_idx[idx] for idx in range(self.n_global_fcell)]
            self._refinement_millers = flex.miller_index(tuple(asu_idx))
            Findices, Fdata = self.S.D.Fhkl_tuple
            vals = [Fdata[self.idx_from_p1[h]] for h in asu_idx]  # TODO am I correct/
            self.fcell_init = deepcopy(vals)  # store the initial values  for rescaling procedure

            if self.Fobs is not None:  # TODO should this ever be None ?
                miller_binner = self.Fobs.binner()
                miller_bin_idx = miller_binner.bin_indices()

                unique_bins = sorted(set(miller_bin_idx))
                sigmas = []
                for i_bin in unique_bins:
                    dmax, dmin = miller_binner.bin_d_range(i_bin)
                    f_selection = self.Fobs.resolution_filter(d_min=dmin, d_max=dmax)
                    sigma = np.sqrt(np.mean(f_selection.data().as_numpy_array() ** 2))
                    sigmas.append(sigma)  # sigma_for_res_id[i_bin] = sigma
                self.sigma_for_res_id = {}
                summed_sigma = 0
                for ii, sigma in enumerate(sigmas):
                    i_bin = unique_bins[ii]
                    if sigma == 0:
                        sigma = nearest_non_zero(sigmas, ii)
                    if sigma == 0:
                        bin_rng = miller_binner.bin_d_range(i_bin)
                        raise ValueError("sigma is being set to 0 for all fcell in range %.4f - %.4f" % bin_rng)
                    if self.rescale_fcell_by_resolution:
                        assert sigma > 0
                        self.sigma_for_res_id[i_bin] = 1. / sigma
                        summed_sigma += 1. / sigma
                    else:
                        self.sigma_for_res_id[i_bin] = 1.
                if self.rescale_fcell_by_resolution:
                    assert summed_sigma > 0
                    for ii in self.sigma_for_res_id.keys():
                        self.sigma_for_res_id[ii] = self.sigma_for_res_id[ii] / summed_sigma

                self.print("SIGMA FOR RES ID:")
                self.print(self.sigma_for_res_id)

                self.res_group_id_from_fcell_index = {}
                for ii, asu_index in enumerate(miller_binner.miller_indices()):
                    if asu_index not in self.idx_from_asu:
                        raise KeyError("something wrong Fobs does not contain the asu indices")
                    i_fcell = self.idx_from_asu[asu_index]
                    self.res_group_id_from_fcell_index[i_fcell] = miller_bin_idx[ii]

                self.resolution_ids_from_i_fcell = np.array([self.res_group_id_from_fcell_index[i_fcell] for i_fcell in range(self.n_global_fcell)])
                self.fcell_sigmas_from_i_fcell = np.array([ self.sigma_for_res_id[res_id]*self.fcell_sigma_scale for res_id in self.resolution_ids_from_i_fcell])
                self.fcell_init_from_i_fcell = np.array(self.fcell_init)

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

    def _get_bg_vals(self, i_shot, i_spot):
        pass

    def _send_ucell_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        pass

    def _run_diffBragg_current(self):
        LOGGER.info("run diffBragg for shot %d" % self._i_shot)
        pfs = self.Modelers[self._i_shot].pan_fast_slow
        nom_h = self.Modelers[self._i_shot].all_nominal_hkl
        self.D.add_diffBragg_spots(pfs, nom_h)
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
        idx, data = self.S.D.Fhkl_tuple
        data = data.as_numpy_array()
        for i_fcell in range(self.n_global_fcell):
            #new_Fcell_amplitude = self._get_fcell_val(i_fcell)
            new_Fcell_amplitude = self._fcell_at_i_fcell[i_fcell]

            # now surgically update the p1 array in nanoBragg with the new amplitudes
            # (need to update each symmetry equivalent)
            p1_indices = self.p1_indices_from_i_fcell[i_fcell]
            data[p1_indices] = new_Fcell_amplitude

        self.S.D.Fhkl_tuple = idx, flex.double(data)  # update nanoBragg again  # TODO: add flag to not re-allocate in nanoBragg!

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
        self._model_pix = self.D.raw_pixels_roi.as_numpy_array()[:npix]

        if self.refine_Fcell:
            dF = self.D.get_derivative_pixels(self._fcell_id)

            self._extracted_fcell_deriv = dF.set_selected(dF != dF, 0)
            self._extracted_fcell_deriv = self._extracted_fcell_deriv.as_numpy_array()[:npix]
            if self.calc_curvatures:
                d2F = self.D.get_second_derivative_pixels(self._fcell_id)
                self._extracted_fcell_second_deriv = d2F.set_selected(d2F != d2F, 0)
                self._extracted_fcell_second_deriv = self._extracted_fcell_second_deriv.as_numpy_array()[:npix]


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
            self.fcell_deriv = SG*(self._extracted_fcell_deriv[self.roi_sel])
            # handles Nan's when Fcell is 0 for whatever reason
            if self.calc_curvatures:
                self.fcell_second_deriv = SG*self._extracted_fcell_second_deriv[self.roi_sel]

    def _get_per_spot_scale(self, i_shot, i_spot):
        pass

    def _extract_pixel_data(self):
        self.model_bragg_spots = self.scale_fac*(self._model_pix[self.roi_sel])
        #self.model_bragg_spots *= self._get_per_spot_scale(self._i_shot, self._i_spot)
        #self._extract_Umatrix_derivative_pixels()
        #self._extract_sausage_derivs()
        #self._extract_Bmatrix_derivative_pixels()
        #self._extract_mosaic_parameter_m_derivative_pixels()
        #self._extract_ncells_def_derivative_pixels()
        #self._extract_detector_distance_derivative_pixels()
        self._extract_Fcell_derivative_pixels()
        #self._extract_spectra_coefficient_derivatives()
        #self._extract_panelRot_derivative_pixels()
        #self._extract_panelXYZ_derivative_pixels()

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
        save_model = self.iterations % self.save_model_freq == 0
        if save_model:
            # output buffers for panel, fast, slow, model, backgorund, data, braggSignal, Zscore, i_fcell
            P,F,S,M,B,D,C,Z,iF = [],[],[],[],[],[],[],[],[]
        for self._i_shot in self.shot_ids:

            self.scale_fac = self._get_spot_scale(self._i_shot)**2

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
            n_spots = len(self.Modelers[self._i_shot].rois)

            self._run_diffBragg_current()

            # CHECK FOR SIGNAL INTERRUPT HERE
            if self.break_signal is not None:
                signal.signal(self.break_signal, self._sig_hand.handle)
                self._MPI_check_for_break_signal()

            # TODO pre-extractions for all parameters
            self._pre_extract_deriv_arrays()
            self._spot_Zscores = []
            for i_spot in range(n_spots):
                self._i_spot = i_spot
                x1, x2, y1, y2 = self.Modelers[self._i_shot].rois[i_spot]
                if x2 - x1 == 0 or y2 - y1 == 0:
                    continue

                self._panel_id = int(self.Modelers[self._i_shot].pids[i_spot])
                self.mod = self.Modelers[self._i_shot]
                self.roi_sel = self.mod.roi_id == i_spot
                self.Imeas = self.mod.all_data[self.roi_sel]
                self._set_background_plane()
                self._extract_pixel_data()
                self._evaluate_averageI()


                if self.poisson_only:
                    self._evaluate_log_averageI()
                else:
                    self._evaluate_log_averageI_plus_sigma_readout()

                self._derivative_convenience_factors()

                if self.iterations % self.saveZ_freq == 0:
                    sigZ = (self._Zscore[self._is_trusted]).std()
                    miller_idx = self.mod.Hi_asu[i_spot]
                    i_fcell = self.idx_from_asu[miller_idx]
                    self._spot_Zscores.append((i_fcell, sigZ))

                if save_model:
                    Npix = len(self.model_Lambda)
                    P += [self._panel_id]*Npix
                    Y, X = np.indices((y2-y1, x2-x1))
                    F += list(X.ravel() + x1)
                    S += list(Y.ravel() + y1)
                    M += list(self.model_Lambda)
                    B += list(self.tilt_plane)
                    D += list(self.Imeas)
                    C += list(self.model_bragg_spots)
                    miller_idx = self.mod.Hi_asu[i_spot]
                    i_fcell = self.idx_from_asu[miller_idx]
                    Z += list(self._Zscore)
                    iF += [i_fcell] * Npix
                self.target_functional += self._target_accumulate()

                # accumulate the per pixel derivatives
                self._spot_scale_derivatives()
                self._Fcell_derivatives(i_spot)
                # Done with derivative accumulation
            if save_model:
                model_info = {"pan": P, "fast": F, "slow": S, "model": M, "background": B, "data": D, "bragg": C,
                              "Zscore": Z, "i_fcell": iF}
                self._save_model(model_info)
            self._shot_Zscores.append(self._spot_Zscores)
        # self.image_corr[self._i_shot] = self.image_corr[self._i_shot] / self.image_corr_norm[self._i_shot]
        tshots = time.time()-tshots
        LOGGER.info("Time rank worked on shots=%.4f" % tshots)
        self._append_global_parameters()
        tmpi = time.time()
        LOGGER.info("MPI aggregation of func and grad")
        self._mpi_aggregation()
        tmpi = tmpi-time.time()
        LOGGER.info("Time for MPIaggregation=%.4f" % tmpi)

        self._f = self.target_functional
        self._g = self.g = self.grad
        self.d = self.curv
        self._curvature_analysis()

        # reset ROI pixels TODO: is this necessary
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
        df = pandas.DataFrame(model_info)
        outdir = os.path.join(self.output_dir, "model_info")
        if not os.path.exists(outdir):
            if COMM.rank == 0:
                os.makedirs(outdir)
            COMM.barrier()
        outname = os.path.join(outdir, "rank%d_iter%d_shot%d.pkl" % (COMM.rank, self.iterations, self._i_shot))
        df.to_pickle(outname)

    def _save_Zscore_data(self):
        if not self.iterations % self.saveZ_freq == 0:
            return
        outdir = os.path.join(self.output_dir, "rank%d_Zscore" % self.rank)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        fname = os.path.join(outdir, "sigZ_iter%d_rank%d" % (self.iterations, self.rank))
        np.save(fname, self._shot_Zscores)

    def _sanity_check_grad(self):
        pass

    def _open_model_output_file(self):
        if not self.I_AM_ROOT:
            return
        if self.output_dir is None:
            self.print("Cannot save without an output dir, continuing without save")
            return
        self.save_time = 0
        from time import time
        t = time()
        from os.path import join as path_join
        from os.path import basename
        from h5py import File as h5_File
        dirpath = basename(self.FNAMES[self._i_shot]) + "-model"
        dirpath = path_join(self.output_dir, "models", dirpath)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        filepath = path_join(dirpath, "model_trial%d.h5" %(self.trial_id+1))
        self.print("Opening model output file %s" % filepath)
        self.model_hout = h5_File(filepath, "w")
        t = time() - t
        self.save_time += t

    def _close_model_output_file(self):
        if not self.I_AM_ROOT:
           return
        self.model_hout.close()
        self.print("Save time took %.4f seconds total" % self.save_time)

    def _Fcell_derivatives(self, i_spot):
        if not self.refine_Fcell:
            return
        # asu index
        miller_idx = self.Modelers[self._i_shot].Hi_asu[i_spot]
        # get multiplicity of this index
        multi = self.hkl_frequency[self.idx_from_asu[miller_idx]]
        # check if we are freezing this index during refinement
        # do the derivative
        if self.refine_Fcell and multi >= self.min_multiplicity:
            hkl_asu = self.Modelers[self._i_shot].Hi_asu[i_spot]
            i_fcell = self.idx_from_asu[hkl_asu]
            xpos = self.fcell_xstart + i_fcell

            self.fcell_dI_dtheta = self.fcell_deriv
            self.fcell_d2I_d2theta2 = self.fcell_second_deriv

            d = self.fcell_dI_dtheta
            d2 = self.fcell_d2I_d2theta2
            if self.rescale_params:
                #fcell = self._get_fcell_val(i_fcell)  # todo: interact with a vectorized object instead
                fcell = self._fcell_at_i_fcell[i_fcell]
                resolution_id = self.res_group_id_from_fcell_index[i_fcell]
                sig = self.sigma_for_res_id[resolution_id] * self.fcell_sigma_scale
                if self.log_fcells:
                    # case 2 rescaling
                    sig_times_fcell = sig*fcell
                    d = sig_times_fcell*self.fcell_dI_dtheta
                    if self.calc_curvatures:
                        d2 = (sig_times_fcell*sig_times_fcell)*self.fcell_d2I_d2theta2 + (sig*sig_times_fcell)*self.fcell_dI_dtheta
                else:
                    # case 1 rescaling
                    d = sig*self.fcell_dI_dtheta
                    if self.calc_curvatures:
                        d2 = (sig*sig)*self.fcell_d2I_d2theta2

            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)

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
        one_over_Lambda = 1. / self.model_Lambda
        self.one_minus_k_over_Lambda = (1. - self.Imeas * one_over_Lambda)
        self.k_over_squared_Lambda = self.Imeas * one_over_Lambda * one_over_Lambda

        self.u = self.Imeas - self.model_Lambda
        self.one_over_v = 1. / (self.model_Lambda + self.sigma_r ** 2)
        self.one_minus_2u_minus_u_squared_over_v = 1 - 2 * self.u - self.u * self.u * self.one_over_v
        self.u_times_one_over_v = self.u*self.one_over_v
        self.u_u_one_over_v = self.u*self.u_times_one_over_v
        self.one_over_v_times_one_minus_2u_minus_u_squared_over_v = self.one_over_v*self.one_minus_2u_minus_u_squared_over_v
        #if self.compute_Z:
        self._Zscore = self.u*np.sqrt(self.one_over_v)

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

    def _append_local_parameters(self):
        if self.parameter_hdf5_path is None:
            return

        exper_name = self.FNAMES[self._i_shot] if self.FNAMES is not None else "shot%d" % self._i_shot
        exper_name = exper_name.replace("/", "+FORWARD+SLASH+")

        crystal_scale = self._get_spot_scale(self._i_shot) ** 2 * self.D.spot_scale
        self.parameters.add_spot_scale(exper_name, crystal_scale)

        ncells_shot = self._get_m_val(self._i_shot)
        if len(ncells_shot) == 1:
            ncells_shot = [ncells_shot[0]]*3
        self.parameters.add_Ncells_abc(exper_name, ncells_shot)

        ncells_def_vals = self._get_ncells_def_vals(self._i_shot)
        self.parameters.add_Ncells_def(exper_name, ncells_def_vals)

        eta_vals = self._get_eta(i_shot=self._i_shot)
        self.parameters.add_eta(exper_name, eta_vals)
        lam01 = self._get_spectra_coefficients()
        if not lam01:
            lam0, lam1 = 0,1
        else:
            lam0, lam1 = lam01
        self.parameters.add_wavelen_offset(exper_name, lam0)
        self.parameters.add_wavelen_scale(exper_name, lam1)

        Bmat = self.get_refined_Bmatrix(self._i_shot)
        self.parameters.add_Bmatrix(exper_name, Bmat)

        rotX = self._get_rotX(self._i_shot)
        rotY = self._get_rotY(self._i_shot)
        rotZ = self._get_rotZ(self._i_shot)
        self.parameters.add_rotXYZ(exper_name, (rotX, rotY, rotZ))

    def _append_global_parameters(self):
        if self.parameter_hdf5_path is None:
            return
        opt_det = self.get_optimized_detector(i_shot=0)
        origs = []
        fasts = []
        slows = []
        for pid in range(len(opt_det)):
            orig = list(opt_det[pid].get_origin())
            fast = list(opt_det[pid].get_fast_axis())
            slow = list(opt_det[pid].get_slow_axis())
            origs.append(orig)
            fasts.append(fast)
            slows.append(slow)
        self.parameters.add_panelOrig(origs)
        self.parameters.add_panelFast(fasts)
        self.parameters.add_panelSlow(slows)
        panX, panY, panZ = zip(*[self._get_panelXYZ_val(pid) for pid in range(len(self.S.detector))])
        panO, panF, panS = zip(*[self._get_panelRot_val(pid) for pid in range(len(self.S.detector))])
        self.parameters.add_panelX(panX)
        self.parameters.add_panelY(panY)
        self.parameters.add_panelZ(panZ)
        self.parameters.add_panelO(panO)
        self.parameters.add_panelF(panF)
        self.parameters.add_panelS(panS)

    def get_refined_Bmatrix(self, i_shot):
        return self.Modelers[i_shot].PAR.ucell_man.B_recipspace

    def curvatures(self):
        return self.curv

    def _MPI_sync_hkl_freq(self):
        pass

    def _MPI_sync_fcell_parameters(self):
        pass

    def _MPI_sync_panel_params(self):
        pass

    def _data_for_write(self, parameter_dict):
        return [parameter_dict]

    def _MPI_aggregate_model_data_correlations(self):
        pass

    def _init_n_bad_shots(self):
        pass

    def _init_gather_ang_off(self):
        pass

    def _get_ang_off(self):
        pass

    def _MPI_reduce_broadcast(self, var):
        return var

    def _MPI_barrier(self):
        pass

    def _MPI_check_for_break_signal(self):
        pass
