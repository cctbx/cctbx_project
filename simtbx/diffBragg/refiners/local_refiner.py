from __future__ import absolute_import, division, print_function

import time
import warnings
import signal
import logging

LOGGER = logging.getLogger("main")
from collections import Iterable
warnings.filterwarnings("ignore")


class SignalHandler:
    """help to exit from HPC jobs before timeout, requires knowing exit signal, for summit its 12 or SIGUSR2 last I checked"""
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

from dxtbx.model import Panel, Detector
import os
import h5py

try:
    import pandas
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# TODO : consider PEP-8 ing these numpy imports, but do a NERSC massively MPI time-test first...
# for now, if it aint broke, dont fix it ...
from numpy import mean, median, unique, percentile, std
from simtbx.diffBragg.utils import makeMoffat_integPSF, convolve_with_psf
from scipy.stats import pearsonr
from numpy import vstack as np_vstack
from numpy import ascontiguousarray as np_ascontiguousarray
from numpy import log as np_log
from numpy.random import choice as np_random_choice
from numpy import sin as SIN
from numpy import indices as np_indices
from numpy import cos as COS
from numpy import arcsin as ASIN
from numpy import exp as np_exp
from numpy import load as np_load
from numpy import abs as ABS
from numpy import save as SAVE
from numpy import savez as SAVEZ
from numpy.linalg import norm
from numpy import ones_like as ONES_LIKE
from numpy import where as WHERE
from numpy import sqrt as SQRT
from numpy import array as ARRAY
from numpy import all as np_all
from numpy import pi as PI
from numpy import zeros as NP_ZEROS
from numpy import nan as NAN

from json import dump as JSON_DUMP
from os import makedirs as MAKEDIRS
from os.path import exists as EXISTS
from os.path import join as PATHJOIN

from simtbx.diffBragg.refiners.parameters import Parameters, RangedParameter
from simtbx.diffBragg.refiners import BreakBecauseSignal, BreakToUseCurvatures
from scitbx.array_family import flex
from cctbx.array_family import flex as cctbx_flex
flex_miller_index = cctbx_flex.miller_index
flex_double = flex.double
FLEX_BOOL = flex.bool
from scitbx.matrix import col
from simtbx.diffBragg.refiners import BaseRefiner
from scipy.optimize import minimize
from collections import Counter
import pylab as plt
from copy import deepcopy
from simtbx.diffBragg.utils import compare_with_ground_truth
from cctbx import miller, sgtbx
import itertools


class LocalRefiner(BaseRefiner):

    def __init__(self, n_total_params, n_local_params, local_idx_start,
                 shot_ucell_managers, shot_rois, shot_nanoBragg_rois,
                 shot_roi_imgs, shot_spectra, shot_crystal_GTs,
                 shot_crystal_models, shot_xrel, shot_yrel, shot_abc_inits, shot_asu=None,
                 global_param_idx_start=None,
                 shot_panel_ids=None,
                 log_of_init_crystal_scales=None,
                 all_crystal_scales=None, init_gain=1, perturb_fcell=False,
                 global_ncells=False, global_ucell=True, global_detector_distance=False,
                 shot_detector_distance_init=None,
                 sgsymbol="P43212",
                 shot_roi_darkRMS=None,
                 omega_kahn=None, selection_flags=None, shot_bg_coef=None, background_estimate=None,
                 verbose=False):
        """
        TODO: parameter x array boundaries should be done in this class, eliminating the need for local_idx_start
        TODO and global_idx_start
        TODO: ROI and nanoBragg ROI should be single variable

        :param n_total_params:  total number of parameters all shots
        :param n_local_params:  total number of parameters for each shot  (e.g. xtal rotation matrices)
        :param local_idx_start: for each shot, this specifies the starting point of its local parameters

        THE FOLLOWING ARE DICTIONARIES WHOSE KEYS ARE SHOT INDICES LOCAL TO THE RANK
        :param shot_ucell_managers: for each shot we have a unit cell manager , from crystal systems folder
        :param shot_rois: for each shot we have the list of ROIs [(x1,x2,y1,y2), .. ]
        :param shot_nanoBragg_rois:  for each shot we have list of ROIS in nanoBRagg format [[(x1,x2+1),(y1,y2+1)] , ...
        :param shot_roi_imgs: for each shot we have the data in each ROI
        :param shot_spectra: for each shot we have the spectra
        :param shot_crystal_GTs:  for each shot we have the ground truht crystal model (or None)
        :param shot_crystal_models: for each shot we have the crystal models
        :param shot_xrel: for each shot we have the relative fast scan coord of each pixel
        :param shot_yrel: for each shot we have relative slow scan coord of each pixel
        :param shot_abc_inits: for each shot we have the a,b,c values of the tilt plane
        :param shot_asu: for each shot we have the ASU miller index for each shoebox
        :param global_param_idx_start: This is position on x array where global parameters begin
        :param shot_panel_ids: for each shot we have panel IDs for each shoebox
        :param log_of_init_crystal_scales: for each shot we have crystal scale factors ( log)
        :param all_crystal_scales: for each shot we have ground truth scale factors (or None)
        :param init_gain: we have init gain correction estimate
        :param perturb_fcell: deprecated, leave as False
        :param global_ncells: do we refine ncells_abc per shot or global for all shots  (global_ncells=True)
        :param global_ucell: do we refine a unitcell per shot or global for all shots (global_ucell=True)
        :param shot_detector_distance_init: per shot origin Z
        :param global_detector_distance: do we refine a single detector Z position for all shots (default is True)
        :param omega_kahn: omega and kahn correction term for each panel
        """
        BaseRefiner.__init__(self)
        assert global_param_idx_start is not None
        assert shot_panel_ids is not None

        if not isinstance(local_idx_start, int):
            raise TypeError("local idx start should be an integer")

        self._sig_hand = SIGHAND
        self.rank = 0
        self.num_negative_model = 0
        self.global_ncells_param = global_ncells
        self.global_ucell_param = global_ucell
        assert not global_detector_distance, "deprecated"
        self.global_detector_distance_param = False
        self.special_flag = None
        self.debug = False
        self.num_Fcell_negative_model = 0
        self.perturb_fcell = perturb_fcell
        # dictionaries whose keys are the shot indices
        self.UCELL_MAN = shot_ucell_managers
        self.CRYSTAL_SCALE_TRUTH = all_crystal_scales  # ground truth of crystal scale factors..
        # cache shot ids and make sure they are identical in all other input dicts
        self.shot_ids = sorted(shot_ucell_managers.keys())
        self.big_dump = False
        self.show_watched = False
        self.n_shots = len(self.shot_ids)
        self.init_image_corr = None
        self.image_corr = {i_shot:None for i_shot in range(self.n_shots)}
        self.image_corr_norm = {i_shot:None for i_shot in range(self.n_shots)}
        self.shot_detector_distance_init = None
        if shot_detector_distance_init is not None:
            if self.verbose:
                self.print("check detector_distance")
            self.shot_detector_distance_init = self._check_keys(shot_detector_distance_init)
        # load the background coefficient  and extracted backgorund estimate
        self.shot_bg_coef = None
        if shot_bg_coef is not None:
            self.shot_bg_coef = self._check_keys(shot_bg_coef)
        self.background_estimate = background_estimate

        self.selection_flags = selection_flags
        if self.selection_flags is not None:
            if verbose:
                self.print("check selection flags")
            self.selection_flags = self._check_keys(self.selection_flags)

        # sanity check: no repeats of the same shot
        assert len(self.shot_ids) == len(set(self.shot_ids))

        self.OMEGA_KAHN = omega_kahn
        if verbose:
            self.print("check ROIS")
        self.ROIS = self._check_keys(shot_rois)
        if verbose:
            self.print("check ASU")
        self.ASU = None
        if shot_asu is not None:
            self.ASU = self._check_keys(shot_asu)
        if verbose:
            self.print("check nanoBRAGG ROIS")
        self.NANOBRAGG_ROIS = self._check_keys(shot_nanoBragg_rois)
        if verbose:
            self.print("check ROI IMAGES")
        self.ROI_IMGS = self._check_keys(shot_roi_imgs)
        self.ROI_DARK_RMS = None
        if shot_roi_darkRMS is not None:
            self.ROI_DARK_RMS = self._check_keys(shot_roi_darkRMS)
        if verbose:
            self.print("check SPECTRA")
        self.SPECTRA = self._check_keys(shot_spectra)
        self.CRYSTAL_GT = None
        if shot_crystal_GTs is not None:
            if verbose:
                self.print("check GT CRYSTAL")
            self.CRYSTAL_GT = self._check_keys(shot_crystal_GTs)
        if verbose:
            self.print("check CRYSTAL")
        self.CRYSTAL_MODELS = self._check_keys(shot_crystal_models)
        if verbose:
            self.print("check XREL YREL")

        self.XREL = self._check_keys(shot_xrel)
        self.YREL = self._check_keys(shot_yrel)
        if verbose:
            self.print("check ABC INIT")
        self.ABC_INIT = self._check_keys(shot_abc_inits)
        if verbose:
            self.print("check PANEL IDS")
        self.PANEL_IDS = self._check_keys(shot_panel_ids)

        # Total number of parameters in the MPI world
        self.n_total_params = n_total_params

        # total number of local parameters
        self.n_local_params = n_local_params


        # here are the indices of the local parameters in the global paramter arrays
        self.local_idx_start = local_idx_start

        self.calc_func = True  # NOTE: leave True, debug flag from older code
        self.multi_panel = True  # we are multi panel for all space and time back to the dawn of creation
        self.f_vals = []  # store the functional over time

        # start with the first shot
        self._i_shot = self.shot_ids[0]

        # These are the per-shot parameters
        self.n_rot_param = 3
        self.n_spot_scale_param = 1
        self.n_ucell_param = len(self.UCELL_MAN[self._i_shot].variables)
        self.n_ncells_param = 1

        self.n_per_shot_detector_distance_param = 1

        self.n_per_shot_ucell_param = self.n_ucell_param
        if global_ucell:
            self.n_per_shot_ucell_param = 0

        self.n_per_shot_ncells_param = 1
        if global_ncells:
            self.n_per_shot_ncells_param = 0

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

        if log_of_init_crystal_scales is None:
            log_of_init_crystal_scales = {s: 0 for s in self.shot_ids}
        else:
            assert sorted(log_of_init_crystal_scales.keys()) == self.shot_ids
        self.log_of_init_crystal_scales = log_of_init_crystal_scales

        self._init_gain = init_gain
        self.num_positive_curvatures = 0
        self._panel_id = None
        self.symbol = sgsymbol
        self.space_group = None
        if self.symbol is not None:
            self.space_group = sgtbx.space_group(sgtbx.space_group_info(symbol=self.symbol).type().hall_symbol())

        self.pid_from_idx = {}
        self.idx_from_pid = {}

        self.idx_from_asu = {}
        self.asu_from_idx = {}

        # For priors, use these ..  experimental
        # TODO: update this for general case (once it is actually working to begin with.. )
        self.ave_ucell = [78.95, 38.12]  ## Angstrom
        self.sig_ucell = [0.025, 0.025]
        self.sig_rot = 0.01  # radian

        # where the global parameters being , initially just gain and detector distance
        self.global_param_idx_start = global_param_idx_start

        self.a = self.b = self.c = None  # tilt plan place holder

        # optional properties
        self.FNAMES = None  # dict where index in shot index, and value is experiment name, optional
        self.PROC_FNAMES = None  # deprecated
        self.init_ang_off = None  # initial misortientations between model Xtals and reference Xtals
        self.current_ang_off = None  # TODO move me to pixel_refinement
        self.I_AM_ROOT = True

    def print(self, s, *args, **kwargs):
        """cheap logger"""
        if self.verbose:
            if isinstance(s, str):
                for line in s.split("\n"):
                    print(line, *args, **kwargs, end=self.print_end)
            else:
                print(s, *args, **kwargs, end=self.print_end)

    def setup_plots(self):
        if self.plot_images:
            if self.plot_residuals:
                from mpl_toolkits.mplot3d import Axes3D
                self.fig = plt.figure()
                self.ax = self.fig.gca(projection='3d')
                self.ax.set_yticklabels([])
                self.ax.set_xticklabels([])
                self.ax.set_zticklabels([])
                self.ax.set_zlabel("model residual")
                self.ax.set_facecolor("gray")
            else:
                self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1, ncols=2)
                self.ax1.imshow([[0, 1, 1], [0, 1, 2]])  # dummie plot
                self.ax2.imshow([[0, 1, 1], [0, 1, 2]])

    def __call__(self, *args, **kwargs):
        _, _ = self.compute_functional_and_gradients()
        return self.x, self._f, self._g, self.d

    @property
    def n(self):
        """LBFGS property"""
        return len(self.x)  # NOTEX

    @property
    def n_global_fcell(self):
        return len(self.idx_from_asu)

    @property
    def image_shape(self):
        panelXdim, panelYdim = self.S.detector[0].get_image_size()
        Npanels = len(self.S.detector)
        return Npanels, panelYdim, panelXdim

    @property
    def x_for_lbfgs(self):
        """LBFGS parameter array"""
        if self.only_pass_refined_x_to_lbfgs:
            return self.Xall.select(self.is_being_refined)
        else:
            return self.Xall

    @property
    def g_for_lbfgs(self):
        """LBFGS parameter array"""
        if self.only_pass_refined_x_to_lbfgs:
            return self.grad.select(self.is_being_refined)
        else:
            return self.grad

    @property
    def d_for_lbfgs(self):
        """LBFGS parameter array"""
        if self.only_pass_refined_x_to_lbfgs:
            return self.curv.select(self.is_being_refined)
        else:
            return self.curv

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
        self.model_Lambda = \
            self.gain_fac * self.gain_fac * (self.tilt_plane + self.model_bragg_spots)
        if self.refine_with_psf:
            self.model_Lambda = convolve_with_psf(self.model_Lambda, psf=self._psf, **self.psf_args)

    def _MPI_make_output_dir(self):
        if self.I_AM_ROOT and self.output_dir is not None and not EXISTS(self.output_dir):
            MAKEDIRS(self.output_dir)

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

        if self.parameters is None:
            self.parameters = Parameters()

        if self.I_AM_ROOT:
            self.print("Setup begins!")
        if self.refine_Fcell and not self.asu_from_idx:
            raise ValueError("Need to supply a non empty asu from idx map")
        if self.refine_Fcell and not self.idx_from_asu:  # # TODO just derive from its inverse
            raise ValueError("Need to supply a non empty idx from asu map")

        self.dummie_detector = deepcopy(self.S.detector)  # need to preserve original detector
        if self.refine_with_psf:
            fwhm_pix = self.psf_args["fwhm"] / self.psf_args["pixel_size"]
            kern_size = self.psf_args["psf_radius"]*2 + 1
            if self.I_AM_ROOT:
                self.print("USING PSF: %f fwhm_pixel and %dx%d kernel size" % (fwhm_pix, kern_size, kern_size))
            self._psf = makeMoffat_integPSF(fwhm_pix, kern_size, kern_size)

        self._MPI_make_output_dir()

        #if self.refine_panelXY or self.refine_panelRotS or self.refine_panelRotF or self.refine_panelRotO:
        if self.panel_group_from_id is not None:
            self.n_panel_groups = len(set(self.panel_group_from_id.values()))
        else:
            self.n_panel_groups = 1

        # get the Fhkl information from P1 array internal to nanoBragg
        if self.I_AM_ROOT:
            self.print("--0 create an Fcell mapping")
        if self.refine_Fcell:
            idx, data = self.S.D.Fhkl_tuple
            self.idx_from_p1 = {h: i for i, h in enumerate(idx)}
            self._make_p1_equiv_mapping()
            # self.p1_from_idx = {i: h for i, h in zip(idx, data)}

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        # determine total number of parameters
        # XYZ per panel group
        #n_global_params = 3*self.n_panel_groups
        # Rotation OFS per panel group
        #n_global_params += 3*self.n_panel_groups

        # Make the global sized parameter array, though here we only update the local portion
        self.Xall = flex_double(self.n_total_params)
        assert isinstance(self.n_total_params, int)
        self.is_being_refined = FLEX_BOOL(self.n_total_params, False)

        # store the starting positions in the parameter array for this shot
        self.rotX_xpos = {}
        self.rotY_xpos = {}
        self.rotZ_xpos = {}
        self.ucell_xstart = {}
        self.ncells_xstart = {}
        self.ncells_def_xstart = {}
        self.detector_distance_xpos = {}
        self.spot_scale_xpos = {}
        self.eta_xstart = {}
        self.eta_params = {}
        self.n_panels = {}
        self.bg_a_xstart = {}
        self.bg_b_xstart = {}
        self.bg_c_xstart = {}
        self.bg_coef_xpos = {}
        self.ucell_params = {}
        self.Ncells_params = {}
        self.Ncells_def_params = {}
        self.per_spot_scale_xpos = {}
        self.sausages_xpos = {}

        self.panelX_params = [RangedParameter()] * self.n_panel_groups
        self.panelY_params = [RangedParameter()] * self.n_panel_groups
        self.panelZ_params = [RangedParameter()] * self.n_panel_groups  # per shot offset-Z to each panel
        self.panelRot_params = [[RangedParameter(), RangedParameter(), RangedParameter()]] * self.n_panel_groups
        self.detector_distance_params = {}  # per shot offset to all panels
        if self.I_AM_ROOT:
            self.print("--1 Setting up per shot parameters")

        # this is a sliding parameter that points to the latest local (per-shot) parameter in the x-array
        _local_pos = self.local_idx_start

        self.panels_fasts_slows = {i_shot: None for i_shot in self.shot_ids}
        for i_shot in self.shot_ids:
            self.pid_from_idx[i_shot] = {i: pid for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.idx_from_pid[i_shot] = {pid: i for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.n_panels[i_shot] = len(self.pid_from_idx[i_shot])

            if self.bg_extracted:
                self.bg_coef_xpos[i_shot] = _local_pos
                self.Xall[self.bg_coef_xpos[i_shot]] = 1
                if not self.rescale_params:
                    raise NotImplementedError("bg coef mode only works in rescale mode")
                if self.refine_background_planes:
                    self.is_being_refined[self.bg_coef_xpos[i_shot]] = True
            else:
                n_spots = len(self.NANOBRAGG_ROIS[i_shot])
                self.bg_a_xstart[i_shot] = []
                self.bg_b_xstart[i_shot] = []
                self.bg_c_xstart[i_shot] = []
                _spot_start = _local_pos

                for i_spot in range(n_spots):
                    self._i_spot = i_spot
                    self.bg_a_xstart[i_shot].append(_spot_start)
                    self.bg_b_xstart[i_shot].append(self.bg_a_xstart[i_shot][i_spot] + 1)
                    self.bg_c_xstart[i_shot].append(self.bg_b_xstart[i_shot][i_spot] + 1)

                    a, b, c = self.ABC_INIT[i_shot][i_spot]
                    if self.bg_offset_only and self.bg_offset_positive:
                        if c < 0:
                            c = np_log(1e-9)
                        else:
                            c = np_log(c)
                    if self.rescale_params:
                        self.Xall[self.bg_a_xstart[i_shot][i_spot]] = 1
                        self.Xall[self.bg_b_xstart[i_shot][i_spot]] = 1
                        self.Xall[self.bg_c_xstart[i_shot][i_spot]] = 1

                    else:
                        self.Xall[self.bg_a_xstart[i_shot][i_spot]] = float(a)
                        self.Xall[self.bg_b_xstart[i_shot][i_spot]] = float(b)
                        self.Xall[self.bg_c_xstart[i_shot][i_spot]] = float(c)

                    _spot_start += 3
                    if self.refine_background_planes:
                        self.is_being_refined[self.bg_c_xstart[i_shot][i_spot]] = True
                        if not self.bg_offset_only:
                            self.is_being_refined[self.bg_a_xstart[i_shot][i_spot]] = True
                            self.is_being_refined[self.bg_b_xstart[i_shot][i_spot]] = True

            if self.bg_extracted:
                self.rotX_xpos[i_shot] = self.bg_coef_xpos[i_shot] + 1
            else:
                self.rotX_xpos[i_shot] = self.bg_c_xstart[i_shot][-1] + 1
            self.rotY_xpos[i_shot] = self.rotX_xpos[i_shot] + 1
            self.rotZ_xpos[i_shot] = self.rotY_xpos[i_shot] + 1

            if self.rescale_params:
                self.Xall[self.rotX_xpos[i_shot]] = 1
                self.Xall[self.rotY_xpos[i_shot]] = 1
                self.Xall[self.rotZ_xpos[i_shot]] = 1
            else:
                self.Xall[self.rotX_xpos[i_shot]] = 0
                self.Xall[self.rotY_xpos[i_shot]] = 0
                self.Xall[self.rotZ_xpos[i_shot]] = 0
            if self.refine_Umatrix:
                if self.refine_rotX:
                    self.is_being_refined[self.rotX_xpos[i_shot]] = True
                if self.refine_rotY:
                    self.is_being_refined[self.rotY_xpos[i_shot]] = True
                if self.refine_rotZ:
                    self.is_being_refined[self.rotZ_xpos[i_shot]] = True

            # continue adding local per shot parameters after rotZ_xpos
            _local_pos = self.rotZ_xpos[i_shot] + 1

            # global always starts here, we have to decide whether to put ncells / unit cell/ detector_distance parameters in global array
            _global_pos = self.global_param_idx_start

            if self.global_ucell_param:
                self.ucell_xstart[i_shot] = _global_pos
                _global_pos += self.n_ucell_param
            else:
                self.ucell_xstart[i_shot] = _local_pos
                _local_pos += self.n_ucell_param
                self.ucell_params[i_shot] = []
                for i_cell in range(self.n_ucell_param):
                    if self.rescale_params:
                        self.Xall[self.ucell_xstart[i_shot] + i_cell] = 1  #self.UCELL_MAN[i_shot].variables[i_cell]
                        uc_param = RangedParameter()
                        if self.use_ucell_ranges:
                            uc_param.init = self.ucell_inits[i_shot][i_cell]
                            uc_param.sigma = self.ucell_sigmas[i_cell]
                            uc_param.maxval = self.ucell_maxs[i_shot][i_cell]
                            uc_param.minval = self.ucell_mins[i_shot][i_cell]
                            self.ucell_params[i_shot].append(uc_param)
                    else:
                        self.Xall[self.ucell_xstart[i_shot] + i_cell] = self.UCELL_MAN[i_shot].variables[i_cell]

            # set refinement flags
            if self.refine_Bmatrix:
                for i_cell in range(self.n_ucell_param):
                    self.is_being_refined[self.ucell_xstart[i_shot] + i_cell] = True

            if self.global_ncells_param:
                self.ncells_xstart[i_shot] = _global_pos
                _global_pos += self.n_ncells_param
            else:
                self.ncells_xstart[i_shot] = _local_pos
                _local_pos += self.n_ncells_param
                self.Ncells_params[i_shot] = []
                for i_ncells in range(self.n_ncells_param):
                    ncells_xval = np_log(self.S.crystal.Ncells_abc[i_ncells]-3)
                    # TODO: each shot gets own starting Ncells
                    if self.use_Ncells_range or self.rescale_params:
                        self.Xall[self.ncells_xstart[i_shot] + i_ncells] = 1
                    else:
                        self.Xall[self.ncells_xstart[i_shot] + i_ncells] = ncells_xval
                    if self.use_Ncells_range:
                        ncells_param = RangedParameter()
                        ncells_param.init = self.m_init[i_shot][i_ncells]
                        ncells_param.sigma = self.m_sigma[i_ncells]
                        if len(self.m_range) == 2:
                            ncells_param.minval = self.m_range[0]
                            ncells_param.maxval = self.m_range[1]
                        elif len(self.m_range) == 6:
                            ncells_param.minval = self.m_range[i_ncells*2]
                            ncells_param.maxval = self.m_range[i_ncells*2 + 1]

                        self.Ncells_params[i_shot].append( ncells_param)
            # set refinement flags
            if self.refine_ncells:
                for i_ncells in range(self.n_ncells_param):
                    self.is_being_refined[self.ncells_xstart[i_shot] + i_ncells] = True

            self.ncells_def_xstart[i_shot] = _local_pos
            _local_pos += 3
            self.Ncells_def_params[i_shot] = []
            for i_ncells in range(3):
                self.Xall[self.ncells_def_xstart[i_shot] + i_ncells] = 1
                ncells_param = RangedParameter()
                ncells_param.init = self.ncells_def_init[i_shot][i_ncells]
                ncells_param.sigma = self.ncells_def_sigma[i_ncells]
                ncells_param.minval = self.ncells_def_range[i_ncells*2]
                ncells_param.maxval = self.ncells_def_range[i_ncells*2 + 1]
                self.Ncells_def_params[i_shot].append(ncells_param)
            if self.refine_ncells_def:
                for i_ncells in range(3):
                    self.is_being_refined[self.ncells_def_xstart[i_shot] + i_ncells] = True

            self.detector_distance_xpos[i_shot] = _local_pos
            detdist_param = RangedParameter()
            detdist_param.init = self.shot_detector_distance_init[i_shot]  # initial offset
            detdist_param.sigma = self.detector_distance_sigma
            detdist_param.minval = self.detector_distance_range[0]
            detdist_param.maxval = self.detector_distance_range[1]
            self.detector_distance_params[i_shot] = detdist_param
            _local_pos += 1
            self.Xall[self.detector_distance_xpos[i_shot]] = 1
            # set refinement flag
            if self.refine_detdist:
                assert self.rescale_params
                self.is_being_refined[self.detector_distance_xpos[i_shot]] = True

            self.spot_scale_xpos[i_shot] = _local_pos
            _local_pos += 1
            if self.rescale_params:
                self.Xall[self.spot_scale_xpos[i_shot]] = 1
            else:
                self.Xall[self.spot_scale_xpos[i_shot]] = self.log_of_init_crystal_scales[i_shot]
            if self.refine_crystal_scale:
                self.is_being_refined[self.spot_scale_xpos[i_shot]] = True

            self.eta_xstart[i_shot] = _local_pos
            self.eta_params[i_shot] = []
            for i_eta in range(3):
                xpos = self.eta_xstart[i_shot] + i_eta
                self.Xall[xpos] = 1
                eta_param = RangedParameter()
                eta_param.init = self.eta_init[i_shot][i_eta]
                eta_param.sigma = self.eta_sigma
                eta_param.minval = self.eta_min
                eta_param.maxval = self.eta_max
                self.eta_params[i_shot].append(eta_param)
                if self.refine_eta:
                    assert self.rescale_params
                    self.is_being_refined[xpos] = True
            _local_pos += 3

            #if self.refine_per_spot_scale:
            n_spots = len(self.NANOBRAGG_ROIS[i_shot])
            self.per_spot_scale_xpos[i_shot] = list(range(_local_pos, _local_pos + n_spots))
            for x in range(_local_pos, _local_pos + n_spots):
                self.Xall[x] = 1
            _local_pos += n_spots

            if self.refine_blueSausages:
                assert self.rescale_params
            self.sausages_xpos[i_shot] = []
            for i_sausage in range(self.num_sausages):
                for i_sausage_param in range(4):
                    #idx = i_sausage_param + 4 * i_sausage
                    self.sausages_xpos[i_shot].append(_local_pos)
                    self.Xall[_local_pos] = 1+i_sausage #idx #*0.1 #i_sausage
                    self.is_being_refined[_local_pos] = True
                    _local_pos += 1

        self.fcell_xstart = _global_pos

        self.spectra_coef_xstart = self.fcell_xstart + self.n_global_fcell
        #if self.refine_lambda0 or self.refine_lambda1 and self.n_total_params == 2
        #    self.spectra_coef_xstart = _local_pos

        self.panelRot_xstart = self.spectra_coef_xstart + self.n_spectra_param

        self.panelXY_xstart = self.panelRot_xstart + 3*self.n_panel_groups
        self.panelZ_xstart = self.panelXY_xstart + 2*self.n_panel_groups

        if self.refine_gain_fac:
            self.gain_xpos = self.n_total_params - 1

        # tally up HKL multiplicity
        if self.I_AM_ROOT:
            self.print("REduction of global data layout")
        self.hkl_totals = []
        fname_totals = []
        panel_id_totals = []
        # img_totals = []
        if self.refine_Fcell:
            for i_shot in self.ASU:
                for i_h, h in enumerate(self.ASU[i_shot]):
                    if self.FNAMES is not None:
                        fname_totals.append(self.FNAMES[i_shot])
                    panel_id_totals.append(self.PANEL_IDS[i_shot][i_h])
                    self.hkl_totals.append(self.idx_from_asu[h])
                    # img_totals.append(self.ROI_IMGS[i_shot][i_h])
            self.hkl_totals = self._MPI_reduce_broadcast(self.hkl_totals)

        import numpy as np

        self._MPI_setup_global_params()

        self._MPI_sync_fcell_parameters()

        self._MPI_sync_panel_params()

        # reduce then broadcast fcell
        if self.I_AM_ROOT == 0:
            self.print("--3 combining parameters across ranks")

        self.Xall = self._MPI_reduce_broadcast(self.Xall)

        # flex bool has no + operator so we convert to numpy
        self.is_being_refined = self._MPI_reduce_broadcast(self.is_being_refined.as_numpy_array())
        assert np.all([isinstance(x, np.bool_) for x in self.is_being_refined]), "problem with is being refined5 %d" % self.rank
        self.is_being_refined = FLEX_BOOL(self.is_being_refined)

        #self._MPI_barrier()
        # set the BFGS parameter array
        self.x = self.x_for_lbfgs

        # make the mapping from x to Xall
        refine_pos = WHERE(self.is_being_refined.as_numpy_array())[0]
        self.x2xall = {xi: xalli for xi, xalli in enumerate(refine_pos)}
        self.xall2x = {xalli: xi for xi, xalli in enumerate(refine_pos)}

        self._MPI_sync_hkl_freq()

        # See if restarting from save state

        if self.x_init is not None: #NOTEX
            print("Initializing with provided x_init array")
            self.Xall = self.x_init
            self.x = self.x_for_lbfgs
        elif self.restart_file is not None:
            print("Restarting from parameter file %s" % self.restart_file)
            self.Xall = flex_double(np_load(self.restart_file)["x"])
            self.x = self.x_for_lbfgs

        if self.I_AM_ROOT:
            self.print("--4 print initial stats")
        rotx, roty, rotz, uc_vals, ncells_vals, scale_vals, _, origZ = self._unpack_internal(self.Xall, lst_is_x=True)
        if self.I_AM_ROOT and self.big_dump and HAS_PANDAS:

            master_data = {"a": uc_vals[0], "c": uc_vals[1],
                           "Ncells": ncells_vals,
                           "scale": scale_vals,
                           "rotx": rotx,
                           "roty": roty,
                           "rotz": rotz,
                           "origZ": origZ}
            master_data = pandas.DataFrame(master_data)
            master_data["gain"] = 1 #self.Xall[self.gain_xpos]
            self.print(master_data.to_string())

        # make the parameter masks for isolating parameters of different types
        self._make_parameter_type_selection_arrays()

        if self.output_dir is not None:
            self._make_x_identifier_array()

        #self._setup_resolution_binner()
        # setup the diffBragg instance
        self.D = self.S.D

        if self.refine_Umatrix:
            if self.refine_rotX:
                self.D.refine(0)  # rotX
            #else:
            #    self.D.fix(0)
            if self.refine_rotY:
                self.D.refine(1)  # rotY
            #else:
            #    self.D.fix(1)
            if self.refine_rotZ:
                self.D.refine(2)  # rotZ
            #else:
            #    self.D.fix(2)
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.D.refine(i + 3)  # unit cell params
            #else:
            #    self.D.fix(i+3)
        if self.refine_ncells:
            self.D.refine(self._ncells_id)
        #else:
        #    self.D.fix(self._ncells_id)
        if self.refine_ncells_def:
            self.D.refine(self._ncells_def_id)
        #else:
        #    self.D.fix(self._ncells_def_id)
        if self.refine_detdist or self.refine_panelZ:
            self.D.refine(self._detector_distance_id)
        #else:
        #    self.D.fix(self._detector_distance_id)
        if self.refine_panelRotO:
            self.D.refine(self._panelRotO_id)
        #else:
        #    self.D.fix(self._panelRotO_id)
        if self.refine_panelRotF:
            self.D.refine(self._panelRotF_id)
        #else:
        #    self.D.fix(self._panelRotF_id)
        if self.refine_panelRotS:
            self.D.refine(self._panelRotS_id)
        #else:
        #    self.D.fix(self._panelRotS_id)
        if self.refine_panelXY:
            self.D.refine(self._panelX_id)
            self.D.refine(self._panelY_id)
        #else:
        #    self.D.fix(self._panelX_id)
        #    self.D.fix(self._panelY_id)
        if self.refine_Fcell:
            self.D.refine(self._fcell_id)
        #else:
        #    self.D.fix(self._fcell_id)
        if self.refine_lambda0:
            self.D.refine(self._lambda0_id)
        #else:
        #    self.D.fix(self._lambda0_id)
        if self.refine_lambda1:
            self.D.refine(self._lambda1_id)
        #else:
        #    self.D.fix(self._lambda1_id)
        if self.refine_eta:
            self.D.refine(self._eta_id)
        #else:
        #    self.D.fix(self._eta_id)
        if self.refine_blueSausages:
            self.D.refine(self._sausage_id)
        #else:
        #    self.D.fix(self._sausage_id)
        self.D.initialize_managers()

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
            # put in estimates for origin vectors
            # TODO: refine at the different hierarchy
            # get te first Z coordinate for now..
            # self.print("Setting origin: %f " % self.S.detector[0].get_local_origin()[2])
            if self.global_ucell_param:
                # TODO have parameter for global init of unit cell, right now its handled in the global_bboxes scripts
                for i_cell in range(self.n_ucell_param):
                    if self.rescale_params:
                        self.Xall[self.ucell_xstart[0] + i_cell] = 1
                    else:
                        self.Xall[self.ucell_xstart[0] + i_cell] = self.UCELL_MAN[0].variables[i_cell]

            if self.global_ncells_param:
                for i_ncells in range(self.n_ncells_param):
                    if self.rescale_params:
                        ncells_xval = 1
                    else:
                        ncells_xval = np_log(self.S.crystal.Ncells_abc[i_ncells] - 3)
                    self.Xall[self.ncells_xstart[0] + i_ncells] = ncells_xval

            # if self.refine_lambda0 or self.refine_lambda1:
            lambda_is_refined = self.refine_lambda0, self.refine_lambda1
            for i_spec_coef in range(self.n_spectra_param):
                xpos = self.spectra_coef_xstart + i_spec_coef
                self.is_being_refined[xpos] = lambda_is_refined[i_spec_coef]
                self.Xall[xpos] = 1
                if not self.rescale_params:
                    raise NotImplementedError("Cant refine spectra without rescale_params=True")

            refine_panelRot = self.refine_panelRotO, self.refine_panelRotF, self.refine_panelRotS
            for i_pan_group in range(self.n_panel_groups):
                for i_pan_rot in range(3):
                    xpos = self.panelRot_xstart + 3*i_pan_group + i_pan_rot

                    self.panelRot_params[i_pan_group][i_pan_rot].init = 0  # self.panelX_init[i_pan_group]
                    self.panelRot_params[i_pan_group][i_pan_rot].sigma = self.panelRot_sigma[i_pan_rot]
                    self.panelRot_params[i_pan_group][i_pan_rot].minval = self.panelRot_range[i_pan_rot][0]
                    self.panelRot_params[i_pan_group][i_pan_rot].maxval = self.panelRot_range[i_pan_rot][1]

                    if refine_panelRot[i_pan_rot]:
                        assert self.rescale_params
                        self.is_being_refined[xpos] = i_pan_group in self.panel_groups_being_refined
                    self.Xall[xpos] = 1

                xpos_X = self.panelXY_xstart + 2*i_pan_group
                self.Xall[xpos_X] = 1
                if self.refine_panelXY:
                    assert self.rescale_params
                    self.is_being_refined[xpos_X] = i_pan_group in self.panel_groups_being_refined #True
                self.panelX_params[i_pan_group].init = 0  # self.panelX_init[i_pan_group]
                self.panelX_params[i_pan_group].sigma = self.panelX_sigma
                self.panelX_params[i_pan_group].minval = self.panelX_range[0]
                self.panelX_params[i_pan_group].maxval = self.panelX_range[1]

                xpos_Y = xpos_X + 1
                self.Xall[xpos_Y] = 1
                if self.refine_panelXY:
                    self.is_being_refined[xpos_Y] = i_pan_group in self.panel_groups_being_refined #True
                self.panelY_params[i_pan_group].init = 0   # self.panelY_init[i_pan_group]
                self.panelY_params[i_pan_group].sigma = self.panelY_sigma
                self.panelY_params[i_pan_group].minval = self.panelY_range[0]
                self.panelY_params[i_pan_group].maxval = self.panelY_range[1]

                xpos_Z = self.panelZ_xstart + i_pan_group
                self.Xall[xpos_Z] = 1
                if self.refine_panelZ:
                    self.is_being_refined[xpos_Z] = i_pan_group in self.panel_groups_being_refined #True
                self.panelZ_params[i_pan_group].init = 0  # self.panelX_init[i_pan_group]
                self.panelZ_params[i_pan_group].sigma = self.panelZ_sigma
                self.panelZ_params[i_pan_group].minval = self.panelZ_range[0]
                self.panelZ_params[i_pan_group].maxval = self.panelZ_range[1]

            if self.output_dir is not None:
                # np.save(os.path.join(self.output_dir, "f_truth"), self.f_truth)  #FIXME by adding in the correct truth from Fref
                SAVE(os.path.join(self.output_dir, "f_asu_map"), self.asu_from_idx)

            # set gain TODO: implement gain dependent statistical model ? Per panel or per gain mode dependent ?
            #self.Xall[self.gain_xpos] = self._init_gain  # gain factor
            self._setup_fcell_params()

    def _setup_fcell_params(self):
        if self.refine_Fcell:
            self.print("----loading fcell data")
            # this is the number of observations of hkl (accessed like a dictionary via global_fcell_index)
            self.print("---- -- counting hkl totes")
            self.hkl_frequency = Counter(self.hkl_totals)
            SAVE(os.path.join(self.output_dir, "f_asu_multi"), self.hkl_frequency)

            # initialize the Fhkl global values
            self.print("--- --- --- inserting the Fhkl array in the parameter array... ")
            asu_idx = [self.asu_from_idx[idx] for idx in range(self.n_global_fcell)]
            self._refinement_millers = flex_miller_index(tuple(asu_idx))
            Findices, Fdata = self.S.D.Fhkl_tuple
            vals = [Fdata[self.idx_from_p1[h]] for h in asu_idx]  # TODO am I correct/
            if self.rescale_params:
                self.fcell_init = deepcopy(vals)  # store the initial values  for rescaling procedure
            if not self.rescale_params and self.log_fcells:
                vals = np_log(vals)
            for i_fcell in range(self.n_global_fcell):
                if self.rescale_params:
                    self.Xall[self.fcell_xstart + i_fcell] = 1
                else:
                    self.Xall[self.fcell_xstart + i_fcell] = vals[i_fcell]
                if self.refine_Fcell:  # TODO only refine if fcell is in the res range
                    self.is_being_refined[self.fcell_xstart + i_fcell] = True

            self.Fref_aligned = self.Fref
            if self.Fref is not None:
                self.Fref_aligned = self.Fref.select_indices(self.Fobs.indices())
                self.init_R1 = self.Fobs_Fref_Rfactor(use_binning=False, auto_scale=self.scale_r1)
                self.print("Initial R1 = %.4f" % self.init_R1)
            else:
                self.init_R1 = -1

            if self.Fobs is not None:  # TODO should this ever be None ?
                miller_binner = self.Fobs.binner()
                miller_bin_idx = miller_binner.bin_indices()

                import numpy as np  # TODO move me to top
                from simtbx.diffBragg.utils import nearest_non_zero

                unique_bins = sorted(set(miller_bin_idx))
                sigmas = []
                for i_bin in unique_bins:
                    dmax, dmin = miller_binner.bin_d_range(i_bin)
                    f_selection = self.Fobs.resolution_filter(d_min=dmin, d_max=dmax)
                    fsel_data = f_selection.data().as_numpy_array()
                    if self.log_fcells:
                        fsel_data = np_log(fsel_data)
                    sigma = SQRT(mean(f_selection.data().as_numpy_array() ** 2))
                    sigmas.append(sigma)  # sigma_for_res_id[i_bin] = sigma
                # min_sigma = min(self.sigma_for_res_id.values())
                # max_sigma = max(self.sigma_for_res_id.values())
                # median_sigma = np.median(self.sigma_for_res_id.values())
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

                if self.rescale_params:
                    self.resolution_ids_from_i_fcell = ARRAY([self.res_group_id_from_fcell_index[i_fcell] for i_fcell in range(self.n_global_fcell)])
                    self.fcell_sigmas_from_i_fcell = ARRAY([ self.sigma_for_res_id[res_id]*self.fcell_sigma_scale for res_id in self.resolution_ids_from_i_fcell])
                    self.fcell_init_from_i_fcell = ARRAY(self.fcell_init)

    def determine_parameter_freeze_order(self):
        param_sels = []
        if self.refine_detdist:
            param_sels.append(self.origin_sel)
        if self.refine_Umatrix:
            param_sels.append(self.umatrix_sel)
        if self.refine_Bmatrix:
            param_sels.append(self.bmatrix_sel)
        if self.refine_Fcell:
            param_sels.append(self.Fcell_sel)
        if self.refine_ncells:
            param_sels.append(self.ncells_sel)
        if self.refine_crystal_scale:
            param_sels.append(self.spot_scale_sel)
        if self.refine_background_planes:
            param_sels.append(self.bg_sel)

        self.param_sels = itertools.cycle(param_sels)

    def _update_Xall_from_x(self):
        """update the master parameter array with values from the parameter array that LBFGS sees"""
        for i, val in enumerate(self.x):
            if self.only_pass_refined_x_to_lbfgs:
                Xall_pos = self.x2xall[i]
            else:
                Xall_pos = i
            self.Xall[Xall_pos] = val
        # TODO maybe its more quick to use set_selected with is_being_refined ?
        # seomthing like self.Xall.set_selected(self.is_being_refined, self.x)

    def _make_x_identifier_array(self):
        """do this in case we need to identify what the parameters in X are at a later time"""
        #TODO update for new parameter types
        parameter_dict = {}
        for i_shot in range(self.n_shots):
            if self.FNAMES is None:
                fname = "rank%d_shot%d" % (self.rank, i_shot)
            else:
                fname = self.FNAMES[i_shot]

            proc_fname = "null"
            if self.PROC_FNAMES is not None:
                proc_fname = self.PROC_FNAMES[i_shot]

            proc_idx = -1
            if self.PROC_IDX is not None:
                proc_idx = int(self.PROC_IDX[i_shot])

            img_fname = fname
            parameter_dict[img_fname] = {}
            parameter_dict[img_fname]["agg_file"] = proc_fname
            parameter_dict[img_fname]["agg_idx"] = proc_idx
            parameter_dict[img_fname]["x_pos"] = {}
            PD = parameter_dict[img_fname]["x_pos"]

            # save the background tilt plane coefficients identifiers
            if not self.bg_extracted:
                nspots_on_shot = len(self.NANOBRAGG_ROIS[i_shot])
                for i_roi in range(nspots_on_shot):
                    i_a = self.bg_a_xstart[i_shot][i_roi]
                    i_b = self.bg_b_xstart[i_shot][i_roi]
                    i_c = self.bg_c_xstart[i_shot][i_roi]
                    if self.BBOX_IDX is not None:
                        bbox_idx = self.BBOX_IDX[i_shot][i_roi]
                    else:
                        bbox_idx = -1
                    PD[i_a] = "t1", bbox_idx
                    PD[i_b] = "t2", bbox_idx
                    PD[i_c] = "t3", bbox_idx
            else:
                i_bg_coef = self.bg_coef_xpos[i_shot]
                PD[i_bg_coef] = "bg_coef"

            # save the rotation angles identifier
            i_rotX = self.rotX_xpos[i_shot]
            i_rotY = self.rotY_xpos[i_shot]
            i_rotZ = self.rotZ_xpos[i_shot]
            PD[i_rotX] = "rX"
            PD[i_rotY] = "rY"
            PD[i_rotZ] = "rZ"

            # save unit cell variables identifier
            ucell_man = self.UCELL_MAN[i_shot]
            names = ucell_man.variable_names
            for i_name, name in enumerate(names):
                i_uc = self.ucell_xstart[i_shot] + i_name
                PD[i_uc] = name

            # save the ncells identifier
            if self.n_ncells_param == 1:
                ncells_names = ("m",)
            elif self.n_ncells_param == 2:
                ncells_names = "N1", "N2"
            else:  # only other choice is n_ncells_param=3
                ncells_names = "Na", "Nb", "Nc"
            for i_nc in range(self.n_ncells_param):
                i_ncells = self.ncells_xstart[i_shot] + i_nc
                PD[i_ncells] = ncells_names[i_nc]

            # save spot scale indentifier
            i_scale = self.spot_scale_xpos[i_shot]
            PD[i_scale] = "Gs"

        self._MPI_write_output(parameter_dict)

    def _MPI_write_output(self, parameter_dict):
        if self.output_dir is not None:
            outdir = PATHJOIN(self.output_dir, "parameter_id")
            if self.I_AM_ROOT and not EXISTS(outdir):
                MAKEDIRS(outdir)
            all_data = self._data_for_write(parameter_dict)
            if self.I_AM_ROOT:
                for i_pd, PD in enumerate(all_data):
                    output_path = PATHJOIN(outdir, "param%d.json" % i_pd)
                    with open(output_path, "w") as out:
                        JSON_DUMP(PD, out)

    def _make_parameter_type_selection_arrays(self):  # experimental , not really used
        from cctbx.array_family import flex
        self.umatrix_sel = flex.bool(len(self.Xall), True)
        self.bmatrix_sel = flex.bool(len(self.Xall), True)
        self.Fcell_sel = flex.bool(len(self.Xall), True)
        self.origin_sel = flex.bool(len(self.Xall), True)
        self.spot_scale_sel = flex.bool(len(self.Xall), True)
        self.ncells_sel = flex.bool(len(self.Xall), True)
        self.bg_sel = flex.bool(len(self.Xall), True)
        for i_shot in range(self.n_shots):
            self.umatrix_sel[self.rotX_xpos[i_shot]] = False
            self.umatrix_sel[self.rotY_xpos[i_shot]] = False
            self.umatrix_sel[self.rotZ_xpos[i_shot]] = False

            for i_uc in range(self.n_ucell_param):
                self.bmatrix_sel[self.ucell_xstart[i_shot] + i_uc] = False

            for i_ncells in range(self.n_ncells_param):
                self.ncells_sel[self.ncells_xstart[i_shot] + i_ncells] = False
            self.spot_scale_sel[self.spot_scale_xpos[i_shot]] = False

            self.origin_sel[self.detector_distance_xpos[i_shot]] = False

            if not self.bg_extracted:
                nspots_on_shot = len(self.NANOBRAGG_ROIS[i_shot])
                for i_spot in range(nspots_on_shot):
                    self.bg_sel[self.bg_a_xstart[i_shot][i_spot]] = False
                    self.bg_sel[self.bg_b_xstart[i_shot][i_spot]] = False
                    self.bg_sel[self.bg_c_xstart[i_shot][i_spot]] = False

        for i_fcell in range(self.n_global_fcell):
            self.Fcell_sel[self.fcell_xstart + i_fcell] = False

    def _get_sausage_parameters(self, i_shot):
        vals = []
        for i_sausage in range(self.num_sausages):
            for i_param in range(4):
                idx = i_sausage*4 + i_param
                xpos = self.sausages_xpos[i_shot][idx]
                sigma = self.sausages_sigma[i_param]
                init = self.sausages_init[i_shot][idx] #_param]
                xval_init = 1 + i_sausage #idx #*0.1
                val = sigma*(self.Xall[xpos] - xval_init) + init
                vals.append(val)
        return vals

    def _get_rotX(self, i_shot):
        if self.rescale_params:
            # FIXME ?
            return self.rotX_sigma*(self.Xall[self.rotX_xpos[i_shot]]-1) + self.rotXYZ_inits[i_shot][0]
        else:
            return self.Xall[self.rotX_xpos[i_shot]]

    def _get_rotY(self, i_shot):
        if self.rescale_params:
            return self.rotY_sigma * (self.Xall[self.rotY_xpos[i_shot]] - 1) + self.rotXYZ_inits[i_shot][1]
        else:
            return self.Xall[self.rotY_xpos[i_shot]]

    def _get_rotZ(self, i_shot):
        if self.rescale_params:
            return self.rotZ_sigma * (self.Xall[self.rotZ_xpos[i_shot]] - 1) + self.rotXYZ_inits[i_shot][2]

        else:
            return self.Xall[self.rotZ_xpos[i_shot]]

    def _get_spectra_coefficients(self):
        vals = []
        if self.rescale_params:
            for i in range(self.n_spectra_param):
                xval = self.Xall[self.spectra_coef_xstart + i]
                sig = self.spectra_coefficients_sigma[i]
                init = self.spectra_coefficients_init[i]
                low, high = self.lambda_coef_ranges[i]
                rng = high-low
                sin_arg = sig*(xval-1) + ASIN(2*(init-low)/rng - 1)
                val = (SIN(sin_arg) + 1)*rng/2 + low
                #val = sig*(xval-1) + init
                vals.append(val)
        else:
            assert NotImplementedError
        return vals

    def _get_ucell_vars(self, i_shot):
        all_p = []
        for i in range(self.n_ucell_param):
            if self.rescale_params:
                x = self.Xall[self.ucell_xstart[i_shot]+i]
                if self.use_ucell_ranges:
                    p = self.ucell_params[i_shot][i].get_val(x)
                else:
                    sig = self.ucell_sigmas[i]
                    init = self.ucell_inits[i_shot][i]
                    p = sig*(self.Xall[self.ucell_xstart[i_shot]+i] - 1) + init
            else:
                p = self.Xall[self.ucell_xstart[i_shot]+i]
            all_p.append(p)
        return all_p

    def _get_panelRot_val(self, panel_id):
        vals = []
        for i_rot in range(3):
            panel_group_id = self.panel_group_from_id[panel_id]
            lbfgs_xval = self.Xall[self.panelRot_xstart + 3*panel_group_id + i_rot]
            val = self.panelRot_params[panel_group_id][i_rot].get_val(lbfgs_xval)
            vals.append(val)
        return vals

    def _get_panelXYZ_val(self, panel_id, i_shot=0):
        panel_group_id = self.panel_group_from_id[panel_id]
        xpos_X = self.panelXY_xstart + 2*panel_group_id
        valX = self.Xall[xpos_X]
        offsetX = self.panelX_params[panel_group_id].get_val(valX)

        xpos_Y = xpos_X + 1
        valY = self.Xall[xpos_Y]
        offsetY = self.panelY_params[panel_group_id].get_val(valY)

        xpos_Z = self.panelZ_xstart + self.panel_group_from_id[panel_id]
        valZ = self.Xall[xpos_Z]
        offsetZ = self.panelZ_params[i_shot].get_val(valZ)

        return offsetX, offsetY, offsetZ

    def _get_detector_distance_val(self, i_shot):
        val = 0
        # TODO will this work with multi-step refinement ?
        if self.refine_detdist:
            xval = self.Xall[self.detector_distance_xpos[i_shot]]
            val = self.detector_distance_params[i_shot].get_val(xval)
        elif self.pershot_detdist_shifts:
            val = self.shot_detector_distance_init[i_shot]
        return val

    def _get_ncells_def_vals(self, i_shot):
        vals = []
        for i_ncell in range(3):
            xval = self.Xall[self.ncells_def_xstart[i_shot] + i_ncell]
            val = self.Ncells_def_params[i_shot][i_ncell].get_val(xval)
            vals.append(val)
        return vals

    def _get_m_val(self, i_shot):
        # NOTE This is stupidly complicated because it supports using 1,2 or 3 free NCELLS parameters, for example
        # NOTE   allowing a refinement where Na=Nb != Nc . This functionality doesnt really make sense,
        # NOTE   and should be dropped for readability
        vals = []
        if self.S.D.isotropic_ncells:
            xval = self.Xall[self.ncells_xstart[i_shot]]
            if self.use_Ncells_range:
                val = self.Ncells_params[i_shot][0].get_val(xval)

            elif self.rescale_params:
                sig = self.m_sigma
                try:
                    sig = sig[0]
                except TypeError:
                    pass
                init = self.m_init[i_shot]
                try:
                    init = init[0]
                except TypeError:
                    pass
                val = np_exp(sig*(xval-1))*(init-3) + 3
            else:
                val = np_exp(xval)+3
            vals.append(val)
        else:
            for i_ncell in range(self.n_ncells_param):
                xval = self.Xall[self.ncells_xstart[i_shot] + i_ncell]
                if self.use_Ncells_range:
                    val = self.Ncells_params[i_shot][i_ncell].get_val(xval)
                elif self.rescale_params:
                    try:
                        sig = self.m_sigma[i_ncell]
                    except (TypeError, IndexError):
                        sig = self.m_sigma
                    init = self.m_init[i_shot][i_ncell]
                    val = np_exp(sig * (xval - 1)) * (init - 3) + 3
                else:
                    val = np_exp(xval) + 3
                vals.append(val)
            if self.ncells_mask is not None and self.n_ncells_param == 2:
                vals = [vals[mask_val] for mask_val in self.ncells_mask]
        return vals

    def _get_eta(self, i_shot):
        vals = [0, 0, 0]
        for i_eta in range(3):
            xval = self.Xall[self.eta_xstart[i_shot]+ i_eta]
            val = self.eta_params[i_shot][i_eta].get_val(xval)
            vals[i_eta] = val
            if not self.S.crystal.has_anisotropic_mosaicity:
                break
        return vals

    def _get_spot_scale(self, i_shot):
        xval = self.Xall[self.spot_scale_xpos[i_shot]]
        if self.rescale_params:
            sig = self.spot_scale_sigma
            init = self.spot_scale_init[i_shot]
            val = sig*(xval-1) + init
            #val = val * val
            #val = np_exp(sig*(xval-1))*init
        else:
            val = np_exp(xval)
        return val

    def _get_bg_coef(self, i_shot):
        assert (self.rescale_params)
        val = self.Xall[self.bg_coef_xpos[i_shot]]
        sig = self.bg_coef_sigma
        init = self.shot_bg_coef[i_shot]
        val = np_exp(sig*(val-1))*init
        return val

    def _set_spot_scale(self, new_val, i_shot):
        """just used in testsing derivatives"""
        if self.rescale_params:
            self.spot_scale_init[0] = new_val
            self.Xall[self.spot_scale_xpos[0]] = 1
        else:
            self.Xall[self.spot_scale_xpos[0]] = np_log(new_val)

    def _get_bg_vals(self, i_shot, i_spot):
        a_val = self.Xall[self.bg_a_xstart[i_shot][i_spot]]
        b_val = self.Xall[self.bg_b_xstart[i_shot][i_spot]]
        c_val = self.Xall[self.bg_c_xstart[i_shot][i_spot]]
        if self.rescale_params:
            a_sig = self.a_sigma
            b_sig = self.b_sigma
            c_sig = self.c_sigma
            a_init, b_init, c_init = self.ABC_INIT[i_shot][i_spot]
            a = a_sig*(a_val-1) + a_init
            b = b_sig*(b_val-1) + b_init
            c = c_sig*(c_val-1) + c_init
            if self.bg_offset_positive:
                c = np_exp(c_sig*(c_val-1))*c_init

        elif self.bg_offset_positive:
            a = a_val
            b = b_val
            c = np_exp(c_val)
        else:
            a = a_val
            b = b_val
            c = c_val

        if self.bg_offset_only:
            a = b = 0

        return a, b, c

    def _unpack_internal(self, lst, lst_is_x=False):
        # x = self..as_numpy_array()
        # note n_shots should be specific for this rank
        if lst_is_x:
            rotx = [self._get_rotX(i_shot) for i_shot in range(self.n_shots)]
            roty = [self._get_rotY(i_shot) for i_shot in range(self.n_shots)]
            rotz = [self._get_rotZ(i_shot) for i_shot in range(self.n_shots)]
        else:
            rotx = [lst[self.rotX_xpos[i_shot]] for i_shot in range(self.n_shots)]
            roty = [lst[self.rotY_xpos[i_shot]] for i_shot in range(self.n_shots)]
            rotz = [lst[self.rotZ_xpos[i_shot]] for i_shot in range(self.n_shots)]

        if self.global_ncells_param:
            if lst_is_x:
                ncells_vals = [self._get_m_val(0)[0]] * len(rotx)
            else:
                ncells_vals = [lst[self.ncells_xstart[0]]] * len(rotx)
        else:
            if lst_is_x:
                ncells_vals = [self._get_m_val(i_shot)[0] for i_shot in range(self.n_shots)]
            else:
                ncells_vals = [lst[self.ncells_xstart[i_shot]] for i_shot in range(self.n_shots)]

        if lst_is_x:
            detector_distance_vals = [self._get_detector_distance_val(i_shot) for i_shot in range(self.n_shots)]
        else:
            detector_distance_vals = [lst[self.detector_distance_xpos[i_shot]] for i_shot in range(self.n_shots)]

        if lst_is_x:
            #ncells_vals = list(np_exp(ncells_vals)+3)
            scale_vals = [self._get_spot_scale(i_shot)**2 for i_shot in range(self.n_shots)]
        else:
            scale_vals = [lst[self.spot_scale_xpos[i_shot]] for i_shot in range(self.n_shots)]

        # this can be used to compare directly
        if self.CRYSTAL_SCALE_TRUTH is not None:
            scale_vals_truths = [self.CRYSTAL_SCALE_TRUTH[i_shot] for i_shot in range(self.n_shots)]
        else:
            scale_vals_truths = None

        if self.global_ucell_param:
            if lst_is_x:
                ucparams = self._get_ucell_vars(0)
            else:
                ucparams = lst[self.ucell_xstart[0]:self.ucell_xstart[0] + self.n_ucell_param]
            ucparams_lsts = []
            for ucp in ucparams:
                ucparams_lsts.append([ucp]*len(rotx))
        else:
            if lst_is_x:
                all_shot_params = [self._get_ucell_vars(i_shot) for i_shot in range(self.n_shots)]
                ucparams_lsts = list(map(list, zip(*all_shot_params)))
            else:
                ucparams_lsts = []
                for i_uc in range(self.n_ucell_param):
                    ucp_lst = [lst[self.ucell_xstart[i_shot] + i_uc] for i_shot in range(self.n_shots)]
                    ucparams_lsts.append(ucp_lst)

        rotx = self._MPI_reduce_broadcast(rotx)
        roty = self._MPI_reduce_broadcast(roty)
        rotz = self._MPI_reduce_broadcast(rotz)
        ncells_vals = self._MPI_reduce_broadcast(ncells_vals)
        scale_vals = self._MPI_reduce_broadcast(scale_vals)
        detector_distance_vals = self._MPI_reduce_broadcast(detector_distance_vals)

        ucparams_all = []
        for ucp in ucparams_lsts:
            ucp = self._MPI_reduce_broadcast(ucp)
            ucparams_all.append(ucp)
        if scale_vals_truths is not None:
            scale_vals_truths = self._MPI_reduce_broadcast(scale_vals_truths)

        return rotx, roty, rotz, ucparams_all, ncells_vals, scale_vals, scale_vals_truths, detector_distance_vals

    def _send_ucell_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i + 3,
                self.UCELL_MAN[self._i_shot].derivative_matrices[i])
            if self.calc_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i + 3, self.UCELL_MAN[self._i_shot].second_derivative_matrices[i])

    def _run_diffBragg_current(self):
        self.pfs, self.roi_ids, _ = self._get_panels_fasts_slows()
        if self.randomize_devices is not None:
            dev = np_random_choice(self.randomize_devices)
            self.D.device_Id = dev
        #p = 136, f = 110, s= 172
        #pfs = [(136, 110, 172)]
        #if self.special_flag is not None:
        #p,f,s = self.pfs[:2]
        #if True: #self.special_flag is not None: #self.FNAMES[self._i_shot]=="/global/cscratch1/sd/dermen/all_3/expers/rank2/stg1_run796_shot655_indexed_0.expt"
        #    p,f,s = (2, 149, 236)
        #    self.D.printout_pixel_fastslow=f,s
        #    self.D.add_diffBragg_spots((p,f,s))
        #    exit()
        #pids = self.pfs[0::3]
        #xs = self.pfs[1::3]
        #ys = self.pfs[2::3]
        #SAVEZ("shit", pids=pids, ys=ys, xs=xs)
        #exit()
        self.D.add_diffBragg_spots(self.pfs)

    def _get_fcell_val(self, i_fcell):
        # TODO vectorize me
        # i_fcell is between 0 and self.n_global_fcell
        # get the asu index and its updated amplitude
        xpos = self.fcell_xstart + i_fcell
        val = self.Xall[xpos]  # new amplitude
        if self.rescale_params:
            resolution_id = self.res_group_id_from_fcell_index[i_fcell]  # TODO
            sig = self.sigma_for_res_id[resolution_id]*self.fcell_sigma_scale  # TODO
            init = self.fcell_init[i_fcell]
            if self.log_fcells:
                val = np_exp(sig*(val - 1))*init
            else:
                if val < 0:  # NOTE this easily happens without the log c.o.v.
                    self.Xall[xpos] = 0
                    val = 0
                    self.num_Fcell_negative_model += 1
                else:
                    val = sig*(val - 1) + init

        else:
            if self.log_fcells:
                val = np_exp(val)
            if val < 0:  # NOTE this easily happens without the log c.o.v.
                self.Xall[xpos] = 0
                val = 0
                self.num_Fcell_negative_model += 1
        return val

    def _store_updated_Fcell(self):
        if not self.refine_Fcell:
            return
        xvals = self.Xall[self.fcell_xstart: self.fcell_xstart+self.n_global_fcell]
        if self.rescale_params and self.log_fcells:
            sigs = self.fcell_sigmas_from_i_fcell
            inits = self.fcell_init_from_i_fcell
            if self.log_fcells:
                vals = np_exp(sigs*(xvals - 1))*inits
            else:
                vals = sigs*(xvals - 1) + inits
                vals[vals < 0] = 0
        else:
            if self.log_fcells:
                vals = np_exp(xvals)
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

        self.S.D.Fhkl_tuple = idx, flex_double(data)  # update nanoBragg again  # TODO: add flag to not re-allocate in nanoBragg!

    def _update_spectra_coefficients(self):
        if self.refine_lambda0 or self.refine_lambda1 or self.update_spectra_during_refinement:
            coeffs = self._get_spectra_coefficients()
            self.D.lambda_coefficients = tuple(coeffs)

    def _update_eta(self):
        if self.update_eta or self.refine_eta:  # TODO think of a better interface
            eta_val = self._get_eta(self._i_shot)
            Umats_dblprime = None
            if self.S.Umats_method == 1:
                assert eta_val[1] == 0 and eta_val[2] == 0
                from simtbx.nanoBragg.tst_gaussian_mosaicity2 import run_uniform
                Umats, Umats_prime = run_uniform(eta_val, self.D.mosaic_domains)
            elif self.S.Umats_method in [2, 3]:
                if self.S.Umats_method == 2:
                    assert eta_val[1] == eta_val[2] == 0
                    eta_val = eta_val[0]
                Umats,Umats_prime,Umats_dblprime = self.S.Umats(eta_val,
                                                            n_mos_doms=self.S.crystal.n_mos_domains,
                                                            crystal=self.crystal_for_mosaicity_model,
                                                            method=self.S.Umats_method,
                                                            angles_per_axis=self.S.crystal.mos_angles_per_axis,
                                                            num_axes=self.S.crystal.num_mos_axes)
            else:
                raise ValueError("Wrong umats method, only support for 1,2,3")
            if isinstance(eta_val, Iterable):
                eta_val_s = "%.4f, %.4f, %.4f" % tuple(eta_val)
            else:
                eta_val_s = "%.4f" % eta_val
            self.print("updating ETA as %s deg sampled from %d blocks from method %d" %
                       (eta_val_s, self.D.mosaic_domains, self.S.Umats_method))

            self.D.set_mosaic_blocks(Umats)

            if self.S.Umats_method==3:
                Umats_prime = Umats_prime[0::3] + Umats_prime[1::3] + Umats_prime[2::3]
            self.D.set_mosaic_blocks_prime(Umats_prime)

            if Umats_dblprime is not None: # and self.compute_curvatures:
                if self.S.Umats_method == 3:
                    Umats_dblprime = Umats_dblprime[0::3] + Umats_dblprime[1::3] + Umats_dblprime[2::3]
                self.D.set_mosaic_blocks_dbl_prime(Umats_dblprime)

            self.D.vectorize_umats()

    def _set_background_plane(self, i_spot):
        saving_model = self.save_model_for_shot is not None and self.save_model_for_shot == self._i_shot
        try:
            if self.bg_extracted:
                self.bg_coef = self._get_bg_coef(self._i_shot)
                (i1, i2), (j1, j2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
                self.tilt_plane = self.bg_coef*self.background_estimate[self._panel_id, j1:j2, i1:i2]
            elif self.background is not None:
                (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][self._i_spot]
                self.tilt_plane = self.background[self._i_shot][self._panel_id, y1:y2, x1:x2]
                assert np_all(self.tilt_plane >= 0)
            else:
                xr = self.XREL[self._i_shot][i_spot]
                yr = self.YREL[self._i_shot][i_spot]
                self.a, self.b, self.c = self._get_bg_vals(self._i_shot, i_spot)
                if self.bg_offset_only:
                    self.tilt_plane = ONES_LIKE(xr)*self.c
                else:
                    self.tilt_plane = xr * self.a + yr * self.b + self.c
                if self.OMEGA_KAHN is not None:
                    (i1, i2), (j1, j2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
                    omega_kahn_correction = self.OMEGA_KAHN[self._panel_id][j1:j2, i1:i2]
                    self.tilt_plane *= omega_kahn_correction
        except Exception as error:
            if not saving_model:
                raise("background error", error)
            else:
                im = self.ROI_IMGS[self._i_shot][i_spot]
                self.tilt_plane = ONES_LIKE(im)
                m = median(im[ im < percentile(im,80)])
                self.tilt_plane *= m

        self.tilt_plane = self.tilt_plane.flatten()

    def _update_sausages(self):
        if self.num_sausages > 0:
            sausage_vals = self._get_sausage_parameters(self._i_shot)
            rotX = flex_double(sausage_vals[0::4])
            rotY = flex_double(sausage_vals[1::4])
            rotZ = flex_double(sausage_vals[2::4])
            scales = flex_double(sausage_vals[3::4])
            self.D.set_sausages(rotX, rotY, rotZ, scales)
            print("Set the sausages!")

    def _update_rotXYZ(self):
        #if self.refine_rotX:
        self.D.set_value(0, self._get_rotX(self._i_shot))
        #if self.refine_rotY:
        self.D.set_value(1, self._get_rotY(self._i_shot))
        #if self.refine_rotZ:
        self.D.set_value(2, self._get_rotZ(self._i_shot))

    def _update_ncells(self):
        vals = self._get_m_val(self._i_shot)
        if self.D.isotropic_ncells:
            self.D.set_value(self._ncells_id, vals[0])
        else:
            self.D.set_ncells_values(tuple(vals))

    def _update_ncells_def(self):
        vals = self._get_ncells_def_vals(self._i_shot)
        self.D.Ncells_def = tuple(vals)

    def _update_dxtbx_detector(self):
        refining_det = any((self.refine_panelRotO, self.refine_panelRotF,
                           self.refine_panelRotS, self.refine_panelXY, self.refine_panelZ))

        if refining_det or self.update_detector_during_refinement:
            npanels = len(self.S.detector)
            for pid in range(npanels):
                new_offsetX, new_offsetY, new_offsetZ = self._get_panelXYZ_val(pid)

                if self.refine_detdist:  # TODO: figure out if this should override the above, or add to it?
                    new_offsetZ = self._get_detector_distance_val(self._i_shot)

                panel_rot_angO, panel_rot_angF, panel_rot_angS = self._get_panelRot_val(pid)

                if self.panel_reference_from_id is not None:
                    self.D.reference_origin = self.panel_reference_from_id[pid]  #self.S.detector[self.panel_reference_from_id[pid]].get_origin()
                else:
                    self.D.reference_origin = self.S.detector[pid].get_origin()

                self.D.update_dxtbx_geoms(self.S.detector, self.S.beam.nanoBragg_constructor_beam, pid,
                                          panel_rot_angO, panel_rot_angF, panel_rot_angS, new_offsetX, new_offsetY, new_offsetZ,
                                          force=False)
        elif self.pershot_detdist_shifts:  # TODO we are hijacking the API to use the shifts from hopper in the legacy code, clean this up
            shiftZ = self._get_detector_distance_val(self._i_shot)
            #self.print("DETZ shift: %1.3e meters" % shiftZ)
            self.S.D.shift_origin_z(self.S.detector,  shiftZ)

    def _extract_spectra_coefficient_derivatives(self):
        self.spectra_derivs = [0]*self.n_spectra_param
        if self.refine_lambda0 or self.refine_lambda1:
            SG = self.scale_fac * self.G2
            if self.refine_lambda0:
                self.spectra_derivs[0] = SG*self._lam0_deriv[self.roi_slice]
            if self.refine_lambda1:
                self.spectra_derivs[1] = SG*self._lam1_deriv[self.roi_slice]

    def _pre_extract_deriv_arrays(self):
        npix = int(len(self.pfs) / 3)
        self._model_pix = self.D.raw_pixels_roi.as_numpy_array()[:npix]

        if self.refine_blueSausages:
            self._sausage_derivs = self.D.get_sausage_derivative_pixels()
            self._sausage_derivs = [d.as_numpy_array()[:npix] for d in self._sausage_derivs]

        ## UMAT
        self._UmatX_deriv_pixels = self._UmatY_deriv_pixels = self._UmatZ_deriv_pixels = 0
        self._UmatX_second_deriv_pixels = self._UmatY_second_deriv_pixels = self._UmatZ_second_deriv_pixels = 0

        if self.refine_Umatrix:

            self._UmatX_deriv_pixels = self.D.get_derivative_pixels(0).as_numpy_array()[:npix]
            self._UmatY_deriv_pixels = self.D.get_derivative_pixels(1).as_numpy_array()[:npix]
            self._UmatZ_deriv_pixels = self.D.get_derivative_pixels(2).as_numpy_array()[:npix]

            if self.calc_curvatures:
                self._UmatX_second_deriv_pixels = self.D.get_second_derivative_pixels(0).as_numpy_array()[:npix]
                self._UmatY_second_deriv_pixels = self.D.get_second_derivative_pixels(1).as_numpy_array()[:npix]
                self._UmatZ_second_deriv_pixels = self.D.get_second_derivative_pixels(2).as_numpy_array()[:npix]

        # BMAT
        self._Bmat_derivs = [0] * self.n_ucell_param
        self._Bmat_second_derivs = [0] * self.n_ucell_param
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self._Bmat_derivs[i] = self.D.get_derivative_pixels(i+3).as_numpy_array()[:npix]
                if self.calc_curvatures:
                    self._Bmat_second_derivs[i] = self.D.get_second_derivative_pixels(i+3).as_numpy_array()[:npix]

        self._lam0_deriv = self._lam1_deriv = 0
        if self.refine_lambda0:
            self._lam0_deriv = self.D.get_derivative_pixels(12).as_numpy_array()[:npix]
        if self.refine_lambda1:
            self._lam1_deriv = self.D.get_derivative_pixels(13).as_numpy_array()[:npix]

        refining_pan_rot = [self.refine_panelRotO, self.refine_panelRotF, self.refine_panelRotS]
        refining_pan_XYZ = [self.refine_panelXY, self.refine_panelXY, self.refine_panelZ or self.refine_detdist]
        rot_manager_ids = self._panelRotO_id, self._panelRotF_id, self._panelRotS_id
        xyz_manager_ids = self._panelX_id, self._panelY_id, self._detector_distance_id

        self._panRot_deriv = [0, 0, 0]
        self._panRot_second_deriv = [0, 0, 0]
        self._panXYZ_deriv = [0, 0, 0]
        self._panXYZ_second_deriv = [0, 0, 0]
        self._eta_deriv = [0,0,0]
        self._eta_second_deriv = [0,0,0]
        for i_pan in range(3):
            if refining_pan_rot[i_pan]:
                rot_man_id = rot_manager_ids[i_pan]
                self._panRot_deriv[i_pan] = self.D.get_derivative_pixels(rot_man_id).as_numpy_array()[:npix]
                if self.calc_curvatures:
                    self._panelRot_second_deriv[i_pan] = self.D.get_second_derivative_pixels(rot_man_id).as_numpy_array()[:npix]

            if refining_pan_XYZ[i_pan]:
                xyz_man_id = xyz_manager_ids[i_pan]
                self._panXYZ_deriv[i_pan] = self.D.get_derivative_pixels(xyz_man_id).as_numpy_array()[:npix]
                if self.calc_curvatures:
                    self._panXYZ_second_deriv[i_pan] = self.D.get_second_derivative_pixels(xyz_man_id).as_numpy_array()[:npix]

        if self.refine_Fcell:
            dF = self.D.get_derivative_pixels(self._fcell_id)

            self._extracted_fcell_deriv = dF.set_selected(dF != dF, 0)
            self._extracted_fcell_deriv = self._extracted_fcell_deriv.as_numpy_array()[:npix]
            if self.calc_curvatures:
                d2F = self.D.get_second_derivative_pixels(self._fcell_id)
                self._extracted_fcell_second_deriv = d2F.set_selected(d2F != d2F, 0)
                self._extracted_fcell_second_deriv = self._extracted_fcell_second_deriv.as_numpy_array()[:npix]

        if self.refine_ncells:
            if self.D.isotropic_ncells:
                self._iso_ncells_deriv = self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()[:npix]
                if self.calc_curvatures:
                    self._iso_ncells_second_deriv = self.D.get_second_derivative_pixels(self._ncells_id).as_numpy_array()[:npix]
            else:
                self._aniso_ncells_derivs = self.D.get_ncells_derivative_pixels()
                self._aniso_ncells_derivs = [d.as_numpy_array()[:npix]for d in self._aniso_ncells_derivs]
                if self.calc_curvatures:
                    self._aniso_ncells_second_derivs = self.D.get_ncells_second_derivative_pixels()
                    self._aniso_ncells_second_derivs = [d.as_numpy_array()[:npix]for d in self._aniso_ncells_second_derivs]

        if self.refine_ncells_def:
            self._ncells_def_derivs = [d.as_numpy_array()[:npix] for d in self.D.get_ncells_def_derivative_pixels()]
            if self.calc_curvatures:
                self._ncells_def_second_derivs = [d.as_numpy_array()[:npix] for d in self.D.get_ncells_def_second_derivative_pixels()]

        if self.refine_eta:
            if not self.S.crystal.has_anisotropic_mosaicity:
                self._eta_deriv[0] = self.D.get_derivative_pixels(self._eta_id).as_numpy_array()[:npix]
            else:
                self._eta_deriv = [d.as_numpy_array()[:npix] for d in self.D.get_aniso_eta_deriv_pixels()]
            if self.calc_curvatures:
                if not self.S.crystal.has_anisotropic_mosaicity:
                    self._eta_second_deriv[0] = self.D.get_second_derivative_pixels(self._eta_id).as_numpy_array()[:npix]
                else:
                    self._eta_second_deriv = [d.as_numpy_array()[:npix] for d in self.D.get_aniso_eta_second_deriv_pixels()]

    def _extract_sausage_derivs(self):
        if self.refine_blueSausages:
            self.sausages_dI_dtheta = [0]*(self.num_sausages*4)
            SG = self.scale_fac*self.G2
            for i_sausage in range(self.num_sausages):
                for i_param in range(4):
                    idx = i_sausage*4 + i_param
                    self.sausages_dI_dtheta[idx] = SG*self._sausage_derivs[idx][self.roi_slice]

    def _extract_Umatrix_derivative_pixels(self):
        self.rotX_dI_dtheta = self.rotY_dI_dtheta = self.rotZ_dI_dtheta = 0
        self.rotX_d2I_dtheta2 = self.rotY_d2I_dtheta2 = self.rotZ_d2I_dtheta2 = 0
        # convenient storage of the gain and scale as a single parameter
        SG = self.scale_fac*self.G2
        if self.refine_Umatrix:
            if self.refine_rotX:
                self.rotX_dI_dtheta = SG*self._UmatX_deriv_pixels[self.roi_slice]
                if self.calc_curvatures:
                    self.rotX_d2I_dtheta2 = SG*self._UmatX_second_deriv_pixels[self.roi_slice]

            if self.refine_rotY:
                self.rotY_dI_dtheta = SG*self._UmatY_deriv_pixels[self.roi_slice]
                if self.calc_curvatures:
                    self.rotY_d2I_dtheta2 = SG*self._UmatY_second_deriv_pixels[self.roi_slice]

            if self.refine_rotZ:
                self.rotZ_dI_dtheta = SG*self._UmatZ_deriv_pixels[self.roi_slice]
                if self.calc_curvatures:
                    self.rotZ_d2I_dtheta2 = SG*self._UmatZ_second_deriv_pixels[self.roi_slice]

    def _extract_Bmatrix_derivative_pixels(self):
        # the Bmatrix derivatives are stored for each unit cell parameter (UcellManager.variables)
        self.ucell_dI_dtheta = [0] * self.n_ucell_param
        self.ucell_d2I_dtheta2 = [0] * self.n_ucell_param
        SG = self.scale_fac*self.G2
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.ucell_dI_dtheta[i] = SG*self._Bmat_derivs[i][self.roi_slice]
                if self.calc_curvatures:
                    self.ucell_d2I_dtheta2[i] = SG*self._Bmat_second_derivs[i][self.roi_slice]

    def _extract_ncells_def_derivative_pixels(self):
        self.ncells_def_dI_dtheta = [0]*3
        self.ncells_def_d2I_dtheta2 = [0]*3
        if self.refine_ncells_def:
            SG = self.scale_fac * self.G2
            for i_ncell in range(3):
                d = self._ncells_def_derivs[i_ncell][self.roi_slice]
                self.ncells_def_dI_dtheta[i_ncell] = SG*d
                if self.calc_curvatures:
                    d2 = self._ncells_def_second_derivs[i_ncell][self.roi_slice]
                    self.ncells_def_d2I_dtheta2[i_ncell] = SG*d2

    def _extract_mosaic_parameter_m_derivative_pixels(self):
        SG = self.scale_fac * self.G2
        if self.D.isotropic_ncells:  # TODO remove need for if/else
            self.m_dI_dtheta = self.m_d2I_dtheta2 = 0
            if self.refine_ncells:
                self.m_dI_dtheta = SG*self._iso_ncells_deriv[self.roi_slice]
                if self.calc_curvatures:
                    self.m_d2I_dtheta2 = SG*self._iso_ncells_second_deriv[self.roi_slice]
            self.m_dI_dtheta = [self.m_dI_dtheta]
            self.m_d2I_dtheta2 = [self.m_d2I_dtheta2]
        else:
            self.m_dI_dtheta = [0] * self.n_ncells_param
            self.m_d2I_dtheta2 = [0] * self.n_ncells_param
            if self.refine_ncells:
                derivs = self._aniso_ncells_derivs
                if self.calc_curvatures:
                    second_derivs = self._aniso_ncells_second_derivs

                # ncells_mask allows specification of two Ncells parameters to be fixed to each other
                # but the diffBragg instance doesnt know about this, rather we control this entirely in the LocalRefiner
                # class here. Hence derivs and second derivs will always be a list of length 3.
                # The attribute `ncells_mask` specifies which Ncells params are fixed, hence  we use it to reduce derivs
                # to length 2, if necessary
                if self.ncells_mask is not None:
                    temp_derivs = [0, 0]
                    temp_second_derivs = [0, 0]
                    for i_mask, mask_val in enumerate(self.ncells_mask):
                        temp_derivs[mask_val] = derivs[i_mask]
                        if self.calc_curvatures:
                            temp_second_derivs[mask_val] = second_derivs[i_mask]
                    derivs = temp_derivs
                    second_derivs = temp_second_derivs

                for i_ncell in range(self.n_ncells_param):
                    d = derivs[i_ncell][self.roi_slice]
                    self.m_dI_dtheta[i_ncell] = SG*d
                    if self.calc_curvatures:
                        d2 = second_derivs[i_ncell][self.roi_slice]
                        self.m_d2I_dtheta2[i_ncell] = SG*d2

    def _extract_detector_distance_derivative_pixels(self):
        self.detdist_dI_dtheta = self.detdist_d2I_dtheta2 = 0
        SG = self.scale_fac*self.G2
        if self.refine_detdist:
            self.detdist_dI_dtheta = SG*self._panXYZ_deriv[2][self.roi_slice]
            if self.calc_curvatures:
                # TODO speed up with pre-extracttion
                self.detdist_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(self._detector_distance_id)[self.roi_slice].as_numpy_array()

    def _extract_panelRot_derivative_pixels(self):
        self.panelRot_dI_dtheta = [0, 0, 0]
        self.panelRot_d2I_dtheta2 = [0, 0, 0]
        SG = self.scale_fac*self.G2
        refining = [self.refine_panelRotO, self.refine_panelRotF, self.refine_panelRotS]
        for i_rot in range(3):
            if refining[i_rot]:
                self.panelRot_dI_dtheta[i_rot] = SG*self._panRot_deriv[i_rot][self.roi_slice]
                if self.calc_curvatures:
                    self.panelRot_d2I_dtheta2[i_rot] = SG*self._panRot_second_deriv[i_rot][self.roi_slice]

    def _extract_panelXYZ_derivative_pixels(self):
        self.panelX_dI_dtheta = self.panelX_d2I_dtheta2 = 0
        self.panelY_dI_dtheta = self.panelY_d2I_dtheta2 = 0
        self.panelZ_dI_dtheta = self.panelZ_d2I_dtheta2 = 0
        SG = self.scale_fac*self.G2
        # TODO: curvatures
        if self.refine_panelXY:
            self.panelX_dI_dtheta = SG*self._panXYZ_deriv[0][self.roi_slice]
            self.panelY_dI_dtheta = SG*self._panXYZ_deriv[1][self.roi_slice]
        if self.refine_panelZ:
            assert not self.refine_detdist
            self.panelZ_dI_dtheta = SG*self._panXYZ_deriv[2][self.roi_slice]

    def _extract_Fcell_derivative_pixels(self):
        # TODO pre-extract
        self.fcell_deriv = self.fcell_second_deriv = 0
        if self.refine_Fcell:
            SG = self.scale_fac*self.G2
            self.fcell_deriv = SG*(self._extracted_fcell_deriv[self.roi_slice])
            # handles Nan's when Fcell is 0 for whatever reason
            if self.calc_curvatures:
                self.fcell_second_deriv = SG*self._extracted_fcell_second_deriv[self.roi_slice]

    def _get_per_spot_scale(self, i_shot, i_spot):
        val = 1
        if self.refine_per_spot_scale:
            assert self.rescale_params

            xpos = self.per_spot_scale_xpos[i_shot][i_spot]
            val = self.Xall[xpos]
            sig = self.per_spot_scale_sigma
            init = 1
            val = np_exp(sig * (val - 1)) * init

        return val

    def _extract_pixel_data(self):
        self.roi_slice = self.roi_ids == self._i_spot
        self.model_bragg_spots = self.scale_fac*(self._model_pix[self.roi_slice])
        self.model_bragg_spots *= self._get_per_spot_scale(self._i_shot, self._i_spot)
        self._extract_Umatrix_derivative_pixels()
        self._extract_sausage_derivs()
        self._extract_Bmatrix_derivative_pixels()
        self._extract_mosaic_parameter_m_derivative_pixels()
        self._extract_ncells_def_derivative_pixels()
        self._extract_detector_distance_derivative_pixels()
        self._extract_Fcell_derivative_pixels()
        self._extract_spectra_coefficient_derivatives()
        self._extract_panelRot_derivative_pixels()
        self._extract_panelXYZ_derivative_pixels()

    def _update_ucell(self):
        if self.rescale_params:
            pars = self._get_ucell_vars(self._i_shot)
        else:
            _s = slice(self.ucell_xstart[self._i_shot], self.ucell_xstart[self._i_shot] + self.n_ucell_param, 1)
            pars = list(self.Xall[_s])
        self.UCELL_MAN[self._i_shot].variables = pars
        self._send_ucell_gradients_to_derivative_managers()
        self.D.Bmatrix = self.UCELL_MAN[self._i_shot].B_recipspace

    def _update_umatrix(self):
        self.D.Umatrix = self.CRYSTAL_MODELS[self._i_shot].get_U()

    def _update_beams(self):
        # sim_data instance has a nanoBragg beam object, which takes spectra and converts to nanoBragg xray_beams
        self.S.beam.spectrum = self.SPECTRA[self._i_shot]
        if self.verbose and self._i_shot == 0:
            self.print("Using %d BEAMZ!" % len(self.SPECTRA[self._i_shot]))
            self.print("Using %d mosaicities" % (self.S.D.mosaic_domains))
            self.print("Using %dx%d oversample rate" % (self.S.D.oversample, self.S.D.oversample))
        self.D.xray_beams = self.S.beam.xray_beams
        assert len(self.D.xray_beams) == len(self.SPECTRA[self._i_shot])

    def _get_panels_fasts_slows(self):
        npan = len(self.S.detector)
        nfast, nslow = self.S.detector[0].get_image_size()
        MASK = NP_ZEROS((npan, nslow, nfast), bool)
        #mx = 65535
        ROI_ID = NP_ZEROS((npan, nslow, nfast), 'uint16')
        #ROI_ID = NP_ONES((npan, nslow, nfast), 'uint16') * mx
        nspots = len(self.NANOBRAGG_ROIS[self._i_shot])
        #if nspots > mx:
        #    raise ValueError("Nspots is too large")
        #DATA = NP_ZEROS((npan, nslow, nfast), "float32")
        #import numpy as np
        #import h5py
        #h = h5py.File("1_again/imgs/rank0/stg1_run796_shot655_indexed_0.expt_0.h5", "r")
        #Z =np.abs(h['Z_data_noise'][()])
        #pids = h['pids'][()]
        #ys = h['ys'][()]
        #xs = h['xs'][()]
        #Zimg = np.zeros((256,254,254))
        #Zimg[pids, ys,xs] = Z

        for i_spot in range(nspots):
            (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
            if x2-x1==0 or y2-y1==0:
                continue

            saving_model = self.save_model_for_shot is not None and self.save_model_for_shot == self._i_shot
            if not saving_model and self.selection_flags is not None and not self.selection_flags[self._i_shot][i_spot]:
                continue
            pid = self.PANEL_IDS[self._i_shot][i_spot]
            MASK[pid, y1:y2, x1:x2] = True
            ROI_ID[pid, y1:y2, x1:x2] = i_spot
            #zbad = Zimg[pid, y1:y2, x1:x2]
            #if np.any(zbad > 4):
            #    from IPython import embed
            #    embed()
            #DATA[pid, y1:y2, x1:x2] = self.ROI_IMGS[self._i_shot][i_spot]
        p,s,f = WHERE(MASK)
        roi_id = ROI_ID[p,s,f]
        #from IPython import embed
        #embed()
        #data = DATA[p,s,f]
        data = None
        pan_fast_slow = np_ascontiguousarray((np_vstack([p,f,s]).T).ravel())
        pan_fast_slow = flex.size_t(pan_fast_slow)

        return pan_fast_slow, roi_id, data

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
        if self.calc_func:
            if self.verbose:
                self._print_iteration_header()

            if not self.only_save_model_for_shot:
                self._MPI_save_state_of_refiner()

            if self.iteratively_freeze_parameters:
                if self.param_sels is None:
                    self.determine_parameter_freeze_order()

            # reset gradient and functional
            self.target_functional = 0
            self._update_Xall_from_x()

            self.grad = flex_double(self.n_total_params)
            if self.calc_curvatures:
                self.curv = flex_double(self.n_total_params)

            # current work has gain_fac at 1 (#TODO gain factor should effect the probability of observing the photons)
            self.gain_fac = 1 #self.Xall[self.gain_xpos]
            self.G2 = self.gain_fac ** 2

            LOGGER.info("start update Fcell")
            self._store_updated_Fcell()
            self._update_Fcell()  # update the structure factor with the new x
            LOGGER.info("done update Fcell")
            self._update_spectra_coefficients()  # updates the diffBragg lambda coefficients if refinining spectra

            if self.CRYSTAL_GT is not None and not self.only_save_model_for_shot:
                self._MPI_initialize_GT_crystal_misorientation_analysis()
            tshots = time.time()

            LOGGER.info("Iterate over %d shots" % len(self.shot_ids))
            self._shot_Zscores = []
            for self._i_shot in self.shot_ids:

                if self.save_model_for_shot is not None:
                    if self.save_model_for_shot == self._i_shot:
                        self.full_image_of_model = NP_ZEROS(self.image_shape)
                        self.full_image_of_spots = NP_ZEROS(self.image_shape)
                        self.full_image_of_sigma_r = NP_ZEROS(self.image_shape)
                    elif self.only_save_model_for_shot:
                        continue

                if self.save_model:
                    self._open_model_output_file()

                if self._i_shot in self.bad_shot_list:
                    continue
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
                n_spots = len(self.NANOBRAGG_ROIS[self._i_shot])
                printed_geom_updates = False

                # CREATE THE PANEL FAST SLOW ARRAY AND RUN DIFFBRAGG
                LOGGER.info("run diffBragg for shot %d" % self._i_shot)
                self._run_diffBragg_current()
                LOGGER.info("finished diffBragg for shot %d" % self._i_shot)

                # CHECK FOR SIGNAL INTERRUPT HERE
                if self.break_signal is not None:
                    signal.signal(self.break_signal, self._sig_hand.handle)
                    self._MPI_check_for_break_signal()

                # TODO pre-extractions for all parameters
                self._append_local_parameters()
                #LOGGER.info("access arrays from diffBragg for shot %d" % self._i_shot)
                self._pre_extract_deriv_arrays()
                #self._per_shot_Z_data = []
                #LOGGER.info("iterate over %d spots and compute gradients for shot %d" % (n_spots, self._i_shot))
                self._spot_Zscores = []
                for i_spot in range(n_spots):
                    self._i_spot = i_spot
                    (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
                    if x2 - x1 == 0 or y2 - y1 == 0:
                        continue

                    saving_model = self.save_model_for_shot is not None and self.save_model_for_shot == self._i_shot
                    #if not saving_model and self.selection_flags is not None:
                    if self.selection_flags is not None:
                        if self._i_shot not in self.selection_flags:
                            continue
                        elif not self.selection_flags[self._i_shot][i_spot]:
                            continue

                    self._panel_id = int(self.PANEL_IDS[self._i_shot][i_spot])

                    if self.verbose and i_spot % self.spot_print_stride == 0:
                        self.print("diffBragg: img %d/%d; spot %d/%d; panel %d" \
                              % (self._i_shot + 1, self.n_shots, i_spot + 1, n_spots, self._panel_id))

                    self.Imeas = self.ROI_IMGS[self._i_shot][i_spot].ravel()
                    self._set_roi_mask()
                    self.sigma_r = self._update_sigma_r_if_possible()

                    if not printed_geom_updates and self._i_shot == 0:
                        if self.refine_panelRotO or self.refine_panelRotF or self.refine_panelRotS:
                            self.print("ROT ANGLES OFS : %.8f %.8f %.8f (degrees)"
                                  % tuple(map(lambda x: x*180/PI, self._get_panelRot_val(self._panel_id))))

                        if self.refine_panelXY or self.refine_panelZ:
                            self.print("XYZ: %.8f %8f %.8f mm" % tuple([val*1000 for val in self._get_panelXYZ_val(self._panel_id)]))
                        if self.refine_spectra:
                            self.print("Spectra Lam0 Lam1: ", self._get_spectra_coefficients())
                        printed_geom_updates = True
                    #self._run_diffBragg_current(i_spot)
                    self._set_background_plane(i_spot)
                    self._extract_pixel_data()
                    self._record_xy_calc()
                    self._evaluate_averageI()

                    if self.save_model:
                        self._save_model_to_disk()
                    if self.save_model_for_shot is not None and self.save_model_for_shot == self._i_shot:
                        self._update_full_image_of_model()

                    # here we can correlate modelLambda with Imeas
                    self._increment_model_data_correlation()

                    if self.poisson_only:
                        self._evaluate_log_averageI()
                    else:
                        self._evaluate_log_averageI_plus_sigma_readout()

                    #self._max_h_sanity_test()
                    self._derivative_convenience_factors()

                    #if self.save_Z_freq is not None:
                    #   self._per_shot_Z_data.append((self.u_times_one_over_v, self._panel_id, self.NANOBRAGG_ROIS[self._i_shot]))
                    sigZ = (self._Zscore[self._is_trusted] ).std()
                    miller_idx = self.ASU[self._i_shot][i_spot]
                    i_fcell = self.idx_from_asu[miller_idx]
                    self._spot_Zscores.append((i_fcell, sigZ))

                    self.target_functional += self._target_accumulate()

                    # make any plots (this only matters if proper flags have been set)
                    self._show_plots(i_spot, n_spots)

                    # accumulate the per pixel derivatives
                    if self.bg_extracted:
                        self._bg_extracted_derivatives(i_spot)
                    else:
                        self._background_derivatives(i_spot)
                    self._Umatrix_derivatives()
                    self._Bmatrix_derivatives()
                    self._mosaic_parameter_m_derivatives()
                    self._ncells_def_derivatives()
                    self._detector_distance_derivatives()
                    self._panelRot_derivatives()
                    self._panelXYZ_derivatives()
                    self._spot_scale_derivatives()
                    self._eta_derivatives()
                    self._per_spot_scale_derivatives()
                    self._gain_factor_derivatives()
                    self._Fcell_derivatives(i_spot)
                    self._spectra_derivatives()
                    self._sausage_derivatives()
                    # Done with derivative accumulation
                self._shot_Zscores.append(self._spot_Zscores)
            #    self.image_corr[self._i_shot] = self.image_corr[self._i_shot] / self.image_corr_norm[self._i_shot]
            tshots =time.time()-tshots
            LOGGER.info("Time rank worked on shots=%.4f" % tshots)
            self._append_global_parameters()
            if not self.only_save_model_for_shot:
                self._MPI_aggregate_model_data_correlations()
            # TODO add in the priors:
            self._priors()
            self._parameter_freezes()
            #if self.save_Z_freq is not None:
            #    self._save_Z()
            if not self.only_save_model_for_shot:
                tmpi = time.time()
                LOGGER.info("MPI aggregation of func and grad")
                self._mpi_aggregation()
                tmpi = tmpi-time.time()
                LOGGER.info("Time for MPIaggregation=%.4f" % tmpi)

                self._f = self.target_functional
                self._g = self.g_for_lbfgs
                self.g = self.g_for_lbfgs  # TODO why all these repeated definitions ?, self.g is needed by _verify_diag

                self._curvature_analysis()

            # reset ROI pixels TODO: is this necessary
            self.D.raw_pixels_roi *= 0
            if self.only_save_model_for_shot:
                return
            self._sanity_check_grad()
            self.gnorm = norm(self.grad)

            if self.verbose:
                # TODO, some of this is becoming crucial, so move out of verobse
                if self.CRYSTAL_GT is not None:
                    self._MPI_print_GT_crystal_misorientation_analysis()
                self._print_image_correlation_analysis()
                self.print_step()
                self.print_step_grads()

            tsave = time.time()
            LOGGER.info("DUMP param and Zscore data")
            self._dump_parameters_to_hdf5()
            self._save_Zscore_data()
            tsave = time.time()-tsave
            LOGGER.info("Time to dump param and Zscore data: %.4f" % tsave)

            self.iterations += 1
            self.f_vals.append(self.target_functional)
            #time.sleep(self.pause_after_iteration)

            if self.calc_curvatures and not self.use_curvatures:
                saving_model = self.save_model_for_shot is not None and self.save_model_for_shot == self._i_shot
                if not saving_model and self.num_positive_curvatures == self.use_curvatures_threshold:
                    raise BreakToUseCurvatures

        if self.save_model:
            self._close_model_output_file()

        LOGGER.info("DONE WITH FUNC GRAD")
        return self._f, self._g

    def _save_Zscore_data(self):
        if not self.iterations % self.saveZ_freq == 0:
            return
        outdir = PATHJOIN(self.output_dir, "rank%d_Zscore" % self.rank)
        if not EXISTS(outdir):
            MAKEDIRS(outdir)
        #i_fcell, sigZ = zip(*self._shot_Zscores)
        fname = PATHJOIN(outdir, "sigZ_iter%d_rank%d" % (self.iterations, self.rank))
        SAVE(fname, self._shot_Zscores)

    def _sanity_check_grad(self):
        if not self.refine_background_planes:
            for xp in self.bg_a_xstart[self._i_shot]:
                assert self.grad[xp] == 0
            for xp in self.bg_b_xstart[self._i_shot]:
                assert self.grad[xp] == 0
            for xp in self.bg_c_xstart[self._i_shot]:
                assert self.grad[xp] == 0
        if not self.refine_blueSausages:
            for i_saus in range(self.num_sausages):
                for i_saus_param in range(4):
                    idx = i_saus * 4 + i_saus_param
                    assert self.grad[self.sausages_xpos[self._i_shot][idx]] == 0
        if not self.refine_eta:
            for i_eta in range(3):
                assert self.grad[self.eta_xstart[self._i_shot]+i_eta] == 0
                if i_eta in [1, 2] and not self.S.crystal.has_anisotropic_mosaicity:
                    assert self.grad[self.eta_xstart[self._i_shot] + i_eta] == 0

        if not self.refine_spectra:
            for i_p in range(2):
                assert self.grad[self.spectra_coef_xstart + i_p] == 0
        if not self.refine_ncells:
            for i_ncells in range(self.n_ncells_param):
                assert self.grad[self.ncells_xstart[self._i_shot] + i_ncells] == 0
        if not self.refine_ncells_def:
            for i_ncells in range(3):
                assert self.grad[self.ncells_def_xstart[self._i_shot] + i_ncells] == 0

        if not self.refine_Umatrix:
            # if not self.refine_rotX:
            assert self.grad[self.rotX_xpos[self._i_shot]] == 0
            # if not self.refine_rotY:
            assert self.grad[self.rotY_xpos[self._i_shot]] == 0
            # if not self.refine_rotZ:
            assert self.grad[self.rotZ_xpos[self._i_shot]] == 0

        if not self.refine_Bmatrix:
            for i_p in range(self.n_ucell_param):
                assert self.grad[self.ucell_xstart[self._i_shot] + i_p] == 0

        if not self.refine_crystal_scale:
            assert self.grad[self.spot_scale_xpos[self._i_shot]] == 0

        for pg in set(list(self.panel_group_from_id.values())):
            if not self.refine_panelRotO:
                assert self.grad[self.panelRot_xstart + 3 * pg] == 0
            if not self.refine_panelRotF:
                assert self.grad[self.panelRot_xstart + 3 * pg + 1] == 0
            if not self.refine_panelRotS:
                assert self.grad[self.panelRot_xstart + 3 * pg + 2] == 0

            if not self.refine_panelXY:
                assert self.grad[self.panelXY_xstart + 2 * pg] == 0
                assert self.grad[self.panelXY_xstart + 2 * pg + 1] == 0
            if not self.refine_panelZ:
                assert self.grad[self.panelZ_xstart] == 0

        if not self.refine_detdist:
            assert self.grad[self.detector_distance_xpos[self._i_shot]] == 0

        if not self.refine_per_spot_scale:
            for xp in self.per_spot_scale_xpos[self._i_shot]:
                assert self.grad[xp] == 0

    def _set_roi_mask(self):
        if self.MASK is not None:
            (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][self._i_spot]
            self._is_trusted = self.MASK[self._panel_id][y1:y2, x1:x2].ravel()

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

    def _save_model_to_disk(self):

        if not self.I_AM_ROOT:
           return
        from time import time
        t = time()
        from numpy import nan_to_num
        #from numpy import str_ as np_str_
        img_file = self.FNAMES[self._i_shot]
        proc_name = self.PROC_FNAMES[self._i_shot]
        proc_idx = self.PROC_IDX[self._i_shot]
        (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][self._i_spot]
        bbox = x1, x2, y1, y2
        miller_idx_asu = self.ASU[self._i_shot][self._i_spot]
        miller_idx = self.Hi[self._i_shot][self._i_spot]
        i_fcell = self.idx_from_asu[self.ASU[self._i_shot][self._i_spot]]
        Famp = self._get_fcell_val(i_fcell)
        tilt = self.tilt_plane
        bbox_idx = self.BBOX_IDX[self._i_shot][self._i_spot]
        bragg = self.model_bragg_spots
        hi_tag = "%+04d%+04d%+04d" % miller_idx
        key = "%s-%s" % (img_file, hi_tag)
        key = key.replace("/", "+-+")
        self.model_hout.create_dataset("%s/data" % key, data=self.Imeas)
        self.model_hout.create_dataset("%s/tilt" % key, data=tilt)
        self.model_hout.create_dataset("%s/bragg" % key, data=bragg)
        self.model_hout.create_dataset("%s/partial" % key, data=nan_to_num(bragg/Famp/Famp))
        self.model_hout.create_dataset("%s/Famp" % key, data=[Famp])
        self.model_hout.create_dataset("%s/bbox" % key, data=bbox)
        self.model_hout.create_dataset("%s/miller_idx" % key, data=miller_idx)
        self.model_hout.create_dataset("%s/miller_idx_asu" % key, data=miller_idx_asu)
        self.model_hout.create_dataset("%s/proc_name" % key, data=ARRAY(proc_name, "S"))
        self.model_hout.create_dataset("%s/proc_idx" % key, data=[proc_idx])
        self.model_hout.create_dataset("%s/bbox_idx" % key, data=[bbox_idx])
        self.model_hout.create_dataset("%s/panel" % key, data=[self._panel_id])
        t = time()-t
        self.save_time += t

    def _increment_model_data_correlation(self):
        if self.image_corr[self._i_shot] is None:
            self.image_corr[self._i_shot] = 0
            self.image_corr_norm[self._i_shot] = 0
        if self.compute_image_model_correlation:
            _overlay_corr, _ = pearsonr(self.Imeas.ravel(), self.model_Lambda.ravel())
        else:
            _overlay_corr = NAN
        self.image_corr[self._i_shot] += _overlay_corr
        self.image_corr_norm[self._i_shot] += 1

    def _background_derivatives(self, i_spot):
        if self.refine_background_planes:
            if self.bg_offset_only:  # option to only refine c (t3 in manuscript) plane
                abc_dI_dtheta = [0, 0, self.G2*self.c]
                abc_d2I_dtheta2 = [0, 0, 0]
            else:
                xr = self.XREL[self._i_shot][i_spot].ravel()  # fast scan pixels
                yr = self.YREL[self._i_shot][i_spot].ravel()  # slow scan pixels
                if self.G2 != 1:
                    abc_dI_dtheta = [xr*self.G2, yr*self.G2, self.G2]
                else:
                    abc_dI_dtheta = [xr, yr, self.G2]
                abc_d2I_dtheta2 = [0, 0, 0]

            if self.rescale_params or self.bg_offset_positive:
                if self.bg_offset_only:
                    # here we apply case2 type reparameterization to the c derivative
                    abc_dI_dx = [0, 0, abc_dI_dtheta[2]*self.c*self.c_sigma]
                    abc_d2I_dx2 = [0, 0, abc_dI_dx[2]*self.c_sigma]
                else:
                    # here we apply case1 type reparameterization to a,b,c derivs
                    abc_dI_dx = [abc_dI_dtheta[0]*self.a_sigma,
                                 abc_dI_dtheta[1]*self.b_sigma,
                                 abc_dI_dtheta[2]*self.c_sigma]
                    abc_d2I_dx2 = [0, 0, 0]
                bg_deriv = abc_dI_dx
                bg_second_deriv = abc_d2I_dx2
            else:
                bg_deriv = abc_dI_dtheta
                bg_second_deriv = abc_d2I_dtheta2

            x_positions = [self.bg_a_xstart[self._i_shot][i_spot],
                           self.bg_b_xstart[self._i_shot][i_spot],
                           self.bg_c_xstart[self._i_shot][i_spot]]
            for ii, xpos in enumerate(x_positions):
                d = bg_deriv[ii]
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2 = bg_second_deriv[ii]
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _get_rotXYZ_first_derivs(self):
        if self.rescale_params:
            # rot XYZ uses a case1 type rescaling
            derivs = [self.rotX_sigma*self.rotX_dI_dtheta,
                      self.rotY_sigma*self.rotY_dI_dtheta,
                      self.rotZ_sigma*self.rotZ_dI_dtheta]
        else:
            derivs = [self.rotX_dI_dtheta, self.rotY_dI_dtheta, self.rotZ_dI_dtheta]

        return derivs

    def _get_rotXYZ_second_derivs(self):
        if self.rescale_params:
            # rot XYZ uses a case1 type rescaling
            second_derivs = [(self.rotX_sigma**2)*self.rotX_d2I_dtheta2,
                             (self.rotY_sigma**2)*self.rotY_d2I_dtheta2,
                             (self.rotZ_sigma**2)*self.rotZ_d2I_dtheta2]
        else:
            second_derivs = [self.rotX_d2I_dtheta2, self.rotY_d2I_dtheta2, self.rotZ_d2I_dtheta2]

        return second_derivs

    def _get_spectra_first_derivs(self):
        if self.refine_lambda0 or self.refine_lambda1:
            assert self.rescale_params
            derivs = []
            for i_coef in range(self.n_spectra_param):
                dI_dtheta = self.spectra_derivs[i_coef]

                # TODO make this part of a `parameter class`, so other parameters can borrow same code when using bounds
                init = self.spectra_coefficients_init[i_coef]
                sigma = self.spectra_coefficients_sigma[i_coef]
                x = self.Xall[self.spectra_coef_xstart + i_coef]
                low, high = self.lambda_coef_ranges[i_coef]
                rng = high - low
                cos_arg = sigma * (x-1) + ASIN(2*(init-low)/rng - 1)
                dtheta_dx = rng/2 * COS(cos_arg) * sigma
                d = dI_dtheta * dtheta_dx
                derivs.append(d)
                #derivs.append(dI_dtheta * sigma)
            return derivs

    def _spectra_derivatives(self):
        if self.refine_lambda0 or self.refine_lambda1:
            derivs = self._get_spectra_first_derivs()
            xstart = self.spectra_coef_xstart
            for i_coef in range(self.n_spectra_param):
                if not self.is_being_refined[xstart+i_coef]:
                    continue
                d = derivs[i_coef]
                self.grad[xstart + i_coef] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    raise NotImplementedError

    def _sausage_derivatives(self):
        if self.refine_blueSausages:
            for i_sausage in range(self.num_sausages):
                for i_param in range(4):
                    idx = i_sausage*4 + i_param
                    xpos = self.sausages_xpos[self._i_shot][idx]
                    d = self.sausages_dI_dtheta[idx] * self.sausages_sigma[i_param]
                    self.grad[xpos] += self._grad_accumulate(d)

    def _Umatrix_derivatives(self):
        if self.refine_Umatrix:
            x_positions = [self.rotX_xpos[self._i_shot],
                           self.rotY_xpos[self._i_shot],
                           self.rotZ_xpos[self._i_shot]]
            derivs = self._get_rotXYZ_first_derivs()
            second_derivs = self._get_rotXYZ_second_derivs()
            for ii, xpos in enumerate(x_positions):
                d = derivs[ii]
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2 = second_derivs[ii]
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _get_ucell_first_derivatives(self):
        derivs = []
        for i_ucell in range(self.n_ucell_param):
            d = self.ucell_dI_dtheta[i_ucell]
            if self.rescale_params:
                if self.use_ucell_ranges:
                    xpos = self.ucell_xstart[self._i_shot] + i_ucell
                    xval = self.Xall[xpos]
                    d = self.ucell_params[self._i_shot][i_ucell].get_deriv(xval, d)
                else:
                    sigma = self.ucell_sigmas[i_ucell]
                    d = d*sigma
            derivs.append(d)
        return derivs

    def _get_ucell_second_derivatives(self):
        second_derivs = []
        for i_ucell in range(self.n_ucell_param):
            d2 = self.ucell_d2I_dtheta2[i_ucell]
            if self.rescale_params:
                if self.use_ucell_ranges:
                    d = self.ucell_dI_dtheta[i_ucell]
                    xpos = self.ucell_xstart[self._i_shot] + i_ucell
                    xval = self.Xall[xpos]
                    d2 = self.ucell_params[self._i_shot][i_ucell].get_second_deriv(xval, d, d2)
                else:
                    sigma_squared = self.ucell_sigmas[i_ucell]**2
                    d2 = d2*sigma_squared
            second_derivs.append(d2)
        return second_derivs

    def _Bmatrix_derivatives(self):
        if self.refine_Bmatrix:
            # unit cell derivative
            derivs = self._get_ucell_first_derivatives()
            second_derivs = self._get_ucell_second_derivatives()#first_deriv=derivs)
            for i_ucell in range(self.n_ucell_param):
                xpos = self.ucell_xstart[self._i_shot] + i_ucell
                d = derivs[i_ucell]
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2 = second_derivs[i_ucell]
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _ncells_def_derivatives(self):
        if self.refine_ncells_def:
            for i_ncell in range(3):
                xpos = self.ncells_def_xstart[self._i_shot] + i_ncell
                xval = self.Xall[xpos]
                d = self.Ncells_def_params[self._i_shot][i_ncell].get_deriv(xval, self.ncells_def_dI_dtheta[i_ncell])
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2 = self.Ncells_def_params[self._i_shot][i_ncell].get_second_deriv(xval,
                                                                                self.ncells_def_dI_dtheta[i_ncell],
                                                                                self.ncells_def_d2I_dtheta2[i_ncell])
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _mosaic_parameter_m_derivatives(self):
        # NOTE This function can be made a lot simpler by dropping support for the case where two Ncells parameters are
        #      constrained to be equal..
        if self.refine_ncells:
            thetas = self._get_m_val(self._i_shot)  # mosaic parameters "m"
            if self.ncells_mask is not None:
                temp_thetas = [0, 0]
                for i_mask, mask_val in enumerate(self.ncells_mask):
                    temp_thetas[mask_val] = thetas[i_mask]
                thetas = temp_thetas
            for i_ncell in range(self.n_ncells_param):
                if self.use_Ncells_range:
                    xpos = self.ncells_xstart[self._i_shot] + i_ncell
                    xval = self.Xall[xpos]
                    d = self.Ncells_params[self._i_shot][i_ncell].get_deriv(xval, self.m_dI_dtheta[i_ncell])
                    d2 = 0
                    if self.calc_curvatures:
                        d2 = self.Ncells_params[self._i_shot][i_ncell].get_second_deriv(xval,
                            self.m_dI_dtheta[i_ncell], self.m_d2I_dtheta2[i_ncell])
                else:
                    theta_minus_three = thetas[i_ncell] - 3
                    try:
                        sig = self.m_sigma[i_ncell]
                    except TypeError:
                        sig = self.m_sigma
                    if self.rescale_params:
                        # case 3 rescaling
                        try:
                            sig = sig[0]
                        except TypeError:
                            pass
                        sig_theta_minus_three = sig*theta_minus_three
                        d = self.m_dI_dtheta[i_ncell]*sig_theta_minus_three
                        d2 = self.m_d2I_dtheta2[i_ncell]*(sig_theta_minus_three*sig_theta_minus_three) + \
                             self.m_dI_dtheta[i_ncell]*(sig*sig_theta_minus_three)
                    else:
                        # case 4 rescaling with theta_o = 3
                        d = self.m_dI_dtheta[i_ncell]*theta_minus_three
                        d2 = self.m_d2I_dtheta2[i_ncell]*(theta_minus_three*theta_minus_three) + self.m_dI_dtheta[i_ncell]*theta_minus_three

                xpos = self.ncells_xstart[self._i_shot] + i_ncell
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _panelRot_derivatives(self):
        refining = self.refine_panelRotO, self.refine_panelRotF, self.refine_panelRotS
        for i_rot in range(3):
            if refining[i_rot]:
                d = self.panelRot_dI_dtheta[i_rot]
                d2 = self.panelRot_d2I_dtheta2[i_rot]

                panel_group_id = self.panel_group_from_id[self._panel_id]
                xpos = self.panelRot_xstart + 3*panel_group_id + i_rot
                xval = self.Xall[xpos]
                d = self.panelRot_params[panel_group_id][i_rot].get_deriv(xval, d) ## apply chain rule for reparameterization

                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    raise NotImplementedError("curvatures and panel rotations not supported")
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _panelXYZ_derivatives(self):
        if self.refine_panelXY or self.refine_panelZ:
            panel_group_id = self.panel_group_from_id[self._panel_id]

            d_X = self.panelX_dI_dtheta
            d_Y = self.panelY_dI_dtheta
            d_Z = self.panelZ_dI_dtheta
            #d2_X = d2_Y = d2_Z = 0 #TODO: curvatires

            xpos_X = self.panelXY_xstart + 2*panel_group_id
            xpos_Y = xpos_X + 1
            xpos_Z = self.panelZ_xstart + panel_group_id

            if self.refine_panelXY:
                d_X = self.panelX_params[panel_group_id].get_deriv(self.Xall[xpos_X], d_X)
                d_Y = self.panelY_params[panel_group_id].get_deriv(self.Xall[xpos_Y], d_Y)
                self.grad[xpos_X] += self._grad_accumulate(d_X)
                self.grad[xpos_Y] += self._grad_accumulate(d_Y)
            if self.refine_panelZ:
                d_Z = self.panelZ_params[panel_group_id].get_deriv(self.Xall[xpos_Z], d_Z)
                self.grad[xpos_Z] += self._grad_accumulate(d_Z)

    def _detector_distance_derivatives(self):
        if self.refine_detdist:
            if self.rescale_params:
                if self.detector_distance_range is not None:
                    init = self.shot_detector_distance_init[self._i_shot]
                    sigma = self.detector_distance_sigma
                    x = self.Xall[self.detector_distance_xpos[self._i_shot]]
                    low, high = self.detector_distance_range
                    rng = high - low
                    cos_arg = sigma * (x - 1) + ASIN(2 * (init - low) / rng - 1)
                    dtheta_dx = rng / 2 * COS(cos_arg) * sigma
                    d = self.detdist_dI_dtheta * dtheta_dx
                    d2 = 0  #TODO implement..

                else:
                    #NOTE old way:
                    # case 1 type of rescaling
                    d = self.detdist_dI_dtheta*self.detector_distance_sigma
                    d2 = self.detdist_d2I_dtheta2*(self.detector_distance_sigma*self.detector_distance_sigma)
            else:
                d = 1*self.detdist_dI_dtheta
                d2 = 1*self.detdist_d2I_dtheta2

            xpos = self.detector_distance_xpos[self._i_shot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                if self.detector_distance_range is not None and self.rescale_params:
                    raise NotImplementedError("You Cannot.")
                self.curv[xpos] += self._curv_accumulate(d, d2)

    def _Fcell_derivatives(self, i_spot):
        if not self.refine_Fcell:
            return
        # asu index
        miller_idx = self.ASU[self._i_shot][i_spot]
        # get multiplicity of this index
        multi = self.hkl_frequency[self.idx_from_asu[miller_idx]]
        # check if we are freezing this index during refinement
        freeze_this_hkl = False
        if self.freeze_idx is not None:
            freeze_this_hkl = self.freeze_idx[miller_idx]

        # do the derivative
        if self.refine_Fcell and multi >= self.min_multiplicity and not freeze_this_hkl:
            i_fcell = self.idx_from_asu[self.ASU[self._i_shot][i_spot]]
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

    def _bg_extracted_derivatives(self, return_derivatives=False):
        if self.refine_background_planes:
            assert self.rescale_params
            dI_dtheta = self.G2 * self.tilt_plane
            sig = self.bg_coef_sigma
            # case 2 type rescaling
            d = dI_dtheta*self.bg_coef * sig
            d2 = dI_dtheta*(self.bg_coef*sig*sig)

            xpos = self.bg_coef_xpos[self._i_shot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)
            if return_derivatives:
                return d, d2

    def _spot_scale_derivatives(self, return_derivatives=False):
        if self.refine_crystal_scale:
            S = SQRT(self.scale_fac)
            dI_dtheta = (2*self.G2/S)*self.model_bragg_spots
            d2I_dtheta2 = (2*self.G2/S/S)*self.model_bragg_spots
            # second derivative is 0 with respect to scale factor
            if self.rescale_params:
                sig = self.spot_scale_sigma
                # case 2 type rescaling
                #d = dI_dtheta*self.scale_fac * sig
                d = dI_dtheta*sig
                d2 = d2I_dtheta2 *(sig**2)
            else:
                # case 4 type rescaling
                d = dI_dtheta*self.scale_fac
                d2 = d  # same as first derivative

            xpos = self.spot_scale_xpos[self._i_shot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)
        if return_derivatives:
            return d, d2

    def _eta_derivatives(self):
        if self.refine_eta:
            for i_eta in range(3):
                xpos = self.eta_xstart[self._i_shot] + i_eta
                dI_dtheta = self.G2*self.scale_fac*self._eta_deriv[i_eta][self.roi_slice]
                xval = self.Xall[xpos]
                d = self.eta_params[self._i_shot][i_eta].get_deriv(xval, dI_dtheta)
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2I_dtheta2 = self.G2 * self.scale_fac * self._eta_second_deriv[i_eta][self.roi_slice]
                    d2 = self.eta_params[self._i_shot][i_eta].get_second_deriv(xval, dI_dtheta, d2I_dtheta2)
                    self.curv[xpos] += self._curv_accumulate(d, d2)
                if not self.S.crystal.has_anisotropic_mosaicity:
                    break

    def _per_spot_scale_derivatives(self):
        if self.refine_per_spot_scale:
            assert self.rescale_params
            per_spot_scale = self._get_per_spot_scale(self._i_shot, self._i_spot)
            dI_dtheta = self.G2 * self.model_bragg_spots / per_spot_scale
            # second derivative is 0 with respect to scale factor
            sig = self.per_spot_scale_sigma
            # case 2 type rescaling
            d = dI_dtheta*self.scale_fac * sig
            d2 = dI_dtheta*(self.scale_fac*sig*sig)

            xpos = self.per_spot_scale_xpos[self._i_shot][self._i_spot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)

    def _gain_factor_derivatives(self):
        if self.refine_gain_fac:
            raise NotImplementedError("gain factor derivatives need more testing")
            d = 2*self.gain_fac*(self.tilt_plane + self.model_bragg_spots)
            self.grad[self.gain_xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                d2 = d / self.gain_fac
                self.curv[self.gain_xpos] += self._curv_accumulate(d, d2)

    def _max_h_sanity_test(self, i_spot):
        return
        max_h = tuple(map(int, self.D.max_I_hkl))
        refinement_h = self.ASU[self._i_shot][i_spot]
        equivs = [i.h() for i in miller.sym_equiv_indices(self.space_group, refinement_h).indices()]
        if not max_h in equivs and self.debug:  # TODO understand this more, how does this effect things
            self.print("Warning max_h  mismatch!!!!!!")

    def _priors(self):
        # experimental, not yet proven to help
        if self.use_ucell_priors and self.refine_Bmatrix:
            for ii in range(self.n_shots):
                for jj in range(self.n_ucell_param):
                    xpos = self.ucell_xstart[ii] + jj
                    ucell_p = self.Xall[xpos]
                    sig_square = self.sig_ucell[jj] ** 2
                    self.target_functional += (ucell_p - self.ave_ucell[jj]) ** 2 / 2 / sig_square
                    self.grad[xpos] += (ucell_p - self.ave_ucell[jj]) / sig_square
                    if self.calc_curvatures:
                        self.curv[xpos] += 1 / sig_square

        if self.use_rot_priors and self.refine_Umatrix:
            for ii in range(self.n_shots):
                x_positions = [self.rotX_xpos[self._i_shot],
                               self.rotY_xpos[self._i_shot],
                               self.rotZ_xpos[self._i_shot]]
                for xpos in x_positions:
                    rot_p = self.Xall[xpos]
                    sig_square = self.sig_rot ** 2
                    self.target_functional += rot_p ** 2 / 2 / sig_square
                    self.grad[xpos] += rot_p / sig_square
                    if self.calc_curvatures:
                        self.curv[xpos] += 1 / sig_square

    def _parameter_freezes(self):
        if self.iteratively_freeze_parameters and self.iterations % self.number_of_frozen_iterations == 0:
            self.print("\n\n\t\tSwitching!!\n\n")
            freeze_sel = next(self.param_sels)
            self.grad.set_selected(freeze_sel, 0)
            if self.calc_curvatures:
                self.curv.set_selected(freeze_sel, 0)

    def _mpi_aggregation(self):
        # reduce the broadcast summed results:
        if self.I_AM_ROOT:
            self.print("\nMPI reduce on functionals and gradients...")
        LOGGER.info("Functional")
        self.target_functional = self._MPI_reduce_broadcast(self.target_functional)
        LOGGER.info("gradients")
        self.grad = self._MPI_reduce_broadcast(self.grad)
        LOGGER.info("unpacks")
        self.rotx, self.roty, self.rotz, self.uc_vals, self.ncells_vals, self.scale_vals, \
        self.scale_vals_truths, self.origZ_vals = self._unpack_internal(self.Xall, lst_is_x=True)
        self.Grotx, self.Groty, self.Grotz, self.Guc_vals, self.Gncells_vals, self.Gscale_vals, _, self.GorigZ_vals = \
            self._unpack_internal(self.grad, lst_is_x=False)
        if self.calc_curvatures:
            self.curv = self._MPI_reduce_broadcast(self.curv)
            self.CUrotx, self.CUroty, self.CUrotz, self.CUuc_vals, self.CUncells_vals, self.CUscale_vals, _, self.CUorigZ_vals = \
                self._unpack_internal(self.curv, lst_is_x=False)
        LOGGER.info("negative fcells")
        self.tot_fcell_negative_model = self._MPI_reduce_broadcast(self.num_Fcell_negative_model)

    def _curvature_analysis(self):
        self.tot_neg_curv = 0
        self.neg_curv_shots = []
        if self.calc_curvatures:
            self.is_negative_curvature = self.curv.as_numpy_array() < 0
            self.tot_neg_curv = sum(self.is_negative_curvature)

        if self.calc_curvatures and not self.use_curvatures:
            if self.tot_neg_curv == 0:
                self.num_positive_curvatures += 1
                self.d = self.d_for_lbfgs #flex_double(self.curv.as_numpy_array())
                self._verify_diag()
            else:
                self.num_positive_curvatures = 0
                self.d = None

        if self.use_curvatures:
            if self.tot_neg_curv == 0:
                self.request_diag_once = False
                self.diag_mode = "always"  # TODO is this proper place to set ?
                self.d = self.d_for_lbfgs #flex_double(self.curv.as_numpy_array())
                self._verify_diag()
            elif self.fix_params_with_negative_curvature:
                self.request_diag_once = False
                self.diag_mode = "always"
                is_ref = self.is_being_refined.as_numpy_array()
                is_ref[self.is_negative_curvature] = False
                self.is_being_refined = FLEX_BOOL(is_ref)
                # set the BFGS parameter array
                self.x = self.x_for_lbfgs
                assert( self.n == len(self.x))
                # make the mapping from x to Xall
                refine_pos = WHERE(self.is_being_refined.as_numpy_array())[0]
                self.x2xall = {xi: xalli for xi, xalli in enumerate(refine_pos)}
                self.xall2x = {xalli: xi for xi, xalli in enumerate(refine_pos)}
                self.g = self.g_for_lbfgs
                self.d = self.d_for_lbfgs
                self._g = self.g_for_lbfgs
                self._verify_diag()
                self.print("Breaking to freeze %d curvatures" % self.tot_neg_curv)
                raise BreakToUseCurvatures
            else:
                if self.debug:
                    self.print("\n\t******************************************")
                    self.print("\tFREEZING THE CURVATURE: DISASTER AVERSION")
                    self.print("*\t*****************************************")
        else:
            self.d = None

    def _MPI_initialize_GT_crystal_misorientation_analysis(self):
        self.all_ang_off = []
        self.current_ang_off = []
        for i in range(self.n_shots):
            try:
                Ctru = self.CRYSTAL_GT[i]
                atru, btru, ctru = Ctru.get_real_space_vectors()
                ang, ax = self.get_correction_misset(as_axis_angle_deg=True, i_shot=i)
                B = self.get_refined_Bmatrix(i)
                C = deepcopy(self.CRYSTAL_MODELS[i])
                C.set_B(B)
                if ang > 0:
                    C.rotate_around_origin(ax, ang)
                ang_off = compare_with_ground_truth(atru, btru, ctru,
                                                    [C],
                                                    symbol=self.symbol)[0]
            except Exception:
                ang_off = -1
            if self.filter_bad_shots and self.iterations == 0:
                if ang_off == -1 or ang_off > 0.015:
                    self.bad_shot_list.append(i)

            self.current_ang_off.append(ang_off)
            self.all_ang_off.append(ang_off)

        if self.init_ang_off is None:
            self.init_ang_off = deepcopy(self.current_ang_off)

        self.bad_shot_list = list(set(self.bad_shot_list))

        self.all_ang_off = self._init_gather_ang_off()
        self.n_bad_shots = self._init_n_bad_shots()

    def get_init_misorientation(self, i_shot):
        ang = NAN
        if self.CRYSTAL_GT is not None:
            ang = self.init_ang_off[i_shot]
        return ang

    def get_current_misorientation(self, i_shot):
        ang = -1 #NAN
        if self.CRYSTAL_GT is not None:
            ang = self.current_ang_off[i_shot]
        return ang

    def _MPI_print_GT_crystal_misorientation_analysis(self):
        all_ang_off = self._get_ang_off()
        self.print(all_ang_off)
        n_broken_misset = sum([1 for aa in all_ang_off if aa == -1])
        n_bad_misset = sum([1 for aa in all_ang_off if aa > 0.1])
        n_misset = len(all_ang_off)
        _pos_misset_vals = [aa for aa in all_ang_off if aa > 0]
        misset_median = median(_pos_misset_vals)
        misset_mean = mean(_pos_misset_vals)
        misset_max = -1
        misset_min = -1
        if _pos_misset_vals:
            misset_max = max(_pos_misset_vals)
            misset_min = min(_pos_misset_vals)
        if self.refine_Umatrix or self.refine_Bmatrix and self.print_all_missets:
            self.print("\nMissets\n========")
            all_ang_off_s = ["%.5f" % aa for aa in all_ang_off]
            self.print(", ".join(all_ang_off_s))
            self.print("N shots deemed bad from missets: %d" % self.n_bad_shots)
        self.print("MISSETTING median: %.4f; mean: %.4f, max: %.4f, min %.4f, num > .1 deg: %d/%d; num broken=%d"
              % (misset_median, misset_mean, misset_max, misset_min, n_bad_misset, n_misset, n_broken_misset))
        self.all_ang_off = all_ang_off

    def _get_image_correlation(self, i_shot):
        corr = NAN
        if self.image_corr[i_shot] is not None:
            corr = self.image_corr[i_shot] / self.image_corr_norm[i_shot]
        return corr

    def _get_init_image_correlation(self, i_shot):
        corr = NAN
        if self.image_corr[i_shot] is not None:
            corr = self.init_image_corr[i_shot]
        return corr

    def _print_image_correlation_analysis(self):
        return
        all_corr_str = ["%.2f" % ic for ic in self.all_image_corr]
        self.print("Correlation stats:")
        if self.print_all_corr:
            self.print(", ".join(all_corr_str))
            self.print("---------------")
        n_bad_corr = sum([1 for ic in self.all_image_corr if ic < 0.25])
        self.print("CORRELATION median: %.4f; mean: %.4f, max: %.4f, min %.4f, num <.25: %d/%d;"
              % (median(self.all_image_corr), mean(self.all_image_corr), max(self.all_image_corr),
                 min(self.all_image_corr), n_bad_corr, len(self.all_image_corr)))

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
            if self.rescale_params:
                fvals = [self._get_fcell_val(i_fcell) for i_fcell in range(self.n_global_fcell)]
                fvals = ARRAY(fvals)
            else:
                fvals = self.Xall[self.fcell_xstart:self.fcell_xstart + self.n_global_fcell].as_numpy_array()
            SAVEZ(outf, fvals=fvals, x=self.Xall.as_numpy_array())

    def _show_plots(self, i_spot, n_spots):
        if self.I_AM_ROOT and self.plot_images and self.iterations % self.plot_stride == 0 and self._i_shot == self.index_of_displayed_image:
            (x1, x2), (y1, y2) = self.NANOBRAGG_ROIS[self._i_shot][self._i_spot]
            img_sh = y2-y1, x2-x1
            if i_spot % self.plot_spot_stride == 0:
                xr = self.XREL[self._i_shot][i_spot]  # fast scan pixels
                yr = self.YREL[self._i_shot][i_spot]  # slow scan pixels
                if self.plot_residuals:
                    self.ax.clear()
                    residual = self.model_Lambda - self.Imeas
                    x = residual.max()
                    #else:
                    #    x = mean([x, residual.max()])
                    self.ax.plot_surface(xr, yr, residual.reshape(img_sh), rstride=2, cstride=2, alpha=0.3, cmap='coolwarm')
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
                    vmax = m + 5 * s
                    vmin = m - s
                    m2 = self.model_Lambda.mean()
                    s2 = self.model_Lambda.std()
                    self.ax1.images[0].set_data(self.model_Lambda.reshape(img_sh))
                    self.ax1.images[0].set_clim(vmin, vmax)
                    self.ax2.images[0].set_data(self.Imeas.reshape(img_sh))
                    self.ax2.images[0].set_clim(vmin, vmax)
                plt.suptitle("Iterations = %d, image %d / %d"
                             % (self.iterations, i_spot + 1, n_spots))
                self.fig.canvas.draw()
                plt.pause(.02)

    def _poisson_target(self):
        fterm = self.model_Lambda - self.Imeas * self.log_Lambda
        if self._is_trusted is not None:
            fterm = fterm[self._is_trusted]
        fterm = fterm.sum()
        return fterm

    def _poisson_d(self, d):
        if self.refine_with_psf:
            d = convolve_with_psf(d, psf=self._psf, **self.psf_args)
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
        if self.refine_with_psf:
            d = convolve_with_psf(d, psf=self._psf, **self.psf_args)
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
        self._Zscore = self.u*SQRT(self.one_over_v)

    def _evaluate_log_averageI(self):  # for Poisson only stats
        try:
            self.log_Lambda = np_log(self.model_Lambda)
        except FloatingPointError:
            pass
        if any((self.model_Lambda <= 0).ravel()):
            self.num_negative_model += 1
            is_bad = self.model_Lambda <= 0
            self.log_Lambda[is_bad] = 1e-6
            self.print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL (negative_models=%d)!!!!!!!!!\n<><><><><><><><><>\n" % self.num_negative_model)
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda[self.model_Lambda <= 0] = 0

    def _evaluate_log_averageI_plus_sigma_readout(self):
        v = self.model_Lambda + self.sigma_r ** 2
        v_is_neg = (v <= 0).ravel()
        if any(v_is_neg):
            self.num_negative_model += 1
            self.print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_v = np_log(v)
        self.log_v[v <= 0] = 0  # but will I ever negative_model ?

    def print_step(self):
        """Deprecated"""
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))

        if self.refine_Umatrix or self.refine_Bmatrix or self.refine_crystal_scale or self.refine_ncells:
            if self.big_dump and HAS_PANDAS:
                master_data = {"Ncells": self.ncells_vals,
                               "scale": self.scale_vals,
                               "rotx": self.rotx,
                               "roty": self.roty,
                               "rotz": self.rotz, "origZ": self.origZ_vals}
                for i_uc in range(self.n_ucell_param):
                    master_data["uc%d" % i_uc] = self.uc_vals[i_uc]


                master_data = pandas.DataFrame(master_data)
                master_data["gain"] = 1 #self.Xall[self.gain_xpos]
                self.print(master_data.to_string(float_format="%2.7g"))

    def print_step_grads(self):
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for i, (n, v) in enumerate(zip(names, vals)):
            grad = self.grad[self.ucell_xstart[self._i_shot] + i]
            ucell_labels.append('G%s=%+2.7g' % (n, grad))

        if self.big_dump and HAS_PANDAS:
            master_data ={"GNcells": self.Gncells_vals,
                          "Gscale": self.Gscale_vals,
                          "Grotx": self.Grotx,
                          "Groty": self.Groty,
                          "Grotz": self.Grotz, "GorigZ": self.GorigZ_vals}
            for i_uc in range(self.n_ucell_param):
                master_data["Guc%d" %i_uc]= self.Guc_vals[i_uc]
            master_data = pandas.DataFrame(master_data)
            master_data["Ggain"] = 0 #self.grad[self.gain_xpos]
            self.print(master_data.to_string(float_format="%2.7g"))

        if self.calc_curvatures:
            if self.big_dump and HAS_PANDAS:
                if self.refine_Umatrix or self.refine_Bmatrix or self.refine_crystal_scale or self.refine_ncells:
                    master_data = {"CUNcells": self.CUncells_vals,
                                   "CUscale": self.CUscale_vals,
                                   "CUrotx": self.CUrotx,
                                   "CUroty": self.CUroty,
                                   "CUrotz": self.CUrotz, "CUorigZ": self.CUorigZ_vals}

                    for i_uc in range(self.n_ucell_param):
                        master_data["CUuc%d" % i_uc] = self.CUuc_vals[i_uc]

                    master_data = pandas.DataFrame(master_data)
                    master_data["CUgain"] = 0 #self.curv[self.gain_xpos]
                    self.print(master_data.to_string(float_format="%2.7g"))

        # Compute the mean, min, max, variance  and median crystal scale
        # Note we must also include the spot_scale in the diffBragg instance if its not unity
        _sv = [self.D.spot_scale*s for s in self.scale_vals]
        stats = (median(_sv),
                 mean(_sv),
                 min(_sv),
                 max(_sv),
                 std(_sv))
        scale_stat_names =["median", "mean", "min", "max", "sigma"]
        scale_stats = ["%s=%.4f" % name_stat for name_stat in zip(scale_stat_names, stats)]
        exper_name = self.FNAMES[self._i_shot] if self.FNAMES is not None else "unspecified"
        scale_stats_string = "Exp %s: SCALE FACTOR STATS: " % exper_name + ", ".join(scale_stats)
        #if len(_sv) == 1:
        scale_stats_string = "Exp %s: spot scale=%1.3g, " % (exper_name, _sv[0])
        if self.scale_vals_truths is not None:
            scale_resid = [ABS(s-stru) for s, stru in zip(_sv, self.scale_vals_truths)]
            scale_stats_string += ", truth_resid=%1.2e" % median(scale_resid)
        #self.parameters.add_spot_scale(exper_name, _sv[0])

        #if self.curv is not None:
        #    if self.curv[self.spot_scale_xpos[self._i_shot]] <=0:
        #        print("Neg spot scale curve!!")
        #    for i in range(self.n_ncells_param):
        #        if self.curv[self.ncells_xstart[self._i_shot]+ i] <= 0:
        #            print("Neg ncells abc curv %d" % i)
        #    for i in range(3):
        #        if self.curv[self.ncells_def_xstart[self._i_shot]+ i] <= 0:
        #            print("Neg_def ncells abc curv %d" % i)
        #    for i in range(self.n_ucell_param):
        #        if self.curv[self.ucell_xstart[self._i_shot] + i] <= 0:
        #            print("Neg ucell curv %d" % i)
        #    if self.curv[self.rotX_xpos[self._i_shot]] <= 0:
        #        print("neg rot X curv!")
        #    if self.curv[self.rotY_xpos[self._i_shot]] <= 0:
        #        print("neg rot Y curv!")
        #    if self.curv[self.rotZ_xpos[self._i_shot]] <= 0:
        #        print("neg rot Z curv!")

        ncells_shot = self._get_m_val(self._i_shot)
        if len(ncells_shot) == 1:
            ncells_shot = [ncells_shot[0]]*3
        scale_stats_string += "Ncells=%.1f, %.1f, %.1f, " % tuple(ncells_shot)
        #self.parameters.add_Ncells_abc(exper_name, ncells_shot)

        ncells_def_vals = self._get_ncells_def_vals(self._i_shot)
        scale_stats_string += "Ncells_def=%.1f, %.1f, %.1f, " % tuple(ncells_def_vals)
        #self.parameters.add_Ncells_def(exper_name, ncells_def_vals)

        degree = "\u00b0"
        angstrom = "\u212b"
        eta_vals = self._get_eta(i_shot=self._i_shot)
        eta_a, eta_b, eta_c = eta_vals
        scale_stats_string += u"Eta=%.3f%s, %.3f%s, %.3f%s, " % (eta_a, degree, eta_b, degree, eta_c, degree)
        #self.parameters.add_eta(exper_name, eta_vals)
        lam01 = self._get_spectra_coefficients()
        if not lam01:
            lam0, lam1 = 0,1
        else:
            lam0, lam1 = lam01
        scale_stats_string += u"wavelen_offset=%1.2g %s, " % (lam0, angstrom)
        scale_stats_string += u"wavelen_scale=%.3f, " % (lam1)
        #self.parameters.add_wavelen_offset(exper_name, lam0)
        #self.parameters.add_wavelen_scale(exper_name, lam1)

        uc_string = ""
        ucparam_names = self.UCELL_MAN[0].variable_names
        for i_ucparam, ucparam_lst in enumerate(self.uc_vals):
            param_val = median(ucparam_lst)
            uc_string += "%s=%.2f, " % (ucparam_names[i_ucparam], param_val)
        scale_stats_string += uc_string

        Bmat = self.get_refined_Bmatrix(self._i_shot)
        #self.parameters.add_Bmatrix(exper_name, Bmat)

        rotX = self._get_rotX(self._i_shot)
        rotY = self._get_rotY(self._i_shot)
        rotZ = self._get_rotZ(self._i_shot)
        #self.parameters.add_rotXYZ(exper_name, (rotX, rotY, rotZ))
        rot_labels = ["rotX=%1.1e" % rotX, "rotY=%1.1e" % rotY, "rotZ=%1.1e" % rotZ]
        for r in rot_labels:
            scale_stats_string += u"%s%s, " % (r, degree)
        scale_stats_string += "detZ_offset=%1.2g m, " % median(self.origZ_vals)

        Xnorm = norm(self.x)  # NOTEX
        R1 = -1
        R1_i = -1
        self.R_overall = -1
        ncurv = 0
        if self.calc_curvatures:
            ncurv = len(self.curv > 0)

        if self.refine_per_spot_scale:
            for i_shot in self.per_spot_scale_xpos:
                nspots = len(self.NANOBRAGG_ROIS[i_shot])
                if self.selection_flags is not None:
                    vals = [self._get_per_spot_scale(i_shot, i_spot) for i_spot in range(nspots) if self.selection_flags[i_shot][i_spot]]
                else:
                    vals = [self._get_per_spot_scale(i_shot, i_spot) for i_spot in range(nspots)]
                m = median(vals)
                mx = max(vals)
                mn = min(vals)
                s = std(vals)
                self.print("Per spot scales shot %d: \n\tmin=%10.7f, \n\tmax=%10.7f, \n\tmedian=%10.7f, \n\tstdev=%10.7f" % (i_shot, mn, mx, m, s))

        if self.Fref is not None and self.iterations % self.merge_stat_frequency == 0 and self.I_AM_ROOT:
            self.R_overall = self.Fobs_Fref_Rfactor(use_binning=False, auto_scale=self.scale_r1)
            self.CC_overall = self.Fobs.correlation(self.Fref_aligned).coefficient()
            self.print("R-factor overall: %.4f, CC overall: %.4f" % (self.R_overall, self.CC_overall))
            #if self.print_resolution_bins:
            #    self.print("R-factor (shells):")
            #    self.print(self.Fobs_Fref_Rfactor(use_binning=True, auto_scale=self.scale_r1).show())
            #    self.print("CC (shells):")
            #    self.Fobs.correlation(self.Fref_aligned, use_binning=True).show()

        self.print(
            "%s%sF=%2.7g, |G|=%2.7g, eps*|X|=%2.7g,%s R1=%2.7g (R1 at start=%2.7g), Fcell negative_models=%d, Neg. Curv.: %d/%d on shots=%s\n"
            % (scale_stats_string, Bcolors.OKBLUE, self._f, self.gnorm, Xnorm * self.trad_conv_eps, Bcolors.ENDC, self.R_overall, self.init_R1,
               self.tot_fcell_negative_model, self.tot_neg_curv, ncurv,
               ", ".join(map(str, self.neg_curv_shots))))
        #self.print("<><><><><><><><> TOP GUN <><><><><><><><>")
        #self.print("                 End of iteration.")
        self._print_sausages()

        if self.testing_mode:
            self.conv_test()

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

    def _print_sausages(self):
        sausage_vals = self._get_sausage_parameters(self._i_shot)
        rotX = map(lambda x: "%1.2e"%x, sausage_vals[0::4])
        rotY = map(lambda x: "%1.2e"%x, sausage_vals[1::4])
        rotZ = map(lambda x: "%1.2e"%x, sausage_vals[2::4])
        scales = map(lambda x: "%1.5e"%x, sausage_vals[3::4])

        self.print("Blue Sausage rotX (rad.): %s" % (", ".join(rotX)))
        self.print("Blue Sausage rotY (rad.): %s" % (", ".join(rotY)))
        self.print("Blue Sausage rotZ (rad.): %s" % (", ".join(rotZ)))
        self.print("Blue Sausage Scale: %s" % (", ".join(scales)))

    def get_refined_Bmatrix(self, i_shot):
        return self.UCELL_MAN[i_shot].B_recipspace

    def curvatures(self):
        return self.curv

    def get_correction_misset(self, as_axis_angle_deg=False, anglesXYZ=None, i_shot=None):
        """
        return the current state of the perturbation matrix
        :return: scitbx.matrix sqr
        """
        if anglesXYZ is None:
            assert i_shot is not None
            rx = self._get_rotX(i_shot)
            ry = self._get_rotY(i_shot)
            rz = self._get_rotZ(i_shot)
            anglesXYZ = rx, ry, rz

        x = col((-1, 0, 0))
        y = col((0, -1, 0))
        z = col((0, 0, -1))
        RX = x.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[0], deg=False)
        RY = y.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[1], deg=False)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[2], deg=False)
        M = RX * RY * RZ
        if as_axis_angle_deg:
            q = M.r3_rotation_matrix_as_unit_quaternion()
            rot_ang, rot_ax = q.unit_quaternion_as_axis_and_angle(deg=True)
            return rot_ang, rot_ax
        else:
            return M

    def Fobs_Fref_Rfactor(self, use_binning=False, auto_scale=False):
        if auto_scale:
            # TODO check for convergence of minimizer and warn if it faile, then set scale factor to 1 ?
            self.r1_scale = minimize(LocalRefiner._rfactor_minimizer_target,
                                     x0=[1], args=(self.Fobs, self.Fref_aligned),
                                     method='Nelder-Mead').x[0]
        else:
            self.r1_scale = 1

        return self.Fobs.r1_factor(self.Fref_aligned,
                                   use_binning=use_binning, scale_factor=self.r1_scale)

    @staticmethod
    def _rfactor_minimizer_target(k, Fobs, Fref):
        return Fobs.r1_factor(Fref, scale_factor=k[0])

    def get_optimized_mtz(self, save_to_file=None, wavelength=1):
        assert self.symbol is not None
        from cctbx import crystal
        # TODO update for non global unit cell case (average over unit cells)

        self._update_Fcell()  # just in case update the Fobs

        um = self.UCELL_MAN[0]
        sym = crystal.symmetry(unit_cell=um.unit_cell_parameters, space_group_symbol=self.symbol)
        mset_obs = miller.set(sym, self.Fobs.indices(), anomalous_flag=True)
        fobs = miller.array(mset_obs, self.Fobs.data()).set_observation_type_xray_amplitude()
        # TODO: what to do in MPI mode when writing ?
        if save_to_file is not None and self.I_AM_ROOT:
            fobs.as_mtz_dataset(column_root_label='fobs', wavelength=wavelength).mtz_object().write(save_to_file)
        return fobs

    def get_corrected_crystal(self, i_shot=0):
        ang, ax = self.get_correction_misset(as_axis_angle_deg=True, i_shot=i_shot)
        B = self.get_refined_Bmatrix(i_shot)
        C = deepcopy(self.CRYSTAL_MODELS[i_shot])
        C.set_B(B)
        if ang > 0:
            C.rotate_around_origin(ax, ang)
        return C

    def conv_test(self):
        assert self.symbol is not None
        err = []
        s = ""
        A = []
        vars = self._get_ucell_vars(0)
        for i in range(self.n_ucell_param):
            if self.rescale_params:
                a = vars[i]
            else:
                a = self.Xall[self.ucell_xstart[0] + i]

            if i == 3:
                a = a * 180 /PI

            s += "%.4f " % a
            A.append(a)
            err.append(ABS(self.gt_ucell[i] - a))

        mn_err = sum(err) / self.n_ucell_param

        shot_refined = []
        all_det_resid = []
        for i_shot in range(self.n_shots):
            Ctru = self.CRYSTAL_GT[i_shot]
            atru, btru, ctru = Ctru.get_real_space_vectors()
            ang, ax = self.get_correction_misset(as_axis_angle_deg=True, i_shot=i_shot)
            B = self.get_refined_Bmatrix(i_shot=i_shot)
            C = deepcopy(self.CRYSTAL_MODELS[i_shot])
            C.set_B(B)
            if ang > 0:
                C.rotate_around_origin(ax, ang)
            try:
                ang_off = compare_with_ground_truth(atru, btru, ctru,
                                                    [C],
                                                    symbol=self.symbol)[0]
            except RuntimeError:
                ang_off = 999

            out_str = "shot %d: MEAN UCELL ERROR=%.4f, ANG OFF %.4f" % (i_shot, mn_err, ang_off)
            ncells_val = self._get_m_val(i_shot)[0]
            ncells_resid = abs(ncells_val - self.gt_ncells)

            if mn_err < 0.01 and ang_off < 0.004 and ncells_resid < 0.1:
                shot_refined.append(True)
            else:
                shot_refined.append(False)

            if self.refine_detdist:
                det_resid = abs(self.detector_distance_gt[i_shot] - self._get_detector_distance_val(i_shot))
                all_det_resid.append(det_resid)
                out_str += ", OrigZ resid = %.4f" % det_resid

            if self.refine_ncells:
                out_str += ", ncells resid=%.4f" % ncells_resid
            self.print(out_str)

        if all(shot_refined):
            if self.refine_detdist:
                if all([det_resid < 0.01 for det_resid in all_det_resid]):
                    self.print("OK")
                    exit()
            else:
                self.print("OK")
                exit()

    # NOTE below are functions which need to be overwritten in the Global MPI class

    def _MPI_sync_hkl_freq(self):
        pass

    def _MPI_sync_fcell_parameters(self):
        pass

    def _MPI_sync_panel_params(self):
        pass

    def _data_for_write(self, parameter_dict):
        return [parameter_dict]

    def _MPI_aggregate_model_data_correlations(self):
        self.all_image_corr = [self._get_image_correlation(i) for i in self.shot_ids]
        if self.init_image_corr is None:
            self.init_image_corr = deepcopy(self.all_image_corr)

    def _init_n_bad_shots(self):
        self.n_bad_shots = len(self.bad_shot_list)
        return self.n_bad_shots

    def _init_gather_ang_off(self):
        return self.all_ang_off

    def _get_ang_off(self):
        return self.all_ang_off

    def _MPI_reduce_broadcast(self, var):
        return var

    def _MPI_combine_data_to_send(self):
        return self.data_to_send

    def _combine_data_to_save(self):
        # Here we can save the refined parameters
        my_shots = self.NANOBRAGG_ROIS.keys()
        x = self.Xall
        self.data_to_send = []
        image_corr = self.image_corr
        if image_corr is None:
            image_corr = [-1] * len(my_shots)

        for i_shot in my_shots:
            rotX = self._get_rotX(i_shot)
            rotY = self._get_rotY(i_shot)
            rotZ = self._get_rotZ(i_shot)

            ang, ax = self.get_correction_misset(as_axis_angle_deg=True, anglesXYZ=(rotX, rotY, rotZ))
            pars = self._get_ucell_vars(i_shot)
            self.UCELL_MAN[i_shot].variables = pars
            Bmat = self.UCELL_MAN[i_shot].B_recipspace

            bg_coef = -1
            C = deepcopy(self.CRYSTAL_MODELS[i_shot])
            a_init, b_init, c_init, al_init, be_init, ga_init = C.get_unit_cell().parameters()
            C.set_B(Bmat)
            try:
                C.rotate_around_origin(ax, ang)
            except RuntimeError:
                pass

            Amat_refined = C.get_A()
            a,b,c,al,be,ga = C.get_unit_cell().parameters()

            panX, panY, panZ = zip(*[self._get_panelXYZ_val(pid) for pid in range(len(self.S.detector))])
            if self.refine_detdist:
                panZ = [self._get_detector_distance_val(i_shot)] * len(self.S.detector)
            panO, panF, panS = zip(*[self._get_panelRot_val(pid) for pid in range(len(self.S.detector))])

            coeffs = self._get_spectra_coefficients()
            if not coeffs:
                lam0, lam1 = -1,-1
            else:
                lam0, lam1 = coeffs

            fcell_xstart = self.fcell_xstart
            ucell_xstart = self.ucell_xstart[i_shot]
            scale_xpos = self.spot_scale_xpos[i_shot]
            ncells_xstart = self.ncells_xstart[i_shot]
            nspots = len(self.NANOBRAGG_ROIS[i_shot])

            bgplane_a_xpos = [self.bg_a_xstart[i_shot][i_spot] for i_spot in range(nspots)]
            bgplane_b_xpos = [self.bg_b_xstart[i_shot][i_spot] for i_spot in range(nspots)]
            bgplane_c_xpos = [self.bg_c_xstart[i_shot][i_spot] for i_spot in range(nspots)]
            bgplane_xpos = list(zip(bgplane_a_xpos, bgplane_b_xpos, bgplane_c_xpos))
            bgplane = [self._get_bg_vals(i_shot, i_spot) for i_spot in range(nspots)]

            crystal_scale = self._get_spot_scale(i_shot)**2 * self.D.spot_scale
            eta = tuple(self._get_eta(i_shot))
            proc_h5_fname = "" #self.all_proc_fnames[i_shot]
            proc_h5_idx = "" #self.all_shot_idx[i_shot]
            proc_bbox_idx = "" #self.all_proc_idx[i_shot]
            if self.FNAMES is not None:
                exp_name = self.FNAMES[i_shot]
            else:
                exp_name = ""

            ncells_val = tuple(self._get_m_val(i_shot))
            ncells_def_vals = tuple(self._get_ncells_def_vals(i_shot))
            init_misori = self.init_ang_off[i_shot] if self.init_ang_off is not None else -1
            img_corr = -1  #image_corr[i_shot]
            init_img_corr = -1
            final_misori = self.get_current_misorientation(i_shot)

            self.data_to_send.append(
                (proc_h5_fname, proc_h5_idx, proc_bbox_idx, exp_name, crystal_scale, eta, Amat_refined, ncells_val,
                 ncells_def_vals, bgplane,
                 img_corr, init_img_corr, fcell_xstart, ucell_xstart, rotX, rotY, rotZ, scale_xpos,
                 ncells_xstart, bgplane_xpos, init_misori, final_misori, a,b,c,al,be,ga,
                 a_init, b_init, c_init, al_init, be_init, ga_init, bg_coef, lam0, lam1, panX, panY, panZ,
                 panO, panF, panS))

    def get_lbfgs_x_array_as_dataframe(self):
        self._combine_data_to_save()

        data = self._MPI_combine_data_to_send()

        if self.I_AM_ROOT:
            import pandas

            fnames, shot_idx, bbox_idx, exp_name, xtal_scales, eta, Amats, ncells_vals, ncells_def_vals,bgplanes, image_corr, init_img_corr, \
                fcell_xstart, ucell_xstart, rotX, rotY, rotZ, scale_xpos, ncells_xstart, bgplane_xpos, \
                init_misori, final_misori, a, b, c, al, be, ga, \
                a_init, b_init, c_init, al_init, be_init, ga_init, bg_coef, lam0, lam1, panX, panY, panZ, \
                panO, panF, panS = zip(*data)

            df = pandas.DataFrame({"proc_fnames": fnames, "proc_shot_idx": shot_idx, "bbox_idx": bbox_idx,
                                   "exp_name": exp_name,
                                   "panX": list(panX), "panY": list(panY), "panZ": list(panZ),
                                   "panO": list(panO), "panF": list(panF), "panS": list(panS),
                                   "spot_scales": xtal_scales, "Amats": Amats, "ncells": ncells_vals,
                                   "ncells_def": ncells_def_vals,
                                   "bgplanes": bgplanes, "image_corr": image_corr,
                                   "init_image_corr": init_img_corr,
                                   "fcell_xstart": fcell_xstart,
                                   "ucell_xstart": ucell_xstart,
                                   "init_misorient": init_misori, "final_misorient": final_misori,
                                   "bg_coef": bg_coef,
                                   "eta": eta,
                                   "rotX": rotX,
                                   "rotY": rotY,
                                   "rotZ": rotZ,
                                   "a": a, "b": b, "c": c, "al": al, "be": be, "ga": ga,
                                   "a_init": a_init, "b_init": b_init, "c_init": c_init, "al_init": al_init,
                                   "lam0": lam0, "lam1": lam1,
                                   "be_init": be_init, "ga_init": ga_init,
                                   "scale_xpos": scale_xpos,
                                   "ncells_xpos": ncells_xstart,
                                   "bgplanes_xpos": bgplane_xpos})
            return df

    def _record_xy_calc(self):
        if self.record_model_predictions:
            I = self.model_bragg_spots.ravel()
            Isum = I.sum()
            x1, x2, y1, y2 = self.ROIS[self._i_shot][self._i_spot]
            roi_shape = y2-y1, x2-x1
            Y, X = np_indices(roi_shape)
            Y = Y.ravel() + y1
            X = X.ravel() + x1
            # this is the calculated (predicted) centroid for this spot
            x_calc_pix = (X*I).sum() / Isum + 0.5
            y_calc_pix = (Y*I).sum() / Isum + 0.5
            if self._i_shot not in self.xy_calc:
                self.xy_calc[self._i_shot] = [(0,0)]*len(self.ROIS[self._i_shot])
            self.xy_calc[self._i_shot][self._i_spot] = x_calc_pix, y_calc_pix

    def get_refined_reflections(self, shot_refls, i_shot=0):
        """
        shot_refls, reflection table corresponding to the input reflections
        i_shot, int specifying the index of the shot corresponding to shot_refls
            usually this is zero in per-shot refinement mode, but in multi shot
            refinement it might be non-zero
        """
        num_spots_in_shot = len(self.ROIS[i_shot])
        if i_shot not in self.xy_calc:
            raise KeyError("No xycalc data for requested shot")

        if len(self.xy_calc[i_shot]) != num_spots_in_shot:
            raise ValueError("shot_refls should be same length as number of spots in refiner")

        selection_flags = flex.bool(self.selection_flags[i_shot])
        refined_refls = deepcopy(shot_refls)
        x, y = zip(*self.xy_calc[i_shot])
        z = [0]*num_spots_in_shot
        xyz_calc = flex.vec3_double(list(zip(x, y, z)))
        refined_refls['dials.xyzcal.px'] = deepcopy(refined_refls['xyzcal.px'])  # make a backup
        refined_refls["xyzcal.px"] = xyz_calc

        refined_refls = refined_refls.select(selection_flags)
        return refined_refls

    def _update_full_image_of_model(self):
        if self.full_image_of_model is None:
            return
        x1, x2, y1, y2 = self.ROIS[self._i_shot][self._i_spot]
        roi_shape = y2-y1, x2-x1
        self.full_image_of_model[self._panel_id, y1:y2, x1:x2] = self.model_Lambda.reshape(roi_shape)
        self.full_image_of_spots[self._panel_id, y1:y2, x1:x2] = self.model_bragg_spots.reshape(roi_shape)
        if isinstance(self.sigma_r, float):
            SIG = self.sigma_r
        else:
            try:
                SIG = self.sigma_r.reshape(roi_shape)
            except Exception as err:
                self.print("something off with sigma_r", err)
                SIG = 0
        self.full_image_of_sigma_r[self._panel_id, y1:y2, x1:x2] = SIG

    def get_model_image(self, i_shot=0, only_save_model=False):
        """
        go through each ROI and execute diffBragg, storing the model pixels along the way
        """
        self.save_model_for_shot = i_shot
        self.only_save_model_for_shot = only_save_model
        self.compute_functional_and_gradients()   # goes thorugh all ROIs
        self.save_model_for_shot = None
        self.only_save_model_for_shot = False
        return self.full_image_of_model, self.full_image_of_spots, self.full_image_of_sigma_r

    def get_optimized_detector(self, i_shot=0):
        n_panels = len(self.S.detector)
        new_det = Detector()
        for pid in range(n_panels):
            new_offsetX, new_offsetY, new_offsetZ = self._get_panelXYZ_val(pid)

            if self.refine_detdist:
                new_offsetZ = self._get_detector_distance_val(i_shot)

            panel_rot_angO, panel_rot_angF, panel_rot_angS = self._get_panelRot_val(pid)

            if self.panel_reference_from_id is not None:
                self.D.reference_origin = self.panel_reference_from_id[pid]
            else:
                self.D.reference_origin = self.S.detector[pid].get_origin()

            self.D.update_dxtbx_geoms(self.S.detector, self.S.beam.nanoBragg_constructor_beam, pid,
                                panel_rot_angO, panel_rot_angF, panel_rot_angS,
                                new_offsetX, new_offsetY, new_offsetZ, force=False)

            fdet = self.D.fdet_vector
            sdet = self.D.sdet_vector
            origin = self.D.get_origin()
            panel_dict = self.S.detector[pid].to_dict()
            panel_dict["fast_axis"] = fdet
            panel_dict["slow_axis"] = sdet
            panel_dict["origin"] = origin
            new_det.add_panel(Panel.from_dict(panel_dict))

        return new_det

    def _update_sigma_r_if_possible(self):
        if self.ROI_DARK_RMS is not None:
            darkRMS = self.ROI_DARK_RMS[self._i_shot][self._i_spot]
            assert darkRMS.min() >= 0
            return darkRMS.ravel()
        else:
            return self.sigma_r

    def _MPI_barrier(self):
        pass

    def _MPI_check_for_break_signal(self):
        pass
