from __future__ import print_function
import time
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    has_mpi = True
except ImportError:
    rank = 0
    size = 1
    has_mpi = False

class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if rank == 0:
    import os
    import pandas
    from numpy import mean, median, unique, std
    from tabulate import tabulate
    from scipy.stats import linregress
    from scipy.stats import pearsonr
    from numpy import log as np_log
    from numpy import exp as np_exp
    from numpy import load as np_load
    from numpy import all as np_all
    from numpy import abs as ABS
    from numpy import std as STD
    from numpy import save as SAVE
    from numpy import savez as SAVEZ
    from numpy.linalg import norm
    from numpy import ones_like as ONES_LIKE
    from numpy import array as ARRAY
    from numpy import pi as PI
    from numpy import allclose as ALL_CLOSE
    
    from simtbx.diffBragg.refiners import BreakToUseCurvatures
    from scitbx.array_family import flex
    from cctbx.array_family import flex as cctbx_flex
    from numpy import nan as NAN
    flex_miller_index = cctbx_flex.miller_index

    flex_double = flex.double
    from scitbx.matrix import col
    from simtbx.diffBragg.refiners import PixelRefinement
    from scipy.optimize import minimize

    from collections import Counter
    import numpy as np

else:
    mean = unique = np_log = np_exp = np_all = norm = median = std = None
    ALL_CLOSE=NAN = ARRAY = SAVE = SAVEZ = ONES_LIKE = STD = ABS = PI= None
    tabulate = None
    flex_miller_index = None
    minimize = None
    np_load = None
    linregress = None
    pearsonr = None
    ARRAY = None
    Counter = None
    # stdout_flush = None
    BreakToUseCurvatures = None
    flex_double = None
    col = None
    PixelRefinement = None

if has_mpi:
    ALL_CLOSE= comm.bcast(ALL_CLOSE)
    ARRAY =comm.bcast(ARRAY) 
    PI =comm.bcast(PI) 
    SAVE = comm.bcast(SAVE) 
    SAVEZ = comm.bcast(SAVEZ) 
    ONES_LIKE = comm.bcast(ONES_LIKE) 
    STD = comm.bcast(STD) 
    ABS = comm.bcast(ABS) 
    mean = comm.bcast(mean, root=0)
    Counter = comm.bcast(Counter, root=0)
    linregress = comm.bcast(linregress, root=0)
    pearsonr = comm.bcast(pearsonr, root=0)
    tabulate = comm.bcast(tabulate, root=0)
    std = comm.bcast(std, root=0)
    flex_miller_index = comm.bcast(flex_miller_index, root=0)
    median = comm.bcast(median, root=0)
    unique = comm.bcast(unique, root=0)
    np_load = comm.bcast(np_load, root=0)
    np_log = comm.bcast(np_log, root=0)
    np_exp = comm.bcast(np_exp, root=0)
    minimize = comm.bcast(minimize, root=0)
    norm = comm.bcast(norm, root=0)
    ARRAY = comm.bcast(ARRAY, root=0)
    np_all = comm.bcast(np_all, root=0)
    # stdout_flush = comm.bcast(stdout_flush)
    BreakToUseCurvatures = comm.bcast(BreakToUseCurvatures)
    flex_double = comm.bcast(flex_double)
    col = comm.bcast(col)
    PixelRefinement = comm.bcast(PixelRefinement)
    NAN = comm.bcast(NAN, root=0)

if rank == 0:
    import pylab as plt

# TODO move me to broadcasts
import sys
import warnings
from copy import deepcopy
from simtbx.diffBragg.utils import compare_with_ground_truth
from cctbx import miller, sgtbx
import itertools
#import numexpr
from IPython import embed

warnings.filterwarnings("ignore")

class FatRefiner(PixelRefinement):

    def __init__(self, n_total_params, n_local_params, n_global_params, local_idx_start,
                 shot_ucell_managers, shot_rois, shot_nanoBragg_rois,
                 shot_roi_imgs, shot_spectra, shot_crystal_GTs,
                 shot_crystal_models, shot_xrel, shot_yrel, shot_abc_inits, shot_asu,
                 global_param_idx_start,
                 shot_panel_ids,
                 log_of_init_crystal_scales=None,
                 all_crystal_scales=None, init_gain=1, perturb_fcell=False,
                 global_ncells=False, global_ucell=True, global_originZ=True,
                 shot_originZ_init=None,
                 sgsymbol="P43212",
                 omega_kahn=None, selection_flags=None):
        """
        TODO: parameter x array boundaries should be done in this class, eliminating the need for local_idx_start
        TODO and global_idx_start
        TODO: ROI and nanoBragg ROI should be single variable

        :param n_total_params:  total number of parameters all shots
        :param n_local_params:  total number of parameters for each shot  (e.g. xtal rotation matrices)
        :param n_global_params: total number of parameters used by all shots (e.g. Fhkl)
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
        :param shot_originZ_init: per shot origin Z
        :param global_originZ: do we refine a single detector Z position for all shots (default is True)
        :param omega_kahn: omega and kahn correction term for each panel
        """
        PixelRefinement.__init__(self)
        # super(GlobalFat, self).__init__()
        self.num_kludge = 0
        self.global_ncells_param = global_ncells
        self.global_ucell_param = global_ucell
        self.global_originZ_param = global_originZ
        self.debug = False
        self.num_Fcell_kludge = 0
        self.perturb_fcell = perturb_fcell
        # dictionaries whose keys are the shot indices
        self.UCELL_MAN = shot_ucell_managers
        self.CRYSTAL_SCALE_TRUTH = all_crystal_scales  # ground truth of crystal scale factors.. 
        # cache shot ids and make sure they are identical in all other input dicts
        self.shot_ids = sorted(shot_ucell_managers.keys())
        self.big_dump = False
        self.show_watched = False
        self.image_corr = None
        self.n_shots = len(self.shot_ids)
        self.shot_originZ_init = None
        if shot_originZ_init is not None:
            self.shot_originZ_init = self._check_keys(shot_originZ_init)
        
        self.selection_flags = selection_flags
        if self.selection_flags is not None: 
            self.selection_flags = self._check_keys(self.selection_flags)
        
        # sanity check: no repeats of the same shot
        assert len(self.shot_ids) == len(set(self.shot_ids))

        self.OMEGA_KAHN = omega_kahn
        self.ROIS = self._check_keys(shot_rois)
        self.ASU = self._check_keys(shot_asu)
        self.NANOBRAGG_ROIS = self._check_keys(shot_nanoBragg_rois)
        self.ROI_IMGS = self._check_keys(shot_roi_imgs)
        self.SPECTRA = self._check_keys(shot_spectra)
        self.CRYSTAL_GT = None
        if shot_crystal_GTs is not None:
            self.CRYSTAL_GT = self._check_keys(shot_crystal_GTs)
        self.CRYSTAL_MODELS = self._check_keys(shot_crystal_models)
        self.XREL = self._check_keys(shot_xrel)
        self.YREL = self._check_keys(shot_yrel)
        self.ABC_INIT = self._check_keys(shot_abc_inits)
        self.PANEL_IDS = self._check_keys(shot_panel_ids)

        # Total number of parameters in the MPI world
        self.n_total_params = n_total_params

        # total number of local parameters
        self.n_local_params = n_local_params

        # total number of params varying globally across all shots
        self.n_global_params = n_global_params

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

        self.n_per_shot_originZ_param = 1
        if global_originZ:
            self.n_per_shot_originZ_param = 0
        else:
            assert shot_originZ_init is not None

        self.n_per_shot_ucell_param = self.n_ucell_param
        if global_ucell:
            self.n_per_shot_ucell_param = 0

        self.n_per_shot_ncells_param = 1
        if global_ncells:
            self.n_per_shot_ncells_param = 0

        self.n_per_shot_params = self.n_rot_param + self.n_spot_scale_param \
            + self.n_per_shot_ncells_param + self.n_per_shot_ucell_param \
            + self.n_per_shot_originZ_param

        self._ncells_id = 9  # diffBragg internal index for Ncells derivative manager
        self._originZ_id = 10  # diffBragg internal index for originZ derivative manager
        self._fcell_id = 11  # diffBragg internal index for Fcell derivative manager

        if log_of_init_crystal_scales is None:
            log_of_init_crystal_scales = {s: 0 for s in self.shot_ids}
        else:
            assert sorted(log_of_init_crystal_scales.keys()) == self.shot_ids
        self.log_of_init_crystal_scales = log_of_init_crystal_scales

        self._init_gain = init_gain
        self.num_positive_curvatures = 0
        self._panel_id = None
        self.symbol = sgsymbol
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
        self.FNAMES = None
        self.PROC_FNAMES = None

    def setup_plots(self):
        if rank == 0:
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
        return len(self.x)
        #return self.n_total_params

    @property
    def n_global_fcell(self):
        return len(self.idx_from_asu)

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
            raise KeyError("input data funky, check GlobalFat inputs")
        return shot_dict

    def _evaluate_averageI(self):
        """model_Lambda means expected intensity in the pixel"""
        self.model_Lambda = \
            self.gain_fac * self.gain_fac * (self.tilt_plane + self.model_bragg_spots)

    def _setup(self):
        # Here we go!  https://youtu.be/7VvkXA6xpqI
        if rank == 0:
            print("Setup begins!")
        if not self.asu_from_idx:
            raise ValueError("Need to supply a non empty asu from idx map")
        if not self.idx_from_asu:  # # TODO just derive from its inverse
            raise ValueError("Need to supply a non empty idx from asu map")

        # get the Fhkl information from P1 array internal to nanoBragg
        if rank == 0:
            print("--0 create an Fcell mapping")
        idx, data = self.S.D.Fhkl_tuple
        self.idx_from_p1 = {h: i for i, h in enumerate(idx)}
        # self.p1_from_idx = {i: h for i, h in zip(idx, data)}

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        # Make the global sized parameter array, though here we only update the local portion
        self.x = flex_double(self.n_total_params)

        # store the starting positions in the parameter array for this shot
        self.rotX_xpos = {}
        self.rotY_xpos = {}
        self.rotZ_xpos = {}
        self.ucell_xstart = {}
        self.ncells_xpos = {}
        self.originZ_xpos = {}
        self.spot_scale_xpos = {}
        self.n_panels = {}
        self.bg_a_xstart = {}
        self.bg_b_xstart = {}
        self.bg_c_xstart = {}
        if rank == 0:
            print("--1 Setting up per shot parameters")

        # this is a sliding parameter that points to the latest local (per-shot) parameter in the x-array
        _local_pos = self.local_idx_start
        for i_shot in self.shot_ids:
            self.pid_from_idx[i_shot] = {i: pid for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.idx_from_pid[i_shot] = {pid: i for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.n_panels[i_shot] = len(self.pid_from_idx[i_shot])

            n_spots = len(self.NANOBRAGG_ROIS[i_shot])
            self.bg_a_xstart[i_shot] = []
            self.bg_b_xstart[i_shot] = []
            self.bg_c_xstart[i_shot] = []
            _spot_start = _local_pos
            for i_spot in range(n_spots):
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
                    self.x[self.bg_a_xstart[i_shot][i_spot]] = 1
                    self.x[self.bg_b_xstart[i_shot][i_spot]] = 1
                    self.x[self.bg_c_xstart[i_shot][i_spot]] = 1

                else:
                    self.x[self.bg_a_xstart[i_shot][i_spot]] = float(a)
                    self.x[self.bg_b_xstart[i_shot][i_spot]] = float(b)
                    self.x[self.bg_c_xstart[i_shot][i_spot]] = float(c)

                _spot_start += 3

            self.rotX_xpos[i_shot] = self.bg_c_xstart[i_shot][-1] + 1
            self.rotY_xpos[i_shot] = self.rotX_xpos[i_shot] + 1
            self.rotZ_xpos[i_shot] = self.rotY_xpos[i_shot] + 1

            if self.rescale_params:
                self.x[self.rotX_xpos[i_shot]] = 1
                self.x[self.rotY_xpos[i_shot]] = 1
                self.x[self.rotZ_xpos[i_shot]] = 1
            else:
                self.x[self.rotX_xpos[i_shot]] = 0
                self.x[self.rotY_xpos[i_shot]] = 0
                self.x[self.rotZ_xpos[i_shot]] = 0

            # continue adding local per shot parameters after rotZ_xpos
            _local_pos = self.rotZ_xpos[i_shot] + 1

            # global always starts here, we have to decide whether to put ncells / unit cell/ originZ parameters in global array
            _global_pos = self.global_param_idx_start

            if self.global_ucell_param:
                self.ucell_xstart[i_shot] = _global_pos
                _global_pos += self.n_ucell_param
            else:
                self.ucell_xstart[i_shot] = _local_pos
                _local_pos += self.n_ucell_param
                for i_cell in range(self.n_ucell_param):
                    if self.rescale_params:
                        self.x[self.ucell_xstart[i_shot] + i_cell] = 1  #self.UCELL_MAN[i_shot].variables[i_cell]
                    else:
                        self.x[self.ucell_xstart[i_shot] + i_cell] = self.UCELL_MAN[i_shot].variables[i_cell]


            if self.global_ncells_param:
                self.ncells_xpos[i_shot] = _global_pos
                _global_pos += self.n_ncells_param
            else:
                self.ncells_xpos[i_shot] = _local_pos
                _local_pos += self.n_ncells_param
                ncells_xval = np_log(self.S.crystal.Ncells_abc[0]-3)
                # TODO: each shot gets own starting Ncells
                if self.rescale_params:
                    self.x[self.ncells_xpos[i_shot]] = 1  #ncells_xval
                else:
                    self.x[self.ncells_xpos[i_shot]] = ncells_xval

            if self.global_originZ_param:
                self.originZ_xpos[i_shot] = _global_pos
                _global_pos += 1
            else:
                self.originZ_xpos[i_shot] = _local_pos
                _local_pos += 1
                if self.rescale_params:
                    self.x[self.originZ_xpos[i_shot]] = 1
                else:
                    self.x[self.originZ_xpos[i_shot]] =self.S.detector[0].get_origin()[2]  #TODO fixme local_origin vs origin

            self.spot_scale_xpos[i_shot] = _local_pos
            _local_pos += 1
            if self.rescale_params:
                self.x[self.spot_scale_xpos[i_shot]] = 1
            else:
                self.x[self.spot_scale_xpos[i_shot]] = self.log_of_init_crystal_scales[i_shot]

        self.fcell_xstart = _global_pos
        self.gain_xpos = self.n_total_params - 1

        # tally up HKL multiplicity
        if rank == 0:
            print("REduction of global data layout")
        hkl_totals = []
        fname_totals = []
        panel_id_totals = []
        # img_totals = []
        for i_shot in self.ASU:
            for i_h, h in enumerate(self.ASU[i_shot]):
                if self.FNAMES is not None:
                    fname_totals.append(self.FNAMES[i_shot])
                panel_id_totals.append(self.PANEL_IDS[i_shot][i_h])
                hkl_totals.append(self.idx_from_asu[h])
                # img_totals.append(self.ROI_IMGS[i_shot][i_h])
        hkl_totals = self._mpi_reduce_broadcast(hkl_totals)

        if rank == 0:
            print("--2 Setting up global parameters")
            # put in estimates for origin vectors
            # TODO: refine at the different hierarchy
            # get te first Z coordinate for now..
            # print("Setting origin: %f " % self.S.detector[0].get_local_origin()[2])
            if self.global_ucell_param:
                # TODO have parameter for global init of unit cell, right now its handled in the global_bboxes scripts
                for i_cell in range(self.n_ucell_param):
                    if self.rescale_params:
                        self.x[self.ucell_xstart[0] + i_cell] = 1
                    else:
                        self.x[self.ucell_xstart[0] + i_cell] = self.UCELL_MAN[0].variables[i_cell]

            if self.global_ncells_param:
                # TODO have parameter for global init of Ncells , right now its handled in the global_bboxes scripts
                if self.rescale_params:
                    ncells_xval = 1
                else:
                    ncells_xval = np_log(self.S.crystal.Ncells_abc[0] - 3)
                self.x[self.ncells_xpos[0]] = ncells_xval

            if self.global_originZ_param:
                # TODO have parameter for global init of originZ param , right now its handled in the global_bboxes scripts
                if self.rescale_params:
                    self.x[self.originZ_xpos[0]] = 1  #self.S.detector[0].get_local_origin()[2]  # NOTE maybe just origin instead?elf.S.detector
                else:
                    self.x[self.originZ_xpos[0]] = self.S.detector[0].get_origin()[2]  # NOTE maybe just origin instead?elf.S.detector
                    #self.x[self.originZ_xpos[0]] = self.S.detector[0].get_local_origin()[2]  # NOTE maybe just origin instead?elf.S.detector

            print("----loading fcell data")
            # this is the number of observations of hkl (accessed like a dictionary via global_fcell_index
            print("---- -- counting hkl totes")
            self.hkl_frequency = Counter(hkl_totals)

            # initialize the Fhkl global values
            print("--- --- --- inserting the Fhkl array in the parameter array... ")
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
                    self.x[self.fcell_xstart + i_fcell] = 1
                else:
                    self.x[self.fcell_xstart + i_fcell] = vals[i_fcell]

            self.Fref_aligned = self.Fref
            if self.Fref is not None:
                self.Fref_aligned = self.Fref.select_indices(self.Fobs.indices())
                self.init_R1 = self.Fobs_Fref_Rfactor(use_binning=False, auto_scale=self.scale_r1)
                print("Initial R1 = %.4f" % self.init_R1)

            if self.Fobs is not None:  # TODO should this ever be None ?
                miller_binner = self.Fobs.binner()
                miller_bin_idx = miller_binner.bin_indices()

                import numpy as np # TODO move me to top
                from simtbx.diffBragg.utils import nearest_non_zero

                unique_bins = sorted(set(miller_bin_idx))
                sigmas = []
                for i_bin in unique_bins:
                    dmax, dmin = miller_binner.bin_d_range(i_bin)
                    f_selection = self.Fobs.resolution_filter(d_min=dmin, d_max=dmax)
                    fsel_data = f_selection.data().as_numpy_array()
                    if self.log_fcells:
                        fsel_data = np_log(fsel_data)
                    sigma = STD(f_selection.data())
                    sigmas.append(sigma) #sigma_for_res_id[i_bin] = sigma
                #min_sigma = min(self.sigma_for_res_id.values())
                #max_sigma = max(self.sigma_for_res_id.values())
                #median_sigma = np.median(self.sigma_for_res_id.values())
                self.sigma_for_res_id = {}
                for ii, sigma in enumerate(sigmas):
                    i_bin = unique_bins[ii]
                    if sigma == 0:
                        sigma = nearest_non_zero(sigmas, ii)
                    if sigma == 0:
                        bin_rng = miller_binner.bin_d_range(i_bin)
                        raise ValueError("sigma is being set to 0 for all fcell in range %.4f - %.4f" % bin_rng)
                    self.sigma_for_res_id[i_bin] = sigma

                self.res_group_id_from_fcell_index = {}
                for ii, asu_index in enumerate(miller_binner.miller_indices()):
                    if asu_index not in self.idx_from_asu:
                        raise KeyError("something wrong Fobs does not contain the asu indices")
                    i_fcell = self.idx_from_asu[asu_index]
                    self.res_group_id_from_fcell_index[i_fcell] = miller_bin_idx[ii]

            if self.output_dir is not None:
                #np.save(os.path.join(self.output_dir, "f_truth"), self.f_truth)  #FIXME by adding in the correct truth from Fref
                SAVE(os.path.join(self.output_dir, "f_asu_map"), self.asu_from_idx)

            # set gain TODO: implement gain dependent statistical model ? Per panel or per gain mode dependent ?
            self.x[self.gain_xpos] = self._init_gain  # gain factor

        # reduce then broadcast self.x
        if not rank==0:
            self.sigma_for_res_id = None 
            self.res_group_id_from_fcell_index = None
            self.fcell_init = None
        if self.rescale_params:
            self.fcell_init = comm.bcast(self.fcell_init)
            self.sigma_for_res_id = comm.bcast(self.sigma_for_res_id)
            self.res_group_id_from_fcell_index = comm.bcast(self.res_group_id_from_fcell_index)

        if rank == 0:
            print("--3 combining parameters across ranks")
        self.x = self._mpi_reduce_broadcast(self.x)

        if rank != 0:
            self.hkl_frequency = None
        if has_mpi:
            self.hkl_frequency = comm.bcast(self.hkl_frequency)

        # See if restarting from save state

        if self.x_init is not None:
            self.x = self.x_init
        elif self.restart_file is not None:
            self.x = flex_double(np_load(self.restart_file)["x"])

        if rank == 0:
            print("--4 print initial stats")
        rotx, roty, rotz, uc_vals, ncells_vals, scale_vals, _, origZ = self._unpack_internal(self.x, lst_is_x=True)
        if rank == 0 and self.big_dump:

            master_data = {"a": uc_vals[0], "c": uc_vals[1],
                           "Ncells": ncells_vals,
                           "scale": scale_vals,
                           "rotx": rotx,
                           "roty": roty,
                           "rotz": rotz,
                           "origZ": origZ}
            master_data = pandas.DataFrame(master_data)
            master_data["gain"] = self.x[self.gain_xpos]
            print(master_data.to_string())

        # make the parameter masks for isolating parameters of different types
        self._make_parameter_type_selection_arrays()

        #self._setup_resolution_binner()
        # setup the diffBragg instance
        self.D = self.S.D

        if self.refine_Umatrix:
            if self.refine_rotX:
                self.D.refine(0)  # rotX
            if self.refine_rotY:
                self.D.refine(1)  # rotY
            if self.refine_rotZ:
                self.D.refine(2)  # rotZ
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.D.refine(i + 3)  # unit cell params
        if self.refine_ncells:
            self.D.refine(self._ncells_id)
        if self.refine_detdist:
            self.D.refine(self._originZ_id)
        if self.refine_Fcell:
            self.D.refine(self._fcell_id)
        self.D.initialize_managers()

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

    def _make_parameter_type_selection_arrays(self):  # experimental , not really used
        from cctbx.array_family import flex
        self.umatrix_sel = flex.bool(len(self.x), True)
        self.bmatrix_sel = flex.bool(len(self.x), True)
        self.Fcell_sel = flex.bool(len(self.x), True)
        self.origin_sel = flex.bool(len(self.x), True)
        self.spot_scale_sel = flex.bool(len(self.x), True)
        self.ncells_sel = flex.bool(len(self.x), True)
        self.bg_sel = flex.bool(len(self.x), True)
        for i_shot in range(self.n_shots):
            self.umatrix_sel[self.rotX_xpos[i_shot]] = False
            self.umatrix_sel[self.rotY_xpos[i_shot]] = False
            self.umatrix_sel[self.rotZ_xpos[i_shot]] = False

            for i_uc in range(self.n_ucell_param):
                self.bmatrix_sel[self.ucell_xstart[i_shot] + i_uc] = False

            self.ncells_sel[self.ncells_xpos[i_shot]] = False
            self.spot_scale_sel[self.spot_scale_xpos[i_shot]] = False

            self.origin_sel[self.originZ_xpos[i_shot]] = False

            nspots_on_shot = len(self.NANOBRAGG_ROIS[i_shot])
            for i_spot in range(nspots_on_shot):
                self.bg_sel[self.bg_a_xstart[i_shot][i_spot]] = False
                self.bg_sel[self.bg_b_xstart[i_shot][i_spot]] = False
                self.bg_sel[self.bg_c_xstart[i_shot][i_spot]] = False

        for i_fcell in range(self.n_global_fcell):
            self.Fcell_sel[self.fcell_xstart + i_fcell] = False

    def _get_rotX(self, i_shot):
        if self.rescale_params:
            # FIXME ? 
            return self.rotX_sigma*(self.x[self.rotX_xpos[i_shot]]-1) + 0.0
        else:
            return self.x[self.rotX_xpos[i_shot]]

    def _get_rotY(self, i_shot):
        if self.rescale_params:
            return self.rotY_sigma * (self.x[self.rotY_xpos[i_shot]] - 1) + 0.0
        else:
            return self.x[self.rotY_xpos[i_shot]]

    def _get_rotZ(self, i_shot):
        if self.rescale_params:
            return self.rotZ_sigma * (self.x[self.rotZ_xpos[i_shot]] - 1) + 0.0
        else:
            return self.x[self.rotZ_xpos[i_shot]]

    def _get_ucell_vars(self, i_shot):
        all_p = []
        for i in range(self.n_ucell_param):
            if self.rescale_params:
                sig = self.ucell_sigmas[i]
                init = self.ucell_inits[i_shot][i]
                p = sig*(self.x[self.ucell_xstart[i_shot]+i] - 1) + init
            else:
                p = self.x[self.ucell_xstart[i_shot]+i]
            all_p.append(p)
        return all_p

    def _get_originZ_val(self, i_shot):
        val = self.x[self.originZ_xpos[i_shot]]
        if self.rescale_params:
            sig = self.originZ_sigma
            init = self.shot_originZ_init[i_shot]
            val = sig*(val-1) + init
        return val

    def _get_m_val(self, i_shot):
        val = self.x[self.ncells_xpos[i_shot]]
        if self.rescale_params:
            sig = self.m_sigma
            init = self.m_init # note all shots start from same mosaic parameter
            val = np_exp(sig*(val-1))*(init-3) + 3
        else:
            val = np_exp(val)+3
        return val

    def _get_spot_scale(self, i_shot):
        val = self.x[self.spot_scale_xpos[i_shot]]
        if self.rescale_params:
            sig = self.spot_scale_sigma
            init = self.spot_scale_init[i_shot]
            val = np_exp(sig*(val-1))*init
        else:
            val = np_exp(val)
        return val

    def _set_spot_scale(self, new_val, i_shot):
        """just used in testsing derivatives"""
        if self.rescale_params:
            self.spot_scale_init[0] = new_val
            self.x[self.spot_scale_xpos[0]] = 1
        else:
            self.x[self.spot_scale_xpos[0]] = np_log(new_val)

    def _get_bg_vals(self, i_shot, i_spot):
        a_val = self.x[self.bg_a_xstart[i_shot][i_spot]]
        b_val = self.x[self.bg_b_xstart[i_shot][i_spot]]
        c_val = self.x[self.bg_c_xstart[i_shot][i_spot]]
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
                ncells_vals = [self._get_m_val(0)] * len(rotx)
            else:
                ncells_vals = [lst[self.ncells_xpos[0]]] * len(rotx)
        else:
            if lst_is_x:
                ncells_vals = [self._get_m_val(i_shot) for i_shot in range(self.n_shots)]
            else:
                ncells_vals = [lst[self.ncells_xpos[i_shot]] for i_shot in range(self.n_shots)]

        if self.global_originZ_param:
            if lst_is_x:
                originZ_vals = [self._get_originZ_val(0)] * len(rotx)
            else:
                originZ_vals = [lst[self.originZ_xpos[0]]] * len(rotx)
        else:
            if lst_is_x:
                originZ_vals = [self._get_originZ_val(i_shot) for i_shot in range(self.n_shots)]
            else:
                originZ_vals = [lst[self.originZ_xpos[i_shot]] for i_shot in range(self.n_shots)]

        if lst_is_x:
            #ncells_vals = list(np_exp(ncells_vals)+3)
            scale_vals = [self._get_spot_scale(i_shot) for i_shot in range(self.n_shots)]
        else:
            scale_vals = [lst[self.spot_scale_xpos[i_shot]] for i_shot in range(self.n_shots)]
        # this can be used to compare directly
        if self.CRYSTAL_SCALE_TRUTH is not None:
            scale_vals_truths = [self.CRYSTAL_SCALE_TRUTH[i_shot] for i_shot in range(self.n_shots)]
        else:
            scale_vals_truths = None

        # TODO generalize for non tetragonal case
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

        rotx = self._mpi_reduce_broadcast(rotx)
        roty = self._mpi_reduce_broadcast(roty)
        rotz = self._mpi_reduce_broadcast(rotz)
        ncells_vals = self._mpi_reduce_broadcast(ncells_vals)
        scale_vals = self._mpi_reduce_broadcast(scale_vals)
        originZ_vals = self._mpi_reduce_broadcast(originZ_vals)

        ucparams_all = []
        for ucp in ucparams_lsts:
            ucp = self._mpi_reduce_broadcast(ucp)
            ucparams_all.append(ucp)
        if scale_vals_truths is not None:
            scale_vals_truths = self._mpi_reduce_broadcast(scale_vals_truths)

        return rotx, roty, rotz, ucparams_all, ncells_vals, scale_vals, scale_vals_truths, originZ_vals

    def _send_ucell_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i + 3,
                self.UCELL_MAN[self._i_shot].derivative_matrices[i])
            if self.calc_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i + 3, self.UCELL_MAN[self._i_shot].second_derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        (i1, i2), (j1, j2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
        self.D.region_of_interest = (int(i1), int(i2)), (int(j1), int(j2))
        self.D.add_diffBragg_spots()

    def _get_fcell_val(self, i_fcell):
        # TODO vectorize me
        # i_fcell is between 0 and self.n_global_fcell
        # get the asu index and its updated amplitude
        xpos = self.fcell_xstart + i_fcell
        val = self.x[xpos]  # new amplitude
        if self.rescale_params:
            resolution_id = self.res_group_id_from_fcell_index[i_fcell]  # TODO
            sig = self.sigma_for_res_id[resolution_id]*self.fcell_sigma_scale  # TODO
            init = self.fcell_init[i_fcell]
            if self.log_fcells:
                val = np_exp(sig*(val - 1))*init
            else:
                if val < 0:  # NOTE this easily happens without the log c.o.v.
                    self.x[xpos] = 0
                    val = 0
                    self.num_Fcell_kludge += 1
                else:
                    val = sig*(val - 1) + init

        else:
            if self.log_fcells:
                val = np_exp(val)
            if val < 0:  # NOTE this easily happens without the log c.o.v.
                self.x[xpos] = 0
                val = 0
                self.num_Fcell_kludge += 1
        return val

    def _update_Fcell(self):
        idx, data = self.S.D.Fhkl_tuple
        for i_fcell in range(self.n_global_fcell):
            # get the asu miller index
            hkl_asu = self.asu_from_idx[i_fcell]

            new_Fcell_amplitude = self._get_fcell_val(i_fcell)

            # now surgically update the p1 array in nanoBragg with the new amplitudes
            # (need to update each symmetry equivalent)
            equivs = [i.h() for i in miller.sym_equiv_indices(self.space_group, hkl_asu).indices()] # todo: speed test.
            for h_equiv in equivs:
                # get the nanoBragg p1 miller table index corresponding to this hkl equivalent
                try:
                    p1_idx = self.idx_from_p1[h_equiv]  # TODO change name to be more specific
                except KeyError as err:
                    if self.debug:
                        print( h_equiv, err)
                    continue
                data[p1_idx] = new_Fcell_amplitude  # set the data with the new value
        self.S.D.Fhkl_tuple = idx, data  # update nanoBragg again  # TODO: add flag to not re-allocate in nanoBragg!

    def _set_background_plane(self, i_spot):
        xr = self.XREL[self._i_shot][i_spot]
        yr = self.YREL[self._i_shot][i_spot]
        self.a, self.b, self.c = self._get_bg_vals(self._i_shot, i_spot)
        if self.bg_offset_only:
            self.tilt_plane = ONES_LIKE(xr)*self.c
        else:
            self.tilt_plane = xr * self.a + yr * self.b + self.c
        if self.OMEGA_KAHN is not None:
            (i1, i2), (j1, j2) = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
            omega_kahn_correction = self.OMEGA_KAHN[self._panel_id][j1:j2+1, i1:i2+1]
            self.tilt_plane *= omega_kahn_correction

    def _update_rotXYZ(self):
        if self.refine_rotX:
            self.D.set_value(0, self._get_rotX(self._i_shot))
        if self.refine_rotY:
            self.D.set_value(1, self._get_rotY(self._i_shot))
        if self.refine_rotZ:
            self.D.set_value(2, self._get_rotZ(self._i_shot))

    def _update_ncells(self):
        val = self._get_m_val(self._i_shot)
        self.D.set_value(self._ncells_id, val)

    def _update_dxtbx_detector(self):
        # TODO: verify that all panels have same local origin to start with..
        det = self.S.detector
        self.S.panel_id = self._panel_id
        # TODO: select hierarchy level at this point
        # NOTE: what does fast-axis and slow-axis mean for the different hierarchy levels?
        node = det[self._panel_id]
        #from IPython import embed
        #embed()
        #orig = node.get_local_origin()
        #new_originZ = self._get_originZ_val(self._i_shot)
        #new_local_origin = orig[0], orig[1], new_originZ
        #node.set_local_frame(node.get_local_fast_axis(),
        #                     node.get_local_slow_axis(),
        #                     new_local_origin)
        orig = node.get_origin()
        new_originZ = self._get_originZ_val(self._i_shot)
        new_origin = orig[0], orig[1], new_originZ
        node.set_frame(node.get_fast_axis(),
                             node.get_slow_axis(),
                             new_origin)
        self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)
        if self.recenter:
            s0 = self.S.beam.nanoBragg_constructor_beam.get_s0()
            assert ALL_CLOSE(node.get_beam_centre(s0), self.D.beam_center_mm)

    def _extract_Umatrix_derivative_pixels(self):
        self.rotX_dI_dtheta = self.rotY_dI_dtheta = self.rotZ_dI_dtheta = 0
        self.rotX_d2I_dtheta2 = self.rotY_d2I_dtheta2 = self.rotZ_d2I_dtheta2 = 0
        # convenient storage of the gain and scale as a single parameter
        SG = self.scale_fac*self.G2
        if self.refine_Umatrix:
            if self.refine_rotX:
                self.rotX_dI_dtheta = SG*self.D.get_derivative_pixels(0).as_numpy_array()
                if self.calc_curvatures:
                    self.rotX_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(0).as_numpy_array()

            if self.refine_rotY:
                self.rotY_dI_dtheta = SG*self.D.get_derivative_pixels(1).as_numpy_array()
                if self.calc_curvatures:
                    self.rotY_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(1).as_numpy_array()

            if self.refine_rotZ:
                self.rotZ_dI_dtheta = SG*self.D.get_derivative_pixels(2).as_numpy_array()
                if self.calc_curvatures:
                    self.rotZ_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(2).as_numpy_array()

    def _extract_Bmatrix_derivative_pixels(self):
        # the Bmatrix derivatives are stored for each unit cell parameter (UcellManager.variables)
        self.ucell_dI_dtheta = [0] * self.n_ucell_param
        self.ucell_d2I_dtheta2 = [0] * self.n_ucell_param
        SG = self.scale_fac*self.G2
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.ucell_dI_dtheta[i] = SG*self.D.get_derivative_pixels(3 + i).as_numpy_array()
                if self.calc_curvatures:
                    self.ucell_d2I_dtheta2[i] = SG*self.D.get_second_derivative_pixels(3 + i).as_numpy_array()

    def _extract_mosaic_parameter_m_derivative_pixels(self):
        self.m_dI_dtheta = self.m_d2I_dtheta2 = 0
        SG = self.scale_fac * self.G2
        if self.refine_ncells:
            self.m_dI_dtheta = SG*self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()
            if self.calc_curvatures:
                self.m_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(self._ncells_id).as_numpy_array()

    def _extract_originZ_derivative_pixels(self):
        self.detdist_dI_dtheta = self.detdist_d2I_dtheta2 = 0
        SG = self.scale_fac*self.G2
        if self.refine_detdist:
            self.detdist_dI_dtheta = SG*self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()
            if self.calc_curvatures:
                self.detdist_d2I_dtheta2 = SG*self.D.get_second_derivative_pixels(self._originZ_id).as_numpy_array()

    def _extract_Fcell_derivative_pixels(self):
        self.fcell_deriv = self.fcell_second_deriv = 0
        if self.refine_Fcell:
            SG = self.scale_fac*self.G2
            self.fcell_deriv = SG*self.D.get_derivative_pixels(self._fcell_id).as_numpy_array()
            if self.calc_curvatures:
                self.fcell_second_deriv = SG*self.D.get_second_derivative_pixels(self._fcell_id).as_numpy_array()

    def _extract_pixel_data(self):
        self.model_bragg_spots = self.scale_fac*self.D.raw_pixels_roi.as_numpy_array()
        self._extract_Umatrix_derivative_pixels()
        self._extract_Bmatrix_derivative_pixels()
        self._extract_mosaic_parameter_m_derivative_pixels()
        self._extract_originZ_derivative_pixels()
        self._extract_Fcell_derivative_pixels()

    def _update_ucell(self):
        if self.rescale_params:
            pars = self._get_ucell_vars(self._i_shot)
        else:
            _s = slice(self.ucell_xstart[self._i_shot], self.ucell_xstart[self._i_shot] + self.n_ucell_param, 1)
            pars = list(self.x[_s])
        self.UCELL_MAN[self._i_shot].variables = pars
        self._send_ucell_gradients_to_derivative_managers()
        self.D.Bmatrix = self.UCELL_MAN[self._i_shot].B_recipspace

    def _update_umatrix(self):
        self.D.Umatrix = self.CRYSTAL_MODELS[self._i_shot].get_U()

    def _update_beams(self):
        # sim_data instance has a nanoBragg beam object, which takes spectra and converts to nanoBragg xray_beams
        self.S.beam.spectra = self.SPECTRA[self._i_shot]
        self.D.xray_beams = self.S.beam.xray_beams

    def compute_functional_gradients_diag(self):
        self.compute_functional_and_gradients()
        return self._f, self._g, self.d

    #@profile
    def compute_functional_and_gradients(self):
        if self.calc_func:
            if self.verbose:
                self._print_iteration_header()

            if rank == 0 and self.output_dir is not None:
                self._save_state_of_refiner()

            if self.iteratively_freeze_parameters:
                if self.param_sels is None:
                    self.determine_parameter_freeze_order()

            # reset gradient and functional
            self.target_functional = 0
            self.grad = flex_double(self.n)
            if self.calc_curvatures:
                self.curv = flex_double(self.n)

            # current work has these all at 1
            self.gain_fac = self.x[self.gain_xpos]
            self.G2 = self.gain_fac ** 2

            self._update_Fcell()  # update the structure factor with the new x

            if self.CRYSTAL_GT is not None:
                self._initialize_GT_crystal_misorientation_analysis()

            self.image_corr = [0] * len(self.shot_ids)
            for self._i_shot in self.shot_ids:
                if self._i_shot in self.bad_shot_list:
                    continue
                self.scale_fac = self._get_spot_scale(self._i_shot)

                # TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
                self._update_beams()
                self._update_umatrix()
                self._update_ucell()
                self._update_ncells()
                self._update_rotXYZ()
                n_spots = len(self.NANOBRAGG_ROIS[self._i_shot])
                for i_spot in range(n_spots):
                    
                    if self.selection_flags is not None:
                        if self._i_shot not in self.selection_flags:
                            continue
                        elif not self.selection_flags[self._i_shot][i_spot]:
                            continue

                    self._panel_id = int(self.PANEL_IDS[self._i_shot][i_spot])

                    if self.verbose and i_spot % self.spot_print_stride == 0:
                        print("diffBragg: img %d/%d; spot %d/%d; panel %d" \
                              % (self._i_shot + 1, self.n_shots, i_spot + 1, n_spots, self._panel_id)) #, flush=True)

                    self.Imeas = self.ROI_IMGS[self._i_shot][i_spot]
                    self._update_dxtbx_detector()
                    self._run_diffBragg_current(i_spot)
                    self._set_background_plane(i_spot)
                    self._extract_pixel_data()
                    self._evaluate_averageI()

                    # here we can correlate modelLambda with Imeas
                    if self.compute_image_model_correlation:
                        _overlay_corr, _ = pearsonr(self.Imeas.ravel(), self.model_Lambda.ravel())
                    else:
                        _overlay_corr = NAN 
                    self.image_corr[self._i_shot] += _overlay_corr

                    if self.poisson_only:
                        self._evaluate_log_averageI()
                    else:
                        self._evaluate_log_averageI_plus_sigma_readout()

                    #self._max_h_sanity_test()
                    self._derivative_convenience_factors()

                    self.target_functional += self._target_accumulate()

                    # make any plots (this only matters if proper flags have been set)
                    self._show_plots(i_spot, n_spots)

                    # accumulate the per pixel derivatives
                    self._background_derivatives(i_spot)
                    self._Umatrix_derivatives()
                    self._Bmatrix_derivatives()
                    self._mosaic_parameter_m_derivatives()
                    self._originZ_derivatives()
                    self._spot_scale_derivatives()
                    self._gain_factor_derivatives()
                    self._Fcell_derivatives(i_spot)
                    # Done with derivative accumulation

                self.image_corr[self._i_shot] = self.image_corr[self._i_shot] / n_spots

            if has_mpi:
                self.all_image_corr = comm.reduce(self.image_corr, MPI.SUM, root=0)
            else:
                self.all_image_corr = self.image_corr
            # TODO add in the priors:
            self._priors()
            self._parameter_freezes()
            self._mpi_aggregation()

            self._f = self.target_functional
            self._g = self.grad
            self.g = self.grad  # TODO why all these repeated definitions ?, self.g is needed by _verify_diag

            self._curvature_analysis()

            # reset ROI pixels TODO: is this necessary
            self.D.raw_pixels *= 0
            self.gnorm = norm(self.grad)

            if self.verbose:
                if self.CRYSTAL_GT is not None:
                    self._print_GT_crystal_misorientation_analysis()
                self._print_image_correlation_analysis()
                self.print_step()
                self.print_step_grads()

            self.iterations += 1
            self.f_vals.append(self.target_functional)
            time.sleep(self.pause_after_iteration)

            if self.calc_curvatures and not self.use_curvatures:
                if self.num_positive_curvatures == self.use_curvatures_threshold:
                    raise BreakToUseCurvatures

        return self._f, self._g

    def _background_derivatives(self, i_spot):
        if self.refine_background_planes:
            if self.bg_offset_only:  # option to only refine c (t3 in manuscript) plane
                abc_dI_dtheta = [0, 0, self.G2*self.c]
                abc_d2I_dtheta2 = [0, 0, 0]
            else:
                xr = self.XREL[self._i_shot][i_spot]  # fast scan pixels
                yr = self.YREL[self._i_shot][i_spot]  # slow scan pixels
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
                sigma = self.ucell_sigmas[i_ucell]
                d = d*sigma
            derivs.append(d)
        return derivs

    def _get_ucell_second_derivatives(self):
        second_derivs = []
        for i_ucell in range(self.n_ucell_param):
            d2 = self.ucell_d2I_dtheta2[i_ucell]
            if self.rescale_params:
                sigma_squared = self.ucell_sigmas[i_ucell]**2
                d2 = d2*sigma_squared
            second_derivs.append(d2)
        return second_derivs

    def _Bmatrix_derivatives(self):
        if self.refine_Bmatrix:
            # unit cell derivative
            derivs = self._get_ucell_first_derivatives()
            second_derivs = self._get_ucell_second_derivatives()
            for i_ucell in range(self.n_ucell_param):
                xpos = self.ucell_xstart[self._i_shot] + i_ucell
                d = derivs[i_ucell]
                self.grad[xpos] += self._grad_accumulate(d)
                if self.calc_curvatures:
                    d2 = second_derivs
                    self.curv[xpos] += self._curv_accumulate(d, d2)

    def _mosaic_parameter_m_derivatives(self):
        if self.refine_ncells:
            theta = self._get_m_val(self._i_shot)  # mosaic parameter "m"
            theta_minus_three = theta - 3
            sig = self.m_sigma
            if self.rescale_params:
                # case 3 rescaling
                sig_theta_minus_three = sig*theta_minus_three
                d = self.m_dI_dtheta*sig_theta_minus_three
                d2 = self.m_d2I_dtheta2*(sig_theta_minus_three*sig_theta_minus_three) + \
                     self.m_dI_dtheta*(sig*sig_theta_minus_three)
            else:
                # case 4 rescaling with theta_o = 3
                d = self.m_dI_dtheta*theta_minus_three
                d2 = self.m_d2I_dtheta2*(theta_minus_three*theta_minus_three) + self.m_dI_dtheta*theta_minus_three

            xpos = self.ncells_xpos[self._i_shot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)

    def _originZ_derivatives(self):
        if self.refine_detdist:
            if self.rescale_params:
                # case 1 type of rescaling
                d = self.detdist_dI_dtheta*self.originZ_sigma
                d2 = self.detdist_d2I_dtheta2*(self.originZ_sigma*self.originZ_sigma)
            else:
                d = 1*self.detdist_dI_dtheta
                d2 = 1*self.detdist_d2I_dtheta2

            xpos = self.originZ_xpos[self._i_shot]
            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)

    def _Fcell_derivatives(self, i_spot):
        miller_idx = self.ASU[self._i_shot][i_spot]
        multi = self.hkl_frequency[self.idx_from_asu[miller_idx]]
        if self.refine_Fcell and multi >= self.min_multiplicity:
            i_fcell = self.idx_from_asu[self.ASU[self._i_shot][i_spot]]
            xpos = self.fcell_xstart + i_fcell

            self.fcell_dI_dtheta = self.fcell_deriv
            self.fcell_d2I_d2theta2 = self.fcell_second_deriv

            if self.rescale_params:
                fcell = self._get_fcell_val(i_fcell)  # todo: interact with a vectorized object instead
                resolution_id = self.res_group_id_from_fcell_index[i_fcell]
                sig = self.sigma_for_res_id[resolution_id] * self.fcell_sigma_scale
                if self.log_fcells:
                    # case 2 rescaling
                    sig_times_fcell = sig*fcell
                    d = sig_times_fcell*self.fcell_dI_dtheta
                    d2 = (sig_times_fcell*sig_times_fcell)*self.fcell_d2I_d2theta2 + (sig*sig_times_fcell)*self.fcell_dI_dtheta
                else:
                    # case 1 rescaling
                    d = sig*self.fcell_dI_dtheta
                    d2 = (sig*sig)*self.fcell_d2I_d2theta2
            else:
                d = self.fcell_dI_dtheta
                d2 = self.fcell_d2I_d2theta2

            self.grad[xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                self.curv[xpos] += self._curv_accumulate(d, d2)

    def _spot_scale_derivatives(self, return_derivatives=False):
        if self.refine_crystal_scale:
            dI_dtheta = (self.G2/self.scale_fac)*self.model_bragg_spots
            # second derivative is 0 with respect to scale factor
            if self.rescale_params:
                sig = self.spot_scale_sigma
                # case 2 type rescaling
                d = dI_dtheta*self.scale_fac * sig
                d2 = dI_dtheta*(self.scale_fac*sig*sig)
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

    def _gain_factor_derivatives(self):
        if self.refine_gain_fac:
            raise NotImplementedError("gain factor derivatives need more testing")
            d = 2*self.gain_fac*(self.tilt_plane + self.model_bragg_spots)
            self.grad[self.gain_xpos] += self._grad_accumulate(d)
            if self.calc_curvatures:
                d2 = d / self.gain_fac
                self.curv[self.gain_xpos] += self._curv_accumulate(d, d2)

    def _max_h_sanity_test(self, i_spot):
        max_h = tuple(map(int, self.D.max_I_hkl))
        refinement_h = self.ASU[self._i_shot][i_spot]
        equivs = [i.h() for i in miller.sym_equiv_indices(self.space_group, refinement_h).indices()]
        if not max_h in equivs and self.debug:  # TODO understand this more, how does this effect things
           print("Warning max_h  mismatch!!!!!!")



    def _priors(self):
        # experimental, not yet proven to help
        if self.use_ucell_priors and self.refine_Bmatrix:
            for ii in range(self.n_shots):
                for jj in range(self.n_ucell_param):
                    xpos = self.ucell_xstart[ii] + jj
                    ucell_p = self.x[xpos]
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
                    rot_p = self.x[xpos]
                    sig_square = self.sig_rot ** 2
                    self.target_functional += rot_p ** 2 / 2 / sig_square
                    self.grad[xpos] += rot_p / sig_square
                    if self.calc_curvatures:
                        self.curv[xpos] += 1 / sig_square

    def _parameter_freezes(self):
        if self.iteratively_freeze_parameters and self.iterations % self.number_of_frozen_iterations == 0:
            print("\n\n\t\tSwitching!!\n\n")
            freeze_sel = next(self.param_sels)
            self.grad.set_selected(freeze_sel, 0)
            if self.calc_curvatures:
                self.curv.set_selected(freeze_sel, 0)

    def _mpi_aggregation(self):
        # reduce the broadcast summed results:
        if rank == 0:
            print("\nMPI reduce on functionals and gradients...")
        self.target_functional = self._mpi_reduce_broadcast(self.target_functional)
        self.grad = self._mpi_reduce_broadcast(self.grad)
        self.rotx, self.roty, self.rotz, self.uc_vals, self.ncells_vals, self.scale_vals, \
        self.scale_vals_truths, self.origZ_vals = self._unpack_internal(self.x, lst_is_x=True)
        self.Grotx, self.Groty, self.Grotz, self.Guc_vals, self.Gncells_vals, self.Gscale_vals, _, self.GorigZ_vals = \
            self._unpack_internal(self.grad, lst_is_x=False)
        if self.calc_curvatures:
            self.curv = self._mpi_reduce_broadcast(self.curv)
            self.CUrotx, self.CUroty, self.CUrotz, self.CUuc_vals, self.CUncells_vals, self.CUscale_vals, _, self.CUorigZ_vals = \
                self._unpack_internal(self.curv, lst_is_x=False)
        self.tot_fcell_kludge = self._mpi_reduce_broadcast(self.num_Fcell_kludge)

    def _curvature_analysis(self):
        self.tot_neg_curv = 0
        self.neg_curv_shots = []
        if self.calc_curvatures:
            self.tot_neg_curv = sum(self.curv < 0)

        if self.calc_curvatures and not self.use_curvatures:
            if self.tot_neg_curv == 0:
                self.num_positive_curvatures += 1
                self.d = flex_double(self.curv.as_numpy_array())
                self._verify_diag()
            else:
                self.num_positive_curvatures = 0
                self.d = None

        if self.use_curvatures:
            if self.tot_neg_curv == 0:
                self.request_diag_once = False
                self.diag_mode = "always"  # TODO is this proper place to set ?
                self.d = flex_double(self.curv.as_numpy_array())
                self._verify_diag()
            else:
                if self.debug:
                    print("\n\t******************************************")
                    print("\tFREEZING THE CURVATURE: DISASTER AVERSION")
                    print("*\t*****************************************")
        else:
            self.d = None

    def _initialize_GT_crystal_misorientation_analysis(self):
        self.all_ang_off = []
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

            self.all_ang_off.append(ang_off)

        self.bad_shot_list = list(set(self.bad_shot_list))
        if has_mpi:
            self.all_ang_off = comm.gather(self.all_ang_off, root=0)
        self.n_bad_shots = len(self.bad_shot_list)
        if has_mpi:
            self.n_bad_shots = comm.bcast(self.n_bad_shots)

    def _print_GT_crystal_misorientation_analysis(self):
        if has_mpi:
            all_ang_off = [s for sl in self.all_ang_off for s in sl]  # flatten the gathered array
        else:
            all_ang_off = self.all_ang_off
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
            print("\nMissets\n========")
            all_ang_off_s = ["%.5f" % aa for aa in all_ang_off]
            print(", ".join(all_ang_off_s))
            print("N shots deemed bad from missets: %d" % self.n_bad_shots)
        print("MISSETTING median: %.4f; mean: %.4f, max: %.4f, min %.4f, num > .1 deg: %d/%d; num broken=%d"
              % (misset_median, misset_mean, misset_max, misset_min, n_bad_misset, n_misset, n_broken_misset))
        self.all_ang_off = all_ang_off

    def _print_image_correlation_analysis(self):
        all_corr_str = ["%.2f" % ic for ic in self.all_image_corr]
        #if self.print_all_corr:
        print("Correlation stats:")
        print(", ".join(all_corr_str))
        print("---------------")
        n_bad_corr = sum([1 for ic in self.all_image_corr if ic < 0.25])
        print("CORRELATION median: %.4f; mean: %.4f, max: %.4f, min %.4f, num <.25: %d/%d;"
              % (median(self.all_image_corr), mean(self.all_image_corr), max(self.all_image_corr),
                 min(self.all_image_corr), n_bad_corr, len(self.all_image_corr)))

    def _get_refinement_string_label(self):
        refine_str = "refining "
        if self.refine_Fcell:
            refine_str += "fcell, "
        if self.refine_ncells:
            refine_str += "Ncells, "
        if self.refine_Bmatrix:
            refine_str += "Bmat, "
        if self.refine_Umatrix:
            refine_str += "Umat, "
        if self.refine_crystal_scale:
            refine_str += "scale, "
        if self.refine_background_planes:
            refine_str += "bkgrnd, "
        if self.refine_detdist:
            refine_str += "originZ, "
        return refine_str

    def _print_iteration_header(self):
        refine_str = self._get_refinement_string_label()
        border = "<><><><><><><><><><><><><><><><>"
        if self.use_curvatures:

            print(
                "%s%s%s%s\nTrial%d (%s): Compute functional and gradients Iter %d %s(Using Curvatures)%s\n%s%s%s%s"
                % (Bcolors.HEADER, border,border,border, self.trial_id + 1, refine_str, self.iterations + 1, Bcolors.OKGREEN, Bcolors.HEADER, border,border,border, Bcolors.ENDC))
        else:
            print("%s%s%s%s\n, Trial%d (%s): Compute functional and gradients Iter %d PosCurva %d\n%s%s%s%s"
                  % (Bcolors.HEADER, border, border, border, self.trial_id + 1, refine_str, self.iterations + 1, self.num_positive_curvatures, border, border,border, Bcolors.ENDC))

    def _save_state_of_refiner(self):
        outf = os.path.join(self.output_dir, "_fcell_trial%d_iter%d" % (self.trial_id, self.iterations))
        if self.rescale_params:
            fvals = [self._get_fcell_val(i_fcell) for i_fcell in range(self.n_global_fcell)]
            fvals = ARRAY(fvals)
        else:
            fvals = self.x[self.fcell_xstart:self.fcell_xstart + self.n_global_fcell].as_numpy_array()
        SAVEZ(outf, fvals=fvals, x=self.x.as_numpy_array())

    def _show_plots(self, i_spot, n_spots):
        if rank == 0 and self.plot_images and self.iterations % self.plot_stride == 0 and self._i_shot == self.index_of_displayed_image:
            if i_spot % self.plot_spot_stride == 0:
                xr = self.XREL[self._i_shot][i_spot]  # fast scan pixels
                yr = self.YREL[self._i_shot][i_spot]  # slow scan pixels
                if self.plot_residuals:
                    self.ax.clear()
                    residual = self.model_Lambda - self.Imeas
                    x = residual.max()
                    #else:
                    #    x = mean([x, residual.max()])
                    self.ax.plot_surface(xr, yr, residual, rstride=2, cstride=2, alpha=0.3, cmap='coolwarm')
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
                    self.ax1.images[0].set_data(self.model_Lambda)
                    self.ax1.images[0].set_clim(vmin, vmax)
                    self.ax2.images[0].set_data(self.Imeas)
                    self.ax2.images[0].set_clim(vmin, vmax)
                plt.suptitle("Iterations = %d, image %d / %d"
                             % (self.iterations, i_spot + 1, n_spots))
                self.fig.canvas.draw()
                plt.pause(.02)

    def _mpi_reduce_broadcast(self, var):
        if has_mpi:
            var = comm.reduce(var, MPI.SUM, root=0)
            var = comm.bcast(var, root=0)
        return var

    def _poisson_target(self):
        fterm = (self.model_Lambda - self.Imeas * self.log_Lambda).sum()
        return fterm

    def _poisson_d(self, d):
        gterm = (d * self.one_minus_k_over_Lambda).sum()
        return gterm

    def _poisson_d2(self, d, d2):
        cterm = d2 * self.one_minus_k_over_Lambda + d * d * self.k_over_squared_Lambda
        return cterm.sum()

    def _gaussian_target(self):
        #fterm = (self.log2pi + 2*self.log_Lambda_plus_sigma_readout + self.u*self.u*self.one_over_v).sum()
        #return .5*fterm
        fterm = (self.log2pi + 2*self.log_Lambda_plus_sigma_readout + self.u_u_one_over_v).sum()
        return .5*fterm

    def _gaussian_d(self, d):
        gterm = (d*self.one_over_v_times_one_minus_2u_minus_u_squared_over_v).sum()
        #a = self.one_over_v_times_one_minus_2u_minus_u_squared_over_v
        #gterm = numexpr.evaluate('sum(d*a)')
        #local_dict={'a': self.one_over_v_times_one_minus_2u_minus_u_squared_over_v, 'd':d})
        #gterm = 0.5*gterm[()]
        return .5*gterm

    def _gaussian_d2(self, d, d2):
        #cterm = self.one_over_v * (d2*self.one_minus_2u_minus_u_squared_over_v -
        #                           d*d*(self.one_over_v*self.one_minus_2u_minus_u_squared_over_v -
        #                                    (2 + 2*self.u*self.one_over_v + self.u*self.u*self.one_over_v*self.one_over_v)))
        #cterm = .5 * (cterm.sum())
        cterm = self.one_over_v * (d2*self.one_minus_2u_minus_u_squared_over_v -
                                   d*d*(self.one_over_v_times_one_minus_2u_minus_u_squared_over_v -
                                        (2 + 2*self.u_times_one_over_v + self.u_u_one_over_v*self.one_over_v)))
        cterm = .5 * (cterm.sum())
        return cterm

    def _derivative_convenience_factors(self):
        one_over_Lambda = 1. / self.model_Lambda
        self.one_minus_k_over_Lambda = (1. - self.Imeas * one_over_Lambda)
        self.k_over_squared_Lambda = self.Imeas * one_over_Lambda * one_over_Lambda

        self.u = self.Imeas - self.model_Lambda
        self.one_over_v = 1. / (self.model_Lambda + self.sigma_r ** 2)
        self.u_times_one_over_v = self.u*self.one_over_v
        self.u_u_one_over_v = self.u*self.u_times_one_over_v
        self.one_minus_2u_minus_u_squared_over_v = 1 - 2 * self.u - self.u_u_one_over_v
        self.one_over_v_times_one_minus_2u_minus_u_squared_over_v = self.one_over_v*self.one_minus_2u_minus_u_squared_over_v

    def _evaluate_log_averageI(self):  # for Poisson only stats
        # fix log(x<=0)
        try:
            self.log_Lambda = np_log(self.model_Lambda)
        except FloatingPointError:
            pass
        if any((self.model_Lambda <= 0).ravel()):
            self.num_kludge += 1
            is_bad = self.model_Lambda <= 0
            self.log_Lambda[is_bad] = 1e-6
            print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL (kludges=%d)!!!!!!!!!\n<><><><><><><><><>\n" % self.num_kludge)
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda[self.model_Lambda <= 0] = 0

    def _evaluate_log_averageI_plus_sigma_readout(self):
        L = self.model_Lambda + self.sigma_r ** 2
        L_is_neg = (L <= 0).ravel()
        if any(L_is_neg):
            self.num_kludge += 1
            print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda_plus_sigma_readout = np_log(L)
        self.log_Lambda_plus_sigma_readout[L <= 0] = 0  # but will I ever kludge ?

    def print_step(self):
        """Deprecated"""
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))
        rotX = self._get_rotX(self._i_shot) #self.rot_scale*self.x[self.rotX_xpos[self._i_shot]]
        rotY = self._get_rotY(self._i_shot) #self.rot_scale*self.x[self.rotY_xpos[self._i_shot]]
        rotZ = self._get_rotZ(self._i_shot)  #self.rot_scale*self.x[self.rotZ_xpos[self._i_shot]]
        rot_labels = ["rotX=%+3.7g" % rotX, "rotY=%+3.7g" % rotY, "rotZ=%+3.4g" % rotZ]

        if self.refine_Umatrix or self.refine_Bmatrix or self.refine_crystal_scale or self.refine_ncells:
            if self.big_dump:
                master_data = {"Ncells": self.ncells_vals,
                               "scale": self.scale_vals,
                               "rotx": self.rotx,
                               "roty": self.roty,
                               "rotz": self.rotz, "origZ": self.origZ_vals}
                for i_uc in range(self.n_ucell_param):
                    master_data["uc%d" % i_uc] = self.uc_vals[i_uc]


                master_data = pandas.DataFrame(master_data)
                master_data["gain"] = self.x[self.gain_xpos]
                print(master_data.to_string(float_format="%2.7g"))

    def print_step_grads(self):
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for i, (n, v) in enumerate(zip(names, vals)):
            grad = self._g[self.ucell_xstart[self._i_shot] + i]
            ucell_labels.append('G%s=%+2.7g' % (n, grad))

        if self.big_dump:
            master_data ={"GNcells": self.Gncells_vals,
                           "Gscale": self.Gscale_vals,
                           "Grotx": self.Grotx,
                           "Groty": self.Groty,
                           "Grotz": self.Grotz, "GorigZ": self.GorigZ_vals}
            for i_uc in range(self.n_ucell_param):
                master_data["Guc%d" %i_uc]= self.Guc_vals[i_uc]
            master_data = pandas.DataFrame(master_data)
            master_data["Ggain"] = self._g[self.gain_xpos]
            print(master_data.to_string(float_format="%2.7g"))

        if self.calc_curvatures:
            if self.big_dump:
                if self.refine_Umatrix or self.refine_Bmatrix or self.refine_crystal_scale or self.refine_ncells:
                    master_data = {"CUNcells": self.CUncells_vals,
                                   "CUscale": self.CUscale_vals,
                                   "CUrotx": self.CUrotx,
                                   "CUroty": self.CUroty,
                                   "CUrotz": self.CUrotz, "CUorigZ": self.CUorigZ_vals}

                    for i_uc in range(self.n_ucell_param):
                        master_data["CUuc%d" % i_uc] = self.CUuc_vals[i_uc]

                    master_data = pandas.DataFrame(master_data)
                    master_data["CUgain"] = self.curv[self.gain_xpos]
                    if self.big_dump:
                        print(master_data.to_string(float_format="%2.7g"))

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
        scale_stats_string = "SCALE FACTOR STATS: " + ", ".join(scale_stats)
        if self.scale_vals_truths is not None:
            scale_resid = [ABS(s-stru) for s, stru in zip(_sv, self.scale_vals_truths)]
            scale_stats_string += ", truth_resid=%.4f" % median(scale_resid)

        # TODO : use a median for a vals if refining Ncells per shot or unit cell per shot
        scale_stats_string += ", Ncells=%.3f" % self.ncells_vals[0]
        uc_string = ", "
        ucparam_names = self.UCELL_MAN[0].variable_names
        for i_ucparam, ucparam_lst in enumerate(self.uc_vals):
            param_val = median(ucparam_lst)
            uc_string += "%s=%.3f, " % (ucparam_names[i_ucparam], param_val)
        scale_stats_string += uc_string
        scale_stats_string += "originZ=%f, " % median(self.origZ_vals)

        Xnorm = norm(self.x)
        R1 = -1
        R1_i = -1
        ncurv = 0
        if self.calc_curvatures:
            ncurv = len(self.curv)


        if self.Fref is not None and self.iterations % self.merge_stat_frequency == 0:
            self.R_overall = self.Fobs_Fref_Rfactor(use_binning=False, auto_scale=self.scale_r1)
            self.CC_overall = self.Fobs.correlation(self.Fref_aligned).coefficient()
            print("R-factor overall: %.4f, CC overall: %.4f" % (self.R_overall, self.CC_overall))
            if self.print_resolution_bins:
                print("R-factor (shells):")
                print(self.Fobs_Fref_Rfactor(use_binning=True, auto_scale=self.scale_r1).show())
                print("CC (shells):")
                self.Fobs.correlation(self.Fref_aligned, use_binning=True).show()

        print(
            "%s\n\t%s|G|=%2.7g, eps*|X|=%2.7g,%s R1=%2.7g (R1 at start=%2.7g), Fcell kludges=%d, Neg. Curv.: %d/%d on shots=%s\n"
            % (scale_stats_string, Bcolors.OKBLUE, self.gnorm, Xnorm * self.trad_conv_eps, Bcolors.ENDC, R1, R1_i, self.tot_fcell_kludge, self.tot_neg_curv, ncurv,
               ", ".join(map(str, self.neg_curv_shots))))
        #print("<><><><><><><><> TOP GUN <><><><><><><><>")
        #print("                 End of iteration.")
        if self.testing_mode:
            self.conv_test()

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
            rx = ry = rz = 0
            if self.refine_Umatrix:
                if self.refine_rotX:
                    rx = self._get_rotX(i_shot)
                if self.refine_rotY:
                    ry = self._get_rotY(i_shot)
                if self.refine_rotZ:
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
            self.r1_scale = minimize(FatRefiner._rfactor_minimizer_target,
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
        from cctbx import crystal
        # TODO update for non global unit cell case (average over unit cells)

        self._update_Fcell()  # just in case update the Fobs

        um = self.UCELL_MAN[0]
        sym = crystal.symmetry(unit_cell=um.unit_cell_parameters, space_group_symbol=self.symbol)
        mset_obs = miller.set(sym, self.Fobs.indices(), anomalous_flag=True)
        fobs = miller.array(mset_obs, self.Fobs.data()).set_observation_type_xray_amplitude()
        # TODO: what to do in MPI mode when writing ?
        if save_to_file is not None and rank == 0:
            fobs.as_mtz_dataset(column_root_label='fobs', wavelength=wavelength).mtz_object().write(save_to_file)
        return fobs

    def conv_test(self):
        err = []
        s = ""
        A = []
        vars = self._get_ucell_vars(0)
        for i in range(self.n_ucell_param):
            if self.rescale_params:
                a = vars[i]
            else:
                a = self.x[self.ucell_xstart[0] + i]

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
            ncells_val = self._get_m_val(i_shot)  #np.exp(self.x[self.ncells_xpos[i_shot]]) + 3
            ncells_resid = abs(ncells_val - self.gt_ncells)

            if mn_err < 0.01 and ang_off < 0.004 and ncells_resid < 0.1:
                shot_refined.append(True)
            else:
                shot_refined.append(False)

            if self.refine_detdist:
                det_resid = abs(self.originZ_gt[i_shot] - self._get_originZ_val(i_shot))
                all_det_resid.append(det_resid)
                out_str += ", OrigZ resid = %.4f" % det_resid

            if self.refine_ncells:
                out_str += ", ncells resid=%.4f" % ncells_resid
            print(out_str)

        if all(shot_refined):
            if self.refine_detdist:
                if all([det_resid < 0.01 for det_resid in all_det_resid]):
                    print("OK")
                    exit()
            else:
                print("OK")
                exit()
