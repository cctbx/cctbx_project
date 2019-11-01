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

if rank == 0:
    import pandas
    from numpy import mean, unique
    from numpy import log as np_log
    from numpy import all as np_all
    from numpy.linalg import norm

    from simtbx.diffBragg.refiners import BreakToUseCurvatures
    from scitbx.array_family import flex

    flex_double = flex.double
    from scitbx.matrix import col
    from simtbx.diffBragg.refiners import PixelRefinement

else:
    mean = unique = np_log = np_all = norm = None
    # stdout_flush = None
    BreakToUseCurvatures = None
    flex_double = None
    col = None
    PixelRefinement = None

if has_mpi:
    mean = comm.bcast(mean, root=0)
    unique = comm.bcast(unique, root=0)
    np_log = comm.bcast(np_log, root=0)
    norm = comm.bcast(norm, root=0)
    np_all = comm.bcast(np_all, root=0)
    # stdout_flush = comm.bcast(stdout_flush)
    BreakToUseCurvatures = comm.bcast(BreakToUseCurvatures)
    flex_double = comm.bcast(flex_double)
    col = comm.bcast(col)
    PixelRefinement = comm.bcast(PixelRefinement)

if rank == 0:
    import pylab as plt
    from IPython import embed

import sys
import warnings
from copy import deepcopy
from cxid9114.helpers import compare_with_ground_truth
from cctbx import miller, sgtbx

warnings.filterwarnings("ignore")


class FatRefiner(PixelRefinement):

    def __init__(self, n_total_params, n_local_params, n_global_params, local_idx_start,
                 shot_ucell_managers, shot_rois, shot_nanoBragg_rois,
                 shot_roi_imgs, shot_spectra, shot_crystal_GTs,
                 shot_crystal_models, shot_xrel, shot_yrel, shot_abc_inits, shot_asu,
                 global_param_idx_start,
                 shot_panel_ids, init_gain=1, init_scale=1):
        PixelRefinement.__init__(self)
        # super(GlobalFat, self).__init__()
        self.num_kludge = 0

        # dictionaries whose keys are the shot indices
        self.UCELL_MAN = shot_ucell_managers
        # cache shot ids and make sure they are identical in all other input dicts
        self.shot_ids = sorted(shot_ucell_managers.keys())
        self.n_shots = len(self.shot_ids)
        # sanity check: no repeats of the same shot
        assert len(self.shot_ids) == len(set(self.shot_ids))

        self.ROIS = self._check_keys(shot_rois)
        self.ASU = self._check_keys(shot_asu)
        self.NANOBRAGG_ROIS = self._check_keys(shot_nanoBragg_rois)
        self.ROI_IMGS = self._check_keys(shot_roi_imgs)
        self.SPECTRA = self._check_keys(shot_spectra)
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
        self.multi_panel = True  # we are multi panel
        self.f_vals = []  # store the functional over time

        # start with the first shot
        self._i_shot = self.shot_ids[0]

        # These are the per-shot parameters
        self.n_rot_param = 3
        self.n_ucell_param = len(self.UCELL_MAN[self._i_shot].variables)
        self.n_ncells_param = 1
        self.n_spot_scale_param = 1
        self.n_per_shot_params = self.n_rot_param + self.n_ncells_param + self.n_spot_scale_param

        assert self.n_per_shot_params * self.n_shots == self.n_local_params

        self._ncells_id = 9  # diffBragg internal index for Ncells derivative manager
        self._originZ_id = 10  # diffBragg internal index for originZ derivative manager
        self._fcell_id = 11  # diffBragg internal index for Fcell derivative manager
        self._init_scale = init_scale
        self._init_gain = init_gain
        self.num_positive_curvatures = 0
        self._panel_id = None
        self.symbol = "P43212"

        self.pid_from_idx = {}
        self.idx_from_pid = {}

        self.idx_from_asu = {}
        self.asu_from_idx = {}

        # For priors, use these ..  experimental
        self.ave_ucell = [78.95, 38.12]  ## Angstrom
        self.sig_ucell = [0.025, 0.025]
        self.sig_rot = 0.01  # radian

        # where the global parameters being , initially just gain and detector distance
        self.global_param_idx_start = global_param_idx_start

        self.a = self.b = self.c = None  # tilt plan place holder

    def setup_plots(self):
        if rank == 0:
            if self.plot_statistics:
                print ("plot_stats")
                assert not self.plot_images  # only one plot at a time for now
                self.fig, self.ax = plt.subplots(nrows=2, ncols=2)
                self.ax1 = self.ax[0][0]
                self.ax2 = self.ax[0][1]
                self.ax3 = self.ax[1][0]
                self.ax4 = self.ax[1][1]

            if self.plot_images:
                if self.plot_residuals:
                    self.fig = plt.figure()
                    self.ax = self.fig.gca(projection='3d')
                    self.ax.set_yticklabels([])
                    self.ax.set_xticklabels([])
                    self.ax.set_zticklabels([])
                    self.ax.set_zlabel("model residual")
                    self.ax.set_facecolor("gray")
                else:
                    self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1, ncols=2)
                    self.ax1.imshow([[0, 1, 1], [0, 1, 2]])
                    self.ax2.imshow([[0, 1, 1], [0, 1, 2]])

    def __call__(self, *args, **kwargs):
        _, _ = self.compute_functional_and_gradients()
        return self.x, self._f, self._g, self.d

    @property
    def n(self):
        """LBFGS property"""
        return self.n_total_params

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
            self.gain_fac * self.gain_fac * (self.tilt_plane + self.scale_fac * self.scale_fac * self.model_bragg_spots)

    def _set_background_plane(self, i_spot):
        xr = self.XREL[self._i_shot][i_spot]
        yr = self.YREL[self._i_shot][i_spot]
        self.tilt_plane = xr * self.a + yr * self.b + self.c

    def _setup(self):
        # Here we go!  https://youtu.be/7VvkXA6xpqI
        if not self.asu_from_idx:
            raise ValueError("Need to supply a non empty asu from idx map")
        if not self.idx_from_asu:  # # TODO just derive from its inverse
            raise ValueError("Need to supply a non empty idx from asu map")

        # get the Fhkl information from P1 array internal to nanoBragg
        if comm.rank==0:
            print ("Making the big ones")
        idx, data = self.S.D.Fhkl_tuple
        self.idx_from_p1 = {h: i for i, h in enumerate(idx)}
        #self.p1_from_idx = {i: h for i, h in zip(idx, data)}
        if comm.rank==0:
            print("Done")

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        # Make te global sized parameter array, though here we only update the local portion
        self.x = flex_double(self.n_total_params)

        self.rotX_xpos = {}
        self.rotY_xpos = {}
        self.rotZ_xpos = {}
        self.ucell_xstart = {}
        self.ncells_xpos = {}
        self.spot_scale_xpos = {}
        self.n_panels = {}

        for i_shot in self.shot_ids:
            self.pid_from_idx[i_shot] = {i: pid for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.idx_from_pid[i_shot] = {pid: i for i, pid in enumerate(unique(self.PANEL_IDS[i_shot]))}
            self.n_panels[i_shot] = len(self.pid_from_idx[i_shot])

            self.rotX_xpos[i_shot] = self.local_idx_start + i_shot * self.n_per_shot_params
            self.rotY_xpos[i_shot] = self.rotX_xpos[i_shot] + 1
            self.rotZ_xpos[i_shot] = self.rotY_xpos[i_shot] + 1

            self.x[self.rotX_xpos[i_shot]] = 0
            self.x[self.rotY_xpos[i_shot]] = 0
            self.x[self.rotZ_xpos[i_shot]] = 0

            self.ucell_xstart[i_shot] = self.global_param_idx_start  # TODO: make global ucell params optional
            self.ncells_xpos[i_shot] = self.rotZ_xpos[i_shot] + 1

            #if self.refine_global_unit_cell:
            #    self.ucell_xstart[i_shot] = self.global_param_idx_start
            #    self.ncells_xpos[i_shot] = self.rotZ_xpos[i_shot] + 1  # self.n_ucell_param

            #else:
            #    self.ucell_xstart[i_shot] = self.rotZ_xpos[i_shot] + 1
            #    # populate the x-array with initial values
            #    for i_uc in range(self.n_ucell_param):
            #        self.x[self.ucell_xstart[i_shot] + i_uc] = self.UCELL_MAN[i_shot].variables[i_uc]

            #    # put in Ncells abc estimate
            #    self.ncells_xpos[i_shot] = self.ucell_xstart[i_shot] + self.n_ucell_param

            self.x[self.ncells_xpos[i_shot]] = self.S.crystal.Ncells_abc[0]  # NOTE: each shot gets own starting Ncells

            self.spot_scale_xpos[i_shot] = self.ncells_xpos[i_shot] + 1
            self.x[self.spot_scale_xpos[i_shot]] = self._init_scale  # TODO: each shot gets own starting scale factor

        #if self.refine_global_unit_cell:
        self.fcell_xstart = self.global_param_idx_start + self.n_ucell_param
        #else:
        #    self.fcell_xstart = self.global_param_idx_start
        self.originZ_xpos = self.n_total_params - 2
        self.gain_xpos = self.n_total_params - 1
        if rank == 0:
            # put in estimates for origin vectors
            # TODO: refine at the different hierarchy
            # get te first Z coordinate for now..
            # print("Setting origin: %f " % self.S.detector[0].get_local_origin()[2])
            # TODO if self.refine_global_unit_cell
            for i_cell in range(self.n_ucell_param):
                self.x[self.ucell_xstart[0] + i_cell] = self.UCELL_MAN[0].variables[i_cell]

            F = self.S.D.Fhkl  # initial values, this table is the high symm table expanded to P1
            # NOTE: dont ever set D.Fhkl property in the refinement setting, it tries to update the unit cell
            for i_fcell in range(self.n_global_fcell):
                asu_hkl = self.asu_from_idx[i_fcell]
                self.x[self.fcell_xstart + i_fcell] = F.value_at_index(asu_hkl)

            # set det dist
            self.x[self.originZ_xpos] = self.S.detector[0].get_local_origin()[2]  # NOTE maybe just origin instead?
            # set gain TODO: remove gain from all of this and never refine it
            self.x[self.gain_xpos] = self._init_gain  # gain factor

            # TODO per panel gain corretion / pedestal correction?
            # n_panels = len(self.S.detector)
            # self.origin_xstart = self.global_param_idx_start
            # for i_pan in range(n_panels):
            #    pid = self.pid_from_idx[i_pan]
            #    self.x[self.origin_xstart + i_pan] = self.S.detector[pid].get_local_origin()[2]
            # lastly, the panel gain correction factor
            # for i_pan in range(n_panels):
            #    self.x[-1] = self._init_gain
            # self.x[-1] = self._init_scale  # initial scale factor
        # reduce then broadcast self.x
        self.x = self._mpi_reduce_broadcast(self.x)
        # self.x = comm.reduce(self.x, MPI.SUM, root=0)
        # self.x = comm.bcast(self.x, root=0)
        if comm.rank == 0:
            rotx, roty, rotz, a_vals, c_vals, ncells_vals, scale_vals = self._unpack_internal(self.x)

            master_data = {"a": a_vals, "c": c_vals,
                           "Ncells": ncells_vals,
                           "scale": scale_vals,
                           "rotx": rotx,
                           "roty": roty,
                           "rotz": rotz}
            master_data = pandas.DataFrame(master_data)

            #print("Avals")
            #print(a_vals)
            #print("Cvals")
            #print(c_vals)
            #print("Ncellsvals")
            #print(ncells_vals)
            #print("Scalevals")
            #print(scale_vals)
            #print("Origin Z")
            #print(self.x[self.originZ_xpos])
            #print ("Gain val")
            #print(self.x[self.gain_xpos])

            master_data["gain"] = self.x[self.gain_xpos]
            master_data["originZ"] = self.x[self.originZ_xpos]
            print(master_data.to_string())
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

    def _unpack_internal(self, lst):
        # NOTE: This is not a generalized method, only works in context of D9114 global refinement
        #x = self..as_numpy_array()
        rotx = lst[:-self.n_global_params:self.n_per_shot_params]
        roty = lst[1:-self.n_global_params:self.n_per_shot_params]
        rotz = lst[2:-self.n_global_params:self.n_per_shot_params]
        #a_vals = lst[3:-self.n_global_params:self.n_per_shot_params]
        #c_vals = lst[4:-self.n_global_params:self.n_per_shot_params]
        ncells_vals = lst[3:-self.n_global_params:self.n_per_shot_params]
        scale_vals = lst[4:-self.n_global_params:self.n_per_shot_params]

        a, c = lst[self.ucell_xstart[0]:self.ucell_xstart[0] + self.n_ucell_param]
        a_vals = [a]*len(rotx)
        c_vals = [c]*len(rotx)

        return rotx, roty, rotz, a_vals, c_vals, ncells_vals, scale_vals

    def _send_gradients_to_derivative_managers(self):
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
        self.D.region_of_interest = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
        self.D.add_diffBragg_spots()

    def _update_Fcell(self):
        idx, data = self.S.D.Fhkl_tuple
        for i_fcell in range(self.n_global_fcell):
            # get the asu index and its updated amplitude
            hkl = self.asu_from_idx[i_fcell]
            xpos = self.fcell_xstart + i_fcell
            new_Fcell = self.x[xpos]  # new amplitude

            #FIXME with a re-parameterization
            if new_Fcell < 0:
                self.x[xpos] = 0
                new_Fcell = 0
                self.num_kludge += 1

            # now surgically update the p1 array in nanoBragg with the new amplitudes
            # need to update each symmetry equivalent
            sg = sgtbx.space_group(sgtbx.space_group_info(symbol=self.symbol).type().hall_symbol())
            equivs = [i.h() for i in miller.sym_equiv_indices(sg, hkl).indices()]
            for h_equiv in equivs:
                # get the nanoBragg p1 miller table index corresponding to this hkl equivalent
                try:
                    p1_idx = self.idx_from_p1[h_equiv]  # TODO change name to be more specific
                except KeyError as err:
                    continue
                data[p1_idx] = new_Fcell  # set the data with the new value
        self.S.D.Fhkl_tuple = idx, data  # update nanoBragg again  # TODO: add flag to not re-allocate in nanoBragg!
        #if comm.rank==0:
        #    print("Done updating Fcell!")

    def _update_rotXYZ(self):
        if self.refine_rotX:
            self.D.set_value(0, self.x[self.rotX_xpos[self._i_shot]])
        if self.refine_rotY:
            self.D.set_value(1, self.x[self.rotY_xpos[self._i_shot]])
        if self.refine_rotZ:
            self.D.set_value(2, self.x[self.rotZ_xpos[self._i_shot]])

    def _update_ncells(self):
        self.D.set_value(self._ncells_id, self.x[self.ncells_xpos[self._i_shot]])

    def _update_dxtbx_detector(self):
        # TODO: verify that all panels have same local origin to start with..
        det = self.S.detector
        self.S.panel_id = self._panel_id
        # TODO: select hierarchy level at this point
        # NOTE: what does fast-axis and slow-axis mean
        # for the different hierarchy levels?
        node = det[self._panel_id]
        orig = node.get_local_origin()
        new_originZ = self.x[self.originZ_xpos]
        new_local_origin = orig[0], orig[1], new_originZ
        node.set_local_frame(node.get_local_fast_axis(),
                             node.get_local_slow_axis(),
                             new_local_origin)
        self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)

    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        self.rot_second_deriv = [0, 0, 0]
        if self.refine_Umatrix:
            if self.refine_rotX:
                self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[0] = self.D.get_second_derivative_pixels(0).as_numpy_array()
            if self.refine_rotY:
                self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[1] = self.D.get_second_derivative_pixels(1).as_numpy_array()
            if self.refine_rotZ:
                self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[2] = self.D.get_second_derivative_pixels(2).as_numpy_array()

        self.ucell_derivatives = [0] * self.n_ucell_param
        self.ucell_second_derivatives = [0] * self.n_ucell_param
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.ucell_derivatives[i] = self.D.get_derivative_pixels(3 + i).as_numpy_array()
                if self.calc_curvatures:
                    self.ucell_second_derivatives[i] = self.D.get_second_derivative_pixels(3 + i).as_numpy_array()

        self.ncells_deriv = self.detdist_deriv = self.fcell_deriv = 0
        self.ncells_second_deriv = self.detdist_second_deriv = self.fcell_second_deriv = 0
        if self.refine_ncells:
            self.ncells_deriv = self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()
            if self.calc_curvatures:
                self.ncells_second_deriv = self.D.get_second_derivative_pixels(self._ncells_id).as_numpy_array()
        if self.refine_detdist:
            self.detdist_deriv = self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()
            if self.calc_curvatures:
                self.detdist_second_deriv = self.D.get_second_derivative_pixels(self._originZ_id).as_numpy_array()
        if self.refine_Fcell:
            self.fcell_deriv = self.D.get_derivative_pixels(self._fcell_id).as_numpy_array()
            if self.calc_curvatures:
                self.fcell_second_deriv = self.D.get_second_derivative_pixels(self._fcell_id).as_numpy_array()

        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a, self.b, self.c = self.ABC_INIT[self._i_shot][i_spot]

    def _update_ucell(self):
        _s = slice(self.ucell_xstart[self._i_shot], self.ucell_xstart[self._i_shot] + self.n_ucell_param, 1)
        pars = list(self.x[_s])
        self.UCELL_MAN[self._i_shot].variables = pars
        self._send_gradients_to_derivative_managers()
        self.D.Bmatrix = self.UCELL_MAN[self._i_shot].B_recipspace

    def _update_umatrix(self):
        self.D.Umatrix = self.CRYSTAL_MODELS[self._i_shot].get_U()

    def _update_beams(self):
        # sim_data instance has a nanoBragg beam object, which takes spectra and converts to nanoBragg xray_beams
        self.S.beam.spectra = self.SPECTRA[self._i_shot]
        self.D.xray_beams = self.S.beam.xray_beams

    def compute_functional_and_gradients(self):
        if self.calc_func:
            if self.verbose:
                if self.use_curvatures:
                    print("Compute functional and gradients Iter %d (Using Curvatures)\n<><><><><><><><><><><><><>"
                          % (self.iterations + 1))
                else:
                    print("Compute functional and gradients Iter %d PosCurva %d\n<><><><><><><><><><><><><>"
                          % (self.iterations + 1, self.num_positive_curvatures))
            #print("Rank=%d, Number of kludges: %d" % (comm.rank, self.num_kludge))
            f = 0
            g = flex_double(self.n)
            if self.calc_curvatures:
                self.curv = flex_double(self.n)

            self.gain_fac = self.x[self.gain_xpos]
            G2 = self.gain_fac ** 2
            self._update_Fcell()  # update the structure factor with the new x  # TODO do I work ?

            for self._i_shot in self.shot_ids:
                self.scale_fac = self.x[self.spot_scale_xpos[self._i_shot]]
                S2 = self.scale_fac ** 2
                # TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
                self._update_beams()
                self._update_umatrix()
                self._update_ucell()
                self._update_ncells()
                self._update_rotXYZ()
                n_spots = len(self.NANOBRAGG_ROIS[self._i_shot])
                for i_spot in range(n_spots):
                    self._panel_id = self.PANEL_IDS[self._i_shot][i_spot]
                    if self.verbose:
                        print "\rdiffBragg: img %d/%d; spot %d/%d; panel %d" \
                              % (self._i_shot + 1, self.n_shots, i_spot + 1, n_spots, self._panel_id),
                        sys.stdout.flush()

                    self._update_dxtbx_detector()
                    self._run_diffBragg_current(i_spot)
                    self._unpack_bgplane_params(i_spot)
                    self._set_background_plane(i_spot)
                    self._extract_pixel_data()
                    self._evaluate_averageI()
                    if self.poisson_only:
                        self._evaluate_log_averageI()
                    else:
                        self._evaluate_log_averageI_plus_sigma_readout()

                    self.Imeas = self.ROI_IMGS[self._i_shot][i_spot]

                    # helper terms for doing derivatives
                    one_over_Lambda = 1. / self.model_Lambda
                    self.one_minus_k_over_Lambda = (1. - self.Imeas * one_over_Lambda)
                    self.k_over_squared_Lambda = self.Imeas * one_over_Lambda * one_over_Lambda

                    self.u = self.Imeas - self.model_Lambda
                    self.one_over_v = 1. / (self.model_Lambda + self.sigma_r ** 2)
                    self.one_minus_2u_minus_u_squared_over_v = 1 - 2 * self.u - self.u * self.u * self.one_over_v

                    f += self._target_accumulate()

                    if self.plot_images and self.iterations % self.plot_stride == 0:
                        xr = self.XREL[self._i_shot][i_spot]  # fast scan pixels
                        yr = self.YREL[self._i_shot][i_spot]  # slow scan pixels
                        if self.plot_residuals:
                            self.ax.clear()
                            residual = self.model_Lambda - self.Imeas
                            if i_spot == 0:
                                x = residual.max()
                            else:
                                x = mean([x, residual.max()])

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
                    if self.refine_Umatrix:
                        x_positions = [self.rotX_xpos[self._i_shot],
                                       self.rotY_xpos[self._i_shot],
                                       self.rotZ_xpos[self._i_shot]]
                        for ii, xpos in enumerate(x_positions):
                            d = S2 * G2 * self.rot_deriv[ii]
                            g[xpos] += self._grad_accumulate(d)
                            if self.calc_curvatures:
                                d2 = S2 * G2 * self.rot_second_deriv[ii]
                                self.curv[xpos] +=self._curv_accumulate(d, d2)

                    if self.refine_Bmatrix:
                        # unit cell derivative
                        for i_ucell_p in range(self.n_ucell_param):
                            xpos = self.ucell_xstart[self._i_shot] + i_ucell_p
                            d = S2 * G2 * self.ucell_derivatives[i_ucell_p]
                            g[xpos] += self._grad_accumulate(d)

                            if self.calc_curvatures:
                                d2 = S2 * G2 * self.ucell_second_derivatives[i_ucell_p]
                                self.curv[xpos] += self._curv_accumulate(d, d2)

                    if self.refine_ncells:
                        d = S2*G2*self.ncells_deriv
                        xpos = self.ncells_xpos[self._i_shot]
                        g[xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = S2 * G2 * self.ncells_second_deriv
                            self.curv[xpos] += self._curv_accumulate(d, d2)

                    if self.refine_detdist:
                        raise NotImplementedError("Cannot refined detdist (yet...)")
                        # if self.calc_curvatures:
                        #    raise NotImplementedError("Cannot use curvatures and refine detdist (yet...)")
                        # origin_xpos = self.origin_xstart + self.idx_from_pid[self._panel_id]
                        # g[origin_xpos] += (S2*G2*self.detdist_deriv*one_minus_k_over_Lambda).sum()

                    if self.refine_Fcell:
                        d = S2*G2*self.fcell_deriv
                        xpos = self.fcell_xstart + self.idx_from_asu[self.ASU[self._i_shot][i_spot]]
                        g[xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = S2*G2*self.fcell_second_deriv
                            self.curv[xpos] += self._curv_accumulate(d, d2)

                    if self.refine_crystal_scale:
                        d = G2 * 2 * self.scale_fac * self.model_bragg_spots
                        xpos = self.spot_scale_xpos[self._i_shot]
                        g[xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = d / self.scale_fac
                            self.curv[xpos] += self._curv_accumulate(d, d2)

                    if self.refine_gain_fac:
                        d = 2 * self.gain_fac * (self.tilt_plane + S2 * self.model_bragg_spots)
                        g[self.gain_xpos] += self._grad_accumulate(d)
                        if self.calc_curvatures:
                            d2 = d / self.gain_fac
                            self.curv[self.gain_xpos] += self._curv_accumulate(d, d2)

            # TODO add in the priors:
            if self.use_ucell_priors and self.refine_Bmatrix:
                for ii in range(self.n_shots):
                    for jj in range(self.n_ucell_param):
                        xpos = self.ucell_xstart[ii] + jj
                        ucell_p = self.x[xpos]
                        sig_square = self.sig_ucell[jj]**2
                        f += (ucell_p - self.ave_ucell[jj])**2 / 2 / sig_square
                        g[xpos] += (ucell_p-self.ave_ucell[jj])/sig_square
                        if self.calc_curvatures:
                            self.curv[xpos] += 1/sig_square

            if self.use_rot_priors and self.refine_Umatrix:
                for ii in range(self.n_shots):
                    x_positions = [self.rotX_xpos[self._i_shot],
                                   self.rotY_xpos[self._i_shot],
                                   self.rotZ_xpos[self._i_shot]]
                    for xpos in x_positions:
                        rot_p = self.x[xpos]
                        sig_square = self.sig_rot**2
                        f += rot_p**2 / 2 / sig_square
                        g[xpos] += rot_p / sig_square
                        if self.calc_curvatures:
                            self.curv[xpos] += 1/sig_square

            # reduce the broadcast summed results:
            f = self._mpi_reduce_broadcast(f)
            g = self._mpi_reduce_broadcast(g)
            if self.calc_curvatures:
                self.curv = self._mpi_reduce_broadcast(self.curv)

            if comm.rank == 0:
                if self.plot_statistics:
                    rotx, roty, rotz, a, c, _, _ = self._unpack_internal(self.x)
                    self.ax1.clear()
                    self.ax2.clear()
                    self.ax3.clear()
                    self.ax4.clear()
                    # plot unit cell a
                    self.ax1.hist(a, bins='auto')
                    # plot unit cell c
                    self.ax2.hist(c, bins='auto')
                    # plot restoring missets
                    self.ax3.hist(rotx, bins='auto', histtype='step', label='x')
                    self.ax3.hist(roty, bins='auto', histtype='step', label='y')
                    self.ax3.hist(rotz, bins='auto', histtype='step', label='z')
                    # plot grad norms per shot as a bar plot
                    assert (self.n_total_params-2) % self.n_per_shot_params == 0
                    total_shots = int((self.n_total_params - 2) / self.n_per_shot_params )
                    g_norm_per_shot = [norm(g[self.n_per_shot_params * i: self.n_per_shot_params * (i + 1)])
                                       for i in range(total_shots)]

                    curv_per_shot = [self.curv[self.n_per_shot_params*i: self.n_per_shot_params * (i + 1)]
                                       for i in range(total_shots)]
                    idx_all_pos, idx_has_neg = [], []
                    for i in range(total_shots):
                        if not np_all(curv_per_shot[i].as_numpy_array() >=0 ):
                            idx_has_neg.append(i)
                        else:
                            idx_all_pos.append(i)

                    self.ax4.bar(idx_has_neg,
                                 height=[g_norm_per_shot[i] for i in idx_has_neg],
                                 width=0.8, color='r')
                    self.ax4.bar(idx_all_pos,
                                 height=[g_norm_per_shot[i] for i in idx_all_pos],
                                 width=0.8, color='g')
                    self.ax4.set_yscale("log")
                    self.fig.canvas.draw()
                    plt.pause(.02)
            self._f = f
            self._g = g
            self.g = g

            if self.calc_curvatures and not self.use_curvatures:
                if np_all(self.curv.as_numpy_array() >= 0):
                    self.num_positive_curvatures += 1
                    self.d = flex_double(self.curv.as_numpy_array())
                    self._verify_diag()
                else:
                    self.num_positive_curvatures = 0
                    self.d = None

            if self.use_curvatures:
                if np_all(self.curv.as_numpy_array() >= 0):
                    self.request_diag_once = False
                    self.d = flex_double(self.curv.as_numpy_array())
                    self._verify_diag()
                else:
                    print("\n**************************************")
                    print("**************************************")
                    print("**************************************")
                    print("**************************************")
                    print("\tFIXING CURVA: ATTEMPTING DISASTER AVERSION")
                    print("*nn***************************************")
                    print("*nn***************************************")
                    print("*nn***************************************")
                    print("*nn***************************************")
            else:
                self.d = None

            #        self.num_positive_curvatures += 1
            #    else:
            #        self.num_positive_curvatures = 0
            #    if self.num_positive_curvatures == self.update_curvatures_every:
            #        print("Updating curva because looks good!!")
            #        self.d = self.curv
            #        self._verify_d()  # derived method
            #        self.num_positive_curvatures = 0

            self.D.raw_pixels *= 0
            gnorm = norm(g)
            all_ang_off = []
            for i in range(self.n_shots):
                try:
                    Ctru = self.CRYSTAL_GT[i]
                    atru, btru, ctru = Ctru.get_real_space_vectors()
                    ang, ax = self.get_correction_misset(as_axis_angle_deg=True, i_shot=i)
                    B = self.get_refined_Bmatrix(i)
                    C = deepcopy(self.CRYSTAL_MODELS[i])
                    C.set_B(B)
                    C.rotate_around_origin(ax, ang)
                    ang_off = compare_with_ground_truth(atru, btru, ctru,
                                                        [C],
                                                        symbol="P43212")[0]
                except Exception:
                    ang_off = 9999
                all_ang_off.append(ang_off)

            all_ang_off = comm.gather(all_ang_off, root=0)
            if self.verbose:
                print("\nMissets\n========")
                all_ang_off = ["%.5f" % s for sl in all_ang_off for s in sl]
                print(", ".join(all_ang_off))
                self.print_step("LBFGS stp", f)
                self.print_step_grads("LBFGS GRADS", gnorm)
            self.iterations += 1
            self.f_vals.append(f)

            if self.calc_curvatures and not self.use_curvatures:
                if self.num_positive_curvatures == self.use_curvatures_threshold:
                    raise BreakToUseCurvatures

        return self._f, self._g

    def _mpi_reduce_broadcast(self, var):
        var = comm.reduce(var, MPI.SUM, root=0)
        var = comm.bcast(var, root=0)
        return var

    def _poisson_target(self):
        fterm = (self.model_Lambda - self.Imeas * self.log_Lambda).sum()
        return fterm

    def _poisson_d(self, d):
        gterm = (d*self.one_minus_k_over_Lambda).sum()
        return gterm

    def _poisson_d2(self, d, d2):
        cterm = d2*self.one_minus_k_over_Lambda + d*d*self.k_over_squared_Lambda
        return cterm.sum()

    def _gaussian_target(self):
        fterm = .5 * (self.log2pi + self.log_Lambda_plus_sigma_readout + self.u * self.u * self.one_over_v).sum()
        return fterm

    def _gaussian_d(self, d):
        gterm = .5*(d*self.one_over_v*self.one_minus_2u_minus_u_squared_over_v).sum()
        return gterm

    def _gaussian_d2(self, d, d2):
        cterm = self.one_over_v*(d2*self.one_minus_2u_minus_u_squared_over_v -
                                 d*d*(self.one_over_v*self.one_minus_2u_minus_u_squared_over_v -
                                      (2+2*self.u*self.one_over_v + self.u*self.u*self.one_over_v*self.one_over_v)))
        cterm = .5*cterm.sum()
        return cterm

    def _evaluate_log_averageI(self):  # for Poisson only stats
        # fix log(x<=0)
        try:
            self.log_Lambda = np_log(self.model_Lambda)
        except FloatingPointError:
            pass
        if any((self.model_Lambda <= 0).ravel()):
            self.num_kludge += 1
            #print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL (kludges=%d)!!!!!!!!!\n<><><><><><><><><>\n" % self.num_kludge)
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda[self.model_Lambda <= 0] = 0  # FIXME with Gaussian model

    def _evaluate_log_averageI_plus_sigma_readout(self):
        L = self.model_Lambda + self.sigma_r**2
        L_is_neg = (L <= 0).ravel()
        if any(L_is_neg):
            self.num_kludge += 1
            print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda_plus_sigma_readout = np_log(L)
        self.log_Lambda_plus_sigma_readout[L <= 0] = 0  # will I even ever kludge ?

    def print_step(self, message, target):
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))
        rotX = self.x[self.rotX_xpos[self._i_shot]]
        rotY = self.x[self.rotY_xpos[self._i_shot]]
        rotZ = self.x[self.rotZ_xpos[self._i_shot]]
        rot_labels = ["rotX=%+3.7g" % rotX, "rotY=%+3.7g" % rotY, "rotZ=%+3.4g" % rotZ]

        rotx, roty, rotz, a_vals, c_vals, ncells_vals, scale_vals = self._unpack_internal(self.x)

        master_data = {"a": a_vals, "c": c_vals,
                       "Ncells": ncells_vals,
                       "scale": scale_vals,
                       "rotx": rotx,
                       "roty": roty,
                       "rotz": rotz}

        master_data = pandas.DataFrame(master_data)
        master_data["gain"] = self.x[self.gain_xpos]
        print("\n")
        print(master_data.to_string(float_format="%2.6g"))

    def print_step_grads(self, message, target):
        names = self.UCELL_MAN[self._i_shot].variable_names
        vals = self.UCELL_MAN[self._i_shot].variables
        ucell_labels = []
        for i, (n, v) in enumerate(zip(names, vals)):
            grad = self._g[self.ucell_xstart[self._i_shot] + i]
            ucell_labels.append('G%s=%+2.7g' % (n, grad))
        rotX = self._g[self.rotX_xpos[self._i_shot]]
        rotY = self._g[self.rotY_xpos[self._i_shot]]
        rotZ = self._g[self.rotZ_xpos[self._i_shot]]
        rot_labels = ["GrotX=%+3.7g" % rotX, "GrotY=%+3.7g" % rotY, "GrotZ=%+3.4g" % rotZ]
        xnorm = norm(self.x)

        rotx, roty, rotz, a_vals, c_vals, ncells_vals, scale_vals = self._unpack_internal(self._g)

        master_data = {"Ga": a_vals, "Gc": c_vals,
                       "GNcells": ncells_vals,
                       "Gscale": scale_vals,
                       "Grotx": rotx,
                       "Groty": roty,
                       "Grotz": rotz}
        master_data = pandas.DataFrame(master_data)
        master_data["Ggain"] = self._g[self.gain_xpos]
        print(master_data.to_string(float_format="%.3g"))

        if self.calc_curvatures:
            rotx, roty, rotz, a_vals, c_vals, ncells_vals, scale_vals = self._unpack_internal(self.curv)

            master_data = {"CUa": a_vals, "CUc": c_vals,
                           "CUNcells": ncells_vals,
                           "CUscale": scale_vals,
                           "CUrotx": rotx,
                           "CUroty": roty,
                           "CUrotz": rotz}
            master_data = pandas.DataFrame(master_data)
            master_data["CUgain"] = self.curv[self.gain_xpos]
            print(master_data.to_string(float_format="%.3g"))
        Xnorm = norm(self.x)
        print("\n\t|G|=%2.7g, eps*|X|=%2.7g" % (target, Xnorm*self.trad_conv_eps))
        print("\n")

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
            anglesXYZ = self.x[self.rotX_xpos[i_shot]], self.x[self.rotY_xpos[i_shot]], self.x[self.rotZ_xpos[i_shot]]
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
