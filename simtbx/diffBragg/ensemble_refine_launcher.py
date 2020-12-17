
from copy import deepcopy
from dials.array_family import flex
import numpy as np
import pandas
from xfel.merging.application.utils.memory_usage import get_memory_usage

from simtbx.diffBragg.utils import map_hkl_list

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner
from simtbx.diffBragg.refine_launcher import LocalRefinerLauncher, ShotData
from simtbx.diffBragg import utils
from dxtbx.model.experiment_list import ExperimentListFactory


def global_refiner_from_parameters(refl_tbl, expt_list, params, pandas_list=None):
    launcher = GlobalRefinerLauncher(params, pandas_list)
    return launcher.launch_refiner(refl_tbl, expt_list, pandas_list)


class GlobalRefinerLauncher(LocalRefinerLauncher):

    def __init__(self, params, pandas_list=None):
        super().__init__(params)
        self.n_shots_on_rank = None
        self.df = None
        if pandas_list is not None:
            self.df = pandas.read_pickle(pandas_list)
        self.WATCH_MISORIENTAION = False   # TODO add a phil

    @property
    def num_shots_on_rank(self):
        return len(self.shot_rois)

    def _alias_refiner(self):
        self._Refiner = GlobalRefiner

    def launch_refiner(self, pandas_table, miller_data=None):
        self._alias_refiner()
        if not "opt_exp_name" in list(pandas_table) or not "predictions" in list(pandas_table):
            raise KeyError("Pandas table needs opt_exp_name and predicictions columns")
        num_exp = len(pandas_table)
        first_exper_file = pandas_table.opt_exp_name.values[0]
        detector = ExperimentListFactory.from_json_file(first_exper_file, check_format=False)[0].detector
        # TODO verify all shots have the same detector ?
        if self.params.refiner.reference_geom is not None:
            detector = ExperimentListFactory.from_json_file(self.params.refiner.reference_geom, check_format=False)[0].detector
            print("Using reference geom from expt %s" % self.params.refiner.reference_geom)

        if COMM.size > num_exp:
            raise ValueError("Requested %d MPI ranks to process %d shots. Reduce number of ranks to %d"
                             % (COMM.size, num_exp, num_exp))
        self._init_panel_group_information(detector)

        shot_idx = 0  # each rank keeps index of the shots local to it
        rank_panel_groups_refined = set()
        rank_local_parameters = []
        exper_names = pandas_table.opt_exp_name.unique()
        # TODO assert all exper are single-file, probably way before this point
        for i_exp, exper_name in enumerate(exper_names):
            if i_exp % COMM.size != COMM.rank:
                continue
            expt = ExperimentListFactory.from_json_file(exper_name, check_format=True)
            expt.detector = detector  # in case of supplied ref geom
            self._check_experiment_integrity(expt)

            exper_dataframe = pandas_table.query("opt_exp_name=='%s'" % exper_name)
            refl_name = exper_dataframe.predictions.values[0]
            refls = flex.reflection_table.from_file(refl_name)

            Hi = list(refls["miller_index"])
            self.Hi[shot_idx] = Hi
            if self.symbol is None:
                self.symbol = expt.crystal.get_space_group().type().lookup_symbol()
            else:
                if self.symbol != expt.crystal.get_space_group().type().lookup_symbol():
                    raise ValueError("Crystals should all have the same space group symmetry")
            self.Hi_asu[shot_idx] = map_hkl_list(self.Hi, True, self.symbol)

            shot_data = self.load_roi_data(refls, expt)
            if shot_data is None:
                raise ValueError("Cannot refine!")

            UcellMan = utils.manager_from_crystal(expt.crystal)

            if shot_idx == 0:  # each rank initializes a simulator only once
                self._init_simulator(expt, miller_data)
                if self.params.refiner.stage_two.Fref_mtzname is not None:
                    self.Fref = utils.open_mtz(self.params.refiner.stage_two.Fref_mtzname,
                                               self.params.refiner.stage_two.Fref_mtzcol)

            self.shot_ucell_managers[shot_idx] = UcellMan
            self.shot_rois[shot_idx] = shot_data.rois
            self.shot_nanoBragg_rois[shot_idx] = shot_data.nanoBragg_rois
            self.shot_roi_imgs[shot_idx] = shot_data.roi_imgs
            if "spectrum_filename" in list(exper_dataframe):
                self.shot_spectra[shot_idx] = utils.load_spectra_from_dataframe(exper_dataframe)
            else:
                self.shot_spectra[shot_idx] = [(expt.beam.get_wavelength(), self.params.simulator.total_flux)]
            self.shot_crystal_models[shot_idx] = expt.crystal
            self.shot_crystal_model_refs[shot_idx] = deepcopy(expt.crystal)
            self.shot_xrel[shot_idx] = shot_data.xrel
            self.shot_yrel[shot_idx] = shot_data.yrel
            self.shot_abc_inits[shot_idx] = shot_data.tilt_abc
            self.shot_panel_ids[shot_idx] = shot_data.pids
            self.shot_originZ_init[shot_idx] = 0
            self.shot_selection_flags[shot_idx] = shot_data.selection_flags
            self.shot_background[shot_idx] = shot_data.background
            shot_idx += 1

            shot_panel_groups_refined = self.determine_refined_panel_groups(shot_data.pids, shot_data.selection_flags)
            rank_panel_groups_refined = rank_panel_groups_refined.union(set(shot_panel_groups_refined))

            # <><><><><><>><><><><><><><><><>
            # determine number of local parameters:
            # <><><><><><>><><><><><><><><><><>
            if not any(self.NCELLS_MASK):
                n_ncells_param = 3
            elif all(self.NCELLS_MASK):
                n_ncells_param = 1
            else:
                n_ncells_param = 2

            nrot_params = 3
            n_unitcell_params = len(UcellMan.variables)  # TODO verify all crystals have same space group sym
            n_spotscale_params = 1
            n_originZ_params = 1
            n_eta_params = 1
            n_tilt_params = 3 * len(shot_data.nanoBragg_rois)
            n_sausage_params = 4*self.params.simulator.crystal.num_sausages
            n_local_unknowns = nrot_params + n_unitcell_params + n_ncells_param + n_spotscale_params + n_originZ_params \
                               + n_tilt_params + n_eta_params + n_sausage_params

            rank_local_parameters.append(n_local_unknowns)

        if not self.shot_rois:
            raise RuntimeError("Cannot refine without shots! Probably Ranks than shots!")

        # TODO warn that per_spot_scale refinement not recommended in ensemble mode

        total_local_param_on_rank = sum(rank_local_parameters)
        local_per_rank = COMM.gather(total_local_param_on_rank)
        local_param_offset_per_rank = total_local_unknowns_all_ranks = None
        if COMM.rank == 0:
            local_param_offset_per_rank = {}
            xpos = 0
            for rank, n_unknown in enumerate(local_per_rank):
                local_param_offset_per_rank[rank] = xpos
                xpos += n_unknown
            total_local_unknowns_all_ranks = sum(local_per_rank)
        total_local_unknowns_all_ranks = COMM.bcast(total_local_unknowns_all_ranks)
        local_param_offset_per_rank = COMM.bcast(local_param_offset_per_rank)[COMM.rank]

        all_refined_groups = COMM.gather(rank_panel_groups_refined)
        panel_groups_refined = None
        if COMM.rank == 0:
            panel_groups_refined = set()
            for set_of_panels in all_refined_groups:
                panel_groups_refined = panel_groups_refined.union(set_of_panels)

        self.panel_groups_refined = list(COMM.bcast(panel_groups_refined))

        self._gather_Hi_information()

        n_spectra_params = 2 if self.params.refiner.refine_spectra is not None else 0
        n_panelRot_params = 3*self.n_panel_groups
        n_panelXYZ_params = 3*self.n_panel_groups
        n_global_params = n_spectra_params + n_panelRot_params + n_panelXYZ_params

        if self.params.refiner.refine_Fcell is not None and any(self.params.refiner.refine_Fcell):
            n_global_params += self.num_hkl_global

        self._init_refiner(n_local_unknowns=total_local_param_on_rank,
                           n_global_unknowns=n_global_params,
                           local_idx_start=local_param_offset_per_rank,
                           global_idx_start=total_local_unknowns_all_ranks)

        self.n_ncells_param = n_ncells_param
        self.n_spectra_params = n_spectra_params

        # in case of GPU
        self.NPIX_TO_ALLOC = self._determine_per_rank_max_num_pix()

        if not self.params.refiner.randomize_devices:
            self.DEVICE_ID = COMM.rank % self.params.refiner.num_devices

        self._mem_usage()

        self._launch(total_local_param_on_rank, n_global_params,
                     local_idx_start=local_param_offset_per_rank,
                     global_idx_start=total_local_unknowns_all_ranks)

        return self.RUC

    def _mem_usage(self):
        memMB = COMM.reduce(get_memory_usage())
        if COMM.rank == 0:
            import socket
            host = socket.gethostname()
            print("Rank 0 reporting memory usage: %f GB on Rank 0 node %s" % (memMB / 1e3, host))

    def _initialize_some_refinement_parameters(self):
        # TODO provide interface for taking inital conditions from all shots (in the event of restarting a simulation)
        # MOSAICBLOCK
        m_init = self.params.simulator.crystal.ncells_abc
        if self.n_ncells_param == 2:
            INVERTED_NCELLS_MASK = [int(not mask_val) for mask_val in self.NCELLS_MASK]
            self.RUC.ncells_mask = INVERTED_NCELLS_MASK
            m_init = [m_init[i_ncell] for i_ncell in sorted(set(INVERTED_NCELLS_MASK))]

        # UNITCELL
        ucell_maxs, ucell_mins = [], []
        if self.params.refiner.ranges.ucell_edge_percentage is not None:
            names = self.shot_ucell_managers[0].variable_names
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
                ucell_mins.append(valmin)
                ucell_maxs.append(valmax)
            self.RUC.use_ucell_ranges = True

        self.RUC.sausages_init = {}
        self.RUC.m_init = {}
        self.RUC.spot_scale_init = {}
        self.RUC.eta_init = {}
        self.RUC.ucell_inits = {}
        self.RUC.ucell_mins = {}
        self.RUC.ucell_maxs = {}
        for i_shot in range(self.num_shots_on_rank):
            self.RUC.sausages_init[i_shot] = [0, 0, 0, 1]
            self.RUC.m_init[i_shot] = m_init
            self.RUC.spot_scale_init[i_shot] = self.params.refiner.init.spot_scale
            self.RUC.eta_init[i_shot] = self.params.simulator.crystal.mosaicity
            self.RUC.ucell_inits[i_shot] = self.shot_ucell_managers[i_shot].variables

            if ucell_maxs and ucell_mins:
                self.RUC.ucell_mins[i_shot] = ucell_mins
                self.RUC.ucell_maxs[i_shot] = ucell_maxs

        #self.RUC.sausages_init = {i: [0, 0, 0, 1] for i in range(self.num_shots_on_rank)}
        #self.RUC.m_init = {i: m_init for i in range(self.num_shots_on_rank)}
        #self.RUC.spot_scale_init = {i: self.params.refiner.init.spot_scale for i in range(self.num_shots_on_rank)}
        #self.RUC.eta_init = {i: self.params.simulator.crystal.mosaicity for i in range(self.num_shots_on_rank)}
        #self.RUC.ucell_inits = {i: self.shot_ucell_managers[i].variables for i in range(self.num_shots_on_rank)}

    def _prep_blue_sausages(self):
        if self.params.refiner.refine_blueSausages:
            from scitbx.array_family import flex
            # TODO params.refiner.init.sausages is not supported for ensemble refinement
            #init = self.params.refiner.init.sausages
            #if init is None:
            init = [0, 0, 0, 1]*self.params.simulator.crystal.num_sausages
            self.SIM.D.update_number_of_sausages(self.params.simulator.crystal.num_sausages)
            x = flex.double(init[0::4])
            y = flex.double(init[1::4])
            z = flex.double(init[2::4])
            scale = flex.double(init[3::4])
            self.SIM.D.set_sausages(x, y, z, scale)

    def determine_refined_panel_groups(self, pids, selection_flags):
        refined_groups = []
        assert len(pids) == len(selection_flags)
        for i, pid in enumerate(pids):
            if selection_flags[i] and self.panel_group_from_id[pid] not in refined_groups:
                refined_groups.append(self.panel_group_from_id[pid])
        return refined_groups

    def _determine_per_rank_max_num_pix(self):
        max_npix = 0
        for i_shot in self.shot_rois:
            rois = self.shot_rois[i_shot]
            x1, x2, y1, y2 = map(np.array, zip(*rois))
            npix = np.sum((x2-x1)*(y2-y1))
            max_npix = max(npix, max_npix)
            print("Rank %d, shot %d has %d pixels" % (COMM.rank, i_shot+1, npix))
        print("Rank %d, max pix to be modeled: %d" % (COMM.rank, max_npix))
        #max_npix = COMM.gather(max_npix)
        #if COMM.rank == 0:
        #    max_npix = np.max(max_npix)
        #max_npix = COMM.bcast(max_npix)
        return max_npix

    def _try_loading_spectrum_filelist(self):
        file_list = None
        fpath = self.params.simulator.spectrum.filename_list
        if fpath is not None:
            file_list = [l.strip() for l in open(fpath, "r").readlines()]
            assert all([len(l.split()) == 1 for l in file_list]), "weird spectrum file %s"% fpath
        return file_list

    def _gather_Hi_information(self):
        nshots_on_this_rank = len(self.Hi)
        # aggregate all miller indices
        self.Hi_all_ranks, self.Hi_asu_all_ranks = [], []
        # TODO assert list types are stored in Hi and Hi_asu
        for i_shot in range(nshots_on_this_rank):
            self.Hi_all_ranks += self.Hi[i_shot]
            self.Hi_asu_all_ranks += self.Hi_asu[i_shot]
        self.Hi_all_ranks = COMM.reduce(self.Hi_all_ranks, root=0)
        self.Hi_all_ranks = COMM.bcast(self.Hi_all_ranks, root=0)
        self.Hi_asu_all_ranks = COMM.reduce(self.Hi_asu_all_ranks, root=0)
        self.Hi_asu_all_ranks = COMM.bcast(self.Hi_asu_all_ranks, root=0)

        marr_unique_h = self._get_unique_Hi()

        # this will map the measured miller indices to their index in the LBFGS parameter array self.x
        self.idx_from_asu = {h: i for i, h in enumerate(set(self.Hi_asu_all_ranks))}
        # we will need the inverse map during refinement to update the miller array in diffBragg, so we cache it here
        self.asu_from_idx = {i: h for i, h in enumerate(set(self.Hi_asu_all_ranks))}

        self.num_hkl_global = len(self.idx_from_asu)

        fres = marr_unique_h.d_spacings()
        self.res_from_asu = {h: res for h, res in zip(fres.indices(), fres.data())}

    def _get_unique_Hi(self):
        if COMM.rank == 0:
            from cctbx.crystal import symmetry
            from cctbx import miller
            from cctbx.array_family import flex as cctbx_flex

            uc = self.all_ucell_mans[0]
            params = uc.a, uc.b, uc.c, uc.al * 180 / np.pi, uc.be * 180 / np.pi, uc.ga * 180 / np.pi
            symm = symmetry(unit_cell=params, space_group_symbol=self.symbol)
            hi_asu_flex = cctbx_flex.miller_index(self.Hi_asu_all_ranks)
            mset = miller.set(symm, hi_asu_flex, anomalous_flag=True)
            marr = miller.array(mset)  # ,data=flex.double(len(hi_asu_felx),0))
            n_bin = 10
            binner = marr.setup_binner(d_max=999, d_min=2.125, n_bins=n_bin)
            from collections import Counter
            print("Average multiplicities:")
            print("<><><><><><><><><><><><>")
            for i_bin in range(n_bin - 1):
                dmax, dmin = binner.bin_d_range(i_bin + 1)
                F_in_bin = marr.resolution_filter(d_max=dmax, d_min=dmin)
                multi_in_bin = np.array(list(Counter(F_in_bin.indices()).values()))
                print("%2.5g-%2.5g : Multiplicity=%.4f" % (dmax, dmin, multi_in_bin.mean()))
                for ii in range(1, 100, 8):
                    print("\t %d refls with multi %d" % (sum(multi_in_bin == ii), ii))

            print("Overall completeness\n<><><><><><><><>")
            symm = symmetry(unit_cell=params, space_group_symbol=self.symbol)
            hi_flex_unique = cctbx_flex.miller_index(list(set(self.Hi_asu_all_ranks)))
            mset = miller.set(symm, hi_flex_unique, anomalous_flag=True)
            self.binner = mset.setup_binner(d_min=self.params.d_min, d_max=self.params.d_max, n_bins=self.params.n_bins)
            mset.completeness(use_binning=True).show()
            marr_unique_h = miller.array(mset)
            print("Rank %d: total miller vars=%d" % (COMM.rank, len(set(self.Hi_asu_all_ranks))))
        else:
            marr_unique_h = None

        marr_unique_h = COMM.bcast(marr_unique_h)
        return marr_unique_h