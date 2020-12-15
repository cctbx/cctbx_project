
from copy import deepcopy
import numpy as np
from simtbx.command_line.stage_one import Script as stage_one_Script
import pandas

import os
from libtbx.mpi4py import MPI
#from mpi4py import MPI
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

    def launch_refiner(self, refl_tbl, expt_list, miller_data=None, pandas_list=None):
        self._alias_refiner()

        file_path_input = False
        num_exp = len(set(refl_tbl["id"]))
        if not isinstance(expt_list, str):
            assert num_exp == len(expt_list)
            detector = expt_list[0].detector  # TODO verify all shots have the same detector
        else:
            detector = ExperimentListFactory.from_json_file(expt_list, check_format=False)[0].detector
            file_path_input = True
        if self.params.refiner.reference_geom is not None:
            detector = ExperimentListFactory.from_json_file(self.params.refiner.reference_geom, check_format=False)[0].detector
            print("USING REFERENCE GEOMZ!!!")

        if COMM.size > num_exp:
            raise ValueError("Requested %d MPI ranks to process %d shots. Reduce number of ranks to %d"
                             % (COMM.size, num_exp, num_exp))
        self._init_panel_group_information(detector)

        spectra_list = self._try_loading_spectrum_filelist()
        if spectra_list is not None:
            assert len(spectra_list) == num_exp

        shot_idx = 0  # each rank keeps index of the shots local to it
        rank_panel_groups_refined = set()
        rank_local_parameters = []
        for i_exp in range(num_exp):
            if i_exp % COMM.size != COMM.rank:
                continue
            if file_path_input:
                expt = stage_one_Script.exper_json_single_file(expt_list, i_exp=i_exp)
            else:
                expt = expt_list[i_exp]
            expt.detector = detector  # in case of supplied ref geom
            refls = refl_tbl.select(refl_tbl['id'] == i_exp)
            self._check_experiment_integrity(expt)

            shot_data = self.load_roi_data(refls, expt)
            if shot_data is None:
                raise ValueError("Cannot refine!")

            UcellMan = utils.manager_from_crystal(expt.crystal)

            if shot_idx == 0:  # each rank initializes a simulator only once
                self._init_simulator(expt, miller_data)

            self.shot_ucell_managers[shot_idx] = UcellMan
            self.shot_rois[shot_idx] = shot_data.rois
            self.shot_nanoBragg_rois[shot_idx] = shot_data.nanoBragg_rois
            self.shot_roi_imgs[shot_idx] = shot_data.roi_imgs
            # TODO spectra ensemble
            #self.shot_spectra[shot_idx] = self.SIM.beam.spectrum
            if spectra_list is not None:
                self.shot_spectra[shot_idx] = utils.load_spectra_file(
                    spectra_list[i_exp],
                    self.params.simulator.total_flux,
                    self.params.simulator.spectrum.stride,
                    as_spectrum=True)
                print("Loaded spectra!!")
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

        assert self.shot_rois, "cannot refine without shots! Probably Ranks than shots!"

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

        n_spectra_params = 2 if self.params.refiner.refine_spectra is not None else 0
        n_panelRot_params = 3*self.n_panel_groups
        n_panelXYZ_params = 3*self.n_panel_groups
        n_global_params = n_spectra_params + n_panelRot_params + n_panelXYZ_params

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

        self._launch(total_local_param_on_rank, n_global_params,
                     local_idx_start=local_param_offset_per_rank,
                     global_idx_start=total_local_unknowns_all_ranks)

        return self.RUC

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

