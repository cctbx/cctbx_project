
from copy import deepcopy

import os
from libtbx.mpi4py import MPI
#from mpi4py import MPI
COMM = MPI.COMM_WORLD

from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner
from simtbx.diffBragg.refine_launcher import LocalRefinerLauncher, ShotData
from simtbx.diffBragg import utils


def global_refiner_from_parameters(refl_tbl, expt_list, params):
    launcher = GlobalRefinerLauncher(params)
    return launcher.launch_refiner(refl_tbl, expt_list)


class GlobalRefinerLauncher(LocalRefinerLauncher):

    def __init__(self, params):
        super().__init__(params)
        self.n_shots_on_rank = None

    @property
    def num_shots_on_rank(self):
        return len(self.shot_rois)

    def _alias_refiner(self):
        self._Refiner = GlobalRefiner

    def launch_refiner(self, refl_tbl, expt_list, miller_data=None):
        self._alias_refiner()

        num_exp = len(expt_list)
        assert num_exp == len(set(refl_tbl["id"]))

        detector = expt_list[0].detector  # TODO verify all shots have the same detector
        self._init_panel_group_information(detector)

        shot_idx = 0  # each rank keeps index of the shots local to it
        rank_panel_groups_refined = set()
        rank_local_parameters = []
        for i_exp in range(num_exp):
            if i_exp % COMM.size != COMM.rank:
                continue
            expt = expt_list[i_exp]
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
            self.shot_spectra[shot_idx] = self.SIM.beam.spectrum
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
        self.RUC.eta_init ={}
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
