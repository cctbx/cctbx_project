from __future__ import absolute_import, division, print_function
from simtbx.diffBragg.stage_two_utils import PAR_from_params
import os
import sys
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
from dials.array_family import flex
import numpy as np
try:
    import pandas
except ImportError:
    print("Pandas is required. Install using 'libtbx.python -m pip install pandas'")
    exit()
from xfel.merging.application.utils.memory_usage import get_memory_usage
from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner
from simtbx.diffBragg.refine_launcher import LocalRefinerLauncher
from simtbx.diffBragg import utils
from simtbx.diffBragg import hopper_utils
from dxtbx.model.experiment_list import ExperimentListFactory
from simtbx.diffBragg.prep_stage2_input import prep_dataframe
import logging

LOGGER = logging.getLogger("main")


def global_refiner_from_parameters(params):
    launcher = GlobalRefinerLauncher(params)
    # TODO read on each rank, or read and broadcast ?
    LOGGER.info("EVENT: read input pickle")
    pandas_table = pandas.read_pickle(params.pandas_table)
    LOGGER.info("EVENT: BEGIN prep dataframe")
    if params.prep_time > 0:
        pandas_table = prep_dataframe(pandas_table, params.prep_time)
    LOGGER.info("EVENT: DONE prep dataframe")
    return launcher.launch_refiner(pandas_table)


class GlobalRefinerLauncher(LocalRefinerLauncher):

    def __init__(self, params):
        super().__init__(params)
        self.n_shots_on_rank = None
        self.df = None
        self.WATCH_MISORIENTAION = False   # TODO add a phil

        self.SCALE_INIT_PER_SHOT = {}
        self.NCELLS_INIT_PER_SHOT = {}

    @property
    def num_shots_on_rank(self):
        return len(self.Modelers)

    def _alias_refiner(self):
        self._Refiner = GlobalRefiner

    def launch_refiner(self, pandas_table, miller_data=None):
        self._alias_refiner()

        if COMM.rank == 0:
            self.create_cache_dir()
        COMM.Barrier()

        num_exp = len(pandas_table)
        #first_exper_file = pandas_table.opt_exp_name.values[0]
        first_exper_file = pandas_table.exp_name.values[0]
        detector = ExperimentListFactory.from_json_file(first_exper_file, check_format=False)[0].detector
        # TODO verify all shots have the same detector ?
        if self.params.refiner.reference_geom is not None:
            detector = ExperimentListFactory.from_json_file(self.params.refiner.reference_geom, check_format=False)[0].detector
            print("Using reference geom from expt %s" % self.params.refiner.reference_geom)

        if COMM.size > num_exp:
            raise ValueError("Requested %d MPI ranks to process %d shots. Reduce number of ranks to %d"
                             % (COMM.size, num_exp, num_exp))
        self._init_panel_group_information(detector)

        self.verbose = False
        if COMM.rank == 0:
            self.verbose = self.params.refiner.verbose > 0
        shot_idx = 0  # each rank keeps index of the shots local to it
        rank_panel_groups_refined = set()
        rank_local_parameters = []
        exper_names = pandas_table.exp_name
        assert len(exper_names) == len(set(exper_names))
        # TODO assert all exper are single-file, probably way before this point
        LOGGER.info("EVENT: begin loading inputs")
        for i_exp, exper_name in enumerate(exper_names):
            if i_exp % COMM.size != COMM.rank:
                continue
            LOGGER.info("EVENT: BEGIN loading experiment list")
            expt_list = ExperimentListFactory.from_json_file(exper_name, check_format=True)
            LOGGER.info("EVENT: DONE loading experiment list")
            if len(expt_list) != 1:
                print("Input experiments need to have length 1, %s does not" % exper_name)
            expt = expt_list[0]
            expt.detector = detector  # in case of supplied ref geom
            self._check_experiment_integrity(expt)

            #exper_dataframe = pandas_table.query("opt_exp_name=='%s'" % exper_name)
            exper_dataframe = pandas_table.query("exp_name=='%s'" % exper_name)
            rotX, rotY, rotZ = exper_dataframe[["rotX", "rotY", "rotZ"]].values[0]
            self.rotXYZ_inits[shot_idx] = rotX, rotY, rotZ

            self._set_initial_model_for_shot(shot_idx, exper_dataframe)

            refl_name = exper_dataframe.predictions.values[0]
            refls = flex.reflection_table.from_file(refl_name)
            # FIXME need to remove (0,0,0) bboxes
            good_sel = flex.bool([h != (0, 0, 0) for h in list(refls["miller_index"])])
            refls = refls.select(good_sel)

            #UcellMan = utils.manager_from_crystal(expt.crystal)
            opt_uc_param = exper_dataframe[["a","b","c","al","be","ga"]].values[0]
            UcellMan = utils.manager_from_params(opt_uc_param)

            if self.symbol is None:
                if self.params.refiner.force_symbol is not None:
                    self.symbol = self.params.refiner.force_symbol
                else:
                    self.symbol = expt.crystal.get_space_group().type().lookup_symbol()
            else:
                if self.params.refiner.force_symbol is None:
                    if expt.crystal.get_space_group().type().lookup_symbol() != self.symbol:
                        raise ValueError("Crystals should all have the same space group symmetry")

            #if "spectrum_filename" in list(exper_dataframe):
            #    self.params.simulator.spectrum.filename = exper_dataframe.spectrum_filename.values[0] #self.input_spectrumnames[self.i_exp]

            if shot_idx == 0:  # each rank initializes a simulator only once
                if self.params.simulator.init_scale != 1:
                    print("WARNING: For stage_two , it is assumed that total scale is stored in the pandas dataframe")
                    print("WARNING: resetting params.simulator.init_scale to 1!")
                    self.params.simulator.init_scale = 1
                self._init_simulator(expt, miller_data)
                if self.params.profile:
                    self.SIM.record_timings = True
                if self.params.refiner.stage_two.Fref_mtzname is not None:
                    self.Fref = utils.open_mtz(self.params.refiner.stage_two.Fref_mtzname,
                                               self.params.refiner.stage_two.Fref_mtzcol)

            #Hi = list(refls["miller_index"])
            #self.Hi[shot_idx] = Hi
            #self.Hi_asu[shot_idx] = map_hkl_list(Hi, True, self.symbol)

            LOGGER.info("EVENT: LOADING ROI DATA")
            shot_modeler = hopper_utils.DataModeler(self.params)
            if not shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=self.symbol):
                raise("Failed to gather data from experiment %s", exper_name)

            self.Hi[shot_idx] = shot_modeler.Hi
            self.Hi_asu[shot_idx] = shot_modeler.Hi_asu

            LOGGER.info("EVENT: DONE LOADING ROI")
            shot_modeler.ucell_man = UcellMan
            #if shot_data.roi_darkRMS is not None:
            #    self.shot_roi_darkRMS[shot_idx] = shot_data.roi_darkRMS
            #if "rlp" in refls:
            #    self.shot_reso[shot_idx] = 1/np.linalg.norm(refls["rlp"], axis=1)

            if self.params.spectrum_from_imageset:
                shot_spectra = hopper_utils.downsamp_spec(self.SIM, self.params, expt, return_and_dont_set=True)

            elif "spectrum_filename" in list(exper_dataframe) and exper_dataframe.spectrum_filename.values[0] is not None:
                shot_spectra = utils.load_spectra_from_dataframe(exper_dataframe)

            else:
                total_flux = exper_dataframe.total_flux.values[0]
                if total_flux is None:
                    total_flux = self.params.simulator.total_flux
                shot_spectra = [(expt.beam.get_wavelength(), total_flux)]

            shot_modeler.spectra = shot_spectra

            if "detz_shift_mm" in list(exper_dataframe):
                shot_modeler.originZ_init = exper_dataframe.detz_shift_mm.values[0]*1e-3
            else:
                shot_modeler.originZ_init = 0
            shot_modeler.exper_name = exper_name

            shot_panel_groups_refined = self.determine_refined_panel_groups(shot_modeler.pids)
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

            n_ncells_def_params = 3
            nrot_params = 3
            n_unitcell_params = len(UcellMan.variables)  # TODO verify all crystals have same space group sym
            n_spotscale_params = 1
            n_originZ_params = 1
            n_eta_params = 3
            n_tilt_params = 3 * len(shot_modeler.rois)
            n_sausage_params = 4*self.params.simulator.crystal.num_sausages
            n_per_spot_scales = len(shot_modeler.rois)
            n_local_unknowns = nrot_params + n_unitcell_params + n_ncells_param + n_ncells_def_params + n_spotscale_params + n_originZ_params \
                               + n_tilt_params + n_eta_params + n_sausage_params + n_per_spot_scales

            rank_local_parameters.append(n_local_unknowns)
            shot_idx += 1
            if COMM.rank == 0:
                self._mem_usage()
                print("Finished loading image %d / %d" % (i_exp+1, len(exper_names)), flush=True)


            shot_modeler.PAR = PAR_from_params(self.params, expt, best=exper_dataframe)
            self.Modelers[i_exp] = shot_modeler

        LOGGER.info("DONE LOADING DATA; ENTER BARRIER")
        COMM.Barrier()
        LOGGER.info("DONE LOADING DATA; EXIT BARRIER")
        #if not self.shot_roi_darkRMS:
        self.shot_roi_darkRMS = None

        # TODO warn that per_spot_scale refinement not intended for ensemble mode

        LOGGER.info("parameter stuff")
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
        LOGGER.info("done with parameter stuff")

        LOGGER.info("EVENT: Gathering global HKL information")
        self._gather_Hi_information()
        LOGGER.info("EVENT: FINISHED gather global HKL information")
        if self.params.roi.cache_dir_only:
            print("Done creating cache directory and cache_dir_only=True, so goodbye.")
            sys.exit()

        n_spectra_params = 2 # if self.params.refiner.refine_spectra is not None else 0
        n_panelRot_params = 3*self.n_panel_groups
        n_panelXYZ_params = 3*self.n_panel_groups
        n_global_params = n_spectra_params + n_panelRot_params + n_panelXYZ_params

        #if self.params.refiner.refine_Fcell is not None and any(self.params.refiner.refine_Fcell):
        n_global_params += self.num_hkl_global

        #TODO why is this init_refiner call here ?
        self._init_refiner(n_local_unknowns=total_local_param_on_rank,
                           n_global_unknowns=n_global_params,
                           local_idx_start=local_param_offset_per_rank,
                           global_idx_start=total_local_unknowns_all_ranks)

        self.n_ncells_param = n_ncells_param
        self.n_spectra_params = n_spectra_params

        # in case of GPU
        LOGGER.info("BEGIN DETERMINE MAX PIX")
        self.NPIX_TO_ALLOC = self._determine_per_rank_max_num_pix()
        # TODO in case of randomize devices, shouldnt this be total max across all ranks?
        n = COMM.gather(self.NPIX_TO_ALLOC)
        if COMM.rank == 0:
            n = max(n)
        self.NPIX_TO_ALLOC = COMM.bcast(n)
        LOGGER.info("DONE DETERMINE MAX PIX")

        if not self.params.refiner.randomize_devices:
            self.DEVICE_ID = COMM.rank % self.params.refiner.num_devices

        self._mem_usage()

        LOGGER.info("EVENT: launch refiner")
        self._launch(total_local_param_on_rank, n_global_params,
                     local_idx_start=local_param_offset_per_rank,
                     global_idx_start=total_local_unknowns_all_ranks)

        return self.RUC

    def _mem_usage(self):
        memMB = get_memory_usage()
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

        self.sausages_init = {}
        self.RUC.eta_init = {}
        self.RUC.ucell_mins = {}
        self.RUC.ucell_maxs = {}

        ncells_def_init = self.params.simulator.crystal.ncells_def
        self.RUC.ncells_def_init = {}
        for i_shot in range(self.num_shots_on_rank):
            self.RUC.ncells_def_init[i_shot] = ncells_def_init # TODO fix
            #TODO sausages from stage 1 ?
            self.sausages_init[i_shot] = [0, 0, 0, 1]*self.params.simulator.crystal.num_sausages
            for i_saus in range(self.params.simulator.crystal.num_sausages):
                self.sausages_init[i_shot][4*i_saus-1] = i_saus + 0.1
            #TODO generalize and allow loading of mosaicity crystal model for anisotropic mosaicity
            if self.params.simulator.crystal.anisotropic_mosaicity is not None:
                raise NotImplemented("Stage 2 doesnt support aniso mosaicity")
            self.RUC.eta_init[i_shot] = [self.params.simulator.crystal.mosaicity, 0, 0]

            if ucell_maxs and ucell_mins:
                self.RUC.ucell_mins[i_shot] = ucell_mins
                self.RUC.ucell_maxs[i_shot] = ucell_maxs

        self.RUC.m_init = self.NCELLS_INIT_PER_SHOT  # enfore these are not None?
        self.RUC.spot_scale_init = self.SCALE_INIT_PER_SHOT
        # TODO per shot ncells def and per shot eta
        #TODO sausages
        #self.sausages_init = {}
        #for i_shot in self.SCALE_INIT_PER_SHOT:
        #    scale_fac = self.SCALE_INIT_PER_SHOT[i_shot]
        #    nsaus = self.params.simulator.crystal.num_sausages
        #    self.sausages_init[i_shot] = []
        #    for i_saus in range(nsaus):
        #        self.sausages_init[i_shot] += [0,0,0,scale_fac / nsaus]

    def _prep_blue_sausages(self):
        if self.params.refiner.refine_blueSausages:
            from scitbx.array_family import flex
            # TODO params.refiner.init.sausages is not supported for ensemble refinement
            init = [0, 0, 0, 1]*self.params.simulator.crystal.num_sausages
            self.SIM.D.update_number_of_sausages(self.params.simulator.crystal.num_sausages)
            x = flex.double(init[0::4])
            y = flex.double(init[1::4])
            z = flex.double(init[2::4])
            scale = flex.double(init[3::4])
            self.SIM.D.set_sausages(x, y, z, scale)

    def determine_refined_panel_groups(self, pids):
        refined_groups = []
        #assert len(pids) == len(selection_flags)
        for i, pid in enumerate(pids):
            if self.panel_group_from_id[pid] not in refined_groups:
                refined_groups.append(self.panel_group_from_id[pid])
        return refined_groups

    def _determine_per_rank_max_num_pix(self):
        max_npix = 0
        for i_shot in self.Modelers:
            modeler = self.Modelers[i_shot]
            x1, x2, y1, y2 = map(np.array, zip(*modeler.rois))
            npix = np.sum((x2-x1)*(y2-y1))
            max_npix = max(npix, max_npix)
            print("Rank %d, shot %d has %d pixels" % (COMM.rank, i_shot+1, npix))
        print("Rank %d, max pix to be modeled: %d" % (COMM.rank, max_npix))
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
        self.Hi_all_ranks = COMM.reduce(self.Hi_all_ranks)
        self.Hi_all_ranks = COMM.bcast(self.Hi_all_ranks)

        self.Hi_asu_all_ranks = COMM.reduce(self.Hi_asu_all_ranks)
        self.Hi_asu_all_ranks = COMM.bcast(self.Hi_asu_all_ranks)

        marr_unique_h = self._get_unique_Hi()

        # this will map the measured miller indices to their index in the LBFGS parameter array self.x
        self.idx_from_asu = {h: i for i, h in enumerate(set(self.Hi_asu_all_ranks))}
        # we will need the inverse map during refinement to update the miller array in diffBragg, so we cache it here
        self.asu_from_idx = {i: h for i, h in enumerate(set(self.Hi_asu_all_ranks))}

        self.num_hkl_global = len(self.idx_from_asu)

        fres = marr_unique_h.d_spacings()
        self.res_from_asu = {h: res for h, res in zip(fres.indices(), fres.data())}

    def _get_unique_Hi(self):
        COMM.barrier()
        if COMM.rank == 0:
            from cctbx.crystal import symmetry
            from cctbx import miller
            from cctbx.array_family import flex as cctbx_flex

            ii = list(self.Modelers.keys())[0]
            uc = self.Modelers[ii].ucell_man
            params = uc.a, uc.b, uc.c, uc.al * 180 / np.pi, uc.be * 180 / np.pi, uc.ga * 180 / np.pi
            if self.params.refiner.force_unit_cell is not None:
                params = self.params.refiner.force_unit_cell
            symm = symmetry(unit_cell=params, space_group_symbol=self.symbol)
            hi_asu_flex = cctbx_flex.miller_index(self.Hi_asu_all_ranks)
            mset = miller.set(symm, hi_asu_flex, anomalous_flag=True)
            marr = miller.array(mset)
            binner = marr.setup_binner(d_max=self.params.refiner.stage_two.d_max, d_min=self.params.refiner.stage_two.d_min,
                                       n_bins=self.params.refiner.stage_two.n_bin)
            from collections import Counter
            print("Average multiplicities:")
            print("<><><><><><><><><><><><>")
            for i_bin in range(self.params.refiner.stage_two.n_bin - 1):
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
            self.binner = mset.setup_binner(d_min=self.params.refiner.stage_two.d_min,
                                            d_max=self.params.refiner.stage_two.d_max,
                                            n_bins=self.params.refiner.stage_two.n_bin)
            mset.completeness(use_binning=True).show()
            marr_unique_h = miller.array(mset)
            print("Rank %d: total miller vars=%d" % (COMM.rank, len(set(self.Hi_asu_all_ranks))))
        else:
            marr_unique_h = None

        marr_unique_h = COMM.bcast(marr_unique_h)
        return marr_unique_h

    def _set_initial_model_for_shot(self, shot_idx, dataframe):
        if shot_idx in self.SCALE_INIT_PER_SHOT or shot_idx in self.NCELLS_INIT_PER_SHOT:
            raise KeyError("Already set initial model for shot %d on rank %d" % (shot_idx, COMM.rank))

        self.SCALE_INIT_PER_SHOT[shot_idx] = np.sqrt(dataframe.spot_scales.values[0])
        # TODO should this be sqrt?

        ncells_init = dataframe.ncells.values[0]  # either a 1-tuple or a 3-tuple
        if len(ncells_init) == 1:
            ncells_init += (ncells_init[0], ncells_init[0])
        self.NCELLS_INIT_PER_SHOT[shot_idx] = ncells_init

    def get_parameter_hdf5_path(self):
        hdf5_path = None
        if self.params.refiner.dump_params:
            if self.params.refiner.parameter_hdf5_path is not None:
                param_folder = self.params.refiner.parameter_hdf5_path
                if not os.path.exists(param_folder):
                    os.makedirs(param_folder)
                hdf5_path = os.path.join(param_folder, "params_rank%d.h5" % COMM.rank)
            elif self.params.refiner.io.output_dir is not None:
                param_dir = os.path.join(self.params.refiner.io.output_dir, "params")
                if COMM.rank == 0:
                    if not os.path.exists(param_dir):
                        os.makedirs(param_dir)
                hdf5_path = os.path.join(param_dir, "params_rank%d.h5" %COMM.rank)
        return hdf5_path
