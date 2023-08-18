from __future__ import absolute_import, division, print_function
from simtbx.diffBragg.stage_two_utils import PAR_from_params
from collections import Counter
from itertools import chain
import os
import socket
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
from simtbx.diffBragg.refiners.stage_two_refiner import StageTwoRefiner
from simtbx.diffBragg import utils
from simtbx.diffBragg import hopper_utils
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from simtbx.diffBragg.prep_stage2_input import prep_dataframe
from cctbx import miller, crystal
import logging

LOGGER = logging.getLogger("diffBragg.main")


def global_refiner_from_parameters(params):
    launcher = RefineLauncher(params)
    # TODO read on each rank, or read and broadcast ?
    LOGGER.info("EVENT: read input pickle")
    pandas_table = pandas.read_pickle(params.pandas_table)
    if params.max_process > 0:
        pandas_table = pandas_table.iloc[:params.max_process]
    LOGGER.info("EVENT: BEGIN prep dataframe")
    if "exp_idx" not in list(pandas_table):
        pandas_table["exp_idx"] = 0
    work_distribution = prep_dataframe(pandas_table)
    LOGGER.info("EVENT: DONE prep dataframe")
    return launcher.launch_refiner(pandas_table, work_distribution=work_distribution)


class RefineLauncher:

    def __init__(self, params):
        self.params = self.check_parameter_integrity(params)
        self.n_shots_on_rank = None
        self.df = None
        self.Modelers = {}
        self.Hi = {}
        self.Hi_asu = {}
        self.symbol = None
        self.DEVICE_ID = 0


    @property
    def NPIX_TO_ALLOC(self):
        return self._NPIX_TO_ALLOC

    @NPIX_TO_ALLOC.setter
    def NPIX_TO_ALLOC(self, val):
        assert val> 0 or val == -1
        self._NPIX_TO_ALLOC = int(val)

    @staticmethod
    def check_parameter_integrity(params):
        if params.refiner.max_calls is None or len(params.refiner.max_calls) == 0:
            raise ValueError("Cannot refine because params.refiner.max_calls is empty")

        if os.environ.get("DIFFBRAGG_CUDA") is not None:
            params.refiner.use_cuda = True

        return params

    @property
    def num_shots_on_rank(self):
        return len(self.Modelers)

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

    def _init_simulator(self, expt, miller_data):
        self.SIM = utils.simulator_for_refinement(expt, self.params)
        # note self.SIM.D is a now diffBragg instance
        # update the miller data ?
        if miller_data is not None:
            self.SIM.crystal.miller_array = miller_data.as_amplitude_array()
            self.SIM.update_Fhkl_tuple()

    @staticmethod
    def _check_experiment_integrity(expt):
        for model in ["crystal", "detector", "beam", "imageset"]:
            if not hasattr(expt, model):
                raise ValueError("No %s in experiment, exiting. " % model)

    def launch_refiner(self, pandas_table, miller_data=None, work_distribution=None):
        self.load_inputs(pandas_table, miller_data=miller_data, work_distribution=work_distribution)
        LOGGER.info("EVENT: launch refiner")
        self._launch()
        return self.RUC

    def load_inputs(self, pandas_table, miller_data=None, work_distribution=None, refls_key='predictions'):
        """

        :param pandas_table: contains path to the experiments (pandas column exp_name) to be loaded
            the pandas table is expected to have been written by diffBragg.hopper or
            diffBragg.hopper_process . See method save_to_pandas in simtbx/command_line/hopper.py
            For example, if the outputdir of diffBragg.hopper was set to `all_shots`, then
            there should be a folder all_shots/pandas created which contains all of the per-shot pandas
            dataframes. They should be concatenated as follows, forming a suitable argument for this method
            >> import glob,pandas
            >> fnames = glob.glob("all_shots/pandas/rank*/*pkl")
            >> df = pandas.concat([ pandas.read_pickle(f) for f in fnames])
            >> df.reset_index(inplace=True, drop=True)
            >> df.to_pickle("all_shots.pkl")
            >> # Then later, as part of an MPI application, the following will load all data:
            >> RefineLauncher_instance.load_inputs(df, refls_key="stage1_refls")

        :param miller_data: Optional miller array for the structure factor component of the model
        :param refls_key: key specifying the reflection tables in the pandas table
            Modeled pixels will lie in shoeboxes centered on each x,y,z in xyzobs.px.value
        :return:
        """
        COMM.Barrier()
        num_exp = len(pandas_table)
        if "exp_idx" not in list(pandas_table):
            pandas_table["exp_idx"] = 0
        first_exper_file = pandas_table.exp_name.values[0]
        detector = ExperimentListFactory.from_json_file(first_exper_file, check_format=False)[0].detector
        if detector is None and self.params.refiner.reference_geom is None:
            raise RuntimeError("No detector in experiment, must provide a reference geom.")
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
            if self.params.refiner.gather_dir is not None and not os.path.exists(self.params.refiner.gather_dir):
              os.makedirs(self.params.refiner.gather_dir)
              LOGGER.info("MADE GATHER DIR %s" % self.params.refiner.gather_dir)
        COMM.barrier()
        shot_idx = 0  # each rank keeps index of the shots local to it
        rank_panel_groups_refined = set()
        exper_names = pandas_table.exp_name
        exper_ids = pandas_table.exp_idx.values
        shot_ids = list(zip(exper_names, exper_ids))
        assert len(shot_ids) == len(set(shot_ids))
        # TODO assert all exper are single-file, probably way before this point
        if work_distribution is None:
            worklist = range(COMM.rank, len(exper_names), COMM.size)
        else:
            worklist = work_distribution[COMM.rank]
        LOGGER.info("EVENT: begin loading inputs")

        # load the Fhkl model once here to check which hkl are missing (and filter from the refls below)
        first_exper = ExperimentList.from_file(exper_names[0], check_format=False)[0]
        Fhkl_model = utils.load_Fhkl_model_from_params_and_expt(self.params, first_exper)
        self.symbol = Fhkl_model.space_group().info().type().lookup_symbol()
        if self.params.refiner.force_symbol is not None:
            self.symbol = self.params.refiner.force_symbol
        LOGGER.info("Will use space group symbol %s" % self.symbol)
        Fhkl_model_p1 = Fhkl_model.expand_to_p1().generate_bijvoet_mates()
        Fhkl_model_p1_indices = set(Fhkl_model_p1.indices())

        for i_work, i_df in enumerate(worklist):
            exper_name = exper_names[i_df]
            exper_id = int(exper_ids[i_df])
            LOGGER.info("EVENT: BEGIN loading experiment list")
            # TODO: test that the diffBragg_benchmarks is not broken
            expt = hopper_utils.DataModeler.exper_json_single_file(exper_name, exper_id)
            expt_list = ExperimentList()
            expt_list.append(expt)
            LOGGER.info("EVENT: DONE loading experiment list")
            expt.detector = detector  # in case of supplied ref geom
            self._check_experiment_integrity(expt)

            exper_dataframe = pandas_table.query("exp_name=='%s'" % exper_name).query("exp_idx==%d" % exper_id)

            refl_name = exper_dataframe[refls_key].values[0]
            refls = flex.reflection_table.from_file(refl_name)
            refls = refls.select(refls['id'] == exper_id)

            try:
                miller_inds = list(refls["miller_index"])
                is_not_000 = [h != (0, 0, 0) for h in miller_inds]
                is_in_Fhkl_model = [h in Fhkl_model_p1_indices for h in miller_inds]
                LOGGER.debug("Only refining %d/%d refls whose HKL are in structure factor model" % (
                    np.sum(is_in_Fhkl_model), len(refls)))
                refl_sel = flex.bool(np.logical_and(is_not_000, is_in_Fhkl_model))
                refls = refls.select(refl_sel)
            except KeyError:
                pass

            opt_uc_param = exper_dataframe[["a","b","c","al","be","ga"]].values[0]
            UcellMan = utils.manager_from_params(opt_uc_param)

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

            LOGGER.info("EVENT: LOADING ROI DATA")
            shot_modeler = hopper_utils.DataModeler(self.params)
            shot_modeler.exper_name = exper_name
            shot_modeler.exper_idx = exper_id
            shot_modeler.refl_name = refl_name
            shot_modeler.rank = COMM.rank
            if self.params.refiner.load_data_from_refl:
                gathered = shot_modeler.GatherFromReflectionTable(expt, refls, sg_symbol=self.symbol)
            else:
                # Note: no need to pass exper_id here because expt and refls have already been sliced out
                gathered = shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=self.symbol)
            if not gathered:
                raise IOError("Failed to gather data from experiment %s", exper_name)
                COMM.abort()

            if self.params.refiner.gather_dir is not None:
                gathered_name = os.path.splitext(os.path.basename(exper_name))[0]
                gathered_name += "_withData.refl"
                gathered_name = os.path.join(self.params.refiner.gather_dir, gathered_name )
                shot_modeler.dump_gathered_to_refl(gathered_name, do_xyobs_sanity_check=False) #True)
                LOGGER.info("SAVED ROI DATA TO %s" % gathered_name)
                if self.params.refiner.test_gathered_file:
                    all_data = shot_modeler.all_data.copy()
                    all_roi_id = shot_modeler.roi_id.copy()
                    all_bg = shot_modeler.all_background.copy()
                    all_trusted = shot_modeler.all_trusted.copy()
                    all_pids = np.array(shot_modeler.pids)
                    all_rois = np.array(shot_modeler.rois)
                    new_Modeler = hopper_utils.DataModeler(self.params)
                    assert new_Modeler.GatherFromReflectionTable(exper_name, gathered_name, sg_symbol=self.symbol)
                    assert np.allclose(new_Modeler.all_data, all_data)
                    assert np.allclose(new_Modeler.all_background, all_bg)
                    assert np.allclose(new_Modeler.rois, all_rois)
                    assert np.allclose(new_Modeler.pids, all_pids)
                    assert np.allclose(new_Modeler.all_trusted, all_trusted)
                    assert np.allclose(new_Modeler.roi_id, all_roi_id)
                    LOGGER.info("Gathered file approved!")

            self.Hi[shot_idx] = shot_modeler.Hi
            self.Hi_asu[shot_idx] = shot_modeler.Hi_asu

            LOGGER.info("EVENT: DONE LOADING ROI")
            shot_modeler.ucell_man = UcellMan
            self.SIM.num_ucell_param = len(shot_modeler.ucell_man.variables)  # for convenience

            loaded_spectra = False
            if self.params.spectrum_from_imageset:
                try:
                    shot_spectra = hopper_utils.downsamp_spec(self.SIM, self.params, expt, return_and_dont_set=True)
                    loaded_spectra = True
                except Exception as err:
                    LOGGER.warning("spectrum_from_imageset is set to True, however failed to load spectra: %s" % err)
                    loaded_spectra = False

            if not loaded_spectra:
                if "spectrum_filename" in list(exper_dataframe) and exper_dataframe.spectrum_filename.values[0] is not None:
                    shot_spectra = utils.load_spectra_from_dataframe(exper_dataframe)
                    LOGGER.debug("Loaded specta from %s" % exper_dataframe.spectrum_filename.values[0])
                    shot_modeler.spec_name = exper_dataframe.spectrum_filename.values[0]

                else:
                    total_flux = exper_dataframe.total_flux.values[0]
                    if total_flux is None:
                        total_flux = self.params.simulator.total_flux
                    shot_spectra = [(expt.beam.get_wavelength(), total_flux)]

            shot_modeler.spectra = shot_spectra
            if self.params.refiner.gather_dir is not None and not self.params.refiner.load_data_from_refl:
                spec_wave, spec_weights = map(np.array, zip(*shot_spectra))
                spec_filename = os.path.splitext(os.path.basename(exper_name))[0]
                spec_filename = os.path.join(self.params.refiner.gather_dir, spec_filename+".lam")
                utils.save_spectra_file(spec_filename, spec_wave, spec_weights)
                LOGGER.info("saved spectra filename %s" % spec_filename)

            LOGGER.info("Will simulate %d energy channels" % len(shot_spectra))

            if "detz_shift_mm" in list(exper_dataframe):
                shot_modeler.originZ_init = exper_dataframe.detz_shift_mm.values[0]*1e-3
            else:
                shot_modeler.originZ_init = 0
            # TODO: is there a reason these 3 attribs are set once more after being set above?
            shot_modeler.exper_name = exper_name
            shot_modeler.exper_idx = exper_id
            shot_modeler.refl_name = refl_name

            shot_panel_groups_refined = self.determine_refined_panel_groups(shot_modeler.pids)
            rank_panel_groups_refined = rank_panel_groups_refined.union(set(shot_panel_groups_refined))

            shot_idx += 1
            self._mem_usage("Process memory usage")
            LOGGER.info("Finished loading image %d / %d (%d / %d)"
                  % (i_df+1, len(exper_names), i_work+1, len(worklist))) #, flush=True)

            shot_modeler.PAR = PAR_from_params(self.params, expt, best=exper_dataframe)
            self.Modelers[i_df] = shot_modeler  # TODO: verify that i_df as a key is ok everywhere

        LOGGER.info("DONE LOADING DATA; ENTER BARRIER")
        COMM.Barrier()
        LOGGER.info("DONE LOADING DATA; EXIT BARRIER")
        self.shot_roi_darkRMS = None

        # TODO warn that per_spot_scale refinement not intended for ensemble mode
        all_refined_groups = COMM.gather(rank_panel_groups_refined)
        panel_groups_refined = None
        if COMM.rank == 0:
            panel_groups_refined = set()
            for set_of_panels in all_refined_groups:
                panel_groups_refined = panel_groups_refined.union(set_of_panels)
        self.panel_groups_refined = list(COMM.bcast(panel_groups_refined))
        self._mem_usage("Mem after panel groups")

        LOGGER.info("EVENT: Gathering global HKL information")
        try:
            self._gather_Hi_information()
        except TypeError:
            pass
        LOGGER.info("EVENT: FINISHED gather global HKL information")
        if self.params.roi.cache_dir_only:
            print("Done creating cache directory and cache_dir_only=True, so goodbye.")
            sys.exit()
        self._mem_usage("Mem after gather Hi info")

        # in case of GPU
        LOGGER.info("BEGIN DETERMINE MAX PIX")
        self.NPIX_TO_ALLOC = self._determine_per_rank_max_num_pix()
        # TODO in case of randomize devices, shouldnt this be total max across all ranks?
        n = COMM.gather(self.NPIX_TO_ALLOC)
        if COMM.rank == 0:
            n = max(n)
        self.NPIX_TO_ALLOC = COMM.bcast(n)
        LOGGER.info("DONE DETERMINE MAX PIX")
        self._mem_usage("Mem after determine max num pix")

        self.DEVICE_ID = COMM.rank % self.params.refiner.num_devices
        self._mem_usage("Mem after load_inputs")

    def _mem_usage(self, msg="Memory usage"):
        memMB = get_memory_usage()
        host = socket.gethostname()
        LOGGER.info("%s: %f GB on node %s" % (msg, memMB / 1024, host))

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
            #print("Rank %d, shot %d has %d pixels" % (COMM.rank, i_shot+1, npix))
        LOGGER.info("Rank %d, max pix to be modeled: %d" % (COMM.rank, max_npix))
        return max_npix

    def _try_loading_spectrum_filelist(self):
        file_list = None
        fpath = self.params.simulator.spectrum.filename_list
        if fpath is not None:
            file_list = [l.strip() for l in open(fpath, "r").readlines()]
            assert all([len(l.split()) == 1 for l in file_list]), "weird spectrum file %s"% fpath
        return file_list

    def _gather_Hi_information(self):
        # aggregate all miller indices

        #self.Hi_asu_all_ranks = []
        # TODO assert list types are stored in Hi and Hi_asu
        #for i_shot in self.Hi: #range(nshots_on_this_rank):
        #    self.Hi_asu_all_ranks += self.Hi_asu[i_shot]

        Hi_asu_possible = self._get_possible_Hi_asu()
        Hi_asu = chain.from_iterable(self.Hi_asu.values())
        Hi_asu_counter = Counter(Hi_asu)
        Hi_asu_possible_counts = [Hi_asu_counter[k] for k in Hi_asu_possible]
        Hi_asu_possible_counts = np.array(Hi_asu_possible_counts, dtype=int)
        Hi_asu_possible_counts_global = np.zeros_like(Hi_asu_possible_counts, dtype=int)
        COMM.Allreduce(Hi_asu_possible_counts, Hi_asu_possible_counts_global, op=MPI.SUM)
        unique_asu_all_ranks = []
        Hi_asu_all_ranks = []
        for p, c in Hi_asu_possible, Hi_asu_possible_counts_global:
            if c:
                unique_asu_all_ranks.append(p)
                Hi_asu_all_ranks.extend([p] * c)

        # TODO: Unused Derek's code, too MPI-intensive - reference, to be removed
        # Hi_asu_all_ranks = COMM.gather(self.Hi_asu)
        # unique_asu_all_ranks = None
        # if COMM.rank==0:
        #     temp = []
        #     for Hi_asu in Hi_asu_all_ranks:
        #         for i_shot in Hi_asu:
        #             temp += Hi_asu[i_shot]
        #     Hi_asu_all_ranks = temp
        #     unique_asu_all_ranks = set(Hi_asu_all_ranks)
        # unique_asu_all_ranks = COMM.bcast(unique_asu_all_ranks)
        # TODO: END of Old Derek's code to be removed - Daniel Tchon

        # Hi_asu_all_ranks should be None on all ranks except rank0
        marr_unique_h = self._get_unique_Hi(Hi_asu_all_ranks)

        # this will map the measured miller indices to their index in the LBFGS parameter array self.x
        self.idx_from_asu = {h: i for i, h in enumerate(unique_asu_all_ranks)}
        # we will need the inverse map during refinement to update the miller array in diffBragg, so we cache it here
        self.asu_from_idx = {i: h for i, h in enumerate(unique_asu_all_ranks)}

        self.num_hkl_global = len(self.idx_from_asu)

        fres = marr_unique_h.d_spacings()
        self.res_from_asu = {h: res for h, res in zip(fres.indices(), fres.data())}

    def _get_first_modeller_symmetry(self):
        ii = list(self.Modelers.keys())[0]
        uc = self.Modelers[ii].ucell_man
        params = uc.a, uc.b, uc.c, uc.al * 180 / np.pi, uc.be * 180 / np.pi, uc.ga * 180 / np.pi
        if self.params.refiner.force_unit_cell is not None:
            params = self.params.refiner.force_unit_cell
        return crystal.symmetry(unit_cell=params, space_group_symbol=self.symbol)

    def _get_possible_Hi_asu(self):
        if COMM.rank == 0:
            symm = self._get_first_modeller_symmetry()
            res_ranges = utils.parse_reso_string(self.params.refiner.res_ranges)
            d_min = min([d_min for d_min, _ in res_ranges]) * 0.8  # accomodate variations in unit cell
            mset_full = symm.build_miller_set(anomalous_flag=True, d_min=d_min)
            Hi_asu_possible = list(mset_full.indices())
        else:
            Hi_asu_possible = None
        return COMM.bcast(Hi_asu_possible, root=0)

    def _get_unique_Hi(self, Hi_asu_all_ranks):
        if COMM.rank == 0:
            from cctbx.array_family import flex as cctbx_flex

            symm = self._get_first_modeller_symmetry()
            hi_asu_flex = cctbx_flex.miller_index(Hi_asu_all_ranks)
            mset = miller.set(symm, hi_asu_flex, anomalous_flag=True)
            marr = miller.array(mset)
            binner = marr.setup_binner(d_max=self.params.refiner.stage_two.d_max, d_min=self.params.refiner.stage_two.d_min,
                                       n_bins=self.params.refiner.stage_two.n_bin)
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
            unique_Hi_asu = set(Hi_asu_all_ranks)
            hi_flex_unique = cctbx_flex.miller_index(list(unique_Hi_asu))
            mset = miller.set(symm, hi_flex_unique, anomalous_flag=True)
            self.binner = mset.setup_binner(d_min=self.params.refiner.stage_two.d_min,
                                            d_max=self.params.refiner.stage_two.d_max,
                                            n_bins=self.params.refiner.stage_two.n_bin)
            mset.completeness(use_binning=True).show()
            marr_unique_h = miller.array(mset)
            print("Rank %d: total miller vars=%d" % (COMM.rank, len(unique_Hi_asu)))
        else:
            marr_unique_h = None

        marr_unique_h = COMM.bcast(marr_unique_h)
        return marr_unique_h

    def _launch(self):
        """
        Usually this method should be modified when new features are added to refinement
        """
        # TODO return None or refiner instance
        LOGGER.info("begin _launch")
        x_init = None
        nmacro = self.params.refiner.num_macro_cycles
        n_trials = len(self.params.refiner.max_calls)

        for i_trial in range(n_trials*nmacro):

            self.RUC = StageTwoRefiner(self.Modelers, self.symbol, self.params)

            if self.will_refine(self.params.refiner.refine_spot_scale):
                self.RUC.refine_crystal_scale = (self.params.refiner.refine_spot_scale*nmacro)[i_trial]

            if self.will_refine(self.params.refiner.refine_Fcell):
                self.RUC.refine_Fcell = (self.params.refiner.refine_Fcell*nmacro)[i_trial]

            self.RUC.panel_group_from_id = self.panel_group_from_id
            self.RUC.panel_reference_from_id = self.panel_reference_from_id
            self.RUC.panel_groups_being_refined = self.panel_groups_refined

            # TODO verify not refining Fcell in case of local refiner
            self.RUC.max_calls = (self.params.refiner.max_calls*nmacro)[i_trial]
            self.RUC.x_init = x_init
            self.RUC.ignore_line_search_failed_step_at_lower_bound = True  # TODO: why was this necessary?

            # plot things
            self.RUC.trial_id = i_trial

            self.RUC.log_fcells = True
            self.RUC.request_diag_once = False
            self.RUC.trad_conv = True
            self.RUC.idx_from_asu = self.idx_from_asu
            self.RUC.asu_from_idx = self.asu_from_idx

            self.RUC.S = self.SIM
            self.RUC.restart_file = self.params.refiner.io.restart_file
            self.RUC.S.update_nanoBragg_instance('update_oversample_during_refinement',
                                                 self.params.refiner.update_oversample_during_refinement)
            self.RUC.S.update_nanoBragg_instance("Npix_to_allocate", self.NPIX_TO_ALLOC)
            self.RUC.S.update_nanoBragg_instance('device_Id', self.DEVICE_ID)
            self.RUC.use_curvatures_threshold = self.params.refiner.use_curvatures_threshold
            if not self.params.refiner.curvatures:
                self.RUC.S.update_nanoBragg_instance('compute_curvatures', False)
            if COMM.rank==0:
                self.RUC.S.update_nanoBragg_instance('verbose', self.params.refiner.verbose)

            LOGGER.info("_launch run setup")
            self._mem_usage("Mem usage before _setup")
            self.RUC.run(setup_only=True)
            self._mem_usage("Mem usage after _setup")
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

            LOGGER.info("_launcher running optimization")

            self.RUC.run(setup=False)
            LOGGER.info("_launcher done running optimization")
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

            if self.params.profile:
                self.RUC.S.D.show_timings(self.RUC.rank)
            if os.environ.get("DIFFBRAGG_USE_CUDA") is not None:
                self.RUC.S.D.gpu_free()

    def will_refine(self, param):
        return param is not None and any(param)
