from __future__ import absolute_import, division, print_function
from scipy.interpolate import interp1d

from simtbx.command_line.hopper import get_data_model_pairs
from cctbx import sgtbx, miller
import time
from collections import Counter
from scipy.optimize import basinhopping
from simtbx.diffBragg import hopper_utils
import h5py
try:
    import pandas
except ImportError:
    print("Please intsall pandas, libtbx.python -m pip install pandas")
    exit()
import pandas
from scitbx.matrix import sqr
try:
    from stream_redirect import Redirect
except ImportError:
    pass

# diffBragg internal parameter indices
FHKL_ID = 11

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.ensemble_more

import numpy as np
np.seterr(invalid='ignore')
import os
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.mpi4py import MPI

COMM = MPI.COMM_WORLD
from libtbx.phil import parse

from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners.parameters import RangedParameter

ensemble_phil = """
include scope simtbx.command_line.hopper.phil_scope
wilson_fat = 1
  .type = float
  .help = scales the width of wilson distribution
use_wilson_restraints = None
  .type = bool
  .help = if True, must provide scaling_reference_mtz information below
scaling_reference_mtz_name = None
  .type = str
  .help = path to mtz file containing amplitudes
  .help = that will be used to restrain the mean squared structure factor
  .help = intensity at each resolution
scaling_reference_mtz_col = None
  .type = str
  .help = mtz column name
sigma_frac = None
  .type = float
  .help = sigma for Fhkl restraints will be some fraction of the starting value
sanity_test_hkl_variation = False
  .type = bool
  .help = measure the variation of each HKL within the shoebox
sanity_test_models = False
  .type = bool
  .help = make sure best models from stage 1 are reproduced at the start
sanity_test_amplitudes = False
  .type = bool
  .help = if True, then quickly run a sanity check ensuring that all h,k,l are predicted
  .help = and/or in the starting miller array
x_write_freq = 25
  .type = int
  .help = save x arrays every x_write_freq iterations
percentile_cut = None
  .type = float
  .help = percentile below which pixels are masked
space_group = P6522
  .type = str
  .help = space group to refine structure factors in
first_n = None
  .type = int
  .help = refine the first n shots only
maxiter = None
  .type = int
  .help = stop refiner after this many iters
disp = False
  .type = bool
  .help = scipy minimize convergence printouts
"""

phil_scope = parse(ensemble_phil, process_includes=True)


class Script:
    def __init__(self):
        from dials.util.options import OptionParser

        self.params = self.parser = None
        if COMM.rank == 0:
            self.parser = OptionParser(
                usage="",  # stage 1 (per-shot) diffBragg refinement",
                sort_options=True,
                phil=phil_scope,
                read_experiments=True,
                read_reflections=True,
                check_format=False,
                epilog="PyCuties")
        self.parser = COMM.bcast(self.parser)
        if COMM.rank == 0:
            self.params, _ = self.parser.parse_args(show_diff_phil=True)
        self.params = COMM.bcast(self.params)

    def run(self):
        assert os.path.exists(self.params.exp_ref_spec_file)
        assert self.params.best_pickle is not None
        input_lines = None
        best_models = None
        if COMM.rank == 0:
            input_lines = open(self.params.exp_ref_spec_file, "r").readlines()
            if self.params.sanity_test_input:
                for line in input_lines:
                    for fname in line.strip().split():
                        if not os.path.exists(fname):
                            raise FileNotFoundError("File %s not there " % fname)
            #if self.params.best_pickle is not None:
            #    if not self.params.quiet: print("reading pickle %s" % self.params.best_pickle)
            #    best_models = pandas.read_pickle(self.params.best_pickle)
        input_lines = COMM.bcast(input_lines)
        #best_models = COMM.bcast(best_models)
        if self.params.best_pickle is not None:
            if not self.params.quiet: print("reading pickle %s" % self.params.best_pickle)
            best_models = pandas.read_pickle(self.params.best_pickle)
        if self.params.skip is not None:
            input_lines = input_lines[self.params.skip:]

        if self.params.first_n is not None:
            input_lines = input_lines[:self.params.first_n]

        shot_roi_dict = count_rois(input_lines, self.params.quiet)
        # gether statistics, e.g. how many total ROIs
        nshots = len(shot_roi_dict)
        nrois = sum([len(shot_roi_dict[s]) for s in shot_roi_dict])
        if not self.params.quiet:print("Rank %d will load %d rois across %d shots" % (COMM.rank, nrois, nshots))

        # make a data modeler for each shot
        Modelers = {}
        bests = {}
        for i_exp in shot_roi_dict:
            if not self.params.quiet: print("COMM.rank %d on shot  %d" % (COMM.rank, i_exp + 1))

            # this corresponds to the expfile, reflfile, specfile for this shot
            line = input_lines[i_exp]
            # which rois to load from this shots refls
            rois_to_load = shot_roi_dict[i_exp]

            # the filenames
            exp, ref, spec = line.strip().split()

            # if there is a starting model, then load it
            best = None
            if best_models is not None:
                best = best_models.query("exp_name=='%s'" % exp)
                if len(best) == 0:
                    best = best_models.query("opt_exp_name=='%s'" % exp)
                if len(best) != 1:
                    raise ValueError("Should be 1 entry for exp %s in best pickle %s" % (exp, self.params.best_pickle))

            # dont think this is necessary, but doesnt matter
            self.params.simulator.spectrum.filename = spec

            # each shot gets a data modeler
            Modeler = DataModeler(self.params)
            # gather the data from the input files
            if not Modeler.GatherFromExperiment(exp, ref, self.params.space_group, rois_to_load):
                continue

            # store the modeler for later use(each rank has one modeler per shot in shot_roi_dict)
            Modeler.exp_name = exp
            Modelers[i_exp] = Modeler
            bests[i_exp] = best  # TODO questionable

        if not Modelers:
            raise ValueError("Need at least 1 modeler per rank, try reducing number of ranks.")

        # count up the total number of pixels being modeled by this rank
        npix = [len(modeler.all_data) for modeler in Modelers.values()]
        if not self.params.quiet: print("Rank %d wil model %d pixels in total" %(COMM.rank, sum(npix)), flush=True)
        COMM.barrier()

        # these are the experient ids corresponding to exper-ref-spectrum input file lines , for this rank
        i_exps = list(Modelers.keys())

        # make a SIM instance, use first Modeler as a template
        self.SIM = get_diffBragg_simulator(Modelers[i_exps[0]].E, self.params)
        setup_Fhkl_attributes(self.SIM, self.params, Modelers)

        # iterate over the shot-modelers on this rank, and set the
        # free parameters and spectrum for each shot
        for i_exp, Modeler in Modelers.items():
            # load spectra for other shots
            if bests[i_exp] is not None:
                total_flux = bests[i_exp].total_flux.values[0]
                spectrum_stride = bests[i_exp].spectrum_stride.values[0]
            else:
                total_flux = self.params.simulator.total_flux
                spectrum_stride = self.params.simulator.spectrum.stride
            spectra_file = input_lines[i_exp].strip().split()[2]
            spectrum = utils.load_spectra_file(spectra_file, total_flux, spectrum_stride, as_spectrum=True)

            # set parameter objects for this shot
            Modeler.SimulatorParamsForExperiment(bests[i_exp])
            #Modeler.spectrum = spectrum  # store the spectrum as part of the modeler
            self.SIM.beam.spectrum = spectrum
            Modeler.spectrum = spectrum
            Modeler.xray_beams = self.SIM.beam.xray_beams

            # define the fcell global index for each modeled pixel
            Modeler.all_fcell_global_idx = np.array([self.SIM.i_fcell_from_asu[h] for h in Modeler.hi_asu_perpix])
            Modeler.unique_i_fcell = set(Modeler.all_fcell_global_idx)
            Modeler.is_i_fcell = {}
            for i_fcell in Modeler.unique_i_fcell:
                sel = Modeler.all_fcell_global_idx == i_fcell
                Modeler.is_i_fcell[i_fcell] = sel

            Modeler.sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon

            Modeler.shiftZ_meters = 0
            if bests[i_exp] is not None and "detz_shift_mm" in list(bests[i_exp]):
                    Modeler.shiftZ_meters = bests[i_exp].detz_shift_mm.values[0]*1e-3

        self.SIM.update_Fhkl = Fhkl_updater(self.SIM, Modelers)

        # there is variable 1 scale factor per shot
        # Bookkeeping:
        # each i_exp in shot_roi_dict should globally point to a single index
        # MORE INFO: i_exp is a line number pointing to an experiment
        # we need a mapping from i_exp to a number on the range 0 to N-1 where N is the total number
        # of shots being modeled. If a shot fails to load in the GatherFromExperiment step, then i_exp
        # will no longer form a uniform range of values from 0 to N-1, hence why we need this mapping
        shot_mapping = {}
        rank_exp_indices = COMM.gather(i_exps, root=0)
        ntimes = None
        if COMM.rank == 0:
            all_indices = [i_exp for indices in rank_exp_indices for i_exp in indices]
            # count how many ranks a shot is divided amongst, this number
            # should be reported back to all ranks to allow ease in computing
            # the restraint terms
            ntimes = Counter(all_indices)  # how many times a shot appears on a ranks list of shots, this is used later to normalize restraint terms
            shot_mapping = {i_exp: ii for ii, i_exp in enumerate(set(all_indices))}
        shot_mapping = COMM.bcast(shot_mapping)  # this now maps each i_exp to a number on [0-1)
        ntimes = COMM.bcast(ntimes)

        # for GPU usage, allocate enough pixels! We only allocate once,
        # and shots are modeled one-at-a-time , so allocate enough pixels to model the shot with the most pixels
        self.NPIX_TO_ALLOC = determine_per_rank_max_num_pix(Modelers)
        # TODO in case of randomize devices, shouldnt this be total max across all ranks?
        n = COMM.gather(self.NPIX_TO_ALLOC)
        if COMM.rank == 0:
            n = max(n)
        self.NPIX_TO_ALLOC = COMM.bcast(n)
        self.SIM.D.Npix_to_allocate = int(self.NPIX_TO_ALLOC)

        global_Nshots = len(shot_mapping)  # total number of shots being modeled
        nparam_per_shot = 7   # 1 scale factor, 3 Nabc terms, 3 RotXYZ terms
        total_params = nparam_per_shot*global_Nshots + self.SIM.n_global_fcell  # total refinement paramters
        if COMM.rank==0:
            print("Refining %d parameters in total" % total_params, flush=True)

        # here we organize the global parameter array which includes per-shot params and Fhkl params shared amongst shots
        total_pershot_params = global_Nshots*nparam_per_shot
        # the term xidx refers to the index (idx) of the parameter in the parameter array (x refers to parameter array)
        # The Fhkl parameters are always stored at the end of the global parameter array
        Fhkl_xidx = list(range(total_pershot_params, total_pershot_params + self.SIM.n_global_fcell))

        # each rank has different modelers, and those will need to each reference different
        # portions of the global parameter array, and rank_xidx holds the referencing information
        rank_xidx = {}
        # Also, in this loop we will get the parameter objects and send em to rank0 for writing to disk
        param_data = []
        for i_exp in i_exps:
            xidx_start = shot_mapping[i_exp]*nparam_per_shot
            xidx = list(range(xidx_start, xidx_start+nparam_per_shot))
            xidx += Fhkl_xidx
            rank_xidx[i_exp] = xidx

            Scale_param = Modelers[i_exp].PAR.Scale
            Nabc_params = Modelers[i_exp].PAR.Nabc_params
            param_data.append([shot_mapping[i_exp], Scale_param, Nabc_params])

        mpi_safe_makedirs(self.params.outdir)
        mpi_safe_makedirs(os.path.join(self.params.outdir, "x"))

        all_param_data = COMM.reduce(param_data)
        if COMM.rank == 0:
            tsave = time.time()

            all_i_shots, all_Scales, all_Nabc_params = zip(*all_param_data)

            # sometimes a shot is distributed across multiple ranks, hence there are duplicates
            # in all_scale_param_data, so we will remove the duplicates here
            output_params = {}
            for i_shot, Scale, Nabc_params in zip(all_i_shots, all_Scales, all_Nabc_params):
                if i_shot not in output_params:
                    output_params[i_shot] = Scale, Nabc_params
                else:
                    assert output_params[i_shot][0].init == Scale.init
                    for i_N in range(3):
                        assert output_params[i_shot][1][i_N].init == Nabc_params[i_N].init

            # at this point there should be a single scale per shot! We verify that here:
            ordered_shot_inds = np.sort(list(output_params.keys()))
            assert np.all(ordered_shot_inds == np.arange(global_Nshots))
            scale_params = [output_params[i_shot][0] for i_shot in range(global_Nshots)]

            # save the parameter objects to pickle files using numpy
            scale_params_file = os.path.join(self.params.outdir, "x", "scale_params.npy")
            print("Saving scale parameter objects to %s" % scale_params_file)
            np.save(scale_params_file, scale_params)
            Nabc_params_file = os.path.join(self.params.outdir, "x", "Nabc_params.npz")
            print("Saving Nabc parameter objects to %s" % Nabc_params_file)
            N_params = []
            for i_N in range(3):
                params_i_N = [output_params[i_shot][1][i_N] for i_shot in range(global_Nshots)]
                N_params.append(params_i_N)
            np.savez(Nabc_params_file, Na=N_params[0], Nb=N_params[1], Nc=N_params[2])

            Fhkl_params_file =os.path.join(self.params.outdir, "x", "Fhkl_params.npy")
            print("Saving Fhkl parameter objects to %s" % Fhkl_params_file)
            np.save(Fhkl_params_file, self.SIM.Fhkl_modelers)
            tsave = time.time()-tsave
            print("Time to save parameter arrays: %f sec" % tsave)
            # These files can be used later in conjunction with the x-array output
        COMM.barrier()

        # initial parameter values, all parameters initialize at 1 due to rescaling
        x0 = np.array([1] * total_params)

        # save the initial model pixels values
        for i_exp in Modelers.keys():
            M = Modelers[i_exp]
            #self.SIM.D.nopolar=True

            if self.params.sanity_test_hkl_variation:
                model(x0[rank_xidx[i_exp]], self.SIM, M, compute_grad=False, sanity_test=0)
                continue
            else:
                _, _, best_mod = model(x0[rank_xidx[i_exp]], self.SIM, M, compute_grad=False)
            #self.SIM.D.nopolar=False
            if self.params.sanity_test_models:
                # NOTE works on one MPI rank only, because shots are divided across ranks
                assert COMM.size==1
                data_subimg, model_subimg, strong_subimg, bragg_subimg = get_data_model_pairs(M.rois, M.pids,
                                                                                              M.roi_id, best_mod-M.all_background,
                                                                                              M.all_data,
                                                                                              background=M.all_background)
                assert bests[i_exp] is not None
                #assert "best_model_image_file" in list(bests[i_exp])
                assert "stage1_output_img" in list(bests[i_exp])
                h = h5py.File(bests[i_exp].stage1_output_img.values[0], 'r')
                nroi = h['rois'].shape[0]
                print("Testing %s" % M.exp_name)
                for i in range(nroi):
                    assert h['pids'][i] == M.pids[i]
                    assert np.all(h['rois'][i] == M.rois[i])
                    ref_img = h['bragg/roi%d' %i][()]
                    current_img = bragg_subimg[i]
                    assert np.allclose(ref_img, current_img)
                    print("Model for %s ROI %d / %d is as expected!" % (M.exp_name, i+1 , nroi ))
                    continue

                    #rois_per_panel = {M.pids[i]: [M.rois[i]]}
                    #model_from_nanoBragg = utils.roi_spots_from_pandas(bests[i_exp], rois_per_panel, quiet=True,
                    #                                    mtz_file=self.params.simulator.structure_factors.mtz_name,
                    #                                    mtz_col=self.params.simulator.structure_factors.mtz_column,
                    #                                    cuda=False, reset_Bmatrix=True, nopolar=True,
                    #                                    force_no_detector_thickness=True,
                    #                                    norm_by_nsource=True)
                    #x1,x2,y1,y2 = M.rois[i]
                    #nanoBragg_img = model_from_nanoBragg[M.pids[i]][y1:y2, x1:x2]
                    #nanoBragg_img *= self.SIM.D.spot_scale
                    ##nanoBragg_img /= len(self.SIM.D.xray_beams)

                    #max_diff = np.abs(nanoBragg_img - current_img)
                    #max_diffs.append(max_diff)
                    #assert np.allclose(nanoBragg_img, current_img)
                #all_m = []
                #for m in max_diffs:
                #    m = m.ravel()
                #    all_m.append(m)
                #all_m = np.hstack(all_m)
                #import pylab as plt
                #plt.hist(all_m, bins=100)
                #ax = plt.gca()
                #ax.set_yscale("log")
                #ax.tick_params(labelsize=12)
                #plt.xlabel("photons", fontsize=14)
                #plt.ylabel("# of pixels", fontsize=14)
                #ax.grid(1, alpha=0.5, which='both')
                #plt.show()

                #print(max_diffs, np.max(max_diffs))
            if self.params.sanity_test_models:
                print("Tested models for sanity, all seem good! exiting...")
                exit()

            img_path = "rank%d_img%d_before.h5" %(COMM.rank, i_exp)
            img_path = os.path.join(self.params.outdir, img_path)
            save_model_Z(img_path, M, best_mod)
        if self.params.sanity_test_hkl_variation:
            COMM.barrier()
            if COMM.rank == 0:
                print("Sanity checked hkl variation")
            exit()

        if COMM.rank==0:
            print("MINIMIZE!", flush=True)
        min_out = Minimize(x0, rank_xidx, self.params, self.SIM, Modelers, ntimes, global_Nshots)
        x = min_out.x
        # TODO analyze the convergence data here

        # save the final models
        for i_exp in Modelers:
            M = Modelers[i_exp]
            _,_,best_mod = model(x[rank_xidx[i_exp]], self.SIM, M, compute_grad=False)
            img_path = "rank%d_img%d_after.h5" %(COMM.rank, i_exp)
            img_path = os.path.join(self.params.outdir, img_path)
            save_model_Z(img_path, M, best_mod)


def get_diffBragg_simulator(expt, params):
    SIM = utils.simulator_from_expt_and_params(expt, params)

    # this works assumes all crystals are of the same crystal system
    SIM.ucell_man = utils.manager_from_crystal(expt.crystal)

    SIM.D.no_Nabc_scale = params.no_Nabc_scale
    SIM.num_xtals = 1
    return SIM


def count_rois(lines, quiet):
    info = []
    for i_line, line in enumerate(lines):
        if i_line % COMM.size != COMM.rank:
            continue
        e,r,s = line.strip().split()
        R = flex.reflection_table.from_file(r)
        info.append((i_line, len(R)))
    info = COMM.reduce(info)
    info = COMM.bcast(info)
    shots, nref = zip(*info)
    if COMM.rank==0:
        print("Input is %d refls on %d shots" %(sum(nref), len(shots)))

    shots_and_rois_to_load = []
    if COMM.rank == 0:
        for i_line, nroi in info:
            shots_to_load = [i_line] * nroi
            rois_to_load = list(range(nroi))
            shots_and_rois_to_load += list(zip(shots_to_load, rois_to_load))

    shots_and_rois_to_load = COMM.bcast(shots_and_rois_to_load)
    shots_and_rois_to_load = np.array_split(shots_and_rois_to_load, COMM.size)[COMM.rank]
    from itertools import groupby
    gb = groupby(sorted(shots_and_rois_to_load, key=lambda x: x[0]), key=lambda x:x[0])
    shot_rois = {shot: [i_roi for _,i_roi in vals] for shot,vals in gb}
    out = "\nRank %d will model\n" %COMM.rank
    for shot in shot_rois:
        roi_s = ",".join(map(str, shot_rois[shot]))
        out += "\tShot %d; rois=%s\n" % (shot, roi_s)
    if not quiet: print(out+"\n")
    return shot_rois


class DataModeler:
    #TODO merge this class with the one in hopper_utils

    def __init__(self, params):
        """ params is a simtbx.diffBragg.hopper phil"""
        self.params = params
        self.SIM = None
        self.PAR = None
        self.spectrum = None
        self.E = None
        self.pan_fast_slow =None
        self.all_background =None
        self.roi_id =None
        self.u_id = None
        self.all_data =None
        self.all_sigmas =None
        self.all_trusted =None
        self.npix_total =None
        self.all_fast =None
        self.all_slow =None
        self.rois=None
        self.pids=None
        self.refls_idx = None
        self.all_refls_idx = None
        self.tilt_abc=None
        self.selection_flags=None
        self.background=None
        self.tilt_cov = None
        self.simple_weights = None
        self.ref_id = None
        self.ref_name = None

        self.Hi_asu = None
        self.hi_asu_perpix = None
        self.all_nominal_hkl = []





    def GatherFromExperiment(self, exp, ref, sg_symbol, ref_indices=None):
        #TODO delete the background attribute because its large(full image)
        """
        :param exp: input experiment filename
        :param ref: intput reflection filename
        :param ref_indices: integer list corresponding to which reflections to load
        :return: True or False depending on success or failure
        """
        self.E = ExperimentListFactory.from_json_file(exp)[0]
        if self.params.opt_det is not None:
            opt_det_E = ExperimentListFactory.from_json_file(self.params.opt_det, False)[0]
            self.E.detector = opt_det_E.detector
        self.ref_name = ref
        refls = flex.reflection_table.from_file(ref)
        refls["refls_idx"] = flex.int(list(range(len(refls))))
        if ref_indices is not None:
            refl_sel = np.zeros(len(refls), bool)
            for i_roi in ref_indices:
                refl_sel[i_roi] = True
            refls = refls.select(flex.bool(refl_sel))
        img_data = utils.image_data_from_expt(self.E)
        img_data /= self.params.refiner.adu_per_photon
        is_trusted = utils.load_mask(self.params.roi.hotpixel_mask)
        hotpix_mask = None
        if is_trusted is not None:
            hotpix_mask = ~is_trusted
        self.sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon
        roi_packet = utils.get_roi_background_and_selection_flags(
            refls, img_data, shoebox_sz=self.params.roi.shoebox_size,
            reject_edge_reflections=self.params.roi.reject_edge_reflections,
            reject_roi_with_hotpix=self.params.roi.reject_roi_with_hotpix,
            background_mask=None, hotpix_mask=hotpix_mask,
            bg_thresh=self.params.roi.background_threshold,
            use_robust_estimation=not self.params.roi.fit_tilt,
            set_negative_bg_to_zero=self.params.roi.force_negative_background_to_zero,
            pad_for_background_estimation=self.params.roi.pad_shoebox_for_background_estimation,
            sigma_rdout=self.sigma_rdout, deltaQ=self.params.roi.deltaQ, experiment=self.E,
            weighted_fit=self.params.roi.fit_tilt_using_weights,
            tilt_relative_to_corner=self.params.relative_tilt, ret_cov=True)

        self.rois, self.pids, self.tilt_abc, self.selection_flags, self.background, self.tilt_cov = roi_packet

        assert len(self.rois) == len(refls)

        if sum(self.selection_flags) == 0:
            if not self.params.quiet: print("No pixels slected, continuing")
            return False

        self.rois = [roi for i_roi, roi in enumerate(self.rois) if self.selection_flags[i_roi]]
        self.tilt_abc = [abc for i_roi, abc in enumerate(self.tilt_abc) if self.selection_flags[i_roi]]
        self.pids = [pid for i_roi, pid in enumerate(self.pids) if self.selection_flags[i_roi]]

        self.refls_idx = [refls[i_roi]["refls_idx"] for i_roi in range(len(refls)) if self.selection_flags[i_roi]]
        self.tilt_cov = [cov for i_roi, cov in enumerate(self.tilt_cov) if self.selection_flags[i_roi]]
        #self.refl_index = [refls[i_roi]["refl_index"] for i_roi in range(len(self.rois)) if self.selection_flags[i_roi]]

        # TODO assert that the order of refls and selection_flags was maintained ?
        refls = refls.select(flex.bool(self.selection_flags))
        Hi = list(refls["miller_index"])
        self.Hi_asu = utils.map_hkl_list(Hi, True, sg_symbol)

        all_data = []
        all_refls_idx = []
        hi_asu_perpix = []
        all_pid = []
        all_fast = []
        all_slow = []
        all_trusted = []
        all_sigmas = []
        all_background = []
        roi_id = []
        self.all_nominal_hkl = []
        #refl_id = []
        for i_roi in range(len(self.rois)):
            pid = self.pids[i_roi]
            x1, x2, y1, y2 = self.rois[i_roi]
            Y, X = np.indices((y2 - y1, x2 - x1))
            data = img_data[pid, y1:y2, x1:x2].copy()

            all_fast += list(X.ravel() + x1)
            all_slow += list(Y.ravel() + y1)

            data = data.ravel()
            trusted = is_trusted[pid, y1:y2, x1:x2].ravel()

            if self.params.percentile_cut is not None:
                lower_cut = np.percentile(data, self.params.percentile_cut)
                trusted[data < lower_cut] = False

            all_background += list(self.background[pid, y1:y2, x1:x2].ravel())
            all_trusted += list(trusted)
            all_sigmas += list(np.sqrt(data + self.sigma_rdout ** 2))
            all_data += list(data)
            npix = len(data)  # np.sum(trusted)
            all_pid += [pid] * npix
            all_refls_idx += [self.refls_idx[i_roi]] * npix
            hi_asu_perpix += [self.Hi_asu[i_roi]] * npix
            roi_id += [i_roi] * npix
            self.all_nominal_hkl += [tuple(Hi[i_roi])]*npix  # this is the nominal l component of the miller index in the P1 setting
            #refl_id += [self.refl_index[i_roi]] * npix
        pan_fast_slow = np.ascontiguousarray((np.vstack([all_pid, all_fast, all_slow]).T).ravel())
        self.pan_fast_slow = flex.size_t(pan_fast_slow)
        self.hi_asu_perpix = hi_asu_perpix
        self.all_background = np.array(all_background)
        self.roi_id = np.array(roi_id)
        #self.refl_id = np.array(refl_id)
        self.all_data = np.array(all_data)
        self.all_sigmas = np.array(all_sigmas)
        # note rare chance for sigmas to be nan if the args of sqrt is below 0
        self.all_trusted = np.logical_and(np.array(all_trusted), ~np.isnan(all_sigmas))
        self.npix_total = len(all_data)
        self.simple_weights = 1/self.all_sigmas**2
        self.u_id = set(self.roi_id)
        self.all_refls_idx = np.array(all_refls_idx)
        return True

    def SimulatorParamsForExperiment(self, best=None):
        """optional best parameter is a single row of a pandas datafame containing the starting
        models, presumably optimized from a previous minimzation using this program"""
        ParameterType = RangedParameter
        PAR = SimParams()

        # per shot Scale factor
        p = ParameterType()
        p.sigma = self.params.sigmas.G
        if best is not None:
            p.init = best.spot_scales.values[0]
        else:
            p.init = self.params.init.G
        p.minval = self.params.mins.G
        p.maxval = self.params.maxs.G
        p.fix = self.params.fix.G
        PAR.Scale = p

        PAR.Nabc_params = []
        for i_N in range(3):
            p = ParameterType()
            p.sigma = self.params.sigmas.Nabc[i_N]
            if best is not None:
                p.init = best.ncells.values[0][i_N]
            else:
                p.init = self.params.init.Nabc[i_N]
            p.minval = self.params.mins.Nabc[i_N]
            p.maxval = self.params.maxs.Nabc[i_N]
            p.fix = self.params.fix.Nabc
            PAR.Nabc_params.append(p)

        PAR.RotXYZ_params = []
        for i_rot in range(3):
            p = ParameterType()
            p.sigma = self.params.sigmas.RotXYZ[i_rot]
            #if best is not None:
            #    p.init = best.values[0][i_rot]
            #else:
            #    p.init = self.params.init.RotXYZ[i_rot]
            p.init = 0  # design choice, rotation parameters are perturbations from the Umat, thus we always init as 0
            p.minval = self.params.mins.RotXYZ[i_rot]
            p.maxval = self.params.maxs.RotXYZ[i_rot]
            p.fix = self.params.fix.RotXYZ
            PAR.RotXYZ_params.append(p)

        # each modeler has an experiment self.E, with a crystal attached
        if best is not None:
            self.E.crystal.set_A(best.Amats.values[0])
        PAR.Umatrix = sqr(self.E.crystal.get_U())
        PAR.Bmatrix = sqr(self.E.crystal.get_B())

        self.PAR = PAR

def Minimize(x0, rank_xidx, params, SIM, Modelers, ntimes, nshots_total):
    if params.refiner.randomize_devices is not None:
        dev = np.random.choice(params.refiner.num_devices)
    else:
        dev = COMM.rank % params.refiner.num_devices
    SIM.D.device_Id = dev

    vary = np.ones(len(x0), bool)
    nparam_per_shot = 7  # 1 scale, 3 Nabc, 3 RotXYZ  # TODO make param order in agreement with hopper_utils
    for i_shot in range(nshots_total):
        if params.fix.G:
            vary[i_shot * nparam_per_shot] = False
        if params.fix.RotXYZ:
            vary[i_shot*nparam_per_shot + 1] = False
            vary[i_shot * nparam_per_shot + 2] = False
            vary[i_shot * nparam_per_shot + 3] = False
        if params.fix.Nabc:
            vary[i_shot * nparam_per_shot + 4] = False
            vary[i_shot * nparam_per_shot + 5] = False
            vary[i_shot * nparam_per_shot + 6] = False

    #n_uc_param = len(SIM.ucell_man.variables)
    #if params.fix.ucell:
    #    for i_uc in range(n_uc_param):
    #        vary[SIM.num_xtals*n_param_per_xtal + i_uc] = False
    #if params.fix.detz_shift:
    #    vary[SIM.num_xtals*n_param_per_xtal + n_uc_param] = False

    target = TargetFunc(params, SIM)

    target.vary = vary  # fixed flags
    target.x0 = np.array(x0, np.float64)  # initial full parameter list
    x0_for_refinement = target.x0[vary]

    niter = params.niter
    SIM.D.refine(hopper_utils.FHKL_ID)
    if not params.fix.Nabc:
        SIM.D.refine(hopper_utils.NCELLS_ID)
    if not params.fix.RotXYZ:
        for ROT_ID in hopper_utils.ROTXYZ_IDS:
            SIM.D.refine(ROT_ID)

    if params.method in ["Nelder-Mead", "Powell"]:
        args = (rank_xidx, SIM, Modelers, True, params, False, ntimes)
        jac = None
    else:
        args = (rank_xidx, SIM, Modelers, True, params, True, ntimes)
        jac = target.jac
    out = basinhopping(target, x0_for_refinement,
                       niter=niter,
                       minimizer_kwargs={'args': args, "method": params.method,
                           "jac": jac,'options':{'maxiter':params.maxiter, 'disp': params.disp},
                                         'hess': params.hess},
                       T=params.temp,
                       callback=None,
                       disp=not params.quiet and COMM.rank==0,
                       stepsize=params.stepsize)

    # save the final value for x
    #target.all_x.append(out.x)
    #target.save_x(optimized=True)
    target.x0[vary] = out.x
    return target.x0

class SimParams:
    def __init__(self):
        self.ucell = None
        self.RotXYZ = None
        self.Scale = None
        self.shift = None
        self.Nabc = None
        self.center = None # fdp edge center
        self.slope = None # fdp edge slope

#@profile
def model(x, SIM, Modeler, compute_grad=True, sanity_test=None):

    pfs = Modeler.pan_fast_slow
    PAR = Modeler.PAR

    # update the simlator crystal model
    SIM.D.Umatrix = PAR.Umatrix
    SIM.D.Bmatrix = PAR.Bmatrix
    #SIM.D.set_ncells_values(PAR.Nabc)
    SIM.D.shift_origin_z(SIM.detector, Modeler.shiftZ_meters)  # set the originZ as SIM.detector.origin plus the shiftZ which should be in meters
    # detector shift code looks like
    # for (int pid=0; pid< detector.size(); pid++)
    #  pix0_vectors[pid*3 + 2] = detector[pid].get_origin()[2]/1000.0 + shiftZ;
    #print(PAR.Nabc)
    #print(PAR.Umatrix)
    #print(PAR.Bmatrix)

    # update the energy spectrum
    #SIM.beam.spectrum = Modeler.spectrum
    #SIM.D.xray_beams = SIM.beam.xray_beams
    #print("Simulating %d energy channels" % len(SIM.D.xray_beams))
    SIM.D.xray_beams = Modeler.xray_beams #SIM.beam.xray_beams
    #SIM.D.update_xray_beams(SIM.beam.xray_beams)

    # how many parameters we simulate
    npix = int(len(pfs) / 3)
    #print("Simulating %d pixels" % npix)
    nparam = len(x)
    grad = np.zeros(nparam)  # gradient

    # set the scale factor for this shot
    scale_reparam = x[0]  # scale is always first
    scale = PAR.Scale.get_val(scale_reparam)

    # set the Ncells for this shot  (nabc are always in x at positions 1,2,3)
    Nabc = [PAR.Nabc_params[i_N].get_val(x[1+i_N]) for i_N in range(3)]
    SIM.D.set_ncells_values(tuple(Nabc))

    # set the Umatrix perturbation for this shot (rotx,roty,rotz are in x at positions 4,5,6)
    RotXYZ = [PAR.RotXYZ_params[i_rot].get_val(x[4+i_rot]) for i_rot in range(3)]
    for i_rot in range(3):
        SIM.D.set_value(hopper_utils.ROTXYZ_IDS[i_rot], RotXYZ[i_rot])

    if sanity_test==0:  # test hkl variation within each shoebox
        SIM.D.track_Fhkl = True
        PFS = np.reshape(pfs, (npix, 3))
        uroi = set(Modeler.roi_id)
        all_good_count_stats = []
        all_bad_count_stats = []
        for ii,i_roi in enumerate(uroi):
            output = Redirect(stdout=True)
            i_roi_sel = Modeler.roi_id==i_roi
            with output:
                sel = np.logical_and(i_roi_sel, Modeler.all_trusted)
                pfs_roi = PFS[sel]
                pfs_roi = np.ascontiguousarray(pfs_roi.ravel())
                pfs_roi = flex.size_t(pfs_roi)
                SIM.D.add_diffBragg_spots(pfs_roi)
            i_fcell = Modeler.all_fcell_global_idx[i_roi_sel][0]
            shoebox_hkl = SIM.asu_from_i_fcell[i_fcell]
            lines = output.stdout.split("\n")
            count_stats = {}
            for l in lines:
                if l.startswith("Pixel"):
                    hkl = l.split()[5]
                    hkl = tuple(map(int, hkl.split(",")))
                    hkl = utils.map_hkl_list([hkl], True, SIM.space_group_symbol)[0]
                    count = int(l.split()[8])
                    if hkl in count_stats:
                        count_stats[hkl] += count
                    else:
                        count_stats[hkl] = count
            ntot = sum(count_stats.values())
            assert shoebox_hkl in count_stats
            print("Shoebox hkl", shoebox_hkl)
            for hkl in count_stats:
                frac = count_stats[hkl] / float(ntot)
                h,k, l = hkl
                print("\tstep hkl %d,%d,%d : frac=%.1f%%" % (h,k,l,frac*100))
                count_stats[hkl] = frac
            if len(count_stats)==1:
                all_good_count_stats.append([shoebox_hkl,count_stats])
            else:
                all_bad_count_stats .append([shoebox_hkl, count_stats])
        if all_bad_count_stats:
            print("Shot %s had %d /  %d rois with HKL variation" %(Modeler.exp_name, len(all_bad_count_stats), len(uroi)))
            percs = [stats[sb_hkl]*100 for sb_hkl, stats in all_bad_count_stats]
            ave_perc = sum(percs) / len(percs)
            min_perc = min(percs)
            nmax = max(len(stats) for _,stats in all_bad_count_stats)
            print("\tMin %.1f%%, Mean=%.1f%%, most variation: %d hkls in a shoebox" %(min_perc, ave_perc, nmax))

            for sb_h,stats in all_bad_count_stats:
                h,k,l = zip(*stats.keys())
                if len(set(h))>1 or len(set(k))>1:
                    print("Weird HK vary: shot %s" % Modeler.exp_name, sb_h)
                if not np.all(np.sort(l) == np.arange(min(l), min(l)+len(l))):
                    print("Weird L sort: shot %s" % Modeler.exp_name, sb_h)
                if len(stats) > 3:
                    print("Weird Nmax: shot %s" % Modeler.exp_name, sb_h)

        return  # end sanity test on hkl variation

    # compute the forward model, and gradients when instructed
    SIM.D.add_diffBragg_spots(pfs, Modeler.all_nominal_hkl)
    bragg_no_scale = SIM.D.raw_pixels_roi[:npix].as_numpy_array()
    bragg = scale*bragg_no_scale

    model_pix = bragg + Modeler.all_background
    resid = (Modeler.all_data - model_pix)
    resid_square = resid ** 2
    V = model_pix + Modeler.sigma_rdout ** 2

    if compute_grad:
        common_grad_term = (0.5 / V * (1 - 2 * resid - resid_square / V))

        # compute the scale factor gradient term, which is related directly to the forward model
        if PAR.Scale.refine:
            scale_grad = bragg_no_scale
            scale_grad = PAR.Scale.get_deriv(scale_reparam, scale_grad)
            grad[0] += (common_grad_term * scale_grad)[Modeler.all_trusted].sum()

        # compute Ncells abc gradient terms
        if PAR.Nabc_params[0].refine:
            Nabc_grad = SIM.D.get_ncells_derivative_pixels()
            for i_N in range(3):
                N_grad = scale * (Nabc_grad[i_N][:npix].as_numpy_array())
                N_grad = PAR.Nabc_params[i_N].get_deriv(x[1+i_N], N_grad)
                grad[1+i_N] += (common_grad_term * N_grad)[Modeler.all_trusted].sum()

        # compute the RotXYZ gradient terms
        if PAR.RotXYZ_params[0].refine:
            for i_rot in range(3):
                rot_grad = scale*SIM.D.get_derivative_pixels(hopper_utils.ROTXYZ_IDS[i_rot]).as_numpy_array()[:npix]
                rot_grad = PAR.RotXYZ_params[i_rot].get_deriv(x[4+i_rot], rot_grad)
                grad[4+i_rot] += (common_grad_term * rot_grad)[Modeler.all_trusted].sum()

        # TODO add a dimension to get_derivative_pixels(FHKL_ID), in case that pixels hold information on multiple HKL
        fcell_grad = SIM.D.get_derivative_pixels(FHKL_ID)
        fcell_grad = scale * (fcell_grad[:npix].as_numpy_array())
        for i_fcell in Modeler.unique_i_fcell:
            sel = Modeler.is_i_fcell[i_fcell]
            this_fcell_grad = fcell_grad[sel]
            rescaled_amplitude = x[-SIM.n_global_fcell+i_fcell]
            #TODO try numexpr to speed up the trig operations in get_deriv
            this_fcell_grad = SIM.Fhkl_modelers[i_fcell].get_deriv(rescaled_amplitude, this_fcell_grad)

            this_common_g = common_grad_term[sel]
            this_trusted = Modeler.all_trusted[sel]
            grad[-SIM.n_global_fcell+i_fcell] += (this_common_g*this_fcell_grad)[this_trusted].sum()
    resid_term = (.5*(np.log(2*np.pi*V) + resid_square / V))[Modeler.all_trusted].sum()   # negative log Likelihood target

    ## sanity check on amplitudes:
    if sanity_test==1:
        bad_rois = []
        for i in set(Modeler.roi_id):
            bragg_i = bragg[Modeler.roi_id == i]
            if np.all(bragg_i == 0):
                bad_rois.append(i)
            else:
                i_fcell = Modeler.all_fcell_global_idx[Modeler.roi_id == i][0]
                if np.isnan(grad[-SIM.n_global_fcell+i_fcell]):
                    bad_rois.append(i)

        if bad_rois:
            #raise ValueError("Oops")
            for b in bad_rois:
                i_fcell = Modeler.all_fcell_global_idx[Modeler.roi_id == b][0]
                hkl_asu = SIM.asu_from_i_fcell[i_fcell]
                try:
                    bad_amp = SIM.crystal.miller_array.value_at_index(hkl_asu)
                except AssertionError:
                    bad_amp = np.nan
                value_of_grad = grad[-SIM.n_global_fcell + i_fcell]
                h,k,l = hkl_asu
                print("No signal in %s; Bragg peak at %d %d %d: Famp=%f, grad=%f" %(Modeler.exp_name, h,k,l,bad_amp, value_of_grad))
            print("NANs above indicate an hkl that is not in the starting structure factor array")
            print("Others indicate that the model has no signal in the shoebox")
            print("Remove the bad ROIS from the input, exiting.")

    #if np.any(np.isnan(grad)):
    #    bad_i = np.where(np.isnan(grad))[0]
    #    scale_is_nan = False
    #    bad_amps = ""
    #    for i in bad_i:
    #        if i==0:
    #            scale_is_nan = True
    #            continue
    #        i_fcell = i - 1
    #        rescaled_amplitude = x[i_fcell]
    #        #TODO try numexpr to speed up the trig operations in get_deriv
    #        amp = SIM.Fhkl_modelers[i_fcell].get_val(rescaled_amplitude)
    #        bad_amps += "%f, " %amp
    #    if scale_is_nan:
    #        bad_amps = "Scale, "+ bad_amps
    #    raise ValueError(bad_amps)

    #    #for i in bad_i:
    #    #    if i==0:
    #    #        raise ValueError("Scale factor gradient is Nan!")


    grad_is_nan = np.isnan(grad)
    assert not np.any(grad_is_nan) , "%s, %s" % \
        (Modeler.exp_name, ",".join( [str(x) for x in np.where(grad_is_nan)[0]] ) )
    #    bad_i = np.where(grad_is_nan)[0]
    #    for i in bad_i:
    #    print("")
    #    print(Modeler.exp_name)


    return resid_term, grad, model_pix


class TargetFunc:
    def __init__(self, params, SIM):
        self.all_x = []
        self.params = params
        self.SIM = SIM
        self.save_count = 0
        self.num_minimum = 0
        self.f_evals = []

        self.vary = None
        self.x0 = None

    def __call__(self, x, *args, **kwargs):
        self.x0[self.vary] = x

        self.all_x.append(self.x0)
        if len(self.all_x) == self.params.x_write_freq:
            self.save_x()
        f, self.g = target_func(self.x0, *args, **kwargs)
        self.f_evals.append(f)
        return f

    def jac(self, x, *args):
        if self.g is not None:
            return self.g[self.vary]

    def save_x(self, optimized=False):
        xdir = os.path.join(self.params.outdir, "x")
        mpi_safe_makedirs(xdir)
        if COMM.rank == 0:
            outpath = os.path.join(xdir, "x_info_hop%d_%d.npz" % (self.num_minimum,self.save_count))
            if optimized:
                outpath = os.path.join(xdir, "x_info_hop%d_final.npz" % self.num_minimum)
            Fidx, Fdata = update_Fhkl(self.SIM, self.all_x[-1])
            # TODO save scale factors
            np.savez(outpath, f_evals=self.f_evals, x=self.all_x[-1], Fidx=Fidx, Fdata=Fdata, Nhkl=self.SIM.n_global_fcell)
        COMM.barrier()
        self.all_x = []
        self.save_count += 1

    def at_minimum(self, x, f, accept):
        # NOTE doesnt seem to hit if niter=0  .. or something else weird happening here .. maybe just pass this meth
        self.x0[self.vary] = x
        self.all_x.append(self.x0)
        self.save_x(optimized=True)
        self.save_count = 0
        self.num_minimum += 1
        self.all_x = []


#@profile
def target_func(x, rank_xidx, SIM, Modelers, verbose=True, params=None, compute_grad=True, ntimes=None, save=None):
    verbose = verbose and COMM.rank==0
    t_start = time.time()
    timestamps = list(Modelers.keys())
    fchi = fG = fNabc = f_RotXYZ = 0
    g = np.zeros_like(x)

    #   do a global update of the Fhkl parameters in the simulator object
    t_update = time.time()
    SIM.update_Fhkl(SIM, x)
    t_update = time.time()-t_update

    all_t_model = 0
    for t in timestamps:
        # data modeler for this expt
        Mod_t = Modelers[t]

        #  x_t should be [G,Fhkl0, Fhkl1, Fhkl2, ..]
        x_t = x[rank_xidx[t]]

        t_model = time.time()

        if params.sanity_test_hkl_variation:
            sanity_test = 0
        elif params.sanity_test_amplitudes:
            sanity_test = 1
        else:
            sanity_test=None

        resid_term, grad, _ = model(x_t, SIM, Mod_t, compute_grad=compute_grad, sanity_test=sanity_test)
        t_model  = time.time()-t_model
        all_t_model += t_model

        # update global target functional
        fchi += resid_term  #(.5*(np.log(2*np.pi*V) + resid_square / V))[trusted_t].sum()   # negative log Likelihood target

        # if a shot is divided across ranks, then the restraint term
        # for that shot will be computed n times, hence we need to reduce
        # restraint terms by that factor
        nn = 1. / ntimes[t]  # n times is the number of ranks which are modeling part of a shot, id imagine its usually a small number like 1 or 2

        # scale factor "G" restraint
        if Mod_t.PAR.Scale.refine:
            G_rescaled = x_t[0]
            G = Mod_t.PAR.Scale.get_val(G_rescaled)
            delG = params.centers.G - G
            G_V = params.betas.G
            fG += nn*.5* delG**2/G_V

        if Mod_t.PAR.Nabc_params[0].refine:
            delNabc = []  # cache delta N for restraints gradient term below
            for i_N in range(3):
                N_V = params.betas.Nabc[i_N]
                N_current = Mod_t.PAR.Nabc_params[i_N].get_val(x_t[1+i_N])  # N always ranges 1-3 in the per-shot param list x_t
                del_N =  params.centers.Nabc[i_N] - N_current
                delNabc.append(del_N)
                fNabc += nn *.5 * del_N**2 / N_V

        if Mod_t.PAR.RotXYZ_params[0].refine:
            delRotXYZ = []
            for i_rot in range(3):
                rot_V = params.betas.RotXYZ # TODO this beta a list ?
                rot_current = Mod_t.PAR.RotXYZ_params[i_rot].get_val(x_t[4+i_rot])
                del_rot = params.centers.RotXYZ[i_rot] - rot_current
                delRotXYZ.append(del_rot)
                f_RotXYZ += nn*.5*del_rot**2 / rot_V

        if compute_grad:
            # copy the per-shot contributions of the gradients to the global gradient array

            # scale gradient term
            if Mod_t.PAR.Scale.refine:
                scale_global_idx = rank_xidx[t][0]
                g[scale_global_idx] += grad[0] # model
                g[scale_global_idx] += nn*Mod_t.PAR.Scale.get_deriv(G_rescaled, -delG / G_V) # restraint

            # Nabc gradient term
            if Mod_t.PAR.Nabc_params[0].refine:
                for i_N in range(3):
                    # 1+i_N is the per-shot position for Nabc
                    Nabc_global_idx = rank_xidx[t][1+i_N]
                    g[Nabc_global_idx] += grad[1+i_N]  # model
                    N_var = params.betas.Nabc[i_N]
                    del_N = delNabc[i_N]
                    g[Nabc_global_idx] += \
                        Mod_t.PAR.Nabc_params[i_N].get_deriv(x_t[1+i_N], -nn*.5*del_N / N_var) # restraint

            # RotXYZ gradient term
            if Mod_t.PAR.RotXYZ_params[0].refine:
                for i_rot in range(3):
                    # 4+i_rot is the per-shot position for Rot_i
                    RotXYZ_global_idx = rank_xidx[t][4+i_rot]
                    g[RotXYZ_global_idx] += grad[4+i_rot]  # model
                    rot_var = params.betas.RotXYZ  #TODO make this beta a list?
                    del_rot = delRotXYZ[i_rot]
                    g[RotXYZ_global_idx] += \
                        Mod_t.PAR.RotXYZ_params[i_rot].get_deriv(x_t[4+i_rot], -nn*.5*del_rot / rot_var) # restraint

            # Fhkl term updates
            g[-SIM.n_global_fcell:] += grad[-SIM.n_global_fcell:]

    # bring in data from all ranks
    if params.sanity_test_amplitudes:
        print("Sanity checked amplitudes, exiting")
        COMM.barrier()
        exit()
    t_mpi_start = time.time()  # This time should have large variance amongst ranks, as some ranks finish ealy and have to wait
    fchi = COMM.bcast(COMM.reduce(fchi))
    fG = COMM.bcast(COMM.reduce(fG))
    fNabc = COMM.bcast(COMM.reduce(fNabc))
    f_RotXYZ = COMM.bcast(COMM.reduce(f_RotXYZ))
    g = COMM.bcast(COMM.reduce(g))
    t_mpi_done = time.time()

    # add the Fhkl restraints
    Fhkl_current = np.array([ \
        SIM.Fhkl_modelers[i_fcell].get_val(x[-SIM.n_global_fcell+i_fcell]) \
        for i_fcell in range(SIM.n_global_fcell)])

    if params.use_wilson_restraints:
        # TODO vectorize or MPI this loop
        f_Fhkl = 0
        for i_fcell in range(SIM.n_global_fcell):
            hkl = SIM.asu_from_i_fcell[i_fcell]
            L = SIM.Fsq_ref_for_asu[hkl]
            F = Fhkl_current[i_fcell]
            if SIM.asu_is_centric[hkl]:
                f_Fhkl += F**2 / L/2
            else:
                #f_Fhkl += -np.log(2 * F / L) + F**2 / L
                f_Fhkl += F**2 / L - np.log(F) 

    else:
        Fhkl_init = np.array([SIM.Fhkl_modelers[i_fcell].init for i_fcell in range(SIM.n_global_fcell)])
        delta_F = Fhkl_init - Fhkl_current
        if params.sigma_frac is None:
            var_F = np.sum(Fhkl_init**2)
        else:
            var_F = (params.sigma_frac*Fhkl_init)**2
        f_Fhkl = np.sum(delta_F**2 / var_F)

    if compute_grad:
        # update the gradient restraint term for structure factor amplitudes
        Fhkl_rescaled = x[-SIM.n_global_fcell:]
        if params.use_wilson_restraints:
            # TODO vectorize or MPI this loop
            for i_fcell in range(SIM.n_global_fcell):
                hkl = SIM.asu_from_i_fcell[i_fcell]
                L = SIM.Fsq_ref_for_asu[hkl]
                F = Fhkl_current[i_fcell]
                if SIM.asu_is_centric[hkl]:
                    g_term = F/L
                else:
                    g_term = 2*F/L - 1/F
                g[-SIM.n_global_fcell+i_fcell] = SIM.Fhkl_modelers[i_fcell].get_deriv(Fhkl_rescaled[i_fcell], g_term)
        else:
            Fhkl_restraint_grad = -2*delta_F / var_F
            g[-SIM.n_global_fcell:] += np.array([ \
                SIM.Fhkl_modelers[i_fcell].get_deriv(Fhkl_rescaled[i_fcell], Fhkl_restraint_grad[i_fcell]) \
                for i_fcell in range(SIM.n_global_fcell)])

    f = fchi + fG + f_Fhkl + fNabc + f_RotXYZ
    chi = fchi / f *100
    gg = fG / f*100
    ff = f_Fhkl / f *100
    nn = fNabc / f * 100
    rr = f_RotXYZ / f * 100
    gnorm = np.linalg.norm(g)

    t_done = time.time()
    t_mpi = t_mpi_done - t_mpi_start
    t_total = t_done - t_start

    frac_mpi = t_mpi / t_total *100.
    frac_model = all_t_model / t_total * 100.
    frac_update = t_update / t_total * 100.

    if verbose:
        print("F=%10.7g (chi: %.1f%%, G: %.1f%%, N: %.1f%%, Rot: %.1f%%, Fhkl: %.1f%%); |g|=%10.7e; Total iter time=%.1f millisec (mpi: %.1f%% , model: %.1f%%, updateFhkl: %.1f%%)" \
              % (f, chi, gg, nn, rr, ff, gnorm, t_total*1000, frac_mpi, frac_model, frac_update))
    return f, g


def save_model_Z(img_path, Modeler, model_pix):
    pfs = Modeler.pan_fast_slow
    pids = pfs[0::3]
    xs = pfs[1::3]
    ys = pfs[2::3]

    sigma = np.sqrt(Modeler.all_data + Modeler.sigma_rdout ** 2)
    sigma2 = np.sqrt(model_pix + Modeler.sigma_rdout ** 2)
    Zdiff = model_pix - Modeler.all_data
    Z = Zdiff / sigma
    Z2 = Zdiff / sigma2
    with h5py.File(img_path, "w") as h5:
        comp = {"compression": "lzf"}
        h5.create_dataset("Z_data_noise", data=Z, **comp)
        h5.create_dataset("Z_model_noise", data=Z2, **comp)
        h5.create_dataset("pids", data=pids, **comp)
        h5.create_dataset("ys", data=ys, **comp)
        h5.create_dataset("xs", data=xs, **comp)
        h5.create_dataset("refls_idx", data=Modeler.all_refls_idx, **comp)
        h5.create_dataset("trusted", data=Modeler.all_trusted, **comp)

        h5.create_dataset("exp_name", data=Modeler.exp_name)
        h5.create_dataset("ref_name", data=Modeler.ref_name)



class Fhkl_updater:

    def __init__(self, SIM, Modelers):
        self.unique_hkl = set()
        for i_exp in Modelers:
            Hi_asu_in_exp = Modelers[i_exp].Hi_asu
            self.unique_hkl = self.unique_hkl.union(set(Hi_asu_in_exp))
        #n_unique = len(self.unique_hkl)

        # occastionally shoeboxes model neighboring spots , so we need to update those as well
        for h, k, l in self.unique_hkl.copy():
            for dl in [-2,0,2]:
                for dk in [-2,0,2]:
                    for dh in [-2,0,2]:
                        hkl_shifted = h+dh, k+dk, l+dl
                        if hkl_shifted in SIM.i_fcell_from_asu:
                            # add the hkls that are being updated potentially on other ranks modeling other shots
                            self.unique_hkl.add(hkl_shifted)
        #if len(self.unique_hkl)> n_unique:
        #    print("rank %d added neighboring spot Fhkl to updater" % COMM.rank)

        self.equiv_hkls = {}
        self.update_idx = flex.miller_index()
        for hkl_asu in self.unique_hkl:
            equivs = [i.h() for i in miller.sym_equiv_indices(SIM.space_group, hkl_asu).indices()]
            self.equiv_hkls[hkl_asu] = equivs
            for hkl_equiv in equivs:
                self.update_idx.append(hkl_equiv)

    #@profile
    def __call__(self, SIM, x):
        #idx, data = SIM.D.Fhkl_tuple
        update_amps = []
        for hkl_asu in self.unique_hkl:
            i_fcell = SIM.i_fcell_from_asu[hkl_asu]
            rescaled_amplitude = x[-SIM.n_global_fcell+i_fcell]
            amp = SIM.Fhkl_modelers[i_fcell].get_val(rescaled_amplitude)
            update_amps += [amp]*len(self.equiv_hkls[hkl_asu])

        update_amps = flex.double(update_amps)
        SIM.D.quick_Fhkl_update((self.update_idx, update_amps))


def update_Fhkl(SIM, x):
    # NOTE, this is the slow version use Fhkl_updater for iterations during ensemble refinement
    #Fidx, Fdata = SIM.D.Fhkl_tuple

    Fidx, Fdata = [],[]
    for i_fcell in range(SIM.n_global_fcell):
        # get the asu miller index
        hkl_asu = SIM.asu_from_i_fcell[i_fcell]

        # get the current amplitude
        xval = x[-SIM.n_global_fcell+i_fcell]
        new_amplitude = SIM.Fhkl_modelers[i_fcell].get_val(xval)

        # now surgically update the p1 array in nanoBragg with the new amplitudes
        # (need to update each symmetry equivalent)
        Fidx.append(hkl_asu)
        Fdata.append(new_amplitude)
        #equivs = [i.h() for i in miller.sym_equiv_indices(SIM.space_group, hkl_asu).indices()]
        #for h_equiv in equivs:
        #    p1_idx = SIM.idx_from_p1[h_equiv]
        #    Fdata[p1_idx] = new_amplitude  # set the data with the new value
    return Fidx,Fdata


def setup_Fhkl_attributes(SIM, params, Modelers):
    SIM.space_group_symbol = params.space_group
    SIM.space_group = sgtbx.space_group(sgtbx.space_group_info(symbol=params.space_group).type().hall_symbol())
    a,b = aggregate_Hi(Modelers)
    SIM.i_fcell_from_asu = a
    SIM.asu_from_i_fcell = b
    SIM.n_global_fcell = len(SIM.i_fcell_from_asu)

    # get the nanoBragg internal p1 positional index from the asu miller-index
    SIM.Fidx, SIM.Fdata = SIM.D.Fhkl_tuple
    SIM.idx_from_p1 = {h: i for i, h in enumerate(SIM.Fidx)}

    # get the initial amplitude value from the asu miller-index
    asu_hi = [SIM.asu_from_i_fcell[i_fcell] for i_fcell in range(SIM.n_global_fcell)]
    SIM.fcell_init_from_asu = {h: SIM.Fdata[SIM.idx_from_p1[h]] for h in asu_hi}

    SIM.Fhkl_modelers = []
    for i_fcell in range(SIM.n_global_fcell):
        p = RangedParameter()
        p.sigma = params.sigmas.Fhkl
        p.maxval = params.maxs.Fhkl
        p.minval = params.mins.Fhkl
        asu = SIM.asu_from_i_fcell[i_fcell]
        p.init = SIM.fcell_init_from_asu[asu]
        SIM.Fhkl_modelers.append(p)

    # sanity test, passes
    #for i in range(SIM.n_global_fcell):
    #    val1 = SIM.Fhkl_modelers[i].get_val(1)
    #    h = SIM.asu_from_i_fcell[i]
    #    val2 = SIM.crystal.miller_array.value_at_index(h)
    #    assert np.allclose(val1, val2)
    # are we using a reference for restraints?
    if params.scaling_reference_mtz_name is not None:

        Fref = utils.open_mtz(params.scaling_reference_mtz_name, params.scaling_reference_mtz_col)
        if not Fref.is_xray_amplitude_array() or np.min(Fref.data()) < 0:
            Fref = Fref.as_amplitude_array()
        dspacing_ref, Fsq_ref = utils.get_ave_FF(Fref)

        temp_mset = Fref.miller_set(flex.miller_index(asu_hi), True)
        temp_data = flex.double([SIM.fcell_init_from_asu[h] for h in asu_hi])
        temp_mill_arr = miller.array(temp_mset, temp_data)
        # this scale factor is applied to Fref to bring it on the same scale as the data
        Fref_scale_factor = utils.compute_scale_to_minmize_r_factor(temp_mill_arr, Fref)
        Fsq_ref = Fsq_ref*Fref_scale_factor**2
        Fsq_ref = Fsq_ref*params.wilson_fat  # does this have any effect at all ?

        dspace_arr = temp_mill_arr.d_spacings()
        d_spacing_map = {h: d for h,d in zip(dspace_arr.indices(), dspace_arr.data())}
        dspacing_at_asu = {h: d_spacing_map[h]  for h in asu_hi}

        # this is variable "L" in my notes
        Fsq_endpt = Fsq_ref[0], Fsq_ref[-1]
        Fsq_ref_at_d = interp1d(dspacing_ref, Fsq_ref, bounds_error=False, fill_value=Fsq_endpt)
        SIM.Fsq_ref_for_asu = {h: Fsq_ref_at_d(d) for h, d in dspacing_at_asu.items()}

        # we will also need the centric flags as well for determining the proper probability distribution
        is_centric = temp_mill_arr.centric_flags()
        is_centric_map = {h: c for h,c in zip(is_centric.indices(), is_centric.data())}
        SIM.asu_is_centric = {h: is_centric_map[h] for h in asu_hi}


def aggregate_Hi(Modelers):
    # aggregate all miller indices
    Hi_asu_all_ranks = []
    for i_exp in Modelers:
        Hi_asu_all_ranks += Modelers[i_exp].Hi_asu
    Hi_asu_all_ranks = COMM.reduce(Hi_asu_all_ranks)
    Hi_asu_all_ranks = COMM.bcast(Hi_asu_all_ranks)

    # this will map the measured miller indices to their index in the LBFGS parameter array self.x
    idx_from_asu = {h: i for i, h in enumerate(set(Hi_asu_all_ranks))}
    # we will need the inverse map during refinement to update the miller array in diffBragg, so we cache it here
    asu_from_idx = {i: h for i, h in enumerate(set(Hi_asu_all_ranks))}

    return idx_from_asu, asu_from_idx


def mpi_safe_makedirs(dname):
    if COMM.rank == 0:
        utils.safe_makedirs(dname)
    COMM.barrier()


def determine_per_rank_max_num_pix(Modelers):
    max_npix = 0
    for i_exp in Modelers:
        rois = Modelers[i_exp].rois
        x1, x2, y1, y2 = map(np.array, zip(*rois))
        npix = np.sum((x2 - x1) * (y2 - y1))
        max_npix = max(npix, max_npix)
        print("Rank %d, shot %d has %d pixels" % (COMM.rank, i_exp + 1, npix))
    print("Rank %d, max pix to be modeled: %d" % (COMM.rank, max_npix))
    return max_npix


if __name__ == '__main__':
    from dials.util import show_mail_on_error

    with show_mail_on_error():
        script = Script()
        script.run()
