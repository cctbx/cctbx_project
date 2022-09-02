from __future__ import absolute_import, division, print_function
import os
from dials.algorithms.shoebox import MaskCode
from copy import deepcopy
from dials.model.data import Shoebox
import numpy as np
from scipy.optimize import dual_annealing, basinhopping
from collections import Counter
from scitbx.matrix import sqr, col
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model import Spectrum
from xfel.util.jungfrau import get_pedestalRMS_from_jungfrau
from simtbx.nanoBragg.utils import downsample_spectrum
from dials.array_family import flex
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners.parameters import RangedParameter, Parameters
from simtbx.diffBragg.attr_list import NB_BEAM_ATTRS, NB_CRYST_ATTRS, DIFFBRAGG_ATTRS

try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = None


import logging
MAIN_LOGGER = logging.getLogger("diffBragg.main")
PROFILE_LOGGER = logging.getLogger("diffBragg.profile")

ROTX_ID = 0
ROTY_ID = 1
ROTZ_ID = 2
ROTXYZ_IDS = ROTX_ID, ROTY_ID, ROTZ_ID
NCELLS_ID = 9
UCELL_ID_OFFSET = 3
DETZ_ID = 10
FHKL_ID = 11
ETA_ID = 19
DIFFUSE_ID = 23
LAMBDA_IDS = 12, 13

DEG = 180 / np.pi


def write_SIM_logs(SIM, log=None, lam=None):
    """
    Logs properties of SIM.D (diffBragg instance), and SIM.crystal, SIM.beam (nanoBragg beam and crystal)
    These are important for reproducing results
    :param SIM: sim_data instance used during hopper refinement, member of the data modeler
    :param log: optional log file to dump attributes of SIM
    :param lam: optional lambda file to dump the spectra to disk, can be read usint diffBragg.utils.load_spectra_file
    """
    if log is not None:
        with open(log, "w") as o:
            print("<><><><><>", file=o)
            print("DIFFBRAGG", file=o)
            print("<><><><><>", file=o)
            for attr in DIFFBRAGG_ATTRS:
                val = getattr(SIM.D, attr)
                print(attr+": ", val, file=o)
            print("\n<><><>", file=o)
            print("BEAM", file=o)
            print("<><><>", file=o)
            for attr in NB_BEAM_ATTRS:
                val = getattr(SIM.beam, attr)
                print(attr+": ", val, file=o)
            print("\n<><><><>", file=o)
            print("CRYSTAL", file=o)
            print("<><><><>", file=o)
            for attr in NB_CRYST_ATTRS:
                val = getattr(SIM.crystal, attr)
                print(attr+": ", val, file=o)
    if lam is not None:
        wavelen, wt = zip(*SIM.beam.spectrum)
        utils.save_spectra_file(lam, wavelen, wt)


def free_SIM_mem(SIM):
    """
    Frees memory allocated to host CPU (and GPU device, if applicable).
    Using this method is critical for serial applications!
    :param SIM: sim_data instance used during hopper refinement, member of the data modeler
    """
    SIM.D.free_all()
    SIM.D.free_Fhkl2()
    try:
        SIM.D.gpu_free()
    except TypeError:
        pass  # occurs on CPU-only builds


def finalize_SIM(SIM, log=None, lam=None):
    """
    thin wrapper to free_SIM_mem and write_SIM_logs
    :param SIM: sim_data instance used during hopper refinement, member of the data modeler
    :param log: optional log file to dump attributes of SIM
    :param lam: optional lambda file to dump the spectra to disk, can be read usint diffBragg.utils.load_spectra_file
    """
    write_SIM_logs(SIM, log, lam)
    free_SIM_mem(SIM)


class DataModeler:

    """
    The data modeler stores information in two ways:
    1- lists whose length is the number of pixels being modeled
    2- lists whose length is the number of shoeboxes being modeled

    for example if one is modeling 3 shoeboxes whose dimensions are 10x10, then
    the objects below like self.all_data will have length 300, and other objects like self.selection_flags
    will have length 3
    """

    def __init__(self, params):
        """ params is a simtbx.diffBragg.hopper phil"""
        self.no_rlp_info = False  # whether rlps are stored in the refls table
        self.params = params  # phil params (see diffBragg/phil.py)
        self._abs_path_params()
        self.SIM = None  # simulator object (instance of nanoBragg.sim_data.SimData
        self.E = None  # placeholder for the dxtbx.model.Experiment instance
        self.pan_fast_slow =None  # (pid, fast, slow) per pixel
        self.all_background =None  # background model per pixel (photon units)
        self.roi_id =None  # region of interest ID per pixel
        self.u_id = None  # set of unique region of interest ids
        self.all_freq = None  # flag for the h,k,l frequency of the observed pixel
        self.best_model = None  # best model value at each pixel
        self.best_model_includes_background = False  # whether the best model includes the background scattering estimate
        self.all_data =None  # data at each pixel (photon units)
        self.all_sigma_rdout = None  # this is either a float or an array. if the phil param use_perpixel_dark_rms=True, then these are different per pixel, per shot
        self.all_gain = None  # gain value per pixel (used during diffBragg/refiners/stage_two_refiner)
        self.all_sigmas =None  # error model for each pixel (photon units)
        self.all_trusted =None  # trusted pixel flags (True is trusted and therefore used during refinement)
        self.npix_total =None  # total number of pixels
        self.all_fast =None  # fast-scan coordinate per pixel
        self.all_slow =None  # slow-scan coordinate per pixel
        self.all_pid = None  # panel id per pixel
        self.rois=None  # region of interest (per spot)
        self.pids=None  # panel id (per spot)
        self.tilt_abc=None  # background plane constants (per spot), a,b are fast,slow scan components, c is offset
        self.selection_flags=None  # whether the spot was selected for refinement (sometimes poorly conditioned spots are rejected)
        self.background=None  # background for entire image (same shape as the detector)
        self.tilt_cov = None  # covariance estimates from background fitting (not used)
        self.simple_weights = None  # not used
        self.refls_idx = None  # position of modeled spot in original refl array
        self.refls = None  # reflection table
        self.nominal_sigma_rdout = None   # the value of the readout noise in photon units
        self.exper_name = None  # optional name specifying where dxtbx.model.Experiment was loaded from
        self.refl_name = None  # optional name specifying where dials.array_family.flex.reflection_table refls were loaded from
        self.spec_name = None  # optional name specifying spectrum file(.lam)

        self.Hi = None  # miller index (P1)
        self.Hi_asu = None  # miller index (high symmetry)

        # which attributes to save when pickling a data modeler
        self.saves = ["all_data", "all_background", "all_trusted", "best_model", "nominal_sigma_rdout",
                      "rois", "pids", "tilt_abc", "selection_flags", "refls_idx", "pan_fast_slow",
                        "Hi", "Hi_asu", "roi_id", "params", "all_pid", "all_fast", "all_slow"]

    def _abs_path_params(self):
        """adds absolute path to certain params"""
        if self.params.simulator.structure_factors.mtz_name is not None:
            self.params.simulator.structure_factors.mtz_name = os.path.abspath(self.params.simulator.structure_factors.mtz_name)

    def __getstate__(self):
        # TODO cleanup/compress
        return {name: getattr(self, name) for name in self.saves}

    def __setstate__(self, state):
        for name in state:
            setattr(self, name, state[name])

    def set_slices(self, attr_name):
        """finds the boundaries for each attr in the 1-D array of per-shot data
        :params: attr_name, str

        For example, in a DataModeler, roi_id is set for every pixel,
        >> Modeler.roi_id returns
        0 0 0 0 1 1 2 2 2 2 3 3 3 3 4 4 4 4 4 ...

        A call Modeler.set_slices("roi_id") adds attributes
            roi_id_slices
            roi_id_unique
        where roi_id_slices is a dictionary whose keys are the unique roi ids and whose values
        are a list of slices corresponding to each connected regeion of the same roi ids.

        In the example,
        >> Modeler.roi_id_slices[0][0]
        slice(0,4,1)
        """
        vals = self.__getattribute__(attr_name)
        # find where the vals change
        splitter = np.where(np.diff(vals) != 0)[0] + 1
        npix = len(self.all_data)
        slices = [slice(V[0], V[-1] + 1, 1) for V in np.split(np.arange(npix), splitter)]
        vals_ids = [V[0] for V in np.split(np.array(vals), splitter)]
        vals_id_slices = {}
        for i_vals, slc in zip(vals_ids, slices):
            if i_vals not in vals_id_slices:
                vals_id_slices[i_vals] = [slc]
            else:
                vals_id_slices[i_vals].append(slc)
        unique_vals = set(vals)
        self.__setattr__("%s_unique" % attr_name, unique_vals)
        logging.debug("Modeler has data on %d unique vals from attribute %s" % (len(unique_vals) , attr_name))
        self.__setattr__("%s_slices" % attr_name, vals_id_slices)

    @property
    def sigma_rdout(self):
        print("WARNING ,this attribute will soon be deprecated, use nominal_sigma_rdout instead!")
        return self.nominal_sigma_rdout

    def clean_up(self):
        free_SIM_mem(self.SIM)

    def set_experiment(self, exp, load_imageset=True):
        if isinstance(exp, str):
            self.E = ExperimentListFactory.from_json_file(exp, load_imageset)[0]
        else:
            self.E = exp
        if self.params.opt_det is not None:
            opt_det_E = ExperimentListFactory.from_json_file(self.params.opt_det, False)[0]
            self.E.detector = opt_det_E.detector
            MAIN_LOGGER.info("Set the optimal detector from %s" % self.params.opt_det)

        if self.params.opt_beam is not None:
            opt_beam_E = ExperimentListFactory.from_json_file(self.params.opt_beam, False)[0]
            self.E.beam = opt_beam_E.beam
            MAIN_LOGGER.info("Set the optimal beam from %s" % self.params.opt_beam)

    def load_refls(self, ref):
        if isinstance(ref, str):
            refls = flex.reflection_table.from_file(ref)
        else:
            # assert is a reflection table. ..
            refls = ref
        return refls

    def is_duplicate_hkl(self, refls):
        nref = len(refls)
        is_duplicate = np.zeros(nref, bool)
        if len(set(refls['miller_index'])) < nref:
            hkls = refls['miller_index']
            dupe_hkl = {h for h, count in Counter(hkls).items() if count > 1}
            for i_ref in range(nref):
                hh = refls[i_ref]['miller_index']
                is_duplicate[i_ref] = hh in dupe_hkl

        return is_duplicate

    def GatherFromReflectionTable(self, exp, ref, sg_symbol=None):
        self.set_experiment(exp, load_imageset=False)
        self.refls = self.load_refls(ref)
        nref = len(self.refls)
        if nref ==0:
            return False
        self.refls_idx = list(range(nref))
        self.rois = [(x1, x2, y1, y2) for x1,x2,y1,y2,_,_ in self.refls["shoebox"].bounding_boxes()]
        self.pids = list(self.refls["panel"])

        npan = len(self.E.detector)
        nfast, nslow = self.E.detector[0].get_image_size()  # NOTE assumes all panels same shape
        img_data = np.zeros((npan, nslow, nfast))
        background = np.zeros_like(img_data)
        is_trusted = np.zeros((npan, nslow, nfast), bool)
        for i_ref in range(nref):
            ref = self.refls[i_ref]
            pid = ref['panel']
            x1, x2, y1, y2 = self.rois[i_ref]

            # these are the in-bounds limits (on the panel)
            x1_onPanel = max(x1,0)
            x2_onPanel = min(x2,nfast)
            y1_onPanel = max(y1,0)
            y2_onPanel = min(y2,nslow)

            xdim = x2_onPanel-x1_onPanel
            ydim = y2_onPanel-y1_onPanel

            sb = ref['shoebox']
            sb_ystart = y1_onPanel - y1
            sb_xstart = x1_onPanel - x1
            sb_sliceY = slice(sb_ystart, sb_ystart+ydim,1)
            sb_sliceX = slice(sb_xstart, sb_xstart+xdim,1)

            dat_sliceY = slice(y1_onPanel, y1_onPanel+ydim,1)
            dat_sliceX = slice(x1_onPanel, x1_onPanel+xdim,1)
            img_data[pid,   dat_sliceY, dat_sliceX] = sb.data.as_numpy_array()[0,sb_sliceY,sb_sliceX]
            sb_bkgrnd = sb.background.as_numpy_array()[0,sb_sliceY,sb_sliceX]
            background[pid, dat_sliceY, dat_sliceX] = sb_bkgrnd
            fg_code = MaskCode.Valid + MaskCode.Foreground  # 5
            bg_code = MaskCode.Valid + MaskCode.Background + MaskCode.BackgroundUsed  # 19
            mask = sb.mask.as_numpy_array()[0,sb_sliceY,sb_sliceX]
            if self.params.refiner.refldata_trusted=="allValid":
                sb_trust = mask > 0
            elif self.params.refiner.refldata_trusted=="fg":
                sb_trust = mask==fg_code
            else:
                sb_trust = np.logical_or(mask==fg_code, mask==bg_code)

            below_zero = sb_bkgrnd <= 0
            if np.any(below_zero):
                nbelow = np.sum(below_zero)
                ntot = sb_bkgrnd.size
                MAIN_LOGGER.debug("background <= zero in %d/%d pixels from shoebox %d! Marking those pixels as untrusted!" %  ( nbelow, ntot, i_ref ))
                sb_trust[below_zero] = False

            is_trusted[pid, dat_sliceY,dat_sliceX] = sb_trust

            self.rois[i_ref] = x1_onPanel, x2_onPanel, y1_onPanel, y2_onPanel


        if self.params.refiner.refldata_to_photons:
            MAIN_LOGGER.debug("Re-scaling reflection data to photon units: conversion factor=%f" % self.params.refiner.adu_per_photon)
            img_data /= self.params.refiner.adu_per_photon
            background /= self.params.refiner.adu_per_photon

        # can be used for Bfactor modeling
        self.Q = np.linalg.norm(self.refls["rlp"], axis=1)
        self.nominal_sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon

        self.Hi = list(self.refls["miller_index"])
        if sg_symbol is not None:
            self.Hi_asu = utils.map_hkl_list(self.Hi, True, sg_symbol)
        else:
            self.Hi_asu = self.Hi

        self.data_to_one_dim(img_data, is_trusted, background)
        return True

    def GatherFromExperiment(self, exp, ref, remove_duplicate_hkl=True, sg_symbol=None):
        self.set_experiment(exp, load_imageset=True)

        refls = self.load_refls(ref)
        if len(refls)==0:
            MAIN_LOGGER.warning("no refls loaded!")
            return False

        if "rlp" not in list(refls[0].keys()):
            try:
                utils.add_rlp_column(refls, self.E)
                assert "rlp" in list(refls[0].keys())
            except KeyError:
                self.no_rlp_info = True
        img_data = utils.image_data_from_expt(self.E)
        img_data /= self.params.refiner.adu_per_photon
        is_trusted = np.ones(img_data.shape, bool)
        hotpix_mask = None
        if self.params.roi.hotpixel_mask is not None:
            is_trusted = utils.load_mask(self.params.roi.hotpixel_mask)
            hotpix_mask = ~is_trusted
        self.nominal_sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon

        roi_packet = utils.get_roi_background_and_selection_flags(
            refls, img_data, shoebox_sz=self.params.roi.shoebox_size,
            reject_edge_reflections=self.params.roi.reject_edge_reflections,
            reject_roi_with_hotpix=self.params.roi.reject_roi_with_hotpix,
            background_mask=None, hotpix_mask=hotpix_mask,
            bg_thresh=self.params.roi.background_threshold,
            use_robust_estimation=not self.params.roi.fit_tilt,
            set_negative_bg_to_zero=self.params.roi.force_negative_background_to_zero,
            pad_for_background_estimation=self.params.roi.pad_shoebox_for_background_estimation,
            sigma_rdout=self.nominal_sigma_rdout, deltaQ=self.params.roi.deltaQ, experiment=self.E,
            weighted_fit=self.params.roi.fit_tilt_using_weights,
            allow_overlaps=self.params.roi.allow_overlapping_spots,
            ret_cov=True, skip_roi_with_negative_bg=self.params.roi.skip_roi_with_negative_bg,
            only_high=self.params.roi.only_filter_zingers_above_mean, centroid=self.params.roi.centroid)

        if roi_packet is None:
            return False

        self.rois, self.pids, self.tilt_abc, self.selection_flags, self.background, self.tilt_cov = roi_packet

        if remove_duplicate_hkl and not self.no_rlp_info:
            is_not_a_duplicate = ~self.is_duplicate_hkl(refls)
            self.selection_flags = np.logical_and( self.selection_flags, is_not_a_duplicate)

        if self.params.refiner.res_ranges is not None:
            # TODO add res ranges support for GatherFromReflectionTable
            if self.no_rlp_info:
                raise NotImplementedError("Cannot set resolution limits when processing refls that are missing the RLP column")
            res_flags = np.zeros(len(refls)).astype(bool)
            res = 1. / np.linalg.norm(refls["rlp"], axis=1)
            for dmin,dmax in utils.parse_reso_string(self.params.refiner.res_ranges):
                MAIN_LOGGER.debug("Parsing res range %.3f - %.3f Angstrom" % (dmin, dmax))
                in_resShell = np.logical_and(res >= dmin, res < dmax)
                res_flags[in_resShell] = True

            MAIN_LOGGER.info("Resolution filter removed %d/%d refls outside of all resolution ranges " \
                              % (sum(~res_flags), len(refls)))
            self.selection_flags[~res_flags] = False

        if "miller_index" in list(refls.keys()):
            self.Hi = list(refls["miller_index"])
            if sg_symbol is not None:
                self.Hi_asu = utils.map_hkl_list(self.Hi, True, sg_symbol)
            else:
                self.Hi_asu = self.Hi

        if sum(self.selection_flags) == 0:
            MAIN_LOGGER.info("No pixels slected, continuing")
            return False
        self.refls = refls
        self.refls_idx = [i_roi for i_roi in range(len(refls)) if self.selection_flags[i_roi]]

        self.rois = [roi for i_roi, roi in enumerate(self.rois) if self.selection_flags[i_roi]]
        self.tilt_abc = [abc for i_roi, abc in enumerate(self.tilt_abc) if self.selection_flags[i_roi]]
        self.pids = [pid for i_roi, pid in enumerate(self.pids) if self.selection_flags[i_roi]]
        self.tilt_cov = [cov for i_roi, cov in enumerate(self.tilt_cov) if self.selection_flags[i_roi]]
        if not self.no_rlp_info:
            self.Q = [np.linalg.norm(refls[i_roi]["rlp"]) for i_roi in range(len(refls)) if self.selection_flags[i_roi]]
        refls = refls.select(flex.bool(self.selection_flags))


        self.data_to_one_dim(img_data, is_trusted, self.background)
        return True

    def data_to_one_dim(self, img_data, is_trusted, background):
        all_data = []
        all_sigma_rdout = []
        all_pid = []
        all_fast = []
        all_slow = []
        all_fast_relative = []
        all_slow_relative = []
        all_trusted = []
        all_sigmas = []
        all_background = []
        roi_id = []
        all_q_perpix = []
        all_refls_idx = []
        pixel_counter = np.zeros_like(img_data)
        self.all_nominal_hkl = []
        self.hi_asu_perpix = []
        numOutOfRange = 0
        perpixel_dark_rms = None
        if self.params.use_perpixel_dark_rms:
            perpixel_dark_rms = get_pedestalRMS_from_jungfrau(self.E)
            perpixel_dark_rms /= self.params.refiner.adu_per_photon

        for i_roi in range(len(self.rois)):
            pid = self.pids[i_roi]
            x1, x2, y1, y2 = self.rois[i_roi]
            Y, X = np.indices((y2 - y1, x2 - x1))
            data = img_data[pid, y1:y2, x1:x2].copy()
            pixel_counter[pid, y1:y2, x1:x2] += 1

            data = data.ravel()
            all_background += list(background[pid, y1:y2, x1:x2].ravel())
            trusted = is_trusted[pid, y1:y2, x1:x2].ravel()
            if perpixel_dark_rms is not None:
                sigma_rdout = perpixel_dark_rms[pid, y1:y2, x1:x2].ravel()
            else:
                sigma_rdout = np.ones_like(data)*self.nominal_sigma_rdout

            if self.params.roi.mask_outside_trusted_range:
                minDat, maxDat = self.E.detector[pid].get_trusted_range()
                data_out_of_range = np.logical_or(data <= minDat, data >= maxDat)
                if self.params.roi.mask_all_if_any_outside_trusted_range:
                    if np.any(data_out_of_range):
                        data_out_of_range[:] = True

                trusted[data_out_of_range] = False
                numOutOfRange +=np.sum(data_out_of_range)

            if self.params.mask_highest_values is not None:
                trusted[np.argsort(data)[-self.params.mask_highest_values:]] = False

            # TODO implement per-shot masking here
            #lower_cut = np.percentile(data, 20)
            #trusted[data < lower_cut] = False

            #d_strong_order = np.argsort(data)
            #trusted[d_strong_order[-1:]] = False
            all_trusted += list(trusted)
            #TODO ignore invalid value warning (handled below), or else mitigate it!
            all_sigmas += list(np.sqrt(data + sigma_rdout ** 2))
            all_sigma_rdout += list(sigma_rdout)
            all_fast += list(X.ravel() + x1)
            all_fast_relative += list(X.ravel())
            all_slow += list(Y.ravel() + y1)
            all_slow_relative += list(Y.ravel())
            all_data += list(data)
            npix = len(data)  # np.sum(trusted)
            all_pid += [pid] * npix
            roi_id += [i_roi] * npix
            all_refls_idx += [self.refls_idx[i_roi]] * npix
            if not self.no_rlp_info:
                all_q_perpix += [self.Q[i_roi]]*npix
            if self.Hi is not None:
                self.all_nominal_hkl += [tuple(self.Hi[i_roi])]*npix
                self.hi_asu_perpix += [self.Hi_asu[i_roi]] * npix

        if self.params.roi.mask_outside_trusted_range:
            MAIN_LOGGER.debug("Found %d pixels outside of trusted range" % numOutOfRange)
        all_freq = []
        for i_roi in range(len(self.rois)):
            pid = self.pids[i_roi]
            x1, x2, y1, y2 = self.rois[i_roi]
            freq = pixel_counter[pid, y1:y2, x1:x2].ravel()
            all_freq += list(freq)
        self.all_freq = np.array(all_freq)  # if no overlapping pixels, this should be an array of 1's
        if not self.params.roi.allow_overlapping_spots:
            if not np.all(self.all_freq==1):
                print(set(self.all_freq))
                raise ValueError("There are overlapping regions of interest, despite the command to not allow overlaps")

        self.all_q_perpix = np.array(all_q_perpix)
        pan_fast_slow = np.ascontiguousarray((np.vstack([all_pid, all_fast, all_slow]).T).ravel())
        self.pan_fast_slow = flex.size_t(pan_fast_slow)
        self.all_background = np.array(all_background)
        self.roi_id = np.array(roi_id)
        self.all_data = np.array(all_data)
        if np.allclose(all_sigma_rdout, self.nominal_sigma_rdout):
            self.all_sigma_rdout = self.nominal_sigma_rdout
        else:
            self.all_sigma_rdout = np.array(all_sigma_rdout)
        self.all_sigmas = np.array(all_sigmas)
        # note rare chance for sigmas to be nan if the args of sqrt is below 0
        self.all_trusted = np.logical_and(np.array(all_trusted), ~np.isnan(all_sigmas))

        if self.params.roi.skip_roi_with_negative_bg:
            # Dont include pixels whose background model is below 0
            self.all_trusted[self.all_background < 0] = False

        self.npix_total = len(all_data)
        self.all_fast = np.array(all_fast)
        self.all_slow = np.array(all_slow)
        self.all_pid = np.array(all_pid)
        #self.simple_weights = 1/self.all_sigmas**2
        self.u_id = set(self.roi_id)
        self.all_refls_idx = np.array(all_refls_idx)

        MAIN_LOGGER.debug("Modeler has %d/ %d trusted pixels" % (self.all_trusted.sum() , self.npix_total))

    def dump_gathered_to_refl(self, output_name, do_xyobs_sanity_check=False):
        """after running GatherFromExperiment, dump the gathered results
        (data, background etc) to a new reflection file which can then be used to run
        diffBragg without the raw data in the experiment (this exists mainly for portability, and
        unit tests)"""
        shoeboxes = []
        R = flex.reflection_table()
        for i_roi, i_ref in enumerate(self.refls_idx):
            roi_sel = self.roi_id==i_roi
            x1, x2, y1, y2 = self.rois[i_roi]
            roi_shape = y2-y1, x2-x1
            roi_img = self.all_data[roi_sel].reshape(roi_shape).astype(np.float32)  #NOTE this has already been converted to photon units
            roi_bg = self.all_background[roi_sel].reshape(roi_shape).astype(np.float32)

            sb = Shoebox((x1, x2, y1, y2, 0, 1))
            sb.allocate()
            sb.data = flex.float(np.ascontiguousarray(roi_img[None]))
            sb.background = flex.float(np.ascontiguousarray(roi_bg[None]))

            dials_mask = np.zeros(roi_img.shape).astype(np.int32)
            mask = self.all_trusted[roi_sel].reshape(roi_shape)
            dials_mask[mask] = dials_mask[mask] + MaskCode.Valid
            sb.mask = flex.int(np.ascontiguousarray(dials_mask[None]))

            # quick sanity test
            if do_xyobs_sanity_check:
                ref = self.refls[i_ref]
                x,y,_ = ref['xyzobs.px.value']
                assert x1 <= x <= x2, "exp %s; refl %d, %f %f %f" % (output_name, i_ref, x1,x,x2)
                assert y1 <= y <= y2, "exp %s; refl %d, %f %f %f" % (output_name, i_ref, y1,y,y2)

            R.extend(self.refls[i_ref: i_ref+1])
            shoeboxes.append(sb)

        R['shoebox'] = flex.shoebox(shoeboxes)
        R.as_file(output_name)

    def SimulatorFromExperiment(self, best=None):
        """
        This instantiates the simulator class and sets the model parameter objects, and attaches them to the
        simulator as an attribute called P.

        Optional best parameter is a single row of a pandas datafame containing the starting
        models, presumably optimized from a previous minimzation using this program
        """

        ParameterType = RangedParameter

        if best is not None:
            if self.params.number_of_xtals > 1:
                raise NotImplementedError("Not configured to load multiple crystals from pandas file")
            # set the crystal Umat (rotational displacement) and Bmat (unit cell)
            # Umatrix
            # NOTE: just set the best Amatrix here
            if self.params.apply_best_crystal_model:
                xax = col((-1, 0, 0))
                yax = col((0, -1, 0))
                zax = col((0, 0, -1))
                rotX,rotY,rotZ = best[["rotX", "rotY", "rotZ"]].values[0]
                RX = xax.axis_and_angle_as_r3_rotation_matrix(rotX, deg=False)
                RY = yax.axis_and_angle_as_r3_rotation_matrix(rotY, deg=False)
                RZ = zax.axis_and_angle_as_r3_rotation_matrix(rotZ, deg=False)
                M = RX * RY * RZ
                U = M * sqr(self.E.crystal.get_U())
                self.E.crystal.set_U(U)

                # Bmatrix:
                ucparam = best[["a","b","c","al","be","ga"]].values[0]
                ucman = utils.manager_from_params(ucparam)
                self.E.crystal.set_B(ucman.B_recipspace)

            ## TODO , currently need this anyway
            ucparam = best[["a","b","c","al","be","ga"]].values[0]
            ucman = utils.manager_from_params(ucparam)
            self.E.crystal.set_B(ucman.B_recipspace)

            # mosaic block
            self.params.init.Nabc = tuple(best.ncells.values[0])
            # scale factor
            self.params.init.G = best.spot_scales.values[0]

            if "detz_shift_mm" in list(best):
                self.params.init.detz_shift = best.detz_shift_mm.values[0]

            # TODO: set best eta_abc params

        # TODO: choose a phil param and remove support for the other: crystal.anositropic_mosaicity, or init.eta_abc
        MAIN_LOGGER.info("Setting initial anisotropic mosaicity from params.init.eta_abc: %f %f %f" % tuple(self.params.init.eta_abc))
        if self.params.simulator.crystal.has_isotropic_mosaicity:
            self.params.simulator.crystal.anisotropic_mosaicity = None
        else:
            self.params.simulator.crystal.anisotropic_mosaicity = self.params.init.eta_abc
        MAIN_LOGGER.info("Number of mosaic domains from params: %d" % self.params.simulator.crystal.num_mosaicity_samples)
        self.SIM = utils.simulator_from_expt_and_params(self.E, self.params)
        if self.SIM.D.mosaic_domains > 1:
            MAIN_LOGGER.info("Will use mosaic models: %d domains" % self.SIM.D.mosaic_domains)
        else:
            MAIN_LOGGER.info("Will not use mosaic models, as simulator.crystal.num_mosaicity_samples=1")

        if not self.params.fix.diffuse_gamma or not self.params.fix.diffuse_sigma:
            assert self.params.use_diffuse_models
        self.SIM.D.use_diffuse = self.params.use_diffuse_models
        if self.params.use_diffuse_models:
            if self.params.symmetrize_diffuse:
                assert self.params.space_group is not None
                self.SIM.D.laue_group_num = utils.get_laue_group_num(self.params.space_group)  # TODO this can also be retrieved from crystal model if params.space_group is None
                MAIN_LOGGER.debug("Set laue group number: %d (for diffuse models)" % self.SIM.D.laue_group_num)
            if self.params.diffuse_stencil_size > 0:
                self.SIM.D.stencil_size = self.params.diffuse_stencil_size
                MAIN_LOGGER.debug("Set diffuse stencil size: %d" % self.SIM.D.stencil_size)
        self.SIM.D.gamma_miller_units = self.params.gamma_miller_units
        self.SIM.isotropic_diffuse_gamma = self.params.isotropic.diffuse_gamma
        self.SIM.isotropic_diffuse_sigma = self.params.isotropic.diffuse_sigma

        if self.params.spectrum_from_imageset:
            downsamp_spec(self.SIM, self.params, self.E)
        elif self.params.gen_gauss_spec:
            set_gauss_spec(self.SIM, self.params, self.E)

        self.SIM.D.no_Nabc_scale = self.params.no_Nabc_scale  # TODO check gradients for this setting
        self.SIM.D.update_oversample_during_refinement = False
        self.SIM.num_xtals = self.params.number_of_xtals

        init = self.params.init
        sigma = self.params.sigmas
        mins = self.params.mins
        maxs = self.params.maxs
        centers = self.params.centers
        betas = self.params.betas
        if not self.params.use_restraints or self.params.fix.ucell:
            centers.ucell = [1,1,1,1,1,1]
            betas.ucell = [1,1,1,1,1,1]
        fix = self.params.fix
        P = Parameters()
        if self.params.init.random_Gs is not None:
            init.G = np.random.choice(self.params.init.random_Gs)
        for i_xtal in range(self.SIM.num_xtals):
            for ii in range(3):

                p = ParameterType(init=0, sigma=sigma.RotXYZ[ii],
                                  minval=mins.RotXYZ[ii], maxval=maxs.RotXYZ[ii],
                                  fix=fix.RotXYZ, name="RotXYZ%d_xtal%d" % (ii,i_xtal),
                                  center=centers.RotXYZ[ii], beta=betas.RotXYZ)
                P.add(p)

            p = RangedParameter(init=init.G + init.G*0.01*i_xtal, sigma=sigma.G,
                              minval=mins.G, maxval=maxs.G,
                              fix=fix.G, name="G_xtal%d" %i_xtal,
                              center=centers.G, beta=betas.G)
            P.add(p)


        # these parameters are equal for all texture-domains within a crystal
        fix_Nabc = [fix.Nabc]*3
        if self.params.simulator.crystal.has_isotropic_ncells:
            fix_Nabc = [fix_Nabc[0], True, True]

        fix_difsig = [fix.diffuse_sigma]*3
        if self.params.isotropic.diffuse_sigma:
            fix_difsig = [fix_difsig[0], True, True]

        fix_difgam = [fix.diffuse_gamma]*3
        if self.params.isotropic.diffuse_gamma:
            fix_difgam = [fix_difgam[0], True, True]

        fix_eta = [fix.eta_abc]*3
        if self.params.simulator.crystal.has_isotropic_mosaicity:
            fix_eta = [fix_eta[0], True, True]

        if self.params.init.random_Nabcs is not None:
            init.Nabc = np.random.choice(self.params.init.random_Nabcs, replace=True, size=3)
        for ii in range(3):
            # Mosaic domain tensor
            p = ParameterType(init=init.Nabc[ii], sigma=sigma.Nabc[ii],
                              minval=mins.Nabc[ii], maxval=maxs.Nabc[ii],
                              fix=fix_Nabc[ii], name="Nabc%d" % (ii,),
                              center=centers.Nabc[ii], beta=betas.Nabc[ii])
            P.add(p)

            # diffuse gamma and sigma
            p = ParameterType(init=init.diffuse_gamma[ii], sigma=sigma.diffuse_gamma[ii],
                              minval=mins.diffuse_gamma[ii], maxval=maxs.diffuse_gamma[ii],
                              fix=fix_difgam[ii], name="diffuse_gamma%d" % (ii,),
                              center=centers.diffuse_gamma[ii], beta=betas.diffuse_gamma[ii])
            P.add(p)

            p = ParameterType(init=init.diffuse_sigma[ii], sigma=sigma.diffuse_sigma[ii],
                              minval=mins.diffuse_sigma[ii], maxval=maxs.diffuse_sigma[ii],
                              fix=fix_difsig[ii], name="diffuse_sigma%d" % (ii,),
                              center=centers.diffuse_sigma[ii], beta=betas.diffuse_sigma[ii])
            P.add(p)

            # mosaic spread (mosaicity)
            p = ParameterType(init=init.eta_abc[ii], sigma=sigma.eta_abc[ii],
                              minval=mins.eta_abc[ii], maxval=maxs.eta_abc[ii],
                              fix=fix_eta[ii], name="eta_abc%d" % (ii,),
                              center=centers.eta_abc[ii], beta=betas.eta_abc[ii])
            P.add(p)

        ucell_man = utils.manager_from_crystal(self.E.crystal)
        ucell_vary_perc = self.params.ucell_edge_perc / 100.
        for i_uc, (name, val) in enumerate(zip(ucell_man.variable_names, ucell_man.variables)):
            if "Ang" in name:
                minval = val - ucell_vary_perc * val
                maxval = val + ucell_vary_perc * val
                if centers.ucell is not None:
                    cent = centers.ucell[i_uc]
                    beta = betas.ucell[i_uc]
                else:
                    if name == 'a_Ang':
                        cent = centers.ucell_a
                        beta = betas.ucell_a
                    elif name== 'b_Ang':
                        cent = centers.ucell_b
                        beta = betas.ucell_b
                    else:
                        cent = centers.ucell_c
                        beta = betas.ucell_c
                    assert cent is not None, "Set the center restraints properly!"
                    assert beta is not None
            else:
                val_in_deg = val * 180 / np.pi
                minval = (val_in_deg - self.params.ucell_ang_abs) * np.pi / 180.
                maxval = (val_in_deg + self.params.ucell_ang_abs) * np.pi / 180.
                if centers.ucell is not None:
                    cent = centers.ucell[i_uc]*np.pi / 180.
                    beta = betas.ucell[i_uc]
                else:
                    if name=='alpha_rad':
                        cent = centers.ucell_alpha
                        beta = betas.ucell_alpha
                    elif name=='beta_rad':
                        cent = centers.ucell_beta
                        beta = betas.ucell_beta
                    else:
                        cent = centers.ucell_gamma
                        beta = betas.ucell_gamma
                    assert cent is not None
                    assert beta is not None
                    cent = cent*np.pi / 180.

            p = ParameterType(init=val, sigma=sigma.ucell[i_uc],
                              minval=minval, maxval=maxval, fix=fix.ucell,
                              name="ucell%d" % (i_uc,),
                              center=cent,
                              beta=beta)
            MAIN_LOGGER.info(
                "Unit cell variable %s (currently=%f) is bounded by %f and %f" % (name, val, minval, maxval))
            P.add(p)

        self.SIM.ucell_man = ucell_man

        p = ParameterType(init=init.detz_shift*1e-3, sigma=sigma.detz_shift,
                          minval=mins.detz_shift*1e-3, maxval=maxs.detz_shift*1e-3,
                          fix=fix.detz_shift,name="detz_shift",
                          center=centers.detz_shift,
                          beta=betas.detz_shift)
        P.add(p)

        self.set_slices("roi_id")  # this creates roi_id_unique
        refls_have_scales = "scale_factor" in list(self.refls.keys())
        for roi_id in self.roi_id_unique:
            slc = self.roi_id_slices[roi_id][0]
            if refls_have_scales:
                refl_idx = int(self.all_refls_idx[slc][0])
                init_scale = self.refls[refl_idx]["scale_factor"]
            else:
                init_scale = 1
            p = RangedParameter(init=init_scale, sigma=self.params.sigmas.roiPerScale,
                              minval=0, maxval=1e12,
                              fix=fix.perRoiScale, name="scale_roi%d" % roi_id,
                              center=1,
                              beta=1e12)
            if isinstance(self.all_q_perpix, np.ndarray) and self.all_q_perpix.size:
                q = self.all_q_perpix[slc][0]
                reso = 1./q
                hkl = self.all_nominal_hkl[slc][0]
                p.misc_data = reso, hkl
            P.add(p)


        # two parameters for optimizing the spectrum
        p = RangedParameter(init=self.params.init.spec[0], sigma=self.params.sigmas.spec[0], minval=-0.01, maxval=0.01, fix=fix.spec,
                            name="lambda_offset", center=0, beta=1e12)
        P.add(p)
        p = RangedParameter(init=self.params.init.spec[1], sigma=self.params.sigmas.spec[1], minval=0.95, maxval=1.05, fix=fix.spec,
                            name="lambda_scale", center=1, beta=1e12)
        P.add(p)
        if not fix.spec:
            self.SIM.D.use_lambda_coefficients = True
            self.SIM.D.lambda_coefficients = tuple(self.params.init.spec)

        self.SIM.P = P
        # TODO , fix this attribute hacking
        self.SIM.roi_id_unique = self.roi_id_unique
        self.SIM.roi_id_slices = self.roi_id_slices

    def get_data_model_pairs(self):
        if self.best_model is None:
            raise ValueError("cannot get the best model, there is no best_model attribute")
        all_dat_img, all_mod_img = [], []
        all_trusted = []
        all_bragg = []
        all_sigma_rdout = []
        for i_roi in range(len(self.rois)):
            x1, x2, y1, y2 = self.rois[i_roi]
            mod = self.best_model[self.roi_id == i_roi].reshape((y2 - y1, x2 - x1))
            if self.all_trusted is not None:
                trusted = self.all_trusted[self.roi_id == i_roi].reshape((y2 - y1, x2 - x1))
                all_trusted.append(trusted)
            else:
                all_trusted.append(None)

            dat = self.all_data[self.roi_id == i_roi].reshape((y2 - y1, x2 - x1))
            all_dat_img.append(dat)
            if isinstance(self.all_sigma_rdout, np.ndarray):
                sig = self.all_sigma_rdout[self.roi_id==i_roi].reshape((y2-y1, x2-x1))
                all_sigma_rdout.append(sig)
            if self.all_background is not None:
                bg = self.all_background[self.roi_id == i_roi].reshape((y2-y1, x2-x1))
                if self.best_model_includes_background:
                    all_bragg.append(mod-bg)
                    all_mod_img.append(mod)
                else:
                    all_bragg.append(mod)
                    all_mod_img.append(mod+bg)
            else:  # assume mod contains background
                all_mod_img.append(mod)
                all_bragg.append(None)
        ret_subimgs = all_dat_img, all_mod_img, all_trusted, all_bragg
        if all_sigma_rdout:
            ret_subimgs += (all_sigma_rdout, )
        return ret_subimgs

    def Minimize(self, x0):
        target = TargetFunc(SIM=self.SIM, niter_per_J=self.params.niter_per_J, profile=self.params.profile)

        # set up the refinement flags
        vary = np.ones(len(x0), bool)
        assert len(x0) == len(self.SIM.P)
        for p in self.SIM.P.values():
            if not p.refine:
                vary[p.xpos] = False

        target.vary = vary  # fixed flags
        target.x0 = np.array(x0, np.float64)  # initial full parameter list
        x0_for_refinement = target.x0[vary]

        if self.params.method is None:
            method = "Nelder-Mead"
        else:
            method = self.params.method

        maxfev = None
        if self.params.nelder_mead_maxfev is not None:
            maxfev = self.params.nelder_mead_maxfev * self.npix_total

        at_min = None
        if self.params.logging.show_params_at_minimum:
            at_min = target.at_minimum
        callback = lambda x: target.callback(x, verbose=self.params.logging.parameters)
        target.terminate_after_n_converged_iterations = self.params.terminate_after_n_converged_iter
        target.percent_change_of_converged = self.params.converged_param_percent_change
        if method in ["L-BFGS-B", "BFGS", "CG", "dogleg", "SLSQP", "Newton-CG", "trust-ncg", "trust-krylov", "trust-exact", "trust-ncg"]:
            if self.SIM.P["lambda_offset"].refine:
                for lam_id in LAMBDA_IDS:
                    self.SIM.D.refine(lam_id)
            if self.SIM.P["RotXYZ0_xtal0"].refine:
                self.SIM.D.refine(ROTX_ID)
                self.SIM.D.refine(ROTY_ID)
                self.SIM.D.refine(ROTZ_ID)
            if self.SIM.P["Nabc0"].refine:
                self.SIM.D.refine(NCELLS_ID)
            if self.SIM.P["ucell0"].refine:
                for i_ucell in range(len(self.SIM.ucell_man.variables)):
                    self.SIM.D.refine(UCELL_ID_OFFSET + i_ucell)
            if self.SIM.P["eta_abc0"].refine:
                self.SIM.D.refine(ETA_ID)
            if self.SIM.P["detz_shift"].refine:
                self.SIM.D.refine(DETZ_ID)
            if self.SIM.D.use_diffuse:
                self.SIM.D.refine(DIFFUSE_ID)

            args = (self.SIM, self.pan_fast_slow, self.all_data,
                    self.all_sigma_rdout, self.all_trusted, self.all_background, True, self.params, True)
            min_kwargs = {'args': args, "method": method, "jac": target.jac,
                          'hess': self.params.hess, 'callback':callback}
            if method=="L-BFGS-B":
                min_kwargs["options"] = {"ftol": self.params.ftol, "gtol": 1e-10, "maxfun":1e5, "maxiter":self.params.lbfgs_maxiter}
        else:
            args = (self.SIM, self.pan_fast_slow, self.all_data,
                    self.all_sigma_rdout, self.all_trusted, self.all_background, True, self.params, False)
            min_kwargs = {'args': args, "method": method,
                          'callback': callback,
                          'options': {'maxfev': maxfev,
                                      'fatol': self.params.nelder_mead_fatol}}

        if self.params.global_method=="basinhopping":
            HOPPER = basinhopping

            try:
                out = HOPPER(target, x0_for_refinement,
                                   niter=self.params.niter,
                                   minimizer_kwargs=min_kwargs,
                                   T=self.params.temp,
                                   callback=at_min,
                                   disp=False,
                                   stepsize=self.params.stepsize)
                target.x0[vary] = out.x
            except StopIteration:
                pass

        else:
            bounds = [(-100,100)] * len(x0_for_refinement)  # TODO decide about bounds, usually x remains close to 1 during refinement
            print("Beginning the annealing process")
            args = min_kwargs.pop("args")
            if self.params.dual.no_local_search:
                compute_grads = args[-1]
                if compute_grads:
                    print("Warning, parameters setup to compute gradients, swicthing off because no_local_search=True")
                args = list(args)
                args[-1] = False  # switch off grad
                args = tuple(args)
            out = dual_annealing(target, bounds=bounds, args=args,
                                 no_local_search=self.params.dual.no_local_search,
                                 x0=x0_for_refinement,
                                 accept=self.params.dual.accept,
                                 visit=self.params.dual.visit,
                                 maxiter=self.params.niter,
                                 local_search_options=min_kwargs,
                                 callback=at_min)
            target.x0[vary] = out.x

        return target.x0


def model(x, SIM, pfs,  compute_grad=True, dont_rescale_gradient=False):

    #params_per_xtal = np.array_split(x[:num_per_xtal_params], SIM.num_xtals)

    # get the unit cell variables
    nucell = len(SIM.ucell_man.variables)
    ucell_params = [SIM.P["ucell%d" % i_uc] for i_uc in range(nucell)]
    ucell_xpos = [p.xpos for p in ucell_params]
    unitcell_var_reparam = [x[xpos] for xpos in ucell_xpos]
    unitcell_variables = [ucell_params[i].get_val(xval) for i, xval in enumerate(unitcell_var_reparam)]
    SIM.ucell_man.variables = unitcell_variables
    Bmatrix = SIM.ucell_man.B_recipspace
    SIM.D.Bmatrix = Bmatrix
    if compute_grad:
        for i_ucell in range(len(unitcell_variables)):
            SIM.D.set_ucell_derivative_matrix(
                i_ucell + UCELL_ID_OFFSET,
                SIM.ucell_man.derivative_matrices[i_ucell])

    # update the mosaicity here
    eta_params = [SIM.P["eta_abc%d" % i_eta] for i_eta in range(3)]
    if SIM.umat_maker is not None:
        # we are modeling mosaic spread
        eta_abc = [p.get_val(x[p.xpos]) for p in eta_params]
        if not SIM.D.has_anisotropic_mosaic_spread:
            eta_abc = eta_abc[0]
        SIM.update_umats_for_refinement(eta_abc)

#   detector parameters
    DetZ = SIM.P["detz_shift"]
    x_shiftZ = x[DetZ.xpos]
    shiftZ = DetZ.get_val(x_shiftZ)
    SIM.D.shift_origin_z(SIM.detector, shiftZ)

    if SIM.P["lambda_offset"].refine:
        p0 = SIM.P["lambda_offset"]
        p1 = SIM.P["lambda_scale"]
        lambda_coef = p0.get_val(x[p0.xpos]), p1.get_val(x[p1.xpos])
        SIM.D.lambda_coefficients = lambda_coef

    # Mosaic block
    Nabc_params = [SIM.P["Nabc%d" % (i_n,)] for i_n in range(3)]
    Na, Nb, Nc = [n_param.get_val(x[n_param.xpos]) for n_param in Nabc_params]
    if SIM.D.isotropic_ncells:
        Nb = Na
        Nc = Na
    SIM.D.set_ncells_values(tuple([Na, Nb, Nc]))

    # diffuse signals
    if SIM.D.use_diffuse:
        diffuse_params_lookup = {}
        iso_flags = {'gamma':SIM.isotropic_diffuse_gamma, 'sigma':SIM.isotropic_diffuse_sigma}
        for diff_type in ['gamma', 'sigma']:
            diff_params = [SIM.P["diffuse_%s%d" % (diff_type,i_gam)] for i_gam in range(3)]
            diffuse_params_lookup[diff_type] = diff_params
            diff_vals = []
            for i_diff, param in enumerate(diff_params):
                val = param.get_val(x[param.xpos])
                if iso_flags[diff_type]:
                    diff_vals = [val]*3
                    break
                else:
                    diff_vals.append(val)
            if diff_type == "gamma":
                SIM.D.diffuse_gamma = tuple(diff_vals)
            else:
                SIM.D.diffuse_sigma = tuple(diff_vals)

    npix = int(len(pfs) / 3)
    nparam = len(x)
    J = None
    if compute_grad:
        J = np.zeros((nparam, npix))  # gradients
    model_pix = None
    model_pix_noRoi = None

    # extract the scale factors per ROI, these might correspond to structure factor intensity scale factors, and quite possibly might result in overfits!
    roiScalesPerPix = 1
    if SIM.P["scale_roi0"].refine:
        perRoiParams = [SIM.P["scale_roi%d" % roi_id] for roi_id in SIM.roi_id_unique]
        perRoiScaleFactors = [p.get_val(x[p.xpos]) for p in perRoiParams]
        roiScalesPerPix = np.zeros(npix)
        for i_roi, roi_id in enumerate(SIM.roi_id_unique):
            slc = SIM.roi_id_slices[roi_id][0]
            roiScalesPerPix[slc] = perRoiScaleFactors[i_roi]

    for i_xtal in range(SIM.num_xtals):
        SIM.D.raw_pixels_roi *= 0

        RotXYZ_params = [SIM.P["RotXYZ%d_xtal%d" % (i_rot, i_xtal)] for i_rot in range(3)]
        rotX,rotY,rotZ = [rot_param.get_val(x[rot_param.xpos]) for rot_param in RotXYZ_params]

        ## update parameters:
        # TODO: if not refining Umat, assert these are 0 , and dont set them here
        SIM.D.set_value(ROTX_ID, rotX)
        SIM.D.set_value(ROTY_ID, rotY)
        SIM.D.set_value(ROTZ_ID, rotZ)

        G = SIM.P["G_xtal%d" % i_xtal]
        scale = G.get_val(x[G.xpos])

        SIM.D.add_diffBragg_spots(pfs)

        pix_noRoiScale = SIM.D.raw_pixels_roi[:npix]
        pix_noRoiScale = pix_noRoiScale.as_numpy_array()

        pix = pix_noRoiScale * roiScalesPerPix

        if model_pix is None:
            model_pix = scale*pix
            model_pix_noRoi = scale*pix_noRoiScale
        else:
            model_pix += scale*pix
            model_pix_noRoi += scale*pix_noRoiScale

        if compute_grad:
            if G.refine:
                scale_grad = pix  # TODO double check multi crystal case
                scale_grad = G.get_deriv(x[G.xpos], scale_grad)
                J[G.xpos] += scale_grad

            if RotXYZ_params[0].refine:
                for i_rot in range(3):
                    rot_grad = scale * SIM.D.get_derivative_pixels(ROTXYZ_IDS[i_rot]).as_numpy_array()[:npix]
                    rot_p = RotXYZ_params[i_rot]
                    rot_grad = rot_p.get_deriv(x[rot_p.xpos], rot_grad)
                    J[rot_p.xpos] += rot_grad

            if Nabc_params[0].refine:
                Nabc_grads = SIM.D.get_ncells_derivative_pixels()
                for i_n in range(3):
                    N_grad = scale*(Nabc_grads[i_n][:npix].as_numpy_array())
                    p = Nabc_params[i_n]
                    N_grad = p.get_deriv(x[p.xpos], N_grad)
                    J[p.xpos] += N_grad
                    if SIM.D.isotropic_ncells:
                        break

            if SIM.D.use_diffuse:
                for t in ['gamma','sigma']:
                    if diffuse_params_lookup[t][0].refine:
                        diffuse_grads = getattr(SIM.D,"get_diffuse_%s_derivative_pixels"%t)()
                        for i_diff in range(3):
                            diff_grad = scale*(diffuse_grads[i_diff][:npix].as_numpy_array())
                            p = diffuse_params_lookup[t][i_diff]
                            diff_grad = p.get_deriv(x[p.xpos], diff_grad)
                            J[p.xpos] += diff_grad

            if eta_params[0].refine:
                if SIM.D.has_anisotropic_mosaic_spread:
                    eta_derivs = SIM.D.get_aniso_eta_deriv_pixels()
                else:
                    eta_derivs = [SIM.D.get_derivative_pixels(ETA_ID)]
                num_eta = 3 if SIM.D.has_anisotropic_mosaic_spread else 1
                for i_eta in range(num_eta):
                    p = eta_params[i_eta]
                    eta_grad = scale * (eta_derivs[i_eta][:npix].as_numpy_array())
                    eta_grad = p.get_deriv(x[p.xpos], eta_grad)
                    J[p.xpos] += eta_grad

            if ucell_params[0].refine:
                for i_ucell in range(nucell):
                    p = ucell_params[i_ucell]
                    deriv = scale*SIM.D.get_derivative_pixels(UCELL_ID_OFFSET+i_ucell).as_numpy_array()[:npix]
                    deriv = p.get_deriv(x[p.xpos], deriv)
                    J[p.xpos] += deriv

            if DetZ.refine:
                d = SIM.D.get_derivative_pixels(DETZ_ID).as_numpy_array()[:npix]
                d = DetZ.get_deriv(x[DetZ.xpos], d)
                J[DetZ.xpos] += d

            if SIM.P["lambda_offset"].refine:
                lambda_derivs = SIM.D.get_lambda_derivative_pixels()
                lambda_param_names = "lambda_offset", "lambda_scale"
                for d,name in zip(lambda_derivs, lambda_param_names):
                    p = SIM.P[name]
                    d = d.as_numpy_array()[:npix]
                    d = p.get_deriv(x[p.xpos], d)
                    J[p.xpos] += d

    if SIM.P["scale_roi0"].refine and compute_grad:
        if compute_grad:
            for p in perRoiParams:
                roi_id = int(p.name.split("scale_roi")[1])
                slc = SIM.roi_id_slices[roi_id][0]
                if dont_rescale_gradient:
                    d = model_pix_noRoi[slc]
                else:
                    d = p.get_deriv(x[p.xpos], model_pix_noRoi[slc])
                J[p.xpos, slc] += d

    return model_pix, J


def look_at_x(x, SIM):
    for name, p in SIM.P.items():
        if name.startswith("scale_roi") and not p.refine:
            continue
        val = p.get_val(x[p.xpos])
        print("%s: %f" % (name, val))


def get_param_from_x(x, SIM, i_xtal=0):
    G = SIM.P['G_xtal%d' %i_xtal]
    scale = G.get_val(x[G.xpos])

    RotXYZ = [SIM.P["RotXYZ%d_xtal%d" % (i, i_xtal)] for i in range(3)]
    rotX, rotY, rotZ = [r.get_val(x[r.xpos]) for r in RotXYZ]

    Nabc = [SIM.P["Nabc%d" % (i, )] for i in range(3)]
    Na, Nb, Nc = [p.get_val(x[p.xpos]) for p in Nabc]

    diff_gam_abc = [SIM.P["diffuse_gamma%d" % i] for i in range(3)]
    diff_gam_a, diff_gam_b, diff_gam_c = [p.get_val(x[p.xpos]) for p in diff_gam_abc]

    diff_sig_abc = [SIM.P["diffuse_sigma%d" % i] for i in range(3)]
    diff_sig_a, diff_sig_b, diff_sig_c = [p.get_val(x[p.xpos]) for p in diff_sig_abc]

    nucell = len(SIM.ucell_man.variables)
    ucell_p = [SIM.P["ucell%d" % i] for i in range(nucell)]
    ucell_var = [p.get_val(x[p.xpos]) for p in ucell_p]
    SIM.ucell_man.variables = ucell_var
    a,b,c,al,be,ga = SIM.ucell_man.unit_cell_parameters

    DetZ = SIM.P["detz_shift"]
    detz = DetZ.get_val(x[DetZ.xpos])

    return scale, rotX, rotY, rotZ, Na, Nb, Nc, diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c, a,b,c,al,be,ga, detz


class TargetFunc:
    def __init__(self, SIM, niter_per_J=1, profile=False):
        self.niter_per_J = niter_per_J
        self.prev_iter_vals = None
        self.global_x = []
        self.percent_change_of_converged = 0.1
        self.all_x = []
        self.terminate_after_n_converged_iterations = None
        self.vary = None #boolean numpy array specifying which params to refine
        self.x0 = None  # 1d array of parameters (should be numpy array, same length as vary)
        self.old_J = None
        self.old_model = None
        self.delta_x = None
        self.iteration = 0
        self.minima = []
        self.all_converged_params = 0
        self.SIM = SIM

    def at_minimum(self, x, f, accept):
        self.iteration = 0
        self.all_x = []
        self.x0[self.vary] = x
        look_at_x(self.x0,self.SIM)
        self.minima.append((f,self.x0,accept))

    def jac(self, x, *args):
        if self.g is not None:
            return self.g[self.vary]

    def __call__(self, x, *args, **kwargs):
        self.x0[self.vary] = x
        if self.all_x:
            self.delta_x = self.x0 - self.all_x[-1]
        update_terms = None
        if not self.iteration % (self.niter_per_J) == 0:
            update_terms = (self.delta_x, self.old_J, self.old_model)
        self.all_x.append(self.x0)
        f, g, modelpix, J = target_func(self.x0, update_terms, *args, **kwargs)
        self.old_model = modelpix
        self.old_J = J
        self.iteration += 1
        self.g = g
        return f

    def callback(self, x, verbose=True):
        xall = self.x0.copy()
        xall[self.vary] = x
        rescaled_vals = np.zeros_like(xall)
        all_perc_change = []
        for name in self.SIM.P:
            p = self.SIM.P[name]

            if not p.refine:
                continue
            xpos = p.xpos
            val = p.get_val(xall[xpos])
            log_s = "Iter %d: %s = %1.2g ."  %(self.iteration, name, val)
            if name.startswith("scale_roi") and p.misc_data is not None:
                reso, (h,k,l) = p.misc_data
                log_s += "reso=%1.3f Ang. (h,k,l)=%d,%d,%d ." %(reso, h,k,l)
            rescaled_vals[xpos] = val
            if self.prev_iter_vals is not None:
                prev_val = self.prev_iter_vals[xpos]
                if val - prev_val == 0 and prev_val == 0:
                    val_percent_diff = 0
                elif prev_val == 0:
                    val_percent_diff = np.abs(val - prev_val) / val *100.
                else:
                    val_percent_diff = np.abs(val - prev_val) / prev_val *100.
                log_s += " Percent change = %1.2f%%" % val_percent_diff
                all_perc_change.append(val_percent_diff)

            if verbose:
                MAIN_LOGGER.debug(log_s)

        if all_perc_change:
            all_perc_change = np.abs(all_perc_change)
            ave_perc_change = np.mean(all_perc_change)
            num_perc_change_small = np.sum(all_perc_change < self.percent_change_of_converged)
            max_perc_change = np.max(all_perc_change)
            if num_perc_change_small == len(all_perc_change):
                self.all_converged_params += 1
            else:
                self.all_converged_params = 0
            if verbose:
                MAIN_LOGGER.info("Iter %d: Mean percent change = %1.2f%%. Num param with %%-change < %1.2f: %d/%d. Max %%-change=%1.2f%%"
                                  % (self.iteration, ave_perc_change, self.percent_change_of_converged,
                                     num_perc_change_small, len(all_perc_change), max_perc_change))

        self.prev_iter_vals = rescaled_vals
        if self.terminate_after_n_converged_iterations is not None and self.all_converged_params >= self.terminate_after_n_converged_iterations:
            # at this point prev_iter_vals are the converged parameters!
            raise StopIteration()  # Refinement has reached convergence!


def target_func(x, udpate_terms, SIM, pfs, data, sigma_rdout, trusted, background, verbose=True, params=None, compute_grad=True):

    if udpate_terms is not None:
        # if approximating the gradients, then fix the parameter refinment managers in diffBragg
        # so we dont waste time computing them
        _compute_grad = False
        SIM.D.fix(NCELLS_ID)
        SIM.D.fix(ROTX_ID)
        SIM.D.fix(ROTY_ID)
        SIM.D.fix(ROTZ_ID)
        for lam_id in LAMBDA_IDS:
            SIM.D.fix(lam_id)
        for i_ucell in range(len(SIM.ucell_man.variables)):
            SIM.D.fix(UCELL_ID_OFFSET + i_ucell)
        SIM.D.fix(DETZ_ID)
        SIM.D.fix(ETA_ID)
        SIM.D.fix(DIFFUSE_ID)
    elif compute_grad:
        # actually compute the gradients
        _compute_grad = True
        if SIM.P["Nabc0"].refine:
            SIM.D.let_loose(NCELLS_ID)
        if SIM.P["RotXYZ0_xtal0"].refine:
            SIM.D.let_loose(ROTX_ID)
            SIM.D.let_loose(ROTY_ID)
            SIM.D.let_loose(ROTZ_ID)
        if SIM.P["ucell0"].refine:
            for i_ucell in range(len(SIM.ucell_man.variables)):
                SIM.D.let_loose(UCELL_ID_OFFSET + i_ucell)
        if SIM.P["detz_shift"].refine:
            SIM.D.let_loose(DETZ_ID)
        if SIM.P["eta_abc0"].refine:
            SIM.D.let_loose(ETA_ID)
        if SIM.P["lambda_offset"].refine:
            for lam_id in LAMBDA_IDS:
                SIM.D.let_loose(lam_id)
    else:
        _compute_grad = False
    model_bragg, Jac = model(x, SIM, pfs, compute_grad=_compute_grad)

    if udpate_terms is not None:
        # try a Broyden update ?
        # https://people.duke.edu/~hpgavin/ce281/lm.pdf  equation 19
        delta_x, prev_J, prev_model_bragg = udpate_terms
        if prev_J is not None:
            delta_y = model_bragg - prev_model_bragg

            delta_J = (delta_y - np.dot(prev_J.T, delta_x))
            delta_J /= np.dot(delta_x,delta_x)
            Jac = prev_J + delta_J
    # Jac has shape of num_param x num_pix

    model_pix = model_bragg + background

    resid = data - model_pix

    # data contributions to target function
    V = model_pix + sigma_rdout**2
    # TODO:what if V is allowed to be negative? The logarithm/sqrt will explore below
    resid_square = resid**2
    fLogLike = (.5*(np.log(2*np.pi*V) + resid_square / V))[trusted].sum()   # negative log Likelihood target

    # width of z-score should decrease as refinement proceeds
    zscore_sigma = np.std( (resid / np.sqrt(V))[trusted])

    restraint_terms = {}
    if params.use_restraints:
        # scale factor restraint
        for name in SIM.P:
            if name.startswith("scale_roi"):
                continue  # No restraints for the per-spot roi scale factors
            p = SIM.P[name]
            val = p.get_restraint_val(x[p.xpos])
            restraint_terms[name] = val

        if params.centers.Nvol is not None:
            Nvol = np.product(SIM.D.Ncells_abc_aniso)
            del_Nvol = params.centers.Nvol - Nvol
            fN_vol = .5*del_Nvol**2/params.betas.Nvol
            restraint_terms["Nvol"] = fN_vol

#   accumulate target function
    f_restraints = 0
    if restraint_terms:
        f_restraints = np.sum(list(restraint_terms.values()))
    f = f_restraints + fLogLike

    restraint_debug_s = "LogLike: %.1f%%; " % (fLogLike / f *100.)
    for name, val in restraint_terms.items():
        if val > 0:
            frac_total = val / f *100.
            restraint_debug_s += "%s: %.2f%%; " % (name, frac_total)

    # fractions of the target function
    g = None  # gradient vector
    gnorm = -1  # norm of gradient vector
    if compute_grad:
        common_grad_term_all = (0.5 /V * (1-2*resid - resid_square / V))
        common_grad_term = common_grad_term_all[trusted]

        g = np.zeros(Jac.shape[0])
        for name ,p in SIM.P.items():
            if not p.refine:
                continue
            Jac_p = Jac[p.xpos]

            if name.startswith("scale_roi"):
                roi_id = int(p.name.split("scale_roi")[1])
                slc = SIM.roi_id_slices[roi_id][0]
                common_grad_slc = common_grad_term_all[slc]
                Jac_slc = Jac_p[slc]
                trusted_slc = trusted[slc]
                g[p.xpos] += (Jac_slc*common_grad_slc)[trusted_slc].sum()
            else:
                g[p.xpos] += (Jac_p[trusted]*common_grad_term).sum()

        # trusted pixels portion of Jacobian
        #Jac_t = Jac[:,trusted]
        # gradient vector
        #g = np.array([np.sum(common_grad_term*Jac_t[param_idx]) for param_idx in range(Jac_t.shape[0])])

        if params.use_restraints:
            # update gradients according to restraints
            for name, p in SIM.P.items():
                g[p.xpos] += p.get_restraint_deriv(x[p.xpos])

            if params.centers.Nvol is not None:
                Na,Nb,Nc = SIM.D.Ncells_abc_aniso
                dNvol_dN = Nb*Nc, Na*Nc, Na*Nb
                for i_N in range(3):
                    p = SIM.P["Nabc%d" % i_N]
                    gterm = -del_Nvol / params.betas.Nvol * dNvol_dN[i_N]
                    g[p.xpos] += p.get_deriv(x[p.xpos], gterm)

        gnorm = np.linalg.norm(g)

    if verbose:
        MAIN_LOGGER.debug("F=%10.7g sigZ=%10.7g (Fracs of F: %s), |g|=%10.7g" \
              % (f, zscore_sigma, restraint_debug_s, gnorm))

    return f, g, model_bragg, Jac


def refine(exp, ref, params, spec=None, gpu_device=None, return_modeler=False, best=None, free_mem=True):
    if gpu_device is None:
        gpu_device = 0
    params.simulator.spectrum.filename = spec
    Modeler = DataModeler(params)
    if params.load_data_from_refls:
        Modeler.GatherFromReflectionTable(exp, ref)
    else:
        assert Modeler.GatherFromExperiment(exp, ref)

    Modeler.SimulatorFromExperiment(best)

    Modeler.SIM.rois = Modeler.rois
    Modeler.SIM.pids = Modeler.pids
    Modeler.SIM.all_data = Modeler.all_data
    Modeler.SIM.E = Modeler.E

    Modeler.SIM.D.device_Id = gpu_device

    nparam = len(Modeler.SIM.P)
    x0 = [1] * nparam

    x = Modeler.Minimize(x0)
    Modeler.best_model, _ = model(x, Modeler.SIM, Modeler.pan_fast_slow, compute_grad=False)
    Modeler.best_model_includes_background = False

    new_crystal = update_crystal_from_x(Modeler.SIM, x)
    new_exp = deepcopy(Modeler.E)
    new_exp.crystal = new_crystal

    try:
        new_exp.beam.set_wavelength(Modeler.SIM.dxtbx_spec.get_weighted_wavelength())
    except Exception: pass
    # if we strip the thickness from the detector, then update it here:
    #new_exp.detector. shift Z mm
    new_det = update_detector_from_x(Modeler.SIM, x)
    new_exp.detector = new_det

    new_refl = get_new_xycalcs(Modeler, new_exp)

    if free_mem:
        Modeler.clean_up()

    if return_modeler:
        return new_exp, new_refl, Modeler, x

    else:
        return new_exp, new_refl


def update_detector_from_x(SIM, x):
    scale, rotX, rotY, rotZ, Na, Nb, Nc, _,_,_,_,_,_,a, b, c, al, be, ga, detz_shift = get_param_from_x(x, SIM)
    detz_shift_mm = detz_shift*1e3
    det = SIM.detector
    det = utils.shift_panelZ(det, detz_shift_mm)
    return det


def new_cryst_from_rotXYZ_and_ucell(rotXYZ, ucparam, orig_crystal):
    """

    :param rotXYZ: tuple of rotation angles about princple axes (radians)
        Crystal will be rotated by rotX about (-1,0,0 ), rotY about (0,-1,0) and rotZ about
        (0,0,-1) . This is the convention used by diffBragg
    :param ucparam: unit cell parameter tuple (a,b,c, alpha, beta, gamma) in Angstrom,degrees)
    :param orig_crystal: dxtbx.model.Crystal object which will be copied and adjusted according to ucell and rotXYZ
    :return: new dxtbx.model.Crystal object
    """
    rotX, rotY, rotZ = rotXYZ
    xax = col((-1, 0, 0))
    yax = col((0, -1, 0))
    zax = col((0, 0, -1))
    ## update parameters:
    RX = xax.axis_and_angle_as_r3_rotation_matrix(rotX, deg=False)
    RY = yax.axis_and_angle_as_r3_rotation_matrix(rotY, deg=False)
    RZ = zax.axis_and_angle_as_r3_rotation_matrix(rotZ, deg=False)
    M = RX * RY * RZ
    U = M * sqr(orig_crystal.get_U())
    new_C = deepcopy(orig_crystal)
    new_C.set_U(U)

    ucman = utils.manager_from_params(ucparam)
    new_C.set_B(ucman.B_recipspace)
    return new_C


def update_crystal_from_x(SIM, x):
    """
    :param SIM: sim_data instance containing a nanoBragg crystal object
    :param x: parameters returned by hopper_utils (instance of simtbx.diffBragg.refiners.parameters.Parameters()
    :return: a new dxtbx.model.Crystal object with updated unit cell and orientation matrix
    """
    scale, rotX, rotY, rotZ, Na, Nb, Nc, _,_,_,_,_,_,a, b, c, al, be, ga, detz_shift = get_param_from_x(x, SIM)
    ucparam = a, b, c, al, be, ga
    return new_cryst_from_rotXYZ_and_ucell((rotX,rotY,rotZ), ucparam, SIM.crystal.dxtbx_crystal)


def get_new_xycalcs(Modeler, new_exp, old_refl_tag="dials"):
    """

    :param Modeler: data modeler instance after refinement
    :param new_exp: post-refinement dxtbx experiment obj
    :param old_refl_tag: exisiting columns will be renamed with this tag as a prefix
    :return: refl table with xyzcalcs derived from the modeling results
    """
    if isinstance(Modeler.all_sigma_rdout, np.ndarray):
        _, _, _, bragg_subimg, _ = Modeler.get_data_model_pairs()
    else:
        _,_,_, bragg_subimg = Modeler.get_data_model_pairs()
    new_refls = deepcopy(Modeler.refls)

    reflkeys = list(new_refls.keys())
    if "xyzcal.px" in reflkeys:
        new_refls['%s.xyzcal.px' % old_refl_tag] = deepcopy(new_refls['xyzcal.px'])
    if "xyzcal.mm" in reflkeys:
        new_refls['%s.xyzcal.mm' % old_refl_tag] = deepcopy(new_refls['xyzcal.mm'])
    if "xyzobs.mm.value" in list(new_refls.keys()):
        new_refls['%s.xyzobs.mm.value' % old_refl_tag] = deepcopy(new_refls['xyzobs.mm.value'])

    new_xycalcs = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    new_xycalcs_mm = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    new_xyobs_mm = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    for i_roi in range(len(bragg_subimg)):

        ref_idx = Modeler.refls_idx[i_roi]

        if np.any(bragg_subimg[i_roi] > 0):
            I = bragg_subimg[i_roi]
            Y, X = np.indices(bragg_subimg[i_roi].shape)
            x1, _, y1, _ = Modeler.rois[i_roi]
            X += x1
            Y += y1
            Isum = I.sum()
            xcom = (X * I).sum() / Isum + .5
            ycom = (Y * I).sum() / Isum + .5
            com = xcom, ycom, 0

            pid = Modeler.pids[i_roi]
            assert pid == new_refls[ref_idx]['panel']
            panel = new_exp.detector[pid]
            xmm, ymm = panel.pixel_to_millimeter((xcom, ycom))
            com_mm = xmm, ymm, 0
            xobs, yobs, _ = new_refls[ref_idx]["xyzobs.px.value"]
            xobs_mm, yobs_mm = panel.pixel_to_millimeter((xobs, yobs))
            obs_com_mm = xobs_mm, yobs_mm, 0

            new_xycalcs[ref_idx] = com
            new_xycalcs_mm[ref_idx] = com_mm
            new_xyobs_mm[ref_idx] = obs_com_mm

    new_refls["xyzcal.px"] = new_xycalcs
    new_refls["xyzcal.mm"] = new_xycalcs_mm
    new_refls["xyzobs.mm.value"] = new_xyobs_mm

    if Modeler.params.filter_unpredicted_refls_in_output:
        sel = [not np.isnan(x) for x,_,_ in new_refls['xyzcal.px']]
        nbefore = len(new_refls)
        new_refls = new_refls.select(flex.bool(sel))
        nafter = len(new_refls)
        MAIN_LOGGER.info("Filtered %d / %d reflections which did not show peaks in model" % (nbefore-nafter, nbefore))

    return new_refls


def get_mosaicity_from_x(x, SIM):
    """
    :param x: refinement parameters
    :param SIM: simulator used during refinement
    :return: float or 3-tuple, depending on whether mosaic spread was modeled isotropically
    """
    eta_params = [SIM.P["eta_abc%d"%i] for i in range(3)]
    eta_abc = [p.get_val(x[p.xpos]) for p in eta_params]
    if not SIM.D.has_anisotropic_mosaic_spread:
        eta_abc = [eta_abc[0]]*3
    return eta_abc


def generate_gauss_spec(central_en=9500, fwhm=10, res=1, nchan=20, total_flux=1e12, as_spectrum=True):
    """
    In cases where raw experimental data do not contain X-ray spectra, one can generate a gaussian spectra
    and use that as part of the diffBragg model.

    :param central_en: nominal beam energy in eV
    :param fwhm: fwhm of desired Gaussian spectrum in eV
    :param res: energy resolution of spectrum in eV
    :param nchan: number of energy channels in spectrum
    :param total_flux: total number of photons across whole spectrum
    :param as_spectrum: return as a spectrum object, a suitable attribute of nanoBragg_beam
    :return: either a nanoBragg_beam.NBbeam spectrum, or a tuple of numpy arrays spec_energies, spec_weights
    """
    sig = fwhm / 2.35482
    min_en = central_en - nchan / 2 * res
    max_en = central_en + nchan / 2 * res
    ens = np.linspace(min_en, max_en, nchan)
    wt = np.exp(-(ens - central_en) ** 2 / 2 / sig ** 2)
    wt /= wt.sum()
    wt *= total_flux
    if as_spectrum:
        wvs = utils.ENERGY_CONV / ens  # wavelengths
        spec = list(zip(list(wvs), list(wt)))
        return spec
    else:
        return ens, wt


def downsamp_spec_from_params(params, expt):
    dxtbx_spec = expt.imageset.get_spectrum(0)
    spec_en = dxtbx_spec.get_energies_eV()
    spec_wt = dxtbx_spec.get_weights()
    if params.downsamp_spec.skip:
        spec_wave = utils.ENERGY_CONV / spec_en.as_numpy_array()
        spectrum = list(zip(spec_wave, spec_wt))
    else:
        spec_en = dxtbx_spec.get_energies_eV()
        spec_wt = dxtbx_spec.get_weights()
        # ---- downsample the spectrum
        method2_param = {"filt_freq": params.downsamp_spec.filt_freq,
                         "filt_order": params.downsamp_spec.filt_order,
                         "tail": params.downsamp_spec.tail,
                         "delta_en": params.downsamp_spec.delta_en}
        downsamp_en, downsamp_wt = downsample_spectrum(spec_en.as_numpy_array(),
                                                       spec_wt.as_numpy_array(),
                                                       method=2, method2_param=method2_param)

        stride = params.simulator.spectrum.stride
        if stride > len(downsamp_en) or stride == 0:
            raise ValueError("Incorrect value for pinkstride")
        downsamp_en = downsamp_en[::stride]
        downsamp_wt = downsamp_wt[::stride]
        tot_fl = params.simulator.total_flux
        if tot_fl is not None:
            downsamp_wt = downsamp_wt / sum(downsamp_wt) * tot_fl

        downsamp_wave = utils.ENERGY_CONV / downsamp_en
        spectrum = list(zip(downsamp_wave, downsamp_wt))
    # the nanoBragg beam has an xray_beams property that is used internally in diffBragg
    starting_wave = expt.beam.get_wavelength()
    waves, specs = map(np.array, zip(*spectrum))
    ave_wave = sum(waves*specs) / sum(specs)
    expt.beam.set_wavelength(ave_wave)
    MAIN_LOGGER.debug("Shifting wavelength from %f to %f" % (starting_wave, ave_wave))
    MAIN_LOGGER.debug("USING %d ENERGY CHANNELS" % len(spectrum))
    return spectrum


# set the X-ray spectra for this shot
def downsamp_spec(SIM, params, expt, return_and_dont_set=False):
    SIM.dxtbx_spec = expt.imageset.get_spectrum(0)
    spec_en = SIM.dxtbx_spec.get_energies_eV()
    spec_wt = SIM.dxtbx_spec.get_weights()
    if params.downsamp_spec.skip:
        spec_wave = utils.ENERGY_CONV / spec_en.as_numpy_array()
        SIM.beam.spectrum = list(zip(spec_wave, spec_wt))
    else:
        spec_en = SIM.dxtbx_spec.get_energies_eV()
        spec_wt = SIM.dxtbx_spec.get_weights()
        # ---- downsample the spectrum
        method2_param = {"filt_freq": params.downsamp_spec.filt_freq,
                         "filt_order": params.downsamp_spec.filt_order,
                         "tail": params.downsamp_spec.tail,
                         "delta_en": params.downsamp_spec.delta_en}
        downsamp_en, downsamp_wt = downsample_spectrum(spec_en.as_numpy_array(),
                                                       spec_wt.as_numpy_array(),
                                                       method=2, method2_param=method2_param)

        stride = params.simulator.spectrum.stride
        if stride > len(downsamp_en) or stride == 0:
            raise ValueError("Incorrect value for pinkstride")
        downsamp_en = downsamp_en[::stride]
        downsamp_wt = downsamp_wt[::stride]
        tot_fl = params.simulator.total_flux
        if tot_fl is not None:
            downsamp_wt = downsamp_wt / sum(downsamp_wt) * tot_fl

        downsamp_wave = utils.ENERGY_CONV / downsamp_en
        SIM.beam.spectrum = list(zip(downsamp_wave, downsamp_wt))
    # the nanoBragg beam has an xray_beams property that is used internally in diffBragg
    starting_wave = expt.beam.get_wavelength()
    waves, specs = map(np.array, zip(*SIM.beam.spectrum))
    ave_wave = sum(waves*specs) / sum(specs)
    expt.beam.set_wavelength(ave_wave)
    MAIN_LOGGER.debug("Shifting wavelength from %f to %f" % (starting_wave, ave_wave))
    if return_and_dont_set:
        return SIM.beam.spectrum
    else:
        SIM.D.xray_beams = SIM.beam.xray_beams


def set_gauss_spec(SIM, params, E):
    """
    Generates a gaussian spectrum for the experiment E and attach it to SIM
    in a way that diffBragg understands

    This updates the xray_beams property of SIM.D and the spectrum property
    of SIM.beam

    :param SIM: simulator object, instance of nanoBragg/sim_data
    :param params: phil params , expected scope is defined in diffBragg/phil.py
    :param E: dxtbx Experiment
    :return: None
    """
    if not hasattr(SIM, "D"):
        raise AttributeError("Cannot set the spectrum until diffBragg has been instantiated!")
    spec_mu = utils.ENERGY_CONV / E.beam.get_wavelength()
    spec_res = params.simulator.spectrum.gauss_spec.res
    spec_nchan = params.simulator.spectrum.gauss_spec.nchannels
    spec_fwhm = params.simulator.spectrum.gauss_spec.fwhm
    spec_ens, spec_wts = generate_gauss_spec(central_en=spec_mu, fwhm=spec_fwhm, res=spec_res,
                                             nchan=spec_nchan, as_spectrum=False)
    SIM.dxtbx_spec = Spectrum(spec_ens, spec_wts)
    spec_wvs = utils.ENERGY_CONV / spec_ens
    SIM.beam.spectrum = list(zip(spec_wvs, spec_wts))
    SIM.D.xray_beams = SIM.beam.xray_beams


def sanity_test_input_lines(input_lines):
    for line in input_lines:
        line_fields = line.strip().split()
        if len(line_fields) not in [2, 3]:
            raise IOError("Input line %s is not formatted properly" % line)
        for fname in line_fields:
            if not os.path.exists(fname):
                raise FileNotFoundError("File %s does not exist" % fname)


def full_img_pfs(img_sh):
    Panel_inds, Slow_inds, Fast_inds = map(np.ravel, np.indices(img_sh) )
    pfs_coords = np.vstack([Panel_inds, Fast_inds, Slow_inds]).T
    pfs_coords_flattened = pfs_coords.ravel()
    pfs_full = flex.size_t( np.ascontiguousarray(pfs_coords_flattened))
    return pfs_full


def print_profile(stats, timed_methods):
    for method in stats.timings.keys():
        filename, header_ln, name = method
        if name not in timed_methods:
            continue
        info = stats.timings[method]
        PROFILE_LOGGER.warning("\n")
        PROFILE_LOGGER.warning("FILE: %s" % filename)
        if not info:
            PROFILE_LOGGER.warning("<><><><><><><><><><><><><><><><><><><><><><><>")
            PROFILE_LOGGER.warning("METHOD %s : Not profiled because never called" % (name))
            PROFILE_LOGGER.warning("<><><><><><><><><><><><><><><><><><><><><><><>")
            continue
        unit = stats.unit

        line_nums, ncalls, timespent = zip(*info)
        fp = open(filename, 'r').readlines()
        total_time = sum(timespent)
        header_line = fp[header_ln-1][:-1]
        PROFILE_LOGGER.warning(header_line)
        PROFILE_LOGGER.warning("TOTAL FUNCTION TIME: %f ms" % (total_time*unit*1e3))
        PROFILE_LOGGER.warning("<><><><><><><><><><><><><><><><><><><><><><><>")
        PROFILE_LOGGER.warning("%5s%14s%9s%10s" % ("Line#", "Time", "%Time", "Line" ))
        PROFILE_LOGGER.warning("%5s%14s%9s%10s" % ("", "(ms)", "", ""))
        PROFILE_LOGGER.warning("<><><><><><><><><><><><><><><><><><><><><><><>")
        for i_l, l in enumerate(line_nums):
            frac_t = timespent[i_l] / total_time * 100.
            line = fp[l-1][:-1]
            PROFILE_LOGGER.warning("%5d%14.2f%9.2f%s" % (l, timespent[i_l]*unit*1e3, frac_t, line))
