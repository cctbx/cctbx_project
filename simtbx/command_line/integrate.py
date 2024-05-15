from __future__ import division, print_function

# LIBTBX_SET_DISPATCHER_NAME diffBragg.integrate

import argparse as ap
parser = ap.ArgumentParser()
parser.add_argument("predPhil", type=str, help="path to a phil config file for diffbragg prediction")
parser.add_argument("procPhil", type=str, help="path to a phil config file for stills process (used for spot finding and integration)")
parser.add_argument("inputGlob", type=str, help="glob of input pandas tables (those that are output by simtbx.diffBragg.hopper or diffBragg.hopper_process")
parser.add_argument("outdir", type=str, help="path to output refls")

parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
parser.add_argument("--dialsInteg", action="store_true", help="Integrate new shoeboxes using dials and write *integrated.expt files")
parser.add_argument("--numdev", type=int, default=1, help="number of GPUs (default=1)")
parser.add_argument("--pklTag", type=str, help="optional suffix for globbing for pandas pickles (default .pkl)", default=".pkl")
parser.add_argument("--loud", action="store_true", help="show lots of screen output")
parser.add_argument("--hopInputName", default="preds_for_hopper", type=str, help="write exp_ref_spec file and best_pickle pointing to the preditction models, such that one can run predicted rois through simtbx.diffBragg.hopper (e.g. to fit per-roi scale factors)")
parser.add_argument("--filterDupes", action="store_true", help="filter refls with same HKL")
parser.add_argument("--keepShoeboxes", action="store_true", help="Optionally keep shoeboxes present in the prediction refl tables (can lead to OOM errors)")
parser.add_argument("--scanWeakFracs", action="store_true", help="optionally stores a variety of inputs for stage2 based filtering different fractions of weak reflections")

args = parser.parse_args()

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

import logging
if not args.loud:
    logging.disable(logging.CRITICAL)
else:
    if COMM.rank==0:
        logger = logging.getLogger("diffBragg.main")
        logger.setLevel(logging.DEBUG)


def printR(*args, **kwargs):
    print("RANK %d" % COMM.rank, *args, **kwargs)
def print0(*args, **kwargs):
    if COMM.rank==0:
        print(*args, **kwargs)

import numpy as np
import json
from simtbx.diffBragg import hopper_utils, utils
from simtbx.modeling import predictions
from simtbx.diffBragg.hopper_utils import downsamp_spec_from_params
import glob
import pandas
import os
from dials.algorithms.integration.stills_significance_filter import SignificanceFilter
from dials.algorithms.indexing.stills_indexer import calc_2D_rmsd_and_displacements
import sys


def filter_weak_reflections(refls, weak_fraction):
    """
    :param pred:  reflection table created by this script
    :param weak_fraction: number from 0-1 (if 0, only strong spots are saved)
    :return: new reflection table with weak reflections filtered according to weak_fraction
    """
    new_refls = None
    for idx in set(refls['id']):
        pred = refls.select(refls['id']==idx)
        weaks = pred.select(pred['is_weak'])
        nweak = len(weaks)
        weaks_sorted = np.argsort(weaks["scatter"])[::-1]
        num_keep = int(nweak * weak_fraction)
        weak_refl_inds_keep = set(np.array(weaks["refl_idx"])[weaks_sorted[:num_keep]])
        weak_sel = flex.bool([i in weak_refl_inds_keep for i in pred['refl_idx']])
        keeps = np.logical_or(pred['is_strong'], weak_sel)
        pred = pred.select(flex.bool(keeps))
        if new_refls is None:
            new_refls = deepcopy(pred)
        else:
            new_refls.extend(pred)
    return new_refls


# Note: these imports and following 3 methods will eventually be in CCTBX/simtbx/diffBragg/utils
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import FilterRunner
from dials.model.data import PixelListLabeller, PixelList
from dials.algorithms.spot_finding.finder import pixel_list_to_reflection_table
from libtbx.phil import parse
from dials.command_line.stills_process import phil_scope
from dials.algorithms.integration.integrator import create_integrator
from dials.algorithms.profile_model.factory import ProfileModelFactory
from dxtbx.model import ExperimentList
from dials.array_family import flex

from copy import deepcopy
from collections import Counter

from simtbx.diffBragg.device import DeviceWrapper

def filter_refls(R):
    vec3_dbl_keys = 'xyzcal.px', 'xyzcal.mm', 'xyzobs.px.value', 'xyzobs.px.value', 'rlp', 's1'

    hkl_dupes = [h for h,count in Counter(R['miller_index']).items() if count > 1]
    print("%d miller indices are duplicates" % len(hkl_dupes))
    hkls = list(R['miller_index'])
    Rnew = None #flex.reflection_table()
    ndupe = 0
    for hkl in hkl_dupes:
        is_h = [h==hkl for h in hkls]
        ndupe += np.sum(is_h)
        Rdupes = R.select(flex.bool(is_h))
        R0 = deepcopy(Rdupes[0:1])
        for k in vec3_dbl_keys:
            xyz = np.mean(Rdupes[k].as_numpy_array(),axis=0)
            R0[k] = flex.vec3_double(1, tuple(xyz))
        if Rnew is None:
            Rnew = R0
        else:
            Rnew.extend(R0)
    print("%d refls belong to duplicates hkls" % ndupe)

    hkl_singles = set(hkls).difference(hkl_dupes)
    for hkl in hkl_singles:
        is_h = [h==hkl for h in hkls]
        i_R = np.where(is_h)[0][0]
        refl = R[i_R: i_R+1]
        if Rnew is None:
            Rnew = refl
        else:
            Rnew.extend(refl)
    print("filtered %d / %d refls" % (len(Rnew), len(R)))
    return Rnew


for i,arg in enumerate(sys.argv):
    if os.path.isfile(arg) or os.path.isdir(arg):
        sys.argv[i] = os.path.abspath(arg)
print0("COMMANDLINE: libtbx.python %s" % " ".join(sys.argv))

def stills_process_params_from_file(phil_file):
    """
    :param phil_file: path to phil file for stills_process
    :return: phil params object
    """
    phil_file = open(phil_file, "r").read()
    user_phil = parse(phil_file)
    phil_sources = [user_phil]
    working_phil, unused = phil_scope.fetch(
        sources=phil_sources, track_unused_definitions=True)
    params = working_phil.extract()
    return params



def process_reference(reference):
    """Load the reference spots."""
    assert "miller_index" in reference
    assert "id" in reference
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(~mask)
    if mask.count(False) > 0:
        reference.del_selected(~mask)
    if len(reference) == 0:
        raise RuntimeError(
            """
    Invalid input for reference reflections.
    Expected > %d indexed spots, got %d
  """
            % (0, len(reference))
        )
    mask = reference["miller_index"] == (0, 0, 0)
    if mask.count(True) > 0:
        rubbish.extend(reference.select(mask))
        reference.del_selected(mask)
    mask = reference["id"] < 0
    if mask.count(True) > 0:
        raise RuntimeError(
            """
    Invalid input for reference reflections.
    %d reference spots have an invalid experiment id
  """
            % mask.count(True)
        )
    return reference, rubbish



def integrate(phil_file, experiments, indexed, predicted):
    """
    integrate a single experiment at the locations specified by the predicted table
    The predicted table should have a column specifying strong reflections
    """
    assert len(experiments)==1

    for refls in [predicted, indexed]:
        refls['id'] = flex.int(len(refls), 0)
        refls['entering'] = flex.bool(len(refls), False)
        eid = refls.experiment_identifiers()
        for k in eid.keys():
            del eid[k]
        eid[0] = '0'
    experiments[0].identifier = '0'

    params = stills_process_params_from_file(phil_file)
    indexed,_ = process_reference(indexed)
    experiments = ProfileModelFactory.create(params, experiments, indexed)

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()
    for expt_id, expt in enumerate(experiments):
        if (
                params.profile.gaussian_rs.parameters.sigma_b_cutoff is None
                or expt.profile.sigma_b()
                < params.profile.gaussian_rs.parameters.sigma_b_cutoff
        ):
            refls = indexed.select(indexed["id"] == expt_id)
            refls["id"] = flex.int(len(refls), len(new_experiments))
            del refls.experiment_identifiers()[expt_id]
            refls.experiment_identifiers()[len(new_experiments)] = expt.identifier
            new_reflections.extend(refls)
            new_experiments.append(expt)

    experiments = new_experiments
    indexed = new_reflections
    if len(experiments) == 0:
        raise RuntimeError("No experiments after filtering by sigma_b")

    predicted.match_with_reference(indexed)
    integrator = create_integrator(params, experiments, predicted)
    integrated = integrator.integrate()

    if params.significance_filter.enable:

        sig_filter = SignificanceFilter(params)
        filtered_refls = sig_filter(experiments, integrated)
        accepted_expts = ExperimentList()
        accepted_refls = flex.reflection_table()
        for expt_id, expt in enumerate(experiments):
            refls = filtered_refls.select(filtered_refls["id"] == expt_id)
            if len(refls) > 0:
                accepted_expts.append(expt)
                refls["id"] = flex.int(len(refls), len(accepted_expts) - 1)
                accepted_refls.extend(refls)

        if len(accepted_refls) == 0:
            raise RuntimeError("No reflections left after applying significance filter")
        experiments = accepted_expts
        integrated = accepted_refls

    # Delete the shoeboxes used for intermediate calculations, if requested
    if params.integration.debug.delete_shoeboxes and "shoebox" in integrated:
        del integrated["shoebox"]


    rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
    log_str = "RMSD indexed (px): %f\n" % rmsd_indexed
    for i in range(6):
        bright_integrated = integrated.select(
            (
                    integrated["intensity.sum.value"]
                    / flex.sqrt(integrated["intensity.sum.variance"])
            )
            >= i
        )
        if len(bright_integrated) > 0:
            rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
        else:
            rmsd_integrated = 0
        log_str += (
                "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"
                % (i, len(bright_integrated), rmsd_integrated)
        )

    for crystal_model in experiments.crystals():
        if hasattr(crystal_model, "get_domain_size_ang"):
            log_str += ". Final ML model: domain size angstroms: {:f}, half mosaicity degrees: {:f}".format(
                crystal_model.get_domain_size_ang(),
                crystal_model.get_half_mosaicity_deg(),
            )

    #print0(log_str)
    return experiments, integrated




def dials_find_spots(data_img, params, trusted_flags=None):
    """
    :param data_img: numpy array image
    :param params: instance of stills_process params.spotfinder
    :param trusted_flags:
    :return:
    """
    if trusted_flags is None:
        trusted_flags = np.ones(data_img.shape, bool)
    thresh = SpotFinderFactory.configure_threshold(params)
    flex_data = flex.double(np.ascontiguousarray(data_img))
    flex_trusted_flags = flex.bool(np.ascontiguousarray(trusted_flags))
    spotmask = thresh.compute_threshold(flex_data, flex_trusted_flags)
    return spotmask.as_numpy_array()


def refls_from_sims(panel_imgs, detector, beam, thresh=0, filter=None, panel_ids=None,
                    max_spot_size=1000, phil_file=None, **kwargs):
    """
    This is for converting the centroids in the noiseless simtbx images
    to a multi panel reflection table
    :param panel_imgs: list or 3D array of detector panel simulations
    :param detector: dxtbx  detector model of a caspad
    :param beam:  dxtxb beam model
    :param thresh: threshol intensity for labeling centroids
    :param filter: optional filter to apply to images before
        labeling threshold, typically one of scipy.ndimage's filters
    :param pids: panel IDS , else assumes panel_imgs is same length as detector
    :param kwargs: kwargs to pass along to the optional filter
    :return: a reflection table of spot centroids
    """
    if panel_ids is None:
        panel_ids = np.arange(len(detector))
    pxlst_labs = []
    badpix_all =None
    min_spot_size=1
    if phil_file is not None:
        params = stills_process_params_from_file(phil_file)
        min_spot_size = params.spotfinder.filter.min_spot_size
    for i, pid in enumerate(panel_ids):
        plab = PixelListLabeller()
        img = panel_imgs[i]
        if phil_file is not None:
            params = stills_process_params_from_file(phil_file)
            badpix = None
            if params.spotfinder.lookup.mask is not None:
                if badpix_all is None:
                    badpix_all = utils.load_mask(params.spotfinder.lookup.mask)
                badpix = badpix_all[pid]
            mask = dials_find_spots(img, params, badpix)
        elif filter is not None:
            mask = filter(img, **kwargs) > thresh
        else:
            mask = img > thresh
        img_sz = detector[int(pid)].get_image_size()  # for some reason the int cast is necessary in Py3
        flex_img = flex.double(img)
        flex_img.reshape(flex.grid(img_sz))

        flex_mask = flex.bool(mask)
        flex_mask.resize(flex.grid(img_sz))
        pl = PixelList(0, flex.double(img), flex.bool(mask))
        plab.add(pl)

        pxlst_labs.append(plab)

    El = utils.explist_from_numpyarrays(panel_imgs, detector, beam)
    iset = El.imagesets()[0]
    refls = pixel_list_to_reflection_table(
        iset, pxlst_labs,
        min_spot_size=min_spot_size,
        max_spot_size=max_spot_size,  # TODO: change this ?
        filter_spots=FilterRunner(),  # must use a dummie filter runner!
        write_hot_pixel_mask=False)[0]
    if phil_file is not None:
        x,y,z = refls['xyzobs.px.value'].parts()
        x -=0.5
        y -=0.5
        refls['xyzobs.px.value'] = flex.vec3_double(x,y,z)

    return refls


if __name__=="__main__":

    if COMM.rank==0:
        if not os.path.exists( args.outdir):
            os.makedirs(args.outdir)
    COMM.barrier()

    #rank_outdir = os.path.join( args.outdir, "rank%d" % COMM.rank)
    #if not os.path.exists(rank_outdir):
    #    os.makedirs(rank_outdir)

    params = utils.get_extracted_params_from_phil_sources(args.predPhil, args.cmdlinePhil)

    # inputGlob can be a glob in strings, a single pandas file, or a hopper output folder
    if os.path.isfile(args.inputGlob) or os.path.isdir(args.inputGlob):
        if os.path.isfile(args.inputGlob):
            fnames = [args.inputGlob]
        else:
            dirname = args.inputGlob
            fnames = glob.glob( os.path.join(dirname, "pandas/hopper_results_rank*.pkl"))
    else:
        fnames = glob.glob(args.inputGlob)

    if not fnames:
        raise OSError("Found no filenames to load!")
    Nf = 0
    shots_per_df = []
    print0("getting total number of shots")
    for i_f, f in enumerate(fnames):
        if i_f % COMM.size != COMM.rank:
            continue
        n = len(pandas.read_pickle(f))
        shots_per_df += [(f, str(x)) for x in range(n)]  # note we cast to string because of mpi reduce
        Nf += n

    shots_per_df = COMM.bcast(COMM.reduce( shots_per_df))
    Nf = COMM.bcast(COMM.reduce(Nf))
    print0("total num shots is %d" % Nf)
    df_rows_per_rank = np.array_split(shots_per_df, COMM.size)[COMM.rank]

    print0("getting dataframe handles")
    df_handles = {}
    dfs = []
    if df_rows_per_rank.size:
        dfs, _ = zip(*df_rows_per_rank)
    for f in set(dfs):
        df_handles[f] = pandas.read_pickle(f).reset_index(drop=True)

    def df_iter():
        for i_df, (df_name, row_idx) in enumerate(df_rows_per_rank):
            row_idx = int(row_idx)
            printR("Opening shot %d / %d" % (i_df+1, len(df_rows_per_rank)))
            df = df_handles[df_name].iloc[row_idx: row_idx+1].copy()
            yield i_df, df

    if params.predictions.verbose:
        params.predictions.verbose = COMM.rank==0

    dev = COMM.rank % args.numdev

    EXPT_DIRS = os.path.join(args.outdir, "expts_and_refls")
    if COMM.rank==0:
        utils.safe_makedirs(EXPT_DIRS)

    if args.scanWeakFracs and params.predictions.weak_fraction != 1:
        print("WARNING: overriding weak_fracion because of scanWeakFracs")
        params.predictions.weak_fraction=1

    print0("Found %d input files" % Nf)
    with DeviceWrapper(dev) as _:
        all_dfs = []
        all_pred_names = []
        exp_ref_spec_lines = []
        all_rank_pred = None
        all_rank_expt = None

        rank_shot_count = 0
        rank_pred_file = os.path.join(EXPT_DIRS, "rank%d_preds.refl" % COMM.rank)
        rank_pred_file = os.path.abspath(rank_pred_file)
        rank_expt_file = rank_pred_file.replace(".refl", ".expt")
        for i_f, df in df_iter():

            expt_name = df.exp_name.values[0]
            expt_idx = df.exp_idx.values[0]
            tag = os.path.splitext(os.path.basename(expt_name))[0]

            data_expt = hopper_utils.DataModeler.exper_json_single_file(expt_name, expt_idx)
            data_exptList = ExperimentList()
            data_exptList.append(data_expt)

            try:
                spectrum_override = None
                if params.spectrum_from_imageset:
                    spectrum_override = downsamp_spec_from_params(params, data_expt)
                pred = predictions.get_predicted_from_pandas(
                    df, params, strong=None, device_Id=dev, spectrum_override=spectrum_override)
                if args.filterDupes:
                    pred = filter_refls(pred)
            except ValueError:
                #os.remove(new_expt_name)
                continue

            data = utils.image_data_from_expt(data_expt)
            Rstrong = refls_from_sims(data, data_expt.detector, data_expt.beam, phil_file=args.procPhil )
            Rstrong['id'] = flex.int(len(Rstrong), 0)
            num_panels = len(data_expt.detector)
            if num_panels > 1:
                assert params.predictions.label_weak_col == "rlp"

            Rstrong.centroid_px_to_mm(data_exptList)
            Rstrong.map_centroids_to_reciprocal_space(data_exptList)
            predictions.label_weak_predictions(pred, Rstrong, q_cutoff=params.predictions.qcut, col=params.predictions.label_weak_col )

            pred['is_strong'] = flex.bool(np.logical_not(pred['is_weak']))
            strong_sel = np.logical_not(pred['is_weak'])

            pred["refl_idx"] = flex.int(np.arange(len(pred)))

            #weaks = pred.select(pred['is_weak'])
            #weaks_sorted = np.argsort(weaks["scatter"])[::-1]
            #nweak = len(weaks)
            #num_keep = int(nweak*params.predictions.weak_fraction)
            #weak_refl_inds_keep = set(np.array(weaks["refl_idx"])[weaks_sorted[:num_keep]])
            #weak_sel = flex.bool([i in weak_refl_inds_keep for i in pred['refl_idx']])
            #keeps = np.logical_or( pred['is_strong'], weak_sel)
            #printR("Sum keeps=%d; num_strong=%d, num_kept_weak=%d" % (sum(keeps), sum(strong_sel), sum(weak_sel)))
            #pred = pred.select(flex.bool(keeps))
            pred = filter_weak_reflections(pred, weak_fraction=params.predictions.weak_fraction)

            nstrong = np.sum(strong_sel)
            printR("Will save %d refls (%d strong, %d weak)" % (len(pred), np.sum(pred["is_strong"]), np.sum(pred["is_weak"])))
            pred['id'] = flex.int(len(pred), rank_shot_count)
            if 'shoebox' in list(pred) and not args.keepShoeboxes:
                del pred['shoebox']
            if all_rank_pred is None:
                all_rank_pred = deepcopy(pred)
            else:
                all_rank_pred.extend(pred)

            # Note, the simple append causes memory leak:
            #all_rank_expt.append(data_expt)
            if all_rank_expt is None:
                all_rank_expt = deepcopy(data_exptList.to_dict())
            else:
                Edict = data_exptList.to_dict()
                for exp_key in 'beam', 'detector', 'crystal', 'imageset':
                    Edict['experiment'][0][exp_key] = rank_shot_count
                for exp_key in 'experiment', 'beam', 'detector', 'crystal', 'imageset':
                    assert len( Edict[exp_key])==1
                    all_rank_expt[exp_key] .append(Edict[exp_key][0])

            Rindexed = Rstrong.select(Rstrong['indexed'])
            if len(Rindexed)==0:
                print("No strong indexed refls for shot %s" % expt_name)
                continue

            utils.refls_to_hkl(Rindexed, data_expt.detector, data_expt.beam, data_expt.crystal, update_table=True)
            if args.dialsInteg:
                # TODO: save these files as multi-shot experiment/refls
                try:
                    int_expt, int_refl = integrate(args.procPhil, data_exptList, Rindexed, pred)
                    int_expt_name = "%s/%s_%d_integrated.expt" % (rank_outdir, tag, i_f)
                    int_expt.as_file(int_expt_name)
                    int_refl['bbox'] = int_refl['shoebox'].bounding_boxes()
                    int_refl_name = int_expt_name.replace(".expt", ".refl")
                    int_refl.as_file(int_refl_name)
                except RuntimeError:
                    print("Integration failed" )

            df['old_exp_name'] = expt_name
            df['old_exp_idx'] = expt_idx
            df['exp_name'] = rank_expt_file
            df['exp_idx'] = rank_shot_count

            df['predictions'] = rank_pred_file
            df['predicted_refs'] = rank_pred_file
            df['num_pred'] = len(pred)

            all_dfs.append(df)
            rank_shot_count += 1

            spec_name = df.spectrum_filename.values[0]
            if spec_name is None:
                spec_name = ""
            exp_ref_spec_lines.append("%s %s %s %d\n" % (rank_expt_file, rank_pred_file, spec_name, rank_shot_count))

        all_rank_pred.as_file(rank_pred_file)
        # NOTE: all_rank_expt is a dictionary to avoid weird OOM, so we write a simple json
        #all_rank_expt.as_file(rank_expt_file)
        with open(rank_expt_file, "w") as file_O:
            json.dump(all_rank_expt, file_O)
        print0("Done with predictions, combining dataframes")
        if all_dfs:
            all_dfs = pandas.concat(all_dfs)
        else:
            all_dfs = None

        if args.scanWeakFracs and all_dfs is not None:
            assert len(all_dfs.predictions.unique()) == 1
            pred_file = all_dfs.predictions.values[0]
            n_total_weak = np.sum(all_rank_pred['is_weak'])
            n_total = len(all_rank_pred)
            weak_fracs = [.11,.22,.33,.44,.55,.66,.77,.88]
            labels = []
            for i_frac, weak_frac in enumerate(weak_fracs):
                filt_refls = filter_weak_reflections(all_rank_pred, weak_frac)
                label="%dperc"%(weak_frac*100,)
                labels.append(label)
                new_pred_file = os.path.splitext(pred_file)[0]+"_%s.refl" % label
                all_dfs['predictions_%s' % label] = new_pred_file
                all_dfs['predicted_refs_%s' %label] = new_pred_file
                num_preds = []
                for exp_id in all_dfs.exp_idx.values:
                    n = np.sum(filt_refls['id'] == int(exp_id))
                    num_preds.append(n)
                all_dfs['num_pred_%s' %label] = num_preds
                filt_refls.as_file(new_pred_file)
                printR("Saved %d/%d refls (%d strong, %d/%d weak) to %s"
                       % (len(filt_refls), n_total, np.sum(filt_refls["is_strong"]), np.sum(filt_refls["is_weak"]), n_total_weak, new_pred_file))
            # note this sanity check below requires that weaK_fracs be sorted
            if sorted(weak_fracs) == weak_fracs:
                # then as weak frac increases, there should be an increasing number of predictions
                num_preds_per_frac = [all_dfs["num_pred_%s" % lab].sum() for lab in labels]
                assert num_preds_per_frac == sorted(num_preds_per_frac)

        print0("MPI gather all_dfs")
        all_dfs = COMM.gather(all_dfs)
        print0("MPI reduce lines")
        exp_ref_spec_lines = COMM.reduce(exp_ref_spec_lines)
        if COMM.rank==0:
            hopper_input_name = os.path.abspath(os.path.join(args.outdir , "%s.txt" % args.hopInputName))
            o = open(hopper_input_name, "w")
            for l in exp_ref_spec_lines:
                o.write(l)
            o.close()
            all_dfs = [df for df in all_dfs if df is not None]
            if not all_dfs:
                raise ValueError("No dataframes to concat: prediction/integration failed for all shots..")

            print("Concat frames")
            all_dfs = pandas.concat([df for df in all_dfs if df is not None])
            all_dfs.reset_index(inplace=True, drop=True)
            best_pkl_name = os.path.abspath(os.path.join(args.outdir , "%s.pkl" % args.hopInputName))
            all_dfs.to_pickle(best_pkl_name)
            print("Wrote %s (best_pickle option for simtbx.diffBragg.hopper) and %s (exp_ref_spec option for simtbx.diffBragg.hopper). Use them to run the predictions through hopper (use phil centroid=cal) or simtbx.diffBragg.stage_two." % (best_pkl_name, hopper_input_name))

            cmd_log_file = os.path.join(args.outdir, "cmdline_execution.txt")
            with open(cmd_log_file, "w") as o:
                o.write("integrate was run from folder: %s\n" % os.getcwd())
                o.write("The command line input was:\n")
                o.write(" ".join(sys.argv) + "\n")
                #TODO: write the diff phils here:
