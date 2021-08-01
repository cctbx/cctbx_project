from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--pkl", type=str, help="input pkl, output of stg1_view_mp.py")
parser.add_argument("--tag", type=str, default="expanded",help="tag name for expanded reflections")
parser.add_argument("--outputExpRef", type=str, help="name for output exper refls file to be used to re-run stage 1")
parser.add_argument("--maxdist", type=float, help="max distance to prediction in pixels", default=3)
parser.add_argument("--minsig", type=float, help="at least one model value needs to be above this threshold", default=-1)
parser.add_argument("--phil", type=str, default=[], nargs="+", help="CCTBX PHIL file containing find_spots and hopper phil params")
args = parser.parse_args()


#from dials.command_line.find_spots import working_phil
#params = working_phil.extract()
phil_str = """
include scope dials.command_line.find_spots.phil_scope
include scope simtbx.command_line.hopper.phil_scope
"""
from libtbx.phil import parse
master_phil = parse(phil_str, process_includes=True)

phil_sources =[]
for phil_file in args.phil:
    phil_sources.append(parse(open(phil_file, "r").read()))

working_phil = master_phil.fetch(sources=phil_sources)
PHIL_PARAMS = working_phil.extract()


import os
from dials.array_family import flex
from simtbx.diffBragg import hopper_utils
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
from scipy.spatial import distance
from itertools import groupby
from simtbx.diffBragg import utils
import numpy as np
from dxtbx.model.experiment_list import ExperimentListFactory
from iotbx.reflection_file_reader import any_reflection_file
import pandas
from copy import deepcopy


def get_pfs(E, rois, pids):
    """
    :param E: experimen -lrt *
    :param rois: list of rois, an roi is (x1,x2,y1,y2) where x is fast-scan, y is slow-scan
    :param pids:  list of panel ids, same length as rois
    :return: pfs object for diffBragg, and roi_id . pfs is 3 times the length of roi_id
    """
    npan = len(E.detector)
    nfast, nslow = E.detector[0].get_image_size()

    MASK = np.zeros((npan, nslow, nfast), bool)
    ROI_ID = np.zeros((npan, nslow, nfast), 'uint16')
    nspots = len(rois)

    for i_spot in range(nspots):
        x1, x2, y1, y2 = rois[i_spot]
        pid = pids[i_spot]
        MASK[pid, y1:y2, x1:x2] = True
        ROI_ID[pid, y1:y2, x1:x2] = i_spot
    p, s, f = np.where(MASK)
    roi_id = ROI_ID[p, s, f]
    pan_fast_slow = np.ascontiguousarray((np.vstack([p, f, s]).T).ravel())
    pan_fast_slow = flex.size_t(pan_fast_slow)
    return pan_fast_slow, roi_id


PHIL_PARAMS.spotfinder.threshold.dispersion.gain=9.481
PHIL_PARAMS.spotfinder.threshold.dispersion.sigma_background=1
PHIL_PARAMS.spotfinder.threshold.dispersion.sigma_strong=1.2
PHIL_PARAMS.spotfinder.threshold.dispersion.kernel_size=[3,3]
PHIL_PARAMS.spotfinder.threshold.dispersion.global_threshold=10
PHIL_PARAMS.spotfinder.filter.min_spot_size=4
PHIL_PARAMS.spotfinder.lookup.mask=PHIL_PARAMS.refiner.mask


REF_NAME = PHIL_PARAMS.simulator.structure_factors.mtz_name
REF_COL = PHIL_PARAMS.simulator.structure_factors.mtz_column
REF = any_reflection_file(REF_NAME).as_miller_arrays()[0]
REF_IDX = set(REF.indices())


deltaQ = 0.015
LOST_SCALE = 15  # TODO eliminate this garbage
cutoff_dist = args.maxdist
quiet = True

nold = 0
n = 0
shots_with_nmiss = 0
nmiss = 0
df = pandas.read_pickle(args.pkl)

# get the average unit cell used
ucell_params = df.a.mean(), df.b.mean(), df.c.mean(), df.al.mean(), df.be.mean(), df.ga.mean()
ucell_man = utils.manager_from_params(ucell_params)

output_df = []
for i_exp_name, exp_name in enumerate(df.opt_exp_name.values):
    if i_exp_name % COMM.size != COMM.rank:
        continue
    df_exp = df.query("opt_exp_name=='%s'" % exp_name)
    exp_name = df_exp.opt_exp_name.values[0]
    print("Loading expt")
    El = ExperimentListFactory.from_json_file(exp_name, True)
    expt = El[0]
    print("Done Loading expt")
    strong_reflections = flex.reflection_table.from_observations(El, PHIL_PARAMS)
    strong_refl_name = exp_name.replace(".expt", "_strong.refl")
    strong_reflections.as_file(strong_refl_name)
    # TODO shfit detZ and strip thickness
    if "detz_shift_mm" in list(df):  # NOTE, this could also be inside expt_name directly
        expt.detector = utils.shift_panelZ(expt.detector, df.detz_shift_mm.values[0])

    if PHIL_PARAMS.simulator.detector.force_zero_thickness:
        expt.detector = utils.strip_thickness_from_detector(expt.detector)


    _=utils.refls_to_q(strong_reflections, expt.detector, expt.beam,update_table=True)

    print("Wrote %d refls for exper pair\n%s %s" %(len(strong_reflections), strong_refl_name, exp_name))
    rois, _ = utils.get_roi_deltaQ(strong_reflections, deltaQ, expt)

    img_data = utils.image_data_from_expt(expt)

    spectrum = None
    if PHIL_PARAMS.spectrum_from_imageset:
        spectrum = hopper_utils.downsamp_spec_from_params(PHIL_PARAMS, expt)

    img_data /= PHIL_PARAMS.refiner.adu_per_photon
    is_trused = np.ones(img_data.shape, bool)
    hotpix_mask = None
    if PHIL_PARAMS.roi.hotpixel_mask is not None:
        is_trusted = utils.load_mask(PHIL_PARAMS.roi.hotpixel_mask)
        hotpix_mask = ~is_trusted
    sigma_rdout = PHIL_PARAMS.refiner.sigma_r / PHIL_PARAMS.refiner.adu_per_photon
    roi_packet = utils.get_roi_background_and_selection_flags(
        strong_reflections, img_data, shoebox_sz=PHIL_PARAMS.roi.shoebox_size,
        reject_edge_reflections=PHIL_PARAMS.roi.reject_edge_reflections,
        reject_roi_with_hotpix=PHIL_PARAMS.roi.reject_roi_with_hotpix,
        background_mask=None, hotpix_mask=hotpix_mask,
        bg_thresh=PHIL_PARAMS.roi.background_threshold,
        use_robust_estimation=not PHIL_PARAMS.roi.fit_tilt,
        set_negative_bg_to_zero=PHIL_PARAMS.roi.force_negative_background_to_zero,
        pad_for_background_estimation=PHIL_PARAMS.roi.pad_shoebox_for_background_estimation,
        sigma_rdout=sigma_rdout, deltaQ=PHIL_PARAMS.roi.deltaQ, experiment=expt,
        weighted_fit=PHIL_PARAMS.roi.fit_tilt_using_weights,
        tilt_relative_to_corner=PHIL_PARAMS.relative_tilt, ret_cov=True)


    rois, pids, tilt_abc, selection_flags, background, tilt_cov = roi_packet

    #pids = strong_reflections['panel']

    pids_rois = sorted(list(zip(pids, rois)), key=lambda x: x[0])
    gb = groupby(pids_rois, key=lambda x: x[0])
    rois_per_panel = {pid: [x[1] for x in list(v)] for pid, v in gb}
    print("Modeling!")
    model = utils.roi_spots_from_pandas(df_exp, rois_per_panel, quiet=quiet,
                                        mtz_file=REF_NAME, mtz_col=REF_COL, cuda=False,
                                        reset_Bmatrix=True,
                                        spectrum_override=spectrum,
                                        force_no_detector_thickness=True,
                                        norm_by_nsource=True)
    print("Done modeling")
    dists = []
    sel = []
    #reso = []
    #sigZ = []
    #ref_inds = []
    #hkls = []
    all_sigZ = []
    all_Z = []

    print("Finding Xcalcs")
    xycalcs = flex.vec3_double()
    for i_ref, (pid, roi) in enumerate(zip(pids, rois)):
        x1,x2,y1,y2 = roi
        img = model[pid][y1:y2, x1:x2]
        img *= LOST_SCALE
        if np.allclose(img,0):
            sel.append(False)
            xycalcs.append((0,0,0))
            dists.append(np.nan)
            all_Z.append(np.nan)
            all_sigZ.append(np.nan)
            continue
        elif np.all(img < args.minsig):
            sel.append(False)
            xycalcs.append((0,0,0))
            dists.append(np.nan)
            all_Z.append(np.nan)
            all_sigZ.append(np.nan)
            continue
        else:
            Y,X = np.indices(img.shape)
            Isum = img.sum()
            xcal = (X*img).sum() / Isum + x1+0.5
            ycal = (Y*img).sum() / Isum + y1+0.5
            xobs, yobs, _ = strong_reflections[i_ref]["xyzobs.px.value"]
            xycalcs.append((xcal, ycal, 0))
            dist = distance.euclidean((xcal, ycal),(xobs,yobs))
            dists.append(dist)
            if dist >cutoff_dist:
                sel.append(False)
                all_Z.append(np.nan)
                all_sigZ.append(np.nan)
                continue
            sel.append(True)
            data_img = img_data[pid,y1:y2, x1:x2 ]
            model_background = background[pid, y1:y2, x1:x2]
            diff_img = data_img - model_background - img
            noise_model = np.sqrt(data_img + sigma_rdout**2)

            Z = diff_img / noise_model
            Zfinite = Z[~np.isnan(Z)]
            sigZ = np.std( Zfinite)
            mnZ = np.mean( Zfinite)
            all_sigZ.append(sigZ)
            all_Z.append(mnZ)
            # HERE compute Z-score
    print("Done Finding Xcalcs")
    if not np.any(sel):
        continue

    strong_reflections["xyzcal.px"] = xycalcs
    strong_reflections["dists"] = flex.double(dists)
    strong_reflections["sigZ"] = flex.double(all_sigZ)
    strong_reflections["Z"] = flex.double(all_Z)
    indexed_refls = strong_reflections.select(flex.bool(sel))
    hkl, hkl_i = utils.refls_to_hkl(indexed_refls, expt.detector, expt.beam, expt.crystal,
                           update_table=True)

    # Shouldn't the hkl be set using the average unit cell ?
    Cave =deepcopy(expt.crystal)
    Cave.set_B(ucell_man.B_recipspace)
    hkl2, hkl_i2 = utils.refls_to_hkl(indexed_refls, expt.detector, expt.beam, Cave,
                                    update_table=True)
    sel_misindexed = np.ones(len(hkl_i), bool)
    if not np.all(hkl_i==hkl_i2):
        shots_with_nmiss+= 1
        for i_h, (h1, h2) in enumerate(zip(hkl_i, hkl_i2)):
            if not np.all(h1 == h2):
                nmiss += 1
                sel_misindexed[i_h] = False

    hkl_i_asu = utils.map_hkl_list(hkl_i2, True, "P6522")
    sel_hkl = np.array([tuple(hi) in REF_IDX for hi in hkl_i_asu])

    SEL = flex.bool(np.logical_and(sel_hkl, sel_misindexed))
    indexed_refls = indexed_refls.select(SEL)

    nidx = len(indexed_refls)
    if nidx==0:
        continue
    md = np.median(indexed_refls["dists"])
    out = "RANK %3d:  shot %d / %d\n------------------\nThe model predicted %d  / %d strong spots with median prediction offset of %f" \
          %(COMM.rank, i_exp_name+1, len(df), nidx, len(rois), md)
    # pred offsets here
    groupA = df_exp.pred_offsets.values[0]
    md_groupA = np.median(groupA)
    out += "\n\tAs comparison, groupA was %d spots, with median prediction offset of %f" %(len(groupA), md_groupA)
    out += "\n\tOriginal exper: %s" % df_exp.exp_name.values[0]

    #old_refl_name = df_exp.exp_name.values[0].replace(".expt",".refl")
    old_refl_name = df_exp.refl_names.values[0]
    oldR= flex.reflection_table.from_file(old_refl_name)
    nold+= len(oldR)
    n += nidx

    #idx_refl_name = exp_name.replace(".expt", "_indexed1.refl")
    idx_refl_name = exp_name.replace(".expt", "_%s.refl" % args.tag)
    indexed_refls.as_file(idx_refl_name)
    out += "\n\tNew indexed refls file: %s" % idx_refl_name
    print(out)

    df_exp.sigmaZ_PoissonDat = [tuple(indexed_refls["sigZ"])]

    hkl_i_asu = utils.map_hkl_list(list(indexed_refls["miller_index"]), True, "P6522")
    df_exp.hkl = [tuple(hkl_i_asu)]

    res = 1 / np.linalg.norm(indexed_refls["rlp"], axis=1)
    df_exp.resolution = [tuple(res)]

    df_exp.refls_idx = [tuple(range(len(indexed_refls)))]
    df_exp.stage1_refls = idx_refl_name
    output_df.append(df_exp)
    #
    #bragg_pix, _ = model(x, SIM, pfs, verbose=True, compute_grad=False):


nold = COMM.reduce(nold)
n = COMM.reduce(n)
nmiss = COMM.reduce(nmiss)
shots_with_nmiss = COMM.reduce(shots_with_nmiss)
all_output_df = COMM.reduce(output_df)
if COMM.rank==0:
    DF = pandas.concat(all_output_df)
    outpkl = os.path.splitext(args.pkl)[0]+"_%s.pkl" % args.tag
    DF.to_pickle(outpkl)
    print("Wrote new pkl %s" % outpkl)
    ofname = args.outputExpRef
    o = open(ofname, "w")
    predictions = []
    for e, eopt, s in zip(df.exp_name, df.opt_exp_name, df.spectrum_filename):
        r = eopt.replace(".expt", "_%s.refl" %args.tag)
        if not os.path.exists(r):  # TODO why dont some exist ?Because of continues in the above script?
            predictions.append(None)
            continue
        predictions.append(r)
        o.write("%s %s %s\n" % (e, r, s))
    df['predictions'] = predictions
    sel = [f is not None for f in df.predictions]
    df = df.loc[sel].reset_index()
    df.to_pickle(args.pkl+".pred")
    o.close()
    print("Saved exper ref file %s" % ofname)
    print("Nold=%d; Nnew=%d" % (nold, n) )
    print("Num mis-indexed by unit cell shift: %d" % nmiss)
    print("Numshots with mis indexed refls by unit cell shift: %d" % shots_with_nmiss)
    print("Average ucell used for indexing:", ucell_params)
    print("Done.")
