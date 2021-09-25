
from simtbx.modeling.forward_models import roi_spots_from_pandas
from simtbx.diffBragg.utils import refls_to_hkl, refls_from_sims
from dials.array_family import flex
from dxtbx.model import ExperimentList
from scipy.spatial import cKDTree, distance
from dials.algorithms.shoebox import MaskCode
SIGNAL_MASK = MaskCode.Valid + MaskCode.Foreground
import numpy as np
from numpy import logical_or as logi_or
from numpy import logical_and as logi_and
from numpy import logical_not as logi_not


def get_predicted_from_pandas(df, params, strong, eid, device_Id=0):
    """
    :param df: pandas dataframe, stage1_df attribute of simtbx.command_line.hopper_process.HopperProcess
    :param params: instance of diffBragg/phil.py phil params
    :param strong: strong (observed) reflections
    :param eid: experiment identifier
    :param device_Id: GPU device Id for simuting forward model
    :return: predicted reflections table , to be passed along to dials.integrate functions
    """

    # returns the images and the experiment including any pre-modeling modifications (e.g. thinning out the detector)
    panel_images, expt = roi_spots_from_pandas(
        df,
        oversample_override=params.predictions.oversample_override,
        Ncells_abc_override=params.predictions.Nabc_override,
        pink_stride_override=params.predictions.pink_stride_override,
        defaultF=params.predictions.default_Famplitude,
        device_Id=device_Id,
        d_max=params.predictions.resolution_range[1],
        d_min=params.predictions.resolution_range[0],
        symbol_override=params.predictions.symbol_override,
        force_no_detector_thickness=params.simulator.detector.force_zero_thickness,
        use_exascale_api=params.predictions.method == "exascale",
        use_db=params.predictions.method == "diffbragg")

    predictions = refls_from_sims(panel_images, expt.detector, expt.beam, thresh=params.predictions.threshold,
                                  max_spot_size=1000)

    print("Found %d Bragg peak predictions above the threshold")

    # TODO: pulled these from comparing to a normal stills_process prediction table, not sure what they imply
    # TODO: multiple experiments per shot
    predictions['flags'] = flex.size_t(len(predictions), 1)
    predictions['id'] = flex.int(len(predictions), 0)
    predictions['entering'] = flex.bool(len(predictions), False)
    predictions['delpsical.rad'] = flex.double(len(predictions), 0)
    if eid:
        predictions.experiment_identifiers()[0] = eid

    El = ExperimentList()
    El.append(expt)
    predictions.centroid_px_to_mm(El)
    predictions.map_centroids_to_reciprocal_space(El)

    strong.centroid_px_to_mm(El)
    strong.map_centroids_to_reciprocal_space(El)

    refls_to_hkl(predictions, expt.detector, expt.beam, expt.crystal, update_table=True)
    predictions['xyzcal.px'] = predictions['xyzobs.px.value']
    predictions['xyzcal.mm'] = predictions['xyzobs.mm.value']
    predictions["num_pixels"] = numpix = predictions["shoebox"].count_mask_values(SIGNAL_MASK)
    predictions['scatter'] = predictions["intensity.sum.value"] / flex.double(np.array(numpix, np.float64))

    # separate out the weak from the strong
    label_weak_predictions(predictions, strong)
    n_weak = sum(predictions["is_weak"])
    n_pred = len(predictions)
    n_strong = n_pred - n_weak
    print("%d / %d predicted refls are near strongs" % (n_strong, n_pred))

    label_weak_spots_for_integration(params.predictions.weak_fraction, predictions)
    print("Will use %d spots for integration" % sum(predictions["is_for_integration"]))
    predictions = predictions.select(predictions["is_for_integration"])

    return predictions


def label_weak_predictions(predictions, strong, q_cutoff=0.005):
    """
    :param predictions: model reflection table
    :param strong: strong observed spots (reflection table)
    :param q_cutoff: distance in RLP space - arbitrary, increasing will bring in more candidates, 0.005 seems reasonable
    """
    strong_tree = cKDTree(strong['rlp'])
    predicted_tree = cKDTree(predictions['rlp'])

    # for each strong refl, find all predictions within q_cutoff of the strong rlp
    pred_idx_candidates = strong_tree.query_ball_tree(predicted_tree, q_cutoff)

    is_weak = flex.bool(len(predictions), True)
    for i_idx, cands in enumerate(pred_idx_candidates):
        if not cands:
            print("WARNING: no predicted refl candidates for strong refl %d - consider reducing changing threshold arg or increasing default_Famplitude" % i_idx)
            continue
        if len(cands) == 1:
            # if 1 spot is within q_cutoff , then its the closest
            pred_idx = cands[0]
        else:
            # in this case there are multiple predictions near the strong refl, we choose the closest one
            dists = []
            for c in cands:
                d = distance.euclidean(strong_tree.data[i_idx], predicted_tree.data[c])
                dists.append(d)
            pred_idx = cands[np.argmin(dists)]
        is_weak[pred_idx] = False
    predictions["is_weak"] = is_weak


def label_weak_spots_for_integration(fraction, predictions, num_res_bins=10):
    """
    :param fraction: fraction of reflections to label as "integratable" within each resolution bin
    :param predictions: dials.flex.reflections_table
    :param num_res_bins: number of resolution bins
    """
    res = 1. / np.linalg.norm(predictions["rlp"], axis=1)
    res_sort = np.sort(res)
    res_bins = [rb[0]-1e-6 for rb in np.array_split( res_sort, num_res_bins)] + [res_sort[-1]+1e-6]
    res_bin_assigments = np.digitize(res, res_bins)
    is_weak_but_integratable = np.zeros(len(predictions)).astype(np.bool)
    for i_res in range(1, num_res_bins+1):
        # grab weak spots in this res bin
        is_weak_and_in_bin = logi_and(res_bin_assigments == i_res, predictions["is_weak"])
        refls_in_bin = predictions.select(flex.bool(is_weak_and_in_bin))

        # determine which weak spots are closer to the Ewald sphere, based on the scatter value
        signal_cutoff_in_bin = np.percentile(refls_in_bin['scatter'], fraction*100)
        is_above_cutoff = predictions["scatter"] > signal_cutoff_in_bin

        is_integratable = logi_and(is_weak_and_in_bin, is_above_cutoff)
        is_weak_but_integratable[is_integratable] = True
    predictions['is_for_integration'] = flex.bool(logi_or(logi_not(predictions["is_weak"]), is_weak_but_integratable))
