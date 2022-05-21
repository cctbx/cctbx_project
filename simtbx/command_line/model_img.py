from __future__ import division, print_function

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("mod", type=str, help="path to the data modeler file")
parser.add_argument("out", type=str, help="output file name")
parser.add_argument("--cudaDevice", type=int, default=None,
                    help="use this gpu device (if not provided, will not use CUDA and will instead use OpenMP)")
parser.add_argument("--plotSpec", action="store_true", help="plot the energy spectrum")
parser.add_argument("--compareData", action="store_true", help="include the data image in the plot")
parser.add_argument("--j", type=int, default=8,
                    help="number of processes to use for background extraction, only relevant if --compareData flag is present (default: 8)")
parser.add_argument("--filtsz", type=int, default=10,
                    help="the median filter used to extract backgroun will have this dimension, pixel units, higher values are slower (default: 10)")
parser.add_argument("--specFile", type=str, default=None,
                    help="Path to the spectrum .lam file. If None, then assume hopper_process output tree and try to autolocate ")
parser.add_argument("--pandaFile", type=str, default=None,
                    help="Path to the pandas model file. If None, then assume hopper_process output tree and try to autolocate ")

args = parser.parse_args()

import glob
import os
import pandas
import numpy as np
from joblib import Parallel, delayed

from simtbx.diffBragg import utils
from simtbx.nanoBragg.utils import H5AttributeGeomWriter
from simtbx.modeling import forward_models
from scipy.ndimage import median_filter
from dxtbx.model import ExperimentList

# LIBTBX_SET_DISPATCHER_NAME diffBragg.model_img


def extract_background_from_data(data, njobs, filt_shape=None):
    """
    :param data: numpy array (3d) of multi panel data
    :param njobs: number of processors to use
    :param filt_shape: shape of median filter used to extract background (default is 10,10, better is 20,20 but its slower)
    :return: background image, same shape as data
    """
    extracted_bg = np.zeros_like(data)
    if filt_shape is None:
        filt_shape = 10,10
    def get_bg(jid, njobs):
        bgs = {}
        for pid, p in enumerate(data):
            if pid % njobs != jid:
                continue
            bg = median_filter(p, filt_shape)
            if jid==0:
                print("processing panel %d / %d" % (pid, len(data)) )
            bgs[pid] = bg
        return bgs
    results = Parallel(n_jobs=njobs)(delayed(get_bg)(jid,njobs) for jid in range(njobs))
    for bgs in results:
        for pid in bgs:
            extracted_bg[pid]+= bgs[pid]
    return extracted_bg


def get_pandas_name_from_mod_name(mod_name):
    """searches the diffBragg.hopper_process output filder for the pandas pickle corresponding to the
    modeler file mod_name"""
    pd_glob = os.path.join(os.path.dirname(mod_name), "../pandas/rank*/*.pkl")
    fnames = glob.glob(pd_glob)
    basename = os.path.basename(mod_name).split("_indexed_modeler.npy")[0]
    pd_name = [f for f in fnames if basename in f]
    assert len(pd_name)==1
    pd_name = pd_name[0]
    return pd_name


if __name__=="__main__":


    # the spectrum file
    spec_file = args.specFile
    if spec_file is None:
        spec_file = args.mod.replace("_modeler.npy", "_spectrum.lam")

    M = np.load(args.mod, allow_pickle=True)[()]
    pd_name = args.pandaFile
    if pd_name is None:
        pd_name = get_pandas_name_from_mod_name(args.mod)
    spec = utils.load_spectra_file(spec_file, as_spectrum=True)
    if args.plotSpec:
        import pylab as plt
        from simtbx.nanoBragg.utils import ENERGY_CONV
        wave,wt = zip(*spec)
        en = ENERGY_CONV/ np.array(wave)
        plt.plot(en, wt)
        plt.xlabel("spectrum energy (ev)")
        plt.ylabel("spectrum fluence")
        plt.suptitle("Close to continue...")
        plt.show()

    print("Will load %s" % pd_name)
    df = pandas.read_pickle(pd_name)
    print(df.iloc[0].to_string())

    Fname = M.params.simulator.structure_factors.mtz_name
    Fcol = M.params.simulator.structure_factors.mtz_column
    no_det_thick = M.params.simulator.detector.force_zero_thickness
    print("Will use mtz file %s, column %s" % (Fname, Fcol))
    print("no detector thickness=%s" % str(no_det_thick))
    print("Loaded data frame:")
    print("Wil simulate %d energy channels loaded from %s (add --plotSpec to display the spectrum)" % (len(spec), spec_file))

    cuda = False
    dev = 0
    if args.cudaDevice is not None:
        cuda = True
        dev = args.cudaDevice
        print("Will use cuda GPU device %d" % dev)
    else:
        print("Will not use cuda, but will use openmp, control with env var OMP_NUM_THREADS")
    imgs, expt = forward_models.model_spots_from_pandas(
            df, mtz_file=Fname, mtz_col=Fcol,
            spectrum_override=spec, cuda=cuda, device_Id=dev,
            force_no_detector_thickness=no_det_thick,
            use_db=True)

    if args.compareData:
        print("Extracting data from experiment list %s..." % df.exp_name.values[0])
        El = ExperimentList.from_file(df.exp_name.values[0])
        iset = El[0].imageset
        # NOTE : assumes a multi-panel detector model, otherwise get_raw_data should have no arg, e.g. iset.get_raw_data()
        data = np.array([a.as_numpy_array() for a in iset.get_raw_data(0)])
        # divide the data by the adu factor
        data /= M.params.refiner.adu_per_photon
        # extract the background:
        print("Extracting background from data ...")
        bg = extract_background_from_data(data, args.j, (args.filtsz, args.filtsz))
        imgs_w_bg  = imgs + bg  # this image is the extracted background and the optimized forward Bragg model

        # this image inclues data except for those pixels that were modeled during stage 1
        imgs2 = data.copy()
        P = M.all_pid
        F = M.all_fast
        S = M.all_slow
        imgs2[P,S,F] = M.best_model + M.all_background

    comp_args = {"compression": "gzip", "compression_opts": 9}
    kwargs ={
        "filename": args.out,
        "image_shape": imgs.shape,
        "num_images": 5 if args.compareData else 1,
        "detector": expt.detector,
        "beam": expt.beam,
        "dtype": np.float32,
        "compression_args": comp_args,
    }

    print("Saving compressed output...")
    with H5AttributeGeomWriter(**kwargs) as writer:
        writer.add_image(imgs)
        if args.compareData:
            writer.add_image(bg)
            writer.add_image(imgs_w_bg)
            writer.add_image(data)
            writer.add_image(imgs2)
    print("Wrote %s" % args.out)
