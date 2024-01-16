from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.make_input_file

from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter as show_defaults_formatter
parser = ArgumentParser(formatter_class=show_defaults_formatter)
parser.add_argument("dirnames", nargs="+", type=str, help="processing dirs")
parser.add_argument("filename", type=str, help="the name of the diffBragg input file written at the end of this script")
parser.add_argument("--splitDir", default=None, type=str, help="optional folder for writing split expts and refls. If None, then "
                                                               "split files will be written in same folders as their sources")
parser.add_argument("--exptSuffix", type=str, default="refined.expt", help="find experiments with this suffix")
parser.add_argument("--reflSuffix", type=str, default="indexed.refl", help="find reflection files with this suffix")
parser.add_argument("--write", action="store_true")


args = parser.parse_args()

import glob
import os
apath = os.path.abspath
import json

from dxtbx.model import ExperimentList
from dials.array_family import flex
from simtbx.diffBragg import hopper_io
import hashlib


from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

if COMM.rank==0 and args.write:
    if args.splitDir is not None and not os.path.exists(args.splitDir):
        os.makedirs(args.splitDir)


def hash_name(name):
    hash_obj = hashlib.md5(name.encode('utf-8'))
    return hash_obj.hexdigest()


def get_idx_path(El):
    """return the idx and path in a single-experiment experimentList"""
    assert len(El)==1
    iset = El[0].imageset
    assert len(iset)==1
    idx = iset.indices()[0]
    path = iset.paths()[0]
    return idx, path


def split_stills_expts(expt_f, refl_f, split_dir, write=False):
    El = ExperimentList.from_file(expt_f, False)
    R = flex.reflection_table.from_file(refl_f)
    expt_names = []
    refl_names = []
    orig_expt_names, orig_refl_names = [],[]  # store the original names for recordkeeping
    seen_isets = {}
    for i_expt in range(len(El)):
        one_exp_El = El[i_expt: i_expt+1]
        idx, path = get_idx_path(one_exp_El)
        iset_id = idx,path
        if iset_id not in seen_isets:
            seen_isets[iset_id] = 0
        else:
            seen_isets[iset_id] += 1
        tag = "%s-%d" % (os.path.basename(os.path.splitext(path)[0]), idx)
        new_expt_name = os.path.splitext(expt_f)[0] + "_%s_xtal%d.expt" % (tag, seen_isets[iset_id])
        if write and split_dir is not None:
            unique_tag = "shot_%s" % hash_name(new_expt_name) + ".expt"
            new_expt_name = os.path.join(split_dir, unique_tag)
        new_refl_name = new_expt_name.replace(".expt", ".refl")
        refls = R.select(R['id'] == i_expt)
        refls.reset_ids()

        if write:
            one_exp_El.as_file(new_expt_name)
            refls.as_file(new_refl_name)
        expt_names.append(new_expt_name)
        refl_names.append(new_refl_name)
        orig_expt_names.append((apath(new_expt_name), (apath(expt_f), i_expt)))
        orig_refl_names.append((apath(new_refl_name), (apath(refl_f), i_expt)))
    return expt_names, refl_names, orig_expt_names, orig_refl_names


exp_names, ref_names = [], []
orig_exp_names, orig_ref_names = [], []
for i_dir, dirname in enumerate(args.dirnames):
    expt_glob = os.path.join( dirname, "*%s" % args.exptSuffix)
    expt_fnames = glob.glob(expt_glob)
    if COMM.rank==0:
        print("dir %s, glob=%s, num files=%d" % (dirname, expt_glob, len(expt_fnames)))

    for i_f, f in enumerate(expt_fnames):
        if i_f % COMM.size != COMM.rank:
            continue
        if COMM.rank == 0:
            print("processing job %d / %d in dir %d / %d"
                  %(i_f+1, len(expt_fnames), i_dir+1, len(args.dirnames) ))
        ref_f = f.replace(args.exptSuffix, args.reflSuffix)
        if not os.path.exists(ref_f):
            raise FileNotFoundError("No matching refl file for expt %s" % f)
        El = ExperimentList.from_file(f, False)
        if len(El.imagesets()) > 1 or len(El.crystals()) > 1:
            exp_fs, ref_fs, orig_exp_fs, orig_ref_fs = split_stills_expts(f, ref_f, args.splitDir)
            exp_names += exp_fs
            ref_names += ref_fs
            orig_exp_names += orig_exp_fs
            orig_ref_names += orig_ref_fs
        else:
            exp_names.append(f)
            ref_names.append(ref_f)
            orig_exp_names.append((apath(f),    (apath(f),0))  )
            orig_ref_names.append((apath(ref_f),(apath(ref_f),0))    )

exp_names = COMM.reduce(exp_names)
ref_names = COMM.reduce(ref_names)
orig_exp_names = COMM.reduce(orig_exp_names)
orig_ref_names = COMM.reduce(orig_ref_names)


if COMM.rank==0:

    if args.write:
        print("Saving the input file for diffBragg")
        hopper_io.save_expt_refl_file(args.filename, exp_names, ref_names, check_exists=True)
        print("Saved %s" % args.filename)
        jname = args.filename + ".json"
        jdat = {"expt": dict(orig_exp_names), "refl": dict(orig_ref_names)}
        with open(jname, "w") as fp:
            json.dump(jdat, fp, indent=1)
        print("Wrote json %s, which maps the hashnames to the original expt files" % jname)
    else:
        _, exp_and_idx = zip(*orig_exp_names)
        _, ref_and_idx = zip(*orig_ref_names)
        exp, exp_idx = zip(*exp_and_idx)
        ref, ref_idx = zip(*ref_and_idx)
        assert exp_idx == ref_idx
        hopper_io.save_expt_refl_file(args.filename, exp, ref, check_exists=True,
                                      indices=exp_idx)
