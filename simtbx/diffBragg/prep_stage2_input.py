from __future__ import division, print_function

import pandas
import numpy as np
from dials.array_family import flex
from libtbx.mpi4py import MPI
from simtbx.diffBragg import utils
COMM = MPI.COMM_WORLD

import logging

LOGGER = logging.getLogger("diffBragg.main")

def get_equal_partition(weights, partitions):
    """
    Use Longest-processing-time-first to schedule the set of samples into partitions.

    :param weights: list of weights
    :param partitions: number of partitions
    :return: list of indices for each partition
    """
    distribution = [[] for _ in range(partitions)]
    load = np.zeros(partitions)
    descending = np.array(weights).argsort()[::-1]
    for idx in descending:
        lightest = load.argmin()
        distribution[lightest].append(idx)
        load[lightest] += weights[idx]
    return distribution

def prep_dataframe(df, refls_key="predictions", res_ranges_string=None):
    """

    :param df: input pandas dataframe for stage2
    :param refls_key: column in df containing the reflection filenames
    :param res_ranges_string: optional res_ranges_string phil param (params.refiner.res_ranges)
    :return:
    """
    # TODO make sure all pred files exist

    res_ranges = None
    if res_ranges_string is not None:
        res_ranges = utils.parse_reso_string(res_ranges_string)

    if refls_key not in list(df):
        raise KeyError("Dataframe has no key %s" % refls_key)
    nshots = len(df)
    df.reset_index(drop=True, inplace=True)
    df['index'] = df.index.values
    refl_info = df[["index", refls_key, "exp_idx"]].values

    # sort and split such that each rank will read many refls from same file
    sorted_names_and_ids = sorted(
        refl_info,
        key=lambda x: x[1])  # sort by name
    df_idx, refl_names, expt_ids = np.array_split(sorted_names_and_ids, COMM.size)[COMM.rank].T

    refls_per_shot = []
    if COMM.rank==0:
        LOGGER.info("Loading nrefls per shot")

    prev_name = ""  # keep track of the most recently read refl table file
    Rall = None
    for (i_shot, name, expt_id) in zip(df_idx, refl_names, expt_ids):
        if Rall is None or name != prev_name:
            Rall = flex.reflection_table.from_file(name)
            prev_name = name

        R = Rall.select(Rall['id'] == int(expt_id))
        if res_ranges is not None:
            num_ref = 0
            if 'rlp' not in set(R.keys()):
                raise KeyError("Cannot filter res ranges if rlp column not in refl tables")
            d = 1. / np.linalg.norm(R["rlp"], axis=1)  # resolution per refl
            for d_fine, d_coarse in res_ranges:
                d_sel = np.logical_and(d >= d_fine, d < d_coarse)
                num_ref += d_sel.sum()
        else:
            num_ref = len(R)

        if num_ref==0:
            LOGGER.critical("Reflection %s id=%d has 0 reflections !" % (name, expt_id, num_ref))
        refls_per_shot.append((i_shot, num_ref))

    refls_per_shot = COMM.reduce(refls_per_shot, root=0)
    work_distribution = None
    if COMM.rank==0:
        print("Found %d shots" % nshots)
        refls_per_shot = sorted(refls_per_shot)
        indices, weights = zip(*refls_per_shot)
        assert list(indices) == list(range(nshots)) # new sanity test

        work_distribution = get_equal_partition(weights, COMM.size)
    work_distribution = COMM.bcast(work_distribution, root=0)
    assert work_distribution is not None, "ERROR! Rank %d has no work_distribution!" % COMM.rank
    return work_distribution


if __name__ == "__main__":
    import sys
    df = pandas.read_pickle(sys.argv[1])
    new_df = prep_dataframe(df)
    new_df.to_pickle(sys.argv[2])
