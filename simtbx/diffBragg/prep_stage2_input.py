from __future__ import division, print_function

import pandas
import numpy as np
from dials.array_family import flex
from libtbx.mpi4py import  MPI
COMM = MPI.COMM_WORLD
import time
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

def prep_dataframe(df):
    # TODO make sure all pred files exist
    nshots = len(df)
    refls_names = df.predictions
    refls_per_shot = []
    if COMM.rank==0:
        LOGGER.info("Loading nrefls per shot")
    for i_shot, name in enumerate(refls_names):
        if i_shot % COMM.size != COMM.rank:
            continue
        R = flex.reflection_table.from_file(name)
        if len(R)==0:
            LOGGER.critical("Reflection %s has 0 reflections !" % (name, len(R)))
        refls_per_shot.append((i_shot, len(R)))

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
