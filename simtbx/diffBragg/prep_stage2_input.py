
import pandas
import numpy as np
from dials.array_family import flex
from libtbx.mpi4py import  MPI
COMM = MPI.COMM_WORLD
import time
import logging

LOGGER = logging.getLogger("main")


def prep_dataframe(df, time_to_try=60., print_time_interval=2):
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
            LOGGER.critical("REflection %s has 0 reflections !" % (name, len(R)))
        refls_per_shot.append((i_shot, len(R)))

    refls_per_shot = COMM.reduce(refls_per_shot)
    if COMM.rank==0:
        refls_per_shot = {i_shot: n for i_shot, n in refls_per_shot}
    refls_per_shot = COMM.bcast(refls_per_shot)

    nloaded_sig = np.inf
    best_order = None
    shot_order = np.arange(nshots)
    #shots_per_rank = {rank: shots for rank, shots in enumerate(np.array_split(shot_order, COMM.size))}
    best_shots_per_rank = None
    if COMM.rank==0:
        print("Found %d shots" % nshots)
        tstart = time.time()
        tlast_print = 0
        while 1:
            nloaded = [0]*COMM.size
            for ii, i_shot in enumerate(shot_order):
                rank = ii%COMM.size
                nloaded[rank] += refls_per_shot[i_shot]
                #if ii % COMM.size != COMM.rank:
                #    continue
                #nloaded += refls_per_shot[i_shot]
            #nloaded = COMM.gather(nloaded)
            #if COMM.rank == 0:
            this_nloaded_sig = np.std(nloaded)
            if this_nloaded_sig < nloaded_sig:
                best_order = shot_order
                #best_shots_per_rank = shots_per_rank
                nloaded_sig = this_nloaded_sig

            shot_order = np.random.permutation(nshots)
            #shots_per_rank = {rank:[] for rank in range(COMM.size)}
            #ranks_assign = np.random.randint(0, COMM.size, nshots )
            #for i_shot, rank in enumerate(ranks_assign):
            #    shots_per_rank[rank].append(i_shot)

            #best_order = COMM.bcast(best_order)
            #nloaded_sig = COMM.bcast(nloaded_sig)
            #shot_order = COMM.bcast(shot_order)
            #shots_per_rank = COMM.bcast(shots_per_rank)
            time_passed = time.time() - tstart
            #if COMM.rank==0 and time.time()-tlast_print > print_time_interval:
            if time.time()-tlast_print > print_time_interval:
                print("nloaded sigma=%.2f refls, time passed= %.2f sec" % (nloaded_sig, time_passed), flush=True)
                tlast_print = time.time()
            if time_passed > time_to_try:
                break

    #return best_shots_per_rank
    best_order = COMM.bcast(best_order)
    # nloaded_sig = COMM.bcast(nloaded_sig)
    assert best_order is not None
    df = df.iloc[best_order]
    return df


if __name__ == "__main__":
    import sys
    df = pandas.read_pickle(sys.argv[1])
    new_df = prep_dataframe(df, time_to_try=float(sys.argv[2]))
    new_df.to_pickle(sys.argv[3])
