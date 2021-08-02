from __future__ import absolute_import, division, print_function

from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

import warnings
from copy import deepcopy
from simtbx.diffBragg.refiners.local_refiner import LocalRefiner
from simtbx.diffBragg.refiners import BreakBecauseSignal

import logging
LOGGER = logging.getLogger("main")

warnings.filterwarnings("ignore")

import time, signal
class MPISignalHandler:
    """help to exit from HPC jobs before timeout, requires knowing exit signal, for summit its 12 or SIGUSR2 last I checked"""
    def __init__(self):
        self.t = time.time()

    def handle(self, signum, frame):
        t = time.time()-self.t
        print("COMM.rank%d Recived signal " % comm.rank, signum," after program running for %f sec" % t)
        for i in range(comm.size):
            if i != comm.rank:
                comm.iend(signum, i, tag=42)
        raise BreakBecauseSignal


MPISIGHAND = MPISignalHandler()


class GlobalRefiner(LocalRefiner):

    def __init__(self, *args, **kwargs):
        super(GlobalRefiner, self).__init__(verbose=rank==0, *args, **kwargs)
        self.rank = rank
        self.I_AM_ROOT = rank == 0
        print("STARTINGREFINEMENT!!! Rank %d Sz %d"% (self.rank, comm.size))
        self._sig_hand = MPISIGHAND  # overwrite single rank version

    def _MPI_sync_hkl_freq(self):
        if self.refine_Fcell:
            if self.rank != 0:
                self.hkl_frequency = None
            self.hkl_frequency = comm.bcast(self.hkl_frequency)

    def _MPI_sync_fcell_parameters(self):
        if not self.I_AM_ROOT:
            self.sigma_for_res_id = None
            self.res_group_id_from_fcell_index = None
            self.fcell_init = None
            self.resolution_ids_from_i_fcell = self.fcell_sigmas_from_i_fcell = self.fcell_init_from_i_fcell = None

        if self.rescale_params:
            if self.refine_Fcell:
                self.fcell_init = comm.bcast(self.fcell_init)
                self.sigma_for_res_id = comm.bcast(self.sigma_for_res_id)
                self.res_group_id_from_fcell_index = comm.bcast(self.res_group_id_from_fcell_index)
                self.resolution_ids_from_i_fcell = comm.bcast(self.resolution_ids_from_i_fcell)
                self.fcell_sigmas_from_i_fcell = comm.bcast(self.fcell_sigmas_from_i_fcell)
                self.fcell_init_from_i_fcell = comm.bcast(self.fcell_init_from_i_fcell)

    def _MPI_sync_panel_params(self):
        if not self.I_AM_ROOT:
            self.panelRot_params = None
            self.panelX_params = None
            self.panelY_params = None
            self.panelZ_params = None
        self.panelRot_params = comm.bcast(self.panelRot_params)
        self.panelX_params = comm.bcast(self.panelX_params)
        self.panelY_params = comm.bcast(self.panelY_params)
        self.panelZ_params = comm.bcast(self.panelZ_params)

    def _data_for_write(self, parameter_dict):
        all_data = comm.gather(parameter_dict)
        return all_data

    def _MPI_aggregate_model_data_correlations(self):
        LOGGER.info("image corr")
        self.all_image_corr = [self._get_image_correlation(i) for i in self.shot_ids]
        if self.init_image_corr is None:
            self.init_image_corr = deepcopy(self.all_image_corr)
        self.all_image_corr = comm.reduce(self.all_image_corr, MPI.SUM, root=0)

    def _init_n_bad_shots(self):
        self.n_bad_shots = len(self.bad_shot_list)
        self.n_bad_shots = comm.bcast(self.n_bad_shots)
        return self.n_bad_shots

    def _init_gather_ang_off(self):
        all_ang_off = self._MPI_reduce_broadcast(list(self.all_ang_off)) #comm.reduce(self.all_ang_off, root=0)
        return all_ang_off

    def _get_ang_off(self):
        return self.all_ang_off

    def _MPI_reduce_broadcast(self, var):
        var = comm.reduce(var, MPI.SUM, root=0)
        var = comm.bcast(var, root=0)
        return var

    def _MPI_combine_data_to_send(self):
        data_to_send = comm.reduce(self.data_to_send, MPI.SUM, root=0)
        return data_to_send

    def _MPI_barrier(self):
        comm.barrier()

    def _MPI_check_for_break_signal(self):
        if self.break_signal is not None:
            for i in range(comm.size):
                if i==comm.rank:
                    continue
                req = comm.irecv(source=i, tag=42)
                received_signal, signal = req.test()
                if received_signal and signal == self.break_signal:
                    raise BreakBecauseSignal
