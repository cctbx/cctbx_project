from __future__ import division
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
from dxtbx.model import DetectorFactory, BeamFactory
from simtbx.nanoBragg.tst_nanoBragg_multipanel import det_descr, beam_descr
from simtbx.diffBragg import diffBragg

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--ngpu", type=int, default=1)
parser.add_argument("--kokkos", action="store_true")
args = parser.parse_args()
NGPU_PER_NODE = args.ngpu
if args.kokkos:
    import os
    os.environ["DIFFBRAGG_USE_KOKKOS"] = "1"

GPU_ID = COMM.rank % NGPU_PER_NODE
beam = BeamFactory.from_dict(beam_descr)
det = DetectorFactory.from_dict(det_descr)

from simtbx.diffBragg.utils import find_diffBragg_instances
from simtbx.diffBragg.device import DeviceWrapper
with DeviceWrapper(GPU_ID) as wrapper:
    for i in range(COMM.size*3):
        if i % COMM.size != COMM.rank:
            continue
        D = diffBragg(det, beam)
        D.verbose = 1
        D.vectorize_umats()
        for _ in range(5):
            D.add_diffBragg_spots()

        D.free_all()
        D.free_Fhkl2()
        print("RANK %d, done with Iter %d" % (COMM.rank, i+1))

    # TODO: the following deletion of the diffBragg instance is necessary to prevent the deallocation error
    # is there someway around this ????
    for name in find_diffBragg_instances(globals()): del globals()[name]
    print("OK!")
