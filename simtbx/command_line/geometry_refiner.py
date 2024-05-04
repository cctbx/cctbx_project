from __future__ import print_function, division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.geometry_refiner
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--phil", type=str, required=True, help="path to a phil string")
parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
progargs = parser.parse_args()

import sys

from libtbx.mpi4py import MPI
from libtbx.phil import parse
from simtbx.diffBragg.device import DeviceWrapper
from simtbx.diffBragg.phil import philz, hopper_phil

COMM = MPI.COMM_WORLD

phil_scope = parse(philz + hopper_phil)
arg_interp = phil_scope.command_line_argument_interpreter(home_scope="")

phil_file = open(progargs.phil, "r").read()
user_phil = parse(phil_file)
phil_sources = [user_phil]

if progargs.cmdlinePhil is not None:
    command_line_phils = [arg_interp.process(phil) for phil in progargs.cmdlinePhil]
    phil_sources += command_line_phils

working_phil, unused = phil_scope.fetch(sources=phil_sources, track_unused_definitions=True)
for loc in unused:
    print("WARNING: unused phil:", loc)
params = working_phil.extract()
device_Id = COMM.rank % params.refiner.num_devices
with DeviceWrapper(device_Id) as _:
    from simtbx.diffBragg.refiners.geometry import geom_min
    from simtbx.diffBragg.utils import find_diffBragg_instances

    if COMM.rank > 0:
        sys.tracebacklimit = 0
    geom_min(params)
    for name in find_diffBragg_instances(globals()): del globals()[name]
