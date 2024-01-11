from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME ens.hopper

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("input", type=str, help="combined pandas pickle")
parser.add_argument("phil", type=str, help="user phil file used to run hopper (see simtbx/diffBragg/phil.py)")
parser.add_argument("--outdir", type=str, default=None, help="output folder")
parser.add_argument("--exp", type=str, default="exp_name", help="column name for input expeirments (default is opt_exp_name)")
parser.add_argument("--refl", type=str, default="stage2_refls", help="column name for refls (default is stage2_refls)")
parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
parser.add_argument("--cell", nargs=6, type=float, default=None, help="unit cell to use when writing MTZ files. If not provided, average will be used")
parser.add_argument("--maxSigma", type=float, default=1e20, help="Fhkls are written to MTZ only if they have a sigma < than this value(default=1e20)")
parser.add_argument("--saveFreq", type=int, default=None, help="save an mtz each N iterations (default is None, i.e. only save after the last iteration)")
parser.add_argument("--preImport", action="store_true", help="convert the data to reflection table format, then exit. Subsequent runs will be quicker, and the data are not portable.")
parser.add_argument("--saveAll", action="store_true", help="save a pandas pickle for each modeled shot whenever an MTZ is written (--saveFreq and/or once at refinement termination)")
parser.add_argument("--saveTag", type=str, default="stage2", help="filename tag (only matters if --saveAll)")
args = parser.parse_args()

import os
import sys
import logging
import pandas

from simtbx.diffBragg.hopper_ensemble_utils import load_inputs
from libtbx.mpi4py import MPI
from simtbx.diffBragg.device import DeviceWrapper

COMM= MPI.COMM_WORLD
LOGGER = logging.getLogger("diffBragg.main")


def write_commandline(params):
    if COMM.rank==0:
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)

        command_fname = os.path.join(params.outdir, "command_line_input.txt")
        with open(command_fname, "w") as o:
            o.write("Command line input:\n")
            o.write(" ".join(sys.argv) + "\n")


if __name__ == "__main__":
    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz, hopper_phil
    from simtbx.diffBragg import mpi_logger

    # phil stuff ==========
    phil_scope = parse(philz+hopper_phil)
    arg_interp = phil_scope.command_line_argument_interpreter(home_scope="")

    phil_file = open(args.phil, "r").read()
    user_phil = parse(phil_file)
    phil_sources = [user_phil]

    if args.cmdlinePhil is not None:
        command_line_phils = [arg_interp.process(phil) for phil in args.cmdlinePhil]
        phil_sources += command_line_phils

    working_phil, unused = phil_scope.fetch(sources=phil_sources, track_unused_definitions=True)
    for loc in unused:
        print("WARNING: unused phil:", loc)
    params = working_phil.extract()
    if args.outdir is not None:
        params.outdir = args.outdir
    params.tag = args.saveTag
    if params.record_device_timings and COMM.rank > 0:
        params.record_device_timings = False  # only record for rank 0 otherwise there's too much output
    # end of phil stuff ========

    write_commandline(params)

    if params.logging.disable:
        logging.disable(level=logging.CRITICAL)  # disables CRITICAL and below
    else:
        mpi_logger.setup_logging_from_params(params)

    df = pandas.read_pickle(args.input)

    if params.skip is not None:
        df = df.iloc[params.skip:]
    if params.max_process is not None:
        df = df.iloc[:params.max_process]
    df.reset_index(inplace=True, drop=True)

    gather_dir=None
    if args.preImport:
        assert args.outdir is not None
        gather_dir = os.path.join(args.outdir, "stage2_imported_data")
        if COMM.rank == 0:
            if not os.path.exists(gather_dir):
                os.makedirs(gather_dir)

    for col in [args.exp, args.refl]:
        if col not in list(df):
            raise KeyError("Col %s is missing from dataframe" % col)

    modelers = load_inputs(df, params, exper_key=args.exp, refls_key=args.refl, gather_dir=gather_dir)
    # note, we only go beyond this point if perImport flag was not passed
    modelers.cell_for_mtz = args.cell
    modelers.max_sigma = args.maxSigma
    modelers.outdir = args.outdir if args.outdir is not None else modelers.params.outdir
    modelers.save_freq = args.saveFreq

    modelers.prep_for_refinement()

    with DeviceWrapper(modelers.SIM.D.device_Id) as _:
        modelers.alloc_max_pix_per_shot()
        modelers.save_modeler_params = args.saveAll

        # do all sanity checks up front before minimization
        modelers.Minimize(save=True)

        LOGGER.debug("Done!")
