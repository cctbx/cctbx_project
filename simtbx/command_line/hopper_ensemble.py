from __future__ import absolute_import, division, print_function

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("input", type=str, help="combined pandas pickle")
parser.add_argument("phil", type=str, help="user phil file used to run hopper (see simtbx/diffBragg/phil.py)")
parser.add_argument("--outdir", type=str, default=None, help="output folder")
parser.add_argument("--exp", type=str, default="opt_exp_name")
parser.add_argument("--refl", type=str, default="predicted_refls")
parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
parser.add_argument("--cell", nargs=6, type=float, default=None, help="unit cell to use when writing MTZ files. If not provided, average will be used")
parser.add_argument("--maxSigma", type=float, default=1e20, help="Fhkls are written to MTZ only if they have a sigma < than this value(default=1e20)")
parser.add_argument("--saveFreq", type=int, default=None, help="save an mtz each N iterations (default is None, i.e. only save after the last iteration)")
args = parser.parse_args()

import socket
import logging
import pandas


from libtbx.mpi4py import MPI
from simtbx.diffBragg.prep_stage2_input import prep_dataframe
from simtbx.diffBragg import hopper_ensemble_utils, hopper_utils
from dials.array_family import flex
from dxtbx.model import ExperimentList
from xfel.merging.application.utils.memory_usage import get_memory_usage


COMM= MPI.COMM_WORLD
LOGGER = logging.getLogger("diffBragg.main")


def mem_usage(rank):
    if COMM.rank == rank:
        memMB = get_memory_usage()
        host = socket.gethostname()
        LOGGER.info("Rank %d reporting memory usage: %f GB on Rank 0 node %s" % (COMM.rank, memMB / 1e3, host))


def load_inputs(pandas_table, params, exper_key="exp_name", refls_key='predictions'):

    print(refls_key)
    work_distribution = prep_dataframe(pandas_table, refls_key)
    COMM.Barrier()
    num_exp = len(pandas_table)
    first_exper_file = pandas_table[exper_key].values[0]
    detector = ExperimentList.from_file(first_exper_file, check_format=False)[0].detector
    if detector is None and params.refiner.reference_geom is None:
        raise RuntimeError("No detector in experiment, must provide a reference geom.")
    # TODO verify all shots have the same detector ?
    if params.refiner.reference_geom is not None:
        detector = ExperimentList.from_file(params.refiner.reference_geom, check_format=False)[
            0].detector
        LOGGER.debug("Using reference geom from expt %s" % params.refiner.reference_geom)

    if COMM.size > num_exp:
        raise ValueError("Requested %d MPI ranks to process %d shots. Reduce number of ranks to %d"
                         % (COMM.size, num_exp, num_exp))

    exper_names = pandas_table[exper_key]
    assert len(exper_names) == len(set(exper_names))
    worklist = work_distribution[COMM.rank]
    LOGGER.info("EVENT: begin loading inputs")

    shot_modelers = hopper_ensemble_utils.DataModelers()
    for ii, i_exp in enumerate(worklist):
        exper_name = exper_names[i_exp]
        LOGGER.info("EVENT: BEGIN loading experiment list")
        expt_list = ExperimentList.from_file(exper_name, check_format=params.refiner.check_expt_format)
        LOGGER.info("EVENT: DONE loading experiment list")
        if len(expt_list) != 1:
            LOGGER.critical("Input experiments need to have length 1, %s does not" % exper_name)
        expt = expt_list[0]
        expt.detector = detector  # in case of supplied ref geom

        exper_dataframe = pandas_table.query("%s=='%s'" % (exper_key, exper_name))

        refl_name = exper_dataframe[refls_key].values[0]
        refls = flex.reflection_table.from_file(refl_name)
        # FIXME need to remove (0,0,0) bboxes

        good_sel = flex.bool([h != (0, 0, 0) for h in list(refls["miller_index"])])
        refls = refls.select(good_sel)

        exp_cry_sym = expt.crystal.get_space_group().type().lookup_symbol()
        # TODO:replace here a cctbx method (is similar symmetry or something)
        if exp_cry_sym.replace(" ", "") != params.space_group:
            raise ValueError("Crystals should all have the same space group symmetry")

        LOGGER.info("EVENT: LOADING ROI DATA")
        shot_modeler = hopper_utils.DataModeler(params)
        shot_modeler.exper_name = exper_name
        shot_modeler.refl_name = refl_name
        shot_modeler.rank = COMM.rank
        if params.refiner.load_data_from_refl:
            gathered = shot_modeler.GatherFromReflectionTable(expt, refls, sg_symbol=params.space_groups)
        else:
            gathered = shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=params.space_group)
        if not gathered:
            raise IOError("Failed to gather data from experiment %s", exper_name)
        LOGGER.info("EVENT: DONE LOADING ROI")

        shot_modeler.set_parameters_for_experiment(best=exper_dataframe)
        shot_modeler.set_spectrum()
        LOGGER.info("Will simulate %d energy channels" % len(shot_modeler.nanoBragg_beam_spectrum))

        shot_modeler.Umatrices = [shot_modeler.E.crystal.get_U()]

        mem_usage(0)
        if COMM.rank==0:
            print("Finished loading image %d / %d" % (ii + 1, len(exper_names)), flush=True)

        shot_modelers.add_modeler(shot_modeler)

    shot_modelers.mpi_set_x_slices()

    assert shot_modelers.num_modelers > 0

    # use the first shot modeler to create a sim data instance:
    shot_modelers.SIM = hopper_utils.get_simulator_for_data_modelers(shot_modelers[0])
    return shot_modelers


if __name__ == "__main__":
    from simtbx.diffBragg.phil import philz, hopper_phil
    from libtbx.phil import parse


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
    # end of phil stuff ========

    from simtbx.diffBragg import mpi_logger

    if params.logging.disable:
        logging.disable(level=logging.CRITICAL)  # disables CRITICAL and below
    else:
        mpi_logger.setup_logging_from_params(params)

    df = pandas.read_pickle(args.input)

    if params.skip is not None:
        df = df.iloc[params.skip:]
    if params.first_n is not None:
        df = df.iloc[:params.first_n]
    df.reset_index(inplace=True, drop=True)
    modelers = load_inputs(df, params, exper_key=args.exp, refls_key=args.refl)
    modelers.cell_for_mtz = args.cell
    modelers.max_sigma = args.maxSigma
    modelers.outdir = args.outdir if args.outdir is not None else modelers.params.outdir
    modelers.save_freq = args.saveFreq
    modelers.prep_for_refinement()

    # do all sanity checks up front before minimization
    modelers.Minimize(save=True)

    LOGGER.debug("Done!")
