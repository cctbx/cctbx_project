from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME ens.hopper

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("input", type=str, help="combined pandas pickle")
parser.add_argument("phil", type=str, help="user phil file used to run hopper (see simtbx/diffBragg/phil.py)")
parser.add_argument("--outdir", type=str, default=None, help="output folder")
parser.add_argument("--exp", type=str, default="opt_exp_name", help="column name for input expeirments (default is opt_exp_name)")
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
import numpy as np
import socket
import logging
import pandas


from cctbx import sgtbx
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


def get_gather_name(exper_name, gather_dir):
    gathered_name = os.path.splitext(os.path.basename(exper_name))[0]
    gathered_name += "_withData.refl"
    gathered_name = os.path.join(gather_dir, gathered_name)
    return os.path.abspath(gathered_name)


def load_inputs(pandas_table, params, exper_key="exp_name", refls_key='predictions',
                gather_dir=None):

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
        if exp_cry_sym.replace(" ", "") != params.space_group:
            gr = sgtbx.space_group_info(params.space_group).group()
            expt.crystal.set_space_group(gr)
            #raise ValueError("Crystals should all have the same space group symmetry")

        LOGGER.info("EVENT: LOADING ROI DATA")
        shot_modeler = hopper_utils.DataModeler(params)
        shot_modeler.exper_name = exper_name
        shot_modeler.refl_name = refl_name
        shot_modeler.rank = COMM.rank
        if params.refiner.load_data_from_refl:
            gathered = shot_modeler.GatherFromReflectionTable(expt, refls, sg_symbol=params.space_group)
            LOGGER.debug("tried loading from reflection table")
        else:
            gathered = shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=params.space_group)
            LOGGER.debug("tried loading data from expt table")
        if not gathered:
            raise IOError("Failed to gather data from experiment %s", exper_name)
        else:
            LOGGER.debug("successfully loaded data")
        LOGGER.info("EVENT: DONE LOADING ROI")

        if gather_dir is not None:
            #gathered_name = os.path.splitext(os.path.basename(exper_name))[0]
            #gathered_name += "_withData.refl"
            #gathered_name = os.path.join(gather_dir, gathered_name)
            gathered_name = get_gather_name(exper_name, gather_dir)
            shot_modeler.dump_gathered_to_refl(gathered_name, do_xyobs_sanity_check=False)
            LOGGER.info("SAVED ROI DATA TO %s" % gathered_name)
            all_data = shot_modeler.all_data.copy()
            all_roi_id = shot_modeler.roi_id.copy()
            all_bg = shot_modeler.all_background.copy()
            all_trusted = shot_modeler.all_trusted.copy()
            all_pids = np.array(shot_modeler.pids)
            all_rois = np.array(shot_modeler.rois)
            new_Modeler = hopper_utils.DataModeler(params)
            assert new_Modeler.GatherFromReflectionTable(exper_name, gathered_name, sg_symbol=params.space_group)
            assert np.allclose(new_Modeler.all_data, all_data)
            assert np.allclose(new_Modeler.all_background, all_bg)
            assert np.allclose(new_Modeler.rois, all_rois)
            assert np.allclose(new_Modeler.pids, all_pids)
            assert np.allclose(new_Modeler.all_trusted, all_trusted)
            assert np.allclose(new_Modeler.roi_id, all_roi_id)
            LOGGER.info("Gathered file approved!")

        if gather_dir is not None:
            continue

        shot_modeler.set_parameters_for_experiment(best=exper_dataframe)
        shot_modeler.set_spectrum()
        LOGGER.info("Will simulate %d energy channels" % len(shot_modeler.nanoBragg_beam_spectrum))

        shot_modeler.Umatrices = [shot_modeler.E.crystal.get_U()]

        mem_usage(0)
        if COMM.rank==0:
            print("Finished loading image %d / %d" % (ii + 1, len(worklist)), flush=True)

        shot_modelers.add_modeler(shot_modeler)

    if gather_dir is not None:
        if COMM.rank==0:
            df['ens.hopper.imported'] = [get_gather_name(f_exp, gather_dir) for f_exp in pandas_table[exper_key]]
            pd_name = os.path.join(params.outdir, "preImport_for_ensemble.pkl")
            df.to_pickle(pd_name)
            print("Wrote file %s to be used to re-run ens.hopper . Use optional ens.hopper arg '--refl ens.hopper.imported', and the phil param load_data_from_refl=True to load the imported data" % pd_name)
        COMM.barrier()
        sys.exit()
    shot_modelers.mpi_set_x_slices()

    assert shot_modelers.num_modelers > 0

    # use the first shot modeler to create a sim data instance:
    shot_modelers.SIM = hopper_utils.get_simulator_for_data_modelers(shot_modelers[0])

    shot_modelers.set_Fhkl_channels()

    return shot_modelers


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
    # end of phil stuff ========

    if COMM.rank==0:
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)

        command_fname = os.path.join(params.outdir, "command_line_input.txt")
        with open(command_fname, "w") as o:
            o.write("Command line input:\n")
            o.write(" ".join(sys.argv) + "\n")


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

    gather_dir=None
    if args.preImport:
        assert args.outdir is not None
        gather_dir = os.path.join(args.outdir, "stage2_imported_data")
        if COMM.rank == 0:
            if not os.path.exists(gather_dir):
                os.makedirs(gather_dir)

    modelers = load_inputs(df, params, exper_key=args.exp, refls_key=args.refl, gather_dir=gather_dir)
    # note, we only go beyond this point if perImport flag was not passed
    modelers.cell_for_mtz = args.cell
    modelers.max_sigma = args.maxSigma
    modelers.outdir = args.outdir if args.outdir is not None else modelers.params.outdir
    modelers.save_freq = args.saveFreq
    modelers.prep_for_refinement()
    modelers.save_modeler_params = args.saveAll

    # do all sanity checks up front before minimization
    modelers.Minimize(save=True)

    LOGGER.debug("Done!")
