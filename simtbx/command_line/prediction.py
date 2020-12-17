from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.prediction

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
import pandas
from simtbx.diffBragg import utils
import time
from libtbx.phil import parse
import os


help_message = "predictions using diffBragg refinement results"

script_phil = """
refl_list = None
  .type = str
  .help = refl file list
exper_list = None
  .type = str
  .help = expt file list
panda_list = None
  .type = str
  .help = pandas file list 
spectra_list = None
  .type = str
  .help = .lam spectra file (precognition) list  
pandas_table = None
  .type = str
  .help = path to a pandas table containing at least one
  .help = entry (data on one shot) output by stage_one
  .help = If present, this will override the other input parameters
oversample_override = None
  .type = int
  .help = force the nanoBragg oversample property to be this value (pixel oversample rate)
Ncells_abc_override = None
  .type = ints(size=3)
  .help = force the nanoBragg Ncells_abc property to be this 3-tuple (Na,Nb,Nc mosaic block size in unit cells)
cuda = False
  .type = bool
  .help = whether to use cuda
d_max = 999
  .type = float
  .help = maximum resolution
d_min = 1.4
  .type = float
  .help = minimum resolution
output_img = None
  .type = str
  .help = a name specifying an output image to write the model
outfile = None
  .type = str
  .help = output reflection file for indexed refls
tolerance = 1
  .type = float
  .help = indexing toleraance for assigning indices to the modeled spots
thresh = 1
  .type = float
  .help = threshold in photons for a modeled pixel to be flagged as part of a Bragg spot
Qdist_cutoff = 0.003
  .type = float
  .help = minimum distance in recip space for a strong spot to a modeled spot
  .help = in order that the strong spot be indexed
ngpu = 1
  .type = int
  .help = number of GPUs to use
njobs = 1
  .type = int
  .help = number of jbs to use, each job will use a randomly assigned gpu, up to ngpu
max_process = None
  .type = int
  .help = maximum number of imgs to predict
pandas_outfile = None
  .type = str
  .help = output file name (this file is suitable input to stage_two refinement, e.g. the pandas_table parameter)
"""

phil_scope = parse(script_phil)

class Script:

    def __init__(self):
        from dials.util.options import OptionParser

        self.parser = OptionParser(
            usage="",
            sort_options=True,
            phil=phil_scope,
            read_experiments=False,
            read_reflections=False,
            check_format=False,
            epilog=help_message)

    def run(self):
        self.params, _ = self.parser.parse_args(show_diff_phil=True)

        if self.params.pandas_table is None:
            explist = self.load_filelist(self.params.exper_list)
            reflist = [None] * len(explist)
            if self.params.refl_list is not None:
                reflist = self.load_filelist(self.params.refl_list)
            speclist = [None]*len(explist)
            if self.params.spectra_list is not None:
                speclist = self.load_filelist(self.params.spectra_list)
            pandalist = self.load_filelist(self.params.panda_list)
            assert (len(reflist) == len(speclist) == len(explist))
            input_iterator = zip(explist, reflist, speclist, pandalist)
            using_pandas_table = False
            NUM_EXPER = len(pandalist)
        else:
            df = pandas.read_pickle(self.params.pandas_table)
            NUM_EXPER = len(df)
            if "opt_exp_name" not in list(df):
                raise KeyError("Pandas dataframe needs to include a path to an optimized experiment")
            def _input_iterator(df):
                exper_names = df.opt_exp_name.unique()
                for i_exp, exper_name in enumerate(exper_names):
                    print("Pandas frame for exper %s" % exper_name)
                    exper_dataframe = df.query('opt_exp_name=="%s"'% exper_name)
                    yield exper_name, None, None, exper_dataframe
            input_iterator = _input_iterator(df)
            using_pandas_table = True

        processed_frames = []
        for i_exp, (exper_file, strong_file, spec_file, panda_file) in enumerate(input_iterator):
            # NOTE panda_file can be a dataframe or a string
            if i_exp % COMM.size != COMM.rank:
                continue

            if self.params.max_process is not None and i_exp >= self.params.max_process:
                break

            print("<><><><><><><><><><><><><><><><><><>")
            print("\tRank %d : iter %d / %d" % (COMM.rank, i_exp+1, NUM_EXPER))
            print("<><><><><><><><><><><><><><><><><><>")

            El = ExperimentListFactory.from_json_file(exper_file, check_format=False)
            exper = El[0]

            if strong_file is not None:
                strong = flex.reflection_table.from_file(strong_file)
            else:
                strong = None

            tstart = time.time()
            dev_id = COMM.rank % self.params.ngpu
            if not using_pandas_table:
                model_imgs = utils.spots_from_pandas_and_experiment(exper, panda_file,
                    spectrum_file=spec_file,
                    cuda=self.params.cuda, d_max=self.params.d_max, d_min=self.params.d_min,
                    output_img=self.params.output_img,
                    njobs=self.params.njobs, device_Id=dev_id, as_numpy_array=True)
            else:
                model_imgs = utils.spots_from_pandas(panda_file,
                    oversample_override=self.params.oversample_override,
                    Ncells_abc_override=self.params.Ncells_abc_override,
                    cuda=self.params.cuda, d_max=self.params.d_max, d_min=self.params.d_min,
                    output_img=self.params.output_img,
                    njobs=self.params.njobs, device_Id=dev_id)

            Rindexed = utils.indexed_from_model(strong, model_imgs, exper, tolerance=self.params.tolerance,
                                                thresh=self.params.thresh, Qdist_cutoff=self.params.Qdist_cutoff)
            Rindexed['id'] = flex.int(len(Rindexed), 0)
            if strong is not None:
                Rindexed = utils.remove_multiple_indexed(Rindexed)
                print("%d / %d are indexed!" % (len(Rindexed), len(strong)))
            prediction_outfile = os.path.splitext(exper_file)[0] + "_diffBragg_prediction.refl"
            Rindexed.as_file(prediction_outfile)
            tdone = time.time() - tstart
            print("Done, saved indexed refls to file %s (took %.4f sec)" % (prediction_outfile, tdone))

            if using_pandas_table:
                panda_file["predictions"] = os.path.abspath(prediction_outfile)
                processed_frames.append(panda_file)

        if using_pandas_table:
            processed_frames = COMM.reduce(processed_frames)
            if COMM.rank == 0:
                if self.params.pandas_outfile is None:
                    master_outfile = os.path.splitext(self.params.pandas_table)[0] + "_predictions.pkl"
                else:
                    master_outfile = self.params.pandas_outfile
                master_frame = pandas.concat(processed_frames)
                master_frame.to_pickle(master_outfile)
                print("Saved predictions dataframe : %s" % master_outfile)


    @staticmethod
    def load_filelist(fname):
        lines = [l.strip() for l in open(fname, "r").readlines()]
        for l in lines:
            if len(l.split()) > 1:
                raise RuntimeError("Input file %s is weird, needs single file path per line" % fname)
            if not os.path.exists(l):
                raise RuntimeError("fpath %s does not exist (from file list in %s)" % (l, fname))
        return lines


if __name__ == '__main__':
    from dials.util import show_mail_on_error
    with show_mail_on_error():
        script = Script()
        script.run()
