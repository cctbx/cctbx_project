from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.prediction

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
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
"""

phil_scope = parse(script_phil)

class Script:

  def __init__(self):
    from dials.util.options import OptionParser

    self.parser = OptionParser(
      usage="",  # stage 1 (per-shot) diffBragg refinement",
      sort_options=True,
      phil=phil_scope,
      read_experiments=False,
      read_reflections=False,
      check_format=False,
      epilog=help_message)

  def run(self):
    self.params, _ = self.parser.parse_args(show_diff_phil=True)
    #from dials.util.options import flatten_experiments, flatten_reflections
    #explist = flatten_experiments(self.params.input.experiments)
    #reflist = flatten_reflections(self.params.input.reflections)

    explist = self.load_filelist(self.params.exper_list)
    reflist = [None] * len(explist)
    if self.params.refl_list is not None:
      reflist = self.load_filelist(self.params.refl_list)
    speclist = [None]*len(explist)
    if self.params.spectra_list is not None:
      speclist = self.load_filelist(self.params.spectra_list)
    pandalist = self.load_filelist(self.params.panda_list)
    assert (len(reflist) == len(speclist) == len(explist))

    for i_exp, (exper_file, strong_file, spec_file, panda_file) in enumerate(zip(explist, reflist, speclist, pandalist)):
      if i_exp % COMM.size != COMM.rank:
        continue

      if self.params.max_process is not None and i_exp >= self.params.max_process:
        break

      El = ExperimentListFactory.from_json_file(exper_file, check_format=False)
      exper = El[0]

      if strong_file is not None:
        strong = flex.reflection_table.from_file(strong_file)
      else:
        strong = None

      tstart = time.time()
      dev_id = COMM.rank % self.params.ngpu
      model_imgs = utils.spots_from_pandas_and_experiment(exper, panda_file,
        spectrum_file=spec_file,
        cuda=self.params.cuda, d_max=self.params.d_max, d_min=self.params.d_min,
        output_img=self.params.output_img,
        njobs=self.params.njobs, device_Id=dev_id, as_numpy_array=True)

      Rindexed = utils.indexed_from_model(strong, model_imgs, exper, tolerance=self.params.tolerance,
                                   thresh=self.params.thresh, Qdist_cutoff=self.params.Qdist_cutoff)
      Rindexed['id'] = flex.int(len(Rindexed), 0)
      if strong is not None:
        Rindexed = utils.remove_multiple_indexed(Rindexed)
        print("%d / %d are indexed!" % (len(Rindexed), len(strong)))
      #outfile = os.path.splitext(strong_file)[0] + "_diffBragg_idx.refl"
      outfile = os.path.splitext(exper_file)[0] + "_diffBragg_idx.refl"
      Rindexed.as_file(outfile)
      tdone = time.time() - tstart
      #print("Done, saved indexed refls to file %s (took %.4f sec)" % (self.params.outfile, tdone))
      print("Done, saved indexed refls to file %s (took %.4f sec)" % (outfile, tdone))

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
