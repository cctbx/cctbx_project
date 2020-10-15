from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.prediction

from dials.array_family import flex
from simtbx.diffBragg import utils
import time
from libtbx.phil import parse


help_message = "predictions using diffBragg refinement results"

script_phil = """
panda_name = None
  .type = str
  .help = pandas file output by diffBRagg
spectrum_file = None
  .type = str
  .help = path to input spectrum file in precognition format (.lam)
  .help = e.g. two columns: wavelength weight
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
"""

phil_scope = parse(script_phil)

class Script:

  def __init__(self):
    from dials.util.options import OptionParser

    self.parser = OptionParser(
      usage="",  # stage 1 (per-shot) diffBragg refinement",
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    self.params, _ = self.parser.parse_args(show_diff_phil=True)
    from dials.util.options import flatten_experiments, flatten_reflections
    explist = flatten_experiments(self.params.input.experiments)
    reflist = flatten_reflections(self.params.input.reflections)

    for El, strong in zip(explist, reflist):
      tstart = time.time()
      model_imgs = utils.spots_from_pandas_and_experiment(El, self.params.panda_name,
        spectrum_file=self.params.spectrum_file,
        cuda=self.params.cuda, d_max=self.params.d_max, d_min=self.params.d_min,
        output_img=self.params.output_img,
        njobs=self.params.njobs, ngpu=self.params.ngpu, as_numpy_array=True)

      Rindexed = utils.indexed_from_model(strong, model_imgs, El[0], tolerance=self.params.tolerance,
                                   thresh=self.params.thresh, Qdist_cutoff=self.paras.Qdist_cutoff)
      Rindexed['id'] = flex.int(len(Rindexed), 0)
      print("%d / %d are indexed!" % (len(Rindexed), len(strong)))
      Rindexed.as_file(self.params.outfile)
      tdone = time.time() - tstart
      print("Done, saved indexed refls to file %s (took %.4f sec)" % (self.params.outfile, tdone))


if __name__ == '__main__':
  from dials.util import show_mail_on_error
  with show_mail_on_error():
    script = Script()
    script.run()
