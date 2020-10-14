from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.prediction

from dials.array_family import flex
from simtbx.diffBragg import utils
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.phil import parse


help_message = "predictions using diffBragg refinement results"

script_phil = """
El_name
  .type = str
  .help = refined experiment
panda_name
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
save_data_too = False
  .type = bool
  .help = whether to store experimental data in `output_img`
strong_refl
  .type = str
  .help = input experiment file
outfile
  .type = str
  .help = output reflection file for indexed refls
tolerance = 1
  .type = float
  .help = indexing toleraance for assigning indices to the `modeled` spots
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
  .help = number of jbs to use, each job will use a randomly assigned gpu, up to `ngpu`
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

    El = ExperimentListFactory.from_json_file(self.params.El_name, check_format=False)
    model_imgs = utils.spots_from_pandas_and_experiment(self.params.El_name, self.params.panda_name,
      spectrum_file=self.params.spectrum_file,
      cuda=self.params.cuda, d_max=self.params.d_max, d_min=self.params.d_min,
      output_img=self.params.output_img,
      njobs=self.params.njobs, ngpu=self.params.ngpu,
      save_expt_data=self.params.save_data_too, as_numpy_array=True)

    strong = flex.reflection_table.from_file(self.params.strong_refl)
    R = utils.indexed_from_model(strong, model_imgs, El[0], tolerance=self.params.tolerance,
                                 thresh=self.params.thresh, Qdist_cutoff=self.paras.Qdist_cutoff)
    R['id'] = flex.int(len(R), 0)
    print("%d / %d are indexed!" % (len(R), len(strong)))
    R.as_file(self.params.outfile)
    print("Done, saved indexed refls to file %s" % self.params.outfile)


if __name__ == '__main__':
  from dials.util import show_mail_on_error
  with show_mail_on_error():
    script = Script()
    script.run()


