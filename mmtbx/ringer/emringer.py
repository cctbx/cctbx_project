from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
import mmtbx.utils
from mmtbx.ringer import iterate_over_residues
from mmtbx.ringer import em_rolling
from mmtbx.ringer import em_scoring

master_params_str = '''
sampling_angle = 5
  .type = int
  .input_size = 64
sampling_method = linear *spline direct
  .help = Method for sampling
  .type = choice(multi=False)
grid_spacing = 1./5
  .type = float
scaling = *sigma volume
  .help = Method for map scaling.
  .type = choice(multi=False)
rolling_window_threshold = 0
  .type = float(value_min=0)
  .help = Threshold for calculating statistics across rolling windows of residues
skip_alt_confs = True
  .type = bool
nproc = 1
  .type = int
  .short_caption = Number of processors
  .input_size = 64
  '''

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class emringer(object):
  def __init__(self, model, miller_array, map_inp, params, out):
    self.model        = model
    self.miller_array = miller_array
    self.map_inp      = map_inp
    self.params       = params
    self.out          = out

  def validate(self):
    assert not None in [self.model, self.params, self.out]
    if (self.model is None):
      raise Sorry("Model is required.")
    if (self.miller_array is None and self.map_inp is None):
      raise Sorry("Map or map coefficients are required.")
    # Sanity check for crystal symmetry
    if (self.map_inp is not None):
      self.cs_consensus = mmtbx.utils.check_and_set_crystal_symmetry(
        models = [self.model], map_inps=[self.map_inp])

  def run(self):
    hierarchy = self.model.get_hierarchy()
    map_data, grid_unit_cell = None, None
    # sanity check for map and model
    if self.map_inp is not None:
      base = map_and_model.input(
        map_data         = self.map_inp.map_data(),
        model            = self.model,
        crystal_symmetry = self.cs_consensus,
        box              = False)

      hierarchy = base.model().get_hierarchy()
      map_data = base.map_data()
      grid_unit_cell = self.map_inp.grid_unit_cell()

    hierarchy.atoms().reset_i_seq()

    self.ringer_result = iterate_over_residues(
      pdb_hierarchy          = hierarchy,
      map_coeffs             = self.miller_array,
      map_data               = map_data,
      unit_cell              = grid_unit_cell,
      params                 = self.params,
      log                    = self.out
      ).results

    if (self.params.output_base is not None):
      plots_dir = self.params.output_base + "_plots"
    else:
      plots_dir = 'emringer_plots'

    import matplotlib
    matplotlib.use("Agg")
    self.scoring_result = em_scoring.main(
      file_name      = self.params.output_base,
      ringer_result  = self.ringer_result,
      out_dir        = plots_dir,
      sampling_angle = self.params.sampling_angle,
      quiet          = self.params.quiet,
      out            = self.out)

    rolling_window_threshold = self.params.rolling_window_threshold
    self.rolling_result = em_rolling.main(
      ringer_results = self.ringer_result,
      dir_name       = plots_dir,
      threshold      = rolling_window_threshold, #scoring.optimal_threshold,
      graph          = False,
      save           = not self.params.quiet,
      out            = self.out)

  def get_results(self):
    return group_args(
      ringer_result  = self.ringer_result,
      scoring_result = self.scoring_result,
      rolling_result = self.rolling_result)
