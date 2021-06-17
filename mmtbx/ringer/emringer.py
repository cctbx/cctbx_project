from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx.map_model_manager import map_model_manager
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
ignore_symmetry_conflicts = False
  .type = bool
  .help = Allows using PDB file with symmetry that does not match map
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
    if (self.map_inp is not None and self.model is not None):
      self.base = map_model_manager(
        map_manager      = self.map_inp,
        model            = self.model,
        ignore_symmetry_conflicts = self.params.ignore_symmetry_conflicts)

      self.cs_consensus = self.base.crystal_symmetry()
    else:
      self.base = None

  def run(self):
    hierarchy = self.model.get_hierarchy()
    map_data, grid_unit_cell = None, None
    if self.base is not None:
      hierarchy = self.base.model().get_hierarchy()
      map_data = self.base.map_manager().map_data()
      grid_unit_cell = self.base.map_manager().grid_unit_cell()

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
