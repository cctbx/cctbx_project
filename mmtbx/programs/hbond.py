"""Tool to find H bonds in an atomic model"""
from __future__ import division, print_function
from libtbx.program_template import ProgramTemplate
import mmtbx.nci.hbond
import mmtbx.nci.skew_kurt_plot
from libtbx.utils import null_out
import os

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.hbond: tool to find H bonds in an atomic model

Usage example:
  phenix.hbond model.pdb
  '''

  datatypes = ['model', 'phil', 'restraint']


  master_phil_str = '''
  include scope mmtbx.nci.hbond.master_phil_str

  output_skew_kurtosis_plot = False
    .type = bool
    .short_caption = Output skew-kurtosis plot with result in png format
  plot_colorblind_friendly = True
    .type = bool
    .short_caption = Use colorblind friendly palette for skew-kurtosis plot
  do_rotate_translate = False
    .type = bool
    .style = hidden
  do_y_log = False
    .type = bool
    .style = hidden
  output_sk_coordinates = False
    .type = bool
  plot_parameters_override
    .help = These parameters will override preset values in plots. The values \
      will be passed directly to matplotlib functions, so should be valid \
      matplotlib color names.
  {
    colormap = None
      .type = str
    contour_left = None
      .type = str
    contour_right = None
      .type = str
    contour_thick_left = None
      .type = float(value_min=0.1, value_max=10)
    contour_thick_right = None
      .type = float(value_min=0.1, value_max=10)
    theta_color = None
      .type = str
    theta_contour = None
      .type = str
    Rha_color = None
      .type = str
    Rha_contour = None
      .type = str
  }
  '''
  # % mmtbx.nci.hbond.master_phil_str

  # ---------------------------------------------------------------------------
  def validate(self):
    self._print('Validating inputs')
    self.data_manager.has_models(raise_sorry=True)

  # ---------------------------------------------------------------------------
  def run(self):
    print('Using model: %s' % self.data_manager.get_default_model_name(),
      file=self.logger)
    inp_models = []
    for model_name in self.data_manager.get_model_names():
      inp_models.append((model_name, self.data_manager.get_model(model_name)))
    theta1_data = []
    Rha_data = []
    if self.params.output_sk_coordinates:
      file_names = []
    for m_name, model in inp_models:
      model.set_log(log = null_out())
      if self.params.hbond.add_hydrogens_if_absent and not model.has_hd():
        from elbow.command_line.ready_set import model_interface as ready_set_model_interface
        model = ready_set_model_interface(
            model=model,
            params=["add_h_to_water=False",
                    "optimise_final_geometry_of_hydrogens=False"],
            )
      model.process(make_restraints=True)
      self.results = mmtbx.nci.hbond.find(model = model)
      if self.params.hbond.show_hbonds:
        self.results.show(log = self.logger)
      self._print("-"*79)
      #self.results.show_summary(log = self.logger)
      prefix=self.params.output.prefix
      if not prefix:
        prefix='%s_hbond' % os.path.basename(m_name).split('.')[0]
      if self.params.hbond.output_pymol_file:
        self.results.as_pymol(prefix=prefix)
      if self.params.hbond.output_restraint_file:
        self.results.as_restraints(file_name='%s.eff' % prefix)
      #
      stats = mmtbx.nci.hbond.stats(model = model, prefix="%s_stats" % prefix,
        output_stats_pdf = self.params.hbond.output_stats_pdf)
      if stats is None:
        self._print('\n\tLimited number of H-bonds so statistics are not calculated.')
        return
      min_data_size=self.params.hbond.min_data_size
      # These are the values to be used for the plot!
      stats.all.show_summary(log = self.logger)
      if stats.all.get_counts(min_data_size=min_data_size):
        theta1_data.append(stats.all.get_counts(min_data_size=min_data_size).theta_1)
        Rha_data.append(   stats.all.get_counts(min_data_size=min_data_size).d_HA)
        if self.params.output_sk_coordinates:
          file_names.append(m_name)

    if self.params.output_skew_kurtosis_plot and len(theta1_data) > 0 and len(Rha_data) > 0:
      # To use other than 'all' type, nci.hbond.find needs to be called with selected model again,
      # like in stats().
      fn = '%s_skew_kurtosis' % prefix
      if self.params.plot_colorblind_friendly:
        fn += "_cbf"
      theta1_c = [(x.skew, x.kurtosis) for x in theta1_data]
      Rha_c = [(x.skew, x.kurtosis) for x in Rha_data]

      op = {}
      for param_name in dir(self.params.plot_parameters_override):
        if not param_name.startswith('__'):
          param_value = getattr(self.params.plot_parameters_override, param_name)
          if param_value != None:
            op[param_name] = param_value

      if self.params.output_sk_coordinates:
        from libtbx import group_args
        from libtbx import easy_pickle
        o = group_args(
          file_names    = file_names,
          theta1_coords = theta1_c,
          Rha_coords    = Rha_c)
        easy_pickle.dump(file_name="%s_coords.pkl"%prefix, obj = o)

      mmtbx.nci.skew_kurt_plot.make_figure(
          file_name=fn,
          theta1_coords=theta1_c,
          Rha_coords=Rha_c,
          dot_size = self.params.hbond.dot_size,
          type='all',
          override_palette = op,
          colorblind_friendly=self.params.plot_colorblind_friendly,
          do_rotate=self.params.do_rotate_translate,
          do_y_log=self.params.do_y_log,
          )
      self._print("\nOutputted plot as %s.png" % fn)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

# =============================================================================
# end
