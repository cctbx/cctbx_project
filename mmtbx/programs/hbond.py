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

  master_phil_str = mmtbx.nci.hbond.master_phil_str

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
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
      print("-"*79, file=self.logger)
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
      min_data_size=self.params.hbond.min_data_size
      # These are the values to be used for the plot!
      stats.all.show_summary(log = self.logger)
      theta1_data.append(stats.all.get_counts(min_data_size=min_data_size).theta_1)
      Rha_data.append(   stats.all.get_counts(min_data_size=min_data_size).d_HA)

    if self.params.hbond.output_skew_kurtosis_plot and self.results.get_counts():
      # To use other than 'all' type, nci.hbond.find needs to be called with selected model again,
      # like in stats().
      fn = '%s_skew_kurtosis' % prefix
      if self.params.hbond.plot_colorblind_friendly:
        fn += "_cbf"
      theta1_c = [(x.skew, x.kurtosis) for x in theta1_data]
      Rha_c = [(x.skew, x.kurtosis) for x in Rha_data]
      mmtbx.nci.skew_kurt_plot.make_figure(
          file_name=fn,
          theta1_coords=theta1_c,
          Rha_coords=Rha_c,
          dot_size = self.params.hbond.dot_size,
          type='all',
          colorblind_friendly=self.params.hbond.plot_colorblind_friendly)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

# =============================================================================
# end
