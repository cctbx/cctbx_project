"""Tool to calculate Rama-Z score. Validation of Ramachandran plot"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import rama_z
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.ramalyze import res_type_labels
from cctbx.maptbx.box import shift_and_box_model

from libtbx.utils import Sorry, null_out
from libtbx import Auto
import os

# =============================================================================

class Program(ProgramTemplate):

  description = '''
mmtbx.rama_z: Tool to calculate Rama-Z score. Validation of Ramachandran plot.

Usage examples:
  mmtbx.rama_z model1.pdb
  '''

  datatypes = ['model', 'phil']
  data_manager_options = ['model_skip_expand_with_mtrix']

  master_phil_str = """\
  write_HSL_models = False
    .type = bool
  write_HSL_plot = False
    .type = bool
  write_HSL_general_only = True
    .type = bool
  write_whole_plot = False
    .type = bool
  write_whole_general_only = True
    .type = bool
"""
  # write everything:
  # write_HSL_models=True write_HSL_plot=True write_HSL_general_only=False write_whole_plot=True write_whole_general_only=False
  # write only general plots:
  # write_HSL_plot=True write_whole_plot=False
  #
  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models()
    m = self.data_manager.get_model()
    print ('Inputs OK', file=self.logger)

  # ---------------------------------------------------------------------------

  def _write_plots_if_needed(self, model, label, type_of_plot='whole'):
    write_plot = getattr(self.params, "write_%s_plot" % type_of_plot)
    write_general_only = getattr(self.params, "write_%s_general_only" % type_of_plot)
    if write_plot:
      self.rama = ramalyze(model.get_hierarchy(), out=null_out())
      self.plots = self.rama.get_plots(
          show_labels=False,
          point_style='.',
          markersize=3,
          markeredgecolor="red",
          dpi=300,
          markerfacecolor="yellow")
      plots_to_write = range(6)
      if write_general_only:
        plots_to_write = [0]
      for i in plots_to_write:
        file_label = res_type_labels[i].replace("/", "_")
        fn = "%s.png" % self.get_default_output_filename(
            prefix='%s_%s_' % (self.inp_fn, label),
            suffix=file_label,
            serial=Auto)
        if os.path.isfile(fn) and not self.params.output.overwrite:
          raise Sorry("%s already exists and overwrite is set to False." % fn)
        print("Saving:", fn, file=self.logger)
        self.plots[i].save_image(fn, dpi=300)

  def run(self):
    self.result = None
    models = []
    for model_name in self.data_manager.get_model_names():
      models.append(self.data_manager.get_model(model_name))

    # model = self.data_manager.get_model()
    self.inp_fn = os.path.basename(self.data_manager.get_default_model_name())[:-4]
    self.rama_z = rama_z.rama_z(
        models = models,
        log = self.logger)
    if len(models) == 1:
      model = models[0]
      cs = model.crystal_symmetry()
      if (cs is None) or (cs.unit_cell() is None):
        model = shift_and_box_model(model)
      self._write_plots_if_needed(model, label='whole', type_of_plot='whole')
      helix_sel, sheet_sel, loop_sel = self.rama_z.get_ss_selections()
      if model.get_hierarchy().models_size() != 1:
        print ("Warning! Outputting partial models and plots are not supported \
  for multi-model files", file=self.logger)
      else:
        for sel, label in [(helix_sel, "helix"),
             (sheet_sel, "sheet"),
             (loop_sel, "loop")]:
          selected_model = model.select(sel)
          if self.params.write_HSL_models:
            fn = "%s" % self.get_default_output_filename(
                prefix='%s_' % self.inp_fn,
                suffix=label,
                serial=Auto)
            print("Writing out partial model: %s" % fn, file=self.logger)
            written_fn = self.data_manager.write_model_file(selected_model, filename=fn)
          self._write_plots_if_needed(selected_model, label, type_of_plot='HSL')
    self.result = self.get_results()
    # This brings 0 value to the user. Per-residue numbers
    # should not be analyzed, therefore no reason to print them.
    # res_info = self.rama_z.get_detailed_values()
    # print ("Individual residues info:", file=self.logger)
    # print ("Residue name, type, SS, (phi, psi), Z", file=self.logger)
    # for i in res_info:
    #   print ('%4s %10s %1s (%7.2f, %7.2f) %7.4f' % (
    #       i[2], res_type_labels[i[1]], i[3], i[4], i[5], i[6]), file=self.logger)

    print(self.result.as_string(prefix=""), file = self.logger)
    # print(self.result.as_json(), file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    if self.result is None:
      self.result = self.rama_z.get_result()
    return self.result

  def get_results_as_JSON(self):
    if self.result is None:
      self.result = self.rama_z.get_result()
    return self.result.as_json()
