# -*- coding: utf-8 -*-
from __future__ import division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import compare_rama
from mmtbx.validation.ramalyze import res_type_labels

import numpy as np
from collections import Counter

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.compare_rama: tool for compare Ramachandran plots, e.g. before-after
  refinement.

Usage examples:
  phenix.compare_rama model1.pdb model2.pdb
  phenix.compare_rama model1.cif model2.cif
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
    include scope mmtbx.validation.compare_rama.master_phil_str
    output
    {
      individual_residues = True
        .type = bool
      sorted_individual_residues = False
        .type = bool
      counts = True
        .type = bool
      prefix = kleywegt
        .type = str
      plots = True
        .type = bool
        .help = output Kleywegt plots - arrows on Rama plot showing where \
          residues moved.
    }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=2, exact_count=True, raise_sorry=True)
    model_1, model_2 = self._get_models()
    assert model_1.get_hierarchy().is_similar_hierarchy(model_2.get_hierarchy())
    for m in [model_1, model_2]:
      assert m.get_hierarchy().models_size() == 1

  # ---------------------------------------------------------------------------
  def run(self):
    # I'm guessing self.data_manager, self.params and self.logger
    # are already defined here...
    # print('Using model: %s' % self.data_manager.get_default_model_name(), file=self.logger)

    # this must be mmtbx.model.manager?
    model_1, model_2 = self._get_models()

    self.rama_comp = compare_rama.rcompare(
        model1 = model_1,
        model2 = model_2,
        params = self.params.compare_rama,
        log = self.logger)

    # outputting results
    results = self.rama_comp.get_results()
    res_columns = zip(*results)
    if self.params.output.individual_residues:
      for r in results:
        self.show_single_result(r)
      print("="*80, file=self.logger)
    if self.params.output.sorted_individual_residues:
      sorted_res = sorted(results, key=lambda tup: tup[1])
      for r in sorted_res:
        self.show_single_result(r)
      print("="*80, file=self.logger)
    print ("mean: %.3f std: %.3f" % (np.mean(res_columns[1]), np.std(res_columns[1])),
        file=self.logger)
    if self.params.output.counts:
      cntr = Counter(res_columns[-2])
      for k, v in cntr.iteritems():
        print("%-20s: %d" % (k,v), file=self.logger)

    if self.params.output.plots:
      rama1, rama2 = self.rama_comp.get_ramalyze_objects()
      plots = rama2.get_plots(
          show_labels=True,
          point_style='bo',
          markersize=3,
          markeredgecolor="black",
          dpi=300,
          markerfacecolor="white")
      for pos, plot in plots.iteritems():
        # prepare data
        got_outliers = [x for x in results if (x[-1]==pos and x[-2].find("-> OUTLIER") > 0)]#.sort(key=lambda x:x[1], reverse=True)
        got_outliers.sort(key=lambda x:x[1], reverse=True)
        print("got_outliers:", len(got_outliers))
        for o in got_outliers:
          self.show_single_result(o)
        got_not_outliers = [x for x in results if (x[-1]==pos and x[-2] == "OUTLIER -> Favored")]#.sort(key=lambda x:x[1], reverse=True)
        got_not_outliers.sort(key=lambda x:x[1], reverse=True)
        print("got_not_outliers:", len(got_not_outliers))
        for o in got_not_outliers:
          self.show_single_result(o)

        for data, color in [(got_outliers, "red"), (got_not_outliers, "lime")]:
          # print (len(data))
          if data and len(data) < 0: continue
          ad = [((x[2], x[3]),(x[4], x[5])) for x in data]
          add_arrows_on_plot(
              plot,
              ad,
              color=color)
        file_label = res_type_labels[pos].replace("/", "_")
        base_fname = "%s--%s" % (self.data_manager.get_model_names()[0].split('.')[0],
            self.data_manager.get_model_names()[1].split('.')[0])
        plot_file_name = "%s_%s_%s_plot.png" % (
            base_fname, self.params.output.prefix, file_label)
        print("saving: '%s'" % plot_file_name)
        plot.save_image(plot_file_name, dpi=300)

  def show_single_result(self, r):
    print("%s %.2f, (%.1f:%.1f), (%.1f:%.1f), %s" % r[:-1], file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.rama_comp.get_results()

  def _get_models(self):
    m_names = self.data_manager.get_model_names()
    model_1 = self.data_manager.get_model(filename=m_names[0])
    model_2 = self.data_manager.get_model(filename=m_names[1])
    return model_1, model_2

def add_arrows_on_plot(p, arrows_data, color='green'):
  """
  p - pyplot
  arrows_data - [((x,y beginning), (x,y end)), ... ((xy),(xy))]
  """
  import matplotlib.patches as patches
  import matplotlib
  style="Simple,head_length=10,head_width=5,tail_width=1"
  for arrow in arrows_data:

    p.plot.add_patch(patches.FancyArrowPatch(
        arrow[0],
        arrow[1],
        arrowstyle=style,
        color = color,
        linewidth=0.5,
        ))
