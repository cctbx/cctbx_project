"""Tool for compare Ramachandran plots, e.g. before-after
  refinement"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from scitbx.array_family import flex

from mmtbx.validation import comparama
from mmtbx.validation.ramalyze import res_type_labels

from matplotlib.backends.backend_pdf import PdfPages

import os
import six

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.comparama: tool for compare Ramachandran plots, e.g. before-after
  refinement.

Usage examples:
  phenix.comparama model1.pdb model2.pdb
  phenix.comparama model1.cif model2.cif
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
    include scope mmtbx.validation.comparama.master_phil_str
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
      plots = False
        .type = bool
        .help = output Kleywegt plots - arrows on Rama plot showing where \
          residues moved.
      pdf = True
        .type = bool
        .help = save the same plots as one pdf file
      wrap_arrows = True
        .type = bool
        .help = wrap long arrows around the edges of plot
    }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=2, exact_count=True, raise_sorry=True)
    model_1, model_2 = self._get_models()
    model_1.add_crystal_symmetry_if_necessary()
    model_2.add_crystal_symmetry_if_necessary()
    # Filter out what you are not interested in
    s1 = model_1.selection(string="protein")
    s2 = model_2.selection(string="protein")
    model_1 = model_1.select(selection = s1)
    model_2 = model_2.select(selection = s2)
    #
    # assert model_1.get_hierarchy().is_similar_hierarchy(model_2.get_hierarchy())
    for m in [model_1, model_2]:
      assert m.get_hierarchy().models_size() == 1

  # ---------------------------------------------------------------------------
  def run(self):
    # I'm guessing self.data_manager, self.params and self.logger
    # are already defined here...
    # print('Using model: %s' % self.data_manager.get_default_model_name(), file=self.logger)

    # this must be mmtbx.model.manager?
    model_1, model_2 = self._get_models()
    # Filter out what you are not interested in
    s1 = model_1.selection(string="protein")
    s2 = model_2.selection(string="protein")
    model_1 = model_1.select(selection = s1)
    model_2 = model_2.select(selection = s2)

    self.rama_comp = comparama.rcompare(
        model1 = model_1,
        model2 = model_2,
        params = self.params.comparama,
        log = self.logger)

    # outputting results
    results = self.rama_comp.get_results()
    if len(results) == 0:
      print("No ramachandran residues found!")
      return
    if self.params.output.individual_residues:
      for r in results:
        self.show_single_result(r)
      print("="*80, file=self.logger)
    if self.params.output.sorted_individual_residues:
      sorted_res = sorted(results, key=lambda tup: tup[1])
      for r in sorted_res:
        self.show_single_result(r)
      print("="*80, file=self.logger)

    skip1, skip2 = self.rama_comp.get_skipped()
    for s in skip1:
      print("WARNING! No match for '%s' in the second model" % s.id_str(), file=self.logger)
    for s in skip2:
      print("WARNING! No match for '%s' in the first model" % s.id_str(), file=self.logger)
    if len(skip1) + len(skip2) > 0:
      print("="*80, file=self.logger)

    nr = self.rama_comp.get_number_results()
    print ("mean: %.3f std: %.3f" % (nr.mean_diff, nr.std_diff), file=self.logger)
    print("Sum of rama scores: \t\t\t %.3f -> %.3f" % (nr.sum_1, nr.sum_2) , file=self.logger)
    print("Sum of rama scores/n_residues:\t\t %.4f -> %.4f (%d residues)" % \
        (nr.sum_1/nr.n_res, nr.sum_2/nr.n_res, nr.n_res), file=self.logger)
    print("Sum of rama scores scaled:\t\t %.3f -> %.3f" % \
        (nr.scaled_sum_1, nr.scaled_sum_2) , file=self.logger)
    print("Sum of rama scores/n_residues scaled:\t %.4f -> %.4f (%d residues)" % \
        (nr.scaled_sum_1/nr.n_res, nr.scaled_sum_2/nr.n_res, nr.n_res), file=self.logger)
    print("Sum of rama scores reverse scaled:\t\t %.3f -> %.3f" % \
        (nr.rev_scaled_sum_1, nr.rev_scaled_sum_2) , file=self.logger)
    print("Sum of rama scores/n_residues reverse scaled:\t %.4f -> %.4f (%d residues)" % \
        (nr.rev_scaled_sum_1/nr.n_res, nr.rev_scaled_sum_2/nr.n_res, nr.n_res), file=self.logger)

    if self.params.output.counts:
      for k, v in six.iteritems(nr.counts):
        print("%-20s: %d" % (k,v), file=self.logger)

    # Pavel's numbers
    s1, s2 = self.rama_comp.get_results_as_vec3()
    pnumber = flex.mean(flex.sqrt((s1-s2).dot()))
    # compare these with log output:
    #  A   2  ASN 25.13, (-60.6:141.2), (-78.7:158.6), Favored, Score: 0.5693 -> 0.2831
    #                      phi1  psi1     phi2  psi2
    # print (list(s1))
    # print (list(s2))
    print ("Pavel's test number: %.4f" % pnumber)

    name1 = os.path.basename(self.data_manager.get_model_names()[0]).split('.')[0]
    name2 = os.path.basename(self.data_manager.get_model_names()[1]).split('.')[0]
    base_fname = "%s--%s" % (name1, name2)
    if self.params.output.plots:
      for pos, plot in six.iteritems(self.rama_comp.get_plots(wrap_arrows=self.params.output.wrap_arrows)):
        file_label = res_type_labels[pos].replace("/", "_")
        plot_file_name = "%s_%s_%s_plot.png" % (
            base_fname, self.params.output.prefix, file_label)
        print("saving: '%s'" % plot_file_name)
        plot.save_image(plot_file_name, dpi=300)

    if self.params.output.pdf:
      pdf_fname = "%s_%s.pdf" % (base_fname, self.params.output.prefix)
      pdfp = PdfPages(pdf_fname)
      for pos, plot in six.iteritems(self.rama_comp.get_plots(wrap_arrows=self.params.output.wrap_arrows)):
        pdfp.savefig(plot.figure)
      print("saving: '%s'" % pdf_fname)
      pdfp.close()

  def show_single_result(self, r):
    print("%s %.2f, (%.1f:%.1f), (%.1f:%.1f), %s, Score: %.4f -> %.4f" % \
        (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[8], r[9]),
        file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.rama_comp.get_results()

  def _get_models(self):
    m_names = self.data_manager.get_model_names()
    model_1 = self.data_manager.get_model(filename=m_names[0])
    model_2 = self.data_manager.get_model(filename=m_names[1])
    return model_1, model_2
