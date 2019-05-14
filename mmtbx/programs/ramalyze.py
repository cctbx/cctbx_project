from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.ramalyze import ramalyze
from libtbx.program_template import ProgramTemplate
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
from libtbx.utils import Sorry

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file
  outliers_only=False   only print outliers
  verbose=False         verbose text output
  plot=False            Create graphics of plots (if Matplotlib is installed)

Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

  # Pavel's style:
  # plot=True show_labels=False markerfacecolor=yellow markeredgecolor=red
  master_phil_str = """
    plot = False
      .type = bool
      .help = Create graphics of plots (if Matplotlib is installed)
    show_labels = True
      .type = bool
      .help = Show labels on outlier residues
    point_style = 'bo'
      .type = str
      .help = choose style of points, use matplotlib format from e.g. here: \
        https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.plot.html \
        very small is ',', little bigger is '.'
    markersize=3
      .type=int
    markerfacecolor = white
      .type = str
    markeredgecolor="black"
      .type = str
    show_filling = True
      .type = bool
    show_contours = True
      .type = bool
    dpi=100
      .type=int
    wxplot = False
      .type = bool
      .help = Display interactive plots (requires wxPython and Matplotlib)
    outliers_only = False
      .type = bool
      .help = "Only display outliers"
    verbose = True
      .type = bool
      .help = '''Verbose'''
    output_prefix = None
      .type = str
      .help = prefix for outputted plots (if plot=True)
"""
  datatypes = ['model','phil']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    results = []
    for model_name in self.data_manager.get_model_names():
      hierarchy = self.data_manager.get_model(model_name).get_hierarchy()
      hierarchy.atoms().reset_i_seq()
      result = ramalyze(
        pdb_hierarchy=hierarchy,
        show_errors=None,
        outliers_only=self.params.outliers_only,
        out=self.logger,
        quiet=False)
      results.append(result)
      if len(self.data_manager.get_model_names()) > 1:
        self.params.verbose=False
        print('\nmodel : %s' % model_name, file=self.logger)
    # combine models
    result = results[0]
    for i in range(1,len(results)):
      result += results[i]
    if self.params.verbose:
      result.show_old_output(out=self.logger, verbose=True)
    if self.params.plot :
      plot_file_base = self.params.output_prefix
      if plot_file_base is None:
        plot_file_base = os.path.splitext(os.path.basename(self.data_manager.get_model_names()[0]))[0]
      result.write_plots(
          plot_file_base=plot_file_base,
          out=self.logger,
          show_labels=self.params.show_labels,
          point_style=self.params.point_style,
          markerfacecolor=self.params.markerfacecolor,
          show_filling=self.params.show_filling,
          show_contours=self.params.show_contours,
          dpi=self.params.dpi,
          markeredgecolor=self.params.markeredgecolor,
          markersize=self.params.markersize)
    if self.params.wxplot :
      try :
        import wxtbx.app
      except ImportError as e :
        raise Sorry("wxPython not available.")
      else :
        app = wxtbx.app.CCTBXApp(0)
        result.display_wx_plots()
        app.MainLoop()

