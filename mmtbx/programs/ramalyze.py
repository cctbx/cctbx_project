"""Report on Ramachandran values for a model"""
from __future__ import absolute_import, division, print_function

import os
import iotbx.phil
from mmtbx.validation.ramalyze import ramalyze
from libtbx.program_template import ProgramTemplate
from six.moves import range
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
from libtbx.utils import Sorry
from datetime import datetime

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
  .type=float
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
json = False
  .type = bool
  .help = "Prints results as JSON format dictionary"
verbose = True
  .type = bool
  .help = '''Verbose'''
output_prefix = None
  .type = str
  .help = prefix for outputted plots (if plot=True)
"""

def master_params():
  return iotbx.phil.parse(master_phil_str)

def compute(hierarchies, params, log, quiet=False, plot_file_base_default=None):
  results = []
  for hierarchy in hierarchies:
    result = ramalyze(
      pdb_hierarchy = hierarchy,
      show_errors   = None,
      outliers_only = params.outliers_only,
      out           = log,
      quiet         = quiet)
    results.append(result)
  # combine models
  result = results[0]
  for i in range(1,len(results)):
    result += results[i]
  if params.verbose and not params.json:
    result.show_old_output(out=log, verbose=True)
  if params.plot:
    plot_file_base = params.output_prefix
    if plot_file_base is None:
      plot_file_base = plot_file_base_default
    result.write_plots(
      plot_file_base  = plot_file_base,
      out             = log,
      show_labels     = params.show_labels,
      point_style     = params.point_style,
      markerfacecolor = params.markerfacecolor,
      show_filling    = params.show_filling,
      show_contours   = params.show_contours,
      dpi             = params.dpi,
      markeredgecolor = params.markeredgecolor,
      markersize      = params.markersize)
  if params.wxplot:
    try:
      import wxtbx.app
    except ImportError as e:
      raise Sorry("wxPython not available.")
    else:
      app = wxtbx.app.CCTBXApp(0)
      result.display_wx_plots()
      app.MainLoop()
  else:
    return result

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file
  outliers_only=False   only print outliers
  json=False            Outputs results as JSON compatible dictionary
  verbose=False         verbose text output
  plot=False            Create graphics of plots (if Matplotlib is installed)

Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

  # Pavel's style:
  # plot=True show_labels=False markerfacecolor=yellow markeredgecolor=red

  master_phil_str = master_phil_str

  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def get_results_as_JSON(self):
    # this calculates the results separately from run() because historically
    # the ramalyze object couldn't handle multi-model files. Multi-model support was
    # added for the JSON code.  Ideally this would get fixed in the future.
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    result = ramalyze(
      pdb_hierarchy = hierarchy,
      outliers_only = self.params.outliers_only,
      out           = self.logger)
    return result.as_JSON(self.info_json)

  def run(self):
    if self.params.json:
      print(self.get_results_as_JSON())
    hierarchies = []
    for model_name in self.data_manager.get_model_names():
      hierarchy = self.data_manager.get_model(model_name).get_hierarchy()
      hierarchies.append(hierarchy)
    fb = os.path.splitext(os.path.basename(
      self.data_manager.get_model_names()[0]))[0]
    self.results = compute(
      hierarchies            = hierarchies,
      params                 = self.params,
      log                    = self.logger,
      quiet                  = False,
      plot_file_base_default = fb)

  def get_results(self):
    return self.results
