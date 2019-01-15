# LIBTBX_SET_DISPATCHER_NAME phenix.ramalyze
# LIBTBX_SET_DISPATCHER_NAME molprobity.ramalyze
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import division
import mmtbx.validation.ramalyze
import iotbx.phil
from libtbx.utils import Usage, Sorry
import os.path
import os, sys

def get_master_phil():
  return iotbx.phil.parse(input_string="""
    include scope mmtbx.validation.molprobity_cmdline_phil_str
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
    model_list = None
      .type = str
      .help = Comma separated file list to accumulate onto one plot
    output_prefix = None
      .type = str
      .help = prefix for outputted plots (if plot=True)
""", process_includes=True)
prog = os.getenv('LIBTBX_DISPATCHER_NAME')
usage_string = """
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file
  outliers_only=False   only print outliers
  verbose=False         verbose text output
  plot=False            Create graphics of plots (if Matplotlib is installed)

Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

def run(args, out=sys.stdout, quiet=False):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None and params.model_list is None):
    raise Usage(usage_string)
  if params.model:
    models = [params.model]
  elif params.model_list:
    if os.path.isfile(params.model_list):
      with open(params.model_list, 'r') as f:
        models = f.read().split('\n')
        models = [m for m in models if m != ""]
    else:
      models = params.model_list.split(',')
    params.verbose=False
    params.model=models[0]
  results = []
  for model in models:
    if not os.path.isfile(model) and params.model_list:
      print "Cannot find '%s', skipping." % model
      continue
    pdb_in = cmdline.get_file(model, force_type="pdb")
    hierarchy = pdb_in.file_object.hierarchy
    hierarchy.atoms().reset_i_seq()
    result = mmtbx.validation.ramalyze.ramalyze(
      pdb_hierarchy=hierarchy,
      show_errors=None,
      outliers_only=params.outliers_only,
      out=out,
      quiet=quiet)
    results.append(result)
    if params.model_list:
      print '\nmodel  : %s' % model
      result.show_summary()
  # combine
  result = results[0]
  for i in range(1,len(results)):
    result += results[i]
  if params.verbose:
    result.show_old_output(out=out, verbose=True)
  if params.plot :
    plot_file_base = params.output_prefix
    if plot_file_base is None:
      plot_file_base = os.path.splitext(os.path.basename(params.model))[0]
    result.write_plots(
        plot_file_base=plot_file_base,
        out=out,
        show_labels=params.show_labels,
        point_style=params.point_style,
        markerfacecolor=params.markerfacecolor,
        show_filling=params.show_filling,
        show_contours=params.show_contours,
        dpi=params.dpi,
        markeredgecolor=params.markeredgecolor,
        markersize=params.markersize)
  if params.wxplot :
    try :
      import wxtbx.app
    except ImportError, e :
      raise Sorry("wxPython not available.")
    else :
      app = wxtbx.app.CCTBXApp(0)
      result.display_wx_plots()
      app.MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
