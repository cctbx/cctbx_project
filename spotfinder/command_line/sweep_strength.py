from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME distl.sweep_strength
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
from six.moves import range
import spotfinder
import libtbx.phil
from libtbx.utils import Sorry
import sys, os
from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs


master_params = libtbx.phil.parse("""\
distl {
  res {
    outer=None
      .type=float
      .help="High resolution limit in angstroms"
    inner=None
      .type=float
      .help="Low resolution limit in angstroms"
  }
  plot {
    file_name = None
      .type = path
      .help = "File path to output a plot of the distl results vs image number."
              "Supported formats include png, jpeg, pdf, ps, eps, svg."
  }
  csv = None
    .type = path
    .help = "Optional file name for output of distl results in CSV format."
  verbosity = 1
    .type = int
    .help="Control the amount of output:"
          "  0: print just the table"
          "  1: also print the parameters used by the program"
          "  2: also print distl summary for each individual image"
  bins {
    verbose = False
      .type = bool
      .help="Additional printout binned by resolution range"
    N = 20
      .type = int
      .help="Maximum number of bins, but fewer can result if there are few spots"
    corner = True
      .type = bool
      .help="Extend the binning all the way to detector corner, otherwise to outermost spot on first image"
  }
  dxtbx = False
    .type = bool
    .help="Switch algorithms to dxtbx models"
}

%s
"""%(spotfinder.labelit_related_commands)+additional_spotfinder_phil_defs)


class Empty(object): pass

def run_sweep_strength(image_file_names, params):
  E = Empty()
  E.argv=['Empty']
  E.argv.extend(image_file_names)
  return run_signal_strength_core(params,E)

def run_signal_strength_core(params,E):
  from spotfinder.applications.wrappers import DistlOrganizer
  verbosity = params.distl.verbosity
  if params.distl.res.inner!=None:
    params.distl_lowres_limit = params.distl.res.inner
  if params.distl.res.outer!=None:
    params.force_method2_resolution_limit = params.distl.res.outer
    params.distl_highres_limit = params.distl.res.outer

  params.distl_force_binning = False
  params.distl_permit_binning = False
  params.wedgelimit = len(E.argv)
  params.spotfinder_header_tests = False
  Org = DistlOrganizer(verbose=(verbosity>1), argument_module=E,
                       phil_params=params)
  Org.printSpots()
  return Org

def table(spotfinder_results, keys):
  headers = ['Image#']
  headers.extend(keys)
  rows = [headers]
  for i in sorted(spotfinder_results.images.keys()):
    spotfinder_one_image = spotfinder_results.images[i]
    row = []
    row.append("%i" %i)
    for k in keys:
      value = spotfinder_one_image[k]
      try:
        as_float = float(value)
        as_int = int(value)
        if float(as_int) == as_float:
          row.append("%i" %as_int)
        else:
          row.append("%.2f" %as_float)
      except (ValueError, TypeError):
        row.append("%s" %value)
    rows.append(row)
  return rows

def print_table(spotfinder_result, keys, out=None):
  if out is None:
    import sys
    out = sys.stdout
  from libtbx import table_utils
  rows = table(spotfinder_result, keys)
  print(table_utils.format(rows, has_header=True), file=out)

def as_columns(spotfinder_results):
  d = {}
  keys = (
    'N_spots_total', 'N_spots_non-ice', 'N_spots_resolution',
    'N_spots_unimodal', 'N_spots_inlier',
    #'intensity', 'area', 'neighbor', #'maxcel',
    'resolution', 'distl_resolution', 'ice-ring_impact')

  for i in sorted(spotfinder_results.images.keys()):
    spotfinder_result = spotfinder_results.images[i]
    for k in keys:
      d.setdefault(k, [])
      d[k].append(spotfinder_result[k])
  return d

def plot(spotfinder_results, file_name):
  try:
    import matplotlib
    matplotlib.use('Agg') # use a non-interactive backend
    # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
    from matplotlib import pyplot
  except ImportError:
    raise Sorry("matplotlib must be installed to generate a plot.")
  columns = as_columns(spotfinder_results)
  n_spots = columns.get('N_spots_inlier')
  resolution = columns.get('resolution')
  i_image = range(1, len(n_spots)+1)
  fig = pyplot.figure()
  ax1 = fig.add_subplot(111)
  sc1 = ax1.scatter(i_image, n_spots, s=20, color='blue', marker='o', alpha=0.5)
  ax1.set_xlabel('Image #')
  ax1.set_ylabel('# spots')
  ax1.set_xlim((0, len(n_spots)))
  ax1.set_ylim(bottom=0)
  ax2 = ax1.twinx()
  sc2 = ax2.scatter(i_image, resolution, s=20, color='red', marker='^', alpha=0.5)
  ax2.set_ylabel(u'resolution (\u00c5)')
  ax2.set_xlim((0, len(n_spots)))
  ax2.invert_yaxis()
  lgd = pyplot.legend(
    (sc1, sc2), ('# good spots', 'resolution (method 2)'), ncol=2,
    loc='upper center',
    mode="expand", borderaxespad=0.,
    bbox_to_anchor=(0.0,-0.22, 1., .102))
  pyplot.savefig(file_name, dpi=600, bbox_extra_artists=(lgd,),
                 bbox_inches='tight')

def as_csv(spotfinder_results, out=None):
  if out is None:
    import sys
    out = sys.stdout
  columns = as_columns(spotfinder_results)
  from iotbx import csv_utils
  csv_utils.writer(out, list(columns.values()), field_names=list(columns.keys()))


def run(args, command_name="distl.sweep_strength"):
  help_str="""\
Similar to distl.signal_strength, but acting on a sweep of images, with
tabulation of the results and optional output of results as CSV file and
plots of number of spots and resolution with image number.
"""

  if (len(args) == 0 or args[0] in ["H","h","-H","-h","help","--help","-help"]):
    print("usage:   %s image_prefix_*.img [parameter=value ...]" % command_name)
    print("example: %s lysozyme_*.img distl.minimum_spot_area=8 plot.file_name=lysozyme.pdf"%command_name)
    master_params.show(attributes_level=1,expert_level=1)
    print(help_str)
    return

  print("%s: characterization of candidate Bragg spots"%command_name)

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")
  image_file_names = []
  moving_pdb_file_name = None

  for arg in args:
    if (os.path.isfile(arg)):
      image_file_names.append(arg)
    else:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  if len(image_file_names) < 2:
    raise RuntimeError(
      "Please provide more than one file. Alternatively use "
      "distl.signal_strength to process a single image file.")

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  if params.distl.verbosity > 0:
    print("#Parameters used:")
    print("#phil __ON__")
    print()
    working_params = master_params.format(python_object=params)
    working_params.show(expert_level=1)
    print()
    print("#phil __OFF__")
    print()

  from spotfinder.applications import signal_strength

  spotfinder_results = run_sweep_strength(image_file_names, params)
  print_table(spotfinder_results.S, keys=["N_spots_inlier", "resolution"])

  csv_file_name = params.distl.csv
  if csv_file_name is not None:
    with open(csv_file_name, 'wb') as f:
      as_csv(spotfinder_results.S, out=f)
  plot_file_name = params.distl.plot.file_name
  if plot_file_name is not None:
    plot(spotfinder_results.S, file_name=plot_file_name)


if (__name__ == "__main__"):
  run(args=sys.argv[1:])
