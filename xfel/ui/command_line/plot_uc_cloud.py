from __future__ import division
from __future__ import print_function

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.db.xfel_db import xfel_db_application
import sys
from xfel.ui import db_phil_str

phil_str = """
  trial = None
    .type = int
  tag = None
    .type = str
    .multiple = True
  tag_selection_mode = *union intersection
    .type = choice
  run = None
    .type = str
    .multiple = True
  rungroup = None
    .type = int
  iqr_ratio = 1.5
    .type = float
    .help = Interquartile range multiplier for outlier rejection. Use None to disable outlier rejection.
"""
phil_scope = parse(phil_str + db_phil_str)

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  app = xfel_db_application(params)

  if params.tag is not None and len(params.tag) > 0:
    tags = []
    for tag in app.get_all_tags():
      for t in params.tag:
        if t == tag.name:
          tags.append(tag)
    extra_title = ",".join([t.name for t in tags])
  else:
    tags = None
    extra_title = None

  if params.run is not None and len(params.run) > 0:
    assert params.rungroup is not None
    runs = [app.get_run(run_number = r) for r in params.run]
    rungroup = app.get_rungroup(rungroup_id = params.rungroup)
  else:
    runs = None
    rungroup = None

  if extra_title is None and runs is not None:
    extra_title = "%s"%(",".join(["%s"%r.run for r in runs]))

  trial = app.get_trial(trial_number=params.trial)
  info = []
  print("Reading data...")
  cells = app.get_stats(trial=trial, tags=tags, isigi_cutoff = 1.0, tag_selection_mode = params.tag_selection_mode, selected_runs = runs, selected_rungroup = rungroup)()
  for cell in cells:
    info.append({'a':cell.cell_a,
                 'b':cell.cell_b,
                 'c':cell.cell_c,
                 'alpha':cell.cell_alpha,
                 'beta':cell.cell_beta,
                 'gamma':cell.cell_gamma,
                 'n_img':0})
  import xfel.ui.components.xfel_gui_plotter as pltr
  plotter = pltr.PopUpCharts()
  plotter.plot_uc_histogram(info_list=[info], extra_title=extra_title, legend_list=[""], iqr_ratio = params.iqr_ratio)
  plotter.plot_uc_3Dplot(info=info, iqr_ratio = params.iqr_ratio)
  plotter.plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
