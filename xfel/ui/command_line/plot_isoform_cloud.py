from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.db.xfel_db import xfel_db_application
import sys
from xfel.ui import db_phil_str

"""
Script to compare whether a given trial is sorting results into the proper isoforms. Requires two trials, one with an isoforms and one without
"""

phil_str = """
  trial_noisoforms = None
    .type = int
  trial_isoforms = None
    .type = int
  tag = None
    .type = str
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

  trial_ni = app.get_trial(trial_number=params.trial_noisoforms)
  trial_i = app.get_trial(trial_number=params.trial_isoforms)

  tags = []
  if params.tag is not None and len(params.tag) > 0:
    for tag in app.get_all_tags():
      if params.tag == tag.name:
        tags.append(tag)
    extra_title = ",".join([t.name for t in tags])
    tag_ids = [t.id for t in tags]
    runs_ni = [r for r in trial_ni.runs if any([t.id in tag_ids for t in r.tags])]
    runs_i = [r for r in trial_i.runs if any([t.id in tag_ids for t in r.tags])]
  else:
    extra_title = None
    runs_ni = trial_ni.runs
    runs_i = trial_i.runs

  runs_ni_str = "(%s)"%(", ".join([str(r.id) for r in runs_ni]))
  exp_tag = params.experiment_tag

  isoforms = app.get_trial_isoforms(trial_i.id)

  print("Reading data...")
  info_list = []
  legend_list = []
  for isoform in isoforms:
    events_i = app.get_all_events(trial = trial_i, runs = runs_i, isoform=isoform)
    if len(events_i) == 0:
      print("No events of isoform %s found"%isoform.name)
      continue
    legend_list.append(isoform.name)
    events_i_str = "(%s)"%(", ".join(["'%s'"%e.timestamp for e in events_i]))

    query = """SELECT crystal.cell_id FROM `%s_crystal` crystal
               JOIN `%s_experiment` exp ON exp.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = exp.imageset_id
               JOIN `%s_imageset_event` ie ON ie.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie.event_id
               JOIN `%s_trial` trial ON trial.id = evt.trial_id
               JOIN `%s_run` run ON run.id = evt.run_id
               WHERE run.id in %s AND trial.id = %d AND evt.timestamp IN %s""" % (
      exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, runs_ni_str, trial_ni.id, events_i_str)
    cell_ids = [str(i[0]) for i in app.execute_query(query).fetchall()]
    if len(cell_ids) == 0:
      cells = []
    else:
      from xfel.ui.db.experiment import Cell
      cells = app.get_all_x(Cell, 'cell', where = "WHERE id IN (%s)"%", ".join(cell_ids))

    print(len(cells))
    info = []
    for cell in cells:
      info.append({'a':cell.cell_a,
                   'b':cell.cell_b,
                   'c':cell.cell_c,
                   'alpha':cell.cell_alpha,
                   'beta':cell.cell_beta,
                   'gamma':cell.cell_gamma,
                   'n_img':0})
    info_list.append(info)

  import xfel.ui.components.xfel_gui_plotter as pltr
  plotter = pltr.PopUpCharts()
  plotter.plot_uc_histogram(info_list=info_list, legend_list=legend_list, extra_title=extra_title)
  plotter.plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
