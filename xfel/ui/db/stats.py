from __future__ import division

class Stats(object):
  def __init__(self, app, trial, tags = None):
    self.app = app
    self.trial = trial
    if tags is None:
      tags = []
    self.tags = tags
    self.tag_ids = [t.id for t in tags]

  def __call__(self):
    runs = []
    run_numbers = []
    for rungroup in self.trial.rungroups:
      for run in rungroup.runs:
        if len(self.tags) > 0:
          found_it = False
          for tag in run.tags:
            if tag.id in self.tag_ids:
              found_it = True
              break
          if not found_it:
            continue

        if run.run not in run_numbers:
          runs.append(run)
          run_numbers.append(run.run)

    if len(runs) == 0:
      return []

    runs_str = "(%s)"%(", ".join([str(r.id) for r in runs]))
    exp_tag = self.app.params.experiment_tag
    query = """SELECT DISTINCT(crystal.cell_id) FROM `%s_crystal` crystal
               JOIN `%s_experiment` exp ON exp.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = exp.imageset_id
               JOIN `%s_imageset_event` ie ON ie.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie.event_id
               JOIN `%s_run` run ON run.id = evt.run_id
               WHERE run.id in %s""" % (
      exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, runs_str)
    cell_ids = self.app.execute_query(query).fetchall()
    cells = [self.app.get_cell(cell_id=i[0]) for i in cell_ids]

    for cell in cells:
      for bin in cell.bins:
        query = """SELECT cell_bin.count FROM `%s_cell_bin` cell_bin
                   JOIN `%s_bin` bin ON bin.id = cell_bin.bin_id
                   JOIN `%s_crystal` crystal ON crystal.id = cell_bin.crystal_id
                   JOIN `%s_experiment` exp ON exp.crystal_id = crystal.id
                   JOIN `%s_imageset` imgset ON imgset.id = exp.imageset_id
                   JOIN `%s_imageset_event` ie ON ie.imageset_id = imgset.id
                   JOIN `%s_event` evt ON evt.id = ie.event_id
                   JOIN `%s_run` run ON run.id = evt.run_id
                   WHERE run.id in %s AND bin.id = %d""" % (
          exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, runs_str, int(bin.id))
        counts = self.app.execute_query(query).fetchall()
        bin.count = sum([c[0] for c in counts])
    return cells
