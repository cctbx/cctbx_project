from __future__ import division
from scitbx.array_family import flex

class Stats(object):
  def __init__(self, app, trial, tags = None, isigi_cutoff = None):
    self.app = app
    self.trial = trial
    if tags is None:
      tags = []
    self.tags = tags
    self.tag_ids = [t.id for t in tags]
    self.isigi_cutoff = isigi_cutoff

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
               JOIN `%s_trial` trial ON trial.id = evt.trial_id
               JOIN `%s_run` run ON run.id = evt.run_id
               WHERE run.id in %s AND trial.id = %d""" % (
      exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, runs_str, self.trial.id)
    cell_ids = [str(i[0]) for i in self.app.execute_query(query).fetchall()]
    if len(cell_ids) == 0:
      cells = []
    else:
      from experiment import Cell
      cells = self.app.get_all_x(Cell, 'cell', where = "WHERE id IN (%s)"%", ".join(cell_ids))

    for cell in cells:
      bin_ids = [str(bin.id) for bin in cell.bins]
      if len(bin_ids) == 0:
        continue

      # Big expensive query to avoid many small queries
      query = """SELECT bin.id, cell_bin.count FROM `%s_cell_bin` cell_bin
                 JOIN `%s_bin` bin ON bin.id = cell_bin.bin_id
                 JOIN `%s_crystal` crystal ON crystal.id = cell_bin.crystal_id
                 JOIN `%s_experiment` exp ON exp.crystal_id = crystal.id
                 JOIN `%s_imageset` imgset ON imgset.id = exp.imageset_id
                 JOIN `%s_imageset_event` ie ON ie.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie.event_id
                 JOIN `%s_run` run ON run.id = evt.run_id
                 WHERE run.id in %s AND bin.id IN (%s)
                       AND cell_bin.avg_intensity > 0
                       """ % (
        exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, exp_tag, runs_str, ", ".join(bin_ids))
      if self.isigi_cutoff is not None and self.isigi_cutoff >= 0:
        query += " AND cell_bin.avg_i_sigi >= %f"%self.isigi_cutoff
      results = self.app.execute_query(query).fetchall()
      ids = flex.int([r[0] for r in results])
      counts = flex.int([r[1] for r in results])
      for bin in cell.bins:
        bin.count = flex.sum(counts.select(ids == bin.id))
    return cells

class HitrateStats(object):
  def __init__(self, app, run_number, trial_number, rungroup_id, d_min = None):
    self.app = app
    self.run = app.get_run(run_number = run_number)
    self.trial = app.get_trial(trial_number = trial_number)
    self.rungroup = app.get_rungroup(rungroup_id = rungroup_id)
    self.d_min = d_min

  def __call__(self):
    from iotbx.detectors.cspad_detector_formats import reverse_timestamp
    run_numbers = [r.run for r in self.trial.runs]
    assert self.run.run in run_numbers
    rungroup_ids = [rg.id for rg in self.trial.rungroups]
    assert self.rungroup.id in rungroup_ids
    isoforms = self.trial.isoforms
    assert len(isoforms) > 0
    low_res_bin_ids = []
    high_res_bin_ids = []
    for isoform in isoforms:
      bins = isoform.cell.bins
      d_mins = [float(b.d_min) for b in bins]
      low_res_bin_ids.append(str(bins[d_mins.index(max(d_mins))].id))
      if self.d_min is None:
        min_bin_index = d_mins.index(min(d_mins))
      else:
        d_maxes = [float(b.d_max) for b in bins]
        qualified_bin_indices = [i for i in xrange(len(bins)) if d_maxes[i] >= self.d_min and d_mins[i] <= self.d_min]
        assert len(qualified_bin_indices) == 1
        min_bin_index = qualified_bin_indices[0]
      high_res_bin_ids.append(str(bins[min_bin_index].id))
    assert len(low_res_bin_ids) > 0
    assert len(high_res_bin_ids) > 0
    assert len(low_res_bin_ids) == len(high_res_bin_ids)

    tag = self.app.params.experiment_tag

    # Get the high and low res avg_i_sigi in one query. Means there will be 2x timestamps retrieved, where each is found twice
    query = """SELECT bin.id, event.timestamp, event.n_strong, cb.avg_i_sigi
               FROM `%s_event` event
               JOIN `%s_imageset_event` is_e ON is_e.event_id = event.id
               JOIN `%s_imageset` imgset ON imgset.id = is_e.imageset_id
               JOIN `%s_experiment` exp ON exp.imageset_id = imgset.id
               JOIN `%s_crystal` crystal ON crystal.id = exp.crystal_id
               JOIN `%s_cell` cell ON cell.id = crystal.cell_id
               JOIN `%s_bin` bin ON bin.cell_id = cell.id
               JOIN `%s_cell_bin` cb ON cb.bin_id = bin.id AND cb.crystal_id = crystal.id
               WHERE event.trial_id = %d AND event.run_id = %d AND event.rungroup_id = %d AND
                     cb.bin_id IN (%s)
            """ % (tag, tag, tag, tag, tag, tag, tag, tag, self.trial.id, self.run.id, self.rungroup.id, ", ".join(low_res_bin_ids + high_res_bin_ids))
    cursor = self.app.execute_query(query)
    timestamps = flex.double()
    n_strong = flex.int()
    average_i_sigi_low = flex.double()
    average_i_sigi_high = flex.double()
    for row in cursor.fetchall():
      b_id, ts, n_s, avg_i_sigi = row
      rts = reverse_timestamp(ts)
      rts = rts[0] + (rts[1]/1000)
      if rts not in timestamps:
        # First time through, figure out which bin is reported (high or low), add avg_i_sigi to that set of results
        timestamps.append(rts)
        n_strong.append(n_s)
        if str(b_id) in low_res_bin_ids:
          if avg_i_sigi is None:
            average_i_sigi_low.append(1e-6)
          else:
            average_i_sigi_low.append(avg_i_sigi)
          average_i_sigi_high.append(0)
        elif str(b_id) in high_res_bin_ids:
          average_i_sigi_low.append(0)
          if avg_i_sigi is None:
            average_i_sigi_high.append(0)
          else:
            average_i_sigi_high.append(avg_i_sigi)
        else:
          assert False
      else:
        # Second time through, already have added to timestamps and n_strong, so fill in missing avg_i_sigi
        index = flex.first_index(timestamps, rts)
        if str(b_id) in low_res_bin_ids:
          average_i_sigi_low[index] = avg_i_sigi
        elif str(b_id) in high_res_bin_ids:
          if avg_i_sigi is None:
            average_i_sigi_high[index] = 0
          else:
            average_i_sigi_high[index] = avg_i_sigi
        else:
          assert False

    # This left join query finds the events with no imageset, meaning they failed to index
    query = """SELECT event.timestamp, event.n_strong
               FROM `%s_event` event
               LEFT JOIN `%s_imageset_event` is_e ON is_e.event_id = event.id
               WHERE is_e.event_id IS NULL AND
                     event.trial_id = %d AND event.run_id = %d AND event.rungroup_id = %d
            """ % (tag, tag, self.trial.id, self.run.id, self.rungroup.id)

    cursor = self.app.execute_query(query)
    for row in cursor.fetchall():
      ts, n_s, = row
      rts = reverse_timestamp(ts)
      timestamps.append(rts[0] + (rts[1]/1000))
      n_strong.append(n_s)
      average_i_sigi_low.append(0)
      average_i_sigi_high.append(0)

    order = flex.sort_permutation(timestamps)
    timestamps = timestamps.select(order)
    n_strong = n_strong.select(order)
    average_i_sigi_low = average_i_sigi_low.select(order)
    average_i_sigi_high = average_i_sigi_high.select(order)

    return timestamps, n_strong, average_i_sigi_low, average_i_sigi_high
