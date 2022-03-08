from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex

class Stats(object):
  def __init__(self, app, trial, tags = None, isigi_cutoff = None, tag_selection_mode="union", selected_runs = None, selected_rungroup = None):
    self.app = app
    self.trial = trial
    if tags is None:
      tags = []
    self.tags = tags
    self.tag_ids = [t.id for t in tags]
    self.isigi_cutoff = isigi_cutoff
    self.tag_selection_mode = tag_selection_mode
    self.selected_runs = selected_runs
    self.selected_rungroup = selected_rungroup

  def __call__(self):
    runs = []
    run_numbers = []
    if self.selected_runs is not None:
      assert self.selected_rungroup is not None
      selected_run_ids = [r.id for r in self.selected_runs]
    for rungroup in self.trial.rungroups:
      if self.selected_rungroup is not None and self.selected_rungroup.id != rungroup.id:
        continue
      for run in rungroup.runs:
        if self.selected_runs is not None:
          if run.id not in selected_run_ids:
            continue
        if len(self.tags) > 0:
          tags_found = []
          for tag in run.tags:
            if tag.id in self.tag_ids:
              tags_found.append(tag.id)
          if self.tag_selection_mode == "union":
            if len(tags_found) == 0:
              continue
          elif self.tag_selection_mode == "intersection":
            if len(tags_found) < len(self.tag_ids):
              continue
          else:
            assert False

        if run.run not in run_numbers:
          runs.append(run)
          run_numbers.append(run.run)

    if len(runs) == 0:
      return []

    runs_str = "(%s)"%(", ".join([str(r.id) for r in runs]))
    tag = self.app.params.experiment_tag

    # Big expensive query to avoid many small queries
    query = """SELECT cell.id, bin.id, SUM(cell_bin.count) FROM `%s_cell_bin` cell_bin
               JOIN `%s_bin` bin ON bin.id = cell_bin.bin_id
               JOIN `%s_cell` cell ON cell.id = bin.cell_id
               JOIN `%s_crystal` crystal ON crystal.id = cell_bin.crystal_id
               JOIN `%s_experiment` exp ON exp.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = exp.imageset_id
               JOIN `%s_imageset_event` ie ON ie.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie.event_id
               JOIN `%s_run` run ON run.id = evt.run_id
               JOIN `%s_rungroup` rg ON rg.id = evt.rungroup_id
               JOIN `%s_trial_rungroup` t_rg on t_rg.rungroup_id = rg.id
               JOIN `%s_trial` trial ON trial.id = t_rg.trial_id
               WHERE run.id IN %s
                     AND cell_bin.avg_intensity > 0
                     AND trial.id = %d
                     AND rg.active = True
                     """ % (
      tag, tag, tag, tag, tag, tag, tag, tag, tag, tag, tag, tag, runs_str, self.trial.id)

    if self.isigi_cutoff is not None and self.isigi_cutoff >= 0:
      query += " AND cell_bin.avg_i_sigi >= %f"%self.isigi_cutoff
    query += " GROUP BY cell.id, bin.id"

    results = self.app.execute_query(query).fetchall()
    if len(results) == 0:
      return []
    cell_ids = set([str(r[0]) for r in results])
    from .experiment import Cell
    cells = self.app.get_all_x(Cell, 'cell', where = "WHERE id IN (%s)"%", ".join(cell_ids))
    self.app.link_cell_bins(cells)

    cell_bins_d = {}
    for cell in cells:
      cell_bins_d[cell.id] = {}
      for bin in cell.bins:
        cell_bins_d[cell.id][bin.id] = bin

    for cell_id, bin_id, count in results:
      cell_bins_d[cell_id][bin_id].count = count

    return cells

class HitrateStats(object):
  def __init__(self, app, run_number, trial_number, rungroup_id, d_min = None, i_sigi_cutoff = 1, raw_data_sampling = 1):
    self.app = app
    self.run = app.get_run(run_number = run_number)
    self.trial = app.get_trial(trial_number = trial_number)
    self.rungroup = app.get_rungroup(rungroup_id = rungroup_id)
    self.d_min = d_min
    self.i_sigi_cutoff = i_sigi_cutoff
    self.sampling = raw_data_sampling

  def __call__(self):
    from iotbx.detectors.cspad_detector_formats import reverse_timestamp
    from xfel.ui.components.timeit import duration
    #import time
    #t1 = time.time()
    run_numbers = [r.run for r in self.trial.runs]
    assert self.run.run in run_numbers
    rungroup_ids = [rg.id for rg in self.trial.rungroups]
    assert self.rungroup.id in rungroup_ids
    if len(self.trial.isoforms) > 0:
      cells = [isoform.cell for isoform in self.trial.isoforms]
    else:
      cells = self.app.get_trial_cells(self.trial.id, self.rungroup.id, self.run.id)

    high_res_bin_ids = []
    for cell in cells:
      bins = cell.bins
      d_mins = [float(b.d_min) for b in bins]
      if len(d_mins) == 0: continue
      if self.d_min is None:
        min_bin_index = d_mins.index(min(d_mins))
      else:
        d_maxes = [float(b.d_max) for b in bins]
        qualified_bin_indices = [i for i in range(len(bins)) if d_maxes[i] >= self.d_min and d_mins[i] <= self.d_min]
        if len(qualified_bin_indices) == 0: continue
        min_bin_index = qualified_bin_indices[0]
      high_res_bin_ids.append(str(bins[min_bin_index].id))

    resolutions = flex.double()
    two_theta_low = flex.double()
    two_theta_high = flex.double()
    tag = self.app.params.experiment_tag
    timestamps, timestamps_s = flex.double(), []
    n_strong = flex.int()
    n_lattices = flex.int()
    if len(high_res_bin_ids) > 0:

      # Get the stats in one query.
      query = """SELECT event.timestamp, event.n_strong, MIN(bin.d_min), event.two_theta_low, event.two_theta_high, COUNT(DISTINCT crystal.id)
                 FROM `%s_event` event
                 JOIN `%s_imageset_event` is_e ON is_e.event_id = event.id
                 JOIN `%s_imageset` imgset ON imgset.id = is_e.imageset_id
                 JOIN `%s_experiment` exp ON exp.imageset_id = imgset.id
                 JOIN `%s_crystal` crystal ON crystal.id = exp.crystal_id
                 JOIN `%s_cell` cell ON cell.id = crystal.cell_id
                 JOIN `%s_bin` bin ON bin.cell_id = cell.id
                 JOIN `%s_cell_bin` cb ON cb.bin_id = bin.id AND cb.crystal_id = crystal.id
                 WHERE event.trial_id = %d AND event.run_id = %d AND event.rungroup_id = %d AND
                       cb.avg_i_sigi >= %f
                 GROUP BY event.id
              """ % (tag, tag, tag, tag, tag, tag, tag, tag, self.trial.id, self.run.id, self.rungroup.id, self.i_sigi_cutoff)
      cursor = self.app.execute_query(query)
      sample = -1
      for row in cursor.fetchall():
        sample += 1
        if sample % self.sampling != 0:
          continue
        ts, n_s, d_min, tt_low, tt_high, n_xtal = row
        try:
          d_min = float(d_min)
        except ValueError:
          d_min = None
        try:
          rts = reverse_timestamp(ts)
          timestamps.append(rts[0] + (rts[1]/1000))
        except ValueError:
          try:
            timestamps.append(float(ts))
          except ValueError:
            timestamps_s.append(ts)
        n_strong.append(n_s)
        two_theta_low.append(tt_low or -1)
        two_theta_high.append(tt_high or -1)
        resolutions.append(d_min or 0)
        n_lattices.append(n_xtal or 0)

    # only get results that are strings or ints, not a mix of both
    assert not (len(timestamps) > 0 and len(timestamps_s) > 0)

    # This left join query finds the events with no imageset, meaning they failed to index
    query = """SELECT event.timestamp, event.n_strong, event.two_theta_low, event.two_theta_high
               FROM `%s_event` event
               LEFT JOIN `%s_imageset_event` is_e ON is_e.event_id = event.id
               WHERE is_e.event_id IS NULL AND
                     event.trial_id = %d AND event.run_id = %d AND event.rungroup_id = %d
            """ % (tag, tag, self.trial.id, self.run.id, self.rungroup.id)

    cursor = self.app.execute_query(query)
    for row in cursor.fetchall():
      ts, n_s, tt_low, tt_high = row
      try:
        rts = reverse_timestamp(ts)
        timestamps.append(rts[0] + (rts[1]/1000))
      except ValueError:
        try:
          rts = float(ts)
          timestamps.append(rts)
        except ValueError:
          timestamps_s.append(ts)
      n_strong.append(n_s)
      two_theta_low.append(tt_low or -1)
      two_theta_high.append(tt_high or -1)
      resolutions.append(0)
      n_lattices.append(0)

    if len(timestamps_s) > 0:
      timestamps = flex.double([i[0] for i in sorted(enumerate(timestamps_s), key=lambda x:x[1])])
      order = flex.size_t([i for i in timestamps.iround()])
      timestamps = flex.sorted(timestamps)
    else:
      order = flex.sort_permutation(timestamps)
      timestamps = timestamps.select(order)
    n_strong = n_strong.select(order)
    two_theta_low = two_theta_low.select(order)
    two_theta_high = two_theta_high.select(order)
    resolutions = resolutions.select(order)
    n_lattices = n_lattices.select(order)

    #t2 = time.time()
    #print "HitrateStats took %s" % duration(t1, t2)
    return timestamps, two_theta_low, two_theta_high, n_strong, resolutions, n_lattices

class SpotfinderStats(object):
  def __init__(self, app, run_number, trial_number, rungroup_id, raw_data_sampling = 1):
    self.app = app
    self.run = app.get_run(run_number = run_number)
    self.trial = app.get_trial(trial_number = trial_number)
    self.rungroup = app.get_rungroup(rungroup_id = rungroup_id)
    self.sampling = raw_data_sampling

  def __call__(self):
    from iotbx.detectors.cspad_detector_formats import reverse_timestamp
    from xfel.ui.components.timeit import duration
    import time
    t1 = time.time()
    run_numbers = [r.run for r in self.trial.runs]
    assert self.run.run in run_numbers
    rungroup_ids = [rg.id for rg in self.trial.rungroups]
    assert self.rungroup.id in rungroup_ids
    if len(self.trial.isoforms) > 0:
      cells = [isoform.cell for isoform in self.trial.isoforms]
    else:
      cells = self.app.get_trial_cells(self.trial.id, self.rungroup.id, self.run.id)

    low_res_bin_ids = []
    for cell in cells:
      bins = cell.bins
      d_mins = [float(b.d_min) for b in bins]
      if len(d_mins) == 0: continue
      low_res_bin_ids.append(str(bins[d_mins.index(max(d_mins))].id))

    tag = self.app.params.experiment_tag
    timestamps = flex.double()
    xtal_ids = flex.double()
    n_strong = flex.int()
    if len(low_res_bin_ids) > 0:

      # Get the spotfinding results from the selected runs
      query = """SELECT bin.id, crystal.id, event.timestamp, event.n_strong
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
              """ % (tag, tag, tag, tag, tag, tag, tag, tag, self.trial.id, self.run.id, self.rungroup.id,
                    ", ".join(low_res_bin_ids))
      cursor = self.app.execute_query(query)
      sample = -1
      for row in cursor.fetchall():
        b_id, xtal_id, ts, n_s = row
        try:
          rts = reverse_timestamp(ts)
          rts = rts[0] + (rts[1]/1000)
        except ValueError:
          rts = float(ts)
        if xtal_id not in xtal_ids:
          sample += 1
          if sample % self.sampling != 0:
            continue
          timestamps.append(rts)
          xtal_ids.append(xtal_id)
          n_strong.append(n_s)

    # This left join query finds the events with no imageset, meaning they failed to index
    query = """SELECT event.timestamp, event.n_strong
               FROM `%s_event` event
               LEFT JOIN `%s_imageset_event` is_e ON is_e.event_id = event.id
               WHERE is_e.event_id IS NULL AND
                     event.trial_id = %d AND event.run_id = %d AND event.rungroup_id = %d
            """ % (tag, tag, self.trial.id, self.run.id, self.rungroup.id)

    cursor = self.app.execute_query(query)
    for row in cursor.fetchall():
      ts, n_s = row
      try:
        rts = reverse_timestamp(ts)
        timestamps.append(rts[0] + (rts[1]/1000))
      except ValueError:
        rts = float(ts)
        timestamps.append(rts)
      n_strong.append(n_s)

    order = flex.sort_permutation(timestamps)
    timestamps = timestamps.select(order)
    n_strong = n_strong.select(order)

    t2 = time.time()
    return timestamps, n_strong
