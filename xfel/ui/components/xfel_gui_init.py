from __future__ import absolute_import, division, print_function
from six.moves import range, zip, map

'''
Author      : Lyubimov, A.Y.
Created     : 06/02/2016
Last Changed: 02/09/2017
Description : XFEL UI Initialization module
'''

import os
import wx
import time
import numpy as np
from threading import Thread
from wx.lib.scrolledpanel import ScrolledPanel
from libtbx import easy_run

try:
  from MySQLdb import OperationalError # test import
except ImportError:
  from libtbx.utils import Sorry
  raise Sorry('Mysql not available')

from xfel.clustering.cluster import Cluster
import xfel.ui.components.xfel_gui_controls as gctr
import xfel.ui.components.xfel_gui_dialogs as dlg
import xfel.ui.components.xfel_gui_plotter as pltr

from xfel.ui import load_cached_settings, save_cached_settings
from xfel.ui.db import get_run_path
from xfel.ui.db.xfel_db import xfel_db_application

from prime.postrefine.mod_gui_frames import PRIMEInputWindow, PRIMERunWindow
from prime.postrefine.mod_input import master_phil
from iota.utils.utils import Capturing, set_base_dir

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')

import getpass
user = getpass.getuser()

import libtbx.load_env
license = os.path.join(libtbx.env.find_in_repositories('cctbx_project'), 'LICENSE.txt')

description = 'The cctbx.xfel UI is developed for use during data collection ' \
              'and initial processing of serial crystallographic data from' \
              'XFELs and synchrotrons.'

class TagSet(object):
  def __init__(self, tag_selection_mode, tags):
    assert tag_selection_mode in ['union', 'intersection']
    self.mode = tag_selection_mode
    self.tags = tags
  def __str__(self):
    if len(self.tags) > 1:
      return ", ".join([t.name for t in self.tags]) + (' (%s)' % self.mode[0])
    else:
      return ", ".join([t.name for t in self.tags])


# ------------------------------- Run Sentinel ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_RUN_REFRESH = wx.NewEventType()
EVT_RUN_REFRESH = wx.PyEventBinder(tp_EVT_RUN_REFRESH, 1)

class RefreshRuns(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class RunSentinel(Thread):
  ''' Worker thread for runs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

    if self.parent.params.facility.name == 'standalone':
      if self.parent.params.facility.standalone.monitor_for == 'folders' and \
         self.parent.params.facility.standalone.folders.method == 'status_file':
        from xfel.ui.db.xfel_db import cheetah_run_finder
        self.finder = cheetah_run_finder(self.parent.params)
      else:
        from xfel.ui.db.xfel_db import standalone_run_finder
        self.finder = standalone_run_finder(self.parent.params)

  def post_refresh(self):
    evt = RefreshRuns(tp_EVT_RUN_REFRESH, -1)
    wx.PostEvent(self.parent.run_window.runs_tab, evt)
    wx.PostEvent(self.parent.run_window.trials_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    db = xfel_db_application(self.parent.params)
    use_ids = self.parent.params.facility.name not in ['lcls']

    while self.active:
      # Find the delta
      known_runs = [r.run for r in db.get_all_runs()]
      if self.parent.params.facility.name == 'lcls':
        unknown_run_runs = [str(run['run']) for run in db.list_lcls_runs() if
                            str(run['run']) not in known_runs]
        unknown_run_paths = [''] * len(unknown_run_runs)
      elif self.parent.params.facility.name == 'standalone':
        standalone_runs = [run for run in self.finder.list_runs() if
                           run[0] not in known_runs]
        unknown_run_runs = [r[0] for r in standalone_runs]
        unknown_run_paths = [r[1] for r in standalone_runs]

      if len(unknown_run_runs) > 0:
        for run_run, run_path in zip(unknown_run_runs, unknown_run_paths):
          db.create_run(run = run_run, path = run_path)
        new_runs = [r for r in db.get_all_runs() if r.run in unknown_run_runs]
        if len(self.parent.run_window.runs_tab.persistent_tags) > 0:
          tags = [t for t in db.get_all_tags() if t.name in self.parent.run_window.runs_tab.persistent_tags]
          for r in new_runs:
            for t in tags:
              r.add_tag(t)
        # Sync new runs to rungroups
        for rungroup in db.get_all_rungroups(only_active=True):
          first_run, last_run = rungroup.get_first_and_last_runs()
          # HACK: to get working -- TODO: make nice
          if use_ids:
            first_run = first_run.id
            last_run = last_run.id if last_run is not None else None
          else:
            first_run = int(first_run.run)
            last_run = int(last_run.run) if last_run is not None else None
          rungroup.sync_runs(first_run, last_run, use_ids=use_ids)

        print("%d new runs" % len(unknown_run_runs))
        self.post_refresh()
      time.sleep(10)

# ------------------------------- Job Monitor ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_JOB_MONITOR = wx.NewEventType()
EVT_JOB_MONITOR = wx.PyEventBinder(tp_EVT_JOB_MONITOR, 1)

class MonitorJobs(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid, trials = None, jobs = None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.trials = trials
    self.jobs = jobs

class JobMonitor(Thread):
  ''' Monitor thread for jobs; generated so that the GUI does not lock up when
      monitoring is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.only_active_jobs = True

  def post_refresh(self, trials = None, jobs = None):
    evt = MonitorJobs(tp_EVT_JOB_MONITOR, -1, trials, jobs)
    wx.PostEvent(self.parent.run_window.jobs_tab, evt)

  def run(self):
    from xfel.ui.components.submission_tracker import TrackerFactory

    # one time post for an initial update
    self.post_refresh()

    db = xfel_db_application(self.parent.params)
    tracker = TrackerFactory.from_params(self.parent.params)

    while self.active:
      self.parent.run_window.jmn_light.change_status('idle')

      trials = db.get_all_trials()
      jobs = db.get_all_jobs(active = self.only_active_jobs)

      for job in jobs:
        if job.status in ['DONE', 'EXIT', 'SUBMIT_FAIL', 'DELETED']:
          continue
        new_status = tracker.track(job.submission_id, job.get_log_path())
        # Handle the case where the job was submitted but no status is available yet
        if job.status == "SUBMITTED" and new_status == "ERR":
          pass
        elif job.status != new_status:
          job.status = new_status

      self.post_refresh(trials, jobs)
      self.parent.run_window.jmn_light.change_status('on')
      time.sleep(5)


# ------------------------------- Job Sentinel ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_JOB_REFRESH = wx.NewEventType()
EVT_JOB_REFRESH = wx.PyEventBinder(tp_EVT_JOB_REFRESH, 1)

class RefreshJobs(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class JobSentinel(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    pass

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    db = xfel_db_application(self.parent.params)

    from xfel.ui.db.job import submit_all_jobs

    while self.active:
      submit_all_jobs(db)
      self.post_refresh()
      time.sleep(2)

# ----------------------------- Progress Sentinel ---------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_PRG_REFRESH = wx.NewEventType()
EVT_PRG_REFRESH = wx.PyEventBinder(tp_EVT_PRG_REFRESH, 1)

class RefreshStats(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class ProgressSentinel(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.output = self.parent.params.output_folder
    self.number_of_pickles = 0
    self.info = {}
    self.noiso_cells = []

    # on initialization (and restart), make sure stats drawn from scratch
    self.parent.run_window.status_tab.redraw_windows = True

  def post_refresh(self):
    evt = RefreshStats(tp_EVT_PRG_REFRESH, -1, self.info)
    wx.PostEvent(self.parent.run_window.status_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    db = xfel_db_application(self.parent.params)

    while self.active:
      self.parent.run_window.prg_light.change_status('idle')

      if len(db.get_all_trials()) > 0:
        trial = db.get_trial(
          trial_number=self.parent.run_window.status_tab.trial_no)

        trial_has_isoforms = len(trial.isoforms) > 0

        tags = self.parent.run_window.status_tab.selected_tags
        tag_ids = [tag.id for tag in tags]
        cells = db.get_stats(trial=trial, tags=tags, isigi_cutoff = self.parent.run_window.status_tab.isigi_cutoff)()

        if self.parent.run_window.status_tab.tag_trial_changed:
          self.parent.run_window.status_tab.redraw_windows = True
          self.parent.run_window.status_tab.tag_trial_changed = False

        run_numbers = []
        runs = []
        for rb in trial.rungroups:
          for run in rb.runs:
            if run.run not in run_numbers:
              if len(tags) > 0:
                for tag in run.tags:
                  if tag.id in tag_ids:
                    run_numbers.append(run.run)
                    runs.append(run)
              else:
                run_numbers.append(run.run)
                runs.append(run)
        if not trial_has_isoforms:
          n_img = len(db.get_all_events(trial, runs))

        for cell in cells:
          # Check for cell isoform
          if cell.isoform is None:
            self.noiso_cells.append({'a':cell.cell_a,
                                     'b':cell.cell_b,
                                     'c':cell.cell_c,
                                     'alpha':cell.cell_alpha,
                                     'beta':cell.cell_beta,
                                     'gamma':cell.cell_gamma,
                                     'n_img':n_img})
          else:
            current_rows = self.parent.run_window.status_tab.rows
            if current_rows != {}:
              if cell.isoform._db_dict['name'] in current_rows:
                bins = cell.bins[
                       :int(current_rows[cell.isoform._db_dict['name']]['high_bin'])]
                highest_bin = cell.bins[int(current_rows[cell.isoform._db_dict['name']]['high_bin'])]
              else:
                assert False, "This isoform is not available yet"
            else:
              bins = cell.bins
              d_mins = [b.d_min for b in bins]
              highest_bin = bins[d_mins.index(min(d_mins))]

            counts_all = [int(i.count) for i in bins]
            totals_all = [int(i.total_hkl) for i in bins]
            counts_highest = int(highest_bin.count)
            totals_highest = int(highest_bin.total_hkl)

            # Apply throttle to multiplicity calculation
            if trial.process_percent is None:
              process_percent = 100
            else:
              process_percent = trial.process_percent

            n_img = len(db.get_all_events(trial, runs, isoform = cell.isoform))

            # Generate multiplicity graph for isoforms
            mult_all = sum(counts_all) / sum(totals_all) / (process_percent / 100)
            mult_highest = counts_highest / totals_highest / (process_percent / 100)
            self.info[cell.isoform._db_dict['name']] = {'multiplicity_all':mult_all,
                                            'multiplicity_highest':mult_highest,
                                            'bins':bins,
                                            'isoform':cell.isoform._db_dict['name'],
                                            'a':cell.cell_a,
                                            'b':cell.cell_b,
                                            'c':cell.cell_c,
                                            'alpha':cell.cell_alpha,
                                            'beta':cell.cell_beta,
                                            'gamma':cell.cell_gamma,
                                            'n_img':n_img}
        #if len(self.noiso_cells) > 0:
        if len(self.info) == 0 and len(self.noiso_cells) > 0:
          sum_n_img = sum([cell['n_img'] for cell in self.noiso_cells])
          mean_a = sum([cell['n_img']*cell['a'] for cell in self.noiso_cells])/sum_n_img
          mean_b = sum([cell['n_img']*cell['b'] for cell in self.noiso_cells])/sum_n_img
          mean_c = sum([cell['n_img']*cell['c'] for cell in self.noiso_cells])/sum_n_img
          mean_alpha = sum([cell['n_img']*cell['alpha'] for cell in self.noiso_cells])/sum_n_img
          mean_beta  = sum([cell['n_img']*cell['beta']  for cell in self.noiso_cells])/sum_n_img
          mean_gamma = sum([cell['n_img']*cell['gamma'] for cell in self.noiso_cells])/sum_n_img
          noiso_entry = {'multiplicity_all':0,
                         'multiplicity_highest':0,
                         'bins':None,
                         'isoform':None,
                         'a':mean_a,
                         'b':mean_b,
                         'c':mean_c,
                         'alpha':mean_alpha,
                         'beta':mean_beta,
                         'gamma':mean_gamma,
                         'n_img':sum_n_img}
          self.info['noiso'] = noiso_entry
      self.post_refresh()
      self.info = {}
      self.parent.run_window.prg_light.change_status('on')
      time.sleep(5)

# ----------------------------- Run Stats Sentinel ---------------------------- #

# Set up events for monitoring hitrate, indexing rate and I/sig(I)
tp_EVT_RUNSTATS_REFRESH = wx.NewEventType()
EVT_RUNSTATS_REFRESH = wx.PyEventBinder(tp_EVT_RUNSTATS_REFRESH, 1)

class RefreshRunStats(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class RunStatsSentinel(Thread):
  ''' Worker thread for run stats; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.output = self.parent.params.output_folder
    self.number_of_pickles = 0
    self.info = {}
    self.run_numbers = []
    self.stats = []
    self.run_tags = []
    self.run_statuses = []

  def post_refresh(self):
    evt = RefreshRunStats(tp_EVT_RUNSTATS_REFRESH, -1, self.info)
    wx.PostEvent(self.parent.run_window.runstats_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    self.db = xfel_db_application(self.parent.params)

    while self.active:
      self.parent.run_window.runstats_light.change_status('idle')
      self.plot_stats()
      self.fetch_timestamps(indexed=True)
      self.fetch_timestamps(indexed=False)
      self.post_refresh()
      self.info = {}
      self.parent.run_window.runstats_light.change_status('on')
      time.sleep(5)

  def refresh_stats(self):
    #from xfel.ui.components.timeit import duration
    from xfel.ui.db.stats import HitrateStats
    import copy, time
    t1 = time.time()
    if self.parent.run_window.runstats_tab.trial_no is not None:
      trial = self.db.get_trial(
        trial_number=self.parent.run_window.runstats_tab.trial_no)
      selected_runs = copy.deepcopy(self.parent.run_window.runstats_tab.selected_runs)
      self.run_numbers = []
      trial_ids = []
      rungroup_ids = []
      self.stats = []
      self.trgr = {}
      self.run_tags = []
      self.run_statuses = []
      for rg in trial.rungroups:
        for run in rg.runs:
          if run.run not in self.run_numbers and run.run in selected_runs:
            self.run_numbers.append(run.run)
            trial_ids.append(trial.id)
            rungroup_ids.append(rg.id)
            self.trgr[run.run] = (trial, rg, run)
            self.stats.append(HitrateStats(self.db, run.run, trial.trial, rg.id,
                                           i_sigi_cutoff=self.parent.run_window.runstats_tab.i_sigi,
                                           d_min=self.parent.run_window.runstats_tab.d_min)())
            self.run_tags.append([tag.name for tag in run.tags])

      jobs = self.db.get_all_jobs()
      for idx in range(len(self.run_numbers)):
        run_no = self.run_numbers[idx]
        rg_id = rungroup_ids[idx]
        t_id = trial_ids[idx]
        found_it = False
        for job in jobs:
          try:
            ok = job.run.run == run_no and job.rungroup.id == rg_id and job.trial.id == t_id
          except AttributeError:
            pass
          else:
            if ok:
              self.run_statuses.append(job.status)
              found_it = True; break
        if not found_it: self.run_statuses.append('UNKWN')
    t2 = time.time()

  def get_xtc_process_params_for_run(self, trial, rg, run):
    params = {}
    params['experiment'] = str(self.parent.db.params.facility.lcls.experiment)
    try:
      params['output_dir'] = os.path.join(str(self.parent.db.params.output_folder),
        "r%04d"%(int(run.run)), "%03d_rg%03d"%(trial.trial, rg.rungroup_id), "all")
    except ValueError:
      params['output_dir'] = os.path.join(str(self.parent.db.params.output_folder),
        "r%s"%(run.run), "%03d_rg%03d"%(trial.trial, rg.rungroup_id), "all")
    params['run'] = run.run
    params['address'] = rg.detector_address
    params['format'] = rg.format
    params['config'] = rg.config_str if hasattr(rg, 'config_str') else None
    params['beamx'] = rg.beamx if hasattr(rg, 'beamx') else None
    params['beamy'] = rg.beamy if hasattr(rg, 'beamy') else None
    params['distance'] = rg.detz_parameter if hasattr(rg, 'detz_parameter') else None
    params['bin_size'] = rg.binning if hasattr(rg, 'binning') else None
    params['energy'] = rg.energy if hasattr(rg, 'energy') else None
    params['gain_mask_level'] = rg.gain_mask_level if hasattr(rg, 'gain_mask_level') else None
    return params

  def fetch_timestamps(self, indexed=False):
    from xfel.ui.components.run_stats_plotter import \
      get_multirun_should_have_indexed_timestamps, get_paths_from_timestamps, get_strings_from_timestamps
    runs, timestamps = \
      get_multirun_should_have_indexed_timestamps(self.stats,
                                                  self.run_numbers,
                                                  self.parent.run_window.runstats_tab.d_min,
                                                  self.parent.run_window.runstats_tab.n_strong,
                                                  indexed=indexed)

    image_paths_by_run = []
    timestamps_and_params_by_run = []
    for i in range(len(runs)):
      run = runs[i]
      ts = timestamps[i]
      outdir = "out" if indexed else "all"
      prepend = os.path.join(get_run_path(self.output, *self.trgr[run]), outdir)
      tag = "idx" if indexed else "shot"
      image_paths = get_paths_from_timestamps(ts, prepend=prepend, tag=tag, ext=self.trgr[run][1].format)
      image_paths_by_run.append(image_paths)
      timestamp_strings = get_strings_from_timestamps(ts, long_form=True)
      timestamps_and_params_by_run.append((self.get_xtc_process_params_for_run(*self.trgr[run]), timestamp_strings))
    if indexed:
      self.parent.run_window.runstats_tab.strong_indexed_image_timestamps = \
        timestamps_and_params_by_run
      self.parent.run_window.runstats_tab.strong_indexed_image_paths = \
        image_paths_by_run
    else:
      self.parent.run_window.runstats_tab.should_have_indexed_timestamps = \
        timestamps_and_params_by_run
      self.parent.run_window.runstats_tab.should_have_indexed_image_paths = \
        image_paths_by_run

  def plot_stats(self):
    from xfel.ui.components.run_stats_plotter import plot_multirun_stats
    self.refresh_stats()
    sizex, sizey = self.parent.run_window.runstats_tab.runstats_panelsize
    figure = self.parent.run_window.runstats_tab.figure
    figure.clear()
    plot_multirun_stats(
      self.stats, self.run_numbers,
      d_min=self.parent.run_window.runstats_tab.d_min,
      n_multiples=self.parent.run_window.runstats_tab.n_multiples,
      interactive=True,
      ratio_cutoff=self.parent.run_window.runstats_tab.ratio,
      n_strong_cutoff=self.parent.run_window.runstats_tab.n_strong,
      i_sigi_cutoff=self.parent.run_window.runstats_tab.i_sigi,
      run_tags=self.run_tags,
      run_statuses=self.run_statuses,
      minimalist=self.parent.run_window.runstats_tab.entire_expt,
      xsize=(sizex-25)/85, ysize=(sizey-25)/95,
      high_vis=self.parent.high_vis,
      figure=figure)
      # convert px to inches with fudge factor for scaling inside borders
    figure.canvas.draw_idle()

# ----------------------------- Spotfinder Sentinel ---------------------------- #

# Set up events for monitoring spotfinder results against a set threshold
tp_EVT_SPOTFINDER_REFRESH = wx.NewEventType()
EVT_SPOTFINDER_REFRESH = wx.PyEventBinder(tp_EVT_SPOTFINDER_REFRESH, 1)

class RefreshSpotfinder(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class SpotfinderSentinel(Thread):
  ''' Worker thread for spotfinder stats; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.output = self.parent.params.output_folder
    self.number_of_pickles = 0
    self.info = {}
    self.run_numbers = []
    self.stats = []
    self.spot_length_stats = []
    self.run_tags = []
    self.run_statuses = []

    # on initialization (and restart), make sure spotfinder stats drawn from scratch
    self.parent.run_window.spotfinder_tab.redraw_windows = True

  def post_refresh(self):
    evt = RefreshSpotfinder(tp_EVT_SPOTFINDER_REFRESH, -1, self.info)
    wx.PostEvent(self.parent.run_window.spotfinder_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    self.db = xfel_db_application(self.parent.params)

    while self.active:
      self.parent.run_window.spotfinder_light.change_status('idle')
      self.plot_stats_static()
      self.post_refresh()
      self.info = {}
      self.parent.run_window.spotfinder_light.change_status('on')
      time.sleep(5)

  def refresh_stats(self):
    from xfel.ui.db.stats import SpotfinderStats
    from xfel.ui.components.spotfinder_scraper import get_spot_length_stats
    import copy
    if self.parent.run_window.spotfinder_tab.trial_no is not None:
      trial = self.db.get_trial(
        trial_number=self.parent.run_window.spotfinder_tab.trial_no)
      selected_runs = copy.deepcopy(self.parent.run_window.spotfinder_tab.selected_runs)
      self.run_numbers = []
      trial_ids = []
      rungroup_ids = []
      self.stats = []
      self.spot_length_stats = []
      self.trgr = {}
      self.run_tags = []
      self.run_statuses = []
      self.output = self.parent.params.output_folder
      for rg in trial.rungroups:
        for run in rg.runs:
          if run.run not in self.run_numbers and run.run in selected_runs:
            self.run_numbers.append(run.run)
            trial_ids.append(trial.id)
            rungroup_ids.append(rg.id)
            self.trgr[run.run] = (trial, rg, run)
            # spot count
            sf_stats = SpotfinderStats(self.db, run.run, trial.trial, rg.id)()
            self.stats.append(sf_stats)
            self.run_tags.append([tag.name for tag in run.tags])
            # spot lengths
            if self.parent.params.dispatcher == "cxi.xtc_process": #LABELIT backend
              outdir = "integration"
            else:
              outdir = "out"
            run_outdir = os.path.join(get_run_path(self.output, trial, rg, run), outdir)
            try:
              self.spot_length_stats.append(get_spot_length_stats(run_outdir, ref_stats=sf_stats))
            except OSError:
              print("Outdir %s no longer accessible." % run_outdir)
            except Exception as e:
              print(e)
              from dials.array_family import flex
              self.spot_length_stats.append((flex.double(), flex.double(), flex.double()))

      jobs = self.db.get_all_jobs()
      for idx in range(len(self.run_numbers)):
        run_no = self.run_numbers[idx]
        rg_id = rungroup_ids[idx]
        t_id = trial_ids[idx]
        for job in jobs:
          if job.run.run == run_no and job.rungroup.id == rg_id and job.trial.id == t_id:
            self.run_statuses.append(job.status)

  def plot_stats_static(self):
    from xfel.ui.components.spotfinder_plotter import plot_multirun_spotfinder_stats
    self.refresh_stats()
    sizex, sizey = self.parent.run_window.spotfinder_tab.spotfinder_panelsize
    self.parent.run_window.spotfinder_tab.png = plot_multirun_spotfinder_stats(
      self.stats, self.run_numbers,
      spot_length_stats=self.spot_length_stats,
      interactive=False,
      run_tags=self.run_tags,
      run_statuses=self.run_statuses,
      n_min=self.parent.run_window.spotfinder_tab.n_min,
      minimalist=self.parent.run_window.spotfinder_tab.entire_expt,
      easy_run=True,
      xsize=(sizex-25)/85, ysize=sizey/95,
      high_vis=self.parent.high_vis)
      # convert px to inches with fudge factor for scaling inside borders
    self.parent.run_window.spotfinder_tab.redraw_windows = True

# ---------------------------- Image Dumping Thread ---------------------------- #

class ImageDumpThread(Thread):
  def __init__(self,
               command):
    Thread.__init__(self)
    self.active = True
    self.command = command

  def run(self):
    print(self.command)
    easy_run.fully_buffered(command=self.command).show_stderr()

# ----------------------------- Unit Cell Sentinel ----------------------------- #

# Set up events for monitoring unit cell statistics
tp_EVT_UNITCELL_REFRESH = wx.NewEventType()
EVT_UNITCELL_REFRESH = wx.PyEventBinder(tp_EVT_UNITCELL_REFRESH, 1)

class RefreshUnitCell(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class UnitCellSentinel(Thread):
  ''' Worker thread for unit cell analysis; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    evt = RefreshUnitCell(tp_EVT_UNITCELL_REFRESH, -1)
    wx.PostEvent(self.parent.run_window.unitcell_tab, evt)

  def run(self):
    import xfel.ui.components.xfel_gui_plotter as pltr

    feature_vectors = {
      "Triclinic": None,
      "Monoclinic": "a,b,c",
      "Orthorhombic": "a,b,c",
      "Tetragonal": "a,c",
      "Hexagonal": "a,c",
      "Cubic": None,
    }

    # one time post for an initial update
    self.post_refresh()
    self.db = xfel_db_application(self.parent.params)

    while self.active:
      self.parent.run_window.unitcell_light.change_status('idle')
      trial = self.parent.run_window.unitcell_tab.trial
      tag_sets = self.parent.run_window.unitcell_tab.tag_sets
      sizex, sizey = self.parent.run_window.unitcell_tab.unit_cell_panelsize

      info_list = []
      legend_list = []
      for tag_set in tag_sets:
        legend_list.append(str(tag_set))
        cells = self.db.get_stats(trial=trial,
                                  tags=tag_set.tags,
                                  isigi_cutoff=1.0,
                                  tag_selection_mode=tag_set.mode)()
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

      iqr_ratio = 1.5 if self.parent.run_window.unitcell_tab.reject_outliers else None

      figure = self.parent.run_window.unitcell_tab.figure
      plotter = pltr.PopUpCharts(interactive=True, figure=figure)

      if not self.parent.run_window.unitcell_tab.plot_clusters:
        figure.clear()
        plotter.plot_uc_histogram(
          info_list=info_list,
          legend_list=legend_list,
          xsize=(sizex-115)/82, ysize=(sizey-115)/82,
          high_vis=self.parent.high_vis,
          iqr_ratio=iqr_ratio)
        figure.canvas.draw_idle()
      elif len(info_list) > 0:
        from uc_metrics.clustering.step1 import phil_scope
        from uc_metrics.clustering.step_dbscan3d import dbscan_plot_manager
        from cctbx.sgtbx import space_group_info

        if len(info_list) > 1:
          print("Warning, only first tag set will be plotted")

        params = phil_scope.extract()
        try:
          sg = self.parent.run_window.unitcell_tab.trial.cell.lookup_symbol
        except AttributeError:
          sg = "P1"
        sg = "".join(sg.split()) # remove spaces
        params.input.space_group = sg

        iterable = ["{a} {b} {c} {alpha} {beta} {gamma} ".format(**c) + sg for c in info_list[0]]
        params.input.__inject__('iterable', iterable)
        params.file_name = None
        params.cluster.dbscan.eps = float(self.parent.run_window.unitcell_tab.plot_eps.eps.GetValue())
        params.show_plot = True
        params.plot.legend = legend_list[0]
        reject_outliers = self.parent.run_window.unitcell_tab.chk_reject_outliers.GetValue()
        params.plot.outliers = not reject_outliers

        sginfo = space_group_info(params.input.space_group)
        cs = sginfo.group().crystal_system()
        params.input.feature_vector = feature_vectors.get(cs)

        if params.input.feature_vector:
          figure = self.parent.run_window.unitcell_tab.figure
          figure.clear()
          plots = dbscan_plot_manager(params)
          plots.wrap_3D_features(fig = figure, embedded = True)
          figure.canvas.draw_idle()
        else:
          print("Unsupported crystal system", cs)

      self.post_refresh()
      self.parent.run_window.unitcell_light.change_status('on')
      time.sleep(5)

# ------------------------------- Frames Sentinel ------------------------------- #

# Set up events for FramesSeninel
tp_EVT_FRAMES_REFRESH = wx.NewEventType()
EVT_FRAMES_REFRESH = wx.PyEventBinder(tp_EVT_FRAMES_REFRESH, 1)

class RefreshFrames(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class FramesSentinel(Thread):
  ''' Worker thread for frames; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    pass

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    db = xfel_db_application(self.parent.parent.params)

    while self.active:
      trial = db.get_trial(trial_number=int(self.parent.trial_number.ctr.GetStringSelection()))
      runs = [db.get_run(run_number=int(r)) for r in self.parent.trial_runs.ctr.GetCheckedStrings()]
      print("Total events in trial", trial.trial, end=' ')
      if len(runs) == 0:
        runs = None
      else:
        print("runs", ", ".join(sorted([str(r.run) for r in runs])), end=' ')
      print(":", len(db.get_all_events(trial, runs)))
      self.post_refresh()
      time.sleep(2)

# ------------------------------- Clustering --------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_CLUSTERING = wx.NewEventType()
EVT_CLUSTERING = wx.PyEventBinder(tp_EVT_CLUSTERING, 1)

class ClusteringResult(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class Clusterer():
  def __init__(self, trial, runblocks, tags, output, sample_size, threshold):
    self.trial = trial
    self.runblocks = runblocks
    self.tags = tags
    self.output = output
    self.sample_size = sample_size
    self.threshold = threshold

  def unit_cell_clustering(self):

    # 1. Get all pickle files, check if new ones arrived
    run_numbers = []
    rb_paths = []
    tag_ids = set([t.id for t in self.tags])
    for rb in self.runblocks:
      for run in rb.runs:
        if run.run not in run_numbers:
          if len(tag_ids) > 0:
            run_tag_ids = set([t.id for t in run.tags])
            if len(tag_ids.intersection(run_tag_ids)) == 0:
              continue
          run_numbers.append(run.run)
          # test for integration folder
          path = os.path.join(get_run_path(self.output, self.trial, rb, run), "integration")
          if not os.path.exists(path):
            path = os.path.join(get_run_path(self.output, self.trial, rb, run), "out")
          rb_paths.append(path)

    all_pickles = []

    for path in rb_paths:
      try:
        pickles = [os.path.join(path, i) for i in os.listdir(path) if
                   i.endswith('pickle') and 'int-' in i]
        all_pickles = all_pickles + pickles
      except OSError as error:
        print('Folder not found!')
        print(error)

    if len(all_pickles) == 0:
      print('No images integrated (yet)')
      return

    # If clustering button was pressed, do clustering
    # 2. Pick subset
    subset = list(np.random.choice(all_pickles, size=int(self.sample_size)))

    # 3. Do clustering
    ucs = Cluster.from_files(subset, use_b=True)
    clusters, _ = ucs.ab_cluster(self.threshold,
                                 log=False, write_file_lists=False,
                                 schnell=False, doplot=False)
    return clusters

class ClusteringWorker(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               trial,
               runblocks,
               tags,
               output,
               sample_size=1000,
               threshold=250,):
    Thread.__init__(self)
    self.parent = parent
    self.trial = trial
    self.runblocks = runblocks
    self.tags = tags
    self.output = output
    self.sample_size = sample_size
    self.threshold = threshold

  def run(self):
    clusterer = Clusterer(self.trial, self.runblocks, self.tags, self.output,
                          self.sample_size, self.threshold)
    self.clusters = clusterer.unit_cell_clustering()

    evt = ClusteringResult(tp_EVT_CLUSTERING, -1, self.clusters)
    wx.PostEvent(self.parent, evt)

# ----------------------------- Merging Stats Sentinel ---------------------------- #

# Set up events for monitoring merging stats
tp_EVT_MERGINGSTATS_REFRESH = wx.NewEventType()
EVT_MERGINGSTATS_REFRESH = wx.PyEventBinder(tp_EVT_MERGINGSTATS_REFRESH, 1)

class RefreshMergingStats(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class MergingStatsSentinel(Thread):
  ''' Worker thread for merging stats; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.output = self.parent.params.output_folder

    # on initialization (and restart), make sure run stats drawn from scratch
    self.parent.run_window.mergingstats_tab.redraw_windows = True

  def post_refresh(self):
    evt = RefreshMergingStats(tp_EVT_MERGINGSTATS_REFRESH, -1)
    wx.PostEvent(self.parent.run_window.mergingstats_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()
    self.db = xfel_db_application(self.parent.params)

    while self.active:
      self.parent.run_window.mergingstats_light.change_status('idle')
      self.plot_stats_static()
      self.post_refresh()
      self.parent.run_window.mergingstats_light.change_status('on')
      time.sleep(5)

  def plot_stats_static(self):
    from xfel.ui.db.merging_log_scraper import Scraper

    if not self.parent.run_window.mergingstats_tab.dataset_versions: return

    sizex, sizey = self.parent.run_window.mergingstats_tab.mergingstats_panelsize
    sizex = (sizex-25)/85
    sizey = (sizey-25)/95

    if len(self.parent.run_window.mergingstats_tab.dataset_versions) > 1:
      all_results = []
      for folder in self.parent.run_window.mergingstats_tab.dataset_versions:
        scraper = Scraper(folder, '#')
        all_results.append(scraper.scrape())
      self.parent.run_window.mergingstats_tab.png = scraper.plot_many_results(all_results,
                                                                              self.parent.run_window.mergingstats_tab.dataset_name,
                                                                              sizex, sizey, interactive=False)
    else:
      scraper = Scraper(self.parent.run_window.mergingstats_tab.dataset_versions[0], '%')
      results = scraper.scrape()
      self.parent.run_window.mergingstats_tab.png = scraper.plot_single_results(results,
                                                                                self.parent.run_window.mergingstats_tab.dataset_name,
                                                                                sizex, sizey, interactive=False)
    self.parent.run_window.mergingstats_tab.redraw_windows = True

# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(200, 200))

    self.run_sentinel = None
    self.job_sentinel = None
    self.job_monitor = None
    self.spotfinder_sentinel = None
    self.runstats_sentinel = None
    self.unitcell_sentinel = None
    self.mergingstats_sentinel = None

    self.params = load_cached_settings()
    self.db = None

    self.high_vis = False

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddTool(wx.ID_EXIT,
                       label='Quit',
                       bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                       bmpDisabled=wx.NullBitmap,
                       shortHelp='Quit',
                       longHelp='Exit CCTBX.XFEL')
    self.toolbar.AddSeparator()
    if not self.params.monitoring_mode:
      self.tb_btn_watch_new_runs = self.toolbar.AddTool(wx.ID_ANY,
                                   label='Watch for new runs',
                                   bitmap=wx.Bitmap('{}/32x32/quick_restart.png'.format(icons)),
                                   bmpDisabled=wx.NullBitmap,
                                   shortHelp='Watch for new runs',
                                   longHelp='Watch for new runs')
      self.tb_btn_auto_submit = self.toolbar.AddTool(wx.ID_ANY,
                                label='Auto-submit jobs',
                                bitmap=wx.Bitmap('{}/32x32/play.png'.format(icons)),
                                bmpDisabled=wx.NullBitmap,
                                shortHelp='Auto-submit jobs',
                                longHelp='Auto-submit all pending jobs')
      self.toolbar.AddSeparator()
      #self.tb_btn_calibrate = self.toolbar.AddTool(wx.ID_ANY,
      #                    label='Calibration',
      #                    bitmap=wx.Bitmap('{}/32x32/calib.png'.format(icons)),
      #                    bmpDisabled=wx.NullBitmap,
      #                    shortHelp='Calibration',
      #                    longHelp='Detector geometry calibration')
      #self.toolbar.AddSeparator()
      self.tb_btn_settings = self.toolbar.AddTool(wx.ID_ANY,
                          label='Settings',
                          bitmap=wx.Bitmap('{}/32x32/settings.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='Settings',
                          longHelp='Database, user and experiment settings')
    self.tb_btn_zoom = self.toolbar.AddCheckTool(wx.ID_ANY,
                                                 label='Large text',
                                                 bitmap1=wx.Bitmap('{}/32x32/search.png'.format(icons)),
                                                 bmpDisabled=wx.NullBitmap,
                                                 shortHelp='Change text size',
                                                 longHelp='Change text size for plots')
    self.toolbar.Realize()

    # Status bar
    self.sb = self.CreateStatusBar()

    # Menu bar
    menubar = wx.MenuBar()
    m_help = wx.Menu()
    self.mb_about = m_help.Append(wx.ID_ANY, '&About')
    self.mb_docs = m_help.Append(wx.ID_ANY, '&Online help')
    menubar.Append(m_help, '&Help')
    self.SetMenuBar(menubar)

    # Place elements in main PRIME window
    main_box = wx.BoxSizer(wx.VERTICAL)

    # Instantiate windows
    self.run_window = RunWindow(self)

    # Single input window
    main_box.Add(self.run_window, 1, flag=wx.ALL | wx.EXPAND, border=10)
    main_box.Add((-1, 20))


    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)
    self.Bind(wx.EVT_MENU, self.OnDocs, self.mb_docs)

    # Bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    if not self.params.monitoring_mode:
      self.Bind(wx.EVT_TOOL, self.onWatchRuns, self.tb_btn_watch_new_runs)
      self.Bind(wx.EVT_TOOL, self.onAutoSubmit, self.tb_btn_auto_submit)
      #self.Bind(wx.EVT_TOOL, self.onCalibration, self.tb_btn_calibrate)
      self.Bind(wx.EVT_TOOL, self.onSettings, self.tb_btn_settings)
    self.Bind(wx.EVT_TOOL, self.onZoom, self.tb_btn_zoom)
    self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onTabChange,
              self.run_window.main_nbook)
    self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.onLeavingTab,
              self.run_window.main_nbook)

    # Draw the main window sizer
    self.SetSizer(main_box)

  def connect_to_db(self, drop_tables = False):
    self.db = xfel_db_application(self.params, drop_tables = drop_tables, verify_tables = True)

    return True

  def stop_sentinels(self):
    if not self.params.monitoring_mode:
      self.stop_run_sentinel()
      self.stop_job_sentinel()
      self.stop_job_monitor()
    #self.stop_spotfinder_sentinel()
    self.stop_runstats_sentinel()
    self.stop_unitcell_sentinel()
    self.stop_mergingstats_sentinel()

  def start_run_sentinel(self):
    self.run_sentinel = RunSentinel(self, active=True)
    self.run_sentinel.start()
    self.run_window.run_light.change_status('on')

    self.toolbar.SetToolNormalBitmap(self.tb_btn_watch_new_runs.Id, wx.Bitmap('{}/32x32/pause.png'.format(icons)))

  def stop_run_sentinel(self, block = True):
    if self.run_sentinel is not None and self.run_sentinel.active:
      self.run_sentinel.active = False
      if block:
        self.run_sentinel.join()
    self.run_window.run_light.change_status('off')
    self.toolbar.SetToolNormalBitmap(self.tb_btn_watch_new_runs.Id, wx.Bitmap('{}/32x32/quick_restart.png'.format(icons)))

  def start_job_monitor(self):
    self.job_monitor = JobMonitor(self, active=True)
    self.job_monitor.start()

  def stop_job_monitor(self, block = True):
    if self.job_monitor is not None and self.job_monitor.active:
      self.job_monitor.active = False
      if block:
        self.job_monitor.join()

  def start_job_sentinel(self):
    self.job_sentinel = JobSentinel(self, active=True)
    self.job_sentinel.start()
    self.run_window.job_light.change_status('on')

    self.toolbar.SetToolNormalBitmap(self.tb_btn_auto_submit.Id, wx.Bitmap('{}/32x32/pause.png'.format(icons)))

  def stop_job_sentinel(self, block = True):
    if self.job_sentinel is not None and self.job_sentinel.active:
      self.job_sentinel.active = False
      if block:
        self.job_sentinel.join()
    self.run_window.job_light.change_status('off')
    self.toolbar.SetToolNormalBitmap(self.tb_btn_auto_submit.Id, wx.Bitmap('{}/32x32/play.png'.format(icons)))

  def start_prg_sentinel(self):
    self.prg_sentinel = ProgressSentinel(self, active=True)
    self.prg_sentinel.start()
    self.run_window.prg_light.change_status('on')

  def stop_prg_sentinel(self, block = True):
    if self.prg_sentinel is not None and self.prg_sentinel.active:
      self.prg_sentinel.active = False
      if block:
        self.prg_sentinel.join()
    self.run_window.prg_light.change_status('off')

  def start_spotfinder_sentinel(self):
    self.spotfinder_sentinel = SpotfinderSentinel(self, active=True)
    self.spotfinder_sentinel.start()
    self.run_window.spotfinder_light.change_status('on')

  def stop_spotfinder_sentinel(self, block = True):
    if self.spotfinder_sentinel is not None and self.spotfinder_sentinel.active:
      self.spotfinder_sentinel.active = False
      if block:
        self.spotfinder_sentinel.join()
    self.run_window.spotfinder_light.change_status('off')

  def start_runstats_sentinel(self):
    self.runstats_sentinel = RunStatsSentinel(self, active=True)
    self.runstats_sentinel.start()
    self.run_window.runstats_light.change_status('on')

  def stop_runstats_sentinel(self, block = True):
    if self.runstats_sentinel is not None and self.runstats_sentinel.active:
      self.runstats_sentinel.active = False
      if block:
        self.runstats_sentinel.join()
    self.run_window.runstats_light.change_status('off')

  def start_unitcell_sentinel(self):
    self.unitcell_sentinel = UnitCellSentinel(self, active=True)
    self.unitcell_sentinel.start()
    self.run_window.unitcell_light.change_status('on')

  def stop_unitcell_sentinel(self, block = True):
    if self.unitcell_sentinel is not None and self.unitcell_sentinel.active:
      self.unitcell_sentinel.active = False
      if block:
        self.unitcell_sentinel.join()
    self.run_window.unitcell_light.change_status('off')

  def start_mergingstats_sentinel(self):
    self.mergingstats_sentinel = MergingStatsSentinel(self, active=True)
    self.mergingstats_sentinel.start()
    self.run_window.mergingstats_light.change_status('on')

  def stop_mergingstats_sentinel(self, block = True):
    if self.mergingstats_sentinel is not None and self.mergingstats_sentinel.active:
      self.mergingstats_sentinel.active = False
      if block:
        self.mergingstats_sentinel.join()
    self.run_window.mergingstats_light.change_status('off')

  def OnAboutBox(self, e):
    ''' About dialog '''
    import wx.adv
    info = wx.adv.AboutDialogInfo()
    info.SetName('cctbx.xfel')
    info.SetLicense(open(license).read())
    info.SetDescription(description)
    info.AddDeveloper('Artem Lyubimov')
    info.AddDeveloper('Aaron Brewster')
    info.AddDeveloper('Iris Young')
    info.AddDeveloper('Asmit Bhowmick')
    info.AddDeveloper('Daniel Paley')
    info.AddDeveloper('Derek A. Mendez')
    info.AddDeveloper('Johannes Blaschke')
    info.AddDeveloper('Robert Bolotovsky')
    info.AddDeveloper('Axel Brunger')
    info.AddDeveloper('Nicholas Sauter')
    wx.adv.AboutBox(info)

  def OnDocs(self, e):
    import webbrowser
    url = 'http://cci.lbl.gov/publications/download/CCN_2019_p22_Brewster.pdf'
    print('Opening', url)
    webbrowser.open(url)

  def onSettings(self, e):
    settings_dlg = dlg.SettingsDialog(self,
                                      params=self.params)
    settings_dlg.db_cred.btn_big.Disable()
    settings_dlg.SetTitle('Settings')

    if (settings_dlg.ShowModal() == wx.ID_OK):
      self.params = settings_dlg.params
      save_cached_settings(self.params)
      if self.params.facility.name == 'lcls':
        self.title = 'CCTBX.XFEL | {} | {}'.format(self.params.experiment_tag,
                                                   self.params.facility.lcls.experiment)
      else:
        self.title = 'CCTBX.XFEL | {}'.format(self.params.experiment_tag)

  def onZoom(self, e):
    self.high_vis = not self.high_vis

  def onCalibration(self, e):
    calib_dlg = dlg.CalibrationDialog(self, db=self.db)
    calib_dlg.Fit()

    calib_dlg.ShowModal()

  def onWatchRuns(self, e):
    ''' Toggle autosubmit '''
    if self.run_sentinel is not None and self.run_sentinel.active:
      self.stop_run_sentinel(block = True)
    else:
      self.start_run_sentinel()

  def onAutoSubmit(self, e):
    ''' Toggle autosubmit '''
    if self.job_sentinel is not None and self.job_sentinel.active:
      self.stop_job_sentinel(block = True)
    else:
      self.start_job_sentinel()

  def onTabChange(self, e):
    name = self.run_window.main_nbook.GetPageText((self.run_window.main_nbook.GetSelection()))
    if name == self.run_window.jobs_tab.name:
      if self.job_monitor is None or not self.job_monitor.active:
        self.start_job_monitor()
        self.run_window.jmn_light.change_status('on')
    # Disabled
    #elif name == self.run_window.spotfinder_tab.name:
    #  if self.job_monitor is None or not self.job_monitor.active:
    #    self.start_job_monitor()
    #    self.run_window.jmn_light.change_status('on')
    #  if self.spotfinder_sentinel is None or not self.spotfinder_sentinel.active:
    #    self.start_spotfinder_sentinel()
    #    self.run_window.spotfinder_light.change_status('on')
    elif name == self.run_window.trials_tab.name:
      self.run_window.trials_tab.refresh_trials()
    elif name == self.run_window.runstats_tab.name:
      if not self.params.monitoring_mode and (self.job_monitor is None or not self.job_monitor.active):
        self.start_job_monitor()
        self.run_window.jmn_light.change_status('on')
      if self.run_window.runstats_tab.auto_update and (self.runstats_sentinel is None or not self.runstats_sentinel.active):
        self.start_runstats_sentinel()
        self.run_window.runstats_light.change_status('on')
    elif name == self.run_window.unitcell_tab.name:
      if self.run_window.unitcell_tab.auto_update and (self.unitcell_sentinel is None or not self.unitcell_sentinel.active):
        self.start_unitcell_sentinel()
        self.run_window.unitcell_light.change_status('on')
    elif name == self.run_window.datasets_tab.name:
      self.run_window.datasets_tab.refresh_datasets()
    elif name == self.run_window.mergingstats_tab.name:
      self.run_window.mergingstats_tab.refresh_datasets()
      if self.mergingstats_sentinel is None or not self.mergingstats_sentinel.active:
        self.start_mergingstats_sentinel()
        self.run_window.mergingstats_light.change_status('on')
    #elif name == self.run_window.merge_tab.name:
    #  self.run_window.merge_tab.find_trials()

  def onLeavingTab(self, e):
    name = self.run_window.main_nbook.GetPageText((self.run_window.main_nbook.GetSelection()))
    if name == self.run_window.jobs_tab.name:
      if self.job_monitor.active:
        self.stop_job_monitor(block = False)
        self.run_window.jmn_light.change_status('off')
    # Disabled
    #elif name == self.run_window.spotfinder_tab.name:
    #  if self.job_monitor.active:
    #    self.stop_job_monitor(block = False)
    #    self.run_window.jmn_light.change_status('off')
    #  if self.spotfinder_sentinel.active:
    #    self.stop_spotfinder_sentinel(block = False)
    #    self.run_window.spotfinder_light.change_status('off')
    elif name == self.run_window.runstats_tab.name:
      if not self.params.monitoring_mode and self.job_monitor.active:
        self.stop_job_monitor(block = False)
        self.run_window.jmn_light.change_status('off')
      if self.runstats_sentinel.active:
        self.stop_runstats_sentinel(block = False)
        self.run_window.runstats_light.change_status('off')
    elif name == self.run_window.unitcell_tab.name:
      if self.unitcell_sentinel.active:
        self.stop_unitcell_sentinel(block = False)
        self.run_window.unitcell_light.change_status('off')
    elif name == self.run_window.mergingstats_tab.name:
      if self.mergingstats_sentinel.active:
        self.stop_mergingstats_sentinel(block = False)
        self.run_window.mergingstats_light.change_status('off')


  def onQuit(self, e):
    self.stop_sentinels()
    save_cached_settings(self.params)
    self.Destroy()


class RunWindow(wx.Panel):
  ''' Window panel that will house all the run tabs '''

  def __init__(self, parent):
    self.parent = parent
    super(RunWindow, self).__init__(self.parent)

    self.main_panel = wx.Panel(self)
    self.main_nbook = wx.Notebook(self.main_panel, style=0)
    self.runs_tab = RunTab(self.main_nbook, main=self.parent)
    self.trials_tab = TrialsTab(self.main_nbook, main=self.parent)
    self.jobs_tab = JobsTab(self.main_nbook, main=self.parent)
    #self.spotfinder_tab = SpotfinderTab(self.main_nbook, main=self.parent) # Disabled
    self.runstats_tab = RunStatsTab(self.main_nbook, main=self.parent)
    self.unitcell_tab = UnitCellTab(self.main_nbook, main=self.parent)
    self.datasets_tab = DatasetTab(self.main_nbook, main=self.parent)
    #self.merge_tab = MergeTab(self.main_nbook, main=self.parent)
    self.mergingstats_tab = MergingStatsTab(self.main_nbook, main=self.parent)
    self.main_nbook.AddPage(self.runs_tab, self.runs_tab.name)
    self.main_nbook.AddPage(self.trials_tab, self.trials_tab.name)
    self.main_nbook.AddPage(self.jobs_tab, self.jobs_tab.name)
    #self.main_nbook.AddPage(self.spotfinder_tab, self.spotfinder_tab.name) # Disabled
    self.main_nbook.AddPage(self.runstats_tab, self.runstats_tab.name)
    self.main_nbook.AddPage(self.unitcell_tab, self.unitcell_tab.name)
    self.main_nbook.AddPage(self.datasets_tab, self.datasets_tab.name)
    #self.main_nbook.AddPage(self.merge_tab, self.merge_tab.name)
    self.main_nbook.AddPage(self.mergingstats_tab, self.mergingstats_tab.name)

    self.sentinel_box = wx.FlexGridSizer(1, 6, 0, 20)
    self.run_light = gctr.SentinelStatus(self.main_panel, label='Run Sentinel')
    self.job_light = gctr.SentinelStatus(self.main_panel, label='Job Sentinel')
    self.jmn_light = gctr.SentinelStatus(self.main_panel, label='Job Monitor')
    #self.spotfinder_light = gctr.SentinelStatus(self.main_panel, label='Spotfinder Sentinel')
    self.runstats_light = gctr.SentinelStatus(self.main_panel, label='Run Stats Sentinel')
    self.unitcell_light = gctr.SentinelStatus(self.main_panel, label='Unit Cell Sentinel')
    self.mergingstats_light = gctr.SentinelStatus(self.main_panel, label='Merging Stats Sentinel')
    self.sentinel_box.Add(self.run_light)
    self.sentinel_box.Add(self.job_light)
    self.sentinel_box.Add(self.jmn_light)
    #self.sentinel_box.Add(self.spotfinder_light)
    self.sentinel_box.Add(self.runstats_light)
    self.sentinel_box.Add(self.unitcell_light)
    self.sentinel_box.Add(self.mergingstats_light)

    nb_sizer = wx.BoxSizer(wx.VERTICAL)
    nb_sizer.Add(self.main_nbook, 1, flag=wx.EXPAND | wx.ALL, border=3)
    nb_sizer.Add((-1, 20))
    nb_sizer.Add(self.sentinel_box, flag=wx.ALIGN_CENTER_HORIZONTAL)
    self.main_panel.SetSizer(nb_sizer)

    main_sizer = wx.BoxSizer(wx.VERTICAL)
    main_sizer.Add(self.main_panel, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.SetSizer(main_sizer)

    if self.parent.params.monitoring_mode:
      self.runs_tab.Hide()
      self.trials_tab.Hide()
      self.jobs_tab.Hide()
      self.datasets_tab.Hide()
      self.run_light.Hide()
      self.job_light.Hide()
      self.jmn_light.Hide()



# --------------------------------- UI Tabs ---------------------------------- #

class BaseTab(wx.Panel):
  ''' Base class for runtime tab '''

  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(200, 200))

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)


class RunTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Runs'
    self.main = main
    self.last_run = 0
    self.all_runs = []
    self.all_tags = []
    self.all_tag_buttons = []
    self.persistent_tags = []

    self.run_panel = ScrolledPanel(self)
    self.run_sizer = wx.BoxSizer(wx.VERTICAL)
    self.run_panel.SetSizer(self.run_sizer)

    self.colname_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    run_label = wx.StaticText(self, label='Run', size=(60, -1))
    tag_label = wx.StaticText(self, label='Sample Tags', size=(620, -1))
    self.colname_sizer.Add(run_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(tag_label, flag=wx.ALIGN_RIGHT | wx.EXPAND)
    self.colname_sizer.AddGrowableCol(1, 1)
    self.main_sizer.Add(self.colname_sizer, flag=wx.ALL | wx.EXPAND, border=10)

    self.btn_multirun_tags = wx.Button(self, label='Change Tags on Multiple Runs', size=(240, -1))
    self.btn_persistent_tags = gctr.Button(self, name='btn_persistent_tags', label='Manage Persistent Tags', size=(240, -1))
    self.btn_manage_tags = gctr.Button(self, name='btn_manage_tags', label='Manage Tags', size=(120, -1))
    self.main_sizer.Add(self.run_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)

    self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.button_sizer.Add(self.btn_multirun_tags,
                          flag=wx.RIGHT | wx.LEFT | wx.BOTTOM,
                          border=10)
    self.button_sizer.Add(self.btn_persistent_tags,
                          flag=wx.RIGHT | wx.LEFT | wx.BOTTOM,
                          border=10)
    self.button_sizer.Add(self.btn_manage_tags,
                          flag=wx.RIGHT | wx.LEFT | wx.BOTTOM,
                          border=10)
    self.main_sizer.Add(self.button_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Bindings
    self.Bind(EVT_RUN_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_BUTTON, self.onMultiRunTags, self.btn_multirun_tags)
    self.Bind(wx.EVT_BUTTON, self.onManagePersistentTags, self.btn_persistent_tags)
    self.Bind(wx.EVT_BUTTON, self.onManageTags, self.btn_manage_tags)

  def onRefresh(self, e):
      self.refresh_rows()

  def onMultiRunTags(self, e):
    """Add or remove tags applied to multiple runs at a time."""
    multirun_tag_dialog = dlg.MultiRunTagDialog(self, db=self.main.db)
    multirun_tag_dialog.Fit()
    multirun_tag_dialog.ShowModal()
    multirun_tag_dialog.Destroy()

  def onManageTags(self, e):
    ''' User can add / remove / edit sample tags '''
    mtag_dlg = dlg.TagDialog(self, db=self.main.db)
    mtag_dlg.Fit()
    mtag_dlg.ShowModal()
    mtag_dlg.Destroy()

    # Update tags on all tag buttons
    for btn in self.all_tag_buttons:
      btn.tags = btn.run.tags
      btn.update_label()

  def onManagePersistentTags(self, e):
    '''Update which tags are applied automatically to new runs'''
    tags = self.main.db.get_all_tags()
    tag_names = [i.name for i in tags]
    mptag_dlg = wx.MultiChoiceDialog(self,
                                     message='Available tags',
                                     caption='Persistent Tags',
                                     choices=tag_names)
    # Get indices of selected items (if any) and set them to checked
    persistent_tags = self.persistent_tags
    indices = [tag_names.index(i) for i in tag_names if i in persistent_tags]
    mptag_dlg.SetSelections(indices)
    mptag_dlg.Fit()
    if (mptag_dlg.ShowModal() == wx.ID_OK):
      indices = mptag_dlg.GetSelections()
      tag_names = [t.name for t in tags if tag_names.index(t.name) in indices]
      self.persistent_tags = tag_names
    mptag_dlg.Destroy()

  def refresh_rows(self, all=False):
    # Get new runs
    old_run_numbers = [run.run for run in self.all_runs]
    all_runs = self.main.db.get_all_runs()
    new_runs = [run for run in all_runs if run.run not in old_run_numbers]

    font = self.GetFont()
    dc = wx.ScreenDC()
    dc.SetFont(font)
    if len(all_runs) > 0:
      max_width = max([dc.GetTextExtent(str(run))[0] for run in all_runs])
    else:
      max_width = None

    # Update either all or only new runs
    if all:
      runs = self.all_runs
    else:
      runs = new_runs
    for run in runs:
      run_row = RunEntry(self.run_panel, run=run, params=self.main.params, label_width = max_width)
      self.all_tag_buttons.append(run_row.tag_button)
      self.run_sizer.Add(run_row, flag=wx.ALL | wx.EXPAND, border=0)

    self.all_runs = all_runs

    # Update labels on all new tag buttons
    self.all_tags = self.main.db.get_all_tags()
    for button in self.all_tag_buttons:
      button.all_tags = self.all_tags
      button.update_label()

    self.run_panel.SetupScrolling(scrollToTop=False)
    self.run_panel.Refresh()

class TrialsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Trials'

    self.main = main
    self.show_active_only = False

    self.trial_panel = ScrolledPanel(self, size=(200, 200))
    self.trial_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.trial_panel.SetSizer(self.trial_sizer)

    self.btn_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.btn_sizer.AddGrowableCol(0)
    self.btn_add_trial = wx.Button(self, label='New Trial', size=(120, -1))
    self.btn_active_only = wx.ToggleButton(self,
                                           label='Show Only Active Trials',
                                    size=(180, self.btn_add_trial.GetSize()[1]))
    self.btn_sizer.Add(self.btn_active_only, flag=wx.ALIGN_RIGHT)
    self.btn_sizer.Add(self.btn_add_trial)

    self.main_sizer.Add(self.trial_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Bindings
    self.Bind(EVT_RUN_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_BUTTON, self.onAddTrial, self.btn_add_trial)
    self.Bind(wx.EVT_TOGGLEBUTTON, self.onActiveOnly, self.btn_active_only)

  def onRefresh(self, e):
    self.refresh_trials()

  def refresh_trials(self):
    self.trial_sizer.Clear(delete_windows=True)
    self.all_trials = self.main.db.get_all_trials()
    for trial in self.all_trials:
      if self.show_active_only:
        if trial.active:
          self.add_trial(trial=trial)
      else:
        self.add_trial(trial=trial)

    self.trial_panel.SetSizer(self.trial_sizer)
    #self.trial_panel.Layout()
    self.trial_sizer.Layout()
    self.trial_panel.SetupScrolling(scrollToTop=False)

  def add_trial(self, trial):
    new_trial = TrialPanel(self.trial_panel,
                           db=self.main.db,
                           trial=trial,
                           box_label='Trial {} {}'.format(trial.trial,
                             trial.comment[:min(len(trial.comment), 10)] if trial.comment is not None else ""))
    new_trial.chk_active.SetValue(trial.active)
    new_trial.refresh_trial()
    self.trial_sizer.Add(new_trial, flag=wx.EXPAND | wx.ALL, border=10)

  def onAddTrial(self, e):
    new_trial_dlg = dlg.TrialDialog(self, db=self.main.db)
    new_trial_dlg.Fit()

    if new_trial_dlg.ShowModal() == wx.ID_OK:
      self.refresh_trials()

  def onActiveOnly(self, e):
    self.show_active_only = self.btn_active_only.GetValue()
    self.refresh_trials()


class JobsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Jobs'

    self.main = main
    self.all_trials = []
    self.all_jobs = None
    self.filter = 'All jobs'
    self.data = {}

    self.job_list = gctr.SortableListCtrl(self, style=wx.LC_REPORT|wx.BORDER_SUNKEN)
    self.job_list.InsertColumn(0, "Job")
    self.job_list.InsertColumn(1, "Type")
    self.job_list.InsertColumn(2, "Dataset")
    self.job_list.InsertColumn(3, "Trial")
    self.job_list.InsertColumn(4, "Run")
    self.job_list.InsertColumn(5, "Block")
    self.job_list.InsertColumn(6, "Task")
    self.job_list.InsertColumn(7, "Subm ID")
    self.job_list.InsertColumn(8, "Status")

    self.job_list_sort_flag = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    self.job_list_col = 0

    self.trial_choice = gctr.ChoiceCtrl(self,
                                        label='Filter by:',
                                        label_size=(60, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])
    self.btn_stop_job = wx.Button(self, label='Stop job', size=(120, -1))
    self.btn_delete_job = wx.Button(self, label='Delete job', size=(120, -1))
    self.btn_restart_job = wx.Button(self, label='Restart job', size=(120, -1))
    self.chk_active = wx.CheckBox(self, label='Only display jobs from active trials/blocks')
    self.chk_active.SetValue(True)
    self.option_sizer = wx.FlexGridSizer(1, 5, 0, 20)
    self.option_sizer.AddMany([(self.trial_choice), (self.chk_active), (self.btn_stop_job), (self.btn_delete_job), (self.btn_restart_job)])

    self.main_sizer.Add(self.job_list, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.option_sizer, flag=wx.EXPAND | wx.ALL, border=10)


    self.Bind(wx.EVT_BUTTON, self.onStopJob, self.btn_stop_job)
    self.Bind(wx.EVT_BUTTON, self.onDeleteJob, self.btn_delete_job)
    self.Bind(wx.EVT_BUTTON, self.onRestartJob, self.btn_restart_job)
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_choice.ctr)
    self.chk_active.Bind(wx.EVT_CHECKBOX, self.onToggleActive)
    self.Bind(EVT_JOB_MONITOR, self.onMonitorJobs)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.onColClick, self.job_list)

  def onToggleActive(self, e):
    self.main.job_monitor.only_active_jobs = self.chk_active.GetValue()

  def onTrialChoice(self, e):
    self.filter = self.trial_choice.ctr.GetString(
                  self.trial_choice.ctr.GetSelection())

  def GetSelectedJobIds(self):
    return [self.job_list.GetItemData(i) for i in range(self.job_list.GetItemCount()) if self.job_list.IsSelected(i)]

  def onStopJob(self, e):
    if self.all_jobs is None:
      return

    jobs_to_stop = self.GetSelectedJobIds()
    if len(jobs_to_stop) == 0:
      return

    if len(jobs_to_stop) == 1:
      message='Are you sure to stop job %d?'%jobs_to_stop[0]
    else:
      message='Are you sure to stop %d jobs?'%len(jobs_to_stop)

    msg = wx.MessageDialog(self,
                           message=message,
                           caption='Warning',
                           style=wx.YES_NO | wx.ICON_EXCLAMATION)

    if (msg.ShowModal() == wx.ID_NO):
      return

    from xfel.ui.components.submission_tracker import JobStopper
    stopper = JobStopper(self.main.params.mp.method)
    for job in self.all_jobs:
      if job.id in jobs_to_stop:
        stopper.stop_job(job.submission_id)

  def onDeleteJob(self, e):
    if self.all_jobs is None:
      return

    jobs_to_delete = self.GetSelectedJobIds()
    if len(jobs_to_delete) == 0:
      return

    if len(jobs_to_delete) == 1:
      message='Are you sure to delete all processing results from job %d?'%jobs_to_delete[0]
    else:
      message='Are you sure to delete all processing results from %d jobs?'%len(jobs_to_delete)

    msg = wx.MessageDialog(self,
                           message=message,
                           caption='Warning',
                           style=wx.YES_NO | wx.ICON_EXCLAMATION)

    if (msg.ShowModal() == wx.ID_NO):
      return

    for job in self.all_jobs:
      if job.id in jobs_to_delete:
        job.delete()

  def onRestartJob(self, e):
    if self.all_jobs is None:
      return

    jobs_to_restart = self.GetSelectedJobIds()
    if len(jobs_to_restart) == 0:
      return

    if len(jobs_to_restart) == 1:
      message='Are you sure to restart job %d? This will delete all processing results from from this job and re-submit it. Be sure the job has been stopped first.'%jobs_to_restart[0]
    else:
      message='Are you sure to restart %d jobs? This will delete all processing results from these jobs and re-submit them. Be sure the jobs have been stopped first.'%len(jobs_to_restart)

    msg = wx.MessageDialog(self,
                           message=message,
                           caption='Warning',
                           style=wx.YES_NO | wx.ICON_EXCLAMATION)

    if (msg.ShowModal() == wx.ID_NO):
      return

    for job in self.all_jobs:
      if job.id in jobs_to_restart:
        job.delete()
        if job.status != "DELETED":
          print("Couldn't restart job", job.id, "job is not deleted")
          continue
        job.remove_from_db()

  def onMonitorJobs(self, e):
    # Find new trials
    if e.trials is not None:
      all_db_trials = [str(i.trial) for i in e.trials]
      new_trials = [i for i in all_db_trials if i not in self.all_trials]
      if len(new_trials) > 0:
        self.find_trials()
        self.all_trials = all_db_trials

    if e.jobs is not None:
      self.all_jobs = e.jobs
      if str(self.filter).lower() == 'all jobs':
        selected_trials = [int(t) for t in self.all_trials]
      else:
        selected_trials = [int(self.filter.split()[-1])]

      selected_jobs = self.GetSelectedJobIds()
      if self.job_list.GetFocusedItem() > 0:
        focused_job_id = self.job_list.GetItemData(self.job_list.GetFocusedItem())
      else:
        focused_job_id = None

      self.data = {} # reset contents of the table, with unique row ids
      for job in e.jobs:
        if job.trial is not None:
          if job.trial.trial not in selected_trials: continue

        # Order: job, type, dataset, trial, run, rungroup, submission id, status
        j = str(job.id)
        jt = job.task.type if job.task is not None else "-"
        ds = job.dataset.name if job.dataset is not None else "-"
        t = "t%03d" % job.trial.trial if job.trial is not None else "-"
        try:
          r = "r%04d" % int(job.run.run) if job.run is not None else "-"
        except ValueError:
          r = "r%s" % job.run.run
        rg = "rg%03d" % job.rungroup.id if job.rungroup is not None else "-"
        tsk = "%d" % job.task.id if job.task is not None else "-"
        sid = str(job.submission_id)

        short_status = str(job.status).strip("'")
        if short_status == "SUBMIT_FAIL":
          short_status = "S_FAIL"
        elif short_status == "SUBMITTED":
          short_status = "SUBMIT"
        s = short_status

        self.data[job.id] = [j, jt, ds, t, r, rg, tsk, sid, s]
        found_it = False
        # Look to see if item already in list
        for i in range(self.job_list.GetItemCount()):
          if self.job_list.GetItemData(i) == job.id:
            for k, item in enumerate(self.data[job.id]):
              self.job_list.SetItem(i, k, item)
            found_it = True
            break
        if found_it: continue

        # Item not present, so deposit the row in the table
        local_job_id = self.job_list.Append(self.data[job.id])
        self.job_list.SetItemData(local_job_id, job.id)

      # Remove items not sent in the event or otherwise filtered out
      # Need to do it in reverse order to avoid list re-ordering when deleting items
      for i in reversed(range(self.job_list.GetItemCount())):
        if self.job_list.GetItemData(i) not in self.data:
          self.job_list.DeleteItem(i)

      # Initialize sortable column mixin
      self.job_list.initialize_sortable_columns(n_col=5, itemDataMap=self.data)
      self.job_list.RestoreSortOrder(self.job_list_col, self.job_list_sort_flag)

      # Restore selected items
      for i in range(self.job_list.GetItemCount()):
        job_id = self.job_list.GetItemData(i)
        if job_id in selected_jobs:
          self.job_list.Select(i)
        if job_id == focused_job_id:
          self.job_list.Focus(i)

  def onColClick(self, e):
    # Note: the sortable list binds the event first and activates this method after.
    # print "Click recorded in column %s" % str(self.job_list._col) ## DEBUG
    # print self.job_list.GetSortState() ## DEBUG
    self.job_list_col = self.job_list._col
    self.job_list_sort_flag = self.job_list._colSortFlag

  def find_trials(self):
    print("Found trials")
    if self.main.db is not None:
      choices = ['All jobs'] + \
                ['trial {}'.format(i.trial) for i in self.main.db.get_all_trials()]
      self.trial_choice.ctr.Clear()
      for choice in choices:
        self.trial_choice.ctr.Append(choice)

      if self.filter == 'All jobs':
        self.trial_choice.ctr.SetSelection(0)
      else:
        self.trial_choice.ctr.SetSelection(int(self.filter[-1]))

class SpotfinderTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Spotfinder'

    self.main = main
    self.all_trials = []
    self.trial_no = None
    self.trial = None
    self.all_runs = []
    self.selected_runs = []
    self.tag_trial_changed = True
    self.tag_runs_changed = True
    self.tag_last_five = False
    self.entire_expt = False
    self.png = None
    self.static_bitmap = None
    self.redraw_windows = True
    self.n_min = 4

    self.spotfinder_panel = wx.Panel(self, size=(100, 100))
    self.spotfinder_panelsize = self.spotfinder_panel.GetSize()
    self.spotfinder_box = wx.StaticBox(self.spotfinder_panel, label='Run Statistics')
    self.spotfinder_sizer = wx.StaticBoxSizer(self.spotfinder_box, wx.HORIZONTAL)
    self.spotfinder_panel.SetSizer(self.spotfinder_sizer)

    self.trial_number = gctr.ChoiceCtrl(self,
                                        label='Trial:',
                                        label_size=(90, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])
    self.last_five_runs =  wx.Button(self,
                                     label='Auto plot last five runs',
                                     size=(200, -1))
    self.plot_entire_expt = wx.Button(self,
                                     label='Auto plot entire experiment',
                                     size=(200,-1))
    self.n_min_selector = gctr.OptionCtrl(self,
                                          label='minimum # spots:',
                                          label_size=(160, -1),
                                          ctrl_size=(30, -1),
                                          items=[('n_min', 4)])
    self.run_numbers =  gctr.CheckListCtrl(self,
                                           label='Selected runs:',
                                           label_size=(200, -1),
                                           label_style='normal',
                                           ctrl_size=(150, 224),
                                           direction='vertical',
                                           choices=[])

    self.bottom_sizer = wx.FlexGridSizer(1, 2, 0, 4)

    options_box = wx.StaticBox(self, label='Display Options')
    self.options_box_sizer = wx.StaticBoxSizer(options_box, wx.VERTICAL)
    self.options_opt_sizer = wx.GridBagSizer(1, 1)

    self.options_opt_sizer.Add(self.trial_number, pos=(0, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.last_five_runs, pos=(1, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.plot_entire_expt, pos=(2, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.n_min_selector, pos=(3, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.run_numbers, pos=(0, 1), span=(8, 1),
                               flag=wx.BOTTOM | wx.TOP | wx.RIGHT | wx.EXPAND,
                               border=10)
    self.options_box_sizer.Add(self.options_opt_sizer)
    self.bottom_sizer.Add(self.options_box_sizer)

    self.main_sizer.Add(self.spotfinder_panel, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.bottom_sizer, 0,
                        flag=wx.EXPAND | wx.ALL, border=10)

    # Bindings
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_number.ctr)
    self.Bind(wx.EVT_BUTTON, self.onLastFiveRuns, self.last_five_runs)
    self.Bind(wx.EVT_BUTTON, self.onEntireExpt, self.plot_entire_expt)
    self.Bind(wx.EVT_TEXT_ENTER, self.onNMin, self.n_min_selector.n_min)
    self.Bind(wx.EVT_CHECKLISTBOX, self.onRunChoice, self.run_numbers.ctr)
    self.Bind(EVT_SPOTFINDER_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_SIZE, self.OnSize)

  def OnSize(self, e):
    self.spotfinder_panelsize = self.spotfinder_panel.GetSize()
    e.Skip()

  def onTrialChoice(self, e):
    trial_idx = self.trial_number.ctr.GetSelection()
    if trial_idx == 0:
      self.trial_no = None
      self.trial = None
      self.run_numbers.ctr.Clear()
      self.all_runs = []
      self.selected_runs = []
    else:
      trial_no = self.trial_number.ctr.GetClientData(trial_idx)
      if trial_no is not None:
        self.trial_no = int(trial_no)
        self.trial = self.main.db.get_trial(trial_number=int(self.trial_no))
        self.spotfinder_box.SetLabel('Spotfinder Results - Trial {}'.format(self.trial_no))
        self.find_runs()

  def onRunChoice(self, e):
    self.tag_last_five = False
    self.entire_expt = False
    run_numbers_selected = self.run_numbers.ctr.GetCheckedStrings()
    if self.trial is not None:
      self.selected_runs = [r.run for r in self.trial.runs if r.run in run_numbers_selected]
      self.main.run_window.spotfinder_light.change_status('idle')

  def find_trials(self):
    all_db_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    new_trials = [i for i in all_db_trials if i not in self.all_trials]
    if len(new_trials) > 0:
      self.trial_number.ctr.Clear()
      self.all_trials = [None] + all_db_trials
      for trial in self.all_trials:
        if trial is None:
          entry = 'None'
          self.trial_number.ctr.Append(entry)
          self.trial_number.ctr.SetClientData(0, None)
        else:
          entry = trial
          self.trial_number.ctr.Append(entry)
          item_idx = self.trial_number.ctr.FindString(entry)
          self.trial_number.ctr.SetClientData(item_idx, trial)

      if self.trial_no is not None:
        self.trial_number.ctr.SetSelection(self.trial_no)
      else:
        self.trial_number.ctr.SetSelection(0)

  def find_runs(self):
    self.run_numbers.ctr.Clear()
    if self.trial is not None:
      self.runs_available = [str(r.run) for r in self.trial.runs]
      if len(self.all_runs) > 0:
        self.run_numbers.ctr.InsertItems(items=self.all_runs, pos=0)

  def onRefresh(self, e):
    self.refresh_trials()
    self.refresh_runs()
    if self.tag_last_five:
      self.select_last_n_runs(5)
    elif self.entire_expt:
      self.select_all()
    if self.redraw_windows:
      self.plot_static_spotfinder_stats()
      self.redraw_windows = False
    if self.trial is not None:
      self.spotfinder_box.SetLabel('Spotfinder Results - Trial {}'.format(self.trial_no))
    else:
      self.spotfinder_box.SetLabel('Spotfinder Results - No trial selected')

  def refresh_trials(self):
    if self.all_trials == []:
      self.find_trials()
    avail_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    for t in avail_trials:
      if t not in self.all_trials:
        self.trial_number.ctr.Append(t)
        self.all_trials.append(t)
        item_idx = self.trial_number.ctr.FindString(t)
        self.trial_number.ctr.SetClientData(item_idx, t)

  def refresh_runs(self):
    if self.all_runs == []:
      self.find_runs()
    if self.trial is not None:
      avail_runs = [str(r.run) for r in self.trial.runs]
      for r in avail_runs:
        if r not in self.all_runs:
          self.run_numbers.ctr.Append(r)
          self.all_runs.append(r)

  def plot_static_spotfinder_stats(self):
    if self.png is not None:
      if self.static_bitmap is not None:
        self.static_bitmap.Destroy()
      img = wx.Image(self.png, wx.BITMAP_TYPE_ANY)
      self.static_bitmap = wx.StaticBitmap(
        self.spotfinder_panel, wx.ID_ANY, wx.Bitmap(img))
      self.spotfinder_sizer.Add(self.static_bitmap, 0, wx.EXPAND | wx.ALL, 3)
      self.spotfinder_panel.SetSizer(self.spotfinder_sizer)
      self.spotfinder_panel.Layout()

  def select_last_n_runs(self, n):
    if self.trial is not None:
      self.selected_runs = [r.run for r in self.trial.runs][-n:]

  def select_all(self):
    if self.trial is not None:
      self.selected_runs = [r.run for r in self.trial.runs]

  def onLastFiveRuns(self, e):
    self.entire_expt = False
    self.tag_last_five = True
    self.select_last_n_runs(5)
    self.main.run_window.spotfinder_light.change_status('idle')

  def onEntireExpt(self, e):
    self.entire_expt = True
    self.tag_last_five = False
    self.select_all()
    self.main.run_window.spotfinder_light.change_status('idle')

  def onNMin(self, e):
    n_min = self.n_min_selector.n_min.GetValue()
    if n_min.isdigit():
      self.n_min = int(n_min)

class RunStatsTab(SpotfinderTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Run Stats'

    self.main = main
    self.all_trials = []
    self.trial_no = None
    self.trial = None
    self.all_runs = []
    self.selected_runs = []
    self.tag_trial_changed = True
    self.tag_runs_changed = True
    self.tag_last_five = False
    self.entire_expt = False
    self.d_min = 2.5
    self.n_multiples = 2
    self.ratio = 1
    self.n_strong = 16
    self.i_sigi = 1
    self.n_dump = 10
    self.should_have_indexed_image_paths = None
    self.should_have_indexed_timestamps = None
    self.strong_indexed_image_paths = None
    self.strong_indexed_image_timestamps = None
    self.auto_update = True

    self.runstats_panel = wx.Panel(self, size=(100, 100))
    self.runstats_box = wx.StaticBox(self.runstats_panel, label='Run Statistics')
    self.runstats_sizer = wx.StaticBoxSizer(self.runstats_box, wx.HORIZONTAL)
    self.runstats_panel.SetSizer(self.runstats_sizer)

    import matplotlib as mpl
    from matplotlib.backends.backend_wxagg import (
      FigureCanvasWxAgg as FigureCanvas,
      NavigationToolbar2WxAgg as NavigationToolbar)

    self.figure = mpl.figure.Figure()
    self.canvas = FigureCanvas(self.runstats_box, -1, self.figure)
    self.toolbar = NavigationToolbar(self.canvas)
    self.toolbar.SetWindowStyle(wx.TB_VERTICAL)
    self.toolbar.Realize()

    self.runstats_sizer.Add(self.canvas, 1, wx.EXPAND)
    self.runstats_sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

    self.options_box = wx.StaticBox(self, label='Statistics Options')

    self.trial_number = gctr.ChoiceCtrl(self.options_box,
                                        label='Trial:',
                                        label_size=(90, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])
    self.last_five_runs =  wx.Button(self.options_box,
                                     label='Auto plot last five runs',
                                     size=(200, -1))
    self.plot_entire_expt = wx.Button(self.options_box,
                                     label='Auto plot entire experiment',
                                     size=(200,-1))
    self.d_min_select = gctr.OptionCtrl(self.options_box,
                                        name='rs_d_min',
                                        label='High resolution limit:',
                                        sub_labels=[''],
                                        label_size=(160, -1),
                                        ctrl_size=(30, -1),
                                        items=[('d_min', 2.5)])
    self.n_multiples_selector = gctr.OptionCtrl(self.options_box,
                                               name='rs_multiples',
                                               label='# multiples threshold:',
                                               sub_labels=[''],
                                               label_size=(160, -1),
                                               ctrl_size=(30, -1),
                                               items=[('multiples', 2)])
    self.ratio_cutoff = gctr.OptionCtrl(self.options_box,
                                        name='rs_ratio',
                                        label='two theta ratio cutoff:',
                                        sub_labels=[''],
                                        label_size=(160, -1),
                                        ctrl_size=(30, -1),
                                        items=[('ratio', 1)])
    self.n_strong_cutoff = gctr.OptionCtrl(self.options_box,
                                           name='rs_n_strong',
                                           label='# strong spots cutoff:',
                                           sub_labels=[''],
                                           label_size=(160, -1),
                                           ctrl_size=(30, -1),
                                           items=[('n_strong', 16)])
    self.i_sigi_cutoff = gctr.OptionCtrl(self.options_box,
                                         name='rs_isigi',
                                         label='I/sig(I) cutoff:',
                                         sub_labels=[''],
                                         label_size=(160, -1),
                                         ctrl_size=(30, -1),
                                         items=[('isigi', 1)])
    self.n_dump_cutoff = gctr.OptionCtrl(self.options_box,
                                         name='rs_n_dump',
                                         label='# images to dump:',
                                         sub_labels=[''],
                                         label_size=(160, -1),
                                         ctrl_size=(30, -1),
                                         items=[('n_dump', 10)])
    self.run_numbers =  gctr.CheckListCtrl(self.options_box,
                                           label='Selected runs:',
                                           label_size=(200, -1),
                                           label_style='normal',
                                           ctrl_size=(150, 224),
                                           direction='vertical',
                                           choices=[])

    self.strong_indexed_box = wx.StaticBox(self, label='Strongest Indexed Images')
    self.strong_indexed_list = wx.TextCtrl(self.strong_indexed_box,
                                           style=wx.TE_MULTILINE | wx.TE_READONLY)
    self.idx_show_images_button = wx.Button(self.strong_indexed_box,
                                            label='Open images',
                                            size=(200, -1))
    self.should_have_indexed_box = wx.StaticBox(self, label='Strong Images that Didn\'t Index')
    self.should_have_indexed_list = wx.TextCtrl(self.should_have_indexed_box,
                                                style=wx.TE_MULTILINE | wx.TE_READONLY)
    self.shi_dump_images_button = wx.Button(self.should_have_indexed_box,
                                            label='Dump images',
                                            size=(200, -1))

    self.bottom_sizer = wx.FlexGridSizer(1, 2, 0, 10)

    self.options_box_sizer = wx.StaticBoxSizer(self.options_box, wx.VERTICAL)
    self.options_opt_sizer = wx.GridBagSizer(1, 1)

    self.options_opt_sizer.Add(self.trial_number, pos=(0, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.last_five_runs, pos=(1, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.plot_entire_expt, pos=(2, 0),
                                   flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.d_min_select, pos=(3, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.n_multiples_selector, pos=(4, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.ratio_cutoff, pos=(5, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.n_strong_cutoff, pos=(7, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.i_sigi_cutoff, pos=(6, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.n_dump_cutoff, pos=(8, 0),
                               flag=wx.ALL, border=2)
    self.options_opt_sizer.Add(self.run_numbers, pos=(0, 1), span=(8, 1),
                               flag=wx.BOTTOM | wx.TOP | wx.RIGHT | wx.EXPAND,
                               border=10)
    self.options_box_sizer.Add(self.options_opt_sizer)
    self.bottom_sizer.Add(self.options_box_sizer)

    self.dump_images_sizer = wx.GridBagSizer(3, 1)

    self.strong_indexed_box_sizer = wx.StaticBoxSizer(self.strong_indexed_box, wx.VERTICAL)

    self.strong_indexed_results_sizer = wx.GridBagSizer(1, 1)
    self.strong_indexed_results_sizer.Add(self.strong_indexed_list, pos=(0, 0),
                                          span=(5, 45),
                                          flag=wx.LEFT | wx.RIGHT | wx.EXPAND,
                                          border=10)
    self.strong_indexed_box_sizer.Add(self.strong_indexed_results_sizer)

    self.strong_indexed_box_sizer.Add(self.idx_show_images_button,
                                      flag=wx.LEFT | wx.RIGHT | wx.ALL | wx.EXPAND,
                                      border=5)

    self.should_have_indexed_box_sizer = wx.StaticBoxSizer(self.should_have_indexed_box, wx.VERTICAL)

    self.should_have_indexed_results_sizer = wx.GridBagSizer(1, 1)
    self.should_have_indexed_results_sizer.Add(self.should_have_indexed_list, pos=(0, 0),
                                               span=(5, 45),
                                               flag=wx.LEFT | wx.RIGHT | wx.EXPAND,
                                               border=10)
    self.should_have_indexed_box_sizer.Add(self.should_have_indexed_results_sizer)

    self.should_have_indexed_box_sizer.Add(self.shi_dump_images_button,
                                           flag=wx.LEFT | wx.RIGHT | wx.ALL | wx.EXPAND,
                                           border=5)

    self.manage_panel = wx.Panel(self)
    self.manage_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.btn_toggle_options = wx.ToggleButton(self.manage_panel,
                                              label='Hide options')
    self.chk_auto_update = wx.CheckBox(self.manage_panel, label='Auto update')
    self.chk_auto_update.SetValue(True)
    self.manage_sizer.Add(self.btn_toggle_options)
    self.manage_sizer.Add(self.chk_auto_update)
    self.manage_panel.SetSizer(self.manage_sizer)

    self.dump_images_sizer.Add(self.manage_panel, pos=(0, 0))
    self.dump_images_sizer.Add(self.strong_indexed_box_sizer, pos=(1, 0))
    self.dump_images_sizer.Add(self.should_have_indexed_box_sizer, pos=(2, 0))

    if self.main.params.dispatcher != "cctbx.xfel.xtc_process":
      self.n_dump_cutoff.Hide()
      self.strong_indexed_box.Hide()
      self.strong_indexed_list.Hide()
      self.idx_show_images_button.Hide()
      self.should_have_indexed_box.Hide()
      self.should_have_indexed_list.Hide()
      self.shi_dump_images_button.Hide()

    # self.bottom_sizer.Add(self.should_have_indexed_box_sizer, flag=wx.EXPAND | wx.ALL)
    # self.bottom_sizer.Add(self.strong_indexed_box_sizer, flag=wx.EXPAND | wx.ALL)
    self.bottom_sizer.Add(self.dump_images_sizer, flag=wx.EXPAND | wx.ALL)

    self.main_sizer.Add(self.runstats_panel, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.bottom_sizer, 0,
                        flag=wx.EXPAND | wx.ALL, border=10)

    self.static_bitmap = wx.StaticBitmap(
      self.runstats_panel, wx.ID_ANY)#, wx.Bitmap(img))
    self.runstats_sizer.Add(self.static_bitmap, 0, wx.EXPAND | wx.ALL, 3)
    self.runstats_panel.SetSizer(self.runstats_sizer)

    # Bindings
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_number.ctr)
    self.Bind(wx.EVT_BUTTON, self.onLastFiveRuns, self.last_five_runs)
    self.Bind(wx.EVT_BUTTON, self.onEntireExpt, self.plot_entire_expt)
    self.Bind(wx.EVT_TEXT_ENTER, self.onDMin, self.d_min_select.d_min)
    self.Bind(wx.EVT_TEXT_ENTER, self.onMultiples, self.n_multiples_selector.multiples)
    self.Bind(wx.EVT_TEXT_ENTER, self.onRatioCutoff, self.ratio_cutoff.ratio)
    self.Bind(wx.EVT_TEXT_ENTER, self.onHitCutoff, self.n_strong_cutoff.n_strong)
    self.Bind(wx.EVT_TEXT_ENTER, self.onIsigICutoff, self.i_sigi_cutoff.isigi)
    self.Bind(wx.EVT_TEXT_ENTER, self.onNDump, self.n_dump_cutoff.n_dump)
    self.Bind(wx.EVT_CHECKLISTBOX, self.onRunChoice, self.run_numbers.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOpenImages, self.idx_show_images_button)
    self.Bind(wx.EVT_BUTTON, self.onDumpImages, self.shi_dump_images_button)
    self.Bind(wx.EVT_TOGGLEBUTTON, self.onToggleOptions, self.btn_toggle_options)
    self.Bind(wx.EVT_CHECKBOX, self.onChkAutoUpdate, self.chk_auto_update)
    self.Bind(EVT_RUNSTATS_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_SIZE, self.OnSize)

    self.Layout()
    self.Fit()
    self.runstats_panelsize = self.runstats_box.GetSize()

  def OnSize(self, e):
    self.runstats_panelsize = self.runstats_box.GetSize()
    e.Skip()

  def onTrialChoice(self, e):
    trial_idx = self.trial_number.ctr.GetSelection()
    self.trial_no = None
    self.trial = None
    self.run_numbers.ctr.Clear()
    self.all_runs = []
    self.selected_runs = []
    if trial_idx > 0:
      trial_no = self.trial_number.ctr.GetClientData(trial_idx)
      if trial_no is not None:
        self.trial_no = int(trial_no)
        self.trial = self.main.db.get_trial(trial_number=int(self.trial_no))
        self.runstats_box.SetLabel('Run Statistics - Trial {}'.format(self.trial_no))
        self.find_runs()

  def onRunChoice(self, e):
    self.tag_last_five = False
    self.entire_expt = False
    run_numbers_selected = self.run_numbers.ctr.GetCheckedStrings()
    if self.trial is not None:
      self.selected_runs = [r.run for r in self.trial.runs if r.run in run_numbers_selected]
      self.main.run_window.runstats_light.change_status('idle')

  def onRefresh(self, e):
    self.refresh_trials()
    self.refresh_runs()
    if self.tag_last_five:
      self.select_last_n_runs(5)
    elif self.entire_expt:
      self.select_all()
    self.print_strong_indexed_paths()
    self.print_should_have_indexed_paths()
    if self.trial is not None:
      self.runstats_box.SetLabel('Run Statistics - Trial {}'.format(self.trial_no))
    else:
      self.runstats_box.SetLabel('Run Statistics - No trial selected')

  def onToggleOptions(self, e):
    if self.btn_toggle_options.GetValue():
      self.options_box.Hide()
      self.strong_indexed_box.Hide()
      self.should_have_indexed_box.Hide()
    else:
      self.options_box.Show()
      if self.main.params.dispatcher == "cctbx.xfel.xtc_process":
        self.strong_indexed_box.Show()
        self.should_have_indexed_box.Show()
    self.Layout()
    self.Fit()

  def onChkAutoUpdate(self, e):
    self.auto_update = self.chk_auto_update.GetValue()

    if self.auto_update and (self.main.runstats_sentinel is None or not self.main.runstats_sentinel.active):
      self.main.start_runstats_sentinel()
    else:
      self.main.stop_runstats_sentinel()

  def print_strong_indexed_paths(self):
    try:
      paths = []
      for p in self.strong_indexed_image_paths:
        paths.extend(p)
      image_paths = '\n'.join(paths)
      self.strong_indexed_list.SetValue(image_paths)
    except TypeError:
      print("Error getting list of best indexed images")
      pass

  def print_should_have_indexed_paths(self):
    if self.trial is not None:
      try:
        paths = []
        for p in self.should_have_indexed_image_paths:
          paths.extend(p)
        image_paths = '\n'.join(paths)
        self.should_have_indexed_list.SetValue(image_paths)
      except TypeError:
        print("Error getting list of images that should have indexed")
        pass

  def onLastFiveRuns(self, e):
    self.entire_expt = False
    self.tag_last_five = True
    self.select_last_n_runs(5)
    self.main.run_window.runstats_light.change_status('idle')

  def onEntireExpt(self, e):
    self.entire_expt = True
    self.tag_last_five = False
    self.select_all()
    self.main.run_window.runstats_light.change_status('idle')

  def onDMin(self, e):
    try:
      d_min = float(self.d_min_select.d_min.GetValue())
      self.d_min = d_min
    except ValueError:
      pass

  def onMultiples(self, e):
    try:
      mult = int(self.n_multiples_selector.multiples.GetValue())
      self.n_multiples = mult
    except ValueError:
      pass

  def onRatioCutoff(self, e):
    try:
      ratio = float(self.ratio_cutoff.ratio.GetValue())
      self.ratio = ratio
    except ValueError:
      pass

  def onHitCutoff(self, e):
    n_strong = self.n_strong_cutoff.n_strong.GetValue()
    if n_strong.isdigit():
      self.n_strong = int(n_strong)

  def onIsigICutoff(self, e):
    try:
      isigi = float(self.i_sigi_cutoff.isigi.GetValue())
      self.i_sigi = isigi
    except ValueError:
      pass

  def onNDump(self, e):
    n_dump = self.n_dump_cutoff.n_dump.GetValue()
    if n_dump.isdigit():
      self.n_dump = int(n_dump)

  def dump_timestamps(self, params, ts_list, img_list):
    if not os.path.isdir(params['output_dir']):
      os.makedirs(params['output_dir'])
    command = ('cctbx.xfel.xtc_dump input.experiment=%s '%params['experiment'])+\
    ('input.run_num=%s input.address=%s '%(str(params['run']), params['address']))+\
    ('format.file_format=%s '%params['format'])+\
    ('output.output_dir=%s '%params['output_dir'])
    if params['format'] == 'cbf':
      command += 'format.cbf.detz_offset=%f '%params['distance']
      if params['energy'] is not None:
        command += 'format.cbf.override_energy=%f '%params['energy']
      if 'Rayonix' in params['address']:
        command += 'format.cbf.mode=rayonix '
        if params['beamx'] is not None:
          command += 'format.cbf.rayonix.override_beam_x=%d '%params['beamx']
        if params['beamy'] is not None:
          command += 'format.cbf.rayonix.override_beam_y=%d '%params['beamy']
        if params['bin_size'] is not None:
          command += 'format.cbf.rayonix.bin_size=%d '%params['bin_size']
      elif 'cspad' in params['address'].lower():
        if params['gain_mask_level'] is not None:
          command += 'format.cbf.cspad.gain_mask_value=%f '% params['gain_mask_level']
    elif params['format'] == 'pickle':
      if params['config'] is not None:
        command += 'input.cfg=%s '%params['config']
    command += 'dispatch.selected_events=True '
    for timestamp_string in ts_list[:self.n_dump]:
      command += 'input.timestamp=%s '%timestamp_string
    command += '&& dials.image_viewer %s'%\
      ' '.join(map(str, img_list[:self.n_dump]))
    thread = ImageDumpThread(command)
    thread.start()

  def onDumpImages(self, e):
    for idx in range(len(self.should_have_indexed_timestamps)):
      params, ts_list = self.should_have_indexed_timestamps[idx]
      imgs = self.should_have_indexed_image_paths[idx]
      self.dump_timestamps(params, ts_list, imgs)

  def onOpenImages(self, e):
    for idx in range(len(self.strong_indexed_image_timestamps)):
      params, ts_list = self.strong_indexed_image_timestamps[idx]
      ext = '.' + params['format']
      image_paths = self.strong_indexed_image_paths[idx][:self.n_dump]
      indexed_paths = [path.split(ext)[0]+'_indexed.refl' for path in image_paths]
      if all([os.path.exists(p) for p in (image_paths + indexed_paths)]):
        command = str('dials.image_viewer ' + ' '.join(image_paths) + \
          ' ' + ' '.join(indexed_paths))
        thread = ImageDumpThread(command)
        thread.start()
      else:
        shot_paths = [p.split('out')[0] + 'all' + p.split('out')[1].replace('idx', 'shot') \
          for p in image_paths]
        self.dump_timestamps(params, ts_list, shot_paths)

class UnitCellTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Unit Cell'

    self.main = main
    self.all_trials = []
    self.trial_no = None
    self.trial = None
    self.tags = []
    self.tag_sets = []
    self.reject_outliers = True
    self.plot_clusters = False
    self.auto_update = True

    # self.tab_panel = wx.Panel(self, size=(300, 300))
    self.tab_sizer = wx.BoxSizer(wx.HORIZONTAL)
    # self.tab_panel.SetSizer(self.tab_sizer)

    self.selection_columns_panel = wx.Panel(self, size=(230, 120))
    self.selection_columns_box = wx.StaticBox(self.selection_columns_panel, label='Select tag sets')
    self.selection_columns_sizer = wx.StaticBoxSizer(self.selection_columns_box, wx.VERTICAL)
    self.selection_columns_panel.SetSizer(self.selection_columns_sizer)

    # self.selection_columns_panel = wx.Panel(self, size=(100, 120))
    # self.selection_columns_sizer = wx.BoxSizer(wx.HORIZONTAL)
    # self.selection_columns_panel.SetSizer(self.selection_columns_sizer)

    self.trial_number = gctr.ChoiceCtrl(self,
                                        label='Trial:',
                                        label_size=(90, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])

    self.tag_checklist = gctr.CheckListCtrl(self,
                                            label='Available tags:',
                                            label_size=(200, -1),
                                            label_style='normal',
                                            ctrl_size=(150, 100),
                                            direction='vertical',
                                            choices=[])

    self.selection_type_radio = gctr.RadioCtrl(self,
                                               name='uc_selection_type',
                                               label='',
                                               label_style='normal',
                                               label_size=(-1, -1),
                                               direction='horizontal',
                                               items={'union':'union',
                                                      'inter':'intersection'})

    self.add_sele_button = wx.Button(self.selection_columns_panel,
                                     label='Add selection',
                                     size=(200, -1))

    self.tag_set_checklist = gctr.CheckListCtrl(self,
                                                label='Tag sets to display:',
                                                label_size=(200, -1),
                                                label_style='normal',
                                                ctrl_size=(150, 100),
                                                direction='vertical',
                                                choices=[])

    self.remove_sele_button = wx.Button(self.selection_columns_panel,
                                        label='Remove selection',
                                        size=(200, -1))

    self.reset_sele_button = wx.Button(self.selection_columns_panel,
                                       label='Reset selections',
                                       size=(200, -1))

    self.chk_reject_outliers = wx.CheckBox(self.selection_columns_panel, label='Reject outliers')
    self.chk_reject_outliers.SetValue(True)

    self.chk_plot_clusters = wx.CheckBox(self.selection_columns_panel, label='Plot clusters')

    self.chk_auto_update = wx.CheckBox(self.selection_columns_panel, label='Auto update')
    self.chk_auto_update.SetValue(True)

    self.plot_eps = gctr.OptionCtrl(self.selection_columns_panel,
                                    name='uc_plot_eps',
                                    label='Cluster epsilon',
                                    sub_labels=[''],
                                    label_size=(160, -1),
                                    ctrl_size=(50, -1),
                                    items=[('eps', 0.8)])
    self.plot_eps.eps.Disable()

    try:
      import uc_metrics # import dependency
    except ImportError:
      self.chk_plot_clusters.Hide()
      self.plot_eps.Hide()

    self.add_sele_sizer = wx.GridBagSizer(4, 1)
    self.add_sele_sizer.Add(self.trial_number, pos=(0, 0),
                            flag=wx.ALL, border=0)
    self.add_sele_sizer.Add(self.tag_checklist, pos=(1, 0),
                            flag=wx.ALL | wx.EXPAND, border=0)
    self.add_sele_sizer.Add(self.selection_type_radio, pos=(2, 0),
                            flag=wx.ALL | wx.ALIGN_CENTER, border=0)
    self.add_sele_sizer.Add(self.add_sele_button, pos=(3, 0),
                            flag=wx.ALL)
    self.selection_columns_sizer.Add(self.add_sele_sizer, flag=wx.ALL | wx.EXPAND, border=10)

    self.remove_sele_sizer = wx.GridBagSizer(3, 1)
    self.remove_sele_sizer.Add(self.tag_set_checklist, pos=(0, 0),
                               flag=wx.ALL | wx.EXPAND, border=0)
    self.remove_sele_sizer.Add(self.remove_sele_button, pos=(1, 0),
                               flag=wx.ALL)
    self.remove_sele_sizer.Add(self.reset_sele_button, pos=(2, 0),
                               flag=wx.ALL)
    self.selection_columns_sizer.Add(self.remove_sele_sizer, flag=wx.ALL | wx.EXPAND, border=10)
    self.selection_columns_sizer.Add(self.chk_reject_outliers, flag=wx.ALL | wx.EXPAND, border=10)
    self.selection_columns_sizer.Add(self.chk_plot_clusters, flag=wx.ALL | wx.EXPAND, border=10)
    self.selection_columns_sizer.Add(self.plot_eps, flag=wx.ALL | wx.EXPAND, border=10)
    self.selection_columns_sizer.Add(self.chk_auto_update, flag=wx.ALL | wx.EXPAND, border=10)

    self.unit_cell_panel = wx.Panel(self, size=(200, 120))
    self.unit_cell_box = wx.StaticBox(self.unit_cell_panel, label='Unit cell analysis')
    self.unit_cell_panelsize = self.unit_cell_box.GetSize()
    self.unit_cell_sizer = wx.StaticBoxSizer(self.unit_cell_box, wx.VERTICAL)
    self.unit_cell_panel.SetSizer(self.unit_cell_sizer)

    import matplotlib as mpl
    from matplotlib.backends.backend_wxagg import (
      FigureCanvasWxAgg as FigureCanvas,
      NavigationToolbar2WxAgg as NavigationToolbar)

    self.figure = mpl.figure.Figure()
    self.canvas = FigureCanvas(self.unit_cell_box, -1, self.figure)
    self.toolbar = NavigationToolbar(self.canvas)
    self.toolbar.Realize()

    self.unit_cell_sizer.Add(self.canvas, 1, wx.EXPAND)
    self.unit_cell_sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

    # self.main_sizer.Add(self.selection_columns_panel, 1,
    #                     flag=wx.EXPAND | wx.ALL, border=10)
    # self.main_sizer.Add(self.unit_cell_panel, 1,
    #                     flag=wx.EXPAND | wx.ALL, border=10)

    self.tab_sizer.Add(self.selection_columns_panel, 0,
                       flag=wx.ALIGN_LEFT | wx.EXPAND, border=10)
    self.tab_sizer.Add(self.unit_cell_panel, 1,
                       flag=wx.EXPAND | wx.ALL, border=0)
    self.main_sizer.Add(self.tab_sizer, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.selection_columns_sizer.Layout()

    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_number.ctr)
    self.Bind(wx.EVT_BUTTON, self.onAddTagSet, self.add_sele_button)
    self.Bind(wx.EVT_BUTTON, self.onRemoveTagSet, self.remove_sele_button)
    self.Bind(wx.EVT_BUTTON, self.onResetTagSets, self.reset_sele_button)
    self.Bind(wx.EVT_CHECKBOX, self.onChkRejectOutliers, self.chk_reject_outliers)
    self.Bind(wx.EVT_CHECKBOX, self.onChkPlotClusters, self.chk_plot_clusters)
    self.Bind(wx.EVT_CHECKBOX, self.onChkAutoUpdate, self.chk_auto_update)
    self.Bind(EVT_UNITCELL_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_SIZE, self.OnSize)

  def OnSize(self, e):
    self.unit_cell_panelsize = self.unit_cell_box.GetSize()
    e.Skip()

  def find_trials(self):
    all_db_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    new_trials = [i for i in all_db_trials if i not in self.all_trials]
    if len(new_trials) > 0:
      self.trial_number.ctr.Clear()
      self.all_trials = [None] + all_db_trials
      for trial in self.all_trials:
        if trial is not None:
          entry = trial
          self.trial_number.ctr.Append(entry)
          item_idx = self.trial_number.ctr.FindString(entry)
          self.trial_number.ctr.SetClientData(item_idx, trial)
        else:
          entry = 'None'
          self.trial_number.ctr.Append(entry)
          self.trial_number.ctr.SetClientData(0, None)
      if self.trial_no is not None:
        self.trial_number.ctr.SetSelection(self.trial_no)
      else:
        self.trial_number.ctr.SetSelection(0)

  def onTrialChoice(self, e):
    trial_idx = self.trial_number.ctr.GetSelection()
    if trial_idx == 0:
      self.trial_no = None
      self.trial = None
    elif self.trial_number.ctr.GetClientData(trial_idx) != self.trial_no:
      self.trial_no = int(self.trial_number.ctr.GetClientData(trial_idx))
      self.trial = self.main.db.get_trial(trial_number=self.trial_no)
    self.find_tags()

  def find_tags(self):
    self.tag_checklist.ctr.Clear()
    if self.trial is not None:
      self.tags = self.trial.tags
      tag_names = [t.name for t in self.tags]
      if tag_names:
        self.tag_checklist.ctr.InsertItems(items=tag_names, pos=0)
    self.refresh_tag_sets()

  def onAddTagSet(self, e):
    checked_items = self.tag_checklist.ctr.GetCheckedStrings()
    selected_tags = [i for i in self.main.db.get_all_tags()
                     if i.name in checked_items]
    selected_tags = selected_tags
    if self.selection_type_radio.union.GetValue() == 1:
      mode = 'union'
    else:
      mode = 'intersection'
    tag_set = TagSet(mode, selected_tags)
    self.tag_sets.append(tag_set)
    self.refresh_tag_sets()

  def refresh_tag_sets(self):
    self.tag_set_checklist.ctr.Clear()
    tag_set_strings = [str(ts) for ts in self.tag_sets]
    if tag_set_strings:
      self.tag_set_checklist.ctr.InsertItems(items = tag_set_strings, pos=0)

  def onRemoveTagSet(self, e):
    all_items = self.tag_set_checklist.ctr.GetStrings()
    checked_items = self.tag_set_checklist.ctr.GetCheckedStrings()
    selected_tag_sets = [ts for ts in self.tag_sets if str(ts) in checked_items]
    for ts in selected_tag_sets:
      idx = all_items.index(str(ts))
      self.tag_set_checklist.ctr.Delete(idx)
      self.tag_sets.remove(ts)

  def onResetTagSets(self, e):
    self.tag_set_checklist.ctr.Clear()
    self.tag_sets = []
    self.selected_tag_sets = []

  def onChkRejectOutliers(self, e):
    self.reject_outliers = self.chk_reject_outliers.GetValue()

  def onChkPlotClusters(self, e):
    self.plot_clusters = self.chk_plot_clusters.GetValue()
    if self.plot_clusters:
      self.plot_eps.eps.Enable()
    else:
      self.plot_eps.eps.Disable()

  def onChkAutoUpdate(self, e):
    self.auto_update = self.chk_auto_update.GetValue()

    if self.auto_update and (self.main.unitcell_sentinel is None or not self.main.unitcell_sentinel.active):
      self.main.start_unitcell_sentinel()
    else:
      self.main.stop_unitcell_sentinel()

  def onRefresh(self, e):
    self.find_trials()

class DatasetTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Datasets'

    self.main = main

    self.show_active_only = False

    self.dataset_panel = ScrolledPanel(self, size=(200, 200))
    self.dataset_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.dataset_panel.SetSizer(self.dataset_sizer)

    self.btn_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.btn_sizer.AddGrowableCol(0)
    self.btn_add_dataset = wx.Button(self, label='New Dataset', size=(120, -1))
    self.btn_active_only = wx.ToggleButton(self,
                                           label='Show Only Active Datasets',
                                    size=(180, self.btn_add_dataset.GetSize()[1]))
    self.btn_sizer.Add(self.btn_active_only, flag=wx.ALIGN_RIGHT)
    self.btn_sizer.Add(self.btn_add_dataset)

    self.main_sizer.Add(self.dataset_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddDataset, self.btn_add_dataset)
    self.Bind(wx.EVT_TOGGLEBUTTON, self.onActiveOnly, self.btn_active_only)

    #self.refresh_datasets()

  def refresh_datasets(self):
    self.dataset_sizer.Clear(delete_windows=True)
    self.all_datasets = self.main.db.get_all_datasets()
    for dataset in self.all_datasets:
      if self.show_active_only:
        if dataset.active:
          self.add_dataset(dataset=dataset)
      else:
        self.add_dataset(dataset=dataset)

    self.dataset_panel.SetSizer(self.dataset_sizer)
    self.dataset_sizer.Layout()
    self.dataset_panel.SetupScrolling(scrollToTop=False)

  def add_dataset(self, dataset):
    new_dataset = DatasetPanel(self.dataset_panel,
                           db=self.main.db,
                           dataset=dataset)
    new_dataset.chk_active.SetValue(dataset.active)
    new_dataset.refresh_dataset()
    self.dataset_sizer.Add(new_dataset, flag=wx.EXPAND | wx.ALL, border=10)

  def onAddDataset(self, e):
    new_dataset_dlg = dlg.DatasetDialog(self, db=self.main.db)
    new_dataset_dlg.Fit()

    if new_dataset_dlg.ShowModal() == wx.ID_OK:
      self.refresh_datasets()

  def onActiveOnly(self, e):
    self.show_active_only = self.btn_active_only.GetValue()
    self.refresh_datasets()

class MergingStatsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Merging stats'

    self.main = main
    self.all_datasets = []
    self.dataset_versions = []
    self.png = None
    self.static_bitmap = None
    self.redraw_windows = True

    self.tab_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.datasets_panel = wx.Panel(self, size=(240, 120))
    self.datasets_box = wx.StaticBox(self.datasets_panel, label='Select dataset')
    self.datasets_sizer = wx.StaticBoxSizer(self.datasets_box, wx.VERTICAL)
    self.datasets_panel.SetSizer(self.datasets_sizer)

    self.datasets = wx.ListBox(self.datasets_panel,
                               size=(220, 100))
    self.datasets_sizer.Add(self.datasets, flag=wx.EXPAND | wx.ALL, border = 5)

    self.chk_active = wx.CheckBox(self.datasets_panel, label='Active only')
    self.chk_active.SetValue(True)
    self.datasets_sizer.Add(self.chk_active, flag=wx.EXPAND | wx.ALL, border = 5)

    self.dataset_version = gctr.ChoiceCtrl(self.datasets_panel,
                                           label='Dataset version:',
                                           label_size=(120, -1),
                                           label_style='normal',
                                           ctrl_size=(100, -1),
                                           choices=[])
    self.datasets_sizer.Add(self.dataset_version, flag=wx.EXPAND | wx.ALL, border = 5)

    self.plots_panel = wx.Panel(self, size=(200, 120))
    self.mergingstats_panelsize = self.plots_panel.GetSize()
    self.plots_box = wx.StaticBox(self.plots_panel, label='Statistics')
    self.plots_sizer = wx.StaticBoxSizer(self.plots_box, wx.VERTICAL)
    self.plots_panel.SetSizer(self.plots_sizer)

    self.tab_sizer.Add(self.datasets_panel, 0,
                       flag=wx.ALIGN_LEFT | wx.EXPAND, border=10)
    self.tab_sizer.Add(self.plots_panel, 1,
                       flag=wx.EXPAND | wx.ALL, border=0)
    self.main_sizer.Add(self.tab_sizer, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)

    self.Bind(wx.EVT_LISTBOX, self.onSelectDataset, self.datasets)
    self.Bind(wx.EVT_CHOICE, self.onVersionChoice, self.dataset_version.ctr)
    self.Bind(EVT_MERGINGSTATS_REFRESH, self.onRefresh)
    self.chk_active.Bind(wx.EVT_CHECKBOX, self.onToggleActivity)
    self.Bind(wx.EVT_SIZE, self.OnSize)

  def OnSize(self, e):
    self.mergingstats_panelsize = self.plots_panel.GetSize()
    e.Skip()

  def onToggleActivity(self, e):
    self.refresh_datasets()

  def refresh_datasets(self):
    self.datasets.Clear()
    self.all_datasets = self.main.db.get_all_datasets()
    if self.chk_active.GetValue():
      self.all_datasets = [d for d in self.all_datasets if d.active]
    for dataset in self.all_datasets:
      self.datasets.Append(dataset.name)
    self.refresh_dataset()

  def onVersionChoice(self, e):
    self.refresh_stats()

  def onSelectDataset(self, e):
    self.refresh_dataset()

  def refresh_dataset(self):
    self.dataset_version.ctr.Clear()
    sel = self.datasets.GetSelection()
    if sel < 0: return
    try:
      dataset = self.all_datasets[sel]
    except IndexError:
      pass
    else:
      self.dataset_version.ctr.Append('All')
      for version in dataset.versions:
        self.dataset_version.ctr.Append(str(version.version))
      self.dataset_version.ctr.SetSelection(0)
      self.refresh_stats()

  def refresh_stats(self):
    sel = self.datasets.GetSelection()
    dataset = self.all_datasets[sel]
    self.dataset_name = dataset.name
    if self.dataset_version.ctr.GetSelection() == 0:
      self.dataset_versions = [version.output_path() for version in dataset.versions]
    else:
      version = dataset.versions[self.dataset_version.ctr.GetSelection()-1]
      self.dataset_name += " v%03d"%version.version
      self.dataset_versions = [version.output_path()]

  def onRefresh(self, e):
    self.plot_merging_stats()

  def plot_merging_stats(self):
    if self.png is not None:
      if self.static_bitmap is not None:
        try:
          self.static_bitmap.Destroy()
        except RuntimeError as e:
          if "StaticBitmap has been deleted" not in str(e):
            raise
      img = wx.Image(self.png, wx.BITMAP_TYPE_ANY)
      self.static_bitmap = wx.StaticBitmap(
        self.plots_panel, wx.ID_ANY, wx.Bitmap(img))
      self.plots_sizer.Add(self.static_bitmap, 0, wx.EXPAND | wx.ALL, 3)
      self.plots_panel.SetSizer(self.plots_sizer)
      self.plots_panel.Layout()

class MergeTab(BaseTab):

  def __init__(self, parent, main, prefix='prime'):
    BaseTab.__init__(self, parent=parent)
    self.name = 'Merge'

    self.main = main
    self.prefix = prefix
    self.prime_filename = '{}.phil'.format(self.prefix)
    self.output = self.main.params.output_folder
    self.run_paths = []
    self.trial_no = None
    self.all_trials = []
    self.all_tags = []
    self.selected_tags = []
    self.run_paths = []

    self.prime_panel = PRIMEInputWindow(self)
    self.toolbar = wx.ToolBar(self, style=wx.TB_HORZ_TEXT | wx.TB_FLAT)
    self.tb_btn_def = self.toolbar.AddTool(wx.ID_ANY, label=' Defaults',
                          bitmap=wx.Bitmap('{}/24x24/def.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='Default Settings',
                          longHelp='Generate default PRIME settings')
    self.tb_btn_load = self.toolbar.AddTool(wx.ID_OPEN, label=' Load PHIL',
                          bitmap=wx.Bitmap('{}/24x24/open.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='Load PHIL file',
                          longHelp='Load PHIL file with PRIME settings')
    self.tb_btn_save = self.toolbar.AddTool(wx.ID_SAVE, label=' Save PHIL',
                          bitmap=wx.Bitmap('{}/24x24/save.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='Save PHIL file',
                          longHelp='Save PHIL file with PRIME settings')
    self.tb_btn_cmd = self.toolbar.AddTool(wx.ID_ANY, label=' Command',
                          bitmap=wx.Bitmap('{}/24x24/term.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='PRIME Command',
                          longHelp='Output PRIME command to stdout')
    self.toolbar.EnableTool(self.tb_btn_cmd.GetId(), False)
    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddTool(wx.ID_ANY, label=' Run PRIME',
                          bitmap=wx.Bitmap('{}/24x24/run.png'.format(icons)),
                          bmpDisabled=wx.NullBitmap,
                          shortHelp='Run PRIME',
                          longHelp='Scale, merge and post-refine with PRIME')
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.Realize()

    # Modify PRIME input window to hide input control
    self.prime_panel.inp_box.Hide()
    self.prime_panel.out_box.ctr.SetValue(self.output)

    # Input box
    self.input_panel = wx.Panel(self)
    input_box = wx.StaticBox(self.input_panel, label='PRIME Input')
    self.input_box_sizer = wx.StaticBoxSizer(input_box, wx.HORIZONTAL)
    self.input_panel.SetSizer(self.input_box_sizer)

    self.trial_number = gctr.ChoiceCtrl(self.input_panel,
                                        label='Trial:',
                                        label_size=(80, -1),
                                        label_style='normal',
                                        ctrl_size=(140, -1),
                                        choices=[])
    self.tag_title = wx.StaticText(self.input_panel, label='Tags:')
    self.tag_list = gctr.CheckListCtrl(self.input_panel,
                                       ctrl_size=(200, 100),
                                       choices=[],
                                       direction='vertical')
    self.opt_prefix = gctr.OptionCtrl(self.input_panel,
                                      label='List prefix:',
                                      label_size=(80, -1),
                                      ctrl_size=(140, -1),
                                      items=[('prefix', 'prime')])
    self.input_number = wx.StaticText(self.input_panel,
                                      label='0 images in 0 folders:')
    self.input_list = wx.TextCtrl(self.input_panel,
                                  style=wx.TE_MULTILINE | wx.TE_READONLY)

    self.trial_tag_sizer = wx.GridBagSizer(2, 3)
    self.trial_tag_sizer.Add(self.opt_prefix, pos=(0, 0))
    self.trial_tag_sizer.Add(self.trial_number, pos=(1, 0),
                             flag=wx.TOP, border=10)
    self.trial_tag_sizer.Add(self.tag_title, pos=(0, 1),
                             flag=wx.LEFT | wx.EXPAND,
                             border=15)
    self.trial_tag_sizer.Add(self.tag_list, pos=(1, 1),
                             flag=wx.LEFT | wx.EXPAND,
                             border=15)
    self.trial_tag_sizer.Add(self.input_number, pos=(0, 2),
                             flag=wx.LEFT | wx.EXPAND | wx.ALIGN_RIGHT,
                             border=15)
    self.trial_tag_sizer.Add(self.input_list, pos=(1, 2),
                             flag=wx.LEFT | wx.EXPAND | wx.ALIGN_RIGHT,
                             border=15)
    self.input_box_sizer.Add(self.trial_tag_sizer, 1, flag=wx.ALL | wx.EXPAND,
                             border=10)
    self.trial_tag_sizer.AddGrowableCol(2)
    self.trial_tag_sizer.AddGrowableRow(1)
    self.main_sizer.Add(self.toolbar, border=10,
                        flag=wx.EXPAND | wx.LEFT | wx.RIGHT)
    self.main_sizer.Add(self.input_panel, proportion=1,
                        flag=wx.ALL | wx.EXPAND, border=10)
    self.main_sizer.Add(self.prime_panel, border=10,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM | wx.EXPAND)


    self.Bind(wx.EVT_TEXT, self.onInput, self.input_list)
    #self.Bind(wx.EVT_BUTTON, self.onIsoRef, self.prime_panel.ref_box.btn_browse)
    #self.Bind(wx.EVT_TEXT, self.onIsoRef, self.prime_panel.ref_box.ctr)
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_number.ctr)
    self.Bind(wx.EVT_CHECKLISTBOX, self.onTagCheck, self.tag_list.ctr)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_cmd)
    self.Bind(wx.EVT_TOOL, self.onLoad, self.tb_btn_load)
    self.Bind(wx.EVT_TOOL, self.onSave, self.tb_btn_save)

  def onTagCheck(self, e):
    checked_items = self.tag_list.ctr.GetCheckedStrings()
    self.selected_tags = [i for i in self.main.db.get_all_tags() if i.name
                          in checked_items]
    self.find_integrated_pickles()

  def onTrialChoice(self, e):
    trial_idx = self.trial_number.ctr.GetSelection()
    if self.trial_number.ctr.GetClientData(trial_idx) == 0:
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
      self.toolbar.EnableTool(self.tb_btn_cmd.GetId(), False)
      self.tag_list.ctr.Clear()
      self.input_list.SetValue('')
    elif self.trial_number.ctr.GetClientData(trial_idx) != self.trial_no:
      self.trial_no = self.trial_number.ctr.GetClientData(trial_idx)
      self.trial = self.main.db.get_trial(trial_number=int(self.trial_no))
      self.find_tags()
      self.find_integrated_pickles()

  def find_tags(self):
    self.tag_list.ctr.Clear()
    self.tags = []
    self.tag_names = []
    tag_ids = []
    for run in self.trial.runs:
      for tag in run.tags:
        if tag.id not in tag_ids:
          self.tags.append(tag)
          tag_ids.append(tag.id)
          self.tag_names.append(tag.name)
    self.tag_title.SetLabel('Tags for trial {}:'.format(self.trial.trial))
    if self.tag_names:
      self.tag_list.ctr.InsertItems(items=self.tag_names, pos=0)

  def find_trials(self):
    all_db_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    new_trials = [i for i in all_db_trials if i not in self.all_trials]
    if len(new_trials) > 0:
      self.trial_number.ctr.Clear()
      self.all_trials = [None] + \
                        [str(i.trial) for i in self.main.db.get_all_trials()]
      for trial in self.all_trials:
        if trial is not None:
          entry = 'Trial {}'.format(trial)
          self.trial_number.ctr.Append(entry)
          item_idx = self.trial_number.ctr.FindString(entry)
          self.trial_number.ctr.SetClientData(item_idx, trial)
        else:
          entry = '-- select a trial --'
          self.trial_number.ctr.Append(entry)
          self.trial_number.ctr.SetClientData(0, None)

      if self.trial_no is not None:
        self.trial_number.ctr.SetSelection(self.trial_no)
      else:
        self.trial_number.ctr.SetSelection(0)

  def find_integrated_pickles(self):

    # Find runblock paths associated with the trial
    run_numbers = []
    run_ids = []
    self.run_paths = []
    if self.main.params.dispatcher == "cxi.xtc_process": #LABELIT backend
      integration_dir = "integration"
    else:
      integration_dir = "out"
    for rb in self.trial.rungroups:
      for run in rb.runs:
        if run.run not in run_numbers:
          if len(self.selected_tags) == 0:
            self.run_paths.append(os.path.join(
              get_run_path(self.output, self.trial, rb, run), integration_dir))
            run_numbers.append(run.run)
          else:
            for tag_id in [int(t.id) for t in self.selected_tags]:
              if tag_id in [int(t.id) for t in run.tags]:
                run_ids.append(int(run.id))
                self.run_paths.append(os.path.join(
                  get_run_path(self.output, self.trial, rb, run), integration_dir))
                break

    # Display paths in input list text control
    input_paths = '\n'.join(self.run_paths)
    self.input_list.SetValue(input_paths)

    # Find appropriate integration pickles in runblock paths
    self.all_pickles = []
    for path in self.run_paths:
      try:
        pickles = [os.path.join(path, i) for i in os.listdir(path) if
                   i.endswith('pickle') and 'int-' in i]
        self.all_pickles = self.all_pickles + pickles
      except OSError as error:
        print('Folder not found: {}'.format(path))
        continue

    self.input_number.SetLabel('{} images in {} folders:'
                               ''.format(len(self.all_pickles),
                                         len(self.run_paths)))

  def onInput(self, e):
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_cmd.GetId(), True)

  def onLoad(self, e):
    # Extract params from file
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      script = load_dlg.GetPaths()[0]
      out_dir = os.path.dirname(script)
      self.prime_filename = os.path.basename(script)
      self.load_script(out_dir=out_dir)
    load_dlg.Destroy()

  def load_script(self, out_dir):
    ''' Loads PRIME script '''
    import iotbx.phil as ip

    script = os.path.join(out_dir, self.prime_filename)
    user_phil = ip.parse(open(script).read())
    self.pparams = master_phil.fetch(sources=[user_phil]).extract()
    self.prime_panel.pparams = self.pparams

    if len(self.pparams.data) > 0:
      self.prime_panel.inp_box.ctr.SetValue(str(self.pparams.data[0]))
    current_dir = os.path.dirname(self.pparams.run_no)
    self.prime_panel.out_box.ctr.SetValue(str(current_dir))
    if str(self.prime_panel.out_box.ctr.GetValue).lower() == '':
      self.prime_panel.out_box.ctr.SetValue(self.out_dir)
    if str(self.pparams.title).lower() != 'none':
      self.prime_panel.title_box.ctr.SetValue(str(self.pparams.title))
    if str(self.pparams.hklisoin).lower() != 'none':
      self.prime_panel.ref_box.ctr.SetValue(str(self.pparams.hklisoin))
    elif str(self.pparams.hklrefin).lower() != 'none':
      self.prime_panel.ref_box.ctr.SetValue(str(self.pparams.hklrefin))
      self.prime_panel.opt_chk_useref.SetValue(True)
    if str(self.pparams.n_residues).lower() == 'none':
      self.prime_panel.opt_spc_nres.SetValue(500)
    else:
      self.prime_panel.opt_spc_nres.SetValue(int(self.pparams.n_residues))
    self.prime_panel.opt_spc_nproc.SetValue(int(self.pparams.n_processors))

  def onSave(self, e):
    self.init_settings()

    # Generate text of params
    final_phil = master_phil.format(python_object=self.pparams)
    with Capturing() as txt_output:
      final_phil.show()
    txt_out = ''
    for one_output in txt_output:
      txt_out += one_output + '\n'

    # Save param file
    save_dlg = wx.FileDialog(self,
                             message="Save PRIME Script",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      with open(save_dlg.GetPath(), 'w') as savefile:
        savefile.write(txt_out)

  def onIsoRef(self, e):
    if self.prime_panel.ref_box.ctr.GetValue() != '':
      self.prime_panel.opt_chk_useref.Enable()
    else:
      self.prime_panel.opt_chk_useref.Disable()

  def init_settings(self):

    # Determine where/what PRIME folders are
    prime_dir = os.path.join(self.output, 'prime')
    self.working_dir = os.path.join(prime_dir, 'trial_{}'.format(self.trial_no))
    if not os.path.exists(prime_dir):
      os.mkdir(prime_dir)
    if not os.path.exists(self.working_dir):
      os.mkdir(self.working_dir)

    # Write list of pickles to file
    list_prefix = self.opt_prefix.prefix.GetValue()
    if list_prefix == None or list_prefix == '':
      list_prefix = 'prime'
    self.pickle_path_file = os.path.join(self.working_dir,
                           '{}_trial_{}.lst'.format(list_prefix, self.trial_no))
    print('Saving list of pickles to ', self.pickle_path_file)

    with open(self.pickle_path_file, 'w') as lfile:
      for pickle in self.all_pickles:
        lfile.write('{}\n'.format(pickle))

    self.pparams = self.prime_panel.pparams
    self.pparams.data = [self.pickle_path_file]
    self.pparams.run_no = set_base_dir(out_dir=self.working_dir)
    self.out_dir = self.prime_panel.out_box.ctr.GetValue()
    self.pparams.title = self.prime_panel.title_box.ctr.GetValue()
    if str(self.prime_panel.ref_box.ctr.GetValue()).lower() != '':
      self.pparams.hklisoin = self.prime_panel.ref_box.ctr.GetValue()
      if self.prime_panel.opt_chk_useref.GetValue():
        self.pparams.hklrefin = self.prime_panel.ref_box.ctr.GetValue()
    self.pparams.n_residues = self.prime_panel.opt_spc_nres.GetValue()
    self.pparams.n_processors = self.prime_panel.opt_spc_nproc.GetValue()

  def onRun(self, e):
    # Run full processing

    from xfel.util.mp import get_lsf_submit_command
    from xfel.ui import settings_dir
    import datetime
    import copy

    params = copy.deepcopy(self.main.params)
    params.mp.nproc = self.prime_panel.opt_spc_nproc.GetValue()

    # Generate script filename (w/ timestamp)
    ts = '{:%Y%m%d_%H%M%S}'.format(datetime.datetime.now())
    script_filename = 'trial_{:03d}_{}.sh'.format(int(self.trial_no), ts)

    self.init_settings()
    prime_phil = master_phil.format(python_object=self.pparams)

    with Capturing() as output:
      prime_phil.show()

    txt_out = ''
    for one_output in output:
      txt_out += one_output + '\n'

    prime_file = os.path.join(settings_dir, self.prime_filename)
    out_file = os.path.join(self.working_dir, 'stdout.log')
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)

    if params.mp.method == 'local':
      command=None
    else:
      job_name = 'prime_t{}'.format(self.trial_no)
      cmd = '-J {} prime.postrefine {}'.format(job_name, prime_file)
      submit_path = os.path.join(settings_dir, script_filename)
      command = str(get_lsf_submit_command(cmd, submit_path, self.working_dir,
                                           params.mp)())

    if e.GetId() == self.tb_btn_run.GetId():
      self.prime_run_window = PRIMERunWindow(self, -1,
                                             title='PRIME Output',
                                             params=self.pparams,
                                             prime_file=prime_file,
                                             # out_file=out_file,
                                             mp_method=params.mp.method,
                                             command=command)
      self.prime_run_window.prev_pids = easy_run.fully_buffered('pgrep -u {} {}'
                                            ''.format(user, 'python')).stdout_lines

      self.prime_run_window.Show(True)

    elif e.GetId() == self.tb_btn_cmd.GetId():
      print('Submission command:')
      print(command)

    # Try and write files to created folder


# ------------------------------- UI Elements -------------------------------- #

class TrialPanel(wx.Panel):
  ''' A scrolled panel that contains run blocks and trial controls '''

  def __init__(self, parent, db, trial, box_label=None):
    wx.Panel.__init__(self, parent=parent, size=(270, 200))

    self.db = db
    self.trial = trial
    self.parent = parent

    trial_box = wx.StaticBox(self, label=box_label)
    self.main_sizer = wx.StaticBoxSizer(trial_box, wx.VERTICAL)

    self.block_panel = ScrolledPanel(self, size=(150, 180))
    self.block_sizer = wx.BoxSizer(wx.VERTICAL)
    self.block_panel.SetSizer(self.block_sizer)
    self.one_block_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.add_panel = wx.Panel(self)
    self.add_sizer = wx.BoxSizer(wx.VERTICAL)
    self.add_panel.SetSizer(self.add_sizer)

    # Add "New Block" button to a separate sizer (so it is always on bottom)
    self.btn_add_block = wx.Button(self.add_panel, label='New Block',
                                   size=(200, -1))
    self.btn_select_blocks = wx.Button(self.add_panel, label='Select Blocks',
                                       size=(200, -1))
    self.btn_view_phil = gctr.BitmapButton(self.add_panel, name='btn_view_phil',
                                           bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    self.chk_active = wx.CheckBox(self.add_panel, label='Active Trial')
    self.view_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.view_sizer.Add(self.btn_view_phil)
    self.view_sizer.Add(self.chk_active, flag=wx.EXPAND)

    self.add_sizer.Add(self.btn_add_block,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.btn_select_blocks,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.view_sizer,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT,
                       border=10)

    self.main_sizer.Add(self.block_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.add_panel, flag=wx.ALL, border=5)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddBlock, self.btn_add_block)
    self.Bind(wx.EVT_BUTTON, self.onSelectBlocks, self.btn_select_blocks)
    self.Bind(wx.EVT_BUTTON, self.onViewPHIL, self.btn_view_phil)
    self.chk_active.Bind(wx.EVT_CHECKBOX, self.onToggleActivity)

    self.SetSizer(self.main_sizer)

  def onViewPHIL(self, e):
    view_dlg = dlg.TrialDialog(self, db=self.db, trial=self.trial, new=False)
    view_dlg.Fit()
    view_dlg.ShowModal()
    view_dlg.Destroy()

  def onToggleActivity(self, e):
    if self.chk_active.GetValue():
      self.trial.active = True
    else:
      self.trial.active = False

  def onAddBlock(self, e):
    rblock_dlg = dlg.RunBlockDialog(self, trial=self.trial,
                                    db=self.db)
    rblock_dlg.Fit()

    if (rblock_dlg.ShowModal() == wx.ID_OK):
      self.refresh_trial()
    rblock_dlg.Destroy()

  def onSelectBlocks(self, e):
    rblocksel_dlg = dlg.SelectRunBlocksDialog(self, trial=self.trial,
                                           db=self.db)
    rblocksel_dlg.Fit()

    if (rblocksel_dlg.ShowModal() == wx.ID_OK):
      self.refresh_trial()
    rblocksel_dlg.Destroy()

  def refresh_trial(self):
    self.block_sizer.Clear(delete_windows=True)
    self.active_blocks = self.trial.rungroups
    for block in self.active_blocks:
      self.draw_block_button(block)
    self.block_panel.Layout()
    self.block_panel.SetupScrolling(scrollToTop=False)

  def draw_block_button(self, block):
    ''' Add new run block button '''
    new_block = gctr.RunBlock(self.block_panel, block=block)
    self.Bind(wx.EVT_BUTTON, self.onRunBlockOptions, new_block.new_runblock)
    self.block_sizer.Add(new_block,
                         flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                         border=5)

  def onRunBlockOptions(self, e):
    ''' Open dialog and change run_block options '''
    run_block = e.GetEventObject().block
    rblock_dlg = dlg.RunBlockDialog(self, block=run_block,
                                    db=self.db)
    rblock_dlg.Fit()

    if (rblock_dlg.ShowModal() == wx.ID_OK):
      wx.CallAfter(self.refresh_trial)
    rblock_dlg.Destroy()

class DatasetPanel(wx.Panel):
  ''' A scrolled panel that contains dataset and task controls '''

  def __init__(self, parent, db, dataset, box_label=""):
    wx.Panel.__init__(self, parent=parent, size=(270, 200))

    self.db = db
    self.dataset = dataset
    self.parent = parent

    self.dataset_box = wx.StaticBox(self, label=box_label)
    self.main_sizer = wx.StaticBoxSizer(self.dataset_box, wx.VERTICAL)

    self.dataset_comment = wx.StaticText(self)
    self.task_panel = ScrolledPanel(self, size=(150, 180))
    self.task_sizer = wx.BoxSizer(wx.VERTICAL)
    self.task_panel.SetSizer(self.task_sizer)
    self.one_task_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.add_panel = wx.Panel(self)
    self.add_sizer = wx.BoxSizer(wx.VERTICAL)
    self.add_panel.SetSizer(self.add_sizer)

    # Add "New task" button to a separate sizer (so it is always on bottom)
    self.btn_add_task = wx.Button(self.add_panel, label='New Task',
                                   size=(200, -1))
    self.btn_select_tasks = wx.Button(self.add_panel, label='Select Tasks',
                                       size=(200, -1))
    self.btn_edit_dataset = wx.BitmapButton(self.add_panel,
                                            bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    self.chk_active = wx.CheckBox(self.add_panel, label='Active Dataset')
    self.chk_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.chk_sizer.Add(self.btn_edit_dataset)
    self.chk_sizer.Add(self.chk_active, flag=wx.EXPAND)

    self.add_sizer.Add(self.btn_add_task,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.btn_select_tasks,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.chk_sizer,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT,
                       border=10)

    self.main_sizer.Add(self.dataset_comment, 0, flag=wx.ALL, border=10)
    self.main_sizer.Add(self.task_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.add_panel, flag=wx.ALL, border=5)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddTask, self.btn_add_task)
    self.Bind(wx.EVT_BUTTON, self.onSelectTasks, self.btn_select_tasks)
    self.Bind(wx.EVT_BUTTON, self.onEditDataset, self.btn_edit_dataset)
    self.chk_active.Bind(wx.EVT_CHECKBOX, self.onToggleActivity)

    self.SetSizer(self.main_sizer)

  def onToggleActivity(self, e):
    if self.chk_active.GetValue():
      self.dataset.active = True
    else:
      self.dataset.active = False

  def onAddTask(self, e):
    task_dlg = dlg.TaskDialog(self, dataset=self.dataset,
                                    db=self.db)
    task_dlg.Fit()

    if (task_dlg.ShowModal() == wx.ID_OK):
      self.refresh_dataset()
    task_dlg.Destroy()

  def onSelectTasks(self, e):
    tasksel_dlg = dlg.SelectTasksDialog(self, dataset=self.dataset,
                                           db=self.db)
    tasksel_dlg.Fit()

    if (tasksel_dlg.ShowModal() == wx.ID_OK):
      self.refresh_dataset()
    tasksel_dlg.Destroy()

  def refresh_dataset(self):
    self.dataset_comment.SetLabel(self.dataset.comment if self.dataset.comment is not None else "")
    self.dataset_box.SetLabel('Dataset {} {}'.format(self.dataset.dataset_id,
                               self.dataset.name[:min(len(self.dataset.name), 10)]
                               if self.dataset.name is not None else ""))
    self.task_sizer.Clear(delete_windows=True)
    tags = self.dataset.tags
    if tags:
      tags_text = "Tags: " + ",".join([t.name for t in tags])
    else:
      tags_text = "No tags selected"
    label = wx.StaticText(self.task_panel, label = tags_text)
    self.task_sizer.Add(label,
                        flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                        border=5)
    for task in self.dataset.tasks:
      self.draw_task_button(task)
    self.task_panel.Layout()
    self.task_panel.SetupScrolling(scrollToTop=False)

  def draw_task_button(self, task):
    ''' Add new run block button '''
    new_task = gctr.TaskCtrl(self.task_panel, task=task)
    self.Bind(wx.EVT_BUTTON, self.onTaskOptions, new_task.new_task)
    self.task_sizer.Add(new_task,
                        flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                        border=5)

  def onTaskOptions(self, e):
    ''' Open dialog and change task options '''
    task = e.GetEventObject().task
    task_dlg = dlg.TaskDialog(self, task=task,
                                    db=self.db)
    task_dlg.Fit()

    if (task_dlg.ShowModal() == wx.ID_OK):
      wx.CallAfter(self.refresh_dataset)
    task_dlg.Destroy()

  def onEditDataset(self, e):
    new_dataset_dlg = dlg.DatasetDialog(self, db=self.db, dataset=self.dataset, new=False)
    new_dataset_dlg.Fit()

    if new_dataset_dlg.ShowModal() == wx.ID_OK:
      self.refresh_dataset()

class RunEntry(wx.Panel):
  ''' Adds run row to table, with average and view buttons'''
  def __init__(self, parent, run, params, label_width = None):
    self.run = run
    self.params = params

    wx.Panel.__init__(self, parent=parent)
    if label_width is None: label_width = 60

    self.sizer = wx.FlexGridSizer(1, 4, 0, 10)
    run_no = wx.StaticText(self, label=str(run),
                           size=(label_width, -1))
    self.tag_button = gctr.TagButton(self, run=run)
    self.avg_button = wx.Button(self, label='Average')
    self.view_button = wx.Button(self, label='View')
    self.view_button.Hide()

    self.sizer.Add(run_no, flag=wx.EXPAND)
    self.sizer.Add(self.tag_button, flag=wx.EXPAND)
    self.sizer.AddGrowableCol(1)
    self.sizer.Add(self.avg_button)
    self.sizer.Add(self.view_button, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    # Button Bindings
    self.Bind(wx.EVT_BUTTON, self.onTagButton, self.tag_button)
    self.Bind(wx.EVT_BUTTON, self.onAvgButton, self.avg_button)
    self.Bind(wx.EVT_BUTTON, self.onViewButton, self.view_button)

    self.SetSizer(self.sizer)

  def onTagButton(self, e):
    self.tag_button.change_tags()

  def onAvgButton(self, e):
    avg = dlg.AveragingDialog(self, self.run, self.params)
    avg.Fit()
    avg.Center()

    if (avg.ShowModal() == wx.ID_OK):
      e.GetEventObject().SetLabel('Running')
      e.GetEventObject().Disable()
      self.view_button.Show()
      # TODO: hook up the calibration app

  def onViewButton(self, e):
    # TODO: hook up view function
    pass
